
import nlopt
import numpy as np
from scipy.interpolate import interp1d
from Engine.classes import fitobjs,inparams
from Engine.rotint import rotint
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv

#-------------------------------------------------------------------------------
def fmodel_chi(par,grad):
    '''
    INPUTS:
       w - The observed wavelength scale (air) in Angstroms.
       x - The array of pixel indices from 0 to npts-1
       par - The model parameters:
          0: The shift of the sunspot spectrum (km/s)
          1: The scale factor for the sunspot spectrum
          2: The shift of the telluric spectrum (km/s)
          3: The scale factor for the telluric spectrum
          4: vsini (km/s)
          5: The instrumental resolution (FWHM) in pixels
          6: Wavelength 0-pt
          7: Wavelength linear component
          8: Wavelength quadratic component
          9: Wavelength cubic component
          10: Continuum zero point
          11: Continuum linear component
          12: Continuum quadratic component

     OUTPUTS:
       The model spectrum on the observed wavelength scale.
    '''

    # Bring in global class of variables needed to generate model.
    # Can't call these directly in function, as NLopt doesn't allow anything to be in the model function call besides par and grad.

    #global fitobj, optimize
    global fitobj_cp, optimize_cp

    watm = fitobj_cp.watm_in;
    satm = fitobj_cp.satm_in;
    mwave = fitobj_cp.mwave_in;
    mflux = fitobj_cp.mflux_in;

    #Make the wavelength scale
    w = par[6] + par[7]*fitobj_cp.x + par[8]*(fitobj_cp.x**2.) + par[9]*(fitobj_cp.x**3.)

    if w[-1] < w[0]:
        #print('Hitting negative wavelength solution for some reason !')
        return 1e7

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792e5
    npts = len(w)

    # Apply velocity shifts and scale
    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    #Verify that new wavelength scale is a subset of old wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        # print('w not subset of watm, w goes from '+str(w[0])+' to '+str(w[-1])+' and watm goes from '+str(watm[0])+' to '+str(watm[-1]))
        return 1e7

    #Now interpolate the spot spectrum onto the telluric wavelength scale
    interpfunc = interp1d(wspot,sspot, kind='linear',bounds_error=False,fill_value='extrapolate')
    sspot2=interpfunc(watm)

    #Handle rotational broadening
    vsini = abs(par[4])
    if vsini != 0:
        rspot = rotint(watm,sspot2,vsini,eps=.4,nr=5,ntheta=25)
    else:
        rspot = sspot2

    #Mutliply rotationally broadened spot by telluric to create total spectrum
    smod = rspot*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    #Handle instrumental broadening
    vhwhm = dw*abs(par[5])/mnw*c/2.
    nsmod = macbro(vel,smod,vhwhm)

    #Rebin continuum to observed wavelength scale
    c2 = rebin_jv(fitobj_cp.a0contwave*1e4,fitobj_cp.continuum,w,False)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    smod *= cont

    # Compute chisq
    chisq = np.sum((fitobj_cp.s - smod)**2. / fitobj_cp.u**2.)

    if optimize_cp == True:
        return chisq
    else:
        return smod,chisq


def fmodel_separate(par):
    global fitobj

    watm = fitobj_cp.watm_in;
    satm = fitobj_cp.satm_in;
    mwave = fitobj_cp.mwave_in;
    mflux = fitobj_cp.mflux_in;

    #Make the wavelength scale
    w = par[6] + par[7]*fitobj_cp.x + par[8]*(fitobj_cp.x**2.) + par[9]*(fitobj_cp.x**3.)

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792e5
    npts = len(w)

    # Apply velocity shifts and scale
    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    #Now interpolate the spot spectrum onto the telluric wavelength scale
    interpfunc = interp1d(wspot,sspot, kind='linear',bounds_error=False,fill_value='extrapolate')
    sspot2=interpfunc(watm)

    #Handle rotational broadening
    vsini = abs(par[4])
    if vsini != 0:
        rspot = rotint(watm,sspot2,vsini,eps=.4,nr=5,ntheta=25)
    else:
        rspot = sspot2

    #Mutliply rotationally broadened spot by telluric to create total spectrum
    smod = rspot*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    #Handle instrumental broadening
    vhwhm = dw*abs(par[5])/mnw*c/2.
    nsmod = macbro(vel,smod,vhwhm)

    #Rebin continuum to observed wavelength scale
    c2 = rebin_jv(fitobj_cp.a0contwave*1e4,fitobj_cp.continuum,w,False)
    # Apply continuum adjustment
    #c2 /= np.median(c2)
    cont1 = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    cont = cont1 * c2

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)


    return w,smod,cont,cont1


def fmodel_chi_plot(par,grad):
    # Bring in global class of variables needed to generate model.
    # Can't call these directly in function, as NLopt doesn't allow anything to be in the model function call besides par and grad.

    #global fitobj, optimize
    global fitobj_cp
    optimize_cp = False

    watm = fitobj_cp.watm_in;
    satm = fitobj_cp.satm_in;
    mwave = fitobj_cp.mwave_in;
    mflux = fitobj_cp.mflux_in;

    #Make the wavelength scale
    w = par[6] + par[7]*fitobj_cp.x + par[8]*(fitobj_cp.x**2.) + par[9]*(fitobj_cp.x**3.)

    if w[-1] < w[0]:
        #print('Hitting negative wavelength solution for some reason !')
        return 1e7

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792e5
    npts = len(w)

    # Apply velocity shifts and scale
    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    #Verify that new wavelength scale is a subset of old wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        # print('w not subset of watm, w goes from '+str(w[0])+' to '+str(w[-1])+' and watm goes from '+str(watm[0])+' to '+str(watm[-1]))
        return 1e7

    #Now interpolate the spot spectrum onto the telluric wavelength scale
    interpfunc = interp1d(wspot,sspot, kind='linear',bounds_error=False,fill_value='extrapolate')
    sspot2=interpfunc(watm)

    #Handle rotational broadening
    vsini = abs(par[4])
    if vsini != 0:
        rspot = rotint(watm,sspot2,vsini,eps=.4,nr=5,ntheta=25)
    else:
        rspot = sspot2

    #Mutliply rotationally broadened spot by telluric to create total spectrum
    smod = rspot*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    #Handle instrumental broadening
    vhwhm = dw*abs(par[5])/mnw*c/2.
    nsmod = macbro(vel,smod,vhwhm)

    #Rebin continuum to observed wavelength scale
    c2 = rebin_jv(fitobj_cp.a0contwave*1e4,fitobj_cp.continuum,w,False)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    smod *= cont

    # Compute chisq
    chisq = np.sum((fitobj_cp.s - smod)**2. / fitobj_cp.u**2.)

    if optimize_cp == True:
        return chisq
    else:
        return smod,chisq


def optimizer(par0, dpar0, fitobj, optimize):
    # NLopt convenience function.
    global fitobj_cp, optimize_cp
    fitobj_cp   = fitobj
    optimize_cp = optimize
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 13)
    opt.set_min_objective(fmodel_chi)
    lows  = par0-dpar0
    highs = par0+dpar0
    # Don't let template powers or vsini be negative
    for frg in [1,3,4]:
        if dpar0[frg] != 0:
            lows[frg] = 0
    for frg in [5]:
        if dpar0[frg] != 0: # Don't even let IP hwhm hit zero (bc throws error)
            lows[frg] = 0.1
    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)
    opt.set_maxtime(600) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_xtol_rel(1e-6)
    parfit = opt.optimize(par0)
    return parfit
