
import nlopt
import numpy as np
from scipy.interpolate import interp1d, splrep,splev
from Engine.classes import fitobjs,inparams
from Engine.rotint import rotint
from Engine.macbro_dynamic    import macbro_dyn
from Engine.rebin_jv import rebin_jv
import time
import sys

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
          13: IP linear component
          14: IP quadratic component

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
    sspot2 = interpfunc(watm)

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

    fwhmraw = par[5] + par[13]*(fitobj_cp.x) + par[14]*(fitobj_cp.x**2)
    try:
        spl = splrep(w,fwhmraw)
    except ValueError:
        return 1e7
    fwhm = splev(watm,spl)
    if min(fwhm) < 1 or max(fwhm) > 7:
        return 1e7

    #Handle instrumental broadening
    vhwhm = dw*abs(fwhm)/mnw*c/2.
    #print(min(fwhm),dw*abs(min(fwhm))/mnw*c/2.,max(fwhm),dw*abs(max(fwhm))/mnw*c/2.)
    nsmod = macbro_dyn(vel,smod,vhwhm)

    #Rebin continuum to observed wavelength scale
    #c2 = rebin_jv(fitobj_cp.a0contwave*1e4,fitobj_cp.continuum,w,False)
    #spl = splrep(fitobj_cp.a0contwave*1e4,fitobj_cp.continuum)
    #c2 = splev(w,spl)
    c2 = fitobj_cp.continuum

    #Rebin model to observed wavelength scale
    #smod = rebin_jv(watm,nsmod,w,False)
    spl = splrep(watm,nsmod)
    smod = splev(w,spl)
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    smod *= cont

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj_cp.s < .05)] = False

    # Compute chisq
    chisq = np.sum((fitobj_cp.s[mask] - smod[mask])**2. / fitobj_cp.u[mask]**2.)
    #chisq = np.sum((fitobj_cp.s - smod)**2. / fitobj_cp.u**2.)

    if optimize_cp == True:
        return chisq
    else:
        return smod,chisq

def fmod(par,fitobj):

    watm = fitobj.watm_in;
    satm = fitobj.satm_in;
    mwave = fitobj.mwave_in;
    mflux = fitobj.mflux_in;

    w = par[6] + par[7]*fitobj.x + par[8]*(fitobj.x**2.) + par[9]*(fitobj.x**3.)

    if w[-1] < w[0]:
        sys.exit('WAVE ERROR 1 {}'.format(par[6:10]))
        return 1e7

    c = 2.99792e5
    npts = len(w)

    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        sys.exit('WAVE ERROR 2 {} {} {} {} {}'.format(par[6:10],watm[0],watm[-1],w[0],w[-1]))
        return 1e7

    interpfunc = interp1d(wspot,sspot, kind='linear',bounds_error=False,fill_value='extrapolate')
    sspot2 = interpfunc(watm)

    vsini = abs(par[4])
    if vsini != 0:
        rspot = rotint(watm,sspot2,vsini,eps=.4,nr=5,ntheta=25)
    else:
        rspot = sspot2

    smod = rspot*satm

    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj.x) + par[14]*(fitobj.x**2)
    if min(fwhmraw) < 1 or max(fwhmraw) > 7:
        sys.exit('IP ERROR 1 {} {} {} {} {}'.format(par[5],par[13],par[14],min(fwhmraw),max(fwhmraw) ))
        return 1e7
    try:
        spl = splrep(w,fwhmraw)
    except ValueError:
        sys.exit('IP ERROR 2 {} {} {}'.format(par[5],par[13],par[14]))
        return 1e7
    fwhm = splev(watm,spl)

    vhwhm = dw*abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel,smod,vhwhm)

    c2 = fitobj.continuum
    spl = splrep(watm,nsmod)
    smod = splev(w,spl)
    smod *= c2/np.median(c2)
    cont = par[10] + par[11]*fitobj.x+ par[12]*(fitobj.x**2)
    smod *= cont

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj.s < .05)] = False
    chisq = np.sum((fitobj.s[mask] - smod[mask])**2. / fitobj.u[mask]**2.)

    return smod,chisq

def optimizer(par0, lows, highs, fitobj, optimize):
    # NLopt convenience function.
    global fitobj_cp, optimize_cp
    fitobj_cp   = fitobj
    optimize_cp = optimize
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 15)
    opt.set_min_objective(fmodel_chi)
    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)
    opt.set_maxtime(600) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_xtol_rel(1e-6)
    parfit = opt.optimize(par0)
    return parfit
