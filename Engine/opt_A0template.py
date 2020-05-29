
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
    smod = macbro_dyn(vel,rspot,vhwhm)

    #Mutliply rotationally and instrumentally broadened spot by A0 to create total spectrum
    nsmod = smod*satm

    #Rebin model and A0 template uncertainty to observed wavelength scale
    spl = splrep(watm,nsmod)
    smod = splev(w,spl)
    spl = splrep(watm,fitobj_cp.uatm_in)
    uatm = splev(w,spl)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    smod *= cont

    # Add uncertainties in A0 template and target spectrum in quadrature
    u = 1/np.sqrt(1/(fitobj_cp.u**2) + 1/(uatm**2))

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj_cp.s < .05)] = False

    # Compute chisq
    chisq = np.sum((fitobj_cp.s[mask] - smod[mask])**2. / u[mask]**2.)
    #chisq = np.sum((fitobj_cp.s - smod)**2. / u**2.)

    if optimize_cp == True:
        return chisq
    else:
        return smod,chisq

def fmod(par,grad):

    global fitobj_cp, optimize_cp

    watm = fitobj_cp.watm_in;
    satm = fitobj_cp.satm_in;
    mwave = fitobj_cp.mwave_in;
    mflux = fitobj_cp.mflux_in;

    w = par[6] + par[7]*fitobj_cp.x + par[8]*(fitobj_cp.x**2.) + par[9]*(fitobj_cp.x**3.)

    if w[-1] < w[0]:
        print('Hitting negative wavelength solution for some reason !')
        return 1e7

    c = 2.99792e5
    npts = len(w)

    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        print('w not subset of watm, w goes from '+str(w[0])+' to '+str(w[-1])+' and watm goes from '+str(watm[0])+' to '+str(watm[-1]))
        return 1e7

    interpfunc = interp1d(wspot,sspot, kind='linear',bounds_error=False,fill_value='extrapolate')
    sspot2 = interpfunc(watm)

    vsini = abs(par[4])
    if vsini != 0:
        rspot = rotint(watm,sspot2,vsini,eps=.4,nr=5,ntheta=25)
    else:
        rspot = sspot2

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

    vhwhm = dw*abs(fwhm)/mnw*c/2.
    smod = macbro_dyn(vel,rspot,vhwhm)

    nsmod = smod*satm

    spl = splrep(watm,nsmod)
    smod = splev(w,spl)
    spl = splrep(watm,fitobj_cp.uatm_in)
    uatm = splev(w,spl)

    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    smod *= cont

    u = 1/np.sqrt(1/(fitobj_cp.u**2) + 1/(uatm**2))

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj_cp.s < .05)] = False

    chisq = np.sum((fitobj_cp.s[mask] - smod[mask])**2. / u[mask]**2.)
    #chisq = np.sum((fitobj_cp.s - smod)**2. / u**2.)

    return smod,chisq

def optimizer(par0,dpar0, hardbounds_v_ip, fitobj, optimize):
    # NLopt convenience function.
    global fitobj_cp, optimize_cp
    fitobj_cp   = fitobj
    optimize_cp = optimize
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 15)
    opt.set_min_objective(fmodel_chi)
    lows  = par0-dpar0
    highs = par0+dpar0
    for frg in [1,3]:
        if dpar0[frg] != 0 and lows[frg] < 0:
            lows[frg] = 0
    if dpar0[4] != 0:
        lows[4] = hardbounds_v_ip[0]; highs[4] = hardbounds_v_ip[1];
    if dpar0[5] != 0:
        lows[5] = hardbounds_v_ip[2]; highs[5] = hardbounds_v_ip[3];
    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)
    opt.set_maxtime(600) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_xtol_rel(1e-6)
    parfit = opt.optimize(par0)
    return parfit
