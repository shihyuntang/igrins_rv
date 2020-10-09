import nlopt
import numpy as np
from scipy.interpolate import interp1d, splrep,splev
from Engine.classes import fitobjs,inparams
from Engine.rotint import rotint
# from Engine.rotint_old import rotint_old
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
        
          If beam is A:

              15: Blaze dip center location
              16: Blaze dip full width
              17: Blaze dip depth
              18: Secondary blaze dip full width
              19: Blaze dip depth

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
        # print(f'{nc_cp}, {nk_cp}, {optkind_cp}: Hitting negative wavelength solution for some reason !')
        return 1e10

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792458e5
    npts = len(w)

    # Apply velocity shifts and scale
    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    #Verify that new wavelength scale is a subset of old wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        # print(f'{nc_cp}, {nk_cp}, {optkind_cp}: w not subset of watm, w goes from '+str(w[0])+' to '+str(w[-1])+' and watm goes from '+str(watm[0])+' to '+str(watm[-1]))
        return 1e10

    vsini = par[4]


    # Rotationally broaden stellar template
    if vsini != 0:
        wspot2,rspot2 = rotint(wspot,sspot,vsini)
    else:
        wspot2 = wspot
        rspot2 = sspot
    '''
    # Rotationally broaden stellar template
    if vsini != 0:
        rspot2 = rotint_old(wspot,sspot,vsini,eps=.6,nr=5,ntheta=25,dif=None)
    else:
        rspot2 = sspot
    wspot2 = wspot
    '''
    #Now rebin the spot spectrum onto the telluric wavelength scale
    sspot2 = rebin_jv(wspot2,rspot2,watm,False)

    #Mutliply rotationally broadened spot by telluric to create total spectrum
    smod = sspot2*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj_cp.x) + par[14]*(fitobj_cp.x**2)
    try:
        spl = splrep(w,fwhmraw)
    except:
        return 1e10
    fwhm = splev(watm,spl)
    if min(fwhm) < 1 or max(fwhm) > 7:
        return 1e10

    #Handle instrumental broadening
    vhwhm = dw*abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel,smod,vhwhm)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)

    # Load saved continuum
    c2 = fitobj_cp.continuum
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x+ par[12]*(fitobj_cp.x**2)
    if fitobj_cp.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj_cp.x >= (par[15]-par[16]/2)) & (fitobj_cp.x <= (par[15]+par[16]/2))] = par[17]
        bucket[(fitobj_cp.x >= (par[15]+par[16]/2-par[18])) & (fitobj_cp.x <= (par[15]+par[16]/2))] += par[19]
        cont -= bucket
    smod *= cont

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj_cp.s < .0)] = False

    if len(fitobj_cp.mask) != 0:
        for maskbounds in fitobj_cp.mask:
            mask[(fitobj_cp.x > maskbounds[0]) & (fitobj_cp.x < maskbounds[1]) ] = False

    # Compute chisq
    chisq = np.sum((fitobj_cp.s[mask] - smod[mask])**2. / fitobj_cp.u[mask]**2.)
    chisq = chisq / (len(smod[mask]) - 15)

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
        return 1e10

    # c = 2.99792e5
    c = 2.99792458e5
    npts = len(w)

    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        sys.exit('WAVE ERROR 2 {} {} {} {} {}'.format(par[6:10],watm[0],watm[-1],w[0],w[-1]))
        return 1e10

    vsini = par[4]



    # Rotationally broaden stellar template
    if vsini != 0:
        wspot2,rspot2 = rotint(wspot,sspot,vsini)
    else:
        wspot2 = wspot
        rspot2 = sspot
    '''
    # Rotationally broaden stellar template
    if vsini != 0:
        rspot2 = rotint_old(wspot,sspot,vsini,eps=.6,nr=5,ntheta=25,dif=None)
    else:
        rspot2 = sspot
    wspot2 = wspot
    '''
    sspot2 = rebin_jv(wspot2,rspot2,watm,False)

    smod = sspot2*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj.x) + par[14]*(fitobj.x**2)
    if min(fwhmraw) < 1 or max(fwhmraw) > 7:
        sys.exit('IP ERROR 1 {} {} {} {} {}'.format(par[5],par[13],par[14],min(fwhmraw),max(fwhmraw) ))
        return 1e10
    try:
        spl = splrep(w,fwhmraw)
    except ValueError:
        sys.exit('IP ERROR 2 {} {} {}'.format(par[5],par[13],par[14]))
        return 1e10
    fwhm = splev(watm,spl)

    vhwhm = dw*abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel,smod,vhwhm)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)

    # Load saved continuum
    c2 = fitobj.continuum
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj.x+ par[12]*(fitobj.x**2)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (par[15]-par[16]/2)) & (fitobj.x <= (par[15]+par[16]/2))] = par[17]
        bucket[(fitobj.x >= (par[15]+par[16]/2-par[18])) & (fitobj.x <= (par[15]+par[16]/2))] += par[19]
        cont -= bucket
    smod *= cont

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj.s < .0)] = False
    chisq = np.sum((fitobj.s[mask] - smod[mask])**2. / fitobj.u[mask]**2.)
    chisq = chisq / (len(smod[mask]) - 15)

    return smod,chisq

def fmod_conti(par,fitobj):

    watm = fitobj.watm_in;
    satm = fitobj.satm_in;
    mwave = fitobj.mwave_in;
    mflux = fitobj.mflux_in;

    w = par[6] + par[7]*fitobj.x + par[8]*(fitobj.x**2.) + par[9]*(fitobj.x**3.)

    if w[-1] < w[0]:
        sys.exit('WAVE ERROR 1 {}'.format(par[6:10]))
        return 1e10

    # c = 2.99792e5
    c = 2.99792458e5
    npts = len(w)

    wspot = mwave*(1.+par[0]/c)
    sspot = mflux**par[1]
    watm = watm*(1.+par[2]/c)
    satm = satm**par[3]

    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        return 1e10

    vsini = par[4]

    # Rotationally broaden stellar template
    if vsini != 0:
        wspot2,rspot2 = rotint(wspot,sspot,vsini)
    else:
        wspot2 = wspot
        rspot2 = sspot

    sspot2 = rebin_jv(wspot2,rspot2,watm,False)

    smod = sspot2*satm

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (watm-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj.x) + par[14]*(fitobj.x**2)
    if min(fwhmraw) < 1 or max(fwhmraw) > 7:
        sys.exit('IP ERROR 1 {} {} {} {} {}'.format(par[5],par[13],par[14],min(fwhmraw),max(fwhmraw) ))
        return 1e10
    try:
        spl = splrep(w,fwhmraw)
    except ValueError:
        sys.exit('IP ERROR 2 {} {} {}'.format(par[5],par[13],par[14]))
        return 1e10
    fwhm = splev(watm,spl)

    vhwhm = dw*abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel,smod,vhwhm)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(watm,nsmod,w,False)

    # Load saved continuum
    c2 = fitobj.continuum
    smod *= c2/np.median(c2)

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj.x+ par[12]*(fitobj.x**2)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (par[15]-par[16]/2)) & (fitobj.x <= (par[15]+par[16]/2))] = par[17]
        bucket[(fitobj.x >= (par[15]+par[16]/2-par[18])) & (fitobj.x <= (par[15]+par[16]/2))] += par[19]
        cont -= bucket
    smod *= cont

    mask = np.ones_like(smod,dtype=bool)
    mask[(fitobj.s < .05)] = False

    return w, smod, cont, c2


# def optimizer(par0,dpar0, hardbounds_v_ip, fitobj, optimize, logger, night, order, tag, optkind, nc, nk):
def optimizer(par0,dpar0, hardbounds_v_ip, fitobj, optimize):
    # NLopt convenience function.
    global fitobj_cp, optimize_cp
    fitobj_cp   = fitobj
    optimize_cp = optimize

    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 20)
    opt.set_min_objective(fmodel_chi)
    lows  = par0-dpar0
    highs = par0+dpar0

    for frg in [1,3]:
        if dpar0[frg] != 0 and lows[frg] < 0:
            lows[frg] = 0

    if dpar0[4] != 0:
        lows[4] = hardbounds_v_ip[0]; highs[4] = hardbounds_v_ip[1];
        if highs[4]-par0[4] < 1e-4:
            par0[4] = par0[4] - 1e-4
        if par0[4] -lows[4] < 1e-4:
            par0[4] = par0[4] + 1e-4

    if dpar0[5] != 0:
        lows[5] = hardbounds_v_ip[2]; highs[5] = hardbounds_v_ip[3];
        if highs[5]-par0[5] < 1e-4:
            par0[5] = par0[5] - 1e-4
        if par0[5] -lows[5] < 1e-4:
            par0[5] = par0[5] + 1e-4

    if fitobj_cp.masterbeam == 'A':
        if dpar0[15] != 0:
            lows[15] = hardbounds_v_ip[4]; highs[15] = hardbounds_v_ip[5];
            if highs[15]-par0[15] < 1e-4:
                par0[15] = par0[15] - 1e-4
            if par0[15] -lows[15] < 1e-4:
                par0[15] = par0[15] + 1e-4

        if dpar0[16] != 0:
            lows[16] = hardbounds_v_ip[6]; highs[16] = hardbounds_v_ip[7];
            if highs[16]-par0[16] < 1e-4:
                par0[16] = par0[16] - 1e-4
            if par0[16] -lows[16] < 1e-4:
                par0[16] = par0[16] + 1e-4

        if dpar0[17] != 0:
            lows[17] = hardbounds_v_ip[8]; highs[17] = hardbounds_v_ip[9];
            if highs[17]-par0[17] < 1e-4:
                par0[17] = par0[17] - 1e-4
            if par0[17] -lows[17] < 1e-4:
                par0[17] = par0[17] + 1e-4

        if dpar0[18] != 0:
            lows[18] = hardbounds_v_ip[10]; highs[18] = hardbounds_v_ip[11];
            if highs[18]-par0[18] < 1e-4:
                par0[18] = par0[18] - 1e-4
            if par0[18] -lows[18] < 1e-4:
                par0[18] = par0[18] + 1e-4  

        if dpar0[19] != 0:  
            lows[19] = hardbounds_v_ip[12]; highs[19] = hardbounds_v_ip[13];
            if highs[19]-par0[19] < 1e-4:
                par0[19] = par0[19] - 1e-4
            if par0[19] -lows[19] < 1e-4:
                par0[19] = par0[19] + 1e-4  

    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)

    maxruntime = 1200 #seconds
    opt.set_maxtime(1200) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_ftol_rel(1e-10)
    parfit = opt.optimize(par0)
    return parfit
