import sys
import matplotlib
import matplotlib.pyplot as plt
import nlopt
import numpy as np
from numpy.polynomial import chebyshev

from scipy.interpolate import splrep,splev #, interp1d
from Engine.classes import FitObjs,InParams
from Engine.rotint import rotint
from Engine.macbro_dynamic    import macbro_dyn
from Engine.rebin_jv import rebin_jv

#-------------------------------------------------------------------------------
def fmodel_chi(par,grad):
    '''
    Function to be optimized. Computes model spectrum and compares it with
    data to calculate reduced chisq.

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
          15: Blaze dip center location        \
          16: Blaze dip full width              |
          17: Blaze dip depth                   | <-- If beam is A
          18: Secondary blaze dip full width    |
          19: Blaze dip depth                  /
          20: Continuum cubic component      \
          21: Continuum quartic component     |  <-- Only enabled for some orders, depending on size of region being fit
          22: Continuum pentic component      |
          23: Continuum hexic component      /
          24: The shift of the secondary stellar spectrum (km/s)
          25: The scale factor for the secondary stellar spectrum
          26: vsini for the secondary stellar spectrum
          27: Flux ratio fo secondary to primary, S2/S1

     OUTPUTS:
       chisq: Reduced chisq
       (If global optimize_cp is False) smod: The model spectrum on the observed
       wavelength scale.
    '''

    # Bring in global class of variables needed to generate model.
    # Can't call these directly in function, as NLopt doesn't allow anything to
    # be in the model function call besides par and grad.

    #global fitobj, optimize
    global fitobj_cp, optimize_cp, binary_cp, rebin2to1_cp

    watm_in = fitobj_cp.watm_in
    satm_in = fitobj_cp.satm_in
    mwave = fitobj_cp.mwave_in
    mflux = fitobj_cp.mflux_in
    if binary_cp:
        mwave2 = fitobj_cp.mwave_in2
        mflux2 = fitobj_cp.mflux_in2

    #Make the wavelength scale
    initwave = fitobj_cp.initwave.copy()
    xgrid = (initwave-np.median(initwave)) / (np.max(initwave)-np.min(initwave))
    dx = chebyshev.chebval(xgrid, par[6:10])
    w = initwave + dx

    if np.all(np.diff(w) > 0) == False:
        # print(f'{nc_cp}, {nk_cp}, {optkind_cp}: Hitting negative wavelength solution for some reason !')
        return 1e10

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792458e5
    npts = len(w)

    # Apply velocity shifts and scale
    watm = watm_in*(1.+par[2]/c)
    satm = satm_in**par[3]

    #Verify that new wavelength scale is a subset of old telluric wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        #print('WAVE ERROR 2: w subset of watm, w goes from ' +
        #            str(w[0]) + ' to '+str(w[-1]) +
        #            ' and watm goes from ' + str(watm[0]) +
        #            ' to ' + str(watm[-1])
        #            )
        return 1e10

    if mwave is not None:

        wspot = mwave*(1.+par[0]/c)
        sspot = mflux**par[1]

        # Verify that new wavelength scale is a subset of stellar wavelength scale.
        if (w[0] < wspot[0]) or (w[-1] > wspot[-1]):
            #print('WAVE ERROR 3:  w not subset of wspot, w goes from '
            #            + str(w[0]) + ' to '+str(w[-1]) +
            #            ' and wspot goes from ' + str(wspot1[0]) +
            #            ' to '+str(wspot1[-1])
            #            )
            return 1e10

        vsini = par[4]

        # Rotationally broaden stellar template
        if vsini >= 0.5:
            wspot1,rspot1 = rotint(wspot, sspot, vsini)
        else:
            wspot1 = wspot
            rspot1 = sspot

        if binary_cp:

            wspot2 = mwave2*(1.+par[24]/c)
            sspot2 = mflux2**par[25]

            vsini2 = par[26]
            # Rotationally broaden stellar template
            if vsini2 >= 0.5:
                wspot2,rspot2 = rotint(wspot2, sspot2, vsini2)
            else:
                wspot2 = wspot2
                rspot2 = sspot2

            rspot2   *=  par[27]

            if fitobj_cp.rebin2to1:
                # Verify that new wavelength scale is a subset of stellar wavelength scale.
                if (wspot2[0] > wspot1[0]) or (wspot2[-1] < wspot1[-1]):
                    return 1e10

                rspot2n = rebin_jv(wspot2,rspot2,wspot1,True)
                rspot   = rspot1 + rspot2n
                wspot   = wspot1.copy()
            else:
                # Verify that new wavelength scale is a subset of stellar wavelength scale.
                if (wspot2[0] < wspot1[0]) or (wspot2[-1] > wspot1[-1]):
                    return 1e10

                rspot1n = rebin_jv(wspot1,rspot1,wspot2,True)
                rspot   = rspot2 + rspot1n
                wspot   = wspot2.copy()

            rspot /= (1+par[27])

        else:
            rspot = rspot1; wspot = wspot1;

        if (wspot[0] < watm[0]) or (wspot[-1] > watm[-1]):
            #print('WAVE ERROR 5: wspot not subset of satm, wspot goes from '
            #            + str(wspot[0]) + ' to ' + str(wspot[-1]) +
            #            ' and watm goes from ' + str(watm[0]) +
            #            ' to ' + str(watm[-1]))
            return 1e10
            
        # Now rebin the telluric spectrum onto the stellar wavelength scale
        satm2 = rebin_jv(watm, satm, wspot, False)

        # Mutliply rotationally broadened spot by telluric to create total spectrum
        smod2 = rspot*satm2
    else:
        smod2 = satm.copy()
        wspot = watm.copy()

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (wspot-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj_cp.x) + par[14]*(fitobj_cp.x**2)
    try:
        spl = splrep(w, fwhmraw)
    except:
        return 1e10
    fwhm = splev(wspot,spl)

    # Have IP extend as constant past wave bounds of data
    fwhm[(wspot < w[0])]  = fwhm[(wspot >= w[0])][0]
    fwhm[(wspot > w[-1])] = fwhm[(wspot <= w[-1])][-1]
    if (np.min(fwhm) < 1) or (np.max(fwhm) > 8):
        return 1e10

    #Handle instrumental broadening
    vhwhm = dw*np.abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel, smod2, vhwhm)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(wspot, nsmod ,w, False)

    # Load saved continuum
    c2 = fitobj_cp.continuum
    smod *= c2

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj_cp.x + par[12]*(fitobj_cp.x**2) \
            + par[20]*(fitobj_cp.x**3) + par[21]*(fitobj_cp.x**4) \
            + par[22]*(fitobj_cp.x**5) + par[23]*(fitobj_cp.x**6)
    if fitobj_cp.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj_cp.x >= (par[15]-par[16]/2))         \
                & (fitobj_cp.x <= (par[15]+par[16]/2))] = par[17]
        bucket[(fitobj_cp.x >= (par[15]+par[16]/2-par[18])) \
                & (fitobj_cp.x <= (par[15]+par[16]/2))] += par[19]
        cont -= bucket
    smod *= cont

    mask = np.ones_like(smod, dtype=bool)
    mask[(fitobj_cp.s < .0)] = False

    if len(fitobj_cp.mask) != 0:
        for maskbounds in fitobj_cp.mask:
            mask[(fitobj_cp.x > maskbounds[0]) \
                 & (fitobj_cp.x < maskbounds[1]) ] = False

    if len(fitobj_cp.molmask) > 0:
        for mb in fitobj_cp.molmask:
            mask[(fitobj_cp.x >= mb[0]) & (fitobj_cp.x <= mb[1])] = False

    if len(fitobj_cp.CRmask[1]) > 0:
        for mb in fitobj_cp.CRmask[1]:
            mask[(fitobj_cp.x >= fitobj_cp.CRmask[0][mb]-1) \
                    & (fitobj_cp.x <= fitobj_cp.CRmask[0][mb]+1)] = False


    # Compute chisq
    chisq = np.sum((fitobj_cp.s[mask] - smod[mask])**2. / fitobj_cp.u[mask]**2.)
    chisq = chisq / (len(smod[mask]) - len(par))

    if optimize_cp == True:
        return chisq
    else:
        return smod,chisq


matplotlib.use('Qt5Agg')

def fmod(par,fitobj,binary=False):
    '''
    Same as fmodel_chi(), but meant to provide best fit model, not for
    optimization. Always returns both smod and chisq.
    '''

    watm_in = fitobj.watm_in
    satm_in = fitobj.satm_in
    mwave = fitobj.mwave_in
    mflux = fitobj.mflux_in
    if binary:
        mwave2 = fitobj.mwave_in2
        mflux2 = fitobj.mflux_in2

    #Make the wavelength scale
    initwave = fitobj.initwave.copy()
    xgrid = (initwave-np.median(initwave)) / (np.max(initwave)-np.min(initwave))
    dx = chebyshev.chebval(xgrid, par[6:10])
    w = initwave + dx

    if np.all(np.diff(w) > 0) == False:
        sys.exit('WAVE ERROR 1 - Hitting negative wavelength solution ',
                    'for some reason - pars: {}'.format(par[6:10]))
        return 1e10

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792458e5
    npts = len(w)

    # Apply velocity shifts and scale
    watm = watm_in*(1.+par[2]/c)
    satm = satm_in**par[3]

    #Verify that new wavelength scale is a subset of old telluric wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        sys.exit('WAVE ERROR 2: w subset of watm, w goes from ' +
                    str(w[0]) + ' to '+str(w[-1]) +
                    ' and watm goes from ' + str(watm[0]) +
                    ' to ' + str(watm[-1])
                    )
        return 1e10

    if mwave is not None:

        wspot1 = mwave*(1.+par[0]/c)
        sspot1 = mflux**par[1]

        #Verify that new wavelength scale is a subset of stellar wavelength scale.
        if (w[0] < wspot1[0]) or (w[-1] > wspot1[-1]):
            sys.exit('WAVE ERROR 3:  w not subset of wspot, w goes from '
                        + str(w[0]) + ' to '+str(w[-1]) +
                        ' and wspot goes from ' + str(wspot1[0]) +
                        ' to '+str(wspot1[-1])
                        )
            return 1e10

        vsini = par[4]

        # Rotationally broaden stellar template
        if vsini >= 0.5:
            wspot1,rspot1 = rotint(wspot1,sspot1,vsini)
        else:
            wspot1 = wspot1
            rspot1 = sspot1


        if binary:

            wspot2 = mwave2*(1.+par[24]/c)
            sspot2 = mflux2**par[25]

            vsini2 = par[26]
            # Rotationally broaden stellar template
            if vsini2 >= 0.5:
                wspot2,rspot2 = rotint(wspot2, sspot2, vsini2)
            else:
                wspot2 = wspot2
                rspot2 = sspot2

            rspot2   *=  par[27]

            if fitobj_cp.rebin2to1:
                # Verify that new wavelength scale is a subset of stellar wavelength scale.
                if (wspot2[0] > wspot1[0]) or (wspot2[-1] < wspot1[-1]):
                    sys.exit('WAVE ERROR 4: wspot1 not subset of wspot2, wspot2 goes from '
                                + str(wspot2[0]) + ' to ' + str(wspot2[-1]) +
                                ' and wspot1 goes from ' + str(wspot1[0]) +
                                ' to ' + str(wspot1[-1]))
                    return 1e10

                rspot2n = rebin_jv(wspot2,rspot2,wspot1,True)
                rspot   = rspot1 + rspot2n
                wspot   = wspot1.copy()
            else:
                # Verify that new wavelength scale is a subset of stellar wavelength scale.
                if (wspot2[0] < wspot1[0]) or (wspot2[-1] > wspot1[-1]):
                    sys.exit('WAVE ERROR 4: wspot2 not subset of wspot1, wspot2 goes from '
                                + str(wspot2[0]) + ' to ' + str(wspot2[-1]) +
                                ' and wspot1 goes from ' + str(wspot1[0]) +
                                ' to ' + str(wspot1[-1]))
                    return 1e10

                rspot1n = rebin_jv(wspot1,rspot1,wspot2,True)
                rspot   = rspot2 + rspot1n
                wspot   = wspot2.copy()

            rspot /= (1+par[27])

        else:
            rspot = rspot1; wspot = wspot1;

        #Verify that stellar wavelength scale is a subset of telluric wavelength scale.
        if (wspot[0] < watm[0]) or (wspot[-1] > watm[-1]):
            sys.exit('WAVE ERROR 5: wspot not subset of satm, wspot goes from '
                        + str(wspot[0]) + ' to ' + str(wspot[-1]) +
                        ' and watm goes from ' + str(watm[0]) +
                        ' to ' + str(watm[-1]))
            return 1e10

        #Now rebin the telluric spectrum onto the stellar wavelength scale
        satm2 = rebin_jv(watm, satm, wspot, False)

        #Mutliply rotationally broadened spot by telluric to create total spectrum
        smod2 = rspot*satm2

    else:
        smod2 = satm.copy(); wspot = watm.copy()

    #Find mean observed wavelength and create a telluric velocity scale
    mnw = np.mean(w)
    dw = (w[-1] - w[0])/(npts-1.)
    vel = (wspot-mnw)/mnw*c

    fwhmraw = par[5] + par[13]*(fitobj.x) + par[14]*(fitobj.x**2)
    try:
        spl = splrep(w, fwhmraw)
    except:
        sys.exit(
            'IP ERROR 1 {} {} {} {} {}'.format(
                par[5], par[13], par[14], np.min(fwhmraw), np.max(fwhmraw)
                ))
        return 1e10
    fwhm = splev(wspot, spl)

    # Have IP extend as constant past wave bounds of data
    fwhm[(wspot < w[0])]  = fwhm[(wspot >= w[0])][0]
    fwhm[(wspot > w[-1])] = fwhm[(wspot <= w[-1])][-1]
    if (np.min(fwhm) < 1) or (np.max(fwhm) > 8):
        sys.exit(
            'IP ERROR 2 {} {} {} {} {}'.format(
                par[5], par[13], par[14], np.min(fwhm), np.max(fwhm)
                ))
        return 1e10

    #Handle instrumental broadening
    vhwhm = dw*np.abs(fwhm)/mnw*c/2.
    nsmod = macbro_dyn(vel, smod2, vhwhm)

    #Rebin model to observed wavelength scale
    smod = rebin_jv(wspot, nsmod, w, False)

    # Load saved continuum
    c2 = fitobj.continuum
    smod *= c2

    # Apply continuum adjustment
    cont = par[10] + par[11]*fitobj.x + par[12]*(fitobj.x**2) \
            + par[20]*(fitobj.x**3) + par[21]*(fitobj.x**4) \
            + par[22]*(fitobj.x**5) + par[23]*(fitobj.x**6)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (par[15]-par[16]/2)) \
                    & (fitobj.x <= (par[15]+par[16]/2))] = par[17]
        bucket[(fitobj.x >= (par[15]+par[16]/2-par[18])) \
                    & (fitobj.x <= (par[15]+par[16]/2))] += par[19]
        cont -= bucket
    smod *= cont

    mask = np.ones_like(smod, dtype=bool)
    mask[(fitobj.s < .0)] = False

    if len(fitobj.mask) != 0:
        for maskbounds in fitobj.mask:
            mask[(fitobj.x > maskbounds[0]) \
                    & (fitobj.x < maskbounds[1]) ] = False

    if len(fitobj.molmask) > 0:
        for mb in fitobj.molmask:
            mask[(fitobj.x >= mb[0]) \
                    & (fitobj.x <= mb[1])] = False

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(fitobj.x >= fitobj.CRmask[0][mb]-1) \
                    & (fitobj.x <= fitobj.CRmask[0][mb]+1)] = False


    # Compute chisq
    chisq = np.sum((fitobj.s[mask] - smod[mask])**2. / fitobj.u[mask]**2.)
    chisq = chisq / (len(smod[mask]) - len(par))

    return smod,chisq, w, cont*c2





# def optimizer(par0,dpar0, hardbounds_v_ip, fitobj, optimize, logger, night, order, tag, optkind, nc, nk):
def optimizer(par0, dpar0, hardbounds_v_ip, fitobj, optimize, binary=False):
    '''
    Prepares and applies NLOpt optimization to find spectral model that
    best fits data

    Inputs:
    par0            : Initial guesses for spectral model parameters
    dpar0           : Amount optimizer can vary spectral model parameters from
                        initial guesses
    hardbounds_v_ip : List of absolute upper and lower boundaries for vsini
                        and IP width
    fitobj          : Class containing spectral data to be fit and templates
                        to be used for fitting
    optimize        : Boolean specifying whether optimization will occur or not

    Outputs:
    parfit   : Best-fit spectral model parameters
    '''

    # NLopt convenience function.
    global fitobj_cp, optimize_cp, binary_cp
    fitobj_cp   = fitobj
    optimize_cp = optimize
    binary_cp   = binary

    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 28)
    opt.set_min_objective(fmodel_chi)
    lows  = par0-dpar0
    highs = par0+dpar0

    for frg in [1,3]:
        if dpar0[frg] != 0 and lows[frg] < 0:
            lows[frg] = 0

    if dpar0[4] != 0:
        lows[4] = hardbounds_v_ip[0]; highs[4] = hardbounds_v_ip[1]
        if highs[4]-par0[4] < 1e-4:
            par0[4] = par0[4] - 1e-4
        if par0[4] -lows[4] < 1e-4:
            par0[4] = par0[4] + 1e-4

    if dpar0[5] != 0:
        lows[5] = hardbounds_v_ip[2]; highs[5] = hardbounds_v_ip[3]
        if highs[5]-par0[5] < 1e-4:
            par0[5] = par0[5] - 1e-4
        if par0[5] -lows[5] < 1e-4:
            par0[5] = par0[5] + 1e-4

    if fitobj_cp.masterbeam == 'A':
        if dpar0[15] != 0:
            lows[15] = hardbounds_v_ip[4]; highs[15] = hardbounds_v_ip[5]
            if highs[15]-par0[15] < 1e-4:
                par0[15] = par0[15] - 1e-4
            if par0[15] -lows[15] < 1e-4:
                par0[15] = par0[15] + 1e-4

        if dpar0[16] != 0:
            lows[16] = hardbounds_v_ip[6]; highs[16] = hardbounds_v_ip[7]
            if highs[16]-par0[16] < 1e-4:
                par0[16] = par0[16] - 1e-4
            if par0[16] -lows[16] < 1e-4:
                par0[16] = par0[16] + 1e-4

        if dpar0[17] != 0:
            lows[17] = hardbounds_v_ip[8]; highs[17] = hardbounds_v_ip[9]
            if highs[17]-par0[17] < 1e-4:
                par0[17] = par0[17] - 1e-4
            if par0[17] -lows[17] < 1e-4:
                par0[17] = par0[17] + 1e-4

        if dpar0[18] != 0:
            lows[18] = hardbounds_v_ip[10]; highs[18] = hardbounds_v_ip[11]
            if highs[18]-par0[18] < 1e-4:
                par0[18] = par0[18] - 1e-4
            if par0[18] -lows[18] < 1e-4:
                par0[18] = par0[18] + 1e-4

        if dpar0[19] != 0:
            lows[19] = hardbounds_v_ip[12]; highs[19] = hardbounds_v_ip[13]
            if highs[19]-par0[19] < 1e-4:
                par0[19] = par0[19] - 1e-4
            if par0[19] -lows[19] < 1e-4:
                par0[19] = par0[19] + 1e-4

    if binary:
        if dpar0[25] != 0 and lows[25] < 0:
            lows[25] = 0
        if dpar0[26] != 0:
            lows[26] = hardbounds_v_ip[-2]; highs[26] = hardbounds_v_ip[-1]
            if highs[26]-par0[26] < 1e-4:
                par0[26] = par0[26] - 1e-4
            if par0[26] -lows[26] < 1e-4:
                par0[26] = par0[26] + 1e-4

    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)

    maxruntime = 1200 #seconds
    opt.set_maxtime(maxruntime)
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_ftol_rel(1e-11)
    parfit = opt.optimize(par0)
    return parfit
