from Engine.importmodule import *
from Engine.set_argparse import _argparse_step1

from Engine.IO_AB import setup_templates_tel, init_fitsread, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import a0cont
from Engine.classes import (FitObjs, InParamsA0,
                            OrderDictCla, _setup_bound_cut)
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.Telfitter import telfitter
from Engine.opt import optimizer, fmod
from Engine.outplotter import outplotter_tel
from Engine.detect_peaks import detect_peaks
from Engine.crmask import cr_masker
from Engine.molmask import h2o_masker
#-------------------------------------------------------------------------------

def a0_fits_write(hdu_0, firstorder, order, outpath, night, masterbeam,
                    band):
    """Output telfit generated synthetic telluric template

    Args:
        hdu_0 (astropy.fits.BinTableHDU): Table to be write out
        firstorder (int): The first order to run in this run
        order (int): Current run order
        outpath (str): Output dir
        night (str): Observation night (or + _tag)
        masterbeam (str): Beam (nodding). Can only be A or B
        band (str): H or K band
    """

    if hdu_0 is None:
        # write a error flag = 1
        c0 = fits.Column(name = f'ERRORFLAG{order}',
                            array = np.array([1]),
                            format='K'
                            )
        cols = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)
    else:
        hdu_1 = hdu_0

    if order == firstorder:
        # If first time writing fits file, make up filler primary hdu
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul = fits.HDUList([primary_hdu, hdu_1])
        hdul.writeto('{}/{}A0_{}treated_{}.fits'.format(outpath, night,
                                                        masterbeam, band
                                                        ),
                                                        overwrite=True)

    else:
        hh = fits.open('{}/{}A0_{}treated_{}.fits'.format(outpath, night,
                                                          masterbeam,
                                                          band)
                                                          )
        hh.append(hdu_1)
        hh.writeto('{}/{}A0_{}treated_{}.fits'.format(outpath, night,
                                                      masterbeam, band),
                                                      overwrite=True)


def setup_fitting_init_pars(inparam, night, band, masterbeam, order):
    """Setup the initial values for the parameters to be optimized (fitted)

    Args:
        inparam (class): [description]
        night (str): Observation night (or + _tag)
        band (str): H or K band
        masterbeam (str): Beam (nodding). Can only be A or B
        order (int): Current run order

    Returns:
        np.array: Initial values for the parameters to be optimized
    """

    # Determine whether IGRINS mounting was loose or
    # the night of interest is in question
    if (int(night) < 20180401) or (int(night) > 20190531):
        IPpars = inparam.ips_tightmount_pars[band][masterbeam][order]
    else:
        IPpars = inparam.ips_loosemount_pars[band][masterbeam][order]

    # start at bucket loc = 1250 +- 100, width = 250 +- 100,
    #  depth = 100 +- 5000 but floor at 0
    centerloc = 1250 if band == 'H' else 1180

    # Initialize parameter array for optimization as well as half-range values
    # for each parameter during the various steps of the optimization.
    # Many of the parameters initialized here will be changed throughout the
    # code before optimization and in between optimization steps.

    parA0 = np.array([
        0.0,       # 0: The shift of the stellar template (km/s)
        0.0,       # 1: The scale factor for the stellar template
        0.0,       # 2: The shift of the telluric template (km/s)
        1.0,       # 3: The scale factor for the telluric template
        0.0,       # 4: vsini (km/s)
        IPpars[2], # 5: The instrumental resolution (FWHM) in pixels
        0.0,       # 6: Wavelength 0-pt
        0.0,       # 7: Wavelength linear component
        0.0,       # 8: Wavelength quadratic component
        0.0,       # 9: Wavelength cubic component
        1.0,       #10: Continuum zero point
        0.0,       #11: Continuum linear component
        0.0,       #12: Continuum quadratic component
        IPpars[1], #13: Instrumental resolution linear component
        IPpars[0], #14: Instrumental resolution quadratic component
        centerloc, #15: Blaze dip center location
        330,       #16: Blaze dip full width
        0.05,      #17: Blaze dip depth
        90,        #18: Secondary blaze dip full width
        0.05,      #19: Blaze dip depth
        0.0,       #20: Continuum cubic component
        0.0,       #21: Continuum quartic component
        0.0,       #22: Continuum quintic component
        0.0,       #23: Continuum hexic component
        0.0,       #24: secondary par
        0.0,       #25: secondary par
        0.0,       #26: secondary par
        0.0        #27: secondary par
    ])

    return parA0


    """

    Parameters
    ----------
    use_sets : list with str


    Returns
    -------
    Dict
        dpars
    """

def base_dpars_dict(use_sets, band, order):
    """Setup basic sets of paramaeter variable ranges

    Args:
        use_sets (list with str): List of dpars_org keys that wish to get
        band (str): H or K band
        order (int): Current run order

    Returns:
        dpars (dict): Optimize parameters' variable ranges
        c_order (int): Continuum fitted orders
    """

    #                     | 0    1    2    3 | | 4 | | 5 | | 6    7    8    9 | |10  11 12| |13 14| |15   16   17   18   19|  |20   21   22    23|
    dpars_org = {
        'cont' : np.array([0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  1e7, 1, 1,  0, 0,   0.0, 0.0, 0.0, 0.0, 0.0,  1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
        'twave': np.array([0.0, 0.0, 0.0, 1.0,  0.0,  0.0,  1.0, 1.0, 1.0, 1.0,  0.0, 0, 0,  0, 0,   0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        'ip'   : np.array([0.0, 0.0, 0.0, 0.0,  0.0,  0.5,  0.0, 0.0, 0.0, 0.0,  0.0, 0, 0,  0, 0,   0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        't'    : np.array([0.0, 0.0, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0, 0,  0, 0,   0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        }

    # blaze fitting order setting
    if band == 'H':
        if order in [13]:
            # fit a quadratic (2) continuum
            dpars_org['cont'][20:] = 0
            c_order = 2
        else:
            # fit a hexic (6) continuum (default)
            c_order = 6
            pass
    elif band == 'K':
        if order in [3]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
            c_order = 3
        elif order in [4, 5, 6]:
            # fit a quartic (4) continuum
            dpars_org['cont'][22:] = 0
            c_order = 4
        else:
            # fit a hexic (6) continuum (default)
            c_order = 6
            pass

    dpars = {k: v for k, v in dpars_org.items() if k in use_sets}

    return dpars, c_order


def main(args, inparam, jerp, orders, masterbeam, i):
    """Main function for A0 fitting that will be threaded over
    by multiprocessing
    """

    order = orders[jerp]            # current looped order
    night = str(inparam.nights[i])  # multiprocess assigned night
    firstorder = orders[0]          # First order that will be analyzed,
                                    # related to file writing

    if args.debug:
        print('Working on order {:02d}/{:02d} ({}), ',
                'night {}/{} ({}) PID:{}...'.format(int(jerp+1),
                                                    len(orders),
                                                    order,
                                                    i+1,
                                                    len(inparam.nights),
                                                    night,
                                                    mp.current_process().pid
                                                    ))
    #-------------------------------------------------------------------------------

    bound_cut = _setup_bound_cut(inparam.bound_cut_dic, args.band, order)

    ### Load relevant A0 spectrum
    x, a0wavelist, a0fluxlist, u = init_fitsread(
                                        inparam.inpath,
                                        'A0',
                                        'combined'+str(masterbeam),
                                        night,
                                        order,
                                        f'{int(inparam.tags[night]):04d}',
                                        args.band,
                                        bound_cut
                                        )
    #-------------------------------------------------------------------------------
    try:
        s2n = a0fluxlist/u
        if np.nanmedian(s2n) < float(args.SN_cut):
            logger.warning('  --> Bad S/N {:1.3f} < {} for {}{}, SKIP'.format(
                np.nanmedian(s2n), args.SN_cut, night, masterbeam)
                )
            pre_err = True
            logger.warning(f'  --> NIGHT {night}, ORDER {order} '
                                'HIT ERROR DURING PRE_OPT')
            a0_fits_write(None, firstorder, order, inparam.outpath, night,
                            masterbeam, args.band)
            return

    except ZeroDivisionError:
        logger.warning('  --> There must be something wrong with '
                            'flux error = 0 for {}{}, SKIP'.format(
                                                    night, masterbeam)
                                                    )
        pre_err = True
        logger.warning(f'  --> NIGHT {night}, ORDER {order} '
                            'HIT ERROR (flux error = 0) DURING PRE_OPT')
        a0_fits_write(None, firstorder, order, inparam.outpath, night,
                        masterbeam, args.band)
        return

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    # Do twice
    nzones = 12
    a0wavelist = basicclip_above(a0wavelist, a0fluxlist, nzones)
    a0x = basicclip_above(x, a0fluxlist, nzones)
    a0u        = basicclip_above(u, a0fluxlist, nzones)
    a0fluxlist = basicclip_above(a0fluxlist, a0fluxlist, nzones)

    a0wavelist = basicclip_above(a0wavelist, a0fluxlist, nzones)
    a0x = basicclip_above(a0x, a0fluxlist, nzones)
    a0u        = basicclip_above(a0u, a0fluxlist, nzones)
    a0fluxlist = basicclip_above(a0fluxlist, a0fluxlist, nzones)

    if masterbeam == 'B':
        # Compute rough blaze function estimate.
        # Better fit will be provided by Telfit later.
        continuum = a0cont(a0wavelist, a0fluxlist, night, order, args.band)
        a0masterwave = a0wavelist.copy()
        a0masterwave *= 1e4 # um --> AA

        # Normalize continuum level of telluric atlas in the given band
        if args.band == 'H':
            contlevel = np.max(
                inparam.satm[ (inparam.watm > 15000) & (inparam.watm < 18000) ]
                )
        else:
            contlevel = np.max(
                inparam.satm[ (inparam.watm > 20000) & (inparam.watm < 24000) ]
                )

        # Trim telluric template to relevant wavelength range
        satm_in = inparam.satm[
            (inparam.watm > np.min(a0wavelist)*1e4 - 11) &
            (inparam.watm < np.max(a0wavelist)*1e4 + 11)
            ]
        watm_in = inparam.watm[
            (inparam.watm > np.min(a0wavelist)*1e4 - 11) &
            (inparam.watm < np.max(a0wavelist)*1e4 + 11)
            ]
        satm_in /= contlevel

        a0x = a0x[ (a0wavelist*1e4 > np.min(watm_in)+5) &
                   (a0wavelist*1e4 < np.max(watm_in)-5)
                ]
        continuum = continuum[ (a0wavelist*1e4 > np.min(watm_in)+5) &
                                (a0wavelist*1e4 < np.max(watm_in)-5)
                            ]

    elif masterbeam == 'A':

        A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_Btreated_{args.band}.fits'

        try:
            hdulist = fits.open(A0loc)
            # Find corresponding table in fits file, given the tables do not go
            # sequentially by order number due to multiprocessing in Step 1
            num_orders = 0
            for i in range(25):
                try:
                    hdulist[i].columns[0].name[9:]
                    num_orders += 1
                except:
                    continue
            fits_layer = [ i for i in np.arange(num_orders)+1 if
                            int(hdulist[i].columns[0].name[9:]) == order ][0]
            tbdata = hdulist[ fits_layer ].data
            flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

            # Check whether Telfit hit critical error in Step 1 for the chosen
            # order with this night. If so, skip.
            if flag == 1:
                logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR '
                                f'IN ORDER: {order} NIGHT: {night}, skipping...')
                a0_fits_write(None, firstorder, order, inparam.outpath,
                                night, masterbeam, args.band)
                return

        except IOError:
            logger.warning(f'  --> No A0-fitted template for night {night}, '
                                'skipping...')
            a0_fits_write(None, firstorder, order, inparam.outpath, night,
                            masterbeam, args.band)
            return

        watm = tbdata['WATM'+str(order)]
        satm = tbdata['SATM'+str(order)]
        a0contx    = tbdata['X'+str(order)]
        continuum  = tbdata['BLAZE'+str(order)]

        # Remove extra rows leftover from having columns of unequal length
        satm = satm[(watm != 0)]
        watm = watm[(watm != 0)]
        # Set very low points to zero so that they don't go to NaN when taken
        # to an exponent by template power in fmodel_chi
        satm[(satm < 1e-4)] = 0.
        a0contx = a0contx[(continuum != 0)]
        continuum = continuum[(continuum != 0)]

        satm_in_fine = satm[(watm > np.min(a0wavelist)*1e4 - 11) &
                       (watm < np.max(a0wavelist)*1e4 + 11)
                       ]
        watm_in_fine = watm[(watm > np.min(a0wavelist)*1e4 - 11) &
                       (watm < np.max(a0wavelist)*1e4 + 11)
                       ]

        coarse_wavesep = 0.10 #AA
        coarse_nstep   = int((watm_in_fine[-1]-watm_in_fine[0])/coarse_wavesep)

        watm_in = np.linspace(watm_in_fine[0],watm_in_fine[-1],coarse_nstep)
        satm_in = rebin_jv(watm_in_fine,satm_in_fine,watm_in,False)
        satm_in[(satm_in < 1e-4)] = 0.

        a0x = a0x[ (a0wavelist*1e4 > np.min(watm_in)+5) &
                   (a0wavelist*1e4 < np.max(watm_in)-5)
                   ]
        continuum = rebin_jv(a0contx, continuum, a0x, False)


    a0fluxlist = a0fluxlist[ (a0wavelist*1e4 > np.min(watm_in)+5) &
                             (a0wavelist*1e4 < np.max(watm_in)-5)
                            ]
    a0u = a0u[ (a0wavelist*1e4 > np.min(watm_in)+5) &
               (a0wavelist*1e4 < np.max(watm_in)-5)
               ]
    a0wavelist = a0wavelist[ (a0wavelist*1e4 > np.min(watm_in)+5) &
                             (a0wavelist*1e4 < np.max(watm_in)-5)
                             ]

    # Initialize parameter array for optimization
    parA0 = setup_fitting_init_pars(inparam, night, args.band,
                                        masterbeam, order)

    # Define main spectrum
    s = a0fluxlist.copy()
    x = a0x.copy()
    u = a0u.copy()

    # Check if we're including part of blaze where dip occurs
    fitdip = False
    if masterbeam == 'A':
        if ( (x[0] < parA0[15]-parA0[16]) and (x[-1] > parA0[15]-parA0[16]) ) \
            or ( (x[0] < parA0[15]+parA0[16]) and (x[-1] > parA0[15]+parA0[16]) ):
            fitdip = True

    # Get initial guess for cubic wavelength solution from reduction pipeline
    f = np.polyfit(a0x, a0wavelist, 3)
    q = np.poly1d(f)
    initwave = q(a0x)*1e4

    # Collect all fit variables into one class
    fitobj = FitObjs(s, x, u, continuum, watm_in, satm_in, None, None, [],
                        masterbeam,
                        [np.array([], dtype=int), np.array([], dtype=int)],
                        initwave, [])

    # Initialize an array that puts hard bounds on vsini and the instrumental
    #  resolution to make sure they do not diverge to unphysical values
    optimize = True
    par_in = parA0.copy()

    # setup fitting boundary
    dpars, c_order = base_dpars_dict(['cont', 'twave', 'ip'],
                                        args.band, order)
    if fitdip:
        # setting for the Blaze dip
        dpars['cont'][15:20] = np.array([10.0, 30.0, 0.2, 50.0, 0.2])

        hardbounds = [
            par_in[4]  - 0,                 par_in[4]  + 0,
            par_in[5]  - dpars['ip'][5],    par_in[5]  + dpars['ip'][5],
            par_in[15] - dpars['cont'][15], par_in[15] + dpars['cont'][15],
            par_in[16] - dpars['cont'][16], par_in[16] + dpars['cont'][16],
            0.,                             par_in[17] + dpars['cont'][17],
            par_in[18] - dpars['cont'][18], par_in[18] + dpars['cont'][18],
            0.,                             par_in[19] + dpars['cont'][19]
            ]
    else:
        hardbounds = [
            par_in[4]  - 0,                 par_in[4]  + 0,
            par_in[5]  - dpars['ip'][5],    par_in[5]  + dpars['ip'][5]
            ]

    if hardbounds[0] < 0.5: hardbounds[0] = 0.5
    if hardbounds[2] < 1:   hardbounds[2] = 1

    # Begin optimization.
    # For every pre-Telfit spectral fit, first fit just template
    # strength/rv/continuum, then just wavelength solution, then template/continuum
    # again, then ip, then finally wavelength. Normally would fit for all but
    # wavelength at the end, but there's no need for the pre-Telfit fit, since
    # all we want is a nice wavelength solution to feed into Telfit.

    optgroup = ['twave', 'cont',
                'twave', 'cont',
                'twave', 'cont']

    # Try 3 different levels of telluric template power to start with. Save
    # the chisqs and best fit pars from each, compare, and then resume
    # optimization from the best pars.
    chisqs = []
    parfitsaves = []

    for telval in [0.5, 1.0, 1.5]:

        parstart = par_in.copy()
        parstart[3] = telval

        nk = 1
        for optkind in optgroup:
            parfit_1 = optimizer(parstart, dpars[optkind], hardbounds,
                                    fitobj, optimize)
            if args.debug == True:
                outplotter_tel(parfit_1, fitobj,
                                '{}_{}_{}_beforeTelfit_telvalstart_{}_{}{}'.format(
                                    order, night, masterbeam, telval, nk, optkind
                                    ), inparam, args, order)
            parstart = parfit_1.copy()
            nk += 1

        smod,chisq,trash,trash2 = fmod(parfit_1,fitobj,False)
        chisqs.append(chisq)
        parfitsaves.append(parfit_1)

    parfit = parfitsaves[np.argmin(chisqs)]

    # If telluric lines too shallow, intentionally cause exception to trigger except
    if parfit[3] < 0.2 or np.min(chisqs) > 1e9:
        logger.warning(f'  --> NIGHT {night}, ORDER {order} HAD TOO LITTLE '
                            'ABSORPTION OR DEVIATED FROM LIVINGSTON TOO MUCH')
        a0_fits_write(None, firstorder, order, inparam.outpath, night,
                        masterbeam, args.band)
        return

     # If A beam and fitting cont dip, do one more cycle of optimization
    if fitdip:
        parstart = parfit.copy()
        for optkind in optgroup:
            parfit_1 = optimizer(parstart, dpars[optkind], hardbounds,
                                    fitobj, optimize)
            if args.debug == True:
                outplotter_tel(parfit_1, fitobj,
                                '{}_{}_{}_beforeTelfit_{}{}'.format(
                                    order, night, masterbeam, nk, optkind
                                    ), inparam, args, order)
            parstart = parfit_1.copy()
            nk += 1
        parfit = parfit_1.copy()

        # If dip present, correct it out of data before running Telfit
        # to enable better fit
        cont = parfit[10] + parfit[11]*fitobj.x+ parfit[12]*(fitobj.x**2) \
                + parfit[20]*(fitobj.x**3) + parfit[21]*(fitobj.x**4) \
                + parfit[22]*(fitobj.x**5) + parfit[23]*(fitobj.x**6)
        cont0 = cont.copy()
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (parfit[15]-parfit[16]/2)) \
            & (fitobj.x <= (parfit[15]+parfit[16]/2))] = parfit[17]

        bucket[(fitobj.x >= (parfit[15]+parfit[16]/2-parfit[18])) \
            & (fitobj.x <= (parfit[15]+parfit[16]/2))] += parfit[19]

        cont -= bucket

        cont *= continuum
        cont0 *= continuum
        justdip = cont/cont0
        a0fluxlistforTelfit = a0fluxlist / justdip
    else:
        a0fluxlistforTelfit = a0fluxlist.copy()

    # Plot results
    if inparam.plotfigs:
        outplotter_tel(parfit, fitobj,
                        f'BeforeTelFit_Order{order}_{night}_{masterbeam}',
                        inparam, args, order)

    # Get best fit wavelength solution
    xgrid = (initwave - np.median(initwave)) / (np.max(initwave) - np.min(initwave))
    dx = chebyshev.chebval(xgrid, parfit[6:10])
    a0w_out_fit = initwave + dx

    fwhmraw = parfit[5] + parfit[13]*(x) + parfit[14]*(x**2)
    resolution_max = np.max(
        a0w_out_fit)/(np.min(fwhmraw)*np.min(np.diff(a0w_out_fit)) )
    # not max because pixels get skipped
    resolution_min = np.min(
        a0w_out_fit)/(np.max(fwhmraw)*np.median(np.diff(a0w_out_fit)) )
    resolution_med = np.median(
        a0w_out_fit)/(np.median(fwhmraw)*np.median(np.diff(a0w_out_fit)) )

    resolutions = [resolution_min, resolution_med, resolution_max]

    # Feed this new wavelength solution into Telfit. Returns high-res synthetic
    # telluric template, parameters of that best fit, and blaze function best fit
    watm1, satm1, telfitparnames, telfitpars, a0contwave, continuum, molnames, molwaves, molfluxes = telfitter(
        a0w_out_fit, a0fluxlistforTelfit, a0u, inparam, night, order, args, masterbeam, c_order, resolutions, logger
        )

    # If Telfit encountered error (details in Telfitter.py), skip night/order combo
    if len(watm1) == 1:
        logger.warning(
            f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} '
                f'NIGHT: {night}')
        a0_fits_write(None, firstorder, order, inparam.outpath, night,
                        masterbeam, args.band)
        return

    else: # If Telfit exited normally, proceed.
        #  Save best blaze function fit
        continuum = rebin_jv(a0contwave, continuum, a0w_out_fit, False)

        # Write out table to fits file with errorflag = 0
        c0  = fits.Column(name='ERRORFLAG'+str(order),      array=np.array([0]),            format='K')
        c1  = fits.Column(name='WAVE'+str(order),           array=a0w_out_fit,              format='D')
        c2  = fits.Column(name='BLAZE'+str(order),          array=continuum,                format='D')
        c3  = fits.Column(name='X'+str(order),              array=a0x,                      format='D')
        c4  = fits.Column(name='INTENS'+str(order),         array=a0fluxlist,               format='D')
        c5  = fits.Column(name='SIGMA'+str(order),          array=a0u,                      format='D')
        c6  = fits.Column(name='WATM'+str(order),           array=watm1,                    format='D')
        c7  = fits.Column(name='SATM'+str(order),           array=satm1,                    format='D')
        c8  = fits.Column(name='TELFITPARNAMES'+str(order), array=np.array(telfitparnames), format='8A')
        c9  = fits.Column(name='TELFITPARS'+str(order),     array=telfitpars,               format='D')
        c10 = fits.Column(name='PARFIT',                    array=parfit,                   format='D')
        c11 = fits.Column(name='MOLNAMES',                  array=np.array(molnames),       format='3A')
        collist = [c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]
        for mol in molnames:
            collist.append(fits.Column(name='WATM'+mol+str(order), array=molwaves[mol],   format='D'))
            collist.append(fits.Column(name='SATM'+mol+str(order), array=molfluxes[mol],  format='D'))
        cols = fits.ColDefs(collist)
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        a0_fits_write(hdu_1, firstorder, order, inparam.outpath,
                        night, masterbeam, args.band)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def use_w(args):
    """Load wavelength regions list file
    """
    try:
        bounddata = Table.read(
            f'./Input/UseWv/WaveRegions_{args.WRegion}_{args.band}.csv',
            format='csv')
    except IOError:
        sys.exit(
            f'WaveRegions FILE "./Input/UseWv/WaveRegions'
                '_{args.WRegion}_{args.band}.csv" NOT FOUND!')

    wavesols = pd.read_csv(f'./Input/UseWv/WaveSolns_{args.band}.csv')
#-------------------------------------------------------------------------------
    XRegion_dir = f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv'
    with open(XRegion_dir,'w') as filew:
        filew.write('order, start,  end, masks\n')

        m_order = np.array(bounddata['order'])
        starts = np.array(bounddata['start'])
        ends = np.array(bounddata['end'])
        ords = list( sorted(OrderDictCla().orderdict[args.band].keys()) )

        Ostarts = [OrderDictCla().orderdict[args.band][k][0] for k in ords]
        Oends = [OrderDictCla().orderdict[args.band][k][1] for k in ords]
        labels = []

        m_orders_unique = np.unique(m_order)

        # For each order specified, find what pixel numbers correspond to the
        # wavelength bounds presented.
        # If multiple wavelength bounds given for a single order, output a
        # pixel mask between the two, as well.
        for o in range(len(m_orders_unique)):

            if len(m_orders_unique) == 9:
                filew.write('9, 150, 1950, []\n')
                continue

            pixs = []
            mini = np.where(m_order == m_orders_unique[o])[0]
            for j in range(len(mini)):
                i = mini[j]

                wavebounds = [starts[i],ends[i]]
                wO   = wavesols['w'+str(m_orders_unique[o])]
                pixO = wavesols['x'+str(m_orders_unique[o])]
                pix  = [pixO[(np.argmin(abs(wO-wavebounds[k])))] for k in [0,1]]
                pixs = pixs + pix

            pixsS = list(sorted(pixs))
            q = pixsS[1:-1]
            if len(pixsS) == 2:
                filew.write('{}, {}, {},[]\n'.format(
                    m_orders_unique[o], pixsS[0], pixsS[-1])
                    )
            else:
                filew.write('{}, {}, {},"{}"\n'.format(
                    m_orders_unique[o], pixsS[0], pixsS[-1],
                    [[first,second] for first, second in zip(q[0::2], q[1::2])]
                    ))

#-------------------------------------------------------------------------------
#-------------------------------------------------  ------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    args = _argparse_step1()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'
#-------------------------------------------------------------------------------
    # Create output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}')

    outpath = f'./Output/{args.targname}_{args.band}/A0Fits'
#-------------------------------------------------------------------------------
    # Handle logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}_A0Fits.log')
    stream_hander= logging.StreamHandler()

    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('####################################################################################\n')
    print(f'Fetching Wavelength Regions to be Analyzed for {args.targname}...')
    time.sleep(2)

    use_w(args)

    print('Fetching Done!')
    print(
        f'File "XRegions_{args.WRegion}_{args.band}.csv" saved under '
            '"./Input/UseWv/"')
    time.sleep(2)

    #-------------------------------------------------------------------------------

    print('###############################################################\n')
    logger.info(f'Using TelFit to create high-resolution, synthetic telluric '
                    'templates based off the telluric standards \nassociated '
                    f'with {args.targname} on a night by night basis...')
    print('This will take a while..........')

    # Read in newly created pixel regions file to get list of orders to analyze.
    # Note that only target star observations will have their fits limited to
    # the wavelength regions specified.
    # For A0 observations, only the orders specified will be analyzed, but each
    # order will be fit as far as there is significant telluric absoprtion.
    bounddata = Table.read(
        f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv', format='csv')
    starts  = np.array(bounddata['start'])
    ends    = np.array(bounddata['end'])
    orders  = np.array(bounddata['order'], dtype=int)
    xbounddict = {orders[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}

    ## Collect relevant file information from Predata files
    A0data = Table.read(f'./Input/Prepdata/Prepdata_A0_{args.targname}.txt',
    format='ascii')

    ind    = [i != 'NA' for i in A0data['humid']]
    humids = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['humid'])}
    tags   = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['tag'])}
    obs    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['obs'])}
    temps  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['temp'])}
    zds    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['zd'])}
    press  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['press'])}
    nightsFinal = np.array(list(sorted(set(A0data[ind]['night']))))

    # Take subset of nights, if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=int)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(f'NIGHT {nnn} EITHER HAS NO CORRESPONDING A0 OR WAS '
                            f'NOT FOUND UNDER "./Input/{args.targname}"')

        nightsFinal = nightstemp
        logger.info(f'Only processing nights: {nightsFinal}')

    logger.info(f'Analyze {len(nightsFinal)} nights')
    intnights = np.array( [np.int(i) for i in nightsFinal] )
    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these nights
                slightly differently. When you run Step 3, RVs will be output
                in two formats: one with the defocus nights separated, and the
                other with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as the
                result of the defocus. IGRINS RV will process defocus nights
                the same way as the others, but when you run Step 3, will still
                output the results in two formats like it does with the K band.
''')

    # orders = np.array([6])
    #--------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates_tel()

    inparam = InParamsA0(inpath, outpath, args.plotfigs, tags, nightsFinal,
                            humids, temps, zds, press, obs, watm, satm,
                            mwave0, mflux0, cdbs_loc, xbounddict, None)

    #--------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    print('Processing the B nods first...')
    for jerp in range(len(orders)):
        if not args.debug:
            print('Working on order {} ({:02d}/{:02d})'.format(
                                orders[jerp], int(jerp+1), len(orders)))
        # main( args, inparam, jerp, orders, 'B',0)
        func = partial(main, args, inparam, jerp, orders, 'B')
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

    print('B nods done! Halfway there! \n Now processing the A nods...')
    for jerp in range(len(orders)):
        if not args.debug:
            print('Working on order {} ({:02d}/{:02d})'.format(
                                orders[jerp], int(jerp+1), len(orders)))
        func = partial(main, args, inparam, jerp, orders, 'A')
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

    warning_r = log_warning_id(
        f'{outpath}/{args.targname}_{args.band}_A0Fits.log', start_time)
    if warning_r:
        print(f'''
**********************************************************************************
WARNING!! you got warning message during this run. Please check the log file under:
          {outpath}/{args.targname}_{args.band}_A0Fits.log
**********************************************************************************
''')
    print('\n')
    if not args.debug: logger.addHandler(stream_hander)
    logger.info('A0 Fitting Done!')

    end_time = datetime.now()
    logger.info(f'A0 Fitting using TelFit finished, Duration: {end_time - start_time}')
    print('The synthetic telluric templates have been saved under {}'.format(outpath))
    print('If you chose to generate plots, they are saved under {}/figs'.format(outpath))
    print('####################################################################################')
    print('You can now run main_step2.py to produce RV and vsini initial guess(es)')
    print('####################################################################################')
