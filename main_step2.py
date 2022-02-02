from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step2

from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import a0cont
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.crmask    import cr_masker
from Engine.molmask    import h2o_masker

#-------------------------------------------------------------------------------
#---------------------------------------------------------------------------
def setup_fitting_init_pars(band, initvsini, order, initvsini2=0, fluxratio=0):
    """Setup the initial values for the parameters to be optimized (fitted)

    Args:
        band (str): H or K band
        initvsini (float): Initial vsini value
        order (int): Current run order

    Returns:
        np.array: Initial values for the parameters to be optimized
    """

    # start at bucket loc = 1250 +- 100, width = 250 +- 100,
    # depth = 100 +- 5000 but floor at 0
    centerloc = 1250 if band == 'H' else 1180

    # Initialize parameter array for optimization as well as half-range values
    # for each parameter during the various steps of the optimization.
    # Many of the parameters initialized here will be changed throughout the
    # code before optimization and in between optimization steps.

    pars0 = np.array([
        np.nan,    # 0: The shift of the stellar template (km/s) [assigned later]
        0.3,       # 1: The scale factor for the stellar template
        0.0,       # 2: The shift of the telluric template (km/s)
        0.6,       # 3: The scale factor for the telluric template
        initvsini, # 4: vsini (km/s)
        np.nan,    # 5: The instrumental resolution (FWHM) in pixels
        0.0,       # 6: Wavelength 0-pt
        0.0,       # 7: Wavelength linear component
        0.0,       # 8: Wavelength quadratic component
        0.0,       # 9: Wavelength cubic component
        1.0,       #10: Continuum zero point
        0.0,       #11: Continuum linear component
        0.0,       #12: Continuum quadratic component
        np.nan,    #13: Instrumental resolution linear component
        np.nan,    #14: Instrumental resolution quadratic component
        centerloc, #15: Blaze dip center location
        330,       #16: Blaze dip full width
        0.05,      #17: Blaze dip depth
        90,        #18: Secondary blaze dip full width
        0.05,      #19: Blaze dip depth
        0.0,       #20: Continuum cubic component
        0.0,       #21: Continuum quartic component
        0.0,       #22: Continuum pentic component
        0.0,       #23: Continuum hexic component
        np.nan,    #24: The shift of the second stellar template (km/s) [assigned later]
        0.3,       #25: The scale factor for the second stellar template
        initvsini2,#26: Secondary vsini (km/s)
        fluxratio  #27: Secondary to primary flux ratio S2/S1 (km/s)
    ])

    if int(order) == 13: pars0[1] = 0.8

    return pars0



def base_dpars_dict(vsinivary, band, order, run_num, vsinivary2=-1):
    """Setup basic sets of paramaeter variable ranges

    Args:
        initvsini (float): Initial vsini value
        band (str): H or K band
        order (int): Current run order
        run_num (int): Number of the optimize sequence that is being running

    Returns:
        dpars_org (dict): Sets of optimize parameters' variable ranges
    """

    #                     | 0    1    2    3 |  | -- 4 -- || 5 | | 6     7     8     9 | |10  11  12| |13 14||15   16   17   18    19 | |20   21   22   23 | 24 25 26
    dpars_org = {
        'cont' : np.array([0.0, 0.0, 0.0, 0.0,  0.0,        0.0,  0.0,  0.0,  0.0,  0.0,  1e7, 1, 1,   0, 0,  10., 30., 0.2, 50.0, 0.2,  1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ]),
        'twave': np.array([0.0, 0.0, 0.0, 1.0,  0.0,        0.0,  1.0,  1.0,  1.0,  1.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]),
        'ip'   : np.array([0.0, 0.0, 0.0, 0.0,  0.0,        0.5,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]),
        's'    : np.array([20.0, 1.0, 0.0, 0.0,  0.0,        0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]),
        'v'    : np.array([0.0, 0.0, 0.0, 0.0,  vsinivary1,  0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]),
        'ts'   : np.array([20.0, 1.0, 0.0, 1.0,  0.0,        0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
        }

    if vsinivary2 != -1:
        dpars_org['s2']   = np.array([0.0, 0.0, 0.0, 0.0,  0.0,        0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 20.0, 1.0, 0.0, 0.0 ])
        dpars_org['v2']   = np.array([0.0, 0.0, 0.0, 0.0,  0.0,         0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, vsinivary2, 0.0 ])
        dpars_org['s1s2'] = np.array([5.0, 1.0, 0.0, 0.0,  0.0,        0.0,  0.0,  0.0,  0.0,  0.0,  0,   0, 0,   0, 0,   0.,  0., 0.0,  0.,  0.0,  0.0, 0.0, 0.0, 0.0, 20.0, 1.0, 0.0, 0.0 ])


    # blaze fitting order setting
    if band == 'H':
        if order in [13]:
            # fit a quadratic (2) continuum
            dpars_org['cont'][20:] = 0
        elif order in [6,14,21]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        else:
            pass
    elif band == 'K':
        if order in [3,5]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        elif order in [4,6]:
            # fit a quartic (4) continuum
            dpars_org['cont'][22:] = 0
        else:
            pass

    if run_num == 2:
        dpars_org['s'][0] = 5.0 # The shift of the stellar template (km/s)
        dpars_org['ts'][0] = 5.0 # The shift of the stellar template (km/s)

    return dpars_org





def main(args, inparam, orders, order_use, trk, step2or3, i):
    """Main function for RV fitting that will be threaded over
    by multiprocessing
    """

    nights   = inparam.nights
    night = nights[i] # current looped night

    order = order_use
    xbounds = inparam.xbounddict[order]

    if args.debug:
        print('Working on order {:02d}, night {:03d}/{:03d} ',
                '({}) PID:{}...'.format(int(order),
                                        i+1,
                                        len(inparam.nights),
                                        night,
                                        mp.current_process().pid) )

    #-------------------------------------------------------------------------------

    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == float:
        initguesses = inparam.initguesses
    else:
        sys.exit(
            'ERROR! EXPECTING SINGLE NUMBER OR FILE FOR INITGUESSES! QUITTING!'
            )

    if np.isnan(initguesses) == True:
        logger.warning(
            f'  --> Previous run of {night} found it inadequate, skipping...'
            )
        return night, np.nan, np.nan

    if args.binary:
        if type(inparam.initguesses2) == dict:
            initguesses2 = inparam.initguesses2[night]
        elif type(inparam.initguesses2) == float:
            initguesses2 = inparam.initguesses2
        else:
            sys.exit('ERROR! EXPECTING SINGLE NUMBER OR FILE FOR '
                        'INITGUESSES 2! QUITTING!')

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = np.array([])
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight = np.append(beamsnight, 'B')

    # Only do B exposures, and just use first B nodding
    masterbeam = 'B'; beam = 'B'
    try:
        tag  = tagsnight[0]
    except IndexError:
        logger.warning(f'  --> No B nodding(frame) for night {night}, skipping...')
        return night, np.nan, np.nan

    if args.binary:
        pars0 = setup_fitting_init_pars(args.band, inparam.initvsini, order, inparam.initvsini2, float(args.fluxratio))
    else:
        pars0 = setup_fitting_init_pars(args.band, inparam.initvsini, order)

    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

    try:
        hdulist = fits.open(A0loc)
    except IOError:
        logger.warning(
            f'  --> No A0-fitted template for night {night}, skipping...'
            )
        return night, np.nan, np.nan

    # Find corresponding table in fits file, given the tables do not go
    # sequentially by order number due to multiprocessing in Step 1
    num_orders = 0
    for i in range(25):
        try:
            hdulist[i].columns[0].name[9:]
            num_orders += 1
        except:
            continue

    fits_layer = [ i for i in np.arange(num_orders)+1 \
                    if np.int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

    # Check whether Telfit hit critical error in Step 1 for the chosen order
    # with this night. If so, try another order. If all hit the error, skip the night.
    nexto = 0
    ordertry = order
    while 1 == 1:
        fits_layer = [ i for i in np.arange(num_orders)+1 \
                        if np.int(hdulist[i].columns[0].name[9:]) == ordertry ][0]

        tbdata = hdulist[ fits_layer ].data
        flag = np.array(tbdata[f'ERRORFLAG{ordertry}'])[0]

        # If Telfit hit unknown critical error in Step 1, this order can't
        # be used for this night. Try another.
        if flag == 1:
            orderbad = ordertry
            ordertry = orders[nexto]
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: '
                                f'{orderbad} NIGHT: {night}, TRYING ORDER '
                                f'{ordertry} INSTEAD...')

        else: # All good, continue
            order = ordertry
            break

        nexto += 1
        if nexto == len(orders):
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ALL '
                                f'ORDERS FOR NIGHT: {night}, skipping...')
            return night, np.nan, np.nan


    # Use instrumental profile dictionary corresponding to whether IGRINS
    # mounting was loose or not
    if np.int(night[:8]) < 20180401 or np.int(night[:8]) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][masterbeam][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][masterbeam][order]

    watm = tbdata['WATM'+str(order)]
    satm = tbdata['SATM'+str(order)]
    a0contx    = tbdata['X'+str(order)]
    continuum  = tbdata['BLAZE'+str(order)]
    molnames   = tbdata['MOLNAMES']

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]

    # set very low points to zero so that they don't go to NaN when taken
    # to an exponent by template power in fmodel_chi
    satm[(satm < 1e-4)] = 0.
    a0contx = a0contx[(continuum != 0)]
    continuum = continuum[(continuum != 0)]

    #-------------------------------------------------------------------------------
    bound_cut = _setup_bound_cut(inparam.bound_cut_dic, args.band, order)

    # Load target spectrum

    x,wave,s,u = init_fitsread(
        f'{inparam.inpath}/',
        'target',
        'combined'+str(masterbeam),
        night,
        order,
        inparam.tagsB[night][0],
        args.band,
        bound_cut
        )

    # Execute S/N cut
    s2n = s/u
    if np.nanmedian(s2n) < np.float(args.SN_cut):
        logger.warning(
            '  --> Bad S/N {:1.3f} < {} for {}{} {}... '.format(
                np.nanmedian(s2n), args.SN_cut, night, beam, tag
                ))
        pass

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    nzones = 5
    x = basicclip_above(x,s,nzones)
    wave = basicclip_above(wave,s,nzones)
    u = basicclip_above(u,s,nzones)
    s = basicclip_above(s,s,nzones)
    x = basicclip_above(x,s,nzones)
    wave = basicclip_above(wave,s,nzones)
    u = basicclip_above(u,s,nzones)
    s = basicclip_above(s,s,nzones)

    # Cut spectrum to within wavelength regions defined in input list
    s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
    x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

    # Save data for second template cutting after optimization cycle 1 done
    s_save = s_piece.copy(); x_save = x_piece.copy(); u_save = u_piece.copy()

    # Trim telluric template to data range +- 15 AA. If telluric template
    # buffer is cut short because A0 lines didn't extend
    # far past data range, cut data range accordingly.
    satm_in = satm[(watm > np.min(wave_piece)*1e4 - 10) \
                        & (watm < np.max(wave_piece)*1e4 + 10)]
    watm_in = watm[(watm > np.min(wave_piece)*1e4 - 10) \
                        & (watm < np.max(wave_piece)*1e4 + 10)]

    s_piece	= s_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                            & (wave_piece*1e4 < np.max(watm_in)-10)]
    u_piece	= u_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                            & (wave_piece*1e4 < np.max(watm_in)-10)]
    x_piece	= x_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                            & (wave_piece*1e4 < np.max(watm_in)-10)]
    wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+10) \
                            & (wave_piece*1e4 < np.max(watm_in)-10)]

    Rstell1 = np.median(np.diff(inparam.mwave0))

    if args.binary:
        Rstell2 = np.median(np.diff(inparam.mwave2))

        if  Rstell1 > Rstell2:
            rebin2to1 = True; extra1 = 0.; extra2 = 10.;
        else:
            rebin2to1 = False; extra1 = 10.; extra2 = 0.;

        mflux_in2 = inparam.mflux2[
                (inparam.mwave2 > np.min(wave_piece)*1e4 - 5 - extra2) \
                & (inparam.mwave2 < np.max(wave_piece)*1e4 + 5 + extra2)
            ]
        mwave_in2 = inparam.mwave2[
                (inparam.mwave2 > np.min(wave_piece)*1e4 - 5 - extra2) \
                & (inparam.mwave2 < np.max(wave_piece)*1e4 + 5 + extra2)
            ]
        Rstell = np.min([Rstell1,Rstell2])

        dstep = Rstell2
        nstep = int((mwave_in2[-1]-mwave_in2[0])/dstep)
        mwave1 = np.linspace(mwave_in2[0],mwave_in2[-1],nstep)
        mflux1 = rebin_jv(mwave_in2,mflux_in2,mwave1,False)
        mwave_in2 = mwave1.copy(); mflux_in2 = mflux1.copy()
        mwave_in2 = mwave_in2[1:-1]
        mflux_in2 = mflux_in2[1:-1]
    else:
        extra1 = 0; extra2 = 0; Rstell = Rstell1;

    # Trim stellar template to data range +- 10 AA
    mflux_in = inparam.mflux0[
            (inparam.mwave0 > np.min(wave_piece)*1e4 - 5 - extra1) \
            & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5 + extra1)
        ]
    mwave_in = inparam.mwave0[
            (inparam.mwave0 > np.min(wave_piece)*1e4 - 5 - extra1) \
            & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5 + extra1)
        ]
    Rtell = np.median(np.diff(watm_in))
    if Rstell < Rtell:
        sys.exit(f'Telluric template resolution ({round(Rtell,4)} AA) '
                    'must be finer than stellar template resolution '
                    '({round(Rstell,4)} AA) !')

    # Rebin stellar template to uniform wavelength scale
    dstep = Rstell1
    nstep = int((mwave_in[-1]-mwave_in[0])/dstep)
    mwave1 = np.linspace(mwave_in[0],mwave_in[-1],nstep)
    mflux1 = rebin_jv(mwave_in,mflux_in,mwave1,False)
    mwave_in = mwave1.copy(); mflux_in = mflux1.copy()
    mwave_in = mwave_in[1:-1]
    mflux_in = mflux_in[1:-1]

    # Normalize continuum from A0 to flux scale of data
    continuum /= np.nanmedian(continuum)
    continuum *= np.nanpercentile(s_piece,99)

    # --------------------------------------------------------------

    par = pars0.copy()

    # Get initial guess for cubic wavelength solution from reduction pipeline
    f = np.polyfit(x_piece, wave_piece, 3)
    q = np.poly1d(f)
    initwave = q(x_piece)*1e4

    # Initial RV with barycentric correction
    par[0] = initguesses-inparam.bvcs[night+tag]
    par[5]  = IPpars[2]
    par[13] = IPpars[1]
    par[14] = IPpars[0]
    if args.binary:
        par[24] = initguesses2-inparam.bvcs[night+tag]
    # setup fitting boundary
    if args.binary:
        dpars1 = base_dpars_dict(inparam.vsinivary, args.band, int(order), run_num=1, inparam.vsinivary2)
        dpars2 = base_dpars_dict(inparam.vsinivary, args.band, int(order), run_num=2, inparam.vsinivary2)
    else:
        dpars1 = base_dpars_dict(inparam.vsinivary, args.band, int(order), run_num=1)
        dpars2 = base_dpars_dict(inparam.vsinivary, args.band, int(order), run_num=2)


    continuum_in = rebin_jv(a0contx, continuum, x_piece, False)
    fitobj = FitObjs(
        s_piece, x_piece, u_piece, continuum_in, watm_in, satm_in,
        mflux_in, mwave_in, ast.literal_eval(inparam.maskdict[order]),
        masterbeam, [np.array([],dtype=int),np.array([],dtype=int)],
        initwave, [])

    if args.binary:
        fitobj.addsecondary(mwave_in2,mflux_in2,rebin2to1)
    #-------------------------------------------------------------------------------

    # Initialize an array that puts hard bounds on vsini and the instrumental
    # resolution to make sure they do not diverge to unphysical values
    optimize = True
    par_in = par.copy()
    hardbounds = [par_in[4] - dpars1['v'][4],  par_in[4] + dpars1['v'][4],
                  par_in[5] - dpars1['ip'][5], par_in[5] + dpars1['ip'][5]
                 ]

    if hardbounds[0] < 0.5:
        hardbounds[0] = 0.5
    if hardbounds[2] < 1:
        hardbounds[2] = 1

    if args.binary:
        hardbounds.append(par_in[26] - dpars['v2'][26])
        hardbounds.append(par_in[26] + dpars['v2'][26])
        if hardbounds[-2] < 0.5:
            hardbounds[-2] = 0.5
    # Begin optimization. Fit the blaze, the wavelength solution, the telluric
    # template power and RV, the stellar template power and RV, the
    # zero point for the instrumental resolution, and the vsini of the star
    # separately, iterating and cycling between each set of parameter fits.

    if args.binary:
        cycles = 4
    else:
        cycles = 2

    optgroup1 = ['cont', 'twave', 'cont', 's',
                'cont', 'twave', 's', 'cont',
                'twave',
                'ip', 'v',
                'ip', 'v',
                'twave',  's',
                'twave',  's']

    optgroup2 = ['cont', 'twave', 'cont', 's1s2',
            'cont', 'twave', 's','s2', 'cont',
            'twave',
            'ip', 'v',
            'ip', 'v',
            'twave',  's','s2',
            'twave',  's','s1s2']

    optgroup = optgroup1.copy()

    nk = 1
    for nc, cycle in enumerate(np.arange(cycles), start=1):
        if cycle == 0:
            parstart = par_in.copy()
            dpars = dpars1
        else:
            dpars = dpars2

        for optkind in optgroup:
            parfit_1 = optimizer(parstart, dpars[optkind], hardbounds, fitobj, optimize)
            parstart = parfit_1.copy()
            if args.debug == True:
                outplotter_23(
                    parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(
                        order,night,tag,nk,optkind),
                        trk, inparam, args, step2or3, order)
                logger.debug(f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
            nk += 1

        if nc == 2:
            optgroup = optgroup2.copy()

    parfit = parfit_1.copy()

    #-------------------------------------------------------------------------------


    # if best fit stellar template power is very low, throw out result
    if parfit[1] < 0.1:
        logger.warning(f'  --> Stellar template power is low for {night}! '
                            'Data likely being misfit! Throwing out result...')
        return night, np.nan, np.nan

    if args.binary and parfit[25] < 0.05:
        logger.warning(f'  --> Secondary stellar template power is low for {night}! '
                            'Data likely being misfit! Throwing out result...')
        continue
    # if best fit stellar or telluric template powers are exactly equal
    # to their starting values, fit failed, throw out result
    if parfit[1] == par_in[1] or parfit[3] == par_in[3] or (args.binary and (parfit[25] == par_in[25])):
        logger.warning(f'  --> Stellar or telluric template powers have not '
                            f'budged from starting values for {night}! Fit is '
                            'broken! Optimizer bounds may be unfeasible, or '
                            'chi-squared may be NaN? Throwing out result...')
        return night, np.nan, np.nan

    # if best fit model dips below zero at any point, we're to close to edge of
    # blaze, fit may be compromised, throw out result
    smod,chisq,trash,trash2 = fmod(parfit,fitobj,args.binary)
    if len(smod[(smod < 0)]) > 0:
        logger.warning(f'  --> Best fit model dips below 0 for {night}! '
                            'May be too close to edge of blaze, throwing '
                            'out result...')
        continue


    #-------------------------------------------------------------------------------
    if args.plotfigs == True:
        parfitS1 = parfit.copy(); parfitS1[3] = 0; parfitS1[24] = 0;
        parfitS2 = parfit.copy(); parfitS2[3] = 0; parfitS2[1] = 0;
        parfitT = parfit.copy(); parfitT[1] = 0; parfitT[24] = 0;
        if args.binary:
            outplotter_23(
                parfitS1, fitobj, 'parfitS1_{}_{}_{}'.format(order,night,tag),
                trk, inparam, args, step2or3,order)
            outplotter_23(
                parfitS2, fitobj, 'parfitS2_{}_{}_{}'.format(order,night,tag),
                trk, inparam, args, step2or3,order)
        else:
            outplotter_23(
                parfitS1, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag),
                trk, inparam, args, step2or3,order)
        outplotter_23(
            parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag),
            trk, inparam, args, step2or3,order)
        outplotter_23(
            parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag),
            trk, inparam, args, step2or3,order)

    rv0 = parfit[0]
    # Barycentric correction
    rvsmini    = rv0 + inparam.bvcs[night+tag] \
                    + rv0*inparam.bvcs[night+tag]/(2.99792458e5**2)
    vsinismini = parfit[4]

    bestguess = np.round(rvsmini,5)
    vsinimini = np.round(vsinismini,5)

    if args.binary:
        rv2 = parfit[24]
        bestguess2 = rv2  + inparam.bvcs[night+tag] \
                    + rv2*inparam.bvcs[night+tag]/(2.99792458e5**2)
        vsinimini2 = parfit[26]
    else:
        bestguess2 = np.nan; vsinimini2 = np.nan;

    return night, bestguess, vsinimini, bestguess2, vsinimini2

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    args = _argparse_step2()
    inpath   = './Input/{}'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    #### Check user inputs

    initvsini = np.float(args.initvsini)
    vsinivary = np.float(args.vsinivary)

    if args.initvsini == '':
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI VALUE, "-i"')

    if (args.guesses == '') & (args.guessesX == ''):
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR RV VALUE(S) '
                    'BY USING "-g" OR "-gX"')

    if (args.temperature == '') & (args.logg == ''):
        sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
                    'STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO '
                    'SEE AVAILABLE TEMPLATES')
    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'./Engine/syn_template')

    if args.binary:
        if args.fluxratio == '':
            sys.exit('ERROR: YOU MUST PROVIDE A FLUX RATIO S2/S1, "-f"')
        if args.initvsini2 == '':
            sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI2 VALUE, "-i2"')
        if (args.temperature2 == '') & (args.logg2 == ''):
            sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
                        'SECONDARY STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE '
                        'AVAILABLE TEMPLATES')
        if args.template2.lower() not in ['synthetic', 'livingston', 'phoenix']:
            sys.exit('ERROR: UNEXPECTED SECONDARY STELLAR TEMPLATE FOR "-t" INPUT!')

        initvsini2 = float(args.initvsini2)
        vsinivary2 = float(args.vsinivary2)

        if args.template2.lower() == 'synthetic':
            #list of all syntheticstellar
            syntemp = [i for i in syntemp if i[:3] == 'syn']
            synT    = [ i.split('_')[2][1:]  for i in syntemp ]
            synlogg = [ i.split('_')[3][4:7] for i in syntemp ]
        elif args.template2.lower() == 'phoenix':
            #list of all phoenix
            syntemp = [i for i in syntemp if i[:3] == 'PHO']
            synT    = [ i.split('-')[1][4:]  for i in syntemp ]
            synlogg = [ i.split('-')[2][:3] for i in syntemp ]
        else:
            synT = [args.temperature2]; synlogg = [args.logg2]

        if args.temperature2 not in synT:
            sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp2" INPUT! '
                        f'{syntemp} AVALIABLE UNDER ./Engine/syn_template/')

        if args.logg2 not in synlogg:
            sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg2" INPUT! {syntemp} '
                        'AVALIABLE UNDER ./Engine/syn_template/')

    syntemp = os.listdir(f'./Engine/syn_template')
    syntemp = [i for i in syntemp if i[:3] == 'syn'] #list of all syntheticstellar

    synT    = [ i.split('_')[2][1:]  for i in syntemp ]
    synlogg = [ i.split('_')[3][4:7] for i in syntemp ]

    if args.temperature not in synT:
        sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp" INPUT! '
                    f'{syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    if args.logg not in synlogg:
        sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg" INPUT! {syntemp} '
                    'AVALIABLE UNDER ./Engine/syn_template/')

    #------------------------------

    if (args.guesses != '') & (args.guessesX != ''):
        sys.exit('ERROR: YOU CAN ONLY CHOOSE EITHER -g OR -gX')

    #------------------------------
    # Specify initial RV guesses as a single value applied to all nights

    if args.guesses != '':
        try:
            initguesses = float(args.guesses)
            initguesses_show = initguesses
        except:
            sys.exit('ERROR: -g ONLY TAKES A NUMBER AS INPUT!')

    #------------------------------

    # Load initial RV guesses from file
    if args.guessesX != '':
        try:
            guessdata = Table.read(
                f'./Output/{args.targname}_{args.band}/'
                    f'Initguesser_results_{args.guessesX}.csv',
                format='csv')

        except:
            sys.exit(
                f'ERROR: "./Output/{args.targname}_{args.band}/'
                    f'Initguesser_results_{args.guessesX}.csv" NOT FOUND!')

        initnights = np.array(guessdata['night'])
        initrvs    = np.array(guessdata['bestguess'])
        initguesses = {}
        initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
        for hrt in range(len(initnights)):
            initguesses[str(initnights[hrt])] = float(initrvs[hrt])
        if args.binary:
        initrvs2    = np.array(guessdata['bestguess2'])
        initguesses2 = {}
        for hrt in range(len(initnights)):
            initguesses2[str(initnights[hrt])] = float(initrvs2[hrt])

    #------------------------------
    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    if int(args.label_use) not in orders:
        sys.exit(
            f'Oops! -l_use INPUT "{args.label_use}" is not in "{orders}" '
                'from the given WRegion list!!')

#-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[37;1;41m {} band \033[0m
    WaveLength file     = \33[37;1;41m WaveRegions_{} \033[0m
    S/N cut             > \33[37;1;41m {} \033[0m
    Order Use           = \33[37;1;41m Order {} \033[0m
    Initial vsini       = \33[37;1;41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[37;1;41m {} km/s \033[0m
    RV initial guess    = \33[37;1;41m {} \033[0m
    Stellar template use= \33[37;1;41m {} \033[0m
    syn template temp   = \33[37;1;41m {} \033[0m
    syn template logg   = \33[37;1;41m {} \033[0m
    Threads use         =  {}
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut,
                args.label_use, initvsini, vsinivary, initguesses_show,
                args.template, args.temperature, args.logg, args.Nthreads))

# Target Spectral Type= \33[41m {} \033[0m             <-------  [late K, M] recommended 'synthetic', [F, G, early K] SpTy recommended 'livingston'
    if not args.skip:
        while True:
            inpp = input("Press [Y]es to continue, [N]o to quit...\n --> ")
            if 'n' in inpp.lower():
                sys.exit('QUIT, PLEASE RE-ENTER YOUR PARAMETERS')
            elif 'y' in inpp.lower():
                break
            else:
                print('I cannot understand what you are saying... TRY AGAIN')
                continue

    if args.binary:
        print(u'''
        PLUS BINARY PARAMETERS:
        Initial vsini #2      = \33[37;1;41m {} km/s \033[0m
        vsini #2 vary range    \u00B1 \33[37;1;41m {} km/s \033[0m
        Stellar template #2 use= \33[37;1;41m {} \033[0m
        syn template temp #2  = \33[37;1;41m {} \033[0m
        syn template logg #2   = \33[37;1;41m {} \033[0m
        syn template B #2   = \33[37;1;41m {} \033[0m
        '''.format(initvsini2, vsinivary2, args.template2,
           args.temperature2, args.logg2, args.B2))
        if not args.skip:
            while True:
                inpp = input("Press [Y]es to continue, [N]o to quite...\n --> ")
                if 'n' in inpp.lower():
                    sys.exit('QUIT, PLEASE RE-ENTER YOUR PARAMETERS')
                elif 'y' in inpp.lower():
                    break
                else:
                    print('I cannot understand what you are saying... TRY AGAIN')
                    continue

    print('---------------------------------------------------------------')
    print('Running Step 2 for {}...'.format(args.targname))
    print('This Will Take a While..........')

    #-------------------------------------------------------------------------------

    # Make output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')
    filesndirs = os.listdir(f'./Output/{args.targname}_{args.band}')

    trk = 1; go = True
    while go == True:
        iniguess_dir = 'Initguesser_results_{}.csv'.format(trk)
        if iniguess_dir not in filesndirs:
            break
        trk += 1

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs')

    step2or3 = 2
    temp_f_dir = f'./Output/{args.targname}_{args.band}/figs/'\
                        f'main_step{step2or3}_{args.band}_{trk}'
    if not os.path.isdir(temp_f_dir):
        os.mkdir(temp_f_dir)

    outpath = f'./Output/{args.targname}_{args.band}'


    #-------------------------------------------------------------------------------

    # Set up logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(
        f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    #-------------------------------------------------------------------------------
    # Create output file to write to
    logger.info(
        f'Writing output to ./Output/{args.targname}_{args.band}/{iniguess_dir}')

    filew = open(
        f'./Output/{args.targname}_{args.band}/{iniguess_dir}','w')
    filew.write('night, bestguess, vsini')
    filew.write('\n')

    #-------------------------------------------------------------------------------


    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(
                    'NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(
                        nnn, args.targname
                        ))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))


    logger.info('Analyze with {} nights'.format(len(nightsFinal)))

    intnights = np.array([int(i[:8]) for i in nightsFinal])
    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these nights
                slightly differently. When you run Step 3, RVs will be output in
                two formats: one with the defocus nights separated, and the other
                with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as the result
                of the defocus. IGRINS RV will process defocus nights the same way
                as the others, but when you run Step 3, will still output the results
                in two formats like it does with the K band.''')

    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(
        logger, args.template, args.band, int(args.temperature),
        float(args.logg), float(args.B)
        )

    # Save pars in class for future use
    inparam = InParams(inpath, outpath, initvsini, vsinivary, args.plotfigs,
                       initguesses, bvcs, tagsA, tagsB, nightsFinal, mwave0,
                       mflux0, None, xbounddict, maskdict)

    if args.binary:
        print('\n Loading secondary stellar template... \n')
        watm,satm, mwave2, mflux2 = setup_templates(
            logger, args.template2, args.band, np.int(args.temperature2),
            np.float(args.logg2), np.float(args.B2)
            )
        inparam.addsecondary(initvsini2,vsinivary2,mwave2,mflux2,initguesses2)

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    func = partial(main, args, inparam, orders,
                        int(args.label_use), trk, step2or3 )
    outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

    # Write outputs to file
    vsinis = []; finalrvs = []; vsinis2 = []; finalrvs2 = [];
    for n in range(len(nightsFinal)):
        nightout = outs[n]
        if args.binary:
            filew.write('{}, {}, {}'.format(nightout[0], nightout[1], nightout[2], nightout[3], nightout[4]))
        else:
            filew.write('{}, {}, {}'.format(nightout[0], nightout[1], nightout[2]))
        filew.write('\n')
        vsinis2.append(nightout[4])
        finalrvs2.append(nightout[3])
        vsinis.append(nightout[2])
        finalrvs.append(nightout[1])

    filew.close()

    warning_r = log_warning_id(f'{outpath}/{args.targname}_{args.band}.log', start_time)
    if warning_r:
        print(f'''
**********************************************************************************
WARNING!! you got warning message during this run. Please check the log file under:
          {outpath}/{args.targname}_{args.band}.log
**********************************************************************************
''')
    print('\n')
    if not args.debug: logger.addHandler(stream_hander)
    print('--------!Initial Guess!--------')
    logger.info('RV results:    mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(finalrvs),
                                                                                                    np.nanmedian(finalrvs),
                                                                                                    np.nanstd(finalrvs)      ))
    logger.info('vsini results: mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(vsinis),
                                                                                                    np.nanmedian(vsinis),
                                                                                                    np.nanstd(vsinis)      ))
    end_time = datetime.now()
    logger.info('RV 2 results:    mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(finalrvs2),
                                                                                                    np.nanmedian(finalrvs2),
                                                                                                    np.nanstd(finalrvs2)      ))
    logger.info('vsini 2 results: mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(vsinis2),
                                                                                                    np.nanmedian(vsinis2),
                                                                                                    np.nanstd(vsinis2)    ))
    logger.info('RV Initial Guess DONE... Duration: {}'.format(end_time - start_time))
    logger.info(f'Output saved under ./Output/{args.targname}_{args.band}/{iniguess_dir}')
    print('---------------------------------------------------------------')
    print('You can now try to get a better RV initial guess with by rerunning Step 2 with -gX set to the run number you just completed.')
    print('OR, you can go on to the full RV analysis in Step 3.')
    print('####################################################################################')
