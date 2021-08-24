from Engine.importmodule import *
from Engine.importmodule import read_prepdata

from Engine.IO_AB      import setup_templates, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs,inparams
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod, fmod_conti
# from Engine.opt_rebintel        import optimizer, fmod, fmod_conti
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.crmask    import CRmasker
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def rv_MPinst(args, inparam, orders, order_use, trk, step2or3, i):

    # Main function for RV fitting that will be threaded over by multiprocessing

    nights = inparam.nights
    night  = nights[i] # current looped night

    order   = orders[order_use]
    xbounds_org = inparam.xbounddict[order]

    if args.debug:
        print('Working on order {:02d}, night {:03d}/{:03d} ({}) PID:{}...'.format(int(order),
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
        sys.exit('ERROR! EXPECING SINGAL NUMBER OR FILE FOR INITGUESSES! QUITTING!')

    if np.isnan(initguesses) == True:
        logger.warning(f'  --> Previous run of {night} found it inadequate, skipping...')
        return night, np.nan, np.nan

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = np.array([]);
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight = np.append(beamsnight, 'B')

    # Only do B exposures, and just use first B nodding
    masterbeam = 'B'; beam = 'B';
    try:
        tag  = tagsnight[0]
    except IndexError:
        logger.warning(f'  --> No B nodding(frame) for night {night}, skipping...')
        return night, np.nan, np.nan


    #-------------------------------------------------------------------------------
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,                                                # 0: The shift of the stellar template (km/s) [assigned later]
                      0.3,                                                   # 1: The scale factor for the stellar template
                      0.0,                                                   # 2: The shift of the telluric template (km/s)
                      0.6,                                                   # 3: The scale factor for the telluric template
                      inparam.initvsini,                                     # 4: vsini (km/s)
                      np.nan,                                                # 5: The instrumental resolution (FWHM) in pixels
                      np.nan,                                                # 6: Wavelength 0-pt
                      np.nan,                                                # 7: Wavelength linear component
                      np.nan,                                                # 8: Wavelength quadratic component
                      np.nan,                                                # 9: Wavelength cubic component
                      1.0,                                                   #10: Continuum zero point
                      0.,                                                    #11: Continuum linear component
                      0.,                                                    #12: Continuum quadratic component
                      np.nan,                                                #13: Instrumental resolution linear component
                      np.nan,                                                #14: Instrumental resolution quadratic component
                      0,                                                     #15: Blaze dip center location
                      0,                                                     #16: Blaze dip full width
                      0,                                                     #17: Blaze dip depth
                      0,                                                     #18: Secondary blaze dip full width
                      0,                                                     #19: Blaze dip depth
                      0.0,                                                   #20: Continuum cubic component
                      0.0,                                                   #21: Continuum quartic component
                      0.0,                                                   #22: Continuum pentic component
                      0.0])                                                  #23: Continuum hexic component


    # This one specific order is small and telluric dominated, start with greater stellar template power to ensure good fits
    if int(order) == 13:
        pars0[1] = 0.8

    # Load synthetic telluric template generated during Step 1
    # [:8] here is to ensure program works under Night_Split mode

    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

    try:
        hdulist = fits.open(A0loc)
    except IOError:
        logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
        return night, np.nan, np.nan

    # Find corresponding table in fits file, given the tables do not go sequentially by order number due to multiprocessing in Step 1
    num_orders = 0
    for i in range(25):
        try:
            hdulist[i].columns[0].name[9:]
            num_orders += 1
        except:
            continue

    fits_layer = [ i for i in np.arange(num_orders)+1 if np.int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

    # Check whether Telfit hit critical error in Step 1 for the chosen order with this night. If so, try another order. If all hit the error, skip the night.
    nexto = 0
    ordertry = order
    while 1 == 1:
        fits_layer = [ i for i in np.arange(num_orders)+1 if np.int(hdulist[i].columns[0].name[9:]) == ordertry ][0]

        tbdata = hdulist[ fits_layer ].data
        flag = np.array(tbdata[f'ERRORFLAG{ordertry}'])[0]

        if flag == 1:  # If Telfit hit unknown critical error in Step 1, this order can't be used for this night. Try another.
            orderbad = ordertry
            ordertry = orders[nexto]
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {orderbad} NIGHT: {night}, TRYING ORDER {ordertry} INSTEAD...')

        else: # All good, continue
            order = ordertry
            break

        nexto += 1
        if nexto == len(orders):
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ALL ORDERS FOR NIGHT: {night}, skipping...')
            return night, np.nan, np.nan


    watm = tbdata['WATM'+str(order)]
    satm = tbdata['SATM'+str(order)]
    a0contx    = tbdata['X'+str(order)]
    continuum  = tbdata['BLAZE'+str(order)]

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]
    satm[(satm < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
    a0contx = a0contx[(continuum != 0)]
    continuum = continuum[(continuum != 0)]


    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if (np.int(night[:8]) < 20180401) or (np.int(night[:8]) > 20190531):
        IPpars = inparam.ips_tightmount_pars[args.band][masterbeam][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][masterbeam][order]

    # Retrieve pixel bounds for where within each other significant telluric absorption is present.
    # If these bounds were not applied, analyzing some orders would give garbage fits.
    if args.band=='K':
        if int(order) in [3, 4, 13, 14]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]

    elif args.band=='H':
        if int(order) in [6, 10, 11, 13, 14, 16, 17, 20, 21, 22]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]

    # Load target spectrum
    x,wave,s,u = init_fitsread(f'{inparam.inpath}/',
                                'target',
                                'combined'+str(masterbeam),
                                night,
                                order,
                                inparam.tagsB[night][0],
                                args.band,
                                bound_cut)

    #-------------------------------------------------------------------------------

    # Execute S/N cut
    s2n = s/u
    if np.nanmedian(s2n) < np.float(args.SN_cut):
        logger.warning('  --> Bad S/N {:1.3f} < {} for {}{} {}... '.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
        pass

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    nzones = 5
    x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
    x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

    # --- xbound re-set -- sytang add ---
    chi2_box = []; par_box = []; xbound_box = []

    for adj_xbound in [0, 10, 20, 30, -10, -20, -30]:
        xbounds = [xbounds_org[0]+adj_xbound, xbounds_org[-1]+adj_xbound]

        # -- check if xbounds is within the bound_cut --
        bound_cut_pixel = [bound_cut[0], 2048-bound_cut[-1]]

        if xbounds[0] < bound_cut_pixel[0]+5:
            xbounds[0] = bound_cut_pixel[0]+5

        if xbounds[-1] > bound_cut_pixel[-1]-5:
            xbounds[-1] = bound_cut_pixel[-1]-5
        # --------------

        # Cut spectrum to within wavelength regions defined in input list
        s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
        x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

        # Trim stellar template to relevant wavelength range
        mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

        # Trim telluric template to relevant wavelength range
        satm_in = satm[(watm > np.min(wave_piece)*1e4 - 11) & (watm < np.max(wave_piece)*1e4 + 11)]
        watm_in = watm[(watm > np.min(wave_piece)*1e4 - 11) & (watm < np.max(wave_piece)*1e4 + 11)]

        # Make sure data is within telluric template range (shouldn't do anything)
        s_piece    = s_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        u_piece    = u_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        x_piece    = x_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]

        # Normalize continuum from A0 to flux scale of data
        continuum /= np.nanmedian(continuum)
        continuum *= np.nanpercentile(s_piece,99)

        # --------------------------------------------------------------

        par = pars0.copy()

        # Get initial guess for cubic wavelength solution from reduction pipeline
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = initguesses-inparam.bvcs[night+tag] # Initial RV with barycentric correction
        par[5] = IPpars[2]; par[13] = IPpars[1]; par[14] = IPpars[0];

        # Arrays defining parameter variations during optimization steps
        # Optimization will cycle twice. In the first cycle, the RVs can vary more than in the second.
        #                             | 0    1    2    3 |  | ------ 4 ------ |  | 5 |   | 6     7     8           9  |  |10  11  12| |13 14|  |15  16  17  18  19|  |20   21   22   23 |
        dpars1 = {'cont' : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    1e7, 1, 1,   0, 0,    0., 0., 0., 0., 0.,   1.0, 1.0, 1.0, 1.0  ]),
                  'twave': np.array([  0.0, 0.0, 0.0, 1.0,   0.0,                 0.0,   10.0, 10.0,  5.00000e-5, 1e-7,   0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  'ip'   : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.5,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  's'    : np.array([ 20.0, 2.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  'v'    : np.array([  0.0, 0.0, 0.0, 0.0,   inparam.vsinivary,   0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ])}

        dpars2 = {'cont' : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    1e7, 1, 1,   0, 0,    0., 0., 0., 0., 0.,  1.0, 1.0, 1.0, 1.0  ]),
                  'twave': np.array([  0.0, 0.0, 0.0, 1.0,   0.0,                 0.0,   10.0, 10.0,  5.00000e-5, 1e-7,   0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  'ip'   : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.5,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  's'    : np.array([  5.0, 2.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ]),
                  'v'    : np.array([  0.0, 0.0, 0.0, 0.0,   inparam.vsinivary,   0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,    0., 0., 0., 0., 0.,   0.0, 0.0, 0.0, 0.0 ])}

        # Use quadratic blaze correction for order 13; cubic for orders 6, 14, 21; quartic for orders 16 and 22
        if args.band == 'H':
            if int(order) in [13]:
                dpars1['cont'][20] = 0.; dpars1['cont'][21] = 0.; dpars1['cont'][22] = 0.; dpars1['cont'][23] = 0.;
                dpars2['cont'][20] = 0.; dpars2['cont'][21] = 0.; dpars2['cont'][22] = 0.; dpars2['cont'][23] = 0.;
            elif int(order) in [6,14,21]:
                dpars1['cont'][21] = 0.; dpars1['cont'][22] = 0.; dpars1['cont'][23] = 0.;
                dpars2['cont'][21] = 0.; dpars2['cont'][22] = 0.; dpars2['cont'][23] = 0.;
            else:
                pass
        else:
            if int(order) in [3,4,5]:
                dpars1['cont'][21] = 0.; dpars1['cont'][22] = 0.; dpars1['cont'][23] = 0.;
                dpars2['cont'][21] = 0.; dpars2['cont'][22] = 0.; dpars2['cont'][23] = 0.;
            elif int(order) in [6]:
                dpars1['cont'][22] = 0.; dpars1['cont'][23] = 0.;
                dpars2['cont'][22] = 0.; dpars2['cont'][23] = 0.;
            else:
                pass

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]),masterbeam,[np.array([],dtype=int),np.array([],dtype=int)])

        #-------------------------------------------------------------------------------

        # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
        optimize = True
        par_in = par.copy()
        hardbounds = [par_in[4] - dpars1['v'][4],  par_in[4] + dpars1['v'][4],
                      par_in[5] - dpars1['ip'][5], par_in[5] + dpars1['ip'][5]
                     ]

        if hardbounds[0] < 0.5:
            hardbounds[0] = 0.5
        if hardbounds[2] < 1:
            hardbounds[2] = 1

        # Begin optimization. Fit the blaze, the wavelength solution, the telluric template power and RV, the stellar template power and RV, the
        # zero point for the instrumental resolution, and the vsini of the star separately, iterating and cycling between each set of parameter fits.

        cycles = 2

        optgroup = ['cont', 'twave', 'cont', 's',
                    'cont', 'twave', 's', 'cont',
                    'twave',
                    'ip', 'v',
                    'ip', 'v',
                    'twave',  's',
                    'twave',  's']

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
                    outplotter_23(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}_xbound_{}-{}'.format(order,night,tag,nk,optkind,xbounds[0],xbounds[-1]), trk, inparam, args, step2or3, order)
                    logger.debug(f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
                nk += 1

        parfit = parfit_1.copy()
        smod,chisq = fmod(parfit,fitobj)

        if args.plotfigs == True:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3, order)
            outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3, order)
            outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}_xbound_{}-{}_chi2-{:1.2f}'.format(order,night,tag,xbounds[0],xbounds[-1],chisq), trk, inparam, args, step2or3, order)

        chi2_box.append(chisq)
        par_box.append(parfit)
        xbound_box.append(xbounds)

    min_chi2_par_idx = np.argmin(chi2_box)
    parfit   = par_box[min_chi2_par_idx]
    xbounds0 = xbound_box[min_chi2_par_idx][0]
    xbounds1 = xbound_box[min_chi2_par_idx][-1]



    #-------------------------------------------------------------------------------


    # if best fit stellar template power is very low, throw out result
    if parfit[1] < 0.1:
        logger.warning(f'  --> Stellar template power is low for {night}! Data likely being misfit! Throwing out result...')
        return night, np.nan, np.nan

    # if best fit stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
    if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
        logger.warning(f'  --> Stellar or telluric template powers have not budged from starting values for {night}! Fit is broken! Optimizer bounds may be unfeasible, or chi-squared may be NaN? Throwing out result...')
        return night, np.nan, np.nan

    # if best fit model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
    smod,chisq = fmod(parfit,fitobj)
    if len(smod[(smod < 0)]) > 0:
        logger.warning(f'  --> Best fit model dips below 0 for {night}! May be too close to edge of blaze, throwing out result...')
        return night, np.nan, np.nan


    #-------------------------------------------------------------------------------
    rv0 = parfit[0]
    rvsmini    = rv0 + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(2.99792458e5**2) # Barycentric correction
    vsinismini = parfit[4]

    bestguess = np.round(rvsmini,5)
    vsinimini = np.round(vsinismini,5)
    return night, bestguess, vsinimini, xbounds0, xbounds1

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 2',
                                     description = '''
                                     Required if the average RV of the target star is unknown to > 5 km/s precision. \n
                                     Performs an abbreviated analysis of the target star observations in order to converge to coarsely accurate RVs, which will be used as starting points for the more precise analysis in the next step; \n
                                     simultaneously does the same for target star's vsini. \n
                                     Only the single most precise wavelength region is used, and all separate observations for a given exposure are combined into one higher S/N spectrum before being fit.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name, no space",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))
    parser.add_argument("-l_use",   dest="label_use",        action="store",
                        help="Specify ORDER used. Default is the first in WRegion list",
                        type=int,   default=int(0))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 ",
                        type=str,   default='50')
    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s)",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",        action="store",
                        help="Range of allowed vsini variation during optimization (float, if set to zero vsini will be held constant), default = 5.0 km/s",
                        type=str,   default='5' )
    parser.add_argument('-g',       dest="guesses",          action="store",
                        help="Initial guess for RV (int or float, km/s) for all nights. Use -gX instead if you want to reference an Initguesser_results file from a previous run of this step, which will have a different initial guess for each night",
                        type=str,   default='')
    parser.add_argument('-gX',      dest="guessesX",         action="store",
                        help="The number, X, that refers to the ./*targname/Initguesser_results_X file you wish to use for initial RV guesses",
                        type=str,   default='')
    parser.add_argument('-t',       dest="template",         action="store",
                        help="Stellar template. Pick from 'synthetic', 'PHOENIX', or 'livingston'. Default = 'synthetic'",
                        type=str,   default='synthetic' )
    parser.add_argument('-temp',      dest="temperature",           action="store",
                        help="The synthetic template temperature used, e.g., 5000",
                        type=str,   default='' )
    parser.add_argument('-logg',      dest="logg",           action="store",
                        help="The synthetic template logg used, e.g., 4.5",
                        type=str,   default='' )
    # parser.add_argument('-sp',      dest="sptype",           action="store",
    #                     help="The spectral type of the *target. (Letter only)",
    #                     type=str,   default='' )
    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If set, will generate plots of the basic fitting results under ./Output/*targname_*band/figs/main_step2_*band_*runnumber")
    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots.")
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 1.0.0')
    args = parser.parse_args()
    inpath   = './Input/{}'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    #### Check user inputs

    initvsini = np.float(args.initvsini)
    vsinivary = np.float(args.vsinivary)

    if args.initvsini == '':
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI VALUE, "-i"')

    if (args.guesses == '') & (args.guessesX == ''):
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR RV VALUE(S) BY USING "-g" OR "-gX"')

    if (args.temperature == '') & (args.logg == ''):
        sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE AVAILABLE TEMPLATES')
    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'./Engine/syn_template')
    syntemp = [i for i in syntemp if i[:3] == 'syn'] #list of all syntheticstellar

    synT    = [ i.split('_')[2][1:]  for i in syntemp ]
    synlogg = [ i.split('_')[3][4:7] for i in syntemp ]

    if args.temperature not in synT:
        sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp" INPUT! {syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    if args.logg not in synlogg:
        sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg" INPUT! {syntemp} AVALIABLE UNDER ./Engine/syn_template/')

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
            guessdata = Table.read(f'./Output/{args.targname}_{args.band}/Initguesser_results_{args.guessesX}.csv', format='csv')

        except:
            sys.exit(f'ERROR: "./Output/{args.targname}_{args.band}/Initguesser_results_{args.guessesX}.csv" NOT FOUND!')

        initnights = np.array(guessdata['night'])
        initrvs    = np.array(guessdata['bestguess'])
        initguesses = {}
        initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
        for hrt in range(len(initnights)):
            initguesses[str(initnights[hrt])] = float(initrvs[hrt])

    #------------------------------
    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    if int(args.label_use) not in orders:
        sys.exit(f'Oops! -l_use INPUT "{args.label_use}" is not in "{orders}" from the given WRegion list!!')

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
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut, args.label_use,
               initvsini, vsinivary, initguesses_show, args.template, args.temperature, args.logg, args.Nthreads))

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

    trk = 1; go = True;
    while go == True:
        iniguess_dir = 'Initguesser_results_{}.csv'.format(trk)
        if iniguess_dir not in filesndirs:
            break
        trk += 1

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs')

    step2or3 = 2
    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs/main_step{step2or3}_{args.band}_{trk}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs/main_step{step2or3}_{args.band}_{trk}')

    outpath = f'./Output/{args.targname}_{args.band}'


    #-------------------------------------------------------------------------------

    # Set up logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    #-------------------------------------------------------------------------------
    # Create output file to write to
    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{iniguess_dir}')

    filew = open(f'./Output/{args.targname}_{args.band}/{iniguess_dir}','w')
    filew.write('night, bestguess, vsini')
    filew.write('\n')

    #-------------------------------------------------------------------------------


    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
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
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, int(args.temperature), float(args.logg))

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):

        func = partial(rv_MPinst, args, inparam, orders, jerp, trk, step2or3 )
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

        # -- advance save with xbound and orders --
        order = orders[jerp]
        # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
        for i in range(len(nightsFinal)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvbox     = outsbox[1]
                vsinibox  = outsbox[2]
                xbound0   = outsbox[3]
                xbound1   = outsbox[4]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvbox     = np.concatenate((rvbox,outsbox[1]))
                vsinibox  = np.concatenate((vsinibox,outsbox[2]))
                xbound0   = np.concatenate((xbound0,outsbox[3]))
                xbound1   = np.concatenate((xbound1,outsbox[4]))

        nightsbox = np.array(nightsbox)
        vsinitags = []

        # Save results to fits file
        c1    = fits.Column(name='NIGHT'+str(order),  array=nightsbox, format='{}A'.format(len(nights[0])) )
        c2    = fits.Column(name='RV'+str(order),     array=rvbox,     format='D')
        c3    = fits.Column(name='VSINI'+str(order),  array=vsinibox,  format='D')
        c4    = fits.Column(name='xBOUND0'+str(order),array=xbound0,    format='D')
        c5    = fits.Column(name='xBOUND1'+str(order),array=xbound1,    format='D')
        cols  = fits.ColDefs([c1,c2,c3,c4,c5])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if jerp == 0: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}'.format(inparam.outpath, iniguess_dir.replace('.csv', '.fits')))
        else:
            hh = fits.open('{}/{}'.format(inparam.outpath, iniguess_dir.replace('.csv', '.fits')))
            hh.append(hdu_1)
            hh.writeto('{}/{}'.format(inparam.outpath, iniguess_dir.replace('.csv', '.fits')),overwrite=True)


        # -- normal save rv and vsini --
        if int(args.label_use) == int(orders[jerp]):
            # Write outputs to file
            vsinis = []; finalrvs = [];
            for n in range(len(nightsFinal)):
                nightout = outs[n]
                filew.write('{}, {}, {}'.format(nightout[0], nightout[1], nightout[2]))
                filew.write('\n')
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
    logger.info('RV Initial Guess DONE... Duration: {}'.format(end_time - start_time))
    logger.info(f'Output saved under ./Output/{args.targname}_{args.band}/{iniguess_dir}')
    print('---------------------------------------------------------------')
    print('You can now try to get a better RV initial guess with by rerunning Step 2 with -gX set to the run number you just completed.')
    print('OR, you can go on to the full RV analysis in Step 3.')
    print('####################################################################################')
