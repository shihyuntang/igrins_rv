from Engine.importmodule import *
from Engine.importmodule import read_prepdata

from Engine.IO_AB      import setup_templates, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs, inparams
from Engine.macbro     import macbro
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod
from Engine.outplotter import outplotter_23
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def ini_MPinst(args, inparam, orders, order_use, trk, step2or3, i):

    # Main function for RV fitting that will be threaded over by multiprocessing

    nights   = inparam.nights
    night    = nights[i] # current looped night

    order   = order_use
    xbounds = inparam.xbounddict[order]
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

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    # Load synthetic telluric template generated during Step 1
    # [:8] here is to ensure program works under Night_Split mode
    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_treated_{args.band}.fits'
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

    # Check whether Telfit hit critical error in Step 1 for the chosen order with this night. If so, try another order. If all hit the error, skip the night.
    nexto = 0
    ordertry = order
    while 1 == 1:
        fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == ordertry ][0]

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
    a0contx   = a0contx[(continuum != 0)]
    continuum = continuum[(continuum != 0)]

    # Use instrumental profile FWHM dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]
#-------------------------------------------------------------------------------
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,                                                # 0: The shift of the stellar template (km/s)
                      0.3,                                                   # 1: The scale factor for the stellar template
                      0.0,                                                   # 2: The shift of the telluric template (km/s)
                      0.6,                                                   # 3: The scale factor for the telluric template
                      inparam.initvsini,                                     # 4: vsini (km/s)
                      IPpars[2],                                             # 5: The instrumental resolution (FWHM) in pixels
                      np.nan,                                                # 6: Wavelength 0-pt
                      np.nan,                                                # 7: Wavelength linear component
                      np.nan,                                                # 8: Wavelength quadratic component
                      np.nan,                                                # 9: Wavelength cubic component
                      1.0,                                                   #10: Continuum zero point
                      0.,                                                    #11: Continuum linear component
                      0.,                                                    #12: Continuum quadratic component
                      IPpars[1],                                             #13: Instrumental resolution linear component
                      IPpars[0]])                                            #14: Instrumental resolution quadratic component

    # rvsmini = []; vsinismini = [];

    # # Iterate over all A/B exposures
    # for t in np.arange(len(tagsnight)):
    tag = tagsnight[0]
    beam = beamsnight[0]

    # Retrieve pixel bounds for where within each other significant telluric absorption is present.
    # If these bounds were not applied, analyzing some orders would give garbage fits.
    if args.band=='K':
        if int(order) in [14]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]

    elif args.band=='H':
        if int(order) in [13, 14, 16, 20]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]

    # Load target spectrum
    x,wave,s,u = init_fitsread(f'{inparam.inpath}/',
                                'target',
                                'combined',
                                night,
                                order,
                                tag,
                                args.band,
                                bound_cut)

    #-------------------------------------------------------------------------------

    # Execute S/N cut
    s2n = s/u
    if np.nanmedian(s2n) < float(args.SN_cut):
        logger.warning('  --> Bad S/N {:1.3f} < {} for {}{} {}... '.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
        pass

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    nzones = 5
    x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones);
    u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
    x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones);
    u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

    # Cut spectrum to within wavelength regions defined in input list
    s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
    x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

    # Trim stellar template to relevant wavelength range
    mwave_in,mflux_in = stellarmodel_setup(wave_piece, inparam.mwave0, inparam.mflux0)

    # Trim telluric template to relevant wavelength range
    satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
    watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

    # Make sure data is within telluric template range (shouldn't do anything)
    s_piece    = s_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    u_piece    = u_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    x_piece    = x_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

    #-------------------------------------------------------------------------------

    par = pars0.copy()

    # Get initial guess for cubic wavelength solution from reduction pipeline
    f = np.polyfit(x_piece,wave_piece,3)
    par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
    par[9] = par9in ;  par[8] = par8in ;  par[7] = par7in ;  par[6] = par6in

    par[0] = initguesses-inparam.bvcs[night+tag] # Initial RV with barycentric correction

    # Arrays defining parameter variations during optimization steps.
    # Optimization will cycle twice. In the first cycle, the RVs can vary more than in the second.
    dpars1 = {'cont' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
              'wave' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0]),
              't'    : np.array([ 0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              'ip'   : np.array([ 0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              's'    : np.array([20.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              'v'    : np.array([ 0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}
    dpars2 = {'cont' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
              'wave' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0]),
              't'    : np.array([ 0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              'ip'   : np.array([ 0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              's'    : np.array([ 5.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
              'v'    : np.array([ 0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}

    continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
    s_piece /= np.median(s_piece)
    fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]))

    #-------------------------------------------------------------------------------
    # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
    optimize = True
    par_in = par.copy()
    hardbounds = [par_in[4]-dpars1['v'][4],  par_in[4]+dpars1['v'][4],
                  par_in[5]-dpars1['ip'][5], par_in[5]+dpars1['ip'][5]]
    if hardbounds[0] < 0:
        hardbounds[0] = 0
    if hardbounds[3] < 0:
        hardbounds[3] = 1

    # Begin optimization. Fit the blaze, the wavelength solution, the telluric template power and RV, the stellar template power and RV, the
    # zero point for the instrumental resolution, and the vsini of the star separately, iterating and cycling between each set of parameter fits.
    cycles = 2

    optgroup = ['cont', 'wave', 't', 'cont', 's',
                'cont', 'wave', 't', 's', 'cont',
                'wave',
                'ip', 'v',
                'ip', 'v',
                't',  's',
                't',  's']

    nk = 1
    for nc, cycle in enumerate(np.arange(cycles), start=1):
        if cycle == 0:
            parstart = par_in.copy()
            dpars = dpars1
        else:
            dpars = dpars2

        for optkind in optgroup:
            parfit_1 = optimizer( parstart, dpars[optkind], hardbounds, fitobj, optimize)
            parstart = parfit_1.copy()
            if args.debug:
                outplotter_23(parfit_1, fitobj, '{}_{}_{}_parfit_{}{}'.format(order,night,tag,nk,optkind), trk, inparam, args, step2or3)
                logger.debug(f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
            nk += 1

    parfit = parfit_1.copy()

    #-------------------------------------------------------------------------------

    # if best fit stellar template power is very low, throw out result
    if parfit[1] < 0.1:
        logger.warning(f'  --> parfit[1] < 0.1, {night} parfit={parfit}...')
        pass

    # if best fit stellar or telluric template powers are exactly equal to their starting values, optimization failed, throw out result
    if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
        logger.warning(f'  --> parfit[1] == par_in[1] or parfit[3] == par_in[3], {night}...')
        pass

    # if best fit model dips below zero at any point, we're too close to edge of blaze, fit may be comrpomised, throw out result
    smod,chisq = fmod(parfit,fitobj)
    if len(smod[(smod < 0)]) > 0:
        logger.warning(f'  --> len(smod[(smod < 0)]) > 0, {night}...')
        pass

#-------------------------------------------------------------------------------
    if args.plotfigs:
        parfitS = parfit.copy(); parfitS[3] = 0
        parfitT = parfit.copy(); parfitT[1] = 0
        outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
        outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
        outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)

    rv0 = parfit[0] - parfit[2]  # Correct for RV of the atmosphere, since we're using that as the basis for the wavelength scale

    rvsmini    = rv0 + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2) # Barycentric correction
    vsinismini = parfit[4]

    bestguess = round(rvsmini,5)
    vsinimini = round(vsinismini,5)
    return night, bestguess, vsinimini

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 2',
                                     description = '''
                                     Only required if the average RV of the target star is unknown to $>$ 5 \kms precision. \n
                                     Performs an abbreviated analysis of the target star observations in order to converge to coarsely accurate RVs,
                                     which will be used as starting points for the more precise analysis in the next step;
                                     simultaneously does the same for target star's \vsini.  \n
                                     Only a single order is used - by default, the first one in the list of wavelength regions, but the user can specify otherwise.  \n
                                     All separate exposures for a given observation are combined into one higher S/N spectrum before fitting occurs.
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
                        help="Stellar template. Pick from 'synthetic', 'livingston', or 'user_defined'. Default = 'synthetic'",
                        type=str,   default='synthetic' )
    parser.add_argument('-Temp',      dest="temperature",           action="store",
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
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots under ./Temp/Debug/*target_*band/main_step2")
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.85')
    args = parser.parse_args()
    inpath   = './Input/{}'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    #### Check user inputs

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)

    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'user_defined']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    # if args.template.lower() in ['synthetic', 'livingston']:
    #     if args.sptype not in ['F','G','K','M']:
    #         sys.exit('ERROR: DEFAULT TEMPLATES ONLY COVER F G K M STARS!')

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
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders = read_prepdata(args)

    if int(args.label_use) not in orders:
        sys.exit(f'Oops! -l_use INPUT "{args.label_use}" is not in "{orders}" from the given WRegion list!!')

#-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[41m {} band \033[0m
    WaveLength file     = \33[41m WaveRegions_{} \033[0m
    S/N cut             > \33[41m {} \033[0m
    Order Use           = \33[41m Order {} \033[0m
    Initial vsini       = \33[41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[41m {} km/s \033[0m
    RV initial guess    = \33[41m {} \033[0m
    Stellar template use= \33[41m {} \033[0m
    syn template temp   = \33[41m {} \033[0m
    syn template logg   = \33[41m {} \033[0m
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

    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, flot(args.temperature), flot(args.logg))

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)

    #-------------------------------------------------------------------------------

    # Run order by order, multiprocessing over nights within an order
    pool = mp.Pool(processes = args.Nthreads)
    func = partial(ini_MPinst, args, inparam, orders, int(args.label_use), trk, step2or3 )
    outs = pool.map(func, np.arange(len(nightsFinal)))
    pool.close()
    pool.join()

    # Write outputs to file
    vsinis = []; finalrvs = [];
    for n in range(len(nightsFinal)):
        nightout = outs[n]
        filew.write('{}, {}, {}'.format(nightout[0], nightout[1], nightout[2]))
        filew.write('\n')
        vsinis.append(nightout[2])
        finalrvs.append(nightout[1])

    filew.close()

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
