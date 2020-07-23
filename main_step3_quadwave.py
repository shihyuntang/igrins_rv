from Engine.importmodule import *
from Engine.importmodule import read_prepdata

from Engine.IO_AB      import setup_templates, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs,inparams
from Engine.macbro     import macbro
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod, fmod_conti
# from Engine.opt_intense import optimizer, fmod, fmod_conti
from Engine.outplotter import outplotter_23
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def rv_MPinst(args, inparam, orders, order_use, trk, step2or3, i):

    # Main function for RV fitting that will be threaded over by multiprocessing

    nights   = inparam.nights
    night = nights[i] # current looped night

    order = orders[order_use]
    xbounds = inparam.xbounddict[order]
    print('Working on order {:02d}/{:02d} ({}), night {:03d}/{:03d} ({}) PID:{}...'.format(int(order_use)+1,
                                                                                           len(orders),
                                                                                           order,
                                                                                           i+1,
                                                                                           len(inparam.nights),
                                                                                           night,
                                                                                           mp.current_process().pid) )

    #-------------------------------------------------------------------------------

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    nightsout = [];
    rvsminibox     = np.ones(len(tagsnight));
    vsiniminibox   = np.ones(len(tagsnight));
    tagsminibox    = np.ones(len(tagsnight));
    parfitminibox  = np.ones((len(tagsnight),15)); # need to match the dpar numbers

    rvsminibox[:]    = np.nan
    vsiniminibox[:]  = np.nan
    tagsminibox[:]   = np.nan
    parfitminibox[:] = np.nan

    for t in tagsnight:
        nightsout.append(night)

#-------------------------------------------------------------------------------
    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == float:
        initguesses = inparam.initguesses
    else:
        sys.exit('ERROR! EXPECING SINGAL NUMBER OR FILE FOR INITGUESSES! QUITTING!')


    # Load synthetic telluric template generated during Step 1
    # [:8] here is to ensure program works under Night_Split mode
    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_treated_{args.band}.fits'
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox, tagsminibox

    # Find corresponding table in fits file, given the tables do not go sequentially by order number due to multiprocessing in Step 1
    num_orders = 0
    for i in range(25):
        try:
            hdulist[i].columns[0].name[9:]
            num_orders += 1
        except:
            continue

    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

    # Check whether Telfit hit critical error in Step 1 for the chosen order with this night. If so, skip.
    if flag == 1:
        logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}, skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox, tagsminibox


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
    if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]
#-------------------------------------------------------------------------------
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,                                                # 0: The shift of the stellar template (km/s) [assigned later]
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


    # Iterate over all A/B exposures
    for t in np.arange(len(tagsnight)):
        tag = tagsnight[t]
        beam = beamsnight[t]

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
        x,wave,s,u = init_fitsread(f'{inparam.inpath}{night}/{beam}/',
                                    'target',
                                    'separate',
                                    night,
                                    order,
                                    tag,
                                    args.band,
                                    bound_cut)

        #-------------------------------------------------------------------------------

        # Execute S/N cut
        s2n = s/u
        if np.nanmedian(s2n) < float(args.SN_cut):
            logger.warning('  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
            continue

        # Trim obvious outliers above the blaze (i.e. cosmic rays)
        nzones = 5
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

        # Cut spectrum to within wavelength regions defined in input list
        s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
        x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

        # Trim stellar template to relevant wavelength range
        mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

        # Trim telluric template to relevant wavelength range
        satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
        watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

        # Make sure data is within telluric template range (shouldn't do anything)
        s_piece    = s_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        u_piece    = u_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        x_piece    = x_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

        # --------------------------------------------------------------

        par = pars0.copy()

        # Get initial guess for cubic wavelength solution from reduction pipeline
        f = np.polyfit(x_piece,wave_piece,2)
        par8in = f[0]*1e4; par7in = f[1]*1e4; par6in = f[2]*1e4;
        par[9] = 0.; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = initguesses-inparam.bvcs[night+tag] # Initial RV with barycentric correction

        # Arrays defining parameter variations during optimization steps
        dpars = {'cont' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                 'wave' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0,    0,   0, 0, 0,    0]),
                 't'    : np.array([0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 'ip'   : np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 's'    : np.array([5.0, 1.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 'v'    : np.array([0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        s_piece /= np.median(s_piece)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]))

        #-------------------------------------------------------------------------------

        # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
        optimize = True
        par_in = par.copy()
        hardbounds = [par_in[4]-dpars['v'][4],  par_in[4]+dpars['v'][4],
                      par_in[5]-dpars['ip'][5], par_in[5]+dpars['ip'][5]]
        if hardbounds[0] < 0:
            hardbounds[0] = 0
        if hardbounds[3] < 0:
            hardbounds[3] = 1

        # Begin optimization. Fit the blaze, the wavelength solution, the telluric template power and RV, the stellar template power and RV, the
        # zero point for the instrumental resolution, and the vsini of the star separately, iterating and cycling between each set of parameter fits.

        cycles = 4

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

            for optkind in optgroup:
                parfit_1 = optimizer(parstart, dpars[optkind], hardbounds, fitobj, optimize)
                parstart = parfit_1.copy()
                if args.debug == True:
                    outplotter_23(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(order,night,tag,nk,optkind), trk, inparam, args, step2or3)
                    logger.debug(f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
                nk += 1

        parfit = parfit_1.copy()

        #-------------------------------------------------------------------------------

        # if best fit stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            logger.warning(f'  --> parfit[1] < 0.1, {night} parfit={parfit}')
            continue

        # if best fit stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
            logger.warning(f'  --> parfit[1] == par_in[1] or parfit[3] == par_in[3], {night}')
            continue

        # if best fit model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq = fmod(parfit,fitobj)
        if len(smod[(smod < 0)]) > 0:
            logger.warning(f'  --> len(smod[(smod < 0)]) > 0, {night}')
            continue

        #-------------------------------------------------------------------------------

        if args.plotfigs == True:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)

        rv0 = parfit[0] - parfit[2]   # Correct for RV of the atmosphere, since we're using that as the basis for the wavelength scale

        rvsminibox[t]   = rv0  + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2) # Barycentric correction
        parfitminibox[t]= parfit
        vsiniminibox[t] = parfit[4]
        tagsminibox[t]  = tag
    return nightsout,rvsminibox,parfitminibox,vsiniminibox,tagsminibox

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 3',
                                     description = '''
                                     Performs a full analysis of each target star observation to produce accurate and precise RVs.\n
                                     All the wavelength regions defined in Step 1 are used, and the code analyzes each exposure that is part of a given observation separately (this allows estimates of the RV uncertainties). \n \n
                                     Unless the target vsini is already known to high accuracy, an initial run of Step 3 is required. \n
                                     This  will provide an estimate of the \vsini, which will then be plugged in as a fixed value in the second run of Step 3. \n \n
                                     Even if vsini is already known, it is recommended to run Step 3 twice, as it allows the user to confirm that RVs between consecutive runs have converged within the calculated uncertainties. \n \n
                                     Additional runs may be necessary depending on the quality of the observed spectra and the \vsini of the target star.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-mode",    dest="mode",             action="store",
                        help="RV standard star (STD) or a normal target (TAR)?",
                        type=str,   default='')
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 ",
                        type=str,   default='50')
    parser.add_argument("-nAB",      dest="nAB",           action="store",
                        help="Minium number of separte A/B exposures within a set for a given observation (ensures accuracy of uncertainy estimates). Default = 2 for STD, 3 for TAR",
                        type=str,   default='')
    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). If no literature value known, use the value given by Step 2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 5.0 km/s. Should be set to 0 for final run.",
                        type=str, default='5.0' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help="For STD star. Initial RV guess for all nights. Given by Step 2 results (float, km/s)",
                        type=str,   default='' )
    parser.add_argument('-gS',       dest="guesses_source",           action="store",
                        help="For TAR star. Source for list of initial RV guesses. 'init' = Initguesser_results_X = past Step 2 result OR 'rvre' = RV_results_X = past Step 3 result",
                        type=str, default='')
    parser.add_argument('-gX',       dest="guessesX",           action="store",
                        help="For TAR star. The number, X, under ./*targname/Initguesser_results_X or ./*targname/RV_results_X, that you wish to use. Prefix determined by -gS",
                        type=str, default='')
    parser.add_argument('-t',       dest="template",         action="store",
                        help="Stellar template. Pick from 'synthetic', 'livingston', or 'user_defined'. Default = 'synthetic'",
                        type=str,   default='synthetic' )
    parser.add_argument('-sp',      dest="sptype",           action="store",
                        help="The spectral type of the *target (letter only).",
                        type=str,   default='' )
    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If set, will generate plots of the fitting results under ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots under ./Temp/Debug/*target_*band/main_step3")
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.85')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    # Check user input

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)

    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'user_defined']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    if args.template.lower() in ['synthetic', 'livingston']:
        if args.sptype not in ['F','G','K','M']:
            sys.exit('ERROR: DEFAULT TEMPLATES ONLY COVER F G K M STARS!')

    #------------------------------

    if (args.mode.lower() == 'std') & (args.guesses_source != ''):
        sys.exit('ERROR: STD CANNOT USE -gS, PLEASE USE -g')
    if (args.mode.lower() == 'tar') & (args.guesses != ''):
        sys.exit('ERROR: TAR CANNOT USE -g, PLEASE USE -gS')

    #------------------------------

    if args.nAB == '' and args.mode.lower() == 'std':
        nAB = 2
    elif args.nAB == '' and args.mode.lower() == 'tar':
        nAB = 3
    else:
        nAB = int(args.nAB)

    #------------------------------

    if args.mode.lower() == 'std': # Specify initial RV guesses as a single value applied to all nights
        initguesses = float(args.guesses)
        initguesses_show = initguesses
    else: # Load initial RV guesses from file
        if args.guesses_source == 'init': # From Step 2 results
            guesses = './Output/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
                                                                         args.band,
                                                                         int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre': # From Step 3 results
            guesses = './Output/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
                                                                      args.band,
                                                                      int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs    = np.array(guessdata['RVfinal'])
            initguesses = {}
            initguesses_show = f'RVresultsSummary_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[41m {} band \033[0m
    WaveLength file     = \33[41m WaveRegions_{} \033[0m
    S/N cut             > \33[41m {} \033[0m
    Minium # of AB sets = \33[41m {} \033[0m             <------- If TAR mode, this should be at least 3. If STD mode, at least 2.
    Initial vsini       = \33[41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[41m {} km/s \033[0m
    RV initial guess    = \33[41m {} \033[0m
    Stellar template use= \33[41m {} \033[0m
    Target Spectral Type= \33[41m {} \033[0m             <-------  [late K, M] recommended 'synthetic', [F, G, early K] SpTy recommended 'livingston'
    Threads use         = {}
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut, nAB,
               initvsini, vsinivary, initguesses_show, args.template, args.sptype, args.Nthreads))
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
    print('Running Step 3 for {}...'.format(args.targname))
    print('This will take a while..........')

    #-------------------------------------------------------------------------------

    # Make output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')
    filesndirs = os.listdir(f'./Output/{args.targname}_{args.band}')

    trk = 1; go = True;
    while go == True:
        name = f'RV_results_{trk}'
        if name not in filesndirs:
            break
        trk += 1

    os.mkdir(f'./Output/{args.targname}_{args.band}/{name}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs')

    step2or3 = 3
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

    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{name}')

    #-------------------------------------------------------------------------------

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders = read_prepdata(args)

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
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, args.sptype)

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)

    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T).
    # All statistical analysis will be performed separately for these two datasets.
    nights    = inparam.nights
    intnights = np.array([int(i[:8]) for i in nights])

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]
    rvmasterboxT  = np.ones((len(nightsT),len(orders)))
    stdmasterboxT = np.ones((len(nightsT),len(orders)))
    rvmasterboxL  = np.ones((len(nightsL),len(orders)))
    stdmasterboxL = np.ones((len(nightsL),len(orders)))
    vsinisT = np.ones((len(nightsT),len(orders)))
    vsinisL  = np.ones((len(nightsL),len(orders)))

    if len(nightsL) > 0:
        nightscomblist = [nightsT,nightsL]
    else:
        nightscomblist = [nightsT]


    #-------------------------------------------------------------------------------

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        pool = mp.Pool(processes = args.Nthreads)
        func = partial(rv_MPinst, args, inparam, orders, jerp, trk, step2or3 )
        outs = pool.map(func, np.arange(len(nightsFinal)))
        pool.close()
        pool.join()

        order = orders[jerp]

        #-------------------------------------------------------------------------------

        # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
        for i in range(len(nights)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvbox     = outsbox[1]
                parfitbox = outsbox[2]
                vsinibox  = outsbox[3]
                tagbox    = outsbox[4]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvbox     = np.concatenate((rvbox,outsbox[1]))
                parfitbox = np.vstack((parfitbox,outsbox[2]))
                vsinibox  = np.concatenate((vsinibox,outsbox[3]))
                tagbox    = np.concatenate((tagbox,outsbox[4]))

        nightsbox = np.array(nightsbox)
        vsinitags = []

        # Save results to fits file
        c1    = fits.Column(name='NIGHT'+str(order),  array=nightsbox, format='{}A'.format(len(nights[0])) )
        c2    = fits.Column(name='RV'+str(order),     array=rvbox,     format='D')
        c3    = fits.Column(name='PARFIT'+str(order), array=parfitbox, format=str(len(parfitbox[0,:]))+'D', dim=(1,len(parfitbox[0,:])))
        c4    = fits.Column(name='VSINI'+str(order),  array=vsinibox,  format='D')
        c5    = fits.Column(name='TAG'+str(order),    array=tagbox,    format='4A')
        cols  = fits.ColDefs([c1,c2,c3,c4,c5])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if jerp == 0: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
        else:
            hh = fits.open('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
            hh.append(hdu_1)
            hh.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name),overwrite=True)

        #-------------------------------------------------------------------------------

        # For each set of nights (tight, loose)...
        T_L = 'T'
        for nights_use in nightscomblist:

            # For each night...
            for i in range(len(nights_use)):
                # Collect the RVs and vsinis determined from different A/B exposures within a night
                indnight  = np.where(nightsbox == nights_use[i])[0]
                rvtags    = rvbox[indnight]
                vsinitags = vsinibox[indnight]

                # Take the mean of the vsinis, and the mean and std of the RVs.
                # If the number of different successfully fitted A/B exposures is less than required, pass NaN instead.
                if T_L == 'T':
                    vsinisT[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxT[i,jerp]  = np.nan
                        stdmasterboxT[i,jerp] = np.nan
                    else:
                        rvmasterboxT[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxT[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))

                else:
                    vsinisL[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxL[i,jerp]  = np.nan
                        stdmasterboxL[i,jerp] = np.nan
                    else:
                        rvmasterboxL[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxL[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))
            T_L = 'L'


    #-------------------------------------------------------------------------------

    # Don't combine Loose and Tight datasets, but make them both easily referenceable
    nightsCombined  = np.array([]); mjdsCombined = np.array([]);
    rvfinalCombined = np.array([]); stdfinalCombined = np.array([]); vsinifinalCombined = np.array([]);

    if len(nightsL) > 0:
        rvboxcomblist  = [rvmasterboxT,rvmasterboxL]
        stdboxcomblist = [stdmasterboxT,stdmasterboxL]
        vsinicomblist  = [vsinisT,vsinisL]
    else:
        rvboxcomblist  = [rvmasterboxT]
        stdboxcomblist = [stdmasterboxT]
        vsinicomblist  = [vsinisT]

    # Iterate over tight and loose mounting data sets...
    for boxind in range(len(rvboxcomblist)):

        rvmasterbox  = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]
        vsinibox     = vsinicomblist[boxind]

        #-------------------------------------------------------------------------------

        if args.mode=='STD': # If RV STD star, calculate uncertainty in method

            # Calculate the precision within an order across nights
            sigma_O2     = np.array([np.nanstd(rvmasterbox[:,ll])**2 for ll in range(len(orders))])
            sigma_ABbar2 = np.ones_like(sigma_O2)

            # Calculate uncertainty in method as difference between variance within an order and mean variance within a night's different exposure RVs
            for ll in range(len(orders)):
                sigma_ABbar2[ll] = np.nanmean(stdmasterbox[:,ll]**2)
            sigma_method2 = sigma_O2 - sigma_ABbar2

        else: # If target star, load the uncertainty in method calculated from our RV STD star runs
            if boxind == 0:
                nights_use = nightsT.copy()
                kind = 'Tight'
                sigma_method2 = inparam.methodvariance_tight[args.band]
            else:
                nights_use = nightsL.copy()
                kind = 'Loose'
                sigma_method2 = inparam.methodvariance_loose[args.band]

        sigma_ON2    = np.ones_like(rvmasterbox)

        #-------------------------------------------------------------------------------

        # Note rvmasterbox indexed as [nights,orders]
        Nnights = len(rvmasterbox[:,0])

        # Calculate the uncertainty in each night/order RV as the sum of the uncertainty in method and the uncertainty in that night's As and Bs RVs
        for ll in range(len(orders)):
            for night in range(Nnights):
                sigma_ON2[night,ll] = sigma_method2[ll] + stdmasterbox[night,ll]**2

        rvfinal    = np.ones(Nnights, dtype=np.float64)
        stdfinal   = np.ones(Nnights, dtype=np.float64)
        vsinifinal = np.ones(Nnights, dtype=np.float64)
        mjds_out   = np.ones(Nnights, dtype=np.float64)

        if boxind == 0:
            nights_use = nightsT.copy(); kind = 'Tight';
        else:
            nights_use = nightsL.copy(); kind = 'Loose';


        # Combine RVs between orders using weights calculated from uncertainties
        for n in range(Nnights):
            weights = (1./sigma_ON2[n,:]) / (np.nansum(1./sigma_ON2[n,:])) # normalized
            stdspre = (1./sigma_ON2[n,:]) #unnormalized weights

            rvfinal[n]  = np.nansum( weights*rvmasterbox[n,:] )
            stdfinal[n] = 1/np.sqrt(np.nansum(stdspre))

            vsinifinal[n] = np.nansum(weights*vsinibox[n,:])
            mjds_out[n]   = mjds[nights_use[n]]

            # if all the RVs going into the observation's final RV calculation were NaN due to any pevious errors, pass NaN
            if np.nansum(weights) == 0:
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan

            # if more than half of the orders going into the observation's final RV calculation were NaN due to any pevious errors, pass NaN
            if np.sum( np.isnan(rvmasterbox[n,:]) ) > np.floor( len(orders) * 0.5 ):
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan

        #-------------------------------------------------------------------------------

        # Plot results
        f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)

        axes.plot(    np.arange(len(rvfinal))+1, rvfinal, '.k', ms=5)
        axes.errorbar(np.arange(len(rvfinal))+1, rvfinal, yerr=stdfinal, ls='none', lw=.5, ecolor='black')
        axes.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinal), np.nanstd(rvfinal)),
                             transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        axes.set_ylim(np.nanmin(rvfinal)-.08,
                     np.nanmax(rvfinal)+.08)
        axes.set_ylabel('RV [km/s]', size=6, style='normal', family='sans-serif' )
        axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(5))
        axes.tick_params(axis='both', which='both', labelsize=5, right=True, top=True, direction='in', width=.6)
        f.savefig('{}/{}/FinalRVs_{}_.png'.format(inparam.outpath, name, kind), format='png', bbox_inches='tight')

        # Save results to fits file separately for each tight/loose dataset
        c1 = fits.Column( name='NIGHT',         array=nights_use,    format='8A')
        c2 = fits.Column( name='JD',            array=mjds_out,      format='D')
        c3 = fits.Column( name='RVBOX',         array=rvmasterbox,   format='{}D'.format(len(orders)))
        c4 = fits.Column( name='STDBOX',        array=stdmasterbox,  format='{}D'.format(len(orders)))
        c7 = fits.Column( name='Sigma_method2', array=sigma_method2, format='D')
        c8 = fits.Column( name='Sigma_ON2',     array=sigma_ON2,     format='{}D'.format(len(orders)))
        c9 = fits.Column( name='RVfinal',       array=rvfinal,       format='D')
        c10 = fits.Column(name='STDfinal',      array=stdfinal,      format='D')

        if args.mode=='STD':
            c5 = fits.Column( name='Sigma_O2',      array=sigma_O2,      format='D')
            c6 = fits.Column( name='Sigma_ABbar2',  array=sigma_ABbar2,  format='D')
            cols  = fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
        else:
            cols  = fits.ColDefs([c1,c2,c3,c4,c7,c8,c9,c10])

        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul        = fits.HDUList([primary_hdu,hdu_1])
        hdul.writeto('{}/{}/RVresultsSummary_{}.fits'.format(inparam.outpath, name, kind), overwrite=True)

        # Combine final RVs from both tight and loose mounting data sets
        nightsCombined     = np.concatenate((nightsCombined,     nights_use))
        mjdsCombined       = np.concatenate((mjdsCombined,       mjds_out))
        rvfinalCombined    = np.concatenate((rvfinalCombined,    rvfinal))
        stdfinalCombined   = np.concatenate((stdfinalCombined,   stdfinal))
        vsinifinalCombined = np.concatenate((vsinifinalCombined, vsinifinal))

        if args.mode=='STD': # If uncertainty in method was calculated, save it
            sigma_method2 = [np.around(float(i), 8) for i in sigma_method2]
            logger.info('sigma_method2 with type = {} is {}'.format(kind, sigma_method2))
        logger.info('During the {} mounting period: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( kind,
                                                                                                        np.nanmean(rvfinal),
                                                                                                        np.nanstd(rvfinal) ))

    #-------------------------------------------------------------------------------

    # Plot combined results
    xscale = np.arange(len(rvfinalCombined))+1

    f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(xscale,rvfinalCombined, '.k', ms=5)
    axes.errorbar(xscale,rvfinalCombined,yerr=stdfinalCombined,ls='none',lw=.5, ecolor='black')
    axes.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=axes.transAxes, size=6, style='normal', family='sans-serif' )

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        axes.text(0.05, 0.1, 'Tight', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        axes.text(0.05, 0.1, 'Loose', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            axes.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Tight', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Loose', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        else:
            axes.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Tight', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Loose', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    axes.set_ylim(np.nanmin(rvfinalCombined)-.08,np.nanmax(rvfinalCombined)+.08)
    axes.set_ylabel('RV (km/s)', size=6, style='normal', family='sans-serif' )
    axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/{}/FinalRVs.png'.format(inparam.outpath, name), format='png', bbox_inches='tight')

    # Output combined final results to fits file
    c1 = fits.Column(name='NIGHT',    array=nightsCombined,         format='{}A'.format(len(nights[0])) )
    c2 = fits.Column(name='JD',       array=mjdsCombined,           format='D')
    c3 = fits.Column(name='RVfinal',  array=rvfinalCombined,        format='D')
    c4 = fits.Column(name='STDfinal', array=stdfinalCombined,       format='D')
    c5 = fits.Column(name='VSINI',    array=vsinifinalCombined,     format='D')

    cols = fits.ColDefs([c1,c2,c3,c4,c5])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    bleh = np.ones((3,3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu,hdu_1])
    hdul.writeto('{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name),overwrite=True)

    tempin = Table.read('{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name), format='fits')
    tempin.write('{}/RVresultsSummary_{}.csv'.format(inparam.outpath, trk), format='csv', overwrite=True)

    logger.info('Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(rvfinalCombined),
                                                                            np.nanstd(rvfinalCombined)))
    logger.info('vsini results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(vsinifinalCombined),
                                                                            np.nanstd(vsinifinalCombined)))
    print('\n')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}/{}'.format(inparam.outpath, name) )
    print('####################################################################################')
