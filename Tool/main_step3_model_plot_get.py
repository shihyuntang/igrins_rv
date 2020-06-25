import sys
sys.path.append("..") # Adds higher directory to python modules path.

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
from Engine.outplotter import outplotter_23
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def rv_MPinst(args, inparam, i, orders, order):
    nights   = inparam.nights
    night = i[0]

    xbounds = inparam.xbounddict[order]
    print('Working on order {:02d}/{:02d} ({}), night {} PID:{}...'.format(int(order)+1,
                                                                               len(orders),
                                                                                order,
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
#-------------------------------------------------------------------------------
    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    wminibox      = np.ones(2048)
    sminibox      = np.ones(2048)
    flminibox_tel = np.ones(2048)
    flminibox_ste = np.ones(2048)
    contiminibox  = np.ones(2048)
    flminibox_mod = np.ones(2048)

    wminibox[:]     = np.nan
    sminibox[:]     = np.nan
    flminibox_tel[:]= np.nan
    flminibox_ste[:]= np.nan
    contiminibox[:] = np.nan
    flminibox_mod[:]  = np.nan

    # Load telluric template from Telfit'd A0
    A0loc = f'../Output/{args.targname}_{args.band}_tool/A0Fits/{night[:8]}A0_treated_{args.band}.fits'
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        sys.exit(f'  --> No A0-fitted template for night {night}, skipping...')

    num_orders = 0
    for i in range(25):
        try:
            hdulist[i].columns[0].name[9:]
            num_orders += 1
        except:
            continue

    # order in A0_treated.fits is no longer sequential...
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

    if flag == 1:  # Telfit hit unknown critical error
        lsys.exit(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}, skipping...')

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
    pars0 = np.array([np.nan,                                                # 0: The shift of the sunspot spectrum (km/s) [assigned later]
                      0.3,                                                   # 1: The scale factor for the sunspot spectrum
                      0.0,                                                   # 2: The shift of the telluric spectrum (km/s)
                      0.6,                                                   # 3: The scale factor for the telluric spectrum
                      inparam.initvsini,                                     # 4: vsini (km/s)
                      IPpars[2],                                             # 5: The instrumental resolution (FWHM) in pixels
                      np.nan,                                                # 6: Wavelength 0-pt
                      np.nan,                                                # 7: Wavelength linear component
                      np.nan,                                                # 8: Wavelength quadratic component
                      np.nan,                                                # 9: Wavelength cubic component
                      1.0,                                                   #10: Continuum zero point
                      0.,                                                    #11: Continuum linear component
                      0.,                                                    #12: Continuum quadratic component
                      IPpars[1],                                             #13: IP linear component
                      IPpars[0]])                                              #14: IP quadratic component


    # Iterate over all A/B exposures
    for t in [0]:
        tag = tagsnight[t]
        beam = beamsnight[t]

        ### Load relevant A0 spectrum
        if args.band=='K':
            if int(order) in [11, 12, 13, 14]:
                bound_cut = inparam.bound_cut_dic[args.band][order]
            else:
                bound_cut = [150, 150]

        elif args.band=='H':
            if int(order) in [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
                bound_cut = inparam.bound_cut_dic[args.band][order]
            else:
                bound_cut = [150, 150]

        x,wave,s,u = init_fitsread(f'{inparam.inpath}{night}/{beam}/',
                                    'target',
                                    'separate',
                                    night,
                                    order,
                                    tag,
                                    args.band,
                                    bound_cut)
#-------------------------------------------------------------------------------
        s2n = s/u
        if np.nanmedian(s2n) < float(args.SN_cut):
            sys.exit('  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))


        nzones = 5
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

        s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
        x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

        mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

        satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
        watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

        # Cut target spec to be within A0 spec wave
        s_piece    = s_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        u_piece    = u_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        x_piece    = x_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
        wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

        # --------------------------------------------------------------
        par = pars0.copy()
        ##Set initial wavelength guess
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = initguesses-inparam.bvcs[night+tag]
        # Arrays defining parameter variations during optimization steps
        dpars = {'cont' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                 'wave' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0]),
                 't'    : np.array([0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 'ip'   : np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 's'    : np.array([5.0, 1.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                 'v'    : np.array([0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        s_piece /= np.median(s_piece)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]))
#-------------------------------------------------------------------------------
        ######## Begin optimization  ########
        optimize = True
        par_in = par.copy()
        hardbounds = [par_in[4]-dpars['v'][4],  par_in[4]+dpars['v'][4],
                      par_in[5]-dpars['ip'][5], par_in[5]+dpars['ip'][5]]
        if hardbounds[0] < 0:
            hardbounds[0] = 0
        if hardbounds[3] < 0:
            hardbounds[3] = 1

        cycles = 4

        optgroup = ['cont', 'wave', 't', 'cont', 's',
                    'cont', 'wave', 't', 's', 'cont',
                    'wave',
                    'ip', 'v',
                    'ip', 'v',
                    't',  's',
                    't',  's']


        for nc, cycle in enumerate(np.arange(cycles), start=1):
            if cycle == 0:
                parstart = par_in.copy()

            for optkind in optgroup:
                parfit_1 = optimizer(parstart, dpars[optkind], hardbounds, fitobj, optimize)
                parstart = parfit_1.copy()
                if args.debug == True:
                    outplotter_23(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(order,night,tag,nc,optkind), trk, inparam, args, step2or3)
                    logger.debug(f'{order}_{tag}_{nc}_{optkind}:\n {parfit_1}')

        parfit = parfit_1.copy()
#-------------------------------------------------------------------------------
        # if stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            sys.exit(f'  --> parfit[1] < 0.1, {night} parfit={parfit}')
            continue
        # if stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
            sys.exit(f'  --> parfit[1] == par_in[1] or parfit[3] == par_in[3], {night}')
            continue
        # if model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq = fmod(parfit,fitobj)
        if len(smod[(smod < 0)]) > 0:
            sys.exit(f'  --> len(smod[(smod < 0)]) > 0, {night}')
            continue
#-------------------------------------------------------------------------------
        # Compute model and divide for residual
        fullmodel,chisq = fmod(parfit,fitobj)

        # Set both stellar and telluric template powers to 0 to compute only continuum
        parcont = parfit.copy();
        parcont[1] = 0.; parcont[3] = 0.;
        contmodel, chisq = fmod(parcont,fitobj)

        # Set stellar tempalte power to 0 to compute only telluric, and vice versa
        parS = parfit.copy(); parT = parfit.copy();
        parT[1] = 0.; parS[3] = 0.;
        stellmodel,chisq = fmod(parS,   fitobj)
        tellmodel, chisq = fmod(parT,   fitobj)

        # Divide everything by continuum model except residual
        dataflat  = fitobj.s/contmodel
        modelflat = fullmodel/contmodel
        stellflat = stellmodel/contmodel
        tellflat  = tellmodel/contmodel

    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)




    wminibox[:n]         = w
    sminibox[:]          = dataflat
    flminibox_mod[:ln]   = modelflat
    flminibox_tel[:n]    = tellflat
    flminibox_ste[:nn]   = stellflat
    contiminibox[:n]     = contmodel
    residualbox[:l]      = residual

    return wminibox,sminibox,flminibox_mod,flminibox_tel,flminibox_ste,contiminibox,residualbox


def mp_run(Nthreads, args, inparam, nights, order0):
    pool = mp.Pool(processes = Nthreads)
    func = partial(rv_MPinst, args, inparam, nights, order0)
    outs = pool.map(func, order0)
    pool.close()
    pool.join()
    return outs

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline',
                                     description = '''
                                     This is a pipeline that helps you to extract radial velocity \n
                                     from IGRINS spectra. \n
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-mode",    dest="mode",             action="store",
                        help="RV standard star (STD) OR a normal target (TAR)??.",
                        type=str,   default='')
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 1",
                        type=int,   default=int(1))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Default = 50 ",
                        type=str,   default='50')
    parser.add_argument("-nAB",      dest="nAB",           action="store",
                        help="Minium request of # of AB sets. Default = 2",
                        type=str,   default='1')

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). Should use the value given by step2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 5.0 km/s",
                        type=str, default='5.0' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help="For STD star, single value given by step2 (float, km/s)",
                        type=str,   default='' )

    parser.add_argument('-gS',       dest="guesses_source",           action="store",
                        help="For TAR star, Source for initial guesses list for RV. Enter init OR rvre (init: Initguesser_results_X, rvre: RV_results_X)",
                        type=str, default='')
    parser.add_argument('-gX',       dest="guessesX",           action="store",
                        help="For TAR star, Please give the number, X, under ./*targname/Initguesser_results_X OR ./*targname/RV_results_X, that you wish to use",
                        type=str, default='')

    parser.add_argument('-t',       dest="template",         action="store",
                        help="Stellar template. Pick from 'synthetic', 'livingston', or 'user_defined'. Default = 'synthetic'",
                        type=str,   default='synthetic' )
    parser.add_argument('-sp',      dest="sptype",           action="store",
                        help="The spectral type of the *target.",
                        type=str,   default='' )
    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If sets, will generate basic fitting result plots")

    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want all process all nights under the ./Input/*target folder, give an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If sets, DeBug logging and extra plots will be given")
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If sets, will skip the input parameters check, handly when run muti-targets line by line ")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.85')
    args = parser.parse_args()
    inpath   = '../Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'
#-------------------------------------------------------------------------------
    # INPUT CHECKING...
    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)
    #------------------------------
    if args.template.lower() not in ['synthetic', 'livingston', 'user_defined']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')
    #------------------------------
    if args.template.lower() in ['synthetic', 'livingston']:
        if args.sptype not in ['F','G','K','M']:
            sys.exit('ERROR: SPECTRAL TYPE FOR DEFAULT TEMPALTE ARE ONLY FOR F G K M! STARS')
    #------------------------------
    if (args.mode.lower() == 'std') & (args.guesses_source != ''):
        sys.exit('ERROR: STD CANNOT USE -gS, PLEASE USE -g')
    if (args.mode.lower() == 'tar') & (args.guesses != ''):
        sys.exit('ERROR: TAR CANNOT USE -g, PLEASE USE -gS')
    #------------------------------

    if args.mode.lower() == 'std':
        initguesses = float(args.guesses)
        initguesses_show = initguesses
    else:
        if args.guesses_source == 'init':
            guesses = '../Output/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
                                                                         args.band,
                                                                         int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre':
            guesses = '../Output/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
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
    Minium # of AB sets = \33[41m {} \033[0m             <------- If is TAR mode, this should be at lease 2
    Initial vsini       = \33[41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[41m {} km/s \033[0m
    RV initial guess    = \33[41m {} \033[0m
    Stellar template use= \33[41m {} \033[0m
    Target Spectral Type= \33[41m {} \033[0m             <-------  [late K, M] recommended 'synthetic', [F, G, early K] SpTy recommended 'livingston'
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut, args.nAB,
               initvsini, vsinivary, initguesses_show, args.template, args.sptype))
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
    print('RV calculation for RV standard star {}...'.format(args.targname))
    print('This will take a while..........')
#-------------------------------------------------------------------------------
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'../Output/{args.targname}_{args.band}_tool'):
        os.mkdir(f'../Output/{args.targname}_{args.band}_tool')
    filesndirs = os.listdir(f'../Output/{args.targname}_{args.band}_tool')

    trk = 1; go = True;
    while go == True:
        name = f'RV_results_{trk}'
        if name not in filesndirs:
            break
        trk += 1

    os.mkdir(f'../Output/{args.targname}_{args.band}_tool/{name}')

    outpath = f'../Output/{args.targname}_{args.band}_tool'
#-------------------------------------------------------------------------------
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}_tool.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
#-------------------------------------------------------------------------------
    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}_tool/{name}')
#-------------------------------------------------------------------------------
    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders = read_prepdata(args)

# ---------------------------------------
#     if args.band == 'K':
#         orders = np.append(np.arange(2, 9), np.array([10, 11, 12, 13, 14, 16]))
#     elif args.band == 'H':
# #        order0 = np.arange(5,11)
#         # order0 = np.arange(2,23)
#         orders = np.array([2, 3, 4, 5, 6, 10, 11, 13, 14, 16, 17, 20, 21, 22])
# #    order0 = np.array([16])
# ---------------------------------------
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))

    if len(nightsFinal)==1:
        print(nightsFinal, len(nightsFinal))
    else:
        sys.exit('only take one night!, we give {}'.format(nightsFinal))

    logger.info('Analyze with {} nights'.format(len(nightsFinal)))
#-------------------------------------------------------------------------------
    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, args.sptype)

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)
#-------------------------------------------------------------------------------
    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T)
    # nights    = inparam.nights
    # intnights = np.array([int(i[:8]) for i in nights])
    #
    # indT = np.where((intnights < 20180401) | (intnights > 20190531))
    # indL = np.where((intnights >= 20180401) & (intnights < 20190531))
    #
    # nightsT = nights[indT]
    # nightsL = nights[indL]
    # rvmasterboxT  = np.ones((len(nightsT),len(orders)))
    # stdmasterboxT = np.ones((len(nightsT),len(orders)))
    # rvmasterboxL  = np.ones((len(nightsL),len(orders)))
    # stdmasterboxL = np.ones((len(nightsL),len(orders)))
    # vsinisT = np.ones((len(nightsT),len(orders)))
    # vsinisL  = np.ones((len(nightsL),len(orders)))
    #
    # if len(nightsL) > 0:
    #     nightscomblist = [nightsT,nightsL]
    # else:
    #     nightscomblist = [nightsT]
    nights    = inparam.nights
#-------------------------------------------------------------------------------
    outs = mp_run(args.Nthreads, args, inparam, nights, orders)

# Collect outputs
for i in range(len(orders)):
    outsbox = outs[i]

    wbox      = outsbox[0]
    stbox     = outsbox[1]
    modbox    = outsbox[2]
    telbox    = outsbox[3]
    stebox    = outsbox[4]
    conti_fl  = outsbox[5]

    # Save results in fits file
    c1 = fits.Column(name='wavelength',    array=wbox,         format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
    c2 = fits.Column(name='s',             array=stbox,        format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
    c3 = fits.Column(name='model_fl',      array=modbox,    format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
    c4 = fits.Column(name='tel_fl',        array=telbox,       format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
    c5 = fits.Column(name='ste_fl',        array=stebox,       format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
    c6 = fits.Column(name='conti_fl',      array=conti_fl,     format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))


    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    if orders[i] == orders[0]:  # If first time writing fits file, make up filler primary hdu
        bleh = np.ones((3, 3))
        primary_hdu1 = fits.PrimaryHDU(bleh)
        hdul = fits.HDUList([primary_hdu1, hdu_1])
        hdul.writeto(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(args.targname, inparam.nights[0], args.band))
    else:
        hh = fits.open(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(args.targname, inparam.nights[0], args.band))
        hh.append(hdu_1)
        hh.writeto(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(args.targname, inparam.nights[0], args.band), overwrite=True)

end_time = datetime.now()
print('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
print('Output saved under {}/{}'.format(args.targname, name) )
print('##########################################################')
