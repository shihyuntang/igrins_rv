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
    nights   = inparam.nights
    night    = nights[i]

    order   = orders[order_use]
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
        sys.exit('ERROR! EXPECTED FILE OR LIST FOR INITGUESSES! QUITTING!')

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    # Load telluric template from Telfit'd A0
    # [:8] here is to ensure program works under Night_Split mode
    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_treated_{args.band}.fits'
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
        return night, np.nan, np.nan

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
        logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}, skipping...')
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

    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]
#-------------------------------------------------------------------------------
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,                                                # 0: The shift of the sunspot spectrum (km/s)
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
                      IPpars[0]])                                            #14: IP quadratic component

    rvsmini = []; vsinismini = [];
    # Iterate over all A/B exposures
    for t in np.arange(len(tagsnight)):
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
            logger.warning('  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
            continue

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

#-------------------------------------------------------------------------------
        par = pars0.copy()
        ##Set initial wavelength guess
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ;  par[8] = par8in ;  par[7] = par7in ;  par[6] = par6in

        par[0] = initguesses-inparam.bvcs[night+tag]
        # Arrays defining parameter variations during optimization steps
        dpars1 = {'cont' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                  'wave' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0.,   0,   0, 0, 0,    0]),
                  't'    : np.array([ 0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  'ip'   : np.array([ 0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  's'    : np.array([20.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  'v'    : np.array([ 0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}
        dpars2 = {'cont' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                  'wave' : np.array([ 0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0.,   0,   0, 0, 0,    0]),
                  't'    : np.array([ 0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  'ip'   : np.array([ 0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  's'    : np.array([ 5.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                  'v'    : np.array([ 0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        s_piece /= np.median(s_piece)

        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]))

        mask = np.ones_like(s_piece,dtype=bool)
        mask[(fitobj.s < .0)] = False
#-------------------------------------------------------------------------------
        ######## Begin optimization  ########

        optimize = True
        par_in = par.copy()
        hardbounds = [par_in[4]-dpars1['v'][4],  par_in[4]+dpars1['v'][4],
                      par_in[5]-dpars1['ip'][5], par_in[5]+dpars1['ip'][5]]
        if hardbounds[0] < 0:
            hardbounds[0] = 0
        if hardbounds[3] < 0:
            hardbounds[3] = 1

        cycles = 2

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
                dpars = dpars1
            else:
                dpars = dpars2

            for optkind in optgroup:
                parfit_1 = optimizer( parstart, dpars[optkind], hardbounds, fitobj, optimize)
                parstart = parfit_1.copy()
                if args.debug:
                    outplotter_23(parfit_1, fitobj, '{}_{}_{}_parfit_{}{}'.format(order,night,tag,nc,optkind), trk, inparam, args, step2or3)
                    logger.debug(f'{order}_{tag}_{nc}_{optkind}:\n {parfit_1}')

        parfit = parfit_1.copy()
#-------------------------------------------------------------------------------
        # if stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            logger.warning(f'  --> parfit[1] < 0.1, {night} parfit={parfit}')
            continue
        # if stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
            logger.warning(f'  --> parfit[1] == par_in[1] or parfit[3] == par_in[3], {night}')
            continue
        # if model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq = fmod(parfit,fitobj)
        if len(smod[(smod < 0)]) > 0:
            logger.warning(f'  --> len(smod[(smod < 0)]) > 0, {night}')
            continue
#-------------------------------------------------------------------------------
        if args.plotfigs:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)

        rv0 = parfit[0] - parfit[2]  # atomosphere velocity correct

        rvsmini.append(rv0 + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2))
        vsinismini.append(parfit[4])

    bestguess = round(np.nanmean(rvsmini),5)
    vsinimini = round(np.nanmean(vsinismini),5)
    return night, bestguess, vsinimini

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline STEP 2',
                                     description = '''
                                     This is a pipeline that helps you to extract radial velocity \n
                                     from IGRINS spectra. \n
                                     Step 2 will help you get a best initial RV guess for starting \n
                                     Step 3.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name, no space",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input/UseWv/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))
    parser.add_argument("-l_use",   dest="label_use",        action="store",
                        help="Only one wavelength range will be used to RV initial guess, pick a label to use, Default is the first label",
                        type=int,   default=int(0))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Default = 50 ",
                        type=str,   default='50')

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s)",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",        action="store",
                        help="Range of allowed vsini variation during optimization (float, if set to zero vsini will be held constant), default = 5.0 km/s",
                        type=str,   default='5' )
    parser.add_argument('-g',       dest="guesses",          action="store",
                        help="Initial guesses for RV (int or float, km/s), or use -gX",
                        type=str,   default='')
    parser.add_argument('-gX',      dest="guessesX",         action="store",
                        help="Please give the number, X, under ./*targname/Initguesser_results_X that you wish to use as initial RV guesses",
                        type=str,   default='')

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
                        help="If you don't want all process all nights under the ./Input/*target folder, give an array of night you wish to process here. e.g., [20181111, 20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If sets, DeBug logging and extra plots will be given")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.85')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)
    #------------------------------
    if args.template not in ['synthetic', 'livingston', 'user_defined']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')
    #------------------------------
    if args.template in ['synthetic', 'livingston']:
        if args.sptype not in ['F','G','K','M']:
            sys.exit('ERROR: SPECTRAL TYPE FOR DEFAULT TEMPALTE ARE ONLY FOR F G K M! STARS')
    #------------------------------
    if (args.guesses != '') & (args.guessesX != ''):
        sys.exit('ERROR: YOU CAN ONLY CHOOSE EITHER -g OR -gX')
    #------------------------------
    if args.guesses != '':
        try:
            initguesses = float(args.guesses)
        except:
            sys.exit('ERROR: -g ONLY TAKE NUMBER AS INPUT!')
    #------------------------------
    if args.guessesX != '':
        try:
            guessdata = Table.read(f'./Output/{args.targname}_{args.band}/Initguesser_results_{args.guessesX}.csv', format='ascii')
        except:
            sys.exit(f'ERROR: "./Output/{args.targname}_{args.band}/Initguesser_results_{args.guessesX}.csv" NOT FOUND!')

        initnights = np.array(guessdata['night'])
        initrvs    = np.array(guessdata['bestguess'])
        initguesses = {}
        for hrt in range(len(initnights)):
            initguesses[str(initnights[hrt])] = [float(initrvs[hrt])]
#-------------------------------------------------------------------------------
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
    start_time = datetime.now()
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    \u00B1 {} km/s
    RV initial guess    = {} km/s
    '''.format(args.targname, initvsini, vsinivary, initguesses))
    print('---------------------------------------------------------------')
    logger.info('RV Initial Guess for {} Per Night...'.format(args.targname))
    print('This Will Take a While..........')
#-------------------------------------------------------------------------------
    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{iniguess_dir}')

    filew = open(f'./Output/{args.targname}_{args.band}/{iniguess_dir}','w')
    filew.write('night, bestguess, vsini')
    filew.write('\n')
#-------------------------------------------------------------------------------
    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders = read_prepdata(args)

    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = np.array(list(nightstemp))
        print('Only processing nights: {}'.format(nightsFinal))

    logger.info('Analyze with {} nights'.format(len(nightsFinal)))
#-------------------------------------------------------------------------------
    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, args.sptype)

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)
#-------------------------------------------------------------------------------
    pool = mp.Pool(processes = args.Nthreads)
    func = partial(ini_MPinst, args, inparam, orders, int(args.label_use), trk, step2or3 )
    outs = pool.map(func, np.arange(len(nightsFinal)))
    pool.close()
    pool.join()

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
    print('You can now try to get a better RV initial guess with by using -gX and rerun main_step2.py')
    print('OR, you can go on to the full RV extractor in main_step3.py')
    print('####################################################################################')
