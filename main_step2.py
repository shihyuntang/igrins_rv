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

def ini_MPinst(label_t, chunk_ind, trk, i):
    nights   = inparam.nights
    night    = nights[i]


    label = label_t[chunk_ind]
    order = label_t[chunk_ind]
    xbounds = inparam.xbounddict[label]

    print('Working on label {}, night {:03d}/{:03d} ({}) PID:{}...'.format(label,
                                                                    i+1,
                                                                    len(inparam.nights),
                                                                    night,
                                                                    mp.current_process().pid) )

#-------------------------------------------------------------------------------
    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]

    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == list:
        initguesses = inparam.initguesses
    else:
        print('EXPECTED FILE OR LIST FOR INITGUESSES! QUITTING!')

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    # Load telluric template from Telfit'd A0
    A0loc = './A0_Fits/A0_Fits_{}/{}A0_treated_{}.fits'.format(args.targname, night[:8], args.band)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('No A0-fitted template for night {}, skipping...'.format(night))
        print(A0loc)
        return night,np.nan,np.nan

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
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]

    if flag == 1:  # Telfit hit unknown critical error
        print('  --> TELFIT RESULT IS BAD, SKIP')
        return night,np.nan,np.nan

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

    ### Load relevant A0 spectra,
    rvcollect = []; vsinicollect = [];
    # Iterate over initial RV guesses
    for initrvguess in initguesses:
        rvsmini = []; vsinismini = [];
        # Iterate over all A/B exposures
        for t in np.arange(len(tagsnight)):
            tag = tagsnight[t]
            beam = beamsnight[t]

            if args.band=='K':
                if order==11:
                    bound_cut = [200, 100]
                elif order==12:
                    bound_cut = [900, 300]
                elif order==13:
                    bound_cut = [200, 400]
                elif order==14:
                    bound_cut = [150, 300]
                else:
                    bound_cut = [150, 150]
            elif args.band=='H':
                if order==10:
                    bound_cut = [250, 150]#ok
                elif order==11:
                    bound_cut = [600, 150]
                elif order==13:
                    bound_cut = [200, 600]#ok
                elif order==14:
                    bound_cut = [700, 100]
                elif order==16:
                    bound_cut = [400, 100]
                elif order==17:
                    bound_cut = [1000, 100]
                elif order==20:
                    bound_cut = [500, 150]
                elif (order==7) or (order==8) or (order==9) or (order==12) or (order==15) or (order==18) or (order==19):
                    bound_cut = [500, 500]
                else:
                    bound_cut = [150, 150]

            x,wave,s,u = init_fitsread('{}{}/{}/'.format(inparam.inpath, night, beam),
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
                print('  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
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
            par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

            par[0] = initrvguess-inparam.bvcs[night+tag]
            # Arrays defining parameter variations during optimization steps
            dpars1 = {'cont' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                     'wave' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0.,   0,   0, 0, 0,    0]),
                     't'    : np.array([0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     'ip'   : np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     's'    : np.array([20.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     'v'    : np.array([0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}
            dpars2 = {'cont' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0]),
                     'wave' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0.,   0,   0, 0, 0,    0]),
                     't'    : np.array([0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     'ip'   : np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     's'    : np.array([5.0, 2.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0]),
                     'v'    : np.array([0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])}

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

            optgroup = ['cont',
                        'wave',
                        't',
                        'cont',
                        's',
                        'cont',
                        'wave',
                        't',
                        's',
                        'cont',
                        'wave',
                        'ip',
                        'v',
                        'ip',
                        'v',
                        't',
                        's',
                        't',
                        's']

            nc = 1
            for cycle in range(cycles):

                if cycle == 0:
                    parstart = par_in.copy()
                    dpars = dpars1
                else:
                    dpars = dpars2

                for optkind in optgroup:
                    parfit_1 = optimizer(parstart,dpars[optkind],hardbounds,fitobj,optimize)
                    parstart = parfit_1.copy()
                    # print('{}: '.format(optkind), parstart)
                    if args.debug == True:
                        outplotter(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(label,night,tag,nc,optkind), trk, 1)
                    nc += 1

            parfit = parfit_1.copy()

            # if stellar template power is very low, throw out result
            if parfit[1] < 0.1:
                print('parfit[1] < 0.1, {} parfit={}'.format(night, parfit))
                continue

            # if stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
            if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
                print('parfit[1] == par_in[1] or parfit[3] == par_in[3], {}'.format(night))
                continue

            # if model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
            smod,chisq = fmod(parfit,fitobj)
            if len(smod[(smod < 0)]) > 0:
                print('len(smod[(smod < 0)]) > 0, {}'.format(night))
                continue

            if args.plotfigs == True:
                # outplotter(parfit, fitobj,'Post_parfit_{}_{}_{}'.format(label,night,tag), trk, 0)
                parfitS = parfit.copy(); parfitS[3] = 0
                parfitT = parfit.copy(); parfitT[1] = 0
                outplotter(parfitS, fitobj,'parfitS_{}_{}_{}'.format(label,night,tag), trk, 0)
                outplotter(parfitT, fitobj,'parfitT_{}_{}_{}'.format(label,night,tag), trk, 0)
                outplotter(parfit, fitobj,'parfit_{}_{}_{}'.format(label,night,tag), trk, 0)

            rv0 = parfit[0] - parfit[2]                         # atomosphere velocity correct

            rvsmini.append(rv0 + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2))
            vsinismini.append(parfit[4])

        rvcollect.append(np.nanmean(rvsmini))
        vsinicollect.append(np.nanmean(vsinismini))

    bestguess = round(np.nanmean(rvsmini),5)
    vsinimini = round(np.nanmean(vsinismini),5)
    return night,bestguess,vsinimini

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
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.5')
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
            guessdata = Table.read(f'./Output/{args.targname}_{args.band}/Initguesser_results_{args.band}.csv', format='ascii')
        except:
            sys.exit(f'ERROR: "./Output/{args.targname}_{args.band}/Initguesser_results_{args.band}.csv" NOT FOUND!')

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

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/main_step2_figs_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/main_step2_figs_{args.band}')

    outpath = f'./Output/{args.targname}_{args.band}/'
#-------------------------------------------------------------------------------
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s: %(message)s-->')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
#-------------------------------------------------------------------------------
    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{iniguess_dir}')

    filew = open(f'./Output/{args.targname}_{args.band}/{iniguess_dir}','w')
    filew.write('night, bestguess, vsini')
    filew.write('\n')
#-------------------------------------------------------------------------------
    start_time = datetime.now()
    print('###############################################################\n')
    print(r'''
Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    $\pm$ {} km/s
    RV initial guess    = {} km/s
    '''.format(args.targname, initvsini, vsinivary, initguesses))
    print('---------------------------------------------------------------')
    print('RV Initial Guess for {} Per Night...'.format(args.targname))
    print('This Will Take a While..........')

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, labels = read_prepdata(args)

    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=int)
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
    pool = mp.Pool(processes = args.Nthreads)
    func = partial(ini_MPinst, labels, int(args.label_use), trk )
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
    print('RV results:    mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(finalrvs),
                                                                                           np.nanmedian(finalrvs),
                                                                                           np.nanstd(finalrvs)      ))
    print('vsini results: mean= {:1.4f} km/s, median= {:1.4f} km/s, std= {:1.4f} km/s'.format(np.nanmean(vsinis),
                                                                                            np.nanmedian(vsinis),
                                                                                            np.nanstd(vsinis)      ))
    end_time = datetime.now()
    print('RV Initial Guess DONE... Duration: {}'.format(end_time - start_time))
    print('Output saved under {}_{}/{}'.format(args.targname, args.band, iniguess_dir) )
    print('---------------------------------------------------------------')
    # print('You can now try to get a better RV initial guess with: ')
    # print('(For RV standards) --> python main_step2.py {} -i {:1.1f} -v [input] -g [{:1.1f}] -c {} -plot'.format(args.targname,
    #                                                                                      np.nanmean(vsinis),
    #                                                                                      np.nanmean(finalrvs),
    #                                                                                      args.Nthreads))
    # print('(For other stars)  --> python main_step2.py {} -i {:1.1f} -v [input] -g {} -c {} -plot'.format(args.targname,
    #                                                                                      np.nanmean(vsinis),
    #                                                                                      trk-1,
    #                                                                                      args.Nthreads))
    # print('OR, you can go on to next analysis step with')
    # print('--> python main_step3tar.py {} -i {:1.1f} -v [input] -g IX{} -c {} -plot'.format(args.targname,
    #                                                                                   np.nanmean(vsinis),
    #                                                                                   trk-1,
    #                                                                                   args.Nthreads))
    # print('###############################################################')
    # print('\n')
