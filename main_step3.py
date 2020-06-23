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
def rv_MPinst(args, inparam, orders, order_use, trk, step2or3, i):
    nights   = inparam.nights
    night = nights[i]

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
    parfitminibox  = np.ones((len(tagsnight),15)); # need to match the dpar numbers

    rvsminibox[:]    = np.nan
    vsiniminibox[:]  = np.nan
    parfitminibox[:] = np.nan

    for t in tagsnight:
        nightsout.append(night)

    # Load telluric template from Telfit'd A0
    A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_treated_{args.band}.fits'
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

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
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

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

        # --------------------------------------------------------------
        par = pars0.copy()
        ##Set initial wavelength guess
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = inparam.initguesses-inparam.bvcs[night+tag]
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
                    outplotter_23(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(label,night,tag,nc,optkind), trk, inparam, args, step2or3)
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
        if args.plotfigs == True:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)
            outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3)

        rv0 = parfit[0] - parfit[2]                         # atomosphere velocity correct

        rvsminibox[t]   = rv0  + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2) # bvcs correct
        parfitminibox[t]= parfit
        vsiniminibox[t] = parfit[4]
    return nightsout,rvsminibox,parfitminibox,vsiniminibox

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
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Default = 50 ",
                        type=str,   default='50')
    parser.add_argument("-nAB",      dest="nAB",           action="store",
                        help="Minium request of # of AB sets. Default = 1",
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
    if (args.mode.lower() == 'std') & (args.guesses_source != ''):
        sys.exit('ERROR: STD CANNOT USE -gS, PLEASE USE -g')
    if (args.mode.lower() == 'tar') & (args.guesses != ''):
        sys.exit('ERROR: TAR CANNOT USE -g, PLEASE USE -gS')
    #------------------------------

    if args.mode.lower() == 'std':
        initguesses = float(args.guesses)
    else:
        if args.guesses_source == 'init':
            guesses = './Output/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
                                                                         args.band,
                                                                         int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre':
            guesses = './Output/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
                                                                      args.band,
                                                                      int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs    = np.array(guessdata['RVfinal'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])
#-------------------------------------------------------------------------------
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
    print('\n')
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
    logger.info('RV calculation for RV standard star {}...'.format(args.targname))
    print('This will take a while..........')
#-------------------------------------------------------------------------------
    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{name}')
#-------------------------------------------------------------------------------
    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders = read_prepdata(args)

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
    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T)
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
    for jerp in range(len(orders)): # Iterate over orders
        pool = mp.Pool(processes = args.Nthreads)
        func = partial(rv_MPinst, args, inparam, orders, jerp, trk, step2or3 )
        outs = pool.map(func, np.arange(len(nightsFinal)))
        pool.close()
        pool.join()

        order = orders[jerp]
#-------------------------------------------------------------------------------
        # Collect outputs
        for i in range(len(nights)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvbox     = outsbox[1]
                parfitbox = outsbox[2]
                vsinibox  = outsbox[3]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvbox     = np.concatenate((rvbox,outsbox[1]))
                parfitbox = np.vstack((parfitbox,outsbox[2]))
                vsinibox  = np.concatenate((vsinibox,outsbox[3]))

        nightsbox = np.array(nightsbox)
        vsinitags = []

        # Save results in fits file
        c1    = fits.Column(name='NIGHT'+str(order),  array=nightsbox, format='{}A'.format(len(nights[0])) )
        c2    = fits.Column(name='RV'+str(order),     array=rvbox,     format='D')
        c3    = fits.Column(name='PARFIT'+str(order), array=parfitbox, format=str(len(parfitbox[0,:]))+'D', dim=(1,len(parfitbox[0,:])))
        c4    = fits.Column(name='VSINI'+str(order),  array=vsinibox,  format='D')
        cols  = fits.ColDefs([c1,c2,c3,c4])
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
            # Iterating over nights, use these weights to combine the RVs from different chunks of one spectrum into a single RV
            for i in range(len(nights_use)):
                indnight  = np.where(nightsbox == nights_use[i])[0]
                rvtags    = rvbox[indnight]
                vsinitags = vsinibox[indnight]

                # Take the mean and std of the RVs determined from different A/B exposures within a night.
                if T_L == 'T':
                    vsinisT[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < int(args.nAB) ):
                        rvmasterboxT[i,jerp]  = np.nan
                        stdmasterboxT[i,jerp] = np.nan
                    else:
                        rvmasterboxT[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxT[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))

                else:
                    vsinisL[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < int(args.nAB) ):
                        rvmasterboxL[i,jerp]  = np.nan
                        stdmasterboxL[i,jerp] = np.nan
                    else:
                        rvmasterboxL[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxL[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))
            T_L = 'L'

#-------------------------------------------------------------------------------
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

    # Iterating over tight and loose mounting nights...
    for boxind in range(len(rvboxcomblist)):

        rvmasterbox  = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]
        vsinibox     = vsinicomblist[boxind]
#-------------------------------------------------------------------------------
        if args.mode=='STD':
            # Calculate the precision within an order across nights
            sigma_O2     = np.array([np.nanstd(rvmasterbox[:,ll])**2 for ll in range(len(orders))])
            sigma_ABbar2 = np.ones_like(sigma_O2)
            sigma_ON2    = np.ones_like(rvmasterbox)

            # Calculate uncertainty in method as difference between variance within an order and mean variance within a night's As and Bs RVs
            for ll in range(len(orders)):
                sigma_ABbar2[ll] = np.nanmean(stdmasterbox[:,ll]**2)
            sigma_method2 = sigma_O2 - sigma_ABbar2

        else:
            # Load the uncertainty from method from GJ 281 analysis
            if boxind == 0:
                nights_use = nightsT.copy()
                kind = 'Tight'
                sigma_method2 = inparam.methodvariance_tight[args.band]
            else:
                nights_use = nightsL.copy()
                kind = 'Loose'
                sigma_method2 = inparam.methodvariance_loose[args.band]
            # Calculate the uncertainty in each night/order RV as the sum of the uncertainty in method and the uncertainty in that night's As and Bs RVs
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

            if np.nansum(weights) == 0:
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan

            if np.sum( np.isnan(rvmasterbox[n,:]) ) > np.floor( len(orders) * 0.5 ):
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan
#-------------------------------------------------------------------------------

        f = plt.figure(figsize=(5,3))
        ax1 = plt.subplot(111)
        ax1.plot(    np.arange(len(rvfinal))+1, rvfinal, '.k', ms=5)
        ax1.errorbar(np.arange(len(rvfinal))+1, rvfinal, yerr=stdfinal, ls='none', lw=.5, ecolor='black')
        ax1.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinal), np.nanstd(rvfinal)),
                             transform=ax1.transAxes)
        ax1.set_ylim(np.nanmin(rvfinal)-.08,
                     np.nanmax(rvfinal)+.08)
        ax1.set_ylabel('RV [km/s]')
        ax1.set_xlabel('Night (#)')
        ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in', width=.6)
        f.savefig('{}/{}/FinalRVs_{}_.png'.format(inparam.outpath, name, kind), format='png', bbox_inches='tight')

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
            cols  = fits.ColDefs([c1,c2,c3,c4,c7,c8,c9,c10])
        else:
            cols  = fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])

        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul        = fits.HDUList([primary_hdu,hdu_1])
        hdul.writeto('{}/{}/RVresultsSummary_{}.fits'.format(inparam.outpath, name, kind), overwrite=True)

        nightsCombined     = np.concatenate((nightsCombined,     nights_use))
        mjdsCombined       = np.concatenate((mjdsCombined,       mjds_out))
        rvfinalCombined    = np.concatenate((rvfinalCombined,    rvfinal))
        stdfinalCombined   = np.concatenate((stdfinalCombined,   stdfinal))
        vsinifinalCombined = np.concatenate((vsinifinalCombined, vsinifinal))

        if args.mode=='STD':
            sigma_method2 = [np.around(float(i), 8) for i in sigma_method2]
            logger.info('sigma_method2 with type = {} is {}'.format(kind, sigma_method2))
        logger.info('Observations when IGRINS is mounting {}: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( kind,
                                                                                                            np.nanmean(rvfinal),
                                                                                                            np.nanstd(rvfinal) ))
#-------------------------------------------------------------------------------
    xscale = np.arange(len(rvfinalCombined))+1

    f = plt.figure(figsize=(5,3))
    ax1 = plt.subplot(111)
    ax1.plot(xscale,rvfinalCombined, '.k', ms=5)
    ax1.errorbar(xscale,rvfinalCombined,yerr=stdfinalCombined,ls='none',lw=.5, ecolor='black')
    ax1.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=ax1.transAxes)

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        ax1.text(0.05, 0.1, 'Loose', transform=ax1.transAxes)
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            ax1.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
            ax1.text(0.9,  0.1, 'Loose', transform=ax1.transAxes)
        else:
            ax1.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
            ax1.text(0.9,  0.1, 'Loose', transform=ax1.transAxes)
    ax1.set_ylim(np.nanmin(rvfinalCombined)-.08,np.nanmax(rvfinalCombined)+.08)
    ax1.set_ylabel('RV (km/s)')
    ax1.set_xlabel('Night (#)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/{}/FinalRVs.png'.format(inparam.outpath, name), format='png', bbox_inches='tight')


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
    logger.info('Output saved under {}/{}'.format(args.outpath, name) )
    print('####################################################################################')
