from Engine.importmodule import *

from Engine.IO_AB import setup_templates_syn, setup_templates_sun, init_fitsread, stellarmodel_setup, setup_outdir, setup_templates
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs, inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer, fmod
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def outplotter(parfit,fitobj,title,trk,debug):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(w,fitobj.s, '-k',  lw=0.5, label='data',  alpha=.6)
    axes.plot(w,fit,      '--r', lw=0.5, label='model',  alpha=.6)

    axes.set_title( title,                 size=5, style='normal', family='sans-serif')
    axes.set_ylabel(r'Normalized Flux',    size=5, style='normal', family='sans-serif')
    axes.set_xlabel(r'Wavelength [$\AA$]', size=5, style='normal', family='sans-serif')

    axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
    axes.legend(fontsize=4, edgecolor='white')
    if debug == 0:
        fig.savefig('{}/figs/main_step2_{}/{}.png'.format(inparam.outpath, trk, title), bbox_inches='tight', format='png', overwrite=True)
    elif debug == 1:
        fig.savefig('./Temp/Debug/{}_{}/main_step2_{}/{}.png'.format(args.targname, args.band, trk, title), bbox_inches='tight', format='png', overwrite=True)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def ini_MPinst(label_t, chunk_ind, trk, i):
    nights   = inparam.nights
    night    = nights[i]

    label = '{}-{}'.format( label_t['0'][chunk_ind], label_t['1'][chunk_ind] )
    order = label_t['0'][chunk_ind]
    chunk = label_t['1'][chunk_ind]
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
        return night,np.nan,np.nan

    num_orders = len( np.unique(label_t['0']) )

    # order in A0_treated.fits is no longer sequential...
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]

    if flag == 1:  # Telfit hit unknown critical error
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
                      IPpars[0],                                            #14: IP quadratic component
                      0.5])                                                 #15: Differential Rotation Coefficient

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
                    bound_cut = [150, 100]
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
                    bound_cut = [150, 100]

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
            dpars = {'cont' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0, 0]),
                     'wave' : np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 0.,   0,   0, 0, 0,    0, 0]),
                     't'    : np.array([0.0, 0.0, 5.0, 1.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0]),
                     'ip'   : np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0]),
                     's'    : np.array([5.0, 5.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0]),
                     'v'    : np.array([0.0, 0.0, 0.0, 0.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0])}

            continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
            s_piece /= np.median(s_piece)
            fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in)
        #-------------------------------------------------------------------------------
                ######## Begin optimization  ########

            optimize = True
            par_in = par.copy()
            # hardbounds = [par_in[4] -dpar[4],   par_in[4]+dpar[4],
            #               par_in[5] -dpar[5],   par_in[5]+dpar[5],
            #               par_in[15]-dpar[15], par_in[15]+dpar[15]]
            hardbounds = [par_in[4]-dpars['v'][4],  par_in[4]+dpars['v'][4],
                          par_in[5]-dpars['ip'][5], par_in[5]+dpars['ip'][5],
                          par_in[15]-dpars['s'][15], par_in[15]+dpars['s'][15]]
            if hardbounds[0] < 0:
                hardbounds[0] = 0
            if hardbounds[3] < 0:
                hardbounds[3] = 1
            if hardbounds[4] < 0.1:
                hardbounds[4] = 0.1
            if hardbounds[5] > 0.9:
                hardbounds[5] = 0.9

            cycles = 4

            optgroup = ['cont',
                        'wave',
                        't',
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

                for optkind in optgroup:
                    parfit_1 = optimizer(parstart,dpars[optkind],hardbounds,fitobj,optimize)
                    parstart = parfit_1.copy()
                    if args.debug == True:
                        outplotter(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(label,night,tag,nc,optkind), trk, 1)
                    nc += 1

            parfit = parfit_1.copy()

            if args.plotfigs == True:
                # outplotter(parfit, fitobj,'Post_parfit_{}_{}_{}'.format(label,night,tag), trk, 0)
                parfitS = parfit.copy(); parfitS[3] = 0
                parfitT = parfit.copy(); parfitT[1] = 0
                outplotter(parfitS, fitobj,'parfitS_{}_{}_{}'.format(label,night,tag), trk, 0)
                outplotter(parfitT, fitobj,'parfitT_{}_{}_{}'.format(label,night,tag), trk, 0)
                outplotter(parfit, fitobj,'parfit_{}_{}_{}'.format(label,night,tag), trk, 0)


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

            #
            # if args.debug == True:
            #     outplotter(parfit_1,fitobj,'Post_parfit_1_{}_{}_{}'.format(label,night,tag), trk, 1)
            #     outplotter(parfit_2,fitobj,'Post_parfit_2_{}_{}_{}'.format(label,night,tag), trk, 1)
            #     outplotter(parfit_3,fitobj,'Post_parfit_3_{}_{}_{}'.format(label,night,tag), trk, 1)
            #     outplotter(parfit_4,fitobj,'Post_parfit_4_{}_{}_{}'.format(label,night,tag), trk, 1)
            #     outplotter(parfit  ,fitobj,'Post_parfit_{}_{}_{}'.format(label,night,tag), trk, 1)

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
                        help="Enter your *target name",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))
    parser.add_argument("-l_use",   dest="label_use",        action="store",
                        help="Only one wavelength range will be used to RV initial guess, pick a label to use, Default is the first label",
                        type=int,   default=int(0))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Default = 50 ",
                        type=str,   default='50')

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s), default = 2.6 km/s",
                        type=str,   default='2.6' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization (float, if set to zero vsini will be held constant), default = 0.5 km/s",
                        type=str,   default='0.5' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help="Initial guesses for RV. May either provide a list (e.g., [0] or [20,25,30], !no space!) that will be uniformly applied to every night, or provide a path stores the previous Initguesses result file",
                        type=str,   default='20.02' )

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If sets, will generate plots")

    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want all process all nights under the Input_Data folder, give an array of night you wish to process here. e.g., [20181111, 20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If sets, will generate files and plots under ./Temp/Debug for debug")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()
    cdbs_loc = '~/cdbs/'
    inpath     = './Input_Data/{}/'.format(args.targname)
    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)
    guesses   = args.guesses


    if guesses[0]=='[':
        initguesses = ast.literal_eval(guesses) #convert str(list) to list
    elif '/' in guesses:
        guessdata  = Table.read(guesses,format='ascii')
        initnights = np.array(guessdata['night'])
        initrvs    = np.array(guessdata['bestguess'])
        initguesses = {}
        for hrt in range(len(initnights)):
            initguesses[str(initnights[hrt])] = [float(initrvs[hrt])]
    else:
        sys.exit('ERROR: INCORRECT "STYLE" INPUT PARAMETER SPECIFIED; EXPECTED EITHER "[list]" OR "file dir"')

#-------------------------------------------------------------------------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('''
Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    = {} km/s
    RV initial guess    = {} km/s
    '''.format(args.targname, initvsini, vsinivary, initguesses))
    print('---------------------------------------------------------------')
    print('RV Initial Guess for {} Per Night...'.format(args.targname))
    print('This Will Take a While..........')

    ## Collect relevant file information from Predata files
    A0data   = Table.read('./Temp/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), format='ascii')
    A0nights = np.array(A0data['night'],dtype='str')
    ams0     = np.array(A0data['airmass'])

    targdata = Table.read('./Temp/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), format='ascii')
    Tnights = np.array(targdata['night'],dtype='str')
    tags0   = np.array(targdata['tag'], dtype='int')
    beams0  = np.array(targdata['beam'],dtype='str')
    mjds0   = np.array(targdata['mjd'])
    bvcs0   = np.array(targdata['bvc'])
    ams     = np.array(targdata['airmass'])

    bounddata = Table.read('./Input_Data/Use_w/XRegions_{}_{}.csv'.format(args.WRegion, args.band), format='csv')
    starts  = np.array(bounddata['start'])
    ends    = np.array(bounddata['end'])
    labels  = np.array(bounddata['label'], dtype=str)
    xbounddict = {labels[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}

    # Attribute A and B exposures to right file numbers
    tagsA = {}; tagsB = {}; mjds = {}; bvcs = {};
    night_orig = Tnights[0]; tagsA0 = []; tagsB0 = [];

    for hrt in range(len(Tnights)):
        tag1 = '{:04d}'.format(tags0[hrt])

        mjds[Tnights[hrt]]                = float(mjds0[hrt])
        bvcs[str(Tnights[hrt])+str(tag1)] = float(bvcs0[hrt])

        if Tnights[hrt] == night_orig:
            if beams0[hrt] == 'A':
                tagsA0.append(tag1)
            else:
                tagsB0.append(tag1)
        else:
            tagsA[Tnights[hrt-1]] = tagsA0
            tagsB[Tnights[hrt-1]] = tagsB0
            tagsA0 = []; tagsB0 = [];
            if beams0[hrt] == 'A':
                tagsA0.append(tag1)
            else:
                tagsB0.append(tag1)
            night_orig = Tnights[hrt].copy()

    tagsA[Tnights[-1]] = tagsA0
    tagsB[Tnights[-1]] = tagsB0

    nightsFinal = np.array(list(sorted(set(Tnights))))
    # nightsFinal = nightsFinal[24:45]

    if args.nights_use != '':
        nightstemp = np.array(args.nights_use, dtype=np.int)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))
#-------------------------------------------------------------------------------
    if not os.path.isdir('./Results/'):
        os.mkdir('./Results/')

    # Create output directory
    try:
        filesndirs = os.listdir('./Results/{}_{}'.format(args.targname, args.band) )
    except:
        os.mkdir('./Results/{}_{}'.format(args.targname, args.band))
        filesndirs = os.listdir( './Results/{}_{}'.format(args.targname, args.band) )
    trk = 1; go = True;
    while go == True:
        iniguess_dir = 'Initguesser_results_{}.csv'.format(trk)
        if iniguess_dir not in filesndirs:
            break
        trk += 1

    if args.debug:
        try:
            os.listdir('./Temp/Debug/{}_{}/main_step2_{}/'.format(args.targname, args.band), trk)
        except OSError:
            os.mkdir('./Temp/Debug/{}_{}/main_step2_{}/'.format(args.targname, args.band), trk)
#-------------------------------------------------------------------------------
    print('Writing output to ./Results/{}_{}/{}'.format(args.targname, args.band, iniguess_dir))
    filew = open('./Results/{}_{}/{}'.format(args.targname, args.band, iniguess_dir),'w')
    filew.write('night, bestguess, vsini')
    filew.write('\n')

    if not os.path.isdir('./Results/{}_{}/figs'.format(args.targname, args.band)):
        os.mkdir('./Results/{}_{}/figs'.format(args.targname, args.band) )

    if not os.path.isdir('./Results/{}_{}/figs/main_step2_{}'.format(args.targname, args.band, trk)):
        os.mkdir('./Results/{}_{}/figs/main_step2_{}'.format(args.targname, args.band, trk) )
    outpath = './Results/{}_{}'.format(args.targname, args.band)
#-------------------------------------------------------------------------------
    # Retrieve stellar and telluric templates

    if (args.targname == 'TauBoo') | (args.targname == 'HD26257'):
        print('Using: SpotAtl_Solar')
        watm,satm, mwave0, mflux0 = setup_templates_sun()
    else:
        if args.band=='K':
            watm,satm, mwave0, mflux0 = setup_templates_syn()
            print('Using: syntheticstellar_kband')
        elif args.band=='H':
            watm,satm, mwave0, mflux0 = setup_templates()
            print('Using: SpotAtl Organized')

    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict)

    # Only use first wavelength region listed
    ### label = labels[0] IF ONLY RV STANDARD, SPECIFY OTHERWISE. OR USE METHOD2

    orders = [ int(labels[i].split('-')[0]) for i in range(len(labels)) ]
    oindex = [ int(labels[i].split('-')[1]) for i in range(len(labels)) ]
    label_t = Table(names=('0', '1'), data=(orders, oindex))
    label_t.sort(['0', '1'])
#-------------------------------------------------------------------------------
    pool = mp.Pool(processes = args.Nthreads)
    func = partial(ini_MPinst, label_t, int(args.label_use), trk )
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
