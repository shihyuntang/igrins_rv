from Engine.importmodule import *

from Engine.IO_AB import setup_templates, setup_templates_syn, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs,inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer, fmod
#-------------------------------------------------------------------------------
def outplotter(parfit,fitobj,title,trk,debug):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(w,fitobj.s, '-k',  lw=0.5, label='data',alpha=.6)
    axes.plot(w,fit,      '--r', lw=0.5, label='model',alpha=.6)

    axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
    axes.set_title(title,   size=5, style='normal' , family='sans-serif' )
    axes.set_ylabel(r'Normalized Flux',   size=5, style='normal' , family='sans-serif' )
    axes.set_xlabel(r'Wavelength [$\AA$]'      size=5, style='normal' , family='sans-serif' )
    axes.legend(fontsize=4, edgecolor='white')
    if debug == 0:
        fig.savefig('{}/figs/main_step3_{}/{}.png'.format(inparam.outpath, trk, title), bbox_inches='tight', format='png', overwrite=True)
    elif debug == 1:
        fig.savefig('./Temp/Debug/{}_{}/main_step3_{}/{}.png'.format(args.targname, args.band, trk, title), bbox_inches='tight', format='png', overwrite=True)


#-------------------------------------------------------------------------------
def rv_MPinst(label_t, chunk_ind, trk, i):
    # Main function to compute RVs for a given night and order
    nights   = inparam.nights
    night = nights[i]

    label = '{}-{}'.format( label_t['0'][chunk_ind], label_t['1'][chunk_ind] )
    order = label_t['0'][chunk_ind]
    chunk = label_t['1'][chunk_ind]

    xbounds = inparam.xbounddict[label]

    print('Working on label {:03d}/{:03d} ({}), night {:03d}/{:03d} ({}) PID:{}...'.format(int(chunk_ind)+1,
                                                                                    len(label_t),
                                                                                    label,
                                                                                    i+1,
                                                                                    len(inparam.nights),
                                                                                    night,
                                                                                    mp.current_process().pid) )

    # --------------------------------------------------------------
    # --------------------------------------------------------------
    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night) < 20180401 or int(night) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]

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
    parfitminibox  = np.ones((len(tagsnight),16));

    rvsminibox[:]    = np.nan
    vsiniminibox[:]  = np.nan
    parfitminibox[:] = np.nan
    
    # Load telluric template from Telfit'd A0
    curdir = os.getcwd()
    A0loc = './A0_Fits/A0_Fits_{}/{}A0_treated_{}.fits'.format(args.targname, night, args.band)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('  --> No A0-fitted template for night '+night+', skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

    num_orders = len( np.unique(label_t['0']) )

    # order in A0_treated.fits is no longer sequential...
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]

    if flag == 1:  # Telfit hit unknown critical error
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

    try:
        if np.isnan(inparam.initguesses[night]):  # Telfit hit unknown critical error
            print('  --> Initial guess for {} is NaN , SKIP...'.format(night))
            return nightsout, rvsminibox, parfitminibox, vsiniminibox
    except:
        print('  --> Initial guess for {} is NaN , SKIP...'.format(night))
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

    watm = tbdata['WATM'+str(order)]
    satm = tbdata['SATM'+str(order)]
    a0contx = tbdata['X'+str(order)]
    continuum  = tbdata['BLAZE'+str(order)]

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]
    satm[(satm < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
    a0contx = a0contx[(continuum != 0)]
    continuum = continuum[(continuum != 0)]

    # --------------------------------------------------------------
    # --------------------------------------------------------------
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
                      IPpars[0],                                              #14: IP quadratic component
                      0.675])                                                #15: Differential Rotation Coefficient

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

        try:
            scomb = np.vstack((scomb,rebin_jv(x_piece,s_piece,xcomb[0,:],False)))
            ucomb = np.vstack((ucomb,rebin_jv(x_piece,u_piece,xcomb[0,:],False)))
            wcomb = np.vstack((wcomb,rebin_jv(x_piece,wave_piece,xcomb[0,:],False)))
        except UnboundLocalError:
            scomb = s_piece.copy()
            ucomb = u_piece.copy()
            wcomb = wave_piece.copy()
            xcomb = x_piece.copy()

    s_piece = np.array([np.nanmean(scomb[:,pip]) for pip in range(len(scomb[0,:]))])
    wave_piece = np.array([np.nanmean(wcomb[:,pip]) for pip in range(len(scomb[0,:]))])
    u_piece = np.array([1/np.nansqrt(np.nansum(1/(ucomb[:,pip]**2))) for pip in range(len(scomb[0,:]))])

    mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

    satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
    watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

    # Cut target spec to be within A0 spec wave
    s_piece    = s_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    u_piece    = u_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    x_piece    = x_piece[   (wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
    wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

    # --------------------------------------------------------------
    # Load initial parameters, assign inital RV guess
    par = pars0.copy()
    f = np.polyfit(x_piece,wave_piece,3)
    par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
    par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

    par[0] = inparam.initguesses[night]-inparam.bvcs[night+tag]
    # Arrays defining parameter variations during optimization steps
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0, 0])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0, 0])
    dpar      = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0])
    dpar_st   = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0])
    dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0])
    
    continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
    s_piece /= np.median(s_piece)
    fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in)
#-------------------------------------------------------------------------------
    ######## Begin optimization  ########

    optimize = True
    par_in = par.copy()
    hardbounds = [par_in[4]-dpar[4],par_in[4]+dpar[4],par_in[5]-dpar[5],par_in[5]+dpar[5]]
    if hardbounds[0] < 0:
        hardbounds[0] = 0
    if hardbounds[3] < 0:
        hardbounds[3] = 1

#        if args.plotfigs == True:#
#            outplotter(targname,par_in,fitobj,'{}_{}_{}_1'.format(label,night,tag))

    parfit_1 = optimizer(par_in,   dpar_cont, hardbounds,fitobj,optimize)
    parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds,fitobj,optimize)
    parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds,fitobj,optimize)
    parfit_4 = optimizer(parfit_3, dpar_wave, hardbounds,fitobj,optimize)
    parfit = optimizer(parfit_4,   dpar,      hardbounds,fitobj,optimize)   # RV fitting

    # if stellar template power is very low, throw out result
    if parfit[1] < 0.1:
        continue

    # if stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
    if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
        continue

    # if model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
    smod,chisq = fmod(parfit,fitobj)
    if len(smod[(smod < 0)]) > 0:
        continue

    if args.plotfigs == True:
        outplotter(parfit, fitobj,'Post_parfit_{}_{}_{}'.format(label,night,tag), trk, 0)

    if args.debug == True:
        outplotter(parfit_1,fitobj,'Post_parfit_1_{}_{}_{}'.format(label,night,tag), trk, 1)
        outplotter(parfit_2,fitobj,'Post_parfit_2_{}_{}_{}'.format(label,night,tag), trk, 1)
        outplotter(parfit_3,fitobj,'Post_parfit_3_{}_{}_{}'.format(label,night,tag), trk, 1)
        outplotter(parfit_4,fitobj,'Post_parfit_4_{}_{}_{}'.format(label,night,tag), trk, 1)
        outplotter(parfit  ,fitobj,'Post_parfit_{}_{}_{}'.format(label,night,tag), trk, 1)

    rv0 = parfit[0] - parfit[2]                         # atomosphere velocity correct

    rvout   = rv0  + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2) # bvcs correct

    return night,rvout, parfit, parfit[4]

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
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Default = 50 ",
                        type=str,   default='50')

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). Should use the value given by step2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 0.0 km/s",
                        type=str, default='0.0' )
    parser.add_argument('-gS',       dest="guesses_source",           action="store",
                        help="Source for initial guesses list for RV. Enter init OR rvre (init: Initguesser_results_X, rvre: RV_results_X)",
                        type=str, default='')
    parser.add_argument('-gX',       dest="guesses",           action="store",
                        help="Please give the number, X, under ./*targname/Initguesser_results_X OR ./*targname/RV_results_X, that you wish to use",
                        type=int, default='')

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
    vsinivary = float(args.vsinivary)

    if args.initvsini != '':
        initvsini = float(args.initvsini)
    else:
        sys.exit('ERROR: EXPECTED FLOAT')

    # Collect init RV guesses
    if args.guesses != '':
        if args.guesses_source == 'init':
            guesses = './Results/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
                                                                          args.band,
                                                                          int(args.guesses))
            guessdata  = Table.read(guesses, format='ascii')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre':
            guesses = './Results/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
                                                                       args.band,
                                                                       int(args.guesses))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs    = np.array(guessdata['RVfinal'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])
    else:
        sys.exit('ERROR: INCORRECT "STYLE" INPUT PARAMETER SPECIFIED; EXPECTED A INT NUMBER')
#-------------------------------------------------------------------------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('''
    Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    = {} km/s
    RV initial guess taken from {}
    '''.format(args.targname, initvsini, vsinivary, guesses))
    print('---------------------------------------------------------------')
    print('RV calculation for target star {}...'.format(args.targname))
    print('This will take a while..........')

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

        mjds[Tnights[hrt]] = float(mjds0[hrt])
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
    if args.nights_use != '':
        nightstemp = np.array(args.nights_use, dtype=np.int)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))
#-------------------------------------------------------------------------------

    # Create output directory
    try:
        filesndirs = os.listdir('./Results/{}_{}_/'.format(args.targname, args.band) )
    except OSError:
        os.mkdir('./Results/{}_{}'.format(args.targname, args.band))
        filesndirs = os.listdir( './Results/{}_{}'.format(args.targname, args.band) )
    trk = 1; go = True;
    while go == True:
        name = 'RV_results_ABCombined_'+str(trk)
        if name not in filesndirs:
            break
        trk += 1
    os.mkdir('./Results/{}_{}/{}'.format(args.targname, args.band, name) )
#-------------------------------------------------------------------------------
    print('Writing output to folder ./Results/{}_{}'.format(args.targname, args.band, name))

    if not os.path.isdir('./Results/{}_{}/figs'.format(args.targname, args.band)):
        os.mkdir('./Results/{}_{}/figs'.format(args.targname, args.band) )

    if not os.path.isdir('./Results/{}_{}/figs/main_step3_{}'.format(args.targname, args.band, trk)):
        os.mkdir('./Results/{}_{}/figs/main_step3_{}'.format(args.targname, args.band, trk) )
    outpath = './Results/{}_{}'.format(args.targname, args.band)
#-------------------------------------------------------------------------------
    # Retrieve stellar and telluric templates
    if args.band=='K':
        watm,satm, mwave0, mflux0 = setup_templates_syn()
    elif args.band=='H':
        watm,satm, mwave0, mflux0 = setup_templates()

    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict)

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T)
    nights    = inparam.nights
    intnights = nights.astype(int)

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]

    rvmasterboxT  = np.ones((len(nightsT),len(labels)))
    stdmasterboxT = np.ones((len(nightsT),len(labels)))
    rvmasterboxL  = np.ones((len(nightsL),len(labels)))
    stdmasterboxL = np.ones((len(nightsL),len(labels)))
    vsinisT = np.ones((len(nightsT),len(labels)))
    vsinisL  = np.ones((len(nightsL),len(labels)))

    if len(nightsL) > 0:
        nightscomblist = [nightsT,nightsL]
    else:
        nightscomblist = [nightsT]

    orders = [ int(labels[i].split('-')[0]) for i in range(len(labels)) ]
    oindex = [ int(labels[i].split('-')[1]) for i in range(len(labels)) ]
    label_t = Table(names=('0', '1'), data=(orders, oindex))
    label_t.sort(['0', '1'])
#-------------------------------------------------------------------------------
    for jerp in range(len(label_t)): # Iterate over orders
        pool = mp.Pool(processes = args.Nthreads)
        func = partial(rv_MPinst, label_t, jerp, trk )
        outs = pool.map(func, np.arange(len(nightsFinal)))
        pool.close()
        pool.join()

        label = '{}-{}'.format( label_t['0'][jerp], label_t['1'][jerp] )
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
        c1    = fits.Column(name='NIGHT'+str(label),  array=nightsbox, format='8A')
        c2    = fits.Column(name='RV'+str(label),     array=rvbox,     format='D')
        c3    = fits.Column(name='PARFIT'+str(label), array=parfitbox, format=str(len(parfitbox[0,:]))+'D', dim=(1,len(parfitbox[0,:])))
        c4    = fits.Column(name='VSINI'+str(label),  array=vsinibox,  format='D')
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
