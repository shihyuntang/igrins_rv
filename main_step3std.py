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
    axes.set_xlabel('Wavelength',       size=5, style='normal' , family='sans-serif' )
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
    parfitminibox  = np.ones((len(tagsnight),15));

    rvsminibox[:]    = np.nan
    vsiniminibox[:]  = np.nan
    parfitminibox[:] = np.nan

    for t in tagsnight:
        nightsout.append(night)

    # Load telluric template from Telfit'd A0
    A0loc = './A0_Fits/A0_Fits_{}/{}A0_treated_{}.fits'.format(args.targname, night, args.band)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('  --> No A0-fitted template for night '+night+', skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

    print(label_t)
    print(label_t['0'])
    print(np.unique(label_t['0']))
    print(len( np.unique(label_t['0']) ))

    num_orders = len( np.unique(label_t['0']) )
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[3].name[1:]) == int(order) ][0]

    try:
        fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[3].name[1:]) == int(order) ][0]
        # same as flag == 1
        # order in A0_treated.fits is no longer sequential...
    except:
        print(A0loc)
        print('  --> {} nights, fits_layer locater ERROR, {} not match order: {}'.format(night, [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[3].name[1:]) == order ],  order))
        return nightsout, rvsminibox, parfitminibox, vsiniminibox

    tbdata = hdulist[ fits_layer ].data


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
                      IPpars[0]])                                            #14: IP quadratic component

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
        SN_cut = 25
        if np.nanmedian(s2n) < SN_cut: # If S/N less than 25, throw out
            print('  --> Bad S/N {:1.1f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), SN_cut, night, beam, tag))
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
        # Load initial parameters, assign inital RV guess
        par = pars0.copy()
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = inparam.initguesses-inparam.bvcs[night+tag]
        # Arrays defining parameter variations during optimization steps
        dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0])
        dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0])
        dpar      = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0])
        dpar_st   = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])
        dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        s_piece /= np.median(s_piece)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in)
#-------------------------------------------------------------------------------
        ######## Begin optimization  ########

        optimize = True
        par_in = par.copy()

#        if args.plotfigs == True:#
#            outplotter(targname,par_in,fitobj,'{}_{}_{}_1'.format(label,night,tag))

        parfit_1 = optimizer(par_in,dpar_cont,fitobj,optimize)
        parfit_2 = optimizer(parfit_1,dpar_st,fitobj,optimize)
        parfit_3 = optimizer(parfit_2,dpar_wave,fitobj,optimize)
        parfit_4 = optimizer(parfit_3,dpar_cont,fitobj,optimize)
        #parfit_5 = optimizer(parfit_4,dpar_ip,fitobj,optimize)
        parfit = optimizer(parfit_4,dpar,fitobj,optimize)   # RV fitting

        if args.plotfigs == True:
            #outplotter(par_in, fitobj,'{}_{}_{}_par_in'.format(label,night,tag), trk, 0)
            outplotter(parfit, fitobj,'{}_{}_{}_parfit'.format(label,night,tag), trk, 0)

        if args.debug == True:
            outplotter(parfit_1,fitobj,'{}_{}_{}_parfit_1'.format(label,night,tag), trk, 1)
            outplotter(parfit_2,fitobj,'{}_{}_{}_parfit_2'.format(label,night,tag), trk, 1)
            outplotter(parfit_3,fitobj,'{}_{}_{}_parfit_3'.format(label,night,tag), trk, 1)
            outplotter(parfit_4,fitobj,'{}_{}_{}_parfit_4'.format(label,night,tag), trk, 1)
            outplotter(parfit_5,fitobj,'{}_{}_{}_parfit_5'.format(label,night,tag), trk, 1)
            outplotter(parfit  ,fitobj,'{}_{}_{}_parfit'.format(label,night,tag), trk, 1)

        rv0 = parfit[0] - parfit[2]                         # atomosphere velocity correct

        rvsminibox[t]   = rv0  + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(3e5**2) # bvcs correct
        parfitminibox[t]= parfit
        vsiniminibox[t] = parfit[4]
    print(nightsout,rvsminibox,parfitminibox,vsiniminibox)
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
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). Should use the value given by step2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 0.0 km/s",
                        type=str, default='0.0' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help=". Should use the single value given by step2 (float, km/s)",
                        type=str,   default='' )
    # parser.add_argument('-gS',       dest="guesses_source",           action="store",
    #                     help="Source for initial guesses list for RV. Enter init OR rvre (init: Initguesser_results_X, rvre: RV_results_X)",
    #                     type=str, default='')
    # parser.add_argument('-gX',       dest="guesses",           action="store",
    #                     help="Please give the number, X, under ./*targname/Initguesser_results_X OR ./*targname/RV_results_X, that you wish to use",
    #                     type=int, default='')

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

    if args.guesses != '':
        guesses = float(args.guesses)
    else:
        sys.exit('ERROR: EXPECTED FLOAT')

    if type(guesses) == float:
        initguesses = guesses
    else:
        sys.exit('ERROR: INCORRECT INITIAL RV GUESSES INPUT PARAMETER SPECIFIED; EXPECTED FLOAT')
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
    print('RV calculation for RV standard star {}...'.format(args.targname))
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
    nightsFinal = nightsFinal[24:45]
    labels      = labels[-2:]

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
        filesndirs = os.listdir('./Results/{}_{}/'.format(args.targname, args.band) )
    except OSError:
        os.mkdir('./Results/{}_{}'.format(args.targname, args.band))
        filesndirs = os.listdir( './Results/{}_{}'.format(args.targname, args.band) )
    trk = 1; go = True;
    while go == True:
        name = 'RV_results_'+str(trk)
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

    print(labels)
    print([ labels[i].split('-')[0] for i in range(len(labels)) ])
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

                    if (len(rvtags) == 0) or (len(rvtags) == 1):
                        rvmasterboxT[i,jerp]  = np.nan
                        stdmasterboxT[i,jerp] = np.nan
                    else:
                        rvmasterboxT[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxT[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))
                else:
                    vsinisL[i,jerp] = np.nanmean(vsinitags)

                    if (len(rvtags) == 0) or (len(rvtags) == 1):
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
        # Calculate the precision within an order across nights
        sigma_O2     = np.array([np.nanstd(rvmasterbox[:,ll])**2 for ll in range(len(labels))])
        sigma_ABbar2 = np.ones_like(sigma_O2)
        sigma_ON2    = np.ones_like(rvmasterbox)

#-------------------------------------------------------------------------------
        # Note rvmasterbox indexed as [nights,orders]
        Nnights = len(rvmasterbox[:,0])

        # Calculate uncertainty in method as difference between variance within an order and mean variance within a night's As and Bs RVs
        for ll in range(len(labels)):
            sigma_ABbar2[ll] = np.nanmean(stdmasterbox[:,ll]**2)
        sigma_method2 = sigma_O2 - sigma_ABbar2

        # Calculate the uncertainty in each night/order RV as the sum of the uncertainty in method and the uncertainty in that night's As and Bs RVs
        for ll in range(len(labels)):
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

        c1 = fits.Column( name='NIGHT',         array=nights_use,        format='8A')
        c2 = fits.Column( name='MJD',           array=mjds_out-2400000.5,format='D')
        c3 = fits.Column( name='RVBOX',         array=rvmasterbox,   format='{}D'.format(len(label_t)))
        c4 = fits.Column( name='STDBOX',        array=stdmasterbox,  format='{}D'.format(len(label_t)))
        c5 = fits.Column( name='Sigma_O2',      array=sigma_O2,      format='D')
        c6 = fits.Column( name='Sigma_ABbar2',  array=sigma_ABbar2,  format='D')
        c7 = fits.Column( name='Sigma_method2', array=sigma_method2, format='D')
        c8 = fits.Column( name='Sigma_ON2',     array=sigma_ON2,     format='{}D'.format(len(label_t)))
        c9 = fits.Column( name='RVfinal',       array=rvfinal,       format='D')
        c10 = fits.Column(name='STDfinal',      array=stdfinal,      format='D')

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

        print('sigma_method2 with type ={} is {}'.format(kind, sigma_method2))
        #print('RV/std for observations when IGRINS mounting was '+kind+': ', np.nanmean(rvfinal),np.nanstd(rvfinal))
        print('Observations when IGRINS is mounting {}: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( kind,
                                                                                                            np.nanmean(rvfinal),
                                                                                                            np.nanstd(rvfinal) ))
#-------------------------------------------------------------------------------
    xscale = np.arange(len(rvfinalCombined))+1

    f = plt.figure(figsize=(5,3))
    ax1 = plt.subplot(111)
    ax1.plot(xscale,rvfinalCombined, '.k', ms=5)
    ax1.errorbar(xscale,rvfinalCombined,yerr=stdfinalCombined,ls='none',lw=.5, ecolor='black')
#    ax1.text(1,np.nanmax(rvfinalCombined)+stdfinalCombined[np.nanargmax(rvfinalCombined)]+.03,
#             'Mean RV: '+str(round(np.nanmean(rvfinalCombined),5))+r'$\ \pm$ '+str(round(np.nanstd(rvfinalCombined),5))+' km/s')
    ax1.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=ax1.transAxes)

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        ax1.text(0.05, 0.1, 'Loose', transform=ax1.transAxes)
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            ax1.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            #ax1.text(xscale[len(nightsT)] -6.5,np.nanmin(rvfinalCombined)-2*stdfinalCombined[np.nanargmax(rvfinalCombined)],'Tight')
            #ax1.text(xscale[len(nightsT)] +1.5,np.nanmin(rvfinalCombined)-2*stdfinalCombined[np.nanargmax(rvfinalCombined)],'Loose')
            ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
            ax1.text(0.9,  0.1, 'Loose', transform=ax1.transAxes)
        else:
            ax1.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            #ax1.text(xscale[len(nightsL)] +1.5,np.nanmin(rvfinalCombined)-2*stdfinalCombined[np.nanargmax(rvfinalCombined)],'Tight')
            #ax1.text(xscale[len(nightsL)] -6.5,np.nanmin(rvfinalCombined)-2*stdfinalCombined[np.nanargmax(rvfinalCombined)],'Loose')
            ax1.text(0.05, 0.1, 'Tight', transform=ax1.transAxes)
            ax1.text(0.9,  0.1, 'Loose', transform=ax1.transAxes)
    ax1.set_ylim(np.nanmin(rvfinalCombined)-.08,np.nanmax(rvfinalCombined)+.08)
    ax1.set_ylabel('RV (km/s)')
    ax1.set_xlabel('Night (#)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/{}/FinalRVs.png'.format(inparam.outpath, name), format='png', bbox_inches='tight')


    c1 = fits.Column(name='NIGHT',    array=nightsCombined,         format='8A')
    c2 = fits.Column(name='MJD',      array=mjdsCombined-2400000.5, format='D')
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

    print('Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(rvfinalCombined),
                                                                            np.nanstd(rvfinalCombined)))
    print('vsini results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(vsinifinalCombined),
                                                                            np.nanstd(vsinifinalCombined)))
    print('\n')
    end_time = datetime.now()
    print('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    print('Output saved under {}/{}'.format(args.targname, name) )
    print('###############################################################')
    print('\n')
