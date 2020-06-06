import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *

from Engine.IO_AB import setup_templates, setup_templates_syn, setup_templates_sun, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs,inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer, fmod

#-------------------------------------------------------------------------------
def rv_main(i, order0, order):
    # Main function to compute RVs for a given night and order
    # global fitobj, optimize;
    nights   = inparam.nights
    targname = args.targname
    night    = i[0]

    print('Working on {} band, order {}/{}, night {} PID:{}...'.format(args.band,
                                                                          order,
                                                                          len(order0),
                                                                          night,
                                                                          mp.current_process().pid) )

    # --------------------------------------------------------------
    # --------------------------------------------------------------
    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
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

    # Number of chunks spectrum will eventually be divided into (two will not be used)
    Nsplit = 8

    wminibox      = np.ones((270, Nsplit))
    sminibox      = np.ones((270, Nsplit))
    flminibox_tel = np.ones((270, Nsplit))
    flminibox_ste = np.ones((270, Nsplit))
    contiminibox  = np.ones((270, Nsplit))
    residualbox   = np.ones((270, Nsplit))
    flminibox_mod = np.ones((270, Nsplit))

    wminibox[:]     = np.nan
    sminibox[:]     = np.nan
    flminibox_tel[:]= np.nan
    flminibox_ste[:]= np.nan
    contiminibox[:] = np.nan
    residualbox[:]  = np.nan
    flminibox_mod[:]  = np.nan

    # Load telluric template from Telfit'd A0
    A0loc = '{}/A0_Fits/{}A0_treated_{}.fits'.format(args.targname, night[:8], args.band)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('No A0-fitted template for night {}, skipping...'.format(night))
        return wminibox,sminibox,flminibox_tel,flminibox_ste,contiminibox,residualbox

    # tbdata = hdulist[order-1].data
    # flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]
    num_orders = len( np.unique(order0) )

    # order in A0_treated.fits is no longer sequential...
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]


    if flag == 1: # Telfit hit unknown critical error
        return  wminibox,sminibox,flminibox_tel,flminibox_ste,contiminibox,residualbox

    if np.isnan(inparam.initguesses): # Telfit hit unknown critical error
        print('Initial guess for {} is NaN , skipping...'.format(night))
        return wminibox,sminibox,flminibox_tel,flminibox_ste,contiminibox,residualbox


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
                      0.5])                                                 #15: Differential Rotation Coefficient

    # Iterate over all A/B exposures
    for t in [0]:
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

        s_piece    = s
        u_piece    = u
        wave_piece = wave
        x_piece    = x

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
        dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0, 0  ])
        dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0, 0  ])
        dpar      = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0.2])
        dpar_st   = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0.2])
        dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0,                 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0  ])
                     #'st'   : np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.0, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        s_piece /= np.median(s_piece)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in)
#-------------------------------------------------------------------------------
        ######## Begin optimization  ########

        optimize = True
        par_in = par.copy()
        hardbounds = [par_in[4] -dpar[4],   par_in[4]+dpar[4],
                      par_in[5] -dpar[5],   par_in[5]+dpar[5],
                      par_in[15]-dpar[15], par_in[15]+dpar[15]]
        if hardbounds[0] < 0:
            hardbounds[0] = 0
        if hardbounds[3] < 0:
            hardbounds[3] = 1
        if hardbounds[4] < 0.1:
            hardbounds[4] = 0.1
        if hardbounds[5] > 0.9:
            hardbounds[5] = 0.9

        parfit_1 = optimizer(par_in,   dpar_cont, hardbounds,fitobj,optimize)
        parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds,fitobj,optimize)
        parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds,fitobj,optimize)
        parfit_4 = optimizer(parfit_3, dpar_wave, hardbounds,fitobj,optimize)
        parfit = optimizer(parfit_4,   dpar,      hardbounds,fitobj,optimize)   # RV fitting
        #
        # # if stellar template power is very low, throw out result
        # if parfit[1] < 0.1:
        #     print('parfit[1] < 0.1')
        #     continue
        #
        # # if stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        # if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
        #     print(' parfit[1] == par_in[1] or parfit[3] == par_in[3]')
        #     continue
        #
        # # # if model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        # smod,chisq = fmod(parfit,fitobj)
        # if len(smod[(smod < 0)]) > 0:
        #     print('len(smod[(smod < 0)]) > 0')
        #     continue

        if args.plotfigs == True:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter(parfitS, fitobj,'{}_{}_{}_parfitS'.format(label,night,tag), trk, 0)
            outplotter(parfitT, fitobj,'{}_{}_{}_parfitT'.format(label,night,tag), trk, 0)
            outplotter(parfit, fitobj,'{}_{}_{}_parfit'.format(label,night,tag), trk, 0)

        if args.debug == True:
            outplotter(parfit_1,fitobj,'{}_{}_{}_parfit_1'.format(label,night,tag), trk, 1)
            outplotter(parfit_2,fitobj,'{}_{}_{}_parfit_2'.format(label,night,tag), trk, 1)
            outplotter(parfit_3,fitobj,'{}_{}_{}_parfit_3'.format(label,night,tag), trk, 1)
            outplotter(parfit_4,fitobj,'{}_{}_{}_parfit_4'.format(label,night,tag), trk, 1)
            outplotter(parfit_5,fitobj,'{}_{}_{}_parfit_5'.format(label,night,tag), trk, 1)
            outplotter(parfit  ,fitobj,'{}_{}_{}_parfit'.format(label,night,tag), trk, 1)

        # Compute model and divide for residual
        fullmodel,chisq = fmod(parfit,fitobj)
        residual = fitobj.s/fullmodel


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

    w_min = np.nanmin(w)
    dw    = (np.nanmax(w) - np.nanmin(w)) / 8
    for nn in range(8):
        wrange = [ (w > (w_min + nn*dw) ) & ((w < (w_min + (nn+1)*dw))) ][0]
        leng_w = sum(wrange)
        wminibox[:leng_w, nn]         = w[wrange]
        sminibox[:leng_w, nn]         = dataflat[wrange]
        flminibox_mod[:leng_w, nn]    = modelflat[wrange]
        flminibox_tel[:leng_w, nn]    = tellflat[wrange]
        flminibox_ste[:leng_w, nn]    = stellflat[wrange]
        contiminibox[:leng_w, nn]     = contmodel[wrange]
        residualbox[:leng_w, nn]      = residual[wrange]

    return wminibox,sminibox,flminibox_mod,flminibox_tel,flminibox_ste,contiminibox,residualbox


def mp_run(Nthreads, nights, order0):
    pool = mp.Pool(processes = Nthreads)
    func = partial(rv_main, nights, order0)
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
                        help="Minium request of # of AB sets. Default = for STD is 1 and TAR is 3 ",
                        type=str,   default='1')

    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). Should use the value given by step2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 0.0 km/s",
                        type=str, default='0.0' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help=". Should use the single value given by step2 (float, km/s)",
                        type=str,   default='' )

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
    inpath    = '../Input_Data/{}/'.format(args.targname)
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
    A0data   = Table.read('../Temp/Prepdata/Prepdata_A0_{}_tool.txt'.format(args.targname), format='ascii')
    A0nights = np.array(A0data['night'],dtype='str')
    ams0     = np.array(A0data['airmass'])

    targdata = Table.read('../Temp/Prepdata/Prepdata_targ_{}_tool.txt'.format(args.targname), format='ascii')
    Tnights = np.array(targdata['night'],dtype='str')
    tags0   = np.array(targdata['tag'], dtype='int')
    beams0  = np.array(targdata['beam'],dtype='str')
    mjds0   = np.array(targdata['mjd'])
    bvcs0   = np.array(targdata['bvc'])
    ams     = np.array(targdata['airmass'])

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

    if len(nightsFinal)==1:
        print(nightsFinal, len(nightsFinal))
    else:
        sys.exit('only take one night!, we give {}'.format(nightsFinal))

    # if args.nights_use != '':
    #     nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
    #     for nnn in nightstemp:
    #         if nnn not in nightsFinal:
    #             sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
    #     nightsFinal = nightstemp
    #     print('Only processing nights: {}'.format(nightsFinal))
#-------------------------------------------------------------------------------

    # Create output directory
    if not os.path.isdir('./{}'.format(args.targname)):
        os.mkdir('./{}'.format(args.targname) )

    filesndirs = os.listdir('./{}'.format(args.targname))
    trk = 1
    go = True
    while go == True:
        name = 'RV_results_'+str(trk)
        if name not in filesndirs:
            break
        trk += 1

    if not os.path.isdir('./{}/{}'.format(args.targname, name)):
        os.mkdir('./{}/{}'.format(args.targname, name) )

    print('Writing output to folder "'+args.targname+'/'+name+'"')

    # Retrieve stellar and telluric templates
    # if args.band=='K':
    #     watm,satm, mwave0, mflux0 = setup_templates_syn()
    # elif args.band=='H':
    #     watm,satm, mwave0, mflux0 = setup_templates_sun()

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

    print('\n')
    outpath = './{}'.format(args.targname)

    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,None)

#    orders = [2,3,4,5,6]
# ---------------------------------------
    if args.band == 'K':
        orders = np.arange(2,17)
#        orders = np.array([6])
    elif args.band == 'H':
        # orders = np.arange(2,23)
        orders = np.array([2, 3, 4, 5, 6, 10, 11, 13, 14, 16, 17, 20, 21, 22])
#    order0 = np.array([16])
# ---------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T)
    nights    = inparam.nights
    intnights = nights.astype(int)


    outs = mp_run(args.Nthreads, nights, orders)
# Collect outputs
for i in range(len(orders)):
    outsbox = outs[i]

    wbox      = outsbox[0]
    stbox     = outsbox[1]
    modbox    = outsbox[2]
    telbox    = outsbox[3]
    stebox    = outsbox[4]
    conti_fl  = outsbox[5]
    residual  = outsbox[6]

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
    c7 = fits.Column(name='residual',      array=residual,    format=str(
        len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))


    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
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
