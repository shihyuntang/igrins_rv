import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *

from Engine.IO_AB import setup_templates, setup_templates_syn, setup_templates_sun, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs,inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer, fmodel_separate

def splitter(master,N):

    # Yield successive length N pieces of an array, except for last one which
    # uses up any excess length leftover after array factored by N

    Npiece = int(len(master)/float(N))
    for ddd in range(N):
        if ddd != N-1:
            yield master[ddd*Npiece:(ddd+1)*Npiece]
        else:
            yield master[ddd*Npiece:-1]

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

    # Number of chunks spectrum will eventually be divided into (two will not be used)
    Nsplit = 8
    # Set up array for output
    rvsminibox   = np.ones((len(tagsnight),Nsplit));
    vsiniminibox = np.ones((len(tagsnight),Nsplit));

    wminibox      = np.ones((270, Nsplit))
    flminibox_tel = np.ones((270, Nsplit))
    flminibox_ste = np.ones((270, Nsplit))
    contiminibox  = np.ones((270, Nsplit))
    stalflatbox   = np.ones((270, Nsplit))
    ubox          = np.ones((270, Nsplit))
    orgfluxbox    = np.ones((270, Nsplit))

    nightsout = [];
    rvsminibox[:]   = np.nan;
    vsiniminibox[:] = np.nan;

    wminibox[:]     = np.nan
    flminibox_tel[:]= np.nan
    flminibox_ste[:]= np.nan
    contiminibox[:] = np.nan
    stalflatbox[:]  = np.nan
    ubox[:]         = np.nan
    orgfluxbox[:]   = np.nan

    for t in tagsnight:
        nightsout.append(night)

    # Load telluric template from Telfit'd A0
    A0loc = './{}/A0_Fits/{}A0_treated_{}.fits'.format(args.targname, night, args.band)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('No A0-fitted template for night , skipping...'.format(night))
        return wminibox,stalflatbox,flminibox_tel,flminibox_ste,ubox,orgfluxbox,contiminibox

    # tbdata = hdulist[order-1].data
    # flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]
    num_orders = len( np.unique(order0) )

    # order in A0_treated.fits is no longer sequential...
    fits_layer = [ i for i in np.arange(num_orders)+1 if int(hdulist[i].columns[0].name[9:]) == order ][0]

    tbdata = hdulist[ fits_layer ].data
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]


    if flag == 1: # Telfit hit unknown critical error
        return  wminibox,stalflatbox,flminibox_tel,flminibox_ste,ubox,orgfluxbox,contiminibox

    if np.isnan(inparam.initguesses): # Telfit hit unknown critical error
        print('Initial guess for {} is NaN , skipping...'.format(night))
        return wminibox,stalflatbox,flminibox_tel,flminibox_ste,ubox,orgfluxbox,contiminibox


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

#    for t in np.arange(len(tagsnight)):
    for t in [0]:
        #print('  PID: {} opt, AB/mode: {}/{}'.format(mp.current_process().pid, t+1, len(tagsnight)))
        tag  = tagsnight[t]
        beam = beamsnight[t]
    # x (list of wavelength used position)
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

        x,wave,s,u = init_fitsread(inparam.inpath+night+'/'+beam+'/',
                                   'target',
                                   'separate',
                                   night,
                                   order,
                                   tag,
                                   args.band,
                                   bound_cut)

        s2n = s/u
        if np.nanmedian(s2n) < 50: # If S/N less than 25, throw out
            print('    Bad S/N ({:01.1f}) for night: {}, order: {}, node: {:04d}{}'.format(np.nanmedian(s2n), night, order, int(tag), beam))
            continue

        nzones = 5
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

        # Split spectrum into 8 ~equal chunks in pixel space, then analyze all chunks but the ones on the ends
        wavegen = splitter(wave.copy(),Nsplit);
        fluxgen = splitter(s.copy(),   Nsplit);
        ugen    = splitter(u.copy(),   Nsplit);
        xgen    = splitter(x.copy(),   Nsplit);
            # Arrays defining parameter variations during optimization steps

        dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0])
        dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0,               0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0])
        dpar_st   = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0])
        dpar      = np.array([5.0, 1.0, 5.0, 3.0, inparam.vsinivary, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0])

        for nn in range(Nsplit):
            wave_piece = next(wavegen);
            s_piece    = next(fluxgen);
            u_piece    = next(ugen);
            x_piece    = next(xgen);

            mwave_in, mflux_in = stellarmodel_setup(wave_piece, inparam.mwave0, inparam.mflux0)

            satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
            watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

            # Load initial IP guess, vsini settings
            par = pars0.copy()
            par[0]  = inparam.initguesses - inparam.bvcs[night+tag]


            f = np.polyfit(x_piece,wave_piece,3)
            par[9] = f[0]*1e4;
            par[8] = f[1]*1e4;
            par[7] = f[2]*1e4;
            par[6] = f[3]*1e4;

            # Cut target spec to be within A0 spec wave
            s_piece = s_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            u_piece = u_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            x_piece = x_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

            continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
            s_piece /= np.median(s_piece)
            fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in,watm_in,satm_in,mflux_in,mwave_in)
#-------------------------------------------------------------------------------
            ######## Begin optimization  ########

            optimize = True
            par_in = par.copy()
            hardbounds = [par_in[4]-dpar[4],
                          par_in[4]+dpar[4],
                          par_in[5]-dpar[5],
                          par_in[5]+dpar[5]]
            if hardbounds[0] < 0:
                hardbounds[0] = 0
            if hardbounds[3] < 0:
                hardbounds[3] = 1

            parfit_1 = optimizer(par_in,   dpar_cont, hardbounds,fitobj,optimize)
            parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds,fitobj,optimize)
            parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds,fitobj,optimize)
            parfit_4 = optimizer(parfit_3, dpar_wave, hardbounds,fitobj,optimize)
            parfit = optimizer(parfit_4,   dpar,      hardbounds,fitobj,optimize)   # RV fitting

#            rv0 = parfit[0] - parfit[2]                         # atomosphere velocity correct

            parfit_tel = parfit.copy() # modified 0503
            parfit_tel[1] = 0
            w,smod_tel,cont,cont1 = fmodel_separate(parfit_tel)

            parfit_ste = parfit.copy() # modified 0503
            parfit_ste[3] = 0
            w,smod_ste,cont,cont1  = fmodel_separate(parfit_ste)

            s2n   = s_piece/u_piece
            sflat = s_piece/cont
            sflat *= np.median(cont)
            u_piece = sflat/s2n

            wminibox[:len(w), nn]                = w
            stalflatbox[:len(sflat), nn]         = sflat
            flminibox_tel[:len(smod_tel), nn]    = smod_tel
            flminibox_ste[:len(smod_ste), nn]    = smod_ste
            ubox[:len(u_piece), nn]              = u_piece
            orgfluxbox[:len(s_piece), nn]        = s_piece
            contiminibox[:len(cont), nn]         = cont1


    return wminibox,stalflatbox,flminibox_tel,flminibox_ste,ubox,orgfluxbox,contiminibox

def mp_run(Nthreads, nights, order0):
    pool = mp.Pool(processes = Nthreads)
    func = partial(rv_main, nights, order0)
    outs = pool.map(func, order0)
    pool.close()
    pool.join()
    return outs

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
    parser.add_argument("-HorK",    dest="band",            action="store",
                        help="Which band to process? H or K?",
                        type=str,   default='K')
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
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()
    cdbs_loc = '~/cdbs/'
    inpath     = '../Input_Data/{}/'.format(args.targname)

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

    #------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('###############################################################')
    print('''
Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    = {} km/s
    RV initial guess    using values in {}
    '''.format(args.targname, initvsini, vsinivary, initguesses))
    print('---------------------------------------------------------------')
    print('RV calculation for target star {}...'.format(args.targname))
    print('This will take a while..........')
    #------------

    ## Collect relevant file information from Predata files
    A0data = Table.read('../Temp/Prepdata/Prepdata_A0_{}_tool.txt'.format(args.targname), format='ascii')
    A0nights = np.array(A0data['night'],dtype='str')
    ams0     = np.array(A0data['airmass'])

    targdata =  Table.read('../Temp/Prepdata/Prepdata_targ_{}_tool.txt'.format(args.targname), format='ascii')
    Tnights  = np.array(targdata['night'],dtype='str')
    tags0    = np.array(targdata['tag'],  dtype='int')
    beams0   = np.array(targdata['beam'], dtype='str')
    mjds0    = np.array(targdata['mjd'])
    bvcs0    = np.array(targdata['bvc'])
    ams      = np.array(targdata['airmass'])

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


    if len(nightsFinal)==1:
        print(nightsFinal, len(nightsFinal))
    else:
        sys.exit('only take one night!, we give {}'.format(nightsFinal))
    # Create output directory

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
    if args.band=='K':
        watm,satm, mwave0, mflux0 = setup_templates_syn()
    elif args.band=='H':
        watm,satm, mwave0, mflux0 = setup_templates_sun()

    # Takes about  seconds to do all 5 orders for a single night, but exact time will vary with number of separate exposures per night
    #print('Will analyze 5 orders of '+str(len(nightsFinal))+' nights, expected runtime: '+str(round(len(nightsFinal)*1000./(3600.*Nthreads),2))+' hours')
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

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]

    rvmasterboxT  = np.ones((len(nightsT),len(orders)))
    stdmasterboxT = np.ones((len(nightsT),len(orders)))
    rvmasterboxL  = np.ones((len(nightsL),len(orders)))
    stdmasterboxL = np.ones((len(nightsL),len(orders)))

    if len(nightsL) > 0:
        nightscomblist = [nightsT,nightsL]
    else:
        nightscomblist = [nightsT]

    vsinis = np.ones(len(orders)); vsinistds = np.ones(len(orders));

#--------------------------
    outs = mp_run(args.Nthreads, nights, orders)
    # Collect outputs
    for i in range(len(orders)):
        outsbox = outs[i]

        wbox      = outsbox[0]
        stbox     = outsbox[1]
        telbox    = outsbox[2]
        stebox    = outsbox[3]
        uubox     = outsbox[4]
        orgbox    = outsbox[5]
        conti_fl  = outsbox[6]

        # Save results in fits file
        c1 = fits.Column(name='wavelength',    array=wbox,         format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c2 = fits.Column(name='sflat',         array=stbox,         format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c3 = fits.Column(name='tel_fl',        array=telbox,    format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c4 = fits.Column(name='ste_fl',        array=stebox,    format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c5 = fits.Column(name='u_piece',        array=uubox,    format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c6 = fits.Column(name='s_piece',        array=orgbox,    format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))
        c7 = fits.Column(name='conti_fl',      array=conti_fl,     format=str(
            len(wbox[0, :]))+'D', dim=(1, len(wbox[0, :])))

        cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if orders[i] == orders[0]:  # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3, 3))
            primary_hdu1 = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu1, hdu_1])
            hdul.writeto(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(targname, inparam.nights[0], args.band))
        else:
            hh = fits.open(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(targname, inparam.nights[0], args.band))
            hh.append(hdu_1)
            hh.writeto(tinparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}.fits'.format(targname, inparam.nights[0], args.band), overwrite=True)

    end_time = datetime.now()
    print('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    print('Output saved under {}/{}'.format(args.targname, name) )
    print('###############################################################')
    print('\n')
