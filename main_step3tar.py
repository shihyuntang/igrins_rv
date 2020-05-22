from Engine.importmodule import *

from Engine.IO_AB import setup_templates_syn, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs, inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer


def splitter(master, N):

    # Yield successive length N pieces of an array, except for last one which
    # uses up any excess length leftover after array factored by N

    Npiece = int(len(master)/float(N))
    for ddd in range(N):
        if ddd != N-1:
            yield master[ddd*Npiece:(ddd+1)*Npiece]
        else:
            yield master[ddd*Npiece:-1]

# -------------------------------------------------------------------------------


def rv_main(order, i):
    # Main function to compute RVs for a given night and order
    # global fitobj, optimize;
    nights = inparam.nights
    targname = args.targname
    night = nights[i]

    print('Working on order {}/5, night {} {}/{} PID:{}...'.format(order-1,
                                                                   night,
                                                                   i+1,
                                                                   len(inparam.nights),
                                                                   mp.current_process().pid))

    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night) < 20180401 or int(night) > 20190531:
        ips = inparam.ips_tightmount[order]
    else:
        ips = inparam.ips_loosemount[order]

    # Collect relevant beam and filenum info
    tagsnight = []
    beamsnight = []
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    # Number of chunks spectrum will eventually be divided into (two will not be used)
    Nsplit = 8
    # Set up array for output
    rvsminibox = np.ones((len(tagsnight), Nsplit-2))
    vsiniminibox = np.ones((len(tagsnight), Nsplit-2))

    nightsout = []
    rvsminibox[:] = np.nan
    vsiniminibox[:] = np.nan

    for t in tagsnight:
        nightsout.append(night)

    # Load telluric template from Telfit'd A0
    curdir = os.getcwd()
    A0loc = curdir+'/A0_Fits_'+targname+'/'+night+'A0_treated.fits'
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('No A0-fitted template for night , skipping...'.format(night))
        return nightsout, rvsminibox, vsiniminibox

    tbdata = hdulist[order-1].data
    flag = np.array(tbdata['ERRORFLAG'+str(order)])[0]

    if flag == 1:  # Telfit hit unknown critical error
        return nightsout, rvsminibox, vsiniminibox

    if np.isnan(inparam.initguesses[night]):  # Telfit hit unknown critical error
        print('Initial guess for {} is NaN , skipping...'.format(night))
        return nightsout, rvsminibox, vsiniminibox

    watm = tbdata['WATM'+str(order)]
    satm = tbdata['SATM'+str(order)]
    a0contwave = tbdata['WAVE'+str(order)]  # sy change
    continuum = tbdata['BLAZE'+str(order)]

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]
    satm[(satm < 1e-4)] = 0.  # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
    continuum = continuum[(a0contwave != 0)]
    a0contwave = a0contwave[(a0contwave != 0)]  # sy change

    # Load relevant A0 spectra,
    # Define initial wavelength guesses (need only be approximate)
    initwave_dict = {
        2: [-3.82127210e-10, -1.06269946e-05, 1.85070280e-01,         2.41718272e+04],
        3: [7.17069776e-10, -1.21967862e-05, 0.18322308314595248, 23850.632494961632],
        4: [-3.90831457e-10, -8.90733387e-06, 0.17881172908636753, 23537.219082651878],
        5: [-2.53480316e-10, -9.18063551e-06, 0.17699954870465473, 23231.746232578],
        6: [3.68743282e-10, -1.10364274e-05, 1.76781196e-01,          2.29338381e+04]
    }
    # Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    # Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan, 1.0, 0.0, 1.0, 0.0, 3.3, 2.29315012e+04,
                      1.75281163e-01, -9.92637874e-06, 0, 1.0, 1e-4, -1e-7])

    # Arrays defining parameter variations during optimization steps
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3,  0.0,  0.0,        0.,   1e7, 1,    1])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 5.00000e-5, 1e-7, 0,   0,    0])
    #dpar = np.array([5.0, 1.0, 5.0, 3.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0,1e4,1e-3,1e-6])
    dpar = np.array([5.0, 1.0, 5.0, 3.0, 0.0, 0.5, 0.0,  0.0,  0.0,        0,    1e4, 1e-3, 1e-6])
    # Iterate over all A/B exposures
    for t in np.arange(len(tagsnight)):
        #print('  PID: {} opt, AB/mode: {}/{}'.format(mp.current_process().pid, t+1, len(tagsnight)))
        tag = tagsnight[t]
        beam = beamsnight[t]
        x, wave, s, u = init_fitsread(inparam.inpath+night+'/'+beam +
                                      '/', 'target', 'separate', night, order, tag, None)

        s2n = s/u
        if np.nanmedian(s2n) < 25:  # If S/N less than 25, throw out
            print('    Bad S/N ({:01.1f}) for night: {}, order: {}, node: {:04d}{}'.format(
                np.nanmedian(s2n), night, order, int(tag), beam))
            continue

        nzones = 5
        x = basicclip_above(x, s, nzones)
        wave = basicclip_above(wave, s, nzones)
        u = basicclip_above(u, s, nzones)
        s = basicclip_above(s, s, nzones)
        x = basicclip_above(x, s, nzones)
        wave = basicclip_above(wave, s, nzones)
        u = basicclip_above(u, s, nzones)
        s = basicclip_above(s, s, nzones)

        # Split spectrum into 8 ~equal chunks in pixel space, then analyze all chunks but the ones on the ends
        wavegen = splitter(wave.copy(), Nsplit)
        fluxgen = splitter(s.copy(),   Nsplit)
        ugen = splitter(u.copy(),   Nsplit)
        xgen = splitter(x.copy(),   Nsplit)

        for nn in range(Nsplit):
            wave_piece = next(wavegen)
            s_piece = next(fluxgen)
            u_piece = next(ugen)
            x_piece = next(xgen)

            if nn == 0 or nn == Nsplit-1:
                continue

            mwave_in, mflux_in = stellarmodel_setup(wave_piece, inparam.mwave0, inparam.mflux0)

            satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
            watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

            # Load initial IP guess, vsini settings
            par = pars0.copy()
            par[0] = inparam.initguesses[night] - inparam.bvcs[night+tag]
            par[5] = ips[nn-1]
            par[4] = inparam.initvsini
            dpar[4] = inparam.vsinivary

            # Cut target spec to be within A0 spec wave
            s_piece = s_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            u_piece = u_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            x_piece = x_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) &
                                    (wave_piece*1e4 < max(watm_in)-5)]

            # Set initial wavelength guess
            try:
                initwavepars = initwave_dict[order]
                par[9] = initwavepars[0]
                par[8] = initwavepars[1]
                par[7] = initwavepars[2]
                par[6] = initwavepars[3]
            except KeyError:
                f = np.polyfit(x_piece, wave_piece, 3)
                par[9] = f[0]*1e4
                par[8] = f[1]*1e4
                par[7] = f[2]*1e4
                par[6] = f[3]*1e4

            s_piece /= np.median(s_piece)

            fitobj = fitobjs(s_piece, x_piece, u_piece, a0contwave,
                             continuum, watm_in, satm_in, mflux_in, mwave_in)

            ######## Begin optimization  ########

            optimize = True
            par_in = par.copy()
            parfit_1 = optimizer(par,     dpar_cont, fitobj, optimize)
            parfit_2 = optimizer(parfit_1, dpar_wave, fitobj, optimize)
            par2 = parfit_2
            parfit_1 = optimizer(par2,    dpar_cont, fitobj, optimize)
            parfit_2 = optimizer(parfit_1, dpar_wave, fitobj, optimize)
            parfit = optimizer(parfit_2,  dpar,      fitobj, optimize)

            rv0 = parfit[0] - parfit[2]
            rvsminibox[t, nn-1] = rv0 + inparam.bvcs[night+tag] + \
                rv0*inparam.bvcs[night+tag]/(3e5**2)
            vsiniminibox[t, nn-1] = parfit[4]




    return nightsout, rvsminibox, vsiniminibox


def mp_run(Nthreads, order, nights):
    pool = mp.Pool(processes=Nthreads)
    func = partial(rv_main, order)
    outs = pool.map(func, np.arange(len(nights)))
    pool.close()
    pool.join()
    return outs

# -------------------------------------------------------------------------------


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='IGRINS Spectra Radial Velocity Pipeline',
        description='''
                                     This is a pipeline that helps you to extract radial velocity \n
                                     from IGRINS spectra. \n
                                     ''',
        epilog="Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)

    parser.add_argument('-dir',     dest="inpath",           action="store",
                        help="Enter path that stores target spectra, default will be under ./*targname",
                        type=str,   default='')
    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). Should use the value given by step2",
                        type=str,   default='')
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 0.0 km/s",
                        type=str, default='0.0')

    parser.add_argument('-gS',       dest="guesses_source",           action="store",
                        help="Source for initial guesses list for RV. Enter init OR rvre (init: Initguesser_results_X, rvre: RV_results_X)",
                        type=str, default='')
    parser.add_argument('-gX',       dest="guesses",           action="store",
                        help="Please give the number, X, under ./*targname/Initguesser_results_X OR ./*targname/RV_results_X, that you wish to use",
                        type=int, default='')

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)" % (
                            mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2))
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If sets, will generate plots")
    parser.add_argument('--version',
                        action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()

    if args.inpath == '':
        args.inpath = './{}/'.format(args.targname)
    if args.inpath[-1] != '/':
        args.inpath += '/'

    vsinivary = float(args.vsinivary)

    if args.initvsini != '':
        initvsini = float(args.initvsini)
    else:
        sys.exit('ERROR: EXPECTED FLOAT')

    if (args.guesses_source != 'init') & (args.guesses_source != 'rvre'):
        sys.exit('ERROR: ONLY "init" OR "rvre" IS ALLOW')
    # Collect init RV guesses
    if args.guesses != '':
        if args.guesses_source == 'init':
            guesses = './{}/Initguesser_results_{}/Initguesser_results_{}'.format(args.targname,
                                                                                  int(args.guesses),
                                                                                  int(args.guesses))
            guessdata = Table.read(guesses, format='ascii')
            initnights = np.array(guessdata['night'])
            initrvs = np.array(guessdata['bestguess'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre':
            guesses = './{}/RV_results_{}/RVresultsSummary.csv'.format(args.targname,
                                                                       int(args.guesses))
            guessdata = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs = np.array(guessdata['RVfinal'])
            initguesses = {}
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])
    else:
        sys.exit('ERROR: INCORRECT "STYLE" INPUT PARAMETER SPECIFIED; EXPECTED A INT NUMBER')

    cdbs_loc = '~/cdbs/'
    targname = args.targname
    inpath = args.inpath
    Nthreads = args.Nthreads

    # ------------
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
    '''.format(targname, initvsini, vsinivary, guesses))
    print('---------------------------------------------------------------')
    print('RV calculation for target star {}...'.format(targname))
    print('This will take a while..........')
    # ------------

    curdir = os.getcwd()
    # Collect relevant file information from Predata files
    A0data = Table.read(curdir+'/Prepdata_A0_'+targname+'.txt', format='ascii')
    A0nights = np.array(A0data['night'], dtype='str')
    ams0 = np.array(A0data['airmass'])

    targdata = Table.read(curdir+'/Prepdata_targ_'+targname+'.txt', format='ascii')
    Tnights = np.array(targdata['night'], dtype='str')
    tags0 = np.array(targdata['tag'],  dtype='int')
    beams0 = np.array(targdata['beam'], dtype='str')
    mjds0 = np.array(targdata['mjd'])
    bvcs0 = np.array(targdata['bvc'])
    ams = np.array(targdata['airmass'])

    # Attribute A and B exposures to right file numbers
    tagsA = {}
    tagsB = {}
    mjds = {}
    bvcs = {}
    night_orig = Tnights[0]
    tagsA0 = []
    tagsB0 = []
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
            tagsA0 = []
            tagsB0 = []
            if beams0[hrt] == 'A':
                tagsA0.append(tag1)
            else:
                tagsB0.append(tag1)
            night_orig = Tnights[hrt].copy()

    tagsA[Tnights[-1]] = tagsA0
    tagsB[Tnights[-1]] = tagsB0

    nightsFinal = np.array(list(sorted(set(Tnights))))

    # Create output directory
    try:
        filesndirs = os.listdir(os.getcwd()+'/'+targname)
    except OSError:
        os.mkdir(targname)
        filesndirs = os.listdir(os.getcwd()+'/'+targname)
    trk = 1
    go = True
    while go == True:
        name = 'RV_results_'+str(trk)
        if name not in filesndirs:
            break
        trk += 1

    os.chdir(os.getcwd()+'/'+targname)
    os.mkdir(name)
    print('\n')
    print('Writing output to folder "'+targname+'/'+name+'"')
    os.chdir(curdir)

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates_syn()

    # Takes about  seconds to do all 5 orders for a single night, but exact time will vary with number of separate exposures per night
    #print('Will analyze 5 orders of '+str(len(nightsFinal))+' nights, expected runtime: '+str(round(len(nightsFinal)*1000./(3600.*Nthreads),2))+' hours')
    print('\n')

    inparam = inparams(inpath, name, initvsini, vsinivary, args.plotfigs,
                       initguesses, bvcs, tagsA, tagsB, nightsFinal, mwave0, mflux0)

    orders = [2, 3, 4, 5, 6]

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T)
    nights = inparam.nights
    intnights = nights.astype(int)

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]

    rvmasterboxT = np.ones((len(nightsT), len(orders)))
    stdmasterboxT = np.ones((len(nightsT), len(orders)))
    rvmasterboxL = np.ones((len(nightsL), len(orders)))
    stdmasterboxL = np.ones((len(nightsL), len(orders)))
    vsinisT = np.ones((len(nightsT), len(orders)))
    vsinisL = np.ones((len(nightsL), len(orders)))

    if len(nightsL) > 0:
        nightscomblist = [nightsT, nightsL]
    else:
        nightscomblist = [nightsT]

    vsinis = np.ones(len(orders))
    vsinistds = np.ones(len(orders))

    for jerp in range(len(orders)):  # Iterate over orders
        order = orders[jerp]
        outs = mp_run(Nthreads, order, nights)

        # Collect outputs
        for i in range(len(nights)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvchunkbox = outsbox[1]
                vsinichunkbox = outsbox[2]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvchunkbox = np.vstack((rvchunkbox, outsbox[1]))
                vsinichunkbox = np.vstack((vsinichunkbox, outsbox[2]))

        nightsbox = np.array(nightsbox)
        vsinitags = []

        # Save results in fits file
        c1 = fits.Column(name='NIGHT', array=nightsbox,    format='8A')
        c2 = fits.Column(name='RV',   array=rvchunkbox,   format=str(
            len(rvchunkbox[0, :]))+'D', dim=(1, len(rvchunkbox[0, :])))
        c3 = fits.Column(name='VSINI', array=vsinichunkbox, format=str(
            len(rvchunkbox[0, :]))+'D', dim=(1, len(rvchunkbox[0, :])))
        cols = fits.ColDefs([c1, c2, c3])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if order == 2:  # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3, 3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu, hdu_1])
            hdul.writeto(targname+'/'+inparam.outpath+'/RVresultsRawBox.fits')
        else:
            hh = fits.open(targname+'/'+inparam.outpath+'/RVresultsRawBox.fits')
            hh.append(hdu_1)
            hh.writeto(targname+'/'+inparam.outpath+'/RVresultsRawBox.fits', overwrite=True)

        # For each set of nights (tight, loose)...
        T_L = 'T'
        for nights_use in nightscomblist:
            # Select relevant indices
            indbox = []
            for i in range(len(nightsbox)):
                if nightsbox[i] in nights_use:
                    indbox.append(i)
            indbox = np.array(indbox)

            # Note rvchunkbox indexed as [t,nn]
            # Use the deviation of RV over time within a chunk position to determine relative weights.
            # These are saved from GJ 281 anaylsis results. 
            if T_L == 'T':
                chunkweights = inparam.chunkweights_tightmount[order]
            else:
                chunkweights = inparam.chunkweights_loosemount[order]

            # Iterating over nights, use these weights to combine the RVs from different chunks of one spectrum into a single RV
            for i in range(len(nights_use)):

                indnight = np.where(nightsbox == nights_use[i])[0]

                rvtags = []
                for a in range(len(indnight)):
                    if np.nansum(rvchunkbox[indnight[a], :]*chunkweights) != 0:
                        rvtags.append(np.nansum(rvchunkbox[indnight[a], :]*chunkweights))
                        vsinitags.append(np.nansum(vsinichunkbox[indnight[a], :]*chunkweights))

                # Take the mean and std of the RVs determined from different A/B exposures within a night.
                if T_L == 'T':
                    vsinisT[i, jerp] = np.nanmean(vsinitags)

                    if len(rvtags) == 0 or len(rvtags) == 1:
                        rvmasterboxT[i, jerp] = np.nan
                        stdmasterboxT[i, jerp] = np.nan
                    else:
                        rvmasterboxT[i, jerp] = np.nanmean(rvtags)
                        stdmasterboxT[i, jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))
                else:
                    vsinisL[i, jerp] = np.nanmean(vsinitags)

                    if len(rvtags) == 0 or len(rvtags) == 1:
                        rvmasterboxL[i, jerp] = np.nan
                        stdmasterboxL[i, jerp] = np.nan
                    else:
                        rvmasterboxL[i, jerp] = np.nanmean(rvtags)
                        stdmasterboxL[i, jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags))

            T_L = 'L'

        # Take the straight mean and std of the vsinis collected, which were combined using the same weightings as the RVs
        #vsinis[jerp] = np.nanmean(vsinitags); vsinistds[jerp] = np.nanstd(vsinitags);

    nightsCombined = np.array([])
    mjdsCombined = np.array([])
    rvfinalCombined = np.array([])
    stdfinalCombined = np.array([])
    vsinifinalCombined = np.array([])

    if len(nightsL) > 0:
        rvboxcomblist = [rvmasterboxT, rvmasterboxL]
        stdboxcomblist = [stdmasterboxT, stdmasterboxL]
        vsinicomblist = [vsinisT, vsinisL]
    else:
        rvboxcomblist = [rvmasterboxT]
        stdboxcomblist = [stdmasterboxT]
        vsinicomblist = [vsinisT]

    # Iterating over tight and loose mounting nights...
    for boxind in range(len(rvboxcomblist)):

        rvmasterbox = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]
        vsinibox = vsinicomblist[boxind]

        # Load the uncertainty from method from GJ 281 analysis (this is from init RV = 20.02, not 20.15)
        # sigma_method2 = np.array([0.00298042, 0.00157407, 0.00109829, 0.00198759, 0.00112364])
        if boxind == 0:
            nights_use = nightsT.copy()
            kind = 'Tight'
            sigma_method2 = inparam.methodvariance_tight
        else:
            nights_use = nightsL.copy()
            kind = 'Loose'
            sigma_method2 = inparam.methodvariance_loose

        # Note rvmasterbox indexed as [nights,orders]
        Nnights = len(rvmasterbox[:, 0])  # sy

        # Calculate the uncertainty in each night/order RV as the sum of the uncertainty in method and the uncertainty in that night's As and Bs RVs
        sigma_ON2 = np.ones_like(rvmasterbox)
        for order in range(len(orders)):
            for night in range(Nnights):
                sigma_ON2[night, order] = sigma_method2[order] + stdmasterbox[night, order]**2

        rvfinal = np.ones(Nnights, dtype=np.float64)
        stdfinal = np.ones(Nnights, dtype=np.float64)
        vsinifinal = np.ones(Nnights, dtype=np.float64)
        mjds_out = np.ones(Nnights, dtype=np.float64)

        if boxind == 0:
            nights_use = nightsT.copy()
            kind = 'Tight'
        else:
            nights_use = nightsL.copy()
            kind = 'Loose'

        # Combine RVs between orders using weights calculated from uncertainties
        for n in range(Nnights):
            weights = (1./sigma_ON2[n, :])/(np.nansum(1./sigma_ON2[n, :]))  # normalized
            stdspre = (1./sigma_ON2[n, :])  # unnormalized weights
            rvfinal[n] = np.nansum(weights*rvmasterbox[n, :])
            stdfinal[n] = 1/np.sqrt(np.nansum(stdspre))

            mjds_out[n] = mjds[nights_use[n]]
            vsinifinal[n] = np.nansum(weights*vsinibox[n, :])

            if np.nansum(weights) == 0:
                rvfinal[n] = np.nan
                stdfinal[n] = np.nan
                vsinifinal[n] = np.nan

        # Print, plot, and save results
        #print('RV/std for observations when IGRINS mounting was '+kind+': ', np.nanmean(rvfinal),np.nanstd(rvfinal))
        print('Observations when IGRINS is mounting {}: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(kind,
                                                                                                           np.nanmean(
                                                                                                               rvfinal),
                                                                                                           np.nanstd(rvfinal)))

        f = plt.figure(figsize=(5, 3))
        ax1 = plt.subplot(111)
        ax1.plot(np.arange(len(rvfinal))+1, rvfinal, '.k', ms=5)
        ax1.errorbar(np.arange(len(rvfinal))+1, rvfinal,
                     yerr=stdfinal, ls='none', lw=0.5, ecolor='black')
#        ax1.text(1,np.nanmax(rvfinal)+stdfinal[np.nanargmax(rvfinal)]+.03,
#                 'Mean RV: '+str(round(np.nanmean(rvfinal),5))+r'$\ \pm$ '+str(round(np.nanstd(rvfinal),5))+' km/s')
        ax1.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinal), np.nanstd(rvfinal)),
                             transform=ax1.transAxes)
        ax1.set_ylim(np.nanmin(rvfinal)-.08, np.nanmax(rvfinal)+.08)
        ax1.set_ylabel('RV (km/s)')
        ax1.set_xlabel('Night (#)')
        ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in', width=.6)
        f.savefig(targname+'/'+inparam.outpath+'/'+'FinalRVs_' +
                  kind+'_.png', format='png', bbox_inches='tight')
        plt.clf()
        plt.close()

        c1 = fits.Column(name='NIGHT',         array=nights_use,    format='8A')
        c2 = fits.Column(name='MJD',           array=mjds_out-2400000.5,      format='D')
        c3 = fits.Column(name='RVBOX',         array=rvmasterbox,   format='5D')
        c4 = fits.Column(name='STDBOX',        array=stdmasterbox,  format='5D')
#        c5 = fits.Column(name='Sigma_O2',      array=sigma_O2,      format='D')
#        c6 = fits.Column(name='Sigma_ABbar2',  array=sigma_ABbar2,  format='D')
        c7 = fits.Column(name='Sigma_method2', array=sigma_method2, format='D')
        c8 = fits.Column(name='Sigma_ON2',     array=sigma_ON2,     format='5D')
        c9 = fits.Column(name='RVfinal',       array=rvfinal,       format='D')
        c10 = fits.Column(name='STDfinal',     array=stdfinal,      format='D')

        cols = fits.ColDefs([c1, c2, c3, c4, c7, c8, c9, c10])
        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3, 3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul = fits.HDUList([primary_hdu, hdu_1])
        hdul.writeto(targname+'/'+inparam.outpath+'/RVresultsSummary_'+kind+'.fits')

        #tempin = Table.read(targname+'/'+inparam.outpath+'/RVresultsSummary_'+kind+'.fits', format='fits')
        #tempin.write(targname+'/'+inparam.outpath+'/RVresultsSummary_'+kind+'.csv', format='csv', overwrite=True)

        nightsCombined = np.concatenate((nightsCombined, nights_use))
        mjdsCombined = np.concatenate((mjdsCombined, mjds_out))
        rvfinalCombined = np.concatenate((rvfinalCombined, rvfinal))
        stdfinalCombined = np.concatenate((stdfinalCombined, stdfinal))
        vsinifinalCombined = np.concatenate((vsinifinalCombined, vsinifinal))

    # Print, plot, and save combined results. These have been analyzed separately, but slapped together at the very end.

    xscale = np.arange(len(rvfinalCombined))+1

    f = plt.figure(figsize=(5, 3))
    ax1 = plt.subplot(111)
    ax1.plot(xscale, rvfinalCombined, '.k', ms=5)
    ax1.errorbar(xscale, rvfinalCombined, yerr=stdfinalCombined, ls='none', lw=0.5, ecolor='black')
#    ax1.text(1,np.nanmax(rvfinalCombined)+stdfinalCombined[np.nanargmax(rvfinalCombined)]+.03,
#             'Mean RV: '+str(round(np.nanmean(rvfinalCombined),5))+r'$\ \pm$ '+str(round(np.nanstd(rvfinalCombined),5))+' km/s')
    ax1.text(0.05, 0.93, r'RV mean= {:1.5f} $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=ax1.transAxes)
    if nightsT[-1] < nightsL[0]:  # if tight epoch precedes loose epoch
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
    ax1.set_ylim(np.nanmin(rvfinalCombined)-.08, np.nanmax(rvfinalCombined)+.08)
    ax1.set_ylabel('RV (km/s)')
    ax1.set_xlabel('Night (#)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig(targname+'/'+inparam.outpath+'/'+'FinalRVs.png', format='png', bbox_inches='tight')
    plt.clf()
    plt.close()

    c1 = fits.Column(name='NIGHT',            array=nightsCombined,   format='8A')
    c2 = fits.Column(name='MJD',               array=mjdsCombined-2400000.5,     format='D')
    c3 = fits.Column(name='RVfinal',          array=rvfinalCombined,  format='D')
    c4 = fits.Column(name='STDfinal',         array=stdfinalCombined, format='D')
    c5 = fits.Column(name='VSINI',            array=vsinifinalCombined,           format='D')

    cols = fits.ColDefs([c1, c2, c3, c4, c5])
    hdu_1 = fits.BinTableHDU.from_columns(cols)
    bleh = np.ones((3, 3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu, hdu_1])
    hdul.writeto(targname+'/'+inparam.outpath+'/RVresultsSummary.fits')

    tempin = Table.read(targname+'/'+inparam.outpath+'/RVresultsSummary.fits', format='fits')
    tempin.write(targname+'/'+inparam.outpath+'/RVresultsSummary.csv', format='csv', overwrite=True)

    print('Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(rvfinalCombined),
                                                                            np.nanstd(rvfinalCombined)))
    print('vsini results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(vsinifinalCombined),
                                                                            np.nanstd(vsinifinalCombined)))
    print('\n')
    end_time = datetime.now()
    print('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    print('Output saved under {}/{}'.format(targname, name))
    print('###############################################################')
    print('\n')
