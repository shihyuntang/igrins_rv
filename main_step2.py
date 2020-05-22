from Engine.importmodule import *

from Engine.IO_AB import setup_templates_syn, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips import basicclip_above
from Engine.contfit import A0cont
from Engine.classes import fitobjs,inparams
from Engine.macbro import macbro
from Engine.rebin_jv import rebin_jv
from Engine.rotint import rotint
from Engine.opt import optimizer
#-------------------------------------------------------------------------------

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

def ini_MPinst(i):
    nights   = inparam.nights
    targname = args.targname
    night    = nights[i]

    print('Working on order 1/1, night {} {}/{} ...'.format(night,
                                                            i+1,
                                                            len(inparam.nights)))
    # Only use the most precise order to hone in on adequate RV initial guesses
    order = 6
    # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
    if int(night) < 20180401 or int(night) > 20190531:
        ips = inparam.ips_tightmount[order]
        chunkweights = inparam.chunkweights_tightmount[order]
    else:
        ips = inparam.ips_loosemount[order]
        chunkweights = inparam.chunkweights_loosemount[order]

    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == list:
        initguesses = inparam.initguesses
    else:
        print('EXPECTED FILE OR LIST FOR INITGUESSES! QUITTING!')

    # Load telluric template from Telfit'd A0
    curdir = os.getcwd()
    A0loc = '{}/A0_Fits_{}/{}A0_treated.fits'.format(curdir, targname, night)
    try:
        hdulist = fits.open(A0loc)
    except IOError:
        print('No A0-fitted template for night '+night+', skipping...')
        return night,np.nan,np.nan

    tbdata = hdulist[order-1].data
    flag   = np.array(tbdata['ERRORFLAG'+str(order)])[0]

    if flag == 1: # Telfit hit unknown critical error
        return night,np.nan,np.nan

    watm = tbdata['WATM'+str(order)]
    satm = tbdata['SATM'+str(order)]
    a0contwave = watm.copy()
    continuum  = tbdata['BLAZE'+str(order)]

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]
    satm[(satm < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
    continuum = continuum[(a0contwave != 0)]

    # For Initialguesser, using combined beams, so take mean of tags' BVCs to get BVC we apply
    minibcs = []
    for tag in inparam.tagsA[night]:
        minibcs.append(inparam.bvcs[night+tag])
    for tag in inparam.tagsB[night]:
        minibcs.append(inparam.bvcs[night+tag])
    bcnight = np.mean(minibcs)

    # Define initial wavelength guesses (need only be approximate)
    initwave_dict = {
                        2: [-3.82127210e-10, -1.06269946e-05,  1.85070280e-01 ,  2.41718272e+04],
                        3: [7.170697761010893e-10,  -1.2196786264552286e-05,  0.18322308314595248, 23850.632494961632],
                        4: [-3.9083145742380996e-10, -8.907333871733518e-06,  0.17881172908636753 , 23537.219082651878],
                        5: [ -2.5348031649790425e-10, -9.180635514849177e-06, 0.17699954870465473, 23231.746232578 ],
                        6: [ 3.68743282e-10, -1.10364274e-05, 1.76781196e-01,  2.29338381e+04]
                    }
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,1.0,0.0,1.0,0.0,3.3,2.29315012e+04,1.75281163e-01,-9.92637874e-06,0,1.0,1e-4,-1e-7])

    # Arrays defining parameter variations during optimization steps
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.,1e7,1,1])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 5.00000e-5, 1e-7,0,0,0])
    dpar      = np.array([10.0, 1.0, 5.0, 3.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0,1e4,1e-3,1e-6]) # sy chnaged 5 -> 60

    vsinimini = []; rvsmini = [];

    # Iterate through initial RV guesses
    for initrvguess in initguesses:

        pars0[0] = initrvguess-bcnight # correct for barycentric velocity

        # Load target star spectrum. Processing similar to that of A0 spectra.
        x,wave,s,u = init_fitsread(inparam.inpath,'target','combined',night,order,None,None)

        nzones = 5
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

        # Split spectrum into 8 ~equal chunks in pixel space, then analyze all chunks but the ones on the ends

        Nsplit = 8
        wavegen = splitter(wave.copy(),Nsplit);
        fluxgen = splitter(s.copy(),   Nsplit);
        ugen    = splitter(u.copy(),   Nsplit);
        xgen    = splitter(x.copy(),   Nsplit);
        rvcollect = []; vsinicollect = [];

        for nn in range(Nsplit):
            wave_piece = next(wavegen);
            s_piece    = next(fluxgen);
            u_piece    = next(ugen);
            x_piece    = next(xgen);

            if nn == 0 or nn == Nsplit-1:
                continue

            mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

            satm_in = satm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]
            watm_in = watm[(watm > min(wave_piece)*1e4 - 11) & (watm < max(wave_piece)*1e4 + 11)]

            # Load initial IP guess, vsini settings
            par     = pars0.copy()
            par[5]  = ips[nn-1]
            par[4]  = inparam.initvsini
            dpar[4] = inparam.vsinivary

            # Cut target spec to be within A0 spec wave
            s_piece = s_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            u_piece = u_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            x_piece = x_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]
            wave_piece = wave_piece[(wave_piece*1e4 > min(watm_in)+5) & (wave_piece*1e4 < max(watm_in)-5)]

            ##Set initial wavelength guess
            try:
                initwavepars = initwave_dict[order]
                par[9] = initwavepars[0]; par[8] = initwavepars[1]; par[7] = initwavepars[2]; par[6] = initwavepars[3];
            except KeyError:
                f = np.polyfit(x_piece,wave_piece,3)
                par[9] = f[0]*1e4; par[8] = f[1]*1e4; par[7] = f[2]*1e4; par[6] = f[3]*1e4;

            s_piece /= np.median(s_piece)

            fitobj = fitobjs(s_piece, x_piece, u_piece, a0contwave,continuum,watm_in,satm_in,mflux_in,mwave_in)

            ######## Begin optimization  ########

            optimize = True
            par_in = par.copy()
            parfit_1 = optimizer(par,     dpar_cont, fitobj, optimize)
            parfit_2 = optimizer(parfit_1,dpar_wave, fitobj, optimize)

            parfit_1 = optimizer(parfit_2,dpar_cont, fitobj, optimize)
            parfit_2 = optimizer(parfit_1,dpar_wave, fitobj, optimize)

            parfit   = optimizer(parfit_2,  dpar,    fitobj, optimize)

            rv0 = parfit[0] - parfit[2]
            rvcollect.append(rv0 + bcnight + rv0*bcnight/(3e5**2))
            vsinicollect.append(parfit[4])

        rvcollect = np.array(rvcollect); vsinicollect = np.array(vsinicollect);
        rvsmini.append(np.nansum(rvcollect*chunkweights))
        vsinimini.append(np.nansum(vsinicollect*chunkweights))

    bestguess   = round(np.nanmean(rvsmini),5)

    return night,bestguess,np.mean(vsinimini)


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

    parser.add_argument('-dir',     dest="inpath",           action="store",
                        help="Enter path that stores target spectra, default will be under ./*targname",
                        type=str,   default='')
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
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()

    if args.inpath == '':
        args.inpath = './{}/'.format(args.targname)
    if args.inpath[-1] != '/':
        args.inpath+='/'

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)
    guesses   = args.guesses
    targname   = args.targname
    inpath     = args.inpath
    Nthreads   = args.Nthreads
    cdbs_loc = '~/cdbs/'
    #guesses = [-20,-10,0,10,20]

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
#        print('ERROR: INCORRECT "STYLE" INPUT PARAMETER SPECIFIED; EXPECTED EITHER "list" OR "file dir"')
#        print(breaker) # forces code to die

    #------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('''
Input Parameters:
    Tartget             = {}
    Initial vsini       = {} km/s
    vsini vary range    = {} km/s
    RV initial guess    = {} km/s
    '''.format(targname, initvsini, vsinivary, initguesses))
    print('---------------------------------------------------------------')
    print('RV Initial Guess for Each Night...')
    print('This will take a while..........')

    curdir = os.getcwd()

    ## Collect relevant file information from Predata files
    A0data   = Table.read('{}/Prepdata_A0_{}.txt'.format(curdir, targname), format='ascii')
    A0nights = np.array(A0data['night'],dtype='str')
    ams0     = np.array(A0data['airmass'])

    targdata = Table.read('{}/Prepdata_targ_{}.txt'.format(curdir, targname), format='ascii')
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

    # Create output directory
    try:
        filesndirs = os.listdir(os.getcwd()+'/'+targname)
    except OSError:
        os.mkdir(targname)
        filesndirs = os.listdir(os.getcwd()+'/'+targname)
    trk = 1; go = True;
    while go == True:
        iniguess_dir = 'Initguesser_results_'+str(trk)
        if iniguess_dir not in filesndirs:
            break
        trk += 1

    os.chdir(curdir+'/'+targname)
    os.mkdir(iniguess_dir)
    print('Writing output to file "'+targname+'/'+iniguess_dir+'"')
    filew = open(iniguess_dir+'/'+iniguess_dir,'w')
    filew.write('night bestguess vsini')
    filew.write('\n')
    os.chdir(curdir)

    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates_syn()

    inparam = inparams(inpath,iniguess_dir,initvsini,vsinivary,args.plotfigs,initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0)

    pool = mp.Pool(processes = Nthreads)
    outs = pool.map(ini_MPinst, np.arange(len(nightsFinal)))
    pool.close()
    pool.join()

    vsinis = []; finalrvs = [];
    for n in range(len(nightsFinal)):
        nightout = outs[n]
        filew.write('{} {} {}'.format(nightout[0], nightout[1], nightout[2]))
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
    print('Output saved under {}/{}'.format(targname, iniguess_dir) )
    print('---------------------------------------------------------------')
    print('You can now try to get a better RV initial guess with')
    print('--> python main_step2.py {} -i {:1.1f} -v {} -g [{:1.1f}] -c {} -plot'.format(targname,
                                                                                         np.nanmean(vsinis),
                                                                                         10,
                                                                                         np.nanmean(finalrvs),
                                                                                         Nthreads))
    print('OR go to next step with')
    print('--> python main_step3tar.py {} -i {:1.1f} -v {} -gX {} -c {} -plot'.format(targname,
                                                                                      np.nanmean(vsinis),
                                                                                      0,
                                                                                      trk,
                                                                                      Nthreads))
    print('###############################################################')
    print('\n')
