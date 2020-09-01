from Engine.importmodule import *

from Engine.IO_AB     import setup_templates_tel, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips     import basicclip_above
from Engine.contfit   import A0cont
from Engine.classes   import fitobjs,inparamsA0,orderdict_cla
from Engine.macbro    import macbro
from Engine.rebin_jv  import rebin_jv
from Engine.rotint    import rotint
from Engine.Telfitter import telfitter
from Engine.opt       import optimizer, fmod
from Engine.outplotter import outplotter_tel
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def MPinst(args, inparam, jerp, orders, i):
    # Main function for A0 fitting that will be threaded over by multiprocessing

    order = orders[jerp]            # current looped order
    night = str(inparam.nights[i])  # multiprocess assigned night
    firstorder = orders[0]          # First order that will be analyzed, related to file writing

    print('Working on order {:02d}/{:02d} ({}), night {}/{} ({}) PID:{}...'.format(int(jerp+1),
                                                                                 len(orders),
                                                                                 order,
                                                                                 i+1,
                                                                                 len(inparam.nights),
                                                                                 night,
                                                                                 mp.current_process().pid) )

    #-------------------------------------------------------------------------------

    # Retrieve pixel bounds for where within each other significant telluric absorption is present.
    # If these bounds were not applied, analyzing some orders would give garbage fits.
    if args.band=='K':
        if int(order) in [14]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]

    elif args.band=='H':
        if int(order) in [13, 14, 16, 20]:
            bound_cut = inparam.bound_cut_dic[args.band][order]
        else:
            bound_cut = [150, 150]
    # bound_cut = [150, 150]


    ### Load relevant A0 spectrum
    x, a0wavelist, a0fluxlist, u = init_fitsread(inparam.inpath,
                                                 'A0',
                                                 'separate',
                                                 night,
                                                 order,
                                                 f'{int(inparam.tags[night]):04d}',
                                                 args.band,
                                                 bound_cut)
    #-------------------------------------------------------------------------------

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    nzones = 12
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(x,a0fluxlist,nzones);
    a0u        = basicclip_above(u,a0fluxlist,nzones);     a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(a0x,a0fluxlist,nzones);
    a0u        = basicclip_above(a0u,a0fluxlist,nzones);   a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);

    # Normalize
    a0fluxlist /= np.median(a0fluxlist)

    # Compute rough blaze function estimate. Better fit will be provided by Telfit later.
    continuum    = A0cont(a0wavelist,a0fluxlist,night,order)
    a0masterwave = a0wavelist.copy()
    a0masterwave *= 1e4

    # Trim stellar template to relevant wavelength range
    mwave_in, mflux_in = stellarmodel_setup(a0wavelist, inparam.mwave0, inparam.mflux0)

    # Trim telluric template to relevant wavelength range
    satm_in = inparam.satm[(inparam.watm > min(a0wavelist)*1e4 - 11) & (inparam.watm < max(a0wavelist)*1e4 + 11)]
    watm_in = inparam.watm[(inparam.watm > min(a0wavelist)*1e4 - 11) & (inparam.watm < max(a0wavelist)*1e4 + 11)]

    # Get initial guess for cubic wavelength solution from reduction pipeline
    f = np.polyfit(a0x, a0wavelist, 3)
    par9in = f[0]*1e4;
    par8in = f[1]*1e4;
    par7in = f[2]*1e4;
    par6in = f[3]*1e4;

    # Determine whether IGRINS mounting was loose or night for the night in question
    if (int(night) < 20180401) or (int(night) > 20190531):
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]

    ### Initialize parameter array for optimization as well as half-range values for each parameter during
    ### the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and
    ### in between optimization steps.
    parA0 = np.array([0.0,           # 0: The shift of the stellar template (km/s)
                      0.0,           # 1: The scale factor for the stellar template
                      0.0,           # 2: The shift of the telluric  template (km/s)
                      1.0,           # 3: The scale factor for the telluric template
                      0.0,           # 4: vsini (km/s)
                      IPpars[2],     # 5: The instrumental resolution (FWHM) in pixels
                      par6in,        # 6: Wavelength 0-pt
                      par7in,        # 7: Wavelength linear component
                      par8in,        # 8: Wavelength quadratic component
                      par9in,        # 9: Wavelength cubic component
                      1.0,           #10: Continuum zero point
                      0.,            #11: Continuum linear component
                      0.,            #12: Continuum quadratic component
                      IPpars[1],     #13: Insrumental resolution linear component
                      IPpars[0]])    #14: Insrumental resolution quadratic component

    # Make sure data is within telluric template range (shouldn't do anything)
    a0fluxlist = a0fluxlist[(a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0u        = a0u[       (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0x        = a0x[       (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    continuum  = continuum[ (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0wavelist = a0wavelist[(a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]

    # Define main spectrum
    s = a0fluxlist.copy(); x = a0x.copy(); u = a0u.copy();

    # Collect all fit variables into one class
    fitobj = fitobjs(s, x, u, continuum, watm_in, satm_in, mflux_in, mwave_in, [])

    # Arrays defining parameter variations during optimization steps
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0])
    dpar      = np.array([0.0, 0.0, 0.0, 3.0, 0.0, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0])
    dpar_st   = np.array([0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0])
    dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0])

    #-------------------------------------------------------------------------------

    # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
    optimize = True
    par_in = parA0.copy()
    hardbounds = [par_in[4] -dpar[4],   par_in[4]+dpar[4],
                  par_in[5] -dpar[5],   par_in[5]+dpar[5]]
    if hardbounds[0] < 0:
        hardbounds[0] = 0
    if hardbounds[3] < 0:
        hardbounds[3] = 1

    # Begin optimization.
    # For every pre-Telfit spectral fit, first fit just template strength/rv/continuum, then just wavelength solution, then template/continuum again, then ip,
    # then finally wavelength. Normally would fit for all but wavelength at the end, but there's no need for the pre-Telfit fit, since all we want
    # is a nice wavelength solution to feed into Telfit.

    pre_err = False
    try:
        parfit_1 = optimizer(par_in,   dpar_st,   hardbounds, fitobj, optimize)
        parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds, fitobj, optimize)
        parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds, fitobj, optimize)
        parfit_4 = optimizer(parfit_3, dpar,      hardbounds, fitobj, optimize)
        parfit = optimizer(parfit_4,   dpar_wave, hardbounds, fitobj, optimize)
    except:
        pre_err = True
        logger.warning(f'  --> NIGHT {night}, ORDER {order} HIT ERROR DURING PRE_OPT')
        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if order == firstorder:
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
        else:
            hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)


    #-------------------------------------------------------------------------------
    if not pre_err:
        # Get best fit wavelength solution
        a0w_out_fit = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)

        # Trim stellar template to new relevant wavelength range
        mwave_in,mflux_in = stellarmodel_setup(a0w_out_fit/1e4, inparam.mwave0, inparam.mflux0)

        # Feed this new wavelength solution into Telfit. Returns high-res synthetic telluric template, parameters of that best fit, and blaze function best fit
        watm1, satm1, telfitparnames, telfitpars, a0contwave, continuum = telfitter(a0w_out_fit,a0fluxlist,a0u,inparam,night,order,args)
    else:
        pass
    #-------------------------------------------------------------------------------

    # If Telfit encountered error (details in Telfitter.py), skip night/order combo
    if pre_err:
        pass
    elif len(watm1) == 1:
        logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}')

        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if order == firstorder:
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
        else:
            hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)

    else: # If Telfit exited normally, proceed.
        #  Save best blaze function fit
        a0contwave /= 1e4
        continuum = rebin_jv(a0contwave,continuum,a0wavelist,False)

        # Fit the A0 again using the new synthetic telluric template.
        # This allows for any tweaks to the blaze function fit that may be necessary.
        fitobj = fitobjs(s, x, u, continuum,watm1,satm1,mflux_in,mwave_in,[])

        try:
            parfit_1 = optimizer(par_in,   dpar_st,   hardbounds, fitobj, optimize)
            parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds, fitobj, optimize)
            parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds, fitobj, optimize)
            parfit_4 = optimizer(parfit_3, dpar_wave, hardbounds, fitobj, optimize)
            parfit   = optimizer(parfit_4, dpar,      hardbounds, fitobj, optimize)
        except:
            pre_err = True
            logger.warning(f'  --> NIGHT {night}, ORDER {order} HIT ERROR DURING POST_OPT')
            # Write out table to fits header with errorflag = 1
            c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
            cols  = fits.ColDefs([c0])
            hdu_1 = fits.BinTableHDU.from_columns(cols)

            # If first time writing fits file, make up filler primary hdu
            if order == firstorder:
                bleh = np.ones((3,3))
                primary_hdu = fits.PrimaryHDU(bleh)
                hdul = fits.HDUList([primary_hdu,hdu_1])
                hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
            else:
                hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
                hh.append(hdu_1)
                hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)

        if not pre_err:
            if inparam.plotfigs: # Plot results
                outplotter_tel(parfit, fitobj, f'FinalFit_Order{order}_{night}', inparam, args)

            if args.debug: # Output debug stuff
                fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)
                axes.plot(fitobj.x,
                          parfit[5] + parfit[13]*(fitobj.x) + parfit[14]*(fitobj.x**2),
                          '-k', lw=0.7, label='data', alpha=.6)

                axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
                axes.set_ylabel('IP',           size=6, style='normal', family='sans-serif')
                axes.set_xlabel('Pixel Number', size=6, style='normal', family='sans-serif')
                axes.legend(fontsize=5, edgecolor='white')
                fig.savefig(f'{inparam.outpath}/figs_{args.band}/IP_{order}_{night}.png',
                            bbox_inches='tight', format='png', overwrite=True)

                outplotter_tel(parfit_1,fitobj, f'Post_parfit1_{order}_{night}', inparam, args)
                outplotter_tel(parfit_2,fitobj, f'Post_parfit2_{order}_{night}', inparam, args)
                outplotter_tel(parfit_3,fitobj, f'Post_parfit3_{order}_{night}', inparam, args)
                outplotter_tel(parfit_4,fitobj, f'Post_parfit4_{order}_{night}', inparam, args)
                outplotter_tel(parfit, fitobj, f'Post_parfit_{order}_{night}', inparam, args)

            logger.debug(f'Post_par_in:\n {par_in}')
            logger.debug(f'Post_parfit1:\n {parfit_1}')
            logger.debug(f'Post_parfit2:\n {parfit_2}')
            logger.debug(f'Post_parfit3:\n {parfit_3}')
            logger.debug(f'Post_parfit4:\n {parfit_4}')
            logger.debug(f'Post_parfit:\n {parfit}')

            #-------------------------------------------------------------------------------

            # Save slightly better wavelength and blaze function solution
            a0w_out  = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)
            cont_adj = parfit[10] + parfit[11]*x + parfit[12]*(x**2.)

            continuum /= np.median(continuum)
            cont_save = continuum*cont_adj

            # Write out table to fits file with errorflag = 0
            c0  = fits.Column(name='ERRORFLAG'+str(order),      array=np.array([0]),            format='K')
            c1  = fits.Column(name='WAVE'+str(order),           array=a0w_out,                  format='D')
            c2  = fits.Column(name='BLAZE'+str(order),          array=cont_save,                format='D')
            c3  = fits.Column(name='X'+str(order),              array=a0x,                      format='D')
            c4  = fits.Column(name='INTENS'+str(order),         array=a0fluxlist,               format='D')
            c5  = fits.Column(name='SIGMA'+str(order),          array=a0u,                      format='D')
            c6  = fits.Column(name='WATM'+str(order),           array=watm1,                    format='D')
            c7  = fits.Column(name='SATM'+str(order),           array=satm1,                    format='D')
            c8  = fits.Column(name='TELFITPARNAMES'+str(order), array=np.array(telfitparnames), format='8A')
            c9  = fits.Column(name='TELFITPARS'+str(order),     array=telfitpars,               format='D')
            c10 = fits.Column(name='PARFIT',                    array=parfit,                   format='D')
            cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
            hdu_1 = fits.BinTableHDU.from_columns(cols)

            if order == firstorder: # If first time writing fits file, make up filler primary hdu
                bleh = np.ones((3,3))
                primary_hdu = fits.PrimaryHDU(bleh)
                hdul = fits.HDUList([primary_hdu,hdu_1])
                hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)
            else:
                hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
                hh.append(hdu_1)
                hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)
        else:
            pass
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def mp_run(args, inparam, Nthreads, jerp, orders, nights):
    # Multiprocessing convenience function
    pool = mp.Pool(processes = Nthreads)
    func = partial(MPinst, args, inparam, jerp, orders)
    outs = pool.map(func, np.arange(len(nights)))
    pool.close()
    pool.join()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def use_w(args):
    # Load wavelength regions list file
    try:
        bounddata = Table.read(f'./Input/UseWv/WaveRegions_{args.WRegion}_{args.band}.csv', format='csv')
    except IOError:
        sys.exit(f'WaveRegions FILE "./Input/UseWv/WaveRegions_{args.WRegion}_{args.band}.csv" NOT FOUND!')

    wavesols = pd.read_csv(f'./Input/UseWv/WaveSolns_{args.band}.csv')
#-------------------------------------------------------------------------------
    with open(f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv','w') as filew:
        filew.write('order, start,  end, masks\n')

        m_order  = np.array(bounddata['order'])
        starts   = np.array(bounddata['start'])
        ends     = np.array(bounddata['end'])
        ords     = list( sorted(orderdict_cla().orderdict[args.band].keys()) )

        Ostarts  = [orderdict_cla().orderdict[args.band][k][0] for k in ords]
        Oends    = [orderdict_cla().orderdict[args.band][k][1] for k in ords]
        labels   = []

        m_orders_unique = np.unique(m_order)

        # For each order specified, find what pixel numbers correspond to the wavelength bounds presented.
        # If multiple wavelength bounds given for a single order, output a pixel mask between the two, as well.
        for o in range(len(m_orders_unique)):
            pixs = [];
            mini = np.where(m_order == m_orders_unique[o])[0]
            for j in range(len(mini)):
                i = mini[j]

                wavebounds = [starts[i],ends[i]]
                wO   = wavesols['w'+str(m_orders_unique[o])]
                pixO = wavesols['x'+str(m_orders_unique[o])];
                pix  = [pixO[(np.argmin(abs(wO-wavebounds[k])))] for k in [0,1]]
                pixs = pixs + pix
                
            pixsS = list(sorted(pixs))
            q = pixsS[1:-1]
            if len(pixsS) == 2:
                filew.write('{}, {}, {},[]\n'.format(m_orders_unique[o],pixsS[0],pixsS[-1]))
            else:
                filew.write('{}, {}, {},"{}"\n'.format(m_orders_unique[o],pixsS[0],pixsS[-1],[[first,second] for first,second in  zip(q[0::2], q[1::2])]))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 1',
                                     description = '''
                                     This step 1) defines the wavelength regions to be analyzed based on user specification \n
                                     2) generates a synthetic, high-resolution telluric
                                     template for use in later model fits on a night by night basis.  \n
                                     Note that only target star observations will have their fits limited to the wavelength regions specified. \n
                                     For A0 observations, only the orders specified will be analyzed, but each order will be fit as far as there is significant telluric absoprtion.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name, no space",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",        action="store_true",
                        help="If set, will generate plots of A0 fitting results under ./Output/A0Fits/*target/fig/")

    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots under ./Temp/Debug/*target_*band/main_step1")
    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.9')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'
#-------------------------------------------------------------------------------
    # Create output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}')

    outpath = f'./Output/{args.targname}_{args.band}/A0Fits'
#-------------------------------------------------------------------------------
    # Handle logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}_A0Fits.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('####################################################################################\n')
    print(f'Fetching Wavelength Regions to be Analyzed for {args.targname}...')
    time.sleep(2)

    use_w(args)

    print('Fetching Done!')
    print(f'File "XRegions_{args.WRegion}_{args.band}.csv" saved under "./Input/UseWv/"')
    time.sleep(5)

    #-------------------------------------------------------------------------------

    print('###############################################################\n')
    logger.info(f'Using TelFit to create high-resolution, synthetic telluric templates based off the telluric standards \nassociated with {args.targname} on a night by night basis...')
    print('This will take a while..........')
    print('\n')

    # Read in newly created pixel regions file to get list of orders to analyze.
    # Note that only target star observations will have their fits limited to the wavelength regions specified.
    # For A0 observations, only the orders specified will be analyzed, but each order will be fit as far as there is significant telluric absoprtion.
    bounddata = Table.read(f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv', format='csv')
    starts  = np.array(bounddata['start'])
    ends    = np.array(bounddata['end'])
    orders  = np.array(bounddata['order'], dtype=int)
    xbounddict = {orders[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}

    ## Collect relevant file information from Predata files
    A0data = Table.read(f'./Input/Prepdata/Prepdata_A0_{args.targname}.txt', format='ascii')

    ind    = [i != 'NA' for i in A0data['humid']]
    humids = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['humid'])}
    tags   = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['tag'])}
    obs    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['obs'])}
    temps  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['temp'])}
    zds    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['zd'])}
    press  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['press'])}
    nightsFinal = np.array(list(sorted(set(A0data[ind]['night']))))

    # Take subset of nights, if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=int)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(f'NIGHT {nnn} EITHER HAS NO CORRESPONDING A0 OR WAS NOT FOUND UNDER "./Input/{args.targname}"')

        nightsFinal = nightstemp
        logger.info(f'Only processing nights: {nightsFinal}')

    logger.info(f'Analyze {len(nightsFinal)} nights')

    time.sleep(6)
    print('\n')

    # print('For paper plot!')
    # if args.band == 'K':
    #     orders = np.array([2, 3, 4, 5, 6, 7, 8, 10, 14, 16])
    # elif args.band=='H':
    #     orders = np.array([2, 3, 4, 6, 13, 14, 16, 20, 21, 22])
    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates_tel()

    inparam = inparamsA0(inpath,outpath,args.plotfigs,tags,nightsFinal,humids,
                         temps,zds,press,obs,watm,satm,mwave0,mflux0,cdbs_loc,xbounddict,None)

    #-------------------------------------------------------------------------------

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        outs = mp_run(args, inparam, args.Nthreads, jerp, orders, nightsFinal)

    print('\n')
    logger.info('A0 Fitting Done!')

    end_time = datetime.now()
    logger.info(f'A0 Fitting using TelFit finished, Duration: {end_time - start_time}')
    print('You can start to run main_step2.py for RV initial guess')
    print('####################################################################################')
