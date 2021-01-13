from Engine.importmodule import *
from Engine.importmodule import read_prepdata

from Engine.IO_AB      import setup_templates, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs,inparams
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod, fmod_conti
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.crmask    import CRmasker

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 3',
                                     description = '''
                                     Performs a full analysis of each target star observation to produce accurate and precise RVs. \n
                                     All the wavelength regions defined in Step 1 are used, and the code analyzes each observation that is part of a given exposure separately. \n
                                     Unless the target vsini is already known to high accuracy, an initial run of Step 3 in which \vsini is allowed to vary is required. \n
                                     This provides an estimate of vsini that can then be plugged into the code as a fixed value in the second run of Step 3. \n
                                     If the user seeks the best possible RV uncertainty estimates, or if their target star has a relatively high \vsini ($>$ 10 \kms), they must run Step 3 once with \vsini held fixed at its estimated value and once with \vsini held fixed at this value plus or minus one sigma. \n
                                     The minor differences in the RVs of the two runs (as low as $<$1 \ms and as high as 7 \ms) can then be incorporated into the final uncertainties. \n
                                     If \vsini is already well-known, it is not necessary to run Step 3 more than once, as the code fully converges to the final RVs (within uncertainty) through just one run.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Run",    dest="run_num",             action="store",
                        help="Which run number to process? For get Telluric_corrected_X",
                        type=str,   default='0')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If set, will generate plots of the fitting results under ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.9')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    # Check user input

    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Band                = \33[41m {} band \033[0m
    WaveLength file     = \33[41m WaveRegions_{} \033[0m
    Run #               = \33[41m {} \033[0m
    '''.format(args.targname, args.band, args.WRegion, args.run_num))
    if not args.skip:
        while True:
            inpp = input("Press [Y]es to continue, [N]o to quite...\n --> ")
            if 'n' in inpp.lower():
                sys.exit('QUIT, PLEASE RE-ENTER YOUR PARAMETERS')
            elif 'y' in inpp.lower():
                break
            else:
                print('I cannot understand what you are saying... TRY AGAIN')
                continue

    print('---------------------------------------------------------------')
    print('Running Step 2.5 for {}...'.format(args.targname))
    print('This will take a while..........')

    #-------------------------------------------------------------------------------

    # Make output directories as needed

    name = f'TelluricCorrected_results_{args.run_num}'

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/{name}/Generated_Template/'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/{name}/Generated_Template/')

    outpath = f'./Output/{args.targname}_{args.band}/{name}/'

    #-------------------------------------------------------------------------------

    # Set up logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)

    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{name}/Generated_Template/')

    #-------------------------------------------------------------------------------

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))

    logger.info('Analyze with {} nights'.format(len(nightsFinal)))

    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T).
    # All statistical analysis will be performed separately for these two datasets.
    nights    = nightsFinal.copy()
    intnights = np.array([int(i[:8]) for i in nights])

    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these
                nights slightly differently.
When you run Step 3, RVs will be output in two formats:
                one with the defocus nights separated,
                and the other with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as
                the result of the defocus. IGRINS RV will process defocus nights
                the same way as the others, but when you run Step 3, will still
                output the results in two formats like it does with the K band.''')

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]

    if len(nightsL) > 0:
        nightscomblist = [nightsT,nightsL]
    else:
        nightscomblist = [nightsT]

    #-------------------------------------------------------------------------------

    T_L = 'T'
    c = 2.99792458e5
    firstorder = orders[0]

    for nights_use in nightscomblist:

        # Run order by order
        for jerp in range(len(orders)):

            order = orders[jerp]
            firstgo = True

            wmins = []; wmaxs = []; ntags = [];

            for night in nights_use:
                
                tags = np.concatenate((np.array(tagsA[night]), np.array(tagsB[night])))
                
                for tag in tags:
                
                    hdu = fits.open('{}/Stellar_Residual_{}_{}_{}.fits'.format(outpath, args.band, night, tag))
                    tbdata = hdu[jerp+1].data

                    wave_in  = np.array(tbdata['WAVE'+str(order)])
                    bvc      = np.array(tbdata['BVC'])[0]

                    wave_corr = wave_in * (1 + bvc/c)

                    wmins.append(wave_corr[0]); wmaxs.append(wave_corr[-1]);
                    ntags.append([night,tag])
            
            # Rebin to the wavelength scale that is middlemost of the possible ranges
            wmins = np.array(wmins); wmaxs = np.array(wmaxs); ntags = np.array(ntags);
            dist = np.sqrt( (wmins-np.min(wmins))**2 + (wmaxs-np.max(wmaxs))**2 )
            masternight = ntags[np.argmin(dist)][0]; mastertag = ntags[np.argmin(dist)][1];
   

            hdu = fits.open('{}/Stellar_Residual_{}_{}_{}.fits'.format(outpath, args.band, masternight, mastertag))
            tbdata = hdu[jerp+1].data
            masterflux  = np.array(tbdata['STELL'+str(order)])
            masterunc   = np.array(tbdata['UNC'+str(order)])
            wave_in     = np.array(tbdata['WAVE'+str(order)])
            bvc         = np.array(tbdata['BVC'])[0]
            masterwave  = wave_in * (1 + bvc/c)

            if args.plotfigs:
                plt.figure(figsize=(16,10))
                plt.plot(masterwave,masterflux,alpha=.2)

            for night in nights_use:
                
                tags = np.concatenate((np.array(tagsA[night]), np.array(tagsB[night])))
                
                for tag in tags:
                
                    if night == masternight and tag == mastertag:
                        continue

                    hdu = fits.open('{}/Stellar_Residual_{}_{}_{}.fits'.format(outpath, args.band, night, tag))
                    tbdata = hdu[jerp+1].data

                    wave_in  = np.array(tbdata['WAVE'+str(order)])
                    stell_in = np.array(tbdata['STELL'+str(order)])
                    unc_in   = np.array(tbdata['UNC'+str(order)])
                    bvc      = np.array(tbdata['BVC'])[0]

                    wave_corr = wave_in * (1 + bvc/c)

                    fluxbinned = rebin_jv(wave_corr,stell_in,masterwave,False)
                    ubinned    = rebin_jv(wave_corr,unc_in,masterwave,False)
      
                    # Throw out extrapolated values
                    fluxbinned[(masterwave < wave_corr[0]) | (masterwave > wave_corr[-1])] = np.nan
                    ubinned[(masterwave < wave_corr[0]) | (masterwave > wave_corr[-1])] = np.nan

                    masterflux = np.vstack((masterflux,fluxbinned))
                    masterunc  = np.vstack((masterunc,ubinned))
        
                    if args.plotfigs:
                        plt.plot(masterwave,fluxbinned,alpha=.2)

            flux_out = np.array([np.nansum( masterflux[:,i] * ((1./(masterunc[:,i]**2)) / (np.nansum(1./(masterunc[:,i]**2)))) ) for i in range(len(masterwave))])
        
            # trim by 5 *MAD ? shouldnt be necessary
            # MAD = np.median(abs(np.median(flux_out)-flux_out)) 

            if args.plotfigs:

                plt.savefig('{}/Generated_Template/StellarTemplate_Separates_{}_{}.png'.format(outpath, T_L, order))
                plt.clf()  
                plt.close()

                plt.figure(figsize=(16,10))
                plt.plot(masterwave,flux_out)           
                plt.savefig('{}/Generated_Template/StellarTemplate_Combined_{}_{}.png'.format(outpath, T_L, order))
                plt.clf()  
                plt.close()
    

            if jerp == 0:
                mwave = masterwave.copy()
                mflux = flux_out.copy()
            else:
                mwave = np.concatenate((masterwave,mwave))
                mflux = np.concatenate((flux_out,mflux))

        c1 = fits.Column(name='MWAVE_'+T_L,       array=mwave,                                 format='D')
        c2 = fits.Column(name='MFLUX_'+T_L,       array=mflux,                                 format='D')
        cols = fits.ColDefs([c1,c2])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if T_L == 'T': # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/Generated_Template/StellarTemplateFromData.fits'.format(outpath), overwrite=True)
        else:
            hh = fits.open('{}/Generated_Template/StellarTemplateFromData.fits'.format(outpath))
            hh.append(hdu_1)
            hh.writeto('{}/Generated_Template/StellarTemplateFromData.fits'.format(outpath), overwrite=True)

        T_L = 'L'

    print('\n')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}/{}'.format(outpath, name) )
    print('####################################################################################')
