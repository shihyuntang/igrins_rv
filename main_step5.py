
from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step5

from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.bisectorcalc    import NightSpecs,BIinst
#from Engine.LS         import LS
from Engine.plot_tool import modtool
from scipy.stats import pearsonr

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    args = _argparse_step5()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    if (args.temperature == '') & (args.logg == ''):
        sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
                    'STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE '
                    'AVAILABLE TEMPLATES')
    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'./Engine/syn_template')

    if args.template.lower() == 'synthetic':
        #list of all syntheticstellar
        syntemp = [i for i in syntemp if i[:3] == 'syn']
        synT    = [ i.split('_')[2][1:]  for i in syntemp ]
        synlogg = [ i.split('_')[3][4:7] for i in syntemp ]
    elif args.template.lower() == 'phoenix':
        #list of all phoenix
        syntemp = [i for i in syntemp if i[:3] == 'PHO']
        synT    = [ i.split('-')[1][4:]  for i in syntemp ]
        synlogg = [ i.split('-')[2][:3] for i in syntemp ]
    else:
        synT = [args.temperature]; synlogg = [args.logg]

    if args.temperature not in synT:
        sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp" INPUT! '
                    f'{syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    if args.logg not in synlogg:
        sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg" INPUT! {syntemp} '
                    'AVALIABLE UNDER ./Engine/syn_template/')

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[37;1;41m {} band \033[0m
    WaveLength file     = \33[37;1;41m WaveRegions_{} \033[0m
    RV run #            = \33[37;1;41m {} \033[0m
    Stellar template use= \33[37;1;41m {} \033[0m
    syn template temp   = \33[37;1;41m {} \033[0m
    syn template logg   = \33[37;1;41m {} \033[0m
    '''.format(args.targname, args.band, args.WRegion, args.run, args.template, args.temperature, args.logg))
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
    print('Running Step 5 for {}...'.format(args.targname))
    print('This will take a while..........')

    #-------------------------------------------------------------------------------

    name = f'Cutouts'

    outpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/{name}'

    if not os.path.isdir(f'{outpath}'):
        os.mkdir(f'{outpath}')
    if not os.path.isdir(f'{outpath}/figs'):
        os.mkdir(f'{outpath}/figs')
    #-------------------------------------------------------------------------------

    # Set up logger
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
    logger.propagate = False

    logger.info(f'Writing output to {outpath}')

    #-------------------------------------------------------------------------------

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(
                    'NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(
                        nnn, args.targname
                        ))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))


    logger.info('Analyze with {} nights'.format(len(nightsFinal)))

    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(
        logger, args.template, args.band, np.int(args.temperature),
        np.float(args.logg), np.float(args.B)
        )

    # Save pars in class for future use
    inparam = InParams(inpath, outpath, None, None, args.plotfigs,
                       None, bvcs, tagsA, tagsB, nightsFinal, mwave0,
                       mflux0, None, xbounddict, maskdict)

    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T).
    # All statistical analysis will be performed separately for these two datasets.
    nights    = inparam.nights
    intnights = np.array([int(i[:8]) for i in nights])

    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these nights
                slightly differently. When you run Step 3, RVs will be output in
                two formats: one with the defocus nights separated, and the other
                with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as the result
                of the defocus. IGRINS RV will process defocus nights the same way
                as the others, but when you run Step 3, will still output the results
                in two formats like it does with the K band.''')

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]

    if len(nightsL) > 0:
        nightscomblist = [nightsT,nightsL]
    else:
        nightscomblist = [nightsT]

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    if args.skipmod:
        pass
    else:
        print('\n \n Constructing residuals...\n \n ')
        rawbox = fits.open(f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/RVresultsRawBox.fits')
        for jerp in range(len(orders)):
            boxdata = rawbox[jerp+1].data
            order = orders[jerp]
            nightsbox = np.array(boxdata['NIGHT'+str(orders[jerp])])
            parfitbox = np.array(boxdata['PARFIT'+str(orders[jerp])])
            tagbox    = np.array(boxdata['TAG'+str(orders[jerp])])

            if not args.debug:
                print('Working on order {} ({:02d}/{:02d})'.format(
                    orders[jerp], int(jerp+1), len(orders)
                    ))

            #modtool(args,jerp,nightsbox,tagbox,parfitbox,inparam,0)
            func = partial(modtool,args,jerp,nightsbox,tagbox,parfitbox,inparam)
            outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

        print('\n \n Done making residuals! Now calculating bisectors...\n \n ')

    name = f'Bisectors'

    outpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/'
    figpath = f'./Output/{args.targname}_{args.band}/figs/main_step5_{args.band}_{args.run}'

    if not os.path.isdir(f'{figpath}'):
        os.mkdir(f'{figpath}')

    inpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/Cutouts'

    s2ns = np.ones((len(nightsFinal)),dtype=float)
    for i in range(len(nightsFinal)):
        nightspec = NightSpecs(inpath, nightsFinal[i], orders, 2)
        if nightspec.flag != 1:
            s2ns[i] = np.nanmedian(nightspec.flux/nightspec.unc)
    refnight = nightsFinal[np.argmax(s2ns)]

    bimasterboxT  = np.ones((len(nightsT),len(orders)))*np.nan
    stdmasterboxT = np.ones((len(nightsT),len(orders)))*np.nan
    bimasterboxL  = np.ones((len(nightsL),len(orders)))*np.nan
    stdmasterboxL = np.ones((len(nightsL),len(orders)))*np.nan

    if len(nightsL) > 0 and len(nightsT) > 0:
        nightscomblist = [nightsT,nightsL]
        T_Ls = ['T','L']
    elif len(nightsL) > 0:
        nightscomblist = [nightsL]
        T_Ls = ['L']
    else:
        nightscomblist = [nightsT]
        T_Ls = ['T']

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        if not args.debug:
            print('Working on order {} ({:02d}/{:02d})'.format(
                orders[jerp], int(jerp+1), len(orders)
                ))

        refspec = NightSpecs(inpath, refnight, orders, jerp)
        order = orders[jerp]

        if refspec.flag == 1: #whole order absent
            nightsbox = nightsFinal.copy()
            bibox  = np.ones_like(nightsFinal,dtype=float)*np.nan
            stdbox = np.ones_like(nightsFinal,dtype=float)*np.nan
            logger.warning(f'WARNING! All of order {order} absent from residuals!')
        else:
            #BIinst(refspec, orders, jerp, nightsFinal, args,inpath,figpath,logger,0)
            func = partial(BIinst, refspec, orders, jerp, nightsFinal, args,inpath,figpath,logger)
            outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

            # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
            outs = np.array(outs)
            nightsbox = outs[:,0]
            bibox = outs[:,1]
            stdbox = outs[:,2]

        # Save results to fits file
        c1    = fits.Column(name='NIGHT'+str(order),  array=nightsbox, format='{}A'.format(len(nights[0])) )
        c2    = fits.Column(name='BI'+str(order),     array=bibox,     format='D')
        c3    = fits.Column(name='BISTD'+str(order),  array=stdbox,  format='D')
        cols  = fits.ColDefs([c1,c2,c3])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if jerp == 0: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/BIresultsRawBox.fits'.format(outpath),
                overwrite=True)
        else:
            hh = fits.open('{}/BIresultsRawBox.fits'.format(outpath))
            hh.append(hdu_1)
            hh.writeto('{}/BIresultsRawBox.fits'.format(outpath),
                overwrite=True)

        #-------------------------------------------------------------------------------

        # For each set of nights (tight, loose)...
        iT_L = 0
        for nights_use in nightscomblist:

            T_L = T_Ls[iT_L]

            # For each night...
            for i in range(len(nights_use)):
                # Collect the RVs and vsinis determined from different A/B
                # exposures within a night
                indnight  = np.where(nightsbox == nights_use[i])[0]
                bitags    =  bibox[indnight]
                stdtags   = stdbox[indnight]

                # Take the mean of the vsinis, and the mean and std of the RVs.
                # If the number of different successfully fitted A/B exposures
                # is less than required, pass NaN instead.
                if T_L == 'T':
                    bimasterboxT[i,jerp]  = bitags[0]
                    stdmasterboxT[i,jerp] = stdtags[0]

                else:
                    bimasterboxL[i,jerp]  = bitags[0]
                    stdmasterboxL[i,jerp] = stdtags[0]
            iT_L += 1


    #-------------------------------------------------------------------------------
    if not args.debug: logger.addHandler(stream_hander)
    print('\n')

    # Don't combine Loose and Tight datasets, but make them both easily referenceable
    jdsCombined  = np.array([])
    nightsCombined  = np.array([])
    bifinalCombined = np.array([])
    stdfinalCombined = np.array([])
    rvfinalCombined = np.array([])
    rvstdfinalCombined = np.array([])

    if T_Ls == ['T','L']:
        rvboxcomblist  = [bimasterboxT,bimasterboxL]
        stdboxcomblist = [stdmasterboxT,stdmasterboxL]
    elif T_Ls == ['L']:
        rvboxcomblist  = [bimasterboxL]
        stdboxcomblist = [stdmasterboxL]
    else:
        rvboxcomblist  = [bimasterboxT]
        stdboxcomblist = [stdmasterboxT]

    hduRV = fits.open(f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/RVresultsSummary.fits')
    tbdataRV = hduRV[1].data
    nightsRV = np.array(tbdataRV['NIGHT'],dtype=str)
    jd0      = np.array(tbdataRV['JD'],dtype=str)
    rvfinal0      = np.array(tbdataRV['RVfinal'],dtype=str)
    rvstdfinal0   = np.array(tbdataRV['STDfinal'],dtype=str)

    # Iterate over tight and loose mounting data sets...
    for boxind in range(len(rvboxcomblist)):

        bimasterbox  = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]

        #-------------------------------------------------------------------------------

        if args.mode=='STD': # If RV STD star, calculate uncertainty in method

            # Calculate the precision within an order across nights
            sigma_O2     = np.array(
                [np.nanstd(bimasterbox[:,ll])**2 for ll in range(len(orders))]
                )
            sigma_ABbar2 = np.ones_like(sigma_O2)

            # Calculate uncertainty in method as difference between variance
            # within an order and mean variance within a night's different exposure RVs
            for ll in range(len(orders)):
                sigma_ABbar2[ll] = np.nanmedian(stdmasterbox[:,ll]**2)
            sigma_method2 = sigma_O2 - sigma_ABbar2

        # If target star, load the uncertainty in method calculated from our
        # RV STD star runs
        else:
            if T_Ls[boxind] == 'T':
                nights_use = nightsT.copy()
                kind = 'Focused'
                sigma_method2 = inparam.BImethodvariance_tight[args.band]
            elif T_Ls[boxind] == 'L':
                nights_use = nightsL.copy()
                kind = 'Defocus'
                sigma_method2 = inparam.BImethodvariance_loose[args.band]

        sigma_ON2    = np.ones_like(bimasterbox)

        #-------------------------------------------------------------------------------

        # Note bimasterbox indexed as [nights,orders]
        Nnights = len(bimasterbox[:,0])

        for ll in range(len(orders)):
            # Calculate the uncertainty in each night/order RV as the sum of the
            # uncertainty in method and the uncertainty in that night's As and Bs RVs
            for night in range(Nnights):
                sigma_ON2[night,ll] = sigma_method2[ll] + stdmasterbox[night,ll]**2

        bifinal    = np.ones(Nnights, dtype=np.float64)
        stdfinal   = np.ones(Nnights, dtype=np.float64)

        if T_Ls[boxind] == 'T':
            nights_use = nightsT.copy(); kind = 'Focused'
        elif T_Ls[boxind] == 'L':
            nights_use = nightsL.copy(); kind = 'Defocused'

        # Combine RVs between orders using weights calculated from uncertainties
        for n in range(Nnights):
            ind = np.where(
                np.isfinite(sigma_ON2[n,:]) & np.isfinite(bimasterbox[n,:]))[0]
            weights = (1./sigma_ON2[n,ind]) / (np.nansum(1./sigma_ON2[n,ind])) # normalized
            stdspre = (1./sigma_ON2[n,ind]) #unnormalized weights

            bifinal[n]  = np.nansum( weights*bimasterbox[n,ind] )
            stdfinal[n] = 1/np.sqrt(np.nansum(stdspre))

            # if all the RVs going into the observation's final RV calculation
            # were NaN due to any pevious errors, pass NaN
            if np.nansum(weights) == 0:
                bifinal[n]    = np.nan
                stdfinal[n]   = np.nan

            # if more than half of the orders going into the observation's final
            # RV calculation were NaN due to any pevious errors, pass NaN
            if np.sum( np.isnan(bimasterbox[n,:]) ) > np.floor( len(orders) * 0.5 ):
                bifinal[n]    = np.nan
                stdfinal[n]   = np.nan

        #-------------------------------------------------------------------------------

        # Plot results
        f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)

        axes.plot(    np.arange(len(bifinal))+1, bifinal, '.k', ms=5)
        axes.errorbar(np.arange(len(bifinal))+1, bifinal, yerr=stdfinal,
            ls='none', lw=.5, ecolor='black')
        axes.text(0.05, 0.93, r'BI mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(
            np.nanmean(bifinal), np.nanstd(bifinal)),
            transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        #axes.set_ylim(np.nanmin(bifinal)-.08, np.nanmax(bifinal)+.08)
        axes.set_ylim(-0.3,0.3)
        axes.set_ylabel('Bisector Span [km/s]', size=6, style='normal', family='sans-serif' )
        axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(5))
        axes.tick_params(axis='both', which='both', labelsize=5, right=True,
            top=True, direction='in', width=.6)
        f.savefig('{}/FinalBIs_{}_.png'.format(outpath, kind),
            format='png', bbox_inches='tight')


        jdfinal = np.ones_like(bifinal)*np.nan
        rvfinal = np.ones_like(bifinal)*np.nan
        rvstdfinal = np.ones_like(bifinal)*np.nan
        for aa in range(len(nights_use)):
            ind = np.where(nightsRV == nights_use[aa])[0]
            rvfinal[aa] =  rvfinal0[ind[0]]
            rvstdfinal[aa] = rvstdfinal0[ind[0]]
            jdfinal[aa] = jd0[ind[0]]


        f, axes = plt.subplots(1, 1, figsize=(8,6), facecolor='white', dpi=300)

        mask = np.ones_like(bifinal,dtype=bool)
        mask[np.isnan(bifinal) | np.isnan(rvfinal)] = False
        pp = pearsonr(rvfinal[mask],bifinal[mask])

        axes.plot(    rvfinal, bifinal, '.k', ms=5)
        axes.errorbar(rvfinal, bifinal, xerr=rvstdfinal,yerr=stdfinal, ls='none', lw=.5, ecolor='black')
        axes.text(0.05, 0.93, 'Correlation Coeff = {} , P = {}'.format(round(pp[0],5),round(pp[1],5)),
                             transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        #axes.set_ylim(np.nanmin(bifinal)-.08,
        #             np.nanmax(bifinal)+.08)
        axes.set_ylim(-0.3,0.3)
        axes.set_xlabel('RV [km/s]', size=6, style='normal', family='sans-serif' )
        axes.set_ylabel('Bisector Span [km/s]', size=6, style='normal', family='sans-serif' )
        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(5))
        axes.tick_params(axis='both', which='both', labelsize=5, right=True, top=True, direction='in', width=.6)
        f.savefig('{}/FinalRVs_vs_BIs_{}_.png'.format(outpath, kind), format='png', bbox_inches='tight')

        # Save results to fits file separately for each tight/loose dataset
        c1 = fits.Column( name='NIGHT',         array=nights_use,    format='8A')
        c2 = fits.Column( name='JD',            array=jdfinal,       format='D')
        c3 = fits.Column( name='RVBOX',         array=bimasterbox,   format='{}D'.format(len(orders)))
        c4 = fits.Column( name='STDBOX',        array=stdmasterbox,  format='{}D'.format(len(orders)))
        c7 = fits.Column( name='Sigma_method2', array=sigma_method2, format='D')
        c8 = fits.Column( name='Sigma_ON2',     array=sigma_ON2,     format='{}D'.format(len(orders)))
        c9 = fits.Column( name='BIfinal',       array=bifinal,       format='D')
        c10 = fits.Column(name='BISTDfinal',      array=stdfinal,      format='D')
        c11 = fits.Column(name='RVfinal',         array=rvfinal,    format='D')
        c12 = fits.Column(name='STDfinal',         array=rvstdfinal,    format='D')

        if args.mode=='STD':
            c5 = fits.Column( name='Sigma_O2',      array=sigma_O2,      format='D')
            c6 = fits.Column( name='Sigma_ABbar2',  array=sigma_ABbar2,  format='D')
            cols  = fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
        else:
            cols  = fits.ColDefs([c1,c2,c3,c4,c7,c8,c9,c10,c11,c12])

        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul        = fits.HDUList([primary_hdu,hdu_1])
        if len(T_Ls) == 2:
            hdul.writeto(
                '{}/BIresultsAdvanced_{}.fits'.format(
                    outpath, kind),
                overwrite=True)
        else:
            hdul.writeto(
                '{}/BIresultsAdvanced.fits'.format(
                    outpath),
                overwrite=True)

        # Combine final RVs from both tight and loose mounting data sets
        nightsCombined     = np.concatenate((nightsCombined,     nights_use))
        jdsCombined        = np.concatenate((jdsCombined,        jdfinal))
        bifinalCombined    = np.concatenate((bifinalCombined,    bifinal))
        stdfinalCombined   = np.concatenate((stdfinalCombined,   stdfinal))
        rvfinalCombined    = np.concatenate((rvfinalCombined,    rvfinal))
        rvstdfinalCombined   = np.concatenate((rvstdfinalCombined,   rvstdfinal))

        if args.mode=='STD': # If uncertainty in method was calculated, save it
            sigma_method2 = [np.around(float(i), 8) for i in sigma_method2]
            logger.info('sigma_method2 during the {} epoch is {}'.format(kind, sigma_method2))
        if len(T_Ls) == 2:
            logger.info('During the {} epoch: BI mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                kind,
                np.nanmean(bifinal),
                np.nanstd(bifinal) ))
        else:
            logger.info('BI mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                np.nanmean(bifinal),
                np.nanstd(bifinal) ))

    #-------------------------------------------------------------------------------

    # Plot combined results
    xscale = np.arange(len(bifinalCombined))+1

    f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(xscale,bifinalCombined, '.k', ms=5)
    axes.errorbar(xscale,bifinalCombined, yerr=stdfinalCombined,
        ls='none', lw=.5, ecolor='black')
    axes.text(0.05, 0.93, r'BI mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(
        np.nanmean(bifinalCombined), np.nanstd(bifinalCombined)),
        transform=axes.transAxes, size=6, style='normal', family='sans-serif' )

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6,
            style='normal', family='sans-serif' )
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        axes.text(0.05, 0.1, 'Defocus', transform=axes.transAxes, size=6,
            style='normal', family='sans-serif' )
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            axes.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6,
                style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6,
                style='normal', family='sans-serif' )
        else:
            axes.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6,
                style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocused', transform=axes.transAxes, size=6,
                style='normal', family='sans-serif' )
    #axes.set_ylim(np.nanmin(bifinalCombined)-.08,np.nanmax(bifinalCombined)+.08)
    axes.set_ylim(-0.3,0.3)
    axes.set_ylabel('Bisector Span (km/s)', size=6, style='normal', family='sans-serif' )
    axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True,
        direction='in', width=.6)
    f.savefig('{}/FinalBIs.png'.format(outpath), format='png',
        bbox_inches='tight')


    f, axes = plt.subplots(1, 1, figsize=(8,6), facecolor='white', dpi=300)

    mask = np.ones_like(bifinalCombined,dtype=bool)
    mask[np.isnan(bifinalCombined) | np.isnan(rvfinalCombined)] = False
    pp = pearsonr(rvfinalCombined[mask],bifinalCombined[mask])

    axes.plot(    rvfinalCombined, bifinalCombined, '.k', ms=5)
    axes.errorbar(rvfinalCombined, bifinalCombined, xerr=rvstdfinalCombined,yerr=stdfinalCombined, ls='none', lw=.5, ecolor='black')
    axes.text(0.05, 0.93, 'Correlation Coeff = {} , P = {}'.format(round(pp[0],5),round(pp[1],5)),
                         transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    #axes.set_ylim(np.nanmin(bifinalCombined)-.08,
    #             np.nanmax(bifinalCombined)+.08)
    axes.set_ylim(-0.3,0.3)
    axes.set_xlabel('RV [km/s]', size=6, style='normal', family='sans-serif' )
    axes.set_ylabel('Bisector Span [km/s]', size=6, style='normal', family='sans-serif' )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=5, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/FinalRVs_vs_BIs_.png'.format(outpath), format='png', bbox_inches='tight')

    # Output combined final results to fits file
    c1 = fits.Column(name='NIGHT',    array=nightsCombined,     format='{}A'.format(len(nights[0])) )
    c2 = fits.Column(name='JD',       array=jdsCombined,        format='D')
    c3 = fits.Column(name='BIfinal',  array=bifinalCombined,    format='D')
    c4 = fits.Column(name='BISTDfinal', array=stdfinalCombined,   format='D')
    c5 = fits.Column(name='RVfinal',  array=rvfinalCombined,    format='D')
    c6 = fits.Column(name='STDfinal', array=rvstdfinalCombined,   format='D')

    cols = fits.ColDefs([c1,c2,c3,c4,c5,c6])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    bleh = np.ones((3,3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu,hdu_1])
    hdul.writeto(
        '{}/BIresultsSummary.fits'.format(outpath),
        overwrite=True)

    tempin = Table.read(
        '{}/BIresultsSummary.fits'.format(outpath),
        format='fits')
    tempin.write(
        '{}/BIresultsSummary.csv'.format(outpath),
        format='csv', overwrite=True)

    if len(T_Ls) == 2:
        logger.info('Combined BI results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(
            np.nanmean(bifinalCombined),
            np.nanstd(bifinalCombined)))

    warning_r = log_warning_id(f'{outpath}/Cutouts/{args.targname}_{args.band}.log', start_time)
    if warning_r:
        print(f'''
    **********************************************************************************
    WARNING!! you got warning message during this run. Please check the log file under:
          {outpath}/Cutouts/{args.targname}_{args.band}_A0Fits.log
    **********************************************************************************
    ''')
    print('\n')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}'.format(outpath) )
    logger.info('The final BI estimates you are looking for are in the BIresultsSummary files!')
    print('####################################################################################')
