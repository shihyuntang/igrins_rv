
from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step5

from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.fakebisectorcalc    import FakeSpec,BIinst
#from Engine.LS         import LS
from Engine.plot_tool import modtool
from scipy.stats import pearsonr

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    #-------------------------------------------------------------------------------

    targname = 'DITau'
    band = 'K'
    mode = 'STD'
    debug = True
    run = 2


    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T).
    # All statistical analysis will be performed separately for these two datasets.
    hduRV = fits.open(f'./Output/DITau_K/RV_results_2/RVresultsSummary.fits')
    tbdataRV = hduRV[1].data
    nightsRV = np.array(tbdataRV['NIGHT'],dtype=str)
    jd0      = np.array(tbdataRV['JD'],dtype=str)
    rvfinal1      = np.array(tbdataRV['RVfinal'],dtype=float)
    rvstdfinal1   = np.array(tbdataRV['STDfinal'],dtype=float)

    tableBI = Table.read(f'./Output/DITau_K/RV_results_2/BIresultsSummary_combined.csv',format='csv')
    bifinal1 = np.array(tableBI['BIfinal'],dtype=float)
    bistdfinal1 = np.array(tableBI['BISTDfinal'],dtype=float)

    rvfinal2  = np.array([11.197363242576982,13.522548777166282,13.122950818323005,
                        1.614368917053188,-0.24032865381311275])
    rvfinal2 += 15.556198659361755

    pow2 = 0.9

    #rvref1 = rvfinal1[1]
    #rvref2 = rvfinal2[1]

    nightsFinal = nightsRV.copy()
    nights = nightsFinal.copy()

    intnights = np.array([int(i[:8]) for i in nights])

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
    print('\n')

    rawbox = fits.open(f'./Output/DITau_K/RV_results_2/RVresultsRawBox.fits')

    #kinds = ['M2obs','M3obs','M5obs','M6obs','3500_4p0','3500_4p5','3000_4p0','3000_4p0_phx','3000_4p5_phx']
    #T2s   = [3600,3500,3200,3100,3500,3500,3000,3000,3000]

    #for k in range(len(kinds)):
    #    kind = kinds[k]; T2 = T2s[k]

    rvfinal20 = rvfinal2.copy()
    print(rvfinal20)
    for tt in np.arange(-25,35,10):
        kind = f'M5obs_RV{tt}'
        T2 = 3200.
        rvfinal2 = rvfinal20.copy() + tt
        print(rvfinal2)

        name = f'FakeBisectors_{kind}'

        outpath = f'./Output/{targname}_{band}/RV_results_{run}/{name}'
        figpath = f'./Output/{targname}_{band}/RV_results_{run}/{name}/figs'


        if not os.path.isdir(f'{outpath}'):
            os.mkdir(f'{outpath}')
        if not os.path.isdir(f'{outpath}/figs'):
            os.mkdir(f'{outpath}/figs')


        # Set up logger
        logger = logging.getLogger(__name__)
        if debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

        file_hander  = logging.FileHandler(f'{outpath}/{targname}_{band}.log')
        stream_hander= logging.StreamHandler()

        # file_hander.setLevel()
        file_hander.setFormatter(formatter)

        logger.addHandler(file_hander)
        logger.addHandler(stream_hander)
        logger.propagate = False

        logger.info(f'Writing output to {outpath}')

        orders = np.array([4,5,6])

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
        if not debug: logger.removeHandler(stream_hander)
        print('\n')

        # Run order by order, multiprocessing over nights within an order
        for jerp in range(len(orders)):
            if not debug:
                print('Working on order {} ({:02d}/{:02d})'.format(
                    orders[jerp], int(jerp+1), len(orders)
                    ))

            order = orders[jerp]

            if order == 4:
                continue

            boxdata = rawbox[jerp+1].data
            order = orders[jerp]
            nightsbox0 = np.array(boxdata['NIGHT'+str(orders[jerp])])
            parfitbox0 = np.array(boxdata['PARFIT'+str(orders[jerp])])
            IPnights1 = {}
            for night in nightsFinal:
                indnight  = np.where(nightsbox0 == night)[0]
                IPnights1[night] = np.nanmean(parfitbox0[indnight][:,5])
                IP2 = np.nanmean(parfitbox0[indnight][:,13])
                IP3 = np.nanmean(parfitbox0[indnight][:,14])


            fakeargs = [kind, 4000., T2, 0.9, pow2, 200,IPnights1,IP2,IP3]

            #BIinst( orders, jerp, nightsFinal,figpath,logger,fakeargs,rvfinal1,rvfinal2,0)
            func = partial(BIinst, orders, jerp, nightsFinal,figpath,logger,fakeargs,rvfinal1,rvfinal2)
            outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=5)

            # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
            outs = np.array(outs)
            nightsbox = outs[:,0]
            bibox = outs[:,1]
            stdbox = outs[:,2]

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
        if not debug: logger.addHandler(stream_hander)
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


        # Iterate over tight and loose mounting data sets...
        for boxind in range(len(rvboxcomblist)):

            bimasterbox  = rvboxcomblist[boxind]
            stdmasterbox = stdboxcomblist[boxind]

            #-------------------------------------------------------------------------------

            if mode=='STD': # If RV STD star, calculate uncertainty in method

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
                    sigma_method2 = inparam.BImethodvariance_tight[band]
                elif T_Ls[boxind] == 'L':
                    nights_use = nightsL.copy()
                    kind = 'Defocus'
                    sigma_method2 = inparam.BImethodvariance_loose[band]

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
            axes.set_ylim(np.nanmin(bifinal)-.08, np.nanmax(bifinal)+.08)
            #axes.set_ylim(-0.3,0.3)
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
                rvfinal[aa] =  rvfinal1[ind[0]]
                rvstdfinal[aa] = rvstdfinal1[ind[0]]
                jdfinal[aa] = jd0[ind[0]]


            f, axes = plt.subplots(1, 1, figsize=(8,6), facecolor='white', dpi=300)

            mask = np.ones_like(bifinal,dtype=bool)
            mask[np.isnan(bifinal) | np.isnan(rvfinal)] = False
            pp = pearsonr(rvfinal[mask],bifinal[mask])

            axes.plot(    rvfinal, bifinal, '.k', ms=5)
            axes.errorbar(rvfinal, bifinal, xerr=rvstdfinal,yerr=stdfinal, ls='none', lw=.5, ecolor='black')
            axes.text(0.05, 0.93, 'Correlation Coeff = {} , P = {}'.format(round(pp[0],5),round(pp[1],5)),
                                 transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.set_ylim(np.nanmin(bifinal)-.08,
                         np.nanmax(bifinal)+.08)
            #axes.set_ylim(-0.3,0.3)
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

            if mode=='STD':
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

            if mode=='STD': # If uncertainty in method was calculated, save it
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
        axes.set_ylim(-0.1,0.25)
        #axes.set_ylim(-0.3,0.3)
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

        axes.plot(    rvfinal1, bifinal1, '.k', ms=5,alpha=0.5,label='Data')
        axes.errorbar(rvfinal1, bifinal1, xerr=rvstdfinal1,yerr=bistdfinal1, ls='none', lw=.5, ecolor='black',alpha=0.5)

        axes.plot(    rvfinalCombined, bifinalCombined, '.r', ms=5,label='Model')
        axes.errorbar(rvfinalCombined, bifinalCombined, xerr=rvstdfinalCombined,yerr=stdfinalCombined, ls='none', lw=.5, ecolor='red')
        axes.text(0.05, 0.93, 'Correlation Coeff = {} , P = {}'.format(round(pp[0],5),round(pp[1],5)),
                             transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        axes.legend()
        axes.set_ylim(-0.1,0.25)
        axes.set_xlim(12.5,15.75)
        axes.set_xlabel('RV [km/s]', size=6, style='normal', family='sans-serif' )
        axes.set_ylabel('Bisector Span [km/s]', size=6, style='normal', family='sans-serif' )
        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(5))
        axes.tick_params(axis='both', which='both', labelsize=5, right=True, top=True, direction='in', width=.6)
        f.savefig('{}/FinalRVs_vs_BIs.png'.format(outpath), format='png', bbox_inches='tight')

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


    print('\n')
    logger.info('The final BI estimates you are looking for are in the BIresultsSummary files!')
    print('####################################################################################')
