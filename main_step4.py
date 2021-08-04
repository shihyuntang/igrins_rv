from Engine.importmodule import *

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 4',
                                     description = '''
                                     Updates RV uncertainty estimates to take into account uncertainty in vsini. Takes two runs of Step 3, one with vsini held fixed at the best guess value and one with vsini held fixed at the best guess value plus or minus one sigma, and uses the difference between the two to produce updated RVs and uncertainties.

                                     For the most part, the uncertainties should change little (~1 m/s), but for high vsini (>~ 15 km/s) objects, it may increase by ~5-7 m/s or so.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-run1",    dest="run1",             action="store",
                        help="First step3 run that will be used, the one with vsini held fixed AT THE BEST GUESS. Takes the string that suffixes 'RV_results_', e.g. for 'RV_results_1', you would set this to '1'.",
                        type=str,   default='')
    parser.add_argument("-run2",    dest="run2",             action="store",
                        help="Second step3 run that will be used, the one with vsini held fixed at the best guess PLUS OR MINUS SIGMA.",
                        type=str,   default='')
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 1.0.0')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    # Check user input

    if args.run1 == '' or args.run2 == '':
        sys.exit('ERROR: REQUIRES TWO RUNS SPECIFIED THROUGH -run1 AND -run2 FLAGS!')


    #-------------------------------------------------------------------------------
    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[41m {} band \033[0m
    Run 1               = \33[41m RV_results_{} \033[0m   <------- Should be run where vsini = best guess
    Run 2               = \33[41m RV_results_{} \033[0m   <------- Should be run where vsini = best guess +- sigma
    '''.format(args.targname, args.band, args.run1,args.run2))
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
    print('Running Step 4 for {}...'.format(args.targname))

    #-------------------------------------------------------------------------------

    name = f'RV_results_{args.run1}_{args.run2}_combined'

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/{name}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/{name}')
    else:
        sys.exit(f'OOPS! \"./Output/{args.targname}_{args.band}/{name}\" ALREADY EXIST!!')

    outpath = f'./Output/{args.targname}_{args.band}'
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

    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{name}')
    #-------------------------------------------------------------------------------

    rvses = []

    for run in [args.run2, args.run1]:
        hdu    = fits.open('{}/{}/RVresultsSummary.fits'.format(outpath, f'RV_results_{run}'))
        tbdata = hdu[1].data
        rvses.append(np.array(tbdata['RVfinal'],dtype=float))
        stds               = np.array(tbdata['STDfinal'], dtype=float)
        jdsCombined        = np.array(tbdata['JD'],       dtype=float)
        nightsCombined     = np.array(tbdata['NIGHT'],    dtype=str)
        vsinifinalCombined = np.array(tbdata['VSINI'],    dtype=float)

    vsini_err        = np.nanstd(rvses[0]-rvses[1])
    stdfinalCombined = np.sqrt(stds**2 + vsini_err**2)
    rvfinalCombined  = rvses[1]

    intnights = np.array([int(i[:8]) for i in nightsCombined])
    indT = np.where((intnights <  20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))
    nightsT = nightsCombined[indT]
    nightsL = nightsCombined[indL]
    #-------------------------------------------------------------------------------

    # Plot combined results
    xscale = np.arange(len(rvfinalCombined))+1

    f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(xscale,rvfinalCombined, '.k', ms=5)
    axes.errorbar(xscale,rvfinalCombined,yerr=stdfinalCombined,ls='none',lw=.5, ecolor='black')
    axes.text(0.05, 0.93, r'RV mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=axes.transAxes, size=6, style='normal', family='sans-serif' )

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        axes.text(0.05, 0.1, 'Defocus', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            axes.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        else:
            axes.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    axes.set_ylim(np.nanmin(rvfinalCombined)-.08,np.nanmax(rvfinalCombined)+.08)
    axes.set_ylabel('RV (km/s)', size=6, style='normal', family='sans-serif' )
    axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/{}/FinalRVs.png'.format(outpath, name), format='png', bbox_inches='tight')
    #-------------------------------------------------------------------------------

    # Output combined final results to fits file
    c1 = fits.Column(name='NIGHT',    array=nightsCombined,         format='{}A'.format(len(nightsCombined[0])) )
    c2 = fits.Column(name='JD',       array=jdsCombined,            format='D')
    c3 = fits.Column(name='RVfinal',  array=rvfinalCombined,        format='D')
    c4 = fits.Column(name='STDfinal', array=stdfinalCombined,       format='D')
    c5 = fits.Column(name='VSINI',    array=vsinifinalCombined,     format='D')

    cols = fits.ColDefs([c1,c2,c3,c4,c5])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    bleh = np.ones((3,3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu,hdu_1])
    hdul.writeto('{}/{}/RVresultsSummary.fits'.format(outpath, name),overwrite=True)

    tempin = Table.read('{}/{}/RVresultsSummary.fits'.format(outpath, name), format='fits')
    tempin.write('{}/{}/RVresultsSummary.csv'.format(outpath, name), format='csv', overwrite=True)
    #-------------------------------------------------------------------------------

    logger.info('Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(rvfinalCombined),
                                                                                  np.nanstd(rvfinalCombined)))
    print('\n')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}/{}'.format(outpath, name) )
    print('####################################################################################')
