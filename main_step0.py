from Engine.importmodule import *

#-------------------------------------------------------------------------------
def DataPrep(args):
# Find all nights of observations of target in master log
    master_log    = pd.read_csv('./Engine/IGRINS_MASTERLOG.csv')
    star_files    = master_log[(master_log['OBJNAME'].str.contains(args.targname, regex=True, na=False)) &
                               (master_log['OBJTYPE'].str.contains('TAR',         regex=True, na=False)) ]
    allnights     = np.array(master_log['CIVIL'],dtype='str')

    # If star input not found in Masterlog, try putting a space in its name somewhere
    n = 1
    while len(star_files['CIVIL']) == 0:
        starnew = args.targname[:n]+' '+args.targname[n:]
        star_files = master_log[(master_log['OBJNAME'].str.contains(starnew, regex=True, na=False)) &
                                (master_log['OBJTYPE'].str.contains('TAR',   regex=True, na=False)) ]
        n += 1
        if n == len(args.targname):
            sys.exit('TARGET NAME NOT FOUND IN CATALOG - CHECK INPUT!')

#-------------------------------------------------------------------------------
# Collect target star information
    fileT = open('./Temp/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), 'w')
    fileT.write('night beam tag mjd facility airmass bvc\n')

    nightsT = [];
    for x in range(len(star_files['CIVIL'])):
        night    = str(  np.array(star_files['CIVIL'])[x]     )
        frame    = str(  np.array(star_files['FRAMETYPE'])[x] )
        tag0     = int(  np.array(star_files['FILENUMBER'])[x])
        airmass  = float(np.array(star_files['AM'])[x]        )
        BVCfile  = float(np.array(star_files['BVC'])[x]       )
        facility = str(  np.array(star_files['FACILITY'])[x]  )
        tag = '{:04d}'.format(tag0)

        try:
            hdulist = fits.open('{}{}/{}/SDC{}_{}_{}.spec.fits'.format(inpath, night, frame, args.band, night, tag))
        except FileNotFoundError:
            continue

        head = hdulist[0].header
        if head['OBSERVAT'] == 'Lowell Observatory':
            obs = 'DCT'
        elif head['OBSERVAT'] == 'McDonald':
            obs = 'McD'
        else:
            print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}, MUST EDIT CODE TO INCLUDE THIS OPTION'.format( head['OBSERVAT'] ))

        try:
            time_midpoint = np.mean([float(head['JD-OBS']),float(head['JD-END'])])
        except KeyError:
            l0 = []
            for nm in ['DATE-OBS','DATE-END']:
                tt1 = head[nm].split('-')
                t1 = Time(tt1[0]+'-'+tt1[1]+'-'+tt1[2]+' '+tt1[3],format='iso')
                l0.append(t1.jd)
            time_midpoint = np.mean(l0)

        mjd = time_midpoint;
        fileT.write(night+' '+frame+' '+str(tag)+' '+str(mjd)+' '+str(facility)+' '+str(airmass)+' '+str(BVCfile))
        fileT.write('\n')
        nightsT.append(night)
    fileT.close()
#-------------------------------------------------------------------------------
    ## Now collect A0 information
    fileA0 = open('./Temp/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), 'w')
    fileA0.write('night tag humid temp zd press obs airmass\n')
    noA0nights = []

    for night0 in np.unique(nightsT):
        night    = str(night0)
        am_stars = star_files['AM'][(np.array(star_files['CIVIL'], dtype=str) == night)]
        am_star  = float(am_stars.values[0])

        std_files = master_log[(allnights == night) & (master_log['OBJTYPE'].str.contains('STD', regex=True, na=False))]
        ams = np.array(std_files['AM'].values,dtype='float')

        # First check whether any A0s observed that night
        if len(ams) == 0:
            noA0nights.append(night)
            tagA = 'NA'; humid = 'NA'; temp = 'NA'; zd = 'NA'; press = 'NA'; obs = 'NA'; AM = 'NA';
        else:
            # Then check whether any A0s files for that night outputted by reduction pipeline.
            # If not, Joe either didn't have the data for them or didn't copy them over.
            anyK = False
            subpath        = '{}std/{}/AB/'.format(inpath, night)
            fullpathprefix = '{}SDC{}_{}_'.format(subpath, args.band, night)

            onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
            for f in onlyfiles:
                q = re.split('_', f)
                if q[0] != 'SDC{}'.format(args.band):
                    continue
                anyK = True

            if anyK == False:
                tagA = 'NA'; humid = 'NA'; temp = 'NA'; zd = 'NA'; press = 'NA'; obs = 'NA'; AM = 'NA';
                noA0nights.append(night)
            else:
                stdname = std_files['OBJNAME'][abs(ams-am_star) == min(abs(ams-am_star))].values[0]
                fileno0 = int(std_files['FILENUMBER'][abs(ams-am_star) == min(abs(ams-am_star))].values[0])
                names   = np.array(std_files['OBJNAME'])
                tagA0s  = np.array(std_files['FILENUMBER'][(names == stdname)].values)
                facs    = np.array(std_files['FACILITY'][(names == stdname)].values)
                am0s    = ams[(names == stdname)]

                firsts = []
                for k, g in groupby(enumerate(tagA0s), lambda ix : ix[0] - ix[1]):
                    q = list(map(itemgetter(1), g))
                    firsts.append(q[0])
                firsts = np.array(firsts)

                tagA0 = firsts[abs(firsts-fileno0) == min(abs(firsts-fileno0))][0]
                am0   = am0s[(tagA0s == tagA0)][0]
                fac   = np.array(facs)[(tagA0s == tagA0)][0]

                if abs(am0-am_star) > float(args.AM_cut):
                    print(night,stdname,am_star,am0,tagA0)
                    sys.exit('WARNING, STD (A0) AIRMASS FOR NIGHT {} HAS A DIFFERENCE LARGER THAN {} FROM TARGET!'.format(night, args.AM_cut))

                tagA = '{:04d}'.format(tagA0)
                subpath = '{}std/{}/AB/SDC{}_{}_{}.spec.fits'.format(inpath, night, args.band, night, tagA)

                try:
                    hdulist = fits.open(subpath)
                    head    = hdulist[0].header

                except FileNotFoundError:
                    # If best airmass match A0 for night not found, check if Joe chose a different A0 instead
                    subpath        = '{}std/{}/AB/'.format(inpath, night)
                    fullpathprefix = '{}SDC{}_{}_'.format(subpath, args.band, night)

                    onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
                    for f in onlyfiles:
                        q = re.split('_', f)
                        if q[0] != 'SDC{}'.format(args.band):
                            continue
                        qr = re.split('\.', q[2])
                        tagA = qr[0]

                        hdulist = fits.open(fullpathprefix+tagA+'.spec.fits')
                        head = hdulist[0].header
                        am0 = np.mean([float(head['AMSTART']),float(head['AMEND'])])

                if head['OBSERVAT'] == 'Lowell Observatory':
                    obs = 'DCT'
                elif (head['OBSERVAT'] == 'McDonald Observatory') or (head['OBSERVAT'] == 'McDonald'):
                    obs = 'McD'
                else:
                    print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}, MUST EDIT CODE TO INCLUDE THIS OPTION'.format( head['OBSERVAT'] ))

                AM = str(am0)
                try:
                    humid = float(head['HUMIDITY'])
                    temp  = float(head['AIRTEMP'])
                    press = float(head['BARPRESS'])
                    zd = np.mean([float(head['ZDSTART']),float(head['ZDEND'])])
                except ValueError: # McDonald headers don't have these quantities :(
                    humid = 'NOINFO'; temp = 'NOINFO'; press = 'NOINFO'; zd = 'NOINFO';

        fileA0.write(night+' '+str(tagA)+' '+str(humid)+' '+str(temp)+' '+str(zd)+' '+str(press)+' '+str(obs)+' '+str(AM))
        fileA0.write('\n')

    fileA0.close()

    print('No reduced A0s found for following nights:')
    for n in noA0nights:
        print(n)
    print('To achieve highest precision, this pipeline defaults to not analyzing target spectra for these nights.\n')





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
                            help="Enter your *target name, no space",
                            type=str)
        parser.add_argument("-HorK",    dest="band",             action="store",
                            help="Which band to process? H or K?. Default = K",
                            type=str,   default='K')

        parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                            help="If sets, will generate files and plots under ./Temp/Debug for debug")
        parser.add_argument('--version',                         action='version',  version='%(prog)s 0.85')
        args   = parser.parse_args()
        inpath = './Input/{}/'.format(args.targname)

        if ' ' in args.targname:
            sys.exit('ERROR! This is Weird... you have a SPACE in your input target name...')
        # make logging dir
        if not os.path.isdir(f'./Runlog/{args.targname}_{args.band}'):
            os.mkdir(f'./Runlog/{args.targname}_{args.band}')

        # make logging dir
        if not os.path.isdir(f'./Input/Prepdata'):
            os.mkdir(f'./Input/Prepdata')

    #-------------------------------------------------------------------------------
        print('\n')
        print('###############################################################')
        print(f'Making your Prepdata for {args.targname}')
        time.sleep(3)

        DataPrep(args)

        print('Finished!')
        print('Prepdata saved under ./Input/Prepdata/')
        print('###############################################################')
        print('\n')
