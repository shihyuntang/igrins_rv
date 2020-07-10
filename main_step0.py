from Engine.importmodule import *
# -------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def DataPrep(args):

    # Collects and organizes all relevant information on target observations and associated telluric standard observations, for ease of use in later steps.
    # Requires that observation be listed in IGRINS_RV_MASTERLOG.csv that comes with this package in the /Engine folder.
    # If your target is not listed, you must construct your own PrepData files. See ReadMe for more details.

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
    # Prepare PrepData file for target star
    fileT = open('./Input/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), 'w')
    fileT.write('night beam tag mjd facility airmass bvc\n')

    # Collect target star information
    nightsT = [];
    for x in range(len(star_files['CIVIL'])):
        night    = str(  np.array(star_files['CIVIL'])[x]     )
        frame    = str(  np.array(star_files['FRAMETYPE'])[x] )
        tag0     = int(  np.array(star_files['FILENUMBER'])[x])
        airmass  = float(np.array(star_files['AM'])[x]        )
        # BVCfile  = float(np.array(star_files['BVC'])[x]       )
        facility = str(  np.array(star_files['FACILITY'])[x]  )
        tag = '{:04d}'.format(tag0)

        # If observation in MASTER_LOG not found in local filesystem, don't list
        try:
            hdulist = fits.open('{}{}/{}/SDC{}_{}_{}.spec.fits'.format(inpath, night, frame, args.band, night, tag))
        except FileNotFoundError:
            print('     --> {}{}/{}/SDC{}_{}_{}.spec.fits NOT FOUND'.format(inpath, night, frame, args.band, night, tag))
            continue

        # Collect observatory
        head = hdulist[0].header
        if head['OBSERVAT'].lower() == 'lowell observatory':
            obs = 'DCT'
        elif (head['OBSERVAT'].lower() == 'mcdonald observatory') or (head['OBSERVAT'].lower()  == 'mcdonald'):
            obs = 'McD'
        else:
            print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}. CODE MUST BE EDITED TO INCLUDE THIS OPTION - CONTACT AUTHORS WITH EXAMPLE OBSERVATION FITS FILE AND THEY WILL UPDATE'.format( head['OBSERVAT'].lower() ))

        # Collect time of mid-exposure
        try:
            time_midpoint = np.mean([float(head['JD-OBS']),float(head['JD-END'])])
        except KeyError:
            l0 = []
            for nm in ['DATE-OBS','DATE-END']:
                tt1 = head[nm].split('-')
                t1 = Time(tt1[0]+'-'+tt1[1]+'-'+tt1[2]+' '+tt1[3],format='iso')
                l0.append(t1.jd)
            time_midpoint = np.mean(l0)

        # BVC calculation
        if args.coord != '':
            print('Calculating BVC base on the input info. ...')
            print(args.coord)
            ra_deg = np.array(ast.literal_eval(args.coord), dtype=float)[0]
            de_deg = np.array(ast.literal_eval(args.coord), dtype=float)[1]

            pmra_deg = np.array(ast.literal_eval(args.pm), dtype=float)[0]
            pmde_deg = np.array(ast.literal_eval(args.pm), dtype=float)[1]

            targ_c = SkyCoord(ra  =  ra_deg                   *units.degree,
                              dec =  de_deg                   *units.degree,
                              pm_ra_cosdec = pmra_deg         *units.mas/units.yr,
                              pm_dec       = pmde_deg         *units.mas/units.yr,
                              distance = float(args.distance) *units.pc,
                              frame='icrs',
                              obstime="J2015.5")

            new_coord = targ_c.apply_space_motion(new_obstime=Time(time_midpoint, format='jd'))


            observatoryN = EarthLocation.of_site('McDonald Observatory')
            new_RA = new_coord.ra
            new_DE = new_coord.dec

            sc = SkyCoord(ra=new_RA, dec=new_DE, frame=ICRS)

            barycorr  = sc.radial_velocity_correction(obstime=Time(time_midpoint, format='jd'), location=observatoryN)
            BVCfile   = barycorr.to(units.km/units.s).value

        else:
            if obs == 'McD':
                print('BVC stright from master log ...')
                observatoryN = EarthLocation.of_site('McDonald Observatory')
                BVCfile  = float(np.array(star_files['BVC'])[x]       ) #BVC in the master log might be wrong, so, re-calculated below...

            elif obs == 'DCT':
                print('Calculating BVC base on the fits header info. ...')
                observatoryN = EarthLocation.of_site('DCT')

                framee = f"{head['RADECSYS'][:2].lower()}{head['RADECSYS'][-1]}"
                sc = SkyCoord(f"{head['TELRA']} {head['TELDEC']}", frame=framee, unit=(units.hourangle, units.deg))
                barycorr = sc.radial_velocity_correction(obstime=Time(time_midpoint, format='jd'), location=observatoryN)
                BVCfile = barycorr.to(units.km/units.s).value

        # Write out collected info
        mjd = time_midpoint;
        fileT.write(night+' '+frame+' '+str(tag)+' '+str(mjd)+' '+str(facility)+' '+str(airmass)+' '+str(BVCfile))
        fileT.write('\n')
        nightsT.append(night)
    fileT.close()

    #-------------------------------------------------------------------------------
    # Prepare PrepData file for A0 stars
    fileA0 = open('./Input/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), 'w')
    fileA0.write('night tag humid temp zd press obs airmass\n')
    noA0nights = []

    ## Now collect A0 information
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
            # Then check whether any A0s files for that night are present in local filesystem.
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
                # If so, review all STD observations for the night and use the one taken soonest before or after the target observation that is also within the set airmass range.
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

                # If best STD observation still has airmass above range limit, throw warning
                if abs(am0-am_star) > float(args.AM_cut):
                    print(night,stdname,am_star,am0,tagA0)
                    sys.exit('WARNING, STD (A0) AIRMASS FOR NIGHT {} HAS A DIFFERENCE LARGER THAN {} FROM TARGET!'.format(night, args.AM_cut))

                tagA = '{:04d}'.format(tagA0)
                subpath = '{}std/{}/AB/SDC{}_{}_{}.spec.fits'.format(inpath, night, args.band, night, tagA)

                # Open the chosen STD observation
                try:
                    hdulist = fits.open(subpath)
                    head    = hdulist[0].header

                # If it is not found locally, use whatever A0 is present locally, in case user selected their STD observations differently
                except FileNotFoundError:
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

                # Collect observatory
                if head['OBSERVAT'] == 'Lowell Observatory':
                    obs = 'DCT'
                elif (head['OBSERVAT'] == 'McDonald Observatory') or (head['OBSERVAT'] == 'McDonald'):
                    obs = 'McD'
                else:
                    print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}. CODE MUST BE EDITED TO INCLUDE THIS OPTION - CONTACT AUTHORS WITH EXAMPLE OBSERVATION FITS FILE AND THEY WILL UPDATE'.format( head['OBSERVAT'] ))

                # Collect other relevant information for use with Telfit
                AM = str(am0)
                try:
                    humid = float(head['HUMIDITY'])
                    temp  = float(head['AIRTEMP'])
                    press = float(head['BARPRESS'])
                    zd = np.mean([float(head['ZDSTART']),float(head['ZDEND'])])
                except ValueError: # McDonald headers don't have these quantities :( They will be fit by Telfit instead of set to their recorded values.
                    humid = 'NOINFO'; temp = 'NOINFO'; press = 'NOINFO'; zd = 'NOINFO';

        # Write to file
        fileA0.write(night+' '+str(tagA)+' '+str(humid)+' '+str(temp)+' '+str(zd)+' '+str(press)+' '+str(obs)+' '+str(AM))
        fileA0.write('\n')

    fileA0.close()

    if len(noA0nights) != 0:
        print('No reduced A0s found for following nights:')
        print(noA0nights)
        print('To achieve highest precision, this pipeline defaults to not analyzing target spectra for these nights.\n')

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 0',
                                     description = '''
                                     This step collects and organizes all relevant information on target observations and associated telluric standard observations, for ease of use in later steps.
                                     It requires that your observations be listed in IGRINS_RV_MASTERLOG.csv, which comes with this package in the /Engine folder.
                                     If your target is not listed, you must construct your own PrepData files. See ReadMe for more details. \n
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your target name, no space",
                        type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-AM",      dest="AM_cut",           action="store",
                        help="AirMass difference allowed between TAR and STD (A0) stars. Default X = 0.25 ",
                        type=str,   default='0.25')

    parser.add_argument("-coord",    dest="coord",            action="store",
                        help="Optional [-XX.xx,-XX.xx] deg, GaiaDR2 coordinates at J2015.5. If give, will calculate BVC base on this info.",
                        type=str,   default='')
    parser.add_argument("-pm",       dest="pm",               action="store",
                        help="Optional [-XX.xx,-XX.xx] [mas/yr], GaiaDR2 proper motion. If give, will calculate BVC base on this info.",
                        type=str,   default='')
    parser.add_argument("-dist",    dest="distance",          action="store",
                        help="Optional (pc), can be from GaiaDR2 parallax [mas] (1/plx), or from Bailer-Jones et al. 2018. If give, will calculate BVC base on this info.",
                        type=str,   default='')

    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.85')
    args   = parser.parse_args()
    inpath = './Input/{}/'.format(args.targname)

    if ' ' in args.targname:
        sys.exit('ERROR! This is weird... you have a SPACE in your input *target name...')

    if not os.path.isdir(f'./Input/Prepdata'):
        os.mkdir(f'./Input/Prepdata')

#-------------------------------------------------------------------------------
    print('####################################################################################\n')
    print(f'Making your *Prepdata* for {args.targname}')
    time.sleep(1)

    DataPrep(args)

    print('Finished!')
    print('Prepdata saved under ./Input/Prepdata/\n')
    print('####################################################################################')
