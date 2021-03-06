import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *
from Engine.IO_AB    import setup_templates, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips    import basicclip_above
from Engine.contfit  import A0cont
from Engine.classes  import fitobjs, inparamsA0
from Engine.rebin_jv import rebin_jv
from Engine.rotint   import rotint
from Engine.opt      import optimizer
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
def DataPrep(args, tar_night, tar_num, tar_frame, file_night_num, std_name, std_num, std_night):
    star = args.targname.replace(' ', '')
    inpath     = '../Input/{}/'.format(star)

    # Find all nights of observations of target in master log
    master_log_fh = '../Engine/IGRINS_MASTERLOG.csv'
    master_log = Table.from_pandas(pd.read_csv(master_log_fh))

    master_log['CIVIL'].format = '%i'
    master_log['FILENUMBER'].format = '%i'

    # Collect target star information
    fileT = open('../Input/Prepdata/Prepdata_targ_{}.txt'.format(star), 'w')
    fileT.write('night beam tag mjd facility airmass bvc\n')

    nightsT = []
#    for x in range(len(star_files['CIVIL'])):
    print('target nights', tar_night)
    for x in range(len(tar_night)):
        night    = str(tar_night[x])
        tag0     = int(tar_num[x])
        frame    = str(tar_frame[x])
        file_nn = file_night_num[x]

        star_files = master_log[ ( master_log['CIVIL']==int(tar_night[x]) ) &
                                 ( master_log['FILENUMBER']==int(tar_num[x]) )  ]

        airmass  = float(np.array(star_files['AM'])[0])
        # BVCfile  = float(np.array(star_files['BVC'])[0])
        facility = str(np.array(star_files['FACILITY'])[0])

        tag = '{:04d}'.format(tag0)
        try:
            temp_dir = 'SDC{}_{}_{}.spec.fits'.format(args.band, night, tag)
            if temp_dir in os.listdir('{}{}_{}/A'.format(inpath, night, tag)):
                print('{}{}_{}/A/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, args.band, night, tag))
                hdulist = fits.open('{}{}_{}/A/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, args.band, night, tag))
                save_yn = 0 ; fram_ty = 'A'
                print('first')
            else:
                print('{}{}_{}/B/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, args.band, night, tag))
                hdulist = fits.open('{}{}_{}/B/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, args.band, night, tag))
                save_yn = 0 ; fram_ty = 'B'
                print('first')

            # hdulist = fits.open('{}{}_{}/{}/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, frame, args.band, night, tag))
            # print('{}{}_{}/{}/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag, frame, args.band, night, tag))
            # print('first')
            # save_yn = 0

        except:
            temp_dir = 'SDC{}_{}_{}.spec.fits'.format(args.band, night, tag)
            if temp_dir in os.listdir('{}{}_{}/A'.format(inpath, night, tag_temp)):
                print('{}{}_{}/A/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag_temp, args.band, night, tag))
                hdulist = fits.open('{}{}_{}/A/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag_temp, args.band, night, tag))
                save_yn = 1 ; fram_ty = 'A'
            else:
                print('{}{}_{}/B/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag_temp, args.band, night, tag))
                hdulist = fits.open('{}{}_{}/B/SDC{}_{}_{}.spec.fits'.format(inpath, night, tag_temp, args.band, night, tag))
                save_yn = 1 ; fram_ty = 'B'


        if save_yn == 0:
            tag_temp = tag

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
                observatoryN = EarthLocation.of_site('McDonald Observatory')
                BVCfile  = float(np.array(star_files['BVC'])[x]       ) #BVC in the master log might be wrong, so, re-calculated below...

            elif obs == 'DCT':
                print('Calculating BVC base on the fits header info. ...')
                observatoryN = EarthLocation.of_site('DCT')

                framee = f"{head['RADECSYS'][:2].lower()}{head['RADECSYS'][-1]}"
                sc = SkyCoord(f"{head['TELRA']} {head['TELDEC']}", frame=framee, unit=(units.hourangle, units.deg))
                barycorr = sc.radial_velocity_correction(obstime=Time(time_midpoint, format='jd'), location=observatoryN)
                BVCfile = barycorr.to(units.km/units.s).value

        mjd = time_midpoint
        fileT.write(file_nn+' '+fram_ty+' '+str(tag)+' '+str(mjd)+' ' +
                    str(facility)+' '+str(airmass)+' '+str(BVCfile))
        fileT.write('\n')
        nightsT.append(night)
    fileT.close()

    # Now collect A0 information
    fileA0 = open('../Input/Prepdata/Prepdata_A0_{}.txt'.format(star), 'w')
    fileA0.write('night tag humid temp zd press obs airmass\n')
    noA0nights = []

    for x in range(len(std_night)):
        night = str(std_night[x])
        name  = str(std_name[x])
        num   = str(std_num[x])

        anyK = False
        subpath        = '{}std/{}/A/'.format(inpath, night)
        fullpathprefix = '{}SDC{}_{}_'.format(subpath, args.band, night)

        onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
        for f in onlyfiles:
            q = re.split('_', f)
            if q[0] != 'SDC{}'.format(args.band):
                continue
            anyK = True

        if anyK == False:
            tagA = 'NA'
            humid = 'NA'
            temp = 'NA'
            zd = 'NA'
            press = 'NA'
            obs = 'NA'
            AM = 'NA'
            noA0nights.append(night)
        else:
            stdname = name
            tagA0 = int(num)

            tagA = '{:04d}'.format(tagA0)
            subpath = '{}std/{}/A/SDC{}_{}_{}.spec.fits'.format(inpath, night, args.band, night, tagA)

            try:
                hdulist = fits.open(subpath)
                head = hdulist[0].header

            except FileNotFoundError:
                # If best airmass match A0 for night not found, check if Joe chose a different A0 instead
                subpath        = '{}std/{}/A/'.format(inpath, night)
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
                    am0 = np.mean([float(head['AMSTART']), float(head['AMEND'])])

            if head['OBSERVAT'] == 'Lowell Observatory':
                obs = 'DCT'
            elif (head['OBSERVAT'] == 'McDonald Observatory') or (head['OBSERVAT'] == 'McDonald'):
                obs = 'McD'
            else:
                sys.exit('EXPECTED LOWELL OR McDonald, GOT ' +
                      str(head['OBSERVAT'])+', MUST EDIT CODE TO INCLUDE THIS OPTION')
            try:
                humid = float(head['HUMIDITY'])
            except:
                humid = 'NOINFO'

            try:
                temp  = float(head['AIRTEMP'])
            except:
                temp  = 'NOINFO'

            try:
                press = float(head['BARPRESS'])
            except:
                press = 'NOINFO'

            try:
                zd    = np.mean([float(head['ZDSTART']), float(head['ZDEND'])])
            except:
                zd    = 'NOINFO'

            try:
                am0   = np.mean([float(head['AMSTART']), float(head['AMEND'])])
            except:
                am0   = 'NOINFO'

            AM = str(am0)

        fileA0.write(night+' '+str(tagA)+' '+str(humid)+' '+str(temp) +
                     ' '+str(zd)+' '+str(press)+' '+str(obs)+' '+str(AM))
        fileA0.write('\n')
    fileA0.close()

    if len(noA0nights) != 0:
        print('No reduced A0s found for following nights:')
        print(n)
        print('To achieve highest precision, this pipeline defaults to not analyzing target spectra for these nights.\n')


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
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

    parser.add_argument("-coord",    dest="coord",            action="store",
                        help="Optional [-XX.xx,-XX.xx] deg, GaiaDR2 coordinates at J2015.5. If give, will calculate BVC base on this info.",
                        type=str,   default='')
    parser.add_argument("-pm",       dest="pm",               action="store",
                        help="Optional [-XX.xx,-XX.xx] [mas/yr], GaiaDR2 proper motion. If give, will calculate BVC base on this info.",
                        type=str,   default='')
    parser.add_argument("-dist",    dest="distance",          action="store",
                        help="Optional (pc), can be from GaiaDR2 parallax [mas] (1/plx), or from Bailer-Jones et al. 2018. If give, will calculate BVC base on this info.",
                        type=str,   default='')

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)" % (
                            mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2))
    parser.add_argument('--version',
                        action='version',  version='%(prog)s 0.9')
    args = parser.parse_args()
    cdbs_loc = '~/cdbs/'
    inpath     = '../Input/{}/'.format(args.targname)

    new_tar_list = os.listdir('./{}_recipes'.format(args.targname.replace(' ', '')))
    target_have  = np.sort([int(dump[:8]) for dump in new_tar_list if dump[-1] == 'p'])

    tar_night = []
    tar_num   = []
    tar_frame = []
    file_night_num = []

    std_night = []
    std_name = []
    std_num  = []
    for i in target_have:
        tempp = Table.read('./{}_recipes/{}.recipes.tmp'.format(args.targname.replace(' ', ''), i), format='ascii')

        temp = tempp[ tempp['OBJNAME'] == '{}'.format(args.targname) ]
        for rows in range(len(temp)):
            nums   = temp[rows]['OBSIDS'].split(' ')
            [tar_num.append(nn) for nn in nums]

            frames = temp[rows]['FRAMETYPES'].split(' ')
            for nn in frames:
                tar_frame.append(nn)
                tar_night.append(i)

                file_night_num.append('{}_{:04d}'.format(i, int(temp[rows]['GROUP1']) ))

        temp = tempp[ tempp['OBJTYPE'] == 'STD' ]
        for rows in range(len(temp)):
            std_name.append(temp[rows]['OBJNAME'])
            std_num.append(temp[rows]['GROUP1'])
            std_night.append(i)

    # ------------
    start_time=datetime.now()
    print('\n')
    print('###############################################################')
    print('Data Preparation for {} (1/2)...'.format(args.targname))
    time.sleep(1)
    DataPrep(args, tar_night, tar_num, tar_frame, file_night_num, std_name, std_num, std_night)

    print('Data Preparation Done!')
