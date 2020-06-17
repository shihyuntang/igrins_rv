import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *
from Engine.rebin_jv import rebin_jv
#-------------------------------------------------------------------------------


def IPval(tar,band):

    filesndirs = os.listdir('./A0_Fits/A0_Fits_{}_IP/'.format(tar))
    filesndirs = [j for j in filesndirs if j[-6:] == '{}.fits'.format(band)]

    nights  = np.array([int(j[:8]) for j in filesndirs ])
    nightsT = np.where((nights < 20180401)  | (nights > 20190531))
    nightsL = np.where((nights >= 20180401) & (nights < 20190531))

    Tdirs = [ './A0_Fits/A0_Fits_{}_IP/{}A0_treated_{}.fits'.format(tar, nn, band) for nn in nights[nightsT] ]
    Ldirs = [ './A0_Fits/A0_Fits_{}_IP/{}A0_treated_{}.fits'.format(tar, nn, band) for nn in nights[nightsL] ]

    filew = open('../Temp/IP_{}.txt'.format(band),'w')

    if len(nightsT) != 0:
        dump1 = 0
        for a0 in Tdirs:
            hdulist = fits.open(a0)

            tt= 1 ; orders = []
            while 1==1:
                try:
                    orders.append( int(hdulist[tt].columns[0].name[9:]) )
                    tt+=1
                except:
                    break

            dump2 = 0
            for o in np.arange(len(orders)):
                try:
                    tbdata = hdulist[o+1].data
                    if dump2 == 0:
                        IP14 = tbdata['PARFIT'][14]
                        IP13 = tbdata['PARFIT'][13]
                        IP5  = tbdata['PARFIT'][5]
                        dump2 += 1
                    else:
                        IP14 = np.append(IP14, tbdata['PARFIT'][14])
                        IP13 = np.append(IP13, tbdata['PARFIT'][13])
                        IP5  = np.append(IP5,  tbdata['PARFIT'][5])
                except:
                        IP14 = np.append(IP14, np.nan)
                        IP13 = np.append(IP13, np.nan)
                        IP5  = np.append(IP5,  np.nan)
            if dump1 == 0:
                IP14box = IP14
                IP13box = IP13
                IP5box  = IP5
                dump1 += 1
            else:
                IP14box = np.vstack((IP14box, IP14))
                IP13box = np.vstack((IP13box, IP13))
                IP5box  = np.vstack((IP5box,  IP5))

        filew.write('Tight \n')
        for o in np.arange(len(orders)):
            filew.write('{}: np.array([{:+1.8f}, {:+1.8f}, {:1.8}]),\n'.format(orders[o], np.nanmean(IP14box[:, o]), np.nanmean(IP13box[:, o]), np.nanmean(IP5box[:, o]) ))

    if len(nightsL) != 0:
        dump1 = 0
        for a0 in Ldirs:
            hdulist = fits.open(a0)

            tt= 1 ; orders = []
            while 1==1:
                try:
                    orders.append( int(hdulist[tt].columns[0].name[9:]) )
                    tt+=1
                except:
                    break

            dump2 = 0
            for o in np.arange(len(orders)):
                try:
                    tbdata = hdulist[o+1].data
                    if dump2 == 0:
                        IP14 = tbdata['PARFIT'][14]
                        IP13 = tbdata['PARFIT'][13]
                        IP5  = tbdata['PARFIT'][5]
                        dump2 += 1
                    else:
                        IP14 = np.append(IP14, tbdata['PARFIT'][14])
                        IP13 = np.append(IP13, tbdata['PARFIT'][13])
                        IP5  = np.append(IP5,  tbdata['PARFIT'][5])
                except:
                        IP14 = np.append(IP14, np.nan)
                        IP13 = np.append(IP13, np.nan)
                        IP5  = np.append(IP5,  np.nan)
            if dump1 == 0:
                IP14box = IP14
                IP13box = IP13
                IP5box  = IP5
                dump1 += 1
            else:
                IP14box = np.vstack((IP14box, IP14))
                IP13box = np.vstack((IP13box, IP13))
                IP5box  = np.vstack((IP5box,  IP5))

        filew.write('Loose \n')
        for o in np.arange(len(orders)):
            filew.write('{}: np.array([{:+1.8f}, {:+1.8f}, {:1.8}]),\n'.format(orders[o], np.nanmean(IP14box[:, o]), np.nanmean(IP13box[:, o]), np.nanmean(IP5box[:, o]) ))

    filew.close()

#-------------------------------------------------------------------------------

def WVsol(band):
    # Enter here the name of all A0_Fits directories you'd like to be used (multiple targets would be nice, for
    # consistency, but only if they use the same orders.
    sourcelist = ['A0_Fits_CITau_TEST']

    basedir = os.getcwd()
    filesndirs = os.listdir(basedir)
    first = True; orders = [];

    xmaster = np.arange(2048)
    wmaster = {}

    for f in sourcelist:
        os.chdir(f)
        for a0 in glob.glob('*.fits'):
            hdulist = fits.open(a0)

            if first == True:
                for o in np.arange(1,24):
                    try:
                        tbdata = hdulist[o]
                    except IndexError:
                        break
                    strhead = tbdata.columns[0].name
                    orders.append(int(strhead[9:]))

            for o in range(len(orders)):
                tbdata = hdulist[o+1].data
                try:
                    flag = np.array(tbdata['ERRORFLAG{}'.format(orders[o])])[0]
                except KeyError:
                    print('Warning, file {} does not have the same set of orders as others!')
                    continue
                if flag  == 1:
                    continue
                x = np.array(tbdata['X{}'.format(orders[o])])
                w = np.array(tbdata['WAVE{}'.format(orders[o])])
                x = x[(w != 0)]
                w = w[(w != 0)]
                w1 = rebin_jv(x,w,xmaster,False)

                if first == True:
                    wmaster[o] = w1
                else:
                    wmaster[o] = np.vstack((wmaster[o],w1))

            first = False

        os.chdir(basedir)


    filew = open('./WaveSolns.csv','w')

    for o in range(len(orders)):

        wmaster[o] = [np.mean(wmaster[o][:,i]) for i in range(len(wmaster[o][0,:]))]
        wmaster[o] = wmaster[o][5:-5]

        test = np.array(wmaster[o])
        diff = test[1:]-test[:-1]
        if len(diff[(diff < 0)]) > 1:
            sys.exit('ERROR! COMBINED WAVELENGTH SOLUTIONS NOT MONOTONIC! QUITTING!')

        if o != len(orders)-1:
            filew.write('x{},w{},'.format(orders[o],orders[o]))
        else:
            filew.write('x{},w{}'.format(orders[o],orders[o]))

    filew.write('\n')
    xmaster = xmaster[5:-5]

    for i in range(len(xmaster)):
        for o in range(len(orders)):
            if o != len(orders)-1:
                filew.write('{},{},'.format(xmaster[i],wmaster[o][i]/1e4))
            else:
                filew.write('{},{}'.format(xmaster[i],wmaster[o][i]/1e4))
                filew.write('\n')
    filew.close()



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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
                        help="Enter target names you wish to use as a list, e.g., [GJ281,CITau]",
                        type=str)
    parser.add_argument("-mode",    dest="mode",             action="store",
                        help="Which mode? 1) get IP & WaveSol, 2) only IP, 3) only WaveSol. Default = 1",
                        type=int,   default=int(1))
    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()

#-------------------------------------------------------------------------------
    print('\n')
    print('###############################################################')
    try:
        tars = args.targname.strip('][').split(',')
    except:
        sys.exit('NO SPACE IS ALLOWED BETWEEN NAMES!' )
#-------------------------------------------------------------------------------
    for i in tars:
        if os.path.isdir( f'./A0_Fits/A0_Fits_{i}_IP' ):
            filesndirs = os.listdir(f'./A0_Fits/A0_Fits_{i}_IP/')
            filesndirs_H = [j for j in filesndirs if j[-6:] == 'H.fits']
            filesndirs_K = [j for j in filesndirs if j[-6:] == 'K.fits']

            print('CONFIRMING... ')
            print('{} of H band & {} of K band under ./A0_Fits/A0_Fits_{}_IP'.format(len(filesndirs_H), len(filesndirs_K), i))
            time.sleep(2)
#-------------------------------------------------------------------------------
            if (args.mode == 1) or (args.mode == 2): #get IP & WaveSol
                print('Getting IP average values...')
                if (len(filesndirs_H) == 0) & (len(filesndirs_K) != 0):
                    IPval(i,'K')
                elif (len(filesndirs_H) != 0) & (len(filesndirs_K) == 0):
                    IPval(i,'H')
                else:
                    IPval(i,'H')
                    IPval(i,'K')
                print('DONE, saving under ../Temp/IP_X.txt')
                time.sleep(1)
        else:
            sys.exit(f'NO FILES FOUND UNDER ./A0_Fits/A0_Fits_{i}_IP/' )
#-------------------------------------------------------------------------------

    for i in tars:
        if os.path.isdir( f'./A0_Fits/A0_Fits_{i}_IP' ):
            continue
        else:
            sys.exit(f'NO FILES FOUND UNDER ./A0_Fits/A0_Fits_{i}_IP/' )

    if (args.mode == 1) or (args.mode == 3):
        print('-------------------------------------------------------------')
        print('Getting Wave Solutions...')
        sys.exit('Oops! not ready for this part yet!' )

    #     WVsol(tars, 'H')
    #     WVsol(tars, 'K')
    #
    #
    #     else:
    #         sys.exit('NO FILES FOUND UNDER ./A0_Fits/A0_Fits_{}/'.format(i) )
    # print('###############################################################')
    #
    #
    # print('Data Preparation Done!')
    # print('Preparation File Saved Under ./Temp/Prepdata/')
    print('###############################################################')
