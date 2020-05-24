
import numpy as np
from astropy.io import fits
import os, glob
from Engine.rebin_jv import rebin_jv

### NOTE: THIS CODE DOES NOT DISTINGUISH BETWEEN H AND K BANDS ###

# Enter here the name of all A0_Fits directories you'd like to be used (multiple targets would be nice, for
# consistency, but only if they use the same orders.
sourcelist = ['A0_Fits_GJ281/']
HorK = 'H'
first = True;

xmaster = np.arange(2048)
wmaster = {}

for f in sourcelist:
    for a0 in glob.glob('./{}*.fits'.format(f)):
        hdulist = fits.open(a0)

        if HorK == 'K':
            orders = np.append(np.arange(2, 9), np.array([10, 11, 12, 13, 14, 16]))
            #orders = [3]
        elif HorK=='H':
            orders = np.append(np.arange(2, 7), np.array([10, 11, 13, 14, 16, 17, 20, 21, 22]))
#        orders = orders[:7]

        for o in range(len(orders)):
            tbdata = hdulist[o+1].data
            try:
                flag = np.array(tbdata['ERRORFLAG{}'.format(orders[o])])[0]
            except KeyError:
                print('Warning, file {} does not have the same set of orders as others!')
                continue

            if flag != 1:
                x = np.array(tbdata['X{}'.format(orders[o])])
                w = np.array(tbdata['WAVE{}'.format(orders[o])])
                x = x[(w != 0)]
                w = w[(w != 0)]
                w1 = rebin_jv(x,w,xmaster,False)

                if first == True:
                    wmaster[orders[o]] = w1
                else:
                    wmaster[orders[o]] = np.vstack((wmaster[orders[o]], w1))
            else:
                if first == True:
                    ww = np.zeros(len(w1)) ; ww[:] = np.nan
                    wmaster[orders[o]] = ww
                else:
                    ww = np.zeros(len(w1)) ; ww[:] = np.nan
                    wmaster[orders[o]] = np.vstack((wmaster[orders[o]], ww))

        first = False


filew = open('WaveSolns.csv','w')

for o in range(len(orders)):

    wmaster[orders[o]] = [np.nanmean(wmaster[orders[o]][:,i]) for i in range(len(wmaster[orders[o]][0,:]))]
    wmaster[orders[o]] = wmaster[orders[o]][5:-5]

    test = np.array(wmaster[orders[o]])
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
            filew.write('{},{},'.format(xmaster[i],wmaster[orders[o]][i]/1e4))
        else:
            filew.write('{},{}'.format(xmaster[i],wmaster[orders[o]][i]/1e4))
            filew.write('\n')
filew.close()



#-------------------
# IP average

for f in sourcelist:
    dump1 = 0
    for a0 in glob.glob('./{}*.fits'.format(f)):
        hdulist = fits.open(a0)

        if HorK == 'K':
            orders = np.append(np.arange(2, 9), np.array([10, 11, 12, 13, 14, 16]))
            #orders = [3]
        elif HorK=='H':
            orders = np.append(np.arange(2, 7), np.array([10, 11, 13, 14, 16, 17, 20, 21, 22]))

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

filew = open('IP.txt','w')
for o in np.arange(len(orders)):

    filew.write('{}: np.array([{:+1.2f}, {:+1.2f}, {:1.2}]),\n'.format(o, np.nanmean(IP14box[:, o]), np.nanmean(IP13box[:, o]), np.nanmean(IP5box[:, o]) ))
filew.close()
