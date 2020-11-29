import numpy as np
from scipy import interpolate
from Engine.detect_peaks import detect_peaks
#import matplotlib.pyplot as plt

def A0cont(a0wavecut,a0vcut,night,order,band):
    # a0vcut is a0fluxlist
    x = np.arange(len(a0vcut))

    # mpd: detect peaks that are at least separated by minimum peak distance
    peaks = detect_peaks(a0vcut, mpd=10)
    if band == 'H':
        xtimes = 1
    else:
        xtimes = 3

    for ii in range(xtimes):
        mask = np.ones(len(peaks), dtype=bool)

        f = np.polyfit(x[peaks],a0vcut[peaks],4)
        q = np.poly1d(f)
        residual = a0vcut[peaks]-q(x[peaks])
        MAD = np.median(np.abs(residual-np.median(residual)))

        '''
        plt.figure(figsize=(20,12))
        plt.plot(x,a0vcut,color='black',alpha=.5)
        plt.scatter(x[peaks],a0vcut[peaks],s=25,color='blue')
        plt.plot(x[peaks],q(x[peaks]),color='red')
        plt.plot(x[peaks],q(x[peaks])+3*MAD,color='orange')
        plt.plot(x[peaks],q(x[peaks])-5*MAD,color='orange')
        plt.savefig(inparam.outpath+'/A0 Contfit_'+str(night)+'_'+str(order)+'_'+str(masterbeam)+'_0')
        plt.clf()
        plt.close()
        '''

        mask[(a0vcut[peaks]/np.nanmedian(a0vcut) < .1)] = False
        mask[(a0vcut[peaks] < q(x[peaks])-5*MAD)] = False
        mask[(a0vcut[peaks] > q(x[peaks])+3*MAD)] = False
        peaks = peaks[mask]

    c = 0

    for smoothing in np.arange(1e6,1e8,1e6):
        f = interpolate.UnivariateSpline(x[peaks], a0vcut[peaks], k=3, s=smoothing)
        continuum = f(x)
        peaks2 = detect_peaks(continuum)
        if len(peaks2) == 1:
            c += 1
        if c == 2:
            break

    if smoothing == 99000000.0:
        for smoothing in np.arange(1e8,1e10,1e8):
            f = interpolate.UnivariateSpline(x[peaks], a0vcut[peaks], k=3, s=smoothing)
            continuum = f(x)
            peaks2 = detect_peaks(continuum)
            if len(peaks2) == 1:
                c += 1
            if c == 2:
                break

    '''
    plt.figure(figsize=(20,12))
    plt.plot(x,a0vcut,color='black',alpha=.5)
    plt.scatter(x[peaks],a0vcut[peaks],s=25,color='blue')
    plt.plot(x,f(x),color='orange',alpha=.5)
    plt.savefig(inparam.outpath+'/A0 Contfit_'+str(night)+'_'+str(order)+'_'+str(masterbeam)+'_1')
    plt.clf()
    plt.close()
    '''

    return continuum
