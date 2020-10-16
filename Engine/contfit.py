import numpy as np
#import matplotlib.pyplot as plt
from scipy import interpolate
from Engine.detect_peaks import detect_peaks

def A0cont(a0wavecut,a0vcut,night,order):
    # a0vcut is a0fluxlist
    x = np.arange(len(a0vcut))

    # mpd: detect peaks that are at least separated by minimum peak distance
    peaks = detect_peaks(a0vcut, mpd=10)
    xtimes = 4
    mask = np.ones(len(peaks), dtype=bool)

    for a in range(len(peaks)):
        # double checking
        if a0vcut[peaks[a]] < .1:
            mask[a] = False
        if order == 6:
            if 1500 > x[peaks[a]] > 1350:
                mask[a] = False
    peaks = peaks[mask]

    for dothisxtimes in range(xtimes):
        mask = np.ones(len(peaks),dtype=bool)
        doer = iter(range(len(peaks)-2))
        for a in doer:
            if a0vcut[peaks[a+1]]-a0vcut[peaks[a]] < -0.2:
                mask[a]   = True
                mask[a+1] = False
                next(doer,None)
                continue
            elif a0vcut[peaks[a+1]]-a0vcut[peaks[a]] < 0.02 and a0vcut[peaks[a+1]]-a0vcut[peaks[a+2]] < 0.02:
                mask[a] = True
                mask[a+1] = False
                mask[a+2] = True
                next(doer,None)
                next(doer,None)
                continue
            else:
                mask[a] = True

        a0vcut = np.array(list(reversed(a0vcut)))
        peaks  = peaks[mask]

    c = 0

    for smoothing in np.arange(1e6,1e8,1e6):
        f = interpolate.UnivariateSpline(x[peaks], a0vcut[peaks], k=3, s=smoothing)
        continuum = f(x)
        peaks2 = detect_peaks(continuum)
        if len(peaks2) == 1:
            c += 1
        if c == 2:
            break

    '''
    if inparam.plotfigs == True:
        plt.figure(figsize=(20,12))
        plt.plot(x,a0vcut,color='black',alpha=.5)
        plt.scatter(x[peaks],a0vcut[peaks],s=25,color='blue')
        plt.plot(x,f(x),color='orange',alpha=.5)
        plt.savefig(inparam.outpath+'/A0 Contfit_'+str(night)+'_'+str(order))
        plt.clf()
        plt.close()
    '''

    return continuum
