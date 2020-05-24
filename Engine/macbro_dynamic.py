
import numpy as np
from scipy.interpolate import interp1d,splev,splrep
from astropy.convolution import Gaussian1DKernel, convolve

def macbro_dyn(w,s,hwhmlist):

    hwhm_ref = min(hwhmlist)

    tmp = hwhmlist[1:]
    hwhmlist = tmp + (hwhmlist[0:-1]-tmp)/2.
    
    tmp = w[1:]
    dw = tmp-w[0:-1]
    dw_stretched = dw/(hwhmlist/hwhm_ref)
    w_stretched = np.concatenate((np.array([w[0]]),w[0]+np.cumsum(dw_stretched)))
    
    dw_fine = min(dw)
    w_fine = np.arange(min(w_stretched),max(w_stretched)+dw_fine,dw_fine)
    ####
       
    spl = splrep(w_stretched,s)
    s_fine = splev(w_fine,spl)
    
    ####
    nw = len(w_fine) # points in spectrum
    nhalf = abs(int(3.3972872 * hwhm_ref/dw_fine)) # points in half gaussian
    ng = 2 * nhalf + 1 # points in gaussian (odd!)
    wg = dw_fine * (np.arange(ng,dtype=float) - (ng-1)/2.0) #wavelength scale of gaussian
    xg = (0.83255461 / hwhm_ref) * wg  #convenient absisca
    gpro = (0.46974832 * dw_fine / hwhm_ref) * np.exp(-xg*xg) #unit area gaussian w/ FWHM
    gpro = gpro / np.sum(gpro)
    ####

    npad = nhalf + 2 # pad pixels on each end
    spad = np.concatenate((s_fine[0]*np.ones(npad),s_fine,np.ones(npad)*s_fine[-1]))
    sout = convolve(spad,gpro)
    sout = sout[npad:npad+nw] #trim to original length
    try:
        spl = splrep(w_fine,sout)
    except TypeError:
        print(min(hwhmlist),min(dw),min(w_stretched),max(w_stretched))
        print(len(w_fine),len(w),len(w_stretched))
        print(breaker)
    sout0 = splev(w_stretched,spl)
    sout0 = sout0 * np.sum(s) / np.sum(sout0)

    return sout0 

