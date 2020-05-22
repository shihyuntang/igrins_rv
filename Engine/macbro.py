
# coding: utf-8


import numpy as np
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian1DKernel, convolve

def macbro(w,s,hwhm):
    #Smooths a spectrum by convolution with a gaussian of specified hwhm.
    # w (input vector) wavelength scale of spectrum to be smoothed
    # s (input vector) spectrum to be smoothed
    # hwhm (input scalar) half width at half maximum of smoothing gaussian.
    #Returns a vector containing the gaussian-smoothed spectrum.
    #Edit History:
    #  -Dec-90 GB,GM Rewrote with fourier convolution algorithm.
    #  -Jul-91 AL	Translated from ANA to IDL.
    #22-Sep-91 JAV	Relaxed constant dispersion check; vectorized, 50% faster.
    #05-Jul-92 JAV	Converted to function, handle nonpositive hwhm.

    #;Warn user and return input argument if hwhm is negative.
    if hwhm <= 0.0:
        print('Warning! Forcing negative smoothing width to zero.')
        return s

    #Calculate (uniform) dispersion.
    nw = len(w) # points in spectrum
    dw = (w[-1] - w[0]) / (nw-1) #wavelength change per pixel

    #Make smoothing gaussian; extend to 4 sigma.
    #Note: 4.0 / sqrt(2.0*alog(2.0)) = 3.3972872 and sqrt(alog(2.0))=0.83255461
    # sqrt(alog(2.0)/pi)=0.46971864 (*1.0000632 to correct for >4 sigma wings)
    nhalf = abs(int(3.3972872 * hwhm/dw)) # points in half gaussian
    ng = 2 * nhalf + 1 # points in gaussian (odd!)
    wg = dw * (np.arange(ng,dtype=float) - (ng-1)/2.0) #wavelength scale of gaussian
    xg = (0.83255461 / hwhm) * wg  #convenient absisca
    gpro = (0.46974832 * dw / hwhm) * np.exp(-xg*xg) #unit area gaussian w/ FWHM
    gpro = gpro / np.sum(gpro)

    #Pad spectrum ends to minimize impact of Fourier ringing.
    npad = nhalf + 2 # pad pixels on each end
    spad = np.concatenate((s[0]*np.ones(npad),s,np.ones(npad)*s[-1]))

    #Convolve and trim.

    #k = Gaussian1DKernel(stddev=hwhm/(np.sqrt(2*np.log(2))))  (bad)
    sout = convolve(spad,gpro)
    #sout = np.convolve(spad,gpro) #convolve with gaussian (gives weird offset to s)
    sout = sout[npad:npad+nw] #trim to original data/length
    sout = sout * np.sum(s) / np.sum(sout)
    return sout #return broadened spectrum.
