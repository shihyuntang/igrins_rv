from scipy.signal import fftconvolve
from Engine.importmodule import *

def rotint_fast(wave_spec,flux_spec,vrot):

    epsilon = 0.6

    wave_ = np.log(wave_spec)
    velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
    flux_ = np.interp(velo_,wave_,flux_spec)
    dvelo = velo_[1]-velo_[0]
    vrot = vrot/(2.99792e5)
    #-- compute the convolution kernel and normalise it
    n = int(2*vrot/dvelo)
    velo_k = np.arange(n)*dvelo
    velo_k -= velo_k[-1]/2.
    y = 1 - (velo_k/vrot)**2 # transformation of velocity
    G = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
    G /= G.sum()
    #-- convolve the flux with the kernel
    flux_conv = fftconvolve(1-flux_,G,mode='same')
    velo_ = np.arange(len(flux_conv))*dvelo+velo_[0]
    wave_conv = np.exp(velo_)
    return wave_conv,1-flux_conv
