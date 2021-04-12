from scipy.signal import fftconvolve
from Engine.importmodule import *

def rotint(wave_spec,flux_spec,vrot):
    
    '''
    Applies rotational broadening to spectrum. This code is from the IvS Python Repository and is referenced as such (Institute for Astronomy at KU Leuven 2018, https://github.com/IvS-KULeuven/IvSPythonRepository)
    
    Inputs:
    wave_spec : Wavelength scale of spectrum
    flux_spec : Corresponding flux of spectrum
    vrot      : vsin(i)
    
    Outputs:
    
    wave_conv   : Rotationally broadened wavelength scale
    1-flux_conv : Rotationally broadened flux
    '''
    
    epsilon = 0.6

    wave_ = np.log(wave_spec)
    velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
    flux_ = np.interp(velo_,wave_,flux_spec)
    dvelo = velo_[1]-velo_[0]
    vrot = vrot/(2.99792e5)
    #-- compute the convolution kernel and normalise it
    n = np.int(2*vrot/dvelo)
    velo_k = np.arange(n)*dvelo
    velo_k -= velo_k[-1]/2.
    y = 1 - (velo_k/vrot)**2 # transformation of velocity
    G = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
    G /= np.sum(G)
    #-- convolve the flux with the kernel
    flux_conv = fftconvolve(1-flux_,G,mode='same')
    velo_ = np.arange(len(flux_conv))*dvelo+velo_[0]
    wave_conv = np.exp(velo_)
    return wave_conv,1-flux_conv
