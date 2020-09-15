

import numpy as np
from scipy.interpolate import interp1d, splev, splrep
import time

def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray

def rebin_jv(Wold,Sold,Wnew,verbose):
    #Interpolates OR integrates a spectrum onto a new wavelength scale, depending
    #  on whether number of pixels per angstrom increases or decreases. Integration
    #  is effectively done analytically under a cubic spline fit to old spectrum.
    # Wold (input vector) old wavelngth scale.
    # Sold (input vector) old spectrum to be binned.
    # Wnew (input vector) new wavelength spectrum.
    # Snew (output vector) newly binned spectrum.
    #Edit History:
    # 10-Oct-90 JAV Create.
    # 22-Sep-91 JAV Translated from IDL to ANA.
    # 27-Aug-93 JAV Fixed bug in endpoint check: the "or" was essentially an "and".
    # 26-Aug-94 JAV	Made endpoint check less restrictive so that identical old and
    # new endpoints are now allowed. Switched to new Solaris library in call_external.

    #Determine spectrum attributes.
    Nold = int(len(Wold)) #number of old points
    Nnew = int(len(Wnew)) #number of new points
    PSold = (Wold[-1] - Wold[0]) / (Nold-1) #old pixel scale
    PSnew = (Wnew[-1] - Wnew[0]) / (Nnew-1) #new pixel scale

    #Verify that new wavelength scale is a subset of old wavelength scale.
    if verbose == True:
        if (Wnew[0] < Wold[0]) or (Wnew[-1] > Wold[-1]):
            print('New wavelength scale not subset of old.')

    #Select integration or interpolation depending on change in dispersion.
    if PSnew <= PSold:
        #pixel scale decreased. Interpolation by cubic spline.
        #Dummy = long(0)
        #Snew = spline(Wold,Sold,Wnew) # inerpolated spectrum
        f = interp1d(Wold,Sold,'cubic',bounds_error=False,fill_value='extrapolate') #interpolated old spectrum
        Snew = f(Wnew)
    else:
        #pixel scale increased. Integration under cubic spline.
        XFac = int(PSnew/PSold + 0.5) #pixel scale expansion factor
        #  Construct another wavelength scale (W) with a pixel scale close to that of
        #    the old wavelength scale (Wold), but with the additional constraint that
        #    every XFac pixels in W will exactly fill a pixel in the new wavelength
        #   scale (Wnew). Optimized for XFac < Nnew.
        
        dW = 0.5 * (Wnew[2:Nnew] - Wnew[0:Nnew-2]) #local pixel scale
        dW = np.concatenate((dW,[2*dW[Nnew-3] - dW[Nnew-4]])) #add trailing endpoint first
        dW = np.concatenate(([2*dW[0] - dW[1]],dW)) #add leading endpoint last
        
        W = np.empty((XFac,Nnew)) #initialize W as array
        for i in range(XFac): #loop thru subpixels
            W[i,:] = Wnew + dW*(float(2*i+1)/(2.0*XFac) - 0.5) #pixel centers in W
            
        W = np.transpose(W) #transpose W before Merging
        nIG = Nnew * XFac #elements in interpolation grid
        W = W.flatten() #make W into 1-dim vector
        #;  Interpolate old spectrum (Sold) onto wavelength scale W to make S. Then
        #;    sum every XFac pixels in S to make a single pixel in the new spectrum
        #;    (Snew). Equivalent to integrating under cubic spline through Sold.
        
        S = splev(W,splrep(Wold, Sold))         
        S /= XFac #take average in each pixel
        Snew = S.reshape(Nnew,XFac).sum(1)

    return Snew


