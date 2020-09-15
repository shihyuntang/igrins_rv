

import numpy as np
from scipy.interpolate import interp1d#, spline

def rebin(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements,
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged,
        if the new shape is bigger array elements are repeated
    See Also
    --------
    resize : Return a new array with the specified shape.
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M//m,n,N//n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m//M, axis=0), n//N, axis=1)

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

    #print '1 pre'
    #Select integration or interpolation depending on change in dispersion.
    if PSnew <= PSold:
        #pixel scale decreased. Interpolation by cubic spline.
        #print '1 post'
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
        W = np.empty((Nnew,XFac)) #initialize W as array
        for i in range(XFac): #loop thru subpixels
            W[:,i] = Wnew + dW*(float(2*i+1)/(2.0*XFac) - 0.5) #pixel centers in W
        W = np.transpose(W) #transpose W before Merging
        nIG = Nnew * XFac #elements in interpolation grid
        W = W.flatten() #make W into 1-dim vector
        #;  Interpolate old spectrum (Sold) onto wavelength scale W to make S. Then
        #;    sum every XFac pixels in S to make a single pixel in the new spectrum
        #;    (Snew). Equivalent to integrating under cubic spline through Sold.
        f = interp1d(Wold,Sold,'cubic',bounds_error=False,fill_value='extrapolate') #interpolated old spectrum
        S = f(W)
        #S = spline(Wold,Sold,W)
        S = S / XFac #take average in each pixel
        Sdummy = S.reshape(XFac,Nnew) #initialize Sdummy as array
        Snew = XFac * rebin(Sdummy,(1,Nnew)) #most efficient pixel sum
        Snew = Snew.flatten() #convert back to vector

    return Snew
