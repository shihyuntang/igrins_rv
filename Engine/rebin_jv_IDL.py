def IDLspline (x,y,t,sigma):

    n = len(x)
    
    xx = x
    yy = y
    tt = t
    yp = np.zeros(2*n,dtype=float)

    delx1 = xx[1] - xx[0]
    dx1 = (yy[1] - yy[0])/delx1

    nm1 = n - 1
    np1 = n + 1

    delx2 = xx[2]-xx[1]
    delx12 = xx[2]-xx[0]
    c1 = -(delx12+delx1)/delx12/delx1
    c2 = delx12/delx1/delx2
    c3 = -delx1/delx12/delx2

    slpp1 = c1*yy[0]+c2*yy[1]+c3*yy[2]
    deln = xx[nm1]-xx[nm1-1]
    delnm1 = xx[nm1-1]-xx[nm1-2]
    delnn = xx[nm1]-xx[nm1-2]
    c1 = (delnn+deln)/delnn/deln
    c2 = -delnn/deln/delnm1
    c3 = deln/delnn/delnm1
    slppn = c3*yy[nm1-2]+c2*yy[nm1-1]+c1*yy[nm1]

    sigmap = sigma*nm1/(xx[nm1]-xx[0])
    dels = sigmap*delx1
    exps = np.exp(dels)
    sinhs = 0.5*(exps-1./exps)
    sinhin = 1./(delx1*sinhs)
    diag1 = sinhin*(dels*0.5*(exps+1./exps)-sinhs)
    diagin = 1./diag1
    yp[0] = diagin*(dx1-slpp1)
    spdiag = sinhin*(sinhs-dels)
    yp[n] = diagin*spdiag

    # Do as much work using vectors as possible.
    delx2 = xx[1:] - xx[:-1]
    dx2 = (yy[1:] - yy[:-1])/delx2
    dels = sigmap*delx2
    exps = np.exp(dels)
    sinhs = 0.5 *(exps-1./exps)
    sinhin = 1./(delx2*sinhs)
    diag2 = sinhin*(dels*(0.5*(exps+1./exps))-sinhs)
    diag2 = [0] + list(diag2[:-1] + diag2[1:])
    dx2nm1 = dx2[nm1-1] #; need to save this to calc yp[nm1]
    dx2 = [0] + list(dx2[1:] - dx2[:-1])
    spdiag = sinhin*(sinhs-dels)

    #; Need to do an iterative loop for this part.
    for i in range(1,nm1): 
        diagin = 1./(diag2[i] - spdiag[i-1]*yp[i+n-1])
        yp[i] = diagin*(dx2[i] - spdiag[i-1]*yp[i-1])
        yp[i+n] = diagin*spdiag[i]


    diagin = 1./(diag1-spdiag[nm1-1]*yp[n+nm1-1])
    yp[nm1] = diagin*(slppn-dx2nm1-spdiag[nm1-1]*yp[nm1-1])
    for i in range(n-2,-1,-1):
        yp[i] = yp[i] - yp[i+n]*yp[i+1]

    m = len(t)
    subs = np.ones(int(m),dtype=int)*nm1 #;subscripts
    s = xx[nm1]-xx[0]
    sigmap = sigma*nm1/s

    j = 0; done = False;
    
    for i in range(1,nm1+1):#;find subscript where xx[subs] > t(j) > xx[subs-1]
        while tt[j] < xx[i]:
            subs[j]=i
            j += 1
            if j == m:
                done = True
                break
        if done == True:
            break


    subs1 = subs - 1
    del1 = tt-xx[subs1]
    del2 = xx[subs] - tt
    dels = xx[subs]-xx[subs1]
    exps1 = np.exp(sigmap*del1)
    sinhd1 = 0.5*(exps1-1./exps1)
    exps = np.exp(sigmap*del2)
    sinhd2 = 0.5*(exps-1./exps)
    exps = exps1*exps
    sinhs = 0.5*(exps-1./exps)
    spl = (yp[subs]*sinhd1+yp[subs1]*sinhd2)/sinhs + \
        ((yy[subs]-yp[subs])*del1+(yy[subs1]-yp[subs1])*del2)/dels

    if m == 1:
        return spl[0]
    else:
        return spl
    
    


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
         
        S = IDLspline(Wold,Sold,W,0.01)     
        S /= XFac #take average in each pixel
        Snew = S.reshape(Nnew,XFac).sum(1)

    return Snew
