
# coding: utf-8

import numpy as np
from scipy.interpolate import interp1d
import time

def rotint(w,s,vsini,eps=None,nr=None,ntheta=None,dif=None):
    #  This routine reads in a spectrum, s, on a wavelength scale, w, and a vsini
    #  with which to rotationally broaden the spectrum.  The rotationally broadened
    #  spectrum is returned in ns.  Parameters that can be set are the coefficient
    #  of the limb darkening law, eps (0.6 default), the number of radial steps on
    #  the disk, nr (default = 10), and the maximum number of steps in angle around
    #  the disk, ntheta (default = 100).  Final optional parameter dif allows for
    #  differential rotation according to the law Omeg(th)/Omeg(eq) = (1. - dif/2
    #  - (dif/2) cos(2 th)).  Dif = .675 nicely reproduces the law proposed by
    #  Smith, 1994, A&A, in press. to unify WTTS and CTTS.  Dif = .23 is similar to
    #  observed solar differential rotation.  Note: the th in the above expression
    #  is the stellar co-latitude, not the same as the integration variable used
    #  below.  This is a disk integration routine.
    #
    #  11-May-1994: Written CMJ.
    #

    if eps == None:
        eps = .6
    if nr == None:
        nr = 10
    if ntheta == None:
        ntheta = 100
    if dif == None:
        dif = 0

    tarea=0.0e00

    dr = 1./nr
    c  = 0
    for j in range(nr):
        r    = dr/2. + j*dr
        area = ( (r+dr/2.)**2 - (r-dr/2.)**2) / int(ntheta*r)*(1.-eps + eps*np.cos(np.arcsin(r)) )
        for k in range(int(ntheta*r)):
            c += 1
            th= np.pi/int(ntheta*r)+k*2.*np.pi/int(ntheta*r)
            if dif != 0:
                vl=vsini*r*np.sin(th)*(1.-dif/2.-dif/2.*np.cos(2.*np.arccos(r*np.cos(th))))
                # ns=ns+area*fspline(w+w*vl/3.e5,s,w)
                f    = interp1d(w+w*vl/3.e5,s, kind='linear',bounds_error=False,fill_value='extrapolate')
                ns   = ns+area*f(w)
                tarea= tarea+area
            else:
                vl = r*vsini*np.sin(th)
                # ns=ns+area*fspline(w+w*vl/3.e5,s,w)
                f = interp1d(w+w*vl/3.e5,s, kind='linear',bounds_error=False,fill_value='extrapolate')
                #if vl != 0:
                #    print '1',min(w+w*vl/3.e5),max(w+w*vl/3.e5)
                #    print '2',min(w),max(w)
                if j == 0 and k == 0:
                    ns = area*f(w)
                else:
                    ns=ns+area*f(w)
                tarea = tarea+area
    #print 'm1',max(ns)
    ns = ns/tarea
    return ns
