
# coding: utf-8

import numpy as np
from scipy.interpolate import interp1d, splrep,splev
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
        eps = .4
    if nr == None:
        nr = 10
    if ntheta == None:
        ntheta = 20
    if dif == None:
        dif = 0

    c = 2.99792458e5
    tarea=0.0e00

    dr = 1./nr

    ns = np.zeros(len(s), dtype=np.float64)

    for j in range(nr):
        r=dr/2.+j*dr # step from dr/2 to 3dr/2, 5dr/2, etc til you hit edge
        area=((r+dr/2.)**2-(r-dr/2.)**2)/np.int(ntheta*r)*(1.-eps+eps*np.cos(np.arcsin(r))) # annulus for 0 to r, r to 2r, etc, times (1 - e + e*cos(arcsin(r))) (limb darkening effectively changes area?
        # divide by int(ntheta*r) because about to iterate through ntheta, higher r means a given ntheta corresponds to larger area, so...what
        for k in range(int(ntheta*r)):
            th=np.pi/int(ntheta*r)+k*2.*np.pi/int(ntheta*r)

            if dif != 0:
                vl = vsini*r*np.sin(th)*(1.-dif/2.-dif/2.*np.cos(2.*np.arccos(r*np.cos(th))))
            else:
                vl=r*vsini*np.sin(th)

            f = interp1d(w+w*vl/c,s, kind='linear',bounds_error=False,fill_value='extrapolate')
            ns = ns + area*f(w)
            tarea += area

    ns /= tarea

    return w, ns


# from numba import njit

# @njit
# def rotint(w,s,vsini,eps=None,nr=None,ntheta=None,dif=None):
#     #  This routine reads in a spectrum, s, on a wavelength scale, w, and a vsini
#     #  with which to rotationally broaden the spectrum.  The rotationally broadened
#     #  spectrum is returned in ns.  Parameters that can be set are the coefficient
#     #  of the limb darkening law, eps (0.6 default), the number of radial steps on
#     #  the disk, nr (default = 10), and the maximum number of steps in angle around
#     #  the disk, ntheta (default = 100).  Final optional parameter dif allows for
#     #  differential rotation according to the law Omeg(th)/Omeg(eq) = (1. - dif/2
#     #  - (dif/2) cos(2 th)).  Dif = .675 nicely reproduces the law proposed by
#     #  Smith, 1994, A&A, in press. to unify WTTS and CTTS.  Dif = .23 is similar to
#     #  observed solar differential rotation.  Note: the th in the above expression
#     #  is the stellar co-latitude, not the same as the integration variable used
#     #  below.  This is a disk integration routine.
#     #
#     #  11-May-1994: Written CMJ.
#     #
#
#     if eps == None:
#         eps = .4
#     if nr == None:
#         nr = 10
#     if ntheta == None:
#         ntheta = 20
#     if dif == None:
#         dif = 0
#
#     c = 2.99792458e5
#     tarea=0.0e00
#
#     dr = 1./nr
#
#     ns = np.zeros(len(s), dtype=np.float64)
#
#     for j in range(nr):
#         # step from dr/2 to 3dr/2, 5dr/2, etc til you hit edge
#         r = dr/2. + j*dr
#         # annulus for 0 to r, r to 2r, etc, times (1 - e + e*cos(arcsin(r))) (limb darkening effectively changes area?
#         area = ( (r+dr/2.)**2-(r-dr/2.)**2 ) / int(ntheta*r)*( 1.-eps + eps*np.cos(np.arcsin(r)) )
#         # divide by int(ntheta*r) because about to iterate through ntheta, higher r means a given ntheta corresponds to larger area, so...what
#         for k in range(int(ntheta*r)):
#             th = np.pi/int(ntheta*r) + k*2.*np.pi/int(ntheta*r)
#
#             if dif != 0:
#                 vl = r*vsini*np.sin(th) * ( 1.-dif/2.-dif/2. * np.cos( 2.*np.arccos(r*np.cos(th)) ) )
#             else:
#                 vl = r*vsini*np.sin(th)
#
#             # f = interp1d(w+w*vl/c,s, kind='linear',bounds_error=False,fill_value='extrapolate')
#             # ns = ns + area*f(w)
#             ns = ns + area * np.interp(w, w+w*vl/c, s)
#             tarea += area
#
#     ns /= tarea
#
#     return w, ns
