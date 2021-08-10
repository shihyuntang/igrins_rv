
from Engine.opt   import fmod
# from Engine.opt_rebintel   import fmod
from Engine.importmodule import *
from Engine.detect_peaks import detect_peaks
from Engine.rebin_jv import rebin_jv


def CRmasker(parfit, fitobj, tel=False):
    '''
    Identify cosmic rays and hot pixels in spectrum, as well as places where the model does not have the ability to reflect the data.

    Inputs:
    parfit    : Best fit spectral model parameters
    fitobj    : Class containing data to be fit and stellar and telluric templates

    Outputs:
    CRmaskF : Pixels to be masked
    '''

    if tel:
        clip_slope_tol = 100
        clip_pixel_tol = 8
    else:
        clip_slope_tol = 300
        clip_pixel_tol = 6

    fit,chi = fmod(parfit, fitobj)

    # Everywhere where data protrudes high above model, check whether slope surrounding protrusion is /\ and mask if sufficiently steep

    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    dstep = np.median(w[1:]-w[:-1])
    nstep = int((w[-1]-w[0])/dstep)
    wreg = np.linspace(w[0],w[-1],nstep)
    sdata = rebin_jv(w,fitobj.s,wreg,False)
    udata = rebin_jv(w,fitobj.u,wreg,False)
    xdata = np.linspace(fitobj.x[0],fitobj.x[-1],nstep)
    w = wreg.copy()

    residual = sdata/fit
    MAD = np.median(np.abs(np.median(residual)-residual))
    CRmask = np.array(np.where(residual > np.median(residual)+2*MAD)[0])

    CRmaskF = []; CRmask = list(CRmask);

    for hit in [0,len(xdata)-1]:
        if hit in CRmask:
            CRmaskF.append(hit)
            CRmask.remove(hit)
    CRmask = np.array(CRmask, dtype=np.int); CRmaskF = np.array(CRmaskF, dtype=np.int);

    for group in mit.consecutive_groups(CRmask):
        group = np.array(list(group))
        if len(group) == 1:
            gL = group-1; gR = group+1;
        else:
            peaks = detect_peaks(sdata[group])
            if len(peaks) < 1:
                group = np.concatenate((np.array([group[0]-1]),group,np.array([group[-1]+1])))
                peaks = detect_peaks(sdata[group])
                if len(peaks) < 1:
                    continue
            if len(peaks) > 1:
                continue
            gL = group[:peaks[0]]; gR = group[peaks[0]+1:];

        slopeL = (sdata[gL+1]-sdata[gL])/(xdata[gL+1]-xdata[gL])
        slopeR = (sdata[gR]-sdata[gR-1])/(xdata[gR]-xdata[gR-1])
        try:
            if (np.min(slopeL) > clip_slope) and (np.max(slopeR) < -clip_slope) and (len(group) < clip_pixel_tol):
                CRmaskF = np.concatenate((CRmaskF,group))
        except ValueError:
            if (slopeL > clip_slope) and (slopeR < -clip_slope):
                CRmaskF = np.concatenate((CRmaskF,group))

    #mask = np.ones_like(xdata,dtype=bool)

    return [xdata,CRmaskF]
