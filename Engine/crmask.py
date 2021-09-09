
from Engine.opt   import fmod, fmod_conti_tell
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
        clip_slope_tol = 40
        clip_pixel_tol = 12
    else:
        clip_slope_tol = 300
        clip_pixel_tol = 6

    fit,chi = fmod(parfit, fitobj)
    w,smod,cont,c2 = fmod_conti_tell(parfit, fitobj)
    continuum = cont*c2

    # Everywhere where data protrudes high above model, check whether slope surrounding protrusion is /\ and mask if sufficiently steep

    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    xdata = fitobj.x.copy(); sdata = fitobj.s.copy(); 

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
            if (np.min(slopeL) > clip_slope_tol) and (np.max(slopeR) < -clip_slope_tol) and (len(group) < clip_pixel_tol):
                CRmaskF = np.concatenate((CRmaskF,group))
        except ValueError:
            if (slopeL > clip_slope) and (slopeR < -clip_slope):
                CRmaskF = np.concatenate((CRmaskF,group))

    sflat = sdata/continuum
    sflat /= np.percentile(sflat,98)
    ind = np.where(sflat < 0.2)[0]
    if len(ind) > 0:
        for group in mit.consecutive_groups(ind):
            group = np.array(list(group))
            CRmaskF = np.concatenate((CRmaskF,group))

    return [xdata,CRmaskF]
