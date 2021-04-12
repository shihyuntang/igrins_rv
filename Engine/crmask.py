
from Engine.opt       import fmod
from Engine.importmodule import *
from Engine.detect_peaks import detect_peaks


def CRmasker(parfit,fitobj):
    '''
 
    Identify cosmic rays and hot pixels in spectrum, as well as places where the model does not have the ability to reflect the data.
    
    Inputs:
    parfit    : Best fit spectral model parameters
    fitobj    : Class containing data to be fit and stellar and telluric templates
    
    Outputs:
    CRmaskF : Pixels to be masked
    
    '''
    
    fit,chi = fmod(parfit, fitobj)

    # Everywhere where data protrudes high above model, check whether slope surrounding protrusion is /\ and mask if sufficiently steep
    residual = fitobj.s/fit
    MAD = np.median(np.abs(np.median(residual)-residual))
    CRmask = np.array(np.where(residual > np.median(residual)+2*MAD)[0])

    CRmaskF = []; CRmask = list(CRmask);

    for hit in [0,len(fitobj.x)-1]:
        if hit in CRmask:
            CRmaskF.append(hit)
            CRmask.remove(hit)
    CRmask = np.array(CRmask, dtype=np.int); CRmaskF = np.array(CRmaskF, dtype=np.int);

    for group in mit.consecutive_groups(CRmask):
        group = np.array(list(group))
        if len(group) == 1:
            gL = group-1; gR = group+1;
        else:
            peaks = detect_peaks(fitobj.s[group])
            if len(peaks) < 1:
                group = np.concatenate((np.array([group[0]-1]),group,np.array([group[-1]+1])))
                peaks = detect_peaks(fitobj.s[group])
                if len(peaks) < 1:
                    continue
            if len(peaks) > 1:
                continue
            gL = group[:peaks[0]]; gR = group[peaks[0]+1:];

        slopeL = (fitobj.s[gL+1]-fitobj.s[gL])/(fitobj.x[gL+1]-fitobj.x[gL])
        slopeR = (fitobj.s[gR]-fitobj.s[gR-1])/(fitobj.x[gR]-fitobj.x[gR-1])
        try:
            if (np.min(slopeL) > 300) and (np.max(slopeR) < -300) and len(group) < 6:
                CRmaskF = np.concatenate((CRmaskF,group))
        except ValueError:
            if (slopeL > 300) and (slopeR < -300):
                CRmaskF = np.concatenate((CRmaskF,group))

    return CRmaskF
