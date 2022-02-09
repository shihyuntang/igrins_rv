
from Engine.importmodule import *
from   Engine.rebin_jv import rebin_jv
import matplotlib.pyplot as plt

def _group_consecutive_numbers_and_get_mask_ranges(ind, waveH2O, halfsep):
    
    maskwave_ranges = []
    for k, g in groupby(enumerate(ind), lambda x:x[0]-x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        maskwave_ranges.append(
            [waveH2O[group[0]]-halfsep, waveH2O[group[-1]]+halfsep]
            )
    return maskwave_ranges

def h2o_masker(
    inparam, args, order, night, watm, satm, molnames, molwaves, 
    molfluxes):

    go = True
    for mol in molnames:

        watmI, satmI = molwaves[mol].copy(), molfluxes[mol].copy()

        if mol =='h2o':#!= dominantmol:
            fluxH2O = np.array(satmI.copy())
            waveH2O = np.array(watmI.copy())
        else:
            if go:
                fluxbox = satmI.copy()
                tempwave = watmI.copy()
                go = False
            else:
                fluxnew = rebin_jv(watmI, satmI, tempwave, False)
                fluxnew[(watmI < tempwave[0]) | (watmI > tempwave[-1])] = np.nan
                fluxbox = np.vstack((fluxbox, fluxnew))

    fluxother = np.array(
        [np.nanmin(fluxbox[:,ii]) for ii in range(len(fluxbox[0,:]))]
        )
    fluxother = rebin_jv(tempwave, fluxother, waveH2O, False)
    fluxother[(tempwave[0] > waveH2O) | (tempwave[-1] < waveH2O)] = 0.

    ind = np.where(fluxother > fluxH2O)[0]
    halfsep = 0.5*np.median(np.diff(waveH2O))
    
    maskwaves = _group_consecutive_numbers_and_get_mask_ranges(
        ind, waveH2O, halfsep
        )


    '''
    if len(maskwaves) > 3:
        maskwaves0 = maskwaves.copy(); maskwaves = [maskwaves0[0]];
        for m in range(1,len(maskwaves0)):
            if abs(maskwaves0[m][0] - maskwaves[-1][1]) < 4:
                maskwaves[-1] = [maskwaves[-1][0],maskwaves0[m][1]]
            else:
                maskwaves.append(maskwaves0[m])
    '''


    return maskwaves
