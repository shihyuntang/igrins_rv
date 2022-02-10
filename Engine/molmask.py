
from Engine.importmodule import *
from  Engine.rebin_jv import rebin_jv
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


def mask_wave2pixel_range(maskwaves, fitobj, w):
    """tansfer masking wave range into pixel range
    e.g., [[234500, 234510],[x, x]] --> [[100, 105],[x, x]]

    Args:
        maskwaves (list): masking wave range
        fitobj (class): fit object par save class
        w (arr): best fit model wavelength 

    Returns:
        list: masking pixel range
    """
    
    molmask_pixels_range = []
    for www in maskwaves:
        ind1 = fitobj.x[
            (np.abs(w-www[0]) == np.min(np.abs(w-www[0])))
            ][0]
        ind2 = fitobj.x[
            (np.abs(w-www[1]) == np.min(np.abs(w-www[1])))
            ][0]
        if ind2 == fitobj.x[0] or ind1 == fitobj.x[-1]:
            continue
        molmask_pixels_range = molmask_pixels_range + [[ind1,ind2]]

    return molmask_pixels_range


# in vac [AA]
kband_nonCO_line_centers = {
    4: {'Fe':[23690.1985, 23573.0969]},
    5: {
        'Na':[23385.5166, 23354.7948],
        'Mg':[23334.4565],
        'Sc':[23411.199],
        'Ti':[23447.8574],
        'Fe':[23331.5547, 23314.833],
        'HF':[23364.7043],
    },
    6: {
        'Si':[23147.9392],
        'Ca':[23059.9758],
        'Sc':[22992.5757],
        'Ti':[22969.5968],
        'Fe':[23150.908, 23051.0722, 23180.0529],
    }
}


def _look_sign_changes(a):

    sign_change_locs = (np.diff(np.sign(a)) != 0)*1
    if np.sum(sign_change_locs) == 0:
        return len(a)
    else:
        first_sign_change_loc = np.where(sign_change_locs==1)[0][:1]
        return first_sign_change_loc

def _apply_rv0_shift(line, rv0):
    c = 2.99792458e5
    return line*(1.+rv0/c)


def nonCO_masker(smod, w, cont, order, rv0, fitobj, flux_cut=0.99):

    flat_smod = smod / cont
    nonCO_line_centers = list(kband_nonCO_line_centers[order].values())
    nonCO_line_centers_flat = list(chain(*nonCO_line_centers))

    nonCO_mask_box = []
    for line_c in nonCO_line_centers_flat:
        
        line_c = _apply_rv0_shift(line_c, rv0)
        
        cen_x = np.nanargmin(np.abs(line_c-w) )

        # check if this line is in the template
        if flat_smod[cen_x] > flux_cut:
            continue

        # flux cut
        flux_over_cut_x = np.where(flat_smod > flux_cut)[0]
        flux_over_cut_x-=cen_x

        if sum(flux_over_cut_x<0) != 0:
            fl_l_bound = flux_over_cut_x[flux_over_cut_x<0][-1] + cen_x
        else: 
            fl_l_bound = 0
        if sum(flux_over_cut_x>0) != 0:
            fl_r_bound = flux_over_cut_x[flux_over_cut_x>0][0] + cen_x
        else:
            fl_r_bound = len(w)-1

        # # slope cut, check where sign changes
        # slopes = np.diff(flat_smod)
        
        # left_part = slopes[:cen_x]
        # right_part = slopes[cen_x+1:]

        # _loc = _look_sign_changes(left_part[::-1])
        # sl_l_bound = cen_x - _loc

        # _loc = _look_sign_changes(right_part)
        # sl_r_bound = cen_x + _loc

        # print('cen_x: ', cen_x)
        # print('fl_l_bound: ', fl_l_bound)
        # print('fl_f_bound: ', fl_r_bound)
        # print('sl_l_bound: ', sl_l_bound)
        # print('sl_f_bound: ', sl_r_bound)

        # use the bound closest to the center
        # l_bound = fl_l_bound if fl_l_bound>sl_l_bound[0] else sl_l_bound[0]
        # r_bound = fl_r_bound if fl_r_bound<sl_r_bound[0] else sl_r_bound[0]
        l_bound = fl_l_bound
        r_bound = fl_r_bound

        z_zero_point = fitobj.x[0]
        # print('zero point: ', z_zero_point)
        l_bound+=z_zero_point
        r_bound+=z_zero_point

        nonCO_mask_box.append([l_bound, r_bound])

    return nonCO_mask_box

def _group_consecutive_numbers(ind):
    
    maskwave_ranges = []
    for k, g in groupby(enumerate(ind), lambda x:x[0]-x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        maskwave_ranges.append(
            [group[0], group[-1]]
            )

    return maskwave_ranges


def merge_pixel_masks(mask1, mask2):

    if len(mask1)==0:
        return mask2
    if len(mask2)==0:
        return mask1

    mask1 = [
        list( range(int(i[0]),int(i[1])+1) ) for i in mask1
        ]
    mask2 = [
        list( range(int(i[0]),int(i[1])+1) ) for i in mask2
        ]

    mask1_flat = list(chain(*mask1))
    mask2_flat = list(chain(*mask2))

    mask = mask1_flat + mask2_flat
    mask = list(np.unique(mask))

    return _group_consecutive_numbers(mask)






