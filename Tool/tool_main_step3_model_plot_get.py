import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *
from Engine.importmodule import read_prepdata

from Engine.IO_AB      import setup_templates, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs,inparams
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod, fmod_conti
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def rv_MPinst(args, inparam, orders, order_use, trk, step2or3, i):

    # Main function for RV fitting that will be threaded over by multiprocessing

    nights   = inparam.nights
    night = nights[i] # current looped night

    order = orders[order_use]
    xbounds = inparam.xbounddict[order]
    firstorder = orders[0]          # First order that will be analyzed, related to file writing

    print('Working on order {:02d}/{:02d} ({}), night {:03d}/{:03d} ({}) PID:{}...'.format(int(order_use)+1,
                                                                                           len(orders),
                                                                                           order,
                                                                                           i+1,
                                                                                           len(inparam.nights),
                                                                                           night,
                                                                                           mp.current_process().pid) )

    #-------------------------------------------------------------------------------

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = [];
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    nightsout = [];

    wminibox      = np.ones(2048)
    sminibox      = np.ones(2048)
    flminibox_tel = np.ones(2048)
    flminibox_ste = np.ones(2048)
    contiminibox  = np.ones(2048)
    flminibox_mod = np.ones(2048)

    wminibox[:]     = np.nan
    sminibox[:]     = np.nan
    flminibox_tel[:]= np.nan
    flminibox_ste[:]= np.nan
    contiminibox[:] = np.nan
    flminibox_mod[:]  = np.nan


    for t in tagsnight:
        nightsout.append(night)

#-------------------------------------------------------------------------------
    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == float:
        initguesses = inparam.initguesses
    else:
        sys.exit('ERROR! EXPECING SINGAL NUMBER OR FILE FOR INITGUESSES! QUITTING!')

    if np.isnan(initguesses) == True:
        logger.warning(f'  --> Previous run of {night} found it inadequate, skipping...')
        return nightsout, rvsminibox, parfitminibox, vsiniminibox, tagsminibox

    # start at bucket loc = 1250 +- 100, width = 250 +- 100, depth = 100 +- 5000 but floor at 0
    if args.band == 'H':
        centerloc = 1250
    else:
        centerloc = 1180

#-------------------------------------------------------------------------------
    ### Initialize parameter array for optimization as well as half-range values for each parameter during the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and in between optimization steps.
    pars0 = np.array([np.nan,                                                # 0: The shift of the stellar template (km/s) [assigned later]
                      0.3,                                                   # 1: The scale factor for the stellar template
                      0.0,                                                   # 2: The shift of the telluric template (km/s)
                      0.6,                                                   # 3: The scale factor for the telluric template
                      inparam.initvsini,                                     # 4: vsini (km/s)
                      np.nan,                                                # 5: The instrumental resolution (FWHM) in pixels
                      np.nan,                                                # 6: Wavelength 0-pt
                      np.nan,                                                # 7: Wavelength linear component
                      np.nan,                                                # 8: Wavelength quadratic component
                      np.nan,                                                # 9: Wavelength cubic component
                      1.0,                                                   #10: Continuum zero point
                      0.,                                                    #11: Continuum linear component
                      0.,                                                    #12: Continuum quadratic component
                      np.nan,                                                #13: Instrumental resolution linear component
                      np.nan,                                                #14: Instrumental resolution quadratic component
                      centerloc,                                             #15: Blaze dip center location
                      330,                                                   #16: Blaze dip full width
                      0.05,                                                  #17: Blaze dip depth
                      90,                                                    #18: Secondary blaze dip full width
                      0.05,                                                  #19: Blaze dip depth
                      0.0,                                                   #20: Continuum cubic component
                      0.0,                                                   #21: Continuum quartic component
                      0.0,                                                   #22: Continuum pentic component
                      0.0])                                                  #23: Continuum hexic component

    # This one specific order is small and telluric dominated, start with greater stellar template power to ensure good fits
    if int(order) == 13:
        pars0[1] = 0.8

    # Iterate over all A/B exposures
    for t in [0]:
        tag = tagsnight[t]
        beam = beamsnight[t]
        masterbeam = beam

        # Load synthetic telluric template generated during Step 1
        # [:8] here is to ensure program works under Night_Split mode

        # Use instrumental profile dictionary corresponding to whether IGRINS mounting was loose or not
        if np.int(night[:8]) < 20180401 or np.int(night[:8]) > 20190531:
            IPpars = inparam.ips_tightmount_pars[args.band][masterbeam][order]
        else:
            IPpars = inparam.ips_loosemount_pars[args.band][masterbeam][order]

        if beam == 'A':
            antibeam = 'B'
        elif beam == 'B':
            antibeam = 'A'
        else:
            sys.exit('uhoh')

        A0loc = f'../Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

        try:
            hdulist = fits.open(A0loc)
        except IOError:
            logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
            return wminibox,sminibox,flminibox_mod,flminibox_tel,flminibox_ste,contiminibox

        # Find corresponding table in fits file, given the tables do not go sequentially by order number due to multiprocessing in Step 1
        num_orders = 0
        for i in range(25):
            try:
                hdulist[i].columns[0].name[9:]
                num_orders += 1
            except:
                continue

        fits_layer = [ i for i in np.arange(num_orders)+1 if np.int(hdulist[i].columns[0].name[9:]) == order ][0]

        tbdata = hdulist[ fits_layer ].data
        flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

        # Check whether Telfit hit critical error in Step 1 for the chosen order with this night. If so, skip.
        if flag == 1:
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}, skipping...')
            return wminibox,sminibox,flminibox_mod,flminibox_tel,flminibox_ste,contiminibox


        watm = tbdata['WATM'+str(order)]
        satm = tbdata['SATM'+str(order)]
        a0contx    = tbdata['X'+str(order)]
        continuum  = tbdata['BLAZE'+str(order)]

        # Remove extra rows leftover from having columns of unequal length
        satm = satm[(watm != 0)]
        watm = watm[(watm != 0)]
        satm[(satm < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
        a0contx = a0contx[(continuum != 0)]
        continuum = continuum[(continuum != 0)]

        # Retrieve pixel bounds for where within each other significant telluric absorption is present.
        # If these bounds were not applied, analyzing some orders would give garbage fits.
        if args.band=='K':
            if int(order) in [13, 14]:
                bound_cut = inparam.bound_cut_dic[args.band][order]
            else:
                bound_cut = [150, 150]

        elif args.band=='H':
            if int(order) in [6, 10, 11, 13, 14, 16, 17, 20, 21, 22]:
                bound_cut = inparam.bound_cut_dic[args.band][order]
            else:
                bound_cut = [150, 150]

        # Load target spectrum
        x,wave,s,u = init_fitsread(f'{inparam.inpath}/{night}/{beam}/',
                                    'target',
                                    'separate',
                                    night,
                                    order,
                                    tag,
                                    args.band,
                                    bound_cut)

        #-------------------------------------------------------------------------------

        # Execute S/N cut
        s2n = s/u
        if np.nanmedian(s2n) < np.float(args.SN_cut):
            logger.warning('  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, beam, tag))
            continue

        # Trim obvious outliers above the blaze (i.e. cosmic rays)
        nzones = 5
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);
        x = basicclip_above(x,s,nzones); wave = basicclip_above(wave,s,nzones); u = basicclip_above(u,s,nzones); s = basicclip_above(s,s,nzones);

        # Cut spectrum to within wavelength regions defined in input list
        s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
        wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
        x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

        # Trim stellar template to relevant wavelength range
        mwave_in,mflux_in = stellarmodel_setup(wave_piece,inparam.mwave0,inparam.mflux0)

        # Trim telluric template to relevant wavelength range
        satm_in = satm[(watm > np.min(wave_piece)*1e4 - 11) & (watm < np.max(wave_piece)*1e4 + 11)]
        watm_in = watm[(watm > np.min(wave_piece)*1e4 - 11) & (watm < np.max(wave_piece)*1e4 + 11)]

        # Make sure data is within telluric template range (shouldn't do anything)
        s_piece    = s_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        u_piece    = u_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        x_piece    = x_piece[   (wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]
        wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+5) & (wave_piece*1e4 < np.max(watm_in)-5)]

        # Normalize continuum from A0 to flux scale of data
        continuum /= np.nanmedian(continuum)
        continuum *= np.nanpercentile(s_piece,99)

        # --------------------------------------------------------------

        par = pars0.copy()

        # Get initial guess for cubic wavelength solution from reduction pipeline
        f = np.polyfit(x_piece,wave_piece,3)
        par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;
        par[9] = par9in ; par[8] = par8in ; par[7] = par7in ; par[6] = par6in

        par[0] = initguesses-inparam.bvcs[night+tag] # Initial RV with barycentric correction
        par[5] = IPpars[2]; par[13] = IPpars[1]; par[14] = IPpars[0];


        # Arrays defining parameter variations during optimization steps
        #                            | 0    1    2    3 |  | ------ 4 ------ |  | 5 |   | 6     7     8           9  |  |10  11  12| |13 14|  |15   16   17   18    19 |  |20   21   22   23 |
        dpars = {'cont' : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    1e7, 1, 1,   0, 0,    10., 20., 0.2, 50.0, 0.2,   1.0, 1.0, 1.0, 1.0 ]),
                 'twave': np.array([  0.0, 0.0, 0.0, 1.0,   0.0,                 0.0,   10.0, 10.0,  5.00000e-5, 1e-7,   0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'ip'   : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,                 0.5,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 's'    : np.array([  5.0, 1.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'v'    : np.array([  0.0, 0.0, 0.0, 0.0,   inparam.vsinivary,   0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'ts'   : np.array([  5.0, 1.0, 0.0, 1.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ])}
        if masterbeam == 'B':
            dpars['cont'] = np.array([0.0, 0.0, 0.0, 0.0,   0.0,                 0.0,    0.0,  0.0,  0.0,        0.0,    1e7, 1, 1,   0, 0,     0.,  0., 0.0,  0.,  0.0,   1.0, 1.0 , 1.0, 1.0 ])

        # Use quadratic blaze correction for order 13; cubic for orders 6, 14, 21; quartic for orders 16 and 22
        if args.band == 'H':
            if np.int(order) in [13]:
                dpars['cont'][20] = 0.; dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
            elif np.int(order) in [6,14,21]:
                dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
            else:
                pass
        else:
            if np.int(order) in [3]:
                dpars['cont'][20] = 0.; dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
            elif np.int(order) in [4,5]:
                dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
            elif np.int(order) in [6]:
                dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
            else:
                pass

        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]),masterbeam,np.array([],dtype=int))

        #-------------------------------------------------------------------------------

        # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
        optimize = True
        par_in = par.copy()
        if masterbeam == 'B':
            hardbounds = [par_in[4] - dpars['v'][4],      par_in[4] + dpars['v'][4],
                          par_in[5] - dpars['ip'][5],     par_in[5]+dpars['ip'][5]]
        else:
            hardbounds = [par_in[4] - dpars['v'][4],      par_in[4] + dpars['v'][4],
                          par_in[5]  - dpars['ip'][5],    par_in[5] + dpars['ip'][5],
                          par_in[15] - dpars['cont'][15], par_in[15] + dpars['cont'][15],
                          par_in[16] - dpars['cont'][16], par_in[16] + dpars['cont'][16],
                          0.,                             par_in[17] + dpars['cont'][17],
                          par_in[18] - dpars['cont'][18], par_in[18] + dpars['cont'][18],
                          0.,                             par_in[19] + dpars['cont'][19]
                          ]
        if hardbounds[0] < 0.5:
            hardbounds[0] = 0.5
        if hardbounds[2] < 1:
            hardbounds[2] = 1

        # Begin optimization. Fit the blaze, the wavelength solution, the telluric template power and RV, the stellar template power and RV, the
        # zero point for the instrumental resolution, and the vsini of the star separately, iterating and cycling between each set of parameter fits.

        cycles = 4

        optgroup = ['cont', 'twave', 'cont', 'ts',
                    'cont', 'twave', 's', 'cont',
                    'twave',
                    'ip', 'v',
                    'ip', 'v',
                    'twave',  's',
                    'twave',  'ts']

        nk = 1
        for nc, cycle in enumerate(np.arange(cycles), start=1):
            if cycle == 0:
                parstart = par_in.copy()

            for optkind in optgroup:
                parfit_1 = optimizer(parstart, dpars[optkind], hardbounds, fitobj, optimize)
                parstart = parfit_1.copy()
                if args.debug == True:
                    outplotter_23(parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(order,night,tag,nk,optkind), trk, inparam, args, step2or3, order)
                    logger.debug(f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
                nk += 1

            ## After first cycle, use best fit model to identify CRs/hot pixels
            if nc == 1:
                parfit = parfit_1.copy()
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

                fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]),masterbeam,CRmaskF)

        parfit = parfit_1.copy()

        #-------------------------------------------------------------------------------

        # if best fit stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            logger.warning(f'  --> parfit[1] < 0.1, {night} parfit={parfit}')
            continue

        # if best fit stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
            logger.warning(f'  --> parfit[1] == par_in[1] or parfit[3] == par_in[3], {night}')
            continue

        # if best fit model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq = fmod(parfit,fitobj)
        if len(smod[(smod < 0)]) > 0:
            logger.warning(f'  --> len(smod[(smod < 0)]) > 0, {night}')
            continue

        #-------------------------------------------------------------------------------

        # Compute model and divide for residual
        fullmodel,chisq = fmod(parfit,fitobj)

        # Set both stellar and telluric template powers to 0 to compute only continuum
        parcont = parfit.copy();
        parcont[1] = 0.; parcont[3] = 0.;
        contmodel, chisq = fmod(parcont,fitobj)

        # Set stellar tempalte power to 0 to compute only telluric, and vice versa
        parS = parfit.copy(); parT = parfit.copy();
        parT[1] = 0.; parS[3] = 0.;
        stellmodel,chisq = fmod(parS,   fitobj)
        tellmodel, chisq = fmod(parT,   fitobj)

        # Divide everything by continuum model except residual
        dataflat  = fitobj.s/contmodel
        modelflat = fullmodel/contmodel
        stellflat = stellmodel/contmodel
        tellflat  = tellmodel/contmodel

    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    wminibox[:len(w)]        = w
    sminibox[:len(w)]        = dataflat
    flminibox_mod[:len(w)]   = modelflat
    flminibox_tel[:len(w)]   = tellflat
    flminibox_ste[:len(w)]   = stellflat
    contiminibox[:len(w)]    = contmodel
    # residualbox[:len(w)]     = residual

    # Save results in fits file
    c1 = fits.Column(name='wavelength',    array=wminibox,         format='D')
    c2 = fits.Column(name='s',             array=sminibox,        format='D')
    c3 = fits.Column(name='model_fl',      array=flminibox_mod,       format='D')
    c4 = fits.Column(name='tel_fl',        array=flminibox_tel,       format='D')
    c5 = fits.Column(name='ste_fl',        array=flminibox_ste,       format='D')
    c6 = fits.Column(name='conti_fl',      array=contiminibox,     format='D')


    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    if order == firstorder::  # If first time writing fits file, make up filler primary hdu
        bleh = np.ones((3, 3))
        primary_hdu1 = fits.PrimaryHDU(bleh)
        hdul = fits.HDUList([primary_hdu1, hdu_1])
        hdul.writeto(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}_{}.fits'.format(args.targname, args.band,night,tag ))
    else:
        hh = fits.open(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}_{}.fits'.format(args.targname, args.band,night,tag ))
        hh.append(hdu_1)
        hh.writeto(inparam.outpath+'/'+name+'/RVresultsRawBox_fit_wl_{}_{}_{}_{}.fits'.format(args.targname,  args.band,night,tag ), overwrite=True)

    return wminibox,sminibox,flminibox_mod,flminibox_tel,flminibox_ste,contiminibox


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 3',
                                     description = '''
                                     Performs a full analysis of each target star observation to produce accurate and precise RVs. \n
                                     All the wavelength regions defined in Step 1 are used, and the code analyzes each observation that is part of a given exposure separately. \n
                                     Unless the target vsini is already known to high accuracy, an initial run of Step 3 in which \vsini is allowed to vary is required. \n
                                     This provides an estimate of vsini that can then be plugged into the code as a fixed value in the second run of Step 3. \n
                                     If the user seeks the best possible RV uncertainty estimates, or if their target star has a relatively high \vsini ($>$ 10 \kms), they must run Step 3 once with \vsini held fixed at its estimated value and once with \vsini held fixed at this value plus or minus one sigma. \n
                                     The minor differences in the RVs of the two runs (as low as $<$1 \ms and as high as 7 \ms) can then be incorporated into the final uncertainties. \n
                                     If \vsini is already well-known, it is not necessary to run Step 3 more than once, as the code fully converges to the final RVs (within uncertainty) through just one run.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-mode",    dest="mode",             action="store",
                        help="RV standard star (STD) or a normal target (TAR)?",
                        type=str,   default='')
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))
    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 ",
                        type=str,   default='50')
    parser.add_argument("-nAB",      dest="nAB",           action="store",
                        help="Minium number of separte A/B exposures within a set for a given observation (ensures accuracy of uncertainy estimates). Default = 2 for STD, 3 for TAR",
                        type=str,   default='')
    parser.add_argument('-i',       dest="initvsini",        action="store",
                        help="Initial vsini (float, km/s). If no literature value known, use the value given by Step 2",
                        type=str,   default='' )
    parser.add_argument('-v',       dest="vsinivary",         action="store",
                        help="Range of allowed vsini variation during optimization, default = 5.0 km/s. Should be set to 0 for final run.",
                        type=str, default='5.0' )
    parser.add_argument('-g',       dest="guesses",           action="store",
                        help="For STD star. Initial RV guess for all nights. Given by Step 2 results (float, km/s)",
                        type=str,   default='' )
    parser.add_argument('-gS',       dest="guesses_source",           action="store",
                        help="For TAR star. Source for list of initial RV guesses. 'init' = Initguesser_results_X = past Step 2 result OR 'rvre' = RV_results_X = past Step 3 result",
                        type=str, default='')
    parser.add_argument('-gX',       dest="guessesX",           action="store",
                        help="For TAR star. The number, X, under ./*targname/Initguesser_results_X or ./*targname/RV_results_X, that you wish to use. Prefix determined by -gS",
                        type=str, default='')
    parser.add_argument('-t',       dest="template",         action="store",
                        help="Stellar template. Pick from 'synthetic', 'PHOENIX', or 'livingston'. Default = 'synthetic'",
                        type=str,   default='synthetic' )
    parser.add_argument('-temp',      dest="temperature",           action="store",
                        help="The synthetic template temperature used, e.g., 5000",
                        type=str,   default='' )
    parser.add_argument('-logg',      dest="logg",           action="store",
                        help="The synthetic template logg used, e.g., 4.5",
                        type=str,   default='' )
    parser.add_argument('-abs_out',   dest="abs",            action="store",
                        help="Take REL and ABS. REL for relative RVs as output, ABS for absolute RVs. Default = REL. Note that ABS mode will have worser precision.",
                        type=str,   default='REL' )
    # parser.add_argument('-sp',      dest="sptype",           action="store",
    #                     help="The spectral type of the *target. (Letter only)",
    #                     type=str,   default='' )
    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
                        help="If set, will generate plots of the fitting results under ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots under ./Temp/Debug/*target_*band/main_step3")
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
                        help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    parser.add_argument('--version',                          action='version',  version='%(prog)s 0.9')
    args = parser.parse_args()
    inpath   = '../Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    # Check user input

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)

    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'../Engine/syn_template')
    syntemp = [i for i in syntemp if i[:3] == 'syn'] #list of all syntheticstellar

    synT    = [ i.split('_')[2][1:]  for i in syntemp ]
    synlogg = [ i.split('_')[3][4:7] for i in syntemp ]

    if args.temperature not in synT:
        sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp" INPUT! {syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    if args.logg not in synlogg:
        sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg" INPUT! {syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    #------------------------------

    if (args.mode.lower() == 'std') & (args.guesses_source != ''):
        sys.exit('ERROR: STD CANNOT USE -gS, PLEASE USE -g')
    if (args.mode.lower() == 'tar') & (args.guesses != ''):
        sys.exit('ERROR: TAR CANNOT USE -g, PLEASE USE -gS')

    #------------------------------

    if args.nAB == '' and args.mode.lower() == 'std':
        nAB = 2
    elif args.nAB == '' and args.mode.lower() == 'tar':
        nAB = 3
    else:
        nAB = int(args.nAB)

    #------------------------------

    if args.abs.lower() not in ['rel', 'abs']:
        sys.exit('ERROR: UNEXPECTED INPUT FOR -abs_out')
    if args.abs.lower() == 'rel' and vsinivary != 0:
        sys.exit('ERROR: -abs_out must be set to "abs" until -v is set to 0!')
    if args.abs.lower() == 'rel':
        print_abs = 'Relative RV'
    else:
        print_abs = 'Absolute RV'


    #------------------------------

    if args.mode.lower() == 'std': # Specify initial RV guesses as a single value applied to all nights
        initguesses = np.float(args.guesses)
        initguesses_show = initguesses
    else: # Load initial RV guesses from file
        if args.guesses_source == 'init': # From Step 2 results
            guesses = '../Output/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
                                                                         args.band,
                                                                         int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])

        elif args.guesses_source == 'rvre': # From Step 3 results
            guesses = '../Output/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
                                                                      args.band,
                                                                      int(args.guessesX))
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs    = np.array(guessdata['RVfinal'])
            initguesses = {}
            initguesses_show = f'RVresultsSummary_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = np.float(initrvs[hrt])

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Tartget             =  {}
    Filter              = \33[41m {} band \033[0m
    WaveLength file     = \33[41m WaveRegions_{} \033[0m
    S/N cut             > \33[41m {} \033[0m
    Minium # of AB sets = \33[41m {} \033[0m             <------- If TAR mode, this should be at least 3. If STD mode, at least 2.
    Initial vsini       = \33[41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[41m {} km/s \033[0m
    RV initial guess    = \33[41m {} \033[0m
    Stellar template use= \33[41m {} \033[0m
    syn template temp   = \33[41m {} \033[0m
    syn template logg   = \33[41m {} \033[0m
    RV Output format    = \33[41m {} \033[0m
    Threads use         = {}
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut, nAB,
               initvsini, vsinivary, initguesses_show, args.template, args.temperature, args.logg, print_abs, args.Nthreads))
    if not args.skip:
        while True:
            inpp = input("Press [Y]es to continue, [N]o to quite...\n --> ")
            if 'n' in inpp.lower():
                sys.exit('QUIT, PLEASE RE-ENTER YOUR PARAMETERS')
            elif 'y' in inpp.lower():
                break
            else:
                print('I cannot understand what you are saying... TRY AGAIN')
                continue

    print('---------------------------------------------------------------')
    print('Running Step 3 for {}...'.format(args.targname))
    print('This will take a while..........')

    #-------------------------------------------------------------------------------

    if not os.path.isdir('../Output'):
        os.mkdir('../Output')

    if not os.path.isdir(f'../Output/{args.targname}_{args.band}_tool'):
        os.mkdir(f'../Output/{args.targname}_{args.band}_tool')
    filesndirs = os.listdir(f'../Output/{args.targname}_{args.band}_tool')

    trk = 1; go = True;
    while go == True:
        name = f'RV_results_{trk}'
        if name not in filesndirs:
            break
        trk += 1

    os.mkdir(f'../Output/{args.targname}_{args.band}_tool/{name}')

    outpath = f'../Output/{args.targname}_{args.band}_tool'
#-------------------------------------------------------------------------------
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}_tool.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
#-------------------------------------------------------------------------------
    logger.info(f'Writing output to ../Output/{args.targname}_{args.band}_tool/{name}')

    #-------------------------------------------------------------------------------

    # Read in the Prepdata under ../Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ../Input_Data/{}'.format(nnn, args.targname))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))

    logger.info('Analyze with {} nights'.format(len(nightsFinal)))

    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm,satm, mwave0, mflux0 = setup_templates(logger, args.template, args.band, np.int(args.temperature), np.float(args.logg))

    # Save pars in class for future use
    inparam = inparams(inpath,outpath,initvsini,vsinivary,args.plotfigs,
                       initguesses,bvcs,tagsA,tagsB,nightsFinal,mwave0,mflux0,None,xbounddict,maskdict)

    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when it was tight (T).
    # All statistical analysis will be performed separately for these two datasets.
    nights    = inparam.nights

    print('For paper plot!')
    if args.band == 'K':
        orders = np.array([3, 4, 5, 6, 8, 10])
    elif args.band=='H':
        orders = np.array([6, 13, 14, 16, 21, 22])
    # orders = np.array([6])
    #-------------------------------------------------------------------------------
    step2or3 = 3
    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        pool = mp.Pool(processes = args.Nthreads)
        func = partial(rv_MPinst, args, inparam, orders, jerp, trk, step2or3 )
        outs = pool.map(func, np.arange(len(nightsFinal)))
        pool.close()
        pool.join()

    logger.info('Output saved under {}/{}'.format(inparam.outpath, name) )
    print('####################################################################################')
