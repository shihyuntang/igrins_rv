from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step3

from Engine.IO_AB      import setup_templates, init_fitsread,stellarmodel_setup, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import A0cont
from Engine.classes    import fitobjs,inparams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod, fmod_conti
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
from Engine.crmask    import CRmasker

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def setup_fitting_init_pars(band, initvsini, order):

    # start at bucket loc = 1250 +- 100, width = 250 +- 100, depth = 100 +- 5000 but floor at 0
    centerloc = 1250 if band == 'H' else 1180
    
    # Initialize parameter array for optimization as well as half-range values for each parameter during
    # the various steps of the optimization.
    # Many of the parameters initialized here will be changed throughout the code before optimization and
    # in between optimization steps.
    pars0 = np.array([np.nan,                   # 0: The shift of the stellar template (km/s) [assigned later]
                      0.3,                      # 1: The scale factor for the stellar template
                      0.0,                      # 2: The shift of the telluric template (km/s)
                      0.6,                      # 3: The scale factor for the telluric template
                      initvsini,                # 4: vsini (km/s)
                      np.nan,                   # 5: The instrumental resolution (FWHM) in pixels
                      0.0,                      # 6: Wavelength 0-pt
                      0.0,                      # 7: Wavelength linear component
                      0.0,                      # 8: Wavelength quadratic component
                      0.0,                      # 9: Wavelength cubic component
                      1.0,                      #10: Continuum zero point
                      0.,                       #11: Continuum linear component
                      0.,                       #12: Continuum quadratic component
                      np.nan,                   #13: Instrumental resolution linear component
                      np.nan,                   #14: Instrumental resolution quadratic component
                      centerloc,                #15: Blaze dip center location
                      330,                      #16: Blaze dip full width
                      0.05,                     #17: Blaze dip depth
                      90,                       #18: Secondary blaze dip full width
                      0.05,                     #19: Blaze dip depth
                      0.0,                      #20: Continuum cubic component
                      0.0,                      #21: Continuum quartic component
                      0.0,                      #22: Continuum pentic component
                      0.0])                     #23: Continuum hexic component
    
    if int(order) == 13: pars0[1] = 0.8
    
    return pars0


def base_dpars_dict(vsinivary, masterbeam, band, order):
    """setup basic sets of par vary range array

    Parameters
    ----------
    use_sets : list with str
        list of dpars_org keys that wish to get

    Returns
    -------
    Dict
        dpars
    """
    
    #                                | 0    1    2    3 |  | -- 4 -- | | 5 |   | 6     7     8           9  |  |10  11  12| |13 14|  |15   16   17   18    19 |  |20   21   22   23 |
    dpars_org = {'cont' : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,        0.0,    0.0,  0.0,  0.0,        0.0,    1e7, 1, 1,   0, 0,    10., 30., 0.2, 50.0, 0.2,   1.0, 1.0, 1.0, 1.0 ]),
                 'twave': np.array([  0.0, 0.0, 0.0, 1.0,   0.0,        0.0,    0.1,  0.1,  0.1,        0.1,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'ip'   : np.array([  0.0, 0.0, 0.0, 0.0,   0.0,        0.5,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 's'    : np.array([  5.0, 1.0, 0.0, 0.0,   0.0,        0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'v'    : np.array([  0.0, 0.0, 0.0, 0.0,   vsinivary,  0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ]),
                 'ts'   : np.array([  5.0, 1.0, 0.0, 1.0,   0.0,        0.0,    0.0,  0.0,  0.0,        0.0,    0,   0, 0,   0, 0,     0.,  0., 0.0,  0.,  0.0,   0.0, 0.0, 0.0, 0.0 ])}
    
    if masterbeam == 'B':
        dpars_org['cont'][15:20] = np.array([0.0, 0.0, 0.0, 0.0, 0.0]) # setting for the Blaze dip
        
    # blaze fitting order setting
    if band == 'H':
        if order in [13]:
            # fit a quadratic (2) continuum
            dpars_org['cont'][20:] = 0
        elif order in [6,14,21]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        else:
            pass
    elif band == 'K':
        if order in [3,4,5]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        elif order in [6]:
            # fit a quartic (4) continuum
            dpars_org['cont'][22:] = 0
        else:
            pass
    
    return dpars_org





def rv_MPinst(args, inparam, orders, order_use, trk, step2or3, i):

    # Main function for RV fitting that will be threaded over by multiprocessing

    nights   = inparam.nights
    night = nights[i] # current looped night

    order = orders[order_use]
    xbounds = inparam.xbounddict[order]

    if args.debug:
        print('Working on order {:02d}/{:02d} ({}), night {:03d}/{:03d} ({}) PID:{}...'.format(int(order_use)+1,
                                                                                               len(orders),
                                                                                               order,
                                                                                               i+1,
                                                                                               len(inparam.nights),
                                                                                               night,
                                                                                               mp.current_process().pid) )

    #-------------------------------------------------------------------------------

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = []
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    nightsout = []
    rvsminibox     = np.ones(len(tagsnight))
    vsiniminibox   = np.ones(len(tagsnight))
    tagsminibox    = np.ones(len(tagsnight))
    chisminibox    = np.ones(len(tagsnight))
    parfitminibox  = np.ones((len(tagsnight),24)) # need to match the dpar numbers

    rvsminibox[:]    = np.nan
    vsiniminibox[:]  = np.nan
    tagsminibox[:]   = np.nan
    parfitminibox[:] = np.nan
    chisminibox[:] = np.nan

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

    pars0 = setup_fitting_init_pars(args.band, inparam.initvsini, order)
    

    # Iterate over all A/B exposures
    for t in np.arange(len(tagsnight)):
        tag = tagsnight[t]
        beam = beamsnight[t]
        masterbeam = beam

        if np.int(night[:8]) == 20170216 and args.targname == 'GJ281' and np.float(tag) == 63:
            continue

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

        A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

        try:
            hdulist = fits.open(A0loc)
        except IOError:
            logger.warning(f'  --> No A0-fitted template for night {night}, skipping...')
            return nightsout, rvsminibox, parfitminibox, vsiniminibox, tagsminibox

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
            return nightsout, rvsminibox, parfitminibox, vsiniminibox, tagsminibox


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


        bound_cut = _setup_bound_cut(inparam.bound_cut_dic, args.band, order)
        
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

        # Save data for second template cutting after optimization cycle 1 done
        s_save = s_piece.copy(); x_save = x_piece.copy(); u_save = u_piece.copy();

        # Trim telluric template to data range +- 15 AA. If telluric template buffer is cut short because A0 lines didn't extend 
        # far past data range, cut data range accordingly.
        satm_in = satm[(watm > np.min(wave_piece)*1e4 - 10) & (watm < np.max(wave_piece)*1e4 + 10)]
        watm_in = watm[(watm > np.min(wave_piece)*1e4 - 10) & (watm < np.max(wave_piece)*1e4 + 10)]

        s_piece	= s_piece[   (wave_piece*1e4 > np.min(watm_in)+10) & (wave_piece*1e4 < np.max(watm_in)-10)]
        u_piece	= u_piece[   (wave_piece*1e4 > np.min(watm_in)+10) & (wave_piece*1e4 < np.max(watm_in)-10)]
        x_piece	= x_piece[   (wave_piece*1e4 > np.min(watm_in)+10) & (wave_piece*1e4 < np.max(watm_in)-10)]
        wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+10) & (wave_piece*1e4 < np.max(watm_in)-10)]

        # Trim stellar template to data range +- 10 AA
        mflux_in = inparam.mflux0[(inparam.mwave0 > np.min(wave_piece)*1e4 - 5) & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5)]
        mwave_in = inparam.mwave0[(inparam.mwave0 > np.min(wave_piece)*1e4 - 5) & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5)]

        Rstell = np.median(np.diff(mwave_in)); Rtell = np.median(np.diff(watm_in));
        if Rstell < Rtell:
            sys.exit(f'Telluric template resolution ({round(Rtell,4)} AA) must be finer than stellar template resolution ({round(Rstell,4)} AA) !')

        # Rebin stellar template to uniform wavelength scale
        dstep = Rstell
        nstep = int((mwave_in[-1]-mwave_in[0])/dstep)
        mwave1 = np.linspace(mwave_in[0],mwave_in[-1],nstep)
        mflux1 = rebin_jv(mwave_in,mflux_in,mwave1,False)
        mwave_in = mwave1.copy(); mflux_in = mflux1.copy();
        mwave_in = mwave_in[1:-1]
        mflux_in = mflux_in[1:-1]
	
        # Normalize continuum from A0 to flux scale of data
        continuum /= np.nanmedian(continuum)
        continuum *= np.nanpercentile(s_piece,99)

        # --------------------------------------------------------------

        par = pars0.copy()

        # Get initial guess for cubic wavelength solution from reduction pipeline
        f = np.polyfit(x_piece,wave_piece,3)
        q = np.poly1d(f)
        initwave = q(x_piece)*1e4
        
        par[0] = initguesses-inparam.bvcs[night+tag] # Initial RV with barycentric correction
        par[5]  = IPpars[2]
        par[13] = IPpars[1]
        par[14] = IPpars[0]
        
        
        # setup fitting boundary
        dpars = base_dpars_dict(inparam.vsinivary, masterbeam, args.band, int(order))
        
        continuum_in = rebin_jv(a0contx,continuum,x_piece,False)
        fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]),masterbeam,[np.array([],dtype=int),np.array([],dtype=int)], initwave)

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

        # added new twave step at start
        optgroup = ['twave', 'cont', 'twave', 'cont', 'ts',
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

            if nc == cycles:
                optgroup = optgroup[:-1] # if is the last run, skip the 'ts'
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
                CRmaskF = CRmasker(parfit,fitobj)

                fitobj = fitobjs(s_piece, x_piece, u_piece, continuum_in, watm_in,satm_in,mflux_in,mwave_in,ast.literal_eval(inparam.maskdict[order]),masterbeam,CRmaskF, initwave)

        parfit = parfit_1.copy()

        #-------------------------------------------------------------------------------

        # if best fit stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            logger.warning(f'  --> Stellar template power is low for {night}! Data likely being misfit! Throwing out result...')
            continue

        # if best fit stellar or telluric template powers are exactly equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3]:
            logger.warning(f'  --> Stellar or telluric template powers have not budged from starting values for {night}! Fit is broken! Optimizer bounds may be unfeasible, or chi-squared may be NaN? Throwing out result...')
            continue

        # if best fit model dips below zero at any point, we're to close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq = fmod(parfit,fitobj)
        if len(smod[(smod < 0)]) > 0:
            logger.warning(f'  --> Best fit model dips below 0 for {night}! May be too close to edge of blaze, throwing out result...')
            continue


        #-------------------------------------------------------------------------------

        if args.plotfigs == True:
            parfitS = parfit.copy(); parfitS[3] = 0
            parfitT = parfit.copy(); parfitT[1] = 0
            outplotter_23(parfitS, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3,order)
            outplotter_23(parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3,order)
            outplotter_23(parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag), trk, inparam, args, step2or3,order)

        rv0 = parfit[0]

        rvsminibox[t]   = rv0  + inparam.bvcs[night+tag] + rv0*inparam.bvcs[night+tag]/(2.99792458e5**2) # Barycentric correction
        parfitminibox[t]= parfit
        vsiniminibox[t] = parfit[4]
        tagsminibox[t]  = tag
        fit,chi = fmod(parfit, fitobj)
        chisminibox[t]  = chi

    # If any tag has a chisq value 10x greater than the best fit, flag it as a misfit and throw out its result
    for t in np.arange(len(tagsnight)):
        if np.isnan(chisminibox[t]) == False:
            if chisminibox[t] > 10.*np.nanmin(chisminibox):
                logger.warning(f'  --> Chi-squared indicates a misfit for observation {night} {order} {tagsnight[t]}')
                rvsminibox[t]   = np.nan
                vsiniminibox[t] = np.nan

    return nightsout,rvsminibox,parfitminibox,vsiniminibox,tagsminibox

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    # parser = argparse.ArgumentParser(
    #                                  prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 3',
    #                                  description = '''
    #                                  Performs a full analysis of each target star observation to produce accurate and precise RVs. \n
    #                                  All the wavelength regions defined in Step 1 are used, and the code analyzes each observation that is part of a given exposure separately. \n
    #                                  Unless the target vsini is already known to high accuracy, an initial run of Step 3 in which \vsini is allowed to vary is required. \n
    #                                  This provides an estimate of vsini that can then be plugged into the code as a fixed value in the second run of Step 3. \n
    #                                  If the user seeks the best possible RV uncertainty estimates, or if their target star has a relatively high \vsini ($>$ 10 \kms), they must run Step 3 once with \vsini held fixed at its estimated value and once with \vsini held fixed at this value plus or minus one sigma. \n
    #                                  The minor differences in the RVs of the two runs (as low as $<$1 \ms and as high as 7 \ms) can then be incorporated into the final uncertainties. \n
    #                                  If \vsini is already well-known, it is not necessary to run Step 3 more than once, as the code fully converges to the final RVs (within uncertainty) through just one run.
    #                                  ''',
    #                                  epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    # parser.add_argument("targname",                          action="store",
    #                     help="Enter your *target name",            type=str)
    # parser.add_argument("-mode",    dest="mode",             action="store",
    #                     help="RV standard star (STD) or a normal target (TAR)?",
    #                     type=str,   default='')
    # parser.add_argument("-HorK",    dest="band",             action="store",
    #                     help="Which band to process? H or K?. Default = K",
    #                     type=str,   default='K')
    # parser.add_argument("-Wr",      dest="WRegion",          action="store",
    #                     help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
    #                     type=int,   default=int(1))
    # parser.add_argument("-SN",      dest="SN_cut",           action="store",
    #                     help="Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 ",
    #                     type=str,   default='50')
    # parser.add_argument("-nAB",      dest="nAB",           action="store",
    #                     help="Minium number of separte A/B exposures within a set for a given observation (ensures accuracy of uncertainy estimates). Default = 2 for STD, 3 for TAR",
    #                     type=str,   default='')
    # parser.add_argument('-i',       dest="initvsini",        action="store",
    #                     help="Initial vsini (float, km/s). If no literature value known, use the value given by Step 2",
    #                     type=str,   default='' )
    # parser.add_argument('-v',       dest="vsinivary",         action="store",
    #                     help="Range of allowed vsini variation during optimization, default = 5.0 km/s. Should be set to 0 for final run.",
    #                     type=str, default='5.0' )
    # parser.add_argument('-g',       dest="guesses",           action="store",
    #                     help="For STD star. Initial RV guess for all nights. Given by Step 2 results (float, km/s)",
    #                     type=str,   default='' )
    # parser.add_argument('-gS',       dest="guesses_source",           action="store",
    #                     help="For TAR star. Source for list of initial RV guesses. 'init' = Initguesser_results_X = past Step 2 result OR 'rvre' = RV_results_X = past Step 3 result",
    #                     type=str, default='')
    # parser.add_argument('-gX',       dest="guessesX",           action="store",
    #                     help="For TAR star. The number, X, under ./*targname/Initguesser_results_X or ./*targname/RV_results_X, that you wish to use. Prefix determined by -gS",
    #                     type=str, default='')
    # parser.add_argument('-t',       dest="template",         action="store",
    #                     help="Stellar template. Pick from 'synthetic', 'PHOENIX', or 'livingston'. Default = 'synthetic'",
    #                     type=str,   default='synthetic' )
    # parser.add_argument('-temp',      dest="temperature",           action="store",
    #                     help="The synthetic template temperature used, e.g., 5000",
    #                     type=str,   default='' )
    # parser.add_argument('-logg',      dest="logg",           action="store",
    #                     help="The synthetic template logg used, e.g., 4.5",
    #                     type=str,   default='' )
    # parser.add_argument('-abs_out',   dest="abs",            action="store",
    #                     help="Take REL and ABS. REL for relative RVs as output, ABS for absolute RVs. Default = REL. Note that ABS mode will have worser precision.",
    #                     type=str,   default='REL' )
    # # parser.add_argument('-sp',      dest="sptype",           action="store",
    # #                     help="The spectral type of the *target. (Letter only)",
    # #                     type=str,   default='' )
    # parser.add_argument('-c',       dest="Nthreads",         action="store",
    #                     help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
    #                     type=int,   default=int(mp.cpu_count()//2) )
    # parser.add_argument('-plot',    dest="plotfigs",          action="store_true",
    #                     help="If set, will generate plots of the fitting results under ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    # parser.add_argument('-n_use',   dest="nights_use",       action="store",
    #                     help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
    #                     type=str,   default='')
    # parser.add_argument('-DeBug',    dest="debug",           action="store_true",
    #                     help="If set, DeBug logging will be output, as well as (lots of) extra plots.")
    # parser.add_argument('-sk_check', dest="skip",           action="store_true",
    #                     help="If set, will skip the input parameters check. Handy when running mutiple targets line by line")
    # parser.add_argument('--version',                          action='version',  version='%(prog)s 1.0.0')
    # args = parser.parse_args()
    
    args = _argparse_step3()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    #-------------------------------------------------------------------------------

    # Check user input

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)

    if args.mode == '':
        sys.exit('ERROR: YOU MUST CHOOSE A MODE, "STD" OR "TAR", for "-mode"')

    if args.initvsini == '':
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI VALUE, "-i"')

    if (args.guesses == '') & (args.guessesX == ''):
        sys.exit('ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR RV VALUE(S) BY USING "-g" OR "-gX"')

    if (args.temperature == '') & (args.logg == ''):
        sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE AVAILABLE TEMPLATES')
    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'./Engine/syn_template')

    if args.template.lower() == 'synthetic':
        syntemp = [i for i in syntemp if i[:3] == 'syn'] #list of all syntheticstellar
        synT    = [ i.split('_')[2][1:]  for i in syntemp ]
        synlogg = [ i.split('_')[3][4:7] for i in syntemp ]
    elif args.template.lower() == 'phoenix':
        syntemp = [i for i in syntemp if i[:3] == 'PHO'] #list of all phoenix
        synT    = [ i.split('-')[1][4:]  for i in syntemp ]
        synlogg = [ i.split('-')[2][:3] for i in syntemp ]
    else:
        synT = [args.temperature]; synlogg = [args.logg];
	
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
            guesses = './Output/{}_{}/Initguesser_results_{}.csv'.format(args.targname,
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
            guesses = './Output/{}_{}/RVresultsSummary_{}.csv'.format(args.targname,
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
    Filter              = \33[37;1;41m {} band \033[0m
    WaveLength file     = \33[37;1;41m WaveRegions_{} \033[0m
    S/N cut             > \33[37;1;41m {} \033[0m
    Minium # of AB sets = \33[37;1;41m {} \033[0m             <------- If TAR mode, this should be at least 3. If STD mode, at least 2.
    Initial vsini       = \33[37;1;41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[37;1;41m {} km/s \033[0m
    RV initial guess    = \33[37;1;41m {} \033[0m
    Stellar template use= \33[37;1;41m {} \033[0m
    syn template temp   = \33[37;1;41m {} \033[0m
    syn template logg   = \33[37;1;41m {} \033[0m
    RV Output format    = \33[37;1;41m {} \033[0m
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

    # Make output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')
    filesndirs = os.listdir(f'./Output/{args.targname}_{args.band}')

    trk = 1; go = True;
    while go == True:
        name = f'RV_results_{trk}'
        if name not in filesndirs:
            break
        trk += 1

    os.mkdir(f'./Output/{args.targname}_{args.band}/{name}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs')

    step2or3 = 3
    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs/main_step{step2or3}_{args.band}_{trk}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs/main_step{step2or3}_{args.band}_{trk}')

    outpath = f'./Output/{args.targname}_{args.band}'

    #-------------------------------------------------------------------------------

    # Set up logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    logger.info(f'Writing output to ./Output/{args.targname}_{args.band}/{name}')

    #-------------------------------------------------------------------------------

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit('NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(nnn, args.targname))
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
    intnights = np.array([int(i[:8]) for i in nights])

    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these nights
                slightly differently. When you run Step 3, RVs will be output in
                two formats: one with the defocus nights separated, and the other
                with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as the result
                of the defocus. IGRINS RV will process defocus nights the same way
                as the others, but when you run Step 3, will still output the results
                in two formats like it does with the K band.''')

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]
    obsT    = obs[indT]
    obsL    = obs[indL]
    rvmasterboxT  = np.ones((len(nightsT),len(orders)))*np.nan
    stdmasterboxT = np.ones((len(nightsT),len(orders)))*np.nan
    rvmasterboxL  = np.ones((len(nightsL),len(orders)))*np.nan
    stdmasterboxL = np.ones((len(nightsL),len(orders)))*np.nan
    vsinisT       = np.ones((len(nightsT),len(orders)))*np.nan
    vsinisL       = np.ones((len(nightsL),len(orders)))*np.nan

    if len(nightsL) > 0 and len(nightsT) > 0:
        nightscomblist = [nightsT,nightsL]
        T_Ls = ['T','L']
    elif len(nightsL) > 0:
        nightscomblist = [nightsL]
        T_Ls = ['L']
    else:
        nightscomblist = [nightsT]
        T_Ls = ['T']

    # orders = np.array([6])
    # print('ONLY process order 6')
    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        if not args.debug: print('Working on order {} ({:02d}/{:02d})'.format(orders[jerp], int(jerp+1), len(orders)))

        #rv_MPinst( args, inparam, orders, jerp, trk, step2or3 , 0)
        func = partial(rv_MPinst, args, inparam, orders, jerp, trk, step2or3 )
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)
        # pool = mp.Pool(processes = args.Nthreads)
        # outs = pool.map(func, np.arange(len(nightsFinal)))
        # pool.close()
        # pool.join()

        order = orders[jerp]

        #-------------------------------------------------------------------------------

        # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
        for i in range(len(nights)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvbox     = outsbox[1]
                parfitbox = outsbox[2]
                vsinibox  = outsbox[3]
                tagbox    = outsbox[4]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvbox     = np.concatenate((rvbox,outsbox[1]))
                parfitbox = np.vstack((parfitbox,outsbox[2]))
                vsinibox  = np.concatenate((vsinibox,outsbox[3]))
                tagbox    = np.concatenate((tagbox,outsbox[4]))

        nightsbox = np.array(nightsbox)
        vsinitags = []

        # Save results to fits file
        c1    = fits.Column(name='NIGHT'+str(order),  array=nightsbox, format='{}A'.format(len(nights[0])) )
        c2    = fits.Column(name='RV'+str(order),     array=rvbox,     format='D')
        c3    = fits.Column(name='PARFIT'+str(order), array=parfitbox, format=str(len(parfitbox[0,:]))+'D', dim=(1,len(parfitbox[0,:])))
        c4    = fits.Column(name='VSINI'+str(order),  array=vsinibox,  format='D')
        c5    = fits.Column(name='TAG'+str(order),    array=tagbox,    format='4A')
        cols  = fits.ColDefs([c1,c2,c3,c4,c5])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if jerp == 0: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
        else:
            hh = fits.open('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
            hh.append(hdu_1)
            hh.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name),overwrite=True)

        #-------------------------------------------------------------------------------

        # For each set of nights (tight, loose)...
        iT_L = 0
        for nights_use in nightscomblist:

            T_L = T_Ls[iT_L]

            # For each night...
            for i in range(len(nights_use)):
                # Collect the RVs and vsinis determined from different A/B exposures within a night
                indnight  = np.where(nightsbox == nights_use[i])[0]
                rvtags    = rvbox[indnight]
                vsinitags = vsinibox[indnight]

                # Take the mean of the vsinis, and the mean and std of the RVs.
                # If the number of different successfully fitted A/B exposures is less than required, pass NaN instead.
                if T_L == 'T':
                    vsinisT[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxT[i,jerp]  = np.nan
                        stdmasterboxT[i,jerp] = np.nan
                    else:
                        rvmasterboxT[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxT[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags[~np.isnan(rvtags)]))

                else:
                    if order == 3 and args.band == 'K': # Don't use vsini estimates from this order during loose epoch
                        vsinisL[i,jerp] = np.nan
                    else:
                        vsinisL[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxL[i,jerp]  = np.nan
                        stdmasterboxL[i,jerp] = np.nan
                    else:
                        rvmasterboxL[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxL[i,jerp] = np.nanstd(rvtags)/np.sqrt(len(rvtags[~np.isnan(rvtags)]))
            iT_L += 1


    #-------------------------------------------------------------------------------
    if not args.debug: logger.addHandler(stream_hander)
    print('\n')

    # Don't combine Loose and Tight datasets, but make them both easily referenceable
    nightsCombined  = np.array([]); jdsCombined = np.array([]);
    rvfinalCombined = np.array([]); stdfinalCombined = np.array([]); vsinifinalCombined = np.array([]);

    if T_Ls == ['T','L']:
        rvboxcomblist  = [rvmasterboxT,rvmasterboxL]
        stdboxcomblist = [stdmasterboxT,stdmasterboxL]
        vsinicomblist  = [vsinisT,vsinisL]
        obscomblist    = [obsT,obsL]
    elif T_Ls == ['L']:
        rvboxcomblist  = [rvmasterboxL]
        stdboxcomblist = [stdmasterboxL]
        vsinicomblist  = [vsinisL]
        obscomblist    = [obsL]
    else:
        rvboxcomblist  = [rvmasterboxT]
        stdboxcomblist = [stdmasterboxT]
        vsinicomblist  = [vsinisT]
        obscomblist    = [obsT]

    # Avg zero-point offset of the different observatory and mounting epoch RVs, to add back in later
    if args.abs.lower() == 'rel':
        rvZeroBase = np.nanmean([np.nanmean(rvboxcomblist[0][(obscomblist[0] == np.unique(obscomblist[0])[0]),ll]) for ll in range(len(orders))])

    # Iterate over tight and loose mounting data sets...
    for boxind in range(len(rvboxcomblist)):

        rvmasterbox  = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]
        vsinibox     = vsinicomblist[boxind]
        obsbox       = obscomblist[boxind]

        #-------------------------------------------------------------------------------

        if args.mode=='STD': # If RV STD star, calculate uncertainty in method

            # Calculate the precision within an order across nights
            sigma_O2     = np.array([np.nanstd(rvmasterbox[:,ll])**2 for ll in range(len(orders))])
            sigma_ABbar2 = np.ones_like(sigma_O2)

            # Calculate uncertainty in method as difference between variance within an order and mean variance within a night's different exposure RVs
            for ll in range(len(orders)):
                sigma_ABbar2[ll] = np.nanmedian(stdmasterbox[:,ll]**2)
            sigma_method2 = sigma_O2 - sigma_ABbar2

        else: # If target star, load the uncertainty in method calculated from our RV STD star runs
            if T_Ls[boxind] == 'T':
                nights_use = nightsT.copy()
                kind = 'Focused'
                sigma_method2 = inparam.methodvariance_tight[args.band]
            elif T_Ls[boxind] == 'L':
                nights_use = nightsL.copy()
                kind = 'Defocus'
                sigma_method2 = inparam.methodvariance_loose[args.band]

        sigma_ON2    = np.ones_like(rvmasterbox)

        #-------------------------------------------------------------------------------

        # Note rvmasterbox indexed as [nights,orders]
        Nnights = len(rvmasterbox[:,0])

        if args.abs.lower() == 'rel':
            rvBase = {}
            for obs_name in np.unique(obsbox):
                #rvBase[obs_name] = np.nanmean([np.nanmean(rvmasterbox[(obsbox == obs_name),ll]) for ll in range(len(orders))])
                rvBase[obs_name] = np.nanmean([np.nansum( ( (1./(stdmasterbox[(obsbox == obs_name),ll]**2)) / (np.nansum(1./(stdmasterbox[(obsbox == obs_name),ll]**2))) )*rvmasterbox[(obsbox == obs_name),ll]) for ll in range(len(orders))])
        
        for ll in range(len(orders)):

	    # Mean-subtract each order's RVs within an observatory epoch
            if args.abs.lower() == 'rel':
                for obs_name in np.unique(obsbox):
                    #rvmasterbox[(obsbox == obs_name),ll] -= np.nanmean(rvmasterbox[(obsbox == obs_name),ll]) - (rvZeroBase - rvBase[obs_name])
                    meanweights = (1./(stdmasterbox[(obsbox == obs_name),ll]**2)) / (np.nansum(1./(stdmasterbox[(obsbox == obs_name),ll]**2))) 
                    rvmasterbox[(obsbox == obs_name),ll] -= np.nansum( meanweights*rvmasterbox[(obsbox == obs_name),ll]) - (rvZeroBase - rvBase[obs_name])
            
            for rvin in rvmasterbox[(obsbox == 'NA'),ll]:
                if np.isnan(rvin) == False:
                    sys.exit('File listed as NA for A0 presence output a RV! Your Prepdata files are not the same version as your A0Fit files!')

            # Calculate the uncertainty in each night/order RV as the sum of the uncertainty in method and the uncertainty in that night's As and Bs RVs
            for night in range(Nnights):
                sigma_ON2[night,ll] = sigma_method2[ll] + stdmasterbox[night,ll]**2

        # If taking absolute RVs, add in uncertainty characterized by scatter between mean RVs of different orders
        if args.abs.lower() == 'abs':
            sigma_order_to_order = np.nanstd([np.nanmean(rvmasterbox[:,ll]) for ll in range(len(orders))])/np.sqrt(len(orders))

        rvfinal    = np.ones(Nnights, dtype=np.float64)
        stdfinal   = np.ones(Nnights, dtype=np.float64)
        vsinifinal = np.ones(Nnights, dtype=np.float64)
        jds_out   = np.ones(Nnights, dtype=np.float64)

        if T_Ls[boxind] == 'T':
            nights_use = nightsT.copy(); kind = 'Focused';
        elif T_Ls[boxind] == 'L':
            nights_use = nightsL.copy(); kind = 'Defocused';


        # Combine RVs between orders using weights calculated from uncertainties
        for n in range(Nnights):
            weights = (1./sigma_ON2[n,:]) / (np.nansum(1./sigma_ON2[n,:])) # normalized
            stdspre = (1./sigma_ON2[n,:]) #unnormalized weights

            rvfinal[n]  = np.nansum( weights*rvmasterbox[n,:] )
            stdfinal[n] = 1/np.sqrt(np.nansum(stdspre))

            vsinifinal[n] = np.nansum(weights*vsinibox[n,:])
            jds_out[n]   = jds[nights_use[n]]

            # if all the RVs going into the observation's final RV calculation were NaN due to any pevious errors, pass NaN
            if np.nansum(weights) == 0:
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan

            # if more than half of the orders going into the observation's final RV calculation were NaN due to any pevious errors, pass NaN
            if np.sum( np.isnan(rvmasterbox[n,:]) ) > np.floor( len(orders) * 0.5 ):
                rvfinal[n]    = np.nan
                stdfinal[n]   = np.nan
                vsinifinal[n] = np.nan

        #-------------------------------------------------------------------------------

        if args.abs.lower() == 'abs':
            stdfinal = np.sqrt(stdfinal**2 + sigma_order_to_order**2)
        else:
            # Correct for zero-point offset between loose and tight epochs
            if T_Ls[boxind] == 'T':
                rvMeanTight = np.nanmean(rvfinal)

            elif (T_Ls[boxind] == 'L') & (len(T_Ls)==2):
                logger.info('Mean RV during the Defocus mounting period, before subtraction = {:1.4f} km/s'.format(np.nanmean(rvfinal)))
                logger.info('Mean RV during the Focused mounting period, before subtraction = {:1.4f} km/s'.format(rvMeanTight))
                logger.info('Value used to correct for this = {:1.4f} km/s'.format(np.nanmean(rvfinal) - rvMeanTight))
                rvfinal -= np.nanmean(rvfinal) - rvMeanTight

        # Plot results
        f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)

        axes.plot(    np.arange(len(rvfinal))+1, rvfinal, '.k', ms=5)
        axes.errorbar(np.arange(len(rvfinal))+1, rvfinal, yerr=stdfinal, ls='none', lw=.5, ecolor='black')
        axes.text(0.05, 0.93, r'RV mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinal), np.nanstd(rvfinal)),
                             transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        axes.set_ylim(np.nanmin(rvfinal)-.08,
                     np.nanmax(rvfinal)+.08)
        axes.set_ylabel('RV [km/s]', size=6, style='normal', family='sans-serif' )
        axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(5))
        axes.tick_params(axis='both', which='both', labelsize=5, right=True, top=True, direction='in', width=.6)
        f.savefig('{}/{}/FinalRVs_{}_.png'.format(inparam.outpath, name, kind), format='png', bbox_inches='tight')

        # Save results to fits file separately for each tight/loose dataset
        c1 = fits.Column( name='NIGHT',         array=nights_use,    format='8A')
        c2 = fits.Column( name='JD',            array=jds_out,      format='D')
        c3 = fits.Column( name='RVBOX',         array=rvmasterbox,   format='{}D'.format(len(orders)))
        c4 = fits.Column( name='STDBOX',        array=stdmasterbox,  format='{}D'.format(len(orders)))
        c7 = fits.Column( name='Sigma_method2', array=sigma_method2, format='D')
        c8 = fits.Column( name='Sigma_ON2',     array=sigma_ON2,     format='{}D'.format(len(orders)))
        c9 = fits.Column( name='RVfinal',       array=rvfinal,       format='D')
        c10 = fits.Column(name='STDfinal',      array=stdfinal,      format='D')
        c11 = fits.Column(name='VSINI',         array=vsinifinal,    format='D')

        if args.mode=='STD':
            c5 = fits.Column( name='Sigma_O2',      array=sigma_O2,      format='D')
            c6 = fits.Column( name='Sigma_ABbar2',  array=sigma_ABbar2,  format='D')
            cols  = fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11])
        else:
            cols  = fits.ColDefs([c1,c2,c3,c4,c7,c8,c9,c10,c11])

        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul        = fits.HDUList([primary_hdu,hdu_1])
        if len(T_Ls) == 2:
            hdul.writeto('{}/{}/RVresultsAdvanced_{}.fits'.format(inparam.outpath, name, kind), overwrite=True)
        else:
            hdul.writeto('{}/{}/RVresultsAdvanced.fits'.format(inparam.outpath, name), overwrite=True)

        # Combine final RVs from both tight and loose mounting data sets
        nightsCombined     = np.concatenate((nightsCombined,     nights_use))
        jdsCombined       = np.concatenate((jdsCombined,       jds_out))
        rvfinalCombined    = np.concatenate((rvfinalCombined,    rvfinal))
        stdfinalCombined   = np.concatenate((stdfinalCombined,   stdfinal))
        vsinifinalCombined = np.concatenate((vsinifinalCombined, vsinifinal))

        if args.mode=='STD': # If uncertainty in method was calculated, save it
            sigma_method2 = [np.around(float(i), 8) for i in sigma_method2]
            logger.info('sigma_method2 during the {} epoch is {}'.format(kind, sigma_method2))
        if len(T_Ls) == 2:
            logger.info('During the {} epoch: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( kind,
                                                                                                        np.nanmean(rvfinal),
                                                                                                        np.nanstd(rvfinal) ))
            logger.info('During the {} epoch: vsini mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( kind,
                                                                                                        np.nanmean(vsinifinal),
                                                                                                        np.nanstd(vsinifinal) ))
        else:
            logger.info('RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( np.nanmean(rvfinal),
                                                                                                        np.nanstd(rvfinal) ))
            logger.info('vsini mean = {:1.4f} km/s, std = {:1.4f} km/s'.format( np.nanmean(vsinifinal),
                                                                                                        np.nanstd(vsinifinal) ))

    #-------------------------------------------------------------------------------

    # Plot combined results
    xscale = np.arange(len(rvfinalCombined))+1

    f, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(xscale,rvfinalCombined, '.k', ms=5)
    axes.errorbar(xscale,rvfinalCombined,yerr=stdfinalCombined,ls='none',lw=.5, ecolor='black')
    axes.text(0.05, 0.93, r'RV mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(np.nanmean(rvfinalCombined), np.nanstd(rvfinalCombined)),
                         transform=axes.transAxes, size=6, style='normal', family='sans-serif' )

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        axes.text(0.05, 0.1, 'Defocus', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            axes.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
        else:
            axes.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
            axes.text(0.9,  0.1, 'Defocused', transform=axes.transAxes, size=6, style='normal', family='sans-serif' )
    axes.set_ylim(np.nanmin(rvfinalCombined)-.08,np.nanmax(rvfinalCombined)+.08)
    axes.set_ylabel('RV (km/s)', size=6, style='normal', family='sans-serif' )
    axes.set_xlabel('Night (#)', size=6, style='normal', family='sans-serif' )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in', width=.6)
    f.savefig('{}/{}/FinalRVs.png'.format(inparam.outpath, name), format='png', bbox_inches='tight')

    # Output combined final results to fits file
    c1 = fits.Column(name='NIGHT',    array=nightsCombined,         format='{}A'.format(len(nights[0])) )
    c2 = fits.Column(name='JD',       array=jdsCombined,           format='D')
    c3 = fits.Column(name='RVfinal',  array=rvfinalCombined,        format='D')
    c4 = fits.Column(name='STDfinal', array=stdfinalCombined,       format='D')
    c5 = fits.Column(name='VSINI',    array=vsinifinalCombined,     format='D')

    cols = fits.ColDefs([c1,c2,c3,c4,c5])
    hdu_1 = fits.BinTableHDU.from_columns(cols)

    bleh = np.ones((3,3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu,hdu_1])
    hdul.writeto('{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name),overwrite=True)

    tempin = Table.read('{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name), format='fits')
    tempin.write('{}/RVresultsSummary_{}.csv'.format(inparam.outpath, trk), format='csv', overwrite=True)

    if len(T_Ls) == 2:
        logger.info('Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(rvfinalCombined),
                                                                                np.nanstd(rvfinalCombined)))
        logger.info('vsini results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(np.nanmean(vsinifinalCombined),
                                                                                np.nanstd(vsinifinalCombined)))

    warning_r = log_warning_id(f'{outpath}/{args.targname}_{args.band}.log', start_time)
    if warning_r:
        print(f'''
**********************************************************************************
WARNING!! you got warning message during this run. Please check the log file under:
          {outpath}/{args.targname}_{args.band}_A0Fits.log
**********************************************************************************
''')
    print('\n')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}/{}'.format(inparam.outpath, name) )
    logger.info('The final RV estimates you are looking for are in the RVresultsSummary files!')
    print('####################################################################################')
