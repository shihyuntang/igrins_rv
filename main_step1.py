from Engine.importmodule import *

from Engine.IO_AB     import setup_templates_tel, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips     import basicclip_above
from Engine.contfit   import A0cont
from Engine.classes   import fitobjs,inparamsA0,orderdict_cla
from Engine.rebin_jv  import rebin_jv
from Engine.rotint    import rotint
from Engine.Telfitter import telfitter
from Engine.opt       import optimizer, fmod
from Engine.outplotter import outplotter_tel
from Engine.detect_peaks import detect_peaks
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def MPinst(args, inparam, jerp, orders, masterbeam, i):
    # Main function for A0 fitting that will be threaded over by multiprocessing

    order = orders[jerp]            # current looped order
    night = str(inparam.nights[i])  # multiprocess assigned night
    firstorder = orders[0]          # First order that will be analyzed, related to file writing

    print('Working on order {:02d}/{:02d} ({}), night {}/{} ({}) PID:{}...'.format(int(jerp+1),
                                                                                 len(orders),
                                                                                 order,
                                                                                 i+1,
                                                                                 len(inparam.nights),
                                                                                 night,
                                                                                 mp.current_process().pid) )

    #-------------------------------------------------------------------------------

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
    # bound_cut = [150, 150]


    ### Load relevant A0 spectrum
    x, a0wavelist, a0fluxlist, u = init_fitsread(inparam.inpath,
                                                 'A0',
                                                 'combined'+str(masterbeam),
                                                 night,
                                                 order,
                                                 f'{int(inparam.tags[night]):04d}',
                                                 args.band,
                                                 bound_cut)
    #-------------------------------------------------------------------------------

    s2n = a0fluxlist/u
    if np.nanmedian(s2n) < float(args.SN_cut):
        logger.warning('  --> Bad S/N {:1.3f} < {} for {}{}, SKIP'.format( np.nanmedian(s2n), args.SN_cut, night, masterbeam))

        pre_err = True
        logger.warning(f'  --> NIGHT {night}, ORDER {order} HIT ERROR DURING PRE_OPT')
        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)
        else:
            hh = fits.open('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)

        return

    # Trim obvious outliers above the blaze (i.e. cosmic rays)
    nzones = 12
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(x,a0fluxlist,nzones);
    a0u        = basicclip_above(u,a0fluxlist,nzones);     a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(a0x,a0fluxlist,nzones);
    a0u        = basicclip_above(a0u,a0fluxlist,nzones);   a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);

    # Compute rough blaze function estimate. Better fit will be provided by Telfit later.
    continuum    = A0cont(a0wavelist,a0fluxlist,night,order,args.band)
    a0masterwave = a0wavelist.copy()
    a0masterwave *= 1e4

    # Trim stellar template to relevant wavelength range
    mwave_in, mflux_in = stellarmodel_setup(a0wavelist, inparam.mwave0, inparam.mflux0)

    # Trim telluric template to relevant wavelength range
    # Normalize continuum level of telluric atlas in the given band
    if args.band == 'H':
        contlevel = np.max(inparam.satm[(inparam.watm > 15000) & (inparam.watm < 18000)])
    else:
        contlevel = np.max(inparam.satm[(inparam.watm > 20000) & (inparam.watm < 24000)])

    # Trim telluric template to relevant wavelength range
    satm_in = inparam.satm[(inparam.watm > np.min(a0wavelist)*1e4 - 11) & (inparam.watm < np.max(a0wavelist)*1e4 + 11)]
    watm_in = inparam.watm[(inparam.watm > np.min(a0wavelist)*1e4 - 11) & (inparam.watm < np.max(a0wavelist)*1e4 + 11)]
    satm_in /= contlevel

    # Get initial guess for cubic wavelength solution from reduction pipeline
    f = np.polyfit(a0x, a0wavelist, 3)
    par9in = f[0]*1e4;
    par8in = f[1]*1e4;
    par7in = f[2]*1e4;
    par6in = f[3]*1e4;

    # Determine whether IGRINS mounting was loose or night for the night in question
    if (int(night) < 20180401) or (int(night) > 20190531):
        IPpars = inparam.ips_tightmount_pars[args.band][masterbeam][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][masterbeam][order]

    # start at bucket loc = 1250 +- 100, width = 250 +- 100, depth = 100 +- 5000 but floor at 0
    if args.band == 'H':
        centerloc = 1250
    else:
        centerloc = 1180

    ### Initialize parameter array for optimization as well as half-range values for each parameter during
    ### the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and
    ### in between optimization steps.
    parA0 = np.array([0.0,           # 0: The shift of the stellar template (km/s)
                      0.0,           # 1: The scale factor for the stellar template
                      0.0,           # 2: The shift of the telluric  template (km/s)
                      1.0,           # 3: The scale factor for the telluric template
                      0.0,           # 4: vsini (km/s)
                      IPpars[2],     # 5: The instrumental resolution (FWHM) in pixels
                      par6in,        # 6: Wavelength 0-pt
                      par7in,        # 7: Wavelength linear component
                      par8in,        # 8: Wavelength quadratic component
                      par9in,        # 9: Wavelength cubic component
                      1.0,           #10: Continuum zero point
                      0.0,           #11: Continuum linear component
                      0.0,           #12: Continuum quadratic component
                      IPpars[1],     #13: Insrumental resolution linear component
                      IPpars[0],     #14: Instrumental resolution quadratic component
                      centerloc,     #15: Blaze dip center location
                      330,           #16: Blaze dip full width
                      0.05,          #17: Blaze dip depth
                      90,            #18: Secondary blaze dip full width
                      0.05,          #19: Blaze dip depth
                      0.0,           #20: Continuum cubic component
                      0.0,           #21: Continuum quartic component
                      0.0,           #22: Continuum quintic component
                      0.0])          #23: Continuum hexic component


    # Make sure data is within telluric template range (shouldn't do anything)
    a0fluxlist = a0fluxlist[(a0wavelist*1e4 > np.min(watm_in)+5) & (a0wavelist*1e4 < np.max(watm_in)-5)]
    a0u        = a0u[       (a0wavelist*1e4 > np.min(watm_in)+5) & (a0wavelist*1e4 < np.max(watm_in)-5)]
    a0x        = a0x[       (a0wavelist*1e4 > np.min(watm_in)+5) & (a0wavelist*1e4 < np.max(watm_in)-5)]
    continuum  = continuum[ (a0wavelist*1e4 > np.min(watm_in)+5) & (a0wavelist*1e4 < np.max(watm_in)-5)]
    a0wavelist = a0wavelist[(a0wavelist*1e4 > np.min(watm_in)+5) & (a0wavelist*1e4 < np.max(watm_in)-5)]

    # Define main spectrum
    s = a0fluxlist.copy(); x = a0x.copy(); u = a0u.copy();

    # Collect all fit variables into one class
    fitobj = fitobjs(s, x, u, continuum, watm_in, satm_in, mflux_in, mwave_in, [], masterbeam, np.array([],dtype=int))

    #                            |0    1    2    3  |  | 4 |  | 5 |   | 6    7    8           9  |    |10 11 12|  |13 14|    |15    16    17   18    19|  |20   21   22    23 |
    dpars = {'cont' :   np.array([0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,  0.0, 0.0,        0.0,    1e7, 1, 1,    0, 0,    10.0, 20.0, 0.2, 50.0, 0.2,   1.0, 1.0, 1.0, 1.0 ]),
             'twave':   np.array([0.0, 0.0, 0.0, 1.0,   0.0,   0.0,   10.0, 10.0, 5.00000e-5, 1e-7,   0.0, 0, 0,    0, 0,     0.0,  0.0, 0.0,  0.0, 0.0,   0.0, 0.0, 0.0, 0.0 ]),
             'ip'   :   np.array([0.0, 0.0, 0.0, 0.0,   0.0,   0.5,    0.0,  0.0, 0.0,        0.0,    0.0, 0, 0,    0, 0,     0.0,  0.0, 0.0,  0.0, 0.0,   0.0, 0.0, 0.0, 0.0 ])}
    if masterbeam == 'B':
        dpars['cont'] = np.array([0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,  0.0, 0.0,        0.,     1e7, 1, 1,    0, 0,     0.0,  0.0, 0.0,  0.0, 0.0,   1.0, 1.0, 1.0, 1.0  ])
    else:
        if (args.band == 'K') and (order == 3):
            parA0[19] = 0.
        #                            |0    1    2    3  |  | 4 |  | 5 |   | 6    7    8           9  |    |10 11 12|  |13 14|    |15    16    17   18    19|  |20   21   22    23 |
            dpars['cont'] = np.array([0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0, 0.0, 0.0,         0.,     1e7, 1, 1,    0, 0,    10.0, 20.0, 0.2, 50.0, 0.0,   1.0, 1.0, 1.0, 1.0  ])
    #-------------------------------------------------------------------------------

    # Use quadratic blaze correction for order 13; cubic for orders 6, 14, 21; quartic for orders 16 and 22
    if args.band == 'H':
        if int(order) in [13]:
            dpars['cont'][20] = 0.; dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
        elif int(order) in [6,14,21]:
            dpars['cont'][21] = 0.; dpars['cont'][22] = 0.; dpars['cont'][23] = 0.;
        else:
            pass

    # Initialize an array that puts hard bounds on vsini and the instrumental resolution to make sure they do not diverge to unphysical values
    optimize = True
    par_in = parA0.copy()
    if masterbeam == 'B':
        hardbounds = [par_in[4]  - 0,                 par_in[4]  + 0,
                      par_in[5]  - dpars['ip'][5],    par_in[5]  + dpars['ip'][5]
                     ]
    else:
        hardbounds = [par_in[4]  - 0,                 par_in[4]  + 0,
                      par_in[5]  - dpars['ip'][5],    par_in[5]  + dpars['ip'][5],
                      par_in[15] - dpars['cont'][15], par_in[15] + dpars['cont'][15],
                      par_in[16] - dpars['cont'][16], par_in[16] + dpars['cont'][16],
                      0.,                             par_in[17] + dpars['cont'][17],
                      par_in[18] - dpars['cont'][18], par_in[18] + dpars['cont'][18],
                      0.,                             par_in[19] + dpars['cont'][19]
                     ]

    if hardbounds[0] < 0:
        hardbounds[0] = 0
    if hardbounds[2] < 0:
        hardbounds[2] = 1

    # Begin optimization.
    # For every pre-Telfit spectral fit, first fit just template strength/rv/continuum, then just wavelength solution, then template/continuum again, then ip,
    # then finally wavelength. Normally would fit for all but wavelength at the end, but there's no need for the pre-Telfit fit, since all we want
    # is a nice wavelength solution to feed into Telfit.

    pre_err = False

    cycles = 3
    optgroup = ['twave', 'cont',
                'cont', 'twave', 'cont',
                'twave',
                'ip', 'twave',  'cont',
                'ip', 'twave',  'cont',
                'twave']
    try:
        
        go = 1; misfit_flag_low = 0; restarted = False;

        while go == 1:

            parstart = par_in.copy()

            if misfit_flag_low == 1:
                parstart[3] = 0.5

            if misfit_flag_low == 2: 
                print(breaker) # deliberately throw error to enter except statement

            nk = 1
            for nc, cycle in enumerate(np.arange(cycles), start=1):

                for optkind in optgroup:
                    start = time.time()
                    parfit_1 = optimizer(parstart, dpars[optkind], hardbounds, fitobj, optimize)

                    if parfit_1[3] < 0.1:
                        misfit_flag_low += 1
                        break

                    if args.debug == True:
                        outplotter_tel(parfit_1,fitobj,'{}_{}_beforeparfit_{}{}'.format(order,night,nk,optkind),inparam, args, order)
                    parstart = parfit_1.copy()
                    nk += 1

                if misfit_flag_low == 1 and restarted == False:
                    restarted = True
                    print('alright')
                    break


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

                    # Redo rough blaze fit in case hot pixels were throwing it off
                    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)
                    mask = np.ones_like(w,dtype=bool)
                    mask[CRmaskF] = False
                    continuum    = A0cont(w[mask]/1e4,s[mask],night,order,args.band)
                    continuum    = rebin_jv(w[mask],continuum,w,False)
                    fitobj = fitobjs(s, x, u, continuum, watm_in, satm_in, mflux_in, mwave_in, [], masterbeam, CRmaskF)

            if misfit_flag_low == 0:
                
                parfit = parfit_1.copy()

                # If dip present, correct it out of data before running Telfit to enable better fit
                if masterbeam == 'A':
                    cont = parfit[10] + parfit[11]*fitobj.x+ parfit[12]*(fitobj.x**2) + parfit[20]*(fitobj.x**3) + parfit[21]*(fitobj.x**4) + parfit[22]*(fitobj.x**5) + parfit[23]*(fitobj.x**6)
                    cont0 = cont.copy()
                    bucket = np.zeros_like(cont)
                    bucket[(fitobj.x >= (parfit[15]-parfit[16]/2)) & (fitobj.x <= (parfit[15]+parfit[16]/2))] = parfit[17]
                    bucket[(fitobj.x >= (parfit[15]+parfit[16]/2-parfit[18])) & (fitobj.x <= (parfit[15]+parfit[16]/2))] += parfit[19]
                    cont -= bucket

                    cont *= continuum
                    cont0 *= continuum
                    justdip = cont/cont0
                    a0fluxlist /= justdip

                go = 0; break;
                
    except:
        pre_err = True
        logger.warning(f'  --> NIGHT {night}, ORDER {order} HIT ERROR DURING PRE_OPT')
        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)
        else:
            hh = fits.open('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)


    #-------------------------------------------------------------------------------
    if not pre_err:

        if inparam.plotfigs: # Plot results
            outplotter_tel(parfit, fitobj, f'BeforeTelFit_Order{order}_{night}_{masterbeam}', inparam, args, order)

        # Get best fit wavelength solution
        a0w_out_fit = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)

        # Trim stellar template to new relevant wavelength range
        mwave_in,mflux_in = stellarmodel_setup(a0w_out_fit/1e4, inparam.mwave0, inparam.mflux0)

        # Feed this new wavelength solution into Telfit. Returns high-res synthetic telluric template, parameters of that best fit, and blaze function best fit
        watm1, satm1, telfitparnames, telfitpars, a0contwave, continuum = telfitter(a0w_out_fit,a0fluxlist,a0u,inparam,night,order,args,masterbeam)
    else:
        pass
    #-------------------------------------------------------------------------------

    # If Telfit encountered error (details in Telfitter.py), skip night/order combo
    if pre_err:
        pass
    elif len(watm1) == 1:
        logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}')

        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name=f'ERRORFLAG{order}', array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        # If first time writing fits file, make up filler primary hdu
        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)
        else:
            hh = fits.open('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)

    else: # If Telfit exited normally, proceed.
        #  Save best blaze function fit
        continuum = rebin_jv(a0contwave,continuum,a0w_out_fit,False)

        # Write out table to fits file with errorflag = 0
        c0  = fits.Column(name='ERRORFLAG'+str(order),      array=np.array([0]),            format='K')
        c1  = fits.Column(name='WAVE'+str(order),           array=a0w_out_fit,              format='D')
        c2  = fits.Column(name='BLAZE'+str(order),          array=continuum,                format='D')
        c3  = fits.Column(name='X'+str(order),              array=a0x,                      format='D')
        c4  = fits.Column(name='INTENS'+str(order),         array=a0fluxlist,               format='D')
        c5  = fits.Column(name='SIGMA'+str(order),          array=a0u,                      format='D')
        c6  = fits.Column(name='WATM'+str(order),           array=watm1,                    format='D')
        c7  = fits.Column(name='SATM'+str(order),           array=satm1,                    format='D')
        c8  = fits.Column(name='TELFITPARNAMES'+str(order), array=np.array(telfitparnames), format='8A')
        c9  = fits.Column(name='TELFITPARS'+str(order),     array=telfitpars,               format='D')
        c10 = fits.Column(name='PARFIT',                    array=parfit,                   format='D')
        cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
        hdu_1 = fits.BinTableHDU.from_columns(cols)


        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)
        else:
            hh = fits.open('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_{}treated_{}.fits'.format(inparam.outpath, night, masterbeam, args.band), overwrite=True)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def mp_run(args, inparam, Nthreads, jerp, orders, nights, masterbeam):
    # Multiprocessing convenience function
    pool = mp.Pool(processes = Nthreads)
    func = partial(MPinst, args, inparam, jerp, orders,masterbeam)
    outs = pool.map(func, np.arange(len(nights)))
    pool.close()
    pool.join()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def use_w(args):
    # Load wavelength regions list file
    try:
        bounddata = Table.read(f'./Input/UseWv/WaveRegions_{args.WRegion}_{args.band}.csv', format='csv')
    except IOError:
        sys.exit(f'WaveRegions FILE "./Input/UseWv/WaveRegions_{args.WRegion}_{args.band}.csv" NOT FOUND!')

    wavesols = pd.read_csv(f'./Input/UseWv/WaveSolns_{args.band}.csv')
#-------------------------------------------------------------------------------
    with open(f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv','w') as filew:
        filew.write('order, start,  end, masks\n')

        m_order  = np.array(bounddata['order'])
        starts   = np.array(bounddata['start'])
        ends     = np.array(bounddata['end'])
        ords     = list( sorted(orderdict_cla().orderdict[args.band].keys()) )

        Ostarts  = [orderdict_cla().orderdict[args.band][k][0] for k in ords]
        Oends    = [orderdict_cla().orderdict[args.band][k][1] for k in ords]
        labels   = []

        m_orders_unique = np.unique(m_order)

        # For each order specified, find what pixel numbers correspond to the wavelength bounds presented.
        # If multiple wavelength bounds given for a single order, output a pixel mask between the two, as well.
        for o in range(len(m_orders_unique)):
            pixs = [];
            mini = np.where(m_order == m_orders_unique[o])[0]
            for j in range(len(mini)):
                i = mini[j]

                wavebounds = [starts[i],ends[i]]
                wO   = wavesols['w'+str(m_orders_unique[o])]
                pixO = wavesols['x'+str(m_orders_unique[o])];
                pix  = [pixO[(np.argmin(abs(wO-wavebounds[k])))] for k in [0,1]]
                pixs = pixs + pix

            pixsS = list(sorted(pixs))
            q = pixsS[1:-1]
            if len(pixsS) == 2:
                filew.write('{}, {}, {},[]\n'.format(m_orders_unique[o],pixsS[0],pixsS[-1]))
            else:
                filew.write('{}, {}, {},"{}"\n'.format(m_orders_unique[o],pixsS[0],pixsS[-1],[[first,second] for first,second in  zip(q[0::2], q[1::2])]))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 1',
                                     description = '''
                                     This step 1) defines the wavelength regions to be analyzed based on user specification \n
                                     2) generates a synthetic, high-resolution telluric
                                     template for use in later model fits on a night by night basis.  \n
                                     Note that only target star observations will have their fits limited to the wavelength regions specified. \n
                                     For A0 observations, only the orders specified will be analyzed, but each order will be fit as far as there is significant telluric absoprtion.
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name, no space",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
                        type=int,   default=int(1))

    parser.add_argument("-SN",      dest="SN_cut",           action="store",
                        help="Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 ",
                        type=str,   default='50')

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",        action="store_true",
                        help="If set, will generate plots of A0 fitting results under ./Output/A0Fits/*target/fig/")

    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want to process all nights under the ./Input/*target/ folder, specify an array of night you wish to process here. e.g., [20181111,20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If set, DeBug logging will be output, as well as (lots of) extra plots under ./Temp/Debug/*target_*band/main_step1")
    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.9')
    args = parser.parse_args()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'
#-------------------------------------------------------------------------------
    # Create output directories as needed
    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/A0Fits/figs_{args.band}')

    outpath = f'./Output/{args.targname}_{args.band}/A0Fits'
#-------------------------------------------------------------------------------
    # Handle logger
    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}_A0Fits.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('####################################################################################\n')
    print(f'Fetching Wavelength Regions to be Analyzed for {args.targname}...')
    time.sleep(2)

    use_w(args)

    print('Fetching Done!')
    print(f'File "XRegions_{args.WRegion}_{args.band}.csv" saved under "./Input/UseWv/"')
    #time.sleep(5)

    #-------------------------------------------------------------------------------

    print('###############################################################\n')
    logger.info(f'Using TelFit to create high-resolution, synthetic telluric templates based off the telluric standards \nassociated with {args.targname} on a night by night basis...')
    print('This will take a while..........')
    print('\n')

    # Read in newly created pixel regions file to get list of orders to analyze.
    # Note that only target star observations will have their fits limited to the wavelength regions specified.
    # For A0 observations, only the orders specified will be analyzed, but each order will be fit as far as there is significant telluric absoprtion.
    bounddata = Table.read(f'./Input/UseWv/XRegions_{args.WRegion}_{args.band}.csv', format='csv')
    starts  = np.array(bounddata['start'])
    ends    = np.array(bounddata['end'])
    orders  = np.array(bounddata['order'], dtype=int)
    xbounddict = {orders[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}

    ## Collect relevant file information from Predata files
    A0data = Table.read(f'./Input/Prepdata/Prepdata_A0_{args.targname}.txt', format='ascii')

    ind    = [i != 'NA' for i in A0data['humid']]
    humids = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['humid'])}
    tags   = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['tag'])}
    obs    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['obs'])}
    temps  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['temp'])}
    zds    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['zd'])}
    press  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['press'])}
    nightsFinal = np.array(list(sorted(set(A0data[ind]['night']))))

    # Take subset of nights, if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=int)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(f'NIGHT {nnn} EITHER HAS NO CORRESPONDING A0 OR WAS NOT FOUND UNDER "./Input/{args.targname}"')

        nightsFinal = nightstemp
        logger.info(f'Only processing nights: {nightsFinal}')

    logger.info(f'Analyze {len(nightsFinal)} nights')

    #time.sleep(6)
    print('\n')

    # print('For paper plot!')
    # if args.band == 'K':
    #     orders = np.array([2, 3, 4, 5, 6, 7, 8, 10, 14, 16])
    # elif args.band=='H':
    #     orders = np.array([2, 3, 4, 6, 13, 14, 16, 20, 21, 22])
    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates_tel()

    inparam = inparamsA0(inpath,outpath,args.plotfigs,tags,nightsFinal,humids,
                         temps,zds,press,obs,watm,satm,mwave0,mflux0,cdbs_loc,xbounddict,None)

    #-------------------------------------------------------------------------------

    # Run order by order, multiprocessing over nights within an order
    print('Processing the A nods first...')
    for jerp in range(len(orders)):
        outs = mp_run(args, inparam, args.Nthreads, jerp, orders, nightsFinal,'A')

    print('A nods done! Halfway there! \n Now processing the B nods...')
    for jerp in range(len(orders)):
        outs = mp_run(args, inparam, args.Nthreads, jerp, orders, nightsFinal,'B')

    print('\n')
    logger.info('A0 Fitting Done!')

    end_time = datetime.now()
    logger.info(f'A0 Fitting using TelFit finished, Duration: {end_time - start_time}')
    print('You can start to run main_step2.py for RV initial guess')
    print('####################################################################################')
