from Engine.importmodule import *
import sys
from numpy.polynomial import chebyshev
from scipy.interpolate import splrep,splev
from Engine.rotint import rotint
from Engine.macbro_dynamic    import macbro_dyn
from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.clips      import basicclip_above
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.opt        import optimizer, fmod, fmod_conti

#-------------------------------------------------------------------------------

matplotlib.use('Qt5Agg')

def templates_semibroadened(par,fitobj,broad="none"):

    watm_in = fitobj.watm_in
    satm_in = fitobj.satm_in
    mwave = fitobj.mwave_in
    mflux = fitobj.mflux_in

    #Make the wavelength scale
    initwave = fitobj.initwave.copy()
    xgrid = (initwave-np.median(initwave)) / (np.max(initwave)-np.min(initwave))
    dx = chebyshev.chebval(xgrid, par[6:10])
    w = initwave + dx

    if np.all(np.diff(w) > 0) == False:
        sys.exit('WAVE ERROR 1 - Hitting negative wavelength solution ',
                    'for some reason - pars: {}'.format(par[6:10]))
        return 1e10

    # Define the speed of light in km/s and other useful quantities
    c = 2.99792458e5
    npts = len(w)

    # Apply velocity shifts and scale
    watm = watm_in*(1.+par[2]/c)
    satm = satm_in**par[3]

    #Verify that new wavelength scale is a subset of old telluric wavelength scale.
    if (w[0] < watm[0]) or (w[-1] > watm[-1]):
        sys.exit('WAVE ERROR 2: w subset of watm, w goes from ' +
                    str(w[0]) + ' to '+str(w[-1]) +
                    ' and watm goes from ' + str(watm[0]) +
                    ' to ' + str(watm[-1])
                    )
        return 1e10

    if mwave is not None:

        wspot = mwave*(1.+par[0]/c)
        sspot = mflux**par[1]

        #Verify that new wavelength scale is a subset of stellar wavelength scale.
        if (w[0] < wspot[0]) or (w[-1] > wspot[-1]):
            sys.exit('WAVE ERROR 3:  w not subset of wspot, w goes from '
                        + str(w[0]) + ' to '+str(w[-1]) +
                        ' and wspot goes from ' + str(wspot[0]) +
                        ' to '+str(wspot[-1])
                        )
            return 1e10

        #Verify that stellar wavelength scale is a subset of telluric wavelength scale.
        if (wspot[0] < watm[0]) or (wspot[-1] > watm[-1]):
            sys.exit('WAVE ERROR 3: wspot not subset of satm, wspot goes from '
                        + str(wspot[0]) + ' to ' + str(wspot[-1]) +
                        ' and watm goes from ' + str(watm[0]) +
                        ' to ' + str(watm[-1]))
            return 1e10

        if broad == "none":
            spot2 = rebin_jv(wspot,sspot,w,False)
            satm2 = rebin_jv(watm,satm,w,False)
            return spot2,satm2

        vsini = par[4]

        # Rotationally broaden stellar template
        if vsini >= 0.5:
            wspot2,rspot2 = rotint(wspot,sspot,vsini)
        else:
            wspot2 = wspot
            rspot2 = sspot

        spot2 = rebin_jv(wspot2,rspot2,w,False)

        return spot2,None



def modtool(args,jerp,nightsbox,tagbox,parfitbox,inparam,index):

    c = 2.99792458e5

    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(args)
    order = orders[jerp]

    night = nightsFinal[index]

    xbounds = xbounddict[order]

    bound_cut = _setup_bound_cut(inparam.bound_cut_dic,
                                    args.band, order)

    indnight  = np.where(nightsbox == night)[0]
    parfittags    =  parfitbox[indnight]
    tags          =     tagbox[indnight]


    tags1 = ['%4.0f' % float(t) for t in tags]
    tags2 = np.array([t.replace(' ','0') for t in tags1])

    firsttag = True; pre_err = True;

    #print(night)
    #plt.figure(figsize=(12,10))

    for beam in ['A','B']:
        if beam == 'A':
            taglist = tagsA[night]
        else:
            taglist = tagsB[night]

        if order == 9:
            A0loc = f'./Output/{args.targname}_{args.band}/A0Fits_order9/{night[:8]}A0_{beam}treated_{args.band}.fits'
        else:
            A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

        try:
            hdulist = fits.open(A0loc)
        except IOError:
            print(f'  --> No A0-fitted template for night {night}, skipping...')
            badtags = taglist.copy()
            continue

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
            print(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: {order} NIGHT: {night}, skipping...')
            badtags = taglist.copy()
            continue

        badtags = []

        for tag in taglist:

            if len(parfittags[(tag == tags2)]) == 0:
                badtags.append(tag)
                #print(night,tag)
                continue
            else:
                parfit = parfittags[(tag == tags2)][0]

            watm = tbdata['WATM'+str(order)]
            satm = tbdata['SATM'+str(order)]
            a0wave = tbdata['WAVE'+str(order)]
            a0flux = tbdata['INTENS'+str(order)]
            a0contx    = tbdata['X'+str(order)]
            continuum  = tbdata['BLAZE'+str(order)]

            # Remove extra rows leftover from having columns of unequal length
            satm = satm[(watm != 0)]
            watm = watm[(watm != 0)]
            a0flux = a0flux[(a0wave != 0)]
            a0wave = a0wave[(a0wave != 0)]
            satm[(satm < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi
            a0contx = a0contx[(continuum != 0)]
            continuum = continuum[(continuum != 0)]

            a0flat = a0flux / continuum

            # Load target spectrum
            x,wave,s,u = init_fitsread(f'./Input/{args.targname}/{night}/{beam}/',
                                'target',
                                'separate',
                                night,
                                order,
                                tag,
                                args.band,
                                bound_cut)

            # Trim obvious outliers above the blaze (i.e. cosmic rays)
            nzones = 5
            x = basicclip_above(x,s,nzones)
            wave = basicclip_above(wave,s,nzones)
            u = basicclip_above(u,s,nzones)
            s = basicclip_above(s,s,nzones)
            x = basicclip_above(x,s,nzones)
            wave = basicclip_above(wave,s,nzones)
            u = basicclip_above(u,s,nzones)
            s = basicclip_above(s,s,nzones)

            # Cut spectrum to within wavelength regions defined in input list
            s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
            u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
            wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
            x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

            # Trim telluric template to data range +- 15 AA. If telluric
            # template buffer is cut short because A0 lines didn't extend
            # far past data range, cut data range accordingly.
            satm_in = satm[(watm > np.min(wave_piece)*1e4 - 10) \
                                & (watm < np.max(wave_piece)*1e4 + 10)]
            watm_in = watm[(watm > np.min(wave_piece)*1e4 - 10) \
                                & (watm < np.max(wave_piece)*1e4 + 10)]

            s_piece	= s_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                                & (wave_piece*1e4 < np.max(watm_in)-10)]
            u_piece	= u_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                                & (wave_piece*1e4 < np.max(watm_in)-10)]
            x_piece	= x_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                                & (wave_piece*1e4 < np.max(watm_in)-10)]
            wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+10) \
                                & (wave_piece*1e4 < np.max(watm_in)-10)]

            # Trim stellar template to data range +- 10 AA
            mflux_in = inparam.mflux0[
                    (inparam.mwave0 > np.min(wave_piece)*1e4 - 5) \
                    & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5)
                ]
            mwave_in = inparam.mwave0[
                    (inparam.mwave0 > np.min(wave_piece)*1e4 - 5) \
                    & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5)
                ]

            Rstell = np.median(np.diff(mwave_in))
            Rtell = np.median(np.diff(watm_in))
            if Rstell < Rtell:
                sys.exit(f'Telluric template resolution ({round(Rtell,4)} AA) '
                            'must be finer than stellar template resolution '
                            '({round(Rstell,4)} AA) !')

            # Rebin stellar template to uniform wavelength scale
            dstep = Rstell
            nstep = int((mwave_in[-1]-mwave_in[0])/dstep)
            mwave1 = np.linspace(mwave_in[0],mwave_in[-1],nstep)
            mflux1 = rebin_jv(mwave_in,mflux_in,mwave1,False)
            mwave_in = mwave1.copy(); mflux_in = mflux1.copy()
            mwave_in = mwave_in[1:-1]
            mflux_in = mflux_in[1:-1]


            # Normalize continuum from A0 to flux scale of data
            continuum /= np.nanmedian(continuum)
            continuum *= np.nanpercentile(s_piece,99)

            continuum_in = rebin_jv(a0contx,continuum,x_piece,False)

            # Get initial guess for cubic wavelength solution from reduction pipeline
            f = np.polyfit(x_piece,wave_piece,3)
            q = np.poly1d(f)
            initwave = q(x_piece)*1e4

            fitobj = FitObjs(s_piece, x_piece, u_piece, continuum_in, watm_in,
                        satm_in, mflux_in, mwave_in,
                        ast.literal_eval(inparam.maskdict[order]),
                        beam, [np.array([], dtype=int),
                        np.array([], dtype=int)],
                        initwave, [])

            parfitS = parfit.copy(); parfitS[3] = 0;
            parfitT = parfit.copy(); parfitT[1] = 0;

            wave_out, smod_tell, cont_dontuse, c2_dontuse = fmod_conti(parfitT,fitobj)
            stell = fitobj.s / smod_tell
            s2n = fitobj.s / fitobj.u

            # Regular model fit
            wave_reg,   mod_out, cont_out, c2_out = fmod_conti(parfit,fitobj)
            wave_reg, stell_reg, cont_out, c2_out = fmod_conti(parfitS,fitobj)
            wave_reg,  tell_reg, cont_out, c2_out = fmod_conti(parfitT,fitobj)
            continuum_out = cont_out*c2_out

            # Unbroadened stellar model
            stell_unbrod,tell_unbrod = templates_semibroadened(parfit,fitobj,broad='none')

            # Rotationally broadened stellar model
            stell_brod,trash = templates_semibroadened(parfit,fitobj,broad='rot')

            # Regular data
            flux_reg = fitobj.s.copy()

            # Flattened data
            flux_flat = fitobj.s/(continuum_out)

            # Flattened, telluric-corrected data
            flux_corr = fitobj.s/tell_reg

            wave_shift = wave_reg*(1 - parfit[0]/c)

            pre_err = False;

            c0 = fits.Column(name = f'ERRORFLAG{order}',
                                array = np.array([0]),
                                format='K'
                                )
            c1 = fits.Column(name='WAVE_RAW'+str(order),       array=wave_reg,                   format='D')
            c2 = fits.Column(name='WAVE_ADJ'+str(order),       array=wave_shift,                   format='D')
            c3 = fits.Column(name='FLUX_RAW',       array=flux_reg,                  format='D')
            c4 = fits.Column(name='FLUX_CORR',       array=flux_corr,                  format='D')
            c5 = fits.Column(name='CONT',       array=continuum_out,                  format='D')
            c50 = fits.Column(name='S2N',       array=s2n,                  format='D')
            c6 = fits.Column(name='STELL',       array=stell_reg,                  format='D')
            c7 = fits.Column(name='STELL_ROTBROADONLY',       array=stell_brod,                  format='D')
            c8 = fits.Column(name='STELL_NOBROAD',       array=stell_unbrod,                  format='D')
            c9 = fits.Column(name='TELL',       array=tell_reg,                  format='D')
            c10 = fits.Column(name='TELL_NOBROAD',       array=tell_unbrod,                  format='D')
            c11 = fits.Column(name='A0WAVE',       array=a0wave,                   format='D')
            c12 = fits.Column(name='A0FLUX',       array=a0flat,                  format='D')
            c13 = fits.Column(name='RV',                    array=np.array([parfit[0]]),     format='D')
            c14 = fits.Column(name='BVC',                    array=np.array([inparam.bvcs[night+tag]]),     format='D')
            cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c50,c6,c7,c8,c9,c10,c11,c12,c13,c14])
            hdu_1 = fits.BinTableHDU.from_columns(cols)


            if jerp == 0:
                bleh = np.ones((3,3))
                primary_hdu = fits.PrimaryHDU(bleh)
                hdul = fits.HDUList([primary_hdu,hdu_1])
                hdul.writeto('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, tag), overwrite=True)
            else:
                hh = fits.open('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, tag))
                hh.append(hdu_1)
                hh.writeto('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, tag), overwrite=True)

            #plt.plot(wave_shift,flux_corr,alpha=0.4)
            if firsttag:
                stellstack = flux_corr.copy(); masterwave = wave_shift.copy(); s2nstack = s2n.copy();
                firsttag = False
            else:
                stellnew = rebin_jv(wave_shift,flux_corr,masterwave,False)
                s2nnew = rebin_jv(wave_shift,s2n,masterwave,False)
                stellnew[(masterwave < wave_shift[0]) | (masterwave > wave_shift[-1])] = np.nan
                s2nnew[(masterwave < wave_shift[0]) | (masterwave > wave_shift[-1])] = np.nan
                stellstack = np.vstack((stellstack,stellnew))
                s2nstack   = np.vstack((s2nstack,s2nnew))

    for badtag in badtags:

        c0 = fits.Column(name = f'ERRORFLAG{order}',
                            array = np.array([1]),
                            format='K'
                            )
        cols = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)
        if jerp == 0:
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, badtag), overwrite=True)
        else:
            hh = fits.open('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, badtag))
            hh.append(hdu_1)
            hh.writeto('{}/Cutout_{}_{}.fits'.format(inparam.outpath, night, badtag), overwrite=True)

    if not pre_err:
        #plt.show()

        try:
            Cstell = np.array([np.nanmean(stellstack[:,jj]) for jj in range(len(stellstack[0,:]))])
            Cs2n   = np.array([np.sqrt(np.nansum(s2nstack[:,jj]**2)) for jj in range(len(s2nstack[0,:]))])
        except IndexError:
            Cstell = stellstack.copy(); Cs2n = s2nstack.copy();

        c0 = fits.Column(name = f'ERRORFLAG{order}',
                            array = np.array([0]),
                            format='K'
                            )
        c1 = fits.Column(name='WAVE'+str(order),       array=masterwave,                   format='D')
        c2 = fits.Column(name='FLUX',       array=Cstell,                  format='D')
        c3 = fits.Column(name='UNC',        array=Cstell/Cs2n,                  format='D')
        cols = fits.ColDefs([c0,c1,c2,c3])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if jerp == 0:
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night), overwrite=True)
        else:
            hh = fits.open('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night))
            hh.append(hdu_1)
            hh.writeto('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night), overwrite=True)

    else:

        c0 = fits.Column(name = f'ERRORFLAG{order}',
                            array = np.array([1]),
                            format='K'
                            )
        cols = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)
        if jerp == 0:
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night), overwrite=True)
        else:
            hh = fits.open('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night))
            hh.append(hdu_1)
            hh.writeto('{}/Cutout_{}Combined.fits'.format(inparam.outpath, night), overwrite=True)
