
from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step6a

from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.outplotter import outplotter_23
from Engine.detect_peaks import detect_peaks
#from Engine.LS         import LS
from Engine.plot_tool import modtool
from scipy.stats import pearsonr
from Engine.rotint import rotint
from Engine.macbro_dynamic    import macbro_dyn
from Engine.rebin_jv import rebin_jv

from scipy.interpolate import splrep,splev #, interp1d
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def makespec(nightsbox0,parfitbox,tagbox,args,indict,jerp,i):

        orders = indict['orders']
        order = orders[jerp]
        night = indict['nights'][i]
        indnight  = np.where(nightsbox0 == night)[0]
        parfittags    =  parfitbox[indnight]
        tags          =     tagbox[indnight]

        tags1 = ['%4.0f' % float(t) for t in tags]
        tags2 = np.array([t.replace(' ','0') for t in tags1])

        firsttag = True
        c = 2.99792458e5
        vsini = indict['vsini']
        fluxratio = float(args.fluxratio)

        for beam in ['A','B']:
            skipnight = False
            if beam == 'A':
                taglist = indict['tagsA'][night]
            else:
                taglist = indict['tagsB'][night]

            A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/{night[:8]}A0_{beam}treated_{args.band}.fits'

            hdulist = fits.open(A0loc)
            num_orders = 0
            for rri in range(25):
                try:
                    hdulist[rri].columns[0].name[9:]
                    num_orders += 1
                except:
                    continue
            fits_layer = [ rri for rri in np.arange(num_orders)+1 \
                            if np.int(hdulist[rri].columns[0].name[9:]) == order ][0]
            tbdata = hdulist[ fits_layer ].data
            flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]
            if flag == 1:
                skipnight = True
            else:
                watm = tbdata['WATM'+str(order)]
                satm = tbdata['SATM'+str(order)]
                a0contx    = tbdata['X'+str(order)]
                continuum  = tbdata['BLAZE'+str(order)]
                satm = satm[(watm != 0)]
                watm = watm[(watm != 0)]
                a0contx = a0contx[(continuum != 0)]
                continuum = continuum[(continuum != 0)]
                satm[(satm < 1e-4)] = 0.

                watm_in = watm.copy(); satm_in = satm.copy();
                mwave1_in = indict['mwave1'].copy(); mflux1_in = indict['mflux1'].copy();
                mwave2_in = indict['mwave2'].copy(); mflux2_in = indict['mflux2'].copy();

                for tag in taglist:

                    skiptag = False
                    par = parfittags[(tag == tags2)][0]

                    if np.isnan(par[0]):
                        skiptag = True
                    else:
                        watm = watm_in*(1.+par[2]/c)
                        satm = satm_in**par[3]

                        filename = f'{indict["inpath"]}/Cutout_{night}_{tag}.fits'
                        hdu2 = fits.open(filename)
                        tbdata2 = hdu2[jerp+1].data
                        flagtag  = np.array(tbdata2['ERRORFLAG'+str(order)])[0]
                        if flagtag == 1:
                            skiptag = True
                        else:
                            w  = np.array(tbdata2['WAVE_RAW'+str(order)],dtype=float)
                            x  = np.array(tbdata2['X'],dtype=float)
                            s2n = np.array(tbdata2['S2N'],dtype=float)
                            x    = x[w!=0]
                            s2n  = s2n[w!=0]
                            w    = w[w!=0]
                            npts = len(w)

                            if np.median(np.diff(mwave1_in)) > np.median(np.diff(mwave2_in)):
                                rebin2to1 = True; extra1 = 0.; extra2 = 10.;
                            else:
                                rebin2to1 = False; extra1 = 10.; extra2 = 0.;

                            mflux1 = mflux1_in[(mwave1_in > w[0]-30-extra1) & (mwave1_in < w[-1]+30+extra1)]
                            mwave1 = mwave1_in[(mwave1_in > w[0]-30-extra1) & (mwave1_in < w[-1]+30+extra1)]
                            mflux2 = mflux2_in[(mwave2_in > w[0]-30-extra2) & (mwave2_in < w[-1]+30+extra2)]
                            mwave2 = mwave2_in[(mwave2_in > w[0]-30-extra2) & (mwave2_in < w[-1]+30+extra2)]

                            # Apply velocity shifts and scale
                            mwave1 = mwave1*(1.+(indict['rv1'][i]-indict['bvcs'][night+tag])/c)
                            mflux1 = mflux1**par[1]
                            mwave2 = mwave2*(1.+(indict['rv2'][i]-indict['bvcs'][night+tag])/c)
                            mflux2 = mflux2**float(args.pow2)

                            wspot1,rspot1 = rotint(mwave1, mflux1, vsini)
                            wspot2,rspot2 = rotint(mwave2, mflux2, vsini)

                            rspot2   *=  fluxratio


                            if rebin2to1:
                                rspot2n = rebin_jv(wspot2,rspot2,wspot1,True)
                                rspot   = rspot1 + rspot2n
                                wspot   = wspot1.copy()
                            else:
                                rspot1n = rebin_jv(wspot1,rspot1,wspot2,True)
                                rspot   = rspot2 + rspot1n
                                wspot   = wspot2.copy()

                            rspot /= (1+fluxratio)

                            dstep0 = np.median(np.diff(wspot))
                            if dstep0 > 0.045:
                                pass
                                #logger.info(f'Stellar template resolution is ~{round(dstep0,4)} '
                                #                'Angstrom, leaving alone...')
                            else:
                                dstep = 0.045
                                nstep = int((wspot[-1]-wspot[0])/dstep)
                                mwave1 = np.linspace(wspot[0],wspot[-1],nstep)
                                mflux1 = rebin_jv(wspot,rspot,mwave1,False)
                                mwave0 = mwave1.copy(); mflux0 = mflux1.copy()
                                mwave0 = mwave0[1:-1]
                                mflux0 = mflux0[1:-1]
                                wspot = mwave0.copy()
                                rspot = mflux0.copy()

                                #logger.info(f'Stellar template resolution is ~{round(dstep0,4)} '
                                #                'Angstrom, rebinning to 0.045 Angstrom...')


                            rspot = rspot[(wspot > watm[0]) & (wspot < watm[-1])]
                            wspot = wspot[(wspot > watm[0]) & (wspot < watm[-1])]
                            satm2 = rebin_jv(watm, satm, wspot, True)
                            rspot *= satm2

                            #Find mean observed wavelength and create a telluric velocity scale
                            mnw = np.mean(w)
                            dw = (w[-1] - w[0])/(npts-1.)
                            vel = (wspot-mnw)/mnw*c

                            if len(x) != len(w):
                                sys.exit('uhoh')

                            fwhmraw = par[5] + par[13]*(x) + par[14]*(x**2)
                            spl = splrep(w, fwhmraw)
                            fwhm = splev(wspot,spl)

                            # Have IP extend as constant past wave bounds of data
                            fwhm[(wspot < w[0])]  = fwhm[(wspot >= w[0])][0]
                            fwhm[(wspot > w[-1])] = fwhm[(wspot <= w[-1])][-1]
                            if (np.min(fwhm) < 1) or (np.max(fwhm) > 8):
                                sys.exit('IP Error!')

                            #Handle instrumental broadening
                            vhwhm = dw*np.abs(fwhm)/mnw*c/2.
                            nsmod = macbro_dyn(vel, rspot, vhwhm)

                            #Rebin model to observed wavelength scale
                            smod = rebin_jv(wspot, nsmod ,w, False)

                            c2 = rebin_jv(a0contx,continuum,x,True)

                            # Apply continuum adjustment
                            cont = par[10] + par[11]*x + par[12]*(x**2) \
                                    + par[20]*(x**3) + par[21]*(x**4) \
                                    + par[22]*(x**5) + par[23]*(x**6)
                            if beam == 'A':
                                bucket = np.zeros_like(cont)
                                bucket[(x >= (par[15]-par[16]/2))         \
                                        & (x <= (par[15]+par[16]/2))] = par[17]
                                bucket[(x >= (par[15]+par[16]/2-par[18])) \
                                        & (x <= (par[15]+par[16]/2))] += par[19]
                                cont -= bucket

                            #matplotlib.use('Qt5Agg')

                            sout = smod*cont*c2

                            sout = sout[s2n > 1]
                            w = w[s2n > 1]
                            x = x[s2n > 1]
                            s2n = s2n[s2n > 1]

                            #print(night,order,w[0],w[-1])
                            #plt.plot(w,sout,color='black')
                            #plt.show()



                    if skiptag:
                        # Save results to fits file separately for each tight/loose dataset
                        c0 = fits.Column( name='FLAG'+str(order),         array=np.array([1]),    format='K')
                        cols  = fits.ColDefs([c0])

                    else:
                        # Save results to fits file separately for each tight/loose dataset
                        c0 = fits.Column( name='FLAG'+str(order),         array=np.array([0]),    format='K')
                        c1 = fits.Column( name='WAVE',         array=w,    format='D')
                        c2 = fits.Column( name='FLUX',            array=sout,       format='D')
                        c3 = fits.Column( name='S2N',         array=s2n,   format='D')
                        c4 = fits.Column( name='X',        array=x,  format='D')
                        cols  = fits.ColDefs([c0,c1,c2,c3,c4])

                    hdu_1 = fits.BinTableHDU.from_columns(cols)

                    if jerp == 0:
                        bleh = np.ones((3,3))
                        primary_hdu = fits.PrimaryHDU(bleh)
                        hdul        = fits.HDUList([primary_hdu,hdu_1])
                        hdul.writeto(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag),
                            overwrite=True)
                    else:
                        hh = fits.open(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag))
                        hh.append(hdu_1)
                        hh.writeto(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag),
                            overwrite=True)

            if skipnight:
                for tag in taglist:
                    # Save results to fits file separately for each tight/loose dataset
                    c0 = fits.Column( name='FLAG',         array=np.array([1]),    format='K')
                    cols  = fits.ColDefs([c0])

                    hdu_1 = fits.BinTableHDU.from_columns(cols)

                    if jerp == 0:
                        bleh = np.ones((3,3))
                        primary_hdu = fits.PrimaryHDU(bleh)
                        hdul        = fits.HDUList([primary_hdu,hdu_1])
                        hdul.writeto(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag),
                            overwrite=True)
                    else:
                        hh = fits.open(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag))
                        hh.append(hdu_1)
                        hh.writeto(
                            '{}/FakeData_{}_{}.fits'.format(
                                indict["outpath"],night,tag),
                            overwrite=True)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    args = _argparse_step6a()

    #-------------------------------------------------------------------------------

    # Check user input

    if args.fluxratio == '':
        sys.exit('ERROR: YOU MUST PROVIDE A GUESS FOR FLUX RATIO VALUE, "-f"')

    if (args.temperature == '') & (args.logg == ''):
        sys.exit('ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
                    'STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE '
                    'AVAILABLE TEMPLATES')
    #------------------------------

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')

    #------------------------------

    syntemp = os.listdir(f'./Engine/syn_template')

    if args.template.lower() == 'synthetic':
        #list of all syntheticstellar
        syntemp = [i for i in syntemp if i[:3] == 'syn']
        synT    = [ i.split('_')[2][1:]  for i in syntemp ]
        synlogg = [ i.split('_')[3][4:7] for i in syntemp ]
    elif args.template.lower() == 'phoenix':
        #list of all phoenix
        syntemp = [i for i in syntemp if i[:3] == 'PHO']
        synT    = [ i.split('-')[1][4:]  for i in syntemp ]
        synlogg = [ i.split('-')[2][:3] for i in syntemp ]
    else:
        synT = [args.temperature]; synlogg = [args.logg]

    if args.temperature not in synT:
        sys.exit(f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp" INPUT! '
                    f'{syntemp} AVALIABLE UNDER ./Engine/syn_template/')

    if args.logg not in synlogg:
        sys.exit(f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg" INPUT! {syntemp} '
                    'AVALIABLE UNDER ./Engine/syn_template/')


    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(args)

    #-------------------------------------------------------------------------------


    tbdataRV = Table.read(f'./Output/{args.targname}_{args.band}/RVresultsSummary_WithRV2_{args.run}.csv',format='csv')
    nightsRV = np.array(tbdataRV['NIGHT'],dtype=str)
    jd0      = np.array(tbdataRV['JD'],dtype=str)
    rvfinal1      = np.array(tbdataRV['RVfinal'],dtype=float)
    rvfinal2      = np.array(tbdataRV['RV2final'],dtype=float)
    vsini      = np.nanmean(np.array(tbdataRV['VSINI'],dtype=float))

    nightsFinal = nightsRV.copy()
    nights = nightsFinal.copy()

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    print('\n')

    #kinds = ['M2obs','M3obs','M5obs','M6obs','3500_4p0','3500_4p5','3000_4p0','3000_4p0_phx','3000_4p5_phx']
    #T2s   = [3600,3500,3200,3100,3500,3500,3000,3000,3000]

    name = f'FakeData_{args.pow2}pow2_fluxratio{args.fluxratio}_{args.template2}_{args.temperature2}_{args.logg2}_{args.B2}kG'

    inpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/Cutouts/'
    outpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/{name}'
    figpath = f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/{name}/figs'


    if not os.path.isdir(f'{outpath}'):
        os.mkdir(f'{outpath}')
    if not os.path.isdir(f'{outpath}/figs'):
        os.mkdir(f'{outpath}/figs')


    # Set up logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    logger.info(f'Writing output to {outpath}')

    logger.removeHandler(stream_hander)
    print('\n')


    watmtrash,satmtrash, mwave1, mflux1 = setup_templates(
        logger, args.template, args.band, np.int(args.temperature),
        np.float(args.logg), np.float(args.B)
        )

    watmtrash,satmtrash, mwave2, mflux2 = setup_templates(
        logger, args.template2, args.band, np.int(args.temperature2),
        np.float(args.logg2), np.float(args.B2)
        )


    indict = {'inpath':inpath,'outpath':outpath,'mwave1':mwave1,'mwave2':mwave2,
              'mflux1':mflux1,'mflux2':mflux2,"vsini":vsini,
              'bvcs':bvcs,'tagsA':tagsA,'tagsB':tagsB,'nights':nights,
              'orders':orders,'rv1':rvfinal1,'rv2':rvfinal2}

    rawbox = fits.open(f'./Output/{args.targname}_{args.band}/RV_results_{args.run}/RVresultsRawBox.fits')

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        print('Working on order {} ({:02d}/{:02d})'.format(
                orders[jerp], int(jerp+1), len(orders)
                ))

        order = orders[jerp]

        boxdata = rawbox[jerp+1].data
        nightsbox0 = np.array(boxdata['NIGHT'+str(orders[jerp])])
        parfitbox = np.array(boxdata['PARFIT'+str(orders[jerp])])
        tagbox    = np.array(boxdata['TAG'+str(orders[jerp])])

        #makespec(nightsbox0,parfitbox,tagbox,args,indict,jerp,0)

        func = partial(makespec, nightsbox0,parfitbox,tagbox,args,indict,jerp)
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)



    print('####################################################################################')
