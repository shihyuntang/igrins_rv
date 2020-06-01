import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *

from Engine.IO_AB     import setup_templates, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips     import basicclip_above
from Engine.contfit   import A0cont
from Engine.classes   import fitobjs,inparamsA0
from Engine.macbro_dynamic    import macbro_dyn
from Engine.rebin_jv  import rebin_jv
from Engine.rotint    import rotint
from Engine.Telfitter import telfitter

from Engine.opt import optimizer, fmod
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def outplotter(parfit,fitobj,title,debug):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
    axes.plot(w,fitobj.s, '-k',  lw=0.5, label='data',alpha=.6)
    axes.plot(w,fit,      '--r', lw=0.5, label='model',alpha=.6)

    axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
    axes.set_title(title,   size=5, style='normal' , family='sans-serif' )
    axes.set_ylabel(r'Normalized Flux',   size=5, style='normal' , family='sans-serif' )
    axes.set_xlabel(r'Wavelength [$\AA$]',       size=5, style='normal' , family='sans-serif' )
    axes.legend(fontsize=4, edgecolor='white')
    if debug == 0:
        fig.savefig('{}/figs/{}.png'.format(inparam.outpath, title), bbox_inches='tight', format='png', overwrite=True)
    elif debug == 1:
        fig.savefig('./Temp/Debug/{}/{}.png'.format(args.targname, title), bbox_inches='tight', format='png', overwrite=True)


def DataPrep(args):
    # Find all nights of observations of target in master log
    master_log_fh = './Engine/IGRINS_MASTERLOG.csv'
    master_log    = pd.read_csv(master_log_fh)

    star_files    = master_log[(master_log['OBJNAME'].str.contains(args.targname, regex=True, na=False)) & (master_log['OBJTYPE'].str.contains('TAR', regex=True, na=False))]
    allnights     = np.array(master_log['CIVIL'],dtype='str')

    # If star input not found in Masterlog, try putting a space in its name somewhere
    n = 1
    while len(star_files['CIVIL']) == 0:
        starnew = args.targname[:n]+' '+args.targname[n:]
        star_files = master_log[(master_log['OBJNAME'].str.contains(starnew, regex=True, na=False)) & (master_log['OBJTYPE'].str.contains('TAR', regex=True, na=False))]
        n += 1
        if n == len(args.targname):
            sys.exit('TARGET NAME NOT FOUND IN CATALOG - CHECK INPUT!')

    #---------------------------------------------------------------------------
    ## Collect target star information
    fileT = open('./Temp/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), 'w')
    fileT.write('night beam tag mjd facility airmass bvc\n')

    nightsT = [];
    for x in range(len(star_files['CIVIL'])):
        night    = str(  np.array(star_files['CIVIL'])[x]     )
        frame    = str(  np.array(star_files['FRAMETYPE'])[x] )
        tag0     = int(  np.array(star_files['FILENUMBER'])[x])
        airmass  = float(np.array(star_files['AM'])[x]        )
        BVCfile  = float(np.array(star_files['BVC'])[x]       )
        facility = str(  np.array(star_files['FACILITY'])[x]  )
        tag = '{:04d}'.format(tag0)

        try:
            hdulist = fits.open('{}{}/{}/SDC{}_{}_{}.spec.fits'.format(inpath, night, frame, args.band, night, tag))

        except FileNotFoundError:
            continue

        head = hdulist[0].header
        if head['OBSERVAT'] == 'Lowell Observatory':
            obs = 'DCT'
        elif head['OBSERVAT'] == 'McDonald':
            obs = 'McD'
        else:
            print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}, MUST EDIT CODE TO INCLUDE THIS OPTION'.format( head['OBSERVAT'] ))

        try:
            time_midpoint = np.mean([float(head['JD-OBS']),float(head['JD-END'])])
        except KeyError:
            l0 = []
            for nm in ['DATE-OBS','DATE-END']:
                tt1 = head[nm].split('-')
                t1 = Time(tt1[0]+'-'+tt1[1]+'-'+tt1[2]+' '+tt1[3],format='iso')
                l0.append(t1.jd)
            time_midpoint = np.mean(l0)

        mjd = time_midpoint;
        fileT.write(night+' '+frame+' '+str(tag)+' '+str(mjd)+' '+str(facility)+' '+str(airmass)+' '+str(BVCfile))
        fileT.write('\n')

        nightsT.append(night)

    fileT.close()


    #---------------------------------------------------------------------------
    ## Now collect A0 information
    fileA0 = open('./Temp/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), 'w')
    fileA0.write('night tag humid temp zd press obs airmass\n')
    noA0nights = []

    for night0 in np.unique(nightsT):
        night    = str(night0)
        am_stars = star_files['AM'][(np.array(star_files['CIVIL'], dtype=str) == night)]
        am_star  = float(am_stars.values[0])

        std_files = master_log[(allnights == night) & (master_log['OBJTYPE'].str.contains('STD', regex=True, na=False))]
        ams = np.array(std_files['AM'].values,dtype='float')

        # First check whether any A0s observed that night
        if len(ams) == 0:
            noA0nights.append(night)
            tagA = 'NA'; humid = 'NA'; temp = 'NA'; zd = 'NA'; press = 'NA'; obs = 'NA'; AM = 'NA';
        else:
            # Then check whether any A0s files for that night outputted by reduction pipeline.
            # If not, Joe either didn't have the data for them or didn't copy them over.
            anyK = False
            subpath        = '{}std/{}/AB/'.format(inpath, night)
            fullpathprefix = '{}SDC{}_{}_'.format(subpath, args.band, night)

            onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
            for f in onlyfiles:
                q = re.split('_', f)
                if q[0] != 'SDC{}'.format(args.band):
                    continue
                anyK = True

            if anyK == False:
                tagA = 'NA'; humid = 'NA'; temp = 'NA'; zd = 'NA'; press = 'NA'; obs = 'NA'; AM = 'NA';
                noA0nights.append(night)
            else:
                stdname = std_files['OBJNAME'][abs(ams-am_star) == min(abs(ams-am_star))].values[0]
                fileno0 = int(std_files['FILENUMBER'][abs(ams-am_star) == min(abs(ams-am_star))].values[0])
                names   = np.array(std_files['OBJNAME'])
                tagA0s  = np.array(std_files['FILENUMBER'][(names == stdname)].values)
                facs    = np.array(std_files['FACILITY'][(names == stdname)].values)
                am0s    = ams[(names == stdname)]

                firsts = []
                for k, g in groupby(enumerate(tagA0s), lambda ix : ix[0] - ix[1]):
                    q = list(map(itemgetter(1), g))
                    firsts.append(q[0])
                firsts = np.array(firsts)

                tagA0 = firsts[abs(firsts-fileno0) == min(abs(firsts-fileno0))][0]
                am0   = am0s[(tagA0s == tagA0)][0]
                fac   = np.array(facs)[(tagA0s == tagA0)][0]

                if abs(am0-am_star) > 0.2:
                    print(night,stdname,am_star,am0,tagA0)
                    sys.exit('WARNING, A0 AIRMASS FOR NIGHT {} HAS A DIFFERENCE LARGER THAN 0.2 FROM TARGET!'.format(night))

                tagA = '{:04d}'.format(tagA0)
                subpath = '{}std/{}/AB/SDC{}_{}_{}.spec.fits'.format(inpath, night, args.band, night, tagA)

                try:
                    hdulist = fits.open(subpath)
                    head    = hdulist[0].header

                except FileNotFoundError:
                    # If best airmass match A0 for night not found, check if Joe chose a different A0 instead
                    subpath        = '{}std/{}/AB/'.format(inpath, night)
                    fullpathprefix = '{}SDC{}_{}_'.format(subpath, args.band, night)

                    onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
                    for f in onlyfiles:
                        q = re.split('_', f)
                        if q[0] != 'SDC{}'.format(args.band):
                            continue
                        qr = re.split('\.', q[2])
                        tagA = qr[0]

                        hdulist = fits.open(fullpathprefix+tagA+'.spec.fits')
                        head = hdulist[0].header
                        am0 = np.mean([float(head['AMSTART']),float(head['AMEND'])])

                if head['OBSERVAT'] == 'Lowell Observatory':
                    obs = 'DCT'
                elif head['OBSERVAT'] == 'McDonald':
                    obs = 'McD'
                else:
                    print('EXPECTED LOWELL OR MCDONALD OBSERVATORY, GOT {}, MUST EDIT CODE TO INCLUDE THIS OPTION'.format( head['OBSERVAT'] ))

                AM = str(am0)
                try:
                    humid = float(head['HUMIDITY'])
                    temp  = float(head['AIRTEMP'])
                    press = float(head['BARPRESS'])
                    zd = np.mean([float(head['ZDSTART']),float(head['ZDEND'])])
                except ValueError: # McDonald headers don't have these quantities :(
                    humid = 'NOINFO'; temp = 'NOINFO'; press = 'NOINFO'; zd = 'NOINFO';

        fileA0.write(night+' '+str(tagA)+' '+str(humid)+' '+str(temp)+' '+str(zd)+' '+str(press)+' '+str(obs)+' '+str(AM))
        fileA0.write('\n')

    fileA0.close()

    print('No reduced A0s found for following nights:')
    for n in noA0nights:
        print(n)
    print('To achieve highest precision, this pipeline defaults to not analyzing target spectra for these nights.\n')


def MPinst(args, chunk_ind, i):
    nights = inparam.nights

    order = inparam.xbounddict[chunk_ind]
    night = str(nights[i])
    firstorder = int(inparam.xbounddict[0])

    print('Working on chunk {}/{}, night {} {}/{} ...'.format(chunk_ind+1,
                                                             len(inparam.xbounddict),
                                                             night,
                                                             i+1,
                                                             len(inparam.nights)) )

    # if int(night) < 20180401 or int(night) > 20190531:
    #     IPpars = inparam.ips_tightmount_pars[order]
    # else:
    #     IPpars = inparam.ips_loosemount_pars[order]

    IPpars = np.array([ 0,  0,  3.3])

    if args.band=='K':
        if order==11:
            bound_cut = [200, 100]
        elif order==12:
            bound_cut = [900, 300]
        elif order==13:
            bound_cut = [200, 400]
        elif order==14:
            bound_cut = [150, 300]
        else:
            bound_cut = [150, 100]
    elif args.band=='H':
        if order==10:
            bound_cut = [250, 150]#ok
        elif order==11:
            bound_cut = [600, 150]
        elif order==13:
            bound_cut = [200, 600]#ok
        elif order==14:
            bound_cut = [700, 100]
        elif order==16:
            bound_cut = [400, 100]
        elif order==17:
            bound_cut = [1000, 100]
        elif order==20:
            bound_cut = [500, 150]
        elif (order==7) or (order==8) or (order==9) or (order==12) or (order==15) or (order==18) or (order==19):
            bound_cut = [500, 500]
        else:
            bound_cut = [150, 100]

    ### Load relevant A0 spectrum
    x, a0wavelist, a0fluxlist, u = init_fitsread(inparam.inpath,
                                                'A0',
                                                'separate',
                                                night,
                                                order,
                                                '{:04d}'.format(int(inparam.tags[night])),
                                                args.band,
                                                bound_cut)

    nzones = 12
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(x,a0fluxlist,nzones);
    a0u        = basicclip_above(u,a0fluxlist,nzones);     a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);

    # do twice?
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(a0x,a0fluxlist,nzones);
    a0u        = basicclip_above(a0u,a0fluxlist,nzones);   a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);

    # Normalize
    a0fluxlist /= np.median(a0fluxlist)

    # Compute rough blaze fn estimate
    continuum    = A0cont(a0wavelist,a0fluxlist,night,order)
    a0masterwave = a0wavelist.copy()
    a0masterwave *= 1e4

    # Trim stellar template to relevant wavelength range
    mwave_in, mflux_in = stellarmodel_setup(a0wavelist, inparam.mwave0, inparam.mflux0)

    # Trim telluric template to relevant wavelength range
    satm_in = inparam.satm[(inparam.watm > min(a0wavelist)*1e4 - 11) & (inparam.watm < max(a0wavelist)*1e4 + 11)]
    watm_in = inparam.watm[(inparam.watm > min(a0wavelist)*1e4 - 11) & (inparam.watm < max(a0wavelist)*1e4 + 11)]

    f = np.polyfit(a0x, a0wavelist, 3)
    par9in = f[0]*1e4; par8in = f[1]*1e4; par7in = f[2]*1e4; par6in = f[3]*1e4;

    ### Initialize parameter array for optimization as well as half-range values for each parameter during
    ### the various steps of the optimization.
    ### Many of the parameters initialized here will be changed throughout the code before optimization and
    ### in between optimization steps.
    parA0 = np.array([0.0,           # 0: The shift of the sunspot spectrum (km/s)
                      0.0,           # 1: The scale factor for the sunspot spectrum
                      0.0,           # 2: The shift of the telluric spectrum (km/s)
                      1.0,           # 3: The scale factor for the telluric spectrum
                      0.0,           # 4: vsini (km/s)
                      IPpars[2],     # 5: The instrumental resolution (FWHM) in pixels
                      par6in,        # 6: Wavelength 0-pt
                      par7in,        # 7: Wavelength linear component
                      par8in,        # 8: Wavelength quadratic component
                      par9in,        # 9: Wavelength cubic component
                      1.0,           #10: Continuum zero point
                      0.,            #11: Continuum linear component
                      0.,            #12: Continuum quadratic component
                      IPpars[1],     #13: IP linear component
                      IPpars[0]])    #14: IP quadratic component

    # Cut target spectrum to be within telluric template wavelengths (should be unncessary)
    a0fluxlist = a0fluxlist[(a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0u        = a0u[       (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0x        = a0x[       (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    continuum  = continuum[ (a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]
    a0wavelist = a0wavelist[(a0wavelist*1e4 > min(watm_in)+5) & (a0wavelist*1e4 < max(watm_in)-5)]

    # Define main spectrum parameters
    s = a0fluxlist.copy(); x = a0x.copy(); u = a0u.copy();

    # Collect all fit variables into one class
    fitobj = fitobjs(s, x, u, continuum, watm_in, satm_in, mflux_in, mwave_in)

    # Arrays defining parameter variations during optimization steps
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.0,  0,    1e7, 1, 1, 0,    0])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 5e-5, 1e-7, 0,   0, 0, 0,    0])
    dpar      = np.array([0.0, 0.0, 5.0, 3.0, 0.0, 0.5, 0.0,   0.0, 0.0,  0,    1e4, 1, 1, 1e-2, 1e-5])
    dpar_st   = np.array([0.0, 0.0, 5.0, 3.0, 0.0, 0.0, 0.0,   0.0, 0.0,  0,    1e4, 1, 1, 0,    0])
    dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0,   0.0, 0.0,  0,    1e4, 1, 1, 1e-2, 1e-5])

    # For every pre-Telfit spectral fit, first fit just template strength/rv/continuum, then just wavelength soln, then template/continuum again, then ip,
    # then finally wavelength. Normally would fit for all but wavelength at the end, but there's no need for the pre-Telfit fit, since all we want
    # is the wavelength solution.


    optimize = True
    par_in = parA0.copy()
    parfit_1 = optimizer(par_in,   dpar_st,   fitobj, optimize)
    parfit_2 = optimizer(parfit_1, dpar_wave, fitobj, optimize)
    parfit_3 = optimizer(parfit_2, dpar_st,   fitobj, optimize)
    parfit_4 = optimizer(parfit_3, dpar_ip,   fitobj, optimize)
    parfit = optimizer(parfit_4, dpar_wave, fitobj, optimize)

    # if inparam.plotfigs == True:
    #     outplotter(parfit, fitobj, '{}_{}_1'.format(order,night), 0)

    # Get fitted wavelength solution
    a0w_out_fit = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)

    # Trim stellar template to relevant wavelength range
    mwave_in,mflux_in = stellarmodel_setup(a0w_out_fit/1e4,inparam.mwave0,inparam.mflux0)

    # Using this new wavelength solution, get Telfit'd telluric template, parameters of that best fit, and blaze fn best fit
    watm1, satm1, telfitparnames, telfitpars, a0contwave, continuum = telfitter(a0w_out_fit,a0fluxlist,a0u,inparam,night,order,args)

    if len(watm1) == 1: # If Telfit encountered error mentioned in Telfitter.py, skip night/order combo

        print('TELFIT ENCOUNTERED CRITICAL ERROR, ORDER '+str(order)+' NIGHT '+str(night))

        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name='ERRORFLAG'+str(order),array=np.array([1]),format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)
        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto(inparam.outpath+'/'+night+'A0_treated.fits')
        else:
            hh = fits.open(inparam.outpath+'/'+night+'A0_treated.fits')
            hh.append(hdu_1)
            hh.writeto(inparam.outpath+'/'+night+'A0_treated.fits', overwrite=True)

    else:

        a0contwave /= 1e4
        continuum = rebin_jv(a0contwave,continuum,a0wavelist,False)

        # Fit whole A0 again to get even better wave soln to use for a0contwave and tweak blaze fn fit as
        # needed with quadratic adjustment
        fitobj = fitobjs(s, x, u, continuum,watm1,satm1,mflux_in,mwave_in)

        parfit_1 = optimizer(par_in,   dpar_st,   fitobj, optimize)
        parfit_2 = optimizer(parfit_1, dpar_wave, fitobj, optimize)
        parfit_3 = optimizer(parfit_2, dpar_st,   fitobj, optimize)
        parfit_4 = optimizer(parfit_3, dpar_ip,   fitobj, optimize)
        parfit_5 = optimizer(parfit_4, dpar_wave, fitobj, optimize)
        parfit   = optimizer(parfit_5, dpar,      fitobj, optimize)

        if inparam.plotfigs == True:
            fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
            axes.plot(fitobj.x,
                      parfit[5] + parfit[13]*(fitobj.x) + parfit[14]*(fitobj.x**2),
                      '-k', lw=0.5, label='data', alpha=.6)

            axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
            axes.set_ylabel(r'IP',   size=5, style='normal' , family='sans-serif' )
            axes.set_xlabel('Pixel Number',       size=5, style='normal' , family='sans-serif' )
            axes.legend(fontsize=4, edgecolor='white')
            fig.savefig(inparam.outpath+'/figs/{}_{}_ip.png'.format(order,night), bbox_inches='tight', format='png', overwrite=True)

            outplotter(parfit,fitobj,'{}_{}_post_parfit'.format(order,night), 0)

        if args.debug == True:
            outplotter(parfit_1,fitobj,'{}_{}_post_parfit_1'.format(order,night), 1)
            outplotter(parfit_2,fitobj,'{}_{}_post_parfit_2'.format(order,night), 1)
            outplotter(parfit_3,fitobj,'{}_{}_post_parfit_3'.format(order,night), 1)
            outplotter(parfit_4,fitobj,'{}_{}_post_parfit_4'.format(order,night), 1)
            outplotter(parfit_5,fitobj,'{}_{}_post_parfit_5'.format(order,night), 1)
            outplotter(parfit  ,fitobj,'{}_{}_post_parfit'.format(order,night), 1)

        a0w_out  = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)
        cont_adj = parfit[10] + parfit[11]*x + parfit[12]*(x**2.)

        continuum /= np.median(continuum)
        cont_save = continuum*cont_adj

        # Write out table to fits file with errorflag = 0
        c0 = fits.Column(name='ERRORFLAG'+str(order),      array=np.array([0]),            format='K')
        c1 = fits.Column(name='WAVE'+str(order),           array=a0w_out,                  format='D')
        c2 = fits.Column(name='BLAZE'+str(order),          array=cont_save,                format='D')
        c3 = fits.Column(name='X'+str(order),              array=a0x,                      format='D')
        c4 = fits.Column(name='INTENS'+str(order),         array=a0fluxlist,               format='D')
        c5 = fits.Column(name='SIGMA'+str(order),          array=a0u,                      format='D')
        c6 = fits.Column(name='WATM'+str(order),           array=watm1,                    format='D')
        c7 = fits.Column(name='SATM'+str(order),           array=satm1,                    format='D')
        c8 = fits.Column(name='TELFITPARNAMES'+str(order), array=np.array(telfitparnames), format='8A')
        c9 = fits.Column(name='TELFITPARS'+str(order),     array=telfitpars,               format='D')
        c10 = fits.Column(name='PARFIT',                    array=parfit,                   format='D')
        cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto(inparam.outpath+'/'+night+'A0_treated.fits',overwrite=True)
        else:
            hh = fits.open(inparam.outpath+'/'+night+'A0_treated.fits')
            hh.append(hdu_1)
            hh.writeto(inparam.outpath+'/'+night+'A0_treated.fits',overwrite=True)


def mp_run(args, Nthreads, jerp, nights):
    pool = mp.Pool(processes = Nthreads)
    func = partial(MPinst, args, jerp)
    outs = pool.map(func, np.arange(len(nights)))
    pool.close()
    pool.join()


#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                                     prog        = 'IGRINS Spectra Radial Velocity Pipeline',
                                     description = '''
                                     This is a pipeline that helps you to extract radial velocity \n
                                     from IGRINS spectra. \n
                                     ''',
                                     epilog = "Contact authors: asa.stahl@rice.edu; sytang@lowell.edu")
    parser.add_argument("targname",                          action="store",
                        help="Enter your *target name",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",        action="store_true",
                        help="If sets, will generate plots")
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If sets, will generate files and plots under ./Temp/Debug for debug")
    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()


    cdbs_loc = '~/cdbs/'
    inpath     = './Input_Data/{}/'.format(args.targname)
    if args.debug:
        try:
            os.listdir('./Temp/Debug/{}/'.format(args.targname))
        except OSError:
            os.mkdir('./Temp/Debug/{}/'.format(args.targname))

    targname   = args.targname
    Nthreads   = args.Nthreads

    #------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('Data Preparation for {} (1/3)...'.format(args.targname))
    time.sleep(1)
    DataPrep(args)

    print('Data Preparation Done!')

    print('\n')

    #time.sleep(5)
    print('---------------------------------------------------------------')
    #------------
    print('A0 Fitting using TelFit for {} (3/3)...'.format(args.targname))
    print('This will take a while..........')
    print('\n')
    curdir = os.getcwd()


    ## Collect relevant file information from Predata files
    A0data = Table.read('./Temp/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), format='ascii')

    ind    = [i != 'NA' for i in A0data['humid']]
    humids = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['humid'])}
    tags   = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['tag'])}
    obs    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['obs'])}
    temps  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['temp'])}
    zds    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['zd'])}
    press  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['press'])}
    nightsFinal = np.array(list(sorted(set(A0data[ind]['night']))))

    nightsFinal = np.append( nightsFinal[1:10:2], nightsFinal[-10:-1:2])

    filesndirs = os.listdir('./A0_Fits/')
    name = 'A0_Fits_'+ args.targname
    if name not in filesndirs:
        os.mkdir('./A0_Fits/{}'.format(name) )
        os.mkdir('./A0_Fits/{}/figs'.format(name) )
    outpath = './A0_Fits/' + name

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates()

    if args.band == 'K':
        orders = np.append(np.arange(2, 9), np.array([10, 11, 12, 13, 14, 16]))
        #orders = [3]
    elif args.band=='H':
        orders = np.append(np.arange(2, 7), np.array([10, 11, 13, 14, 16, 17, 20, 21, 22]))

    inparam = inparamsA0(inpath,outpath,args.plotfigs,tags,nightsFinal,humids,
                         temps,zds,press,obs,watm,satm,mwave0,mflux0,cdbs_loc,orders)

    print('Processing Orders: {}'.format(orders))
    for jerp in range(len(orders)): # Iterate over regions
        outs = mp_run(args, args.Nthreads, jerp, nightsFinal)

    print('\n')
    print('A0 Fitting Done!')

    end_time = datetime.now()
    print('A0 Fitting using TelFit finished, Duration: {}'.format(end_time - start_time))
    print('You can start to run main_step2.py for RV initial guess')
    print('###############################################################')
    print('\n')