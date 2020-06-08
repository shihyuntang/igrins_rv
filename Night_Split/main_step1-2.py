import sys
sys.path.append("..") # Adds higher directory to python modules path.

from Engine.importmodule import *

from Engine.IO_AB     import setup_templates, init_fitsread, stellarmodel_setup, setup_outdir
from Engine.clips     import basicclip_above
from Engine.contfit   import A0cont
from Engine.classes   import fitobjs,inparamsA0,orderdict_cla
from Engine.macbro    import macbro
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
    axes.plot(w,fitobj.s, '-',  c = 'k',        lw=0.5, label='data',  alpha=.6)
    axes.plot(w,fit,      '--', c = 'tab:red',  lw=0.5, label='model', alpha=.6)

    axes.set_title( title,                 size=5, style='normal', family='sans-serif')
    axes.set_ylabel(r'Normalized Flux',    size=5, style='normal', family='sans-serif')
    axes.set_xlabel(r'Wavelength [$\AA$]', size=5, style='normal', family='sans-serif')

    axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
    axes.legend(fontsize=4, edgecolor='white')
    if debug == 0:
        fig.savefig('{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title), bbox_inches='tight', format='png', overwrite=True)
    elif debug == 1:
        fig.savefig('../Temp/Debug/{}/{}.png'.format(args.targname, title), bbox_inches='tight', format='png', overwrite=True)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# def MPinst(order, i):
def MPinst(args, chunk_ind, orders, i):
    order = orders[chunk_ind]
    night = str(inparam.nights[i])
    firstorder = orders[0]
    print('Working on order {}/{}, night {} {}/{} PID:{}...'.format(chunk_ind+1,
                                                             len(orders),
                                                             night,
                                                             i+1,
                                                             len(inparam.nights),
                                                             mp.current_process().pid) )

    if int(night) < 20180401 or int(night) > 20190531:
        IPpars = inparam.ips_tightmount_pars[args.band][order]
    else:
        IPpars = inparam.ips_loosemount_pars[args.band][order]
#-------------------------------------------------------------------------------
    ### Load relevant A0 spectrum

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

    x, a0wavelist, a0fluxlist, u = init_fitsread(inparam.inpath,
                                                'A0',
                                                'separate',
                                                night,
                                                order,
                                                '{:04d}'.format(int(inparam.tags[night])),
                                                args.band,
                                                bound_cut)
#-------------------------------------------------------------------------------
    nzones = 12
    a0wavelist = basicclip_above(a0wavelist,a0fluxlist,nzones);   a0x = basicclip_above(x,a0fluxlist,nzones);
    a0u        = basicclip_above(u,a0fluxlist,nzones);     a0fluxlist = basicclip_above(a0fluxlist,a0fluxlist,nzones);
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
                      IPpars[0],     #14: IP quadratic component
                      0.0])          #15: Differential Rotation Coefficient

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
    dpar_cont = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,        0.,   1e7, 1, 1, 0,    0, 0])
    dpar_wave = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,  10.0, 5.00000e-5, 1e-7, 0,   0, 0, 0,    0, 0])
    dpar      = np.array([0.0, 0.0, 5.0, 3.0, 0.0, 0.5, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0])
    dpar_st   = np.array([0.0, 0.0, 5.0, 3.0, 0.0, 0.0, 0.0,   0.0,  0.0,        0,    1e4, 1, 1, 0,    0, 0])
    dpar_ip   = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0,   0.0,  0.0,        0,    0,   0, 0, 0,    0, 0])

#-------------------------------------------------------------------------------
    # For every pre-Telfit spectral fit, first fit just template strength/rv/continuum, then just wavelength soln, then template/continuum again, then ip,
    # then finally wavelength. Normally would fit for all but wavelength at the end, but there's no need for the pre-Telfit fit, since all we want
    # is the wavelength solution.

    optimize = True
    par_in = par.copy()
    hardbounds = [par_in[4] -dpar[4],   par_in[4]+dpar[4],
                  par_in[5] -dpar[5],   par_in[5]+dpar[5],
                  par_in[15]-dpar[15], par_in[15]+dpar[15]]
    if hardbounds[0] < 0:
        hardbounds[0] = 0
    if hardbounds[3] < 0:
        hardbounds[3] = 1

#        if args.plotfigs == True:#
#            outplotter(targname,par_in,fitobj,'{}_{}_{}_1'.format(label,night,tag))

    parfit_1 = optimizer(par_in,   dpar_st,   hardbounds,fitobj,optimize)
    parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds,fitobj,optimize)
    parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds,fitobj,optimize)
    parfit_4 = optimizer(parfit_3, dpar,      hardbounds,fitobj,optimize)
    parfit = optimizer(parfit_4,   dpar_wave, hardbounds,fitobj,optimize)

    # if inparam.plotfigs == True:
    #     outplotter(parfit, fitobj, '{}_{}_1'.format(label,night), 0)
#-------------------------------------------------------------------------------
    # Get fitted wavelength solution
    a0w_out_fit = parfit[6] + parfit[7]*x + parfit[8]*(x**2.) + parfit[9]*(x**3.)

    # Trim stellar template to relevant wavelength range
    mwave_in,mflux_in = stellarmodel_setup(a0w_out_fit/1e4, inparam.mwave0, inparam.mflux0)

    # Using this new wavelength solution, get Telfit'd telluric template, parameters of that best fit, and blaze fn best fit
    watm1, satm1, telfitparnames, telfitpars, a0contwave, continuum = telfitter(a0w_out_fit,a0fluxlist,a0u,inparam,night,order,args)
    # watm1, satm1 from Telfit, fack one.
#-------------------------------------------------------------------------------
    if len(watm1) == 1: # If Telfit encountered error mentioned in Telfitter.py, skip night/order combo
        print('TELFIT ENCOUNTERED CRITICAL ERROR, ORDER '+str(order)+' NIGHT '+str(night))
        # Write out table to fits header with errorflag = 1
        c0    = fits.Column(name='ERRORFLAG'+str(order), array=np.array([1]), format='K')
        cols  = fits.ColDefs([c0])
        hdu_1 = fits.BinTableHDU.from_columns(cols)

        if order == firstorder: # If first time writing fits file, make up filler primary hdu
            bleh = np.ones((3,3))
            primary_hdu = fits.PrimaryHDU(bleh)
            hdul = fits.HDUList([primary_hdu,hdu_1])
            hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band) )
        else:
            hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)
    else:

        a0contwave /= 1e4
        continuum = rebin_jv(a0contwave,continuum,a0wavelist,False)

        # Fit whole A0 again to get even better wave soln to use for a0contwave and tweak blaze fn fit as
        # needed with quadratic adjustment
        fitobj = fitobjs(s, x, u, continuum,watm1,satm1,mflux_in,mwave_in)

        parfit_1 = optimizer(par_in,   dpar_st,   hardbounds, fitobj, optimize)
        parfit_2 = optimizer(parfit_1, dpar_wave, hardbounds, fitobj, optimize)
        parfit_3 = optimizer(parfit_2, dpar_st,   hardbounds, fitobj, optimize)
        parfit_4 = optimizer(parfit_3, dpar_wave, hardbounds, fitobj, optimize)
        parfit   = optimizer(parfit_4, dpar,      hardbounds, fitobj, optimize)

        if inparam.plotfigs == True:
            fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
            axes.plot(fitobj.x,
                      parfit[5] + parfit[13]*(fitobj.x) + parfit[14]*(fitobj.x**2),
                      '-k', lw=0.5, label='data', alpha=.6)

            axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
            axes.set_ylabel(r'Normalized Flux',   size=5, style='normal' , family='sans-serif' )
            axes.set_xlabel('Wavelength',       size=5, style='normal' , family='sans-serif' )
            axes.legend(fontsize=4, edgecolor='white')
            fig.savefig('{}/figs_{}/IP_{}_{}.png'.format(inparam.outpath, args.band, order, night), bbox_inches='tight', format='png', overwrite=True)

            outplotter(parfit,fitobj,'Post_parfit_{}_{}'.format(order,night), 0)

        if args.debug == True:
            outplotter(parfit_1,fitobj,'Post_parfit_1_{}_{}'.format(order,night), 1)
            outplotter(parfit_2,fitobj,'Post_parfit_2_{}_{}'.format(order,night), 1)
            outplotter(parfit_3,fitobj,'Post_parfit_3_{}_{}'.format(order,night), 1)
            outplotter(parfit_4,fitobj,'Post_parfit_4_{}_{}'.format(order,night), 1)
            outplotter(parfit  ,fitobj,'Post_parfit_{}_{}'.format(order,night), 1)

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
            hdul.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band) ,overwrite=True)
        else:
            hh = fits.open('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band))
            hh.append(hdu_1)
            hh.writeto('{}/{}A0_treated_{}.fits'.format(inparam.outpath, night, args.band), overwrite=True)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def mp_run(args, Nthreads, jerp, orders, nights):
    pool = mp.Pool(processes = Nthreads)
    func = partial(MPinst, args, jerp, orders)
    outs = pool.map(func, np.arange(len(nights)))
    pool.close()
    pool.join()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def use_w(args):
    try:
        bounddata = Table.read('../Input_Data/Use_w/WaveRegions_{}_{}.csv'.format(args.WRegion, args.band), format='csv')
    except IOError:
        sys.exit('WaveRegions FILE ../Input_Data/Use_w/WaveRegions_{}_{}.csv NOT FOUND!'.format(args.WRegion, args.band))
    wavesols = pd.read_csv('../Input_Data/Use_w/WaveSolns_{}.csv'.format(args.band))
#-------------------------------------------------------------------------------
    filew = open('../Input_Data/Use_w/XRegions_{}_{}.csv'.format(args.WRegion, args.band),'w')
    filew.write('label, start,  end\n')

    m_order  = np.array(bounddata['order'])
    starts   = np.array(bounddata['start'])
    ends     = np.array(bounddata['end'])
    ords     = list( sorted(orderdict_cla().orderdict[args.band].keys()) )

    Ostarts  = [orderdict_cla().orderdict[args.band][k][0] for k in ords]
    Oends    = [orderdict_cla().orderdict[args.band][k][1] for k in ords]
    labels   = []
    for i in range(len(starts)):
        indS = list(np.where((starts[i] > Ostarts) & (starts[i] < Oends))[0])
        indE = list(np.where((ends[i]   > Ostarts) & (ends[i]   < Oends))[0])
        indboth = indS
        indboth.extend(x for x in indE if x not in indboth)
        for ind in indboth:
            wavebounds = [max([starts[i],Ostarts[ind]]),min([ends[i],Oends[ind]])]
            wO   = wavesols['w'+str(ords[ind])]
            pixO = wavesols['x'+str(ords[ind])];
            pix  = [pixO[(np.argmin(abs(wO-wavebounds[k])))] for k in [0,1]]

            p = 1
            if ords[ind] == m_order[i]:
                while 1 == 1:
                    lab = '{}-{}'.format(ords[ind],p)
                    if lab not in labels:
                        filew.write('{}, {}, {}\n'.format(lab,pix[0],pix[1]))
                        labels.append(lab)
                        break
                    else:
                        p += 1
    filew.close()

#-------------------------------------------------------------------------------
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
                        help="Enter your *target name, no space",            type=str)
    parser.add_argument("-HorK",    dest="band",             action="store",
                        help="Which band to process? H or K?. Default = K",
                        type=str,   default='K')
    parser.add_argument("-Wr",      dest="WRegion",          action="store",
                        help="Which ./Input_Data/Use_w/WaveRegions_X to use, Default X = 0",
                        type=int,   default=int(0))
    parser.add_argument("-AM",      dest="AM_cut",           action="store",
                        help="AirMass difference allowed between TAR and STD (A0) stars. Default X = 0.25 ",
                        type=str,   default='0.25')

    parser.add_argument('-c',       dest="Nthreads",         action="store",
                        help="Number of cpu (threads) to use, default is 1/2 of avalible ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
                        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot',    dest="plotfigs",        action="store_true",
                        help="If sets, will generate basic plots of A0 model fitting under ./A0_Fits/")

    parser.add_argument('-n_use',   dest="nights_use",       action="store",
                        help="If you don't want all process all nights under the Input_Data folder, give an array of night you wish to process here. e.g., [20181111, 20181112]",
                        type=str,   default='')
    parser.add_argument('-DeBug',    dest="debug",           action="store_true",
                        help="If sets, will generate files and plots under ./Temp/Debug for debug")
    parser.add_argument('--version',                         action='version',  version='%(prog)s 0.5')
    args = parser.parse_args()
    cdbs_loc = '~/cdbs/'
    inpath     = '../Input_Data/{}/'.format(args.targname.replace(' ', ''))

    if args.debug:
        try:
            os.listdir('../Temp/Debug/{}/'.format(args.targname))
        except OSError:
            os.mkdir('../Temp/Debug/{}/'.format(args.targname))

#-------------------------------------------------------------------------------
    start_time = datetime.now()
    print('\n')
    print('###############################################################')
    print('Fetching Wavelength Regions to be Analyzed for {} (2/3)...'.format(args.targname.replace(' ', '')))

    use_w(args)

    print('Fetching Done!')
    print('###############################################################')
    print('---------------------------------------------------------------')
#-------------------------------------------------------------------------------
    print('A0 Fitting using TelFit for {} (2/2)...'.format(args.targname.replace(' ', '')))
    print('This will take a while..........')
    print('\n')

    bounddata = Table.read('../Input_Data/Use_w/XRegions_{}_{}.csv'.format(args.WRegion, args.band), format='csv')
    starts  = np.array(bounddata['start'])
    ends    = np.array(bounddata['end'])
    labels  = np.array(bounddata['label'], dtype=str)
    xbounddict = {labels[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}

    targname   = args.targname.replace(' ', '')

    ## Collect relevant file information from Predata files
    A0data = Table.read('../Temp/Prepdata/Prepdata_A0_{}.txt'.format(targname), format='ascii')

    ind    = [i != 'NA' for i in A0data['humid']]
    humids = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['humid'])}
    tags   = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['tag'])}
    obs    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['obs'])}
    temps  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['temp'])}
    zds    = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['zd'])}
    press  = {str(k):str(v) for k,v in zip(A0data[ind]['night'],A0data[ind]['press'])}
    nightsFinal = np.array(list(sorted(set(A0data[ind]['night']))))

    print('Analyze with {} nights'.format(len(nightsFinal)))
    time.sleep(6)
    print('\n')
#-------------------------------------------------------------------------------
    if not os.path.isdir('../A0_Fits/'):
        os.mkdir('../A0_Fits/')

    filesndirs = os.listdir('../A0_Fits/')
    name = 'A0_Fits_'+ args.targname
    if name not in filesndirs:
        os.mkdir('../A0_Fits/{}'.format(name) )

    if not os.path.isdir('../A0_Fits/{}/figs_{}'.format(name, args.band)):
        os.mkdir('../A0_Fits/{}/figs_{}'.format(name, args.band) )

    outpath = '../A0_Fits/' + name

    # Retrieve stellar and telluric templates
    watm, satm, mwave0, mflux0 = setup_templates()

    inparam = inparamsA0(inpath,outpath,args.plotfigs,tags,nightsFinal,humids,
                         temps,zds,press,obs,watm,satm,mwave0,mflux0,cdbs_loc,xbounddict)

#    global order

    orders = [ int(labels[i].split('-')[0]) for i in range(len(labels)) ]
    orders = np.unique(orders)
    orders = np.sort(orders)
    for jerp in range(len(orders)):
#    for jerp in range(1):
        outs = mp_run(args, args.Nthreads, jerp, orders, nightsFinal)

    print('\n')
    print('A0 Fitting Done!')

    end_time = datetime.now()
    print('A0 Fitting using TelFit finished, Duration: {}'.format(end_time - start_time))
    print('You can start to run main_step2.py for RV initial guess')
    print('###############################################################')
    print('\n')
