from Engine.importmodule import *
from Engine.opt   import fmod
# from Engine.opt_rebintel   import fmod
from Engine.rebin_jv import rebin_jv

def outplotter_tel(parfit, fitobj, title, inparam, args, order):
    '''
    Plots model fit to telluric standard observation.

    Inputs:
    parfit     : Best fit parameters
    fitobj     : Class containing spectral data to be fit and templates for use in fit
    title      : Title of plot file
    inparam    : Class containing variety of information (e.g. on observing conditions)
    order      : Echelle order, as characterized by file index (as opposed to m number; for conversion between the two, see Stahl et al. 2021)
    args       : Information as input by user from command line
    '''

    fit,chi = fmod(parfit, fitobj)

    #--- to match the wavelengh scale of fit from fmod ---
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    xdata = fitobj.x.copy(); sdata = fitobj.s.copy(); 
    #---
    npars = len(parfit)

    mask = np.ones_like(sdata,dtype=bool)
    mask[(sdata < .0)] = False

    if len(fitobj.mask) != 0:
        for maskbounds in fitobj.mask:
            mask[(xdata > maskbounds[0]) & (xdata < maskbounds[1]) ] = False

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(xdata >= fitobj.CRmask[0][mb]-1) & (xdata <= fitobj.CRmask[0][mb]+1)] = False


    if args.band == 'H':
        if np.int(order) in [13]:
            npars -= 4
        elif np.int(order) in [6,14,21]:
            npars -= 3
        else:
            pass
    else:
        # print("We haven't determined what polynomial orders for K band yet and hardcoded this!")
        # if np.int(order) in [3]:
        #     npars -= 4
        if np.int(order) in [3,4,5]:
            npars -= 3
        else:
            pass

    if fitobj.masterbeam == 'B':
        npars -= 5

    npars -= 6 # subtract 6 from npars total: 2 for linear/quadratic IP, 1 for RV_telluric, 2 fot stellar template power and RV, 1 for vsini

    chi_new = chi*(len(sdata[mask]) - len(parfit))/(len(sdata[mask]) - npars) # correct reduce chisq

    w = parfit[6] + parfit[7]*xdata + parfit[8]*(xdata**2.) + parfit[9]*(xdata**3.)

    cont = parfit[10] + parfit[11]*xdata+ parfit[12]*(xdata**2) + parfit[20]*(xdata**3) + parfit[21]*(xdata**4) + parfit[22]*(xdata**5) + parfit[23]*(xdata**6)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(xdata >= (parfit[15]-parfit[16]/2)) & (xdata <= (parfit[15]+parfit[16]/2))] = parfit[17]
        bucket[(xdata >= (parfit[15]+parfit[16]/2-parfit[18])) & (xdata <= (parfit[15]+parfit[16]/2))] += parfit[19]
        cont -= bucket

    c2 = fitobj.continuum
    cont *= c2

    mask2 = np.ones_like(xdata,dtype=bool)

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask2[(xdata >= fitobj.CRmask[0][mb]-1) & (xdata <= fitobj.CRmask[0][mb]+1)] = False


    fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=250)

    axes.plot(w        ,sdata, '-',  c = 'k',        lw=0.7, label='data',  alpha=.3)
    axes.plot(w[mask2] ,sdata[mask2], '-',  c='k',       lw=0.7, label='data (emission removed)',  alpha=.8)
    axes.plot(w[mask2] ,fit[mask2],      '--', c='tab:red', lw=0.7, label='model', alpha=.8)
    axes.plot(w[mask2] ,cont[mask2],     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.8)
    axes.set_title( title,                 size=6, style='normal', family='sans-serif')
    axes.set_ylabel(r'Flux',    size=6, style='normal', family='sans-serif')
    axes.set_xlabel(r'Wavelength [$\AA$]', size=6, style='normal', family='sans-serif')
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(2))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
    axes.legend(fontsize=5, edgecolor='white')
    fig.text(0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(chi_new),
                        size=6, style='normal', family='sans-serif')

    fig.savefig('{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title),
                bbox_inches='tight', format='png', overwrite=True)


def outplotter_23(parfit, fitobj, title, trk, inparam, args, step2or3, order):
    '''
    Plots model fit to science target observation.

    Inputs:
    parfit     : Best fit parameters
    fitobj     : Class containing spectral data to be fit and templates for use in fit
    title      : Title of plot file
    trk        : Number of run (e.g. RV_results_1, RV_results_2)
    inparam    : Class containing variety of information (e.g. on observing conditions)
    args       : Information as input by user from command line
    step2or3   : Whether run is Step 2 or Step 3
    order      : Echelle order, as characterized by file index (as opposed to m number; for conversion between the two, see Stahl et al. 2021)
    '''

    fit,chi = fmod(parfit, fitobj)

    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    xdata = fitobj.x.copy(); sdata = fitobj.s.copy(); 

    npars = len(parfit)

    mask = np.ones_like(sdata,dtype=bool)
    mask[(sdata < .0)] = False

    if len(fitobj.mask) != 0:
        for maskbounds in fitobj.mask:
            mask[(xdata > maskbounds[0]) & (xdata < maskbounds[1]) ] = False

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(xdata >= fitobj.CRmask[0][mb]-1) & (xdata <= fitobj.CRmask[0][mb]+1)] = False

    if args.band == 'H':
        if np.int(order) in [13]:
            npars -= 4
        elif np.int(order) in [6,14,21]:
            npars -= 3
        else:
            pass
    else:
        # print("We haven't determined what polynomial orders for K band yet and hardcoded this!")
        # if np.int(order) in [3]:
        #     npars -= 4
        if np.int(order) in [3,4,5]:
            npars -= 3
        else:
            pass

    if fitobj.masterbeam == 'B':
        npars -= 5

    npars -= 3 # subtract 3 from npars total, 2 for linear/quadratic IP and 1 for RV_telluric

    chi_new = chi*(len(sdata[mask]) - len(parfit))/(len(sdata[mask]) - npars)

    # Apply continuum adjustment
    cont = parfit[10] + parfit[11]*xdata+ parfit[12]*(xdata**2) + parfit[20]*(xdata**3) + parfit[21]*(xdata**4) + parfit[22]*(xdata**5) + parfit[23]*(xdata**6)


    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(xdata >= (parfit[15]-parfit[16]/2))         & (xdata <= (parfit[15]+parfit[16]/2))] = parfit[17]
        bucket[(xdata >= (parfit[15]+parfit[16]/2-parfit[18])) & (xdata <= (parfit[15]+parfit[16]/2))] += parfit[19]
        cont -= bucket

    c2 = fitobj.continuum
    cont *= c2

    fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=250)

    mask2 = np.ones_like(xdata,dtype=bool)

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask2[(xdata >= fitobj.CRmask[0][mb]-1) & (xdata <= fitobj.CRmask[0][mb]+1)] = False

    n = len(fitobj.mask)

    if n > 0:
        widths = [fitobj.mask[0][0]-xdata[0]]
        for m in range(n-1):
            widths.append(fitobj.mask[m+1][0]-fitobj.mask[m][1])
        widths.append(xdata[-1]-fitobj.mask[n-1][1])
        gs = gridspec.GridSpec(1, n+1, width_ratios=widths)
        for m in range(n+1):
            ax0 = plt.subplot(gs[m])

            ax0.plot(w,       sdata,           '--', c='k',       lw=0.7, label='data',  alpha=.3)
            ax0.plot(w[mask2],sdata[mask2],    '-',  c='k',       lw=0.7, label='data (emission removed)',  alpha=.8)
            ax0.plot(w[mask2],fit[mask2],      '--', c='tab:red', lw=0.7, label='model', alpha=.8)
            ax0.plot(w[mask2],cont[mask2],     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.8)

            if title[6]=='_':
                # only plot residual on the "parfit"
                ax0.plot(w[mask2],fit[mask2] - sdata[mask2], 's', c='k', ms=0.3, mew=0.3, label='residual', alpha=1)
                ax0.axhline(0, color='tab:grey', lw=0.2, zorder=0, alpha=1)

            kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False,lw= 0.6)
            if m == 0:
                ax0.tick_params(axis='both', labelsize=6, right=False, top=True, direction='in')
                left = w[0]
                right = parfit[6] + parfit[7]*fitobj.mask[m][0] + parfit[8]*(fitobj.mask[m][0]**2.) + parfit[9]*(fitobj.mask[m][0]**3.)
                ax0.plot([right,right],[min(sdata),max(sdata)],'--k',lw=0.75)
            elif m == n:
                ax0.tick_params(axis='both', labelsize=6, left=False, right=True, top=True, direction='in')
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                left = parfit[6] + parfit[7]*fitobj.mask[m-1][1] + parfit[8]*(fitobj.mask[m-1][1]**2.) + parfit[9]*(fitobj.mask[m-1][1]**3.)
                right = w[-1]
            else:
                ax0.tick_params(axis='both', labelsize=6, right=False, left=False,top=True, direction='in')
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                ax0.plot([right,right],[min(sdata),max(sdata)],'--k',lw=0.75)
                left = parfit[6] + parfit[7]*fitobj.mask[m-1][1] + parfit[8]*(fitobj.mask[m-1][1]**2.) + parfit[9]*(fitobj.mask[m-1][1]**3.)
                right = parfit[6] + parfit[7]*fitobj.mask[m][0] + parfit[8]*(fitobj.mask[m][0]**2.) + parfit[9]*(fitobj.mask[m][0]**3.)

            ax0.set_xlim(left,right)

            if m != 0:
                ax0.spines['left'].set_visible(False)
            if m != n:
                ax0.spines['right'].set_visible(False)

        fig.tight_layout(pad=0.0)
        fig.suptitle( title,     x=0.5,y=1.05, size=6, style='normal', family='sans-serif')
        fig.text(0.5, -0.04, r'Wavelength [$\rm\AA$]', ha='center', size=6, style='normal', family='sans-serif')
        fig.text(-0.04, 0.5, r'Flux',       va='center', rotation='vertical', size=6, style='normal', family='sans-serif')
        fig.text(0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(chi_new),
                            size=6, style='normal', family='sans-serif')
        ax0.legend(fontsize=5, edgecolor='white', markerscale=2.5)

    else:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)
        axes.plot(w,       sdata,        '--', c='k',       lw=0.7, label='data',  alpha=.3)
        axes.plot(w[mask2],sdata[mask2], '-',  c='k',       lw=0.7, label='data (emission removed)',  alpha=.8)
        axes.plot(w[mask2],fit[mask2],      '--', c='tab:red', lw=0.7, label='model', alpha=.8)
        axes.plot(w[mask2],cont[mask2],     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.8)

        if title[6]=='_':
            # only plot residual on the "parfit"
            axes.plot(w[mask2],fit[mask2] - sdata[mask2], 's', c='k', ms=0.3, mew=0.3, label='residual', alpha=0.8)
            axes.axhline(0, color='tab:grey', lw=0.2, zorder=0, alpha=1)

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.set_title(title,  size=6, style='normal' , family='sans-serif' )
        axes.set_ylabel(r'Flux',        size=6, style='normal', family='sans-serif' )
        axes.set_xlabel(r'Wavelength [$\rm\AA$]',  size=6, style='normal', family='sans-serif' )
        fig.text(0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(chi_new),
                            size=6, style='normal', family='sans-serif')

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.legend(fontsize=5, edgecolor='white', markerscale=2.5)

    fig.savefig(f'{inparam.outpath}/figs/main_step{step2or3}_{args.band}_{trk}/{title}.png', bbox_inches='tight', format='png', overwrite=True)
