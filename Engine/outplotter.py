import matplotlib as mpl

from Engine.importmodule import *
from Engine.opt   import fmod
from Engine.rebin_jv import rebin_jv

mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.family'] = 'sans-serif'

def outplotter_tel(parfit, fitobj, title, inparam, args, order, chi_new):
    '''
    Plots model fit to telluric standard observation.

    Inputs:
    parfit     : Best fit parameters
    fitobj     : Class containing spectral data to be fit and templates for use in fit
    title      : Title of plot file
    inparam    : Class containing variety of information (e.g. on observing conditions)
    order      : Echelle order, as characterized by file index (as opposed to
                    m number; for conversion between the two, see Stahl et al. 2021)
    args       : Information as input by user from command line
    '''

    fit,chi,w,cont = fmod(parfit, fitobj)

    #--- to match the wavelengh scale of fit from fmod ---

    xdata = fitobj.x.copy()
    sdata = fitobj.s.copy()
    #---
    npars = len(parfit)

    mask = np.ones_like(sdata, dtype=bool)
    mask[(sdata < .0)] = False

    if len(fitobj.mask) != 0:
        for maskbounds in fitobj.mask:
            mask[(xdata > maskbounds[0]) \
                    & (xdata < maskbounds[1]) ] = False

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(xdata >= fitobj.CRmask[0][mb]-1) \
                    & (xdata <= fitobj.CRmask[0][mb]+1)] = False

    mask = np.ones_like(xdata, dtype=bool)

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(xdata >= fitobj.CRmask[0][mb]-1) \
                    & (xdata <= fitobj.CRmask[0][mb]+1)] = False

    if len(fitobj.molmask) > 0:
        for mb in fitobj.molmask:
            mask[(xdata >= mb[0]) & (xdata <= mb[1])] = False

    fig, axes = plt.subplots(1, 1, figsize=(6,3), dpi=200)

    axes.plot(w ,sdata, '-', c = 'k', lw=0.7, label='data', alpha=.3)
    axes.plot(
        w[mask] ,sdata[mask], '-',  c='k', lw=0.7,
        label='data (emission removed)', alpha=.8
        )
    axes.plot(
        w[mask] ,fit[mask], '--', c='tab:red', lw=0.7, label='model', 
        alpha=.8
        )
    axes.plot(
        w[mask] ,cont[mask], '--', c='tab:blue',  lw=0.7, label='cont', 
        alpha=.8
        )

    axes.set_title(title, size=6)
    axes.set_ylabel(r'Flux', size=6)
    axes.set_xlabel(r'Wavelength [$\AA$]', size=6)
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(2))
    axes.tick_params(
        axis='both', which='both', labelsize=6, right=True, top=True, 
        direction='in'
        )
    axes.legend(fontsize=5, edgecolor='white')
    fig.text(
        0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(chi_new),
        size=6
        )

    fig.savefig(
        '{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title),
        bbox_inches='tight', format='png'
        )


def outplotter_23(parfit, fitobj, title, trk, inparam, args, step2or3, order, chi_new):
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
    order      : Echelle order, as characterized by file index (as opposed to m
                    number; for conversion between the two, see Stahl et al. 2021)
    '''

    fit,chi,w,cont = fmod(parfit, fitobj,args.binary)

    xdata = fitobj.x.copy()
    sdata = fitobj.s.copy()


    mask = np.ones_like(sdata, dtype=bool)
    mask[(sdata < .0)] = False

    if len(fitobj.mask) != 0:
        for maskbounds in fitobj.mask:
            mask[(xdata > maskbounds[0]) & (xdata < maskbounds[1]) ] = False

    if len(fitobj.CRmask[1]) > 0:
        for mb in fitobj.CRmask[1]:
            mask[(xdata >= fitobj.CRmask[0][mb]-1) \
                    & (xdata <= fitobj.CRmask[0][mb]+1)] = False

    if len(fitobj.molmask) > 0:
        for mb in fitobj.molmask:
            mask[(xdata >= mb[0]) & (xdata <= mb[1])] = False


    fig, axes = plt.subplots(1, 1, figsize=(6,3), dpi=250)

    n = len(fitobj.mask)

    # masked_w = w.copy()
    # masked_w[~mask] = np.nan

    masked_sdata = sdata.copy()
    masked_sdata[~mask] = np.nan

    masked_fit = fit.copy()
    masked_fit[~mask] = np.nan

    masked_cont = cont.copy()
    masked_cont[~mask] = np.nan

    if n > 0:
        widths = [fitobj.mask[0][0]-xdata[0]]
        for m in range(n-1):
            widths.append(fitobj.mask[m+1][0]-fitobj.mask[m][1])
        widths.append(xdata[-1]-fitobj.mask[n-1][1])
        gs = gridspec.GridSpec(1, n+1, width_ratios=widths)
        for m in range(n+1):
            ax0 = plt.subplot(gs[m])

            ax0.plot(w, sdata, '--', c='k', lw=0.7, label='data', alpha=.3)
            ax0.plot(
                w, masked_sdata, '-',  c='k', lw=0.7, 
                label='data (emission removed)', alpha=.8
                )
            ax0.plot(
                w, masked_fit, '--', c='tab:red', lw=0.7, label='model', 
                alpha=.8
                )
            ax0.plot(
                w, masked_cont, '--', c='tab:blue',  lw=0.7, label='cont', 
                alpha=.8
                )

            if title[6]=='_':
                # only plot residual on the "parfit"
                ax0.plot(
                    w, masked_fit - masked_sdata, 's', c='k', ms=0.3,
                    mew=0.3, label='residual', alpha=1
                    )
                ax0.axhline(
                    0, color='tab:grey', lw=0.2, zorder=0, alpha=1
                    )

            kwargs = dict(
                transform=ax0.transAxes, color='k', clip_on=False,
                lw= 0.6
                )
            if m == 0:
                ax0.tick_params(
                    axis='both', labelsize=6, right=False, top=True,
                    direction='in'
                    )
                left = w[0]
                right = w[(xdata >= fitobj.mask[m][0])][0]
                ax0.plot(
                    [right,right], [np.min(sdata), np.max(sdata)],
                    '--k',lw=0.75
                    )
            elif m == n:
                ax0.tick_params(
                    axis='both', labelsize=6, left=False, right=True,
                    top=True, direction='in'
                    )
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                left = w[(xdata >= fitobj.mask[m-1][1])][0]
                right = w[-1]
            else:
                ax0.tick_params(
                    axis='both', labelsize=6, right=False, left=False,
                    top=True, direction='in'
                    )
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                ax0.plot(
                    [right,right], [np.min(sdata), np.max(sdata)],
                    '--k',  lw=0.75
                    )
                left = w[(xdata >= fitobj.mask[m-1][1])][0]
                right = w[(xdata >= fitobj.mask[m][0])][0]

            ax0.set_xlim(left, right)

            if m != 0:
                ax0.spines['left'].set_visible(False)
            if m != n:
                ax0.spines['right'].set_visible(False)

        fig.tight_layout(pad=0.0)
        fig.suptitle( 
            title, x=0.5, y=1.05, size=6)
        fig.text(
            0.5, -0.04, r'Wavelength [$\rm\AA$]', ha='center', size=6
            )
        fig.text(
            -0.04, 0.5, 'Flux', va='center', rotation='vertical', size=6
            )
        fig.text(
            0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(
                chi_new), size=6
                )
        ax0.legend(fontsize=5, edgecolor='white', markerscale=2.5)

    else:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), dpi=200)
        axes.plot(w, sdata, '--', c='k', lw=0.7, label='data', alpha=.3)
        axes.plot(
            w, masked_sdata, '-',  c='k', lw=0.7, 
            label='data (emission removed)', alpha=.8
            )
        axes.plot(
            w, masked_fit, '--', c='tab:red', lw=0.7, label='model', 
            alpha=.8
            )
        axes.plot(
            w, masked_cont, '--', c='tab:blue', lw=0.7, label='cont', 
            alpha=.8
            )

        if title[6]=='_':
            # only plot residual on the "parfit"
            axes.plot(
                w, masked_fit - masked_sdata, 's', c='k', ms=0.3, mew=0.3, 
                label='residual', alpha=0.8
                )
            axes.axhline(0, color='tab:grey', lw=0.2, zorder=0, alpha=1)

        axes.tick_params(
            axis='both', labelsize=6, right=True, top=True, direction='in'
            )
        axes.set_title(title, size=6)
        axes.set_ylabel('Flux', size=6)
        axes.set_xlabel(r'Wavelength [$\rm\AA$]',  size=6)
        fig.text(
            0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(
                chi_new), size=6
                )

        axes.tick_params(
            axis='both', labelsize=6, right=True, top=True, direction='in'
            )
        axes.legend(fontsize=5, edgecolor='white', markerscale=2.5)

    fig.savefig(
        f'{inparam.outpath}/figs/main_step{step2or3}_{args.band}_{trk}/{title}.png',
        bbox_inches='tight', format='png', overwrite=True)

def outplotter_rv(rvp, stdp, args, inparam, name, ara, kind):

    f, axes = plt.subplots(1, 1, figsize=(5,3), dpi=300)

    # axes.plot(np.arange(len(rvp))+1, rvp, '.k', ms=5)
    axes.errorbar(
        np.arange(len(rvp))+1, rvp, yerr=stdp, fmt='.', ls='none', lw=.5, 
        ms=5, color='k', ecolor='k'
        )
    
    axes.text(
        0.05, 0.93, r'RV mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(
        np.nanmean(rvp), np.nanstd(rvp)), transform=axes.transAxes, 
        size=6 
        )
    axes.set_ylim(np.nanmin(rvp)-.08, np.nanmax(rvp)+.08)
    axes.set_ylabel('RV [km/s]', size=6)
    axes.set_xlabel('Night (#)', size=6)
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(
        axis='both', which='both', labelsize=5, right=True, top=True, 
        direction='in', width=.6
        )
    if args.binary:
        f.savefig(
            '{}/{}/FinalRVs{}_{}.png'.format(
                inparam.outpath, name, ara+1, kind
                ), format='png', bbox_inches='tight')
    else:
        f.savefig(
            '{}/{}/FinalRVs_{}.png'.format(
                inparam.outpath, name, kind
                ), format='png', bbox_inches='tight')

def outplotter_rv_combind(
    xscale, rvp, stdp, args, inparam, name, ara, kind, nightsT, nightsL):

    f, axes = plt.subplots(1, 1, figsize=(5,3), dpi=300)
    
    axes.errorbar(
        xscale, rvp, yerr=stdp, fmt='.', ls='none', lw=.5, 
        ms=5, color='k', ecolor='k'
        )

    axes.text(
        0.05, 0.93, r'RV mean= ${:1.5f}$ $\pm$ {:1.5f} km/s'.format(
        np.nanmean(rvp), np.nanstd(rvp)),
        transform=axes.transAxes, size=6 )

    if (len(nightsT) != 0) & (len(nightsL) == 0):
        axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6)
    elif (len(nightsT) == 0) & (len(nightsL) != 0):
        axes.text(0.05, 0.1, 'Defocus', transform=axes.transAxes, size=6)
    else:
        if nightsT[-1] < nightsL[0]: # if tight epoch precedes loose epoch #sy
            axes.axvline(xscale[len(nightsT)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6)
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6)
        else:
            axes.axvline(xscale[len(nightsL)] - 0.5, linewidth=.7, color='black')
            axes.text(0.05, 0.1, 'Focused', transform=axes.transAxes, size=6)
            axes.text(0.9,  0.1, 'Defocus', transform=axes.transAxes, size=6)
    
    axes.set_ylim(np.nanmin(rvp)-.08,np.nanmax(rvp)+.08)
    axes.set_ylabel('RV (km/s)', size=6 )
    axes.set_xlabel('Night (#)', size=6 )
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(5))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True,
        direction='in', width=.6)
    if args.binary:
        f.savefig(
            '{}/{}/FinalRVs{}.png'.format(
                inparam.outpath, name, ara+1), 
                format='png', bbox_inches='tight'
                )
    else:
        f.savefig(
            '{}/{}/FinalRVs.png'.format(
                inparam.outpath, name), 
                format='png', bbox_inches='tight'
                )
