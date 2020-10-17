from Engine.importmodule import *
from Engine.optwithtwodip   import fmod

def outplotter_tel(parfit, fitobj, title, inparam, args):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    c2 = fitobj.continuum
    c2 = c2#/np.median(c2)
    cont = parfit[10] + parfit[11]*fitobj.x+ parfit[12]*(fitobj.x**2)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (parfit[15]-parfit[16]/2)) & (fitobj.x <= (parfit[15]+parfit[16]/2))] = parfit[17]
        bucket[(fitobj.x >= (parfit[15]+parfit[16]/2-parfit[18])) & (fitobj.x <= (parfit[15]+parfit[16]/2))] += parfit[19]
        cont -= bucket
    cont *= c2

    fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)

    axes.plot(w,fitobj.s, '-',  c = 'k',        lw=0.7, label='data',  alpha=.6)
    axes.plot(w,fit,      '--', c = 'tab:red',  lw=0.7, label='model', alpha=.6)
    axes.plot(w,cont,     '--', c = 'tab:blue',  lw=0.7, label='cont', alpha=.6)
    axes.set_title( title,                 size=6, style='normal', family='sans-serif')
    axes.set_ylabel(r'Flux',    size=6, style='normal', family='sans-serif')
    axes.set_xlabel(r'Wavelength [$\AA$]', size=6, style='normal', family='sans-serif')
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(2))
    axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
    axes.legend(fontsize=5, edgecolor='white')
    fig.savefig('{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title),
                bbox_inches='tight', format='png', overwrite=True)

    # fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)
    #
    # axes.plot(fitobj.x,fitobj.s, '-',  c = 'k',        lw=0.7, label='data',  alpha=.6)
    # axes.plot(fitobj.x,fit,      '--', c = 'tab:red',  lw=0.7, label='model', alpha=.6)
    # axes.plot(fitobj.x,cont,     '--', c = 'tab:blue',  lw=0.7, label='cont', alpha=.6)
    # axes.set_title( title,                 size=6, style='normal', family='sans-serif')
    # axes.set_ylabel(r'Flux',    size=6, style='normal', family='sans-serif')
    # axes.set_xlabel(r'Wavelength [$\AA$]', size=6, style='normal', family='sans-serif')
    # axes.xaxis.set_minor_locator(AutoMinorLocator(5))
    # axes.yaxis.set_minor_locator(AutoMinorLocator(2))
    # axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
    # axes.legend(fontsize=5, edgecolor='white')
    # fig.savefig('{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title+str('_X')),
    #             bbox_inches='tight', format='png', overwrite=True)


def outplotter_23(parfit, fitobj, title, trk, inparam, args, step2or3):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    c2 = fitobj.continuum
    c2 = c2#/np.median(c2)
    cont = parfit[10] + parfit[11]*fitobj.x+ parfit[12]*(fitobj.x**2)
    if fitobj.masterbeam == 'A':
        bucket = np.zeros_like(cont)
        bucket[(fitobj.x >= (parfit[15]-parfit[16]/2)) & (fitobj.x <= (parfit[15]+parfit[16]/2))] = parfit[17]
        bucket[(fitobj.x >= (parfit[15]+parfit[16]/2-parfit[18])) & (fitobj.x <= (parfit[15]+parfit[16]/2))] += parfit[19]
        cont -= bucket
    cont *= c2

    fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)

    n = len(fitobj.mask)

    if n > 0:
        widths = [fitobj.mask[0][0]-fitobj.x[0]]
        for m in range(n-1):
            widths.append(fitobj.mask[m+1][0]-fitobj.mask[m][1])
        widths.append(fitobj.x[-1]-fitobj.mask[n-1][1])
        gs = gridspec.GridSpec(1, n+1, width_ratios=widths)
        for m in range(n+1):
            ax0 = plt.subplot(gs[m])
            ax0.plot(w,fitobj.s, '-',  c='k',       lw=0.7, label='data',  alpha=.6)
            ax0.plot(w,fit,      '--', c='tab:red', lw=0.7, label='model', alpha=.6)
            ax0.plot(w,cont,     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.6)
            kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False,lw= 0.6)
            if m == 0:
                ax0.tick_params(axis='both', labelsize=6, right=False, top=True, direction='in')
                left = w[0]
                right = parfit[6] + parfit[7]*fitobj.mask[m][0] + parfit[8]*(fitobj.mask[m][0]**2.) + parfit[9]*(fitobj.mask[m][0]**3.)
                ax0.plot([right,right],[min(fitobj.s),max(fitobj.s)],'--k',lw=0.75)
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
                ax0.plot([right,right],[min(fitobj.s),max(fitobj.s)],'--k',lw=0.75)
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
        ax0.legend(fontsize=5, edgecolor='white')

    else:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)
        axes.plot(w,fitobj.s, '-',  c='k',       lw=0.7, label='data',  alpha=.6)
        axes.plot(w,fit,      '--', c='tab:red', lw=0.7, label='model', alpha=.6)
        axes.plot(w,cont,     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.6)

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.set_title(title,  size=6, style='normal' , family='sans-serif' )
        axes.set_ylabel(r'Flux',        size=6, style='normal', family='sans-serif' )
        axes.set_xlabel(r'Wavelength [$\rm\AA$]',  size=6, style='normal', family='sans-serif' )

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.legend(fontsize=5, edgecolor='white')

    fig.savefig(f'{inparam.outpath}/figs/main_step{step2or3}_{args.band}_{trk}/{title}.png', bbox_inches='tight', format='png', overwrite=True)


    fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)

    n = len(fitobj.mask)

    if n > 0:
        widths = [fitobj.mask[0][0]-fitobj.x[0]]
        for m in range(n-1):
            widths.append(fitobj.mask[m+1][0]-fitobj.mask[m][1])
        widths.append(fitobj.x[-1]-fitobj.mask[n-1][1])
        gs = gridspec.GridSpec(1, n+1, width_ratios=widths)
        for m in range(n+1):
            ax0 = plt.subplot(gs[m])
            ax0.plot(fitobj.x,fitobj.s, '-',  c='k',       lw=0.7, label='data',  alpha=.6)
            ax0.plot(fitobj.x,fit,      '--', c='tab:red', lw=0.7, label='model', alpha=.6)
            ax0.plot(fitobj.x,cont,     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.6)
            kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False,lw= 0.6)
            if m == 0:
                ax0.tick_params(axis='both', labelsize=6, right=False, top=True, direction='in')
                left = fitobj.x[0]
                right = fitobj.mask[m][0]
                ax0.plot([right,right],[min(fitobj.s),max(fitobj.s)],'--k',lw=0.75)
            elif m == n:
                ax0.tick_params(axis='both', labelsize=6, left=False, right=True, top=True, direction='in')
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                left = fitobj.mask[m-1][1]
                right = fitobj.x[-1]
            else:
                ax0.tick_params(axis='both', labelsize=6, right=False, left=False,top=True, direction='in')
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                ax0.plot([right,right],[min(fitobj.s),max(fitobj.s)],'--k',lw=0.75)
                left = fitobj.mask[m-1][1]
                right = fitobj.mask[m][0]

            ax0.set_xlim(left,right)

            if m != 0:
                ax0.spines['left'].set_visible(False)
            if m != n:
                ax0.spines['right'].set_visible(False)

        fig.tight_layout(pad=0.0)
        fig.suptitle( title,     x=0.5,y=1.05, size=6, style='normal', family='sans-serif')
        fig.text(0.5, -0.04, r'Pixel', ha='center', size=6, style='normal', family='sans-serif')
        fig.text(-0.04, 0.5, r'Flux',       va='center', rotation='vertical', size=6, style='normal', family='sans-serif')
        ax0.legend(fontsize=5, edgecolor='white')

    else:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)
        axes.plot(fitobj.x,fitobj.s, '-',  c='k',       lw=0.7, label='data',  alpha=.6)
        axes.plot(fitobj.x,fit,      '--', c='tab:red', lw=0.7, label='model', alpha=.6)
        axes.plot(fitobj.x,cont,     '--', c='tab:blue',  lw=0.7, label='cont', alpha=.6)

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.set_title(title,  size=6, style='normal' , family='sans-serif' )
        axes.set_ylabel(r'Flux',        size=6, style='normal', family='sans-serif' )
        axes.set_xlabel(r'Pixel',  size=6, style='normal', family='sans-serif' )

        axes.tick_params(axis='both', labelsize=6, right=True, top=True, direction='in')
        axes.legend(fontsize=5, edgecolor='white')

    fig.savefig(f'{inparam.outpath}/figs/main_step{step2or3}_{args.band}_{trk}/{title}_X.png', bbox_inches='tight', format='png', overwrite=True)



#     print('''
# ************************************************************************************
# ___    ____    ___    ___    ___   __     _    ____       ___   __      __
#  |    /     \   |    |    \   |    | \    |  /           |    \  \      /
#  |   |  ____    |    |___ /   |    |  \   |  | ___       |____/   \    /
#  |   |      |   |    |    \   |    |   \  |        |     |    \    \  /
# ---   \____/   ---   |    |  ---   |    \_|   ____ /     |    |     \/
# ************************************************************************************
#     ''')
