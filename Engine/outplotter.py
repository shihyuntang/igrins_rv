from Engine.importmodule import *
from Engine.opt       import fmod

def outplotter_tel(parfit, fitobj, title, inparam, args):
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
    fig.savefig('{}/figs_{}/{}.png'.format(inparam.outpath, args.band, title), bbox_inches='tight', format='png', overwrite=True)


def outplotter_23(parfit,fitobj,title,trk):
    fit,chi = fmod(parfit, fitobj)
    w = parfit[6] + parfit[7]*fitobj.x + parfit[8]*(fitobj.x**2.) + parfit[9]*(fitobj.x**3.)

    fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)


    n = len(fitobj.mask)

    if n > 0:

        widths = [fitobj.mask[0][0]-fitobj.x[0]]
        for m in range(n-1):
            widths.append(fitobj.mask[m+1][0]-fitobj.mask[m][1])
        widths.append(fitobj.x[-1]-fitobj.mask[n-1][1])
        gs = gridspec.GridSpec(1, n+1, width_ratios=widths)
        for m in range(n+1):
            ax0 = plt.subplot(gs[m])
            ax0.plot(w,fitobj.s, '-k',  lw=0.5, label='data',  alpha=.6)
            ax0.plot(w,fit,      '--r', lw=0.5, label='model',  alpha=.6)
            kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False,lw= 0.6)
            if m == 0:
                ax0.tick_params(axis='both', labelsize=4.5, right=False, top=True, direction='in')
                left = w[0]
                right = parfit[6] + parfit[7]*fitobj.mask[m][0] + parfit[8]*(fitobj.mask[m][0]**2.) + parfit[9]*(fitobj.mask[m][0]**3.)
                ax0.plot([right,right],[min(fitobj.s),max(fitobj.s)],'--k',lw=0.75)
            elif m == n:
                ax0.tick_params(axis='both', labelsize=4.5, left=False, right=True, top=True, direction='in')
                ax0.set_yticklabels([])
                #ax0.vlines([left+1],'--k',lw=0.75)
                left = parfit[6] + parfit[7]*fitobj.mask[m-1][1] + parfit[8]*(fitobj.mask[m-1][1]**2.) + parfit[9]*(fitobj.mask[m-1][1]**3.)
                right = w[-1]
            else:
                ax0.tick_params(axis='both', labelsize=4.5, right=False, left=False,top=True, direction='in')
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
        fig.suptitle( title,     x=0.5,y=1.05,size=5, style='normal', family='sans-serif')
        fig.text(0.5, -0.04, r'Wavelength [$\rm\AA$]', ha='center',size=5, style='normal', family='sans-serif')
        fig.text(-0.04, 0.5, r'Normalized Flux', va='center', rotation='vertical',size=5, style='normal', family='sans-serif')
        ax0.legend(fontsize=4, edgecolor='white')

    else:
        fig, axes = plt.subplots(1, 1, figsize=(5,3), facecolor='white', dpi=300)
        axes.plot(w,fitobj.s, '-k',  lw=0.5, label='data',alpha=.6)
        axes.plot(w,fit,      '--r', lw=0.5, label='model',alpha=.6)

        axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
        axes.set_title(title,   size=5, style='normal' , family='sans-serif' )
        axes.set_ylabel(r'Normalized Flux',   size=5, style='normal' , family='sans-serif' )
        axes.set_xlabel(r'Wavelength [$\AA$]',       size=5, style='normal' , family='sans-serif' )

        axes.tick_params(axis='both', labelsize=4.5, right=True, top=True, direction='in')
        axes.legend(fontsize=4, edgecolor='white')

    fig.savefig(f'{inparam.outpath}/main_step2_figs_{args.band}_{trk}/{title}.png', bbox_inches='tight', format='png', overwrite=True)
