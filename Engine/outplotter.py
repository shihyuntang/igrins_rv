from Engine.importmodule import *
from Engine.opt       import fmod

def outplotter_tel(parfit, fitobj, title):
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
