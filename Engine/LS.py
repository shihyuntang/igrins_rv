from Engine.importmodule import *
from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.rebin_jv   import rebin_jv
from Engine.detect_peaks import detect_peaks
from astropy.timeseries import LombScargle
from PyAstronomy.pyasl import foldAt


def LSandfold(x,y,u,args,out,title):

    Protdict = {'V830Tau':2.741,'V827Tau':3.76,'Hubble4':1.55,'V1075Tau':2.43,
                'DKTau':8.24,'AATau':8.2,'DHTau':7.0,'IQTau':6.25, 'LkCa15':5.7,
                 'GITau':7.03}
    Porbdict = {"CITau":8.989}

    x = x[np.isfinite(y)];     u = u[np.isfinite(y)];    y = y[np.isfinite(y)];

    freq0,pow0 = LombScargle(x,y,u).autopower(normalization='standard',
        maximum_frequency = 1/args.minper,
        minimum_frequency = 1/args.maxper)

    per0 = 1/freq0
    periodfound = per0[np.argmax(pow0)]

    faps = LombScargle(x,y,u).false_alarm_level([0.1,0.05,0.01])

    nmpd = 20
    peaks,extra = find_peaks(pow0,distance=50,height=0.4)

    bestpeaks = [j for _,j in sorted(zip(pow0[peaks],per0[peaks]),reverse=True)] #x is sorted by y

    for alias in [2., 3., 4., 5.]:
        for peak in bestpeaks:
            if abs(bestpeaks[0]*alias - peak) < .05:
                bestpeaks.remove(bestpeaks[0])
                break
            elif abs(bestpeaks[0]/alias - peak) < .05:
                bestpeaks.remove(peak)
                break
            else:
                pass
    for peak in bestpeaks[1:]:
        if abs(bestpeaks[0]- peak) < .05:
            bestpeaks.remove(peak)

    if len(bestpeaks) == 0:
        fig,axs = plt.subplots(len(bestpeaks)+1,1,figsize=(12,8))

        plt.subplots_adjust(hspace=.25)
        axs.plot(per0,pow0,color='black')
        if Prot0 is not None:
            axs.axvline(Prot0,ls='-',color='red',alpha=.6)
            for alias in [1/5.,1/4.,1/3.,1/2.]:
                if Prot0*alias > 1:
                    ypred = pow0[(abs(per0-Prot0*alias) == np.min(abs(per0-Prot0*alias)))]
                    axs.vlines([Prot0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='red',alpha=.4)
        if Porb0 is not None:
            axs.axvline(Porb0,ls='-',color='blue',alpha=.6)
            for alias in [1/5.,1/4.,1/3.,1/2.]:
                if Porb0*alias > 1:
                    ypred = pow0[(abs(per0-Porb0*alias) == np.min(abs(per0-Porb0*alias)))]
                    axs.vlines([Porb0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='blue',alpha=.4)
        ys = axs.get_ylim()
        for fap in faps:
            if fap < ys[-1]:
                axs.axhline(fap,color='black',alpha=.25,ls='--')
        axs.tick_params(axis='x',labelsize=15)
        axs.tick_params(axis='y',labelsize=15)
        axs.set_xlabel('Period (d)',labelpad=15,fontsize=25,weight='normal')
        axs.set_ylabel('Power',labelpad=30,fontsize=25,weight='normal')
        plt.savefig(out+title+'.png',format='png',bbox_inches='tight')
        plt.clf()
        plt.close()
        return


    markers = ['x','d','s','*','^']
    if len(bestpeaks) > len(markers):
        bestpeaks = bestpeaks[:len(markers)]
        print('WARNING! TOO MANY PEAKS IDENTIFIED, SHOWING ONLY TOP FIVE!')


    for kind in ['','_timecolor']:

        fig,axs = plt.subplots(len(bestpeaks)+1,1,figsize=(12,8+8*len(bestpeaks)))

        plt.subplots_adjust(hspace=.25)
        axs[0].plot(per0,pow0,color='black')
        try:
            Prot0 = Protdict[args.targname]
            axs[0].axvline(Prot0,ls='-',color='red',alpha=.6)
        except:
            pass
        try:
            Porb0 = Porbdict[args.targname]
            axs[0].axvline(Porb0,ls='-',color='blue',alpha=.6)
        except:
            pass

        for alias in [1/5.,1/4.,1/3.,1/2.]:
            if Prot0 is not None:
                if Prot0*alias > 1:
                    ypred = pow0[(abs(per0-Prot0*alias) == np.min(abs(per0-Prot0*alias)))]
                    axs[0].vlines([Prot0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='red',alpha=.4)
            if Porb0 is not None:
                if Porb0*alias > 1:
                    ypred = pow0[(abs(per0-Porb0*alias) == np.min(abs(per0-Porb0*alias)))]
                    axs[0].vlines([Porb0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='blue',alpha=.4)

        ypreds = sorted(pow0[peaks],reverse=True)[:len(bestpeaks)]
        markers = ['x','d','s','*','^']
        for j in range(len(bestpeaks)):

            axs[0].scatter(bestpeaks[j],ypreds[j]*1.1,color='black',marker= markers[j] ,label='P = '+str(round(bestpeaks[j],3)))
            for alias in [1/5.,1/4.,1/3.,1/2.]:
                if bestpeaks[j]*alias > 1:
                    ypred = pow0[(abs(per0-bestpeaks[j]*alias) == np.min(abs(per0-bestpeaks[j]*alias)))]
                    axs[0].scatter(bestpeaks[j]*alias,ypred*1.1,color='black',marker=markers[j],s=5,alpha=.4)

        axs[0].set_xlabel('Period (d)',labelpad=15,fontsize=25,weight='normal')
        axs[0].set_ylabel('Power',labelpad=30,fontsize=25,weight='normal')
        axs[0].legend()
        axs[0].tick_params(axis='x',labelsize=15)
        axs[0].tick_params(axis='y',labelsize=15)

        ys = axs[0].get_ylim()
        for fap in faps:
            if fap < ys[-1]:
                axs[0].axhline(fap,color='black',alpha=.25,ls='--')

        if len(bestpeaks) > 0:

            for jj in range(len(bestpeaks)):
                phases     = foldAt(x,  bestpeaks[jj], centralzero=True)
                sortIndi1  = np.argsort(phases);
                ff         = y[sortIndi1]
                uf         = u[sortIndi1]
                jf         = x[sortIndi1]
                xf         = phases[sortIndi1]
                diff = np.nanpercentile(ff,95)-np.nanpercentile(ff,5)
                if kind == '':
                    axs[jj+1].errorbar(xf,ff,yerr=uf,color='black',ls='none',alpha=1)
                    axs[jj+1].scatter(xf,ff,color='black',s=7,alpha=1)
                else:
                    axs[jj+1].scatter(xf,ff,c=jf,cmap=cm.cool,s=7,alpha=1)
                    clb = axs[jj+1].colorbar(q)
                    norm = matplotlib.colors.Normalize(vmin=min(jf), vmax=max(jf), clip=True)
                    mapper = cm.ScalarMappable(norm=norm, cmap='cool')
                    time_color = np.array([(mapper.to_rgba(v)) for v in jf])
                    for xx, yy, ee, color in zip(xf, ff, uf, time_color):
                        axs[jj+1].plot(xx, yy, 'o', color=color)
                        axs[jj+1].errorbar(xx, yy, ee, lw=1, capsize=3, color=color)
                axs[jj+1].set_ylim(-diff,diff)
                axs[jj+1].tick_params(axis='x',labelsize=15)
                axs[jj+1].tick_params(axis='y',labelsize=15)
                axs[jj+1].set_ylabel('Flux (Normalized)',labelpad=30,fontsize=25,weight='normal')
                axs[jj+1].set_xlabel('Phase (Folded at {} d)'.format(round(bestpeaks[jj],3)),labelpad=15,fontsize=25,weight='normal')

        axs[0].set_title(title,fontsize=35,pad=15)

        fig.align_ylabels(axs)
        plt.savefig(out+title+kind+'.png',format='png',bbox_inches='tight')
        plt.clf()
        plt.close()



    bounds = [[2456990,2457600],[2457600,2457801],[2457801,2458200],[2458200,2460000]]

    plt.figure(figsize=(16,10))

    snum = 1
    colors = ['red','orange','green','blue']

    phases     = foldAt(x, periodfound, centralzero=True)
    sortIndi1  = np.argsort(phases);
    ff         = y[sortIndi1]
    uf         = u[sortIndi1]
    jf         = x[sortIndi1]
    xf         = phases[sortIndi1]

    for b in range(len(bounds)):

        bound = bounds[b]

        xfI   = xf[(jf > bound[0]) & (jf < bound[1])]
        ffI   = ff[(jf > bound[0]) & (jf < bound[1])]
        ufI   =  uf[(jf > bound[0]) & (jf < bound[1])]

        plt.scatter(xfI,ffI,color=colors[b],s=20,label='S'+str(snum))
        plt.errorbar(xfI,ffI,yerr=ufI,color=colors[b],ls='none')
        snum += 1

    plt.ylabel('RV (km/s)')
    plt.xlabel('Phase (Folded at {} d)'.format(round(periodfound,5)))
    plt.legend()
    plt.savefig(out+title+'_folded_seasonstogether.png')
    plt.clf()
    plt.close()


    fig, axs = plt.subplots(2,3,figsize=(16,10))
    n = 0; m = 0;
    snum = 1

    for b in range(len(bounds)):

        bound = bounds[b]

        xfI   = xf[(jf > bound[0]) & (jf < bound[1])]
        ffI   = ff[(jf > bound[0]) & (jf < bound[1])]
        ufI   =  uf[(jf > bound[0]) & (jf < bound[1])]
        jfI   =  jf[(jf > bound[0]) & (jf < bound[1])]

        try:
        	r1 = Time(np.min(jfI),format='jd')
        	r2 = Time(np.max(jfI),format='jd')
        	axs[m,n].scatter(xfI,ffI,color=colors[b],s=20,label='S'+str(snum))
        	axs[m,n].errorbar(xfI,ffI,yerr=ufI,color=colors[b],ls='none')
        	axs[m,n].set_title('S{}: {} - {}'.format(snum,r1.isot[:7],r2.isot[:7]))
        	axs[m,n].set_ylabel('RV (km/s)')
        	axs[m,n].set_xlim(0,1)
        	axs[m,n].set_xlabel('Phase (Folded at {} d)'.format(round(periodfound,5)))
        except ValueError:
        	fig.delaxes(axs[m,n])
        	pass

        snum += 1

        n += 1;
        if n == 3:
        	n = 0; m = 1;

    fig.delaxes(axs[-1,-1])
    axs[m,n].scatter(xf,ff,color='black',s=20,label='S'+str(snum))
    axs[m,n].errorbar(xf,ff,yerr=uf,color='black',ls='none')
    axs[m,n].set_title('All')
    axs[m,n].set_ylabel('RV (km/s)')
    axs[m,n].set_xlim(0,1)
    axs[m,n].set_xlabel('Phase (Folded at {} d)'.format(round(periodfound,5)))

    ys = axs[m,n].get_ylim()
    n = 0; m = 0;
    for b in range(len(bounds)):
    	axs[m,n].set_ylim(ys)
    	n += 1;
    	if n == 3:
    		n = 0; m = 1;
    plt.tight_layout()
    plt.savefig(out+title+'_folded_seasonsseparate.png')
    plt.clf()
    plt.close()

    return


def LS_fringe(x,y,u,title,out,Prot0=None,Porb0=None):

    x = x[np.isfinite(y)];     u = u[np.isfinite(y)];    y = y[np.isfinite(y)];

    freq0,pow0 = LombScargle(x,y,u).autopower(normalization='standard',maximum_frequency = 1)
    per0 = 1/freq0
    pow0 = pow0[(per0 < 20)]
    per0 = per0[(per0 < 20)]
    periodfound = per0[np.argmax(pow0)]

    faps = LombScargle(x,y,u).false_alarm_level([0.1,0.05,0.01])

    nmpd = 20
    peaks,extra = find_peaks(pow0,distance=50,height=0.4)

    bestpeaks = [j for _,j in sorted(zip(pow0[peaks],per0[peaks]),reverse=True)] #x is sorted by y

    for alias in [2., 3., 4., 5.]:
        for peak in bestpeaks:
            if abs(bestpeaks[0]*alias - peak) < .05:
                bestpeaks.remove(bestpeaks[0])
                break
            elif abs(bestpeaks[0]/alias - peak) < .05:
                bestpeaks.remove(peak)
                break
            else:
                pass
    for peak in bestpeaks[1:]:
        if abs(bestpeaks[0]- peak) < .05:
            bestpeaks.remove(peak)

    #plt.plot(per0,pow0)
    #plt.scatter(per0[peaks],pow0[peaks],color='red')
    #plt.show()

    if len(bestpeaks) == 0:
        fig,axs = plt.subplots(len(bestpeaks)+1,1,figsize=(12,8))

        plt.subplots_adjust(hspace=.25)
        axs.plot(per0,pow0,color='black')
        if Prot0 is not None:
            axs.axvline(Prot0,ls='-',color='red',alpha=.6)
            for alias in [1/5.,1/4.,1/3.,1/2.]:
                if Prot0*alias > 1:
                    ypred = pow0[(abs(per0-Prot0*alias) == np.min(abs(per0-Prot0*alias)))]
                    axs.vlines([Prot0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='red',alpha=.4)
        if Porb0 is not None:
            axs.axvline(Porb0,ls='-',color='blue',alpha=.6)
            for alias in [1/5.,1/4.,1/3.,1/2.]:
                if Porb0*alias > 1:
                    ypred = pow0[(abs(per0-Porb0*alias) == np.min(abs(per0-Porb0*alias)))]
                    axs.vlines([Porb0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='blue',alpha=.4)
        ys = axs.get_ylim()
        for fap in faps:
            if fap < ys[-1]:
                axs.axhline(fap,color='black',alpha=.25,ls='--')
        axs.tick_params(axis='x',labelsize=15)
        axs.tick_params(axis='y',labelsize=15)
        axs.set_xlabel('Period (Angstrom)',labelpad=15,fontsize=25,weight='normal')
        axs.set_ylabel('Power',labelpad=30,fontsize=25,weight='normal')
        plt.savefig(out+title+'.png',format='png',bbox_inches='tight')
        plt.clf()
        plt.close()
        return


    fig,axs = plt.subplots(len(bestpeaks)+1,1,figsize=(12,8+8*len(bestpeaks)))

    plt.subplots_adjust(hspace=.25)
    axs[0].plot(per0,pow0,color='black')
    if Prot0 is not None:
        axs[0].axvline(Prot0,ls='-',color='red',alpha=.6)
    if Porb0 is not None:
        axs[0].axvline(Porb0,ls='-',color='blue',alpha=.6)

    for alias in [1/5.,1/4.,1/3.,1/2.]:
        if Prot0 is not None:
            if Prot0*alias > 1:
                ypred = pow0[(abs(per0-Prot0*alias) == np.min(abs(per0-Prot0*alias)))]
                axs[0].vlines([Prot0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='red',alpha=.4)
        if Porb0 is not None:
            if Porb0*alias > 1:
                ypred = pow0[(abs(per0-Porb0*alias) == np.min(abs(per0-Porb0*alias)))]
                axs[0].vlines([Porb0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='blue',alpha=.4)

    ypreds = sorted(pow0[peaks],reverse=True)[:len(bestpeaks)]
    markers = ['x','d','s','*','^']
    for j in range(len(bestpeaks)):

        axs[0].scatter(bestpeaks[j],ypreds[j]*1.1,color='black',marker= markers[j] ,label='P = '+str(round(bestpeaks[j],3)))
        for alias in [1/5.,1/4.,1/3.,1/2.]:
            if bestpeaks[j]*alias > 1:
                ypred = pow0[(abs(per0-bestpeaks[j]*alias) == np.min(abs(per0-bestpeaks[j]*alias)))]
                axs[0].scatter(bestpeaks[j]*alias,ypred*1.1,color='black',marker=markers[j],s=5,alpha=.4)

    axs[0].set_xlabel('Period (Angstroms)',labelpad=15,fontsize=25,weight='normal')
    axs[0].set_ylabel('Power',labelpad=30,fontsize=25,weight='normal')
    axs[0].legend()
    axs[0].tick_params(axis='x',labelsize=15)
    axs[0].tick_params(axis='y',labelsize=15)

    ys = axs[0].get_ylim()
    for fap in faps:
        if fap < ys[-1]:
            axs[0].axhline(fap,color='black',alpha=.25,ls='--')

    if len(bestpeaks) > 0:

        for jj in range(len(bestpeaks)):
            xf,ff = folder(x,bestpeaks[jj],y)
            xf,uf = folder(x,bestpeaks[jj],u)
            diff = np.nanpercentile(ff,95)-np.nanpercentile(ff,5)
            axs[jj+1].errorbar(xf/bestpeaks[jj],ff,yerr=uf,color='black',ls='none',alpha=1)
            axs[jj+1].scatter(xf/bestpeaks[jj],ff,color='black',s=7,alpha=1)
            axs[jj+1].set_ylim(-diff,diff)
            axs[jj+1].tick_params(axis='x',labelsize=15)
            axs[jj+1].tick_params(axis='y',labelsize=15)
            axs[jj+1].set_ylabel('Flux (Normalized)',labelpad=30,fontsize=25,weight='normal')
            axs[jj+1].set_xlabel('Phase (Folded at {} d)'.format(round(bestpeaks[jj],3)),labelpad=15,fontsize=25,weight='normal')

    #for jj in range(len(bestpeaks)):
    #    axs[jj+1].legend()

    axs[0].set_title(title,fontsize=35,pad=15)

    fig.align_ylabels(axs)
    plt.savefig(out+title+'.png',format='png',bbox_inches='tight')
    plt.clf()
    plt.close()


def LS2(x,y,u,x2,y2,u2,label1,label2,title,out,Prot0=None,Porb0=None):

    freq0,pow0 = LombScargle(x,y,u).autopower(normalization='standard',maximum_frequency = 1)
    per0 = 1/freq0
    pow0 = pow0[(per0 < 20)]
    per0 = per0[(per0 < 20)]
    periodfound = per0[np.argmax(pow0)]

    faps = LombScargle(x,y,u).false_alarm_level([0.1,0.05,0.01])

    freq02,pow02 = LombScargle(x2,y2,u2).autopower(normalization='standard',maximum_frequency = 1)
    per02 = 1/freq02
    pow02 = pow02[(per02 < 20)]
    per02 = per02[(per02 < 20)]
    periodfound2 = per02[np.argmax(pow02)]

    faps2 = LombScargle(x2,y2,u2).false_alarm_level([0.1,0.05,0.01])

    fig,axs = plt.subplots(1,1,figsize=(12,8))

    plt.subplots_adjust(hspace=.25)
    axs.plot(per0,pow0,color='blue',alpha=.6,label=label1)
    axs.plot(per02,pow02,color='orange',alpha=.6,label=label2)
    if Prot0 is not None:
        axs.axvline(Prot0,ls='-',color='red',alpha=.1)
        for alias in [1/5.,1/4.,1/3.,1/2.]:
            if Prot0*alias > 1:
                ypred = pow0[(abs(per0-Prot0*alias) == np.min(abs(per0-Prot0*alias)))]
                axs.vlines([Prot0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='red',alpha=.1)
    if Porb0 is not None:
        axs.axvline(Porb0,ls='-',color='blue',alpha=.1)
        if Porb0*alias > 1:
            ypred = pow0[(abs(per0-Porb0*alias) == np.min(abs(per0-Porb0*alias)))]
            axs.vlines([Porb0*alias],ymin=ypred*1.1,ymax=ypred*1.2,color='blue',alpha=.1)
    ys = axs.get_ylim()
    for fap in faps:
        if fap < ys[-1]:
            axs.axhline(fap,color='black',alpha=.25,ls='--')
    axs.tick_params(axis='x',labelsize=15)
    axs.tick_params(axis='y',labelsize=15)
    axs.set_xlabel('Period (d)',labelpad=15,fontsize=25,weight='normal')
    axs.set_ylabel('Power',labelpad=30,fontsize=25,weight='normal')
    plt.legend()
    plt.savefig(out+title+'.png',format='png',bbox_inches='tight')
    plt.clf()
    plt.close()
