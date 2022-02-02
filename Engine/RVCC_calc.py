

from Engine.importmodule import *
from Engine.rebin_jv   import rebin_jv
from Engine.detect_peaks import detect_peaks
from scipy.optimize import curve_fit
import numpy as np
import glob

c = 2.99792458e5


# Gaussian fit convience function
def gauss_fit(x):
    def innerfit(*p):
        #print('p',p)
        mu = p[1]
        sigma = p[2]
        offset = p[3]
        scale = p[4]
        slope = p[5]
        kurt = p[6]
        #print('went')
        return(offset + slope*x + kurt*(x**2) + ((scale*np.sqrt(2*np.pi))*np.exp(-.5*((x-mu)/sigma)**2)))
    return innerfit

# Gaussian convience function
def gauss(x,mu,sigma,offset,scale,slope,kurt):
    return(offset+ slope*x + kurt*(x**2) +((scale*np.sqrt(2*np.pi))*np.exp(-.5*((x-mu)/sigma)**2)))


# Cross-correlation function
def c_correlate(x,y,lag):
    x = np.array(x,dtype=float); y = np.array(y,dtype=float);
    nX = len(x)
    Xd = x-np.nansum(x)/nX
    Yd = y-np.nansum(y)/nX
    nLag = len(lag)
    Cross = np.ones(nLag,dtype=float)
    for k in range(nLag):
        if lag[k] >= 0.:
            Cross[k] = np.nansum(Xd[0:nX - lag[k]] * Yd[lag[k]:])
        else:
            Cross[k] = np.nansum(Yd[0:nX + lag[k]] * Xd[-lag[k]:])
    return Cross/(np.sqrt(np.nansum(Xd**2)*np.nansum(Yd**2)))


# Cross-correlation deployment
def correlator(flux,fluxT,rv_sep):

    c = 2.99792458e5

    lags = np.arange(-100,100)
    cc = c_correlate(flux,fluxT,lags)

    rvs = lags*rv_sep

    lags = lags[(abs(rvs) < 200)]
    cc = cc[(abs(rvs) < 200)]
    rvs = rvs[(abs(rvs) < 200)]

    return rvs,cc


# Convenience class for storing stellar residual data
class NightSpecs:

    def __init__(self,inpath, night, orders, jerp,args):
        self.filename = '{}/Cutout_{}Combined.fits'.format(inpath, night)
        hdu = fits.open(self.filename)
        tbdata = hdu[jerp+1].data
        order = orders[jerp]
        self.flag = 1
        try:
            self.flag  = np.array(tbdata['ERRORFLAG'+str(order)])[0]
        except KeyError:
            for nnn in range(1,10):
                try:
                    tbdata = hdu[jerp+1+nnn].data
                    self.flag  = np.array(tbdata['ERRORFLAG'+str(order)])[0]
                except:
                    continue
                break

        if self.flag == 0:
            self.wave  = np.array(tbdata['WAVE'+str(order)],dtype=float)
            if args.component == 'primary':
                self.flux  = np.array(tbdata['FLUX'],dtype=float)
            else:
                self.flux  = np.array(tbdata['RESID'],dtype=float)
            self.unc   = np.array(tbdata['UNC'],dtype=float)
            self.night = str(night)


class NightTagSpecs:

    def __init__(self,inpath, night, orders, tag, jerp,args):
        self.filename = '{}/Cutout_{}_{}.fits'.format(inpath, night,tag)
        hdu = fits.open(self.filename)
        tbdata = hdu[jerp+1].data
        order = orders[jerp]
        self.flag = 1
        try:
            self.flag  = np.array(tbdata['ERRORFLAG'+str(order)])[0]
        except KeyError:
            for nnn in range(1,10):
                try:
                    tbdata = hdu[jerp+1+nnn].data
                    self.flag  = np.array(tbdata['ERRORFLAG'+str(order)])[0]
                except:
                    continue
                break

        if self.flag == 0:
            self.wave  = np.array(tbdata['WAVE_ADJ'+str(order)],dtype=float)
            if args.component == 'primary':
                self.flux  = np.array(tbdata['FLUX_CORR'],dtype=float)
            else:
                self.flux  = np.array(tbdata['FLUX_RESID'],dtype=float)
            S2N        = np.array(tbdata['S2N'],dtype=float)
            self.unc   = self.flux/S2N
            self.flux = self.flux[self.wave!=0]
            self.unc = self.unc[self.wave!=0]
            self.wave = self.wave[self.wave!=0]
            self.night = str(night)

class FakeNightSpecs:

    def __init__(self,wave,flux):
        self.wave  = wave
        self.flux  = flux
        self.night = 'template'
        self.flag  = 0

# Get CC curve and clean and trim it
def peakfinderSpec(x,y,night,order,outpath,rv_sep,plot,tag,debug):

    x0 = x[np.argmax(y)]

    dist = 5*rv_sep

    xtofit = x[(x > x0-dist) & (x < x0+dist)]
    ytofit = y[(x > x0-dist) & (x < x0+dist)]

    xbute = np.linspace(xtofit[0],xtofit[-1],10000)

    minguess = min(ytofit)
    rvrange = 2

    # mu,sigma,offset,scale,slope,kurt
    init =  [x0,                4., minguess,    1e-1,  1e-6,   1e-4]
    low  =  [x0-rvrange*rv_sep, .1, minguess-.5, 1e-5, -1e-1,  -1e-1]
    high =  [x0+rvrange*rv_sep, 30.,minguess+.5, 1e2,   1e-1,   1e-1]

    popt, pcov = curve_fit(gauss_fit(xtofit), xtofit,ytofit, p0=init, bounds=(low,high),
                           sigma = np.ones_like(xtofit)*.001,ftol=1e-3)

    fitcurve = gauss(xbute,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    rv_out = xbute[np.argmax(fitcurve)]

    valleys = detect_peaks(fitcurve,valley=True)
    valleys = valleys[(10000-250 > valleys) & (valleys > 250)]

    if len(valleys) > 0: # do fitting again but expand to include more in direction(s) of inversion(s)
        if np.min(valleys) < 5000:
            distL = 8*rv_sep
        else:
            distL = 5*rv_sep
        if np.max(valleys) > 5000:
            distR = 8*rv_sep
        else:
            distR = 5*rv_sep

        xtofit = x[(x > x0-distL) & (x < x0+distR)]
        ytofit = y[(x > x0-distL) & (x < x0+distR)]

        xbute = np.linspace(xtofit[0],xtofit[-1],10000)
        minguess = min(ytofit)
        rvrange = 7

        init =  [x0,4.,minguess,1e-1,1e-6,1e-4]
        low = [x0-rvrange*rv_sep,.1,minguess-.5,1e-5,-1e-1,-1e-1]
        high = [x0+rvrange*rv_sep,30.,minguess+.5,1e2,1e-1,1e-1]
        popt, pcov = curve_fit(gauss_fit(xtofit), xtofit,ytofit, p0=init, bounds=(low,high),
                           sigma = np.ones_like(xtofit)*.001,ftol=1e-3)

        fitcurve = gauss(xbute,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
        rv_out = xbute[np.argmax(fitcurve)]

    if debug == True:
        plt.figure(figsize=(16,10))
        plt.plot(x,y,color='black',alpha=0.5)
        plt.plot(xbute,fitcurve,color='red',alpha=.5)
        plt.scatter(xbute[valleys],fitcurve[valleys],color='blue',s=60)
        #plt.axvline(popt[0])
        plt.axvline(rv_out)
        plt.savefig('{}/Crosscorr_{}_{}_{}.png'.format(outpath,night,order,tag))
        plt.clf()
        plt.close()

    dist = 20*rv_sep
    rv = x[(x > x0-dist-50) & (x < x0+dist+50)]
    cc = y[(x > x0-dist-50) & (x < x0+dist+50)]

    cc -= np.min(cc)
    cc /= np.max(cc)

    return rv_out, rv, cc



def RVinst(refspec,orders,jerp,nights_use,args,inpath,outpath,logger,kind,i):

    night = nights_use[i];    order = orders[jerp];

    if kind == 'self':
        if night == refspec.night:
            return night, np.nan, np.nan
    else:
        if kind[0] == 'M':
            mininamedict = {'M2obs':'HD95735_M2_V_K','M3obs':'BD+443567_M3_V_K',
                            'M5obs':'BD-074003_M5_V_K','M6obs':'HD1326B_M6_V_K'}
            tempdat = Table.read(f'./Engine/spec_library/{mininamedict[kind]}.txt',format='ascii')
            rwave   = np.array(tempdat.columns[0],dtype=float)*1e4
            rflux   = np.array(tempdat.columns[1],dtype=float)
            #runc    = np.array(tempdat.columns[2],dtype=float)
            #rs2n    = rflux/runc
            rflux  /= np.median(rflux)
            #runc    = rflux/rs2n
        else:
            mininamedict = {'3500_4p0':'syntheticstellar_kband_T3500_logg4.0_0.0kG',
                            '3500_4p5':'syntheticstellar_kband_T3500_logg4.5_0.0kG',
                            '3000_4p0':'syntheticstellar_kband_T3000_logg4.0_0.0kG',
                            '4000_4p5':'syntheticstellar_kband_T4000_logg4.5_0.0kG',
                            '3000_4p0_phx':'PHOENIX-lte03000-4.00-0.0_contadjK',
                            '3000_4p5_phx':'PHOENIX-lte03000-4.50-0.0_contadjK'}
            tempdat = Table.read(f'./Engine/syn_template/{mininamedict[kind]}.txt',format='ascii')
            rwave   = np.array(tempdat['wave'],dtype=float)
            rflux   = np.array(tempdat['flux'],dtype=float)
            # Remove duplicate wavelength values from stellar template
            # (doesn't affect steps 1-3, but needed for bisectors)
            ind = []
            maxwave = rwave[0]
            for j in range(1,len(rwave)-1):
                if rwave[j] > maxwave:
                    maxwave = rwave[j]
                else:
                    ind.append(j)
            ind = np.array(ind)
            mask = np.ones(len(rwave), dtype=bool)
            if len(ind) > 0:
                mask[ind] = False
                rwave = rwave[mask]
                rflux = rflux[mask]

        refspec  =  FakeNightSpecs(rwave,rflux)


    if np.int(night[:8]) < 20180401 or np.int(night[:8]) > 20190531:
        T_L = 'T'
    else:
        T_L = 'L'

    rvtags = []

    if args.ABkind == 'separate':
        filelist = glob.glob(f'{inpath}/Cutout_{night}_*.fits')
    else:
        filelist = [f'{inpath}/Cutout_{night}Combined.fits']

    for ff in filelist:

        if args.ABkind == 'separate':
            tag = ff[-9:-5]
            newspec = NightTagSpecs(inpath, night,    orders, tag, jerp, args)
        else:
            newspec = NightSpecs(inpath, night,    orders, jerp, args)
            tag = 'combined'

        if newspec.flag == 1:
            logger.warning('WARNING! {} order {} is flagged!'.format(newspec.filename,order))
            return night,np.nan,np.nan

        # Rebin refspec to uniform wave scale
        dw = np.median(np.diff(newspec.wave))
        nw = int( (newspec.wave[-1] - newspec.wave[0])/dw )
        uniwave  = np.linspace(newspec.wave[0],newspec.wave[-1],nw)
        uniflux  = rebin_jv(newspec.wave,newspec.flux,uniwave,False)
        uniwave = uniwave[1:-1]
        uniflux = uniflux[1:-1]
        newspec.wave = uniwave; newspec.flux = uniflux;


        # Rebin spectra to same wavelength scale
        rvrange = 100.

        # Now split treatment:
        refspec1 = deepcopy(refspec); refspec2 = deepcopy(refspec);
        newspec1 = deepcopy(newspec); newspec2 = deepcopy(newspec);
        newflux = rebin_jv(refspec2.wave,refspec2.flux,newspec2.wave,False)


        # for pixel-wise correlation, cut both to be on same non-overinterpolated wave scale
        refspec2.flux  =       newflux[(newspec.wave > refspec.wave[0]) & (newspec.wave < refspec.wave[-1])]
        newspec2.flux  = newspec2.flux[(newspec.wave > refspec.wave[0]) & (newspec.wave < refspec.wave[-1])]
        refspec2.wave  = newspec2.wave[(newspec.wave > refspec.wave[0]) & (newspec.wave < refspec.wave[-1])]
        newspec2.wave  = newspec2.wave[(newspec.wave > refspec.wave[0]) & (newspec.wave < refspec.wave[-1])]

        refspec2.flux[(refspec2.flux > 1.03)] = np.nan
        newspec2.flux[(newspec2.flux > 1.03)] = np.nan

        if args.debug:
            plt.figure(figsize=(16,10))
            plt.plot(newspec2.wave,newspec2.flux,color='black',alpha=.5)
            plt.plot(refspec2.wave,refspec2.flux,color='red',alpha=.5)
            plt.savefig('{}/Spec_{}_{}_{}.png'.format(outpath,newspec.night,order,tag))
            plt.clf()
            plt.close()


        ### COMMENTED OUT ALL THE PYASTRO STUFF

        # For wavelength-wise correlation, cut edges of one spectrum for shifts from the other
        #When template left edge is shifted 100 km/s to the right, it must still be less than left edge of data
        #When template right edge is shifted 100 km/s to the left, it must still be more than right edge of data
        #newspec1.flux = newspec1.flux[(refspec1.wave > newspec1.wave[0]) & (refspec1.wave < newspec1.wave[-1])]
        #newspec1.wave = refspec1.wave[(refspec1.wave > newspec1.wave[0]) & (refspec1.wave < newspec1.wave[-1])]
        #newspec1.flux = newspec1.flux[(newspec1.wave > refspec1.wave[0]*(1+rvrange/c)) & (newspec1.wave < refspec1.wave[-1]*(1-rvrange/c))]
        #newspec1.wave = newspec1.wave[(newspec1.wave > refspec1.wave[0]*(1+rvrange/c)) & (newspec1.wave < refspec1.wave[-1]*(1-rvrange/c))]

        rv_sep  = c*np.median(np.diff(newspec2.wave))/np.median(newspec2.wave)

        #rv10, cc10 = crosscorrRV( newspec1.wave, newspec1.flux, refspec1.wave, refspec1.flux,  -100, 100, 0.01, skipedge=0)
        #rvguess = rv10[np.argmax(cc10)]
        #rv1, cc1 = crosscorrRV(  newspec1.wave, newspec1.flux, refspec1.wave, refspec1.flux, rvguess-10, rvguess+10, 0.001, skipedge=0)
        rv2, cc2 = correlator(refspec2.flux,newspec2.flux,rv_sep)

        #rv_out0,rv0,cc0 = peakfinderSpec(rv10,cc10,night,order,outpath,rv_sep,args.plotfigs,'pyAwide',args.debug)
        #rv_out1,rv1,cc1 = peakfinderSpec(rv1,cc1,night,order,outpath,rv_sep,args.plotfigs,'pyA',args.debug)
        rv_out2,rv2,cc2 = peakfinderSpec(rv2,cc2,night,order,outpath,rv_sep,args.plotfigs,tag,args.debug)

        #if abs(rv_out1-rv_out2) > 1:
        #    logger.info(f'WARNING! {night} order {order} shows CC RV disagreement of {round(abs(rv_out1-rv_out2),2)} km/s!')
        #    return night,np.nan,np.nan

        #rv_out = rv_out1
        #rv = rv1.copy(); cc = cc1.copy(); # ONLY USE PYA CC
        rvtags.append(rv_out2)

    rvtags = np.array(rvtags)

    return night,np.nanmean(rvtags),np.nanstd(rvtags)/np.sqrt(len(rvtags[np.isfinite(rvtags)]))
