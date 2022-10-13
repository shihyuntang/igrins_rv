

from Engine.importmodule import *
from Engine.rebin_jv   import rebin_jv
from Engine.detect_peaks import detect_peaks
from scipy.optimize import curve_fit
import numpy as np

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

    def __init__(self,inpath, night, orders, jerp):
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
            self.flux  = np.array(tbdata['FLUX'],dtype=float)
            self.unc   = np.array(tbdata['UNC'],dtype=float)
            self.night = str(night)


# Get CC curve and clean and trim it
def peakfinderSpec(x,y,night,order,outpath,rv_sep,plot,tag,debug):

    x0 = x[np.argmax(y)]

    dist = 5*rv_sep

    xtofit = x[(x > x0-dist) & (x < x0+dist)]
    ytofit = y[(x > x0-dist) & (x < x0+dist)]

    xbute = np.linspace(xtofit[0],xtofit[-1],10000)

    minguess = min(ytofit)
    rvrange = 7

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



def bisector(rv,cc,rv_out,outpath,night,order,debug):

    points0 = np.concatenate((np.linspace(.4,.6,30),np.linspace(.7,.90,30)))

    valleys = detect_peaks(cc,valley=True)
    valleysL = valleys[(rv[valleys] < rv_out - 7) & (cc[valleys] < 0.2)]
    valleysR = valleys[(rv[valleys] > rv_out + 7) & (cc[valleys] < 0.2)]
    if len(valleysL) > 0:
        ccleft  = cc[(rv[valleysL][-1] < rv) & (rv < rv_out)]
    else:
        ccleft  = cc[(rv < rv_out)]
    if len(valleysR) > 0:
        ccright = cc[(rv[valleysR][0] > rv) & (rv > rv_out)]
    else:
        ccright = cc[(rv > rv_out)]

    lowbis = []; highbis = []; bis = []; points = []; errbis = [];

    for point in points0:
        if len(ccleft[(ccleft < point)]) == 0 or len(ccright[(ccright < point)]) == 0:
            continue
        if len(ccleft[(ccleft > point)]) == 0 or len(ccright[(ccright > point)]) == 0:
            continue

        points.append(point)
        ccleftover = np.min(ccleft[(ccleft > point)])
        ccleftunder = np.max(ccleft[(ccleft < point)])
        fleft = interp1d([ccleftunder,ccleftover],[rv[(cc == ccleftunder)][0],rv[(cc == ccleftover)][0]])

        ccrightover = np.min(ccright[(ccright > point)])
        ccrightunder = np.max(ccright[(ccright < point)])
        fright = interp1d([ccrightover,ccrightunder],[rv[(cc == ccrightover)][0],rv[(cc == ccrightunder)][0]])

        sloperight = (rv[(cc == ccrightover)][0]-rv[(cc == ccrightunder)][0])/(ccrightover-ccrightunder)
        #eright = np.mean([errs[(cc == ccrightover)],errs[(cc == ccrightunder)]])*sloperight
        slopeleft = (rv[(cc == ccleftover)][0]-rv[(cc == ccleftunder)][0])/(ccleftover-ccleftunder)
        #eleft = np.mean([errs[(cc == ccleftover)],errs[(cc == ccleftunder)]])*slopeleft

        if point < .5:
            lowbis.append((fright(point)+fleft(point))/2.0)
        else:
            highbis.append((fright(point)+fleft(point))/2.0)
        bis.append((fright(point)+fleft(point))/2.0)
        #errbis.append(np.mean([eleft,eright])/np.sqrt(2))

    if debug == True:
        plt.plot(rv,cc,alpha=1,color='blue')
        plt.scatter(rv[valleysL],cc[valleysL],color='blue',s=15)
        plt.scatter(rv[valleysR],cc[valleysR],color='blue',s=15)
        plt.scatter(bis,points,color='blue',s=15)
        plt.text(np.median(rv)+5,0.5,round(np.mean(lowbis)-np.mean(highbis),3))
        plt.savefig('{}/Bisector_{}_{}.png'.format(outpath,night,order))
        plt.clf()
        plt.close()

    bi = np.mean(lowbis)-np.mean(highbis)
    bistd = np.sqrt( (np.std(lowbis)/np.sqrt(len(lowbis)))**2 + (np.std(highbis)/np.sqrt(len(highbis)))**2 )
    #bistd = np.sqrt(  np.mean(errbis)**2 + (np.std(lowbis)/np.sqrt(len(lowbis)))**2 + (np.std(highbis)/np.sqrt(len(highbis)))**2 )

    return bi,bistd


def BIinst(refspec,orders,jerp,nights_use,args,inpath,outpath,logger,i):

    night = nights_use[i];    order = orders[jerp];

    if night == refspec.night:
        return night, np.nan, np.nan

    if np.int(night[:8]) < 20180401 or np.int(night[:8]) > 20190531:
        T_L = 'T'
    else:
        T_L = 'L'

    newspec = NightSpecs(inpath, night,    orders, jerp)

    if newspec.flag == 1:
        logger.warning('WARNING! {} order {} is flagged!'.format(newspec.filename,order))
        return night,np.nan,np.nan

    if args.debug:
        plt.figure(figsize=(16,10))
        plt.plot(newspec.wave,newspec.flux,color='black',alpha=.5)
        plt.plot(refspec.wave,refspec.flux,color='red',alpha=.5)
        plt.savefig('{}/Spec_{}_{}.png'.format(outpath,newspec.night,order))
        plt.clf()
        plt.close()


    # Rebin refspec to uniform wave scale
    dw = np.median(np.diff(refspec.wave))
    nw = int( (refspec.wave[-1] - refspec.wave[0])/dw )
    uniwave  = np.linspace(refspec.wave[0],refspec.wave[-1],nw)
    uniflux  = rebin_jv(refspec.wave,refspec.flux,uniwave,False)
    uniwave = uniwave[1:-1]
    uniflux = uniflux[1:-1]
    refspec.wave = uniwave; refspec.flux = uniflux;

    # Rebin spectra to same wavelength scale
    rvrange = 100.
    newspec.flux = rebin_jv(newspec.wave,newspec.flux,refspec.wave,False)

    # Now split treatment:
    refspec1 = deepcopy(refspec); refspec2 = deepcopy(refspec);
    newspec1 = deepcopy(newspec); newspec2 = deepcopy(newspec);

    # for pixel-wise correlation, cut both to be on same non-overinterpolated wave scale
    newspec2.flux = newspec2.flux[(refspec.wave > newspec.wave[0]) & (refspec.wave < newspec.wave[-1])]
    refspec2.flux = refspec2.flux[(refspec.wave > newspec.wave[0]) & (refspec.wave < newspec.wave[-1])]
    newspec2.wave = refspec2.wave[(refspec.wave > newspec.wave[0]) & (refspec.wave < newspec.wave[-1])]
    refspec2.wave = refspec2.wave[(refspec.wave > newspec.wave[0]) & (refspec.wave < newspec.wave[-1])]
    refspec2.flux[(refspec2.flux > 1.03)] = np.nan
    newspec2.flux[(newspec2.flux > 1.03)] = np.nan

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
    rv_out2,rv2,cc2 = peakfinderSpec(rv2,cc2,night,order,outpath,rv_sep,args.plotfigs,'uni',args.debug)

    #if abs(rv_out1-rv_out2) > 1:
    #    logger.info(f'WARNING! {night} order {order} shows CC RV disagreement of {round(abs(rv_out1-rv_out2),2)} km/s!')
    #    return night,np.nan,np.nan

    #rv_out = rv_out1
    #rv = rv1.copy(); cc = cc1.copy(); # ONLY USE PYA CC

    bi,bistd = bisector(rv2,cc2,rv_out2,outpath,night,order,args.debug)

    return night,bi,bistd
