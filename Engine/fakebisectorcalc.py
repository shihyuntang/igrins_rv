

from Engine.importmodule import *
from Engine.rebin_jv   import rebin_jv
from Engine.detect_peaks import detect_peaks
from scipy.interpolate import splrep,splev #, interp1d
from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from scipy.optimize import curve_fit
from Engine.rotint import rotint
from Engine.macbro_dynamic    import macbro_dyn
from Engine.rebin_jv import rebin_jv
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


class FakeSpec:

    def __init__(self,kind, T1, T2, pow1, pow2, rv1,rv2, s2n, night,order,jerp,IPnights1,IP2,IP3,logger):

        IP1 = IPnights1[night]

        watm,satm, mwave1, mflux1 = setup_templates(
            logger, 'synthetic', 'K', np.int(T1),
            np.float(4.5), np.float(0)
            )

        if 'RV' in kind:
            kind = kind[:5]

        if kind[0] == 'M':
            mininamedict = {'M2obs':'HD95735_M2_V_K','M3obs':'BD+443567_M3_V_K',
                            'M5obs':'BD-074003_M5_V_K','M6obs':'HD1326B_M6_V_K'}
            tempdat = Table.read(f'./Engine/spec_library/{mininamedict[kind]}.txt',format='ascii')
            mwave2   = np.array(tempdat.columns[0],dtype=float)*1e4
            mflux2   = np.array(tempdat.columns[1],dtype=float)
            #runc    = np.array(tempdat.columns[2],dtype=float)
            #rs2n    = mflux2/runc
            mflux2  /= np.median(mflux2)
            #runc    = mflux2/rs2n
        else:
            mininamedict = {'3500_4p0':'syntheticstellar_kband_T3500_logg4.0_0.0kG',
                            '3500_4p5':'syntheticstellar_kband_T3500_logg4.5_0.0kG',
                            '3000_4p0':'syntheticstellar_kband_T3000_logg4.0_0.0kG',
                            '4000_4p5':'syntheticstellar_kband_T4000_logg4.5_0.0kG',
                            '3000_4p0_phx':'PHOENIX-lte03000-4.00-0.0_contadjK',
                            '3000_4p5_phx':'PHOENIX-lte03000-4.50-0.0_contadjK'}
            tempdat = Table.read(f'./Engine/syn_template/{mininamedict[kind]}.txt',format='ascii')
            mwave2   = np.array(tempdat['wave'],dtype=float)
            mflux2   = np.array(tempdat['flux'],dtype=float)
            # Remove duplicate wavelength values from stellar template
            # (doesn't affect steps 1-3, but needed for bisectors)
            ind = []
            maxwave = mwave2[0]
            for j in range(1,len(mwave2)-1):
                if mwave2[j] > maxwave:
                    maxwave = mwave2[j]
                else:
                    ind.append(j)
            ind = np.array(ind)
            mask = np.ones(len(mwave2), dtype=bool)
            if len(ind) > 0:
                mask[ind] = False
                mwave2 = mwave2[mask]
                mflux2 = mflux2[mask]

        vsini = 12.3774
        fluxratio = 1/8.

        filename = './Output/DITau_K/RV_results_2/Cutouts/Cutout_{}Combined.fits'.format(night)
        hdu = fits.open(filename)
        tbdata = hdu[jerp+1].data
        w  = np.array(tbdata['WAVE'+str(order)],dtype=float)
        c = 2.99792458e5
        npts = len(w)

        if np.median(np.diff(mwave1)) > np.median(np.diff(mwave2)):
            rebin2to1 = True; extra1 = 0.; extra2 = 10.;
        else:
            rebin2to1 = False; extra1 = 10.; extra2 = 0.;

        mflux1 = mflux1[(mwave1 > w[0]-30-extra1) & (mwave1 < w[-1]+30+extra1)]
        mwave1 = mwave1[(mwave1 > w[0]-30-extra1) & (mwave1 < w[-1]+30+extra1)]
        mflux2 = mflux2[(mwave2 > w[0]-30-extra2) & (mwave2 < w[-1]+30+extra2)]
        mwave2 = mwave2[(mwave2 > w[0]-30-extra2) & (mwave2 < w[-1]+30+extra2)]

        # Apply velocity shifts and scale
        mwave1 = mwave1*(1.+rv1/c)
        mflux1 = mflux1**pow1
        mwave2 = mwave2*(1.+rv2/c)
        mflux2 = mflux2**pow2

        wspot1,rspot1 = rotint(mwave1, mflux1, vsini)
        wspot2,rspot2 = rotint(mwave2, mflux2, vsini)

        rspot2   *=  fluxratio

        matplotlib.use('Qt5Agg')

        if rebin2to1:
            rspot2n = rebin_jv(wspot2,rspot2,wspot1,False)
            rspot   = rspot1 + rspot2n
            wspot   = wspot1.copy()
        else:
            rspot1n = rebin_jv(wspot1,rspot1,wspot2,False)
            rspot   = rspot2 + rspot1n
            wspot   = wspot2.copy()

        rspot /= (1+fluxratio)

        #Find mean observed wavelength and create a telluric velocity scale
        mnw = np.mean(w)
        dw = (w[-1] - w[0])/(npts-1.)
        vel = (wspot-mnw)/mnw*c

        xdict = {4:np.arange(159,940),5:np.arange(164,1477),6:np.arange(180,1842)}
        x = xdict[order]
        if len(x) > len(w):
            x = x[:len(w)-len(x)]

        fwhmraw = IP1 + IP2*(x) + IP3*(x**2)
        spl = splrep(w, fwhmraw)
        fwhm = splev(wspot,spl)

        # Have IP extend as constant past wave bounds of data
        fwhm[(wspot < w[0])]  = fwhm[(wspot >= w[0])][0]
        fwhm[(wspot > w[-1])] = fwhm[(wspot <= w[-1])][-1]
        if (np.min(fwhm) < 1) or (np.max(fwhm) > 8):
            sys.exit('IP Error!')

        #Handle instrumental broadening
        vhwhm = dw*np.abs(fwhm)/mnw*c/2.
        nsmod = macbro_dyn(vel, rspot, vhwhm)

        #Rebin model to observed wavelength scale
        smod = rebin_jv(wspot, nsmod ,w, False)

        '''
        plt.figure(figsize=(12,10))
        plt.plot(wspot,rspot,color='red',alpha=0.5)
        plt.plot(wspot1,rspot1,color='blue',alpha=0.5)
        plt.plot(w,smod,color='black',alpha=0.5)
        plt.show()
        '''

        sout = np.random.normal(smod,smod/s2n)

        self.x  = x
        self.wave  = w
        self.flux  = sout
        self.night = night
        self.unc   = sout/s2n
        self.flag  = 0



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


def BIinst(orders,jerp,nights_use,outpath,logger,fakeargs,rvs1in,rvs2in,i):

    night = nights_use[i];    order = orders[jerp];

    rv1in = rvs1in[i]
    rv2in = rvs2in[i]

    newspec = FakeSpec(fakeargs[0],fakeargs[1],fakeargs[2],fakeargs[3],fakeargs[4],
                        rv1in,rv2in,fakeargs[5],night,order,jerp,fakeargs[6],
                        fakeargs[7],fakeargs[8],logger)

    if night == '20170928':
        refspec = FakeSpec(fakeargs[0],fakeargs[1],fakeargs[2],fakeargs[3],fakeargs[4],
                        rvs1in[4],rvs2in[4], fakeargs[5], '20171122',order,jerp,fakeargs[6],
                        fakeargs[7],fakeargs[8],logger)
    else:
        refspec = FakeSpec(fakeargs[0],fakeargs[1],fakeargs[2],fakeargs[3],fakeargs[4],
                        rvs1in[3],rvs2in[3], fakeargs[5], '20170928',order,jerp,fakeargs[6],
                        fakeargs[7],fakeargs[8],logger)

    if np.int(night[:8]) < 20180401 or np.int(night[:8]) > 20190531:
        T_L = 'T'
    else:
        T_L = 'L'



    if newspec.flag == 1:
        logger.warning('WARNING! {} order {} is flagged!'.format(newspec.filename,order))
        return night,np.nan,np.nan

    if True:
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

    #rv_out0,rv0,cc0 = peakfinderSpec(rv10,cc10,night,order,outpath,rv_sep,True,'pyAwide',True)
    #rv_out1,rv1,cc1 = peakfinderSpec(rv1,cc1,night,order,outpath,rv_sep,True,'pyA',True)
    rv_out2,rv2,cc2 = peakfinderSpec(rv2,cc2,night,order,outpath,rv_sep,True,'uni',True)

    #if abs(rv_out1-rv_out2) > 1:
    #    logger.info(f'WARNING! {night} order {order} shows CC RV disagreement of {round(abs(rv_out1-rv_out2),2)} km/s!')
    #    return night,np.nan,np.nan

    #rv_out = rv_out1
    #rv = rv1.copy(); cc = cc1.copy(); # ONLY USE PYA CC

    bi,bistd = bisector(rv2,cc2,rv_out2,outpath,night,order,True)

    return night,bi,bistd
