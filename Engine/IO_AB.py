
import numpy as np
import os

from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord, solar_system, EarthLocation, ICRS
from astropy import units
import os
from os import listdir
from os.path import isfile, join, isdir
import re
import sys

def partial_loader(inpath0,order):

    hdulist = fits.open(inpath0+'.spec.fits')
    tbdata_w = hdulist[1].data
    tbdata_f = hdulist[0].data

    hdulist_sn = fits.open(inpath0+'.sn.fits')
    tbdata_sn = hdulist_sn[0].data

    wavelist = np.array(tbdata_w[order],dtype=np.float64)
    fluxlist = np.array(tbdata_f[order],dtype=np.float64)
    s2nlist = np.array(tbdata_sn[order],dtype=np.float64)

    return wavelist,fluxlist,s2nlist


def init_fitsread(path,kind,beam,night,order,tag,band,Ncuts=None):

    if beam not in ['combined','separate']:
        print('beam MUST BE "combined" OR "separate",  FORCE QUITTING!')
        print(breaker)

    if kind not in ['A0','target']:
        print('kind MUST BE "target" OR "A0", FORCE QUITTING!')
        print(breaker)

    # Initguesser takes combined AB versions of target spectra, of which there may be multiple, and combines /those/ together.
    if beam == 'combined':
        if kind != 'target':
            sys.exit('IO_AB ERROR: KIND SHOULD BE TARGET, THIS ERROR SHOULD ONLY THROW FROM INITGUESSER')


        subpath        = '{}{}/AB/'.format(path, night)
        fullpathprefix = '{}SDC{}_{}_'.format(subpath, band, night[:8])

        onlyfiles = [f for f in listdir(subpath) if isfile(join(subpath, f))]
        contmeds = []
        for f in onlyfiles:
            q = re.split('_', f)
            if q[0] != 'SDC{}'.format(band):
                continue
            qr = re.split('\.', q[2])
            tag = qr[0]

            # Readin fits
            wavelist0, fluxlist0, s2nlist0 = partial_loader(fullpathprefix+tag,order)

            if Ncuts != None:
                Nstartcut = Ncuts[0]; Nendcut = Ncuts[-1];
                ia = np.ones_like(wavelist0,dtype=bool)
                ia[0:Nstartcut] = False; ia[-Nendcut:len(ia)] = False;
                wavelist1 = wavelist0[ia]; fluxlist1 = fluxlist0[ia]; s2nlist1 = s2nlist0[ia];
                xlist = np.arange(len(wavelist1),dtype=float) + float(Nstartcut)
            else:
                wavelist1 = wavelist0; fluxlist1 = fluxlist0; s2nlist1 = s2nlist0;
                xlist = np.arange(len(wavelist1),dtype=float)

            try:
                fluxstack = np.vstack((fluxstack,fluxlist1))
                wavestack = np.vstack((wavestack,wavelist1))
                s2nstack    = np.vstack((s2nstack,s2nlist1))
            except UnboundLocalError:
                fluxstack = fluxlist1.copy()
                wavestack = wavelist1.copy()
                s2nstack     = s2nlist1.copy()

        try:
            wavelist = np.array([np.nanmean(wavestack[:,i]) for i in range(len(wavelist1))])
            fluxlist = np.array([np.nanmean(fluxstack[:,i]) for i in range(len(wavelist1))])
            s2nlist = np.array([np.sqrt(np.nansum(s2nstack[:,i]**2)) for i in range(len(wavelist1))])
            #ulist    = np.array([1/np.sqrt(np.nansum(1/(ustack[:,i]**2))) for i in range(len(wavelist1))])
        except UnboundLocalError:
            return 0,0,0,0

    # For A0 loading, use combined AB and load from AB/ folder. Shouldn't be more than one AB file, so shouldn't need to do any extra combining, unlike above.
    # For RVProc procedure, use As and Bs separately and load from A/ or B/ folders.
    else: # if separate beams
        if kind == 'target':
            pass
        else:
            path = '{}std/{}/AB/'.format(path, night)

        wavelist0,fluxlist0,s2nlist0 = partial_loader('{}SDC{}_{}_{}'.format(path, band, night[:8], tag), order)

        if Ncuts != None:
            Nstartcut = Ncuts[0]; Nendcut = Ncuts[-1];
            ia = np.ones_like(wavelist0,dtype=bool)
            ia[0:Nstartcut] = False; ia[-Nendcut:len(ia)] = False;
            wavelist = wavelist0[ia]; fluxlist = fluxlist0[ia]; s2nlist = s2nlist0[ia];
            xlist = np.arange(len(wavelist),dtype=float) + float(Nstartcut)
        else:
            wavelist = wavelist0; fluxlist = fluxlist0; s2nlist = s2nlist0;
            xlist = np.arange(len(wavelist),dtype=float)

    MAD = np.nanmedian(abs(np.nanmedian(s2nlist)-s2nlist))
    x3    = xlist[   (np.isnan(fluxlist) == False) & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    wave3 = wavelist[(np.isnan(fluxlist) == False) & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    s2n3    = s2nlist[   (np.isnan(fluxlist) == False) & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    s3    = fluxlist[(np.isnan(fluxlist) == False) & (np.isnan(s2nlist) == False) & (fluxlist > 0)]

    # Cut absurdities that will mess with fit
    #wave = wave3[(10 < s2n3) & (s2n3 < np.nanmedian(s2n3)+MAD*5)]
    #s = s3[(10 < s2n3) & (s2n3 < np.nanmedian(s2n3)+MAD*5)]
    #x = x3[(10 < s2n3) & (s2n3 < np.nanmedian(s2n3)+MAD*5)]
    #s2n = s2n3[(10 < s2n3) & (s2n3 < np.nanmedian(s2n3)+MAD*5)]
    wave = wave3[(s2n3 < np.nanmedian(s2n3)+MAD*5)]
    s = s3[(s2n3 < np.nanmedian(s2n3)+MAD*5)]
    x = x3[(s2n3 < np.nanmedian(s2n3)+MAD*5)]
    s2n = s2n3[(s2n3 < np.nanmedian(s2n3)+MAD*5)]

    if Ncuts == None:
        wave = wave[5:-5]
        s = s[5:-5]
        s2n = s2n[5:-5]
        x = x[5:-5]

    u = s/s2n
    u[(s2n < 10)] = 1e3*max(u)

    return x,wave,s,u

def airtovac(wave):
    """Converts air wavelengths to vaccuum wavelengths returns a float or array of wavelengths
    INPUTS:
      wave - Wavelengths in air, in Angstroms, float or array
    OUTPUTS:
      newwave - Wavelengths in a vacuum, in Angstroms, float or array
    NOTES:
      1. The procedure uses the IAU standard conversion formula
         from Morton (1991 Ap. J. Suppl. 77, 119) """
    wave = np.array(wave,float)
    # Converting to Wavenumber squared
    sigma2 = (1.0e4/wave)**2
    # Computing conversion factor
    fact = 1.0 + 6.4328e-5 + 2.94981e-2/(146.0 - sigma2) + 2.5540e-4/( 41.0 - sigma2)
    fact = fact*(wave >= 2000.0) + 1.0*(wave < 2000.0)
    newwave = wave*fact
    return newwave

def setup_templates(logger, kind='synthetic', band='K', sptype='M'):
    if (kind == 'synthetic') and (band == 'K'):
        if sptype not in ['K','M']:
            sys.exit('Pipeline does not have a stellar template for early type stars in K band! Upload your own?')
        logger.info('Using synthetic stellar template...')
        logger.info('!!!!!!!! INTERNAL TEST!!!! T5000 logg4.5!!!!!')
        if os.getcwd()[-1]=='v':
            stelldata = Table.read('./Engine/syntheticstellar_kband_T5000_logg4.5.txt',format='ascii')
            # stelldata = Table.read('./Engine/PHOENIX-lte03700-4.50-0.0_contadj.txt',format='ascii')
        else:
            # stelldata = Table.read('../Engine/PHOENIX-lte03700-4.50-0.0_contadj.txt',format='ascii')
            stelldata = Table.read('../Engine/syntheticstellar_kband_T5000_logg4.5.txt',format='ascii')
        mwave0 = np.array(stelldata['wave'])#*10000.0
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0
        mwave0 = airtovac(mwave0)
    elif kind == 'livingston' and band == 'K':
        if sptype not in ['K','M']:
            sys.exit('Pipeline does not have a stellar template for early type stars in K band! Upload your own?')
        logger.info('Using sunspot for stellar template...')
        if os.getcwd()[-1]=='v':
            stelldata = Table.read('./Engine/SpotAtl_contadjusted.txt',format='ascii')
        else:
            stelldata = Table.read('../Engine/SpotAtl_contadjusted.txt',format='ascii')
        mwave0 = np.array(stelldata['wave'])*10000.0
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0
    elif kind == 'synthetic' and band == 'H':
        if sptype not in ['F','G','K']:
            sys.exit('Pipeline does not have a stellar template for late type stars in H band! Upload your own?')
        logger.info('Using synthetic stellar template...')
        logger.info('!!!!!!!! INTERNAL TEST!!!! T4000 logg4.5!!!!!')
        if os.getcwd()[-1]=='v':
            stelldata = Table.read('./Engine/syntheticstellar_hband_T4000_logg4.5.txt',format='ascii')
        else:
            stelldata = Table.read('../Engine/syntheticstellar_hband_T4000_logg4.5.txt',format='ascii')
        mwave0 = np.array(stelldata['wave'])
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0
        mwave0 = airtovac(mwave0)
    elif kind == 'livingston' and band == 'H':
        if sptype not in ['F','G','K']:
            sys.exit('Pipeline does not have a stellar template for late type stars in H band! Upload your own?')
        logger.info('Using quiet sun for stellar template...')
        if os.getcwd()[-1]=='v':
            spotdata = Table.read('./Engine/PhotoAtl_Solar_contadjusted.txt',format='ascii')
        else:
            spotdata = Table.read('../Engine/PhotoAtl_Solar_contadjusted.txt',format='ascii')
        mwave0 = np.array(spotdata['wave'])*10000.0
        mflux0 = np.array(spotdata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0

    if os.getcwd()[-1]=='v':
        telluricdata = Table.read('./Engine/PhotoAtl Organized.txt',format='ascii')
    else:
        telluricdata = Table.read('../Engine/PhotoAtl Organized.txt',format='ascii')
    watm = np.array(telluricdata['wave'])*10000.0
    satm = np.array(telluricdata['flux'])
    watm = watm[(np.isfinite(satm))]
    satm = satm[(np.isfinite(satm))]
    satm[(satm < 0)] = 0
    return watm, satm, mwave0, mflux0


def setup_templates_tel():

    if os.getcwd()[-1]=='v':
        spotdata = Table.read('./Engine/SpotAtl Organized.txt',format='ascii')
    else:
        spotdata = Table.read('../Engine/SpotAtl Organized.txt',format='ascii')

    mwave0 = np.array(spotdata['wave'])*10000.0
    mflux0 = np.array(spotdata['flux'])
    mwave0 = mwave0[(np.isfinite(mflux0))]
    mflux0 = mflux0[(np.isfinite(mflux0))]
    mflux0[(mflux0 < 0)] = 0

    if os.getcwd()[-1]=='v':
        telluricdata = Table.read('./Engine/PhotoAtl Organized.txt',format='ascii')
    else:
        telluricdata = Table.read('../Engine/PhotoAtl Organized.txt',format='ascii')

    watm = np.array(telluricdata['wave'])*10000.0
    satm = np.array(telluricdata['flux'])
    watm = watm[(np.isfinite(satm))]
    satm = satm[(np.isfinite(satm))]
    satm[(satm < 0)] = 0
    return watm, satm, mwave0, mflux0


def stellarmodel_setup(wave,mwave0,mflux0):
    mflux = mflux0[(mwave0/1e4 >= min(wave) - .003) & (mwave0/1e4 <= max(wave) + .002)]
    mwave = mwave0[(mwave0/1e4 >= min(wave) - .003) & (mwave0/1e4 <= max(wave) + .002)]

    return mwave, mflux


def setup_outdir(prefix):
    filesndirs = os.listdir(os.getcwd())
    trk = 1; go = True;
    while go == True:
        name = prefix+'_results_'+str(trk)
        if name not in filesndirs:
            break
        trk += 1

    print('Writing outputs to folder "'+name+'"')
    os.mkdir(name)
    return name
