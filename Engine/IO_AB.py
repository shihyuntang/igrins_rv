import os
from os import listdir
from os.path import isfile, join#, isdir
import re
import sys

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits

from Engine.rebin_jv import rebin_jv

def partial_loader(inpath0, order):
    '''
    Access data arrays from reduced IGRINS fits files.

    Inputs:
    inpath0 : Path to file, including filename up to ".spec.fits" suffix
    order   : Echelle order, as characterized by file index (as opposed to m 
                number; for conversion between the two, see Stahl et al. 2021)

    Outputs:
    wavelist : Wavelength solution from plp reduction (known to be imprecise)
    fluxlist : Corresponding flux
    s2nlist  : Corresponding signal to noise ratio
    '''

    if os.path.exists(inpath0+'.spec.fits'):
        hdulist = fits.open(inpath0+'.spec.fits')
    elif os.path.exists(inpath0+'.spec.fits.gz'):
        hdulist = fits.open(inpath0+'.spec.fits.gz')
    else:
        sys.exit(f'ERROR: data not found: {inpath0}.spec.fits')
    tbdata_w = hdulist[1].data
    tbdata_f = hdulist[0].data

    if os.path.exists(inpath0+'.sn.fits'):
        hdulist_sn = fits.open(inpath0+'.sn.fits')
    elif os.path.exists(inpath0+'.sn.fits.gz'):
        hdulist_sn = fits.open(inpath0+'.sn.fits.gz')
    else:
        sys.exit(f'ERROR: data not found: {inpath0}.sn.fits')

    tbdata_sn = hdulist_sn[0].data

    wavelist = np.array(tbdata_w[order],dtype=np.float64)
    fluxlist = np.array(tbdata_f[order],dtype=np.float64)
    s2nlist = np.array(tbdata_sn[order],dtype=np.float64)

    return wavelist,fluxlist,s2nlist


def init_fitsread(path,kind,beam,night,order,tag,band,Ncuts=None):
    '''
    Fetches reduced spectrum from file, performs some basic quality checks and 
    cleaning, and returns the results in a convenient format.

    Inputs:
    path  : Path to fits file up to some subdirectory (folder organization 
            varies depending on whether loading target or telluric standard, 
            separate A and Bs or combined
    kind  : Specifies whether loading a telluric standard or target
    beam  : Frame type (A, B, or combinedAB)
    night : Date of observation in YYYYMMDD
    order : Echelle order, as characterized by file index (as opposed to m 
            number; for conversion between the two, see Stahl et al. 2021)
    tag   : Number of observation, represented by last four digits in fits filename
    band  : H or K band
    Ncuts : List of format [L,R], where L specifies number of pixels to trim 
            off left (low wavelength) end of spectrum, and R vice versa

    Outputs:
    x     : Absolute pixel value (if spectra has edges trimmed, will begin at 
            a nonzero number)
    wave  : Wavelength solution from plp reduction (known to be imprecise)
    s     : Corresponding flux
    u     :  Corresponding uncertainty
    '''

    if kind not in ['A0', 'target']:
        sys.exit('kind MUST BE "target" OR "A0", FORCE QUITTING!')

    # Initguesser takes combined AB versions of target spectra, of which 
    # there may be multiple, and combines /those/ together.
    if beam[:8] == 'combined':

        if kind == 'target':
            subpath = '{}{}/{}/'.format(path, night,beam[8:])
        else:
            subpath = '{}std/{}/{}/'.format(path, night,beam[8:])

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
            wavelist0, fluxlist0, s2nlist0 = partial_loader(
                fullpathprefix+tag,order
                )

            if Ncuts != None:
                Nstartcut = Ncuts[0]
                Nendcut = Ncuts[-1]

                ia = np.ones_like(wavelist0, dtype=bool)
                ia[0:Nstartcut] = False
                ia[-Nendcut:len(ia)] = False
                
                wavelist1 = wavelist0[ia]
                fluxlist1 = fluxlist0[ia]
                s2nlist1 = s2nlist0[ia]
                xlist = np.arange(len(wavelist1), dtype=float) \
                            + np.float(Nstartcut)
            else:
                wavelist1 = wavelist0
                fluxlist1 = fluxlist0
                s2nlist1 = s2nlist0
                xlist = np.arange(len(wavelist1), dtype=float)

            try:
                fluxstack = np.vstack((fluxstack, fluxlist1))
                wavestack = np.vstack((wavestack, wavelist1))
                s2nstack  = np.vstack((s2nstack, s2nlist1))
            except UnboundLocalError:
                fluxstack = fluxlist1.copy()
                wavestack = wavelist1.copy()
                s2nstack = s2nlist1.copy()

        try:
            wavelist = np.array(
                [np.nanmean(wavestack[:,i]) for i in range(len(wavelist1))]
                )
            fluxlist = np.array(
                [np.nanmean(fluxstack[:,i]) for i in range(len(wavelist1))]
                )
            s2nlist = np.array(
                [np.sqrt(
                    np.nansum(s2nstack[:,i]**2)) for i in range(len(wavelist1))
                    ]
                )
            #ulist    = np.array([1/np.sqrt(np.nansum(1/(ustack[:,i]**2))) for i in range(len(wavelist1))])
        except UnboundLocalError:
            return 0,0,0,0

    # For A0 loading, use combined AB and load from AB/ folder. Shouldn't be 
    # more than one AB file, so shouldn't need to do any extra combining, 
    # unlike above.
    # For RVProc procedure, use As and Bs separately and load from 
    # A/ or B/ folders.
    else: # if separate beams
        if kind == 'target':
            pass
        else:
            path = '{}std/{}/AB/'.format(path, night)

        wavelist0, fluxlist0, s2nlist0 = partial_loader(
            '{}SDC{}_{}_{}'.format(path, band, night[:8], tag),
            order)

        if Ncuts != None:
            Nstartcut = Ncuts[0]
            Nendcut = Ncuts[-1]

            ia = np.ones_like(wavelist0, dtype=bool)
            ia[0:Nstartcut] = False
            ia[-Nendcut:len(ia)] = False

            wavelist = wavelist0[ia]
            fluxlist = fluxlist0[ia]
            s2nlist = s2nlist0[ia]
            xlist = np.arange(len(wavelist), dtype=float) + np.float(Nstartcut)
        else:
            wavelist = wavelist0
            fluxlist = fluxlist0
            s2nlist = s2nlist0
            xlist = np.arange(len(wavelist), dtype=float)

    MAD = np.nanmedian(
            np.abs(np.nanmedian(s2nlist)-s2nlist)
            )
    x3    = xlist[ (np.isnan(fluxlist) == False) \
                    & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    wave3 = wavelist[(np.isnan(fluxlist) == False) \
                    & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    s2n3  = s2nlist[ (np.isnan(fluxlist) == False) \
                    & (np.isnan(s2nlist) == False) & (fluxlist > 0)]
    s3    = fluxlist[(np.isnan(fluxlist) == False) \
                    & (np.isnan(s2nlist) == False) & (fluxlist > 0)]

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
    u[(s2n < 10)] = 1e3*np.max(u)
    s[(s < 1e-6)] = 1e-6 # Provide absolute floor to flux values of spectra

    return x, wave, s, u

def air2vac(wave):
    """
    This code is from the PyAstronomy repository and is referenced as such 
    (Czela et al. 2019, https://pyastronomy.readthedocs.io/en/latest/index.html)
    Converts air wavelengths to vaccuum wavelengths returns a float or array 
    of wavelengths

    Inputs:
      wave - Wavelengths in air, in Angstroms, float or array

    Outputs:
      newwave - Wavelengths in a vacuum, in Angstroms, float or array

    NOTES:
      1. The procedure uses the IAU standard conversion formula
         from Morton (1991 Ap. J. Suppl. 77, 119)
    """

    wave = np.array(wave, float)
    # Converting to Wavenumber squared
    sigma2 = (1.0e4/wave)**2
    # Computing conversion factor
    fact = 1.0 + 6.4328e-5 + 2.94981e-2 / (146.0-sigma2) \
            + 2.5540e-4 / (41.0-sigma2)
    fact = fact*(wave >= 2000.0) + 1.0*(wave < 2000.0)
    newwave = wave*fact
    return newwave

def _pdread2astrotable(csvgzdir):
    """As pandas can read .csv.gz direactally, use pd read then trans to 
    astropy table.
    """
    df = pd.read_csv(csvgzdir)
    tb = Table.from_pandas(df)
    return tb

def setup_templates(logger, kind='synthetic', band='K', 
        temperature=5000, logg=4.5, B=0):
    '''
    Fetches static stellar and/or telluric templates from file.

    Inputs:
    logger      : Mechanism to keep updating log.txt
    kind        : Specifies kind of stellar template ('synthetic' = IGRINS RV 
                    team generated models, 'phoenix' = for IGRINS RV team only, 
                    'user' = user provided template)
    band        : H or K band
    temperature : Effective temperature corresponding to stellar template
    logg        : log(g) corresponding to stellar template

    Outputs:
    watm    : Wavelength scale of static telluric template
    satm    : Corresponding flux of static telluric template
    mwave0  : Wavelength scale of stellar template
    mflux0  :  Corresponding flux of stellar template
    '''

    if (kind == 'synthetic'):
        logger.info(f'Using {band}-band synthetic stellar template...')
        logger.info(
            f'synthetic stellar template with T{temperature} logg{logg}!!!!!')

        if str(B) == '0':
            temploc = f'syntheticstellar_{band.lower()}band_T{temperature}_logg{logg}_0.0kG.csv.gz'
        else:
            temploc = f'syntheticstellar_{band.lower()}band_T{temperature}_logg{logg}_{B}kG.csv.gz'
        if 'igrins' in os.getcwd().split('/')[-1]:
            if os.path.exists(f'./Engine/syn_template/{temploc}'):
                stelldata = _pdread2astrotable(f'./Engine/syn_template/{temploc}')
            else:
                stelldata = _pdread2astrotable(f'./Engine/syn_template/{temploc[:-1]}')
        else:
            if os.path.exists(f'../Engine/syn_template/{temploc}'):
                stelldata = _pdread2astrotable(f'../Engine/syn_template/{temploc}')
            else:
                stelldata = _pdread2astrotable(f'../Engine/syn_template/{temploc[:-1]}')

        mwave0 = np.array(stelldata['wave'])
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0

    elif (kind.lower() == 'phoenix'):
        logger.info(f'Using {band}-band PHOENIX stellar template...')
        logger.info(f'PHOENIX stellar template with T{temperature} logg{logg}!!!!!')

        temploc = f'PHOENIX-lte0{temperature}-{logg}0-0.0_contadj.csv.gz'
        if 'igrins' in os.getcwd().split('/')[-1]:
            if os.path.exists(f'./Engine/syn_template/{temploc}'):
                stelldata = _pdread2astrotable(f'./Engine/syn_template/{temploc}')
            else:
                stelldata = _pdread2astrotable(f'./Engine/syn_template/{temploc[:-1]}')
        else:
            if os.path.exists(f'../Engine/syn_template/{temploc}'):
                stelldata = _pdread2astrotable(f'../Engine/syn_template/{temploc}')
            else:
                stelldata = _pdread2astrotable(f'../Engine/syn_template/{temploc[:-1]}')
        
        mwave0 = np.array(stelldata['wave'])
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0

    elif (kind.lower() == 'user'):
        logger.info(f'Using user-provided stellar template with T{temperature} '
                        'logg{logg}...')
        logger.warning(f'WARNING! PRECISION AND ACCURACY OF THIS PIPELINE IS '
                            'NOT GUARANTEED WITH USER-PROVIDED STELLAR '
                            'TEMPLATES!!!!!')
        logger.warning(f'Make sure to characterize your errors appropriately '
                            '- see IGRINS RV paper for details.')
        logger.warning(f'Also be sure your templates follow the naming and '
                            'formatting conventions described in the github '
                            'wiki AND are placed under ./Engine/user_templates.')

        if 'igrins' in os.getcwd().split('/')[-1]:
            stelldata = Table.read(
                f'./Engine/user_templates/user_T{temperature}_logg{logg}_{band}band.txt',
                format='ascii')
        else:
            stelldata = Table.read(
                f'../Engine/user_templates/user_T{temperature}_logg{logg}_{band}band.txt',
                format='ascii')

        mwave0 = np.array(stelldata['wave'])
        mflux0 = np.array(stelldata['flux'])
        mwave0 = mwave0[(np.isfinite(mflux0))]
        mflux0 = mflux0[(np.isfinite(mflux0))]
        mflux0[(mflux0 < 0)] = 0

    else:
        logger.info(f'Input kind is {kind}, but must be either "synthetic", '
                        '"phoenix" (for IGRINS RV team usage only), or "user"!')

    if 'igrins' in os.getcwd().split('/')[-1]:
        telluricdata = Table.read(
            './Engine/PhotoAtl_Organized.csv', format='csv')  
    else:
        telluricdata = Table.read(
            '../Engine/PhotoAtl_Organized.csv', format='csv') 


    watm = np.array(telluricdata['wave'])*10000.0
    satm = np.array(telluricdata['flux'])
    watm = watm[(np.isfinite(satm))]
    satm = satm[(np.isfinite(satm))]
    satm[(satm < 0)] = 0
    
    # Make sure wavelength is ascending always
    mflux0 = np.array([x for _, x in sorted(zip(mwave0, mflux0))])
    mwave0 = np.array(list(sorted(mwave0)))

    # Remove duplicate wavelength values from stellar template 
    # (doesn't affect steps 1-3, but needed for bisectors)
    ind = []
    maxwave = mwave0[0]
    for i in range(1,len(mwave0)-1):
        if mwave0[i] > maxwave:
            maxwave = mwave0[i]
        else:
            ind.append(i)
    ind = np.array(ind)
    mask = np.ones(len(mwave0), dtype=bool)
    if len(ind) > 0:
        mask[ind] = False
        mwave0 = mwave0[mask]
        mflux0 = mflux0[mask]
        
    dstep0 = np.median(np.diff(mwave0))
    if dstep0 > 0.045:
        logger.info(f'Stellar template resolution is ~{round(dstep,4)} '
                        'Angstrom, leaving alone...')
    else:
        dstep = 0.045
        nstep = int((mwave0[-1]-mwave0[0])/dstep)
        mwave1 = np.linspace(mwave0[0],mwave0[-1],nstep)
        mflux1 = rebin_jv(mwave0,mflux0,mwave1,False)
        mwave0 = mwave1.copy(); mflux0 = mflux1.copy()
        mwave0 = mwave0[1:-1]
        mflux0 = mflux0[1:-1]

        logger.info(f'Stellar template resolution is ~{round(dstep0,4)} '
                        'Angstrom, rebinning to 0.045 Angstrom...')
        
    return watm, satm, mwave0, mflux0



def setup_templates_tel():
    '''
    Stripped version of setup_templates() for Step 1, where only telluric 
    template needed, so the stellar template is chosen by default.

    Outputs:
    watm    : Wavelength scale of static telluric template
    satm    : Corresponding flux of static telluric template
    mwave0  : Wavelength scale of stellar template
    mflux0  : Corresponding flux of stellar template
    '''

    if 'igrins' in os.getcwd().split('/')[-1]:
        spotdata = Table.read(
            './Engine/SpotAtl_Organized.csv', format='csv')
    else:
        spotdata = Table.read(
            '../Engine/SpotAtl_Organized.csv', format='csv')

    mwave0 = np.array(spotdata['wave'])*10000.0
    mflux0 = np.array(spotdata['flux'])
    mwave0 = mwave0[(np.isfinite(mflux0))]
    mflux0 = mflux0[(np.isfinite(mflux0))]
    mflux0[(mflux0 < 0)] = 0

    if 'igrins' in os.getcwd().split('/')[-1]:
        telluricdata = Table.read(
            './Engine/PhotoAtl_Organized.csv', format='csv')  
    else:
        telluricdata = Table.read(
            '../Engine/PhotoAtl_Organized.csv', format='csv') 

    watm = np.array(telluricdata['wave'])*10000.0
    satm = np.array(telluricdata['flux'])
    watm = watm[(np.isfinite(satm))]
    satm = satm[(np.isfinite(satm))]
    satm[(satm < 0)] = 0
    return watm, satm, mwave0, mflux0



def setup_outdir(prefix):
    '''
    Checks what number of times the code has been run on a given target/band 
    and makes a new output directory accordingly.

    Inputs:
    prefix : Directory of overall output for target/band

    Outputs:
    name : Directory of output for this run of code

    '''
    filesndirs = os.listdir(os.getcwd())
    trk = 1; go = True
    while go == True:
        name = prefix+'_results_'+str(trk)
        if name not in filesndirs:
            break
        trk += 1

    print('Writing outputs to folder "'+name+'"')
    os.mkdir(name)
    return name
