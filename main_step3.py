from Engine.importmodule import *
from Engine.importmodule import read_prepdata
from Engine.set_argparse import _argparse_step3

from Engine.IO_AB      import setup_templates, init_fitsread, setup_outdir
from Engine.clips      import basicclip_above
from Engine.contfit    import a0cont
from Engine.classes    import FitObjs,InParams,_setup_bound_cut
from Engine.rebin_jv   import rebin_jv
from Engine.rotint     import rotint
from Engine.opt        import optimizer, fmod
from Engine.outplotter import (outplotter_23, outplotter_rv,
                                    outplotter_rv_combind)
from Engine.detect_peaks import detect_peaks
from Engine.crmask    import cr_masker
from Engine.molmask    import (h2o_masker, mask_wave2pixel_range,
                                    merge_pixel_masks, nonCO_masker)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def setup_fitting_init_pars(band, initvsini, order, initvsini2=0, fluxratio=0):
    """Setup the initial values for the parameters to be optimized (fitted)

    Args:
        band (str): H or K band
        initvsini (float): initial vsini value
        order (int): Current run order
        initvsini2 (float): initial vsini value for the secondary
        fluxratio (float): flux ration between the secondary and the promary

    Returns:
        np.array: initial values for the parameters to be optimized
    """

    # start at bucket loc = 1250 +- 100, width = 250 +- 100,
    # depth = 100 +- 5000 but floor at 0
    centerloc = 1250 if band == 'H' else 1180

    # Initialize parameter array for optimization as well as half-range values
    # for each parameter during the various steps of the optimization.
    # Many of the parameters initialized here will be changed throughout the
    # code before optimization and in between optimization steps.

    pars0 = np.array([
        np.nan,    # 0: The shift of the stellar template (km/s) [assigned later]
        0.3,       # 1: The scale factor for the stellar template
        0.0,       # 2: The shift of the telluric template (km/s)
        0.6,       # 3: The scale factor for the telluric template
        initvsini, # 4: vsini (km/s)
        np.nan,    # 5: The instrumental resolution (FWHM) in pixels
        0.0,       # 6: Wavelength 0-pt
        0.0,       # 7: Wavelength linear component
        0.0,       # 8: Wavelength quadratic component
        0.0,       # 9: Wavelength cubic component
        1.0,       #10: Continuum zero point
        0.0,       #11: Continuum linear component
        0.0,       #12: Continuum quadratic component
        np.nan,    #13: Instrumental resolution linear component
        np.nan,    #14: Instrumental resolution quadratic component
        centerloc, #15: Blaze dip center location
        330,       #16: Blaze dip full width
        0.05,      #17: Blaze dip depth
        90,        #18: Secondary blaze dip full width
        0.05,      #19: Blaze dip depth
        0.0,       #20: Continuum cubic component
        0.0,       #21: Continuum quartic component
        0.0,       #22: Continuum pentic component
        0.0,       #23: Continuum hexic component
        np.nan,    #24: The shift of the second stellar template (km/s) [assigned later]
        0.3,       #25: The scale factor for the second stellar template
        initvsini2,#26: Secondary vsini (km/s)
        fluxratio  #27: Secondary to primary flux ratio S2/S1 (km/s)
    ])

    if int(order) == 13: pars0[1] = 0.8

    return pars0

def _make_dpars(key_name, locs, dpar, numofpars, dpars_org):
    """Conventional func for make new dpar array.
    Not meant to call by the user.

    Args:
        key_name (srt): dict key name
        locs (list): location where dpar needs to change
        dpar (list): dpar values for the locs
        numofpars (int): number of parameters
        dpars_org (dict): dict for the dpars_org

    Returns:
        [dict]: dpars_org
    """

    init_dpars = np.zeros(numofpars)
    init_dpars[locs] = dpar

    dpars_org[key_name] = init_dpars

    return dpars_org


def base_dpars_dict(
    vsini_v1, masterbeam, band, order, numofpars, vsini_v2=-1):
    """Setup basic sets of paramaeter variable ranges

    Args:
        vsini_v1 (float): initial vsini value
        masterbeam (str): A or B beam (nodding) that is running
        band (str): H or K band
        order (int): Current run order
        vsini_v2 (float): initial vsini value for the secondary

    Returns:
        dpars_org (dict): Sets of optimize parameters' variable ranges
    """

    dpars_org = {}
    dpars_org = _make_dpars('cont',
                            [10, 11, 12, 15, 16, 17,  18,  19, 20, 21, 22, 23],
                            [1e7, 1,  1, 10, 30, 0.2, 50, 0.2,  1,  1,  1,  1],
                            numofpars, dpars_org
                            )
    dpars_org = _make_dpars('twave',
                            [3, 6, 7, 8, 9],
                            [1, 1, 1, 1, 1],
                            numofpars, dpars_org
                            )
    dpars_org = _make_dpars('ip',
                            [5],
                            [0.5],
                            numofpars, dpars_org
                            )
    dpars_org = _make_dpars('s',
                            [0, 1],
                            [5, 1],
                            numofpars, dpars_org
                            )
    dpars_org = _make_dpars('v',
                            [4],
                            [vsini_v1],
                            numofpars, dpars_org
                            )
    dpars_org = _make_dpars('ts',
                            [0, 1, 3],
                            [5, 1, 1],
                            numofpars, dpars_org
                            )

    if vsini_v2 != -1:
        dpars_org = _make_dpars('s2',
                            [24, 25],
                            [20, 1],
                            numofpars, dpars_org
                            )
        dpars_org = _make_dpars('v2',
                            [26],
                            [vsini_v2],
                            numofpars, dpars_org
                            )
        dpars_org = _make_dpars('s1s2',
                            [0, 1, 24, 25],
                            [5, 1, 20, 1],
                            numofpars, dpars_org
                            )

    if masterbeam == 'B':
        # setting for the Blaze dip
        dpars_org['cont'][15:20] = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    # blaze fitting order setting
    if band == 'H':
        if order in [13]:
            # fit a quadratic (2) continuum
            dpars_org['cont'][20:] = 0
        elif order in [6,14,21]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        else:
            pass
    elif band == 'K':
        if order in [3,5]:
            # fit a cubic (3) continuum
            dpars_org['cont'][21:] = 0
        elif order in [4,6]:
            # fit a quartic (4) continuum
            dpars_org['cont'][22:] = 0
        else:
            pass

    return dpars_org

def setup_tel(args, night, beam, order, tbdata):
    """Load synthetic telluric template generated during Step 1.
    Also loc molecule line loc that need to be masked

    Args:
        args: information specified by user at command line
        night (str): current run night
        beam (str): A or B beam (nodding) that is running
        order (int): current run order

    Returns:
        watm (array): model telluric wavelength
        satm (array): model telluric flux
        a0contx (array): model telluric x
        continuum (array): model telluric continuum (blaze)
        molnames (list): molecule present in select wave region
        watmmols (array): molecule line wavelength in select wave region
        satmmols (array): molecule line flux in select wave region
    """

    watm = tbdata[f'WATM{order}']
    satm = tbdata[f'SATM{order}']
    a0contx    = tbdata[f'X{order}']
    continuum  = tbdata[f'BLAZE{order}']
    molnames   = tbdata['MOLNAMES']

    # Remove extra rows leftover from having columns of unequal length
    satm = satm[(watm != 0)]
    watm = watm[(watm != 0)]

    # set very low points to zero so that they don't go to NaN when taken
    # to an exponent by template power in fmodel_chi
    satm[(satm < 1e-4)] = 0.
    a0contx = a0contx[(continuum != 0)]
    continuum = continuum[(continuum != 0)]
    molnames  =  molnames[(molnames != '')]

    watmmols = {}; satmmols = {}
    for mol in molnames:
        watm1mol = tbdata['WATM'+mol+str(order)]
        satm1mol = tbdata['SATM'+mol+str(order)]
        satm1mol = satm1mol[(watm1mol != 0)]
        watm1mol = watm1mol[(watm1mol != 0)]
        watmmols[mol] = watm1mol; satmmols[mol] = satm1mol

    return watm, satm, a0contx, continuum, molnames, watmmols, satmmols

def trim_obs_data(x, wave, s, u, xbounds, maskwaves):
    """ Trim obvious outliers above the blaze (i.e. cosmic rays)
    """
    nzones = 5
    x = basicclip_above(x,s,nzones)
    wave = basicclip_above(wave,s,nzones)
    u = basicclip_above(u,s,nzones)
    s = basicclip_above(s,s,nzones)
    x = basicclip_above(x,s,nzones)
    wave = basicclip_above(wave,s,nzones)
    u = basicclip_above(u,s,nzones)
    s = basicclip_above(s,s,nzones)

    # Cut spectrum to within wavelength regions defined in input list
    s_piece    = s[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    u_piece    = u[    (x > xbounds[0]) & (x < xbounds[-1]) ]
    wave_piece = wave[ (x > xbounds[0]) & (x < xbounds[-1]) ]
    x_piece    = x[    (x > xbounds[0]) & (x < xbounds[-1]) ]

    testmask = np.ones_like(s_piece, dtype=bool)
    if len(maskwaves) != 0:
        for maskbounds in maskwaves:
            testmask[(wave_piece*1e4 > maskbounds[0]) \
                        & (wave_piece*1e4 < maskbounds[1]) ] = False

    return s_piece, u_piece, wave_piece, x_piece, testmask

def trim_tel_data(watm, satm, wave_piece, s_piece, u_piece, x_piece):
    """Trim telluric template to data range +- 15 AA. If telluric
    template buffer is cut short because A0 lines didn't extend
    far past data range, cut data range accordingly.
    """

    satm_in = satm[(watm > np.min(wave_piece)*1e4 - 10) \
                        & (watm < np.max(wave_piece)*1e4 + 10)]
    watm_in = watm[(watm > np.min(wave_piece)*1e4 - 10) \
                        & (watm < np.max(wave_piece)*1e4 + 10)]

    s_piece	= s_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                        & (wave_piece*1e4 < np.max(watm_in)-10)]
    u_piece	= u_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                        & (wave_piece*1e4 < np.max(watm_in)-10)]
    x_piece	= x_piece[ (wave_piece*1e4 > np.min(watm_in)+10) \
                        & (wave_piece*1e4 < np.max(watm_in)-10)]
    wave_piece = wave_piece[(wave_piece*1e4 > np.min(watm_in)+10) \
                        & (wave_piece*1e4 < np.max(watm_in)-10)]

    return satm_in, watm_in, wave_piece, s_piece, u_piece, x_piece

def check_if_template_exist(args, singleORdouble=1):

    syntemp = os.listdir(f'./Engine/syn_template')

    if singleORdouble == 1:
        template_kind = args.template.lower()
        temperature = args.temperature
        logg = args.logg
        ll = ''
    elif singleORdouble == 2:
        template_kind = args.template2.lower()
        temperature = args.temperature2
        logg = args.logg2
        ll = '2'
    else:
        sys.exit(f'singleORdouble gives {singleORdouble}, can onlu be 1 or 2.')

    if template_kind == 'synthetic':
        #list of all syntheticstellar
        syntemp = [i for i in syntemp if i[:3] == 'syn']
        synT    = [ i.split('_')[2][1:]  for i in syntemp ]
        synlogg = [ i.split('_')[3][4:7] for i in syntemp ]
    elif template_kind == 'phoenix':
        #list of all phoenix
        syntemp = [i for i in syntemp if i[:3] == 'PHO']
        synT    = [ i.split('-')[1][4:]  for i in syntemp ]
        synlogg = [ i.split('-')[2][:3] for i in syntemp ]
    else:
        synT = [temperature]
        synlogg = [logg]

    if temperature not in synT:
        sys.exit(
            f'ERROR: UNEXPECTED STELLAR TEMPERATURE FOR "-temp{ll}" INPUT! '
            f'{syntemp} AVALIABLE UNDER ./Engine/syn_template/'
            )

    if logg not in synlogg:
        sys.exit(
            f'ERROR: UNEXPECTED STELLAR LOGG FOR "-logg{ll}" INPUT! {syntemp} '
            'AVALIABLE UNDER ./Engine/syn_template/'
            )

def check_user_input(args, singleORdouble=1):

    if args.mode == '':
        sys.exit(
            'ERROR: YOU MUST CHOOSE A MODE, "STD" OR "TAR", for "-mode"'
            )

    if args.initvsini == '':
        sys.exit(
            'ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI VALUE, "-i"'
            )

    if (args.guesses == '') & (args.guessesX == ''):
        sys.exit(
            'ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR RV VALUE(S) BY '
            'USING "-g" OR "-gX"'
            )

    if (args.temperature == '') & (args.logg == ''):
        sys.exit(
            'ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
            'STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE '
            'AVAILABLE TEMPLATES'
            )

    if args.template.lower() not in ['synthetic', 'livingston', 'phoenix']:
        sys.exit('ERROR: UNEXPECTED STELLAR TEMPLATE FOR "-t" INPUT!')


    if singleORdouble == 2:

        if args.fluxratio == '':
            sys.exit('ERROR: YOU MUST PROVIDE A FLUX RATIO S2/S1, "-f"')

        if args.initvsini2 == '':
            sys.exit(
                'ERROR: YOU MUST PROVIDE AN INITIAL GUESS FOR VSINI2 VALUE, "-i2"'
                )

        if (args.temperature2 == '') & (args.logg2 == ''):
            sys.exit(
                'ERROR: YOU MUST PROVIDE THE TEMPERATURE AND LOGG VALUE FOR '
                'SECONDARY STELLAR TEMPLATE. GO TO "./Engine/syn_template/" TO SEE '
                'AVAILABLE TEMPLATES'
                )

        if args.template2.lower() not in ['synthetic', 'livingston', 'phoenix']:
            sys.exit(
                'ERROR: UNEXPECTED SECONDARY STELLAR TEMPLATE FOR "-t" INPUT!'
                )


def setup_init_rv_guess(args):

    if args.mode.lower() == 'std':
        initguesses = float(args.guesses)
        initguesses_show = initguesses

    else: # Load initial RV guesses from file
        if args.guesses_source == 'init': # From Step 2 results
            guesses = './Output/{}_{}/Initguesser_results_{}.csv'.format(
                args.targname,
                args.band,
                int(args.guessesX)
                )
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['night'])
            initrvs    = np.array(guessdata['bestguess'])
            initguesses = {}
            initguesses_show = f'Initguesser_results_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])
            if args.binary:
                initrvs2    = np.array(guessdata['bestguess2'])
                initguesses2 = {}
                for hrt in range(len(initnights)):
                    initguesses2[str(initnights[hrt])] = float(initrvs2[hrt])

                return initguesses, initguesses_show, initguesses2

        elif args.guesses_source == 'rvre': # From Step 3 results
            guesses = './Output/{}_{}/RVresultsSummary_{}.csv'.format(
                args.targname,
                args.band,
                int(args.guessesX)
                )
            guessdata  = Table.read(guesses, format='csv')
            initnights = np.array(guessdata['NIGHT'])
            initrvs    = np.array(guessdata['RVfinal'])
            initguesses = {}
            initguesses_show = f'RVresultsSummary_{args.guessesX}.csv'
            for hrt in range(len(initnights)):
                initguesses[str(initnights[hrt])] = float(initrvs[hrt])
            if args.binary:
                initrvs2    = np.array(guessdata['RV2final'])
                initguesses2 = {}
                for hrt in range(len(initnights)):
                    initguesses2[str(initnights[hrt])] = float(initrvs2[hrt])

                return initguesses, initguesses_show, initguesses2

    return initguesses, initguesses_show

def mkdir_output_dic(args):

    if not os.path.isdir('./Output'):
        os.mkdir('./Output')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}'):
        os.mkdir(f'./Output/{args.targname}_{args.band}')
    filesndirs = os.listdir(f'./Output/{args.targname}_{args.band}')

    trk = 1; go = True
    while go == True:
        name = f'RV_results_{trk}'
        if name not in filesndirs:
            break
        trk += 1

    os.mkdir(f'./Output/{args.targname}_{args.band}/{name}')

    if not os.path.isdir(f'./Output/{args.targname}_{args.band}/figs'):
        os.mkdir(f'./Output/{args.targname}_{args.band}/figs')

    step2or3 = '3'
    temp_dir = f'./Output/{args.targname}_{args.band}/figs/' \
                    f'main_step{step2or3}_{args.band}_{trk}'
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    outpath = f'./Output/{args.targname}_{args.band}'

    return trk, outpath, step2or3, name

def setup_logger(args, outpath, name):

    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')

    file_hander  = logging.FileHandler(f'{outpath}/{args.targname}_{args.band}.log')
    stream_hander= logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    logger.info(
        f'Writing output to ./Output/{args.targname}_{args.band}/{name}')

    return logger, stream_hander

def save_raw_box(args, nights, inparam, name, order,
                    nightsbox, rvbox, parfitbox, vsinibox, tagbox,
                    chibox, rvbox2=None, vsinibox2=None):

    c1 = fits.Column(name=f'NIGHT{order}',  array=nightsbox, format=f'{len(nights[0])}A' )
    c2 = fits.Column(name=f'RV{order}',     array=rvbox,     format='D')
    c3 = fits.Column(name=f'PARFIT{order}', array=parfitbox, format=f'{len(parfitbox[0,:])}D',
                        dim=(1,len(parfitbox[0,:])))
    c4 = fits.Column(name=f'VSINI{order}',  array=vsinibox,  format='D')
    c5 = fits.Column(name=f'TAG{order}',    array=tagbox,    format='4A')
    c6 = fits.Column(name=f'CHI{order}',    array=chibox,    format='D')
    if args.binary:
        c7    = fits.Column(name=f'RV2{order}',     array=rvbox2,     format='D')
        c8    = fits.Column(name=f'VSINI2{order}',  array=vsinibox2,  format='D')
        cols  = fits.ColDefs([c1,c2,c3,c4,c5, c6,c7,c8])
    else:
        cols  = fits.ColDefs([c1,c2,c3,c4,c5,c6])

    hdu_1 = fits.BinTableHDU.from_columns(cols)

    if jerp == 0: # If first time writing fits file, make up filler primary hdu
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul = fits.HDUList([primary_hdu,hdu_1])
        hdul.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
    else:
        hh = fits.open('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name))
        hh.append(hdu_1)
        hh.writeto('{}/{}/RVresultsRawBox.fits'.format(inparam.outpath, name),
            overwrite=True)


def combine_rvs_between_orders(
        n, sigma_ON2, rvmasterbox, vsinibox, jds,
        rvfinal, stdfinal, vsinifinal, orders, Nind,
        offsets, 
        std1=None,nights_use=None):
            
    ind = np.where(
        np.isfinite(sigma_ON2[n,:]) & np.isfinite(rvmasterbox[n,:]))[0]
    weights = (1./sigma_ON2[n,ind]) / (np.nansum(1./sigma_ON2[n,ind])) # normalized
    stdspre = (1./sigma_ON2[n,ind]) #unnormalized weights

    stdbtworders   = np.nanstd(rvmasterbox[n,ind])/np.sqrt(len(rvmasterbox[n,ind]))
    stdofallorders = np.sqrt(np.nansum(sigma_ON2[n,ind]))
    std2 = np.sqrt(stdbtworders**2 - stdofallorders**2)

    rvfinal[n]  = np.nansum( weights*rvmasterbox[n,ind] )
    std0        = 1/np.sqrt(np.nansum(stdspre))
    if offsets:
      stdfinal[n] = np.sqrt(std0**2 + std1**2 + std2**2)
    else:
      stdfinal[n] = std0

    vsinifinal[n] = np.nansum(weights*vsinibox[n,ind])

    # if all the RVs going into the observation's final RV calculation
    # were NaN due to any pevious errors, pass NaN
    if np.nansum(weights) == 0:
        rvfinal[n]    = np.nan
        stdfinal[n]   = np.nan
        vsinifinal[n] = np.nan

    # if more than half of the orders going into the observation's final
    # RV calculation were NaN due to any pevious errors, pass NaN
    if np.sum( np.isnan(rvmasterbox[n,:]) ) > np.floor( len(orders) * 0.5 ):
        rvfinal[n]    = np.nan
        stdfinal[n]   = np.nan
        vsinifinal[n] = np.nan

    if nights_use is None:
        return rvfinal, stdfinal, vsinifinal, None
    else:
        jds_out[n]   = jds[nights_use[n]]
        return rvfinal, stdfinal, vsinifinal, jds_out

def _make_fits_ColDefs(col_name, save_data, save_format):

    return fits.Column(name=col_name, array=save_data, format=save_format)

def _add_npar(par, optgroup, dpars, fitobj):
    """add npar (number if var pars) in fitobj
    Returns:
        class: fitobj
    """

    parmask = np.ones_like(par, dtype=bool)
    parmask[:] = False
    for optkind in optgroup:
        parmask[(dpars[optkind] != 0)] = True
    fitobj.npar = len(par[parmask])

    return fitobj

def main(args, inparam, orders, order_use, trk, step2or3, i):
    """Main function for RV fitting that will be threaded over
    by multiprocessing
    """

    nights = inparam.nights
    night = nights[i] # current looped night

    order = orders[order_use]
    xbounds = inparam.xbounddict[order]

    if args.debug:
        print('Working on order {:02d}/{:02d} ({}), night '
                '{:03d}/{:03d} ({}) PID:{}...'.format(int(order_use)+1,
                                                    len(orders),
                                                    order,
                                                    i+1,
                                                    len(inparam.nights),
                                                    night,
                                                    mp.current_process().pid)
                                                    )

    if args.binary:
        pars0 = setup_fitting_init_pars(
            args.band, inparam.initvsini, order,
            inparam.initvsini2, float(args.fluxratio)
            )
    else:
        pars0 = setup_fitting_init_pars(
            args.band, inparam.initvsini, order
            )

    # Collect relevant beam and filenum info
    tagsnight = []; beamsnight = []
    for tag in inparam.tagsA[night]:
        tagsnight.append(tag)
        beamsnight.append('A')
    for tag in inparam.tagsB[night]:
        tagsnight.append(tag)
        beamsnight.append('B')

    rvsminibox    = np.ones(len(tagsnight))*np.nan
    vsiniminibox  = np.ones(len(tagsnight))*np.nan
    rvsminibox2   = np.ones(len(tagsnight))*np.nan
    vsiniminibox2 = np.ones(len(tagsnight))*np.nan
    tagsminibox   = np.ones(len(tagsnight))*np.nan
    chisminibox   = np.ones(len(tagsnight))*np.nan
    # need to match the dpar numbers
    parfitminibox = np.ones(
        (len(tagsnight), len(pars0))
        )*np.nan

    nightsout = []
    for _ in tagsnight:
        nightsout.append(night)

#-------------------------------------------------------------------------------
    # Collect initial RV guesses
    if type(inparam.initguesses) == dict:
        initguesses = inparam.initguesses[night]
    elif type(inparam.initguesses) == float:
        initguesses = inparam.initguesses
    else:
        sys.exit('ERROR! EXPECTING SINGLE NUMBER OR FILE FOR '
                    'INITGUESSES! QUITTING!')

    if args.binary:
        if type(inparam.initguesses2) == dict:
            initguesses2 = inparam.initguesses2[night]
        elif type(inparam.initguesses2) == float:
            initguesses2 = inparam.initguesses2
        else:
            sys.exit('ERROR! EXPECTING SINGLE NUMBER OR FILE FOR '
                        'INITGUESSES 2! QUITTING!')

    if np.isnan(initguesses) == True:
        logger.warning(f'  --> Previous run of {night} found it '
                                'inadequate (nan), skipping...')
        return (nightsout, rvsminibox, parfitminibox, vsiniminibox,
                    tagsminibox, rvsminibox2, vsiniminibox2, chisminibox)

    if order == 4 and args.band == 'K' and (20180401 < int(night[:8]) < 20190531):
        return (nightsout, rvsminibox, parfitminibox, vsiniminibox,
                    tagsminibox, rvsminibox2, vsiniminibox2, chisminibox)

    # Iterate over all A/B exposures
    for t in range(len(tagsnight)):
        tag = tagsnight[t]
        beam = beamsnight[t]
        masterbeam = beam

        if int(night[:8]) == 20170216 \
            and args.targname == 'GJ281' \
            and float(tag) == 63:
            continue


        # Use instrumental profile dictionary corresponding to whether
        # IGRINS mounting was loose or not
        # [:8] here is to ensure program works under Night_Split mode
        if int(night[:8]) < 20180401 or int(night[:8]) > 20190531:
            IPpars = inparam.ips_tightmount_pars[args.band][masterbeam][order]
        else:
            IPpars = inparam.ips_loosemount_pars[args.band][masterbeam][order]

        if beam == 'A':
            antibeam = 'B'
        elif beam == 'B':
            antibeam = 'A'
        else:
            sys.exit(
                f'EXIT, beam (nodding) can only be A or B, getting {beam}...')

        #-------------------------------------------------------------------------------
        # Load synthetic telluric template generated during Step 1.
        A0loc = f'./Output/{args.targname}_{args.band}/A0Fits/'\
                    f'{night[:8]}A0_{beam}treated_{args.band}.fits'

        hdulist = fits.open(A0loc)
        num_orders = 0
        for i in range(25):
            try:
                hdulist[i].columns[0].name[9:]
                num_orders += 1
            except:
                continue
        fits_layer = [ i for i in np.arange(num_orders)+1 \
                        if int(hdulist[i].columns[0].name[9:]) == order ][0]
        tbdata = hdulist[ fits_layer ].data
        flag = np.array(tbdata[f'ERRORFLAG{order}'])[0]

        # Check whether Telfit hit critical error in Step 1 for the chosen
        # order with this night. If so, skip.
        if flag == 1:
            logger.warning(f'  --> TELFIT ENCOUNTERED CRITICAL ERROR IN ORDER: '
                                f'{order} NIGHT: {night}, skipping...')

            return (nightsout, rvsminibox, parfitminibox, vsiniminibox,
                        tagsminibox, rvsminibox2, vsiniminibox2, chisminibox)

        watm, satm, a0contx, continuum, molnames, watmmols, satmmols = setup_tel(
            args, night, beam, order, tbdata
            )

        maskwaves = h2o_masker(
            inparam, args, order, night, watm, satm, molnames, watmmols, 
            satmmols
            )
        #-------------------------------------------------------------------------------

        bound_cut = _setup_bound_cut(
            inparam.bound_cut_dic, args.band, order
            )

        # Load target spectrum
        x,wave,s,u = init_fitsread(
            f'{inparam.inpath}/{night}/{beam}/',
            'target',
            'separate',
            night,
            order,
            tag,
            args.band,
            bound_cut,
            )

        # Execute S/N cut
        s2n = s/u
        if np.nanmedian(s2n) < float(args.SN_cut):
            logger.warning(
                '  --> Bad S/N {:1.3f} < {} for {}{} {}, SKIP'.format(
                    np.nanmedian(s2n), args.SN_cut, night, beam, tag
                    )
                )
            continue

        s_piece, u_piece, wave_piece, x_piece, testmask = trim_obs_data(
            x, wave, s, u, xbounds, maskwaves
            )

        if (len(s_piece[testmask]) / len(s_piece) < 0.3) \
            or (len(s_piece[testmask]) < 350):
            logger.warning(
                'Only {} unmasked pixels for {} order {}, SKIP'.format(
                    len(s_piece[testmask]), night, order
                    ))
            return (nightsout, rvsminibox, parfitminibox, vsiniminibox,
                        tagsminibox, rvsminibox2, vsiniminibox2, chisminibox)

        # Save data for second template cutting after optimization cycle 1 done
        s_save = s_piece.copy()
        x_save = x_piece.copy()
        u_save = u_piece.copy()

        satm_in, watm_in, wave_piece, s_piece, u_piece, x_piece = trim_tel_data(
            watm, satm, wave_piece, s_piece, u_piece, x_piece
            )

        Rstell1 = np.median(np.diff(inparam.mwave0))

        if args.binary:
            Rstell2 = np.median(np.diff(inparam.mwave2))

            if  Rstell1 > Rstell2:
                rebin2to1 = True; extra1 = 0.; extra2 = 10.
            else:
                rebin2to1 = False; extra1 = 10.; extra2 = 0.

            mflux_in2 = inparam.mflux2[
                    (inparam.mwave2 > np.min(wave_piece)*1e4 - 5 - extra2) \
                    & (inparam.mwave2 < np.max(wave_piece)*1e4 + 5 + extra2)
                ]
            mwave_in2 = inparam.mwave2[
                    (inparam.mwave2 > np.min(wave_piece)*1e4 - 5 - extra2) \
                    & (inparam.mwave2 < np.max(wave_piece)*1e4 + 5 + extra2)
                ]
            Rstell = np.min([Rstell1,Rstell2])

            dstep = Rstell2
            nstep = int((mwave_in2[-1]-mwave_in2[0])/dstep)
            mwave1 = np.linspace(mwave_in2[0],mwave_in2[-1],nstep)
            mflux1 = rebin_jv(mwave_in2,mflux_in2,mwave1,False)
            mwave_in2 = mwave1.copy(); mflux_in2 = mflux1.copy()
            mwave_in2 = mwave_in2[1:-1]
            mflux_in2 = mflux_in2[1:-1]
        else:
            extra1 = 0
            extra2 = 0
            Rstell = Rstell1

        # Trim stellar template to data range +- 10 AA
        mflux_in = inparam.mflux0[
                (inparam.mwave0 > np.min(wave_piece)*1e4 - 5 - extra1) \
                & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5 + extra1)
            ]
        mwave_in = inparam.mwave0[
                (inparam.mwave0 > np.min(wave_piece)*1e4 - 5 - extra1) \
                & (inparam.mwave0 < np.max(wave_piece)*1e4 + 5 + extra1)
            ]
        Rtell = np.median(np.diff(watm_in))
        if Rstell < Rtell:
            sys.exit(f'Telluric template resolution ({np.round(Rtell,4)} AA) '
                        'must be finer than stellar template resolution '
                        '({round(Rstell,4)} AA) !')

        # Rebin stellar template to uniform wavelength scale
        dstep = Rstell1
        nstep = int((mwave_in[-1]-mwave_in[0])/dstep)
        mwave1 = np.linspace(mwave_in[0], mwave_in[-1], nstep)
        mflux1 = rebin_jv(mwave_in, mflux_in, mwave1, False)
        mwave_in = mwave1.copy(); mflux_in = mflux1.copy()
        mwave_in = mwave_in[1:-1]
        mflux_in = mflux_in[1:-1]

        # Normalize continuum from A0 to flux scale of data
        continuum /= np.nanmedian(continuum)
        continuum *= np.nanpercentile(s_piece, 99)

        # --------------------------------------------------------------

        par = pars0.copy()

        # Get initial guess for cubic wavelength solution from reduction pipeline
        f = np.polyfit(x_piece, wave_piece, 3)
        q = np.poly1d(f)
        initwave = q(x_piece)*1e4

        # Initial RV with barycentric correction
        par[0] = initguesses-inparam.bvcs[night+tag]
        par[5]  = IPpars[2]
        par[13] = IPpars[1]
        par[14] = IPpars[0]

        if args.binary:
            par[24] = initguesses2 - inparam.bvcs[night+tag]


        # setup fitting boundary
        if args.binary:
            dpars = base_dpars_dict(inparam.vsinivary, masterbeam,
                                    args.band, int(order), len(pars0),
                                    inparam.vsinivary2
                                    )
        else:
            dpars = base_dpars_dict(inparam.vsinivary, masterbeam,
                                    args.band, int(order), len(pars0)
                                    )

        continuum_in = rebin_jv(a0contx, continuum, x_piece, False)
        fitobj = FitObjs(s_piece, x_piece, u_piece, continuum_in, watm_in,
                            satm_in, mflux_in, mwave_in,
                            ast.literal_eval(inparam.maskdict[order]),
                            masterbeam, [np.array([], dtype=int),
                            np.array([], dtype=int)],
                            initwave, [])

        if args.binary:
            fitobj.addsecondary(mwave_in2, mflux_in2, rebin2to1)

        #-------------------------------------------------------------------------------

        # Initialize an array that puts hard bounds on vsini and the instrumental
        # resolution to make sure they do not diverge to unphysical values
        optimize = True
        par_in = par.copy()
        if masterbeam == 'B':
            hardbounds = [
                par_in[4] - dpars['v'][4],      par_in[4] + dpars['v'][4],
                par_in[5] - dpars['ip'][5],     par_in[5]+dpars['ip'][5]
                ]
        else:
            hardbounds = [
                par_in[4] - dpars['v'][4],      par_in[4] + dpars['v'][4],
                par_in[5]  - dpars['ip'][5],    par_in[5] + dpars['ip'][5],
                par_in[15] - dpars['cont'][15], par_in[15] + dpars['cont'][15],
                par_in[16] - dpars['cont'][16], par_in[16] + dpars['cont'][16],
                0.,                             par_in[17] + dpars['cont'][17],
                par_in[18] - dpars['cont'][18], par_in[18] + dpars['cont'][18],
                0.,                             par_in[19] + dpars['cont'][19]
                ]
        if hardbounds[0] < 0.5:
            hardbounds[0] = 0.5
        if hardbounds[2] < 1:
            hardbounds[2] = 1

        if args.binary:
            hardbounds.append(par_in[26] - dpars['v2'][26])
            hardbounds.append(par_in[26] + dpars['v2'][26])
            if hardbounds[-2] < 0.5:
                hardbounds[-2] = 0.5

        # Begin optimization. Fit the blaze, the wavelength solution, the telluric
        # template power and RV, the stellar template power and RV, the
        # zero point for the instrumental resolution, and the vsini of the star
        # separately, iterating and cycling between each set of parameter fits.

        cycles = 4

        # added new twave step at start
        optgroup1 = ['twave', 'cont', 'twave', 'cont', 'ts',
                    'cont', 'twave', 's', 'cont',
                    'twave',
                    'ip', 'v',
                    'ip', 'v',
                    'twave',  's',
                    'twave',  'ts']

        optgroup2 = ['twave', 'cont', 'twave', 'cont', 'ts', 's1s2',
                    'cont', 'twave', 's', 's2', 'cont',
                    'twave',
                    'ip', 'v', 'v2',
                    'ip', 'v', 'v2',
                    'twave',  's', 's2',
                    'twave',  'ts', 's1s2']

        optgroup = optgroup1.copy()
        initstellpow2 = par_in[25]
        par_in[25] = 0.
        
        nk = 1

        for nc, cycle in enumerate(np.arange(cycles), start=1):
            if cycle == 0:
                parstart = par_in.copy()

            if nc == cycles:
                optgroup = optgroup[:-1] # if is the last run, skip the 'ts'
            
            fitobj = _add_npar(parstart, optgroup, dpars, fitobj)

            for optkind in optgroup:
                parfit_1 = optimizer(parstart, dpars[optkind], hardbounds,
                                        fitobj, optimize, binary=args.binary)
                parstart = parfit_1.copy()
                if args.debug == True:
                    outplotter_23(
                        parfit_1,fitobj,'{}_{}_{}_parfit_{}{}'.format(
                            order,night,tag,nk,optkind),
                            trk, inparam, args, step2or3, order)
                    logger.debug(
                        f'{order}_{tag}_{nk}_{optkind}:\n {parfit_1}')
                nk += 1

            ## After first cycle, use best fit model to identify CRs/hot pixels
            if nc == 1:
                parfit = parfit_1.copy()
                CRmaskF = cr_masker(parfit, fitobj, args.binary)

                _, _, w, _ = fmod(parfit, fitobj, binary=args.binary)
                
                molmask = mask_wave2pixel_range(maskwaves, fitobj, w)
                
                # get st template
                # parfit[3] = 0; parfit[24] = 0
                # smod, _, w, cont = fmod(parfit, fitobj, binary=args.binary)
    
                # nonCO_mask_box = nonCO_masker(
                #     smod, w, cont, int(order), parfit[0], 
                #     fitobj, flux_cut=0.96
                #     )

                # template_mask = merge_pixel_masks(molmask, nonCO_mask_box)
                template_mask = molmask

                fitobj = FitObjs(
                    s_piece, x_piece, u_piece, continuum_in, watm_in,
                    satm_in, mflux_in, mwave_in,
                    ast.literal_eval(inparam.maskdict[order]),
                    masterbeam, CRmaskF, initwave, template_mask)

                if args.binary:
                    optgroup = optgroup2.copy()
                    parstart[25] = initstellpow2
                    fitobj.addsecondary(mwave_in2, mflux_in2, rebin2to1)

                    fitobj = _add_npar(parfit, optgroup, dpars, fitobj)
                
            # if nc > 1:
            #     parfit = parfit_1.copy()
            #     # update the st nonCO mask region based on new rv0
            #     parfit[3] = 0; parfit[24] = 0
            #     smod, _, w, cont = fmod(parfit, fitobj, binary=args.binary)
                
            #     nonCO_mask_box = nonCO_masker(
            #         smod, w, cont, int(order), parfit[0], 
            #         fitobj, flux_cut=0.96
            #         )
                
            #     template_mask = merge_pixel_masks(molmask, nonCO_mask_box)

            #     fitobj.molmask = template_mask


        parfit = parfit_1.copy()

        #-------------------------------------------------------------------------------

        # if best fit stellar template power is very low, throw out result
        if parfit[1] < 0.1:
            logger.warning(f'  --> Stellar template power is low for {night}! '
                                'Data likely being misfit! Throwing out result...')
            continue
        if args.binary and parfit[25] < 0.05:
            logger.warning(f'  --> Secondary stellar template power is low for {night}! '
                                'Data likely being misfit! Throwing out result...')
            continue

        # if best fit stellar or telluric template powers are exactly
        # equal to their starting values, fit failed, throw out result
        if parfit[1] == par_in[1] or parfit[3] == par_in[3] \
            or (args.binary and (parfit[25] == par_in[25])):
            logger.warning(f'  --> Stellar or telluric template powers have '
                                f'not budged from starting values for {night}! '
                                'Fit is broken! Optimizer bounds may be '
                                'unfeasible, or chi-squared may be NaN? '
                                'Throwing out result...')
            continue

        # if best fit model dips below zero at any point, we're to
        # close to edge of blaze, fit may be comrpomised, throw out result
        smod,chisq,trash,trash2 = fmod(parfit, fitobj, args.binary)
        if len(smod[(smod < 0)]) > 0:
            logger.warning(f'  --> Best fit model dips below 0 for {night}! '
                                'May be too close to edge of blaze, throwing '
                                'out result...')
            continue

           
        #-------------------------------------------------------------------------------

        if args.plotfigs == True:
            parfitS1 = parfit.copy(); parfitS1[3] = 0; parfitS1[25] = 0
            parfitS2 = parfit.copy(); parfitS2[3] = 0; parfitS2[1] = 0
            parfitT = parfit.copy(); parfitT[1] = 0; parfitT[25] = 0
            if args.binary:
                outplotter_23(
                    parfitS1, fitobj, 'parfitS1_{}_{}_{}'.format(order,night,tag),
                    trk, inparam, args, step2or3,order)
                outplotter_23(
                    parfitS2, fitobj, 'parfitS2_{}_{}_{}'.format(order,night,tag),
                    trk, inparam, args, step2or3,order)
            else:
                outplotter_23(
                    parfitS1, fitobj, 'parfitS_{}_{}_{}'.format(order,night,tag),
                    trk, inparam, args, step2or3,order)
            outplotter_23(
                parfitT, fitobj, 'parfitT_{}_{}_{}'.format(order,night,tag),
                trk, inparam, args, step2or3,order)
            outplotter_23(
                parfit, fitobj,  'parfit_{}_{}_{}'.format(order,night,tag),
                trk, inparam, args, step2or3,order, chi_new=chisq)

        rv0 = parfit[0]

        # Barycentric correction
        rvsminibox[t]   = rv0  + inparam.bvcs[night+tag] \
                            + rv0*inparam.bvcs[night+tag]/(2.99792458e5**2)
        parfitminibox[t]= parfit
        vsiniminibox[t] = parfit[4]
        tagsminibox[t]  = tag
        fit,chi,trash,trash2 = fmod(parfit, fitobj,args.binary)
        chisminibox[t]  = chi
        if args.binary:
            rv2 = parfit[24]
            rvsminibox2[t]   = rv2  + inparam.bvcs[night+tag] \
                        + rv2*inparam.bvcs[night+tag]/(2.99792458e5**2)
            vsiniminibox2[t] = parfit[26]

    # If any tag has a chisq value 10x greater than the best fit, flag it
    # as a misfit and throw out its result
    for t in np.arange(len(tagsnight)):
        if np.isnan(chisminibox[t]) == False:
            if chisminibox[t] > 10.*np.nanmin(chisminibox):
                logger.warning(
                    f'  --> Chi-squared indicates a misfit for observation '
                            f'{night} {order} {tagsnight[t]}')
                rvsminibox[t]   = np.nan
                vsiniminibox[t] = np.nan

    return (nightsout, rvsminibox, parfitminibox, vsiniminibox,
                tagsminibox, rvsminibox2, vsiniminibox2, chisminibox)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    args = _argparse_step3()
    inpath   = './Input/{}/'.format(args.targname)
    cdbs_loc = '~/cdbs/'

    check_user_input(args, singleORdouble=1)
    check_if_template_exist(args, singleORdouble=1)

    initvsini = float(args.initvsini)
    vsinivary = float(args.vsinivary)

    if args.binary:
        initvsini2 = float(args.initvsini2)
        vsinivary2 = float(args.vsinivary2)

        check_user_input(args, singleORdouble=2)
        check_if_template_exist(args, singleORdouble=2)

    #------------------------------

    if (args.mode.lower() == 'std') & (args.guesses_source != ''):
        sys.exit('ERROR: STD CANNOT USE -gS, PLEASE USE -g')
    if (args.mode.lower() == 'tar') & (args.guesses != ''):
        sys.exit('ERROR: TAR CANNOT USE -g, PLEASE USE -gS')

    #------------------------------

    if args.nAB == '' and args.mode.lower() == 'std':
        nAB = 2
    elif args.nAB == '' and args.mode.lower() == 'tar':
        nAB = 3
    else:
        nAB = int(args.nAB)

    #------------------------------

    if args.binary:
        initguesses, initguesses_show, initguesses2 = setup_init_rv_guess(args)
    else:
        initguesses, initguesses_show = setup_init_rv_guess(args)

    #-------------------------------------------------------------------------------

    start_time = datetime.now()
    print('\n')
    print('####################################################################################\n')
    print('---------------------------------------------------------------')
    print(u'''
Input Parameters:
    Target             =  {}
    Filter              = \33[37;1;41m {} band \033[0m
    WaveLength file     = \33[37;1;41m WaveRegions_{} \033[0m
    S/N cut             > \33[37;1;41m {} \033[0m
    Minium # of AB sets = \33[37;1;41m {} \033[0m  <------- If TAR mode, this should be at least 3. If STD mode, at least 2.
    Initial vsini       = \33[37;1;41m {} km/s \033[0m
    vsini vary range    \u00B1 \33[37;1;41m {} km/s \033[0m
    RV initial guess    = \33[37;1;41m {} \033[0m
    Stellar template use= \33[37;1;41m {} \033[0m
    syn template temp   = \33[37;1;41m {} \033[0m
    syn template logg   = \33[37;1;41m {} \033[0m
    syn template B      = \33[37;1;41m {} \033[0m
    Threads use         = {}
    '''.format(args.targname, args.band, args.WRegion, args.SN_cut, nAB,
               initvsini, vsinivary, initguesses_show, args.template,
               args.temperature, args.logg, args.B, args.Nthreads))
    if args.binary:
        print(u'''
PLUS BINARY PARAMETERS:
    Initial vsini #2       = \33[37;1;41m {} km/s \033[0m
    vsini #2 vary range    \u00B1 \33[37;1;41m {} km/s \033[0m
    Stellar template #2 use= \33[37;1;41m {} \033[0m
    syn template temp #2   = \33[37;1;41m {} \033[0m
    syn template logg #2   = \33[37;1;41m {} \033[0m
    syn template B #2      = \33[37;1;41m {} \033[0m
       '''.format(initvsini2, vsinivary2, args.template2,
                  args.temperature2, args.logg2, args.B2))
    if not args.skip:
        while True:
            inpp = input("Press [Y]es to continue, [N]o to quite...\n --> ")
            if 'n' in inpp.lower():
                sys.exit('QUIT, PLEASE RE-ENTER YOUR PARAMETERS')
            elif 'y' in inpp.lower():
                break
            else:
                print('I cannot understand what you are saying... TRY AGAIN')
                continue

    print('---------------------------------------------------------------')
    print('Running Step 3 for {}...'.format(args.targname))
    print('This will take a while..........')

    trk, outpath, step2or3, name = mkdir_output_dic(args)
    logger, stream_hander = setup_logger(args, outpath, name)

    # Read in the Prepdata under ./Input/Prpedata/
    xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs = read_prepdata(
        args
        )

    # Use subset of nights if specified
    if args.nights_use != '':
        nightstemp = np.array(ast.literal_eval(args.nights_use), dtype=str)
        for nnn in nightstemp:
            if nnn not in nightsFinal:
                sys.exit(
                    'NIGHT {} NOT FOUND UNDER ./Input_Data/{}'.format(
                        nnn, args.targname
                        ))
        nightsFinal = nightstemp
        print('Only processing nights: {}'.format(nightsFinal))

    logger.info('Analyze with {} nights'.format(len(nightsFinal)))

    #-------------------------------------------------------------------------------

    # Retrieve stellar and telluric templates
    print('\n Loading stellar template... \n')
    watm,satm, mwave0, mflux0 = setup_templates(
        logger, args.template, args.band, int(args.temperature),
        float(args.logg), float(args.B)
        )

    # Save pars in class for future use
    inparam = InParams(
        inpath, outpath, initvsini, vsinivary, args.plotfigs, initguesses,
        bvcs, tagsA, tagsB, nightsFinal, mwave0, mflux0, None, xbounddict,
        maskdict
        )

    if args.binary:
        print('\n Loading secondary stellar template... \n')
        watm,satm, mwave2, mflux2 = setup_templates(
            logger, args.template2, args.band, int(args.temperature2),
            float(args.logg2), float(args.B2)
            )
        inparam.addsecondary(
            initvsini2, vsinivary2, mwave2, mflux2, initguesses2
            )

    #-------------------------------------------------------------------------------

    # Divide between nights where IGRINS mounting was loose (L) and when
    # it was tight (T). All statistical analysis will be performed separately
    # for these two datasets.
    nights    = inparam.nights
    intnights = np.array([int(i[:8]) for i in nights])

    if len(intnights[(intnights >= 20180401) & (intnights < 20190531)]) > 0:
        logger.info('''
WARNING: Some of these nights were when the IGRINS K band was defocused!
For K band RVs: IGRINS RV will take this into account and process these nights
                slightly differently. When you run Step 3, RVs will be output
                in two formats: one with the defocus nights separated, and the
                other with all nights together.
For H band RVs: We do not expect any systematic changes in the H band as the
                result of the defocus. IGRINS RV will process defocus nights
                the same way as the others, but when you run Step 3, will still
                output the results in two formats like it does with the K band.
''')

    indT = np.where((intnights < 20180401) | (intnights > 20190531))
    indL = np.where((intnights >= 20180401) & (intnights < 20190531))

    nightsT = nights[indT]
    nightsL = nights[indL]
    obsT    = obs[indT]
    obsL    = obs[indL]
    rvmasterboxT  = np.ones((len(nightsT),len(orders)))*np.nan
    stdmasterboxT = np.ones((len(nightsT),len(orders)))*np.nan
    rvmasterboxL  = np.ones((len(nightsL),len(orders)))*np.nan
    stdmasterboxL = np.ones((len(nightsL),len(orders)))*np.nan
    vsinisT       = np.ones((len(nightsT),len(orders)))*np.nan
    vsinisL       = np.ones((len(nightsL),len(orders)))*np.nan
    rvmasterboxT2  = np.ones((len(nightsT),len(orders)))*np.nan
    stdmasterboxT2 = np.ones((len(nightsT),len(orders)))*np.nan
    rvmasterboxL2  = np.ones((len(nightsL),len(orders)))*np.nan
    stdmasterboxL2 = np.ones((len(nightsL),len(orders)))*np.nan
    vsinisT2       = np.ones((len(nightsT),len(orders)))*np.nan
    vsinisL2       = np.ones((len(nightsL),len(orders)))*np.nan

    if len(nightsL) > 0 and len(nightsT) > 0:
        nightscomblist = [nightsT,nightsL]
        T_Ls = ['T','L']
    elif len(nightsL) > 0:
        nightscomblist = [nightsL]
        T_Ls = ['L']
    else:
        nightscomblist = [nightsT]
        T_Ls = ['T']

    #-------------------------------------------------------------------------------
    # if not in debug mode than enter quite mode, i.e., all message saved in log file
    if not args.debug: logger.removeHandler(stream_hander)
    print('\n')

    # Run order by order, multiprocessing over nights within an order
    for jerp in range(len(orders)):
        if not args.debug:
            print('Working on order {} ({:02d}/{:02d})'.format(
                orders[jerp], int(jerp+1), len(orders)
                ))

        func = partial(main, args, inparam, orders, jerp, trk, step2or3 )
        outs = pqdm(np.arange(len(nightsFinal)), func, n_jobs=args.Nthreads)

        # for ii in np.arange(len(nightsFinal)):
        #     main(args, inparam, orders, jerp, trk, step2or3, ii)

        #-------------------------------------------------------------------------------

        # Collect outputs: the reference night, the best fit RV, vsini, and other parameters
        for i in range(len(nights)):
            outsbox = outs[i]
            if i == 0:
                nightsbox = outsbox[0]
                rvbox     = outsbox[1]
                parfitbox = outsbox[2]
                vsinibox  = outsbox[3]
                tagbox    = outsbox[4]
                rvbox2    = outsbox[5]
                vsinibox2 = outsbox[6]
                chibox    = outsbox[7]
            else:
                nightsbox = nightsbox + outsbox[0]
                rvbox     = np.concatenate((rvbox,outsbox[1]))
                parfitbox = np.vstack((parfitbox,outsbox[2]))
                vsinibox  = np.concatenate((vsinibox,outsbox[3]))
                tagbox    = np.concatenate((tagbox,outsbox[4]))
                rvbox2    = np.concatenate((rvbox2,outsbox[5]))
                vsinibox2 = np.concatenate((vsinibox2,outsbox[6]))
                chibox    = np.concatenate((chibox,outsbox[7]))

        order = orders[jerp]
        nightsbox = np.array(nightsbox)
        vsinitags = []

        if args.binary:
            save_raw_box(
                args, nights, inparam, name, order,
                nightsbox, rvbox, parfitbox, vsinibox, tagbox,
                chibox, rvbox2=rvbox2, vsinibox2=vsinibox2
                )
        else:
            save_raw_box(
                args, nights, inparam, name, order,
                nightsbox, rvbox, parfitbox, vsinibox, tagbox, chibox
                )

        #-------------------------------------------------------------------------------
        # For each set of nights (tight, loose)...
        iT_L = 0
        for nights_use in nightscomblist:

            T_L = T_Ls[iT_L]

            # For each night...
            for i in range(len(nights_use)):
                # Collect the RVs and vsinis determined from different A/B
                # exposures within a night
                indnight  = np.where(nightsbox == nights_use[i])[0]
                rvtags    = rvbox[indnight]
                vsinitags = vsinibox[indnight]

                # Take the mean of the vsinis, and the mean and std of the RVs.
                # If the number of different successfully fitted A/B exposures
                # is less than required, pass NaN instead.
                if T_L == 'T':
                    vsinisT[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxT[i,jerp]  = np.nan
                        stdmasterboxT[i,jerp] = np.nan
                    else:
                        rvmasterboxT[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxT[i,jerp] = np.nanstd(rvtags)\
                            /np.sqrt(len(rvtags[~np.isnan(rvtags)]))

                else:
                    # Don't use vsini estimates from this order during loose epoch
                    if order == 3 and args.band == 'K':
                        vsinisL[i,jerp] = np.nan
                    else:
                        vsinisL[i,jerp] = np.nanmean(vsinitags)

                    if (np.sum(~np.isnan(rvtags)) < nAB ):
                        rvmasterboxL[i,jerp]  = np.nan
                        stdmasterboxL[i,jerp] = np.nan
                    else:
                        rvmasterboxL[i,jerp]  = np.nanmean(rvtags)
                        stdmasterboxL[i,jerp] = np.nanstd(rvtags)\
                            /np.sqrt(len(rvtags[~np.isnan(rvtags)]))

                if args.binary:
                    rvtags    = rvbox2[indnight]
                    vsinitags = vsinibox2[indnight]
                    if T_L == 'T':
                        vsinisT2[i,jerp] = np.nanmean(vsinitags)
                        if (np.sum(~np.isnan(rvtags)) < nAB ):
                            rvmasterboxT2[i,jerp]  = np.nan
                            stdmasterboxT2[i,jerp] = np.nan
                        else:
                            rvmasterboxT2[i,jerp]  = np.nanmean(rvtags)
                            stdmasterboxT2[i,jerp] = np.nanstd(rvtags)\
                                /np.sqrt(len(rvtags[~np.isnan(rvtags)]))
                    else:
                        if order == 3 and args.band == 'K':
                            vsinisL2[i,jerp] = np.nan
                        else:
                            vsinisL2[i,jerp] = np.nanmean(vsinitags)

                        if (np.sum(~np.isnan(rvtags)) < nAB ):
                            rvmasterboxL2[i,jerp]  = np.nan
                            stdmasterboxL2[i,jerp] = np.nan
                        else:
                            rvmasterboxL2[i,jerp]  = np.nanmean(rvtags)
                            stdmasterboxL2[i,jerp] = np.nanstd(rvtags)\
                                /np.sqrt(len(rvtags[~np.isnan(rvtags)]))

            iT_L += 1


    #-------------------------------------------------------------------------------
    if not args.debug: logger.addHandler(stream_hander)
    print('\n')

    # Don't combine Loose and Tight datasets, but make them both easily referenceable
    nightsCombined  = np.array([])
    jdsCombined = np.array([])
    rvfinalCombined = np.array([])
    stdfinalCombined = np.array([])
    vsinifinalCombined = np.array([])
    rvfinalCombined2 = np.array([])
    stdfinalCombined2 = np.array([])
    vsinifinalCombined2 = np.array([])

    if T_Ls == ['T','L']:
        rvboxcomblist  = [rvmasterboxT,rvmasterboxL]
        stdboxcomblist = [stdmasterboxT,stdmasterboxL]
        vsinicomblist  = [vsinisT,vsinisL]
        rvboxcomblist2  = [rvmasterboxT2,rvmasterboxL2]
        stdboxcomblist2 = [stdmasterboxT2,stdmasterboxL2]
        vsinicomblist2  = [vsinisT2,vsinisL2]
        obscomblist    = [obsT,obsL]
    elif T_Ls == ['L']:
        rvboxcomblist  = [rvmasterboxL]
        stdboxcomblist = [stdmasterboxL]
        vsinicomblist  = [vsinisL]
        rvboxcomblist2  = [rvmasterboxL2]
        stdboxcomblist2 = [stdmasterboxL2]
        vsinicomblist2  = [vsinisL2]
        obscomblist    = [obsL]
    else:
        rvboxcomblist  = [rvmasterboxT]
        stdboxcomblist = [stdmasterboxT]
        vsinicomblist  = [vsinisT]
        rvboxcomblist2  = [rvmasterboxT2]
        stdboxcomblist2 = [stdmasterboxT2]
        vsinicomblist2  = [vsinisT2]
        obscomblist    = [obsT]

    # Iterate over tight and loose mounting data sets...
    for boxind in range(len(rvboxcomblist)):

        rvmasterbox  = rvboxcomblist[boxind]
        stdmasterbox = stdboxcomblist[boxind]
        vsinibox     = vsinicomblist[boxind]
        rvmasterbox2  = rvboxcomblist2[boxind]
        stdmasterbox2 = stdboxcomblist2[boxind]
        vsinibox2     = vsinicomblist2[boxind]
        obsbox       = obscomblist[boxind]

        # For every order, check what portion of observations did not
        # return RVs. If that portion is less than half those from the
        # most successful order, then it was masked too much to be
        # useable for this target and it would be inconsistent to
        # consider it for such a small fraction of the data --
        # so we throw it out.
        goodcounts = np.array(
            [len(np.where(np.isfinite(rvmasterbox[:,ll]))[0]) for ll in range(len(orders))]
            )

        for ll in range(len(orders)):
            if goodcounts[ll] < 0.5*np.max(goodcounts):
                rvmasterbox[:,ll] = np.nan
                rvmasterbox2[:,ll] = np.nan

        #-------------------------------------------------------------------------------

        if args.mode=='STD': # If RV STD star, calculate uncertainty in method

            # Calculate the precision within an order across nights
            sigma_O2     = np.array(
                [np.nanstd(rvmasterbox[:,ll])**2 for ll in range(len(orders))]
                )
            sigma_ABbar2 = np.ones_like(sigma_O2)

            # Calculate uncertainty in method as difference between variance
            # within an order and mean variance within a night's different exposure RVs
            for ll in range(len(orders)):
                sigma_ABbar2[ll] = np.nanmedian(stdmasterbox[:,ll]**2)
            sigma_method2 = sigma_O2 - sigma_ABbar2

        # If target star, load the uncertainty in method calculated from our
        # RV STD star runs
        else:
            if T_Ls[boxind] == 'T':
                nights_use = nightsT.copy()
                kind = 'Focused'
                sigma_method2 = inparam.methodvariance_tight[args.band]
            elif T_Ls[boxind] == 'L':
                nights_use = nightsL.copy()
                kind = 'Defocus'
                sigma_method2 = inparam.methodvariance_loose[args.band]

        sigma_ON2      = np.ones_like(rvmasterbox)
        sigma_ON2bi    = np.ones_like(rvmasterbox)

        #-------------------------------------------------------------------------------

        # Note rvmasterbox indexed as [nights,orders]
        Nnights = len(rvmasterbox[:,0])

        for ll in range(len(orders)):
            # Calculate the uncertainty in each night/order RV as the sum of the
            # uncertainty in method and the uncertainty in that night's As and Bs RVs
            for night in range(Nnights):
                sigma_ON2[night,ll]   = sigma_method2[ll] + stdmasterbox[night,ll]**2
                sigma_ON2bi[night,ll] = sigma_method2[ll] + stdmasterbox2[night,ll]**2

        rvfinal    = np.ones(Nnights, dtype=float)
        stdfinal   = np.ones(Nnights, dtype=float)
        vsinifinal = np.ones(Nnights, dtype=float)
        rvfinal2    = np.ones(Nnights, dtype=float)
        stdfinal2   = np.ones(Nnights, dtype=float)
        vsinifinal2 = np.ones(Nnights, dtype=float)
        jds_out   = np.ones(Nnights, dtype=float)

        if T_Ls[boxind] == 'T':
            nights_use = nightsT.copy()
            kind = 'Focused'
        elif T_Ls[boxind] == 'L':
            nights_use = nightsL.copy()
            kind = 'Defocused'

        Nord = len(orders)
        ordermeans = np.ones(Nord)*np.nan
        orderstds  = np.ones(Nord)*np.nan

        for jerp in range(Nord):
            ind = np.where(
                np.isfinite(sigma_ON2[:,jerp]) & np.isfinite(rvmasterbox[:,jerp]))[0]
            weights = (1./sigma_ON2[ind,jerp]) / (np.nansum(1./sigma_ON2[ind,jerp])) # normalized
            stdspre = (1./sigma_ON2[ind,jerp]) #unnormalized weights
            ordermeans[jerp] = np.nansum( weights*rvmasterbox[ind,jerp] )
            orderstds[jerp]  = 1/np.sqrt(np.nansum(stdspre))

        std1 = np.sqrt(np.nansum([orderstds[jerp]**2 for jerp in range(Nord)]))

        offsets = False
        for jerp1 in range(Nord):
            for jerp2 in range(Nord):
                if jerp1 != jerp2 and np.isfinite(ordermeans[jerp1]) and np.isfinite(ordermeans[jerp2]):
                    if abs(ordermeans[jerp1]-ordermeans[jerp2]) > np.sqrt(orderstds[jerp1]**2 + orderstds[jerp2]**2):
                        offsets = True

        masterorder = 6
        masterjerp  = np.arange(Nord)[orders == masterorder]
        if offsets == True:
            for jerp in range(Nord):
                rvmasterbox[:,jerp] = rvmasterbox[:,jerp] + ordermeans[masterjerp] - ordermeans[jerp]

        if args.binary and offsets:
            ordermeans2 = np.ones(Nord)*np.nan
            orderstds2  = np.ones(Nord)*np.nan
            
            for jerp in range(Nord):
                ind = np.where(
                    np.isfinite(sigma_ON2bi[:,jerp]) & np.isfinite(rvmasterbox2[:,jerp]))[0]
                weights = (1./sigma_ON2bi[ind,jerp]) / (np.nansum(1./sigma_ON2bi[ind,jerp])) # normalized
                stdspre = (1./sigma_ON2bi[ind,jerp]) #unnormalized weights
                ordermeans2[jerp] = np.nansum( weights*rvmasterbox2[ind,jerp] )
                orderstds2[jerp]  = 1/np.sqrt(np.nansum(stdspre))

            std1_2 = np.sqrt(np.nansum([orderstds2[jerp]**2 for jerp in range(Nord)]))

            for jerp in range(Nord):
                rvmasterbox2[:,jerp] = rvmasterbox2[:,jerp] + ordermeans2[masterjerp] - ordermeans2[jerp]

            
        # Combine RVs between orders using weights calculated from uncertainties
        for n in range(Nnights):

            Nind = np.where(intnights == int(nightsFinal[n][:8]))[0]

            rvfinal, stdfinal, vsinifinal, jds_out = combine_rvs_between_orders(
                n, sigma_ON2, rvmasterbox, vsinibox, jds, rvfinal,
                stdfinal, vsinifinal, orders, Nind, offsets, std1=std1, nights_use
                )

            if args.binary:
                rvfinal2, stdfinal2, vsinifinal2, _ = combine_rvs_between_orders(
                    n, sigma_ON2bi, rvmasterbox2, vsinibox2, jds, rvfinal2,
                    stdfinal2, vsinifinal2, orders, Nind, offsets, std1=std1_2
                    )

        #-------------------------------------------------------------------------------

        if args.binary:
            rvps = [rvfinal,rvfinal2]
            stdps = [stdfinal,stdfinal2]
        else:
            rvps = [rvfinal]
            stdps = [stdfinal]

        for ara in range(len(rvps)):
            rvp = rvps[ara]
            stdp = stdps[ara]
            # Plot results

            outplotter_rv(rvp, stdp, args, inparam, name, ara, kind)

        #-------------------------------------------------------------------------------
        # Save results to fits file separately for each tight/loose dataset
        name_list = [
            'NIGHT', 'JD', 'RVBOX', 'STDBOX', 'Sigma_method2', 'Sigma_ON2',
            'RVfinal', 'STDfinal', 'VSINI'
            ]
        data_list = [
            nights_use, jds_out, rvmasterbox, stdmasterbox, sigma_method2,
            sigma_ON2, rvfinal, stdfinal, vsinifinal
            ]
        format_list = [
            '8A', 'D', f'{len(orders)}D', f'{len(orders)}D', 'D',
            f'{len(orders)}D', 'D', 'D', 'D'
            ]

        if args.mode=='STD':
            name_list.extend(['Sigma_O2', 'Sigma_ABbar2'])
            data_list.extend([sigma_O2, sigma_ABbar2])
            format_list.extend(['D', 'D'])

        elif args.binary:
            name_list.extend(
                ['RVBOX2', 'STDBOX2', 'Sigma_ON2bi', 'RV2final', 'STD2final',
                'VSINI2']
                )
            data_list.extend(
                [rvmasterbox2, stdmasterbox2, sigma_ON2bi, rvfinal2,
                stdfinal2, vsinifinal2]
                )
            format_list.extend(
                [f'{len(orders)}D', f'{len(orders)}D', f'{len(orders)}D',
                'D', 'D', 'D']
                )

        cc = []
        for nn, dd, ff in zip(name_list, data_list, format_list):
            cc.append(
                _make_fits_ColDefs(nn, dd, ff)
                )
        cols  = fits.ColDefs(cc)

        hdu_1 = fits.BinTableHDU.from_columns(cols)
        bleh = np.ones((3,3))
        primary_hdu = fits.PrimaryHDU(bleh)
        hdul        = fits.HDUList([primary_hdu,hdu_1])
        if len(T_Ls) == 2:
            hdul.writeto(
                '{}/{}/RVresultsAdvanced_{}.fits'.format(
                    inparam.outpath, name, kind),
                overwrite=True)
        else:
            hdul.writeto(
                '{}/{}/RVresultsAdvanced.fits'.format(
                    inparam.outpath, name),
                overwrite=True)

        #-------------------------------------------------------------------------------
        # Combine final RVs from both tight and loose mounting data sets
        nightsCombined     = np.concatenate((nightsCombined, nights_use))
        jdsCombined        = np.concatenate((jdsCombined, jds_out))
        rvfinalCombined    = np.concatenate((rvfinalCombined, rvfinal))
        stdfinalCombined   = np.concatenate((stdfinalCombined, stdfinal))
        vsinifinalCombined = np.concatenate((vsinifinalCombined, vsinifinal))
        rvfinalCombined2   = np.concatenate((rvfinalCombined2, rvfinal2))
        stdfinalCombined2  = np.concatenate((stdfinalCombined2, stdfinal2))
        vsinifinalCombined2= np.concatenate((vsinifinalCombined2,vsinifinal2))

        if args.mode=='STD': # If uncertainty in method was calculated, save it
            sigma_method2 = [np.around(float(i), 8) for i in sigma_method2]
            logger.info(
                'sigma_method2 during the {} epoch is {}'.format(
                    kind, sigma_method2)
                    )
        if len(T_Ls) == 2:
            logger.info(
                'During the {} epoch: RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                kind,
                np.nanmean(rvfinal),
                np.nanstd(rvfinal) )
                )
            logger.info(
                'During the {} epoch: vsini mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                kind,
                np.nanmean(vsinifinal),
                np.nanstd(vsinifinal) )
                )
        else:
            logger.info(
                'RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                np.nanmean(rvfinal),
                np.nanstd(rvfinal) )
                )
            logger.info(
                'vsini mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                np.nanmean(vsinifinal),
                np.nanstd(vsinifinal) )
                )
            if args.binary:
                logger.info(
                    'Secondary RV mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                    np.nanmean(rvfinal2),
                    np.nanstd(rvfinal2) )
                    )
                logger.info(
                    'Secondary vsini mean = {:1.4f} km/s, std = {:1.4f} km/s'.format(
                    np.nanmean(vsinifinal2),
                    np.nanstd(vsinifinal2) )
                    )

    #-------------------------------------------------------------------------------

    # Plot combined results
    if args.binary:
        rvps = [rvfinalCombined, rvfinalCombined2]
        stdps = [stdfinalCombined, stdfinalCombined2]
    else:
        rvps = [rvfinalCombined]
        stdps = [stdfinalCombined]

    for ara in range(len(rvps)):
        rvp = rvps[ara]
        stdp = stdps[ara]

        xscale = np.arange(len(rvp))+1

        outplotter_rv_combind(
            xscale, rvp, stdp, args, inparam, name, ara, kind, nightsT, nightsL)

    #-------------------------------------------------------------------------------
    # Output combined final results to fits file
    name_list = ['NIGHT', 'JD', 'RVfinal', 'STDfinal', 'VSINI']
    data_list = [
        nightsCombined, jdsCombined, rvfinalCombined, stdfinalCombined,
        vsinifinalCombined
        ]
    format_list = [f'{len(nights[0])}A', 'D', 'D', 'D', 'D']

    if args.binary:
        name_list.extend(['RV2final', 'STD2final', 'VSINI2'])
        data_list.extend(
            [rvfinalCombined2, stdfinalCombined2, vsinifinalCombined2]
            )
        format_list.extend(['D', 'D', 'D'])

    cc = []
    for nn, dd, ff in zip(name_list, data_list, format_list):
        cc.append(
            _make_fits_ColDefs(nn, dd, ff)
            )
    cols  = fits.ColDefs(cc)

    hdu_1 = fits.BinTableHDU.from_columns(cols)

    bleh = np.ones((3,3))
    primary_hdu = fits.PrimaryHDU(bleh)
    hdul = fits.HDUList([primary_hdu,hdu_1])
    hdul.writeto(
        '{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name),
        overwrite=True)

    tempin = Table.read(
        '{}/{}/RVresultsSummary.fits'.format(inparam.outpath, name),
        format='fits')
    tempin.write(
        '{}/RVresultsSummary_{}.csv'.format(inparam.outpath, trk),
        format='csv', overwrite=True)
#-------------------------------------------------------------------------------

    if len(T_Ls) == 2:
        logger.info(
            'Combined RV results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(
            np.nanmean(rvfinalCombined),
            np.nanstd(rvfinalCombined))
            )
        logger.info(
            'vsini results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(
            np.nanmean(vsinifinalCombined),
            np.nanstd(vsinifinalCombined))
            )

        if args.binary:
            logger.info(
                'Combined RV 2 results: mean={:1.4f} km/s, std={:1.4f} km/s'.format(
                np.nanmean(rvfinalCombined2),
                np.nanstd(rvfinalCombined2))
                )
            logger.info(
                'vsini 2 results:       mean={:1.4f} km/s, std={:1.4f} km/s'.format(
                np.nanmean(vsinifinalCombined2),
                np.nanstd(vsinifinalCombined2))
                )

    warning_r = log_warning_id(
        f'{outpath}/{args.targname}_{args.band}.log', start_time
        )
    if warning_r:
        print(f'''
**********************************************************************************
WARNING!! you got warning message during this run. Please check the log file under:
          {outpath}/{args.targname}_{args.band}_A0Fits.log
**********************************************************************************
''')
    end_time = datetime.now()
    logger.info('Whole process DONE!!!!!!, Duration: {}'.format(end_time - start_time))
    logger.info('Output saved under {}/{}'.format(inparam.outpath, name) )
    logger.info('The final RV estimates you are looking for are in the RVresultsSummary files!')
    print('####################################################################################')
