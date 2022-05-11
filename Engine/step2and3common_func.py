"""
collection of functions that used in both step2 and 3
"""

import numpy as np

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


def trim_obs_data(x, wave, s, u, xbounds, maskwaves=None, step2or3=2):
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

    if step2or3 == 2:
        return s_piece, u_piece, wave_piece, x_piece
    
    elif step2or3 == 3:
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

def check_user_input(args, singleORdouble=1, step2or3=2):

    if step2or3 == 3:
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


def setup_logger(args, outpath, name=None, step2or3=2):

    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        '%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s'
        )

    file_hander = logging.FileHandler(
        f'{outpath}/{args.targname}_{args.band}.log'
        )
    stream_hander = logging.StreamHandler()

    # file_hander.setLevel()
    file_hander.setFormatter(formatter)

    logger.addHandler(file_hander)
    logger.addHandler(stream_hander)
    logger.propagate = False

    if step2or3 == 3:
        logger.info(
            f'Writing output to ./Output/{args.targname}_{args.band}/{name}'
            )

    return logger, stream_hander

def _add_npar(par_in, optgroup, dpars, fitobj):
    """add npar (number if var pars) in fitobj
    Returns:
        class: fitobj
    """

    parmask = np.ones_like(par_in, dtype=bool)
    parmask[:] = False
    for optkind in optgroup:
        parmask[(dpars[optkind] != 0)] = True
    fitobj.npar = len(par_in[parmask])

    return fitobj