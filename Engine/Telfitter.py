from Engine.importmodule import *
import os, nlopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from   astropy import units
from   Engine.rebin_jv import rebin_jv
from   telfit import TelluricFitter, DataStructures

#-------------------------------------------------------------------------------
def round_decimals_down(number, decimals=2):
    """
    Returns a value rounded down to a specific number of decimal places.
    Borrowed from https://kodify.net/python/math/round-decimals/
    """
    factor = 10 ** decimals
    return np.floor(number * factor) / factor


def wavefunc(par,grad):
    '''
    Takes Telfitted template and uses an input wavelength solution to rebin it for direct comparison with Livingston.

    Inputs:
    par  : array of polynomial coefficients specifying wavelength solution
    grad : Always "None" (has to be this way for NLOpt)

    Outputs reduced chisq of model fit.
    '''

    global watm_Liv, satm_Liv, satmLivGen, x;
    #Make the wavelength scale
    f = np.poly1d(par)
    w = f(x)

    if (w[-1] < w[0]) or (w[-1] > watm_Liv[-1]+5):
        return 1e3

    satmTel2 = rebin_jv(w,satmLivGen,watm_Liv,False)
    return np.sum((satm_Liv - satmTel2)**2) / (len(satmTel2) - len(par))


def wavefit(par0, dpar0):
    '''
    NLopt convenience function for fitting wavelength solution of Telfitted Livingston Atlas such that it matches with Livingston Atlas.

    Inputs:
    par0  : Initial guesses for polynomial coefficients specifying wavelength solution
    dpar0 : Amount each initial guess can vary (higher or lower)

    Outputs:
    parfit : Best fit polynomial coefficients specifying wavelength solution
    '''

    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 7)
    opt.set_min_objective(wavefunc)
    lows  = par0-dpar0
    highs = par0+dpar0
    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)
    # opt.set_maxtime(1200) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_ftol_rel(1e-14)
    parfit = opt.optimize(par0)

    # parfit_floor = round_decimals_down(parfit, decimals=14)
    return parfit

#------------
# to suppress print out from Telfit
@suppress_stdout
def suppress_Fit(fitter, data):
    model = fitter.Fit(data=data, resolution_fit_mode="SVD", adjust_wave="model", air_wave=False)
    return model

@suppress_stdout
def suppress_GenerateModel(fitter, parfit, args):
    model = fitter.GenerateModel(parfit, nofit=True, air_wave=False)
    return model

#------------


def telfitter(watm_in, satm_in, a0ucut, inparam, night, order, args, masterbeam, logger):
    '''
    Produce synthetic telluric template from fit to telluric standard observation. How and why it works is detailed in comments throughout the code.

    Inputs:
    watm_in    : Wavelength scale of telluric standard spectrum
    satm_in    : Corresponding flux of telluric standard spectrum
    a0ucut     : Corresponding uncertainty of telluric standard spectrum
    inparam    : Class containing variety of information (e.g. on observing conditions)
    night      : Date of observation in YYYYMMDD
    order      : Echelle order, as characterized by file index (as opposed to m number; for conversion between the two, see Stahl et al. 2021)
    args       : Information as input by user from command line
    masterbeam : A or B frame

    Outputs:
    wavefitted : Wavelength scale of synethetic telluric spectrum
    satmTel    : Corresponding flux of synthetic telluric spectrum
    names      : Descriptors of Telfit parameters
    parfitted  : Values of best-fit Telfit parameters
    wcont1     : Wavelength scale corresponding to best-fit continuum (from intermediate Telfit step)
    cont1      : Flux corresponding to best-fit continuum (from intermediate Telfit step)
    '''

    os.environ['PYSYN_CDBS'] = inparam.cdbsloc
    fitter = TelluricFitter(debug=False, print_lblrtm_output=args.debug)

    #Set the observatory location with a keyword
    DCT_props     = {"latitude": 34.744, "altitude": 2.36} #altitude in km
    McD_props     = {"latitude": 30.710, "altitude": 2.07}
    GaminiS_props = {"latitude": -30.241, "altitude": 2.72}

    if inparam.obses[night] == 'DCT':
        fitter.SetObservatory(DCT_props)
    elif inparam.obses[night] == 'McD':
        fitter.SetObservatory(McD_props)
    elif inparam.obses[night] == 'GeminiS':
        fitter.SetObservatory(GaminiS_props)
    else:
        sys.exit('TELFIT OBSERVATORY ERROR, OLNY SUPPORT DCT, McD & GeminiS IN THIS VERSION!')

    # Read in data
    watm_in = watm_in/10 # AA --> nm
    data = DataStructures.xypoint(x=watm_in, y=satm_in, cont=None, err=a0ucut) # input wavelength in nm

    # DCT data has parameters describing night of observation that the McDonald data does not.
    if inparam.temps[night] != 'NOINFO': # If such information is available:

        angle       = np.float(inparam.zds[night])           #Zenith distance
        pressure    = np.float(inparam.press[night])         #Pressure, in hPa
        humidity    = np.float(inparam.humids[night])        #Percent humidity, at the observatory altitude
        temperature = np.float(inparam.temps[night])+273.15  #Temperature in Kelvin

        if (order <= 4):
            resolution  = 55000.0                             #Resolution lambda/delta-lambda
        else:
            resolution  = 45000.0                             #Resolution lambda/delta-lambda

        # Ideally, we'd fit resolution as well since that varies across the detector.
        # But in practice the Telfit's resolution fits often diverge to unphysically high values.
        # Ultimately, we only want accurate estimates for the chemical compositions, which are unaffacted
        # by fixing the resolution at 45000. The final telluric template we'll be synthesizing from this
        # will be at a set high resolution, anyway, and when we need/have needed to estimate resolution
        # across the detector, we (have) done so via the fits of the telluric and stellar templates to the
        # observed A0 and GJ281 spectra.

        # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
        if (3 < order < 9) & (args.band =='K'):
            num_fit = 4
            # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co": 5e-3, "n2o":5e-2})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "pressure": pressure,\
                                "temperature": temperature,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co2": 3.675e2,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-5,1e2],\
                              "co": [ 1e-6,1e2]})
        elif (order >= 9 or order <= 3) & (args.band =='K'):
            num_fit = 4
            # Only molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, N2O, and CO2.
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co2": 3.675e2, "n2o":5e-2})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "pressure": pressure,\
                                "temperature": temperature,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co": 5e-3,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-5,1e2],\
                              "co2": [1.0, 1e4]})
        elif args.band =='H':
            num_fit = 3
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co": 5e-3,"co2": 3.675e2,"n2o" : 5e-2})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "pressure": pressure,\
                                "temperature": temperature,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-6,1e2],\
                              "co": [1e-6,1e2],\
                              "co2": [1.0, 1e4]})

    elif inparam.zds[night] != 'NOINFO': # If GeminiS data, some but not all parameters are in fits file.
          # If parameters are not in fits file, use initial guesses and letting them vary.
          # Guesses are taken from mean of parameters from DCT GJ281 data.

        angle       = np.float(inparam.zds[night])           #Zenith distance
        humidity    = np.float(inparam.humids[night])        #Percent humidity, at the observatory altitude

        if (order <= 4):
            resolution  = 55000.0                             #Resolution lambda/delta-lambda
        else:
            resolution  = 45000.0                             #Resolution lambda/delta-lambda

        # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
        if (3 < order < 9) & (args.band =='K'):
            num_fit = 4
            # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co": 5e-3, "n2o":5e-2, "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co2": 3.675e2,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-5,1e2],\
                              "temperature": [265.,300.],\
                              "pressure": [1010.,1035.],\
                              "co": [ 1e-6,1e2]})

        elif (order >= 9 or order <= 3) & (args.band =='K'):
            num_fit = 4
            # Only molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, N2O, and CO2.
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co2": 3.675e2, "n2o":5e-2, "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co": 5e-3,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-5,1e2],\
                              "temperature": [265.,300.],\
                              "pressure": [1010.,1035.],\
                              "co2": [1.0, 1e4]})
        elif args.band =='H':
            num_fit = 3
            fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co": 5e-3,"co2": 3.675e2,"n2o" : 5e-2, "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"angle": angle,\
                                "resolution": resolution,
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,  10.0],\
                              "n2o": [1e-6,1e2],\
                              "co": [1e-6,1e2],\
                              "temperature": [265.,300.],\
                              "pressure": [1010.,1035.],\
                              "co2": [1.0, 1e4]})

    else: # If parameters are not in fits file, use initial guesses and letting them vary.
          # Guesses are taken from mean of parameters from DCT GJ281 data.

        if (order <= 4):
            resolution  = 55000.0                             #Resolution lambda/delta-lambda
        else:
            resolution  = 45000.0                             #Resolution lambda/delta-lambda

        # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
        if (3 < order < 9) & (args.band =='K'):
            num_fit = 6
            fitter.FitVariable({"h2o": 43.,"ch4": 1.8,"co": 5e-3, "n2o":5e-2,
                                "angle": 39., "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"resolution": resolution,\
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co2": 3.675e2,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,10.0],\
                              "n2o": [1e-5,1e2],\
                              "temperature": [265.,300.],\
                              "angle": [1.,75.],\
                              "pressure": [1010.,1035.],\
                              "co": [ 1e-6,1e2]})

        elif (order >= 9 or order <= 3) & (args.band =='K'):
            num_fit = 7
            # Only molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, N2O, and CO2.
            fitter.FitVariable({"h2o": 43.,"ch4": 1.8,"co2": 3.675e2, "n2o": 5e-2,
                                "angle": 39., "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"resolution": resolution,\
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "co": 5e-3,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,10.0],\
                              "temperature": [265.,300.],\
                              "angle": [1.,75.],\
                              "n2o":[1e-5,1e2],\
                              "pressure": [1010.,1035.],\
                              "co2": [1.0, 1e4]})

        elif args.band =='H':
            num_fit = 6
            fitter.FitVariable({"h2o": 43.,"ch4": 1.8,"co": 5e-3,"co2": 3.675e2,"n2o" : 5e-2,
                                "angle": 39., "pressure":1023., "temperature":280.87})

            #Adjust parameters that will not be fit, but are important
            fitter.AdjustValue({"resolution": resolution,\
                                "wavestart": data.x[0]-0.001,\
                                "waveend": data.x[-1]+0.001,\
                                "o3": 7.6e-4,\
                                "o2": 2.1e5,\
                                "no": 0.,\
                                "so2": 5e-9,\
                                "no2": 5e-9,\
                                "nh3": 5e-9,\
                                "hno3": 1.56e-7})

            #Set bounds on the variables being fit
            fitter.SetBounds({"h2o": [1.0, 99.0],\
                              "ch4": [.1,10.0],\
                              "n2o": [1e-6,1e2],\
                              "co": [1e-6,1e2],\
                              "temperature": [265.,300.],\
                              "angle": [1.,75.],\
                              "pressure": [1010.,1035.],\
                              "co2": [ 1,1e4]})

    try:
        if args.debug:
            model = fitter.Fit(data=data, resolution_fit_mode="SVD", adjust_wave="model",air_wave=False)
        else:
            model = suppress_Fit(fitter, data)
    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan],[np.nan],[np.nan]

    '''
      resolution_fit_mode = SVD ought to give faster, more accurate fits for the deep telluric lines we mostly see in K band
      air_wave = False because data in vacuum wavelengths
      adjust_wave =
                    From Telfit comments: "Can be set to either 'data' or 'model'. To wavelength calibrate the
                    data to the telluric lines, set to 'data'. If you think the wavelength
                    calibration is good on the data (such as Th-Ar lines in the optical),
                    then set to 'model' Note that currently, the vacuum --> air conversion
                    for the telluric model is done in a very approximate sense, so
                    adjusting the data wavelengths may introduce a small (few km/s) offset
                    from what it should be. That is fine for relative RVs, but probably not
                    for absolute RVs."

                    As it turns out, the model internal to Telfit is still not very precise in wavelength space, since it relies on the HITRAN
                    database. Some lines are accurate to 1 m/s, but some less so.

                    Hence, our procedure is as follows:

                    1) Fit the A0 spectrum using the Livingston telluric template. Get out a precisely calibrated wavelength solution for the spectrum.
                        (This is done in A0Fitter, not Telfitter)
                    2) Use that wavelength solution was input for Telfit, and let the wavelength scale of the model vary with respect to it.
                    3) Using Telfit's best fit parameters, generate a telluric template at high resolution.
                       Note: Telfit's default when generating a model is to employ a vacuum/air conversion. In order to avoid that,
                             I have manually edited one of Telfit's files.
                    4) To properly calibrate this template in wavelength space, we fit it to a Telfit'd version of the Livingston telluric template, only allowing
                       the wavelength solution to vary.
    '''

    #Get the improved continuum from the fitter
    cont1  = fitter.data.cont
    wcont1 = model.x*10 # nm-->AA

    # chi_new = np.sum((satm_in - model.y*cont1)**2. / model.u**2.)
    # chi_new = chisq / (len(model.y) - num_fit)

    if args.plotfigs:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=300)

        axes.plot(10*watm_in, satm_in,       color='black',    alpha=.8, label='data',      lw=0.7)
        axes.plot(10*model.x, model.y*cont1, color='tab:red',  alpha=.8, label='model fit', lw=0.7)
        axes.plot(10*model.x, cont1,         color='tab:blue', alpha=.8, label='blaze fit', lw=0.7)

        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(2))
        axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
        axes.set_ylabel(r'Flux',       size=6, style='normal' , family='sans-serif' )
        axes.set_xlabel(r'Wavelength [$\AA$]', size=6, style='normal' , family='sans-serif' )
        axes.legend(fontsize=5, edgecolor='white')
        axes.set_title('A0Telfit_Order{}_{}_{}.png'.format(order, night, masterbeam),
                         size=6, style='normal', family='sans-serif')
        # fig.text(0.65, 0.2, r'$\rm \chi^{{2}}_{{\nu}}$ = {:1.2f}'.format(chi_new),
        #                     size=6, style='normal', family='sans-serif')
        fig.savefig('{}/figs_{}/A0Telfit_Order{}_{}_{}.png'.format(inparam.outpath, args.band, order, night, masterbeam),
                    format='png', bbox_inches='tight', overwrite=True)

    ############### Generate template with these parameters but at higher resolution

    names = ["pressure", "temperature", "angle", "resolution",'wavestart','waveend',
                         "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "no",
                         "so2", "no2", "nh3", "hno3"]

    parfitted = np.ones_like(names, dtype=float)
    for k in range(len(names)):
        parfitted[k] = np.float(fitter.GetValue(names[k]) )

    fitter2 = TelluricFitter(debug=False, print_lblrtm_output=args.debug)

    if inparam.obses[night] == 'DCT':
        fitter2.SetObservatory(DCT_props)
    elif inparam.obses[night] == 'McD':
        fitter2.SetObservatory(McD_props)
    elif inparam.obses[night] == 'GeminiS':
        fitter2.SetObservatory(GaminiS_props)

    # Compute telluric template with highest resolution of Livingston template.
    # Add extra space at ends to make sure template covers wider range than data.
    Livingston_minimum_wsep = .035/10
    IGRINS_minimum_wsep     = .130 # <-- This would compute template with IGRINS resolution, sensibly coarser than Livingston

    newwave = np.arange(np.min(watm_in)-2.5, np.max(watm_in)+2.5, Livingston_minimum_wsep) #in nm

    data2 = DataStructures.xypoint(x=newwave,
                                   y=None,
                                   cont=None,
                                   err=None)
    params = {}
    for k in range(len(names)):
        params[names[k]] = np.float(parfitted[k])

    params['wavestart'] = data2.x[0] -0.001
    params['waveend']   = data2.x[-1]+0.001

    fitter2.AdjustValue(params)
    fitter2.ImportData(data2)

    # Call the modeller. On rare occasions, this returns an error. I have no idea what is causing this error, as the
    # FORTRAN readout is quite unhelpful and anyone else who apepars to have experienced this problem had it randomly go away at some point.
    # If this happens, simply deliver NAN arrays, and in later parts of the RV analysis A0 fits from the nearest compatible observation will be used.
    try:
        if args.debug:
            model2 = fitter2.GenerateModel(parfitted, nofit=True, air_wave=False)
        else:
            model2 = suppress_GenerateModel(fitter2, parfitted, args)

    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan],[np.nan],[np.nan]

    watm_save = watm_in.copy(); satm_save = satm_in.copy();
    newwave1 = newwave[(newwave > watm_in[0]-1.0) & (newwave < watm_in[-1]+1.0)]

    # Parameters for reproducing Livingston template with Telfit
    if args.band == 'K':
        telparsdict = { '1':np.array([1.01469894e+03, 2.86278974e+02, 2.57160505e+01, 6.00000000e+05,
               2.42400187e+03, 2.50999640e+03, 1.44631214e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 7.06751837e+00, 2.21814724e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '2':np.array([1.01000565e+03, 2.86022221e+02, 1.02518785e+01, 6.00000000e+05,
               2.39100388e+03, 2.47599903e+03, 1.60984192e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 9.51814513e-02, 2.80273169e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '3': np.array([1.01208206e+03, 2.83372443e+02, 1.30994741e+01, 6.00000000e+05,
               2.36000083e+03, 2.44399739e+03, 2.03888166e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 1.78694998e-01, 2.63207260e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '4': np.array([1.01044914e+03, 3.09376975e+02, 5.00457250e+01, 6.00000000e+05,
               2.32700464e+03, 2.41199998e+03, 2.23826025e+00, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 1.44749523e-01, 1.78797014e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '5': np.array([1.01028238e+03, 2.80000941e+02, 5.02034212e+01, 6.00000000e+05,
               2.29700026e+03, 2.38099471e+03, 2.00001482e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 5.04544799e-03, 1.79978608e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '6': np.array([1.01003245e+03, 2.82986081e+02, 4.69226227e+01, 6.00000000e+05,
               2.26700430e+03, 2.35099941e+03, 1.83825870e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 1.33618090e-01, 1.75621792e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '7': np.array([1.01676398e+03, 2.80004945e+02, 5.00497774e+01, 6.00000000e+05,
               2.23800235e+03, 2.32199758e+03, 1.99995159e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 5.00089318e-03, 1.79986683e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '8': np.array([1.02000000e+03, 2.80000000e+02, 5.00000000e+01, 6.00000000e+05,
               2.21000474e+03, 2.29299818e+03, 2.00000000e+01, 3.67500000e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '9': np.array([1.01010028e+03, 2.83521347e+02, 4.94584722e+01, 6.00000000e+05,
               2.20299994e+03, 2.24499864e+03, 1.62628503e+01, 2.09026148e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.84204544e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '10': np.array([1.02057631e+03, 2.80065296e+02, 4.97608855e+01, 6.00000000e+05,
               2.17600055e+03, 2.21800093e+03, 1.99862859e+01, 3.75841930e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.77986802e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '11': np.array([1.01175621e+03, 2.80258981e+02, 5.01757000e+01, 6.00000000e+05,
               2.15000090e+03, 2.19200022e+03, 1.94041560e+01, 3.17393616e+02,
               7.60000000e-04, 1.54427423e-01, 5.00000000e-03, 1.69487644e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '12': np.array([1.01939759e+03, 3.03543029e+02, 5.16166925e+01, 6.00000000e+05,
               2.12399934e+03, 2.16599710e+03, 4.88527175e+00, 3.62082287e+02,
               7.60000000e-04, 4.24872185e-01, 5.00000000e-03, 1.76259774e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '13': np.array([1.01289532e+03, 2.97875731e+02, 3.68795863e+01, 6.00000000e+05,
               2.09899925e+03, 2.14099888e+03, 8.74545147e+00, 4.70227887e+02,
               7.60000000e-04, 4.02655042e-01, 5.00000000e-03, 4.78062748e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '14': np.array([1.01088512e+03, 2.82651527e+02, 4.97893186e+01, 6.00000000e+05,
               2.07500079e+03, 2.11599814e+03, 1.70464824e+01, 3.75172445e+02,
               7.60000000e-04, 3.21242442e-01, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '15': np.array([1.01015623e+03, 3.09948799e+02, 3.32690854e+01, 6.00000000e+05,
               2.05100270e+03, 2.09199736e+03, 4.38724419e+00, 3.98166105e+02,
               7.60000000e-04, 4.99995596e-02, 5.00000000e-03, 1.79687720e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '16': np.array([1.01020405e+03, 2.99035342e+02, 6.85450975e+01, 6.00000000e+05,
               2.02800148e+03, 2.06800000e+03, 3.16442794e+00, 1.73470928e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
       }
    elif args.band == 'H':
        telparsdict = { '1': np.array([1.01099787e+03, 2.77972664e+02, 1.74812183e+01, 6.00000000e+05,
                1.76500270e+03, 1.84299940e+03, 2.59940438e+01, 3.67500000e+02,
                7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.43429720e+00,
                2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
                5.00000000e-09, 1.56000000e-07]),
                        '2': np.array([1.01488875e+03, 2.81026692e+02, 5.96311953e+01, 6.00000000e+05,
               1.76800166e+03, 1.80500074e+03, 1.22804054e+01, 2.59351413e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.33456803e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '3': np.array([1.01303120e+03, 2.85356979e+02, 4.55001103e+01, 6.00000000e+05,
               1.75099969e+03, 1.78799814e+03, 1.24741562e+01, 2.10505343e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.74585431e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '4': np.array([1.01010048e+03, 2.79724074e+02, 4.64849179e+01, 6.00000000e+05,
               1.73399959e+03, 1.77100071e+03, 1.89163537e+01, 2.06761848e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.04884875e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '5': np.array([1.01201335e+03, 2.80571184e+02, 5.12883875e+01, 6.00000000e+05,
               1.71800032e+03, 1.75399961e+03, 1.61496968e+01, 2.30222114e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.88354988e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '6': np.array([1.01002060e+03, 2.79834802e+02, 4.92317644e+01, 6.00000000e+05,
               1.70199977e+03, 1.73800017e+03, 1.94448331e+01, 1.44128653e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.78512125e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '7': np.array([1.02349233e+03, 2.82841942e+02, 4.84055036e+01, 6.00000000e+05,
               1.68600077e+03, 1.72200035e+03, 1.56440323e+01, 2.29711325e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.75546431e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '8': np.array([1.02535183e+03, 2.84690731e+02, 6.38909473e+01, 6.00000000e+05,
               1.67000105e+03, 1.70600005e+03, 8.62895867e+00, 2.93112111e+03,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.08243949e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '9': np.array([1.01005569e+03, 2.85368493e+02, 4.69600789e+01, 6.00000000e+05,
               1.65500012e+03, 1.68999935e+03, 1.34436019e+01, 6.84155129e+01,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.76576528e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '10': np.array([1.01048567e+03, 2.97058150e+02, 4.33608311e+01, 6.00000000e+05,
               1.64000122e+03, 1.67500013e+03, 6.53445899e+00, 4.08655037e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.90785867e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '11': np.array([1.02491139e+03, 2.85559174e+02, 4.08532668e+01, 6.00000000e+05,
               1.62600054e+03, 1.65999850e+03, 2.00545756e+01, 4.11449098e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.88659154e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '12': np.array([1.01003290e+03, 3.09888667e+02, 4.54061158e+01, 6.00000000e+05,
               1.61100114e+03, 1.64599985e+03, 3.64022137e+00, 3.81309830e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.91039536e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '13': np.array([1.01342132e+03, 3.08696226e+02, 4.20101993e+01, 6.00000000e+05,
               1.59700075e+03, 1.63200043e+03, 1.87270207e+00, 3.92945895e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.57637068e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '14': np.array([1.01388819e+03, 3.08306291e+02, 3.84989496e+01, 6.00000000e+05,
               1.58299931e+03, 1.61799908e+03, 3.67128385e+00, 4.27356912e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.23969000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '15': np.array([1.01014958e+03, 3.09907179e+02, 4.01340902e+01, 6.00000000e+05,
               1.57000020e+03, 1.60399922e+03, 3.71266195e+00, 4.13591102e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 4.45900885e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '16': np.array([1.02110710e+03, 2.95717915e+02, 5.10020728e+01, 6.00000000e+05,
               1.55600027e+03, 1.59000003e+03, 3.33875074e+01, 3.23470908e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 4.38118995e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '17': np.array([1.01016410e+03, 3.07744028e+02, 4.29708519e+01, 6.00000000e+05,
               1.54300110e+03, 1.57699913e+03, 3.10294397e+00, 3.89081633e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '18': np.array([1.02523331e+03, 2.80364698e+02, 4.91323141e+01, 6.00000000e+05,
               1.52999979e+03, 1.56400054e+03, 2.00474401e+01, 3.58065516e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '19': np.array([1.01003349e+03, 2.81022630e+02, 4.79869931e+01, 6.00000000e+05,
               1.51799946e+03, 1.55100021e+03, 2.00658592e+01, 3.75629220e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.80000000e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '20': np.array([1.01002788e+03, 2.79946092e+02, 4.98150515e+01, 6.00000000e+05,
               1.50500004e+03, 1.53799890e+03, 1.98874188e+01, 3.63105569e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.08959964e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '21': np.array([1.01198904e+03, 2.75314098e+02, 4.99899615e+01, 6.00000000e+05,
               1.49300048e+03, 1.52599904e+03, 2.58987815e+01, 5.00528381e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.74730881e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '22': np.array([1.01503380e+03, 2.77013372e+02, 5.99130176e+01, 6.00000000e+05,
               1.48099946e+03, 1.51400039e+03, 1.72525774e+01, 2.29210179e+01,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 2.36183449e+00,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
                        '23': np.array([1.01136214e+03, 2.86167482e+02, 5.33964947e+01, 6.00000000e+05,
               1.47015399e+03, 1.50200068e+03, 9.83566652e+00, 2.71708854e+02,
               7.60000000e-04, 5.00000000e-02, 5.00000000e-03, 1.06492323e-01,
               2.10000000e+05, 0.00000000e+00, 5.00000000e-09, 5.00000000e-09,
               5.00000000e-09, 1.56000000e-07]),
}

    fitterL = TelluricFitter(debug=False, print_lblrtm_output=args.debug)

    NSO_props = {"latitude": 31.958, "altitude":2.096} #alt in km
    fitterL.SetObservatory(NSO_props)

    dataL = DataStructures.xypoint(x=newwave1, y=None, cont=None, err=None)

    parfittedL = telparsdict[str(order)]
    paramsL = {}
    for k in range(len(names)):
        paramsL[names[k]] = np.float(parfittedL[k])

    paramsL['wavestart'] = dataL.x[0] # nm
    paramsL['waveend']   = dataL.x[-1] # nm
    fitterL.AdjustValue(paramsL)
    fitterL.ImportData(data2)

    try:
        if args.debug:
            modelL = fitterL.GenerateModel(parfittedL, nofit=True, air_wave=False)
        else:
            modelL = suppress_GenerateModel(fitterL, parfittedL, args)


    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan],[np.nan],[np.nan]

    global x, satmLivGen, watm_Liv,satm_Liv;

    satmTel    = rebin_jv(model2.x*10, model2.y, newwave1*10, True, logger=logger) # nm --> AA
    satmLivGen = rebin_jv(modelL.x*10, modelL.y, newwave1*10, True, logger=logger) # nm --> AA
    watmLivGen = newwave1.copy() ; watmLivGen*=10 # nm --> AA

    # Fit wavelength scale to Telfit'd Livingston
    x = np.arange(len(satmLivGen))
    initguess = np.polyfit(x, watmLivGen, 6)

#    print('watmLivGen=\n', watmLivGen)
#    print('x=\n', x)

    watm_Liv  = inparam.watm[ (inparam.watm > watmLivGen[0]+1) & (inparam.watm < watmLivGen[-1]-1) ]
    satm_Liv  = inparam.satm[ (inparam.watm > watmLivGen[0]+1) & (inparam.watm < watmLivGen[-1]-1 )]
    dpar = np.abs(initguess)*10
    dpar[-1] = 5

    waveparfit = wavefit(initguess, dpar)
    f = np.poly1d(waveparfit)
    wavefitted = f(x)

    satmTel[(satmTel < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi

    return wavefitted, satmTel, names, parfitted, wcont1, cont1
