from Engine.importmodule import *
import os, nlopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from   Engine.rebin_jv import rebin_jv
from   telfit import TelluricFitter, DataStructures
from scipy.optimize import leastsq


# From Telfit:
def Poly( pars, middle, low, high, x):
        """
        Generates a polynomial with the given parameters
        for all of the x-values.
        x is assumed to be a np.ndarray!
         Not meant to be called directly by the user!
        """
        xgrid = (x - middle) / (high - low)  # re-scale
        return chebyshev.chebval(xgrid, pars)

# From Telfit:
def WavelengthErrorFunctionNew(pars, datax,datay, telluricx, telluricy, maxdiff=0.05):
        """
        Cost function for the new wavelength fitter.
        Not meant to be called directly by the user!
        """
        dx = Poly(pars, np.median(datax), min(datax), max(datax), datax)
        penalty = np.sum(np.abs(dx[np.abs(dx) > maxdiff]))
        #retval = (datay - model(datax + dx)) + penalty
        retval = (datay - rebin_jv(telluricx,telluricy,datax + dx,False)) + penalty
        return retval

# From Telfit:
def FitWavelengthNew( data_originalx, data_originaly, telluricx, telluricy, fitorder=3, be_safe=True):
        """
        This is a vastly simplified version of FitWavelength.
        It takes the same inputs and returns the same thing,
        so is a drop-in replacement for the old FitWavelength.

        Instead of finding the lines, and generating a polynomial
        to apply to the axis as x --> f(x), it fits a polynomial
        to the delta-x. So, it fits the function for x --> x + f(x).
        This way, we can automatically penalize large deviations in
        the wavelength.
        """
        #modelfcn = UnivariateSpline(telluricx, telluricy, s=0)
        pars = np.zeros(fitorder + 1)
        if be_safe:
            #args = (data_originalx, data_originaly, modelfcn, 0.05)
            args = (data_originalx, data_originaly, telluricx, telluricy, 0.05)
        else:
            #args = (data_originalx, data_originaly, modelfcn, 100)
            args = (data_originalx, data_originaly, telluricx, telluricy, 100)
        output = leastsq(WavelengthErrorFunctionNew, pars, args=args, full_output=True, xtol=1e-12, ftol=1e-12)
        pars = output[0]

        return partial(Poly, pars, np.median(data_originalx), min(data_originalx), max(data_originalx)), 0.0



#------------
# to suppress print out from Telfit
@suppress_stdout
def suppress_Fit(fitter, data, c_order):
    model = fitter.Fit(data=data, resolution_fit_mode="SVD", adjust_wave="data", air_wave=False, continuum_fit_order=c_order)
    return model

@suppress_stdout
def suppress_GenerateModel(fitter, parfit, args):
    model = fitter.GenerateModel(parfit, nofit=True, air_wave=False)
    return model

# some share default values and fitting boundaries for telfit
_telfit_default_values_dic = {"h2o" : 43.0,
                              "ch4" :  1.8,
                              "co"  :  5e-3,
                              "n2o" :  5e-2,
                              "co2" :  3.675e2,
                              "o3"  :  7.6e-4,
                              "o2"  :  2.1e5,
                              "no"  :  0.0,
                              "so2" :  5e-9,
                              "no2" :  5e-9,
                              "nh3" :  5e-9,
                              "hno3":  1.56e-7,
                              "pressure"   : 1023.0, 
                              "temperature":  280.87,
                              "angle"      :   39.0 
                              }

_telfit_default_vary_bound_dic = {"h2o": [1.0, 99.0],
                                  "ch4": [0.1, 10.0],
                                  "n2o": [1e-5, 1e2],
                                  "co" : [1e-6, 1e2],
                                  "co2": [1.0,  1e4],
                                  "temperature": [ 265.0,  300.0],
                                  "pressure"   : [1010.0, 1035.0],
                                  "angle"      : [   1.0,   75.0]
                                  }


def telfitter(watm_in, satm_in, a0ucut, inparam, night, order, args, masterbeam, c_order,resolutions, logger):
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

    # Define to fitter objects, one for actual fit and one for high res model generation
    fitter  = TelluricFitter(debug=False, print_lblrtm_output=args.debug)
    fitter2 = TelluricFitter(debug = False, print_lblrtm_output = args.debug)

    #Set the observatory location with a keyword
    DCT_props     = {"latitude":  34.744, "altitude": 2.36} # altitude in km
    McD_props     = {"latitude":  30.710, "altitude": 2.07}
    GeminiS_props = {"latitude": -30.241, "altitude": 2.72}

    if inparam.obses[night] == 'DCT':
        fitter.SetObservatory(DCT_props)
        fitter2.SetObservatory(DCT_props)
    elif inparam.obses[night] == 'McD':
        fitter.SetObservatory(McD_props)
        fitter2.SetObservatory(McD_props)
    elif inparam.obses[night] == 'GeminiS':
        fitter.SetObservatory(GeminiS_props)
        fitter2.SetObservatory(GaminiS_props)
    else:
        sys.exit('TELFIT OBSERVATORY ERROR, OLNY SUPPORT DCT, McD & GeminiS IN THIS VERSION!')

    # Read in data
    # do not use astropy units in telfit to speed things up. telfit default wavelength units is in nm
    watm_in = watm_in/10 # AA --> nm
    data = DataStructures.xypoint(x=watm_in, y=satm_in, cont=None, err=a0ucut) # input wavelength in nm
    
    # set spectrum resolution
    # Ideally, we'd fit resolution as well since that varies across the detector.
    # But in practice the Telfit's resolution fits often diverge to unphysically high values.
    # Ultimately, we only want accurate estimates for the chemical compositions, which are unaffacted
    # by fixing the resolution at 45000. The final telluric template we'll be synthesizing from this
    # will be at a set high resolution, anyway, and when we need/have needed to estimate resolution
    # across the detector, we (have) done so via the fits of the telluric and stellar templates to the
    # observed A0 and GJ281 spectra.
    
    resolution_min, resolution_med, resolution_max = resolutions

    # Define dicts of parameters to fit, adjust to set values, and bounds of fit params
    tofit = {}; toadjust = {}; tobound = {};

    names = ["pressure", "temperature", "angle", "resolution", "wavestart", "waveend",
             "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "no", "so2", "no2", "nh3", "hno3"]
        
    # DCT data has parameters describing night of observation that the McDonald and GeminiS data does not.
    if inparam.temps[night] != 'NOINFO': # If temperature info. is available:

        toadjust["angle"]       = np.float(inparam.zds[night]    )        # Zenith distance
        toadjust["pressure"]    = np.float(inparam.press[night]  )        # Pressure, in hPa
        toadjust["temperature"] = np.float(inparam.temps[night]  )+273.15 # Temperature in Kelvin
        tofit["h2o"]            = np.float(inparam.humids[night] )        # Percent humidity, at the observatory altitude
        tobound["h2o"]          = _telfit_default_vary_bound_dic["h2o"]
        
    elif inparam.zds[night] != 'NOINFO': 
        # If GeminiS data, some but not all parameters are in fits file.
        # If parameters are not in fits file, use initial guesses and letting them vary.
        # Guesses are taken from mean of parameters from DCT GJ281 data.

        toadjust["angle"]       = np.float(inparam.zds[night]    )        # Zenith distance
        tofit["h2o"]            = np.float(inparam.humids[night] )        # Percent humidity, at the observatory altitude
        tobound["h2o"]          = _telfit_default_vary_bound_dic["h2o"]
    else:
        pass
        
    fitthis = False # For each of these parameters, check if they haven't been fetched from fits file. If not, set to default values.
    for par in ["pressure", "temperature", "angle","h2o"]:
        if par != "h2o":
            try:
                toadjust[par]
            except KeyError:
                fitthis = True
        else:
            try:
                tofit[par]
            except KeyError:
                fitthis = True
        if fitthis:
            tofit[par]   = _telfit_default_values_dic[par]
            tobound[par] = _telfit_default_vary_bound_dic[par]
        fitthis = False
        
    # Depending on wavelength region, set different molecules to be fit
    if (3 < order < 9) & (args.band == 'K'):
        molstofit = ["ch4","co","n2o"]   # Only 4 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, N2O, and CO.
    elif (order >= 9 or order <= 3) & (args.band == 'K'): 
        molstofit = ["ch4","co2","n2o"]   # Only 4 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, N2O and CO2.
    elif args.band =='H':
        molstofit = ["ch4","co","co2","n2o"]
        
    for mol in molstofit:
        tofit[mol]   = _telfit_default_values_dic[mol]
        tobound[mol] = _telfit_default_vary_bound_dic[mol]
        
    tofit["resolution"]     = resolution_med
    tobound["resolution"] = [resolution_min,resolution_max]

    # Rest of molecules are set to defaults
    for par in names[4:]:
        if par not in molstofit:
            if par == 'wavestart':
                toadjust[par] = data.x[0] - 0.001
            elif par == 'waveend':
                toadjust[par] = data.x[-1] + 0.001
            else:
                toadjust[par] = _telfit_default_values_dic[par]
            
    fitter.FitVariable(tofit)
    fitter.AdjustValue(toadjust) 
    fitter.SetBounds(tobound) 
            
    try:
        if args.debug:
            model = fitter.Fit(data=data, resolution_fit_mode="SVD", adjust_wave="data", air_wave=False, continuum_fit_order=c_order)
        else:
            model = suppress_Fit(fitter, data, c_order)
    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]

    '''
    Note
    ----
    resolution_fit_mode = SVD ought to give faster, more accurate fits for the deep telluric lines we mostly see in K band
    air_wave = False, because data in vacuum wavelengths
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
                4) To properly calibrate this template in wavelength space, we fit it to a Telfit'd version of the Livingston telluric template, only allowing
                    the wavelength solution to vary.
    '''

    #Get the improved continuum from the fitter
    cont1  = fitter.data.cont
    wcont1 = model.x*10 # nm-->AA

    if args.plotfigs:
        fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=250)

        axes.plot(10*watm_in, satm_in,       color='black',    alpha=.8, label='data',      lw=0.7)
        axes.plot(10*model.x, model.y*cont1, color='tab:red',  alpha=.8, label='model fit', lw=0.7)
        axes.plot(10*model.x, cont1,         color='tab:blue', alpha=.8, label='blaze fit', lw=0.7)

        axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        axes.yaxis.set_minor_locator(AutoMinorLocator(2))
        axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
        axes.set_ylabel(r'Flux',       size=6, style='normal' , family='sans-serif' )
        axes.set_xlabel(r'Wavelength [$\rm \AA$]', size=6, style='normal' , family='sans-serif' )
        axes.legend(fontsize=5, edgecolor='white')
        axes.set_title('A0Telfit_Order{}_{}_{}.png'.format(order, night, masterbeam),
                         size=6, style='normal', family='sans-serif')

        fig.savefig('{}/figs_{}/A0Telfit_Order{}_{}_{}.png'.format(inparam.outpath, args.band, order, night, masterbeam),
                    format='png', bbox_inches='tight', overwrite=True)

    ############### Generate telluric template with these parameters but at higher resolution

    # Get best fit params in list and dict
    parfittedL = np.ones_like(names, dtype=float); parfittedD = {};
    for k in range(len(names)):
        parfittedL[k]        = np.float(fitter.GetValue(names[k]) )
        parfittedD[names[k]] = np.float(fitter.GetValue(names[k]) )
    parfittedD0 = deepcopy(parfittedD)

    # Compute telluric template with highest resolution of Livingston template.
    # Add extra space at ends to make sure template covers wider range than data.
    Livingston_minimum_wsep = .035/10 # <-- Set this to be lower to increase resolution of output template!
    IGRINS_minimum_wsep     = .130 # This would compute template with IGRINS resolution, sensibly coarser than Livingston

    newwave = np.arange(np.min(watm_in)-2.5, np.max(watm_in)+2.5, Livingston_minimum_wsep) # in nm

    data2 = DataStructures.xypoint(x=newwave,
                                   y=None,
                                   cont=None,
                                   err=None)

    parfittedD['wavestart'] = data2.x[0] 
    parfittedD['waveend']   = data2.x[-1]

    fitter2.AdjustValue(parfittedD)
    fitter2.ImportData(data2)

    # Call the modeler. On rare occasions, this returns an error as noted on https://github.com/kgullikson88/Telluric-Fitter/issues/9#issuecomment-123074731
    # The cause of this issue is with the gfortran compiler.
    # If this happens, we deliver NAN arrays, and in later parts of the RV analysis, A0 fits from the nearest compatible observation will be used.
    # HOWEVER, the better solution is to use the ifort compiler.
    try:
        if args.debug:
            model2 = fitter2.GenerateModel(parfittedD, nofit=True, air_wave=False)
        else:
            model2 = suppress_GenerateModel(fitter2, parfittedD, args)

    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]

    wavefitted = model2.x.copy()*10.
    satmTel = model2.y.copy()
    satmTel[(satmTel < 1e-6)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi

    mastermols = names[6:]; molwaves = {}; molfluxes = {};

    # For every molecule being fit, copy the fit params, set all other molecules to 0, and produce a model with only that molecule
    molsout = molstofit+["h2o"]

    for mol in molsout:
        parfittedDmol = deepcopy(parfittedD0)
        for othermol in mastermols:
            if othermol != mol:
                parfittedDmol[othermol] = 0.

        fitter2.AdjustValue(parfittedDmol)
    
        if args.debug:
            modelM = fitter2.GenerateModel(parfittedDmol, nofit=True, air_wave=False)
        else:
            modelM = suppress_GenerateModel(fitter2, parfittedDmol, args)

        molwaves[mol]  = modelM.x.copy()*10.
        molfluxes[mol] = modelM.y.copy()

        if args.debug:
            fig, axes = plt.subplots(1, 1, figsize=(6,3), facecolor='white', dpi=250)
            axes.plot(10*watm_in, satm_in,       color='black',    alpha=.8, label='data',      lw=0.7)
            axes.plot(10*model.x, rebin_jv(modelM.x,modelM.y,model.x,False)*cont1, color='tab:red',  alpha=.8, label='model fit', lw=0.7)
            axes.xaxis.set_minor_locator(AutoMinorLocator(5))
            axes.yaxis.set_minor_locator(AutoMinorLocator(2))
            axes.tick_params(axis='both', which='both', labelsize=6, right=True, top=True, direction='in')
            axes.set_ylabel(r'Flux',       size=6, style='normal' , family='sans-serif' )
            axes.set_xlabel(r'Wavelength [$\rm \AA$]', size=6, style='normal' , family='sans-serif' )
            axes.legend(fontsize=5, edgecolor='white')
            axes.set_title('A0Telfit_Order{}_{}_{}_{}.png'.format(order, night, masterbeam, mol),
                         size=6, style='normal', family='sans-serif')
            fig.savefig('{}/figs_{}/A0Telfit_Order{}_{}_{}_{}.png'.format(inparam.outpath, args.band, order, night, masterbeam, mol),
                    format='png', bbox_inches='tight', overwrite=True)

    return wavefitted, satmTel, names, parfittedL, wcont1, cont1, molsout, molwaves, molfluxes
