
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy import units
from Engine.rebin_jv import rebin_jv
from telfit import TelluricFitter, DataStructures
import nlopt
import matplotlib.patches as mpatches


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
        return(offset + slope*x + kurt*((x-mu)**2) + ((scale*np.sqrt(2*np.pi))*np.exp(-.5*((x-mu)/sigma)**2)))
    return innerfit

def gauss(x,mu,sigma,offset,scale,slope,kurt):
    return(offset+ slope*x + kurt*((x-mu)**2) + ((scale*np.sqrt(2*np.pi))*np.exp(-.5*((x-mu)/sigma)**2)))


def wavefunc(par,grad):

    global watm_Liv, satm_Liv, satmLivGen, x;

    # This function takes Telfitted template and uses an input wavelength solution to rebin it for direct comparison with Livingston.


    #Make the wavelength scale
    f = np.poly1d(par)
    w = f(x)

    if w[-1] < w[0]:
        return 1e20

    satmTel2 = rebin_jv(w,satmLivGen,watm_Liv,False)

    return np.sum((satm_Liv - satmTel2)**2.)


def wavefit(par0, dpar0):
    # NLopt convenience function.
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, 7)
    opt.set_min_objective(wavefunc)
    lows = par0-dpar0
    highs = par0+dpar0
    opt.set_lower_bounds(lows)
    opt.set_upper_bounds(highs)
    opt.set_maxtime(600) #seconds
    # Quit optimization based on relative change in output fit parameters between iterations.
    # Choosing smaller change tolerance than 1e-6 has demonstrated no improvement in precision.
    opt.set_xtol_rel(1e-6)
    parfit = opt.optimize(par0)
    return parfit



def telfitter(watm_in, satm_in, a0ucut, inparam, night, order):

    # Code to produced fitted telluric template. How and why it works is detailed in comments throughout the code.

    os.environ['PYSYN_CDBS'] = inparam.cdbsloc

    fitter = TelluricFitter(debug=False)

    #Set the observatory location with a keyword
    DCT_props = {"latitude": 34.744, "altitude":2.36} #altitude in km
    McD_props = {"latitude": 30.71, "altitude": 2.07}
    if inparam.obses[night] == 'DCT':
        fitter.SetObservatory(DCT_props)
    elif inparam.obses[night] == 'McD':
        fitter.SetObservatory(McD_props)
    else:
        sys.exit('TELFIT OBSERVATORY ERROR')
#        print('TELFIT OBSERVATORY ERROR')
#        print(breaker) #force quit)

    # Read in data
    data = DataStructures.xypoint(x=watm_in*units.angstrom, y=satm_in, cont=None, err=a0ucut)

    # DCT data has parameters describing night of observation that the McDonald data does not.
    if inparam.zds[night] != 'NOINFO': # If such information is available:

        angle       = float(inparam.zds[night])  #Zenith distance
        pressure    = float(inparam.press[night])  #Pressure, in hPa
        humidity    = float(inparam.humids[night])  #Percent humidity, at the observatory altitude
        temperature = float(inparam.temps[night])+273.15  #Temperature in Kelvin
        resolution  = 45000.0                          #Resolution lambda/delta-lambda

        # Ideally, we'd fit resolution as well since that varies across the detector.
        # But in practice the Telfit's resolution fits often diverge to unphysically high values.
        # Ultimately, we only want accurate estimates for the chemical compositions, which are unaffacted
        # by fixing the resolution at 45000. The final telluric template we'll be synthesizing from this
        # will be at a set high resolution, anyway, and when we need/have needed to estimate resolution
        # across the detector, we (have) done so via the fits of the telluric and stellar templates to the
        # observed A0 and GJ281 spectra.

        # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
        fitter.FitVariable({"h2o": humidity,"ch4": 1.8,"co": 5e-3})

        #Adjust parameters that will not be fit, but are important
        fitter.AdjustValue({"angle": angle,\
                            "pressure": pressure,\
                            "temperature": temperature,\
                            "resolution": resolution,
                            "wavestart": data.x[0],\
                            "waveend": data.x[-1],\
                            "co2": 3.675e2,\
                            "o3": 7.6e-4,\
                            "n2o": 5e-2,\
                            "o2": 2.1e5,\
                            "no": 0.,\
                            "so2": 5e-9,\
                            "no2": 5e-9,\
                            "nh3": 5e-9,\
                            "hno3": 1.56e-7})

        #Set bounds on the variables being fit
        fitter.SetBounds({"h2o": [1.0, 99.0],\
                          "ch4": [.1,  10.0],\
                          "co": [ 1e-6,1e2]})

    else: # If parameters are not in fits file, use initial guesses and letting them vary.
          # Guesses are taken from mean of parameters from DCT GJ281 data.

        resolution = 45000.0                          #Resolution lambda/delta-lambda

        # Only 3 molecules present in chosen IGRINS orders' wavelength range are H2O, CH4, and CO.
        fitter.FitVariable({"h2o": 43.,"ch4": 1.8,"co": 5e-3,
                            "angle": 39., "pressure":1023., "temperature":280.87})

        #Adjust parameters that will not be fit, but are important
        fitter.AdjustValue({"resolution": resolution,\
                            "wavestart": data.x[0],\
                            "waveend": data.x[-1],\
                            "co2": 3.675e2,\
                            "o3": 7.6e-4,\
                            "n2o": 5e-2,\
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
                          "pressure": [1010.,1035.],\
                          "co": [ 1e-6,1e2]})

    try:
        model = fitter.Fit(data=data, resolution_fit_mode="SVD", adjust_wave="model",air_wave=False)
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
    wcont1 = model.x

    if inparam.plotfigs == True:
        plt.figure(figsize=(20,12))
        plt.plot(watm_in,satm_in,color='black',alpha=.6,label='data')
        plt.plot(model.x,model.y*cont1,color='red',alpha=.6,label='model fit')
        plt.plot(model.x,cont1,color='blue',alpha=.6,label='blaze fit')
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Flux (Normalized)')
        plt.legend()
        plt.savefig(inparam.outpath+'/A0 Telfit '+str(night)+'_'+str(order)+'.png')
        plt.clf()
        plt.close()

    ############### Generate template with these parameters but at higher resolution

    names = ["pressure", "temperature", "angle", "resolution",'wavestart','waveend',
                         "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "no",
                         "so2", "no2", "nh3", "hno3"]

    parfitted = np.ones_like(names, dtype=float)
    for k in range(len(names)):
        parfitted[k] = float(fitter.GetValue(names[k]) )


    fitter2 = TelluricFitter(debug=False)

    if inparam.obses[night] == 'DCT':
        fitter2.SetObservatory(DCT_props)
    elif inparam.obses[night] == 'McD':
        fitter2.SetObservatory(McD_props)

    # Compute telluric template with highest resolution of Livingston template.
    # Add extra space at ends to make sure template covers wider range than data.
    Livingston_minimum_wsep = .035
    IGRINS_minimum_wsep     = .130 # <-- This would compute template with IGRINS resolution, sensibly coarser than Livingston

    newwave = np.arange(min(watm_in)-25, max(watm_in)+25, Livingston_minimum_wsep)

    data2 = DataStructures.xypoint(x=newwave*units.angstrom,
                                   y=None,
                                   cont=None,
                                   err=None)
    params = {}
    for k in range(len(names)):
        params[names[k]] = float(parfitted[k])

    params['wavestart'] = data2.x[0]
    params['waveend']   = data2.x[-1]

    fitter2.AdjustValue(params)
    fitter2.ImportData(data2)

    # Call the modeller. On rare occasions, this returns an error. I have no idea what is causing this error, as the
    # FORTRAN readout is quite unhelpful and anyone else who apepars to have experienced this problem had it randomly go away at some point.
    # If this happens, simply deliver NAN arrays, and in later parts of the RV analysis A0 fits from the nearest compatible observation will be used.
    try:
        model2 = fitter2.GenerateModel(parfitted,nofit=True)
    except TypeError:
        return [np.nan], [np.nan], [np.nan], [np.nan],[np.nan],[np.nan]

    watm_save = watm_in.copy(); satm_save = satm_in.copy();
    newwave1 = newwave[(newwave > watm_in[0]-10) & (newwave < watm_in[-1]+10)]

    # Parameters for reproducing Livingston template with Telfit
    telparsdict = {'2':np.array([1.01000565e+03, 2.86022221e+02, 1.02518785e+01, 6.00000000e+05,
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
           5.00000000e-09, 1.56000000e-07])}

    fitterL = TelluricFitter(debug=False)

    NSO_props = {"latitude": 31.958, "altitude":2.096} #alt in km
    fitterL.SetObservatory(NSO_props)

    dataL = DataStructures.xypoint(x=newwave1*units.angstrom, y=None, cont=None, err=None)

    parfittedL = telparsdict[str(order)]
    paramsL = {}
    for k in range(len(names)):
        paramsL[names[k]] = float(parfittedL[k])

    #paramsL['wavestart'] = dataL.x[0]
    #paramsL['waveend']   = dataL.x[-1]
    #fitterL.AdjustValue(paramsL)
    #fitterL.ImportData(dataL)

    params['wavestart'] = dataL.x[0]
    params['waveend']   = dataL.x[-1]
    fitterL.AdjustValue(paramsL)
    fitterL.ImportData(data2)

    modelL = fitterL.GenerateModel(parfittedL,nofit=True)

    global x, satmLivGen, watm_Liv,satm_Liv;

    satmTel    = rebin_jv(model2.x, model2.y, newwave1,True)
    satmLivGen = rebin_jv(modelL.x, modelL.y, newwave1,True)
    watmLivGen = newwave1.copy()

    # Fit wavelength scale to Telfit'd Livingston
    x = np.arange(len(satmLivGen))
    initguess = np.polyfit(x, watmLivGen, 6)

#    print('watmLivGen=\n', watmLivGen)
#    print('x=\n', x)

    watm_Liv  = inparam.watm[ (inparam.watm > watmLivGen[0]+1) & (inparam.watm < watmLivGen[-1]-1) ]
    satm_Liv  = inparam.satm[ (inparam.watm > watmLivGen[0]+1) & (inparam.watm < watmLivGen[-1]-1 )]
    dpar = abs(initguess)*100
    dpar[-1] = 5

    waveparfit = wavefit(initguess, dpar)
    f = np.poly1d(waveparfit)
    wavefitted = f(x)

    satmTel[(satmTel < 1e-4)] = 0. # set very low points to zero so that they don't go to NaN when taken to an exponent by template power in fmodel_chi

    return wavefitted, satmTel, names, parfitted, wcont1, cont1
