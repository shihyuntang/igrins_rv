# Step1: Telluric Modelling
***
*Defines the wavelength regions to be analyzed; generates a synthetic, high-resolution telluric template for use in later model fits on a night by night basis.*
***

### Input data
If you got 1D spectra from running `igrins plp`, please make sure you've moved your folder \
`plp-master/final_A_B_spec/[target_name]` to `./igrins_rv-master/input/`.

The folder structure should look like:
```
igrins_rv(-master)
├── Engine
│   └── ...
│
├── Input
│   ├── Prepdata
│   │   ├── Prepdata_A0_[targname].txt
│   │   └── Prepdata_targ_[targname].txt
│   ├── [targname]
│   │   ├── YYYYMMDD
│   │   ├── YYYYMMDD
│   │   └── ...
│   └── UseWv
│       └── ... 
├── Output
│   └── ... 
└── ... 
```

###  Defining analysis regions
By default, **IGRINS RV** analyzes portions of 5 echelle orders in the H band and 4 in the K band. These regions were handpicked based off visual inspections of preliminary fits to target star spectra.

**IGRINS RV** comes with these regions pre-loaded, but users can also customize their own lists of wavelength regions they would like the code to fit; this is recommended if the user is studying target stars are different spectral types than those presented in the paper linked at the start of this wiki (F,G,K,M). As part of Step 1, **IGRINS RV** will automatically take the input list of wavelength ranges and convert it into the echelle orders and pixel ranges that will be fit as part of RV estimation. 

If you do choose to select your own wavelength regions, be careful! There are a variety of criteria that should be considered, and you will need to monitor a RV standard to determine the code's precision with the specified regions.

To use the default wavelength regions, do not specify any value for the keyword `-Wr`, or specify the default value of `1`.

To use your own wavelength regions, go to `./Input/UseWv/` and create your own `WaveRegions_X_Y.csv` file, where `X` is a greater number than already present and `Y` is the band (H or K) you wish to analyze. In the `.csv` file, list the orders and wavelength ranges (in microns) you'd like to analyze, using a default WaveRegions file as a template. Then specify `-Wr X`, where `X` is the number of the file you created, when you run Step 1.


###  Generating telluric templates
The bus mostly drives itself here, with nothing required beyond the PrepData files and the wavelength regions. **IGRINS RV** will process all the B frames first for all orders and nights, then all the A frames.

### The terminal command
Call as:
```shell
(igrins_rv) ~$ python main_step1.py [target name] [keywords]
```
NOTE: In Step 1, the keywords all have values that they will automatically default to, and so the code will technically run without any specified. However, the user must be sure to specify keywords as needed if their desired run of the code varies from the defaults (e.g. if they want to run H band instead of K, or use a different number of CPU threads).

*&dagger; For example*:
```shell
(igrins_rv) ~$ python main_step1.py GJ281 -HorK K -Wr 1 -c 4 -plot
```

With keywords 

* -HorK : Specifies which band to process, default = K 
* -Wr : Specifies which wavelength regions file in `./Input/UseWv/WaveRegions_X` to use, default is those used in methods paper (-Wr 1)
* -c : Specifies number of CPU threads to use, default is 1/2 of available ones
* -SN: Specifies spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50
* -n_use : Specifies an array of nights you wish to process here (e.g. [20181111,20181112]) if you don't want to process all the nights in the `./Input/*target/` folder
* -plot : If set, will generate plots of A0 fitting results under `./Output/[target_name]_[band]/A0Fits/figs_[band]/`
* -DeBug : If set, DeBug logging will be output, as well as (lots of) extra plots


### Outputs

The synthetic telluric templates are saved to `./Output/[target_name]_[band]/A0Fits` as `[NIGHT]A0_[A/B nodding]treated_[band].fits`

Visualization of the `Telfit` fitting results can be found under `./Output/[target_name]_[band]/A0Fits/figs_[band]`. Users will mostly (or only) need to check the `A0Telfit_OrderXX_[NIGHT]_[A/B nodding].png`, which shows the final fitting result.
