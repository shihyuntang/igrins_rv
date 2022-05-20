# Step2: Initial Convergence

*Required if the average RV of the target star is unknown to > 5 km/s precision. Performs an abbreviated analysis of the target star observations in order to converge to coarsely accurate RVs, which will be used as starting points for the more precise analysis in the next step; simultaneously does the same for target star's vsini. Only a single echelle region is used, and only one B observation from a given exposure (e.g. AB or ABBA) is analyzed.*
***

It is recommended that the user choose what they think may be the most precise of the wavelength regions defined in Step 1 for Step 2 (presumably whichever has the greatest number of strong stellar and telluric lines). This can be specified in the `-l_use` keyword.

This code can be run on both RV standard stars and stars where RV variability is expected. The difference appears in how you specify the keywords `-g` and/or `-gX`. Note that these keywords also allow you to easily use the RV results of a previous Step 2 run as the input RV guesses for a new run.


## The terminal command for step2
Call as:
```shell
(igrins_rv) ~$ python main_step2.py [target name] [keywords]
```

> NOTE: Not all keywords have default values, meaning that **the code will not run without certain keywords explicitly specified**. These keywords are **-i**, **-temp**, **-logg**, and either one of **-g** or **-gX**. The user must be sure to specify keywords as needed if their desired run of the code varies from the defaults, as well (e.g. if they want to run H band instead of K, or use a different number of CPU threads).

*&dagger; An example*:
```shell
(igrins_rv) ~$ python main_step2.py GJ281 -HorK K -Wr 1 -i 10 -v 5 -g 25 -t synthetic -temp 4000 -logg 4.5 -c 4 -plot -l_use 6
```
> When this code is done executing, it would print the mean RV and vsin(i) of the input observations. The measurements of each individual observation would be saved in a `Initguesser_results_[run_number].csv` file, under `./Output/[target_name]_[band]/`.

> Say in this case, for the sake of example, the code output a mean RV of 20.2 km/s and vsin(i) of 2.5 km/s. If you were then comfortable with these estimates (perhaps std of these estimates are to your liking, or the fit plots all look nice), you would proceed to Step 3. On the other hand, if you wanted to refine these estimates, you would run Step 2 again, but this time with its starting guesses for RV and vsin(i) given by the estimates just produced:

```shell
(igrins_rv) ~$ python main_step2.py GJ281 -HorK K -Wr 1 -i 2.5 -v 3 (-g 20.2 OR -gX 1) -t synthetic -temp 4000 -logg 4.5 -c 36 -plot -l_use 6
```

> where -g 20.2 would apply a singular initial guess RV uniformly to all observations, while -gX 1 would apply a different starting guess for each observation based on the contents of `Initguesser_results_1.csv`. The former is used for targets where little intrinsic variance is expected (like RV standards); the latter is used for targets where RV variance is expected. The -g keyword is also used the first time the user runs Step 2, since presumably even if the RV variance is expected, the user only has an idea of the average absolute RV of the star.

### Keywords for step2

* -HorK : Specifies which band to process, default = K 
* -Wr : Specifies which wavelength regions file in `./Input/UseWv/WaveRegions_X` to use, default is those used in methods paper (-Wr 1)
* -l_use : Specifies which single echelle order will be analyzed. Default is the first in WRegion list.
* -SN : Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 
* -i : Initial vsini guess (km/s)
* -v : Range of allowed vsini variation during optimization (if set to zero vsini will be held constant), default = 5.0 km/s
* -g : Initial guess for RV (km/s) that will be uniformly applied to all observations. Use -gX instead if you want to reference an Initguesser_results file from a previous run of Step 2, which will have a different initial guess for each observation. 
* -gX : The number, X, that refers to the `./*targname/Initguesser_results_X` file you wish to use for initial RV guesses
* -t : Stellar template kind. Pick from 'synthetic' or 'user' (which then needs to be supplied in `./Engine/user_templates`). Default = 'synthetic'"
* -temp : Synthetic stellar template Teff (K). See `./Engine/syn_template` for all those available.
* -logg : Synthetic stellar template log(g) (cgs). See `./Engine/syn_template` for all those available.
* -c : Specifies number of CPU threads to use, default is 1/2 of available ones
* -n_use : Specifies an array of observations you wish to process here (e.g. [20181111,20181112]) if you don't want to process all the observations in the `./Input/*target/` folder
* -plot : If set, will generate plots of fitting results under `./Output/[target_name][band]/figs/main_step2[band]_[run_number]`
* -DeBug : If set, DeBug logging will be output, as well as (lots of) extra plots


## Outputs from step2
### --RVs and vsin(i)s--
Output files are saved as `./Output/[target_name]_[band]/Initguesser_results_[run_number].csv`. The mean RV and vsin(i) of the observations are printed from the terminal.

### --Plots from step2--

If the -plot keyword was specified, plots of individual spectral fits will be under `./Output/[target_name]_[band]/figs/main_step3_[band]_[run_number]`. The plots with the prefix:
* "parfit" show the full spectral fit to the data
* "parfitS" shows solely the stellar template contribution
* "parfitT" shows solely the telluric template contribution