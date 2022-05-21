# Step3: Full Analysis

*Performs a full analysis of each target star observation to produce accurate and precise RVs.*
***

Now, all the wavelength regions defined in Step 1 are used, and the code analyzes each exposure that is part of a given observation separately (this allows estimates of the RV uncertainties). In other words, if a star is observed with the ABBA beam pattern, Step 2 analyzes just one B observation, while Step 3 analyzes each A and B separately. 

Unless the target vsini is already known to high accuracy, an initial run of Step 3 in which vsini is allowed to vary is required. This will provide an estimate of vsini, which will then be plugged in as a fixed value in the second run of Step 3. 

If the user seeks the best possible RV uncertainty estimates, or if their target star has a relatively high vsini (> 10 km/s), they must run Step 3 once with vsini held fixed at its estimated value and once with vsini held fixed at this value plus or minus one sigma. The minor differences in the RVs of the two runs (as low as < 1 m/s and as high as 5 m/s) can then be incorporated into the final uncertainties through Step 4. If vsini is already well-known, it is not necessary to run Step 3 more than once, as the code fully converges to the final RVs (within uncertainty) through just one run. 

Like Step 2, this code can be run on both RV standard stars (use "-mode STD") and stars where RV variability is expected (use "-mode TAR"). The difference appears in how you specify the keywords `-mode`, `-g`, `-gS`, and/or `-gX`. Note that these keywords also allow you to easily use the RV results of a previous Step 3 run as the input RV guesses for a new run.

## The terminal command for step3
Call as:
```shell
(igrins_rv) ~$ python main_step3.py [target name] [keywords]
```

> NOTE: Not all keywords have default values, meaning that **the code will not run without certain keywords explicitly specified**. These keywords are **-i**, **-temp**, **-logg**, and either one of **-g** or **-gX**. The user must be sure to specify keywords as needed if their desired run of the code varies from the defaults, as well (e.g. if they want to run H band instead of K, or use a different number of CPU threads).

*&dagger; An example for STD mode (for the example data, GJ281)*
```shell
(igrins_rv) ~$ python main_step3.py GJ281 -mode STD -HorK K -Wr 1 -nAB 1 -i 2.7561 -v 5 -g 20.1478 -t synthetic -temp 4000 -logg 4.5 -c 4 -plot
```

An example for TAR mode
```shell
(igrins_rv) ~$ python main_step3.py HD189733 -mode TAR -HorK K -Wr 1 -i 4.0276 -v 5 -gS init -gX 1 -c 4 -plot -t synthetic -temp 5000 -logg 4.5
```
> The keyword `-gS init -gX 1` will tell **IGRINS RV** to look for RV initial guesses in the output of Step2, run #1, `Initguesser_results_1.csv`.

> If you want or need to run Step 3 multiple times on the same data set (see top of page) then you can instruct **IGRINS RV** to look for RV initial guesses in the output of Step3. The keyword `-gS rvre -gX 1` would fetch initial guesses from `RVresultsSummary_1.csv`, for example: 
```shell
(igrins_rv) ~$ python main_step3.py HD189733 -mode TAR -HorK K -Wr 1 -i 4.3789 -v 0 -gS rvre -gX 1 -c 36 -plot -t synthetic -temp 5000 -logg 4.5
```
<!-- >  In the above command, we have also set the keyword `-abs_out rel`, which will give relative RVs (more precise) rather than absolute RVs (less precise). 

*&dagger; An example for target GJ281 to get the relative RVs*
```shell
(igrins_rv) ~$ python main_step3.py GJ281 -mode STD -HorK K -Wr 1 -nAB 1 -i 3.0755 -v 0 -g 19.9611 -t synthetic -temp 4000 -logg 4.5 -c 4 -plot
``` -->

### Keywords for step3

* -mode : Specifies whether analysis is of RV standard (STD) or a normal target (TAR). 
* -g : For STD star (-mode STD). Initial guess for RV (km/s) that will be uniformly applied to all observations. Given by Step 2 results. Use -gX instead if you want to reference an Initguesser_results file from a previous run of Step 2, which will have a different initial guess for each observation (-mode TAR). 
* -gS : Takes `init` or `rvre`. Source for list of initial RV guesses. `init` means initial RV guesses come from the file `Initguesser_results_X`, where X is specified by -gX. This is a past Step 2 result. OR, specify `rvre` to have initial RV guesses come from the file `RV_results_X`, a past Step 3 result. Use for "-mode TAR" only. 
* -gX : Takes number. For "-mode TAR" only. The number, X, of `./*targname/Initguesser_results_X` OR `./*targname/RV_results_X`, that you wish to use. Prefix determined by -gS.
* -nAB : Specifies minimum number of separate A/B exposures within a set for a given observation (ensures accuracy of uncertainy estimates). Default = 2 for STD, 3 for TAR

* -HorK : Specifies which band to process, default = K 
* -Wr : Specifies which wavelength regions file in `./Input/UseWv/WaveRegions_X` to use, default is those used in methods paper (-Wr 1)
* -SN : Spectrum S/N quality cut. Spectra with median S/N below this will not be analyzed. Default = 50 
* -i : Initial vsini guess (km/s)
* -v : Range of allowed vsini variation during optimization (if set to zero vsini will be held constant), default = 5.0 km/s
* -t : Stellar template kind. Pick from 'synthetic' or 'user' (which then needs to be supplied in `./Engine/user_templates`). Default = 'synthetic'
* -temp : Synthetic stellar template Teff (K). See `./Engine/syn_template` for all those available.
* -logg : Synthetic stellar template log(g) (cgs). See `./Engine/syn_template` for all those available.
<!-- * -abs_out : Set to 'REL' for more precise, relative RVs (default); set to 'ABS' for slightly less precise absolute RVs. -->
* -c : Specifies number of CPU threads to use, default is 1/2 of available ones
* -n_use : Specifies an array of observations you wish to process here (e.g. [20181111,20181112]) if you don't want to process all the observations in the `./Input/*target/` folder
* -plot : If set, will generate plots of fitting results under `./Output/[target_name]_[band]/figs/main_step3_[band]_[run_number]`
* -DeBug : If set, DeBug logging will be output, as well as (lots of) extra plots

## Outputs from step3
### --RVs (step3)--
Output files are saved to `./Output/[target_name]_[band]/RV_results_[run_number]`. The primary output is the `RVresultsSummary.fits` file, which contains the following columns in its table:
* **NIGHT** : Designation of observation based on PrepData files, either in format YYYYMMDD or YYYYMMDD_XXXX, where XXXX is for Multi-RVs per Night mode (see Setup: PrepData Files page)
* **JD** : Julian date
* **RVfinal** - Final RV estimate for each 'NIGHT'
* **STDfinal** - Final uncertainty estimates of RVs
* **VSINI** - Final vsin(i) estimate for each 'NIGHT'

This is likely all you are looking for. 

> the csv file under `./Output/[target_name]_[band]/`, `RVresultsSummary_[run_number].csv` has the same content as the `RVresultsSummary.fits`. 
>
> `RVresultsSummary.fits` file is a fits table, therefore, you can read it by
```python
from astropy.table import Table
data = Table.read('RVresultsSummary.fits', format='fits')
```

### --Plots from step3--

A plot summarizing the final RV measurements will be in the same directory as the RVresultsSummary file. If some of your observations were taken during the time when IGRINS experienced the defocus and some during other times, IGRINS RV will also produce separate plots for each of these epochs.

Additionally, if the -plot keyword was specified, plots of individual spectral fits will be under `./Output/[target_name]_[band]/figs/main_step3_[band]_[run_number]`. The plots with the prefix:
* **parfit** show the full spectral fit to the data
* **parfitS** shows solely the stellar template contribution
* **parfitT** shows solely the telluric template contribution

### --Data-- (Advanced)

If you are using a different means of generating stellar templates than those provided, or have defined your own wavelength regions for analysis, then it is recommended that you apply your version of the pipeline to an RV standard star in order to characterize these variances under the new conditions. In this case, you would run Step 3 in STD mode on an RV standard, and take the output `sigma_method2` array and plug it back into IGRINS RV under `./Engine/classes.py`. Your uncertainties would then be properly characterized when you run science targets with the new setup.

This `sigma_method2` array is printed automatically when Step 3 is done processing in STD mode. It can also be found in the `RVresultsAdvanced.fits` file in the column named `Sigma_method2`. 

The `RVresultsAdvanced.fits` file contains additional columns, but these should not be needed by the user:
* **RVBOX** - 2D array of the RV estimates from different orders for each 'NIGHT'
* **STDBOX** - 2D array of the uncertainty estimates in the RVs from different orders for each 'NIGHT', without method uncertainty added in
* **Sigma_ON2**- 2D array of the variance estimates in the RVs from different orders for each 'NIGHT', with method uncertainty added in

If Step 3 was run in STD mode, there will also be:
* **Sigma_O2** - The variance of RVs within a given order across all exposures
* **Sigma_ABbar2** - 2D array of the median  of the variances of each set of  observations

> `RVresultsAdvanced.fits` file is not a single level fits table, therefore, you will need to open it via
```python
from astropy.io import fits
h = fits.open('RVresultsAdvanced.fits')
```

The output directory should also contain `RVresultsRawBox.fits`, which provides the results of the spectral model fits without any additional processing. Users should not need this.



