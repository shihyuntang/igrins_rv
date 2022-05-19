# Step4: vsini Uncertainty Calculation

*Updates RV uncertainty estimates to take into account uncertainty in vsini.*
***

This step takes two runs of Step 3, one with vsini held fixed at the best guess value and one with vsini held fixed at the best guess value plus or minus one sigma, and uses the difference between the two to produce updated RVs and uncertainties. For the most part, the uncertainties should change little (1 m/s), but for high vsini (> 15 km/s) objects, they may increase by ~5 m/s or so.
                                     

## The terminal command for step4
Call as:
```shell
(igrins_rv) ~$ python main_step4.py [target name] [keywords]
```

With keywords

* -run1 : First step3 run that will be used, the one with vsini held fixed AT THE BEST GUESS. Takes the string that suffixes `RV_results_`, e.g. for `RV_results_1`, you would set this to `1`.
* -run2 : Second step3 run that will be used, the one with vsini held fixed at the best guess PLUS OR MINUS SIGMA.
* -HorK : Which band to process? H or K?. Default = K

Keywords -run1 and -run2 are mandatory.


## Outputs from step4
### --RVs (step4)--
Output files are saved to `./Output/[target_name]_[band]/RV_results_[run_number_1]_[run_number_2]_combined`. The primary output is the `RVresultsSummary.fits/csv` file, which contains the following columns in its table:
* 'NIGHT' : Designation of observation based on PrepData files, either in format YYYYMMDD or YYYYMMDD_XXXX, where XXXX is for Multi-RVs per Night mode (see Setup: PrepData Files page)
* 'JD' : Julian date
* 'RVfinal' - Final RV estimate for each 'NIGHT'
* 'STDfinal' - Final uncertainty estimates of RVs
* 'VSINI' - Final vsin(i) estimate for each 'NIGHT'

This is the same as in Step 3.