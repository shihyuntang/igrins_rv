# Setup: Input 1D spectra (plp)

To enable error estimation, **IGRINS RV** **cannot** get RVs from the standard 1D spectra reduced by [igrins plp v2.2.0](https://github.com/igrins/plp). **You will need the raw fits file, and to perform a specific (but semi-automated) reduction procedure**.

To derive uncertainty estimates for its radial velocity measurements, **IGRINS RV** analyzes each observation from an exposure (a nodding sequence) separately and calculates an RV from each. When it comes to data reduction, this is done by adding the ``--frac-slit`` command flag for A0 standards and the science targets during reduction.

***


## Running IGRINS pipeline version 2.2.0 (plp v2.2.0)
First, open a new terminal and cd into the plp directory. Make sure you are not in **ANY** conda environment. All plp_main_stepXs will open/close the environment for you.

Because we are only modifying **input of the recipe**, the first step is to make sure that ``plp v2.2.0`` is able to run on your machine. See [https://github.com/igrins/plp/wiki/How-to-run-pipeline](https://github.com/igrins/plp/wiki/How-to-run-pipeline) for more detail. 

> If you already have ``igrins plp v2.2.0`` on your machine, you can just download `plp_main_step1.py`, `plp_main_step2.py`, and `plp_main_step2_s.py`. Put these three `.py` files at the same level as the `igr_pipe.py`.

> You can follow the step by step tutorials below with example data of the *GJ281* given under the `./indata/20170216` 

## Modify the raw data and make new recipe
  * Place your `.tmp` recipe file under `./recipes/[target_name]_recipes/`. 
> The `.tmp` recipe file can be obtained via the simple command `python igr_pipe.py prepare-recipe-logs [NIGHT]` built into ``igrins plp`` **IF** you have the observation log under ``indata/[NIGHT]/IGRINS_DT_Log_[NIGHT]-1_[band].txt`` (For more detail, please visit [How to prepare recipe logs](https://github.com/igrins/plp/wiki/How-to-prepare-recipe-logs)). Without the observation log, one will need to make the `.tmp` recipe file manually.

An example of the `.tmp` recipe file:
```python
OBJNAME, OBJTYPE, GROUP1, GROUP2, EXPTIME, RECIPE, OBSIDS, FRAMETYPES
# Available recipes : FLAT, THAR, SKY, A0V_AB, A0V_ONOFF, STELLAR_AB, STELLAR_ONOFF, EXTENDED_AB, EXTENDED_ONOFF
FLAT ON/OFF, FLAT, 11, 1, 30.000000, FLAT,11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30, OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF ON ON ON ON ON ON ON ON ON ON
 , TAR, 59sky, 1, 300.000000, SKY,59, A
HR 1482, STD, 39, 1, 60.000000, A0V_AB, 39 40 41 42, A B B A
GJ 281, TAR, 60, 54, 600.000000, STELLAR_AB,60 61 62 63, A B B A
```
> Example .tmp for GJ281 is provided under: `./recipes/GJ281_recipes/20170216.recipes.tmp`
* !!! There is a "space" in the first column of the SKY !!!

* Run plp_main_step1.py via `--> python plp_main_stepX.py [target name] [-c -h keywords]`

*&dagger; Example command for the example data*:
```shell
~$ python plp_main_step1.py GJ\ 281 -c 4 
```
This step will generate a new recipe file `20170216.recipes` that separate the ABBA nodding into two pairs of AB and BA:
```python
OBJNAME, OBJTYPE, GROUP1, GROUP2, EXPTIME, RECIPE, OBSIDS, FRAMETYPES
# Available recipes : FLAT, THAR, SKY, A0V_AB, A0V_ONOFF, STELLAR_AB, STELLAR_ONOFF, EXTENDED_AB, EXTENDED_ONOFF
FLAT ON/OFF, FLAT, 11, 1, 30.000000, FLAT,11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30, OFF OFF OFF OFF OFF OFF OFF OFF OFF OFF ON ON ON ON ON ON ON ON ON ON
 , TAR, 59sky, 1, 300.000000, SKY,59, A
HR 1482, STD, 39, 1, 60.000000, A0V_AB, 39 40 41 42, A B B A
GJ 281, TAR, 60, 1, 600.000000, STELLAR_AB, 60 61, A B
GJ 281, TAR, 62, 1, 600.000000, STELLAR_AB, 62 63, B A
```
* There will be one `[NIGHT].recipes` file saved under `./recipes/[target_name]_recipes/`, and it will also be copied to `./recipes_logs`. 
  The `igrins plp` only monitors (reads) recipes under `./recipes_logs`. An example file structure:
```
plp(-master)
├── functions
│   ├── run_IGRINS.py
│   └── ...
│
├── recipe_logs
│   └── 20170216.recipes
│
├── recipes                  
│   ├── GJ281_recipes
│   │   ├── 20170216.recipes
│   │   └── 20170216.recipes.tmp   
│   ├── [Target]_recipes   
│   └── ... 
│
├── igr_pipe.py     
├── plp_main_step1.py
├── plp_main_step2.py
└── ...
```

* Note that depending on how you setup your environment, line 10 in `./functions/run_IGRINS.py` **might need to be changed**. The default is `fh.write("source activate igr-pipe\n")`, one might need to change `source` to `conda`, or change `igr-pipe` to their environment names.

  > One or more `.sh` bash scripts are generated under `./run_sh/`. Run them via `bash GJ281_run_igrins1.sh`.\
They contain commands to obtain the A and B separated spectra, e.g:
```bash
  python igr_pipe.py a0v-ab $UTDATE --basename-postfix="_A" --frac-slit="0,0.5"
  python igr_pipe.py a0v-ab $UTDATE --basename-postfix="_B" --frac-slit="0.5,1"
  python igr_pipe.py stellar-ab $UTDATE --basename-postfix="_A" --frac-slit="0,0.5"
  python igr_pipe.py stellar-ab $UTDATE --basename-postfix="_B" --frac-slit="0.5,1"
```

* Tip: Open/login with muti-terminal and run `bash GJ281_run_igrinsX.sh` one by one in each terminal to save time.


### Analyzing Multiple Observations Per Night (i.e., Multiple RVs)
If you have intensely observed a target such that you have enough data to constitute multiple observations per night, you can change the reduction accordingly. In the `./tmp` file with:
```python
GJ 281, TAR, 60, 54, 600.000000, STELLAR_AB,60 61 62 63 64 65 66 67 68 69 70 71, A B B A A B B A A B B A
```
you can modify it to:
```python
GJ 281, TAR, 60, 54, 600.000000, STELLAR_AB,60 61 62 63, A B B A
GJ 281, TAR, 64, 54, 600.000000, STELLAR_AB,64 65 66 67, A B B A
GJ 281, TAR, 68, 54, 600.000000, STELLAR_AB,68 69 70 71, A B B A
```
This would set up **IGRINS RV** to ultimately generate three RVs instead of one for the night.

  * Then you can run `plp_main_step1.py` normally.

> Note that for the **Muti-RVs Per Night** analysis, the files under `./recipes/[target_name]_recipes/` are important for **IGRINS RV** to track the split you did in the later processing.


## Getting the **input** data for **IGRINS RV**
After seeing the `IGRINS Pipeline Finished` print out, we are ready for `plp_main_step2.py` or `plp_main_step2_s.py`.
* `plp_main_step2.py` for **normal** usage
* `plp_main_step2_s.py` for **Muti-RVs Per Night** usage

*&dagger; Example command*:
```shell
~$ python plp_main_step2.py GJ\ 281
```

Run step2 like step1, and you will get final spectra under `./final_A_B_spec/[target_name]`.\
Now, move the entire folder `plp-master/final_A_B_spec/[target_name]` to `./igrins_rv-master/input/`. 

Voilà, you are almost ready to get RVs out from your data with **IGRINS RV**!\
Next you need to {doc}`Setup PrepData Files`,
and choose (see {doc}`Setup Stellar Templates`) the stellar template you would like to use.


