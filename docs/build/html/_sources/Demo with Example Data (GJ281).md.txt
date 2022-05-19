# Demo with Example Data (GJ281)

***
This is a step by step walkthrough of **IGRINS RV** using the example data provided with the package installation.

For the sake of efficiency, **IGRINS RV** is broken up into several different steps. The purpose of this page is to provide an example of how one user might run these steps to go from raw data to final RV measurements. Input commands and their resultant printouts are provided for illustration. 

As you go through this demo, navigate to the other wiki pages as they are linked to get a more detailed description of how each step works and why you're running it. Once you've finished the demo and would like to begin analyzing your real science targets, these other pages will also provide explain what options you have for modifying the code's default inputs to best suit the needs of your science.

Make sure you have all the packages installed before you begin. If not, please visit {doc}`Installation` first.


**Example data**:
1. Under `plp/indata/`. Contain raw spectral images for testing all `plp` related steps.
2. Under `igrins_rv/Example/`. Contain 1D reduced spectra for testing all **IGRINS RV**  related steps.
> If you are only testing the function of **IGRINS RV**, and not the `plp` package it makes use of, please skip all steps related to `plp`.


## AN IMPORTANT NOTE
GJ281 is an RV standard (STD) in reality, but for the sake of demonstration, we treat GJ281 as a normal target (TAR) in this demo! If you want to see how the input commands would be different if we were treating GJ281 as an RV standard, that information can be found throughout the wiki pages on Step1 -- Step4.


If you are only testing the function of **IGRINS RV** -- and not that of the `plp` package it makes use of -- please go to [Demo for IGRINS RV](https://github.com/shihyuntang/igrins_rv/wiki/Demo-with-Example-Data-(GJ281)#demo-for-igrins-rv) section directly.


***

## Demo for getting the 1D spectra with `plp v2.2.0`
### Step 1: Making the bash run script

Make sure you are in the /plp(-master) directory.

Input:
```shell
~/plp$ python plp_main_step1.py GJ\ 281 -c 6
```

Output:
```
Step 1...
-------------------------------------
process dates:
 20170216 
-------------------------------------
Making A/B/AB recipes under plp read folder: ./recipe_log
Clearing dir ./recipe_log/ 
done
 processing: 20170216, 100.00% 1/1 -------------------------------------
Making .sh file
Step 2 Ended, please run --> bash GJ281_run_igrins.sh manually
```

### Step 2: run `plp v2.2.0` with the bash script
Now `cd` to `./run_sh`, and run
```shell
~/plp/run_sh $ bash GJ281_run_igrins1.sh
...
```
> While the code is running, visit [Setup: Input 1D spectra (plp)](https://github.com/shihyuntang/igrins_rv/wiki/Setup:-Input-1D-spectra-(plp)) to learn more about what is going on here. You might want to pay special attention to the section about making the **recipe** for `plp`.

### Step 3: 
Now `cd ..` to get back to the directory `./run_sh/`. Then run:
```shell
~/plp $ python plp_main_step2.py GJ\ 281
```
Output:
```
Step 2...
-------------------------------------
process dates:
 20170216 
20170216
-------------------------------------
Now overwrite the fits header from RAW to get the right JD time for each nodding...

Step 2 Ended. You are now ready for running IGRINS RV
Please copy the reduced 1D spectra target folder under "./final_A_B_spec" to the "input" folder in the "igrins_rv-master"
```
You can find the reduced 1D spectra under `./plp/final_A_B_spec/GJ281`. Copy this folder to `./igrins_rv/Input/`, where **IGRINS RV** will look for it. 
***
## Demo for IGRINS RV
### Input data
For this example, we have 6 nights of reduced data ready for you to use. Simply move `igrins_rv(-master)/Example/GJ281` to `igrins_rv/Input/`. If you completed the `plp` demo above, you should already done this step.

All users now move `igrins_rv/Example/Prepdata_A0_GJ281.txt` and `igrins_rv/Example/Prepdata_targ_GJ281.txt` to `igrins_rv/Input/Prepdata/`. 

> These files tell **IGRINS RV** what data are available for the science target and its associated A0 standards, as well as some information about the observations (e.g. time, barycentric velocity correction, observing conditions). Read [Setup: PrepData Files](https://github.com/shihyuntang/igrins_rv/wiki/Setup:-PrepData-Files) to learn how to make your own `PrepData Files` for your science targets (you may want to wait until you begin running Step 1, below). 

The folder structure should now look like:
```
igrins_rv(-master)
├── Input
│   ├── Prepdata
│   │   ├── Prepdata_A0_GJ281.txt
│   │   └── Prepdata_targ_GJ281.txt
│   ├── GJ281
│   │   ├── 20170106
│   │   ├── 20170128
│   │   └── ...
│   └── UseWv
│       └── ... 
├── Output
│   └── ... 
└── ... 
```

Activate the **igrins_rv** environment.
```shell
~/igrins_rv $ source activate igrins_rv
```
or 
```shell
~/igrins_rv $ conda activate igrins_rv
```

### Step 1: Telluric Modelling
```shell
(igrins_rv) ~/igrins_rv $ python main_step1.py GJ281 -HorK K -Wr 1 -c 6 -plot
```
While the code is running, visit [Step 1: Telluric Modelling](https://github.com/shihyuntang/igrins_rv/wiki/Step-1:-Telluric-Modelling) to learn about what you've just started. And if you're curious about how what you've done so far fits into the larger overall structure of the pipeline, check out [Overview and Workflow](https://github.com/shihyuntang/igrins_rv/wiki/Overview-and-Workflow) (you may want to skim this and return to it once you've finished the demo).

After a while, you should get this final info:

```
A0 Fitting Done!
A0 Fitting using TelFit finished, Duration: 0:45:41.780903
The synthetic telluric templates have been saved under ./Output/GJ281_K/A0Fits
If you chose to generate plots, they are saved under ./Output/GJ281_K/A0Fits/figs
####################################################################################
You can now run main_step2.py to produce RV and vsini initial guess(es)
####################################################################################
```

> Tip: All important printouts are also saved under the `.log` file, e.g., `/igrins_rv-master/Output/GJ281_K/A0Fits/GJ281_K_A0Fits.log`


### Step 2: Initial Convergence
#### Step 2 - Run No. 1

Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step2.py GJ281 -HorK K -Wr 1 -i 10 -v 5 -g 25 -t synthetic -temp 4000 -logg 4.5 -c 6 -plot -l_use 6
```
Output:
```
####################################################################################

---------------------------------------------------------------

Input Parameters:
    Tartget             =  GJ281
    Filter              =  K band 
    WaveLength file     =  WaveRegions_1 
    S/N cut             >  50 
    Order Use           =  Order 6 
    Initial vsini       =  10.0 km/s 
    vsini vary range    ±  5.0 km/s 
    RV initial guess    =  25.0 
    Stellar template use=  synthetic 
    syn template temp   =  4000 
    syn template logg   =  4.5 
    Threads use         =  6
    
Press [Y]es to continue, [N]o to quit...
 --> y
---------------------------------------------------------------
Running Step 2 for GJ281...
This Will Take a While..........
Writing output to ./Output/GJ281_K/Initguesser_results_1.csv
Analyze with 6 nights
Using K-band synthetic stellar template...
synthetic stellar template with T4000 logg4.5!!!!!

SUBMITTING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 462.79it/s]
PROCESSING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [02:14<00:00, 22.37s/it]
COLLECTING | : 100%|███████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 87078.98it/s]

```
While the code is running, visit [Step 2: Initial Convergence](https://github.com/shihyuntang/igrins_rv/wiki/Step-2:-Initial-Convergence) to learn more about what is going on here.

After a while, you should get this final printout:
```
--------!Initial Guess!--------
RV results:    mean= 20.1224 km/s, median= 20.1315 km/s, std= 0.0404 km/s
vsini results: mean= 5.0000 km/s, median= 5.0000 km/s, std= 0.0000 km/s
RV Initial Guess DONE... Duration: 0:01:42.228705
Output saved under ./Output/GJ281_K/Initguesser_results_1.csv
---------------------------------------------------------------
You can now try to get a better RV initial guess with by rerunning Step 2 with -gX set to the run number you just completed.
OR, you can go on to the full RV analysis in Step 3.
####################################################################################
```

#### Step 2 - Run No. 2 (This is optional and we do it here for demonstration purposes)
Let's rerun `Step 2: Initial Convergence` again, but with each observation's initial RV guess taken from `./Output/GJ281_K/Initguesser_results_1.csv`.

> The command will be different from the previous one in that:
> 1. -i 10 --> -i 5
> 2. -g 25 --> -gX 1

> Tip: add `-sk_check` to your command to skip the [Y]es and [N]o Input Parameters checking process. This is useful when running a sequence of commands in a script.

Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step2.py GJ281 -HorK K -Wr 1 -i 5 -v 5 -gX 1 -t synthetic -temp 4000 -logg 4.5 -c 6 -plot -l_use 6 -sk_check
```

Output:
```
####################################################################################

---------------------------------------------------------------
Running Step 2 for GJ281...
This Will Take a While..........
Writing output to ./Output/GJ281_K/Initguesser_results_2.csv
Analyze with 6 nights
Using K-band synthetic stellar template...
synthetic stellar template with T4000 logg4.5!!!!!

SUBMITTING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 22.99it/s]
PROCESSING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [02:14<00:00, 22.47s/it]
COLLECTING | : 100%|██████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 133152.51it/s]

--------!Initial Guess!--------
RV results:    mean= 20.1373 km/s, median= 20.1465 km/s, std= 0.0390 km/s
vsini results: mean= 3.1850 km/s, median= 3.2428 km/s, std= 0.2046 km/s
RV Initial Guess DONE... Duration: 0:01:33.994598
Output saved under ./Output/GJ281_K/Initguesser_results_2.csv
---------------------------------------------------------------
You can now try to get a better RV initial guess with by rerunning Step 2 with -gX set to the run number you just completed.
OR, you can go on to the full RV analysis in Step 3.
####################################################################################
```


### Step 3: Full Analysis
#### Step 3 - Run No. 1 (getting precise absolute RVs)
The input parameters for Step 3 are based on the results from step 2. \
For the starting vsini guess (`-i`), use the printout from Step2: `vsini results: mean= 3.1850 km/s` from step2 (command: `-i 3.1850`).\
For the starting guesses for the RVs, use the `./Output/GJ281_K/Initguesser_results_2.csv` from Step 2 (command: `-gS init -gX 2`).

Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step3.py GJ281 -mode TAR -HorK K -Wr 1 -nAB 1 -i 3.1850 -v 5 -gS init -gX 2 -t synthetic -temp 4000 -logg 4.5 -c 6 -plot -abs_out abs -sk_check
```
Output:
```
####################################################################################

---------------------------------------------------------------
Running Step 3 for GJ281...
This will take a while..........
Writing output to ./Output/GJ281_K/RV_results_1
Analyze with 6 nights
Using K-band synthetic stellar template...
synthetic stellar template with T4000 logg4.5!!!!!

Working on order 3 (01/04)
SUBMITTING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 19.18it/s]
PROCESSING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [05:25<00:00, 54.30s/it]
COLLECTING | : 100%|██████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 139037.70it/s]
Working on order 4 (02/04)
SUBMITTING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 659.03it/s]
PROCESSING | : 100%|██████████████████████████████████████████████████████████████████| 6/6 [07:11<00:00, 71.91s/it]
COLLECTING | : 100%|██████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 142987.64it/s]
Working on order 5 (03/04)
SUBMITTING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 726.92it/s]
PROCESSING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [12:50<00:00, 128.48s/it]
COLLECTING | : 100%|██████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 142987.64it/s]
Working on order 6 (04/04)
SUBMITTING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 733.59it/s]
PROCESSING | : 100%|█████████████████████████████████████████████████████████████████| 6/6 [17:32<00:00, 175.37s/it]
COLLECTING | : 100%|██████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 139810.13it/s]
...
```
While the code is running, visit [Step 3: Full Analysis](https://github.com/shihyuntang/igrins_rv/wiki/Step-3:-Full-Analysis) to learn more about what is going on here.

After a while, you should get these final info.
```
RV mean = 19.9862 km/s, std = 0.0364 km/s
vsini mean = 3.3244 km/s, std = 0.4313 km/s

**********************************************************************************
WARNING!! you got warning message during this run. Please check the log file under:
          ./Output/GJ281_K/GJ281_K_A0Fits.log
**********************************************************************************

Whole process DONE!!!!!!, Duration: 0:32:15.134722
Output saved under ./Output/GJ281_K/RV_results_1
The final RV estimates you are looking for are in the RVresultsSummary files!
```
You can visually check how the fitting went by looking at the plots saved under `./Output/GJ281_K/figs/main_step3_K_1/`


#### Step 3 - Run No. 2 (getting relative RVs that are more precise than the absolute RVs)
> The command differences from the previous one are:
> 1. -i 3.1850 --> 3.3244 (start the mean vsini output by the previous run)
> 2. -gS init -gX 2 --> -gS rvre -gX 1 (start from the output RVs from the previous run, e.g., `./Output/GJ281_K/RVresultsSummary_1.csv`)
> 3. -abs_out abs --> -abs_out rel (get relative, not absolute, RVs)

Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step3.py GJ281 -mode TAR -HorK K -Wr 1 -nAB 1 -i 3.3244 -v 0 -gS rvre -gX 1 -t synthetic -temp 4000 -logg 4.5 -c 6 -plot -abs_out rel -sk_check
```
Output:
```
...
RV mean = -0.0071 km/s, std = 0.0291 km/s
vsini mean = 3.3244 km/s, std = 0.0000 km/s

Whole process DONE!!!!!!, Duration: 0:32:08.984956
Output saved under ./Output/GJ281_K/RV_results_2
The final RV estimates you are looking for are in the RVresultsSummary files!
####################################################################################
```

#### Step 3 - Run No. 3 (getting the relative RVs with vsini fixed at its best fit plus 1 sigma)
> Difference in the terminal command from the previous one is:
> 1. -i 3.3244 --> 3.7557 (Use RV_results_1's mean vsini plus the standard deviation of that estimate, 3.3244 + 0.4313 = 3.7557)

Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step3.py GJ281 -mode TAR -HorK K -Wr 1 -nAB 1 -i 3.7557 -v 0 -gS rvre -gX 1 -t synthetic -temp 4000 -logg 4.5 -c 6 -plot -abs_out rel -sk_check
```
Output:
```
...
RV mean = -0.0075 km/s, std = 0.0316 km/s
vsini mean = 3.7557 km/s, std = 0.0000 km/s

Whole process DONE!!!!!!, Duration: 0:30:54.705114
Output saved under ./Output/GJ281_K/RV_results_3
The final RV estimates you are looking for are in the RVresultsSummary files!
####################################################################################
```

Why would you purposefully try to calculate RVs with vsini offset from your best estimate? Read about [Step 4](https://github.com/shihyuntang/igrins_rv/wiki/Step-4:-vsini-Uncertainty-Calculation).

### Step 4: Calculating the impact of vsini uncertainty on your RVs
Input:
```shell
(igrins_rv) ~/igrins_rv $ python main_step4.py GJ281 -run1 2 -run2 3 -HorK K
```
Output:
```
####################################################################################

---------------------------------------------------------------
Input Parameters:
    Tartget             =  GJ281
    Filter              =  K band 
    Run 1               =  RV_results_2    <------- Should be run where vsini = best guess
    Run 2               =  RV_results_3    <------- Should be run where vsini = best guess +- sigma
    
Press [Y]es to continue, [N]o to quite...
 --> y
---------------------------------------------------------------
Running Step 4 for GJ281...
Writing output to ./Output/GJ281_K/RV_results_2_3_combined
Combined RV results: mean=-0.0071 km/s, std=0.0291 km/s


Whole process DONE!!!!!!, Duration: 0:00:02.393307
Output saved under ./Output/GJ281_K/RV_results_2_3_combined
####################################################################################
```

Voilà, you are ready to run **IGRINS RV** with your science data!
> It might be useful to take a look at [Overview and Workflow](https://github.com/shihyuntang/igrins_rv/wiki/Overview-and-Workflow) again before you try to run the pipeline with your science data.


