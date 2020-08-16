# IGRINS Radial Relocity Pipeline [IGRINS RV](https://github.com/shihyuntang/igrins_rv)

**This package is NOT yet ready for public use.**

``IGRINS RV`` if a ``python`` open source pipeline for extracting radial velocity (RV) for the IGRINS instrument. It's core technique use a modified forward modeling technique that utilize the telluric absorption lines as the common-path wavelength calibrator. With RV precision at K band around 25 m/s and H band 50 m/s, ``IGRINS RV`` had successfully retrieval of planet-induced RV signals for both HD 189733 and &tau; Boo A. Visit [arXiv]() to see the published paper.

``IGRINS RV`` is now under intense development. \
If you wish to join the internal review and testing, please contact us (asa.stahl@rice.edu or sytang@lowell.edu).

News:\
**2020/08/04: igrins_rv v0.9-beta.1 public beta is under internal review and testing!!**\
2020/06/23: igrins_rv v0.85 is under internal review.. will come to public soon!!

***
Please visit the [GitHub wiki page](https://github.com/shihyuntang/igrins_rv/wiki) for a more detailed documentation and tutorials.
Here only contains a short XXX

Setting up and running igrins_rv is easy if you have Conda installed.
NOTE: ``IGRINS RV`` has only be tested on Linux machine...

## Packages Installation (part 1)
**Installation with conda**\
Single command to setup environment with all pkg needed (using the with the environment.yml file):
```
conda env create
```

This will create an environment called ``igrins_rv``. You can use
```
conda info --envs
```
to check all available environments.

Use
```
conda (source) activate igrins_rv
```
to enter the environment.\
If you do it correctly, your command line will now looks like:
```
(igrins_rv) -->
```
**If you do not have conda**, the basic requirement for running ``IGRINS RV``\
is python.3.7, and the following packages/versions:
>    - python=3.7
>    - matplotlib=3.1.3
>    - numpy=1.18.1
>    - pandas=1.0.3
>    - scipy=1.4.1
>    - astropy=4.0
>    - multiprocess
>    - cython=0.29.15
>    - requests=2.23.0
>    - nlopt=2.6.1
>    - pip:
>      - pysynphot=0.9.14

## Packages installation (part 2)
**``IGRINS RV`` comes with the files for Telfit, but they need to be installed.**

Note that the original Telfit pkg is available on https://github.com/kgullikson88/Telluric-Fitter, by @kgullikson88.\
But the files that come with ``IGRINS RV`` are **modified**, such that one line in ``./src/TelluricFitter.py`` is different:
```
line 672: resolution=None, vac2air=self.air_wave)
           --> resolution=None, vac2air=False)
```
To install Telfit, First, check your mechine's ``gcc`` version by
```
gcc --version
```
**Telfit will not work with version later than v9.X.X**\
Not sure for v8.X.X, but v5.5.0 & v7.5.0 work fine.

Also, Telfit only recognize "gcc", so please make sure the new installed gcc **is callable** by
```
gcc
```
If you intend to use more the 6 threads to run this pipeline, we recommended you go into the ``setup.py`` file, and\
modify ``num_rundirs = `` to a number that is **THREE** times the threads you intend to use.

Now, enter the ``igrins_rv`` environment (within which Telfit must be installed)\
and cd to Telluric-Fitter-master, then run:
```
(igrins_rv) --> python setup.py build
(igrins_rv) --> python setup.py install
```

## Running igrins_rv

For a detailed understanding of how ``IGRINS RV`` works and how it should be run, it is HIGHLY recommended that the user reads the associated [paper]().\
``IGRINS RV`` is divided into four main steps: Setup, Telluric Modelling, Initial Convergence, and Analysis.\
Each is provided in the package as a separate file and # can be run from the command line with keywords specifying\
all relevant information.

**All steps are run as**
```
(igrins_rv) --> python main_stepX.py [target name] [keywords]
```
where X is the step number (0-4) and information on the keywords for the step can be found with
```
(igrins_rv) --> python main_stepX.py [target name] -h
```
(The -h command here means --help.)

### Step 0 - Setup: Collects and organizes the available data for the target spectra and their associated telluric (A0) standards. 

**If the user have** a copy of the IGRINS Observation MASTERLOG, Step 0 will automatically consult a master observing log file\
to retrieve the relevant information on filenames, coordinates, observing times, and conditions.
```
(igrins_rv) --> python main_step0.py [target name] -
```
**If the user do not have** a copy of the IGRINS Observation MASTERLOG, the output files it generates must be compiled by the user themselves. 


### Step 1 - Telluric Modelling: \
Defines the wavelength regions to be analyzed; generates a synthetic, high-resolution telluric template for use in later model fits on a night by night basis. 

### Step 2 - Initial Convergence: \
Required if the average RV of the target star is unknown to $>$ 5 \kms precision. Performs an abbreviated analysis of the target star observations in order to converge to coarsely accurate RVs, which will be used as starting points for the more precise analysis in the next step; simultaneously does the same for target star's \vsini. Only the single most precise echelle region is used, and all separate exposures for a given observation are combined into one higher S/N spectrum before fitting occurs. 

### Step 3 - Analysis: \
Performs a full analysis of each target star observation to produce accurate and precise RVs. All the wavelength regions defined in Step 1 are used, and the code analyzes each exposure that is part of a given observation separately (this allows estimates of the RV uncertainties - see Section \ref{ssec:err}). In other words, if a star is observed with the ABBA beam pattern, Step 2 analyzes one spectrum made from the combination of the four exposures, while Step 3 analyzes each A and B separately. 

Unless the target \vsini is already known to high accuracy, an initial run of Step 3 in which \vsini is allowed to vary is required. This will provide an estimate of \vsini, which will then be plugged in as a fixed value in the second run of Step 3. 

Even if \vsini is already known, it is recommended to run Step 3 twice, as it allows the user to confirm that RVs between consecutive runs have converged within the calculated uncertainties. Additional runs may be necessary depending on the quality of the observed spectra and the \vsini of the target star. 

All above three steps can be ran by:
```
(igrins_rv) --> python main_step1.py [-h][keywords]
```
The `-h` command will show the basic information of main_step1.py, and tell you what parameters you can put in.


