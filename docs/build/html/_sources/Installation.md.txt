# Installation

*Setting up and running **IGRINS RV** is easy if you have **Anaconda(Conda)** installed!\
NOTE:The current versions (v1.0.0) of **IGRINS RV** has only been tested on Linux and MacOS...*
***

## Package installation (part 1)

First, download/clone **igrins_rv** from github to the directory of your choice. Then `cd` to that directory in terminal. 

### Package installation with conda
You can use a single command to setup an environment with all needed packages (except the ``Telfit`` package):
```
conda env create
```
(this reads the **environment.yml** file and installs all the packages listed inside)

The above command will create an environment called ``igrins_rv``. You can use
```
conda info --envs
```
to check all available environments.

Use
```
conda activate igrins_rv
```
or
```
source activate igrins_rv
```
to enter the environment.

If you are successful, your command line will look something like
```
(igrins_rv) ~$
```

### Manual package installation
If you do not have **Conda**, the basic requirement for running **IGRINS RV** 
is with python.3.7 or later, and the following packages/versions:
* python=3.7
* matplotlib=3.1.3
* numpy=1.18.1
* pandas=1.0.3
* scipy=1.4.1
* astropy=4.0
* multiprocess
* cython=0.29.15
* requests=2.23.0
* nlopt=2.6.1
* via pip: pysynphot=0.9.14
* via pip: more_itertools
* via pip: pqdm


***

## Packages installation (part 2) - ``Telfit``

The most up to date version of `Telfit v1.4.0` is still under the beta test stage of pip installation; thus, 
```diff
---PLEASE install it via source---
```
Go to  [`Telfit`](https://github.com/kgullikson88/Telluric-Fitter), download the `master branch`, and install it from source (If you've never installed a pkg from source, no worries, we will walk you through it). 

**To install ``Telfit`` from source:** Enter the `igrins_rv` environment (within which Telfit must be installed) and `cd` into `Telluric-Fitter(-master)`, then run
```
(igrins_rv) ~$ python setup.py build
(igrins_rv) ~$ python setup.py install
```

After that, you can use
```
conda list
```
to check if `Telfit` installed successfully.


***
## Packages installation (part 3) - ``igrins plp``

> If you are only testing the functionality of **IGRINS RV**, you can skip steps related to `plp`, and go to the {ref}`Demo for IGRINS RV` section under the {doc}`Demo with Example Data (GJ281)`. 

To enable error estimation, **IGRINS RV** cannot get RVs from the 1D spectra reduced by ``igrins plp v2.2.0`` from combined ABs/ABBAs/etc. Users will always need to re-run ``igrins plp`` to get separated A and B frame 1D spectra.
Download or clone ``igrins plp v2.2.0`` at [https://github.com/shihyuntang/plp](https://github.com/shihyuntang/plp).
See {doc}`Setup Input 1D spectra (plp)` for more details.

> If you already have ``igrins plp v2.2.0`` on your machine, you can just download `plp_main_step1.py`, `plp_main_step2.py`, and `plp_main_step2_s.py`. Put these three `.py` files at the same level as the `igr_pipe.py`.

