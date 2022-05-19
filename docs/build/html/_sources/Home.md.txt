Welcome to the **IGRINS RV** wiki page!\

If you have any question, idea or wish to report any bug, please let us know by either open an issue or contact us (asa.stahl@rice.edu or sytang@lowell.edu).


News:
* **2021/02/18: **IGRINS RV** v0.9.6-beta.3 public beta is ready for JOSS review!!\
(minor bugs related to example data fixed)**
* 2021/01/29: **IGRINS RV** v0.9.6-beta.1 public beta\
(update the ability to process IGRINS data at Gemini South)
* 2020/12/05: **IGRINS RV** v0.9.5-beta.1 public beta is under internal review and testing!!\
(update with modeling the flux suppression effect)
* 2020/08/04: **IGRINS RV** v0.9-beta.1 public beta is under internal review and testing!!
* 2020/06/23: **IGRINS RV** v0.85 is under internal review.. will come to public soon!!


***

For the sake of efficiency, **IGRINS RV** is broken up into several different steps. After [installing](https://github.com/shihyuntang/igrins_rv/wiki/Installation) all the required packages, the best way to learn how these steps work is to follow our [demo](https://github.com/shihyuntang/igrins_rv/wiki/Demo-with-Example-Data-(GJ281)) using the example data provided. The demo will demonstrate how one user might run the steps of **IGRINS RV** to go from raw data to final RV measurements. 

As you go through the demo, navigate to the other wiki pages as they are linked to get a more detailed description of how each step works and why you're running it. Once you've finished the demo and would like to begin analyzing your real science targets, these other pages will also explain what options you have for modifying the code's default inputs to best suit the needs of your science:

* [Installation](https://github.com/shihyuntang/igrins_rv/wiki/Installation)
* [Overview and Workflow](https://github.com/shihyuntang/igrins_rv/wiki/Overview-and-Workflow)
* Setup
   * [Setup: Input 1D spectra (plp)](https://github.com/shihyuntang/igrins_rv/wiki/Setup:-Input-1D-spectra-(plp))
   * [Setup: PrepData Files (Step 0)](https://github.com/shihyuntang/igrins_rv/wiki/Setup:-PrepData-Files)
   * [Setup: Stellar Templates](https://github.com/shihyuntang/igrins_rv/wiki/Setup:-Stellar-Templates)
* Main RV program
   * [Step 1: Telluric Modelling](https://github.com/shihyuntang/igrins_rv/wiki/Step-1:-Telluric-Modelling)
   * [Step 2: Initial Convergence](https://github.com/shihyuntang/igrins_rv/wiki/Step-2:-Initial-Convergence)
   * [Step 3: Full Analysis](https://github.com/shihyuntang/igrins_rv/wiki/Step-3:-Full-Analysis)
   * [Step 4: vsini Uncertainty Calculation](https://github.com/shihyuntang/igrins_rv/wiki/Step-4:-vsini-Uncertainty-Calculation)
* [FAQ](https://github.com/shihyuntang/igrins_rv/wiki/FAQ)

**Example data**:
1. Under `plp/indata/`. Contain raw spectral images for testing all `plp` related steps.
2. Under `igrins_rv/Example/`. Contain 1D reduced spectra for testing all **IGRINS RV**  related steps.
> If you are only testing the function of **IGRINS RV**, and not the `plp` package it makes use of, please skip all steps related to `plp`.


***
Please cite these papers if you use **IGRINS RV**:
* ``igrins plp v2.2.0`` https://zenodo.org/record/845059#.Xzlp5C0wJQI
* ``Telfit`` paper on AJ https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract
* **IGRINS RV** paper on AJ https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract
* **IGRINS RV** paper on JOSS https://joss.theoj.org/papers/10.21105/joss.03095#


