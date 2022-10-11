% igrins_rv documentation master file, created by
% sphinx-quickstart on Thu May 19 13:50:00 2022.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.
<!-- 
```{include} ../../README.md
``` -->

[![status](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd/status.svg)](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd)
[![DOI](https://zenodo.org/badge/266670787.svg)](https://zenodo.org/badge/latestdoi/266670787)

# IGRINS RV: A Radial Velocity Pipeline for IGRINS

**IGRINS RV** is a ``python`` open source pipeline for extracting radial velocities (RVs) from spectra taken with the Immersion GRating INfrared Spectrometer (IGRINS). It uses a modified forward modeling technique that leverages telluric absorption lines as a common-path wavelength calibrator. **IGRINS RV** achieves an RV precision in the H and K bands of around 25-30 m/s for narrow-line stars, and it has successfully recovered the planet-induced RV signals of both HD 189733 and &tau; Boo A. Visit [Stahl et al. 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract) to see the published paper.

If you have any questions, suggestions, or wish to report a bug, please let us know by either opening an issue or contacting us: Asa G. Stahl (asa.stahl@rice.edu) or Shih-Yun Tang (sytang@lowell.edu).
More on how to contribute can be found in {ref}`Q: How do I help make **IGRINS RV** better?` page.

> Starting from v1.5.1, the example data & synthetic stellar templates are not shipped with **IGRINS RV** as to reduce the repo size. Synthetic stellar templates are now stored on google drive [link](https://drive.google.com/drive/folders/1WRiQ3PKCbhueQi6htd0zusq_1ieKXgHP?usp=sharing).

**News:**
* **2022/10/11: IGRINS RV v1.5.1!!**
    * **Who might be effected?**
      * People with only a few observations (<~5) that use K band.
      * People who need more *precise* ABSOLUTE RV.
    * **What's new?:**
      1. Significantly reduced the K band RV shift between orders.
      2. Removed choice between *relative* and *absolute* RV modes (for both H and K band). Now only absolute RVs are calculated, but an order-to-order RV correction is still applied, and the uncertainties provided are updated to take into account this correction. (for more info, see {doc}`What to Know About v1.5.1 Update`)
      3. Dropped the usage of Order 3 for K band RV calculation, as well as Order 4 during the period when IGRINS was slightly defocused as the result of a mounting issue.

***

For the sake of efficiency, **IGRINS RV** is broken up into several different steps. After installed (see {ref}`Installation`) all the required packages, the best way to learn how these steps work is to follow our {doc}`Demo with Example Data (GJ281)` using the example data provided. The demo will demonstrate how one user might run the steps of **IGRINS RV** to go from raw data to final RV measurements. 

```{toctree}
:caption: 'Quick Start'
:maxdepth: 1

Installation
Demo with Example Data (GJ281)
Overview and Workflow
```

As you go through the demo, navigate to the other document pages as they are linked to get a more detailed description of how each step works and why you're running it. Once you've finished the demo and would like to begin analyzing your real science targets, these other pages will also explain what options you have for modifying the code's default inputs to best suit the needs of your science:

```{toctree}
:caption: 'Setup'
:maxdepth: 1

Setup Input 1D spectra (plp)
Setup PrepData Files
Setup Stellar Templates
```

```{toctree}
:caption: 'Main RV program'
:maxdepth: 1

Step 1 Telluric Modelling
Step 2 Initial Convergence
Step 3 Full Analysis
Step 4 vsini Uncertainty Calculation
```

```{toctree}
:caption: 'Other'
:maxdepth: 1

FAQ
Use latest lblrtm
Version history
What to Know About v1.5.1 Update
```


***
Please cite these papers if you use **IGRINS RV**:
* **igrins plp v2.2.0**: [https://zenodo.org/record/845059#.Xzlp5C0wJQI](https://zenodo.org/record/845059#.Xzlp5C0wJQI)
* **Telfit** paper on AJ: [https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract](https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract)
* **IGRINS RV** paper on AJ: [https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract](https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract)
* **IGRINS RV** paper on JOSS: [https://joss.theoj.org/papers/10.21105/joss.03095#](https://joss.theoj.org/papers/10.21105/joss.03095#)


<!-- ```{toctree}
:caption: 'Contents:'
:maxdepth: 2
:hidden:

Home
Installation
Overview-and-Workflow
FAQ
``` -->

<!-- # Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search` -->
