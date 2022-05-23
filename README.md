[![status](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd/status.svg)](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd)
[![DOI](https://zenodo.org/badge/266670787.svg)](https://zenodo.org/badge/latestdoi/266670787)

# [IGRINS RV](https://github.com/shihyuntang/igrins_rv): A Radial Velocity Pipeline for IGRINS


``IGRINS RV`` is a ``python`` open source pipeline for extracting radial velocities (RVs) from spectra taken with the Immersion GRating INfrared Spectrometer (IGRINS). It uses a modified forward modeling technique that leverages telluric absorption lines as a common-path wavelength calibrator. ``IGRINS RV`` achieves an RV precision in the H and K bands of around 25-30 m/s for narrow-line stars, and it has successfully recovered the planet-induced RV signals of both HD 189733 and &tau; Boo A. Visit [Stahl et al. 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract) to see the published paper.

If you have any questions, suggestions, or wish to report a bug, please let us know by either opening an issue or contacting us (asa.stahl@rice.edu or sytang@lowell.edu).
More on how to contribute can be found in the [FAQ](https://igrins-rv.readthedocs.io/en/latest/FAQ.html#q-how-do-i-report-bugs) page.

## Installation
A detailed installation guide can be found [here](https://igrins-rv.readthedocs.io/en/latest/Installation.html).

## How to run
The best way to learn how IGRINS RV works is to first play with the example data, following steps in [Demo Run With Example Data of GJ281](https://igrins-rv.readthedocs.io/en/latest/Demo%20with%20Example%20Data%20%28GJ281%29.html). While you are waiting for the example code to finish, read the detailed documentation on [https://igrins-rv.readthedocs.io](https://igrins-rv.readthedocs.io) so you know how to alter the example commands to fit your OWN science targets.

***
If you use ``IGRINS RV`` in your work, please go to the [FAQ page](https://igrins-rv.readthedocs.io/en/latest/FAQ.html) to see how to properly cite ``IGRINS RV``.

Acknowledgements:\
Many thanks to Dr. Gregory Mace for helping improve the user experence with `IGRINS RV`!
