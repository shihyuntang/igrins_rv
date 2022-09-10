# FAQ

**Still under construction...**

<!-- ## Installation 
* [How do I install `gfortran` for `Telfit`?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-do-i-install-gfortran-for-telfit)
* [Why do I see error message with `ifort` while installing `Telfit`?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-why-do-i-see-an-error-message-with-ifort-while-installing-telfit)

## Running **IGRINS RV** 
* [What is the typical run time of **IGRINS RV**?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-what-is-the-typical-run-time-of-igrins-rv)
* [How much RAM do I need to run **IGRINS RV**?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-much-ram-do-i-need-to-run-igrins-rv)
* [Can **IGRINS RV** do any better in RV precision?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-can-igrins-rv-do-any-better-in-rv-precision)

## Other
* [How do I report bugs?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ/_edit#q-how-do-i-report-bugs) 
* [How do I help make **IGRINS RV** better?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ/_edit#q-how-do-i-help-make-igrins-rv-better) 
* [How do I properly cite **IGRINS RV**?](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-do-i-properly-cite-igrins-rv)

*** -->

## Q: How do I install `gfortran` for `Telfit`?
**Answer:** 
* On Mac: Use [Homebrew](https://brew.sh/index_ja) with `brew install gcc`
* On Linux: With `sudo apt install gfortran`

If successful, try `gfortran --version` to see if you can see the versions number.
If so, you are good!

--

## Q: How can I tell `Telfit` to use `ifort` compiler?
**Answer:** Change the order of compilers in `/Telluric-Fitter(-master)/setup.py` line [117](https://github.com/kgullikson88/Telluric-Fitter/blob/7ae98db278525e157d2d0abaf4697e2fe778d6bc/setup.py#L117) to 122
from 
```
    compilers = ["gfortran",
                 "ifort",
                 "g95"]
    comp_strs = ["GNU",
                 "INTEL",
                 "G95"]
```
to
```
    compilers = ["ifort",
                 "gfortran",
                 "g95"]
    comp_strs = ["INTEL",
                 "GNU",
                 "G95"]
```

--
## Q: Why do I see an error message with `ifort` while installing `Telfit`?
**Answer:** Well...lucky you. This issue occurs if you have a new version of `ifort` but running with old versions of lblrtm/lnfl. Updating the lblrtm/lnfl (you might as well also update the aer_line_files) to the latest version will fix this. See {doc}`Use latest lblrtm` for detail instructions.

--
## Q: What is the typical run time of **IGRINS RV**?
**Answer:** This really depends on the setup of your machine. With an Intel Core i9-9980XE CPU (18 cores/36 threads) running GJ281 K band on 64 nights and 4 orders, it took about 6.6 hours for Step 1 (Telluric Modelling) and about 1.2 hours for Step 3 (Analysis) to finish. The runtime will be about 1.5 times longer in the H band because there are more orders to process.

--
## Q: How many memory (RAM) do I need to run **IGRINS RV**?
**Answer:** The most RAM consuming task is in step 3, and it also depends on how many threads you use. The more threads you use, the more RAM you need. Running GJ281 K band on 64 nights with 36 threads, the maximum RAM used are about 8 GB.

--
## Q: Can **IGRINS RV** do any better in RV precision?
**Answer:** **IGRINS RV** uses telluric (atmospheric) absorption lines as a common-path wavelength calibrator, so its precision is limited by the internal RV stability of Earth's atmosphere. This has been estimated to be ~10-20 m/s, so it is unlikely **IGRINS RV** will achieve better  than the current ~25m/s precision. If you are obtaining worse precisions than this for your science targets and want to try to push toward better precisions, improving the accuracy of your stellar templates or increasing the signal quality of your data could help (note the precision achievable by **IGRINS RV** also depends on the vsin(i) of the target in question).

--
## Q: What happened to absolute vs relative RV modes?
**Answer:** All RVs are now absolute, not relative. See Section "What's new in v1.5.1?".

--
## Q: How do I report bugs?
**Answer:** Please open an issue at [https://github.com/shihyuntang/igrins_rv/issues](https://github.com/shihyuntang/igrins_rv/issues). When doing so, please let us know:
1. your operating system, and the OS version.
2. detailed bug description
2. detailed steps on how to reproduce the bug reported.

--
## Q: How do I help make **IGRINS RV** better?
**Answer:** Please open an issue at [https://github.com/shihyuntang/igrins_rv/issues](https://github.com/shihyuntang/igrins_rv/issues) and describe your ideas to us. Currently, any help on making **IGRINS RV** run on Windows is welcome. Go to [this issue](https://github.com/shihyuntang/igrins_rv/issues/7) and join the discussion with us.

--
## Q: How do I properly cite **IGRINS RV**?
**Answer:** Please cite both the AJ and the JOSS papers for **IGRINS RV**:
* AJ 161:283: [Stahl et al. 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract)
* JOSS 6:62: [Tang et al. 2021](https://joss.theoj.org/papers/10.21105/joss.03095#)

And cite:
* ``igrins plp v2.2.0`` [https://zenodo.org/record/845059#.Xzlp5C0wJQI](https://zenodo.org/record/845059#.Xzlp5C0wJQI)
for 1d spectra reduction, and
* ``Telfit`` paper on AJ [https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract](https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract)
for telluric line fitting.

BibTex:
```bibtex
@ARTICLE{2021AJ....161..283S,
       author = {{Stahl}, Asa G. and {Tang}, Shih-Yun and {Johns-Krull}, Christopher M. and {Prato}, L. and {Llama}, Joe and {Mace}, Gregory N. and {Joon Lee}, Jae and {Oh}, Heeyoung and {Luna}, Jessica and {Jaffe}, Daniel T.},
        title = "{IGRINS RV: A Precision Radial Velocity Pipeline for IGRINS Using Modified Forward Modeling in the Near-infrared}",
      journal = {\aj},
     keywords = {Exoplanet detection methods, Exoplanet formation, Radial velocity, Near infrared astronomy, Open source software, Astronomy software, Starspots, Young stellar objects, 489, 492, 1332, 1093, 1866, 1855, 1572, 1834, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2021,
        month = jun,
       volume = {161},
       number = {6},
          eid = {283},
        pages = {283},
          doi = {10.3847/1538-3881/abf5e7},
archivePrefix = {arXiv},
       eprint = {2104.02082},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@article{Tang2021,
  doi = {10.21105/joss.03095},
  url = {https://doi.org/10.21105/joss.03095},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {62},
  pages = {3095},
  author = {Shih-Yun Tang and Asa G. Stahl and Christopher M. Johns-Krull and L. Prato and Joe Llama},
  title = {IGRINS RV: A Python Package for Precision Radial Velocities with Near-Infrared Spectra},
  journal = {Journal of Open Source Software}
}

@software{jae_joon_lee_2017_845059,
  author       = {Jae-Joon Lee and
                  Kevin Gullikson and
                  Kyle Kaplan},
  title        = {igrins/plp 2.2.0},
  month        = aug,
  year         = 2017,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.845059},
  url          = {https://doi.org/10.5281/zenodo.845059}
}

@ARTICLE{2014AJ....148...53G,
       author = {{Gullikson}, Kevin and {Dodson-Robinson}, Sarah and {Kraus}, Adam},
        title = "{Correcting for Telluric Absorption: Methods, Case Studies, and Release of the TelFit Code}",
      journal = {\aj},
     keywords = {atmospheric effects, instrumentation: spectrographs, techniques: spectroscopic, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2014,
        month = sep,
       volume = {148},
       number = {3},
          eid = {53},
        pages = {53},
          doi = {10.1088/0004-6256/148/3/53},
archivePrefix = {arXiv},
       eprint = {1406.6059},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
