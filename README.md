[![status](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd/status.svg)](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd)


# IGRINS Radial Relocity Pipeline [IGRINS RV](https://github.com/shihyuntang/igrins_rv)


``IGRINS RV`` is a ``python`` open source pipeline for extracting radial velocity (RV) for the Immersion GRating INfrared Spectrometer (IGRINS) instrument. It's core technique use a modified forward modeling technique that utilize the telluric absorption lines as the common-path wavelength calibrator. With RV precision at K band around 25 m/s and H band 50 m/s, ``IGRINS RV`` had successfully retrieval of planet-induced RV signals for both HD 189733 and &tau; Boo A. Visit [Stahl et al. 2021 arXiv](https://ui.adsabs.harvard.edu/abs/2021arXiv210402082S/abstract) to see the published paper.

``IGRINS RV`` is now under intense development. \
If you have any question, idea or wish to report any bug, please let us know by either open an issue or contact us (asa.stahl@rice.edu or sytang@lowell.edu).
More on how to contribute can be found in the [FAQ](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-do-i-report-bugs) page.

News:
* **2021/02/18: `IGRINS RV` v0.9.6-beta.3 public beta is ready for JOSS review!!\
(minor bugs related to example data fixed)**
* 2021/01/29: `IGRINS RV` v0.9.6-beta.1 public beta\
(update the ability to process IGRINS data at Gemini South)
* 2020/12/05: `IGRINS RV` v0.9.5-beta.1 public beta is under internal review and testing!!\
(update with modeling the flux suppression effect)
* 2020/08/04: `IGRINS RV` v0.9-beta.1 public beta is under internal review and testing!!
* 2020/06/23: `IGRINS RV` v0.85 is under internal review.. will come to public soon!!

***
The best way to learn how IGRINS RV works is to first play with the example data following steps in [Demo Run With Example Data of GJ281](https://github.com/shihyuntang/igrins_rv/wiki/Demo-Run-With-Example-Data-of-GJ281). While you are waiting for the example code to finish, please read all detailed documentation on the [GitHub wiki page](https://github.com/shihyuntang/igrins_rv/wiki) so you know how to change the command to fit your OWN science targets:

Please Cite these papers if you use ``IGRINS RV``:
* `` igrins plp v2.2.0`` paper https://zenodo.org/record/845059#.Xzlp5C0wJQI
* ``Telfit`` paper https://ui.adsabs.harvard.edu/abs/2014AJ....148...53G/abstract
* ``IGRINS RV paper on arXiv`` https://ui.adsabs.harvard.edu/abs/2021arXiv210402082S/abstract

Acknowledgements:\
Many thanks to Dr. Gregory Mace for helping on the improvment of the user experence with `IGRINS RV`!
