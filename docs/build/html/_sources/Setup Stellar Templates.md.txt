# Setup: Stellar Templates

***IGRINS RV** comes with a handful of synthetic stellar templates, but the user may wish to supply their own.*\
This page describes those that come with the code as well as how to use your own templates.
Details on how the supplied templates were constructed, as well as how sub-optimal choices of template Teff and log(g) may affect RVs, are given in the paper.
***


## SYNTHMAG
The following templates are from the `SYNTHMAG IDL code` (Kochukhov 2007) with atmospheres profile from the PHOENIX NextGen model (Hauschildt et al. 1999) and the stellar lines from the VALD database (Ryabchikova et al. 2015). 
Call these templates via for example:
```shell 
(igrins_rv) ~$ python main_stepX.py [target name] -HorK H -t synthetic -temp 5000 -logg 4.5
```
|  Band  |   Teff (K)    |  log(g) (cgs) |
| :----: | :-----------: | ------------- |
| H & K | 3000 | 4.0 |
| K     | 3500 | 4.0 |
| H & K | 3500 | 4.5 |
| K     | 3600 | 4.5 |
| K     | 3700 | 4.5 |
| K     | 3900 | 4.5 |
| K     | 4000 | 3.5 |
| H & K | 4000 | 4.5 |
| K     | 4200 | 5.0 |
| K     | 4400 | 4.0 & 4.5 |
| K     | 4600 | 4.5 |
| K     | 4800 | 5.0 |
| H & K | 5000 | 4.5 |
| H     | 5800 | 4.5 |
| H     | 6200 | 3.5 & 4.5 |
| H & K | 6400 | 4.5 |
| H     | 6600 | 4.5 |

> These can generally be applied to stars of the same log(g) and a Teff within a few hundred K of the template. 

If you used any of the above templates, please cite these papers:
* [Kochukhov 2007 (Spectrum synthesis for magnetic, chemically stratified stellar atmospheres)](https://ui.adsabs.harvard.edu/abs/2007pms..conf..109K/abstract)
* [Hauschildt et al. 1999 (The NEXTGEN Model Atmosphere Grid. II. Spherically Symmetric Model Atmospheres for Giant Stars with Effective Temperatures between 3000 and 6800 K)](https://ui.adsabs.harvard.edu/abs/1999ApJ...525..871H/abstract)
* [Ryabchikova et al. 2015 (A major upgrade of the VALD database)](https://iopscience.iop.org/article/10.1088/0031-8949/90/5/054005)

## PHOENIX template from the Gottingen spectral library
A few sets of [Gottingen spectral templates](http://phoenix.astro.physik.uni-goettingen.de/?page_id=15) are available.
Call these templates with (for example):
```shell 
(igrins_rv) ~$ python main_stepX.py [target name] -HorK K -t PHOENIX -temp 3700 -logg 4.5
```
|  Band  |   Teff (K)    |  log(g) (cgs) |
| :----: | :-----------: | ------------- |
| H & K | 3000 | 4.0 & 4.5 |
| K     | 3700 | 4.5 |
| K     | 4000 | 4.5 |
| H     | 6200 | 4.5 |


## User supply templates
If you would like to supply your own template, the filename must follow the convention 
`user_T[temperature]_logg[logg]_[band]band.txt`, the file must be placed in ``./Engine/user_templates/`` , and the file itself must follow the format:

```python
 #      wave           flux
       19380.000      0.99825212
       19380.036      0.99824445
       19380.073      0.99823262
```

with the same first line listed here, followed by the (vacuum) wavelength in angstroms and the normalized flux. There should be no trailing blank lines at the end of the file.

```diff
---All stellar templates need to be **flattened** before use!!---
```
