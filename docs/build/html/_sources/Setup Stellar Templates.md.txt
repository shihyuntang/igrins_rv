# Setup: Stellar Templates

***IGRINS RV** comes with a handful of synthetic stellar templates, but the user may wish to supply their own.*\
This page describes those that come with the code as well as how to use your own templates.
Details on how the supplied templates were constructed, as well as how sub-optimal choices of template Teff and log(g) may affect RVs, are given in the paper.
***

> Starting from v1.5.1, the synthetic stellar templates are not shipped with **IGRINS RV** as to reduce the repo size. Synthetic stellar templates are now stored on google drive [link](https://drive.google.com/drive/folders/1WRiQ3PKCbhueQi6htd0zusq_1ieKXgHP?usp=sharing) under `./syn_template/`. 
> 
> Also, Starting from v1.5.1, IGRINS RV can read in syn spectra as .gz compress format, so NOT NEED TO **unzip** them!

```diff
---!! It is recommenced to use SYNTHMAG templates for target >~3500K, 
---!!                  and use PHOENIX templates with target <~3500K!!
```
## SYNTHMAG
The following templates are all available on the google drive. Template named `syntheticstellar*` are generated from the `SYNTHMAG IDL code` (Kochukhov 2007) with atmospheres profile from the PHOENIX NextGen model (Hauschildt et al. 1999) and the stellar lines from the VALD database (Ryabchikova et al. 2015). 

```none
syntheticstellar_hband_T3000_logg4.0_0.0kG.csv.gz	syntheticstellar_kband_T3800_logg3.5_0.5kG.csv.gz
syntheticstellar_hband_T3500_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T3800_logg4.0_0.0kG.csv.gz
syntheticstellar_hband_T4000_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T3800_logg4.0_0.5kG.csv.gz
syntheticstellar_hband_T5000_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T3800_logg4.5_0.0kG.csv.gz
syntheticstellar_hband_T5800_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T3900_logg3.5_0.0kG.csv.gz
syntheticstellar_hband_T6200_logg3.5_0.0kG.csv.gz	syntheticstellar_kband_T3900_logg4.0_0.0kG.csv.gz
syntheticstellar_hband_T6200_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T3900_logg4.5_0.0kG.csv.gz
syntheticstellar_hband_T6400_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4000_logg3.5_0.0kG.csv.gz
syntheticstellar_hband_T6600_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4000_logg3.5_2.2kG.csv.gz
syntheticstellar_kband_T3000_logg4.0_0.0kG.csv.gz	syntheticstellar_kband_T4000_logg4.0_0.0kG.csv.gz
syntheticstellar_kband_T3200_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4000_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3300_logg5.0_0.0kG.csv.gz	syntheticstellar_kband_T4000_logg4.5_2.5kG.csv.gz
syntheticstellar_kband_T3400_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4100_logg3.5_0.0kG.csv.gz
syntheticstellar_kband_T3500_logg3.5_0.0kG.csv.gz	syntheticstellar_kband_T4100_logg4.0_0.0kG.csv.gz
syntheticstellar_kband_T3500_logg4.0_0.0kG.csv.gz	syntheticstellar_kband_T4100_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3500_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4200_logg4.0_0.0kG.csv.gz
syntheticstellar_kband_T3600_logg3.5_0.0kG.csv.gz	syntheticstellar_kband_T4200_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3600_logg3.5_0.5kG.csv.gz	syntheticstellar_kband_T4200_logg5.0_0.0kG.csv.gz
syntheticstellar_kband_T3600_logg4.0_0.0kG.csv.gz	syntheticstellar_kband_T4300_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3600_logg4.0_0.5kG.csv.gz	syntheticstellar_kband_T4400_logg4.0_0.0kG.csv.gz
syntheticstellar_kband_T3600_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T4400_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3700_logg3.5_0.0kG.csv.gz	syntheticstellar_kband_T4600_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3700_logg3.5_0.5kG.csv.gz	syntheticstellar_kband_T4800_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3700_logg4.0_0.0kG.csv.gz	syntheticstellar_kband_T4800_logg5.0_0.0kG.csv.gz
syntheticstellar_kband_T3700_logg4.0_0.5kG.csv.gz	syntheticstellar_kband_T5000_logg4.5_0.0kG.csv.gz
syntheticstellar_kband_T3700_logg4.5_0.0kG.csv.gz	syntheticstellar_kband_T5000_logg5.0_0.0kG.csv.gz
syntheticstellar_kband_T3800_logg3.5_0.0kG.csv.gz	syntheticstellar_kband_T6400_logg4.5_0.0kG.csv.gz
```

> These can generally be applied to stars of the same log(g) and a Teff within a few hundred K of the template. 

Call these templates via for example:
```shell 
(igrins_rv) ~$ python main_stepX.py [target name] -HorK H -t synthetic -temp 5000 -logg 4.5
```

If you used any of the above templates, please cite these papers:
* [Kochukhov 2007 (Spectrum synthesis for magnetic, chemically stratified stellar atmospheres)](https://ui.adsabs.harvard.edu/abs/2007pms..conf..109K/abstract)
* [Hauschildt et al. 1999 (The NEXTGEN Model Atmosphere Grid. II. Spherically Symmetric Model Atmospheres for Giant Stars with Effective Temperatures between 3000 and 6800 K)](https://ui.adsabs.harvard.edu/abs/1999ApJ...525..871H/abstract)
* [Ryabchikova et al. 2015 (A major upgrade of the VALD database)](https://iopscience.iop.org/article/10.1088/0031-8949/90/5/054005)

## PHOENIX template from the Gottingen spectral library
Several sets of [Gottingen spectral templates](http://phoenix.astro.physik.uni-goettingen.de/?page_id=15) are available.

```none
PHOENIX-lte02600-4.00-0.0_contadj.csv.gz	PHOENIX-lte03800-3.50-0.0_contadj.csv.gz
PHOENIX-lte02700-4.00-0.0_contadj.csv.gz	PHOENIX-lte03800-4.50-0.0_contadj.csv.gz
PHOENIX-lte02800-4.00-0.0_contadj.csv.gz	PHOENIX-lte03900-3.50-0.0_contadj.csv.gz
PHOENIX-lte02900-4.00-0.0_contadj.csv.gz	PHOENIX-lte03900-4.50-0.0_contadj.csv.gz
PHOENIX-lte03000-4.00-0.0_contadj.csv.gz	PHOENIX-lte04000-3.50-0.0_contadj.csv.gz
PHOENIX-lte03000-4.50-0.0_contadj.csv.gz	PHOENIX-lte04000-4.50-0.0_contadj.csv.gz
PHOENIX-lte03100-4.50-0.0_contadj.csv.gz	PHOENIX-lte04100-3.50-0.0_contadj.csv.gz
PHOENIX-lte03200-3.50-0.0_contadj.csv.gz	PHOENIX-lte04100-4.50-0.0_contadj.csv.gz
PHOENIX-lte03200-4.00-0.0_contadj.csv.gz	PHOENIX-lte04200-3.50-0.0_contadj.csv.gz
PHOENIX-lte03200-4.50-0.0_contadj.csv.gz	PHOENIX-lte04200-4.50-0.0_contadj.csv.gz
PHOENIX-lte03300-3.50-0.0_contadj.csv.gz	PHOENIX-lte04300-3.50-0.0_contadj.csv.gz
PHOENIX-lte03300-4.00-0.0_contadj.csv.gz	PHOENIX-lte04300-4.50-0.0_contadj.csv.gz
PHOENIX-lte03300-4.50-0.0_contadj.csv.gz	PHOENIX-lte04400-3.50-0.0_contadj.csv.gz
PHOENIX-lte03400-3.50-0.0_contadj.csv.gz	PHOENIX-lte04400-4.50-0.0_contadj.csv.gz
PHOENIX-lte03400-4.00-0.0_contadj.csv.gz	PHOENIX-lte04500-3.50-0.0_contadj.csv.gz
PHOENIX-lte03400-4.50-0.0_contadj.csv.gz	PHOENIX-lte04500-4.50-0.0_contadj.csv.gz
PHOENIX-lte03500-3.50-0.0_contadj.csv.gz	PHOENIX-lte04600-3.50-0.0_contadj.csv.gz
PHOENIX-lte03500-4.50-0.0_contadj.csv.gz	PHOENIX-lte04600-4.50-0.0_contadj.csv.gz
PHOENIX-lte03600-3.50-0.0_contadj.csv.gz	PHOENIX-lte04700-3.50-0.0_contadj.csv.gz
PHOENIX-lte03600-4.50-0.0_contadj.csv.gz	PHOENIX-lte04700-4.50-0.0_contadj.csv.gz
PHOENIX-lte03700-3.50-0.0_contadj.csv.gz	PHOENIX-lte04800-4.50-0.0_contadj.csv.gz
PHOENIX-lte03700-4.50-0.0_contadj.csv.gz	PHOENIX-lte05000-4.50-0.0_contadj.csv.gz
```


Call these templates with (for example):
```shell 
(igrins_rv) ~$ python main_stepX.py [target name] -HorK K -t PHOENIX -temp 3700 -logg 4.5
```
## User supply templates
If you would like to supply your own template, the filename must follow the convention 
`user_T[temperature]_logg[logg]_[band]band.txt`, the file must be placed in ``./Engine/user_templates/`` , and the file itself must follow the format:

```python
 #      wave           flux
       19380.000      0.99825212
       19380.036      0.99824445
       19380.073      0.99823262
```

with the same first line listed here, followed by the **(vacuum) wavelength** in **angstroms** and the normalized flux. There should be no trailing blank lines at the end of the file.

```diff
---All stellar templates need to be **flattened** before use!!---
```
