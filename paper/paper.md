---
title: "IGRINS RV: A Python package for precision radial velocities with Near-Infrared Spectra"
tags:
  - Python
  - astronomy
  - radial velocity
  - IGRINS
authors:
 - name: Shih-Yun Tang
   orcid: 0000-0003-4247-1401
   affiliation: "1, 2"
 - name: Asa G. Stahl
   orcid: 0000-0002-0848-6960
   affiliation: "3"
 - name: Christopher M. Johns-Krull
   affiliation: "3"
 - name: L. Prato
   orcid: 0000-0001-7998-226X
   affiliation: "1, 2"
 - name: Joe Llama
   orcid: 0000-0003-4450-0368
   affiliation: "1"
affiliations:
 - name: Lowell Observatory, 1400 W. Mars Hill Road, Flagstaff, AZ 86001, USA
   index: 1
 - name: Department of Astronomy and Planetary Sciences, Northern Arizona University, Flagstaff, AZ 86011, USA
   index: 2
 - name: Department of Physics and Astronomy, Rice University, 6100 Main Street, Houston, TX 77005, USA
   index: 3
aas-doi:
aas-journal: Astronomical Journal
date:
bibliography: master.bib
---

# Summary

The relative velocity of a star compare to the Sun can be calculated from the spectrum using the Doppler Effect. This line-of-site relative velocity is called the Radial Velocity (RV). RV has been an essential parameter for astrophysicists not only to detect exoplanets, but also to study the orbit of binaries, and plays a key role in star clusters and moving groups member identification.

In the past decay, RV was mostly measured from spectra in the optical wavelength regime. This is partly because of the more advanced development in the detector, but mostly because of the Earth's atmosphere absorption features (telluric lines) are far less in the optical wavelength. Yet, when astronomer gains their interest in the much fainter objects, stars with mass less than half of the Sun (M type star), RV measurement in a longer wavelength, the Near-Infrared (NIR) regime, is needed. Moreover, a fainter, cooler stellar object emits more energy in the NIR. This means we can have a less exposure time when observing an M type star in the NIR than that in the optical. Moreover, the fake planet like RV signal, e.g., dusk accretion on to the stellar surface or sun spots rotation period, is shown to be less server in the NIR compared to optical.

<!--However, the required precision in RVs is quite different from studies to studies. Only a few km/s precision is needed to determine the orbit of a binary system and for the identification for cluster members. Nevertheless, a 10 m/s precision is required for finding Earth-size exoplanets.-->

``IGRINS RV`` is a pipeline build for extracting precision RVs for spectra of the Immersion GRating INfrared Spectrometer (IGRINS) spectrograph [@yuk10, @park14, @mace16, @mace18]. This pipeline is built on the forward-modeling methodology that performed well with the CSHELL and PHOENIX spectra [@croc12]. However, the RV precision had improved nearly 10 times to about 25--50 m/s with a more robust approach on wavelength calibration. Users can choose to obtain an absolute RV solution or a relative RV solution. The later comes with higher precision by subtracting the mean of individual RV from each echelle order before combined. ``IGRINS RV`` requires `` igrins plp v2.2`` [@leegul17] and ``Telfit`` [@gull14] package pre-installed. Detailed documentation and tutorials can be found on the GitHub wiki page.


# Acknowledgements

If any.

# References
