---
title: "IGRINS RV: A Python Package for Precision Radial Velocities with Near-Infrared Spectra"
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

The relative radial velocity of a star with respect to the Sun can be calculated from its electromagnetic spectrum using the Doppler Effect. This line-of-sight motion, called the Radial Velocity (RV), is an essential tool for astrophysicists. RVs are not only used to detect and characterize exoplanets, but also play a key role in studies of binary stars, star clusters, and moving group member identification.

In the past decade, RVs have primarily been measured from spectra in the optical wavelength regime. This is partly because of advancements in detector technology, but also because of the paucity of Earth's atmospheric absorption features (telluric lines) in the optical. Yet for a fainter, cooler, smaller stellar object like an M type star (stars with mass less than half of the Sun), which emits more energy in the Near-Infrared (NIR), observations in the NIR can save a considerable amount of exposure time. Also, the M type star is the most common type of star. This along with its size increases the detectability of Earth-like planets around them. Moreover, the stellar activity that can drive false positive exoplanet detections, e.g., star spots carried into view by stellar rotation or gas accretion from the circumstellar disk in young star systems, is shown to be less severe in the NIR compared to optical.

``IGRINS RV`` is a pipeline built for extracting precision RVs from spectra taken with the Immersion GRating INfrared Spectrometer (IGRINS) spectrograph [@yuk10; @park14; @mace16; @mace18]. This pipeline is built on the forward-modeling methodology that was successfully applied to CSHELL and PHOENIX spectra in the past [@croc12]. Comparing to obtain RVs via cross-correlation with stellar templates used by past studies, ``IGRINS RV`` gives three times better RV precision to about 25--30 m/s around narrowline stars in both H and K bands shown by years of monitoring on two RV standard stars, GJ\ 281 and HD\ 26257. This improvement is because of the use of a more robust approach to wavelength calibration and better telluric modeling through the generation of high resolution synthetic templates from fits to A0 standard observations. ``IGRINS RV`` has also demonstrated its effectiveness in identifying orbiting companions by successfully recovering the planet-induced RV signals of HD\ 189733 and Tau\ Boo\ A. ``IGRINS RV`` lets users choose to obtain absolute RVs or relative RVs, depending on whether their priority is precise RV monitoring or more coarse RV characterization. ``IGRINS RV`` requires that the ``igrins plp v2.2.0`` [@leegul17] and ``Telfit`` [@gull14] packages be pre-installed. Detailed documentation and tutorials can be found on the GitHub wiki page.


# Acknowledgements

Partial support for this work was provided by NASA Exoplanet Research Program grant 80-NSSC19K-0289 to L. Prato.  We are grateful for the generous donations of John and Ginger Giovale, the BF Foundation, and others which made the IGRINS-DCT program possible. Additional funding was provided by the Mt. Cuba Astronomical Foundation and the Orr Family Foundation. IGRINS was developed under a collaboration between the University of Texas at Austin and the Korea Astronomy and Space Science Institute (KASI) with the financial support of the US National Science Foundation under grant AST-1229522 and AST-1702267, of the University of Texas at Austin, and of the Korean GMT Project of KASI.

# References
