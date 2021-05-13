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
aas-doi: 10.3847/1538-3881/abf5e7
aas-journal: Astronomical Journal
date:
bibliography: master.bib
---

# Summary

The relative radial velocity of a star with respect to the Sun can be calculated from its electromagnetic spectrum using the Doppler Effect. This line-of-sight motion, called the Radial Velocity (RV), is an essential tool for astrophysicists. RVs are not only used to detect and characterize exoplanets, but also play a key role in studies of binary stars, star clusters, and moving group member identification.

In the past decade, RVs have primarily been measured from spectra in the optical wavelength regime. This is partly because of advancements in detector technology, but also because of the paucity of Earth's atmospheric absorption features (telluric lines) in the optical. Yet for a fainter, cooler, smaller stellar object like an M type star (stars with mass less than half of the Sun), which emits more energy in the Near-Infrared (NIR), observations in the NIR can save a considerable amount of exposure time. Also, the M type star is the most common type of star. This along with its size increases the detectability of Earth-like planets around them. Moreover, the stellar activity that can drive false positive exoplanet detections, e.g., star spots carried into view by stellar rotation is shown to be less severe in the NIR compared to optical.


# Statement of need

Current RV pipelines and techniques that can deliver RV precision in tens of m/s (or beyond) in the NIR, e.g., ``PySHELL`` [@cale19], ``wobble`` [@bede19], ``SERVAL`` [@zech18] or the PCA-based cross-correlating method used for the SPIRou spectrograph [@mout20], all require instruments that are highly stabilized and with well-characterized wavelength solutions. For example, the iSHELL spectrograph can be equipped with the methane isotopologue gas cell, and the SPIRou and the CARMENES (NIR channel) spectrograph come with uranium-neon hollow-cathode lamps and stabilized Fabry-Perot etalons. The Immersion GRating INfrared Spectrometer (IGRINS) spectrograph [@yuk10; @park14; @mace16; @mace18], on the other hand, was not designed to be RV-stable and comes with no means of wavelength calibration accurate enough to achieve RVs precise to tens of m/s using existing techniques.  A new approach to extract precision RVs is needed for an echelle spectrograph like IGRINS, which still offers fertile ground for RV science with its high resolution (R ~ 45,000) and broad spectral grasp (the full H and K bands).

``IGRINS RV`` is a pipeline tailored for extracting precision RVs from spectra taken with IGRINS on different facilities. This pipeline is built on the forward-modeling methodology that was successfully applied to CSHELL and PHOENIX spectra [@croc12] that utilized telluric lines as a common-path wavelength calibrator. Compared to RVs obtained by cross-correlation with stellar templates adopted by past studies, ``IGRINS RV`` gives three times higher RV precision, about 25--30 m/s, around narrow-line stars in both H and K bands, shown by years of monitoring on two RV standard stars, GJ\ 281 and HD\ 26257. ``IGRINS RV`` also pushes this technique, using telluric lines as wavelength calibrator for RV calculation, to its limits as studies found the stability of the telluric lines is about 10--20 m/s. Moreover, ``IGRINS RV`` is also tailored to take into account specific aspects of the IGRINS instrument, like the variations in spectral resolution across the detector and the year-long K band detector defocus.

``IGRINS RV`` has demonstrated its effectiveness in identifying orbiting companions by successfully recovering the planet-induced RV signals of HD\ 189733 and Tau\ Boo\ A. ``IGRINS RV`` lets users choose to obtain absolute RVs or relative RVs, depending on whether their priority is coarse RV characterization or more precise RV monitoring. The code extends the science capabilities of an already powerful spectrograph, which lacked a publicly available RV pipeline until now. It facilitates the detection and/or characterization of exoplanets, binary stars, star clusters, and moving group members, and it enables such studies to be done in a more precise and uniform way.

``IGRINS RV`` makes use of the ``astropy`` [@astropy2013; @astropy2018] on handing sky coordinates and barycentric velocity correction, ``scipy`` [@2020SciPy] and ``numpy`` [@2020NumPy] on mathmatical calculation, ``nlopt`` [@john08; @box1965new] on the optimization process, ``pandas`` [@reback2020pandas; @mckinney2010] on data management, and ``matplotlib`` [@Hunter2007] on plotting. We also used a part of code from `BMC` [@marcos2021] for peak detection. ``IGRINS RV`` requires that the ``igrins plp v2.2.0`` [@leegul17] and ``Telfit`` [@gull14] packages be pre-installed. Detailed documentation and tutorials can be found on the GitHub wiki page.


# Acknowledgements

Partial support for this work was provided by NASA Exoplanet Research Program grant 80-NSSC19K-0289 to L. Prato. CMJ would like to acknowledge partial support for this work provided through grants to Rice University provided by NASA (award 80-NSSC18K-0828) and the NSF (awards AST-2009197 and AST-1461918). We are grateful for the generous donations of John and Ginger Giovale, the BF Foundation, and others which made the IGRINS-DCT program possible. Additional funding was provided by the Mt. Cuba Astronomical Foundation and the Orr Family Foundation. IGRINS was developed under a collaboration between the University of Texas at Austin and the Korea Astronomy and Space Science Institute (KASI) with the financial support of the US National Science Foundation under grant AST-1229522 and AST-1702267, of the University of Texas at Austin, and of the Korean GMT Project of KASI.

# References
