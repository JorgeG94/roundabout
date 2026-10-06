.. _references:

-----------------------------------
Acknowledgements and references
-----------------------------------

.. contents::
   :local:


NOAA-GFDL, MOM6 and SIS2
========================

The numerical methods in Roundabout were developed at NOAA's Geophysical
Fluid Dynamics Laboratory (NOAA-GFDL) by the developers of MOM6 and SIS2.
Roundabout adapts them to GPUs. Its own work is the port: rewriting those
methods around Fortran ``do concurrent``, with OpenACC for data movement and
reductions, so the same source runs on NVIDIA GPUs and on CPUs. Full credit
for the numerics goes to the MOM6 and SIS2 developers.

* **MOM6**, the Modular Ocean Model version 6,
  https://github.com/NOAA-GFDL/MOM6. The ocean dynamical core follows it:
  the split-explicit barotropic solver, continuity-PPM, the finite-volume
  pressure gradient, ALE remapping, the vertical and lateral closures, and
  the vertical-friction and bottom-boundary-layer treatment.
  ``tools/om_topo_limit.f90`` is derived from MOM6 code.
* **SIS2**, the Sea Ice Simulator version 2,
  https://github.com/NOAA-GFDL/SIS2. The sea-ice component follows it:
  thermodynamics, the ice-thickness distribution, category transport, the
  C-grid elastic-viscous-plastic rheology and the ice-ocean coupling.


The physics oracle
==================

MOM6 is Roundabout's physics oracle. Wherever Roundabout implements an
algorithm MOM6 also has, MOM6's behaviour is the reference it is checked
against. Divergences from it are deliberate and documented in
``docs/CLOSURE_MATRIX.md`` and ``docs/CAPABILITIES_AND_LIMITATIONS.md``. MOM6
and SIS2 runtime parameters (for example ``BT_STRONG_DRAG`` or
``HVEL_SCHEME``) are named where a Roundabout option corresponds to one. The
code itself cites the physics and the paper behind each algorithm.


Licensing
=========

MOM6 and SIS2 are licensed under the Apache License, Version 2.0. The code
ported from them is used and modified under that licence. The repository
carries a copy of it in ``LICENSES/Apache-2.0.txt``, and a ``NOTICE`` file
that states what was ported and that it was modified. Roundabout's own code
is released under the MIT License (``LICENSE``). Neither NOAA-GFDL nor the
MOM6 or SIS2 developers endorse Roundabout, and they are not responsible for
it.


References
==========

These entries carry the bibliographic details given in the module docstrings
that use them.

* Adcroft, A. (2013): "Representation of topography by porous barriers and
  objective interpolation of topographic data." *Ocean Modelling* 67, 13-27.
  (Porous barriers.)
* Adcroft, A., Hallberg, R. and Harrison, M. (2008): *Ocean Modelling* 24,
  1-2. Analytic finite-volume pressure gradient.
* Asay-Davis, X. S. et al. (2016): "Experimental design for three interrelated
  marine ice sheet and ocean model intercomparison projects: MISMIP v. 3
  (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1)." *Geosci.
  Model Dev.* 9, 2471-2497, doi:10.5194/gmd-9-2471-2016. (ISOMIP+.)
* Beckmann, A. and Haidvogel, D. B. (1993): *J. Phys. Oceanogr.* 23,
  1736-1753. (Sigma-coordinate slope criterion.)
* Burchard, H. et al. (2022): *Ocean Modelling* 179, 102119. (Ice-shelf
  cavity fluxes.)
* Chelton, D. B., deSzoeke, R. A., Schlax, M. G., El Naggar, K. and
  Siwertz, N. (1998): *J. Phys. Oceanogr.* 28, 433-460. (Deformation
  radius.)
* Eady, E. T. (1949): *Tellus* 1, 33-52. (Baroclinic growth rate.)
* Gill, A. E. (1982): *Atmosphere-Ocean Dynamics*.
* Hallberg, R. (2013): *Ocean Modelling* 72, 92-103. (Resolution-aware
  eddy parameterisation.)
* Hallberg, R. and Adcroft, A. (2014): "An order-invariant real-to-integer
  conversion sum." *Parallel Computing* 40(5-6),
  doi:10.1016/j.parco.2014.04.007. (Reproducing sums.)
* Haney, R. L. (1991): *J. Phys. Oceanogr.* 21, 610-619. (Hydrostatic
  consistency.)
* Holland, D. M. and Jenkins, A. (1999): *J. Phys. Oceanogr.* 29,
  1787-1800. (Ice-shelf melt thermodynamics.)
* Jenkins, A., Nicholls, K. W. and Corr, H. F. J. (2010): "Observation and
  parameterization of ablation at the base of Ronne Ice Shelf, Antarctica."
  *J. Phys. Oceanogr.* 40, 2298-2312.
* Killworth, P. D. and Edwards, N. R. (1999): *J. Phys. Oceanogr.* 29,
  1221-1238. (Bottom boundary layer of prescribed thickness.)
* McPhee, M. G., Maykut, G. A. and Morison, J. H. (1987): "Dynamics and
  thermodynamics of the ice/upper ocean system in the marginal ice zone of the
  Greenland Sea." *J. Geophys. Res.* 92(C7), 7017-7031.
* Murray, R. J. (1996): *J. Comput. Phys.* 126, 251-273. (Tripolar grid.)
* Reichl, B. G. and Hallberg, R. (2018): *Ocean Modelling* 132. (ePBL.)
* Roquet, F., Madec, G., McDougall, T. J. and Barker, P. M. (2015):
  "Accurate polynomial expressions for the density and specific volume of
  seawater using the TEOS-10 standard." *Ocean Modelling* 90, 29-43.
* Simmons, H. L., Jayne, S. R., St Laurent, L. C. and Weaver, A. J. (2004):
  *Ocean Modelling* 6. (Tidal mixing.)
* Smith, R. D. and McWilliams, J. C. (2003): "Anisotropic horizontal
  viscosity for ocean models." *Ocean Modelling* 5(2).
* Visbeck, M., Marshall, J., Haine, T. and Spall, M. (1997): *J. Phys.
  Oceanogr.* 27, 381-402. (Eddy diffusivity.)
* White, L., Adcroft, A. and Hallberg, R. (2009): *J. Comput. Phys.* 228.
  (High-order in-layer reconstruction for the pressure gradient.)
* Yung, C. K. et al. (2025): "Sensitivity of Antarctic ice shelf melt to the
  ice-ocean boundary layer parameterisation." *The Cryosphere* 19, 5827-5861.
* Yung, C. K. et al. (2026): *The Cryosphere* 20, 2053-2088.


Also cited by author and year
=============================

The following works are cited by author and year in the module docstrings,
next to the algorithm each one describes. Their full bibliographic entries
have not yet been added to this page.

Arakawa & Hsu (1990); Balsara & Shu (2000); Bleck (2002); Bodner et al.
(2023); Borges et al. (2008); Bryan & Lewis (1979); Colella & Woodward
(1984); Egbert & Ray (2001); Ferreira, Marshall & Campin (2010); Flather
(1976); Fox-Kemper, Ferrari & Hallberg (2008); Gent & McWilliams (1990);
Gent et al. (1995); Griffies (1998); Griffies & Hallberg (2000); Hallberg
(1997); Hallberg & Adcroft (2009); Harrison & Hallberg (2008); Henyey, Wright
& Flatté (1986); Hibler (1979); Jackson, Hallberg & Legg (2008); Jayne & St
Laurent (2001); Large et al. (1994); Leith (1968); Lipscomb (2001); Losch
(2008); Millero (1978); Orlanski (1976); Pacanowski & Philander (1981); Redi
(1982); Reichl & Li (2019); Sadourny (1975); Shao (2016); Smagorinsky
(1963); St Laurent et al. (2002); Vreugdenhil & Taylor (2019); White &
Adcroft (2008); Winton (2000); Wright (1997).
