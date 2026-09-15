.. Roundabout documentation master file.
   You can adapt this file completely to your liking, but it should at
   least contain the root `toctree` directive.

========================
User guide to Roundabout
========================

Roundabout is a GPU-native ocean solver written in modern Fortran. It solves
the hydrostatic, Boussinesq primitive equations in layer form on a
structured Arakawa C-grid, with continuity-PPM layer transport, a
split-explicit time integration that sub-cycles the fast barotropic mode,
and ALE vertical coordinates. It targets regional and global hydrostatic
configurations down to submesoscale-resolving resolution.

Two things distinguish it. The first is that **layer thickness is a
prognostic variable**: continuity is a transport equation rather than a
constraint, so there is no Poisson solve and no elliptic stage anywhere in
the dynamical core. The second is that the parallelism is expressed in
standard Fortran — all data-parallel loops are ``do concurrent``, with
OpenACC used only for data movement and reductions — so one source builds
for NVIDIA GPUs, for CPU multicore, and for plain serial CPU across NVHPC,
gfortran and ifx.

Roundabout is named for Canberra, the one large Australian city with no
ocean, and its roundabouts, which look a lot like eddies.

The API documentation for the code itself — per-module and per-procedure
reference generated from the source — is hosted separately:
https://jorgeg94.github.io/roundabout/

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   getting_started
   theory
   discretisation
   vertical_coordinates
   pressure_gradient
   parameterizations
   grids_and_boundaries
   running_simulations
   python_interface
   validation

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide:

   extending
