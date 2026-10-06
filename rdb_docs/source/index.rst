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

**The numerical methods in Roundabout were developed at NOAA's Geophysical
Fluid Dynamics Laboratory (NOAA-GFDL). Roundabout adapts them to GPUs.** The
ocean dynamical core follows `MOM6 <https://github.com/NOAA-GFDL/MOM6>`_ and
the sea-ice component follows `SIS2 <https://github.com/NOAA-GFDL/SIS2>`_.
Roundabout's own work is rewriting those methods around ``do concurrent`` and
OpenACC. Full credit for the numerics goes to the MOM6 and SIS2 developers.
MOM6 is also Roundabout's **physics oracle**: where Roundabout implements an
algorithm MOM6 also has, MOM6's behaviour is the reference it is checked
against, and deliberate divergences are documented. MOM6 and SIS2 are
licensed under the Apache License 2.0, and the ported code is used and
modified under it. See :doc:`references` for the full acknowledgement,
licensing and the papers the physics comes from.

Roundabout is a hobby project, written for fun and as the author's
learning project. The author is part of
`MOM6-GPU <https://github.com/MOM6-GPU>`_, the effort to port MOM6 to
GPUs. Roundabout runs in parallel to, and separately from, that work and
the GPU porting efforts at NOAA-GFDL and
`ACCESS-NRI <https://www.access-nri.org.au/>`_, and it is not part of
them.

**Authorship.** Most of Roundabout's code was written by Claude,
Anthropic's AI model, working through Claude Code. The baseline design is the
author's: the architecture, the API, the memory model and the parallelism
model, all designed by Jorge Luis Gálvez Vallejo.

Roundabout makes no claim to be a novel implementation of a general
circulation model. It is an LLM-driven, human-supervised experiment in
taking a general circulation model to GPUs.

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
   references

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide:

   extending
