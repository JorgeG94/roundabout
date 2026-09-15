.. _running_simulations:

-------------------
Running simulations
-------------------

.. contents::
   :local:

Everything about a Roundabout run is in its namelist. This page covers how
that file is organised, then walks through three real cases.


Anatomy of a namelist
=====================

A Roundabout namelist is a sequence of Fortran namelist groups. Order does
not matter, and any group you leave out takes its defaults:

.. code-block:: fortran

   &sim_nml
      sim_type = "ocean"
   /

   &grid_nml
      nx     = 50
      ny     = 50
      dx     = 2000.0
      dy     = 2000.0
      nghost = 2
   /

Groups fall into three families.

**Core groups** describe the run itself: ``&sim_nml`` (``sim_type``, whose
only accepted value is ``"ocean"``), ``&grid_nml`` (cell counts, spacing,
halo width), ``&time_nml``, ``&physics_nml``, ``&nonhydrostatic_nml``
(which, despite the name, is where ``nz_layers`` lives), ``&vcoord_nml``,
``&tracer_nml``, ``&output_nml`` and ``&logging_nml``.

**Ocean groups** are named ``&ocean_<concern>_nml`` and hold the physics:
``&ocean_grid_nml``, ``&ocean_topo_nml``, ``&ocean_ic_nml``,
``&ocean_bc_nml``, ``&ocean_coriolis_nml``, ``&ocean_pgf_nml``,
``&ocean_eos_nml``, ``&ocean_hvisc_nml``, ``&ocean_vmix_nml``,
``&ocean_vdiff_nml``, ``&ocean_bdrag_nml``, ``&ocean_bt_nml``,
``&ocean_thermo_nml``, ``&ocean_continuity_nml``, ``&ocean_diag_nml``, and
one per optional closure — ``&ocean_epbl_nml``, ``&ocean_gm_nml``,
``&ocean_meke_nml``, ``&ocean_foxkemper_nml``, and so on.

The naming rule is that the ``ocean_`` prefix is dropped from each key
inside its group: the viscosity knob is ``&ocean_hvisc_nml nu_h``, not
``ocean_nu_h``.

**Groups you will rarely touch** exist for machinery — the schema carries
57 groups in total.

.. note::

   The complete knob-to-default map is ``docs/generated_nml_knobs.md`` in
   the repository, produced by the ``rdb_nml_doc`` binary directly from the
   configuration schema. It is generated, so it cannot drift. Consult it
   rather than any hand-written list, this page included.

Units and time
--------------

``&time_nml time_unit`` — accepting ``"s"`` (the default), ``"min"``,
``"hr"``, ``"day"`` and ``"year"`` — scales several time-valued knobs:
``t_end``, ``&logging_nml status_interval`` and ``&ocean_diag_nml dt_out``.

.. warning::

   ``dt_fixed`` is **always seconds**, regardless of ``time_unit``. A
   namelist with ``time_unit = "day"`` and ``t_end = 10.0`` runs for ten
   days with a ``dt_fixed = 1200.0`` step of twenty minutes. Mixing the two
   conventions up is the most common namelist error.

The driver requires ``dt_fixed > 0``; there is no adaptive-CFL helper.

Validation
----------

The configuration is checked before anything is allocated, and the checks
are **fail-loud by design**. An unknown key, a value outside an enum, or a
combination the code cannot honour stops the run with a message naming the
problem and, usually, the fix. Nothing is silently ignored and nothing
silently falls back to a default.

This is worth leaning on. If a configuration runs, the combination of
options it names is one the code actually supports.


Case 1: an at-rest test
=======================

``validation_examples/ocean/seamount/seamount.nml``

The most informative kind of test: one where the right answer is *nothing*.

.. code-block:: fortran

   &grid_nml
      nx     = 50
      ny     = 50
      dx     = 2000.0
      dy     = 2000.0
      nghost = 2
   /

   &time_nml
      t_end        = 2592000.0   ! 30 days
      dt_fixed     = 300.0       ! 5 min outer
      cfl_interval = 1
   /

   &physics_nml
      coriolis_f    = 0.0
      wind_stress_x = 0.0
      wind_stress_y = 0.0
   /

   &nonhydrostatic_nml
      nz_layers = 15
   /

   &vcoord_nml
      vcoord_type = "zstar_sigma"
   /

   &tracer_nml
      initial_temperature = 15.0
      initial_salinity    = 35.0
      T_init_surface      = 0.0
      T_init_bottom       = 0.0
   /

A 100 km closed basin with a single Gaussian seamount rising from 4000 m to
200 m over a 25 km half-width — a maximum slope of about 150 m/km,
comparable to a real shelf break. Fifteen layers.

Every source of motion has been removed deliberately. Temperature and
salinity are uniform, so density is uniform and there is no available
potential energy. :math:`f = 0`, so nothing can resonate. No wind, no
surface flux, no sponges, closed walls.

Therefore **any velocity the model produces is a numerical artefact**,
almost certainly in the pressure gradient over the sloping coordinate
surfaces (:ref:`pressure_gradient`). There is nothing to tune and no way
for a wrong answer to look plausible. The suite gates this case at a mean
kinetic energy corresponding to about 1 cm/s of spurious current — the
Beckmann & Haidvogel (1993) bar for sigma-coordinate pressure-gradient
error — and the current code sits about an order of magnitude under it.

The two settings that make it a *dynamical-core* test rather than a physics
test are worth noting:

.. code-block:: fortran

   &ocean_hvisc_nml
      nu_h = 0.0
      nu_4 = 0.0
   /

   &ocean_vmix_nml
      use_closure = .false.
      use_kpp     = .false.
   /

With viscosity and mixing off, nothing can damp an artefact into
invisibility.


Case 2: an analytical benchmark
===============================

``validation_examples/ocean/eady/eady.nml``

The classical Eady (1949) front: a 250 km × 250 km reentrant channel,
50 × 50 cells at 5 km, ten layers over a 1000 m flat bottom, on an
:math:`f`-plane at :math:`f = 10^{-4}\ \mathrm{s^{-1}}`.

.. code-block:: fortran

   &ocean_bc_nml
      west = "periodic"
      east = "periodic"
   /

   &ocean_ic_nml
      ic_config      = "eady"
      alpha_T        = 0.17
      rho_0          = 1035.0
      eady_dT_dz     = 0.01
      eady_dT_dy     = -2.0e-5
      eady_T_ref     = 10.0
      eady_pert_amp  = 1.0e-3
      eady_pert_seed = 12345
   /

The initial condition is a meridional temperature gradient in thermal-wind
balance with a vertically sheared zonal jet, seeded with
:math:`\pm 0.5\ \mathrm{mK}` of white noise to break the symmetry.

What makes this case valuable is that **linear theory predicts the answer**.
The stratification and shear give a deformation radius
:math:`R_d = NH/f = 40.1\ \mathrm{km}` and a Richardson number of 155. In a
channel the zonal wavenumber is quantised, and the mode that grows is the
125 km one, at a predicted rate of :math:`1.98\times10^{-6}\
\mathrm{s^{-1}}`. The model produces :math:`1.93\times10^{-6}`, about 3 %
low.

Three points of practice come out of this case.

**Measure the right thing.** The cross-channel velocity :math:`v` is zero
in the basic state, so :math:`\max|v|` is a *pure* perturbation measure.
``En`` is not — it is dominated by the :math:`\pm 0.16\ \mathrm{m\,s^{-1}}`
background jet, and it actually drifts *down* by about 6 % over the first
45 days as the jet drains against the side walls, turning up only once the
instability is large. Reading ``En`` as the growth signal here gives
exactly the wrong conclusion.

**Viscosity competes with the mode you are measuring.** This file runs
``nu_h = 20``. At ``nu_h = 100`` the growth rate is only slightly lower,
but the jet drains faster and the mode needs about ten extra days to clear
the noise floor. A closure strong enough to stabilise a run is also strong
enough to suppress the physics you are trying to validate.

**Check that growth is wavenumber-selective.** Genuine baroclinic
instability picks a mode: here :math:`k_x = 2` grows while
:math:`k_x \ge 5` decay. Growth that is the *same* at every wavenumber is
not baroclinic instability, whatever its rate.

The run stops at 60 days, deliberately: the linear phase is clean to about
day 62, and the front's nonlinear collapse past that point is outside the
validated envelope at this resolution and time step. This is a
**linear-growth benchmark by construction**.


Case 3: a reference reproduction
================================

``validation_examples/ocean/double_gyre/double_gyre_mom6.nml``

A wind-driven double gyre mirroring the MOM6 ``ocean_only/double_gyre``
example: 44 × 40 cells over 22° × 20° at mid-latitude, two layers,
2000 m deep, spoon bathymetry, a two-gyre wind stress of 0.1 Pa,
:math:`\Delta t = 1200\ \mathrm{s}`, ten days.

.. code-block:: fortran

   &ocean_grid_nml
      axis_units = "degrees"
      len_lon    = 22.0
      len_lat    = 20.0
      rad_earth  = 6.378e6
   /

   &time_nml
      t_end        = 10.0
      time_unit    = "day"
      dt_fixed     = 1200.0
      cfl_interval = 1
   /

   &nonhydrostatic_nml
      nz_layers = 2
   /

   &ocean_pgf_nml
      form        = "gprime"
      gprime_gfs  = 0.98
      gprime_gint = 0.0098
      maxvel      = 6.0
   /

Note the grid. ``dx`` and ``dy`` are **derived** from the degree extent, not
stated — pinning them by hand once produced a 16 700 km basin that ran
perfectly happily. Note also that ``nz_layers = 2`` is not a choice: the
``gprime`` pressure-gradient form is a two-layer kernel
(:ref:`pressure_gradient`), and the configuration is refused with any other
layer count.

Reproducing another model exactly is rarely possible, and the header is
candid about where this one diverges:

* MOM6 runs this case adiabatically with thermodynamics disabled.
  Roundabout's multilayer kernel always advects :math:`T` and :math:`S`, so
  the equivalent is reached by making the linear EOS insensitive to them —
  the density then does not depend on the tracers, and any tracer drift
  cannot reach the dynamics.
* MOM6's linear drag is distributed over the bottom 10 m. The effective
  Rayleigh rate is matched,
  :math:`r = C_d \cdot U_{bg}/H_{bbl} = 2.5\times10^{-3} \times 0.1/10 =
  2.5\times10^{-5}\ \mathrm{s^{-1}}`, but the bottom-boundary-layer
  footprint differs, so bottom-current structure will not match.
* The density contrast of :math:`1.0\ \mathrm{kg\,m^{-3}}` gives
  :math:`g' \approx 0.00948` against MOM6's 0.0098 — a 3 % difference.

The general lesson is that a reference reproduction is only as good as its
list of known divergences. A case that claims to match another model and
does not enumerate where it cannot is not a validation.


Output and diagnostics
======================

``&ocean_diag_nml`` controls what is written.

.. code-block:: fortran

   &ocean_diag_nml
      enabled  = .true.
      filename = "double_gyre"
      dt_out   = 1.0        ! time_unit = "day" above
   /

   &output_nml
      output_dir = "./out_double_gyre_mom6"
   /

Files are written **per rank** as
``<output_dir>/<filename>_rank_NNNNNN.nc``, with no gather;
``tools/merge_output.py`` stitches them offline.

Selecting fields
----------------

``&ocean_diag_nml diags`` takes a token list naming what to output. Each
token may carry suffixes selecting a time operation and an output vertical
grid, so a single knob covers instantaneous values, time means, extrema and
integrals, on layers or remapped to :math:`z`, sigma, z\* or density bins.

``vgrid`` sets the default output vertical grid for the whole stream —
``layer`` (the default, native model layers), ``z_fixed``, ``sigma``,
``zstar`` or ``density``. Remapping happens at write time; the model always
computes on its own coordinate.

Cadence and cost
----------------

Diagnostics are filled **on the device** and pulled back to the host only
when the cadence fires. Pulling every step instead costs about three times
as much.

.. warning::

   ``dt_out`` is easy to underestimate. A frame of a
   :math:`600 \times 600 \times 50` case is roughly 0.9 GB; eighty frames
   of that fill a disk. Set ``dt_out`` from how much output you can store,
   not from how much you would like.

``&ocean_diag_nml output_precision = "single"`` halves the bytes per frame
on the diagnostic stream, which is write-bandwidth-bound. It is deliberately
scoped to that one stream: restarts and the console conservation totals stay
at working precision, so a lossy restart is not requestable. Opting in means
regenerating any baseline that byte-compares diagnostic NetCDF.

Compression is available as ``&output_nml compress_output`` with
``compress_level`` from 1 to 9, typically giving three- to five-fold
reduction.
