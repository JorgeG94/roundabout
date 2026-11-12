.. _sea_ice:

-------
Sea ice
-------

.. contents::
   :local:

Roundabout carries an optional sea-ice component that sits on top of the
ocean surface and exchanges heat, salt and momentum with it. It is a port
of `SIS2 <https://github.com/NOAA-GFDL/SIS2>`_, the Sea Ice Simulator
version 2 from NOAA-GFDL, rewritten around ``do concurrent`` like the rest
of the model. It lives under ``src/core/ice/`` and runs on the ocean's own
C-grid and domain decomposition.

The short version of what it is: **a sea-ice dynamical core and column
model, not yet a complete sea-ice model.** What is there — the column
thermodynamics, the ice-thickness distribution, transport and the
elastic-viscous-plastic rheology — is tested against analytical cases.
What is missing is mostly on the mass budget: there is no ridging, no melt
ponds and no lateral melt (see `Limitations`_).

Everything on this page is **off by default**. A namelist without
``&ocean_ice_nml enable = .true.`` never initialises, maps or steps the
ice, and gives byte-identical results.

.. warning::

   The sea ice is tuned for global, ice-covered configurations and is
   **not** performance-tuned. It runs on every outer step, and the EVP
   dynamics subcycles hundreds of times per step. Expect a large
   throughput hit, and a very slow run on a regional or coastal domain
   where most of the grid is ice-free. The model prints a warning saying
   so whenever the ice is enabled.


What the model carries
======================

Each ocean column holds a sea-ice state split into **thickness
categories**: category 0 is open water and categories ``1..ncat`` hold ice
of increasing thickness. Per category the model stores the area fraction
(summing to one over ``0..ncat``), the ice and snow mass per unit
ice-covered area, the ice and snow **specific enthalpy** (the
thermodynamic prognostic is enthalpy, never temperature) and the ice bulk
salinity. Ice is split into ``nk_ice`` layers (default 2) under one snow
layer. The ice velocity lives on the ocean's C-grid u- and v-faces.

**A single category is a special mode.** With ``ncat = 1`` the model runs
in a legacy *lumped* mode: there is no thickness distribution, masses are
per unit cell area, and the concentration is binary (1 wherever there is
ice, 0 elsewhere). Category transport is refused in this mode. Use
``ncat > 1`` (the default is 5) for anything realistic.


Turning it on
=============

The smallest useful configuration is thermodynamics only. The block below
is the sea-ice part of
``validation_examples/ocean/sea_ice/polar_freezeup_thermo.nml``, a 40 × 30
cell, 10 km, 50 m deep f-plane basin at about 80°N that starts ice-free
just above the freezing point and is cooled at the surface:

.. code-block:: fortran

   &ocean_thermo_nml
      enable_thermodynamics = .true.   ! the ice thermo chain needs this
      q_heat = -60.0                    ! W/m^2, cools the surface to T_f
   /

   &ocean_vmix_nml
      dt_therm_ratio = 1                ! ice thermodynamic cadence
   /

   &ocean_ice_nml
      enable         = .true.
      ncat           = 5
      nk_ice         = 2
      air_temp       = -20.0            ! degC, slab atmosphere
      restore_lambda = 20.0             ! W/m^2/K
      sw_down        = 0.0              ! polar night
      transport      = .false.
      dynamics       = .false.
   /

Ice nucleates as frazil once the surface supercools, and the column model
then grows it over the 90-day run. To add motion, the companion file
``polar_freezeup_dynamics.nml`` applies a steady 0.05 Pa eastward wind
(``&physics_nml wind_stress_x``) and adds, in the same group:

.. code-block:: fortran

      transport      = .true.
      dynamics       = .true.
      evp_sub_steps  = 432
      cfl_trunc      = 0.5     ! SIS2's own default; Roundabout's is 0
      project_ci     = .true.  ! SIS2's own default; Roundabout's is .false.

To start from an existing pack instead of growing one, see
`Initial conditions`_ and ``validation_examples/ocean/sea_ice_pack/``.

.. note::

   Open water never grows ice on its own. The column thermodynamics only
   steps categories that already carry ice, so the **only** way new ice
   appears is the frazil path: the ocean surface has to supercool first.
   A run that starts ice-free and never reaches the freezing point stays
   ice-free.


What happens in one step
========================

The ice runs after the ocean's dynamics and tracer step, on the ocean's
final state for that step. Two parts run on **every** outer step:

1. **Frazil banking.** Any supercooling in the ocean surface is clamped
   back to the freezing point and the heat deficit is banked for the ice.
2. **EVP dynamics** (when ``dynamics = .true.``), followed immediately by
   the ice-to-ocean **stress blend** that the next ocean momentum step
   feels.

The rest runs only on **thermodynamic steps**, the cadence set by
``&ocean_vmix_nml dt_therm_ratio`` (every outer step at the default 1), and
only when ``&ocean_thermo_nml enable_thermodynamics`` is on. The order is
fixed and load-bearing:

3. Category **transport** (when ``transport = .true.``).
4. The **atmospheric forcing** seam is filled: from the scalar restoring
   knobs, or from a prescribed atmosphere (see `Atmospheric forcing`_).
5. The ocean-to-ice **basal heat flux**.
6. **Frazil uptake**: the banked heat is spent as new ice.
7. The **column thermodynamics** (conduction, growth, melt, snow).
8. The open-water share of **snowfall** is delivered to the ocean.
9. The **brine, heat and shortwave couplers** write the ocean's surface
   salt and heat fluxes.
10. The **ice-thickness distribution** is restored.

The fluxes written in step 9 are applied to the ocean over the *next*
thermodynamic window, and the stress blend in step 2 is felt on the *next*
ocean momentum step. The coupling is therefore **lagged by one step** (the
model's convention for frazil and MEKE as well). A fresh run's first step
drives the ocean with the pure wind stress.

With ``enable_thermodynamics = .false.`` only steps 1 and 2 run: frazil is
banked but never spent.


Physics
=======

Frazil
------

After each outer step, wherever the surface layer is colder than its
freezing point :math:`T_f(S)`, the layer is raised back to :math:`T_f` and
the heat spent is added to a per-column frazil bank,

.. math::

   \Delta E \;=\; \rho\, c_p\, h_{\mathrm{sfc}}\,
     \bigl(T_f - T\bigr)^{+} \qquad [\mathrm{J\,m^{-2}}] ,

which the ice model later spends as latent heat of new ice. Energy is
banked, never discarded, and the ocean heat budget includes the bank so it
still closes with ice on. The freezing point is the ocean equation of
state's (``&ocean_eos_nml tfreeze_set`` / ``tfreeze_form``), evaluated at
zero pressure by default.

The "surface layer" here is the **first live layer** counting down from
the top, not literally ``k = nz``: if the top layer has collapsed to a
vanished thickness at run time, frazil (and the frazil uptake and basal
flux below) read and write the first layer that still holds water.

Three knobs change the frazil scheme, following ACCESS-OM3:

* ``frazil_column`` sweeps the whole column from the bed upward with a
  rising heat carry instead of clamping the surface only. Cold deep layers
  freeze, warmer layers above re-melt the carry, and only what survives to
  the top reaches the bank.
* ``pressure_dependent_frazil`` evaluates :math:`T_f` at each layer's
  hydrostatic centre pressure :math:`g \rho_0 z`. It needs
  ``frazil_column``.
* ``reclaim_frazil`` (default ``.true.``, used only with
  ``frazil_column``) spends the standing bank on a warm surface layer
  before the sweep.

On thermodynamic steps the bank is spent as new ice in the thinnest
category, annexing open water if needed. The salt the ice does not keep is
rejected to the ocean as brine (see `Coupling to the ocean`_).

Column thermodynamics
---------------------

Each category with ice runs a one-dimensional conduction and resize step
ported from SIS2: the **Winton (2000)** two-layer formulation, with
enthalpy as the prognostic. It uses one simplification, a brine heat
capacity equal to the ice heat capacity, which turns every
temperature–enthalpy inversion into an exact closed-form quadratic with no
iteration. The Stefan melt problem is reproduced to under 1 %.

The column step absorbs shortwave with the SIS2 **CSIM4** optics (snow and
ice albedo with a melting reduction, and Beer's-law penetration into the
ice and through to the ocean; the delta-Eddington branch is not ported),
solves conduction against the surface energy balance, then grows or melts
ice at the top and base and rebalances the layers to equal mass.

The ocean warms the ice from below through the **basal heat flux**,

.. math::

   F_b \;=\; \rho\, c_p\, \max\!\bigl(0,\; T_{\mathrm{sfc}} - T_f\bigr)\,
     \frac{h_{\mathrm{sfc}}}{\Delta t_{\mathrm{therm}}} ,

the full heat available above freezing in the surface layer, with no cap
(the column's melt is self-limiting and spills any excess back to the
ocean). It is zero on ice-free columns.

**Snow.** ``snowfall`` (kg/m²/s) is a uniform frozen-precipitation rate.
On ice it accumulates as snow. The share that falls on open water melts on
contact: its latent heat and its freshening are delivered to the ocean.
New snow is added after the optics and conduction of the same window, so
its albedo and insulation take effect on the next window.

**Snow-ice flooding** (``snow_ice``) converts snow to ice when its weight
pushes the snow–ice interface below the waterline (an Archimedes
freeboard step, SIS2's ``SN2IC``). It exchanges nothing with the ocean.
The converted ice keeps the snow's enthalpy and zero salinity, so it is
somewhat too cold and too fresh.

Ice-thickness distribution
--------------------------

With ``ncat > 1`` the ice is binned by thickness. After the thermodynamics
each thermodynamic step, a category whose mean thickness has crossed a bin
edge moves **all** of its area and mass into the neighbouring category,
merged by an area-weighted average. This is the scheme SIS2 actually ships
(not a Lipscomb (2001) linear remap). The thinnest category instead
shrinks in area when it thins below its lower edge.

The category lower edges default to SIS2's table,
1e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 m, extended by constant width
when ``ncat`` needs more. The top category has no upper bound. ``hlim``
overrides the table: give between 2 and ``ncat + 1`` strictly increasing
values with the first one positive, and the rest are extrapolated by
constant width, for example ``hlim = 1.0e-10, 0.2, 0.5, 1.0, 2.0``.

Transport
---------

With ``transport = .true.`` the ice and snow are advected between columns
on each thermodynamic step, in ``adv_substeps`` sub-iterations. The scheme
is SIS2's default velocity-driven category transport: a PPM reconstruction
of the category-summed mass gives one face flux, which is split between
categories in proportion to their share, with snow on its own PPM solve.
Enthalpy and salinity ride along upwind. A ``compress_ice`` step then
removes any over-full cell by compacting categories; this is area
compaction, not ridging. ``roll_factor`` sets SIS2's thin-ice rolling
floor (0 disables it).

The velocity the transport uses comes from the EVP dynamics when
``dynamics = .true.``. Without dynamics it is **sampled from the ocean's
top-layer face velocities**, so the ice drifts with the surface current.
A conservation or positivity failure in the transport stops the run.

EVP dynamics
------------

``dynamics = .true.`` turns on the C-grid **elastic-viscous-plastic (EVP)**
rheology of Hunke and Dukowicz (1997), ported from SIS2's C-grid solver,
with the Hibler (1979) ice strength and elliptical yield curve. It runs on
every outer step, independently of the thermodynamic cadence, in
``evp_sub_steps`` subcycles. The ice feels the wind stress, Coriolis, its
internal stress and a quadratic ice-ocean drag against the ocean's
**top-layer** (``k = nz``) velocity.

The strength and rheology knobs are SIS2's: ``p0`` (ice-strength
:math:`P^*`), ``c0`` (:math:`C^*`), ``ec`` (yield-curve axis ratio, 0 for a
cavitating fluid), ``cdw`` and ``rho_ocean`` (the drag), ``del_sh_min_scale``
(the viscosity floor) and ``tdamp`` (elastic damping time: positive in
seconds, 0 for :math:`0.2\,\Delta t`, negative as a fraction of the outer
step). An ice-free cell simply has zero strength; there is no ice-edge
special case.

Three options default **off** in Roundabout to keep earlier results
bit-identical, although SIS2 turns the first two on:

* ``cfl_trunc`` clips the final ice velocity to a transport CFL bound
  (SIS2 ``CFL_TRUNCATE``, 0.5 there). Each clip is counted and logged as a
  warning, never an abort. ``cfl_trunc_dyn_its`` also clips on every
  subcycle.
* ``project_ci`` projects the ice concentration forward within the
  subcycle loop and recomputes the strength (SIS2
  ``PROJECT_ICE_CONCENTRATION``).
* ``a_face_stress`` weights both the wind stress and the ice-ocean drag in
  the ice momentum balance by the face ice concentration. This makes the
  ice–ocean momentum budget close at every fractional cover. Without it,
  the default form leaks :math:`(1-a)(\tau_a - \tau_{io})` per face during
  transients at partial cover. SIS2 itself weights neither term.

Running ``dynamics`` without ``transport`` is accepted with a warning: the
velocity has no CFL backstop in that mode, and it is unvalidated.

The EVP solver has no sea-surface-tilt term (the ice assumes a flat sea
surface), and landfast ice is not ported.


Atmospheric forcing
===================

The ice surface sees the atmosphere through a linearised surface energy
balance, :math:`SF(T) = sf_0 + \tfrac{\partial SF}{\partial T}\,T` (upward
positive), plus a downwelling shortwave and a snowfall rate. Two fillers
supply it.

**Scalar restoring (the default).** A uniform, time-invariant slab
atmosphere:

.. math::

   SF(T) \;=\; \lambda\,\bigl(T - T_{\mathrm{air}}\bigr) ,

with :math:`\lambda` = ``restore_lambda`` (W/m²/K) and
:math:`T_{\mathrm{air}}` = ``air_temp`` (°C); ``sw_down`` and ``snowfall``
are spread uniformly. With ``restore_lambda = 0`` (the default) the column
feels only shortwave and the basal flux. A uniform, constant forcing cannot
produce a seasonal cycle.

**A prescribed atmosphere.** When ``&ocean_atm_forcing_nml enable`` is on
(together with ``&ocean_bulk_flux_nml enable``), the ice seam is filled
from the JRA55-do fields instead, and the four scalar knobs above are not
used. The surface energy balance is linearised at 0 °C with the Large and
Yeager (2004, 2009) constant over-ice transfer coefficients
(:math:`C_D = C_H = C_E = 1.63\times10^{-3}`, no stability iteration);
the shortwave and snowfall are the ``rsds`` and ``prsn`` fields, and the
longwave emissivity is ``&ocean_bulk_flux_nml emissivity``. The ocean's own
bulk fluxes are then weighted by the open-water fraction. This path needs a
non-cartesian grid and is single-rank for now.

**Wind stress.** In either case the ice feels a snapshot of the ocean's
wind stress taken once at configure. It does not follow a time-varying
wind, and there is no ice-specific bulk drag law.


Coupling to the ocean
=====================

The ice changes the ocean surface fluxes; the ocean's own interior
equations are untouched.

**Salt.** Freezing seawater of salinity :math:`S` into ice that keeps a
bulk salinity of 4 psu rejects the difference as brine; melting reverses
it. The rate is applied as a virtual salt flux on the surface
layer,

.. math::

   Q_{\mathrm{salt}} \;=\; Q_{\mathrm{salt}}^{\mathrm{const}}
     \;+\; \frac{m_{\mathrm{net}}\,\bigl(S - S_{\mathrm{ice}}\bigr)}
                {\Delta t_{\mathrm{therm}}} ,

with :math:`m_{\mathrm{net}}` the net mass frozen (positive) or melted. The
ocean's water mass is not changed: like MOM6 and SIS2's default, the ocean
is volume-conserving with a virtual salt flux.

**Heat.** The column's growth and melt exchange, net of the basal flux, is
added to the surface heat flux, and the shortwave transmitted through the
ice is added once on top. The frazil latent heat is not added here: the
frazil clamp already credited it to the ocean.

**Momentum** (only with ``dynamics = .true.``). On each face the ocean sees
a blend of the wind and the ice drag, weighted by the face ice
concentration :math:`a`:

.. math::

   \tau \;=\; (1 - a)\,\tau_a \;+\; a\,\tau_{io},
   \qquad a_u = \tfrac12\bigl(c_{i-1,j} + c_{i,j}\bigr) ,

and the KPP and EPBL friction velocity is re-derived from this blended
stress. Without ``dynamics`` the ocean sees the unmodified wind stress.

The constant ``&ocean_thermo_nml`` surface fluxes remain in
:math:`Q_{\mathrm{heat}}` and :math:`Q_{\mathrm{salt}}`; the ice terms are
added to them and are not weighted by ice cover. With
``&ocean_forcing_nml enable_components = .true.`` the same numbers arrive
through the surface-flux component set instead.

**Which layer.** The frazil, frazil-uptake and basal-flux kernels, and the
surface flux deposits, use the first live layer from the top. The ocean
velocity the ice drag and the transport sampler read is still the
``k = nz`` slice, a known gap on columns whose top layer has collapsed.

``bound_salinity`` (off by default) floors the top-layer salinity at
``min_salinity`` after the ice salt exchange. It is **non-conserving**: the
salt it creates is not in the salt budget.

Under ``&ocean_eos_nml ts_convention = "cons_abs"`` the ice keeps
potential temperature and practical salinity internally and the
conversion is made once at the coupler boundary.


Initial conditions
==================

By default the run starts ice-free. ``&ocean_ice_ic_nml`` seeds an analytic
pack before the first step:

.. code-block:: fortran

   &ocean_ice_ic_nml
      conc_config = "uniform"   ! "zero" (default), "uniform" or "latitudes"
      h_ice       = 2.0         ! m, required > 0
      conc        = 1.0         ! "uniform" only
   /

``"latitudes"`` covers everything poleward of ``arctic_edge`` and
``antarctic_edge`` (both default ±91°, i.e. no ice) and needs a
non-cartesian grid. The pack is binned straight into the thickness
categories, so the first ITD restore is a no-op. ``t_ice`` must be below
the freezing point at ``s_ice``, and snow needs ice under it. With
``ncat = 1`` a uniform concentration must be exactly 1. There is no file
reader. On a warm restart the ice IC is skipped (the restart replaces it),
with a warning.


Namelist reference
==================

``&ocean_ice_nml``:

=============================  ===========  ========================================
Knob                           Default      Meaning
=============================  ===========  ========================================
``enable``                     ``.false.``  Master switch for the whole component.
``ncat``                       ``5``        Thickness categories (category 0 is open water); 1 is the lumped mode.
``hlim``                       unset        Category lower edges (m); unset uses the SIS2 table.
``nk_ice``                     ``2``        Ice layers per category (2 = Winton two-layer).
``air_temp``                   ``0.0``      Slab-atmosphere air temperature (°C), scalar forcing.
``restore_lambda``             ``0.0``      Restoring coefficient :math:`\partial SF/\partial T` (W/m²/K).
``sw_down``                    ``0.0``      Downwelling shortwave onto the ice (W/m²).
``snowfall``                   ``0.0``      Uniform snowfall rate (kg/m²/s).
``snow_ice``                   ``.false.``  Snow-ice flooding.
``frazil_column``              ``.false.``  Full-column frazil sweep (OM3 ``FRAZIL``).
``pressure_dependent_frazil``  ``.false.``  Freezing point at layer pressure; needs ``frazil_column``.
``reclaim_frazil``             ``.true.``   Spend the bank on a warm surface; only with ``frazil_column``.
``bound_salinity``             ``.false.``  Floor the top-layer salinity (non-conserving).
``min_salinity``               ``0.0``      Floor for ``bound_salinity`` (g/kg).
``transport``                  ``.false.``  Horizontal category transport.
``adv_substeps``               ``1``        Transport sub-iterations per call (SIS2 ``NSTEPS_ADV``).
``roll_factor``                ``1.0``      Thin-ice rolling floor (SIS2 ``SEA_ICE_ROLL_FACTOR``); 0 disables.
``dynamics``                   ``.false.``  C-grid EVP dynamics.
``evp_sub_steps``              ``432``      EVP subcycles per outer step (SIS2 ``NSTEPS_DYN``).
``p0``                         ``2.75e4``   Ice-strength :math:`P^*` (Pa).
``c0``                         ``20.0``     Ice-strength :math:`C^*`.
``ec``                         ``2.0``      Yield-curve axis ratio; 0 = cavitating fluid.
``cdw``                        ``3.24e-3``  Ice-ocean drag coefficient.
``rho_ocean``                  ``1030.0``   Drag reference density (kg/m³).
``del_sh_min_scale``           ``2.0``      Viscosity-floor scale.
``tdamp``                      ``-0.2``     Elastic damping time (>0 s, 0 auto, <0 fraction of the step).
``a_face_stress``              ``.false.``  Concentration-weighted wind and drag (momentum-conserving).
``cfl_trunc``                  ``0.0``      CFL ceiling on the final ice velocity; 0 disables.
``cfl_trunc_dyn_its``          ``.false.``  Also clip every subcycle; needs ``cfl_trunc > 0``.
``project_ci``                 ``.false.``  Project concentration within the subcycles.
=============================  ===========  ========================================

``&ocean_ice_ic_nml``: ``conc_config`` (``"zero"``), ``conc`` (0.0),
``h_ice`` (0.0 m), ``h_snow`` (0.0 m), ``t_ice`` (−4.0 °C), ``s_ice``
(4.0 psu), ``arctic_edge`` (91.0), ``antarctic_edge`` (−91.0).

The generated map in ``docs/generated_nml_knobs.md`` is the authority if
this table and the code ever disagree.

What the configuration check refuses
------------------------------------

Like the rest of the model, the sea ice fails loud at configure rather than
silently ignoring a knob. In particular:

* ``transport`` needs ``enable``, ``ncat > 1`` and
  ``&ocean_thermo_nml enable_thermodynamics``.
* ``dynamics`` needs ``enable`` and every ``&ocean_bc_nml`` edge to be
  ``"wall"`` or ``"periodic"`` (the north edge may also be
  ``"tripolar_fold"``, which needs ``nghost >= 2``). Open, tidal, sponge,
  clamped and Chapman edges are refused. ``p0``, ``c0``, ``cdw``,
  ``rho_ocean`` and ``del_sh_min_scale`` must be positive, ``ec``
  non-negative.
* ``a_face_stress``, ``cfl_trunc > 0``, ``cfl_trunc_dyn_its`` and
  ``project_ci`` need ``dynamics``.
* ``snowfall > 0`` and ``snow_ice`` need ``enable`` and
  ``enable_thermodynamics``; ``frazil_column`` and ``hlim`` need
  ``enable``.
* Sea ice is mutually exclusive with the ice-shelf cavity
  (``&ocean_cavity_dyn_nml``, ``&ocean_cavity_melt_nml``), with
  ``&ocean_psurf_nml in_eos``, and with the pseudo-salt tracer
  (``&ocean_tracers_nml enable_pseudo_salt``).

``./rdb --validate-only case.nml`` lists every refusal without running.


Parallel runs and restarts
==========================

The sea ice runs on the ocean's domain decomposition and halo width. Its
category state, flux fields, EVP velocity and the blended surface stress
are halo-exchanged, and under a tripolar north fold they cross the fold
through the same machinery as the ocean prognostics, on any process grid.
A decomposed run is **bit-identical** to the single-rank run, including
under a tripolar fold split east–west; the decomposition test covers sea
ice with and without EVP and transport, with a mid-run checkpoint and
resume.

Restarts are bitwise, on one rank and decomposed, including the EVP
stresses and the exact blended surface stress. As for the ocean, a restart
must resume on the same decomposition.


Diagnostics
===========

With ``&ocean_ice_nml enable`` the default output gains a set of
two-dimensional sea-ice fields named after CICE's history fields, so
notebooks written for ACCESS-OM3 output read the same names. They are
time means over the output interval. Add or remove fields with
``&ocean_diag_nml diags`` (``name`` to add, ``name:off`` to drop); the full
catalogue and the selection syntax are on :doc:`diagnostics`.

=============  ======  ============  ========================================
Name           Units   Registered    Meaning
=============  ======  ============  ========================================
``aice``       1       ice           Ice area fraction (concentration).
``hi``         m       ice           Grid-mean (effective) ice thickness: ice mass over ice density.
``hs``         m       ice           Grid-mean snow thickness.
``Tsfc``       °C      ice           Ice-area-weighted surface temperature; open water counts at the freezing point.
``congel``     cm/day  ice + thermo  Basal (congelation) growth rate.
``frazil``     cm/day  ice + thermo  Frazil growth rate.
``snoice``     cm/day  ice + thermo  Snow-ice formation rate.
``meltt``      cm/day  ice + thermo  Top ice melt rate (snow melt excluded).
``meltb``      cm/day  ice + thermo  Basal ice melt rate.
``meltl``      cm/day  ice + thermo  Lateral melt rate. Always zero: there is no lateral melt.
``dvidtt``     cm/day  ice + thermo  Thermodynamic volume tendency, ``congel + frazil + snoice - meltt - meltb``.
``ice_conc``   1       opt-in        Alias of ``aice`` (a former default name).
``ice_thick``  m       opt-in        Alias of ``hi`` (a former default name).
``ice_speed``  m/s     opt-in        Ice drift speed at the cell centre.
``ice_u``      m/s     opt-in        Ice velocity, grid x component, at the cell centre.
``ice_v``      m/s     opt-in        Ice velocity, grid y component, at the cell centre.
``uvel``       m/s     opt-in        Alias of ``ice_u``.
``vvel``       m/s     opt-in        Alias of ``ice_v``.
=============  ======  ============  ========================================

"ice" means registered by default whenever ``&ocean_ice_nml enable`` is
on; "ice + thermo" also needs ``&ocean_thermo_nml enable_thermodynamics``;
"opt-in" fields appear only when named in ``diags``. All are 2-D. The rates are the masses
formed or lost over each thermodynamic window divided by the window length
and the ice density, then averaged over the output interval. The
dynamical volume tendency (CICE ``dvidtd``) is not provided.

.. note::

   ``ice_u``, ``ice_v``, ``uvel`` and ``vvel`` are the face velocities
   averaged to the cell centre along the model grid's own axes. CICE's
   ``uvel``/``vvel`` are B-grid corner values, so the staggering differs.

The console status line also gains an ice line with the ocean-area mean
concentration and thickness, for example:

.. code-block:: text

       Ice  : conc  0.4213  thick    0.8127 m  (ocean-area mean)


Limitations
===========

* **No ridging.** Under convergence ``compress_ice`` thickens ice within
  its own category (SIS2 with ``DO_RIDGING = .false.``), so the thickness
  distribution and the EVP strength are biased in convergent regimes.
* **No melt ponds.** The CSIM4 melting-albedo ramp stands in, so summer
  melt is biased low and the pond albedo feedback is absent.
* **No lateral melt or floe-size effects.** Marginal-ice-zone retreat is
  biased.
* **No freshwater or mass exchange.** Salt is exchanged as a virtual flux;
  the ice puts no load on the sea surface.
* **Simple forcing.** Without a prescribed atmosphere the forcing is a
  uniform, constant slab. The wind stress the ice feels is a configure-time
  snapshot in every case.
* **Momentum leak at partial cover** unless ``a_face_stress`` is on.
* **The ice drag reads the ocean's top-layer velocity** (``k = nz``)
  even when that layer has collapsed.
* **Analytic initial conditions only.**

``docs/CAPABILITIES_AND_LIMITATIONS.md`` (section "Sea ice") and
``docs/CLOSURE_MATRIX.md`` (section "Sea ice") carry the per-feature status
and the tests behind each piece.


Provenance and references
=========================

The sea ice follows SIS2, and its numerics are the SIS2 developers' work;
Roundabout's contribution is the port. Where SIS2 has a runtime parameter
for a knob, the table above names it. Deliberate divergences from SIS2 are
documented in the module docstrings, rendered in the API documentation.

The physics comes from:

* Winton (2000) — the two-layer column thermodynamics and resize step.
* Hunke and Dukowicz (1997) — the elastic-viscous-plastic rheology.
* Hibler (1979) — the ice strength, elliptical yield curve and momentum
  balance.
* Thorndike, Rothrock, Maykut and Colony (1975), *J. Geophys. Res.* 80,
  4501-4513, and Bitz, Holland, Weaver and Eby (2001), *J. Geophys. Res.*
  106, 2441-2463 — the thickness-distribution discretisation.
* Hunke et al. (2015), the CICE v5 user guide — melt-energy accounting and
  layer rebalancing.
* Leppäranta (1983), *Nordic Hydrology* 14, 59-70, and Fichefet and
  Morales Maqueda (1997), *J. Geophys. Res.* 102, 12609-12646 — snow-ice
  flooding.
* Large and Yeager (2004, 2009) — the over-ice bulk transfer coefficients.

The frazil options and the salinity floor follow ACCESS-OM3. See
:doc:`references` for the acknowledgement of MOM6 and SIS2 and their
licensing, and :doc:`parameterizations` for the ocean-side closures the
ice feeds (KPP, EPBL, shortwave penetration).
