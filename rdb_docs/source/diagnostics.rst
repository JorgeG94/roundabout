.. _diagnostics:

-----------
Diagnostics
-----------

.. contents::
   :local:

This page lists every diagnostic field Roundabout can write and explains how
to choose them. :ref:`running_simulations` has the short version, under
"Output and diagnostics". This page is the full catalogue.

Roundabout reports in two places:

* the **diagnostic stream**: NetCDF files, one per rank, configured by
  ``&ocean_diag_nml``;
* the **console status line**: conservation and stability numbers printed to
  the log every ``&logging_nml status_interval``.

There are also some log-only forensic probes, covered at the end.


Quick start
===========

.. code-block:: fortran

   &time_nml
      time_unit = "day"
   /

   &ocean_diag_nml
      enabled  = .true.
      filename = "tasman"
      dt_out   = 1.0                     ! in time_unit: one frame per day
      diags    = "KE:off  vorticity_z  temperature:6h:max  transport_x:mean"
   /

   &output_nml
      output_dir      = "./out_tasman"
      compress_output = .true.
      compress_level  = 1
   /

This writes the default set without ``KE``. It adds a daily snapshot of
``vorticity_z``, a 6-hourly maximum of ``temperature`` and a daily mean of
``transport_x``. Rank ``r`` writes ``./out_tasman/tasman_rank_<r>.nc``, with
the rank number zero-padded to six digits.

.. note::

   Diagnostics are **on by default** (``enabled = .true.``), and the default
   cadence is ``dt_out = 3600`` **in your** ``time_unit``. With
   ``time_unit = "day"`` that is one frame every 3600 days. Set ``dt_out``
   explicitly. The model warns when ``dt_out`` is longer than ``t_end``.


Choosing what is written
========================

The default set
---------------

When ``diags`` is empty (the default), the model registers the **canonical
set**. It always includes:

``SSH``, ``temperature``, ``salinity``, ``u``, ``v``, ``KE``

It also includes every canonical field whose feature is switched on:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Fields added
     - When this is on
   * - ``age``
     - ``&ocean_tracers_nml enable_ideal_age``
   * - ``pseudo_salt``, ``pseudo_salt_diff``
     - ``&ocean_tracers_nml enable_pseudo_salt``
   * - ``MLD_EPBL``, ``Kd_EPBL``
     - ``&ocean_epbl_nml enable``
   * - ``Kd_KSHEAR``
     - ``&ocean_kappa_shear_nml enable``
   * - ``Kd_TIDAL``
     - ``&ocean_tidal_mixing_nml enable``
   * - ``Kd_BBL``
     - ``&ocean_bbl_mixing_nml enable``
   * - ``aice``, ``hi``, ``hs``, ``Tsfc``
     - ``&ocean_ice_nml enable``
   * - ``congel``, ``frazil``, ``snoice``, ``meltt``, ``meltb``, ``meltl``,
       ``dvidtt``
     - ``&ocean_ice_nml enable`` and ``&ocean_thermo_nml
       enable_thermodynamics`` (default on)
   * - ``atm_tas`` … ``atm_cosz`` (ten fields)
     - ``&ocean_atm_forcing_nml enable`` (NetCDF build)

Every other field in the catalogue is **opt-in**: it is written only if you
name it in ``diags``. That covers the derived fields, the ``ice_conc`` and
``ice_thick`` aliases, and the sea-ice drift fields.

The ``diags`` selection list
----------------------------

``&ocean_diag_nml diags`` is a single string of **entries**. Separate entries
with spaces or commas. Each entry has the form::

   name[:attribute[:attribute...]]

* ``name`` is a catalogue name. It is **case-sensitive**: ``SSH`` and ``KE``
  are upper case, ``vorticity_z`` is lower case.
* Attributes are case-insensitive and can come in any order. Each one sets one
  of the following:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Effect
   * - ``off``
     - Do not register this field. Use it to drop a field from the default
       set. On an opt-in name that was not requested, it does nothing.
   * - ``instant`` / ``mean`` / ``max`` / ``min`` / ``integral``
     - Time operation (see `Time operations and cadence`_).
   * - ``<n>s`` / ``<n>m`` / ``<n>h`` / ``<n>d``
     - Cadence for this field only, overriding ``dt_out``. ``<n>`` is a
       positive **integer**, and ``m`` means **minutes**. The cadence is in
       absolute units and does not depend on ``time_unit``. ``1.5d`` is
       rejected; write ``36h``.
   * - ``layer`` / ``z`` / ``zstar`` (or ``z*``) / ``sigma`` / ``density``
       (or ``rho``)
     - Output vertical grid for this field, overriding ``vgrid``. 2-D fields
       ignore it. Note the per-field spelling ``z``: the global ``vgrid``
       knob spells the same grid ``z_fixed``.

Examples:

.. code-block:: fortran

   diags = "SSH:1d  temperature:1d  u:off  v:off  KE:off  vorticity_z:1d"
   diags = "temperature:z:mean  salinity:z:mean  umo:density:1d"
   diags = "aice:6h  hi:6h  ice_speed:6h:mean"

What happens to each entry:

* A **canonical** name changes that field's time operation, cadence or grid,
  or removes it with ``off``.
* A **derived** name (one not in the canonical set) registers that field.
  Derived fields default to ``instant``. Unless you give attributes, they use
  ``dt_out`` and ``vgrid``.
* A canonical name whose feature is off is **not** registered. The run goes
  ahead and logs a warning that names the switch to check.
* A derived name whose prerequisite is off (the ice-shelf cavity fields) is a
  **configure-time error**. A plane of missing values could otherwise end up
  in a figure looking like a real result.
* An unknown name, or an unknown attribute, stops the run at configure time.
  The error lists the valid names.

Time operations and cadence
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 14 50 36

   * - Operation
     - What is written each frame
     - Time stamp
   * - ``instant``
     - A snapshot, taken only on the step the cadence fires.
     - End of the window
   * - ``mean``
     - Average over the window, weighted by ``dt``. Sampled every step.
     - Centre of the window
   * - ``max`` / ``min``
     - Running maximum or minimum for each cell over the window. Sampled
       every step.
     - End of the window
   * - ``integral``
     - :math:`\sum(\text{sample} \cdot \Delta t)` over the window, not divided
       by the window length. Use it to close a budget. The ``units``
       attribute is **not** multiplied by seconds.
     - Centre of the window

Default operations:

* ``SSH`` is ``instant``.
* Every other canonical field is ``mean``.
* Every derived field is ``instant``.

Every field keeps its own window. A window starts at registration and closes
when the time accumulated in it reaches the field's cadence. The first frame
is therefore at :math:`t = \text{cadence}`, and the initial condition is not
written.

Fields with a time operation other than ``instant`` are filled and folded in
**every step** on the device. That costs GPU time even between frames.
``instant`` fields cost nothing between frames.


Output vertical grids
=====================

3-D fields are computed on the model layers. You can write them on that grid
or remap them conservatively onto another vertical grid when they are output.
For a ``mean``, ``max``, ``min`` or ``integral`` field, the remap runs every
step, so the average is formed on the output grid.

.. list-table::
   :header-rows: 1
   :widths: 14 30 56

   * - Grid
     - Levels knob
     - Meaning of the levels
   * - ``layer``
     - None
     - Native model layers. **Index 1 is the bed and index nz is the
       surface** (the model's bottom-up convention).
   * - ``z_fixed`` (``:z``)
     - ``z_levels`` + ``n_z_levels``, which must be set
     - Bottom-interface depths in metres, positive down, shallow to deep,
       measured from the top of the water column. The surface interface at
       0 is implicit, so do not list it.
   * - ``sigma``
     - ``sigma_levels`` + ``n_sigma_levels``. Without them, ``nz`` uniform
       levels are used.
     - Cumulative fractions of the column depth, 0 to 1, shallow to deep.
   * - ``zstar``
     - ``zstar_levels`` + ``n_zstar_levels``. Without them, ``nz`` uniform
       levels are used.
     - Reference interface depths in metres, shallow to deep. Each column's
       grid is stretched by the ratio of the column depth to the deepest
       level. With uniform levels this grid is identical to ``sigma``. It only
       differs when you give a non-uniform (for example, fine-near-surface)
       reference.
   * - ``density``
     - ``rho_levels`` + ``n_rho_levels``, which must be set and must
       increase strictly
     - Target potential densities in kg m\ :sup:`-3`, light to dense. The
       potential density is taken at ``&vcoord_nml rho_ref_pressure``
       (default :math:`2\times10^7` Pa).

On every **remapped** grid, output index 1 is the **shallowest** (or
lightest) cell. That is the opposite of native layers. The file has no
vertical coordinate variable, so keep your level list with the output.

You can give at most ``NZ_STACK_MAX`` levels. This is a CMake option,
``RDB_NZ_STACK_MAX``, and defaults to 128.

The remap redistributes thickness-integrated fields (``umo``, ``vmo``) so the
column sum is preserved. It averages every other field. The reconstruction is
set by ``diag_remap_scheme``:

* ``pcm``, ``plm``, ``ppm`` (the default) and ``ppm_h4`` work.
* ``pqm`` is accepted by the namelist but has no diagnostic mapping. It
  silently falls back to PLM, so do not use it here.

Cells below the bottom of a shallow column read **0** by default. Set
``mask_vanished_layers = .true.`` to write the missing value there instead.
That applies to every grid except ``density``.


Catalogue
=========

Every field is written on the horizontal grid of the cell centres (T cells).
The **Shape** column uses these terms:

* **2-D**: one value per column.
* **layer**: one value per model layer, before any remap.
* **interface**: one value per layer, taken at the layer's **bottom**
  interface. The surface interface is dropped.
* **u face** / **v face**: the value on the cell's west / south face, written
  as it is. The eastmost / northmost face of the domain has no slot and is
  not written.

Other fields that live on faces or corners are averaged to the cell centre,
and the table says so.

Core ocean state (always registered)
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 10 10 64

   * - Name
     - Units
     - Shape
     - Meaning
   * - ``SSH``
     - m
     - 2-D
     - Sea-surface height: total column thickness minus bathymetry depth.
       Default ``instant``.
   * - ``temperature``
     - degC
     - layer
     - Potential temperature. Under ``&ocean_eos_nml ts_convention =
       "cons_abs"`` it is conservative temperature, with the CF metadata to
       match.
   * - ``salinity``
     - psu
     - layer
     - Practical salinity. Under ``ts_convention = "cons_abs"`` it is
       absolute salinity in ``g kg-1``.
   * - ``u``
     - m s-1
     - layer
     - Eastward velocity, the 2-point average of the u faces.
   * - ``v``
     - m s-1
     - layer
     - Northward velocity, the 2-point average of the v faces.
   * - ``KE``
     - m2 s-2
     - layer
     - Kinetic energy per unit mass,
       :math:`\tfrac14(u_W^2+u_E^2+v_S^2+v_N^2)`. This is the C-grid
       stencil the Coriolis kinetic-energy term uses.

Passive tracers
---------------

.. list-table::
   :header-rows: 1
   :widths: 16 8 8 46 22

   * - Name
     - Units
     - Shape
     - Meaning
     - Requires
   * - ``age``
     - s
     - layer
     - Ideal age of sea water.
     - ``&ocean_tracers_nml enable_ideal_age``
   * - ``pseudo_salt``
     - psu
     - layer
     - The pseudo-salt verification tracer (Shao 2016).
     - ``&ocean_tracers_nml enable_pseudo_salt``
   * - ``pseudo_salt_diff``
     - psu
     - layer
     - ``pseudo_salt`` minus salinity: how far the passive-tracer transport
       has drifted from the active-tracer transport.
     - ``&ocean_tracers_nml enable_pseudo_salt``

Vertical-mixing closures
------------------------

Each ``Kd_*`` field below is a **single closure's own contribution**. For the
total diffusivity that is actually applied to temperature, see
``Kd_interface`` in `Derived: transports, stress and total mixing`_.

.. list-table::
   :header-rows: 1
   :widths: 14 9 11 40 26

   * - Name
     - Units
     - Shape
     - Meaning
     - Requires
   * - ``MLD_EPBL``
     - m
     - 2-D
     - Depth of the EPBL active mixing layer.
     - ``&ocean_epbl_nml enable``
   * - ``Kd_EPBL``
     - m2 s-1
     - interface
     - EPBL diffusivity.
     - ``&ocean_epbl_nml enable``
   * - ``Kd_KSHEAR``
     - m2 s-1
     - interface
     - Kappa-shear diffusivity.
     - ``&ocean_kappa_shear_nml enable``
   * - ``Kd_TIDAL``
     - m2 s-1
     - interface
     - Internal-tide (tidal mixing) diffusivity.
     - ``&ocean_tidal_mixing_nml enable``
   * - ``Kd_BBL``
     - m2 s-1
     - interface
     - Law-of-the-wall bottom-boundary-layer diffusivity.
     - ``&ocean_bbl_mixing_nml enable``

Sea ice
-------

These fields need ``&ocean_ice_nml enable``. They follow CICE's history
names, and all are 2-D at the cell centre. See :doc:`sea_ice` for the model
behind them.

.. list-table::
   :header-rows: 1
   :widths: 13 10 49 28

   * - Name
     - Units
     - Meaning
     - Registered
   * - ``aice``
     - 1
     - Total ice concentration, 0 to 1.
     - Default when ice is on
   * - ``hi``
     - m
     - Grid-cell mean ice thickness: ice mass divided by ice density. This is
       the effective thickness, not the floe thickness, so the field has no
       CF ``standard_name``.
     - Default when ice is on
   * - ``hs``
     - m
     - Grid-cell mean snow thickness.
     - Default when ice is on
   * - ``Tsfc``
     - degC
     - Aggregate surface temperature: the ice-area-weighted mean of the skin
       temperatures for each category. The open-water share of the cell is
       counted at the seawater freezing point.
     - Default when ice is on
   * - ``ice_conc``
     - 1
     - Alias of ``aice``.
     - Only when named
   * - ``ice_thick``
     - m
     - Alias of ``hi``.
     - Only when named
   * - ``congel``
     - cm/day
     - Basal congelation growth rate.
     - Default, needs thermodynamics
   * - ``frazil``
     - cm/day
     - Frazil growth rate.
     - Default, needs thermodynamics
   * - ``snoice``
     - cm/day
     - Rate of snow-ice formation by flooding.
     - Default, needs thermodynamics
   * - ``meltt``
     - cm/day
     - Top melt rate of ice, not counting snow melt.
     - Default, needs thermodynamics
   * - ``meltb``
     - cm/day
     - Basal melt rate of ice.
     - Default, needs thermodynamics
   * - ``meltl``
     - cm/day
     - Lateral melt rate. The model has **no lateral-melt scheme**, so this is
       always 0. It exists so that tools expecting CICE names find it.
     - Default, needs thermodynamics
   * - ``dvidtt``
     - cm/day
     - Thermodynamic ice-volume tendency:
       ``congel + frazil + snoice - meltt - meltb``. This is not the dynamic
       tendency.
     - Default, needs thermodynamics

The rates are window masses from the most recent thermodynamic interval,
converted to cm/day. Land and halo cells are NaN.

The **ice drift** fields are opt-in derived fields. They are not gated: with
ice off they are written as zeros.

.. list-table::
   :header-rows: 1
   :widths: 14 10 76

   * - Name
     - Units
     - Meaning
   * - ``ice_speed``
     - m s-1
     - Ice drift speed :math:`|\mathbf{u}_{ice}|`, with face pairs averaged
       to the cell centre.
   * - ``ice_u``
     - m s-1
     - Eastward ice velocity, averaged to the cell centre.
   * - ``ice_v``
     - m s-1
     - Northward ice velocity, averaged to the cell centre.
   * - ``uvel``
     - m s-1
     - Alias of ``ice_u``. CICE's ``uvel`` is on B-grid corners, but this one
       is at the cell centre.
   * - ``vvel``
     - m s-1
     - Alias of ``ice_v``, staggered the same way as ``uvel``.

The console also prints a sea-ice line (see `Console status line`_).

Prescribed atmosphere
---------------------

These fields need ``&ocean_atm_forcing_nml enable`` and a NetCDF build. They
are the prescribed forcing fields copied through as they are, all 2-D.

.. list-table::
   :header-rows: 1
   :widths: 16 14 70

   * - Name
     - Units
     - Meaning
   * - ``atm_tas``
     - K
     - 2 m air temperature.
   * - ``atm_huss``
     - 1
     - 2 m specific humidity.
   * - ``atm_uas``
     - m s-1
     - 10 m eastward wind.
   * - ``atm_vas``
     - m s-1
     - 10 m northward wind.
   * - ``atm_psl``
     - Pa
     - Sea-level pressure.
   * - ``atm_rlds``
     - W m-2
     - Downwelling longwave radiation.
   * - ``atm_rsds``
     - W m-2
     - Downwelling shortwave radiation: the final field after any shaping by
       the solar zenith angle.
   * - ``atm_prra``
     - kg m-2 s-1
     - Rainfall flux.
   * - ``atm_prsn``
     - kg m-2 s-1
     - Snowfall flux.
   * - ``atm_cosz``
     - 1
     - Cosine of the solar zenith angle. It stays 0 under
       ``sw_interp = "linear"``, which never computes it.

Derived: ocean state and dynamics
---------------------------------

These fields are opt-in and need no other switch.

.. list-table::
   :header-rows: 1
   :widths: 16 10 10 64

   * - Name
     - Units
     - Shape
     - Meaning
   * - ``h_layer``
     - m
     - layer
     - Layer thickness.
   * - ``rho_layer``
     - kg m-3
     - layer
     - Layer in-situ density from the equation of state. Vanished layers are
       NaN.
   * - ``vorticity_z``
     - s-1
     - layer
     - Relative vorticity :math:`\zeta = \partial_x v - \partial_y u`. It is
       the model's own corner vorticity, including the free-slip or no-slip
       treatment, averaged from the four corners to the cell centre. Land
       cells carry the missing value.
   * - ``ke_total``
     - m3 s-2
     - 2-D
     - Depth-integrated kinetic energy, :math:`\sum_k \tfrac12 h_k
       (u_c^2+v_c^2)`.
   * - ``transport_x``
     - m2 s-1
     - 2-D
     - Depth-integrated eastward transport, :math:`\sum_k h_k u_c`, at the
       cell centre.
   * - ``transport_y``
     - m2 s-1
     - 2-D
     - Depth-integrated northward transport, at the cell centre.
   * - ``mld_density``
     - m
     - 2-D
     - Mixed-layer depth from a density threshold of 0.03 kg m\ :sup:`-3`
       relative to the top live layer (de Boyer Montégut 2004). A column with
       no crossing reports its full depth.
   * - ``temperature_pt``
     - degC
     - layer
     - Potential temperature :math:`\theta`. Under ``cons_abs`` it is
       converted from conservative temperature and absolute salinity.
       Otherwise it is a copy of ``temperature``.
   * - ``salinity_sp``
     - 1e-3
     - layer
     - Practical salinity. Under ``cons_abs`` it is converted from absolute
       salinity. Otherwise it is a copy of ``salinity``.

Derived: transports, stress and total mixing
--------------------------------------------

These fields are opt-in and need no other switch. Their names follow
MOM6/CMOR so that output can be compared with ACCESS-OM3.

.. list-table::
   :header-rows: 1
   :widths: 14 9 12 65

   * - Name
     - Units
     - Shape
     - Meaning
   * - ``umo``
     - kg s-1
     - layer, u face
     - Layer mass transport through the west face:
       :math:`\rho_0 \times` the volume flux that continuity uses. This is
       the flux from the **last** continuity call of the step, not a mean
       over the step's stages. When remapped, the column sum is preserved.
   * - ``vmo``
     - kg s-1
     - layer, v face
     - The same through the south face.
   * - ``tauuo``
     - N m-2
     - 2-D, u face
     - Zonal wind stress as applied to the momentum equation, after any
       masking by ice cover.
   * - ``tauvo``
     - N m-2
     - 2-D, v face
     - Meridional wind stress as applied.
   * - ``Kd_interface``
     - m2 s-1
     - interface
     - Total diffusivity for heat after every interior and boundary-layer
       closure has contributed (MOM6 ``difvho``).

Derived: ice-shelf cavity
-------------------------

All of these are 2-D. Outside the ice-shelf cover they are written as
missing values, so a domain mean of the output is a mean over the cavity.
Asking for one without its prerequisite is a configure-time error.

.. list-table::
   :header-rows: 1
   :widths: 19 10 45 26

   * - Name
     - Units
     - Meaning
     - Requires
   * - ``melt``
     - kg m-2 s-1
     - Basal melt mass flux, positive for melting.
     - ``&ocean_cavity_melt_nml enable``
   * - ``melt_m_per_yr``
     - m yr-1
     - Basal melt rate in the ISOMIP+ reporting unit: metres of ice-equivalent
       freshwater per year, using a freshwater density of
       1000 kg m\ :sup:`-3` and a 365-day year.
     - ``&ocean_cavity_melt_nml enable``
   * - ``thermal_driving``
     - degC
     - Far-field temperature above the in-situ freezing point at the
       interface pressure. A positive value drives melting.
     - ``&ocean_cavity_melt_nml enable``
   * - ``haline_driving``
     - g kg-1
     - Far-field salinity minus the interface salinity.
     - ``&ocean_cavity_melt_nml enable``
   * - ``tbdry``
     - degC
     - Interface temperature, on the liquidus at the interface salinity.
     - ``&ocean_cavity_melt_nml enable``
   * - ``sbdry``
     - g kg-1
     - Interface salinity.
     - ``&ocean_cavity_melt_nml enable``
   * - ``tfreeze_ib``
     - degC
     - In-situ freezing point of the far field at the ice base.
     - ``&ocean_cavity_melt_nml enable``
   * - ``exch_vel_t``
     - m s-1
     - Thermal exchange velocity the melt solve converged on.
     - ``&ocean_cavity_melt_nml enable``
   * - ``exch_vel_s``
     - m s-1
     - Haline exchange velocity the melt solve converged on.
     - ``&ocean_cavity_melt_nml enable``
   * - ``ustar_shelf``
     - m s-1
     - Friction velocity in the melt law. It drives no drag, and KPP and EPBL
       do not read it.
     - ``&ocean_cavity_melt_nml enable``
   * - ``cavity_melt_status``
     - 1
     - Melt-solver status code for each column. 0 is OK. Any other value
       means the column fell back to zero melt.
     - ``&ocean_cavity_melt_nml enable``
   * - ``z_draft``
     - m
     - Ice-shelf draft below the geoid, positive down.
     - ``&ocean_cavity_dyn_nml enable``
   * - ``water_column``
     - m
     - Thickness of the water column under the shelf.
     - ``&ocean_cavity_dyn_nml enable``

``&ocean_cavity_melt_nml`` itself requires ``&ocean_cavity_dyn_nml``.

Counting names, not distinct quantities, the catalogue has **37 canonical**
and **32 derived** names. Four of them are aliases: ``ice_conc``,
``ice_thick``, ``uvel`` and ``vvel``.


Output files
============

Layout
------

Each rank writes ``<output_dir>/<filename>_rank_NNNNNN.nc``. There is no
online gather. If ``output_dir`` does not exist, the model creates it, but
only one level deep: a missing parent directory is still an error.

Each file is CF-1.8 and is laid out as follows:

* Every field has **its own time axis**: an unlimited dimension and a
  variable named ``time_<name>``, with units ``seconds since simulation
  start``. Fields on different cadences are therefore independent.
* A 2-D field is ``(x, y, time_<name>)``. A 3-D field is
  ``(x, y, z_<name>, time_<name>)``.
* ``x`` and ``y`` are the rank's **full local extent, halo included**, with
  ``nghost`` ghost cells on each side. There are no ``x`` or ``y`` coordinate
  variables.
* Each data variable carries ``units``, ``long_name`` and, where one
  applies, ``standard_name``.
* In a run with more than one rank, the global attributes ``px``, ``py``,
  ``i_start``, ``j_start``, ``nx_global``, ``ny_global``, ``nx_local``,
  ``ny_local`` and ``nghost`` describe where the rank's tile sits.
* The file is synced after every frame, so frames already written survive
  an aborted run.

**Missing data.** Treat both of the following as missing:

* **NaN**: vanished layers in tracer and density fields, land and halo cells
  in the default sea-ice fields (``aice`` through ``dvidtt``), and cells
  outside the cover in the cavity fields;
* **1.0e20**: advertised as ``_FillValue`` and ``missing_value`` on
  ``vorticity_z``, on the cavity fields, and on remapped fields when
  ``mask_vanished_layers = .true.``.

Precision and compression
-------------------------

``output_precision = "single"`` writes the data variables as 32-bit floats
and halves the bytes per frame. Time axes, restarts and console totals stay
at working precision.

``&output_nml compress_output = .true.`` turns on deflate compression, with
``compress_level`` from 1 to 9. Level 1 is usually the best trade-off.

Merging per-rank files
----------------------

.. code-block:: bash

   python3 tools/merge_output.py ./out_tasman tasman_merged.nc --prefix tasman

The tool trims each rank's halo and places it in the global array. Pass
``--prefix`` with your ``filename``: the tool's default prefix is ``output``,
which does not match the diagnostic default ``ocean_diag``. The tool imports
``netCDF4`` and ``numpy``.

The log line for each frame
---------------------------

Every frame also prints one log line, so a run without NetCDF post-processing
still shows the range of each field::

   [diag] t=<time> <name> [<units>]  min=<v>  max=<v>  mean=<v>  missing=<n>/<total>

The min, max and mean cover finite cells only. The ``missing=`` suffix
appears only when some cells were excluded.


Console status line
===================

Separately from the NetCDF stream, rank 0 prints a status block every
``&logging_nml status_interval``. The interval is in ``time_unit``, and 0
means every 100 steps. The block contains:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Item
     - Meaning
   * - ``Day``, ``step``
     - Simulation time in days, and the step count.
   * - ``En``
     - Total kinetic energy per unit mass, in m\ :sup:`2` s\ :sup:`-2`. This
       is comparable to MOM6's ``En``.
   * - ``MaxCFL``
     - The largest advective CFL number, :math:`|u|\Delta t/\Delta x +
       |v|\Delta t/\Delta y`.
   * - ``Mass``
     - Total mass :math:`\sum h\,A\,\rho_0`, in kg.
   * - ``Salt``, ``Temp``
     - Volume-mean salinity and temperature. Shown when thermodynamics is on.
   * - ``Mass :`` / ``Salt :`` / ``Heat :`` lines
     - The total, and ``Error``: the change since :math:`t=0` relative to the
       initial value. When that quantity's budget is tracked, the cumulative
       boundary outflux ``out`` and surface source ``src`` are also shown, and
       ``Error`` becomes the closed-budget residual. It should stay near
       round-off.
   * - ``Age :``
     - Volume-mean ideal age, in days. Shown when ideal age is on.
   * - ``En :``
     - Absolute kinetic energy in J, and its growth relative to the initial
       KE. This tracks spin-up and is not a conservation check.
   * - ``Ice :``
     - Ocean-area mean ice concentration and thickness. Shown when ice is on.

The run **aborts** if any console total is NaN, or if ``MaxCFL`` exceeds 0.9.
The last good line is printed before the abort.

With ``&ocean_diag_nml reproducing_sums = .true.`` (the default), the totals
and budget terms are summed in fixed point. The console then prints the
**same digits on any rank count**, and a non-finite value anywhere shows up
as NaN rather than a plausible-looking number. Set it to ``.false.`` to use
ordinary floating-point sums.


Forensic probes (log only)
==========================

``&ocean_debug_nml`` holds three heavy debugging probes, all off by default.
They print greppable rows to the log and write nothing to NetCDF:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Knob
     - What it prints
   * - ``budget``
     - ``BT_BUDGET`` rows: for each region, the power each slow tendency
       term (PGF, Coriolis/advection, viscosity, bottom drag, surface stress)
       puts into the barotropic mode.
   * - ``ke_attr``
     - ``KE_ATTR`` rows: the change in layer KE across each segment of the
       split stage, for the whole domain, a rim band and an east band.
       Limit the steps sampled with ``ke_attr_start_step`` and
       ``ke_attr_end_step``. It is a single-rank tool.
   * - ``chksum``
     - ``CHKSUM`` rows: sum, min, max, non-finite count and a
       decomposition-invariant bit count for each field at each phase.
       Diff two logs to find the first operator where they diverge. Limit
       the steps sampled with ``chksum_start_step`` and ``chksum_end_step``.
       ``chksum_interior`` restricts the reduction to physical cells, which
       makes the bit count comparable between 1 rank and N ranks.


Knob reference
==============

``&ocean_diag_nml``:

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Knob
     - Default
     - Meaning
   * - ``enabled``
     - ``.true.``
     - Master switch for the diagnostic stream.
   * - ``filename``
     - ``"ocean_diag"``
     - File basename. The suffix ``_rank_NNNNNN.nc`` is appended.
   * - ``dt_out``
     - ``3600``
     - Default cadence, **in** ``time_unit``.
   * - ``diags``
     - ``""``
     - The selection list (see `The diags selection list`_). Empty means the
       canonical set.
   * - ``vgrid``
     - ``"layer"``
     - Default output grid for layered fields: ``layer``, ``z_fixed``,
       ``sigma``, ``zstar`` or ``density``.
   * - ``z_levels``, ``n_z_levels``
     - none, ``0``
     - Depths in metres, positive down. Required for ``z_fixed``.
   * - ``sigma_levels``, ``n_sigma_levels``
     - none, ``0``
     - Fractions of the column depth. If unset, ``nz`` uniform levels are
       used.
   * - ``zstar_levels``, ``n_zstar_levels``
     - none, ``0``
     - Reference depths in metres. If unset, ``nz`` uniform levels are used.
   * - ``rho_levels``, ``n_rho_levels``
     - none, ``0``
     - Densities in kg m\ :sup:`-3`, strictly increasing. Required for
       ``density``.
   * - ``mask_vanished_layers``
     - ``.false.``
     - Write the missing value, instead of 0, into remapped cells below the
       bottom.
   * - ``diag_remap_scheme``
     - ``"ppm"``
     - ``pcm``, ``plm``, ``ppm`` or ``ppm_h4``.
   * - ``output_precision``
     - ``"double"``
     - ``"single"`` halves the bytes per frame.
   * - ``reproducing_sums``
     - ``.true.``
     - Fixed-point sums for the console totals.

Related knobs in other groups:

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Knob
     - Default
     - Meaning
   * - ``&output_nml output_dir``
     - ``"./output"``
     - Directory for the per-rank files.
   * - ``&output_nml compress_output``
     - ``.false.``
     - Deflate compression for NetCDF.
   * - ``&output_nml compress_level``
     - ``1``
     - Deflate level, 1 to 9.
   * - ``&logging_nml status_interval``
     - ``0``
     - Console cadence in ``time_unit``. 0 means every 100 steps.
   * - ``&time_nml time_unit``
     - ``"s"``
     - Unit for ``dt_out``, ``status_interval`` and ``t_end``.

The diagnostic stream does not read ``&output_nml output_to_file``. It is
controlled only by ``&ocean_diag_nml enabled``. ``use_io_server`` is not
supported on the ocean path, and the model warns if it is set.

From Python, :doc:`python_interface` covers ``rdb.available_diagnostics()``
and the ``Diagnostics`` object, which writes this namelist group for you.


Adding a diagnostic
===================

To add a new field to the catalogue, follow the developer recipe in
``docs/howto/add_diagnostic.md``.
