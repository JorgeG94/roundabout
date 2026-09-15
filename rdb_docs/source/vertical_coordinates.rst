.. _vertical_coordinates:

--------------------
Vertical coordinates
--------------------

.. contents::
   :local:


The bottom-up convention
========================

.. warning::

   **Layers are indexed bottom-up.** :math:`k = 1` is the layer resting on
   the bed; :math:`k = N_z` is the surface layer.

This is the single most important thing to know about the vertical in
Roundabout, and it is the opposite of the convention several other ocean
models use. It follows ROMS. It is not negotiable and it is not
configurable.

Everything downstream assumes it:

* Surface forcing — wind stress, heat and salt flux, shortwave deposition,
  restoring — lands at ``(:, :, nz)``.
* Bed forcing — bottom drag, geothermal heat — lands at ``(:, :, 1)``.
* The hydrostatic pressure recursion of :ref:`theory` sums over layers
  *above* :math:`k`, which are the ones with the **larger** index.
* Vertical-coordinate target construction, per-layer drag, the PGF, the
  reference density profile, vertical advection and vertical diffusion all
  walk the column in this direction.

If you write a new layered kernel, follow it. The convention is guarded by
regression tests — ``test_ocean_conservation_salt_heat``,
``test_ocean_sw_penetration`` (positive :math:`Q` must warm :math:`k = N_z`
and preserve the stratification) and the ``test_ocean_vcoord*`` family will
all fail loudly on a flipped column, but they will fail in ways that take a
while to diagnose.


ALE: advance Lagrangian, then remap
===================================

Roundabout does not hold the layer interfaces fixed during a step. The
approach is **Arbitrary Lagrangian-Eulerian**:

#. **Advance.** The dynamics run with the layer interfaces free to move.
   Layer thickness is prognostic — :eq:`eq-continuity` transports it — so
   the layers simply go where the flow takes them. Within a step the
   coordinate is Lagrangian.
#. **Remap.** The column is then conservatively remapped from wherever the
   layers ended up onto the *target* arrangement that the chosen coordinate
   asks for.

The whole of the vertical-coordinate machinery therefore reduces to one
question: **what target thicknesses does this coordinate want, given the
column's total depth and its free surface?** Every family below answers that
question differently, and then goes through the identical remap.

The remap itself is a conservative reconstruct-and-reintegrate operation
selected by ``&vcoord_nml remap_method``, whose accepted values are ``pcm``,
``plm``, ``ppm`` (the default), ``ppm_h4`` and ``pqm``. It conserves
:math:`h\theta` for every registered tracer to round-off; the regression
gate is drift :math:`\le 10^{-15}` over ten remaps.

The remap is gated on the thermodynamic cadence
(``&ocean_vmix_nml dt_therm_ratio``) rather than running every dynamics
step.

.. note::

   ``sigma`` is the one coordinate that needs no remap at all — its target
   is by construction the arrangement the Lagrangian step already produced,
   up to the column rescaling. The code short-circuits accordingly.


The coordinate families
=======================

``&vcoord_nml vcoord_type`` selects the family. The accepted spellings are
``sigma`` (the default), ``zsigma``, ``zstar``, ``zstar_full``,
``zstar_sigma``, ``z``, ``eulerian_z``, ``isopycnal``, ``lagrangian``,
``gprime``, ``z_fixed``, ``rho`` and ``hycom`` — thirteen spellings over ten
distinct coordinate families, since ``z``/``eulerian_z``,
``isopycnal``/``lagrangian`` and ``gprime``/``z_fixed`` are pairs of aliases.

Terrain-following
-----------------

``sigma``
   Pure terrain-following. Every layer is a fixed fraction of the local
   water column: :math:`h_k = \mathrm{d}\sigma_k \, (H + \eta)`. Layers
   never vanish, however shallow the water, so nothing ever divides by a
   near-zero thickness. This is the default and the most robust choice.

   Its weakness is the **pressure-gradient error**: over steep bathymetry,
   surfaces of constant :math:`\sigma` cut steeply across surfaces of
   constant density, and the two large terms of the horizontal pressure
   gradient nearly cancel, leaving a spurious residual. Where that matters,
   move to a z-like family or improve the PGF form
   (:ref:`pressure_gradient`).

   **Use it for** shallow, gently sloped, or strongly tidal domains, for
   wetting-and-drying, and as the first thing to try when a run is
   misbehaving and you suspect the vertical coordinate.

Fixed-depth
-----------

``z`` / ``eulerian_z``
   The classic Eulerian z-coordinate: interfaces sit at fixed depths and do
   not move with the free surface. No pressure-gradient error by
   construction, because coordinate surfaces are geopotential surfaces.

   .. warning::

      ``eulerian_z`` is outside the ``pred_corr`` outer-scheme envelope.
      Since ``pred_corr`` is the default, a ``eulerian_z`` configuration
      must also set ``&ocean_bt_nml split_scheme = "ssp_rk2"``, and
      ``validate_config`` will refuse the run fail-loud otherwise, naming
      the fix.

``z_fixed`` / ``gprime``
   A prescribed stack of interface depths. This is the reduced-gravity
   setup: a small number of thick layers whose interfaces you specify,
   typically paired with ``&ocean_pgf_nml form = "gprime"``. The
   MOM6-reference double gyre is configured this way.

z-star
------

The z-star families stretch a reference level arrangement by the free
surface, so that the coordinate tracks the SSH instead of being cut by it.

``zstar``
   *z\*-lite.* A single global reference profile :math:`z_{\mathrm{ref}}`,
   stretched per column by the local :math:`(H + \eta)/H`. SSH-tracking,
   conservative, and robust: layers do not vanish. The pragmatic
   general-purpose choice when ``sigma``'s pressure-gradient error is a
   problem.

``zstar_full``
   *Full per-column z\*.* Each column builds its own
   :math:`z_{\mathrm{ref}}(0{:}N_z)` from its own bathymetry, and the
   surface layer is anchored at ``&vcoord_nml zstar_h_surf_target``
   regardless of the total depth — so near-surface resolution stays
   uniform across a domain with wildly varying depth, which is exactly what
   you want when the boundary layer is the thing you are resolving.

   The price is that **bed-side layers can vanish** where the water is
   shallower than the reference stack is deep. They are floored at
   ``&vcoord_nml zstar_h_min`` (default :math:`10^{-4}\ \mathrm{m}`), and
   every operator that divides by a layer thickness gates on the vanished-
   layer tolerance.

   .. warning::

      In intertidal domains with frequent wetting and drying,
      ``zstar_full`` still leaks on the order of 1–2 % of salt per
      wet/dry cycle. Use ``sigma`` or ``zstar`` there.

``zstar_sigma``
   A hybrid: sigma in shallow water, z\*-lite in deep. Conserves by
   construction. The intent is to get z\*'s clean pressure gradient in the
   deep interior without inheriting its vanishing-layer problem on the
   shelf.

``zsigma``
   A smoothstep blend from sigma near the surface to z at depth,
   conservatively remapped. The blend depth and width are internal
   constants on the coordinate object, not namelist knobs.

Isopycnal and hybrid
--------------------

``isopycnal`` / ``lagrangian``
   No remap at all — the layer *is* the coordinate, and it goes wherever
   the flow puts it. This is a true isopycnal model in the stacked
   shallow-water sense. Pair it with the ``&ocean_isopycnal_nml`` grounding
   controls.

``rho``
   Places layer interfaces on prescribed potential-density surfaces
   :math:`\rho_{\mathrm{target}}(0{:}N_z)` — lightest at the surface,
   densest at the bed — referenced to ``&vcoord_nml rho_ref_pressure``
   (default :math:`2\times10^{7}\ \mathrm{Pa}`). The column density profile
   is PPM-reconstructed from :math:`T` and :math:`S` through the EOS and
   then inverted, by a bracketing sweep followed by a fixed eight-iteration
   Newton solve, for the depth at which :math:`\rho` equals each interior
   target. Conserves :math:`h\theta` to round-off.

   .. warning::

      Weakly stratified columns collapse: with no density contrast there is
      no depth that matches the target, and the layers pile up. ``rho``
      alone is **validation-grade**; for production use ``hycom``.

``hycom``
   The production hybrid, and the answer to ``rho``'s collapse. It runs the
   same density-space inversion, then applies a z\* **surface-resolution
   floor**: a top-down sweep enforces
   :math:`z(k) \ge \sum \mathrm{d}\sigma \, (H + \eta)`, so the near-surface
   band keeps fixed z\* resolution while the deep interior tracks
   isopycnals. Reuses ``rho_target`` and ``rho_ref_pressure`` plus the
   existing :math:`\mathrm{d}\sigma` for the z\* band. Note the floor is a
   *surface-side minimum* only — the deep interior may still collapse, by
   design.


Choosing one
============

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - If you are
     - Start with
   * - Getting a new domain to run at all
     - ``sigma``. It cannot produce a vanishing layer, so it removes a
       whole class of failure from the picture while you debug everything
       else.
   * - Running shallow, tidal, or wetting-and-drying
     - ``sigma``, or ``zstar`` if you need the deep column too.
       ``zstar_full`` leaks across wet/dry cycles.
   * - Fighting spurious currents over steep topography
     - ``zstar``, then ``zstar_sigma``. Also revisit
       :ref:`pressure_gradient` — the coordinate and the PGF form are two
       handles on the same problem.
   * - Resolving a surface boundary layer over varying depth
     - ``zstar_full``, with ``zstar_h_surf_target`` set to the resolution
       you need.
   * - Running a reduced-gravity or idealised layered setup
     - ``z_fixed`` with ``&ocean_pgf_nml form = "gprime"``.
   * - Doing interior water-mass work where spurious diapycnal mixing is
       the enemy
     - ``hycom``. (``rho`` only if you can guarantee stratification
       everywhere.)

.. note::

   ``&vcoord_nml thickness_config`` selects the **initial** layer-thickness
   profile — accepted values ``sigma`` (default) and ``uniform_z`` — and is
   a separate thing from the running coordinate. It is easy to conflate the
   two; they answer different questions, one about the state at
   :math:`t = 0` and one about the target the remap aims at every
   thermodynamic step.


Vanishing layers
================

Two constants govern thin layers, and picking the wrong one is a recurring
source of bugs:

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Constant
     - Value
     - Role
   * - ``H_VANISHED``
     - :math:`1.5\times10^{-4}\ \mathrm{m}`
     - The *dynamic* vanish threshold. A layer this thin carries negligible
       mass and should be **skipped or merged**, not clamped. Kernels that
       would otherwise divide by it gate on this.
   * - ``H_DIV_EPS``
     - :math:`10^{-20}`
     - Pure divide-by-zero armour. Not a physical threshold — it exists
       only so that :math:`1/0` cannot happen.

Both live in ``rdb_constants`` with their roles documented. A new kernel
picks whichever matches its intent: if the question is "does this layer
matter physically?", it is ``H_VANISHED``; if the question is "could this
denominator be exactly zero?", it is ``H_DIV_EPS``.
