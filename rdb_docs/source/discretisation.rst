.. _discretisation:

--------------
Discretisation
--------------

.. contents::
   :local:

This page describes how the equations of :ref:`theory` become an algorithm:
where the variables live, how continuity is solved, how the fast and slow
modes are separated, and how the two outer time schemes differ. The
pressure-gradient term is involved enough to get its own page,
:ref:`pressure_gradient`.


The Arakawa C-grid
==================

Variables are staggered on an Arakawa C-grid. Four positions exist in a
cell, and the code names them after the MOM6 convention:

.. list-table::
   :header-rows: 1
   :widths: 12 16 30 42

   * - Tag
     - Position
     - Array shape
     - What lives there
   * - ``T``
     - Cell centre
     - ``(nx, ny)``
     - :math:`\eta`, :math:`h`, :math:`\rho`, every tracer
       (:math:`S`, :math:`T`, …), kinetic energy :math:`K`, cell area
   * - ``Cu``
     - x-face
     - ``(nx+1, ny)``
     - :math:`u`, :math:`hu`, the zonal mass flux, ``dy_cu``
   * - ``Cv``
     - y-face
     - ``(nx, ny+1)``
     - :math:`v`, :math:`hv`, the meridional mass flux, ``dx_cv``
   * - ``Bu``
     - Cell corner
     - ``(nx+1, ny+1)``
     - Relative vorticity :math:`\zeta`, the Coriolis parameter
       :math:`f`, potential vorticity

Adding the vertical index gives the layered arrays: ``h_layer`` is
``(nx, ny, nz)``, ``u_face_x_layer`` is ``(nx+1, ny, nz)``,
``v_face_y_layer`` is ``(nx, ny+1, nz)``, and the diagnostic vertical
velocity ``w_interface`` sits at cell centres on layer *interfaces*,
``(nx, ny, nz+1)``.

Face index :math:`i` sits between cells :math:`i-1` and :math:`i`, so a
cell's divergence reads faces :math:`i` and :math:`i+1`.

The C-grid is chosen because it puts the divergence and the pressure
gradient on the same two-point stencil — there is no averaging between the
:math:`\eta` gradient and the velocity it accelerates, so the grid supports
gravity waves without the computational mode that an A-grid produces. The
price is that the Coriolis term now needs the *other* velocity component,
which is not co-located, so it must be averaged. That averaging is where
the Coriolis discretisation earns its complexity.

.. note::

   Arrays are stored **structure-of-arrays** — one allocatable per field,
   indexed ``(i, j, k)`` — so that a ``do concurrent`` loop over the
   fastest-varying index reads contiguous memory. This is a GPU
   requirement, not a stylistic choice; see :ref:`extending`.

All metric arrays are filled **including ghost rows and columns**. An
unfilled ghost is a recurring bug class: it lets the EOS fall back to
:math:`\rho = \rho_0` at a wall-adjacent face, which produces a spurious
density jump and a blow-up with an e-folding time of hours.


Continuity: PPM transport
=========================

Layer continuity :eq:`eq-continuity` is solved by a **Lin & Rood
directionally split continuity-PPM** scheme. Per layer, per step:

.. code-block:: text

   1. zonal_flux                 Phi_x from h^n
   2. tracer_advect_zonal        hTr <- hTr - dt d(Phi_x T)/dx   at h^n
   3. apply_zonal                h   <- h^n - dt dPhi_x/dx       (= h*)
   4. meridional_flux            Phi_y from h*
   5. tracer_advect_meridional   hTr <- hTr - dt d(Phi_y T)/dy   at h*
   6. apply_meridional           h   <- h* - dt dPhi_y/dy        (= h^{n+1})

The mass flux through an x-face is

.. math::

   \Phi_x(i,j,k) \;=\; u(i,j,k) \; h_{\mathrm{face}}(i,j,k) \;
                       \Delta y_{cu}(i,j) ,

and the thickness update that consumes it is

.. math::

   h(i,j,k) \;\leftarrow\; h(i,j,k)
     \;-\; \Delta t \,
       \bigl[ \Phi_x(i{+}1,j,k) - \Phi_x(i,j,k) \bigr] \,
       \frac{1}{A_T(i,j)} .

The face thickness :math:`h_{\mathrm{face}}` is the interesting part. It is
a **swept-average**: a piecewise-parabolic reconstruction of :math:`h`
(Colella & Woodward 1984, eq. 1.6, with the eq. 1.10 monotonicity limiter)
is integrated over the region the flow sweeps across the face during
:math:`\Delta t`, Godunov-style, and the upwind side supplies the
reconstruction. Integrating the reconstruction over the swept region rather
than evaluating it at a point is what makes the scheme single-stage stable;
the classical point-value form is weakly unstable in this integrator.

Notice what is *absent*. There is no velocity projection, no pressure
Poisson solve, no FFT. Thickness is advanced by a flux divergence exactly
as a tracer would be. That single property is what lets the vertical
coordinate be Lagrangian within a step, which in turn is what makes the ALE
remap of :ref:`vertical_coordinates` a drop-in stage rather than a rewrite.

**Consistency with continuity.** Steps 2 and 5 above use the *same*
:math:`\Phi` that steps 3 and 6 use. This is not an optimisation, it is a
correctness requirement: the discrete tracer equation must reduce to the
discrete continuity equation when :math:`\theta \equiv 1`, or a uniform
tracer field stops being uniform. The layer fluxes are additionally
renormalised so that :math:`\sum_k \Phi_{x,k}` equals the barotropic
transport the fast loop produced, by a uniform velocity correction
:math:`\delta u = (\overline{uh} - \sum_k \Phi_k) / \sum_k
h_{\mathrm{face}}`. That keeps the slow thickness advance consistent with
the fast loop's end-of-step :math:`\eta` with no post-hoc rescaling.


Coriolis and advection
======================

The vector-invariant combination :math:`(f + \zeta)\hat{\mathbf{z}} \times
\mathbf{u}` is discretised as a **flux of potential vorticity** —
Sadourny's scheme. The corner vorticity

.. math::

   \zeta(i,j) \;=\; \Bigl[
     \bigl( v(i,j)\,\Delta y_{cv}(i,j) - v(i{-}1,j)\,\Delta y_{cv}(i{-}1,j) \bigr)
     - \bigl( u(i,j)\,\Delta x_{cu}(i,j) - u(i,j{-}1)\,\Delta x_{cu}(i,j{-}1) \bigr)
   \Bigr] \frac{1}{A_{Bu}(i,j)}

is combined with :math:`f` at the same corner and interpolated onto the
faces where the momentum tendency is needed. ``&ocean_coriolis_nml form=``
selects the variant:

``sadourny`` (default)
   The enstrophy-conserving velocity form.

``sadourny_energy``
   The energy-conserving transport form.

``sadourny_hk``
   Adds the Hollingsworth-Källén correction, which suppresses a
   computational instability the plain enstrophy form admits at high
   resolution. Carried as a guard option rather than the production path.

An orthogonal knob, ``&ocean_coriolis_nml pv_adv_scheme=``, chooses how the
corner vorticity is interpolated onto the face. The default ``centered`` is
a two-point average and is bit-identical to the historical behaviour.
``weno3``, ``weno5`` and ``weno7`` instead use an upwind-biased WENO-Z
reconstruction, which sharpens submesoscale PV fronts and supplies
scale-selective dissipation that controls grid-scale vorticity noise where
explicit viscosity is low. The wider stencils need more ghost cells —
``weno5`` requires ``nghost >= 3`` and ``weno7`` requires ``nghost >= 4``,
checked fail-loud at configure.


The split-explicit mode split
=============================

The system carries waves whose speeds differ by two orders of magnitude.
External (barotropic) gravity waves travel at :math:`\sqrt{gH} \approx
200\ \mathrm{m\,s^{-1}}` in the deep ocean; the internal waves, eddies and
currents that carry the physics move at :math:`\mathcal{O}(1)\
\mathrm{m\,s^{-1}}`. Resolving the fast wave with a single explicit step
would force :math:`\Delta t \sim \Delta x / \sqrt{gH}` — seconds — on
everything.

The mode split fixes this by treating the two separately:

**The slow (baroclinic) part** — Coriolis-advection, the pressure gradient,
lateral viscosity, vertical mixing, bottom drag, surface stress, the tracer
equations — is evaluated once per outer step at :math:`\Delta t`.

**The fast (barotropic) part** — the depth-integrated
:math:`(\eta, \overline{u}, \overline{v})` system, which is where the
gravity wave lives — is sub-cycled ``n_inner`` times at
:math:`\Delta t / n_{\mathrm{inner}}`, driven by a forcing
:math:`\mathbf{F}_{bt}` that is the depth-mean of the slow tendencies, held
fixed across the substeps.

**The correction** then folds the barotropic result back into the layers.

Consequently **the outer time step is limited by advective and baroclinic
CFL, not by** :math:`\Delta x/\sqrt{gH}`. In practice that is minutes:
:math:`\Delta t = 1200\ \mathrm{s}` in the MOM6-reference double gyre.

.. note::

   This is a **sub-cycled explicit** fast mode, not a semi-implicit free
   surface. Roundabout does not treat the gravity-wave term implicitly in
   the Casulli/SCHISM sense; the substeps resolve it explicitly, just at
   their own smaller step.

The barotropic substep
----------------------

Each substep is **forward-backward Euler**: update :math:`\eta` first,
consuming :math:`u^n` and :math:`v^n`, then update :math:`u` and :math:`v`
consuming the *just-updated* :math:`\eta`. Schematically, per substep, with
:math:`\delta t = \Delta t/n_{\mathrm{inner}}`:

.. math::

   \eta^{m+1} &= \eta^{m}
     - \delta t \, \nabla \!\cdot\! \bigl( h_{\mathrm{face}}\,
       \overline{\mathbf{u}}^{\,m} \bigr) \\[4pt]
   \overline{u}^{\,m+1} &= r_u \Bigl[ \overline{u}^{\,m} + \delta t \bigl(
       (\zeta + f) \overline{v}
       - g \, \partial_x \eta^{m+1}
       - \partial_x K
       + F_{bt,u} \bigr) \Bigr]

and the mirror expression for :math:`\overline{v}`. The factor :math:`r_u`
is a per-face remainder that is unity unless substep drag, barotropic wave
drag or land masking is active.

The :math:`\zeta` and :math:`K` terms carry the nonlinearity of the
barotropic mode, and ``&ocean_bt_nml substep_zeta_ke`` (default ``.true.``)
controls whether they are recomputed each substep or frozen into the slow
forcing. Freezing them is the MOM6-parity behaviour.

The substep loop reports both the **time-mean** :math:`\eta`, velocity and
transport over the substeps, and the **end-of-substep** values. Both are
used, for different things: layer continuity consumes the time-mean
transports, and the momentum correction consumes the end-of-step velocity.

Choosing ``n_inner``
--------------------

``&ocean_bt_nml n_inner`` sets the substep count directly. Setting
``&ocean_bt_nml auto_n_inner = .true.`` derives it instead from the
gravity-wave CFL:

.. math::

   \ell_{\mathrm{CFL}} \;=\;
     \min_{\mathrm{cells}}
     \frac{1}{\sqrt{\dfrac{1}{\Delta x^2} + \dfrac{1}{\Delta y^2}}} ,
   \qquad
   c_{\mathrm{ext}} \;=\; \sqrt{g \, H_{\max}} ,

.. math::

   \delta t_{\mathrm{safe}} \;=\;
     \frac{\texttt{cfl\_bt\_safety} \cdot \ell_{\mathrm{CFL}}}
          {c_{\mathrm{ext}}} ,
   \qquad
   n_{\mathrm{inner}} \;=\; \max\!\left(1,
     \left\lceil \frac{\Delta t}{\delta t_{\mathrm{safe}}} \right\rceil
   \right) .

``cfl_bt_safety`` defaults to ``0.65``.

The :math:`\ell_{\mathrm{CFL}}` expression is the **two-dimensional** CFL
length, and the cross-direction term matters: on a uniform Cartesian grid
it evaluates to :math:`\Delta x/\sqrt{2}`, not :math:`\Delta x`. A
one-dimensional estimate under-counts :math:`n_{\mathrm{inner}}` by a
factor of :math:`\sqrt{2}` and leaves the effective two-dimensional CFL at
:math:`\approx 0.92` — right on the edge of the forward-backward scheme's
stability, which is exactly how a 2 km configuration was once made to NaN.

.. warning::

   ``n_inner`` is derived **once, at configure time**, and the resolved
   value is written back into the configuration. It is *not* re-derived
   per step, and there is no adaptive controller. The truncation counter
   ``ntrunc_total`` that the console reports counts CFL clips; it does not
   feed back into the time step. If your bathymetry or your forcing
   changes the effective wave speed materially mid-run, you must choose
   :math:`\Delta t` and ``n_inner`` to cover the worst case.

The correction back into the layers
-----------------------------------

After the substeps, the barotropic increment is distributed to every layer:

.. math::

   \Delta \overline{u} \;=\;
     \overline{u}_{\mathrm{end}} \;-\; \overline{u}^{\,n}
     \;-\; \Delta t \, F_{bt,u} ,
   \qquad
   u_k \;\leftarrow\; u_k + \Delta \overline{u} .

Using the **end-of-step** barotropic velocity here, while layer continuity
used the **time-mean** transports, is the Hallberg (2009) convention and is
deliberate. Optionally the increment can be distributed with thickness
weighting (``correction_h_weighted``) or further biased by the viscous
remnant of the vertical-friction solve (``correction_visc_rem``), both
matching MOM6 options.

A non-finite increment is skipped and counted rather than applied — see the
note on NaN-laundering clamps in :ref:`extending`.


The outer time scheme
=====================

``&ocean_bt_nml split_scheme`` selects the outer time integration. Two
values are accepted: ``pred_corr``, the default, and ``ssp_rk2``. The two
schemes assemble the stages of an outer step differently, and they differ
in how they amplify or damp certain numerical artifacts, so results are not
bit-identical between them.

.. note::

   The spellings ``split_rk2`` and ``mom6_pc`` are **retired**. A namelist
   that still uses them fails at configure with a hint naming the
   replacement, rather than silently selecting a default.

``pred_corr``
-------------

A predictor-corrector. One outer step is a predictor followed by a
corrector; both run the same stage routine, and what differs is a set of
role gates.

**Predictor.** Advance with an off-centred fraction of the step,
:math:`\Delta t_{\mathrm{pred}} = \texttt{pc\_be} \cdot \Delta t`:

.. math::

   u^{p} \;=\; u^{n}
     \;+\; \mathrm{BE}\,\Delta t \,
       \bigl( \mathrm{CorAd} + \mathrm{PGF}
              + \mathrm{diffu}^{\,n-1} + \mathrm{drag}
              + \mathrm{stress} \bigr)
     \;+\; \mathrm{BE}\,\Delta \overline{u} ,

with :math:`\mathrm{BE} = \texttt{pc\_be}`, default **0.6**. The predictor
reuses the previous step's lateral-viscosity tendency rather than
recomputing it, and does no tracer physics.

**Corrector.** The provisional :math:`u^p`, :math:`v^p`, :math:`h^p` are
discarded; the state is restored to :math:`u^n` and the same expression is
evaluated at the *full* :math:`\Delta t`, with the slow tendencies taken on
the step time-means :math:`\overline{u}_{\mathrm{av}}`,
:math:`\overline{v}_{\mathrm{av}}`, :math:`\overline{h}_{\mathrm{av}}` that
the predictor produced. The corrector's single full-step update *is* the
step — there is no average of two stage outputs.

.. note::

   ``pc_be`` is an off-centring fraction of the step, not a weight in an
   average of stage outputs.

The other structural feature is **where continuity runs**: under
``pred_corr`` the continuity-and-tracer chain is deferred until *after* the
velocity update and the implicit friction, so that the thickness advances
with the updated velocities. That ordering is the forward-backward
gravity-wave pairing.

**The v1 envelope.** ``pred_corr`` is refused **fail-loud** at configure in
three cases, each of which names the fix in the error message:

* ``&vcoord_nml vcoord_type = "eulerian_z"``
* ``&ocean_wetdry_nml enable = .true.``
* ``&ocean_vmix_nml dt_tracer_advect_ratio > 1``

In each case the remedy is to set
``&ocean_bt_nml split_scheme = "ssp_rk2"`` explicitly.

``ssp_rk2``
-----------

A two-stage Heun-type scheme: two identical stages, then an arithmetic
average.

.. math::

   \phi^{(1)} &= \phi^{n} + \Delta t \, L(\phi^{n}) \\
   \phi^{(2)} &= \phi^{(1)} + \Delta t \, L(\phi^{(1)}) \\
   \phi^{n+1} &= \tfrac{1}{2} \left( \phi^{n} + \phi^{(2)} \right)

applied to :math:`h`, :math:`u`, :math:`v` and every tracer. Here the
continuity-and-tracer chain runs *before* the velocity update.

Its configuration envelope is the wider of the two: every vertical
coordinate including ``eulerian_z``, wetting and drying, and windowed
tracer advection. The six shipped namelists that use one of those features
pin it explicitly.
