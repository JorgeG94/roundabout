.. _pressure_gradient:

---------------------
The pressure gradient
---------------------

.. contents::
   :local:


Why this term is hard
=====================

In a layered model the horizontal pressure-gradient force is the one term
that cannot be discretised naively, and it is the classic source of
spurious currents in coordinate ocean models.

The difficulty is geometric. :eq:`eq-pressure` gives the pressure at layer
:math:`k` in a column, and the momentum equation wants
:math:`-\rho_0^{-1}\nabla_h p` **at constant depth**. But layer :math:`k` in
the column to the left of a face and layer :math:`k` in the column to the
right sit at *different depths* whenever the layer is tilted — over
bathymetry, under a sloping free surface, or simply because a z\*-layer
follows the SSH. Differencing :math:`p_k` between the two columns therefore
mixes the horizontal gradient you want with a vertical gradient you do not.

The vertical gradient is large. Hydrostatic pressure in the deep ocean is
:math:`\mathcal{O}(2\times10^{7}\ \mathrm{Pa})`, while the horizontal
pressure difference that actually drives the flow is
:math:`\mathcal{O}(10^{3}\ \mathrm{Pa})`. Recovering the small number as
the difference of two large ones carries a condition number of order
:math:`10^{4}`. Over steep topography in a terrain-following coordinate,
the residual of that near-cancellation can exceed the real signal — the
**pressure-gradient error**, which shows up as a steady spurious current
over a seamount that ought to be motionless.

Roundabout ships five discretisations of this term, and the choice matters.
Every one of them is measured **from the free surface downward**, so that
none of them carries the barotropic :math:`-g\nabla\eta` term: that is the
barotropic substep's job, and including it here would double-count it.


The five forms
==============

``&ocean_pgf_nml form=`` selects the variant. The accepted values are
``mont`` (the **default**), ``fv_lite``, ``fv_wright``, ``gprime`` and
``fv_mom6``.

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - ``form``
     - Use it when
   * - ``mont``
     - **Default.** General purpose, valid over any bathymetry and any
       vertical coordinate. Start here.
   * - ``fv_lite``
     - You want the straightforward finite-volume form, or you are
       reproducing an older result. Algebraically equivalent to ``mont``
       on aligned columns, but carrying the cancellation error.
   * - ``fv_wright``
     - Compressibility matters — deep water, large pressure range — and
       you are running the Wright EOS.
   * - ``gprime``
     - Two-layer reduced-gravity idealised setups. **Requires exactly
       two layers** and is refused fail-loud otherwise.
   * - ``fv_mom6``
     - You are reproducing MOM6 closely, or you need the barotropic
       ``pbce`` coupling, which only this form populates.


``mont`` — the Montgomery potential
===================================

The default, and the one worth understanding properly.

The idea
--------

Define the **Boussinesq Montgomery potential**

.. math::
   :label: eq-montgomery

   M \;=\; \frac{p}{\rho_0} \;+\; \rho_{*} \, z ,
   \qquad
   \rho_{*} \;=\; \frac{g \, \rho_{\mathrm{layer}}}{\rho_0} .

The point of this combination is that **inside a layer of horizontally
uniform density,** :math:`M` **is constant with depth.** Moving up by
:math:`\mathrm{d}z` within the layer costs
:math:`-g\rho\,\mathrm{d}z/\rho_0` of :math:`p/\rho_0` — that is exactly
hydrostatic balance :eq:`eq-hydrostatic` — and gains precisely
:math:`\rho_{*}\,\mathrm{d}z` from the second term. The two cancel
identically.

That constancy is what licenses a single horizontal difference. Because
:math:`M` does not vary through the layer, differencing it between two
columns cannot pick up a vertical contamination, and

.. math::
   :label: eq-mont-identity

   -\frac{1}{\rho_0} \left. \frac{\partial p}{\partial x} \right|_{z}
     \;=\; -\,\left. \frac{\partial M}{\partial x} \right|_{\mathrm{layer}} .

The left-hand side is evaluated at constant depth; the right-hand side is
evaluated *along the layer*, which is where the model's data actually is.
**One horizontal difference of** :math:`M` **is legitimate because**
:math:`M` **already carries the geopotential** :math:`\rho_{*} z`. Note also
that no :math:`\rho_0^{-1}` appears on the right-hand side: :math:`M`
already has units of geopotential, :math:`\mathrm{m^2\,s^{-2}}`, so its
horizontal gradient *is* an acceleration. The :math:`1/\rho_0` lives inside
:math:`\rho_{*}`.

Contrast this with the hydrostatic-stack forms below, which build :math:`p`
and then have to *correct* the difference with an explicit
:math:`g\rho\,\Delta z` term. The Montgomery form never builds a pressure
stack at all — the hydrostatic relation is folded analytically into the
definition of :math:`M`.

The vertical recursion
----------------------

:math:`M` is built by a recursion that runs **from the surface downward**.
In the bottom-up convention (:ref:`vertical_coordinates`) that means
starting at :math:`k = N_z` and marching to :math:`k = 1`:

.. math::
   :label: eq-mont-recursion

   M_{N_z} &= 0 , \\[4pt]
   M_{k}   &= M_{k+1}
              \;+\; \bigl( \rho_{*,k} - \rho_{*,k+1} \bigr) \, e_{k+1} ,
              \qquad k = N_z - 1, \ldots, 1 ,

where :math:`e_{k+1}` is the height of the interface **shared** by layers
:math:`k` and :math:`k+1` — the *top* of layer :math:`k`, recovered from
the layer centre as :math:`e_{k+1} = z_k + h_k/2`.

Two things to note.

**The seed.** :math:`M_{N_z} = 0`, not :math:`p_{\mathrm{surf}}/\rho_0`.
Heights are measured relative to the free surface — :math:`z = 0` at the
surface, negative below — so the surface-relative interface height and the
surface pressure both vanish there. This is what keeps the barotropic
:math:`-g\nabla\eta` out of the PGF and in the barotropic substep where it
belongs.

**The recursion step.** At the shared interface, :math:`p` and :math:`z`
agree between the two layers by construction, so the entire jump in
:math:`M` across the interface is the jump in :math:`\rho_{*}`. That is why
the step is the density *difference* times the interface height, and
nothing else.

The horizontal-density term
---------------------------

The acceleration on an x-face is then

.. math::
   :label: eq-mont-face

   \mathrm{PGF}_x(i,j,k) \;=\;
     \Bigl[
       -\bigl( M_{i,j,k} - M_{i-1,j,k} \bigr)
       \;+\; \Delta\rho_{*} \, z_{\mathrm{eff}}
     \Bigr] \frac{1}{\Delta x_{cu}(i,j)} ,

with

.. math::
   :label: eq-mont-zeff

   \Delta\rho_{*} \;=\;
     \frac{g}{\rho_0}\bigl( \rho_{i,j,k} - \rho_{i-1,j,k} \bigr) ,
   \qquad
   z_{\mathrm{eff}} \;=\;
     \frac{e_L h_R + e_R h_L - h_L h_R}{h_L + h_R} ,

where :math:`h_L, h_R` are the layer thicknesses either side of the face
and :math:`e_L, e_R` the corresponding top-interface heights. The y-face is
the exact mirror, using :math:`\Delta y_{cv}`.

**The second term is load-bearing, not a refinement.** This is worth
stating plainly, because the shape of :eq:`eq-mont-face` invites the
opposite reading — one big term plus one small correction.

The constancy argument behind :eq:`eq-mont-identity` assumed the layer
density was horizontally uniform *within the layer*. Where it is not — which
is to say, wherever there is any baroclinicity at all, which is the entire
point of running a layered ocean model — the exact relation picks up an
additional term :math:`+\,z\,\partial_x \rho_{*}`, and
:math:`z_{\mathrm{eff}}` is the thickness-weighted height at which to
evaluate it. Drop it and the scheme gets **the sign of a horizontal density
contrast wrong**. Not the magnitude — the sign.

The history here is instructive. Before this form was corrected, ``mont``
computed nothing but a horizontal difference of layer-centre pressure, with
no geopotential term at all:

.. math::

   \mathrm{PGF}_x \;=\;
     -\frac{1}{\rho_0} \,
      \frac{p^{\,\mathrm{centre}}_{R} - p^{\,\mathrm{centre}}_{L}}
           {\Delta x} .

That expression is exact only where :math:`\Delta z = 0` across the face —
a flat, aligned column — and it produced NaN on every sloping case tried.
It was named for a Montgomery potential it did not contain.

On **aligned columns**, where :math:`h_L = h_R`, :math:`z_{\mathrm{eff}}`
collapses to the arithmetic mean layer-centre height and
:eq:`eq-mont-face` reduces *algebraically* to ``fv_lite``. The two forms
are not rivals; ``mont`` is the generalisation that stays correct when the
columns are not aligned.

What it buys
------------

``mont`` is exact at rest in an isopycnal column over any bathymetry where
the layer is present on both sides of the face: flat isopycnals make every
interface height uniform and every :math:`\Delta\rho_{*}` zero, so both
terms vanish identically rather than nearly cancelling.

That property shows up in the measurements. On three quiescent
sigma-coordinate seamount cases — configurations whose correct answer is
*no motion* — the mean kinetic energy dropped by three to five orders of
magnitude when ``mont`` became the default, for example
:math:`\mathrm{En}\; 6.00\times10^{-4} \rightarrow 7.82\times10^{-7}\
\mathrm{m^2\,s^{-2}}`. On an adversarial at-rest stratified seamount it
agrees with ``fv_mom6`` to six digits. Eleven of the twenty-eight unpinned
namelists were bit-identical under the change, and stratified z\*-sigma
seamounts moved by less than :math:`5\times10^{-5}` relative.

The application, finally, is a forward-Euler accumulation onto the face
velocities:

.. math::

   u_k \;\leftarrow\; u_k + \Delta t \cdot \mathrm{PGF}_x ,
   \qquad
   v_k \;\leftarrow\; v_k + \Delta t \cdot \mathrm{PGF}_y .


``fv_lite`` — the finite-volume form
====================================

Builds a layer-edge pressure stack per column, from the surface down, using
the layer-mean density:

.. math::

   p_{N_z+1} = 0 ,
   \qquad
   p_{k} \;=\; p_{k+1} + g \, \rho_k \, h_k ,
   \qquad
   p^{\,\mathrm{centre}}_{k} \;=\; \tfrac{1}{2}\left( p_k + p_{k+1} \right) ,

then differences the cell-centred pressure across the face *and corrects it
back to constant* :math:`z`:

.. math::

   \mathrm{PGF}_x \;=\;
     -\frac{1}{\rho_0}
      \left[
        \frac{p^{\,\mathrm{centre}}_{R} - p^{\,\mathrm{centre}}_{L}}
             {\Delta x}
        \;+\; g \, \overline{\rho}_{\mathrm{face}} \,
              \frac{z_R - z_L}{\Delta x}
      \right] ,

with :math:`\overline{\rho}_{\mathrm{face}}` the two-point density average.
The second term is the explicit :math:`z`-correction that ``mont`` gets for
free from the structure of :math:`M`.

This form is correct, and it reduces to ``mont`` on aligned columns. Its
weakness is the cancellation error described at the top of this page: it
recovers an :math:`\mathcal{O}(10^{3})\ \mathrm{Pa}` signal by differencing
an :math:`\mathcal{O}(2\times10^{7})\ \mathrm{Pa}` stack.


``fv_wright`` — in-situ density
===============================

Identical face algebra to ``fv_lite``, but the pressure stack *and* the
face density are built from **in-situ** density, re-evaluated with the
Wright (1997) EOS at the pressure obtained from a single Picard iteration:

.. math::

   p^{\,\mathrm{centre}}_{k} \;=\;
     p^{\,\mathrm{above}} + \tfrac{1}{2} g \,
       \rho\!\left(T_k, S_k, p^{\,\mathrm{centre}}_{k}\right) h_k ,

after which the column pressure advances by
:math:`g \rho_k^{\,\mathrm{in\text{-}situ}} h_k`.

This adds the compressibility contribution that ``fv_lite`` misses, and
reduces to ``fv_lite`` for incompressible water. It needs the Wright EOS
for a sensible Picard seed; combining ``fv_wright`` with
``&ocean_eos_nml eos = "roquet_spv"`` is refused fail-loud at configure.


``gprime`` — two-layer reduced gravity
======================================

No EOS, no pressure stack, no iteration. The two densities are fixed by two
reduced gravities: ``&ocean_pgf_nml gprime_gfs`` (:math:`g_{FS}`, default
:math:`9.81\ \mathrm{m\,s^{-2}}`) and ``gprime_gint`` (:math:`g'`, default
:math:`0.0098\ \mathrm{m\,s^{-2}}`).

The free surface is recovered from the thickness sum and the bathymetry,
:math:`\eta = (h_1 + h_2) - b`, and the two accelerations are

.. math::

   a_{\mathrm{top}} \;(k = 2) &= -\, g_{FS} \, \nabla \eta \\[4pt]
   a_{\mathrm{bot}} \;(k = 1) &= -\, g_{FS} \, \nabla \eta
                                 \;-\; g' \, \nabla h_1

— the top layer feels the free-surface slope alone, the bottom layer feels
that plus the interface slope weighted by the reduced gravity.

.. warning::

   ``gprime`` requires exactly **two layers**, checked fail-loud at
   configure. Selecting it with a deeper stack would silently leave layers
   :math:`k > 2` with no pressure gradient at all.

When this form is selected, the barotropic fast loop also switches to
:math:`g_{FS}` for its gravity-wave speed.


``fv_mom6`` — the MOM6 port
===========================

A faithful port of MOM6's ``PressureForce_FV_Bouss`` on the Boussinesq
per-layer-density path. Instead of layer-centre pressure plus a
:math:`z`-correction, it works with **layer-integrated** pressure
differences divided by the face-averaged thickness. It builds interface
heights from the bed up, a pressure *anomaly* stack relative to
:math:`\rho_{\mathrm{ref}}\, g \, z`,

.. math::

   p^{a}_{N_z+1} = \rho_{\mathrm{ref}} \, g \, \eta ,
   \qquad
   p^{a}_{k} = p^{a}_{k+1}
     + \left( \rho_k - \rho_{\mathrm{ref}} \right) g \, h_k ,

together with the first moment
:math:`\tfrac{1}{2} (\rho_k - \rho_{\mathrm{ref}}) g h_k^2`, and assembles
the face force from those integrals. Working with the anomaly rather than
the full pressure is what controls the cancellation error.

Options specific to this form:

``mass_weight``
   Default ``.false.``. Enables the MOM6 shelf-break mass-weighting blend
   in the face integrals.

``gfs_scale``
   Default ``1.0``, at which the corresponding correction is a no-op. Below
   one it subtracts a Montgomery-style
   :math:`(1 - \texttt{gfs\_scale})\,(g/\rho_0)\,\rho_{\mathrm{surf}}
   \nabla\eta` from every layer, matching MOM6's ``GFS_SCALE``.

``reconstruct_for_pressure``
   Default ``.false.``, and **only valid with this form** — any other
   ``form`` is refused fail-loud. It replaces the piecewise-constant
   layer-mean density in the pressure stack with a five-point Boole
   quadrature of the EOS along a monotone PLM or PPM sub-layer
   :math:`T`/:math:`S` profile:

   .. math::

      \overline{\rho}^{\,a} \;=\; \frac{1}{90}
        \Bigl[ 7 (r_1 + r_5) + 32 (r_2 + r_4) + 12 \, r_3 \Bigr] ,

   with :math:`r_n` the density anomaly at the five equally spaced
   sub-points. The companion ``recon_scheme`` chooses the reconstruction —
   ``1`` for PLM (the default) or ``2`` for PPM. The face assembly is
   untouched; only the in-layer integration changes. It removes the
   spurious-PGF moment error on thick, sloped layers.

.. note::

   ``fv_mom6`` is the only form that populates the interface-height field
   the barotropic ``pbce`` coupling requires. Selecting that coupling with
   any other form is an ``error stop``.
