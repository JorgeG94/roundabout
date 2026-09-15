.. _theory:

-------------------
Governing equations
-------------------

.. contents::
   :local:

Roundabout integrates the **hydrostatic, Boussinesq primitive equations in
layer form**. This page sets out the continuous system the solver
discretises, and defines every symbol that appears elsewhere in this
documentation. How the system is discretised is the subject of
:ref:`discretisation`; the pressure-gradient term in particular gets its own
page, :ref:`pressure_gradient`.

Two properties of the formulation drive almost every design decision in the
code, so they are worth stating before the algebra:

#. **Layer thickness is a prognostic variable.** Continuity is a transport
   equation for :math:`h`, not a constraint on the velocity field. There is
   no Poisson solve and no FFT projection anywhere in the dynamical core.
#. **The vertical coordinate is generalised.** The equations below never
   commit to a particular set of layer interfaces. The interfaces are free
   to move with the flow within a step, and are restored to a target
   arrangement afterwards — the ALE approach described in
   :ref:`vertical_coordinates`.


The layer decomposition
=======================

The water column is divided into :math:`N_z` layers. Layers are indexed
**bottom-up**: :math:`k = 1` is the layer resting on the bed, and
:math:`k = N_z` is the surface layer. This convention is load-bearing
throughout the code and is discussed in :ref:`vertical_coordinates`.

Layer :math:`k` has thickness :math:`h_k > 0`, and the thicknesses sum to
the total water depth,

.. math::

   \sum_{k=1}^{N_z} h_k \;=\; H + \eta ,

where :math:`H(x,y)` is the resting depth of the bed below the reference
geoid :math:`z = 0` and :math:`\eta(x,y,t)` is the free-surface elevation
above it. The interface depths follow by summation from either end; the
free surface sits at :math:`z = \eta` and the bed at :math:`z = -H`.

Each layer carries a horizontal velocity :math:`\mathbf{u}_k = (u_k, v_k)`,
which is understood as the layer-average of the three-dimensional
horizontal velocity, and a set of tracer concentrations :math:`\theta_k`.


Layer continuity
================

Mass conservation for a layer, with the vertical coordinate free to move, is

.. math::
   :label: eq-continuity

   \frac{\partial h_k}{\partial t}
     \;+\; \nabla_h \!\cdot\! \left( h_k \, \mathbf{u}_k \right)
     \;=\; 0 ,

where :math:`\nabla_h = (\partial_x, \partial_y)` is the horizontal gradient
taken **along the layer**, not along a surface of constant :math:`z`.

Equation :eq:`eq-continuity` is the reason the solver has no elliptic stage.
In a fixed-coordinate model, the analogous statement is
:math:`\nabla \!\cdot\! \mathbf{u} = 0`, a *constraint* that must be
projected onto — typically by solving a Poisson problem for pressure. Here
the same physics appears as an *evolution equation* for a prognostic field.
Summing :eq:`eq-continuity` over :math:`k` recovers the free-surface
equation,

.. math::

   \frac{\partial \eta}{\partial t}
     \;+\; \nabla_h \!\cdot\! \sum_{k} h_k \mathbf{u}_k \;=\; 0 ,

which is the equation the barotropic sub-cycle advances; see
:ref:`discretisation`.

Vertical motion is not a prognostic variable. Where a diagnostic vertical
velocity is needed — for the vertical advection of momentum and tracers
across moving interfaces — it is recovered from the divergence of the layer
transports, as the residual between the layer's actual thickness tendency
and the thickness tendency the coordinate demands.


Momentum
========

Horizontal momentum is carried in **vector-invariant form**, which
re-expresses the nonlinear advection term using the identity

.. math::

   \left( \mathbf{u} \!\cdot\! \nabla_h \right) \mathbf{u}
     \;=\; \zeta \, \hat{\mathbf{z}} \times \mathbf{u}
         \;+\; \nabla_h \! \left( \tfrac{1}{2} |\mathbf{u}|^2 \right) ,

so that advection and rotation appear in the same term. The layer momentum
equation is

.. math::
   :label: eq-momentum

   \frac{\partial \mathbf{u}_k}{\partial t}
     \;+\; \left( f + \zeta_k \right) \hat{\mathbf{z}} \times \mathbf{u}_k
     \;+\; \nabla_h K_k
     \;=\;
     -\frac{1}{\rho_0} \nabla_h p_k
     \;+\; \nabla_h \!\cdot\! \left( \nu_h \nabla_h \mathbf{u}_k \right)
     \;+\; \frac{\partial}{\partial z}
           \left( \nu_v \frac{\partial \mathbf{u}_k}{\partial z} \right)
     \;+\; \mathbf{F}_k ,

with

.. math::

   \zeta_k \;=\; \frac{\partial v_k}{\partial x}
                  - \frac{\partial u_k}{\partial y} ,
   \qquad
   K_k \;=\; \tfrac{1}{2} \left( u_k^2 + v_k^2 \right) .

The symbols are:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Symbol
     - Meaning
   * - :math:`\mathbf{u}_k`
     - Layer horizontal velocity, :math:`\mathrm{m\,s^{-1}}`. Stored on
       C-grid faces: :math:`u` on x-faces, :math:`v` on y-faces.
   * - :math:`f`
     - Coriolis parameter, :math:`2\Omega \sin\phi`, :math:`\mathrm{s^{-1}}`.
       On a Cartesian grid this is the :math:`f`- or :math:`\beta`-plane
       value; on a spherical or curvilinear grid it is evaluated from the
       cell latitude.
   * - :math:`\zeta_k`
     - Relative vorticity, :math:`\mathrm{s^{-1}}`. Lives at cell corners.
   * - :math:`K_k`
     - Kinetic energy per unit mass, :math:`\mathrm{m^2\,s^{-2}}`. Lives at
       cell centres.
   * - :math:`p_k`
     - Pressure in the layer, :math:`\mathrm{Pa}`.
   * - :math:`\rho_0`
     - Boussinesq reference density, :math:`\mathrm{kg\,m^{-3}}`.
   * - :math:`\nu_h, \nu_v`
     - Horizontal and vertical eddy viscosity,
       :math:`\mathrm{m^2\,s^{-1}}`. Both are supplied by closures — see
       :ref:`parameterizations`.
   * - :math:`\mathbf{F}_k`
     - Everything else: surface wind stress at :math:`k = N_z`, bottom drag
       at :math:`k = 1`, tidal body forcing, eddy parameterization
       tendencies, sponge relaxation.

The quantity :math:`(f + \zeta_k)` is the absolute vorticity, and dividing
it by :math:`h_k` gives the **potential vorticity**
:math:`q_k = (f + \zeta_k)/h_k`. The Coriolis-plus-advection term is
discretised as a flux of :math:`q`, which is what makes the scheme
enstrophy- or energy-conserving; see :ref:`discretisation`.

.. note::

   The Boussinesq approximation appears in :eq:`eq-momentum` as the constant
   :math:`\rho_0` in the pressure-gradient term: density variations are
   retained where they produce buoyancy and neglected where they only
   change inertia. Roundabout is Boussinesq throughout — there is no
   non-Boussinesq branch.


Hydrostatic balance
===================

The vertical momentum equation is reduced to hydrostatic balance,

.. math::
   :label: eq-hydrostatic

   \frac{\partial p}{\partial z} \;=\; -\,\rho \, g ,

which is to say the vertical acceleration :math:`\mathrm{D}w/\mathrm{D}t`
and the vertical friction and Coriolis terms acting on :math:`w` are all
discarded. This is an excellent approximation whenever the horizontal scale
of the motion greatly exceeds the vertical — true of essentially everything
at the basin, mesoscale and submesoscale range that Roundabout targets, and
false for convective plumes and internal-wave breaking, which must instead
be parameterized (:ref:`parameterizations`).

Integrating :eq:`eq-hydrostatic` downward from the surface gives the
pressure at layer :math:`k`,

.. math::
   :label: eq-pressure

   p_k \;=\; p_{\mathrm{surf}}
             \;+\; g \!\! \sum_{m > k} \rho_m h_m
             \;+\; \tfrac{1}{2} \, g \, \rho_k h_k ,

the sum running over the layers *above* :math:`k` — which, in the bottom-up
convention, are the ones with the larger index. The half-layer term
evaluates the pressure at the mid-depth of layer :math:`k` itself.
:math:`p_{\mathrm{surf}}` is the load at the free surface: atmospheric
pressure when ``&ocean_psurf_nml`` is enabled, sea-ice mass loading
when that is wired, and zero otherwise.

Because Roundabout is a **layered** model, :math:`-\nabla_h p` cannot be
evaluated by simply differencing :eq:`eq-pressure` between neighbouring
columns: the layers in the two columns sit at different depths, so a naive
difference contains a spurious contribution from the tilt of the layer
itself. Handling this correctly is the whole content of the
pressure-gradient discretisation, and it is treated separately in
:ref:`pressure_gradient`.


Tracers
=======

Each registered tracer :math:`\theta` — salinity, temperature, and any
passive tracer added to the registry — obeys a **flux-form** conservation
law written for the thickness-weighted quantity :math:`h_k \theta_k`:

.. math::
   :label: eq-tracer

   \frac{\partial \left( h_k \theta_k \right)}{\partial t}
     \;+\; \nabla_h \!\cdot\! \left( h_k \mathbf{u}_k \theta_k \right)
     \;=\;
     \nabla_h \!\cdot\! \left( \kappa_h \, h_k \nabla_h \theta_k \right)
     \;+\; \frac{\partial}{\partial z}
           \left( \kappa_v \frac{\partial \theta_k}{\partial z} \right)
     \;+\; Q_\theta ,

with :math:`\kappa_h` and :math:`\kappa_v` the horizontal and vertical
tracer diffusivities and :math:`Q_\theta` the source term — surface heat and
salt flux at :math:`k = N_z`, penetrating shortwave radiation distributed
through the column, geothermal heat at :math:`k = 1`, restoring, and sponge
relaxation.

Writing :eq:`eq-tracer` in flux form for :math:`h\theta` rather than as an
advection equation for :math:`\theta` is what makes conservation exact
rather than approximate. The two conserved quantities the solver reports
each status line are precisely the column sums of :math:`h S` and
:math:`h T`; see :ref:`getting_started`.

**Consistency with continuity (CWC).** The tracer flux in
:eq:`eq-tracer` must use the *same* mass flux :math:`h_k \mathbf{u}_k` that
:eq:`eq-continuity` used, face by face and layer by layer. If it does not,
a spatially uniform tracer field stops being uniform — the model invents
structure out of nothing. Every advection path in the code is written to
satisfy this identity exactly, and the property is regression-tested.

Vertical tracer diffusion is solved as a backward-Euler tridiagonal system
per column, which is unconditionally stable and therefore imposes no time-step
restriction of its own, however thin the layers become.


Equation of state
=================

Density closes the system:

.. math::
   :label: eq-eos

   \rho_k \;=\; \rho\!\left( S_k,\, T_k,\, p \right) ,

with :math:`S` the practical salinity (PSU) and :math:`T` the potential
temperature (:math:`{}^\circ\mathrm{C}`). Four variants ship, selected by
``&ocean_eos_nml eos``, whose accepted values are ``linear`` (the
**default**), ``wright``, ``roquet_spv`` and ``teos10``:

**Wright (1997)** (``wright``) — the recommended production path, and what
every realistic configuration in ``validation_examples/ocean/`` selects. A
rational-function fit to the
equation of state of seawater, of the form

.. math::

   \rho(S,T,p) \;=\;
     \frac{p \;+\; p_0(S,T)}
          {\lambda(S,T) \;+\; \alpha_0(S,T)\,\bigl(p + p_0(S,T)\bigr)} ,

in which :math:`p_0`, :math:`\lambda` and :math:`\alpha_0` are low-order
polynomials in :math:`S` and :math:`T`. The rational form is what makes it
cheap to also obtain the compressibility and the thermal and haline
expansion coefficients analytically, which several closures need.

**Roquet et al. (2015)** (``roquet_spv``) — a specific-volume polynomial
fitted to TEOS-10, available for validation. It is mutually exclusive with
``&ocean_pgf_nml form="fv_wright"``, and the combination is refused at
configure time.

**TEOS-10** (``teos10``) — a further TEOS-10-flavoured variant.

**Linear** (``linear``, the default) — for analytical test cases and for
reduced-gravity configurations:

.. math::
   :label: eq-eos-linear

   \rho_k \;=\; \rho_0
                \;+\; \beta_S \left( S_k - S_{\mathrm{ref}} \right)
                \;-\; \alpha_T \left( T_k - T_{\mathrm{ref}} \right) .

Here :math:`\beta_S` is the haline contraction coefficient and
:math:`\alpha_T` the thermal expansion coefficient. Setting
:math:`\alpha_T = \beta_S = 0` decouples the tracers from the dynamics
entirely, giving an adiabatic reduced-gravity configuration — this is how
the MOM6-reference double gyre is run.

.. note::

   On the ocean path the linear-EOS coefficients live on the EOS object and
   default to :math:`\rho_0 = 1035\ \mathrm{kg\,m^{-3}}`,
   :math:`\alpha_T = 1.7\times10^{-4}\ \mathrm{kg\,m^{-3}\,{}^\circ C^{-1}}`
   and :math:`\beta_S = 7.6\times10^{-4}\ \mathrm{kg\,m^{-3}\,PSU^{-1}}`.
   Only :math:`\alpha_T` and :math:`\rho_0` are namelist-settable, as
   ``&ocean_ic_nml alpha_T`` and ``&ocean_ic_nml rho_0``;
   :math:`\beta_S`, :math:`S_{\mathrm{ref}}` and :math:`T_{\mathrm{ref}}`
   are compiled-in. The similarly-named ``&tracer_nml alpha_T`` and
   ``&tracer_nml beta_S`` keys are **coastal-era leftovers and are dead on
   the ocean path** — the configuration schema marks them as such, and
   setting them changes nothing.


What is deliberately absent
===========================

Stating the omissions is as useful as stating the equations:

* **No non-hydrostatic branch.** :eq:`eq-hydrostatic` is the vertical
  momentum equation, full stop. There is no :math:`\mathrm{D}w/\mathrm{D}t`
  and no pressure-correction Poisson solve on this path.
* **No semi-implicit free surface.** The fast barotropic mode is
  *sub-cycled* explicitly, not treated implicitly in the Casulli sense.
* **No online bulk flux formulae.** :math:`\mathbf{F}` and
  :math:`Q_\theta` must be supplied as stress and flux on the model grid;
  the model does not convert from :math:`(U_{10}, T_{\mathrm{air}},
  q_{\mathrm{air}}, \mathrm{SST}, \mathrm{SLP})`.
* **No biogeochemistry or sediment.** The tracer registry will carry
  passive tracers, but no source terms beyond the surface/bottom flux set
  are provided.
