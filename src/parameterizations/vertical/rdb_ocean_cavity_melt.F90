!! Ice-shelf basal-melt interface thermodynamics — the three-equation
!! kernel for an ice-shelf cavity.  **KERNEL ONLY: nothing in this module
!! is called by the ocean step yet** (no namelist group, no state slot, no
!! engine wiring — that is the coupling PR).  What ships here is the
!! scalar, `pure`, `!$acc routine seq` physics, under unit test against a
!! 17-significant-digit golden oracle.
!!
!! SCOPE.  One water column at a time: given the far-field (T, S), the
!! interface pressure and a friction velocity, return the interface state
!! (T_b, S_b) and the melt MASS flux, plus the heat/salt fluxes that a
!! future coupling seam will hand to the tracer budgets.
!!
!! UNIT AND SIGN CONVENTIONS (the one table that matters — get a sign or a
!! density wrong here and everything downstream is plausible and wrong):
!!
!! | symbol | meaning | unit | sign |
!! |---|---|---|---|
!! | `m_mass` | **canonical** ICE mass leaving the ice per unit area | kg/m^2/s | **>0 MELTING**, <0 freezing |
!! | `T_b`, `S_b` | interface temperature / salinity (on the liquidus) | degC, g/kg | — |
!! | `T_w`, `S_w` | far-field ("mixed layer") temperature / salinity | degC, g/kg | — |
!! | `p_b` | interface pressure `g*rho_i*h_ice` (the undiluted ice load) | **Pa** | >= 0 |
!! | `gamma_t`, `gamma_s` | exchange VELOCITIES | m/s | > 0 |
!! | `u_star` | interface friction velocity | m/s | > 0 |
!! | `T_star` | thermal driving `T_w - T_f(S_w, p_b)` | degC | > 0 ⇒ melting |
!! | `q_ocean` | turbulent heat flux ocean → interface | W/m^2 | > 0 warms the interface (the ocean cools) |
!! | `q_ice` | conductive + ice-warming flux interface → ice | W/m^2 | > 0 into the ice |
!! | `q_latent` | latent heat consumed by the phase change | W/m^2 | > 0 when melting |
!! | `b_flux` | interfacial buoyancy flux | m^2/s^3 | **< 0 STABILISING** (melting) |
!! | `l_plus` | viscous Obukhov scale `L/delta_nu` | — | **> 0 stabilising** |
!!
!! `m_mass` is canonical because the literature uses FOUR incompatible
!! densities for the same physical melt rate (ice 916/918/920, freshwater
!! 1000, seawater 1025/1027).  A 30 m/yr ISOMIP+ melt rate is 32.7 m/yr of
!! solid ice and 29.2 m/yr of seawater-equivalent.  The mass flux is the
!! only convention-free number, so it is what this kernel returns;
!! `cavity_m_ice_from_mass` / `cavity_m_weq_from_mass` are reporting-only.
!!
!! **Pressure is in Pa, always.**  Three of the source papers work in dbar
!! and one in metres of depth; the liquidus pressure coefficient differs by
!! 1e4 between Pa and dbar, which is the kind of error that produces a
!! plausible-looking answer.
!!
!! THE THREE EQUATIONS solved here, in `m_mass` form:
!!
!!   (E1) liquidus  `T_b = lambda1*S_b + lambda2 + lambda3*p_b`
!!   (E2) heat      `m_mass*L_eff(T_b) = q_ocean - q_ice`,
!!                  `q_ocean = rho_w*c_w*gamma_t*(T_w - T_b)`
!!   (E3) salt      `m_mass*(S_b - S_i) = rho_w*gamma_s*(S_w - S_b)`
!!
!! (E1) Holland & Jenkins (1999) eq. (1) p. 1788; Jenkins, Nicholls & Corr
!! (2010) eq. (4) p. 2300; Asay-Davis et al. (2016) eq. (25) p. 2485.
!! (E2) Holland & Jenkins (1999) eqs. (2), (3), (9) pp. 1788-1791;
!! Asay-Davis et al. (2016) eq. (24) p. 2485.
!! (E3) Holland & Jenkins (1999) eqs. (4), (5), (10) pp. 1788-1791;
!! Asay-Davis et al. (2016) eq. (26) p. 2485.
!!
!! **THE LIQUIDUS IS NOT DUPLICATED HERE.**  (E1) is evaluated by
!! `eos_freezing_point(eos, S, p)` and its three coefficients are read off
!! the shared `eos_t` handle (`tfr_s`/`tfr_0`/`tfr_p`), selected as a named
!! set by `&ocean_eos_nml tfreeze_set` (`"seaice"` — the SIS2 form — or
!! `"isomip"`).  The two sets are 0.03 degC apart at S = 34.5, which is
!! enough to flip the SIGN of the melt rate over a 0.03 degC band of
!! far-field temperature, so a cavity run configures the set deliberately.
!! No liquidus coefficient is hard-coded in this module.
!!
!! FAIL-LOUD WITHOUT SPEAKING.  Every solver returns a `CAVITY_MELT_*`
!! status instead of `error stop`: a `pure` `!$acc routine seq` procedure
!! can neither log nor abort, and a column kernel must not take the whole
!! run down for one bad column.  On any non-OK status the outputs are set
!! to the documented SAFE STATE — `m_mass = 0` exactly, `S_b = S_w`,
!! `T_b = T_f(S_w, p_b)` — and the status is returned so the caller can
!! count the failures and fail loud itself.  Zero melt is safe; a
!! finite-but-wrong melt is not.
!!
!! NaN IS NEVER LAUNDERED.  Under nvfortran's relaxed-FP default (`-fast`,
!! no `-Kieee`) an `if (x > hi) x = hi` clamp lowers to a NaN-blind min/max
!! select, so a NaN comes OUT as the clamp bound — corruption silently
!! becomes a plausible extreme value.  Every clamp in this module
!! (`u_star`'s floor, the Yung et al. (2025) `min()` caps, the H&J99
!! `eta* := 1` switch) is therefore preceded by an `ieee_is_finite` test of
!! its own input.  Comparisons with NaN are FALSE, so `if (x > thresh)`
!! guards also silently SKIP a NaN — every one here is paired with the
!! finite test.
!!
!! PROVENANCE.  Transliterated from the Python prototype
!! `python_prototypes/ice_shelf_melt/melt.py` (which is the golden oracle's
!! generator), equation by equation, with every `# FORTRAN-GUARD:` marker
!! in that file realised as an `ieee_is_finite` test here.  House rule: the
!! citations below are to the PAPERS — no other model's source was opened.
!!
!!   * Holland, D. M. and Jenkins, A. (1999): "Modeling thermodynamic
!!     ice-ocean interactions at the base of an ice shelf."  J. Phys.
!!     Oceanogr. 29, 1787-1800.
!!   * Jenkins, A., Nicholls, K. W. and Corr, H. F. J. (2010):
!!     "Observation and parameterization of ablation at the base of Ronne
!!     Ice Shelf, Antarctica."  J. Phys. Oceanogr. 40, 2298-2312.
!!   * Asay-Davis, X. S. et al. (2016): "Experimental design for three
!!     interrelated marine ice sheet and ocean model intercomparison
!!     projects."  Geosci. Model Dev. 9, 2471-2497.  (ISOMIP+)
!!   * McPhee, M. G., Maykut, G. A. and Morison, J. H. (1987): "Dynamics
!!     and thermodynamics of the ice/upper ocean system in the marginal ice
!!     zone of the Greenland Sea."  J. Geophys. Res. 92(C7), 7017-7031.
!!   * Yung, C. K. et al. (2025): "Sensitivity of Antarctic ice shelf melt
!!     to the ice-ocean boundary layer parameterisation."  The Cryosphere
!!     19, 5827-5861.
module rdb_ocean_cavity_melt
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_eos, only: eos_t, eos_freezing_point
   implicit none
   private

   public :: ocean_cavity_const_t
   public :: ocean_cavity_exchange_t
   public :: ocean_cavity_ice_t
   public :: ocean_cavity_solution_t
   public :: parse_cavity_exchange_law
   public :: parse_cavity_ice_mode
   public :: cavity_ustar
   public :: cavity_exchange_velocities
   public :: cavity_three_equation
   public :: cavity_two_equation
   public :: cavity_solve_melt
   public :: cavity_melt_point
   public :: cavity_melt_columns
   public :: cavity_heat_fluxes
   public :: cavity_salt_fluxes
   public :: cavity_buoyancy_flux
   public :: cavity_obukhov_length
   public :: cavity_l_plus_from_state
   public :: cavity_l_plus_is_neutral
   public :: cavity_m_ice_from_mass
   public :: cavity_m_weq_from_mass

   ! ======================================================================
   ! Solver status codes.  Returned, never logged — see the module header.
   ! ======================================================================

   integer, parameter, public :: CAVITY_MELT_OK = 0
      !! Solved.  The outputs are physics.
   integer, parameter, public :: CAVITY_MELT_NONFINITE_INPUT = 1
      !! An INPUT was NaN or +/-Inf.  Every guard fires BEFORE any
      !! min/max/clamp, so a non-finite input can never come back out as a
      !! plausible clamp bound (the nvfortran relaxed-FP hazard).
   integer, parameter, public :: CAVITY_MELT_NONFINITE_STATE = 2
      !! An INTERMEDIATE went non-finite (overflow in the discriminant, a
      !! non-finite root or melt rate) from finite inputs.
   integer, parameter, public :: CAVITY_MELT_BAD_INPUT = 3
      !! Finite but outside the solver's domain: a negative drag or
      !! friction-velocity floor, a negative exchange velocity, a negative
      !! ice salinity, `u_star <= 0`, `gamma_s <= 0` in the
      !! three-equation form, or `S_w <= S_i` (which is what the
      !! root-bracketing argument rests on).
   integer, parameter, public :: CAVITY_MELT_NO_PHYSICAL_ROOT = 4
      !! The quadratic has no admissible root: a negative discriminant
      !! (impossible for a well-posed set, so it is treated as corruption
      !! and never clamped to zero), a fully degenerate `A = B = 0`, a
      !! root that failed `S_b > S_i`, or a non-positive effective latent
      !! heat in the two-equation form.
   integer, parameter, public :: CAVITY_MELT_NOT_CONVERGED = 5
      !! The outer stratification bisection did not bracket or did not
      !! converge within `CAVITY_MAX_ITER`.
   integer, parameter, public :: CAVITY_MELT_NOT_IMPLEMENTED = 6
      !! A RESERVED exchange law or ice-conduction mode — the enum value
      !! exists so the dispatch seam is stable, the physics does not ship
      !! yet.  See the per-enum `!!` notes for what the prototype has.
   integer, parameter, public :: CAVITY_MELT_LAW_INVALID = 7
      !! The law / ice-mode code is not one of the enum values at all
      !! (e.g. `CAVITY_LAW_INVALID` straight out of a mistyped namelist
      !! string).  Distinct from NOT_IMPLEMENTED: this one is a typo.
   integer, parameter, public :: CAVITY_MELT_NO_CORIOLIS = 8
      !! The Holland & Jenkins (1999) law was asked for at `f = 0`.  Their
      !! eq. (15) p. 1792 takes `ln(u* xi_N eta*^2 / (|f| h_nu))` and
      !! eq. (18) divides by `f L_O`, so the law simply does not exist on
      !! the equator — the same fail-loud stance this repository already
      !! takes for the Henyey background-mixing latitude factor on a
      !! cartesian grid.
   integer, parameter, public :: CAVITY_MELT_LAW_DOMAIN = 9
      !! An exchange law was evaluated outside its own domain: a
      !! non-positive logarithm argument, a non-positive `eta*` argument,
      !! or a non-positive `Gamma_Turb + Gamma_Mole` denominator (H&J99's
      !! `Gamma_Turb` can go strongly negative at very small `u*`; the
      !! molecular terms normally dominate, but a negative exchange
      !! velocity is refused rather than returned).

   ! ======================================================================
   ! Exchange-law enum.  Codes 1-2 and 4 SHIP; the rest are RESERVED —
   ! their values are nailed down now so that adding one later is not a
   ! renumbering, and `cavity_exchange_velocities` returns
   ! `CAVITY_MELT_NOT_IMPLEMENTED` for them.  Every reserved law IS
   ! implemented in the Python prototype (`ice_shelf_melt/melt.py`), whose
   ! `CITATIONS.md` records the equation numbers and the published
   ! ambiguities each one carries.
   ! ======================================================================

   integer, parameter, public :: CAVITY_LAW_INVALID = 0
      !! Unrecognised law string.  `parse_cavity_exchange_law` returns
      !! this rather than falling back to a default — a mistyped exchange
      !! law is a silent physics change, exactly like a mistyped liquidus.
   integer, parameter, public :: CAVITY_LAW_CONST_GAMMA = 1
      !! **SHIPS.**  `gamma = Gamma*u*` — Jenkins, Nicholls & Corr (2010)
      !! eqs. (1), (2), (5) p. 2300; the ISOMIP+ form, Asay-Davis et al.
      !! (2016) eqs. (24), (26) p. 2485.  The recommended default.
   integer, parameter, public :: CAVITY_LAW_HJ99 = 2
      !! **SHIPS.**  Holland & Jenkins (1999) eqs. (14)-(18) p. 1792 —
      !! turbulent + molecular sublayer with the McPhee (1981) stability
      !! parameter `eta*`.  Stratification-dependent, hence implicit in
      !! the melt rate (outer bisection).
   integer, parameter, public :: CAVITY_LAW_JENKINS91 = 3
      !! RESERVED.  Kader & Yaglom smooth-wall form, Holland & Jenkins
      !! (1999) eqs. (11)-(12) p. 1792.  Explicit in the melt rate, but it
      !! carries a free boundary-layer thickness `h` for which the paper
      !! gives no value.
   integer, parameter, public :: CAVITY_LAW_YUNG25 = 4
      !! **SHIPS.**  Yung et al. (2025) "StratFeedback", eqs. (7)-(8)
      !! p. 5832 — two power laws in the viscous Obukhov scale `L+`,
      !! capped at the Vreugdenhil & Taylor (2019) passive-scalar maxima,
      !! so it reduces EXACTLY to a constant-Gamma law in the
      !! shear-dominated limit.  Stratification-dependent (implicit).
   integer, parameter, public :: CAVITY_LAW_ROSEVEAR22 = 5
      !! RESERVED.  Rosevear, Gayen & Galton-Fenzi (2022) eqs. (27)-(28)
      !! p. 2601.  Needs two documented choices the paper does not make
      !! (natural vs base-10 logarithm; behaviour below `L+ = 2500`).
   integer, parameter, public :: CAVITY_LAW_VT19 = 6
      !! RESERVED.  Vreugdenhil & Taylor (2019) Monin-Obukhov recipe.
      !! Carries a free reference height `z_inf` with no defensible
      !! default, and its authors disown their own salt branch.
   integer, parameter, public :: CAVITY_LAW_MK18 = 7
      !! RESERVED.  McConnochie & Kerr convective floor, reachable only
      !! through Yung et al. (2025) eqs. (10)-(13) p. 5836.  It is the
      !! ONLY law that reads the interface state `(T_b, S_b)` — which is
      !! why those two arguments are in
      !! `cavity_exchange_velocities`' argument list already.
   integer, parameter, public :: CAVITY_LAW_BURCHARD22 = 8
      !! RESERVED.  Burchard et al. (2022) resolution-robust log-layer
      !! law — the only law in the set that treats the far-field sampling
      !! depth as a physical input, and the natural control for a
      !! vertical-coordinate study.
   integer, parameter, public :: CAVITY_LAW_JENKINS21 = 9
      !! RESERVED, and probably permanently: Jenkins (2021) is a 1-D
      !! boundary-CURRENT model, not a per-column transfer law.  Its only
      !! separable piece is the two-equation constant-Stanton law already
      !! available as `cavity_two_equation`.

   ! ======================================================================
   ! Ice-side heat-conduction enum.
   ! ======================================================================

   integer, parameter, public :: CAVITY_ICE_INVALID = 0
      !! Unrecognised ice-conduction string.
   integer, parameter, public :: CAVITY_ICE_INSULATING = 1
      !! **SHIPS.**  Perfect insulator, `q_ice = 0` — Holland & Jenkins
      !! (1999) section 2d(1) p. 1793.  PRESCRIBED by the ISOMIP+ protocol
      !! (Asay-Davis et al. 2016 Table 4 p. 2483 sets `kappa_i = 0`, and
      !! p. 2485 explicitly instructs participants NOT to use the H&J99
      !! advection-diffusion scheme).  `T_ice` is then unread.
   integer, parameter, public :: CAVITY_ICE_ADV_DIFF = 2
      !! **SHIPS.**  Holland & Jenkins (1999) constant-vertical-advection
      !! + diffusion in the linearised eq. (31) p. 1794 form, where the
      !! amplification factor is replaced by its asymptote `Pi = w_I
      !! H_I/kappa_I` for MELTING and `Pi = 0` for FREEZING.  Substituting
      !! (31) into (26)+(6) cancels `H_I` and `kappa_I` identically and
      !! collapses the whole ice column to
      !!
      !!   `q_ice = m_mass*c_i*(T_b - T_ice)`  (melting),  `0` (freezing)
      !!
      !! i.e. the melt flux must additionally warm the incorporated ice
      !! from `T_ice` to `T_b`.  It is LINEAR in `m_mass`, so the
      !! closed-form quadratic survives — it merely replaces `L_f` by an
      !! effective latent heat `L_eff = L_f + c_i*(T_b - T_ice)`.  The ice
      !! thickness and diffusivity are NOT needed in this mode.  The
      !! zeroing on freezing IS what the paper prescribes (eq. 31), not an
      !! approximation added here: p. 1796, "The model with constant
      !! vertical heat advection in the ice shelf has no effect unless the
      !! mixed layer is warmer than the freezing point."
   integer, parameter, public :: CAVITY_ICE_DIFFUSIVE = 3
      !! RESERVED.  Holland & Jenkins (1999) section 2d(2) eq. (21)
      !! p. 1793, `q_ice = k_ice*(T_b - T_ice)/h_ice`: steady, no
      !! advection, linear ice profile.  Not shipped because it also
      !! changes the BRANCH LOGIC — `q_ice` is then independent of
      !! `m_mass`, so `sign(m) = sign(T*)` no longer holds and a strictly
      !! POSITIVE thermal driving is required for zero melt (their
      !! p. 1796: "causes a net shift toward freezing").  Wiring it up
      !! means revisiting the pre-solve melt/freeze branch AND adding
      !! `k_ice`/`h_ice` to the parameter bundles.

   ! ======================================================================
   ! Named constants.  House rule: cite the PAPER.
   ! ======================================================================

   real(wp), parameter, public :: CAVITY_GAMMA_T_ISOMIP = 2.2e-2_wp
      !! ISOMIP+ heat-transfer coefficient `Gamma_T`, dimensionless —
      !! Asay-Davis et al. (2016) section 3.2.1 p. 2487, which derives it
      !! from the Stanton number `sqrt(C_D,top)*Gamma_T = 1.1e-3`
      !! "suggesting that Gamma_T = 2.2 x 10^-2 might be a good initial
      !! guess".  **A STARTING GUESS, NOT A CONSTANT OF NATURE:** the
      !! protocol has participants TUNE it until the Ocean0 mean melt
      !! lands in the prescribed band, and Yung et al. (2026) Table 2
      !! p. 2058 shows the twelve ISOMIP+ submissions landing anywhere
      !! from 0.011 to 0.2.  A namelist knob in the coupling PR, never a
      !! hard-wired number.
   real(wp), parameter, public :: CAVITY_GAMMA_S_ISOMIP = CAVITY_GAMMA_T_ISOMIP/35.0_wp
      !! ISOMIP+ salt-transfer coefficient, `Gamma_S = Gamma_T/35` —
      !! Asay-Davis et al. (2016) Table 4 p. 2483.  The 35 is Jenkins,
      !! Nicholls & Corr (2010) p. 2309: the ratio "should lie somewhere
      !! in the range 35-70.  Adopting a value at the lower end of this
      !! range...".
   real(wp), parameter, public :: CAVITY_CD_ISOMIP = 2.5e-3_wp
      !! ISOMIP+ top drag coefficient `C_D,top` — Asay-Davis et al. (2016)
      !! Table 4 p. 2483.  **The least constrained number in the whole
      !! subject:** the literature spans 1.5e-3 (Holland & Jenkins 1999
      !! Table 1 p. 1790) to 9.7e-3 (Jenkins, Nicholls & Corr 2010 Table 2
      !! p. 2309) — a factor of 6.5 — and Yung et al. (2025) p. 5831
      !! reports order-of-magnitude variation within a single crevasse.
      !! Must be a namelist knob.
   real(wp), parameter, public :: CAVITY_U_TIDE_ISOMIP = 1.0e-2_wp
      !! ISOMIP+ RMS tidal velocity `u_tidal` (m/s) entering the MELT
      !! friction velocity — Asay-Davis et al. (2016) Table 4 p. 2483 and
      !! eq. (27) p. 2485, after Jenkins, Nicholls & Corr (2010) eq. (10)
      !! p. 2309 `u*^2 = C_d (U^2 + <U_T^2>)`.  TRAP: the protocol applies
      !! it to the melt `u*` ONLY, not to the momentum drag (p. 2486, "The
      !! computation of top and bottom drag do not incorporate utidal").
   real(wp), parameter, public :: CAVITY_USTAR_MIN_YUNG25 = 1.0e-4_wp
      !! Friction-velocity floor (m/s) — Yung et al. (2025) eq. (14)
      !! p. 5836 with the value from their Table 2 p. 5838.  It exists
      !! because "a friction velocity of zero (perhaps created by
      !! initialising the model at rest) will result in identically zero
      !! melt ... which would be inconsistent with the presence of heat
      !! available for melting" (p. 5836).

   real(wp), parameter, public :: CAVITY_L_PLUS_NEUTRAL = huge(1.0_wp)
      !! Sentinel for "neutral / unsuppressed", i.e. the `L+ -> +infinity`
      !! limit, carried in the same scalar as a real `L+` so that the
      !! exchange-law seam is ONE argument list.  Any `l_plus` that is
      !! non-finite, non-positive or `>= CAVITY_L_PLUS_NEUTRAL` is treated
      !! as neutral (see `cavity_l_plus_is_neutral`), which is exactly
      !! what every stratification-dependent law prescribes for a
      !! destabilising or vanishing buoyancy flux.

   ! ---- Yung et al. (2025) "StratFeedback" fit, Table 1 + eqs. (7)-(8) p. 5832 ----
   real(wp), parameter :: Y25_A_T = -3.21_wp
      !! `Gamma_T = 10^A_T * (L+)^n_T` prefactor exponent, eq. (7) p. 5832.
   real(wp), parameter :: Y25_N_T = 0.322_wp
      !! `Gamma_T` power-law slope in `L+`, eq. (7) p. 5832.
   real(wp), parameter :: Y25_A_S = -4.30_wp
      !! `Gamma_S = 10^A_S * (L+)^n_S` prefactor exponent, eq. (8) p. 5832.
   real(wp), parameter :: Y25_N_S = 0.223_wp
      !! `Gamma_S` power-law slope in `L+`, eq. (8) p. 5832.
   real(wp), parameter :: Y25_GAMMA_T_CC = 0.012_wp
      !! Constant-coefficient cap on `Gamma_T` (Yung et al. 2025 Table 1
      !! p. 5832) — the Vreugdenhil & Taylor (2019) passive-scalar
      !! maximum.  NOTE this law's neutral limit is its OWN cap, not the
      !! configured `Gamma_T`; `par%gamma_t_coeff` is ignored by design.
   real(wp), parameter :: Y25_GAMMA_S_CC = 3.9e-4_wp
      !! Constant-coefficient cap on `Gamma_S`, Yung et al. (2025) Table 1
      !! p. 5832.  The forced crossover to the caps is APPROXIMATE: with
      !! the published exponents the power laws reach 0.011967 and
      !! 3.876e-4 at `L+ = 1e4`, so the `min()` actually engages at
      !! `L+ = 1.04e4` (heat) and `1.13e4` (salt).  Immaterial
      !! physically; material to any test that asserts equality AT 1e4.

   ! ---- Holland & Jenkins (1999) molecular sublayer, eq. (16) p. 1792 ----
   real(wp), parameter :: HJ99_MOLE_SLOPE = 12.5_wp
      !! `Gamma_Mole = 12.5*(Pr,Sc)^(2/3) - 6`, attributed by Holland &
      !! Jenkins (1999) eq. (16) p. 1792 to Kader & Yaglom (1972).
   real(wp), parameter :: HJ99_MOLE_OFFSET = 6.0_wp
      !! The `- 6` of the same equation.  (Malyarenko et al. (2020)
      !! Table B.1 records this additive constant appearing as -10.1,
      !! -8.68, -9 and -6 across the lineage; we use H&J99's own
      !! rendering, which is what the oracle was generated with.)
   real(wp), parameter :: HJ99_H_NU_COEFF = 5.0_wp
      !! Viscous sublayer thickness `h_nu = 5*nu/u*`, Holland & Jenkins
      !! (1999) eq. (17) p. 1792 — their hydraulically-smooth replacement
      !! for McPhee, Maykut & Morison (1987) eq. (10) p. 7029's roughness
      !! length `z_0`.

   ! ---- Outer stratification iteration.  OURS, and defensible:  no paper ----
   ! gives a tolerance, an iteration cap, a relaxation factor or a
   ! non-convergence fallback (not Yung et al. 2025, not Rosevear et al.
   ! 2022, not Vreugdenhil & Taylor 2019 — every one says "iterate" and
   ! stops).  Bisection on `x = ln(L+)` is used rather than Newton or
   ! secant because the map has a KINK where the buoyancy flux changes
   ! sign (`eta* := 1` in H&J99; the ConstCoeff fallback in Yung 2025) and
   ! another wherever a `min()` cap engages, and no paper proves
   ! uniqueness.  Bisection on a bracket with a sign change is
   ! unconditionally convergent; Newton on a kinked map is not.
   real(wp), parameter :: CAVITY_LP_X_LO = -18.420680743952367_wp
      !! Lower bracket, `ln(1e-8)`.  As `x -> x_lo` the exchange
      !! velocities go to zero, so `|B_b| -> 0` and the re-diagnosed
      !! `L+ -> +infinity`: the residual `G = ln(L+_new) - x` is `> 0`.
      !! (Spelled as a literal, not `log(1.0e-8_wp)`, because a
      !! transcendental is not a Fortran constant expression.)
   real(wp), parameter :: CAVITY_LP_X_HI = 46.051701859880914_wp
      !! Upper bracket, `ln(1e20)`, and the NEUTRAL evaluation point: at
      !! or above it the trial `L+` is taken to be `+infinity`, every law
      !! sits at its own neutral limit, `|B_b|` is maximal and `L+_new` is
      !! finite, so `G < 0`.
   real(wp), parameter :: CAVITY_LP_X_TOL = 1.0e-13_wp
      !! Bracket width in `ln(L+)` at which the bisection is declared
      !! converged.  At `x ~ 46` the double-precision spacing is ~7e-15,
      !! so this is ~15 representable steps above the floor — tight
      !! enough that the returned `L+` is good to ~1e-13 relative, loose
      !! enough to terminate on every toolchain.
   integer, parameter :: CAVITY_MAX_ITER = 200
      !! Hard cap on the outer bisection.  A 64.5-wide bracket halved 200
      !! times is far beyond exhausting double precision, so the loop
      !! normally exits on `CAVITY_LP_X_TOL` or on `mid <= lo`; the cap
      !! exists so a corrupted residual cannot spin forever on device.

   ! ======================================================================
   ! Parameter bundles.  All three are FLAT PODs — no allocatable
   ! component — so they copy into registers when passed BY VALUE into a
   ! `!$acc routine seq` kernel, exactly like `eos_t`.  Under
   ! `mem:separate` that means a host assignment at configure needs no
   ! `!$acc update device`.
   ! ======================================================================

   type :: ocean_cavity_const_t
      !! Thermodynamic + turbulence constants.  All SI.
      !!
      !! DENSITY TRAP: `rho_w` multiplies the OCEAN-side turbulent fluxes
      !! in (E2)/(E3); `rho_i` and `rho_fw` are used ONLY to convert the
      !! canonical mass flux into a thickness rate for REPORTING.  They
      !! are three different numbers in three different papers and must
      !! never be interchanged.
      !!
      !! Defaults are the ISOMIP+ protocol set (Asay-Davis et al. (2016)
      !! Table 4 p. 2483) plus the turbulence constants that protocol does
      !! not specify, taken from Holland & Jenkins (1999) Table 1 p. 1790
      !! and Yung et al. (2025) Table A1 p. 5849.
      real(wp) :: L_f = 3.34e5_wp
         !! Latent heat of fusion (J/kg) — Holland & Jenkins (1999)
         !! Table 1 p. 1790; Asay-Davis et al. (2016) Table 4 p. 2483;
         !! Yung et al. (2025) Table 1 p. 5832.  (Burchard et al. (2022)
         !! Table 1 p. 8 uses 3.335e5.)
      real(wp) :: c_w = 3974.0_wp
         !! Specific heat capacity of seawater (J/kg/K) — unanimous across
         !! Holland & Jenkins (1999), Jenkins et al. (2010), Asay-Davis et
         !! al. (2016) and Yung et al. (2025).
      real(wp) :: c_i = 2009.0_wp
         !! Specific heat capacity of ice (J/kg/K) — Holland & Jenkins
         !! (1999) Table 1 p. 1790.  Read only by `CAVITY_ICE_ADV_DIFF`.
      real(wp) :: rho_w = 1028.0_wp
         !! Seawater density multiplying the turbulent fluxes (kg/m^3) —
         !! Asay-Davis et al. (2016) p. 2479.  The papers span 1025 (H&J99)
         !! to 1030 (Jenkins et al. 2010), a 0.5% spread that lands
         !! directly on the melt rate.
      real(wp) :: rho_i = 918.0_wp
         !! Ice density (kg/m^3), REPORTING ONLY — Asay-Davis et al.
         !! (2016) p. 2479 and Yung et al. (2025).
      real(wp) :: rho_fw = 1000.0_wp
         !! Freshwater density (kg/m^3), REPORTING ONLY — the density
         !! ISOMIP+ reports its `m_w` melt rate with (Asay-Davis et al.
         !! (2016) eq. (24) p. 2485).
      real(wp) :: alpha_T = 3.733e-5_wp
         !! FRACTIONAL thermal expansion coefficient (1/degC) of the
         !! ISOMIP+ linear EOS — Asay-Davis et al. (2016) Table 4 p. 2483.
         !! Feeds the interfacial buoyancy flux only.  NOTE
         !! `eos_t%alpha_T` in this repository is DIMENSIONAL (kg/m^3 per
         !! degC) = `rho0` times this; convert at the coupling seam.
         !! OPEN QUESTION, recorded not resolved: near -2 degC the true
         !! thermal expansion under a nonlinear EOS is near zero or
         !! NEGATIVE, which flips the sign of the (small) temperature term
         !! in the buoyancy flux and therefore of the stratification
         !! feedback near neutrality.  No paper in the set addresses it.
      real(wp) :: beta_S = 7.843e-4_wp
         !! FRACTIONAL haline contraction coefficient (1/(g/kg)) of the
         !! ISOMIP+ linear EOS — Asay-Davis et al. (2016) Table 4 p. 2483.
      real(wp) :: g = 9.81_wp
         !! Gravitational acceleration (m/s^2), as used by the buoyancy
         !! flux — Asay-Davis et al. (2016) Table 4 p. 2483.
      real(wp) :: nu = 1.95e-6_wp
         !! Kinematic viscosity of seawater (m^2/s) — Holland & Jenkins
         !! (1999) Table 1 p. 1790.
      real(wp) :: Pr = 13.8_wp
         !! Molecular Prandtl number — Holland & Jenkins (1999) Table 1
         !! p. 1790, confirmed by McPhee, Maykut & Morison (1987) p. 7029.
      real(wp) :: Sc = 2432.0_wp
         !! Molecular Schmidt number — same two sources.
      real(wp) :: kappa_vk = 0.40_wp
         !! Von Karman constant — Holland & Jenkins (1999) Table 1
         !! p. 1790.  Also the `kappa` of the Obukhov scale.
      real(wp) :: xi_N = 0.052_wp
         !! McPhee stability constant `xi_N` — Holland & Jenkins (1999)
         !! Table 1 p. 1790; McPhee, Maykut & Morison (1987) p. 7029; Yung
         !! et al. (2025) Table A1 p. 5849.  (NOT 0.13: no paper in the
         !! set contains a `zeta_N` or that value.)
      real(wp) :: R_c = 0.20_wp
         !! Critical flux Richardson number — Holland & Jenkins (1999)
         !! Table 1 p. 1790.
   end type ocean_cavity_const_t

   type :: ocean_cavity_exchange_t
      !! Exchange-law selector + its parameters.  ONE bundle for every
      !! law, so the dispatch is a single `select case` and a later
      !! `do concurrent` kernel can call it without reshaping.  Members a
      !! given law does not read are ignored.
      integer :: law = CAVITY_LAW_CONST_GAMMA
         !! `CAVITY_LAW_*`.  Default = the recommended `const_gamma`.
      real(wp) :: gamma_t_coeff = CAVITY_GAMMA_T_ISOMIP
         !! Dimensionless `Gamma_T` of `gamma_t = Gamma_T*u*`.  Named
         !! `*_coeff` to keep it distinct from the exchange VELOCITY
         !! `gamma_t` (m/s) the laws return.
      real(wp) :: gamma_s_coeff = CAVITY_GAMMA_S_ISOMIP
         !! Dimensionless `Gamma_S`.
      real(wp) :: f_cor = -1.4e-4_wp
         !! Coriolis parameter (1/s) — read by `CAVITY_LAW_HJ99` only,
         !! and only as `|f|`.  Holland & Jenkins (1999) Table 1 p. 1790
         !! prints `f = -1.0e-4` (Southern Hemisphere), but their eq. (15)
         !! takes `ln(.../f h_nu)` and eq. (18) needs `f L_O > 0`: both
         !! only make sense with the magnitude.  **Neither that paper nor
         !! McPhee, Maykut & Morison (1987) says `|f|`** — it is inferred
         !! here, and recorded as inferred.
   end type ocean_cavity_exchange_t

   type :: ocean_cavity_ice_t
      !! Ice-side conduction selector + its parameter.
      integer :: mode = CAVITY_ICE_INSULATING
         !! `CAVITY_ICE_*`.  Default = insulating, which is what the
         !! ISOMIP+ protocol prescribes.
      real(wp) :: T_ice = -25.0_wp
         !! Ice interior / surface temperature (degC), read by
         !! `CAVITY_ICE_ADV_DIFF` only — Holland & Jenkins (1999) Table 1
         !! p. 1790 uses `T_S ~ -25.0`.  Under `CAVITY_ICE_INSULATING`
         !! this member is IGNORED and the kernel substitutes exactly
         !! zero, so a stale or absent ice temperature cannot leak into an
         !! insulating run.
   end type ocean_cavity_ice_t

   type :: ocean_cavity_solution_t
      !! Everything `cavity_solve_melt` returns: the interface state, the
      !! fluxes a coupling seam needs, and the solver diagnostics.
      !!
      !! DELIBERATELY WITHOUT DEFAULT INITIALISERS.  A derived type that
      !! carries them cannot be a `do concurrent` `local(...)` variable on
      !! gfortran 15 ("LOCAL specifier ... of derived type with default
      !! initializer is not yet supported"), which would bar the coupling
      !! kernel — and this suite's device test — from declaring one per
      !! column.  `cavity_solve_melt` therefore DEFINES every component on
      !! every path (`cavity_solution_reset` first, before any early
      !! return), so `intent(out)` never leaves one undefined.
      real(wp) :: T_b
         !! Interface temperature (degC), on the liquidus by construction.
      real(wp) :: S_b
         !! Interface salinity (g/kg).
      real(wp) :: m_mass
         !! **Canonical** melt mass flux (kg/m^2/s of ice), > 0 melting.
      real(wp) :: gamma_t
         !! Heat exchange velocity actually used (m/s).
      real(wp) :: gamma_s
         !! Salt exchange velocity actually used (m/s).
      real(wp) :: u_star
         !! Friction velocity the solve was given (m/s).
      real(wp) :: b_flux
         !! Interfacial buoyancy flux (m^2/s^3), < 0 stabilising.
      real(wp) :: l_plus
         !! Viscous Obukhov scale (dimensionless), > 0 stabilising;
         !! `CAVITY_L_PLUS_NEUTRAL` for a vanishing buoyancy flux.
      real(wp) :: L_obukhov
         !! Dimensional Obukhov length (m), > 0 stabilising.
      real(wp) :: T_star
         !! Thermal driving `T_w - T_f(S_w, p_b)` (degC).
      real(wp) :: S_star
         !! Haline driving `S_w - S_b` (g/kg).
      real(wp) :: q_ocean
         !! Turbulent heat flux ocean → interface (W/m^2).
      real(wp) :: q_ice
         !! Conductive + ice-warming flux interface → ice (W/m^2).
      real(wp) :: q_latent
         !! Latent heat consumed by the phase change (W/m^2).
      integer :: n_iter
         !! Outer bisection iterations taken (0 for an explicit law or a
         !! destabilising short-circuit).
      logical :: converged
         !! Outer iteration converged (always `.true.` for an explicit
         !! law).
   end type ocean_cavity_solution_t

contains

   ! ======================================================================
   ! Namelist-string parsing.  Both mirror `parse_tfreeze_set`: an
   ! unrecognised string becomes the INVALID code and the CALLER fails
   ! loud, because silently defaulting a mistyped closure is a physics
   ! change with no run-time symptom.
   ! ======================================================================

   pure function parse_cavity_exchange_law(name) result(code)
      !! Translate an exchange-law string into a `CAVITY_LAW_*` code.
      !! RESERVED laws parse successfully — the refusal belongs to
      !! `cavity_exchange_velocities`, which returns
      !! `CAVITY_MELT_NOT_IMPLEMENTED` — so that a typo
      !! (`CAVITY_LAW_INVALID`) and an honest request for unwritten
      !! physics stay distinguishable.
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: key
      integer :: code
      key = trim(adjustl(name))
      select case (key)
      case ("const_gamma", "CONST_GAMMA", "const", "isomip")
         code = CAVITY_LAW_CONST_GAMMA
      case ("hj99", "HJ99", "holland_jenkins99")
         code = CAVITY_LAW_HJ99
      case ("jenkins91", "JENKINS91")
         code = CAVITY_LAW_JENKINS91
      case ("yung25", "YUNG25", "stratfeedback")
         code = CAVITY_LAW_YUNG25
      case ("rosevear22", "ROSEVEAR22")
         code = CAVITY_LAW_ROSEVEAR22
      case ("vt19", "VT19")
         code = CAVITY_LAW_VT19
      case ("mk18", "MK18")
         code = CAVITY_LAW_MK18
      case ("burchard22", "BURCHARD22")
         code = CAVITY_LAW_BURCHARD22
      case ("jenkins21", "JENKINS21")
         code = CAVITY_LAW_JENKINS21
      case default
         code = CAVITY_LAW_INVALID
      end select
   end function parse_cavity_exchange_law

   pure function parse_cavity_ice_mode(name) result(code)
      !! Translate an ice-conduction string into a `CAVITY_ICE_*` code.
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: key
      integer :: code
      key = trim(adjustl(name))
      select case (key)
      case ("insulating", "INSULATING", "none")
         code = CAVITY_ICE_INSULATING
      case ("adv_diff", "ADV_DIFF", "hj99")
         code = CAVITY_ICE_ADV_DIFF
      case ("diffusive", "DIFFUSIVE")
         code = CAVITY_ICE_DIFFUSIVE
      case default
         code = CAVITY_ICE_INVALID
      end select
   end function parse_cavity_ice_mode

   ! ======================================================================
   ! Friction velocity
   ! ======================================================================

   pure subroutine cavity_ustar(u, v, cd, u_tide, ustar_min, u_star, ierr)
      !! Interface friction velocity (m/s),
      !!
      !!   `u* = max( sqrt( cd*(u^2 + v^2 + u_tide^2) ), ustar_min )`
      !!
      !! Jenkins, Nicholls & Corr (2010) eq. (10) p. 2309 introduces the
      !! tidal variance term `u*^2 = C_d (U^2 + <U_T^2>)`; the ISOMIP+
      !! protocol adopts it as eq. (27) p. 2485 with `u_tidal = 0.01 m/s`
      !! RMS.  The floor `ustar_min` is Yung et al. (2025) eq. (14)
      !! p. 5836.  Both exist for the same stated reason: an ocean at rest
      !! under an ice shelf would otherwise melt exactly nothing even with
      !! heat available.  With `u = v = 0`, `cd = 2.5e-3` and
      !! `u_tide = 0.01` this returns 5e-4 m/s, the operative ISOMIP+
      !! floor.
      !!
      !! TRAP: the tidal term belongs to the MELT friction velocity only —
      !! ISOMIP+ p. 2486, "The computation of top and bottom drag do not
      !! incorporate utidal".  Do not reuse this `u*` for the momentum
      !! drag.
      !!
      !! GUARD: the `max()` is reached only after every input has been
      !! PROVEN finite.  Reversing that order is the documented
      !! NaN-laundering bug — a NaN velocity would come back out as
      !! `ustar_min` and produce a small, plausible melt rate forever.
      !$acc routine seq
      real(wp), intent(in) :: u
         !! Far-field velocity component (m/s).
      real(wp), intent(in) :: v
         !! Far-field velocity component (m/s).
      real(wp), intent(in) :: cd
         !! Top drag coefficient (dimensionless), >= 0.
      real(wp), intent(in) :: u_tide
         !! RMS tidal velocity (m/s), >= 0 by squaring.
      real(wp), intent(in) :: ustar_min
         !! Friction-velocity floor (m/s), >= 0.
      real(wp), intent(out) :: u_star
         !! Friction velocity (m/s).  Exactly zero on any non-OK status.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: us

      ierr = CAVITY_MELT_OK
      u_star = 0.0_wp
      if (.not. (ieee_is_finite(u) .and. ieee_is_finite(v) .and. &
                 ieee_is_finite(cd) .and. ieee_is_finite(u_tide) .and. &
                 ieee_is_finite(ustar_min))) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      if (cd < 0.0_wp .or. ustar_min < 0.0_wp) then
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if
      us = sqrt(cd*(u*u + v*v + u_tide*u_tide))
      if (us > ustar_min) then
         u_star = us
      else
         u_star = ustar_min
      end if
   end subroutine cavity_ustar

   ! ======================================================================
   ! Exchange laws
   ! ======================================================================

   pure function cavity_l_plus_is_neutral(l_plus) result(is_neutral)
      !! Is this trial `L+` the neutral / unsuppressed limit?
      !!
      !! TRUE for a non-finite `L+`, for `L+ <= 0` and for the
      !! `CAVITY_L_PLUS_NEUTRAL` sentinel.  A non-positive `L+` means the
      !! interfacial buoyancy flux is DESTABILISING (or zero), which every
      !! stratification-dependent law treats as its neutral limit: Holland
      !! & Jenkins (1999) p. 1792, "If the Obukhov length is negative
      !! (i.e., the buoyancy flux is destabilizing) the stability
      !! parameter is set to 1"; Yung et al. (2025) pp. 5832-5833, "Since
      !! the LES studies we follow do not explore freezing conditions, we
      !! use the ConstCoeff transfer coefficients ... when L+ < 0".
      !!
      !! It also keeps `log()` off a negative argument in the outer
      !! residual — a real excursion, not a hypothetical: for `hj99` at
      !! strong suppression `Gamma_Turb -> infinity` makes the heat and
      !! salt denominators converge, so `gamma_s/gamma_t -> 1`,
      !! `S_b -> S_w`, the stabilising salt term vanishes and the
      !! destabilising temperature term is left holding the sign.
      !$acc routine seq
      real(wp), intent(in) :: l_plus
         !! Viscous Obukhov scale, or `CAVITY_L_PLUS_NEUTRAL`.
      logical :: is_neutral
      is_neutral = (.not. ieee_is_finite(l_plus)) .or. (l_plus <= 0.0_wp) &
                   .or. (l_plus >= CAVITY_L_PLUS_NEUTRAL)
   end function cavity_l_plus_is_neutral

   pure subroutine cavity_gamma_hj99(u_star, l_plus, par, const, gamma_t, gamma_s, ierr)
      !! Holland & Jenkins (1999) eqs. (14)-(18) p. 1792:
      !!
      !!   `gamma_{T,S} = u* / (Gamma_Turb + Gamma_Mole^{T,S})`        (14)
      !!   `Gamma_Turb  = (1/k)*ln( u* xi_N eta*^2/(|f| h_nu) )
      !!                  + 1/(2 xi_N eta*) - 1/k`                     (15)
      !!   `Gamma_Mole  = 12.5*(Pr,Sc)^(2/3) - 6`                      (16)
      !!   `h_nu        = 5 nu/u*`                                     (17)
      !!   `eta*        = ( 1 + xi_N u*/(|f| L_O R_c) )^(-1/2) <= 1`   (18)
      !!
      !! Identical to McPhee, Maykut & Morison (1987) eq. (10) p. 7029
      !! with the roughness length replaced by the viscous-sublayer
      !! thickness (the hydraulically-smooth assumption, H&J99 p. 1792).
      !!
      !! Check value: `u* = 1e-2 m/s`, `eta* = 1`, `|f| = 1e-4` gives
      !! `gamma_T = 1.06e-4 m/s` and `gamma_S/gamma_T = 0.041`, against
      !! H&J99 Table 1's `gamma_T ~ 1.0e-4` and their p. 1797 ratio of
      !! 0.04.  (Table 1's companion `gamma_S ~ 5.05e-7` is a Hellmer &
      !! Olbers constant-coefficient value, NOT this formulation — the two
      !! are not a self-consistent check pair.)
      !!
      !! BRANCH (must be reproduced exactly for oracle agreement): a
      !! destabilising or vanishing buoyancy flux sets `eta* := 1`.
      !!
      !! GUARDS: `|f| = 0` is refused outright (`CAVITY_MELT_NO_CORIOLIS`)
      !! — the law has `|f|` inside a logarithm and divides by it, so it
      !! does not exist on the equator.  A non-positive logarithm
      !! argument, a non-positive `eta*` argument and a non-positive
      !! denominator are all refused rather than clamped.
      !$acc routine seq
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), > 0.
      real(wp), intent(in) :: l_plus
         !! Trial viscous Obukhov scale, or `CAVITY_L_PLUS_NEUTRAL`.
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Exchange bundle; reads `f_cor` only.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: gamma_t
         !! Heat exchange velocity (m/s).
      real(wp), intent(out) :: gamma_s
         !! Salt exchange velocity (m/s).
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: fa, eta, arg, l_obukhov, h_nu, g_turb, den_t, den_s

      ierr = CAVITY_MELT_OK
      gamma_t = 0.0_wp
      gamma_s = 0.0_wp

      if (.not. ieee_is_finite(par%f_cor)) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      fa = abs(par%f_cor)
      if (fa <= 0.0_wp) then
         ierr = CAVITY_MELT_NO_CORIOLIS
         return
      end if

      ! eq. (18), with the H&J99 p. 1792 destabilising branch.  `L+` is
      ! converted back to the DIMENSIONAL Obukhov length the paper uses,
      ! `L_O = L+ * delta_nu = L+ * nu/u*`.
      if (cavity_l_plus_is_neutral(l_plus)) then
         eta = 1.0_wp
      else
         l_obukhov = l_plus*const%nu/u_star
         if ((.not. ieee_is_finite(l_obukhov)) .or. l_obukhov <= 0.0_wp) then
            eta = 1.0_wp
         else
            arg = 1.0_wp + const%xi_N*u_star/(fa*l_obukhov*const%R_c)
            if (.not. ieee_is_finite(arg)) then
               ierr = CAVITY_MELT_NONFINITE_STATE
               return
            end if
            if (arg <= 0.0_wp) then
               ierr = CAVITY_MELT_LAW_DOMAIN
               return
            end if
            eta = arg**(-0.5_wp)
         end if
      end if

      ! eqs. (15)+(17)
      h_nu = HJ99_H_NU_COEFF*const%nu/u_star
      arg = u_star*const%xi_N*eta*eta/(fa*h_nu)
      if (.not. ieee_is_finite(arg)) then
         ierr = CAVITY_MELT_NONFINITE_STATE
         return
      end if
      if (.not. (arg > 0.0_wp)) then
         ierr = CAVITY_MELT_LAW_DOMAIN
         return
      end if
      g_turb = log(arg)/const%kappa_vk + 1.0_wp/(2.0_wp*const%xi_N*eta) &
               - 1.0_wp/const%kappa_vk

      ! eq. (16) + eq. (14)
      den_t = g_turb + (HJ99_MOLE_SLOPE*const%Pr**(2.0_wp/3.0_wp) - HJ99_MOLE_OFFSET)
      den_s = g_turb + (HJ99_MOLE_SLOPE*const%Sc**(2.0_wp/3.0_wp) - HJ99_MOLE_OFFSET)
      if (.not. (ieee_is_finite(den_t) .and. ieee_is_finite(den_s))) then
         ierr = CAVITY_MELT_NONFINITE_STATE
         return
      end if
      if (den_t <= 0.0_wp .or. den_s <= 0.0_wp) then
         ! `Gamma_Turb` can go strongly negative at very small `u*`.  The
         ! molecular terms (65.9 heat, 2254.6 salt) normally dominate;
         ! refuse rather than return a negative exchange velocity.
         ierr = CAVITY_MELT_LAW_DOMAIN
         return
      end if
      gamma_t = u_star/den_t
      gamma_s = u_star/den_s
   end subroutine cavity_gamma_hj99

   pure subroutine cavity_gamma_yung25(u_star, l_plus, gamma_t, gamma_s, ierr)
      !! Yung et al. (2025) "StratFeedback", eqs. (7)-(8) p. 5832:
      !!
      !!   `Gamma_T = min( 10^A_T*(L+)^n_T , Gamma_T,CC )`
      !!   `Gamma_S = min( 10^A_S*(L+)^n_S , Gamma_S,CC )`
      !!
      !! with `gamma = Gamma*u*`.  The caps are the Vreugdenhil & Taylor
      !! (2019) passive-scalar maxima, and the fit is forced toward them
      !! at `L+ = 1e4` so that the law reduces to a constant-Gamma law in
      !! the well-mixed / shear-dominated limit (p. 5834, "StratFeedback
      !! limits to ConstCoeff at high friction velocities and lower
      !! thermal driving").
      !!
      !! BRANCH: `L+ <= 0` (freezing / destabilising) takes the caps —
      !! the paper's own prescription, pp. 5832-5833, because the LES
      !! studies it follows do not explore freezing.
      !!
      !! NOTE this law IGNORES `par%gamma_t_coeff` by design: its neutral
      !! limit is its own published cap, not the configured `Gamma_T`.  A
      !! gate that asserts "reduces to the configured constant-Gamma law"
      !! is therefore the wrong gate; "reduces to its own documented
      !! neutral limit" is the right one.
      !!
      !! GUARD: the two `min()` caps are reached only after the power laws
      !! have been proven finite.
      !$acc routine seq
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), > 0.
      real(wp), intent(in) :: l_plus
         !! Trial viscous Obukhov scale, or `CAVITY_L_PLUS_NEUTRAL`.
      real(wp), intent(out) :: gamma_t
         !! Heat exchange velocity (m/s).
      real(wp), intent(out) :: gamma_s
         !! Salt exchange velocity (m/s).
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: big_gamma_t, big_gamma_s

      ierr = CAVITY_MELT_OK
      if (cavity_l_plus_is_neutral(l_plus)) then
         gamma_t = Y25_GAMMA_T_CC*u_star
         gamma_s = Y25_GAMMA_S_CC*u_star
         return
      end if
      big_gamma_t = 10.0_wp**Y25_A_T*l_plus**Y25_N_T
      big_gamma_s = 10.0_wp**Y25_A_S*l_plus**Y25_N_S
      if (.not. (ieee_is_finite(big_gamma_t) .and. ieee_is_finite(big_gamma_s))) then
         gamma_t = 0.0_wp
         gamma_s = 0.0_wp
         ierr = CAVITY_MELT_NONFINITE_STATE
         return
      end if
      if (big_gamma_t >= Y25_GAMMA_T_CC) big_gamma_t = Y25_GAMMA_T_CC
      if (big_gamma_s >= Y25_GAMMA_S_CC) big_gamma_s = Y25_GAMMA_S_CC
      gamma_t = big_gamma_t*u_star
      gamma_s = big_gamma_s*u_star
   end subroutine cavity_gamma_yung25

   pure subroutine cavity_exchange_velocities(par, const, u_star, l_plus, &
                                              T_w, S_w, T_b, S_b, gamma_t, gamma_s, ierr)
      !! Exchange-velocity dispatch — ONE argument list for every law, so
      !! a later `do concurrent` kernel dispatches with a single
      !! `select case` and no reshaping.  `l_plus` is the viscous Obukhov
      !! scale of the CURRENT iterate: the single scalar that carries the
      !! stratification feedback for every implicit law.  Pass
      !! `CAVITY_L_PLUS_NEUTRAL` for the neutral / unsuppressed
      !! evaluation.
      !!
      !! `T_w`, `S_w`, `T_b`, `S_b` are in the list for `CAVITY_LAW_MK18`
      !! (reserved), the only law whose exchange velocities depend on the
      !! interface state; every law that ships ignores them.  They are
      !! carried now precisely so that adding MK18 later does not
      !! re-signature this seam or any kernel that calls it.
      !$acc routine seq
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Law selector + parameters.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), strictly positive.
      real(wp), intent(in) :: l_plus
         !! Trial viscous Obukhov scale, or `CAVITY_L_PLUS_NEUTRAL`.
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC) — reserved-law argument.
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg) — reserved-law argument.
      real(wp), intent(in) :: T_b
         !! Trial interface temperature (degC) — reserved-law argument.
      real(wp), intent(in) :: S_b
         !! Trial interface salinity (g/kg) — reserved-law argument.
      real(wp), intent(out) :: gamma_t
         !! Heat exchange velocity (m/s).  Zero on any non-OK status.
      real(wp), intent(out) :: gamma_s
         !! Salt exchange velocity (m/s).  Zero on any non-OK status.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.

      ierr = CAVITY_MELT_OK
      gamma_t = 0.0_wp
      gamma_s = 0.0_wp

      if (.not. (ieee_is_finite(u_star) .and. ieee_is_finite(T_w) .and. &
                 ieee_is_finite(S_w) .and. ieee_is_finite(T_b) .and. &
                 ieee_is_finite(S_b))) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      if (u_star <= 0.0_wp) then
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if

      select case (par%law)
      case (CAVITY_LAW_CONST_GAMMA)
         if (.not. (ieee_is_finite(par%gamma_t_coeff) .and. &
                    ieee_is_finite(par%gamma_s_coeff))) then
            ierr = CAVITY_MELT_NONFINITE_INPUT
            return
         end if
         gamma_t = par%gamma_t_coeff*u_star
         gamma_s = par%gamma_s_coeff*u_star
      case (CAVITY_LAW_HJ99)
         call cavity_gamma_hj99(u_star, l_plus, par, const, gamma_t, gamma_s, ierr)
      case (CAVITY_LAW_YUNG25)
         call cavity_gamma_yung25(u_star, l_plus, gamma_t, gamma_s, ierr)
      case (CAVITY_LAW_JENKINS91, CAVITY_LAW_ROSEVEAR22, CAVITY_LAW_VT19, &
            CAVITY_LAW_MK18, CAVITY_LAW_BURCHARD22, CAVITY_LAW_JENKINS21)
         ! Reserved.  Each is implemented in the Python prototype
         ! `ice_shelf_melt/melt.py`, with the published ambiguities it
         ! carries recorded in that directory's CITATIONS.md.
         ierr = CAVITY_MELT_NOT_IMPLEMENTED
      case default
         ierr = CAVITY_MELT_LAW_INVALID
      end select
   end subroutine cavity_exchange_velocities

   ! ======================================================================
   ! The three-equation system
   ! ======================================================================
   !
   ! DERIVATION (ours; no paper states the root selection — Yung et al.
   ! (2025) p. 5830 says only "this system of equations reduces to a
   ! quadratic equation", Burchard et al. (2022) p. 7 says "the positive
   ! root", and the others do not mention the quadratic at all).
   !
   ! Write the liquidus as `T_b = a*S_b + b`, `a = lambda1`,
   ! `b = lambda2 + lambda3*p`, and let `D = T_w - b`, so
   ! `T_w - T_b = D - a*S_b`.
   !
   ! Heat (E2), with `L_eff = E0 + E1*S_b` and `q_ice = kh*(T_b - T_ice)`:
   !     E0 = L_f + c_i_eff*(b - T_ice)        E1 = c_i_eff*a
   !     kh = k_ice/h_ice  (reserved DIFFUSIVE mode)  else 0
   !     c_i_eff = c_i     (ADV_DIFF and melting)     else 0
   !     m*(E0 + E1*S_b) = rho_w*c_w*gamma_t*(D - a*S_b)
   !                       - kh*(a*S_b + b - T_ice)  =  H0 + H1*S_b
   !     H0 = rho_w*c_w*gamma_t*D - kh*(b - T_ice)
   !     H1 = -a*(rho_w*c_w*gamma_t + kh)
   ! Salt (E3):  m*(S_b - S_i) = P1*(S_w - S_b),  P1 = rho_w*gamma_s
   ! Eliminating m (valid for S_b /= S_i) gives A*S_b^2 + B*S_b + C = 0:
   !     A = -P1*E1 - H1
   !     B =  P1*(S_w*E1 - E0) - H0 + H1*S_i
   !     C =  P1*S_w*E0 + H0*S_i
   !
   ! WHICH ROOT IS PHYSICAL — the LARGER one, always.  With
   ! `f(S) = A S^2 + B S + C`:
   !   * `A = a*(rho_w*c_w*gamma_t + kh - rho_w*gamma_s*c_i_eff)`, and
   !     `a = lambda1 < 0` while the bracket is positive for any sane
   !     parameter set, so A < 0 — the parabola opens DOWNWARD;
   !   * `f(S_i) = P1*(S_w - S_i)*L_eff(S_i) > 0` whenever `S_w > S_i`,
   !     so `S_i` lies strictly BETWEEN the two roots.
   !   Hence `r_minus < S_i <= r_plus` and the larger root is the only one
   !   with `S_b >= S_i >= 0`.  This holds for melting (`S_i < S_b < S_w`),
   !   for freezing (`S_b > S_w`) and for `S_i = 0` alike.
   !
   ! MELT/FREEZE BRANCH, decided BEFORE the solve.  For the two shipped
   ! ice modes (`kh = 0`),
   !     f(S_w) = -rho_w*c_w*gamma_t*T_star*(S_w - S_i)
   ! so `sign(S_w - S_b) = sign(T_star) = sign(m)` exactly.  Deciding the
   ! conduction branch from `T_star` up front makes it impossible for the
   ! branch to disagree with its own answer.  (This breaks for the
   ! reserved DIFFUSIVE mode — see `CAVITY_ICE_DIFFUSIVE`.)
   !
   ! NUMERICAL STABILITY — why not `(-B - sqrt(disc))/(2A)`.  The larger
   ! root of a downward parabola is exactly that expression, and it
   ! cancels catastrophically whenever `|B| >> sqrt(4AC)` — which is the
   ! strong-double-diffusion limit `gamma_s/gamma_t -> 0`, i.e. where this
   ! model spends most of its life (`Gamma_S = Gamma_T/35`).  The standard
   ! two-branch form is used instead:
   !     q = -0.5*(B + sign(sqrt(disc), B));  roots are q/A and C/q
   !     B >= 0: q < 0 and q/A is the LARGER root
   !     B <  0: q > 0 and C/q is the LARGER root
   !
   ! A NONLINEAR LIQUIDUS would break the closed form: `T_b = T_f(S_b, p)`
   ! would no longer be affine in `S_b`, so the quadratic would have to be
   ! replaced by a bracketed bisection of
   !     f(S) = rho_w*gamma_s*(S_w - S)*L_eff(S) - (q_ocean(S) - q_ice(S))*(S - S_i)
   ! over `[S_i, S_hi]` (`f(S_i) > 0` by the argument above; grow `S_hi`
   ! geometrically from `S_w` for the freezing branch).  That solver is
   ! implemented in the prototype as `three_equation_bracketed` and agrees
   ! with the closed form to ~1e-14 on a linear liquidus.  It is NOT
   ! ported here because `eos_freezing_point` is linear by construction —
   ! this is where it would go if `&ocean_eos_nml` ever grows a
   ! Millero/TEOS-10 `tfreeze_form`.
   ! ======================================================================

   pure function cavity_t_ice(ice) result(T_ice)
      !! Ice temperature actually used: exactly zero when the mode ignores
      !! it, so an unset or stale `T_ice` cannot leak into an insulating
      !! run through `L_eff`.
      !$acc routine seq
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      real(wp) :: T_ice
      if (ice%mode == CAVITY_ICE_INSULATING) then
         T_ice = 0.0_wp
      else
         T_ice = ice%T_ice
      end if
   end function cavity_t_ice

   pure subroutine cavity_ice_terms(ice, melting, const, c_i_eff, kh, ierr)
      !! `(c_i_eff, kh)` for the requested ice-conduction mode — the two
      !! numbers through which every mode enters the quadratic.
      !!
      !! `CAVITY_ICE_ADV_DIFF` is Holland & Jenkins (1999) eq. (31)
      !! p. 1794: the heat-advection term is present for MELTING and
      !! exactly zero for FREEZING, which is what the paper prescribes
      !! rather than an approximation made here.
      !$acc routine seq
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      logical, intent(in) :: melting
         !! Melt/freeze branch, decided from the thermal driving.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: c_i_eff
         !! Effective ice heat capacity (J/kg/K) in `L_eff = L_f +
         !! c_i_eff*(T_b - T_ice)`.
      real(wp), intent(out) :: kh
         !! Conductance `k_ice/h_ice` (W/m^2/K).  Exactly zero for both
         !! shipped modes; carried so the derivation above and the code
         !! stay the same expression when the reserved DIFFUSIVE mode
         !! lands.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.

      ierr = CAVITY_MELT_OK
      c_i_eff = 0.0_wp
      kh = 0.0_wp
      select case (ice%mode)
      case (CAVITY_ICE_INSULATING)
         ! q_ice == 0 identically.
      case (CAVITY_ICE_ADV_DIFF)
         if (melting) c_i_eff = const%c_i
      case (CAVITY_ICE_DIFFUSIVE)
         ierr = CAVITY_MELT_NOT_IMPLEMENTED
      case default
         ierr = CAVITY_MELT_LAW_INVALID
      end select
   end subroutine cavity_ice_terms

   pure subroutine cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
      !! The documented SAFE STATE returned on every non-OK status: zero
      !! melt, interface salinity equal to the far field, interface
      !! temperature on the liquidus there.  `m_mass` is EXACTLY zero, so
      !! a failed column contributes nothing to a heat/salt budget rather
      !! than contributing a plausible wrong number.
      !!
      !! `T_b`/`S_b` inherit whatever `S_w`/`p_b` were: if the caller
      !! handed in a NaN, a NaN comes back.  That is deliberate — this
      !! routine must not INVENT a finite interface state out of
      !! corrupted input.  The melt flux, the one output that would
      !! silently poison a budget, is the one pinned to zero.
      !$acc routine seq
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — carries the liquidus coefficient set.
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa).
      real(wp), intent(out) :: T_b
         !! Interface temperature (degC) = `T_f(S_w, p_b)`.
      real(wp), intent(out) :: S_b
         !! Interface salinity (g/kg) = `S_w`.
      real(wp), intent(out) :: m_mass
         !! Melt mass flux (kg/m^2/s) = exactly 0.
      T_b = eos_freezing_point(eos, S_w, p_b)
      S_b = S_w
      m_mass = 0.0_wp
   end subroutine cavity_safe_state

   pure subroutine cavity_three_equation(T_w, S_w, p_b, gamma_t, gamma_s, S_i, &
                                         ice, eos, const, T_b, S_b, m_mass, ierr)
      !! Closed-form solve of (E1)-(E3) on the linear liquidus carried by
      !! `eos`.  Returns the interface state and the canonical melt mass
      !! flux (kg/m^2/s, > 0 melting).  See the derivation block above for
      !! the root selection, the cancellation-safe quadratic and the
      !! pre-solve melt/freeze branch.
      !!
      !! `gamma_s > 0` is REQUIRED: the `gamma_s -> infinity` limit is the
      !! two-equation form and has its own entry point
      !! (`cavity_two_equation`), while `gamma_s = 0` is not a limit of
      !! this system at all.
      !$acc routine seq
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).  Must exceed `S_i`.
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa) — `g*rho_i*h_ice`, the undiluted ice
         !! load.
      real(wp), intent(in) :: gamma_t
         !! Heat exchange velocity (m/s), >= 0.
      real(wp), intent(in) :: gamma_s
         !! Salt exchange velocity (m/s), > 0.
      real(wp), intent(in) :: S_i
         !! Ice salinity (g/kg), >= 0.  Holland & Jenkins (1999) p. 1789
         !! treats marine ice as fresh ("we can treat S_I as zero
         !! always"); the argument is kept because the root-bracketing
         !! proof is stated for general `S_i < S_w`.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.  By value, flat POD.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: T_b
         !! Interface temperature (degC).
      real(wp), intent(out) :: S_b
         !! Interface salinity (g/kg).
      real(wp), intent(out) :: m_mass
         !! Melt mass flux (kg/m^2/s), > 0 melting.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.  On anything but OK the three outputs
         !! are the safe state (see `cavity_safe_state`).
      real(wp) :: T_ice, T_star, c_i_eff, kh
      real(wp) :: a, b, d_w, e0, e1, rc, h0, h1, p1, qa, qb, qc, disc, sq, q
      logical :: melting

      call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
      T_ice = cavity_t_ice(ice)

      ierr = CAVITY_MELT_OK
      if (.not. (ieee_is_finite(T_w) .and. ieee_is_finite(S_w) .and. &
                 ieee_is_finite(p_b) .and. ieee_is_finite(gamma_t) .and. &
                 ieee_is_finite(gamma_s) .and. ieee_is_finite(S_i) .and. &
                 ieee_is_finite(T_ice))) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      if (gamma_t < 0.0_wp .or. gamma_s <= 0.0_wp .or. S_i < 0.0_wp) then
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if
      if (.not. (S_w > S_i)) then
         ! The whole root-bracketing argument rests on `f(S_i) > 0`, which
         ! needs `S_w > S_i`.
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if

      T_star = T_w - eos_freezing_point(eos, S_w, p_b)
      melting = T_star > 0.0_wp
      call cavity_ice_terms(ice, melting, const, c_i_eff, kh, ierr)
      if (ierr /= CAVITY_MELT_OK) return

      a = eos%tfr_s
      b = eos%tfr_0 + eos%tfr_p*p_b
      d_w = T_w - b
      e0 = const%L_f + c_i_eff*(b - T_ice)
      e1 = c_i_eff*a
      rc = const%rho_w*const%c_w*gamma_t
      h0 = rc*d_w - kh*(b - T_ice)
      h1 = -a*(rc + kh)
      p1 = const%rho_w*gamma_s

      qa = -p1*e1 - h1
      qb = p1*(S_w*e1 - e0) - h0 + h1*S_i
      qc = p1*S_w*e0 + h0*S_i

      if (qa == 0.0_wp) then
         ! Degenerate: a salinity-independent liquidus (`lambda1 = 0`), or
         ! the pathological `c_w*gamma_t + kh/rho_w == c_i*gamma_s`.  The
         ! system is then linear in `S_b`.
         if (qb == 0.0_wp) then
            ierr = CAVITY_MELT_NO_PHYSICAL_ROOT
            return
         end if
         S_b = -qc/qb
      else
         disc = qb*qb - 4.0_wp*qa*qc
         if (.not. ieee_is_finite(disc)) then
            ierr = CAVITY_MELT_NONFINITE_STATE
            return
         end if
         if (disc < 0.0_wp) then
            ! Cannot happen for a well-posed set (`S_i` lies between the
            ! roots, so `disc > 0`).  Treat it as corruption and refuse —
            ! NEVER clamp it to zero, which would manufacture a double
            ! root and a plausible melt rate out of a broken column.
            ierr = CAVITY_MELT_NO_PHYSICAL_ROOT
            return
         end if
         sq = sqrt(disc)
         if (qb >= 0.0_wp) then
            q = -0.5_wp*(qb + sq)
            S_b = q/qa
         else
            q = -0.5_wp*(qb - sq)
            ! `q == 0` requires `B == 0` and `disc == 0`, excluded above.
            S_b = qc/q
         end if
      end if

      if (.not. ieee_is_finite(S_b)) then
         call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
         ierr = CAVITY_MELT_NONFINITE_STATE
         return
      end if
      if (S_b <= S_i) then
         call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
         ierr = CAVITY_MELT_NO_PHYSICAL_ROOT
         return
      end if

      T_b = eos_freezing_point(eos, S_b, p_b)
      m_mass = const%rho_w*gamma_s*(S_w - S_b)/(S_b - S_i)
      if (.not. (ieee_is_finite(m_mass) .and. ieee_is_finite(T_b))) then
         call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
         ierr = CAVITY_MELT_NONFINITE_STATE
      end if
   end subroutine cavity_three_equation

   pure subroutine cavity_two_equation(T_w, S_w, p_b, gamma_t, ice, eos, const, &
                                       T_b, S_b, m_mass, ierr)
      !! Two-equation variant: the interface salinity is the FAR-FIELD
      !! salinity, `S_b = S_w` exactly, i.e. the `gamma_s -> infinity`
      !! limit of the three-equation form (NOT the `gamma_s -> 0` limit,
      !! which sends `T_b -> T_w` and the melt rate to zero).  Holland &
      !! Jenkins (1999) section 2b(2) p. 1791; Jenkins, Nicholls & Corr
      !! (2010) eq. (6) p. 2302, where the single transfer coefficient is
      !! `Gamma_TS ~ 0.006` and the freezing point is evaluated at the
      !! far-field salinity.
      !!
      !! **IDENTITY WITH THE SHIPPED SEA-ICE BASAL FLUX.**  This form, at
      !! `gamma_t*dt = h`, IS `rdb_ice_basal_flux`'s
      !! "relax the top layer to the freezing point in one thermo step"
      !! law: that kernel forms
      !! `fb = rho*c_p*max(0, SST - T_f)*h/dt_therm`, and setting
      !! `gamma_t = h/dt` in `q_ocean = rho_w*c_w*gamma_t*(T_w - T_b)`
      !! with `T_b = T_f(S_w, p)` reproduces it term for term.  It is an
      !! exact identity, not an asymptote — the two differences are that
      !! the sea-ice kernel clamps the supercooled branch away with a
      !! `max(0, .)` (frazil owns that side) and evaluates the liquidus at
      !! `p = 0`.  The practical reading: the shipped sea-ice law already
      !! has an implied exchange velocity, and it is resolution- and
      !! timestep-dependent by construction, which is exactly the
      !! far-field-sampling problem the cavity literature is about.
      !$acc routine seq
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa).
      real(wp), intent(in) :: gamma_t
         !! Heat exchange velocity (m/s), >= 0.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: T_b
         !! Interface temperature (degC) = `T_f(S_w, p_b)`.
      real(wp), intent(out) :: S_b
         !! Interface salinity (g/kg) = `S_w` exactly.
      real(wp), intent(out) :: m_mass
         !! Melt mass flux (kg/m^2/s), > 0 melting.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: T_ice, c_i_eff, kh, l_eff, q_oc, q_ic
      logical :: melting

      call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
      T_ice = cavity_t_ice(ice)

      ierr = CAVITY_MELT_OK
      if (.not. (ieee_is_finite(T_w) .and. ieee_is_finite(S_w) .and. &
                 ieee_is_finite(p_b) .and. ieee_is_finite(gamma_t) .and. &
                 ieee_is_finite(T_ice))) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      if (gamma_t < 0.0_wp) then
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if

      melting = (T_w - T_b) > 0.0_wp
      call cavity_ice_terms(ice, melting, const, c_i_eff, kh, ierr)
      if (ierr /= CAVITY_MELT_OK) return

      l_eff = const%L_f + c_i_eff*(T_b - T_ice)
      if (l_eff <= 0.0_wp) then
         ierr = CAVITY_MELT_NO_PHYSICAL_ROOT
         return
      end if
      q_oc = const%rho_w*const%c_w*gamma_t*(T_w - T_b)
      q_ic = kh*(T_b - T_ice)
      m_mass = (q_oc - q_ic)/l_eff
      if (.not. ieee_is_finite(m_mass)) then
         call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
         ierr = CAVITY_MELT_NONFINITE_STATE
      end if
   end subroutine cavity_two_equation

   ! ======================================================================
   ! Buoyancy flux, Obukhov scales
   ! ======================================================================

   pure function cavity_buoyancy_flux(const, T_w, S_w, T_b, S_b, gamma_t, gamma_s) &
      result(b_flux)
      !! Interfacial buoyancy flux `B_b` (m^2/s^3), NEGATIVE = stabilising
      !! — Yung et al. (2025) eq. (6) p. 5831,
      !!
      !!   `B_b = -g*( beta*(S_w - S_b)*gamma_s - alpha*(T_w - T_b)*gamma_t )`
      !!
      !! (their form carries `gamma = Gamma*u*` explicitly).  Equivalent to
      !! the McPhee, Maykut & Morison (1987) p. 7029 expression once their
      !! dimensional expansion coefficients are written as `rho_0*beta`,
      !! `rho_0*alpha` and the kinematic fluxes as `gamma*(far - interface)`.
      !!
      !! SIGN CHECK: melting FRESHENS the interface (`S_b < S_w`), which is
      !! stabilising, and COOLS it (`T_b < T_w`), which is destabilising;
      !! the salt term wins by about an order of magnitude at seawater
      !! salinities, so melting gives `B_b < 0` and therefore a POSITIVE
      !! Obukhov length.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `g`, `alpha_T`, `beta_S`.
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: T_b
         !! Interface temperature (degC).
      real(wp), intent(in) :: S_b
         !! Interface salinity (g/kg).
      real(wp), intent(in) :: gamma_t
         !! Heat exchange velocity (m/s).
      real(wp), intent(in) :: gamma_s
         !! Salt exchange velocity (m/s).
      real(wp) :: b_flux
      b_flux = -const%g*(const%beta_S*(S_w - S_b)*gamma_s &
                         - const%alpha_T*(T_w - T_b)*gamma_t)
   end function cavity_buoyancy_flux

   pure function cavity_obukhov_length(const, u_star, b_flux) result(l_obukhov)
      !! Dimensional Obukhov length `L = -u*^3/(kappa*B_b)` (m), POSITIVE
      !! for a stabilising (melting) buoyancy flux.  McPhee, Maykut &
      !! Morison (1987) p. 7029; the same scale appears as Yung et al.
      !! (2025) eq. (5) p. 5831.  (Holland & Jenkins (1999) uses `L_O` in
      !! their eq. (18) but never defines it.)
      !!
      !! Returns `CAVITY_L_PLUS_NEUTRAL` for a vanishing buoyancy flux —
      !! the `L -> +infinity` neutral limit, carried as a finite sentinel
      !! so no downstream arithmetic has to handle an actual infinity.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `kappa_vk`.
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), > 0.
      real(wp), intent(in) :: b_flux
         !! Interfacial buoyancy flux (m^2/s^3).
      real(wp) :: l_obukhov
      if (b_flux == 0.0_wp) then
         l_obukhov = CAVITY_L_PLUS_NEUTRAL
      else
         l_obukhov = -(u_star*u_star*u_star)/(const%kappa_vk*b_flux)
      end if
   end function cavity_obukhov_length

   pure function cavity_l_plus_from_state(const, u_star, b_flux) result(l_plus)
      !! Viscous Obukhov scale `L+ = L/delta_nu` with `delta_nu = nu/u*`,
      !! i.e. `L+ = -u*^4/(nu*kappa*B_b)` — Yung et al. (2025) eq. (5)
      !! p. 5831; the same definition in Vreugdenhil & Taylor (2019)
      !! eq. (27) and Rosevear et al. (2022) eqs. (6)+(8) p. 2592.
      !! POSITIVE for melting.  Returns `CAVITY_L_PLUS_NEUTRAL` for a
      !! vanishing buoyancy flux.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `nu`, `kappa_vk`.
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), > 0.
      real(wp), intent(in) :: b_flux
         !! Interfacial buoyancy flux (m^2/s^3).
      real(wp) :: l_plus
      if (b_flux == 0.0_wp) then
         l_plus = CAVITY_L_PLUS_NEUTRAL
      else
         l_plus = -(u_star**4)/(const%nu*const%kappa_vk*b_flux)
      end if
   end function cavity_l_plus_from_state

   ! ======================================================================
   ! Diagnostic fluxes and conversions
   ! ======================================================================

   pure subroutine cavity_heat_fluxes(T_w, T_b, m_mass, gamma_t, ice, const, &
                                      q_ocean, q_ice, q_latent)
      !! The three heat fluxes of (E2), W/m^2:
      !!
      !!   `q_ocean  = rho_w*c_w*gamma_t*(T_w - T_b)`  ocean → interface
      !!   `q_ice    = kh*(T_b - T_ice) + m_mass*c_i_eff*(T_b - T_ice)`
      !!   `q_latent = m_mass*L_f`
      !!
      !! The closure `q_ocean - q_ice - q_latent = 0` is what the unit
      !! tests assert to round-off, independently of any golden number.
      !!
      !! The melt/freeze branch is taken here from `sign(m_mass)`, which
      !! for both shipped ice modes is identical to the `sign(T_star)`
      !! branch the solve used (see the derivation block) — the two agree
      !! by construction, including at exactly zero melt.
      !!
      !! This is a DIAGNOSTIC re-evaluation, so it carries no status: a
      !! reserved or invalid ice mode leaves `c_i_eff = kh = 0`, i.e. the
      !! insulating fluxes, and the refusal is reported by the SOLVER
      !! that produced `m_mass` in the first place.
      !$acc routine seq
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: T_b
         !! Interface temperature (degC).
      real(wp), intent(in) :: m_mass
         !! Melt mass flux (kg/m^2/s).
      real(wp), intent(in) :: gamma_t
         !! Heat exchange velocity (m/s).
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: q_ocean
         !! Turbulent heat flux ocean → interface (W/m^2).
      real(wp), intent(out) :: q_ice
         !! Conductive + ice-warming flux interface → ice (W/m^2).
      real(wp), intent(out) :: q_latent
         !! Latent heat consumed by the phase change (W/m^2).
      real(wp) :: T_ice, c_i_eff, kh
      integer :: ice_stat

      T_ice = cavity_t_ice(ice)
      call cavity_ice_terms(ice, m_mass > 0.0_wp, const, c_i_eff, kh, ice_stat)
      q_ocean = const%rho_w*const%c_w*gamma_t*(T_w - T_b)
      q_ice = kh*(T_b - T_ice) + m_mass*c_i_eff*(T_b - T_ice)
      q_latent = m_mass*const%L_f
   end subroutine cavity_heat_fluxes

   pure subroutine cavity_salt_fluxes(const, S_w, S_b, m_mass, gamma_s, S_i, &
                                      f_turb, f_phase)
      !! The two sides of (E3), in (g/kg)*kg/m^2/s:
      !!
      !!   `f_turb  = rho_w*gamma_s*(S_w - S_b)`  turbulent salt flux
      !!                                          toward the interface
      !!   `f_phase = m_mass*(S_b - S_i)`         salt rejected/diluted by
      !!                                          the phase change
      !!
      !! (E3) says these are equal; the unit tests assert it to round-off.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `rho_w`.
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: S_b
         !! Interface salinity (g/kg).
      real(wp), intent(in) :: m_mass
         !! Melt mass flux (kg/m^2/s).
      real(wp), intent(in) :: gamma_s
         !! Salt exchange velocity (m/s).
      real(wp), intent(in) :: S_i
         !! Ice salinity (g/kg).
      real(wp), intent(out) :: f_turb
         !! Turbulent salt flux toward the interface.
      real(wp), intent(out) :: f_phase
         !! Phase-change salt flux.
      f_turb = const%rho_w*gamma_s*(S_w - S_b)
      f_phase = m_mass*(S_b - S_i)
   end subroutine cavity_salt_fluxes

   pure function cavity_m_ice_from_mass(const, m_mass) result(m_ice)
      !! Solid-ice thickness rate (m/s) from the canonical mass flux,
      !! `m_ice = m_mass/rho_i`.  REPORTING ONLY — this is Jenkins,
      !! Nicholls & Corr (2010)'s `a_b` convention.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `rho_i`.
      real(wp), intent(in) :: m_mass
         !! Melt mass flux (kg/m^2/s).
      real(wp) :: m_ice
      m_ice = m_mass/const%rho_i
   end function cavity_m_ice_from_mass

   pure function cavity_m_weq_from_mass(const, m_mass) result(m_weq)
      !! Freshwater-equivalent thickness rate (m/s) from the canonical
      !! mass flux, `m_weq = m_mass/rho_fw`.  REPORTING ONLY — this is
      !! ISOMIP+'s `m_w` (Asay-Davis et al. (2016) eq. (24) p. 2485), the
      !! number their figures are in.
      !$acc routine seq
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle — `rho_fw`.
      real(wp), intent(in) :: m_mass
         !! Melt mass flux (kg/m^2/s).
      real(wp) :: m_weq
      m_weq = m_mass/const%rho_fw
   end function cavity_m_weq_from_mass

   ! ======================================================================
   ! The full solve, with the outer stratification iteration
   ! ======================================================================

   pure subroutine cavity_solution_reset(sol)
      !! Define every component of the solution bundle.  Called FIRST by
      !! `cavity_solve_melt`, before any early return, because the type
      !! carries no default initialisers (see its docstring: they would
      !! bar it from a `do concurrent` `local(...)` clause on gfortran).
      !$acc routine seq
      type(ocean_cavity_solution_t), intent(out) :: sol
         !! Solution bundle, zeroed.
      sol%T_b = 0.0_wp
      sol%S_b = 0.0_wp
      sol%m_mass = 0.0_wp
      sol%gamma_t = 0.0_wp
      sol%gamma_s = 0.0_wp
      sol%u_star = 0.0_wp
      sol%b_flux = 0.0_wp
      sol%l_plus = 0.0_wp
      sol%L_obukhov = 0.0_wp
      sol%T_star = 0.0_wp
      sol%S_star = 0.0_wp
      sol%q_ocean = 0.0_wp
      sol%q_ice = 0.0_wp
      sol%q_latent = 0.0_wp
      sol%n_iter = 0
      sol%converged = .false.
   end subroutine cavity_solution_reset

   pure function cavity_law_is_implicit(law) result(is_implicit)
      !! Does this law's `(gamma_t, gamma_s)` depend on the interfacial
      !! buoyancy flux — and therefore on the melt rate it produces?  Yung
      !! et al. (2025) p. 5833 states the consequence: "Since the transfer
      !! coefficients depend on L+, which in turn depends on melt rate via
      !! surface buoyancy forcing, iteration is required for convergence
      !! of the three-equation parameterisation solution."
      !$acc routine seq
      integer, intent(in) :: law
         !! `CAVITY_LAW_*` code.
      logical :: is_implicit
      is_implicit = (law == CAVITY_LAW_HJ99) .or. (law == CAVITY_LAW_YUNG25)
   end function cavity_law_is_implicit

   pure function cavity_outer_residual(lp_new, x) result(g)
      !! Outer-iteration residual `G(x) = ln(L+_new(x)) - x`, with the
      !! destabilising branch folded in: a non-positive, non-finite or
      !! sentinel `L+_new` means the buoyancy flux at this iterate is
      !! destabilising (or zero), which every law treats as
      !! `L+ = +infinity`, so `G = +infinity`.  Mapping it that way keeps
      !! the bisection bracket valid instead of taking `log()` of a
      !! negative number.
      !$acc routine seq
      real(wp), intent(in) :: lp_new
         !! `L+` re-diagnosed from the state this iterate produced.
      real(wp), intent(in) :: x
         !! The trial `ln(L+)` it was produced at.
      real(wp) :: g
      if (cavity_l_plus_is_neutral(lp_new)) then
         g = huge(1.0_wp)
      else
         g = log(lp_new) - x
      end if
   end function cavity_outer_residual

   pure subroutine cavity_state_at_x(x, par, ice, eos, const, u_star, T_w, S_w, p_b, S_i, &
                                     T_b, S_b, m_mass, gamma_t, gamma_s, b_flux, lp_new, ierr)
      !! Evaluate the whole interface at a trial `x = ln(L+)`: exchange
      !! velocities at that stratification, the closed-form three-equation
      !! solve, the buoyancy flux the answer implies and the `L+` it
      !! re-diagnoses.  The fixed point of `x -> ln(L+_new)` is the
      !! solution of the implicit system.
      !$acc routine seq
      real(wp), intent(in) :: x
         !! Trial `ln(L+)`.  At or above `CAVITY_LP_X_HI` the trial `L+`
         !! is the neutral sentinel.
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Exchange-law bundle.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s).
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa).
      real(wp), intent(in) :: S_i
         !! Ice salinity (g/kg).
      real(wp), intent(out) :: T_b
         !! Interface temperature (degC).
      real(wp), intent(out) :: S_b
         !! Interface salinity (g/kg).
      real(wp), intent(out) :: m_mass
         !! Melt mass flux (kg/m^2/s).
      real(wp), intent(out) :: gamma_t
         !! Heat exchange velocity (m/s) at this iterate.
      real(wp), intent(out) :: gamma_s
         !! Salt exchange velocity (m/s) at this iterate.
      real(wp), intent(out) :: b_flux
         !! Interfacial buoyancy flux (m^2/s^3) the answer implies.
      real(wp), intent(out) :: lp_new
         !! `L+` re-diagnosed from that buoyancy flux.
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: l_plus

      b_flux = 0.0_wp
      lp_new = CAVITY_L_PLUS_NEUTRAL
      if (x >= CAVITY_LP_X_HI) then
         l_plus = CAVITY_L_PLUS_NEUTRAL
      else
         l_plus = exp(x)
      end if

      ! `T_b`/`S_b` are the RESERVED-law arguments here (only MK18 reads
      ! them); the trial interface state is not known yet, so the
      ! prototype's convention of passing `(0, S_w)` is kept.
      call cavity_exchange_velocities(par, const, u_star, l_plus, T_w, S_w, &
                                      0.0_wp, S_w, gamma_t, gamma_s, ierr)
      if (ierr /= CAVITY_MELT_OK) then
         call cavity_safe_state(eos, S_w, p_b, T_b, S_b, m_mass)
         return
      end if

      call cavity_three_equation(T_w, S_w, p_b, gamma_t, gamma_s, S_i, ice, eos, &
                                 const, T_b, S_b, m_mass, ierr)
      if (ierr /= CAVITY_MELT_OK) return

      b_flux = cavity_buoyancy_flux(const, T_w, S_w, T_b, S_b, gamma_t, gamma_s)
      lp_new = cavity_l_plus_from_state(const, u_star, b_flux)
   end subroutine cavity_state_at_x

   pure subroutine cavity_solve_melt(T_w, S_w, p_b, u_star, S_i, par, ice, eos, const, &
                                     sol, ierr)
      !! Solve the three-equation system with any implemented exchange
      !! law, and return the full interface state plus the fluxes a
      !! coupling seam will consume.
      !!
      !! For an EXPLICIT law (`const_gamma`) this is one call to the
      !! closed-form quadratic.
      !!
      !! For a STRATIFICATION-DEPENDENT law (`hj99`, `yung25`) the system
      !! is implicit — the exchange coefficients depend on the interfacial
      !! buoyancy flux, which depends on the melt rate they produce — and
      !! the outer iteration is BISECTION on `x = ln(L+)` over
      !! `[CAVITY_LP_X_LO, CAVITY_LP_X_HI]`.  It is guaranteed convergent
      !! because the residual `G(x) = ln(L+_new(x)) - x` provably changes
      !! sign across that bracket:
      !!
      !!   * `x -> x_lo` (maximal suppression): `gamma -> 0`, so
      !!     `|B_b| -> 0` and `L+_new -> +infinity`, giving `G > 0`;
      !!   * `x -> x_hi` (neutral): `gamma` sits at its cap, `|B_b|` is
      !!     maximal and `L+_new` is finite, giving `G < 0`.
      !!
      !! UNLESS the interface is DESTABILISING (`B_b >= 0`, i.e. freezing
      !! or a strongly cooling interface), in which case `L+_new` is
      !! `+infinity` everywhere and the neutral limit IS the fixed point —
      !! short-circuited, not iterated.  That is the branch KINK: H&J99
      !! sets `eta* := 1` there (p. 1792) and Yung et al. (2025) falls
      !! back to the ConstCoeff values (pp. 5832-5833).  The melt rate is
      !! continuous across it; its derivative is not, which is exactly why
      !! this is bisection and not Newton.
      !!
      !! After convergence the whole column is re-evaluated ONCE at the
      !! converged scalar, so the returned exchange velocities and the
      !! returned `(T_b, S_b, m_mass)` belong to ONE value of `L+`.
      !$acc routine seq
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).  Must exceed `S_i`.
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa).
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s), strictly positive — from
         !! `cavity_ustar`.
      real(wp), intent(in) :: S_i
         !! Ice salinity (g/kg), >= 0.
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Exchange-law bundle.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      type(ocean_cavity_solution_t), intent(out) :: sol
         !! Interface state, fluxes and solver diagnostics.  On any non-OK
         !! status this carries the safe state (`m_mass` exactly zero).
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      real(wp) :: T_b, S_b, m_mass, gamma_t, gamma_s, b_flux, lp_new, l_plus
      real(wp) :: lo, hi, mid, xstar, g_lo, g_hi, g_mid
      integer :: it
      logical :: converged, settled

      ierr = CAVITY_MELT_OK
      call cavity_solution_reset(sol)
      sol%u_star = u_star
      sol%l_plus = CAVITY_L_PLUS_NEUTRAL
      sol%L_obukhov = CAVITY_L_PLUS_NEUTRAL
      call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)

      if (.not. (ieee_is_finite(T_w) .and. ieee_is_finite(S_w) .and. &
                 ieee_is_finite(p_b) .and. ieee_is_finite(u_star) .and. &
                 ieee_is_finite(S_i))) then
         ierr = CAVITY_MELT_NONFINITE_INPUT
         return
      end if
      if (u_star <= 0.0_wp) then
         ierr = CAVITY_MELT_BAD_INPUT
         return
      end if

      converged = .true.
      it = 0
      l_plus = CAVITY_L_PLUS_NEUTRAL

      ! The neutral evaluation is both the answer for an explicit law and
      ! the upper bracket for an implicit one.
      call cavity_state_at_x(CAVITY_LP_X_HI, par, ice, eos, const, u_star, T_w, S_w, &
                             p_b, S_i, T_b, S_b, m_mass, gamma_t, gamma_s, b_flux, &
                             lp_new, ierr)
      if (ierr /= CAVITY_MELT_OK) then
         call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
         return
      end if
      l_plus = lp_new

      if (cavity_law_is_implicit(par%law)) then
         g_hi = cavity_outer_residual(lp_new, CAVITY_LP_X_HI)
         ! `cavity_l_plus_is_neutral(lp_new)` is the destabilising
         ! short-circuit and makes `g_hi = +huge >= 0`; the explicit
         ! `g_hi >= 0` test below therefore covers both accept branches.
         if (g_hi < 0.0_wp) then
            call cavity_state_at_x(CAVITY_LP_X_LO, par, ice, eos, const, u_star, T_w, &
                                   S_w, p_b, S_i, T_b, S_b, m_mass, gamma_t, gamma_s, &
                                   b_flux, lp_new, ierr)
            if (ierr /= CAVITY_MELT_OK) then
               call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
               return
            end if
            g_lo = cavity_outer_residual(lp_new, CAVITY_LP_X_LO)
            if (g_lo <= 0.0_wp) then
               ! The bracket argument failed — refuse rather than return
               ! whichever end happens to look plausible.
               call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
               ierr = CAVITY_MELT_NOT_CONVERGED
               return
            end if

            lo = CAVITY_LP_X_LO
            hi = CAVITY_LP_X_HI
            converged = .false.
            do it = 1, CAVITY_MAX_ITER
               mid = 0.5_wp*(lo + hi)
               settled = (mid <= lo) .or. (mid >= hi)
               if (.not. settled) then
                  call cavity_state_at_x(mid, par, ice, eos, const, u_star, T_w, S_w, &
                                         p_b, S_i, T_b, S_b, m_mass, gamma_t, gamma_s, &
                                         b_flux, lp_new, ierr)
                  if (ierr /= CAVITY_MELT_OK) then
                     call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
                     return
                  end if
                  g_mid = cavity_outer_residual(lp_new, mid)
                  if (g_mid == 0.0_wp) then
                     settled = .true.
                  else if (g_mid > 0.0_wp) then
                     lo = mid
                  else
                     hi = mid
                  end if
                  if ((.not. settled) .and. (hi - lo < CAVITY_LP_X_TOL)) settled = .true.
               end if
               if (settled) then
                  converged = .true.
                  exit
               end if
            end do
            if (.not. converged) then
               call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
               ierr = CAVITY_MELT_NOT_CONVERGED
               return
            end if

            ! One final evaluation at the converged scalar, so the
            ! returned gammas and the returned interface state belong to
            ! the SAME `L+`.
            xstar = 0.5_wp*(lo + hi)
            call cavity_state_at_x(xstar, par, ice, eos, const, u_star, T_w, S_w, p_b, &
                                   S_i, T_b, S_b, m_mass, gamma_t, gamma_s, b_flux, &
                                   lp_new, ierr)
            if (ierr /= CAVITY_MELT_OK) then
               call cavity_safe_state(eos, S_w, p_b, sol%T_b, sol%S_b, sol%m_mass)
               return
            end if
            l_plus = exp(xstar)
         end if
      end if

      sol%T_b = T_b
      sol%S_b = S_b
      sol%m_mass = m_mass
      sol%gamma_t = gamma_t
      sol%gamma_s = gamma_s
      sol%b_flux = b_flux
      sol%l_plus = l_plus
      sol%L_obukhov = cavity_obukhov_length(const, u_star, b_flux)
      sol%T_star = T_w - eos_freezing_point(eos, S_w, p_b)
      sol%S_star = S_w - S_b
      call cavity_heat_fluxes(T_w, T_b, m_mass, gamma_t, ice, const, &
                              sol%q_ocean, sol%q_ice, sol%q_latent)
      sol%n_iter = it
      sol%converged = converged
   end subroutine cavity_solve_melt

   pure subroutine cavity_melt_point(T_w, S_w, p_b, u_star, S_i, par, ice, eos, const, &
                                     T_b, S_b, m_mass, q_ocean, ierr)
      !! Scalar entry point returning only what a coupling seam consumes:
      !! the interface state, the canonical melt mass flux and the ocean
      !! -> interface heat flux.  A thin wrapper over `cavity_solve_melt`
      !! that keeps the solution BUNDLE inside the callee, so a
      !! `do concurrent` over columns needs no derived-type `local(...)`
      !! clause at all — each iteration writes its own array elements.
      !$acc routine seq
      real(wp), intent(in) :: T_w
         !! Far-field temperature (degC).
      real(wp), intent(in) :: S_w
         !! Far-field salinity (g/kg).
      real(wp), intent(in) :: p_b
         !! Interface pressure (Pa).
      real(wp), intent(in) :: u_star
         !! Friction velocity (m/s).
      real(wp), intent(in) :: S_i
         !! Ice salinity (g/kg).
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Exchange-law bundle.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle.
      real(wp), intent(out) :: T_b
         !! Interface temperature (degC).
      real(wp), intent(out) :: S_b
         !! Interface salinity (g/kg).
      real(wp), intent(out) :: m_mass
         !! Melt mass flux (kg/m^2/s), > 0 melting; EXACTLY zero on any
         !! non-OK status.
      real(wp), intent(out) :: q_ocean
         !! Turbulent heat flux ocean -> interface (W/m^2).
      integer, intent(out) :: ierr
         !! `CAVITY_MELT_*` status.
      type(ocean_cavity_solution_t) :: sol

      call cavity_solve_melt(T_w, S_w, p_b, u_star, S_i, par, ice, eos, const, sol, ierr)
      T_b = sol%T_b
      S_b = sol%S_b
      m_mass = sol%m_mass
      q_ocean = sol%q_ocean
   end subroutine cavity_melt_point

   pure subroutine cavity_melt_columns(n, T_w, S_w, p_b, u_star, S_i, par, ice, eos, &
                                       const, T_b, S_b, m_mass, q_ocean, ierr_col)
      !! Data-parallel driver: solve `n` independent columns.
      !!
      !! WHY THIS EXISTS, and it is not only convenience.  On the GPU
      !! build the whole kernel lives in `librdb_core.so`, and nvlink
      !! cannot resolve a `!$acc routine seq` device symbol from a SHARED
      !! library into a `do concurrent` compiled in a DIFFERENT
      !! translation unit — the test that first tried it failed with
      !! `nvlink error: Undefined reference to
      !! 'rdb_ocean_cavity_melt_cavity_solve_melt_'`, on the GPU
      !! toolchain only (gfortran linked and ran it happily).  So the
      !! column loop must live in the same object as the routine it
      !! calls: HERE, not in the caller.  The coupling PR's cavity kernel
      !! should therefore extend this routine (or sit beside it in this
      !! module), not write its own `do concurrent` over
      !! `cavity_solve_melt` in the engine.
      !!
      !! No locality clause is needed: every iteration writes only its
      !! own `i`-th elements, and `cavity_melt_point` keeps the solution
      !! bundle inside the callee.  Explicit-shape dummies throughout,
      !! per the repo's `do concurrent` descriptor rule.
      !!
      !! `mem:separate` contract: this routine moves NOTHING.  The caller
      !! owns the mapping — every one of the ten arrays must already be
      !! device-present (`!$acc enter data copyin(...)` for the five
      !! inputs, `create(...)` for the five outputs) and the results must
      !! be pulled back with `!$acc update self(...)`.
      !!
      !! HOST routine — it LAUNCHES the kernel, so it carries no
      !! `!$acc routine seq` of its own; `cavity_melt_point` and
      !! everything below it do.
      integer, intent(in) :: n
         !! Number of columns.
      real(wp), intent(in) :: T_w(n)
         !! Far-field temperature per column (degC).
      real(wp), intent(in) :: S_w(n)
         !! Far-field salinity per column (g/kg).
      real(wp), intent(in) :: p_b(n)
         !! Interface pressure per column (Pa).
      real(wp), intent(in) :: u_star(n)
         !! Friction velocity per column (m/s).
      real(wp), intent(in) :: S_i(n)
         !! Ice salinity per column (g/kg).
      type(ocean_cavity_exchange_t), intent(in) :: par
         !! Exchange-law bundle, shared by every column.
      type(ocean_cavity_ice_t), intent(in) :: ice
         !! Ice-conduction bundle, shared by every column.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.  Flat POD, by value.
      type(ocean_cavity_const_t), intent(in) :: const
         !! Constants bundle, shared by every column.
      real(wp), intent(out) :: T_b(n)
         !! Interface temperature per column (degC).
      real(wp), intent(out) :: S_b(n)
         !! Interface salinity per column (g/kg).
      real(wp), intent(out) :: m_mass(n)
         !! Melt mass flux per column (kg/m^2/s), > 0 melting.
      real(wp), intent(out) :: q_ocean(n)
         !! Ocean -> interface heat flux per column (W/m^2).
      integer, intent(out) :: ierr_col(n)
         !! `CAVITY_MELT_*` status PER COLUMN.  A column kernel must not
         !! take the run down for one bad column, so the failures are
         !! counted by the caller, not raised here.
      integer :: i

      do concurrent(i=1:n)
         call cavity_melt_point(T_w(i), S_w(i), p_b(i), u_star(i), S_i(i), par, ice, &
                                eos, const, T_b(i), S_b(i), m_mass(i), q_ocean(i), &
                                ierr_col(i))
      end do
   end subroutine cavity_melt_columns

end module rdb_ocean_cavity_melt
