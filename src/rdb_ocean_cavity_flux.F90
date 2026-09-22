!! Ice-shelf basal melt, COUPLED: the slot, the far-field sampler and the
!! once-per-thermo-step driver that turns `rdb_ocean_cavity_melt`'s scalar
!! three-equation kernel into two surface-flux components the ocean
!! integrates.
!!
!! `&ocean_cavity_melt_nml enable` (default `.false.`) is the gate; the
!! geometry it stands on (`z_draft`, `cover_frac`, `p_ice_ref`) is
!! `&ocean_cavity_dyn_nml`'s and is required.  Knob off ⇒ the slot's
!! arrays stay at their `(1,1)` placeholders, no kernel is launched, and
!! every path is byte-identical.
!!
!! ## What this module owns
!!
!! 1. `ocean_cavity_flux_t` — the 2-D device-resident state of the
!!    interface: the sampled far field, `u*`, the interface `(T_b, S_b)`,
!!    the canonical melt mass flux, the ocean → interface heat flux, and
!!    the per-column solver status.
!! 2. `cavity_far_field_impl` — the far-field sample, in METRES below the
!!    ice base, never "layer `nz`".
!! 3. `ocean_cavity_flux_step` — the driver: sample, solve (through
!!    `cavity_melt_columns_2d`, which must own the `do concurrent`; see
!!    its docstring for the nvlink reason), deliver, count, fail loud.
!!
!! ## Far-field sampling (why metres)
!!
!! The melt rate a three-equation law returns is roughly linear in the
!! thermal driving it is handed, and the thermal driving depends on HOW
!! FAR FROM THE ICE the model sampled it.  That is the dominant
!! resolution artefact in the subject, not a refinement: Gwyther et al.
!! (2020), Burchard et al. (2022) Table 2 p. 15 (the all-bulk melt-rate
!! error GROWS under refinement — -8 % at 6.66 m, -30 % at 0.015 m) and
!! Yung et al. (2026) p. 2074 all identify it.  Sampling "the top layer"
!! would therefore make the melt rate a function of the vertical
!! coordinate's cell thickness — which is precisely the quantity a
!! coordinate study must hold fixed.  So `far_field_depth` is in metres
!! and the sampler is a thickness-weighted mean over the layers that
!! span it, with a PARTIAL last layer.
!!
!! Vanished layers (`h <= H_VANISHED`, the D4 skip/merge marker) carry no
!! mass and are skipped rather than clamped.  A column whose sample finds
!! no mass at all is dropped from the cover mask — zero melt, `OK`
!! status — not solved on a fabricated state.
!!
!! ## The three pressures, and which one this is
!!
!! The liquidus is evaluated at `multilayer_state_t%p_top` — THE
!! interface pressure, the same field the FV_MOM6 pressure-stack
!! boundary condition and the in-situ EOS read (`src/core/ocean/README.md`,
!! the `p_top` seam contract).  It is NOT `eos%p_ref` (a scalar potential
!! -density reference, deliberately horizontally uniform) and NOT the
!! `eta_forcing` seam (a gradient).  One pressure, three consumers, and
!! this module is a CONSUMER — it writes nothing to `p_top`.  The sole
!! producer is `configure_ocean_cavity`, which assembles
!! `p_top = metrics%p_ice_ref + sf%p_surf` (re-assembled per outer step
!! in `ocean_dyn_step_split` when the psurf seam makes `sf%p_surf`
!! live); `configure_ocean_cavity_melt` asserts that it ran rather than
!! seeding a second copy, because a second writer running after it would
!! silently drop `sf%p_surf`.  A consequence worth stating: an
!! atmospheric load under the shelf depresses the freezing point exactly
!! as the ice load does, with no extra wiring.
!!
!! ## Sign and units — the table that decides whether the answer is right
!!
!! | quantity | symbol | unit | sign |
!! |---|---|---|---|
!! | melt | `melt` (`m_mass`) | kg/m^2/s | **> 0 melting** |
!! | kernel heat flux | `q_ocean` | W/m^2 | **> 0 warms the INTERFACE** (the ocean cools) |
!! | delivered heat | `sf%heat_cavity` | W/m^2 | positive DOWN into the ocean |
!! | delivered salt | `sf%salt_cavity` | same as `sf%salt_flux` | positive SALINIFIES |
!!
!! **Heat.**  `q_ocean = rho_w*c_w*gamma_t*(T_w - T_b) > 0` is the
!! turbulent flux the interface takes FROM the ocean, so the component
!! the ocean must be given is its negative:
!!
!!     heat_cavity = -q_ocean
!!
!! The assembler adds it to `Q_heat` with a plain `+`, and
!! `apply_surface_src_2d_impl` deposits `dt/(rho0*cp) * Q_heat` at
!! `k = nz` — positive warms.  Warm water under a shelf therefore COOLS
!! the top of the column, which is the whole point.
!!
!! **Salt — and why it is not the same construction.**  Melting adds
!! freshwater MASS the Boussinesq column does not yet carry (that is
!! Phase 3).  Emulate the dilution with a salt flux and the fixed-mass
!! equivalent is exact, not approximate.  For a column of mass `M` per
!! unit area receiving mass `m` at salinity `S_i`:
!!
!!     d(M*S)/dt = m*S_i ,   dM/dt = m
!!  => M*dS/dt = m*S_i - S*m = -m*(S - S_i)
!!
!! so the virtual flux that reproduces `dS/dt` at FIXED `M` is
!!
!!     salt_cavity = -m_mass*(S_far - s_ice)
!!
!! referenced to the SAME far-field salinity the solve was handed.
!! Melting (`m_mass > 0`, `S_far > s_ice`) freshens.
!!
!! The heat twin of that `-S*dM/dt` dilution term is `-m*c_w*(T_w - T_b)`
!! and it is DELIBERATELY DROPPED.  The asymmetry is physical, not
!! sloppiness: `S_far - s_ice ~ 34 g/kg` is O(1) of the salinity itself,
!! while `T_w - T_b` is a few hundredths of a degree, so the dropped heat
!! term is `m/(rho_w*gamma_t) ~ 1e-3` of `q_ocean` — a 0.1 % correction
!! that belongs with the real-mass PR, whereas the salt term is the
!! entire meltwater buoyancy signal and cannot wait for it.
!!
!! NOTE the prototype's `column_demo.py` removes `f_turb = rho_w*gamma_s*
!! (S_w - S_b) = m_mass*(S_b - s_ice)` instead.  That is the salt flux
!! ACROSS the interface, not the dilution equivalent — it omits the
!! meltwater that in reality returns to the column carrying `S_b`, and so
!! under-freshens by `(S_b - s_ice)/(S_far - s_ice)`.  That demo is
!! explicitly not part of the oracle (its README: "for reading"); the
!! identity above is what the ocean's fixed-mass salinity equation
!! requires, and it is what the budget test asserts.
!!
!! ## Phase 3 — REAL freshwater MASS (`&ocean_cavity_melt_nml freshwater`)
!!
!! `freshwater = "virtual"` is the default and everything above is what
!! runs.  `freshwater = "mass"` moves the meltwater as a real Boussinesq
!! VOLUME instead.  The derivation, in full, because every sign and
!! every omitted term below is load-bearing.
!!
!! ### Which density converts kg to m
!!
!! `m_mass` is a MASS flux (kg m^-2 s^-1).  The model is BOUSSINESQ: it
!! conserves VOLUME, `h_layer` is a thickness, and every extensive total
!! the console prints is a volume integral scaled by a constant density.
!! The one density that is allowed to convert a mass flux into a
!! thickness tendency is therefore the Boussinesq reference `rho_0` —
!! the SAME `rho_0` the surface-flux apply already divides by
!! (`dS/dt = Q_salt/(rho_0 h)`), and the same one the EOS carries.
!! Using the meltwater's own freshwater density (~1000) instead would
!! make a kilogram of meltwater occupy more volume than a kilogram of
!! the sea water it joins, which in a Boussinesq model is not a
!! refinement but an inconsistency: the column's mass total is
!! `rho_0 * volume` by construction, so a second density would put the
!! tracked source and the tracked total on different scales.  So
!!
!!     w  ==  m_mass / rho_0          [m/s, > 0 melting]
!!     dh ==  w * dt                  added to `h_layer(:,:,nz)`
!!
!! (The console's `total_mass = sum(h*areaT)*RHO_WATER` uses the MODULE
!! constant `RHO_WATER`, which a run may configure `rho_0` away from.
!! The tracked mass source is therefore `RHO_WATER * sum(dh*areaT)`, i.e.
!! `(RHO_WATER/rho_0) * sum(m*area*dt)` — exactly `sum(m*area*dt)` at the
!! default `rho_0 = RHO_WATER`.  The accumulator must measure the same
!! mass the total does; it is not free to pick its own density.)
!!
!! ### The three top-layer tendencies
!!
!! Take the ocean column as the control volume.  Its top boundary is the
!! ice base, and the only things that cross it are the meltwater (mass
!! `m`, salinity `s_ice`, temperature `T_b`) and the turbulent fluxes.
!!
!! **Mass.**  `d h/dt = w`.
!!
!! **Salt.**  Across the ocean's top boundary, salt enters ADVECTIVELY
!! with the water leaving the interface (`m*S_b`) and leaves
!! DIFFUSIVELY into the interface (`rho_w*gamma_S*(S_w - S_b)`).  The
!! interface itself stores nothing, so its own balance — the third of
!! the three equations — is
!!
!!     m*s_ice + rho_w*gamma_S*(S_w - S_b) = m*S_b
!!
!! and the NET salt the ocean gains is
!!
!!     m*S_b - rho_w*gamma_S*(S_w - S_b) = m*S_b - m*(S_b - s_ice)
!!                                       = m*s_ice
!!
!! — the ice is the only thing on the other side of the boundary, and it
!! carries `s_ice` per kilogram.  The advective and diffusive halves
!! cancel EXACTLY; there is nothing to approximate.  So
!!
!!     d(h*S)/dt = w*s_ice
!!
!! and for the ISOMIP+ `s_ice = 0` the salt content of the column does
!! not change at all — the salinity falls purely because `h` grows.
!!
!! **Heat.**  The same control volume loses the turbulent flux
!! `q_ocean` and gains the meltwater's enthalpy.  The model's heat
!! content is `rho_0*c_p*integral(T dz)` with an implicit 0 degC
!! reference (`total_heat = sum(h*T*areaT)*RHO_WATER` carries no
!! offset), so in `hTr = h*T` space the enthalpy of the added water is
!! simply `dh*T_b`: NO heat capacity enters, and the model's `c_p` and
!! the melt law's own `c_w` never have to agree.
!!
!!     d(h*T)/dt = -q_ocean/(rho_0*c_p) + w*T_b
!!
!! ### Why the virtual salt flux must NOT also be applied — and why it
!! ### is nevertheless still assembled
!!
!! Applying both would double-count: the dilution would happen once
!! through the growing `h` and once through the removed salt.  But the
!! virtual flux is ALSO the entire meltwater buoyancy signal that KPP
!! and EPBL read out of `Q_salt` to build `B_0` (`beta*Q_salt/rho_0` is
!! ~6x the `alpha*Q_heat/(rho_0*c_p)` term at these temperatures), and
!! the real volume source does not reach `B_0` at all.  Dropping it from
!! `Q_salt` would leave both boundary-layer schemes blind to the
!! dominant, destabilising half of the surface buoyancy flux.
!!
!! So under `"mass"` the assembled `Q_salt` is UNCHANGED — `salt_cavity`
!! keeps carrying `-m*(S_far - s_ice)` and `B_0` is the same under both
!! forms, which is the point — and `ocean_cavity_mass_step` removes that
!! increment again from the TRACER in the same stage, adding the real
!! `w*s_ice` in its place.  The undo is the exact negation of what
!! `apply_surface_src_2d_impl` stamped (same `dt/rho_0`, same field, and
!! `cover_frac = 1` makes the atmospheric group exactly zero on a
!! covered column), so `x + (-x)` is zero in floating point and the
!! salinity tracer load is left EXACTLY unchanged when `s_ice = 0`.
!! Both increments are mirrored into `salt_budget_surface` in the same
!! order, so the budget sees the identical pair and closes to round-off.
!! The pseudo-salt mirror receives the same correction, so it keeps
!! receiving "exactly salinity's surface salt flux" — its documented
!! purpose — under either form.
!!
!! Passive tracers need NOTHING: `h` grows, `h*C` does not, so a tracer
!! with zero concentration in the meltwater is diluted automatically.
!! The ideal-age tracer is diluted the same way, which is the right
!! answer (meltwater is new water).
!!
!! ### Agreement with the virtual form, and the size of the difference
!!
!! With `eps = w*dt/h0` the two forms differ EXACTLY by
!!
!!     dS_mass - dS_virtual = +(S0 - s_ice) * eps^2/(1 + eps)
!!     dT_mass - dT_virtual = -eps*(A + T0 - T_b)/(1 + eps),
!!        A = -q_ocean*dt/(rho_0*c_p*h0) = dT_virtual
!!
!! so the salinity difference is second order in `eps` and the
!! temperature difference is `eps` times the (tiny) dropped dilution
!! term plus `eps` times the virtual increment itself.  `test_ocean_
!! cavity_freshwater` asserts those two identities rather than asserting
!! "the two agree", because agreeing to zero is what a double count
!! would also do.
!!
!! ### Sea level (`volume_compensation`)
!!
!! Real mass into a CLOSED domain raises it.  ISOMIP+ Sect. 3.1.3 states
!! the protocol's position: the freshwater flux "is not compensated by a
!! corresponding removal of water elsewhere in the domain" for the
!! restored Ocean0-2 configurations — whose only open seam is the
!! northern restoring band, which relaxes T and S and moves no volume —
!! whereas the closed Ocean3/4 boxes may compensate.  Scale, for Ocean0:
!! ~30 m/yr of melt over ~1e10 m^2 of shelf into ~4e10 m^2 of surface is
!! metres per year of sea level.  Hence the knob.
!!
!! `volume_compensation = "uniform_open_ocean"` removes the
!! domain-integrated melt volume again each thermo step, spread
!! uniformly PER UNIT AREA over the wet cells the ice does not cover,
!! each parcel carrying that cell's own `T` and `S` so no concentration
!! there changes.  It is tracked as a sink in all three budgets (a
!! parcel that carries `S` does change the domain SALT total, so the
!! salt and heat budgets must name it too).
!!
!! ### Budgets
!!
!! Mass gains a tracked SOURCE (`ms%mass_src`), accumulated with the same
!! per-stage weight `ms%mass_out` uses, so the console residual
!! `(M - M0) + mass_out - mass_src` stays at round-off while the total
!! legitimately grows.  Salt and heat need no new accumulator: every
!! increment this path makes is mirrored into the EXISTING
!! `salt_budget_surface` / `heat_budget_surface` contributors, by the
!! same in-stage kernel and at the same cadence as the surface fluxes,
!! so they carry the right `ocean_budget_stage_weight` by construction.
!!
!! ## Budgets
!!
!! Both components ride `Q_heat`/`Q_salt`, so they are integrated by
!! `ocean_surface_flux_apply_tracers` and land in the EXISTING
!! `ms%heat_budget_surface` / `ms%salt_budget_surface` contributors —
!! which the console already folds in with the correct
!! `ocean_budget_stage_weight`.  No separate frazil-style accumulator is
!! needed, and no full-weight/half-weight decision is taken here: the
!! source is applied by the same per-stage kernel as every other surface
!! flux, so it carries the same weight by construction.
!!
!! ## `mem:separate`
!!
!! Every array this module touches is device-resident and mapped by the
!! slot's `enter_data`.  Nothing is copied host↔device per step except
!! the four status COUNTS, which come back as `do concurrent ... reduce`
!! results (host scalars), exactly like `continuity_t%n_limited_step`.
!!
!! ## Citations (papers, never another model's source)
!!
!!   * Holland, D. M. and Jenkins, A. (1999): J. Phys. Oceanogr. 29,
!!     1787-1800.
!!   * Jenkins, A., Nicholls, K. W. and Corr, H. F. J. (2010): J. Phys.
!!     Oceanogr. 40, 2298-2312.
!!   * Asay-Davis, X. S. et al. (2016): Geosci. Model Dev. 9, 2471-2497
!!     (ISOMIP+).
!!   * Burchard, H. et al. (2022): Ocean Modelling 179, 102119.
!!   * Yung, C. K. et al. (2025): The Cryosphere 19, 5827-5861.
!!   * Yung, C. K. et al. (2026): The Cryosphere 20, 2053-2088.
module rdb_ocean_cavity_flux
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS, RHO_WATER
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t
   use rdb_mem_report, only: arr_bytes
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_cavity_melt, only: ocean_cavity_const_t, ocean_cavity_exchange_t, &
                                    ocean_cavity_ice_t, cavity_melt_columns_2d, &
                                    CAVITY_MELT_OK, CAVITY_MELT_NONFINITE_INPUT, &
                                    CAVITY_MELT_NONFINITE_STATE, &
                                    CAVITY_MELT_NOT_CONVERGED, &
                                    CAVITY_MELT_NO_PHYSICAL_ROOT, &
                                    CAVITY_FW_VIRTUAL, CAVITY_FW_MASS, &
                                    CAVITY_VC_NONE, CAVITY_VC_UNIFORM_OPEN
   implicit none
   private

   public :: ocean_cavity_flux_t
   public :: ocean_cavity_flux_step
   public :: ocean_cavity_mass_step
   public :: cavity_far_field_impl
   public :: cavity_flux_fill_impl
   public :: cavity_status_counts_impl
   public :: cavity_melt_status_is_fatal
   public :: cavity_mass_apply_impl
   public :: cavity_mass_salt_mirror_impl
   public :: cavity_mass_totals_impl
   public :: cavity_comp_apply_impl
   public :: cavity_comp_scale_tracer_impl
   public :: cavity_mass_thin_is_fatal
   public :: cavity_comp_withdrawal
   public :: H_CAVITY_FLOOR

   real(wp), parameter :: H_CAVITY_FLOOR = 2.0_wp*H_VANISHED
      !! Thickness (m) a clamped top-layer WITHDRAWAL is allowed to leave
      !! behind — one full marker ABOVE `H_VANISHED`, not on it.
      !!
      !! `H_VANISHED` is the D4 skip/merge MARKER, and every gate in the
      !! tree tests it with a strict `>`: a layer sitting exactly ON the
      !! marker reads as VANISHED.  So clamping a melt or compensation
      !! withdrawal to `H_VANISHED` exactly — which is what these kernels
      !! used to do, with an in-line comment arguing for landing on the
      !! marker — hands the next ALE remap a layer it treats as empty:
      !! `ocean_remap_tracer_column`'s `c = hTr/h` guard returns 0 there,
      !! so the layer's heat and salt content is DELETED, silently and at
      !! the next thermo step.  One ulp above the marker it is preserved.
      !! A clamped withdrawal is already a fail-loud condition
      !! (`cavity_mass_thin_is_fatal`), but the abort happens AFTER the
      !! state is written, so the value written has to survive being read.
      !!
      !! `2*H_VANISHED` is the convention the rest of the tree already
      !! uses where a producer must land ABOVE the marker on purpose —
      !! `compute_target_h_rho_impl`'s `h_floor_eff = max(zstar_h_min,
      !! 2*H_VANISHED)` and `seed_land_h_floor`'s `max(angstrom_h,
      !! 2*H_VANISHED)`.  Written as a multiple of the constant of record
      !! rather than as a literal, so it cannot drift from it.

   type :: ocean_cavity_flux_t
      !! Ice-shelf basal-melt slot: the 2-D interface state plus the
      !! parameter bundles the kernel is called with.
      !!
      !! Every array is `(nx_total, ny_total)` when `enable`, `(1,1)`
      !! otherwise — the `z_draft` gating convention, latched in
      !! `ocean_state_init_from_config` BEFORE `init` so the allocation
      !! gate can read it.
      logical :: is_init = .false.
         !! Set by `init` after every allocation succeeds; cleared by
         !! `destroy` first.  Always test this, never `allocated(...)`.
      logical :: enable = .false.
         !! `&ocean_cavity_melt_nml enable`, latched before `init`.

      ! ---- Knobs the kernel is called with (host scalars) ----
      real(wp) :: far_field_depth = 10.0_wp
         !! Thickness (m) below the ice base the far field is averaged
         !! over.  See the module docstring: metres, not layers.
      real(wp) :: cdrag_top = 2.5e-3_wp
         !! Top drag coefficient for the MELT friction velocity only.
      real(wp) :: u_tide = 1.0e-2_wp
         !! RMS tidal velocity (m/s) in the melt `u*`.
      real(wp) :: ustar_min = 1.0e-4_wp
         !! Friction-velocity floor (m/s).
      real(wp) :: s_ice = 0.0_wp
         !! Ice salinity (g/kg).
      integer :: freshwater = CAVITY_FW_VIRTUAL
         !! `&ocean_cavity_melt_nml freshwater`, parsed.
         !! `CAVITY_FW_VIRTUAL` (default) ⇒ `ocean_cavity_mass_step` is
         !! an immediate return and the path is byte-identical.
      integer :: volume_comp = CAVITY_VC_NONE
         !! `&ocean_cavity_melt_nml volume_compensation`, parsed.
         !! Requires `freshwater = CAVITY_FW_MASS` (validate_config).
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3) — THE density that
         !! converts the melt MASS flux into a thickness tendency, and
         !! the same one `ocean_surface_flux_t%rho0` divides the surface
         !! tracer sources by.  Assigned at configure from the single
         !! `rho_0` of record; the literal is only the pre-configure type
         !! default.  The two MUST agree or the salt undo is not the
         !! exact negation of what the apply kernel stamped, which is
         !! asserted at configure.
      type(ocean_cavity_exchange_t) :: par
         !! Exchange-law selector + coefficients.  Flat POD; its
         !! `f_cor` member is unused here (the driver substitutes the
         !! per-column `f_cor` array).
      type(ocean_cavity_ice_t) :: ice
         !! Ice-conduction selector + `T_ice`.  Flat POD.
      type(ocean_cavity_const_t) :: const
         !! Thermodynamic + turbulence constants.  Flat POD, left at the
         !! ISOMIP+ protocol values (Asay-Davis et al. (2016) Table 4):
         !! `rho_w`, `c_w`, `alpha_T`, `beta_S` here are the MELT LAW's
         !! own calibrated constants, NOT copies of the model's
         !! Boussinesq `rho_0` — the same "out of scope on purpose"
         !! category as `RHO_WATER` and `&ocean_ice_nml rho_ocean` in
         !! `src/core/ocean/README.md`'s reference-density table.
         !! Overriding them would break ISOMIP+ comparability, which is
         !! the reason this path exists.

      ! ---- Static, filled once at configure ----
      real(wp), allocatable :: f_cor(:, :)
         !! Cell-centred Coriolis parameter (1/s), from the same
         !! `fill_coriolis_centre` every other slot uses.  Read by
         !! `hj99` only, as `|f|`; the law does not exist at `f = 0`.

      ! ---- Per-step interface state (all device-resident) ----
      real(wp), allocatable :: active(:, :)
         !! Composed solve mask (0/1): `cover_frac` AND wet AND "the
         !! far-field sample found mass".  This is what the kernel's
         !! `cover` argument is given.
      real(wp), allocatable :: t_far(:, :)
         !! Far-field temperature (degC), thickness-weighted over
         !! `far_field_depth`.
      real(wp), allocatable :: s_far(:, :)
         !! Far-field salinity (g/kg), same average.
      real(wp), allocatable :: u_far(:, :)
         !! Far-field x velocity at the cell CENTRE (m/s), same average.
      real(wp), allocatable :: v_far(:, :)
         !! Far-field y velocity at the cell CENTRE (m/s), same average.
      real(wp), allocatable :: ustar(:, :)
         !! Melt friction velocity (m/s).
      real(wp), allocatable :: t_b(:, :)
         !! Interface temperature (degC), on the liquidus.
      real(wp), allocatable :: s_b(:, :)
         !! Interface salinity (g/kg).
      real(wp), allocatable :: melt(:, :)
         !! **Canonical** melt mass flux (kg/m^2/s), > 0 melting.
      real(wp), allocatable :: q_ocean(:, :)
         !! Turbulent heat flux ocean → interface (W/m^2), > 0 cools the
         !! ocean.
      real(wp), allocatable :: gamma_t(:, :)
         !! Thermal exchange velocity (m/s) the column's solve converged
         !! on — NOT `par%gamma_t_coeff`, which is the dimensionless
         !! coefficient.  Under `hj99` / `yung25` it is an implicit
         !! function of the converged interface state, so it is stored
         !! rather than re-derived: a diagnostic that rebuilt it from
         !! `u*` alone would report the NEUTRAL value instead of the
         !! stratification-suppressed one.  Zero where inactive.
      real(wp), allocatable :: gamma_s(:, :)
         !! Haline exchange velocity (m/s), same convention.
      integer, allocatable :: status(:, :)
         !! `CAVITY_MELT_*` per column; `CAVITY_MELT_OK` where inactive.
      real(wp), allocatable :: comp_scale(:, :)
         !! Per-cell top-layer thickness RATIO `h_new/h_old` left behind
         !! by the `volume_compensation` sink, exactly `1` everywhere the
         !! sink did not act.  Every PASSIVE tracer's top-layer load is
         !! multiplied by it (`cavity_comp_scale_tracer_impl`) so the
         !! removed parcel carries that cell's own concentration and
         !! changes none of them.  Allocated with the rest of the slot
         !! (full size iff `enable`); left at `1` when
         !! `volume_comp = CAVITY_VC_NONE`.

      ! ---- Status counters (host scalars, reduced off the device) ----
      integer :: n_nonfinite_step = 0
         !! Active columns returning `NONFINITE_INPUT`/`NONFINITE_STATE`
         !! this step.  **FATAL** — the driver fails loud.
      integer :: n_not_converged_step = 0
         !! Active columns returning `NOT_CONVERGED` this step (zero melt
         !! applied there — the kernel's documented safe state).
      integer :: n_no_root_step = 0
         !! Active columns returning `NO_PHYSICAL_ROOT` this step.
      integer :: n_other_step = 0
         !! Active columns returning any other non-OK status
         !! (`BAD_INPUT`, `NO_CORIOLIS`, `LAW_DOMAIN`, ...).
      integer(int64) :: n_not_converged_total = 0_int64
         !! Running total, drained to the console like
         !! `continuity_t%n_limited_total`.  int64 because a
         !! long run with one stubborn column can exceed 2e9.
      integer(int64) :: n_no_root_total = 0_int64
         !! Running total of `NO_PHYSICAL_ROOT`.
      integer(int64) :: n_other_total = 0_int64
         !! Running total of every other non-OK, non-fatal status.
      integer :: n_thin_step = 0
         !! Columns this thermo step whose top layer could NOT give up
         !! the thickness the real-mass path asked of it without falling
         !! through `H_CAVITY_FLOOR` — a FREEZING column (`m < 0`) whose
         !! top layer is thinner than `|m|*dt/rho_0`, or a compensation
         !! sink deeper than the open-ocean top layer.  The withdrawal is
         !! CLAMPED to what is there ABOVE the floor (so `h` stays
         !! strictly above the vanish marker, never on or below it) and
         !! the count is FATAL: a clamped withdrawal is a
         !! withdrawal the tracked mass source no longer matches, so
         !! continuing would print a budget that silently stopped
         !! closing.  Same stance as the non-finite melt column.
      integer(int64) :: n_thin_total = 0_int64
         !! Running total of `n_thin_step`, for the record.
      real(wp) :: melt_volume_step = 0.0_wp
         !! Domain-integrated meltwater VOLUME (m^3) added this thermo
         !! step over the INTERIOR (ghosts excluded, globally combined).
         !! Diagnostic + the input the compensation sink spends.
      real(wp) :: open_area = 0.0_wp
         !! Interior wet area (m^2) the ice does NOT cover, globally
         !! combined — the denominator of the uniform sink.
      real(wp) :: comp_withdrawal_step = 0.0_wp
         !! Thickness (m) the compensation sink removed from each open
         !! -ocean column this thermo step; zero when compensation is off.
   contains
      procedure, non_overridable :: init => ocean_cavity_flux_init
      procedure, non_overridable :: destroy => ocean_cavity_flux_destroy
      procedure, non_overridable :: enter_data => ocean_cavity_flux_enter_data
      procedure, non_overridable :: exit_data => ocean_cavity_flux_exit_data
      procedure, non_overridable :: bytes => ocean_cavity_flux_bytes
   end type ocean_cavity_flux_t

contains

   ! ======================================================================
   ! Lifecycle
   ! ======================================================================

   subroutine ocean_cavity_flux_init(this, grid)
      !! Allocate the slot.  Gated on `enable` (latched by
      !! `ocean_state_init_from_config` before this runs), so a run
      !! without a cavity pays fourteen `(1,1)` placeholders.
      class(ocean_cavity_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny

      if (this%enable) then
         nx = grid%nx_total
         ny = grid%ny_total
      else
         nx = 1
         ny = 1
      end if

      allocate (this%f_cor(nx, ny), source=0.0_wp)
      allocate (this%active(nx, ny), source=0.0_wp)
      allocate (this%t_far(nx, ny), source=0.0_wp)
      allocate (this%s_far(nx, ny), source=0.0_wp)
      allocate (this%u_far(nx, ny), source=0.0_wp)
      allocate (this%v_far(nx, ny), source=0.0_wp)
      allocate (this%ustar(nx, ny), source=0.0_wp)
      allocate (this%t_b(nx, ny), source=0.0_wp)
      allocate (this%s_b(nx, ny), source=0.0_wp)
      allocate (this%melt(nx, ny), source=0.0_wp)
      allocate (this%q_ocean(nx, ny), source=0.0_wp)
      allocate (this%gamma_t(nx, ny), source=0.0_wp)
      allocate (this%gamma_s(nx, ny), source=0.0_wp)
      allocate (this%status(nx, ny), source=CAVITY_MELT_OK)
      ! `comp_scale` is a RATIO, so its inert value is 1, not 0 — a
      ! zero-initialised plane would annihilate every passive tracer's
      ! top layer the first time it was applied.
      allocate (this%comp_scale(nx, ny), source=1.0_wp)
      this%is_init = .true.
   end subroutine ocean_cavity_flux_init

   subroutine ocean_cavity_flux_destroy(this)
      !! Release the slot.  `is_init` is cleared FIRST.
      class(ocean_cavity_flux_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%f_cor)) deallocate (this%f_cor)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%t_far)) deallocate (this%t_far)
      if (allocated(this%s_far)) deallocate (this%s_far)
      if (allocated(this%u_far)) deallocate (this%u_far)
      if (allocated(this%v_far)) deallocate (this%v_far)
      if (allocated(this%ustar)) deallocate (this%ustar)
      if (allocated(this%t_b)) deallocate (this%t_b)
      if (allocated(this%s_b)) deallocate (this%s_b)
      if (allocated(this%melt)) deallocate (this%melt)
      if (allocated(this%q_ocean)) deallocate (this%q_ocean)
      if (allocated(this%gamma_t)) deallocate (this%gamma_t)
      if (allocated(this%gamma_s)) deallocate (this%gamma_s)
      if (allocated(this%status)) deallocate (this%status)
      if (allocated(this%comp_scale)) deallocate (this%comp_scale)
   end subroutine ocean_cavity_flux_destroy

   subroutine ocean_cavity_flux_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so
      !! the device-attach map base is the heap object, not a
      !! polymorphic stack box.
      class(ocean_cavity_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_cavity_flux_t)
         call ocean_cavity_flux_enter_data_impl(this)
      end select
   end subroutine ocean_cavity_flux_enter_data

   subroutine ocean_cavity_flux_enter_data_impl(this)
      !! `copyin` (not `create`) throughout: `f_cor` carries a
      !! configure-time host fill that MUST reach the device, and the
      !! rest carry the zero `init` promised on both toolchains.
      type(ocean_cavity_flux_t), intent(inout) :: this
      !$acc enter data copyin(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&                  this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&                  this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&                  this%gamma_s, this%status, this%comp_scale)
      !$acc update device(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&              this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&              this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&              this%gamma_s, this%status, this%comp_scale)
   end subroutine ocean_cavity_flux_enter_data_impl

   subroutine ocean_cavity_flux_exit_data(this)
      class(ocean_cavity_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_cavity_flux_t)
         call ocean_cavity_flux_exit_data_impl(this)
      end select
   end subroutine ocean_cavity_flux_exit_data

   subroutine ocean_cavity_flux_exit_data_impl(this)
      type(ocean_cavity_flux_t), intent(inout) :: this
      !$acc exit data delete(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&                 this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&                 this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&                 this%gamma_s, this%status, this%comp_scale)
   end subroutine ocean_cavity_flux_exit_data_impl

   pure function ocean_cavity_flux_bytes(this) result(nbytes)
      !! Counted allocatable footprint (0 when unallocated).  One
      !! `arr_bytes` term per array — add one here when an array joins
      !! the type.
      class(ocean_cavity_flux_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_cor) &
               + arr_bytes(this%active) &
               + arr_bytes(this%t_far) &
               + arr_bytes(this%s_far) &
               + arr_bytes(this%u_far) &
               + arr_bytes(this%v_far) &
               + arr_bytes(this%ustar) &
               + arr_bytes(this%t_b) &
               + arr_bytes(this%s_b) &
               + arr_bytes(this%melt) &
               + arr_bytes(this%q_ocean) &
               + arr_bytes(this%gamma_t) &
               + arr_bytes(this%gamma_s) &
               + arr_bytes(this%status) &
               + arr_bytes(this%comp_scale)
   end function ocean_cavity_flux_bytes

   ! ======================================================================
   ! Far-field sampling
   ! ======================================================================

   pure subroutine cavity_far_field_impl(nx, ny, nz, far_depth, cover, wet_mask, &
                                         h_layer, hTr_T, hTr_S, u_face_x, v_face_y, &
                                         active, t_far, s_far, u_far, v_far)
      !! Thickness-weighted mean of `(T, S, u, v)` over `far_depth`
      !! METRES below the ice base, with a PARTIAL last layer.
      !!
      !! Bottom-up stack (`k = nz` is the surface layer, i.e. the one
      !! against the ice base), so the walk runs `k = nz, 1, -1` and
      !! stops when the budget is spent.  `far_depth >= column
      !! thickness` ⇒ the whole column, which is the correct limit and
      !! not an error.
      !!
      !! Vanished layers (`h <= H_VANISHED`) are SKIPPED, not clamped —
      !! the D4 taxonomy's skip/merge marker.  They carry no mass, so
      !! including them would be dividing a zero tracer load by a
      !! near-zero thickness.
      !!
      !! Velocities are centred from the C-grid faces per layer BEFORE
      !! the vertical average (`u_c = 1/2 (u_{i} + u_{i+1})`), not after:
      !! the two commute only for a uniform column, and the order that
      !! matches "the speed the ice base feels" is the one that averages
      !! the centred, per-layer velocity.
      !!
      !! `active` is the composed solve mask handed to the melt kernel:
      !! covered AND wet AND the sample found mass.  A covered column
      !! whose entire sample is vanished gets `active = 0` and zeroed
      !! outputs — no melt, no status noise.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count; `k = nz` is the surface / ice-base layer.
      real(wp), intent(in) :: far_depth
         !! Sampling thickness (m), > 0.
      real(wp), intent(in) :: cover(nx, ny)
         !! Ice-cover fraction (v1 binary).
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Static wet (1) / land (0) mask.
      real(wp), intent(in) :: h_layer(nx, ny, nz)
         !! Layer thickness (m).
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
         !! `h*T` (degC m).
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
         !! `h*S` ((g/kg) m).
      real(wp), intent(in) :: u_face_x(nx + 1, ny, nz)
         !! Zonal face velocity (m/s).
      real(wp), intent(in) :: v_face_y(nx, ny + 1, nz)
         !! Meridional face velocity (m/s).
      real(wp), intent(out) :: active(nx, ny)
         !! Composed solve mask, 0 or 1.
      real(wp), intent(out) :: t_far(nx, ny)
         !! Sampled temperature (degC); 0 where inactive.
      real(wp), intent(out) :: s_far(nx, ny)
         !! Sampled salinity (g/kg); 0 where inactive.
      real(wp), intent(out) :: u_far(nx, ny)
         !! Sampled centred x velocity (m/s); 0 where inactive.
      real(wp), intent(out) :: v_far(nx, ny)
         !! Sampled centred y velocity (m/s); 0 where inactive.
      integer :: i, j, k
      real(wp) :: acc_h, acc_t, acc_s, acc_u, acc_v, remain, h, w, inv_h

      do concurrent(j=1:ny, i=1:nx) &
         local(k, acc_h, acc_t, acc_s, acc_u, acc_v, remain, h, w, inv_h)
         acc_h = 0.0_wp
         acc_t = 0.0_wp
         acc_s = 0.0_wp
         acc_u = 0.0_wp
         acc_v = 0.0_wp
         remain = far_depth
         if (cover(i, j) > 0.5_wp .and. wet_mask(i, j) > 0.5_wp) then
            do k = nz, 1, -1
               if (remain <= 0.0_wp) exit
               h = h_layer(i, j, k)
               if (h <= H_VANISHED) cycle
               w = min(h, remain)
               inv_h = 1.0_wp/h
               acc_h = acc_h + w
               acc_t = acc_t + w*hTr_T(i, j, k)*inv_h
               acc_s = acc_s + w*hTr_S(i, j, k)*inv_h
               acc_u = acc_u + w*0.5_wp*(u_face_x(i, j, k) + u_face_x(i + 1, j, k))
               acc_v = acc_v + w*0.5_wp*(v_face_y(i, j, k) + v_face_y(i, j + 1, k))
               remain = remain - w
            end do
         end if
         if (acc_h > 0.0_wp) then
            active(i, j) = 1.0_wp
            t_far(i, j) = acc_t/acc_h
            s_far(i, j) = acc_s/acc_h
            u_far(i, j) = acc_u/acc_h
            v_far(i, j) = acc_v/acc_h
         else
            active(i, j) = 0.0_wp
            t_far(i, j) = 0.0_wp
            s_far(i, j) = 0.0_wp
            u_far(i, j) = 0.0_wp
            v_far(i, j) = 0.0_wp
         end if
      end do
   end subroutine cavity_far_field_impl

   ! ======================================================================
   ! Flux delivery
   ! ======================================================================

   pure subroutine cavity_flux_fill_impl(nx, ny, s_ice, active, melt, q_ocean, s_far, &
                                         heat_cavity, salt_cavity)
      !! Fill the two OWNED surface-flux components from the solved
      !! interface.  FULL OVERWRITE of the whole plane, never `+=` —
      !! the same rule the sea-ice coupler's fillers follow, and for the
      !! same reason (a `+=` ratchets across outer steps with no bound).
      !!
      !! See the module docstring for the sign derivation.  In one line:
      !! the ocean loses the heat the interface takes (`-q_ocean`), and
      !! the dilution by meltwater of salinity `s_ice` is emulated at
      !! fixed column mass by `-m_mass*(S_far - s_ice)`.
      !!
      !! No wet-mask factor here: `active` already carries it (the
      !! sampler composes it), and the assembler multiplies `Q_heat` /
      !! `Q_salt` by `ms%wet_mask` again.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: s_ice
         !! Ice salinity (g/kg).
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      real(wp), intent(in) :: melt(nx, ny)
         !! Melt mass flux (kg/m^2/s), > 0 melting.
      real(wp), intent(in) :: q_ocean(nx, ny)
         !! Ocean → interface heat flux (W/m^2).
      real(wp), intent(in) :: s_far(nx, ny)
         !! Far-field salinity the solve used (g/kg).
      real(wp), intent(out) :: heat_cavity(nx, ny)
         !! Heat component (W/m^2), positive down into the ocean.
      real(wp), intent(out) :: salt_cavity(nx, ny)
         !! Virtual salt component, positive salinifies.
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            heat_cavity(i, j) = -q_ocean(i, j)
            salt_cavity(i, j) = -melt(i, j)*(s_far(i, j) - s_ice)
         else
            heat_cavity(i, j) = 0.0_wp
            salt_cavity(i, j) = 0.0_wp
         end if
      end do
   end subroutine cavity_flux_fill_impl

   ! ======================================================================
   ! Status accounting
   ! ======================================================================

   pure subroutine cavity_status_counts_impl(nx, ny, active, status, n_nonfinite, &
                                             n_not_converged, n_no_root, n_other)
      !! Reduce the per-column status plane into four counts, ON DEVICE.
      !!
      !! `do concurrent ... reduce(+:)` and not `count(...)`: the array
      !! is device-resident under `mem:separate`, and a host intrinsic
      !! would silently read the stale host shadow — the same trap
      !! `continuity_t%n_limited_step` documents.
      !!
      !! Only ACTIVE columns are counted.  An inactive column always
      !! carries `CAVITY_MELT_OK`, but gating the count on the mask
      !! keeps that a property of the caller rather than of the kernel.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      integer, intent(in) :: status(nx, ny)
         !! `CAVITY_MELT_*` per column.
      integer, intent(out) :: n_nonfinite
         !! Count of `NONFINITE_INPUT` + `NONFINITE_STATE`.
      integer, intent(out) :: n_not_converged
         !! Count of `NOT_CONVERGED`.
      integer, intent(out) :: n_no_root
         !! Count of `NO_PHYSICAL_ROOT`.
      integer, intent(out) :: n_other
         !! Count of every other non-OK status.
      integer :: i, j, c_nf, c_nc, c_nr, c_ot

      c_nf = 0
      c_nc = 0
      c_nr = 0
      c_ot = 0
      do concurrent(j=1:ny, i=1:nx) reduce(+:c_nf, c_nc, c_nr, c_ot)
         if (active(i, j) > 0.5_wp) then
            select case (status(i, j))
            case (CAVITY_MELT_OK)
               continue
            case (CAVITY_MELT_NONFINITE_INPUT, CAVITY_MELT_NONFINITE_STATE)
               c_nf = c_nf + 1
            case (CAVITY_MELT_NOT_CONVERGED)
               c_nc = c_nc + 1
            case (CAVITY_MELT_NO_PHYSICAL_ROOT)
               c_nr = c_nr + 1
            case default
               c_ot = c_ot + 1
            end select
         end if
      end do
      n_nonfinite = c_nf
      n_not_converged = c_nc
      n_no_root = c_nr
      n_other = c_ot
   end subroutine cavity_status_counts_impl

   pure function cavity_melt_status_is_fatal(n_nonfinite) result(is_fatal)
      !! Is this step's status tally a FAIL-LOUD condition?
      !!
      !! A non-finite input or intermediate on a covered column means the
      !! state feeding the interface is already corrupt; the kernel's
      !! safe state (zero melt) would hide it behind a plausible run.
      !! `NOT_CONVERGED` and `NO_PHYSICAL_ROOT` are NOT fatal — they are
      !! counted and warned, and those columns take the documented
      !! zero-melt safe state.
      !!
      !! A `pure` predicate so the test suite can assert the DECISION
      !! without provoking `error stop`, the way the repo tests every
      !! other fail-loud rule.
      integer, intent(in) :: n_nonfinite
         !! `n_nonfinite_step` from `cavity_status_counts_impl`.
      logical :: is_fatal
      is_fatal = (n_nonfinite > 0)
   end function cavity_melt_status_is_fatal

   ! ======================================================================
   ! Phase 3 — real freshwater MASS
   ! ======================================================================

   pure subroutine cavity_mass_apply_impl(nx, ny, nz, dt_over_rho0, inv_rho0_dt, s_ice, &
                                          active, melt, s_far, t_b, &
                                          h_layer, hTr_S, hTr_T, &
                                          salt_budget, heat_budget, k_top, n_thin)
      !! The real-mass top-layer source: add the meltwater VOLUME to
      !! `h_layer` at the first LIVE layer `k_top(i,j)`, replace the
      !! virtual salt flux the assembler
      !! stamped with the real advective salt `w*s_ice`, and add the
      !! enthalpy `dh*T_b` — mirroring both tracer increments into the
      !! existing surface budget contributors.
      !!
      !! See the module docstring for the derivation.  In one line:
      !! `dh = m*dt/rho_0`, `d(hS) = dh*s_ice`, `d(hT) = dh*T_b`, on top
      !! of the `-q_ocean` the assembler already delivers through
      !! `Q_heat`.
      !!
      !! `dt_over_rho0` and `inv_rho0_dt` are the SAME scalar, passed
      !! twice on purpose: the first multiplies `melt`, the second is the
      !! `dt/rho_0` the surface-flux apply used, and the salt undo is
      !! only the EXACT negation of that stamp if the two are the same
      !! bit pattern.  Passing one argument and reusing it makes that
      !! structural rather than a comment.
      !!
      !! FREEZING (`m < 0`) withdraws.  A column whose top layer cannot
      !! give up `|dh|` without falling through `H_CAVITY_FLOOR` is
      !! CLAMPED to what is there above that floor and COUNTED; the
      !! caller treats a non-zero count as fatal
      !! (`cavity_mass_thin_is_fatal`), because a clamped withdrawal no
      !! longer matches the tracked mass source.  The clamp is a plain
      !! `max` on a quantity that has already been range-tested, not a
      !! NaN-laundering `if/else` chain: a non-finite `melt` is caught
      !! upstream by the solver's fatal status.
      !!
      !! The floor is `H_CAVITY_FLOOR`, not `H_VANISHED`, and that gap is
      !! load-bearing: the abort fires AFTER this kernel has written the
      !! state, so whatever it leaves behind gets read at least once more.
      !! A layer left exactly ON the marker reads as VANISHED to every
      !! strict-`>` gate downstream, and `ocean_remap_tracer_column` then
      !! returns `hTr = 0` for it — the heat and salt the clamp was
      !! protecting are deleted by the very next remap.  See
      !! `H_CAVITY_FLOOR`.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count; `k = nz` is the surface / ice-base layer.
      real(wp), intent(in) :: dt_over_rho0
         !! `dt/rho_0` (m^3 s / kg) — multiplies `melt` to give `dh`.
      real(wp), intent(in) :: inv_rho0_dt
         !! The same `dt/rho_0`, used to undo the surface-flux stamp.
      real(wp), intent(in) :: s_ice
         !! Ice salinity (g/kg).
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask; 1 only on covered, wet, sampled columns.
      real(wp), intent(in) :: melt(nx, ny)
         !! Melt mass flux (kg/m^2/s), > 0 melting.
      real(wp), intent(in) :: s_far(nx, ny)
         !! Far-field salinity (g/kg) the solve used.  The VIRTUAL salt
         !! component `-melt*(s_far - s_ice)` is REBUILT from it here,
         !! with the same expression `cavity_flux_fill_impl` uses, rather
         !! than read back off `sf`: on a covered column the assembler's
         !! open-water factor is exactly zero, so `Q_salt` IS that
         !! component bit for bit, and rebuilding it keeps this kernel
         !! free of the surface-flux slot.
      real(wp), intent(in) :: t_b(nx, ny)
         !! Interface temperature (degC): the temperature the meltwater
         !! joins the column at.
      real(wp), intent(inout) :: h_layer(nx, ny, nz)
         !! Layer thickness (m); only `k = k_top(i,j)` is touched.
      real(wp), intent(inout) :: hTr_S(nx, ny, nz)
         !! `h*S` ((g/kg) m).
      real(wp), intent(inout) :: hTr_T(nx, ny, nz)
         !! `h*T` (degC m).
      real(wp), intent(inout) :: salt_budget(nx, ny, nz)
         !! `ms%salt_budget_surface`.
      real(wp), intent(inout) :: heat_budget(nx, ny, nz)
         !! `ms%heat_budget_surface`.
      integer, intent(in) :: k_top(nx, ny)
         !! `ms%k_top` — the first LIVE layer counting down from the
         !! top, `nz` wherever nothing vanishes against the top.  Under
         !! a quasi-geopotential coordinate beneath the shelf `k = nz` is
         !! an inert filler on every covered column: growing IT would put
         !! the meltwater where the remap drain deletes it, and shrinking
         !! it would trip the fatal thin-withdrawal clamp on every column
         !! at once (a filler is already AT the marker).
      integer, intent(out) :: n_thin
         !! Columns whose withdrawal had to be clamped.
      integer :: i, j, kt, c_thin
      real(wp) :: dh, h0, hn, cell_s, cell_t, virt

      c_thin = 0
      do concurrent(j=1:ny, i=1:nx) local(kt, dh, h0, hn, cell_s, cell_t, virt) reduce(+:c_thin)
         if (active(i, j) > 0.5_wp) then
            kt = k_top(i, j)
            h0 = h_layer(i, j, kt)
            dh = melt(i, j)*dt_over_rho0
            hn = h0 + dh
            if (dh < 0.0_wp .and. hn < H_CAVITY_FLOOR) then
               ! Pin the RESULT at the floor and derive the applied `dh`
               ! from it, rather than pinning `dh` and adding: `h0 +
               ! (H_CAVITY_FLOOR - h0)` rounds and can land a ulp BELOW
               ! the floor, which is the one place a thin-layer gate must
               ! never be.
               !
               ! `min(h0, ...)` takes only what sits ABOVE the floor: a
               ! column already at or below it has its withdrawal REFUSED
               ! (`hn = h0`, `dh = 0`) rather than turned into a deposit.
               ! Clamping `hn` UP to the floor there would invent mass the
               ! tracked source does not name — a worse failure than the
               ! one being prevented, and it is what the unconditional
               ! `hn = marker` form used to do.
               !
               ! Gated on `dh < 0` for the same reason: a MELTING column
               ! (`dh >= 0`) is adding mass, and a thin result then means
               ! the layer was already thin on arrival — not this kernel's
               ! doing, and not this kernel's to fabricate away.
               hn = min(h0, H_CAVITY_FLOOR)
               dh = hn - h0
               c_thin = c_thin + 1
            end if
            virt = -melt(i, j)*(s_far(i, j) - s_ice)
            cell_s = -inv_rho0_dt*virt + dh*s_ice
            cell_t = dh*t_b(i, j)
            h_layer(i, j, kt) = hn
            hTr_S(i, j, kt) = hTr_S(i, j, kt) + cell_s
            salt_budget(i, j, kt) = salt_budget(i, j, kt) + cell_s
            hTr_T(i, j, kt) = hTr_T(i, j, kt) + cell_t
            heat_budget(i, j, kt) = heat_budget(i, j, kt) + cell_t
         end if
      end do
      n_thin = c_thin
   end subroutine cavity_mass_apply_impl

   pure subroutine cavity_mass_salt_mirror_impl(nx, ny, nz, dt_over_rho0, inv_rho0_dt, &
                                                s_ice, active, melt, s_far, hTr, k_top)
      !! The pseudo-salt mirror of `cavity_mass_apply_impl`'s SALT
      !! increment, with no budget accumulation — `budget_id = NONE`, so
      !! it must not touch `salt_budget_surface`.
      !!
      !! Pseudo-salt's documented purpose is to be a PASSIVE tracer given
      !! exactly salinity's surface salt flux, so the deviation measures
      !! the passive-vs-active transport-path error.  Under the real-mass
      !! form salinity's surface salt flux is `w*s_ice` (and the dilution
      !! by the added volume, which pseudo-salt gets for free because it
      !! shares `h`), so the mirror must receive the same correction — or
      !! it would keep the virtual flux salinity no longer has and the
      !! deviation would measure the bookkeeping instead of the transport.
      !!
      !! The arithmetic is the SAME EXPRESSION as the salt branch above,
      !! so the increment is bit-identical to salinity's.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count.
      real(wp), intent(in) :: dt_over_rho0
         !! `dt/rho_0`.
      real(wp), intent(in) :: inv_rho0_dt
         !! The same `dt/rho_0` (see the sibling's docstring).
      real(wp), intent(in) :: s_ice
         !! Ice salinity (g/kg).
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      real(wp), intent(in) :: melt(nx, ny)
         !! Melt mass flux (kg/m^2/s).
      real(wp), intent(in) :: s_far(nx, ny)
         !! Far-field salinity (g/kg); the virtual component being
         !! undone is rebuilt from it, exactly as in the sibling.
      real(wp), intent(inout) :: hTr(nx, ny, nz)
         !! Pseudo-salt `h*C`.
      integer, intent(in) :: k_top(nx, ny)
         !! The SAME first-live-layer index the salt branch used — the
         !! mirror's whole contract is that its increment is bit-identical
         !! to salinity's, which includes landing on the same row.
      integer :: i, j
      real(wp) :: cell_s, virt

      do concurrent(j=1:ny, i=1:nx) local(cell_s, virt)
         if (active(i, j) > 0.5_wp) then
            virt = -melt(i, j)*(s_far(i, j) - s_ice)
            cell_s = -inv_rho0_dt*virt + melt(i, j)*dt_over_rho0*s_ice
            hTr(i, j, k_top(i, j)) = hTr(i, j, k_top(i, j)) + cell_s
         end if
      end do
   end subroutine cavity_mass_salt_mirror_impl

   pure subroutine cavity_mass_totals_impl(nx, ny, nghost, dt_over_rho0, active, &
                                           wet_mask, cover_frac, melt, areaT, &
                                           vol_melt, area_open)
      !! The two INTERIOR integrals the mass budget and the compensation
      !! sink need, reduced on device in one pass:
      !!
      !!   `vol_melt`  — the meltwater volume (m^3) added this step,
      !!                 `sum(melt*dt/rho_0 * areaT)` over active columns;
      !!   `area_open` — the wet area (m^2) the ice does NOT cover.
      !!
      !! GHOSTS ARE EXCLUDED (`nghost+1 .. n-nghost`), exactly as
      !! `ocean_accumulate_mass_out` excludes them: the console's totals
      !! are interior integrals, so a source term that counted the halo
      !! would be measuring a different domain than the total it has to
      !! close.  The per-cell APPLY still runs over the whole array, like
      !! every other surface kernel — the halo exchange owns the ghosts.
      !!
      !! `melt*dt_over_rho0` is formed with the same scalar and in the
      !! same order as the apply kernel's `dh`, so the two agree to the
      !! last bit per cell and the residual is reduction round-off only.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nghost
         !! Halo width to exclude on every side.
      real(wp), intent(in) :: dt_over_rho0
         !! `dt/rho_0`.
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Static wet (1) / land (0) mask.
      real(wp), intent(in) :: cover_frac(nx, ny)
         !! Ice-cover fraction (v1 binary).
      real(wp), intent(in) :: melt(nx, ny)
         !! Melt mass flux (kg/m^2/s).
      real(wp), intent(in) :: areaT(nx, ny)
         !! Cell area (m^2).
      real(wp), intent(out) :: vol_melt
         !! Meltwater volume this step (m^3).
      real(wp), intent(out) :: area_open
         !! Uncovered wet area (m^2).
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: acc_v, acc_a

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      acc_v = 0.0_wp
      acc_a = 0.0_wp
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) reduce(+:acc_v, acc_a)
         if (active(i, j) > 0.5_wp) then
            acc_v = acc_v + melt(i, j)*dt_over_rho0*areaT(i, j)
         end if
         if (wet_mask(i, j) > 0.5_wp .and. cover_frac(i, j) < 0.5_wp) then
            acc_a = acc_a + areaT(i, j)
         end if
      end do
      vol_melt = acc_v
      area_open = acc_a
   end subroutine cavity_mass_totals_impl

   pure function cavity_comp_withdrawal(vol_melt, area_open) result(dw)
      !! The uniform per-unit-area thickness the `uniform_open_ocean`
      !! sink removes: `vol_melt/area_open`, and exactly zero when there
      !! is no open ocean to remove it from (a fully ice-covered domain —
      !! the sink then does nothing rather than dividing by zero, and the
      !! volume stays in, which the mass budget reports honestly as a
      !! growing total).
      !!
      !! `pure`, separate from the kernels, so the test suite can assert
      !! the arithmetic and the no-open-ocean fallback directly.
      real(wp), intent(in) :: vol_melt
         !! Meltwater volume this step (m^3); either sign.
      real(wp), intent(in) :: area_open
         !! Uncovered wet area (m^2), >= 0.
      real(wp) :: dw
      if (area_open > H_DIV_EPS) then
         dw = vol_melt/area_open
      else
         dw = 0.0_wp
      end if
   end function cavity_comp_withdrawal

   pure subroutine cavity_comp_apply_impl(nx, ny, nz, dw, wet_mask, cover_frac, &
                                          h_layer, hTr_S, hTr_T, &
                                          salt_budget, heat_budget, scale, n_thin)
      !! The `volume_compensation = "uniform_open_ocean"` sink: remove
      !! `dw` metres of the top layer from every wet cell the ice does
      !! NOT cover, the removed parcel carrying that cell's own `T` and
      !! `S` so no concentration there changes.
      !!
      !! "Carries its own concentration" is implemented as a RATIO:
      !! `hTr *= h_new/h_old`, which leaves `Tr = hTr/h` algebraically
      !! identical.  The ratio is published in `scale` so every PASSIVE
      !! tracer can be given exactly the same treatment
      !! (`cavity_comp_scale_tracer_impl`) — a tracer left alone here
      !! would be CONCENTRATED by the sink, which is a different physical
      !! statement from the one the knob makes.
      !!
      !! The salt and heat increments ARE mirrored into the surface
      !! budget contributors: a parcel that carries `S` out of the domain
      !! changes the domain salt total, so the budget has to name it.
      !!
      !! Same clamp (`H_CAVITY_FLOOR`, strictly above the vanish marker)
      !! + count + fatal policy as the mass source.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count; only `k = nz` is touched.
      real(wp), intent(in) :: dw
         !! Thickness (m) to remove from each open-ocean column.
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Static wet (1) / land (0) mask.
      real(wp), intent(in) :: cover_frac(nx, ny)
         !! Ice-cover fraction (v1 binary).
      real(wp), intent(inout) :: h_layer(nx, ny, nz)
         !! Layer thickness (m).
      real(wp), intent(inout) :: hTr_S(nx, ny, nz)
         !! `h*S`.
      real(wp), intent(inout) :: hTr_T(nx, ny, nz)
         !! `h*T`.
      real(wp), intent(inout) :: salt_budget(nx, ny, nz)
         !! `ms%salt_budget_surface`.
      real(wp), intent(inout) :: heat_budget(nx, ny, nz)
         !! `ms%heat_budget_surface`.
      real(wp), intent(out) :: scale(nx, ny)
         !! `h_new/h_old` where the sink acted, exactly 1 elsewhere.
      integer, intent(out) :: n_thin
         !! Columns whose withdrawal had to be clamped.
      integer :: i, j, c_thin
      real(wp) :: h0, hn, d, f, cell_s, cell_t

      c_thin = 0
      do concurrent(j=1:ny, i=1:nx) local(h0, hn, d, f, cell_s, cell_t) reduce(+:c_thin)
         if (wet_mask(i, j) > 0.5_wp .and. cover_frac(i, j) < 0.5_wp) then
            h0 = h_layer(i, j, nz)
            d = dw
            hn = h0 - d
            if (hn < H_CAVITY_FLOOR) then
               ! Same floor and same refusal rule as the mass source: pin
               ! the RESULT (never the increment — `h0 - (h0 - floor)`
               ! rounds and can land a ulp below), take only what sits
               ! ABOVE the floor, and take nothing at all from a column
               ! already at or below it (`min(h0, ...)`, so `hn <= h0`
               ! always and the sink can never become a deposit).
               hn = min(h0, H_CAVITY_FLOOR)
               c_thin = c_thin + 1
            end if
            if (h0 > H_DIV_EPS) then
               f = hn/h0
            else
               f = 1.0_wp
            end if
            cell_s = hTr_S(i, j, nz)*(f - 1.0_wp)
            cell_t = hTr_T(i, j, nz)*(f - 1.0_wp)
            h_layer(i, j, nz) = hn
            hTr_S(i, j, nz) = hTr_S(i, j, nz) + cell_s
            salt_budget(i, j, nz) = salt_budget(i, j, nz) + cell_s
            hTr_T(i, j, nz) = hTr_T(i, j, nz) + cell_t
            heat_budget(i, j, nz) = heat_budget(i, j, nz) + cell_t
            scale(i, j) = f
         else
            scale(i, j) = 1.0_wp
         end if
      end do
      n_thin = c_thin
   end subroutine cavity_comp_apply_impl

   pure subroutine cavity_comp_scale_tracer_impl(nx, ny, nz, scale, hTr)
      !! Apply the compensation sink's top-layer ratio to one PASSIVE
      !! tracer's load.  Unconditional (`scale = 1` off the sink), so
      !! there is no mask branch inside the kernel.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count.
      real(wp), intent(in) :: scale(nx, ny)
         !! `h_new/h_old` from `cavity_comp_apply_impl`.
      real(wp), intent(inout) :: hTr(nx, ny, nz)
         !! The tracer's `h*C`.
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         hTr(i, j, nz) = hTr(i, j, nz)*scale(i, j)
      end do
   end subroutine cavity_comp_scale_tracer_impl

   pure function cavity_mass_thin_is_fatal(n_thin) result(is_fatal)
      !! Is this step's clamped-withdrawal tally a FAIL-LOUD condition?
      !!
      !! Yes, on any count.  A clamped withdrawal is one the tracked mass
      !! source no longer matches, so from that step on the console's
      !! `Mass Error` would stop being a leak measurement while still
      !! being printed as one — the exact failure mode the budget gate
      !! exists to make impossible.  `h` is nevertheless clamped first,
      !! so the state handed to the abort path is finite and the last
      !! printed line is meaningful.
      !!
      !! A `pure` predicate, like `cavity_melt_status_is_fatal`, so the
      !! suite can assert the DECISION without provoking `error stop`.
      integer, intent(in) :: n_thin
         !! `n_thin_step`.
      logical :: is_fatal
      is_fatal = (n_thin > 0)
   end function cavity_mass_thin_is_fatal

   subroutine ocean_cavity_mass_step(grid, metrics, cav, ms, dt, weight, active)
      !! The real-freshwater MASS update — one call, at the THERMO
      !! cadence, from inside the RK2 stage immediately after
      !! `ocean_surface_flux_apply_tracers`.
      !!
      !! WHY THERE, and not in `engine_step_finalize` next to the melt
      !! solve.  Three things have to line up:
      !!
      !!   1. **The same melt value.**  The salt and heat halves of the
      !!      melt ride `Q_salt`/`Q_heat`, which the stage's surface-flux
      !!      apply spends.  Putting the volume anywhere else would spend
      !!      `melt(n)` for mass and `melt(n-1)` for salt.
      !!   2. **The same stage weight, by construction.**  Under
      !!      `ssp_rk2` both stages apply and `rk2_average` halves the
      !!      pair; under `pred_corr` `therm_active` is false in the
      !!      predictor and the corrector's single application IS the
      !!      step.  The caller passes the matching `weight` for the
      !!      budget accumulator, the same one `ocean_accumulate_mass_out`
      !!      takes.
      !!   3. **The barotropic mode sees it.**  `derive_bt_from_layers`
      !!      rebuilds `bt_eta = sum_k h - bt_H_ref` at the TOP of every
      !!      stage, so a thickness source applied here is in `bt_eta`
      !!      one stage later with no separate barotropic forcing term —
      !!      the free surface under the cavity datum simply rises.  The
      !!      BT substep's transport renormalisation (`bt_uhbt`) is a
      !!      constraint on the layer TRANSPORTS within a stage and never
      !!      reads a thickness source, so it is undisturbed.  This is
      !!      the documented `F_slow`-style operator split, and it is what
      !!      MOM6 does with its own surface mass fluxes.
      !!
      !! It runs BEFORE the vertical-mixing block in the same stage, so
      !! the implicit vdiff/drag solves see the thickened top layer —
      !! which is the right order: the added volume is part of the column
      !! before the column is mixed.  The ALE remap then redistributes it
      !! to the coordinate's target after `rk2_average`.
      !!
      !! NON-`pure`: it logs, it fails loud, it mutates `ms%mass_src`,
      !! and (with compensation on) it issues one collective.  Its six
      !! kernels are all `pure`.
      !!
      !! `mem:separate`: every array is mapped by its owning slot's
      !! `enter_data`; the tracer-registry derefs happen on the HOST
      !! before each kernel (the outer-shim rule).
      use pic_logger, only: global_logger
      use pic_strings, only: to_string
      use rdb_error_ring, only: fail
      use rdb_ocean_status, only: OCEAN_STATUS_ERR_SETUP
      use rdb_halo, only: halo_allreduce_sum
      type(hgrid_t), intent(in) :: grid
         !! Horizontal grid (`nx_total`/`ny_total`/`nghost`).
      type(ocean_metrics_t), intent(in) :: metrics
         !! Reads `cover_frac` and `areaT`.
      type(ocean_cavity_flux_t), intent(inout) :: cav
         !! The cavity-melt slot; reads `melt`/`active`/`t_b`, writes
         !! `comp_scale` and the thin counters.
      type(multilayer_state_t), intent(inout) :: ms
         !! Writes `h_layer(:,:,nz)`, the S/T/passive tracer loads, the
         !! two surface budget contributors and `mass_src`.
      real(wp), intent(in) :: dt
         !! Thermo timestep (s) — the same `therm_dt` the surface-flux
         !! apply was given.
      real(wp), intent(in) :: weight
         !! Per-stage weight for the `mass_src` accumulator (0.5 per
         !! SSP-RK2 stage; 0 / 1 for the pred_corr predictor /
         !! corrector), matching `ocean_accumulate_mass_out`.
      logical, intent(in), optional :: active
         !! Thermo-cadence gate.  Present-and-false ⇒ early return;
         !! absent ⇒ run.
      integer :: nx, ny, nz, idx_t, idx_s, idx_ps, it, n_thin_src, n_thin_comp
      real(wp) :: dt_over_rho0, vol_melt, area_open, dw, tmp

      if (.not. cav%is_init) return
      if (.not. cav%enable) return
      if (cav%freshwater /= CAVITY_FW_MASS) return
      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      if (idx_t <= 0 .or. idx_s <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      idx_ps = ms%idx_pseudo_salt
      dt_over_rho0 = dt/cav%rho0

      ! (1) Interior integrals FIRST: the melt volume the budget tracks
      ! and the open-ocean area the sink would spend it over.  Combined
      ! across ranks so a decomposed domain removes the same total the
      ! whole domain gained.
      call cavity_mass_totals_impl(nx, ny, grid%nghost, dt_over_rho0, cav%active, &
                                   ms%wet_mask, metrics%cover_frac, cav%melt, &
                                   metrics%areaT, vol_melt, area_open)
      tmp = vol_melt
      call halo_allreduce_sum(tmp, vol_melt)
      tmp = area_open
      call halo_allreduce_sum(tmp, area_open)
      cav%melt_volume_step = vol_melt
      cav%open_area = area_open

      ! (2) The source.
      call cavity_mass_apply_impl(nx, ny, nz, dt_over_rho0, dt_over_rho0, cav%s_ice, &
                                  cav%active, cav%melt, cav%s_far, cav%t_b, &
                                  ms%h_layer, ms%tracers(idx_s)%hTr, ms%tracers(idx_t)%hTr, &
                                  ms%salt_budget_surface, ms%heat_budget_surface, &
                                  ms%k_top, n_thin_src)
      if (idx_ps > 0) then
         call cavity_mass_salt_mirror_impl(nx, ny, nz, dt_over_rho0, dt_over_rho0, &
                                           cav%s_ice, cav%active, cav%melt, &
                                           cav%s_far, ms%tracers(idx_ps)%hTr, ms%k_top)
      end if
      ms%mass_src = ms%mass_src + weight*RHO_WATER*vol_melt

      ! (3) The sink, when the knob asks for it.
      n_thin_comp = 0
      cav%comp_withdrawal_step = 0.0_wp
      if (cav%volume_comp == CAVITY_VC_UNIFORM_OPEN) then
         dw = cavity_comp_withdrawal(vol_melt, area_open)
         cav%comp_withdrawal_step = dw
         call cavity_comp_apply_impl(nx, ny, nz, dw, ms%wet_mask, metrics%cover_frac, &
                                     ms%h_layer, ms%tracers(idx_s)%hTr, &
                                     ms%tracers(idx_t)%hTr, ms%salt_budget_surface, &
                                     ms%heat_budget_surface, cav%comp_scale, n_thin_comp)
         ! NOTE the sink stays on `k = nz` and is NOT routed through
         ! `k_top`, deliberately: its own gate is `cover_frac < 0.5`, so
         ! it only ever acts on OPEN-OCEAN columns, and an open-ocean
         ! column has `z_top = 0` ⇒ no top-side filler ⇒ `k_top ≡ nz`.
         ! Routing it would be a provable no-op; leaving it spells out
         ! that "the top layer" and "the first live layer" are the same
         ! row wherever this kernel runs.  `cavity_comp_scale_tracer_impl`
         ! below rides the same argument (it is unconditional, but
         ! `comp_scale` is exactly 1 off the sink).
         do it = 1, size(ms%tracers)
            if (it == idx_s .or. it == idx_t) cycle
            call cavity_comp_scale_tracer_impl(nx, ny, nz, cav%comp_scale, &
                                               ms%tracers(it)%hTr)
         end do
         ! The removed interior volume is exactly `dw*area_open` — the
         ! same product the withdrawal was derived from — so with no
         ! clamping this cancels the source term to the last bit and the
         ! domain mass is constant.
         ms%mass_src = ms%mass_src - weight*RHO_WATER*dw*area_open
      end if

      ! (4) Accounting, then fail loud.
      cav%n_thin_step = n_thin_src + n_thin_comp
      cav%n_thin_total = cav%n_thin_total + int(cav%n_thin_step, int64)
      if (cavity_mass_thin_is_fatal(cav%n_thin_step)) then
         call global_logger%error("cavity real freshwater: "// &
                                  to_string(n_thin_src)//" melt column(s) and "// &
                                  to_string(n_thin_comp)//" compensation column(s) "// &
                                  "could not give up the requested thickness")
         call fail("&ocean_cavity_melt_nml freshwater='mass': "// &
                   to_string(cav%n_thin_step)//" column(s) would have been driven "// &
                   "below H_CAVITY_FLOOR (2*H_VANISHED) by the top-layer "// &
                   "withdrawal (a freezing "// &
                   "column thinner than |m|*dt/rho_0, or a compensation sink "// &
                   "deeper than the open-ocean top layer).  The withdrawal was "// &
                   "clamped so the state stays finite, but a clamped withdrawal "// &
                   "no longer matches the tracked mass source, so the console's "// &
                   "Mass Error would stop being a leak measurement while still "// &
                   "being printed as one.  Use a thicker top layer (a coordinate "// &
                   "with a larger surface target), a shorter dt, or "// &
                   "volume_compensation='none'.", code=OCEAN_STATUS_ERR_SETUP)
      end if
   end subroutine ocean_cavity_mass_step

   ! ======================================================================
   ! The driver
   ! ======================================================================

   subroutine ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf, active)
      !! One cavity basal-melt update: sample the far field, solve the
      !! three-equation interface on every covered column, deliver the
      !! two owned surface-flux components, and account for the solver
      !! status.
      !!
      !! CADENCE + PLACEMENT.  Called once per outer step at THERMO
      !! cadence from `engine_step_finalize`, immediately BEFORE
      !! `ocean_surface_flux_assemble` — the assembler must see the
      !! components this routine writes.  The tracers then integrate the
      !! assembled `Q_heat`/`Q_salt` on the NEXT outer step, which is the
      !! same one-step lag the sea-ice coupler documents.
      !!
      !! NON-`pure`, and it is the only non-pure procedure this module
      !! adds: it logs (`n_not_converged` warnings) and it FAILS LOUD on
      !! a non-finite column.  Its three kernels are all `pure`.
      !!
      !! `mem:separate`: every array read or written here is mapped by
      !! its owning slot's `enter_data`.  `ms%tracers(idx)%hTr` is
      !! dereferenced on the HOST before the kernels (the outer-shim
      !! rule for the array-of-derived-types registry).
      use pic_logger, only: global_logger
      use pic_strings, only: to_string
      use rdb_error_ring, only: fail
      use rdb_ocean_status, only: OCEAN_STATUS_ERR_SETUP
      type(hgrid_t), intent(in) :: grid
         !! Horizontal grid (for `nx_total`/`ny_total`).
      type(ocean_cavity_flux_t), intent(inout) :: cav
         !! The cavity-melt slot.
      type(ocean_metrics_t), intent(in) :: metrics
         !! Reads `cover_frac` only.
      type(multilayer_state_t), intent(in) :: ms
         !! Reads `h_layer`, the S/T tracer loads, the face velocities,
         !! `p_top` (THE interface pressure) and `wet_mask`.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.  Flat POD, by value.
      type(ocean_surface_flux_t), intent(inout) :: sf
         !! Writes `heat_cavity`/`salt_cavity` and latches
         !! `has_heat`/`has_salt` host-side.
      logical, intent(in), optional :: active
         !! Thermo-cadence gate.  Present-and-false ⇒ early return;
         !! absent ⇒ run (the `ocean_surface_flux_assemble` convention).
      integer :: nx, ny, nz, idx_t, idx_s

      if (.not. cav%is_init) return
      if (.not. cav%enable) return
      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      if (idx_t <= 0 .or. idx_s <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call cavity_far_field_impl(nx, ny, nz, cav%far_field_depth, &
                                 metrics%cover_frac, ms%wet_mask, ms%h_layer, &
                                 ms%tracers(idx_t)%hTr, ms%tracers(idx_s)%hTr, &
                                 ms%u_face_x_layer, ms%v_face_y_layer, &
                                 cav%active, cav%t_far, cav%s_far, cav%u_far, cav%v_far)

      call cavity_melt_columns_2d(nx, ny, cav%active, cav%t_far, cav%s_far, ms%p_top, &
                                  cav%u_far, cav%v_far, cav%s_ice, cav%f_cor, &
                                  cav%cdrag_top, cav%u_tide, cav%ustar_min, &
                                  cav%par, cav%ice, eos, cav%const, &
                                  cav%ustar, cav%t_b, cav%s_b, cav%melt, cav%q_ocean, &
                                  cav%gamma_t, cav%gamma_s, cav%status)

      call cavity_flux_fill_impl(nx, ny, cav%s_ice, cav%active, cav%melt, cav%q_ocean, &
                                 cav%s_far, sf%heat_cavity, sf%salt_cavity)

      ! The filler's own obligations (PR-12 fill contract): latch the
      ! has_* flags HOST-side, never from a device reduction.  They stay
      ! latched for the rest of the run — the melt rate can legitimately
      ! pass through zero, and un-latching there would change the
      ! apply-path operand count mid-run.
      sf%has_heat = .true.
      sf%has_salt = .true.

      call cavity_status_counts_impl(nx, ny, cav%active, cav%status, &
                                     cav%n_nonfinite_step, cav%n_not_converged_step, &
                                     cav%n_no_root_step, cav%n_other_step)
      cav%n_not_converged_total = cav%n_not_converged_total &
                                  + int(cav%n_not_converged_step, int64)
      cav%n_no_root_total = cav%n_no_root_total + int(cav%n_no_root_step, int64)
      cav%n_other_total = cav%n_other_total + int(cav%n_other_step, int64)

      if (cavity_melt_status_is_fatal(cav%n_nonfinite_step)) then
         call fail("&ocean_cavity_melt_nml: "//to_string(cav%n_nonfinite_step)// &
                   " ice-covered column(s) fed the basal-melt solver a non-finite "// &
                   "far-field temperature, salinity, velocity or interface pressure "// &
                   "(CAVITY_MELT_NONFINITE_*).  The kernel's safe state is zero melt, "// &
                   "which would hide an already-corrupt column behind a plausible "// &
                   "run, so this is fatal.  NOT_CONVERGED / NO_PHYSICAL_ROOT are "// &
                   "counted and warned instead.", code=OCEAN_STATUS_ERR_SETUP)
      end if

      if (cav%n_not_converged_step > 0 .or. cav%n_no_root_step > 0 .or. &
          cav%n_other_step > 0) then
         call global_logger%warning("cavity melt: zero melt applied on "// &
                                    to_string(cav%n_not_converged_step)// &
                                    " non-converged, "//to_string(cav%n_no_root_step)// &
                                    " no-physical-root and "// &
                                    to_string(cav%n_other_step)// &
                                    " otherwise-refused column(s) this step")
      end if
   end subroutine ocean_cavity_flux_step

end module rdb_ocean_cavity_flux
