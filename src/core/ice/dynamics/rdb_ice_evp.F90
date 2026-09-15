!! C-grid elastic-viscous-plastic (EVP) sea-ice rheology (SIS2 port, PR 5).
module rdb_ice_evp
   !! Mechanical transliteration of the validated Python prototype
   !! `tmp_local_artifacts/proto_evp_core.py` (`Channel1D`) +
   !! `tmp_local_artifacts/proto_evp_validate.py` (`Box2D`), grounded
   !! against `SIS_dyn_cgrid.F90` (`SIS_C_dynamics` :603-1614,
   !! `limit_stresses` :1619-1741, `SIS_C_dyn_init` :209-270) per
   !! `SPEC_ice-pr5-evp.md`.  Field names mirror SIS2 1:1 (`str_d`,
   !! `str_t`, `str_s`, `sh_Dd`, `sh_Dt`, `sh_Ds`, `zeta`, `del_sh`,
   !! `mi_ratio_A_q`, `Tdamp`, `EC`, ...).
   !!
   !! **Index convention** (SPEC §1): rdb u-face `(i,j)` is the WEST
   !! face of cell `(i,j)` (SIS2 face `I` = east of cell `i`); rdb
   !! v-face `(i,j)` is the SOUTH face; rdb corner `(i,j)` is the SW
   !! corner of cell `(i,j)` (SIS2 corner `(I,J)` = NE of cell `(i,j)`).
   !! The 4 T-cells around rdb corner `(ic,jc)` are `(ic-1,jc-1)
   !! (ic,jc-1) (ic-1,jc) (ic,jc)` — the `wet_q` convention.
   !!
   !! **Ice margins (Lens C, critical — no ice-edge code path anywhere).**
   !! `mi=0`/`ci=0` cells are ORDINARY wet T-cells: `pres_mice*mice=0 =>
   !! zeta=0 => str_d` decays geometrically toward 0 via the `I_1pdt_T`
   !! relaxation — emergent, never a Dirichlet special case.  Masks built
   !! here (`mask_t_w`/`mask_u_w`/`mask_v_w`/`mask_q_w`) encode LAND +
   !! non-periodic-boundary policy ONLY; they never test ice presence.
   !! The momentum-solve denominator's 0/0 guard is `m_neglect`; it is not
   !! the module's only division guard — `dxharm>0`, `denom/=0`, and the
   !! Adcroft `i_htot` reciprocal each guard their own quotient.  With PR 62
   !! `a_face_stress=.true.` the guard alone is NOT sufficient at an
   !! ice-free face: `drag_eff = a_u*drag_u` is explicitly BRANCHED around
   !! (`a_fac > 0.0`, not floored) in `evp_u_momentum_impl`/
   !! `evp_v_momentum_impl`, because `m_neglect` alone against a generally
   !! nonzero `dt*fxic_now` blows up to `O(1e30)` on the first substep.
   !!
   !! **Persistent workspace** (never local-allocate scratch in a
   !! per-substep kernel on `-stdpar=gpu`): the EVP scratch lives on the
   !! `evp_workspace_t` slot (`ice%evp_ws`, in `rdb_ice_state`), eagerly
   !! allocated in `ocean_sea_ice_t%init` and GPU-mapped via
   !! `ocean_sea_ice_t%enter_data` (which rides `ocean_state_enter_data`) —
   !! matching the rest of `src/core/ocean/` (zero module-level `save`
   !! allocatables).  `ice_evp_dynamics` is a thin shim that unpacks the
   !! slot's components into the explicit-shape flat-impl args of
   !! `ice_evp_dynamics_impl`; the whole substep loop is device-resident:
   !! no H<->D inside `do n = 1, evp_sub_steps`.
   !!
   !! **Documented divergences (D-list, SPEC §4.8):**
   !!   D1  no sea-surface-tilt term (SIS2 `PFu`/`PFv`) — v1 assumes flat
   !!       eta for the ice; a future PR wires the ocean SSH.
   !!   D2  CLOSED by PR 36 — `PROJECT_ICE_CONCENTRATION` is ported as
   !!       `&ocean_ice_nml project_ci` (SIS2 default `.true.`; Roundabout
   !!       default `.false.` ⇒ byte-identical, the house bit-identity
   !!       rule).  When on, `evp_project_ci_impl` projects `ci` (and hence
   !!       `pres_mice`) forward each subcycle from the CALL's initial
   !!       concentration and cumulative elapsed time —
   !!       `ci_proj = ci*exp(-t_cum*sh_dd)`, `t_cum = n*dt`.
   !!   D3  not ported: landfast (Lemieux/ITD, SIS2 default off),
   !!       `drag_max`/`MIN_OCN_INTERTIAL_H` (default off),
   !!       `vel_underflow`/`str_underflow` (default 0), `weak_low_shear`
   !!       (default off), `DT_RHEOLOGY` (NSTEPS_DYN only), hi-freq/sigI/
   !!       sigII diagnostics, `drag_bg_vel2` (SIS2 hardwires 0).  CFL
   !!       truncation's CFL half is CLOSED by PR 36 — `&ocean_ice_nml
   !!       cfl_trunc` (SIS2 `CFL_TRUNCATE`, default 0.5 there, `0.0` here
   !!       ⇒ byte-identical) clips the FINAL transport velocity to
   !!       `0.95*cfl_trunc*areaT(donor)/(dt_transport*dy_cu)`
   !!       (`evp_truncate_final_impl`), counting ice-bearing faces
   !!       touched into a driver-logged warning (NOT an abort — SIS2
   !!       pairs its counter with `MAXTRUNC=0`, a run-stopper Roundabout does
   !!       not port; PR-4b transport's conservation/positivity check
   !!       stays the backstop of last resort).  `cfl_trunc_dyn_its` (SIS2
   !!       `CFL_TRUNC_DYN_ITS`, default off, matches) additionally clips
   !!       to the EXACT bound at the bottom of every subcycle.  Four
   !!       documented divergences from SIS2's port of this feature: (i)
   !!       the bound uses the dt TRANSPORT will consume
   !!       (`dt_transport`, an optional argument threaded through
   !!       `ice_evp_dynamics`/`ice_evp_step`), NOT this call's `dt_slow`
   !!       — SIS2 assumes the two are the same dt, Roundabout decouples EVP
   !!       (every outer step) from transport (thermo cadence); (ii)
   !!       `ci_proj` (D2) is a `local()` scalar, not an array — SIS2
   !!       materialises it only for sigI/sigII/find_ice_strength
   !!       diagnostics Roundabout does not have; (iii) the truncation count
   !!       (`n_trunc`) drives a rank-0 driver WARNING, never an abort —
   !!       no `MAXTRUNC`; (iv) SIS2 defaults `CFL_TRUNCATE=0.5` /
   !!       `PROJECT_ICE_CONCENTRATION=.true.`, Roundabout defaults both off
   !!       (house bit-identity rule) — the shipped
   !!       `polar_freezeup_dynamics.nml` example carries SIS2's defaults
   !!       instead.  Landfast/`drag_max`/underflows/`weak_low_shear`/
   !!       `DT_RHEOLOGY`/diagnostics remain unported.
   !!   D4  v-momentum reads `u_tmp` (the PRE-update u), per SIS2
   !!       :1258-1273 — resolved toward SIS2 (the prototype's Box2D used
   !!       the updated u, a defect with no effect on any analytic gate).
   !!   D5  `ncat==1` lumped concentration: `ci = 1` where `m_ice > 0`
   !!       (SIS2 has no lumped mode) — see `ice_cell_concentration_impl`
   !!       in `rdb_ice_state` (shared with the tau coupler).
   !!   D6  domain edges wall-or-periodic only; no OBC, no tripolar fold.
   !!   D7  one atmospheric stress field: the ice feels the FULL wind
   !!       stress snapshot (`tau_a_x`/`tau_a_y`); no ice-specific bulk
   !!       drag law (that part is unchanged — a future ice-specific bulk
   !!       drag law is PR-55/RESUME #8's problem).  PR 62's
   !!       `&ocean_ice_nml a_face_stress` (default OFF) weights BOTH this
   !!       wind stress AND the ice-ocean drag in the momentum balance by
   !!       the face ice concentration `a_u = 0.5*(ci(i-1,j)+ci(i,j))`
   !!       (`evp_u_momentum_impl`/`evp_v_momentum_impl`), giving the
   !!       textbook `m du/dt = grad.sigma + a*(tau_a - tau_w)` (Hibler 1979
   !!       eq. 1) and an EXACTLY closing ice<->ocean momentum budget at
   !!       every fractional cover `a`, not just steady free drift.
   !!       `a_face_stress=.false.` (default, byte-identical) is the legacy
   !!       form: the ice absorbs the FULL wind and sheds the FULL drag
   !!       while `ice_ocean_stress_flux` hands the ocean
   !!       `(1-a)*tau_a + a*fxoc` — leaking `(1-a)*(tau_a-fxoc)` per face
   !!       per step at fractional cover (zero at `a in {0,1}` and at
   !!       steady free drift `fxoc==tau_a`, but nonzero in a generic
   !!       transient — see the F5 caveat in `ice_ocean_stress_flux`'s
   !!       docstring, `rdb_ice_ocean_coupler`).  DELIBERATE DIVERGENCE FROM
   !!       SIS2: `SIS_C_dynamics` weights NEITHER term (`fxat`/`drag_u` are
   !!       bare against a per-TOTAL-area `mis`) and carries the same leak;
   !!       the nearest SIS2 analogue, `set_wind_stresses_C`'s ice-cover-
   !!       weighted interpolation, degenerates under Roundabout's single wind
   !!       field to a pure `a>0` presence gate (kills the ghost-drift
   !!       artefact below, does not close the budget).  `a_face_stress` is
   !!       therefore more correct than the reference on Hibler (1979)/CICE
   !!       conservation grounds — the same call D4 already made toward
   !!       SIS2.  Weighting the wind ALONE (without the drag) is NOT a
   !!       valid alternative: it converts today's leak (zero at steady
   !!       free drift) into a PERMANENT one, `-(1-a)*a*tau_a`, nonzero at
   !!       steady state forever — do not "simplify" to a single term.
   !!       `a_face_stress=.true.` also kills the ghost free-drift artefact
   !!       at ice-free wet faces (today's unweighted form converges a
   !!       massless slab to the full Nansen free-drift speed
   !!       `sqrt(tau_a/(rho_o*Cdw))`, regenerated every substep; with the
   !!       knob on, `a_fac == 0` branches `uio_c` to exactly 0 => `ui==uo`).
   !!       PRE-EXISTING, NOT fixed by PR 62: the EVP's `a_u` is built from
   !!       `ci_w` (masked + periodic-wrapped); the coupler's `a_u`
   !!       (`ice_ocean_stress_flux_impl`) is built from the raw halo
   !!       (`ice_cell_concentration_impl`, no wrap).  They agree at
   !!       physical faces; at the first physical face of a periodic domain
   !!       they can differ, and the budget then closes only to that
   !!       difference there.  Lives in the coupler's gather — a future PR's
   !!       problem, not this one's.
   !!   D8  tau mediation is one-step-lagged (MEKE/frazil convention); a
   !!       fresh run's first outer step drives the ocean with pure wind
   !!       — literally true as of PR 63 (previously the resume fold's
   !!       reconstruct gave a fresh run with ice at configure
   !!       `(1-a)*tau_a` on step 1, an artefact of `fxoc` initialising to
   !!       0, not this claim; see PR 63's F4 note below and its own
   !!       plan §11.3/§14 Q5).  The ELASTIC stress state
   !!       (`str_d`/`str_t`/`str_s`) + `u_ice`/`v_ice` + `fxoc`/`fyoc`
   !!       round-trip bit-exact through the restart; `tau_x`/`tau_y`
   !!       (the field the ocean actually consumes) round-trips
   !!       bit-exact too, as of PR 63 — `ocean_sea_ice_t%tau_ocn_x/y`
   !!       mirror the exact blended value `ice_ocean_stress_flux` last
   !!       wrote and `ice_ocean_stress_resume_apply`
   !!       (`rdb_ice_ocean_coupler`) COPIES it back at configure,
   !!       formula-agnostic.  F4 (now historical): the PRE-PR-63
   !!       resume fold RECONSTRUCTED `tau_x`/`tau_y` from the
   !!       checkpoint's POST-thermo `ci`, which differs from the
   !!       PRE-thermo/pre-transport `ci` the uninterrupted run's blend
   !!       actually used whenever a checkpoint step's thermo/transport
   !!       changed `ci` after the blend — closed by carrying the
   !!       blend's output instead of recomputing it.  `tau_a_x`/`tau_a_y`
   !!       remain a fresh configure-time snapshot, not restart-carried
   !!       (D8 is otherwise unchanged: the mediation is still one-step
   !!       lagged, only the RESUME path changed).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t, ice_cell_concentration_impl, evp_workspace_t
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_y_2d
   implicit none
   private

   public :: ice_evp_params_t
   public :: ice_evp_params_from_config
   public :: ice_evp_dynamics
   public :: ice_evp_step
   public :: ice_evp_mi_ratio_point
   public :: evp_truncate_final_impl

   real(wp), parameter :: M_NEGLECT_FACTOR = 1.0e-30_wp
      !! SIS2 `H_subroundoff` (:757) — `m_neglect = ICE_RHO_ICE*1e-30`.
   real(wp), parameter :: EVP_DRAG_LINEARIZE_THRESHOLD = 1.0e8_wp
      !! Large-`b_vel0` cutover in the semi-implicit ice-ocean drag solve
      !! (SIS2 `SIS_C_dynamics`): when `b_vel0**2 >
      !! EVP_DRAG_LINEARIZE_THRESHOLD * I_cdRhoDt * |m_*io_explicit|` the
      !! quadratic drag term is negligible against the linear one, so the
      !! predicted relative velocity linearizes to `m_*io_explicit *
      !! I_cdRhoDt / b_vel0`. Same value in the u- and v-momentum kernels.
   real(wp), parameter :: TRUNC_BACKOFF = 0.95_wp
      !! PR 36: back-off factor on the FINAL CFL velocity clip (SIS2
      !! `SIS_dyn_cgrid.F90:1456`) — the clipped value is set to
      !! `0.95*bound`, not the bound itself, so a re-check cannot be
      !! marginal. The in-loop (`cfl_trunc_dyn_its`) clip uses the exact
      !! bound (backoff = 1.0) instead — it is not the last word on the
      !! velocity before transport reads it.

   type :: ice_evp_params_t
      !! EVP physical + numerical parameters (SIS2 `SIS_C_dyn_CS` subset).
      !! Passed `intent(in)` into the core; every scalar is hoisted to a
      !! local before the substep loop (no derived-type deref inside a
      !! `do concurrent`).  Defaults + meanings mirror the `&ocean_ice_nml`
      !! block in `rdb_config.F90` (the config is the authoritative knob
      !! set; keep the two in sync).
      real(wp) :: p0 = 2.75e4_wp
         !! SIS2 `ICE_STRENGTH_PSTAR` — ice-strength pressure constant [Pa].
      real(wp) :: c0 = 20.0_wp
         !! SIS2 `ICE_STRENGTH_CSTAR` — ice-strength exponent constant [nondim].
      real(wp) :: ec = 2.0_wp
         !! SIS2 `ICE_YIELD_ELLIPTICITY` — yield-curve axis ratio [nondim].
         !! 0 => cavitating-fluid rheology (`str_t`/`str_s` stay exactly 0).
      real(wp) :: cdw = 3.24e-3_wp
         !! SIS2 `ICE_CDRAG_WATER` — ice-ocean drag coefficient [nondim].
      real(wp) :: rho_ocean = 1030.0_wp
         !! SIS2 `RHO_OCEAN` — ice-drag reference density [kg/m^3].
         !! Deliberately independent of the ocean's `rho0` (usually 1035).
      real(wp) :: del_sh_min_scale = 2.0_wp
         !! SIS2 `ICE_DEL_SH_MIN_SCALE` — viscosity-floor scale [nondim].
      real(wp) :: tdamp = -0.2_wp
         !! SIS2 `ICE_TDAMP_ELASTIC` — elastic damping timescale selector.
         !! `> 0` => seconds; `== 0` => `max(0.2*dt_slow, 3*dt)`; `< 0` =>
         !! the special case `max(|tdamp|*dt_slow, 3*dt)` (i.e. `|tdamp|` is
         !! a fraction of the slow step). Sign-free — no positivity guard.
      integer :: evp_sub_steps = 432
         !! SIS2 `NSTEPS_DYN` — EVP subcycles per slow (outer) step.
      logical :: a_face_stress = .false.
         !! PR 62: weight the atmospheric stress AND the ice-ocean drag in
         !! the momentum balance by the face ice concentration `a_u`
         !! (`&ocean_ice_nml a_face_stress`). Default off ⇒ byte-identical.
      real(wp) :: cfl_trunc = 0.0_wp
         !! PR 36: SIS2 `CFL_TRUNCATE` (SIS2 default 0.5). Transport-CFL
         !! ceiling on the final ice velocity; `0` disables the clip.
         !! Type-level default is the bit-identity mechanism for every
         !! test call site that does not set it.
      logical :: cfl_trunc_dyn_its = .false.
         !! PR 36: SIS2 `CFL_TRUNC_DYN_ITS` (SIS2 default `.false.`, matches).
         !! Also clip at the bottom of every EVP subcycle.
      logical :: project_ci = .false.
         !! PR 36: SIS2 `PROJECT_ICE_CONCENTRATION` (SIS2 default `.true.`).
         !! Project `ci` forward along the current divergence each subcycle
         !! and recompute `pres_mice` from it.
   end type ice_evp_params_t

contains

   ! =====================================================================
   ! Params constructor
   ! =====================================================================

   pure function ice_evp_params_from_config(p0, c0, ec, cdw, rho_ocean, &
                                            del_sh_min_scale, tdamp, &
                                            evp_sub_steps, a_face_stress, &
                                            cfl_trunc, cfl_trunc_dyn_its, project_ci) result(par)
      !! Small constructor — build once from `&ocean_ice_nml` config.
      !! 11 args (> the style guide's 6): pre-existing deviation, sanctioned
      !! by the derived-type-grouping escape hatch (`FORTRAN_STYLE.md`
      !! §Public procedure arguments) — the whole point of `ice_evp_params_t`
      !! is to be this constructor's one-shot host. Positional (not
      !! `optional`): a knob threaded to the config/schema but not to this
      !! constructor reads from the namelist, validates, and does nothing —
      !! the dead-knob class the audit indicts. The compiler catches the
      !! omission; `optional` would not.
      real(wp), intent(in) :: p0, c0, ec, cdw, rho_ocean, del_sh_min_scale, tdamp
      integer, intent(in) :: evp_sub_steps
      logical, intent(in) :: a_face_stress
      real(wp), intent(in) :: cfl_trunc
      logical, intent(in) :: cfl_trunc_dyn_its, project_ci
      type(ice_evp_params_t) :: par

      par%p0 = p0
      par%c0 = c0
      par%ec = ec
      par%cdw = cdw
      par%rho_ocean = rho_ocean
      par%del_sh_min_scale = del_sh_min_scale
      par%tdamp = tdamp
      par%evp_sub_steps = evp_sub_steps
      par%a_face_stress = a_face_stress
      par%cfl_trunc = cfl_trunc
      par%cfl_trunc_dyn_its = cfl_trunc_dyn_its
      par%project_ci = project_ci
   end function ice_evp_params_from_config

   ! =====================================================================
   ! Driver shim: gather category state, pull ocean velocity, call core
   ! =====================================================================

   subroutine ice_evp_step(grid, metrics, f_corner, ice, ms, dt_slow, par, &
                           periodic_x, periodic_y, dt_transport, n_trunc)
      !! Gathers `mis`/`mice`/`ci` from the category state (mode-branched,
      !! mirrors PR 4b's IST->CAS dispatch), pulls the one-step-lagged
      !! ocean surface velocity, and calls `ice_evp_dynamics` on
      !! `ice%u_ice/v_ice/str_d/str_t/str_s/fxoc/fyoc`.  No-op when the
      !! ice slot is not live or `dynamics` is off (defence-in-depth; the
      !! driver already gates this call on `ice%dynamics`).
      !!
      !! DEVIATION from SPEC §4.1's literal signature: `par` is an
      !! explicit argument here (the spec's shim signature omits it, but
      !! the shim has no other route to `&ocean_ice_nml` — the driver
      !! builds `par` ONCE via `ice_evp_params_from_config` and passes it
      !! into every call, which is both cheaper and clearer than a hidden
      !! module-level singleton).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
         !! Coriolis parameter at corners, shape (nx+1,ny+1) (`coriolis_adv_t%f_corner`).
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt_slow
      type(ice_evp_params_t), intent(in) :: par
      logical, intent(in) :: periodic_x, periodic_y
      real(wp), intent(in), optional :: dt_transport
         !! PR 36: the dt the TRANSPORT step will actually consume
         !! (`ocean_dyn%therm_dt(dt)`), NOT this call's `dt_slow` — EVP
         !! runs every outer step, transport at thermo cadence. Absent =>
         !! `dt_slow` (read only when `par%cfl_trunc > 0`).
      integer, intent(out), optional :: n_trunc
         !! PR 36: count of ice-bearing faces the final CFL clip touched
         !! (`0` when `par%cfl_trunc <= 0`). Mirrors `ice_transport_step`'s
         !! `ok` idiom — the caller (driver) logs, this module does not.

      integer :: nx, ny, nz

      if (.not. ice%is_init .or. .not. ice%dynamics) then
         if (present(n_trunc)) n_trunc = 0
         return
      end if

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! F3: gather into the dedicated INPUT buffers on the workspace slot,
      ! distinct from the `mis_w`/`mice_w`/`ci_w` that `ice_evp_dynamics`
      ! fills from these — no dummy-argument aliasing.
      call ice_cell_concentration_impl(metrics%wet_T, ice%part_size, ice%m_ice, &
                                       ice%m_snow, ice%evp_ws%mis_in_w, ice%evp_ws%mice_in_w, &
                                       ice%evp_ws%ci_in_w, ice%ncat, nx, ny)

      call ice_evp_dynamics(grid, metrics, f_corner, ice%evp_ws%mis_in_w, &
                            ice%evp_ws%mice_in_w, ice%evp_ws%ci_in_w, &
                            ms%u_face_x_layer(:, :, nz), ms%v_face_y_layer(:, :, nz), &
                            ice%tau_a_x, ice%tau_a_y, ice%u_ice, ice%v_ice, &
                            ice%str_d, ice%str_t, ice%str_s, ice%fxoc, ice%fyoc, &
                            dt_slow, par, periodic_x, periodic_y, ice%evp_ws, &
                            dt_transport, n_trunc)
   end subroutine ice_evp_step

   ! =====================================================================
   ! Public core (test seam)
   ! =====================================================================

   subroutine ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, dt_slow, par, periodic_x, periodic_y, ws, &
                               dt_transport, n_trunc)
      !! One outer (slow) EVP call: `evp_sub_steps` subcycles advancing
      !! `ui`/`vi`/`str_d`/`str_t`/`str_s`, plus the subcycle-averaged
      !! ice->ocean stress `fxoc`/`fyoc`.
      !!
      !! Thin public shim (test seam): unpacks the caller-supplied
      !! `evp_workspace_t` slot's components into the explicit-shape
      !! flat-impl args of `ice_evp_dynamics_impl` — the outer-shim +
      !! flat-impl pattern.  NEVER read `ws%component` inside a
      !! `do concurrent` (per-launch descriptor copies); the impl takes the
      !! scratch as explicit-shape `(nx,ny)` dummies instead.  `ws` must be
      !! `init`+`enter_data`'d for this `(grid%nx_total, grid%ny_total)`
      !! (the driver does this in `ocean_sea_ice_t%init`/`enter_data`; tests
      !! build a local `evp_workspace_t`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
      real(wp), intent(in) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), intent(in) :: uo(:, :), vo(:, :)
      real(wp), intent(in) :: tau_ax(:, :), tau_ay(:, :)
      real(wp), intent(inout) :: ui(:, :), vi(:, :)
      real(wp), intent(inout) :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), intent(inout) :: fxoc(:, :), fyoc(:, :)
      real(wp), intent(in) :: dt_slow
      type(ice_evp_params_t), intent(in) :: par
      logical, intent(in) :: periodic_x, periodic_y
      type(evp_workspace_t), intent(inout) :: ws
      real(wp), intent(in), optional :: dt_transport
         !! PR 36: dt TRANSPORT will use for the CFL bound. Absent => `dt_slow`.
      integer, intent(out), optional :: n_trunc
         !! PR 36: count of ice-bearing faces the final clip touched.
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      call ice_evp_dynamics_impl(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                 tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                 fxoc, fyoc, dt_slow, par, periodic_x, periodic_y, nx, ny, &
                                 ws%mis_w, ws%mice_w, ws%ci_w, ws%pres_mice_w, &
                                 ws%del_sh_min_pr_w, ws%sh_dd_w, ws%sh_dt_w, ws%zeta_w, &
                                 ws%del_sh_w, ws%mask_t_w, ws%mi_u_w, ws%mask_u_w, &
                                 ws%u_tmp_w, ws%mi_v_w, ws%mask_v_w, ws%a_u_w, ws%a_v_w, &
                                 ws%sh_ds_w, ws%mi_ratio_a_q_w, ws%q_w, ws%mask_q_w, &
                                 dt_transport, n_trunc)
   end subroutine ice_evp_dynamics

   subroutine ice_evp_dynamics_impl(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                    tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                    fxoc, fyoc, dt_slow, par, periodic_x, periodic_y, nx, ny, &
                                    mis_w, mice_w, ci_w, pres_mice_w, del_sh_min_pr_w, &
                                    sh_dd_w, sh_dt_w, zeta_w, del_sh_w, mask_t_w, &
                                    mi_u_w, mask_u_w, u_tmp_w, mi_v_w, mask_v_w, a_u_w, a_v_w, &
                                    sh_ds_w, mi_ratio_a_q_w, q_w, mask_q_w, &
                                    dt_transport, n_trunc)
      !! Flat-impl core of `ice_evp_dynamics`: the EVP subcycle body with
      !! the persistent scratch passed as EXPLICIT-SHAPE dummies (memory:
      !! never assumed-shape into a `do concurrent` feeder — NVHPC would
      !! emit descriptor-walk memcpys per launch).  The `*_w` scratch names
      !! mirror the retired module workspace 1:1, so the body below is
      !! unchanged from the pre-slot version.
      !!
      !! **Ghost policy is SPLIT — not "all owned here".**  This routine
      !! owns the ghost policy (periodic wrap or zero) for the workspace
      !! copies of `mis`/`mice`/`ci`, for `ui`/`vi`, and for
      !! `str_d`/`str_t`/`str_s` — it wraps those below.  The `intent(in)`
      !! `uo`/`vo`/`tau_ax`/`tau_ay` are consumed DIRECTLY (the momentum
      !! kernel reads their ghost rows, e.g. `uo(i,j)`, `vo(i-1,j+1)`), so
      !! their ghosts are the CALLER's responsibility: a test in a periodic
      !! config must fill those four with their ghosts already wrapped.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
      real(wp), intent(in) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), intent(in) :: uo(:, :), vo(:, :)
      real(wp), intent(in) :: tau_ax(:, :), tau_ay(:, :)
      real(wp), intent(inout) :: ui(:, :), vi(:, :)
      real(wp), intent(inout) :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), intent(inout) :: fxoc(:, :), fyoc(:, :)
      real(wp), intent(in) :: dt_slow
      type(ice_evp_params_t), intent(in) :: par
      logical, intent(in) :: periodic_x, periodic_y
      integer, intent(in) :: nx, ny
         !! T-cell extents (declared before the explicit-shape scratch that
         !! uses them — decl-order rule).
      real(wp), intent(inout) :: mis_w(nx, ny), mice_w(nx, ny), ci_w(nx, ny)
      real(wp), intent(inout) :: pres_mice_w(nx, ny), del_sh_min_pr_w(nx, ny)
      real(wp), intent(inout) :: sh_dd_w(nx, ny), sh_dt_w(nx, ny)
      real(wp), intent(inout) :: zeta_w(nx, ny), del_sh_w(nx, ny)
      real(wp), intent(inout) :: mask_t_w(nx, ny)
      real(wp), intent(inout) :: mi_u_w(nx + 1, ny), mask_u_w(nx + 1, ny), u_tmp_w(nx + 1, ny)
      real(wp), intent(inout) :: mi_v_w(nx, ny + 1), mask_v_w(nx, ny + 1)
      real(wp), intent(inout) :: a_u_w(nx + 1, ny), a_v_w(nx, ny + 1)
         !! PR 62: face ice concentration, valid ONLY when `par%a_face_stress`
         !! (uninitialised device memory otherwise — never read off-gate).
      real(wp), intent(inout) :: sh_ds_w(nx + 1, ny + 1), mi_ratio_a_q_w(nx + 1, ny + 1)
      real(wp), intent(inout) :: q_w(nx + 1, ny + 1), mask_q_w(nx + 1, ny + 1)
      real(wp), intent(in), optional :: dt_transport
         !! PR 36: the dt TRANSPORT will actually consume
         !! (`ocean_dyn%therm_dt(dt)`), NOT this call's `dt_slow` — EVP
         !! runs every outer step, transport at thermo cadence
         !! (`rdb_driver.F90`). Absent => `dt_slow` (SIS2's own assumption:
         !! `SIS_C_dynamics` and `SIS_transport` share `dt_slow`). Read
         !! only when `par%cfl_trunc > 0`.
      integer, intent(out), optional :: n_trunc
         !! PR 36: count of ice-bearing faces (`mi_u`/`mi_v > m_neglect`)
         !! the FINAL clip touched. `0` when `par%cfl_trunc <= 0`.

      integer :: nx_phys, ny_phys, nghost, n
      real(wp) :: dt, tdamp_eff, dt_2tdamp, i_1pdt_t, ec2, i_ec2
      real(wp) :: cdrho, i_cdrhodt, p0_rho, m_neglect, m_neglect2, m_neglect4
      real(wp) :: dt_tr, dt_cum
      logical :: a_face_on, do_trunc_its, do_trunc_fin
      integer :: n_out

      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      ! ---- Scalar precompute (hoisted before the substep loop) ----
      a_face_on = par%a_face_stress
      dt = dt_slow/real(par%evp_sub_steps, wp)
      if (par%tdamp > 0.0_wp) then
         tdamp_eff = par%tdamp
      else if (par%tdamp == 0.0_wp) then
         tdamp_eff = max(0.2_wp*dt_slow, 3.0_wp*dt)
      else
         tdamp_eff = max(-par%tdamp*dt_slow, 3.0_wp*dt)
      end if
      dt_2tdamp = dt/(2.0_wp*tdamp_eff)
      ec2 = par%ec*par%ec
      i_ec2 = 0.0_wp
      if (ec2 > 0.0_wp) i_ec2 = 1.0_wp/ec2
      i_1pdt_t = 1.0_wp/(1.0_wp + dt_2tdamp)
      cdrho = par%cdw*par%rho_ocean
      i_cdrhodt = 1.0_wp/(par%cdw*par%rho_ocean*dt)
      p0_rho = par%p0/ICE_RHO_ICE
      m_neglect = ICE_RHO_ICE*M_NEGLECT_FACTOR
      m_neglect2 = m_neglect*m_neglect
      m_neglect4 = m_neglect2*m_neglect2

      ! ---- PR 36: CFL-truncation gates + the dt the bound must use.
      ! The bound is against the dt TRANSPORT will consume, NOT this call's
      ! `dt_slow` -- Roundabout decouples EVP (every outer step) from transport
      ! (thermo cadence); SIS2's structure assumes they are the same dt.
      ! Absent `dt_transport` => `dt_slow`, matching SIS2's own assumption. ----
      dt_tr = dt_slow
      if (present(dt_transport)) dt_tr = dt_transport
      do_trunc_its = par%cfl_trunc_dyn_its .and. (par%cfl_trunc > 0.0_wp) .and. (dt_tr > 0.0_wp)
      do_trunc_fin = (par%cfl_trunc > 0.0_wp) .and. (dt_tr > 0.0_wp)

      ! ---- Effective masks (SIS2 mask2dT/Cu/Cv/Bu), built once ----
      call evp_build_masks_impl(metrics%wet_T, mask_t_w, mask_u_w, mask_v_w, mask_q_w, &
                                nx_phys, ny_phys, nghost, periodic_x, periodic_y, nx, ny)

      ! ---- Category fields into the workspace, ghost-wrapped/zeroed ----
      call evp_fill_cell_fields_impl(mask_t_w, mis, mice, ci, mis_w, mice_w, ci_w, nx, ny)
      call ocean_periodic_wrap_centre_2d(mis_w, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)
      call ocean_periodic_wrap_centre_2d(mice_w, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)
      call ocean_periodic_wrap_centre_2d(ci_w, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)

      ! ---- Zero ice velocities with no mass (SIS2 :899-907) ----
      call evp_zero_massless_velocity_impl(mask_u_w, mask_v_w, mis_w, ui, vi, nx, ny)

      ! ---- pres_mice + del_sh_min_pr precompute (:877-890) ----
      call evp_pres_mice_impl(metrics%dxT, metrics%dyT, ci_w, p0_rho, par%c0, &
                              par%del_sh_min_scale, tdamp_eff, dt, pres_mice_w, &
                              del_sh_min_pr_w, nx, ny)

      ! ---- mi_u / mi_v (:967-974) ----
      call evp_mi_face_impl(mis_w, mi_u_w, mi_v_w, nx, ny)

      ! ---- PR 62: face ice concentration a_u/a_v, ONLY when a_face_stress.
      ! Reuses evp_mi_face_impl's 0.5-face-average — same expression/edge
      ! convention `ice_ocean_stress_flux_impl`'s `a_u` already uses, so the
      ! momentum budget it closes matches the coupler bit-for-bit (§5.4). ----
      if (a_face_on) call evp_mi_face_impl(ci_w, a_u_w, a_v_w, nx, ny)

      ! ---- q + mi_ratio_A_q (:926-982) ----
      call evp_q_and_mi_ratio_impl(metrics%areaT, f_corner, mask_t_w, mask_u_w, mask_v_w, &
                                   mask_q_w, mis_w, m_neglect, m_neglect2, m_neglect4, &
                                   q_w, mi_ratio_a_q_w, nx, ny)

      ! ---- limit_stresses ONCE before the substep loop (:896) — req (2) ----
      call ice_limit_stresses(metrics%areaT, mask_t_w, pres_mice_w, mice_w, &
                              str_d, str_t, str_s, par%ec, nx, ny)
      call ocean_periodic_wrap_centre_2d(str_d, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)
      call ocean_periodic_wrap_centre_2d(str_t, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)
      call evp_wrap_corner_impl(str_s, nx + 1, ny + 1, nx_phys, ny_phys, nghost, &
                                periodic_x, periodic_y)

      ! ---- Zero the subcycle-averaged ice->ocean stress (F1: an explicit
      ! device kernel, NOT a bare host whole-array assignment — fxoc/fyoc
      ! are copyin-mapped, so a host `= 0.0_wp` would leave the DEVICE copy
      ! stale and each outer call would accumulate onto the prior call's
      ! average). ----
      call evp_zero_stress_impl(fxoc, fyoc, nx, ny)

      ! ---- The EVP subcycle loop — device-resident, no per-substep H<->D ----
      dt_cum = 0.0_wp
      do n = 1, par%evp_sub_steps
         ! PR 36: cumulative elapsed time within THIS call, at the TOP of the
         ! loop (SIS2 :1023,1028) => subcycle n sees t_cum = n*dt. Host
         ! scalar by value -- the substep loop is host-driven, so this costs
         ! nothing and breaks no device residency.
         dt_cum = dt_cum + dt

         call ocean_periodic_wrap_face_x_2d(ui, nx + 1, ny, nx_phys, ny_phys, nghost, &
                                            periodic_x, periodic_y)
         call ocean_periodic_wrap_face_y_2d(vi, nx, ny + 1, nx_phys, ny_phys, nghost, &
                                            periodic_x, periodic_y)

         call evp_sh_ds_impl(metrics%dx_dyBu, metrics%dy_dxBu, metrics%idxCu, metrics%idyCv, &
                             mask_q_w, ui, vi, sh_ds_w, nx, ny)
         call evp_sh_dd_dt_impl(metrics%dy_dxT, metrics%dx_dyT, metrics%iareaT, &
                                metrics%idyCu, metrics%idxCv, metrics%dyCu, metrics%dxCv, &
                                ui, vi, sh_dd_w, sh_dt_w, nx, ny)

         ! ---- PR 36: PROJECT_ICE_CONCENTRATION -- SIS2's position exactly
         ! (:1077 -> :1082): after sh_Dd, before zeta. ci_w is the CALL's
         ! initial (entry-gathered) concentration, never mutated by this --
         ! keep it that way (the projection is `ci*exp(-t_cum*sh_dd)`, not
         ! an incremental accumulator). del_sh_min_pr_w is NOT recomputed
         ! (no ci dependence). ----
         if (par%project_ci) then
            call evp_project_ci_impl(ci_w, sh_dd_w, dt_cum, p0_rho, par%c0, pres_mice_w, nx, ny)
         end if

         call evp_zeta_impl(sh_dd_w, sh_dt_w, sh_ds_w, i_ec2, pres_mice_w, mice_w, &
                            del_sh_min_pr_w, del_sh_w, zeta_w, nx, ny)
         call evp_stress_relax_impl(zeta_w, sh_dd_w, sh_dt_w, pres_mice_w, mice_w, &
                                    i_1pdt_t, dt_2tdamp, i_ec2, str_d, str_t, nx, ny)
         call evp_str_s_relax_impl(metrics%areaT, zeta_w, sh_ds_w, mi_ratio_a_q_w, &
                                   i_1pdt_t, dt_2tdamp, i_ec2, str_s, nx, ny)

         call evp_copy_u_impl(ui, u_tmp_w, nx, ny)

         call evp_u_momentum_impl(metrics%idxCu, metrics%idyCu, metrics%dy2h, &
                                  metrics%dx2q, metrics%iareaCu, mask_u_w, mi_u_w, mi_v_w, &
                                  q_w, str_d, str_t, str_s, uo, vo, tau_ax, ui, vi, &
                                  fxoc, m_neglect, i_cdrhodt, cdrho, dt, nx_phys, ny_phys, &
                                  nghost, nx, ny, a_u_w, a_face_on)
         call evp_v_momentum_impl(metrics%idyCv, metrics%idxCv, metrics%dx2h, &
                                  metrics%dy2q, metrics%iareaCv, mask_v_w, mi_v_w, mi_u_w, &
                                  q_w, str_d, str_t, str_s, uo, vo, tau_ay, u_tmp_w, vi, &
                                  fyoc, m_neglect, i_cdrhodt, cdrho, dt, nx_phys, ny_phys, &
                                  nghost, nx, ny, a_v_w, a_face_on)

         ! ---- PR 36: in-loop CFL clip (cfl_trunc_dyn_its) -- SIS2's
         ! position (:1338, bottom of the loop, after both momentum
         ! solves). Exact bound, no count. No re-wrap needed: the next
         ! iteration's wrap at the top of the loop does it. ----
         if (do_trunc_its) then
            call evp_truncate_velocity_impl(metrics%areaT, metrics%dy_cu, metrics%dx_cv, &
                                            ui, vi, par%cfl_trunc, dt_tr, 1.0_wp, &
                                            nghost, nx_phys, ny_phys, nx, ny)
         end if
      end do

      ! ---- PR 36: FINAL CFL clip (cfl_trunc) -- always on when cfl_trunc>0
      ! and dt_tr>0. 0.95*bound back-off; counts ice-bearing faces touched.
      ! The clip walks PHYSICAL faces only -- the periodic re-wrap after it
      ! is MANDATORY (without it a periodic ghost face keeps its unclipped
      ! value, transport's PPM reads it as a donor velocity, and the run
      ! aborts anyway for a reason no test would name). ----
      n_out = 0
      if (do_trunc_fin) then
         call evp_truncate_final_impl(metrics%areaT, metrics%dy_cu, metrics%dx_cv, &
                                      mi_u_w, mi_v_w, ui, vi, par%cfl_trunc, dt_tr, &
                                      m_neglect, nghost, nx_phys, ny_phys, nx, ny, n_out)
         call ocean_periodic_wrap_face_x_2d(ui, nx + 1, ny, nx_phys, ny_phys, nghost, &
                                            periodic_x, periodic_y)
         call ocean_periodic_wrap_face_y_2d(vi, nx, ny + 1, nx_phys, ny_phys, nghost, &
                                            periodic_x, periodic_y)
      end if
      if (present(n_trunc)) n_trunc = n_out

      ! ---- fxoc/fyoc average + mask (:1415-1441) ----
      call evp_average_stress_impl(mask_u_w, mask_v_w, fxoc, fyoc, par%evp_sub_steps, nx, ny)
   end subroutine ice_evp_dynamics_impl

   ! =====================================================================
   ! mi_ratio_A_q — full SIS2 harmonic-mean form (requirement 5)
   ! =====================================================================

   pure function ice_evp_mi_ratio_point(mis_sw, mis_se, mis_nw, mis_ne, &
                                        mask_u_below, mask_u_above, &
                                        mask_v_left, mask_v_right, &
                                        mask_q, area_sw, area_se, area_nw, area_ne, &
                                        mask_t_sw, mask_t_se, mask_t_nw, mask_t_ne, &
                                        m_neglect2, m_neglect4) result(mi_ratio)
      !! `mi_ratio_A_q` at a single corner (SIS2 :926-964), FULL form —
      !! all four branches (interior / corner-coast / straight-coast /
      !! land). `weak_coast_stress=.false.` hardwired (SIS2 default):
      !! `sum_area` is the MASKED area sum of the 4 surrounding T-cells.
      !! Factored out of the fill kernel so a unit test can pin it
      !! directly (SPEC §7 gate 8).
      real(wp), intent(in) :: mis_sw, mis_se, mis_nw, mis_ne
         !! Ice+snow mass per cell area at the 4 T-cells around the
         !! corner (SW/SE/NW/NE, rdb `wet_q` convention).
      real(wp), intent(in) :: mask_u_below, mask_u_above
         !! u-face masks below/above the corner (SIS2
         !! `mask2dCu(I,j)`/`mask2dCu(I,j+1)`).
      real(wp), intent(in) :: mask_v_left, mask_v_right
         !! v-face masks left/right of the corner (SIS2
         !! `mask2dCv(i,J)`/`mask2dCv(i+1,J)`).
      real(wp), intent(in) :: mask_q
         !! `mask2dBu` at this corner (1 = genuinely interior ocean point).
      real(wp), intent(in) :: area_sw, area_se, area_nw, area_ne
         !! T-cell areas at the 4 surrounding cells.
      real(wp), intent(in) :: mask_t_sw, mask_t_se, mask_t_nw, mask_t_ne
         !! T-cell wet masks at the 4 surrounding cells (land => 0).
      real(wp), intent(in) :: m_neglect2, m_neglect4
      real(wp) :: mi_ratio

      real(wp) :: sum_area, muq2, mvq2, muq, mvq

      sum_area = (mask_t_sw*area_sw + mask_t_ne*area_ne) + &
                 (mask_t_nw*area_nw + mask_t_se*area_se)

      if (sum_area <= 0.0_wp) then
         mi_ratio = 0.0_wp
      else if (mask_q > 0.0_wp) then
         muq2 = 0.25_wp*(mis_sw + mis_se)*(mis_nw + mis_ne)
         mvq2 = 0.25_wp*(mis_sw + mis_nw)*(mis_se + mis_ne)
         mi_ratio = 32.0_wp*muq2*mvq2/((m_neglect4 + (muq2 + mvq2)* &
                                        ((mis_sw + mis_ne) + (mis_nw + mis_se))**2)*sum_area)
      else if ((mask_u_below + mask_u_above) + (mask_v_left + mask_v_right) > 1.5_wp) then
         muq = 0.5_wp*(mask_u_below*(mis_sw + mis_se) + mask_u_above*(mis_nw + mis_ne))
         mvq = 0.5_wp*(mask_v_left*(mis_sw + mis_nw) + mask_v_right*(mis_se + mis_ne))
         mi_ratio = 4.0_wp*muq*mvq/((m_neglect2 + (muq + mvq)**2)*sum_area)
      else
         mi_ratio = 1.0_wp/sum_area
      end if
   end function ice_evp_mi_ratio_point

   ! =====================================================================
   ! Flat-impl kernels
   ! =====================================================================

   pure subroutine evp_build_masks_impl(wet_t, mask_t, mask_u, mask_v, mask_q, &
                                        nx_phys, ny_phys, nghost, periodic_x, periodic_y, &
                                        nx, ny)
      !! `mask_t`: `wet_T` inside the physical domain; ghosts = periodic
      !! wrap or 0 (SIS2 `mask2dT` semantics — pins non-periodic ghosts
      !! to land, matching SIS2's own domain-edge convention even when a
      !! driver run's `wet_T` ghost happens to read 1).
      !! `mask_u(i,j) = mask_t(i-1,j)*mask_t(i,j)`, `mask_v` ditto in y,
      !! `mask_q` = product of the 4 surrounding `mask_t` (SIS2
      !! `mask2dBu`).
      integer, intent(in) :: nx_phys, ny_phys, nghost, nx, ny
      real(wp), intent(in) :: wet_t(nx, ny)
      logical, intent(in) :: periodic_x, periodic_y
      real(wp), intent(out) :: mask_t(nx, ny)
      real(wp), intent(out) :: mask_u(nx + 1, ny)
      real(wp), intent(out) :: mask_v(nx, ny + 1)
      real(wp), intent(out) :: mask_q(nx + 1, ny + 1)
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: mt_sw, mt_se, mt_nw, mt_ne

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      do concurrent(j=1:ny, i=1:nx)
         if (i >= i_lo .and. i <= i_hi .and. j >= j_lo .and. j <= j_hi) then
            mask_t(i, j) = merge(1.0_wp, 0.0_wp, wet_t(i, j) > 0.5_wp)
         else
            mask_t(i, j) = 0.0_wp
         end if
      end do
      call ocean_periodic_wrap_centre_2d(mask_t, nx, ny, nx_phys, ny_phys, nghost, &
                                         periodic_x, periodic_y)

      do concurrent(j=1:ny, i=1:nx + 1)
         if (i == 1 .or. i == nx + 1) then
            mask_u(i, j) = 0.0_wp
         else
            mask_u(i, j) = mask_t(i - 1, j)*mask_t(i, j)
         end if
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         if (j == 1 .or. j == ny + 1) then
            mask_v(i, j) = 0.0_wp
         else
            mask_v(i, j) = mask_t(i, j - 1)*mask_t(i, j)
         end if
      end do
      do concurrent(j=1:ny + 1, i=1:nx + 1) local(mt_sw, mt_se, mt_nw, mt_ne)
         if (i == 1 .or. i == nx + 1 .or. j == 1 .or. j == ny + 1) then
            mask_q(i, j) = 0.0_wp
         else
            mt_sw = mask_t(i - 1, j - 1)
            mt_se = mask_t(i, j - 1)
            mt_nw = mask_t(i - 1, j)
            mt_ne = mask_t(i, j)
            mask_q(i, j) = mt_sw*mt_se*mt_nw*mt_ne
         end if
      end do
   end subroutine evp_build_masks_impl

   pure subroutine evp_fill_cell_fields_impl(mask_t, mis_in, mice_in, ci_in, &
                                             mis_out, mice_out, ci_out, nx, ny)
      !! Interior copy of the gathered `mis`/`mice`/`ci`, masked to `mask_t`
      !! (defence-in-depth beyond the caller's own `wet_T` gate); ghost
      !! rows/cols zeroed (the periodic wrap that follows fills them).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: mask_t(nx, ny)
      real(wp), intent(in) :: mis_in(nx, ny), mice_in(nx, ny), ci_in(nx, ny)
      real(wp), intent(out) :: mis_out(nx, ny), mice_out(nx, ny), ci_out(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         mis_out(i, j) = mask_t(i, j)*mis_in(i, j)
         mice_out(i, j) = mask_t(i, j)*mice_in(i, j)
         ci_out(i, j) = mask_t(i, j)*ci_in(i, j)
      end do
   end subroutine evp_fill_cell_fields_impl

   pure subroutine evp_zero_massless_velocity_impl(mask_u, mask_v, mis, ui, vi, nx, ny)
      !! SIS2 :899-907 — zero ice velocities where BOTH neighbouring
      !! cells are massless (or the face is masked/land).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: mask_u(nx + 1, ny)
      real(wp), intent(in) :: mask_v(nx, ny + 1)
      real(wp), intent(in) :: mis(nx, ny)
      real(wp), intent(inout) :: ui(nx + 1, ny)
      real(wp), intent(inout) :: vi(nx, ny + 1)
      integer :: i, j
      real(wp) :: mleft, mright

      ! F2: explicit `if` branches — `merge(mis(i-1,j), 0, i>1)` would
      ! still EVALUATE the OOB `mis(0,j)` reference at i=1 (merge does not
      ! conditionally evaluate its args), which is UB and traps under
      ! -Mbounds. Guard the read itself (mirrors `evp_mi_face_impl`).
      do concurrent(j=1:ny, i=1:nx + 1) local(mleft, mright)
         if (i > 1) then
            mleft = mis(i - 1, j)
         else
            mleft = 0.0_wp
         end if
         if (i <= nx) then
            mright = mis(i, j)
         else
            mright = 0.0_wp
         end if
         if (mask_u(i, j)*(mleft + mright) == 0.0_wp) ui(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx) local(mleft, mright)
         if (j > 1) then
            mleft = mis(i, j - 1)
         else
            mleft = 0.0_wp
         end if
         if (j <= ny) then
            mright = mis(i, j)
         else
            mright = 0.0_wp
         end if
         if (mask_v(i, j)*(mleft + mright) == 0.0_wp) vi(i, j) = 0.0_wp
      end do
   end subroutine evp_zero_massless_velocity_impl

   pure subroutine evp_pres_mice_impl(dxT, dyT, ci, p0_rho, c0, del_sh_min_scale, &
                                      tdamp_eff, dt, pres_mice, del_sh_min_pr, nx, ny)
      !! `pres_mice = p0_rho*exp(-c0*max(1-ci,0))` (:878); `dxharm =
      !! 2*dxT*dyT/(dxT+dyT)`; `del_sh_min_pr = 2*del_sh_min_scale*dt^2 /
      !! (Tdamp*dxharm^2)` guarded on `dxharm > 0` (:880-890).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), ci(nx, ny)
      real(wp), intent(in) :: p0_rho, c0, del_sh_min_scale, tdamp_eff, dt
      real(wp), intent(out) :: pres_mice(nx, ny), del_sh_min_pr(nx, ny)
      integer :: i, j
      real(wp) :: dxharm

      do concurrent(j=1:ny, i=1:nx) local(dxharm)
         pres_mice(i, j) = p0_rho*exp(-c0*max(1.0_wp - ci(i, j), 0.0_wp))
         dxharm = 2.0_wp*dxT(i, j)*dyT(i, j)/(dxT(i, j) + dyT(i, j))
         if (dxharm > 0.0_wp) then
            del_sh_min_pr(i, j) = (2.0_wp*del_sh_min_scale*dt**2)/(tdamp_eff*dxharm**2)
         else
            del_sh_min_pr(i, j) = 0.0_wp
         end if
      end do
   end subroutine evp_pres_mice_impl

   pure subroutine evp_project_ci_impl(ci, sh_dd, dt_cum, p0_rho, c0, pres_mice, nx, ny)
      !! PR 36: `PROJECT_ICE_CONCENTRATION` (SIS2 `SIS_dyn_cgrid.F90:1064-
      !! 1077`). `ci_proj = ci*exp(-dt_cum*sh_dd)` then `pres_mice =
      !! p0_rho*exp(-c0*max(1-ci_proj, 0))`. `ci_proj` is a `local()`
      !! scalar, NOT an array: SIS2 materialises it only for the
      !! sigI/sigII/find_ice_strength diagnostics Roundabout does not have
      !! (documented divergence). `del_sh_min_pr` is NOT recomputed here
      !! (it has no `ci` dependence, `evp_pres_mice_impl` above). `ci_proj`
      !! is deliberately unclamped above 1 -- `max(1-ci_proj, 0)` already
      !! saturates the effect at `p0_rho`, and for `dt_cum*|sh_dd| > 709`
      !! (an unreachable regime in any sane run) `exp` overflows to `+Inf`,
      !! `max(1-Inf, 0) = 0`, `exp(0) = 1` -- IEEE launders the overflow to
      !! exactly the correct saturated value, so no guard is needed (SIS2
      !! has none either).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: ci(nx, ny), sh_dd(nx, ny)
      real(wp), intent(in) :: dt_cum, p0_rho, c0
      real(wp), intent(inout) :: pres_mice(nx, ny)
      integer :: i, j
      real(wp) :: ci_proj

      do concurrent(j=1:ny, i=1:nx) local(ci_proj)
         ci_proj = ci(i, j)*exp(-dt_cum*sh_dd(i, j))
         pres_mice(i, j) = p0_rho*exp(-c0*max(1.0_wp - ci_proj, 0.0_wp))
      end do
   end subroutine evp_project_ci_impl

   pure subroutine evp_truncate_velocity_impl(areaT, dy_cu, dx_cv, ui, vi, cfl_trunc, dt_tr, &
                                              backoff, nghost, nx_phys, ny_phys, nx, ny)
      !! PR 36: the shared CFL-clip algebra -- the transport-CFL bound on
      !! the ice velocity (SIS2 `SIS_dyn_cgrid.F90:839-870`, the in-loop
      !! half at `:1338-1361`, the final half at `:1443-1500`; this
      !! routine is the counting-free, caller-chosen-backoff form both
      !! reuse; `evp_truncate_final_impl` below wraps it with the 0.95
      !! back-off and the `mi > m_neglect` count).
      !!
      !! `u_max(face) = +cfl_trunc*areaT(donor for u>0)/(dt_tr*dy_cu(face))`,
      !! `u_min(face) = -cfl_trunc*areaT(donor for u<0)/(dt_tr*dy_cu(face))`
      !! -- "the flux out of a cell in one slow step cannot exceed
      !! cfl_trunc of its volume". The donor asymmetry is load-bearing:
      !! `+u` at rdb u-face `(i,j)` (the WEST face of cell `(i,j)`,
      !! module docstring §1) drains the WEST cell `(i-1,j)`; `-u` drains
      !! the EAST cell `(i,j)`. v-mirror: `+v` drains the SOUTH cell
      !! `(i,j-1)`, `-v` drains the NORTH cell `(i,j)`.
      !!
      !! `dy_cu`/`dx_cv` (NOT the unmasked `dyCu`/`dxCv`) are the
      !! topography-aware OPEN face widths -- zero at a closed/land face,
      !! which is why the bound is guarded `> 0.0`: a closed face gets
      !! `u_hi = u_lo = 0` (forces `ui = 0` there), not a finite spurious
      !! bound from dividing by a nonzero length at land.
      !!
      !! Loop ranges are copied VERBATIM from `evp_u_momentum_impl` (u) and
      !! `evp_v_momentum_impl` (v) -- physical faces only. Over that range
      !! `i-1 >= nghost >= 1` (resp. `j-1 >= nghost >= 1`) always, so no
      !! array-edge branch is needed (unlike `evp_mi_face_impl`, which
      !! loops the full `1:nx+1`/`1:ny+1` and does need one).
      integer, intent(in) :: nghost, nx_phys, ny_phys, nx, ny
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp), intent(inout) :: ui(nx + 1, ny), vi(nx, ny + 1)
      real(wp), intent(in) :: cfl_trunc, dt_tr, backoff
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: u_hi, u_lo, v_hi, v_lo, loc_scale

      i_lo = nghost + 1
      i_hi = nghost + nx_phys + 1
      j_lo = nghost + 1
      j_hi = nghost + ny_phys
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(u_hi, u_lo, loc_scale)
         u_hi = 0.0_wp
         u_lo = 0.0_wp
         if (dy_cu(i, j) > 0.0_wp) then
            loc_scale = cfl_trunc/(dt_tr*dy_cu(i, j))
            u_hi = backoff*loc_scale*areaT(i - 1, j)
            u_lo = -backoff*loc_scale*areaT(i, j)
         end if
         if (ui(i, j) > u_hi) then
            ui(i, j) = u_hi
         else if (ui(i, j) < u_lo) then
            ui(i, j) = u_lo
         end if
      end do

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys + 1
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(v_hi, v_lo, loc_scale)
         v_hi = 0.0_wp
         v_lo = 0.0_wp
         if (dx_cv(i, j) > 0.0_wp) then
            loc_scale = cfl_trunc/(dt_tr*dx_cv(i, j))
            v_hi = backoff*loc_scale*areaT(i, j - 1)
            v_lo = -backoff*loc_scale*areaT(i, j)
         end if
         if (vi(i, j) > v_hi) then
            vi(i, j) = v_hi
         else if (vi(i, j) < v_lo) then
            vi(i, j) = v_lo
         end if
      end do
   end subroutine evp_truncate_velocity_impl

   pure subroutine evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                           cfl_trunc, dt_tr, m_neglect, nghost, nx_phys, &
                                           ny_phys, nx, ny, n_trunc)
      !! PR 36: the FINAL CFL clip (SIS2 `:1443-1500`) -- `TRUNC_BACKOFF`
      !! (0.95) back-off instead of the exact bound, PLUS a count of the
      !! ice-bearing faces it touched (`mi > m_neglect`, SIS2 `:1466,1469`
      !! -- massless faces clip silently, matching SIS2: counting them
      !! would flood the driver's warning with meaningless ice-free clips
      !! at every margin). Not a `do concurrent`: reductions use
      !! `!$acc parallel loop reduction(...)` (`ice_compress_impl` is the
      !! local precedent for a reduction that also mutates the arrays it
      !! walks). Same bound algebra as `evp_truncate_velocity_impl`,
      !! duplicated rather than shared: the in-loop variant runs
      !! `evp_sub_steps` (432 by default) times per outer step and must
      !! NOT carry a reduction (each would be a device->host sync); this
      !! variant runs once and must. `CLAUDE.md`'s "duplicate explicitly"
      !! rule -- merging the two costs 432 syncs per outer step.
      integer, intent(in) :: nghost, nx_phys, ny_phys, nx, ny
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp), intent(in) :: mi_u(nx + 1, ny), mi_v(nx, ny + 1)
      real(wp), intent(inout) :: ui(nx + 1, ny), vi(nx, ny + 1)
      real(wp), intent(in) :: cfl_trunc, dt_tr, m_neglect
      integer, intent(out) :: n_trunc
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: u_hi, u_lo, v_hi, v_lo, loc_scale
      integer :: n_acc

      n_acc = 0

      i_lo = nghost + 1
      i_hi = nghost + nx_phys + 1
      j_lo = nghost + 1
      j_hi = nghost + ny_phys
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(u_hi, u_lo, loc_scale) reduce(+:n_acc)
         u_hi = 0.0_wp
         u_lo = 0.0_wp
         if (dy_cu(i, j) > 0.0_wp) then
            loc_scale = cfl_trunc/(dt_tr*dy_cu(i, j))
            u_hi = TRUNC_BACKOFF*loc_scale*areaT(i - 1, j)
            u_lo = -TRUNC_BACKOFF*loc_scale*areaT(i, j)
         end if
         if (ui(i, j) > u_hi .or. ui(i, j) < u_lo) then
            if (mi_u(i, j) > m_neglect) n_acc = n_acc + 1
            ui(i, j) = merge(u_hi, u_lo, ui(i, j) > u_hi)
         end if
      end do

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys + 1
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(v_hi, v_lo, loc_scale) reduce(+:n_acc)
         v_hi = 0.0_wp
         v_lo = 0.0_wp
         if (dx_cv(i, j) > 0.0_wp) then
            loc_scale = cfl_trunc/(dt_tr*dx_cv(i, j))
            v_hi = TRUNC_BACKOFF*loc_scale*areaT(i, j - 1)
            v_lo = -TRUNC_BACKOFF*loc_scale*areaT(i, j)
         end if
         if (vi(i, j) > v_hi .or. vi(i, j) < v_lo) then
            if (mi_v(i, j) > m_neglect) n_acc = n_acc + 1
            vi(i, j) = merge(v_hi, v_lo, vi(i, j) > v_hi)
         end if
      end do

      n_trunc = n_acc
   end subroutine evp_truncate_final_impl

   pure subroutine evp_mi_face_impl(mis, mi_u, mi_v, nx, ny)
      !! `mi_u(i,j) = 0.5*(mis(i-1,j)+mis(i,j))`; `mi_v(i,j) =
      !! 0.5*(mis(i,j-1)+mis(i,j))` (SIS2 :967-974, rdb index
      !! translation §1). Array-edge faces (`i=1`/`i=nx+1`, `j=1`/
      !! `j=ny+1`) have no neighbour on one side; `mis` at those ghost
      !! rows/cols was already periodic-wrapped or zeroed, so a naive
      !! `mis(i-1,j)`/`mis(i,j)` read is always in-bounds here EXCEPT at
      !! the two hard array edges themselves — those faces are handled
      !! explicitly.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: mis(nx, ny)
      real(wp), intent(out) :: mi_u(nx + 1, ny)
      real(wp), intent(out) :: mi_v(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         if (i == 1) then
            mi_u(i, j) = 0.5_wp*mis(1, j)
         else if (i == nx + 1) then
            mi_u(i, j) = 0.5_wp*mis(nx, j)
         else
            mi_u(i, j) = 0.5_wp*(mis(i - 1, j) + mis(i, j))
         end if
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         if (j == 1) then
            mi_v(i, j) = 0.5_wp*mis(i, 1)
         else if (j == ny + 1) then
            mi_v(i, j) = 0.5_wp*mis(i, ny)
         else
            mi_v(i, j) = 0.5_wp*(mis(i, j - 1) + mis(i, j))
         end if
      end do
   end subroutine evp_mi_face_impl

   pure subroutine evp_q_and_mi_ratio_impl(areaT, f_corner, mask_t, mask_u, mask_v, mask_q, &
                                           mis, m_neglect, m_neglect2, m_neglect4, &
                                           q, mi_ratio_a_q, nx, ny)
      !! `q(ic,jc) = f_corner*tot_area / (Σ areaT*mis over the 4 cells +
      !! tot_area*m_neglect)` (:977-982); `mi_ratio_A_q` via
      !! `ice_evp_mi_ratio_point` (requirement 5).  4 T-cells around
      !! corner `(ic,jc)`: `(ic-1,jc-1) (ic,jc-1) (ic-1,jc) (ic,jc)`
      !! (SW/SE/NW/NE, §1). Array-edge corners (no T-cell on one side)
      !! get `q=0`/`mi_ratio=0` (land-corner convention — consistent with
      !! `mask_t=0` beyond the array edge).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1)
      real(wp), intent(in) :: mask_t(nx, ny)
      real(wp), intent(in) :: mask_u(nx + 1, ny)
      real(wp), intent(in) :: mask_v(nx, ny + 1)
      real(wp), intent(in) :: mask_q(nx + 1, ny + 1)
      real(wp), intent(in) :: mis(nx, ny)
      real(wp), intent(in) :: m_neglect, m_neglect2, m_neglect4
      real(wp), intent(out) :: q(nx + 1, ny + 1)
      real(wp), intent(out) :: mi_ratio_a_q(nx + 1, ny + 1)
      integer :: ic, jc
      real(wp) :: tot_area, mass_sum
      real(wp) :: a_sw, a_se, a_nw, a_ne
      real(wp) :: m_sw, m_se, m_nw, m_ne
      real(wp) :: mt_sw, mt_se, mt_nw, mt_ne

      do concurrent(jc=1:ny + 1, ic=1:nx + 1) &
         local(tot_area, mass_sum, a_sw, a_se, a_nw, a_ne, m_sw, m_se, m_nw, m_ne, &
               mt_sw, mt_se, mt_nw, mt_ne)
         if (ic == 1 .or. ic == nx + 1 .or. jc == 1 .or. jc == ny + 1) then
            q(ic, jc) = 0.0_wp
            mi_ratio_a_q(ic, jc) = 0.0_wp
         else
            a_sw = areaT(ic - 1, jc - 1)
            a_se = areaT(ic, jc - 1)
            a_nw = areaT(ic - 1, jc)
            a_ne = areaT(ic, jc)
            m_sw = mis(ic - 1, jc - 1)
            m_se = mis(ic, jc - 1)
            m_nw = mis(ic - 1, jc)
            m_ne = mis(ic, jc)
            mt_sw = mask_t(ic - 1, jc - 1)
            mt_se = mask_t(ic, jc - 1)
            mt_nw = mask_t(ic - 1, jc)
            mt_ne = mask_t(ic, jc)

            tot_area = (a_sw + a_ne) + (a_nw + a_se)
            mass_sum = (a_sw*m_sw + a_ne*m_ne) + (a_nw*m_nw + a_se*m_se)
            q(ic, jc) = f_corner(ic, jc)*tot_area/(mass_sum + tot_area*m_neglect)

            mi_ratio_a_q(ic, jc) = ice_evp_mi_ratio_point( &
                                   m_sw, m_se, m_nw, m_ne, &
                                   mask_u(ic, jc - 1), mask_u(ic, jc), &
                                   mask_v(ic - 1, jc), mask_v(ic, jc), &
                                   mask_q(ic, jc), a_sw, a_se, a_nw, a_ne, &
                                   mt_sw, mt_se, mt_nw, mt_ne, m_neglect2, m_neglect4)
         end if
      end do
   end subroutine evp_q_and_mi_ratio_impl

   ! ---------------------------------------------------------------------
   ! limit_stresses (requirements 2 + 3)
   ! ---------------------------------------------------------------------

   pure subroutine ice_limit_stresses(areaT, mask_t, pres_mice, mice, str_d, str_t, str_s, &
                                      ec, nx, ny)
      !! SIS2 `limit_stresses` (:1619-1684), `lim=1` (no optional arg).
      !! Called ONCE per `ice_evp_dynamics` call, BEFORE the substep loop
      !! — requirement (2). Corner clamp uses the MASKED-area-weighted
      !! mean pressure of the <=4 wet neighbours — requirement (3).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: mask_t(nx, ny)
      real(wp), intent(in) :: pres_mice(nx, ny)
      real(wp), intent(in) :: mice(nx, ny)
      real(wp), intent(inout) :: str_d(nx, ny)
      real(wp), intent(inout) :: str_t(nx, ny)
      real(wp), intent(inout) :: str_s(nx + 1, ny + 1)
      real(wp), intent(in) :: ec
      integer :: i, j, ic, jc
      real(wp) :: pressure, i_2ec, lim_2
      real(wp) :: sum_area, pres_avg

      i_2ec = 0.0_wp
      if (ec > 0.0_wp) i_2ec = 0.5_wp/ec
      lim_2 = 0.5_wp

      do concurrent(j=1:ny, i=1:nx) local(pressure)
         pressure = pres_mice(i, j)*mice(i, j)
         if (str_d(i, j) < -pressure) str_d(i, j) = -pressure
         if (ec*str_t(i, j) > lim_2*pressure) str_t(i, j) = i_2ec*pressure
         if (ec*str_t(i, j) < -lim_2*pressure) str_t(i, j) = -i_2ec*pressure
      end do

      do concurrent(jc=1:ny + 1, ic=1:nx + 1) local(sum_area, pres_avg)
         if (ic == 1 .or. ic == nx + 1 .or. jc == 1 .or. jc == ny + 1) then
            ! Array-edge corner: no 4th neighbour exists; leave str_s
            ! untouched (these are always ghost/land corners under the
            ! periodic-or-wall ghost policy — never read by the momentum
            ! solve at a physical interior face).
            continue
         else
            sum_area = (mask_t(ic - 1, jc - 1)*areaT(ic - 1, jc - 1) + &
                        mask_t(ic, jc)*areaT(ic, jc)) + &
                       (mask_t(ic - 1, jc)*areaT(ic - 1, jc) + &
                        mask_t(ic, jc - 1)*areaT(ic, jc - 1))
            pres_avg = 0.0_wp
            if (sum_area > 0.0_wp) then
               pres_avg = ((mask_t(ic - 1, jc - 1)*areaT(ic - 1, jc - 1)* &
                            (pres_mice(ic - 1, jc - 1)*mice(ic - 1, jc - 1)) + &
                            mask_t(ic, jc)*areaT(ic, jc)* &
                            (pres_mice(ic, jc)*mice(ic, jc))) + &
                           (mask_t(ic - 1, jc)*areaT(ic - 1, jc)* &
                            (pres_mice(ic - 1, jc)*mice(ic - 1, jc)) + &
                            mask_t(ic, jc - 1)*areaT(ic, jc - 1)* &
                            (pres_mice(ic, jc - 1)*mice(ic, jc - 1))))/sum_area
            end if
            if (ec*str_s(ic, jc) > lim_2*pres_avg) str_s(ic, jc) = i_2ec*pres_avg
            if (ec*str_s(ic, jc) < -lim_2*pres_avg) str_s(ic, jc) = -i_2ec*pres_avg
         end if
      end do
   end subroutine ice_limit_stresses

   ! ---------------------------------------------------------------------
   ! Subcycle-loop kernels (SIS2 :1026-1406)
   ! ---------------------------------------------------------------------

   pure subroutine evp_wrap_corner_impl(fld, nx_face, ny_face, nx_phys, ny_phys, nghost, &
                                        wrap_x, wrap_y)
      !! Periodic ghost-wrap for a corner-staggered field (e.g. `str_s`),
      !! shape (nx_total+1, ny_total+1). No corner-wrap helper exists in
      !! `rdb_ocean_periodic` (only centre/face_x/face_y) — this is the
      !! EVP-local twin, same two-pass (x-then-y) structure.
      integer, intent(in) :: nx_face, ny_face, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_face, ny_face)
      logical, intent(in) :: wrap_x, wrap_y
      integer :: i, j
      integer :: i_w, i_e, j_s, j_n

      i_w = nghost + 1
      i_e = nghost + nx_phys + 1
      j_s = nghost + 1
      j_n = nghost + ny_phys + 1

      if (wrap_x) then
         do concurrent(j=1:ny_face, i=1:nx_face)
            if (i <= nghost) then
               fld(i, j) = fld(i + nx_phys, j)
            end if
            if (i > nx_phys + nghost + 1) then
               fld(i, j) = fld(i - nx_phys, j)
            end if
            if (i == i_e) fld(i, j) = fld(i_w, j)
         end do
      end if
      if (wrap_y) then
         do concurrent(j=1:ny_face, i=1:nx_face)
            if (j <= nghost) then
               fld(i, j) = fld(i, j + ny_phys)
            end if
            if (j > ny_phys + nghost + 1) then
               fld(i, j) = fld(i, j - ny_phys)
            end if
            if (j == j_n) fld(i, j) = fld(i, j_s)
         end do
      end if
   end subroutine evp_wrap_corner_impl

   pure subroutine evp_sh_ds_impl(dx_dyBu, dy_dxBu, idxCu, idyCv, mask_q, ui, vi, sh_ds, nx, ny)
      !! sh_Ds at corners (:1045-1050) — requirement (4): the SINGLE
      !! scalar no-slip factor `(2-mask_q)` on the WHOLE combined strain.
      !! Computed over the interior+1 ring (ic,jc in [1,nx+1]x[1,ny+1] —
      !! the full corner array; out-of-band neighbours contribute 0 via
      !! zero ghost velocities at the hard array edges, never per-term
      !! mirroring).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dx_dyBu(nx + 1, ny + 1), dy_dxBu(nx + 1, ny + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: mask_q(nx + 1, ny + 1)
      real(wp), intent(in) :: ui(nx + 1, ny), vi(nx, ny + 1)
      real(wp), intent(out) :: sh_ds(nx + 1, ny + 1)
      integer :: ic, jc
      real(wp) :: du_term, dv_term

      do concurrent(jc=1:ny + 1, ic=1:nx + 1) local(du_term, dv_term)
         du_term = 0.0_wp
         if (jc <= ny .and. jc >= 1) then
            du_term = ui(ic, jc)*idxCu(ic, jc)
         end if
         if (jc - 1 >= 1 .and. jc - 1 <= ny) then
            du_term = du_term - ui(ic, jc - 1)*idxCu(ic, jc - 1)
         end if
         dv_term = 0.0_wp
         if (ic <= nx .and. ic >= 1) then
            dv_term = vi(ic, jc)*idyCv(ic, jc)
         end if
         if (ic - 1 >= 1 .and. ic - 1 <= nx) then
            dv_term = dv_term - vi(ic - 1, jc)*idyCv(ic - 1, jc)
         end if
         sh_ds(ic, jc) = (2.0_wp - mask_q(ic, jc))*(dx_dyBu(ic, jc)*du_term + &
                                                    dy_dxBu(ic, jc)*dv_term)
      end do
   end subroutine evp_sh_ds_impl

   pure subroutine evp_sh_dd_dt_impl(dy_dxT, dx_dyT, iareaT, idyCu, idxCv, dyCu, dxCv, &
                                     ui, vi, sh_dd, sh_dt, nx, ny)
      !! sh_Dt / sh_Dd at cells (:1053-1061).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dy_dxT(nx, ny), dx_dyT(nx, ny), iareaT(nx, ny)
      real(wp), intent(in) :: idyCu(nx + 1, ny), idxCv(nx, ny + 1)
      real(wp), intent(in) :: dyCu(nx + 1, ny), dxCv(nx, ny + 1)
      real(wp), intent(in) :: ui(nx + 1, ny), vi(nx, ny + 1)
      real(wp), intent(out) :: sh_dd(nx, ny), sh_dt(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         sh_dt(i, j) = dy_dxT(i, j)*(idyCu(i + 1, j)*ui(i + 1, j) - idyCu(i, j)*ui(i, j)) - &
                       dx_dyT(i, j)*(idxCv(i, j + 1)*vi(i, j + 1) - idxCv(i, j)*vi(i, j))
         sh_dd(i, j) = iareaT(i, j)*((dyCu(i + 1, j)*ui(i + 1, j) - dyCu(i, j)*ui(i, j)) + &
                                     (dxCv(i, j + 1)*vi(i, j + 1) - dxCv(i, j)*vi(i, j)))
      end do
   end subroutine evp_sh_dd_dt_impl

   pure subroutine evp_zeta_impl(sh_dd, sh_dt, sh_ds, i_ec2, pres_mice, mice, &
                                 del_sh_min_pr, del_sh, zeta, nx, ny)
      !! del_sh / zeta (:1082-1095). `shear_at_T` averages the 4
      !! surrounding corner sh_Ds values.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sh_dd(nx, ny), sh_dt(nx, ny)
      real(wp), intent(in) :: sh_ds(nx + 1, ny + 1)
      real(wp), intent(in) :: i_ec2
      real(wp), intent(in) :: pres_mice(nx, ny), mice(nx, ny)
      real(wp), intent(in) :: del_sh_min_pr(nx, ny)
      real(wp), intent(out) :: del_sh(nx, ny), zeta(nx, ny)
      integer :: i, j
      real(wp) :: shear_at_t, denom

      do concurrent(j=1:ny, i=1:nx) local(shear_at_t, denom)
         shear_at_t = 0.25_wp*((sh_ds(i, j) + sh_ds(i + 1, j + 1)) + &
                               (sh_ds(i, j + 1) + sh_ds(i + 1, j)))
         del_sh(i, j) = sqrt(sh_dd(i, j)**2 + i_ec2*(sh_dt(i, j)**2 + shear_at_t**2))
         denom = max(del_sh(i, j), del_sh_min_pr(i, j)*pres_mice(i, j))
         if (denom /= 0.0_wp) then
            zeta(i, j) = 0.5_wp*pres_mice(i, j)*mice(i, j)/denom
         else
            zeta(i, j) = 0.0_wp
         end if
      end do
   end subroutine evp_zeta_impl

   pure subroutine evp_stress_relax_impl(zeta, sh_dd, sh_dt, pres_mice, mice, &
                                         i_1pdt_t, dt_2tdamp, i_ec2, str_d, str_t, nx, ny)
      !! str_d/str_t semi-implicit relax (:1124-1134), non-weak_low_shear
      !! branch only.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: zeta(nx, ny), sh_dd(nx, ny), sh_dt(nx, ny)
      real(wp), intent(in) :: pres_mice(nx, ny), mice(nx, ny)
      real(wp), intent(in) :: i_1pdt_t, dt_2tdamp, i_ec2
      real(wp), intent(inout) :: str_d(nx, ny), str_t(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         str_d(i, j) = i_1pdt_t*(str_d(i, j) + dt_2tdamp* &
                                 (zeta(i, j)*sh_dd(i, j) - 0.5_wp*pres_mice(i, j)*mice(i, j)))
         str_t(i, j) = i_1pdt_t*(str_t(i, j) + (i_ec2*dt_2tdamp)*(zeta(i, j)*sh_dt(i, j)))
      end do
   end subroutine evp_stress_relax_impl

   pure subroutine evp_str_s_relax_impl(areaT, zeta, sh_ds, mi_ratio_a_q, &
                                        i_1pdt_t, dt_2tdamp, i_ec2, str_s, nx, ny)
      !! str_s relax (:1137-1143). Corners in [1,nx+1]x[1,ny+1]; the 4
      !! surrounding T-cells at an array-edge corner are handled by
      !! `zeta`'s own ghost values (zero-mass ghost cells => zeta=0
      !! there, contributing nothing) — no special-case branch needed.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: zeta(nx, ny)
      real(wp), intent(in) :: sh_ds(nx + 1, ny + 1)
      real(wp), intent(in) :: mi_ratio_a_q(nx + 1, ny + 1)
      real(wp), intent(in) :: i_1pdt_t, dt_2tdamp, i_ec2
      real(wp), intent(inout) :: str_s(nx + 1, ny + 1)
      integer :: ic, jc
      real(wp) :: zeta_sw, zeta_se, zeta_nw, zeta_ne, a_sw, a_se, a_nw, a_ne
      real(wp) :: weighted_zeta

      do concurrent(jc=1:ny + 1, ic=1:nx + 1) &
         local(zeta_sw, zeta_se, zeta_nw, zeta_ne, a_sw, a_se, a_nw, a_ne, weighted_zeta)
         if (ic == 1 .or. ic == nx + 1 .or. jc == 1 .or. jc == ny + 1) then
            ! Array-edge corner: zeta/areaT have no defined 4th neighbour;
            ! these are ghost/land corners under the ghost policy and are
            ! never read by a physical-interior momentum face. Leave
            ! str_s untouched.
            continue
         else
            zeta_sw = zeta(ic - 1, jc - 1)
            zeta_se = zeta(ic, jc - 1)
            zeta_nw = zeta(ic - 1, jc)
            zeta_ne = zeta(ic, jc)
            a_sw = areaT(ic - 1, jc - 1)
            a_se = areaT(ic, jc - 1)
            a_nw = areaT(ic - 1, jc)
            a_ne = areaT(ic, jc)
            weighted_zeta = ((a_sw*zeta_sw + a_ne*zeta_ne) + (a_se*zeta_se + a_nw*zeta_nw))
            str_s(ic, jc) = i_1pdt_t*(str_s(ic, jc) + (i_ec2*dt_2tdamp)* &
                                      (weighted_zeta*mi_ratio_a_q(ic, jc)*sh_ds(ic, jc)))
         end if
      end do
   end subroutine evp_str_s_relax_impl

   pure subroutine evp_copy_u_impl(ui, u_tmp, nx, ny)
      !! `u_tmp = ui` (full array — the v-momentum MUST read pre-update
      !! u, D4). Explicit `do concurrent` element copy rather than a bare
      !! whole-array assignment (repo convention for device-resident
      !! arrays — see `rdb_ml_dynamics`'s `h_layer0` save pattern).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: ui(nx + 1, ny)
      real(wp), intent(out) :: u_tmp(nx + 1, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         u_tmp(i, j) = ui(i, j)
      end do
   end subroutine evp_copy_u_impl

   pure subroutine evp_zero_stress_impl(fxoc, fyoc, nx, ny)
      !! Zero the subcycle-averaged ice->ocean stress accumulators via an
      !! explicit `do concurrent` device kernel (F1): `fxoc`/`fyoc` are
      !! copyin-mapped device-resident arrays, so a host `= 0.0_wp` would
      !! zero only the HOST copy and leave the device copy carrying the
      !! prior call's average (the accumulate below would then converge to
      !! S/(N-1) instead of S/N). Runs on the device-present arrays; inert
      !! no-op on host builds.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: fxoc(nx + 1, ny)
      real(wp), intent(out) :: fyoc(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         fxoc(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         fyoc(i, j) = 0.0_wp
      end do
   end subroutine evp_zero_stress_impl

   pure subroutine evp_u_momentum_impl(idxCu, idyCu, dy2h, dx2q, iareaCu, mask_u, mi_u, mi_v, &
                                       q, str_d, str_t, str_s, uo, vo, tau_ax, ui, vi, &
                                       fxoc, m_neglect, i_cdrhodt, cdrho, dt, &
                                       nx_phys, ny_phys, nghost, nx, ny, a_u, a_face_on)
      !! u-momentum (:1172-1231, requirement 1: fxic_now carries the FULL
      !! str_t force term). Loop over u-faces `ng+1..ng+nxp+1` x
      !! `ng+1..ng+nyp` — each iteration writes only its own face.
      !!
      !! PR 62 (`a_face_on`): weights BOTH the wind (`tau_ax`) AND the
      !! ice-ocean drag (`drag_u`) by the face ice concentration `a_u`, in
      !! the momentum balance ONLY — `fxoc` stays unweighted (per unit ice
      !! area) so `ice_ocean_stress_flux_impl`'s `a_u*fxoc` on the coupler
      !! side is the ocean's share (see the module docstring D7 + the
      !! `ice_ocean_stress_flux` F5 caveat). Weighting the wind alone would
      !! convert today's leak (zero at steady free drift) into a permanent
      !! one — do not "simplify" this to a single weighted term.
      !!
      !! The `a_fac > 0.0` branch is a MANDATORY 0/0 guard, not defensive
      !! tidiness: at an ice-free face `mi_u = 0`, so a naive
      !! `a_fac*drag_u` collapses the denominator to `m_neglect` alone
      !! against a generally-nonzero `dt*fxic_now`, producing `O(1e30)` on
      !! the first substep. The `else` branch (`uio_c = 0` => `ui = uo`) is
      !! the SIS2 limit (`set_wind_stresses_C`'s `else WindStr_x_Cu = 0.0`)
      !! reached without the division hazard. Do NOT floor `a_fac` instead
      !! of branching — a `max(a_fac, eps)` floor reintroduces a
      !! (much smaller but nonzero) ghost-drift artefact.
      !!
      !! The drag PREDICTOR (`b_vel0`/`uio_pred`, below) is deliberately NOT
      !! folded by `a_u`: conservation depends only on `drag_u`'s use in the
      !! `uio_c`/`fxoc` pair (§3.3 of the PR-62 plan), not on the predictor's
      !! accuracy, and `drag_u`'s own `max(uio_init**2, ...)` converges to
      !! the exact quadratic drag as the substep loop converges regardless.
      integer, intent(in) :: nx_phys, ny_phys, nghost, nx, ny
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in) :: dy2h(nx, ny), dx2q(nx + 1, ny + 1)
      real(wp), intent(in) :: iareaCu(nx + 1, ny)
      real(wp), intent(in) :: mask_u(nx + 1, ny)
      real(wp), intent(in) :: mi_u(nx + 1, ny), mi_v(nx, ny + 1)
      real(wp), intent(in) :: q(nx + 1, ny + 1)
      real(wp), intent(in) :: str_d(nx, ny), str_t(nx, ny), str_s(nx + 1, ny + 1)
      real(wp), intent(in) :: uo(nx + 1, ny), vo(nx, ny + 1)
      real(wp), intent(in) :: tau_ax(nx + 1, ny)
      real(wp), intent(inout) :: ui(nx + 1, ny)
      real(wp), intent(in) :: vi(nx, ny + 1)
      real(wp), intent(inout) :: fxoc(nx + 1, ny)
      real(wp), intent(in) :: m_neglect, i_cdrhodt, cdrho, dt
      real(wp), intent(in) :: a_u(nx + 1, ny)
         !! PR 62: face ice concentration. Valid ONLY when `a_face_on`.
      logical, intent(in) :: a_face_on
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: cor, f2dt_u, i1_f2dt2_u
      real(wp) :: azon, bzon, czon, dzon
      real(wp) :: fxic_now, v2_at_u, uio_init
      real(wp) :: m_uio_explicit, b_vel0, uio_pred, drag_u, uio_c
      real(wp) :: a_fac, tau_eff, drag_eff

      i_lo = nghost + 1
      i_hi = nghost + nx_phys + 1
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(cor, f2dt_u, i1_f2dt2_u, azon, bzon, czon, dzon, fxic_now, v2_at_u, &
               uio_init, m_uio_explicit, b_vel0, uio_pred, drag_u, uio_c, &
               a_fac, tau_eff, drag_eff)
         azon = 0.25_wp*mi_v(i, j + 1)*q(i, j + 1)
         bzon = 0.25_wp*mi_v(i - 1, j + 1)*q(i, j + 1)
         czon = 0.25_wp*mi_v(i - 1, j)*q(i, j)
         dzon = 0.25_wp*mi_v(i, j)*q(i, j)

         cor = 0.25_wp*(q(i, j + 1)*(mi_v(i, j + 1)*vi(i, j + 1) + mi_v(i - 1, j + 1)*vi(i - 1, j + 1)) + &
                        q(i, j)*(mi_v(i - 1, j)*vi(i - 1, j) + mi_v(i, j)*vi(i, j)))

         f2dt_u = dt*4.0_wp*((azon**2 + czon**2) + (bzon**2 + dzon**2))
         i1_f2dt2_u = 1.0_wp/(1.0_wp + dt*f2dt_u)

         fxic_now = idxCu(i, j)*(str_d(i, j) - str_d(i - 1, j)) + &
                    (idyCu(i, j)*(dy2h(i, j)*str_t(i, j) - dy2h(i - 1, j)*str_t(i - 1, j)) + &
                     idxCu(i, j)*(dx2q(i, j + 1)*str_s(i, j + 1) - dx2q(i, j)*str_s(i, j)))* &
                    iareaCu(i, j)

         v2_at_u = 0.25_wp*(((vi(i - 1, j + 1) - vo(i - 1, j + 1))**2 + &
                             (vi(i, j) - vo(i, j))**2) + &
                            ((vi(i, j + 1) - vo(i, j + 1))**2 + &
                             (vi(i - 1, j) - vo(i - 1, j))**2))

         uio_init = ui(i, j) - uo(i, j)

         ! TWO FULLY SEPARATE ARMS, not `a_fac = 1.0_wp` feeding one shared
         ! expression: the off arm below is TEXTUALLY UNCHANGED from the
         ! pre-PR kernel, predictor included.  A shared `tau_eff`/`drag_eff`
         ! computed via `a_fac = 1.0_wp` is mathematically exact (1.0*x==x
         ! in IEEE) but NVHPC's GPU codegen does not guarantee identical FMA
         ! contraction/rounding across ~5000 chained substeps for two
         ! syntactically different expression trees that merely evaluate to
         ! the same VALUE — confirmed empirically (`a_face_full_cover_
         ! bitident` drifted ~1e-14 rel under exactly that construction
         ! before this fix; the `7d283fe8` GPU FMA precedent, CLAUDE.md
         ! Gotchas).  Keep the off arm untouched, full stop.
         if (a_face_on) then
            a_fac = a_u(i, j)
            tau_eff = a_fac*tau_ax(i, j)

            drag_u = 0.0_wp
            if (mask_u(i, j) > 0.0_wp) then
               m_uio_explicit = uio_init*mi_u(i, j) + dt*(cor*mi_u(i, j) + (fxic_now + tau_eff))
               b_vel0 = mi_u(i, j)*i_cdrhodt + (sqrt(uio_init**2 + v2_at_u) - abs(uio_init))
               if (b_vel0**2 > EVP_DRAG_LINEARIZE_THRESHOLD*i_cdrhodt*abs(m_uio_explicit)) then
                  if (b_vel0 /= 0.0_wp) then
                     uio_pred = m_uio_explicit*i_cdrhodt/b_vel0
                  else
                     uio_pred = 0.0_wp
                  end if
               else
                  uio_pred = 0.5_wp*(sqrt(b_vel0**2 + 4.0_wp*i_cdrhodt*abs(m_uio_explicit)) - b_vel0)
               end if
               drag_u = cdrho*sqrt(max(uio_init**2, uio_pred**2) + v2_at_u)
            end if

            drag_eff = a_fac*drag_u
            if (a_fac > 0.0_wp) then
               uio_c = mask_u(i, j)*(mi_u(i, j)*((ui(i, j) + dt*cor)*i1_f2dt2_u - uo(i, j)) + &
                                     dt*(fxic_now + tau_eff))/ &
                       (mi_u(i, j) + m_neglect + dt*drag_eff)
            else
               uio_c = 0.0_wp     ! no ice at either neighbour: no ice momentum here
            end if
         else
            drag_u = 0.0_wp
            if (mask_u(i, j) > 0.0_wp) then
               m_uio_explicit = uio_init*mi_u(i, j) + dt*(cor*mi_u(i, j) + (fxic_now + tau_ax(i, j)))
               b_vel0 = mi_u(i, j)*i_cdrhodt + (sqrt(uio_init**2 + v2_at_u) - abs(uio_init))
               if (b_vel0**2 > EVP_DRAG_LINEARIZE_THRESHOLD*i_cdrhodt*abs(m_uio_explicit)) then
                  if (b_vel0 /= 0.0_wp) then
                     uio_pred = m_uio_explicit*i_cdrhodt/b_vel0
                  else
                     uio_pred = 0.0_wp
                  end if
               else
                  uio_pred = 0.5_wp*(sqrt(b_vel0**2 + 4.0_wp*i_cdrhodt*abs(m_uio_explicit)) - b_vel0)
               end if
               drag_u = cdrho*sqrt(max(uio_init**2, uio_pred**2) + v2_at_u)
            end if

            uio_c = mask_u(i, j)*(mi_u(i, j)*((ui(i, j) + dt*cor)*i1_f2dt2_u - uo(i, j)) + &
                                  dt*(fxic_now + tau_ax(i, j)))/ &
                    (mi_u(i, j) + m_neglect + dt*drag_u)
         end if

         ui(i, j) = (uio_c + uo(i, j))*mask_u(i, j)
         fxoc(i, j) = fxoc(i, j) + drag_u*uio_c          ! UNWEIGHTED: the coupler applies a_u
      end do
   end subroutine evp_u_momentum_impl

   pure subroutine evp_v_momentum_impl(idyCv, idxCv, dx2h, dy2q, iareaCv, mask_v, mi_v, mi_u, &
                                       q, str_d, str_t, str_s, uo, vo, tau_ay, u_tmp, vi, &
                                       fyoc, m_neglect, i_cdrhodt, cdrho, dt, &
                                       nx_phys, ny_phys, nghost, nx, ny, a_v, a_face_on)
      !! v-momentum (:1257-1334, mirror of u). D4: reads `u_tmp` (the
      !! PRE-update u), never the just-updated `ui`. **Minus** on the
      !! str_t divergence term (:1263-1267).
      !!
      !! PR 62 (`a_face_on`): exact mirror of `evp_u_momentum_impl`'s
      !! weighting — see that kernel's docstring for the full rationale
      !! (both terms weighted, `fyoc` unweighted, the `a_fac > 0.0` 0/0
      !! guard, and the untouched drag predictor).
      integer, intent(in) :: nx_phys, ny_phys, nghost, nx, ny
      real(wp), intent(in) :: idyCv(nx, ny + 1), idxCv(nx, ny + 1)
      real(wp), intent(in) :: dx2h(nx, ny), dy2q(nx + 1, ny + 1)
      real(wp), intent(in) :: iareaCv(nx, ny + 1)
      real(wp), intent(in) :: mask_v(nx, ny + 1)
      real(wp), intent(in) :: mi_v(nx, ny + 1), mi_u(nx + 1, ny)
      real(wp), intent(in) :: q(nx + 1, ny + 1)
      real(wp), intent(in) :: str_d(nx, ny), str_t(nx, ny), str_s(nx + 1, ny + 1)
      real(wp), intent(in) :: uo(nx + 1, ny), vo(nx, ny + 1)
      real(wp), intent(in) :: tau_ay(nx, ny + 1)
      real(wp), intent(in) :: u_tmp(nx + 1, ny)
      real(wp), intent(inout) :: vi(nx, ny + 1)
      real(wp), intent(inout) :: fyoc(nx, ny + 1)
      real(wp), intent(in) :: m_neglect, i_cdrhodt, cdrho, dt
      real(wp), intent(in) :: a_v(nx, ny + 1)
         !! PR 62: face ice concentration. Valid ONLY when `a_face_on`.
      logical, intent(in) :: a_face_on
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: cor, f2dt_v, i1_f2dt2_v
      real(wp) :: amer, bmer, cmer, dmer
      real(wp) :: fyic_now, u2_at_v, vio_init
      real(wp) :: m_vio_explicit, b_vel0, vio_pred, drag_v, vio_c
      real(wp) :: a_fac, tau_eff, drag_eff

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys + 1

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(cor, f2dt_v, i1_f2dt2_v, amer, bmer, cmer, dmer, fyic_now, u2_at_v, &
               vio_init, m_vio_explicit, b_vel0, vio_pred, drag_v, vio_c, &
               a_fac, tau_eff, drag_eff)
         amer = 0.25_wp*mi_u(i, j - 1)*q(i, j)
         bmer = 0.25_wp*mi_u(i + 1, j - 1)*q(i + 1, j)
         cmer = 0.25_wp*mi_u(i + 1, j)*q(i + 1, j)
         dmer = 0.25_wp*mi_u(i, j)*q(i, j)

         cor = -0.25_wp*(q(i, j)*(mi_u(i, j - 1)*u_tmp(i, j - 1) + mi_u(i, j)*u_tmp(i, j)) + &
                         q(i + 1, j)*(mi_u(i + 1, j - 1)*u_tmp(i + 1, j - 1) + &
                                      mi_u(i + 1, j)*u_tmp(i + 1, j)))

         f2dt_v = dt*4.0_wp*((amer**2 + cmer**2) + (bmer**2 + dmer**2))
         i1_f2dt2_v = 1.0_wp/(1.0_wp + dt*f2dt_v)

         fyic_now = idyCv(i, j)*(str_d(i, j) - str_d(i, j - 1)) + &
                    (-idxCv(i, j)*(dx2h(i, j)*str_t(i, j) - dx2h(i, j - 1)*str_t(i, j - 1)) + &
                     idyCv(i, j)*(dy2q(i + 1, j)*str_s(i + 1, j) - dy2q(i, j)*str_s(i, j)))* &
                    iareaCv(i, j)

         u2_at_v = 0.25_wp*(((u_tmp(i + 1, j - 1) - uo(i + 1, j - 1))**2 + &
                             (u_tmp(i, j) - uo(i, j))**2) + &
                            ((u_tmp(i + 1, j) - uo(i + 1, j))**2 + &
                             (u_tmp(i, j - 1) - uo(i, j - 1))**2))

         vio_init = vi(i, j) - vo(i, j)

         ! TWO FULLY SEPARATE ARMS -- see evp_u_momentum_impl's comment at
         ! the mirror site for why the off arm (predictor included) must
         ! stay textually untouched rather than routed through a shared
         ! `a_fac = 1.0_wp` expression.
         if (a_face_on) then
            a_fac = a_v(i, j)
            tau_eff = a_fac*tau_ay(i, j)

            drag_v = 0.0_wp
            if (mask_v(i, j) > 0.0_wp) then
               m_vio_explicit = vio_init*mi_v(i, j) + dt*(cor*mi_v(i, j) + (fyic_now + tau_eff))
               b_vel0 = mi_v(i, j)*i_cdrhodt + (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
               if (b_vel0**2 > EVP_DRAG_LINEARIZE_THRESHOLD*i_cdrhodt*abs(m_vio_explicit)) then
                  if (b_vel0 /= 0.0_wp) then
                     vio_pred = m_vio_explicit*i_cdrhodt/b_vel0
                  else
                     vio_pred = 0.0_wp
                  end if
               else
                  vio_pred = 0.5_wp*(sqrt(b_vel0**2 + 4.0_wp*i_cdrhodt*abs(m_vio_explicit)) - b_vel0)
               end if
               drag_v = cdrho*sqrt(max(vio_init**2, vio_pred**2) + u2_at_v)
            end if

            drag_eff = a_fac*drag_v
            if (a_fac > 0.0_wp) then
               vio_c = mask_v(i, j)*(mi_v(i, j)*((vi(i, j) + dt*cor)*i1_f2dt2_v - vo(i, j)) + &
                                     dt*(fyic_now + tau_eff))/ &
                       (mi_v(i, j) + m_neglect + dt*drag_eff)
            else
               vio_c = 0.0_wp     ! no ice at either neighbour: no ice momentum here
            end if
         else
            drag_v = 0.0_wp
            if (mask_v(i, j) > 0.0_wp) then
               m_vio_explicit = vio_init*mi_v(i, j) + dt*(cor*mi_v(i, j) + (fyic_now + tau_ay(i, j)))
               b_vel0 = mi_v(i, j)*i_cdrhodt + (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
               if (b_vel0**2 > EVP_DRAG_LINEARIZE_THRESHOLD*i_cdrhodt*abs(m_vio_explicit)) then
                  if (b_vel0 /= 0.0_wp) then
                     vio_pred = m_vio_explicit*i_cdrhodt/b_vel0
                  else
                     vio_pred = 0.0_wp
                  end if
               else
                  vio_pred = 0.5_wp*(sqrt(b_vel0**2 + 4.0_wp*i_cdrhodt*abs(m_vio_explicit)) - b_vel0)
               end if
               drag_v = cdrho*sqrt(max(vio_init**2, vio_pred**2) + u2_at_v)
            end if

            vio_c = mask_v(i, j)*(mi_v(i, j)*((vi(i, j) + dt*cor)*i1_f2dt2_v - vo(i, j)) + &
                                  dt*(fyic_now + tau_ay(i, j)))/ &
                    (mi_v(i, j) + m_neglect + dt*drag_v)
         end if

         vi(i, j) = (vio_c + vo(i, j))*mask_v(i, j)
         fyoc(i, j) = fyoc(i, j) + drag_v*vio_c          ! UNWEIGHTED: the coupler applies a_v
      end do
   end subroutine evp_v_momentum_impl

   pure subroutine evp_average_stress_impl(mask_u, mask_v, fxoc, fyoc, evp_sub_steps, nx, ny)
      !! `fxoc *= mask_u/evp_sub_steps`, `fyoc *= mask_v/evp_sub_steps`
      !! (:1415-1441).
      integer, intent(in) :: evp_sub_steps, nx, ny
      real(wp), intent(in) :: mask_u(nx + 1, ny)
      real(wp), intent(in) :: mask_v(nx, ny + 1)
      real(wp), intent(inout) :: fxoc(nx + 1, ny)
      real(wp), intent(inout) :: fyoc(nx, ny + 1)
      integer :: i, j
      real(wp) :: i_sub_steps

      i_sub_steps = 1.0_wp/real(evp_sub_steps, wp)
      do concurrent(j=1:ny, i=1:nx + 1)
         fxoc(i, j) = fxoc(i, j)*(mask_u(i, j)*i_sub_steps)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         fyoc(i, j) = fyoc(i, j)*(mask_v(i, j)*i_sub_steps)
      end do
   end subroutine evp_average_stress_impl

end module rdb_ice_evp
