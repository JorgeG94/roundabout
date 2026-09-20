!! Split-explicit RK2 dynamics state + driver step.
module rdb_ocean_dyn
   !! Driver-level state + step routines for the outer (baroclinic) +
   !! inner (barotropic) split-explicit RK2 integrator of the ocean
   !! dynamical core.  Carries the substep ratio `n_inner`, the fast-mode
   !! time-averaging accumulators, and energy/CFL diagnostic slots.
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t, &
                                       barotropic_workstate_enter_data_impl, &
                                       barotropic_workstate_exit_data_impl
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear_interior
   use rdb_ocean_bt_wide, only: bt_wide_t, bt_wide_substep
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, OBC_WALL
   use rdb_ocean_periodic, only: ocean_periodic_wrap_state, &
                                 ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_centre_3d, &
                                 ocean_periodic_wrap_face_x_3d, &
                                 ocean_periodic_wrap_face_y_3d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_state
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_ocean_halo, only: ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y, &
                             ocean_halo_bt_group_2d, ocean_halo_centre, &
                             ocean_halo_face_x, ocean_halo_face_y
   use rdb_ocean_ghost_poison, only: ocean_poison_ghost_bands
   use rdb_ocean_sponge, only: ocean_sponge_apply, ocean_sponge_apply_tracers, &
                               ocean_sponge_apply_maps, ocean_sponge_t
   use rdb_ocean_obc_baroclinic, only: ocean_obc_apply_baroclinic, &
                                       ocean_obc_fill_ghosts, &
                                       ocean_obc_refill_ghost_ssh, &
                                       ocean_obc_update_reservoirs
   use rdb_barotropic_coupling, only: derive_bt_from_layers, &
                                      compute_h_face_upstream, &
                                      sum_slow_tendencies_into_F_slow, &
                                      subtract_fast_cor_ref, &
                                      set_cor_ref_velocity, &
                                      face_depth_mean_u, face_depth_mean_v, &
                                      face_depth_mean_rem_u, face_depth_mean_rem_v, &
                                      apply_bt_correction, &
                                      snapshot_eta_PF, compute_pbce, &
                                      compute_gtot_faces, compute_e_anom, &
                                      compute_bt_rem, reset_bt_rem, &
                                      compute_bt_rem_wave_drag, mask_bt_rem, &
                                      set_local_BT_cont_types
   use rdb_ocean_bt_budget_probe, only: print_bt_budget
   use rdb_ocean_ke_probe, only: ke_probe_t, ke_probe_sample, ke_probe_coradv_split
   use rdb_ocean_chksum, only: chksum_probe_t, chksum_state, chksum_bt, chksum_hotface
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_porous, only: porous_update_face_areas
   use rdb_ocean_metrics, only: ocean_metrics_t, &
                                metrics_fill_cartesian, metrics_fill_spherical, &
                                metrics_finalize, metrics_fill_coriolis, &
                                GRID_CONFIG_CARTESIAN, GRID_CONFIG_SPHERICAL, &
                                CORIOLIS_SCHEME_BETA_PLANE, &
                                parse_grid_config, parse_coriolis_scheme
   use rdb_continuity, only: continuity_t, &
                             continuity_compute_fluxes_barotropic, &
                             continuity_apply_fluxes_barotropic, &
                             continuity_tracer_step_split, &
                             continuity_tracer_drain, &
                             TR_MODE_ADVECT, TR_MODE_ACCUMULATE, TR_MODE_NONE
   use rdb_coriolis_adv, only: coriolis_adv_t, &
                               coriolis_adv_compute_tendencies_barotropic, &
                               coriolis_adv_apply_tendencies_barotropic, &
                               coriolis_adv_compute_tendencies, &
                               coriolis_adv_apply_tendencies
   use rdb_eos, only: eos_t
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       ocean_pressure_force_apply
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies, &
                                             ocean_horizontal_viscosity_apply_tendencies, &
                                             ocean_horizontal_viscosity_compute_ke_diss
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, ocean_lateral_mix_compute
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_mle, only: ocean_mle_t, mle_compute_transports
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports
   use rdb_ocean_redi, only: ocean_redi_t, redi_calc_coeffs, redi_apply_flux
   use rdb_ocean_varmix, only: ocean_varmix_t, varmix_compute
   use rdb_ocean_meke, only: ocean_meke_t, meke_step, meke_backscatter_apply
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t, wavespeed_compute
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, &
                                    ocean_bottom_drag_compute_tendencies, &
                                    ocean_bottom_drag_apply_tendencies, &
                                    ocean_channel_drag_compute_tendencies, &
                                    ocean_channel_drag_apply_tendencies
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_compute_tendencies, &
                                       ocean_surface_stress_apply_tendencies
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_apply_sw_penetration, &
                                     ocean_surface_restore_apply_tracers
   use rdb_ocean_geothermal, only: ocean_geothermal_t, &
                                   ocean_geothermal_apply_tracers
   use rdb_ocean_ideal_age, only: ocean_ideal_age_apply, ocean_ideal_age_reset_surface, &
                                  ocean_ideal_age_young_val
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t, &
                                           compute_w_from_continuity, &
                                           tracer_advect_vertical
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t, &
                                     tracer_hdiff
   use rdb_ocean_vdiff, only: ocean_vdiff_t, &
                              vdiff_apply_momentum, &
                              vdiff_apply_tracers
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, &
                             vmix_apply_kpp_overlay, vmix_add_kv_ml_invz2, &
                             vmix_apply_nonlocal_tendencies, vmix_assemble, &
                             vmix_apply_convection, vmix_split_kd_heat_salt, &
                             VMIX_INTERIOR_PP81
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_compute, epbl_merge_into_kv_kt
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t, kappa_shear_compute, &
                                    kappa_shear_merge_into_kv_kt
   use rdb_ocean_tidal_mixing, only: ocean_tidal_mixing_t, tidal_mixing_compute, &
                                     tidal_mixing_merge_into_kt
   use rdb_ocean_tides, only: ocean_tides_t, tides_update_eta_eq, &
                              tides_update_eta_sal
   use rdb_ocean_p_surf, only: ocean_p_surf_t, p_surf_update_seam
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN
   use rdb_ocean_remap, only: ocean_apply_ale_remap_step
   use rdb_ocean_min_thickness, only: ocean_apply_conservative_min_thickness
   use rdb_profiler, only: profiler_start, profiler_stop
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_dyn_t
   public :: ocean_dyn_step_barotropic
   public :: ocean_dyn_step
   public :: ocean_dyn_step_split
   public :: ocean_porous_refresh
   public :: apply_velocity_truncation
   public :: ocean_dt_tracer_advect_ratios_ok
   public :: ocean_dyn_flush_tracer_window
   public :: isopycnal_vanish_tol
   public :: ocean_dyn_enable_bt_wide
   public :: accel_visc_rem_snapshot
   public :: accel_visc_rem_reweight

#ifdef RDB_ENABLE_TESTING
   public :: reset_vanished_layer_velocities
   public :: mask_layer_velocities
#endif

   ! ---- Baroclinic-stability diagnostic (debug-gated, default OFF) ----
   ! Probes in `run_stage_split` print max|hTr_S/h_layer - bcdiag_S_ref|
   ! after each kernel call.  All probes early-return when disabled, so
   ! production paths are unaffected; tests set `bcdiag_enabled = .true.`.
   logical, public, save :: bcdiag_enabled = .false.
   real(wp), public, save :: bcdiag_S_ref = 35.0_wp
   integer, public, save :: bcdiag_step_limit = 2
      !! Probe only the first `bcdiag_step_limit` outer steps (keeps trace bounded).

   ! Outer split-explicit time-scheme selector (`&ocean_bt_nml
   ! split_scheme`, SPEC §4 S3/S4).
   integer, parameter, public :: SPLIT_SCHEME_SSP_RK2 = 0
      !! Two identical stages + SSP average (`&ocean_bt_nml split_scheme
      !! = "ssp_rk2"`).  **EXPERIMENTAL** — a label on the ANSWER, not a
      !! deprecation: the path is fully supported, fully tested, and the
      !! stability suite runs an `ssp_rk2` twin of every case whose
      !! namelist does not pin a scheme.  It has the widest envelope
      !! (every vcoord including `eulerian_z`, wet/dry,
      !! `dt_tracer_advect_ratio > 1`) and is the only scheme wired
      !! through the windowed tracer-advection path, which is why six
      !! shipped namelists pin it.
      !!
      !! **What it costs you, measured.**  The two-stage average amplifies
      !! an internal gravity wave by `√(1 + (ω·dt)⁴/4)` per step, so it
      !! MANUFACTURES energy from a motionless stratified state.  On
      !! `validation_examples/ocean/eady/resting_stratified_channel.nml`
      !! — a flat-bottomed periodic channel at rest, stably stratified,
      !! seeded with ±0.5 mK of noise and then left alone, with NO energy
      !! source of any kind — En reaches **2.992E-05 m²/s² by day 25**
      !! (7.7 mm/s of current out of nothing) and is still climbing on a
      !! 2.5-day e-folding.  `pred_corr` on the identical file holds
      !! **1.739E-09** (17 000× less, 83-day e-folding).  The cause is
      !! the OUTER time splitting and nothing else: substituting the
      !! Coriolis form (`sadourny_hk`, `sadourny_energy`), the ALE remap
      !! (`vcoord = sigma`) and the PGF form (`fv_lite`) each moves the
      !! answer by < 0.1 % — all three EXONERATED; `eady_dT_dz = 0` drops
      !! En 119×, and REMOVING the lateral viscosity RAISES it —
      !! viscosity damps the mode, it is not its source.
      !!
      !! Practically: the error is a gravity-wave-scale numerical noise
      !! floor that grows with `(ω·dt)⁴`, so it is set by how hard you
      !! push `dt` against the internal-wave period.  A forced,
      !! energetic, viscous configuration runs decades above that floor
      !! and never notices it; a quiescent, weakly-damped or long
      !! spin-up one does, and there the manufactured energy IS the
      !! signal.  Carried as a scoped XFAIL on
      !! `resting_stratified_channel__ssp_rk2`.
   integer, parameter, public :: SPLIT_SCHEME_PRED_CORR = 1
      !! MOM6 predictor-corrector: off-centred predictor at `pc_be·dt`,
      !! slow tendencies on the `u_av`/`h_av` step time-means, ONE
      !! prognostic update in the corrector, forward-backward
      !! gravity-wave pairing.  **The DEFAULT** (`&ocean_bt_nml
      !! split_scheme = "pred_corr"`, since 2026-09-14).
      !!
      !! Neutrally stable to `ω·dt = 2`, so it does not have the
      !! resting-state growth described above (17 000× less on the same
      !! file) and it lifts the internal-wave `dt` ceiling.  It costs
      !! about 7 % more wall time per step than `ssp_rk2` — not a
      !! different order.  `validate_config` refuses it fail-loud
      !! outside its v1 envelope (`eulerian_z`, wet/dry,
      !! `dt_tracer_advect_ratio > 1`) rather than degrading silently;
      !! those configurations must pin `ssp_rk2`.
      !!
      !! **Was renamed twice.**  This scheme shipped as `"mom6_pc"`,
      !! then briefly as `"split_rk2"` — a name that collided with
      !! `dynamics/split_rk2/`, the directory holding the outer-loop
      !! machinery of BOTH schemes, so it named the family as well as
      !! one member.  A namelist still carrying EITHER dead spelling
      !! fails loud naming `pred_corr` (`nml_enum`'s `retired` list),
      !! never silently falls back to a default.

   type :: ocean_dyn_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      integer :: split_scheme = SPLIT_SCHEME_PRED_CORR
         !! Outer time-scheme selector (see the module constants).  Kept in
         !! step with the `&ocean_bt_nml split_scheme` default in
         !! `rdb_config.F90` — `ocean_setup` assigns this field from the
         !! config on BOTH branches (`rdb_ocean_setup.F90`, "SPEC S3"), so no
         !! production run reads this value, but a bare `ocean_dyn_t` built by
         !! a unit test inherits it and a divergence would have the Fortran
         !! suite exercising a scheme the model does not ship by default.
         !! The two have diverged before; change them in the same commit.
         !!
         !! **Why the default moved (2026-09-14).**  `ssp_rk2` grows
         !! internal gravity waves out of a stratified REST state — En
         !! 2.992E-05 against pred_corr's 1.739E-09 at day 25 on
         !! `resting_stratified_channel.nml`, 17 000× — and the three
         !! things that kept the predictor-corrector off the default are
         !! all closed: the land-mask NaN (the fast-loop Coriolis
         !! reference was evaluated on the stage-entry `u^n` while the
         !! slow Coriolis used `u_av`; see `set_cor_ref_velocity`), the
         !! `nz = 1` vdiff tridiagonal defect, and SEVEN GPU-only
         !! unit-test failures — `periodic`, `obc_baroclinic`,
         !! `dyn_split`, `ice_restart`, `wetdry` (non-finite within 1-5
         !! steps), plus `restart` bit-exactness and uniform-`p_surf`
         !! gauge invariance — which had ONE cause, and it was not scheme
         !! semantics: `scratch_3d_buffer_t` attached its payload with
         !! `!$acc enter data create`, which carries no host value, and
         !! the pred_corr PREDICTOR deliberately reads
         !! `hv%du_visc`/`dv_visc` without recomputing them (MOM6's
         !! `diffu(u[n-1])` reuse).  At step 1 there is no previous
         !! producer, so the host read the zero `init` promised and the
         !! device read the allocator's leftovers.  `enter_data` now
         !! device-zeroes the payload (`rdb_test_scratch_3d_device` is the
         !! gate), and `ctest -R rdb` is 187/187 on BOTH toolchains.
      real(wp) :: pc_be = 0.6_wp
         !! pred_corr predictor fraction (MOM6 `BE`, 0.6 in the control
         !! run): the predictor advances the provisional velocity to
         !! `dt_pred = pc_be·dt` (SPEC §2 P8); only the corrector takes
         !! the full step.  Unused under ssp_rk2.
      integer :: n_inner = 0
         !! Number of barotropic substeps per outer baroclinic step
         !! (0 = "auto, derived from CFL").
      integer :: outer_step_count = 0
         !! Outer (baroclinic) step counter.
      real(wp) :: dt_inner = 0.0_wp
         !! Inner barotropic substep length (s).
      real(wp) :: dt_outer = 0.0_wp
         !! Outer baroclinic timestep length (s).

      ! ---- Barotropic working state ----
      ! 2D barotropic-substep fields, time-mean accumulators, nonlinear
      ! scratch, and slow-tendency accumulators.  The substep kernels
      ! (`rdb_barotropic_substep`) + coupling kernels
      ! (`rdb_barotropic_coupling`) operate on `dyn%bt_work` directly.
      type(barotropic_workstate_t) :: bt_work

      ! ---- Diagnostic accumulators ----
      real(wp) :: ke_outer = 0.0_wp
         !! Total kinetic energy after the last outer step (diagnostic).
      real(wp) :: cfl_outer = 0.0_wp
         !! Realised outer CFL (diagnostic).
      real(wp) :: cfl_inner_max = 0.0_wp
         !! Max realised inner CFL within the last outer step.

      ! ---- Velocity-truncation clamp (MOM6 MAXVEL) ----
      real(wp) :: maxvel = 0.0_wp
         !! When > 0, face velocities are clipped to `[-maxvel, +maxvel]`
         !! after each outer step (safety net, not a fix).  0 = off.
      real(wp) :: cfl_trunc = 0.0_wp
         !! Advective-CFL truncation threshold (nondim).  When > 0, any
         !! face whose `|u|·dt/dx` exceeds it is clipped to
         !! `0.9·cfl_trunc·dx/dt` (sign preserved), then the `maxvel`
         !! cap is applied.  0 = off (bit-identical).
      integer :: ntrunc_step = 0
         !! Face components CFL-truncated in the most recent outer step
         !! (host-scalar reduction; reset each `apply_velocity_truncation`).
      integer :: ntrunc_total = 0
         !! Cumulative CFL-truncation count over the run (driver logs it).
      integer :: n_nanzero_step = 0
         !! Non-finite (NaN/Inf) face velocities zeroed by the truncation's
         !! NaN-catch in the most recent outer step.  MUST be 0 on a
         !! healthy run; any non-zero is a producer 0/0 upstream (loud —
         !! NaNs used to launder to ±maxvel silently, pdc 8c2fd674).
      integer(int64) :: n_nanzero_total = 0_int64
         !! Cumulative NaN-catch count over the run.
      integer :: bt_nonfin_step = 0
         !! Faces whose BT-correction Δ was non-finite this stage (fold
         !! write skipped; pdc 2d5338f9).  0 on a healthy run.

      ! ---- DT_THERM ratio (MOM6) ----
      integer :: dt_therm_ratio = 1
         !! Tracer/thermo step runs every `dt_therm_ratio` outer steps
         !! with effective dt = ratio · dt.  1 (default) = every step
         !! (bit-identical).  Applies only when `enable_thermodynamics`.
      integer :: dt_tracer_advect_ratio = 1
         !! MOM6 DT_TRACER_ADVECT analogue.  Horizontal tracer advection
         !! fires every `dt_tracer_advect_ratio` outer steps over the
         !! accumulated face transports (`continuity_t%uhtr/vhtr`).
         !! 1 (default) ⇒ every-step path verbatim (bit-identical).
         !! Configure requires `dt_therm_ratio` to be an integer multiple
         !! so the ALE remap never fires mid-accumulation-window.

      ! ---- Debug probe ----
      logical :: debug_bt_budget = .false.
         !! When `.true.`, calls the BT-budget probe once per RK2 stage
         !! after the slow tendencies.  Heavy (D->H transfers + prints).
         !! Default `.false.`; driven by `ocean_debug_bt_budget` namelist.
      type(ke_probe_t) :: ke_probe
         !! Per-segment layer-KE attribution meter (`&ocean_debug_nml
         !! ke_attr`).  Default off ⇒ bit-identical; when on,
         !! serialises the async velocity-apply chain (device waits +
         !! one reduction per segment) — debug only.
      type(chksum_probe_t) :: chksum_probe
         !! MOM6-style per-phase field checksums + HOTFACE argmax rows
         !! (`&ocean_debug_nml chksum` + step window).  Default off ⇒
         !! bit-identical; when on, waits the device + reduces at each
         !! phase seam — the first-diverging-operator attribution probe.
      logical :: accel_visc_rem = .false.
         !! MOM6 `MOM_dynamics_split_RK2` parity: attenuate the slow
         !! EXPLICIT accelerations by the per-layer viscous remnant —
         !! `u_new = u_entry + visc_rem·(u_applied − u_entry)` after the
         !! CorAdv/PGF/hvisc/drag applies, before the BT correction —
         !! so friction-dominated near-massless layers cannot receive a
         !! full-strength dt·F kick (the 2026-07-28 forensics: explicit
         !! force × dt on outcropped mm-layers is the dt=800 blow-up
         !! injector; the fold/vdiff mopped ±30 m/s per stage until
         !! escape).  `&ocean_vdiff_nml accel_visc_rem`; requires
         !! `&ocean_bt_nml correction_visc_rem` (the producer).  Default
         !! off ⇒ bit-identical.  Split path only (v1).
      real(wp), allocatable :: avr_u0(:, :, :)
      real(wp), allocatable :: avr_v0(:, :, :)
         !! accel_visc_rem stage-entry velocity snapshots.  Eagerly
         !! allocated at setup (`configure_ocean_bt`) and device-mapped in
         !! `ocean_dyn_enter_data_impl` — but ONLY when the `accel_visc_rem`
         !! knob is on (unallocated / zero-footprint on the default-off
         !! path); released in destroy/exit_data.

      ! ---- Ghost-band poison (debug) ----
      logical :: poison_ghosts = .false.
         !! When `.true.`, sentinel-NaN the exchange-covered ghost bands at
         !! each outer-step start (before any exchange or kernel).  Any kernel
         !! consuming an unexchanged ghost produces a loud NaN at that step.
         !! Default `.false.` = bit-identical (untaken branch per step).
         !! Driven by `&ocean_mpi_nml poison_ghosts` via `rdb_ocean_setup`.
         !! Requires the `bc` optional argument (edge topology); no-op when
         !! `bc` is absent.

      ! ---- Thermodynamic switch ----
      logical :: enable_thermodynamics = .true.
         !! Mirrors MOM6 `ENABLE_THERMODYNAMICS`.  When `.false.` the step
         !! skips EOS, tracer hdiff, tracer vertical advection/diffusion,
         !! KPP non-local transport, and surface tracer flux.  Tracers
         !! still advect horizontally but with `alpha_T = beta_S = 0` can't
         !! influence dynamics ⇒ adiabatic.  Wire via
         !! `&ocean_setup_nml ocean_enable_thermodynamics`.

      ! ---- Phase-2/3 Lagrangian grounding-stability knobs ----
      real(wp) :: angstrom_h = 0.0_wp
         !! Phase-2/3 copy of `cfg%ocean%isopycnal%angstrom_h` (m).
         !! Used to compute `isopycnal_vanish_tol` for the reset and CFL
         !! gates; 0 ⇒ off ⇒ bit-identical.
      logical :: reset_vanished_u = .false.
         !! Phase 2: when `.true.` AND vcoord is VCOORD_LAGRANGIAN, call
         !! `reset_vanished_layer_velocities` after each
         !! `mask_layer_velocities` (per-stage + RK2-averaged site) to
         !! zero face velocities where BOTH adjacent cell-thicknesses are
         !! at or below `isopycnal_vanish_tol(angstrom_h)`.
      logical :: cfl_ignore_vanished = .false.
         !! Phase 3: when `.true.` AND vcoord is VCOORD_LAGRANGIAN, pass
         !! `vanish_tol` to `compute_max_cfl` and `apply_velocity_truncation`
         !! so vanished-layer face spikes are excluded from MaxCFL /
         !! panic and are zeroed (not CFL-clipped) by truncation.
      logical :: check_h_positive = .false.
         !! DEBUG copy of `cfg%ocean%isopycnal%check_h_positive`.  When
         !! `.true.`, `check_h_positive_or_die` runs after each `h_layer`-
         !! writing stage and aborts on the first negative thickness, naming
         !! the stage.  `.false.` ⇒ never called ⇒ bit-identical.
      ! ---- Wide-halo BT march-in (Phase 3c) ----
      integer :: bt_halo = 0
         !! Wide-halo march-in width (0 = per-substep v1, bit-identical).
         !! Set by `ocean_dyn_enable_bt_wide` after `init` and `enter_data`.
      type(bt_wide_t), allocatable :: bt_wide
         !! Wide shadow state; allocated only when `bt_halo > 0`.

      ! ---- Ideal-age tracer (PR-7) ----
      real(wp) :: ideal_age_young_val = 0.0_wp
         !! Surface-band age value (s), set by `configure_ocean_tracers`
         !! from `&ocean_tracers_nml ideal_age_young_val`.  0 (default) =
         !! today's hard-coded zero reset ⇒ bit-identical.
      real(wp) :: ideal_age_sfc_growth_rate = 0.0_wp
         !! Exponential growth rate of the surface value (1/s), set by
         !! `configure_ocean_tracers` from
         !! `&ocean_tracers_nml ideal_age_sfc_growth_rate`.  0 (default)
         !! ⇒ young_val constant ⇒ no `exp()` ⇒ bit-identical.
   contains
      procedure, non_overridable :: init => ocean_dyn_init
      procedure, non_overridable :: destroy => ocean_dyn_destroy
      procedure, non_overridable :: enter_data => ocean_dyn_enter_data
      procedure, non_overridable :: exit_data => ocean_dyn_exit_data
      procedure, non_overridable :: is_thermo_step => ocean_dyn_is_thermo_step
      procedure, non_overridable :: therm_dt => ocean_dyn_therm_dt
      procedure, non_overridable :: is_tracer_advect_step => ocean_dyn_is_tracer_advect_step
      procedure, non_overridable :: bytes => ocean_dyn_bytes
   end type ocean_dyn_t

contains

   subroutine ocean_porous_refresh(grid, metrics, ms)
      !! Recompute the porous-barrier layer-averaged open-area fractions
      !! from the current layer thicknesses.  No-op (and untouched
      !! placeholder arrays) when `&ocean_porous_nml enable` is off.
      !!
      !! Cadence: ONCE PER OUTER STEP, before the RK2 stages — the
      !! fractions are then held fixed across both stages.  That is
      !! MOM6's cadence (`porous_widths_layer` runs in `step_MOM` ahead
      !! of the dynamics call and the porous type enters the split-RK2
      !! driver `intent(in)`), and it keeps the widths consistent between
      !! the two stages' continuity and Coriolis transports.
      !!
      !! Device-side `do concurrent`: every array it touches is already
      !! mapped by `ocean_metrics_enter_data` / the multilayer state, so
      !! there is no host round-trip and nothing to allocate.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(in) :: ms

      if (.not. metrics%use_porous) return

      call porous_update_face_areas(grid%nx_total, grid%ny_total, ms%nz_ml, &
                                    metrics%porous_eta_interp, &
                                    metrics%porous_mask_depth, &
                                    metrics%por_bed, ms%h_layer, &
                                    metrics%por_dmin_u, metrics%por_dmax_u, &
                                    metrics%por_davg_u, &
                                    metrics%por_dmin_v, metrics%por_dmax_v, &
                                    metrics%por_davg_v, &
                                    metrics%dy_cu, metrics%dx_cv, &
                                    metrics%por_face_area_u, &
                                    metrics%por_face_area_v, &
                                    metrics%dy_cu_bt, metrics%dx_cv_bt)
   end subroutine ocean_porous_refresh

   subroutine ocean_dyn_init(this, grid, nz_ml)
      !! Initialise the barotropic working-state slot.  Pass `nz_ml`
      !! to also allocate the split-driver slow-tendency
      !! accumulators (`F_slow_u/v`, `F_bt_u/v`, `ubt_at_n/vbt_at_n`).
      !! Barotropic-substep unit tests can skip the optional argument since
      !! they don't exercise the split-driver coupling.
      class(ocean_dyn_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      if (present(nz_ml)) then
         call this%bt_work%init(grid, nz_ml=nz_ml)
      else
         call this%bt_work%init(grid)
      end if
      this%is_init = .true.
   end subroutine ocean_dyn_init

   pure function ocean_dyn_is_thermo_step(this) result(yes)
      !! Returns `.true.` if the thermodynamic / tracer kernels should
      !! fire on this outer step.  `dt_therm_ratio <= 1` (default) →
      !! always true (every step is a thermo step, bit-identical to
      !! prior behaviour).  `ratio >= 2` → fires every Nth step,
      !! aligned to `outer_step_count = 0` for the IC snapshot.
      class(ocean_dyn_t), intent(in) :: this
      logical :: yes
      if (this%dt_therm_ratio <= 1) then
         yes = .true.
      else
         yes = (mod(this%outer_step_count, this%dt_therm_ratio) == 0)
      end if
   end function ocean_dyn_is_thermo_step

   pure function ocean_dyn_therm_dt(this, dt) result(dt_th)
      !! Effective `dt` for the thermo / tracer kernels.  When
      !! `dt_therm_ratio = 1` returns `dt` exactly; otherwise returns
      !! `ratio · dt`, since the kernels only fire every Nth step and
      !! must advance by that aggregate interval.
      class(ocean_dyn_t), intent(in) :: this
      real(wp), intent(in) :: dt
      real(wp) :: dt_th
      if (this%dt_therm_ratio <= 1) then
         dt_th = dt
      else
         dt_th = real(this%dt_therm_ratio, wp)*dt
      end if
   end function ocean_dyn_therm_dt

   pure function ocean_dyn_is_tracer_advect_step(this) result(yes)
      !! Returns `.true.` when the windowed horizontal tracer-advect drain
      !! should fire on the current outer step, evaluated with the
      !! PRE-increment `outer_step_count` so it aligns bit-for-bit with
      !! `is_thermo_step()` (which gates the ALE remap on the same count).
      !!
      !! Cadence reasoning: at `outer_step_count = 0` the IC-aligned first
      !! window has not accumulated yet, so the drain must NOT fire then.
      !! It fires when `mod(outer_step_count, ratio) == 0 .and.
      !! outer_step_count > 0` — i.e. exactly `ratio` accumulate-steps have
      !! completed since the last reset.  Because the config constraint
      !! forces `dt_therm_ratio` to be an integer multiple of
      !! `dt_tracer_advect_ratio`, every thermo step (where ALE fires) is
      !! also a tracer-advect step, so the drain always precedes the ALE
      !! remap and no accumulation window is ever cut by a remap.
      !! `ratio <= 1` ⇒ never (the every-step bypass owns that case).
      class(ocean_dyn_t), intent(in) :: this
      logical :: yes
      if (this%dt_tracer_advect_ratio <= 1) then
         yes = .false.
      else
         yes = (this%outer_step_count > 0) .and. &
               (mod(this%outer_step_count, this%dt_tracer_advect_ratio) == 0)
      end if
   end function ocean_dyn_is_tracer_advect_step

   pure function ocean_dt_tracer_advect_ratios_ok(dt_therm_ratio, &
                                                  dt_tracer_advect_ratio) result(ok)
      !! Configure-time validity of the (`dt_therm_ratio`,
      !! `dt_tracer_advect_ratio`) pair.  Both must be >= 1 and
      !! `dt_therm_ratio` must be an integer multiple of
      !! `dt_tracer_advect_ratio` so the ALE remap (which fires at the
      !! DT_THERM cadence) never lands inside an open tracer-flux
      !! accumulation window.  The setup layer (`configure_ocean_vmix`)
      !! calls this and `error stop`s with a descriptive message on
      !! `.false.`; exposed as a pure predicate so the validation logic
      !! is unit-testable without constructing a full ocean state.
      integer, intent(in) :: dt_therm_ratio
      integer, intent(in) :: dt_tracer_advect_ratio
      logical :: ok
      ok = .true.
      if (dt_tracer_advect_ratio < 1) ok = .false.
      if (dt_therm_ratio < 1) ok = .false.
      if (ok) then
         if (mod(dt_therm_ratio, dt_tracer_advect_ratio) /= 0) ok = .false.
      end if
   end function ocean_dt_tracer_advect_ratios_ok

   subroutine ocean_dyn_flush_tracer_window(grid, metrics, dyn, ct, ms, bc)
      !! Mandatory end-of-segment flush of the windowed tracer-advect
      !! accumulators (spec §(c) — MOM6's `n == n_max`).  Drains any OPEN
      !! accumulation window (`t_dyn_rel_adv > 0`) so no Lagrangian-advanced
      !! `h_layer` is ever paired with FROZEN `hTr` at an output write, a
      !! restart checkpoint, or the end of a run segment whose length is not
      !! an exact multiple of `dt_tracer_advect_ratio`.
      !!
      !! Hard no-op at `dt_tracer_advect_ratio <= 1` (ratio = 1 never
      !! accumulates ⇒ `t_dyn_rel_adv` stays 0 ⇒ bit-identical).  Idempotent:
      !! the drain resets the accumulators + clock, so a redundant flush over
      !! an already-empty window reconstructs `hprev = h_end` (div of the
      !! zeroed accumulators) and drains zero transport — an exact no-op.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(in) :: dyn
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_bc_state_t), intent(in), optional :: bc

      if (dyn%dt_tracer_advect_ratio <= 1) return
      if (ct%t_dyn_rel_adv <= 0.0_wp) return

      if (present(bc)) then
         call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio, bc=bc)
      else
         call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio)
      end if
   end subroutine ocean_dyn_flush_tracer_window

   pure function isopycnal_vanish_tol(angstrom_h, pd_floor) result(tol)
      !! Shared vanish-tolerance for Phase-2/3 kernels: the layer is
      !! considered vanished when its thickness is at or below this value.
      !! Defined as max(angstrom_h, H_VANISHED) so both the physical
      !! floor (Phase 1) and the skip/merge marker are covered by a
      !! single threshold.  Providing this as a `pure` helper ensures
      !! Phase 2 (`reset_vanished_layer_velocities`) and Phase 3
      !! (`apply_velocity_truncation` gate + `compute_max_cfl` gate)
      !! cannot drift in their vanish definition.
      !!
      !! **`pd_floor` — the floor/tolerance collision fix (pdc 2c7ad0a1).**
      !! Positive-definite continuity holds every floored layer at
      !! `h >= h_lim = angstrom_h`, so a floored (dynamically dead) layer
      !! sits AT or just-above the plain `max(angstrom_h, H_VANISHED) =
      !! angstrom_h` tolerance — and the `h <= tol` guards then classify it
      !! LIVE, leaving phantom maxvel velocities that the CFL clip can't
      !! reach and that pin the panic MaxCFL (600² spoon, day-16 2.4696).
      !! When `pd_floor = .true.` (PD active, `angstrom_h > H_VANISHED`)
      !! the tolerance is lifted to `angstrom_h + H_VANISHED` — strictly
      !! ABOVE the floor it now rests on — so the at-floor band
      !! `[h_lim, h_lim + H_VANISHED]` reads as vanished.  `pd_floor`
      !! absent / .false. ⇒ unchanged ⇒ bit-identical.
      real(wp), intent(in) :: angstrom_h
         !! Phase-1 floor (m); 0.0 when Phase 1 is off.
      logical, intent(in), optional :: pd_floor
         !! When .true. AND `angstrom_h > H_VANISHED`: lift the tolerance
         !! to `angstrom_h + H_VANISHED` (above the PD floor).
      real(wp) :: tol
      logical :: pd
      pd = .false.
      if (present(pd_floor)) pd = pd_floor
      if (pd .and. angstrom_h > H_VANISHED) then
         tol = angstrom_h + H_VANISHED
      else
         tol = max(angstrom_h, H_VANISHED)
      end if
   end function isopycnal_vanish_tol

   subroutine ocean_dyn_destroy(this)
      class(ocean_dyn_t), intent(inout) :: this
      this%is_init = .false.
      call this%bt_work%destroy()
      if (allocated(this%bt_wide)) then
         call this%bt_wide%destroy()
         deallocate (this%bt_wide)
      end if
      if (allocated(this%avr_u0)) deallocate (this%avr_u0)
      if (allocated(this%avr_v0)) deallocate (this%avr_v0)
   end subroutine ocean_dyn_destroy

   subroutine ocean_dyn_enter_data(this)
      class(ocean_dyn_t), intent(inout) :: this
      select type (this)
      type is (ocean_dyn_t)
         call ocean_dyn_enter_data_impl(this)
      end select
   end subroutine ocean_dyn_enter_data

   subroutine ocean_dyn_enter_data_impl(this)
      type(ocean_dyn_t), intent(inout) :: this
      call barotropic_workstate_enter_data_impl(this%bt_work)
      ! accel_visc_rem stage-entry snapshots: eagerly allocated at setup
      ! (host source=0) only when the knob is on, so `copyin` (not `create`)
      ! and gate on allocation — symmetric with the setup allocation and the
      ! exit_data delete.
      if (allocated(this%avr_u0)) then
         !$acc enter data copyin(this%avr_u0, this%avr_v0)
      end if
      ! bt_wide is attached by ocean_dyn_enable_bt_wide (AFTER enter_data
      ! for the rest of dyn); no attach here — see ocean_dyn_enable_bt_wide.
   end subroutine ocean_dyn_enter_data_impl

   subroutine ocean_dyn_exit_data(this)
      class(ocean_dyn_t), intent(inout) :: this
      select type (this)
      type is (ocean_dyn_t)
         call ocean_dyn_exit_data_impl(this)
      end select
   end subroutine ocean_dyn_exit_data

   subroutine ocean_dyn_exit_data_impl(this)
      type(ocean_dyn_t), intent(inout) :: this
      if (allocated(this%bt_wide)) call this%bt_wide%exit_data()
      ! avr_* are eagerly allocated + mapped at setup (only when the
      ! accel_visc_rem knob is on), so guard the delete on allocation —
      ! symmetric with the enter_data copyin.
      if (allocated(this%avr_u0)) then
         !$acc exit data delete(this%avr_u0, this%avr_v0)
      end if
      call barotropic_workstate_exit_data_impl(this%bt_work)
   end subroutine ocean_dyn_exit_data_impl

   subroutine ocean_dyn_enable_bt_wide(dyn, grid, dx, dy, lon_west, lat_south, &
                                       rad_earth, grid_config_str, &
                                       f_0, beta, y_ref, coriolis_scheme_str)
      !! Allocate, initialise, and GPU-attach the wide-halo shadow state
      !! from `dyn%bt_halo` (already set by the caller).
      !!
      !! Preconditions:
      !!   - `dyn%init` has been called.
      !!   - `dyn%enter_data` has been called (bt_work is on the GPU; we add
      !!     bt_wide to the same GPU context).
      !!   - `dyn%bt_halo > 0`.
      !!
      !! Odd bt_halo is rounded DOWN to even (with a warning).
      !! ng_wide = nghost + bt_halo must be <= min(nx_phys, ny_phys)
      !! (fail-loud: the wide halo must fit inside the local domain).
      type(ocean_dyn_t), intent(inout) :: dyn
         !! Dynamics state; `bt_halo` must be set before calling.
      type(hgrid_t), intent(in) :: grid
         !! Normal-width grid for this subdomain.
      real(wp), intent(in) :: dx, dy
         !! Cell spacing (m for Cartesian; deg for spherical).
      real(wp), intent(in) :: lon_west, lat_south, rad_earth
         !! Spherical-grid parameters (ignored for Cartesian).
      character(len=*), intent(in) :: grid_config_str
         !! Grid-config string (e.g. "cartesian", "spherical").
      real(wp), intent(in) :: f_0, beta, y_ref
         !! Beta-plane Coriolis parameters.
      character(len=*), intent(in) :: coriolis_scheme_str
         !! Coriolis-scheme string (e.g. "beta_plane").

      integer :: bh, ng_w, grid_cfg, cor_scheme

      bh = dyn%bt_halo
      if (bh <= 0) return

      ! Round odd bt_halo down to even.
      if (mod(bh, 2) /= 0) then
         call logger%warning("ocean_dyn_enable_bt_wide: bt_halo="//to_string(bh)// &
                             " is odd; rounding down to "//to_string(bh - 1))
         bh = bh - 1
         dyn%bt_halo = bh
      end if

      ng_w = grid%nghost + bh
      if (ng_w > min(grid%nx_phys, grid%ny_phys)) then
         error stop "ocean_dyn_enable_bt_wide: bt_halo too large — ng_wide exceeds "// &
            "min(nx_phys,ny_phys); reduce bt_halo or increase the domain"
      end if

      grid_cfg = parse_grid_config(grid_config_str)
      cor_scheme = parse_coriolis_scheme(coriolis_scheme_str)

      allocate (dyn%bt_wide)
      dyn%bt_wide%bt_halo = bh
      call dyn%bt_wide%init(grid, dx, dy, lon_west, lat_south, rad_earth, &
                            grid_cfg, f_0, beta, y_ref, cor_scheme)
      call dyn%bt_wide%enter_data()
   end subroutine ocean_dyn_enable_bt_wide

   pure subroutine ocean_dyn_step_barotropic(grid, metrics, dyn, cor, ct, bs, dt)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Unsplit SSP-RK2 (Heun's method) outer step on the barotropic
      !! C-grid state.  Couples continuity-PPM (h-update) and the
      !! Sadourny Coriolis-advection tendency (u, v update) into one
      !! second-order-accurate step.
      !!
      !! SSP-RK2 (Shu-Osher form):
      !!   u^(1)    = u^n + dt * L(u^n)
      !!   u^(n+1)  = 1/2 * u^n + 1/2 * (u^(1) + dt * L(u^(1)))
      !!
      !! In our code:
      !!   1. Save u^n into the *_0 buffers on `bs`.
      !!   2. Stage 1: compute continuity flux divergence and the
      !!      Sadourny momentum tendency at u^n, apply both as a
      !!      forward-Euler step.  State is now u^(1).
      !!   3. Stage 2: compute tendencies at u^(1), apply another
      !!      forward-Euler step.  State is u^(1) + dt * L(u^(1)).
      !!   4. RK2 average: state <- 1/2 * (u_0 + state).
      !!
      !! Errors per step are O(dt^3); compared to plain FE, the
      !! inertial-oscillator amplitude growth drops from O(dt^2) per
      !! step to O(dt^4) per step.  At f*dt = 0.01 over 157 steps the
      !! magnitude error drops from ~0.8% (FE) to ~2e-7 (RK2).
      !!
      !! `dyn` carries diagnostic state (step counter, CFL, KE) for
      !! the outer step; Phase 4a only bumps `outer_step_count`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      type(continuity_t), intent(inout) :: ct
      type(barotropic_state_t), intent(inout) :: bs
      real(wp), intent(in) :: dt

      integer :: i, j, nx, ny, nx_face, ny_uface, nx_vface, ny_face

      nx = grid%nx_total
      ny = grid%ny_total
      nx_face = size(bs%u_face_x, 1)
      ny_uface = size(bs%u_face_x, 2)
      nx_vface = size(bs%v_face_y, 1)
      ny_face = size(bs%v_face_y, 2)

      ! ---- 1. Save u^n into h0 / u_face_x0 / v_face_y0 ----
      do concurrent(j=1:ny, i=1:nx)
         bs%h0(i, j) = bs%h(i, j)
      end do
      do concurrent(j=1:ny_uface, i=1:nx_face)
         bs%u_face_x0(i, j) = bs%u_face_x(i, j)
      end do
      do concurrent(j=1:ny_face, i=1:nx_vface)
         bs%v_face_y0(i, j) = bs%v_face_y(i, j)
      end do

      ! ---- 2. Stage 1: tendencies at u^n, FE step -> u^(1) ----
      call continuity_compute_fluxes_barotropic(grid, metrics, ct, bs)
      call coriolis_adv_compute_tendencies_barotropic(grid, metrics, cor, bs)
      call continuity_apply_fluxes_barotropic(bs, dt)
      call coriolis_adv_apply_tendencies_barotropic(cor, bs, dt)

      ! ---- 3. Stage 2: tendencies at u^(1), FE step -> u^(1) + dt*L(u^(1)) ----
      call continuity_compute_fluxes_barotropic(grid, metrics, ct, bs)
      call coriolis_adv_compute_tendencies_barotropic(grid, metrics, cor, bs)
      call continuity_apply_fluxes_barotropic(bs, dt)
      call coriolis_adv_apply_tendencies_barotropic(cor, bs, dt)

      ! ---- 4. RK2 average: u^(n+1) = 1/2 * (u^n + (u^(1) + dt*L(u^(1)))) ----
      do concurrent(j=1:ny, i=1:nx)
         bs%h(i, j) = 0.5_wp*(bs%h0(i, j) + bs%h(i, j))
      end do
      do concurrent(j=1:ny_uface, i=1:nx_face)
         bs%u_face_x(i, j) = 0.5_wp*(bs%u_face_x0(i, j) + bs%u_face_x(i, j))
      end do
      do concurrent(j=1:ny_face, i=1:nx_vface)
         bs%v_face_y(i, j) = 0.5_wp*(bs%v_face_y0(i, j) + bs%v_face_y(i, j))
      end do

      dyn%outer_step_count = dyn%outer_step_count + 1
   end subroutine ocean_dyn_step_barotropic

   subroutine ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, dt, sf, geo, lateral_mix, epbl, kshear, slopes, vmix_tidal, bc, t)
      !! Multilayer extension of `ocean_dyn_step_barotropic`.  One
      !! SSP-RK2 outer step that orchestrates the full per-layer
      !! dynamical core:
      !!
      !!   EOS  → continuity-PPM  → Sadourny Coriolis-adv
      !!        → Mont pressure-force  → Laplacian hvisc
      !!        → bottom drag  → surface stress
      !!        → horizontal tracer advection
      !!        → diagnose w  → vertical tracer advection (Eulerian z)
      !!        → applies
      !!
      !! Each stage runs every compute kernel in sequence, then the
      !! apply step.  The two SSP-RK2 stages save the prognostic
      !! state into the *_0 buffers (h_layer0, u_face_x_layer0,
      !! v_face_y_layer0, plus hTr0 on every registered tracer),
      !! run two full FE substeps, and average with the saved
      !! state to recover the second-order-accurate result.
      !!
      !! The Coriolis parameter lives on the `coriolis_adv_t` slot
      !! as `f_corner(:, :)` (C-grid corners).  Defaults to a uniform
      !! `f_0` at init; call `cor%set_beta_plane(grid, f_0, beta, y_ref)`
      !! for the beta-plane variant `f(y) = f_0 + beta*(y - y_ref)`.
      !!
      !! Caller passes the individual slots rather than the full
      !! `ocean_state_t` to avoid a circular module dependency
      !! (rdb_ocean_state already `use`s rdb_ocean_dyn).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics — forwarded to the PGF (and,
         !! in slice 2, the other geometry-aware kernels).
      type(ocean_dyn_t), intent(inout) :: dyn
      type(eos_t), intent(in) :: eos
      type(coriolis_adv_t), intent(inout) :: cor
      type(continuity_t), intent(inout) :: ct
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      type(ocean_surface_flux_t), intent(in), optional :: sf
      type(ocean_geothermal_t), intent(in), optional :: geo
         !! Geothermal bottom-heat-flux slot.  Absent or `enable=.false.`
         !! preserves the historical no-geothermal path bit-identically.
      type(ocean_lateral_mix_t), intent(inout), optional :: lateral_mix
      type(ocean_epbl_t), intent(inout), optional :: epbl
         !! Energetics-based PBL slot.  Absent or `enable=.false.`
         !! preserves the historical PP81/KPP path bit-identically.
      type(ocean_kappa_shear_t), intent(inout), optional :: kshear
         !! Kappa-shear interior closure slot.  Absent or
         !! `enable=.false.` preserves the historical path bit-identically.
      type(ocean_slopes_t), intent(inout), optional :: slopes
         !! Isopycnal-slope diagnostic slot.  Absent or `enable=.false.`
         !! preserves the historical path bit-identically (diagnostic
         !! only — never feeds back into prognostics).
      type(ocean_tidal_mixing_t), intent(inout), optional :: vmix_tidal
         !! Tidal-mixing interior closure slot.  Absent or
         !! `enable=.false.` preserves the historical path bit-identically.
      type(ocean_bc_state_t), intent(in), optional :: bc
         !! Boundary state — forwarded ONLY to the windowed tracer-advect
         !! drain so its single-rank periodic-x / north-fold seam wraps
         !! fire for periodic unsplit runs at `dt_tracer_advect_ratio > 1`.
         !! Absent ⇒ wall seams (the historical unsplit default); ratio = 1
         !! never drains so this is inert there (bit-identical).
      real(wp), intent(in), optional :: t
         !! Model time (s) since run start, used to evaluate the ideal-age
         !! vintage-mode surface value (PR-7).  Mirrors the split driver's
         !! `t` dummy (`ocean_dyn_step_split`).  Absent ⇒ `t = 0`.

      integer :: it
      real(wp) :: t_now

      ! ---- Save u^n into the *_0 buffers ----
      call save_state(ms)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            call copy_field_3d(ms%tracers(it)%hTr, ms%tracers(it)%hTr0, &
                               size(ms%tracers(it)%hTr, 1), &
                               size(ms%tracers(it)%hTr, 2), &
                               size(ms%tracers(it)%hTr, 3))
         end do
      end if

      ! ---- Stage 1: tendencies at u^n, FE step -> u^(1) ----
      call run_stage(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, &
                     vd, vmix, ms, dt, 1, sf=sf, geo=geo, lateral_mix=lateral_mix, &
                     epbl=epbl, kshear=kshear, slopes=slopes, vmix_tidal=vmix_tidal)

      ! ---- Stage 2: tendencies at u^(1), FE step -> u^(1) + dt*L(u^(1)) ----
      call run_stage(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, &
                     vd, vmix, ms, dt, 2, sf=sf, geo=geo, lateral_mix=lateral_mix, &
                     epbl=epbl, kshear=kshear, slopes=slopes, vmix_tidal=vmix_tidal)

      ! ---- RK2 average: u^(n+1) = 0.5 * (u^n + stage2 result) ----
      call rk2_average(ms)
      if (dyn%check_h_positive) then
         call check_h_positive_or_die(grid, ms, "after rk2_average", 3, &
                                      dyn%outer_step_count + 1, check_layers=.true.)
      end if
      ! Land-face velocity reset on the averaged state (C4 / R4b.2).
      call mask_layer_velocities(grid, metrics, ms, bt_work=dyn%bt_work)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            call rk2_average_field_3d(ms%tracers(it)%hTr0, ms%tracers(it)%hTr, &
                                      size(ms%tracers(it)%hTr, 1), &
                                      size(ms%tracers(it)%hTr, 2), &
                                      size(ms%tracers(it)%hTr, 3))
         end do
      end if

      ! Velocity housekeeping (E7): advective-CFL truncation then the
      ! absolute maxvel cap.  Both no-op when their knob <= 0.
      call apply_velocity_truncation(ms, metrics, dt, dyn%cfl_trunc, dyn%maxvel, dyn%ntrunc_step)
      dyn%ntrunc_total = dyn%ntrunc_total + dyn%ntrunc_step

      ! Phase 2 (6b) windowed tracer-advect drain (unsplit reference path).
      ! No-op at ratio = 1 (never accumulated).  No vcoord here (the unsplit
      ! path does not remap), so the drain fires on the window-full predicate
      ! only.  `bc` IS forwarded (when present) so periodic-x / north-fold
      ! seam wraps fire for periodic unsplit runs — matching the split path.
      if (dyn%dt_tracer_advect_ratio > 1) then
         ct%t_dyn_rel_adv = ct%t_dyn_rel_adv + dt
         if (dyn%is_tracer_advect_step()) then
            if (present(bc)) then
               call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio, bc=bc)
            else
               call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio)
            end if
         end if
      end if

      ! Ideal-age surface reset (PR-7): once per outer step, after
      ! rk2_average and the windowed-drain block above (last operator to
      ! touch k=nz — there is no ALE remap on the unsplit path, so
      ! post-average + post-drain is already last).  See the ordering
      ! comment in `ocean_dyn_step_split` for the full rationale.
      if (dyn%is_thermo_step()) then
         t_now = 0.0_wp
         if (present(t)) t_now = t
         call ocean_ideal_age_reset_surface(grid, ms, &
                                            ocean_ideal_age_young_val(dyn%ideal_age_young_val, &
                                                                      dyn%ideal_age_sfc_growth_rate, t_now))
      end if

      dyn%outer_step_count = dyn%outer_step_count + 1
   end subroutine ocean_dyn_step

   subroutine run_stage(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, dt, stage, sf, geo, lateral_mix, epbl, kshear, slopes, vmix_tidal)
      !! One FE stage of the multilayer step.  Order of operations:
      !!
      !!   1. EOS: rho_layer <- linear(T, S)
      !!   2. Velocity-tendency computes — Coriolis, PGF, hvisc,
      !!      bottom drag, surface stress.  All read h_old, u, v
      !!      and write to their own tendency buffers; none touch
      !!      h yet.
      !!   3. continuity_tracer_step_split — interleaved
      !!      Lie-split horizontal step: zonal flux + tracer zonal +
      !!      apply zonal + meridional flux + tracer meridional +
      !!      apply meridional.  Updates h AND hTr together; the
      !!      CWC discrete theorem holds (uniform T stays uniform).
      !!   4. tracer_hdiff — Laplacian diffusion on hTr.
      !!   5. compute_w_from_continuity — diagnose w from the total
      !!      horizontal divergence stored in `flux_h_layer`.
      !!   6. tracer_advect_vertical — upwind-in-z + h update that
      !!      cancels the horizontal h change (Eulerian-z mode).
      !!   7. Velocity-tendency applies (Coriolis, PGF, hvisc, drag,
      !!      surface stress) — all additive.
      !!   8. Surface tracer fluxes (heat, salt).
      !!   9. Vertical mixing closure + vdiff.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics — forwarded to the PGF.
      type(ocean_dyn_t), intent(in) :: dyn
      type(eos_t), intent(in) :: eos
      type(coriolis_adv_t), intent(inout) :: cor
      type(continuity_t), intent(inout) :: ct
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      integer, intent(in) :: stage
         !! RK2 stage (1 or 2); see vmix_apply_in_stage.
      type(ocean_surface_flux_t), intent(in), optional :: sf
      type(ocean_geothermal_t), intent(in), optional :: geo
         !! Geothermal bottom-heat-flux slot.  See `ocean_dyn_step`.
      type(ocean_lateral_mix_t), intent(inout), optional :: lateral_mix
      type(ocean_epbl_t), intent(inout), optional :: epbl
      type(ocean_kappa_shear_t), intent(inout), optional :: kshear
      type(ocean_slopes_t), intent(inout), optional :: slopes
      type(ocean_tidal_mixing_t), intent(inout), optional :: vmix_tidal

      logical :: therm_active
      real(wp) :: therm_dt

      therm_active = dyn%enable_thermodynamics .and. dyn%is_thermo_step()
      therm_dt = dyn%therm_dt(dt)

      ! EOS is DYNAMICS, not thermodynamics: rho_layer feeds the baroclinic
      ! PGF every step. Gating it on `is_thermo_step()` zero-order-holds the
      ! restoring force of the internal-gravity-wave oscillator for
      ! tau = (dt_therm_ratio-1)*dt, which is a delayed-restoring-force
      ! instability: sigma ~ omega^2*tau/2, unstable for ANY tau > 0, maximal
      ! at the grid scale. Measured on eady: sigma ∝ alpha_T ∝ N^2, sigma ∝ k^2,
      ! the mode is the diagonal 2-delta checkerboard, and the rate recovers
      ! the mode-2 internal wave speed. Survival was only ever by viscosity
      ! margin (nu*k^2 > sigma) — eady at nu_h=100 is marginal and blows up;
      ! double_gyre at nu_h=10000 merely looks fine.
      !
      ! So recompute rho EVERY dynamics step whenever thermodynamics is on.
      ! The EXPENSIVE thermo (mixing, ALE remap, tracer advection/hdiff,
      ! surface fluxes) stays on the slow `therm_active` cadence below, which
      ! is where the cost actually is. Bit-identical at dt_therm_ratio = 1,
      ! where is_thermo_step() is always true.
      call ocean_eos_compute(eos, ms, active=dyn%enable_thermodynamics)
      call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      call ocean_lateral_mix_compute(grid, metrics, lateral_mix, ms)
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, &
                                                         lateral_mix=lateral_mix, dt=dt)
      ! KE dissipation rate for the MEKE frictional source — captured here
      ! (du_visc fresh, u_face still pre-viscous) before the apply below.
      ! No-op unless MEKE's frictional source is enabled.
      call ocean_horizontal_viscosity_compute_ke_diss(hv, ms)
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, dt)
      call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
      call ocean_surface_stress_compute_tendencies(grid, ss, ms)

      ! Horizontal step + tracer chain.  Continuity-tracer is
      ! unconditional (the h advection lives here regardless of
      ! thermodynamics); the tracer-only kernels self-gate on
      ! `therm_active`.
      !
      ! Phase 2 (6b) horizontal-tracer-advect cadence dispatch
      ! (`dt_tracer_advect_ratio`):
      !   ratio == 1 (default) → the existing every-step fused
      !     continuity+tracer kernel runs verbatim ⇒ bit-identical.
      !   ratio  > 1 → TR_MODE_ACCUMULATE: advance h + accumulate
      !     0.5·mass_flux·dt into ct%uhtr/vhtr per RK2 stage (hTr frozen);
      !     the boundary drain in `ocean_dyn_step` spends them.
      if (dyn%dt_tracer_advect_ratio <= 1) then
         call continuity_tracer_step_split(grid, metrics, ct, ms, dt)
      else
         call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                           tracer_mode=TR_MODE_ACCUMULATE)
      end if
      call tracer_hdiff(grid, metrics, hd, ms, therm_dt, active=therm_active)
      ! Mass budget: flux_h_layer now holds the total horizontal divergence
      ! (the same field compute_w_from_continuity reads).  Accumulate the
      ! boundary mass outflux for this RK2 stage (weight 0.5, since rk2_average
      ! halves each stage's contribution); vertical advection + ALE remap only
      ! redistribute within a column, so they don't affect the column-mass
      ! budget.  Closes the ocean mass Error to round-off with open BCs.
      call ocean_accumulate_mass_out(ms, ms%flux_h_layer, metrics%areaT, &
                                     grid%nghost, dt, 0.5_wp)
      call compute_w_from_continuity(grid, va, ms)
      call tracer_advect_vertical(grid, va, ms, therm_dt, active=therm_active)

      ! Velocity-side applies (no thermodynamics gate).  Async-chained on
      ! OpenACC queue 1 (additive accumulations onto u_face/v_face, FIFO-
      ! ordered).  The interleaved tracer applies below run on the default
      ! queue but touch disjoint arrays (hTr), so there is no cross-queue
      ! race.  ONE `!$acc wait(1)` before vmix — the first device consumer
      ! that reads the applied u_face/v_face (shear).
      call coriolis_adv_apply_tendencies(cor, ms, dt, no_wait=.true.)
      call ocean_pressure_force_apply(pgf, ms, dt, no_wait=.true.)
      call ocean_horizontal_viscosity_apply_tendencies(hv, ms, dt, no_wait=.true.)
      ! Double-count guard: when the bottom drag / wind stress are folded
      ! into the implicit vdiff tridiagonal (`&ocean_vdiff_nml implicit_*`),
      ! their explicit pre-solve applies are SKIPPED here so the forcing is
      ! not applied twice.  Channel (side-wall) drag is a distinct lateral
      ! term and always applies.  Defaults (both off) ⇒ both applies run ⇒
      ! bit-identical to the prior path.
      if (.not. vd%implicit_drag) then
         call ocean_bottom_drag_apply_tendencies(bd, ms, dt, no_wait=.true.)
      end if
      call ocean_channel_drag_apply_tendencies(bd, ms, dt, no_wait=.true.)
      if (.not. vd%implicit_stress) then
         call ocean_surface_stress_apply_tendencies(ss, ms, dt, no_wait=.true.)
      end if

      ! Surface tracer fluxes, geothermal bottom flux, ideal-age, vmix.
      call ocean_surface_flux_apply_tracers(grid, sf, ms, therm_dt, active=therm_active)
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, therm_dt, active=therm_active)
      call ocean_surface_restore_apply_tracers(grid, sf, ms, therm_dt, active=therm_active)
      call ocean_geothermal_apply_tracers(grid, geo, ms, therm_dt, active=therm_active)
      ! Ideal-age interior aging only (PR-7): thermo-cadence gated, mirrors
      ! its tracer-kernel neighbours above.  The surface Dirichlet reset is
      ! NOT here — it runs once per outer step, after rk2_average, in
      ! `ocean_dyn_step` (a per-stage reset is halved by the RK2 average).
      call ocean_ideal_age_apply(grid, ms, therm_dt, active=dyn%is_thermo_step())

      ! Isopycnal-slope diagnostics — purely diagnostic, refreshed at
      ! thermo cadence (same gate as the lateral closures it feeds).
      ! No-op when absent / disabled (bit-identical).
      if (present(slopes) .and. therm_active) then
         if (slopes%enable) then
            call ocean_slopes_compute(grid, metrics, eos, slopes, ms, therm_dt)
         end if
      end if

      !$acc wait(1)
      ! No `bt_work` here: the unsplit path has no BT correction at all,
      ! so there is no consumer for visc_rem — do_remnant stays .false.
      call vmix_apply_in_stage(grid, dyn, vmix, vd, ss, bd, ms, dt, stage, sf, epbl=epbl, kshear=kshear, &
                               vmix_tidal=vmix_tidal, metrics=metrics)
   end subroutine run_stage

   subroutine vmix_apply_in_stage(grid, dyn, vmix, vd, ss, bd, ms, dt, stage, sf, epbl, kshear, vmix_tidal, bt_work, &
                                  apply_tracers, metrics)
      !! Bundle the per-stage vmix closure / KPP overlay / KV_ML_INVZ2 /
      !! assembly gate / vdiff dispatch into one routine so the run_stage
      !! drivers can call `vmix_apply_in_stage(grid, dyn, vmix, vd, ss,
      !! ms, dt, sf)` instead of carrying 30 lines of nested if-branching.
      !!
      !! Closure chain (all upstream CONTRIBUTORS into kv/kt):
      !!   PP81 interior → KPP overlay XOR EPBL merge → kappa-shear
      !!   additive merge → tidal-mixing additive merge → KV_ML_INVZ2
      !!   surface band → convective adjustment (Brunt-Vaisala trigger,
      !!   interior-only, masked below the active KPP/EPBL boundary
      !!   layer) → **vmix_split_kd_heat_salt** (derives ks from kt; NOT a
      !!   contributor, must stay last) → **vmix_assemble** (the single
      !!   downstream gate: background floors, kv_max/kd_max ceilings,
      !!   optional smoothing, optional guard — applies to kv, kt, AND ks)
      !!   → vdiff.
      !!
      !! Logic preserved verbatim from the prior in-driver dispatch:
      !!   - When `vmix%use_closure` is true: optional PP81 / KPP
      !!     overlay populate `vmix%kv`/`vmix%kt`, then vdiff reads
      !!     them; KPP non-local γ is applied after the tracer vdiff
      !!     solve when thermodynamics are active.
      !!   - When `vmix%use_closure` is false: vdiff uses its scalar
      !!     `K_v_*` defaults — kv_source not passed.
      !!   - Tracer vdiff + KPP non-local fire only when
      !!     `enable_thermodynamics .and. is_thermo_step()`.
      !!
      !! Profiler regions are emitted from here uniformly (the prior
      !! run_stage path had no profiler regions around vmix; the
      !! split path did — now both share the same labels).
      type(hgrid_t), intent(in) :: grid
      type(ocean_dyn_t), intent(in) :: dyn
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_surface_stress_t), intent(in) :: ss
      type(ocean_bottom_drag_t), intent(in) :: bd
         !! Bottom-drag slot — supplies the bed-layer Rayleigh-rate field
         !! `lambda_bot_u/v` for the implicit-drag vdiff fold.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      integer, intent(in) :: stage
         !! RK2 stage (1 or 2): the EPBL/kappa-shear COLUMN SOLVES fire
         !! only at stage 1 of a thermo step (their kd is
         !! stage-invariant by design — recomputing at stage 2 doubled
         !! the closure cost for no accuracy: 19.6%% of GPU time was
         !! kappa-shear at 2x cadence).  The merges still run EVERY
         !! stage (PP81 rewrites kv/kt per stage).
      type(ocean_surface_flux_t), intent(in), optional :: sf
      type(ocean_epbl_t), intent(inout), optional :: epbl
         !! EPBL slot.  When present and enabled, replaces the KPP
         !! overlay (configure enforces the mutual exclusion):
         !! `epbl_compute` refreshes `kd_int` at thermo cadence and
         !! the merge into `vmix%kv` / `vmix%kt` runs every stage
         !! (PP81 rewrites those arrays each stage).
      type(ocean_kappa_shear_t), intent(inout), optional :: kshear
         !! Kappa-shear interior closure slot.  When present and enabled,
         !! `kappa_shear_compute` refreshes `kd_int` at thermo cadence and
         !! the additive merge into `vmix%kv` / `vmix%kt` runs every stage.
         !! Coexists with KPP / EPBL (no mutual exclusion).
      type(ocean_tidal_mixing_t), intent(inout), optional :: vmix_tidal
         !! St-Laurent/Simmons tidal-mixing interior closure slot.  When
         !! present and enabled, `tidal_mixing_compute` refreshes `kd_int`
         !! at thermo cadence and the additive merge into `vmix%kv` /
         !! `vmix%kt` runs every stage.  Coexists with KPP / EPBL /
         !! kappa-shear (no mutual exclusion).
      type(barotropic_workstate_t), intent(inout), optional :: bt_work
         !! The BT-corrector workstate — supplies `visc_rem_u/v` as the
         !! vdiff kernel's OUTPUT.  Present only from `run_stage_split`
         !! (the unsplit path has no BT correction to consume it).  When
         !! present AND `bt_work%bt_correction_visc_rem`, the viscous
         !! remnant γ is (re)computed here, at step 9 of the CURRENT
         !! stage — the BT corrector at step 7 of the NEXT stage reads
         !! it, a one-stage (Δt/2) lag (see the step-9 call site below).
         !! Absent, or the knob off, ⇒ no remnant work ⇒ bit-identical.

      logical, intent(in), optional :: apply_tracers
         !! `.false.` = momentum-only: skip the tracer vdiff + KPP
         !! non-local applies regardless of the thermo gate.  The pred_corr
         !! PREDICTOR passes this — MOM6's predictor applies
         !! `vertvisc(up, dt_pred)` to the provisional velocity (so the
         !! spurious grounded-layer accelerations are absorbed BEFORE
         !! continuity forms `u_av`) but never touches tracers.
         !! Default `.true.` = historical.
      type(ocean_metrics_t), intent(in), optional :: metrics
         !! Metrics slot — supplies the halo-valid wet masks
         !! (`wet_T`/`wet_u`/`wet_v`) the kappa-shear VERTEX form's
         !! corner gather needs, and `metrics%geolatT` — the C7 Henyey
         !! latitude factor's only spatial input, forwarded to
         !! `vmix_assemble` (unread when `bkgnd_henyey` is off).
         !! Optional so the split_rk2/legacy call shapes stay valid; both
         !! consumers fail loud when their knob is on and the slot is
         !! absent — `kappa_shear_compute` for `at_vertex`, and
         !! `vmix_assemble` for `bkgnd_henyey`.

      logical :: epbl_active, kshear_active, tidal_active, do_remnant
      logical :: do_tracers, vertex_kv

      do_tracers = .true.
      if (present(apply_tracers)) do_tracers = apply_tracers
      epbl_active = .false.
      if (present(epbl)) epbl_active = epbl%enable
      kshear_active = .false.
      if (present(kshear)) kshear_active = kshear%enable
      ! Vertex kappa-shear routes its viscosity corner->face: the
      ! momentum vdiff gets `kd_corner` as `kv_corner_source` (scaled by
      ! prandtl_turb) and the cell-centred kv merge is suppressed inside
      ! `kappa_shear_merge_into_kv_kt` (MOM6 zeroes the tracer-point
      ! Kv_shear when VERTEX_SHEAR is on — no double-count).
      vertex_kv = .false.
      if (present(kshear)) vertex_kv = kshear_active .and. kshear%at_vertex
      tidal_active = .false.
      if (present(vmix_tidal)) tidal_active = vmix_tidal%enable
      do_remnant = .false.
      if (present(bt_work)) do_remnant = bt_work%bt_correction_visc_rem

      if (vmix%use_closure) then
         call profiler_start("ocean_vmix_compute")
         if (vmix%interior_closure == VMIX_INTERIOR_PP81) then
            call vmix_compute_pp81(grid, vmix, ms)
         else
            ! PR-6 fail-loud: VMIX_INTERIOR_LARGE94 / VMIX_INTERIOR_CVMIX are
            ! reserved-but-unwired (no kernel).  Selecting one would leave
            ! kv/kt with NO interior mixing.  There is no namelist key for
            ! interior_closure today, so this is defence-in-depth for the
            ! next code/config consumer that sets it.
            call logger%error("vmix_apply_in_stage: the selected interior closure "// &
                              "has no kernel (only VMIX_INTERIOR_PP81 is implemented; "// &
                              "VMIX_INTERIOR_LARGE94/CVMIX are reserved-but-unwired) — "// &
                              "running it would leave kv/kt with no interior mixing")
            error stop "vmix_apply_in_stage: unimplemented vmix interior closure"
         end if
         if (vmix%use_kpp .and. .not. epbl_active) then
            if (.not. present(sf)) then
               call logger%error("vmix_apply_in_stage: KPP requires the surface-flux slot (sf)")
               error stop "vmix_apply_in_stage: use_kpp requires the surface-flux slot (sf)"
            end if
            call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
         end if
         if (epbl_active) then
            if (.not. present(sf)) then
               call logger%error("vmix_apply_in_stage: EPBL requires the surface-flux slot (sf)")
               error stop "vmix_apply_in_stage: EPBL requires the surface-flux slot (sf)"
            end if
            if (dyn%enable_thermodynamics .and. dyn%is_thermo_step() .and. stage == 1) then
               call epbl_compute(grid, epbl, ms, ss, dyn%therm_dt(dt), sf)
            end if
            call epbl_merge_into_kv_kt(epbl, size(vmix%kt, 1), size(vmix%kt, 2), &
                                       size(vmix%kt, 3), vmix%kv, vmix%kt)
         end if
         if (kshear_active) then
            if (dyn%enable_thermodynamics .and. dyn%is_thermo_step() .and. stage == 1) then
               if (kshear%at_vertex) then
                  ! Vertex form needs the halo-valid wet masks for the
                  ! corner gather; fail-loud inside if metrics is absent.
                  if (.not. present(metrics)) then
                     call logger%error("vmix_apply_in_stage: kappa-shear at_vertex "// &
                                       "requires the metrics slot (wet masks)")
                     error stop "vmix_apply_in_stage: at_vertex needs metrics"
                  end if
                  call kappa_shear_compute(grid, kshear, ms, dyn%therm_dt(dt), &
                                           wet_t=metrics%wet_T, wet_u=metrics%wet_u, &
                                           wet_v=metrics%wet_v)
               else
                  call kappa_shear_compute(grid, kshear, ms, dyn%therm_dt(dt))
               end if
            end if
            call kappa_shear_merge_into_kv_kt(kshear, size(vmix%kt, 1), size(vmix%kt, 2), &
                                              size(vmix%kt, 3), vmix%kv, vmix%kt)
         end if
         if (tidal_active) then
            if (dyn%enable_thermodynamics .and. dyn%is_thermo_step() .and. stage == 1) then
               call tidal_mixing_compute(grid, vmix_tidal, ms, dyn%therm_dt(dt))
            end if
            call tidal_mixing_merge_into_kt(vmix_tidal, size(vmix%kt, 1), size(vmix%kt, 2), &
                                            size(vmix%kt, 3), vmix%kv, vmix%kt)
         end if
         call vmix_add_kv_ml_invz2(grid, vmix, ms)
         if (vmix%conv_enable) then
            ! Brunt-Vaisala-triggered convective adjustment.  A CONTRIBUTOR
            ! (max() floor), so it must run BEFORE vmix_assemble -- see the
            ! rdb_ocean_vmix module docstring for the D1-D4 divergences from
            ! MOM6's MOM_CVMix_conv.  Masks against the live BL depth of
            ! whichever surface scheme is active this stage.
            if (epbl_active) then
               call vmix_apply_convection(grid, vmix, ms, epbl%mld)
            else
               call vmix_apply_convection(grid, vmix, ms, vmix%bl_depth)
            end if
         end if
         ! `vmix_split_kd_heat_salt` is NOT a contributor -- it must be
         ! the LAST statement before `vmix_assemble`, always.  Any future
         ! PR adding a kv/kt contributor (e.g. convective adjustment)
         ! inserts ABOVE this line, never below it -- ks is derived from
         ! whatever kt holds at this point, so a contributor placed after
         ! the split silently never reaches ks.
         call vmix_split_kd_heat_salt(grid, vmix, ms)
         ! Single downstream assembly gate: background floors, kv_max/kd_max
         ! ceilings, optional 1-2-1 smoothing, optional negative/NaN guard.
         ! Defaults reproduce the pre-assembly chain bit-for-bit.  Gates kv,
         ! kt, AND ks.  `geolat` rides on the optional metrics slot (the C7
         ! Henyey latitude factor's only spatial input) — forwarded when the
         ! caller threaded metrics through, omitted otherwise.  Omitting it
         ! is safe: vmix_assemble fails loud if `bkgnd_henyey` is on and
         ! geolat is absent, so a caller that forgot cannot silently lose
         ! the latitude factor.
         if (present(metrics)) then
            call vmix_assemble(grid, vmix, ms, geolat=metrics%geolatT)
         else
            call vmix_assemble(grid, vmix, ms)
         end if
         call profiler_stop("ocean_vmix_compute")
         call profiler_start("ocean_vdiff_apply")
         if (vertex_kv) then
            ! Corner Kv seam: same calls + the corner viscosity source.
            if (do_remnant) then
               call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                         tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                         lambda_bot_u=bd%lambda_bot_u, &
                                         lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                         visc_rem_u=bt_work%visc_rem_u, visc_rem_v=bt_work%visc_rem_v, &
                                         kv_corner_source=kshear%kd_corner, &
                                         kv_corner_prandtl=kshear%prandtl_turb)
            else
               call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                         tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                         lambda_bot_u=bd%lambda_bot_u, &
                                         lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                         kv_corner_source=kshear%kd_corner, &
                                         kv_corner_prandtl=kshear%prandtl_turb)
            end if
         else if (do_remnant) then
            call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                      visc_rem_u=bt_work%visc_rem_u, visc_rem_v=bt_work%visc_rem_v)
         else
            call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0)
         end if
         if (do_tracers .and. dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
            call vdiff_apply_tracers(grid, vd, ms, dyn%therm_dt(dt), &
                                     kt_source=vmix%kt, ks_source=vmix%ks)
            if (vmix%use_kpp) then
               call vmix_apply_nonlocal_tendencies(grid, vmix, ms, dyn%therm_dt(dt))
            end if
         end if
         call profiler_stop("ocean_vdiff_apply")
      else
         call profiler_start("ocean_vdiff_apply")
         if (do_remnant) then
            call vdiff_apply_momentum(grid, vd, ms, dt, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                      visc_rem_u=bt_work%visc_rem_u, visc_rem_v=bt_work%visc_rem_v)
         else
            call vdiff_apply_momentum(grid, vd, ms, dt, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0)
         end if
         if (do_tracers .and. dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
            call vdiff_apply_tracers(grid, vd, ms, dyn%therm_dt(dt))
         end if
         call profiler_stop("ocean_vdiff_apply")
      end if
   end subroutine vmix_apply_in_stage

   subroutine visc_rem_precompute(grid, dyn, vmix, vd, ss, bd, ms, dt, kshear)
      !! Refresh `bt_work%visc_rem_u/v` from the CURRENT stage state
      !! BEFORE the barotropic forcing assembly (PGF_BUG.md §9) — the
      !! MOM6-order parity (`vertvisc_coef` runs before `btstep` every
      !! stage).  The stage-end producer alone leaves visc_rem at its
      !! init value (≡ 1) for the whole first stage, so the rem-weighted
      !! `F_bt` degenerates to the plain mean exactly when the spurious
      !! grounded-layer PGF is at its ballistic worst, and the Δu
      !! corrector then deposits the spurious column-mean into wet
      !! layers.  Remnant-only vdiff call: builds the same matrix the
      !! stage-end solve will build (kv one stage stale — benign; the
      !! thickness field, which the BBL glue keys on, is current) and
      !! does NOT touch the velocities.  Vertex kappa-shear: the corner
      !! Kv source enters this matrix too (same operator as the
      !! stage-end momentum solve — a remnant built without it would
      !! weight the BT corrector with a different friction operator
      !! than the one actually applied).
      type(hgrid_t), intent(in) :: grid
      type(ocean_dyn_t), intent(inout) :: dyn
      type(ocean_vmix_t), intent(in) :: vmix
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_surface_stress_t), intent(in) :: ss
      type(ocean_bottom_drag_t), intent(in) :: bd
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      type(ocean_kappa_shear_t), intent(in), optional :: kshear
         !! Kappa-shear slot; only read when enabled + vertex mode
         !! (supplies the corner Kv source).

      logical :: vertex_kv

      vertex_kv = .false.
      if (present(kshear)) vertex_kv = kshear%enable .and. kshear%at_vertex

      if (vmix%use_closure) then
         if (vertex_kv) then
            call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                      visc_rem_u=dyn%bt_work%visc_rem_u, &
                                      visc_rem_v=dyn%bt_work%visc_rem_v, &
                                      remnant_only=.true., &
                                      kv_corner_source=kshear%kd_corner, &
                                      kv_corner_prandtl=kshear%prandtl_turb)
         else
            call vdiff_apply_momentum(grid, vd, ms, dt, kv_source=vmix%kv, &
                                      tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                      lambda_bot_u=bd%lambda_bot_u, &
                                      lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                      visc_rem_u=dyn%bt_work%visc_rem_u, &
                                      visc_rem_v=dyn%bt_work%visc_rem_v, &
                                      remnant_only=.true.)
         end if
      else
         call vdiff_apply_momentum(grid, vd, ms, dt, &
                                   tau_u=ss%tau_x, tau_v=ss%tau_y, &
                                   lambda_bot_u=bd%lambda_bot_u, &
                                   lambda_bot_v=bd%lambda_bot_v, rho0=ss%rho0, &
                                   visc_rem_u=dyn%bt_work%visc_rem_u, &
                                   visc_rem_v=dyn%bt_work%visc_rem_v, &
                                   remnant_only=.true.)
      end if
   end subroutine visc_rem_precompute

   subroutine accel_visc_rem_snapshot(n1, n2, n3, vel, snap)
      !! accel_visc_rem stage-entry snapshot: `snap = vel`, device-side.
      !! Public only for the unit-test suite.  Under `mem:separate` this
      !! `do concurrent` runs on the device-resident arrays in production
      !! (both are mapped on the ocean state); unit tests must map their
      !! own arrays explicitly (`!$acc enter data copyin` / `update self`).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: vel(n1, n2, n3)
      real(wp), intent(out) :: snap(n1, n2, n3)
      integer :: i, j, k

      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         snap(i, j, k) = vel(i, j, k)
      end do
   end subroutine accel_visc_rem_snapshot

   subroutine accel_visc_rem_reweight(n1, n2, n3, snap, rem, vel)
      !! accel_visc_rem post-apply reweight:
      !! `vel = snap + rem·(vel − snap)` — the whole explicit-tendency
      !! sum accumulated since the snapshot is attenuated by the
      !! per-layer viscous remnant (linearity ⇒ identical to weighting
      !! each tendency individually, MOM6 `u = u_init + dt·visc_rem·
      !! (CAu + PFu + diffu)`).  Faces with `rem == 1` (the unconditioned
      !! `source=1.0` init, and every face before the first vdiff fills
      !! the producer) are SKIPPED, not rewritten — `snap + 1·(vel−snap)`
      !! is not an FP identity, and the skip keeps rem≡1 bitwise inert.
      !! Friction-dominated near-massless layers (`rem → 0`) keep their
      !! entry velocity.  Public only for the unit-test suite.  Under
      !! `mem:separate` this `do concurrent` runs on the device-resident
      !! arrays in production (all mapped on the ocean state); unit tests
      !! must map their own arrays explicitly.  The masked write stays
      !! inside the `do concurrent` body (legal — no cross-iteration dep).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: snap(n1, n2, n3)
      real(wp), intent(in) :: rem(n1, n2, n3)
      real(wp), intent(inout) :: vel(n1, n2, n3)
      integer :: i, j, k

      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         if (rem(i, j, k) /= 1.0_wp) then
            vel(i, j, k) = snap(i, j, k) &
                           + rem(i, j, k)*(vel(i, j, k) - snap(i, j, k))
         end if
      end do
   end subroutine accel_visc_rem_reweight

   pure subroutine apply_maxvel_clamp(ms, maxvel)
      !! Truncate face velocities to `|u| ≤ maxvel`.  MOM6's MAXVEL
      !! analogue.  Called once at the end of each outer step (after
      !! the RK2 average so the clipped state is what gets carried
      !! into the next stage).  No-op when `maxvel <= 0`.
      !!
      !! Clip preserves the sign and the velocity direction — it's
      !! a per-component clip, not a magnitude clip.  Matches MOM6's
      !! `if (abs(u) > maxvel) u = sign(maxvel, u)` form so the WBC
      !! jet that pushes through 6 m/s gets clipped to ±6 m/s without
      !! introducing direction-reversal artifacts.
      !!
      !! Not conservative: clipping a face velocity from 10 m/s to
      !! 6 m/s loses momentum.  That's acceptable as a safety net —
      !! when the clamp is active, momentum conservation is already
      !! broken by whatever generated the runaway velocity.  In
      !! production runs the clamp should fire rarely or never.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: maxvel
      integer :: i, j, k, nx_face, ny_uface, nx_vface, ny_face, nz

      if (maxvel <= 0.0_wp) return

      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      nz = ms%nz_ml

      ! NaN-safe (pdc 8c2fd674 + 8675b0d4): the ieee_is_finite guard both
      ! (a) leaves a non-finite value untouched (NEVER laundering NaN →
      ! ±maxvel) and (b) blocks nvfortran -fast from lowering the
      ! `if(u>hi)…else if(u<lo)` pair to a NaN-blind min/max clamp on GPU.
      ! Non-finite velocities are zeroed upstream by the truncation's
      ! NaN-catch; this is the second line of defence.
      do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
         if (ieee_is_finite(ms%u_face_x_layer(i, j, k))) then
            if (ms%u_face_x_layer(i, j, k) > maxvel) then
               ms%u_face_x_layer(i, j, k) = maxvel
            else if (ms%u_face_x_layer(i, j, k) < -maxvel) then
               ms%u_face_x_layer(i, j, k) = -maxvel
            end if
         end if
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         if (ieee_is_finite(ms%v_face_y_layer(i, j, k))) then
            if (ms%v_face_y_layer(i, j, k) > maxvel) then
               ms%v_face_y_layer(i, j, k) = maxvel
            else if (ms%v_face_y_layer(i, j, k) < -maxvel) then
               ms%v_face_y_layer(i, j, k) = -maxvel
            end if
         end if
      end do
   end subroutine apply_maxvel_clamp

   subroutine apply_velocity_truncation(ms, metrics, dt, cfl_trunc, maxvel, ntrunc_step, &
                                        vanish_tol, n_nanzero, clip_cell_metric, &
                                        nan_i, nan_j, nan_k, nan_is_u, nan_dx, nan_visc_cfl, nu_h)
      !! Post-RK2 velocity housekeeping: the advective-CFL truncation
      !! (E7) followed by the absolute `maxvel` cap.  Called once at the
      !! end of each outer step (after the RK2 average, before the ALE
      !! remap) so the carried-forward / remapped velocity field is
      !! bounded.  Replaces the bare `apply_maxvel_clamp` call at both
      !! `ocean_dyn_step` / `ocean_dyn_step_split` sites.
      !!
      !! Order (well-defined, idempotent):
      !!   0. **Vanished-layer zero** (Phase 3, gated `vanish_tol > 0`):
      !!      a u-face where BOTH adjacent centre thicknesses are at or
      !!      below `vanish_tol` is set to 0 and NOT counted as a CFL
      !!      truncation (it is a diagnostic-clean zero, not a real clip).
      !!      A face with at least one massive side is untouched (R1).
      !!      Runs BEFORE the CFL clip so a zeroed face cannot re-trigger
      !!      the clip.  Idempotent with Phase-2 reset (R6).
      !!   1. **CFL clip** (gated `cfl_trunc > 0`): any face whose local
      !!      advective CFL `|u|·dt·idx` exceeds `cfl_trunc` is reset to
      !!      `sign(0.9·cfl_trunc/(dt·idx), u_old)` — i.e. magnitude
      !!      `0.9·cfl_trunc·dx/dt`, sign preserved.  The `0.9`
      !!      relaxation (`CFL_TRUNC_RELAX`) keeps the clipped face below
      !!      threshold so it does not re-trip on round-off next step.
      !!      idx is the stored reciprocal metric (`idxCu`/`idyCv`).
      !!   2. **maxvel cap** (gated `maxvel > 0`, via `apply_maxvel_clamp`):
      !!      the coarser absolute physical backstop.  Runs second, so a
      !!      face clipped to `0.9·cfl_trunc·dx/dt > maxvel` is further
      !!      bounded to `±maxvel`.
      !!
      !! The CFL clip + count run in a single `!$acc parallel loop
      !! reduction(+:n)` per face component — *not* a bare `do concurrent`
      !! + `sum()` over scratch (which silently returns 0 on GPU under
      !! stdpar; see the no-managed-memory gotcha).  `ntrunc_step` is the
      !! host-scalar reduction result (reset to 0 here each call).
      !!
      !! Not conservative — clipping a runaway face loses momentum; when
      !! the truncation fires, conservation is already broken by whatever
      !! produced the runaway.  PointAccel-style storage of the worst
      !! offender's (i,j,k) + acceleration breakdown is a deferred future
      !! extension; v1 ships the count only.
      !!
      !! `vanish_tol` absent / <= 0 ⇒ step 0 skipped ⇒ bit-identical.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: cfl_trunc
      real(wp), intent(in) :: maxvel
      integer, intent(out) :: ntrunc_step
      real(wp), intent(in), optional :: vanish_tol
         !! When present (> 0): zero both-sided-vanished faces before the
         !! CFL clip.  Absent / 0 ⇒ skip ⇒ bit-identical.

      real(wp), parameter :: CFL_TRUNC_RELAX = 0.9_wp
         !! MOM6 relaxation factor: clip to 90% of the threshold so the
         !! clipped face stays sub-threshold and does not re-trip.
      integer :: nx_uface, ny_uface, nx_vface, ny_vface, nz
      integer :: nx_centre, ny_centre
      integer, intent(out), optional :: n_nanzero
         !! Count of non-finite (NaN/Inf) face velocities zeroed by the
         !! step -1 NaN-catch this call.  A non-finite velocity is a
         !! producer bug (a 0/0 upstream); it MUST be caught here so it
         !! can never launder to ±maxvel (nvfortran lowers the
         !! `if(u>hi)…else if(u<lo)` pair to a NaN-skipping min/max clamp
         !! on GPU ⇒ NaN → -maxvel).  Loud counter — 0 on a healthy run.
      logical, intent(in), optional :: clip_cell_metric
         !! When present AND .true.: the CFL clip bounds each face on the
         !! CELL metric `max(idxT(i-1,j), idxT(i,j))` (v: `idyT`) — the SAME
         !! metric the console panic / `compute_max_cfl` uses — instead of
         !! the face metric `idxCu`/`idyCv` (pdc 6991988f).  Guarantees the
         !! panic-visible CFL is bounded even where the face metric is
         !! anomalous at a grounding face (`|u|·dt·idxCu ≤ cfl_trunc` while
         !! `|u|·dt·idxT > cfl_trunc` — the θ-edge escape; measured at
         !! 1024²/dt=800 as MaxCFL 0.689 sailing past a 0.5 ceiling).
         !! Bit-identical on grids where `idxCu == idxT`.  Absent /
         !! .false. ⇒ face metric ⇒ bit-identical.
      integer, intent(out), optional :: nan_i, nan_j, nan_k
         !! Grid location (local, including ghosts) of the FIRST non-finite
         !! face this call caught, in (k,j,i)-ascending scan order —
         !! actionable in place of the old bare count ("producer 0/0
         !! upstream — investigate"; FINDINGS.md's global-tripolar-aquaplanet
         !! debugging session had nothing better to go on).  0 when
         !! `n_nanzero` is absent/0 (nothing to report) or when none of
         !! these outputs were requested (the search is skipped entirely —
         !! see `nan_dx`/`nan_visc_cfl`).
      logical, intent(out), optional :: nan_is_u
         !! `.true.` if the first non-finite face was a u-face
         !! (`u_face_x_layer`), `.false.` if a v-face. Only meaningful
         !! when `nan_i`/`nan_j`/`nan_k` were actually found (`n_nan > 0`
         !! AND the location search ran).
      real(wp), intent(out), optional :: nan_dx
         !! Local cell size (m) at `(nan_i, nan_j)` — `1/max(idxT(i-1,j),
         !! idxT(i,j))` (v-face: idyT analogue), the SAME cell metric
         !! `clip_cell_metric` uses. Cheap (one extra metrics read at the
         !! single located cell, not a scan) — always computed alongside
         !! the location when a location was found.
      real(wp), intent(out), optional :: nan_visc_cfl
         !! Local viscous CFL `nu_h*dt/nan_dx^2` at the located cell —
         !! only computed when `nu_h` is supplied (the caller's
         !! `ocean_horizontal_viscosity_t%nu_h`); 0 otherwise. Mirrors
         !! `rdb_ocean_stability_audit`'s configure-time check, evaluated
         !! HERE at the exact runtime location the NaN was caught, so a
         !! user reading the message can immediately see whether an
         !! under-resolved viscous CFL is the likely producer.
      real(wp), intent(in), optional :: nu_h
         !! Constant horizontal viscosity (m^2/s), for `nan_visc_cfl`
         !! only. Absent ⇒ `nan_visc_cfl` (if requested) is 0.
      integer :: i, j, k, n_u, n_v, n_nan
      real(wp) :: idx_use, idy_use
      logical :: do_clip_cell
      real(wp) :: vt
      logical :: do_vanish_zero
      logical :: want_nan_loc
      integer(int64) :: lin_u_min, lin_v_min, lin
      integer(int64) :: nxu64, nyu64, nxv64, nyv64

      do_clip_cell = .false.
      if (present(clip_cell_metric)) do_clip_cell = clip_cell_metric
      ntrunc_step = 0
      do_vanish_zero = .false.
      vt = 0.0_wp
      if (present(vanish_tol)) then
         if (vanish_tol > 0.0_wp) then
            do_vanish_zero = .true.
            vt = vanish_tol
         end if
      end if

      ! Step -1: NaN/Inf CATCH (loud, unconditional; pdc 8c2fd674).  A
      ! non-finite face velocity — born of a 0/0 in an upstream producer
      ! (implicit vdiff on an all-at-floor column, a per-column BT-fold/PGF
      ! division at a grounding θ-edge) — is NaN-false to the CFL clip's
      ! `abs(u)…>cfl_trunc` (no clip) and LAUNDERS to ±maxvel in
      ! apply_maxvel_clamp (nvfortran lowers `if(u>hi)…else if(u<lo)…` to a
      ! NaN-skipping min/max clamp on GPU).  It then transports mass at
      ! maxvel.  Zero it here — BEFORE clip + maxvel — and count it loudly
      ! (max-reductions are NaN-blind, so the count is explicit).  Split
      ! count/zero (write out of the reduction loop, per the clip fix).
      ! Finite fields ⇒ 0 non-finite ⇒ bit-identical.
      nz = ms%nz_ml
      n_nan = 0
      nx_uface = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_vface = size(ms%v_face_y_layer, 2)
      !$acc parallel loop collapse(3) reduction(+:n_nan) present(ms%u_face_x_layer)
      do k = 1, nz
         do j = 1, ny_uface
            do i = 1, nx_uface
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) n_nan = n_nan + 1
            end do
         end do
      end do
      !$acc parallel loop collapse(3) reduction(+:n_nan) present(ms%v_face_y_layer)
      do k = 1, nz
         do j = 1, ny_vface
            do i = 1, nx_vface
               if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) n_nan = n_nan + 1
            end do
         end do
      end do
      ! Locate the FIRST non-finite face (min-reduction over an encoded
      ! linear index, k-major) — ONLY when n_nan>0 (something is already
      ! broken) AND the caller actually asked for a location (any of
      ! nan_i/nan_j/nan_k/nan_is_u/nan_dx/nan_visc_cfl present). This is
      ! the per-step-path cost gate: a healthy run (n_nan==0, the
      ! overwhelming common case) never runs this block at all, and a
      ! caller not asking for a location (the legacy call sites) pays
      ! nothing either. MUST run BEFORE the zero-out below — the zeroed
      ! array has nothing left to locate.
      want_nan_loc = present(nan_i) .or. present(nan_j) .or. present(nan_k) .or. &
                     present(nan_is_u) .or. present(nan_dx) .or. present(nan_visc_cfl)
      if (present(nan_i)) nan_i = 0
      if (present(nan_j)) nan_j = 0
      if (present(nan_k)) nan_k = 0
      if (present(nan_is_u)) nan_is_u = .true.
      if (present(nan_dx)) nan_dx = 0.0_wp
      if (present(nan_visc_cfl)) nan_visc_cfl = 0.0_wp
      if (n_nan > 0 .and. want_nan_loc) then
         nxu64 = int(nx_uface, int64)
         nyu64 = int(ny_uface, int64)
         nxv64 = int(nx_vface, int64)
         nyv64 = int(ny_vface, int64)
         lin_u_min = huge(1_int64)
         !$acc parallel loop collapse(3) reduction(min:lin_u_min) present(ms%u_face_x_layer)
         do k = 1, nz
            do j = 1, ny_uface
               do i = 1, nx_uface
                  if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) then
                     lin_u_min = min(lin_u_min, &
                                     (int(k - 1, int64)*nyu64 + int(j - 1, int64))*nxu64 + int(i - 1, int64))
                  end if
               end do
            end do
         end do
         lin_v_min = huge(1_int64)
         !$acc parallel loop collapse(3) reduction(min:lin_v_min) present(ms%v_face_y_layer)
         do k = 1, nz
            do j = 1, ny_vface
               do i = 1, nx_vface
                  if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) then
                     lin_v_min = min(lin_v_min, &
                                     (int(k - 1, int64)*nyv64 + int(j - 1, int64))*nxv64 + int(i - 1, int64))
                  end if
               end do
            end do
         end do

         block
            integer :: fi, fj, fk
            logical :: found_u
            real(wp) :: idx_local, dx_local
            integer(int64) :: lin_dec, lin_rem

            ! n_nan>0 guarantees at least one of lin_u_min/lin_v_min is a
            ! real (non-huge) encoded index; a tie prefers u.
            found_u = lin_u_min <= lin_v_min
            if (found_u) then
               lin_dec = lin_u_min
               fi = int(mod(lin_dec, nxu64), kind(fi)) + 1
               lin_rem = lin_dec/nxu64
               fj = int(mod(lin_rem, nyu64), kind(fj)) + 1
               fk = int(lin_rem/nyu64, kind(fk)) + 1
            else
               lin_dec = lin_v_min
               fi = int(mod(lin_dec, nxv64), kind(fi)) + 1
               lin_rem = lin_dec/nxv64
               fj = int(mod(lin_rem, nyv64), kind(fj)) + 1
               fk = int(lin_rem/nyv64, kind(fk)) + 1
            end if

            if (present(nan_i)) nan_i = fi
            if (present(nan_j)) nan_j = fj
            if (present(nan_k)) nan_k = fk
            if (present(nan_is_u)) nan_is_u = found_u

            if (present(nan_dx) .or. present(nan_visc_cfl)) then
               ! Host scalar read of metrics%idxT/idyT: SAFE under
               ! mem:separate despite metrics being device-resident —
               ! these are static, `copyin`-once metric arrays (never
               ! re-written after configure_ocean_metrics), so the host
               ! and device copies stay identical for the life of the
               ! run; no `!$acc update self` needed (contrast with a
               ! per-step state array, which would be stale here).
               nx_centre = size(ms%h_layer, 1)
               ny_centre = size(ms%h_layer, 2)
               if (found_u) then
                  idx_local = max(metrics%idxT(max(fi - 1, 1), fj), &
                                  metrics%idxT(min(fi, nx_centre), fj))
               else
                  idx_local = max(metrics%idyT(fi, max(fj - 1, 1)), &
                                  metrics%idyT(fi, min(fj, ny_centre)))
               end if
               dx_local = 0.0_wp
               if (idx_local > 0.0_wp) dx_local = 1.0_wp/idx_local
               if (present(nan_dx)) nan_dx = dx_local
               if (present(nan_visc_cfl)) then
                  if (present(nu_h) .and. dx_local > 0.0_wp) then
                     nan_visc_cfl = nu_h*dt/dx_local**2
                  end if
               end if
            end if
         end block
      end if

      if (n_nan > 0) then
         do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_uface)
            if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) ms%u_face_x_layer(i, j, k) = 0.0_wp
         end do
         do concurrent(k=1:nz, j=1:ny_vface, i=1:nx_vface)
            if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) ms%v_face_y_layer(i, j, k) = 0.0_wp
         end do
      end if
      if (present(n_nanzero)) n_nanzero = n_nan

      ! Step 0: zero both-sided-vanished faces (Phase 3, BEFORE CFL clip).
      ! A face zeroed here does NOT increment ntrunc_step (it is not a CFL
      ! clip — it is a diagnostic-clean zero of a face with no real mass).
      ! Data-parallel masked write (no reduction) ⇒ do concurrent with the
      ! both-sided-vanished test folded into the loop mask; mirrors
      ! reset_vanished_layer_velocities.  ms arrays are device-resident.
      if (do_vanish_zero) then
         nx_uface = size(ms%u_face_x_layer, 1)
         ny_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         ny_vface = size(ms%v_face_y_layer, 2)
         nx_centre = size(ms%h_layer, 1)
         ny_centre = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ! u-faces: i straddles centres i-1 and i; safe range i=2..nx_uface-1
         ! (both-sided-vanished test as an inner `if`, not a DC mask — masked
         ! headers are not auto-convertible by tools/dc_to_omp.py)
         do concurrent(k=1:nz, j=1:ny_uface, i=2:nx_uface - 1)
            if (max(ms%h_layer(i - 1, j, k), ms%h_layer(i, j, k)) <= vt) then
               ms%u_face_x_layer(i, j, k) = 0.0_wp
            end if
         end do
         ! v-faces: j straddles centres j-1 and j; safe range j=2..ny_vface-1
         do concurrent(k=1:nz, j=2:ny_vface - 1, i=1:nx_vface)
            if (max(ms%h_layer(i, j - 1, k), ms%h_layer(i, j, k)) <= vt) then
               ms%v_face_y_layer(i, j, k) = 0.0_wp
            end if
         end do
      end if

      if (cfl_trunc > 0.0_wp) then
         nx_uface = size(ms%u_face_x_layer, 1)
         ny_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         ny_vface = size(ms%v_face_y_layer, 2)
         nz = ms%nz_ml
         n_u = 0
         n_v = 0

         ! u-faces (Cu): clip in TWO passes — a read-only COUNT reduction,
         ! then a masked-write `do concurrent` CLIP.  A conditional array
         ! WRITE fused into an `!$acc parallel loop reduction` silently
         ! fails to land on the cc70 600² build (the 6 m/s face survived
         ! an unconditionally-true clip — the truncation was partially
         ! inert, "truncations climb but don't save").  Both
         ! passes read the ORIGINAL u ⇒ same faces, same ntrunc_step ⇒
         ! bit-identical to a correct fused loop.
         ! Clip metric: face idxCu by default; the CELL metric
         ! max(idxT(i-1), idxT(i)) — the panic basis — under
         ! clip_cell_metric.  On idxCu==idxT grids the two are identical
         ! ⇒ bit-identical.
         nx_centre = size(ms%h_layer, 1)
         ny_centre = size(ms%h_layer, 2)
         !$acc parallel loop collapse(3) reduction(+:n_u) private(idx_use) &
         !$acc   present(ms%u_face_x_layer, metrics%idxCu, metrics%idxT)
         do k = 1, nz
            do j = 1, ny_uface
               do i = 1, nx_uface
                  if (do_clip_cell) then
                     idx_use = max(metrics%idxT(max(i - 1, 1), j), &
                                   metrics%idxT(min(i, nx_centre), j))
                  else
                     idx_use = metrics%idxCu(i, j)
                  end if
                  if (abs(ms%u_face_x_layer(i, j, k))*dt*idx_use > cfl_trunc) then
                     n_u = n_u + 1
                  end if
               end do
            end do
         end do
         do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_uface) local(idx_use)
            if (do_clip_cell) then
               idx_use = max(metrics%idxT(max(i - 1, 1), j), &
                             metrics%idxT(min(i, nx_centre), j))
            else
               idx_use = metrics%idxCu(i, j)
            end if
            if (abs(ms%u_face_x_layer(i, j, k))*dt*idx_use > cfl_trunc) then
               ms%u_face_x_layer(i, j, k) = sign( &
                                            CFL_TRUNC_RELAX*cfl_trunc/(dt*idx_use), &
                                            ms%u_face_x_layer(i, j, k))
            end if
         end do

         ! v-faces (Cv): same two-pass split as the u-faces above.
         !$acc parallel loop collapse(3) reduction(+:n_v) private(idy_use) &
         !$acc   present(ms%v_face_y_layer, metrics%idyCv, metrics%idyT)
         do k = 1, nz
            do j = 1, ny_vface
               do i = 1, nx_vface
                  if (do_clip_cell) then
                     idy_use = max(metrics%idyT(i, max(j - 1, 1)), &
                                   metrics%idyT(i, min(j, ny_centre)))
                  else
                     idy_use = metrics%idyCv(i, j)
                  end if
                  if (abs(ms%v_face_y_layer(i, j, k))*dt*idy_use > cfl_trunc) then
                     n_v = n_v + 1
                  end if
               end do
            end do
         end do
         do concurrent(k=1:nz, j=1:ny_vface, i=1:nx_vface) local(idy_use)
            if (do_clip_cell) then
               idy_use = max(metrics%idyT(i, max(j - 1, 1)), &
                             metrics%idyT(i, min(j, ny_centre)))
            else
               idy_use = metrics%idyCv(i, j)
            end if
            if (abs(ms%v_face_y_layer(i, j, k))*dt*idy_use > cfl_trunc) then
               ms%v_face_y_layer(i, j, k) = sign( &
                                            CFL_TRUNC_RELAX*cfl_trunc/(dt*idy_use), &
                                            ms%v_face_y_layer(i, j, k))
            end if
         end do

         ntrunc_step = n_u + n_v
      end if

      ! Absolute physical backstop (MOM6 MAXVEL); no-op when maxvel <= 0.
      call apply_maxvel_clamp(ms, maxvel)
   end subroutine apply_velocity_truncation

   pure subroutine save_state(ms)
      !! Copy h_layer, u_face_x_layer, v_face_y_layer into the *_0
      !! save buffers on the multilayer state.  Per-tracer hTr is
      !! saved by the caller via copy_field_3d (the registry has to
      !! be walked on the host).
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j, k, nx, ny, nz, nx_face, ny_uface, nx_vface, ny_face

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%h_layer0(i, j, k) = ms%h_layer(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
         ms%u_face_x_layer0(i, j, k) = ms%u_face_x_layer(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer0(i, j, k) = ms%v_face_y_layer(i, j, k)
      end do
   end subroutine save_state

   pure subroutine restore_state(ms)
      !! Copy the *_0 save buffers back into h_layer / u_face_x_layer /
      !! v_face_y_layer — the pred_corr between-stage reset (SPEC §2): the
      !! predictor's provisional up/vp/hp are DISCARDED (only u_av/v_av/
      !! h_av survive it), and the corrector advances from u^n / h^n.
      !! Tracers are untouched by the predictor (TR_MODE_NONE + no
      !! thermodynamics), so no tracer restore is needed.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j, k, nx, ny, nz, nx_face, ny_uface, nx_vface, ny_face

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%h_layer(i, j, k) = ms%h_layer0(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
         ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer0(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer0(i, j, k)
      end do
   end subroutine restore_state

   subroutine ocean_accumulate_mass_out(ms, flux_h_layer, areaT, nghost, dt, weight)
      !! Accumulate the net mass (kg) that left the domain this RK stage into
      !! `ms%mass_out`.  `flux_h_layer` is the total horizontal divergence
      !! (`h_layer -= dt·flux_h_layer`), so `-Σ_interior(−flux_h_layer)·areaT`
      !! is the boundary outflux (interior faces telescope — divergence
      !! theorem), and it is the SAME field the thickness update consumes, so
      !! with the RK2 stage weight this closes the mass budget to round-off.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: flux_h_layer(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      real(wp), intent(in) :: dt, weight
      integer :: i, j, k, nz, nx, ny, i_lo, i_hi, j_lo, j_hi
      real(wp) :: acc

      nx = min(size(flux_h_layer, 1), size(areaT, 1))
      ny = min(size(flux_h_layer, 2), size(areaT, 2))
      nz = size(flux_h_layer, 3)
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      acc = 0.0_wp
      !$acc parallel loop collapse(3) reduction(+:acc) present(flux_h_layer, areaT)
      do k = 1, nz
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               acc = acc + flux_h_layer(i, j, k)*areaT(i, j)
            end do
         end do
      end do
      ! h_layer -= dt·flux_h_layer ⇒ mass change = −dt·ρ·Σ; outflux (positive
      ! = leaving) is its negative.
      ms%mass_out = ms%mass_out + weight*dt*RHO_WATER*acc
      ms%mass_out_tracked = .true.
   end subroutine ocean_accumulate_mass_out

   pure subroutine rk2_average(ms)
      !! State <- 0.5 * (state_0 + state) for h_layer, u_face_x_layer,
      !! v_face_y_layer.  Tracer averages are done by the caller.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j, k, nx, ny, nz, nx_face, ny_uface, nx_vface, ny_face

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%h_layer(i, j, k) = 0.5_wp*(ms%h_layer0(i, j, k) + ms%h_layer(i, j, k))
      end do
      do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
         ms%u_face_x_layer(i, j, k) = 0.5_wp*(ms%u_face_x_layer0(i, j, k) + &
                                              ms%u_face_x_layer(i, j, k))
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = 0.5_wp*(ms%v_face_y_layer0(i, j, k) + &
                                              ms%v_face_y_layer(i, j, k))
      end do
   end subroutine rk2_average

   pure subroutine mask_layer_velocities(grid, metrics, ms, bt_work)
      !! Zero the per-layer face velocities at land faces (spec §14 C4 /
      !! R4(b.2) / MOM6 `up = mask2dCu·(u+dt·accel)`).  Runs once per
      !! RK2 stage AFTER
      !! `apply_bt_correction` and inside the RK2 average so the
      !! additive layer tendencies + the BT correction cannot leave a
      !! re-ingested land-face velocity (a slow conservation leak).  The
      !! transports themselves already ride zeroed metrics; this resets
      !! the prognostic velocity so `derive_bt_from_layers` next step
      !! sees zero there.  All-wet ⇒ `wet_u/v≡1` ⇒ literal no-op.
      !!
      !! When `bt_work` is passed AND wet/dry is enabled
      !! (docs/ocean_wetdry_plan.md §4), the DYNAMIC bed-blocking face
      !! mask (`wd_open_u/v`, filled by the last BT substep's Pass 2)
      !! composes multiplicatively with the static one — layer momentum
      !! at a dynamically blocked (drying-front) face resets to zero, so
      !! rewetting starts from rest.  Knob off / absent ⇒ the ORIGINAL
      !! loops run untouched (byte-identical).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(in), optional :: bt_work
      integer :: i, j, k, nz, nx_face, ny_uface, nx_vface, ny_face
      logical :: wd_on

      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      wd_on = .false.
      if (present(bt_work)) wd_on = bt_work%wetdry_enable

      if (wd_on) then
         do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
            ms%u_face_x_layer(i, j, k) = metrics%wet_u(i, j)* &
                                         bt_work%wd_open_u(i, j)* &
                                         ms%u_face_x_layer(i, j, k)
         end do
         do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
            ms%v_face_y_layer(i, j, k) = metrics%wet_v(i, j)* &
                                         bt_work%wd_open_v(i, j)* &
                                         ms%v_face_y_layer(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
            ms%u_face_x_layer(i, j, k) = metrics%wet_u(i, j)*ms%u_face_x_layer(i, j, k)
         end do
         do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
            ms%v_face_y_layer(i, j, k) = metrics%wet_v(i, j)*ms%v_face_y_layer(i, j, k)
         end do
      end if
   end subroutine mask_layer_velocities

   pure subroutine reset_vanished_layer_velocities(ms, vanish_tol)
      !! Zero the per-layer face velocity at any face where BOTH adjacent
      !! centre-cell thicknesses are at or below `vanish_tol`
      !! (`isopycnal_vanish_tol(angstrom_h)` = `max(angstrom_h, H_VANISHED)`).
      !!
      !! **Both-sided rule (R1/R6):** a face is reset only when BOTH
      !! neighbouring centre cells are vanished.  A face with one massive
      !! side (`h >> vanish_tol`) is LEFT UNTOUCHED — that flux is
      !! legitimate (it may be re-wetting the thin layer at a grounding
      !! FRONT) and zeroing it would create a spurious vorticity dipole.
      !!
      !! **Idempotency (R6):** faces zeroed here are already 0 when
      !! Phase-3 CFL truncation runs (`apply_velocity_truncation`);
      !! the composition is harmless.
      !!
      !! Mirroring `mask_layer_velocities`: plain `do concurrent`, no
      !! `!$acc` annotations, pure, explicit dims before arrays.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: vanish_tol
      integer :: i, j, k, nz, nx_face, ny_uface, nx_vface, ny_face
      integer :: nx_centre, ny_centre

      nz = ms%nz_ml
      nx_centre = size(ms%h_layer, 1)
      ny_centre = size(ms%h_layer, 2)
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)

      ! u-face at (i,j,k) straddles centres (i-1,j,k) and (i,j,k).
      ! Skip outermost columns (i=1 and i=nx_face) which are ghost/wall
      ! faces — their adjacent interior index would be out of range.
      ! Safe inner range: i=2..nx_face-1 (both i-1 and i are in 1..nx_centre).
      ! Both-sided-vanished test as an inner `if`, not a DC mask — masked
      ! headers are not auto-convertible by tools/dc_to_omp.py.
      do concurrent(k=1:nz, j=1:ny_uface, i=2:nx_face - 1)
         if (max(ms%h_layer(i - 1, j, k), ms%h_layer(i, j, k)) <= vanish_tol) then
            ms%u_face_x_layer(i, j, k) = 0.0_wp
         end if
      end do

      ! v-face at (i,j,k) straddles centres (i,j-1,k) and (i,j,k).
      ! Safe inner range: j=2..ny_face-1.
      do concurrent(k=1:nz, j=2:ny_face - 1, i=1:nx_vface)
         if (max(ms%h_layer(i, j - 1, k), ms%h_layer(i, j, k)) <= vanish_tol) then
            ms%v_face_y_layer(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine reset_vanished_layer_velocities

   pure subroutine copy_field_3d(src, dst, nx, ny, nz)
      !! src -> dst flat copy on the device.  Bare-array shim avoids
      !! the deep struct deref inside do-concurrent (tracer hTr lives
      !! in an array-of-derived-types registry).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)  :: src(nx, ny, nz)
      real(wp), intent(out) :: dst(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         dst(i, j, k) = src(i, j, k)
      end do
   end subroutine copy_field_3d

   pure subroutine rk2_average_field_3d(saved, current, nx, ny, nz)
      !! current <- 0.5 * (saved + current) on the device.  Same
      !! bare-array shim rationale as copy_field_3d.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)    :: saved(nx, ny, nz)
      real(wp), intent(inout) :: current(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         current(i, j, k) = 0.5_wp*(saved(i, j, k) + current(i, j, k))
      end do
   end subroutine rk2_average_field_3d

   subroutine ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, dt, n_inner, sf, geo, vcoord, bc, sp, t, &
                                   lateral_mix, epbl, kshear, mle, slopes, gm, varmix, wavespeed, &
                                   redi, meke, vmix_tidal, tides, psurf)
      !! Split-explicit SSP-RK2 outer step on the multilayer state.
      !! Parallel to `ocean_dyn_step` (the unsplit driver
      !! still ships for tests + reference).  Phase 4b-MVP scope:
      !! gravity-wave-stability fix on momentum only.  ALE-aware
      !! `h_layer` redistribution and full nonlinear-bt corrections
      !! are deferred to a follow-up branch.
      !!
      !! Algorithm per SSP-RK2 stage:
      !!   1. Derive (η^n, u_bt^n, v_bt^n) from current `ms` and stash
      !!      u_bt^n, v_bt^n in `dyn%bt_work%ubt_at_n` / `dyn%bt_work%vbt_at_n`.
      !!   2. Run all slow computes (writes the per-kernel scratch
      !!      buffers).
      !!   3. Sum the per-face slow tendency from those scratches
      !!      into `dyn%bt_work%F_slow_u` / `dyn%bt_work%F_slow_v`.
      !!   4. Depth-mean → `dyn%bt_work%F_bt_u` / `dyn%bt_work%F_bt_v`.  (Weighted by
      !!      face thickness.)
      !!   5. Existing applies for momentum + continuity (uses the
      !!      same per-kernel scratches — produces u^*, v^*, h^* the
      !!      unsplit driver would).
      !!   6. Run the barotropic substep with F_bt as constant forcing for
      !!      `n_inner` substeps at dt_inner = dt / n_inner — fills
      !!      `dyn%bt_work%bt_eta` / `dyn%bt_work%bt_ubt` / `dyn%bt_work%bt_vbt` with the
      !!      time-mean.
      !!   7. Correction step: add (⟨u_bt⟩ - u_bt^n - dt·F_bt_u) to
      !!      every layer's u, v.  This replaces the unsplit bt mode
      !!      (which is FE-amplified on the gravity wave) with the
      !!      barotropic-substep's resolved bt mode (stable under FBE for
      !!      dt_inner·c·k < 2).
      !!
      !! Two RK2 stages, then RK2 average.  Per-tracer save/average
      !! is identical to the unsplit driver.
      !!
      !! `dyn%bt_work%bt_H_ref` must be set BEFORE the first split step (the
      !! caller initialises it from the time-mean total H of the
      !! initial multilayer state; for Eulerian-z with a flat
      !! bathymetry it's just the total column depth).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics — forwarded to the PGF (and,
         !! in slice 2, the other geometry-aware kernels).
      type(ocean_dyn_t), intent(inout) :: dyn
      type(eos_t), intent(in) :: eos
      type(coriolis_adv_t), intent(inout) :: cor
      type(continuity_t), intent(inout) :: ct
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      integer, intent(in) :: n_inner
      type(ocean_surface_flux_t), intent(in), optional :: sf
      type(ocean_geothermal_t), intent(in), optional :: geo
         !! Geothermal bottom-heat-flux slot.  Absent or `enable=.false.`
         !! preserves the historical no-geothermal path bit-identically.
      type(ocean_vcoord_t), intent(inout), optional :: vcoord
         !! Vertical-coordinate state.  When present and
         !! `vcoord%coord_type /= VCOORD_EULERIAN_Z`, the orchestrator
         !! runs the ALE remap step after the RK2 average (Lagrangian-
         !! then-remap pattern, MOM6-style).  Absent or EULERIAN_Z
         !! preserves the pre-existing Eulerian-z behaviour bit-
         !! identically.
      type(ocean_bc_state_t), intent(inout), optional :: bc
         !! Open-boundary config.  Absent or all-OBC_WALL preserves
         !! the closed-wall behaviour bit-identically (the barotropic substep's
         !! tag dispatch falls through to hard-zero).  When any edge
         !! is non-WALL, the corresponding BC variant fires inside
         !! the barotropic substep's per-substep wall closure.
      type(ocean_sponge_t), intent(in), optional :: sp
         !! Map-driven sponge slot (PR-23).  Absent or `enable=.false.`
         !! (default) preserves the legacy `bc`-band sponge path bit-
         !! identically; `enable=.true.` supersedes it (exactly one of the
         !! two runs — see `run_stage_split`'s dispatch).
      real(wp), intent(in), optional :: t
         !! Wall-clock time at the start of this outer step (s),
         !! used to evaluate the OBC_TIDAL constituent table.
      type(ocean_lateral_mix_t), intent(inout), optional :: lateral_mix
         !! Flow-aware lateral-viscosity closure (Leith / Smagorinsky).
         !! Absent or `closure = LMIX_NONE` falls through to the
         !! scalar `nu_h` in `hv`, bit-identically.
      type(ocean_epbl_t), intent(inout), optional :: epbl
         !! Energetics-based PBL slot.  Absent or `enable=.false.`
         !! preserves the historical PP81/KPP path bit-identically.
      type(ocean_kappa_shear_t), intent(inout), optional :: kshear
         !! Kappa-shear interior closure slot.  Absent or
         !! `enable=.false.` preserves the historical path bit-identically.
      type(ocean_tidal_mixing_t), intent(inout), optional :: vmix_tidal
         !! Tidal-mixing interior closure slot.  Absent or
         !! `enable=.false.` preserves the historical path bit-identically.
      type(ocean_mle_t), intent(inout), optional :: mle
         !! Fox-Kemper MLE slot (B5).  Absent or `enable=.false.`
         !! preserves bit-identity.  Transports are computed once per
         !! outer step at THERMO cadence (from the prior step's
         !! `epbl%mld`) and folded into the continuity mass fluxes in
         !! both RK2 stages.
      type(ocean_slopes_t), intent(inout), optional :: slopes
         !! Isopycnal-slope diagnostics slot.  Refreshed at THERMO cadence
         !! (the split driver otherwise never calls `ocean_slopes_compute`)
         !! so the GM slot has a fresh slope to consume.  Absent or
         !! disabled ⇒ no-op.
      type(ocean_gm_t), intent(inout), optional :: gm
         !! Gent-McWilliams thickness-diffusion slot (capability [2]).
         !! Absent or `enable=.false.` preserves bit-identity.  Transports
         !! are computed once per outer step at THERMO cadence (after the
         !! slope refresh) and folded into the continuity mass fluxes in
         !! both RK2 stages, exactly like `mle`.
      type(ocean_varmix_t), intent(inout), optional :: varmix
         !! VarMix slot (capability [4]): spatially-varying GM/Redi
         !! coefficients.  When present + `enable=.true.` (and wavespeed
         !! present) `varmix_compute` fills the pre-CFL base `khth_u/v` at
         !! THERMO cadence BEFORE GM, and GM consumes them as its external
         !! base.  Absent or disabled ⇒ GM uses its scalar `khth` ⇒
         !! bit-identity.
      type(ocean_wave_speed_t), intent(inout), optional :: wavespeed
         !! Wave-speed slot supplying `cg1` to the VarMix resolution
         !! function.  Only read when `varmix` is active.
      type(ocean_redi_t), intent(inout), optional :: redi
         !! Redi continuous neutral-diffusion slot (capability [3]).  Absent
         !! or `enable=.false.` preserves bit-identity.  Phase-A neutral-
         !! surface coefficients are computed ONCE per outer step at THERMO
         !! cadence here (tracer-independent geometry); Phase B applies the
         !! rotated flux per tracer inside `run_stage_split` after the
         !! along-coordinate `tracer_hdiff`.
      type(ocean_meke_t), intent(inout), optional :: meke
         !! MEKE prognostic eddy-energy slot (capability [5]).  Stepped once
         !! per outer step at THERMO cadence BETWEEN `varmix_compute` and
         !! `gm_compute_transports`: it reads `gm%gm_src` from the PREVIOUS
         !! thermo step (one-step lag) and feeds the geom-mean of its derived
         !! `kh` into `varmix%khth_u/v` (+ khtr) BEFORE GM's CFL clamp.
         !! Absent or `enable=.false.` ⇒ no-op (bit-identical).
      type(ocean_tides_t), intent(inout), optional :: tides
         !! Equilibrium body-force tide slot (C1) + scalar SAL (C2).  When
         !! present and `enable`, `eta_eq` is refreshed ONCE per outer step
         !! (held static across the inner substep loop); scalar SAL then
         !! folds `beta_sal*eta` (lagged barotropic SSH) into the combined
         !! `eta_forcing = eta_eq + eta_sal` forwarded to the barotropic
         !! PGF.  Absent / disabled / `use_sal=.false.` ⇒ bit-identical.
      type(ocean_p_surf_t), intent(inout), optional :: psurf
         !! Atmospheric surface-pressure loading / inverse barometer slot
         !! (PR-17).  When present and `enable`, the assembled surface
         !! pressure `sf%p_surf` (read only; Pa) is converted to
         !! `eta_ib = -p_surf/(rho0 g_bt)` ONCE per outer step and folded
         !! into `eta_seam = eta_ib [+ tide eta_forcing]`, which then feeds
         !! the barotropic PGF in place of the tide-only seam.  Requires the
         !! PR-12 component set (`sf%p_surf`) — guarded at configure.
         !! Absent / disabled ⇒ bit-identical.

      integer :: it, stage
      integer :: i, j, nx_ptop, ny_ptop
         !! Loop indices + extents for the E3 `ms%p_top` refresh below.
      logical :: tide_on, psurf_on, p_top_live
      real(wp) :: t_now
         !! Model time (s) for `ocean_ideal_age_young_val`; `t` fallback (PR-7).
      integer :: nan_i, nan_j, nan_k
      logical :: nan_is_u
      real(wp) :: nan_dx, nan_visc_cfl
         !! Location + local diagnostics for the FIRST non-finite face the
         !! velocity-truncation NaN-catch finds this outer step (see the
         !! "[nan-catch] outer step" report below) — always requested from
         !! `apply_velocity_truncation` (cheap: only actually searched
         !! inside that call when `n_nan>0`).

      ! ---- Ghost-band poison (debug knob) ----
      ! Sentinel-NaN every exchange-covered ghost band BEFORE any exchange
      ! or kernel.  Any kernel that reads a ghost that was not overwritten by
      ! the subsequent exchange will encounter a NaN, failing the run loudly at
      ! the offending step.  Default off (dyn%poison_ghosts = .false.) = the
      ! branch is never taken and the run is bit-identical to the unpoisoned path.
      ! Requires bc (edge topology); bc absent means no topology to gate on.
      if (dyn%poison_ghosts .and. present(bc)) then
         call ocean_poison_ghost_bands( &
            grid, ms, dyn%bt_work, ss, &
            (.not. bc%has_west) .or. bc%periodic_x, &
            (.not. bc%has_east) .or. bc%periodic_x, &
            (.not. bc%has_south) .or. bc%periodic_y, &
            (.not. bc%has_north) .or. bc%periodic_y)
      end if

      tide_on = .false.
      if (present(tides)) then
         if (tides%enable) tide_on = .true.
      end if
      ! Refresh the equilibrium-tide elevation ONCE per outer step, before
      ! the RK2 stages (tidal periods are O(hr) >> inner dt).
      if (tide_on) then
         if (present(t)) then
            call tides_update_eta_eq(tides, t)
         else
            call tides_update_eta_eq(tides, 0.0_wp)
         end if
         ! Scalar SAL (C2): fold `beta_sal*eta` into the combined seam field.
         ! Uses the LAGGED barotropic surface elevation `bt_work%bt_eta`
         ! (the SSH the barotropic substep re-derives each stage; here it
         ! still holds the previous outer step's value — standard
         ! one-step-lagged scalar SAL).  `use_sal=.false.` ⇒
         ! `eta_forcing == eta_eq`, bit-identical to C1.
         call tides_update_eta_sal(tides, dyn%bt_work%bt_eta)
      end if

      ! Atmospheric surface-pressure loading / inverse barometer (PR-17).
      ! Refresh the combined seam field ONCE per outer step (held static
      ! across the inner substep loop, exactly like the tide): fold
      ! `eta_ib = -p_surf/(rho0 g_bt)` into `eta_seam = eta_ib [+ tide
      ! eta_forcing]`.  `g_bt` is the barotropic substep gravity (what the
      ! seam consumer runs on).  `sf%p_surf` is the assembled total (read
      ! only here; seeded from p_surf_atm at configure, and — once PR-18
      ! lands — overwritten with the ice load via the ice path's own inout
      ! access before this step).  The seam is differenced in the ghost
      ! band, so a halo exchange follows the fill (single-rank no-op;
      ! correct multi-rank).  `sf%p_surf` is mapped under the PR-12
      ! component set; `validate_config` requires it (and `present(sf)`)
      ! when psurf is enabled, so the read is always device-present.
      psurf_on = .false.
      if (present(psurf) .and. present(sf)) then
         if (psurf%enable) psurf_on = .true.
      end if
      if (psurf_on) then
         if (tide_on) then
            call p_surf_update_seam(psurf, sf%p_surf, dyn%bt_work%g_bt, &
                                    eta_tide=tides%eta_forcing)
         else
            call p_surf_update_seam(psurf, sf%p_surf, dyn%bt_work%g_bt)
         end if
         call ocean_halo_centre(psurf%eta_seam, device_resident=.true.)
      end if

      ! E3: the SAME assembled `sf%p_surf` also becomes the top-of-column
      ! pressure the EOS's IN-SITU builders measure down from (`ms%p_top`,
      ! Pa) when `&ocean_psurf_nml in_eos` is set — the ice-shelf-cavity
      ! seam, where 1e6-2e7 Pa of overburden makes the historical
      ! "in-situ pressure starts at 0 at the free surface" a systematic
      ! ~4-5 kg/m^3 density error under a nonlinear EOS.  It does NOT
      ! touch the POTENTIAL density `ms%rho_layer`, which stays at the
      ! horizontally uniform `eos%p_ref` by design.  Refreshed here, once
      ! per outer step, before the PGF of the step and held static across
      ! the stages.
      !
      ! P5.0: the SAME refresh also covers `&ocean_pgf_nml p_top_in_bc`,
      ! the second consumer of `ms%p_top` — it puts the load in the
      ! FV_MOM6 pressure-stack surface BC (`pa(nz+1)`), independently of
      ! whether it also reaches the EOS arguments.  The gate is the
      ! DISJUNCTION so neither consumer can ever read a p_top that the
      ! configure-time seed left behind while `sf%p_surf` moved on.
      !
      ! Written INLINE as a `do concurrent` rather than as a call: a
      ! host-gated call handing a state array to an EXTERNAL subroutine
      ! makes nvfortran treat that array as escaping and pessimises EVERY
      ! `do concurrent` in this routine, whether or not the branch is
      ! taken (CLAUDE.md, measured at +4.8 % on `ocean_continuity`).
      ! The copy spans the WHOLE array, ghosts included, so `p_top`
      ! inherits exactly the halo validity `p_surf` has and needs no
      ! exchange of its own.
      p_top_live = .false.
      if (psurf_on) then
         if (psurf%in_eos .or. pgf%p_top_in_bc) p_top_live = .true.
      end if
      if (p_top_live) then
         nx_ptop = size(ms%p_top, 1)
         ny_ptop = size(ms%p_top, 2)
         do concurrent(j=1:ny_ptop, i=1:nx_ptop)
            ms%p_top(i, j) = sf%p_surf(i, j)
         end do
      end if

      call probe_dS(grid, ms, "outer step entry", 0, dyn%outer_step_count + 1)

      ! Fox-Kemper MLE restratification (B5): compute the ML-confined
      ! overturning transports ONCE per outer step at THERMO cadence (it
      ! is a slow buoyancy-driven process; per-stage recompute doubles
      ! cost for no benefit).  Uses the MLD diagnosed by EPBL at the end
      ! of the previous step.  `run_stage_split` then folds the same
      ! uhml/vhml into the mass fluxes in both stages.  No-op when
      ! absent / disabled.
      if (present(mle) .and. present(epbl)) then
         if (dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
            call profiler_start("ocean_foxkemper")
            ! Pass the thermo window so the per-layer FK availability cap
            ! keeps the windowed tracer-drain hprev reconstruction positive
            ! in thin (EPBL-MLD) surface layers (dt_tracer_advect_ratio > 1).
            ! At dt_therm_ratio = 1 this is dt — a tighter-but-harmless cap.
            if (present(bc)) then
               call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss, &
                                           dt_limit=dyn%therm_dt(dt), bc=bc)
            else
               call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss, &
                                           dt_limit=dyn%therm_dt(dt))
            end if
            call profiler_stop("ocean_foxkemper")
         end if
      end if

      ! Wave speed (B1): refresh the first-baroclinic gravity-wave speed
      ! `cg1` + the Rossby deformation radius `rd`/`rd_over_dx` ONCE per
      ! outer step at THERMO cadence (further gated by `n_wavespeed`,
      ! `mod(outer_step_count, n_wavespeed) == 0`; effective cadence is
      ! lcm(dt_therm_ratio, n_wavespeed) when both are > 1).  Placed at the
      ! TOP of the outer step — before both `varmix_compute` call sites
      ! and `run_meke_step` — so `cg1` is fresh for all three consumers
      ! this step; `rho_layer` is whatever the PREVIOUS outer step (or,
      ! at `outer_step_count == 0`, the initial condition) left, the same
      ! one-step lag GM/MEKE already accept and the same placement MOM6's
      ! `calc_resoln_function` uses relative to `step_MOM`.  At step 0
      ! this means `cg1 = 0` (if `rho_layer` is uniform-rho0 IC) so
      ! `Res_fn = 1` for one step — a documented decision, not a bug.
      ! No-op when absent / disabled ⇒ bit-identical.  `wavespeed_compute`
      ! is `pure` (called from `do concurrent`-adjacent code), so the
      ! profiler wrapper lives here, not inside it.
      if (present(wavespeed)) then
         if (wavespeed%is_init .and. wavespeed%enable) then
            if (dyn%enable_thermodynamics .and. dyn%is_thermo_step() .and. &
                mod(dyn%outer_step_count, wavespeed%n_wavespeed) == 0) then
               call profiler_start("ocean_wavespeed")
               call wavespeed_compute(grid, metrics, wavespeed, ms)
               call profiler_stop("ocean_wavespeed")
            end if
         end if
      end if

      ! Resolution-scaled momentum viscosity (Gap 1): the lateral-mix
      ! kernels read the VarMix resolution function `res_fn_u/v`.  When GM
      ! is active it owns the VarMix refresh (below), so this standalone
      ! pass runs ONLY when GM is absent / disabled but `resoln_scaled_visc`
      ! still needs a fresh `res_fn` — fired once per outer step at THERMO
      ! cadence (`res_fn` is slowly varying; it stays device-resident
      ! between thermo steps).  No-op when VarMix is disabled or the knob
      ! is off ⇒ bit-identical.  `res_fn` depends only on the static grid
      ! terms + cg1 (not on the slopes GM refreshes), so a no-GM run still
      ! yields the correct resolution function.
      if (present(varmix) .and. present(wavespeed) .and. &
          present(lateral_mix) .and. present(slopes)) then
         if (varmix%enable .and. lateral_mix%resoln_scaled_visc) then
            if (dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
               if (.not. gm_refreshes_varmix(gm, slopes, dyn)) then
                  ! `slopes` is guaranteed present here; `varmix_compute`
                  ! self-gates on `slopes%is_init` (the configure invariant
                  ! ensures VarMix-enabled ⇒ slopes allocated) and uses the
                  ! slope only for the (here-zero) Eady SN term — the
                  ! resolution function it fills needs only cg1.
                  call varmix_compute(grid, metrics, varmix, slopes, &
                                      wavespeed, ms)
               end if
            end if
         end if
      end if

      ! Gent-McWilliams thickness diffusion (capability [2]): refresh the
      ! isopycnal slopes then compute the bolus thickness transports ONCE
      ! per outer step at THERMO cadence (a slow buoyancy-driven process).
      ! `run_stage_split` then folds the same uhD/vhD into the mass fluxes
      ! in both stages.  The split driver does not otherwise call
      ! `ocean_slopes_compute`, so GM owns the slope refresh.  No-op when
      ! absent / disabled.
      if (present(gm)) then
         if (dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
            if (gm%enable .and. present(slopes)) then
               call profiler_start("ocean_gm")
               call ocean_slopes_compute(grid, metrics, eos, slopes, ms, dyn%therm_dt(dt))
               ! VarMix (capability [4]): fill the spatially-varying pre-CFL
               ! base KhTh face field from the fresh slopes + cg1, then hand
               ! it to GM as the external base.  When VarMix is absent /
               ! disabled GM falls back to its scalar khth (bit-identical).
               if (present(varmix) .and. present(wavespeed)) then
                  if (varmix%enable) then
                     call varmix_compute(grid, metrics, varmix, slopes, &
                                         wavespeed, ms)
                     ! MEKE (capability [5]): step the prognostic eddy-energy
                     ! field AFTER VarMix and BEFORE GM so it reads the prior
                     ! step's gm_src (one-step lag) and adds the geom-mean kh
                     ! into varmix%khth before GM's CFL clamp.  No-op when
                     ! absent / disabled.
                     call run_meke_step(grid, metrics, gm, varmix, wavespeed, hv, &
                                        ms, dyn%therm_dt(dt), meke)
                     call gm_compute_transports(grid, metrics, gm, slopes, ms, &
                                                dyn%therm_dt(dt), &
                                                khth_ext_u=varmix%khth_u, &
                                                khth_ext_v=varmix%khth_v)
                  else
                     ! VarMix off ⇒ MEKE still evolves E (feedback inert, no
                     ! face accumulator to feed) — pass varmix so meke_step
                     ! can still no-op cleanly on the seam.
                     call run_meke_step(grid, metrics, gm, varmix, wavespeed, hv, &
                                        ms, dyn%therm_dt(dt), meke)
                     call gm_compute_transports(grid, metrics, gm, slopes, ms, &
                                                dyn%therm_dt(dt))
                  end if
               else
                  call run_meke_step(grid, metrics, gm, varmix, wavespeed, hv, &
                                     ms, dyn%therm_dt(dt), meke)
                  call gm_compute_transports(grid, metrics, gm, slopes, ms, &
                                             dyn%therm_dt(dt))
               end if
               call profiler_stop("ocean_gm")
            end if
         end if
      end if

      ! Redi neutral diffusion (capability [3]): build the tracer-
      ! INDEPENDENT neutral-surface coefficients ONCE per outer step at
      ! THERMO cadence (a slow geometry — recomputing per stage/tracer is
      ! wasteful).  `run_stage_split` then applies the rotated flux to each
      ! tracer (Phase B) after `tracer_hdiff`.  No-op when absent / disabled.
      if (present(redi)) then
         if (dyn%enable_thermodynamics .and. dyn%is_thermo_step()) then
            if (redi%enable) then
               call profiler_start("ocean_redi")
               call redi_calc_coeffs(grid, metrics, eos, redi, ms)
               call profiler_stop("ocean_redi")
            end if
         end if
      end if

      ! SPEC S1 seed (pred_corr): before the very first step the time-mean
      ! family must hold the initial state — the first predictor's
      ! CorAd/hvisc read u_av/h_av, and the allocation default (0) would
      ! hand them a zero-thickness field.
      if (dyn%split_scheme == SPLIT_SCHEME_PRED_CORR .and. dyn%outer_step_count == 0) then
         call copy_field_3d(ms%u_face_x_layer, ms%u_av_layer, &
                            size(ms%u_face_x_layer, 1), size(ms%u_face_x_layer, 2), &
                            size(ms%u_face_x_layer, 3))
         call copy_field_3d(ms%v_face_y_layer, ms%v_av_layer, &
                            size(ms%v_face_y_layer, 1), size(ms%v_face_y_layer, 2), &
                            size(ms%v_face_y_layer, 3))
         call copy_field_3d(ms%h_layer, ms%h_av_layer, &
                            size(ms%h_layer, 1), size(ms%h_layer, 2), &
                            size(ms%h_layer, 3))
      end if

      call save_state(ms)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            call copy_field_3d(ms%tracers(it)%hTr, ms%tracers(it)%hTr0, &
                               size(ms%tracers(it)%hTr, 1), &
                               size(ms%tracers(it)%hTr, 2), &
                               size(ms%tracers(it)%hTr, 3))
         end do
      end if

      ! Forward all optionals straight through to `run_stage_split` —
      ! Fortran 2008+ propagates `present(...)` correctly across
      ! optional dummy arguments, so an absent optional passed as the
      ! actual arg stays absent in the callee.  This collapses the
      ! previous 4-branch cartesian-product dispatch (vcoord × bc)
      ! into a single stage loop.
      do stage = 1, 2
         if (psurf_on) then
            ! `eta_seam` already includes the tide when it is on (folded in
            ! p_surf_update_seam above), so this single branch subsumes the
            ! tide-on case; the two branches below are the pre-PR-17 code.
            call run_stage_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                 va, hd, vd, vmix, ms, dt, n_inner, &
                                 sf=sf, geo=geo, stage=stage, vcoord=vcoord, bc=bc, sp=sp, t=t, &
                                 lateral_mix=lateral_mix, epbl=epbl, kshear=kshear, mle=mle, gm=gm, &
                                 redi=redi, varmix=varmix, vmix_tidal=vmix_tidal, meke=meke, &
                                 eta_forcing=psurf%eta_seam)
         else if (tide_on) then
            call run_stage_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                 va, hd, vd, vmix, ms, dt, n_inner, &
                                 sf=sf, geo=geo, stage=stage, vcoord=vcoord, bc=bc, sp=sp, t=t, &
                                 lateral_mix=lateral_mix, epbl=epbl, kshear=kshear, mle=mle, gm=gm, &
                                 redi=redi, varmix=varmix, vmix_tidal=vmix_tidal, meke=meke, &
                                 eta_forcing=tides%eta_forcing)
         else
            call run_stage_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                 va, hd, vd, vmix, ms, dt, n_inner, &
                                 sf=sf, geo=geo, stage=stage, vcoord=vcoord, bc=bc, sp=sp, t=t, &
                                 lateral_mix=lateral_mix, epbl=epbl, kshear=kshear, mle=mle, gm=gm, &
                                 redi=redi, varmix=varmix, vmix_tidal=vmix_tidal, meke=meke)
         end if
         ! pred_corr between-stage reset (SPEC §2): the predictor's
         ! provisional up/vp/hp are discarded — only u_av/v_av/h_av carry
         ! forward — and the corrector advances from u^n / h^n.  The
         ! stage-2 entry then re-derives the barotropic state from the
         ! restored layers (derive_bt_from_layers), so the corrector's
         ! btstep also starts from the step-entry state, as MOM6's does.
         if (dyn%split_scheme == SPLIT_SCHEME_PRED_CORR .and. stage == 1) then
            call restore_state(ms)
         end if
      end do

      if (dyn%split_scheme == SPLIT_SCHEME_PRED_CORR) then
         ! No SSP average: the corrector's single full-dt update IS the
         ! step (SPEC §1 fact 1, §5 trap 4 — averaging FB stages
         ! annihilates the internal wave, |R| = 0.089 at θ = 1.35).
      else
         call rk2_average(ms)
         if (dyn%check_h_positive) then
            call check_h_positive_or_die(grid, ms, "after rk2_average", 3, &
                                         dyn%outer_step_count + 1, check_layers=.true.)
         end if
         ! Land-face velocity reset on the averaged state (C4 / R4b.2).
         call mask_layer_velocities(grid, metrics, ms, bt_work=dyn%bt_work)
         ! Phase-2 vanished-layer velocity reset on the RK2-averaged state.
         ! Only for VCOORD_LAGRANGIAN; no-op otherwise (bit-identical).
         if (dyn%reset_vanished_u .and. present(vcoord)) then
            if (vcoord%coord_type == VCOORD_LAGRANGIAN) then
               call reset_vanished_layer_velocities(ms, &
                                                    isopycnal_vanish_tol(dyn%angstrom_h, pd_floor=ct%positive_definite))
            end if
         end if
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               call rk2_average_field_3d(ms%tracers(it)%hTr0, ms%tracers(it)%hTr, &
                                         size(ms%tracers(it)%hTr, 1), &
                                         size(ms%tracers(it)%hTr, 2), &
                                         size(ms%tracers(it)%hTr, 3))
            end do
         end if
      end if
      call probe_dS(grid, ms, "after rk2_average", 3, dyn%outer_step_count + 1)

      ! Phase 2 (6b): advance the windowed-advect clock once per OUTER step
      ! (not per RK2 stage) when accumulating.  Reset to 0 by the drain.
      if (dyn%dt_tracer_advect_ratio > 1) ct%t_dyn_rel_adv = ct%t_dyn_rel_adv + dt

      ! Velocity housekeeping (E7): advective-CFL truncation then the
      ! absolute maxvel cap.  Applied after the RK2 average but before
      ! the ALE remap — so the remap sees a bounded velocity field.
      ! Both no-op when their knob <= 0.
      ! Phase-3: pass vanish_tol when cfl_ignore_vanished + VCOORD_LAGRANGIAN.
      if (dyn%cfl_ignore_vanished .and. present(vcoord)) then
         if (vcoord%coord_type == VCOORD_LAGRANGIAN) then
            call apply_velocity_truncation(ms, metrics, dt, dyn%cfl_trunc, dyn%maxvel, &
                                           dyn%ntrunc_step, &
                                           vanish_tol=isopycnal_vanish_tol(dyn%angstrom_h, &
                                                                           pd_floor=ct%positive_definite), &
                                           n_nanzero=dyn%n_nanzero_step, &
                                           clip_cell_metric=.true., &
                                           nan_i=nan_i, nan_j=nan_j, nan_k=nan_k, nan_is_u=nan_is_u, &
                                           nan_dx=nan_dx, nan_visc_cfl=nan_visc_cfl, nu_h=hv%nu_h)
         else
            call apply_velocity_truncation(ms, metrics, dt, dyn%cfl_trunc, dyn%maxvel, &
                                           dyn%ntrunc_step, n_nanzero=dyn%n_nanzero_step, &
                                           clip_cell_metric=.true., &
                                           nan_i=nan_i, nan_j=nan_j, nan_k=nan_k, nan_is_u=nan_is_u, &
                                           nan_dx=nan_dx, nan_visc_cfl=nan_visc_cfl, nu_h=hv%nu_h)
         end if
      else
         call apply_velocity_truncation(ms, metrics, dt, dyn%cfl_trunc, dyn%maxvel, &
                                        dyn%ntrunc_step, n_nanzero=dyn%n_nanzero_step, &
                                        clip_cell_metric=.true., &
                                        nan_i=nan_i, nan_j=nan_j, nan_k=nan_k, nan_is_u=nan_is_u, &
                                        nan_dx=nan_dx, nan_visc_cfl=nan_visc_cfl, nu_h=hv%nu_h)
      end if
      dyn%ntrunc_total = dyn%ntrunc_total + dyn%ntrunc_step
      ! Loud NaN-catch accounting: any non-zero is a producer 0/0 upstream.
      ! Actionable (FINDINGS.md, the global-tripolar-aquaplanet debugging
      ! session that had nothing but this bare count to go on): name the
      ! first non-finite face's location, its local cell size, and — since
      ! `hv%nu_h` was cheap to thread through — the local viscous CFL a
      ! plain scalar-nu_h Laplacian would see there, so an under-resolved
      ! viscous CFL (the actual root cause that day) is visible immediately
      ! instead of requiring a separate bisection.
      if (dyn%n_nanzero_step > 0) then
         dyn%n_nanzero_total = dyn%n_nanzero_total + int(dyn%n_nanzero_step, int64)
         write (output_unit, '("[nan-catch] outer step ", i0, ": zeroed ", i0, &
            &" non-finite face velocities; first at ", a, "-face (i=", i0, ", j=", i0, &
            &", k=", i0, "), local dx=", es10.3, " m, local nu_h*dt/dx^2=", es10.3, &
            &" (nu_h=", es10.3, " m2/s)")') &
            dyn%outer_step_count + 1, dyn%n_nanzero_step, merge("u", "v", nan_is_u), &
            nan_i, nan_j, nan_k, nan_dx, nan_visc_cfl, hv%nu_h
         flush (output_unit)
      end if
      ! Outer-step seam (stage 0): the RK2-averaged, truncation-bounded
      ! state the ALE remap + next step will consume.
      call chksum_state(grid, ms, dyn%chksum_probe, "post_trunc", 0, dyn%outer_step_count + 1)

      ! Phase 2 (6b) windowed horizontal tracer-advect DRAIN.  Only active
      ! when `dt_tracer_advect_ratio > 1`; the every-step path (ratio = 1)
      ! never accumulated, so this is a hard no-op there (bit-identical).
      ! Runs AFTER the RK2 average (hTr is the window-start frozen mass,
      ! h_layer is the RK2-averaged window-end thickness) and BEFORE the
      ! ALE remap (spec §(c) ordering: advect-then-reset → remap).
      !
      ! Fires when the window just filled (`is_tracer_advect_step`), OR as
      ! a mandatory FLUSH before any ALE remap if the window is non-empty
      ! (`t_dyn_rel_adv > 0`) — so no Lagrangian-then-remapped state is ever
      ! emitted with un-drained accumulated transport.  The integer-multiple
      ! config constraint makes the natural window-close coincide with the
      ! thermo step, so the flush path is the safety net, not the norm.
      ! TODO(flush): wire the same drain into the output/restart cadence in
      ! the driver for the end-of-run segment (here it is remap-aligned).
      if (dyn%dt_tracer_advect_ratio > 1) then
         if (dyn%is_tracer_advect_step() .or. &
             (present(vcoord) .and. dyn%is_thermo_step() .and. ct%t_dyn_rel_adv > 0.0_wp)) then
            call profiler_start("ocean_tracer_drain")
            if (present(bc)) then
               call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio, bc=bc)
            else
               call continuity_tracer_drain(grid, metrics, ct, ms, dyn%dt_tracer_advect_ratio)
            end if
            call profiler_stop("ocean_tracer_drain")
            call probe_dS(grid, ms, "after tracer drain", 3, dyn%outer_step_count + 1)
         end if
      end if

      ! Lagrangian-then-remap: dynamics has advanced on the existing
      ! grid; now relayer onto vcoord%target_h.  Orchestrator no-ops
      ! for VCOORD_EULERIAN_Z, so the call is bit-safe when vcoord
      ! is in its default configuration.
      !
      ! Gated on the THERMO cadence (`is_thermo_step()`): with the default
      ! `dt_therm_ratio = 1` every outer step is a thermo step, so the
      ! remap fires every step exactly as before (bit-identical).  With
      ! `dt_therm_ratio > 1` the layers run Lagrangian (h_layer is
      ! prognostic; continuity advances it every step) and relayer onto
      ! `target_h` once per thermo interval — the MOM6 DT_THERM design
      ! (Adcroft & Hallberg 2006).  The Lagrangian state between remaps is
      ! self-consistent: diag output remaps to fixed-z independent of the
      ! layer grid, and restart checkpoints carry `outer_step_count` so a
      ! mid-Lagrangian resume re-aligns the cadence bit-exactly.
      if (present(vcoord) .and. dyn%is_thermo_step()) then
         call profiler_start("ocean_ale_remap")
         ! The grid time-filter relaxes once per thermo interval, so it
         ! must see the aggregate thermo dt (= dt at dt_therm_ratio=1, the
         ! bit-identity default).
         call ocean_apply_ale_remap_step(grid, vcoord, ms, dyn%bt_work%bt_eta, dyn%bt_work%bt_H_ref, &
                                         method=vcoord%remap_method, eos=eos, dt=dyn%therm_dt(dt))
         call profiler_stop("ocean_ale_remap")
         call probe_dS(grid, ms, "after ALE remap", 3, dyn%outer_step_count + 1)
      end if

      ! Ideal-age surface reset (PR-7): the Dirichlet BC `age = young_val`
      ! on k=nz must be the LAST operator to touch the top layer this outer
      ! step — applied per-RK2-stage it is halved by `rk2_average_field_3d`
      ! above (a Dirichlet condition inside an averaged sub-step is not a
      ! Dirichlet condition); applied before the ALE remap it is overwritten
      ! by the remap's vertical redistribution of subsurface age into the
      ! new top cell. So: after rk2_average, after continuity_tracer_drain,
      ! after the ALE remap, before the ghost-SSH refill (SSH-only, doesn't
      ! touch tracers). Thermo-cadence gated to compose with the remap;
      ! self-gates internally on `ms%idx_age <= 0`.
      if (dyn%is_thermo_step()) then
         t_now = 0.0_wp
         if (present(t)) t_now = t
         call ocean_ideal_age_reset_surface(grid, ms, &
                                            ocean_ideal_age_young_val(dyn%ideal_age_young_val, &
                                                                      dyn%ideal_age_sfc_growth_rate, t_now))
      end if

      ! Re-establish the zero-gradient free surface in the open-edge GHOST
      ! columns before the diagnostic manager reads the state.  The slow
      ! continuity drifts the ghost h_layer (array-edge fluxes) and the
      ! conservative remap preserves that drift, so `SSH = Σh_layer − b`
      ! blows up in the halo (worst at open corners) while the physical
      ! interior is healthy.  Physics-neutral: next step's stage-1
      ! `ocean_obc_fill_ghosts` re-fills these ghosts before any kernel
      ! reads them.  No-op for WALL/PERIODIC edges (bit-identical).
      if (present(bc)) call ocean_obc_refill_ghost_ssh(grid, bc, ms, dyn%bt_work%bt_H_ref)

      dyn%outer_step_count = dyn%outer_step_count + 1
   end subroutine ocean_dyn_step_split

   subroutine run_continuity_chain(grid, metrics, dyn, ct, hd, va, redi, varmix, ms, &
                                   dt, therm_dt, therm_active, is_lagrangian, &
                                   h_min_floor, mass_out_weight, h_only, &
                                   stage_id, step_id, bc, mle, gm)
      !! The slow horizontal continuity + tracer chain (ghost fills →
      !! constrained continuity+tracer split → reservoirs → halo/wrap →
      !! tracer hdiff → Redi → vertical advection), extracted verbatim
      !! from `run_stage_split` so the pred_corr path can run it AFTER the
      !! velocity update + implicit friction (the forward-backward
      !! pairing, SPEC §2 C8/§1 fact 5) while the historical ssp_rk2
      !! path keeps it before the applies (bit-identical).
      !! `h_only` selects the predictor's TR_MODE_NONE (SPEC §2 P9);
      !! `mass_out_weight` is the per-call budget weight (0.5 per SSP
      !! stage; 0 for the discarded predictor state, 1 for the
      !! corrector).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(continuity_t), intent(inout) :: ct
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_redi_t), intent(inout), optional :: redi
      type(ocean_varmix_t), intent(inout), optional :: varmix
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt, therm_dt
      logical, intent(in) :: therm_active, is_lagrangian
      real(wp), intent(in) :: h_min_floor, mass_out_weight
      logical, intent(in) :: h_only
      integer, intent(in) :: stage_id, step_id
      type(ocean_bc_state_t), intent(inout), optional :: bc
      type(ocean_mle_t), intent(inout), optional :: mle
      type(ocean_gm_t), intent(inout), optional :: gm

      integer :: tr_mode
      real(wp) :: h_min_pass

      ! ---- 5. Ghost fills at open edges (before continuity PPM) ----
      ! Fill h_layer and tracer hTr ghost cells at any open-ish edge
      ! using zero-gradient (h) and upwind-aware (hTr) logic.
      ! This generalises the CLAMPED-only ghost fill inside
      ! rdb_continuity::tracer_advect_zonal/meridional to cover OPEN /
      ! TIDAL / CHAPMAN / NESTED as well.
      !
      ! Authority split (no conflict): on CLAMPED edges the old fill inside
      ! rdb_continuity runs LAST (once per Lie-split sub-flux, after h_layer
      ! has been updated mid-split) and is the authority — it unconditionally
      ! clamps the ghost to clamped_tracer*h, matching what the BT mode
      ! imposes.  This fill's CLAMPED writes here are simply overwritten by
      ! it.  On the open-class edges (OPEN/TIDAL/CHAPMAN/NESTED) the old fill
      ! is inert, so THIS fill's upwind-aware values survive and govern.
      ! No-op when bc is absent or all edges are WALL.
      if (present(bc)) call ocean_obc_fill_ghosts(grid, bc, ms)

      ! ---- 5b. Slow horizontal continuity + tracer, constrained ----
      ! Pass `bt_uhbt, bt_vhbt` so the per-layer mass fluxes are
      ! renormalised to vertically sum to the barotropic-substep's transport.
      ! After apply_zonal + apply_meridional, `sum_k(h_layer)` equals
      ! `H_ref + bt_eta_end` to machine precision — no h rescale
      ! needed afterwards.  Tracer rides the same renormalised
      ! fluxes ⇒ per-column `T = hTr/h` stays exact.
      call profiler_start("ocean_continuity")
      ! `mle_fold_active = is_thermo_step()` gates the Fox-Kemper fold to
      ! the thermo cadence: the FK transports are computed once per thermo
      ! interval (mle_compute_transports, above), so folding them on the
      ! intervening non-thermo steps would re-apply stale transports.  At
      ! the default dt_therm_ratio = 1 this is .true. every step ⇒ bit-identical.
      !
      ! Phase 2 (6b) horizontal-tracer-advect cadence: ratio == 1 ⇒
      ! TR_MODE_ADVECT (fused every-step advect, bit-identical); ratio > 1
      ! ⇒ TR_MODE_ACCUMULATE (advance h + accumulate ½·flux·dt into
      ! ct%uhtr/vhtr per RK2 stage, hTr frozen).  The boundary drain runs
      ! once per window in ocean_dyn_step_split.
      if (h_only) then
         ! pred_corr PREDICTOR (SPEC §2 P9): continuity advances h (into the
         ! provisional hp, discarded after the stage) and produces u_av via
         ! u_cor, but must not touch tracers or the advect window.
         tr_mode = TR_MODE_NONE
      else if (dyn%dt_tracer_advect_ratio <= 1) then
         tr_mode = TR_MODE_ADVECT
      else
         tr_mode = TR_MODE_ACCUMULATE
      end if
      ! Conservative min-thickness mode: pass h_min=0 so continuity performs
      ! the RAW (non-injecting) h-update; the sub-floor layers are then repaired
      ! conservatively by ocean_apply_conservative_min_thickness below.  The
      ! legacy injecting floor stays the default (h_min=h_min_floor).
      h_min_pass = h_min_floor
      if (ct%conservative_floor) h_min_pass = 0.0_wp
      ! MOM6 time-mean fields (SPEC §2 C7/C9): stash h_in in h_av, average with
      ! h_out after continuity.  NOTE: under the present two-stage scheme this
      ! yields a PER-STAGE mean, not the per-step mean MOM6 forms; it becomes
      ! exact when S4 restructures the outer loop.  Nothing reads it until S3.
      call copy_field_3d(ms%h_layer, ms%h_av_layer, &
                         size(ms%h_layer, 1), size(ms%h_layer, 2), size(ms%h_layer, 3))
      ! SPEC S2b (`&ocean_bt_nml renorm_visc_rem`): forward `visc_rem_u/v`
      ! so the transport-matching inversion runs in the γ-weighted MOM6
      ! form (`u_cor = u + du·γ_k`, Jacobian `dy·h_marg·γ_k`) — a
      ! heavily-frictioned layer receives a smaller share of `du`.  Off ⇒
      ! the historical uniform-`du` renormaliser, bit-identical.
      if (dyn%bt_work%bt_renorm_visc_rem) then
         if (present(bc)) then
            call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                              uhbt=dyn%bt_work%bt_uhbt, &
                                              vhbt=dyn%bt_work%bt_vhbt, bc=bc, mle=mle, &
                                              mle_fold_active=dyn%is_thermo_step(), &
                                              tracer_mode=tr_mode, gm=gm, &
                                              h_min=h_min_pass, &
                                              visc_rem_u=dyn%bt_work%visc_rem_u, &
                                              visc_rem_v=dyn%bt_work%visc_rem_v, &
                                              u_cor=ms%u_av_layer, v_cor=ms%v_av_layer)
         else
            call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                              uhbt=dyn%bt_work%bt_uhbt, &
                                              vhbt=dyn%bt_work%bt_vhbt, mle=mle, &
                                              mle_fold_active=dyn%is_thermo_step(), &
                                              tracer_mode=tr_mode, gm=gm, &
                                              h_min=h_min_pass, &
                                              visc_rem_u=dyn%bt_work%visc_rem_u, &
                                              visc_rem_v=dyn%bt_work%visc_rem_v, &
                                              u_cor=ms%u_av_layer, v_cor=ms%v_av_layer)
         end if
      else if (present(bc)) then
         call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                           uhbt=dyn%bt_work%bt_uhbt, &
                                           vhbt=dyn%bt_work%bt_vhbt, bc=bc, mle=mle, &
                                           mle_fold_active=dyn%is_thermo_step(), &
                                           tracer_mode=tr_mode, gm=gm, &
                                           h_min=h_min_pass, &
                                           u_cor=ms%u_av_layer, v_cor=ms%v_av_layer)
      else
         call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                           uhbt=dyn%bt_work%bt_uhbt, &
                                           vhbt=dyn%bt_work%bt_vhbt, mle=mle, &
                                           mle_fold_active=dyn%is_thermo_step(), &
                                           tracer_mode=tr_mode, gm=gm, &
                                           h_min=h_min_pass, &
                                           u_cor=ms%u_av_layer, v_cor=ms%v_av_layer)
      end if
      ! MOM6 C9: h_av = 0.5*(h_in + h_out).
      call rk2_average_field_3d(ms%h_layer, ms%h_av_layer, &
                                size(ms%h_layer, 1), size(ms%h_layer, 2), size(ms%h_layer, 3))
      ! Conservative minimum-thickness borrow (isopycnal grounding stability).
      ! Only when the knob is on AND we have a positive floor (⇒ VCOORD_LAGRANGIAN
      ! by the setup fail-loud guard).  No-op on any ungrounded column ⇒ the
      ! isopycnal interface structure is untouched where all layers meet floor.
      ! NOT checked under `conservative_floor`: that mode deliberately passes
      ! `h_min = 0`, and the three clamp sites are gated on `h_min_use > 0`, so
      ! the h-update here is genuinely RAW — a transiently negative layer is
      ! legal and is repaired by the borrow below (which flags any
      ! `h < h_floor`, negatives included, and inflates it back). Checking here
      ! would abort on the first benign intermediate. With the injecting floor
      ! the clamp DOES hold, so a negative there is a real defect.
      if (dyn%check_h_positive) then
         call check_h_positive_or_die(grid, ms, "after continuity_tracer_step_split", &
                                      0, dyn%outer_step_count + 1, &
                                      check_layers=.not. ct%conservative_floor)
      end if
      if (ct%conservative_floor .and. h_min_floor > 0.0_wp) then
         call profiler_start("ocean_min_thickness")
         call ocean_apply_conservative_min_thickness(grid, ms, ct%mt_h_new%data, &
                                                     ct%mt_grounded%data, h_min_floor)
         call profiler_stop("ocean_min_thickness")
         if (dyn%check_h_positive) then
            call check_h_positive_or_die(grid, ms, "after conservative_min_thickness", &
                                         0, dyn%outer_step_count + 1, check_layers=.true.)
         end if
      end if
      ! Mass budget: flux_h_layer now holds the total horizontal divergence;
      ! accumulate the boundary outflux for this RK2 stage (weight 0.5).
      call ocean_accumulate_mass_out(ms, ms%flux_h_layer, metrics%areaT, &
                                     grid%nghost, dt, mass_out_weight)
      call profiler_stop("ocean_continuity")
      call probe_dS(grid, ms, "after continuity_tracer_split", stage_id, step_id)

      ! Reservoir update (§1, v2): evolve tres toward T_int / T_data using
      ! the stage's wall-face mass_flux_*_layer values.  Called while those
      ! arrays still hold this stage's fluxes (before any overwrite).
      ! No-op when res_lscale_out == 0 and res_lscale_in == 0 (default).
      if (present(bc)) call ocean_obc_update_reservoirs(grid, bc, ms, dt)

      ! Post-continuity periodic wrap (design §1.5): re-wrap h_layer and
      ! tracers after the Lie-split continuity step and before hdiff, so
      ! the horizontal diffusion kernel reads correct ghost-zone values.
      ! No-op when neither periodic axis is set.
      ! O2: unconditional multi-rank ghost exchange (D0 — no-op on 1 rank).
      call ocean_halo_exchange_ml_state(ms)
      if (present(bc)) call ocean_periodic_wrap_state(grid, bc, ms, &
                                                      skip_x=ocean_halo_is_decomposed_x(), skip_y=ocean_halo_is_decomposed_y())
      ! Tripolar north-fold after the post-continuity periodic wrap, so
      ! hdiff reads fold-consistent ghost values at the seam.
      if (present(bc)) call ocean_fold_wrap_state(grid, bc, ms)

      call profiler_start("ocean_tracer_hdiff")
      call tracer_hdiff(grid, metrics, hd, ms, therm_dt, active=therm_active, bc=bc)
      call profiler_stop("ocean_tracer_hdiff")
      call probe_dS(grid, ms, "after tracer_hdiff", stage_id, step_id)

      ! Redi neutral diffusion (capability [3]) — Phase B.  AUGMENTS the
      ! along-coordinate `tracer_hdiff` above with the rotated (along-
      ! isopycnal) flux, using the Phase-A coefficients built once this
      ! outer step.  THERMO cadence (self-gated inside on `therm_active`
      ! via the same `is_thermo_step` window the coeffs were built under).
      ! No-op when absent / disabled / khtr<=0.
      if (present(redi) .and. therm_active) then
         call profiler_start("ocean_redi")
         ! VarMix seam: feed Redi the spatially-varying KhTr when VarMix is
         ! enabled; otherwise Redi falls back to its scalar khtr.
         if (present(varmix)) then
            if (varmix%enable) then
               call redi_apply_flux(grid, metrics, redi, ms, therm_dt, &
                                    khtr_u_ext=varmix%khtr_u, khtr_v_ext=varmix%khtr_v, bc=bc)
            else
               call redi_apply_flux(grid, metrics, redi, ms, therm_dt, bc=bc)
            end if
         else
            call redi_apply_flux(grid, metrics, redi, ms, therm_dt, bc=bc)
         end if
         call profiler_stop("ocean_redi")
      end if
      ! ---- Vertical advection (Eulerian-z only) ----
      ! In Eulerian-z mode the vertical w-divergence cancels the
      ! horizontal flux divergence per layer, pinning `h_layer` to
      ! reference.  In Lagrangian mode (`vcoord` present and
      ! coord_type /= EULERIAN_Z) we skip this — h_layer is allowed
      ! to evolve under the MOM6-constrained horizontal continuity
      ! alone, and the ALE remap at end of outer step relayers
      ! conservatively onto `vcoord%target_h`.  Skipping here
      ! removes the surface F(nz+1) = 0 vs w(nz+1) ≠ 0 CWC
      ! inconsistency that seeds the salt-baroclinic instability
      ! under realistic β_S.
      if (.not. is_lagrangian) then
         call profiler_start("ocean_vertical_advect")
         call compute_w_from_continuity(grid, va, ms)
         call probe_dS(grid, ms, "after compute_w_from_cont", stage_id, step_id)
         call tracer_advect_vertical(grid, va, ms, therm_dt, active=therm_active)
         call probe_dS(grid, ms, "after tracer_advect_vertical", stage_id, step_id)
         call profiler_stop("ocean_vertical_advect")
      end if
   end subroutine run_continuity_chain

   subroutine run_meke_step(grid, metrics, gm, varmix, wavespeed, hv, ms, dt, meke)
      !! Thin dispatcher: call `meke_step` only when the MEKE slot is present
      !! AND enabled, forwarding the (also-optional) VarMix + wavespeed slots
      !! so the feedback seam + length scales engage when those are on.
      !! `meke_step` itself no-ops on `enable=.false.`; this guard avoids the
      !! call (and the present-propagation noise) when the slot is absent.
      !! `hv%ke_diss` (the lateral-viscosity KE dissipation rate) is forwarded
      !! for the frictional source; it is 0 unless `hv%compute_ke_diss` is set
      !! (and `meke%frcoeff<0` ignores it) ⇒ inert by default.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_gm_t), intent(in) :: gm
      type(ocean_varmix_t), intent(inout), optional :: varmix
      type(ocean_wave_speed_t), intent(in), optional :: wavespeed
      type(ocean_horizontal_viscosity_t), intent(in) :: hv
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
      type(ocean_meke_t), intent(inout), optional :: meke

      if (.not. present(meke)) return
      if (.not. meke%enable) return
      call meke_step(grid, metrics, meke, gm, varmix, wavespeed, ms, dt, &
                     ke_diss_ext=hv%ke_diss)
   end subroutine run_meke_step

   pure function gm_refreshes_varmix(gm, slopes, dyn) result(refreshes)
      !! `.true.` when the GM block below this step already calls
      !! `varmix_compute` (GM present + enabled + slopes present, at a
      !! thermo step) — so the standalone Gap-1 resolution-function refresh
      !! must NOT run it again (avoids a double compute and preserves the
      !! GM+VarMix slope ordering bit-for-bit).
      type(ocean_gm_t), intent(in), optional :: gm
      type(ocean_slopes_t), intent(in), optional :: slopes
      type(ocean_dyn_t), intent(in) :: dyn
      logical :: refreshes
      refreshes = .false.
      if (.not. present(gm)) return
      if (.not. present(slopes)) return
      if (.not. gm%enable) return
      refreshes = dyn%enable_thermodynamics .and. dyn%is_thermo_step()
   end function gm_refreshes_varmix

   pure function lateral_mix_uses_resoln(lateral_mix, varmix) result(uses)
      !! `.true.` when the lateral-mix compute call should be handed the
      !! VarMix resolution-function face fields (Gap 1): the lateral-mix
      !! slot is present + initialised with `resoln_scaled_visc`, AND VarMix
      !! is present + enabled (so `res_fn_u/v` carry a valid, up-to-date
      !! resolution function).  Either absent / off ⇒ unscaled coefficients
      !! (bit-identical).
      type(ocean_lateral_mix_t), intent(in), optional :: lateral_mix
      type(ocean_varmix_t), intent(in), optional :: varmix
      logical :: uses
      uses = .false.
      if (.not. present(lateral_mix)) return
      if (.not. present(varmix)) return
      if (.not. lateral_mix%is_init) return
      uses = lateral_mix%resoln_scaled_visc .and. varmix%enable
   end function lateral_mix_uses_resoln

   subroutine check_h_positive_or_die(grid, ms, label, stage, outer_step, check_layers)
      !! `&ocean_isopycnal_nml check_h_positive` guard.  Abort on the FIRST
      !! negative layer thickness, naming the pipeline stage that produced it,
      !! the offending `(i,j,k)`, and that column's full thickness profile.
      !!
      !! Why this exists: a negative `h` surfaces much later and somewhere else
      !! — as a console-stats NaN, or as an hourly diagnostic minimum — by
      !! which point the producing kernel cannot be identified.  Under
      !! `conservative_floor` it should not be reachable at all: continuity is
      !! handed `h_min = 0`, so its `max(h - dt*div, h_min)` clamps at zero,
      !! and the conservative borrow only redistributes WITHIN a column.  A
      !! trip here therefore localises a real defect rather than a tuning
      !! problem.
      !!
      !! Cheap on the healthy path: one device-side min-reduction, no H<-D
      !! copy.  The host walk + column dump run only when already aborting.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      character(len=*), intent(in) :: label
         !! Pipeline stage that last wrote `h_layer` (e.g. "after continuity").
      integer, intent(in) :: stage, outer_step
      logical, intent(in) :: check_layers
         !! `.true.` => also abort on a negative single LAYER.  Pass `.false.`
         !! between the raw continuity update and the conservative borrow,
         !! where a transiently negative layer is legal; the column-total
         !! check runs unconditionally either way.
      integer :: i, j, k, ig, i0, i1, j0, j1, ib, jb, kb
      real(wp) :: h_min, col_sum, col_min, worst, acc

      ig = grid%nghost
      i0 = ig + 1
      i1 = grid%nx_total - ig
      j0 = ig + 1
      j1 = grid%ny_total - ig

      ! Column TOTAL is the invariant that always holds.  A negative single
      ! LAYER is legal between the raw continuity update and the conservative
      ! borrow (see the call-site comment), but no column may ever reach a
      ! non-positive total — and if one does, `min_thickness_target_column`'s
      ! degenerate branch (`total <= nz*h_floor`) spreads that negative total
      ! uniformly over every layer, poisoning the whole column silently.
      ! Catching the total here names the stage that evacuated the column.
      col_min = huge(1.0_wp)
      !$acc parallel loop collapse(2) reduction(min:col_min) &
      !$acc   private(k, col_sum) present(ms%h_layer)
      do j = j0, j1
         do i = i0, i1
            col_sum = 0.0_wp
            do k = 1, ms%nz_ml
               col_sum = col_sum + ms%h_layer(i, j, k)
            end do
            col_min = min(col_min, col_sum)
         end do
      end do

      h_min = huge(1.0_wp)
      if (check_layers) then
         !$acc parallel loop collapse(3) reduction(min:h_min) present(ms%h_layer)
         do k = 1, ms%nz_ml
            do j = j0, j1
               do i = i0, i1
                  h_min = min(h_min, ms%h_layer(i, j, k))
               end do
            end do
         end do
      end if
      if (h_min >= 0.0_wp .and. col_min > 0.0_wp) return

      ! Abort path.  Pull the COMPONENT array, never the aggregate `ms` — an
      ! aggregate D->H copy overwrites the host allocatable descriptors with
      ! DEVICE addresses and the next host read segfaults.
      !$acc update self(ms%h_layer)

      ! Locate the worst column: if a column total went non-positive that is
      ! the primary failure, so report the most-negative TOTAL.  Otherwise
      ! report the column holding the most-negative single layer.
      ib = i0
      jb = j0
      worst = huge(1.0_wp)
      if (col_min <= 0.0_wp) then
         do j = j0, j1
            do i = i0, i1
               acc = 0.0_wp
               do k = 1, ms%nz_ml
                  acc = acc + ms%h_layer(i, j, k)
               end do
               if (acc < worst) then
                  worst = acc
                  ib = i
                  jb = j
               end if
            end do
         end do
      else
         do k = 1, ms%nz_ml
            do j = j0, j1
               do i = i0, i1
                  if (ms%h_layer(i, j, k) < worst) then
                     worst = ms%h_layer(i, j, k)
                     ib = i
                     jb = j
                  end if
               end do
            end do
         end do
      end if
      ! Thinnest layer within the reported column.
      kb = 1
      do k = 1, ms%nz_ml
         if (ms%h_layer(ib, jb, k) < ms%h_layer(ib, jb, kb)) kb = k
      end do

      write (output_unit, "(a)") repeat("=", 68)
      if (col_min <= 0.0_wp) then
         write (output_unit, '("[h-guard] NON-POSITIVE COLUMN TOTAL (column evacuated)")')
      else
         write (output_unit, '("[h-guard] NEGATIVE LAYER THICKNESS")')
      end if
      write (output_unit, '("[h-guard]   stage      : ", a)') trim(label)
      write (output_unit, '("[h-guard]   outer step : ", i0, "   rk2 stage : ", i0)') &
         outer_step, stage
      write (output_unit, '("[h-guard]   at (i,j,k) : ", i0, ", ", i0, ", ", i0, &
                          &"   (physical i0/j0 = ", i0, "/", i0, ")")') &
         ib, jb, kb, i0, j0
      write (output_unit, '("[h-guard]   thinnest k : ", i0, "   h = ", es13.5)') &
         kb, ms%h_layer(ib, jb, kb)
      write (output_unit, '("[h-guard]   min layer h over interior  : ", es13.5)') h_min
      write (output_unit, '("[h-guard]   min column total over interior : ", es13.5)') col_min
      col_sum = 0.0_wp
      do k = 1, ms%nz_ml
         col_sum = col_sum + ms%h_layer(ib, jb, k)
      end do
      write (output_unit, '("[h-guard]   column sum : ", es13.5)') col_sum
      write (output_unit, '("[h-guard]   column profile (k, h):")')
      do k = 1, ms%nz_ml
         write (output_unit, '("[h-guard]     ", i4, "  ", es14.6)') k, ms%h_layer(ib, jb, k)
      end do
      write (output_unit, "(a)") repeat("=", 68)
      flush (output_unit)
      error stop "h-guard: negative layer thickness (see [h-guard] block above)"
   end subroutine check_h_positive_or_die

   subroutine probe_dS(grid, ms, label, stage, outer_step)
      !! Debug-gated diagnostic.  Pulls `h_layer` + `hTr_S` from device
      !! and prints `max|hTr_S/h_layer - bcdiag_S_ref|` over interior
      !! columns.  Early-returns when `bcdiag_enabled = .false.` or
      !! after `bcdiag_step_limit` outer steps have completed.
      !!
      !! Used to localise where the first 9e-5 dS jump enters during
      !! one outer step of `geostrophic_adjust` under realistic β_S.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      character(len=*), intent(in) :: label
      integer, intent(in) :: stage, outer_step

      integer :: i, j, k, iS, ig, i0, i1, j0, j1
      real(wp) :: max_dS, dS

      if (.not. bcdiag_enabled) return
      if (outer_step > bcdiag_step_limit) return
      if (ms%idx_salinity <= 0) return
      if (.not. allocated(ms%tracers)) return

      iS = ms%idx_salinity
      !$acc update self(ms%h_layer, ms%tracers(iS)%hTr)

      ig = grid%nghost
      i0 = ig + 1
      i1 = grid%nx_total - ig
      j0 = ig + 1
      j1 = grid%ny_total - ig

      max_dS = 0.0_wp
      do k = 1, ms%nz_ml
         do j = j0, j1
            do i = i0, i1
               if (ms%h_layer(i, j, k) > 0.0_wp) then
                  dS = abs(ms%tracers(iS)%hTr(i, j, k)/ms%h_layer(i, j, k) - bcdiag_S_ref)
                  if (dS > max_dS) max_dS = dS
               end if
            end do
         end do
      end do

      write (output_unit, '("[bcdiag] step=", i0, " stage=", i0, " ", a40, " max|dS|=", es12.4)') &
         outer_step, stage, label, max_dS
      flush (output_unit)
   end subroutine probe_dS

   subroutine probe_h_vs_eta_residual(grid, ms, bt_work, stage, outer_step)
      !! Debug probe: print max|sum_k(h_layer) - (H_ref + bt_eta_end)|
      !! over interior cells.  In Eulerian-z this is forced to zero
      !! by `apply_bt_correction`'s h-rescale.  In Lagrangian mode
      !! the rescale is skipped, so this residual is what the slow
      !! continuity actually drifts to (expected to be FP).  Watching
      !! how it grows over time across many outer steps is what
      !! pins down whether the long-run momentum NaN is FP
      !! accumulation in `sum_k(h_layer) - H - bt_eta_end`.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(in) :: bt_work
      integer, intent(in) :: stage, outer_step

      integer :: i, j, k, i0, i1, j0, j1, ig
      real(wp) :: max_res, sum_h, target_h, res

      if (.not. bcdiag_enabled) return
      if (outer_step > bcdiag_step_limit) return

      !$acc update self(ms%h_layer, bt_work%bt_eta_end, bt_work%bt_H_ref)

      ig = grid%nghost
      i0 = ig + 1
      i1 = grid%nx_total - ig
      j0 = ig + 1
      j1 = grid%ny_total - ig

      max_res = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            sum_h = 0.0_wp
            do k = 1, ms%nz_ml
               sum_h = sum_h + ms%h_layer(i, j, k)
            end do
            target_h = bt_work%bt_H_ref(i, j) + bt_work%bt_eta_end(i, j)
            res = abs(sum_h - target_h)
            if (res > max_res) max_res = res
         end do
      end do

      write (output_unit, &
             '("[bcdiag] step=", i0, " stage=", i0, " ", a40, " max|h_sum - (H+eta_end)|=", es12.4)') &
         outer_step, stage, "h_vs_eta_residual", max_res
      flush (output_unit)
   end subroutine probe_h_vs_eta_residual

   subroutine run_stage_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                              va, hd, vd, vmix, ms, dt, n_inner, sf, geo, stage, vcoord, bc, sp, t, &
                              lateral_mix, epbl, kshear, mle, gm, redi, varmix, vmix_tidal, meke, &
                              eta_forcing)
      !! One FE stage of the split-explicit step.  See the
      !! `ocean_dyn_step_split` header for the design.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics — forwarded to the PGF (and,
         !! in slice 2, the other geometry-aware kernels).
      type(ocean_dyn_t), intent(inout) :: dyn
      type(eos_t), intent(in) :: eos
      type(coriolis_adv_t), intent(inout) :: cor
      type(continuity_t), intent(inout) :: ct
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      integer, intent(in) :: n_inner
      type(ocean_surface_flux_t), intent(in), optional :: sf
      type(ocean_geothermal_t), intent(in), optional :: geo
         !! Geothermal bottom-heat-flux slot.  See `ocean_dyn_step_split`.
      integer, intent(in), optional :: stage
         !! Outer SSP-RK2 stage index (1 or 2) — only used by the
         !! debug probe so trace lines self-identify which half-step
         !! they came from.  Production paths leave it unset.
      type(ocean_vcoord_t), intent(in), optional :: vcoord
         !! Vertical-coordinate state.  When present and
         !! `coord_type /= VCOORD_EULERIAN_Z` the stage runs in
         !! Lagrangian mode: vertical advection is skipped (no
         !! `compute_w_from_continuity` / `tracer_advect_vertical`)
         !! and `apply_bt_correction`'s h-rescale is skipped.  The
         !! MOM6-constrained slow continuity makes `sum_k(h_layer)`
         !! self-consistent with the barotropic-substep η; the ALE remap at
         !! the end of the outer step relayers (h, hTr) onto
         !! `vcoord%target_h` conservatively.  Absent or EULERIAN_Z
         !! ⇒ the historical Eulerian-z code path.
      type(ocean_bc_state_t), intent(inout), optional :: bc
         !! Open-boundary config forwarded straight to the barotropic substep.
         !! Absent / all-OBC_WALL keeps the closed-wall path.
      type(ocean_sponge_t), intent(in), optional :: sp
         !! Map-driven sponge slot (PR-23).  See `ocean_dyn_step_split`.
      real(wp), intent(in), optional :: t
         !! Wall-clock time for OBC_TIDAL constituent evaluation.
      type(ocean_lateral_mix_t), intent(inout), optional :: lateral_mix
         !! Flow-aware lateral closure.  See the matching arg on
         !! `ocean_dyn_step_split`.
      type(ocean_epbl_t), intent(inout), optional :: epbl
         !! EPBL slot.  See the matching arg on `ocean_dyn_step_split`.
      type(ocean_kappa_shear_t), intent(inout), optional :: kshear
         !! Kappa-shear slot.  See `ocean_dyn_step_split`.
      type(ocean_tidal_mixing_t), intent(inout), optional :: vmix_tidal
         !! Tidal-mixing slot.  See `ocean_dyn_step_split`.
      type(ocean_mle_t), intent(inout), optional :: mle
         !! Fox-Kemper MLE slot (B5).  Forwarded to
         !! `continuity_tracer_step_split` to fold the precomputed
         !! `uhml`/`vhml` into the mass fluxes.  See `ocean_dyn_step_split`.
      type(ocean_gm_t), intent(inout), optional :: gm
         !! GM thickness-diffusion slot (capability [2]).  Forwarded to
         !! `continuity_tracer_step_split` to fold the precomputed
         !! `uhD`/`vhD` into the mass fluxes.  See `ocean_dyn_step_split`.
      type(ocean_redi_t), intent(inout), optional :: redi
         !! Redi neutral-diffusion slot (capability [3]).  The Phase-A
         !! coefficients are precomputed in `ocean_dyn_step_split`; here
         !! Phase B (`redi_apply_flux`) adds the rotated tracer flux after
         !! the along-coordinate `tracer_hdiff`, at THERMO cadence.
      type(ocean_varmix_t), intent(inout), optional :: varmix
         !! VarMix coefficient slot (capability [4]).  When enabled, its
         !! per-face `khtr_u`/`khtr_v` feed the Redi flux (the Visbeck KhTr
         !! seam); absent / disabled ⇒ Redi uses its scalar `khtr`.
      type(ocean_meke_t), intent(in), optional :: meke
         !! MEKE slot (capability [5]).  Read-only here: when
         !! `meke%backscatter` is on, its `ku` field is injected into the
         !! resolved per-face harmonic viscosity right after
         !! `ocean_lateral_mix_compute` (the energy-return seam, Gap 2).
         !! Absent / `backscatter` off ⇒ bit-identical.
      real(wp), intent(in), optional :: eta_forcing(grid%nx_total, grid%ny_total)
         !! Equilibrium-tide elevation (C1), held static across the inner
         !! substep loop.  Forwarded to `barotropic_substep_nonlinear`'s
         !! PGF; absent ⇒ bit-identical.

      integer :: i, j, nx_face, ny_uface, nx_vface, ny_face
      integer :: stage_id, step_id
      integer :: bc_w_drv, bc_e_drv, bc_s_drv, bc_n_drv
         !! Per-edge BC tags used for the slow-path transport wall-zeroing.
         !! Default OBC_WALL; overridden from bc%<edge>%bc_type when bc is present.
      logical :: has_w_drv, has_e_drv, has_s_drv, has_n_drv
         !! Physical-edge flags for the transport wall-zeroing.  .true. = physical
         !! domain edge (WALL zero applies); .false. = MPI seam (halo owns it, O0).
      logical :: is_lagrangian, therm_active
      logical :: is_pc, is_pred
         !! pred_corr stage roles (SPEC §2/§4 S4): `is_pred` = the predictor
         !! (stage 1 under pred_corr) — provisional `up = u + BE·dt·accel`,
         !! h-only continuity into a discarded hp, coefficients/remnant
         !! only, no tracer physics.  Stage 2 is the corrector: the single
         !! full-dt prognostic update.  Both .false. under ssp_rk2 ⇒ every
         !! gate below is untaken ⇒ bit-identical.
      logical :: sponge_maps_on
         !! .true. when the map-driven sponge (PR-23) supersedes the legacy
         !! band path. A local logical because Fortran does not guarantee
         !! `.and.` short-circuits past `present()`.
      real(wp) :: dt_inner, therm_dt
      real(wp) :: h_min_floor
      real(wp) :: chain_weight, dt_vel

      stage_id = 0
      if (present(stage)) stage_id = stage
      ! Probes label the OUTER step we're INSIDE — i.e. the one
      ! `outer_step_count` increments to at the end of this call.
      step_id = dyn%outer_step_count + 1
      is_lagrangian = .false.
      if (present(vcoord)) is_lagrangian = vcoord%coord_type /= VCOORD_EULERIAN_Z
      ! Phase-1 Lagrangian h-floor: only non-zero for VCOORD_LAGRANGIAN
      ! specifically (NOT the broad is_lagrangian = /= EULERIAN_Z — R3).
      h_min_floor = 0.0_wp
      if (present(vcoord)) then
         if (vcoord%coord_type == VCOORD_LAGRANGIAN) h_min_floor = ct%angstrom_h
      end if
      is_pc = dyn%split_scheme == SPLIT_SCHEME_PRED_CORR
      is_pred = is_pc .and. stage_id == 1
      ! Predictor: no tracer physics of any kind (MOM6's predictor never
      ! touches thermodynamics — SPEC §2; the corrector owns the full
      ! tracer chain at the thermo cadence).
      therm_active = dyn%enable_thermodynamics .and. dyn%is_thermo_step() &
                     .and. .not. is_pred
      therm_dt = dyn%therm_dt(dt)
      ! Velocity-apply dt: the predictor's provisional velocity advances
      ! to dt_pred = BE·dt (SPEC §2 P8); everything else (substep,
      ! continuity) stays at the full dt as MOM6 does.
      dt_vel = dt
      if (is_pred) dt_vel = dyn%pc_be*dt
      ! Mass-budget weight for the continuity chain: 0.5 per SSP stage;
      ! under pred_corr the predictor's h advance is DISCARDED (weight 0)
      ! and the corrector's is the whole step (weight 1).
      chain_weight = 0.5_wp
      if (is_pc) chain_weight = merge(0.0_wp, 1.0_wp, is_pred)

      call probe_dS(grid, ms, "entry", stage_id, step_id)

      ! ---- 0. Stage-entry periodic wrap ----
      ! For periodic axes: fill ghost cells of h_layer, u/v layer faces,
      ! and all tracers before derive_bt_from_layers, so every slow
      ! tendency kernel (EOS, Coriolis-adv, PGF, hvisc, bdrag, surface
      ! stress) sees physically-correct values at the seam.  No-op when
      ! neither periodic flag is set.
      ! O2: unconditional multi-rank ghost exchange (D0 — no-op on 1 rank).
      call ocean_halo_exchange_ml_state(ms)
      if (present(bc)) call ocean_periodic_wrap_state(grid, bc, ms, &
                                                      skip_x=ocean_halo_is_decomposed_x(), skip_y=ocean_halo_is_decomposed_y())
      ! Tripolar north-fold seam — periodic-FIRST-fold-SECOND (Appendix A):
      ! folds h/u/v/tracers AFTER the periodic wrap so it reads the already
      ! cyclically-wrapped corner columns.  No-op when bc%north_fold is off.
      if (present(bc)) call ocean_fold_wrap_state(grid, bc, ms)

      ! ---- 0b. pred_corr: the SAME seam fill for the step time-means ----
      ! Under `split_scheme = "pred_corr"` the Coriolis-advection and the
      ! horizontal-viscosity tendencies are evaluated on `u_av`/`v_av`/`h_av`,
      ! NOT on the prognostic fields that step 0 just wrapped.  Those means
      ! come out of the continuity solve (`u_cor`) and the h_in/h_out average,
      ! both of which write the INTERIOR only — so without this their ghost
      ! band holds stale values (zeros on the first step) and the seam
      ! tendency is wrong.
      !
      ! Symptom, and why this is not cosmetic: `periodic_shifted_domain_
      ! identity` (tests/test_ocean_periodic.F90) runs the same physical IC
      ! twice with run B's domain circularly shifted half a period and
      ! demands shift(B) == A BIT-FOR-BIT.  With the ghosts unfilled the
      ! answer depends on where the seam falls: max_diff = 8.6e-07 on h ~ 50
      ! after 8 steps, from exactly 0 under ssp_rk2.  ssp_rk2 never reads
      ! these arrays, so the `is_pc` gate keeps it bit-identical.
      if (is_pc) then
         if (allocated(ms%u_av_layer)) then
            ! Multi-rank seam first, exactly as step 0 does for the
            ! prognostics.  GATED on an actually-decomposed axis: the halo
            ! specifics take EXPLICIT-SHAPE dummies sized from the module's
            ! `oh_nx_total`/`oh_ny_total`, which a unit test that never calls
            ! `ocean_halo_init` leaves at 0 -- a mis-shaped device dummy, and
            ! on the GPU build that reads as garbage rather than as an error.
            if (ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()) then
               call ocean_halo_face_x(ms%u_av_layer, size(ms%u_av_layer, 3))
               call ocean_halo_face_y(ms%v_av_layer, size(ms%v_av_layer, 3))
               call ocean_halo_centre(ms%h_av_layer, size(ms%h_av_layer, 3))
            end if
            if (present(bc)) then
               if (bc%periodic_x .or. bc%periodic_y) then
                  call ocean_periodic_wrap_face_x_3d( &
                     ms%u_av_layer, size(ms%u_av_layer, 1), size(ms%u_av_layer, 2), &
                     size(ms%u_av_layer, 3), grid%nx_phys, grid%ny_phys, grid%nghost, &
                     bc%periodic_x .and. .not. ocean_halo_is_decomposed_x(), &
                     bc%periodic_y .and. .not. ocean_halo_is_decomposed_y())
                  call ocean_periodic_wrap_face_y_3d( &
                     ms%v_av_layer, size(ms%v_av_layer, 1), size(ms%v_av_layer, 2), &
                     size(ms%v_av_layer, 3), grid%nx_phys, grid%ny_phys, grid%nghost, &
                     bc%periodic_x .and. .not. ocean_halo_is_decomposed_x(), &
                     bc%periodic_y .and. .not. ocean_halo_is_decomposed_y())
                  call ocean_periodic_wrap_centre_3d( &
                     ms%h_av_layer, size(ms%h_av_layer, 1), size(ms%h_av_layer, 2), &
                     size(ms%h_av_layer, 3), grid%nx_phys, grid%ny_phys, grid%nghost, &
                     bc%periodic_x .and. .not. ocean_halo_is_decomposed_x(), &
                     bc%periodic_y .and. .not. ocean_halo_is_decomposed_y())
               end if
            end if
         end if
      end if

      ! KE attribution: stage-entry baseline (debug_ke_attr; no-op when off).
      call ke_probe_sample(grid, ms, dyn%ke_probe, "entry", stage_id, step_id)
      call chksum_state(grid, ms, dyn%chksum_probe, "entry", stage_id, step_id)

      ! ---- 1. Snapshot u_bt^n, v_bt^n at start of stage ----
      call derive_bt_from_layers(grid, dyn%bt_work, ms)
      ! Build the per-face upstream-PPM column-sum thickness on the
      ! same snapshot.  No-op when `use_upstream_h_face = .false.`;
      ! otherwise feeds the BT substep + corrector with the same
      ! face-thickness convention slow continuity uses, eliminating
      ! the centred-vs-upstream mismatch at slopes.
      call compute_h_face_upstream(grid, dyn%bt_work, ms)
      ! Build the per-face BT_cont_type flux closure from the same ML
      ! snapshot.  No-op when `use_bt_cont_type = .false.`; otherwise
      ! BTCL_u/v feed the BT substep's flux paths instead of the
      ! naive `uh = u·h_face`.
      call set_local_BT_cont_types(grid, metrics, dyn%bt_work, ms, dt)
      nx_face = size(dyn%bt_work%bt_ubt, 1)
      ny_uface = size(dyn%bt_work%bt_ubt, 2)
      nx_vface = size(dyn%bt_work%bt_vbt, 1)
      ny_face = size(dyn%bt_work%bt_vbt, 2)
      do concurrent(j=1:ny_uface, i=1:nx_face)
         dyn%bt_work%ubt_at_n(i, j) = dyn%bt_work%bt_ubt(i, j)
      end do
      do concurrent(j=1:ny_face, i=1:nx_vface)
         dyn%bt_work%vbt_at_n(i, j) = dyn%bt_work%bt_vbt(i, j)
      end do

      ! Snapshot η before the slow PGF runs — used by the bc-PGF
      ! correction below to form `e_anom = (η_avg − eta_PF)`.
      ! No-op when bt_correction_bc_pgf is off; safe to always run.
      if (dyn%bt_work%bt_correction_bc_pgf) then
         call snapshot_eta_PF(dyn%bt_work)
      end if

      ! ---- 2. Compute all slow tendencies ----
      ! All read h^n / u^n / v^n into per-kernel scratch buffers.
      ! No writes to h_layer, hTr, u_face, or v_face yet.
      call profiler_start("ocean_eos")
      ! See the note at the other ocean_eos_compute call site: the EOS is a
      ! dynamics term and must NOT be held across the thermo window. Note the
      ! predictor exclusion (`.not. is_pred`) also drops away here — the
      ! predictor's PGF needs a current rho just as much as the corrector's.
      call ocean_eos_compute(eos, ms, active=dyn%enable_thermodynamics)
      call profiler_stop("ocean_eos")
      call probe_dS(grid, ms, "after eos", stage_id, step_id)
      call profiler_start("ocean_coriolis_adv")
      ! SPEC S3 (`split_scheme = "pred_corr"`): the Coriolis-advection
      ! tendency reads the `u_av`/`h_av` step time-means, never the
      ! prognostic (MOM6 `CorAdCalc(u_av, v_av, h_av, ...)`).
      ! Default ⇒ prognostic, bit-identical.
      !
      ! Mass-consistent CorAdCalc (`&ocean_coriolis_nml use_state_fluxes`):
      ! in the pred_corr CORRECTOR, `ms%mass_flux_*_layer` still hold the
      ! PREDICTOR chain's renormalised transports — the uh/vh from the same
      ! solve that produced this `u_av` evaluation state.  Consuming them
      ! (instead of the kernel's own `u·h_face` recompute) closes the
      ! evaluate-at-u_av / transport-mismatch energy leak the KE-attribution
      ! meter pinned on coriolis_adv (pdc M16, the rim-Kelvin slow
      ! exponential).  Predictor stages keep the recompute.
      if (dyn%split_scheme == SPLIT_SCHEME_PRED_CORR) then
         call coriolis_adv_compute_tendencies(grid, metrics, cor, ms, &
                                              u_src=ms%u_av_layer, &
                                              v_src=ms%v_av_layer, &
                                              h_src=ms%h_av_layer, &
                                              use_state_fluxes=(cor%state_fluxes &
                                                                .and. stage_id == 2))
      else
         call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
      end if
      call profiler_stop("ocean_coriolis_adv")
      call profiler_start("ocean_pgf")
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      call profiler_stop("ocean_pgf")
      ! bc-PGF correction: compute per-layer `pbce` + face-centred
      ! `gtot_*` now (PGF has populated `pgf%e_face`).  Used later
      ! by `apply_bt_correction` with `use_bc_pgf=.true.`.  Gated;
      ! no-op when bt_correction_bc_pgf is off.
      if (dyn%bt_work%bt_correction_bc_pgf) then
         call compute_pbce(grid, dyn%bt_work, pgf, ms)
         call compute_gtot_faces(grid, dyn%bt_work, ms)
      end if
      call profiler_start("ocean_hvisc")
      ! pred_corr predictor: SKIP the viscous recompute — MOM6's predictor
      ! uses the PREVIOUS step's diffu (SPEC §2 P3, "diffu(u[n-1])"), and
      ! the du_visc/dv_visc buffers persist from the last corrector's C1.
      if (.not. is_pred) then
         ! Flow-aware lateral closure first (Leith / Smagorinsky) so the
         ! hvisc kernel reads the freshly-computed per-face viscosity.
         ! Both calls self-gate on `lateral_mix` absence / `LMIX_NONE`.
         ! Gap 1: when `resoln_scaled_visc` is on AND VarMix is active, forward
         ! the resolution-function face fields so the dynamic coefficients are
         ! scaled by `Res_fn` before the clamps.  Absent / off ⇒ the plain call
         ! ⇒ bit-identical.
         if (lateral_mix_uses_resoln(lateral_mix, varmix)) then
            call ocean_lateral_mix_compute(grid, metrics, lateral_mix, ms, &
                                           res_fn_u=varmix%res_fn_u, &
                                           res_fn_v=varmix%res_fn_v)
         else
            call ocean_lateral_mix_compute(grid, metrics, lateral_mix, ms)
         end if
         ! MEKE backscatter (capability [5], Gap 2): subtract the eddy-energy
         ! harmonic backscatter `ku` (CFL-floored) from the freshly-computed
         ! per-face resolved viscosity, so the hvisc Laplacian returns energy to
         ! the resolved flow.  Reads the prior thermo step's `ku` (one-step lag).
         ! No-op unless `meke%backscatter`; needs the face fields (LMIX active).
         if (present(meke) .and. present(lateral_mix)) then
            if (lateral_mix%is_init) then
               call meke_backscatter_apply(grid, metrics, meke, dt, &
                                           lateral_mix%ah_face_x, lateral_mix%ah_face_y)
            end if
         end if
         ! SPEC S3: under pred_corr the viscous tendency also reads the
         ! time-means (MOM6 `horizontal_viscosity(u_av, v_av, h_av, ...)`).
         ! The lateral-mix COEFFICIENT
         ! computation above still reads the prognostic (documented
         ! deviation: the closure coefficient is one flow-state stale,
         ! second-order in dt; full parity is a follow-on).
         if (dyn%split_scheme == SPLIT_SCHEME_PRED_CORR) then
            call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, &
                                                               lateral_mix=lateral_mix, dt=dt, &
                                                               u_src=ms%u_av_layer, &
                                                               v_src=ms%v_av_layer, &
                                                               h_src=ms%h_av_layer)
         else
            call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, &
                                                               lateral_mix=lateral_mix, dt=dt)
         end if
         ! KE dissipation rate for the MEKE frictional source (no-op unless on).
         call ocean_horizontal_viscosity_compute_ke_diss(hv, ms)
      end if
      call profiler_stop("ocean_hvisc")
      call profiler_start("ocean_bdrag")
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, dt)
      call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
      call profiler_stop("ocean_bdrag")
      call profiler_start("ocean_surfstress")
      call ocean_surface_stress_compute_tendencies(grid, ss, ms)
      call profiler_stop("ocean_surfstress")
      call probe_dS(grid, ms, "after slow computes", stage_id, step_id)

      ! Per-edge OBC tags + physical-edge flags — consumed by
      ! subtract_fast_cor_ref (below) and the post-substep transport
      ! wall-zeroing.  Defaults: WALL tags, .true. flags ⇒ single-rank
      ! closed-basin behaviour.
      bc_w_drv = OBC_WALL
      bc_e_drv = OBC_WALL
      bc_s_drv = OBC_WALL
      bc_n_drv = OBC_WALL
      has_w_drv = .true.
      has_e_drv = .true.
      has_s_drv = .true.
      has_n_drv = .true.
      if (present(bc)) then
         bc_w_drv = bc%west%bc_type
         bc_e_drv = bc%east%bc_type
         bc_s_drv = bc%south%bc_type
         bc_n_drv = bc%north%bc_type
         has_w_drv = bc%has_west
         has_e_drv = bc%has_east
         has_s_drv = bc%has_south
         has_n_drv = bc%has_north
      end if

      ! ---- 3. Sum slow tendencies into F_slow_u/v + depth-mean to F_bt ----
      ! Builds the constant-per-substep forcing for the barotropic substep.
      ! Uses h^n (which is still untouched) for the depth weighting.
      call profiler_start("ocean_F_slow_assembly")
      ! BT-budget probe: per-region per-term BT power decomposition.
      ! Pull device buffers down to host first (the probe is host-side
      ! Fortran loops + writes to stdout).  Gated on the namelist knob
      ! so production runs pay nothing.
      if (dyn%debug_bt_budget) then
         !$acc update self(ms%h_layer)
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         !$acc update self(dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt)
         !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
         !$acc update self(cor%pv_flux_x%data, cor%pv_flux_y%data)
         !$acc update self(hv%du_visc%data, hv%dv_visc%data)
         !$acc update self(bd%du_drag%data, bd%dv_drag%data)
         !$acc update self(ss%du_stress%data, ss%dv_stress%data)
         call print_bt_budget(grid, ms, dyn%bt_work, pgf, cor, hv, bd, ss, &
                              real(step_id, wp), "step", "S"//achar(48 + stage_id), &
                              header=(step_id == 1 .and. stage_id == 1))
      end if
      ! Pre-substep visc_rem refresh (MOM6 order parity): the rem-weighted
      ! forcing below, the γ-weighted continuity renormaliser, and the Δu
      ! corrector all read visc_rem_u/v; without this refresh they lag one
      ! stage — fatally at step 1 stage 1, where the init value (≡ 1)
      ! makes F_bt the plain mean of the ballistic spurious tendencies
      ! (PGF_BUG.md §9).  MOM6 computes vertvisc_coef + remnant in the
      ! predictor BEFORE btstep/continuity (SPEC §2 P5).
      if (dyn%bt_work%bt_forcing_visc_rem .or. dyn%bt_work%bt_renorm_visc_rem &
          .or. is_pc) then
         call visc_rem_precompute(grid, dyn, vmix, vd, ss, bd, ms, dt, kshear=kshear)
      end if
      call sum_slow_tendencies_into_F_slow(dyn%bt_work, pgf, cor, hv, bd, ss, ms)
      ! MOM6 wt_u parity (`&ocean_bt_nml forcing_visc_rem`): weight the
      ! forcing depth-mean by h·visc_rem so layers the implicit friction
      ! will immediately damp (grounded stacks under the vdiff BBL glue)
      ! do not force the fast loop (PGF_BUG.md §9).  The PGF-projection
      ! subtraction below MUST use the
      ! same weights — a plain-h-weighted subtraction would re-inject the
      ! very depth-mean the weighting removed, with the opposite sign.
      ! visc_rem lags one stage (the stage-end vdiff producer), same as
      ! the corrector's documented convention; it initializes to 1 so the
      ! first stage degenerates to the plain h-mean.
      if (dyn%bt_work%bt_forcing_visc_rem) then
         call face_depth_mean_rem_u(grid, dyn%bt_work%F_slow_u, ms%h_layer, &
                                    dyn%bt_work%visc_rem_u, dyn%bt_work%F_bt_u, ms%nz_ml)
         call face_depth_mean_rem_v(grid, dyn%bt_work%F_slow_v, ms%h_layer, &
                                    dyn%bt_work%visc_rem_v, dyn%bt_work%F_bt_v, ms%nz_ml)
      else
         call face_depth_mean_u(grid, dyn%bt_work%F_slow_u, ms%h_layer, dyn%bt_work%F_bt_u, ms%nz_ml)
         call face_depth_mean_v(grid, dyn%bt_work%F_slow_v, ms%h_layer, dyn%bt_work%F_bt_v, ms%nz_ml)
      end if
      ! Subtract the bt projection of the PGF: the barotropic substep has
      ! its own `-G·∂η/∂x` term, so without this subtraction the
      ! bt mode would integrate the PGF twice and √(gH) inflates to
      ! √(2gH).
      if (dyn%bt_work%bt_forcing_visc_rem) then
         call face_depth_mean_rem_u(grid, pgf%dpdx_face%data, ms%h_layer, &
                                    dyn%bt_work%visc_rem_u, dyn%bt_work%F_bt_u_fast, ms%nz_ml)
         call face_depth_mean_rem_v(grid, pgf%dpdy_face%data, ms%h_layer, &
                                    dyn%bt_work%visc_rem_v, dyn%bt_work%F_bt_v_fast, ms%nz_ml)
      else
         call face_depth_mean_u(grid, pgf%dpdx_face%data, ms%h_layer, dyn%bt_work%F_bt_u_fast, ms%nz_ml)
         call face_depth_mean_v(grid, pgf%dpdy_face%data, ms%h_layer, dyn%bt_work%F_bt_v_fast, ms%nz_ml)
      end if
      do concurrent(j=1:ny_uface, i=1:nx_face)
         dyn%bt_work%F_bt_u_fast(i, j) = dyn%bt_work%F_bt_u(i, j) - dyn%bt_work%F_bt_u_fast(i, j)
      end do
      do concurrent(j=1:ny_face, i=1:nx_vface)
         dyn%bt_work%F_bt_v_fast(i, j) = dyn%bt_work%F_bt_v(i, j) - dyn%bt_work%F_bt_v_fast(i, j)
      end do
      ! Subtract the fast-loop Coriolis + advection evaluated at the
      ! stage-entry bt state — the Coriolis/advection analogue of the PGF
      ! projection subtraction just above (MOM6 `Cor_ref_u/v`).  Without
      ! it the barotropic Coriolis is integrated twice (frozen in F_bt_u
      ! via `cor%pv_flux_*`'s depth mean AND live in the substep), which
      ! pumps an exponential wall-trapped barotropic mode on shelf rims
      ! (the 600² Lagrangian double-gyre h-guard trap).  `bt_ubt/bt_vbt`
      ! still hold the stage-entry state here (`derive_bt_from_layers`
      ! filled them; the substep's copy-in below re-sets them).  The
      ! per-edge tags/flags are derived early (above) for this call.
      !
      ! Wet/dry exclusion (v1 composition deferral, same class as the
      ! wetdry × use_cont_type/upstream_h_face exclusions): the reference
      ! is evaluated on the ungated stage-entry velocities, but the
      ! substep's live terms are zeroed at dry faces by `wd_open_u/v` —
      ! subtracting an ungated reference there injects net forcing at
      ! dry faces and drives `h_layer` negative
      ! (test_ocean_wetdry_driver).  Skipping keeps wetdry runs on the
      ! pre-fix forcing (the double-count is mild at coastal scales).
      if (.not. dyn%bt_work%wetdry_enable) then
         ! Reference velocity first.  Under `ssp_rk2` this is a copy of
         ! the stage-entry `bt_ubt/bt_vbt` the reference used to read
         ! directly (bit-identical); under `pred_corr` it is the depth
         ! mean of `u_av/v_av` — the velocity step 2 actually evaluated
         ! `cor%pv_flux_*` on — under the same weights the forcing
         ! depth-mean above used.  Without this the uncancelled
         ! `f × (v̄_av − v̄^n)` forces every substep and pumps the basin's
         ! gravest Poincaré seiche (see `set_cor_ref_velocity`).
         call set_cor_ref_velocity(grid, dyn%bt_work, ms, is_pc)
         call subtract_fast_cor_ref(grid, metrics, dyn%bt_work, cor%f_corner, &
                                    bc_w_drv, bc_e_drv, bc_s_drv, bc_n_drv, &
                                    has_w_drv, has_e_drv, has_s_drv, has_n_drv)
      end if
      call profiler_stop("ocean_F_slow_assembly")

      ! ---- 4. Run barotropic substep ----
      ! Outputs bt_eta_end, bt_ubt_end, bt_vbt_end (Hallberg end-step
      ! anchors for the Δu correction), bt_eta / bt_ubt / bt_vbt
      ! (time-mean for diagnostics), and — new — bt_uhbt / bt_vhbt
      ! (time-mean depth-integrated transport).  The slow continuity
      ! below renormalises its per-layer mass fluxes to vertically
      ! sum to those, so the layer h evolution is consistent with
      ! the barotropic-substep η evolution by construction (MOM6 split-RK2).
      do concurrent(j=1:ny_uface, i=1:nx_face)
         dyn%bt_work%bt_ubt(i, j) = dyn%bt_work%ubt_at_n(i, j)
      end do
      do concurrent(j=1:ny_face, i=1:nx_vface)
         dyn%bt_work%bt_vbt(i, j) = dyn%bt_work%vbt_at_n(i, j)
      end do
      ! bt_eta is fresh from derive_bt_from_layers — leave it.

      dt_inner = dt/real(n_inner, wp)
      ! MOM6 bt_rem_u: populate the per-face multiplicative damping
      ! factor read by the substep loop.  When the knob is off, leave
      ! `bt_rem_u/v` at their init value of 1 ⇒ multiplication is a
      ! no-op (bit-identical to pre-knob path).
      if (dyn%bt_work%bt_substep_drag) then
         call compute_bt_rem(grid, dyn%bt_work, ms, bd%r_linear, bd%hbbl, dt_inner)
      else if (dyn%bt_work%lwd_enable) then
         ! `bt_rem_u/v` is reset ONLY by `compute_bt_rem` above; when
         ! `substep_drag` is off but wave drag is on, nothing else resets
         ! it, and `compute_bt_rem_wave_drag` below MULTIPLIES into it —
         ! without this reset bt_rem would compound geometrically across
         ! outer steps (bt_rem = R^n after n stages), silently annihilating
         ! the barotropic mode. See `src/core/ocean/README.md`.
         call reset_bt_rem(grid, dyn%bt_work)
      end if
      if (dyn%bt_work%lwd_enable) then
         ! Barotropic linear wave drag (Egbert & Ray 2001; Jayne & St
         ! Laurent 2001): MULTIPLIES the static piston-velocity map into
         ! `bt_rem_u/v`, composing with `substep_drag` exactly as MOM6
         ! composes `lin_drag_u` with the viscous remnant.
         ! Must run AFTER the base fill
         ! above and BEFORE `mask_bt_rem` (land masking must be last).
         call compute_bt_rem_wave_drag(grid, dyn%bt_work, ms, dt_inner)
      end if
      ! Fold the static land face masks into bt_rem (C4 / R4a): zeroes the
      ! BT-substep velocity update across land faces.  No-op for all-wet.
      call mask_bt_rem(grid, metrics, dyn%bt_work)
      call profiler_start("ocean_barotropic_solver")
      ! `eta_forcing` is itself optional here: passing an absent optional
      ! as the actual for the substep's optional dummy propagates absence
      ! (F2018 15.5.2.13) ⇒ bit-identical when the C1 tide is off.
      if (dyn%bt_halo > 0) then
         ! ---- Wide-halo march-in path (Phase 3c) ----
         ! Tides × bt_halo is a configure-time exclusion.
         if (present(eta_forcing)) then
            error stop "run_stage_split: bt_halo > 0 with eta_forcing — excluded at configure time"
         end if
         call dyn%bt_wide%copy_in(grid, &
                                  dyn%bt_work%bt_eta, dyn%bt_work%bt_H_ref, &
                                  dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt, &
                                  dyn%bt_work%bt_ubt_prev, dyn%bt_work%bt_vbt_prev, &
                                  dyn%bt_work%bt_rem_u, dyn%bt_work%bt_rem_v, &
                                  dyn%bt_work%F_bt_u_fast, dyn%bt_work%F_bt_v_fast)
         call dyn%bt_wide%entry_exchange()
         call bt_wide_substep(dyn%bt_wide, dyn%bt_work, n_inner, dt_inner, bc)
         ! Copy wide outputs back to normal-width bt_work fields.
         ! The copy_out DC loops are synchronous (no acc async); drain async(1)
         ! first so copy_out reads the completed time-mean arrays.
         !$acc wait(1)
         call dyn%bt_wide%copy_out(grid, &
                                   dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt, &
                                   dyn%bt_work%bt_uhbt, dyn%bt_work%bt_vhbt, &
                                   dyn%bt_work%bt_eta_end, dyn%bt_work%bt_ubt_end, &
                                   dyn%bt_work%bt_vbt_end)
         ! Exit seam freshen: leaves downstream consumers exactly what v1's
         ! last in-loop exchange left (fresh normal-width seam ghosts).
         call profiler_start("ocean_comms_bt")
         call ocean_halo_bt_group_2d(dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, &
                                     dyn%bt_work%bt_vbt)
         call profiler_stop("ocean_comms_bt")
      else
         call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                    cor%f_corner, n_inner, dt_inner, &
                                                    bc=bc, t=t, eta_forcing=eta_forcing)
      end if
      call profiler_stop("ocean_barotropic_solver")
      call probe_dS(grid, ms, "after barotropic substep", stage_id, step_id)
      ! BT in/out chksums at the substep exit: when the fold's loud count
      ! fires, these rows name which fold INPUT (ubt_end/ubt_at_n/F_bt)
      ! went non-finite — the producer the nan-catch cannot see.
      call chksum_bt(grid, dyn%bt_work, dyn%chksum_probe, "post_bt", stage_id, step_id)

      ! Physical-wall reconciliation for the transport constraint.
      ! The barotropic substep closes `bt_ubt` at array-edge faces (i=1,
      ! i=nx+1) — convenient for its own integration but several
      ! cells away from where the slow continuity closes (the
      ! physical walls at i = nghost+1 and i = nghost+nx_phys+1).
      ! Passing the raw barotropic-substep transport into continuity would
      ! ask for nonzero mass flux through the physical wall, which
      ! the wall-zero step then erases — leaving the constraint
      ! unsatisfied at those faces and a small leak in
      ! `sum_k(h_layer)` vs `H_ref + bt_eta_end`.  Zero the wall
      ! faces of `bt_uhbt / bt_vhbt` so the constraint is
      ! self-consistent with the closed-wall slow continuity.
      !
      ! OBC dispatch: for non-WALL edges, the barotropic substep already
      ! computed a physical transport at the wall face (Flather,
      ! Chapman, clamped data) — DON'T zero it.  The slow continuity
      ! will honour the same `bc` tag and let the matching mass
      ! flux through.  (Tags/flags bc_*_drv / has_*_drv are derived
      ! earlier, before the F_slow assembly, for subtract_fast_cor_ref.)
      ! Seam faces carry real transport; the halo owns them (O0, plan D4).
      do concurrent(j=1:ny_uface)
         if (bc_w_drv == OBC_WALL .and. has_w_drv) dyn%bt_work%bt_uhbt(grid%nghost + 1, j) = 0.0_wp
         if (bc_e_drv == OBC_WALL .and. has_e_drv) dyn%bt_work%bt_uhbt(grid%nghost + grid%nx_phys + 1, j) = 0.0_wp
      end do
      do concurrent(i=1:nx_vface)
         if (bc_s_drv == OBC_WALL .and. has_s_drv) dyn%bt_work%bt_vhbt(i, grid%nghost + 1) = 0.0_wp
         if (bc_n_drv == OBC_WALL .and. has_n_drv) dyn%bt_work%bt_vhbt(i, grid%nghost + grid%ny_phys + 1) = 0.0_wp
      end do

      ! Continuity + tracer chain.  ssp_rk2: here (historical position,
      ! bit-identical).  pred_corr: deferred to AFTER the velocity update +
      ! implicit friction — the forward-backward pairing (SPEC §2 C8).
      if (.not. is_pc) then
         call run_continuity_chain(grid, metrics, dyn, ct, hd, va, redi, varmix, ms, &
                                   dt, therm_dt, therm_active, is_lagrangian, &
                                   h_min_floor, chain_weight, is_pred, &
                                   stage_id, step_id, bc=bc, mle=mle, gm=gm)
      end if

      ! ---- 6. Velocity-tendency applies ----
      ! Async-chained on OpenACC queue 1: these five applies are all
      ! additive forward-Euler accumulations onto u_face/v_face with no
      ! intervening default-queue or host-reading op between them (the
      ! profiler calls are host-side wall-clock timers + NVTX markers; they
      ! neither sync the device nor read device data).  Same queue ⇒ FIFO
      ! ⇒ the additive sequence is order-preserved.  ONE `!$acc wait(1)`
      ! before `apply_bt_correction` — the first device consumer that reads
      ! the freshly-applied u_face/v_face — completes the chain.  This
      ! eliminates the per-launch host sync gap (~half the per-stage apply
      ! cost in nsys).
      ! accel_visc_rem (MOM6 parity): snapshot the pre-apply velocity so
      ! the post-apply reweight below can attenuate the whole explicit
      ! tendency sum by the per-layer viscous remnant.  The `!$acc wait`
      ! orders the snapshot against any still-in-flight async producer of
      ! u_face/v_face (the q1 chain convention).
      if (dyn%accel_visc_rem) then
         !$acc wait
         call accel_visc_rem_snapshot(grid%nx_total + 1, grid%ny_total, ms%nz_ml, &
                                      ms%u_face_x_layer, dyn%avr_u0)
         call accel_visc_rem_snapshot(grid%nx_total, grid%ny_total + 1, ms%nz_ml, &
                                      ms%v_face_y_layer, dyn%avr_v0)
      end if
      call profiler_start("ocean_velocity_apply")
      ! `dt_vel` = dt (historical / corrector) or BE·dt (pred_corr
      ! predictor, SPEC §2 P8).
      ! chksum seams (&ocean_debug_nml chksum, inert when off): each sample
      ! waits the device, so the async queue-1 chain is serialised ONLY inside
      ! the chksum window — the per-tendency attribution that named the
      ! 30 m/s injector.  FIFO order is unaffected (same queue).
      call coriolis_adv_apply_tendencies(cor, ms, dt_vel, no_wait=.true.)
      call chksum_state(grid, ms, dyn%chksum_probe, "post_coradv", stage_id, step_id)
      call chksum_hotface(grid, ms, dyn%bt_work%visc_rem_u, dyn%bt_work%visc_rem_v, &
                          dyn%chksum_probe, "post_coradv", stage_id, step_id)
      call ocean_pressure_force_apply(pgf, ms, dt_vel, no_wait=.true.)
      call chksum_state(grid, ms, dyn%chksum_probe, "post_pgf", stage_id, step_id)
      call chksum_hotface(grid, ms, dyn%bt_work%visc_rem_u, dyn%bt_work%visc_rem_v, &
                          dyn%chksum_probe, "post_pgf", stage_id, step_id)
      call ocean_horizontal_viscosity_apply_tendencies(hv, ms, dt_vel, no_wait=.true.)
      call chksum_state(grid, ms, dyn%chksum_probe, "post_hvisc", stage_id, step_id)
      ! Double-count guard (see run_stage): skip the explicit per-layer
      ! drag/stress apply when it is folded into the implicit vdiff
      ! tridiagonal.  NOTE (split path): the drag/stress tendency buffers
      ! still feed F_slow → the barotropic mode above; the implicit fold
      ! here acts on the post-bt-correction per-layer (baroclinic)
      ! velocity.  PR-19 supplies the visc_rem QUANTITY (produced by
      ! vdiff_apply_momentum, consumed by apply_bt_correction's
      ! h-weighted branch below); feeding it into F_slow / the
      ! barotropic substep itself — the actual barotropic-coupling flip
      ! — is PR-56's territory (changes the barotropic mode, gated on
      ! the Bleck/Hallberg instability test).  For strict split-path use
      ! today, keep the explicit / &ocean_bdrag_nml implicit split-apply.
      ! Defaults (folds off) ⇒ both applies run ⇒ bit-identical.
      if (.not. vd%implicit_drag) then
         call ocean_bottom_drag_apply_tendencies(bd, ms, dt_vel, no_wait=.true.)
      end if
      call ocean_channel_drag_apply_tendencies(bd, ms, dt_vel, no_wait=.true.)
      if (.not. vd%implicit_stress) then
         call ocean_surface_stress_apply_tendencies(ss, ms, dt_vel, no_wait=.true.)
      end if
      call chksum_state(grid, ms, dyn%chksum_probe, "post_drag", stage_id, step_id)
      call profiler_stop("ocean_velocity_apply")
      call probe_dS(grid, ms, "after velocity applies", stage_id, step_id)

      ! ---- 7. bt correction ----
      ! In Eulerian-z mode: Δu correction + h-rescale (the rescale
      ! brings any small FP drift in `sum_k(h_layer)` back to
      ! `H_ref + bt_eta_end`, important for gravity-wave temporal
      ! continuity).  In Lagrangian mode: only the Δu correction —
      ! the slow continuity's `h_layer` is authoritative, and any
      ! FP residual gets absorbed by the ALE remap at end of outer
      ! step rather than redistributed across layers (which would
      ! corrupt S = hTr/h).
      ! bc-PGF correction: form `e_anom` from the post-substep
      ! η state and the pre-PGF snapshot.  Then apply_bt_correction
      ! consumes it via the optional `use_bc_pgf` path.  Gated;
      ! no-op when bt_correction_bc_pgf is off.
      if (dyn%bt_work%bt_correction_bc_pgf) then
         call compute_e_anom(dyn%bt_work)
      end if
      ! Close the async velocity-apply chain (queue 1) before the first
      ! device consumer reads the applied u_face/v_face.  apply_bt_correction
      ! reads u_face_x_layer / v_face_y_layer on the default queue, so the
      ! whole coriolis→pgf→hvisc→bdrag→surfstress chain must have landed.
      !$acc wait(1)
      ! accel_visc_rem reweight: `u = u_entry + visc_rem·(u − u_entry)` —
      ! friction-dominated near-massless layers cannot keep a
      ! full-strength explicit dt·F kick (MOM6 MOM_dynamics_split_RK2:
      ! `u = u_init + dt·visc_rem_u·(CAu + PFu + diffu)`; the 2026-07-28
      ! forensics injector-#2 fix).  visc_rem is the LAGGED stage
      ! producer (vdiff_apply_momentum) — ≡ 1 before the first vdiff ⇒
      ! no-op, matching MOM6's semantics.
      if (dyn%accel_visc_rem) then
         call accel_visc_rem_reweight(grid%nx_total + 1, grid%ny_total, ms%nz_ml, &
                                      dyn%avr_u0, dyn%bt_work%visc_rem_u, ms%u_face_x_layer)
         call accel_visc_rem_reweight(grid%nx_total, grid%ny_total + 1, ms%nz_ml, &
                                      dyn%avr_v0, dyn%bt_work%visc_rem_v, ms%v_face_y_layer)
      end if
      call profiler_start("ocean_bt_correction")
      call apply_bt_correction(dyn%bt_work, ms, dt, &
                               skip_h_rescale=is_lagrangian, &
                               use_h_weighted=dyn%bt_work%bt_correction_h_weighted, &
                               grid=grid, &
                               use_bc_pgf=dyn%bt_work%bt_correction_bc_pgf, &
                               use_visc_rem=dyn%bt_work%bt_correction_visc_rem, &
                               metrics=metrics, &
                               scale=merge(dyn%pc_be, 1.0_wp, is_pred), &
                               n_nonfin=dyn%bt_nonfin_step)
      if (dyn%bt_nonfin_step > 0) then
         write (output_unit, '("[nan-catch] stage ", i0, " step ", i0, ": ", i0, &
            &" non-finite BT-correction faces (fold skipped — BT loop blown up?)")') &
            stage_id, step_id, dyn%bt_nonfin_step
         flush (output_unit)
      end if
      call chksum_state(grid, ms, dyn%chksum_probe, "post_bt_fold", stage_id, step_id)
      ! Land-face velocity reset (C4 / R4b.2): zero re-ingested land-face
      ! layer velocity after the BT correction.  No-op for all-wet.
      call mask_layer_velocities(grid, metrics, ms, bt_work=dyn%bt_work)
      ! Phase-2 vanished-layer velocity reset (R5): runs AFTER apply_bt_correction
      ! so it sees the recombined per-layer velocity.  Only for VCOORD_LAGRANGIAN.
      if (dyn%reset_vanished_u .and. present(vcoord)) then
         if (vcoord%coord_type == VCOORD_LAGRANGIAN) then
            call reset_vanished_layer_velocities(ms, &
                                                 isopycnal_vanish_tol(dyn%angstrom_h, pd_floor=ct%positive_definite))
         end if
      end if
      call profiler_stop("ocean_bt_correction")
      call probe_dS(grid, ms, "after apply_bt_correction", stage_id, step_id)
      call probe_h_vs_eta_residual(grid, ms, dyn%bt_work, stage_id, step_id)

      ! Baroclinic OBC: set per-layer normal velocity at open faces to
      ! Flather mean + zero-gradient baroclinic anomaly (OPEN/TIDAL/
      ! CHAPMAN/NESTED) or uniform clamped_u/v (CLAMPED).  Runs AFTER
      ! apply_bt_correction so it sees the recombined per-layer velocity.
      ! No-op when bc is absent or no edge is open-ish.
      if (present(bc)) call ocean_obc_apply_baroclinic(grid, bc, dyn%bt_work, ms, dt)

      ! Sponge.  Map-driven path (&ocean_sponge_nml enable, PR-23)
      ! supersedes the legacy per-edge band; exactly one of the two runs.
      ! Both run AFTER bt-correction so they see the recombined per-layer
      ! velocity. `sponge_maps_on` is a local logical because Fortran does
      ! not guarantee `.and.` short-circuits past `present()`.
      sponge_maps_on = .false.
      if (present(sp)) sponge_maps_on = sp%enable
      if (sponge_maps_on) then
         call ocean_sponge_apply_maps(grid, sp, ms, dt)
      else if (present(bc)) then
         call ocean_sponge_apply(grid, bc, ms, dt)
         call ocean_sponge_apply_tracers(grid, bc, ms, dt)
      end if

      ! ---- 8. Surface tracer fluxes (heat / salt) ----
      ! Same ordering as the unsplit driver: after horizontal +
      ! vertical tracer transport, before vmix/vdiff.  Both `sf`
      ! and `active` are optional — kernel no-ops when either is
      ! absent / false.
      ! Wet/dry: the dynamic cell mask gates the flux (a dry column's
      ! mm-scale sliver must not be heated, plan §4.4).  sw_pen /
      ! restore / geothermal are configure-time excluded with wetdry
      ! (their additive/piston structure mis-composes with a masked
      ! main deposit); knob off ⇒ the original call, byte-identical.
      if (dyn%bt_work%wetdry_enable) then
         call ocean_surface_flux_apply_tracers(grid, sf, ms, therm_dt, &
                                               active=therm_active, &
                                               wet_dyn=dyn%bt_work%wd_wet_dyn)
      else
         call ocean_surface_flux_apply_tracers(grid, sf, ms, therm_dt, active=therm_active)
      end if
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, therm_dt, active=therm_active)
      call ocean_surface_restore_apply_tracers(grid, sf, ms, therm_dt, active=therm_active)
      call ocean_geothermal_apply_tracers(grid, geo, ms, therm_dt, active=therm_active)
      call probe_dS(grid, ms, "after surface flux", stage_id, step_id)
      ! Ideal-age tracer: interior aging only (1 s/s), thermo-cadence
      ! gated (PR-7).  The surface Dirichlet reset runs once per outer
      ! step, after rk2_average + the ALE remap, in
      ! `ocean_dyn_step_split` — see the ordering comment there.
      !
      ! `.not. is_pred` is what makes the source rate SCHEME-INDEPENDENT.
      ! Under ssp_rk2 this fires on BOTH identical stages and `rk2_average`
      ! turns the two `+therm_dt` into exactly one — the contract
      ! `rdb_ocean_ideal_age`'s header documents.  Under pred_corr there is no
      ! average, so a predictor application would simply survive into the
      ! corrector's and age the ocean at 2x real time (caught by
      ! `split_ideal_age_grows_linearly_via_driver`, which asserts
      ! age == N*DT).  Skipping the predictor leaves ONE application per
      ! outer step, which is the same +therm_dt.  ssp_rk2 is unaffected:
      ! `is_pred` is false on both of its stages.
      call ocean_ideal_age_apply(grid, ms, therm_dt, &
                                 active=dyn%is_thermo_step() .and. .not. is_pred)

      ! ---- 9. Vertical mixing (implicit, stable under any dt) ----
      ! `bt_work=dyn%bt_work` is the visc_rem PRODUCER call: when
      ! `bt_correction_visc_rem` is on, this (re)fills `bt_work%visc_rem_u/v`
      ! from THIS stage's momentum solve.  The BT correction that CONSUMES
      ! it already ran at step 7 above, so the corrector always reads the
      ! PREVIOUS stage's γ — a one-stage (Δt/2) lag, accepted for v1 (see
      ! `vmix_apply_in_stage`'s `bt_work` docstring and PLAN_PR19 §11.1).
      ! At stage 1 of step 1, γ is still at its `source=1.0` init, so the
      ! very first correction is h-only ⇒ benign.
      ! pred_corr predictor: MOM6 applies the implicit friction to the
      ! provisional velocity too — `vertvisc(up, dt_pred)` (SPEC §1 fact 8's
      ! coef-only claim is WRONG) — so the
      ! spurious grounded-layer ballistics are absorbed BEFORE the
      ! predictor continuity forms u_av.  Momentum-only there (tracers
      ! untouched); dt_vel = BE·dt matches MOM6's dt_pred.
      if (is_pred) then
         call vmix_apply_in_stage(grid, dyn, vmix, vd, ss, bd, ms, dt_vel, stage, sf, epbl=epbl, &
                                  kshear=kshear, vmix_tidal=vmix_tidal, bt_work=dyn%bt_work, &
                                  apply_tracers=.false., metrics=metrics)
      else
         call vmix_apply_in_stage(grid, dyn, vmix, vd, ss, bd, ms, dt, stage, sf, epbl=epbl, &
                                  kshear=kshear, vmix_tidal=vmix_tidal, bt_work=dyn%bt_work, &
                                  metrics=metrics)
      end if
      ! KE attribution: implicit vertical friction (+ folded drag/stress
      ! when implicit_*) — the stage-close segment (debug_ke_attr).
      call ke_probe_sample(grid, ms, dyn%ke_probe, "vdiff_vmix", stage_id, step_id)
      call chksum_state(grid, ms, dyn%chksum_probe, "stage_close", stage_id, step_id)
      call probe_dS(grid, ms, "after vdiff/vmix (stage end)", stage_id, step_id)

      ! pred_corr: the continuity + tracer chain runs HERE — after the
      ! single prognostic update (corrector) / provisional up (predictor)
      ! and after the implicit friction — so the thickness advances with
      ! the UPDATED velocities via uhbt + u_cor.  This is the
      ! forward-backward gravity-wave pairing that lifts the
      ! internal-wave dt ceiling (SPEC §1 fact 5, §2 C8).
      if (is_pc) then
         call run_continuity_chain(grid, metrics, dyn, ct, hd, va, redi, varmix, ms, &
                                   dt, therm_dt, therm_active, is_lagrangian, &
                                   h_min_floor, chain_weight, is_pred, &
                                   stage_id, step_id, bc=bc, mle=mle, gm=gm)
      end if

      ! ---- 10. Stage-end seam reconciliation (tripolar only) ----
      ! The velocity applies + bt-correction + sponge updated v/u/h/tracers
      ! AFTER the post-continuity fold; the duplicated-DOF v/corner row must
      ! be re-projected and the north ghosts re-folded so the stage output
      ! is fully seam-consistent (the v on-row antisymmetric projection is
      ! "applied after each update of v" — Appendix A).  Periodic-first-
      ! fold-second: re-wrap periodic ghosts, then fold.  No-op when not
      ! folding (and when no periodic axis is set the periodic wrap no-ops too).
      ! O2: unconditional multi-rank ghost exchange (D0 — no-op on 1 rank).
      call ocean_halo_exchange_ml_state(ms)
      if (present(bc)) then
         if (bc%north_fold) then
            call ocean_periodic_wrap_state(grid, bc, ms, &
                                           skip_x=ocean_halo_is_decomposed_x(), skip_y=ocean_halo_is_decomposed_y())
            call ocean_fold_wrap_state(grid, bc, ms)
         end if
      end if
   end subroutine run_stage_split

   pure function ocean_dyn_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the split-RK2 driver (BT work state
      !! + wide-halo shadow state) slot (0 when unallocated).
      class(ocean_dyn_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%bt_work%bytes()
      ! accel_visc_rem stage-entry snapshots (allocated at setup when the
      ! knob is on; 0 on the default-off path where they stay unallocated).
      nbytes = nbytes + arr_bytes(this%avr_u0) + arr_bytes(this%avr_v0)
      ! Wide-halo BT march-in shadow state (`&ocean_bt_nml bt_halo > 0`).
      ! Allocated + device-mapped by `ocean_dyn_enable_bt_wide`, which the
      ! driver runs AFTER `ocean_state_enter_data`; the driver therefore
      ! latches `mem_set_counted_budget` after that call so this term is in
      ! the snapshot.  Unallocated at bt_halo = 0 ⇒ 0.
      if (allocated(this%bt_wide)) nbytes = nbytes + this%bt_wide%bytes()
   end function ocean_dyn_bytes

end module rdb_ocean_dyn
