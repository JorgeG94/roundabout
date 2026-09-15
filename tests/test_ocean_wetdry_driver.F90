!! Full-driver + surface-flux-gate coverage for ocean dynamic
!! wetting/drying (docs/ocean_wetdry_plan.md §4.2, §4.4, §8).
!!
!! The sibling `test_ocean_wetdry` drives `barotropic_substep_nonlinear`
!! in ISOLATION — no tracers, no layers, no ALE remap / vdiff / clamp.
!! That leaves the plan's most consequential downstream claims
!! (§4.2 "hTr → 0 conservatively", §4.4 "surface fluxes masked on dry
!! columns") with NO test coverage.  This module closes that gap:
!!
!!   1. `wetdry_tracer_conserv_cycle` — drives the FULL split-explicit
!!      driver (`ocean_dyn_step_split`) on a flat shallow basin
!!      (VCOORD_SIGMA, S/T tracers, ppm_limit_pos, closed walls) whose
!!      rest depth sits in the wet/dry band.  A gravest-mode cosine
!!      seiche sweeps each column's depth D across dry_depth and
!!      rewet_depth every half period, so `wd_wet_dyn` toggles
!!      1 -> 0 -> 1 (the operational dry->wet->dry cycle) while D stays
!!      strictly positive.  Asserts global tracer CONTENT (`Σ hS`,
!!      `Σ hT`) is conserved to round-off, no NaN/Inf anywhere, and
!!      positivity (`h_layer ≥ 0`, `D ≥ 0`) throughout — the
!!      coastal-ZSTAR_FULL salt-leak failure mode reachable on sigma.
!!      NOTE: the stricter collapse to D -> 0 / below-H_VANISHED is NOT
!!      stably reachable through this driver (it NaNs in `h_layer`);
!!      see the review report.
!!
!!   2. `wetdry_sflux_gate` — the surface-flux gate (§4.4, plan §8
!!      `test_ocean_wetdry_sflux`).  A dynamically-dry column
!!      (`wd_wet_dyn = 0`) must receive NO surface heat/salt flux (its
!!      top-layer T/S untouched).  Separately CHARACTERISES a held-wet
!!      thin-band column (`dry_depth < D < rewet_depth`): the binary
!!      `wet_dyn` gate applies the FULL flux into a mm-scale top layer
!!      with NO thickness-aware throttle, so its ΔT dwarfs a deep
!!      column's — the review's S1 concern.  The magnitudes are
!!      asserted against the analytic `Q·dt/(ρ·cp·h_top)` and reported.
module test_ocean_wetdry_driver
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, VCOORD_SIGMA
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_LITE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, SEAWATER_CP
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private

   public :: collect_ocean_wetdry_driver_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_wetdry_driver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wetdry_tracer_conserv_cycle", test_tracer_conserv_cycle), &
                  new_unittest("wetdry_sflux_gate", test_sflux_gate) &
                  ]
   end subroutine collect_ocean_wetdry_driver_tests

   ! -----------------------------------------------------------------
   ! TEST 1 — tracer conservation through a dry->wet->dry cycle,
   !          full split driver.
   ! -----------------------------------------------------------------

   subroutine test_tracer_conserv_cycle(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc

      integer, parameter :: NXP = 48, NYP = 4, NZ = 4
      integer, parameter :: N_STEPS = 200
      integer, parameter :: N_INNER = 15
      integer, parameter :: PROBE = 1
      real(wp), parameter :: DX = 200.0_wp
      real(wp), parameter :: H0 = 0.09_wp
         !! Flat shallow basin whose rest depth sits in the wet/dry
         !! hysteresis band region.  A gravest-mode cosine seiche
         !! (amplitude AMP) then sweeps each column's total depth D
         !! across dry_depth (dries, wd_wet_dyn -> 0) and rewet_depth
         !! (rewets, wd_wet_dyn -> 1) every half period — the operational
         !! dry->wet->dry cycle — while D stays strictly positive
         !! (AMP < H0 - a margin), so the full tracer pipeline runs on
         !! genuinely thin (but finite) columns.  See the report note on
         !! the D -> 0 / below-H_VANISHED extreme, which is NOT reachable
         !! stably through this driver (a separate finding).
      real(wp), parameter :: AMP = 0.045_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DRY_D = 0.05_wp, REWET_D = 0.10_wp
      real(wp), parameter :: CONSERV_TOL = 1.0e-12_wp
         !! Round-off gate on Σ hS / Σ hT relative drift.  Continuity-PPM
         !! + sigma ALE remap + vdiff are all flux-form/conservative, so
         !! the expectation is FP round-off; a drift above this measures a
         !! real leak in the tracer-under-drying path (the coastal
         !! ZSTAR_FULL failure mode on sigma).

      integer :: i, j, k, ig, i0, i1, j0, j1, step
      real(wp) :: frac, d_ic
      real(wp) :: hS0, hT0, hS1, hT1, drift_S, drift_T
      real(wp) :: min_h, min_d, worst_S, worst_T
      logical :: hit_nan, saw_dry, saw_rewet
      logical, allocatable :: was_wet_rest(:), dried(:)
      character(len=256) :: msg
      character(len=32) :: msg_field

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 0.0_wp
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call sf%init(grid)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      pgf%variant = OPGF_VARIANT_FV_LITE

      ! Quiescent physics: no wind, no surface flux, no bottom drag,
      ! no horizontal viscosity/diffusion.  Vertical tracer diffusion
      ! ON (small K_v) so the vdiff tridiagonal actually runs on the
      ! vanishing-then-rewetting columns (the vdiff mass-drop bug
      ! pattern, plan §4.2) — it must stay conservative.
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      sf%Q_heat_const = 0.0_wp
      sf%Q_salt_const = 0.0_wp
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 1.0e-4_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      ! wet/dry requires the layer positivity limiter.
      ct%use_ppm_limit_pos = .true.
      ! Wet/dry composition (mirrors rdb_ocean_setup wiring from
      ! cfg%ocean%wetdry%enable): keep the legacy single-step uhbt
      ! renormalisation — the Newton donor re-pick + CFL bracket
      ! destabilise the drying front (h_layer goes negative).
      ct%renorm_legacy_single_step = .true.

      allocate (was_wet_rest(grid%nx_total), source=.false.)
      allocate (dried(grid%nx_total), source=.false.)

      ! Flat shallow basin (rest depth H0, in the wet/dry band region)
      ! + gravest-mode cosine seiche.  Each column's D = H0 + eta sweeps
      ! across dry_depth and rewet_depth every half period, toggling
      ! wd_wet_dyn, while D stays > 0.  A cross-basin salinity gradient
      ! (uniform in z) makes the horizontal tracer flux non-trivial;
      ! sigma keeps the ALE tracer remap a vertical no-op so any Σ hS
      ! drift is a genuine continuity/vdiff/clamp leak.  All columns are
      ! wet at rest (H0 > rewet_depth).
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
            dyn%bt_work%bt_H_ref(i, j) = H0
            d_ic = max(H0 + AMP*cos(PI*frac), 0.0_wp)
            do k = 1, NZ
               ms%h_layer(i, j, k) = d_ic/real(NZ, wp)
               ! S: west saltier (S_ref+1) -> east fresher (S_ref-1);
               ! T uniform.  hTr = conc * h_layer.
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  (eos%S_ref + 2.0_wp*(0.5_wp - frac))*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  eos%T_ref*ms%h_layer(i, j, k)
               ms%rho_layer(i, j, k) = eos%rho0
            end do
            was_wet_rest(i) = .true.
         end do
      end do

      ! Enable wet/dry on the BT workstate the way configure_ocean_wetdry
      ! does: knobs + wd_* workspaces + hysteresis seed from D (eta=0
      ! reference here uses bt_H_ref).  Allocated BEFORE enter_data so
      ! dyn%enter_data attaches them to the device.
      call enable_wetdry_full(dyn, grid, DRY_D, REWET_D)

      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP

      hS0 = tracer_content(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)
      hT0 = tracer_content(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)

      saw_dry = .false.
      saw_rewet = .false.
      hit_nan = .false.
      min_h = huge(1.0_wp)
      min_d = huge(1.0_wp)

      call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call vc%enter_data()

      checks: block
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT_of(H0 + AMP, DX, N_INNER), N_INNER, &
                                      sf=sf, vcoord=vc)

            if (mod(step, PROBE) == 0 .or. step == N_STEPS) then
               !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
               !$acc&            ms%rho_layer, &
               !$acc&            ms%tracers(ms%idx_salinity)%hTr, &
               !$acc&            ms%tracers(ms%idx_temperature)%hTr, &
               !$acc&            dyn%bt_work%wd_wet_dyn)

               ! No NaN/Inf anywhere.
               if (.not. all_finite(ms, NZ, msg_field)) then
                  hit_nan = .true.
                  min_h = minval(ms%h_layer(i0:i1, j0:j1, :))
                  write (msg, '("wetdry cycle: NaN/Inf in ", a, " at step ", i0, &
                                &" (min_h=", es10.3, ")")') trim(msg_field), step, min_h
                  call check(error, .false., trim(msg)); exit checks
               end if

               ! Positivity: h_layer >= 0 and column depth D >= 0.
               min_h = min(min_h, minval(ms%h_layer(i0:i1, j0:j1, :)))
               do j = j0, j1
                  do i = i0, i1
                     min_d = min(min_d, column_depth(ms, i, j, NZ))
                  end do
               end do

               ! Hysteresis engaged: an interior column that was wet at
               ! rest dynamically dried (wd_wet_dyn -> 0, i.e. D fell
               ! below dry_depth) and later rewetted (wd_wet_dyn -> 1).
               ! This is the operational dry->wet->dry cycle the mask
               ! implements; the stricter below-H_VANISHED collapse is
               ! NOT stably reachable through this driver (report note).
               do i = i0, i1
                  if (dyn%bt_work%wd_wet_dyn(i, j0) < 0.5_wp .and. was_wet_rest(i)) then
                     saw_dry = .true.
                     dried(i) = .true.
                  else if (dyn%bt_work%wd_wet_dyn(i, j0) > 0.5_wp .and. dried(i)) then
                     saw_rewet = .true.
                  end if
               end do
            end if
         end do

         hS1 = tracer_content(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)
         hT1 = tracer_content(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)
         drift_S = abs(hS1 - hS0)/max(abs(hS0), tiny(1.0_wp))
         drift_T = abs(hT1 - hT0)/max(abs(hT0), tiny(1.0_wp))
         worst_S = drift_S
         worst_T = drift_T

         ! Machinery must have demonstrably engaged, else the test is
         ! blind to the drying path.
         write (msg, '("wetdry cycle: no interior column dried (min_D reached=", es10.3, ")")') min_d
         call check(error, saw_dry, trim(msg))
         if (allocated(error)) exit checks
         call check(error, saw_rewet, &
                    "wetdry cycle: no dried column ever rewetted (hysteresis stuck)")
         if (allocated(error)) exit checks

         call check(error, min_h >= 0.0_wp, "wetdry cycle: h_layer went negative")
         if (allocated(error)) exit checks
         call check(error, min_d >= -1.0e-12_wp, "wetdry cycle: column depth D went negative")
         if (allocated(error)) exit checks

         write (msg, '("wetdry cycle: Sigma hS drift_rel=", es12.5, " > ", es9.2, &
                       &" (min_h=", es10.3, ", min_D=", es10.3, ")")') &
            worst_S, CONSERV_TOL, min_h, min_d
         call check(error, drift_S < CONSERV_TOL, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("wetdry cycle: Sigma hT drift_rel=", es12.5, " > ", es9.2)') &
            worst_T, CONSERV_TOL
         call check(error, drift_T < CONSERV_TOL, trim(msg))
      end block checks

      call vc%exit_data()
      call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

      deallocate (was_wet_rest, dried)
      call vc%destroy()
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_tracer_conserv_cycle

   ! -----------------------------------------------------------------
   ! TEST 2 — surface-flux gate on thin / dry columns (plan §8
   !          `test_ocean_wetdry_sflux`, §4.4).
   ! -----------------------------------------------------------------

   subroutine test_sflux_gate(error)
      !! Three columns, one surface-flux apply:
      !!   * DRY   (wd_wet_dyn = 0): must be untouched (the gate).
      !!   * DEEP  (wd_wet_dyn = 1, D = 100 m): normal, tiny ΔT.
      !!   * THIN  (wd_wet_dyn = 1, held-wet band D in (dry,rewet)):
      !!     the binary gate applies the FULL flux into a mm-scale top
      !!     layer -> ΔT dwarfs DEEP.  Asserted against the analytic
      !!     `Q·dt/(ρ·cp·h_top)` and reported (review S1 gate gap).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      real(wp), allocatable :: wet_dyn(:, :)
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      integer, parameter :: NXP = 6, NYP = 4, NZ = 4
      real(wp), parameter :: Q_HEAT = 100.0_wp, Q_SALT = 1.0e-4_wp
      real(wp), parameter :: DT = 1200.0_wp
      real(wp), parameter :: D_DEEP = 100.0_wp
      real(wp), parameter :: D_THIN = 0.07_wp   ! held-wet band: dry(.05) < D < rewet(.10)
      integer :: nx, ny, ig, i, j, k, idT, idS
      integer :: i_dry, i_deep, i_thin
      real(wp) :: htop_deep, htop_thin
      real(wp) :: dT_dry, dS_dry, dT_deep, dT_thin, dS_thin
      real(wp) :: dT_thin_ana, dT_deep_ana, inv_heat, inv_salt
      character(len=320) :: msg

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      call sf%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      idT = ms%idx_temperature
      idS = ms%idx_salinity

      ! Three probe columns in the physical interior.
      i_dry = ig + 1
      i_deep = ig + 2
      i_thin = ig + 3
      htop_deep = D_DEEP/real(NZ, wp)
      htop_thin = D_THIN/real(NZ, wp)

      allocate (wet_dyn(nx, ny), source=1.0_wp)
      ! Column depths per sigma partition; wd_wet_dyn per column.
      do j = 1, ny
         do i = 1, nx
            do k = 1, NZ
               ms%h_layer(i, j, k) = htop_deep   ! default: deep everywhere
               ms%tracers(idT)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               ms%tracers(idS)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
            end do
         end do
      end do
      ! Override the DRY and THIN columns (all physical rows share the
      ! per-i configuration for a clean per-column average).
      do j = 1, ny
         call set_column(ms, i_dry, j, NZ, htop_deep, eos%T_ref, eos%S_ref)  ! depth irrelevant (masked)
         call set_column(ms, i_thin, j, NZ, htop_thin, eos%T_ref, eos%S_ref)
         wet_dyn(i_dry, j) = 0.0_wp   ! dynamically dry -> gate must block
      end do

      allocate (hT_ic, source=ms%tracers(idT)%hTr)
      allocate (hS_ic, source=ms%tracers(idS)%hTr)

      call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

      call run_apply_dyn(grid, sf, ms, wet_dyn, DT)

      ! Per-column top-layer (k = nz) concentration change, averaged
      ! over physical rows.
      dT_dry = col_dtop(ms, hT_ic, idT, i_dry, ig, NYP, NZ)
      dS_dry = col_dtop(ms, hS_ic, idS, i_dry, ig, NYP, NZ)
      dT_deep = col_dtop(ms, hT_ic, idT, i_deep, ig, NYP, NZ)
      dT_thin = col_dtop(ms, hT_ic, idT, i_thin, ig, NYP, NZ)
      dS_thin = col_dtop(ms, hS_ic, idS, i_thin, ig, NYP, NZ)

      ! Analytic ΔT = Q·dt/(ρ·cp·h_top) for a wet column.
      inv_heat = DT/(sf%rho0*sf%cp)
      inv_salt = DT/sf%rho0
      dT_deep_ana = Q_HEAT*inv_heat/htop_deep
      dT_thin_ana = Q_HEAT*inv_heat/htop_thin

      checks: block
         ! GATE: a dynamically-dry column receives no flux at all.
         call check(error, abs(dT_dry) == 0.0_wp .and. abs(dS_dry) == 0.0_wp, &
                    "sflux gate: dry column (wd_wet_dyn=0) received surface flux")
         if (allocated(error)) exit checks

         ! DEEP wet column warms by the analytic amount (gate passes it).
         write (msg, '("sflux gate: deep column ΔT=", es12.5, " != analytic ", es12.5)') &
            dT_deep, dT_deep_ana
         call check(error, abs(dT_deep - dT_deep_ana) <= 1.0e-12_wp*abs(dT_deep_ana), &
                    trim(msg))
         if (allocated(error)) exit checks

         ! THIN held-wet band: the binary gate applies the FULL flux
         ! into the mm-scale top layer — matches the UNTHROTTLED
         ! analytic value to round-off (NO thickness-aware protection).
         write (msg, '("sflux gate: thin-band ΔT=", es12.5, " != analytic ", es12.5, &
                       &" (h_top=", es10.3, " m)")') dT_thin, dT_thin_ana, htop_thin
         call check(error, abs(dT_thin - dT_thin_ana) <= 1.0e-12_wp*abs(dT_thin_ana), &
                    trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("sflux gate: thin-band ΔS=", es12.5, " != analytic ", es12.5)') &
            dS_thin, Q_SALT*inv_salt/htop_thin
         call check(error, abs(dS_thin - Q_SALT*inv_salt/htop_thin) <= &
                    1.0e-12_wp*abs(Q_SALT*inv_salt/htop_thin), trim(msg))
         if (allocated(error)) exit checks

         ! The gate gap is real and exercised: a held-wet band column
         ! warms by orders of magnitude more per step than a deep one
         ! (D_DEEP/D_THIN ~ 1400x here).  This documents S1 — if a
         ! thickness-aware throttle is ever added, this assertion (and
         ! the analytic-identity ones above) must be revisited.
         write (msg, '("sflux gate: thin/deep ΔT ratio=", es10.3, " (expected ~", es10.3, ")")') &
            dT_thin/dT_deep, htop_deep/htop_thin
         call check(error, dT_thin > 100.0_wp*dT_deep, trim(msg))
      end block checks

      deallocate (wet_dyn, hT_ic, hS_ic)
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_sflux_gate

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   pure real(wp) function DT_of(h_deep, dx, n_inner) result(dt)
      !! Outer dt sized so the inner barotropic substep sits well within
      !! the gravity-wave CFL for the deepest column.
      real(wp), intent(in) :: h_deep, dx
      integer, intent(in) :: n_inner
      dt = 0.4_wp*dx/sqrt(GRAVITY*h_deep)*real(n_inner, wp)
   end function DT_of

   subroutine enable_wetdry_full(dyn, grid, dry_depth, rewet_depth)
      !! Mirror `configure_ocean_wetdry` on a bare `ocean_dyn_t`: set the
      !! knobs, allocate the wd_* workspaces, seed the hysteresis mask
      !! from the seeded depth (eta = 0 reference => D0 = bt_H_ref).
      type(ocean_dyn_t), intent(inout) :: dyn
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: dry_depth, rewet_depth
      integer :: nx, ny, i, j
      nx = grid%nx_total
      ny = grid%ny_total
      dyn%bt_work%wetdry_enable = .true.
      dyn%bt_work%wd_dry_depth = dry_depth
      dyn%bt_work%wd_rewet_depth = rewet_depth
      allocate (dyn%bt_work%wd_wet_dyn(nx, ny), source=1.0_wp)
      allocate (dyn%bt_work%wd_theta(nx, ny), source=1.0_wp)
      allocate (dyn%bt_work%wd_flux_x(nx + 1, ny), source=0.0_wp)
      allocate (dyn%bt_work%wd_flux_y(nx, ny + 1), source=0.0_wp)
      allocate (dyn%bt_work%wd_open_u(nx + 1, ny), source=1.0_wp)
      allocate (dyn%bt_work%wd_open_v(nx, ny + 1), source=1.0_wp)
      do j = 1, ny
         do i = 1, nx
            if (dyn%bt_work%bt_H_ref(i, j) < dry_depth) then
               dyn%bt_work%wd_wet_dyn(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine enable_wetdry_full

   pure real(wp) function tracer_content(ms, idx, i0, i1, j0, j1, nz) result(c)
      !! Σ hTr over interior physical cells (all layers) — the conserved
      !! tracer CONTENT (concentration·thickness integrated).
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx, i0, i1, j0, j1, nz
      integer :: i, j, k
      c = 0.0_wp
      do k = 1, nz
         do j = j0, j1
            do i = i0, i1
               c = c + ms%tracers(idx)%hTr(i, j, k)
            end do
         end do
      end do
   end function tracer_content

   pure real(wp) function column_depth(ms, i, j, nz) result(d)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: i, j, nz
      integer :: k
      d = 0.0_wp
      do k = 1, nz
         d = d + ms%h_layer(i, j, k)
      end do
   end function column_depth

   logical function all_finite(ms, nz, culprit) result(ok)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: nz
      character(len=*), intent(out), optional :: culprit
      integer :: t
      ok = .true.
      if (present(culprit)) culprit = ""
      if (any(ieee_is_nan(ms%h_layer)) .or. .not. all(ieee_is_finite(ms%h_layer))) then
         ok = .false.; if (present(culprit)) culprit = "h_layer"; return
      end if
      if (any(ieee_is_nan(ms%rho_layer))) then
         ok = .false.; if (present(culprit)) culprit = "rho_layer"; return
      end if
      if (any(ieee_is_nan(ms%u_face_x_layer))) then
         ok = .false.; if (present(culprit)) culprit = "u_face_x_layer"; return
      end if
      if (any(ieee_is_nan(ms%v_face_y_layer))) then
         ok = .false.; if (present(culprit)) culprit = "v_face_y_layer"; return
      end if
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%hTr)) cycle
            if (any(ieee_is_nan(ms%tracers(t)%hTr))) then
               ok = .false.
               if (present(culprit)) write (culprit, '("tracer ", i0, " hTr")') t
               return
            end if
         end do
      end if
   end function all_finite

   subroutine set_column(ms, i, j, nz, htop, T_ref, S_ref)
      !! Uniform-sigma column of NZ layers each `htop` thick, tracers at
      !! reference concentration.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: i, j, nz
      real(wp), intent(in) :: htop, T_ref, S_ref
      integer :: k
      do k = 1, nz
         ms%h_layer(i, j, k) = htop
         ms%tracers(ms%idx_temperature)%hTr(i, j, k) = T_ref*htop
         ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S_ref*htop
      end do
   end subroutine set_column

   pure real(wp) function col_dtop(ms, ic, idx, i, ig, nyp, nz) result(dtop)
      !! Row-averaged top-layer (k = nz) concentration change for column
      !! `i`: mean over physical rows of (hTr - hTr_ic)/h_layer.
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: ic(:, :, :)
      integer, intent(in) :: idx, i, ig, nyp, nz
      integer :: j
      real(wp) :: acc, h
      acc = 0.0_wp
      do j = ig + 1, ig + nyp
         h = ms%h_layer(i, j, nz)
         if (h > 0.0_wp) then
            acc = acc + (ms%tracers(idx)%hTr(i, j, nz) - ic(i, j, nz))/h
         end if
      end do
      dtop = acc/real(nyp, wp)
   end function col_dtop

   subroutine run_apply_dyn(grid, sf, ms, wet_dyn, dt)
      !! Device round-trip for the wet/dry-gated surface-flux apply.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(inout) :: wet_dyn(:, :)
      real(wp), intent(in) :: dt
      !$acc enter data copyin(ms, sf, wet_dyn)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_tracers(grid, sf, ms, dt, wet_dyn=wet_dyn)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr)
         !$acc update self(hT, hS)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf, wet_dyn)
   end subroutine run_apply_dyn

   ! ---- slot map/destroy (trimmed from test_ocean_baroclinic_longrun) ----

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call sf%enter_data(); call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call destroy_cartesian_metrics(metrics)
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data(); call sf%exit_data()
      call ss%exit_data(); call bd%exit_data(); call hv%exit_data()
      call pgf%exit_data(); call cor%exit_data(); call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call sf%destroy(); call ss%destroy()
      call bd%destroy(); call hv%destroy(); call pgf%destroy(); call cor%destroy()
      call ct%destroy(); call ms%destroy()
   end subroutine destroy_all

end module test_ocean_wetdry_driver
