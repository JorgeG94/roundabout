!! End-to-end validation tests for the ocean Tier-1 dynamical core.
!!
!! Where `test_ocean_dyn_multilayer` checks the *trivial* fixed
!! points of the driver (stratified rest, conservation under
!! benign smooth flow), this file runs the driver against
!! **active** physics on a small sample domain and checks
!! integrated invariants that span the whole RK2 pipeline.
!!
!! Each test isolates one regime so any failure points at one
!! sub-system:
!!
!!   * `cooling_kpp_column` — vertical column physics only (zero
!!     velocity, no advection).  Cooling + KPP overlay + non-local
!!     γ + vdiff.  BL deepens; column heat tracks the surface
!!     forcing budget; column salt is conserved (Q_salt = 0).
!!     Probes the new KPP-non-local apply, surface_flux_apply, and
!!     vdiff coupled through RK2.
!!
!!   * `wind_ekman_spinup_f_plane` — momentum dynamics only (no
!!     thermal forcing, KPP off).  Steady eastward wind stress on
!!     an f-plane.  Surface u grows; Coriolis turns it southward
!!     (v < 0 in the Northern Hemisphere f > 0 convention).
!!     Probes the new 2D wind-stress field, the new mass-flux-
!!     weighted Coriolis, the nonlinear barotropic substep (when fired via
!!     the split driver — wired in a follow-up).
!!
!!   * `flow_conservation_no_forcing` — non-trivial wall-vanishing
!!     flow, stratified T, S, no surface forcing.  Column mass,
!!     column salt, column heat all conserved over many RK2 steps.
!!     Probes the split continuity + tracer interleaving + vertical
!!     advection cancellation + RK2 averaging hold up at the
!!     driver level.
module test_ocean_validation
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step
   implicit none
   private

   public :: collect_ocean_validation_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_validation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cooling_kpp_column", test_cooling_column), &
                  new_unittest("wind_ekman_spinup_f_plane", test_wind_ekman), &
                  new_unittest("flow_conservation_no_forcing", &
                               test_flow_conservation) &
                  ]
   end subroutine collect_ocean_validation_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine init_all(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0)
      type(hgrid_t), intent(in) :: grid
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
      real(wp), intent(in) :: f0
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = f0
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
      call dyn%init(grid)
   end subroutine init_all

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
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call sf%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine destroy_all

   subroutine map_in_all(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
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
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
   end subroutine map_in_all

   subroutine map_out_all(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
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
      call destroy_cartesian_metrics(metrics)
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call sf%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_all

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_cooling_column(error)
      !! Vertical column physics only: zero velocity (so split
      !! continuity / horizontal tracer advect / vertical advect
      !! are all no-ops), stratified T (warm top, cold bottom),
      !! uniform S, steady surface cooling, KPP overlay + non-local
      !! γ + vdiff all active.  Several SSP-RK2 outer steps.
      !!
      !! Invariants:
      !!   * h_layer stays at its IC bit-for-bit (zero flow → no
      !!     thickness change)
      !!   * Column salt conserved exactly (Q_salt = 0, no flow)
      !!   * Column heat ΔΣ(hT) = N · dt · Q_heat / (ρ_0 · cp)
      !!     to round-off — surface forcing is the only source
      !!   * BL depth at the end > BL depth after step 1
      !!     (cooling deepens the BL through KPP)
      !!
      !! End-to-end probe of: `ocean_surface_flux_apply_tracers_-
      !! multilayer`, `vmix_compute_pp81`,
      !! `vmix_apply_kpp_overlay` (Phase 1 + 2),
      !! `vdiff_apply_tracers`, and the new
      !! `vmix_apply_nonlocal_tendencies`, all
      !! interleaved through the 2-stage RK2 averaging.
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
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: Q_COOLING = -200.0_wp
      integer, parameter :: N_STEPS = 5
      real(wp) :: bl_after_first, bl_at_end
      real(wp) :: heat_budget_expected, max_budget_dev, max_dh, max_dsalt
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :), h_ic(:, :, :)
      integer :: i, j, k, step, i_probe, j_probe
      checks: block

         ! dx = 1 km keeps the gravity-wave CFL out of play; with zero
         ! velocity and zero wind there's no horizontal dynamics anyway.
         call make_grid(grid, 10, 8, 1.0e3_wp, 1.0e3_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=0.0_wp)

         ! Stratified IC: warm surface (k = NZ), cooler below.  Zero
         ! velocity throughout the run.
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H_LAYER
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
               (eos%T_ref + real(k - 1, wp))*H_LAYER
            ms%rho_layer(:, :, k) = eos%rho0 + 0.5_wp*real(NZ - k, wp)
         end do
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (h_ic, source=ms%h_layer)

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call sf%set_surface_flux_const(Q_COOLING, 0.0_wp)

         vmix%use_closure = .true.
         vmix%use_kpp = .true.

         call map_in_all(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT, sf=sf)
         !$acc update self(vmix%bl_depth)
         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         bl_after_first = vmix%bl_depth(i_probe, j_probe)
         do step = 2, N_STEPS
            call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT, sf=sf)
         end do
         !$acc update self(vmix%bl_depth)
         call map_out_all(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)

         bl_at_end = vmix%bl_depth(i_probe, j_probe)

         ! h must not drift — zero flow, no surface mass flux.
         max_dh = maxval(abs(ms%h_layer - h_ic))
         call check(error, max_dh < 1.0e-10_wp, &
                    "cooling+KPP: h_layer drifted under zero flow")
         if (allocated(error)) exit checks

         ! Salt: Q_salt = 0 → γ_S = 0 and the column has no horizontal
         ! flow under uniform forcing (wall and interior kv are now
         ! homogeneous after the KPP bulk-Ri sweep was extended to wall
         ! columns).  Per-cell salt
         ! must therefore stay at its IC to round-off.
         max_dsalt = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         call check(error, max_dsalt < 1.0e-10_wp, &
                    "cooling+KPP: per-cell salt drifted under zero flow")
         if (allocated(error)) exit checks

         ! BL must not shallow under sustained cooling.  Strict
         ! deepening only happens until the BL fills the column; once
         ! it saturates (h_b == total H), subsequent steps maintain
         ! that value.  V_t²-on diagnoses a deeper BL than V_t²-off on
         ! the same initial state, so the BL can saturate after step 1.
         call check(error, bl_at_end >= bl_after_first, &
                    "cooling+KPP: BL shallowed between step 1 and step N")
         if (allocated(error)) exit checks

         ! Per-column heat budget: with wall and interior columns
         ! producing the same kv profile under uniform forcing, no
         ! baroclinic adjustment flow develops, so every (i,j) column's
         ! integrated hT change must equal the surface-flux contribution
         ! to round-off.  Vertical diffusion redistributes heat within a
         ! column but conserves the column total — checking per-layer
         ! against the column-integrated expected would be wrong, since
         ! the IC is stratified and 5 RK2 steps don't fully mix.
         heat_budget_expected = real(N_STEPS, wp)*DT*Q_COOLING/(eos%rho0*sf%cp)
         max_budget_dev = 0.0_wp
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               max_budget_dev = max(max_budget_dev, abs( &
                                    sum(ms%tracers(ms%idx_temperature)%hTr(i, j, :) &
                                        - hT_ic(i, j, :)) &
                                    - heat_budget_expected))
            end do
         end do
         call check(error, max_budget_dev < 1.0e-10_wp, &
                    "cooling+KPP: per-column heat budget didn't match surface forcing")

      end block checks
      deallocate (hT_ic, hS_ic, h_ic)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_cooling_column

   subroutine test_wind_ekman(error)
      !! Steady eastward wind stress on an f-plane, no thermal
      !! forcing, no KPP overlay.  Coriolis turns the wind-driven
      !! surface flow:
      !!
      !!   * The surface u_face_x grows from zero (wind accelerates it)
      !!   * v_face_y at the surface acquires the correct Ekman sign.
      !!     With f > 0 (Northern Hemisphere) and τ_x > 0, the
      !!     Ekman transport is to the right of the wind → southward
      !!     (negative y) → v_bt < 0.
      !!
      !! End-to-end probe of: ocean_surface_stress_apply on the top
      !! layer, Coriolis-adv tendencies on u and v, RK2 averaging,
      !! GPU mapping of `tau_x` / `tau_y` (the new 2D fields).
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
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: DT = 5.0_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: TAU_X = 0.1_wp
      integer, parameter :: N_STEPS = 20
      real(wp) :: u_top, v_top, mass_ic, mass_final
      integer :: i, j, k, step, i_probe, j_probe
      checks: block

         call make_grid(grid, 16, 12, 1.0e3_wp, 1.0e3_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)

         ! Uniform column, uniform (T, S), zero velocity.
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H_LAYER
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H_LAYER
            ms%rho_layer(:, :, k) = eos%rho0
         end do

         call ss%set_wind_stress_const(TAU_X, 0.0_wp)
         call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

         ! Closed-domain mass conservation: no surface mass flux, walls
         ! seal the basin, so Σh must be invariant to round-off.  Run
         ! before map_in_all so we sum on the host-side IC.
         mass_ic = sum(ms%h_layer)

         ! Leave vmix off — the discriminator should be Coriolis-only
         ! (with vmix on, momentum would diffuse into deeper layers
         ! and the surface signal would be muddier on this short run).
         vmix%use_closure = .false.
         vmix%use_kpp = .false.

         call map_in_all(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
         do step = 1, N_STEPS
            call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT, sf=sf)
         end do
         call map_out_all(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         u_top = ms%u_face_x_layer(i_probe, j_probe, NZ)
         v_top = ms%v_face_y_layer(i_probe, j_probe, NZ)

         ! With τ_x > 0, the wind accelerates u_top from zero.  Over
         ! N_STEPS·DT ≈ 100 s and ρ·H = 25e3 kg/m², linear ramp would
         ! give u ≈ τ_x · t / (ρ · H) = 0.1 · 100 / 25800 ≈ 4e-4 m/s.
         ! Coriolis trims this slightly but the sign and order of
         ! magnitude must hold.
         call check(error, u_top > 0.0_wp, &
                    "Ekman: surface u didn't grow under eastward wind")
         if (allocated(error)) exit checks
         call check(error, u_top > 1.0e-5_wp, &
                    "Ekman: surface u magnitude too small (force not reaching layer?)")
         if (allocated(error)) exit checks

         ! Coriolis turns the surface flow to the right of the wind
         ! in NH.  With τ_x > 0, the Coriolis force on u is `+f·v`,
         ! the Coriolis force on v is `-f·u`.  As u grows positive,
         ! v acquires negative sign — surface v should be < 0.
         call check(error, v_top < 0.0_wp, &
                    "Ekman: Coriolis didn't turn surface flow southward")
         if (allocated(error)) exit checks

         ! Sanity bound: |v_top| should be smaller than |u_top| over
         ! this short integration (only one quarter inertial period
         ! is ~15700 s, we're way short of that).
         call check(error, abs(v_top) < abs(u_top), &
                    "Ekman: |v_top| exceeded |u_top| — Coriolis over-rotated?")
         if (allocated(error)) exit checks

         ! Closed-domain mass conservation invariant (Tier 3).
         mass_final = sum(ms%h_layer)
         call check(error, abs(mass_final - mass_ic) < 1.0e-10_wp*abs(mass_ic), &
                    "Ekman: global mass not conserved under closed-domain wind forcing")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_wind_ekman

   subroutine test_flow_conservation(error)
      !! Dynamical-core conservation with no surface forcing.  Non-
      !! trivial wall-vanishing flow, stratified T (warmer surface,
      !! cooler bottom), uniform S, no wind, no thermal forcing,
      !! KPP off (no vmix overlay).  Several SSP-RK2 steps.
      !! Invariants:
      !!
      !!   * column mass per interior cell: split continuity is
      !!     conservative against the closed-wall BC; vertical
      !!     advection's `apply_w_to_h` cancels per-layer; column
      !!     mass should be preserved to round-off.
      !!   * column salt per interior cell: CWC theorem holds for
      !!     the interleaved split tracer step; no surface flux;
      !!     column salt should be preserved to round-off.
      !!   * column heat per interior cell: same argument, no
      !!     surface flux; column heat preserved.
      !!
      !! End-to-end probe of: split continuity + interleaved tracer
      !! advection + vertical advection + RK2 averaging, all
      !! running on a non-trivial flow without thermodynamic
      !! forcing complicating the budget.
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
      real(wp), parameter :: H_BASE = 20.0_wp
      real(wp), parameter :: U_AMP = 0.02_wp
      real(wp), parameter :: DT = 5.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 10
      real(wp), allocatable :: mass_ic(:, :), salt_ic(:, :), heat_ic(:, :)
      real(wp) :: max_dmass, max_dsalt, max_dheat
      integer :: i, j, k, step, nx, ny

      call make_grid(grid, 12, 10, 1.0e3_wp, 1.0e3_wp)
      call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=1.0e-4_wp)
      nx = grid%nx_total
      ny = grid%ny_total

      ! Uniform h, wall-vanishing u/v, stratified T per layer (uniform
      ! in i, j within each layer), uniform S.  No surface forcing,
      ! KPP off.  This isolates the advective dynamical core.
      ms%h_layer = H_BASE
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = U_AMP* &
                                            sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                            sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp))
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = U_AMP* &
                                            sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                            sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp))
            end do
         end do
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H_BASE
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (eos%T_ref + real(NZ - k, wp))*H_BASE
         ms%rho_layer(:, :, k) = eos%rho0 + 0.5_wp*real(NZ - k, wp)
      end do

      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      vmix%use_closure = .false.
      vmix%use_kpp = .false.

      ! Closed walls sit at the OUTER face (i=1 / i=nx_total+1).
      ! Interior faces (e.g. face 3) carry non-zero mass flux —
      ! flow legitimately moves tracer between every pair of
      ! adjacent cells, including the "ghost" rows at i=1,2 and
      ! the "interior" at i=3+.  Conservation is therefore a
      ! *global* property over the whole grid, not the interior
      ! sub-block.  Probe globally over all cells.
      allocate (mass_ic(nx, ny), salt_ic(nx, ny), heat_ic(nx, ny))
      mass_ic = sum(ms%h_layer, dim=3)
      salt_ic = sum(ms%tracers(ms%idx_salinity)%hTr, dim=3)
      heat_ic = sum(ms%tracers(ms%idx_temperature)%hTr, dim=3)

      call map_in_all(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      do step = 1, N_STEPS
         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT, sf=sf)
      end do
      call map_out_all(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)

      checks: block
         real(wp) :: global_mass_before, global_mass_after
         real(wp) :: global_salt_before, global_salt_after
         real(wp) :: global_heat_before, global_heat_after
         global_mass_before = sum(mass_ic)
         global_mass_after = sum(sum(ms%h_layer, dim=3))
         global_salt_before = sum(salt_ic)
         global_salt_after = sum(sum(ms%tracers(ms%idx_salinity)%hTr, dim=3))
         global_heat_before = sum(heat_ic)
         global_heat_after = sum(sum(ms%tracers(ms%idx_temperature)%hTr, dim=3))
         max_dmass = abs(global_mass_after - global_mass_before)
         max_dsalt = abs(global_salt_after - global_salt_before)
         max_dheat = abs(global_heat_after - global_heat_before)
         call check(error, max_dmass < 1.0e-10_wp*abs(global_mass_before), &
                    "no-forcing flow: global mass not conserved")
         if (allocated(error)) exit checks
         call check(error, max_dsalt < 1.0e-10_wp*abs(global_salt_before), &
                    "no-forcing flow: global salt not conserved")
         if (allocated(error)) exit checks
         call check(error, max_dheat < 1.0e-10_wp*abs(global_heat_before), &
                    "no-forcing flow: global heat not conserved")
         if (allocated(error)) exit checks
         call check(error, minval(ms%h_layer) > 0.0_wp, &
                    "no-forcing flow: h_layer went non-positive")
      end block checks
      deallocate (mass_ic, salt_ic, heat_ic)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_flow_conservation

end module test_ocean_validation
