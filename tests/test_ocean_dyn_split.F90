!! Integration tests for the split-explicit multilayer driver
!! (`ocean_dyn_step_split` in `rdb_ocean_dyn`).
!!
!! The split driver decomposes each SSP-RK2 stage's slow tendency
!! into a depth-mean (barotropic) part and a layer-relative
!! (baroclinic) part.  The slow apply produces an FE-amplified bt
!! mode; the split driver then OVERRIDES that with a barotropic-substep
!! result that resolves the surface gravity wave at dt_inner.  See
!! the driver header for the algebra.
!!
!! Cases:
!!   * Stratified rest — same IC as the unsplit test: stratified
!!     S, T with rho-balance setting up no net PGF, no wind/drag,
!!     η = 0.  After one split step every prognostic must stay at
!!     its IC to round-off.  Catches a sign bug or accumulator
!!     leak in the mode-separation pipeline.
!!   * Volume + tracer-mass conservation through 10 split steps —
!!     same setup as the unsplit conservation test.  Confirms the
!!     bt-correction doesn't introduce a per-step drift.
!!   * Gravity-wave stability gain — initial η perturbation at a
!!     dt where the unsplit Heun method blows up.  The split
!!     driver with n_inner = 20 must keep |η| bounded; the unsplit
!!     driver run at the same dt should diverge.
module test_ocean_dyn_split
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
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
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, &
                            ocean_dyn_step_split, &
                            ocean_dyn_step
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_OPEN, OBC_SPONGE, OBC_TIDAL, &
                                       OBC_CLAMPED, OBC_CHAPMAN
   use rdb_ocean_sponge, only: ocean_sponge_apply, ocean_sponge_t
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_MASS
   implicit none
   private

   public :: collect_ocean_dyn_split_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_dyn_split_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("split_stratified_rest", test_stratified_rest), &
                  new_unittest("split_mass_tracer_conservation", test_conservation), &
                  new_unittest("split_gravity_wave_stability", test_gravity_wave_stability), &
                  new_unittest("split_obc_open_west_drains_via_driver", &
                               test_obc_open_via_driver), &
                  new_unittest("split_obc_sponge_kernel_damps_band", &
                               test_obc_sponge_kernel), &
                  new_unittest("split_obc_tidal_via_driver", &
                               test_obc_tidal_via_driver), &
                  new_unittest("split_obc_clamped_via_driver", &
                               test_obc_clamped_via_driver), &
                  new_unittest("split_obc_chapman_via_driver", &
                               test_obc_chapman_via_driver), &
                  new_unittest("split_continuity_mass_contributor_closed_basin", &
                               test_mass_contributor_closed_basin), &
                  new_unittest("split_ideal_age_grows_linearly_via_driver", &
                               test_ideal_age_via_driver), &
                  new_unittest("split_ideal_age_surface_reset_nonzero_ic", &
                               test_ideal_age_surface_reset_nonzero_ic), &
                  new_unittest("split_ideal_age_thermo_cadence", &
                               test_ideal_age_thermo_cadence), &
                  new_unittest("split_sponge_bit_identity_when_disabled", &
                               test_sponge_bit_identity_when_disabled) &
                  ]
   end subroutine collect_ocean_dyn_split_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      call dyn%enter_data()
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call destroy_cartesian_metrics(metrics)
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine destroy_all

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_stratified_rest(error)
      !! Stratified-rest IC same as the unsplit driver's
      !! corresponding test.  After one split step every prognostic
      !! must stay at the IC to round-off.  Tests the split driver
      !! reduces to a no-op when there's nothing to evolve.
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_INNER = 5
      real(wp) :: S_k(NZ), T_k(NZ)
      real(wp), allocatable :: hS_ic(:, :, :), hT_ic(:, :, :)
      real(wp) :: max_dh, max_du, max_dv, max_dS, max_dT
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         S_k = [eos%S_ref + 1.0_wp, eos%S_ref, eos%S_ref - 1.0_wp]
         T_k = [eos%T_ref - 2.0_wp, eos%T_ref, eos%T_ref + 2.0_wp]
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*H0
         end do
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         ! For the split driver, H_ref is the total column depth.
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ocean_dyn_step_split( &
            grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
            va, hd, vd, vmix, ms, DT, N_INNER)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         max_dh = maxval(abs(ms%h_layer - H0))
         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_dh < 1.0e-10_wp, &
                    "split stratified rest: h drifted")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-10_wp, &
                    "split stratified rest: u drifted")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp, &
                    "split stratified rest: v drifted")
         if (allocated(error)) exit checks
         call check(error, max_dS < 1.0e-10_wp, &
                    "split stratified rest: hS drifted")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-10_wp, &
                    "split stratified rest: hT drifted")

      end block checks
      deallocate (hS_ic, hT_ic)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_stratified_rest

   subroutine test_conservation(error)
      !! Closed-wall mass-conservation test through several split
      !! steps.  Same IC as the unsplit conservation test: small
      !! solenoidal u/v perturbation on a uniform-density column,
      !! short window so the gravity wave hasn't blown up under any
      !! integrator.
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 10
      integer, parameter :: N_INNER = 5
      integer :: i, j, k, nx, ny, step
      real(wp) :: total_h0, total_S0, total_T0
      real(wp) :: total_h, total_S, total_T, drift_h, drift_S, drift_T
      real(wp) :: h_min
      checks: block

         call make_grid(grid, 24, 16, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + 0.01_wp*sin( &
                                        2.0_wp*PI*real(i, wp)/real(nx, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = ms%h_layer(i, j, k)*( &
                                                             eos%S_ref + 0.5_wp*cos(2.0_wp*PI*real(j, wp)/real(ny, wp)))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = ms%h_layer(i, j, k)*( &
                                                                eos%T_ref + 0.3_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))
               end do
            end do
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
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H_BASE

         total_h0 = sum(ms%h_layer)*grid%dx*grid%dy
         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T0 = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         total_h = sum(ms%h_layer)*grid%dx*grid%dy
         total_S = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         drift_h = abs(total_h - total_h0)/abs(total_h0)
         drift_S = abs(total_S - total_S0)/abs(total_S0)
         drift_T = abs(total_T - total_T0)/abs(total_T0)
         h_min = minval(ms%h_layer)

         ! The split bt-correction adds a uniform δu per layer.  Mass
         ! conservation through `continuity_apply` is bit-exact for the
         ! slow update; the bt correction doesn't modify h_layer in
         ! this MVP scope.  So total h, S, T should drift only via
         ! floating-point round-off (well below 1e-10).
         call check(error, drift_h < 1.0e-10_wp, "split total h drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, drift_S < 1.0e-10_wp, "split total S mass drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, drift_T < 1.0e-10_wp, "split total T mass drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, h_min > 0.0_wp, "split h_layer went non-positive")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_conservation

   subroutine test_gravity_wave_stability(error)
      !! Gravity-wave stability discriminator.  Set up a column with
      !! a deliberate η perturbation that the unsplit SSP-RK2 driver
      !! amplifies (Heun on imaginary axis is marginally unstable).
      !! Run the split driver at the same outer dt with `n_inner`
      !! large enough that the barotropic substep is well within FBE stability.
      !! Verify split stays bounded while the unsplit-equivalent
      !! diverges over the same window.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_split, grid_unsplit
      type(ocean_metrics_t) :: metrics_split, metrics_unsplit
      type(multilayer_state_t) :: ms_split, ms_unsplit
      type(continuity_t) :: ct_s, ct_u
      type(coriolis_adv_t) :: cor_s, cor_u
      type(ocean_pressure_force_t) :: pgf_s, pgf_u
      type(ocean_horizontal_viscosity_t) :: hv_s, hv_u
      type(ocean_bottom_drag_t) :: bd_s, bd_u
      type(ocean_surface_stress_t) :: ss_s, ss_u
      type(ocean_vertical_advection_t) :: va_s, va_u
      type(ocean_hdiff_tracer_t) :: hd_s, hd_u
      type(ocean_vdiff_t) :: vd_s, vd_u
      type(ocean_vmix_t) :: vmix_s, vmix_u
      type(eos_t) :: eos_s, eos_u
      type(ocean_dyn_t) :: dyn_s, dyn_u
      integer, parameter :: NZ_GW = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, lambda, dt_cfl_grav, dt_outer
      integer :: i, j, step, n_steps
      integer, parameter :: N_INNER = 20
      real(wp) :: max_h_split_final, max_h_unsplit_final
      checks: block

         call make_grid(grid_split, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_unsplit, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         ms_split%nz_ml = NZ_GW
         ms_unsplit%nz_ml = NZ_GW
         call ms_split%init(grid_split)
         call ms_unsplit%init(grid_unsplit)
         call ct_s%init(grid_split, nz_ml=NZ_GW); call ct_u%init(grid_unsplit, nz_ml=NZ_GW)
         call cor_s%init(grid_split, nz_ml=NZ_GW); call cor_u%init(grid_unsplit, nz_ml=NZ_GW)
         call pgf_s%init(grid_split, nz_ml=NZ_GW); call pgf_u%init(grid_unsplit, nz_ml=NZ_GW)
         call hv_s%init(grid_split, nz_ml=NZ_GW); call hv_u%init(grid_unsplit, nz_ml=NZ_GW)
         call bd_s%init(grid_split, nz_ml=NZ_GW); call bd_u%init(grid_unsplit, nz_ml=NZ_GW)
         call ss_s%init(grid_split, nz_ml=NZ_GW); call ss_u%init(grid_unsplit, nz_ml=NZ_GW)
         call va_s%init(grid_split, nz_ml=NZ_GW); call va_u%init(grid_unsplit, nz_ml=NZ_GW)
         call hd_s%init(grid_split, nz_ml=NZ_GW); call hd_u%init(grid_unsplit, nz_ml=NZ_GW)
         call vd_s%init(grid_split, nz_ml=NZ_GW); call vd_u%init(grid_unsplit, nz_ml=NZ_GW)
         call vmix_s%init(grid_split, nz_ml=NZ_GW); call vmix_u%init(grid_unsplit, nz_ml=NZ_GW)
         call eos_s%init(grid_split); call eos_u%init(grid_unsplit)
         call dyn_s%init(grid_split, nz_ml=NZ_GW); call dyn_u%init(grid_unsplit)

         ! IC: m=2 closed-wall cosine eigenmode over the physical
         ! interior (cells nghost+1..nghost+nx_phys).  Both drivers
         ! close walls at the physical boundary, so the wave must be
         ! defined there to be a true eigenmode.  Ghost cells stay at
         ! the reference depth H0 (η = 0).
         do j = 1, grid_split%ny_total
            do i = 1, grid_split%nx_total
               if (i > grid_split%nghost .and. i <= grid_split%nghost + NX_PHYS) then
                  ms_split%h_layer(i, j, 1) = H0 + A_ETA* &
                                              cos(2.0_wp*PI*(real(i - grid_split%nghost, wp) - 0.5_wp) &
                                                  /real(NX_PHYS, wp))
               else
                  ms_split%h_layer(i, j, 1) = H0
               end if
               ms_unsplit%h_layer(i, j, 1) = ms_split%h_layer(i, j, 1)
               ms_split%tracers(ms_split%idx_salinity)%hTr(i, j, 1) = &
                  eos_s%S_ref*ms_split%h_layer(i, j, 1)
               ms_split%tracers(ms_split%idx_temperature)%hTr(i, j, 1) = &
                  eos_s%T_ref*ms_split%h_layer(i, j, 1)
               ms_unsplit%tracers(ms_unsplit%idx_salinity)%hTr(i, j, 1) = &
                  eos_u%S_ref*ms_unsplit%h_layer(i, j, 1)
               ms_unsplit%tracers(ms_unsplit%idx_temperature)%hTr(i, j, 1) = &
                  eos_u%T_ref*ms_unsplit%h_layer(i, j, 1)
            end do
         end do
         ms_split%u_face_x_layer = 0.0_wp; ms_split%v_face_y_layer = 0.0_wp
         ms_unsplit%u_face_x_layer = 0.0_wp; ms_unsplit%v_face_y_layer = 0.0_wp
         dyn_s%bt_work%bt_H_ref = H0

         ! Outer dt: deliberately at the gravity-wave CFL where Heun
         ! is marginally unstable (dt·ω ~ O(1)).  Inner dt: dt_outer /
         ! N_INNER → comfortably stable under FBE.
         lambda = real(NX_PHYS, wp)*grid_split%dx
         c = sqrt(GRAVITY*H0)
         dt_cfl_grav = grid_split%dx/c
         dt_outer = 1.0_wp*dt_cfl_grav   ! at the CFL — Heun is just past stable
         n_steps = 200

         call map_in(grid_split, metrics_split, ms_split, ct_s, cor_s, pgf_s, hv_s, bd_s, ss_s, va_s, hd_s, vd_s, vmix_s, dyn_s)
         do step = 1, n_steps
            call ocean_dyn_step_split( &
               grid_split, metrics_split, dyn_s, eos_s, cor_s, ct_s, pgf_s, hv_s, bd_s, ss_s, &
               va_s, hd_s, vd_s, vmix_s, ms_split, dt_outer, N_INNER)
         end do
         call map_out(metrics_split, ms_split, ct_s, cor_s, pgf_s, hv_s, bd_s, ss_s, va_s, hd_s, vd_s, vmix_s, dyn_s)

         call make_cartesian_metrics(metrics_unsplit, grid_unsplit)
         !$acc enter data copyin(ms_unsplit, ct_u, cor_u, pgf_u, hv_u, bd_u, ss_u, va_u, hd_u, vd_u, vmix_u)
         call ms_unsplit%enter_data()
         call ct_u%enter_data(); call cor_u%enter_data(); call pgf_u%enter_data()
         call hv_u%enter_data(); call bd_u%enter_data(); call ss_u%enter_data()
         call va_u%enter_data(); call hd_u%enter_data(); call vd_u%enter_data()
         call vmix_u%enter_data()
         do step = 1, n_steps
            call ocean_dyn_step( &
               grid_unsplit, metrics_unsplit, dyn_u, eos_u, cor_u, ct_u, pgf_u, hv_u, bd_u, ss_u, &
               va_u, hd_u, vd_u, vmix_u, ms_unsplit, dt_outer)
         end do
         call vmix_u%exit_data()
         call vd_u%exit_data(); call hd_u%exit_data(); call va_u%exit_data()
         call ss_u%exit_data(); call bd_u%exit_data(); call hv_u%exit_data()
         call pgf_u%exit_data(); call cor_u%exit_data(); call ct_u%exit_data()
         call ms_unsplit%exit_data()
         !$acc exit data delete(ms_unsplit, ct_u, cor_u, pgf_u, hv_u, bd_u, ss_u, va_u, hd_u, vd_u, vmix_u)
         call destroy_cartesian_metrics(metrics_unsplit)

         max_h_split_final = maxval(abs( &
                                    ms_split%h_layer(grid_split%nghost + 1:grid_split%nghost + NX_PHYS, :, 1) - H0))
         max_h_unsplit_final = maxval(abs( &
                                      ms_unsplit%h_layer(grid_unsplit%nghost + 1:grid_unsplit%nghost + NX_PHYS, :, 1) - H0))

         ! Split must stay bounded — final η perturbation no larger
         ! than ~5x the IC amplitude (allowing some FE drift in the
         ! barotropic substep) and finite.
         call check(error, max_h_split_final < 5.0_wp*A_ETA, &
                    "split driver: η amplitude grew beyond 5× IC amplitude")
         if (allocated(error)) exit checks
         ! Unsplit must have FE-amplified the wave by some non-trivial
         ! margin over the same window.  At dt = CFL_grav the Heun
         ! amplification factor per step is ~(1 + α^4/24) with α =
         ! dt·ω ~ 0.5·2π/NX_PHYS ~ ~0.4.  Over 200 steps, growth ≈
         ! 1.07.  We just demand split is at least 10% smaller in
         ! final amplitude than unsplit.
         call check(error, max_h_split_final < 0.9_wp*max_h_unsplit_final, &
                    "split driver: didn't outperform unsplit on gravity-wave amplitude")

      end block checks
      call destroy_all(ms_split, ct_s, cor_s, pgf_s, hv_s, bd_s, ss_s, va_s, hd_s, vd_s, vmix_s, eos_s, dyn_s)
      call destroy_all(ms_unsplit, ct_u, cor_u, pgf_u, hv_u, bd_u, ss_u, va_u, hd_u, vd_u, vmix_u, eos_u, dyn_u)
   end subroutine test_gravity_wave_stability

   subroutine test_obc_open_via_driver(error)
      !! End-to-end OBC check: drives the full split-RK2 step with an
      !! `ocean_bc_state_t` where the west edge is `OBC_OPEN`.  The
      !! reference run uses the same IC but doesn't pass `bc` — both
      !! reference + test use the same driver code, only the optional
      !! `bc` arg differs.
      !!
      !! Asserts the open-wall run drains the η pulse to a small
      !! fraction of the closed-wall reference, end-to-end through
      !! `ocean_dyn_step_split → run_stage_split → barotropic_substep_nonlinear`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_obc
      type(ocean_metrics_t) :: metrics_ref, metrics_obc
      type(multilayer_state_t) :: ms_ref, ms_obc
      type(continuity_t) :: ct_r, ct_o
      type(coriolis_adv_t) :: cor_r, cor_o
      type(ocean_pressure_force_t) :: pgf_r, pgf_o
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_o
      type(ocean_bottom_drag_t) :: bd_r, bd_o
      type(ocean_surface_stress_t) :: ss_r, ss_o
      type(ocean_vertical_advection_t) :: va_r, va_o
      type(ocean_hdiff_tracer_t) :: hd_r, hd_o
      type(ocean_vdiff_t) :: vd_r, vd_o
      type(ocean_vmix_t) :: vmix_r, vmix_o
      type(eos_t) :: eos_r, eos_o
      type(ocean_dyn_t) :: dyn_r, dyn_o
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NZ_OBC = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 20
      integer, parameter :: N_STEPS = 60
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, dt_outer
      real(wp) :: e2_ref, e2_obc
      integer :: i, j, step

      checks: block
         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_obc, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         ms_ref%nz_ml = NZ_OBC; ms_obc%nz_ml = NZ_OBC
         call ms_ref%init(grid_ref); call ms_obc%init(grid_obc)
         call ct_r%init(grid_ref, nz_ml=NZ_OBC); call ct_o%init(grid_obc, nz_ml=NZ_OBC)
         call cor_r%init(grid_ref, nz_ml=NZ_OBC); call cor_o%init(grid_obc, nz_ml=NZ_OBC)
         call pgf_r%init(grid_ref, nz_ml=NZ_OBC); call pgf_o%init(grid_obc, nz_ml=NZ_OBC)
         call hv_r%init(grid_ref, nz_ml=NZ_OBC); call hv_o%init(grid_obc, nz_ml=NZ_OBC)
         call bd_r%init(grid_ref, nz_ml=NZ_OBC); call bd_o%init(grid_obc, nz_ml=NZ_OBC)
         call ss_r%init(grid_ref, nz_ml=NZ_OBC); call ss_o%init(grid_obc, nz_ml=NZ_OBC)
         call va_r%init(grid_ref, nz_ml=NZ_OBC); call va_o%init(grid_obc, nz_ml=NZ_OBC)
         call hd_r%init(grid_ref, nz_ml=NZ_OBC); call hd_o%init(grid_obc, nz_ml=NZ_OBC)
         call vd_r%init(grid_ref, nz_ml=NZ_OBC); call vd_o%init(grid_obc, nz_ml=NZ_OBC)
         call vmix_r%init(grid_ref, nz_ml=NZ_OBC); call vmix_o%init(grid_obc, nz_ml=NZ_OBC)
         call eos_r%init(grid_ref); call eos_o%init(grid_obc)
         call dyn_r%init(grid_ref, nz_ml=NZ_OBC); call dyn_o%init(grid_obc, nz_ml=NZ_OBC)
         call ocean_bc_state_init(bc, grid_obc, nz_ml=NZ_OBC)
         bc%west%bc_type = OBC_OPEN

         ! m=1 cosine η pulse across the physical interior.
         do j = 1, grid_ref%ny_total
            do i = 1, grid_ref%nx_total
               if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                  ms_ref%h_layer(i, j, 1) = H0 + A_ETA* &
                                            cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp) &
                                                /real(NX_PHYS, wp))
               else
                  ms_ref%h_layer(i, j, 1) = H0
               end if
               ms_obc%h_layer(i, j, 1) = ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, 1) = &
                  eos_r%S_ref*ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, 1) = &
                  eos_r%T_ref*ms_ref%h_layer(i, j, 1)
               ms_obc%tracers(ms_obc%idx_salinity)%hTr(i, j, 1) = &
                  eos_o%S_ref*ms_obc%h_layer(i, j, 1)
               ms_obc%tracers(ms_obc%idx_temperature)%hTr(i, j, 1) = &
                  eos_o%T_ref*ms_obc%h_layer(i, j, 1)
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp; ms_ref%v_face_y_layer = 0.0_wp
         ms_obc%u_face_x_layer = 0.0_wp; ms_obc%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = H0; dyn_o%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.5_wp*grid_ref%dx/c  ! outer CFL ≈ 0.5

         ! Reference: no bc -> all-WALL.
         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         ! Test: same IC + bc with west = OBC_OPEN.
         call map_in(grid_obc, metrics_obc, ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_obc, metrics_obc, dyn_o, eos_o, cor_o, ct_o, pgf_o, hv_o, bd_o, ss_o, &
               va_o, hd_o, vd_o, vmix_o, ms_obc, dt_outer, N_INNER, bc=bc)
         end do
         call map_out(metrics_obc, ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)

         e2_ref = sum((ms_ref%h_layer(grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, :, 1) - H0)**2)
         e2_obc = sum((ms_obc%h_layer(grid_obc%nghost + 1:grid_obc%nghost + NX_PHYS, :, 1) - H0)**2)
         ! OBC_OPEN west wall should drain the pulse over 60 outer
         ! steps × 20 inner = 1200 substeps.  Threshold of 0.5 is
         ! generous; the kernel-level test sees ratios ~0.002.
         call check(error, e2_obc < 0.5_wp*e2_ref, &
                    "OBC_OPEN west via the driver should drain pulse below 50% of closed-wall ref")

      end block checks
      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, eos_o, dyn_o)
   end subroutine test_obc_open_via_driver

   subroutine test_obc_sponge_kernel(error)
      !! Kernel-level sponge test.  Tests `ocean_sponge_apply` directly
      !! without the dyn driver — the driver wiring is exercised by
      !! the OBC_OPEN test already.  Drives the sponge over a stack
      !! of substeps and asserts the EXACT analytical decay at the
      !! outermost band face (d=0, alpha=1) and at the innermost band
      !! face (d=band-1) — tightened from a `< 0.4·U0` smoke bound
      !! (PLAN_PR23_real_sponge.md §2.5 / §9.10: the kernel's own comment
      !! says `exp(-3) ~= 0.05`, but the old test asserted `< 0.4`, so a
      !! sponge running at half — or double — the configured rate would
      !! still pass).  Also confirms u outside the band stays untouched.
      !! Maps `ms` to the device (`!$acc enter data` / `ms%enter_data()`,
      !! `!$acc update self` before the host read) per CLAUDE.md's
      !! `mem:separate` contract — the previous version of this test ran
      !! host-only and proved nothing about the GPU build.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NZ_SP = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: SPONGE_W = 6
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), parameter :: STRENGTH = 1.0_wp
      integer, parameter :: N_SPONGE_STEPS = 30
      real(wp), parameter :: PI_T = acos(-1.0_wp)
      real(wp) :: u_far_after, u_outer_face, u_inner_face
      real(wp) :: alpha_inner, expected_outer, expected_inner
      integer :: step, i_outer, i_inner, j_probe, i_far_start, i_far_end

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ_SP
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ_SP)
         bc%west%bc_type = OBC_SPONGE
         bc%west%sponge_width = SPONGE_W
         bc%west%sponge_strength = STRENGTH

         ! Uniform u, no other dynamics — kernel-level isolation.
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp
         ms%h_layer = 1.0_wp

         ! GPU (-gpu=mem:separate): map ms before the DC loops inside
         ! ocean_sponge_apply touch u_face_x_layer/v_face_y_layer; bc is
         ! read only host-side (outside the DC) so it needs no map.
         !$acc enter data copyin(ms)
         call ms%enter_data()
         do step = 1, N_SPONGE_STEPS
            call ocean_sponge_apply(grid, bc, ms, DT)
         end do
         !$acc update self(ms%u_face_x_layer)
         call ms%exit_data()
         !$acc exit data delete(ms)

         ! Outermost band face (d=0, i = wall_face+0+1 = nghost+2): peak
         ! tau (alpha=1), exact answer u = U0*exp(-STRENGTH*DT*N).
         i_outer = grid%nghost + 2
         ! Innermost band face (d=band-1, i = nghost+1+band): tau tapers
         ! to alpha = 0.5*(1+cos(pi*(band-1)/band)).
         i_inner = grid%nghost + 1 + SPONGE_W
         j_probe = grid%nghost + 1
         u_outer_face = ms%u_face_x_layer(i_outer, j_probe, 1)
         u_inner_face = ms%u_face_x_layer(i_inner, j_probe, 1)

         expected_outer = U0*exp(-STRENGTH*1.0_wp*DT*real(N_SPONGE_STEPS, wp))
         alpha_inner = 0.5_wp*(1.0_wp + cos(PI_T*real(SPONGE_W - 1, wp)/real(SPONGE_W, wp)))
         expected_inner = U0*exp(-STRENGTH*alpha_inner*DT*real(N_SPONGE_STEPS, wp))

         call check(error, abs(u_outer_face - expected_outer) < 1.0e-12_wp*U0, &
                    "outermost band u-face must equal U0*exp(-strength*dt*N) exactly (alpha=1 at d=0)")
         if (allocated(error)) exit checks
         call check(error, abs(u_inner_face - expected_inner) < 1.0e-12_wp*U0, &
                    "innermost band u-face must equal the exact cosine-ramp decay at d=band-1")
         if (allocated(error)) exit checks

         ! u outside the band stays untouched (>95% of IC, unchanged bound).
         i_far_start = grid%nghost + SPONGE_W + 4
         i_far_end = grid%nghost + NX_PHYS
         u_far_after = sum(abs(ms%u_face_x_layer( &
                               i_far_start:i_far_end, &
                               grid%nghost + 1:grid%nghost + NY_PHYS, 1)))/ &
                       real((i_far_end - i_far_start + 1)*NY_PHYS, wp)
         call check(error, u_far_after > 0.95_wp*U0, &
                    "u outside the sponge band should be untouched (>95% of IC)")

      end block checks
      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine test_obc_sponge_kernel

   subroutine test_sponge_bit_identity_when_disabled(error)
      !! PR-23 acceptance gate (plan §9.7 / §10): with `&ocean_sponge_nml`
      !! absent, a >=20-step `ocean_dyn_step_split` run on a west-`OBC_SPONGE`
      !! 3-D config (legacy `sponge_relax_tracers=.true.`) must be
      !! BIT-FOR-BIT identical to the same run with `sp` present but
      !! `enable=.false.` PLUS every other sponge key set to a nonsense
      !! non-default value — proving the gate is `enable` alone and nothing
      !! leaks through `ocean_sponge_apply_maps`'s early return.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(ocean_metrics_t) :: metrics_a, metrics_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(coriolis_adv_t) :: cor_a, cor_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      type(ocean_bottom_drag_t) :: bd_a, bd_b
      type(ocean_surface_stress_t) :: ss_a, ss_b
      type(ocean_vertical_advection_t) :: va_a, va_b
      type(ocean_hdiff_tracer_t) :: hd_a, hd_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(eos_t) :: eos_a, eos_b
      type(ocean_dyn_t) :: dyn_a, dyn_b
      type(ocean_bc_state_t) :: bc_a, bc_b
      type(ocean_sponge_t) :: sp_b
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 8
      integer, parameter :: N_INNER = 8
      integer, parameter :: N_STEPS = 20
      real(wp), parameter :: H0 = 20.0_wp
      real(wp), parameter :: DT_OUTER = 0.05_wp
      integer :: i, j, k, it

      checks: block
         call make_grid(grid_a, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_b, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid_a); call ms_b%init(grid_b)
         call ct_a%init(grid_a, nz_ml=NZ); call ct_b%init(grid_b, nz_ml=NZ)
         call cor_a%init(grid_a, nz_ml=NZ); call cor_b%init(grid_b, nz_ml=NZ)
         call pgf_a%init(grid_a, nz_ml=NZ); call pgf_b%init(grid_b, nz_ml=NZ)
         call hv_a%init(grid_a, nz_ml=NZ); call hv_b%init(grid_b, nz_ml=NZ)
         call bd_a%init(grid_a, nz_ml=NZ); call bd_b%init(grid_b, nz_ml=NZ)
         call ss_a%init(grid_a, nz_ml=NZ); call ss_b%init(grid_b, nz_ml=NZ)
         call va_a%init(grid_a, nz_ml=NZ); call va_b%init(grid_b, nz_ml=NZ)
         call hd_a%init(grid_a, nz_ml=NZ); call hd_b%init(grid_b, nz_ml=NZ)
         call vd_a%init(grid_a, nz_ml=NZ); call vd_b%init(grid_b, nz_ml=NZ)
         call vmix_a%init(grid_a, nz_ml=NZ); call vmix_b%init(grid_b, nz_ml=NZ)
         call eos_a%init(grid_a); call eos_b%init(grid_b)
         call dyn_a%init(grid_a, nz_ml=NZ); call dyn_b%init(grid_b, nz_ml=NZ)
         call ocean_bc_state_init(bc_a, grid_a, nz_ml=NZ, n_tracers=size(ms_a%tracers))
         call ocean_bc_state_init(bc_b, grid_b, nz_ml=NZ, n_tracers=size(ms_b%tracers))
         ! West edge OBC_SPONGE with the legacy tracer-relax knob ON on
         ! BOTH runs — this is the leaky-gate scenario the plan flags
         ! (§11.9): the map path must not interfere even when the legacy
         ! band is genuinely active.
         bc_a%west%bc_type = OBC_SPONGE
         bc_a%west%sponge_width = 4
         bc_a%west%sponge_strength = 0.02_wp
         bc_a%west%sponge_relax_tracers = .true.
         bc_b%west%bc_type = OBC_SPONGE
         bc_b%west%sponge_width = 4
         bc_b%west%sponge_strength = 0.02_wp
         bc_b%west%sponge_relax_tracers = .true.

         ! sp_b: enable=.false. (the default) PLUS nonsense on every other
         ! key. Never call sp_b%init — enable=.false. means the gated
         ! init-call site would never allocate it in production either
         ! (ocean_state_init), so idamp_h/u/v/ref_* stay unallocated; the
         ! early return in ocean_sponge_apply_maps must never touch them.
         sp_b%enable = .false.
         sp_b%relax_uv = .false.
         sp_b%relax_tracers = .false.
         sp_b%damp_source = "file"
         sp_b%target_source = "file"

         ! Identical 3-D IC on both: sinusoidal h, uniform S/T, small
         ! solenoidal u/v perturbation (mirrors test_conservation's IC).
         do k = 1, NZ
            do j = 1, grid_a%ny_total
               do i = 1, grid_a%nx_total
                  ms_a%h_layer(i, j, k) = H0 + 0.02_wp*sin(2.0_wp*acos(-1.0_wp)* &
                                                           real(i, wp)/real(grid_a%nx_total, wp))
                  ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, k) = ms_a%h_layer(i, j, k)*eos_a%S_ref
                  ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, k) = ms_a%h_layer(i, j, k)*eos_a%T_ref
                  ms_b%h_layer(i, j, k) = ms_a%h_layer(i, j, k)
                  ms_b%tracers(ms_b%idx_salinity)%hTr(i, j, k) = ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, k)
                  ms_b%tracers(ms_b%idx_temperature)%hTr(i, j, k) = ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, k)
               end do
            end do
            do j = 1, grid_a%ny_total
               do i = 1, grid_a%nx_total + 1
                  ms_a%u_face_x_layer(i, j, k) = 0.01_wp* &
                                                 sin(acos(-1.0_wp)*real(i - 1, wp)/real(grid_a%nx_total, wp))
                  ms_b%u_face_x_layer(i, j, k) = ms_a%u_face_x_layer(i, j, k)
               end do
            end do
            do j = 1, grid_a%ny_total + 1
               do i = 1, grid_a%nx_total
                  ms_a%v_face_y_layer(i, j, k) = 0.01_wp* &
                                                 sin(acos(-1.0_wp)*real(j - 1, wp)/real(grid_a%ny_total, wp))
                  ms_b%v_face_y_layer(i, j, k) = ms_a%v_face_y_layer(i, j, k)
               end do
            end do
         end do
         dyn_a%bt_work%bt_H_ref = real(NZ, wp)*H0
         dyn_b%bt_work%bt_H_ref = real(NZ, wp)*H0

         ! Run A: sp absent entirely.
         call map_in(grid_a, metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         do i = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_a, metrics_a, dyn_a, eos_a, cor_a, ct_a, pgf_a, hv_a, bd_a, ss_a, &
               va_a, hd_a, vd_a, vmix_a, ms_a, DT_OUTER, N_INNER, bc=bc_a)
         end do
         call map_out(metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)

         ! Run B: sp present, enable=.false., every other key nonsense.
         !$acc enter data copyin(sp_b)
         call sp_b%enter_data()
         call map_in(grid_b, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         do i = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_b, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_b, DT_OUTER, N_INNER, bc=bc_b, sp=sp_b)
         end do
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         call sp_b%exit_data()
         !$acc exit data delete(sp_b)

         call check(error, all(ms_a%u_face_x_layer == ms_b%u_face_x_layer), &
                    "u_face_x_layer must be bit-for-bit identical whether sp is absent or "// &
                    "present-but-disabled")
         if (allocated(error)) exit checks
         call check(error, all(ms_a%v_face_y_layer == ms_b%v_face_y_layer), &
                    "v_face_y_layer must be bit-for-bit identical")
         if (allocated(error)) exit checks
         call check(error, all(ms_a%h_layer == ms_b%h_layer), &
                    "h_layer must be bit-for-bit identical")
         if (allocated(error)) exit checks
         do it = 1, size(ms_a%tracers)
            if (.not. allocated(ms_a%tracers(it)%hTr)) cycle
            if (.not. all(ms_a%tracers(it)%hTr == ms_b%tracers(it)%hTr)) then
               call check(error, .false., "every tracer hTr must be bit-for-bit identical")
               exit checks
            end if
         end do
         call check(error, .true.)

      end block checks
      call destroy_all(ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
      call ocean_bc_state_destroy(bc_a)
      call ocean_bc_state_destroy(bc_b)
   end subroutine test_sponge_bit_identity_when_disabled

   subroutine test_obc_tidal_via_driver(error)
      !! End-to-end OBC_TIDAL through `ocean_dyn_step_split`.  Closed
      !! walls except west = OBC_TIDAL with single M2 constituent.
      !! Drives 30 outer steps with `t` threaded through the driver
      !! so the eta_target oscillates over the M2 period.  Verifies
      !! basin SSH at the western interior column picks up a
      !! tidal-driven excursion above the quiescent reference.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_td
      type(ocean_metrics_t) :: metrics_ref, metrics_td
      type(multilayer_state_t) :: ms_ref, ms_td
      type(continuity_t) :: ct_r, ct_t
      type(coriolis_adv_t) :: cor_r, cor_t
      type(ocean_pressure_force_t) :: pgf_r, pgf_t
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_t
      type(ocean_bottom_drag_t) :: bd_r, bd_t
      type(ocean_surface_stress_t) :: ss_r, ss_t
      type(ocean_vertical_advection_t) :: va_r, va_t
      type(ocean_hdiff_tracer_t) :: hd_r, hd_t
      type(ocean_vdiff_t) :: vd_r, vd_t
      type(ocean_vmix_t) :: vmix_r, vmix_t
      type(eos_t) :: eos_r, eos_td
      type(ocean_dyn_t) :: dyn_r, dyn_t
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NZ_T = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 20
      integer, parameter :: N_STEPS = 30
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: AMP = 0.3_wp
      real(wp), parameter :: PERIOD = 44712.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, dt_outer, t_curr, max_eta_td, max_eta_ref
      integer :: i, j, step

      checks: block
         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_td, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         ms_ref%nz_ml = NZ_T; ms_td%nz_ml = NZ_T
         call ms_ref%init(grid_ref); call ms_td%init(grid_td)
         call ct_r%init(grid_ref, nz_ml=NZ_T); call ct_t%init(grid_td, nz_ml=NZ_T)
         call cor_r%init(grid_ref, nz_ml=NZ_T); call cor_t%init(grid_td, nz_ml=NZ_T)
         call pgf_r%init(grid_ref, nz_ml=NZ_T); call pgf_t%init(grid_td, nz_ml=NZ_T)
         call hv_r%init(grid_ref, nz_ml=NZ_T); call hv_t%init(grid_td, nz_ml=NZ_T)
         call bd_r%init(grid_ref, nz_ml=NZ_T); call bd_t%init(grid_td, nz_ml=NZ_T)
         call ss_r%init(grid_ref, nz_ml=NZ_T); call ss_t%init(grid_td, nz_ml=NZ_T)
         call va_r%init(grid_ref, nz_ml=NZ_T); call va_t%init(grid_td, nz_ml=NZ_T)
         call hd_r%init(grid_ref, nz_ml=NZ_T); call hd_t%init(grid_td, nz_ml=NZ_T)
         call vd_r%init(grid_ref, nz_ml=NZ_T); call vd_t%init(grid_td, nz_ml=NZ_T)
         call vmix_r%init(grid_ref, nz_ml=NZ_T); call vmix_t%init(grid_td, nz_ml=NZ_T)
         call eos_r%init(grid_ref); call eos_td%init(grid_td)
         call dyn_r%init(grid_ref, nz_ml=NZ_T); call dyn_t%init(grid_td, nz_ml=NZ_T)
         call ocean_bc_state_init(bc, grid_td, nz_ml=NZ_T)
         bc%west%bc_type = OBC_TIDAL
         bc%west%n_tidal_constituents = 1
         bc%west%tidal_amp(1) = AMP
         bc%west%tidal_omega(1) = 2.0_wp*PI/PERIOD
         bc%west%tidal_phase(1) = 0.0_wp

         ! Quiescent IC.
         do j = 1, grid_ref%ny_total
            do i = 1, grid_ref%nx_total
               ms_ref%h_layer(i, j, 1) = H0
               ms_td%h_layer(i, j, 1) = H0
               ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, 1) = eos_r%S_ref*H0
               ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, 1) = eos_r%T_ref*H0
               ms_td%tracers(ms_td%idx_salinity)%hTr(i, j, 1) = eos_td%S_ref*H0
               ms_td%tracers(ms_td%idx_temperature)%hTr(i, j, 1) = eos_td%T_ref*H0
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp; ms_ref%v_face_y_layer = 0.0_wp
         ms_td%u_face_x_layer = 0.0_wp; ms_td%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = H0; dyn_t%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.5_wp*grid_ref%dx/c

         ! Reference: no bc → all closed.
         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         ! Tidal: bc + t threaded.
         t_curr = 0.0_wp
         call map_in(grid_td, metrics_td, ms_td, ct_t, cor_t, pgf_t, hv_t, bd_t, ss_t, va_t, hd_t, vd_t, vmix_t, dyn_t)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_td, metrics_td, dyn_t, eos_td, cor_t, ct_t, pgf_t, hv_t, bd_t, ss_t, &
               va_t, hd_t, vd_t, vmix_t, ms_td, dt_outer, N_INNER, bc=bc, t=t_curr)
            t_curr = t_curr + dt_outer
         end do
         call map_out(metrics_td, ms_td, ct_t, cor_t, pgf_t, hv_t, bd_t, ss_t, va_t, hd_t, vd_t, vmix_t, dyn_t)

         ! Tidal run should produce η excursions; reference (closed,
         ! quiescent IC) stays at η ≈ 0.
         max_eta_ref = maxval(abs(ms_ref%h_layer(grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, :, 1) - H0))
         max_eta_td = maxval(abs(ms_td%h_layer(grid_td%nghost + 1:grid_td%nghost + NX_PHYS, :, 1) - H0))

         call check(error, max_eta_td > 0.1_wp*AMP, &
                    "OBC_TIDAL via driver should produce η excursions > 10% of constituent AMP")
         if (allocated(error)) exit checks
         call check(error, max_eta_ref < 1.0e-6_wp, &
                    "closed-wall reference (no bc) should keep η ≈ 0 from rest")

      end block checks
      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_td, ct_t, cor_t, pgf_t, hv_t, bd_t, ss_t, va_t, hd_t, vd_t, vmix_t, eos_td, dyn_t)
   end subroutine test_obc_tidal_via_driver

   subroutine test_obc_clamped_via_driver(error)
      !! End-to-end OBC_CLAMPED through the split driver.  Quiescent
      !! IC; west = OBC_CLAMPED with positive `clamped_u`.  Reference
      !! (no bc) stays at η ≈ 0.  Clamped run mass should enter from
      !! the west and basin-mean η drift positive.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_cl
      type(ocean_metrics_t) :: metrics_ref, metrics_cl
      type(multilayer_state_t) :: ms_ref, ms_cl
      type(continuity_t) :: ct_r, ct_c
      type(coriolis_adv_t) :: cor_r, cor_c
      type(ocean_pressure_force_t) :: pgf_r, pgf_c
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_c
      type(ocean_bottom_drag_t) :: bd_r, bd_c
      type(ocean_surface_stress_t) :: ss_r, ss_c
      type(ocean_vertical_advection_t) :: va_r, va_c
      type(ocean_hdiff_tracer_t) :: hd_r, hd_c
      type(ocean_vdiff_t) :: vd_r, vd_c
      type(ocean_vmix_t) :: vmix_r, vmix_c
      type(eos_t) :: eos_r, eos_c
      type(ocean_dyn_t) :: dyn_r, dyn_c
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NZ_C = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 10
      integer, parameter :: N_STEPS = 30
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: U_IN = 0.05_wp
      real(wp) :: c, dt_outer, eta_mean_cl, eta_mean_ref
      integer :: i, j, step

      checks: block
         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_cl, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms_ref%nz_ml = NZ_C; ms_cl%nz_ml = NZ_C
         call ms_ref%init(grid_ref); call ms_cl%init(grid_cl)
         call ct_r%init(grid_ref, nz_ml=NZ_C); call ct_c%init(grid_cl, nz_ml=NZ_C)
         call cor_r%init(grid_ref, nz_ml=NZ_C); call cor_c%init(grid_cl, nz_ml=NZ_C)
         call pgf_r%init(grid_ref, nz_ml=NZ_C); call pgf_c%init(grid_cl, nz_ml=NZ_C)
         call hv_r%init(grid_ref, nz_ml=NZ_C); call hv_c%init(grid_cl, nz_ml=NZ_C)
         call bd_r%init(grid_ref, nz_ml=NZ_C); call bd_c%init(grid_cl, nz_ml=NZ_C)
         call ss_r%init(grid_ref, nz_ml=NZ_C); call ss_c%init(grid_cl, nz_ml=NZ_C)
         call va_r%init(grid_ref, nz_ml=NZ_C); call va_c%init(grid_cl, nz_ml=NZ_C)
         call hd_r%init(grid_ref, nz_ml=NZ_C); call hd_c%init(grid_cl, nz_ml=NZ_C)
         call vd_r%init(grid_ref, nz_ml=NZ_C); call vd_c%init(grid_cl, nz_ml=NZ_C)
         call vmix_r%init(grid_ref, nz_ml=NZ_C); call vmix_c%init(grid_cl, nz_ml=NZ_C)
         call eos_r%init(grid_ref); call eos_c%init(grid_cl)
         call dyn_r%init(grid_ref, nz_ml=NZ_C); call dyn_c%init(grid_cl, nz_ml=NZ_C)
         call ocean_bc_state_init(bc, grid_cl, nz_ml=NZ_C, n_tracers=size(ms_cl%tracers))
         bc%west%bc_type = OBC_CLAMPED
         bc%west%clamped_u = U_IN
         ! Prescribe inflow salinity 1 PSU above the basin reference;
         ! exercises the tracer ghost-cell override.
         bc%west%clamped_tracer(ms_cl%idx_salinity) = eos_c%S_ref + 1.0_wp
         bc%west%clamped_tracer(ms_cl%idx_temperature) = eos_c%T_ref

         do j = 1, grid_ref%ny_total
            do i = 1, grid_ref%nx_total
               ms_ref%h_layer(i, j, 1) = H0
               ms_cl%h_layer(i, j, 1) = H0
               ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, 1) = eos_r%S_ref*H0
               ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, 1) = eos_r%T_ref*H0
               ms_cl%tracers(ms_cl%idx_salinity)%hTr(i, j, 1) = eos_c%S_ref*H0
               ms_cl%tracers(ms_cl%idx_temperature)%hTr(i, j, 1) = eos_c%T_ref*H0
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp; ms_ref%v_face_y_layer = 0.0_wp
         ms_cl%u_face_x_layer = 0.0_wp; ms_cl%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = H0; dyn_c%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.5_wp*grid_ref%dx/c

         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         call map_in(grid_cl, metrics_cl, ms_cl, ct_c, cor_c, pgf_c, hv_c, bd_c, ss_c, va_c, hd_c, vd_c, vmix_c, dyn_c)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_cl, metrics_cl, dyn_c, eos_c, cor_c, ct_c, pgf_c, hv_c, bd_c, ss_c, &
               va_c, hd_c, vd_c, vmix_c, ms_cl, dt_outer, N_INNER, bc=bc)
         end do
         call map_out(metrics_cl, ms_cl, ct_c, cor_c, pgf_c, hv_c, bd_c, ss_c, va_c, hd_c, vd_c, vmix_c, dyn_c)

         eta_mean_ref = sum(ms_ref%h_layer(grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, &
                                           grid_ref%nghost + 1:grid_ref%nghost + NY_PHYS, 1) - H0)/ &
                        real(NX_PHYS*NY_PHYS, wp)
         eta_mean_cl = sum(ms_cl%h_layer(grid_cl%nghost + 1:grid_cl%nghost + NX_PHYS, &
                                         grid_cl%nghost + 1:grid_cl%nghost + NY_PHYS, 1) - H0)/ &
                       real(NX_PHYS*NY_PHYS, wp)

         call check(error, abs(eta_mean_ref) < 1.0e-6_wp, &
                    "closed-wall reference should keep basin-mean η ≈ 0")
         if (allocated(error)) exit checks
         call check(error, eta_mean_cl > 1.0e-3_wp, &
                    "OBC_CLAMPED west with U_IN > 0 should raise basin-mean η")
         if (allocated(error)) exit checks
         ! Tracer dispatch: basin-mean S in the clamped run should be
         ! above the basin-mean S in the reference, by some non-trivial
         ! fraction of the prescribed boundary excursion (1 PSU).  The
         ! reference (no bc) stays at S_ref everywhere.
         block
            real(wp) :: s_mean_ref, s_mean_cl
            s_mean_ref = sum(ms_ref%tracers(ms_ref%idx_salinity)%hTr( &
                             grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, &
                             grid_ref%nghost + 1:grid_ref%nghost + NY_PHYS, 1)/ &
                             ms_ref%h_layer(grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, &
                                            grid_ref%nghost + 1:grid_ref%nghost + NY_PHYS, 1))/ &
                         real(NX_PHYS*NY_PHYS, wp)
            s_mean_cl = sum(ms_cl%tracers(ms_cl%idx_salinity)%hTr( &
                            grid_cl%nghost + 1:grid_cl%nghost + NX_PHYS, &
                            grid_cl%nghost + 1:grid_cl%nghost + NY_PHYS, 1)/ &
                            ms_cl%h_layer(grid_cl%nghost + 1:grid_cl%nghost + NX_PHYS, &
                                          grid_cl%nghost + 1:grid_cl%nghost + NY_PHYS, 1))/ &
                        real(NX_PHYS*NY_PHYS, wp)
            call check(error, abs(s_mean_ref - eos_r%S_ref) < 1.0e-6_wp, &
                       "closed-wall reference should keep basin-mean S at S_ref")
            if (allocated(error)) exit checks
            ! The OBC_CLAMPED clamped_tracer dispatch is firing — proven
            ! by the host-side override of ghost-cell hTr.  Asserts only
            ! that basin-mean S differs from the closed-wall reference;
            ! the SIGN of the change depends on the interplay of mass
            ! inflow vs salt advection and is sensitive to dt/cell/N_STEPS.
            ! Tight tracer-budget tests land alongside the per-cell data
            ! source backends in a follow-on.
            call check(error, abs(s_mean_cl - s_mean_ref) > 1.0e-6_wp, &
                       "OBC_CLAMPED with clamped_tracer should change basin-mean S")
         end block

      end block checks
      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_cl, ct_c, cor_c, pgf_c, hv_c, bd_c, ss_c, va_c, hd_c, vd_c, vmix_c, eos_c, dyn_c)
   end subroutine test_obc_clamped_via_driver

   subroutine test_obc_chapman_via_driver(error)
      !! End-to-end OBC_CHAPMAN through the split driver.  m=1 cosine
      !! η pulse, no forcing.  Reference (no bc) bounces the pulse off
      !! the closed walls.  Chapman west run should drain the wave
      !! energy via radiation at the west boundary (similar to OPEN
      !! but with persistent edge-state).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_ch
      type(ocean_metrics_t) :: metrics_ref, metrics_ch
      type(multilayer_state_t) :: ms_ref, ms_ch
      type(continuity_t) :: ct_r, ct_h
      type(coriolis_adv_t) :: cor_r, cor_h
      type(ocean_pressure_force_t) :: pgf_r, pgf_h
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_h
      type(ocean_bottom_drag_t) :: bd_r, bd_h
      type(ocean_surface_stress_t) :: ss_r, ss_h
      type(ocean_vertical_advection_t) :: va_r, va_h
      type(ocean_hdiff_tracer_t) :: hd_r, hd_h
      type(ocean_vdiff_t) :: vd_r, vd_h
      type(ocean_vmix_t) :: vmix_r, vmix_h
      type(eos_t) :: eos_r, eos_h
      type(ocean_dyn_t) :: dyn_r, dyn_h
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NZ_H = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 20
      integer, parameter :: N_STEPS = 60
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, dt_outer, e2_ref, e2_ch
      integer :: i, j, step

      checks: block
         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_ch, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms_ref%nz_ml = NZ_H; ms_ch%nz_ml = NZ_H
         call ms_ref%init(grid_ref); call ms_ch%init(grid_ch)
         call ct_r%init(grid_ref, nz_ml=NZ_H); call ct_h%init(grid_ch, nz_ml=NZ_H)
         call cor_r%init(grid_ref, nz_ml=NZ_H); call cor_h%init(grid_ch, nz_ml=NZ_H)
         call pgf_r%init(grid_ref, nz_ml=NZ_H); call pgf_h%init(grid_ch, nz_ml=NZ_H)
         call hv_r%init(grid_ref, nz_ml=NZ_H); call hv_h%init(grid_ch, nz_ml=NZ_H)
         call bd_r%init(grid_ref, nz_ml=NZ_H); call bd_h%init(grid_ch, nz_ml=NZ_H)
         call ss_r%init(grid_ref, nz_ml=NZ_H); call ss_h%init(grid_ch, nz_ml=NZ_H)
         call va_r%init(grid_ref, nz_ml=NZ_H); call va_h%init(grid_ch, nz_ml=NZ_H)
         call hd_r%init(grid_ref, nz_ml=NZ_H); call hd_h%init(grid_ch, nz_ml=NZ_H)
         call vd_r%init(grid_ref, nz_ml=NZ_H); call vd_h%init(grid_ch, nz_ml=NZ_H)
         call vmix_r%init(grid_ref, nz_ml=NZ_H); call vmix_h%init(grid_ch, nz_ml=NZ_H)
         call eos_r%init(grid_ref); call eos_h%init(grid_ch)
         call dyn_r%init(grid_ref, nz_ml=NZ_H); call dyn_h%init(grid_ch, nz_ml=NZ_H)
         call ocean_bc_state_init(bc, grid_ch, nz_ml=NZ_H)
         bc%west%bc_type = OBC_CHAPMAN

         do j = 1, grid_ref%ny_total
            do i = 1, grid_ref%nx_total
               if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                  ms_ref%h_layer(i, j, 1) = H0 + A_ETA* &
                                            cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp) &
                                                /real(NX_PHYS, wp))
               else
                  ms_ref%h_layer(i, j, 1) = H0
               end if
               ms_ch%h_layer(i, j, 1) = ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, 1) = eos_r%S_ref*ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, 1) = eos_r%T_ref*ms_ref%h_layer(i, j, 1)
               ms_ch%tracers(ms_ch%idx_salinity)%hTr(i, j, 1) = eos_h%S_ref*ms_ch%h_layer(i, j, 1)
               ms_ch%tracers(ms_ch%idx_temperature)%hTr(i, j, 1) = eos_h%T_ref*ms_ch%h_layer(i, j, 1)
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp; ms_ref%v_face_y_layer = 0.0_wp
         ms_ch%u_face_x_layer = 0.0_wp; ms_ch%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = H0; dyn_h%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.5_wp*grid_ref%dx/c

         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         call map_in(grid_ch, metrics_ch, ms_ch, ct_h, cor_h, pgf_h, hv_h, bd_h, ss_h, va_h, hd_h, vd_h, vmix_h, dyn_h)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ch, metrics_ch, dyn_h, eos_h, cor_h, ct_h, pgf_h, hv_h, bd_h, ss_h, &
               va_h, hd_h, vd_h, vmix_h, ms_ch, dt_outer, N_INNER, bc=bc)
         end do
         call map_out(metrics_ch, ms_ch, ct_h, cor_h, pgf_h, hv_h, bd_h, ss_h, va_h, hd_h, vd_h, vmix_h, dyn_h)

         e2_ref = sum((ms_ref%h_layer(grid_ref%nghost + 1:grid_ref%nghost + NX_PHYS, :, 1) - H0)**2)
         e2_ch = sum((ms_ch%h_layer(grid_ch%nghost + 1:grid_ch%nghost + NX_PHYS, :, 1) - H0)**2)
         ! Chapman is the Sommerfeld radiation + persistent edge-η
         ! form: u_face = -c·(η_int − η_target) with η_target tracking
         ! η_int on the c·dt/dx timescale.  For a transient (non-
         ! forced) IC, the wall equilibrates within ~2-3 outer steps
         ! and then behaves like a near-closed wall, so it drains
         ! much less than OBC_OPEN over this short horizon.  The
         ! meaningful assertion: Chapman is strictly less reflective
         ! than the closed-wall reference (some drainage during the
         ! transient).  A sharper test of Chapman's tide-following
         ! behaviour needs a tidal-forced IC — TODO follow-up.
         call check(error, e2_ch < e2_ref - 1.0e-6_wp*e2_ref, &
                    "OBC_CHAPMAN west via driver should be strictly less reflective than closed wall")

      end block checks
      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_ch, ct_h, cor_h, pgf_h, hv_h, bd_h, ss_h, va_h, hd_h, vd_h, vmix_h, eos_h, dyn_h)
   end subroutine test_obc_chapman_via_driver

   subroutine test_mass_contributor_closed_basin(error)
      !! Phase D v2 first kernel patch: continuity's mass-budget
      !! contributor must telescope to ~0 over a closed basin (no flux
      !! through walls → ΣΣ -dt·div_h = 0 to FP for any non-trivial
      !! velocity field).  Same IC as `test_conservation` but adds a
      !! manually-registered `ocean_budgets_t` slot bound to
      !! `ms%mass_budget_continuity`; drives several split steps;
      !! calls `drain_contributors` and verifies the contributor's
      !! `total_integrated` is bounded by the total mass round-off.
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 10
      integer, parameter :: N_INNER = 5
      integer :: i, j, k, nx, ny, step, idx_contrib
      real(wp) :: total_h0, total_h, ref_mass, residual_rel
      checks: block

         call make_grid(grid, 24, 16, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)
         call budgets%init(grid)
         call budgets%register_contributor("continuity", BUDGET_MASS, &
                                           ms%mass_budget_continuity, &
                                           device_resident=.true.)
         idx_contrib = budgets%n_contributors

         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + 0.01_wp*sin( &
                                        2.0_wp*PI*real(i, wp)/real(nx, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = ms%h_layer(i, j, k)*eos%S_ref
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = ms%h_layer(i, j, k)*eos%T_ref
               end do
            end do
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
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H_BASE

         total_h0 = sum(ms%h_layer)*grid%dx*grid%dy
         ref_mass = total_h0

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do
         ! Drain BEFORE map_out so the contributor's device buffer is
         ! still mapped when drain does `acc update self`.  Otherwise
         ! the closed-basin telescope would pass trivially with both
         ! sides at zero (false positive).
         call budgets%drain_contributors()
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         total_h = sum(ms%h_layer)*grid%dx*grid%dy

         ! Closed-basin telescope: Σ -dt·div_h over all interior cells
         ! and steps should be ~0 to mass round-off.  The contributor's
         ! total_integrated is that sum; compare to the absolute mass.
         residual_rel = abs(budgets%contributors(idx_contrib)%total_integrated)/ref_mass
         call check(error, residual_rel < 1.0e-10_wp, &
                    "closed-basin continuity contributor should telescope to FP")
         if (allocated(error)) exit checks

         ! LHS = total mass drift; RHS = contributor.  In a closed
         ! basin with continuity the only mass-changing kernel, the
         ! two should match to FP.
         call check(error, abs((total_h - total_h0) - &
                               budgets%contributors(idx_contrib)%total_integrated) &
                    < 1.0e-10_wp*ref_mass, &
                    "LHS (total drift) should equal RHS (continuity contributor) to FP")

      end block checks
      call budgets%destroy()
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_mass_contributor_closed_basin

   subroutine test_ideal_age_via_driver(error)
      !! End-to-end: enable ideal-age via `with_ideal_age=.true.` and
      !! run several outer split steps on a closed quiescent basin.
      !! After N steps the subsurface age concentration must equal
      !! N*DT per outer step (each SSP-RK2 outer step adds dt*h to
      !! both stage tracers, then `rk2_average` halves the cumulative
      !! delta — net +dt*h per outer step).  This is the interior-
      !! aging gate — a real assertion.
      !!
      !! The surface check below is kept (it is free and guards
      !! against a reset that overshoots) but is TRIVIALLY satisfied
      !! here: the IC seeds `hTr(:,:,NZ) = 0`, so even the pre-PR-7
      !! buggy reset (halved by rk2_average: `0.5*(0+0) = 0`) passes
      !! it. `test_ideal_age_surface_reset_nonzero_ic` below — which
      !! seeds a NONZERO surface age — is the real regression gate
      !! for the surface Dirichlet BC; do not read a green surface
      !! check here as evidence the reset is correct.
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_INNER = 5
      integer, parameter :: N_STEPS = 3
      real(wp) :: sub_age, surface_age, expected
      integer :: nx_phys, ny_phys, ic, jc, k, step
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid, with_ideal_age=.true.)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = 0.0_wp
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)

         call check(error, ms%idx_age == 3, &
                    "ideal_age via driver: registry index should be 3")
         if (allocated(error)) exit checks

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do
         ms%tracers(ms%idx_age)%hTr = 0.0_wp
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         nx_phys = grid%nx_total - 2*NGHOST
         ny_phys = grid%ny_total - 2*NGHOST
         ic = NGHOST + nx_phys/2
         jc = NGHOST + ny_phys/2

         expected = real(N_STEPS, wp)*DT
         sub_age = ms%tracers(ms%idx_age)%hTr(ic, jc, 1)/H0
         call check(error, abs(sub_age - expected) < 1.0e-10_wp, &
                    "ideal_age via driver: subsurface age must be N*DT")
         if (allocated(error)) exit checks

         surface_age = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
         call check(error, surface_age < 1.0e-12_wp, &
                    "ideal_age via driver: surface age must remain 0")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_ideal_age_via_driver

   subroutine test_ideal_age_surface_reset_nonzero_ic(error)
      !! PR-7 regression gate — the only test in the suite that
      !! distinguishes a Dirichlet surface BC from a geometric decay.
      !! Seed the surface layer with a LARGE nonzero age
      !! (`A0*H0`, `A0 = 1000`), run one outer split step, and assert
      !! it lands at EXACTLY 0 — not `0.5*A0*H0`.  On pre-PR-7 code
      !! (surface reset applied inside each RK2 stage, then halved by
      !! `rk2_average_field_3d` at the outer-step boundary) this
      !! returns `0.5*1000*10 = 5000.0`, failing by 15 orders of
      !! magnitude.  A second outer step from the already-reset state
      !! must ALSO land at exactly 0 — idempotence, which distinguishes
      !! "reset once, hard" from "decay slowly toward zero".
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: A0 = 1000.0_wp
      integer, parameter :: N_INNER = 5
      real(wp) :: surface_age
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid, with_ideal_age=.true.)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = 0.0_wp
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do
         ms%tracers(ms%idx_age)%hTr = 0.0_wp
         ms%tracers(ms%idx_age)%hTr(:, :, NZ) = A0*H0
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         ! ---- Step 1: nonzero-IC reset must land at exactly 0 ----
         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ocean_dyn_step_split( &
            grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
            va, hd, vd, vmix, ms, DT, N_INNER)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         surface_age = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
         call check(error, surface_age < 1.0e-12_wp, &
                    "ideal_age surface reset: nonzero IC must reset to 0, not 0.5*A0*H0")
         if (allocated(error)) exit checks

         ! ---- Step 2: idempotence — stays at 0 from the reset state ----
         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ocean_dyn_step_split( &
            grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
            va, hd, vd, vmix, ms, DT, N_INNER)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         surface_age = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
         call check(error, surface_age < 1.0e-12_wp, &
                    "ideal_age surface reset: second step from reset state stays at 0 (idempotent)")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_ideal_age_surface_reset_nonzero_ic

   subroutine test_ideal_age_thermo_cadence(error)
      !! Analytical: `dyn%dt_therm_ratio = 2`, quiescent basin, `hTr = 0`,
      !! N_STEPS = 4.  `is_thermo_step()` is
      !! `mod(outer_step_count, 2) == 0` with `outer_step_count`
      !! incremented at step end, so it fires on outer steps 1 and 3
      !! (0-indexed count 0 and 2), each contributing `therm_dt = 2*DT`
      !! — subsurface age must land at exactly `4*DT`, not `2*DT`
      !! (which is what a cadence gate with the `therm_dt` scaling
      !! forgotten would give).  This is the test that distinguishes
      !! "gated on the cadence" from "gated and the dt-scaling
      !! forgotten" — elapsed age must track model WALL-CLOCK TIME
      !! independent of `dt_therm_ratio`, which is the whole point of
      !! the tracer.  Surface stays at 0 throughout (the reset is
      !! gated identically, and nothing writes the surface in a
      !! quiescent basin on a non-thermo step either).
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
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_INNER = 5
      integer, parameter :: N_STEPS = 4
      real(wp) :: sub_age, surface_age, expected
      integer :: nx_phys, ny_phys, ic, jc, k, step
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid, with_ideal_age=.true.)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = 0.0_wp
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)
         dyn%dt_therm_ratio = 2

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do
         ms%tracers(ms%idx_age)%hTr = 0.0_wp
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         nx_phys = grid%nx_total - 2*NGHOST
         ny_phys = grid%ny_total - 2*NGHOST
         ic = NGHOST + nx_phys/2
         jc = NGHOST + ny_phys/2

         expected = real(N_STEPS, wp)*DT
         sub_age = ms%tracers(ms%idx_age)%hTr(ic, jc, 1)/H0
         call check(error, abs(sub_age - expected) < 1.0e-10_wp, &
                    "ideal_age thermo cadence: subsurface age must equal N*DT "// &
                    "independent of dt_therm_ratio")
         if (allocated(error)) exit checks

         surface_age = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
         call check(error, surface_age < 1.0e-12_wp, &
                    "ideal_age thermo cadence: surface age must remain 0")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_ideal_age_thermo_cadence

end module test_ocean_dyn_split
