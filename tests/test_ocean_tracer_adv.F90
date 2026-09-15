!! Unit tests for the per-layer PPM tracer-advection kernel
!! (rdb_continuity::tracer_advect).  Iterates over the
!! tracer registry on `multilayer_state_t%tracers(:)` and
!! advects each enabled tracer using the same mass fluxes that
!! continuity-PPM produced — the consistency-with-continuity (CWC)
!! property — so a uniform tracer field stays uniform under any
!! divergent flow.
!!
!! Cases:
!!   * Zero-velocity preservation — non-trivial Tr, zero u, v.
!!     mass_flux is zero everywhere, so hTr cannot change.  Bit-
!!     for-bit guard.
!!   * CWC constancy preservation — uniform Tr = Tr_const, smooth
!!     non-zero velocity, run continuity + tracer advect + apply
!!     for several steps.  Tr_new = hTr_new/h_new must still equal
!!     Tr_const everywhere to round-off.  The defining test of the
!!     advection scheme.
!!   * Tracer mass conservation — non-trivial Tr, closed-basin
!!     velocity, 50 steps.  Total tracer mass `sum(hTr) * dx*dy`
!!     must drift by less than 1e-10 relative.
!!   * `do_horizontal_advection` toggle respected — disable
!!     advection on temperature, leave salinity enabled, step.
!!     Salinity's hTr changes; temperature's hTr stays at its
!!     initial value.
module test_ocean_tracer_adv
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, &
                             continuity_compute_fluxes, &
                             continuity_apply_fluxes, &
                             tracer_advect, &
                             continuity_tracer_step_split
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_tracer_adv_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_tracer_adv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zero_velocity_preserves_tracer", &
                               test_zero_velocity), &
                  new_unittest("cwc_constancy_preservation", test_cwc_constancy), &
                  new_unittest("tracer_mass_conservation", test_mass_conservation), &
                  new_unittest("do_horizontal_advection_toggle", &
                               test_toggle_advection), &
                  new_unittest("split_cwc_constancy_preservation", &
                               test_split_cwc_constancy), &
                  new_unittest("split_tracer_mass_conservation", &
                               test_split_mass_conservation) &
                  ]
   end subroutine collect_ocean_tracer_adv_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in

   subroutine map_out(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine run_step(grid, metrics, ct, ms, dt)
      !! One coupled continuity + tracer step.
      !!   compute_fluxes -> tracer_advect (uses h_old + mass_flux)
      !!   -> continuity_apply (updates h)
      !! Both updates consume the SAME mass_flux, which is the CWC
      !! property.  Order between tracer_advect and continuity_apply
      !! doesn't matter algebraically (tracer_advect reads h_old and
      !! mass_flux, neither changes between the two applies).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      call continuity_compute_fluxes(grid, metrics, ct, ms)
      call tracer_advect(grid, metrics, ct, ms, dt)
      call continuity_apply_fluxes(ms, dt)
   end subroutine run_step

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_zero_velocity(error)
      !! Zero u, v means mass_flux = 0, which means no tracer mass
      !! moves.  hTr must stay bit-for-bit equal to its IC.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 12.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: max_dS, max_dT
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = S0*H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0

         call map_in(ms, ct)
         call run_step(grid, metrics, ct, ms, DT)
         call map_out(ms, ct)

         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - S0*H0))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - T0*H0))
         call check(error, max_dS < 1.0e-14_wp, &
                    "salinity hTr drifted under zero velocity")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-14_wp, &
                    "temperature hTr drifted under zero velocity")

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_zero_velocity

   subroutine test_cwc_constancy(error)
      !! Uniform Tr field + non-trivial velocity: after coupled
      !! continuity + tracer advect + continuity-apply, Tr_new =
      !! hTr_new / h_new must still equal the initial constant.
      !! This is the defining CWC test for tracer advection.  Any
      !! mass_flux mismatch (e.g. tracer reads a different flux than
      !! continuity used to update h) breaks it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: S_CONST = 34.5_wp
      real(wp), parameter :: T_CONST = 18.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 30
      integer :: i, j, k, nx, ny, step
      real(wp) :: max_dS, max_dT, h_min
      real(wp), allocatable :: S_new(:, :, :), T_new(:, :, :)
      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H_BASE
         ms%tracers(ms%idx_salinity)%hTr = S_CONST*H_BASE
         ms%tracers(ms%idx_temperature)%hTr = T_CONST*H_BASE
         ! Wall-vanishing velocity field per layer
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
         end do

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call run_step(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         h_min = minval(ms%h_layer)
         call check(error, h_min > 0.0_wp, &
                    "h_layer went non-positive — kernel unstable, can't probe CWC")
         if (allocated(error)) exit checks

         allocate (S_new(nx, ny, NZ), T_new(nx, ny, NZ))
         S_new = ms%tracers(ms%idx_salinity)%hTr/ms%h_layer
         T_new = ms%tracers(ms%idx_temperature)%hTr/ms%h_layer
         max_dS = maxval(abs(S_new - S_CONST))
         max_dT = maxval(abs(T_new - T_CONST))
         deallocate (S_new, T_new)

         call check(error, max_dS < 1.0e-10_wp, &
                    "CWC broken: uniform S deviated under non-zero flow")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-10_wp, &
                    "CWC broken: uniform T deviated under non-zero flow")

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_cwc_constancy

   subroutine test_mass_conservation(error)
      !! Non-trivial S, T fields + wall-vanishing velocity, 50
      !! coupled steps.  Total tracer mass sum(hTr)*dx*dy must
      !! drift by less than 1e-10 (round-off floor; walls force
      !! zero tracer flux through them since mass_flux at walls is
      !! zero by the continuity kernel).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: S_BASE = 35.0_wp
      real(wp), parameter :: T_BASE = 15.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 50
      integer :: i, j, k, nx, ny, step
      real(wp) :: total_S0, total_T0, total_S, total_T, drift_S, drift_T
      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H_BASE
         ! Spatially varying S, T
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = H_BASE*(S_BASE + &
                                                                     0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                                                                     cos(2.0_wp*PI*real(j, wp)/real(ny, wp)))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = H_BASE*(T_BASE + &
                                                                        0.3_wp*cos(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                                                                        sin(2.0_wp*PI*real(j, wp)/real(ny, wp)))
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

         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T0 = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call run_step(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         total_S = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         drift_S = abs(total_S - total_S0)/abs(total_S0)
         drift_T = abs(total_T - total_T0)/abs(total_T0)

         call check(error, drift_S < 1.0e-10_wp, &
                    "total salinity mass drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, drift_T < 1.0e-10_wp, &
                    "total temperature mass drift > 1e-10")

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_mass_conservation

   subroutine test_toggle_advection(error)
      !! Disable advection on temperature, leave salinity enabled.
      !! After a step with non-zero velocity, salinity changes
      !! (because flow is non-trivial); temperature's hTr stays
      !! bit-for-bit at its initial value.  Confirms the registry's
      !! per-tracer `do_horizontal_advection` flag is honoured —
      !! prerequisite for the extendable design (passive
      !! diagnostics, externally-forced tracers).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      integer :: i, j, k, nx, ny
      real(wp) :: max_dT, max_dS
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      checks: block

         call make_grid(grid, 16, 8, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ! Spatially varying S so PPM has gradients to operate on
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = H0*(35.0_wp + &
                                                                 0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = H0*(15.0_wp + &
                                                                    0.3_wp*cos(2.0_wp*PI*real(i, wp)/real(nx, wp)))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp))
               end do
            end do
            ms%v_face_y_layer(:, :, k) = 0.0_wp
         end do

         ! Snapshot initial hT, hS for the post-step comparison
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         ! Disable temperature advection
         ms%tracers(ms%idx_temperature)%do_horizontal_advection = .false.

         call map_in(ms, ct)
         call run_step(grid, metrics, ct, ms, DT)
         call map_out(ms, ct)

         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))

         call check(error, max_dT < 1.0e-14_wp, &
                    "temperature hTr changed despite do_horizontal_advection = .false.")
         if (allocated(error)) exit checks
         call check(error, max_dS > 1.0e-6_wp, &
                    "salinity hTr did not change — flow setup is insufficient")

      end block checks
      deallocate (hT_ic, hS_ic)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_toggle_advection

   ! -----------------------------------------------------------------
   ! Direction-split (Lie) form
   ! -----------------------------------------------------------------

   subroutine test_split_cwc_constancy(error)
      !! CWC under `continuity_tracer_step_split`: same
      !! defining test as the unsplit form, lifted to the
      !! interleaved Lie split.  Uniform Tr must stay uniform.  If
      !! the interleaving is wrong (e.g. tracer reads post-zonal h
      !! against mass_flux_x computed from pre-zonal h), Tr drifts
      !! and the test fails.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: S_CONST = 34.5_wp
      real(wp), parameter :: T_CONST = 18.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 30
      integer :: i, j, k, nx, ny, step
      real(wp) :: max_dS, max_dT, h_min
      real(wp), allocatable :: S_new(:, :, :), T_new(:, :, :)
      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H_BASE
         ms%tracers(ms%idx_salinity)%hTr = S_CONST*H_BASE
         ms%tracers(ms%idx_temperature)%hTr = T_CONST*H_BASE
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
         end do

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         h_min = minval(ms%h_layer)
         call check(error, h_min > 0.0_wp, &
                    "h_layer went non-positive — split kernel unstable")
         if (allocated(error)) exit checks

         allocate (S_new(nx, ny, NZ), T_new(nx, ny, NZ))
         S_new = ms%tracers(ms%idx_salinity)%hTr/ms%h_layer
         T_new = ms%tracers(ms%idx_temperature)%hTr/ms%h_layer
         max_dS = maxval(abs(S_new - S_CONST))
         max_dT = maxval(abs(T_new - T_CONST))
         deallocate (S_new, T_new)

         call check(error, max_dS < 1.0e-10_wp, &
                    "split CWC broken: uniform S deviated")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-10_wp, &
                    "split CWC broken: uniform T deviated")

      end block checks
      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_split_cwc_constancy

   subroutine test_split_mass_conservation(error)
      !! Tracer mass per layer conserved to round-off under the
      !! interleaved split step.  Same IC pattern as the unsplit
      !! `test_mass_conservation`; the split form must hit the
      !! 1e-10 floor too — each substep is a face-flux update on
      !! a closed-wall domain, so total mass is conserved by
      !! construction.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_STEPS = 50
      real(wp) :: total_S_init(NZ), total_T_init(NZ)
      real(wp) :: total_S_final(NZ), total_T_final(NZ), drift, layer_phase
      integer :: i, j, k, nx, ny, step
      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            layer_phase = real(k - 1, wp)*PI/real(NZ, wp)
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = H_BASE*( &
                                                             35.0_wp + 0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp) + layer_phase))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = H_BASE*( &
                                                             15.0_wp + 2.0_wp*cos(2.0_wp*PI*real(j, wp)/real(ny, wp) + layer_phase))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp) + layer_phase)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                               sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp) + layer_phase)
               end do
            end do
         end do

         do k = 1, NZ
            total_S_init(k) = sum(ms%tracers(ms%idx_salinity)%hTr(:, :, k))*grid%dx*grid%dy
            total_T_init(k) = sum(ms%tracers(ms%idx_temperature)%hTr(:, :, k))*grid%dx*grid%dy
         end do

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         do k = 1, NZ
            total_S_final(k) = sum(ms%tracers(ms%idx_salinity)%hTr(:, :, k))*grid%dx*grid%dy
            total_T_final(k) = sum(ms%tracers(ms%idx_temperature)%hTr(:, :, k))*grid%dx*grid%dy
            drift = abs(total_S_final(k) - total_S_init(k))/abs(total_S_init(k))
            call check(error, drift < 1.0e-10_wp, &
                       "split-form salinity mass drift exceeded 1e-10")
            if (allocated(error)) exit checks
            drift = abs(total_T_final(k) - total_T_init(k))/abs(total_T_init(k))
            call check(error, drift < 1.0e-10_wp, &
                       "split-form temperature mass drift exceeded 1e-10")
            if (allocated(error)) exit checks
         end do

      end block checks
      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_split_mass_conservation

end module test_ocean_tracer_adv
