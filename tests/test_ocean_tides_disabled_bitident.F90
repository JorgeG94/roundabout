!! Bit-identity guard for the C1 equilibrium-tide barotropic seam.
!!
!! The tide enters `barotropic_substep_nonlinear` through a trailing
!! `optional` `eta_forcing` argument gated inside the existing PGF
!! `do concurrent`.  This test pins three invariants:
!!   1. absent `eta_forcing`  ==  present `eta_forcing = 0`  (byte-for-byte);
!!   2. a spatially CONSTANT `eta_forcing` (zero gradient) is also
!!      byte-for-byte identical to the absent path (a uniform tide exerts
!!      no body force);
!!   3. a spatially VARYING `eta_forcing` DOES change the solution (the
!!      seam is actually live, not dead code).
module test_ocean_tides_disabled_bitident
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_tides_disabled_bitident_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 12, NY_PHYS = 10
   real(wp), parameter :: H_REF = 100.0_wp
   integer, parameter :: N_STEPS = 8
   real(wp), parameter :: DT_INNER = 0.1_wp

contains

   subroutine collect_ocean_tides_disabled_bitident_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("absent_equals_zero_forcing", test_absent_eq_zero), &
                  new_unittest("uniform_tide_is_noop", test_uniform_noop), &
                  new_unittest("varying_tide_is_live", test_varying_live)]
   end subroutine collect_ocean_tides_disabled_bitident_tests

   subroutine seed_ic(dyn, grid, fu, fv)
      !! A non-trivial rest-perturbed IC so the substep produces a
      !! non-zero, forcing-sensitive solution.
      type(ocean_dyn_t), intent(inout) :: dyn
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: fu(:, :), fv(:, :)
      integer :: i, j
      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            dyn%bt_work%bt_eta(i, j) = 0.3_wp*cos(real(i, wp)*0.5_wp) &
                                       *sin(real(j, wp)*0.4_wp)
         end do
      end do
      fu = 1.0e-4_wp
      fv = -5.0e-5_wp
   end subroutine seed_ic

   subroutine run_variant(grid, metrics, dyn, cor, fu, fv, mode, eta)
      !! mode 0 = no eta_forcing; 1 = pass `eta`.  Manages device data.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      real(wp), intent(in) :: fu(:, :), fv(:, :)
      integer, intent(in) :: mode
      real(wp), intent(in) :: eta(:, :)
      !$acc enter data copyin(dyn, cor, fu, fv, eta)
      call dyn%enter_data()
      call cor%enter_data()
      if (mode == 1) then
         call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                           fu, fv, &
                                           N_STEPS, DT_INNER, &
                                           bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                           bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                           eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                           bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                           bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                           uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                           bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                           bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                           bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                           vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                           bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                           bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           eta_forcing=eta, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      else
         call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                           fu, fv, &
                                           N_STEPS, DT_INNER, &
                                           bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                           bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                           eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                           bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                           bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                           uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                           bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                           bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                           bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                           vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                           bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                           bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      end if
      !$acc update self(dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt)
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv, eta)
   end subroutine run_variant

   subroutine test_absent_eq_zero(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :), eta(:, :)
      real(wp), allocatable :: e0(:, :), u0(:, :), v0(:, :)
      real(wp), allocatable :: e1(:, :), u1(:, :), v1(:, :)

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      call cor%init(grid)
      allocate (fu(grid%nx_total + 1, grid%ny_total))
      allocate (fv(grid%nx_total, grid%ny_total + 1))
      allocate (eta(grid%nx_total, grid%ny_total), source=0.0_wp)

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 0, eta)
      e0 = dyn%bt_work%bt_eta
      u0 = dyn%bt_work%bt_ubt
      v0 = dyn%bt_work%bt_vbt

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 1, eta)   ! eta = 0
      e1 = dyn%bt_work%bt_eta
      u1 = dyn%bt_work%bt_ubt
      v1 = dyn%bt_work%bt_vbt
      call destroy_cartesian_metrics(metrics)
      call dyn%destroy()
      call cor%destroy()

      call check(error, maxval(abs(e1 - e0)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(u1 - u0)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(v1 - v0)), 0.0_wp, thr=0.0_wp)
   end subroutine test_absent_eq_zero

   subroutine test_uniform_noop(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :), eta(:, :)
      real(wp), allocatable :: u0(:, :), v0(:, :)
      real(wp), allocatable :: u1(:, :), v1(:, :)

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      call cor%init(grid)
      allocate (fu(grid%nx_total + 1, grid%ny_total))
      allocate (fv(grid%nx_total, grid%ny_total + 1))
      allocate (eta(grid%nx_total, grid%ny_total), source=1.234_wp)   ! uniform

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 0, eta)
      u0 = dyn%bt_work%bt_ubt
      v0 = dyn%bt_work%bt_vbt

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 1, eta)
      u1 = dyn%bt_work%bt_ubt
      v1 = dyn%bt_work%bt_vbt
      call destroy_cartesian_metrics(metrics)
      call dyn%destroy()
      call cor%destroy()

      call check(error, maxval(abs(u1 - u0)), 0.0_wp, thr=0.0_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(v1 - v0)), 0.0_wp, thr=0.0_wp)
   end subroutine test_uniform_noop

   subroutine test_varying_live(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :), eta(:, :)
      real(wp), allocatable :: u0(:, :), u1(:, :)
      integer :: i, j

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      call cor%init(grid)
      allocate (fu(grid%nx_total + 1, grid%ny_total))
      allocate (fv(grid%nx_total, grid%ny_total + 1))
      allocate (eta(grid%nx_total, grid%ny_total))
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            eta(i, j) = 0.2_wp*real(i, wp)   ! non-zero x-gradient
         end do
      end do

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 0, eta)
      u0 = dyn%bt_work%bt_ubt

      call seed_ic(dyn, grid, fu, fv)
      call run_variant(grid, metrics, dyn, cor, fu, fv, 1, eta)
      u1 = dyn%bt_work%bt_ubt
      call destroy_cartesian_metrics(metrics)
      call dyn%destroy()
      call cor%destroy()

      call check(error, maxval(abs(u1 - u0)) > 1.0e-8_wp, &
                 "varying eta_forcing must change u_bt (seam is live)")
   end subroutine test_varying_live

end module test_ocean_tides_disabled_bitident
