!! Fast-loop Coriolis/advection reference subtraction for the
!! split-explicit barotropic forcing (`subtract_fast_cor_ref`).
!!
!! `F_bt_u/v` (the depth-mean slow forcing handed to the barotropic
!! substep) contains the depth mean of the layer Coriolis-advection
!! tendencies, and the substep integrates its own live `(ζ+f)·v − ∇KE`
!! on top — so without a reference subtraction the barotropic Coriolis
!! is integrated twice (MOM6 removes it with `Cor_ref_u/v`,
!! `MOM_barotropic.F90:1526-1535`).  `subtract_fast_cor_ref` subtracts
!! the fast-loop terms evaluated at the reference barotropic velocity
!! `bt_work%cor_ref_u/v` from `F_bt_u_fast`/`F_bt_v_fast`.  (That
!! reference is a named slot filled by `set_cor_ref_velocity`: the
!! stage-entry `bt_ubt/bt_vbt` under `ssp_rk2`, the depth mean of
!! `u_av/v_av` under `pred_corr`.  These cases drive it directly.)
!!
!! Covers:
!!   1. Uniform zonal flow `u0` on an f-plane: the fast-loop reference
!!      is exactly `−f·u0` at every interior v-face (ζ = 0, ∇KE = 0),
!!      so the subtraction ADDS `+f·u0` to `F_bt_v_fast` and leaves
!!      `F_bt_u_fast` untouched (v = 0).  Analytic, uniform grid.
!!   2. Rest state: `u = v = 0` ⇒ both forcing fields byte-unchanged
!!      (the rest-preservation gate — the subtraction cannot perturb a
!!      resting ocean).
module test_ocean_bt_cor_ref
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_barotropic_coupling, only: subtract_fast_cor_ref
   use rdb_ocean_boundary_types, only: OBC_WALL
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_bt_cor_ref_tests

   integer, parameter :: NX = 12, NY = 10, NGHOST = 3
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: F0 = 1.0e-4_wp   !! f-plane Coriolis (1/s)
   real(wp), parameter :: U0 = 0.3_wp      !! uniform zonal flow (m/s)

contains

   subroutine collect_ocean_bt_cor_ref_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cor_ref_uniform_flow_fplane", test_uniform_flow), &
                  new_unittest("cor_ref_rest_state_no_op", test_rest_no_op) &
                  ]
   end subroutine collect_ocean_bt_cor_ref_tests

   subroutine run_subtract(grid, metrics, dyn, f_corner)
      !! Map, run the kernel, pull the forcing fields back.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      real(wp), intent(in) :: f_corner(:, :)  ! assumed-shape-ok: host-side test driver
      !$acc enter data copyin(dyn, f_corner)
      call dyn%enter_data()
      call subtract_fast_cor_ref(grid, metrics, dyn%bt_work, f_corner, &
                                 OBC_WALL, OBC_WALL, OBC_WALL, OBC_WALL, &
                                 .true., .true., .true., .true.)
      !$acc update self(dyn%bt_work%F_bt_u_fast, dyn%bt_work%F_bt_v_fast)
      call dyn%exit_data()
      !$acc exit data delete(dyn, f_corner)
   end subroutine run_subtract

   subroutine test_uniform_flow(error)
      !! Uniform `u = U0`, `v = 0`, f-plane: reference `−f·u0` at every
      !! interior v-face ⇒ `F_bt_v_fast` gains exactly `+f·u0`;
      !! `F_bt_u_fast` stays zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: f_corner(:, :)
      real(wp) :: worst
      integer :: i, j

      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid, nz_ml=2)
         dyn%bt_work%cor_ref_u = U0
         dyn%bt_work%cor_ref_v = 0.0_wp
         dyn%bt_work%F_bt_u_fast = 0.0_wp
         dyn%bt_work%F_bt_v_fast = 0.0_wp
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1), source=F0)

         call run_subtract(grid, metrics, dyn, f_corner)

         ! Interior v-faces (j = 2..ny_total): expect exactly +f·u0.
         worst = 0.0_wp
         do j = 2, grid%ny_total
            do i = 1, grid%nx_total
               worst = max(worst, abs(dyn%bt_work%F_bt_v_fast(i, j) - F0*U0))
            end do
         end do
         call check(error, worst < 1.0e-14_wp, &
                    "uniform flow: F_bt_v_fast /= +f*u0 at an interior v-face")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%F_bt_u_fast)) < 1.0e-14_wp, &
                    "uniform flow: F_bt_u_fast perturbed (v = 0 so the u-reference is 0)")
      end block checks
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_uniform_flow

   subroutine test_rest_no_op(error)
      !! Rest state: the subtraction must leave both forcing fields
      !! byte-unchanged (guards the flat-bottom rest-preservation gate).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: f_corner(:, :)
      real(wp), parameter :: F_SEED = 0.123_wp

      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid, nz_ml=2)
         dyn%bt_work%cor_ref_u = 0.0_wp
         dyn%bt_work%cor_ref_v = 0.0_wp
         dyn%bt_work%F_bt_u_fast = F_SEED
         dyn%bt_work%F_bt_v_fast = F_SEED
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1), source=F0)

         call run_subtract(grid, metrics, dyn, f_corner)

         call check(error, maxval(abs(dyn%bt_work%F_bt_u_fast - F_SEED)) == 0.0_wp, &
                    "rest state: F_bt_u_fast changed")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%F_bt_v_fast - F_SEED)) == 0.0_wp, &
                    "rest state: F_bt_v_fast changed")
      end block checks
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_rest_no_op

end module test_ocean_bt_cor_ref
