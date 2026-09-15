!! Analytical tests for the 2-D barotropic gravity-wave CFL used by
!! `auto_n_inner` (`bt_auto_n_inner` + `metrics_bt_cfl_length`).
!!
!! The barotropic substep must satisfy the 2-D external-gravity-wave CFL
!!   c_ext * dt_bt * sqrt(1/dx^2 + 1/dy^2) <= 1,
!! i.e. the CFL length is `l = 1/sqrt(1/dx^2 + 1/dy^2)`, which on square cells
!! is `dx/sqrt(2)` — NOT the single grid length `dx` the legacy estimate used.
!! These tests pin the formula (the floor, the sqrt(2) factor vs the old 1-D
!! length) and the metrics length helper on uniform + anisotropic grids.
module test_ocean_bt_cfl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_setup, only: bt_auto_n_inner, metrics_bt_cfl_length
   use ocean_test_metrics, only: make_cartesian_metrics, make_anisotropic_metrics, &
                                 destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_bt_cfl_tests

   integer, parameter :: NGHOST = 3
   real(wp), parameter :: SQRT2 = sqrt(2.0_wp)

contains

   subroutine collect_ocean_bt_cfl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bt_n_inner_formula", test_n_inner_formula), &
                  new_unittest("bt_n_inner_floor", test_n_inner_floor), &
                  new_unittest("bt_n_inner_2d_is_sqrt2_of_1d", test_2d_vs_1d), &
                  new_unittest("bt_cfl_length_uniform", test_length_uniform), &
                  new_unittest("bt_cfl_length_anisotropic", test_length_anisotropic) &
                  ]
   end subroutine collect_ocean_bt_cfl_tests

   subroutine test_n_inner_formula(error)
      !! n_inner = max(1, ceiling(dt_outer / (safety*l_cfl/c_ext))).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 300.0_wp, SAFETY = 0.65_wp, C = 198.0_wp
      real(wp) :: l_cfl, dt_bt
      integer :: n, n_expect

      l_cfl = 5000.0_wp/SQRT2                  ! square 5 km cell
      dt_bt = SAFETY*l_cfl/C
      n_expect = max(1, ceiling(DT/dt_bt))
      n = bt_auto_n_inner(DT, SAFETY, C, l_cfl)
      call check(error, n == n_expect, "bt_auto_n_inner: wrong n_inner")
      if (allocated(error)) return
      ! Substep must actually satisfy the CFL: dt_outer/n <= safety*l/c.
      call check(error, DT/real(n, wp) <= SAFETY*l_cfl/C + 1.0e-9_wp, &
                 "bt_auto_n_inner: substep violates CFL")
   end subroutine test_n_inner_formula

   subroutine test_n_inner_floor(error)
      !! Very long CFL length / tiny c => n_inner floored at 1, never 0.
      type(error_type), allocatable, intent(out) :: error
      integer :: n
      n = bt_auto_n_inner(300.0_wp, 0.65_wp, 1.0e-3_wp, 1.0e6_wp)
      call check(error, n == 1, "bt_auto_n_inner: floor not 1")
   end subroutine test_n_inner_floor

   subroutine test_2d_vs_1d(error)
      !! The 2-D length (dx/sqrt2) yields ~sqrt(2)x more substeps than the
      !! legacy 1-D length (dx) — and always at least as many (the bug was
      !! under-counting).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 1200.0_wp, SAFETY = 0.65_wp, C = 198.0_wp, DX = 5000.0_wp
      integer :: n_1d, n_2d
      real(wp) :: ratio

      n_1d = bt_auto_n_inner(DT, SAFETY, C, DX)          ! legacy 1-D length
      n_2d = bt_auto_n_inner(DT, SAFETY, C, DX/SQRT2)    ! correct 2-D length
      call check(error, n_2d >= n_1d, "2-D n_inner must be >= 1-D")
      if (allocated(error)) return
      ratio = real(n_2d, wp)/real(n_1d, wp)
      call check(error, abs(ratio - SQRT2) < 0.05_wp, "2-D/1-D n_inner ratio not ~sqrt(2)")
   end subroutine test_2d_vs_1d

   subroutine test_length_uniform(error)
      !! Uniform Cartesian dx=dy => l_cfl = dx/sqrt(2).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp) :: l_cfl

      call grid%init(16, 12, NGHOST, 5000.0_wp, 5000.0_wp)
      call make_cartesian_metrics(metrics, grid)
      l_cfl = metrics_bt_cfl_length(metrics, grid)
      call check(error, abs(l_cfl - 5000.0_wp/SQRT2) < 1.0e-6_wp, &
                 "uniform l_cfl /= dx/sqrt(2)")
      call destroy_cartesian_metrics(metrics)
   end subroutine test_length_uniform

   subroutine test_length_anisotropic(error)
      !! dxT varies with j (min dxT = dx0 at j=1), dy uniform =>
      !! l_cfl = 1/sqrt(1/dx0^2 + 1/dy^2).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX0 = 4000.0_wp, DY = 6000.0_wp, AMP = 1.0_wp
      real(wp) :: l_cfl, expect, dx_min_interior
      integer :: j_min, ny_total

      call grid%init(16, 12, NGHOST, DX0, DY)
      call make_anisotropic_metrics(metrics, grid, DX0, DY, AMP)
      l_cfl = metrics_bt_cfl_length(metrics, grid)
      ! Helper fills dxT(i,j) = dx0*(1 + amp*(j-1)/ny_total) with j the TOTAL
      ! index; dxT grows with j, so the smallest interior dxT (=> the limiting
      ! 2-D CFL cell) is at the first interior row j = ng+1.
      ny_total = grid%ny_total
      j_min = NGHOST + 1
      dx_min_interior = DX0*(1.0_wp + AMP*real(j_min - 1, wp)/real(ny_total, wp))
      expect = 1.0_wp/sqrt(1.0_wp/dx_min_interior**2 + 1.0_wp/DY**2)
      call check(error, abs(l_cfl - expect) < 1.0e-6_wp, &
                 "anisotropic l_cfl /= 1/sqrt(1/min(dxT)^2+1/dy^2)")
      call destroy_cartesian_metrics(metrics)
   end subroutine test_length_anisotropic

end module test_ocean_bt_cfl
