!! Analytical tests for the 2-D barotropic gravity-wave CFL used by
!! `auto_n_inner` (`bt_auto_n_inner` + `metrics_bt_cfl_length`).
!!
!! The barotropic substep must satisfy the 2-D external-gravity-wave CFL
!!   c_ext * dt_bt * sqrt(1/dx^2 + 1/dy^2) <= 1,
!! i.e. the CFL length is `l = 1/sqrt(1/dx^2 + 1/dy^2)`, which on square cells
!! is `dx/sqrt(2)` — NOT the single grid length `dx` the legacy estimate used.
!! These tests pin the formula (the floor, the sqrt(2) factor vs the old 1-D
!! length) and the metrics length helper on uniform + anisotropic grids.
!!
!! The per-wet-cell limit (`bt_cfl_dt_wet`, MOM6 `set_dtbt`) is pinned
!! against the global-extremes estimate it replaced: bit-identical where
!! the deepest and the smallest wet cell coincide, land excluded, and a
!! deep-coarse / shallow-fine pair limited by the WORSE point, not by the
!! combination of the two extremes.
module test_ocean_bt_cfl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_setup, only: bt_auto_n_inner, metrics_bt_cfl_length, &
                              bt_auto_n_inner_from_dt, bt_cfl_dt_wet
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
                  new_unittest("bt_cfl_length_anisotropic", test_length_anisotropic), &
                  new_unittest("bt_cfl_wet_coincident_is_bitwise_legacy", &
                               test_wet_coincident_bitwise), &
                  new_unittest("bt_cfl_wet_ignores_land", test_wet_ignores_land), &
                  new_unittest("bt_cfl_wet_is_per_point", test_wet_per_point), &
                  new_unittest("bt_cfl_wet_no_wet_cell", test_wet_none) &
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

   ! -----------------------------------------------------------------
   ! Per-wet-cell limit (`bt_cfl_dt_wet`)
   ! -----------------------------------------------------------------

   pure subroutine fill_uniform(nx, ny, b, wet, dxT, dyT, h, dx, dy)
      !! Flat bed `h`, all wet, uniform `dx` x `dy`.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: b(nx, ny), wet(nx, ny), dxT(nx, ny), dyT(nx, ny)
      real(wp), intent(in) :: h, dx, dy
      b = h
      wet = 1.0_wp
      dxT = dx
      dyT = dy
   end subroutine fill_uniform

   subroutine test_wet_coincident_bitwise(error)
      !! Flat bed, uniform grid, all wet: the deepest and the smallest cell
      !! are the SAME cell, so the per-point minimum must equal the legacy
      !! `safety*l_cfl/c_ext` from the global extremes BIT FOR BIT — and so
      !! must `n_inner`.  This is the byte-identity promise for every
      !! all-wet flat configuration (the double gyre among them).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 10, NY = 8
      real(wp), parameter :: DT = 1200.0_wp, SAFETY = 0.65_wp
      real(wp) :: b(NX, NY), wet(NX, NY), dxT(NX, NY), dyT(NX, NY)
      real(wp) :: dt_bt, h_at, l_at, l_old, c_old, dt_old
      integer :: n_wet

      call fill_uniform(NX, NY, b, wet, dxT, dyT, 4000.0_wp, 20000.0_wp, 17000.0_wp)
      call bt_cfl_dt_wet(NX, NY, 2, NX - 1, 2, NY - 1, b, wet, dxT, dyT, SAFETY, &
                         dt_bt, h_at, l_at, n_wet)
      ! The legacy arithmetic, spelled exactly as configure_ocean_bt_split
      ! spelled it before the per-point rewrite.
      l_old = 1.0_wp/sqrt(maxval(1.0_wp/dxT(2:NX - 1, 2:NY - 1)**2 &
                                 + 1.0_wp/dyT(2:NX - 1, 2:NY - 1)**2))
      c_old = sqrt(GRAVITY*max(maxval(b(2:NX - 1, 2:NY - 1)), 1.0_wp))
      dt_old = SAFETY*l_old/c_old
      call check(error, n_wet == (NX - 2)*(NY - 2), "wet cells not all scanned")
      if (allocated(error)) return
      call check(error, dt_bt == dt_old, "coincident extremes: dt_bt not bit-identical to legacy")
      if (allocated(error)) return
      call check(error, bt_auto_n_inner_from_dt(DT, dt_bt) == &
                 bt_auto_n_inner(DT, SAFETY, c_old, l_old), &
                 "coincident extremes: n_inner differs from legacy")
   end subroutine test_wet_coincident_bitwise

   subroutine test_wet_ignores_land(error)
      !! The 1-degree tripolar defect in miniature: 6000 m of ocean on
      !! 100 km cells, plus ONE 362 m land cell (a land-locked bipole) that
      !! also carries a deep bed value.  Land-inclusive extremes buy
      !! thousands of substeps; the wet-only limit is set by the ocean alone.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 10, NY = 8
      real(wp), parameter :: DT = 1800.0_wp, SAFETY = 0.65_wp
      real(wp) :: b(NX, NY), wet(NX, NY), dxT(NX, NY), dyT(NX, NY)
      real(wp) :: dt_bt, h_at, l_at, dt_ocean, dt_land_incl
      integer :: n_wet, n_new, n_old

      call fill_uniform(NX, NY, b, wet, dxT, dyT, 6000.0_wp, 1.0e5_wp, 1.0e5_wp)
      wet(5, 4) = 0.0_wp
      dxT(5, 4) = 362.0_wp
      dyT(5, 4) = 362.0_wp
      b(5, 4) = 7000.0_wp
      call bt_cfl_dt_wet(NX, NY, 2, NX - 1, 2, NY - 1, b, wet, dxT, dyT, SAFETY, &
                         dt_bt, h_at, l_at, n_wet)
      dt_ocean = SAFETY*(1.0e5_wp/sqrt(2.0_wp))/sqrt(GRAVITY*6000.0_wp)
      dt_land_incl = SAFETY*(362.0_wp/sqrt(2.0_wp))/sqrt(GRAVITY*7000.0_wp)
      n_new = bt_auto_n_inner_from_dt(DT, dt_bt)
      n_old = bt_auto_n_inner_from_dt(DT, dt_land_incl)
      call check(error, n_wet == (NX - 2)*(NY - 2) - 1, "land cell was scanned")
      if (allocated(error)) return
      call check(error, abs(dt_bt - dt_ocean) <= 1.0e-12_wp*dt_ocean, &
                 "wet-only dt_bt is not the ocean cells' limit")
      if (allocated(error)) return
      call check(error, h_at == 6000.0_wp, "limiting cell is not an ocean cell")
      if (allocated(error)) return
      call check(error, n_new < n_old/100, "land cell still drives n_inner")
   end subroutine test_wet_ignores_land

   subroutine test_wet_per_point(error)
      !! A deep COARSE cell and a shallow FINE cell: the limit is the worse
      !! of the two per-point limits, and strictly looser than pairing the
      !! deep cell's depth with the fine cell's length (which no point has).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 6, NY = 6
      real(wp), parameter :: SAFETY = 0.65_wp
      real(wp) :: b(NX, NY), wet(NX, NY), dxT(NX, NY), dyT(NX, NY)
      real(wp) :: dt_bt, h_at, l_at, dt_deep, dt_shelf, dt_combined
      integer :: n_wet

      ! Background: shallow + coarse everywhere, so it never limits.
      call fill_uniform(NX, NY, b, wet, dxT, dyT, 100.0_wp, 5.0e4_wp, 5.0e4_wp)
      b(3, 3) = 5000.0_wp                  ! deep, coarse
      b(4, 4) = 50.0_wp                    ! shallow, fine
      dxT(4, 4) = 2.0e3_wp
      dyT(4, 4) = 2.0e3_wp
      call bt_cfl_dt_wet(NX, NY, 2, NX - 1, 2, NY - 1, b, wet, dxT, dyT, SAFETY, &
                         dt_bt, h_at, l_at, n_wet)
      dt_deep = SAFETY*(5.0e4_wp/sqrt(2.0_wp))/sqrt(GRAVITY*5000.0_wp)
      dt_shelf = SAFETY*(2.0e3_wp/sqrt(2.0_wp))/sqrt(GRAVITY*50.0_wp)
      dt_combined = SAFETY*(2.0e3_wp/sqrt(2.0_wp))/sqrt(GRAVITY*5000.0_wp)
      call check(error, abs(dt_bt - min(dt_deep, dt_shelf)) <= 1.0e-12_wp*dt_bt, &
                 "per-point limit is not min over points")
      if (allocated(error)) return
      call check(error, dt_bt > 5.0_wp*dt_combined, &
                 "per-point limit collapsed to the global-extremes combination")
   end subroutine test_wet_per_point

   subroutine test_wet_none(error)
      !! All land: no constraint (`huge`, the min-reduction identity), and
      !! `n_inner` floors at 1.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 5
      real(wp) :: b(NX, NY), wet(NX, NY), dxT(NX, NY), dyT(NX, NY)
      real(wp) :: dt_bt, h_at, l_at
      integer :: n_wet

      call fill_uniform(NX, NY, b, wet, dxT, dyT, 100.0_wp, 1.0e3_wp, 1.0e3_wp)
      wet = 0.0_wp
      call bt_cfl_dt_wet(NX, NY, 2, NX - 1, 2, NY - 1, b, wet, dxT, dyT, 0.65_wp, &
                         dt_bt, h_at, l_at, n_wet)
      call check(error, n_wet == 0 .and. dt_bt == huge(1.0_wp), "all-land window not inert")
      if (allocated(error)) return
      call check(error, bt_auto_n_inner_from_dt(1200.0_wp, dt_bt) == 1, "all-land n_inner /= 1")
   end subroutine test_wet_none

end module test_ocean_bt_cfl
