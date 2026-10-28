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
   use rdb_barotropic_dtbt, only: bt_dtbt_dt_max, bt_dtbt_n_inner, bt_dtbt_safe_step
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
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
                  new_unittest("bt_cfl_wet_no_wet_cell", test_wet_none), &
                  new_unittest("bt_dtbt_uniform_f0", test_dtbt_uniform_f0), &
                  new_unittest("bt_dtbt_coriolis_4f2", test_dtbt_coriolis), &
                  new_unittest("bt_dtbt_anisotropic", test_dtbt_aniso), &
                  new_unittest("bt_dtbt_n_inner_ceiling", test_dtbt_n_inner), &
                  new_unittest("bt_dtbt_reset_follows_depth", test_dtbt_depth_change), &
                  new_unittest("bt_dtbt_land_faces_closed", test_dtbt_land), &
                  new_unittest("bt_dtbt_device_branch_matches_host", test_dtbt_device), &
                  new_unittest("bt_dtbt_config_refusals", test_dtbt_refusals) &
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

   subroutine dtbt_grid(nx, ny, dxv, dyv, h, wet, dyc, dxc, idxc, idyc, ia, fc, depth, fval)
      !! Uniform `nx x ny` grid, `nz = 2` layers of `depth/2`, all wet.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dxv, dyv, depth, fval
      real(wp), intent(out) :: h(nx, ny, 2), wet(nx, ny), dyc(nx + 1, ny), dxc(nx, ny + 1)
      real(wp), intent(out) :: idxc(nx + 1, ny), idyc(nx, ny + 1), ia(nx, ny)
      real(wp), intent(out) :: fc(nx + 1, ny + 1)
      h = 0.5_wp*depth
      wet = 1.0_wp
      dyc = dyv
      dxc = dxv
      idxc = 1.0_wp/dxv
      idyc = 1.0_wp/dyv
      ia = 1.0_wp/(dxv*dyv)
      fc = fval
   end subroutine dtbt_grid

   subroutine test_dtbt_uniform_f0(error)
      !! Flat bottom, f = 0, dx = dy: 1/dt_max^2 = 1/2 (1 + 2 bebt) * 4 g D / dx^2.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 7
      real(wp), parameter :: DX = 12000.0_wp, D = 3000.0_wp, BEBT = 0.2_wp
      real(wp) :: h(NX, NY, 2), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_max, expect

      call dtbt_grid(NX, NY, DX, DX, h, wet, dyc, dxc, idxc, idyc, ia, fc, D, 0.0_wp)
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_max)
      expect = 1.0_wp/sqrt(0.7_wp*4.0_wp*GRAVITY*D/DX**2)
      call check(error, abs(dt_max - expect) <= 1.0e-14_wp*expect, &
                 "uniform f=0 dt_max /= 1/sqrt(0.7*4*g*D/dx^2)")
   end subroutine test_dtbt_uniform_f0

   subroutine test_dtbt_coriolis(error)
      !! The four corner f^2 terms enter as 4 f^2 inside the same 1/2(1+2 bebt).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 7
      real(wp), parameter :: DX = 12000.0_wp, D = 3000.0_wp, BEBT = 0.2_wp, F0 = 1.0e-3_wp
      real(wp) :: h(NX, NY, 2), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_max, expect

      call dtbt_grid(NX, NY, DX, DX, h, wet, dyc, dxc, idxc, idyc, ia, fc, D, F0)
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_max)
      expect = 1.0_wp/sqrt(0.7_wp*(4.0_wp*GRAVITY*D/DX**2 + 4.0_wp*F0**2))
      call check(error, abs(dt_max - expect) <= 1.0e-14_wp*expect, &
                 "f-plane dt_max /= 1/sqrt(0.7*(4gD/dx^2 + 4f^2))")
      if (allocated(error)) return
      ! Quiescent-pressure limit: f alone (g = 0) gives 1/(f*sqrt(0.7*4)).
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, 0.0_wp, BEBT, .false., dt_max)
      expect = 1.0_wp/(F0*sqrt(0.7_wp*4.0_wp))
      call check(error, abs(dt_max - expect) <= 1.0e-14_wp*expect, "g=0 dt_max /= 1/(2 f sqrt(0.7))")
   end subroutine test_dtbt_coriolis

   subroutine test_dtbt_aniso(error)
      !! dx /= dy: E+W faces give 2 g D/dx^2, N+S give 2 g D/dy^2.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 6, NY = 6
      real(wp), parameter :: DX = 20000.0_wp, DY = 5000.0_wp, D = 800.0_wp, BEBT = 0.2_wp
      real(wp) :: h(NX, NY, 2), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_max, expect

      call dtbt_grid(NX, NY, DX, DY, h, wet, dyc, dxc, idxc, idyc, ia, fc, D, 0.0_wp)
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_max)
      expect = 1.0_wp/sqrt(0.7_wp*(2.0_wp*GRAVITY*D/DX**2 + 2.0_wp*GRAVITY*D/DY**2))
      call check(error, abs(dt_max - expect) <= 1.0e-14_wp*expect, &
                 "anisotropic dt_max /= 1/sqrt(0.7*(2gD/dx^2 + 2gD/dy^2))")
   end subroutine test_dtbt_aniso

   subroutine test_dtbt_n_inner(error)
      !! n_inner = ceiling(dt / step); dtbt = 0 means -0.98; dtbt > 0 is seconds.
      type(error_type), allocatable, intent(out) :: error
      call check(error, bt_dtbt_n_inner(900.0_wp, -0.9_wp, 100.0_wp) == 10, "step 90 s: 900/90 -> 10")
      if (allocated(error)) return
      call check(error, bt_dtbt_n_inner(900.0_wp, -0.9_wp, 99.0_wp) == 11, "step 89.1 s: -> 11")
      if (allocated(error)) return
      call check(error, abs(bt_dtbt_safe_step(0.0_wp, 100.0_wp) - 98.0_wp) < 1.0e-12_wp, "dtbt=0 /= -0.98")
      if (allocated(error)) return
      call check(error, bt_dtbt_n_inner(900.0_wp, 0.0_wp, 100.0_wp) == 10, "dtbt=0: 900/98 -> 10")
      if (allocated(error)) return
      call check(error, bt_dtbt_n_inner(900.0_wp, 300.0_wp, 1.0_wp) == 3, "fixed 300 s: -> 3")
      if (allocated(error)) return
      call check(error, bt_dtbt_n_inner(900.0_wp, 250.0_wp, 1.0_wp) == 4, "fixed 250 s: -> 4")
      if (allocated(error)) return
      call check(error, bt_dtbt_n_inner(900.0_wp, -0.9_wp, huge(1.0_wp)) == 1, "no wet cell -> 1")
   end subroutine test_dtbt_n_inner

   subroutine test_dtbt_depth_change(error)
      !! Period 0 recomputes from the live depths: with f = 0, quadrupling the
      !! depth halves dt_max exactly (dt_max ~ 1/sqrt(D)).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 7
      real(wp), parameter :: DX = 12000.0_wp, BEBT = 0.2_wp
      real(wp) :: h(NX, NY, 2), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_a, dt_b

      call dtbt_grid(NX, NY, DX, DX, h, wet, dyc, dxc, idxc, idyc, ia, fc, 1000.0_wp, 0.0_wp)
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_a)
      h = 4.0_wp*h
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_b)
      call check(error, abs(dt_b - 0.5_wp*dt_a) <= 1.0e-14_wp*dt_a, "4x depth did not halve dt_max")
      if (allocated(error)) return
      ! One deep cell sets the domain minimum, and a deep cell's neighbour
      ! faces carry the average depth.
      call dtbt_grid(NX, NY, DX, DX, h, wet, dyc, dxc, idxc, idyc, ia, fc, 1000.0_wp, 0.0_wp)
      h(4, 4, :) = 2500.0_wp
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_b)
      ! cell (4,4): four faces of depth (5000+1000)/2 = 3000 each
      dt_a = 1.0_wp/sqrt(0.7_wp*4.0_wp*GRAVITY*3000.0_wp/DX**2)
      call check(error, dt_b <= dt_a*(1.0_wp + 1.0e-14_wp), "deep cell does not set the minimum")
   end subroutine test_dtbt_depth_change

   subroutine test_dtbt_land(error)
      !! A land neighbour closes the face; an all-land window is inert (huge).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 6, NY = 6
      real(wp), parameter :: DX = 10000.0_wp, D = 1000.0_wp, BEBT = 0.2_wp
      real(wp) :: h(NX, NY, 2), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_max, expect

      call dtbt_grid(NX, NY, DX, DX, h, wet, dyc, dxc, idxc, idyc, ia, fc, D, 0.0_wp)
      wet = 0.0_wp
      wet(3, 3) = 1.0_wp                   ! a lone wet cell: all four faces closed
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_max)
      call check(error, dt_max == huge(1.0_wp), "isolated wet cell (no open face, f=0) is not unconstrained")
      if (allocated(error)) return
      wet(4, 3) = 1.0_wp                   ! a two-cell pond: one open face each
      call bt_dtbt_dt_max(NX, NY, 2, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_max)
      expect = 1.0_wp/sqrt(0.7_wp*GRAVITY*D/DX**2)
      call check(error, abs(dt_max - expect) <= 1.0e-14_wp*expect, "one open face: dt_max /= 1/sqrt(0.7 g D/dx^2)")
   end subroutine test_dtbt_land

   subroutine test_dtbt_device(error)
      !! The `on_device` branch (the `do concurrent` reduction `engine_step`
      !! runs) must reproduce the host loop on a bathymetry with a deep cell,
      !! a land block and an f gradient.  Every array is mapped by hand
      !! (`mem:separate`: no implicit copies); a no-op on host builds.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 12, NY = 10, NZ = 2
      real(wp), parameter :: DX = 15000.0_wp, DY = 9000.0_wp, BEBT = 0.2_wp
      real(wp) :: h(NX, NY, NZ), wet(NX, NY), dyc(NX + 1, NY), dxc(NX, NY + 1)
      real(wp) :: idxc(NX + 1, NY), idyc(NX, NY + 1), ia(NX, NY), fc(NX + 1, NY + 1)
      real(wp) :: dt_host, dt_dev
      integer :: i, j

      call dtbt_grid(NX, NY, DX, DY, h, wet, dyc, dxc, idxc, idyc, ia, fc, 2000.0_wp, 0.0_wp)
      h(5, 4, :) = 2800.0_wp
      h(6, 4, 1) = 90.0_wp
      wet(8:9, 5:7) = 0.0_wp
      do j = 1, NY + 1
         do i = 1, NX + 1
            fc(i, j) = 1.0e-5_wp*real(j, wp)
         end do
      end do
      call bt_dtbt_dt_max(NX, NY, NZ, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .false., dt_host)
      !$acc enter data copyin(h, wet, dyc, dxc, idxc, idyc, ia, fc)
      call bt_dtbt_dt_max(NX, NY, NZ, 2, NX - 1, 2, NY - 1, h, wet, dyc, dxc, idxc, idyc, ia, &
                          fc, GRAVITY, BEBT, .true., dt_dev)
      !$acc exit data delete(h, wet, dyc, dxc, idxc, idyc, ia, fc)
      call check(error, abs(dt_dev - dt_host) <= 1.0e-13_wp*dt_host, &
                 "device dt_max differs from the host loop")
   end subroutine test_dtbt_device

   subroutine test_dtbt_refusals(error)
      !! dtbt replaces auto_n_inner/cfl_bt_safety; a reset period needs dtbt.
      type(error_type), allocatable, intent(out) :: error
      call expect_dtbt(error, "dtbt = -0.9", .true., "dtbt alone")
      if (allocated(error)) return
      call expect_dtbt(error, "dtbt = -0.9, dtbt_reset_period = 0.0", .true., "dtbt + reset 0")
      if (allocated(error)) return
      call expect_dtbt(error, "dtbt = -0.9, auto_n_inner = .true.", .false., "dtbt + auto_n_inner")
      if (allocated(error)) return
      call expect_dtbt(error, "dtbt = -0.9, cfl_bt_safety = 0.5", .false., "dtbt + cfl_bt_safety")
      if (allocated(error)) return
      call expect_dtbt(error, "dtbt = -1.5", .false., "dtbt < -1")
      if (allocated(error)) return
      call expect_dtbt(error, "dtbt_reset_period = 0.0", .false., "reset period without dtbt")
   end subroutine test_dtbt_refusals

   subroutine expect_dtbt(error, bt_body, valid, what)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: bt_body, what
      logical, intent(in) :: valid
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ierr
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 400.0 /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true. /"//new_line("a")// &
            "&ocean_bt_nml n_inner = 8, "//bt_body//" /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "namelist must PARSE: "//what)
      if (allocated(error)) return
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      if (valid) then
         call check(error, ierr == OCEAN_STATUS_OK, "must CONFIGURE: "//what)
      else
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, "must be REFUSED: "//what)
      end if
   end subroutine expect_dtbt

end module test_ocean_bt_cfl
