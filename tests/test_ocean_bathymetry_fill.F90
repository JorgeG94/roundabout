!! Regression tests for the ghost-cell bathymetry bug (commit 20b1b84).
!!
!! Both `set_bathymetry_spoon` and `set_bathymetry_seamount` originally
!! restricted their loops to interior cells (`ng+1 : ng+nx_phys`),
!! leaving the ghost rows at the alloc-time zero value.  Downstream:
!!
!!   1. `seed_h_layer_uniform_impl` does `h_layer = b/nz` → ghost
!!      cells get `h_layer = 0`.
!!   2. The EOS hits its vanishing-layer fallback (`rho_layer = rho_0`).
!!   3. Interior cells get the actual EOS-computed density.  Whenever
!!      the IC has `T_init ≠ T_ref` or `S_init ≠ S_ref`, that density
!!      differs from `rho_0`, putting a step jump at every wall-
!!      adjacent face.
!!   4. The BPG operator picks up the spurious density gradient and
!!      injects horizontal pressure → exponential 12-hr e-fold
!!      instability over ~10 outer steps (caught by the seamount
!!      Tier-1.5 test as NaN by day 1).
!!
!! These tests assert that both bathymetry routines fill the full
!! (1 : nx_total, 1 : ny_total) extent of `b`.
module test_ocean_bathymetry_fill
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, DEG2RAD
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: set_bathymetry_spoon, set_bathymetry_seamount, &
                              topo_length_to_grid_units
   implicit none
   private

   public :: collect_ocean_bathymetry_fill_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RAD_EARTH = 6.378e6_wp

contains

   subroutine collect_ocean_bathymetry_fill_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("seamount_fills_ghost_rows", test_seamount_ghost_fill), &
                  new_unittest("spoon_fills_ghost_rows", test_spoon_ghost_fill), &
                  new_unittest("length_to_grid_units_converts", test_length_conversion), &
                  new_unittest("spherical_seamount_not_flat", test_spherical_seamount) &
                  ]
   end subroutine collect_ocean_bathymetry_fill_tests

   subroutine test_length_conversion(error)
      !! `topo_length_to_grid_units` passes metres through unchanged on
      !! Cartesian, and converts to degrees of latitude on spherical/
      !! curvilinear grids (metres-per-degree = rad_earth·π/180).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m_per_deg
      m_per_deg = RAD_EARTH*DEG2RAD
      checks: block
         ! Cartesian: identity.
         call check(error, abs(topo_length_to_grid_units(25000.0_wp, "cartesian", RAD_EARTH) &
                               - 25000.0_wp) < 1.0e-9_wp, "cartesian passes metres through")
         if (allocated(error)) exit checks
         ! Spherical: one metre-per-degree of length => 1 degree.
         call check(error, abs(topo_length_to_grid_units(m_per_deg, "spherical", RAD_EARTH) &
                               - 1.0_wp) < 1.0e-9_wp, "spherical metres->degrees")
         if (allocated(error)) exit checks
         call check(error, abs(topo_length_to_grid_units(m_per_deg, "tripolar", RAD_EARTH) &
                               - 1.0_wp) < 1.0e-9_wp, "tripolar metres->degrees")
      end block checks
   end subroutine test_length_conversion

   subroutine test_spherical_seamount(error)
      !! Regression for the spherical-grid seamount unit bug: a metres
      !! `slope_scale` against a degrees grid position made
      !! `exp(-r²/L²) ≈ 1` everywhere, collapsing the basin to a flat
      !! `peak_depth`.  With the metres->degrees conversion the seamount is
      !! a real bump: the basin must span [peak_depth, ~max_depth].  The
      !! test also confirms the UNCONVERTED (buggy) call is ~flat, so it
      !! fails if the conversion is dropped.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: MAX_DEPTH = 4000.0_wp, PEAK_DEPTH = 200.0_wp
      real(wp), parameter :: HALF_WIDTH_M = 5.0e5_wp   ! 500 km
      real(wp) :: hw_grid
      integer :: nx, ny
      ! 40x40 lon-lat sector at 1 degree spacing (positions in DEGREES).
      call grid%init(40, 40, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (b(nx, ny), source=0.0_wp)
      checks: block
         ! Fixed path: convert the metres half-width to grid degrees.
         hw_grid = topo_length_to_grid_units(HALF_WIDTH_M, "spherical", RAD_EARTH)
         call set_bathymetry_seamount(b, grid, MAX_DEPTH, PEAK_DEPTH, hw_grid)
         call check(error, maxval(b) - minval(b) > 1000.0_wp, &
                    "spherical seamount must vary by > 1000 m (not a flat basin)")
         if (allocated(error)) exit checks
         call check(error, minval(b) >= PEAK_DEPTH - 1.0e-6_wp .and. &
                    maxval(b) <= MAX_DEPTH + 1.0e-6_wp, "depth within [peak, max]")
         if (allocated(error)) exit checks
         ! Buggy path (raw metres half-width on a degrees grid) is ~flat:
         ! guards that the test would catch a regression that drops the fix.
         call set_bathymetry_seamount(b, grid, MAX_DEPTH, PEAK_DEPTH, HALF_WIDTH_M)
         call check(error, maxval(b) - minval(b) < 1.0_wp, &
                    "unconverted metres half-width collapses the basin flat (the bug)")
      end block checks
   end subroutine test_spherical_seamount

   subroutine test_seamount_ghost_fill(error)
      !! After `set_bathymetry_seamount`, every cell of `b` (interior
      !! + ghost rows on all four sides) must be positive and within
      !! [peak_depth, max_depth].  Pre-fix the ghost rows held the
      !! zero from the alloc-time `b(:,:) = 0` left by ms%init —
      !! that's the smoking-gun condition the EOS fallback bites on.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: MAX_DEPTH = 4000.0_wp
      real(wp), parameter :: PEAK_DEPTH = 200.0_wp
      real(wp), parameter :: HALF_WIDTH = 25.0_wp
      integer :: nx, ny

      call grid%init(20, 16, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (b(nx, ny), source=0.0_wp)

      call set_bathymetry_seamount(b, grid, MAX_DEPTH, PEAK_DEPTH, HALF_WIDTH)

      ! All cells, including the four ghost bands, must hold a real
      ! bathymetry value (> 0) within the formula's range.  Pre-fix
      ! the first/last NGHOST rows and columns stayed at zero.
      call check(error, minval(b) > 0.0_wp, &
                 "set_bathymetry_seamount: some cell left at zero (ghost not filled)")
      if (allocated(error)) return
      call check(error, minval(b) >= PEAK_DEPTH - 1.0e-6_wp, &
                 "seamount: b dips below peak_depth")
      if (allocated(error)) return
      call check(error, maxval(b) <= MAX_DEPTH + 1.0e-6_wp, &
                 "seamount: b exceeds max_depth")
      deallocate (b)
   end subroutine test_seamount_ghost_fill

   subroutine test_spoon_ghost_fill(error)
      !! Mirror of the seamount test for the MOM6-style spoon
      !! bathymetry — same ghost-row coverage bug, same fix.  Spoon
      !! parameters: edge_depth at the basin perimeter, max_depth in
      !! the south interior, exponential decay toward the north wall.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: MAX_DEPTH = 4000.0_wp
      real(wp), parameter :: EDGE_DEPTH = 100.0_wp
      real(wp), parameter :: SLOPE_SCALE = 5.0_wp
      integer :: nx, ny

      call grid%init(16, 12, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (b(nx, ny), source=0.0_wp)

      call set_bathymetry_spoon(b, grid, MAX_DEPTH, EDGE_DEPTH, SLOPE_SCALE)

      call check(error, minval(b) > 0.0_wp, &
                 "set_bathymetry_spoon: some cell left at zero (ghost not filled)")
      if (allocated(error)) return
      ! Spoon clamps to [edge_depth, max_depth] inside the formula.
      call check(error, minval(b) >= EDGE_DEPTH - 1.0e-6_wp, &
                 "spoon: b dips below edge_depth")
      if (allocated(error)) return
      call check(error, maxval(b) <= MAX_DEPTH + 1.0e-6_wp, &
                 "spoon: b exceeds max_depth")
      deallocate (b)
   end subroutine test_spoon_ghost_fill

end module test_ocean_bathymetry_fill
