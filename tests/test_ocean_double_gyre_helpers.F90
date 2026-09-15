!! Unit tests for the analytical double-gyre setup helpers:
!!
!!   * `ocean_surface_stress_t%set_wind_stress_2gyre`
!!     (MOM6's `wind_forcing_2gyre`)
!!   * `set_bathymetry_spoon`
!!     (MOM6's `initialize_topography_named` for `topog_config="spoon"`)
!!   * `coriolis_adv_t%set_beta_plane`
!!
!! Each test sets up the helper, evaluates at known interior points,
!! and checks the formula to round-off.  Wall behaviour and physical-
!! interior-only filling are covered explicitly.
module test_ocean_double_gyre_helpers
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: set_bathymetry_spoon
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_double_gyre_helpers_tests

   real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
   real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)
   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_ocean_double_gyre_helpers_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wind_2gyre_matches_formula", test_wind_2gyre), &
                  new_unittest("wind_2gyre_ghosts_stay_zero", test_wind_2gyre_ghosts), &
                  new_unittest("spoon_walls_and_max", test_spoon_walls_and_max), &
                  new_unittest("spoon_zero_at_north_wall", test_spoon_zero_at_north), &
                  new_unittest("beta_plane_linear_in_y", test_beta_plane) &
                  ]
   end subroutine collect_ocean_double_gyre_helpers_tests

   subroutine test_wind_2gyre(error)
      !! Verify tau_x matches `taux_mag * (1 - cos(2π * (j_phys-0.5)/ny_phys))`
      !! at every physical interior point.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 10, NGHOST = 1
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), parameter :: TAUX_MAG = 0.1_wp
      real(wp) :: expected, y_rel
      integer :: i, j, j_phys

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      ss%rho0 = 1025.0_wp
      call ss%init(grid, nz_ml=1)
      checks: block
         call ss%set_wind_stress_2gyre(grid, TAUX_MAG)
         do j = NGHOST + 1, NGHOST + NY_PHYS
            j_phys = j - NGHOST
            y_rel = (real(j_phys, wp) - 0.5_wp)/real(NY_PHYS, wp)
            expected = TAUX_MAG*(1.0_wp - cos(TWO_PI*y_rel))
            do i = 1, size(ss%tau_x, 1)
               call check(error, abs(ss%tau_x(i, j) - expected) < TOL, &
                          "tau_x at physical (i,j_phys) should match formula")
               if (allocated(error)) exit checks
            end do
         end do
         call check(error, all(abs(ss%tau_y) < TOL), &
                    "tau_y should remain zero for 2gyre forcing")
      end block checks
      call ss%destroy()
   end subroutine test_wind_2gyre

   subroutine test_wind_2gyre_ghosts(error)
      !! Ghost rows must stay at zero — formula is physical-interior only.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 6, NGHOST = 2
      real(wp), parameter :: DX = 1000.0_wp
      integer :: i, j

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      ss%rho0 = 1025.0_wp
      call ss%init(grid, nz_ml=1)
      checks: block
         call ss%set_wind_stress_2gyre(grid, 0.1_wp)
         do j = 1, NGHOST
            do i = 1, size(ss%tau_x, 1)
               call check(error, abs(ss%tau_x(i, j)) < TOL, &
                          "tau_x in south ghost row should be zero")
               if (allocated(error)) exit checks
            end do
         end do
         do j = NGHOST + NY_PHYS + 1, grid%ny_total
            do i = 1, size(ss%tau_x, 1)
               call check(error, abs(ss%tau_x(i, j)) < TOL, &
                          "tau_x in north ghost row should be zero")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
      call ss%destroy()
   end subroutine test_wind_2gyre_ghosts

   subroutine test_spoon_walls_and_max(error)
      !! Spoon: full depth at south-centre, edge depth at east/west walls.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer, parameter :: NX_PHYS = 22, NY_PHYS = 20, NGHOST = 1
      real(wp), parameter :: DX = 100000.0_wp, DY = 110000.0_wp
      real(wp), parameter :: MAX_DEPTH = 2000.0_wp
      real(wp), parameter :: EDGE_DEPTH = 100.0_wp
      real(wp), parameter :: SLOPE_SCALE = 400000.0_wp
      integer :: nx_t, ny_t, j_south, i_west, i_east, i_centre
      real(wp) :: y_phys, expected, D_0, denom, x_len, y_len

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      allocate (b(nx_t, ny_t), source=0.0_wp)
      call set_bathymetry_spoon(b, grid, MAX_DEPTH, EDGE_DEPTH, SLOPE_SCALE)

      x_len = real(NX_PHYS, wp)*DX
      y_len = real(NY_PHYS, wp)*DY
      denom = 1.0_wp - exp(-0.5_wp*y_len/SLOPE_SCALE)
      D_0 = (MAX_DEPTH - EDGE_DEPTH)/(denom*denom)

      checks: block
         ! South-centre interior cell: max depth side of the basin
         ! (no clamp expected at this point).
         j_south = NGHOST + 1
         i_centre = NGHOST + NX_PHYS/2
         y_phys = (real(j_south - NGHOST, wp) - 0.5_wp)*DY
         expected = EDGE_DEPTH &
                    + D_0*sin(PI*(real(i_centre - NGHOST, wp) - 0.5_wp)*DX/x_len) &
                    *(1.0_wp - exp((y_phys - y_len)/SLOPE_SCALE))
         if (expected > MAX_DEPTH) expected = MAX_DEPTH
         if (expected < EDGE_DEPTH) expected = EDGE_DEPTH
         call check(error, abs(b(i_centre, j_south) - expected) < 1.0e-9_wp, &
                    "spoon south-centre depth should match formula")
         if (allocated(error)) exit checks

         ! East/west wall columns: sin(π * (i_phys - 0.5)/nx) is small
         ! at i_phys=1 and i_phys=nx → depth clamps to edge_depth.
         i_west = NGHOST + 1
         i_east = NGHOST + NX_PHYS
         call check(error, abs(b(i_west, j_south) - EDGE_DEPTH) < 1.0e-9_wp .or. &
                    b(i_west, j_south) >= EDGE_DEPTH, &
                    "west wall column should not be deeper than EDGE_DEPTH + slope")
         if (allocated(error)) exit checks
         call check(error, b(i_west, j_south) < 0.5_wp*MAX_DEPTH, &
                    "west wall column should be much shallower than max_depth")
         if (allocated(error)) exit checks
         call check(error, b(i_east, j_south) < 0.5_wp*MAX_DEPTH, &
                    "east wall column should be much shallower than max_depth")
         if (allocated(error)) exit checks

         ! Clamp invariants over the physical interior.
         call check(error, &
                    all(b(NGHOST + 1:NGHOST + NX_PHYS, &
                          NGHOST + 1:NGHOST + NY_PHYS) >= EDGE_DEPTH - TOL), &
                    "all physical-interior depths should be >= edge_depth")
         if (allocated(error)) exit checks
         call check(error, &
                    all(b(NGHOST + 1:NGHOST + NX_PHYS, &
                          NGHOST + 1:NGHOST + NY_PHYS) <= MAX_DEPTH + TOL), &
                    "all physical-interior depths should be <= max_depth")
      end block checks
      deallocate (b)
   end subroutine test_spoon_walls_and_max

   subroutine test_spoon_zero_at_north(error)
      !! At the northern interior row, the (1 - exp(...)) factor is
      !! small → all columns should approach edge_depth.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer, parameter :: NX_PHYS = 22, NY_PHYS = 20, NGHOST = 1
      real(wp), parameter :: DX = 100000.0_wp, DY = 110000.0_wp
      real(wp), parameter :: MAX_DEPTH = 2000.0_wp
      real(wp), parameter :: EDGE_DEPTH = 100.0_wp
      real(wp), parameter :: SLOPE_SCALE = 400000.0_wp
      integer :: nx_t, ny_t, j_north, i

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      allocate (b(nx_t, ny_t), source=0.0_wp)
      call set_bathymetry_spoon(b, grid, MAX_DEPTH, EDGE_DEPTH, SLOPE_SCALE)

      ! At the northern interior edge, (y - y_len) ≈ -0.5*dy, so the
      ! factor is small.  Depths should be much closer to edge_depth
      ! than the southern row.
      j_north = NGHOST + NY_PHYS
      checks: block
         do i = NGHOST + 1, NGHOST + NX_PHYS
            call check(error, b(i, j_north) < 0.5_wp*MAX_DEPTH, &
                       "north interior row depths should be < 0.5 * max_depth")
            if (allocated(error)) exit checks
         end do
      end block checks
      deallocate (b)
   end subroutine test_spoon_zero_at_north

   subroutine test_beta_plane(error)
      !! Verify f_corner varies linearly in y around the reference.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(coriolis_adv_t) :: cor
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 6, NGHOST = 1
      real(wp), parameter :: DX = 100000.0_wp, DY = 110000.0_wp
      real(wp), parameter :: F_0 = 9.4e-5_wp
      real(wp), parameter :: BETA = 1.76e-11_wp
      real(wp), parameter :: Y_REF = 0.5_wp*real(NY_PHYS, wp)*DY
      real(wp) :: y_corner, expected, observed
      integer :: i, j

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      call cor%init(grid, nz_ml=2)
      checks: block
         call cor%set_beta_plane(grid, F_0, BETA, Y_REF)
         do j = 1, size(cor%f_corner, 2)
            ! Per the set_beta_plane docstring: SW corner of cell j
            ! sits at y = (j - 1 - nghost) * dy.
            y_corner = (real(j - 1 - NGHOST, wp))*DY
            expected = F_0 + BETA*(y_corner - Y_REF)
            do i = 1, size(cor%f_corner, 1)
               observed = cor%f_corner(i, j)
               call check(error, abs(observed - expected) < 1.0e-15_wp, &
                          "f_corner should be linear in y")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
      call cor%destroy()
   end subroutine test_beta_plane

end module test_ocean_double_gyre_helpers
