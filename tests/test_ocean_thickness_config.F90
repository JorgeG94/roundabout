!! Analytical tests for `&vcoord_nml thickness_config = "uniform_z"` —
!! the MOM6 `initialize_thickness_uniform` initial-thickness seed.
!!
!! Why the knob exists.  Under the linear density-coordinate IC (MOM6
!! `COORD_CONFIG="linear"`, Roundabout `rho_lightest`/`rho_range`) the layer
!! densities are horizontally uniform, so the layer interfaces ARE the
!! isopycnals.  The legacy `"sigma"` seed (`h_layer = b/nz`) therefore makes
!! every isopycnal follow the bathymetry, standing the whole `rho_range`
!! contrast up across each shelf break at t=0.  On the MOM6 double-gyre spoon
!! (edge_depth = 100 m, max_depth = 2000 m, rho_range = 2 kg/m³) that is a
!! ~1.9 kg/m³ horizontal density jump — g' ≈ 0.018 m/s² — over layers only
!! `edge_depth/nz` thick, which drains the rim layers and drives the day-1
!! grounding blow-up.  `"uniform_z"` seeds flat resting isopycnals instead.
!!
!! The invariant that matters most is the LAST one: `sum(h) == b` exactly.
!! The seed floor must not inject mass the way the runtime `angstrom_h`
!! continuity clamp does.
module test_ocean_thickness_config
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_state, only: seed_h_layer_uniform_z_impl, seed_h_layer_uniform_impl
   implicit none
   private

   public :: collect_ocean_thickness_config_tests

   ! MOM6 double-gyre spoon envelope.
   real(wp), parameter :: MAX_DEPTH = 2000.0_wp
   real(wp), parameter :: EDGE_DEPTH = 100.0_wp
   real(wp), parameter :: ANGSTROM = 1.0e-2_wp
   integer, parameter :: NZ = 100

contains

   subroutine collect_ocean_thickness_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_z_deep_column_is_even_z", test_deep_column), &
                  new_unittest("uniform_z_rim_column_matches_mom6", test_rim_column), &
                  new_unittest("uniform_z_conserves_column_depth", test_conservation), &
                  new_unittest("uniform_z_flattens_isopycnals", test_flat_interfaces), &
                  new_unittest("uniform_z_degenerate_column_floors", test_degenerate), &
                  new_unittest("uniform_z_thick_floor_still_conserves", test_thick_floor) &
                  ]
   end subroutine collect_ocean_thickness_config_tests

   pure function column_sum(h, i, j) result(s)
      !! Sum a single column of the seeded (nx, ny, nz) array.
      real(wp), intent(in) :: h(:, :, :)
      integer, intent(in) :: i, j
      real(wp) :: s
      integer :: k
      s = 0.0_wp
      do k = 1, size(h, 3)
         s = s + h(i, j, k)
      end do
   end function column_sum

   subroutine test_deep_column(error)
      !! A column at exactly `max_depth` gets `nz` layers of `max_depth/nz`:
      !! every z target clears the floor, so nothing collapses.  This is the
      !! branch where "uniform_z" and "sigma" agree, and it must agree exactly.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(1, 1, NZ), h_sigma(1, 1, NZ), b(1, 1)
      real(wp), parameter :: DZ = MAX_DEPTH/real(NZ, wp)
      integer :: k

      b(1, 1) = MAX_DEPTH
      call seed_h_layer_uniform_z_impl(h, b, NZ, MAX_DEPTH, ANGSTROM)
      call seed_h_layer_uniform_impl(h_sigma, b, NZ)

      do k = 1, NZ
         call check(error, abs(h(1, 1, k) - DZ) < 1.0e-10_wp, &
                    "uniform_z: deep column layer is not max_depth/nz")
         if (allocated(error)) return
         ! At b == max_depth the two seeds must coincide.
         call check(error, abs(h(1, 1, k) - h_sigma(1, 1, k)) < 1.0e-10_wp, &
                    "uniform_z: deep column disagrees with the sigma seed")
         if (allocated(error)) return
      end do
   end subroutine test_deep_column

   subroutine test_rim_column(error)
      !! The 100 m spoon rim, hand-derived against MOM6's sweep.  With
      !! max_depth = 2000, nz = 100 the z targets are 20 m apart, so a 100 m
      !! column holds four full 20 m layers plus a partial one, and collapses
      !! the remaining 95 to the floor:
      !!
      !!   k =  1..95  ->  angstrom            (bed side, densest, grounded)
      !!   k = 96      ->  20 - 95*angstrom    (absorbs the collapsed stack)
      !!   k = 97..100 ->  20                  (surface side)
      !!
      !! Bottom-up indexing: k = 1 is the bed, k = nz is the surface.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(1, 1, NZ), b(1, 1)
      integer :: k

      b(1, 1) = EDGE_DEPTH
      call seed_h_layer_uniform_z_impl(h, b, NZ, MAX_DEPTH, ANGSTROM)

      do k = 1, 95
         call check(error, abs(h(1, 1, k) - ANGSTROM) < 1.0e-12_wp, &
                    "uniform_z rim: collapsed layer is not at the angstrom floor")
         if (allocated(error)) return
      end do
      call check(error, abs(h(1, 1, 96) - (20.0_wp - 95.0_wp*ANGSTROM)) < 1.0e-10_wp, &
                 "uniform_z rim: partial layer does not absorb the collapsed stack")
      if (allocated(error)) return
      do k = 97, NZ
         call check(error, abs(h(1, 1, k) - 20.0_wp) < 1.0e-10_wp, &
                    "uniform_z rim: surface-side layer is not 20 m")
         if (allocated(error)) return
      end do
   end subroutine test_rim_column

   subroutine test_conservation(error)
      !! `sum(h) == b` to round-off for every column across the full spoon
      !! depth range, including the shelf break where layers collapse.  The
      !! seed floor must borrow from the partial layer, never inject.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCOL = 40
      real(wp) :: h(NCOL, 1, NZ), b(NCOL, 1)
      integer :: i

      do i = 1, NCOL
         ! Sweep edge_depth -> max_depth.
         b(i, 1) = EDGE_DEPTH + (MAX_DEPTH - EDGE_DEPTH)* &
                   real(i - 1, wp)/real(NCOL - 1, wp)
      end do
      call seed_h_layer_uniform_z_impl(h, b, NZ, MAX_DEPTH, ANGSTROM)

      do i = 1, NCOL
         call check(error, abs(column_sum(h, i, 1) - b(i, 1)) < 1.0e-9_wp, &
                    "uniform_z: column sum does not equal the bathymetry")
         if (allocated(error)) return
         call check(error, minval(h(i, 1, :)) >= ANGSTROM - 1.0e-14_wp, &
                    "uniform_z: layer seeded below the angstrom floor")
         if (allocated(error)) return
      end do
   end subroutine test_conservation

   subroutine test_flat_interfaces(error)
      !! The point of the knob: interfaces are FLAT in z across a shelf break,
      !! whereas the sigma seed tilts them with the bathymetry.  Compare the
      !! depth of one interface (the top of layer k) between a rim column and
      !! a deep column.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(2, 1, NZ), hs(2, 1, NZ), b(2, 1)
      real(wp) :: z_rim, z_deep, z_rim_s, z_deep_s
      integer, parameter :: KTOP = NZ - 1   ! interface one layer below the surface
      integer :: k

      b(1, 1) = EDGE_DEPTH     ! rim
      b(2, 1) = MAX_DEPTH      ! interior
      call seed_h_layer_uniform_z_impl(h, b, NZ, MAX_DEPTH, ANGSTROM)
      call seed_h_layer_uniform_impl(hs, b, NZ)

      ! Interface depth measured DOWN from the surface: sum the layers above it.
      z_rim = 0.0_wp; z_deep = 0.0_wp; z_rim_s = 0.0_wp; z_deep_s = 0.0_wp
      do k = KTOP + 1, NZ
         z_rim = z_rim + h(1, 1, k)
         z_deep = z_deep + h(2, 1, k)
         z_rim_s = z_rim_s + hs(1, 1, k)
         z_deep_s = z_deep_s + hs(2, 1, k)
      end do

      ! uniform_z: same z on both sides of the shelf break => flat isopycnal.
      call check(error, abs(z_rim - z_deep) < 1.0e-9_wp, &
                 "uniform_z: interface is not flat across the shelf break")
      if (allocated(error)) return
      ! sigma: the same interface tilts by the full depth ratio.  Guards the
      ! test itself — if this ever stops holding, the comparison is vacuous.
      call check(error, abs(z_rim_s - z_deep_s) > 1.0_wp, &
                 "sigma seed unexpectedly produced a flat interface")
   end subroutine test_flat_interfaces

   subroutine test_degenerate(error)
      !! A column too shallow to hold `nz` floored layers (or a dry/negative
      !! land bed) falls back to the floored even split rather than producing
      !! negative or zero thicknesses.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(2, 1, NZ), b(2, 1)

      b(1, 1) = 0.5_wp*real(NZ, wp)*ANGSTROM   ! below nz*floor
      b(2, 1) = -3.0_wp                        ! emerged / land
      call seed_h_layer_uniform_z_impl(h, b, NZ, MAX_DEPTH, ANGSTROM)

      call check(error, minval(h(1, 1, :)) >= ANGSTROM - 1.0e-14_wp, &
                 "uniform_z: shallow column seeded below the floor")
      if (allocated(error)) return
      call check(error, minval(h(2, 1, :)) > 0.0_wp, &
                 "uniform_z: dry column seeded a non-positive thickness")
      if (allocated(error)) return
      call check(error, abs(h(2, 1, 1) - ANGSTROM) < 1.0e-14_wp, &
                 "uniform_z: dry column did not fall back to the floor")
   end subroutine test_degenerate

   subroutine test_thick_floor(error)
      !! Guards the room clamp.  When the floor is NOT small compared with
      !! `max_depth/nz`, a naive sweep runs past z = 0 and the telescoping sum
      !! stops equalling `b`.  Here `max_depth/nz = 1 m` against a 0.9 m floor.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZS = 20
      real(wp) :: h(3, 1, NZS), b(3, 1)
      real(wp), parameter :: MD = 20.0_wp     ! => 1 m z targets
      real(wp), parameter :: FLOOR_H = 0.9_wp
      integer :: i

      b(1, 1) = 19.0_wp   ! just under max_depth
      b(2, 1) = 20.0_wp   ! exactly max_depth
      b(3, 1) = 25.0_wp   ! deeper than max_depth (bed layer absorbs excess)
      call seed_h_layer_uniform_z_impl(h, b, NZS, MD, FLOOR_H)

      do i = 1, 3
         call check(error, abs(column_sum(h, i, 1) - b(i, 1)) < 1.0e-9_wp, &
                    "uniform_z: thick floor broke column conservation")
         if (allocated(error)) return
         call check(error, minval(h(i, 1, :)) >= FLOOR_H - 1.0e-12_wp, &
                    "uniform_z: thick floor seeded below the floor")
         if (allocated(error)) return
      end do
      ! Deeper than max_depth: the excess lands in the bed-most layer.
      call check(error, h(3, 1, 1) > 5.0_wp, &
                 "uniform_z: bed layer did not absorb the sub-max_depth excess")
   end subroutine test_thick_floor

end module test_ocean_thickness_config
