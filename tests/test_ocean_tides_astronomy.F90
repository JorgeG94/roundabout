!! Unit tests for the equilibrium-tide astronomy generator
!! (`rdb_ocean_tide_astro`).  Asserts the mean longitudes, equilibrium
!! arguments V_c (incl. the load-bearing ±pi/2 diurnal signs), the nodal
!! f/u corrections, and the Gregorian day-number helper against the
!! clean-room golden numbers (Kowalik-Luick / Schureman), computed by
!! the module and cross-checked in `tools/prototype/tide_astro.py`.
module test_ocean_tides_astronomy
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_tide_astro, only: gregorian_day_number, days_since_1900, &
                                   mean_longitudes, equilibrium_arguments, &
                                   nodal_fu, tide_name_index, TIDE_DEG2RAD, &
                                   TIDES_CATALOG_SIZE
   implicit none
   private

   public :: collect_ocean_tides_astronomy_tests

   real(wp), parameter :: TOL_DEG = 1.0e-5_wp
   real(wp), parameter :: TOL_F = 1.0e-6_wp

contains

   subroutine collect_ocean_tides_astronomy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("day_number_helper", test_day_number), &
                  new_unittest("mean_longitudes_1900", test_mean_long_1900), &
                  new_unittest("mean_longitudes_2000", test_mean_long_2000), &
                  new_unittest("equilibrium_arguments_1900", test_v_1900), &
                  new_unittest("diurnal_pm_half_pi_signs", test_diurnal_signs), &
                  new_unittest("nodal_fu_m2_1900", test_nodal_m2), &
                  new_unittest("nodal_off_is_identity", test_nodal_off), &
                  new_unittest("name_index_lookup", test_name_index)]
   end subroutine collect_ocean_tides_astronomy_tests

   pure function vdeg(v_rad) result(d)
      !! Equilibrium argument in canonical degrees [0,360).
      real(wp), intent(in) :: v_rad
      real(wp) :: d
      d = mod(v_rad/TIDE_DEG2RAD, 360.0_wp)
      if (d < 0.0_wp) d = d + 360.0_wp
   end function vdeg

   subroutine test_day_number(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, days_since_1900(1900, 1, 1), 0.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
      ! 100 yrs x 365 + 24 leap days (1900 is NOT a leap year) = 36524.
      call check(error, days_since_1900(2000, 1, 1), 36524.0_wp, thr=1.0e-9_wp)
      if (allocated(error)) return
      ! One ordinary day forward.
      call check(error, days_since_1900(1900, 1, 2), 1.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
      ! Known JDN: 2000-01-01 (proleptic Gregorian) = 2451545.
      call check(error, real(gregorian_day_number(2000, 1, 1), wp), &
                 2451545.0_wp, thr=1.0e-6_wp)
   end subroutine test_day_number

   subroutine test_mean_long_1900(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: s, h, p, n
      call mean_longitudes(0.0_wp, s, h, p, n)
      call check(error, s, 277.0248_wp, thr=TOL_DEG)
      if (allocated(error)) return
      call check(error, h, 280.1895_wp, thr=TOL_DEG)
      if (allocated(error)) return
      call check(error, p, 334.3853_wp, thr=TOL_DEG)
      if (allocated(error)) return
      call check(error, n, 259.1568_wp, thr=TOL_DEG)
   end subroutine test_mean_long_1900

   subroutine test_mean_long_2000(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: s, h, p, n
      ! D = 36524.  N raw is ≈ -1674.93°; the module folds to [0,360) so we
      ! assert the canonical residue directly (no bare-mod sign trap).
      call mean_longitudes(36524.0_wp, s, h, p, n)
      call check(error, s, 211.740103_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, h, 279.973056_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, p, 83.297596_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, n, 125.069854_wp, thr=1.0e-4_wp)
   end subroutine test_mean_long_2000

   subroutine test_v_1900(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(TIDES_CATALOG_SIZE)
      call equilibrium_arguments(0.0_wp, v)
      call check(error, vdeg(v(1)), 6.32940_wp, thr=1.0e-4_wp)
      if (allocated(error)) return   ! M2
      call check(error, vdeg(v(2)), 0.0_wp, thr=1.0e-6_wp)
      if (allocated(error)) return       ! S2
      call check(error, vdeg(v(3)), 63.68990_wp, thr=1.0e-4_wp)
      if (allocated(error)) return  ! N2
      call check(error, vdeg(v(4)), 200.37900_wp, thr=1.0e-3_wp)
      if (allocated(error)) return  ! K2
      call check(error, vdeg(v(5)), 10.18950_wp, thr=1.0e-4_wp)
      if (allocated(error)) return  ! K1
      call check(error, vdeg(v(6)), 356.13990_wp, thr=1.0e-3_wp)
      if (allocated(error)) return  ! O1
      call check(error, vdeg(v(7)), 349.81050_wp, thr=1.0e-3_wp)
      if (allocated(error)) return  ! P1
      call check(error, vdeg(v(8)), 53.50040_wp, thr=1.0e-4_wp)
      if (allocated(error)) return  ! Q1
      call check(error, vdeg(v(9)), 194.04960_wp, thr=1.0e-3_wp)
      if (allocated(error)) return  ! MF
      call check(error, vdeg(v(10)), 302.63950_wp, thr=1.0e-3_wp)                              ! MM
   end subroutine test_v_1900

   subroutine test_diurnal_signs(error)
      !! K1 carries +pi/2, O1/P1/Q1 carry -pi/2.  Isolate the sign by
      !! subtracting the longitude combination and reading the residual.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: v(TIDES_CATALOG_SIZE)
      real(wp) :: s, h, p, n, srad, hrad, prad
      call equilibrium_arguments(0.0_wp, v)
      call mean_longitudes(0.0_wp, s, h, p, n)
      srad = s*TIDE_DEG2RAD
      hrad = h*TIDE_DEG2RAD
      prad = p*TIDE_DEG2RAD
      call check(error, (v(5) - hrad)/TIDE_DEG2RAD, 90.0_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, (v(6) - (-2.0_wp*srad + hrad))/TIDE_DEG2RAD, -90.0_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, (v(7) - (-hrad))/TIDE_DEG2RAD, -90.0_wp, thr=1.0e-4_wp)
      if (allocated(error)) return
      call check(error, (v(8) - (-3.0_wp*srad + hrad + prad))/TIDE_DEG2RAD, -90.0_wp, thr=1.0e-4_wp)
   end subroutine test_diurnal_signs

   subroutine test_nodal_m2(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: f(TIDES_CATALOG_SIZE), u(TIDES_CATALOG_SIZE)
      call nodal_fu(0.0_wp, .true., f, u)
      ! M2 golden: f=1.00696051, u=+3.59975196e-2 rad (+2.06250595 deg).
      call check(error, f(1), 1.00696051_wp, thr=TOL_F)
      if (allocated(error)) return
      call check(error, u(1), 3.59975196e-2_wp, thr=1.0e-7_wp)
      if (allocated(error)) return
      ! O1 carries the +10.8 deg coefficient, K1 the -8.9 deg coefficient:
      ! their nodal phases must have OPPOSITE signs (the sin(N) factor is
      ! shared, so the coefficient sign difference shows through as a sign
      ! flip between them — independent of the epoch's sin(N)).
      call check(error, u(5)*u(6) < 0.0_wp, "K1/O1 nodal phases must be opposite-signed")
      if (allocated(error)) return
      ! Concrete values at 1900 (sin(N) < 0): O1 u < 0, K1 u > 0.
      call check(error, u(6) < 0.0_wp, "O1 nodal u sign wrong at 1900")
      if (allocated(error)) return
      call check(error, u(5) > 0.0_wp, "K1 nodal u sign wrong at 1900")
   end subroutine test_nodal_m2

   subroutine test_nodal_off(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: f(TIDES_CATALOG_SIZE), u(TIDES_CATALOG_SIZE)
      integer :: c
      call nodal_fu(0.0_wp, .false., f, u)
      do c = 1, TIDES_CATALOG_SIZE
         call check(error, f(c), 1.0_wp, thr=1.0e-12_wp)
         if (allocated(error)) return
         call check(error, u(c), 0.0_wp, thr=1.0e-12_wp)
         if (allocated(error)) return
      end do
   end subroutine test_nodal_off

   subroutine test_name_index(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, tide_name_index("M2"), 1)
      if (allocated(error)) return
      call check(error, tide_name_index("m2"), 1)
      if (allocated(error)) return
      call check(error, tide_name_index("Q1"), 8)
      if (allocated(error)) return
      call check(error, tide_name_index("MM"), 10)
      if (allocated(error)) return
      call check(error, tide_name_index("ZZ"), -1)
   end subroutine test_name_index

end module test_ocean_tides_astronomy
