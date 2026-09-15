!! Unit tests for the accel_visc_rem kernels (`rdb_ocean_dyn`):
!! the MOM6-parity visc_rem attenuation of the slow explicit velocity
!! applies.  Analytical: `u = u0 + rem·(u − u0)` exact at rem = 0.5,
!! bitwise-inert at rem ≡ 1 (the skip branch), entry-velocity restore
!! at rem = 0, and snapshot exactness.  Arrays follow the chksum
!! mem:separate device-map pattern (`!$acc enter data copyin` before
!! the kernel, `update self` after, `exit data delete` at the end) so
!! the test exercises real device data motion on the GPU build; the
!! directives are inert on host/multicore builds.
module test_ocean_accel_visc_rem
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_dyn, only: accel_visc_rem_snapshot, accel_visc_rem_reweight
   implicit none
   private

   public :: collect_ocean_accel_visc_rem_tests

   integer, parameter :: N1 = 5, N2 = 4, N3 = 3

contains

   subroutine collect_ocean_accel_visc_rem_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("snapshot_exact", test_snapshot_exact), &
                  new_unittest("reweight_half_exact", test_reweight_half_exact), &
                  new_unittest("reweight_rem1_bitwise", test_reweight_rem1_bitwise), &
                  new_unittest("reweight_rem0_restores_entry", test_reweight_rem0) &
                  ]
   end subroutine collect_ocean_accel_visc_rem_tests

   pure subroutine fill_fields(vel0, vel)
      !! Deterministic non-trivial entry velocity + post-apply velocity:
      !! exact in wp (small integers), so every check is equality.
      real(wp), intent(out) :: vel0(N1, N2, N3), vel(N1, N2, N3)
      integer :: i, j, k
      do k = 1, N3
         do j = 1, N2
            do i = 1, N1
               vel0(i, j, k) = real(i + 10*j + 100*k, wp)
               ! Explicit-apply delta of ±(i+j+k): sign-alternating so
               ! the reweight has to preserve direction per face.
               vel(i, j, k) = vel0(i, j, k) &
                              + real((-1)**(i + j + k)*(i + j + k), wp)
            end do
         end do
      end do
   end subroutine fill_fields

   subroutine test_snapshot_exact(error)
      !! snapshot: snap == vel bitwise after the device-side copy.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vel(N1, N2, N3), vel0(N1, N2, N3), snap(N1, N2, N3)
      integer :: i, j, k

      call fill_fields(vel0, vel)
      snap = -777.0_wp
      !$acc enter data copyin(vel, snap)
      call accel_visc_rem_snapshot(N1, N2, N3, vel, snap)
      !$acc update self(snap)
      !$acc exit data delete(vel, snap)

      do k = 1, N3
         do j = 1, N2
            do i = 1, N1
               call check(error, snap(i, j, k) == vel(i, j, k), &
                          "snapshot must copy the velocity bitwise")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_snapshot_exact

   subroutine test_reweight_half_exact(error)
      !! rem = 0.5 everywhere: u = u0 + 0.5·Δ exactly (Δ is an integer
      !! in wp, so 0.5·Δ is exact and the check is equality).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vel0(N1, N2, N3), vel(N1, N2, N3), rem(N1, N2, N3)
      real(wp) :: expect
      integer :: i, j, k

      call fill_fields(vel0, vel)
      rem = 0.5_wp
      !$acc enter data copyin(vel0, vel, rem)
      call accel_visc_rem_reweight(N1, N2, N3, vel0, rem, vel)
      !$acc update self(vel)
      !$acc exit data delete(vel0, vel, rem)

      do k = 1, N3
         do j = 1, N2
            do i = 1, N1
               expect = vel0(i, j, k) &
                        + 0.5_wp*real((-1)**(i + j + k)*(i + j + k), wp)
               call check(error, vel(i, j, k) == expect, &
                          "rem=0.5 must give u0 + 0.5*delta exactly")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_reweight_half_exact

   subroutine test_reweight_rem1_bitwise(error)
      !! rem ≡ 1: the skip branch leaves the applied velocity BITWISE
      !! unchanged (irrational entries — no FP-identity accident: the
      !! rewrite path `u0 + 1·(u − u0)` would NOT reproduce them).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vel0(N1, N2, N3), vel(N1, N2, N3), rem(N1, N2, N3)
      real(wp) :: before(N1, N2, N3)
      integer :: i, j, k

      call fill_fields(vel0, vel)
      ! Perturb into "awkward" FP territory where x0 + (x - x0) /= x.
      vel0 = vel0*1.0e12_wp
      vel = vel + 1.0e-7_wp*acos(-1.0_wp)
      before = vel
      rem = 1.0_wp
      !$acc enter data copyin(vel0, vel, rem)
      call accel_visc_rem_reweight(N1, N2, N3, vel0, rem, vel)
      !$acc update self(vel)
      !$acc exit data delete(vel0, vel, rem)

      do k = 1, N3
         do j = 1, N2
            do i = 1, N1
               call check(error, vel(i, j, k) == before(i, j, k), &
                          "rem=1 must be bitwise inert (skip branch)")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_reweight_rem1_bitwise

   subroutine test_reweight_rem0(error)
      !! rem = 0 at ONE face restores that face's entry velocity exactly
      !! (the friction-dominated-layer limit); every other face carries
      !! rem = 0.25 and must NOT collapse to u0.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vel0(N1, N2, N3), vel(N1, N2, N3), rem(N1, N2, N3)
      real(wp) :: expect
      integer :: i, j, k

      call fill_fields(vel0, vel)
      rem = 0.25_wp
      rem(2, 3, 1) = 0.0_wp
      !$acc enter data copyin(vel0, vel, rem)
      call accel_visc_rem_reweight(N1, N2, N3, vel0, rem, vel)
      !$acc update self(vel)
      !$acc exit data delete(vel0, vel, rem)

      call check(error, vel(2, 3, 1) == vel0(2, 3, 1), &
                 "rem=0 must restore the entry velocity exactly")
      if (allocated(error)) return
      do k = 1, N3
         do j = 1, N2
            do i = 1, N1
               if (i == 2 .and. j == 3 .and. k == 1) cycle
               expect = vel0(i, j, k) &
                        + 0.25_wp*real((-1)**(i + j + k)*(i + j + k), wp)
               call check(error, vel(i, j, k) == expect, &
                          "rem=0.25 faces must get the quarter-strength delta")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_reweight_rem0

end module test_ocean_accel_visc_rem
