!! Unit tests for the MOM6-style per-phase chksum probe
!! (`rdb_ocean_chksum`).  Analytical: exact sum/min/max on a hand-built
!! field, non-finite counting (NaN + Inf), and the step-window gate.
module test_ocean_chksum
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                                               ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_chksum, only: chksum_probe_t, chksum_stats_t, &
                               chksum_active, chksum_stats_3d, chksum_stats_2d, &
                               chksum_argmax, chksum_loc_extents, &
                               LOC_H, LOC_U, LOC_V, LOC_Q
   implicit none
   private

   public :: collect_ocean_chksum_tests

contains

   subroutine collect_ocean_chksum_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("stats_exact_3d", test_stats_exact_3d), &
                  new_unittest("stats_nonfinite_count", test_stats_nonfinite), &
                  new_unittest("window_gate", test_window_gate), &
                  new_unittest("hotface_argmax_interior", test_argmax_interior), &
                  new_unittest("bits_reorder_invariant", test_bits_reorder_invariant), &
                  new_unittest("loc_extents", test_loc_extents) &
                  ]
   end subroutine collect_ocean_chksum_tests

   subroutine test_stats_exact_3d(error)
      !! arr(i,j,k) = i + 10j + 100k on 4x3x2: sum/min/max are exact
      !! integers in wp, so the checks are equality, not tolerance.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: arr(4, 3, 2)
      type(chksum_stats_t) :: st
      integer :: i, j, k
      real(wp) :: expect_sum

      expect_sum = 0.0_wp
      do k = 1, 2
         do j = 1, 3
            do i = 1, 4
               arr(i, j, k) = real(i + 10*j + 100*k, wp)
               expect_sum = expect_sum + arr(i, j, k)
            end do
         end do
      end do
      !$acc enter data copyin(arr)
      call chksum_stats_3d(arr, 4, 3, 2, st)
      !$acc exit data delete(arr)

      call check(error, st%total == expect_sum, "3D sum must be exact")
      if (allocated(error)) return
      call check(error, st%minv == 111.0_wp, "3D min must be arr(1,1,1)")
      if (allocated(error)) return
      call check(error, st%maxv == 234.0_wp, "3D max must be arr(4,3,2)")
      if (allocated(error)) return
      call check(error, st%nonfin == 0, "finite field must count 0 nonfin")
   end subroutine test_stats_exact_3d

   subroutine test_stats_nonfinite(error)
      !! Seed one NaN + two Infs: nonfin == 3 (the corruption signal —
      !! min/max are NaN-blind, sum goes NaN; only nonfin is load-bearing).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: arr(5, 4)
      type(chksum_stats_t) :: st

      arr = 1.0_wp
      arr(2, 2) = ieee_value(1.0_wp, ieee_quiet_nan)
      arr(3, 3) = ieee_value(1.0_wp, ieee_positive_inf)
      arr(4, 1) = -ieee_value(1.0_wp, ieee_positive_inf)
      !$acc enter data copyin(arr)
      call chksum_stats_2d(arr, 5, 4, st)
      !$acc exit data delete(arr)

      call check(error, st%nonfin == 3, &
                 "nonfin must count NaN + both Infs (the corruption signal)")
   end subroutine test_stats_nonfinite

   subroutine test_window_gate(error)
      !! chksum_active: enable gate + [start, end] window semantics,
      !! including the 0 = unbounded conventions.
      type(error_type), allocatable, intent(out) :: error
      type(chksum_probe_t) :: p

      call check(error,.not. chksum_active(p, 1), "default-off must gate")
      if (allocated(error)) return
      p%enable = .true.
      call check(error, chksum_active(p, 1), "enabled + no window = always on")
      if (allocated(error)) return
      p%start_step = 10
      p%end_step = 20
      call check(error,.not. chksum_active(p, 9), "step below window must gate")
      if (allocated(error)) return
      call check(error, chksum_active(p, 10) .and. chksum_active(p, 20), &
                 "window bounds are inclusive")
      if (allocated(error)) return
      call check(error,.not. chksum_active(p, 21), "step above window must gate")
      if (allocated(error)) return
      p%end_step = 0
      call check(error, chksum_active(p, 99999), "end_step 0 = unbounded")
   end subroutine test_window_gate

   subroutine test_argmax_interior(error)
      !! HOTFACE argmax core: |·| semantics (a large NEGATIVE interior
      !! value wins over smaller positives), and the caller's interior
      !! bounds exclude a LARGER value planted in the ghost ring.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: arr(7, 6, 3)
      integer :: im, jm, km
      real(wp) :: amax

      arr = 0.5_wp
      arr(4, 3, 2) = -9.0_wp     ! interior max by magnitude, negative
      arr(6, 4, 1) = 7.0_wp      ! interior runner-up
      arr(1, 1, 3) = 99.0_wp     ! ghost-ring decoy — bounds must exclude
      arr(7, 6, 2) = -88.0_wp    ! ghost-ring decoy, negative side

      call chksum_argmax(arr, 7, 6, 3, 2, 6, 2, 5, im, jm, km, amax)

      call check(error, amax == 9.0_wp, "argmax must report |max|, sign-blind")
      if (allocated(error)) return
      call check(error, im == 4 .and. jm == 3 .and. km == 2, &
                 "argmax must land on the interior magnitude max, not the ghost decoys")
   end subroutine test_argmax_interior

   subroutine test_bits_reorder_invariant(error)
      !! The core guarantee of the decomposition-invariant checksum: the
      !! `bits` POPCNT reduction is IDENTICAL when the SAME multiset of
      !! values is visited in a different order (integer addition is
      !! exactly associative + commutative), whereas the FP `total` is
      !! order-sensitive — here a huge magnitude next to tiny ones makes
      !! the forward and reversed sums round differently.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 4*3*2
      real(wp) :: vals(n), arr(4, 3, 2), rev(4, 3, 2)
      type(chksum_stats_t) :: sf, sr
      integer :: i, j, k, p

      ! A spread of magnitudes so the FP sum is genuinely reorder-sensitive.
      vals = [1.0e16_wp, 1.0_wp, -1.0_wp, 3.0_wp, -7.5_wp, 2.0e-8_wp, &
              1.25_wp, -4.0e15_wp, 8.0_wp, 0.5_wp, -0.25_wp, 6.0_wp, &
              9.0e14_wp, -2.0_wp, 42.0_wp, 1.0e-3_wp, -1.0e12_wp, 5.0_wp, &
              7.0_wp, -3.5_wp, 11.0_wp, -0.125_wp, 4.0_wp, 100.0_wp]
      p = 0
      do k = 1, 2
         do j = 1, 3
            do i = 1, 4
               p = p + 1
               arr(i, j, k) = vals(p)        ! forward fill
               rev(i, j, k) = vals(n + 1 - p) ! same multiset, reversed order
            end do
         end do
      end do

      !$acc enter data copyin(arr, rev)
      call chksum_stats_3d(arr, 4, 3, 2, sf)
      call chksum_stats_3d(rev, 4, 3, 2, sr)
      !$acc exit data delete(arr, rev)

      ! Guard against a no-op permutation: the two arrangements must really
      ! differ (a huge value where a small one now sits).  The FP `total`
      ! NEED NOT match across the two orders (documented, not asserted —
      ! its rounding is compiler-dependent); only `bits` is invariant.
      call check(error, arr(1, 1, 1) /= rev(1, 1, 1), &
                 "the two fills must be a genuine reordering")
      if (allocated(error)) return
      call check(error, sf%bits == sr%bits, &
                 "bits POPCNT reduction must be order-invariant (decomp-invariant)")
   end subroutine test_bits_reorder_invariant

   subroutine test_loc_extents(error)
      !! Grid-location API bounds: each LOC tag selects the correct
      !! Arakawa-C extents from the grid — u/v faces carry the extra
      !! wall-normal row/column, the corner both.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nx, ny

      grid%nx_total = 10
      grid%ny_total = 7

      call chksum_loc_extents(grid, LOC_H, nx, ny)
      call check(error, nx == 10 .and. ny == 7, "LOC_H = cell-centre (nx, ny)")
      if (allocated(error)) return
      call chksum_loc_extents(grid, LOC_U, nx, ny)
      call check(error, nx == 11 .and. ny == 7, "LOC_U = (nx+1, ny)")
      if (allocated(error)) return
      call chksum_loc_extents(grid, LOC_V, nx, ny)
      call check(error, nx == 10 .and. ny == 8, "LOC_V = (nx, ny+1)")
      if (allocated(error)) return
      call chksum_loc_extents(grid, LOC_Q, nx, ny)
      call check(error, nx == 11 .and. ny == 8, "LOC_Q = (nx+1, ny+1)")
   end subroutine test_loc_extents

end module test_ocean_chksum
