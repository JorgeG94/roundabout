!! Unit tests for the fused diagnostic statistics reduction,
!! `rdb_ocean_diag::diag_reduce_stats`.
!!
!! The emit path of `ocean_diag_t%step` used to make three host passes over
!! every diagnostic buffer (`minval` / `maxval` / `sum`).  Those are now one
!! `do concurrent(k, j, i) reduce(min:) reduce(max:) reduce(+:)`, gated on
!! `ocean_diag_t%on_device` so it only runs where the buffer actually lives.
!! This file pins down both halves of that contract:
!!
!!   * NUMERICS — the fused reduction agrees with the intrinsics it replaced,
!!     across ordinary and degenerate shapes and across every sign pattern.
!!   * RESIDENCY — the reduction really does read DEVICE memory on a GPU
!!     build, not a stale host copy.
!!
!! Exactness policy, applied consistently below:
!!   * `min` / `max` are associative AND commutative, so no partitioning of the
!!     index space can change the answer.  They are asserted BIT-IDENTICAL to
!!     `minval` / `maxval`; a mismatch is a miscompile or a logic bug, never
!!     round-off.
!!   * `+` on floating point is NOT associative.  A parallel reduction sums a
!!     different tree than the intrinsic, so the result may differ in the last
!!     few ulps (the commit that introduced the kernel measured exactly this:
!!     an emitted mean moving -2.23986E-17 -> -2.23987E-17 on a
!!     cancellation-dominated field).  The general sum is therefore compared
!!     within a relative tolerance scaled by `sum(abs(...))`, which stays
!!     meaningful even when the signed sum cancels to near zero.
!!   * The one exception is `reduce_sum_exact_on_integer_valued`: for
!!     integer-valued data whose partial sums never leave the exactly
!!     representable range (|x| < 2**53), EVERY addition is exact whatever the
!!     order, so the sum is asserted bit-identical there too.  That case would
!!     catch a dropped or double-counted element that a tolerance hides.
module test_ocean_diag_reduce
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_diag, only: ocean_diag_t, diag_reduce_stats, diag_field_stats, &
                             DIAG_MISSING_VALUE
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_diag_reduce_tests

   !! Relative tolerance for the reassociated sum, scaled by sum(abs(buf)).
   !! Generous by ulp standards (~1e-4 of it) and still ~4 orders tighter than
   !! any realistic bug in a whole-array reduction.
   real(wp), parameter :: SUM_RTOL = 1.0e-12_wp

   !! Host-poison sentinel for the device-residency tests.  Chosen well outside
   !! the range of every pattern used here so a host read is unambiguous.
   real(wp), parameter :: POISON = -7.77e5_wp

contains

   subroutine collect_ocean_diag_reduce_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("reduce_matches_intrinsics_over_shapes", test_matches_intrinsics), &
                  new_unittest("reduce_sum_exact_on_integer_valued", test_sum_exact), &
                  new_unittest("reduce_n3_one_is_2d_diagnostic", test_n3_one), &
                  new_unittest("reduce_single_element", test_single_element), &
                  new_unittest("reduce_all_equal_values", test_all_equal), &
                  new_unittest("reduce_all_negative_data", test_all_negative), &
                  new_unittest("reduce_all_positive_data", test_all_positive), &
                  new_unittest("reduce_mixed_sign_data", test_mixed_sign), &
                  new_unittest("reduce_skips_nan_missing_cells", test_reduce_skips_nan), &
                  new_unittest("stats_mean_excludes_nan_missing_cells", test_stats_skips_nan), &
                  new_unittest("stats_all_finite_matches_intrinsics", test_stats_all_finite), &
                  new_unittest("stats_all_missing_reports_sentinel", test_stats_all_missing), &
                  new_unittest("device_resident_buffer_is_read_on_device", test_device_resident), &
                  new_unittest("device_resident_all_negative_buffer", test_device_resident_negative), &
                  new_unittest("device_resident_nan_cells_are_skipped", test_device_resident_nan), &
                  new_unittest("device_resident_stats_shim_stays_on_device", test_device_resident_stats), &
                  new_unittest("on_device_flag_lifecycle", test_on_device_flag) &
                  ]
   end subroutine collect_ocean_diag_reduce_tests

   ! ---------------------------------------------------------------------
   ! Data patterns
   ! ---------------------------------------------------------------------

   pure function pattern(i, j, k) result(v)
      !! Deterministic, mixed-sign, NON-integer-valued pseudo-random pattern in
      !! roughly [-143, +143].  Non-integer on purpose: the fractional mantissa
      !! is what makes the reassociated sum differ from `sum()` at all, so the
      !! tolerance branch is genuinely exercised rather than accidentally exact.
      integer, intent(in) :: i, j, k
      real(wp) :: v
      v = real(mod(i*7919 + j*104729 + k*1299709, 2003), wp)/7.0_wp - 143.0_wp
   end function pattern

   pure function pattern_int(i, j, k) result(v)
      !! Integer-valued in [-500, 499].  Every partial sum stays an exact
      !! float64 integer, so the reduction result is order-independent.
      integer, intent(in) :: i, j, k
      real(wp) :: v
      v = real(mod(i*7919 + j*631 + k*97, 1000) - 500, wp)
   end function pattern_int

   subroutine fill_pattern(buf, use_int, shift)
      !! Stamp `buf` from one of the two patterns, optionally shifted so the
      !! whole array lands on one side of zero.  Plain sequential loops (never
      !! `do concurrent`) so this is unambiguously host code that runs before
      !! any device mapping.
      ! assumed-shape-ok: host-only setup helper, no `do concurrent` here.
      real(wp), intent(out) :: buf(:, :, :)
      logical, intent(in) :: use_int
      real(wp), intent(in) :: shift
      integer :: i, j, k
      do k = 1, size(buf, 3)
         do j = 1, size(buf, 2)
            do i = 1, size(buf, 1)
               if (use_int) then
                  buf(i, j, k) = pattern_int(i, j, k) + shift
               else
                  buf(i, j, k) = pattern(i, j, k) + shift
               end if
            end do
         end do
      end do
   end subroutine fill_pattern

   ! ---------------------------------------------------------------------
   ! Assertion helper
   ! ---------------------------------------------------------------------

   subroutine check_against_intrinsics(error, buf, tag)
      !! Run `diag_reduce_stats` over `buf` and compare to the intrinsics it
      !! replaced.  min/max exact, sum within SUM_RTOL * sum(abs(buf)).
      ! assumed-shape-ok: host-only assertion helper, no `do concurrent` here.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: buf(:, :, :)
      character(len=*), intent(in) :: tag
      real(wp) :: vmin, vmax, vsum, tol

      call diag_reduce_stats(buf, size(buf, 1), size(buf, 2), size(buf, 3), &
                             vmin, vmax, vsum)

      checks: block
         ! min/max are associative + commutative => order-invariant => exact.
         call check(error, vmin == minval(buf), &
                    "["//tag//"] vmin must be BIT-IDENTICAL to minval (min is associative)")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), &
                    "["//tag//"] vmax must be BIT-IDENTICAL to maxval (max is associative)")
         if (allocated(error)) exit checks
         ! Floating-point + is not associative: allow the reassociation ulps.
         tol = SUM_RTOL*sum(abs(buf))
         call check(error, abs(vsum - sum(buf)) <= tol, &
                    "["//tag//"] vsum must match sum() to within the reassociation tolerance")
         if (allocated(error)) exit checks
      end block checks
   end subroutine check_against_intrinsics

   ! ---------------------------------------------------------------------
   ! Numerics
   ! ---------------------------------------------------------------------

   subroutine test_matches_intrinsics(error)
      !! Several buffer shapes, including flat and pencil-shaped ones, against
      !! the intrinsics.  Catches an index-order or bounds slip in the kernel
      !! (a transposed `(i, j, k)` walk over a non-cubic buffer reads out of
      !! range or misses elements).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NSHAPE = 6
      integer, parameter :: S1(NSHAPE) = [3, 8, 17, 32, 1, 64]
      integer, parameter :: S2(NSHAPE) = [5, 4, 2, 16, 9, 1]
      integer, parameter :: S3(NSHAPE) = [7, 3, 11, 8, 1, 2]
      real(wp), allocatable :: buf(:, :, :)
      character(len=32) :: tag
      integer :: s

      do s = 1, NSHAPE
         allocate (buf(S1(s), S2(s), S3(s)))
         call fill_pattern(buf, use_int=.false., shift=0.0_wp)
         write (tag, "(A,I0,A,I0,A,I0)") "shape ", S1(s), "x", S2(s), "x", S3(s)
         call check_against_intrinsics(error, buf, trim(tag))
         deallocate (buf)
         if (allocated(error)) return
      end do
   end subroutine test_matches_intrinsics

   subroutine test_sum_exact(error)
      !! Integer-valued data: the sum is order-independent, so demand exact
      !! equality.  This is the assertion that would catch an element visited
      !! twice or not at all, which SUM_RTOL would happily absorb.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 13, N2 = 6, N3 = 5
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum

      call fill_pattern(buf, use_int=.true., shift=0.0_wp)
      call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)

      checks: block
         call check(error, vmin == minval(buf), "integer-valued vmin must be exact")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), "integer-valued vmax must be exact")
         if (allocated(error)) exit checks
         ! |partial sums| <= 500 * 390 << 2**53, so every addition is exact in
         ! float64 regardless of the association order the reduction picks.
         call check(error, vsum == sum(buf), &
                    "integer-valued vsum must be BIT-IDENTICAL to sum() (all additions exact)")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_sum_exact

   subroutine test_n3_one(error)
      !! A 2-D diagnostic (SSH, mixed-layer depth, ...) is registered with
      !! `n3 = 1`.  The innermost `do concurrent` extent is then a single
      !! iteration, which is exactly where an off-by-one collapses to zero
      !! trips and leaves the seeds `+huge` / `-huge` / `0` untouched.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 11, N2 = 7
      real(wp) :: buf(N1, N2, 1)
      real(wp) :: vmin, vmax, vsum

      call fill_pattern(buf, use_int=.false., shift=0.0_wp)
      call diag_reduce_stats(buf, N1, N2, 1, vmin, vmax, vsum)

      checks: block
         call check(error, vmin == minval(buf), "n3=1 vmin must equal minval")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), "n3=1 vmax must equal maxval")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - sum(buf)) <= SUM_RTOL*sum(abs(buf)), &
                    "n3=1 vsum must match sum()")
         if (allocated(error)) exit checks
         ! Explicitly rule out "loop never ran, seeds survived".
         call check(error, vmin < huge(1.0_wp), "n3=1 vmin still at the +huge seed")
         if (allocated(error)) exit checks
         call check(error, vmax > -huge(1.0_wp), "n3=1 vmax still at the -huge seed")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_n3_one

   subroutine test_single_element(error)
      !! 1x1x1 — the smallest legal buffer.  min, max and sum must all be the
      !! single value, which pins the seeds down completely.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: V = -3.25_wp
      real(wp) :: buf(1, 1, 1)
      real(wp) :: vmin, vmax, vsum

      buf(1, 1, 1) = V
      call diag_reduce_stats(buf, 1, 1, 1, vmin, vmax, vsum)

      checks: block
         call check(error, vmin == V, "1x1x1 vmin must be the single element")
         if (allocated(error)) exit checks
         call check(error, vmax == V, "1x1x1 vmax must be the single element")
         if (allocated(error)) exit checks
         call check(error, vsum == V, "1x1x1 vsum must be the single element")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_single_element

   subroutine test_all_equal(error)
      !! A constant field — what a freshly filled or masked-off diagnostic
      !! looks like.  min == max == the value, and the sum is exact because
      !! 2.5 * 2**m is representable for the counts used here.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 8, N2 = 4, N3 = 2
      real(wp), parameter :: V = 2.5_wp
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum

      buf = V
      call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)

      checks: block
         call check(error, vmin == V, "constant field vmin must be the constant")
         if (allocated(error)) exit checks
         call check(error, vmax == V, "constant field vmax must be the constant")
         if (allocated(error)) exit checks
         call check(error, vsum == V*real(N1*N2*N3, wp), &
                    "constant field vsum must be exactly n*value")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_all_equal

   subroutine test_all_negative(error)
      !! Every element strictly negative.  Seeding `vmax` with `0` instead of
      !! `-huge` returns 0 here — a value that is not in the buffer at all —
      !! so this test is the guard on the max seed.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 9, N2 = 5, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum

      ! pattern() spans [-143, +143]; shift down by 200 to clear zero.
      call fill_pattern(buf, use_int=.false., shift=-200.0_wp)

      checks: block
         call check(error, maxval(buf) < 0.0_wp, "test setup: buffer must be all-negative")
         if (allocated(error)) exit checks

         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)

         call check(error, vmax == maxval(buf), &
                    "all-negative vmax must equal maxval (a 0 seed would leak through)")
         if (allocated(error)) exit checks
         call check(error, vmax < 0.0_wp, "all-negative vmax must itself be negative")
         if (allocated(error)) exit checks
         call check(error, vmin == minval(buf), "all-negative vmin must equal minval")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - sum(buf)) <= SUM_RTOL*sum(abs(buf)), &
                    "all-negative vsum must match sum()")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_all_negative

   subroutine test_all_positive(error)
      !! Mirror image: every element strictly positive, so a `0` seed on
      !! `vmin` would leak a value that is not in the buffer.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 9, N2 = 5, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum

      call fill_pattern(buf, use_int=.false., shift=200.0_wp)

      checks: block
         call check(error, minval(buf) > 0.0_wp, "test setup: buffer must be all-positive")
         if (allocated(error)) exit checks

         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)

         call check(error, vmin == minval(buf), &
                    "all-positive vmin must equal minval (a 0 seed would leak through)")
         if (allocated(error)) exit checks
         call check(error, vmin > 0.0_wp, "all-positive vmin must itself be positive")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), "all-positive vmax must equal maxval")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - sum(buf)) <= SUM_RTOL*sum(abs(buf)), &
                    "all-positive vsum must match sum()")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_all_positive

   subroutine test_mixed_sign(error)
      !! Straddles zero — the ordinary case for u, v, w, SSH anomaly.  Both
      !! extremes must be strictly on their own side of zero, so neither seed
      !! can be `0` and pass.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 12, N2 = 6, N3 = 4
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum

      call fill_pattern(buf, use_int=.false., shift=0.0_wp)

      checks: block
         call check(error, minval(buf) < 0.0_wp .and. maxval(buf) > 0.0_wp, &
                    "test setup: buffer must straddle zero")
         if (allocated(error)) exit checks

         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)

         call check(error, vmin == minval(buf), "mixed-sign vmin must equal minval")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), "mixed-sign vmax must equal maxval")
         if (allocated(error)) exit checks
         call check(error, vmin < 0.0_wp .and. vmax > 0.0_wp, &
                    "mixed-sign extremes must land on opposite sides of zero")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - sum(buf)) <= SUM_RTOL*sum(abs(buf)), &
                    "mixed-sign vsum must match sum()")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_mixed_sign

   ! ---------------------------------------------------------------------
   ! Device residency
   !
   ! `-gpu=..., mem:separate` — the GPU build has NO managed/unified memory,
   ! so host and device memory are two distinct copies and the only thing that
   ! ties them together is an explicit map.  A device test that merely maps a
   ! buffer, reduces it, and compares to the host values proves NOTHING on its
   ! own: delete the `enter data` and it still passes on every build, because
   ! the two copies were identical the whole time.  It would be a host test
   ! wearing a device costume.
   !
   ! The construction below makes the two copies DISAGREE.  The buffer is
   ! copied in, and then, ON A GPU BUILD ONLY, the HOST copy is overwritten
   ! with POISON:
   !
   !     device : the real pattern      host : POISON everywhere
   !
   ! The assertions demand the pattern's statistics, so only a reduction that
   ! read the DEVICE copy can pass.  Concretely, for this kernel `-Minfo`
   ! reports
   !
   !     Generating implicit copyin(buf(:n1,:n2,:n3)) [if not already present]
   !
   ! — i.e. with the map in place the kernel binds to the resident device copy,
   ! and without it nvfortran implicitly stages the HOST buffer across, which
   ! is now POISON.  VERIFIED by deleting the two directives below and
   ! rebuilding: exactly these two tests go red (`vmin must be the mapped
   ! device data, not the host copy`) while the nine host-side numerics tests
   ! stay green.  Whether an unmapped buffer degrades to an implicit stage-in
   ! or to raw garbage is a compiler-version detail; either way the answer is
   ! not the device copy, and either way these assertions catch it.
   !
   ! The poison is `#ifdef RDB_GPU_OFFLOAD` because that macro is defined
   ! exactly when the `!$acc enter data` above it is a real device transfer
   ! (cmake/compiler_flags.cmake).  On a host / multicore build the directive
   ! is an inert comment, there is only one copy of the buffer, and poisoning
   ! it would be poisoning the data under test.  The `!$acc exit data delete`
   ! is mandatory, not tidiness: `buf` is a stack array, and leaving a present
   ! -table entry pointing at a dead frame corrupts a later test in the same
   ! binary.
   ! ---------------------------------------------------------------------

   subroutine test_device_resident(error)
      !! Mixed-sign buffer, reduced where it lives.  See the block comment
      !! above for why failure to map produces a wrong answer rather than a
      !! silent pass.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 16, N2 = 8, N3 = 5
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: ref(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs

      call fill_pattern(buf, use_int=.false., shift=0.0_wp)
      ! `ref` is a plain host array that is never mapped and never poisoned:
      ! the reference the device answer is judged against.
      ref = buf
      ref_min = minval(ref)
      ref_max = maxval(ref)
      ref_sum = sum(ref)
      ref_abs = sum(abs(ref))

      checks: block
         ! The discriminator only has teeth if POISON is outside the data.
         call check(error, POISON < ref_min, &
                    "test setup: POISON must sit below the buffer minimum")
         if (allocated(error)) exit checks

         !$acc enter data copyin(buf)
#ifdef RDB_GPU_OFFLOAD
         call poison_host(buf)
#endif
         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)
         !$acc exit data delete(buf)

         call check(error, vmin == ref_min, &
                    "device reduction vmin must be the mapped device data, not the host copy")
         if (allocated(error)) exit checks
         call check(error, vmax == ref_max, &
                    "device reduction vmax must be the mapped device data, not the host copy")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - ref_sum) <= SUM_RTOL*ref_abs, &
                    "device reduction vsum must be the mapped device data, not the host copy")
         if (allocated(error)) exit checks
         ! Belt and braces: name the failure mode explicitly.
         call check(error, vmin /= POISON, &
                    "device reduction returned POISON -- it read the HOST buffer")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_device_resident

   subroutine test_device_resident_negative(error)
      !! Same contract, all-negative data: `vmax` must come back negative.  On
      !! the device the seeds are set inside the kernel prologue rather than by
      !! a host store, so this re-checks the `-huge` seed in the place it
      !! actually matters.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 10, N2 = 6, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: ref(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs

      call fill_pattern(buf, use_int=.false., shift=-200.0_wp)
      ref = buf
      ref_min = minval(ref)
      ref_max = maxval(ref)
      ref_sum = sum(ref)
      ref_abs = sum(abs(ref))

      checks: block
         call check(error, ref_max < 0.0_wp, "test setup: buffer must be all-negative")
         if (allocated(error)) exit checks
         call check(error, POISON < ref_min, &
                    "test setup: POISON must sit below the buffer minimum")
         if (allocated(error)) exit checks

         !$acc enter data copyin(buf)
#ifdef RDB_GPU_OFFLOAD
         call poison_host(buf)
#endif
         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum)
         !$acc exit data delete(buf)

         call check(error, vmax == ref_max, "device all-negative vmax must equal maxval")
         if (allocated(error)) exit checks
         call check(error, vmax < 0.0_wp, "device all-negative vmax must be negative")
         if (allocated(error)) exit checks
         call check(error, vmin == ref_min, "device all-negative vmin must equal minval")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - ref_sum) <= SUM_RTOL*ref_abs, &
                    "device all-negative vsum must match sum()")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_device_resident_negative

   subroutine test_device_resident_nan(error)
      !! The missing-data skip must happen where the buffer LIVES.  Same
      !! discriminator as the two tests above -- the host copy is poisoned on
      !! a GPU build, so only a reduction that read the mapped device copy
      !! can return the pattern's statistics -- with NaN sentinels folded in,
      !! so it also proves the `ieee_is_finite` guard survives device codegen
      !! (an unguarded device `+` reduction hands back NaN; a relaxed-FP
      !! `min`/`max` lowered to a NaN-blind select hands back the sentinel).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 12, N2 = 6, N3 = 4
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: ref(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs
      integer :: n_valid, n_finite

      call fill_pattern_missing(buf, every=3)
      ! `ref` is a plain host array that is never mapped and never poisoned.
      ref = buf
      call finite_reference(ref, ref_min, ref_max, ref_sum, ref_abs, n_finite)

      checks: block
         call check(error, n_finite > 0 .and. n_finite < N1*N2*N3, &
                    "test setup: the buffer must be partly, not wholly, missing")
         if (allocated(error)) exit checks
         call check(error, POISON < ref_min, &
                    "test setup: POISON must sit below the buffer minimum")
         if (allocated(error)) exit checks

         !$acc enter data copyin(buf)
#ifdef RDB_GPU_OFFLOAD
         call poison_host(buf)
#endif
         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum, n_valid)
         !$acc exit data delete(buf)

         call check(error,.not. ieee_is_nan(vsum), &
                    "the device reduction must not let the NaN sentinel poison the sum")
         if (allocated(error)) exit checks
         call check(error, n_valid == n_finite, &
                    "the device reduction must count the finite cells")
         if (allocated(error)) exit checks
         call check(error, vmin == ref_min, &
                    "device vmin must be the mapped device data's finite minimum")
         if (allocated(error)) exit checks
         call check(error, vmax == ref_max, &
                    "device vmax must be the mapped device data's finite maximum")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - ref_sum) <= SUM_RTOL*ref_abs, &
                    "device vsum must be the mapped device data's finite sum")
         if (allocated(error)) exit checks
         call check(error, vmin /= POISON, &
                    "device reduction returned POISON -- it read the HOST buffer")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_device_resident_nan

   subroutine test_device_resident_stats(error)
      !! `diag_field_stats` is the emit path's entry point and it hands the
      !! buffer on to the device kernel through ONE extra call hop, from an
      !! assumed-shape dummy to an explicit-shape one.  That hop is exactly
      !! where a compiler that cannot prove contiguity would sequence-
      !! associate via a host copy — whose address misses the device present
      !! table, so the kernel would reduce stale host memory with no error.
      !! Same poisoned-host discriminator as above, so a shim that lost
      !! device residency comes back as POISON rather than passing quietly.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 9, N2 = 5, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: ref(N1, N2, N3)
      real(wp) :: vmin, vmax, vmean
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs
      integer :: n_valid, n_total, n_finite

      call fill_pattern_missing(buf, every=4)
      ref = buf
      call finite_reference(ref, ref_min, ref_max, ref_sum, ref_abs, n_finite)

      checks: block
         call check(error, n_finite > 0 .and. n_finite < N1*N2*N3, &
                    "test setup: the buffer must be partly, not wholly, missing")
         if (allocated(error)) exit checks
         call check(error, POISON < ref_min, &
                    "test setup: POISON must sit below the buffer minimum")
         if (allocated(error)) exit checks

         !$acc enter data copyin(buf)
#ifdef RDB_GPU_OFFLOAD
         call poison_host(buf)
#endif
         call diag_field_stats(buf, on_device=.true., vmin=vmin, vmax=vmax, &
                               vmean=vmean, n_valid=n_valid, n_total=n_total)
         !$acc exit data delete(buf)

         call check(error, vmin /= POISON, &
                    "the stats shim returned POISON -- it lost device residency")
         if (allocated(error)) exit checks
         call check(error, n_total == N1*N2*N3, "n_total must be the whole buffer size")
         if (allocated(error)) exit checks
         call check(error, n_valid == n_finite, &
                    "the stats shim must count the device buffer's finite cells")
         if (allocated(error)) exit checks
         call check(error, vmin == ref_min, "shim vmin must be the device data's finite minimum")
         if (allocated(error)) exit checks
         call check(error, vmax == ref_max, "shim vmax must be the device data's finite maximum")
         if (allocated(error)) exit checks
         call check(error, abs(vmean - ref_sum/real(n_finite, wp)) <= SUM_RTOL*ref_abs, &
                    "shim vmean must average the device data's finite cells")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_device_resident_stats

   ! ---------------------------------------------------------------------
   ! Missing-data (NaN sentinel) handling
   !
   ! A diagnostic buffer legitimately carries IEEE NaN: `fill_tracer_impl`
   ! writes one into every land column and every dynamically vanished layer,
   ! because 0 degC / 0 PSU are legal ocean values and must not be confused
   ! with "no water here" (`test_ocean_diag/fill_temperature_vanished_layer_is_nan`
   ! pins that sentinel).  The emit path therefore has to reduce over the
   ! FINITE cells only.
   !
   ! The bug these tests lock out: a plain `minval`/`maxval`/`sum` over the
   ! whole buffer.  Comparisons with NaN are FALSE, so min/max silently
   ! skipped the sentinel cells while `sum` propagated them -- the console
   ! line came out as `min= 1.5E+01  max= 1.5E+01  mean= NaN` for runs whose
   ! state was entirely healthy (island / coastline / vanishing-layer
   ! configurations), and the regression suite's NaN gate failed them.
   ! Reproduced identically under gfortran 15.1 and nvfortran 25.5.
   ! ---------------------------------------------------------------------

   subroutine fill_pattern_missing(buf, every)
      !! Stamp the mixed-sign pattern, then overwrite every `every`-th cell
      !! (linear order) with a quiet NaN -- the missing-data sentinel.
      !! Plain sequential loops, never `do concurrent`: host setup only.
      ! assumed-shape-ok: host-only setup helper, no `do concurrent` here.
      real(wp), intent(out) :: buf(:, :, :)
      integer, intent(in) :: every
      integer :: i, j, k, n
      real(wp) :: qnan
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      n = 0
      do k = 1, size(buf, 3)
         do j = 1, size(buf, 2)
            do i = 1, size(buf, 1)
               n = n + 1
               if (mod(n, every) == 0) then
                  buf(i, j, k) = qnan
               else
                  buf(i, j, k) = pattern(i, j, k)
               end if
            end do
         end do
      end do
   end subroutine fill_pattern_missing

   subroutine finite_reference(buf, ref_min, ref_max, ref_sum, ref_abs, n_finite)
      !! Host reference over the finite cells of `buf`, built with plain
      !! sequential loops and an explicit `ieee_is_nan` test -- never the
      !! intrinsics under test, and never a NaN comparison.
      ! assumed-shape-ok: host-only assertion helper, no `do concurrent` here.
      real(wp), intent(in) :: buf(:, :, :)
      real(wp), intent(out) :: ref_min, ref_max, ref_sum, ref_abs
      integer, intent(out) :: n_finite
      integer :: i, j, k
      ref_min = huge(1.0_wp)
      ref_max = -huge(1.0_wp)
      ref_sum = 0.0_wp
      ref_abs = 0.0_wp
      n_finite = 0
      do k = 1, size(buf, 3)
         do j = 1, size(buf, 2)
            do i = 1, size(buf, 1)
               if (ieee_is_nan(buf(i, j, k))) cycle
               ref_min = min(ref_min, buf(i, j, k))
               ref_max = max(ref_max, buf(i, j, k))
               ref_sum = ref_sum + buf(i, j, k)
               ref_abs = ref_abs + abs(buf(i, j, k))
               n_finite = n_finite + 1
            end do
         end do
      end do
   end subroutine finite_reference

   subroutine test_reduce_skips_nan(error)
      !! The fused reduction must fold only the finite cells and report how
      !! many it folded.  Pre-fix, `vsum` came back NaN here while `vmin` /
      !! `vmax` looked perfectly healthy.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 7, N2 = 5, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vsum
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs
      integer :: n_valid, n_finite

      call fill_pattern_missing(buf, every=4)
      call finite_reference(buf, ref_min, ref_max, ref_sum, ref_abs, n_finite)

      checks: block
         call check(error, n_finite > 0 .and. n_finite < N1*N2*N3, &
                    "test setup: the buffer must be partly, not wholly, missing")
         if (allocated(error)) exit checks

         call diag_reduce_stats(buf, N1, N2, N3, vmin, vmax, vsum, n_valid)

         call check(error,.not. ieee_is_nan(vsum), &
                    "the NaN missing sentinel must not poison the sum")
         if (allocated(error)) exit checks
         call check(error, n_valid == n_finite, "n_valid must count the finite cells")
         if (allocated(error)) exit checks
         call check(error, vmin == ref_min, "vmin must be the minimum over finite cells")
         if (allocated(error)) exit checks
         call check(error, vmax == ref_max, "vmax must be the maximum over finite cells")
         if (allocated(error)) exit checks
         call check(error, abs(vsum - ref_sum) <= SUM_RTOL*ref_abs, &
                    "vsum must be the sum over finite cells")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_reduce_skips_nan

   subroutine test_stats_skips_nan(error)
      !! The emitted `[diag]` statistic: mean over the finite cells, divided
      !! by the VALID count -- not by the array size, which would silently
      !! drag the mean toward zero even once the NaN itself was handled.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 9, N2 = 6, N3 = 4
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vmean
      real(wp) :: ref_min, ref_max, ref_sum, ref_abs, ref_mean, wrong_mean
      integer :: n_valid, n_total, n_finite

      call fill_pattern_missing(buf, every=5)
      call finite_reference(buf, ref_min, ref_max, ref_sum, ref_abs, n_finite)
      ref_mean = ref_sum/real(n_finite, wp)
      wrong_mean = ref_sum/real(N1*N2*N3, wp)

      checks: block
         call check(error, n_finite > 0 .and. n_finite < N1*N2*N3, &
                    "test setup: the buffer must be partly, not wholly, missing")
         if (allocated(error)) exit checks
         call check(error, abs(ref_mean - wrong_mean) > SUM_RTOL*ref_abs, &
                    "test setup: the two divisors must give distinguishable means")
         if (allocated(error)) exit checks

         call diag_field_stats(buf, on_device=.false., vmin=vmin, vmax=vmax, &
                               vmean=vmean, n_valid=n_valid, n_total=n_total)

         call check(error,.not. ieee_is_nan(vmean), &
                    "the emitted mean must not be NaN when the field itself is healthy")
         if (allocated(error)) exit checks
         call check(error, n_total == N1*N2*N3, "n_total must be the whole buffer size")
         if (allocated(error)) exit checks
         call check(error, n_valid == n_finite, "n_valid must count the finite cells")
         if (allocated(error)) exit checks
         call check(error, vmin == ref_min, "vmin must be the minimum over finite cells")
         if (allocated(error)) exit checks
         call check(error, vmax == ref_max, "vmax must be the maximum over finite cells")
         if (allocated(error)) exit checks
         call check(error, abs(vmean - ref_mean) <= SUM_RTOL*ref_abs, &
                    "the mean must divide the finite sum by the VALID cell count")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_skips_nan

   subroutine test_stats_all_finite(error)
      !! No missing cell anywhere -- the overwhelmingly common case, and the
      !! one that must not move: every case in the regression corpus without
      !! land or a vanishing layer has to keep emitting the console numbers
      !! it always did.  The statistic must take the UNMASKED intrinsic
      !! branch and divide by the full size.
      !!
      !! Exactness follows this file's stated policy: `min`/`max` are
      !! order-independent and asserted bit-identical to the intrinsics; the
      !! mean carries a float `+` and so is compared within SUM_RTOL, scaled
      !! by the mean of |buf|.  (Asserting the mean bit-identical to a
      !! separately-compiled `sum(buf)/N` expression fails on nvfortran,
      !! which reassociates the two sums differently -- a property of the
      !! reference expression, not of the code under test, which runs the
      !! same single `sum(buf)` it always ran.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 11, N2 = 7, N3 = 3
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vmean, ref_mean, tol
      integer :: n_valid, n_total

      call fill_pattern(buf, use_int=.false., shift=0.0_wp)
      ref_mean = sum(buf)/real(N1*N2*N3, wp)
      tol = SUM_RTOL*sum(abs(buf))/real(N1*N2*N3, wp)

      checks: block
         call diag_field_stats(buf, on_device=.false., vmin=vmin, vmax=vmax, &
                               vmean=vmean, n_valid=n_valid, n_total=n_total)

         call check(error, n_valid == N1*N2*N3 .and. n_total == N1*N2*N3, &
                    "an all-finite buffer must report every cell valid")
         if (allocated(error)) exit checks
         call check(error, vmin == minval(buf), "vmin must equal minval bit-for-bit")
         if (allocated(error)) exit checks
         call check(error, vmax == maxval(buf), "vmax must equal maxval bit-for-bit")
         if (allocated(error)) exit checks
         call check(error, abs(vmean - ref_mean) <= tol, &
                    "the all-finite mean must divide the whole-buffer sum by the whole size")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_all_finite

   subroutine test_stats_all_missing(error)
      !! A wholly dry / wholly masked field: nothing to average.  Reporting
      !! the reduction's untouched `+huge` / `-huge` seeds would be nonsense
      !! and a 0 would read as a legal ocean value, so all three statistics
      !! come back as the missing-value sentinel.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 4, N2 = 3, N3 = 2
      real(wp) :: buf(N1, N2, N3)
      real(wp) :: vmin, vmax, vmean
      integer :: n_valid, n_total

      call fill_pattern_missing(buf, every=1)

      checks: block
         call diag_field_stats(buf, on_device=.false., vmin=vmin, vmax=vmax, &
                               vmean=vmean, n_valid=n_valid, n_total=n_total)

         call check(error, n_valid == 0, "an all-NaN buffer must report zero valid cells")
         if (allocated(error)) exit checks
         call check(error, n_total == N1*N2*N3, "n_total must be the whole buffer size")
         if (allocated(error)) exit checks
         call check(error, vmean == DIAG_MISSING_VALUE, &
                    "an all-missing field's mean must be the missing-value sentinel")
         if (allocated(error)) exit checks
         call check(error, vmin == DIAG_MISSING_VALUE .and. vmax == DIAG_MISSING_VALUE, &
                    "an all-missing field's min/max must be the missing-value sentinel")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_all_missing

   subroutine poison_host(buf)
      !! Overwrite the HOST copy of an already-mapped buffer.  A plain
      !! sequential loop, never `do concurrent` and never array syntax, so
      !! there is no chance of the store being offloaded to the device (which
      !! would poison the very copy under test).
      ! assumed-shape-ok: host-only helper, no `do concurrent` here.
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: i, j, k
      do k = 1, size(buf, 3)
         do j = 1, size(buf, 2)
            do i = 1, size(buf, 1)
               buf(i, j, k) = POISON
            end do
         end do
      end do
   end subroutine poison_host

   ! ---------------------------------------------------------------------
   ! on_device flag lifecycle
   ! ---------------------------------------------------------------------

   subroutine test_on_device_flag(error)
      !! `ocean_diag_t%on_device` is what routes the emit path between the
      !! device reduction and the host fallback, so its lifecycle is part of
      !! the contract: false before `enter_data`, true after, false again
      !! after `exit_data`.  A manager that never had `init` called must stay
      !! false through both, because both impls return early on `is_init`.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 6, NY = 4
      type(hgrid_t) :: grid
      type(ocean_diag_t) :: diag
      type(ocean_diag_t) :: uninit

      checks: block
         ! An uninitialised manager: enter_data/exit_data are no-ops.
         call uninit%enter_data()
         call check(error,.not. uninit%on_device, &
                    "enter_data on an uninitialised manager must not set on_device")
         if (allocated(error)) exit checks

         call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
         call diag%init(grid)
         call diag%register("probe", units="m", fill=dummy_fill, n1=NX, n2=NY, n3=1)

         call check(error,.not. diag%on_device, &
                    "on_device must be .false. before enter_data")
         if (allocated(error)) exit checks

         call diag%enter_data()
         call check(error, diag%on_device, &
                    "on_device must be .true. after enter_data")
         if (allocated(error)) exit checks

         call diag%exit_data()
         call check(error,.not. diag%on_device, &
                    "on_device must be .false. after exit_data")
         if (allocated(error)) exit checks

         call diag%destroy()
      end block checks
   end subroutine test_on_device_flag

   subroutine dummy_fill(state_handle, buf)
      !! Minimal `diag_fill_proc` so `register` has something to bind; the
      !! lifecycle test never fires a cadence, so this is never called.
      ! assumed-shape-ok: matches the diag_fill_proc abstract interface.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      buf = 0.0_wp
      if (.false.) then
         select type (state_handle)
         class default
         end select
      end if
   end subroutine dummy_fill

end module test_ocean_diag_reduce
