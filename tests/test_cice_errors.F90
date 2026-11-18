module test_cice_errors
   !! Tests of `rdb_ice_errors` on a seeded 64x64 grid mapped to the device:
   !! the OR-mask, every per-bit count, the first-offender cell and word, the
   !! FATAL/COUNTED classification, the clean-grid zero, and the clear kernel.
   !! Expected values come from an independent host brute-force over the seeded
   !! grid and from hand-planted cells, never from `ice_err_reduce` itself.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_ice_errors, only: ICE_NERR, ICE_ERR_FATAL_MASK, ICE_ERR_COUNTED_MASK, &
                             ICE_ERR_ITD_HBNEW_ORDER, ICE_ERR_TIN_RANGE, ICE_ERR_NONFINITE, &
                             ICE_ERR_FLUX_DROPPED, ICE_ERR_REMAP_FALLBACK, ICE_ERR_RIDGE_PASSES, &
                             ICE_ERR_TSFC_NOCONV, ice_err_summary_t, ice_err_name, &
                             ice_err_fatal_set, ice_err_clear, ice_err_reduce, ice_err_is_fatal, &
                             ice_err_dump, ice_err_counted_line, ice_err_accumulate
   implicit none
   private

   public :: collect_cice_errors_tests

   integer, parameter :: NX = 64, NY = 64
   integer, parameter :: NG = 2
      !! Ghost width of the padded grid in the seeded test.
   integer, parameter :: NXA = NX + 2*NG, NYA = NY + 2*NG

contains

   subroutine collect_cice_errors_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("registry_classification", test_registry_classification), &
                  new_unittest("clean_grid_reduces_to_zero", test_clean_grid), &
                  new_unittest("seeded_grid_mask_counts_first", test_seeded_grid), &
                  new_unittest("counted_only_grid_has_no_first_offender", test_counted_only), &
                  new_unittest("fatal_on_error_promotes_itd_bit", test_fatal_on_error), &
                  new_unittest("clear_kernel_zeroes_device_grid", test_clear), &
                  new_unittest("totals_and_console_line", test_totals_line) &
                  ]
   end subroutine collect_cice_errors_tests

   subroutine test_registry_classification(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: b
      integer(int32) :: all
      all = 0_int32
      do b = 0, ICE_NERR - 1
         call check(error, ice_err_name(b) /= "UNKNOWN", "every registered bit is named")
         if (allocated(error)) return
      end do
      call check(error, ice_err_name(ICE_NERR) == "UNKNOWN")
      if (allocated(error)) return
      ! counted by default: bits 10, 22, 23, 24 only
      call check(error, ICE_ERR_COUNTED_MASK == &
                 ior(ior(ishft(1_int32, 10), ishft(1_int32, 22)), ior(ishft(1_int32, 23), ishft(1_int32, 24))))
      if (allocated(error)) return
      call check(error, iand(ICE_ERR_FATAL_MASK, ICE_ERR_COUNTED_MASK) == 0_int32, "disjoint")
      if (allocated(error)) return
      call check(error, ior(ICE_ERR_FATAL_MASK, ICE_ERR_COUNTED_MASK) == ishft(1_int32, ICE_NERR) - 1_int32, &
                 "union covers exactly the registered bits")
      if (allocated(error)) return
      call check(error, ICE_ERR_TSFC_NOCONV == 1_int32 .and. ICE_ERR_FLUX_DROPPED == 16777216_int32)
      if (allocated(error)) return
      call check(error, popcnt(ICE_ERR_FATAL_MASK) == 21, "twenty-one fatal bits by default")
   end subroutine test_registry_classification

   subroutine test_clean_grid(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int32) :: err(NX, NY)
      type(ice_err_summary_t) :: s
      err = 0_int32
      !$acc enter data copyin(err)
      call ice_err_reduce(err, NX, NY, 1, NX, 1, NY, 0, 0, NX, ice_err_fatal_set(.false.), s)
      !$acc exit data delete(err)
      call check(error, s%mask == 0_int32 .and. s%first_idx == 0_int64 .and. s%first_word == 0_int32)
      if (allocated(error)) return
      call check(error, all(s%count == 0_int32) .and. .not. ice_err_is_fatal(s, ICE_ERR_FATAL_MASK))
   end subroutine test_clean_grid

   subroutine seed_full(err)
      !! Owned cells (offset NG) get every one of the 25 bits on a distinct,
      !! bit-dependent lattice.  Fatal bits are kept out of the first four owned rows
      !! and out of owned cells (i<17, j=5) so the lowest fatal owned cell is the
      !! planted (17,5); counted bits start at owned (1,1).  The ghost ring carries
      !! EVERY bit and must never be counted nor named the first offender.
      integer(int32), intent(out) :: err(NXA, NYA)
      integer :: i, j, b
      err = ishft(1_int32, ICE_NERR) - 1_int32
      do j = 1, NY
         do i = 1, NX
            err(i + NG, j + NG) = 0_int32
            do b = 0, ICE_NERR - 1
               if (mod(i*(b + 3) + j*(2*b + 1), b + 2) /= 0) cycle
               if (btest(ICE_ERR_FATAL_MASK, b)) then
                  if (j < 5 .or. (j == 5 .and. i < 17)) cycle
               end if
               err(i + NG, j + NG) = ibset(err(i + NG, j + NG), b)
            end do
         end do
      end do
      err(1 + NG, 1 + NG) = ior(err(1 + NG, 1 + NG), ICE_ERR_ITD_HBNEW_ORDER)
      err(17 + NG, 5 + NG) = ior(err(17 + NG, 5 + NG), ior(ICE_ERR_TIN_RANGE, ICE_ERR_NONFINITE))
   end subroutine seed_full

   subroutine test_seeded_grid(error)
      type(error_type), allocatable, intent(out) :: error
      call seeded_case(error, 0, 0, NX)
      if (allocated(error)) return
      call seeded_case(error, 100, 200, 1000)
   end subroutine test_seeded_grid

   subroutine seeded_case(error, ioff, joff, nxg)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ioff, joff, nxg
      integer(int32) :: err(NXA, NYA), expect_mask, fmask
      integer :: expect_count(0:ICE_NERR - 1), i, j, b
      integer(int64) :: first
      type(ice_err_summary_t) :: s
      call seed_full(err)
      fmask = ice_err_fatal_set(.false.)
      expect_mask = 0_int32
      expect_count = 0
      first = 0_int64
      do j = 1, NY
         do i = 1, NX
            expect_mask = ior(expect_mask, err(i + NG, j + NG))
            do b = 0, ICE_NERR - 1
               if (btest(err(i + NG, j + NG), b)) expect_count(b) = expect_count(b) + 1
            end do
            if (first == 0_int64 .and. iand(err(i + NG, j + NG), fmask) /= 0_int32) &
               first = int(j + joff - 1, int64)*nxg + (i + ioff)
         end do
      end do
      ! the seed must exercise every bit with several owned cells, or the test is vacuous
      call check(error, expect_mask == ishft(1_int32, ICE_NERR) - 1_int32 .and. minval(expect_count) >= 2, &
                 "seed sets all 25 bits")
      if (allocated(error)) return
      call check(error, first == int(4 + joff, int64)*nxg + 17 + ioff, "seed puts the first fatal cell at (17,5)")
      if (allocated(error)) return

      !$acc enter data copyin(err)
      call ice_err_reduce(err, NXA, NYA, 1 + NG, NX + NG, 1 + NG, NY + NG, ioff, joff, nxg, fmask, s)
      !$acc exit data delete(err)

      call check(error, s%mask == expect_mask, "OR-mask")
      if (allocated(error)) return
      do b = 0, ICE_NERR - 1
         call check(error, int(s%count(b)) == expect_count(b), "per-bit count (ghosts excluded) "//ice_err_name(b))
         if (allocated(error)) return
      end do
      call check(error, s%first_idx == first, "first-offender GLOBAL linear index")
      if (allocated(error)) return
      call check(error, s%first_word == err(17 + NG, 5 + NG), "first-offender word")
      if (allocated(error)) return
      call check(error, ice_err_is_fatal(s, fmask))
      if (allocated(error)) return
      call ice_err_dump(s, nxg, fmask)   ! must log, not abort
   end subroutine seeded_case

   subroutine test_counted_only(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int32) :: err(NX, NY)
      type(ice_err_summary_t) :: s
      err = 0_int32
      err(3, 7) = ior(ICE_ERR_REMAP_FALLBACK, ICE_ERR_RIDGE_PASSES)
      err(60, 60) = ICE_ERR_REMAP_FALLBACK
      !$acc enter data copyin(err)
      call ice_err_reduce(err, NX, NY, 1, NX, 1, NY, 0, 0, NX, ice_err_fatal_set(.false.), s)
      !$acc exit data delete(err)
      call check(error, s%mask == ior(ICE_ERR_REMAP_FALLBACK, ICE_ERR_RIDGE_PASSES))
      if (allocated(error)) return
      call check(error, s%count(22) == 2_int32 .and. s%count(23) == 1_int32 .and. s%count(24) == 0_int32)
      if (allocated(error)) return
      call check(error, s%first_idx == 0_int64 .and. .not. ice_err_is_fatal(s, ice_err_fatal_set(.false.)))
   end subroutine test_counted_only

   subroutine test_fatal_on_error(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int32) :: err(NX, NY)
      type(ice_err_summary_t) :: s
      call check(error, iand(ice_err_fatal_set(.false.), ICE_ERR_ITD_HBNEW_ORDER) == 0_int32)
      if (allocated(error)) return
      call check(error, ice_err_fatal_set(.true.) == ior(ICE_ERR_FATAL_MASK, ICE_ERR_ITD_HBNEW_ORDER))
      if (allocated(error)) return
      err = 0_int32
      err(10, 2) = ICE_ERR_ITD_HBNEW_ORDER
      !$acc enter data copyin(err)
      call ice_err_reduce(err, NX, NY, 1, NX, 1, NY, 0, 0, NX, ice_err_fatal_set(.true.), s)
      !$acc exit data delete(err)
      call check(error, s%first_idx == NX + 10 .and. s%first_word == ICE_ERR_ITD_HBNEW_ORDER)
      if (allocated(error)) return
      call check(error, ice_err_is_fatal(s, ice_err_fatal_set(.true.)) .and. &
                 .not. ice_err_is_fatal(s, ice_err_fatal_set(.false.)))
   end subroutine test_fatal_on_error

   subroutine test_clear(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int32) :: err(NXA, NYA)
      call seed_full(err)
      !$acc enter data copyin(err)
      call ice_err_clear(err, NXA, NYA)
      !$acc update self(err)
      !$acc exit data delete(err)
      call check(error, all(err == 0_int32))
   end subroutine test_clear

   subroutine test_totals_line(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_err_summary_t) :: s
      integer(int64) :: tot(0:ICE_NERR - 1)
      tot = 0_int64
      s%count(22) = 3_int32
      s%count(0) = 5_int32
      call ice_err_accumulate(tot, s)
      call ice_err_accumulate(tot, s)
      call check(error, tot(22) == 6_int64 .and. tot(0) == 10_int64)
      if (allocated(error)) return
      call check(error, ice_err_counted_line(tot, ice_err_fatal_set(.false.)) == "REMAP_FALLBACK=6", &
                 "fatal bits never appear on the counted line")
   end subroutine test_totals_line

end module test_cice_errors
