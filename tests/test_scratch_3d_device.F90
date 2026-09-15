!! `scratch_3d_buffer_t`: what a kernel reads out of a freshly attached
!! buffer must be ZERO, on the device as well as on the host.
!!
!! `init` allocates `source = 0.0_wp`, so every consumer is entitled to
!! treat an untouched buffer as zero — and some deliberately do: a
!! producer may be SKIPPED for a step and leave its buffer to be read as
!! "whatever the previous step produced" (MOM6's predictor reuses
!! `diffu(u[n-1])`; `&ocean_bt_nml split_scheme = "pred_corr"` skips the
!! viscous recompute in stage 1, and on step 1 there is no previous
!! producer).  `enter_data` attaches with `create`, which allocates
!! device memory and copies NOTHING, so before this contract was
!! enforced that first read returned whatever the device allocator
!! happened to hand back.
!!
!! The test poisons the device pool on purpose: a same-shaped buffer is
!! attached, filled with a large value BY A DEVICE KERNEL, and released,
!! so the block most likely to be recycled into the buffer under test
!! holds non-zero.  A host build cannot see the difference (there is one
!! copy of everything and `init` zeroed it) — this is a GPU-path gate,
!! and it must run on the actual `-stdpar=gpu -gpu=...,mem:separate`
!! build to mean anything.
module test_scratch_3d_device
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   implicit none
   private

   public :: collect_scratch_3d_device_tests

   integer, parameter :: N1 = 13
   integer, parameter :: N2 = 11
   integer, parameter :: N3 = 7
   real(wp), parameter :: POISON = 1.0e30_wp

contains

   subroutine collect_scratch_3d_device_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("scratch_3d_attaches_zeroed", test_attaches_zeroed) &
                  ]
   end subroutine collect_scratch_3d_device_tests

   subroutine test_attaches_zeroed(error)
      type(error_type), allocatable, intent(out) :: error
      type(scratch_3d_buffer_t) :: poisoner, buf
      real(wp) :: total

      ! 1. Poison: attach a same-shaped buffer, write POISON into it with
      !    a DEVICE kernel (so the device copy, not the host copy, is
      !    dirtied), then release it back to the allocator.
      call poisoner%init(N1, N2, N3, "scratch_3d_poisoner")
      call scratch_3d_buffer_enter_data_impl(poisoner)
      call fill(poisoner%data, N1, N2, N3, POISON)
      call scratch_3d_buffer_exit_data_impl(poisoner)
      call poisoner%destroy()

      ! 2. The buffer under test: attached and never written.
      call buf%init(N1, N2, N3, "scratch_3d_under_test")
      call scratch_3d_buffer_enter_data_impl(buf)

      ! 3. Read it back where the kernels read it.  A reduction, not an
      !    `update self`: on the host build the directive is inert and
      !    this reduces the one and only copy, which is the same
      !    question asked of the same memory.
      total = sum_abs(buf%data, N1, N2, N3)

      call scratch_3d_buffer_exit_data_impl(buf)
      call buf%destroy()

      call check(error, total == 0.0_wp, &
                 "scratch_3d_buffer_t must attach ZEROED: a kernel read of a "// &
                 "freshly enter_data'd buffer returned non-zero (device memory "// &
                 "from `create` carries no host value)")
   end subroutine test_attaches_zeroed

   subroutine fill(arr, n1, n2, n3, val)
      !! Device-side write, so the poison lands in the device copy.
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(inout) :: arr(n1, n2, n3)
      real(wp), intent(in) :: val
      integer :: i, j, k
      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         arr(i, j, k) = val
      end do
   end subroutine fill

   function sum_abs(arr, n1, n2, n3) result(total)
      !! Device-side reduction (OpenACC: `do concurrent` reductions are
      !! not reliable on this toolchain — see test_dc_reduce).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: arr(n1, n2, n3)
      real(wp) :: total
      integer :: i, j, k
      total = 0.0_wp
      !$acc parallel loop collapse(3) reduction(+:total) present(arr)
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               total = total + abs(arr(i, j, k))
            end do
         end do
      end do
   end function sum_abs

end module test_scratch_3d_device
