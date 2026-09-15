!! Unit tests for Phase-2 (6a) flux-accumulated tracer-transport
!! scaffold: the `dt_tracer_advect_ratio` configure-time validation
!! predicate and the `continuity_t` flux-accumulator slots
!! (`uhtr`/`vhtr`/`t_dyn_rel_adv`).
!!
!! 6a ships the plumbing + a bit-identical bypass only — the windowed
!! accumulate+drain numerics land in 6b.  These tests therefore cover:
!!   (1) the host-side fail-loud constraint `dt_therm_ratio` must be an
!!       integer multiple of `dt_tracer_advect_ratio` (and both >= 1),
!!       exercised through the real predicate the setup layer calls; and
!!   (2) the accumulator slots allocate at the expected C-grid face
!!       shapes and are zeroed on init.
module test_ocean_dt_tracer_advect
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_dyn, only: ocean_dt_tracer_advect_ratios_ok
   use rdb_continuity, only: continuity_t
   implicit none
   private

   integer, parameter :: NGHOST = 3

   public :: collect_ocean_dt_tracer_advect_tests

contains

   subroutine collect_ocean_dt_tracer_advect_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ratios_reject_non_multiple", test_reject_non_multiple), &
                  new_unittest("ratios_reject_below_one", test_reject_below_one), &
                  new_unittest("ratios_accept_valid_multiple", test_accept_valid_multiple), &
                  new_unittest("ratios_accept_identity", test_accept_identity), &
                  new_unittest("accum_slots_alloc_shape", test_accum_alloc_shape), &
                  new_unittest("accum_slots_zeroed", test_accum_zeroed) &
                  ]
   end subroutine collect_ocean_dt_tracer_advect_tests

   ! ----------------------------------------------------------------
   ! Configure-time ratio validation (the real predicate that
   ! configure_ocean_vmix error-stops on).
   ! ----------------------------------------------------------------

   subroutine test_reject_non_multiple(error)
      !! dt_therm_ratio = 3, dt_tracer_advect_ratio = 2 is INVALID:
      !! 3 is not an integer multiple of 2, so the ALE remap could land
      !! mid-accumulation-window.  Must be rejected.
      type(error_type), allocatable, intent(out) :: error
      call check(error,.not. ocean_dt_tracer_advect_ratios_ok(3, 2), &
                 "dt_therm_ratio=3 / dt_tracer_advect_ratio=2 must be rejected (non-multiple)")
   end subroutine test_reject_non_multiple

   subroutine test_reject_below_one(error)
      !! dt_tracer_advect_ratio < 1 is INVALID regardless of the other
      !! ratio.
      type(error_type), allocatable, intent(out) :: error
      logical :: both_rejected
      both_rejected = (.not. ocean_dt_tracer_advect_ratios_ok(6, 0)) .and. &
                      (.not. ocean_dt_tracer_advect_ratios_ok(6, -1))
      call check(error, both_rejected, &
                 "dt_tracer_advect_ratio < 1 must be rejected")
   end subroutine test_reject_below_one

   subroutine test_accept_valid_multiple(error)
      !! dt_therm_ratio = 6, dt_tracer_advect_ratio = 2 is VALID:
      !! 6 = 3 * 2, so every accumulation window closes before a remap.
      type(error_type), allocatable, intent(out) :: error
      call check(error, ocean_dt_tracer_advect_ratios_ok(6, 2), &
                 "dt_therm_ratio=6 / dt_tracer_advect_ratio=2 must be accepted (6 = 3*2)")
   end subroutine test_accept_valid_multiple

   subroutine test_accept_identity(error)
      !! The default pair (1, 1) — advect every step, thermo every step
      !! — must be accepted (the bit-identical bypass configuration).
      type(error_type), allocatable, intent(out) :: error
      call check(error, ocean_dt_tracer_advect_ratios_ok(1, 1), &
                 "default ratios (1, 1) must be accepted")
   end subroutine test_accept_identity

   ! ----------------------------------------------------------------
   ! Accumulator slot allocation / shape / zeroing.
   ! ----------------------------------------------------------------

   subroutine test_accum_alloc_shape(error)
      !! `continuity_t%init` must allocate the flux accumulators at the
      !! C-grid face shapes: uhtr(nx+1,ny,nz) east-faces,
      !! vhtr(nx,ny+1,nz) north-faces (nx/ny including ghosts =
      !! nx_total/ny_total).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(continuity_t) :: ct
      integer :: nx, ny, nz

      nz = 4
      call grid%init(8, 6, NGHOST, 1000.0_wp, 1000.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
      call ct%init(grid, nz_ml=nz)

      call check(error, allocated(ct%uhtr) .and. allocated(ct%vhtr), &
                 "uhtr/vhtr must be allocated after continuity init")
      if (allocated(error)) then
         call ct%destroy()
         return
      end if
      call check(error, all(shape(ct%uhtr) == [nx + 1, ny, nz]), &
                 "uhtr must be shaped (nx+1, ny, nz)")
      if (.not. allocated(error)) &
         call check(error, all(shape(ct%vhtr) == [nx, ny + 1, nz]), &
                    "vhtr must be shaped (nx, ny+1, nz)")
      call ct%destroy()
   end subroutine test_accum_alloc_shape

   subroutine test_accum_zeroed(error)
      !! The accumulators (and the scalar clock) must be zeroed on init
      !! — the windowed drain assumes an empty window at the start.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(continuity_t) :: ct

      call grid%init(8, 6, NGHOST, 1000.0_wp, 1000.0_wp)
      ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
      call ct%init(grid, nz_ml=3)

      call check(error, maxval(abs(ct%uhtr)) == 0.0_wp, &
                 "uhtr must be zeroed on init")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(ct%vhtr)) == 0.0_wp, &
                    "vhtr must be zeroed on init")
      if (.not. allocated(error)) &
         call check(error, ct%t_dyn_rel_adv == 0.0_wp, &
                    "t_dyn_rel_adv must be zeroed on init")
      call ct%destroy()
   end subroutine test_accum_zeroed

end module test_ocean_dt_tracer_advect
