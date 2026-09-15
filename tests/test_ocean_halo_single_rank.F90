!! Unit tests for the ocean halo exchange at ONE rank.
!!
!! There is no stub backend any more -- these drive the real
!! `rdb_ocean_halo`, which at px==1 / py==1 takes its local path instead
!! of entering a collective.  That path is what a single-rank build runs
!! in production, so it is worth pinning directly.
!!
!! Two test cases:
!!   1. periodic_local_wrap: calls ocean_halo_centre with periodic=.true.
!!      and verifies that ghost cells are filled from the opposite physical edge.
!!   2. nonperiodic_noop: calls ocean_halo_centre with periodic=.false.;
!!      ghost cells are unchanged (no-op).
module test_ocean_halo_single_rank
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy, &
                             ocean_halo_is_init, &
                             ocean_halo_centre, ocean_halo_face_x, &
                             ocean_halo_face_y
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_halo_single_rank_tests

   integer, parameter :: NX_PHYS = 8
      !! Physical cells in x (small grid for fast tests)
   integer, parameter :: NY_PHYS = 6
      !! Physical cells in y
   integer, parameter :: NGHOST = 3
      !! Ghost cell width (required >= 3 for periodic)
   logical :: comm_inited = .false.
      !! Guard: only initialise comm-env once across both test subroutines

contains

   subroutine collect_ocean_halo_single_rank_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("halo_single_rank_periodic_local_wrap", test_periodic_local_wrap), &
                  new_unittest("halo_single_rank_nonperiodic_noop", test_nonperiodic_noop) &
                  ]
   end subroutine collect_ocean_halo_single_rank_tests

   ! -----------------------------------------------------------------
   ! Test 1: periodic_local_wrap
   !
   ! Fill centre-2D field with unique i+j values in physical cells.
   ! Ghost cells are initialised to -1.  Call ocean_halo_centre
   ! with periodic_x=.true., periodic_y=.true.  Verify that:
   !   - west ghost  fld(k, j) == fld(k + nx_phys, j)  for k=1..ng
   !   - east ghost  fld(ng+nx_phys+k, j) == fld(ng+k, j)  k=1..ng
   !   - south ghost fld(i, k) == fld(i, k + ny_phys)  for k=1..ng
   !   - north ghost fld(i, ng+ny_phys+k) == fld(i, ng+k)  k=1..ng
   ! -----------------------------------------------------------------

   subroutine test_periodic_local_wrap(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = NX_PHYS + 2*NGHOST
      integer, parameter :: ny = NY_PHYS + 2*NGHOST
      real(wp) :: fld(nx, ny)
      type(decomp_t) :: decomp
      integer :: i, j, k, ng

      ng = NGHOST

      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if

      ! Single rank: px=py=1
      call decomp_init(decomp, NX_PHYS, NY_PHYS, 1, 1, 0)
      call ocean_halo_init(decomp, ng, .true., .true.)

      checks: block
         ! Fill physical cells
         fld = -1.0_wp
         do j = ng + 1, ng + NY_PHYS
            do i = ng + 1, ng + NX_PHYS
               fld(i, j) = real(i + j, wp)
            end do
         end do

         call ocean_halo_centre(fld)

         ! West ghosts: fld(k, j) == fld(k + nx_phys, j)  for k=1..ng
         do j = ng + 1, ng + NY_PHYS
            do k = 1, ng
               call check(error, fld(k, j), fld(k + NX_PHYS, j), &
                          "west ghost", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! East ghosts: fld(ng+nx_phys+k, j) == fld(ng+k, j)  for k=1..ng
         do j = ng + 1, ng + NY_PHYS
            do k = 1, ng
               call check(error, fld(ng + NX_PHYS + k, j), fld(ng + k, j), &
                          "east ghost", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! South ghosts: fld(i, k) == fld(i, k + ny_phys)  for k=1..ng
         do k = 1, ng
            do i = ng + 1, ng + NX_PHYS
               call check(error, fld(i, k), fld(i, k + NY_PHYS), &
                          "south ghost", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! North ghosts: fld(i, ng+ny_phys+k) == fld(i, ng+k)  for k=1..ng
         do k = 1, ng
            do i = ng + 1, ng + NX_PHYS
               call check(error, fld(i, ng + NY_PHYS + k), fld(i, ng + k), &
                          "north ghost", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do
      end block checks

      call ocean_halo_destroy()
   end subroutine test_periodic_local_wrap

   ! -----------------------------------------------------------------
   ! Test 2: nonperiodic_noop
   !
   ! Ghost cells initialised to sentinel value -1.  Call
   ! ocean_halo_centre with periodic_x=.false., periodic_y=.false.
   ! Verify that no ghost cell was modified (undecomposed + non-periodic
   ! is a true no-op).
   ! -----------------------------------------------------------------

   subroutine test_nonperiodic_noop(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = NX_PHYS + 2*NGHOST
      integer, parameter :: ny = NY_PHYS + 2*NGHOST
      real(wp) :: fld(nx, ny)
      type(decomp_t) :: decomp
      integer :: i, j, k, ng
      real(wp), parameter :: SENTINEL = -999.0_wp

      ng = NGHOST

      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if

      call decomp_init(decomp, NX_PHYS, NY_PHYS, 1, 1, 0)
      call ocean_halo_init(decomp, ng, .false., .false.)

      checks: block
         ! Fill everything (physical + ghost) with a known sentinel
         fld = SENTINEL
         ! Overwrite physical cells with different values
         do j = ng + 1, ng + NY_PHYS
            do i = ng + 1, ng + NX_PHYS
               fld(i, j) = real(i*j, wp)
            end do
         end do

         call ocean_halo_centre(fld)

         ! West ghosts must still be SENTINEL
         do j = ng + 1, ng + NY_PHYS
            do k = 1, ng
               call check(error, fld(k, j), SENTINEL, "west ghost unchanged", &
                          thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! East ghosts must still be SENTINEL
         do j = ng + 1, ng + NY_PHYS
            do k = 1, ng
               call check(error, fld(ng + NX_PHYS + k, j), SENTINEL, &
                          "east ghost unchanged", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! South ghosts must still be SENTINEL
         do k = 1, ng
            do i = ng + 1, ng + NX_PHYS
               call check(error, fld(i, k), SENTINEL, "south ghost unchanged", &
                          thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do

         ! North ghosts must still be SENTINEL
         do k = 1, ng
            do i = ng + 1, ng + NX_PHYS
               call check(error, fld(i, ng + NY_PHYS + k), SENTINEL, &
                          "north ghost unchanged", thr=1.0e-12_wp)
               if (allocated(error)) exit checks
            end do
         end do
      end block checks

      call ocean_halo_destroy()
   end subroutine test_nonperiodic_noop

end module test_ocean_halo_single_rank
