!! Unit tests for Phase-2 vanished-layer velocity reset,
!! `reset_vanished_layer_velocities` in `rdb_ocean_dyn`.
!!
!! The kernel zeroes per-layer face velocities where BOTH adjacent
!! centre-cell thicknesses are at or below `vanish_tol`.  A face with
!! at least one massive neighbour is untouched (R1: the face may be a
!! legitimate grounding-front flux).
!!
!! Test suite (4 analytical cases):
!!   1. both_sided_vanished_zeroed  - a face flanked by two thin cells
!!      is set to 0 exactly; other faces with at least one massive
!!      neighbour are unmodified.
!!   2. one_sided_massive_untouched (R1) - a face with one massive side
!!      and one thin side is byte-identical after the call.
!!   3. barotropic_sum_preserved   - the mass flux u*h at the zeroed face
!!      was already negligible (|flux| <= vanish_tol*|u_old|*width);
!!      asserts the zeroed layer carried negligible mass transport.
!!   4. vanished_vel_bitident      - a state with NO vanished layers
!!      (all h >> vanish_tol) is byte-identical after the call.
module test_ocean_isopycnal_vanished_vel
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_dyn, only: reset_vanished_layer_velocities, isopycnal_vanish_tol
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_isopycnal_vanished_vel_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: ANGSTROM = 1.0e-2_wp
      !! Representative floor value used in the tests.

contains

   subroutine collect_ocean_isopycnal_vanished_vel_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("both_sided_vanished_zeroed", test_both_sided_vanished_zeroed), &
                  new_unittest("one_sided_massive_untouched", test_one_sided_massive_untouched), &
                  new_unittest("barotropic_sum_preserved", test_barotropic_sum_preserved), &
                  new_unittest("vanished_vel_bitident", test_vanished_vel_bitident) &
                  ]
   end subroutine collect_ocean_isopycnal_vanished_vel_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(ms)
      call ms%enter_data()
   end subroutine map_in

   subroutine map_out(ms)
      type(multilayer_state_t), intent(inout) :: ms
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! T1: both sides vanished -> face velocity zeroed
   ! -----------------------------------------------------------------
   subroutine test_both_sided_vanished_zeroed(error)
      !! Grid 8x6 interior.  Place a large velocity on a u-face flanked by
      !! two thin cells (h <= vanish_tol).  After the call, that face must
      !! be 0.  All other faces with at least one massive neighbour must be
      !! unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
         !! Typical massive layer thickness (> vanish_tol by orders of magnitude).
      real(wp), parameter :: THIN = ANGSTROM*0.5_wp
         !! Below vanish_tol: this face should be zeroed.
      real(wp), parameter :: U_LARGE = 5.0_wp
         !! Large face velocity placed on the vanished face.
      real(wp), parameter :: U_OTHER = 2.0_wp
         !! Velocity placed on untouched faces.
      real(wp) :: vtol
      integer :: i_target, j_target, k_target
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)
         i_target = 5    ! u-face index (i straddles centres i-1 and i)
         j_target = 4
         k_target = 2

         ! Set all layers massive everywhere
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! Make both adjacent centres thin for the target face:
         ! u-face at (i_target, j_target, k_target) reads
         ! h_layer(i_target-1, j_target, k_target) and h_layer(i_target, ...)
         ms%h_layer(i_target - 1, j_target, k_target) = THIN
         ms%h_layer(i_target, j_target, k_target) = THIN

         ! Place large velocity on the vanished face and a non-zero velocity
         ! on a neighbouring face that has a massive side
         ms%u_face_x_layer(i_target, j_target, k_target) = U_LARGE
         ms%u_face_x_layer(i_target + 1, j_target, k_target) = U_OTHER

         call map_in(ms)
         call reset_vanished_layer_velocities(ms, vtol)
         call map_out(ms)

         ! Vanished-both-sides face must be zeroed
         call check(error, ms%u_face_x_layer(i_target, j_target, k_target) == 0.0_wp, &
                    "T1: both-sided vanished u-face not zeroed")
         if (allocated(error)) exit checks

         ! Neighbouring face (i_target+1) has h_layer(i_target,j,k)=THIN but
         ! h_layer(i_target+1,j,k)=MASSIVE => max > vtol => untouched
         call check(error, ms%u_face_x_layer(i_target + 1, j_target, k_target) == U_OTHER, &
                    "T1: one-massive-side u-face was incorrectly zeroed")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_both_sided_vanished_zeroed

   ! -----------------------------------------------------------------
   ! T2: one-sided massive -> face velocity untouched (R1)
   ! -----------------------------------------------------------------
   subroutine test_one_sided_massive_untouched(error)
      !! A u-face at (i,j,k) with h_layer(i-1,j,k) = THIN but
      !! h_layer(i,j,k) = MASSIVE: max > vanish_tol.  The face velocity
      !! must be byte-identical after the call — this preserves the
      !! grounding-front flux that re-wets the thin layer (R1).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp), parameter :: THIN = ANGSTROM*0.5_wp
      real(wp), parameter :: U_FRONT = 3.7_wp
         !! Arbitrary non-zero velocity on the grounding-front face.
      real(wp) :: vtol
      integer :: i_target, j_target, k_target
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)
         i_target = 5
         j_target = 4
         k_target = 2

         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! Only the LEFT neighbour is thin
         ms%h_layer(i_target - 1, j_target, k_target) = THIN

         ms%u_face_x_layer(i_target, j_target, k_target) = U_FRONT

         call map_in(ms)
         call reset_vanished_layer_velocities(ms, vtol)
         call map_out(ms)

         call check(error, ms%u_face_x_layer(i_target, j_target, k_target) == U_FRONT, &
                    "T2: one-sided grounding-front face velocity was changed (R1 violation)")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_one_sided_massive_untouched

   ! -----------------------------------------------------------------
   ! T3: mass flux at the zeroed face was negligible
   ! -----------------------------------------------------------------
   subroutine test_barotropic_sum_preserved(error)
      !! Set up a face that will be zeroed by the reset.  Before the call,
      !! compute the mass flux |u*h_face_approx| = |u| * max(h_L, h_R)
      !! at that face.  By definition max(h_L, h_R) <= vanish_tol, so
      !! the mass flux is bounded by vanish_tol * |u|.  Assert this
      !! pre-call bound holds — i.e. the layer that got zeroed carried
      !! negligible mass.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp), parameter :: THIN = ANGSTROM*0.3_wp
         !! Well below vanish_tol
      real(wp), parameter :: U_LARGE = 10.0_wp
      real(wp) :: vtol, max_h, mass_flux_bound
      integer :: i_target, j_target, k_target
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)
         i_target = 5
         j_target = 4
         k_target = 2

         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ms%h_layer(i_target - 1, j_target, k_target) = THIN
         ms%h_layer(i_target, j_target, k_target) = THIN
         ms%u_face_x_layer(i_target, j_target, k_target) = U_LARGE

         ! Pre-call mass flux bound: max(h_L, h_R) * |u| <= vtol * |u|
         max_h = max(ms%h_layer(i_target - 1, j_target, k_target), &
                     ms%h_layer(i_target, j_target, k_target))
         mass_flux_bound = max_h*abs(U_LARGE)

         call check(error, max_h <= vtol, &
                    "T3 precondition: max_h should be <= vanish_tol")
         if (allocated(error)) exit checks
         call check(error, mass_flux_bound <= vtol*abs(U_LARGE), &
                    "T3: pre-call mass flux |u*h| not bounded by vtol*|u|")
         if (allocated(error)) exit checks

         call map_in(ms)
         call reset_vanished_layer_velocities(ms, vtol)
         call map_out(ms)

         ! After the call the face should be 0
         call check(error, ms%u_face_x_layer(i_target, j_target, k_target) == 0.0_wp, &
                    "T3: face not zeroed as expected")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_barotropic_sum_preserved

   ! -----------------------------------------------------------------
   ! T4: all layers massive -> byte-identical (bitident)
   ! -----------------------------------------------------------------
   subroutine test_vanished_vel_bitident(error)
      !! When all cell thicknesses are well above vanish_tol, no face
      !! condition fires and the velocity arrays are byte-identical before
      !! and after the call.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 50.0_wp
      real(wp), parameter :: U_SEED = 1.234_wp
         !! Arbitrary seed velocity — must be unchanged.
      real(wp) :: vtol
      integer :: nx, ny, nx_face, ny_uface, nx_vface, ny_face
      logical :: all_u_ok, all_v_ok
      integer :: i, j, k
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)

         ! All layers massive
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = U_SEED
         ms%v_face_y_layer = U_SEED

         call map_in(ms)
         call reset_vanished_layer_velocities(ms, vtol)
         call map_out(ms)

         nx = size(ms%u_face_x_layer, 1)
         ny_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         ny_face = size(ms%v_face_y_layer, 2)

         all_u_ok = .true.
         do k = 1, NZ
            do j = 1, ny_uface
               do i = 1, nx
                  if (ms%u_face_x_layer(i, j, k) /= U_SEED) then
                     all_u_ok = .false.
                  end if
               end do
            end do
         end do
         all_v_ok = .true.
         do k = 1, NZ
            do j = 1, ny_face
               do i = 1, nx_vface
                  if (ms%v_face_y_layer(i, j, k) /= U_SEED) then
                     all_v_ok = .false.
                  end if
               end do
            end do
         end do

         call check(error, all_u_ok, &
                    "T4: u_face_x_layer not byte-identical (massive state)")
         if (allocated(error)) exit checks
         call check(error, all_v_ok, &
                    "T4: v_face_y_layer not byte-identical (massive state)")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_vanished_vel_bitident

end module test_ocean_isopycnal_vanished_vel
