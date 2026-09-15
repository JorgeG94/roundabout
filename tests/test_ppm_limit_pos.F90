!! Unit tests for `ppm_limit_pos` (rdb_continuity).
!!
!! `ppm_limit_pos` is MOM6's `PPM_limit_pos` analogue: it adjusts the
!! PPM face values (h_left, h_right) so the parabolic reconstruction's
!! interior minimum sits at exactly `h_min`.  When `h_centre <= h_min`
!! it FORCES the parabola to a constant by setting `h_left = h_right =
!! h_centre` — equivalent to upwind for that cell.  This is the
!! load-bearing mechanism preventing mass-flux overshoot at thin /
!! vanishing layers.
!!
!! Tests:
!!   1. No-op when h_centre is well above h_min and curvature small
!!      (max interior in parabola; the limiter shouldn't touch it).
!!   2. Force constant when h_centre <= h_min (the "vanishing layer"
!!      branch).
!!   3. Curvature scaled down when h_centre is marginally above h_min
!!      and the parabola predicts an interior minimum below h_min.
module test_ppm_limit_pos
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_continuity, only: ppm_limit_pos
   implicit none
   private

   public :: collect_ppm_limit_pos_tests

contains

   subroutine collect_ppm_limit_pos_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("noop_when_above_h_min", test_noop), &
                  new_unittest("force_constant_at_vanishing", test_force_constant), &
                  new_unittest("scale_curvature_marginal", test_scale_curvature), &
                  new_unittest("noop_when_interior_max", test_noop_interior_max) &
                  ]
   end subroutine collect_ppm_limit_pos_tests

   subroutine test_noop(error)
      !! Smooth thick column — h_centre = 1000, h_left = h_right = 1010
      !! (gentle parabolic dip).  Interior min = ?  Let's pick values
      !! where the interior min is well above h_min, so the limiter
      !! must NOT touch the edges.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hL, hR, hL_in, hR_in
      real(wp), parameter :: h_centre = 1000.0_wp, h_min = 1.0_wp
      checks: block
         ! Set h_L = h_R = 1010 (symmetric parabola, interior min at centre = 1000)
         hL_in = 1010.0_wp; hR_in = 1010.0_wp
         hL = hL_in; hR = hR_in
         call ppm_limit_pos(h_centre, hL, hR, h_min)
         ! curv = 3·(2010 - 2000) = 30; dh = 0; |dh| < curv → interior min branch
         ! interior min at h_centre = 1000 > h_min = 1 → no force constant
         ! 12·30·(1000-1) = 359640;  curv² + 3·dh² = 900 + 0 = 900
         ! 359640 > 900 → no scaling
         ! → no-op
         call check(error, abs(hL - hL_in) < 1.0e-14_wp .and. abs(hR - hR_in) < 1.0e-14_wp, &
                    "smooth thick column: limit_pos should be a no-op")
      end block checks
   end subroutine test_noop

   subroutine test_force_constant(error)
      !! Vanishing layer — h_centre = 0.05 m < h_min = 1.0 m.
      !! Limiter MUST force h_L = h_R = h_centre.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hL, hR
      real(wp), parameter :: h_centre = 0.05_wp, h_min = 1.0_wp
      checks: block
         hL = 0.1_wp; hR = 0.08_wp  ! some interior-min-prone edges
         call ppm_limit_pos(h_centre, hL, hR, h_min)
         call check(error, abs(hL - h_centre) < 1.0e-14_wp, &
                    "vanishing layer: h_left should be forced to h_centre")
         if (allocated(error)) exit checks
         call check(error, abs(hR - h_centre) < 1.0e-14_wp, &
                    "vanishing layer: h_right should be forced to h_centre")
      end block checks
   end subroutine test_force_constant

   subroutine test_scale_curvature(error)
      !! Marginal case: h_centre = 2, h_min = 1, but PPM parabola has
      !! strong curvature that would predict an interior min < h_min.
      !! Pick h_L = h_R = 50 — large symmetric edges → strong curv.
      !! Interior min at h_centre = 2 (by symmetry), but the
      !! parabolic profile away from centre dips lower... wait, with
      !! h_L = h_R the interior extremum is AT h_centre.  Hmm.
      !!
      !! Actually rethink: with h_L = h_R > h_centre, the parabola
      !! curves UP — interior is a MIN at h_centre.  curv > 0,
      !! |dh| = 0 < curv, h_centre > h_min, so we check the scaling
      !! condition: 12·curv·(h_centre - h_min) vs curv² + 3·dh².
      !! curv = 3·(100 - 4) = 288; 12·288·(2-1) = 3456;
      !! curv² + 3·dh² = 82944 + 0 = 82944 → 3456 < 82944 → SCALE.
      !! After scaling, edges shrink toward centre.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hL, hR, hL_in, hR_in
      real(wp), parameter :: h_centre = 2.0_wp, h_min = 1.0_wp
      checks: block
         hL_in = 50.0_wp; hR_in = 50.0_wp
         hL = hL_in; hR = hR_in
         call ppm_limit_pos(h_centre, hL, hR, h_min)
         ! Both edges should have shrunk toward h_centre = 2
         call check(error, hL < hL_in .and. hR < hR_in, &
                    "marginal: edges should shrink toward h_centre")
         if (allocated(error)) exit checks
         call check(error, hL > h_centre .and. hR > h_centre, &
                    "marginal: edges should stay above h_centre (curv direction preserved)")
         if (allocated(error)) exit checks
         ! By symmetry h_L should equal h_R
         call check(error, abs(hL - hR) < 1.0e-12_wp, &
                    "marginal: symmetric input should give symmetric output")
      end block checks
   end subroutine test_scale_curvature

   subroutine test_noop_interior_max(error)
      !! Interior maximum (h_centre > h_L, h_centre > h_R) — the
      !! limiter checks `curv > 0` (which is false for interior max),
      !! so it must be a no-op regardless of h_min.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hL, hR, hL_in, hR_in
      real(wp), parameter :: h_centre = 10.0_wp, h_min = 1.0_wp
      checks: block
         hL_in = 8.0_wp; hR_in = 9.0_wp  ! both below centre → curv < 0 (max)
         hL = hL_in; hR = hR_in
         call ppm_limit_pos(h_centre, hL, hR, h_min)
         call check(error, abs(hL - hL_in) < 1.0e-14_wp .and. abs(hR - hR_in) < 1.0e-14_wp, &
                    "interior max: limit_pos should not touch the edges")
      end block checks
   end subroutine test_noop_interior_max

end module test_ppm_limit_pos
