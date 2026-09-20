!! Analytical tests for the ALE remap's boundary-cell closure.
module test_remap_boundary_extrap
   !! `&vcoord_nml remap_boundary_extrap` (MOM6 `BOUNDARY_EXTRAPOLATION`)
   !! replaces the PCM flatten at `k=1` / `k=nz` with the linear-exact
   !! one-sided edge pair, so the whole column reproduces a profile that is
   !! LINEAR IN z exactly — for any source and any target thickness set.
   !!
   !! Why that is worth an analytical gate.  A stratified ocean at REST has
   !! a tracer profile that is linear in z to leading order, and under a
   !! terrain-following coordinate the ALE remap runs on it every thermo
   !! step.  With the default closure the remap is FIRST-ORDER in the two
   !! cells next to the bed and the surface (identically so for PLM, PPM,
   !! PPM_H4 and PQM — they share the closure), so every step injects a
   !! spurious diapycnal tracer flux there.  Over a slope the injection
   !! differs between neighbouring columns, which is a horizontal density
   !! gradient, which is a pressure-gradient force; with rotation it feeds
   !! a growing grid mode trapped in exactly those layers.  Measured on a
   !! 48x6x15 sigma-over-slope rest case (bed 226 -> 709 m, f-plane 75S,
   !! zero viscosity, exact FV pressure gradient): `sigma_En` = 0.355 /day
   !! with the default closure against 0.047 /day — the no-remap floor —
   !! with extrapolation on, and day-45 En 1.7E-20 against 4.7E-25.
   !!
   !! Four of the cases below FAIL with `bnd_extrap` absent or `.false.`,
   !! which is the point: they measure the closure, not the scheme.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, REMAP_PCM, REMAP_PLM, REMAP_PPM, &
                            REMAP_PPM_H4, REMAP_PQM
   use rdb_remap_column, only: remap_column
   implicit none
   private

   public :: collect_remap_boundary_extrap_tests

   integer, parameter :: NZ = 8
   real(wp), parameter :: Q0 = 3.5_wp       !! Tracer value at z = 0.
   real(wp), parameter :: DQDZ = -0.025_wp  !! Linear gradient (per metre).
   real(wp), parameter :: TOL_EXACT = 1.0e-12_wp
      !! Relative-to-range tolerance for "exact for a linear profile".  The
      !! column spans ~1 unit of q over ~40 m, so this is ~1e-12 absolute —
      !! several decades above double round-off and several decades BELOW
      !! the O(h) default-closure error the test is separating from.

contains

   subroutine collect_remap_boundary_extrap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("linear_exact_plm", test_linear_exact_plm), &
                  new_unittest("linear_exact_ppm", test_linear_exact_ppm), &
                  new_unittest("linear_exact_ppm_h4", test_linear_exact_ppm_h4), &
                  new_unittest("linear_exact_pqm", test_linear_exact_pqm), &
                  new_unittest("default_closure_is_first_order", test_default_is_first_order), &
                  new_unittest("default_off_is_bit_identical", test_default_off_bit_identical), &
                  new_unittest("pcm_ignores_the_knob", test_pcm_unaffected), &
                  new_unittest("conservation_with_extrap", test_conservation), &
                  new_unittest("uniform_profile_untouched", test_uniform_profile), &
                  new_unittest("bounded_by_the_local_gradient", test_bounded) &
                  ]
   end subroutine collect_remap_boundary_extrap_tests

   pure subroutine linear_column(dz, q)
      !! Cell averages of `q(z) = Q0 + DQDZ*z` on the thicknesses `dz`,
      !! with z measured upward from the bed (k=1).  For a linear profile
      !! the cell AVERAGE equals the value at the cell CENTRE, which is
      !! what makes the expected answer closed-form.
      real(wp), intent(in) :: dz(NZ)
      real(wp), intent(out) :: q(NZ)
      integer :: k
      real(wp) :: z
      z = 0.0_wp
      do k = 1, NZ
         q(k) = Q0 + DQDZ*(z + 0.5_wp*dz(k))
         z = z + dz(k)
      end do
   end subroutine linear_column

   subroutine check_linear_exact(error, method, label, uniform_src)
      !! Remap a linear-in-z profile onto a DIFFERENT, non-uniform target
      !! column of the same total, with the boundary closure on, and
      !! require the exact cell averages back.
      !!
      !! `uniform_src` selects the source column, and the distinction is
      !! real rather than cosmetic: PLM's minmod half-difference and PPM's
      !! `(7/12,-1/12)` edge estimate are both written for UNIFORM layers,
      !! so their INTERIORS are linear-exact only there.  PPM_H4 and PQM
      !! carry the thickness-weighted stencils and are linear-exact on any
      !! thickness set.  The boundary closure under test is exact for any
      !! thickness pair in all four — so a uniform source isolates it for
      !! PLM/PPM, and the non-uniform source additionally exercises the
      !! interior for the two schemes that claim it.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: method
      character(len=*), intent(in) :: label
      logical, intent(in) :: uniform_src
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_new(NZ), q_exact(NZ)
      integer :: k

      if (uniform_src) then
         dz_old = 4.5_wp
      else
         dz_old = [3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp]
      end if
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      call linear_column(dz_old, q_old)
      call linear_column(dz_new, q_exact)

      call remap_column(method, NZ, dz_old, dz_new, q_old, q_new, bnd_extrap=.true.)

      do k = 1, NZ
         call check(error, abs(q_new(k) - q_exact(k)) < TOL_EXACT, &
                    label//": linear profile must remap exactly with boundary "// &
                    "extrapolation (cell k)")
         if (allocated(error)) return
      end do
   end subroutine check_linear_exact

   subroutine test_linear_exact_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PLM, "PLM", .true.)
   end subroutine test_linear_exact_plm

   subroutine test_linear_exact_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PPM, "PPM", .true.)
   end subroutine test_linear_exact_ppm

   subroutine test_linear_exact_ppm_h4(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PPM_H4, "PPM_H4", .true.)
      if (allocated(error)) return
      call check_linear_exact(error, REMAP_PPM_H4, "PPM_H4 (non-uniform source)", .false.)
   end subroutine test_linear_exact_ppm_h4

   subroutine test_linear_exact_pqm(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PQM, "PQM", .true.)
      if (allocated(error)) return
      call check_linear_exact(error, REMAP_PQM, "PQM (non-uniform source)", .false.)
   end subroutine test_linear_exact_pqm

   subroutine test_default_is_first_order(error)
      !! The complement of the four cases above, and the reason the knob
      !! exists: with the default closure the error on a linear profile is
      !! NOT round-off, it is O(h*dq/dz), and it sits in the boundary cells
      !! while the interior stays exact.  Asserting the defect explicitly
      !! keeps the four exactness tests honest — without this, a remap that
      !! silently became exact everywhere would leave them passing for the
      !! wrong reason.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_def(NZ), q_exact(NZ)
      real(wp) :: err_bnd, err_int
      integer :: k

      dz_old = 4.5_wp
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      call linear_column(dz_old, q_old)
      call linear_column(dz_new, q_exact)

      call remap_column(REMAP_PPM, NZ, dz_old, dz_new, q_old, q_def)

      err_bnd = max(abs(q_def(1) - q_exact(1)), abs(q_def(NZ) - q_exact(NZ)))
      err_int = 0.0_wp
      do k = 3, NZ - 2
         err_int = max(err_int, abs(q_def(k) - q_exact(k)))
      end do

      call check(error, err_bnd > 1.0e-4_wp, &
                 "default closure: the boundary-cell error on a linear profile "// &
                 "is first-order, not round-off")
      if (allocated(error)) return
      call check(error, err_int < TOL_EXACT, &
                 "default closure: the INTERIOR is already exact for a linear "// &
                 "profile — the defect is the boundary closure alone")
   end subroutine test_default_is_first_order

   subroutine test_default_off_bit_identical(error)
      !! `bnd_extrap` absent, `.false.`, and every shipped scheme: the
      !! result must be bit-for-bit what the pre-knob kernel produced.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_a(NZ), q_b(NZ)
      integer :: k, m
      integer :: methods(5)

      methods = [REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM]
      dz_old = [3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp]
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      q_old = [12.0_wp, 11.0_wp, 9.5_wp, 9.5_wp, 7.0_wp, 4.0_wp, 3.0_wp, 2.5_wp]

      do m = 1, size(methods)
         call remap_column(methods(m), NZ, dz_old, dz_new, q_old, q_a)
         call remap_column(methods(m), NZ, dz_old, dz_new, q_old, q_b, bnd_extrap=.false.)
         do k = 1, NZ
            call check(error, q_a(k) == q_b(k), &
                       "absent and .false. bnd_extrap must be bit-identical")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_default_off_bit_identical

   subroutine test_pcm_unaffected(error)
      !! PCM has no reconstruction to close, so the knob is inert there —
      !! the documented behaviour, and the reason a PCM remap is NOT a
      !! workaround for the rest-state mode.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_a(NZ), q_b(NZ)
      integer :: k

      dz_old = [3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp]
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      call linear_column(dz_old, q_old)

      call remap_column(REMAP_PCM, NZ, dz_old, dz_new, q_old, q_a, bnd_extrap=.false.)
      call remap_column(REMAP_PCM, NZ, dz_old, dz_new, q_old, q_b, bnd_extrap=.true.)
      do k = 1, NZ
         call check(error, q_a(k) == q_b(k), "PCM must ignore bnd_extrap")
         if (allocated(error)) return
      end do
   end subroutine test_pcm_unaffected

   subroutine test_conservation(error)
      !! Extrapolation changes the reconstruction, never the integral: the
      !! column total must still be conserved to round-off for every scheme.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_new(NZ)
      real(wp) :: mass_old, mass_new
      integer :: m
      integer :: methods(4)

      methods = [REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM]
      dz_old = [3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp]
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      q_old = [12.0_wp, 11.0_wp, 9.5_wp, 9.5_wp, 7.0_wp, 4.0_wp, 3.0_wp, 2.5_wp]
      mass_old = sum(q_old*dz_old)

      do m = 1, size(methods)
         call remap_column(methods(m), NZ, dz_old, dz_new, q_old, q_new, bnd_extrap=.true.)
         mass_new = sum(q_new*dz_new)
         call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                    "boundary extrapolation must not break conservation")
         if (allocated(error)) return
      end do
   end subroutine test_conservation

   subroutine test_uniform_profile(error)
      !! A constant column has `dq_up = 0`, so the half-jump is zero and
      !! the extrapolated closure degenerates to the PCM one: a uniform
      !! tracer stays uniform exactly.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_new(NZ)
      integer :: k

      dz_old = [3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp]
      dz_new = [5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 3.0_wp, 3.0_wp]
      q_old = 34.75_wp

      call remap_column(REMAP_PPM, NZ, dz_old, dz_new, q_old, q_new, bnd_extrap=.true.)
      do k = 1, NZ
         call check(error, abs(q_new(k) - 34.75_wp) < 1.0e-13_wp, &
                    "uniform profile must be preserved under extrapolation")
         if (allocated(error)) return
      end do
   end subroutine test_uniform_profile

   subroutine test_bounded(error)
      !! Extrapolation is not free monotonicity: by construction the
      !! boundary cell MAY leave the source range (that is what makes it
      !! exact on a linear profile).  What it may not do is run away — the
      !! `|d| <= |dq_up|` clamp bounds the excursion by ONE cell-mean
      !! difference.  Checked on a step profile, the adversarial case.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(NZ), dz_new(NZ), q_old(NZ), q_new(NZ)
      real(wp) :: q_lo, q_hi, slack
      integer :: k

      dz_old = 4.0_wp
      dz_new = [2.0_wp, 6.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 2.0_wp]
      q_old = [0.0_wp, 0.0_wp, 0.0_wp, 30.0_wp, 30.0_wp, 30.0_wp, 30.0_wp, 30.0_wp]

      q_lo = minval(q_old)
      q_hi = maxval(q_old)
      ! One cell-mean difference of slack at each end — the clamp's bound.
      slack = max(abs(q_old(2) - q_old(1)), abs(q_old(NZ) - q_old(NZ - 1)))

      call remap_column(REMAP_PPM, NZ, dz_old, dz_new, q_old, q_new, bnd_extrap=.true.)
      do k = 1, NZ
         call check(error, q_new(k) >= q_lo - slack - 1.0e-10_wp .and. &
                    q_new(k) <= q_hi + slack + 1.0e-10_wp, &
                    "extrapolated boundary cell must stay within one cell-mean "// &
                    "difference of the source range")
         if (allocated(error)) return
      end do
   end subroutine test_bounded

end module test_remap_boundary_extrap
