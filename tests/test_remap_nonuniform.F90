!! Analytical tests for the ALE remap's non-uniform-grid reconstruction weights.
module test_remap_nonuniform
   !! `&vcoord_nml remap_nonuniform_weights` replaces PLM's minmod
   !! half-difference and PPM's `(7/12, -1/12)` edge estimate — both of them
   !! EQUAL-THICKNESS specialisations of Colella & Woodward (1984) — with the
   !! thickness-weighted (1.6)-(1.8) forms.
   !!
   !! Why that is worth an analytical gate.  Linear exactness is the property
   !! that makes a remap harmless on a stratified column at REST: the tracer
   !! profile is linear in z to leading order, so an exact scheme returns it
   !! unchanged and a scheme that is merely consistent injects a diapycnal
   !! flux every thermo step.  `remap_boundary_extrap` bought that property
   !! back at `k=1` and `k=nz` — but only for a UNIFORM source column, which
   !! is what the existing gate
   !! (`tests/test_remap_boundary_extrap.F90 :: check_linear_exact`) uses and
   !! says so in its own comment.  Every geometric vcoord family but sigma on
   !! a flat bed hands the remap a STRETCHED source column, and there PLM and
   !! PPM carry an O(dh/h) interior error on the same linear profile.
   !! Measured on random columns (`h` uniform on [0.5, 6.5] m, `nz = 12`,
   !! `dq/dz = 0.25`): max|q_new - q_exact| is 1.5E-01 (PLM) and 1.7E-01
   !! (PPM) with the knob off, and 2.8E-14 / 4.3E-14 with it on.
   !!
   !! PPM_H4 and PQM already carry thickness-weighted stencils; the cases
   !! below assert that too, because it is the reason the knob is scoped to
   !! PLM and PPM rather than applied to all four.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, REMAP_PCM, REMAP_PLM, REMAP_PPM, &
                            REMAP_PPM_H4, REMAP_PQM
   use rdb_remap_column, only: remap_column, remap_column_preconditions_ok
   implicit none
   private

   public :: collect_remap_nonuniform_tests

   integer, parameter :: NZ = 12
   real(wp), parameter :: Q0 = 4.25_wp       !! Tracer value at z = 0.
   real(wp), parameter :: DQDZ = -0.037_wp   !! Linear gradient (per metre).
   real(wp), parameter :: TOL_EXACT = 1.0e-13_wp
      !! Absolute tolerance for "exact for a linear profile".  The column
      !! spans ~1.5 units of q over ~40 m, so this is ~1e-13 relative to the
      !! range — two decades above double round-off on these operand
      !! magnitudes, and twelve decades BELOW the knob-off error.

   ! Two deliberately UNEQUAL thickness sets.  `SRC` is monotonically
   ! stretched (thickest at the bed, thinnest at the surface — a z*-like
   ! column over a shelf); `TGT` is stretched the other way and has a
   ! different distribution entirely, so every target cell straddles a
   ! different number of source cells.  Totals match to the last bit by
   ! construction: both are integers summing to 78.
   real(wp), parameter :: SRC(NZ) = [ &
                          11.0_wp, 10.0_wp, 9.0_wp, 8.0_wp, 7.0_wp, 6.0_wp, &
                          6.0_wp, 5.0_wp, 5.0_wp, 4.0_wp, 4.0_wp, 3.0_wp]
   real(wp), parameter :: TGT(NZ) = [ &
                          2.0_wp, 3.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 7.0_wp, &
                          8.0_wp, 8.0_wp, 9.0_wp, 9.0_wp, 9.0_wp, 9.0_wp]

contains

   subroutine collect_remap_nonuniform_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("linear_exact_plm_stretched", test_linear_plm), &
                  new_unittest("linear_exact_ppm_stretched", test_linear_ppm), &
                  new_unittest("already_exact_ppm_h4_and_pqm", test_linear_h4_pqm), &
                  new_unittest("knob_off_is_first_order", test_off_is_first_order), &
                  new_unittest("knob_off_is_bit_identical", test_off_bit_identical), &
                  new_unittest("uniform_source_ppm_reduces", test_uniform_source_ppm), &
                  new_unittest("uniform_source_plm_swaps_limiter", test_uniform_source_plm), &
                  new_unittest("conservation_to_round_off", test_conservation), &
                  new_unittest("step_profile_stays_bounded", test_monotone_step), &
                  new_unittest("vanished_layers_stay_finite", test_vanished_layers), &
                  new_unittest("preconditions_predicate", test_preconditions) &
                  ]
   end subroutine collect_remap_nonuniform_tests

   pure subroutine linear_column(dz, q)
      !! Cell averages of `q(z) = Q0 + DQDZ*z`, z upward from the bed.  For a
      !! linear profile the cell AVERAGE is the value at the cell CENTRE,
      !! which is what makes the expected answer closed-form on any dz.
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

   subroutine check_linear_exact(error, method, label)
      !! Remap a linear-in-z profile from a STRETCHED source column onto a
      !! differently stretched target of the same total, with BOTH closures
      !! on, and require the exact cell averages back.  Both are needed: the
      !! knob under test fixes the interior, `bnd_extrap` the two outermost
      !! cells, and a column is exact only with the pair.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: method
      character(len=*), intent(in) :: label
      real(wp) :: q_old(NZ), q_new(NZ), q_exact(NZ)
      integer :: k

      call linear_column(SRC, q_old)
      call linear_column(TGT, q_exact)
      call remap_column(method, NZ, SRC, TGT, q_old, q_new, &
                        bnd_extrap=.true., nonunif=.true.)
      do k = 1, NZ
         call check(error, abs(q_new(k) - q_exact(k)) < TOL_EXACT, &
                    label//": a linear profile on a STRETCHED source column must "// &
                    "remap exactly under the non-uniform weights")
         if (allocated(error)) return
      end do
   end subroutine check_linear_exact

   subroutine test_linear_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PLM, "PLM")
   end subroutine test_linear_plm

   subroutine test_linear_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call check_linear_exact(error, REMAP_PPM, "PPM")
   end subroutine test_linear_ppm

   subroutine test_linear_h4_pqm(error)
      !! PPM_H4's (White & Adcroft 2008) H4/H3 edge stencils and PQM's
      !! implicit-h4/h3 solves are thickness-weighted already, so they are
      !! exact on a stretched source with the knob OFF.  Asserting that is
      !! what justifies scoping the knob to PLM and PPM — if it ever stops
      !! being true, this case says so rather than the knob quietly growing.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q_old(NZ), q_new(NZ), q_exact(NZ)
      integer :: k, m
      integer :: methods(2)

      methods = [REMAP_PPM_H4, REMAP_PQM]
      call linear_column(SRC, q_old)
      call linear_column(TGT, q_exact)
      do m = 1, size(methods)
         call remap_column(methods(m), NZ, SRC, TGT, q_old, q_new, bnd_extrap=.true.)
         do k = 1, NZ
            call check(error, abs(q_new(k) - q_exact(k)) < TOL_EXACT, &
                       "PPM_H4/PQM are already linear-exact on a stretched source")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_linear_h4_pqm

   subroutine test_off_is_first_order(error)
      !! The complement, and the reason the knob exists: with the shipped
      !! equal-thickness weights the error on a linear profile is NOT
      !! round-off but O(dh/h * dq/dz * h), and it sits in the INTERIOR — so
      !! the boundary closure alone does not buy exactness on a stretched
      !! column.  Without this case the four exactness assertions above would
      !! keep passing for the wrong reason.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q_old(NZ), q_def(NZ), q_exact(NZ)
      real(wp) :: err_int
      integer :: k

      call linear_column(SRC, q_old)
      call linear_column(TGT, q_exact)
      do k = 1, 2
         if (k == 1) then
            call remap_column(REMAP_PLM, NZ, SRC, TGT, q_old, q_def, bnd_extrap=.true.)
         else
            call remap_column(REMAP_PPM, NZ, SRC, TGT, q_old, q_def, bnd_extrap=.true.)
         end if
         err_int = maxval(abs(q_def(3:NZ - 2) - q_exact(3:NZ - 2)))
         call check(error, err_int > 1.0e-5_wp, &
                    "equal-thickness weights: the INTERIOR error on a linear "// &
                    "profile over a stretched source is first-order, not round-off")
         if (allocated(error)) return
      end do
   end subroutine test_off_is_first_order

   subroutine test_off_bit_identical(error)
      !! `nonunif` absent and `.false.` must agree bit-for-bit, for every
      !! shipped method and on a profile with no structure to hide behind.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q_old(NZ), q_a(NZ), q_b(NZ)
      integer :: k, m, b
      integer :: methods(5)
      logical :: be

      methods = [REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM]
      do k = 1, NZ
         q_old(k) = 12.0_wp - 0.7_wp*real(k, wp) + 0.4_wp*cos(1.7_wp*real(k, wp))
      end do

      do b = 1, 2
         be = (b == 2)
         do m = 1, size(methods)
            call remap_column(methods(m), NZ, SRC, TGT, q_old, q_a, bnd_extrap=be)
            call remap_column(methods(m), NZ, SRC, TGT, q_old, q_b, &
                              bnd_extrap=be, nonunif=.false.)
            do k = 1, NZ
               call check(error, q_a(k) == q_b(k), &
                          "absent and .false. nonunif must be bit-identical")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_off_bit_identical

   subroutine test_uniform_source_ppm(error)
      !! On an EQUAL-thickness source column CW84 (1.6) with the UNLIMITED
      !! (1.7) jump reduces algebraically to the shipped `(7/12, -1/12)`
      !! estimate — the whole correction term collapses because its two
      !! thickness weights become equal.  It must therefore agree to
      !! round-off there for ANY profile, monotone or not, which is what
      !! bounds how far a shipped PPM answer can move when the knob is
      !! switched on: only by the non-uniformity it exists to fix.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_u(NZ), dz_t(NZ), q_old(NZ), q_a(NZ), q_b(NZ)
      integer :: k

      dz_u = 6.5_wp
      dz_t = TGT*(sum(dz_u)/sum(TGT))
      do k = 1, NZ
         q_old(k) = 12.0_wp - 0.7_wp*real(k, wp) + 0.4_wp*cos(1.7_wp*real(k, wp))
      end do
      call remap_column(REMAP_PPM, NZ, dz_u, dz_t, q_old, q_a)
      call remap_column(REMAP_PPM, NZ, dz_u, dz_t, q_old, q_b, nonunif=.true.)
      do k = 1, NZ
         call check(error, abs(q_a(k) - q_b(k)) < 1.0e-12_wp, &
                    "on a uniform source the non-uniform PPM edge weights must "// &
                    "reduce to the shipped (7/12, -1/12) estimate")
         if (allocated(error)) return
      end do
   end subroutine test_uniform_source_ppm

   subroutine test_uniform_source_plm(error)
      !! PLM is the deliberate exception, and this case pins it so nobody
      !! reads the knob as "inert on a uniform column".  CW84 (1.7) on equal
      !! thicknesses is the CENTRED difference under the (1.8) bound, where
      !! the shipped kernel uses the strictly more diffusive
      !! `0.5*minmod(dq_l, dq_r)` — so the knob swaps PLM's limiter as well
      !! as its weighting, and the answer moves even with no stretching.
      !! What must hold is that the swap is a REFINEMENT: still monotone
      !! (no new extremum on a step), still conservative, and now
      !! linear-exact.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_u(NZ), dz_t(NZ), q_old(NZ), q_a(NZ), q_b(NZ), q_exact(NZ)
      real(wp) :: diff
      integer :: k

      dz_u = 6.5_wp
      dz_t = TGT*(sum(dz_u)/sum(TGT))
      do k = 1, NZ
         q_old(k) = 12.0_wp - 0.7_wp*real(k, wp) + 0.4_wp*cos(1.7_wp*real(k, wp))
      end do
      call remap_column(REMAP_PLM, NZ, dz_u, dz_t, q_old, q_a)
      call remap_column(REMAP_PLM, NZ, dz_u, dz_t, q_old, q_b, nonunif=.true.)
      diff = maxval(abs(q_a - q_b))
      call check(error, diff > 1.0e-6_wp, &
                 "PLM: the knob swaps minmod for the CW84 limiter, so a uniform "// &
                 "source is NOT bit-identical — pin that rather than assume it")
      if (allocated(error)) return

      ! Still monotone on a step.
      q_old = 0.0_wp
      q_old(7:NZ) = 30.0_wp
      call remap_column(REMAP_PLM, NZ, dz_u, dz_t, q_old, q_b, nonunif=.true.)
      do k = 1, NZ
         call check(error, q_b(k) >= -1.0e-12_wp .and. q_b(k) <= 30.0_wp + 1.0e-12_wp, &
                    "PLM under the CW84 limiter must stay monotone on a uniform source")
         if (allocated(error)) return
      end do

      ! And linear-exact, which minmod already was on a uniform column.
      dz_t = 6.5_wp*0.5_wp
      dz_t(1:NZ/2) = 6.5_wp*1.5_wp
      dz_t = dz_t*(sum(dz_u)/sum(dz_t))
      call linear_column(dz_u, q_old)
      call linear_column(dz_t, q_exact)
      call remap_column(REMAP_PLM, NZ, dz_u, dz_t, q_old, q_b, &
                        bnd_extrap=.true., nonunif=.true.)
      do k = 1, NZ
         call check(error, abs(q_b(k) - q_exact(k)) < TOL_EXACT, &
                    "PLM under the CW84 limiter must stay linear-exact")
         if (allocated(error)) return
      end do
   end subroutine test_uniform_source_plm

   subroutine test_conservation(error)
      !! The weights change the reconstruction, never the integral: the
      !! column total survives to round-off for every method.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q_old(NZ), q_new(NZ), mass_old, mass_new
      integer :: k, m
      integer :: methods(4)

      methods = [REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM]
      do k = 1, NZ
         q_old(k) = 34.0_wp + 0.5_wp*real(k, wp) - 0.9_wp*cos(0.8_wp*real(k, wp))
      end do
      mass_old = sum(q_old*SRC)
      do m = 1, size(methods)
         call remap_column(methods(m), NZ, SRC, TGT, q_old, q_new, &
                           bnd_extrap=.true., nonunif=.true.)
         mass_new = sum(q_new*TGT)
         call check(error, abs(mass_new - mass_old) < 1.0e-11_wp*abs(mass_old), &
                    "the non-uniform weights must not break conservation")
         if (allocated(error)) return
      end do
   end subroutine test_conservation

   subroutine test_monotone_step(error)
      !! Linear exactness must not be bought with monotonicity: CW84's (1.8)
      !! bound is retained, so a step profile gains no new extremum.  Checked
      !! with `bnd_extrap` OFF, so the whole column is under the limiter
      !! (boundary extrapolation deliberately leaves the two outermost cells
      !! free — that is its own test's subject, not this one's).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q_old(NZ), q_new(NZ)
      integer :: k, m
      integer :: methods(2)

      methods = [REMAP_PLM, REMAP_PPM]
      q_old = 0.0_wp
      q_old(7:NZ) = 30.0_wp
      do m = 1, size(methods)
         call remap_column(methods(m), NZ, SRC, TGT, q_old, q_new, nonunif=.true.)
         do k = 1, NZ
            call check(error, q_new(k) >= -1.0e-12_wp .and. q_new(k) <= 30.0_wp + 1.0e-12_wp, &
                       "a step profile must gain no new extremum under the "// &
                       "non-uniform weights")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_monotone_step

   subroutine test_vanished_layers(error)
      !! Every new denominator is a SUM of thicknesses, armoured with
      !! `H_DIV_EPS` — the pure 1/0 role, since the caller's precondition
      !! already guarantees the summands are non-negative.  An interior layer
      !! at exactly zero and another at 1e-12 m (the `zstar_full` /
      !! `z_fixed` filler geometry) must therefore leave the result finite,
      !! conservative, and free of new extrema rather than producing a NaN
      !! or an Inf from a 0/0.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_o(NZ), dz_n(NZ), q_old(NZ), q_new(NZ)
      real(wp) :: mass_old, mass_new, slack
      integer :: k, m
      integer :: methods(2)

      methods = [REMAP_PLM, REMAP_PPM]
      dz_o = 50.0_wp
      dz_o(5) = 0.0_wp
      dz_o(6) = 1.0e-12_wp
      dz_n = sum(dz_o)/real(NZ, wp)
      do k = 1, NZ
         q_old(k) = 10.0_wp + real(k, wp)
      end do
      mass_old = sum(q_old*dz_o)
      slack = max(abs(q_old(2) - q_old(1)), abs(q_old(NZ) - q_old(NZ - 1))) + 1.0e-9_wp
      do m = 1, size(methods)
         call remap_column(methods(m), NZ, dz_o, dz_n, q_old, q_new, &
                           bnd_extrap=.true., nonunif=.true.)
         mass_new = sum(q_new*dz_n)
         call check(error, abs(mass_new - mass_old) < 1.0e-11_wp*abs(mass_old), &
                    "vanished interior layers must stay conservative")
         if (allocated(error)) return
         do k = 1, NZ
            call check(error, ieee_is_finite(q_new(k)), &
                       "vanished interior layers must not produce a non-finite value")
            if (allocated(error)) return
            ! `bnd_extrap` deliberately lets `k=1`/`k=nz` leave the source
            ! range by up to one cell-mean difference — that excursion is
            ! what makes it exact on a linear profile.  The bound here is
            ! that the DEGENERATE column does not widen it further, i.e. no
            ! 1/0 has turned a vanished layer into a runaway edge.
            call check(error, q_new(k) >= minval(q_old) - slack .and. &
                       q_new(k) <= maxval(q_old) + slack, &
                       "vanished interior layers must not produce a runaway value")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_vanished_layers

   subroutine test_preconditions(error)
      !! `remap_column_preconditions_ok` is the `pure` predicate the
      !! cadence-bounded caller asserts (audit V5, V6).  It must accept a
      !! well-posed column, reject a negative source thickness — which makes
      !! the interface stack non-monotone and lets the overlap sweep CREATE
      !! mass — and reject a column-total mismatch, which silently deletes
      !! the non-overlapping tail.  A land column (both totals zero) must
      !! pass rather than trip on a 0/0 relative test.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_o(NZ), dz_n(NZ)
      real(wp), parameter :: RTOL = 1.0e-10_wp

      call check(error, remap_column_preconditions_ok(NZ, SRC, TGT, RTOL), &
                 "a matched, non-negative column pair must pass")
      if (allocated(error)) return

      dz_o = SRC
      dz_o(3) = -2.0_wp
      dz_o(4) = SRC(4) + 2.0_wp
      call check(error,.not. remap_column_preconditions_ok(NZ, dz_o, TGT, RTOL), &
                 "a negative SOURCE thickness must be refused even when the "// &
                 "column total still matches")
      if (allocated(error)) return

      dz_n = TGT
      dz_n(5) = -1.0_wp
      dz_n(6) = TGT(6) + 1.0_wp
      call check(error,.not. remap_column_preconditions_ok(NZ, SRC, dz_n, RTOL), &
                 "a negative TARGET thickness must be refused")
      if (allocated(error)) return

      dz_n = TGT*0.9_wp
      call check(error,.not. remap_column_preconditions_ok(NZ, SRC, dz_n, RTOL), &
                 "a 10% short target column must be refused")
      if (allocated(error)) return

      dz_o = 0.0_wp
      dz_n = 0.0_wp
      call check(error, remap_column_preconditions_ok(NZ, dz_o, dz_n, RTOL), &
                 "a land column (both totals zero) must pass")
   end subroutine test_preconditions

end module test_remap_nonuniform
