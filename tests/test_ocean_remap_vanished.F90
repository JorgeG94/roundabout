!! The vanished-layer content rule around the ALE tracer remap.
!!
!! Invariant **I1′**: `h <= H_VANISHED ⇒ hTr = h·c_live` for every registered
!! tracer — `c_live` the concentration of the filler's DONOR, the nearest live
!! layer above it (the topmost live layer for fillers above it; `hTr = 0` in a
!! column with no live layer) — with the per-column content sum conserved to
!! round-off.  (It replaced I1, `hTr = 0`, which conserved the content but not
!! a uniform concentration; see `tests/test_ocean_vanished_constancy.F90`.)
!! `i1p_holds` below is an INDEPENDENT oracle for the donor map: a search up
!! then down, not a copy of the included rule.
!!
!! **Who keeps I1′ has changed.**  The tracer remap now computes on
!! CONCENTRATION (White & Adcroft 2008): its input is the donor-consistent
!! `c_old`, its output `c_new·h_new`, with no read-side pool merge, no
!! write-side pool merge and no defect fold.  I1′ is therefore an INPUT
!! convention of the remap (a filler's `c_old` is its donor's) and is
!! re-established, for the output, only by the tail sweep
!! `multilayer_state_t%enforce_vanished_content` — the tests of the sweep, the
!! scan and the host twin are unchanged.  What the remap itself guarantees,
!! and what is asserted here:
!!   - on an I1′ input, the column content is conserved to round-off, every
!!     target concentration (filler targets included) lies inside the range
!!     of the source column, and an exactly empty target holds the mean of
!!     the source cell containing its bottom interface;
!!   - content a filler holds ABOVE `h·c_donor` (parked there by this step's
!!     writers before the sweep) is not read: the result is independent of
!!     it, and the column content falls by exactly that excess, which is
!!     what the remap budget slot records;
!!   - a column whose target is a few ulp short of the source loses exactly
!!     `c·δh` of water, nothing else;
!!   - a UNIFORM column stays uniform through the remap, fillers included;
!!   - a column with no sub-threshold layer equals the plain
!!     `c = hTr/h` → `remap_column` → `c·h` pipeline exactly, except where
!!     the bounding acted, and there by exactly the recorded amount;
!!   - the degenerate column with no live layer at all (land / grounded);
!!   - end-to-end through `ocean_apply_ale_remap_centres` with a
!!     `VCOORD_ZSTAR_FULL` filler column: the `*_budget_remap` contributor
!!     telescopes to zero per column on an I1′ input.
!! And, unchanged: the I1′ enforcement point, the tripwire scan, the merge
!! helper and the cavity donor map.
module test_ocean_remap_vanished
   use rdb_constants, only: wp, REMAP_PPM, REMAP_PLM, H_VANISHED, NZ_STACK_MAX, &
                            VCOORD_ZSTAR_FULL
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_remap_column, only: remap_column
   use rdb_ocean_remap, only: ocean_remap_tracer_column, &
                              ocean_remap_merge_vanished_content, &
                              ocean_apply_ale_remap_centres
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_remap_vanished_tests

contains

   subroutine collect_ocean_remap_vanished_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("section_g_parked_excess_not_read_and_measured", test_section_g_reproducer), &
                  new_unittest("fillers_bounded_conserved_after_remap", test_i1_after_remap), &
                  new_unittest("random_columns_with_fillers", test_random_columns), &
                  new_unittest("no_filler_equals_plain_pipeline_but_recorded_clamp", test_bit_identity), &
                  new_unittest("dead_column_stays_dead", test_dead_column), &
                  new_unittest("merge_helper_properties", test_merge_helper), &
                  new_unittest("remap_budget_telescopes_with_fillers", test_budget_telescopes), &
                  new_unittest("state_enforcement_point_establishes_i1", test_enforcement_point), &
                  new_unittest("state_scan_detects_and_clears", test_i1_scan), &
                  new_unittest("state_scan_zero_crossing_roundoff", test_i1_scan_zero_crossing), &
                  new_unittest("uniform_column_stays_uniform_through_remap", test_remap_constancy), &
                  new_unittest("cavity_top_fillers_take_k_top", test_top_fillers_k_top), &
                  new_unittest("host_enforcement_matches_device_sweep", test_host_twin), &
                  new_unittest("mismatched_column_loses_exactly_c_dh", test_mismatch_fold) &
                  ]
   end subroutine collect_ocean_remap_vanished_tests

   ! ---------------------------------------------------------------- helpers

   function i1p_holds(n, h, q) result(ok)
      !! Independent I1′ oracle for one column: every layer at or below
      !! `H_VANISHED` must hold `h·c_d`, `c_d = q/h` of its donor — the first
      !! live layer found searching UP from it, else the first found searching
      !! DOWN — to `1e-12` relative; a column with no live layer must hold 0.
      integer, intent(in) :: n
      real(wp), intent(in) :: h(:), q(:)
      logical :: ok
      integer :: k, kd, kk
      real(wp) :: c_d
      ok = .true.
      do k = 1, n
         if (h(k) > H_VANISHED) cycle
         kd = 0
         do kk = k + 1, n
            if (h(kk) > H_VANISHED) then
               kd = kk
               exit
            end if
         end do
         if (kd == 0) then
            do kk = k - 1, 1, -1
               if (h(kk) > H_VANISHED) then
                  kd = kk
                  exit
               end if
            end do
         end if
         if (kd == 0) then
            if (q(k) /= 0.0_wp) ok = .false.
         else
            c_d = q(kd)/h(kd)
            if (abs(q(k) - h(k)*c_d) > 1.0e-12_wp*abs(h(k)*c_d)) ok = .false.
         end if
      end do
   end function i1p_holds

   subroutine impose_i1p(n, h, q)
      !! Put a column in the I1′ state: every layer at or below `H_VANISHED`
      !! holds `h·c_d`, `c_d` its donor's concentration (the first live layer
      !! searching UP from it, else DOWN), `0` with no live layer.  Independent
      !! of the included rule (same search as `i1p_holds`).
      integer, intent(in) :: n
      real(wp), intent(in) :: h(:)
      real(wp), intent(inout) :: q(:)
      integer :: k, kd, kk
      real(wp) :: qn(n)
      qn = q(1:n)
      do k = 1, n
         if (h(k) > H_VANISHED) cycle
         kd = 0
         do kk = k + 1, n
            if (h(kk) > H_VANISHED) then
               kd = kk
               exit
            end if
         end do
         if (kd == 0) then
            do kk = k - 1, 1, -1
               if (h(kk) > H_VANISHED) then
                  kd = kk
                  exit
               end if
            end do
         end if
         if (kd == 0) then
            qn(k) = 0.0_wp
         else
            qn(k) = h(k)*(q(kd)/h(kd))
         end if
      end do
      q(1:n) = qn
   end subroutine impose_i1p

   pure function col_sum(n, a) result(s)
      integer, intent(in) :: n
      real(wp), intent(in) :: a(:)
      real(wp) :: s
      integer :: k
      s = 0.0_wp
      do k = 1, n
         s = s + a(k)
      end do
   end function col_sum

   pure subroutine lcg_next(seed, r)
      !! Reproducible pseudo-random in [0,1) — a 31-bit Lehmer generator, so the
      !! "random columns" case is the same on every toolchain.
      integer, intent(inout) :: seed
      real(wp), intent(out) :: r
      integer, parameter :: A = 16807, M = 2147483647, Q = 127773, R_C = 2836
      integer :: hi, lo, t
      hi = seed/Q
      lo = seed - hi*Q
      t = A*lo - R_C*hi
      if (t > 0) then
         seed = t
      else
         seed = t + M
      end if
      r = real(seed, wp)/real(M, wp)
   end subroutine lcg_next

   ! ------------------------------------------------------------------ tests

   subroutine test_section_g_reproducer(error)
      !! §G.4 of the diagnosis: `nz = 4`, bed filler at `h = 1e-4` holding
      !! 34.5 PSU of parked content (3.45e-3 PSU.m) plus the step's advective
      !! increment (8.889e-4), remapped onto a target where the bed layer has
      !! just become live (1.758e-3 m).
      !!
      !! CHANGED with the concentration-form remap.  The read-side pool merge
      !! used to carry the parked excess `hTr − h·c_donor` into the live
      !! layers, so the remap conserved the column.  The remap now reads
      !! `c_old` (the filler's is its donor's, 34.5) and never the filler's
      !! content: the result is BIT-IDENTICAL to the same remap fed the
      !! I1′-consistent column, and the column content falls by exactly the
      !! excess (the advective increment, 8.889e-4), which is the amount the
      !! remap budget slot (`hTr_new − hTr_old`) records.  Carrying that
      !! increment is the tail sweep's job at the end of the step, and
      !! disappears with the content storage.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), hTr_clean(NZ)
      real(wp) :: before, after, excess

      h_old = [1.0e-4_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      h_new = [1.758272e-3_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      h_new(2) = h_new(2) + (sum(h_old) - sum(h_new))
      hTr = [3.45e-3_wp + 8.888926e-4_wp, 690.0_wp, 690.0_wp, 0.0_wp]
      hTr_clean = hTr
      call impose_i1p(NZ, h_old, hTr_clean)
      excess = (hTr(1) - hTr_clean(1)) + (hTr(4) - hTr_clean(4))

      before = col_sum(NZ, hTr)
      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      after = col_sum(NZ, hTr)
      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr_clean, REMAP_PPM)

      checks: block
         call check(error, all(hTr == hTr_clean), &
                    "section-G column: the remap must not read the filler's parked content")
         if (allocated(error)) exit checks
         call check(error, abs((before - after) - excess) <= 1.0e-12_wp*abs(before), &
                    "section-G column: the content falls by exactly the parked excess")
      end block checks
   end subroutine test_section_g_reproducer

   subroutine test_i1_after_remap(error)
      !! Fillers at the bed, in the interior and at the top, and a target whose
      !! filler set differs from the source's, fed an I1′ column.
      !!
      !! CHANGED: the old assertion was I1′ on the OUTPUT, which the
      !! write-side pool merge established.  The remap no longer merges; the
      !! output is `c_new·h_new` with `c_new` a bounded concentration, so what
      !! is asserted is (a) conservation to round-off, (b) every target —
      !! filler targets included — carries a concentration inside the source
      !! range (not `0`, not a made-up value), (c) the bounding amount is
      !! round-off.  I1′ itself is the tail sweep's, tested below.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), c_new(NZ), clamp(NZ)
      real(wp) :: before, after, c_lo, c_hi
      integer :: k
      logical :: bounded

      !            bed filler  live   filler  live   live   filler  live   top filler
      h_old = [1.0e-4_wp, 12.0_wp, 1.0e-4_wp, 30.0_wp, 25.0_wp, 5.0e-5_wp, 18.0_wp, 1.0e-4_wp]
      h_new = [3.0_wp, 9.0_wp, 1.0e-4_wp, 28.0_wp, 27.0_wp, 1.0e-4_wp, 18.9_wp, 1.0e-4_wp]
      h_new(7) = sum(h_old) - (sum(h_new) - h_new(7))
      do k = 1, NZ
         hTr(k) = (34.0_wp + 0.1_wp*real(k, wp))*h_old(k)
      end do
      call impose_i1p(NZ, h_old, hTr)
      c_lo = huge(1.0_wp)
      c_hi = -huge(1.0_wp)
      do k = 1, NZ
         c_lo = min(c_lo, hTr(k)/h_old(k))
         c_hi = max(c_hi, hTr(k)/h_old(k))
      end do
      before = col_sum(NZ, hTr)

      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM, c_new_out=c_new, clamp_out=clamp)
      after = col_sum(NZ, hTr)

      bounded = .true.
      do k = 1, NZ
         if (c_new(k) < c_lo .or. c_new(k) > c_hi) bounded = .false.
      end do

      checks: block
         call check(error, abs(after - before) <= 1.0e-12_wp*abs(before), &
                    "mixed-filler column: content conserved")
         if (allocated(error)) exit checks
         call check(error, bounded, &
                    "every target concentration, filler targets included, lies in the source range")
         if (allocated(error)) exit checks
         call check(error, abs(col_sum(NZ, clamp)) <= 1.0e-12_wp*abs(before), &
                    "the recorded bounding amount is round-off")
      end block checks
   end subroutine test_i1_after_remap

   subroutine test_random_columns(error)
      !! Pseudo-random I1′ columns, filler runs placed at the bed, at the top and
      !! in the interior in turn.  Every column: conservation to 1e-12 relative
      !! and every target concentration (filler targets included) inside the
      !! source range.  CHANGED from "I1′ on the target grid": the remap no
      !! longer pools; the source fillers enter in the I1′ state instead of
      !! carrying parked content.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 12, NTRIAL = 60
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), c_new(NZ)
      real(wp) :: before, after, worst_rel, r, live_new, defect, c_lo, c_hi
      integer :: seed, trial, k, place, k0, k1
      logical :: i1_ok

      seed = 20260921
      worst_rel = 0.0_wp
      i1_ok = .true.
      do trial = 1, NTRIAL
         do k = 1, NZ
            call lcg_next(seed, r)
            h_old(k) = 1.0_wp + 40.0_wp*r
            call lcg_next(seed, r)
            h_new(k) = 1.0_wp + 40.0_wp*r
            call lcg_next(seed, r)
            hTr(k) = (30.0_wp + 8.0_wp*r)*h_old(k)
         end do
         place = mod(trial, 3)          ! 0 bed, 1 top, 2 interior
         select case (place)
         case (0)
            k0 = 1; k1 = 1 + mod(trial, 4)
         case (1)
            k0 = NZ - mod(trial, 4); k1 = NZ
         case default
            k0 = 4; k1 = 4 + mod(trial, 3)
         end select
         do k = k0, k1
            ! Source filler (its content is set to the I1′ value below).
            h_old(k) = 1.0e-4_wp
            call lcg_next(seed, r)
            hTr(k) = 34.5_wp*1.0e-4_wp*(1.0_wp + r)
            ! Target: half the trials keep it a filler, half wake it up.
            if (mod(trial, 2) == 0) then
               h_new(k) = 1.0e-4_wp
            else
               call lcg_next(seed, r)
               h_new(k) = 1.0e-3_wp + r
            end if
         end do
         ! The target grid must span the SAME column as the source — that is
         ! what `compute_target_h` guarantees in production, and `remap_column`
         ! is conservative only on a matched column.  Absorb the mismatch into
         ! the live target layers, proportionally.
         live_new = 0.0_wp
         do k = 1, NZ
            if (h_new(k) > H_VANISHED) live_new = live_new + h_new(k)
         end do
         defect = col_sum(NZ, h_old) - col_sum(NZ, h_new)
         do k = 1, NZ
            if (h_new(k) > H_VANISHED) h_new(k) = h_new(k)*(1.0_wp + defect/live_new)
         end do

         call impose_i1p(NZ, h_old, hTr)
         c_lo = huge(1.0_wp)
         c_hi = -huge(1.0_wp)
         do k = 1, NZ
            c_lo = min(c_lo, hTr(k)/h_old(k))
            c_hi = max(c_hi, hTr(k)/h_old(k))
         end do
         before = col_sum(NZ, hTr)
         call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM, c_new_out=c_new)
         after = col_sum(NZ, hTr)
         worst_rel = max(worst_rel, abs(after - before)/abs(before))
         if (any(c_new < c_lo) .or. any(c_new > c_hi)) i1_ok = .false.
      end do

      checks: block
         call check(error, worst_rel <= 1.0e-12_wp, &
                    "random filler columns: content conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, i1_ok, "random filler columns: every target concentration inside the source range")
      end block checks
   end subroutine test_random_columns

   subroutine test_bit_identity(error)
      !! A column with NO sub-threshold layer equals the plain pipeline
      !! `c = hTr/h` -> `remap_column` -> `c*h`, spelled out here, EXACTLY
      !! (`==`) wherever the bounding did not act; where it did (an overshoot
      !! of a few ulp on a target that coincides with a source extremum) the
      !! difference is the recorded `clamp_out` amount to round-off.
      !! CHANGED from bit-identity to the one-sided remap: the bounding is
      !! new, so exact equality is only owed where it is inactive.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10, NTRIAL = 40
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), hTr_ref(NZ), clamp(NZ)
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX)
      real(wp) :: r
      integer :: seed, trial, k, n_clamped
      logical :: identical, clamp_recorded

      seed = 991
      identical = .true.
      clamp_recorded = .true.
      n_clamped = 0
      do trial = 1, NTRIAL
         do k = 1, NZ
            call lcg_next(seed, r)
            h_old(k) = 0.5_wp + 30.0_wp*r
            call lcg_next(seed, r)
            h_new(k) = 0.5_wp + 30.0_wp*r
            call lcg_next(seed, r)
            hTr(k) = (25.0_wp + 12.0_wp*r)*h_old(k)
         end do
         hTr_ref = hTr

         ! --- the plain pipeline, verbatim
         do k = 1, NZ
            c_old(k) = hTr_ref(k)/h_old(k)
         end do
         call remap_column(REMAP_PPM, NZ, h_old, h_new, c_old(1:NZ), c_new(1:NZ))
         do k = 1, NZ
            hTr_ref(k) = c_new(k)*h_new(k)
         end do

         call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM, clamp_out=clamp)
         do k = 1, NZ
            if (clamp(k) == 0.0_wp) then
               if (hTr(k) /= hTr_ref(k)) identical = .false.
            else
               n_clamped = n_clamped + 1
               if (abs((hTr(k) - hTr_ref(k)) - clamp(k)) > 4.0_wp*epsilon(1.0_wp)*abs(hTr(k))) &
                  clamp_recorded = .false.
            end if
         end do
      end do

      checks: block
         call check(error, identical, &
                    "no sub-threshold layer ⇒ equal to the plain pipeline wherever the bounding is inactive")
         if (allocated(error)) exit checks
         call check(error, clamp_recorded, &
                    "where the bounding acted, hTr moved by exactly the recorded amount")
      end block checks
   end subroutine test_bit_identity

   subroutine test_dead_column(error)
      !! Land / fully grounded column: every layer a filler, content zero by
      !! construction.  Nothing must be created, and I1 is trivially met.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ)
      integer :: k
      logical :: ok

      h_old = 1.0e-4_wp
      h_new = 1.0e-4_wp
      hTr = 0.0_wp
      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      ok = .true.
      do k = 1, NZ
         if (hTr(k) /= 0.0_wp) ok = .false.
      end do
      call check(error, ok, "dead column: remap must not manufacture content")
   end subroutine test_dead_column

   subroutine test_merge_helper(error)
      !! `rdb_vl_merge_content` on its own (through the remap's test shim):
      !! conserves the column sum, leaves every sub-threshold layer holding
      !! its donor's concentration, pools each filler run with the RIGHT donor,
      !! is idempotent to round-off, and is a textual no-op on a column with
      !! no filler.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp) :: h(NZ_STACK_MAX), q(NZ_STACK_MAX), q2(NZ_STACK_MAX)
      real(wp) :: before, after, c_a, c_b
      integer :: k
      logical :: ok

      h = 0.0_wp
      q = 0.0_wp
      h(1:NZ) = [1.0e-4_wp, 1.0e-4_wp, 10.0_wp, 1.0e-4_wp, 20.0_wp, 1.0e-4_wp]
      q(1:NZ) = [1.0_wp, 2.0_wp, 300.0_wp, 4.0_wp, 600.0_wp, 6.0_wp]
      before = col_sum(NZ, q)

      call ocean_remap_merge_vanished_content(NZ, h, q)
      after = col_sum(NZ, q)

      checks: block
         call check(error, abs(after - before) <= 1.0e-13_wp*abs(before), &
                    "merge helper conserves the column sum")
         if (allocated(error)) exit checks
         call check(error, i1p_holds(NZ, h(1:NZ), q(1:NZ)), &
                    "merge helper: every vanished layer holds h*c_live")
         if (allocated(error)) exit checks

         ! Pools: {1, 2, 3} (bed fillers + the live layer ABOVE them) and
         ! {4, 5, 6} (interior filler 4 goes UP into 5; top filler 6 has no
         ! live layer above, so it pools with the topmost live layer, 5).
         c_a = 303.0_wp/(10.0_wp + 2.0e-4_wp)
         c_b = 610.0_wp/(20.0_wp + 2.0e-4_wp)
         call check(error, abs(q(3) - 10.0_wp*c_a) <= 1.0e-12_wp*303.0_wp, &
                    "bed fillers pool with the live layer above")
         if (allocated(error)) exit checks
         call check(error, abs(q(5) - 20.0_wp*c_b) <= 1.0e-12_wp*610.0_wp, &
                    "interior + top fillers pool with the live layer between them")
         if (allocated(error)) exit checks
         call check(error, abs(q(1)/1.0e-4_wp - c_a) <= 1.0e-12_wp*c_a .and. &
                    abs(q(6)/1.0e-4_wp - c_b) <= 1.0e-12_wp*c_b, &
                    "a filler reads its pool's concentration")
         if (allocated(error)) exit checks

         q2 = q
         call ocean_remap_merge_vanished_content(NZ, h, q2)
         ok = .true.
         do k = 1, NZ
            if (abs(q2(k) - q(k)) > 1.0e-13_wp*abs(q(k))) ok = .false.
         end do
         call check(error, ok, "merge helper is idempotent (to round-off)")
         if (allocated(error)) exit checks

         h(1:NZ) = [5.0_wp, 5.0_wp, 10.0_wp, 5.0_wp, 20.0_wp, 5.0_wp]
         q(1:NZ) = [1.0_wp, 2.0_wp, 300.0_wp, 4.0_wp, 600.0_wp, 6.0_wp]
         q2 = q
         call ocean_remap_merge_vanished_content(NZ, h, q2)
         ok = .true.
         do k = 1, NZ
            if (q2(k) /= q(k)) ok = .false.
         end do
         call check(error, ok, "no filler ⇒ merge helper is a textual no-op")
      end block checks
   end subroutine test_merge_helper

   subroutine test_budget_telescopes(error)
      !! End-to-end through the orchestrator on a ZSTAR_FULL column whose bed
      !! layers vanish, the fillers in the I1′ state.  The `*_budget_remap`
      !! contributor must telescope to zero per column — nothing leaves the
      !! column, so the accumulator stays a leak detector — the column content
      !! must be conserved, and (uniform S) every layer, filler included, must
      !! read S0.  CHANGED: the fillers used to enter with PARKED content and
      !! the test asserted I1′ on the output; the remap no longer pools, so
      !! parked content is not carried (see `test_section_g_reproducer`) and
      !! the output assertion is the uniform concentration instead.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 6
      real(wp), parameter :: H_BED = 120.0_wp, S0 = 34.5_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :), h_bed_2d(:, :)
      real(wp), allocatable :: col_before(:, :)
      real(wp) :: s, b, worst_content, worst_budget
      integer :: nx_tot, ny_tot, i, j, k
      logical :: i1_ok

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! Two bed fillers in the I1′ state, four live layers above.
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               if (k <= 2) then
                  ms%h_layer(i, j, k) = 1.0e-4_wp
                  ms%tracers(ms%idx_salinity)%c(i, j, k) = S0*1.0e-4_wp
               else
                  ms%h_layer(i, j, k) = (H_BED - 2.0e-4_wp)/4.0_wp
                  ms%tracers(ms%idx_salinity)%c(i, j, k) = S0
               end if
            end do
         end do
      end do

      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = 20.0_wp
      vc%zstar_n_surf = 2
      allocate (h_bed_2d(nx_tot, ny_tot), source=H_BED)
      call vc%build_zref_full(h_bed_2d)
      allocate (bt_eta(nx_tot, ny_tot), source=0.0_wp)
      allocate (bt_H_ref(nx_tot, ny_tot), source=H_BED)
      allocate (col_before(nx_tot, ny_tot), source=0.0_wp)

      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               col_before(i, j) = col_before(i, j) + ms%tracers(ms%idx_salinity)%c(i, j, k)
            end do
         end do
      end do

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref, method=REMAP_PLM)

      worst_content = 0.0_wp
      worst_budget = 0.0_wp
      i1_ok = .true.
      do j = 1, ny_tot
         do i = 1, nx_tot
            s = 0.0_wp
            b = 0.0_wp
            do k = 1, NZ
               s = s + ms%tracers(ms%idx_salinity)%c(i, j, k)
               b = b + ms%salt_budget_remap(i, j, k)
            end do
            do k = 1, NZ
               if (abs(ms%tracers(ms%idx_salinity)%c(i, j, k) - S0*ms%h_layer(i, j, k)) > &
                   1.0e-12_wp*S0*ms%h_layer(i, j, k)) i1_ok = .false.
            end do
            worst_content = max(worst_content, abs(s - col_before(i, j)))
            worst_budget = max(worst_budget, abs(b))
         end do
      end do

      checks: block
         call check(error, worst_content <= 1.0e-10_wp*abs(col_before(1, 1)), &
                    "ZSTAR_FULL filler column: content conserved through the orchestrator")
         if (allocated(error)) exit checks
         call check(error, worst_budget <= 1.0e-10_wp*abs(col_before(1, 1)), &
                    "remap budget contributor telescopes to zero per column")
         if (allocated(error)) exit checks
         call check(error, i1_ok, "uniform S stays uniform in every layer, fillers included")
      end block checks

      deallocate (bt_eta, bt_H_ref, h_bed_2d, col_before)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_budget_telescopes

   subroutine test_enforcement_point(error)
      !! `multilayer_state_t%enforce_vanished_content` — THE I1′ enforcement
      !! point.  Park inconsistent content in bed and top fillers of every
      !! registered tracer (including a passive one, to prove the registry
      !! loop carries it), sweep, and require: I1′ everywhere, the column
      !! content conserved, and the live column with no filler left
      !! BIT-identical (the sweep is a textual no-op there).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NX = 4, NY = 3, NZ = 6
      real(wp), parameter :: S0 = 34.5_wp
      real(wp), allocatable :: before_s(:, :), live_ref(:, :, :)
      integer :: nx_tot, ny_tot, i, j, k, t, idx_pass
      real(wp) :: s, worst_content
      logical :: i1_ok, live_identical

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ms%register_passive_tracer(grid, "dye", "1", "a passive dye", idx_pass)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! j = 1 columns are all live (bit-identity control); every other column
      ! carries a bed filler at k = 1 and a top filler at k = NZ, both holding
      ! parked content far off their donor's concentration.
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               ms%h_layer(i, j, k) = 20.0_wp
               if (j > 1 .and. (k == 1 .or. k == NZ)) ms%h_layer(i, j, k) = 1.0e-4_wp
               do t = 1, size(ms%tracers)
                  ms%tracers(t)%c(i, j, k) = S0*20.0_wp + real(t, wp)
               end do
            end do
         end do
      end do

      allocate (before_s(nx_tot, ny_tot), source=0.0_wp)
      allocate (live_ref(nx_tot, ny_tot, NZ), source=0.0_wp)
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               before_s(i, j) = before_s(i, j) + ms%tracers(ms%idx_salinity)%c(i, j, k)
               live_ref(i, j, k) = ms%tracers(ms%idx_salinity)%c(i, j, k)
            end do
         end do
      end do

      call ms%enforce_vanished_content(nx_tot, ny_tot)

      i1_ok = .true.
      live_identical = .true.
      worst_content = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            s = 0.0_wp
            do k = 1, NZ
               s = s + ms%tracers(ms%idx_salinity)%c(i, j, k)
               if (j == 1) then
                  if (ms%tracers(ms%idx_salinity)%c(i, j, k) /= live_ref(i, j, k)) &
                     live_identical = .false.
               end if
            end do
            do t = 1, size(ms%tracers)
               if (.not. i1p_holds(NZ, ms%h_layer(i, j, :), ms%tracers(t)%c(i, j, :))) &
                  i1_ok = .false.
            end do
            worst_content = max(worst_content, abs(s - before_s(i, j)))
         end do
      end do

      checks: block
         call check(error, idx_pass > 0, "the passive tracer must register")
         if (allocated(error)) exit checks
         call check(error, i1_ok, "enforcement point: I1' over the whole registry")
         if (allocated(error)) exit checks
         call check(error, worst_content <= 1.0e-12_wp*abs(before_s(1, 1)), &
                    "enforcement point: column content conserved")
         if (allocated(error)) exit checks
         call check(error, live_identical, &
                    "enforcement point: an all-live column is BIT-identical")
      end block checks

      deallocate (before_s, live_ref)
      call ms%destroy()
   end subroutine test_enforcement_point

   subroutine test_i1_scan(error)
      !! `multilayer_state_t%scan_vanished_content` — the tripwire's pure
      !! half.  It must SEE a violation (count + worst `|hTr − h·c_live|`) —
      !! in particular the EMPTY filler the previous rule left behind — pass a
      !! filler that holds its donor's concentration, and report clean once
      !! the enforcement point has run.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      integer :: nx_tot, ny_tot, i, j, k, n_bad
      real(wp) :: worst

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               ms%h_layer(i, j, k) = 10.0_wp
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 345.0_wp
               ms%tracers(ms%idx_temperature)%c(i, j, k) = 10.0_wp
            end do
            ! One bed filler, EMPTY — the I1 state, which I1′ forbids.
            ms%h_layer(i, j, 1) = 1.0e-4_wp
            ms%tracers(ms%idx_salinity)%c(i, j, 1) = 0.0_wp
            ms%tracers(ms%idx_temperature)%c(i, j, 1) = 0.0_wp
         end do
      end do

      checks: block
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == 2*nx_tot*ny_tot, &
                    "the scan must count every violating (cell, tracer)")
         if (allocated(error)) exit checks
         ! Salinity is the worse: |0 - 1e-4 * 34.5|.
         call check(error, abs(worst - 3.45e-3_wp) <= 1.0e-15_wp, &
                    "the scan must report the worst |hTr - h*c_live| found in a filler")
         if (allocated(error)) exit checks

         ! A filler holding exactly its donor's concentration is clean.
         do j = 1, ny_tot
            do i = 1, nx_tot
               ms%tracers(ms%idx_salinity)%c(i, j, 1) = 1.0e-4_wp*34.5_wp
               ms%tracers(ms%idx_temperature)%c(i, j, 1) = 1.0e-4_wp*1.0_wp
            end do
         end do
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == 0, "a filler at h*c_live must pass the scan")
         if (allocated(error)) exit checks

         ! Break it again, sweep, and require clean.
         do j = 1, ny_tot
            do i = 1, nx_tot
               ms%tracers(ms%idx_salinity)%c(i, j, 1) = 3.0e-3_wp
            end do
         end do
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == nx_tot*ny_tot, "a drifted filler must fail the scan")
         if (allocated(error)) exit checks
         call ms%enforce_vanished_content(nx_tot, ny_tot)
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == 0, "the scan must be clean after enforcement")
         if (allocated(error)) exit checks
         call check(error, worst == 0.0_wp, "clean scan reports worst = 0")
      end block checks

      call ms%destroy()
   end subroutine test_i1_scan

   subroutine test_i1_scan_zero_crossing(error)
      !! The tripwire must not flag ROUND-OFF where the pooled concentration
      !! passes through zero.  Bed filler + live layer pool to a mixed
      !! temperature of ~1e-17 degC out of terms of magnitude 1 (the OM3 1 deg
      !! z* run tripped on exactly this: |hTr - h*c_live| = 5e-20, tolerance
      !! relative to |h*c_live| ~ 0).  After the sweep the scan must be clean
      !! while a REAL violation (filler at 0 under a 1 degC donor) is still seen.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      integer :: nx_tot, ny_tot, i, j, k, n_bad
      real(wp) :: worst

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               ms%h_layer(i, j, k) = 10.0_wp
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 345.0_wp
               ms%tracers(ms%idx_temperature)%c(i, j, k) = 10.0_wp
            end do
            ms%h_layer(i, j, 1) = 1.0e-4_wp
            ms%tracers(ms%idx_salinity)%c(i, j, 1) = 1.0e-4_wp*34.5_wp
            ! Pool {k=1 filler, k=2 live} sums to one ulp of 1 (~1.1e-16).
            ms%tracers(ms%idx_temperature)%c(i, j, 2) = 1.0_wp
            ms%tracers(ms%idx_temperature)%c(i, j, 1) = &
               -(1.0_wp - epsilon(1.0_wp)/2.0_wp)
         end do
      end do

      checks: block
         call ms%enforce_vanished_content(nx_tot, ny_tot)
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == 0, &
                    "a pool mixing through zero must pass the scan after the sweep")
         if (allocated(error)) exit checks

         ! A real violation under the same column scale is still caught.
         do j = 1, ny_tot
            do i = 1, nx_tot
               ms%tracers(ms%idx_temperature)%c(i, j, 1) = 0.0_wp
               ms%tracers(ms%idx_temperature)%c(i, j, 2) = 1.0e1_wp
            end do
         end do
         call ms%scan_vanished_content(nx_tot, ny_tot, n_bad, worst)
         call check(error, n_bad == nx_tot*ny_tot, &
                    "an empty filler under a 1 degC donor must still fail the scan")
      end block checks

      call ms%destroy()
   end subroutine test_i1_scan_zero_crossing

   subroutine test_remap_constancy(error)
      !! THE property I1′ exists for, at the level of one remap: a UNIFORM
      !! column stays uniform through the remap — fillers at the bed, in the
      !! interior and at the top, a source filler set that differs from the
      !! target's, and the source fillers already carrying `h·c` (the I1′
      !! state).  Every target layer, filler or live, must read `c0` to
      !! round-off.  Now structural: the concentration form never leaves the
      !! (here one-point) source range; unchanged assertion.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp), parameter :: C0 = 35.0_wp
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), worst
      integer :: k

      h_old = [1.0e-4_wp, 12.0_wp, 1.0e-4_wp, 30.0_wp, 25.0_wp, 5.0e-5_wp, 18.0_wp, 1.0e-4_wp]
      h_new = [3.0_wp, 9.0_wp, 1.0e-4_wp, 28.0_wp, 27.0_wp, 1.0e-4_wp, 18.9_wp - 5.0e-5_wp, &
               1.0e-4_wp]
      ! Make the two grids span the same column exactly.
      h_new(2) = h_new(2) + (sum(h_old) - sum(h_new))
      hTr = C0*h_old

      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      worst = 0.0_wp
      do k = 1, NZ
         worst = max(worst, abs(hTr(k)/h_new(k)/C0 - 1.0_wp))
      end do
      call check(error, worst <= 1.0e-13_wp, &
                 "a uniform column must stay uniform through the remap, fillers included")
   end subroutine test_remap_constancy

   subroutine test_top_fillers_k_top(error)
      !! An ice-shelf column under `Z_FIXED`: fillers inside the ice ABOVE the
      !! live column and fillers below the bed partial cell.  The top run has
      !! no live layer above it, so its donor is the topmost live layer
      !! (`k_top`); the bed run's donor is the lowest live layer.  Enforce,
      !! then check the donors explicitly — the top fillers must carry the
      !! `k_top` concentration, not the bed's.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NX = 2, NY = 2, NZ = 7
      real(wp), parameter :: HF = 1.0e-4_wp
      real(wp) :: c_top, c_bot, qs(NZ)
      integer :: nx_tot, ny_tot, i, j
      logical :: ok

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      !   k:   1 bed filler, 2 bed partial, 3-4 live, 5 = k_top, 6-7 in the ice
      do j = 1, ny_tot
         do i = 1, nx_tot
            ms%h_layer(i, j, :) = [HF, 3.0_wp, 40.0_wp, 40.0_wp, 17.0_wp, HF, HF]
            ! Salinity: 34.6 deep, 34.2 at k_top; fillers empty (the I1 state).
            ms%tracers(ms%idx_salinity)%c(i, j, :) = [0.0_wp, 3.0_wp*34.6_wp, &
                                                      40.0_wp*34.5_wp, 40.0_wp*34.4_wp, &
                                                      17.0_wp*34.2_wp, 0.0_wp, 0.0_wp]
         end do
      end do

      call ms%enforce_vanished_content(nx_tot, ny_tot)

      ok = .true.
      do j = 1, ny_tot
         do i = 1, nx_tot
            qs = ms%tracers(ms%idx_salinity)%c(i, j, :)
            c_top = qs(5)/17.0_wp
            c_bot = qs(2)/3.0_wp
            if (abs(qs(6)/HF - c_top) > 1.0e-12_wp*c_top) ok = .false.
            if (abs(qs(7)/HF - c_top) > 1.0e-12_wp*c_top) ok = .false.
            if (abs(qs(1)/HF - c_bot) > 1.0e-12_wp*c_bot) ok = .false.
            ! The donors moved by at most the filler share of an empty pool.
            if (abs(c_top - 34.2_wp) > 34.2_wp*2.0_wp*HF/17.0_wp) ok = .false.
            if (abs(c_bot - 34.6_wp) > 34.6_wp*HF/3.0_wp) ok = .false.
            ! The interior is untouched.
            if (qs(3) /= 40.0_wp*34.5_wp .or. qs(4) /= 40.0_wp*34.4_wp) ok = .false.
         end do
      end do
      call check(error, ok, "top fillers take k_top's concentration, bed fillers the bed layer's")
      call ms%destroy()
   end subroutine test_top_fillers_k_top

   subroutine test_host_twin(error)
      !! `enforce_vanished_content_host` — the setup-time twin the seed calls
      !! before `enter_data` — must produce the same content as the sweep
      !! (both run the one included `rdb_vl_merge_content`), and I1′.  To
      !! ROUND-OFF, not bitwise: on the offload build the sweep runs on the
      !! device, whose FMA contraction of `h*c` / `Σq − Σh*c` differs from
      !! the host's (measured: bitwise equal on gfortran, NOT bitwise equal
      !! on nvfortran cc70) — the "never assert bit-zero across an
      !! FMA-contractible expression" rule.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      integer, parameter :: NX = 4, NY = 3, NZ = 6
      integer :: nx_tot, ny_tot, i, j, k, t, seed
      real(wp) :: r
      logical :: same, i1_ok

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      ms_a%nz_ml = NZ
      ms_b%nz_ml = NZ
      call ms_a%init(grid)
      call ms_b%init(grid)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      seed = 4242
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               call lcg_next(seed, r)
               if (r < 0.35_wp) then
                  ms_a%h_layer(i, j, k) = 1.0e-4_wp*r
               else
                  ms_a%h_layer(i, j, k) = 1.0_wp + 30.0_wp*r
               end if
               do t = 1, size(ms_a%tracers)
                  call lcg_next(seed, r)
                  ms_a%tracers(t)%c(i, j, k) = 10.0_wp + 25.0_wp*r
               end do
            end do
         end do
      end do
      ms_b%h_layer = ms_a%h_layer
      do t = 1, size(ms_a%tracers)
         ms_b%tracers(t)%c = ms_a%tracers(t)%c
      end do

      call ms_a%enforce_vanished_content(nx_tot, ny_tot)
      call ms_b%enforce_vanished_content_host(nx_tot, ny_tot)

      same = .true.
      i1_ok = .true.
      do t = 1, size(ms_a%tracers)
         do j = 1, ny_tot
            do i = 1, nx_tot
               do k = 1, NZ
                  if (abs(ms_a%tracers(t)%c(i, j, k) - ms_b%tracers(t)%c(i, j, k)) > &
                      1.0e-13_wp*abs(ms_a%tracers(t)%c(i, j, k))) same = .false.
               end do
               if (.not. i1p_holds(NZ, ms_b%h_layer(i, j, :), ms_b%tracers(t)%c(i, j, :))) &
                  i1_ok = .false.
            end do
         end do
      end do
      checks: block
         call check(error, same, "host twin must match the sweep to round-off")
         if (allocated(error)) exit checks
         call check(error, i1_ok, "host twin must establish I1'")
      end block checks
      call ms_a%destroy()
      call ms_b%destroy()
   end subroutine test_host_twin

   subroutine test_mismatch_fold(error)
      !! The target grid a few ulp SHORTER than the source (what the target
      !! builders hand the remap every step, here exaggerated to 1e-9
      !! relative, the precondition tolerance, so the assertion cannot pass by
      !! accident), on a column whose top fillers carry `c_live`.
      !!
      !! CHANGED with the concentration-form remap, which has no
      !! `remap_fold_filler_defect`: the sliver at the top of the source is
      !! simply not covered by any target, and the column loses exactly the
      !! water it cannot place, at ITS concentration — `c_top·δh`, the
      !! consistent answer for a tracer whose concentration is what is
      !! conserved — and nothing else.  The loss is part of `hTr_new −
      !! hTr_old`, hence of the remap budget slot.  (The fold used to hand
      !! the sliver's content to the top live layer; with the tracer stored as
      !! concentration there is no content to hand over.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), before, after, dh

      !        bed      live      live      k_top     filler     filler
      h_old = [40.0_wp, 48.0_wp, 48.0_wp, 36.0_wp, 1.0e-4_wp, 1.0e-4_wp]
      h_new = h_old
      h_new(1) = h_old(1) - 1.72e-7_wp       ! the bed absorbs the short fall
      dh = sum(h_new) - sum(h_old)
      hTr = [40.0_wp*34.6_wp, 48.0_wp*34.5_wp, 48.0_wp*34.4_wp, 36.0_wp*34.2_wp, &
             1.0e-4_wp*34.2_wp, 1.0e-4_wp*34.2_wp]
      before = col_sum(NZ, hTr)
      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      after = col_sum(NZ, hTr)
      call check(error, abs((after - before) - 34.2_wp*dh) <= 64.0_wp*epsilon(1.0_wp)*abs(before), &
                 "a short target column loses exactly c_top*dh of water")
   end subroutine test_mismatch_fold

end module test_ocean_remap_vanished
