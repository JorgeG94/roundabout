!! The vanished-layer content rule for the ALE tracer remap.
!!
!! Invariant **I1**: after a remap, `h <= H_VANISHED ⇒ hTr = 0` for every
!! registered tracer, with the per-column content sum conserved to round-off.
!!
!! The guard used to be ONE-SIDED — the READ side floored (`h_old <= H_VANISHED
!! ⇒ c_old = 0`) while the WRITE side did not (`hTr_new = c_new·h_new` for
!! EVERY target layer, filler included). So the remap parked content in a
!! sub-threshold layer and deleted it one step later with no budget
!! contributor: the day-16 salt/heat break on
!! `validation_examples/ocean/isomip_plus/ocean0_idealised_zfixed.nml`.
!!
!! Coverage here:
!!   - the section-G reproducer: content parked in a 1e-4 filler plus an
!!     advective increment, one remap, column sum conserved (this FAILS on the
!!     one-sided guard — it loses exactly the parked content);
!!   - I1 holds after a remap, for fillers at the bed, at the top and in the
!!     interior;
!!   - pseudo-random columns with filler runs in all three places: conservation
!!     + I1 + no content ever left in a filler;
!!   - **bit-identity**: a column with no sub-threshold layer takes textually
!!     the old code path — asserted as EXACT equality against the old
!!     expression, not a tolerance;
!!   - the degenerate column with no live layer at all (land / grounded);
!!   - the merge helper itself: conservation, idempotence, no-op on a live
!!     column;
!!   - end-to-end through `ocean_apply_ale_remap_centres` with a
!!     `VCOORD_ZSTAR_FULL` filler column: the `*_budget_remap` contributor
!!     telescopes to zero per column, which is what makes it a valid leak
!!     detector rather than a record of the leak.
module test_ocean_remap_vanished
   use rdb_constants, only: wp, REMAP_PPM, REMAP_PLM, H_VANISHED, NZ_STACK_MAX, &
                            VCOORD_ZSTAR_FULL
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_remap_column, only: remap_column
   use rdb_ocean_remap, only: ocean_remap_tracer_column, &
                              remap_merge_vanished_content, &
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
                  new_unittest("section_g_parked_content_conserved", test_section_g_reproducer), &
                  new_unittest("i1_holds_after_remap", test_i1_after_remap), &
                  new_unittest("random_columns_with_fillers", test_random_columns), &
                  new_unittest("bit_identical_without_fillers", test_bit_identity), &
                  new_unittest("dead_column_stays_dead", test_dead_column), &
                  new_unittest("merge_helper_properties", test_merge_helper), &
                  new_unittest("remap_budget_telescopes_with_fillers", test_budget_telescopes) &
                  ]
   end subroutine collect_ocean_remap_vanished_tests

   ! ---------------------------------------------------------------- helpers

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
      !! §G.4 of the diagnosis, verbatim: `nz = 4`, bed filler at `h = 1e-4`
      !! holding 34.5 PSU of parked content (3.45e-3 PSU.m), plus the step's
      !! advective increment (8.889e-4), remapped onto a target where the bed
      !! layer has just become live (1.758e-3 m).
      !!
      !! One-sided guard: loses 4.3413e-3 PSU.m — the whole parked content.
      !! Two-sided: conserved to round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ)
      real(wp) :: before, after

      h_old = [1.0e-4_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      h_new = [1.758272e-3_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      hTr = [3.45e-3_wp + 8.888926e-4_wp, 690.0_wp, 690.0_wp, 0.0_wp]

      before = col_sum(NZ, hTr)
      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      after = col_sum(NZ, hTr)

      checks: block
         call check(error, abs(after - before) <= 1.0e-12_wp*abs(before), &
                    "section-G column: parked content must survive the remap")
         if (allocated(error)) exit checks
         ! And it must not be parked AGAIN: the top filler is still 1e-4.
         call check(error, hTr(NZ) == 0.0_wp, &
                    "section-G column: no content may be written into a filler")
      end block checks
   end subroutine test_section_g_reproducer

   subroutine test_i1_after_remap(error)
      !! I1 after a remap, with fillers at the bed, in the interior and at the
      !! top, and a target whose filler set differs from the source's.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ)
      real(wp) :: before, after
      integer :: k
      logical :: ok

      !            bed filler  live   filler  live   live   filler  live   top filler
      h_old = [1.0e-4_wp, 12.0_wp, 1.0e-4_wp, 30.0_wp, 25.0_wp, 5.0e-5_wp, 18.0_wp, 1.0e-4_wp]
      h_new = [3.0_wp, 9.0_wp, 1.0e-4_wp, 28.0_wp, 27.0_wp, 1.0e-4_wp, 18.9_wp, 1.0e-4_wp]
      do k = 1, NZ
         hTr(k) = (34.0_wp + 0.1_wp*real(k, wp))*h_old(k)
      end do
      before = col_sum(NZ, hTr)

      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
      after = col_sum(NZ, hTr)

      checks: block
         call check(error, abs(after - before) <= 1.0e-12_wp*abs(before), &
                    "mixed-filler column: content conserved")
         if (allocated(error)) exit checks
         ok = .true.
         do k = 1, NZ
            if (h_new(k) <= H_VANISHED .and. hTr(k) /= 0.0_wp) ok = .false.
         end do
         call check(error, ok, "I1: a vanished target layer must hold no content")
      end block checks
   end subroutine test_i1_after_remap

   subroutine test_random_columns(error)
      !! Pseudo-random columns, filler runs placed at the bed, at the top and in
      !! the interior in turn.  Every column: conservation to 1e-12 relative and
      !! I1 on the target grid.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 12, NTRIAL = 60
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ)
      real(wp) :: before, after, worst_rel, r, live_new, defect
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
            ! Source filler carrying PARKED content — the failure mode.
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

         before = col_sum(NZ, hTr)
         call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
         after = col_sum(NZ, hTr)
         worst_rel = max(worst_rel, abs(after - before)/abs(before))
         do k = 1, NZ
            if (h_new(k) <= H_VANISHED .and. hTr(k) /= 0.0_wp) i1_ok = .false.
         end do
      end do

      checks: block
         call check(error, worst_rel <= 1.0e-12_wp, &
                    "random filler columns: content conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, i1_ok, "random filler columns: I1 holds on every target")
      end block checks
   end subroutine test_random_columns

   subroutine test_bit_identity(error)
      !! A column with NO sub-threshold layer must take textually the old code
      !! path.  Reference = the pre-fix expression, spelled out here; equality
      !! is asserted EXACTLY (`==`), not to a tolerance.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10, NTRIAL = 40
      real(wp) :: h_old(NZ), h_new(NZ), hTr(NZ), hTr_ref(NZ)
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX)
      real(wp) :: r
      integer :: seed, trial, k
      logical :: identical

      seed = 991
      identical = .true.
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

         ! --- the ONE-SIDED (pre-fix) algorithm, verbatim
         do k = 1, NZ
            if (h_old(k) > H_VANISHED) then
               c_old(k) = hTr_ref(k)/h_old(k)
            else
               c_old(k) = 0.0_wp
            end if
         end do
         call remap_column(REMAP_PPM, NZ, h_old, h_new, c_old(1:NZ), c_new(1:NZ))
         do k = 1, NZ
            hTr_ref(k) = c_new(k)*h_new(k)
         end do

         call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PPM)
         do k = 1, NZ
            if (hTr(k) /= hTr_ref(k)) identical = .false.
         end do
      end do

      call check(error, identical, &
                 "no sub-threshold layer ⇒ bit-identical to the one-sided remap")
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
      !! `remap_merge_vanished_content` on its own: conserves the column sum,
      !! zeroes every sub-threshold layer, is idempotent, and is a textual
      !! no-op on a column with no filler.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp) :: h(NZ_STACK_MAX), q(NZ_STACK_MAX), q2(NZ_STACK_MAX)
      real(wp) :: before, after
      integer :: k
      logical :: ok

      h = 0.0_wp
      q = 0.0_wp
      h(1:NZ) = [1.0e-4_wp, 1.0e-4_wp, 10.0_wp, 1.0e-4_wp, 20.0_wp, 1.0e-4_wp]
      q(1:NZ) = [1.0_wp, 2.0_wp, 300.0_wp, 4.0_wp, 600.0_wp, 6.0_wp]
      before = col_sum(NZ, q)

      call remap_merge_vanished_content(NZ, h, q)
      after = col_sum(NZ, q)

      checks: block
         call check(error, abs(after - before) <= 1.0e-13_wp*abs(before), &
                    "merge helper conserves the column sum")
         if (allocated(error)) exit checks
         ok = .true.
         do k = 1, NZ
            if (h(k) <= H_VANISHED .and. q(k) /= 0.0_wp) ok = .false.
         end do
         call check(error, ok, "merge helper zeroes every vanished layer")
         if (allocated(error)) exit checks

         ! Nearest live layer: 1 + 2 go UP into layer 3; 4 goes up into 5;
         ! 6 has nothing above, so it comes back DOWN into 5.
         call check(error, q(3) == 303.0_wp, "bed fillers merge into the layer above")
         if (allocated(error)) exit checks
         call check(error, q(5) == 610.0_wp, "top filler merges back into the layer below")
         if (allocated(error)) exit checks

         q2 = q
         call remap_merge_vanished_content(NZ, h, q2)
         ok = .true.
         do k = 1, NZ
            if (q2(k) /= q(k)) ok = .false.
         end do
         call check(error, ok, "merge helper is idempotent")
         if (allocated(error)) exit checks

         h(1:NZ) = [5.0_wp, 5.0_wp, 10.0_wp, 5.0_wp, 20.0_wp, 5.0_wp]
         q(1:NZ) = [1.0_wp, 2.0_wp, 300.0_wp, 4.0_wp, 600.0_wp, 6.0_wp]
         q2 = q
         call remap_merge_vanished_content(NZ, h, q2)
         ok = .true.
         do k = 1, NZ
            if (q2(k) /= q(k)) ok = .false.
         end do
         call check(error, ok, "no filler ⇒ merge helper is a textual no-op")
      end block checks
   end subroutine test_merge_helper

   subroutine test_budget_telescopes(error)
      !! End-to-end through the orchestrator on a ZSTAR_FULL column whose bed
      !! layers vanish, with content deliberately PARKED in a bed filler (the
      !! state the one-sided guard leaves behind).  The `*_budget_remap`
      !! contributor must telescope to zero per column — nothing leaves the
      !! column, so the accumulator stays a leak detector — and the column
      !! content must be conserved.
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

      ! Two bed fillers holding PARKED content, four live layers above.
      do j = 1, ny_tot
         do i = 1, nx_tot
            do k = 1, NZ
               if (k <= 2) then
                  ms%h_layer(i, j, k) = 1.0e-4_wp
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*1.0e-3_wp
               else
                  ms%h_layer(i, j, k) = (H_BED - 2.0e-4_wp)/4.0_wp
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*ms%h_layer(i, j, k)
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
               col_before(i, j) = col_before(i, j) + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
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
               s = s + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
               b = b + ms%salt_budget_remap(i, j, k)
               if (ms%h_layer(i, j, k) <= H_VANISHED .and. &
                   ms%tracers(ms%idx_salinity)%hTr(i, j, k) /= 0.0_wp) i1_ok = .false.
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
         call check(error, i1_ok, "I1 after the orchestrator's tracer remap")
      end block checks

      deallocate (bt_eta, bt_H_ref, h_bed_2d, col_before)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_budget_telescopes

end module test_ocean_remap_vanished
