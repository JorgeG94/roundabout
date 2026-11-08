!! Vanishing / refilling layers through the ALE tracer remap, CONCENTRATION form.
!!
!! The tracer prognostic `c` is a concentration (White & Adcroft 2008).  A
!! layer at or below `H_VANISHED` (a "filler": z-fixed / z* bed layers, ice
!! shelf lids, land) carries NO content of its own and just keeps its last
!! concentration: the remap reads that `c` as a cell value of the stencil
!! (limited by its neighbours) and a waking filler takes a concentration from
!! the water it sits in.  There is no donor map, no pool merge and no
!! tail sweep any more; the old invariant I1′ (`hTr = h·c_donor` on a filler)
!! and its enforcement point, scan, merge helper and host twin are gone, and
!! so are the tests of them.  What the remap itself guarantees, asserted here:
!!   - the content `Σ h·c` of a column is conserved to round-off once the
!!     recorded bounding amount is added back, and every target concentration
!!     (filler targets included) lies inside the source range;
!!   - a layer that wakes up (bed layer going live) takes its neighbour's
!!     concentration, never `0`;
!!   - a UNIFORM column stays uniform, fillers included;
!!   - a column with no sub-threshold layer equals the plain `remap_column`
!!     pipeline exactly, except where the bounding acted, and there by the
!!     recorded amount;
!!   - the degenerate column with no live layer (land / grounded) invents
!!     nothing;
!!   - a target a few ulp short of the source loses exactly `c·δh` of water;
!!   - end-to-end through `ocean_apply_ale_remap_centres` with a
!!     `VCOORD_ZSTAR_FULL` filler column: the live-layer content `h·c` is
!!     conserved, the `*_budget_remap` contributor telescopes to zero per
!!     column, and a uniform tracer stays uniform in every layer.
module test_ocean_remap_vanished
   use rdb_constants, only: wp, REMAP_PPM, REMAP_PLM, H_VANISHED, NZ_STACK_MAX, &
                            VCOORD_ZSTAR_FULL
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_remap_column, only: remap_column
   use rdb_ocean_remap, only: ocean_remap_tracer_column, &
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
                  new_unittest("waking_bed_layer_takes_neighbour_concentration", test_section_g_reproducer), &
                  new_unittest("fillers_bounded_conserved_after_remap", test_i1_after_remap), &
                  new_unittest("random_columns_with_fillers", test_random_columns), &
                  new_unittest("no_filler_equals_plain_pipeline_but_recorded_clamp", test_bit_identity), &
                  new_unittest("dead_column_invents_nothing", test_dead_column), &
                  new_unittest("remap_budget_telescopes_with_fillers", test_budget_telescopes), &
                  new_unittest("uniform_column_stays_uniform_through_remap", test_remap_constancy), &
                  new_unittest("mismatched_column_loses_exactly_c_dh", test_mismatch_fold) &
                  ]
   end subroutine collect_ocean_remap_vanished_tests

   ! ---------------------------------------------------------------- helpers

   pure function content(n, h, c) result(s)
      !! Column content `Σ h·c`.
      integer, intent(in) :: n
      real(wp), intent(in) :: h(:), c(:)
      real(wp) :: s
      integer :: k
      s = 0.0_wp
      do k = 1, n
         s = s + h(k)*c(k)
      end do
   end function content

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

   pure subroutine fillers_hold_neighbour_c(n, h, c)
      !! Give every filler (`h <= H_VANISHED`) the concentration of the
      !! nearest live layer ABOVE it, else BELOW (what a filler that has just
      !! been vanished by the remap would still hold).  A stale value that is
      !! physically sensible, so the bounds asserted below are meaningful.
      integer, intent(in) :: n
      real(wp), intent(in) :: h(:)
      real(wp), intent(inout) :: c(:)
      integer :: k, kd, kk
      real(wp) :: cn(n)
      cn = c(1:n)
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
         if (kd /= 0) cn(k) = c(kd)
      end do
      c(1:n) = cn
   end subroutine fillers_hold_neighbour_c

   pure subroutine lcg_next(seed, r)
      !! Reproducible pseudo-random in [0,1) -- a 31-bit Lehmer generator, so the
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
      !! `nz = 4` uniform S = 34.5 with a bed filler at `h = 1e-4` that WAKES
      !! UP (target 1.758e-3 m) and a top filler that stays one.  The waking
      !! layer must come out at 34.5 (the water it joins), not `0`, and the
      !! column content must be conserved.  (Replaces the section-G
      !! "parked excess is not read" case, which only existed because fillers
      !! carried content.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp), parameter :: C0 = 34.5_wp
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), before, after

      h_old = [1.0e-4_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      h_new = [1.758272e-3_wp, 20.0_wp, 20.0_wp, 1.0e-4_wp]
      h_new(2) = h_new(2) + (sum(h_old) - sum(h_new))
      c = C0
      before = content(NZ, h_old, c)
      call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM)
      after = content(NZ, h_new, c)

      checks: block
         call check(error, abs(c(1) - C0) <= 1.0e-13_wp*C0, &
                    "a waking bed layer takes the neighbour's concentration, not 0")
         if (allocated(error)) exit checks
         call check(error, all(abs(c - C0) <= 1.0e-13_wp*C0), &
                    "uniform column with a waking filler stays uniform")
         if (allocated(error)) exit checks
         call check(error, abs(after - before) <= 1.0e-12_wp*abs(before), &
                    "waking-layer column: content conserved")
      end block checks
   end subroutine test_section_g_reproducer

   subroutine test_i1_after_remap(error)
      !! Fillers at the bed, in the interior and at the top, and a target whose
      !! filler set differs from the source's.  Asserted on the OUTPUT of the
      !! concentration-form remap: (a) content `Σ h·c` conserved to round-off,
      !! (b) every target -- filler targets included -- carries a
      !! concentration inside the source range (not `0`, not a made-up value),
      !! (c) the bounding amount is round-off.  (Replaces the old assertion of
      !! invariant I1′, `hTr = h·c_donor`, on the output.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), c_new(NZ), clamp(NZ)
      real(wp) :: before, after, c_lo, c_hi
      integer :: k

      !            bed filler  live   filler  live   live   filler  live   top filler
      h_old = [1.0e-4_wp, 12.0_wp, 1.0e-4_wp, 30.0_wp, 25.0_wp, 5.0e-5_wp, 18.0_wp, 1.0e-4_wp]
      h_new = [3.0_wp, 9.0_wp, 1.0e-4_wp, 28.0_wp, 27.0_wp, 1.0e-4_wp, 18.9_wp, 1.0e-4_wp]
      h_new(7) = sum(h_old) - (sum(h_new) - h_new(7))
      do k = 1, NZ
         c(k) = 34.0_wp + 0.1_wp*real(k, wp)
      end do
      call fillers_hold_neighbour_c(NZ, h_old, c)
      c_lo = minval(c)
      c_hi = maxval(c)
      before = content(NZ, h_old, c)

      call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM, c_new_out=c_new, clamp_out=clamp)
      after = content(NZ, h_new, c)

      checks: block
         call check(error, abs(after - before) <= 1.0e-12_wp*abs(before), &
                    "mixed-filler column: content conserved")
         if (allocated(error)) exit checks
         call check(error, all(c_new >= c_lo) .and. all(c_new <= c_hi), &
                    "every target concentration, filler targets included, lies in the source range")
         if (allocated(error)) exit checks
         call check(error, abs(col_sum(NZ, clamp)) <= 1.0e-12_wp*abs(before), &
                    "the recorded bounding amount is round-off")
      end block checks
   end subroutine test_i1_after_remap

   subroutine test_random_columns(error)
      !! Pseudo-random columns, filler runs placed at the bed, at the top and
      !! in the interior in turn.  Every column: content conserved to 1e-12
      !! relative and every target concentration (filler targets included)
      !! inside the source range.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 12, NTRIAL = 60
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), c_new(NZ)
      real(wp) :: before, after, worst_rel, r, live_new, defect, c_lo, c_hi
      integer :: seed, trial, k, place, k0, k1
      logical :: in_range

      seed = 20260921
      worst_rel = 0.0_wp
      in_range = .true.
      do trial = 1, NTRIAL
         do k = 1, NZ
            call lcg_next(seed, r)
            h_old(k) = 1.0_wp + 40.0_wp*r
            call lcg_next(seed, r)
            h_new(k) = 1.0_wp + 40.0_wp*r
            call lcg_next(seed, r)
            c(k) = 30.0_wp + 8.0_wp*r
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
            h_old(k) = 1.0e-4_wp
            ! Target: half the trials keep it a filler, half wake it up.
            if (mod(trial, 2) == 0) then
               h_new(k) = 1.0e-4_wp
            else
               call lcg_next(seed, r)
               h_new(k) = 1.0e-3_wp + r
            end if
         end do
         ! The target grid must span the SAME column as the source -- what
         ! `compute_target_h` guarantees in production.  Absorb the mismatch
         ! into the live target layers, proportionally.
         live_new = 0.0_wp
         do k = 1, NZ
            if (h_new(k) > H_VANISHED) live_new = live_new + h_new(k)
         end do
         defect = col_sum(NZ, h_old) - col_sum(NZ, h_new)
         do k = 1, NZ
            if (h_new(k) > H_VANISHED) h_new(k) = h_new(k)*(1.0_wp + defect/live_new)
         end do

         call fillers_hold_neighbour_c(NZ, h_old, c)
         c_lo = minval(c)
         c_hi = maxval(c)
         before = content(NZ, h_old, c)
         call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM, c_new_out=c_new)
         after = content(NZ, h_new, c)
         worst_rel = max(worst_rel, abs(after - before)/abs(before))
         if (any(c_new < c_lo) .or. any(c_new > c_hi)) in_range = .false.
      end do

      checks: block
         call check(error, worst_rel <= 1.0e-12_wp, &
                    "random filler columns: content conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, in_range, "random filler columns: every target concentration inside the source range")
      end block checks
   end subroutine test_random_columns

   subroutine test_bit_identity(error)
      !! A column with NO sub-threshold layer equals the plain pipeline
      !! `remap_column` on the stored concentration, spelled out here, EXACTLY
      !! (`==`) wherever the bounding did not act; where it did (an overshoot
      !! of a few ulp on a target that coincides with a source extremum) the
      !! content difference is the recorded `clamp_out` amount to round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10, NTRIAL = 40
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), c_ref(NZ), clamp(NZ)
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX)
      real(wp) :: r
      integer :: seed, trial, k
      logical :: identical, clamp_recorded

      seed = 991
      identical = .true.
      clamp_recorded = .true.
      do trial = 1, NTRIAL
         do k = 1, NZ
            call lcg_next(seed, r)
            h_old(k) = 0.5_wp + 30.0_wp*r
            call lcg_next(seed, r)
            h_new(k) = 0.5_wp + 30.0_wp*r
            call lcg_next(seed, r)
            c(k) = 25.0_wp + 12.0_wp*r
         end do

         ! --- the plain pipeline, verbatim
         c_old(1:NZ) = c
         call remap_column(REMAP_PPM, NZ, h_old, h_new, c_old(1:NZ), c_new(1:NZ))
         c_ref = c_new(1:NZ)

         call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM, clamp_out=clamp)
         do k = 1, NZ
            if (clamp(k) == 0.0_wp) then
               if (c(k) /= c_ref(k)) identical = .false.
            else
               if (abs(h_new(k)*(c(k) - c_ref(k)) - clamp(k)) > 4.0_wp*epsilon(1.0_wp)*abs(h_new(k)*c(k))) &
                  clamp_recorded = .false.
            end if
         end do
      end do

      checks: block
         call check(error, identical, &
                    "no sub-threshold layer => equal to the plain pipeline wherever the bounding is inactive")
         if (allocated(error)) exit checks
         call check(error, clamp_recorded, &
                    "where the bounding acted, the content moved by exactly the recorded amount")
      end block checks
   end subroutine test_bit_identity

   subroutine test_dead_column(error)
      !! Land / fully grounded column: every layer a filler.  Nothing must be
      !! invented: a uniform stale concentration stays exactly that (finite, no
      !! `0`, no NaN), and a non-uniform one stays inside its own range.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ)

      h_old = 1.0e-4_wp
      h_new = 1.0e-4_wp
      c = 7.0_wp
      call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM)
      checks: block
         call check(error, all(abs(c - 7.0_wp) <= 1.0e-13_wp*7.0_wp), &
                    "dead column: a uniform stale concentration is kept")
         if (allocated(error)) exit checks
         c = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
         call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM)
         call check(error, all(c >= 1.0_wp - 1.0e-13_wp) .and. all(c <= 5.0_wp + 1.0e-13_wp), &
                    "dead column: the remap must not manufacture concentration outside the stale range")
      end block checks
   end subroutine test_dead_column

   subroutine test_budget_telescopes(error)
      !! End-to-end through the orchestrator on a ZSTAR_FULL column whose bed
      !! layers vanish, the two bed fillers holding their donor's
      !! concentration.  (1) With a vertical gradient the column content
      !! `Σ h·S` must be conserved and the `*_budget_remap` contributor must
      !! telescope to zero per column -- nothing leaves the column, so the
      !! accumulator stays a leak detector.  (2) With a uniform S every layer,
      !! filler included, must read S0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 6
      real(wp), parameter :: H_BED = 120.0_wp, S0 = 34.5_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :), h_bed_2d(:, :)
      real(wp), allocatable :: col_before(:, :)
      real(wp) :: s, b, worst_content, worst_budget, worst_uniform, c_lo, c_hi
      integer :: nx_tot, ny_tot, i, j, k, kk, kb
      logical :: in_range, gradient

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      allocate (h_bed_2d(nx_tot, ny_tot), source=H_BED)
      allocate (bt_eta(nx_tot, ny_tot), source=0.0_wp)
      allocate (bt_H_ref(nx_tot, ny_tot), source=H_BED)
      allocate (col_before(nx_tot, ny_tot), source=0.0_wp)

      worst_content = 0.0_wp
      worst_budget = 0.0_wp
      worst_uniform = 0.0_wp
      in_range = .true.
      c_lo = S0
      c_hi = S0 + 0.05_wp*real(NZ, wp)
      do k = 1, 2     ! pass 1: vertical gradient; pass 2: uniform
         gradient = (k == 1)
         ms%nz_ml = NZ
         call ms%init(grid)
         ! Two bed fillers, four live layers above.
         do j = 1, ny_tot
            do i = 1, nx_tot
               do kk = 1, NZ
                  if (kk <= 2) then
                     ms%h_layer(i, j, kk) = 1.0e-4_wp
                  else
                     ms%h_layer(i, j, kk) = (H_BED - 2.0e-4_wp)/4.0_wp
                  end if
                  if (gradient) then
                     ms%tracers(ms%idx_salinity)%c(i, j, kk) = S0 + 0.05_wp*real(max(kk, 3), wp)
                  else
                     ms%tracers(ms%idx_salinity)%c(i, j, kk) = S0
                  end if
               end do
            end do
         end do

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 20.0_wp
         vc%zstar_n_surf = 2
         call vc%build_zref_full(h_bed_2d)

         do j = 1, ny_tot
            do i = 1, nx_tot
               col_before(i, j) = content(NZ, ms%h_layer(i, j, :), ms%tracers(ms%idx_salinity)%c(i, j, :))
            end do
         end do

         call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref, method=REMAP_PLM)

         do j = 1, ny_tot
            do i = 1, nx_tot
               s = content(NZ, ms%h_layer(i, j, :), ms%tracers(ms%idx_salinity)%c(i, j, :))
               b = 0.0_wp
               do kb = 1, NZ
                  b = b + ms%salt_budget_remap(i, j, kb)
                  if (ms%tracers(ms%idx_salinity)%c(i, j, kb) < c_lo - 1.0e-12_wp .or. &
                      ms%tracers(ms%idx_salinity)%c(i, j, kb) > c_hi + 1.0e-12_wp) in_range = .false.
                  if (.not. gradient) worst_uniform = max(worst_uniform, &
                                                          abs(ms%tracers(ms%idx_salinity)%c(i, j, kb) - S0))
               end do
               if (gradient) then
                  worst_content = max(worst_content, abs(s - col_before(i, j)))
                  worst_budget = max(worst_budget, abs(b))
               end if
            end do
         end do
         call vc%destroy()
         call ms%destroy()
      end do

      checks: block
         call check(error, worst_content <= 1.0e-10_wp*abs(col_before(1, 1)), &
                    "ZSTAR_FULL filler column: content conserved through the orchestrator")
         if (allocated(error)) exit checks
         call check(error, worst_budget <= 1.0e-10_wp*abs(col_before(1, 1)), &
                    "remap budget contributor telescopes to zero per column")
         if (allocated(error)) exit checks
         call check(error, in_range, "gradient column: every layer's concentration stays inside the source range")
         if (allocated(error)) exit checks
         call check(error, worst_uniform <= 1.0e-12_wp*S0, &
                    "uniform S stays uniform in every layer, fillers included")
      end block checks

      deallocate (bt_eta, bt_H_ref, h_bed_2d, col_before)
   end subroutine test_budget_telescopes

   subroutine test_remap_constancy(error)
      !! THE property that motivates the concentration form, at the level of one
      !! remap: a UNIFORM column stays uniform -- fillers at the bed, in the
      !! interior and at the top, a source filler set that differs from the
      !! target's.  Every target layer, filler or live, must read `c0` to
      !! round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp), parameter :: C0 = 35.0_wp
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), worst

      h_old = [1.0e-4_wp, 12.0_wp, 1.0e-4_wp, 30.0_wp, 25.0_wp, 5.0e-5_wp, 18.0_wp, 1.0e-4_wp]
      h_new = [3.0_wp, 9.0_wp, 1.0e-4_wp, 28.0_wp, 27.0_wp, 1.0e-4_wp, 18.9_wp - 5.0e-5_wp, &
               1.0e-4_wp]
      ! Make the two grids span the same column exactly.
      h_new(2) = h_new(2) + (sum(h_old) - sum(h_new))
      c = C0

      call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM)
      worst = maxval(abs(c/C0 - 1.0_wp))
      call check(error, worst <= 1.0e-13_wp, &
                 "a uniform column must stay uniform through the remap, fillers included")
   end subroutine test_remap_constancy

   subroutine test_mismatch_fold(error)
      !! The target grid a few ulp SHORTER than the source (what the target
      !! builders hand the remap every step, here exaggerated), on a column
      !! whose top fillers hold `c_top`.  The column loses exactly the water it
      !! cannot place, at ITS concentration -- `c_top·δh` -- and nothing else.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp) :: h_old(NZ), h_new(NZ), c(NZ), before, after, dh

      !        bed      live      live      k_top     filler     filler
      h_old = [40.0_wp, 48.0_wp, 48.0_wp, 36.0_wp, 1.0e-4_wp, 1.0e-4_wp]
      h_new = h_old
      h_new(1) = h_old(1) - 1.72e-7_wp       ! the bed absorbs the short fall
      dh = sum(h_new) - sum(h_old)
      c = [34.6_wp, 34.5_wp, 34.4_wp, 34.2_wp, 34.2_wp, 34.2_wp]
      before = content(NZ, h_old, c)
      call ocean_remap_tracer_column(NZ, h_old, h_new, c, REMAP_PPM)
      after = content(NZ, h_new, c)
      call check(error, abs((after - before) - 34.2_wp*dh) <= 64.0_wp*epsilon(1.0_wp)*abs(before), &
                 "a short target column loses exactly c_top*dh of water")
   end subroutine test_mismatch_fold

end module test_ocean_remap_vanished
