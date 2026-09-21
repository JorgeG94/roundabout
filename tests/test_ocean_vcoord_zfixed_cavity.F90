!! `VCOORD_Z_FIXED` under a RIGID TOP — the quasi-geopotential target
!! grid an ice-shelf cavity needs (P6.2).
!!
!! Every assertion here is on the absolute GEOPOTENTIAL depth of an
!! interface, `e_K = Σ_{k >= K} target_h(k)` measured UP from the bed,
!! not merely on the column sum: a sum-only test cannot see a stack
!! placed in the wrong half of the column, which is the failure mode
!! this slice exists to remove (Yung, Hallberg, Adcroft & Morrison 2026,
!! JAMES 18, e2025MS005645, Fig. 1b — quasi-z layers are geopotential and
!! VANISH where they outcrop into the ice base).
!!
!! The bit-identity case is the hinge of the whole phase and is asserted
!! with `==`, against an independent transcription of the pre-cavity
!! algorithm held in `reference_pre_cavity_column` below.
module test_ocean_vcoord_zfixed_cavity
   use, intrinsic :: iso_fortran_env, only: real64
   use rdb_constants, only: wp, H_VANISHED, VCOORD_Z_FIXED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_zfixed_cavity_tests

   ! Geometry shared by most cases: a 1000 m bed, nz = 20 ⇒ 50 m nominal
   ! spacing, a 500 m flat lid ⇒ a 500 m water column.
   integer, parameter :: NZ = 20
   real(wp), parameter :: H_REF = 1000.0_wp
   real(wp), parameter :: H_NOM = H_REF/real(NZ, wp)
   real(wp), parameter :: H_MIN = 1.0e-4_wp
   real(wp), parameter :: TOP_FRAC = 0.1_wp
      !! Must track `Z_FIXED_TOP_PARTIAL_FRAC` in `rdb_ocean_vcoord`.

contains

   subroutine collect_ocean_vcoord_zfixed_cavity_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("z_top_zero_is_bit_identical", test_bit_identity), &
                  new_unittest("flat_lid_interface_depths", test_flat_lid_depths), &
                  new_unittest("sloping_lid_first_live_steps", test_sloping_lid), &
                  new_unittest("partial_top_cell_sliver_merges_below", test_sliver_merge), &
                  new_unittest("eta_lands_in_a_live_layer", test_eta_placement), &
                  new_unittest("both_ends_vanish_sum_exact", test_both_ends_vanish) &
                  ]
   end subroutine collect_ocean_vcoord_zfixed_cavity_tests

   subroutine build(vc, grid, z_top_val, total_h_val, eta_val, target_col, nz_use)
      !! Build a one-column-per-test target stack: uniform `z_top`,
      !! `total_h` and `eta` over the whole (ghosted) array, then read
      !! the interior column back.
      type(ocean_vcoord_t), intent(inout) :: vc
      type(hgrid_t), intent(inout) :: grid
      real(wp), intent(in) :: z_top_val, total_h_val, eta_val
      integer, intent(in) :: nz_use
      real(wp), intent(out) :: target_col(nz_use)
      real(wp), allocatable :: total_h(:, :), eta(:, :)
      integer :: k

      call grid%init(4, 4, 1, 1000.0_wp, 1000.0_wp)
      call vc%init(grid, nz_ml=nz_use)
      vc%coord_type = VCOORD_Z_FIXED
      vc%z_fixed_h_ref = H_REF
      vc%zstar_h_min = H_MIN
      vc%z_top = z_top_val
      allocate (total_h(grid%nx_total, grid%ny_total), source=total_h_val)
      allocate (eta(grid%nx_total, grid%ny_total), source=eta_val)
      call vc%compute_target_h(total_h, eta)
      do k = 1, nz_use
         target_col(k) = vc%target_h(2, 2, k)
      end do
      deallocate (total_h, eta)
   end subroutine build

   pure subroutine reference_pre_cavity_column(column, nz, h_nominal, h_min, target_col)
      !! Verbatim transcription of the Z_FIXED target walk as it stood
      !! BEFORE the rigid top existed (`rdb_ocean_vcoord.F90`, the
      !! `case (VCOORD_Z_FIXED)` body of `compute_target_h_impl`).  The
      !! bit-identity test asserts the live code reproduces this to the
      !! last bit when `z_top = 0`.
      integer, intent(in) :: nz
      real(wp), intent(in) :: column, h_nominal, h_min
      real(wp), intent(out) :: target_col(nz)
      integer :: k
      real(wp) :: z_below_loc, z_above_nominal_loc
      z_below_loc = column
      do k = 1, nz
         z_above_nominal_loc = real(nz - k, wp)*h_nominal
         if (z_above_nominal_loc > z_below_loc - h_min) then
            target_col(k) = h_min
            z_below_loc = z_below_loc - h_min
         else
            target_col(k) = z_below_loc - z_above_nominal_loc
            z_below_loc = z_above_nominal_loc
         end if
      end do
   end subroutine reference_pre_cavity_column

   pure function interface_depth(target_col, nz, k, bed_depth) result(depth)
      !! Geopotential depth (m, positive down) of the interface ABOVE
      !! layer `k`: the bed less everything stacked below that interface.
      integer, intent(in) :: nz, k
      real(wp), intent(in) :: target_col(nz), bed_depth
      real(wp) :: depth
      integer :: kk
      depth = bed_depth
      do kk = 1, k
         depth = depth - target_col(kk)
      end do
   end function interface_depth

   ! ------------------------------------------------------------------
   ! (1) THE BIT-IDENTITY HINGE.  `z_top = 0` must reproduce the
   !     pre-cavity arithmetic exactly — same expressions, same operands
   !     — over a deep column (nothing vanishes), a shallow one (the bed
   !     side vanishes) and a non-zero free surface.
   ! ------------------------------------------------------------------
   subroutine test_bit_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: got(NZ), want(NZ)
      real(wp) :: cases_h(4), cases_eta(4)
      integer :: c, k
      cases_h = [1000.0_wp, 400.0_wp, 137.3_wp, 1000.0_wp]
      cases_eta = [0.0_wp, 0.0_wp, 0.0_wp, 0.37_wp]
      cases: block
         do c = 1, 4
            call build(vc, grid, 0.0_wp, cases_h(c), cases_eta(c), got, NZ)
            call vc%destroy()
            call reference_pre_cavity_column(cases_h(c) + cases_eta(c), NZ, &
                                             H_NOM, H_MIN, want)
            do k = 1, NZ
               call check(error, got(k) == want(k), &
                          "z_top=0 must be BIT-identical to the pre-cavity target")
               if (allocated(error)) exit cases
            end do
         end do
      end block cases
   end subroutine test_bit_identity

   ! ------------------------------------------------------------------
   ! (2) A FLAT 500 m LID OVER A 1000 m BED, nz = 20, h_ref = 1000.
   !     The upper 500 m of the nominal stack is inside the ice: layers
   !     k = 11..20 are inert fillers, k = 1..10 are live and sit at
   !     their OPEN-OCEAN geopotential depths, and the column sums to the
   !     500 m of water exactly.
   ! ------------------------------------------------------------------
   subroutine test_flat_lid_depths(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: col(NZ), total, tol, filler_debt
      integer :: k
      checks: block
         call build(vc, grid, 500.0_wp, 500.0_wp, 0.0_wp, col, NZ)
         call vc%destroy()

         ! Tolerance: the interface depth is an nz-term running sum of
         ! O(h_ref) terms, so a few ulp of h_ref per term.
         tol = 4.0_wp*real(NZ, wp)*epsilon(1.0_wp)*H_REF

         ! The ten layers inside the ice are inert fillers — at, not
         ! merely below, the filler thickness, and <= H_VANISHED so every
         ! operator that gates on the marker skips them.
         do k = 11, NZ
            call check(error, abs(col(k) - H_MIN) <= tol, &
                       "layers 11..20 must be inert fillers at zstar_h_min")
            if (allocated(error)) exit checks
            call check(error, col(k) <= H_VANISHED, &
                       "a top-side filler must be <= H_VANISHED")
            if (allocated(error)) exit checks
         end do

         ! The nine fully-submerged layers carry the nominal spacing and
         ! their interfaces sit at the OPEN-OCEAN geopotential depths
         ! (nz - k)*h_nominal, independent of the ice.
         do k = 1, 9
            call check(error, abs(col(k) - H_NOM) <= tol, &
                       "a fully submerged layer must carry h_nominal")
            if (allocated(error)) exit checks
            call check(error, abs(interface_depth(col, NZ, k, H_REF) &
                                  - real(NZ - k, wp)*H_NOM) <= tol, &
                       "live interfaces must sit at their geopotential depths")
            if (allocated(error)) exit checks
         end do

         ! k = 10 is the partial TOP cell: the nominal 50 m less the ten
         ! fillers' h_min debt, so the stack closes on the water column.
         filler_debt = real(NZ - 10, wp)*H_MIN
         call check(error, abs(col(10) - (H_NOM - filler_debt)) <= tol, &
                    "the partial top cell pays for the fillers above it")
         if (allocated(error)) exit checks

         total = sum(col)
         call check(error, abs(total - 500.0_wp) <= tol, &
                    "the target must sum to the 500 m water column")
      end block checks
   end subroutine test_flat_lid_depths

   ! ------------------------------------------------------------------
   ! (3) A SLOPING DRAFT: the first live layer index STEPS across the
   !     domain as the draft crosses nominal levels, the live interfaces
   !     stay at their geopotential depths on every column, and every
   !     column conserves.
   ! ------------------------------------------------------------------
   subroutine test_sloping_lid(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: col(NZ), draft, tol
      integer :: c, k, k_live, k_live_prev, n_step
      logical :: seen_step
      checks: block
         tol = 4.0_wp*real(NZ, wp)*epsilon(1.0_wp)*H_REF
         k_live_prev = 0
         n_step = 0
         seen_step = .false.
         ! Sweep the draft down through two nominal levels.
         do c = 0, 20
            draft = 380.0_wp + 6.0_wp*real(c, wp)      ! 380 .. 500 m
            call build(vc, grid, draft, H_REF - draft, 0.0_wp, col, NZ)
            call vc%destroy()

            ! Column conservation on every column of the slope.
            call check(error, abs(sum(col) - (H_REF - draft)) <= tol, &
                       "a sloping-draft column must sum to its water column")
            if (allocated(error)) exit checks

            ! First live layer, counting down from the top.
            k_live = 1
            do k = NZ, 1, -1
               if (col(k) > H_VANISHED) then
                  k_live = k
                  exit
               end if
            end do
            ! Every LIVE interface strictly below the first live layer is
            ! geopotential: it does not move as the ice base moves.
            do k = 1, k_live - 1
               call check(error, abs(interface_depth(col, NZ, k, H_REF) &
                                     - real(NZ - k, wp)*H_NOM) <= tol, &
                          "interior interfaces must not follow the ice base")
               if (allocated(error)) exit checks
            end do
            ! The partial top cell never falls below the minimum
            ! partial-cell fraction, and never exceeds the nominal
            ! spacing by more than one merged sliver.
            call check(error, col(k_live) <= (1.0_wp + TOP_FRAC)*H_NOM + tol, &
                       "the partial top cell cannot exceed h_nominal + a sliver")
            if (allocated(error)) exit checks
            call check(error, col(k_live) >= TOP_FRAC*H_NOM - tol, &
                       "the partial top cell must clear the sliver floor")
            if (allocated(error)) exit checks

            if (c > 0) then
               if (k_live /= k_live_prev) then
                  n_step = n_step + 1
                  seen_step = .true.
                  call check(error, k_live == k_live_prev - 1, &
                             "the first live layer must step DOWN by one")
                  if (allocated(error)) exit checks
               end if
            end if
            k_live_prev = k_live
         end do
         call check(error, seen_step, &
                    "a 120 m draft sweep over a 50 m spacing must step the index")
         if (allocated(error)) exit checks
         call check(error, n_step >= 2, "the sweep must cross at least two levels")
      end block checks
   end subroutine test_sloping_lid

   ! ------------------------------------------------------------------
   ! (4) THE SLIVER RULE.  A draft placed just above a nominal level
   !     would cut a partial top cell thinner than
   !     Z_FIXED_TOP_PARTIAL_FRAC*h_nominal.  It must NOT ship as a
   !     sliver: the index becomes another filler and the water merges
   !     into the layer BELOW — the mirror of the bed side, where a
   !     would-be sub-h_min cell collapses and hands its water UP.
   ! ------------------------------------------------------------------
   subroutine test_sliver_merge(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: col_sliver(NZ), col_clean(NZ), tol, draft_sliver, draft_clean
      integer :: k, k_live_sliver, k_live_clean
      checks: block
         tol = 4.0_wp*real(NZ, wp)*epsilon(1.0_wp)*H_REF
         ! Layer 11 spans the nominal depth band [450, 500] m.  A draft
         ! of 480 m cuts it to 20 m — a legal partial top cell.  A draft
         ! of 499 m would cut it to 1 m, below the 0.1*50 = 5 m floor, so
         ! layer 11 becomes a filler too and its 1 m merges into layer
         ! 10, which then carries 50 + 1 m.
         draft_sliver = 499.0_wp
         draft_clean = 480.0_wp

         call build(vc, grid, draft_sliver, H_REF - draft_sliver, 0.0_wp, col_sliver, NZ)
         call vc%destroy()
         call build(vc, grid, draft_clean, H_REF - draft_clean, 0.0_wp, col_clean, NZ)
         call vc%destroy()

         k_live_sliver = 1
         do k = NZ, 1, -1
            if (col_sliver(k) > H_VANISHED) then
               k_live_sliver = k
               exit
            end if
         end do
         k_live_clean = 1
         do k = NZ, 1, -1
            if (col_clean(k) > H_VANISHED) then
               k_live_clean = k
               exit
            end if
         end do

         call check(error, k_live_clean == 11, &
                    "a 480 m draft must leave layer 11 as the partial top cell")
         if (allocated(error)) exit checks
         call check(error, abs(col_clean(11) - (20.0_wp - real(NZ - 11, wp)*H_MIN)) <= tol, &
                    "the legal partial top cell is the cut less the filler debt")
         if (allocated(error)) exit checks
         call check(error, k_live_sliver == 10, &
                    "a 499 m draft must merge the 1 m sliver and drop to layer 10")
         if (allocated(error)) exit checks
         ! The merged layer is the nominal spacing PLUS the sliver, less
         ! the fillers' debt: it is thicker than nominal, not thinner.
         call check(error, col_sliver(10) > H_NOM, &
                    "the sliver must be absorbed by the layer below")
         if (allocated(error)) exit checks
         call check(error, abs(col_sliver(10) - (H_NOM + 1.0_wp &
                                                 - real(NZ - 10, wp)*H_MIN)) <= tol, &
                    "the merged layer must carry exactly h_nominal + the sliver")
         if (allocated(error)) exit checks
         call check(error, abs(sum(col_sliver) - (H_REF - draft_sliver)) <= tol, &
                    "the merge must not break column conservation")
      end block checks
   end subroutine test_sliver_merge

   ! ------------------------------------------------------------------
   ! (5) WHERE η GOES.  The free-surface anomaly must land in a LIVE
   !     layer — never in a vanished k = nz filler.  The closing rule is
   !     the pre-cavity one, unchanged: the LOWEST live layer (the
   !     partial BOTTOM cell) absorbs η and the bed-side filler debt,
   !     while the top end pays for its own fillers.
   ! ------------------------------------------------------------------
   subroutine test_eta_placement(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: col0(NZ), col1(NZ), tol, d_eta
      integer :: k
      checks: block
         tol = 4.0_wp*real(NZ, wp)*epsilon(1.0_wp)*H_REF
         d_eta = 0.75_wp
         call build(vc, grid, 500.0_wp, 500.0_wp, 0.0_wp, col0, NZ)
         call vc%destroy()
         call build(vc, grid, 500.0_wp, 500.0_wp, d_eta, col1, NZ)
         call vc%destroy()

         ! Everything except the lowest live layer is untouched by η.
         do k = 2, NZ
            call check(error, abs(col1(k) - col0(k)) <= tol, &
                       "eta must not disturb any layer but the lowest live one")
            if (allocated(error)) exit checks
         end do
         call check(error, abs((col1(1) - col0(1)) - d_eta) <= tol, &
                    "eta must land in the lowest LIVE layer, in full")
         if (allocated(error)) exit checks
         ! And in particular NOT in the top filler.
         call check(error, col1(NZ) <= H_VANISHED, &
                    "k = nz must stay a vanished filler when eta /= 0")
         if (allocated(error)) exit checks
         call check(error, abs(sum(col1) - (500.0_wp + d_eta)) <= tol, &
                    "the column must still sum to H + eta")
      end block checks
   end subroutine test_eta_placement

   ! ------------------------------------------------------------------
   ! (6) BOTH ENDS VANISH.  A near-grounding-line column: a 500 m draft
   !     over a 520 m bed leaves 20 m of water, so the stack is fillers
   !     below the bed, fillers inside the ice, and ONE live layer cut at
   !     both ends.  The single closing rule must still conserve.
   ! ------------------------------------------------------------------
   subroutine test_both_ends_vanish(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: col(NZ), tol
      integer :: k, n_live
      checks: block
         tol = 4.0_wp*real(NZ, wp)*epsilon(1.0_wp)*H_REF
         call build(vc, grid, 500.0_wp, 20.0_wp, 0.0_wp, col, NZ)
         call vc%destroy()

         n_live = 0
         do k = 1, NZ
            if (col(k) > H_VANISHED) n_live = n_live + 1
            call check(error, col(k) > 0.0_wp, &
                       "no target thickness may be zero or negative")
            if (allocated(error)) exit checks
         end do
         call check(error, n_live == 1, &
                    "a 20 m column on a 50 m spacing carries ONE live layer")
         if (allocated(error)) exit checks
         call check(error, col(10) > H_VANISHED, &
                    "the live layer must be the one the ice base cuts, k = 10")
         if (allocated(error)) exit checks
         call check(error, abs(sum(col) - 20.0_wp) <= tol, &
                    "both ends vanishing must still conserve the column")
      end block checks
   end subroutine test_both_ends_vanish

end module test_ocean_vcoord_zfixed_cavity
