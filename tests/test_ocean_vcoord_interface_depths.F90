!! Interface-DEPTH gate for every vertical-coordinate family the ocean
!! path can build a geometric target for.
!!
!! **Why this suite exists.**  Every other vcoord test in the tree asserts
!! that `target_h` is positive and that `Σ_k target_h = H + η`.  Both
!! properties hold *identically* for a stack laid in the wrong half of the
!! column — so a sum-only test cannot see a coordinate that puts its fine
!! resolution 500 m too deep, nor one whose whole column has collapsed
!! into the bed layer.  This suite asserts the **absolute geopotential
!! depth of every target interface**,
!!
!!     e(K) = −b + Σ_{k ≤ K} target_h(k),      e(0) = −b (the bed),
!!
!! against a hand-derived analytic table, for each family and for three
!! column geometries:
!!
!!   1. open ocean, flat bed, `η = 0`   — interfaces at their analytic depths
!!   2. open ocean, sloping bed          — bed-side vanishing where the family
!!                                         has it; live layers at the right depths
!!   3. the SAME columns with the column TOP displaced to `z = −z_top`
!!      (a flat lid at 500 m over a 1000 m bed: the water column is 500 m
!!      thick and occupies `−1000 ≤ z ≤ −500`).
!!
!! Case 3 is the one that matters for the ice-shelf-cavity work.  For a
!! **sigma-like** family the analytic expectation is simply "the live
!! column divided proportionally", and the families meet it.  For a
!! **z-like** family the analytic expectation is that layers whose nominal
!! range lies ABOVE `−z_top` vanish against the TOP while the live layers
!! keep their open-ocean geopotential depths:
!!
!!     e(K) = −max((NZ − K)·h_nominal, z_top)
!!
!! — geopotential depth, clipped at the rigid top.
!!
!! **`VCOORD_Z_FIXED` now meets that (P6.2).**  The slot carries a
!! per-column top depth `ocean_vcoord_t%z_top` (`metrics%z_draft` under a
!! cavity, `0` otherwise), so the builder finally knows where the column
!! starts, and `z_fixed_under_a_lid_vanishes_against_the_top` asserts the
!! table above — it is the flipped twin of what was a `documents_*` case
!! pinning the old, wrong placement.  Its companion
!! `z_fixed_z_top_zero_reproduces_the_old_placement` keeps the other half
!! of the contract: at `z_top = 0` every target thickness is reproduced
!! BIT-for-bit, so the family's pre-cavity answers cannot move.
!!
!! `VCOORD_ZSTAR_FULL` is still carried as a `documents_*` case: its
!! per-column reference table is built from the TRUE bed and has not been
!! taught the draft (that is P6.11).  It asserts the CURRENT, wrong
!! placement, prints the measured depth beside the analytic one in the
!! failure message, and is flagged in-line.  The suite is green; the
!! remaining defect is pinned and loud.
!!
!! Families covered — every branch of `ocean_vcoord_compute_target_h`:
!! LAGRANGIAN (no target), EULERIAN_Z, SIGMA, ZSTAR (shares the SIGMA
!! branch), ZSIGMA, ZSTAR_SIGMA, ZSTAR_FULL, Z_FIXED.  The two
!! density-space families (`VCOORD_RHO`, `VCOORD_HYCOM`) come in through
!! the sibling `compute_target_h_rho` wrapper and place their interfaces
!! on prescribed potential densities, not at geometric depths — their
!! placement gate is `test_ocean_vcoord_rho` /
!! `test_ocean_vcoord_hycom`, which assert exactly that.
!!
!! Tolerances: `TOL_EXACT` (1e-9 m) where the arithmetic is exact, and
!! `TOL_FILLER` (2e-3 m) wherever a column carries inert filler layers —
!! twice the whole-column filler budget `nz · zstar_h_min = 1e-3 m`.  The
!! defects this suite guards against are 100–500 m.
module test_ocean_vcoord_interface_depths
   use rdb_constants, only: wp, H_VANISHED, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_interface_depths_tests

   integer, parameter :: NZ = 10
      !! Layer count for every case.  Bottom-up: `k = 1` is the bed,
      !! `k = NZ` the surface layer (the repo-wide ROMS convention).
   integer, parameter :: NX = 2, NY = 1
      !! Interior extents; `nghost = 1` ⇒ `nx_total = 4`, `ny_total = 3`.
      !! The four i-columns carry the sloping-bed profile `B_SLOPE`.
   real(wp), parameter :: Z_FIXED_H_REF = 1000.0_wp
      !! `VCOORD_Z_FIXED` reference stack depth ⇒ `h_nominal = 100 m`.
   real(wp), parameter :: H_NOMINAL = Z_FIXED_H_REF/real(NZ, wp)
   real(wp), parameter :: B_FLAT = 1000.0_wp
      !! Flat-bed depth for cases 1 and 3.
   real(wp), parameter :: B_SLOPE(NX + 2) = [1000.0_wp, 650.0_wp, 250.0_wp, 120.0_wp]
      !! Case-2 bed profile.  Chosen so the ZSTAR_FULL reference table is
      !! exactly representable: with `zstar_h_surf_target = 100` and the
      !! auto `n_surf = max(1, nz/3) = 3`, a bed of 1000 m gives a uniform
      !! 100 m table and 650 m gives 100/100/100 then 50 m coarse layers,
      !! while 250 m and 120 m fall in the too-shallow reserve branch.
   real(wp), parameter :: Z_TOP_LID = 500.0_wp
      !! Case-3 rigid-top depth: the lid sits at `z = −500` over the
      !! `B_FLAT = 1000 m` bed, so the live column is 500 m thick.
   real(wp), parameter :: TOL_EXACT = 1.0e-9_wp
   real(wp), parameter :: TOL_FILLER = 2.0e-3_wp

contains

   subroutine collect_ocean_vcoord_interface_depths_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("lagrangian_builds_no_target", test_lagrangian), &
                  new_unittest("eulerian_z_open_ocean_depths", test_eulerian_open), &
                  new_unittest("eulerian_z_under_a_lid_divides_live_column", test_eulerian_lid), &
                  new_unittest("sigma_open_ocean_depths", test_sigma_open), &
                  new_unittest("sigma_under_a_lid_divides_live_column", test_sigma_lid), &
                  new_unittest("zstar_lite_depths_match_sigma", test_zstar_lite), &
                  new_unittest("zstar_sigma_open_ocean_depths", test_zstar_sigma_open), &
                  new_unittest("zstar_sigma_under_a_lid_divides_live_column", test_zstar_sigma_lid), &
                  new_unittest("documents_zsigma_dimensionless_zref_collapse", test_zsigma_collapse), &
                  new_unittest("zstar_full_open_ocean_depths", test_zstar_full_open), &
                  new_unittest("documents_zstar_full_band_anchored_at_z0", test_zstar_full_lid), &
                  new_unittest("z_fixed_open_ocean_depths", test_z_fixed_open), &
                  new_unittest("z_fixed_under_a_lid_vanishes_against_the_top", &
                               test_z_fixed_lid), &
                  new_unittest("z_fixed_z_top_zero_reproduces_the_old_placement", &
                               test_z_fixed_z_top_zero), &
                  new_unittest("inert_fillers_stay_at_or_below_h_vanished", test_filler_contract) &
                  ]
   end subroutine collect_ocean_vcoord_interface_depths_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   pure subroutine column_depths(target_h_col, bed_depth, e_iface)
      !! Geopotential depth of every interface of one column, measured
      !! UP from the bed: `e(0) = −b`, `e(K) = e(K−1) + target_h(K)`.
      !! `e(NZ)` is the top of the water column (`z = η` in the open
      !! ocean, `z = −z_top + η` under a rigid lid).
      real(wp), intent(in) :: target_h_col(NZ)
      real(wp), intent(in) :: bed_depth
         !! Bed depth (m, positive-down), so the bed is at `z = −bed_depth`.
      real(wp), intent(out) :: e_iface(0:NZ)
      integer :: k
      e_iface(0) = -bed_depth
      do k = 1, NZ
         e_iface(k) = e_iface(k - 1) + target_h_col(k)
      end do
   end subroutine column_depths

   pure subroutine proportional_depths(bed_depth, column_total, e_want)
      !! Analytic interfaces of a SIGMA-like family with uniform `dsig`:
      !! the live column is divided proportionally, so the stack runs
      !! from the bed up to `−bed_depth + column_total`.
      real(wp), intent(in) :: bed_depth, column_total
      real(wp), intent(out) :: e_want(0:NZ)
      integer :: k
      do k = 0, NZ
         e_want(k) = -bed_depth + column_total*real(k, wp)/real(NZ, wp)
      end do
   end subroutine proportional_depths

   pure subroutine z_fixed_depths_from_column_top(bed_depth, column_total, e_want)
      !! Analytic interfaces of TODAY's `VCOORD_Z_FIXED`: the nominal
      !! stack hangs from the COLUMN TOP (`z = −bed_depth + column_total`)
      !! and the leftover layers vanish against the BED.
      !!
      !!   e(NZ − m) = (column top) − m·h_nominal,   m = 0 … n_live − 1
      !!   e(K)      = −bed_depth                     for K ≤ NZ − n_live
      !!
      !! with `n_live = ceiling(column_total / h_nominal)`.
      real(wp), intent(in) :: bed_depth, column_total
      real(wp), intent(out) :: e_want(0:NZ)
      integer :: k, m, n_live
      real(wp) :: z_column_top
      z_column_top = -bed_depth + column_total
      n_live = min(NZ, ceiling(column_total/H_NOMINAL - TOL_EXACT))
      do m = 0, n_live - 1
         e_want(NZ - m) = z_column_top - real(m, wp)*H_NOMINAL
      end do
      do k = 0, NZ - n_live
         e_want(k) = -bed_depth
      end do
   end subroutine z_fixed_depths_from_column_top

   pure subroutine z_fixed_depths_under_rigid_top(bed_depth, z_top, e_want)
      !! Analytic interfaces of `VCOORD_Z_FIXED` under a RIGID TOP (P6.2):
      !! the nominal stack is GEOPOTENTIAL and is clipped at the top,
      !!
      !!   e(K) = -max((NZ - K)*h_nominal, z_top)
      !!
      !! so every interface deeper than the lid keeps its open-ocean depth
      !! and every one that would sit inside the ice collapses onto the
      !! lid.  The fillers stacked under the lid displace the topmost
      !! interfaces by at most `NZ*zstar_h_min = 1e-3 m`, which is why
      !! this table is checked at `TOL_FILLER`.
      real(wp), intent(in) :: bed_depth, z_top
      real(wp), intent(out) :: e_want(0:NZ)
      integer :: k
      do k = 0, NZ
         e_want(k) = -max(real(NZ - k, wp)*H_NOMINAL, z_top)
      end do
      ! The bed is the bed whatever the lid does.
      e_want(0) = -bed_depth
   end subroutine z_fixed_depths_under_rigid_top

   pure subroutine z_fixed_reference_pre_cavity(column_total, target_col)
      !! Verbatim transcription of the `VCOORD_Z_FIXED` target walk as it
      !! stood BEFORE the rigid top existed.  Used by the `z_top = 0`
      !! twin to assert bit-for-bit reproduction with `==`; an analytic
      !! depth table cannot make that statement.
      real(wp), intent(in) :: column_total
      real(wp), intent(out) :: target_col(NZ)
      integer :: k
      real(wp) :: z_below_loc, z_above_nominal_loc
      z_below_loc = column_total
      do k = 1, NZ
         z_above_nominal_loc = real(NZ - k, wp)*H_NOMINAL
         if (z_above_nominal_loc > z_below_loc - 1.0e-4_wp) then
            target_col(k) = 1.0e-4_wp
            z_below_loc = z_below_loc - 1.0e-4_wp
         else
            target_col(k) = z_below_loc - z_above_nominal_loc
            z_below_loc = z_above_nominal_loc
         end if
      end do
   end subroutine z_fixed_reference_pre_cavity

   subroutine check_depths(error, family, e_got, e_want, tol)
      !! Compare a measured interface-depth column against its analytic
      !! table, naming the family and printing BOTH depths for the first
      !! interface that misses.  Not `pure`: `check` allocates `error`.
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: family
      real(wp), intent(in) :: e_got(0:NZ), e_want(0:NZ)
      real(wp), intent(in) :: tol
      character(len=256) :: msg
      integer :: k
      do k = 0, NZ
         if (abs(e_got(k) - e_want(k)) > tol) then
            write (msg, '(a,i0,a,f14.6,a,f14.6,a,es10.3,a)') &
               family//": interface e(", k, ") at z = ", e_got(k), &
               " m, analytic z = ", e_want(k), " m (miss ", &
               abs(e_got(k) - e_want(k)), " m)"
            call check(error, .false., trim(msg))
            return
         end if
      end do
   end subroutine check_depths

   subroutine check_column_sum(error, family, target_h_col, column_total)
      !! The invariant every other vcoord test already covers, kept here
      !! as a cheap secondary so a depth failure can be told apart from a
      !! conservation failure.  Not `pure`: `check` allocates `error`.
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: family
      real(wp), intent(in) :: target_h_col(NZ)
      real(wp), intent(in) :: column_total
      call check(error, abs(sum(target_h_col) - column_total) <= TOL_FILLER, &
                 family//": sum(target_h) should equal H + eta")
   end subroutine check_column_sum

   subroutine make_slot(vc, grid, coord_type)
      !! Fresh `ocean_vcoord_t` on an `NX x NY` grid with one ghost ring.
      !! Not `pure`: `init` allocates.
      type(ocean_vcoord_t), intent(out) :: vc
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: coord_type
      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = coord_type
      vc%zstar_h_min = 1.0e-4_wp
      vc%zstar_h_surf_target = 100.0_wp
      vc%z_fixed_h_ref = Z_FIXED_H_REF
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 100.0_wp
   end subroutine make_slot

   ! ------------------------------------------------------------------
   ! LAGRANGIAN — no geometric target at all
   ! ------------------------------------------------------------------

   subroutine test_lagrangian(error)
      !! `VCOORD_LAGRANGIAN` returns before touching `target_h` (the ALE
      !! remap is a no-op for it), so it has no interface depths to gate:
      !! the target array keeps whatever it held.  Datum-free by
      !! construction — it is the one family a displaced column top
      !! cannot mis-place.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      call make_slot(vc, grid, VCOORD_LAGRANGIAN)
      vc%target_h = -7.0_wp
      total_h = B_FLAT
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))
      call check(error, all(abs(vc%target_h + 7.0_wp) < TOL_EXACT), &
                 "LAGRANGIAN: compute_target_h must not write target_h")
      call vc%destroy()
   end subroutine test_lagrangian

   ! ------------------------------------------------------------------
   ! EULERIAN_Z — H·dsig, eta dropped
   ! ------------------------------------------------------------------

   subroutine test_eulerian_open(error)
      !! Cases 1 + 2.  `H·dsig` with uniform `dsig` is a proportional
      !! division of the column, so the interfaces run bed → surface in
      !! equal steps on a flat bed AND on a sloping one.  Note the family
      !! drops `η` — asserted here with `η = 0` so the two agree.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: i
      checks: block
         call make_slot(vc, grid, VCOORD_EULERIAN_Z)
         eta = 0.0_wp
         ! --- case 1: flat bed
         total_h = B_FLAT
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         call proportional_depths(B_FLAT, B_FLAT, e_want)
         call check_depths(error, "EULERIAN_Z flat bed", e_got, e_want, TOL_EXACT)
         if (allocated(error)) exit checks
         ! --- case 2: sloping bed
         do i = 1, grid%nx_total
            total_h(i, :) = B_SLOPE(i)
         end do
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         do i = 1, grid%nx_total
            call column_depths(vc%target_h(i, 2, :), B_SLOPE(i), e_got)
            call proportional_depths(B_SLOPE(i), B_SLOPE(i), e_want)
            call check_depths(error, "EULERIAN_Z sloping bed", e_got, e_want, TOL_EXACT)
            if (allocated(error)) exit checks
            call check_column_sum(error, "EULERIAN_Z", vc%target_h(i, 2, :), B_SLOPE(i))
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_eulerian_open

   subroutine test_eulerian_lid(error)
      !! Case 3.  EULERIAN_Z is a stretched SIGMA with `η` dropped — it is
      !! NOT a geopotential coordinate despite the name — so under a rigid
      !! lid it divides the live column proportionally and lands exactly
      !! where the analytic sigma-like expectation puts it.  Nothing for
      !! P6.2 to change here.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      call make_slot(vc, grid, VCOORD_EULERIAN_Z)
      total_h = B_FLAT - Z_TOP_LID
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))
      call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
      call proportional_depths(B_FLAT, B_FLAT - Z_TOP_LID, e_want)
      call check_depths(error, "EULERIAN_Z under a lid", e_got, e_want, TOL_EXACT)
      call vc%destroy()
   end subroutine test_eulerian_lid

   ! ------------------------------------------------------------------
   ! SIGMA / ZSTAR (one branch)
   ! ------------------------------------------------------------------

   subroutine test_sigma_open(error)
      !! Cases 1 + 2 for `VCOORD_SIGMA`: `(H + η)·dsig`, a proportional
      !! division of the live column on any bed.  Run with `η = 1.5 m` so
      !! the free-surface anchor is exercised: the whole stack lifts by η.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp), parameter :: ETA0 = 1.5_wp
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: i
      checks: block
         call make_slot(vc, grid, VCOORD_SIGMA)
         eta = ETA0
         total_h = B_FLAT
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         call proportional_depths(B_FLAT, B_FLAT + ETA0, e_want)
         call check_depths(error, "SIGMA flat bed", e_got, e_want, TOL_EXACT)
         if (allocated(error)) exit checks
         call check(error, abs(e_got(NZ) - ETA0) < TOL_EXACT, &
                    "SIGMA: the top interface must sit at z = eta")
         if (allocated(error)) exit checks
         do i = 1, grid%nx_total
            total_h(i, :) = B_SLOPE(i)
         end do
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         do i = 1, grid%nx_total
            call column_depths(vc%target_h(i, 2, :), B_SLOPE(i), e_got)
            call proportional_depths(B_SLOPE(i), B_SLOPE(i) + ETA0, e_want)
            call check_depths(error, "SIGMA sloping bed", e_got, e_want, TOL_EXACT)
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_sigma_open

   subroutine test_sigma_lid(error)
      !! Case 3.  SIGMA is terrain-following at BOTH ends: handed a 500 m
      !! live column it divides it proportionally, which IS the analytic
      !! expectation for a sigma-like family under a rigid lid.  This is
      !! the row that makes sigma the only coordinate the cavity accepts
      !! today, and the control leg of the coordinate study.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      call make_slot(vc, grid, VCOORD_SIGMA)
      total_h = B_FLAT - Z_TOP_LID
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))
      call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
      call proportional_depths(B_FLAT, B_FLAT - Z_TOP_LID, e_want)
      call check_depths(error, "SIGMA under a lid", e_got, e_want, TOL_EXACT)
      if (allocated(error)) then
         call vc%destroy()
         return
      end if
      call check(error, abs(e_got(NZ) + Z_TOP_LID) < TOL_EXACT, &
                 "SIGMA under a lid: the top interface must sit at z = -z_top")
      call vc%destroy()
   end subroutine test_sigma_lid

   subroutine test_zstar_lite(error)
      !! `VCOORD_ZSTAR` shares the SIGMA branch on the ocean path, so its
      !! interface depths must be BIT-IDENTICAL to sigma's in every one of
      !! the three geometries.  Pinning that here keeps the "z*-lite is
      !! sigma in this (H, η) form" claim in `CLAUDE.md` honest.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc_s, vc_z
      type(hgrid_t) :: grid_s, grid_z
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      integer :: i
      checks: block
         call make_slot(vc_s, grid_s, VCOORD_SIGMA)
         call make_slot(vc_z, grid_z, VCOORD_ZSTAR)
         eta = 0.75_wp
         do i = 1, grid_s%nx_total
            total_h(i, :) = B_SLOPE(i)
         end do
         call vc_s%compute_target_h(total_h(1:grid_s%nx_total, 1:grid_s%ny_total), &
                                    eta(1:grid_s%nx_total, 1:grid_s%ny_total))
         call vc_z%compute_target_h(total_h(1:grid_z%nx_total, 1:grid_z%ny_total), &
                                    eta(1:grid_z%nx_total, 1:grid_z%ny_total))
         call check(error, all(vc_s%target_h == vc_z%target_h), &
                    "ZSTAR-lite must be bit-identical to SIGMA (sloping bed)")
         if (allocated(error)) exit checks
         total_h = B_FLAT - Z_TOP_LID
         call vc_s%compute_target_h(total_h(1:grid_s%nx_total, 1:grid_s%ny_total), &
                                    eta(1:grid_s%nx_total, 1:grid_s%ny_total))
         call vc_z%compute_target_h(total_h(1:grid_z%nx_total, 1:grid_z%ny_total), &
                                    eta(1:grid_z%nx_total, 1:grid_z%ny_total))
         call check(error, all(vc_s%target_h == vc_z%target_h), &
                    "ZSTAR-lite must be bit-identical to SIGMA (under a lid)")
      end block checks
      call vc_s%destroy()
      call vc_z%destroy()
   end subroutine test_zstar_lite

   ! ------------------------------------------------------------------
   ! ZSTAR_SIGMA — fractional rescale of the global reference table
   ! ------------------------------------------------------------------

   subroutine test_zstar_sigma_open(error)
      !! Cases 1 + 2 on the PRODUCTION `z_ref_global`.  The branch
      !! rescales the table by `column_total / z_ref_global(nz)`, so it is
      !! datum-free — and, because the only writer of `z_ref_global` on
      !! the ocean path is the uniform `k/nz` init in `ocean_vcoord_init`,
      !! the rescaled table is uniform and the family is *numerically
      !! indistinguishable from SIGMA in production*.  That is recorded
      !! here deliberately: the fractional arithmetic is sound, but with
      !! no metre-valued table ever reaching the slot the "z*-lite in
      !! deep water" half of the family does nothing.  See the ZSIGMA
      !! case below, which consumes the same table ABSOLUTELY and is
      !! therefore broken rather than merely inert.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: i
      checks: block
         call make_slot(vc, grid, VCOORD_ZSTAR_SIGMA)
         eta = 0.0_wp
         total_h = B_FLAT
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         call proportional_depths(B_FLAT, B_FLAT, e_want)
         call check_depths(error, "ZSTAR_SIGMA flat bed", e_got, e_want, TOL_EXACT)
         if (allocated(error)) exit checks
         do i = 1, grid%nx_total
            total_h(i, :) = B_SLOPE(i)
         end do
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         do i = 1, grid%nx_total
            call column_depths(vc%target_h(i, 2, :), B_SLOPE(i), e_got)
            call proportional_depths(B_SLOPE(i), B_SLOPE(i), e_want)
            call check_depths(error, "ZSTAR_SIGMA sloping bed", e_got, e_want, TOL_EXACT)
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_zstar_sigma_open

   subroutine test_zstar_sigma_lid(error)
      !! Case 3.  Purely fractional ⇒ datum-invariant: the live column is
      !! divided proportionally and nothing is mis-placed.  ZSTAR_SIGMA is
      !! therefore SAFE under a rigid top — and useless for the cavity
      !! work, because it follows both boundaries.  A second control leg.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      call make_slot(vc, grid, VCOORD_ZSTAR_SIGMA)
      total_h = B_FLAT - Z_TOP_LID
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))
      call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
      call proportional_depths(B_FLAT, B_FLAT - Z_TOP_LID, e_want)
      call check_depths(error, "ZSTAR_SIGMA under a lid", e_got, e_want, TOL_EXACT)
      call vc%destroy()
   end subroutine test_zstar_sigma_lid

   ! ------------------------------------------------------------------
   ! ZSIGMA — the dimensionless-table defect
   ! ------------------------------------------------------------------

   subroutine test_zsigma_collapse(error)
      !! `documents_*` — this asserts BROKEN behaviour on purpose.
      !!
      !! ZSIGMA's deep branch reads `z_ref_global` as a table of absolute
      !! depths **in metres** (`min(z_ref_global(k), column_total)`), but
      !! the only writer of that array on the ocean path is the
      !! dimensionless `z_ref_global(k) = k/nz` init in
      !! `ocean_vcoord_init`.  Nothing on the namelist path ever replaces
      !! it, so every z-level interval is `1/nz` **metres** and the whole
      !! column lands in the bed layer through the deficit line.
      !!
      !! Measured here on a 1000 m flat bed with `nz = 10`: nine layers of
      !! 0.1 m stacked in the top 90 cm and a 999.1 m bed layer — the
      !! interfaces `e(1) … e(9)` sit between `z = −0.9 m` and `z = 0`
      !! where a metre-valued table would have spread them through the
      !! column.  `Σ target_h = H + η` holds throughout, which is exactly
      !! why no existing test caught it.
      !!
      !! `validate_config` now refuses `vcoord_type = "zsigma"` on the
      !! ocean path for this reason.  The follow-up that fixes it (fill
      !! `z_ref_global` in metres, then make ZSIGMA the hybrid seat) must
      !! DELETE this case rather than update it.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: k
      checks: block
         call make_slot(vc, grid, VCOORD_ZSIGMA)
         total_h = B_FLAT
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         ! Today: e(0) = -b; the bed layer swallows H - (nz-1)/nz, and the
         ! nine remaining interfaces are 1/nz METRES apart under z = 0.
         e_want(0) = -B_FLAT
         do k = 1, NZ
            e_want(k) = -real(NZ - k, wp)/real(NZ, wp)
         end do
         call check_depths(error, "ZSIGMA (documented defect)", e_got, e_want, TOL_EXACT)
         if (allocated(error)) exit checks
         call check(error, vc%target_h(2, 2, 1) > 0.99_wp*B_FLAT, &
                    "ZSIGMA (documented defect): the bed layer should hold "// &
                    "essentially the whole column")
         if (allocated(error)) exit checks
         call check_column_sum(error, "ZSIGMA", vc%target_h(2, 2, :), B_FLAT)
      end block checks
      call vc%destroy()
   end subroutine test_zsigma_collapse

   ! ------------------------------------------------------------------
   ! ZSTAR_FULL — per-column reference table
   ! ------------------------------------------------------------------

   subroutine test_zstar_full_open(error)
      !! Cases 1 + 2.  The per-column table is built by `build_zref_full`
      !! from the LOCAL bed with `zstar_h_surf_target = 100 m` and the
      !! auto `n_surf = max(1, nz/3) = 3`, so the analytic interfaces are
      !!
      !!   b = 1000: 100 m fine band then a 100 m coarse fill — uniform
      !!   b =  650: 100/200/300 then 50 m coarse layers to the bed
      !!   b =  250: 100/200 then the reserve branch — two full layers,
      !!             one partial, and the rest inert fillers ON THE BED
      !!   b =  120: 100 then a 20 m partial, the rest fillers on the bed
      !!
      !! measured DOWN from `z = 0`, which here is both the column top and
      !! the geopotential datum — the two are indistinguishable in the
      !! open ocean, which is the whole reason case 3 exists.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2), h_bed(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ), z_ref_want(0:NZ)
      integer :: i, k
      checks: block
         call make_slot(vc, grid, VCOORD_ZSTAR_FULL)
         eta = 0.0_wp
         do i = 1, grid%nx_total
            total_h(i, :) = B_SLOPE(i)
            h_bed(i, :) = B_SLOPE(i)
         end do
         call vc%build_zref_full(h_bed(1:grid%nx_total, 1:grid%ny_total))
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         do i = 1, grid%nx_total
            call zstar_full_table(B_SLOPE(i), z_ref_want)
            ! Bottom-up interfaces are the top-down table reversed.
            do k = 0, NZ
               e_want(k) = -z_ref_want(NZ - k)
            end do
            call column_depths(vc%target_h(i, 2, :), B_SLOPE(i), e_got)
            call check_depths(error, "ZSTAR_FULL open ocean", e_got, e_want, TOL_FILLER)
            if (allocated(error)) exit checks
            call check_column_sum(error, "ZSTAR_FULL", vc%target_h(i, 2, :), B_SLOPE(i))
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_zstar_full_open

   pure subroutine zstar_full_table(bed_depth, z_ref_want)
      !! Analytic `z_ref(0:NZ)` (top-down, metres below `z = 0`) for
      !! `zstar_h_surf_target = 100`, `n_surf = 3`, `zstar_h_min = 1e-4`.
      !! Deep columns get a 3 x 100 m fine band plus a uniform coarse
      !! fill; columns shallower than the 300 m band fill top-down at
      !! 100 m and pin the remainder at the bed (the reserve branch),
      !! which is why TOL_FILLER rather than TOL_EXACT applies there.
      real(wp), intent(in) :: bed_depth
      real(wp), intent(out) :: z_ref_want(0:NZ)
      integer :: k
      real(wp) :: dz_coarse
      z_ref_want(0) = 0.0_wp
      if (bed_depth > 300.0_wp) then
         do k = 1, 3
            z_ref_want(k) = 100.0_wp*real(k, wp)
         end do
         dz_coarse = (bed_depth - 300.0_wp)/real(NZ - 3, wp)
         do k = 4, NZ
            z_ref_want(k) = 300.0_wp + dz_coarse*real(k - 3, wp)
         end do
      else
         do k = 1, NZ
            z_ref_want(k) = min(100.0_wp*real(k, wp), bed_depth)
         end do
      end if
      z_ref_want(NZ) = bed_depth
   end subroutine zstar_full_table

   subroutine test_zstar_full_lid(error)
      !! `documents_*` — this asserts BROKEN placement on purpose.
      !!
      !! Case 3 with a genuinely NON-uniform table: `h_surf_target = 20 m`
      !! over the 1000 m bed gives a 3 x 20 m fine band at the top of the
      !! table and 134.2857 m coarse layers below it.  The table is built
      !! from the TRUE bed (as production does,
      !! `rdb_ocean_state :: build_zref_full(state%barotropic%b)`), while
      !! the target walk is handed the 500 m live column, so
      !! `eta_loc = 500 − 1000 = −500` and the `η < 0` branch keeps the
      !! table's SHALLOW entries and vanishes its DEEP ones.
      !!
      !! Measured: the 20 m fine band lands at `z = −500 … −560`, hard
      !! against the lid, the three 134 m coarse layers fill `−560 …
      !! −962.9`, and three inert fillers sit on the bed at `−1000`.  The
      !! family does not merely anchor at `z = 0` — it INVERTS which half
      !! of the column is resolved: the table's near-surface refinement is
      !! applied 500 m too deep, and the deep water is represented by
      !! shallow-ocean spacing.
      !!
      !! P6.2 must flip this: with a per-column top depth the fine band
      !! (`z = 0 … −60`) lies entirely above the lid and should vanish
      !! against the TOP, leaving the live layers at their open-ocean
      !! geopotential depths.  The analytic target is printed below.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp), parameter :: H_SURF_FINE = 20.0_wp
      real(wp), parameter :: DZ_COARSE = (B_FLAT - 3.0_wp*H_SURF_FINE)/real(NZ - 3, wp)
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2), h_bed(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: k
      checks: block
         call make_slot(vc, grid, VCOORD_ZSTAR_FULL)
         vc%zstar_h_surf_target = H_SURF_FINE
         h_bed = B_FLAT                      ! the TRUE bed, as production passes it
         total_h = B_FLAT - Z_TOP_LID        ! the live column under the lid
         eta = 0.0_wp
         call vc%build_zref_full(h_bed(1:grid%nx_total, 1:grid%ny_total))
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         ! TODAY (wrong): three fillers on the bed, one partial layer,
         ! three coarse layers, then the fine band hard against the lid.
         e_want(0) = -B_FLAT
         e_want(1) = -B_FLAT
         e_want(2) = -B_FLAT
         e_want(3) = -B_FLAT
         e_want(4) = -(B_FLAT - (Z_TOP_LID - 3.0_wp*DZ_COARSE - 3.0_wp*H_SURF_FINE))
         do k = 5, 7
            e_want(k) = e_want(4) + DZ_COARSE*real(k - 4, wp)
         end do
         do k = 8, NZ
            e_want(k) = e_want(7) + H_SURF_FINE*real(k - 7, wp)
         end do
         call check_depths(error, "ZSTAR_FULL under a lid (documented defect)", &
                           e_got, e_want, TOL_FILLER)
         if (allocated(error)) exit checks
         ! The headline: the FINE band sits at the bottom of the ice, not
         ! at z = 0.  P6.2 must move it above the lid and vanish it there.
         call check(error, abs(vc%target_h(2, 2, NZ) - H_SURF_FINE) < TOL_FILLER, &
                    "ZSTAR_FULL under a lid (documented defect): the fine band "// &
                    "should still be hanging from the lid before P6.2")
         if (allocated(error)) exit checks
         call check(error, vc%target_h(2, 2, 1) <= H_VANISHED, &
                    "ZSTAR_FULL under a lid (documented defect): the vanished "// &
                    "layers should still be on the BED before P6.2")
         if (allocated(error)) exit checks
         call check_column_sum(error, "ZSTAR_FULL", vc%target_h(2, 2, :), B_FLAT - Z_TOP_LID)
      end block checks
      call vc%destroy()
   end subroutine test_zstar_full_lid

   ! ------------------------------------------------------------------
   ! Z_FIXED — the first family P6.2 makes quasi-z
   ! ------------------------------------------------------------------

   subroutine test_z_fixed_open(error)
      !! Cases 1 + 2.  `h_nominal = z_fixed_h_ref/nz = 100 m`.  In the
      !! OPEN ocean, "depth below the column top" and "depth below z = 0"
      !! coincide at `η = 0`, so the live interfaces land on exact
      !! multiples of 100 m below `z = 0` and the leftover layers vanish
      !! against the bed — a quasi-z coordinate with a free-surface
      !! anchor and bed-side vanishing.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: i, k
      checks: block
         call make_slot(vc, grid, VCOORD_Z_FIXED)
         eta = 0.0_wp
         total_h = B_FLAT
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         do k = 0, NZ
            e_want(k) = -real(NZ - k, wp)*H_NOMINAL
         end do
         call check_depths(error, "Z_FIXED flat bed", e_got, e_want, TOL_EXACT)
         if (allocated(error)) exit checks
         do i = 1, grid%nx_total
            total_h(i, :) = B_SLOPE(i)
         end do
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         do i = 1, grid%nx_total
            call column_depths(vc%target_h(i, 2, :), B_SLOPE(i), e_got)
            call z_fixed_depths_from_column_top(B_SLOPE(i), B_SLOPE(i), e_want)
            call check_depths(error, "Z_FIXED sloping bed", e_got, e_want, TOL_FILLER)
            if (allocated(error)) exit checks
            call check_column_sum(error, "Z_FIXED", vc%target_h(i, 2, :), B_SLOPE(i))
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_z_fixed_open

   subroutine test_z_fixed_lid(error)
      !! Case 3, FLIPPED by P6.2.  This assertion used to be a
      !! `documents_*` case pinning the broken placement: handed a 500 m
      !! column under a lid at `z = -500`, `VCOORD_Z_FIXED` hung five
      !! 100 m layers from the LID and put five inert fillers on the BED,
      !! a 500 m per-index error at `e(5)`.
      !!
      !! It now asserts the intended table.  The slot carries the
      !! column-top depth (`vc%z_top`), the nominal interface depths stay
      !! GEOPOTENTIAL, the five layers whose nominal range lies inside
      !! the ice vanish against the TOP and the five live ones keep their
      !! open-ocean depths:
      !!
      !!   interface     before P6.2    now (= the analytic table)
      !!     e(0)         -1000.0 m       -1000.0 m   (the bed)
      !!     e(1)         -1000.0 m        -900.0 m
      !!     e(2)         -1000.0 m        -800.0 m
      !!     e(3)         -1000.0 m        -700.0 m
      !!     e(4)         -1000.0 m        -600.0 m
      !!     e(5)         -1000.0 m        -500.0 m   (the ice base)
      !!     e(6..10)      -900 .. -500    -500.0 m   (the filler stack)
      !!
      !! The two stacks span the same interval only because the lid depth
      !! is a whole multiple of `h_nominal`; the LAYER each interface
      !! belongs to is reversed, which is what every `k`-indexed consumer
      !! sees, and reversing it is the point of the slice.
      !!
      !! Here the lid falls exactly on a nominal level, so `k = 5` is a
      !! FULL 100 m layer rather than a partial cut; it is still debited
      !! the five fillers' `h_min`, which is the `5e-4 m` that keeps
      !! `Sum target_h` at exactly the 500 m water column.  The partial
      !! cut, its minimum thickness and the sliver merge are gated in
      !! `test_ocean_vcoord_zfixed_cavity`.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ)
      integer :: k, n_top_fill
      checks: block
         call make_slot(vc, grid, VCOORD_Z_FIXED)
         vc%z_top = Z_TOP_LID
         total_h = B_FLAT - Z_TOP_LID
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         call z_fixed_depths_under_rigid_top(B_FLAT, Z_TOP_LID, e_want)
         call check_depths(error, "Z_FIXED under a lid", e_got, e_want, TOL_FILLER)
         if (allocated(error)) exit checks
         ! The inversion, stated as the thing that had to change: the
         ! fillers are now under the TOP and the live stack is anchored to
         ! z = 0, not to the lid.
         n_top_fill = NZ - nint((B_FLAT - Z_TOP_LID)/H_NOMINAL)
         do k = NZ - n_top_fill + 1, NZ
            call check(error, vc%target_h(2, 2, k) <= H_VANISHED, &
                       "Z_FIXED under a lid: the layers that outcrop into the "// &
                       "ice must be inert fillers")
            if (allocated(error)) exit checks
         end do
         do k = 1, NZ - n_top_fill
            call check(error, vc%target_h(2, 2, k) > H_VANISHED, &
                       "Z_FIXED under a lid: every layer below the ice base "// &
                       "must be LIVE")
            if (allocated(error)) exit checks
         end do
         call check(error, abs(e_got(NZ - n_top_fill) + Z_TOP_LID) < TOL_FILLER, &
                    "Z_FIXED under a lid: e(5) must sit at the ICE BASE (-500 m), "// &
                    "not at the bed")
         if (allocated(error)) exit checks
         call check(error, abs(e_got(1) + (B_FLAT - H_NOMINAL)) < TOL_EXACT, &
                    "Z_FIXED under a lid: e(1) must keep its open-ocean "// &
                    "geopotential depth (-900 m)")
         if (allocated(error)) exit checks
         call check_column_sum(error, "Z_FIXED", vc%target_h(2, 2, :), B_FLAT - Z_TOP_LID)
      end block checks
      call vc%destroy()
   end subroutine test_z_fixed_lid

   subroutine test_z_fixed_z_top_zero(error)
      !! The other half of the P6.2 contract: with NO rigid top the
      !! family's pre-cavity answers must not move.  Same 500 m column,
      !! `z_top` left at its init-time zero — the target must reproduce
      !! the old walk BIT-for-bit (asserted with `==`, against a verbatim
      !! transcription held in this file), and the old analytic placement
      !! with it: the stack hangs from the column top at `z = -500` and
      !! the fillers sit on the BED.
      !!
      !! That is not a contradiction with the case above.  `z_top = 0` IS
      !! the open ocean, where "the column top" and "z = 0" are the same
      !! place, and a stack hung from the column top is then a
      !! geopotential stack.  The lid case is the only one that can tell
      !! the two readings apart, which is why it is the one that flipped.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: e_got(0:NZ), e_want(0:NZ), want_h(NZ)
      integer :: k
      checks: block
         call make_slot(vc, grid, VCOORD_Z_FIXED)
         total_h = B_FLAT - Z_TOP_LID
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call z_fixed_reference_pre_cavity(B_FLAT - Z_TOP_LID, want_h)
         do k = 1, NZ
            call check(error, vc%target_h(2, 2, k) == want_h(k), &
                       "Z_FIXED at z_top = 0 must be BIT-identical to the "// &
                       "pre-cavity target")
            if (allocated(error)) exit checks
         end do
         call column_depths(vc%target_h(2, 2, :), B_FLAT, e_got)
         call z_fixed_depths_from_column_top(B_FLAT, B_FLAT - Z_TOP_LID, e_want)
         call check_depths(error, "Z_FIXED at z_top = 0", e_got, e_want, TOL_FILLER)
         if (allocated(error)) exit checks
         call check_column_sum(error, "Z_FIXED", vc%target_h(2, 2, :), B_FLAT - Z_TOP_LID)
      end block checks
      call vc%destroy()
   end subroutine test_z_fixed_z_top_zero

   ! ------------------------------------------------------------------
   ! The inert-filler contract
   ! ------------------------------------------------------------------

   subroutine test_filler_contract(error)
      !! On the GEOMETRIC vanishing families the filler thickness must
      !! stay at or below `H_VANISHED = 1.5e-4 m`, because every
      !! downstream vanish test is a STRICT `> H_VANISHED`: a filler on
      !! the marker reads as vanished, one above it becomes dynamically
      !! live (EOS / PGF / remap-drain / vdiff) while the coordinate still
      !! treats it as throwaway.  Asserted on the shallowest column of the
      !! sloping bed for both families that produce fillers.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2), h_bed(NX + 2, NY + 2)
      real(wp) :: h_shallow
      checks: block
         h_shallow = B_SLOPE(NX + 2)
         call make_slot(vc, grid, VCOORD_Z_FIXED)
         total_h = h_shallow
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call check(error, vc%target_h(2, 2, 1) <= H_VANISHED, &
                    "Z_FIXED: the bed filler must be <= H_VANISHED")
         if (allocated(error)) exit checks
         call check(error, vc%target_h(2, 2, 1) > 0.0_wp, &
                    "Z_FIXED: the bed filler must be strictly positive")
         if (allocated(error)) exit checks
         call vc%destroy()

         call make_slot(vc, grid, VCOORD_ZSTAR_FULL)
         h_bed = h_shallow
         call vc%build_zref_full(h_bed(1:grid%nx_total, 1:grid%ny_total))
         call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                                  eta(1:grid%nx_total, 1:grid%ny_total))
         call check(error, vc%target_h(2, 2, 1) <= H_VANISHED, &
                    "ZSTAR_FULL: the bed filler must be <= H_VANISHED")
         if (allocated(error)) exit checks
         call check(error, vc%target_h(2, 2, 1) > 0.0_wp, &
                    "ZSTAR_FULL: the bed filler must be strictly positive")
      end block checks
      call vc%destroy()
   end subroutine test_filler_contract

end module test_ocean_vcoord_interface_depths
