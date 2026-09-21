!! Property-style INVARIANT suite for the vertical-coordinate target builder
!! and for the ALE remap column, run over an ADVERSARIAL set of columns.
!!
!! **How this differs from its sibling.**
!! `test_ocean_vcoord_interface_depths` asserts the absolute geopotential
!! DEPTH of every target interface against a hand-derived analytic table,
!! for three hand-picked geometries.  It answers *"is the stack in the right
!! place?"*.  This suite answers a different question — *"does the builder
!! ever produce a column that is not a column?"* — by sweeping properties
!! that must hold for EVERY family on EVERY column, including the ones no
!! shipped namelist has ever handed it:
!!
!!   * `H` from `2·H_VANISHED` (3e-4 m) to 5000 m, plus a land column;
!!   * `η ∈ {0, ±1e-10, ±0.01·H, ±0.3·H}` — the `±1e-10` pair is the
!!     FLIP-FLOP probe (P5 below);
!!   * `nz ∈ {1, 2, 15, 75}`, i.e. far more layers than a shallow column can
!!     hold;
!!   * a displaced column top `z_top` (a rigid lid / ice base).
!!
!! **The column properties.**
!!
!!   P1  finite          — no NaN/Inf anywhere in `target_h`.
!!   P2  conservation    — `Σ_k target_h = H + η` to round-off.  `EULERIAN_Z`
!!                         is the documented exception: it drops `η` by
!!                         construction, so its sum is `H`.
!!   P3  positivity      — every thickness is STRICTLY positive.
!!   P4  monotonicity    — the interface stack `e(0) = −b`,
!!                         `e(K) = e(K−1) + target_h(K)` is non-decreasing.
!!                         P3 implies it, but it is asserted separately so a
!!                         sign failure is reported as a stack inversion
!!                         rather than as a negative number.
!!   P5  no flip-flop    — the LIVE/VANISHED pattern (`target_h > H_VANISHED`)
!!                         is identical for `η` and `η ± 1e-10`.  A coordinate
!!                         whose pattern flickers under a 0.1 nm change in sea
!!                         level regrids a layer in and out of existence every
!!                         step, and the ALE drain then mints and destroys
!!                         tracer at the flicker rate.  P2 and P3 cannot see
!!                         that: a flickering layer is positive at every
!!                         sample and conserves at every sample.
!!
!! **The thin-column envelope, which P2/P3 are SCOPED to.**  The two families
!! that lay inert fillers (`VCOORD_ZSTAR_FULL`, `VCOORD_Z_FIXED`) give every
!! layer a floor of `zstar_h_min`, so a column with `H + η < nz·zstar_h_min`
!! cannot hold its own stack: the floor wins and the target column MINTS
!! thickness.  Measured here, and pinned in
!! `documents_filler_floor_beats_a_thin_column`:
!!
!!   * `ZSTAR_FULL`, `H = 3e-4 m`, `nz = 15` — 14 fillers of `1e-4` plus one
!!     layer of **exactly 0.0**, summing to `1.400000e-03` against a 3e-4 m
!!     column: conservation out by **4.67x**, and a zero-thickness layer is
!!     below the family's own `zstar_h_min` promise.
!!   * `Z_FIXED`, `H = 3e-4 m`, `nz = 75` — 75 fillers of `1e-4`, summing to
!!     `7.500000e-03`: conservation out by **25x**.
!!
!! Nothing downstream currently hands the builder such a column (the wet
!! mask drops a cell long before `H` reaches `nz·zstar_h_min`), which is why
!! no shipped case sees it — and exactly why it is worth a gate: the
!! condition is `nz·zstar_h_min`, so it moves with the layer count, not with
!! the physics.
!!
!! **Remap-column properties**, for every `remap_method`:
!!
!!   R1  identity        — `dz_new == dz_old` must reproduce `q_old`.  It
!!                         does NOT do so bit-for-bit (pinned below).
!!   R2  conservation    — `Σ q·dz` is preserved to round-off, including when
!!                         vanished layers sit at the top, at the bed, or in
!!                         the interior.  HOLDS for all five methods.
!!   R3  no new extrema  — a limited scheme must not overshoot the old
!!                         column's range, on a step profile.  HOLDS for all
!!                         five methods.
!!   R4  linear exactness— a profile linear in `z` remaps exactly.  This is
!!                         the property that separates the methods, and the
!!                         answer depends on whether the OLD grid is uniform:
!!
!!         old grid      pcm        plm        ppm        ppm_h4     pqm
!!         uniform       2.2e-02    5.7e-14    5.0e-14    5.0e-14    5.0e-14
!!         stretched     4.3e-02    2.3e-03    6.3e-04    5.7e-14    1.1e-13
!!
!!                         (worst absolute error over cells 3…nz−2 on a
!!                         `q = 34.2 + 0.013·z` profile, 20 layers, gfortran
!!                         15.1 Release.)  `pcm` is first order by definition
!!                         and is not gated on R4.  `plm` and `ppm` are
!!                         exact on a UNIFORM column and lose it on a
!!                         stretched one — their slope / edge stencils are
!!                         the uniform-grid forms — while `ppm_h4` and `pqm`
!!                         carry the non-uniform h4 stencil (White & Adcroft
!!                         2008) and stay exact.  Every ALE step in this
!!                         model hands the remap a STRETCHED column (a sigma
!!                         column over a slope, a z* column with fillers), so
!!                         the shipped default `remap_method = ppm` is the
!!                         6.3e-04 cell of that table, not the 5.0e-14 one.
!!
!! Anything this suite pins as `documents_*` is a defect it has MEASURED, not
!! a bar it has chosen.  Nothing here was tuned to pass.
module test_ocean_vcoord_invariants
   use rdb_constants, only: wp, H_VANISHED, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED, &
                            REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, &
                            REMAP_PQM
   use rdb_grid, only: hgrid_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z
   use rdb_remap_column, only: remap_column
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_invariants_tests

   ! ---------------------------------------------------------------------
   ! The adversarial column set
   ! ---------------------------------------------------------------------
   integer, parameter :: N_FAM = 7
      !! Every GEOMETRIC family `ocean_vcoord_compute_target_h` dispatches on.
      !! `VCOORD_LAGRANGIAN` is excluded: it returns without writing
      !! `target_h` at all (gated in the interface-depth suite).  The two
      !! density-space families come in through `compute_target_h_rho` and
      !! need per-layer T/S + the EOS, so they are gated by
      !! `test_ocean_vcoord_rho` / `test_ocean_vcoord_hycom`.
   integer, parameter :: FAM_CODE(N_FAM) = [VCOORD_EULERIAN_Z, VCOORD_SIGMA, &
                                            VCOORD_ZSTAR, VCOORD_ZSTAR_SIGMA, &
                                            VCOORD_ZSTAR_FULL, VCOORD_Z_FIXED, &
                                            VCOORD_ZSIGMA]
   character(len=12), parameter :: FAM_NAME(N_FAM) = [character(len=12) :: &
                                                      "EULERIAN_Z", "SIGMA", "ZSTAR", "ZSTAR_SIGMA", &
                                                      "ZSTAR_FULL", "Z_FIXED", "ZSIGMA"]
   logical, parameter :: FAM_FILLS(N_FAM) = [.false., .false., .false., .false., &
                                             .true., .true., .false.]
      !! Does the family lay `zstar_h_min` fillers?  Only these two have a
      !! thin-column envelope (see the header), so only these two are scoped
      !! out of P2/P3 below `nz·zstar_h_min`.

   integer, parameter :: N_BED = 5
   real(wp), parameter :: BED(N_BED) = [3.0e-4_wp, 1.0_wp, 50.0_wp, &
                                        500.0_wp, 5000.0_wp]
      !! `2·H_VANISHED` (the thinnest column the solver still calls wet) up
      !! to an abyssal 5000 m.

   integer, parameter :: N_NZ = 4
   integer, parameter :: NZ_SET(N_NZ) = [1, 2, 15, 75]
      !! Includes `nz = 75` against a 3e-4 m bed: 250 000 times more layers
      !! than the column can hold at any sane thickness.

   real(wp), parameter :: ETA_TINY = 1.0e-10_wp
      !! The FLIP-FLOP probe.  0.1 nm of sea level: nothing physical, and far
      !! below every threshold any family documents.

   real(wp), parameter :: TOL_REL = 1.0e-12_wp
      !! Relative round-off bar for the conservation sums.
   real(wp), parameter :: H_MIN_SLOT = 1.0e-4_wp
      !! `zstar_h_min` as this suite sets it — at or below `H_VANISHED`, the
      !! rule `VCOORD_ZSTAR_FULL` and `VCOORD_Z_FIXED` fail loud on.

   integer, parameter :: N_METHOD = 5
   integer, parameter :: M_CODE(N_METHOD) = [REMAP_PCM, REMAP_PLM, REMAP_PPM, &
                                             REMAP_PPM_H4, REMAP_PQM]
   character(len=6), parameter :: M_NAME(N_METHOD) = [character(len=6) :: &
                                                      "pcm", "plm", "ppm", "ppm_h4", "pqm"]

contains

   subroutine collect_ocean_vcoord_invariants_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("column_invariants_over_every_family", &
                               test_column_invariants), &
                  new_unittest("eta_perturbation_does_not_flip_flop", &
                               test_no_flip_flop), &
                  new_unittest("displaced_top_column_invariants", &
                               test_displaced_top), &
                  new_unittest("documents_uniform_stack_on_the_vanish_marker", &
                               test_vanish_marker), &
                  new_unittest("documents_filler_floor_beats_a_thin_column", &
                               test_thin_column), &
                  new_unittest("land_column_is_finite_and_non_negative", test_land), &
                  new_unittest("remap_conserves_with_vanished_layers", &
                               test_remap_conservation), &
                  new_unittest("remap_adds_no_new_extrema", test_remap_extrema), &
                  new_unittest("remap_linear_exact_on_a_uniform_column", &
                               test_remap_linear_uniform), &
                  new_unittest("remap_h4_stencils_linear_exact_when_stretched", &
                               test_remap_linear_stretched), &
                  new_unittest("documents_remap_identity_is_not_bit_exact", &
                               test_remap_identity), &
                  new_unittest("documents_plm_ppm_lose_linearity_when_stretched", &
                               test_remap_linear_stretched_defect) &
                  ]
   end subroutine collect_ocean_vcoord_invariants_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   pure function is_finite(x) result(ok)
      !! `.true.` iff `x` is neither NaN nor Inf, written without
      !! `ieee_arithmetic` so the helper stays usable from a `pure` context
      !! on every toolchain in the matrix.
      real(wp), intent(in) :: x
      logical :: ok
      ok = (x == x) .and. (abs(x) <= huge(1.0_wp))
   end function is_finite

   pure function sum_tol(column_total, nz) result(tol)
      !! Conservation bar for one column: relative round-off on the column
      !! total, but never tighter than the absolute round-off of the filler
      !! budget a vanishing family lays down (`nz·zstar_h_min`).
      real(wp), intent(in) :: column_total
      integer, intent(in) :: nz
      real(wp) :: tol
      tol = max(TOL_REL*max(column_total, 1.0_wp), &
                1.0e-9_wp*real(nz, wp)*H_MIN_SLOT)
   end function sum_tol

   pure function on_vanish_marker(column_total, nz) result(onit)
      !! `.true.` when a UNIFORM stack of `nz` layers lands every interface
      !! exactly on `H_VANISHED` — `column_total = nz·H_VANISHED`.  This is
      !! the ONE documented threshold in the P5 sweep: the live test is a
      !! strict `> H_VANISHED`, so on that column an arbitrarily small `η`
      !! decides the answer, for every family whose stack is uniform.  It
      !! is not a defect of a family, it is what "strictly greater" means at
      !! the marker, and it is asserted directly by
      !! `documents_uniform_stack_on_the_vanish_marker` rather than skipped
      !! silently here.
      real(wp), intent(in) :: column_total
      integer, intent(in) :: nz
      logical :: onit
      onit = abs(column_total - real(nz, wp)*H_VANISHED) &
             <= 1.0e-9_wp*real(nz, wp)*H_VANISHED
   end function on_vanish_marker

   pure function thin_for_family(ifam, column_total, nz) result(thin)
      !! `.true.` when this column falls inside the documented thin-column
      !! envelope of a FILLER family — `H + η` below the whole-column filler
      !! budget `nz·zstar_h_min`, with a factor of 2 of headroom because the
      !! reserve branch needs room for a partial cell on top of the fillers.
      !! Those rows are asserted by `documents_filler_floor_beats_a_thin_
      !! column` instead, with the measured numbers.
      integer, intent(in) :: ifam, nz
      real(wp), intent(in) :: column_total
      logical :: thin
      thin = FAM_FILLS(ifam) .and. &
             (column_total < 2.0_wp*real(nz, wp)*H_MIN_SLOT)
   end function thin_for_family

   subroutine build_column(vc, grid, family, nz, bed_depth, eta_val, z_top_val)
      !! Fresh slot on a 1x1 interior grid with one ghost ring, configured
      !! the way production configures it, then one `compute_target_h` call
      !! for the (bed, eta) column.  Not `pure`: `init` allocates.
      type(ocean_vcoord_t), intent(out) :: vc
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: family
      integer, intent(in) :: nz
      real(wp), intent(in) :: bed_depth
      real(wp), intent(in) :: eta_val
      real(wp), intent(in) :: z_top_val
      real(wp) :: total_h(3, 3), eta(3, 3), h_bed(3, 3)
      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=nz)
      vc%coord_type = family
      vc%zstar_h_min = H_MIN_SLOT
      vc%zstar_h_surf_target = 100.0_wp
      vc%z_fixed_h_ref = 1000.0_wp
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 100.0_wp
      vc%z_top = z_top_val
      total_h = bed_depth
      eta = eta_val
      h_bed = bed_depth
      if (family == VCOORD_ZSTAR_FULL) then
         call vc%build_zref_full(h_bed(1:grid%nx_total, 1:grid%ny_total))
      end if
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))
   end subroutine build_column

   pure function expected_total(family, bed_depth, eta_val) result(total)
      !! What `Σ_k target_h` must equal.  `EULERIAN_Z` is `H·dsig` — it drops
      !! `η` by construction (that is the family's definition, not a defect),
      !! so its column total is the bed depth alone.
      integer, intent(in) :: family
      real(wp), intent(in) :: bed_depth, eta_val
      real(wp) :: total
      if (family == VCOORD_EULERIAN_Z) then
         total = bed_depth
      else
         total = bed_depth + eta_val
      end if
   end function expected_total

   ! ------------------------------------------------------------------
   ! P1..P4 over the whole sweep
   ! ------------------------------------------------------------------

   subroutine test_column_invariants(error)
      !! P1 finite, P2 `Σ target_h = H + η`, P3 strict positivity, P4 a
      !! monotone interface stack — asserted for every family on every
      !! (bed, η, nz) column of the adversarial set, outside the documented
      !! thin-column envelope.
      !!
      !! 7 families x 5 beds x 7 sea levels x 4 layer counts = 980 columns.
      !! The first miss stops the sweep and names the family, the column and
      !! the property, so a failure reads as a coordinate defect rather than
      !! as "a vcoord test failed".
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      character(len=220) :: msg
      real(wp) :: eta_set(7)
      real(wp) :: bed_depth, eta_val, want, got, e_prev, e_here
      integer :: ifam, ibed, ieta, inz, nz, k
      sweep: block
         do ifam = 1, N_FAM
            do ibed = 1, N_BED
               bed_depth = BED(ibed)
               eta_set = [0.0_wp, ETA_TINY, -ETA_TINY, &
                          0.01_wp*bed_depth, -0.01_wp*bed_depth, &
                          0.30_wp*bed_depth, -0.30_wp*bed_depth]
               do ieta = 1, size(eta_set)
                  eta_val = eta_set(ieta)
                  do inz = 1, N_NZ
                     nz = NZ_SET(inz)
                     want = expected_total(FAM_CODE(ifam), bed_depth, eta_val)
                     if (thin_for_family(ifam, want, nz)) cycle
                     call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                       bed_depth, eta_val, 0.0_wp)
                     got = 0.0_wp
                     e_prev = -bed_depth
                     do k = 1, nz
                        got = got + vc%target_h(2, 2, k)
                        ! --- P1 finite
                        if (.not. is_finite(vc%target_h(2, 2, k))) then
                           write (msg, '(a,i0,a,es10.3,a,es10.3,a,i0)') &
                              trim(FAM_NAME(ifam))//" P1 finite: target_h(k=", k, &
                              ") is not finite; H = ", bed_depth, ", eta = ", &
                              eta_val, ", nz = ", nz
                           call check(error, .false., trim(msg))
                           exit sweep
                        end if
                        ! --- P3 strict positivity
                        if (.not. (vc%target_h(2, 2, k) > 0.0_wp)) then
                           write (msg, '(a,i0,a,es13.6,a,es10.3,a,es10.3,a,i0)') &
                              trim(FAM_NAME(ifam))//" P3 positivity: target_h(k=", k, &
                              ") = ", vc%target_h(2, 2, k), " <= 0; H = ", bed_depth, &
                              ", eta = ", eta_val, ", nz = ", nz
                           call check(error, .false., trim(msg))
                           exit sweep
                        end if
                        ! --- P4 monotone interface stack
                        e_here = e_prev + vc%target_h(2, 2, k)
                        if (e_here < e_prev) then
                           write (msg, '(a,i0,a,es13.6,a,es13.6)') &
                              trim(FAM_NAME(ifam))//" P4 monotone: e(", k, ") = ", &
                              e_here, " is below e(k-1) = ", e_prev
                           call check(error, .false., trim(msg))
                           exit sweep
                        end if
                        e_prev = e_here
                     end do
                     ! --- P2 conservation
                     if (abs(got - want) > sum_tol(want, nz)) then
                        write (msg, '(a,es13.6,a,es13.6,a,es10.3,a,es10.3,a,i0)') &
                           trim(FAM_NAME(ifam))//" P2 conservation: sum(target_h) = ", &
                           got, " but H + eta = ", want, "; H = ", bed_depth, &
                           ", eta = ", eta_val, ", nz = ", nz
                        call check(error, .false., trim(msg))
                        exit sweep
                     end if
                     call vc%destroy()
                  end do
               end do
            end do
         end do
      end block sweep
      call vc%destroy()
   end subroutine test_column_invariants

   ! ------------------------------------------------------------------
   ! P5 -- the flip-flop probe
   ! ------------------------------------------------------------------

   subroutine test_no_flip_flop(error)
      !! P5.  Build each column three times — at `η`, at `η + 1e-10` and at
      !! `η − 1e-10` — and require the LIVE/VANISHED pattern
      !! (`target_h > H_VANISHED`) to be identical across the three.
      !!
      !! **Documented thresholds.**  A family MAY legitimately change its
      !! pattern when a nominal interface coincides with the free surface —
      !! `Z_FIXED` at `H + η = m·h_nominal` and `ZSTAR_FULL` at a table
      !! entry.  The sweep deliberately avoids those coincidences (the bed
      !! set is 3e-4 / 1 / 50 / 500 / 5000 m against `h_nominal = 1000/nz`,
      !! and the `±0.01 H` / `±0.3 H` offsets move every column off a round
      !! number), so any pattern change this test sees is a flicker at an
      !! arbitrary sea level and is a defect.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      character(len=220) :: msg
      logical :: live0(80), live_p(80), live_m(80)
      real(wp) :: eta_set(3)
      real(wp) :: bed_depth, eta_val
      integer :: ifam, ibed, ieta, inz, nz, k
      sweep: block
         do ifam = 1, N_FAM
            do ibed = 1, N_BED
               bed_depth = BED(ibed)
               eta_set = [0.0_wp, 0.01_wp*bed_depth, -0.30_wp*bed_depth]
               do ieta = 1, size(eta_set)
                  eta_val = eta_set(ieta)
                  do inz = 1, N_NZ
                     nz = NZ_SET(inz)
                     if (on_vanish_marker(bed_depth + eta_val, nz)) cycle
                     call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                       bed_depth, eta_val, 0.0_wp)
                     live0(1:nz) = vc%target_h(2, 2, 1:nz) > H_VANISHED
                     call vc%destroy()
                     call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                       bed_depth, eta_val + ETA_TINY, 0.0_wp)
                     live_p(1:nz) = vc%target_h(2, 2, 1:nz) > H_VANISHED
                     call vc%destroy()
                     call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                       bed_depth, eta_val - ETA_TINY, 0.0_wp)
                     live_m(1:nz) = vc%target_h(2, 2, 1:nz) > H_VANISHED
                     do k = 1, nz
                        if ((live0(k) .neqv. live_p(k)) .or. &
                            (live0(k) .neqv. live_m(k))) then
                           write (msg, '(a,i0,a,es10.3,a,es10.3,a,i0,a)') &
                              trim(FAM_NAME(ifam))//" P5 flip-flop: layer ", k, &
                              " changes live/vanished state under a 1e-10 m sea-level "// &
                              "perturbation; H = ", bed_depth, ", eta = ", eta_val, &
                              ", nz = ", nz, " (no family documents a threshold here)"
                           call check(error, .false., trim(msg))
                           exit sweep
                        end if
                     end do
                     call vc%destroy()
                  end do
               end do
            end do
         end do
      end block sweep
      call vc%destroy()
   end subroutine test_no_flip_flop

   ! ------------------------------------------------------------------
   ! The displaced column top
   ! ------------------------------------------------------------------

   subroutine test_displaced_top(error)
      !! The same invariants with the column top displaced to `z = −z_top`
      !! (a rigid lid / ice base).  `VCOORD_Z_FIXED` is the one family that
      !! READS `vcoord%z_top`; the rest are handed it and must be
      !! indifferent to it, which is itself a property worth pinning — a
      !! family that silently started consuming the slot would change every
      !! cavity answer without any test noticing.
      !!
      !! The live column here is `H − z_top`, so the conservation target is
      !! `H − z_top`: `compute_target_h` is handed the LIVE column, exactly
      !! as `rdb_ocean_remap` hands it one under a cavity.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      character(len=220) :: msg
      real(wp) :: z_top_set(3)
      real(wp) :: bed_depth, live_col, want, got, z_top_val
      integer :: ifam, ibed, itop, inz, nz, k
      sweep: block
         do ifam = 1, N_FAM
            do ibed = 2, N_BED          ! the 3e-4 m bed has no room for a lid
               bed_depth = BED(ibed)
               z_top_set = [0.1_wp*bed_depth, 0.5_wp*bed_depth, 0.9_wp*bed_depth]
               do itop = 1, size(z_top_set)
                  z_top_val = z_top_set(itop)
                  live_col = bed_depth - z_top_val
                  do inz = 1, N_NZ
                     nz = NZ_SET(inz)
                     want = live_col
                     if (thin_for_family(ifam, want, nz)) cycle
                     call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                       live_col, 0.0_wp, z_top_val)
                     got = 0.0_wp
                     do k = 1, nz
                        got = got + vc%target_h(2, 2, k)
                        if (.not. is_finite(vc%target_h(2, 2, k)) .or. &
                            .not. (vc%target_h(2, 2, k) > 0.0_wp)) then
                           write (msg, '(a,i0,a,es13.6,a,es10.3,a,es10.3,a,i0)') &
                              trim(FAM_NAME(ifam))//" displaced top: target_h(k=", k, &
                              ") = ", vc%target_h(2, 2, k), " is not a live thickness; "// &
                              "H = ", bed_depth, ", z_top = ", z_top_val, ", nz = ", nz
                           call check(error, .false., trim(msg))
                           exit sweep
                        end if
                     end do
                     if (abs(got - want) > sum_tol(want, nz)) then
                        write (msg, '(a,es13.6,a,es13.6,a,es10.3,a,i0)') &
                           trim(FAM_NAME(ifam))//" displaced top: sum(target_h) = ", &
                           got, " but the live column is ", want, "; z_top = ", &
                           z_top_val, ", nz = ", nz
                        call check(error, .false., trim(msg))
                        exit sweep
                     end if
                     call vc%destroy()
                  end do
               end do
            end do
         end do
      end block sweep
      call vc%destroy()
   end subroutine test_displaced_top

   subroutine test_vanish_marker(error)
      !! `documents_*` — the ONE flip-flop threshold the P5 sweep skips,
      !! asserted here directly so the list is machine-checked rather than
      !! written down in a comment.
      !!
      !! On a column of exactly `nz·H_VANISHED` a uniform-stack family puts
      !! every layer AT the vanish marker, and because every live test in
      !! the solver is a strict `> H_VANISHED`, an `η` of 1e-10 m decides
      !! whether the whole column is live or the whole column is throwaway.
      !! Measured for `VCOORD_SIGMA` at `H = 3.0e-4`, `nz = 2`
      !! (`h = 1.5e-4 = H_VANISHED` exactly): `η = 0` and `η = −1e-10` read
      !! VANISHED, `η = +1e-10` reads LIVE.
      !!
      !! This is a property of the marker, not of a family — but it is also
      !! the reason `zstar_h_min` must stay strictly BELOW `H_VANISHED`
      !! (which `validate_config` enforces for `ZSTAR_FULL` and `Z_FIXED`):
      !! a filler laid AT the marker would sit on this knife edge on every
      !! column, every step.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp), parameter :: H_MARK = 2.0_wp*H_VANISHED
      logical :: live_zero, live_plus
      checks: block
         call build_column(vc, grid, VCOORD_SIGMA, 2, H_MARK, 0.0_wp, 0.0_wp)
         live_zero = vc%target_h(2, 2, 1) > H_VANISHED
         call vc%destroy()
         call build_column(vc, grid, VCOORD_SIGMA, 2, H_MARK, ETA_TINY, 0.0_wp)
         live_plus = vc%target_h(2, 2, 1) > H_VANISHED
         call check(error, (.not. live_zero) .and. live_plus, &
                    "SIGMA on the vanish marker (H = 2*H_VANISHED, nz = 2): "// &
                    "the documented knife edge is gone. If the live test is no "// &
                    "longer a strict `> H_VANISHED`, delete this case and drop "// &
                    "`on_vanish_marker` from the P5 sweep.")
         if (allocated(error)) exit checks
      end block checks
      call vc%destroy()
   end subroutine test_vanish_marker

   ! ------------------------------------------------------------------
   ! The thin-column envelope -- pinned, not fixed
   ! ------------------------------------------------------------------

   subroutine test_thin_column(error)
      !! `documents_*` — this asserts BROKEN behaviour on purpose.
      !!
      !! Below `nz·zstar_h_min` the filler floor of `VCOORD_ZSTAR_FULL` and
      !! `VCOORD_Z_FIXED` is thicker than the column it is filling, and the
      !! target grid stops being a partition of the water column.  Measured
      !! on gfortran 15.1 Release, `zstar_h_min = 1e-4`, `H = 3e-4 m`,
      !! `η = 0`:
      !!
      !!   family      nz    Σ target_h      should be   ratio   worst layer
      !!   ZSTAR_FULL  15    1.400000e-03    3.0e-04     4.67x   0.0 (k=15)
      !!   Z_FIXED     75    7.500000e-03    3.0e-04     25.0x   1.0e-04
      !!
      !! The `ZSTAR_FULL` row is the worse of the two: it mints thickness
      !! AND emits a layer of EXACTLY ZERO, which is below the family's own
      !! `zstar_h_min` promise and is what every `1/h` guard in the solver
      !! is armoured against seeing.
      !!
      !! The intended behaviour is that the builder either (a) conserves and
      !! lets the layers go below `zstar_h_min`, or (b) refuses the column
      !! fail-loud at configure the way `zstar_h_min > H_VANISHED` is
      !! refused.  Whichever lands, it must DELETE this case rather than
      !! update it.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      character(len=220) :: msg
      real(wp), parameter :: H_THIN = 3.0e-4_wp
      real(wp) :: got
      checks: block
         ! --- ZSTAR_FULL, nz = 15: 1.4e-3 against a 3e-4 column, one zero layer
         call build_column(vc, grid, VCOORD_ZSTAR_FULL, 15, H_THIN, 0.0_wp, 0.0_wp)
         got = sum(vc%target_h(2, 2, :))
         if (abs(got - 1.4e-3_wp) > 1.0e-9_wp) then
            write (msg, '(a,es13.6,a)') &
               "ZSTAR_FULL thin column: sum(target_h) = ", got, &
               " -- was 1.400000e-03. The thin-column envelope moved; "// &
               "re-measure and either delete this case or restate it."
            call check(error, .false., trim(msg))
            exit checks
         end if
         call check(error, minval(vc%target_h(2, 2, :)) == 0.0_wp, &
                    "ZSTAR_FULL thin column: the exactly-zero layer is gone -- "// &
                    "if the builder now floors every layer, delete this case.")
         if (allocated(error)) exit checks
         call vc%destroy()
         ! --- Z_FIXED, nz = 75: the whole filler budget, 25x the column
         call build_column(vc, grid, VCOORD_Z_FIXED, 75, H_THIN, 0.0_wp, 0.0_wp)
         got = sum(vc%target_h(2, 2, :))
         if (abs(got - 75.0_wp*H_MIN_SLOT) > 1.0e-9_wp) then
            write (msg, '(a,es13.6,a)') &
               "Z_FIXED thin column: sum(target_h) = ", got, &
               " -- was 7.500000e-03 (= nz*zstar_h_min). The thin-column "// &
               "envelope moved; re-measure or delete this case."
            call check(error, .false., trim(msg))
            exit checks
         end if
      end block checks
      call vc%destroy()
   end subroutine test_thin_column

   ! ------------------------------------------------------------------
   ! Land
   ! ------------------------------------------------------------------

   subroutine test_land(error)
      !! A land column (`H = 0`, `η = 0`).  Nothing may be NaN and nothing
      !! may be negative: a land column reaches the EOS, the PGF and the
      !! remap through the same arrays as a wet one, and a NaN seeded here
      !! propagates into the interior through the halo on the first
      !! exchange.  That part HOLDS for all seven families.
      !!
      !! The column total does NOT come out at zero for the two filler
      !! families: `ZSTAR_FULL` and `Z_FIXED` emit a full stack of
      !! `zstar_h_min`, i.e. `nz·1e-4` m of target column on a dry cell
      !! (measured: 1.5e-03 m at `nz = 15`).  That is the `H = 0` corner of
      !! the same thin-column envelope pinned above, and it is bounded here
      !! rather than asserted away: the target may not exceed the filler
      !! budget, which is what keeps a dry cell from reaching the remap with
      !! a metre-scale column.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      character(len=220) :: msg
      real(wp) :: got
      integer :: ifam, inz, nz, k
      sweep: block
         do ifam = 1, N_FAM
            do inz = 1, N_NZ
               nz = NZ_SET(inz)
               call build_column(vc, grid, FAM_CODE(ifam), nz, &
                                 0.0_wp, 0.0_wp, 0.0_wp)
               got = 0.0_wp
               do k = 1, nz
                  got = got + vc%target_h(2, 2, k)
                  if (.not. is_finite(vc%target_h(2, 2, k)) .or. &
                      vc%target_h(2, 2, k) < 0.0_wp) then
                     write (msg, '(a,i0,a,es13.6,a,i0)') &
                        trim(FAM_NAME(ifam))//" land column: target_h(k=", k, &
                        ") = ", vc%target_h(2, 2, k), " on H = 0, nz = ", nz
                     call check(error, .false., trim(msg))
                     exit sweep
                  end if
               end do
               if (got > real(nz, wp)*H_MIN_SLOT*(1.0_wp + TOL_REL)) then
                  write (msg, '(a,es13.6,a,es13.6,a,i0)') &
                     trim(FAM_NAME(ifam))//" land column: sum(target_h) = ", got, &
                     " exceeds the filler budget ", real(nz, wp)*H_MIN_SLOT, &
                     " on a dry column, nz = ", nz
                  call check(error, .false., trim(msg))
                  exit sweep
               end if
               call vc%destroy()
            end do
         end do
      end block sweep
      call vc%destroy()
   end subroutine test_land

   ! ------------------------------------------------------------------
   ! Remap-column properties
   ! ------------------------------------------------------------------

   pure subroutine rebalance(nz, dz, total)
      !! Scale `dz` so it spans the same interval as the old column: the
      !! remap is only defined between two partitions of ONE interval.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: dz(nz)
      real(wp), intent(in) :: total
      dz = dz*(total/sum(dz))
   end subroutine rebalance

   pure subroutine vanish_layers(nz, dz, k_lo, k_hi, total)
      !! Collapse `k_lo…k_hi` onto `H_VANISHED` and give the freed thickness
      !! back to the remaining layers, so the column still spans `total`.
      integer, intent(in) :: nz, k_lo, k_hi
      real(wp), intent(inout) :: dz(nz)
      real(wp), intent(in) :: total
      real(wp) :: freed, rest
      integer :: k
      freed = 0.0_wp
      do k = k_lo, k_hi
         freed = freed + dz(k) - H_VANISHED
         dz(k) = H_VANISHED
      end do
      rest = 0.0_wp
      do k = 1, nz
         if (k < k_lo .or. k > k_hi) rest = rest + dz(k)
      end do
      do k = 1, nz
         if (k < k_lo .or. k > k_hi) dz(k) = dz(k)*(rest + freed)/rest
      end do
      dz = dz*(total/sum(dz))
   end subroutine vanish_layers

   pure subroutine linear_column(nz, uniform_old, dz_old, dz_new, q_old, z_new, &
                                 slope, base)
      !! Build the shared linear-profile column: an old grid (uniform or
      !! stretched), a DIFFERENT stretched new grid spanning the same
      !! interval, and `q_old` set to the exact CELL AVERAGE of
      !! `base + slope·z` (which, for a linear profile, is the value at the
      !! cell centre).
      integer, intent(in) :: nz
      logical, intent(in) :: uniform_old
      real(wp), intent(out) :: dz_old(nz), dz_new(nz), q_old(nz)
      real(wp), intent(out) :: z_new(0:nz)
      real(wp), intent(in) :: slope, base
      real(wp) :: z_old(0:nz), total
      integer :: k
      do k = 1, nz
         if (uniform_old) then
            dz_old(k) = 20.0_wp
         else
            dz_old(k) = 20.0_wp + 4.0_wp*sin(0.8_wp*real(k, wp))
         end if
         dz_new(k) = 20.0_wp + 5.0_wp*cos(0.55_wp*real(k, wp))
      end do
      total = sum(dz_old)
      dz_new = dz_new*(total/sum(dz_new))
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
         q_old(k) = base + slope*0.5_wp*(z_old(k - 1) + z_old(k))
      end do
   end subroutine linear_column

   subroutine test_remap_identity(error)
      !! `documents_*` — this asserts a KNOWN SHORTFALL on purpose.
      !!
      !! R1.  `dz_new == dz_old` OUGHT to reproduce `q_old` bit-for-bit: an
      !! identity remap that is merely accurate to 1e-15 is a per-step
      !! diffusion of the whole tracer field, and the ALE remap fires every
      !! thermo step for the life of a run.
      !!
      !! It does not.  Every method rebuilds the interface stack by
      !! cumulative sum and recomputes `∫q dz / dz_new` even where the two
      !! stacks coincide, so the answer comes back one to three ulp off.
      !! Measured (20 layers, stretched column, gfortran 15.1 Release):
      !!
      !!   method   cells not bit-exact (of 20)   worst |Δq|
      !!   pcm      16                            2.84e-14
      !!   plm      16                            2.84e-14
      !!   ppm      16                            2.84e-14
      !!   ppm_h4   16                            2.84e-14
      !!   pqm      18                            8.53e-14
      !!
      !! against `q ~ 34`, i.e. 8e-16 to 2.5e-15 relative — 4 to 12 ulp of
      !! `real64`.  The fix is a short circuit on `dz_new == dz_old`, which
      !! this branch does not carry.  What IS gated here is that the error
      !! stays at that scale: the bound below is 20 ulp relative, so a
      !! reconstruction change that turned the identity into a real
      !! diffusion would fail it.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_MAX = 60
      integer, parameter :: NZR(4) = [1, 2, 15, NZ_MAX]
      character(len=220) :: msg
      real(wp) :: dz(NZ_MAX), q_old(NZ_MAX), q_new(NZ_MAX)
      real(wp) :: worst, bar
      integer :: im, inz, nz, k
      sweep: block
         do im = 1, N_METHOD
            do inz = 1, size(NZR)
               nz = NZR(inz)
               do k = 1, nz
                  dz(k) = 5.0_wp + 0.37_wp*real(k, wp)
                  q_old(k) = 34.0_wp + 0.11_wp*real(k, wp) &
                             + 0.4_wp*sin(0.7_wp*real(k, wp))
               end do
               call remap_column(M_CODE(im), nz, dz(1:nz), dz(1:nz), &
                                 q_old(1:nz), q_new(1:nz))
               worst = maxval(abs(q_new(1:nz) - q_old(1:nz)))
               bar = 20.0_wp*epsilon(1.0_wp)*maxval(abs(q_old(1:nz)))
               if (.not. is_finite(worst) .or. worst > bar) then
                  write (msg, '(a,es12.5,a,es12.5,a,i0)') &
                     "remap "//trim(M_NAME(im))// &
                     " R1 identity: worst |dq| = ", worst, " over the 20-ulp bar ", &
                     bar, " at nz = ", nz
                  call check(error, .false., trim(msg))
                  exit sweep
               end if
            end do
         end do
      end block sweep
   end subroutine test_remap_identity

   subroutine test_remap_conservation(error)
      !! R2.  `Σ q·dz` preserved to round-off, on four column shapes that
      !! between them put a vanished layer everywhere it can go:
      !!
      !!   1. a plain regrid (no vanished layers) — the control;
      !!   2. vanished layers at the BED (`k = 1…3`), the `Z_FIXED` /
      !!      `ZSTAR_FULL` shallow-column shape;
      !!   3. vanished layers at the TOP (`k = nz−2…nz`), the shape
      !!      `Z_FIXED` produces under an ice base;
      !!   4. a vanished layer in the INTERIOR, which no shipped geometric
      !!      coordinate builds today but which an isopycnal outcrop
      !!      (`VCOORD_RHO`, `VCOORD_HYCOM`) makes routine.
      !!
      !! The vanished thickness used is `H_VANISHED` itself, the marker the
      !! rest of the solver tests against.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_C = 15
      character(len=220) :: msg
      real(wp) :: dz_old(NZ_C), dz_new(NZ_C), q_old(NZ_C), q_new(NZ_C)
      real(wp) :: total_old, total_new, tol
      integer :: im, ishape, k
      sweep: block
         do ishape = 1, 4
            do k = 1, NZ_C
               dz_old(k) = 20.0_wp
               q_old(k) = 5.0_wp + 0.25_wp*real(k, wp)
            end do
            do k = 1, NZ_C
               dz_new(k) = 20.0_wp + 3.0_wp*sin(0.9_wp*real(k, wp))
            end do
            call rebalance(NZ_C, dz_new, sum(dz_old))
            select case (ishape)
            case (2)
               call vanish_layers(NZ_C, dz_new, 1, 3, sum(dz_old))
            case (3)
               call vanish_layers(NZ_C, dz_new, NZ_C - 2, NZ_C, sum(dz_old))
            case (4)
               call vanish_layers(NZ_C, dz_new, 8, 8, sum(dz_old))
            end select
            do im = 1, N_METHOD
               call remap_column(M_CODE(im), NZ_C, dz_old, dz_new, q_old, q_new)
               total_old = sum(q_old*dz_old)
               total_new = sum(q_new*dz_new)
               tol = 1.0e-11_wp*abs(total_old)
               if (.not. is_finite(total_new) .or. &
                   abs(total_new - total_old) > tol) then
                  write (msg, '(a,i0,a,es22.15,a,es22.15)') &
                     "remap "//trim(M_NAME(im))// &
                     " R2 conservation (shape ", ishape, "): sum(q*dz) went ", &
                     total_old, " -> ", total_new
                  call check(error, .false., trim(msg))
                  exit sweep
               end if
            end do
         end do
      end block sweep
   end subroutine test_remap_conservation

   subroutine test_remap_extrema(error)
      !! R3.  No method may create a value outside the old column's range.
      !! Asserted on a step profile — the hardest case for a
      !! reconstruction, and the one an overflow or a lock exchange puts in
      !! front of the remap on every step.
      !!
      !! The bound carries a round-off allowance of `1e-12·range`, not a
      !! tuned slack: an exact limiter clips AT the neighbouring cell
      !! average, and the overlap integral that follows is a sum of products
      !! whose last bit is not the limiter's.
      !!
      !! `REMAP_PQM` is INCLUDED rather than excused as "4th order, so
      !! overshoot is expected": White & Adcroft (2008) present PQM with a
      !! full monotonicity limiter, so an overshoot here would be a defect
      !! and is gated as one.  Measured: all five methods hold.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_C = 20
      character(len=220) :: msg
      real(wp) :: dz_old(NZ_C), dz_new(NZ_C), q_old(NZ_C), q_new(NZ_C)
      real(wp) :: lo, hi, slack
      integer :: im, k
      sweep: block
         do k = 1, NZ_C
            dz_old(k) = 25.0_wp
            if (k <= NZ_C/2) then
               q_old(k) = 34.0_wp
            else
               q_old(k) = 35.0_wp
            end if
            dz_new(k) = 25.0_wp + 6.0_wp*sin(1.3_wp*real(k, wp))
         end do
         call rebalance(NZ_C, dz_new, sum(dz_old))
         lo = minval(q_old)
         hi = maxval(q_old)
         slack = 1.0e-12_wp*(hi - lo)
         do im = 1, N_METHOD
            call remap_column(M_CODE(im), NZ_C, dz_old, dz_new, q_old, q_new)
            do k = 1, NZ_C
               if (q_new(k) < lo - slack .or. q_new(k) > hi + slack) then
                  write (msg, '(a,i0,a,es22.15,a,f8.4,a,f8.4,a)') &
                     "remap "//trim(M_NAME(im))//" R3 new extremum: q_new(", k, &
                     ") = ", q_new(k), " is outside the old range [", lo, ", ", &
                     hi, "]"
                  call check(error, .false., trim(msg))
                  exit sweep
               end if
            end do
         end do
      end block sweep
   end subroutine test_remap_extrema

   subroutine test_remap_linear_uniform(error)
      !! R4 on a UNIFORM old column: every method from PLM up reproduces a
      !! profile linear in `z` to round-off, in the interior AND in the
      !! boundary cells.  This is the order-of-accuracy check the ALE
      !! literature states (Adcroft & Hallberg 2006; White & Adcroft 2008).
      !!
      !! `REMAP_PCM` is excluded and that is not a loosening: a
      !! piecewise-CONSTANT reconstruction is first order everywhere by
      !! definition (measured here at 2.2e-02), so asserting linear
      !! exactness on it would be asserting that the method is not what its
      !! name says.  Its conservation and monotonicity are gated above.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_C = 20
      real(wp), parameter :: SLOPE = 0.013_wp, BASE = 34.2_wp
      real(wp), parameter :: TOL_LIN = 1.0e-10_wp
      character(len=220) :: msg
      real(wp) :: dz_old(NZ_C), dz_new(NZ_C), q_old(NZ_C), q_new(NZ_C)
      real(wp) :: z_new(0:NZ_C), worst
      integer :: im, k_worst
      sweep: block
         call linear_column(NZ_C, .true., dz_old, dz_new, q_old, z_new, SLOPE, BASE)
         do im = 2, N_METHOD
            call remap_column(M_CODE(im), NZ_C, dz_old, dz_new, q_old, q_new)
            call worst_linear_error(NZ_C, q_new, z_new, SLOPE, BASE, 1, NZ_C, &
                                    worst, k_worst)
            if (worst > TOL_LIN) then
               write (msg, '(a,es12.5,a,i0,a,es12.5)') &
                  "remap "//trim(M_NAME(im))// &
                  " R4 (uniform old column): worst error ", worst, " at k = ", &
                  k_worst, ", bar ", TOL_LIN
               call check(error, .false., trim(msg))
               exit sweep
            end if
         end do
      end block sweep
   end subroutine test_remap_linear_uniform

   subroutine test_remap_linear_stretched(error)
      !! R4 on a STRETCHED old column — which is what every ALE step in this
      !! model actually hands the remap (a sigma column over a slope, a z*
      !! column with fillers, a `pred_corr` corrector after the layers have
      !! moved).  `ppm_h4` and `pqm` carry the non-uniform h4 edge stencil
      !! of White & Adcroft (2008) and stay exact here; `plm` and `ppm` do
      !! not, and are pinned separately in
      !! `documents_plm_ppm_lose_linearity_when_stretched`.
      !!
      !! Boundary cells are INCLUDED in the bound: measured 8.5e-14 for
      !! `ppm_h4` and 1.1e-13 for `pqm`, i.e. the PCM boundary closure does
      !! not cost these two their linear exactness on this column.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_C = 20
      real(wp), parameter :: SLOPE = 0.013_wp, BASE = 34.2_wp
      real(wp), parameter :: TOL_LIN = 1.0e-10_wp
      character(len=220) :: msg
      real(wp) :: dz_old(NZ_C), dz_new(NZ_C), q_old(NZ_C), q_new(NZ_C)
      real(wp) :: z_new(0:NZ_C), worst
      integer :: im, k_worst
      sweep: block
         call linear_column(NZ_C, .false., dz_old, dz_new, q_old, z_new, SLOPE, BASE)
         do im = 4, N_METHOD          ! ppm_h4, pqm
            call remap_column(M_CODE(im), NZ_C, dz_old, dz_new, q_old, q_new)
            call worst_linear_error(NZ_C, q_new, z_new, SLOPE, BASE, 1, NZ_C, &
                                    worst, k_worst)
            if (worst > TOL_LIN) then
               write (msg, '(a,es12.5,a,i0,a,es12.5)') &
                  "remap "//trim(M_NAME(im))// &
                  " R4 (stretched old column): worst error ", worst, " at k = ", &
                  k_worst, ", bar ", TOL_LIN
               call check(error, .false., trim(msg))
               exit sweep
            end if
         end do
      end block sweep
   end subroutine test_remap_linear_stretched

   subroutine test_remap_linear_stretched_defect(error)
      !! `documents_*` — this asserts a KNOWN DEFECT on purpose.
      !!
      !! `remap_column_plm` limits with
      !! `slope(k) = 0.5·minmod(q(k+1) − q(k), q(k) − q(k−1))`, and
      !! `remap_column_ppm` interpolates its edge values with the matching
      !! 4-point formula.  Both are the UNIFORM-GRID forms: on a stretched
      !! column the two one-sided differences of a LINEAR profile are not
      !! equal, minmod takes the smaller, and the reconstruction
      !! under-estimates the slope.  A linear profile is then not
      !! reproduced, so the scheme is not second order on the grids it is
      !! actually run on.
      !!
      !! Measured (20 layers, `q = 34.2 + 0.013·z`, old `dz` stretched
      !! ±20 %, gfortran 15.1 Release):
      !!
      !!   method   worst interior (3..nz-2)   worst boundary (1, nz)
      !!   plm      2.25538e-03                9.48601e-04
      !!   ppm      6.28944e-04                8.27615e-04
      !!
      !! against `ppm_h4`'s 5.7e-14 on the identical column.  In salinity
      !! units that is 2.3 mPSU of error per remap on a 0.013 PSU/m
      !! gradient, on a coordinate that regrids every thermo step — and the
      !! shipped default is `remap_method = ppm`.
      !!
      !! The intended value is 0 for both.  A fix (the non-uniform PLM slope
      !! of Colella & Woodward 1984, and the h4 edge stencil `ppm_h4`
      !! already has) must DELETE this case and move `plm` / `ppm` into
      !! `remap_h4_stencils_linear_exact_when_stretched`.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_C = 20
      real(wp), parameter :: SLOPE = 0.013_wp, BASE = 34.2_wp
      real(wp), parameter :: TOL_LIN = 1.0e-10_wp
      character(len=220) :: msg
      real(wp) :: dz_old(NZ_C), dz_new(NZ_C), q_old(NZ_C), q_new(NZ_C)
      real(wp) :: z_new(0:NZ_C), worst
      integer :: im, k_worst
      sweep: block
         call linear_column(NZ_C, .false., dz_old, dz_new, q_old, z_new, SLOPE, BASE)
         do im = 2, 3                 ! plm, ppm
            call remap_column(M_CODE(im), NZ_C, dz_old, dz_new, q_old, q_new)
            call worst_linear_error(NZ_C, q_new, z_new, SLOPE, BASE, 3, NZ_C - 2, &
                                    worst, k_worst)
            if (worst <= TOL_LIN) then
               write (msg, '(a,es12.5,a)') &
                  "remap "//trim(M_NAME(im))//" is now linear-exact on a "// &
                  "stretched column (worst ", worst, &
                  ") -- the uniform-grid stencil is fixed. DELETE this case."
               call check(error, .false., trim(msg))
               exit sweep
            end if
            ! ...and it must not get WORSE than what was measured.
            if (worst > 1.0e-2_wp) then
               write (msg, '(a,es12.5,a)') &
                  "remap "//trim(M_NAME(im))//" stretched-column error ", worst, &
                  " is an order above the pinned 2.3e-03/6.3e-04 -- a "// &
                  "regression, not the documented defect."
               call check(error, .false., trim(msg))
               exit sweep
            end if
         end do
      end block sweep
   end subroutine test_remap_linear_stretched_defect

   pure subroutine worst_linear_error(nz, q_new, z_new, slope, base, k_lo, k_hi, &
                                      worst, k_worst)
      !! Worst absolute departure of `q_new` from the exact cell average of
      !! `base + slope·z` over `k_lo … k_hi`.
      integer, intent(in) :: nz, k_lo, k_hi
      real(wp), intent(in) :: q_new(nz)
      real(wp), intent(in) :: z_new(0:nz)
      real(wp), intent(in) :: slope, base
      real(wp), intent(out) :: worst
      integer, intent(out) :: k_worst
      real(wp) :: q_want
      integer :: k
      worst = 0.0_wp
      k_worst = 0
      do k = k_lo, k_hi
         q_want = base + slope*0.5_wp*(z_new(k - 1) + z_new(k))
         if (abs(q_new(k) - q_want) > worst) then
            worst = abs(q_new(k) - q_want)
            k_worst = k
         end if
      end do
   end subroutine worst_linear_error

end module test_ocean_vcoord_invariants
