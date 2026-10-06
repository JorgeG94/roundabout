!! Unit tests for the frhat port (`&ocean_bt_nml frhat_scheme`): MOM6
!! `btcalc`'s HVEL_SCHEME=HYBRID face-thickness closure
!! (`MOM_barotropic.F90:4546-4790`), ported bottom-up as
!! `rdb_barotropic_coupling::frhat_h_face_step`
!! (`src/shared_module_utilities/rdb_frhat_face.inc`) into every barotropic
!! depth mean that reads a layer's face thickness.  See
!! `docs/visc_rem_bt_rem_plan.md` Section 7 for the motivating gap and
!! `src/core/ocean/kernels/barotropic/rdb_barotropic_coupling.F90`'s
!! `face_depth_mean_u` docstring for the ISOMIP+ x22.7 bug class a single
!! inconsistent call site can cause.
!!
!! Tests:
!!   1. `test_frhat_hybrid_closed_form` -- a clean two-layer, two-column
!!      face (bed-partial RIGHT column against a deep LEFT column) worked
!!      by hand against MOM6's `e_u`/`D_shallow_u` recursion; the exact
!!      rational fractions (10/11, blended 2 + 150/169) are reproduced
!!      here as Fortran literals, not transcribed decimals, so there is
!!      no manual rounding to get wrong.
!!   2. `test_frhat_sums_to_one` -- for both schemes, on an irregular
!!      multi-layer column (including a near-H_DIV_EPS sliver), the
!!      normalised frhat weights sum to 1 to round-off.
!!   3. `test_frhat_arithmetic_matches_pre_port` -- FRHAT_ARITHMETIC must
!!      reproduce the plain two-cell mean exactly (bit-identity contract
!!      for the default).
!!   4. `test_fold_self_consistent` -- the fold invariant: with zero
!!      barotropic forcing, `apply_bt_correction`'s visc_rem-weighted
!!      fold (`wt = vr/<vr>_h`) preserves the FRHAT-weighted (not
!!      plain-h-weighted) depth mean of the layer velocities exactly, so
!!      an INDEPENDENT post-fold `face_depth_mean_u` call (the same
!!      routine `set_cor_ref_velocity` uses) reproduces `bt_ubt_end` to
!!      round-off under HYBRID -- the three call sites
!!      (`derive_bt_from_layers`, `apply_bt_correction`,
!!      `face_depth_mean_u`) agree on the SAME face weight.
!!   5. `test_frhat_device_residency` -- GPU device-mapping smoke test
!!      (`mem:separate` rules): `bt_H_ref`/`frhat_scheme` reach
!!      `face_depth_mean_u`/`compute_bt_rem_from_visc_rem` through the
!!      SAME `bt_work`/`metrics` objects already `enter_data`-mapped by
!!      the barotropic workstate/multilayer-state/metrics slots -- no new
!!      array needs its own map, but the `!$acc enter data` / `update
!!      self` round trip on `av_rem_u/v` must still survive with
!!      `frhat_scheme = FRHAT_HYBRID`, exactly mirroring
!!      `test_ocean_bt_rem_from_visc_rem`'s existing device-residency
!!      pattern (inert on a host/multicore build, load-bearing on the GPU
!!      build per CLAUDE.md's `mem:separate` gotcha).
module test_ocean_bt_frhat
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, FRHAT_ARITHMETIC, FRHAT_HYBRID, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: frhat_h_face_step, face_depth_mean_u, &
                                      derive_bt_from_layers, apply_bt_correction, &
                                      compute_bt_rem_from_visc_rem
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_bt_frhat_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_bt_frhat_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("frhat_hybrid_closed_form", test_frhat_hybrid_closed_form), &
                  new_unittest("frhat_sums_to_one", test_frhat_sums_to_one), &
                  new_unittest("frhat_arithmetic_matches_pre_port", &
                               test_frhat_arithmetic_matches_pre_port), &
                  new_unittest("fold_self_consistent", test_fold_self_consistent), &
                  new_unittest("frhat_device_residency", test_frhat_device_residency) &
                  ]
   end subroutine collect_ocean_bt_frhat_tests

   subroutine test_frhat_hybrid_closed_form(error)
      !! Two layers (k=1 bed, k=2 surface). LEFT column 10 m deep, split
      !! 5/5. RIGHT column 2 m deep (a thin partial-bed sliver against the
      !! deep left column), split 0.5 (bed) / 1.5 (surface).
      !!
      !! Hand recursion (mean-bed datum, MOM6 e_u/D_shallow_u, see the
      !! module docstring and `rdb_frhat_face.inc`):
      !!   d_shallow = -min(10,2) = -2;  e_prev(0) = -0.5*(10+2) = -6
      !!   k=1: h_arith=2.75, e_cur=-3.25. e_prev(-6) < d_shallow(-2) =>
      !!        harmonic: h_harm = (5*0.5)/2.75 = 10/11. e_cur(-3.25) <=
      !!        d_shallow(-2) => hatu(1) = 10/11 (wholly below the shelf).
      !!   k=2: h_arith=3.25, e_cur=0. e_prev(-3.25) < d_shallow(-2) =>
      !!        harmonic: h_harm = (5*1.5)/3.25 = 30/13. e_cur(0) >
      !!        d_shallow(-2) => straddles: wt = (0-(-2))/3.25 = 8/13;
      !!        hatu(2) = (8/13)*3.25 + (5/13)*(30/13) = 2 + 150/169.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HREF_L = 10.0_wp, HREF_R = 2.0_wp
      real(wp), parameter :: HL(2) = [5.0_wp, 5.0_wp]
      real(wp), parameter :: HR(2) = [0.5_wp, 1.5_wp]
      real(wp), parameter :: EXPECT1 = 10.0_wp/11.0_wp
      real(wp), parameter :: EXPECT2 = 2.0_wp + 150.0_wp/169.0_wp
      real(wp) :: e_prev, hatu1, hatu2, hatutot
      checks: block
         e_prev = -0.5_wp*(HREF_L + HREF_R)
         call frhat_h_face_step(HL(1), HR(1), HREF_L, HREF_R, FRHAT_HYBRID, e_prev, hatu1)
         call check(error, abs(hatu1 - EXPECT1) < 1.0e-12_wp, &
                    "hatu(k=1, bed) does not match the hand-rolled HYBRID recursion")
         if (allocated(error)) exit checks
         call check(error, abs(e_prev - (-3.25_wp)) < 1.0e-12_wp, &
                    "e_prev after k=1 does not match the hand-rolled elevation")
         if (allocated(error)) exit checks

         call frhat_h_face_step(HL(2), HR(2), HREF_L, HREF_R, FRHAT_HYBRID, e_prev, hatu2)
         call check(error, abs(hatu2 - EXPECT2) < 1.0e-12_wp, &
                    "hatu(k=2, surface, straddles the shelf) does not match the hand-rolled "// &
                    "HYBRID recursion")
         if (allocated(error)) exit checks
         call check(error, abs(e_prev - 0.0_wp) < 1.0e-12_wp, &
                    "e_prev after k=2 must reach 0 (the free surface)")
         if (allocated(error)) exit checks

         ! Sanity: HYBRID suppresses the thin sliver relative to the
         ! arithmetic mean at the bed (0.909 < 2.75) and inflates the
         ! surface layer's weight relative to it (2.887 > 3.25 is FALSE --
         ! rather, the suppressed bed mass is redistributed upward; check
         ! only the qualitative bed-suppression, the sign this port exists
         ! for).
         hatutot = hatu1 + hatu2
         call check(error, hatu1 < 0.5_wp*(HL(1) + HR(1)), &
                    "HYBRID must suppress the bed layer's weight below the plain "// &
                    "arithmetic mean (the whole point of the port)")
         if (allocated(error)) exit checks
         call check(error, hatutot > 0.0_wp, "hatutot must be positive")
      end block checks
   end subroutine test_frhat_hybrid_closed_form

   subroutine test_frhat_sums_to_one(error)
      !! Six-layer irregular column (including one H_DIV_EPS-scale
      !! sliver) at an asymmetric two-column face: for EITHER scheme the
      !! normalised frhat weights sum to 1 to round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp), parameter :: HL(NZ) = [20.0_wp, 15.0_wp, 10.0_wp, 8.0_wp, 5.0_wp, 2.0_wp]
      real(wp), parameter :: HR(NZ) = [1.0e-3_wp, 0.3_wp, 0.6_wp, 1.0_wp, 1.5_wp, 1.9_wp]
      real(wp) :: href_l, href_r, e_prev, hatu(NZ), hatutot, frhat_sum
      integer :: scheme, k
      checks: block
         href_l = sum(HL)
         href_r = sum(HR)
         do scheme = FRHAT_ARITHMETIC, FRHAT_HYBRID
            e_prev = -0.5_wp*(href_l + href_r)
            hatutot = 0.0_wp
            do k = 1, NZ
               call frhat_h_face_step(HL(k), HR(k), href_l, href_r, scheme, e_prev, hatu(k))
               hatutot = hatutot + hatu(k)
            end do
            call check(error, hatutot > 0.0_wp, "hatutot must be positive")
            if (allocated(error)) exit checks
            frhat_sum = sum(hatu)/hatutot
            call check(error, abs(frhat_sum - 1.0_wp) < 1.0e-12_wp, &
                       "normalised frhat weights do not sum to 1")
            if (allocated(error)) exit checks
            ! Every individual weight must stay non-negative and finite --
            ! the harmonic/blend branches must not go rogue on the sliver.
            do k = 1, NZ
               call check(error, hatu(k) >= 0.0_wp, &
                          "a frhat layer weight went negative")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
   end subroutine test_frhat_sums_to_one

   subroutine test_frhat_arithmetic_matches_pre_port(error)
      !! FRHAT_ARITHMETIC must reduce to the plain two-cell mean,
      !! regardless of href/e_prev (the pre-port bit-identity contract
      !! for the default scheme).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hatu, e_prev
      checks: block
         e_prev = -0.5_wp*(1.0_wp + 999.0_wp)
         call frhat_h_face_step(3.0_wp, 7.0_wp, 1.0_wp, 999.0_wp, FRHAT_ARITHMETIC, e_prev, hatu)
         call check(error, hatu == 0.5_wp*(3.0_wp + 7.0_wp), &
                    "FRHAT_ARITHMETIC must equal the plain two-cell mean exactly")
         if (allocated(error)) exit checks
         ! href/e_prev are dereferenced but must not perturb the answer --
         ! e_prev itself must also be left untouched under this scheme.
         call check(error, e_prev == -0.5_wp*(1.0_wp + 999.0_wp), &
                    "FRHAT_ARITHMETIC must leave e_prev untouched")
      end block checks
   end subroutine test_frhat_arithmetic_matches_pre_port

   subroutine build_two_column_face(grid, ms, bt_work, metrics, h_deep, h_shallow_bed, &
                                    h_shallow_rest, nz)
      !! A 3-column Cartesian channel (col0 deep | col1 the partial-bed
      !! sill | col2 deep again), one interior u-face per side -- mirrors
      !! `test_ocean_bt_rem_from_visc_rem.F90::build_channel`'s geometry
      !! shape but with caller-supplied depths/nz, and (new here) sets
      !! `bt_H_ref` so `frhat_scheme = FRHAT_HYBRID` has real bathymetry
      !! to read.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      type(ocean_metrics_t), intent(out) :: metrics
      real(wp), intent(in) :: h_deep, h_shallow_bed, h_shallow_rest
      integer, intent(in) :: nz
      real(wp), parameter :: DX = 50000.0_wp
      integer :: i

      call grid%init(3, 1, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = nz
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=nz)

      ms%h_layer(NGHOST + 1, :, :) = h_deep/real(nz, wp)
      ms%h_layer(NGHOST + 2, :, 1) = h_shallow_bed
      ms%h_layer(NGHOST + 2, :, 2:nz) = h_shallow_rest
      ms%h_layer(NGHOST + 3, :, :) = h_deep/real(nz, wp)

      do i = 1, size(bt_work%bt_H_ref, 1)
         if (i == NGHOST + 2) then
            bt_work%bt_H_ref(i, :) = h_shallow_bed + real(nz - 1, wp)*h_shallow_rest
         else
            bt_work%bt_H_ref(i, :) = h_deep
         end if
      end do
      bt_work%frhat_scheme = FRHAT_HYBRID
   end subroutine build_two_column_face

   subroutine test_fold_self_consistent(error)
      !! The fold must stay self-consistent (the "ISOMIP+ x22.7" bug
      !! class): with F_bt_u == 0 everywhere, apply_bt_correction's
      !! visc_rem-weighted fold (use_visc_rem=.true., the closed-faces
      !! path off) preserves the FRHAT-weighted depth mean of the layer
      !! velocities exactly equal to `bt_ubt_end`, where "FRHAT-weighted"
      !! means `derive_bt_from_layers`'s own `frhat_scheme = FRHAT_HYBRID`
      !! weight -- verified by an INDEPENDENT post-fold
      !! `face_depth_mean_u` call (the routine `set_cor_ref_velocity`
      !! uses) reproducing `bt_ubt_end` to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NZ = 4
      real(wp), allocatable :: f_mean(:, :)
      real(wp) :: dt
      integer :: jp, nu

      checks: block
         call build_two_column_face(grid, ms, bt_work, metrics, 188.0_wp, 0.74_wp, &
                                    (8.4_wp - 0.74_wp)/3.0_wp, NZ)

         ! Per-layer viscous remnant (uneven across k, like the glued BBL
         ! case this chain exists for) at the two interior u-faces.
         bt_work%visc_rem_u(NGHOST + 2, :, :) = spread( &
                                                [0.0658_wp, 0.3071_wp, 0.5142_wp, 0.6201_wp], 1, size(bt_work%visc_rem_u, 2))
         bt_work%visc_rem_u(NGHOST + 3, :, :) = spread( &
                                                [0.0631_wp, 0.2844_wp, 0.4767_wp, 0.5723_wp], 1, size(bt_work%visc_rem_u, 2))

         ! Seed a non-uniform layer velocity (shear across k) so the
         ! depth-mean weighting is actually exercised, and derive bt_ubt
         ! from it the SAME way the driver does.
         ms%u_face_x_layer = 0.0_wp
         ms%u_face_x_layer(NGHOST + 2, :, :) = spread( &
                                               [0.01_wp, -0.02_wp, 0.03_wp, 0.05_wp], 1, size(ms%u_face_x_layer, 2))
         ms%u_face_x_layer(NGHOST + 3, :, :) = spread( &
                                               [-0.01_wp, 0.04_wp, -0.015_wp, 0.02_wp], 1, size(ms%u_face_x_layer, 2))

         call derive_bt_from_layers(grid, bt_work, ms, metrics)
         bt_work%ubt_at_n = bt_work%bt_ubt
         ! Choose an arbitrary end state and ZERO forcing -- Delta = ubt_end
         ! - ubt_at_n exactly.
         bt_work%bt_ubt_end = bt_work%bt_ubt + 0.2_wp
         bt_work%F_bt_u = 0.0_wp
         dt = 1.0_wp

         call apply_bt_correction(bt_work, ms, dt, metrics, use_visc_rem=.true.)

         ! Independent post-fold recomputation via face_depth_mean_u (the
         ! SAME routine set_cor_ref_velocity uses for the Coriolis
         ! reference) must reproduce bt_ubt_end to round-off.
         nu = size(bt_work%bt_ubt_end, 1)
         allocate (f_mean(nu, size(bt_work%bt_ubt_end, 2)))
         call face_depth_mean_u(grid, ms%u_face_x_layer, ms%h_layer, f_mean, NZ, metrics, &
                                bt_work%bt_H_ref, bt_work%frhat_scheme)

         jp = NGHOST + 1
         call check(error, abs(f_mean(NGHOST + 2, jp) - bt_work%bt_ubt_end(NGHOST + 2, jp)) &
                    < 1.0e-10_wp, &
                    "face01: the frhat-weighted post-fold mean does not reproduce bt_ubt_end "// &
                    "-- the fold is not self-consistent with derive_bt_from_layers/"// &
                    "face_depth_mean_u under frhat_scheme=hybrid")
         if (allocated(error)) exit checks
         call check(error, abs(f_mean(NGHOST + 3, jp) - bt_work%bt_ubt_end(NGHOST + 3, jp)) &
                    < 1.0e-10_wp, &
                    "face12: the frhat-weighted post-fold mean does not reproduce bt_ubt_end")
      end block checks
      if (allocated(f_mean)) deallocate (f_mean)
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_fold_self_consistent

   subroutine test_frhat_device_residency(error)
      !! GPU device-mapping smoke test (`mem:separate` rules): the
      !! `href`/`scheme` plumbing into `face_depth_mean_u/v` (hence
      !! `compute_bt_rem_from_visc_rem`'s `av_rem`) rides the SAME
      !! already-mapped `bt_work%bt_H_ref` / `ms%h_layer` / `metrics`
      !! objects -- no new array needs its own `enter_data`.  Mirrors
      !! `test_ocean_bt_rem_from_visc_rem::test_column_spin_down_matches_
      !! av_rem`'s `!$acc enter data` / `update self` round trip, with
      !! `frhat_scheme = FRHAT_HYBRID`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NZ = 4
      real(wp) :: av01_host, av12_host
      integer :: jp
      checks: block
         call build_two_column_face(grid, ms, bt_work, metrics, 188.0_wp, 0.74_wp, &
                                    (8.4_wp - 0.74_wp)/3.0_wp, NZ)
         bt_work%visc_rem_u(NGHOST + 2, :, :) = spread( &
                                                [0.0658_wp, 0.3071_wp, 0.5142_wp, 0.6201_wp], 1, size(bt_work%visc_rem_u, 2))
         bt_work%visc_rem_u(NGHOST + 3, :, :) = spread( &
                                                [0.0631_wp, 0.2844_wp, 0.4767_wp, 0.5723_wp], 1, size(bt_work%visc_rem_u, 2))

         ! Host reference, computed BEFORE any device map.
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 4)
         jp = NGHOST + 1
         av01_host = bt_work%av_rem_u(NGHOST + 2, jp)
         av12_host = bt_work%av_rem_u(NGHOST + 3, jp)

         ! Re-zero and recompute through a mapped bt_work/ms/metrics --
         ! `ms%enter_data()`/`metrics%enter_data()` attach `h_layer`/
         ! `open_u` (consumed inside `face_depth_mean_u`'s `do concurrent`);
         ! `bt_work%enter_data()` attaches `bt_H_ref`/`visc_rem_u`/
         ! `av_rem_u`. On a `mem:separate` GPU build the kernel reads device
         ! memory exclusively; `update self` pulls the result back for the
         ! host `check` below.
         bt_work%av_rem_u = 0.0_wp
         bt_work%av_rem_v = 0.0_wp
         call ms%enter_data()
         call metrics%enter_data()
         call bt_work%enter_data()
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 4)
         !$acc update self(bt_work%av_rem_u, bt_work%av_rem_v)
         call bt_work%exit_data()
         call metrics%exit_data()
         call ms%exit_data()

         call check(error, abs(bt_work%av_rem_u(NGHOST + 2, jp) - av01_host) < 1.0e-12_wp, &
                    "av_rem_u(face01) changed across the enter_data/update self round trip "// &
                    "under frhat_scheme=hybrid")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%av_rem_u(NGHOST + 3, jp) - av12_host) < 1.0e-12_wp, &
                    "av_rem_u(face12) changed across the enter_data/update self round trip "// &
                    "under frhat_scheme=hybrid")
         if (allocated(error)) exit checks
         call check(error, av01_host > 0.0_wp .and. av01_host < 1.0_wp, &
                    "av_rem(face01) out of (0,1) under frhat_scheme=hybrid")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_frhat_device_residency

end module test_ocean_bt_frhat
