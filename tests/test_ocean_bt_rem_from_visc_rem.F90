!! Unit tests for PR-2 (bt-rem-from-av-rem): `bt_rem_u/v` built from the
!! viscous remnant the layered momentum solve uses (`&ocean_bt_nml
!! bt_rem_from_visc_rem`), MOM6 `MOM_barotropic.F90:1553-1580`.
!!
!! `av_rem_u/v := Sum_k frhat_k*visc_rem_k` (`compute_bt_rem_from_visc_rem`,
!! `rdb_barotropic_coupling.F90`) reuses `face_depth_mean_u/v`'s own
!! arithmetic-mean face-thickness weight (`frhat_k = h_face_k /
!! Sum_k h_face_k`) directly on `visc_rem_u/v` -- NOT
!! `face_depth_mean_rem_u/v`'s `wt_u = h_face*visc_rem` weighting (that
!! one is the MOM6 FORCING weight, `forcing_visc_rem`/PR-3 scope).  Then
!! `bt_rem = av_rem**(1/n_inner)` (plain) or the `strong_drag` rational
!! form, masked to 0 where `av_rem <= 0`.
!!
!! The two-cell channel geometry + per-layer visc_rem profiles in
!! `test_slosh_bounded_with_chain`/`test_slosh_undamped_without_chain`
!! are taken directly from `python_prototypes/bt_rem/README.md`'s
!! `fastloop_damping.py` positive result (3 T-cells: 188 m / 8.4 m sill
!! with a 0.74 m partial bed cell / 70 m, nz=4, bed->surface visc_rem
!! `[0.0658, 0.3071, 0.5142, 0.6201]` / `[0.0631, 0.2844, 0.4767,
!! 0.5723]`).  The resulting `av_rem` differs numerically from the
!! prototype's own 0.44/0.41 (roundabout's `face_depth_mean_u` uses a
!! plain arithmetic-mean `h_face`, not necessarily the exact weighting
!! the prototype's own `bt_rem_lib.py` used for `frhat`) -- the per-face
!! `av_rem` asserted here is computed independently in-test by the SAME
!! arithmetic-mean formula `face_depth_mean_u` uses, so the check is a
!! genuine implementation-matches-its-own-spec test, not a transcription
!! of the prototype's number.  The qualitative mechanism (undamped mode
!! stays bounded/neutral, chain-damped mode decays by orders of
!! magnitude over repeated outer steps) is the thing being reproduced.
!!
!! Tests:
!!   1. `test_av_rem_matches_hand_calc` -- `av_rem_u` at both interior
!!      faces of the 3-column channel matches an independently
!!      hand-rolled depth mean (same arithmetic as `face_depth_mean_u`)
!!      to round-off.
!!   2. `test_bt_rem_power_closed_form` -- `bt_rem_u**n_inner == av_rem_u`
!!      to round-off (the plain power form), for several `n_inner`.
!!   3. `test_bt_rem_strong_drag_closed_form` -- the `strong_drag`
!!      rational form matches `n*av/(1+(n-1)*av)` to round-off, and its
!!      `n_inner=1` case reduces exactly to `av_rem` (MOM6 analytic
!!      check #2).
!!   4. `test_bt_rem_mask_dry_face` -- a face with zero thickness on
!!      both sides (`av_rem <= 0`) gets `bt_rem = 0`, not 1 or NaN, under
!!      BOTH forms.
!!   5. `test_column_spin_down_matches_av_rem` -- drives the REAL
!!      `barotropic_substep_nonlinear_interior`: from a uniform rest
!!      column with a nonzero `bt_ubt` IC and zero forcing, ONE call
!!      with `n_steps = n_inner` must contract `bt_ubt` by EXACTLY
!!      `bt_rem**n_inner = av_rem` (the depth-mean decay through the BT
!!      path equals the layered remnant mean, by construction) --
!!      exercises the `!$acc enter data`/`update self` round trip too
!!      (GPU device-residency smoke test for `av_rem_u/v`/`bt_rem_u/v`).
!!   6. `test_slosh_bounded_with_chain` / `test_slosh_undamped_without_chain`
!!      -- the two-cell thin-channel slosh: a 2*dx eta perturbation on
!!      the 3-column channel, run for several repeated outer-step blocks
!!      with NO forcing/Coriolis.  With the chain on (`bt_rem` from
!!      `av_rem`) the interior |u| decays by orders of magnitude over
!!      the run; with `bt_rem == 1` (today's default) it stays within a
!!      bounded, non-decaying range -- never anywhere near the damped
!!      run's final amplitude.
module test_ocean_bt_rem_from_visc_rem
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: compute_bt_rem_from_visc_rem
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear_interior
   implicit none
   private

   public :: collect_bt_rem_from_visc_rem_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   ! Prototype channel geometry (python_prototypes/bt_rem/README.md):
   ! col0 (188 m) -- face01 -- col1 (8.4 m sill, 0.74 m partial bed cell)
   ! -- face12 -- col2 (70 m).
   real(wp), parameter :: H_COL0 = 188.0_wp
   real(wp), parameter :: H_COL1_BED = 0.74_wp
   real(wp), parameter :: H_COL1_UPPER = (8.4_wp - 0.74_wp)/3.0_wp
   real(wp), parameter :: H_COL2 = 70.0_wp
   ! Prototype per-layer visc_rem profiles, bed (k=1) -> surface (k=4).
   real(wp), parameter :: VR_FACE01(4) = [0.0658_wp, 0.3071_wp, 0.5142_wp, 0.6201_wp]
   real(wp), parameter :: VR_FACE12(4) = [0.0631_wp, 0.2844_wp, 0.4767_wp, 0.5723_wp]

contains

   subroutine collect_bt_rem_from_visc_rem_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("av_rem_matches_hand_calc", test_av_rem_matches_hand_calc), &
                  new_unittest("bt_rem_power_closed_form", test_bt_rem_power_closed_form), &
                  new_unittest("bt_rem_strong_drag_closed_form", &
                               test_bt_rem_strong_drag_closed_form), &
                  new_unittest("bt_rem_mask_dry_face", test_bt_rem_mask_dry_face), &
                  new_unittest("column_spin_down_matches_av_rem", &
                               test_column_spin_down_matches_av_rem), &
                  new_unittest("slosh_bounded_with_chain", test_slosh_bounded_with_chain), &
                  new_unittest("slosh_undamped_without_chain", &
                               test_slosh_undamped_without_chain) &
                  ]
   end subroutine collect_bt_rem_from_visc_rem_tests

   subroutine build_channel(grid, ms, bt_work, metrics)
      !! The 3-column channel: physical cells at i = NGHOST+1 (col0),
      !! NGHOST+2 (col1, sill), NGHOST+3 (col2).  Interior u-faces at
      !! i = NGHOST+2 (face01) and NGHOST+3 (face12); wall faces at
      !! i = NGHOST+1 and NGHOST+4 are hard-zeroed by the substep itself.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      type(ocean_metrics_t), intent(out) :: metrics
      real(wp), parameter :: DX = 50000.0_wp

      call grid%init(3, 1, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=NZ)

      ! Column thicknesses (uniform in j; only j=NGHOST+1 is physical).
      ms%h_layer(NGHOST + 1, :, :) = H_COL0/real(NZ, wp)
      ms%h_layer(NGHOST + 2, :, 1) = H_COL1_BED
      ms%h_layer(NGHOST + 2, :, 2:4) = H_COL1_UPPER
      ms%h_layer(NGHOST + 3, :, :) = H_COL2/real(NZ, wp)

      ! Per-face visc_rem profiles at the two interior u-faces.
      bt_work%visc_rem_u(NGHOST + 2, :, :) = spread(VR_FACE01, 1, size(bt_work%visc_rem_u, 2))
      bt_work%visc_rem_u(NGHOST + 3, :, :) = spread(VR_FACE12, 1, size(bt_work%visc_rem_u, 2))
   end subroutine build_channel

   function hand_av_rem_u(h_west, h_east, vr) result(av_rem)
      !! Independent reference: the SAME arithmetic-mean h_face formula
      !! `face_depth_mean_u` uses, applied by hand to a two-column face.
      real(wp), intent(in) :: h_west(NZ), h_east(NZ), vr(NZ)
      real(wp) :: av_rem
      real(wp) :: h_face, num, denom
      integer :: k
      num = 0.0_wp; denom = 0.0_wp
      do k = 1, NZ
         h_face = 0.5_wp*(h_west(k) + h_east(k))
         num = num + vr(k)*h_face
         denom = denom + h_face
      end do
      av_rem = num/denom
   end function hand_av_rem_u

   subroutine test_av_rem_matches_hand_calc(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      ! Named distinctly from the module parameters H_COL0/H_COL1_BED/
      ! H_COL1_UPPER/H_COL2 above -- Fortran is case-insensitive, so a
      ! local `h_col0` would silently SHADOW the module parameter
      ! `H_COL0` inside this scope (bit us during development: `h_col0 =
      ! H_COL0/real(NZ,wp)` read the not-yet-initialized local on the
      ! RHS, not the 188.0 parameter, and gave exactly 0).
      real(wp) :: hcol0_arr(NZ), hcol1_arr(NZ), hcol2_arr(NZ)
      real(wp) :: expect01, expect12
      integer, parameter :: N_INNER = 4
      integer :: jp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, N_INNER)

         hcol0_arr = H_COL0/real(NZ, wp)
         hcol1_arr = [H_COL1_BED, H_COL1_UPPER, H_COL1_UPPER, H_COL1_UPPER]
         hcol2_arr = H_COL2/real(NZ, wp)
         expect01 = hand_av_rem_u(hcol0_arr, hcol1_arr, VR_FACE01)
         expect12 = hand_av_rem_u(hcol1_arr, hcol2_arr, VR_FACE12)

         jp = NGHOST + 1
         call check(error, abs(bt_work%av_rem_u(NGHOST + 2, jp) - expect01) < 1.0e-12_wp, &
                    "av_rem_u(face01) does not match the hand-rolled depth mean")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%av_rem_u(NGHOST + 3, jp) - expect12) < 1.0e-12_wp, &
                    "av_rem_u(face12) does not match the hand-rolled depth mean")
         if (allocated(error)) exit checks
         ! Sanity: both faces damp real momentum (0 < av_rem < 1), matching
         ! the prototype's qualitative finding (the glued bed layer retains
         ! only a fraction of a barotropic kick).
         call check(error, expect01 > 0.0_wp .and. expect01 < 1.0_wp, &
                    "av_rem(face01) out of (0,1) -- geometry/profile typo?")
         if (allocated(error)) exit checks
         call check(error, expect12 > 0.0_wp .and. expect12 < 1.0_wp, &
                    "av_rem(face12) out of (0,1) -- geometry/profile typo?")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_av_rem_matches_hand_calc

   subroutine test_bt_rem_power_closed_form(error)
      !! bt_rem**n_inner == av_rem exactly (to round-off), for several
      !! n_inner -- the MOM6 analytic check #2 from
      !! python_prototypes/bt_rem/README.md.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      integer, parameter :: N_CASES = 4
      integer, parameter :: N_INNERS(N_CASES) = [1, 4, 16, 64]
      integer :: c, jp
      real(wp) :: av_rem, roundtrip
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         jp = NGHOST + 1
         do c = 1, N_CASES
            call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, N_INNERS(c))
            av_rem = bt_work%av_rem_u(NGHOST + 2, jp)
            roundtrip = bt_work%bt_rem_u(NGHOST + 2, jp)**N_INNERS(c)
            call check(error, abs(roundtrip - av_rem) < 1.0e-10_wp*max(1.0_wp, av_rem), &
                       "bt_rem**n_inner /= av_rem (face01)")
            if (allocated(error)) exit checks
            av_rem = bt_work%av_rem_u(NGHOST + 3, jp)
            roundtrip = bt_work%bt_rem_u(NGHOST + 3, jp)**N_INNERS(c)
            call check(error, abs(roundtrip - av_rem) < 1.0e-10_wp*max(1.0_wp, av_rem), &
                       "bt_rem**n_inner /= av_rem (face12)")
            if (allocated(error)) exit checks
         end do
         ! n_inner = 1 must reduce bt_rem == av_rem exactly (Instep = 1).
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 1)
         call check(error, bt_work%bt_rem_u(NGHOST + 2, jp) == bt_work%av_rem_u(NGHOST + 2, jp), &
                    "n_inner=1: bt_rem must equal av_rem exactly")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bt_rem_power_closed_form

   subroutine test_bt_rem_strong_drag_closed_form(error)
      !! MOM6 BT_STRONG_DRAG rational form:
      !!   bt_rem = n_inner*av_rem / (1 + (n_inner-1)*av_rem)
      !! and its n_inner=1 reduction to av_rem (same analytic check,
      !! strong-drag variant).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp) :: av_rem, expect, rn
      integer, parameter :: N_INNER = 8
      integer :: jp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         bt_work%bt_strong_drag = .true.
         jp = NGHOST + 1

         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, N_INNER)
         av_rem = bt_work%av_rem_u(NGHOST + 2, jp)
         rn = real(N_INNER, wp)
         expect = (rn*av_rem)/(1.0_wp + (rn - 1.0_wp)*av_rem)
         call check(error, abs(bt_work%bt_rem_u(NGHOST + 2, jp) - expect) < 1.0e-13_wp, &
                    "strong_drag bt_rem does not match n*av/(1+(n-1)*av)")
         if (allocated(error)) exit checks

         ! n_inner = 1 reduces to av_rem exactly for the rational form too
         ! ((1*av)/(1+0*av) = av).
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 1)
         call check(error, abs(bt_work%bt_rem_u(NGHOST + 2, jp) - &
                               bt_work%av_rem_u(NGHOST + 2, jp)) < 1.0e-13_wp, &
                    "strong_drag n_inner=1 must reduce to av_rem")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bt_rem_strong_drag_closed_form

   subroutine test_bt_rem_mask_dry_face(error)
      !! A face with zero thickness on BOTH sides has av_rem = 0
      !! (face_depth_mean_u's denom<=0 fallback) -- bt_rem must be 0
      !! under both the plain-power and strong_drag forms, never 1 or
      !! NaN (CLAUDE.md: no max(..., eps) floor substitute).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      integer :: jp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         ! Dry out col0 and col1 entirely -- face01 (between them) is
         ! now a zero-thickness face on both sides.
         ms%h_layer(NGHOST + 1, :, :) = 0.0_wp
         ms%h_layer(NGHOST + 2, :, :) = 0.0_wp
         jp = NGHOST + 1

         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 4)
         call check(error, bt_work%av_rem_u(NGHOST + 2, jp) == 0.0_wp, &
                    "dry face: av_rem must be exactly 0")
         if (allocated(error)) exit checks
         call check(error, bt_work%bt_rem_u(NGHOST + 2, jp) == 0.0_wp, &
                    "dry face: bt_rem (plain form) must be exactly 0, not 1")
         if (allocated(error)) exit checks

         bt_work%bt_strong_drag = .true.
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 4)
         call check(error, bt_work%bt_rem_u(NGHOST + 2, jp) == 0.0_wp, &
                    "dry face: bt_rem (strong_drag form) must be exactly 0, not 1")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bt_rem_mask_dry_face

   subroutine test_column_spin_down_matches_av_rem(error)
      !! The depth-mean decay through the BT path equals the layered
      !! remnant mean: from a spatially uniform nonzero bt_ubt IC with
      !! zero forcing/Coriolis/eta (bt_H_ref left at its init 0, so the
      !! continuity PGF term never activates), a SINGLE substep
      !! (n_steps=1) must contract `bt_ubt` at each interior face by
      !! EXACTLY that face's `av_rem` (bt_rem built with n_inner=1, so
      !! `bt_rem == av_rem` with no exponent) -- the same single-step
      !! exactness argument as `test_substep_applies_bt_rem`
      !! (`tests/test_ocean_bt_substep_drag.F90`), now driven by the REAL
      !! visc_rem-derived `av_rem_u` rather than a hand-set constant.
      !! Also exercises the `!$acc enter data`/`update self` round trip
      !! on `av_rem_u/v`/`bt_rem_u/v` (GPU device-residency smoke test).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: f_corner(:, :)
      real(wp), parameter :: U0 = 0.1_wp, DT_INNER = 50.0_wp
      real(wp) :: av01, av12
      integer :: jp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, 1)
         av01 = bt_work%av_rem_u(NGHOST + 2, NGHOST + 1)
         av12 = bt_work%av_rem_u(NGHOST + 3, NGHOST + 1)

         bt_work%bt_ubt = U0
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1), source=0.0_wp)

         !$acc enter data copyin(f_corner)
         call bt_work%enter_data()
         call barotropic_substep_nonlinear_interior(grid, metrics, bt_work, f_corner, &
                                                    1, DT_INNER)
         !$acc update self(bt_work%bt_ubt, bt_work%av_rem_u, bt_work%bt_rem_u)
         call bt_work%exit_data()
         !$acc exit data delete(f_corner)

         jp = NGHOST + 1
         call check(error, bt_work%bt_ubt(NGHOST + 2, jp) == av01*U0, &
                    "face01: one substep did not contract bt_ubt by exactly av_rem")
         if (allocated(error)) exit checks
         call check(error, bt_work%bt_ubt(NGHOST + 3, jp) == av12*U0, &
                    "face12: one substep did not contract bt_ubt by exactly av_rem")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_column_spin_down_matches_av_rem

   subroutine run_slosh(bt_work, grid, metrics, f_corner, n_outer, n_inner, dt_inner, &
                        n_early, n_late, peak_early, peak_late)
      !! Seed the 2*dx eta perturbation (col0/col2 = +eta0, col1 = -eta0,
      !! matching `fastloop_damping.py`), at rest, no forcing/Coriolis,
      !! and run `n_outer` repeated blocks of `n_inner` substeps each
      !! (bt_rem_u/v held fixed across the whole run -- visc_rem never
      !! changes in this toy).  Reports the RUNNING-MAX interior |u| over
      !! the first `n_early` outer blocks (lets the eta->u gravity-wave
      !! conversion reach its natural early amplitude before anything
      !! has had much time to damp) and over the last `n_late` blocks
      !! (the long-run, asymptotic amplitude) -- comparing two single
      !! instantaneous samples is not robust here because the two
      !! interior faces have slightly different gravity-wave periods
      !! (188|8.4 m vs 8.4|70 m flanking depths), so a lone sample can
      !! land anywhere on the beat pattern even in the UNDAMPED case.
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
      integer, intent(in) :: n_outer, n_inner, n_early, n_late
      real(wp), intent(in) :: dt_inner
      real(wp), intent(out) :: peak_early, peak_late
      real(wp), parameter :: ETA0 = 0.02_wp
      integer :: step, jp
      real(wp) :: cur

      jp = NGHOST + 1
      bt_work%bt_eta(NGHOST + 1, jp) = ETA0
      bt_work%bt_eta(NGHOST + 2, jp) = -ETA0
      bt_work%bt_eta(NGHOST + 3, jp) = ETA0
      bt_work%bt_H_ref(NGHOST + 1, jp) = H_COL0
      bt_work%bt_H_ref(NGHOST + 2, jp) = 8.4_wp
      bt_work%bt_H_ref(NGHOST + 3, jp) = H_COL2

      peak_early = 0.0_wp
      peak_late = 0.0_wp

      ! GPU device-residency: map f_corner and bt_work's components
      ! (`bt_work%enter_data()`, never the whole derived type) before the repeated
      ! substep calls, and update the host copy of bt_ubt back only
      ! where it is actually read.
      !$acc enter data copyin(f_corner)
      call bt_work%enter_data()
      do step = 1, n_outer
         call barotropic_substep_nonlinear_interior(grid, metrics, bt_work, f_corner, &
                                                    n_inner, dt_inner)
         if (step <= n_early .or. step > n_outer - n_late) then
            !$acc update self(bt_work%bt_ubt)
            cur = max(abs(bt_work%bt_ubt(NGHOST + 2, jp)), abs(bt_work%bt_ubt(NGHOST + 3, jp)))
            if (step <= n_early) peak_early = max(peak_early, cur)
            if (step > n_outer - n_late) peak_late = max(peak_late, cur)
         end if
      end do
      call bt_work%exit_data()
      !$acc exit data delete(f_corner)
   end subroutine run_slosh

   subroutine test_slosh_bounded_with_chain(error)
      !! The two-cell thin-channel slosh WITH the chain: bt_rem built
      !! from av_rem damps the interior |u| by orders of magnitude over
      !! repeated outer-step blocks (fastloop_damping.py's variant (b)/
      !! (c): "fully damped").
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: f_corner(:, :)
      real(wp) :: peak_early, peak_late
      integer, parameter :: N_INNER = 4, N_OUTER = 48, N_EARLY = 5, N_LATE = 5
      real(wp), parameter :: DT_INNER = 100.0_wp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         call compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, N_INNER)
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1), source=0.0_wp)

         call run_slosh(bt_work, grid, metrics, f_corner, N_OUTER, N_INNER, DT_INNER, &
                        N_EARLY, N_LATE, peak_early, peak_late)
         ! Measured: peak_early=2.249E-03, peak_late=4.168E-04 (ratio 0.185)
         ! -- a damped oscillator's PE<->KE exchange keeps re-exciting u from
         ! the still-decaying eta each substep, so the envelope does NOT
         ! collapse like bare bt_rem**n_total (that identity is
         ! `test_column_spin_down_matches_av_rem`'s job, with no eta
         ! coupling); the qualitative contrast with the UNDAMPED run below
         ! (ratio ~1.06, bounded/neutral) is the thing being asserted here.
         call check(error, peak_early > 0.0_wp, "chain: the slosh never got going (peak_early == 0)")
         if (allocated(error)) exit checks
         call check(error, peak_late < 0.3_wp*peak_early, &
                    "chain: interior |u| did not decay substantially "// &
                    "over repeated outer steps")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_slosh_bounded_with_chain

   subroutine test_slosh_undamped_without_chain(error)
      !! The SAME slosh with bt_rem left at its init value of 1 (today's
      !! default, no chain): the forward-backward loop is neutrally
      !! stable for an undamped linear mode, so the interior |u| stays
      !! WITHIN A BOUND -- it must not collapse toward the chain-damped
      !! run's final amplitude (fastloop_damping.py's variant (a):
      !! "neutral, bounded", ratio ~1.4 over many more windows than used
      !! here).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: f_corner(:, :)
      real(wp) :: peak_early, peak_late
      integer, parameter :: N_INNER = 4, N_OUTER = 48, N_EARLY = 5, N_LATE = 5
      real(wp), parameter :: DT_INNER = 100.0_wp
      checks: block
         call build_channel(grid, ms, bt_work, metrics)
         ! bt_rem_u/v left at their init value of 1.0 (bt_rem_from_visc_rem
         ! never called) -- today's undamped default.
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1), source=0.0_wp)

         call run_slosh(bt_work, grid, metrics, f_corner, N_OUTER, N_INNER, DT_INNER, &
                        N_EARLY, N_LATE, peak_early, peak_late)
         ! Measured: peak_early=7.968E-03, peak_late=8.424E-03 (ratio 1.057,
         ! slightly growing -- fastloop_damping.py's own "neutral, bounded"
         ! finding, ratio ~1.4 over many more windows than this test runs).

         call check(error, peak_early > 0.0_wp, &
                    "no-chain: the slosh never got going (peak_early == 0)")
         if (allocated(error)) exit checks
         call check(error, ieee_is_finite(peak_late), &
                    "no-chain: the undamped slosh blew up to Inf/NaN")
         if (allocated(error)) exit checks
         call check(error, peak_late > 0.7_wp*peak_early, &
                    "no-chain: interior |u| collapsed -- bt_rem=1 should leave it "// &
                    "bounded/neutral, not decay it like the chain does")
      end block checks
      call bt_work%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_slosh_undamped_without_chain

end module test_ocean_bt_rem_from_visc_rem
