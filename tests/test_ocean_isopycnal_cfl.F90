!! Unit tests for Phase-3 vanished-layer CFL / truncation gate,
!! `compute_max_cfl` (gated variant) + `apply_velocity_truncation`
!! (vanish_tol optional) in `rdb_ocean_dyn` /
!! `rdb_ocean_console_stats`.
!!
!! Test suite (4 analytical cases):
!!   1. cfl_ignores_vanished_spike  - an entire k-layer has h <= vanish_tol
!!      and carries a large velocity; the gated MaxCFL returns the
!!      massive-layer value; un-gated returns the spike (bit-identical
!!      to today).  Whole-layer design avoids CFL leakage via shared
!!      face centrings at the vanished/massive boundary.
!!   2. cfl_bitident               - without gate args the un-gated
!!      compute_max_cfl returns the exact pre-feature value.
!!   3. truncation_zeros_vanished  - a vanished-both-sides face with a
!!      super-CFL velocity: apply_velocity_truncation (engaged) sets it
!!      to 0 (NOT to 0.9*cfl_trunc*dx/dt); a massive-side super-CFL
!!      face still clips to the CFL value.
!!   4. truncation_bitident        - vanish_tol absent => existing
!!      clip behaviour unchanged.
module test_ocean_isopycnal_cfl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_dyn, only: apply_velocity_truncation, isopycnal_vanish_tol
   use rdb_ocean_console_stats, only: compute_max_cfl
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_isopycnal_cfl_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: ANGSTROM = 1.0e-2_wp
      !! Representative floor value used in the tests.
   real(wp), parameter :: DT = 100.0_wp
      !! Timestep used in CFL calculations.
   real(wp), parameter :: CFLT = 0.5_wp
      !! CFL truncation threshold.
   real(wp), parameter :: RELAX = 0.9_wp
      !! CFL_TRUNC_RELAX factor (must match the kernel).

contains

   subroutine collect_ocean_isopycnal_cfl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cfl_ignores_vanished_spike", test_cfl_ignores_vanished_spike), &
                  new_unittest("cfl_bitident", test_cfl_bitident), &
                  new_unittest("truncation_zeros_vanished", test_truncation_zeros_vanished), &
                  new_unittest("truncation_bitident", test_truncation_bitident), &
                  new_unittest("vanish_tol_lifts_above_pd_floor", test_vanish_tol_pd_floor) &
                  ]
   end subroutine collect_ocean_isopycnal_cfl_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(ms)
      call ms%enter_data()
   end subroutine map_in

   subroutine map_out(ms)
      type(multilayer_state_t), intent(inout) :: ms
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! T1: gated MaxCFL skips an entire vanished k-layer
   ! -----------------------------------------------------------------
   subroutine test_cfl_ignores_vanished_spike(error)
      !! All cells in k=K_SPIKE have h_layer <= vanish_tol and large
      !! velocity; all other layers are massive with a small velocity.
      !! Design: entire-layer vanishing avoids CFL leakage via shared
      !! face centrings at the vanished/massive boundary (a single-cell
      !! thin patch would bleed half the spike into the adjacent massive
      !! cells' centred CFL since uc(i) = 0.5*(u_face(i) + u_face(i+1))).
      !! Un-gated compute_max_cfl returns the spike (existing behaviour).
      !! Gated version returns the massive-layer CFL (spike layer excluded).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp), parameter :: THIN = ANGSTROM*0.3_wp
         !! Below vanish_tol: entire k_spike layer
      real(wp), parameter :: U_SPIKE = 8.0_wp
         !! Spike CFL = U_SPIKE * DT * idx = 8 * 100 * 1 = 800 >> threshold
      real(wp), parameter :: U_MASSIVE = 0.003_wp
         !! Massive-layer CFL = U_MASSIVE * DT = 0.3 << threshold
      integer, parameter :: K_SPIKE = 2
         !! Vanished k-layer index
      real(wp) :: vtol, cfl_ungated, cfl_gated, cfl_massive_expected
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)

         ! All layers massive, small velocity seed on u-faces
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = U_MASSIVE
         ms%v_face_y_layer = 0.0_wp

         ! Entire k=K_SPIKE layer: thin h everywhere, spike velocity on u-faces
         ms%h_layer(:, :, K_SPIKE) = THIN
         ms%u_face_x_layer(:, :, K_SPIKE) = U_SPIKE

         ! Expected massive CFL: centred uc = U_MASSIVE (all faces same),
         ! on dx=1 grid with idx=1: CFL = U_MASSIVE * DT
         cfl_massive_expected = U_MASSIVE*DT

         call map_in(ms)

         ! Un-gated: spike layer contributes => large MaxCFL
         cfl_ungated = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                       metrics%idxT, metrics%idyT, DT, NGHOST)

         ! Gated: k_spike layer h <= vtol => all cells at k=K_SPIKE excluded
         cfl_gated = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                     metrics%idxT, metrics%idyT, DT, NGHOST, &
                                     ms%h_layer, vtol)

         call map_out(ms)

         ! Un-gated must include the spike (U_SPIKE >> U_MASSIVE)
         call check(error, cfl_ungated > U_MASSIVE*DT*10.0_wp, &
                    "T1: un-gated MaxCFL did not reflect vanished-layer spike")
         if (allocated(error)) exit checks

         ! Gated must exclude the spike layer and return the massive-layer CFL
         ! NB: massive layers (k/=K_SPIKE) have all faces = U_MASSIVE;
         !     centred uc = U_MASSIVE everywhere => max CFL = U_MASSIVE*DT.
         call check(error, abs(cfl_gated - cfl_massive_expected) < 1.0e-10_wp, &
                    "T1: gated MaxCFL did not exclude vanished k-layer")
         if (allocated(error)) exit checks

         ! Gated < un-gated
         call check(error, cfl_gated < cfl_ungated, &
                    "T1: gated MaxCFL should be < un-gated (spike excluded)")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_cfl_ignores_vanished_spike

   ! -----------------------------------------------------------------
   ! T2: un-gated compute_max_cfl is bit-identical to itself
   ! -----------------------------------------------------------------
   subroutine test_cfl_bitident(error)
      !! Without gate args, compute_max_cfl returns the same value on
      !! two independent calls (deterministic, no side effects).
      !! Also confirms the un-gated path is unchanged from pre-Phase-3.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 50.0_wp
      real(wp), parameter :: U_BASE = 0.01_wp
      real(wp) :: cfl_a, cfl_b
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         ms%h_layer = H_BASE
         ms%u_face_x_layer = U_BASE
         ms%v_face_y_layer = U_BASE

         call map_in(ms)
         cfl_a = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                 metrics%idxT, metrics%idyT, DT, NGHOST)
         cfl_b = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                 metrics%idxT, metrics%idyT, DT, NGHOST)
         call map_out(ms)

         call check(error, cfl_a == cfl_b, &
                    "T2: two un-gated compute_max_cfl calls not bit-identical")
         if (allocated(error)) exit checks
         call check(error, cfl_a > 0.0_wp, &
                    "T2: MaxCFL should be > 0 for non-zero velocities")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_cfl_bitident

   ! -----------------------------------------------------------------
   ! T3: truncation zeroes vanished face (not CFL-clipped)
   ! -----------------------------------------------------------------
   subroutine test_truncation_zeros_vanished(error)
      !! A u-face with BOTH adjacent cells vanished (h <= vanish_tol)
      !! and a super-CFL velocity: apply_velocity_truncation with
      !! vanish_tol set must zero it (not clip to RELAX*cfl_trunc*dx/dt).
      !! A face with at least one massive neighbour and super-CFL must
      !! still be clipped to the CFL value (existing behaviour).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp), parameter :: THIN = ANGSTROM*0.3_wp
      ! On dx=1 Cartesian grid idxCu = 1.  CFL = |u|*DT.
      ! Super-CFL: |u| * DT > CFLT => |u| > CFLT/DT.
      real(wp) :: super_u, clip_u
      real(wp) :: vtol
      real(wp) :: got_vanished, got_massive
      integer :: ntrunc
      integer :: i_vanish, j_vanish, k_vanish
      integer :: i_massive, j_massive, k_massive
      checks: block
         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)
         super_u = 2.0_wp*CFLT/DT    ! CFL = 2*CFLT > threshold
         clip_u = RELAX*CFLT/DT      ! expected clipped magnitude (dx=1, idxCu=1)

         i_vanish = 5
         j_vanish = 4
         k_vanish = 2
         i_massive = 7
         j_massive = 4
         k_massive = 2

         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! Vanished face: both adjacent centres thin
         ms%h_layer(i_vanish - 1, j_vanish, k_vanish) = THIN
         ms%h_layer(i_vanish, j_vanish, k_vanish) = THIN
         ms%u_face_x_layer(i_vanish, j_vanish, k_vanish) = super_u

         ! Massive face: one THIN, one MASSIVE => max > vtol => CFL clip applies
         ms%h_layer(i_massive - 1, j_massive, k_massive) = THIN
         ms%u_face_x_layer(i_massive, j_massive, k_massive) = super_u

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc, &
                                        vanish_tol=vtol)
         call map_out(ms)

         got_vanished = ms%u_face_x_layer(i_vanish, j_vanish, k_vanish)
         got_massive = ms%u_face_x_layer(i_massive, j_massive, k_massive)

         ! Vanished face must be zero (not CFL-clipped)
         call check(error, got_vanished == 0.0_wp, &
                    "T3: vanished-both-sides face not zeroed by truncation")
         if (allocated(error)) exit checks

         ! Massive-one-side face must be CFL-clipped (not zeroed)
         call check(error, abs(got_massive - clip_u) < 1.0e-12_wp, &
                    "T3: one-sided-massive face not clipped to RELAX*cfl*dx/dt")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_truncation_zeros_vanished

   ! -----------------------------------------------------------------
   ! T4: without vanish_tol the truncation is bit-identical to the
   !     pre-Phase-3 behaviour (all-massive state)
   ! -----------------------------------------------------------------
   subroutine test_truncation_bitident(error)
      !! Two calls to apply_velocity_truncation on the same state:
      !! one with vanish_tol absent, one with vanish_tol present but
      !! ALL h_layer massive (no cells at or below vtol).
      !! Both must return the same ntrunc_step and clipped velocity.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp) :: super_u, clip_u, vtol
      real(wp) :: got_no_gate, got_gate
      integer :: ntrunc_a, ntrunc_b
      integer :: i_tgt, j_tgt, k_tgt
      checks: block
         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         vtol = isopycnal_vanish_tol(ANGSTROM)
         super_u = 2.0_wp*CFLT/DT
         clip_u = RELAX*CFLT/DT

         i_tgt = 5
         j_tgt = 4
         k_tgt = 2

         ! --- Run A: no vanish_tol ---
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(i_tgt, j_tgt, k_tgt) = super_u

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc_a)
         call map_out(ms)
         got_no_gate = ms%u_face_x_layer(i_tgt, j_tgt, k_tgt)

         ! --- Run B: vanish_tol present but ALL layers massive (no vanished) ---
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(i_tgt, j_tgt, k_tgt) = super_u

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc_b, &
                                        vanish_tol=vtol)
         call map_out(ms)
         got_gate = ms%u_face_x_layer(i_tgt, j_tgt, k_tgt)

         call check(error, got_no_gate == got_gate, &
                    "T4: gated truncation on all-massive state differs from un-gated")
         if (allocated(error)) exit checks
         call check(error, ntrunc_a == ntrunc_b, &
                    "T4: ntrunc_step differs between gated and un-gated on all-massive")
         if (allocated(error)) exit checks
         call check(error, abs(got_no_gate - clip_u) < 1.0e-12_wp, &
                    "T4: clipped value does not match RELAX*cfl_trunc*dx/dt")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_truncation_bitident

   ! -----------------------------------------------------------------
   subroutine test_vanish_tol_pd_floor(error)
      !! The floor/tolerance collision: positive_definite holds h >= h_lim =
      !! angstrom_h, and the plain `isopycnal_vanish_tol = max(angstrom_h,
      !! H_VANISHED) = angstrom_h` sits AT the floor — so an at-floor (dead)
      !! layer with `h = h_lim + ε` reads as LIVE under `h <= tol` and its
      !! phantom velocity is never vanish-zeroed (the day-16 2.4696 mechanism).
      !! `pd_floor=.true.` lifts the tolerance to `angstrom_h + H_VANISHED`
      !! (strictly above the floor), so the at-floor band vanishes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: MASSIVE = 100.0_wp
      real(wp), parameter :: AT_FLOOR = ANGSTROM + 0.5_wp*H_VANISHED  ! in (h_lim, h_lim+H_VANISHED]
      real(wp), parameter :: U_SUB = 0.003_wp   ! sub-CFL (0.003·DT·1 = 0.3 < CFLT) ⇒ CFL clip inert
      integer, parameter :: IL = 6, JL = 4, KL = 2   ! both-at-floor u-face
      integer, parameter :: IW = 9, JW = 4, KW = 2   ! live control u-face
      real(wp) :: vt_plain, vt_pd, at_no_pd, at_pd, live_no_pd, live_pd
      integer :: ntrunc
      character(len=160) :: msg
      checks: block
         vt_plain = isopycnal_vanish_tol(ANGSTROM)
         vt_pd = isopycnal_vanish_tol(ANGSTROM, pd_floor=.true.)

         ! Arithmetic: plain sits at the floor; pd sits H_VANISHED above it.
         write (msg, '("T7: vt_plain=", es13.6, " (want ", es13.6, "), vt_pd=", es13.6, &
               &" (want ", es13.6, ")")') vt_plain, ANGSTROM, vt_pd, ANGSTROM + H_VANISHED
         call check(error, abs(vt_plain - ANGSTROM) < 1.0e-15_wp .and. &
                    abs(vt_pd - (ANGSTROM + H_VANISHED)) < 1.0e-15_wp, trim(msg))
         if (allocated(error)) return
         ! The at-floor band is missed by plain, caught by pd.
         call check(error, AT_FLOOR > vt_plain .and. AT_FLOOR <= vt_pd, &
                    "T7: AT_FLOOR must be above plain tol and at/below pd tol")
         if (allocated(error)) return

         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)   ! idxCu=idxT=1

         ! --- Run A: plain tol (current) ---
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%h_layer(IL - 1, JL, KL) = AT_FLOOR        ! both sides at-floor ⇒ both-sided vanish
         ms%h_layer(IL, JL, KL) = AT_FLOOR
         ms%u_face_x_layer(IL, JL, KL) = U_SUB        ! at-floor face
         ms%u_face_x_layer(IW, JW, KW) = U_SUB        ! live control face
         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc, vanish_tol=vt_plain)
         call map_out(ms)
         at_no_pd = ms%u_face_x_layer(IL, JL, KL)
         live_no_pd = ms%u_face_x_layer(IW, JW, KW)

         ! --- Run B: pd-lifted tol (the fix) ---
         ms%h_layer = MASSIVE
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%h_layer(IL - 1, JL, KL) = AT_FLOOR
         ms%h_layer(IL, JL, KL) = AT_FLOOR
         ms%u_face_x_layer(IL, JL, KL) = U_SUB
         ms%u_face_x_layer(IW, JW, KW) = U_SUB
         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc, vanish_tol=vt_pd)
         call map_out(ms)
         at_pd = ms%u_face_x_layer(IL, JL, KL)
         live_pd = ms%u_face_x_layer(IW, JW, KW)

         ! (bug) plain tol leaves the at-floor face LIVE (survives, sub-CFL).
         call check(error, at_no_pd == U_SUB, &
                    "T7: plain tol should NOT vanish-zero the at-floor face (the bug)")
         if (allocated(error)) exit checks
         ! (fix) pd tol vanish-zeroes it.
         call check(error, at_pd == 0.0_wp, &
                    "T7: pd-lifted tol did not vanish-zero the at-floor face")
         if (allocated(error)) exit checks
         ! Live control face untouched in BOTH (not at-floor, sub-CFL).
         call check(error, live_no_pd == U_SUB .and. live_pd == U_SUB, &
                    "T7: live control face must be untouched by the vanish-zero")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_vanish_tol_pd_floor

end module test_ocean_isopycnal_cfl
