!! #178: a free-surface drawdown under `z_fixed` can thin the nominal bed
!! sliver `k_bot_u/v` (the STATIC index, filled once at configure from the
!! `η = 0` target) points at below `H_VANISHED` without moving the index
!! itself. Before the fix the bottom-drag kernels and the vdiff bed-BC row
!! stayed glued to that now-massless row; `zlevel_faces` then decoupled it
!! from the live layer above, so `visc_rem` sat at exactly `1.0` and NO
!! bottom drag reached the water (measured on a global 1-degree ice run: a
!! grounding shelf jet ran away and the run died by `ERROR STOP` at step
!! 8296 — see `tmp_local_artifacts/ice_hunt/RESULT_A.md` /
!! `RESULT_B2.md`). `kb_live` — the first row at/above the static index
!! whose face thickness is still live by the shared `rdb_blf_is_live`
!! criterion (`src/shared_module_utilities/rdb_bed_live_face.inc`) — fixes
!! this in every drag/BBL consumer.
!!
!! Each column here has ONE live layer surviving a drawdown that vanished
!! the static `k_bot_u`'s own row (and, for the "drawdown" case, a couple
!! of configure-time fillers below that too): `kb_live = nz`, which
!! degenerates the vdiff column solve to the documented single-layer case
!! (no interior coupling at all), giving a closed-form answer this test
!! can check to near machine precision instead of just a sign/trend.
!!
!! Cases:
!!   1. explicit path (`bbl_glue = .false.`): `ocean_bottom_drag_t`'s own
!!      implicit quadratic drag, called directly with no vdiff solve.
!!   2. BBL-glue path (`bbl_glue = .true.`, the production default): the
!!      MOM6 bed piston folded into the vdiff tridiagonal via
!!      `vdiff_apply_momentum`, read back through `visc_rem_u`.
!! Each case runs a DRAWDOWN column (static `k_bot_u` vanished, kb_live
!! one row up) against a CONTROL column (no drawdown, kb_live = k_bot_u
!! already) built from the exact same formula, which is the bit-identity
!! guard: the control's analytic answer is identical in form, only with
!! `kb_live = k_bot_u` baked in from the start.
module test_ocean_bottom_drag_live_kbot
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, &
                                    ocean_bottom_drag_compute_tendencies, &
                                    ocean_bottom_drag_apply_tendencies, &
                                    BDRAG_QUADRATIC
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_momentum
   implicit none
   private

   public :: collect_ocean_bottom_drag_live_kbot_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 4, NYP = 3
      !! Physical extent; total = +2*NGHOST.
   integer, parameter :: NFILL = 2
      !! Configure-time bed fillers below the static `k_bot` (drawdown
      !! case only) — exercises the "several vanished rows, not just one"
      !! shape, matching the real `z_fixed` staircase.
   real(wp), parameter :: H_FILL = 1.0e-4_wp
      !! `zstar_h_min`-style filler thickness; strictly below `H_VANISHED`.
   real(wp), parameter :: H_LIVE = 40.0_wp
   real(wp), parameter :: U0 = 0.3_wp
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: CD = 2.5e-3_wp
   real(wp), parameter :: DT = 600.0_wp
   real(wp), parameter :: PISTON = 1.0e-3_wp
      !! `&ocean_bdrag_nml` constant-piston test value (MOM6 `bbl_piston`);
      !! paired with `HBBL_THICK` below gives `kv_bbl = PISTON*HBBL_THICK`.
   real(wp), parameter :: HBBL_THICK = 5.0_wp
      !! `< 0.5*H_LIVE`, so the piston's `min(0.5*hvel, bbl_thick)` floor
      !! resolves to this value exactly — keeps the closed form simple.

contains

   subroutine collect_ocean_bottom_drag_live_kbot_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("explicit_quadratic_drag_on_kb_live", test_explicit), &
                  new_unittest("bbl_glue_piston_on_kb_live", test_bbl_glue) &
                  ]
   end subroutine collect_ocean_bottom_drag_live_kbot_tests

   ! ------------------------------------------------------------------
   ! Builder
   ! ------------------------------------------------------------------

   subroutine build_column(grid, ms, drawdown)
      !! `drawdown = .true.`: `NFILL` configure-time fillers, then the
      !! static `k_bot`'s OWN row ALSO thinned to `H_FILL` (the free-
      !! surface drawdown #178 is about), then ONE live layer at
      !! `kb_live = k_bot + 1`, carrying the uniform jet `U0`.
      !! `drawdown = .false.`: the control — `k_bot = 1`, immediately
      !! live, `kb_live = k_bot`, same jet.  Both: `v_face_y_layer = 0`
      !! throughout, so `|U| = U0` exactly (clean quadratic-law algebra).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      logical, intent(in) :: drawdown
      integer :: k, kb, nzl

      call grid%init(NXP, NYP, NGHOST, 1000.0_wp, 1000.0_wp)
      if (drawdown) then
         kb = NFILL + 1
         nzl = kb + 1
            !! Rows 1..NFILL config fillers, row kb (=NFILL+1) the
            !! drawdown-vanished static bed, row kb+1 = kb_live the sole
            !! live survivor -- single live layer, same as the control.
      else
         kb = 1
         nzl = 1
            !! No filler, no drawdown: `k_bot = 1` is already live, so
            !! `kb_live = k_bot = nz` trivially -- the SAME single-layer
            !! shape the drawdown case reduces to after the fix.
      end if
      ms%nz_ml = nzl
      call ms%init(grid)
      do k = 1, nzl
         if (drawdown .and. k <= kb) then
            ms%h_layer(:, :, k) = H_FILL
            ms%u_face_x_layer(:, :, k) = 0.0_wp
         else
            ms%h_layer(:, :, k) = H_LIVE
            ms%u_face_x_layer(:, :, k) = U0
         end if
      end do
      ms%v_face_y_layer = 0.0_wp
      ms%k_bot = kb
      ms%k_bot_u = kb
      ms%k_bot_v = kb
   end subroutine build_column

   subroutine map_in(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(ms)
      call ms%enter_data()
   end subroutine map_in

   subroutine map_out(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! 1. Explicit path (bbl_glue off)
   ! ------------------------------------------------------------------

   subroutine test_explicit(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_drawdown, u_control
      real(wp) :: expect_drawdown, expect_control

      call run_explicit(.true., u_drawdown)
      call run_explicit(.false., u_control)

      ! Closed form: implicit quadratic drag, one step, lambda = Cd*|U0|/H_LIVE,
      ! u_new = U0/(1 + dt*lambda) -- identical in both cases because kb_live
      ! sits on the SAME H_LIVE face thickness and the SAME U0 jet either way.
      expect_drawdown = U0/(1.0_wp + DT*CD*abs(U0)/H_LIVE)
      expect_control = expect_drawdown

      call check(error, abs(u_drawdown - expect_drawdown) < 1.0e-12_wp*abs(expect_drawdown), &
                 "drawdown: kb_live jet decelerates at the implicit quadratic rate")
      if (allocated(error)) return
      call check(error, abs(u_control - expect_control) < 1.0e-12_wp*abs(expect_control), &
                 "control (no drawdown): identical closed form at kb_live = k_bot")
      if (allocated(error)) return
      ! The #178 bug signature: on unfixed code the tendency is computed
      ! AND applied at the static (vanished) k_bot, so kb_live's jet is
      ! left at EXACTLY U0 -- the discriminator this assertion catches.
      call check(error, abs(u_drawdown - U0) > 1.0e-6_wp, &
                 "drawdown jet actually moved -- not left at the untouched U0 (#178 signature)")
   end subroutine test_explicit

   subroutine run_explicit(drawdown, u_live)
      logical, intent(in) :: drawdown
      real(wp), intent(out) :: u_live
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      integer :: nzl

      call build_column(grid, ms, drawdown)
      nzl = ms%nz_ml
      call bd%init(grid, nz_ml=nzl)
      bd%variant = BDRAG_QUADRATIC
      bd%c_drag = CD
      bd%implicit = .true.
      bd%zlevel_faces = .true.
         !! Exercise the #178 zlevel_faces branch of the shared kb_live
         !! search; the column is uniform in (i,j) so min() == arith mean
         !! here and cannot change the analytic answer above.

      call map_in(ms)
      !$acc enter data copyin(bd)
      call bd%enter_data()
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
      call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
      call bd%exit_data()
      !$acc exit data delete(bd)
      call map_out(ms)

      u_live = ms%u_face_x_layer(3, 2, nzl)
      call bd%destroy()
      call ms%destroy()
   end subroutine run_explicit

   ! ------------------------------------------------------------------
   ! 2. BBL-glue path (production default)
   ! ------------------------------------------------------------------

   subroutine test_bbl_glue(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_drawdown, u_control, rem_drawdown, rem_control
      real(wp) :: expect

      call run_bbl_glue(.true., u_drawdown, rem_drawdown)
      call run_bbl_glue(.false., u_control, rem_control)

      ! Single-live-layer column => kb_live == nz => the column solve
      ! degenerates to the documented no-interior-coupling case: a pure
      ! diagonal piston, u_new = U0/(1 + dt*kv_bbl/(H_LIVE*HBBL_THICK)).
      expect = U0/(1.0_wp + DT*(PISTON*HBBL_THICK)/(H_LIVE*HBBL_THICK))

      call check(error, abs(u_drawdown - expect) < 1.0e-10_wp*abs(expect), &
                 "drawdown: BBL-glue piston decelerates kb_live at the MOM6 closed form")
      if (allocated(error)) return
      call check(error, abs(u_control - expect) < 1.0e-10_wp*abs(expect), &
                 "control: identical closed form, guards bit-identity off drawdown")
      if (allocated(error)) return
      call check(error, rem_drawdown < 0.99_wp, &
                 "#178 signature: visc_rem at kb_live is NOT pinned to 1.0")
      if (allocated(error)) return
      call check(error, rem_control < 0.99_wp, &
                 "control: visc_rem also responds to the piston (same physics)")
      if (allocated(error)) return
      call check(error, abs(u_drawdown - U0) > 1.0e-6_wp, &
                 "drawdown jet actually moved -- not left at the untouched U0 (#178 signature)")
   end subroutine test_bbl_glue

   subroutine run_bbl_glue(drawdown, u_live, rem_live)
      logical, intent(in) :: drawdown
      real(wp), intent(out) :: u_live, rem_live
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      integer :: nzl, nx, ny
      real(wp), allocatable :: visc_rem_u(:, :, :), visc_rem_v(:, :, :)

      call build_column(grid, ms, drawdown)
      nzl = ms%nz_ml
      nx = grid%nx_total
      ny = grid%ny_total
      call vd%init(grid, nz_ml=nzl)
      vd%K_v_momentum = 0.0_wp
      vd%bbl_glue = .true.
      vd%hvel_mom6 = .true.
      vd%zlevel_faces = .true.
      vd%bbl_piston = PISTON
      vd%hbbl_visc = HBBL_THICK
      allocate (visc_rem_u(nx + 1, ny, nzl), source=1.0_wp)
      allocate (visc_rem_v(nx, ny + 1, nzl), source=1.0_wp)

      call map_in(ms)
      !$acc enter data copyin(vd)
      call vd%enter_data()
      !$acc enter data copyin(visc_rem_u, visc_rem_v)
      call vdiff_apply_momentum(grid, vd, ms, DT, &
                                visc_rem_u=visc_rem_u, visc_rem_v=visc_rem_v)
      !$acc update self(visc_rem_u, visc_rem_v)
      !$acc exit data delete(visc_rem_u, visc_rem_v)
      call vd%exit_data()
      !$acc exit data delete(vd)
      call map_out(ms)

      u_live = ms%u_face_x_layer(3, 2, nzl)
      rem_live = visc_rem_u(3, 2, nzl)
      call vd%destroy()
      call ms%destroy()
      deallocate (visc_rem_u, visc_rem_v)
   end subroutine run_bbl_glue

end module test_ocean_bottom_drag_live_kbot
