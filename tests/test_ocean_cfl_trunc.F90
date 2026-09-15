!! Unit tests for the ocean advective-CFL velocity truncation (E7),
!! `apply_velocity_truncation` in `rdb_ocean_dyn`.
!!
!! The kernel runs two stages at the end of each outer dyn step:
!!   1. CFL clip (gated `cfl_trunc > 0`): any face whose local
!!      advective CFL `|u|·dt·idx` exceeds `cfl_trunc` is reset to
!!      `sign(0.9·cfl_trunc/(dt·idx), u_old)` and counted.
!!   2. maxvel cap (gated `maxvel > 0`): the absolute physical
!!      backstop, applied second.
!!
!! Cases (analytical, high-leverage):
!!   * T1 clip correctness — a super-CFL +u face and a super-CFL −u
!!     face clip to ±0.9·cfl_trunc·dx/dt (sign preserved); a sub-CFL
!!     face is byte-identical.  Same for a v-face with idyCv.
!!   * T2 count — seed exactly N faces (mix u + v) above threshold,
!!     assert ntrunc_step == N.
!!   * T4 maxvel interaction — both knobs set.  (a) CFL clip lands a
!!     face above maxvel ⇒ final = maxvel (maxvel runs second, wins).
!!     (b) maxvel large, cfl_trunc binds ⇒ final = 0.9·cfl_trunc·dx/dt.
!!   * T5 NaN-catch (step −1; pdc 8c2fd674/366dd5f4, adapted 25d6f911)
!!     — NaN/±Inf faces are zeroed + counted in `n_nanzero` BEFORE the
!!     clip and the maxvel cap, with maxvel ON: proves a non-finite
!!     velocity can never launder to ±maxvel.
!!   * T6 cell-metric clip (pdc 6991988f, adapted 74966912) — the
!!     θ-edge escape: a face benign on idxCu but super-threshold on the
!!     adjacent idxT passes the face-metric clip (mode a, escape
!!     documented) and is clipped on max(idxT−,idxT+) under
!!     clip_cell_metric=.true. (mode b, the panic basis).
!!
!! The double_gyre byte-gate (T3, bit-identity at the default
!! cfl_trunc = 0) is a run-time regression, not a unit test — it is
!! exercised by the build/validate harness, not here.
module test_ocean_cfl_trunc
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                                               ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_dyn, only: apply_velocity_truncation
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_cfl_trunc_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: RELAX = 0.9_wp
      !! Must match CFL_TRUNC_RELAX inside apply_velocity_truncation.

contains

   subroutine collect_ocean_cfl_trunc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("clip_correctness", test_clip_correctness), &
                  new_unittest("truncation_count", test_count), &
                  new_unittest("maxvel_interaction", test_maxvel_interaction), &
                  new_unittest("nan_catch", test_nan_catch), &
                  new_unittest("cell_metric_clip", test_cell_metric_clip) &
                  ]
   end subroutine collect_ocean_cfl_trunc_tests

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
   ! T1: clip correctness
   ! -----------------------------------------------------------------
   subroutine test_clip_correctness(error)
      !! Cartesian dx=dy=1 ⇒ idxCu = idyCv = 1.  Seed one +u super-CFL
      !! face, one −u super-CFL face, one sub-CFL u face, and one +v
      !! super-CFL face.  After the call the super-CFL faces must equal
      !! ±RELAX·cfl_trunc/(dt·idx) (sign preserved); sub-CFL untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: CFLT = 0.5_wp
      real(wp), parameter :: MAXVEL = 0.0_wp   ! maxvel off
      real(wp) :: idx, super_pos, super_neg, sub_val, clip_mag
      real(wp) :: got_pos, got_neg, got_sub, got_v
      integer :: ntrunc
      checks: block
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         idx = 1.0_wp                              ! 1/dx with dx=1
         clip_mag = RELAX*CFLT/(DT*idx)            ! expected clipped magnitude
         super_pos = 2.0_wp*CFLT/(DT*idx)          ! CFL = 2*cfl_trunc > threshold
         super_neg = -super_pos
         sub_val = 0.5_wp*CFLT/(DT*idx)            ! CFL = 0.5*cfl_trunc < threshold

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(5, 4, 2) = super_pos
         ms%u_face_x_layer(6, 4, 2) = super_neg
         ms%u_face_x_layer(7, 4, 2) = sub_val
         ms%v_face_y_layer(5, 5, 3) = super_pos

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, MAXVEL, ntrunc)
         call map_out(ms)

         got_pos = ms%u_face_x_layer(5, 4, 2)
         got_neg = ms%u_face_x_layer(6, 4, 2)
         got_sub = ms%u_face_x_layer(7, 4, 2)
         got_v = ms%v_face_y_layer(5, 5, 3)

         call check(error, abs(got_pos - clip_mag) < 1.0e-12_wp, &
                    "T1: +u super-CFL face not clipped to +RELAX*cfl*dx/dt")
         if (allocated(error)) exit checks
         call check(error, abs(got_neg + clip_mag) < 1.0e-12_wp, &
                    "T1: -u super-CFL face not clipped to -RELAX*cfl*dx/dt")
         if (allocated(error)) exit checks
         call check(error, got_sub == sub_val, &
                    "T1: sub-CFL u face was modified (must be bitwise unchanged)")
         if (allocated(error)) exit checks
         call check(error, abs(got_v - clip_mag) < 1.0e-12_wp, &
                    "T1: +v super-CFL face not clipped to +RELAX*cfl*dy/dt")
         if (allocated(error)) exit checks
         call check(error, ntrunc == 3, "T1: expected 3 truncations (2 u + 1 v)")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_clip_correctness

   ! -----------------------------------------------------------------
   ! T2: truncation count
   ! -----------------------------------------------------------------
   subroutine test_count(error)
      !! Seed exactly N faces (mix u + v) above threshold and assert
      !! ntrunc_step == N.  Faces left at zero (sub-CFL) must not count.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: CFLT = 0.5_wp
      real(wp) :: super
      integer :: ntrunc, n_u, n_v
      checks: block
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         super = 3.0_wp*CFLT/(DT*1.0_wp)   ! super-CFL magnitude

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! 4 u-faces over threshold
         ms%u_face_x_layer(3, 3, 1) = super
         ms%u_face_x_layer(4, 3, 1) = -super
         ms%u_face_x_layer(5, 6, 2) = super
         ms%u_face_x_layer(6, 7, 3) = -super
         ! 3 v-faces over threshold
         ms%v_face_y_layer(3, 3, 1) = super
         ms%v_face_y_layer(4, 4, 2) = -super
         ms%v_face_y_layer(5, 5, 3) = super
         n_u = 4
         n_v = 3

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc)
         call map_out(ms)

         call check(error, ntrunc == n_u + n_v, &
                    "T2: ntrunc_step does not equal seeded over-threshold count")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_count

   ! -----------------------------------------------------------------
   ! T4: maxvel interaction (documented ordering)
   ! -----------------------------------------------------------------
   subroutine test_maxvel_interaction(error)
      !! (a) maxvel binds: CFL clip lands a face at RELAX*cfl*dx/dt,
      !!     then maxvel (smaller) caps it ⇒ final = maxvel.
      !! (b) cfl binds: maxvel is large, so CFL clip wins ⇒ final =
      !!     RELAX*cfl*dx/dt.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: CFLT = 0.5_wp
      real(wp) :: clip_mag, super, maxvel_small, maxvel_big, got
      integer :: ntrunc
      checks: block
         clip_mag = RELAX*CFLT/(DT*1.0_wp)   ! = 0.0045
         super = 5.0_wp*CFLT/(DT*1.0_wp)

         ! ---- (a) maxvel runs second and wins ----
         maxvel_small = 0.5_wp*clip_mag
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(5, 4, 2) = super
         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, maxvel_small, ntrunc)
         call map_out(ms)
         got = ms%u_face_x_layer(5, 4, 2)
         call check(error, abs(got - maxvel_small) < 1.0e-12_wp, &
                    "T4a: maxvel (smaller) must win after the CFL clip")
         call destroy_cartesian_metrics(metrics)
         call ms%destroy()
         if (allocated(error)) exit checks

         ! ---- (b) cfl_trunc binds, maxvel does not ----
         maxvel_big = 100.0_wp
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(5, 4, 2) = super
         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, maxvel_big, ntrunc)
         call map_out(ms)
         got = ms%u_face_x_layer(5, 4, 2)
         call check(error, abs(got - clip_mag) < 1.0e-12_wp, &
                    "T4b: CFL clip must bind when maxvel is large")
         call destroy_cartesian_metrics(metrics)
         call ms%destroy()
      end block checks
   end subroutine test_maxvel_interaction

   ! -----------------------------------------------------------------
   ! T5: NaN-catch — non-finite faces zeroed + counted, never laundered
   ! -----------------------------------------------------------------
   subroutine test_nan_catch(error)
      !! Step −1 of apply_velocity_truncation.  Comparisons with NaN are
      !! FALSE, so the CFL clip skips a NaN face, and nvfortran -fast
      !! lowers the maxvel if/else clamp to a NaN-blind min/max — a NaN
      !! would come out as ±maxvel and transport mass at maxvel forever.
      !! The catch must zero + count every non-finite face BEFORE the
      !! clip and the cap.  maxvel is deliberately ON here: the check
      !! `face == 0` (not ±maxvel) is the anti-laundering assertion.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: CFLT = 0.5_wp
      real(wp), parameter :: MAXVEL = 1.0_wp   ! ON — the laundering target
      real(wp) :: nan_val, inf_val, sub_val, super, clip_mag
      integer :: ntrunc, nnan
      checks: block
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         nan_val = ieee_value(1.0_wp, ieee_quiet_nan)
         inf_val = ieee_value(1.0_wp, ieee_positive_inf)
         clip_mag = RELAX*CFLT/(DT*1.0_wp)
         super = 2.0_wp*CFLT/(DT*1.0_wp)
         sub_val = 0.5_wp*CFLT/(DT*1.0_wp)

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(5, 4, 2) = nan_val
         ms%u_face_x_layer(6, 5, 1) = inf_val
         ms%v_face_y_layer(4, 4, 3) = -inf_val
         ms%u_face_x_layer(8, 6, 2) = super     ! real clip must still fire
         ms%u_face_x_layer(9, 7, 3) = sub_val   ! finite bystander

         call map_in(ms)
         call apply_velocity_truncation(ms, metrics, DT, CFLT, MAXVEL, ntrunc, &
                                        n_nanzero=nnan)
         call map_out(ms)

         call check(error, nnan == 3, &
                    "T5: n_nanzero must count all 3 non-finite faces")
         if (allocated(error)) exit checks
         call check(error, ms%u_face_x_layer(5, 4, 2) == 0.0_wp, &
                    "T5: NaN u face must be zeroed, not laundered to +-maxvel")
         if (allocated(error)) exit checks
         call check(error, ms%u_face_x_layer(6, 5, 1) == 0.0_wp, &
                    "T5: +Inf u face must be zeroed")
         if (allocated(error)) exit checks
         call check(error, ms%v_face_y_layer(4, 4, 3) == 0.0_wp, &
                    "T5: -Inf v face must be zeroed")
         if (allocated(error)) exit checks
         call check(error, abs(ms%u_face_x_layer(8, 6, 2) - clip_mag) < 1.0e-12_wp, &
                    "T5: real super-CFL face must still be clipped")
         if (allocated(error)) exit checks
         call check(error, ms%u_face_x_layer(9, 7, 3) == sub_val, &
                    "T5: finite sub-CFL face must be bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, ntrunc == 1, &
                    "T5: only the real super-CFL face counts as a truncation")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_nan_catch

   ! -----------------------------------------------------------------
   ! T6: cell-metric clip — the theta-edge escape is closed
   ! -----------------------------------------------------------------
   subroutine test_cell_metric_clip(error)
      !! The console panic / compute_max_cfl measure CFL on the CELL
      !! metric idxT; the clip historically bounded on the FACE metric
      !! idxCu.  Where the two diverge (grounding theta-edge, anomalous
      !! cell) a face can satisfy `|u|·dt·idxCu ≤ cfl_trunc` while the
      !! panic sees `|u|·dt·idxT > cfl_trunc` — truncations climb but
      !! the panic still fires (measured: MaxCFL 0.689 past a 0.5
      !! ceiling at 1024²/dt=800).  Uniform Cartesian metrics with one
      !! idxT (and one idyT) spiked ×4 by hand reproduce the divergence:
      !! (a) face-metric mode leaves the face untouched (the escape),
      !! (b) clip_cell_metric=.true. clips on max(idxT−,idxT+) — and a
      !! uniform-metric face already clipped in (a) stays sub-threshold
      !! in (b) (the 0.9 relax + metric equality ⇒ bit-identity there).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: CFLT = 0.5_wp
      real(wp), parameter :: SPIKE = 4.0_wp    ! anomalous idxT/idyT factor
      real(wp) :: u0, super, clip_uniform, clip_spiked
      integer :: ntrunc
      checks: block
         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)

         ! Face CFL = 0.8*CFLT (benign); cell CFL = SPIKE*0.8*CFLT = 1.6 (super)
         u0 = 0.8_wp*CFLT/DT
         super = 2.0_wp*CFLT/DT
         clip_uniform = RELAX*CFLT/(DT*1.0_wp)
         clip_spiked = RELAX*CFLT/(DT*SPIKE)

         ! Spike the cell metric next to the target faces (host), then
         ! push to the device (mapped arrays get no implicit updates).
         metrics%idxT(5, 4) = SPIKE          ! u-face (5,4): max(idxT(4,4), idxT(5,4)) = SPIKE
         metrics%idyT(6, 6) = SPIKE          ! v-face (6,6): max(idyT(6,5), idyT(6,6)) = SPIKE
         !$acc update device(metrics%idxT, metrics%idyT)

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(5, 4, 2) = u0      ! theta-edge u face
         ms%v_face_y_layer(6, 6, 1) = u0      ! theta-edge v face
         ms%u_face_x_layer(9, 7, 1) = super   ! uniform-metric control face

         call map_in(ms)

         ! (a) face-metric mode: the escape — theta-edge faces untouched,
         !     only the uniform control face clips.
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc)
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         call check(error, ms%u_face_x_layer(5, 4, 2) == u0, &
                    "T6a: face-metric mode must NOT clip the theta-edge u face (the escape)")
         if (allocated(error)) then
            call map_out(ms)
            exit checks
         end if
         call check(error, ntrunc == 1, &
                    "T6a: only the uniform control face clips under the face metric")
         if (allocated(error)) then
            call map_out(ms)
            exit checks
         end if

         ! (b) cell-metric mode: both theta-edge faces clip on the spiked
         !     idxT/idyT; the control face (already at 0.9 relax) stays put.
         call apply_velocity_truncation(ms, metrics, DT, CFLT, 0.0_wp, ntrunc, &
                                        clip_cell_metric=.true.)
         call map_out(ms)

         call check(error, abs(ms%u_face_x_layer(5, 4, 2) - clip_spiked) < 1.0e-12_wp, &
                    "T6b: theta-edge u face must clip to RELAX*cfl/(dt*idxT_spike)")
         if (allocated(error)) exit checks
         call check(error, abs(ms%v_face_y_layer(6, 6, 1) - clip_spiked) < 1.0e-12_wp, &
                    "T6b: theta-edge v face must clip to RELAX*cfl/(dt*idyT_spike)")
         if (allocated(error)) exit checks
         call check(error, abs(ms%u_face_x_layer(9, 7, 1) - clip_uniform) < 1.0e-12_wp, &
                    "T6b: uniform-metric face clipped in (a) must stay sub-threshold in (b)")
         if (allocated(error)) exit checks
         call check(error, ntrunc == 2, &
                    "T6b: exactly the two theta-edge faces clip under the cell metric")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine test_cell_metric_clip

end module test_ocean_cfl_trunc
