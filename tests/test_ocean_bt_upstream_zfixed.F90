!! `&ocean_bt_nml upstream_h_face` (and `substep_drag` / `wave_drag`)
!! under the partial-step z-level face closure (`&vcoord_nml
!! zfixed_closed_faces`) — the barotropic face depth must be the OPEN
!! column's.
!!
!! ### The defect
!!
!! Under closed faces the barotropic substep transports on
!! `h_face_up·ubt·dy_cu_bt`, where `dy_cu_bt = dy_cu·φ_c` carries the
!! CENTRED open fraction `φ_c = Σ h_c·open / Σ h_c`, while `ubt` is the
!! OPEN-column upstream mean and the renormaliser hands `uhbt` to the
!! open layers only.  `compute_h_face_upstream` used to sum the FULL
!! upstream column `H_up`, so the fast loop transported `H_up·φ_c` per
!! unit `ubt` where the layers carry `Σ_k h_up,k·open_k`.  At a staircase
!! face between a deep column `H_D` and a shallow one `H_S` that ratio is
!! `2·H_D/(H_D+H_S)` for flow off the deep side and `2·H_S/(H_D+H_S)` off
!! the shallow side — O(1), not round-off (`python_prototypes/
!! bt_upstream_zfixed/`).  On the 1-degree Southern Ocean it took the
!! barotropic velocity to the `maxvel` clamp and to NaN at step 309.
!!
!! ### The assertions
!!
!! **`upstream_face_depth_is_open_column`** — the DIRECT statement, on a
!! hand-built staircase face with `dy_cu_bt` from the production
!! `closed_faces_update_bt_widths`: the fast-loop face depth per unit
!! width, `h_face_up·dy_cu_bt/dy_cu`, must equal the open upstream column
!! `Σ_k h_up,k·open_k` for BOTH flow directions, to a derived round-off
!! bound.  It also asserts the full-column product DIFFERS here (by the
!! factor above), so the case cannot pass vacuously.  Fails on the
!! pre-port code by `4/3` and `2/3`.
!!
!! **`bt_rem_is_open_column`** — `compute_bt_rem` (`substep_drag`) and
!! `compute_bt_rem_wave_drag` (`wave_drag`) damp on the OPEN-column face
!! depth `Σ_k h_face·open` — the closed form, to round-off — and the
!! full-column answer is shown to differ on the same face.
!!
!! **`resting_stratified_staircase_upstream`** — a stably stratified
!! staircase basin at rest, `upstream_h_face` on, `pred_corr`: it stays
!! at rest.  The bed is aligned with the nominal layers so there is no
!! partial cell and the open-face PGF is exactly zero; the bound is
!! round-off on the velocity.
!!
!! **`staircase_seiche_upstream`** — the unforced, undamped staircase
!! seiche with `upstream_h_face` on, the production per-outer-step
!! `ocean_porous_refresh` of the narrowed BT widths, and `pred_corr`:
!! closed faces carry exactly zero velocity, `KE+PE` does not grow, and
!! the period is the EXACT stepped-basin one (`stepped_basin_period`).
!!
!! **`rotating_staircase_upstream`** — the same basin at `f = 1e-3`:
!! closed faces exactly zero, no energy growth.
!!
!! What discriminates, and what does not.  The two DIRECT statements fail
!! on the pre-port code (face depth `266.67` against `200`, the `4/3` of
!! the 400 | 200 m face; `bt_rem` at the full-column value).  The three
!! integration legs are regression gates for the ported path and do NOT
!! discriminate at this size: MEASURED (gfortran 15.1) `KE+PE` ratios
!! 0.896 / 0.947 (seiche / rotating) ported against 0.896 / 0.942 with
!! the full-column sum restored — an 80 km, 4-layer, unforced basin
!! damps the inconsistency faster than it amplifies.  The amplification
!! needs the production envelope; it is measured on the 1-degree Southern
!! Ocean (`zfixed_audit/bt_upstream_h_face`: day-1 MaxCFL 0.196 and NaN at
!! step 309 before, day-10 `En 5.709E-04` / MaxCFL `0.02935` after,
!! against the knob-off `5.493E-04` / `0.02326`).
module test_ocean_bt_upstream_zfixed
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_LITE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, ocean_porous_refresh, &
                            SPLIT_SCHEME_PRED_CORR
   use rdb_barotropic_coupling, only: compute_h_face_upstream, compute_bt_rem, &
                                      compute_bt_rem_wave_drag
   use rdb_ocean_porous, only: closed_faces_update_bt_widths
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_Z_FIXED, &
                               ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks
   use testdrive, only: error_type, check, new_unittest, unittest_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_bt_upstream_zfixed_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 40
   integer, parameter :: NYP = 4
   integer, parameter :: NZ = 4
   integer, parameter :: N_INNER = 20
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: GRAV = 9.80665_wp
   real(wp), parameter :: H_NOM = 100.0_wp
   real(wp), parameter :: H_MIN = 1.0e-4_wp
   real(wp), parameter :: H_DEEP = 400.0_wp
      !! Western column: 4 full nominal layers.
   real(wp), parameter :: H_SHELF = 100.0_wp
      !! Eastern column: THREE nominal layers shallower, so only the
      !! surface layer is live there and every deeper layer closes at the
      !! faces touching the shelf.  The full-column vs open-column face
      !! depth then differs at the step by `2·H_D/(H_D+H_S) = 1.6` off the
      !! deep side and `0.4` off the shelf — the ratios of the prototype.
   real(wp), parameter :: ETA0 = 0.05_wp
   real(wp), parameter :: PI_L = 3.14159265358979324_wp
   real(wp), parameter :: DT = 60.0_wp
   integer, parameter :: N_PERIODS = 22

   real(wp), parameter :: ENERGY_GROWTH_BAR = 1.30_wp
      !! Bar for `(KE+PE)_end/(KE+PE)_0` over `N_PERIODS` seiche periods,
      !! the same one-sided bar `test_ocean_zfixed_bt_seiche` uses for
      !! its centred-face `pred_corr` leg.  Unforced and undamped ⇒ any
      !! growth is manufactured.
   real(wp), parameter :: PERIOD_TOL = 0.03_wp
      !! Fractional tolerance on the seiche period against the EXACT
      !! gravest-mode period of the two-step basin (`stepped_basin_period`;
      !! the width-weighted `2L/√(g·H_eff)` estimate the sibling test uses
      !! is far off for a 4:1 step).  MEASURED 1.01 % apart on gfortran
      !! (4243 s against 4201 s; the residual is the C-grid dispersion and
      !! the discrete step), so 3 % is ~3x of headroom.  A barotropic mode that saw the FULL
      !! upstream column would be off by the depth error, not by this.
   real(wp), parameter :: REST_U_BAR = 1.0e-12_wp
      !! Bound on `max|u|` (m/s) for the resting stratified staircase.
      !! The open-face PGF is identically zero (layers aligned with the
      !! bed, identical `ρ` per layer in every column), so the only
      !! velocity is round-off in the `η`-gradient of a flat free surface:
      !! `g·Δη·dt/dx` with `Δη ~ 1e-16·H` is `~1e-15` per step.
   real(wp), parameter :: REST_PERIODS = 3.0_wp
      !! Length of the resting leg (seiche periods, ~200 steps): it has no
      !! mode to resolve, only round-off that must stay round-off.
   real(wp), parameter :: ROT_PERIODS = 8.0_wp
      !! Length of the rotating leg (~5 inertial periods at `F_ROT`).
   real(wp), parameter :: F_ROT = 1.0e-3_wp
      !! Coriolis parameter of the rotating leg (1/s).  Ten times a
      !! mid-latitude `f` so the barotropic deformation radius
      !! `√(gH)/f ≈ 30-60 km` fits inside the 80 km basin and rotation
      !! actually shapes the mode (at `1e-4` the leg is the f = 0 seiche
      !! to three figures).
   real(wp), parameter :: DRHO = 2.0_wp
      !! Bed-to-surface density contrast of the resting case (kg/m³).

contains

   subroutine collect_ocean_bt_upstream_zfixed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("upstream_face_depth_is_open_column", test_face_depth), &
                  new_unittest("bt_rem_is_open_column", test_bt_rem_open), &
                  new_unittest("resting_stratified_staircase_upstream", test_resting), &
                  new_unittest("staircase_seiche_upstream", test_seiche), &
                  new_unittest("rotating_staircase_upstream", test_rotating) &
                  ]
   end subroutine collect_ocean_bt_upstream_zfixed_tests

   ! ---------------------------------------------------------------
   !  1. The face-depth identity, on one hand-built staircase face.
   ! ---------------------------------------------------------------

   subroutine build_step_face(grid, metrics, ms, bt_work, h_w, h_e, u_k)
      !! Two columns (west `i=1`, east `i=2`), one interior u-face `I=2`,
      !! `nghost = 0`.  Masks from the production builder on `h` itself
      !! (a layer live on both sides is open), `dy_cu_bt`/`dx_cv_bt` from
      !! the production `closed_faces_update_bt_widths`.  Everything is
      !! device-mapped and the host inputs are pushed before any kernel
      !! runs (`mem:separate` contract).
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(inout) :: bt_work
      real(wp), intent(in) :: h_w(NZ), h_e(NZ), u_k(NZ)
      integer :: k

      call grid%init(2, 1, 0, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      do k = 1, NZ
         ms%h_layer(1, 1, k) = h_w(k)
         ms%h_layer(2, 1, k) = h_e(k)
         ms%u_face_x_layer(:, 1, k) = u_k(k)
      end do
      call bt_work%init(grid, nz_ml=NZ)
      bt_work%use_upstream_h_face = .true.
      allocate (bt_work%h_face_up_x(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
      allocate (bt_work%h_face_up_y(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
      metrics%use_closed_faces = .true.

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call bt_work%enter_data()
      ! Both builders are device kernels on device-present arrays.
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, ms%h_layer, &
                                          grid%nx_total, grid%ny_total, NZ, H_VANISHED)
      call closed_faces_update_bt_widths(grid%nx_total, grid%ny_total, NZ, .false., &
                                         metrics%dy_cu, metrics%dx_cv, ms%h_layer, &
                                         metrics%open_u, metrics%open_v, &
                                         metrics%open_u, metrics%open_v, &
                                         metrics%dy_cu_bt, metrics%dx_cv_bt)
      !$acc update self(metrics%open_u, metrics%open_v, metrics%dy_cu_bt, metrics%dx_cv_bt)
   end subroutine build_step_face

   subroutine teardown_step_face(metrics, ms, bt_work)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(inout) :: bt_work
      call bt_work%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call bt_work%destroy()
      call ms%destroy()
   end subroutine teardown_step_face

   subroutine one_direction(error, u_sign, label)
      !! One flow direction across the 400 m | 200 m step.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: u_sign
      character(len=*), intent(in) :: label
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp) :: h_w(NZ), h_e(NZ), h_up(NZ), u_k(NZ)
      real(wp) :: d_open, d_full, d_bt, tol
      integer :: k
      character(len=400) :: msg

      ! West = deep (4 live layers), east = shelf (bed + next layer are
      ! fillers).  Bottom-up: k = 1 is the bed.
      h_w = [H_NOM, H_NOM, H_NOM, H_NOM]
      h_e = [H_MIN, H_MIN, H_NOM, H_NOM]
      u_k = u_sign*[0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp]
      call build_step_face(grid, metrics, ms, bt_work, h_w, h_e, u_k)

      call compute_h_face_upstream(grid, bt_work, ms, metrics)
      !$acc update self(bt_work%h_face_up_x)

      if (u_sign > 0.0_wp) then
         h_up = h_w
      else
         h_up = h_e
      end if
      d_open = 0.0_wp
      d_full = 0.0_wp
      do k = 1, NZ
         d_open = d_open + h_up(k)*metrics%open_u(2, 1, k)
         d_full = d_full + h_up(k)
      end do
      ! What the fast loop multiplies `ubt` by, per unit of the UN-narrowed
      ! width: `h_face_up·dy_cu_bt/dy_cu`.
      d_bt = bt_work%h_face_up_x(2, 1)*metrics%dy_cu_bt(2, 1)/metrics%dy_cu(2, 1)

      ! Anti-vacuity: the face is a partial one (2 of 4 layers open) and the
      ! FULL-column product differs from the open depth by O(1).
      write (msg, '(A,": the step face must be PARTIAL (2 of 4 layers open), got ", &
            &F4.1," open")') label, sum(metrics%open_u(2, 1, :))
      call check(error, sum(metrics%open_u(2, 1, :)) == 2.0_wp, trim(msg))
      if (allocated(error)) then
         call teardown_step_face(metrics, ms, bt_work)
         return
      end if
      write (msg, '(A,": full-column product ",es12.5," vs open depth ",es12.5, &
            &" -- they must DIFFER here or this case cannot see the defect")') &
         label, d_full*metrics%dy_cu_bt(2, 1)/metrics%dy_cu(2, 1), d_open
      call check(error, abs(d_full*metrics%dy_cu_bt(2, 1)/metrics%dy_cu(2, 1) - d_open) &
                 > 0.1_wp*d_open, trim(msg))
      if (allocated(error)) then
         call teardown_step_face(metrics, ms, bt_work)
         return
      end if

      ! The identity.  `d_bt` is `s·(dy/dy_bt)·dy_bt/dy` with `s` the
      ! open sum: four round-offs on an O(d_open) value.
      tol = 8.0_wp*epsilon(1.0_wp)*d_open
      write (msg, '(A,": fast-loop face depth h_face_up*dy_cu_bt/dy_cu = ",es20.12, &
            &" vs the open upstream column sum = ",es20.12,". The barotropic ", &
            &"transport must equal the renormalised open-layer transport.")') &
         label, d_bt, d_open
      call check(error, abs(d_bt - d_open) <= tol, trim(msg))
      call teardown_step_face(metrics, ms, bt_work)
   end subroutine one_direction

   subroutine test_face_depth(error)
      type(error_type), allocatable, intent(out) :: error
      call one_direction(error, 1.0_wp, "flow off the DEEP side")
      if (allocated(error)) return
      call one_direction(error, -1.0_wp, "flow off the SHELF side")
   end subroutine test_face_depth

   ! ---------------------------------------------------------------
   !  2. substep_drag / wave_drag damp on the open-column depth.
   ! ---------------------------------------------------------------

   subroutine test_bt_rem_open(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp), parameter :: R_LIN = 2.5e-4_wp, HBBL = 10.0_wp, DT_IN = 30.0_wp
      real(wp), parameter :: R_H = 1.0e-2_wp
      real(wp) :: h_w(NZ), h_e(NZ), u_k(NZ), h_open, h_full, expect, full_ans, rem_bt
      integer :: k
      character(len=400) :: msg

      h_w = [H_NOM, H_NOM, H_NOM, H_NOM]
      h_e = [H_MIN, H_MIN, H_NOM, H_NOM]
      u_k = 0.1_wp
      call build_step_face(grid, metrics, ms, bt_work, h_w, h_e, u_k)

      h_open = 0.0_wp
      h_full = 0.0_wp
      do k = 1, NZ
         h_open = h_open + 0.5_wp*(h_w(k) + h_e(k))*metrics%open_u(2, 1, k)
         h_full = h_full + 0.5_wp*(h_w(k) + h_e(k))
      end do

      ! ---- substep_drag ----
      call compute_bt_rem(grid, bt_work, ms, metrics, R_LIN, HBBL, DT_IN)
      !$acc update self(bt_work%bt_rem_u)
      rem_bt = bt_work%bt_rem_u(2, 1)
      expect = h_open/(h_open + R_LIN*HBBL*DT_IN)
      full_ans = h_full/(h_full + R_LIN*HBBL*DT_IN)
      write (msg, '("substep_drag: bt_rem_u = ",es20.12," vs open-column ",es20.12, &
            &" (full-column would be ",es20.12,")")') rem_bt, expect, full_ans
      call check(error, abs(rem_bt - expect) <= 4.0_wp*epsilon(1.0_wp) .and. &
                 abs(full_ans - expect) > 1.0e-6_wp, trim(msg))
      if (allocated(error)) then
         call teardown_step_face(metrics, ms, bt_work)
         return
      end if

      ! ---- wave_drag (MULTIPLIES into bt_rem: reset to 1 first) ----
      allocate (bt_work%lwd_drag_u(grid%nx_total + 1, grid%ny_total), source=R_H)
      allocate (bt_work%lwd_drag_v(grid%nx_total, grid%ny_total + 1), source=R_H)
      !$acc enter data copyin(bt_work%lwd_drag_u, bt_work%lwd_drag_v)
      bt_work%bt_rem_u = 1.0_wp
      bt_work%bt_rem_v = 1.0_wp
      !$acc update device(bt_work%bt_rem_u, bt_work%bt_rem_v)
      call compute_bt_rem_wave_drag(grid, bt_work, ms, metrics, DT_IN)
      !$acc update self(bt_work%bt_rem_u)
      !$acc exit data delete(bt_work%lwd_drag_u, bt_work%lwd_drag_v)
      rem_bt = bt_work%bt_rem_u(2, 1)
      expect = h_open/(h_open + R_H*DT_IN)
      full_ans = h_full/(h_full + R_H*DT_IN)
      write (msg, '("wave_drag: bt_rem_u = ",es20.12," vs open-column ",es20.12, &
            &" (full-column would be ",es20.12,")")') rem_bt, expect, full_ans
      call check(error, abs(rem_bt - expect) <= 4.0_wp*epsilon(1.0_wp) .and. &
                 abs(full_ans - expect) > 1.0e-6_wp, trim(msg))
      deallocate (bt_work%lwd_drag_u, bt_work%lwd_drag_v)
      call teardown_step_face(metrics, ms, bt_work)
   end subroutine test_bt_rem_open

   ! ---------------------------------------------------------------
   !  3. Integration: the staircase basin, upstream_h_face on.
   ! ---------------------------------------------------------------

   subroutine run_basin(stratified, seiche, f0, n_periods, finite, energy_ratio, &
                        period_meas, n_closed, closed_u_max, u_max)
      !! The `test_ocean_zfixed_bt_seiche` staircase basin (400 m west,
      !! 200 m shelf east) with `upstream_h_face` on, the production
      !! per-outer-step `ocean_porous_refresh`, `pred_corr`, f = 0, no
      !! wind / drag / viscosity / diffusion / mixing / thermodynamics.
      logical, intent(in) :: stratified
         !! Layer density `ρ0 + DRHO·(NZ−k)/(NZ−1)` (dense at the bed).
      logical, intent(in) :: seiche
         !! Seed the gravest seiche (`η = ETA0·cos(πx/L)`); `.false.` ⇒ rest.
      real(wp), intent(in) :: f0
         !! Coriolis parameter (1/s).
      real(wp), intent(in) :: n_periods
         !! Run length in units of the f = 0 stepped-basin seiche period.
      logical, intent(out) :: finite
      real(wp), intent(out) :: energy_ratio, period_meas, closed_u_max, u_max
      integer, intent(out) :: n_closed

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc

      integer :: i, j, k, ig, i0, i1, j0, j1, step, n_steps, n_cross, nx_t, ny_t
      real(wp) :: energy0, energy1, x_rel, eta_seed, period_analytic
      real(wp) :: t_first_cross, t_last_cross, d_prev, d_now
      real(wp), allocatable :: tgt(:, :, :), tot_h(:, :), eta0f(:, :), z_top(:, :)

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = f0
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_FV_LITE
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_Z_FIXED
      vc%z_fixed_h_ref = H_NOM*real(NZ, wp)
      vc%zstar_h_min = H_MIN
      vc%zfixed_closed_faces = .true.
      dyn%split_scheme = SPLIT_SCHEME_PRED_CORR

      ! THE knob under test.  Allocated exactly as `configure_ocean_bt_split`
      ! does, BEFORE `dyn%enter_data` so the device map carries it.
      dyn%bt_work%use_upstream_h_face = .true.
      allocate (dyn%bt_work%h_face_up_x(nx_t + 1, ny_t), source=0.0_wp)
      allocate (dyn%bt_work%h_face_up_y(nx_t, ny_t + 1), source=0.0_wp)

      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      vd%zlevel_faces = .true.
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP

      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)

      allocate (tot_h(nx_t, ny_t), source=H_DEEP)
      allocate (eta0f(nx_t, ny_t), source=0.0_wp)
      allocate (z_top(nx_t, ny_t), source=0.0_wp)
      allocate (tgt(nx_t, ny_t, NZ), source=0.0_wp)
      do j = 1, ny_t
         do i = 1, nx_t
            if (i > ig + NXP/2) tot_h(i, j) = H_SHELF
         end do
      end do
      call ocean_vcoord_z_fixed_target_uniform(tgt, tot_h, eta0f, z_top, &
                                               nx_t, ny_t, NZ, H_NOM, H_MIN)
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                          tgt, nx_t, ny_t, NZ, H_VANISHED)
      metrics%use_closed_faces = .true.
      ! The mask kernel wrote the DEVICE copy (masks are mapped); pull it
      ! back for the host census and scans.  Inert on host builds.
      !$acc update self(metrics%open_u, metrics%open_v)

      n_closed = 0
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1 + 1
               if (metrics%open_u(i, j, k) == 0.0_wp) n_closed = n_closed + 1
            end do
         end do
      end do

      do j = 1, ny_t
         do i = 1, nx_t
            eta_seed = 0.0_wp
            if (seiche) then
               x_rel = (real(i - ig, wp) - 0.5_wp)/real(NXP, wp)
               eta_seed = ETA0*cos(PI_L*x_rel)
            end if
            do k = 1, NZ
               ms%h_layer(i, j, k) = tgt(i, j, k)
               if (stratified) then
                  ms%rho_layer(i, j, k) = eos%rho0 + DRHO*real(NZ - k, wp)/real(NZ - 1, wp)
               else
                  ms%rho_layer(i, j, k) = eos%rho0
               end if
            end do
            ms%h_layer(i, j, NZ) = ms%h_layer(i, j, NZ) + eta_seed
            do k = 1, NZ
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = 0.0_wp
            dyn%bt_work%bt_H_ref(i, j) = tot_h(i, j)
         end do
      end do
      ms%u_face_x_layer(nx_t + 1, :, :) = 0.0_wp
      ms%v_face_y_layer(:, ny_t + 1, :) = 0.0_wp

      period_analytic = stepped_basin_period()
      n_steps = nint(n_periods*period_analytic/DT)

      energy0 = basin_energy(ms, dyn, i0, i1, j0, j1)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      closed_u_max = 0.0_wp
      u_max = 0.0_wp
      n_cross = 0
      t_first_cross = 0.0_wp
      t_last_cross = 0.0_wp
      d_prev = ssh_tilt(ms, dyn, i0, i1, j0, j1)

      do step = 1, n_steps
         ! Production cadence: the engine refreshes the narrowed BT widths
         ! from the live `h` once per outer step.  Without it `dy_cu_bt`
         ! stays at the un-narrowed seed and the fast loop runs on the
         ! FULL column whatever `h_face_up` says.
         call ocean_porous_refresh(grid, metrics, ms)
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         call scan_faces(ms, metrics, i0, i1, j0, j1, closed_u_max, u_max)
         d_now = ssh_tilt(ms, dyn, i0, i1, j0, j1)
         if (d_prev*d_now < 0.0_wp) then
            n_cross = n_cross + 1
            if (n_cross == 1) t_first_cross = real(step, wp)*DT
            t_last_cross = real(step, wp)*DT
         end if
         d_prev = d_now
      end do

      finite = .true.
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) finite = .false.
            end do
         end do
      end do

      energy1 = basin_energy(ms, dyn, i0, i1, j0, j1)
      if (finite .and. energy0 > 0.0_wp .and. ieee_is_finite(energy1)) then
         energy_ratio = energy1/energy0
      else if (finite .and. energy0 == 0.0_wp) then
         energy_ratio = energy1
      else
         energy_ratio = huge(1.0_wp)
      end if
      if (n_cross >= 3) then
         period_meas = 2.0_wp*(t_last_cross - t_first_cross)/real(n_cross - 1, wp)
      else
         period_meas = 0.0_wp
      end if

      call vc%exit_data()
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data()
      call ss%exit_data(); call bd%exit_data(); call hv%exit_data()
      call pgf%exit_data(); call cor%exit_data(); call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
      deallocate (tgt, tot_h, eta0f, z_top)
   end subroutine run_basin

   pure subroutine scan_faces(ms, metrics, i0, i1, j0, j1, closed_u_max, u_max)
      !! Worst `|u|` on a CLOSED u-face (must be exactly 0) and anywhere.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: i0, i1, j0, j1
      real(wp), intent(inout) :: closed_u_max, u_max
      integer :: i, j, k

      do j = j0, j1
         do i = i0 + 1, i1
            do k = 1, NZ
               u_max = max(u_max, abs(ms%u_face_x_layer(i, j, k)))
               if (metrics%open_u(i, j, k) == 0.0_wp) then
                  closed_u_max = max(closed_u_max, abs(ms%u_face_x_layer(i, j, k)))
               end if
            end do
         end do
      end do
      do j = j0 + 1, j1
         do i = i0, i1
            do k = 1, NZ
               u_max = max(u_max, abs(ms%v_face_y_layer(i, j, k)))
            end do
         end do
      end do
   end subroutine scan_faces

   pure function ssh_tilt(ms, dyn, i0, i1, j0, j1) result(tilt)
      !! `eta(west end) − eta(east end)`, the gravest-mode amplitude.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_dyn_t), intent(in) :: dyn
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: tilt
      integer :: j, k
      real(wp) :: w, e

      w = 0.0_wp
      e = 0.0_wp
      do j = j0, j1
         do k = 1, NZ
            w = w + ms%h_layer(i0, j, k)
            e = e + ms%h_layer(i1, j, k)
         end do
         w = w - dyn%bt_work%bt_H_ref(i0, j)
         e = e - dyn%bt_work%bt_H_ref(i1, j)
      end do
      tilt = w - e
   end function ssh_tilt

   pure function basin_energy(ms, dyn, i0, i1, j0, j1) result(energy)
      !! `Σ_k h·(u² + v²)/2 + g·η²/2` per unit area and density.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_dyn_t), intent(in) :: dyn
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: energy
      integer :: i, j, k
      real(wp) :: eta, h_col

      energy = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            h_col = 0.0_wp
            do k = 1, NZ
               h_col = h_col + ms%h_layer(i, j, k)
               energy = energy + 0.5_wp*ms%h_layer(i, j, k)* &
                        (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
            end do
            eta = h_col - dyn%bt_work%bt_H_ref(i, j)
            energy = energy + 0.5_wp*GRAV*eta*eta
         end do
      end do
   end function basin_energy

   subroutine test_resting(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: fin
      real(wp) :: ratio, period_meas, closed_u_max, u_max
      integer :: n_closed
      character(len=400) :: msg

      call run_basin(.true., .false., 0.0_wp, REST_PERIODS, fin, ratio, period_meas, n_closed, closed_u_max, u_max)
      call check(error, n_closed > 0, "resting staircase: the mask closed nothing -- vacuous")
      if (allocated(error)) return
      call check(error, fin, "resting staircase: went non-finite")
      if (allocated(error)) return
      write (msg, '("resting stratified staircase, upstream_h_face on: max|u| = ",es11.3, &
            &" m/s over the run (bar ",es9.2,"). The open-face PGF is identically ", &
            &"zero here, so any motion is manufactured.")') u_max, REST_U_BAR
      call check(error, u_max <= REST_U_BAR, trim(msg))
   end subroutine test_resting

   subroutine test_seiche(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: fin
      real(wp) :: ratio, period_meas, closed_u_max, u_max, period_analytic, rel
      integer :: n_closed
      character(len=400) :: msg

      call run_basin(.false., .true., 0.0_wp, real(N_PERIODS, wp), fin, ratio, period_meas, n_closed, closed_u_max, u_max)
      write (msg, '("seiche: the staircase closed ",I0," interior u-face entries; zero ", &
            &"means every assertion below is vacuous")') n_closed
      call check(error, n_closed > 0, trim(msg))
      if (allocated(error)) return
      call check(error, fin, "seiche, upstream_h_face on: went non-finite")
      if (allocated(error)) return
      write (msg, '("seiche: max |u| on a CLOSED face = ",es11.3," (must be exactly 0)")') &
         closed_u_max
      call check(error, closed_u_max == 0.0_wp, trim(msg))
      if (allocated(error)) return
      write (msg, '("seiche, upstream_h_face on: KE+PE ratio = ",es11.3," over ",I0, &
            &" periods (bar ",F5.2,"). Unforced and undamped: growth is manufactured ", &
            &"-- the fast loop must transport on the OPEN upstream column.")') &
         ratio, N_PERIODS, ENERGY_GROWTH_BAR
      call check(error, ratio < ENERGY_GROWTH_BAR, trim(msg))
      if (allocated(error)) return
      period_analytic = stepped_basin_period()
      rel = abs(period_meas - period_analytic)/period_analytic
      write (msg, '("seiche: period ",es11.3," s vs the exact stepped-basin ",es11.3, &
            &" s (",f6.2,"% apart)")') period_meas, period_analytic, 100.0_wp*rel
      call check(error, period_meas > 0.0_wp .and. rel < PERIOD_TOL, trim(msg))
   end subroutine test_seiche

   subroutine test_rotating(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: fin
      real(wp) :: ratio, period_meas, closed_u_max, u_max
      integer :: n_closed
      character(len=400) :: msg

      call run_basin(.false., .true., F_ROT, ROT_PERIODS, fin, ratio, period_meas, &
                     n_closed, closed_u_max, u_max)
      call check(error, n_closed > 0, "rotating: the mask closed nothing -- vacuous")
      if (allocated(error)) return
      call check(error, fin, "rotating staircase, upstream_h_face on: went non-finite")
      if (allocated(error)) return
      write (msg, '("rotating: max |u| on a CLOSED face = ",es11.3," (must be exactly 0)")') &
         closed_u_max
      call check(error, closed_u_max == 0.0_wp, trim(msg))
      if (allocated(error)) return
      write (msg, '("rotating staircase, upstream_h_face on: KE+PE ratio = ",es11.3, &
            &" (bar ",F5.2,")")') ratio, ENERGY_GROWTH_BAR
      call check(error, ratio < ENERGY_GROWTH_BAR, trim(msg))
   end subroutine test_rotating

   pure function stepped_basin_period() result(period)
      !! Exact gravest-mode period of a closed 1-D basin of length `L`,
      !! depth `H1 = H_DEEP` on the western half and `H2 = H_SHELF` on the
      !! eastern half.  `η = A·cos(k1·x)` west, `B·cos(k2·(L−x))` east,
      !! `k = ω/√(gH)`; continuity of `η` and of the transport `H·∂η/∂x`
      !! at `x = L/2` gives
      !!
      !! ```
      !! F(ω) = √H1·tan(ω·L/(2c1)) + √H2·tan(ω·L/(2c2)) = 0
      !! ```
      !!
      !! whose smallest positive root lies in `(π·c2/L, π·c1/L)` (the
      !! shelf argument past `π/2`, the deep one below it), where `F` runs
      !! monotonically from `−∞` to `+∞`.  Bisection to round-off.
      real(wp) :: period
      real(wp) :: c1, c2, el, lo, hi, mid, f_mid
      integer :: it

      c1 = sqrt(GRAV*H_DEEP)
      c2 = sqrt(GRAV*H_SHELF)
      el = real(NXP, wp)*DX
      lo = PI_L*c2/el*(1.0_wp + 1.0e-12_wp)
      hi = PI_L*c1/el*(1.0_wp - 1.0e-12_wp)
      do it = 1, 200
         mid = 0.5_wp*(lo + hi)
         f_mid = sqrt(H_DEEP)*tan(mid*el/(2.0_wp*c1)) + sqrt(H_SHELF)*tan(mid*el/(2.0_wp*c2))
         if (f_mid < 0.0_wp) then
            lo = mid
         else
            hi = mid
         end if
      end do
      period = 2.0_wp*PI_L/(0.5_wp*(lo + hi))
   end function stepped_basin_period

end module test_ocean_bt_upstream_zfixed
