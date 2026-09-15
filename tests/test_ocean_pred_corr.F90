!! Dedicated coverage for the `split_scheme = "pred_corr"`
!! (`SPLIT_SCHEME_PRED_CORR`) predictor-corrector split integrator in
!! `ocean_dyn_step_split`.  The scheme is otherwise only exercised
!! incidentally (`test_ocean_bt_substep_zeta_ke` drives it to probe the
!! fast-loop knob), leaving the predictor/corrector trajectory itself
!! without a direct guard.
!!
!! Three properties, in ascending sensitivity:
!!
!!   1. REST INVARIANT (the high-leverage gate) — a quiescent basin
!!      (zero velocity, flat free surface, uniform density, no forcing)
!!      integrated a dozen pred_corr steps must stay at rest to round-off.
!!      The predictor's discarded-h advance, the between-stage
!!      barotropic reset (SPEC §2), and the S1 time-mean seed each get a
!!      chance to inject a spurious acceleration into an ocean that
!!      should not move; any of those bugs shows up here as `|u|` growing
!!      far above round-off.  This is the pattern that has repeatedly
!!      caught real bugs in this core.
!!
!!   2. MASS CONSERVATION under a NON-trivial trajectory — a depth-
!!      uniform barotropic shear jet (so the corrector integrates a
!!      genuine `∇·(hu)`) must conserve total column mass `Σh` to
!!      round-off over the same run.  A predictor whose h advance leaks
!!      into the corrector (instead of being weighted out) breaks this.
!!
!!   3. CROSS-SCHEME STABILITY — the same benign jet run under BOTH
!!      `pred_corr` and `ssp_rk2` stays finite and conserves mass.  The two
!!      integrators are NOT expected to agree bit-for-bit (different
!!      schemes); the guard is only that pred_corr is not categorically
!!      less stable than the reference SSP-RK2 path on a case the latter
!!      handles cleanly.
module test_ocean_pred_corr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, &
                            SPLIT_SCHEME_PRED_CORR, SPLIT_SCHEME_SSP_RK2
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   use rdb_config, only: config_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_pred_corr_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 32
   integer, parameter :: NYP = 8
   integer, parameter :: NZ = 2
   integer, parameter :: N_STEPS = 12
   integer, parameter :: N_INNER = 24
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: H0 = 200.0_wp        !! flat total depth (m)
   real(wp), parameter :: V0 = 0.3_wp          !! barotropic jet amplitude (m/s)
   real(wp), parameter :: F0 = 1.0e-4_wp       !! f-plane Coriolis (1/s)
   real(wp), parameter :: DT = 60.0_wp
   real(wp), parameter :: JET_CELLS = 8.0_wp   !! jet wavelength (grid cells)
   real(wp), parameter :: PI_L = 3.14159265358979324_wp

contains

   subroutine collect_ocean_pred_corr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("pred_corr_rest_invariant", test_rest_invariant), &
                  new_unittest("pred_corr_mass_conserved_active", test_mass_conserved), &
                  new_unittest("pred_corr_vs_ssp_rk2_both_stable", test_cross_scheme_stable), &
                  new_unittest("pred_corr_is_the_default", test_default_in_step) &
                  ]
   end subroutine collect_ocean_pred_corr_tests

   subroutine test_default_in_step(error)
      !! The two halves of the default live in different modules -- the
      !! `&ocean_bt_nml split_scheme` string in `rdb_config.F90` and the
      !! `ocean_dyn_t%split_scheme` enum here -- and `ocean_setup` only
      !! ever writes the second from the first, so a divergence is
      !! invisible to every production run and silently hides a scheme
      !! from the whole Fortran suite (which builds bare `ocean_dyn_t`
      !! objects).  They have diverged before.  Assert them equal.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      type(config_t) :: cfg

      call check(error, dyn%split_scheme == SPLIT_SCHEME_PRED_CORR, &
                 "ocean_dyn_t%split_scheme default must be SPLIT_SCHEME_PRED_CORR")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%bt%split_scheme) == "pred_corr", &
                 "&ocean_bt_nml split_scheme default must be 'pred_corr'")
   end subroutine test_default_in_step

   subroutine run_case(split_scheme, at_rest, finite, mass_rel_err, max_speed)
      !! Drive `N_STEPS` of `ocean_dyn_step_split` on a flat two-layer
      !! f-plane basin under `split_scheme`.  `at_rest = .true.` seeds a
      !! quiescent state (all velocities zero); `.false.` seeds a
      !! depth-uniform (barotropic) sinusoidal v-jet in x.  Undamped and
      !! unforced either way — Lagrangian vcoord (no remap) so the layers
      !! are material and a rest state is a fixed point by construction.
      integer, intent(in) :: split_scheme
      logical, intent(in) :: at_rest
      logical, intent(out) :: finite
      real(wp), intent(out) :: mass_rel_err
      real(wp), intent(out) :: max_speed
         !! Max |velocity| over the interior at the end of the run.

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

      integer :: i, j, k, ig, i0, i1, j0, j1, step
      real(wp) :: x_f, mass0, mass1, vjet

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = F0
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
      vc%coord_type = VCOORD_LAGRANGIAN
      dyn%split_scheme = split_scheme

      ! Fully undamped + unforced: the only trajectory driver is the
      ! split integrator itself, so a rest state stays exactly at rest.
      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP

      ! Flat two-layer column; uniform density (=> zero PGF); optional
      ! depth-uniform (barotropic) v-jet in x.
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            ms%h_layer(i, j, 1) = 0.5_wp*H0
            ms%h_layer(i, j, 2) = 0.5_wp*H0
            ms%rho_layer(i, j, :) = eos%rho0
            do k = 1, NZ
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            if (at_rest) then
               vjet = 0.0_wp
            else
               x_f = real(i - ig, wp)
               vjet = V0*sin(2.0_wp*PI_L*x_f/JET_CELLS)
            end if
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = vjet
            dyn%bt_work%bt_H_ref(i, j) = H0
         end do
      end do
      ms%u_face_x_layer(grid%nx_total + 1, :, :) = 0.0_wp
      ms%v_face_y_layer(:, grid%ny_total + 1, :) = 0.0_wp

      mass0 = 0.0_wp
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1
               mass0 = mass0 + ms%h_layer(i, j, k)
            end do
         end do
      end do

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
      end do

      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      finite = .true.
      mass1 = 0.0_wp
      max_speed = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               mass1 = mass1 + ms%h_layer(i, j, k)
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) finite = .false.
               max_speed = max(max_speed, abs(ms%u_face_x_layer(i, j, k)))
               max_speed = max(max_speed, abs(ms%v_face_y_layer(i, j, k)))
            end do
         end do
      end do
      mass_rel_err = abs(mass1 - mass0)/mass0

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
   end subroutine run_case

   subroutine test_rest_invariant(error)
      !! Quiescent rest under pred_corr must stay at rest to round-off:
      !! no spurious acceleration, mass exactly conserved, all finite.
      type(error_type), allocatable, intent(out) :: error
      logical :: finite
      real(wp) :: mass_rel_err, max_speed
      character(len=256) :: msg

      call run_case(SPLIT_SCHEME_PRED_CORR, at_rest=.true., finite=finite, &
                    mass_rel_err=mass_rel_err, max_speed=max_speed)

      call check(error, finite, "pred_corr rest: state went non-finite")
      if (allocated(error)) return

      write (msg, '("pred_corr rest: max |vel| = ", es11.3, &
            &" m/s (must be round-off; nonzero => spurious acceleration)")') max_speed
      call check(error, max_speed < 1.0e-10_wp, trim(msg))
      if (allocated(error)) return

      write (msg, '("pred_corr rest: mass rel err = ", es11.3, " (must be round-off)")') &
         mass_rel_err
      call check(error, mass_rel_err < 1.0e-12_wp, trim(msg))
   end subroutine test_rest_invariant

   subroutine test_mass_conserved(error)
      !! A live barotropic shear jet under pred_corr: the corrector
      !! integrates a genuine divergence, yet total column mass must be
      !! conserved to round-off (predictor h advance weighted out).
      type(error_type), allocatable, intent(out) :: error
      logical :: finite
      real(wp) :: mass_rel_err, max_speed
      character(len=256) :: msg

      call run_case(SPLIT_SCHEME_PRED_CORR, at_rest=.false., finite=finite, &
                    mass_rel_err=mass_rel_err, max_speed=max_speed)

      call check(error, finite, "pred_corr jet: state went non-finite")
      if (allocated(error)) return

      write (msg, '("pred_corr jet: mass rel err = ", es11.3, " (must be round-off)")') &
         mass_rel_err
      call check(error, mass_rel_err < 1.0e-12_wp, trim(msg))
   end subroutine test_mass_conserved

   subroutine test_cross_scheme_stable(error)
      !! The SAME benign jet under pred_corr and ssp_rk2 both stay finite
      !! and conserve mass — not a bit-equality check (different
      !! schemes), just a guard that pred_corr is no less stable than the
      !! reference on a case ssp_rk2 handles cleanly.
      type(error_type), allocatable, intent(out) :: error
      logical :: fin_pc, fin_rk2
      real(wp) :: merr_pc, merr_rk2, spd_pc, spd_rk2

      call run_case(SPLIT_SCHEME_PRED_CORR, at_rest=.false., finite=fin_pc, &
                    mass_rel_err=merr_pc, max_speed=spd_pc)
      call run_case(SPLIT_SCHEME_SSP_RK2, at_rest=.false., finite=fin_rk2, &
                    mass_rel_err=merr_rk2, max_speed=spd_rk2)

      call check(error, fin_rk2, "ssp_rk2 jet: reference went non-finite")
      if (allocated(error)) return
      call check(error, fin_pc, "pred_corr jet: went non-finite where ssp_rk2 stayed finite")
      if (allocated(error)) return
      call check(error, merr_pc < 1.0e-12_wp, "pred_corr jet: mass not conserved")
      if (allocated(error)) return
      call check(error, merr_rk2 < 1.0e-12_wp, "ssp_rk2 jet: mass not conserved")
   end subroutine test_cross_scheme_stable

end module test_ocean_pred_corr
