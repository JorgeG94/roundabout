!! `&ocean_bt_nml substep_zeta_ke` — MOM6-parity planetary-only BT fast
!! loop (live ζ_bt/∇KE dropped from the substeps, frozen copies retained
!! in `F_bt_*_fast`, Cor_ref reduced to `f·v̄`).
!!
!! Two guards:
!!   1. WIRING — the knob must change the trajectory of a state with
!!      nonzero barotropic ζ and KE gradient (a sinusoidal barotropic
!!      shear jet on an f-plane).  This is the parsed-but-unwired
!!      regression class (`cont_corr_bounds` shipped consuming nothing).
!!   2. SANITY — the planetary-only mode must stay finite and conserve
!!      column mass to round-off over the same run.
!!
!! Default `.true.` bit-identity needs no test here: the gate multiplies
!! by exactly 1.0 (IEEE identity), and the full suite runs the default.
module test_ocean_bt_substep_zeta_ke
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_bt_substep_zeta_ke_tests

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

   subroutine collect_bt_substep_zeta_ke_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bt_planetary_substep_wired_and_sane", &
                               test_planetary_wired_and_sane) &
                  ]
   end subroutine collect_bt_substep_zeta_ke_tests

   subroutine run_jet(live_zeta_ke, u_out, v_out, finite, mass_rel_err)
      !! Barotropic sinusoidal jet v(x) on an f-plane, closed flat basin:
      !! ζ_bt = ∂v/∂x ≠ 0 and KE varies with x, so the live and
      !! planetary-only fast loops integrate different equations.
      logical, intent(in) :: live_zeta_ke
      real(wp), intent(out) :: u_out(NXP, NYP), v_out(NXP, NYP)
      logical, intent(out) :: finite
      real(wp), intent(out) :: mass_rel_err

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
      real(wp) :: x_f, mass0, mass1

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
      dyn%split_scheme = SPLIT_SCHEME_PRED_CORR
      dyn%bt_work%substep_zeta_ke = live_zeta_ke

      ! Undamped: only the fast-loop structure may differ between runs.
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

      ! Flat two-layer column; depth-uniform (barotropic) v jet in x.
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
            x_f = real(i - ig, wp)
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = V0*sin(2.0_wp*PI_L*x_f/JET_CELLS)
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
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               mass1 = mass1 + ms%h_layer(i, j, k)
            end do
            if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, 1))) finite = .false.
            if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, 1))) finite = .false.
            u_out(i - ig, j - ig) = ms%u_face_x_layer(i, j, 1)
            v_out(i - ig, j - ig) = ms%v_face_y_layer(i, j, 1)
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
   end subroutine run_jet

   subroutine test_planetary_wired_and_sane(error)
      !! Live vs planetary-only runs of the same barotropic shear jet
      !! must diverge (knob wired), and the planetary-only run must stay
      !! finite with round-off mass conservation.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_live(NXP, NYP), v_live(NXP, NYP)
      real(wp) :: u_plan(NXP, NYP), v_plan(NXP, NYP)
      real(wp) :: merr_live, merr_plan, dmax
      logical :: fin_live, fin_plan
      character(len=256) :: msg

      call run_jet(.true., u_live, v_live, fin_live, merr_live)
      call run_jet(.false., u_plan, v_plan, fin_plan, merr_plan)

      dmax = max(maxval(abs(u_live - u_plan)), maxval(abs(v_live - v_plan)))
      write (msg, '("live vs planetary max |du| = ", es11.3, &
            &" (zero => knob parsed but unwired)")') dmax
      call check(error, dmax > 1.0e-12_wp, trim(msg))
      if (allocated(error)) return

      write (msg, '("planetary run finite = ", l1)') fin_plan
      call check(error, fin_plan, trim(msg))
      if (allocated(error)) return

      write (msg, '("planetary mass rel err = ", es11.3, &
            &" (must be round-off)")') merr_plan
      call check(error, merr_plan < 1.0e-12_wp, trim(msg))
   end subroutine test_planetary_wired_and_sane

end module test_ocean_bt_substep_zeta_ke
