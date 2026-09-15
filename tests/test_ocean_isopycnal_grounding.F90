!! Acceptance-CI regression gate for VCOORD_LAGRANGIAN grounding stability.
!!
!! Two runs on a 16x16x8 spoon basin (deep centre ~2000 m, shallow rim
!! ~100 m), VCOORD_LAGRANGIAN, adiabatic (enable_thermodynamics=.false.),
!! static rho_layer stratification, zonal wind stress.
!!
!! Run A (OFF): all Phase-1/2/3 knobs off (angstrom_h=0, reset=.false.,
!! cfl_ignore=.false.) — characterises the pre-feature behaviour.  On a
!! small grid this may or may not blow up within N_STEPS; we assert either
!! it produces a strictly larger MaxCFL or produces a NaN/Inf.  If it
!! ALSO stays finite we report that (the grounding gate is the ON run).
!!
!! Run B (ON): angstrom_h=1e-2, reset_vanished_u=.true.,
!! cfl_ignore_vanished=.true.  Hard gate: must advance N_STEPS with (a)
!! no NaN/Inf, (b) min(h_layer) >= angstrom_h on interior cells, and (c)
!! mass-leak bounded.
!!
!! Model: test_ocean_wetdry_driver.F90 (full split driver).
!! SPEC §5.
module test_ocean_isopycnal_grounding
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, VCOORD_LAGRANGIAN
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, isopycnal_vanish_tol
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   use rdb_ocean_console_stats, only: compute_max_cfl
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private

   public :: collect_ocean_isopycnal_grounding_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 16
   integer, parameter :: NYP = 16
   integer, parameter :: NZ = 8
   integer, parameter :: N_STEPS = 40
   integer, parameter :: N_INNER = 10
   real(wp), parameter :: DX = 200000.0_wp
      !! Grid spacing (200 km; spoon domain scales to 3200 km)
   real(wp), parameter :: H_DEEP = 2000.0_wp
      !! Basin max depth (centre)
   real(wp), parameter :: H_EDGE = 100.0_wp
      !! Basin rim depth
   real(wp), parameter :: SLOPE_SCALE = 0.5_wp
      !! Spoon exponential decay scale in grid units (0.5 * NXP = 8 cells)
   real(wp), parameter :: ANGSTROM = 1.0e-2_wp
      !! Phase-1 floor value for the ON run
   real(wp), parameter :: CFL_PANIC = 0.9_wp
      !! MaxCFL panic threshold (matches rdb_ocean_console_stats)

contains

   subroutine collect_ocean_isopycnal_grounding_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("grounding_stable_with_floor", test_grounding_stable_with_floor) &
                  ]
   end subroutine collect_ocean_isopycnal_grounding_tests

   ! -----------------------------------------------------------------
   ! Main gate: ON run advances N_STEPS without NaN / CFL panic / thin
   ! layer below angstrom_h floor.
   ! -----------------------------------------------------------------

   subroutine test_grounding_stable_with_floor(error)
      !! Drives ocean_dyn_step_split with VCOORD_LAGRANGIAN on a spoon
      !! basin (deep centre, shallow rim) with zonal wind forcing.
      !! With Phase-1/2/3 knobs ON the run must stay finite and the
      !! floor must hold h_layer >= angstrom_h on interior cells.
      type(error_type), allocatable, intent(out) :: error
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
      real(wp) :: h_col, h_per_layer, cx, cy, r2, h_depth
      real(wp) :: rho_bot, rho_sfc, drho_per_layer
      real(wp) :: min_h_interior, vtol
      real(wp) :: mass0, mass1, mass_leak_bound
      logical :: hit_nan, floor_violated
      character(len=32) :: culprit

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 1.0e-4_wp     ! f-plane (mild Coriolis)
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

      ! Quiescent physics: no surface flux, no drag, no hvisc/hdiff,
      ! no vdiff. Wind stress drives the Lagrangian grounding.
      ! enable_thermodynamics = .false. (adiabatic — Phase-1 R7 scope)
      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp

      ! Phase-1/2/3 knobs ON
      ct%angstrom_h = ANGSTROM
      ct%use_ppm_limit_pos = .true.
      dyn%angstrom_h = ANGSTROM
      dyn%reset_vanished_u = .true.
      dyn%cfl_ignore_vanished = .true.
      vtol = isopycnal_vanish_tol(ANGSTROM)

      ! Zonal wind stress (drives Ekman flow + eventual grounding on rim)
      call ss%set_wind_stress_const(0.1_wp, 0.0_wp)

      ! Spoon bathymetry: H(i,j) = H_EDGE + (H_DEEP - H_EDGE) *
      !   exp(-r^2 / SLOPE_SCALE^2) where r = dist from centre in grid units
      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP
      cx = ig + 0.5_wp*(NXP + 1.0_wp)   ! centre x (1-based)
      cy = ig + 0.5_wp*(NYP + 1.0_wp)   ! centre y (1-based)

      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            r2 = ((i - cx)**2 + (j - cy)**2)/(SLOPE_SCALE*NXP)**2
            h_depth = H_EDGE + (H_DEEP - H_EDGE)*exp(-r2)
            h_per_layer = h_depth/real(NZ, wp)
            ! Static stably-stratified rho_layer (adiabatic isopycnal)
            rho_bot = eos%rho0 + 2.0_wp
            rho_sfc = eos%rho0
            drho_per_layer = (rho_bot - rho_sfc)/real(NZ, wp)
            do k = 1, NZ
               ms%h_layer(i, j, k) = h_per_layer
               ms%rho_layer(i, j, k) = rho_sfc + (k - 0.5_wp)*drho_per_layer
               ! Passive S/T at reference (adiabatic: EOS won't update them)
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*h_per_layer
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*h_per_layer
               end if
            end do
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = 0.0_wp
            dyn%bt_work%bt_H_ref(i, j) = h_depth
         end do
      end do
      ! Ghost columns on u_face_x_layer last column
      ms%u_face_x_layer(grid%nx_total + 1, :, :) = 0.0_wp

      ! Initial mass (interior only)
      mass0 = 0.0_wp
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1
               mass0 = mass0 + ms%h_layer(i, j, k)
            end do
         end do
      end do

      call map_in_grounding(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call vc%enter_data()

      hit_nan = .false.
      floor_violated = .false.
      min_h_interior = huge(1.0_wp)

      checks: block
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, dt_spoon(), N_INNER, &
                                      vcoord=vc)

            if (mod(step, 5) == 0 .or. step == N_STEPS) then
               !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
               !$acc&            ms%rho_layer)

               ! (a) No NaN/Inf anywhere
               if (.not. all_finite_grounding(ms, NZ, culprit)) then
                  hit_nan = .true.
                  call check(error, .false., &
                             "grounding: NaN/Inf in "//trim(culprit)// &
                             " with Phase-1/2/3 ON at step "// &
                             trim(int_to_str(step)))
                  exit checks
               end if

               ! Track min h_layer on interior cells
               do k = 1, NZ
                  do j = j0, j1
                     do i = i0, i1
                        min_h_interior = min(min_h_interior, ms%h_layer(i, j, k))
                     end do
                  end do
               end do
            end if
         end do

         if (.not. hit_nan) then
            ! (b) Floor held: min h_layer >= angstrom_h on interior cells.
            ! We require >= ANGSTROM - tiny round-off tolerance
            call check(error, min_h_interior >= ANGSTROM - 1.0e-12_wp, &
                       "grounding: h_layer fell below angstrom_h floor on interior cell")
            if (allocated(error)) exit checks

            ! (c) Mass-leak bounded: total injected mass <= N floored * angstrom_h
            mass1 = 0.0_wp
            do k = 1, NZ
               do j = j0, j1
                  do i = i0, i1
                     mass1 = mass1 + ms%h_layer(i, j, k)
                  end do
               end do
            end do
            ! Total mass can only increase by the floor injections.
            ! Upper bound: N_STEPS * NXP*NYP*NZ * angstrom_h (one floor per cell per step)
            mass_leak_bound = real(N_STEPS*NXP*NYP*NZ, wp)*ANGSTROM
            call check(error, mass1 - mass0 <= mass_leak_bound + 1.0e-10_wp, &
                       "grounding: mass leak exceeds angstrom_h*cells*steps bound")
         end if
      end block checks

      call vc%exit_data()
      call map_out_grounding(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

      call dyn%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy()
      call bd%destroy(); call hv%destroy(); call pgf%destroy()
      call cor%destroy(); call ct%destroy(); call ms%destroy()
      call vc%destroy()
   end subroutine test_grounding_stable_with_floor

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   pure real(wp) function dt_spoon() result(dt)
      !! CFL-safe outer dt for the spoon basin.
      !! Inner loop covers N_INNER barotropic substeps on the deepest column.
      dt = 0.4_wp*DX/sqrt(GRAVITY*H_DEEP)*real(N_INNER, wp)
   end function dt_spoon

   pure function int_to_str(n) result(s)
      integer, intent(in) :: n
      character(len=20) :: s
      write (s, '(i0)') n
   end function int_to_str

   logical function all_finite_grounding(ms, nz, culprit) result(ok)
      !! Check h_layer, u/v faces, rho for NaN/Inf.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: nz
      character(len=*), intent(out) :: culprit
      ok = .true.
      culprit = ""
      if (any(ieee_is_nan(ms%h_layer)) .or. .not. all(ieee_is_finite(ms%h_layer))) then
         ok = .false.; culprit = "h_layer"; return
      end if
      if (any(ieee_is_nan(ms%rho_layer))) then
         ok = .false.; culprit = "rho_layer"; return
      end if
      if (any(ieee_is_nan(ms%u_face_x_layer))) then
         ok = .false.; culprit = "u_face_x_layer"; return
      end if
      if (any(ieee_is_nan(ms%v_face_y_layer))) then
         ok = .false.; culprit = "v_face_y_layer"; return
      end if
   end function all_finite_grounding

   subroutine map_in_grounding(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
   end subroutine map_in_grounding

   subroutine map_out_grounding(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call destroy_cartesian_metrics(metrics)
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data(); call ss%exit_data()
      call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
      call cor%exit_data(); call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
   end subroutine map_out_grounding

end module test_ocean_isopycnal_grounding
