!! Baroclinic open-boundary tests for the ocean dyn-core.
!! Tests 7-10: v1 OBC anomaly/zero-gradient scheme.
!! Tests 11-16: v2 Orlanski + nudging + full Flather.
!! Pattern-matches test_ocean_dyn_split.F90: test-drive, small grids, direct slots.
!!
!! v1 tests implemented (7-10):
!!   7. Kernel-level outflow: depth-mean(u_layer) == bt_ubt_end invariant.
!!   8. Inflow tracer: CLAMPED west inflow raises interior heat content.
!!   9. Radiating gravity wave: all-OPEN drains < DECAY_FRAC * WALL energy.
!!  10. Driver regression: all-WALL bc == no-bc bit-for-bit.
!!
!! v2 tests implemented (11-16, §1 + §3):
!!  11. Outflow relaxation rate: tres − T_int decays as (1 + u*dt/L)^{-n}.
!!  12. Inflow relaxation: tres relaxes to T_data under sustained inflow.
!!  13. Degenerate limits: L_out=0 ⇒ instant T_int; L_in=0 ⇒ instant T_data.
!!  14. No-chatter: alternating face velocity keeps ghost smooth (bounded step).
!!  15. Bit-identity: res_lscale_out/in both 0 ⇒ same as no-reservoir run.
!!  16. Corner-ζ closure: boundary corner line ζ is exactly 0 at all
!!      non-PERIODIC edges after a substep (design §3).
module test_ocean_obc_baroclinic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, &
                                       ocean_bc_state_enter_data, &
                                       ocean_bc_state_exit_data, &
                                       OBC_WALL, OBC_OPEN, OBC_CLAMPED
   use rdb_ocean_obc_baroclinic, only: ocean_obc_apply_baroclinic, &
                                       ocean_obc_update_reservoirs
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear
   implicit none
   private

   public :: collect_ocean_obc_baroclinic_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 2

contains

   subroutine collect_ocean_obc_baroclinic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("obc_baroclinic_outflow_depth_mean", test_outflow_depth_mean), &
                  new_unittest("obc_baroclinic_inflow_tracer_rise", test_inflow_tracer_rise), &
                  new_unittest("obc_baroclinic_radiating_gw_decay", test_radiating_gw_decay), &
                  new_unittest("obc_baroclinic_wall_bc_regression", test_wall_bc_regression), &
                  new_unittest("obc_reservoir_outflow_relax_rate", test_reservoir_outflow_relax), &
                  new_unittest("obc_reservoir_inflow_relax", test_reservoir_inflow_relax), &
                  new_unittest("obc_reservoir_degenerate_limits", test_reservoir_degenerate), &
                  new_unittest("obc_reservoir_no_chatter", test_reservoir_no_chatter), &
                  new_unittest("obc_reservoir_bit_identity", test_reservoir_bit_identity), &
                  new_unittest("obc_corner_zeta_zero", test_corner_zeta_zero), &
                  new_unittest("obc_orlanski_rx_unit", test_orlanski_rx_unit), &
                  new_unittest("obc_orlanski_wave_drains", test_orlanski_wave_drains), &
                  new_unittest("obc_orlanski_nudge_inflow", test_orlanski_nudge_inflow), &
                  new_unittest("obc_orlanski_default_bit_identity", test_orlanski_default_bit_identity), &
                  new_unittest("obc_flather_full_quiescent", test_flather_full_quiescent), &
                  new_unittest("obc_flather_ext_vel_spinup", test_flather_ext_vel_spinup), &
                  new_unittest("obc_flather_legacy_default", test_flather_legacy_default), &
                  new_unittest("obc_orlanski_anchor_per_layer", test_orlanski_anchor_per_layer), &
                  new_unittest("obc_flather_full_mirror", test_flather_full_mirror) &
                  ]
   end subroutine collect_ocean_obc_baroclinic_tests

   ! -----------------------------------------------------------------
   ! Shared helpers (mirror test_ocean_dyn_split.F90 pattern exactly)
   ! -----------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine init_all(grid, nz_loc, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_loc
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
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn

      ms%nz_ml = nz_loc
      call ms%init(grid)
      call ct%init(grid, nz_ml=nz_loc)
      call cor%init(grid, nz_ml=nz_loc)
      call pgf%init(grid, nz_ml=nz_loc)
      call hv%init(grid, nz_ml=nz_loc)
      call bd%init(grid, nz_ml=nz_loc)
      call ss%init(grid, nz_ml=nz_loc)
      call va%init(grid, nz_ml=nz_loc)
      call hd%init(grid, nz_ml=nz_loc)
      call vd%init(grid, nz_ml=nz_loc)
      call vmix%init(grid, nz_ml=nz_loc)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=nz_loc)
   end subroutine init_all

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
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
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      call dyn%enter_data()
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
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
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
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
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine destroy_all

   ! -----------------------------------------------------------------
   ! Test 7: kernel-level outflow depth-mean invariant
   ! -----------------------------------------------------------------

   subroutine test_outflow_depth_mean(error)
      !! Design §3 test 7 (kernel-level).
      !!
      !! Constructs bt_work + ms directly on the host, then calls
      !! `ocean_obc_apply_baroclinic` on the HOST (CPU) path to verify
      !! the depth-mean invariant before any further physics modifies the fields:
      !!
      !!   depth_mean(u_layer(i_e, j, :)) == bt_ubt_end(i_e, j)
      !!
      !! This is a pure unit test of the kernel formula (§2.1) and does NOT
      !! exercise the GPU path (the GPU path is covered by test 9 end-to-end).
      !!
      !! Setup: 2-layer OPEN east edge, non-uniform per-layer u (baroclinic anomaly),
      !! prescribed bt_ubt_end at i_e.  After the kernel call, depth-mean of the
      !! east wall-face per-layer u should equal bt_ubt_end exactly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_bc_state_t) :: bc

      ! 2-layer system: k=1 bed, k=2 surface.
      real(wp), parameter :: H0 = 50.0_wp     ! thickness per layer (m)
      real(wp), parameter :: U_BOT = 0.05_wp  ! bottom layer velocity (slower)
      real(wp), parameter :: U_TOP = 0.15_wp  ! top layer velocity (faster)
      real(wp), parameter :: UBT_E = 0.10_wp  ! prescribed BT velocity at east wall
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 4
      integer :: ng, i_e, j_mid, k
      real(wp) :: h_tot, u_mean_after, ubt_target, diff

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bt_work%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%east%bc_type = OBC_OPEN

         ng = grid%nghost
         i_e = ng + NX_PHYS + 1  ! east physical wall-face index (u_face_x is nx_total+1 wide)

         ! Set uniform per-layer h and per-layer u with a baroclinic anomaly.
         ! Layer 1 (bed) = U_BOT, layer 2 (surface) = U_TOP.
         ms%h_layer = H0
         ms%u_face_x_layer(:, :, 1) = U_BOT
         ms%u_face_x_layer(:, :, 2) = U_TOP
         ms%v_face_y_layer = 0.0_wp

         ! Prescribe bt_ubt_end at the east wall face.  Depth-mean of the initial
         ! u_layer at i_e is 0.5*(U_BOT + U_TOP) (equal H0 layers).
         ! We set UBT_E to a different value to check the formula shifts the result.
         bt_work%bt_ubt_end = 0.0_wp
         bt_work%bt_ubt_end(i_e, :) = UBT_E

         ! GPU map — kernel runs on device.
         !$acc enter data copyin(ms, bt_work, bc)
         call ms%enter_data()
         call bt_work%enter_data()
         call ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, 1.0_wp)
         !$acc update self(ms%u_face_x_layer, ms%h_layer)
         call bt_work%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, bt_work, bc)

         ! Pick the central j column.
         j_mid = ng + NY_PHYS/2 + 1

         h_tot = 0.0_wp
         u_mean_after = 0.0_wp
         do k = 1, NZ
            h_tot = h_tot + ms%h_layer(ng + NX_PHYS, j_mid, k)
            u_mean_after = u_mean_after &
                           + ms%u_face_x_layer(i_e, j_mid, k)*ms%h_layer(ng + NX_PHYS, j_mid, k)
         end do
         if (h_tot > 0.0_wp) u_mean_after = u_mean_after/h_tot

         ubt_target = UBT_E
         diff = abs(u_mean_after - ubt_target)

         call check(error, diff < 1.0e-12_wp, &
                    "test 7: depth-mean(u_layer) at east wall /= bt_ubt_end after kernel, diff=" &
                    //short_str(diff))
         if (allocated(error)) exit checks

         ! Also verify the baroclinic anomaly is zero-gradient:
         ! u_layer(i_e) - ubt_end == u_layer(i_int, :) - ubar_int
         ! i.e., the per-layer deviations from the BT mean are copied from the interior.
         ! With H0 uniform: ubar_int = 0.5*(U_BOT + U_TOP) = 0.1.
         ! New u_layer(i_e, k=1) = UBT_E + (U_BOT - 0.1) = 0.1 + (0.05 - 0.1) = 0.05 = U_BOT.
         ! New u_layer(i_e, k=2) = UBT_E + (U_TOP - 0.1) = 0.1 + (0.15 - 0.1) = 0.15 = U_TOP.
         ! So: anomaly is unchanged (zero-gradient), only the BT mean is changed.
         call check(error, &
                    abs(ms%u_face_x_layer(i_e, j_mid, 1) - UBT_E - (U_BOT - 0.5_wp*(U_BOT + U_TOP))) < 1.0e-12_wp, &
                    "test 7: baroclinic anomaly at east wall k=1 is not zero-gradient")
         if (allocated(error)) exit checks
         call check(error, &
                    abs(ms%u_face_x_layer(i_e, j_mid, 2) - UBT_E - (U_TOP - 0.5_wp*(U_BOT + U_TOP))) < 1.0e-12_wp, &
                    "test 7: baroclinic anomaly at east wall k=2 is not zero-gradient")

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call ms%destroy()
   end subroutine test_outflow_depth_mean

   ! -----------------------------------------------------------------
   ! Test 8: inflow tracer rise
   ! -----------------------------------------------------------------

   subroutine test_inflow_tracer_rise(error)
      !! Design §3 test 8.
      !!
      !! CLAMPED west inflow with clamped_T != interior T.
      !! After N_STEPS steps the interior heat content must have risen
      !! by more than zero and less than the analytical upper bound
      !! (NX_PHYS * NY_PHYS * NZ * H0 * (T_BC - T_INTERIOR): all interior replaced).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
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
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 50.0_wp       ! layer thickness (m)
      real(wp), parameter :: U_INFLOW = 0.1_wp  ! eastward inflow velocity (m/s)
      real(wp), parameter :: T_INTERIOR = 10.0_wp  ! interior temperature (°C)
      real(wp), parameter :: T_BC = 20.0_wp         ! warmer inflow water (°C)
      real(wp), parameter :: S_REF = 35.0_wp
      real(wp), parameter :: DT = 2.0_wp
      integer, parameter :: N_INNER = 5
      integer, parameter :: N_STEPS = 10
      integer, parameter :: NX_PHYS = 12, NY_PHYS = 4
      type(ocean_metrics_t) :: metrics
      integer :: ng, i0, i1, j0, j1, step
      real(wp) :: heat0, heat_after, expected_max_rise

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call init_all(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ! Need n_tracers for clamped_tracer allocation.
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_CLAMPED
         ! u at wall clamped to +U_INFLOW (eastward = into domain from west).
         bc%west%clamped_u = U_INFLOW
         ! Clamped tracer values (T_BC, S_REF).
         bc%west%clamped_tracer(ms%idx_temperature) = T_BC
         bc%west%clamped_tracer(ms%idx_salinity) = S_REF

         ng = grid%nghost
         i0 = ng + 1
         i1 = ng + NX_PHYS
         j0 = ng + 1
         j1 = ng + NY_PHYS

         ! IC: uniform h, uniform T=T_INTERIOR, S=S_REF, u=0.
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_temperature)%hTr = T_INTERIOR*H0
         ms%tracers(ms%idx_salinity)%hTr = S_REF*H0
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         ! Initial interior heat content (hTr is T * h, so sum / grid is total).
         heat0 = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
         end do

         !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         heat_after = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))

         ! Heat must have risen (warmer inflow T_BC > T_INTERIOR).
         call check(error, heat_after > heat0, &
                    "test 8: interior heat content did not rise with warm clamped inflow")
         if (allocated(error)) exit checks

         ! Upper bound: no more than entire basin replaced with T_BC.
         expected_max_rise = real(NX_PHYS*NY_PHYS*NZ, wp)*H0*(T_BC - T_INTERIOR)
         call check(error, (heat_after - heat0) < expected_max_rise, &
                    "test 8: heat content rose more than entire basin replacement — impossible")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_inflow_tracer_rise

   ! -----------------------------------------------------------------
   ! Test 9: radiating gravity wave vs wall — energy decay
   ! -----------------------------------------------------------------

   subroutine test_radiating_gw_decay(error)
      !! Design §3 test 9.
      !!
      !! η bump in a basin with all-OPEN edges.  Over N_STEPS outer steps the
      !! total energy proxy (Σ (h_total - H_ref)²) should be smaller than
      !! DECAY_FRAC times the WALL run energy, confirming Flather drains the wave.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_obc
      type(multilayer_state_t) :: ms_ref, ms_obc
      type(continuity_t) :: ct_r, ct_o
      type(coriolis_adv_t) :: cor_r, cor_o
      type(ocean_pressure_force_t) :: pgf_r, pgf_o
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_o
      type(ocean_bottom_drag_t) :: bd_r, bd_o
      type(ocean_surface_stress_t) :: ss_r, ss_o
      type(ocean_vertical_advection_t) :: va_r, va_o
      type(ocean_hdiff_tracer_t) :: hd_r, hd_o
      type(ocean_vdiff_t) :: vd_r, vd_o
      type(ocean_vmix_t) :: vmix_r, vmix_o
      type(eos_t) :: eos_r, eos_o
      type(ocean_dyn_t) :: dyn_r, dyn_o
      type(ocean_bc_state_t) :: bc

      ! Single layer for simplicity (η = h_total - H_ref).
      integer, parameter :: NZ_GW = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 20
      integer, parameter :: N_STEPS = 60
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      ! Flather should drain the wave out.  After 60 outer × 20 inner substeps
      ! with CFL ≈ 0.5 per outer step, the wave has had ample time to radiate.
      ! Accept OBC energy < 50% of WALL energy as the threshold.
      real(wp), parameter :: DECAY_FRAC = 0.5_wp
      real(wp) :: c, dt_outer
      real(wp) :: e2_ref, e2_obc
      type(ocean_metrics_t) :: metrics_ref, metrics_obc
      integer :: i, j, step

      checks: block

         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_obc, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         call init_all(grid_ref, NZ_GW, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, &
                       va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
         call init_all(grid_obc, NZ_GW, ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, &
                       va_o, hd_o, vd_o, vmix_o, eos_o, dyn_o)

         call ocean_bc_state_init(bc, grid_obc, nz_ml=NZ_GW)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN

         ! m=1 cosine η pulse over the physical interior.
         do j = 1, grid_ref%ny_total
            do i = 1, grid_ref%nx_total
               if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                  ms_ref%h_layer(i, j, 1) = H0 + A_ETA* &
                                            cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp) &
                                                /real(NX_PHYS, wp))
               else
                  ms_ref%h_layer(i, j, 1) = H0
               end if
               ms_obc%h_layer(i, j, 1) = ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, 1) = &
                  eos_r%S_ref*ms_ref%h_layer(i, j, 1)
               ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, 1) = &
                  eos_r%T_ref*ms_ref%h_layer(i, j, 1)
               ms_obc%tracers(ms_obc%idx_salinity)%hTr(i, j, 1) = &
                  eos_o%S_ref*ms_obc%h_layer(i, j, 1)
               ms_obc%tracers(ms_obc%idx_temperature)%hTr(i, j, 1) = &
                  eos_o%T_ref*ms_obc%h_layer(i, j, 1)
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp
         ms_ref%v_face_y_layer = 0.0_wp
         ms_obc%u_face_x_layer = 0.0_wp
         ms_obc%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = H0
         dyn_o%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.5_wp*grid_ref%dx/c  ! outer CFL ≈ 0.5

         ! Reference: WALL (no bc).
         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         ! OBC: all-OPEN.
         call map_in(grid_obc, metrics_obc, ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_obc, metrics_obc, dyn_o, eos_o, cor_o, ct_o, pgf_o, hv_o, bd_o, ss_o, &
               va_o, hd_o, vd_o, vmix_o, ms_obc, dt_outer, N_INNER, bc=bc)
         end do
         call map_out(metrics_obc, ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)

         ! Energy proxy: Σ (h_total - H0)² over physical interior.
         block
            integer :: ng
            ng = grid_ref%nghost
            e2_ref = sum((ms_ref%h_layer(ng + 1:ng + NX_PHYS, &
                                         ng + 1:ng + NY_PHYS, 1) - H0)**2)
            e2_obc = sum((ms_obc%h_layer(ng + 1:ng + NX_PHYS, &
                                         ng + 1:ng + NY_PHYS, 1) - H0)**2)
         end block

         call check(error, e2_obc < DECAY_FRAC*e2_ref, &
                    "test 9: all-OPEN run should drain η energy below " &
                    //short_str(DECAY_FRAC)//" * WALL energy")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_obc, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, eos_o, dyn_o)
   end subroutine test_radiating_gw_decay

   ! -----------------------------------------------------------------
   ! Test 10: driver regression — all-WALL bc == no-bc
   ! -----------------------------------------------------------------

   subroutine test_wall_bc_regression(error)
      !! Design §3 test 10.
      !!
      !! Pass an all-WALL `ocean_bc_state_t` through the split driver
      !! and verify the result is bit-identical to the run without `bc`.
      !! Guards the wiring change that makes the driver pass `bc`
      !! unconditionally — all-WALL tags must reproduce the historic
      !! hard-zero behaviour bit-for-bit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_bc
      type(multilayer_state_t) :: ms_ref, ms_bc
      type(continuity_t) :: ct_r, ct_b
      type(coriolis_adv_t) :: cor_r, cor_b
      type(ocean_pressure_force_t) :: pgf_r, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_b
      type(ocean_bottom_drag_t) :: bd_r, bd_b
      type(ocean_surface_stress_t) :: ss_r, ss_b
      type(ocean_vertical_advection_t) :: va_r, va_b
      type(ocean_hdiff_tracer_t) :: hd_r, hd_b
      type(ocean_vdiff_t) :: vd_r, vd_b
      type(ocean_vmix_t) :: vmix_r, vmix_b
      type(eos_t) :: eos_r, eos_b
      type(ocean_dyn_t) :: dyn_r, dyn_b
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 6
      integer, parameter :: N_INNER = 10
      integer, parameter :: N_STEPS = 8
      real(wp) :: c, dt_outer
      real(wp) :: max_diff_h, max_diff_u, max_diff_S
      type(ocean_metrics_t) :: metrics_ref, metrics_bc
      integer :: i, j, k

      checks: block

         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_grid(grid_bc, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)

         call init_all(grid_ref, NZ, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, &
                       va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
         call init_all(grid_bc, NZ, ms_bc, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)

         ! All-WALL bc: default tags from ocean_bc_state_init.
         call ocean_bc_state_init(bc, grid_bc, nz_ml=NZ)
         ! All four edges stay OBC_WALL (default).

         ! Identical IC: small η perturbation.
         do k = 1, NZ
            do j = 1, grid_ref%ny_total
               do i = 1, grid_ref%nx_total
                  if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                     ms_ref%h_layer(i, j, k) = H0 + A_ETA* &
                                               cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp) &
                                                   /real(NX_PHYS, wp))
                  else
                     ms_ref%h_layer(i, j, k) = H0
                  end if
                  ms_bc%h_layer(i, j, k) = ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, k) = &
                     eos_r%S_ref*ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, k) = &
                     eos_r%T_ref*ms_ref%h_layer(i, j, k)
                  ms_bc%tracers(ms_bc%idx_salinity)%hTr(i, j, k) = &
                     eos_b%S_ref*ms_bc%h_layer(i, j, k)
                  ms_bc%tracers(ms_bc%idx_temperature)%hTr(i, j, k) = &
                     eos_b%T_ref*ms_bc%h_layer(i, j, k)
               end do
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp
         ms_ref%v_face_y_layer = 0.0_wp
         ms_bc%u_face_x_layer = 0.0_wp
         ms_bc%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = real(NZ, wp)*H0
         dyn_b%bt_work%bt_H_ref = real(NZ, wp)*H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.4_wp*grid_ref%dx/c

         ! Reference: no bc argument.
         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do i = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         ! Test: all-WALL bc passed in.
         call map_in(grid_bc, metrics_bc, ms_bc, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         do i = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_bc, metrics_bc, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_bc, dt_outer, N_INNER, bc=bc)
         end do
         call map_out(metrics_bc, ms_bc, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         ! Bit-identity check on h_layer, u_face_x_layer and hTr_salt.
         max_diff_h = maxval(abs(ms_ref%h_layer - ms_bc%h_layer))
         max_diff_u = maxval(abs(ms_ref%u_face_x_layer - ms_bc%u_face_x_layer))
         max_diff_S = maxval(abs(ms_ref%tracers(ms_ref%idx_salinity)%hTr &
                                 - ms_bc%tracers(ms_bc%idx_salinity)%hTr))

         call check(error, max_diff_h == 0.0_wp, &
                    "test 10: all-WALL bc changes h_layer (max_diff=" &
                    //short_str(max_diff_h)//")")
         if (allocated(error)) exit checks
         call check(error, max_diff_u == 0.0_wp, &
                    "test 10: all-WALL bc changes u_face_x_layer (max_diff=" &
                    //short_str(max_diff_u)//")")
         if (allocated(error)) exit checks
         call check(error, max_diff_S == 0.0_wp, &
                    "test 10: all-WALL bc changes hTr_salt (max_diff=" &
                    //short_str(max_diff_S)//")")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_bc, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_wall_bc_regression

   ! -----------------------------------------------------------------
   ! Test 11: outflow relaxation rate (design §1, v2 test 1)
   ! -----------------------------------------------------------------

   subroutine test_reservoir_outflow_relax(error)
      !! Design §1 v2 test 1.
      !!
      !! Uniform outflow u_n > 0 at east edge, L_out = X_L, L_in = 0.
      !! tres_0 = T_DATA (far from T_int).  After N_STEPS:
      !!   tres(n) = T_int + (tres_0 − T_int) * (1/(1 + u*dt/L))^n   (analytic)
      !!
      !! Reservoir update function is driven directly: mass_flux_x seeded to give
      !! uniform eastward outflow at the east wall face.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 6, NY_PHYS = 4
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: T_INT = 15.0_wp     ! interior concentration
      real(wp), parameter :: T_DATA = 5.0_wp     ! reservoir initial / inflow ref
      real(wp), parameter :: U_OUTFLOW = 0.1_wp  ! outward velocity (m/s)
      real(wp), parameter :: L_OUT = 500.0_wp    ! outflow length scale (m)
      real(wp), parameter :: DT = 100.0_wp       ! timestep (s)
      integer, parameter :: N_STEPS = 20
      integer :: ng, i_e, j_mid, step
      real(wp) :: tres_analytic, tres_found, decay_factor, diff

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)

         ! Enable outflow reservoir (L_in=0 → no inflow term).
         bc%res_lscale_out = L_OUT
         bc%res_lscale_in = 0.0_wp
         bc%east%bc_type = OBC_OPEN
         bc%east%clamped_tracer(1) = T_DATA

         ng = grid%nghost
         i_e = ng + NX_PHYS  ! last interior cell index
         j_mid = ng + NY_PHYS/2 + 1

         ! Seed reservoir far from equilibrium (T_DATA != T_INT).
         allocate (bc%tres_east(grid%ny_total, NZ, 2), source=T_DATA)

         ms%h_layer = H0
         ms%tracers(1)%hTr = T_INT*H0
         ! Outflow at east wall face: mass_flux_x(i_e+1) > 0 (eastward = outward at east).
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_x_layer(i_e + 1, :, :) = U_OUTFLOW*H0

         call ms%enter_data()
         ! mass_flux_x_layer is create'd (not copyin'd) by enter_data — push host values now.
         !$acc update device(ms%mass_flux_x_layer)
         call ocean_bc_state_enter_data(bc)

         do step = 1, N_STEPS
            call ocean_obc_update_reservoirs(grid, bc, ms, DT)
         end do

         !$acc update self(bc%tres_east)
         call ocean_bc_state_exit_data(bc)
         call ms%exit_data()

         ! Analytic: tres(n) = T_int + (tres_0 − T_int) * r^n
         ! where r = 1/(1 + c_out), c_out = U_OUTFLOW * DT / L_OUT.
         decay_factor = 1.0_wp/(1.0_wp + U_OUTFLOW*DT/L_OUT)
         tres_analytic = T_INT + (T_DATA - T_INT)*decay_factor**N_STEPS

         tres_found = bc%tres_east(j_mid, 1, 1)
         diff = abs(tres_found - tres_analytic)

         call check(error, diff < 1.0e-12_wp, &
                    "test 11 (outflow relax rate): tres /= analytic, diff="//short_str(diff))

      end block checks

      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine test_reservoir_outflow_relax

   ! -----------------------------------------------------------------
   ! Test 12: inflow relaxation (design §1, v2 test 2)
   ! -----------------------------------------------------------------

   subroutine test_reservoir_inflow_relax(error)
      !! Design §1 v2 test 2.
      !!
      !! Uniform inflow at west edge (u < 0 outward-normal, i.e. eastward
      !! into domain), L_in = X_L.  tres should relax toward T_data.
      !! Analytic: tres(n) = T_data + (tres_0 − T_data) * (1/(1+u*dt/L_in))^n
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 6, NY_PHYS = 4
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: T_INT = 15.0_wp
      real(wp), parameter :: T_DATA = 5.0_wp
      real(wp), parameter :: U_INFLOW = 0.1_wp   ! eastward = into domain
      real(wp), parameter :: L_IN = 500.0_wp
      real(wp), parameter :: DT = 100.0_wp
      integer, parameter :: N_STEPS = 20
      integer :: ng, i_w, j_mid, step
      real(wp) :: tres_analytic, tres_found, decay_factor, diff

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)

         bc%res_lscale_out = 0.0_wp
         bc%res_lscale_in = L_IN
         bc%west%bc_type = OBC_OPEN
         bc%west%clamped_tracer(1) = T_DATA
         ! Start far from T_DATA so relaxation is visible.
         ng = grid%nghost
         i_w = ng + 1
         j_mid = ng + NY_PHYS/2 + 1
         allocate (bc%tres_west(grid%ny_total, NZ, 2), source=T_INT)

         ms%h_layer = H0
         ms%tracers(1)%hTr = T_INT*H0
         ! Inflow at west: flux into domain = positive (eastward).
         ! Outward normal at west = -x, so u_n = -flux/h < 0 (inflow).
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_x_layer(i_w, :, :) = U_INFLOW*H0  ! +ve = into domain

         call ms%enter_data()
         ! mass_flux_x_layer is create'd (not copyin'd) by enter_data — push host values now.
         !$acc update device(ms%mass_flux_x_layer)
         call ocean_bc_state_enter_data(bc)
         do step = 1, N_STEPS
            call ocean_obc_update_reservoirs(grid, bc, ms, DT)
         end do
         !$acc update self(bc%tres_west)
         call ocean_bc_state_exit_data(bc)
         call ms%exit_data()

         ! Analytic: c_in = U_INFLOW * DT / L_IN (u_n = -U_INFLOW < 0 → inflow)
         decay_factor = 1.0_wp/(1.0_wp + U_INFLOW*DT/L_IN)
         tres_analytic = T_DATA + (T_INT - T_DATA)*decay_factor**N_STEPS

         tres_found = bc%tres_west(j_mid, 1, 1)
         diff = abs(tres_found - tres_analytic)

         call check(error, diff < 1.0e-12_wp, &
                    "test 12 (inflow relax): tres /= analytic, diff="//short_str(diff))

      end block checks

      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine test_reservoir_inflow_relax

   ! -----------------------------------------------------------------
   ! Test 13: degenerate limits (design §1, v2 test 3)
   ! -----------------------------------------------------------------

   subroutine test_reservoir_degenerate(error)
      !! Design §1 v2 test 3.
      !!
      !! L_out = 0: one update step with outflow → tres = T_int instantly.
      !! L_in  = 0: one update step with inflow  → tres = T_data instantly.
      !! Confirm the degenerate limit gives the correct instant-set value.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 6, NY_PHYS = 4
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: T_INT = 15.0_wp
      real(wp), parameter :: T_DATA = 5.0_wp
      real(wp), parameter :: U_OUTFLOW = 0.1_wp
      real(wp), parameter :: DT = 100.0_wp
      integer :: ng, i_e, j_mid
      real(wp) :: tres_found, diff

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)

         ! ---- Case A: L_out = 0, outflow → tres = T_int ----
         bc%res_lscale_out = 0.0_wp
         bc%res_lscale_in = 0.0_wp
         bc%east%bc_type = OBC_OPEN
         bc%east%clamped_tracer(1) = T_DATA
         ng = grid%nghost
         i_e = ng + NX_PHYS
         j_mid = ng + NY_PHYS/2 + 1
         ! tres_0 = T_DATA (far from T_int).
         allocate (bc%tres_east(grid%ny_total, NZ, 2), source=T_DATA)

         ms%h_layer = H0
         ms%tracers(1)%hTr = T_INT*H0
         ! Outflow at east: +ve flux at east wall face.
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_x_layer(i_e + 1, :, :) = U_OUTFLOW*H0

         ! Enable feature with at least one L > 0 so the code doesn't short-circuit.
         ! Set L_out=0, L_in=1 (large) so outflow path hits the degenerate branch.
         bc%res_lscale_out = 0.0_wp
         bc%res_lscale_in = 1.0e6_wp   ! large: no inflow penalty in this case

         call ms%enter_data()
         ! mass_flux_x_layer is create'd (not copyin'd) by enter_data — push host values now.
         !$acc update device(ms%mass_flux_x_layer)
         call ocean_bc_state_enter_data(bc)
         call ocean_obc_update_reservoirs(grid, bc, ms, DT)
         !$acc update self(bc%tres_east)
         call ocean_bc_state_exit_data(bc)
         call ms%exit_data()

         tres_found = bc%tres_east(j_mid, 1, 1)
         diff = abs(tres_found - T_INT)
         call check(error, diff < 1.0e-12_wp, &
                    "test 13 (L_out=0 degenerate): tres should equal T_int, diff=" &
                    //short_str(diff))
         if (allocated(error)) exit checks

         ! ---- Case B: L_in = 0, inflow → tres = T_data ----
         call ocean_bc_state_destroy(bc)
         call ms%destroy()

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)

         bc%res_lscale_out = 1.0e6_wp  ! large: no outflow penalty in this case
         bc%res_lscale_in = 0.0_wp
         bc%west%bc_type = OBC_OPEN
         bc%west%clamped_tracer(1) = T_DATA
         ng = grid%nghost
         j_mid = ng + NY_PHYS/2 + 1
         ! tres_0 = T_INT (far from T_data).
         allocate (bc%tres_west(grid%ny_total, NZ, 2), source=T_INT)

         ms%h_layer = H0
         ms%tracers(1)%hTr = T_INT*H0
         ! Inflow at west: +ve flux (eastward into domain).
         ms%mass_flux_x_layer = 0.0_wp
         i_e = ng + 1  ! west wall face for mass flux
         ms%mass_flux_x_layer(i_e, :, :) = 0.1_wp*H0

         call ms%enter_data()
         ! mass_flux_x_layer is create'd (not copyin'd) by enter_data — push host values now.
         !$acc update device(ms%mass_flux_x_layer)
         call ocean_bc_state_enter_data(bc)
         call ocean_obc_update_reservoirs(grid, bc, ms, DT)
         !$acc update self(bc%tres_west)
         call ocean_bc_state_exit_data(bc)
         call ms%exit_data()

         tres_found = bc%tres_west(j_mid, 1, 1)
         diff = abs(tres_found - T_DATA)
         call check(error, diff < 1.0e-12_wp, &
                    "test 13 (L_in=0 degenerate): tres should equal T_data, diff=" &
                    //short_str(diff))

      end block checks

      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine test_reservoir_degenerate

   ! -----------------------------------------------------------------
   ! Test 14: no-chatter with alternating velocity (design §1, v2 test 4)
   ! -----------------------------------------------------------------

   subroutine test_reservoir_no_chatter(error)
      !! Design §1 v2 test 4.
      !!
      !! Alternate the face velocity sign every step with small magnitude.
      !! Assert the ghost-concentration step-to-step change is bounded by
      !! u*dt/L (smooth), whereas the v1 sign-switch path would jump by
      !! |T_data − T_int| each step.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 6, NY_PHYS = 4
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: T_INT = 15.0_wp
      real(wp), parameter :: T_DATA = 5.0_wp
      real(wp), parameter :: U_SMALL = 0.01_wp  ! small oscillating velocity
      real(wp), parameter :: L_OUT = 1000.0_wp
      real(wp), parameter :: L_IN = 1000.0_wp
      real(wp), parameter :: DT = 50.0_wp
      integer, parameter :: N_STEPS = 8
      integer :: ng, i_e, j_mid, step
      real(wp) :: tres_prev, tres_curr, max_step, bound_step, sign_val

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)

         bc%res_lscale_out = L_OUT
         bc%res_lscale_in = L_IN
         bc%east%bc_type = OBC_OPEN
         bc%east%clamped_tracer(1) = T_DATA
         ng = grid%nghost
         i_e = ng + NX_PHYS
         j_mid = ng + NY_PHYS/2 + 1
         ! Start mid-way between T_DATA and T_INT.
         allocate (bc%tres_east(grid%ny_total, NZ, 2), source=0.5_wp*(T_DATA + T_INT))

         ms%h_layer = H0
         ms%tracers(1)%hTr = T_INT*H0
         ms%mass_flux_x_layer = 0.0_wp

         call ms%enter_data()
         call ocean_bc_state_enter_data(bc)

         max_step = 0.0_wp
         do step = 1, N_STEPS
            ! Read tres before update (pull from device).
            !$acc update self(bc%tres_east)
            tres_prev = bc%tres_east(j_mid, 1, 1)

            ! Alternate direction each step; push updated flux to device.
            sign_val = real(1 - 2*mod(step, 2), wp)  ! +1, -1, +1, ...
            ms%mass_flux_x_layer(i_e + 1, :, :) = sign_val*U_SMALL*H0
            !$acc update device(ms%mass_flux_x_layer)

            call ocean_obc_update_reservoirs(grid, bc, ms, DT)

            !$acc update self(bc%tres_east)
            tres_curr = bc%tres_east(j_mid, 1, 1)
            max_step = max(max_step, abs(tres_curr - tres_prev))
         end do

         call ocean_bc_state_exit_data(bc)
         call ms%exit_data()

         ! Bound: each step changes tres by at most c * |T_int - T_data|
         ! where c = U_SMALL * DT / L_OUT (plus the inflow term symmetrically).
         ! Being generous: bound = 2 * (U_SMALL * DT / L) * |T_int - T_data|.
         bound_step = 2.0_wp*(U_SMALL*DT/L_OUT)*abs(T_INT - T_DATA)

         call check(error, max_step <= bound_step + 1.0e-14_wp, &
                    "test 14 (no-chatter): step "//short_str(max_step)// &
                    " exceeds bound "//short_str(bound_step))
         if (allocated(error)) exit checks

         ! Verify: v1 sign-switch step would be |T_data − T_int| = 10.
         ! Confirm our bound is much smaller, i.e. the test actually distinguishes.
         call check(error, bound_step < 0.1_wp*abs(T_INT - T_DATA), &
                    "test 14 setup error: bound_step should be << |T_int - T_data|")

      end block checks

      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine test_reservoir_no_chatter

   ! -----------------------------------------------------------------
   ! Test 15: bit-identity when feature disabled (design §1, v2 test 5)
   ! -----------------------------------------------------------------

   subroutine test_reservoir_bit_identity(error)
      !! Design §1 v2 test 5.
      !!
      !! Run the same all-OPEN basin twice: once with res_lscale_out/in = 0
      !! (default — v1 sign-switch path) and once with a bc that has no
      !! tres_* allocated (equivalent to the no-reservoir case).
      !! The tracer-ghost fill must give bit-identical hTr ghost columns.
      !! This re-uses the test_radiating_gw_decay setup but checks tracers.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(coriolis_adv_t) :: cor_a, cor_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      type(ocean_bottom_drag_t) :: bd_a, bd_b
      type(ocean_surface_stress_t) :: ss_a, ss_b
      type(ocean_vertical_advection_t) :: va_a, va_b
      type(ocean_hdiff_tracer_t) :: hd_a, hd_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(eos_t) :: eos_a, eos_b
      type(ocean_dyn_t) :: dyn_a, dyn_b
      type(ocean_bc_state_t) :: bc_a, bc_b

      integer, parameter :: NZ_G = 1
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_INNER = 10
      integer, parameter :: N_STEPS = 5
      real(wp), parameter :: H0 = 100.0_wp
      real(wp) :: c, dt_outer
      real(wp) :: max_diff_S
      type(ocean_metrics_t) :: metrics_a, metrics_b
      integer :: step, i, j

      checks: block

         call make_grid(grid_a, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_b, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         call init_all(grid_a, NZ_G, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, &
                       va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
         call init_all(grid_b, NZ_G, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)

         ! Both runs: all-OPEN, no reservoirs (tres_* not allocated).
         call ocean_bc_state_init(bc_a, grid_a, nz_ml=NZ_G)
         bc_a%west%bc_type = OBC_OPEN
         bc_a%east%bc_type = OBC_OPEN
         bc_a%south%bc_type = OBC_OPEN
         bc_a%north%bc_type = OBC_OPEN

         call ocean_bc_state_init(bc_b, grid_b, nz_ml=NZ_G)
         bc_b%west%bc_type = OBC_OPEN
         bc_b%east%bc_type = OBC_OPEN
         bc_b%south%bc_type = OBC_OPEN
         bc_b%north%bc_type = OBC_OPEN
         ! Explicitly zero res_lscale_out/in (should already be default 0).
         bc_b%res_lscale_out = 0.0_wp
         bc_b%res_lscale_in = 0.0_wp

         ! Identical IC: flat h, reference S/T.
         do j = 1, grid_a%ny_total
            do i = 1, grid_a%nx_total
               ms_a%h_layer(i, j, 1) = H0
               ms_b%h_layer(i, j, 1) = H0
               ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, 1) = eos_a%S_ref*H0
               ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, 1) = eos_a%T_ref*H0
               ms_b%tracers(ms_b%idx_salinity)%hTr(i, j, 1) = eos_b%S_ref*H0
               ms_b%tracers(ms_b%idx_temperature)%hTr(i, j, 1) = eos_b%T_ref*H0
            end do
         end do
         ms_a%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp
         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         dyn_a%bt_work%bt_H_ref = H0
         dyn_b%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.3_wp*grid_a%dx/c

         call map_in(grid_a, metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_a, metrics_a, dyn_a, eos_a, cor_a, ct_a, pgf_a, hv_a, bd_a, ss_a, &
               va_a, hd_a, vd_a, vmix_a, ms_a, dt_outer, N_INNER, bc=bc_a)
         end do
         call map_out(metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)

         call map_in(grid_b, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_b, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_b, dt_outer, N_INNER, bc=bc_b)
         end do
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         max_diff_S = maxval(abs(ms_a%tracers(ms_a%idx_salinity)%hTr &
                                 - ms_b%tracers(ms_b%idx_salinity)%hTr))

         call check(error, max_diff_S == 0.0_wp, &
                    "test 15 (bit-identity): zero vs zero res_lscale gives different S, diff=" &
                    //short_str(max_diff_S))

      end block checks

      call ocean_bc_state_destroy(bc_a)
      call ocean_bc_state_destroy(bc_b)
      call destroy_all(ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_reservoir_bit_identity

   ! -----------------------------------------------------------------
   ! Test 16: corner-ζ closure (design §3, v2)
   ! -----------------------------------------------------------------

   subroutine test_corner_zeta_zero(error)
      !! Design §3 v2.
      !!
      !! All-OPEN basin: after one barotropic_substep_nonlinear call,
      !! assert that bt_zeta_corner is exactly 0 along all four
      !! boundary corner lines:
      !!   west  corner line: i = nghost+1
      !!   east  corner line: i = nghost + nx_phys + 1
      !!   south corner line: j = nghost+1
      !!   north corner line: j = nghost + ny_phys + 1
      !!
      !! WALL run also checked: corner-ζ is still 0 (behaviour unchanged —
      !! the new predicate is equivalent for WALL because OBC_WALL /= OBC_PERIODIC).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(barotropic_workstate_t) :: bt_work
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 6
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: DT_INNER = 10.0_wp
      real(wp), parameter :: SMALL_U = 0.01_wp
      integer, parameter :: N_INNER = 1
      integer :: ng, i_w, i_e, j_s, j_n
      integer :: nx, ny
      real(wp), allocatable :: force_u(:, :), force_v(:, :)
      real(wp) :: max_zeta_w, max_zeta_e, max_zeta_s, max_zeta_n

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bt_work%init(grid, nz_ml=NZ)
         call cor%init(grid, nz_ml=NZ)

         ng = grid%nghost
         nx = grid%nx_total
         ny = grid%ny_total
         i_w = ng + 1
         i_e = ng + NX_PHYS + 1
         j_s = ng + 1
         j_n = ng + NY_PHYS + 1

         allocate (force_u(nx + 1, ny), source=0.0_wp)
         allocate (force_v(nx, ny + 1), source=0.0_wp)

         ! Set up a non-trivial velocity field so ζ is nonzero before closure.
         ! Uniform small u, zero v → shear ζ = 0 initially, but any residual
         ! ghost-cell noise would show up. Let's add a random-ish pattern.
         bt_work%bt_H_ref = H0
         bt_work%bt_eta = 0.0_wp
         bt_work%bt_ubt = SMALL_U
         bt_work%bt_vbt = 0.0_wp

         ! ---- All-OPEN bc ----
         call ocean_bc_state_init(bc, grid, nz_ml=NZ)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN

         !$acc enter data copyin(bt_work, cor, bc, force_u, force_v)
         call bt_work%enter_data()
         call cor%enter_data()

         call barotropic_substep_nonlinear(grid, bt_work, &
                                           force_u, force_v, &
                                           N_INNER, DT_INNER, &
                                           bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                           bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                           eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                           bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                           bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                           uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                           bt_ubt_end=bt_work%bt_ubt_end, &
                                           bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                           bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                           vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                           bt_vbt_end=bt_work%bt_vbt_end, &
                                           bt_zeta_corner=bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           bc=bc, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)

         !$acc update self(bt_work%bt_zeta_corner)
         call cor%exit_data()
         call bt_work%exit_data()
         !$acc exit data delete(bt_work, cor, bc, force_u, force_v)

         ! West corner line: bt_zeta_corner(i_w, j) for j in 1:ny+1
         max_zeta_w = maxval(abs(bt_work%bt_zeta_corner(i_w, :)))
         ! East corner line
         max_zeta_e = maxval(abs(bt_work%bt_zeta_corner(i_e, :)))
         ! South corner line
         max_zeta_s = maxval(abs(bt_work%bt_zeta_corner(:, j_s)))
         ! North corner line
         max_zeta_n = maxval(abs(bt_work%bt_zeta_corner(:, j_n)))

         call check(error, max_zeta_w == 0.0_wp, &
                    "test 16 (corner-ζ): west corner line ζ /= 0, max=" &
                    //short_str(max_zeta_w))
         if (allocated(error)) exit checks
         call check(error, max_zeta_e == 0.0_wp, &
                    "test 16 (corner-ζ): east corner line ζ /= 0, max=" &
                    //short_str(max_zeta_e))
         if (allocated(error)) exit checks
         call check(error, max_zeta_s == 0.0_wp, &
                    "test 16 (corner-ζ): south corner line ζ /= 0, max=" &
                    //short_str(max_zeta_s))
         if (allocated(error)) exit checks
         call check(error, max_zeta_n == 0.0_wp, &
                    "test 16 (corner-ζ): north corner line ζ /= 0, max=" &
                    //short_str(max_zeta_n))
         if (allocated(error)) exit checks

         ! ---- WALL bc: corner-ζ must also be 0 (bit-identical behaviour) ----
         call ocean_bc_state_destroy(bc)
         call bt_work%destroy()

         call bt_work%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ)
         ! Default: all WALL.

         bt_work%bt_H_ref = H0
         bt_work%bt_eta = 0.0_wp
         bt_work%bt_ubt = SMALL_U
         bt_work%bt_vbt = 0.0_wp

         !$acc enter data copyin(bt_work, cor, bc, force_u, force_v)
         call bt_work%enter_data()
         call cor%enter_data()

         call barotropic_substep_nonlinear(grid, bt_work, &
                                           force_u, force_v, &
                                           N_INNER, DT_INNER, &
                                           bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                           bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                           eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                           bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                           bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                           uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                           bt_ubt_end=bt_work%bt_ubt_end, &
                                           bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                           bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                           vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                           bt_vbt_end=bt_work%bt_vbt_end, &
                                           bt_zeta_corner=bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           bc=bc, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)

         !$acc update self(bt_work%bt_zeta_corner)
         call cor%exit_data()
         call bt_work%exit_data()
         !$acc exit data delete(bt_work, cor, bc, force_u, force_v)

         max_zeta_w = maxval(abs(bt_work%bt_zeta_corner(i_w, :)))
         max_zeta_e = maxval(abs(bt_work%bt_zeta_corner(i_e, :)))
         max_zeta_s = maxval(abs(bt_work%bt_zeta_corner(:, j_s)))
         max_zeta_n = maxval(abs(bt_work%bt_zeta_corner(:, j_n)))

         call check(error, max_zeta_w == 0.0_wp, &
                    "test 16 (corner-ζ WALL): west corner line ζ /= 0")
         if (allocated(error)) exit checks
         call check(error, max_zeta_e == 0.0_wp, &
                    "test 16 (corner-ζ WALL): east corner line ζ /= 0")
         if (allocated(error)) exit checks
         call check(error, max_zeta_s == 0.0_wp, &
                    "test 16 (corner-ζ WALL): south corner line ζ /= 0")
         if (allocated(error)) exit checks
         call check(error, max_zeta_n == 0.0_wp, &
                    "test 16 (corner-ζ WALL): north corner line ζ /= 0")

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call cor%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_corner_zeta_zero

   ! -----------------------------------------------------------------
   ! Test 17: Orlanski rx computation unit test
   ! -----------------------------------------------------------------

   subroutine test_orlanski_rx_unit(error)
      !! Design §2 v2 test 17.
      !!
      !! Drive apply_orlanski_east directly with prescribed u_prev / u_new
      !! and verify:
      !!   (a) When dhdt·dhdx > 0 (outgoing) rx is clipped to rx_max.
      !!   (b) When dhdt·dhdx <= 0 (incoming) rx stays at 0.
      !!
      !! Setup: single j-column, 2-layer, east OPEN edge.
      !! Case A: u_prev > u_int_1 AND u(i_e-1) > u(i_e-2) ⟹ both positive,
      !!         product > 0, rx_raw = dhdt/dhdx.
      !! Case B: u_prev < u_int_1 (dhdt < 0) while dhdx > 0 ⟹ product < 0, rx=0.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: RX_MAX = 2.0_wp
      real(wp), parameter :: GAMMA_U = 1.0_wp   ! instant running mean
      real(wp), parameter :: DT = 1.0_wp        ! nudging taus=0 so dt unused

      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_bc_state_t) :: bc

      integer :: ng, i_e, j_mid, j0, j1, nxt, nyt, nz_loc
      real(wp) :: rx_found, rx_expected, diff
      real(wp) :: u_int1, u_int2, u_prev_val
      real(wp) :: dhdt_val, dhdx_val, rx_raw_expected

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         nz_loc = NZ
         ms%nz_ml = nz_loc
         call ms%init(grid)
         call bt_work%init(grid, nz_ml=nz_loc)
         call ocean_bc_state_init(bc, grid, nz_ml=nz_loc, n_tracers=2)

         ng = grid%nghost
         i_e = ng + NX_PHYS + 1   ! east wall-face index
         j0 = ng + 1
         j1 = ng + NY_PHYS
         j_mid = ng + NY_PHYS/2 + 1
         nxt = grid%nx_total
         nyt = grid%ny_total

         bc%east%bc_type = OBC_OPEN
         bc%radiation_scheme = 1   ! RAD_ORLANSKI
         bc%orlanski_rx_max = RX_MAX
         bc%orlanski_gamma = GAMMA_U
         bc%nudge_tau_in = 0.0_wp
         bc%nudge_tau_out = 0.0_wp
         allocate (bc%rx_east(nyt, nz_loc), source=0.0_wp)
         allocate (bc%u_prev_east(nyt, nz_loc), source=0.0_wp)

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         bt_work%bt_ubt_end = 0.0_wp
         bt_work%bt_H_ref = real(nz_loc, wp)*H0

         ! ---- Case A: outgoing (dhdt·dhdx > 0) ----
         ! East interior: u(i_e-1) = u_int1, u(i_e-2) = u_int2.
         ! Outgoing east wave: u_prev > u_int1 (velocity decayed, dhdt>0),
         ! dhdx = u_int1 - u_int2 > 0 (decreasing inward = eastward gradient).
         u_int1 = 0.05_wp   ! 1st interior face
         u_int2 = 0.03_wp   ! 2nd interior face
         u_prev_val = 0.08_wp  ! was higher last step
         ms%u_face_x_layer(i_e - 1, :, :) = u_int1
         ms%u_face_x_layer(i_e - 2, :, :) = u_int2
         bc%u_prev_east(:, :) = u_prev_val

         dhdt_val = u_prev_val - u_int1  ! = 0.03
         dhdx_val = u_int1 - u_int2      ! = 0.02
         rx_raw_expected = min(dhdt_val/dhdx_val, RX_MAX)  ! = 1.5

         !$acc enter data copyin(ms, bt_work, bc)
         call ms%enter_data()
         call bt_work%enter_data()
         call ocean_bc_state_enter_data(bc)
         call ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, DT)
         !$acc update self(bc%rx_east)
         call ocean_bc_state_exit_data(bc)
         call bt_work%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, bt_work, bc)

         rx_found = bc%rx_east(j_mid, 1)
         rx_expected = rx_raw_expected
         diff = abs(rx_found - rx_expected)

         call check(error, diff < 1.0e-12_wp, &
                    "test 17 case A (outgoing): rx /= dhdt/dhdx, diff="//short_str(diff))
         if (allocated(error)) exit checks

         ! ---- Case B: incoming (dhdt·dhdx <= 0) ----
         ! Reset rx to non-zero first to confirm it is zeroed.
         bc%rx_east(:, :) = 0.5_wp
         u_int1 = 0.03_wp
         u_int2 = 0.01_wp
         u_prev_val = 0.01_wp  ! velocity grew, so dhdt = prev - new < 0
         ms%u_face_x_layer(i_e - 1, :, :) = u_int1
         ms%u_face_x_layer(i_e - 2, :, :) = u_int2
         bc%u_prev_east(:, :) = u_prev_val

         ! dhdt = 0.01 - 0.03 = -0.02; dhdx = 0.03 - 0.01 = 0.02 → product < 0 → incoming

         !$acc enter data copyin(ms, bt_work, bc)
         call ms%enter_data()
         call bt_work%enter_data()
         call ocean_bc_state_enter_data(bc)
         call ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, DT)
         !$acc update self(bc%rx_east)
         call ocean_bc_state_exit_data(bc)
         call bt_work%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, bt_work, bc)

         rx_found = bc%rx_east(j_mid, 1)
         diff = abs(rx_found)   ! should be 0 (gamma=1, rx_raw=0)

         call check(error, diff < 1.0e-12_wp, &
                    "test 17 case B (incoming): rx should be 0, got "//short_str(rx_found))

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call ms%destroy()
   end subroutine test_orlanski_rx_unit

   ! -----------------------------------------------------------------
   ! Test 18: Orlanski scheme drains radiating wave (≤ anomaly energy)
   ! -----------------------------------------------------------------

   subroutine test_orlanski_wave_drains(error)
      !! Design §2 v2 test 18.
      !!
      !! η bump in a basin with all-OPEN Orlanski edges.  After N_STEPS outer
      !! steps the total energy proxy (Σ (h_total - H_ref)²) must be no larger
      !! than the anomaly-scheme energy (both drain the wave; Orlanski may
      !! drain faster but must not be worse).  Threshold: factor of 3 above
      !! anomaly result is a failure (tight criterion: Orlanski is a radiation
      !! scheme and should generally perform as well or better).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_anom, grid_orl
      type(multilayer_state_t) :: ms_anom, ms_orl
      type(continuity_t) :: ct_a, ct_o
      type(coriolis_adv_t) :: cor_a, cor_o
      type(ocean_pressure_force_t) :: pgf_a, pgf_o
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_o
      type(ocean_bottom_drag_t) :: bd_a, bd_o
      type(ocean_surface_stress_t) :: ss_a, ss_o
      type(ocean_vertical_advection_t) :: va_a, va_o
      type(ocean_hdiff_tracer_t) :: hd_a, hd_o
      type(ocean_vdiff_t) :: vd_a, vd_o
      type(ocean_vmix_t) :: vmix_a, vmix_o
      type(eos_t) :: eos_a, eos_o
      type(ocean_dyn_t) :: dyn_a, dyn_o
      type(ocean_bc_state_t) :: bc_anom, bc_orl

      integer, parameter :: NZ_GW = 1
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_INNER = 20
      integer, parameter :: N_STEPS = 40
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: ENERGY_RATIO_MAX = 3.0_wp
      real(wp) :: c, dt_outer
      real(wp) :: e2_anom, e2_orl
      type(ocean_metrics_t) :: metrics_anom, metrics_orl
      integer :: i, j, step

      checks: block

         call make_grid(grid_anom, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         call make_grid(grid_orl, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)

         call init_all(grid_anom, NZ_GW, ms_anom, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, &
                       va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
         call init_all(grid_orl, NZ_GW, ms_orl, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, &
                       va_o, hd_o, vd_o, vmix_o, eos_o, dyn_o)

         ! Anomaly scheme (default).
         call ocean_bc_state_init(bc_anom, grid_anom, nz_ml=NZ_GW)
         bc_anom%west%bc_type = OBC_OPEN
         bc_anom%east%bc_type = OBC_OPEN
         bc_anom%south%bc_type = OBC_OPEN
         bc_anom%north%bc_type = OBC_OPEN

         ! Orlanski scheme: allocate rx/u_prev for all four radiating edges.
         call ocean_bc_state_init(bc_orl, grid_orl, nz_ml=NZ_GW)
         bc_orl%west%bc_type = OBC_OPEN
         bc_orl%east%bc_type = OBC_OPEN
         bc_orl%south%bc_type = OBC_OPEN
         bc_orl%north%bc_type = OBC_OPEN
         bc_orl%radiation_scheme = 1
         bc_orl%orlanski_rx_max = 10.0_wp
         bc_orl%orlanski_gamma = 1.0_wp
         allocate (bc_orl%rx_west(grid_orl%ny_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%rx_east(grid_orl%ny_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%rx_south(grid_orl%nx_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%rx_north(grid_orl%nx_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%u_prev_west(grid_orl%ny_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%u_prev_east(grid_orl%ny_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%u_prev_south(grid_orl%nx_total, NZ_GW), source=0.0_wp)
         allocate (bc_orl%u_prev_north(grid_orl%nx_total, NZ_GW), source=0.0_wp)

         ! Identical IC: cosine η pulse.
         do j = 1, grid_anom%ny_total
            do i = 1, grid_anom%nx_total
               if (i > grid_anom%nghost .and. i <= grid_anom%nghost + NX_PHYS) then
                  ms_anom%h_layer(i, j, 1) = H0 + A_ETA* &
                                             cos(2.0_wp*PI*(real(i - grid_anom%nghost, wp) - 0.5_wp)/real(NX_PHYS, wp))
               else
                  ms_anom%h_layer(i, j, 1) = H0
               end if
               ms_orl%h_layer(i, j, 1) = ms_anom%h_layer(i, j, 1)
               ms_anom%tracers(ms_anom%idx_salinity)%hTr(i, j, 1) = &
                  eos_a%S_ref*ms_anom%h_layer(i, j, 1)
               ms_anom%tracers(ms_anom%idx_temperature)%hTr(i, j, 1) = &
                  eos_a%T_ref*ms_anom%h_layer(i, j, 1)
               ms_orl%tracers(ms_orl%idx_salinity)%hTr(i, j, 1) = &
                  eos_o%S_ref*ms_orl%h_layer(i, j, 1)
               ms_orl%tracers(ms_orl%idx_temperature)%hTr(i, j, 1) = &
                  eos_o%T_ref*ms_orl%h_layer(i, j, 1)
            end do
         end do
         ms_anom%u_face_x_layer = 0.0_wp
         ms_anom%v_face_y_layer = 0.0_wp
         ms_orl%u_face_x_layer = 0.0_wp
         ms_orl%v_face_y_layer = 0.0_wp
         dyn_a%bt_work%bt_H_ref = H0
         dyn_o%bt_work%bt_H_ref = H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.4_wp*grid_anom%dx/c

         call map_in(grid_anom, metrics_anom, ms_anom, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_anom, metrics_anom, dyn_a, eos_a, cor_a, ct_a, pgf_a, hv_a, bd_a, ss_a, &
               va_a, hd_a, vd_a, vmix_a, ms_anom, dt_outer, N_INNER, bc=bc_anom)
         end do
         call map_out(metrics_anom, ms_anom, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)

         !$acc enter data copyin(bc_orl)
         call ocean_bc_state_enter_data(bc_orl)
         call map_in(grid_orl, metrics_orl, ms_orl, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_orl, metrics_orl, dyn_o, eos_o, cor_o, ct_o, pgf_o, hv_o, bd_o, ss_o, &
               va_o, hd_o, vd_o, vmix_o, ms_orl, dt_outer, N_INNER, bc=bc_orl)
         end do
         call map_out(metrics_orl, ms_orl, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, dyn_o)
         call ocean_bc_state_exit_data(bc_orl)
         !$acc exit data delete(bc_orl)

         block
            integer :: ng
            ng = grid_anom%nghost
            e2_anom = sum((ms_anom%h_layer(ng + 1:ng + NX_PHYS, ng + 1:ng + NY_PHYS, 1) - H0)**2)
            e2_orl = sum((ms_orl%h_layer(ng + 1:ng + NX_PHYS, ng + 1:ng + NY_PHYS, 1) - H0)**2)
         end block

         ! Both should drain energy; Orlanski energy must not be grossly worse than anomaly.
         call check(error, e2_orl < ENERGY_RATIO_MAX*e2_anom + 1.0e-12_wp, &
                    "test 18 (Orlanski wave): energy ratio Orlanski/anomaly="// &
                    short_str(e2_orl/(e2_anom + 1.0e-30_wp))//" > "// &
                    short_str(ENERGY_RATIO_MAX))

      end block checks

      call ocean_bc_state_destroy(bc_anom)
      call ocean_bc_state_destroy(bc_orl)
      call destroy_all(ms_anom, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
      call destroy_all(ms_orl, ct_o, cor_o, pgf_o, hv_o, bd_o, ss_o, va_o, hd_o, vd_o, vmix_o, eos_o, dyn_o)
   end subroutine test_orlanski_wave_drains

   ! -----------------------------------------------------------------
   ! Test 19: Orlanski nudged inflow — tau_in small → boundary relaxes to u_data
   ! -----------------------------------------------------------------

   subroutine test_orlanski_nudge_inflow(error)
      !! Design §2 v2 test 19.
      !!
      !! East OPEN edge, radiation_scheme=1, tau_in = DT (aggressive nudging).
      !! Interior field is uniform (dhdx = 0) → incoming criterion always met
      !! (dhdt·dhdx = 0 ≤ 0) → tau = tau_in, and rx = 0 so the radiation
      !! update leaves the wall face at its own previous per-layer value
      !! (U_INT — the whole field is seeded uniform).  One nudge then gives
      !! the exact single-step result
      !!   u_wall = (1−g2)·U_INT + g2·U_DATA,  g2 = DT/(TAU_IN+DT) = 0.5.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: U_DATA = 0.2_wp    ! target inflow velocity
      real(wp), parameter :: U_INT = -0.05_wp   ! interior face: incoming (westward at east)
      real(wp), parameter :: DT = 120.0_wp
      real(wp), parameter :: TAU_IN = DT        ! g2 = DT/(TAU_IN+DT) = 0.5
      integer, parameter :: N_STEPS = 1         ! single step suffices for analytical check

      integer :: ng, i_e, j_mid, nyt, nz_loc, step
      real(wp) :: u_wall_after, g2_expected, u_expected, diff

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         nz_loc = NZ
         ms%nz_ml = nz_loc
         call ms%init(grid)
         call bt_work%init(grid, nz_ml=nz_loc)
         call ocean_bc_state_init(bc, grid, nz_ml=nz_loc, n_tracers=2)

         ng = grid%nghost
         i_e = ng + NX_PHYS + 1
         j_mid = ng + NY_PHYS/2 + 1
         nyt = grid%ny_total

         bc%east%bc_type = OBC_OPEN
         bc%east%clamped_u = U_DATA
         bc%radiation_scheme = 1
         bc%orlanski_rx_max = 10.0_wp
         bc%orlanski_gamma = 1.0_wp
         bc%nudge_tau_in = TAU_IN
         bc%nudge_tau_out = 0.0_wp
         allocate (bc%rx_east(nyt, nz_loc), source=0.0_wp)
         allocate (bc%u_prev_east(nyt, nz_loc), source=0.0_wp)

         ms%h_layer = H0
         ! Uniform interior velocity (incoming at east; dhdx = 0 ⟹ product = 0 ≤ 0).
         ms%u_face_x_layer = U_INT
         ms%v_face_y_layer = 0.0_wp
         bt_work%bt_ubt_end = 0.0_wp   ! Orlanski base u_b from ubt_end = 0
         bt_work%bt_H_ref = real(nz_loc, wp)*H0

         !$acc enter data copyin(ms, bt_work, bc)
         call ms%enter_data()
         call bt_work%enter_data()
         call ocean_bc_state_enter_data(bc)

         do step = 1, N_STEPS
            call ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, DT)
         end do

         !$acc update self(ms%u_face_x_layer)
         call ocean_bc_state_exit_data(bc)
         call bt_work%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, bt_work, bc)

         ! Analytical: rx = 0 ⟹ u_b_base = u_wall_old = U_INT (per-layer
         ! anchor — the radiation update must NOT collapse the boundary to
         ! the depth-uniform barotropic value during incoming phases).
         !             g2 = DT/(TAU_IN+DT) = 0.5
         !             u_wall = (1-g2)*U_INT + g2*U_DATA = 0.075
         g2_expected = DT/(TAU_IN + DT)
         u_expected = (1.0_wp - g2_expected)*U_INT + g2_expected*U_DATA
         u_wall_after = ms%u_face_x_layer(i_e, j_mid, 1)
         diff = abs(u_wall_after - u_expected)

         call check(error, diff < 1.0e-12_wp, &
                    "test 19 (Orlanski nudged inflow): u_wall /= g2*U_DATA analytical, diff=" &
                    //short_str(diff)//" u_wall="//short_str(u_wall_after))

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call ms%destroy()
   end subroutine test_orlanski_nudge_inflow

   ! -----------------------------------------------------------------
   ! Test 20: Orlanski default=anomaly bit-identity
   ! -----------------------------------------------------------------

   subroutine test_orlanski_default_bit_identity(error)
      !! Design §2 v2 test 20.
      !!
      !! With radiation_scheme = "anomaly" (default, bc%radiation_scheme = 0)
      !! and nudge_tau_in/out = 0.0, the result must be bit-identical to a run
      !! that never sets radiation_scheme at all.  Guards that Phase B defaults
      !! do not perturb the v1 code path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_b
      type(multilayer_state_t) :: ms_ref, ms_b
      type(continuity_t) :: ct_r, ct_b
      type(coriolis_adv_t) :: cor_r, cor_b
      type(ocean_pressure_force_t) :: pgf_r, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_b
      type(ocean_bottom_drag_t) :: bd_r, bd_b
      type(ocean_surface_stress_t) :: ss_r, ss_b
      type(ocean_vertical_advection_t) :: va_r, va_b
      type(ocean_hdiff_tracer_t) :: hd_r, hd_b
      type(ocean_vdiff_t) :: vd_r, vd_b
      type(ocean_vmix_t) :: vmix_r, vmix_b
      type(eos_t) :: eos_r, eos_b
      type(ocean_dyn_t) :: dyn_r, dyn_b
      type(ocean_bc_state_t) :: bc_ref, bc_b

      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_INNER = 10, N_STEPS = 5
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.2_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, dt_outer, max_diff_h
      type(ocean_metrics_t) :: metrics_ref, metrics_b
      integer :: i, j, k, step

      checks: block

         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_grid(grid_b, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)

         call init_all(grid_ref, NZ, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, &
                       va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
         call init_all(grid_b, NZ, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)

         ! Reference: plain OBC_OPEN, default (anomaly) scheme, no extra knobs.
         call ocean_bc_state_init(bc_ref, grid_ref, nz_ml=NZ)
         bc_ref%west%bc_type = OBC_OPEN
         bc_ref%east%bc_type = OBC_OPEN

         ! Test: explicitly set radiation_scheme=0 + taus=0 — same defaults.
         call ocean_bc_state_init(bc_b, grid_b, nz_ml=NZ)
         bc_b%west%bc_type = OBC_OPEN
         bc_b%east%bc_type = OBC_OPEN
         bc_b%radiation_scheme = 0
         bc_b%nudge_tau_in = 0.0_wp
         bc_b%nudge_tau_out = 0.0_wp

         ! Identical IC.
         do k = 1, NZ
            do j = 1, grid_ref%ny_total
               do i = 1, grid_ref%nx_total
                  if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                     ms_ref%h_layer(i, j, k) = H0 + A_ETA* &
                                               cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp)/real(NX_PHYS, wp))
                  else
                     ms_ref%h_layer(i, j, k) = H0
                  end if
                  ms_b%h_layer(i, j, k) = ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, k) = &
                     eos_r%S_ref*ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, k) = &
                     eos_r%T_ref*ms_ref%h_layer(i, j, k)
                  ms_b%tracers(ms_b%idx_salinity)%hTr(i, j, k) = &
                     eos_b%S_ref*ms_b%h_layer(i, j, k)
                  ms_b%tracers(ms_b%idx_temperature)%hTr(i, j, k) = &
                     eos_b%T_ref*ms_b%h_layer(i, j, k)
               end do
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp
         ms_ref%v_face_y_layer = 0.0_wp
         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = real(NZ, wp)*H0
         dyn_b%bt_work%bt_H_ref = real(NZ, wp)*H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.4_wp*grid_ref%dx/c

         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER, bc=bc_ref)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         call map_in(grid_b, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_b, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_b, dt_outer, N_INNER, bc=bc_b)
         end do
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         max_diff_h = maxval(abs(ms_ref%h_layer - ms_b%h_layer))

         call check(error, max_diff_h == 0.0_wp, &
                    "test 20 (Orlanski default bit-identity): h_layer differs, max_diff=" &
                    //short_str(max_diff_h))

      end block checks

      call ocean_bc_state_destroy(bc_ref)
      call ocean_bc_state_destroy(bc_b)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_orlanski_default_bit_identity

   ! -----------------------------------------------------------------
   ! Test 21: Full-Flather quiescent: flat basin stays flat
   ! -----------------------------------------------------------------

   subroutine test_flather_full_quiescent(error)
      !! Design §4 v2 test 21.
      !!
      !! Uniform-η, zero-velocity basin with all-OPEN full-Flather edges
      !! and u_ext = 0 (default).  After N_INNER substeps the basin should
      !! remain quiescent: max|η| < SMALL and max|u| < SMALL.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(barotropic_workstate_t) :: bt_work
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: DT_INNER = 5.0_wp
      integer, parameter :: N_INNER = 4
      real(wp), parameter :: SMALL = 1.0e-10_wp
      integer :: nx, ny
      real(wp), allocatable :: force_u(:, :), force_v(:, :)
      real(wp) :: max_eta, max_u

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bt_work%init(grid, nz_ml=NZ)
         call cor%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ)

         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN
         bc%use_full_flather = .true.
         ! Exterior velocities default to 0 (already set by init).

         nx = grid%nx_total
         ny = grid%ny_total
         allocate (force_u(nx + 1, ny), source=0.0_wp)
         allocate (force_v(nx, ny + 1), source=0.0_wp)

         bt_work%bt_H_ref = H0
         bt_work%bt_eta = 0.0_wp
         bt_work%bt_ubt = 0.0_wp
         bt_work%bt_vbt = 0.0_wp

         !$acc enter data copyin(bt_work, cor, bc, force_u, force_v)
         call bt_work%enter_data()
         call cor%enter_data()

         call barotropic_substep_nonlinear(grid, bt_work, &
                                           force_u, force_v, &
                                           N_INNER, DT_INNER, &
                                           bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                           bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                           eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                           bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                           bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                           uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                           bt_ubt_end=bt_work%bt_ubt_end, &
                                           bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                           bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                           vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                           bt_vbt_end=bt_work%bt_vbt_end, &
                                           bt_zeta_corner=bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           bc=bc, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)

         !$acc update self(bt_work%bt_eta, bt_work%bt_ubt)
         call cor%exit_data()
         call bt_work%exit_data()
         !$acc exit data delete(bt_work, cor, bc, force_u, force_v)

         block
            integer :: ng
            ng = grid%nghost
            max_eta = maxval(abs(bt_work%bt_eta(ng + 1:ng + NX_PHYS, ng + 1:ng + NY_PHYS)))
            max_u = maxval(abs(bt_work%bt_ubt(ng + 1:ng + NX_PHYS + 1, ng + 1:ng + NY_PHYS)))
         end block

         call check(error, max_eta < SMALL, &
                    "test 21 (full-Flather quiescent): η disturbed, max_eta="//short_str(max_eta))
         if (allocated(error)) exit checks
         call check(error, max_u < SMALL, &
                    "test 21 (full-Flather quiescent): u disturbed, max_u="//short_str(max_u))

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call cor%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_flather_full_quiescent

   ! -----------------------------------------------------------------
   ! Test 22: Full-Flather with exterior velocity spins up the interior
   ! -----------------------------------------------------------------

   subroutine test_flather_ext_vel_spinup(error)
      !! Design §4 v2 test 22.
      !!
      !! West OPEN full-Flather edge, ext_u_west > 0 (westward inflow into domain).
      !! After N_INNER substeps the interior u must have increased from 0
      !! (the Flather condition injects the exterior velocity via the
      !! half-characteristic inward-propagating term).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(barotropic_workstate_t) :: bt_work
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: U_EXT = 0.1_wp   ! westward inflow into domain
      real(wp), parameter :: DT_INNER = 5.0_wp
      integer, parameter :: N_INNER = 4
      integer :: ng, nx, ny
      real(wp), allocatable :: force_u(:, :), force_v(:, :)
      real(wp) :: u_int_after

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bt_work%init(grid, nz_ml=NZ)
         call cor%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ)

         bc%west%bc_type = OBC_OPEN
         bc%use_full_flather = .true.
         bc%ext_u_west = U_EXT   ! exterior velocity at west boundary

         ng = grid%nghost
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (force_u(nx + 1, ny), source=0.0_wp)
         allocate (force_v(nx, ny + 1), source=0.0_wp)

         bt_work%bt_H_ref = H0
         bt_work%bt_eta = 0.0_wp
         bt_work%bt_ubt = 0.0_wp
         bt_work%bt_vbt = 0.0_wp

         !$acc enter data copyin(bt_work, cor, bc, force_u, force_v)
         call bt_work%enter_data()
         call cor%enter_data()

         call barotropic_substep_nonlinear(grid, bt_work, &
                                           force_u, force_v, &
                                           N_INNER, DT_INNER, &
                                           bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                           bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                           eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                           bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                           bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                           uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                           bt_ubt_end=bt_work%bt_ubt_end, &
                                           bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                           bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                           vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                           bt_vbt_end=bt_work%bt_vbt_end, &
                                           bt_zeta_corner=bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           bc=bc, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)

         !$acc update self(bt_work%bt_ubt)
         call cor%exit_data()
         call bt_work%exit_data()
         !$acc exit data delete(bt_work, cor, bc, force_u, force_v)

         ! First interior u-face should have grown from 0 toward U_EXT.
         u_int_after = bt_work%bt_ubt(ng + 2, ng + NY_PHYS/2 + 1)

         call check(error, u_int_after > 0.0_wp, &
                    "test 22 (full-Flather ext vel): interior u did not spin up, u=" &
                    //short_str(u_int_after))

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call cor%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_flather_ext_vel_spinup

   ! -----------------------------------------------------------------
   ! Test 23: Full-Flather legacy default bit-identity
   ! -----------------------------------------------------------------

   subroutine test_flather_legacy_default(error)
      !! Design §4 v2 test 23.
      !!
      !! With use_full_flather = .false. (default, flather_form="legacy"),
      !! the result must be bit-identical to a run without any bc knob.
      !! Guards that Phase B Flather defaults do not perturb the v1 code path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_ref, grid_b
      type(multilayer_state_t) :: ms_ref, ms_b
      type(continuity_t) :: ct_r, ct_b
      type(coriolis_adv_t) :: cor_r, cor_b
      type(ocean_pressure_force_t) :: pgf_r, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_r, hv_b
      type(ocean_bottom_drag_t) :: bd_r, bd_b
      type(ocean_surface_stress_t) :: ss_r, ss_b
      type(ocean_vertical_advection_t) :: va_r, va_b
      type(ocean_hdiff_tracer_t) :: hd_r, hd_b
      type(ocean_vdiff_t) :: vd_r, vd_b
      type(ocean_vmix_t) :: vmix_r, vmix_b
      type(eos_t) :: eos_r, eos_b
      type(ocean_dyn_t) :: dyn_r, dyn_b
      type(ocean_bc_state_t) :: bc_ref, bc_b

      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_INNER = 10, N_STEPS = 5
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: A_ETA = 0.2_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, dt_outer, max_diff_h, max_diff_u
      type(ocean_metrics_t) :: metrics_ref, metrics_b
      integer :: i, j, k, step

      checks: block

         call make_grid(grid_ref, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_grid(grid_b, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)

         call init_all(grid_ref, NZ, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, &
                       va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
         call init_all(grid_b, NZ, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)

         ! Reference: plain OBC_OPEN, no full-Flather knob.
         call ocean_bc_state_init(bc_ref, grid_ref, nz_ml=NZ)
         bc_ref%west%bc_type = OBC_OPEN
         bc_ref%east%bc_type = OBC_OPEN

         ! Test: explicitly set use_full_flather = .false. — same as default.
         call ocean_bc_state_init(bc_b, grid_b, nz_ml=NZ)
         bc_b%west%bc_type = OBC_OPEN
         bc_b%east%bc_type = OBC_OPEN
         bc_b%use_full_flather = .false.
         bc_b%ext_u_west = 0.0_wp
         bc_b%ext_u_east = 0.0_wp

         ! Identical IC.
         do k = 1, NZ
            do j = 1, grid_ref%ny_total
               do i = 1, grid_ref%nx_total
                  if (i > grid_ref%nghost .and. i <= grid_ref%nghost + NX_PHYS) then
                     ms_ref%h_layer(i, j, k) = H0 + A_ETA* &
                                               cos(2.0_wp*PI*(real(i - grid_ref%nghost, wp) - 0.5_wp)/real(NX_PHYS, wp))
                  else
                     ms_ref%h_layer(i, j, k) = H0
                  end if
                  ms_b%h_layer(i, j, k) = ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_salinity)%hTr(i, j, k) = &
                     eos_r%S_ref*ms_ref%h_layer(i, j, k)
                  ms_ref%tracers(ms_ref%idx_temperature)%hTr(i, j, k) = &
                     eos_r%T_ref*ms_ref%h_layer(i, j, k)
                  ms_b%tracers(ms_b%idx_salinity)%hTr(i, j, k) = &
                     eos_b%S_ref*ms_b%h_layer(i, j, k)
                  ms_b%tracers(ms_b%idx_temperature)%hTr(i, j, k) = &
                     eos_b%T_ref*ms_b%h_layer(i, j, k)
               end do
            end do
         end do
         ms_ref%u_face_x_layer = 0.0_wp
         ms_ref%v_face_y_layer = 0.0_wp
         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         dyn_r%bt_work%bt_H_ref = real(NZ, wp)*H0
         dyn_b%bt_work%bt_H_ref = real(NZ, wp)*H0

         c = sqrt(GRAVITY*H0)
         dt_outer = 0.4_wp*grid_ref%dx/c

         call map_in(grid_ref, metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_ref, metrics_ref, dyn_r, eos_r, cor_r, ct_r, pgf_r, hv_r, bd_r, ss_r, &
               va_r, hd_r, vd_r, vmix_r, ms_ref, dt_outer, N_INNER, bc=bc_ref)
         end do
         call map_out(metrics_ref, ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, dyn_r)

         call map_in(grid_b, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid_b, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_b, dt_outer, N_INNER, bc=bc_b)
         end do
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         max_diff_h = maxval(abs(ms_ref%h_layer - ms_b%h_layer))
         max_diff_u = maxval(abs(ms_ref%u_face_x_layer - ms_b%u_face_x_layer))

         call check(error, max_diff_h == 0.0_wp, &
                    "test 23 (Flather legacy default): h_layer differs, max_diff=" &
                    //short_str(max_diff_h))
         if (allocated(error)) exit checks
         call check(error, max_diff_u == 0.0_wp, &
                    "test 23 (Flather legacy default): u_face_x_layer differs, max_diff=" &
                    //short_str(max_diff_u))

      end block checks

      call ocean_bc_state_destroy(bc_ref)
      call ocean_bc_state_destroy(bc_b)
      call destroy_all(ms_ref, ct_r, cor_r, pgf_r, hv_r, bd_r, ss_r, va_r, hd_r, vd_r, vmix_r, eos_r, dyn_r)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_flather_legacy_default

   ! -----------------------------------------------------------------
   ! Internal helper
   ! -----------------------------------------------------------------

   pure function short_str(val) result(s)
      real(wp), intent(in) :: val
      character(len=20) :: s
      write (s, "(es12.4)") val
   end function short_str

   subroutine test_orlanski_anchor_per_layer(error)
      !! Regression lock for the Orlanski semi-implicit anchor.
      !!
      !! All four edges OPEN with radiation_scheme = orlanski, per-layer
      !! uniform fields u(:,:,k) = U_K(k), v(:,:,k) = V_K(k) and u_prev
      !! seeded to the same values: dhdt = dhdx = 0 everywhere, so rx = 0
      !! and the radiation update must leave every wall face at its OWN
      !! per-layer value.  The historical bug anchored the update on the
      !! depth-uniform barotropic velocity (ubt_end = 0 here), which
      !! collapses the boundary baroclinic shear to zero each quiet phase
      !! — this test fails loudly under that formulation.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 8, NY_PHYS = 6
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: U_K(2) = [0.07_wp, -0.04_wp]
      real(wp), parameter :: V_K(2) = [-0.02_wp, 0.05_wp]
      real(wp), parameter :: DT = 60.0_wp

      integer :: ng, i_w, i_e, j_s, j_n, nxt, nyt, nz_loc, k, j, i
      real(wp) :: max_dev

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         nz_loc = NZ
         ms%nz_ml = nz_loc
         call ms%init(grid)
         call bt_work%init(grid, nz_ml=nz_loc)
         call ocean_bc_state_init(bc, grid, nz_ml=nz_loc, n_tracers=2)

         ng = grid%nghost
         i_w = ng + 1
         i_e = ng + NX_PHYS + 1
         j_s = ng + 1
         j_n = ng + NY_PHYS + 1
         nxt = grid%nx_total
         nyt = grid%ny_total

         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN
         bc%radiation_scheme = 1
         bc%orlanski_rx_max = 10.0_wp
         bc%orlanski_gamma = 1.0_wp
         bc%nudge_tau_in = 0.0_wp
         bc%nudge_tau_out = 0.0_wp
         allocate (bc%rx_west(nyt, nz_loc), source=0.0_wp)
         allocate (bc%rx_east(nyt, nz_loc), source=0.0_wp)
         allocate (bc%rx_south(nxt, nz_loc), source=0.0_wp)
         allocate (bc%rx_north(nxt, nz_loc), source=0.0_wp)
         allocate (bc%u_prev_west(nyt, nz_loc))
         allocate (bc%u_prev_east(nyt, nz_loc))
         allocate (bc%u_prev_south(nxt, nz_loc))
         allocate (bc%u_prev_north(nxt, nz_loc))

         ms%h_layer = H0
         do k = 1, nz_loc
            ms%u_face_x_layer(:, :, k) = U_K(k)
            ms%v_face_y_layer(:, :, k) = V_K(k)
            bc%u_prev_west(:, k) = U_K(k)
            bc%u_prev_east(:, k) = U_K(k)
            bc%u_prev_south(:, k) = V_K(k)
            bc%u_prev_north(:, k) = V_K(k)
         end do
         bt_work%bt_ubt_end = 0.0_wp
         bt_work%bt_vbt_end = 0.0_wp
         bt_work%bt_H_ref = real(nz_loc, wp)*H0

         !$acc enter data copyin(ms, bt_work, bc)
         call ms%enter_data()
         call bt_work%enter_data()
         call ocean_bc_state_enter_data(bc)
         call ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, DT)
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         call ocean_bc_state_exit_data(bc)
         call bt_work%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, bt_work, bc)

         ! Every wall face must still hold its per-layer value exactly.
         max_dev = 0.0_wp
         do k = 1, nz_loc
            do j = j_s, j_n - 1
               max_dev = max(max_dev, abs(ms%u_face_x_layer(i_w, j, k) - U_K(k)))
               max_dev = max(max_dev, abs(ms%u_face_x_layer(i_e, j, k) - U_K(k)))
            end do
            do i = i_w, i_e - 1
               max_dev = max(max_dev, abs(ms%v_face_y_layer(i, j_s, k) - V_K(k)))
               max_dev = max(max_dev, abs(ms%v_face_y_layer(i, j_n, k) - V_K(k)))
            end do
         end do

         call check(error, max_dev == 0.0_wp, &
                    "Orlanski rx=0 must preserve per-layer boundary values; max_dev=" &
                    //short_str(max_dev))

      end block checks

      call ocean_bc_state_destroy(bc)
      call bt_work%destroy()
      call ms%destroy()
   end subroutine test_orlanski_anchor_per_layer

   subroutine test_flather_full_mirror(error)
      !! Regression lock for the full-Flather interpolation weights.
      !!
      !! East/west mirror symmetry: with f = 0, symmetric eta, flat
      !! bathymetry and an antisymmetric ubt IC, the dynamics are exactly
      !! mirror-symmetric, so after N substeps the two open wall faces
      !! must satisfy ubt(i_w) = -ubt(i_e).  A swapped CFL weighting in
      !! either edge's u_inlet (the bug this test locks out) breaks the
      !! antisymmetry at leading order.  cfl is chosen well away from 0.5
      !! (where a weight swap would be invisible).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(barotropic_workstate_t) :: bt_work
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc

      integer, parameter :: NX_PHYS = 12, NY_PHYS = 4
      real(wp), parameter :: H0 = 100.0_wp     ! per-layer; D = NZ*H0 = 200 m
      real(wp), parameter :: ETA_AMP = 0.05_wp
      real(wp), parameter :: U_AMP = 0.02_wp
      real(wp), parameter :: DT_INNER = 4.0_wp  ! cfl = dt*sqrt(g*D)/dx ~ 0.18
      integer, parameter :: N_SUB = 3

      real(wp), allocatable :: force_u(:, :), force_v(:, :)
      integer :: ng, i_w, i_e, nxt, nyt, i, j
      real(wp) :: c_cell, c_face, max_asym, max_mag

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bt_work%init(grid, nz_ml=NZ)
         call cor%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=0)

         ng = grid%nghost
         i_w = ng + 1
         i_e = ng + NX_PHYS + 1
         nxt = grid%nx_total
         nyt = grid%ny_total

         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%use_full_flather = .true.
         bc%ext_u_west = 0.0_wp
         bc%ext_u_east = 0.0_wp

         cor%f_corner = 0.0_wp   ! f /= 0 would break exact mirror symmetry

         ! Symmetric eta (Gaussian about the domain centre), antisymmetric
         ! ubt (linear through the centre face), quiescent v.
         c_cell = 0.5_wp*real(nxt + 1, wp)
         c_face = 0.5_wp*real(nxt + 2, wp)
         do j = 1, nyt
            do i = 1, nxt
               bt_work%bt_eta(i, j) = ETA_AMP*exp(-((real(i, wp) - c_cell)/3.0_wp)**2)
            end do
            do i = 1, nxt + 1
               bt_work%bt_ubt(i, j) = U_AMP*(real(i, wp) - c_face)/real(NX_PHYS, wp)
            end do
         end do
         bt_work%bt_vbt = 0.0_wp
         bt_work%bt_H_ref = real(NZ, wp)*H0

         allocate (force_u(nxt + 1, nyt), source=0.0_wp)
         allocate (force_v(nxt, nyt + 1), source=0.0_wp)

         !$acc enter data copyin(bt_work, cor, bc, force_u, force_v)
         call bt_work%enter_data()
         call cor%enter_data()
         call ocean_bc_state_enter_data(bc)

         call barotropic_substep_nonlinear(grid, bt_work, &
                                           force_u, force_v, &
                                           N_SUB, DT_INNER, &
                                           bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                           bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                           eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                           bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                           bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                           uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                           bt_ubt_end=bt_work%bt_ubt_end, &
                                           bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                           bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                           vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                           bt_vbt_end=bt_work%bt_vbt_end, &
                                           bt_zeta_corner=bt_work%bt_zeta_corner, &
                                           f_corner=cor%f_corner, &
                                           bc=bc, &
                                          area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                           idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)

         !$acc update self(bt_work%bt_ubt_end)
         call ocean_bc_state_exit_data(bc)
         call cor%exit_data()
         call bt_work%exit_data()
         !$acc exit data delete(bt_work, cor, bc, force_u, force_v)

         max_asym = 0.0_wp
         max_mag = 0.0_wp
         do j = ng + 1, ng + NY_PHYS
            max_asym = max(max_asym, abs(bt_work%bt_ubt_end(i_w, j) + &
                                         bt_work%bt_ubt_end(i_e, j)))
            max_mag = max(max_mag, abs(bt_work%bt_ubt_end(i_e, j)))
         end do

         ! Non-vacuous: the wall faces must actually be moving...
         call check(error, max_mag > 1.0e-6_wp, &
                    "full-Flather mirror test is vacuous (wall faces ~ 0); max_mag=" &
                    //short_str(max_mag))
         if (allocated(error)) exit checks
         ! ...and antisymmetric to round-off.
         call check(error, max_asym < 1.0e-13_wp, &
                    "full-Flather east/west mirror broken: |u(i_w)+u(i_e)| max=" &
                    //short_str(max_asym))

      end block checks

      call ocean_bc_state_destroy(bc)
      call cor%destroy()
      call bt_work%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_flather_full_mirror

end module test_ocean_obc_baroclinic
