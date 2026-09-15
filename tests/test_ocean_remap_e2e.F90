!! End-to-end validation for the ALE remap pipeline wired through
!! `ocean_dyn_step_split`.  Drives one or more outer steps
!! of the production split driver with a non-trivial `vcoord%coord_type`
!! and asserts that the resting-state guarantees hold to round-off:
!!
!!   - Lake-at-rest under VCOORD_SIGMA: stratified S/T + zero velocity
!!     + flat bath + η=0 → state unchanged across one full RK2 step
!!     that includes the remap.  Identity-remap case: target_h equals
!!     h_layer to the bit, so PLM column kernel collapses to copy.
!!   - Lake-at-rest under VCOORD_ZSTAR_FULL with a stretched per-column
!!     z_ref: exercises the most subtle target_h path against the
!!     trivial dynamics tendency.  Initial h_layer is set to the
!!     z_ref intervals so the first remap is also identity.
!!
!! These tests prove the e2e wire-up is correct: the optional `vcoord`
!! reaches the orchestrator, GPU mapping (vcoord%enter_data → device
!! arrays) survives a full driver step, and the remap is genuinely
!! a no-op when the target grid matches the Lagrangian grid.
!!
!! A third intended test — multi-step conservation under VCOORD_SIGMA
!! with non-trivial flow — is deferred to Layer 3c follow-up: the
!! orchestrator alone is per-column-conservative (verified in
!! `test_ocean_remap`), but the interaction with the split driver's
!! barotropic-substep bt_eta accumulator produces a small non-round-off drift
!! (~6e-6 per step) that needs separate diagnosis.  Documented in
!! docs/ROADMAP_OCEAN.md Phase 5g.
!!
!! Runs on the same harness as `test_ocean_dyn_split` (the pre-existing
!! split-driver suite) — same map_in/map_out shape, plus the new
!! vcoord argument.
module test_ocean_remap_e2e
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, VCOORD_SIGMA, VCOORD_ZSTAR, VCOORD_ZSTAR_FULL
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
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
   use rdb_ocean_vcoord, only: ocean_vcoord_t, STRETCH_UNIFORM
   implicit none
   private

   public :: collect_ocean_remap_e2e_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_remap_e2e_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("e2e_lake_at_rest_sigma", test_lake_at_rest_sigma), &
                  new_unittest("e2e_lake_at_rest_zstar_full", test_lake_at_rest_zstar_full), &
                  new_unittest("e2e_mass_tracer_conservation_sigma", test_conservation_sigma), &
                  new_unittest("e2e_remap_cadence_zstar", test_remap_cadence_zstar), &
                  new_unittest("e2e_remap_cadence_ratio1", test_remap_cadence_ratio1) &
                  ]
   end subroutine collect_ocean_remap_e2e_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c, h_per_layer)
      type(hgrid_t), intent(inout) :: grid
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
      type(ocean_vcoord_t), intent(inout) :: vc
      real(wp), intent(in) :: f_c, h_per_layer
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = f_c
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
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
      dyn%bt_work%bt_H_ref = real(NZ, wp)*h_per_layer
      ms%h_layer = h_per_layer
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine init_all

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
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
      type(ocean_vcoord_t), intent(inout) :: vc
      call vc%destroy()
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

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
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
      type(ocean_vcoord_t), intent(inout) :: vc
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
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
      call vc%enter_data()
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
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
      type(ocean_vcoord_t), intent(inout) :: vc
      call destroy_cartesian_metrics(metrics)
      call vc%exit_data()
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
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
   end subroutine map_out

   subroutine stamp_stratification(ms, eos, S_top, T_top, dSdk, dTdk)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: S_top, T_top, dSdk, dTdk
      real(wp) :: S_k(NZ), T_k(NZ)
      integer :: k
      ! k=NZ is the surface, k=1 the bed.  Top-anchored.
      do k = 1, NZ
         S_k(k) = S_top + dSdk*real(NZ - k, wp)
         T_k(k) = T_top + dTdk*real(NZ - k, wp)
      end do
      if (.false.) S_k(1) = eos%S_ref   ! reference EOS access; silence unused
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*ms%h_layer(:, :, k)
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*ms%h_layer(:, :, k)
      end do
   end subroutine stamp_stratification

   subroutine test_lake_at_rest_sigma(error)
      !! Stratified S/T, zero velocity, flat bath, eta=0, VCOORD_SIGMA.
      !! After 1 RK2 step including the remap, every prognostic should
      !! stay at its IC to round-off.  The remap target should equal
      !! the initial h_layer exactly (uniform → uniform under SIGMA at
      !! eta=0), so the remap is an identity transformation.
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
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_INNER = 5
      real(wp), allocatable :: hS_ic(:, :, :), hT_ic(:, :, :)
      real(wp) :: max_dh, max_du, max_dv, max_dS, max_dT
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c=1.0_wp, h_per_layer=H0)
         vc%coord_type = VCOORD_SIGMA
         call stamp_stratification(ms, eos, S_top=eos%S_ref - 1.0_wp, T_top=eos%T_ref + 2.0_wp, &
                                   dSdk=0.5_wp, dTdk=-1.0_wp)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         max_dh = maxval(abs(ms%h_layer - H0))
         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_dh < 1.0e-10_wp, "lake-at-rest SIGMA: h drifted")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-10_wp, "lake-at-rest SIGMA: u drifted")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp, "lake-at-rest SIGMA: v drifted")
         if (allocated(error)) exit checks
         call check(error, max_dS < 1.0e-10_wp, "lake-at-rest SIGMA: salt drifted")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-10_wp, "lake-at-rest SIGMA: heat drifted")

      end block checks
      deallocate (hS_ic, hT_ic)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
   end subroutine test_lake_at_rest_sigma

   subroutine test_lake_at_rest_zstar_full(error)
      !! Lake-at-rest under VCOORD_ZSTAR_FULL.  The wrinkle is that the
      !! initial h_layer must be set to the reference z_ref intervals,
      !! not uniform — otherwise the first remap shuffles mass and the
      !! "no drift" guarantee won't hold.
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
      real(wp), parameter :: H_BED = 200.0_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_INNER = 5
      real(wp), allocatable :: h_bed_2d(:, :), hS_ic(:, :, :), hT_ic(:, :, :)
      real(wp) :: max_dh, max_du, max_dv, max_dS, max_dT
      integer :: i, j, k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c=1.0_wp, h_per_layer=H_BED/real(NZ, wp))
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 20.0_wp
         vc%zstar_n_surf = 2
         vc%zstar_stretching = STRETCH_UNIFORM
         allocate (h_bed_2d(grid%nx_total, grid%ny_total), source=H_BED)
         call vc%build_zref_full(h_bed_2d)

         ! Initialise h_layer to the z_ref intervals so the first remap
         ! is identity.  ROMS-ordered: k=1 bed, k=NZ surface.
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  ms%h_layer(i, j, k) = vc%z_ref(i, j, NZ - k + 1) - vc%z_ref(i, j, NZ - k)
               end do
            end do
         end do
         dyn%bt_work%bt_H_ref = H_BED
         call stamp_stratification(ms, eos, S_top=eos%S_ref - 1.0_wp, T_top=eos%T_ref + 2.0_wp, &
                                   dSdk=0.5_wp, dTdk=-1.0_wp)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         max_dh = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  max_dh = max(max_dh, abs(ms%h_layer(i, j, k) - &
                                           (vc%z_ref(i, j, NZ - k + 1) - vc%z_ref(i, j, NZ - k))))
               end do
            end do
         end do
         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_dh < 1.0e-9_wp, "lake-at-rest ZSTAR_FULL: h drifted")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-9_wp, "lake-at-rest ZSTAR_FULL: u drifted")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-9_wp, "lake-at-rest ZSTAR_FULL: v drifted")
         if (allocated(error)) exit checks
         call check(error, max_dS < 1.0e-9_wp, "lake-at-rest ZSTAR_FULL: salt drifted")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-9_wp, "lake-at-rest ZSTAR_FULL: heat drifted")

      end block checks
      deallocate (h_bed_2d, hS_ic, hT_ic)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
   end subroutine test_lake_at_rest_zstar_full

   subroutine test_conservation_sigma(error)
      !! Total volume + tracer mass conservation through 10 split-driver
      !! steps under VCOORD_SIGMA, with a small solenoidal u/v
      !! perturbation.  Same IC pattern as `test_conservation` in
      !! test_ocean_dyn_split (which clears 1e-10 without vcoord);
      !! this is the vcoord-on counterpart.
      !!
      !! With PPM as the default remap method the spurious vertical-
      !! mixing artefact from a 1st-order limiter is eliminated, so the
      !! tracer drift drops to round-off scale.
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
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 10
      integer, parameter :: N_INNER = 5
      integer :: i, j, k, nx, ny, step
      real(wp) :: total_h0, total_S0, total_T0
      real(wp) :: total_h, total_S, total_T, drift_h, drift_S, drift_T
      checks: block

         call make_grid(grid, 24, 16, 1.0_wp, 1.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c=F_C, h_per_layer=H_BASE)
         vc%coord_type = VCOORD_SIGMA
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + 0.01_wp*sin( &
                                        2.0_wp*PI*real(i, wp)/real(nx, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = ms%h_layer(i, j, k)*( &
                                                             eos%S_ref + 0.5_wp*cos(2.0_wp*PI*real(j, wp)/real(ny, wp)))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = ms%h_layer(i, j, k)*( &
                                                                eos%T_ref + 0.3_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                               sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp))
               end do
            end do
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H_BASE

         total_h0 = sum(ms%h_layer)*grid%dx*grid%dy
         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T0 = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         total_h = sum(ms%h_layer)*grid%dx*grid%dy
         total_S = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         drift_h = abs(total_h - total_h0)/abs(total_h0)
         drift_S = abs(total_S - total_S0)/abs(total_S0)
         drift_T = abs(total_T - total_T0)/abs(total_T0)

         call check(error, drift_h < 1.0e-10_wp, "SIGMA: total h drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, drift_S < 1.0e-10_wp, "SIGMA: total S mass drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, drift_T < 1.0e-10_wp, "SIGMA: total T mass drift > 1e-10")
         if (allocated(error)) exit checks
         call check(error, minval(ms%h_layer) > 0.0_wp, "SIGMA: h_layer went non-positive")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
   end subroutine test_conservation_sigma

   subroutine test_remap_cadence_zstar(error)
      !! DT_THERM cadence gate: the ALE remap fires on every EVEN outer step
      !! (`dt_therm_ratio = 2`, thermo steps at `outer_step_count = 0, 2, 4, …`)
      !! and is SKIPPED on ODD steps.
      !!
      !! Setup: VCOORD_ZSTAR with a small tilt in h_layer and a non-zero
      !! barotropic velocity so continuity provably moves mass each step.
      !! After the thermo step (step 1, outer_step_count=0) the remap fires
      !! and h_layer is set to target_h by the remap kernel — verified by
      !! pulling both arrays from device.  After the non-thermo step (step 2,
      !! outer_step_count=1) the remap is skipped, continuity has evolved
      !! h_layer, but device target_h is stale (unchanged by the remap) —
      !! so max|h_layer - target_h| must exceed the round-off floor.
      !! Conservation: total h and tracer mass are conserved to 1e-10 across
      !! all steps (Lagrangian continuity is still conservative).
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
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_INNER = 5
      integer :: i, j, k, nx, ny
      real(wp) :: total_h0, total_S0, total_T0
      real(wp) :: total_h, total_S, total_T
      real(wp) :: max_diff_after_thermo, max_diff_after_nonthermo
      checks: block

         call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c=1.0_wp, h_per_layer=H_BASE)
         vc%coord_type = VCOORD_ZSTAR
         nx = grid%nx_total
         ny = grid%ny_total

         ! Small tilt in h_layer so continuity moves mass each step and
         ! h_layer diverges from target_h on a non-thermo step.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + 0.02_wp*sin( &
                                        2.0_wp*PI*real(i, wp)/real(nx, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = ms%h_layer(i, j, k)* &
                                                             (eos%S_ref + 0.3_wp*sin( &
                                                              PI*real(j, wp)/real(ny, wp)))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = ms%h_layer(i, j, k)* &
                                                                (eos%T_ref + 0.2_wp*cos( &
                                                                 PI*real(i, wp)/real(nx, wp)))
               end do
            end do
            ! Non-zero face velocities so continuity changes h_layer.
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))
               end do
            end do
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H_BASE
         dyn%dt_therm_ratio = 2

         total_h0 = sum(ms%h_layer)*grid%dx*grid%dy
         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T0 = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         ! Step 1: outer_step_count=0 => thermo step => remap fires.
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         ! Pull both h_layer and target_h to host.
         !$acc update self(ms%h_layer, vc%target_h)
         ! After the remap, h_layer == target_h (the remap kernel writes
         ! ms%h_layer = vcoord%target_h — step 6 of ocean_apply_ale_remap_step).
         max_diff_after_thermo = maxval(abs(ms%h_layer - vc%target_h))
         call check(error, max_diff_after_thermo < 1.0e-10_wp, &
                    "cadence ZSTAR: h_layer /= target_h after thermo step")
         if (allocated(error)) exit checks

         ! Step 2: outer_step_count=1 => non-thermo step => remap skipped.
         ! Continuity still advances h_layer (Lagrangian), so h_layer drifts
         ! away from the stale target_h that was set during step 1.
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         !$acc update self(ms%h_layer, vc%target_h)
         max_diff_after_nonthermo = maxval(abs(ms%h_layer - vc%target_h))
         call check(error, max_diff_after_nonthermo > 1.0e-8_wp, &
                    "cadence ZSTAR: h_layer == target_h after non-thermo step (remap ran when it should not)")
         if (allocated(error)) exit checks

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         ! Conservation: both steps are conservative regardless of whether
         ! remap fires (Lagrangian continuity conserves mass; remap conserves
         ! by its own discrete accounting).
         total_h = sum(ms%h_layer)*grid%dx*grid%dy
         total_S = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         call check(error, abs(total_h - total_h0)/abs(total_h0) < 1.0e-10_wp, &
                    "cadence ZSTAR: total h not conserved")
         if (allocated(error)) exit checks
         call check(error, abs(total_S - total_S0)/abs(total_S0) < 1.0e-10_wp, &
                    "cadence ZSTAR: total S not conserved")
         if (allocated(error)) exit checks
         call check(error, abs(total_T - total_T0)/abs(total_T0) < 1.0e-10_wp, &
                    "cadence ZSTAR: total T not conserved")

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
   end subroutine test_remap_cadence_zstar

   subroutine test_remap_cadence_ratio1(error)
      !! Control: dt_therm_ratio=1 (default) remaps every step.
      !! After EACH outer step h_layer == target_h to round-off.
      !! This verifies the bit-identical default path still remaps every step.
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
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_INNER = 5
      integer :: i, j, k, nx, ny, step
      real(wp) :: max_diff
      checks: block

         call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc, &
                       f_c=1.0_wp, h_per_layer=H_BASE)
         vc%coord_type = VCOORD_ZSTAR
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + 0.02_wp*sin( &
                                        2.0_wp*PI*real(i, wp)/real(nx, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = ms%h_layer(i, j, k)*eos%S_ref
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = ms%h_layer(i, j, k)*eos%T_ref
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))
               end do
            end do
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H_BASE
         ! dt_therm_ratio = 1 is the default; set explicitly for clarity.
         dyn%dt_therm_ratio = 1

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

         do step = 1, 4
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
            !$acc update self(ms%h_layer, vc%target_h)
            max_diff = maxval(abs(ms%h_layer - vc%target_h))
            call check(error, max_diff < 1.0e-10_wp, &
                       "cadence ratio=1: h_layer /= target_h (remap skipped on every-step path)")
            if (allocated(error)) exit
         end do

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, vc)

      end block checks
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn, vc)
   end subroutine test_remap_cadence_ratio1

end module test_ocean_remap_e2e
