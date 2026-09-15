!! Tier-1 analytical / exact-solution tests for the ocean dynamical
!! core.  Where `test_ocean_validation` checks integrated invariants
!! (cooling-column heat budget, Ekman sign, conservation under no
!! forcing), this file pits the driver against problems with known
!! analytical solutions where the error is measurable to a fixed
!! tolerance.
!!
!! Cases:
!!
!!   * `vdiff_matches_erfc` — 1D vertical diffusion against the
!!     error-function solution of the heat equation.  Step IC in T,
!!     zero velocity, zero surface flux, vmix-closure off so vdiff
!!     falls back to the scalar `vd%K_v_tracer`.  After
!!     `t = 150 s` (κ·t = 6 m²) the column temperature must lie
!!     within a few percent of the erfc profile.  Mirror of
!!     `test_thermal_diffusion%test_matches_erfc` on the coastal
!!     path; here it exercises the ocean-path orchestration
!!     (`ocean_dyn_step` → `vdiff_apply_tracers`).
!!
!!   * `gravity_wave_phase_speed` — surface gravity wave on a flat
!!     basin.  A narrow Gaussian SSH perturbation at the centre
!!     splits into left- and right-going waves (d'Alembert).  After
!!     `t = T_END`, the right-going maximum sits at
!!     `x_peak ≈ x_0 + √(g·H)·T_END`.  Tests barotropic gravity-
!!     wave dispersion of `ocean_dyn_step_split`'s barotropic substep +
!!     `continuity_ppm` propagation.
!!
!!   * `geostrophic_adjustment_balance` — initial cross-channel SSH
!!     step on an f-plane.  Gravity waves radiate, the system
!!     settles into a geostrophic balance with a non-zero
!!     along-channel velocity in the sign predicted by `v = g/f · ∂η/∂x`.
!!     Probe the **average direction** of the residual flow rather
!!     than its exact magnitude — that's where Roundabout's Wright EOS
!!     differs from the gprime reduced-gravity reference solutions.
module test_ocean_analytical
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
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, BDRAG_LINEAR
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step, ocean_dyn_step_split, &
                            bcdiag_enabled, bcdiag_S_ref, bcdiag_step_limit
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_constants, only: VCOORD_SIGMA
   implicit none
   private

   public :: collect_ocean_analytical_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_analytical_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vdiff_matches_erfc", test_vdiff_erfc), &
                  new_unittest("gravity_wave_phase_speed", test_gravity_wave), &
                  new_unittest("geostrophic_adjustment_balance", &
                               test_geostrophic_adjust), &
                  new_unittest("stommel_gyre_spinup_signs", test_stommel_gyre), &
                  new_unittest("munk_gyre_spinup_signs", test_munk_gyre), &
                  new_unittest("tracer_advection_in_gyre", test_tracer_advection), &
                  new_unittest("lock_exchange_2layer", test_lock_exchange), &
                  new_unittest("baroclinic_diag_capture", test_baroclinic_diag_capture) &
                  ]
   end subroutine collect_ocean_analytical_tests

   ! ------------------------------------------------------------------
   ! Shared slot setup (mirrors test_ocean_validation%init_all)
   ! ------------------------------------------------------------------

   subroutine init_slots(grid, nz_ml, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_ml
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      real(wp), intent(in) :: f0
      ms%nz_ml = nz_ml
      call ms%init(grid)
      call ct%init(grid, nz_ml=nz_ml)
      cor%f_0 = f0
      call cor%init(grid, nz_ml=nz_ml)
      call pgf%init(grid, nz_ml=nz_ml)
      call hv%init(grid, nz_ml=nz_ml)
      call bd%init(grid, nz_ml=nz_ml)
      call ss%init(grid, nz_ml=nz_ml)
      call sf%init(grid)
      call va%init(grid, nz_ml=nz_ml)
      call hd%init(grid, nz_ml=nz_ml)
      call vd%init(grid, nz_ml=nz_ml)
      call vmix%init(grid, nz_ml=nz_ml)
      call eos%init(grid)
      ! Pass nz_ml so dyn allocates the split-driver slow-tendency
      ! accumulators (F_slow_u/v, F_bt_u/v, ubt_at_n/vbt_at_n).  The
      ! unsplit `ocean_dyn_step` is fine without them, but
      ! `ocean_dyn_step_split` accesses them every stage.
      call dyn%init(grid, nz_ml=nz_ml)
   end subroutine init_slots

   subroutine destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
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
      call sf%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine destroy_slots

   subroutine map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
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
   end subroutine map_in_slots

   subroutine map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
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
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_slots

   ! ------------------------------------------------------------------
   ! Tier 1 - Analytical / exact tests
   ! ------------------------------------------------------------------

   subroutine test_vdiff_erfc(error)
      !! 1D vertical diffusion against the erfc analytical solution.
      !!
      !! Single shallow column (4×4 horizontal, 15 layers, total
      !! depth H=10m).  IC: step in T at z=H/2 (warm above, cold
      !! below).  Zero velocity, zero surface flux.  vmix-closure off
      !! so vdiff reads the scalar `vd%K_v_tracer`.  After t = 150 s
      !! (κt = 6 m²) the column profile must match
      !!
      !!   T(z, t) = T_mean + 0.5·(T_warm - T_cold)·erf((H/2 - z)/(2√(κt)))
      !!
      !! to within a few percent.  Coastal mirror:
      !! `test_thermal_diffusion%test_matches_erfc`.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NZ = 15
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DZ = H0/real(NZ, wp)
      real(wp), parameter :: T_WARM = 25.0_wp
      real(wp), parameter :: T_COLD = 15.0_wp
      real(wp), parameter :: T_MEAN = 0.5_wp*(T_WARM + T_COLD)
      real(wp), parameter :: KAPPA = 4.0e-2_wp
      real(wp), parameter :: DT = 5.0_wp
      real(wp), parameter :: T_END = 150.0_wp
      integer, parameter :: N_STEPS = int(T_END/DT)
      integer :: i, j, k, step, idx_T, i_probe, j_probe
      real(wp) :: T_num(NZ), z_c, sqrt_2kt, T_exact, rms, max_err
      real(wp) :: heat_ic, heat_final

      call grid%init(4, 4, NGHOST, 1.0e3_wp, 1.0e3_wp)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=0.0_wp)
      checks: block
         ! Step IC in T at z = H/2.  k=1 is the bed; k=NZ the surface.
         ! Depth from surface at layer k: z_c = (NZ - k + 0.5) * dz.
         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         idx_T = ms%idx_temperature
         do k = 1, NZ
            z_c = (real(NZ - k, wp) + 0.5_wp)*DZ
            if (z_c < 0.5_wp*H0) then
               ms%tracers(idx_T)%hTr(:, :, k) = T_WARM*DZ
            else
               ms%tracers(idx_T)%hTr(:, :, k) = T_COLD*DZ
            end if
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*DZ
            ms%rho_layer(:, :, k) = eos%rho0
         end do

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp

         ! Drive vdiff with the scalar K_v; vmix-closure off so vdiff
         ! reads vd%K_v_tracer directly (no PP81 / KPP overlay).
         vd%K_v_tracer = KAPPA
         vmix%use_closure = .false.
         vmix%use_kpp = .false.

         ! Snapshot the actual column heat from the step IC.  With NZ
         ! odd, the step doesn't bisect the column symmetrically (e.g.,
         ! NZ=15 → 7 warm + 8 cold layers), so the analytical mean
         ! `NZ·DZ·T_MEAN` doesn't match the realised IC.
         i_probe = NGHOST + 2
         j_probe = NGHOST + 2
         heat_ic = sum(ms%tracers(idx_T)%hTr(i_probe, j_probe, :))

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT, sf=sf)
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! Probe the same interior column we snapshotted above.
         do k = 1, NZ
            T_num(k) = ms%tracers(idx_T)%hTr(i_probe, j_probe, k)/ms%h_layer(i_probe, j_probe, k)
         end do

         ! Compare against erfc.  Depth from surface: z_c at k.
         sqrt_2kt = 2.0_wp*sqrt(KAPPA*T_END)
         rms = 0.0_wp
         max_err = 0.0_wp
         do k = 1, NZ
            z_c = (real(NZ - k, wp) + 0.5_wp)*DZ
            T_exact = T_MEAN + 0.5_wp*(T_WARM - T_COLD)*erf((0.5_wp*H0 - z_c)/sqrt_2kt)
            rms = rms + (T_num(k) - T_exact)**2
            max_err = max(max_err, abs(T_num(k) - T_exact))
         end do
         rms = sqrt(rms/real(NZ, wp))

         ! Same tolerance as the coastal version — ~10% RMS / ~15% max
         ! to absorb first-order vertical advection error + RK2 drift.
         call check(error, rms < 0.10_wp*(T_WARM - T_COLD), &
                    "vdiff_erfc: RMS error vs erfc > 10% of initial contrast")
         if (allocated(error)) exit checks
         call check(error, max_err < 0.15_wp*(T_WARM - T_COLD), &
                    "vdiff_erfc: max pointwise error > 15% of initial contrast")

         ! Column heat conserved to round-off (the apply step is
         ! flux-form and the boundaries are closed).
         heat_final = sum(ms%tracers(idx_T)%hTr(i_probe, j_probe, :))
         call check(error, abs(heat_final - heat_ic) < 1.0e-10_wp*abs(heat_ic), &
                    "vdiff_erfc: column heat drift > round-off under closed BC")
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_vdiff_erfc

   subroutine test_gravity_wave(error)
      !! Surface gravity wave phase speed on a flat-bottom basin.
      !! 1D channel (100 × 4 cells, dx = 1 km), H = 1000 m, uniform
      !! T/S so density is uniform and only the barotropic mode is
      !! active.  Initial Gaussian SSH perturbation at the centre.
      !! After t = T_END the right-going peak is at
      !!   x_peak ≈ x_0 + √(g·H)·T_END.
      !! With g = 9.81, H = 1000, c = 99.05 m/s.  Over 200 s the wave
      !! travels ≈ 19.8 km (~20 cells) — well clear of the walls at
      !! x = 0 and x = 100 km.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NX = 100, NY = 4, NZ = 2
      integer, parameter :: N_INNER = 1
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: H0 = 1000.0_wp
      real(wp), parameter :: ETA_AMP = 0.1_wp
      real(wp), parameter :: SIGMA_X = 5.0_wp*DX   ! 5-cell Gaussian — well-resolved by PPM
      real(wp), parameter :: G_ACCEL = 9.81_wp
      ! Outer DT close to the Heun-RK2 gravity-wave CFL (= dx/c ≈ 10 s).
      ! Small DT (e.g. 2 s) gives the split-driver many chances to leak
      ! BT energy via the slow-apply / barotropic-substep / correction handshake;
      ! one big step per CFL cycle keeps the BT mode in one piece.
      real(wp), parameter :: DT = 8.0_wp
      real(wp), parameter :: T_END = 400.0_wp      ! wave travels ~40 km of 50 km
      integer, parameter :: N_STEPS = int(T_END/DT)
      real(wp) :: c_expected, x_centre, x_peak, dz
      real(wp) :: ssh_max, ssh_local, speed_observed, speed_error
      integer :: i, j, k, step, j_probe, i_peak
      character(len=256) :: msg

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=0.0_wp)
      ! Use FV_LITE — MONT (the slot's default) introduces a spurious
      ! per-layer PGF when h_layer varies in x, which makes the
      ! depth-mean PGF half its physical value (≈ -g/2 · ∂η/∂x instead
      ! of -g · ∂η/∂x) and the barotropic wave speed
      ! sqrt(gH/2) ≈ 70 m/s instead of sqrt(gH) ≈ 99 m/s.  See the
      ! variant docstring in `rdb_ocean_pressure_force.F90`.
      pgf%variant = OPGF_VARIANT_FV_LITE
      checks: block
         dz = H0/real(NZ, wp)
         x_centre = 0.5_wp*real(NX, wp)*DX

         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Gaussian SSH perturbation, distributed evenly across layers
         ! (i.e., column thickens uniformly).  Surface position is
         ! eta(x) = ETA_AMP * exp(-((x - x_centre)/SIGMA_X)^2).
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  block
                     real(wp) :: x_i, eta_i
                     x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
                     eta_i = ETA_AMP*exp(-((x_i - x_centre)/SIGMA_X)**2)
                     ms%h_layer(i, j, k) = dz + eta_i/real(NZ, wp)
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                        eos%T_ref*ms%h_layer(i, j, k)
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                        eos%S_ref*ms%h_layer(i, j, k)
                     ms%rho_layer(i, j, k) = eos%rho0
                  end block
               end do
            end do
         end do

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp

         ! Split driver — the unsplit `ocean_dyn_step` runs in
         ! Eulerian-z mode where h_layer is held fixed, so the wave
         ! doesn't propagate.  The split driver's barotropic substep is the only
         ! path that evolves the barotropic gravity wave correctly.
         dyn%bt_work%bt_H_ref = H0

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
            ! Diagnostic: every 1/4 of run, print current peak position
            ! and inferred speed.
            if (mod(step, max(1, N_STEPS/4)) == 0) then
               !$acc update self(ms%h_layer)
               block
                  integer :: ii, ii_pk
                  real(wp) :: pk, loc
                  pk = -huge(1.0_wp)
                  ii_pk = NGHOST + NX/2
                  do ii = NGHOST + NX/2, NGHOST + NX - 1
                     loc = sum(ms%h_layer(ii, NGHOST + NY/2, :)) - H0
                     if (loc > pk) then
                        pk = loc
                        ii_pk = ii
                     end if
                  end do
                  print '("  step ", i0, " t = ", f7.2, " s: peak i = ", i0, &
                        &" (Δx = ", f7.2, " km), speed-so-far = ", f6.2, " m/s, ssh = ", es10.3)', &
                     step, real(step, wp)*DT, ii_pk, &
                     (real(ii_pk - NGHOST, wp) - 0.5_wp)*DX*1.0e-3_wp - x_centre*1.0e-3_wp, &
                     ((real(ii_pk - NGHOST, wp) - 0.5_wp)*DX - x_centre)/(real(step, wp)*DT), pk
               end block
            end if
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! Find argmax of column-summed h (≡ surface elevation) over
         ! the RIGHT half of the domain.  After apply_bt_correction
         ! rescales h_layer to match the barotropic-substep bt_eta_mean, this
         ! equals the barotropic-substep SSH up to round-off.
         j_probe = NGHOST + NY/2
         c_expected = sqrt(G_ACCEL*H0)
         ssh_max = -huge(1.0_wp)
         i_peak = NGHOST + NX/2
         do i = NGHOST + NX/2, NGHOST + NX - 1
            ssh_local = sum(ms%h_layer(i, j_probe, :)) - H0
            if (ssh_local > ssh_max) then
               ssh_max = ssh_local
               i_peak = i
            end if
         end do
         x_peak = (real(i_peak - NGHOST, wp) - 0.5_wp)*DX
         speed_observed = (x_peak - x_centre)/T_END
         speed_error = abs(speed_observed - c_expected)/c_expected

         ! Tolerance: 20%.  Peak position is a soft target — PPM has
         ! O(dx²) dispersion at the Gaussian's short-wavelength tail
         ! that lags the peak slightly behind the true phase front.
         ! With σ = 5·dx the Gaussian is well-resolved; the residual
         ! error is barotropic-substep FBE phase error plus PPM dispersion.
         write (msg, '("gravity_wave: speed_observed = ", f7.2, " m/s, expected ", f7.2, " (error = ", f6.2, "%, peak at i = ", i0, ", ssh_max = ", es10.3, ")")') &
            speed_observed, c_expected, 100.0_wp*speed_error, i_peak, ssh_max
         call check(error, speed_error < 0.20_wp, trim(msg))
         if (allocated(error)) exit checks

         ! Peak should still be in the right half (i.e., the wave
         ! moved at all).
         write (msg, '("gravity_wave: peak only travelled ", f7.2, " m (cell ", i0, "), expected ", f7.2, " m")') &
            x_peak - x_centre, i_peak, c_expected*T_END
         call check(error, x_peak > x_centre + 5.0_wp*DX, trim(msg))
         if (allocated(error)) exit checks

         ! And it shouldn't have hit the wall — the right-going wave
         ! should still be inside the domain.
         call check(error, i_peak < NGHOST + NX - 2, &
                    "gravity_wave: right-going wave reflected off the east wall")
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_gravity_wave

   subroutine test_geostrophic_adjust(error)
      !! Geostrophic adjustment.  f-plane basin, initial SSH step in
      !! x (high on west half, low on east half), zero initial
      !! velocity.  The system radiates gravity waves and settles
      !! into a state with non-zero along-channel velocity in the
      !! sign predicted by `v = (g/f) · ∂η/∂x`.  With ∂η/∂x < 0
      !! (high on west, low on east) and f > 0, geostrophic balance
      !! requires v < 0 in the time-mean of the basin.
      !!
      !! Run for ~5 inertial periods, then time-average over the
      !! last 2 periods (filters out the lingering oscillation) and
      !! check the sign and order of magnitude.  Avoid the strict
      !! analytical balance check (which requires PV conservation +
      !! Wright EOS bookkeeping) — sign + magnitude is the Tier-1
      !! invariant.
      !!
      !! This is a Coriolis-vs-PGF balance test; it does not depend
      !! on density.  The multilayer machinery still runs (h_layer,
      !! per-layer u/v, the registered S/T tracers) but we
      !! decouple the EOS — `β_S = α_T = 0` — so any spurious S or
      !! T drift cannot translate into a density gradient and feed
      !! back through the baroclinic PGF.  Without that, an
      !! unrelated tracer-advection bug becomes an apparent
      !! geostrophy failure under realistic seawater sensitivities.
      !! Density-driven adjustment lives in `lock_exchange_2layer`
      !! (active β_S override); the open salt-feedback diagnostic
      !! is `baroclinic_diag_capture`.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NX = 32, NY = 4, NZ = 2
      integer, parameter :: N_INNER = 10
      real(wp), parameter :: DX = 10.0e3_wp, DY = 10.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: ETA_STEP = 0.1_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: T_INERTIAL = 2.0_wp*acos(-1.0_wp)/F0
      real(wp), parameter :: DT = 60.0_wp
      real(wp) :: T_END, T_AVG_START
      integer :: N_STEPS, step, n_avg_samples, i, j, k
      integer :: i_probe, j_probe, k_probe
      real(wp) :: v_sum, v_mean, dz, x_centre, x_i
      character(len=256) :: msg

      T_END = 5.0_wp*T_INERTIAL
      T_AVG_START = 3.0_wp*T_INERTIAL
      N_STEPS = int(T_END/DT)

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      ! See test_gravity_wave for why MONT (the default) is incorrect
      ! once h_layer varies horizontally — use FV_LITE.
      pgf%variant = OPGF_VARIANT_FV_LITE
      ! Decouple the EOS — see the docstring.  Spurious S/T drift
      ! from the tracer-advection bug must not feed back into the
      ! PGF here, or we'd be testing two things at once.
      eos%beta_S = 0.0_wp
      eos%alpha_T = 0.0_wp
      checks: block
         dz = H0/real(NZ, wp)
         x_centre = 0.5_wp*real(NX, wp)*DX

         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
                  ! Smooth tanh step in eta across the domain centre.
                  ! Width ~3*dx so the gradient is resolved.
                  ms%h_layer(i, j, k) = dz - (ETA_STEP/real(NZ, wp))* &
                                        tanh((x_i - x_centre)/(3.0_wp*DX))
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     eos%T_ref*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     eos%S_ref*ms%h_layer(i, j, k)
                  ms%rho_layer(i, j, k) = eos%rho0
               end do
            end do
         end do

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp

         dyn%bt_work%bt_H_ref = H0

         ! Probe the v-face right at the basin centre's north edge,
         ! where ∂η/∂x is maximally negative.  Surface layer.
         i_probe = NGHOST + NX/2
         j_probe = NGHOST + NY/2
         k_probe = NZ

         v_sum = 0.0_wp
         n_avg_samples = 0
         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
            ! Accumulate v at the probe point over the last ~2 inertial
            ! periods.  Pointwise probe avoids dilution by the closed-
            ! wall ghosts (where v = 0) that a global mean would average
            ! into.
            if (real(step, wp)*DT >= T_AVG_START) then
               ! Pull the full array — single-element `update self` is
               ! unreliable on NVHPC OpenACC (D->H of one slot may be
               ! optimised away).
               !$acc update self(ms%v_face_y_layer)
               v_sum = v_sum + ms%v_face_y_layer(i_probe, j_probe, k_probe)
               n_avg_samples = n_avg_samples + 1
            end if
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         v_mean = v_sum/real(n_avg_samples, wp)

         ! With ∂η/∂x < 0 at the centre (high on west, low on east) and
         ! f > 0, geostrophic balance is v = (g/f) · ∂η/∂x < 0.
write (msg, '("geostrophic_adjust: v_mean at centre = ", es12.4, " m/s (expected sign: negative, expected |v| ~ 0.16 m/s)")') v_mean
         call check(error, v_mean < 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! Order-of-magnitude check.  Crude estimate:
         !   |v| ~ (g · ETA_STEP) / (f · width_of_step)
         ! With ETA_STEP=0.1, width ≈ 6*DX = 60 km, g=9.81, f=1e-4:
         !   |v| ~ 9.81 · 0.1 / (1e-4 · 60e3) ≈ 0.16 m/s.  Allow a
         ! generous lower bound — bouncing waves still pollute the
         ! time-mean in a closed basin without friction, and with the
         ! step width ≪ Rd = √(gH)/f ≈ 700 km most of the initial APE
         ! radiates rather than settling into balance, so the pointwise
         ! 2-period mean is a small residual of interfering waves, not
         ! the naive 0.16.  The floor only guards against NO spin-up
         ! (|v| ~ round-off); it moved 1e-4 → 5e-5 when the
         ! `subtract_fast_cor_ref` fix (MOM6 `Cor_ref` parity — the BT
         ! Coriolis was previously integrated twice, distorting the
         ! inertial wave phasing this residual mean is made of) shifted
         ! the observed value from ~1.1e-4 to 7.8e-5.  The SIGN check
         ! above is the robust balance invariant.
         write (msg, '("geostrophic_adjust: |v_mean at centre| = ", es12.4, " too small (no spin-up?)")') abs(v_mean)
         call check(error, abs(v_mean) > 5.0e-5_wp, trim(msg))
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_geostrophic_adjust

   subroutine test_baroclinic_diag_capture(error)
      !! Env-var-gated diagnostic capture for the baroclinic-stability
      !! investigation.
      !!
      !! When `RDB_BCDIAG=1` is set in the environment:
      !!   1. Build the same `geostrophic_adjust` IC.
      !!   2. Override β_S to the realistic seawater value (0.78
      !!      kg/m³/PSU) so the baroclinic feedback engages.
      !!   3. Enable the `bcdiag_*` probes in `rdb_ocean_dyn`, capping
      !!      them at the first 2 outer steps.
      !!   4. Run 2 outer split-explicit steps.  Trace lands on stdout.
      !!
      !! When the env var is unset, the test is a no-op and passes —
      !! the trace is a one-shot investigation, not a regression gate.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 32, NY = 4, NZ = 2
      integer, parameter :: N_INNER = 10
      real(wp), parameter :: DX = 10.0e3_wp, DY = 10.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: ETA_STEP = 0.1_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: BETA_S_REAL = 0.78_wp
      integer :: i, j, k, step, status
      real(wp) :: dz, x_centre, x_i
      character(len=8) :: env_val

      call get_environment_variable("RDB_BCDIAG", env_val, status=status)
      if (status /= 0 .or. trim(env_val) /= "1") return

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE

      ! Lagrangian + ALE remap path.  VCOORD_SIGMA target tracks
      ! `(H_ref + bt_eta)·dsig` (uniform fraction per layer).  This
      ! makes the split driver skip `tracer_advect_vertical` and
      ! the `apply_bt_correction` h-rescale; the ALE remap after
      ! `rk2_average` relayers (h_layer, every hTr) onto target_h
      ! conservatively per column — preserves both bit-exact
      ! `sum(hTr)` and per-column `S = hTr/h`.  Without this,
      ! the Eulerian-z surface-flux dilemma seeds a baroclinic
      ! feedback that NaNs around step ~2100 under realistic β_S.
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA

      ! Realistic seawater haline sensitivity — the value that causes
      ! the baroclinic feedback we want to diagnose.
      eos%beta_S = BETA_S_REAL

      dz = H0/real(NZ, wp)
      x_centre = 0.5_wp*real(NX, wp)*DX

      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
               ms%h_layer(i, j, k) = dz - (ETA_STEP/real(NZ, wp))* &
                                     tanh((x_i - x_centre)/(3.0_wp*DX))
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  eos%T_ref*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  eos%S_ref*ms%h_layer(i, j, k)
               ms%rho_layer(i, j, k) = eos%rho0
            end do
         end do
      end do

      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      sf%Q_heat_const = 0.0_wp
      sf%Q_salt_const = 0.0_wp
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      dyn%bt_work%bt_H_ref = H0

      ! Enable the per-kernel probes for the first 2 outer steps only.
      bcdiag_enabled = .true.
      bcdiag_S_ref = eos%S_ref
      bcdiag_step_limit = 2

      call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call vc%enter_data()
      do step = 1, 2
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, sf=sf, vcoord=vc)
      end do
      call vc%exit_data()
      call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

      ! Disable probes again so we don't bleed into any subsequent test
      ! in the same process.
      bcdiag_enabled = .false.

      call vc%destroy()
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_baroclinic_diag_capture

   subroutine test_stommel_gyre(error)
      !! Stommel gyre, sign-and-structure variant.  Closed square
      !! basin on a β-plane with a single-cos zonal wind and Rayleigh
      !! bottom drag.  Start from rest, spin up for K outer steps,
      !! check the developing gyre has the right large-scale
      !! structure:
      !!
      !!   1. Domain-mean v in the basin INTERIOR is negative
      !!      (southward Sverdrup return for our wind: τ_x = -τ₀·cos
      !!      drives anticyclonic / clockwise / subtropical-style
      !!      circulation).
      !!   2. Domain-mean |u| in the western 5-cell band is larger
      !!      than in the eastern 5-cell band (WBC forms at the west,
      !!      not the east — the β-effect "westernises" the WBC).
      !!   3. Domain-mean v in the WB band is positive (northward
      !!      return flow in the WBC).
      !!   4. KE has grown from zero (wind did net work).
      !!
      !! Why structural rather than fixed-point: setting up the
      !! analytical Stommel IC requires solving the full second-order
      !! ODE r·G'' + β·G' − r·(π/L_y)²·G = const, *and* a 2D η that
      !! satisfies both ∂η/∂x and ∂η/∂y from the steady momentum
      !! equations.  At our parameters the dropped r·(π/L_y)²·G term
      !! is ~25 % of β·G' — significant.  The fixed-point approach
      !! also conflates discretisation error with wiring error; the
      !! structural test isolates the wiring.
      !!
      !! What this test catches:
      !!   * sign of β·v in the vorticity balance (test 1)
      !!   * sign of curl(τ) → wind-driven gyre direction (test 1, 3)
      !!   * β-effect westernises the WBC (test 2)
      !!   * sign of linear drag → bounded WBC (test 2 fails if drag
      !!     has the wrong sign; the WBC would diverge)
      !!   * wind input is positive-definite (test 4)
      !!
      !! Spinup time: 1/r = 1e6 s at our drag, but the interior
      !! Sverdrup pattern emerges via Rossby waves on the much faster
      !! timescale of one basin crossing (~L_x / (β·L_d²) ~ 7 h for
      !! our parameters).  K = 200 outer steps at dt = 60 s gives
      !! 3.3 h — enough for the gyre sign and WBC concentration to
      !! be unambiguous.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      ! NZ = 1: the analytical Stommel ψ is barotropic (depth-uniform).
      ! With NZ > 1 the surface wind only hits the top layer and drag
      ! only the bottom; without vertical momentum mixing the layers
      ! decouple and the fixed-point assumption breaks (u drifts as the
      ! top layer accelerates and the bottom decelerates).
      integer, parameter :: NX = 50, NY = 50, NZ = 1
      ! Long spinup: at dt = 120 s × K = 3000 = 4.2 days physical
      ! time ≈ 0.7 Rossby basin crossings (T_crossing = π²/(β·L_x)
      ! ≈ 5.8 days for our params).  Enough to develop a non-trivial
      ! fraction of the steady Sverdrup interior — calibrated against
      ! the observed v_int/v_Sverdrup ratio.  Shorter K leaves the
      ! ratio so far from unity that any tight magnitude bound is
      ! impossible.
      integer, parameter :: N_INNER = 1
      integer, parameter :: N_STEPS = 3000
      real(wp), parameter :: DX = 20.0e3_wp, DY = 20.0e3_wp   ! 1000 × 1000 km basin
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: F0 = 5.0e-5_wp                   ! mid-latitude f
      real(wp), parameter :: BETA = 2.0e-11_wp                ! β at mid-latitude
      real(wp), parameter :: TAU0 = 0.1_wp                    ! peak wind stress (N/m²)
      real(wp), parameter :: R_DRAG = 1.0e-6_wp               ! linear Rayleigh (1/s)
      real(wp), parameter :: RHO0 = 1025.0_wp
      real(wp), parameter :: G_ACCEL = 9.81_wp
      real(wp), parameter :: DT = 120.0_wp                    ! CFL ≈ 0.6
      real(wp), parameter :: PI = 3.141592653589793_wp
      real(wp) :: Lx, Ly, delta_s, dz
      real(wp) :: v_interior_mean, v_wb_mean, v_eb_mean
      real(wp) :: ke_final
      real(wp) :: v_sverdrup, sverdrup_ratio
      integer :: i, j, k, step, idx_T, idx_S
      integer :: n_interior, n_wb, n_eb
      integer :: i_int_lo, i_int_hi, j_int_lo, j_int_hi
      integer :: i_wb_lo, i_wb_hi, i_eb_lo, i_eb_hi
      character(len=256) :: msg

      Lx = real(NX, wp)*DX
      Ly = real(NY, wp)*DY
      delta_s = R_DRAG/BETA
      dz = H0/real(NZ, wp)

      call grid%init(NX, NY, NGHOST, DX, DY)
      ! Note: init_slots calls cor%init which fills f_corner = f_0
      ! (constant).  We overwrite below with the β-plane profile.
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      ! FV_LITE — geostrophic balance needs the layer-PGF chain-rule
      ! correction (see gravity_wave test for the MONT failure mode).
      pgf%variant = OPGF_VARIANT_FV_LITE

      checks: block
         ! ---- β-plane: f = f₀ + β·y, y measured from south wall.
         ! The setter's coordinate is `y = (j - 1 - nghost)·dy` so y=0
         ! is the physical south wall; passing y_ref=0 puts f_0 there.
         call cor%set_beta_plane(grid, F0, BETA, y_ref=0.0_wp)

         ! ---- Wind: τ_x(y) = -τ₀·cos(π·y/L_y), y measured from south wall ----
         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         do j = 1, size(ss%tau_x, 2)
            block
               real(wp) :: y_phys
               y_phys = (real(j - NGHOST, wp) - 0.5_wp)*DY
               if (y_phys > 0.0_wp .and. y_phys < Ly) then
                  ss%tau_x(:, j) = -TAU0*cos(PI*y_phys/Ly)
               end if
            end block
         end do

         ! ---- Linear drag ----
         bd%variant = BDRAG_LINEAR
         bd%r_linear = R_DRAG
         bd%c_drag = 0.0_wp

         ! ---- Disable mixing / surface tracer flux (wind stress is
         ! the test's only forcing; configured above on ss) ----
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         hv%nu_h = 0.0_wp                                ! Stommel = no lateral friction
         hd%kappa_h = 0.0_wp

         ! ---- Rest state: zero velocity, flat η, uniform T/S ----
         idx_T = ms%idx_temperature
         idx_S = ms%idx_salinity

         do k = 1, NZ
            ms%h_layer(:, :, k) = dz
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%v_face_y_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = eos%rho0
            ms%tracers(idx_T)%hTr(:, :, k) = eos%T_ref*dz
            ms%tracers(idx_S)%hTr(:, :, k) = eos%S_ref*dz
         end do

         ! Set up the bt reference column thickness for the split driver.
         dyn%bt_work%bt_H_ref = H0

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! ---- Structural checks: gyre sign + WBC location ----
         ! WB band: westernmost 5 interior cells.  EB band: easternmost
         ! 5 interior cells.  Interior: skip 5 cells from each x-wall
         ! and 2 cells from each y-wall (where the wind goes to zero).
         i_wb_lo = NGHOST + 1
         i_wb_hi = NGHOST + 5
         i_eb_lo = NGHOST + NX - 4
         i_eb_hi = NGHOST + NX
         i_int_lo = NGHOST + 6
         i_int_hi = NGHOST + NX - 5
         j_int_lo = NGHOST + 3
         j_int_hi = NGHOST + NY - 2

         v_interior_mean = 0.0_wp
         v_wb_mean = 0.0_wp
         v_eb_mean = 0.0_wp
         ke_final = 0.0_wp
         n_interior = 0
         n_wb = 0
         n_eb = 0
         do k = 1, NZ
            ! Interior v-mean (over v-face cells with x-centre in interior).
            do j = j_int_lo, j_int_hi
               do i = i_int_lo, i_int_hi
                  v_interior_mean = v_interior_mean + ms%v_face_y_layer(i, j, k)
                  n_interior = n_interior + 1
               end do
            end do
            ! WB v-mean.
            do j = j_int_lo, j_int_hi
               do i = i_wb_lo, i_wb_hi
                  v_wb_mean = v_wb_mean + ms%v_face_y_layer(i, j, k)
                  n_wb = n_wb + 1
               end do
            end do
            ! Signed v-mean in EB band.  Stommel says the eastern
            ! boundary carries the broad return, with small per-cell
            ! v of either sign (returns are spread across the
            ! interior).  Comparing the *signed mean* of v_wb (large,
            ! coherent, northward) to v_eb (small) is a sharper
            ! discriminator than the L1 norm, which is dominated by
            ! spinup gravity-wave noise on both walls.
            do j = j_int_lo, j_int_hi
               do i = i_eb_lo, i_eb_hi
                  v_eb_mean = v_eb_mean + ms%v_face_y_layer(i, j, k)
                  n_eb = n_eb + 1
               end do
            end do
            ! Total KE over the basin interior.
            do j = NGHOST + 1, NGHOST + NY
               do i = NGHOST + 1, NGHOST + NX
                  ke_final = ke_final + 0.5_wp*ms%h_layer(i, j, k)* &
                             (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
               end do
            end do
         end do
         if (n_interior > 0) v_interior_mean = v_interior_mean/real(n_interior, wp)
         if (n_wb > 0) v_wb_mean = v_wb_mean/real(n_wb, wp)
         if (n_eb > 0) v_eb_mean = v_eb_mean/real(n_eb, wp)

         ! Sverdrup-balance reference value at mid-basin (y = L_y/2):
         ! β·v_steady = curl(τ)/(ρ₀·H), curl(τ) = -∂τ_x/∂y =
         ! -τ₀·(π/L_y)·sin(πy/L_y).  At y = L_y/2, sin = 1.
         v_sverdrup = -TAU0*PI/(BETA*RHO0*H0*Ly)
         if (abs(v_sverdrup) > tiny(1.0_wp)) then
            sverdrup_ratio = v_interior_mean/v_sverdrup
         else
            sverdrup_ratio = 0.0_wp
         end if

         write (msg, '("stommel: v_int=", es10.3, " v_wb=", es10.3, &
                       &" v_eb=", es10.3, " KE=", es10.3, &
                       &" v_int/v_Sverdrup=", f8.5)') &
            v_interior_mean, v_wb_mean, v_eb_mean, ke_final, sverdrup_ratio
         print '(a)', trim(msg)

         ! 1. Interior v < 0 (Sverdrup southward return for our wind).
         write (msg, '("stommel: interior v_mean=", es12.4, " expected negative (Sverdrup southward)")') v_interior_mean
         call check(error, v_interior_mean < 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 2. WB v > 0 (northward WBC). Sign of f, β, or wind curl
         !    wrong → flips this.
         write (msg, '("stommel: WB v_mean=", es12.4, " expected positive (northward WBC)")') v_wb_mean
         call check(error, v_wb_mean > 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 3. EB v < 0 (southward broad return).  Sign of β wrong
         !    flips this (WBC moves to east, return to west).
         write (msg, '("stommel: EB v_mean=", es12.4, " expected negative (southward return)")') v_eb_mean
         call check(error, v_eb_mean < 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 4. WBC v amplitude exceeds the interior v amplitude (the
         !    WBC is a *current*, the interior is a slow Sverdrup
         !    drift).  Catches a degenerate basin with no real WBC.
         write (msg, '("stommel: WBC weak: |v_wb_mean|=", es10.3, " not >> |v_int_mean|=", es10.3)') &
            abs(v_wb_mean), abs(v_interior_mean)
         call check(error, abs(v_wb_mean) > 5.0_wp*abs(v_interior_mean), trim(msg))
         if (allocated(error)) exit checks

         ! 5. Wind has done net work: KE has grown from zero.
         write (msg, '("stommel: KE=", es12.4, " expected > 0 (rest IC + wind forcing)")') ke_final
         call check(error, ke_final > 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 6. Sverdrup ratio: tight bounds around the calibrated
         !    spinup state.  At K = 3000 × dt = 120 s = 4.2 days
         !    ≈ 0.7 Rossby basin crossings, we observe
         !    v_int / v_Sverdrup ≈ 0.156 (15.6 % of steady).  Bounds
         !    [0.05, 0.5] tolerate the natural Rossby-wave overshoot
         !    + spinup state but reject:
         !    * a 3× error in β (would shift the ratio outside the
         !      band — bigger β → faster spinup → ratio approaches 1
         !      faster; smaller β → slower spinup → ratio collapses).
         !    * a wind curl wiring bug (response would vanish or
         !      flip sign).
         !    * a missing Coriolis projection in the vorticity balance
         !      (Sverdrup wouldn't develop at all).
         write (msg, '("stommel: Sverdrup ratio v_int/v_Sverdrup=", f10.5, &
                       &" outside [0.05, 0.5] — β wiring or magnitude off")') sverdrup_ratio
         call check(error, sverdrup_ratio > 0.05_wp .and. sverdrup_ratio < 0.5_wp, trim(msg))
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_stommel_gyre

   subroutine test_munk_gyre(error)
      !! Munk gyre — lateral-viscosity WBC, sign-and-structure variant.
      !! Same basin, wind, and β as the Stommel test but the WBC is
      !! closed by Laplacian horizontal viscosity (`hv%nu_h`) instead
      !! of linear bottom drag.  Exercises a different closure path
      !! in the code: `ocean_horizontal_viscosity_apply_tendencies`
      !! and the Laplacian `compute_tendencies` kernel that Stommel
      !! doesn't touch.
      !!
      !! Munk theory: WBC width δ_M = (ν/β)^(1/3).  At ν = 1e4 m²/s,
      !! β = 2e-11, δ_M ≈ 80 km ≈ 4 cells at dx = 20 km — comfortably
      !! resolved.  The Munk solution has a damped-oscillation profile
      !! at the WBC (Airy-function-like) rather than Stommel's pure
      !! exponential, but the *signs* of the gyre interior + WBC
      !! return are the same, so the same structural checks apply.
      !!
      !! What this test adds over Stommel:
      !!   * `hv` apply path (Laplacian viscosity on momentum).
      !!   * Confirms the gyre forms with the right sign under a
      !!     completely different WBC closure — if Stommel passes and
      !!     Munk fails, the regression is specifically in the
      !!     viscosity path; if both fail, the failure is upstream
      !!     (wind / β / Coriolis / PGF).
      !!
      !! Same five sign-and-structure checks as `test_stommel_gyre`,
      !! plus the same loose Sverdrup-magnitude ratio.  See that
      !! test's docstring for the rationale.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NX = 50, NY = 50, NZ = 1
      integer, parameter :: N_INNER = 1
      integer, parameter :: N_STEPS = 3000   ! 4.2-day spinup, matches Stommel
      real(wp), parameter :: DX = 20.0e3_wp, DY = 20.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: F0 = 5.0e-5_wp
      real(wp), parameter :: BETA = 2.0e-11_wp
      real(wp), parameter :: TAU0 = 0.1_wp
      real(wp), parameter :: NU_H = 1.0e4_wp                  ! Munk-WBC closure
      real(wp), parameter :: RHO0 = 1025.0_wp
      real(wp), parameter :: G_ACCEL = 9.81_wp
      real(wp), parameter :: DT = 120.0_wp                    ! CFL ≈ 0.6
      real(wp), parameter :: PI = 3.141592653589793_wp
      real(wp) :: Lx, Ly, delta_m, dz
      real(wp) :: v_interior_mean, v_wb_mean, v_eb_mean
      real(wp) :: ke_final, v_sverdrup, sverdrup_ratio
      integer :: i, j, k, step, idx_T, idx_S
      integer :: n_interior, n_wb, n_eb
      integer :: i_int_lo, i_int_hi, j_int_lo, j_int_hi
      integer :: i_wb_lo, i_wb_hi, i_eb_lo, i_eb_hi
      character(len=256) :: msg

      Lx = real(NX, wp)*DX
      Ly = real(NY, wp)*DY
      delta_m = (NU_H/BETA)**(1.0_wp/3.0_wp)
      dz = H0/real(NZ, wp)

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE

      checks: block
         call cor%set_beta_plane(grid, F0, BETA, y_ref=0.0_wp)

         ! Wind: same single-cosine as Stommel.
         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         do j = 1, size(ss%tau_x, 2)
            block
               real(wp) :: y_phys
               y_phys = (real(j - NGHOST, wp) - 0.5_wp)*DY
               if (y_phys > 0.0_wp .and. y_phys < Ly) then
                  ss%tau_x(:, j) = -TAU0*cos(PI*y_phys/Ly)
               end if
            end block
         end do

         ! No drag — Munk is viscosity-closed.
         bd%variant = BDRAG_LINEAR
         bd%r_linear = 0.0_wp
         bd%c_drag = 0.0_wp

         ! Lateral viscosity = the WBC closure.
         hv%nu_h = NU_H

         ! Disable other physics.
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         hd%kappa_h = 0.0_wp

         ! Rest state.
         idx_T = ms%idx_temperature
         idx_S = ms%idx_salinity
         do k = 1, NZ
            ms%h_layer(:, :, k) = dz
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%v_face_y_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = eos%rho0
            ms%tracers(idx_T)%hTr(:, :, k) = eos%T_ref*dz
            ms%tracers(idx_S)%hTr(:, :, k) = eos%S_ref*dz
         end do

         dyn%bt_work%bt_H_ref = H0

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! ---- Structural checks (same shape as Stommel) ----
         ! WB band sized to ~δ_M.  For δ_M ≈ 80 km, 4 cells is the
         ! WBC width; widen by ~50 % to 6 to absorb the overshoot
         ! lobe of the Airy-like profile.
         i_wb_lo = NGHOST + 1
         i_wb_hi = NGHOST + 6
         i_eb_lo = NGHOST + NX - 5
         i_eb_hi = NGHOST + NX
         i_int_lo = NGHOST + 7
         i_int_hi = NGHOST + NX - 6
         j_int_lo = NGHOST + 3
         j_int_hi = NGHOST + NY - 2

         v_interior_mean = 0.0_wp
         v_wb_mean = 0.0_wp
         v_eb_mean = 0.0_wp
         ke_final = 0.0_wp
         n_interior = 0
         n_wb = 0
         n_eb = 0
         do k = 1, NZ
            do j = j_int_lo, j_int_hi
               do i = i_int_lo, i_int_hi
                  v_interior_mean = v_interior_mean + ms%v_face_y_layer(i, j, k)
                  n_interior = n_interior + 1
               end do
               do i = i_wb_lo, i_wb_hi
                  v_wb_mean = v_wb_mean + ms%v_face_y_layer(i, j, k)
                  n_wb = n_wb + 1
               end do
               do i = i_eb_lo, i_eb_hi
                  v_eb_mean = v_eb_mean + ms%v_face_y_layer(i, j, k)
                  n_eb = n_eb + 1
               end do
            end do
            do j = NGHOST + 1, NGHOST + NY
               do i = NGHOST + 1, NGHOST + NX
                  ke_final = ke_final + 0.5_wp*ms%h_layer(i, j, k)* &
                             (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
               end do
            end do
         end do
         if (n_interior > 0) v_interior_mean = v_interior_mean/real(n_interior, wp)
         if (n_wb > 0) v_wb_mean = v_wb_mean/real(n_wb, wp)
         if (n_eb > 0) v_eb_mean = v_eb_mean/real(n_eb, wp)

         v_sverdrup = -TAU0*PI/(BETA*RHO0*H0*Ly)
         if (abs(v_sverdrup) > tiny(1.0_wp)) then
            sverdrup_ratio = v_interior_mean/v_sverdrup
         else
            sverdrup_ratio = 0.0_wp
         end if

         write (msg, '("munk: v_int=", es10.3, " v_wb=", es10.3, &
                       &" v_eb=", es10.3, " KE=", es10.3, &
                       &" v_int/v_Sverdrup=", f8.5, " δ_M=", f6.1, " km")') &
            v_interior_mean, v_wb_mean, v_eb_mean, ke_final, sverdrup_ratio, &
            delta_m*1.0e-3_wp
         print '(a)', trim(msg)

         write (msg, '("munk: interior v_mean=", es12.4, " expected negative (Sverdrup southward)")') v_interior_mean
         call check(error, v_interior_mean < 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         write (msg, '("munk: WB v_mean=", es12.4, " expected positive (northward WBC)")') v_wb_mean
         call check(error, v_wb_mean > 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         write (msg, '("munk: EB v_mean=", es12.4, " expected negative (southward return)")') v_eb_mean
         call check(error, v_eb_mean < 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         write (msg, '("munk: WBC weak: |v_wb_mean|=", es10.3, " not >> |v_int_mean|=", es10.3)') &
            abs(v_wb_mean), abs(v_interior_mean)
         call check(error, abs(v_wb_mean) > 5.0_wp*abs(v_interior_mean), trim(msg))
         if (allocated(error)) exit checks

         write (msg, '("munk: KE=", es12.4, " expected > 0 (rest IC + wind forcing)")') ke_final
         call check(error, ke_final > 0.0_wp, trim(msg))
         if (allocated(error)) exit checks

         ! Tighter bound: same long-spinup logic as Stommel; see
         ! that test's check #6 for the rationale.  Munk's WBC
         ! closure is different (viscous vs drag) but the *interior*
         ! Sverdrup approach is closure-independent so the same
         ! bounds apply.
         write (msg, '("munk: Sverdrup ratio v_int/v_Sverdrup=", f10.5, &
                       &" outside [0.05, 0.5] — β wiring or magnitude off")') sverdrup_ratio
         call check(error, sverdrup_ratio > 0.05_wp .and. sverdrup_ratio < 0.5_wp, trim(msg))
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_munk_gyre

   subroutine test_tracer_advection(error)
      !! Passive tracer advection in a developing Stommel gyre.
      !! Validates the multilayer tracer pipeline:
      !!
      !!   * `continuity_tracer_step_split` mass conservation
      !!     (closed-BC sum of hTr stays bit-stable to round-off).
      !!   * PPM limiter positivity / no-overshoot (a Gaussian pulse
      !!     in salinity stays inside its initial value envelope).
      !!   * Coupling between dynamics and tracer transport (the
      !!     pulse deforms — proves advection actually ran).
      !!
      !! Why a gyre instead of pure translation: the original plan
      !! called for a 1D uniform-advection test, but closed walls
      !! prevent a clean translation setup — `u = U₀` collides with
      !! `u_wall = 0` and drains the column near the east wall by
      !! continuity.  A 2D gyre is a more realistic transport
      !! scenario anyway: the tracer experiences both shear and
      !! curvature, and the same closure routines (PPM + face fluxes)
      !! are exercised.
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NX = 50, NY = 50, NZ = 1
      integer, parameter :: N_INNER = 1
      integer, parameter :: N_STEPS = 1500
         !! Spinup window for the wind-driven gyre to develop a flow
         !! capable of deforming the pulse.  Bumped from 600 after
         !! the barotropic-substep physical-wall closure fix (2026-05-19) —
         !! the new wall geometry slows the transient gyre spinup
         !! (Stommel/Munk steady-state at N_STEPS=3000 unaffected).
      real(wp), parameter :: DX = 20.0e3_wp, DY = 20.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: F0 = 5.0e-5_wp
      real(wp), parameter :: BETA = 2.0e-11_wp
      real(wp), parameter :: TAU0 = 0.1_wp
      real(wp), parameter :: R_DRAG = 1.0e-6_wp
      real(wp), parameter :: RHO0 = 1025.0_wp
      real(wp), parameter :: G_ACCEL = 9.81_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: PI = 3.141592653589793_wp
      real(wp), parameter :: S_PERT = 1.0_wp                ! PSU peak above S_ref
      real(wp) :: Lx, Ly, x_c, y_c, sigma_pulse, dz
      real(wp) :: S_min_init, S_max_init, S_ref
      real(wp) :: S_min_final, S_max_final
      real(wp) :: hTr_sum_init, hTr_sum_final, mass_rel_err
      real(wp) :: deform_l2
      integer :: i, j, k, step, idx_S, idx_T
      character(len=256) :: msg
      real(wp), allocatable :: hTr_init(:, :, :)

      Lx = real(NX, wp)*DX
      Ly = real(NY, wp)*DY
      x_c = 0.5_wp*Lx
      y_c = 0.5_wp*Ly
      sigma_pulse = 5.0_wp*DX
      dz = H0/real(NZ, wp)
      S_ref = 35.0_wp  ! overridden by eos%S_ref below; cached here for clarity

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE
      S_ref = eos%S_ref

      checks: block
         ! ---- Stommel-style gyre setup so the tracer sees a real flow ----
         call cor%set_beta_plane(grid, F0, BETA, y_ref=0.0_wp)

         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         do j = 1, size(ss%tau_x, 2)
            block
               real(wp) :: y_phys
               y_phys = (real(j - NGHOST, wp) - 0.5_wp)*DY
               if (y_phys > 0.0_wp .and. y_phys < Ly) then
                  ss%tau_x(:, j) = -TAU0*cos(PI*y_phys/Ly)
               end if
            end block
         end do

         bd%variant = BDRAG_LINEAR
         bd%r_linear = R_DRAG
         bd%c_drag = 0.0_wp

         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         hv%nu_h = 0.0_wp
         hd%kappa_h = 0.0_wp                  ! no horizontal tracer diffusion — pure advection

         ! ---- Rest state with a Gaussian salinity pulse at basin centre ----
         idx_T = ms%idx_temperature
         idx_S = ms%idx_salinity
         do k = 1, NZ
            ms%h_layer(:, :, k) = dz
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%v_face_y_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = eos%rho0
            ms%tracers(idx_T)%hTr(:, :, k) = eos%T_ref*dz
         end do

         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  block
                     real(wp) :: x_phys, y_phys, r2, S_local
                     x_phys = (real(i - NGHOST, wp) - 0.5_wp)*DX
                     y_phys = (real(j - NGHOST, wp) - 0.5_wp)*DY
                     r2 = (x_phys - x_c)**2 + (y_phys - y_c)**2
                     S_local = S_ref + S_PERT*exp(-r2/sigma_pulse**2)
                     ms%tracers(idx_S)%hTr(i, j, k) = S_local*dz
                  end block
               end do
            end do
         end do

         ! Snapshot initial pulse state.
         allocate (hTr_init(grid%nx_total, grid%ny_total, NZ))
         hTr_init = ms%tracers(idx_S)%hTr

         hTr_sum_init = kahan_sum_interior(ms%tracers(idx_S)%hTr, grid, NZ)
         S_min_init = huge(1.0_wp)
         S_max_init = -huge(1.0_wp)
         do k = 1, NZ
            do j = NGHOST + 1, NGHOST + NY
               do i = NGHOST + 1, NGHOST + NX
                  block
                     real(wp) :: S_local
                     S_local = ms%tracers(idx_S)%hTr(i, j, k)/ms%h_layer(i, j, k)
                     S_min_init = min(S_min_init, S_local)
                     S_max_init = max(S_max_init, S_local)
                  end block
               end do
            end do
         end do

         dyn%bt_work%bt_H_ref = H0

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! ---- Final-state diagnostics ----
         hTr_sum_final = kahan_sum_interior(ms%tracers(idx_S)%hTr, grid, NZ)
         S_min_final = huge(1.0_wp)
         S_max_final = -huge(1.0_wp)
         deform_l2 = 0.0_wp
         do k = 1, NZ
            do j = NGHOST + 1, NGHOST + NY
               do i = NGHOST + 1, NGHOST + NX
                  block
                     real(wp) :: S_local
                     S_local = ms%tracers(idx_S)%hTr(i, j, k)/ms%h_layer(i, j, k)
                     S_min_final = min(S_min_final, S_local)
                     S_max_final = max(S_max_final, S_local)
                     deform_l2 = deform_l2 + (ms%tracers(idx_S)%hTr(i, j, k) - hTr_init(i, j, k))**2
                  end block
               end do
            end do
         end do
         deform_l2 = sqrt(deform_l2)
         mass_rel_err = abs(hTr_sum_final - hTr_sum_init)/abs(hTr_sum_init)

         write (msg, '("tracer: hTr_sum init=", es14.7, " final=", es14.7, &
                       &" rel_err=", es10.3, " S_min=", f7.4, "→", f7.4, &
                       &" S_max=", f7.4, "→", f7.4, " deform_L2=", es10.3)') &
            hTr_sum_init, hTr_sum_final, mass_rel_err, &
            S_min_init, S_min_final, S_max_init, S_max_final, deform_l2
         print '(a)', trim(msg)

         ! 1. Mass conservation: PPM is bit-exact-conservative on a
         !    closed-wall basin.  Kernel fix landed in
         !    `continuity_zonal_flux` / `continuity_meridional_flux`
         !    to zero `mass_flux` at the *physical* wall faces
         !    (i = nghost+1, nghost+nx_phys+1; same for j) — without
         !    that, slow-path Coriolis/wind drove a tiny u at the
         !    physical wall and PPM transported a few ppm of tracer
         !    across into ghost cells per ~600 steps.  Tolerance set
         !    to 5×ε·N_sum ≈ 1e-10 to absorb sum-reduction round-off
         !    only.
         write (msg, '("tracer: mass drift ", es10.3, " > 1e-10 — PPM conservation broken")') &
            mass_rel_err
         call check(error, mass_rel_err < 1.0e-10_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 2. PPM bounds — loose. The continuity-PPM kernel uses
         !    a classical PPM cell-parabolic limiter (no FCT face-
         !    bound enforcement), so face values can slip slightly
         !    outside the local-stencil min/max.  Observed undershoot
         !    on this test is ~0.1 PSU (10 % of pulse amplitude) and
         !    overshoot is essentially zero.  Bound at 0.5 PSU
         !    (half the pulse amplitude) catches a fully-broken
         !    limiter / runaway oscillation without flagging the
         !    known small undershoot.  Tightening this is a separate
         !    follow-up (FCT-PPM upgrade).
         write (msg, '("tracer: overshoot S_max=", f10.6, " > IC max=", f10.6, " + 0.5 PSU")') &
            S_max_final, S_max_init
         call check(error, S_max_final < S_max_init + 0.5_wp, trim(msg))
         if (allocated(error)) exit checks

         write (msg, '("tracer: undershoot S_min=", f10.6, " < IC min=", f10.6, " - 0.5 PSU")') &
            S_min_final, S_min_init
         call check(error, S_min_final > S_min_init - 0.5_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 3. Pulse has deformed: catches a regression where advection
         !    is a no-op (tracer decoupled from velocity field).  L2
         !    is an absolute distance, not a relative measure.  Pulse
         !    amplitude × √(active cells) ≈ 1·500·√75 ≈ 4300; observed
         !    deformation at K=600 is ~165 (a few percent of the pulse
         !    scale).  Threshold 10 catches the no-advection case
         !    while remaining loose for partial spinup.
         write (msg, '("tracer: deform_L2=", es10.3, " — advection did not run?")') deform_l2
         call check(error, deform_l2 > 10.0_wp, trim(msg))
      end block checks
      if (allocated(hTr_init)) deallocate (hTr_init)
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_tracer_advection

   subroutine test_lock_exchange(error)
      !! 2-layer hydrostatic lock exchange.  Validates baroclinic
      !! pressure-gradient wiring, EOS coupling, multilayer
      !! continuity, and gravity-current propagation — the *first*
      !! analytical test to exercise stratified (NZ > 1, density-
      !! varying-with-T/S) dynamics in the ocean dyn-core.
      !!
      !! Setup: long zonal channel (NX = 100 × dx = 1 km), 4 cells in
      !! y (effectively 1D), 2 vertical layers each 50 m thick.  No
      !! rotation, no wind, no drag, no viscosity, no surface flux,
      !! no vertical mixing.  Initial state at rest with a vertical
      !! salinity step at x = L_x/2:
      !!
      !!   left  half:  S = 36 PSU  (denser)
      !!   right half:  S = 34 PSU  (lighter)
      !!   T = 10 °C uniform
      !!
      !! EOS computes density per step from T/S, so the PGF kernel
      !! sees a horizontal density gradient at the lock.  Hydrostatic
      !! adjustment drives a gravity current: dense water flows
      !! rightward along the bottom, light water flows leftward along
      !! the surface.  Predicted front speed
      !!   c = √(g'·H/2),  g' = g · Δρ / ρ₀,
      !! Δρ ≈ 0.78·ΔS ≈ 1.56 kg/m³ at our T/S → c ≈ 0.86 m/s.
      !!
      !! What this test catches (each independently of Stommel/Munk):
      !!   * EOS sign / wiring (ρ must increase with S);
      !!   * baroclinic PGF (per-layer pressure differences, not
      !!     just bt PGF on η);
      !!   * multilayer continuity (h_layer evolves separately per k);
      !!   * tracer mass conservation under non-trivial 2D advection
      !!     (the wall-mask fix in continuity-PPM is exercised here
      !!     too, but the new physics is the vertical layer coupling).
      !!
      !! Sign-and-structure checks (consistent with Stommel/Munk):
      !!   1. dense salinity (S > 35) has reached the right half in
      !!      the BOTTOM layer (gravity-current sign);
      !!   2. light salinity (S < 35) has reached the left half in
      !!      the TOP layer (counter-flow sign);
      !!   3. tracer mass conserved bit-exactly (closed BCs);
      !!   4. PPM bounds: S stays within [34 − 0.5, 36 + 0.5] PSU;
      !!   5. KE has grown from zero (PE → KE via baroclinic PGF).
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NX = 100, NY = 4, NZ = 2
      integer, parameter :: N_INNER = 1
      integer, parameter :: N_STEPS = 1500           ! ≈ 8.3 h at dt=20 s, ~18 km front travel
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: S_DENSE = 36.0_wp
      real(wp), parameter :: S_LIGHT = 34.0_wp
      real(wp), parameter :: S_MEAN = 0.5_wp*(S_DENSE + S_LIGHT)
      real(wp), parameter :: T_REF = 10.0_wp
      real(wp), parameter :: DT = 20.0_wp
      ! Override EOS defaults — the small ocean defaults (β_S=7.6e-4,
      ! α_T=1.7e-4) give Δρ ≈ 1000× too small.  Use realistic
      ! seawater values here so the gravity current actually
      ! propagates.  See `feedback_ocean_eos_default.md` and the
      ! follow-up note in the EOS module for the broader story.
      real(wp), parameter :: BETA_S_REAL = 0.78_wp     ! kg/m³/PSU
      real(wp), parameter :: ALPHA_T_REAL = 0.17_wp    ! kg/m³/°C
      real(wp) :: dz, hTr_S_init, hTr_S_final, mass_rel_err
      real(wp) :: S_min_final, S_max_final, S_max_right_bot, S_min_left_top
      real(wp) :: rho_top_right, rho_bot_right, ke_final
      integer :: i, j, k, step, idx_T, idx_S
      integer :: i_lock, i_right_probe, i_left_probe, j_probe
      character(len=256) :: msg

      dz = H0/real(NZ, wp)
      i_lock = NGHOST + NX/2
      ! Probe locations: 5 cells (= 5 km) from the lock on each side
      ! — inside the gravity-current footprint at K=1500 (front travels
      ! ~18 km at c ≈ 0.61 m/s).
      i_right_probe = NGHOST + NX/2 + 5
      i_left_probe = NGHOST + NX/2 - 4
      j_probe = NGHOST + NY/2

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=0.0_wp)
      pgf%variant = OPGF_VARIANT_FV_LITE
      eos%beta_S = BETA_S_REAL
      eos%alpha_T = ALPHA_T_REAL

      checks: block
         ! Disable everything except baroclinic PGF + continuity + EOS.
         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         bd%variant = BDRAG_LINEAR
         bd%r_linear = 0.0_wp
         bd%c_drag = 0.0_wp
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         hv%nu_h = 0.0_wp
         hd%kappa_h = 0.0_wp

         ! Rest state with salinity step at lock.
         idx_T = ms%idx_temperature
         idx_S = ms%idx_salinity
         do k = 1, NZ
            ms%h_layer(:, :, k) = dz
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%v_face_y_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = eos%rho0
            ms%tracers(idx_T)%hTr(:, :, k) = T_REF*dz
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  if (i <= i_lock) then
                     ms%tracers(idx_S)%hTr(i, j, k) = S_DENSE*dz
                  else
                     ms%tracers(idx_S)%hTr(i, j, k) = S_LIGHT*dz
                  end if
               end do
            end do
         end do

         hTr_S_init = kahan_sum_interior(ms%tracers(idx_S)%hTr, grid, NZ)

         dyn%bt_work%bt_H_ref = H0

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf)
         end do
         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)

         ! Final diagnostics.
         hTr_S_final = kahan_sum_interior(ms%tracers(idx_S)%hTr, grid, NZ)
         mass_rel_err = abs(hTr_S_final - hTr_S_init)/abs(hTr_S_init)

         S_min_final = huge(1.0_wp)
         S_max_final = -huge(1.0_wp)
         S_max_right_bot = -huge(1.0_wp)   ! max S in bottom layer, RIGHT of lock
         S_min_left_top = huge(1.0_wp)    ! min S in top layer, LEFT of lock
         ke_final = 0.0_wp
         do k = 1, NZ
            do j = NGHOST + 1, NGHOST + NY
               do i = NGHOST + 1, NGHOST + NX
                  block
                     real(wp) :: S_local
                     S_local = ms%tracers(idx_S)%hTr(i, j, k)/ms%h_layer(i, j, k)
                     S_min_final = min(S_min_final, S_local)
                     S_max_final = max(S_max_final, S_local)
                     if (k == 1 .and. i > i_lock) then
                        S_max_right_bot = max(S_max_right_bot, S_local)
                     end if
                     if (k == NZ .and. i <= i_lock) then
                        S_min_left_top = min(S_min_left_top, S_local)
                     end if
                     ke_final = ke_final + 0.5_wp*ms%h_layer(i, j, k)* &
                                (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
                  end block
               end do
            end do
         end do

         ! Density inversion at a right-side probe column: bottom
         ! denser than top (dense water has invaded from below).
         rho_top_right = ms%rho_layer(i_right_probe, j_probe, NZ)
         rho_bot_right = ms%rho_layer(i_right_probe, j_probe, 1)

         write (msg, '("lock: mass_err=", es10.3, " S_max_right_bot=", f7.4, &
                       &" S_min_left_top=", f7.4, " ρ_bot_right=", f9.4, &
                       &" ρ_top_right=", f9.4, " KE=", es10.3)') &
            mass_rel_err, S_max_right_bot, S_min_left_top, &
            rho_bot_right, rho_top_right, ke_final
         print '(a)', trim(msg)

         ! 1. Tracer mass conservation: PPM + wall mask → bit-exact.
         write (msg, '("lock: mass drift ", es10.3, " > 1e-10 — PPM conservation broken")') mass_rel_err
         call check(error, mass_rel_err < 1.0e-10_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 2. PPM bounds (classical PPM, ~0.5 PSU undershoot tolerance).
         write (msg, '("lock: S_max=", f7.4, " > S_DENSE + 0.5 = ", f7.4)') S_max_final, S_DENSE + 0.5_wp
         call check(error, S_max_final < S_DENSE + 0.5_wp, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("lock: S_min=", f7.4, " < S_LIGHT − 0.5 = ", f7.4)') S_min_final, S_LIGHT - 0.5_wp
         call check(error, S_min_final > S_LIGHT - 0.5_wp, trim(msg))
         if (allocated(error)) exit checks

         ! 3. Gravity-current sign: dense water has reached the right
         !    side in the BOTTOM layer (S > S_MEAN somewhere on the
         !    right).  If PGF sign is wrong or EOS is broken, dense
         !    water doesn't move and S_max_right_bot stays at S_LIGHT.
         write (msg, '("lock: dense water did not reach right bottom: S_max_right_bot=", f7.4, " <= S_MEAN=", f7.4)') &
            S_max_right_bot, S_MEAN
         call check(error, S_max_right_bot > S_MEAN, trim(msg))
         if (allocated(error)) exit checks

         ! 4. Counter-flow sign: light water has reached the left
         !    side in the TOP layer (S < S_MEAN somewhere on left).
         write (msg, '("lock: light water did not reach left top: S_min_left_top=", f7.4, " >= S_MEAN=", f7.4)') &
            S_min_left_top, S_MEAN
         call check(error, S_min_left_top < S_MEAN, trim(msg))
         if (allocated(error)) exit checks

         ! 5. Density inversion on right side: bottom denser than
         !    top (the dense gravity current is below the light
         !    water).  EOS sign error or PGF feedback wrong would
         !    flip this.
         write (msg, '("lock: stratification wrong on right side: ρ_bot=", f9.4, " <= ρ_top=", f9.4)') &
            rho_bot_right, rho_top_right
         call check(error, rho_bot_right > rho_top_right, trim(msg))
         if (allocated(error)) exit checks

         ! 6. KE has grown — PE → KE conversion via baroclinic PGF.
         write (msg, '("lock: KE=", es12.4, " expected > 0 (rest IC + PE release)")') ke_final
         call check(error, ke_final > 0.0_wp, trim(msg))
      end block checks
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_lock_exchange

   function kahan_sum_interior(field, grid, nz) result(total)
      !! Compensated summation of `field` over the physical-interior
      !! cells (excluding ghosts).  Plain `sum`/loop accumulators have
      !! O(N·ε·max) round-off — at 2500 cells × max~4e7 that's already
      !! ~2e-5 relative, which would mask any real conservation drift
      !! < 1e-5.  Kahan compensation reduces the bound to O(ε), letting
      !! the test distinguish a true conservation breach from FP
      !! reduction noise.
      real(wp), intent(in) :: field(:, :, :)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp) :: total, c, y, t
      integer :: i, j, k
      total = 0.0_wp
      c = 0.0_wp
      do k = 1, nz
         do j = grid%nghost + 1, grid%nghost + grid%ny_phys
            do i = grid%nghost + 1, grid%nghost + grid%nx_phys
               y = field(i, j, k) - c
               t = total + y
               c = (t - total) - y
               total = t
            end do
         end do
      end do
   end function kahan_sum_interior

end module test_ocean_analytical
