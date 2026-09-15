!! Long-run baroclinic-stability regression.  Carved out from
!! `test_ocean_analytical` so we can iterate on it without
!! re-running the seven analytical cases each rebuild.
!!
!! Scope: drive the split-explicit RK2 driver under realistic
!! seawater haline sensitivity (`β_S = 0.78 kg/m³/PSU`) with
!! the Lagrangian + ALE remap path (`VCOORD_SIGMA`).  The
!! salt-feedback seed that previously NaN'd the `geostrophic_
!! adjust` IC around step ~2615 has been closed by Phase 1 (MOM6
!! transport constraint on slow continuity) + Phase 2
!! (Lagrangian + ALE remap), so `max|hTr_S/h_layer - S_ref|`
!! should stay at FP round-off indefinitely.  This test asserts
!! that, plus the basic invariants (no NaN, h_layer stays
!! positive, global mass / tracer conservation to round-off).
!!
!! Long-run NaN around step ~1400 — separate (non-tracer)
!! corner-case still under investigation.  The test runs to a
!! step count BEFORE that NaN so it's a stable gate; the
!! step count can be pushed up as the corner case gets fixed.
module test_ocean_baroclinic_longrun
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, VCOORD_SIGMA
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, &
                            bcdiag_enabled, bcdiag_S_ref, bcdiag_step_limit
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_budgets, only: ocean_budgets_t, &
                                BUDGET_MASS, BUDGET_SALT_TOTAL, BUDGET_HEAT_TOTAL, BUDGET_KE
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private

   ! Split into one collector per subtest so each registers as its own
   ! ctest binary (rdb_test_ocean_baroclinic_{geostrophic,stommel,budget}) —
   ! lets the heavy 20000-step geostrophic run be focused/skipped in isolation.
   public :: collect_ocean_baroclinic_geostrophic_tests
   public :: collect_ocean_baroclinic_stommel_tests
   public :: collect_ocean_baroclinic_budget_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_baroclinic_geostrophic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("baroclinic_geostrophic_realistic_betaS_lagrangian", &
                               test_baroclinic_longrun) &
                  ]
   end subroutine collect_ocean_baroclinic_geostrophic_tests

   subroutine collect_ocean_baroclinic_stommel_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("baroclinic_stommel_realistic_betaS_lagrangian", &
                               test_stommel_realistic_betaS) &
                  ]
   end subroutine collect_ocean_baroclinic_stommel_tests

   subroutine collect_ocean_baroclinic_budget_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("baroclinic_budget_conservation_lagrangian", &
                               test_budget_conservation) &
                  ]
   end subroutine collect_ocean_baroclinic_budget_tests

   subroutine test_baroclinic_longrun(error)
      !! Drive the same geostrophic_adjust IC used by the
      !! `baroclinic_diag_capture` harness but with realistic
      !! `β_S = 0.78`, `VCOORD_SIGMA` for the Lagrangian + ALE
      !! remap path, and N_STEPS_LONG outer steps.  Assert
      !! no-NaN + positive h_layer + bounded `max|dS|` along
      !! the way.
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
      integer, parameter :: N_STEPS_LONG = 20000
         !! Long-run regression. With the barotropic-substep physical-wall
         !! closure fix (2026-05-19) `max|dS|` stays at FP round-off
         !! and `max|h_sum - (H+eta_end)|` stays at ~6e-14.
      integer, parameter :: PROBE_INTERVAL = 50
      real(wp), parameter :: DX = 10.0e3_wp, DY = 10.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: ETA_STEP = 0.1_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: BETA_S_REAL = 0.78_wp
      real(wp), parameter :: DS_TOLERANCE = 1.0e-6_wp
         !! Generous bound on `max|hTr_S/h - S_ref|`.  Lagrangian
         !! + ALE remap delivers ~1e-12 in practice; any drift
         !! above 1e-6 indicates the salt-feedback bug is back.
      integer :: i, j, k, step
      real(wp) :: dz, x_centre, x_i
      real(wp) :: max_dS, h_min
      logical :: hit_nan
      character(len=256) :: msg, field_msg

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE

      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA

      eos%beta_S = BETA_S_REAL
      ! Debug toggle: set this to 0.0_wp to decouple the EOS and
      ! see whether the long-run blowup is salt-feedback (then the
      ! blowup goes away) or pure-momentum (still blows up).
      ! eos%beta_S = 0.0_wp; eos%alpha_T = 0.0_wp

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

      call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call vc%enter_data()

      ! Diagnostic harness — when RDB_BCDIAG=1 is set in the
      ! environment, the in-driver probes print per-stage state.
      ! Otherwise this is a no-op (`bcdiag_enabled` stays false).
      bcdiag_S_ref = eos%S_ref
      bcdiag_step_limit = N_STEPS_LONG

      checks: block
         do step = 1, N_STEPS_LONG
            ! Cadence-gate the in-driver probes so the trace doesn't
            ! flood stdout on a healthy long run.
            bcdiag_enabled = mod(step, PROBE_INTERVAL) == 0 .or. step == 1
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf, vcoord=vc)
            bcdiag_enabled = .false.

            if (mod(step, PROBE_INTERVAL) == 0 .or. step == 1 .or. step == N_STEPS_LONG) then
               ! Pull the state from device for inspection.  Velocity
               ! and tracers are also pulled so the NaN sweep covers
               ! every prognostic field.
               !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
               !$acc&             ms%tracers(ms%idx_salinity)%hTr, &
               !$acc&             ms%tracers(ms%idx_temperature)%hTr)

               h_min = minval(ms%h_layer)

               call locate_nan(ms, field_msg, hit_nan)
               write (msg, '("baroclinic_longrun: step=", i0, " — NaN in ", a, " (prev h_min=", es10.3, ")")') &
                  step, trim(field_msg), h_min
               call check(error,.not. hit_nan, trim(msg))
               if (allocated(error)) exit checks

               write (msg, '("baroclinic_longrun: step=", i0, " — h_layer went non-positive (min = ", es10.3, ")")') step, h_min
               call check(error, h_min > 0.0_wp, trim(msg))
               if (allocated(error)) exit checks

               max_dS = max_dS_in_interior(grid, ms, eos%S_ref)
               write (msg, '("baroclinic_longrun: step=", i0, " — max|dS|=", es10.3, " > ", es10.3)') &
                  step, max_dS, DS_TOLERANCE
               call check(error, max_dS < DS_TOLERANCE, trim(msg))
               if (allocated(error)) exit checks
            end if
         end do
      end block checks

      call vc%exit_data()
      call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      call vc%destroy()
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_baroclinic_longrun

   ! ----- helpers (mirror test_ocean_analytical layout) -----

   subroutine locate_nan(ms, field_msg, has_nan)
      !! Cheap host-side sweep across every prognostic field.  Sets
      !! `has_nan = .true.` and writes the field name into
      !! `field_msg` if any NaN is found.  Caller has already pulled
      !! the fields to host via `!$acc update self`.
      type(multilayer_state_t), intent(in) :: ms
      character(len=*), intent(out) :: field_msg
      logical, intent(out) :: has_nan
      integer :: t

      has_nan = .false.
      field_msg = ""
      if (any(ieee_is_nan(ms%h_layer))) then
         has_nan = .true.; field_msg = "h_layer"; return
      end if
      if (any(ieee_is_nan(ms%u_face_x_layer))) then
         has_nan = .true.; field_msg = "u_face_x_layer"; return
      end if
      if (any(ieee_is_nan(ms%v_face_y_layer))) then
         has_nan = .true.; field_msg = "v_face_y_layer"; return
      end if
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%hTr)) cycle
            if (any(ieee_is_nan(ms%tracers(t)%hTr))) then
               has_nan = .true.
               write (field_msg, '("tracers(", i0, ")%hTr")') t
               return
            end if
         end do
      end if
   end subroutine locate_nan

   subroutine test_stommel_realistic_betaS(error)
      !! Wind-driven Stommel gyre with realistic `β_S = 0.78`,
      !! 2-layer uniform S/T, Lagrangian + ALE-remap path.
      !! Complements `test_baroclinic_longrun` (geostrophic-
      !! adjust IC) by stressing the salt feedback under a real
      !! circulation — wind-driven Ekman pumping interacting
      !! with the multilayer tracer flux.
      !!
      !! Note: layer S is uniform (= eos%S_ref) so the sigma-coord
      !! ALE remap is a no-op for tracer.  Any drift in
      !! `max|hTr_S/h - S_ref|` then comes from the multilayer
      !! continuity pipeline — exactly the path the salt-feedback
      !! bug lived in.
      !!
      !! Catches:
      !!   * Phase 1/2 regression (per-column S preservation
      !!     under multilayer continuity).
      !!   * Phase 3 regression (barotropic-substep wall closure under
      !!     a developing wind-driven barotropic mode).
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
      integer, parameter :: NX = 32, NY = 32, NZ = 2
      integer, parameter :: N_INNER = 1
      integer, parameter :: N_STEPS = 1500
      integer, parameter :: PROBE_INTERVAL = 250
      real(wp), parameter :: DX = 20.0e3_wp, DY = 20.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: F0 = 5.0e-5_wp
      real(wp), parameter :: BETA = 2.0e-11_wp
      real(wp), parameter :: TAU0 = 0.1_wp
      real(wp), parameter :: R_DRAG = 1.0e-6_wp
      real(wp), parameter :: DT = 120.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: BETA_S_REAL = 0.78_wp
      real(wp), parameter :: DS_TOLERANCE = 1.0e-6_wp
         !! Bound on `max|hTr_S/h - eos%S_ref|` across layers.  With
         !! uniform IC and Lagrangian + ALE remap, observed FP-level
         !! drift ~1e-12.  The runaway bug grew dS by orders of
         !! magnitude per few hundred steps, so 1e-6 is a wide
         !! regression bound that fires only on the actual bug.
      real(wp), parameter :: KE_MAX = 1.0e8_wp
         !! Loose bound on total KE — runaway momentum blow up
         !! produces O(1e10+).
      real(wp) :: Lx, Ly, dz, y_phys
      real(wp) :: max_dS, ke_total, h_min
      integer :: i, j, k, step, idx_S, idx_T
      logical :: hit_nan
      character(len=256) :: msg, field_msg

      Lx = real(NX, wp)*DX
      Ly = real(NY, wp)*DY

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE

      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA

      eos%beta_S = BETA_S_REAL

      checks: block
         call cor%set_beta_plane(grid, F0, BETA, y_ref=0.0_wp)

         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         do j = 1, size(ss%tau_x, 2)
            y_phys = (real(j - NGHOST, wp) - 0.5_wp)*DY
            if (y_phys > 0.0_wp .and. y_phys < Ly) then
               ss%tau_x(:, j) = -TAU0*cos(PI*y_phys/Ly)
            end if
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
         hd%kappa_h = 0.0_wp

         ! Uniform S/T per layer — the salt feedback bug surfaces
         ! through the multilayer continuity pipeline, not through
         ! IC stratification.  ALE remap is a no-op for tracer in
         ! this configuration so any drift means the bug is back.
         dz = H0/real(NZ, wp)
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
         call vc%enter_data()

         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf, vcoord=vc)

            if (mod(step, PROBE_INTERVAL) == 0 .or. step == N_STEPS) then
               !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
               !$acc&            ms%tracers(idx_S)%hTr, ms%tracers(idx_T)%hTr)

               call locate_nan(ms, field_msg, hit_nan)
               if (hit_nan) then
                  write (msg, '("stommel β_S real: NaN at step ", i0, " in ", a)') &
                     step, trim(field_msg)
                  call check(error, .false., trim(msg))
                  exit checks
               end if

               h_min = minval(ms%h_layer(NGHOST + 1:NGHOST + NX, &
                                         NGHOST + 1:NGHOST + NY, :))
               if (h_min <= 0.0_wp) then
                  write (msg, '("stommel β_S real: h_layer went non-positive at step ", &
                                &i0, " (h_min=", es10.3, ")")') step, h_min
                  call check(error, .false., trim(msg))
                  exit checks
               end if

               max_dS = max_dS_in_interior(grid, ms, eos%S_ref)
               if (max_dS > DS_TOLERANCE) then
                  write (msg, '("stommel β_S real: max|S - S_ref|=", &
                                &es10.3, " > ", es10.3, " at step ", i0)') &
                     max_dS, DS_TOLERANCE, step
                  call check(error, .false., trim(msg))
                  exit checks
               end if

               ke_total = 0.5_wp*( &
                          sum(ms%u_face_x_layer(NGHOST + 1:NGHOST + NX + 1, &
                                                NGHOST + 1:NGHOST + NY, :)**2) + &
                          sum(ms%v_face_y_layer(NGHOST + 1:NGHOST + NX, &
                                                NGHOST + 1:NGHOST + NY + 1, :)**2))
               if (ke_total > KE_MAX) then
                  write (msg, '("stommel β_S real: KE ", es10.3, " > ", es10.3, &
                                &" at step ", i0)') ke_total, KE_MAX, step
                  call check(error, .false., trim(msg))
                  exit checks
               end if
            end if
         end do

         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         call vc%exit_data()

      end block checks
      call vc%destroy()
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_stommel_realistic_betaS

   subroutine test_budget_conservation(error)
      !! Phase D v1 end-to-end check: closed-wall stratified rest IC
      !! under realistic β_S + Lagrangian + ALE remap.  The dyn-core
      !! conserves mass, salt and heat by construction (continuity is
      !! flux-form, ALE remap is conservative); the budget machinery
      !! should detect FP-level drift, no more.
      !!
      !! IC: uniform layers, stable S stratification (denser bottom).
      !! No wind, no surface flux, no horizontal viscosity.  Vertical
      !! diffusion off so tracer-T/S can't trade between layers.
      !! Bottom drag off.
      !!
      !! The test directly calls `state%budgets%init_snapshot(ms)`
      !! before the loop and `state%budgets%evaluate(ms)` after, so
      !! the cadence-fire path is bypassed (cadence is for
      !! production logging; here we just want a final residual).
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
      type(ocean_budgets_t) :: budgets
      integer, parameter :: NX = 16, NY = 8, NZ = 2
      integer, parameter :: N_INNER = 10
      integer, parameter :: N_STEPS = 500
      real(wp), parameter :: DX = 10.0e3_wp, DY = 10.0e3_wp
      real(wp), parameter :: H0 = 500.0_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: ETA_STEP = 0.1_wp
      real(wp), parameter :: BETA_S_REAL = 0.78_wp
      real(wp), parameter :: DRIFT_TOL = 1.0e-10_wp
         !! Relative drift tolerance.  Closed-wall + Lagrangian +
         !! ALE remap is bit-exact to FP for mass / salt / heat;
         !! 1e-10 leaves headroom for sum-reduction order.
      integer :: i, j, k, step, idx_S, idx_T
      real(wp) :: dz, x_centre, x_i
      real(wp) :: drift_mass, drift_S, drift_T, ref
      character(len=256) :: msg
      real(wp), allocatable :: S_layer(:)

      call grid%init(NX, NY, NGHOST, DX, DY)
      call init_slots(grid, NZ, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn, f0=F0)
      pgf%variant = OPGF_VARIANT_FV_LITE

      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA

      eos%beta_S = BETA_S_REAL
      call budgets%init(grid)

      checks: block
         ss%tau_x = 0.0_wp
         ss%tau_y = 0.0_wp
         sf%Q_heat_const = 0.0_wp
         sf%Q_salt_const = 0.0_wp
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         hv%nu_h = 0.0_wp
         hd%kappa_h = 0.0_wp
         bd%c_drag = 0.0_wp
         bd%r_linear = 0.0_wp

         dz = H0/real(NZ, wp)
         idx_T = ms%idx_temperature
         idx_S = ms%idx_salinity
         allocate (S_layer(NZ))
         S_layer = eos%S_ref
         x_centre = 0.5_wp*real(NX, wp)*DX
         ! Mild η-step (same shape as `test_baroclinic_longrun`) keeps
         ! the dyn-core out of the perfectly-uniform corner case.  The
         ! perturbation is tiny relative to H0 so conservation can be
         ! checked at FP either way; if mass / salt / heat survive
         ! 500 outer steps of evolution to FP, the budget machinery is
         ! wired correctly through the dyn driver.
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
                  ms%h_layer(i, j, k) = dz - (ETA_STEP/real(NZ, wp))* &
                                        tanh((x_i - x_centre)/(3.0_wp*DX))
                  ms%rho_layer(i, j, k) = eos%rho0
                  ms%tracers(idx_T)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
                  ms%tracers(idx_S)%hTr(i, j, k) = S_layer(k)*ms%h_layer(i, j, k)
               end do
            end do
         end do
         dyn%bt_work%bt_H_ref = H0

         ! Snapshot BEFORE map_in so the snapshot is read from host
         ! arrays (matches the evaluate path that pulls from host).
         call budgets%init_snapshot(ms)

         call map_in_slots(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         call vc%enter_data()

         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, sf=sf, vcoord=vc)
         end do

         call map_out_slots(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
         call vc%exit_data()

         call budgets%evaluate(ms)
         ref = max(budgets%values_init(BUDGET_MASS), tiny(0.0_wp))
         drift_mass = abs(budgets%values(BUDGET_MASS) - budgets%values_init(BUDGET_MASS))/ref
         ref = max(budgets%values_init(BUDGET_SALT_TOTAL), tiny(0.0_wp))
         drift_S = abs(budgets%values(BUDGET_SALT_TOTAL) - budgets%values_init(BUDGET_SALT_TOTAL))/ref
         ref = max(budgets%values_init(BUDGET_HEAT_TOTAL), tiny(0.0_wp))
         drift_T = abs(budgets%values(BUDGET_HEAT_TOTAL) - budgets%values_init(BUDGET_HEAT_TOTAL))/ref

         write (msg, '("mass drift_rel=", es10.3, " > ", es10.3, " after ", i0, " outer steps")') &
            drift_mass, DRIFT_TOL, N_STEPS
         call check(error, drift_mass < DRIFT_TOL, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("salt drift_rel=", es10.3, " > ", es10.3, " after ", i0, " outer steps")') &
            drift_S, DRIFT_TOL, N_STEPS
         call check(error, drift_S < DRIFT_TOL, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("heat drift_rel=", es10.3, " > ", es10.3, " after ", i0, " outer steps")') &
            drift_T, DRIFT_TOL, N_STEPS
         call check(error, drift_T < DRIFT_TOL, trim(msg))

      end block checks
      if (allocated(S_layer)) deallocate (S_layer)
      call vc%destroy()
      call budgets%destroy()
      call destroy_slots(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine test_budget_conservation

   real(wp) function max_dS_in_interior(grid, ms, S_ref) result(max_dS)
      !! `max|hTr_S/h_layer - S_ref|` over interior cells.  Mirrors
      !! the `probe_dS` helper in `rdb_ocean_dyn` but lives here so
      !! the test owns the assertion path.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: S_ref
      integer :: i, j, k, iS, ig, i0, i1, j0, j1
      real(wp) :: dS

      max_dS = 0.0_wp
      if (ms%idx_salinity <= 0) return
      if (.not. allocated(ms%tracers)) return

      iS = ms%idx_salinity
      ig = grid%nghost
      i0 = ig + 1
      i1 = grid%nx_total - ig
      j0 = ig + 1
      j1 = grid%ny_total - ig

      do k = 1, ms%nz_ml
         do j = j0, j1
            do i = i0, i1
               if (ms%h_layer(i, j, k) > 0.0_wp) then
                  dS = abs(ms%tracers(iS)%hTr(i, j, k)/ms%h_layer(i, j, k) - S_ref)
                  if (dS > max_dS) max_dS = dS
               end if
            end do
         end do
      end do
   end function max_dS_in_interior

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
      call ms%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, dyn)
      !$acc exit data delete(ms)
   end subroutine map_out_slots

end module test_ocean_baroclinic_longrun
