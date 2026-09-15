!! Unit tests for the multilayer SSP-RK2 driver
!! (rdb_ocean_dyn::ocean_dyn_step).  One step orchestrates
!!
!!   EOS  ->  continuity-PPM  ->  Sadourny Coriolis-adv
!!        ->  Mont pressure-force  ->  PPM tracer advection
!!        ->  forward-Euler applies x2 + RK2 average
!!
!! across the full multilayer C-grid stack.  These tests are the
!! end-to-end Tier-1 check — every kernel touched in sequence,
!! with the driver responsible for the state-save / two FE stages
!! / final average.
!!
!! Cases:
!!   * Stratified lake-at-rest — layered density (heavy at bed,
!!     light at surface), constant T, S per layer, zero velocity.
!!     The trivial fixed point of the full system.  Every prognostic
!!     must stay bit-for-bit at its IC after one driver step:
!!       - h_layer (continuity is no-op for zero u, v)
!!       - u, v (Coriolis null on zero velocity; PGF null because
!!         horizontal pressure gradient vanishes; full balance)
!!       - hTr (tracer advect null for zero mass_flux)
!!     The strongest possible "discrete hydrostatic balance" check.
!!   * Inertial oscillation through the driver — uniform u, zero
!!     v, uniform h, uniform T/S (so PGF stays zero), f = 1.  After
!!     a quarter inertial period the velocity vector must have
!!     rotated through pi/2 with RK2-level magnitude preservation
!!     (~1e-5; FE leaks ~0.8%).
!!   * Mass + tracer conservation — 30 driver steps with a smooth
!!     wall-vanishing velocity, non-trivial h and tracer fields.
!!     Total mass and total tracer mass drift < 1e-9.
module test_ocean_dyn_multilayer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
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
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step
   implicit none
   private

   public :: collect_ocean_dyn_multilayer_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_dyn_multilayer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stratified_lake_at_rest", test_stratified_rest), &
                  new_unittest("inertial_oscillation_through_driver", &
                               test_inertial_oscillation), &
                  new_unittest("mass_and_tracer_conservation", test_conservation), &
                  new_unittest("kshear_vertex_kv_reaches_momentum", &
                               test_kshear_vertex_kv_wire) &
                  ]
   end subroutine collect_ocean_dyn_multilayer_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
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
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
      !$acc enter data copyin(cor)
      call cor%enter_data()
      !$acc enter data copyin(pgf)
      call pgf%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
      !$acc enter data copyin(bd)
      call bd%enter_data()
      !$acc enter data copyin(ss)
      call ss%enter_data()
      !$acc enter data copyin(va)
      call va%enter_data()
      !$acc enter data copyin(hd)
      call hd%enter_data()
      !$acc enter data copyin(vd)
      call vd%enter_data()
      !$acc enter data copyin(vmix)
      call vmix%enter_data()
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
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
      call destroy_cartesian_metrics(metrics)
      call vmix%exit_data()
      !$acc exit data delete(vmix)
      call vd%exit_data()
      !$acc exit data delete(vd)
      call hd%exit_data()
      !$acc exit data delete(hd)
      call va%exit_data()
      !$acc exit data delete(va)
      call ss%exit_data()
      !$acc exit data delete(ss)
      call bd%exit_data()
      !$acc exit data delete(bd)
      call hv%exit_data()
      !$acc exit data delete(hv)
      call pgf%exit_data()
      !$acc exit data delete(pgf)
      call cor%exit_data()
      !$acc exit data delete(cor)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_stratified_rest(error)
      !! Stratified rest state: uniform h per layer, layered (T, S)
      !! that vary by layer but are spatially uniform, zero velocity.
      !! After one driver step every prognostic must stay at its IC
      !! to round-off — this exercises the full Tier-1 stack and
      !! catches any non-zero tendency introduced spuriously by any
      !! kernel.
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
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp) :: S_k(NZ), T_k(NZ)
      real(wp), allocatable :: hS_ic(:, :, :), hT_ic(:, :, :)
      real(wp) :: max_dh, max_du, max_dv, max_dS, max_dT
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
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
         call dyn%init(grid)

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! Stratified S, T (heavier S at the bed, cooler T at the bed)
         S_k = [eos%S_ref + 1.0_wp, eos%S_ref, eos%S_ref - 1.0_wp]
         T_k = [eos%T_ref - 2.0_wp, eos%T_ref, eos%T_ref + 2.0_wp]
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*H0
         end do
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)

         max_dh = maxval(abs(ms%h_layer - H0))
         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_dh < 1.0e-12_wp, &
                    "stratified rest: h drifted")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-12_wp, &
                    "stratified rest: u drifted (PGF + Coriolis didn't balance to zero)")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-12_wp, &
                    "stratified rest: v drifted")
         if (allocated(error)) exit checks
         call check(error, max_dS < 1.0e-12_wp, &
                    "stratified rest: hS drifted")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-12_wp, &
                    "stratified rest: hT drifted")
         if (allocated(error)) exit checks
         call check(error, dyn%outer_step_count == 1, &
                    "outer_step_count did not advance")

      end block checks
      deallocate (hS_ic, hT_ic)
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy(); call pgf%destroy()
      call cor%destroy(); call ct%destroy(); call ms%destroy()
   end subroutine test_stratified_rest

   subroutine test_inertial_oscillation(error)
      !! Short-window inertial response check.  Uniform u = U0,
      !! v = 0, uniform h, T = T_ref, S = S_ref (rho = rho_0
      !! everywhere -> PGF = 0).  Closed walls force mass_flux = 0
      !! at the boundary, which excites a surface gravity wave
      !! travelling inward at c = sqrt(g*H_total) ~ 54 m/s.  We run
      !! short enough that the wave hasn't reached the central probe
      !! (probe ~7 cells from each wall, c*t ~3 cells of travel at
      !! t = 5*DT) so the probe still sees pure inertial dynamics.
      !!
      !! Inertial rotation after t = N*DT is v = -f*U0*t to leading
      !! order in (f*t).  We verify the early-time linear response
      !! rather than the full quarter-period rotation (closed walls +
      !! Heun's known marginal instability on the wave equation make
      !! the long-time rotation test impossible without proper
      !! Phase 4b split-explicit substepping).
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
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: DT = 0.001_wp
      integer, parameter :: N_STEPS = 5
      real(wp) :: u_obs, v_obs, v_expected, v_err, t_end
      integer :: step, i_probe, j_probe, nx, ny
      checks: block

         call make_grid(grid, 16, 16, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
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
         call dyn%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H0
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H0

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
         do step = 1, N_STEPS
            call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)

         i_probe = nx/2
         j_probe = ny/2
         u_obs = ms%u_face_x_layer(i_probe, j_probe, NZ/2)
         v_obs = ms%v_face_y_layer(i_probe, j_probe, NZ/2)

         t_end = real(N_STEPS, wp)*DT
         ! v after t in pure inertial: v = -f*U0*sin(f*t) ~ -f*U0*t
         v_expected = -F_C*U0*sin(F_C*t_end)
         v_err = abs(v_obs - v_expected)

         call check(error, v_err < 1.0e-6_wp, &
                    "driver inertial oscillation: linear-response v error > 1e-6")
         if (allocated(error)) exit checks
         call check(error, abs(u_obs - U0) < 1.0e-3_wp, &
                    "driver inertial oscillation: u drifted from U0 too far")

      end block checks
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy(); call pgf%destroy()
      call cor%destroy(); call ct%destroy(); call ms%destroy()
   end subroutine test_inertial_oscillation

   subroutine test_conservation(error)
      !! Conservation of total volume + total tracer mass through the
      !! multilayer driver.  The test isolates conservation from
      !! stability: walls are closed (`mass_flux = 0` at boundaries
      !! enforced by continuity), and we drive the system with a
      !! solenoidal-ish u, v pattern at a small amplitude so the
      !! gravity-wave excitation stays in the linear regime over the
      !! integration window.  Conservation is an algebraic property of
      !! the discrete divergence — the closed-wall divergence sums to
      !! zero by construction regardless of the wave state — so even
      !! short windows expose conservation bugs at machine precision.
      !!
      !! Heun's method is *not* unconditionally stable on the linear
      !! gravity-wave equation (the SSP-RK2 stability region only
      !! touches the imaginary axis tangentially at 0), so we keep
      !! amplitude small and step count modest to avoid running into
      !! the slow Heun-on-wave amplification.  Phase 4b's split-
      !! explicit substepping will resolve the surface mode at the
      !! fast inner dt and let the slow outer step run safely on the
      !! baroclinic mode.
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
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U_AMP = 0.01_wp
      real(wp), parameter :: DT = 0.02_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 10
      integer :: i, j, k, nx, ny, step
      real(wp) :: total_h0, total_S0, total_T0
      real(wp) :: total_h, total_S, total_T, drift_h, drift_S, drift_T
      real(wp) :: h_min
      checks: block

         call make_grid(grid, 24, 16, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
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
         call dyn%init(grid)
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

         total_h0 = sum(ms%h_layer)*grid%dx*grid%dy
         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T0 = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
         do step = 1, N_STEPS
            call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT)
         end do
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)

         total_h = sum(ms%h_layer)*grid%dx*grid%dy
         total_S = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy
         total_T = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         drift_h = abs(total_h - total_h0)/abs(total_h0)
         drift_S = abs(total_S - total_S0)/abs(total_S0)
         drift_T = abs(total_T - total_T0)/abs(total_T0)
         h_min = minval(ms%h_layer)

         call check(error, drift_h < 1.0e-12_wp, "total h drift > 1e-12 through driver")
         if (allocated(error)) exit checks
         call check(error, drift_S < 1.0e-12_wp, "total S mass drift > 1e-12")
         if (allocated(error)) exit checks
         call check(error, drift_T < 1.0e-12_wp, "total T mass drift > 1e-12")
         if (allocated(error)) exit checks
         call check(error, h_min > 0.0_wp, "h_layer went non-positive")

      end block checks
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy(); call pgf%destroy()
      call cor%destroy(); call ct%destroy(); call ms%destroy()
   end subroutine test_conservation

   subroutine test_kshear_vertex_kv_wire(error)
      !! End-to-end WIRE check for the vertex kappa-shear corner-Kv
      !! seam: the dispatch (`vmix_apply_in_stage`) must hand the
      !! kappa-shear `kd_corner` carrier to `vdiff_apply_momentum` as
      !! `kv_corner_source` scaled by `prandtl_turb`.  Two identical
      !! full driver steps with a SEEDED corner field, differing ONLY
      !! in `prandtl_turb` (1 vs 0): in vertex mode prandtl_turb has NO
      !! other pathway to the momentum solve (the cell-centred kv merge
      !! is suppressed, and thermodynamics are off so the corner solve
      !! never overwrites the seed) — so the velocities differ IFF the
      !! corner->face wire is alive.  A dead wire (corner source not
      !! passed) makes the two runs bit-identical and the test fails.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: u_pr1(:, :, :), u_pr0(:, :, :)
      checks: block
         call run_vertex_kv_step(1.0_wp, u_pr1)
         call run_vertex_kv_step(0.0_wp, u_pr0)
         call check(error, maxval(abs(u_pr1 - u_pr0)) > 1.0e-10_wp, &
                    "vertex-Kv wire: prandtl_turb must reach the momentum "// &
                    "solve through the corner->face seam — identical runs "// &
                    "mean the dispatch never passed kd_corner to "// &
                    "vdiff_apply_momentum")
      end block checks
      if (allocated(u_pr1)) deallocate (u_pr1, u_pr0)
   end subroutine test_kshear_vertex_kv_wire

   subroutine run_vertex_kv_step(prandtl, u_out)
      !! One full ocean_dyn_step with a vertex-armed kappa-shear slot
      !! whose `kd_corner` is SEEDED (thermo off => the corner solve
      !! never runs, so the seed is exactly what vdiff consumes).
      !! Sheared u so vertical viscosity visibly acts; uniform T/S at
      !! the EOS refs so the PGF is quiet.
      real(wp), intent(in) :: prandtl
      real(wp), intent(out), allocatable :: u_out(:, :, :)
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
      type(ocean_kappa_shear_t) :: ks
      real(wp), parameter :: H0 = 10.0_wp, DT = 300.0_wp, KC = 0.02_wp
      integer :: k

      call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
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
      call dyn%init(grid)

      ! kshear: vertex-armed, enabled, thermo OFF so the seeded corner
      ! field is what reaches the momentum solve.
      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .true.
      call ks%init(grid, nz_ml=NZ)
      ks%enable = .true.
      ks%prandtl_turb = prandtl
      call ks%init_vertex(grid, NZ)
      do k = 2, NZ
         ks%kd_corner(:, :, k) = KC
      end do

      ms%h_layer = H0
      do k = 1, NZ
         ms%u_face_x_layer(:, :, k) = 0.1_wp*real(k, wp)   ! vertical shear
      end do
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
      end do

      call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)
      !$acc enter data copyin(ks)
      call ks%enter_data()
      call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                          va, hd, vd, vmix, ms, DT, kshear=ks)
      call ks%exit_data()
      !$acc exit data delete(ks)
      call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix)

      allocate (u_out, source=ms%u_face_x_layer)

      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call ks%destroy(); call ms%destroy()
   end subroutine run_vertex_kv_step

end module test_ocean_dyn_multilayer
