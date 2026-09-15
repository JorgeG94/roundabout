!! Atmospheric surface-pressure loading / inverse barometer (PR-17).
!!
!! Three properties, two tiers:
!!  (a) SEAM FOLD (unit, no time-stepping) — `p_surf_update_seam` gives
!!      `eta_ib = -p_surf/(rho0 g)` to round-off with the CORRECT (negative)
!!      sign, `eta_seam = eta_ib` (no tide) / `eta_ib + eta_tide` (tide), and
!!      the ~1 cm/hPa magnitude.  Carries its own NEGATIVE CONTROL: the
!!      positive-sign reference must NOT match (the minus sign is the one
!!      thing the roadmap got wrong).
!!  (b) INVERSE-BAROMETER EQUILIBRIUM (integration) — on a closed Cartesian
!!      basin with uniform density, initialising the SSH to the analytical
!!      response `eta = eta_ib = -p_surf/(rho0 g)` is a discrete FIXED POINT:
!!      the seam forcing `-g grad(eta - eta_seam)` cancels to round-off, so
!!      the flow stays quiescent.  The SAME p_surf started from a FLAT SSH is
!!      strongly forced.  The contrast (equilibrium quiet << flat loud) is
!!      the defining inverse-barometer property AND its own negative control:
!!      with the wrong sign the equilibrium IC would be the LOUD one.
!!  (c) GAUGE INVARIANCE (integration) — a spatially UNIFORM p_surf has no
!!      gradient, so it is byte-for-byte identical to p_surf = 0 (only
!!      grad(p_surf) is physical).
!!
!! GPU/mem:separate: every state object AND its scratch companion is
!! enter_data'd before the first kernel; host-set inputs (sf%p_surf) are
!! update device'd via sf%enter_data; fields read back are update self'd.
!! All directives are inert on the multicore build — a green multicore run
!! proves nothing about device data motion; verify on the GPU build.
!! Single-rank: the ghost-band halo exchange of eta_seam is a no-op here, so
!! this test cannot see a missing exchange (it is a multi-rank-only bug).
module test_ocean_p_surf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics
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
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_p_surf, only: ocean_p_surf_t, p_surf_configure, p_surf_update_seam
   implicit none
   private

   public :: collect_ocean_p_surf_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 20, NY_PHYS = 20, NZ = 2
   real(wp), parameter :: DX = 50000.0_wp, DY = 50000.0_wp
   real(wp), parameter :: H_TOTAL = 4000.0_wp
   real(wp), parameter :: DT = 300.0_wp
   integer, parameter :: N_INNER = 16
   integer, parameter :: N_STEPS = 20
   real(wp), parameter :: RHO0 = 1035.0_wp
   integer, parameter :: N_RESP = 3               ! steps for the ramp response
   ! Linear ramp: total 2000 Pa (20 hPa) across the physical domain width.
   real(wp), parameter :: GRAD_P = 2000.0_wp/(real(NX_PHYS, wp)*DX)

contains

   subroutine collect_ocean_p_surf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("psurf_seam_fold", test_seam_fold), &
                  new_unittest("inverse_barometer_response", test_ib_response), &
                  new_unittest("gauge_invariance_uniform_p_surf", test_gauge_invariance)]
   end subroutine collect_ocean_p_surf_tests

   ! ---- (a) unit-level seam fold: sign + magnitude + composition ----
   subroutine test_seam_fold(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_p_surf_t) :: psurf
      integer :: nx, ny, i, j
      real(wp), allocatable :: p_surf(:, :), eta_tide(:, :)
      real(wp) :: g_bt, i_rho0_g, ref, maxdev, maxdev_wrongsign, one_hpa_ib

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      g_bt = 9.81_wp
      i_rho0_g = 1.0_wp/(RHO0*g_bt)

      psurf%enable = .true.
      psurf%rho0 = RHO0
      call p_surf_configure(psurf, nx, ny)

      allocate (p_surf(nx, ny), eta_tide(nx, ny))
      do j = 1, ny
         do i = 1, nx
            ! A spatially varying, strictly positive (>=0) surface pressure.
            p_surf(i, j) = 1.0e5_wp + 50.0_wp*real(i, wp) - 30.0_wp*real(j, wp)
            eta_tide(i, j) = 0.001_wp*real(i - j, wp)
         end do
      end do

      call psurf%enter_data()
      !$acc enter data copyin(p_surf, eta_tide)
      !$acc update device(p_surf, eta_tide)

      ! --- no tide: eta_seam == eta_ib == -p_surf/(rho0 g) ---
      call p_surf_update_seam(psurf, p_surf, g_bt)
      !$acc update self(psurf%eta_ib, psurf%eta_seam)
      maxdev = 0.0_wp
      maxdev_wrongsign = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            ref = -p_surf(i, j)*i_rho0_g
            maxdev = max(maxdev, abs(psurf%eta_ib(i, j) - ref))
            ! NEGATIVE CONTROL: the positive-sign (roadmap) reference.
            maxdev_wrongsign = max(maxdev_wrongsign, &
                                   abs(psurf%eta_ib(i, j) - (-ref)))
            maxdev = max(maxdev, abs(psurf%eta_seam(i, j) - psurf%eta_ib(i, j)))
         end do
      end do
      call check(error, maxdev < 1.0e-14_wp, &
                 "eta_ib must equal -p_surf/(rho0*g_bt) and eta_seam==eta_ib")
      if (allocated(error)) then
         call cleanup_fold(psurf, p_surf, eta_tide)
         return
      end if
      ! Teeth: the wrong (positive) sign must be VISIBLY off (p_surf ~ 1e5 Pa
      ! => eta_ib ~ 10 m, so a sign flip is a ~20 m discrepancy).
      call check(error, maxdev_wrongsign > 1.0_wp, &
                 "NEGATIVE CONTROL: +p_surf/(rho0 g) sign must NOT match")
      if (allocated(error)) then
         call cleanup_fold(psurf, p_surf, eta_tide)
         return
      end if

      ! --- ~1 cm per hPa magnitude, explicit ---
      one_hpa_ib = 100.0_wp*i_rho0_g          ! |eta_ib| for 1 hPa = 100 Pa
      call check(error, abs(one_hpa_ib - 0.00985_wp) < 5.0e-5_wp, &
                 "inverse barometer must be ~1 cm (9.85 mm) per hPa")
      if (allocated(error)) then
         call cleanup_fold(psurf, p_surf, eta_tide)
         return
      end if

      ! --- with tide: eta_seam == eta_ib + eta_tide (exact additive) ---
      call p_surf_update_seam(psurf, p_surf, g_bt, eta_tide=eta_tide)
      !$acc update self(psurf%eta_ib, psurf%eta_seam)
      maxdev = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            ref = psurf%eta_ib(i, j) + eta_tide(i, j)
            maxdev = max(maxdev, abs(psurf%eta_seam(i, j) - ref))
         end do
      end do
      call check(error, maxdev < 1.0e-14_wp, &
                 "eta_seam must equal eta_ib + eta_tide (additive composition)")

      call cleanup_fold(psurf, p_surf, eta_tide)
   end subroutine test_seam_fold

   subroutine cleanup_fold(psurf, p_surf, eta_tide)
      type(ocean_p_surf_t), intent(inout) :: psurf
      real(wp), allocatable, intent(inout) :: p_surf(:, :), eta_tide(:, :)
      !$acc exit data delete(p_surf, eta_tide)
      call psurf%exit_data()
      call psurf%destroy()
      deallocate (p_surf, eta_tide)
   end subroutine cleanup_fold

   ! ---- (b) inverse-barometer response: acceleration = -(1/rho0) grad p ----
   subroutine test_ib_response(error)
      !! From REST under a linear surface-pressure ramp `p_surf = P0 +
      !! GRAD_P*x`, the depth-uniform inverse-barometer acceleration
      !! `a = -(1/rho0) dp/dx` drives the barotropic mode.  Over N_RESP outer
      !! steps from rest (before the SSH adjusts enough to oppose it) the
      !! net velocity impulse is the analytical `u = -N_RESP*dt*GRAD_P/rho0`,
      !! independent of basin geometry (a constant gradient => a uniform
      !! response, no seiche).  Pins the SIGN (a high-pressure gradient drives
      !! flow DOWN-gradient, toward the low) and the 1/rho0 magnitude
      !! end-to-end through the real split solver.  A UNIFORM p_surf control
      !! (zero gradient) must stay at rest — the built-in negative control.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_ramp, u_ctrl, u_expect

      call run_ramp(error, grad_p=GRAD_P, mean_u=u_ramp)
      if (allocated(error)) return
      call run_ramp(error, grad_p=0.0_wp, mean_u=u_ctrl)
      if (allocated(error)) return

      u_expect = -real(N_RESP, wp)*DT*GRAD_P/RHO0

      ! Negative control: no pressure GRADIENT => no motion.
      call check(error, abs(u_ctrl) < 1.0e-6_wp, &
                 "uniform p_surf (zero gradient) must not drive flow")
      if (allocated(error)) return
      ! The ramp must drive flow with the CORRECT sign (down-gradient).
      call check(error, u_ramp*u_expect > 0.0_wp .and. abs(u_ramp) > 1.0e-6_wp, &
                 "surface-pressure ramp drove flow the WRONG way (sign error "// &
                 "in the inverse-barometer fold) or not at all")
      if (allocated(error)) return
      ! The 1/rho0 magnitude, end-to-end: the mean response matches the
      ! analytical impulse -N*dt*grad(p)/rho0 to 30%.  The response is
      ! slightly BELOW the free impulse (the SSH's early back-reaction
      ! opposes the forcing over N_RESP steps); the bound still rejects any
      ! gross error (a wrong sign is 100% off, a missing 1/rho0 is ~rho0 x).
      call check(error, abs(u_ramp - u_expect) < 0.30_wp*abs(u_expect), &
                 "inverse-barometer response magnitude != -dt*grad(p)/rho0")
   end subroutine test_ib_response

   ! ---- (c) uniform p_surf is byte-for-byte identical to p_surf = 0 ----
   subroutine test_gauge_invariance(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: eta_u(:, :), u_u(:, :, :), v_u(:, :, :)
      real(wp), allocatable :: eta_0(:, :), u_0(:, :, :), v_0(:, :, :)
      real(wp) :: dmax

      call run_uniform(error, p_const=5.0e4_wp, eta_out=eta_u, u_out=u_u, v_out=v_u)
      if (allocated(error)) return
      call run_uniform(error, p_const=0.0_wp, eta_out=eta_0, u_out=u_0, v_out=v_0)
      if (allocated(error)) return

      dmax = maxval(abs(eta_u - eta_0))
      dmax = max(dmax, maxval(abs(u_u - u_0)))
      dmax = max(dmax, maxval(abs(v_u - v_0)))
      call check(error, dmax == 0.0_wp, &
                 "a spatially UNIFORM p_surf must be byte-for-byte identical "// &
                 "to p_surf = 0 (gauge invariance: only grad(p_surf) is physical)")

      deallocate (eta_u, u_u, v_u, eta_0, u_0, v_0)
   end subroutine test_gauge_invariance

   ! ============================================================
   ! Linear-ramp response driver.
   ! ============================================================
   subroutine run_ramp(error, grad_p, mean_u)
      !! Cartesian basin at REST, uniform density, no Coriolis, with a
      !! linear surface-pressure ramp `p_surf = P0 + grad_p*(x - x0)`.
      !! Integrates N_RESP outer steps and returns the mean signed interior
      !! `u_face_x_layer` — the barotropic response to the depth-uniform
      !! inverse-barometer acceleration.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: grad_p
      real(wp), intent(out) :: mean_u
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_p_surf_t) :: psurf
      integer :: nx, ny, i, j, k, step, i0, ni
      real(wp) :: usum

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      i0 = grid%nghost + NX_PHYS/2
      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 0.0_wp                     ! no Coriolis: clean 1-D response
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
      dyn%bt_work%bt_H_ref = H_TOTAL

      ! Linear surface-pressure ramp in x (p0 arbitrary; only the constant
      ! gradient is physical).  Strictly positive => a valid load.
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      do j = 1, ny
         do i = 1, nx
            sf%p_surf(i, j) = 1.0e5_wp + grad_p*real(i - i0, wp)*DX
         end do
      end do

      ! Uniform density, uniform S/T, flat SSH, at REST.
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%h_layer(:, :, k) = H_TOTAL/real(NZ, wp)
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*(H_TOTAL/real(NZ, wp))
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*(H_TOTAL/real(NZ, wp))
      end do

      psurf%enable = .true.
      psurf%rho0 = eos%rho0
      call p_surf_configure(psurf, nx, ny)

      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, sf, psurf)
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
      call sf%enter_data()          ! maps + update device p_surf*
      call psurf%enter_data()

      do step = 1, N_RESP
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, &
                                   sf=sf, psurf=psurf)
      end do

      !$acc update self(ms%u_face_x_layer)
      ! Mean signed u over interior u-faces (away from the walls).
      usum = 0.0_wp
      ni = 0
      do k = 1, NZ
         do j = grid%nghost + 2, grid%nghost + NY_PHYS - 1
            do i = grid%nghost + 2, grid%nghost + NX_PHYS - 1
               usum = usum + ms%u_face_x_layer(i, j, k)
               ni = ni + 1
            end do
         end do
      end do
      mean_u = usum/real(max(ni, 1), wp)

      call sf%exit_data()
      call psurf%exit_data()
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
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, sf, psurf)
      call metrics%exit_data()

      call check(error, mean_u == mean_u, "mean_u is NaN")

      call psurf%destroy()
      call sf%destroy()
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
   end subroutine run_ramp

   subroutine run_uniform(error, p_const, eta_out, u_out, v_out)
      !! Same basin but a UNIFORM p_surf = p_const and a non-trivial IC
      !! (a tilted SSH + a sinusoidal velocity), integrated N_STEPS.  Returns
      !! the final barotropic SSH + layer velocities for a byte comparison.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: p_const
      real(wp), allocatable, intent(out) :: eta_out(:, :), u_out(:, :, :), v_out(:, :, :)
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
      type(ocean_surface_flux_t) :: sf
      type(ocean_p_surf_t) :: psurf
      integer :: nx, ny, i, j, k, step
      real(wp) :: tilt

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 1.0e-4_wp
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
      dyn%bt_work%bt_H_ref = H_TOTAL

      call sf%init(grid)
      call sf%set_components(grid, .true.)
      sf%p_surf = p_const

      ! Non-trivial IC: a gently tilted SSH + a sinusoidal u so the run is
      ! dynamically active (a trivial rest state could hide a leak).
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               tilt = 0.05_wp*real(i - nx/2, wp)/real(nx, wp)
               ms%h_layer(i, j, k) = (H_TOTAL + tilt)/real(NZ, wp)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  eos%S_ref*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  eos%T_ref*ms%h_layer(i, j, k)
               ms%u_face_x_layer(i, j, k) = &
                  0.02_wp*sin(2.0_wp*acos(-1.0_wp)*real(j, wp)/real(ny, wp))
               ms%v_face_y_layer(i, j, k) = 0.0_wp
            end do
         end do
      end do

      psurf%enable = .true.
      psurf%rho0 = eos%rho0
      call p_surf_configure(psurf, nx, ny)

      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, sf, psurf)
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
      call sf%enter_data()
      call psurf%enter_data()

      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, &
                                   sf=sf, psurf=psurf)
      end do

      !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer, dyn%bt_work%bt_eta)
      ! C-grid faces are staggered (u has nx+1 in dim 1, v has ny+1 in
      ! dim 2), so take the arrays' own shapes via source-only allocation.
      allocate (eta_out, source=dyn%bt_work%bt_eta)
      allocate (u_out, source=ms%u_face_x_layer)
      allocate (v_out, source=ms%v_face_y_layer)

      call sf%exit_data()
      call psurf%exit_data()
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
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, sf, psurf)
      call metrics%exit_data()

      call check(error, all(eta_out == eta_out), "eta_out is NaN")

      call psurf%destroy()
      call sf%destroy()
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
   end subroutine run_uniform

end module test_ocean_p_surf
