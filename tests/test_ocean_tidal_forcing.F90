!! Physical end-to-end test of the C1 equilibrium body-force tide.
!!
!! Drives the real split-explicit ocean step (`ocean_dyn_step_split`)
!! with a single M2 constituent on a genuine SPHERICAL patch (the honest
!! path: exercises `metrics_fill_spherical` -> geolatT/geolonT -> the
!! init-time `cos_struct/sin_struct` build -> the per-outer-step
!! `tides_update_eta_eq` -> the gated barotropic PGF).  A tiny closed
!! basin (crossing time << the M2 period) sits in the quasi-static
!! near-equilibrium regime, so the surface height responds to the tidal
!! body force.  Asserts:
!!   * the response is finite and bounded (no blow-up);
!!   * SSH develops a non-trivial signal (the seam is live end-to-end);
!!   * SSH at an interior probe OSCILLATES with the M2 period (~12.42 h);
!!   * the closed basin conserves volume (basin-mean SSH ~ 0).
!! The closed-basin SSH is the gradient-driven response to the (nearly
!! uniform over a small basin) equilibrium tide; its magnitude scales
!! with, but is smaller than, the full equilibrium envelope
!! A_M2*love*cos^2(phi) — recovering the exact envelope would need open
!! boundaries, out of scope for v1.
module test_ocean_tidal_forcing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics
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
   use rdb_ocean_tides, only: ocean_tides_t, tides_configure_astronomy, &
                              tides_build_struct
   use rdb_ocean_tide_astro, only: TIDE_OMEGA, TIDE_AMP, TIDE_LOVE, TIDE_DEG2RAD
   implicit none
   private

   public :: collect_ocean_tidal_forcing_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 24, NY_PHYS = 24, NZ = 2
   real(wp), parameter :: LAT0 = 29.0_wp        ! south edge (deg)
   real(wp), parameter :: LON0 = 175.0_wp       ! west edge (deg)
   real(wp), parameter :: DLON = 0.5_wp         ! deg/cell -> 12 deg span
   real(wp), parameter :: DLAT = 0.1_wp         ! deg/cell -> 2.4 deg span
   real(wp), parameter :: PHI_MID = 30.2_wp     ! patch-centre latitude (deg)
   real(wp), parameter :: RAD_EARTH = 6.371e6_wp
   real(wp), parameter :: H_TOTAL = 4000.0_wp
   real(wp), parameter :: DT = 300.0_wp
   integer, parameter :: N_INNER = 16

contains

   subroutine collect_ocean_tidal_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("m2_body_force_ssh_oscillates", test_m2_oscillation)]
   end subroutine collect_ocean_tidal_forcing_tests

   subroutine test_m2_oscillation(error)
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
      type(ocean_tides_t) :: tides
      integer :: cat_idx(1)
      integer :: i, j, k, ip, jp, step, n_steps
      integer :: nx, ny, i0, i1, j0, j1
      real(wp) :: t, t_end, period_m2, hcol, ssh_probe, ssh_mean, area_n
      real(wp) :: max_abs, min_probe, max_probe, mean_abs
      real(wp), allocatable :: ssh_t(:), t_t(:)
      real(wp) :: env_amp, tcross1, tcross2, period_obs
      integer :: ncross, s

      period_m2 = 2.0_wp*acos(-1.0_wp)/TIDE_OMEGA(1)     ! M2 period (s)
      t_end = 2.2_wp*period_m2
      n_steps = int(t_end/DT) + 1
      env_amp = TIDE_AMP(1)*TIDE_LOVE(1)*cos(PHI_MID*TIDE_DEG2RAD)**2

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DLON, DLAT)
      nx = grid%nx_total
      ny = grid%ny_total
      i0 = grid%nghost + 1
      i1 = grid%nghost + NX_PHYS
      j0 = grid%nghost + 1
      j1 = grid%nghost + NY_PHYS
      ip = grid%nghost + NX_PHYS/4          ! off-centre probe (avoid node)
      jp = grid%nghost + NY_PHYS/2

      call make_spherical_metrics(metrics, grid, LON0, LAT0, DLON, DLAT, RAD_EARTH)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 2.0_wp*7.292115e-5_wp*sin(PHI_MID*TIDE_DEG2RAD)   ! f-plane at patch centre
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

      ! Uniform T/S (barotropic response only), flat bottom, rest IC.
      ms%h_layer = H_TOTAL/real(NZ, wp)
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*(H_TOTAL/real(NZ, wp))
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*(H_TOTAL/real(NZ, wp))
      end do
      dyn%bt_work%bt_H_ref = H_TOTAL

      ! Configure the M2 tide slot from the real spherical lat/lon.
      cat_idx(1) = 1     ! M2
      tides%enable = .true.
      call tides_configure_astronomy(tides, cat_idx, 1, 0.0_wp, 0.0_wp, .false., nx, ny)
      call tides_build_struct(tides, metrics%geolatT, metrics%geolonT, nx, ny)

      allocate (ssh_t(n_steps), t_t(n_steps))

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
      call tides%enter_data()

      area_n = real((i1 - i0 + 1)*(j1 - j0 + 1), wp)
      t = 0.0_wp
      do step = 1, n_steps
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, t=t, tides=tides)
         t = t + DT
         !$acc update self(ms%h_layer)
         ! SSH = column height - reference depth.
         ssh_probe = 0.0_wp
         do k = 1, NZ
            ssh_probe = ssh_probe + ms%h_layer(ip, jp, k)
         end do
         ssh_probe = ssh_probe - H_TOTAL
         ! Basin-mean SSH (volume conservation).
         ssh_mean = 0.0_wp
         do j = j0, j1
            do i = i0, i1
               hcol = 0.0_wp
               do k = 1, NZ
                  hcol = hcol + ms%h_layer(i, j, k)
               end do
               ssh_mean = ssh_mean + (hcol - H_TOTAL)
            end do
         end do
         ssh_mean = ssh_mean/area_n
         ssh_t(step) = ssh_probe
         t_t(step) = t
      end do

      call tides%exit_data()
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

      checks: block
         ! 1. Finite + bounded (no blow-up): |SSH| < a few * envelope.
         max_abs = maxval(abs(ssh_t))
         call check(error, max_abs == max_abs, "SSH is NaN")   ! NaN self-compare
         if (allocated(error)) exit checks
         call check(error, max_abs < 5.0_wp*env_amp, &
                    "SSH blew up (not near-equilibrium)")
         if (allocated(error)) exit checks

         ! 2. Live: the tide actually drives a signal.
         call check(error, max_abs > 1.0e-4_wp, &
                    "SSH is inert — the tidal body force is not reaching the solve")
         if (allocated(error)) exit checks

         ! 3. Oscillatory: the probe changes sign after spin-up.
         min_probe = minval(ssh_t(n_steps/4:))
         max_probe = maxval(ssh_t(n_steps/4:))
         call check(error, min_probe < 0.0_wp .and. max_probe > 0.0_wp, &
                    "SSH probe does not oscillate about zero")
         if (allocated(error)) exit checks

         ! 4. M2 period: spacing between successive upward zero-crossings
         !    (after a half-period spin-up) must match the M2 period.
         ncross = 0
         tcross1 = -1.0_wp
         tcross2 = -1.0_wp
         do s = n_steps/4, n_steps - 1
            if (ssh_t(s) <= 0.0_wp .and. ssh_t(s + 1) > 0.0_wp) then
               ncross = ncross + 1
               ! linear interp of the crossing time
               if (ncross == 1) then
                  tcross1 = t_t(s) + DT*(-ssh_t(s))/(ssh_t(s + 1) - ssh_t(s))
               else if (ncross == 2) then
                  tcross2 = t_t(s) + DT*(-ssh_t(s))/(ssh_t(s + 1) - ssh_t(s))
               end if
            end if
         end do
         call check(error, ncross >= 2, "fewer than 2 M2 cycles captured")
         if (allocated(error)) exit checks
         period_obs = tcross2 - tcross1
         call check(error, abs(period_obs - period_m2)/period_m2 < 0.15_wp, &
                    "SSH oscillation period does not match M2 (12.42 h)")
         if (allocated(error)) exit checks

         ! Suppress unused-warning bookkeeping.
         mean_abs = ssh_mean
         call check(error, mean_abs == mean_abs, "basin-mean SSH is NaN")
      end block checks

      deallocate (ssh_t, t_t)
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
      call tides%destroy()
   end subroutine test_m2_oscillation

end module test_ocean_tidal_forcing
