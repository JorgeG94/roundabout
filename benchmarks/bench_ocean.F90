!! Standalone NetCDF-free ocean throughput benchmark.
!!
!! Drives the production ocean dyn-core (`ocean_dyn_step_split`) on a
!! namelist-configured basin, reusing the EXACT setup path the driver
!! uses — `ocean_state%init_from_config` + `ocean_state_seed_from_cfg` +
!! `register_default_tracers` + the shared `configure_ocean_*` helpers
!! (`rdb_ocean_setup`) — but skipping the diag-manager / NetCDF I/O.  So
!! it builds + runs with `-DRDB_ENABLE_NETCDF=OFF` (Intel / AMD
!! portability, profiling) while staying bit-faithful to the stress-tested
!! driver config (e.g. double_gyre_mom6.nml).
!!
!! Output is console only: a config echo, a periodic progress line with a
!! NaN/CFL guard, an in-terminal ASCII map of the free surface (so you can
!! see the gyres without NetCDF), a Mcells/s throughput number, and the
!! per-kernel profiler breakdown.
!!
!! Usage:
!!   ./build/benchmarks/bench_ocean [namelist]
!!   (default: validation_examples/ocean/double_gyre/double_gyre_mom6.nml)
program bench_ocean
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use pic_timer, only: timer_type
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t, read_config, validate_config
   use rdb_state, only: register_default_tracers
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, &
                              ocean_state_exit_data, ocean_state_seed_from_cfg
   use rdb_ocean_setup, only: configure_ocean_metrics, configure_ocean_land_mask, &
                              configure_ocean_forcing, configure_ocean_drag, &
                              configure_ocean_vmix, configure_ocean_lateral, &
                              configure_ocean_pgf, configure_ocean_bt, &
                              configure_ocean_bt_split
   use rdb_ocean_vcoord, only: parse_ocean_vcoord_type
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_dyn, only: ocean_dyn_step_split
   use rdb_profiler, only: profiler_init, profiler_report
   implicit none

   character(len=*), parameter :: DEFAULT_NML = &
                                  "validation_examples/ocean/double_gyre/double_gyre_mom6.nml"

   type(config_t) :: cfg
   type(hgrid_t) :: grid
   type(ocean_state_t) :: os
   type(ocean_surface_flux_t) :: sf
   type(timer_type) :: wall_timer

   character(len=512) :: nml_path
   integer :: nstep, n_inner, ng
   real(wp) :: t_current, dt, wall_sec, mcells_per_sec
   real(wp) :: next_status, status_interval

   ! ---- Namelist path (CLI arg #1, else the double-gyre MOM6 reference) ----
   if (command_argument_count() >= 1) then
      call get_command_argument(1, nml_path)
   else
      nml_path = DEFAULT_NML
   end if

   call read_config(trim(nml_path), cfg)
   call validate_config(cfg)   ! aborts on invalid config

   if (trim(cfg%sim_type) /= "ocean") then
      write (error_unit, "(a)") "bench_ocean: namelist sim_type must be 'ocean' (got '"// &
         trim(cfg%sim_type)//"')"
      stop 1
   end if
   if (cfg%dt_fixed <= 0.0_wp) then
      write (error_unit, "(a)") "bench_ocean: requires dt_fixed > 0"
      stop 1
   end if

   ! ---- Setup — the SAME sequence as driver_run_ocean, minus diag/IO ----
   call grid%init(cfg%nx, cfg%ny, cfg%nghost, cfg%dx, cfg%dy)
   call os%init_from_config(cfg, grid)

   ! vcoord params BEFORE the seed (the seed's build_zref_full reads them).
   if (os%use_multilayer) then
      os%vcoord%coord_type = parse_ocean_vcoord_type(cfg%vcoord_type)
      os%vcoord%zstar_h_surf_target = cfg%zstar_h_surf_target
      os%vcoord%zstar_h_min = cfg%zstar_h_min
   end if

   call ocean_state_seed_from_cfg(os, grid, cfg)
   call register_default_tracers( &
      os%multilayer%tracers(os%multilayer%idx_salinity), &
      os%multilayer%tracers(os%multilayer%idx_temperature), cfg)

   call sf%init(grid)
   call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

   call configure_ocean_metrics(cfg, os, grid, 0)
   call configure_ocean_forcing(cfg, os, grid, 0)
   call configure_ocean_drag(cfg, os, 0)
   call configure_ocean_vmix(cfg, os, 0)
   call configure_ocean_lateral(cfg, os, grid, 0)
   call configure_ocean_pgf(cfg, os, 0)
   call configure_ocean_bt(cfg, os, grid, 0)
   call configure_ocean_bt_split(cfg, os, grid, 0)   ! may auto-set cfg%ocean%bt%n_inner
   call configure_ocean_land_mask(cfg, os, grid, 0)  ! static land masking (CHUNK A)

   n_inner = cfg%ocean%bt%n_inner
   if (n_inner < 1) then
      write (error_unit, "(a)") "bench_ocean: needs split-explicit barotropic "// &
         "(ocean_bt n_inner >= 1 or auto_n_inner = .true.)"
      stop 1
   end if

   call print_banner()
   write (output_unit, '(a, i6, " x ", i6, " x ", i4, "  (dx = ", f9.1, " m)")') &
      " grid       :", cfg%nx, cfg%ny, cfg%nz_layers, cfg%dx
   write (output_unit, '(a, f10.1, " s  dt = ", f8.1, " s   n_inner = ", i4)') &
      " duration   :", cfg%t_end, cfg%dt_fixed, n_inner
   write (output_unit, "(a)") " namelist   : "//trim(nml_path)
   write (output_unit, "(a)") ""

   ! ---- Device placement (same as the driver) ----
   call ocean_state_enter_data(os)
   !$acc enter data copyin(sf)
   call sf%enter_data()

   call profiler_init(.true.)

   ! ---- Time loop ----
   t_current = 0.0_wp
   nstep = 0
   status_interval = cfg%status_interval
   if (status_interval <= 0.0_wp) status_interval = max(cfg%t_end/20.0_wp, cfg%dt_fixed)
   next_status = status_interval

   call wall_timer%start()
   do while (t_current < cfg%t_end)
      dt = min(cfg%dt_fixed, cfg%t_end - t_current)
      call ocean_dyn_step_split(grid, os%metrics, os%dyn, os%eos, &
                                os%coriolis_adv, os%continuity, &
                                os%pressure_force, os%hvisc, &
                                os%bdrag, os%surface_stress, &
                                os%vert_advect, os%hdiff_tracer, &
                                os%vdiff, os%vmix, &
                                os%multilayer, dt, n_inner, sf=sf, &
                                vcoord=os%vcoord, t=t_current, &
                                lateral_mix=os%lateral_mix)
      t_current = t_current + dt
      nstep = nstep + 1
      if (t_current >= next_status .or. t_current >= cfg%t_end) then
         call check_finite_or_die(nstep, t_current)
         call report_progress(nstep, t_current)
         next_status = next_status + status_interval
      end if
   end do
   call wall_timer%stop()

   wall_sec = wall_timer%get_elapsed_time()
   mcells_per_sec = real(nstep, wp)*real(cfg%nx, wp)*real(cfg%ny, wp)* &
                    real(cfg%nz_layers, wp)/(wall_sec*1.0e6_wp)

   call print_eta_map()
   write (output_unit, "(a)") ""
   write (output_unit, '(a, i8, " steps in ", f9.2, " s   (", f9.3, " Mcells·layer/s)")') &
      " throughput :", nstep, wall_sec, mcells_per_sec
   write (output_unit, "(a)") ""
   call profiler_report()

   call sf%exit_data()
   !$acc exit data delete(sf)
   call ocean_state_exit_data(os)

contains

   subroutine print_banner()
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "  ╔══════════════════════════════════════════════════╗"
      write (output_unit, "(a)") "  ║   Roundabout · ocean throughput benchmark            ║"
      write (output_unit, "(a)") "  ║   C-grid split-RK2 dyn-core · NetCDF-free        ║"
      write (output_unit, "(a)") "  ╚══════════════════════════════════════════════════╝"
   end subroutine print_banner

   subroutine check_finite_or_die(step_idx, t_secs)
      !! Stop on the first non-finite prognostic field so a blown run
      !! doesn't keep burning wall time producing NaNs.
      integer, intent(in) :: step_idx
      real(wp), intent(in) :: t_secs
      character(len=64) :: bad
      ! associate-leaf: ifx/flang OpenMP-target silently no-op a `update`
      ! (target update from) whose operand is a deep DT chain, leaving the
      ! host copy stale on Intel/AMD (NVHPC tolerates the chain).  Hoist each
      ! to a one-deep name so the D->H sync lands.  Same rule as
      ! solver_enter_data's enter/exit-data clauses.
      associate (h => os%multilayer%h_layer, &
                 u => os%multilayer%u_face_x_layer, &
                 v => os%multilayer%v_face_y_layer, &
                 eta => os%dyn%bt_work%bt_eta)
         !$acc update self(h, u, v, eta)
         bad = ""
         if (.not. all(ieee_is_finite(h))) then
            bad = "h_layer"
         else if (.not. all(ieee_is_finite(u))) then
            bad = "u_face_x_layer"
         else if (.not. all(ieee_is_finite(v))) then
            bad = "v_face_y_layer"
         else if (.not. all(ieee_is_finite(eta))) then
            bad = "bt_eta"
         end if
      end associate
      if (len_trim(bad) > 0) then
         write (error_unit, "(a)") ""
         write (error_unit, "(a)") "*** bench_ocean: non-finite state ("//trim(bad)//") ***"
         write (error_unit, '(a, i8, "   t = ", f12.1, " s")') "   step ", step_idx, t_secs
         stop 1
      end if
   end subroutine check_finite_or_die

   subroutine report_progress(step_idx, t_secs)
      integer, intent(in) :: step_idx
      real(wp), intent(in) :: t_secs
      real(wp) :: eta_mx, u_mx, elapsed, eta_est
      ! associate-leaf so the device->host sync lands on ifx/flang (deep DT
      ! chains in `update` are a silent no-op there; see check_finite_or_die).
      associate (u => os%multilayer%u_face_x_layer, &
                 v => os%multilayer%v_face_y_layer, &
                 eta => os%dyn%bt_work%bt_eta)
         !$acc update self(u, v, eta)
         eta_mx = maxval(abs(eta))
         u_mx = max(maxval(abs(u)), maxval(abs(v)))
      end associate
      elapsed = wall_timer%get_elapsed_time()
      eta_est = 0.0_wp
      if (t_secs > 0.0_wp) eta_est = elapsed*(cfg%t_end - t_secs)/t_secs
      write (output_unit, &
         '(a, i7, "  t=", f10.1, "s  |eta|=", f6.3, "m  |u|=", f6.3, &
           &"m/s  wall=", f7.1, "s  ETA=", f7.1, "s")') &
         " step ", step_idx, t_secs, eta_mx, u_mx, elapsed, eta_est
   end subroutine report_progress

   subroutine print_eta_map()
      !! Coarse ASCII heat-map of the free surface η — the gyres in-terminal.
      integer, parameter :: MAP_COLS = 70, MAP_ROWS = 24
      character(len=*), parameter :: RAMP = " .:-=+*o#%@"
      integer :: cols, rows, r, c, ii, jj, lev, nxp, nyp
      real(wp) :: emin, emax, val, frac
      character(len=MAP_COLS) :: line
      ! associate-leaf so the D->H sync lands on ifx/flang (see report_progress).
      associate (eta => os%dyn%bt_work%bt_eta)
         !$acc update self(eta)
      end associate
      ng = cfg%nghost
      nxp = cfg%nx
      nyp = cfg%ny
      emin = huge(1.0_wp)
      emax = -huge(1.0_wp)
      do jj = ng + 1, ng + nyp
         do ii = ng + 1, ng + nxp
            emin = min(emin, os%dyn%bt_work%bt_eta(ii, jj))
            emax = max(emax, os%dyn%bt_work%bt_eta(ii, jj))
         end do
      end do
      write (output_unit, "(a)") ""
      write (output_unit, '(a, f8.3, " m  to  ", f8.3, " m   (north up)")') &
         " free surface η:  ", emin, emax
      if (emax - emin <= 0.0_wp) then
         write (output_unit, "(a)") "   (flat)"
         return
      end if
      cols = min(MAP_COLS, nxp)
      rows = min(MAP_ROWS, nyp)
      do r = rows, 1, -1
         jj = ng + 1 + int((real(r, wp) - 0.5_wp)/real(rows, wp)*real(nyp, wp))
         do c = 1, cols
            ii = ng + 1 + int((real(c, wp) - 0.5_wp)/real(cols, wp)*real(nxp, wp))
            val = os%dyn%bt_work%bt_eta(ii, jj)
            frac = (val - emin)/(emax - emin)
            lev = 1 + int(frac*real(len(RAMP) - 1, wp) + 0.5_wp)
            lev = max(1, min(len(RAMP), lev))
            line(c:c) = RAMP(lev:lev)
         end do
         write (output_unit, "(a, a)") "   ", line(1:cols)
      end do
   end subroutine print_eta_map

end program bench_ocean
