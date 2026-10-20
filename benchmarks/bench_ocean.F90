!! Standalone NetCDF-free ocean throughput benchmark.
!!
!! Drives the production ocean dyn-core on a namelist-configured basin
!! through the SAME lifecycle the driver runs: `engine_setup` (config
!! validation of the process grid, IC/bathymetry seed, every
!! `configure_ocean_*` stage in driver order, halo/periodic/fold seeding,
!! `configure_ocean_bt_split`, land mask, porous, stability audit),
!! `engine_enter_data`, then `engine_step` + `engine_step_ice` +
!! `engine_step_finalize` per outer step, `engine_exit_data`,
!! `engine_teardown`.  Nothing is re-copied here, so the benchmark cannot
!! drift from `driver_run` again: the only differences are that the diag
!! manager is forced off (no NetCDF stream) and only the time loop is timed.
!! It therefore builds + runs with `-DRDB_ENABLE_NETCDF=OFF` (Intel / AMD
!! portability, profiling).
!!
!! Output is console only: a config echo, a periodic progress line with a
!! NaN/CFL guard, an in-terminal ASCII map of the free surface (so you can
!! see the gyres without NetCDF), a Mcells/s throughput number, and the
!! per-kernel profiler breakdown.
!!
!! Usage:
!!   ./build/benchmarks/bench_ocean [namelist] [max_steps]
!!   (default: validation_examples/ocean/double_gyre/double_gyre_mom6.nml;
!!    `max_steps` > 0 truncates the run to that many outer steps — the CTest
!!    smoke guard uses it.)
program bench_ocean
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use pic_timer, only: timer_type
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_dyn, only: ocean_dyn_flush_tracer_window
   use rdb_profiler, only: profiler_init, profiler_report
   implicit none

   character(len=*), parameter :: DEFAULT_NML = &
                                  "validation_examples/ocean/double_gyre/double_gyre_mom6.nml"

   type(config_t) :: cfg
   type(ocean_engine_t) :: engine
   type(timer_type) :: wall_timer

   character(len=512) :: nml_path
   character(len=32) :: arg2
   integer :: nstep, n_inner, ng, max_steps, ios, ierr
   real(wp) :: t_current, dt, wall_sec, mcells_per_sec
   real(wp) :: next_status, status_interval

   ! ---- Namelist path (CLI arg #1, else the double-gyre MOM6 reference) ----
   if (command_argument_count() >= 1) then
      call get_command_argument(1, nml_path)
   else
      nml_path = DEFAULT_NML
   end if
   max_steps = 0
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg2)
      read (arg2, *, iostat=ios) max_steps
      if (ios /= 0) max_steps = 0
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

   ! ---- Setup: the driver's own engine path (may auto-set cfg n_inner) ----
   ! Diagnostics are the one deliberate difference: no NetCDF stream.
   cfg%ocean%diag%enabled = .false.
   call engine_setup(engine, cfg, ierr)
   if (ierr /= OCEAN_STATUS_OK) then
      write (error_unit, "(a)") "bench_ocean: engine_setup failed"
      stop 1
   end if

   n_inner = engine%n_inner
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

   ! ---- Device placement (the driver's engine_enter_data) ----
   call engine_enter_data(engine, cfg)

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
      call engine_step(engine, dt, t_current, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) stop 1
      call engine_step_ice(engine, cfg, dt, t_current, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) stop 1
      call engine_step_finalize(engine, dt, t_current, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) stop 1
      t_current = t_current + dt
      nstep = nstep + 1
      if (max_steps > 0 .and. nstep >= max_steps) then
         call check_finite_or_die(nstep, t_current)
         call report_progress(nstep, t_current)
         exit
      end if
      if (t_current >= next_status .or. t_current >= cfg%t_end) then
         call check_finite_or_die(nstep, t_current)
         call report_progress(nstep, t_current)
         next_status = next_status + status_interval
      end if
   end do
   call ocean_dyn_flush_tracer_window(engine%grid, engine%state%metrics, &
                                      engine%state%dyn, engine%state%continuity, &
                                      engine%state%multilayer, bc=engine%state%bc)
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

   call engine_exit_data(engine)
   call engine_teardown(engine)

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
      associate (h => engine%state%multilayer%h_layer, &
                 u => engine%state%multilayer%u_face_x_layer, &
                 v => engine%state%multilayer%v_face_y_layer, &
                 eta => engine%state%dyn%bt_work%bt_eta)
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
         else if (cfg%ocean%pgf%maxvel > 0.0_wp .and. &
                  max(maxval(abs(u)), maxval(abs(v))) >= 0.9_wp*cfg%ocean%pgf%maxvel) then
            ! The velocity-truncation clamp keeps a blown run FINITE (phantom
            ! +-maxvel currents), so finiteness alone cannot catch it.
            bad = "|u| at the maxvel clamp (blown run)"
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
      associate (u => engine%state%multilayer%u_face_x_layer, &
                 v => engine%state%multilayer%v_face_y_layer, &
                 eta => engine%state%dyn%bt_work%bt_eta)
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
      associate (eta => engine%state%dyn%bt_work%bt_eta)
         !$acc update self(eta)
      end associate
      ng = cfg%nghost
      nxp = cfg%nx
      nyp = cfg%ny
      emin = huge(1.0_wp)
      emax = -huge(1.0_wp)
      do jj = ng + 1, ng + nyp
         do ii = ng + 1, ng + nxp
            emin = min(emin, engine%state%dyn%bt_work%bt_eta(ii, jj))
            emax = max(emax, engine%state%dyn%bt_work%bt_eta(ii, jj))
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
            val = engine%state%dyn%bt_work%bt_eta(ii, jj)
            frac = (val - emin)/(emax - emin)
            lev = 1 + int(frac*real(len(RAMP) - 1, wp) + 0.5_wp)
            lev = max(1, min(len(RAMP), lev))
            line(c:c) = RAMP(lev:lev)
         end do
         write (output_unit, "(a, a)") "   ", line(1:cols)
      end do
   end subroutine print_eta_map

end program bench_ocean
