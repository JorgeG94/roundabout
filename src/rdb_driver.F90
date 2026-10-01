!! Driver: full compute-rank lifecycle for the Roundabout solver.
module rdb_driver
   !! Owns the end-to-end compute-rank path that used to live inline in
   !! `app/main.F90`: setup (decomp / state allocation / scatter / BC /
   !! forcing / nesting / output / GPU enter_data), the time-stepping
   !! loop (forcing update, solver_step, status print, snapshot output,
   !! restart, gauges), and finalize (exit_data, final snapshot,
   !! diagnostics, profiler report, resource cleanup).
   !!
   !! Exposed as `driver_run(cfg)` so:
   !!   * `app/main.F90` stays thin (MPI init, argv parsing, dispatch).
   !!   * Unit tests can drive a full compute lifecycle from an in-code
   !!     `config_t`, catching integration-level bugs (OpenACC partial
   !!     presence, init ordering, missing attaches) that per-kernel
   !!     unit tests miss.
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_restart_read
   use rdb_ocean_boundary_types, only: ocean_bc_has_tracer_open_edge
   use rdb_ocean_dyn, only: ocean_dyn_flush_tracer_window, SPLIT_SCHEME_PRED_CORR
   ! P2.4: the setup/step/teardown sequence itself (21 configure_ocean_*-family
   ! stages, device placement, the dyn-core advance) lives in ONE shared engine
   ! module now, so driver_run_ocean, rdb_ocean_api and bench_ocean stop
   ! carrying independent copies. See src/core/ocean/rdb_ocean_engine.F90.
   ! P2.4b: the sea-ice per-step block (frazil/EVP/thermo-driver/transport)
   ! moved into the same shared engine (`engine_step_ice`) — the driver no
   ! longer `use`s the individual rdb_ice_* kernel modules directly.
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown, diag_vgrid_from_name
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_console_stats, only: console_stats_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, &
                                      ocean_budget_stage_weight
   use rdb_io_netcdf, only: output_rank_filename
   use rdb_halo, only: halo_allreduce_sum
   use pic_logger, only: logger => global_logger, &
                         debug_level, verbose_level, info_level, &
                         performance_level, warning_level, error_level
   use pic_strings, only: to_string
   use pic_io, only: to_char_count
   use pic_types, only: int64
   use pic_timer, only: timer_type
   use rdb_comm_env, only: comm_env_rank, comm_env_size, comm_env_bcast_real, &
                           comm_env_abort, &
                           comm_env_io_server_rank, comm_env_compute_rank, &
                           comm_env_compute_size
   use rdb_profiler, only: profiler_init, profiler_start, profiler_stop, &
                           profiler_report, profiler_end
   use rdb_banner, only: print_banner
   use rdb_mem_report, only: mem_host_rss_bytes, mem_log_state_budget, &
                             mem_log_device_actuals, mem_log_device_growth, &
                             mem_set_counted_budget, mem_log_computed_budget, &
                             mem_log_computed_line
   use pic_knowledge, only: get_knowledge
   implicit none
   private

   public :: driver_run
   public :: configure_log_level
   public :: diag_vgrid_from_name

contains

   subroutine driver_run(cfg)
      !! Compute-rank entry point.  Runs the Arakawa C-grid +
      !! continuity-PPM split-RK2 ocean dynamical core in `src/core/ocean/`.
      type(config_t), intent(inout) :: cfg

      call driver_run_ocean(cfg)
   end subroutine driver_run

   subroutine report_throughput(global_cells, n_steps, elapsed, compute_size, &
                                compute_rank, total_mcells)
      !! Log per-GPU and total horizontal throughput on the root compute rank.
      !!
      !! `global_cells` is the total horizontal cells across all ranks -- the
      !! caller does any reduction needed to obtain it (the ocean path
      !! already knows the global grid). Per-GPU is the balanced-partition
      !! average
      !! (`total / compute_size`). Throughput is horizontal cells x steps, NOT
      !! x nz (the codebase convention). Optionally returns the total so a
      !! caller can thread it into a later summary.
      real(wp), intent(in) :: global_cells
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: elapsed
      integer, intent(in) :: compute_size
      integer, intent(in) :: compute_rank
      real(wp), intent(out), optional :: total_mcells

      real(wp) :: per_gpu, total

      total = 0.0_wp
      per_gpu = 0.0_wp
      if (elapsed > 0.0_wp .and. compute_size > 0) then
         total = global_cells*real(n_steps, wp)/elapsed/1.0e6_wp
         per_gpu = total/real(compute_size, wp)
      end if
      if (present(total_mcells)) total_mcells = total
      if (compute_rank == 0) then
         call logger%info("  throughput: "//to_string(per_gpu)//" Mcells/s per GPU, "// &
                          to_string(total)//" Mcells/s total ("// &
                          to_string(compute_size)//" rank(s))")
      end if
   end subroutine report_throughput

   subroutine driver_run_ocean(cfg)
      !! Run the full compute-rank lifecycle for the ocean regime.
      !!
      !! Phase 6 scope (single-rank, serial NetCDF, scalar surface
      !! forcing, stub IC from cfg scalars).  Sets up the ocean god
      !! state, registers the default tracer + diag-manager variables,
      !! opens a per-rank NetCDF output stream, runs `ocean_dyn_step`
      !! per outer step at `cfg%dt_fixed`, calls `diag%step` after each
      !! advance, and closes the stream at the end.  MPI scatter, ocean
      !! halo exchanges, restart, gauges, and the I/O server hand-off
      !! are deferred to Phase 7+.
      type(config_t), intent(inout) :: cfg

      type(ocean_engine_t) :: engine
         !! P2.4: owns the god state + the setup products that used to
         !! be driver stack locals (grid, geo, bc_source, decomp,
         !! evp_params, ic_par, cfl_vtol) — see
         !! `src/core/ocean/rdb_ocean_engine.F90`. Bound to the local
         !! names `ocean_state`/`grid`/`geo`/`decomp`/`evp_params`/
         !! `cfl_vtol` via `associate` below so the (unchanged) time-loop
         !! + teardown body reads exactly as it did before this refactor.
      type(timer_type) :: wall_timer, setup_timer, t_enter_data
      type(console_stats_t) :: console_stats
         !! MOM6-style status-line scalars (mass / KE / mean S, T /
         !! max-CFL) emitted at the driver's `status_interval` cadence.

      real(wp) :: t_current, dt, wall_elapsed
      real(wp) :: next_status_time
      real(wp) :: next_restart_time, t_restart
      integer :: n_steps, step_restart
      logical :: use_restart_output
      character(len=512) :: restart_filename
      integer :: mpi_rank, compute_rank, compute_size
      integer :: ntrunc_last_report
      integer(int64) :: n_limited_last_report
      integer(int64) :: rss0, rss1
      integer :: setup_ierr, step_ierr
         !! `rdb_ocean_status` codes from `engine_setup`/`engine_step`/
         !! `engine_step_ice`/`engine_step_finalize` — the driver still
         !! aborts on failure (via `comm_env_abort`, which is MPI-safe,
         !! rather than the bare `error stop` the engine falls back to
         !! when `ierr` is absent), it just does so through the
         !! returnable-status path
         !! now instead of engine-internal `error stop`s.

      mpi_rank = comm_env_rank()
      compute_rank = comm_env_compute_rank()
      compute_size = comm_env_compute_size()

      if (cfg%use_io_server) then
         if (compute_rank == 0) then
            call logger%warning("I/O server is not supported on the ocean path; "// &
                                "diag manager writes serial NetCDF per rank.")
         end if
      end if

      call setup_timer%start()

      ! Host RSS before state allocation — host growth across init is an
      ! ADVISORY host-side figure + a conservative upper bound for the
      ! pre-flight OOM warning; the true device footprint is measured
      ! across enter_data (rdb_mem_report).
      rss0 = mem_host_rss_bytes()

      ! P2.4: the full 21-stage configure_ocean_* sequence (decomp -> grid
      ! -> god state -> IC seed -> restart -> diag -> surface-flux seed ->
      ! metrics/tides/psurf/forcing/dataovr/drag/hdiff/vmix/tracers/
      ! lateral -> ice EVP params -> pgf/bt/bt_split/bc/sponge -> ghost
      ! wraps -> halo init -> land mask -> wave drag/porous -> ice IC),
      ! in the driver's exact original order, now lives in ONE place —
      ! see `engine_setup`'s docstring for why order is load-bearing here.
      call engine_setup(engine, cfg, setup_ierr, compute_rank=compute_rank, &
                        compute_size=compute_size, mpi_rank=mpi_rank, &
                        restart_file=cfg%restart_file, t_restart=t_restart, &
                        step_restart=step_restart)
      if (setup_ierr /= OCEAN_STATUS_OK) call comm_env_abort(1)
      ! engine_setup itself logs the "Ocean warm restart from ..." info line
      ! (gated on restart_file being non-blank) — nothing to do here.

      ! Place ocean state on the GPU.  `sf` is a scalar-only derived
      ! type (no allocatables), so a plain copyin is enough.  First
      ! OpenACC directive in the process — implicitly does CUDA context
      ! init + CUBIN load on top of the H→D memcpy.
      ! Memory budget: log the host-growth estimate + device free/total
      ! BEFORE the mapping (the motivating 10M-cell failure died here).
      rss1 = mem_host_rss_bytes()
      if (compute_rank == 0 .and. rss0 >= 0_int64 .and. rss1 >= 0_int64) then
         call mem_log_state_budget("ocean state", max(rss1 - rss0, 0_int64))
      end if

      ! Initialise profiler BEFORE engine_enter_data so that per-slot
      ! NVTX ranges inside the orchestrator land in the same profiler context.
      ! The old placement (post-banner) dropped every NVTX range in enter_data.
      call profiler_init()

      call t_enter_data%start()
      call profiler_start("ocean_enter_data")
      call engine_enter_data(engine, cfg)
      call profiler_stop("ocean_enter_data")
      call t_enter_data%stop()

      ! Latch the counted state-array footprint for the enter_data
      ! reconciliation (mem_log_device_actuals reads it); the human-readable
      ! budget prints in the post-setup banner block below.
      if (compute_rank == 0) call mem_set_counted_budget(engine%state%bytes())

      if (compute_rank == 0) call mem_log_device_actuals("ocean state mapped")

      call setup_timer%stop()

      if (compute_rank == 0) then
         call print_banner("ocean")
         call logger%info("")
         call logger%info("  Regime:   ocean (C-grid + continuity-PPM + split-RK2)")
         call logger%info("  Grid:     "//to_string(cfg%nx)//" x "//to_string(cfg%ny)// &
                          " ("//to_char_count(int(cfg%nx, int64)*int(cfg%ny, int64))//")")
         call logger%info("  Spacing:  dx = "//to_string(cfg%dx)//" m, dy = "// &
                          to_string(cfg%dy)//" m")
         call logger%info("  Layers:   "//to_string(cfg%nz_layers))
         call logger%info("  Duration: "//to_string(cfg%t_end)//" s")
         call logger%info("  Fixed dt: "//to_string(cfg%dt_fixed)//" s")
         if (cfg%ocean%diag%enabled) then
            call logger%info("  Diag:     "//trim(cfg%ocean%diag%filename)// &
                             " every "//to_string(cfg%ocean%diag%dt_out)//" s")
         end if
         ! Counted state-array footprint (exact allocatable sum; the gated
         ! default-off closures count 0, so the breakdown reflects the
         ! conditional allocation).  Latched earlier for the reconciliation.
         block
            integer(int64) :: b_total, b_baro, b_layers
            b_total = engine%state%bytes()
            b_baro = engine%state%barotropic%bytes()
            b_layers = 0_int64
            if (engine%state%use_multilayer) b_layers = engine%state%multilayer%bytes()
            call mem_log_computed_budget("ocean state", b_total)
            call mem_log_computed_line("barotropic (C-grid)", b_baro)
            call mem_log_computed_line("layers + tracers", b_layers)
            call mem_log_computed_line("metrics + dyn-core + closures", &
                                       max(b_total - b_baro - b_layers, 0_int64))
         end block
         call logger%info("  Setup time: "// &
                          to_string(setup_timer%get_elapsed_time())//" s")
         ! P2.4: the setup phase is now one engine_setup call (was
         ! alloc/ic_seed/diag_reg/forcing_cfg sub-timers around distinct
         ! pieces of what is now a single sequence) — enter_data stays
         ! separately timed since it is genuinely a separate engine call.
         call logger%info("    enter_data  : "// &
                          to_string(t_enter_data%get_elapsed_time())//" s")
         call logger%info("")
         call logger%info("     Step         t (s)         dt (s)     Wall (s)  Remaining (s)")
         call logger%info("  -------  ------------  ------------  -----------  -------------")
      end if

      t_current = t_restart
      n_steps = step_restart
      ntrunc_last_report = 0
      n_limited_last_report = 0
      if (cfg%status_interval > 0.0_wp) then
         next_status_time = t_current + cfg%status_interval
      end if

      ! Restart-write cadence (regime-agnostic knob, mirrors coastal).
      use_restart_output = (cfg%restart_interval > 0.0_wp)
      if (use_restart_output) then
         next_restart_time = t_current + cfg%restart_interval
      end if

      ! P2.4: the rest of this subroutine (loop-prep snapshot, time loop,
      ! teardown) is UNCHANGED from before the engine extraction — it reads
      ! `ocean_state`/`grid`/`geo`/`decomp`/`evp_params`/`cfl_vtol` exactly
      ! as it did when those were driver stack locals; they are now
      ! `associate`-bound to the matching `engine` components so no call
      ! site below had to change.
      associate (ocean_state => engine%state, grid => engine%grid, geo => engine%geo, &
                 decomp => engine%decomp, evp_params => engine%evp_params, &
                 cfl_vtol => engine%cfl_vtol)

         ! Capture conservation-reference snapshot — the first call latches
         ! `console_stats%mass0 / salt0 / heat0 / ke0`, so every subsequent
         ! fire reports drift relative to this reference rather than
         ! post-first-step.  On a cold start this is the true t=0 IC; on a
         ! warm restart it is the resumed state at `t_restart` (drift is then
         ! measured from the restart point, not the original t=0 — the
         ! pre-restart baseline lives in the run that wrote the checkpoint).
         ! Skip the Salt/Heat lines when
         ! thermodynamics is disabled — those tracers stay at IC by
         ! construction and the columns are noise.
         console_stats%report_thermodynamics = cfg%ocean%thermo%enable_thermodynamics
         ! COLLECTIVE: called on ALL ranks (allreduce inside); is_root gates the
         ! printed lines to rank 0. compute_rank is keyword-passed (it follows the
         ! optional budget args in the signature).
         ! `heat_budget_frazil` rides the F2008 unallocated-actual-to-optional
         ! rule: with ice off the slot array is never allocated, the dummy is
         ! absent, and the report is bit-identical to the pre-ice path.
         !
         ! `horiz_adv_budget_valid` is deliberately NOT passed (absent ⇒
         ! `.true.`, the closed budget stays on).  It used to carry
         ! `dt_tracer_advect_ratio == 1`, because the windowed tracer-advect
         ! path left `*_budget_horiz_adv` unfilled — which demoted the whole
         ! Salt/Heat `Error` column to raw drift, and raw drift cannot
         ! subtract a surface source, so a conservative run with `q_heat`
         ! printed the flux itself as a ~5e-5 "leak".  The windowed path is
         ! instrumented end to end now (`continuity_tracer_drain` + both
         ! halves of the concentration hold); see
         ! `ocean_budget_is_active`'s docstring and the regression test
         ! `windowed_drain_with_surface_flux`.
         if (cfl_vtol > 0.0_wp) then
            call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                            ocean_state%multilayer, &
                                            t_current, cfg%dt_fixed, n_steps, &
                                            redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                 ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                            cfl_vanish_tol=cfl_vtol, &
                                            heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                            ice_part_size=ocean_state%ice%part_size, &
                                            ice_m_ice=ocean_state%ice%m_ice, &
                                            ice_ncat=ocean_state%ice%ncat, &
                                            compute_rank=compute_rank, &
                                            budget_stage_weight= &
                                            ocean_budget_stage_weight( &
                                            ocean_state%dyn%split_scheme == &
                                            SPLIT_SCHEME_PRED_CORR), &
                                            reproducing_sums=cfg%ocean%diag%reproducing_sums)
         else
            call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                            ocean_state%multilayer, &
                                            t_current, cfg%dt_fixed, n_steps, &
                                            redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                 ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                            heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                            ice_part_size=ocean_state%ice%part_size, &
                                            ice_m_ice=ocean_state%ice%m_ice, &
                                            ice_ncat=ocean_state%ice%ncat, &
                                            compute_rank=compute_rank, &
                                            budget_stage_weight= &
                                            ocean_budget_stage_weight( &
                                            ocean_state%dyn%split_scheme == &
                                            SPLIT_SCHEME_PRED_CORR), &
                                            reproducing_sums=cfg%ocean%diag%reproducing_sums)
         end if

         if (compute_rank == 0) call wall_timer%start()
         call profiler_start("time_loop")

         do while (t_current < cfg%t_end)
            dt = min(cfg%dt_fixed, cfg%t_end - t_current)

            ! P2.4: file-forcing update/apply, boundary-data refresh,
            ! porous-area refresh, then the dyn-core advance (split-RK2 or
            ! the legacy unsplit path, exactly as `engine%n_inner` dictates)
            ! — all now one `engine_step` call; see its docstring for why
            ! `ocean_surface_flux_assemble`/the diag step are a SEPARATE
            ! `engine_step_finalize` call made after the sea-ice block below
            ! rather than folded in here.
            call engine_step(engine, dt, t_current, ierr=step_ierr)
            if (step_ierr /= OCEAN_STATUS_OK) call comm_env_abort(1)

            ! P2.4b: sea-ice per-step physics (frazil accumulation, EVP
            ! dynamics + ice->ocean stress, and, at thermo cadence, transport
            ! + the atmospheric-forcing/basal/frazil-uptake/column-thermo/
            ! snowfall/brine/heat/shortwave/ITD chain) now lives on the same
            ! shared engine as the dyn-core advance -- see
            ! `engine_step_ice`'s docstring (the "MANDATED ORDER" contract
            ! moved there verbatim, unchanged). No-op when
            ! `&ocean_ice_nml enable = .false.`.
            call engine_step_ice(engine, cfg, dt, t_current, ierr=step_ierr)
            if (step_ierr /= OCEAN_STATUS_OK) call comm_env_abort(1)

            ! PR-12: derive Q_heat/Q_salt from the component set (no-op
            ! unless &ocean_forcing_nml enable_components), then (when
            ! configured) one diag-manager step.  MUST sit here —
            ! immediately after the ice block closes, before t_current
            ! advances — to preserve the ice coupler's documented one-step
            ! lag (the ice writes its components at the end of outer step
            ! N; the ocean integrates the assembled net field on step N+1,
            ! same convention as the pre-PR-12 Q_heat/Q_salt full-overwrite).
            ! This is why `engine_step_finalize` is a call SEPARATE from
            ! `engine_step` above — see its docstring.
            call engine_step_finalize(engine, dt, t_current, ierr=step_ierr)
            if (step_ierr /= OCEAN_STATUS_OK) call comm_env_abort(1)

            t_current = t_current + dt
            n_steps = n_steps + 1

            ! Per-rank restart checkpoint at cadence.  Fires at the top of a
            ! completed outer step (step-aligned — RK saves + bt
            ! accumulators are scratch by then).  D->H pulls happen inside.
            if (use_restart_output .and. t_current >= next_restart_time) then
               ! FIX-1 (spec §(c) mandatory flush): a restart checkpoint is the
               ! prognostic snapshot a run resumes from — hTr and h_layer are
               ! both read back as state.  If this write lands mid-accumulation-
               ! window (run cadence not an exact multiple of
               ! dt_tracer_advect_ratio), drain first so frozen hTr is never
               ! checkpointed against advanced h_layer.  No-op at ratio = 1 and
               ! idempotent (window already empty ⇒ exact no-op).
               call ocean_dyn_flush_tracer_window(grid, ocean_state%metrics, &
                                                  ocean_state%dyn, ocean_state%continuity, &
                                                  ocean_state%multilayer, bc=ocean_state%bc)
               restart_filename = output_rank_filename(trim(cfg%output_dir), "restart", mpi_rank)
               call ocean_state_restart_write(ocean_state, grid, decomp, &
                                              trim(restart_filename), t_current, n_steps)
               next_restart_time = next_restart_time + cfg%restart_interval
            end if

            if (cfg%status_interval > 0.0_wp .and. t_current >= next_status_time) then
               ! COLLECTIVE: ocean_console_stats_report allreduces inside, so it
               ! runs on ALL ranks; is_root (compute_rank) gates its printing.
               ! print_status_line / ntrunc / mem_log stay rank-0-only below.
               ! compute_rank is keyword-passed (it follows the optional args).
               if (cfl_vtol > 0.0_wp) then
                  call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                                  ocean_state%multilayer, &
                                                  t_current, dt, n_steps, &
                                                  redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                       ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                                  cfl_vanish_tol=cfl_vtol, &
                                                  heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                                  ice_part_size=ocean_state%ice%part_size, &
                                                  ice_m_ice=ocean_state%ice%m_ice, &
                                                  ice_ncat=ocean_state%ice%ncat, &
                                                  compute_rank=compute_rank, &
                                                  budget_stage_weight= &
                                                  ocean_budget_stage_weight( &
                                                  ocean_state%dyn%split_scheme == &
                                                  SPLIT_SCHEME_PRED_CORR), &
                                                  reproducing_sums=cfg%ocean%diag%reproducing_sums)
               else
                  call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                                  ocean_state%multilayer, &
                                                  t_current, dt, n_steps, &
                                                  redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                       ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                                  heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                                  ice_part_size=ocean_state%ice%part_size, &
                                                  ice_m_ice=ocean_state%ice%m_ice, &
                                                  ice_ncat=ocean_state%ice%ncat, &
                                                  compute_rank=compute_rank, &
                                                  budget_stage_weight= &
                                                  ocean_budget_stage_weight( &
                                                  ocean_state%dyn%split_scheme == &
                                                  SPLIT_SCHEME_PRED_CORR), &
                                                  reproducing_sums=cfg%ocean%diag%reproducing_sums)
               end if
               ! ntrunc: collective sum so the printed total is global.
               block
                  real(wp) :: nt_l, nt_g
                  nt_l = real(ocean_state%dyn%ntrunc_total, wp)
                  call halo_allreduce_sum(nt_l, nt_g)
                  if (compute_rank == 0) then
                     call print_status_line(n_steps, t_current, dt, &
                                            real(wall_timer%get_elapsed_time(), wp), &
                                            cfg%t_end)
                     if (nint(nt_g) > ntrunc_last_report) then
                        call logger%info("CFL truncations:  "// &
                                         to_string(nint(nt_g) - ntrunc_last_report)// &
                                         " this report ("//to_string(nint(nt_g))// &
                                         " total)")
                        ntrunc_last_report = nint(nt_g)
                     end if
                     ! Quiet drift check: one-shot warning if device memory grew
                     ! > 15% past the post-enter_data baseline (lazy workspaces /
                     ! leak) — no info lines at status cadence.
                     call mem_log_device_growth("ocean status", quiet=.true.)
                  end if
               end block
               ! positive-definite limiter: collective sum, per-report delta log
               ! (mirrors the ntrunc drain above).  The running total grows only
               ! when the knob is on, so the `> last_report` gate alone suffices.
               block
                  real(wp) :: nl_l, nl_g
                  nl_l = real(ocean_state%continuity%n_limited_total, wp)
                  call halo_allreduce_sum(nl_l, nl_g)
                  if (compute_rank == 0) then
                     if (nint(nl_g, int64) > n_limited_last_report) then
                        call logger%info("positive-definite limiter:  "// &
                                         to_string(nint(nl_g, int64) - n_limited_last_report)// &
                                         " faces limited this report ("//to_string(nint(nl_g, int64))// &
                                         " total)")
                        n_limited_last_report = nint(nl_g, int64)
                     end if
                  end if
               end block
               next_status_time = next_status_time + cfg%status_interval
            else if (cfg%status_interval <= 0.0_wp .and. mod(n_steps, 100) == 0) then
               ! COLLECTIVE: same allreduce-inside pattern; the cadence condition is
               ! replicated across ranks (t_current, n_steps identical) so all ranks
               ! enter together.
               if (cfl_vtol > 0.0_wp) then
                  call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                                  ocean_state%multilayer, &
                                                  t_current, dt, n_steps, &
                                                  redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                       ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                                  cfl_vanish_tol=cfl_vtol, &
                                                  heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                                  ice_part_size=ocean_state%ice%part_size, &
                                                  ice_m_ice=ocean_state%ice%m_ice, &
                                                  ice_ncat=ocean_state%ice%ncat, &
                                                  compute_rank=compute_rank, &
                                                  budget_stage_weight= &
                                                  ocean_budget_stage_weight( &
                                                  ocean_state%dyn%split_scheme == &
                                                  SPLIT_SCHEME_PRED_CORR), &
                                                  reproducing_sums=cfg%ocean%diag%reproducing_sums)
               else
                  call ocean_console_stats_report(console_stats, grid, ocean_state%metrics, &
                                                  ocean_state%multilayer, &
                                                  t_current, dt, n_steps, &
                                                  redi_with_open_edge=(ocean_state%redi%enable .and. &
                                                                       ocean_bc_has_tracer_open_edge(ocean_state%bc)), &
                                                  heat_budget_frazil=ocean_state%ice%heat_budget_frazil, &
                                                  ice_part_size=ocean_state%ice%part_size, &
                                                  ice_m_ice=ocean_state%ice%m_ice, &
                                                  ice_ncat=ocean_state%ice%ncat, &
                                                  compute_rank=compute_rank, &
                                                  budget_stage_weight= &
                                                  ocean_budget_stage_weight( &
                                                  ocean_state%dyn%split_scheme == &
                                                  SPLIT_SCHEME_PRED_CORR), &
                                                  reproducing_sums=cfg%ocean%diag%reproducing_sums)
               end if
               ! ntrunc: collective sum.
               block
                  real(wp) :: nt_l, nt_g
                  nt_l = real(ocean_state%dyn%ntrunc_total, wp)
                  call halo_allreduce_sum(nt_l, nt_g)
                  if (compute_rank == 0) then
                     call print_status_line(n_steps, t_current, dt, &
                                            real(wall_timer%get_elapsed_time(), wp), cfg%t_end)
                     if (nint(nt_g) > ntrunc_last_report) then
                        call logger%info("CFL truncations:  "// &
                                         to_string(nint(nt_g) - ntrunc_last_report)// &
                                         " this report ("//to_string(nint(nt_g))// &
                                         " total)")
                        ntrunc_last_report = nint(nt_g)
                     end if
                  end if
               end block
               ! positive-definite limiter: collective sum, per-report delta log
               ! (mirrors the ntrunc drain above).
               block
                  real(wp) :: nl_l, nl_g
                  nl_l = real(ocean_state%continuity%n_limited_total, wp)
                  call halo_allreduce_sum(nl_l, nl_g)
                  if (compute_rank == 0) then
                     if (nint(nl_g, int64) > n_limited_last_report) then
                        call logger%info("positive-definite limiter:  "// &
                                         to_string(nint(nl_g, int64) - n_limited_last_report)// &
                                         " faces limited this report ("//to_string(nint(nl_g, int64))// &
                                         " total)")
                        n_limited_last_report = nint(nl_g, int64)
                     end if
                  end if
               end block
            end if
         end do
         call profiler_stop("time_loop")

         if (compute_rank == 0) then
            wall_elapsed = wall_timer%get_elapsed_time()
            call logger%info("")
            call logger%info("Solver complete: "//to_string(n_steps)//" steps in "// &
                             to_string(wall_elapsed)//" s")
            call report_throughput(real(cfg%nx, wp)*real(cfg%ny, wp), n_steps, &
                                   wall_elapsed, compute_size, compute_rank)
            ! End-of-run device-memory truth: state is still device-resident
            ! here (exit_data below) — the growth since enter_data attributes
            ! any lazily allocated device workspaces even on short runs.
            call mem_log_device_growth("ocean end of run")
         end if

         ! FIX-1 (spec §(c) mandatory end-of-run flush): a run whose length is
         ! not an exact multiple of dt_tracer_advect_ratio ends mid-window with
         ! frozen hTr and advanced h_layer.  Drain before the clean-end
         ! checkpoint (and before the diag stream closes) so neither the final
         ! restart nor any pending state is emitted with an un-drained window.
         ! No-op at ratio = 1 (bit-identical) and idempotent.
         call ocean_dyn_flush_tracer_window(grid, ocean_state%metrics, &
                                            ocean_state%dyn, ocean_state%continuity, &
                                            ocean_state%multilayer, bc=ocean_state%bc)

         ! Clean-end checkpoint: one final restart so a completed run can be
         ! continued (and the bit-exact gate has a top-of-step snapshot).
         ! State is still device-resident here; the write pulls D->H.
         if (use_restart_output) then
            restart_filename = output_rank_filename(trim(cfg%output_dir), "restart", mpi_rank)
            call ocean_state_restart_write(ocean_state, grid, decomp, &
                                           trim(restart_filename), t_current, n_steps)
         end if

         ! P2.4: diag NetCDF stream close + ocean-halo module teardown + host
         ! deallocation (geo%destroy/state%destroy) all now live in
         ! `engine_teardown`, called below (after exit_data, mirroring the
         ! original relative order of ocean_state_exit_data before
         ! ocean_halo_destroy/geo%destroy/ocean_state%destroy — the ONE
         ! deliberate reordering here is the diag stream close moving from
         ! before exit_data to inside engine_teardown, i.e. after it; that is
         ! a NetCDF-only I/O finalization with no device-state dependency, so
         ! it does not affect the simulation trajectory, only where its call
         ! lands relative to a couple of NVTX profiler range boundaries).
         call profiler_start("ocean_exit_data")
         call engine_exit_data(engine)
         call profiler_stop("ocean_exit_data")

         if (compute_rank == 0) then
            call wall_timer%stop()
            call logger%info("Total steps: "//to_string(n_steps))
            call logger%info("Total wall time: "// &
                             to_string(wall_timer%get_elapsed_time())//" s")
         end if

         call profiler_start("xd_ocean_halo", nvtx_only=.true.)
         call engine_teardown(engine)
         call profiler_stop("xd_ocean_halo")

         ! profiler_report / profiler_end after all D->H teardown so that
         ! ocean_exit_data + xd_ocean_halo ranges appear in the nsys timeline.
         if (compute_rank == 0) then
            call profiler_report("Compute", root_region="time_loop")
         end if
         call profiler_end()

         if (compute_rank == 0) call get_knowledge()

      end associate
   end subroutine driver_run_ocean

   ! ================================================================
   ! Helper subroutines (private to the driver module)
   ! ================================================================

   ! configure_ocean_diag / diag_vgrid_from_name / uniform_diag_levels moved
   ! to `rdb_ocean_engine` (P2.4) — see `engine_configure_diag` there.
   ! `diag_vgrid_from_name` is re-exported below (`use rdb_ocean_engine,
   ! only: diag_vgrid_from_name` + the `public ::` at the top of this
   ! module) so any existing `use rdb_driver, only: diag_vgrid_from_name`
   ! caller is unaffected.

   subroutine print_status_line(step, t, dt, wall, t_end, mass, mean_S, mean_T)
      !! Print a table-formatted status line with estimated remaining time.
      !! When mean_S / mean_T are present, two extra columns are written
      !! between Mass and Wall so the multilayer header lines up.
      integer, intent(in) :: step
      real(wp), intent(in) :: t, dt, wall, t_end
      real(wp), intent(in), optional :: mass
      real(wp), intent(in), optional :: mean_S, mean_T

      real(wp) :: remaining
      character(len=256) :: line

      if (t > 0.0_wp) then
         remaining = wall*(t_end - t)/t
      else
         remaining = 0.0_wp
      end if

      if (present(mass) .and. present(mean_S) .and. present(mean_T)) then
         write (line, "(I9,F14.2,F14.4,ES13.3,F13.4,F13.4,F13.2,F15.2)") &
            step, t, dt, mass, mean_S, mean_T, wall, remaining
      else if (present(mass)) then
         write (line, "(I9,F14.2,F14.4,ES13.3,F13.2,F15.2)") &
            step, t, dt, mass, wall, remaining
      else
         write (line, "(I9,F14.2,F14.4,13X,F13.2,F15.2)") &
            step, t, dt, wall, remaining
      end if
      call logger%info("  "//trim(line))
   end subroutine print_status_line

   subroutine configure_log_level(level_str)
      !! Set logger verbosity from string
      character(len=*), intent(in) :: level_str

      select case (trim(level_str))
      case ("debug")
         call logger%configure(debug_level)
      case ("verbose")
         call logger%configure(verbose_level)
      case ("info")
         call logger%configure(info_level)
      case ("performance")
         call logger%configure(performance_level)
      case ("warning")
         call logger%configure(warning_level)
      case ("error")
         call logger%configure(error_level)
      case default
         call logger%configure(info_level)
      end select
   end subroutine configure_log_level

end module rdb_driver
