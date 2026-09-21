!! Ocean setup / step / teardown engine (P2.4) — the ONE shared lifecycle
!! path for the ocean dyn-core, called by `driver_run_ocean`
!! (`src/driver/rdb_driver.F90`), the C ABI (`src/api/rdb_ocean_api.F90`)
!! and, where it would not distort the benchmark's own measurement,
!! `benchmarks/bench_ocean.F90`.
module rdb_ocean_engine
   !! Before this module, there were THREE copies of the ocean setup
   !! sequence: `driver_run_ocean` ran 21 `configure_ocean_*`-family
   !! stages, `rdb_ocean_create_from_string` ran 9, and `bench_ocean`
   !! ran 9 (mirroring the API). Twelve stages — `configure_ocean_bc`,
   !! `configure_ocean_diag`, `configure_ocean_hdiff`,
   !! `configure_ocean_porous`, `configure_ocean_p_surf`,
   !! `configure_ocean_sponge`, `configure_ocean_tides`,
   !! `configure_ocean_tracers`, `configure_ocean_wave_drag`,
   !! `configure_ocean_wetdry`, `ocean_data_forcing_configure`,
   !! `ocean_halo_init` — were parsed, validated and silently INERT on
   !! the API path: a Python caller could set a boundary condition, a
   !! tide, a sponge, wet/dry, porous barriers, tracer hdiff, surface
   !! pressure, wave drag or file-backed forcing and none of it would
   !! take effect. `engine_setup` below is `driver_run_ocean`'s exact
   !! 21-stage sequence, in its exact order — order is load-bearing
   !! and heavily commented at each site; see the per-stage comments
   !! carried over from the driver.
   !!
   !! `ocean_engine_t` owns the god state (`ocean_state_t`) plus the
   !! setup products that are NOT state slots and used to be driver
   !! stack locals: `grid`, `geo` (geothermal), `bc_source` (constant
   !! boundary-data backend), `decomp`, `evp_params` (sea-ice EVP),
   !! `ic_par` (sea-ice IC), `cfl_vtol` (console-stats CFL-vanish
   !! tolerance) — exactly the things the API path silently lacked a
   !! place to put. `n_inner` (resolved barotropic substep count) and
   !! `diag_enabled` (latched `cfg%ocean%diag%enabled`, so `engine_step`
   !! does not need `cfg`) are new fields the engine itself needs to
   !! function; `device_mapped`/`is_setup` are idempotency guards
   !! mirroring `rdb_handle`'s `ocean_handle_t`.
   !!
   !! `engine_step` + `engine_step_finalize` are intentionally more
   !! than a bare `ocean_dyn_step_split` call: because `engine_setup`
   !! now reaches `configure_ocean_bc`/`_porous`/`_diag`/
   !! `ocean_data_forcing_configure` for every caller, the matching
   !! PER-STEP halves of those stages (`bc_source%update`,
   !! `ocean_porous_refresh`, the file-forcing `update_all`/`apply`
   !! pair, `ocean_surface_flux_assemble`, `ocean_diag_t%step`) must
   !! run too, or those stages would be "configured but not executed"
   !! for API/bench callers — precisely the class of bug this module
   !! exists to close. The two calls are split (not one) because the
   !! driver's own sea-ice block sits, in `driver_run_ocean`'s
   !! sequence, between the dyn-core advance and
   !! `ocean_surface_flux_assemble` — see `engine_step`'s docstring.
   !! What stays OUT of both, by design: `t_current`/`n_steps`
   !! bookkeeping, status/restart cadence, and `console_stats` (the
   !! caller latches its own `mass0`/`salt0`/`heat0`/`ke0` reference on
   !! first call — see `driver_run_ocean`).
   !!
   !! **P2.4b closed the sea-ice gap P2.4 left open.** `engine_step_ice`
   !! below is `driver_run_ocean`'s former inline sea-ice per-step block
   !! (frazil/EVP/thermo-driver/transport), transcribed verbatim —
   !! same call order, same gates (`ice%enable`, `ice%dynamics`, the
   !! thermo cadence via `is_thermo_step()`, `&ocean_ice_nml transport`),
   !! same "MANDATED ORDER" contract documented at its own call sites.
   !! A caller with sea ice enabled calls `engine_step`, then
   !! `engine_step_ice`, then `engine_step_finalize` — exactly
   !! `driver_run_ocean`'s sequence (both the driver and the C ABI's
   !! `rdb_ocean_step` do this now). It is `engine_step_ice`, not
   !! `engine_step`, precisely because it must run BETWEEN the dyn-core
   !! advance and `ocean_surface_flux_assemble`: the assembler reads the
   !! salt/heat components the ice block writes. `engine_step_ice`
   !! itself is a no-op (returns immediately) when
   !! `&ocean_ice_nml enable = .false.` — bit-identical to before this
   !! phase for every non-ice caller/config. A caller with sea ice
   !! enabled via `engine_setup` who never calls `engine_step_ice` is
   !! back to the P2.4 gap (configured, not advanced) — that is now a
   !! caller bug, not a library one.
   !!
   !! NetCDF is optional (`RDB_ENABLE_NETCDF`/`RDB_NO_NETCDF`).
   !! This module compiles unconditionally (the API's `core_objs`
   !! target requires it even in a NetCDF-free configure), so every
   !! NetCDF-only piece — the diag NetCDF stream, file-backed forcing,
   !! and restart I/O's *filename resolution* (`output_rank_filename`,
   !! `rdb_io_netcdf`) — is `#ifndef RDB_NO_NETCDF`-guarded, with an
   !! `#else` `fail()` if a caller actually asks for one of those
   !! features in a build that cannot provide it (previously an
   !! impossible combination: `rdb_driver.F90` itself was excluded from
   !! the NetCDF-free build, so no code path could reach this
   !! contradiction; reachable via the API today).
   use rdb_constants, only: wp
   use rdb_config, only: config_t, resolve_bt_halo, bt_halo_auto_exclusion, &
                         BT_HALO_AUTO_SENTINEL
   use rdb_vcoord, only: parse_remap_method
   use rdb_state, only: register_default_tracers
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init_from_config, decomp_log_summary
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data, &
                              ocean_state_seed_from_cfg, ocean_state_restart_read
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy, ocean_halo_reserve, &
                             ocean_halo_centre, &
                             ocean_halo_bt_group_2d, ocean_halo_face_x, &
                             ocean_halo_bt_group_2d_wide
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_ocean_periodic, only: ocean_periodic_wrap_state, ocean_periodic_wrap_centre_2d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_state, ocean_fold_wrap_eta_2d
   use rdb_ocean_boundary_data, only: ocean_boundary_data_constant_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_set_edges, ocean_bc_state_set_topology
   use rdb_ocean_metrics, only: metrics_assemble_from_supergrid_arrays, metrics_finalize
   use rdb_ocean_dyn, only: ocean_dyn_step, ocean_dyn_step_split, ocean_porous_refresh, &
                            ocean_dyn_enable_bt_wide, isopycnal_vanish_tol
   use rdb_ocean_surface_flux, only: ocean_surface_flux_assemble
   use rdb_ocean_surface_stress, only: ocean_surface_stress_set_shelf_from_ustar
   use rdb_ocean_cavity_flux, only: ocean_cavity_flux_step
   use rdb_ocean_vcoord, only: parse_ocean_vcoord_type, VCOORD_LAGRANGIAN
   use rdb_ocean_setup, only: configure_ocean_metrics, configure_ocean_land_mask, &
                              configure_ocean_forcing, &
                              configure_ocean_drag, &
                              configure_ocean_hdiff, &
                              configure_ocean_vmix, configure_ocean_tracers, configure_ocean_lateral, &
                              configure_ocean_reference_density, &
                              configure_ocean_pgf, configure_ocean_bt, &
                              configure_ocean_bt_split, configure_ocean_bc, &
                              configure_ocean_tides, configure_ocean_p_surf, &
                              configure_ocean_wave_drag, configure_ocean_porous, &
                              configure_ocean_closed_faces, &
                              configure_ocean_k_top, &
                              configure_ocean_cavity, &
                              configure_ocean_cavity_melt, &
                              configure_ocean_top_drag, &
                              configure_ocean_wetdry, &
                              configure_ocean_sponge
   use rdb_ocean_stability_audit, only: ocean_stability_audit
   use rdb_ocean_sponge, only: ocean_sponge_snapshot_reference, ocean_sponge_refresh_target
   use rdb_ocean_geothermal, only: ocean_geothermal_t
   use rdb_ocean_diag_fills, only: set_diag_remap_method, parse_diag_remap_scheme, &
                                   set_diag_mask_vanished
   use rdb_ocean_diag_derived, only: apply_diag_selection
   use rdb_ocean_diag, only: DIAG_VGRID_LAYER, DIAG_VGRID_Z_FIXED, &
                             DIAG_VGRID_SIGMA, DIAG_VGRID_ZSTAR, DIAG_VGRID_DENSITY
#ifndef RDB_NO_NETCDF
   use rdb_ocean_diag_netcdf, only: open_stream, close_stream
   use rdb_ocean_data_input, only: ocean_data_input_update_all
   use rdb_ocean_data_forcing, only: ocean_data_forcing_configure, ocean_data_forcing_apply
   use rdb_io_netcdf, only: output_rank_filename, ensure_directory_exists
#endif
   use rdb_ice_evp, only: ice_evp_params_t, ice_evp_params_from_config, ice_evp_step
   use rdb_ice_init, only: ice_ic_params_t, ice_ic_params_from_config, ice_init_apply, &
                           ICE_IC_CONC_ZERO
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_resume_apply, ice_ocean_stress_cleanup, &
                                    ice_ocean_brine_flux, ice_ocean_heat_flux, &
                                    ice_ocean_sw_flux, ice_ocean_stress_flux
   use rdb_ice_frazil, only: ice_frazil_accumulate
   use rdb_ice_frazil_uptake, only: ice_frazil_uptake
   use rdb_ice_atm_forcing, only: ice_atm_forcing_restoring
   use rdb_ice_basal_flux, only: ice_compute_basal_flux
   use rdb_ice_thermo_driver, only: ice_thermo_driver_step
   use rdb_ice_snow, only: ice_snowfall_ocean_share
   use rdb_ice_itd, only: ice_adjust_categories
   use rdb_ice_transport, only: ice_transport_step
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP, OCEAN_STATUS_ERR_BAD_SHAPE
   use rdb_error_ring, only: fail
   use rdb_profiler, only: profiler_start, profiler_stop
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_engine_t
   public :: engine_setup
   public :: engine_enter_data
   public :: engine_step
   public :: engine_step_ice
   public :: engine_step_finalize
   public :: engine_exit_data
   public :: engine_teardown
   public :: diag_vgrid_from_name

   type :: ocean_engine_t
      !! One ocean simulation's setup products, minus the caller-owned
      !! time-loop bookkeeping (t_current/n_steps/status-restart
      !! cadence/console_stats stay with the caller — see the module
      !! docstring).
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      type(ocean_geothermal_t) :: geo
      type(ocean_boundary_data_constant_t) :: bc_source
         !! Boundary data source, refreshed once per outer step via
         !! `engine_step` -> `bc_source%update(t, state%bc)`.
      type(decomp_t) :: decomp
      type(ice_evp_params_t) :: evp_params
         !! Sea-ice EVP physical + numerical parameters, built once from
         !! `&ocean_ice_nml` at setup.
      type(ice_ic_params_t) :: ic_par
         !! Sea-ice initial-condition parameters, built once from
         !! `&ocean_ice_ic_nml` at setup.
      real(wp) :: cfl_vtol = 0.0_wp
         !! Console-stats MaxCFL vanish tolerance; 0 when off. Consumed
         !! by the caller's own `ocean_console_stats_report`, not by
         !! the engine itself.
      integer :: n_inner = 0
         !! Resolved barotropic fast-loop substep count
         !! (`cfg%ocean%bt%n_inner`, possibly auto-derived by
         !! `configure_ocean_bt_split`). `engine_step` dispatches to
         !! the split-RK2 path when `n_inner >= 1`, else the legacy
         !! unsplit `ocean_dyn_step` — mirrors `driver_run_ocean`'s own
         !! dispatch exactly; `engine_setup` does not itself require
         !! `n_inner >= 1` (a caller wanting that guarantee, such as
         !! the C ABI, checks it after `engine_setup` returns).
      logical :: diag_enabled = .false.
         !! Latched `cfg%ocean%diag%enabled` at setup, so `engine_step`
         !! can gate its `state%diag%step` call without needing `cfg`.
      logical :: device_mapped = .false.
         !! True between `engine_enter_data` and `engine_exit_data`.
      logical :: is_setup = .false.
         !! True once `engine_setup` has completed (host-side only;
         !! device mapping is a separate step).

      ! ---- P2.5: pre-create geometry injection staging ----
      ! Populated (by the C ABI's rdb_ocean_stage_* calls, or directly by
      ! a Fortran caller) BEFORE engine_setup runs, and consumed inside it —
      ! geometry must land before ocean_state_enter_data, exactly like
      ! passive-tracer registration (registry_locked closes at enter_data).
      ! Unset (not `allocated`/`.false.`) ⇒ engine_setup falls through to
      ! its ordinary namelist-driven path, byte-identical to before P2.5.
      real(wp), allocatable :: staged_bathymetry(:, :)
         !! Interior-sized (nx_phys, ny_phys) bathymetry, RAW in the
         !! caller's own sign convention (`staged_bathymetry_convention`).
         !! Consumed by `ocean_state_seed_from_cfg`'s `injected_b` argument,
         !! which overrides `cfg%ocean%topo%topo_config` entirely.
      integer :: staged_bathymetry_convention = 0
         !! `BATHY_CONVENTION_*` (`rdb_ocean_bathymetry_inject`). Only
         !! meaningful when `staged_bathymetry` is allocated — no default
         !! (D6.2: sign is the single most dangerous argument here).
      real(wp), allocatable :: staged_metrics_x(:, :), staged_metrics_y(:, :)
      real(wp), allocatable :: staged_metrics_dx(:, :), staged_metrics_dy(:, :)
      real(wp), allocatable :: staged_metrics_area(:, :)
         !! In-memory MOM6-style supergrid arrays
         !! (`metrics_assemble_from_supergrid_arrays`'s `sg_x/sg_y/sg_dx/
         !! sg_dy/sg_area`). When `staged_metrics_x` is allocated,
         !! `engine_setup` assembles metrics directly from these instead of
         !! dispatching on `cfg%ocean%grid%grid_config`.
      logical :: has_staged_topology = .false.
         !! True after a topology injection call; `staged_periodic_x/y`
         !! are only meaningful when this is true.
      logical :: staged_periodic_x = .false.
      logical :: staged_periodic_y = .false.
         !! Oceananigans-style "the grid owns periodicity": applied via
         !! `ocean_bc_state_set_topology` right after `configure_ocean_bc`,
         !! overriding whatever the namelist edge tags derived.
   end type ocean_engine_t

contains

   ! ================================================================
   ! Setup — driver_run_ocean's exact 21-stage sequence
   ! ================================================================

   subroutine engine_setup(engine, cfg, ierr, compute_rank, compute_size, mpi_rank, &
                           restart_file, t_restart, step_restart)
      !! Host-side setup: decomposition -> grid -> god state -> IC seed
      !! -> restart (optional) -> the 21 `configure_ocean_*`-family
      !! stages -> ghost wraps -> halo init -> land mask -> wave
      !! drag/porous -> sea-ice IC. Does NOT map the state onto the
      !! device — call `engine_enter_data` next. Preserves
      !! `driver_run_ocean`'s exact stage order; see the per-stage
      !! comments below (carried over from the driver almost verbatim).
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero (`rdb_ocean_status` code) on a bad config/setup/IC
         !! failure when present; absent behaves as today — the
         !! offending `configure_ocean_*` stage (or this routine's own
         !! pre-flight checks) `error stop`s.
      integer, intent(in), optional :: compute_rank
         !! This rank's 0-based compute-rank index. Default 0
         !! (single-rank: the API/bench callers).
      integer, intent(in), optional :: compute_size
         !! Total compute-rank count. Default 1.
      integer, intent(in), optional :: mpi_rank
         !! Rank used for the restart filename convention
         !! (`output_rank_filename`). Default = `compute_rank` (the
         !! ocean path has no separate I/O-server rank).
      character(len=*), intent(in), optional :: restart_file
         !! Warm-restart source (file, or a run-directory — see
         !! `driver_run_ocean`'s `.nc`-suffix convention). Absent/blank
         !! => cold start (matches the API/bench callers today).
      real(wp), intent(out), optional :: t_restart
         !! Restored simulation time (0 on a cold start).
      integer, intent(out), optional :: step_restart
         !! Restored outer-step count (0 on a cold start).

      integer :: rank, csize, mrank
      character(len=512) :: restart_filename
      real(wp) :: t_restart_local
      integer :: step_restart_local
      logical :: did_restart
      logical :: bt_excluded
      character(len=:), allocatable :: bt_excl_reason
      integer :: bt_halo_req, bt_halo_res

      rank = 0
      if (present(compute_rank)) rank = compute_rank
      csize = 1
      if (present(compute_size)) csize = compute_size
      mrank = rank
      if (present(mpi_rank)) mrank = mpi_rank

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      t_restart_local = 0.0_wp
      step_restart_local = 0

      ! ---- pre-flight (was comm_env_abort in the driver; a library
      ! cannot abort its host process, so these are now fail()-based) ----
      if (cfg%dt_fixed <= 0.0_wp) then
         call fail("sim_type='ocean' requires dt_fixed > 0 (no adaptive CFL helper yet).", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%px*cfg%py /= csize) then
         call fail("Process grid px*py = "//to_string(cfg%px*cfg%py)// &
                   " does not match number of compute ranks = "//to_string(csize)// &
                   " (ocean path)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (csize > 1) then
         ! D6: tripolar north fold requires the north rank-row undecomposed
         ! in x; allow it only when px == 1.
         if (trim(cfg%ocean%bc%north) == "tripolar_fold" .and. cfg%px > 1) then
            call fail("Tripolar north fold requires the north rank-row "// &
                      "undecomposed in x (px must be 1) — deferred D6.", &
                      ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! E1: windowed tracer advect drain halo not yet wired for multi-rank.
         if (cfg%ocean%vmix%dt_tracer_advect_ratio > 1) then
            call fail("dt_tracer_advect_ratio > 1 multi-rank drain halo "// &
                      "is deferred (E1); run single-rank or set ratio = 1.", &
                      ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if

      call decomp_init_from_config(engine%decomp, cfg, csize, rank)
      if (rank == 0 .and. csize > 1) call decomp_log_summary(engine%decomp, csize)

      ! Resolve the BT wide-halo march-in sentinel now that compute_size and
      ! the feature flags are both known (see rdb_config::resolve_bt_halo).
      bt_halo_req = cfg%ocean%bt%bt_halo
      call bt_halo_auto_exclusion(cfg, bt_excluded, bt_excl_reason)
      bt_halo_res = resolve_bt_halo(bt_halo_req, csize, bt_excluded)
      cfg%ocean%bt%bt_halo = bt_halo_res
      if (rank == 0 .and. bt_halo_req == BT_HALO_AUTO_SENTINEL) then
         if (bt_halo_res > 0) then
            call logger%info("BT march-in: bt_halo auto -> "// &
                             to_string(bt_halo_res)//" (multi-rank)")
         else if (csize <= 1) then
            call logger%info("BT march-in: bt_halo auto -> 0 (serial)")
         else
            call logger%info("BT march-in: bt_halo auto -> 0 "// &
                             "(auto disabled: "//trim(bt_excl_reason)//" active)")
         end if
      end if

      call engine%grid%init(engine%decomp%nx_local, engine%decomp%ny_local, &
                            cfg%nghost, cfg%dx, cfg%dy)
      engine%grid%i_offset_global = engine%decomp%i_start - 1
      engine%grid%j_offset_global = engine%decomp%j_start - 1
      engine%grid%nx_global = engine%decomp%nx_global
      engine%grid%ny_global = engine%decomp%ny_global
      call engine%state%init_from_config(cfg, engine%grid)

      ! Wire vcoord parameters BEFORE the seed runs — the seed calls
      ! `vcoord%build_zref_full(b)` at its tail, which reads
      ! `zstar_h_surf_target` etc.
      if (engine%state%use_multilayer) then
         engine%state%vcoord%coord_type = parse_ocean_vcoord_type(cfg%vcoord_type)
         engine%state%vcoord%remap_method = parse_remap_method(cfg%remap_method)
         engine%state%vcoord%zstar_h_surf_target = cfg%zstar_h_surf_target
         engine%state%vcoord%zstar_h_min = cfg%zstar_h_min
      end if

      ! Analytical IC from cfg scalars — or, when P2.5 geometry injection
      ! staged a bathymetry array (rdb_ocean_stage_bathymetry), that
      ! array overrides cfg%ocean%topo%topo_config entirely.
      if (allocated(engine%staged_bathymetry)) then
         call ocean_state_seed_from_cfg(engine%state, engine%grid, cfg, ierr=ierr, &
                                        injected_b=engine%staged_bathymetry, &
                                        injected_b_convention=engine%staged_bathymetry_convention)
      else
         call ocean_state_seed_from_cfg(engine%state, engine%grid, cfg, ierr=ierr)
      end if
      if (setup_failed(ierr)) return

      ! PR-23 sponge: snapshot the reference state (target_source="ic")
      ! from the JUST-SEEDED initial condition, BEFORE any warm-restart
      ! read below (load-bearing ordering).
      call ocean_sponge_snapshot_reference(engine%state%sponge, engine%grid, engine%state%multilayer)

      ! Overlay tracer metadata + scalar IC values from cfg.
      call register_default_tracers( &
         engine%state%multilayer%tracers(engine%state%multilayer%idx_salinity), &
         engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature), &
         cfg)

      ! Dynamic wetting/drying: MUST precede the restart read (persistent
      ! hysteresis mask must be allocated + registered before the restart
      ! registry walk). Default off => no-op, byte-identical.
      call configure_ocean_wetdry(cfg, engine%state, engine%grid, rank)

      ! Warm restart (optional — absent/blank restart_file => cold start,
      ! matching the API/bench callers today). NetCDF-only: filename
      ! resolution needs `output_rank_filename`.
      did_restart = .false.
      if (present(restart_file)) then
         if (len_trim(restart_file) > 0) then
            did_restart = .true.
#ifndef RDB_NO_NETCDF
            if (len_trim(restart_file) >= 3 .and. &
                restart_file(len_trim(restart_file) - 2:len_trim(restart_file)) == ".nc") then
               restart_filename = trim(restart_file)
            else
               restart_filename = output_rank_filename(trim(restart_file), "restart", mrank)
            end if
            call ocean_state_restart_read(engine%state, engine%grid, engine%decomp, &
                                          trim(restart_filename), t_restart_local, &
                                          step_restart_local, ierr=ierr)
            if (setup_failed(ierr)) return
            if (rank == 0) then
               call logger%info("Ocean warm restart from "//trim(restart_filename)// &
                                " at t = "//to_string(t_restart_local)//" s (step "// &
                                to_string(step_restart_local)//")")
            end if
#else
            call fail("engine_setup: restart_file was given but this build has no "// &
                      "NetCDF support (RDB_ENABLE_NETCDF=OFF)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
#endif
         end if
      end if
      if (present(t_restart)) t_restart = t_restart_local
      if (present(step_restart)) step_restart = step_restart_local

      ! Diag-manager: optional z-levels, default variable set, per-rank
      ! NetCDF stream. NetCDF-only (moved from driver_run_ocean's private
      ! helper of the same name).
      engine%diag_enabled = cfg%ocean%diag%enabled
      call engine_configure_diag(cfg, engine%state, mrank, engine%decomp, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Surface heat/salt/p_surf/sw-penetration/restore seeding — 2D
      ! fields device-mapped by ocean_state_enter_data; re-seeded on
      ! resume (configure-time static, not in the restart registry).
      call engine%state%surface_flux%set_components(engine%grid, &
                                                    cfg%ocean%forcing%enable_components)
      if (rank == 0 .and. cfg%ocean%forcing%enable_components) then
         call logger%info("Forcing components: ON")
      end if
      call engine%state%surface_flux%set_surface_flux_const( &
         cfg%ocean%thermo%q_heat, cfg%ocean%thermo%q_salt)
      call engine%state%surface_flux%set_p_surf_const( &
         cfg%ocean%psurf%p_surf_const)
      ! E3 top-of-column IN-SITU EOS pressure: seed `ms%p_top` from the
      ! assembled `sf%p_surf` HERE, on the host and before `enter_data`, so
      ! the very first PGF of the run already sees the load.  The
      ! outer-step driver refreshes it every step from the same source,
      ! which is what keeps it live once the ice mass-loading PR makes
      ! `p_surf` dynamic.  Gated: with `in_eos = .false.` `p_top` stays the
      ! zero array it was allocated as and every EOS evaluation is
      ! bit-identical.
      !
      ! P5.0 joins `&ocean_pgf_nml p_top_in_bc` — the load in the FV_MOM6
      ! pressure-stack surface BC — to the SAME seed, so the unsplit
      ! driver (which has no per-step refresh) still gets the configure
      ! value rather than a zero array, and the split driver's step-1 PGF
      ! is already loaded.
      if ((cfg%ocean%psurf%in_eos .or. cfg%ocean%pgf%p_top_in_bc) .and. &
          allocated(engine%state%surface_flux%p_surf)) then
         engine%state%multilayer%p_top = engine%state%surface_flux%p_surf
      end if
      call engine%state%surface_flux%set_sw_penetration( &
         cfg%ocean%thermo%sw_pen_frac, cfg%ocean%thermo%sw_band_ratio, &
         cfg%ocean%thermo%sw_zeta1, cfg%ocean%thermo%sw_zeta2, &
         sw_source=cfg%ocean%thermo%sw_source)
      call engine%state%surface_flux%set_restore( &
         cfg%ocean%restore%enable_restore_temp, &
         cfg%ocean%restore%enable_restore_salt, &
         cfg%ocean%restore%piston_t, cfg%ocean%restore%piston_s, &
         cfg%ocean%restore%restore_sst, cfg%ocean%restore%restore_sss)
      if (rank == 0 .and. engine%state%surface_flux%has_restore_T) then
         call logger%info("Restore SST:      piston = "// &
                          to_string(cfg%ocean%restore%piston_t)//" m/day, target = "// &
                          to_string(cfg%ocean%restore%restore_sst)//" degC")
      end if
      if (rank == 0 .and. engine%state%surface_flux%has_restore_S) then
         call logger%info("Restore SSS:      piston = "// &
                          to_string(cfg%ocean%restore%piston_s)//" m/day, target = "// &
                          to_string(cfg%ocean%restore%restore_sss)//" PSU")
      end if
      if (rank == 0 .and. &
          (cfg%ocean%thermo%q_heat /= 0.0_wp .or. cfg%ocean%thermo%q_salt /= 0.0_wp)) then
         call logger%info("Surface flux:     Q_heat = "// &
                          to_string(cfg%ocean%thermo%q_heat)//" W/m2, Q_salt = "// &
                          to_string(cfg%ocean%thermo%q_salt)//" kg/m2/s")
      end if

      ! Sea-ice: resume-fold the restart-carried brine/heat/shortwave
      ! contributions back into the just-reseeded Q_salt/Q_heat (cold
      ! start: exact +0.0).
      if (engine%state%ice%enable) then
         engine%state%surface_flux%Q_salt = engine%state%surface_flux%Q_salt &
                                            + engine%state%ice%salt_flux_diag
         engine%state%surface_flux%has_salt = .true.
         engine%state%surface_flux%Q_heat = engine%state%surface_flux%Q_heat &
                                            + engine%state%ice%heat_flux_diag &
                                            + engine%state%ice%sw_thru_diag
         engine%state%surface_flux%has_heat = .true.
         if (engine%state%surface_flux%use_components) then
            engine%state%surface_flux%has_q_sw = .true.
         end if
      end if

      ! Geothermal bottom heat flux (held on the engine, not a state slot;
      ! the split-driver's `geo` arg is optional).
      call engine%geo%init(engine%grid)
      engine%geo%enable = cfg%ocean%geothermal%enable
      engine%geo%q_geo_const = cfg%ocean%geothermal%q_geo
      if (rank == 0 .and. engine%geo%enable .and. cfg%ocean%geothermal%q_geo /= 0.0_wp) then
         call logger%info("Geothermal flux:  Q_geo = "// &
                          to_string(cfg%ocean%geothermal%q_geo)//" W/m2")
      end if

      ! P2.5 geometry injection: in-memory supergrid arrays
      ! (rdb_ocean_stage_metrics) bypass the cfg%ocean%grid%grid_config
      ! dispatch entirely — metrics_assemble_from_supergrid_arrays is the
      ! same battle-tested index-sum path metrics_fill_from_supergrid /
      ! the tripolar generator use, just fed in-memory arrays instead of a
      ! mosaic file.
      if (allocated(engine%staged_metrics_x)) then
         if (size(engine%staged_metrics_x, 1) /= 2*engine%grid%nx_phys + 1 .or. &
             size(engine%staged_metrics_x, 2) /= 2*engine%grid%ny_phys + 1) then
            call fail("engine_setup: staged metrics x/y shape mismatch — expected "// &
                      "(2*nx_phys+1, 2*ny_phys+1) = ("//to_string(2*engine%grid%nx_phys + 1)// &
                      ","//to_string(2*engine%grid%ny_phys + 1)//")", ierr, OCEAN_STATUS_ERR_BAD_SHAPE)
            return
         end if
         call metrics_assemble_from_supergrid_arrays(engine%state%metrics, engine%grid, &
                                                     engine%staged_metrics_x, engine%staged_metrics_y, &
                                                     engine%staged_metrics_dx, engine%staged_metrics_dy, &
                                                     engine%staged_metrics_area)
         call metrics_finalize(engine%state%metrics)
         if (rank == 0) then
            call logger%info("Grid config:      injected (in-memory supergrid arrays)")
         end if
         if (present(ierr)) ierr = OCEAN_STATUS_OK
      else
         call configure_ocean_metrics(cfg, engine%state, engine%grid, rank, ierr=ierr)
         if (setup_failed(ierr)) return
      end if

      ! Equilibrium body-force tide (C1): needs the filled geolatT/geolonT
      ! from configure_ocean_metrics.
      call configure_ocean_tides(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Atmospheric surface-pressure loading / inverse barometer (PR-17).
      call configure_ocean_p_surf(cfg, engine%state, engine%grid, rank)

      call configure_ocean_forcing(cfg, engine%state, engine%grid, rank, decomp=engine%decomp)

      ! File-backed surface forcing (PR-15): must run after set_components
      ! (heat/salt destination) and after configure_ocean_forcing (whose
      ! formula wind it overrides), before enter_data (registration
      ! allocates the reader's bracket buffers). NetCDF-only.
#ifndef RDB_NO_NETCDF
      call ocean_data_forcing_configure(cfg%ocean%dataovr, engine%state%data_input, engine%grid, &
                                        engine%state%surface_stress, engine%state%surface_flux, &
                                        engine%state%bc, engine%state%data_forcing, ierr=ierr)
      if (setup_failed(ierr)) return
#else
      if (cfg%ocean%dataovr%enable) then
         call fail("engine_setup: &ocean_dataovr_nml enable = .true. but this build has no "// &
                   "NetCDF support (RDB_ENABLE_NETCDF=OFF)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
#endif

      call configure_ocean_drag(cfg, engine%state, rank, ierr=ierr)
      if (setup_failed(ierr)) return
      call configure_ocean_hdiff(cfg, engine%state, rank)
      call configure_ocean_vmix(cfg, engine%state, rank, ierr=ierr)
      if (setup_failed(ierr)) return
      call configure_ocean_tracers(cfg, engine%state, rank)
      call configure_ocean_lateral(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Sea-ice PR 5: EVP params + the atmospheric-stress snapshot ice
      ! feels (D7). MANDATED ORDER: the tau_a snapshot must read the
      ! pristine wind BEFORE the resume apply overwrites
      ! surface_stress%tau_x/y with the restart-carried blend.
      if (engine%state%ice%enable .and. engine%state%ice%dynamics) then
         engine%evp_params = ice_evp_params_from_config(cfg%ocean%ice%p0, cfg%ocean%ice%c0, &
                                                        cfg%ocean%ice%ec, cfg%ocean%ice%cdw, &
                                                        cfg%ocean%ice%rho_ocean, &
                                                        cfg%ocean%ice%del_sh_min_scale, &
                                                        cfg%ocean%ice%tdamp, cfg%ocean%ice%evp_sub_steps, &
                                                        cfg%ocean%ice%a_face_stress, &
                                                        cfg%ocean%ice%cfl_trunc, &
                                                        cfg%ocean%ice%cfl_trunc_dyn_its, &
                                                        cfg%ocean%ice%project_ci)
         engine%state%ice%tau_a_x = engine%state%surface_stress%tau_x
         engine%state%ice%tau_a_y = engine%state%surface_stress%tau_y
         call ice_ocean_stress_resume_apply(engine%state%surface_stress, engine%state%ice)
      end if

      call configure_ocean_pgf(cfg, engine%state, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! The single configured ρ₀ (`&ocean_ic_nml rho_0` -> `eos%rho0`) out to
      ! the slots that keep their own copy: surface heat/salt flux, wind
      ! stress, KPP-vmix and the engine-held geothermal slot.  Runs after
      ! every slot-specific configure (nothing downstream re-derives a
      ! reference density) and well before `ocean_state_enter_data`, which is
      ! what makes the on-device `vmix%rho0` reads correct without an
      ! explicit `!$acc update device` — see the routine's docstring.
      call configure_ocean_reference_density(engine%state, geo=engine%geo)

      call configure_ocean_bt(cfg, engine%state, engine%grid, rank)

      call configure_ocean_bt_split(cfg, engine%state, engine%grid, rank)   ! may auto-set n_inner
      ! BC config: edge tags + tidal constituents + tracer inflow values.
      call configure_ocean_bc(cfg, engine%state, rank, ierr=ierr)
      if (setup_failed(ierr)) return
      ! Physical-domain-edge flags: at px*py=1 all four are .true.
      call ocean_bc_state_set_edges(engine%state%bc, engine%decomp%has_west, &
                                    engine%decomp%has_east, engine%decomp%has_south, &
                                    engine%decomp%has_north)
      ! P2.5 geometry injection: Oceananigans-style "the grid owns
      ! periodicity" (rdb_ocean_stage_topology) — overrides whatever
      ! configure_ocean_bc just derived from the namelist edge tags. Must
      ! run BEFORE the periodic ghost wrap / ocean_halo_init /
      ! configure_ocean_land_mask below, all of which read periodic_x/y.
      if (engine%has_staged_topology) then
         call ocean_bc_state_set_topology(engine%state%bc, engine%staged_periodic_x, &
                                          engine%staged_periodic_y, ierr=ierr)
         if (setup_failed(ierr)) return
      end if
      ! PR-23 sponge: build the per-cell idamp maps from the just-configured
      ! edge tags. Must run AFTER configure_ocean_bc, BEFORE enter_data.
      call configure_ocean_sponge(cfg, engine%state, engine%grid, rank)
      ! Seed the boundary data source from the same config so the per-step
      ! update is idempotent with configure_ocean_bc.
      engine%bc_source%u_west = cfg%ocean%bc%west_clamped_u
      engine%bc_source%u_east = cfg%ocean%bc%east_clamped_u
      engine%bc_source%v_south = cfg%ocean%bc%south_clamped_v
      engine%bc_source%v_north = cfg%ocean%bc%north_clamped_v
      engine%bc_source%eta_west = cfg%ocean%bc%west_clamped_eta
      engine%bc_source%eta_east = cfg%ocean%bc%east_clamped_eta
      engine%bc_source%eta_south = cfg%ocean%bc%south_clamped_eta
      engine%bc_source%eta_north = cfg%ocean%bc%north_clamped_eta

      ! Init-time periodic ghost wrap: after topo+IC are seeded and BC
      ! tags are configured, wrap bathymetry + multilayer state ghost
      ! cells so all kernels see correct periodic seam values on the
      ! first step. Run on the HOST here (before enter_data).
      if (engine%state%bc%periodic_x .or. engine%state%bc%periodic_y .or. &
          engine%state%bc%north_fold) then
         call ocean_periodic_wrap_centre_2d( &
            engine%state%barotropic%b, &
            engine%grid%nx_total, engine%grid%ny_total, &
            engine%grid%nx_phys, engine%grid%ny_phys, engine%grid%nghost, &
            engine%state%bc%periodic_x, engine%state%bc%periodic_y)
         call ocean_periodic_wrap_state(engine%grid, engine%state%bc, engine%state%multilayer)
         call ocean_fold_wrap_eta_2d(engine%grid, engine%state%bc, engine%state%barotropic%b)
         call ocean_fold_wrap_state(engine%grid, engine%state%bc, engine%state%multilayer)
         ! The ice draft is bathymetry-class static geometry, so it takes
         ! the bathymetry's ghost treatment VERBATIM: the analytic setter
         ! already filled the ghosts, and the wrap/fold then overwrites
         ! them with the seam-correct values on a periodic or folded edge.
         ! `bt_H_ref = b - z_draft` was latched from the UNWRAPPED pair,
         ! which is exactly why all three are re-wrapped here (and why
         ! they must be re-wrapped TOGETHER — a draft whose seam disagreed
         ! with the datum's would count the ice load twice at that face).
         if (engine%state%metrics%use_cavity) then
            call ocean_periodic_wrap_centre_2d( &
               engine%state%metrics%z_draft, &
               engine%grid%nx_total, engine%grid%ny_total, &
               engine%grid%nx_phys, engine%grid%ny_phys, engine%grid%nghost, &
               engine%state%bc%periodic_x, engine%state%bc%periodic_y)
            call ocean_fold_wrap_eta_2d(engine%grid, engine%state%bc, &
                                        engine%state%metrics%z_draft)
            call ocean_periodic_wrap_centre_2d( &
               engine%state%metrics%cover_frac, &
               engine%grid%nx_total, engine%grid%ny_total, &
               engine%grid%nx_phys, engine%grid%ny_phys, engine%grid%nghost, &
               engine%state%bc%periodic_x, engine%state%bc%periodic_y)
            call ocean_fold_wrap_eta_2d(engine%grid, engine%state%bc, &
                                        engine%state%metrics%cover_frac)
         end if
         ! bt_H_ref was snapshotted from the UNWRAPPED b inside
         ! configure_ocean_bt_split (above) — re-wrap it too.
         if (cfg%ocean%bt%n_inner >= 1) then
            call ocean_periodic_wrap_centre_2d( &
               engine%state%dyn%bt_work%bt_H_ref, &
               engine%grid%nx_total, engine%grid%ny_total, &
               engine%grid%nx_phys, engine%grid%ny_phys, engine%grid%nghost, &
               engine%state%bc%periodic_x, engine%state%bc%periodic_y)
            call ocean_fold_wrap_eta_2d(engine%grid, engine%state%bc, &
                                        engine%state%dyn%bt_work%bt_H_ref)
         end if
      end if

      ! Initialise the ocean halo module. Must run BEFORE any ocean_halo_*
      ! call, with the BC periodicity flags already set.
      call ocean_halo_init(engine%decomp, engine%grid%nghost, &
                           engine%state%bc%periodic_x, engine%state%bc%periodic_y, ierr=ierr)
      if (setup_failed(ierr)) return
      call ocean_halo_reserve(cfg%nz_layers, &
                              merge(engine%grid%nghost + cfg%ocean%bt%bt_halo, 0, &
                                    cfg%ocean%bt%bt_halo > 0), ierr=ierr)
      if (setup_failed(ierr)) return

      ! Host-side seam ghost fill (D0 init-halo, O2): single-rank
      ! non-periodic ⇒ no-op, periodic ⇒ local wrap.
      call ocean_halo_centre(engine%state%barotropic%b, device_resident=.false.)
      if (engine%state%metrics%use_cavity) then
         call ocean_halo_centre(engine%state%metrics%z_draft, device_resident=.false.)
         call ocean_halo_centre(engine%state%metrics%cover_frac, device_resident=.false.)
      end if
      call ocean_halo_exchange_ml_state(engine%state%multilayer, device_resident=.false.)
      if (cfg%ocean%bt%n_inner >= 1) then
         call ocean_halo_centre(engine%state%dyn%bt_work%bt_H_ref, device_resident=.false.)
      end if

      ! Static land masking: derive the C-grid face/corner masks from the
      ! seeded wet_mask + zero the 6 face metrics at land faces.
      call configure_ocean_land_mask(cfg, engine%state, engine%grid, rank)

      ! Configure-time stability audit (viscous CFL / kappa_h diffusive
      ! number / Munk-layer resolution / ah_max-clamps-nu_h): MUST run
      ! after configure_ocean_metrics + configure_ocean_land_mask (needs
      ! the real per-cell metric arrays, not nominal &grid_nml dx/dy) and
      ! before enter_data. See rdb_ocean_stability_audit.F90 for the
      ! motivating failure (a global tripolar aquaplanet NaN, diagnosed
      ! only after the fact — this audit is the fix).
      ! `bt_H_ref` (b - z_draft afloat, 0 where grounded) is the COLUMN the
      ! vertical coordinate divides, so it — not the bathymetry — is what
      ! the terrain-following stiffness check must see under an ice shelf.
      ! It exists only once the barotropic datum has been built.
      if (cfg%ocean%bt%n_inner >= 1) then
         call ocean_stability_audit(cfg, engine%state%metrics, engine%grid, rank, &
                                    ierr=ierr, &
                                    column=engine%state%dyn%bt_work%bt_H_ref)
      else
         call ocean_stability_audit(cfg, engine%state%metrics, engine%grid, rank, ierr=ierr)
      end if
      if (setup_failed(ierr)) return

      ! Barotropic linear wave drag: host-side r_H map + h->face average.
      ! AFTER bathymetry + land masking, BEFORE enter_data.
      call configure_ocean_wave_drag(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Porous barriers: grow + fill the static along-face subgrid
      ! statistics. Same ordering constraints as wave drag.
      call configure_ocean_porous(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Ice-shelf cavity: build the static isostatic load from the draft
      ! the IC seed already laid down, and assert the counted-once datum
      ! invariant. AFTER configure_ocean_pgf (it needs the PGF reference
      ! density) and configure_ocean_bt_split (it checks that latch),
      ! BEFORE enter_data.
      call configure_ocean_cavity(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Ice-shelf basal melt (P2b): copy the thermodynamic knobs onto the
      ! melt slot, resolve the gamma_s sentinel, build the per-column
      ! Coriolis array the hj99 law needs, and SEED ms%p_top from the
      ! isostatic load configure_ocean_cavity just built.  Immediately
      ! after it (that is where p_ice_ref comes from), BEFORE enter_data.
      call configure_ocean_cavity_melt(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! Ice-shelf TOP drag (Phase 4a): coefficients + the static FACE
      ! cover masks projected from `metrics%cover_frac`.  AFTER the
      ! cover_frac halo exchange and after configure_ocean_cavity_melt
      ! (it owns the one-C_d rule against the melt slot's cdrag_top),
      ! BEFORE enter_data (the face masks are host-filled and reach the
      ! device on the slot's `copyin` map).
      call configure_ocean_top_drag(cfg, engine%state, engine%grid, rank)

      ! Partial-step z-level face closure (&vcoord_nml
      ! zfixed_closed_faces).  LAST of the static-geometry builders and
      ! still BEFORE enter_data: it needs `vcoord%z_top` (configure_ocean
      ! _cavity), `bt_work%bt_H_ref` (configure_ocean_bt_split), the
      ! land-masked `dy_cu`/`dx_cv` (configure_ocean_land_mask) and the
      ! periodic-wrap + halo pass above — the mask is built from the
      ! z_fixed target at eta = 0, and its ghost-band correctness IS the
      ! seam correctness of those two inputs.  Knob off => literal no-op.
      call configure_ocean_closed_faces(cfg, engine%state, engine%grid, rank, ierr=ierr)
      if (setup_failed(ierr)) return

      ! The shared FIRST-LIVE-LAYER index `ms%k_top` (+ its two face
      ! twins) — what every top-side consumer reads instead of spelling
      ! `nz`, so that a `z_fixed` column whose top layers are inert
      ! fillers inside the ice draft is forced on the ice-adjacent LIVE
      ! layer and not on the filler.  Same inputs and same ordering
      ! constraints as the closed-face mask above (it is built from the
      ! same `z_fixed` target at eta = 0), still before enter_data.
      ! Literal no-op on every coordinate but `z_fixed` under a cavity:
      ! the arrays already hold the `nz` fallback.
      call configure_ocean_k_top(cfg, engine%state, engine%grid, rank)

      ! Sea-ice PR 24: analytic IC path. Host-side, run once, AFTER
      ! wet_mask/geolatT/wet_T are valid, BEFORE enter_data. Skips on a
      ! warm restart (the restart read already replaced the IC).
      if (engine%state%ice%enable) then
         if (.not. did_restart) then
            engine%ic_par = ice_ic_params_from_config(cfg%ocean%ice_ic%conc_config, &
                                                      cfg%ocean%ice_ic%conc, &
                                                      cfg%ocean%ice_ic%h_ice, cfg%ocean%ice_ic%h_snow, &
                                                      cfg%ocean%ice_ic%t_ice, cfg%ocean%ice_ic%s_ice, &
                                                      cfg%ocean%ice_ic%arctic_edge, &
                                                      cfg%ocean%ice_ic%antarctic_edge)
            call ice_init_apply(engine%grid, engine%state%multilayer, engine%state%metrics, &
                                engine%state%ice, engine%ic_par)
            if (rank == 0 .and. engine%ic_par%conc_config /= ICE_IC_CONC_ZERO) then
               call logger%info("Sea-ice IC:       conc_config='"// &
                                trim(cfg%ocean%ice_ic%conc_config)// &
                                "', h_ice = "//to_string(cfg%ocean%ice_ic%h_ice)// &
                                " m, conc = "//to_string(cfg%ocean%ice_ic%conc))
            end if
         end if
      end if

      ! Resolve + store n_inner (does NOT fail loud here — see the
      ! n_inner docstring on ocean_engine_t; a caller wanting the C
      ! ABI's stricter "n_inner must be >= 1" contract checks it itself).
      engine%n_inner = cfg%ocean%bt%n_inner

      ! Phase-3 vanish_tol for console MaxCFL gating.
      engine%cfl_vtol = 0.0_wp
      if (engine%state%dyn%cfl_ignore_vanished) then
         if (engine%state%vcoord%coord_type == VCOORD_LAGRANGIAN) then
            engine%cfl_vtol = isopycnal_vanish_tol(engine%state%dyn%angstrom_h)
         end if
      end if

      engine%is_setup = .true.
   end subroutine engine_setup

   pure function setup_failed(ierr) result(bad)
      !! `.true.` iff `ierr` is present and non-OK — the "shall I
      !! return early" test used after every ierr-threaded call in
      !! `engine_setup`. When `ierr` is absent the callee already
      !! `error stop`ped on failure, so this is always `.false.` here.
      integer, intent(in), optional :: ierr
      logical :: bad
      bad = .false.
      if (present(ierr)) bad = (ierr /= OCEAN_STATUS_OK)
   end function setup_failed

   subroutine engine_configure_diag(cfg, state, mpi_rank, decomp, ierr)
      !! Register the default diag-manager variable set and open the
      !! per-rank NetCDF stream. No-op when diagnostics are disabled.
      !! Moved verbatim from `driver_run_ocean`'s private
      !! `configure_ocean_diag` helper (P2.4) so the API/bench setup
      !! paths can reach it too; F5 residual — now takes optional
      !! `ierr` instead of `error stop`ping (see `rdb_ocean_diag`'s
      !! `set_output_*_levels` for the leaf-level conversions this
      !! threads through).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: mpi_rank
      type(decomp_t), intent(in) :: decomp
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. cfg%ocean%diag%enabled) return

#ifdef RDB_NO_NETCDF
      call fail("engine_configure_diag: &ocean_diag_nml enabled = .true. but this build "// &
                "has no NetCDF support (RDB_ENABLE_NETCDF=OFF)", ierr, OCEAN_STATUS_ERR_SETUP)
#else
      block
         character(len=512) :: rank_filename
         logical :: multi_rank_diag

         call set_diag_remap_method(parse_diag_remap_scheme(cfg%ocean%diag%diag_remap_scheme))
         call set_diag_mask_vanished(cfg%ocean%diag%mask_vanished_layers)

         if (cfg%ocean%diag%n_z_levels > 0) then
            call state%diag%set_output_z_levels( &
               cfg%ocean%diag%z_levels(1:cfg%ocean%diag%n_z_levels), ierr=ierr)
            if (setup_failed(ierr)) return
         end if
         if (cfg%ocean%diag%n_sigma_levels > 0) then
            call state%diag%set_output_sigma_levels( &
               cfg%ocean%diag%sigma_levels(1:cfg%ocean%diag%n_sigma_levels), ierr=ierr)
         else
            call state%diag%set_output_sigma_levels( &
               uniform_diag_levels(state%multilayer%nz_ml), ierr=ierr)
         end if
         if (setup_failed(ierr)) return
         if (cfg%ocean%diag%n_zstar_levels > 0) then
            call state%diag%set_output_zstar_levels( &
               cfg%ocean%diag%zstar_levels(1:cfg%ocean%diag%n_zstar_levels), ierr=ierr)
         else
            call state%diag%set_output_zstar_levels( &
               uniform_diag_levels(state%multilayer%nz_ml), ierr=ierr)
         end if
         if (setup_failed(ierr)) return
         if (cfg%ocean%diag%n_rho_levels > 0) then
            call state%diag%set_output_density_levels( &
               cfg%ocean%diag%rho_levels(1:cfg%ocean%diag%n_rho_levels), ierr=ierr)
            if (setup_failed(ierr)) return
         end if
         call apply_diag_selection(state, trim(cfg%ocean%diag%diags), &
                                   dt_out=cfg%ocean%diag%dt_out, &
                                   default_coord=diag_vgrid_from_name(cfg%ocean%diag%vgrid))

         ! P7 F1: a fresh checkout's default `output_dir = "./output"`
         ! (rdb_config.F90:2316) does not exist until something creates
         ! it, and `nc_create_file` used to `error stop` the whole host
         ! process on the resulting ENOENT — the single most likely
         ! real-world `create()` failure for a first-time caller.
         ! Best-effort `mkdir`; a still-missing directory (e.g. a missing
         ! GRANDPARENT, or an unwritable parent) surfaces as the ordinary
         ! `open_stream` ierr/error-stop below, naming the exact path via
         ! the underlying `nf90_create` failure message.
         call ensure_directory_exists(trim(cfg%output_dir))

         rank_filename = output_rank_filename(trim(cfg%output_dir), &
                                              trim(cfg%ocean%diag%filename), mpi_rank)
         multi_rank_diag = decomp%px*decomp%py > 1
         if (multi_rank_diag) then
            if (cfg%compress_output) then
               call open_stream(state%diag, trim(rank_filename), &
                                deflate_level=cfg%compress_level, &
                                px=decomp%px, py=decomp%py, &
                                i_start=decomp%i_start, j_start=decomp%j_start, &
                                nx_global=decomp%nx_global, ny_global=decomp%ny_global, &
                                nx_local=decomp%nx_local, ny_local=decomp%ny_local, &
                                nghost=cfg%nghost, &
                                output_precision=trim(cfg%ocean%diag%output_precision), &
                                ierr=ierr)
            else
               call open_stream(state%diag, trim(rank_filename), &
                                px=decomp%px, py=decomp%py, &
                                i_start=decomp%i_start, j_start=decomp%j_start, &
                                nx_global=decomp%nx_global, ny_global=decomp%ny_global, &
                                nx_local=decomp%nx_local, ny_local=decomp%ny_local, &
                                nghost=cfg%nghost, &
                                output_precision=trim(cfg%ocean%diag%output_precision), &
                                ierr=ierr)
            end if
         else
            if (cfg%compress_output) then
               call open_stream(state%diag, trim(rank_filename), &
                                deflate_level=cfg%compress_level, &
                                output_precision=trim(cfg%ocean%diag%output_precision), &
                                ierr=ierr)
            else
               call open_stream(state%diag, trim(rank_filename), &
                                output_precision=trim(cfg%ocean%diag%output_precision), &
                                ierr=ierr)
            end if
         end if
         if (setup_failed(ierr)) return
      end block
#endif
   end subroutine engine_configure_diag

   pure function diag_vgrid_from_name(name) result(vgrid)
      !! Map the `&ocean_diag_nml vgrid` string to a `DIAG_VGRID_*` tag.
      !! Unrecognised (schema-validated upstream) falls back to LAYER.
      !! Moved from `rdb_driver` (P2.4) so `engine_configure_diag` can
      !! use it without a driver dependency; re-exported by
      !! `rdb_driver` for any existing caller of that name.
      character(len=*), intent(in) :: name
      integer :: vgrid
      select case (trim(name))
      case ("z_fixed")
         vgrid = DIAG_VGRID_Z_FIXED
      case ("sigma")
         vgrid = DIAG_VGRID_SIGMA
      case ("zstar")
         vgrid = DIAG_VGRID_ZSTAR
      case ("density")
         vgrid = DIAG_VGRID_DENSITY
      case default
         vgrid = DIAG_VGRID_LAYER
      end select
   end function diag_vgrid_from_name

   pure function uniform_diag_levels(nz) result(levels)
      !! `nz` uniform cumulative fractions (m/nz, m=1..nz) — the auto
      !! sigma / z* output grid when no explicit levels are configured.
      !! Moved from `rdb_driver` (P2.4).
      integer, intent(in) :: nz
      real(wp) :: levels(max(nz, 1))
      integer :: m, n
      n = max(nz, 1)
      do m = 1, n
         levels(m) = real(m, wp)/real(n, wp)
      end do
   end function uniform_diag_levels

   ! ================================================================
   ! Device placement
   ! ================================================================

   subroutine engine_enter_data(engine, cfg)
      !! Map the ocean state onto the device, then (when `bt_halo > 0`)
      !! allocate + attach the wide-halo BT march-in workspace, then
      !! warm up the device-path halo exchanges (UCX cuda_ipc handle
      !! open, outside any timed region — see `driver_run_ocean`'s
      !! "Content-safety argument" comment, carried over unchanged).
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(in) :: cfg

      call ocean_state_enter_data(engine%state)
      engine%state%dyn%bt_halo = cfg%ocean%bt%bt_halo
      if (cfg%ocean%bt%bt_halo > 0) then
         call ocean_dyn_enable_bt_wide(engine%state%dyn, engine%grid, &
                                       cfg%dx, cfg%dy, &
                                       cfg%ocean%grid%lon_west, &
                                       cfg%ocean%grid%lat_south, &
                                       cfg%ocean%grid%rad_earth, &
                                       cfg%ocean%grid%grid_config, &
                                       cfg%coriolis_f, &
                                       cfg%ocean%topo%coriolis_beta, &
                                       cfg%ocean%topo%coriolis_y_ref, &
                                       cfg%ocean%grid%coriolis_scheme)
      end if

      call ocean_halo_exchange_ml_state(engine%state%multilayer)
      if (cfg%ocean%bt%n_inner >= 1) then
         call ocean_halo_bt_group_2d(engine%state%dyn%bt_work%bt_eta, &
                                     engine%state%dyn%bt_work%bt_ubt, &
                                     engine%state%dyn%bt_work%bt_vbt)
         call ocean_halo_face_x(engine%state%dyn%bt_work%bt_ubt)
         if (cfg%ocean%bt%bt_halo > 0) then
            call ocean_halo_bt_group_2d_wide( &
               engine%state%dyn%bt_wide%w_eta, &
               engine%state%dyn%bt_wide%w_ubt, &
               engine%state%dyn%bt_wide%w_vbt, &
               engine%state%dyn%bt_wide%ng_wide)
         end if
      end if

      engine%device_mapped = .true.
   end subroutine engine_enter_data

   ! ================================================================
   ! One outer step
   ! ================================================================

   subroutine engine_step(engine, dt, t, ierr)
      !! Advance one outer step at fixed `dt`, starting from simulation
      !! time `t` (consumed by tide/astro forcing inside
      !! `ocean_dyn_step_split`; the caller advances its own
      !! `t_current` by `dt` afterwards — see the module docstring for
      !! what stays caller-side). Wraps `ocean_dyn_step_split`
      !! (`rdb_ocean_dyn.F90:1986`) plus the per-step halves of the
      !! setup stages that must run BEFORE the dyn-core advance:
      !! file-forcing update/apply, boundary-data refresh, porous-area
      !! refresh.
      !!
      !! Split from `engine_step_finalize` (surface-flux component
      !! assembly + the diag-manager step) because
      !! `driver_run_ocean`'s sea-ice block sits, in the driver's own
      !! per-step sequence, BETWEEN the dyn-core advance and
      !! `ocean_surface_flux_assemble` — the assembler reads the
      !! ice-written salt/heat components, so it must run AFTER the ice
      !! block closes (the driver's own "PR-12 ... MUST sit here"
      !! comment). Since `engine_step` does not itself run the
      !! sea-ice per-step physics, a caller calls `engine_step`, then
      !! `engine_step_ice` (P2.4b — no-op when ice is disabled), THEN
      !! `engine_step_finalize` — exactly `driver_run_ocean`'s order
      !! (both `driver_run_ocean` and the C ABI's `rdb_ocean_step`
      !! do this now).
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: t
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK

#ifndef RDB_NO_NETCDF
      call ocean_data_input_update_all(engine%state%data_input, t)
      call ocean_data_forcing_apply(engine%state%data_forcing, engine%state%data_input, &
                                    engine%grid, engine%state%surface_stress, &
                                    engine%state%surface_flux, engine%state%bc, t)
#endif

      call engine%bc_source%update(t, engine%state%bc)

      call ocean_porous_refresh(engine%grid, engine%state%metrics, engine%state%multilayer)

      ! PR-23b: rebuild the ANALYTIC sponge target on the live layer
      ! geometry, once per outer step, beside the porous refresh and
      ! before the dyn step -- so the whole step relaxes toward the
      ! geopotential profile that was asked for rather than toward the
      ! t = 0 layer positions.  No-op unless `&ocean_sponge_nml
      ! target_source = "linear_z"`, so every other configuration is
      ! bit-identical.
      call ocean_sponge_refresh_target(engine%grid, engine%state%sponge, &
                                       engine%state%multilayer)

      if (engine%n_inner >= 1) then
         call ocean_dyn_step_split(engine%grid, engine%state%metrics, engine%state%dyn, &
                                   engine%state%eos, &
                                   engine%state%coriolis_adv, engine%state%continuity, &
                                   engine%state%pressure_force, engine%state%hvisc, &
                                   engine%state%bdrag, engine%state%surface_stress, &
                                   engine%state%vert_advect, engine%state%hdiff_tracer, &
                                   engine%state%vdiff, engine%state%vmix, &
                                   engine%state%multilayer, dt, engine%n_inner, &
                                   sf=engine%state%surface_flux, &
                                   geo=engine%geo, &
                                   vcoord=engine%state%vcoord, t=t, &
                                   bc=engine%state%bc, sp=engine%state%sponge, &
                                   lateral_mix=engine%state%lateral_mix, &
                                   epbl=engine%state%epbl, kshear=engine%state%kshear, &
                                   mle=engine%state%mle, slopes=engine%state%slopes, &
                                   gm=engine%state%gm, varmix=engine%state%varmix, &
                                   wavespeed=engine%state%wavespeed, &
                                   redi=engine%state%redi, meke=engine%state%meke, &
                                   vmix_tidal=engine%state%vmix_tidal, &
                                   tides=engine%state%tides, &
                                   psurf=engine%state%p_surf, &
                                   td=engine%state%tdrag, &
                                   cav=engine%state%cavity_flux)
      else
         call profiler_start("ocean_dyn_step")
         call ocean_dyn_step(engine%grid, engine%state%metrics, engine%state%dyn, engine%state%eos, &
                             engine%state%coriolis_adv, engine%state%continuity, &
                             engine%state%pressure_force, engine%state%hvisc, &
                             engine%state%bdrag, engine%state%surface_stress, &
                             engine%state%vert_advect, engine%state%hdiff_tracer, &
                             engine%state%vdiff, engine%state%vmix, &
                             engine%state%multilayer, dt, sf=engine%state%surface_flux, &
                             geo=engine%geo, &
                             lateral_mix=engine%state%lateral_mix, &
                             epbl=engine%state%epbl, kshear=engine%state%kshear, &
                             slopes=engine%state%slopes, &
                             vmix_tidal=engine%state%vmix_tidal, &
                             bc=engine%state%bc, t=t, td=engine%state%tdrag, &
                             cav=engine%state%cavity_flux)
         call profiler_stop("ocean_dyn_step")
      end if
   end subroutine engine_step

   subroutine engine_step_ice(engine, cfg, dt, t, ierr)
      !! P2.4b: sea-ice per-step physics — ocean-side frazil
      !! accumulation, EVP dynamics (every outer step, independent of
      !! the thermo cadence: the ice's own fast/slow split) plus the
      !! resulting ice->ocean stress coupling, and, at thermo cadence
      !! (`engine%state%dyn%is_thermo_step()`), category transport
      !! followed by the atmospheric-forcing / basal-flux /
      !! frazil-uptake / column-thermo / snowfall / brine / heat /
      !! shortwave / ITD chain. Transcribed VERBATIM from
      !! `driver_run_ocean`'s former inline sea-ice block — the
      !! internal call order is load-bearing (see the "MANDATED ORDER"
      !! comment below, itself carried over unchanged) and is NOT
      !! reordered here. Returns immediately, a no-op, when
      !! `&ocean_ice_nml enable = .false.` — bit-identical to before
      !! this phase.
      !!
      !! Call between `engine_step` and `engine_step_finalize` — the
      !! surface-flux assembler the finalize call runs reads the
      !! salt/heat components this routine writes (the driver's own
      !! "PR-12 ... MUST sit here" comment), so it must run AFTER this
      !! returns. `cfg` is threaded through (unlike `engine_step`)
      !! because the thermo-forcing knobs this block reads every
      !! thermo-cadence step (`air_temp`/`restore_lambda`/`sw_down`/
      !! `snowfall`/`transport`/`adv_substeps`/`roll_factor`) are read
      !! straight from `&ocean_ice_nml` in the driver too, not cached
      !! at setup (unlike `evp_params`/`ic_par`, which the config
      !! builds once, up front).
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: t
         !! Unused today (the sea-ice block has no direct time
         !! dependence — only `dt` and the `therm_dt`/thermo-cadence
         !! machinery derived from it); present for signature symmetry
         !! with `engine_step`/`engine_step_finalize`.
      integer, intent(out), optional :: ierr
         !! Present for signature symmetry; no failure path in this
         !! routine sets it today — `ice_transport_step`'s
         !! conservation/positivity check still `error stop`s,
         !! matching the pre-existing driver behaviour verbatim.

      logical :: ice_ok
      integer :: ice_n_trunc

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. engine%state%ice%enable) return

      ! Sea-ice PR 1/3b/3c: ocean-side frazil accumulation, then (at thermo
      ! cadence) the atmospheric-forcing seam + basal flux + column
      ! thermodynamics + frazil-bank spend, mediated back to the ocean via
      ! the salt/heat couplers.  ice_frazil_accumulate runs on the FINAL
      ! (post-RK2-average, post-ALE-remap) tracer state of the outer step —
      ! clamping inside the RK2 stages would not bound the averaged result —
      ! so T(k=nz) >= T_f(S) holds at every step boundary the ice model
      ! observes.
      call profiler_start("ice_frazil")
      call ice_frazil_accumulate(engine%grid, engine%state%eos, &
                                 engine%state%multilayer, &
                                 engine%state%ice%frazil_heat, &
                                 engine%state%ice%heat_budget_frazil)
      call profiler_stop("ice_frazil")

      ! Sea-ice PR 5: C-grid EVP dynamics, every outer step (NOT
      ! thermo-cadence gated — EVP is the ice's own fast/slow split,
      ! analogous to the ocean's barotropic/baroclinic split, and runs
      ! every outer step regardless of the thermo cadence).  ice_evp_step
      ! writes ice%u_ice/v_ice directly (replacing the PR-4b transport
      ! sampler, gated off below); ice_ocean_stress_flux immediately
      ! mediates the resulting drag into the ocean's surface-stress field so
      ! the NEXT ocean momentum step feels it (one-step-lagged, same
      ! convention as the frazil/heat couplers).
      if (engine%state%ice%dynamics) then
         call profiler_start("ice_evp")
         ! PR 36: the CFL bound must use the dt TRANSPORT will actually
         ! consume, not this call's outer `dt` -- EVP runs every outer
         ! step, transport at thermo cadence (`engine%state%dyn%therm_dt(dt)`,
         ! the same bound function the transport call below already uses).
         call ice_evp_step(engine%grid, engine%state%metrics, &
                           engine%state%coriolis_adv%f_corner, &
                           engine%state%ice, engine%state%multilayer, dt, &
                           engine%evp_params, &
                           engine%state%bc%periodic_x, engine%state%bc%periodic_y, &
                           dt_transport=engine%state%dyn%therm_dt(dt), &
                           n_trunc=ice_n_trunc)
         if (ice_n_trunc > 0 .and. engine%decomp%rx == 0 .and. engine%decomp%ry == 0) then
            call logger%warning("ice_evp_step: ice velocity CFL-truncated at N faces "// &
                                "(&ocean_ice_nml cfl_trunc); the ice dynamics is "// &
                                "unstable or the transport step is too long")
         end if
         call ice_ocean_stress_flux(engine%state%metrics, engine%state%surface_stress, &
                                    engine%state%ice)
         call profiler_stop("ice_evp")
      end if

      ! Sea-ice PR 3b/3c: at thermo cadence, drive the column live and spend
      ! the bank.  Post-step, outer_step_count has already incremented, so
      ! is_thermo_step() fires at the END of each thermo window — the rates
      ! written here are integrated by exactly ONE apply (therm_dt) in the
      ! NEXT window (MEKE/frazil one-step-lag convention).  Gated on
      ! enable_thermodynamics too: with thermo off the apply path never
      ! fires, and spending the bank / stepping the column would strand
      ! un-mediated salt/heat.
      !
      ! MANDATED ORDER (PLAN_ICE_PR3c "Driver wiring" + the ordering caution
      ! therein; PR 4a appends the ITD restore; PR 26 inserts the snowfall
      ! ocean-share contributor; PR 31 inserts the shortwave coupler):
      ! forcing -> basal -> frazil uptake -> column driver -> snowfall ocean
      ! share -> brine coupler -> heat coupler -> shortwave coupler ->
      ! adjust categories (ITD restore; ncat=1 short-circuits inside).
      ! Frazil uptake MUST precede the column driver: `ice_frazil_uptake_impl`
      ! OVERWRITES `salt_flux_diag` (unconditional zero, then a gated write),
      ! so if the column ran first its net-melt salt contribution would be
      ! clobbered.  Running frazil first (fresh window) then having the
      ! column driver ADD its net-melt term on top composes correctly — a
      ! cell with no frazil but melting ice starts at 0 and the column adds
      ! its (negative) contribution.  `ice_snowfall_ocean_share` (PR 26,
      ! gated on `has_snowfall`) MUST run AFTER the column driver (which
      ! unconditionally zeroes `heat_flux_diag` — running snowfall first
      ! would have its contribution clobbered) and BEFORE the brine/heat
      ! couplers (which overwrite `Q_salt`/`Q_heat` from
      ! `salt_flux_diag`/`heat_flux_diag` — running it after would never
      ! reach the ocean).  It ADDS to both diags, same contract as the
      ! column driver's net-melt term.  The couplers run next, after every
      ! contributor has written salt_flux_diag/heat_flux_diag for this
      ! window.  `ice_ocean_sw_flux` (PR 31) runs immediately AFTER
      ! `ice_ocean_heat_flux`: in the components-off default it ADDS
      ! `sw_thru_diag` onto the `Q_heat` the heat coupler just
      ! full-overwrote, so it must observe that overwrite first (the
      ! shortwave is deliberately NOT in heat_flux_diag — see
      ! `ice_ocean_sw_flux`, no double count in either component mode).
      ! `ice_adjust_categories` runs LAST: it only reshuffles category
      ! area/mass/enthalpy/salt (no thermodynamics, no diag writes) after
      ! every thermodynamic thickness change this window, so the couplers
      ! (which only READ the per-cell diags written above) are unaffected
      ! by it running after them.
      if (engine%state%dyn%enable_thermodynamics .and. &
          engine%state%dyn%is_thermo_step()) then
         ! Sea-ice PR 4b: category ice/snow transport + compress_ice,
         ! BEFORE the thermo forcing chain (SIS2 slow sequence:
         ! dynamics+transport, then slow thermo).  Cadence = the thermo step
         ! (the ice slow step); gated on `&ocean_ice_nml transport` (default
         ! off => byte-identical).
         if (cfg%ocean%ice%transport) then
            call profiler_start("ice_transport")
            call ice_transport_step(engine%grid, engine%state%metrics, &
                                    engine%state%multilayer, &
                                    engine%state%ice, engine%state%dyn%therm_dt(dt), &
                                    cfg%ocean%ice%adv_substeps, cfg%ocean%ice%roll_factor, &
                                    ice_ok)
            call profiler_stop("ice_transport")
            if (.not. ice_ok) then
               call logger%error("ice_transport_step: conservation/positivity "// &
                                 "violation (negative mass, orphan snow, or a "// &
                                 "compress-time consistency failure)")
               error stop "ice_transport_step: conservation/positivity violation"
            end if
         end if
         call profiler_start("ice_thermo")
         call ice_atm_forcing_restoring(engine%state%ice, &
                                        cfg%ocean%ice%air_temp, &
                                        cfg%ocean%ice%restore_lambda, &
                                        cfg%ocean%ice%sw_down, &
                                        cfg%ocean%ice%snowfall)
         call ice_compute_basal_flux(engine%grid, engine%state%eos, &
                                     engine%state%multilayer, engine%state%ice, &
                                     engine%state%dyn%therm_dt(dt))
         call ice_frazil_uptake(engine%grid, engine%state%eos, engine%state%multilayer, &
                                engine%state%ice, engine%state%dyn%therm_dt(dt))
         call ice_thermo_driver_step(engine%grid, engine%state%eos, &
                                     engine%state%multilayer, engine%state%ice, &
                                     engine%state%dyn%therm_dt(dt))
         if (engine%state%ice%has_snowfall) then
            call ice_snowfall_ocean_share(engine%grid, engine%state%ice, &
                                          engine%state%dyn%therm_dt(dt))
         end if
         call ice_ocean_brine_flux(engine%state%surface_flux, engine%state%ice)
         call ice_ocean_heat_flux(engine%state%surface_flux, engine%state%ice)
         call ice_ocean_sw_flux(engine%state%surface_flux, engine%state%ice)
         call ice_adjust_categories(engine%grid, engine%state%multilayer, engine%state%ice)
         call profiler_stop("ice_thermo")
      end if
   end subroutine engine_step_ice

   subroutine engine_step_finalize(engine, dt, t, ierr)
      !! Second half of one outer step: derive `Q_heat`/`Q_salt` from
      !! the surface-flux component set (no-op unless
      !! `&ocean_forcing_nml enable_components`) and, when diagnostics
      !! are configured, run one `ocean_diag_t%step`. Call immediately
      !! after `engine_step` and, when sea ice is enabled,
      !! `engine_step_ice` too — see `engine_step`'s docstring for why
      !! this is a separate call.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: t
         !! Simulation time at the START of the step this finalizes
         !! (the same `t` passed to the matching `engine_step` call) —
         !! `t + dt` is what reaches the diag manager, matching
         !! `driver_run_ocean`'s post-advance `t_current`.
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      ! Ice-shelf basal melt (P2b): solve the three-equation interface on
      ! every covered column and fill the OWNED heat_cavity/salt_cavity
      ! components.  MUST precede the assembler, which folds them into
      ! Q_heat/Q_salt.  Same thermo cadence, and a no-op (immediate
      ! return) when &ocean_cavity_melt_nml enable=.false.  Cavity x sea
      ! ice is refused at configure, so the ordering against
      ! engine_step_ice's fillers is not a live question.
      call ocean_cavity_flux_step(engine%grid, engine%state%cavity_flux, &
                                  engine%state%metrics, engine%state%multilayer, &
                                  engine%state%eos, engine%state%surface_flux, &
                                  active=engine%state%dyn%enable_thermodynamics &
                                  .and. engine%state%dyn%is_thermo_step())

      ! Phase 4b — the MELT-ONLY fallback for the boundary-layer `u_*`.
      !
      ! When `&ocean_tdrag_nml` is on, the RK2 stage drivers publish
      ! `ss%stress_shelf` from the top drag's own `stress_top`, in-stage
      ! and unlagged, and this branch stays out of the way.  When it is
      ! OFF but basal melt is on, nothing else would give KPP/EPBL an
      ! under-ice `u_*` at all — they would mix a covered column on the
      ! masked (exactly zero) wind stress.  The melt slot already solved
      ! for a friction velocity with the SAME `C_d` (the one-drag-
      ! coefficient rule: `&ocean_cavity_melt_nml cdrag_top` must equal
      ! `&ocean_tdrag_nml cd`), so re-deriving `|tau_top| = rho_0*u_*^2`
      ! from it is the consistent answer, not a second drag law.
      !
      ! HONEST LIMIT: `cav%ustar` is refreshed at the THERMO cadence, at
      ! the END of the outer step, so this path reaches the boundary-layer
      ! schemes ONE OUTER STEP LATE.  The top-drag path has no such lag.
      ! `cav%ustar` is exactly 0 on every uncovered column, so the
      ! published field keeps `stress_shelf`'s "zero off the cover"
      ! invariant.
      !
      ! A CALL, not an inline `do concurrent`, and the reason is a GPU
      ! rule rather than a style preference.  A `do concurrent` written
      ! here would reference `engine%state%surface_stress%stress_shelf`,
      ! i.e. it would walk the ENGINE, and `ocean_engine_t` is not a
      ! mapped object -- only `engine%state` is.  nvfortran then emits a
      ! data clause for the whole `engine` and aborts at run time with
      ! "variable in data clause is partially present on the device:
      ! name=engine".  Measured, cc70, 2026-09-20: it did exactly that.
      ! Host-dereferencing the two component arrays AT the call site and
      ! handing them to a flat explicit-shape kernel is the fix -- the
      ! outer-shim + flat-impl pattern.
      !
      ! CLAUDE.md's "write a host-gated pass INLINE" rule does not apply:
      ! that rule exists because an escaping state array pessimises the
      ! OTHER `do concurrent` loops in the calling routine, and
      ! `engine_step_finalize` has none -- it is a four-call orchestrator.
      !
      ! Placeholder safety is STRUCTURAL, not a runtime branch: the gate
      ! is `cavity_flux%enable`, and the melt slot allocates `ustar` at
      ! `(nx, ny)` exactly when that is set (`(1,1)` otherwise), while
      ! `stress_shelf` is always full size.  So the explicit-shape dummies
      ! below can only ever be reached with matching, full-size actuals.
      if (engine%state%cavity_flux%enable .and. .not. engine%state%tdrag%enable) then
         call ocean_surface_stress_set_shelf_from_ustar( &
            engine%state%surface_stress%stress_shelf, &
            engine%state%cavity_flux%ustar, &
            engine%state%surface_stress%rho0, &
            size(engine%state%surface_stress%stress_shelf, 1), &
            size(engine%state%surface_stress%stress_shelf, 2))
      end if

      ! Ice-shelf cover: the assembler is where the atmospheric bands and
      ! the cavity's own heat_cavity/salt_cavity are still separable, so
      ! it is where `1 - cover_frac` is applied — and applying it there
      ! (rather than at apply time) is also what makes the `Q_heat` /
      ! `Q_salt` that KPP and EPBL read for `B_0` the MASKED values.
      ! Cavity off ⇒ the original call, byte-identical.
      if (engine%state%metrics%use_cavity) then
         call ocean_surface_flux_assemble(engine%grid, engine%state%surface_flux, &
                                          engine%state%multilayer, &
                                          active=engine%state%dyn%enable_thermodynamics &
                                          .and. engine%state%dyn%is_thermo_step(), &
                                          cover_frac=engine%state%metrics%cover_frac)
      else
         call ocean_surface_flux_assemble(engine%grid, engine%state%surface_flux, &
                                          engine%state%multilayer, &
                                          active=engine%state%dyn%enable_thermodynamics &
                                          .and. engine%state%dyn%is_thermo_step())
      end if

      if (engine%diag_enabled) then
         call engine%state%diag%step(engine%state, dt, t + dt)
      end if
   end subroutine engine_step_finalize

   ! ================================================================
   ! Teardown
   ! ================================================================

   subroutine engine_exit_data(engine)
      !! Unwind device residency. Idempotent-adjacent: caller checks
      !! `device_mapped` (mirrors `rdb_handle`'s `ocean_handle_t`).
      type(ocean_engine_t), intent(inout) :: engine

      call ocean_state_exit_data(engine%state)
      call ice_ocean_stress_cleanup()
      engine%device_mapped = .false.
   end subroutine engine_exit_data

   subroutine engine_teardown(engine)
      !! Close the diag NetCDF stream (if one was opened), release the
      !! process-global ocean-halo module state, then release host-side
      !! allocations. Call after `engine_exit_data`.
      type(ocean_engine_t), intent(inout) :: engine

#ifndef RDB_NO_NETCDF
      if (engine%diag_enabled) call close_stream(engine%state%diag)
#endif
      call ocean_halo_destroy()
      call engine%geo%destroy()
      call engine%state%destroy()
      engine%is_setup = .false.
   end subroutine engine_teardown

end module rdb_ocean_engine
