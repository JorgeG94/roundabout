# Core infrastructure sources
set(CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_constants.F90
    # Leaf module (no deps): named ierr status codes for the ocean
    # solver-creation path (Python API plan P0). Used by rdb_eos,
    # rdb_ocean_boundary_types, rdb_config, rdb_ocean_metrics, rdb_ocean_setup,
    # rdb_ocean_state, rdb_ocean_z_init.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_status.F90
    # Leaf module (only deps: pic_logger): ring buffer of specific error
    # messages behind the coarse ierr/status codes above (Python API plan P2).
    # Used by the same create()-path modules as rdb_ocean_status.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_error_ring.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_console_stats.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_nml_schema.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_array_utils.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ALE/rdb_vcoord.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_decomp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_config.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_config_schema.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_grid.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_tracer.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rdb_scratch_3d.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rdb_safe_math.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rdb_efp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rdb_mem_report.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_barotropic_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_barotropic_workstate.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_multilayer_state.F90
    # Ocean-regime god state + per-kernel state types.  All Phase 0d empty
    # shells; declared here (not under solver/kernel sources) because rdb_state
    # composes ocean_state_t.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/dynamics/split_rk2/rdb_ocean_bt_wide.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/dynamics/split_rk2/rdb_ocean_dyn.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/barotropic/rdb_bt_cont_type.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/barotropic/rdb_barotropic_substep.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/barotropic/rdb_barotropic_coupling.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/continuity_ppm/rdb_continuity.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/continuity_ppm/rdb_ocean_min_thickness.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/coriolis_adv/rdb_coriolis_adv.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/kernels/vertical_advection/rdb_ocean_vertical_advection.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/pressure_force/rdb_ocean_pgf_reconstruct.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/pressure_force/rdb_ocean_pressure_force.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_vmix.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_epbl.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_wave_speed.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_kappa_shear.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_tidal_mixing.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_bottom_drag.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_surface_stress.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_surface_flux.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_geothermal.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_cavity_melt.F90
    # P2b ice-shelf basal-melt COUPLING: the slot, the far-field sampler (metres
    # below the ice base, not "layer nz") and the per-thermo-step driver that
    # fills the owned heat_cavity/salt_cavity components.
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_cavity_flux.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rdb_ocean_vdiff.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_lateral_mix.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_isopycnal_slopes.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_mle.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_gm.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_varmix.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_meke.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_redi.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_horizontal_viscosity.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/rdb_ocean_hdiff_tracer.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/equation_of_state/rdb_eos.F90
    # Ocean-side EOS shim over the C-grid multilayer state — keeps rdb_eos
    # itself free of ocean-state dependencies; depends on rdb_eos +
    # rdb_multilayer_state (both declared above).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_eos_compute.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/forcing/rdb_ocean_tide_astro.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/forcing/rdb_ocean_tides.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/forcing/rdb_ocean_p_surf.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_obc.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_boundary_types.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_boundary_data.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_sponge.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_periodic.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_fold.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_fold_apply.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/rdb_ocean_halo_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/rdb_ocean_ghost_poison.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/boundary/rdb_ocean_obc_baroclinic.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_diag_mask.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_diag.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_budgets.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_bt_budget_probe.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_chksum.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_ke_probe.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/io/rdb_ocean_restart.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/vcoord/rdb_ocean_vcoord.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_massless.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_bipolar.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_metrics.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_porous.F90
    # P5.1 static ice-shelf cavity geometry: the prescribed draft setters, the
    # grounding/cross-validation helpers and the isostatic load fill. Leaf
    # module (rdb_constants + rdb_grid only) — deliberately free of rdb_config
    # so the formula dispatch, which needs the namelist, stays in the state
    # builder next to the bathymetry dispatch it mirrors.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_cavity.F90
    # P2.5 pre-create geometry injection (Python runtime API plan): sign
    # normalisation/validation + ghost fill for an in-memory bathymetry array.
    # Deliberately NetCDF-free (unlike rdb_bathymetry, which requires
    # RDB_ENABLE_NETCDF=ON purely because it also contains the NetCDF reader) —
    # array injection must work with no filesystem touch.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_bathymetry_inject.F90
    # P2.5: one required_halo() query folding the six scattered minimum- nghost
    # rules. Depends on rdb_coriolis_adv (CORE_SOURCES, above) and
    # rdb_recon_weno (KERNEL_SOURCES, below) — CMake's Fortran module dependency
    # scanning orders the actual compilation, not this list.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_halo_width.F90
    # Sea-ice PR 2: pure device-callable enthalpy library (T<->E inversions);
    # leaf module. Precedes rdb_ice_optics/_mass/_column (PR 3a), which all
    # depend on it.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_enthalpy.F90
    # Sea-ice PR 3a: shortwave optics (CSIM4 branch) + column mass bookkeeping
    # (bottom-freeze, melt peel, rebalance) — both leaf modules `rdb_ice_column`
    # depends on.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_optics.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_mass.F90
    # Sea-ice PR 3a: single-column Winton NkIce-layer thermodynamics (SEB +
    # conduction solve, per-column orchestrator, do-concurrent driver).
    # Exercised only by the new test suite — nothing in the driver calls it yet.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_column.F90
    # Sea-ice slot scaffold (SIS2 port, PR 0) — composed by rdb_ocean_state, so
    # it must precede it in the list. PR 3a: now depends on rdb_ice_column for
    # ICE_BULK_SALINITY/ICE_NK_MAX (the Winton column prognostics' init values).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/state/rdb_ice_state.F90
    # Sea-ice PR 24: analytic initial-condition path (uniform/latitudes
    # concentration seeding). Host-side setup code, no device kernel; depends on
    # rdb_ice_state, rdb_ice_column (ICE_RHO_ICE/ICE_RHO_SNOW/
    # ICE_BULK_SALINITY), rdb_ice_enthalpy (ice_enth_from_ts), rdb_grid,
    # rdb_multilayer_state, rdb_ocean_metrics (all declared above).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/state/rdb_ice_init.F90
    # Sea-ice PR 4a: multi-category ITD restore (adjust_ice_categories); depends
    # on rdb_ice_state, rdb_ice_column (ICE_RHO_ICE et al, declared above),
    # rdb_constants, rdb_grid, rdb_multilayer_state.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/itd/rdb_ice_itd.F90
    # Sea-ice PR 4b: category ice/snow transport + compress_ice; depends on
    # rdb_ice_state, rdb_ice_itd (ice_adjust_categories, Phase 5),
    # rdb_continuity (the 5 promoted PPM helpers, declared earlier in this file
    # at the ocean/kernels/continuity_ppm entry), rdb_multilayer_state,
    # rdb_ocean_metrics.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/transport/rdb_ice_transport.F90
    # Sea-ice PR 5: C-grid EVP ice dynamics; depends on rdb_ice_state
    # (ice_cell_concentration_impl), rdb_ice_column (ICE_RHO_ICE),
    # rdb_ocean_metrics, rdb_multilayer_state, rdb_ocean_periodic (the
    # face/centre ghost-wrap helpers, declared earlier in this file at the
    # ocean/boundary entry).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/dynamics/rdb_ice_evp.F90
    # Sea-ice PR 1: ocean-side frazil accumulator (freezing-point clamp), called
    # from the driver after each ocean outer step.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_frazil.F90
    # Sea-ice PR 3c: coupleable atmospheric-forcing seam filler (v1 scalar
    # restoring); depends only on rdb_ice_state (already declared above).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_atm_forcing.F90
    # Sea-ice PR 3c: ocean -> ice basal heat flux + the sst/s_surf/tfw sample
    # seam; depends on rdb_ice_state + rdb_ocean_surface_flux (SEAWATER_CP,
    # declared earlier in this file) + rdb_ice_column (ICE_RHO_ICE, for the
    # ice-presence gate — declared above at rdb_ice_column.F90).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_basal_flux.F90
    # Sea-ice PR 3b: frazil-bank spend (new-ice formation) + the ice->ocean
    # brine-rejection coupler (Q_salt refresh); both called from the driver at
    # thermo cadence, after rdb_ice_frazil in the loop.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_frazil_uptake.F90
    # Sea-ice PR 3c: driver-facing column step + ocean-coupling-diag reduction;
    # depends on rdb_ice_column (Winton column) + rdb_ice_state.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_thermo_driver.F90
    # Sea-ice PR 26: snowfall's ice-free-share ocean contributor (heat +
    # virtual-salt); depends on rdb_ice_state + rdb_ice_enthalpy (declared
    # above). Runs AFTER rdb_ice_thermo_driver, BEFORE rdb_ice_ocean_coupler in
    # the driver's thermo-cadence block (mandated order).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/thermo/rdb_ice_snow.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ice/coupling/rdb_ice_ocean_coupler.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_setup.F90
    # Configure-time numerical-stability audit (viscous CFL, tracer diffusive
    # number, Munk-layer resolution, ah_max/nu_h clamp): depends on rdb_config +
    # rdb_ocean_metrics (both declared above).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/state/rdb_ocean_stability_audit.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_diag_fills.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_diag_derived.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_console_stats.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rdb_profiler.F90)

# The z-level initialisation + data-override readers load NetCDF inputs.  Skip
# them when the I/O subsystem is disabled so the kernels/benchmarks still build.
if(RDB_ENABLE_NETCDF)
  list(
    APPEND
    CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/io/rdb_ocean_z_init.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/io/rdb_ocean_data_input.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/forcing/rdb_ocean_data_forcing.F90)
endif()

# Communication layer — ONE implementation, always compiled.
#
# There is no `single/` stub twin any more. Every call site talks to
# `pic_mpi_lib`, and pic-mpi decides at ITS build time (`PIC_ENABLE_MPI`)
# whether that is a real MPI library or its serial backend. See `comm/README.md`
# for the single-rank contract that makes this safe.
set(COMM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/comm/rdb_comm_env.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/comm/rdb_halo.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/comm/rdb_ocean_halo_counters.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/comm/rdb_ocean_halo.F90)

if(RDB_ENABLE_MPI)
  message(STATUS "Comm backend: pic-mpi over MPI (multi-rank)")
else()
  message(STATUS "Comm backend: pic-mpi serial backend (single rank)")
endif()

# High-level driver — owns the full compute-rank lifecycle.  Dispatches straight
# into the split-RK2 dyn-core (core/ocean/dynamics/split_rk2/); there is no
# separate solver-driver module.  Lives below `app/main.F90`.  Test-callable so
# unit tests can run a full lifecycle and catch integration-level bugs. Pulls in
# output/gauges/forcing/nesting, so it requires NetCDF.
set(DRIVER_SOURCES)
if(RDB_ENABLE_NETCDF)
  list(APPEND DRIVER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/driver/rdb_driver.F90)
endif()

set(KERNEL_SOURCES
    # Structured tracer reconstruction (PLM/WENO face-swept helpers) + the
    # shipped passive-tracer packages.
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/rdb_recon_weno.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/rdb_ocean_ideal_age.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/rdb_ocean_pseudo_salt.F90
    # Column remap kernel + the ocean-side orchestrator that wires it across the
    # tracer registry — both live under ALE/ since they're vcoord machinery, not
    # split-RK2-specific.
    ${CMAKE_CURRENT_SOURCE_DIR}/ALE/rdb_remap_column.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ALE/rdb_ocean_remap.F90)
if(RDB_PARALLEL_BACKEND STREQUAL "openmp")
  set(_kernel_directive_label "OpenMP target")
else()
  set(_kernel_directive_label "OpenACC")
endif()
if(RDB_ENABLE_GPU)
  message(
    STATUS
      "Kernel backend: GPU (do concurrent + ${_kernel_directive_label} offload)"
  )
else()
  message(
    STATUS
      "Kernel backend: CPU (do concurrent serial, ${_kernel_directive_label} directives ignored)"
  )
endif()

# I/O sources (NetCDF)
set(IO_SOURCES)
if(RDB_ENABLE_NETCDF)
  list(
    APPEND
    IO_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rdb_io_netcdf.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rdb_banner.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rdb_bathymetry.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/diag/rdb_ocean_diag_netcdf.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/io/rdb_ocean_restart_io.F90)
endif()

# Ocean setup/step/teardown engine (P2.4) — the ONE shared lifecycle path used
# by driver_run_ocean, rdb_ocean_api and (partially) bench_ocean. ALWAYS
# compiled: API_SOURCES below needs it even when RDB_ENABLE_NETCDF=OFF, so it
# internally #ifndef RDB_NO_NETCDF-guards the NetCDF-only stages (diag stream,
# file-backed forcing, restart filename resolution) rather than being excluded
# from the source list the way DRIVER_SOURCES is. Depends on rdb_ocean_state,
# rdb_ocean_setup, rdb_ocean_dyn, rdb_ice_evp, rdb_ice_init,
# rdb_ice_ocean_coupler, rdb_ocean_diag*, rdb_ocean_sponge, rdb_ocean_geothermal
# (declared above in CORE_SOURCES) plus rdb_ocean_halo/rdb_ocean_boundary_data/
# rdb_ocean_periodic/rdb_ocean_fold_apply/rdb_ocean_dyn (KERNEL_SOURCES/
# COMM_SOURCES, declared above too).
set(ENGINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/ocean/rdb_ocean_engine.F90)

# C ABI (Python runtime API plan, P1) — compiles INTO core_objs, never linked
# against it: RDC cannot cross a shared-library boundary (see the core_objs
# comment block in the root CMakeLists.txt), so the shim has to be part of the
# same object library that backs both the static and shared cores. Depends on
# rdb_config, rdb_state, rdb_ocean_state, rdb_ocean_setup, rdb_ocean_vcoord,
# rdb_ocean_dyn, rdb_ocean_status, rdb_ocean_engine — all declared above in
# CORE_SOURCES/ENGINE_SOURCES — so this list must be appended to target_sources
# after them.
set(API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/api/rdb_handle.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/api/rdb_ocean_api.F90)

target_sources(
  ${core_objs}
  PRIVATE ${CORE_SOURCES}
          ${COMM_SOURCES}
          ${IO_SOURCES}
          ${KERNEL_SOURCES}
          ${ENGINE_SOURCES}
          ${DRIVER_SOURCES}
          ${API_SOURCES})
