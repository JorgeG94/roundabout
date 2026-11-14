# Sources of the standalone ice core library `rdb_ice_core_objs` (OBJECT).
# Infrastructure-only dependencies (pic, safe-math includes): nothing here may
# `use` the ocean state, engine, config or halo.  One path per line,
# alphabetical.
set(RDB_ICE_CORE_SOURCES
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_aggregate.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_atmo.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_basal_heat.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_bl99.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_dims.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_errors.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_flux.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_itd_params.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_itd_shift.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_lateral_melt.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_ledger_sums.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_params.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_restart_convert.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_ridge.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_ridge_itd.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_seed.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_shortwave.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_thermo_shared.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_thickness.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/cice/rdb_ice_vprofile.F90
    ${PROJECT_SOURCE_DIR}/src/core/ice/thermo/rdb_ice_enthalpy.F90)
