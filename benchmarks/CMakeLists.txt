# Benchmark suite for the ocean dyn-core

# Benchmarks link the SHARED core by default (fast link), mirroring the test
# suite (RDB_TESTS_LINK_SHARED, on by default). A self-contained `do concurrent`
# links against the .so fine; benches that CALL a core `!$acc routine seq` from
# their own `do concurrent` (a device→device reference that can't cross a .so
# boundary) were dropped rather than kept on the static core. When the shared
# core wasn't built, `bench_core` falls back to static for everything.
if(TARGET ${main_lib}_shared)
  set(bench_core ${main_lib}_shared)
else()
  set(bench_core ${main_lib})
endif()

# Full ocean dyn-core throughput bench — reads an ocean namelist (default
# double_gyre_mom6.nml) and drives the engine (`engine_setup` / `engine_step`)
# -- the SAME lifecycle the driver runs (rdb_ocean_engine) -- but NetCDF-free,
# so it builds with -DRDB_ENABLE_NETCDF=OFF for Intel/AMD portability +
# profiling.
add_executable(bench_ocean bench_ocean.F90)
target_include_directories(bench_ocean PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                               "${RDB_SHARED_INCLUDE_DIR}")
target_link_libraries(bench_ocean PRIVATE ${bench_core})

# Drift guard: bench_ocean runs the driver's own engine_setup/engine_step
# lifecycle, so a short run of the canonical double-gyre namelist must stay
# finite AND below the velocity clamp (the old hand-copied setup skipped steps
# the driver ran and blew up here, at |u| = maxvel by step ~70). Only exists
# when benchmarks are built (-DRDB_BUILD_BENCHMARKS=ON).
add_test(
  NAME rdb_bench_ocean_smoke
  COMMAND
    bench_ocean
    "${CMAKE_SOURCE_DIR}/validation_examples/ocean/double_gyre/double_gyre_mom6.nml"
    150)
set_tests_properties(rdb_bench_ocean_smoke
                     PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
