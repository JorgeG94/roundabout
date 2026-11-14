cmake_minimum_required(VERSION 3.25)
project(
  rdb
  VERSION 0.1.0
  LANGUAGES Fortran)

set(main_lib core_rdb)
set(exe_name rdb)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# GNUInstallDirs picks lib vs lib64 from CMAKE_SIZEOF_VOID_P. In a Fortran-only
# configure that comes from the Fortran ABI probe
# (CMAKE_Fortran_SIZEOF_DATA_PTR), which some toolchains do not report (seen on
# GH200 / neoverse_v2): an author warning and a `lib` default where the platform
# wants `lib64`. Every supported target is 64-bit, so assume it only when the
# probe is silent; a later enable_language(C) still overrides it.
if(NOT DEFINED CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P 8)
endif()
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Build options + cache variables (16 options, the backend enum check, the
# build-type and GPU-arch defaults).
include(options)

# Compiler-identity guards. AFTER options (it forces some of them off), BEFORE
# compiler_flags.
include(toolchain_guards)

# A shared library needs all code (dependencies included) built -fPIC; so does
# the shared-core test/benchmark path. Reads the two knobs *after* the guards
# above may have forced them off.
if(RDB_BUILD_SHARED OR (RDB_TESTS_LINK_SHARED AND (RDB_ENABLE_TESTING
                                                   OR RDB_BUILD_BENCHMARKS)))
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Preprocessor definitions, collected here and applied to the core target once
# it exists (see target_compile_definitions below).
set(RDB_COMPILE_DEFS)
if(RDB_ENABLE_DOUBLE)
  list(APPEND RDB_COMPILE_DEFS RDB_DOUBLE_PRECISION)
endif()
if(RDB_ENABLE_NVTX)
  list(APPEND RDB_COMPILE_DEFS USE_NVTX)
  message(STATUS "NVTX profiling enabled")
endif()

# Compiler-specific flags (mutates CMAKE_Fortran_FLAGS* and RDB_COMPILE_DEFS)
include(compiler_flags)

include(summary)
rdb_print_build_summary()

# Dependencies: test-drive, pic, NetCDF-Fortran, MPI, pic-mpi (appends to
# RDB_COMPILE_DEFS and sets RDB_MPI_LANG).
include(dependencies)

# Benchmark-only API surface.  Kernels guarded by `#ifdef RDB_BENCHMARK_API`
# expose internal helpers (e.g. `ml_equation_of_state_impl`) only when the
# benchmark suite is being built, so production deployments don't get a wider
# public ABI than they need.
if(RDB_BUILD_BENCHMARKS)
  list(APPEND RDB_COMPILE_DEFS RDB_BENCHMARK_API)
endif()

# Per-column stack-workspace bound (NZ_STACK_MAX in rdb_constants). Plumbed as a
# preprocessor define so it is build-configurable without editing source.
list(APPEND RDB_COMPILE_DEFS RDB_NZ_STACK_MAX=${RDB_NZ_STACK_MAX})

# Core library — built via an OBJECT library so the same compiled objects back
# BOTH an always-static core AND, on request, a shared core, without a second
# compile.
#
# Why the split: tests and benchmarks have their own `do concurrent` kernels
# that call core `!$acc routine seq` device functions.  Relocatable device code
# (RDC) can only be device-linked WITHIN a single executable/archive — it cannot
# cross a shared-library boundary.  So tests/benches must link a STATIC core
# (device code device-links into each exe).  The app has no device kernels of
# its own into the core, so it can link the SHARED core (.so) when
# RDB_BUILD_SHARED=ON.
set(core_objs rdb_ocean_core_objs)
add_library(${core_objs} OBJECT)
set_target_properties(${core_objs} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(${core_objs} PUBLIC ${RDB_COMPILE_DEFS})
target_include_directories(${core_objs} PUBLIC ${CMAKE_BINARY_DIR}/modules)

# Shared include bodies (`src/shared_module_utilities/*.inc`) — small `pure` +
# `!$acc routine seq` helpers `#include`d into the `contains` section of each
# consuming module so NVHPC gets a module-LOCAL copy it can inline into a `do
# concurrent` kernel, with ONE source of truth.  See that directory's README.
# PUBLIC here so it propagates to every core consumer; tests and benchmarks that
# link the STATIC core do NOT inherit include dirs from an OBJECT library
# through the archive, so they add it themselves.
set(RDB_SHARED_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/src/shared_module_utilities)
target_include_directories(${core_objs} PUBLIC ${RDB_SHARED_INCLUDE_DIR})

# Three object-library layers, each compiled ONCE: rdb_core_objs      shared
# infrastructure (rdb_constants, rdb_safe_math, rdb_efp) rdb_ocean_core_objs the
# ocean (${core_objs}, above; legacy sea-ice sources under src/core/ice stay
# here: they use ocean types) rdb_ice_core_objs  the CICE-port ice core
# (cmake/ice_core_sources.cmake); depends only on rdb_core_objs + pic All three
# end up in the same final libraries (below); the `rdb_test_cice_*` tests link
# only the ice and core layers, so they build without the ocean.
add_library(
  rdb_core_objs OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/core/rdb_constants.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/src/framework/rdb_efp.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/src/framework/rdb_safe_math.F90)
set_target_properties(
  rdb_core_objs PROPERTIES POSITION_INDEPENDENT_CODE ON
                           Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_compile_definitions(rdb_core_objs PUBLIC ${RDB_COMPILE_DEFS})
target_include_directories(rdb_core_objs PUBLIC ${CMAKE_BINARY_DIR}/modules
                                                ${RDB_SHARED_INCLUDE_DIR})
target_link_libraries(rdb_core_objs PUBLIC pic::pic)

include(ice_core_sources)
add_library(rdb_ice_core_objs OBJECT ${RDB_ICE_CORE_SOURCES})
set_target_properties(
  rdb_ice_core_objs
  PROPERTIES POSITION_INDEPENDENT_CODE ON Fortran_MODULE_DIRECTORY
                                          ${CMAKE_BINARY_DIR}/modules)
target_compile_definitions(rdb_ice_core_objs PUBLIC ${RDB_COMPILE_DEFS})
target_include_directories(rdb_ice_core_objs PUBLIC ${CMAKE_BINARY_DIR}/modules
                                                    ${RDB_SHARED_INCLUDE_DIR})
target_link_libraries(rdb_ice_core_objs PUBLIC rdb_core_objs pic::pic)
# REAL link edges (not add_dependencies): Ninja's Fortran dyndep scan only sees
# the module directories of targets that are LINKED, so the ocean objects need
# these to find rdb_constants.mod etc. and to recompile when they change.
target_link_libraries(${core_objs} PUBLIC rdb_core_objs rdb_ice_core_objs)

set(libraries_to_link)

# Source directories (populate the OBJECT library)
add_subdirectory(src)

# Module output directory — compilation happens in the OBJECT library
set_target_properties(${core_objs} PROPERTIES Fortran_MODULE_DIRECTORY
                                              ${CMAKE_BINARY_DIR}/modules)

# Link dependencies (PUBLIC on the objects ⇒ propagate to every core consumer)
list(APPEND libraries_to_link pic::pic)
if(RDB_ENABLE_NETCDF)
  list(APPEND libraries_to_link netCDF::netcdff)
endif()
list(APPEND libraries_to_link pic-mpi::pic-mpi)
if(RDB_ENABLE_MPI)
  list(APPEND libraries_to_link MPI::MPI_${RDB_MPI_LANG})
endif()
target_link_libraries(${core_objs} PUBLIC ${libraries_to_link})
if(RDB_ENABLE_NVTX AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  target_compile_options(${core_objs} PUBLIC "-cudalib=nvtx")
  target_link_options(${core_objs} PUBLIC "-cudalib=nvtx")
endif()

# STATIC core — ALWAYS built.  This is what tests + benchmarks link, so their
# device kernels device-link within the executable.
add_library(${main_lib} STATIC)
target_link_libraries(${main_lib} PUBLIC ${core_objs})
set_target_properties(${main_lib} PROPERTIES POSITION_INDEPENDENT_CODE ON
                                             LINKER_LANGUAGE Fortran)
# The infrastructure and ice objects are compiled once and archived here (an
# OBJECT library does not propagate objects added to another OBJECT library).
target_sources(${main_lib} PRIVATE $<TARGET_OBJECTS:rdb_core_objs>
                                   $<TARGET_OBJECTS:rdb_ice_core_objs>)

# SHARED core (.so) — same object library as the static core, so no second
# compile.  Two consumers: the app when RDB_BUILD_SHARED is on (it crosses no
# device boundary into the core), and the test suite + benchmarks when
# RDB_TESTS_LINK_SHARED is on — every test/bench except the few with an
# on-device call into a core `!$acc routine seq` links this instead of the
# static archive, skipping the per-exe static pull + RDC device link.
if(RDB_BUILD_SHARED OR (RDB_TESTS_LINK_SHARED AND (RDB_ENABLE_TESTING
                                                   OR RDB_BUILD_BENCHMARKS)))
  add_library(${main_lib}_shared SHARED)
  target_link_libraries(${main_lib}_shared PUBLIC ${core_objs})
  target_sources(${main_lib}_shared PRIVATE $<TARGET_OBJECTS:rdb_core_objs>
                                            $<TARGET_OBJECTS:rdb_ice_core_objs>)
  set_target_properties(
    ${main_lib}_shared
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE Fortran
               OUTPUT_NAME rdb_core
               VERSION ${PROJECT_VERSION}
               SOVERSION ${PROJECT_VERSION_MAJOR})

  # NVHPC's `pgacclnk` two-stage device link picks up a PRE-EXISTING
  # librdb_core.so from the build tree and hands it to `ld` as an input, which
  # then refuses with "input file 'librdb_core.so' is the same as output file"
  # and exits 2 -- while leaving the stale .so in place, so a test run
  # immediately afterwards silently executes old code.  Deleting the output
  # before linking is the reliable fix.  (Not a CMake dependency cycle: the .so
  # is nowhere on the link line, and it reproduces under Ninja too.)
  add_custom_command(
    TARGET ${main_lib}_shared
    PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E rm -f $<TARGET_FILE:${main_lib}_shared>
    COMMENT "Removing stale librdb_core.so before relink (NVHPC pgacclnk)")
endif()

# Install rules. Only the DISTRIBUTABLE shared core (RDB_BUILD_SHARED=ON) is
# installed -- a .so that exists solely to speed up RDB_TESTS_LINK_SHARED test
# builds is a dev-only artifact. The C ABI header ships with it so the install
# is link-ready for a C/C++ or ctypes/cffi caller without this source tree.
if(RDB_BUILD_SHARED)
  install(
    TARGETS ${main_lib}_shared
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${PROJECT_SOURCE_DIR}/include/rdb_ocean.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# Distributable core linked by the app + FFI: a SHARED .so when requested, else
# the static core.  (Tests pick their core per-target in tests/CMakeLists.txt.)
if(RDB_BUILD_SHARED)
  set(core_dist ${main_lib}_shared)
  message(
    STATUS
      "librdb_core: SHARED (.so) for app/FFI; STATIC core for tests/benchmarks")
else()
  set(core_dist ${main_lib})
endif()

# Executable — depends on rdb_driver, which is part of the NetCDF I/O subsystem.
# Without NetCDF the namelist-driven simulation has no way to write output, so
# we skip the executable entirely; the library + benchmarks are still useful for
# portability testing.
if(RDB_ENABLE_NETCDF)
  add_executable(${exe_name} ${PROJECT_SOURCE_DIR}/app/main.F90)
  target_include_directories(${exe_name} PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                                 "${RDB_SHARED_INCLUDE_DIR}")
  target_link_libraries(${exe_name} PRIVATE ${core_dist})
endif()

# Namelist knob-reference generator (host-only). Walks the strict schema and
# renders every knob: `./rdb_nml_doc docs/generated_nml_knobs.md`.
add_executable(rdb_nml_doc ${PROJECT_SOURCE_DIR}/app/nml_doc.F90)
target_include_directories(rdb_nml_doc PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                               "${RDB_SHARED_INCLUDE_DIR}")
target_link_libraries(rdb_nml_doc PRIVATE ${core_dist})

# Same schema walk, JSON sink. Feeds tools/gen_python_config.py, which generates
# the checked-in python/rdb/_config_generated.py: `./rdb_nml_json <out.json>`.
add_executable(rdb_nml_json ${PROJECT_SOURCE_DIR}/app/nml_json.F90)
target_include_directories(rdb_nml_json PRIVATE "${CMAKE_BINARY_DIR}/modules"
                                                "${RDB_SHARED_INCLUDE_DIR}")
target_link_libraries(rdb_nml_json PRIVATE ${core_dist})

if(RDB_ENABLE_TESTING)
  add_subdirectory(tests)
endif()
if(RDB_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# Coverage: configure -DRDB_ENABLE_COVERAGE=ON, then `cmake --build build
# --target coverage` -> build/coverage_report/index.html.
if(RDB_ENABLE_COVERAGE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include(coverage)
  rdb_add_coverage_target()
endif()
