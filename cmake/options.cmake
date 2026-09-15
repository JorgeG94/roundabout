# Build options and cache variables. Plain include() (not a function) so every
# entry lands in the caller's scope.

# GPU offload is OPT-IN. An explicit ON with a compiler that has no GPU offload
# path is a configure-time FATAL_ERROR (see cmake/compiler_flags.cmake) — never
# a silent CPU-only downgrade. Reported on the summary's "GPU offload:" line.
option(RDB_ENABLE_GPU "Build with GPU offloading (do concurrent + OpenACC)" OFF)

option(RDB_ENABLE_THREADS
       "Enable multicore threading (do concurrent + OpenACC multicore)" OFF)

# Setting `openmp` on `main` only changes compiler flags — the source still
# carries `!$acc` directives until acc_to_omp.py runs. The OpenMP variant is
# consumed via the auto-regenerated `auto/dc-openmp` branch, which sets this and
# applies the overlay patches under patches/openmp/.
set(RDB_PARALLEL_BACKEND
    "openacc"
    CACHE STRING "Parallel backend: openacc | openmp (advanced; leave default)")
set_property(CACHE RDB_PARALLEL_BACKEND PROPERTY STRINGS openacc openmp)
mark_as_advanced(RDB_PARALLEL_BACKEND)
if(NOT RDB_PARALLEL_BACKEND MATCHES "^(openacc|openmp)$")
  message(FATAL_ERROR "RDB_PARALLEL_BACKEND must be one of: openacc, openmp "
                      "(got '${RDB_PARALLEL_BACKEND}')")
endif()

option(RDB_ENABLE_DOUBLE "Use double precision (real64)" ON)
option(
  RDB_ENABLE_MPI
  "Build multi-rank: link an MPI library and build pic-mpi against it.\nOFF still uses pic-mpi -- its serial backend -- and builds single-rank."
  OFF)
option(RDB_CUDA_AWARE_MPI "Use CUDA-aware MPI for GPU-direct halo exchange" OFF)
option(
  RDB_ENABLE_LEGACY_MPI
  "Use legacy 'use mpi' instead of 'use mpi_f08' (must match pic-mpi build)"
  OFF)
option(
  RDB_USE_VAPAA
  "pic-mpi was built with PIC_USE_VAPAA: link MPI::MPI_C and let vapaa supply mpi_f08"
  OFF)
option(RDB_ENABLE_TESTING "Build test suite" ON)
option(
  RDB_BUILD_BENCHMARKS
  "Build benchmark executables (dev-only; some drive core GPU device routines directly, so they cannot link against a shared core — keep OFF with RDB_BUILD_SHARED)"
  OFF)
option(RDB_VERBOSE_COMPILE "Show detailed compiler diagnostics (-Minfo=all)"
       OFF)
option(RDB_ENABLE_NVTX "Enable NVTX profiling annotations for Nsight Systems"
       OFF)
option(RDB_BUILD_SHARED
       "Build librdb as a shared (.so) library instead of a static archive" OFF)
option(
  RDB_TESTS_LINK_SHARED
  "Link tests against a shared core (.so) for faster test builds. ON by default: only the tests in tests/CMakeLists.txt's RDB_STATIC_CORE_TESTS list (device code in their own source) keep the static core; everything else links the .so. Set OFF to force every test back to the static core."
  ON)
option(
  RDB_ENABLE_NETCDF
  "Build the NetCDF-backed I/O subsystem (driver, output, gauges, restart, forcing, nesting). Disable for portability testing on AMD/Intel where bringing along NetCDF is inconvenient — kernels and benchmarks still build."
  ON)
option(
  RDB_ENABLE_COVERAGE
  "Build with gcov instrumentation (-O0 -g --coverage). GNU/gfortran only; requires lcov + genhtml on PATH for the `coverage` target."
  OFF)

if(RDB_ENABLE_COVERAGE AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(
    WARNING "RDB_ENABLE_COVERAGE=ON requires GNU/gfortran (current: "
            "${CMAKE_Fortran_COMPILER_ID}). Coverage flags NOT applied; the "
            "`coverage` target will not be created.")
endif()

# LFortran 0.64 (alpha) needs a handful of source-level workarounds, all guarded
# by the LFORTRAN_PASSING macro so the default build compiles the REAL code.
# Leave OFF to test a newer lfortran against the unmodified source.
option(
  RDB_LFORTRAN_PASSING
  "Activate the LFortran 0.64 source workarounds (defines LFORTRAN_PASSING)"
  OFF)
if(RDB_LFORTRAN_PASSING)
  add_compile_definitions(LFORTRAN_PASSING)
endif()

# NZ_STACK_MAX in rdb_constants: the per-column fixed-size stack workspace bound
# inside do-concurrent kernels. Must stay >= nz + 1 (`validate_config` refuses
# the run fail-loud otherwise). NOT 2*nz+2 — the Redi neutral-surface locals are
# declared 2*NZ_STACK_MAX+2 and scale with this. 128 covers the largest shipped
# case (nz=90) and test (nz=100).
set(RDB_NZ_STACK_MAX
    "128"
    CACHE
      STRING
      "Max vertical layers for per-column stack workspace (lower for AMD/flang if the device link hits stack-frame limits)"
)

# Default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type (Debug/Release)" FORCE)
endif()

# GPU architecture — compiler-dependent default. AMD/LLVM flang targets gfx90a
# (`--offload-arch`); every other compiler (notably NVHPC, which consumes it as
# `-gpu=ccNN`) defaults to cc70.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(_rdb_default_gpu_arch "gfx90a")
else()
  set(_rdb_default_gpu_arch "cc70")
endif()
set(RDB_GPU_ARCH
    "${_rdb_default_gpu_arch}"
    CACHE STRING "GPU arch: AMD gfxNNN (LLVM flang) or NVIDIA ccNN (NVHPC)")
