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
    "openmp"
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
# CUDA-aware MPI. DECLARES that the MPI library is CUDA-aware (a GPU build
# always hands device buffers to MPI; without the declaration a multi-rank GPU
# run is refused at startup). Default: ON exactly when MPI and GPU offload are
# both requested (RDB_ENABLE_GPU=ON already fails loud on a compiler with no GPU
# path), OFF otherwise. Provenance: RDB_CUDA_AWARE_MPI_AUTO (INTERNAL) holds the
# value this file auto-chose. While the cache still equals it the value is
# re-derived on every configure (so toggling MPI/GPU in the same build dir
# follows); any other value is explicit and never touched. An explicit -D equal
# to the previously auto-chosen value is indistinguishable and keeps following
# the default; a re-derivation that changes the value is reported. ON without
# MPI + GPU is a configure FATAL_ERROR, like RDB_ENABLE_GPU.
if(RDB_ENABLE_MPI AND RDB_ENABLE_GPU)
  set(_rdb_cuda_aware_default ON)
else()
  set(_rdb_cuda_aware_default OFF)
endif()
if(NOT DEFINED RDB_CUDA_AWARE_MPI)
  set(RDB_CUDA_AWARE_MPI_AUTO
      "${_rdb_cuda_aware_default}"
      CACHE INTERNAL "")
elseif(DEFINED RDB_CUDA_AWARE_MPI_AUTO AND "${RDB_CUDA_AWARE_MPI}" STREQUAL
                                           "${RDB_CUDA_AWARE_MPI_AUTO}")
  if(NOT "${RDB_CUDA_AWARE_MPI}" STREQUAL "${_rdb_cuda_aware_default}")
    message(
      STATUS
        "RDB_CUDA_AWARE_MPI: auto default re-derived ${RDB_CUDA_AWARE_MPI} -> "
        "${_rdb_cuda_aware_default} (MPI/GPU changed); to pin a value use a "
        "fresh build dir with -DRDB_CUDA_AWARE_MPI=<value>, or "
        "-U RDB_CUDA_AWARE_MPI first, then -D...")
  endif()
  set(RDB_CUDA_AWARE_MPI
      "${_rdb_cuda_aware_default}"
      CACHE BOOL "" FORCE)
  set(RDB_CUDA_AWARE_MPI_AUTO
      "${_rdb_cuda_aware_default}"
      CACHE INTERNAL "")
else()
  unset(RDB_CUDA_AWARE_MPI_AUTO CACHE)
endif()
option(
  RDB_CUDA_AWARE_MPI
  "Declare the MPI library CUDA-aware (device buffers go to MPI; required for multi-rank GPU runs). Default ON with MPI + GPU"
  ${_rdb_cuda_aware_default})
if(RDB_CUDA_AWARE_MPI AND NOT (RDB_ENABLE_MPI AND RDB_ENABLE_GPU))
  message(
    FATAL_ERROR
      "RDB_CUDA_AWARE_MPI=ON requires both RDB_ENABLE_MPI=ON and "
      "RDB_ENABLE_GPU=ON (got MPI=${RDB_ENABLE_MPI}, GPU=${RDB_ENABLE_GPU}). "
      "Pass -DRDB_CUDA_AWARE_MPI=OFF (or -U RDB_CUDA_AWARE_MPI), or enable both."
  )
endif()
option(
  RDB_ENABLE_LEGACY_MPI
  "Use legacy 'use mpi' instead of 'use mpi_f08' (must match pic-mpi build)"
  OFF)
option(
  RDB_USE_VAPAA
  "pic-mpi was built with PIC_USE_VAPAA: link MPI::MPI_C and let vapaa supply mpi_f08"
  OFF)

# Escape hatch for flags this file has no opinion about. The motivating case is
# `-tp` on NVHPC: nvfortran targets the BUILD HOST's CPU by default -- nothing
# in compiler_flags.cmake sets it, and neither `-fast` nor the build type
# changes that -- so a binary built on one machine can die with SIGILL on
# another. CI pins it to a baseline on the NVHPC legs; a cluster user building
# on a login node for differently-specced compute nodes wants the same.
set(RDB_EXTRA_FORTRAN_FLAGS
    ""
    CACHE STRING
          "Extra Fortran flags, appended after every flag this project sets")
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
