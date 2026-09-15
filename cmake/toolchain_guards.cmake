# Compiler-identity guards. Include AFTER options.cmake (it overrides some of
# them) and BEFORE compiler_flags.cmake. Plain include(): the overrides below
# must land in the caller's scope.

# GCC implements the F2018 `do concurrent` locality specifiers that the ocean
# kernels are written with only from GCC 15; the FATAL_ERROR spells out why.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION
                                                VERSION_LESS 15)
  message(
    FATAL_ERROR
      "gfortran ${CMAKE_Fortran_COMPILER_VERSION} is too old to build Roundabout.\n"
      "\n"
      "  Required: GNU Fortran >= 15.0\n"
      "  Found:    ${CMAKE_Fortran_COMPILER_VERSION}\n"
      "            (${CMAKE_Fortran_COMPILER})\n"
      "\n"
      "WHY: the ocean dynamical core is written with F2018 `do concurrent` "
      "locality specifiers -- `local(...)`, `local_init(...)`, `shared(...)` "
      "-- which are how the kernels declare per-iteration private state so "
      "they parallelise on CPU and GPU from the SAME source. GCC implements "
      "those clauses only from GCC 15. Older gfortran accepts `do concurrent` "
      "itself but rejects the locality clauses, so the build would fail with "
      "a wall of syntax errors inside the kernels instead of this message.\n"
      "\n"
      "HOW TO FIX, pick one:\n"
      "  * install GCC >= 15 and point CMake at it:\n"
      "      cmake -B build -S . -DCMAKE_Fortran_COMPILER=/path/to/gfortran-15\n"
      "  * use the NVIDIA HPC SDK (nvfortran), which supports the locality\n"
      "    specifiers and is the toolchain for GPU builds:\n"
      "      cmake -B build -S . -DCMAKE_Fortran_COMPILER=nvfortran\n"
      "  * use Intel ifx (also supported).\n"
      "\n"
      "See FORTRAN_STYLE.md for why the locality specifiers are load-bearing "
      "rather than stylistic.")
endif()

# CMake's Apple platform rules inject `-install_name @rpath/lib….dylib`, which
# lfortran's CLI rejects. The static-core path is fully supported, so force it.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  set(RDB_BUILD_SHARED OFF)
  set(RDB_TESTS_LINK_SHARED OFF)
  message(
    STATUS
      "LFortran: forcing the static core (no shared .dylib — lfortran lacks the macOS -install_name flag)"
  )
endif()

# Flang/AMD offload advisory: the per-work-item device stack frame of the heavy
# per-column solvers (e.g. kappa-shear) scales with NZ_STACK_MAX and can exceed
# the AMDGPU static stack-frame limit on some targets (MI250X / gfx90a).
if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang") # matches "Flang" and "LLVMFlang"
  message(
    STATUS
      "Flang detected: NZ_STACK_MAX=${RDB_NZ_STACK_MAX}. If a GPU-offload device "
      "link reports 'stack frame size exceeds limit', reconfigure with a smaller "
      "-DRDB_NZ_STACK_MAX=<n> (must stay >= nz + 1 for your vertical layer count)."
  )
endif()
