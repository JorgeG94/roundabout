# External dependencies: test-drive, pic, NetCDF-Fortran, MPI, pic-mpi. Plain
# include() so the RDB_COMPILE_DEFS appends and RDB_MPI_LANG reach the caller.

if(RDB_ENABLE_TESTING)
  enable_testing()
  # Expose unit-test-only module symbols (e.g. the unsplit reference routines +
  # directionally-split building blocks in rdb_continuity) to the test suite via
  # `#ifdef RDB_ENABLE_TESTING`, without widening the production API.
  list(APPEND RDB_COMPILE_DEFS RDB_ENABLE_TESTING)
  if(NOT TARGET test-drive::test-drive)
    find_package("test-drive" REQUIRED)
  endif()
endif()

if(NOT TARGET pic::pic)
  find_package("pic" REQUIRED)
endif()

# NetCDF-Fortran: CMake config first (modern installs), then pkg-config
# (autotools installs).
if(RDB_ENABLE_NETCDF)
  find_package(netCDF-Fortran QUIET HINTS "$ENV{NETCDFF_DIR}"
               "$ENV{NETCDF_FORTRAN_DIR}")

  if(netCDF-Fortran_FOUND)
    message(
      STATUS "NetCDF-Fortran found via CMake config: ${netCDF-Fortran_DIR}")
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NETCDFF REQUIRED IMPORTED_TARGET netcdf-fortran)
    add_library(netCDF::netcdff ALIAS PkgConfig::NETCDFF)
    message(
      STATUS "NetCDF-Fortran found via pkg-config: ${NETCDFF_LINK_LIBRARIES}")
  endif()
else()
  list(APPEND RDB_COMPILE_DEFS RDB_NO_NETCDF)
  message(
    STATUS "NetCDF disabled: skipping I/O subsystem (driver, output, gauges, "
           "restart, forcing, nesting, main executable). "
           "Kernels and benchmarks still build.")
endif()

# pic-mpi is ALWAYS a dependency: rdb has ONE comm implementation and it always
# calls `pic_mpi_lib`.  RDB_ENABLE_MPI selects which pic-mpi BACKEND gets built
# (real MPI vs pic-mpi's serial one), whether rdb links an MPI library itself,
# and whether the multi-rank test suite is built.
set(PIC_ENABLE_MPI
    ${RDB_ENABLE_MPI}
    CACHE BOOL "Build pic-mpi against MPI (OFF selects its serial backend)"
          FORCE)

if(RDB_ENABLE_MPI)
  if(RDB_ENABLE_LEGACY_MPI AND RDB_USE_VAPAA)
    message(
      FATAL_ERROR
        "RDB_ENABLE_LEGACY_MPI and RDB_USE_VAPAA are mutually exclusive (vapaa supplies mpi_f08 only)."
    )
  endif()
  if(RDB_USE_VAPAA)
    enable_language(C)
    set(RDB_MPI_LANG C)
    # Propagate to pic-mpi: when fetched as a subproject its option() sees this
    # pre-set cache value and skips the default OFF.  Findpic-mpi.cmake also
    # FORCEs it for safety in case the consumer set it on the command line.
    set(PIC_USE_VAPAA
        ON
        CACHE BOOL "Use vapaa to link to a C MPI library" FORCE)
  else()
    set(RDB_MPI_LANG Fortran)
  endif()
  find_package(MPI REQUIRED COMPONENTS ${RDB_MPI_LANG})
  # Use mpirun from PATH rather than CMake's MPIEXEC_EXECUTABLE, which can pick
  # up the wrong MPI installation.
  set(MPIEXEC_EXECUTABLE
      "mpirun"
      CACHE FILEPATH "MPI launcher" FORCE)
  set(MPIEXEC_NUMPROC_FLAG
      "-np"
      CACHE STRING "MPI numproc flag" FORCE)
  list(APPEND RDB_COMPILE_DEFS RDB_ENABLE_MPI)
  if(RDB_CUDA_AWARE_MPI)
    list(APPEND RDB_COMPILE_DEFS RDB_CUDA_AWARE_MPI)
    message(STATUS "CUDA-aware MPI enabled: GPU-direct halo exchange")
  endif()
  if(RDB_ENABLE_LEGACY_MPI)
    list(APPEND RDB_COMPILE_DEFS USE_LEGACY_MPI)
    message(STATUS "MPI interface: legacy (use mpi)")
  elseif(RDB_USE_VAPAA)
    message(STATUS "MPI interface: modern (use mpi_f08 via vapaa, MPI::MPI_C)")
  else()
    message(STATUS "MPI interface: modern (use mpi_f08)")
  endif()
endif()

if(NOT TARGET pic-mpi::pic-mpi)
  find_package("pic-mpi" REQUIRED)
endif()
