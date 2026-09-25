# External dependencies: test-drive, pic, NetCDF-Fortran, MPI, pic-mpi. Plain
# include() so the RDB_COMPILE_DEFS appends and RDB_MPI_LANG reach the caller.

# Do NOT build our dependencies' own test suites. They cost build time, they
# register in our ctest (86 pic + 6 test-drive tests -- the reason every
# instruction here says `ctest -R rdb`), and a failure in one of them fails OUR
# build for a reason that is not ours: test-drive v0.6.1's `test-drive-tester`
# does not link under LLVM Flang, which took the whole build down.
#
# This project's own tests are unaffected -- they go through `enable_testing()`
# and `add_test` directly and never consult BUILD_TESTING.
foreach(
  _rdb_dep_tests
  BUILD_TESTING # test-drive's fallback
  TESTDRIVE_BUILD_TESTING # test-drive
  PIC_ENABLE_TESTING # pic
  ENABLE_TESTING) # pic-mpi
  set(${_rdb_dep_tests}
      OFF
      CACHE BOOL "Build a dependency's own test suite" FORCE)
endforeach()
unset(_rdb_dep_tests)

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
  # C has to be enabled IN THIS SCOPE, before the search. A netcdf-fortran CMake
  # config built against a netcdf-c that itself ships a config drags the whole C
  # chain in behind it:
  #
  # netCDF-FortranConfig -> find_dependency(netCDF) -> netCDFConfig     ->
  # find_dependency(HDF5) -> hdf5-config    -> find_dependency(Threads) ->
  # FindThreads: "only works if either C or CXX language is enabled"
  #
  # Note this is a SCOPE problem, not a missing-compiler one. `pic` and
  # `pic-mpi` both declare `LANGUAGES Fortran C`, so a C compiler is already
  # found and configured by the time we get here -- but they are added as
  # subprojects, and `CMAKE_C_COMPILER_LOADED` propagates DOWN into a
  # subdirectory, never back UP into ours. FindThreads tests that variable in
  # the calling scope, so it fails even though the compiler is right there.
  #
  # conda-forge installs are what expose it, because they ship a config for BOTH
  # netcdf-c and HDF5. Where either lacks one, netcdf-fortran bakes an absolute
  # path into its own config and never re-enters the chain -- which is why an
  # apt-style netcdf-c configures fine without this line.
  enable_language(C)

  find_package(netCDF-Fortran QUIET HINTS "$ENV{NETCDFF_DIR}"
               "$ENV{NETCDF_FORTRAN_DIR}")

  if(netCDF-Fortran_FOUND)
    message(
      STATUS "NetCDF-Fortran found via CMake config: ${netCDF-Fortran_DIR}")
  else()
    find_package(PkgConfig)
    if(PkgConfig_FOUND)
      pkg_check_modules(NETCDFF IMPORTED_TARGET netcdf-fortran)
    endif()
    if(NETCDFF_FOUND)
      add_library(netCDF::netcdff ALIAS PkgConfig::NETCDFF)
      message(
        STATUS "NetCDF-Fortran found via pkg-config: ${NETCDFF_LINK_LIBRARIES}")
    else()
      # Worth a real message: this is the single most common way a first build
      # fails, and `netcdf-fortran` is the ONE dependency that has to match the
      # Fortran compiler (`.mod` files are compiler-specific).
      message(
        FATAL_ERROR
          "netcdf-fortran not found, and it must be a build made by THIS "
          "Fortran compiler (${CMAKE_Fortran_COMPILER_ID}) -- Fortran .mod "
          "files are not portable between compilers.\n"
          "Ways to get one:\n"
          "  * Build just the wrapper (~30 s) on top of any netcdf-c "
          "you already have.\n"
          "    The netcdf-c does NOT have to match your compiler:\n"
          "      tools/build_netcdf_fortran.sh --fc <your compiler>\n"
          "    then configure with the prefix it prints:\n"
          "      -DCMAKE_PREFIX_PATH=<prefix>\n"
          "  * No netcdf-c at all, only HDF5? Build both into one prefix:\n"
          "      tools/build_netcdf.sh --fc <your compiler>\n"
          "    (docs/howto/deploy_without_netcdf.md)\n"
          "  * Your site's module tree, if it ships one per compiler.\n"
          "  * conda-forge's netcdf-fortran -- GFORTRAN BUILDS ONLY.\n"
          "  * environments/spack.yaml for the underlying C libraries.\n"
          "Or drop the I/O subsystem entirely and build kernels + tests "
          "only:\n"
          "      cmake -B build -S . -DRDB_ENABLE_NETCDF=OFF\n"
          "CMake looked at the netCDF-Fortran CMake config package (hinted "
          "by $NETCDFF_DIR / $NETCDF_FORTRAN_DIR) and at pkg-config's "
          "netcdf-fortran.")
    endif()
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
