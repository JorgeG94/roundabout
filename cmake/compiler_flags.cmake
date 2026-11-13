# Compiler-specific flags for Roundabout.
#
# Mutates CMAKE_Fortran_FLAGS / CMAKE_Fortran_FLAGS_{DEBUG,RELEASE} and may
# append to RDB_COMPILE_DEFS.  Kept as a global-flag include (rather than
# target-scoped INTERFACE libraries) because GPU offload (-stdpar=gpu / -acc)
# must be applied uniformly at compile *and* link time across every Fortran TU
# and dependent target — easiest to guarantee via the global flag set.

# Compilers that expose the `openacc` Fortran module + acc_set_device_num API
# regardless of whether OpenACC is the active codegen backend. Used by
# rdb_comm_env to bind the stdpar/OpenACC device under multi-GPU MPI even on the
# OpenMP backend (do concurrent on these compilers still dispatches through the
# OpenACC runtime). Extend the list as new toolchains gain support (AMD flang
# OpenACC, etc).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID
                                                 STREQUAL "Cray")
  list(APPEND RDB_COMPILE_DEFS RDB_HAS_OPENACC_RUNTIME)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  # Common NVHPC flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree -Mbackslash")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Mbounds -Mchkptr -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  if(RDB_VERBOSE_COMPILE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Minfo=all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=all")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "${CMAKE_Fortran_FLAGS_RELEASE} -Minfo=accel,opt")
  endif()

  if(RDB_ENABLE_GPU)
    # NVHPC consumes RDB_GPU_ARCH as `-gpu=ccNN`; the default (gfx90a, for
    # AMD/LLVM-flang) is invalid here, so fail fast with guidance rather than
    # letting NVHPC emit a cryptic codegen error.
    if(NOT RDB_GPU_ARCH MATCHES "^cc[0-9]+$")
      message(
        FATAL_ERROR
          "NVHPC GPU build needs RDB_GPU_ARCH=ccNN (e.g. cc70/cc80/cc90); got "
          "'${RDB_GPU_ARCH}'. The default is gfx90a for AMD/LLVM-flang — pass "
          "-DRDB_GPU_ARCH=cc70 for an NVIDIA NVHPC build.")
    endif()
    if(RDB_PARALLEL_BACKEND STREQUAL "openacc")
      # -acc=gpu is needed so OpenACC directives like `!$acc routine seq` on the
      # inlined `*_impl` pure helpers (e.g. in `rdb_kernel_flux_unstr` and
      # `rdb_kernel_extrapolate_unstr`) are explicitly enabled. Without it NVHPC
      # may parse the directives but not generate device-callable versions of
      # the helpers, leading to runtime CUDA_ERROR_ILLEGAL_ADDRESS inside `do
      # concurrent` GPU kernels that span module boundaries.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=gpu -acc=gpu -gpu=${RDB_GPU_ARCH},mem:separate"
      )
      message(
        STATUS
          "GPU enabled (do concurrent + OpenACC, NVHPC): targeting ${RDB_GPU_ARCH}"
      )
    else() # openmp
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=gpu -mp=gpu -gpu=${RDB_GPU_ARCH},mem:separate"
      )
      message(
        STATUS
          "GPU enabled (do concurrent + OpenMP target, NVHPC): targeting ${RDB_GPU_ARCH}"
      )
    endif()
    # NVHPC 26.5 reads loop trip-counts from the device copy, inserting an extra
    # per-kernel data refresh in multi-loop `!$acc kernels` regions (the
    # barotropic fast loop) — ~2x on the barotropic solver. Revert to host
    # trip-counts until the compiler regression (NVIDIA TPR #38714) is fixed.
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 26.5)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gpu=tripcount:host")
    endif()
    list(APPEND RDB_COMPILE_DEFS RDB_GPU_OFFLOAD)
  elseif(RDB_ENABLE_THREADS)
    if(RDB_PARALLEL_BACKEND STREQUAL "openacc")
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=multicore -acc=multicore")
      message(STATUS "CPU multicore (do concurrent + OpenACC threaded, NVHPC)")
    else() # openmp
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=multicore -mp=multicore")
      message(STATUS "CPU multicore (do concurrent + OpenMP host, NVHPC)")
    endif()
  else()
    message(STATUS "CPU serial (no threading)")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} -ffree-form -std=gnu -O3 -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -Wall -Wextra")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops")
  if(RDB_ENABLE_COVERAGE)
    # gcov instrumentation: disable optimisation so line counts map cleanly,
    # apply --coverage to both compile and link (gcc driver expands it to
    # -fprofile-arcs -ftest-coverage and links libgcov).
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    message(STATUS "Coverage enabled: -O0 -g --coverage (GNU/gfortran)")
  endif()
  if(RDB_VERBOSE_COMPILE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopt-info-loop")
    message(
      STATUS
        "GNU verbose compile: loop optimisation diagnostics enabled (-fopt-info-loop)"
    )
  endif()
  if(RDB_ENABLE_GPU)
    # Fail loud: RDB_ENABLE_GPU=ON here is always an explicit user request (the
    # auto-detect default is OFF for GNU) — a silent CPU-only downgrade would
    # misreport what was built.
    message(
      FATAL_ERROR
        "RDB_ENABLE_GPU=ON but gfortran has no `do concurrent` GPU offload "
        "path. Use the NVHPC toolchain (nvfortran) for GPU builds, or "
        "configure with -DRDB_ENABLE_GPU=OFF for a gfortran CPU build.")
  endif()
  if(RDB_ENABLE_THREADS)
    if(RDB_PARALLEL_BACKEND STREQUAL "openmp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
      message(STATUS "CPU multicore (do concurrent + OpenMP host, gfortran)")
    else() # openacc / default
      set(RDB_NPROC
          "0"
          CACHE STRING
                "Thread count for -ftree-parallelize-loops (0 = auto-detect)")
      if(RDB_NPROC GREATER 0)
        set(_gnu_nproc ${RDB_NPROC})
      else()
        include(ProcessorCount)
        ProcessorCount(_gnu_nproc)
        if(_gnu_nproc EQUAL 0)
          set(_gnu_nproc 4)
        endif()
      endif()
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -ftree-parallelize-loops=${_gnu_nproc}")
      message(
        STATUS
          "CPU multicore (do concurrent, gfortran): -ftree-parallelize-loops=${_gnu_nproc}"
      )
    endif()
  else()
    message(STATUS "CPU serial (no threading)")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Classic ifort (CMAKE_Fortran_COMPILER_ID == "Intel"). Kept distinct from ifx
  # ("IntelLLVM") with STREQUAL — MATCHES would conflate the two. No OpenMP here
  # on purpose: ifort's `do concurrent` codegen is mature and the analytical
  # suite passes serially at -O3.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -stand f18")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model=precise")
  if(RDB_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RDB_ENABLE_GPU=ON but classic ifort has no GPU offload path. Use "
        "NVHPC (nvfortran) for GPU builds, or -DRDB_ENABLE_GPU=OFF.")
  endif()
  message(STATUS "ifort host CPU (serial do concurrent)")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # -qopenmp is mandatory for ifx, not optional. At -O2+ ifx miscompiles
  # multi-index `do concurrent` bodies that write via `intent(out)` scalars
  # through pure-subroutine calls — the writes silently become zero. The warning
  # is "Locality information is ignored without one of these command line
  # qualifiers '-qopenmp or -parallel'", and the analytical tests fail wholesale
  # (compute_flux_hll returns zeros for non-trivial gradients). See
  # tmp/ifx_repro/flux_repro.f90 for the standalone repro. Since `do concurrent`
  # is the project's primary parallel construct, we add -qopenmp
  # unconditionally; RDB_PARALLEL_BACKEND only gates GPU target offload here.
  #
  # -fp-model=precise pins away from the default fast-math at -O2+ (FMA
  # contraction + reassociation + reciprocal approx) so symmetry and
  # mass-conservation tests match gfortran / nvhpc IEEE behavior.
  #
  # -xHost is intentionally omitted: Gadi login and compute nodes can have
  # different ISA, so a login-node build with -xHost may emit instructions
  # unsupported on the compute node.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -stand f18 -qopenmp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model=precise")
  if(RDB_PARALLEL_BACKEND STREQUAL "openmp" AND RDB_ENABLE_GPU)
    # -fopenmp-target-do-concurrent lowers `do concurrent` onto the OpenMP
    # target region — ifx's do-concurrent GPU-offload path. Without it the DC
    # loops run on the host even with the target backend selected.
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -fopenmp-targets=spir64 -fopenmp-target-do-concurrent"
    )
    message(STATUS "ifx OpenMP target offload enabled "
                   "(-fopenmp-targets=spir64 -fopenmp-target-do-concurrent)")
  elseif(RDB_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RDB_ENABLE_GPU=ON with ifx requires the OpenMP backend "
        "(-DRDB_PARALLEL_BACKEND=openmp -> -fopenmp-targets=spir64); the "
        "'${RDB_PARALLEL_BACKEND}' backend has no ifx GPU path. Either add "
        "-DRDB_PARALLEL_BACKEND=openmp or drop -DRDB_ENABLE_GPU=ON.")
  else()
    message(STATUS "ifx host CPU (-qopenmp, do concurrent)")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ")
  if(RDB_PARALLEL_BACKEND STREQUAL "openmp")
    if(RDB_ENABLE_GPU)
      # OpenMP target offload to an AMD/NVIDIA GPU.  `--offload-arch` selects
      # the device (RDB_GPU_ARCH, default gfx90a = MI200/MI250); applied at
      # compile AND link (global flags) so the offload images are produced +
      # bundled.
      #
      # -fPIC is REQUIRED on the AMD offload path: LLVM Flang emits
      # position-dependent relocations (R_X86_64_32/32S) by default, which
      # ld.lld rejects ("cannot be used against local symbol; recompile with
      # -fPIC") when the offload objects are bundled/device-linked. Scoped to
      # the GPU-offload build only (this branch is already gated on
      # RDB_ENABLE_GPU) so host-only builds are unaffected.
      #
      # -fdo-concurrent-to-openmp=device lowers `do concurrent` onto an OpenMP
      # target region (the GPU); =host would map to host multicore, =none leaves
      # it serial. This is the Flang DC-offload switch, mirroring ifx's
      # -fopenmp-target-do-concurrent above.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -fopenmp --offload-arch=${RDB_GPU_ARCH} -fPIC -fdo-concurrent-to-openmp=device"
      )
      message(
        STATUS
          "LLVM Flang OpenMP target offload "
          "(--offload-arch=${RDB_GPU_ARCH} -fdo-concurrent-to-openmp=device)")
      list(APPEND RDB_COMPILE_DEFS RDB_GPU_OFFLOAD)
    else()
      # Host multicore: map `do concurrent` onto host OpenMP threads.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -fopenmp -fdo-concurrent-to-openmp=host")
      message(STATUS "LLVM Flang OpenMP host "
                     "(-fopenmp -fdo-concurrent-to-openmp=host)")
    endif()
  elseif(RDB_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RDB_ENABLE_GPU=ON with LLVM Flang requires the OpenMP backend "
        "(-DRDB_PARALLEL_BACKEND=openmp -> --offload-arch); the "
        "'${RDB_PARALLEL_BACKEND}' backend has no Flang GPU path. Either "
        "add -DRDB_PARALLEL_BACKEND=openmp or drop -DRDB_ENABLE_GPU=ON.")
  else()
    # No threading: leave `do concurrent` serial (no OpenMP mapping).
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -fdo-concurrent-to-openmp=none")
    message(STATUS "LLVM Flang host CPU "
                   "(serial do concurrent, -fdo-concurrent-to-openmp=none)")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  # Only -O3. Do NOT add `--cpp` or `--std=`: CMake's own LFortran module
  # already compiles with `--cpp-infer`, and LFortran rejects a second `--std`
  # outright ("--std: At Most 1 required but received 2"), which kills the
  # dependency builds before the project is reached. pic likewise sets no
  # standard flags for LFortran.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "")
  set(CMAKE_Fortran_FLAGS_RELEASE "")
endif()

# FMA contraction policy (CPU builds): OFF, for decomposition reproducibility.
#
# The contract is "a decomposed run equals the 1-rank run BITWISE on owned
# cells".  On CPU builds the vectoriser FMA-contracts a loop's vector body and
# its scalar remainder/peel differently, so a cell's rounding depends on its
# position in the TILE, hence on the tile width.  Almost every kernel loops over
# full storage (writing open-edge ghosts that are never exchanged) and owned
# cells have the same body-vs-remainder exposure, and WHICH loops are exposed
# depends on the compiler version -- so a per-file list is not an honest scope.
# One blanket flag per compiler replaces the per-file `-Mnofma` /
# `-ffp-contract=off` lines this file used to carry.  Applied at directory scope
# (project AND fetched dependencies: pic, pic-mpi, test-drive), after the config
# flags, like RDB_EXTRA_FORTRAN_FLAGS below.  The GPU build is NOT touched: its
# device kernels have no vector remainder and the GPU decomposition gate passes
# without it (its HOST code -- setup, host loops -- still contracts). This does
# NOT make vectorised libm calls (glibc libmvec `_ZGV*`) position-independent;
# that is a separate, tracked follow-up (CLAUDE.md, Gotchas).
#
# History: the per-file cases this blanket flag replaces (each caught by the
# gate named).  Entries are independent; GNU and NVHPC cases are separate.
#
# cmake-format: off
#
# NVHPC, rdb_ocean_horizontal_viscosity.F90 (2026-10): the vector body and the
#   scalar remainder contract differently; pred_corr carries `hvisc_du_visc`
#   into the prognostics, and z* + visc_rem on 4x1 drifted 1e-8 in eta in 24
#   steps (decomp_bitid `visc_rem_zstar`, `hvisc_thin_layer`).
#
# NVHPC, rdb_barotropic_coupling.F90 (2026-10): `compute_gtot_faces`'s
#   `p_sum = p_sum + h_face*pbce` in a `do k` inside the outer face loop,
#   ~7e-15 on `bt_visc_rem_u/v`, `epbl_kd_int` (compat matrix cell c023,
#   DECOMP leg).
#
# NVHPC, rdb_remap_column.F90 (2026-10): under `-fast -tp=haswell` one face's
#   remapped `v_face_y_layer` in the first owned column of the 4x1 split's last
#   tile was 1 ulp off the 2x2 / 1x4 value; the ulp reached the mass flux and
#   the EFP `mass_out` budget kept it (decomp_bitid `gm_ebt_struct/ssp_rk2`,
#   np 4).  A per-file bisection of `-Mnofma` isolated this file alone.
#
# GNU, rdb_barotropic_coupling.F90 (2026-10): same `compute_gtot_faces` loop
#   under `-O3 -march=x86-64-v3`, same symptom and gate as the NVHPC case
#   above (compat cell c023, 1 rank vs 4x1, ~7e-15 relative).
#
# GNU, rdb_ocean_horizontal_viscosity.F90 (2026-10): the stress-form path with
#   `bound_thin_layer` left `v_face_y_layer` / `u_face_x_layer` 1 ulp apart
#   between splits (decomp_bitid `hvisc_thin_layer/ssp_rk2`, 2x1 and 4x1).
#
# GNU, rdb_remap_column.F90 (2026-10): compat cell c023 (sigma, ssp_rk2,
#   EPBL + tidal + Bryan-Lewis, Leith, MLE, OBC) was bitwise only by luck once
#   the hvisc line moved Leith's rounding; the 4x1 split's `bt_visc_rem_u/v`
#   and `epbl_kd_int` came out a few ulp off the 1-rank run.
#
# NVHPC, rdb_continuity (found on this branch, not a per-file line on beta):
#   with a Flather east edge the open-edge ghost columns evolved by the
#   continuity chain reach the owned boundary face through the horizontal
#   viscosity; first owned difference at step 12 (6.7e-16), 3.6e-13 by step 24
#   (decomp_bitid `obc_epbl`, np 4).
# cmake-format: on
if(NOT RDB_ENABLE_GPU)
  set(RDB_FMA_POLICY "OFF (CPU build, decomposition reproducibility)")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Mnofma>")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Contraction is already off at -O0 (Debug); harmless there.
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ffp-contract=off>")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID
                                                       STREQUAL "IntelLLVM")
    # UNTESTED (no ifx/ifort here).  `-fp-model=precise` alone does not stop FMA
    # formation; `-no-fma` is the documented Intel switch for it.
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-no-fma>")
    set(RDB_FMA_POLICY
        "OFF, -no-fma (UNTESTED; CPU build, decomposition reproducibility)")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
    # UNTESTED (no flang here).  flang-new accepts clang's -ffp-contract.
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ffp-contract=off>")
    set(RDB_FMA_POLICY
        "OFF, -ffp-contract=off (UNTESTED; CPU build, decomposition reproducibility)"
    )
  else()
    set(RDB_FMA_POLICY
        "compiler default (no contraction switch known for ${CMAKE_Fortran_COMPILER_ID})"
    )
  endif()
else()
  set(RDB_FMA_POLICY "ON (GPU build)")
endif()

# Applied LAST, and deliberately via add_compile_options rather than by
# appending to CMAKE_Fortran_FLAGS. CMake emits CMAKE_Fortran_FLAGS, then
# CMAKE_Fortran_FLAGS_<CONFIG>, then COMPILE_OPTIONS -- so only this form can
# override something set on the Release line, and overriding is the whole point
# of an escape hatch.
#
# Directory scope, set before any target or fetched dependency exists, so the
# subprojects compile with it too. That matters for the `-tp` case: a core built
# for one CPU target and a dependency built for another is exactly the mismatch
# this is meant to prevent.
if(RDB_EXTRA_FORTRAN_FLAGS)
  separate_arguments(_rdb_extra_flags NATIVE_COMMAND
                     "${RDB_EXTRA_FORTRAN_FLAGS}")
  foreach(_rdb_flag IN LISTS _rdb_extra_flags)
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${_rdb_flag}>")
  endforeach()
  unset(_rdb_flag)
  unset(_rdb_extra_flags)
endif()
