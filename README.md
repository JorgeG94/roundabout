# Roundabout

I had to learn about ocean dynamics so I wrote a code to do it. It runs on CPUs and GPUs
with `do concurrent` being the way to access parallelism. MPI is availbale too. 

It is also a way to find bugs on Intel, AMD, and NVIDIA and check for LFortran portability.

The code is called Roundabout in honour of the city of Canberra, which is where I work. It is 
an inland city and has no ocea. I thought it funny to write an ocean code from a non-ocean city,
and one can think that a roundabout is like an Eddie.

## Features

- **C-grid dynamical core** — Arakawa C-grid, continuity-PPM, PV-conserving Sadourny Coriolis with optional WENO PV advection, split-explicit RK2 with a nonlinear barotropic fast loop.
- **Thermodynamics** — Wright (1997) nonlinear EOS (linear + Roquet available), salinity/temperature plus a growable passive-tracer registry (ideal age, pseudo-salt), shortwave penetration, surface buoyancy restoring.
- **Vertical mixing** — PP81 interior + KPP boundary overlay (default on), EPBL, kappa-shear, tidal mixing, convective adjustment, double diffusion, Bryan-Lewis / Henyey backgrounds — all assembled through one `vmix_assemble` floors/ceilings gate.
- **Lateral closures** — Leith, Smagorinsky KH+AH, Leith-biharmonic, anisotropic viscosity, MEKE backscatter, GM/Redi, Fox-Kemper MLE.
- **Vertical coordinates** — sigma, z-sigma hybrid, and the z-star family (z\*-lite, z\*/sigma hybrid, full per-column z\*) through a conservative ALE remap (PCM/PLM/PPM/PQM).
- **Horizontal grids** — Cartesian, spherical lon-lat, MOM6 supergrid mosaic, and tripolar (Murray 1996 bipolar Arctic cap with a north fold).
- **Boundaries** — wall, open (Flather), tidal, clamped, Chapman, sponge, periodic, Orlanski radiation; interior land masking + dynamic wet/dry.
- **Sea ice** — Winton column thermodynamics, multi-category ITD, category transport, C-grid EVP dynamics, frazil + brine coupling.
- **I/O** — per-rank NetCDF (CF-1.8) diagnostics with cadence dispatch and per-diag vcoord remap, restart.
- **MPI multi-GPU** — domain decomposition + halo exchange (incl. CUDA-aware GPU-direct).

See [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md) for the full feature/limitation list, [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) for the enabled-closure ground truth + tunable knobs.

## Backend Strategy

The idea is to use standard Fortran parallelism as much as possible and OpenACC where it can't be helped. We have
a CI and a linter to make sure any openacc directive is portable to OpenMP for AMD and Intel GPUs.

- **`do concurrent`** for all data-parallel loops — offloads to GPU with `-stdpar=gpu`, runs threaded with `-stdpar=multicore`, or runs serially without flags
- **`!$acc parallel loop reduction(...)`** for reduction loops (CFL timestep, CG dot products) — inert comments on non-OpenACC compilers, falls back to sequential
- **`!$acc enter/exit data`** for GPU memory management — written directly with no `#ifdef` wrapping (inert comments on CPU builds)

No vendor-specific extensions — the same source compiles with NVHPC (GPU + multicore), gfortran, and ifx.

| Build mode | Compiler | Flags | `do concurrent` | Reductions |
|------------|----------|-------|-----------------|------------|
| **GPU** | nvfortran | `-stdpar=gpu -acc=gpu` | GPU | `!$acc parallel loop` |
| **Multicore** | nvfortran | `-stdpar=multicore -acc=multicore` | threaded | sequential |
| **Serial** | nvfortran | (none) | sequential | sequential |
| **CPU (gcc)** | gfortran | `-ftree-parallelize-loops=N` | threaded | sequential |

You can use the alternative backends if you checkout the branches:

- auto/dc-openmp
- auto/openmp

### auto/dc-openmp

Automatically generated from the main branch upon Pull-Request. It uses OpenMP of data mvoement and do concurrent for compute.

### auto/openmp

Automatically generated from the main branch upon Pull-Request. It uses OpenMP of both data movement and compute.


## Build

```bash
# GPU build (NVHPC, default)
module load cmake nvhpc
cmake -B build -S . && cmake --build build
cd build && ctest --output-on-failure

# NVHPC multicore CPU (threaded)
cmake -B build -S . -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_THREADS=ON

# NVHPC serial CPU (no threading)
cmake -B build -S . -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_THREADS=OFF

# gfortran (CPU, auto-parallel via -ftree-parallelize-loops)
cmake -B build -S . -DCMAKE_Fortran_COMPILER=gfortran -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_THREADS=ON

# gfortran with explicit thread count (this is very weird to work with with GNU)
cmake -B build -S . -DCMAKE_Fortran_COMPILER=gfortran -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_THREADS=ON -DRDB_NPROC=8

# MPI multi-GPU build
cmake -B build -S . -DRDB_ENABLE_MPI=ON && cmake --build build

# NetCDF-free build (portability/CI testing on AMD/Intel where pulling in NetCDF is inconvenient or annoying).
# Skips the I/O subsystem (driver, diagnostics, restart, forcing, main executable);
# kernels and benchmarks still build, and unit tests that don't need I/O run.
cmake -B build -S . -DRDB_ENABLE_NETCDF=OFF && cmake --build build
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `RDB_ENABLE_GPU` | `OFF` | GPU offloading (`-stdpar=gpu -acc=gpu`); opt-in on the NVHPC toolchain. Explicit `ON` with a compiler that has no GPU path is a configure-time error (no silent CPU downgrade); the value prints on the `GPU offload:` configure line. Sets `RDB_GPU_OFFLOAD` preprocessor define (used only for the CUDA-aware MPI selection and the multi-GPU build-config FATAL guard) |
| `RDB_ENABLE_THREADS` | `OFF` | Multicore threading (`-stdpar=multicore -acc=multicore` for NVHPC; `-ftree-parallelize-loops=N` for GCC). Ignored when `RDB_ENABLE_GPU=ON` |
| `RDB_ENABLE_MPI` | `OFF` | Multi-rank: link an MPI library and build pic-mpi against it. `OFF` still uses pic-mpi — its serial backend — and builds single-rank |
| `RDB_CUDA_AWARE_MPI` | `OFF` | Use GPU-direct MPI for halo exchange |
| `RDB_ENABLE_DOUBLE` | `ON` | Double precision (`real64`) |
| `RDB_GPU_ARCH` | `cc70` | GPU compute capability (e.g. `cc70`, `cc80`, `cc90`) |
| `RDB_ENABLE_NETCDF` | `ON` | Build the NetCDF-backed I/O subsystem (driver, diagnostics, restart, forcing, main executable). Disable for portability testing on platforms where NetCDF is inconvenient — kernels and benchmarks still build |
| `RDB_BUILD_SHARED` | `OFF` | Build `librdb_core` as a shared library instead of a static archive |

### Dependencies

- **Required**: [pic](https://github.com/JorgeG94/pic/) (types, strings, timers, logger) — fetched automatically via CMake
- **Required**: [pic-mpi](https://github.com/JorgeG94/pic-mpi/) (MPI wrappers) — fetched automatically via CMake
- **Required**: [test-drive](https://github.com/JorgeG94/test-drive) (unit testing) — fetched automatically
- **Optional**: NetCDF-Fortran (bathymetry I/O, diagnostics, restart, forcing) — required for the main executable; skip with `-DRDB_ENABLE_NETCDF=OFF` for kernel-only portability builds
- **Optional**: MPI (multi-GPU, via `mpi_f08`)

## Project Layout

```
src/
  core/
    rdb_{grid,tracer,barotropic_state,multilayer_state}.F90   shared state types
    ocean/         the C-grid dyn-core — god-state slot map (see its README)
    ice/           sea ice — thermo/, itd/, transport/, dynamics/, coupling/
  ALE/  tracer/  equation_of_state/  pressure_force/  parameterizations/  framework/   shared physics
  driver/  io/  comm/    lifecycle, NetCDF I/O, MPI facade
app/main.F90       thin entry — calls driver_run
tests/  benchmarks/  validation_examples/  docs/
```

## Validation Tests

| Test | What it validates |
|------|-------------------|
| Double gyre (MOM6 reference) | Full split-RK2 dyn-core |
| Geostrophic adjustment | Coriolis / PGF balance |
| Baroclinic channel, Eady | Baroclinic instability growth rates |
| Lock exchange / seamount | ALE remap conservation under stratification |
| Island at rest | Well-balanced property with land masking |
| Salt/heat conservation | Closed-domain tracer budgets |

Run all tests:
```bash
cd build && ctest --output-on-failure
```

Cross-backend validation:
```bash
./validate.sh                    # all backends
./validate.sh gcc-serial nvhpc-gpu  # selected backends
```

## Documentation

- [`CLAUDE.md`](CLAUDE.md) — project overview, architecture, build, conventions, + the doc map
- [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md) — what the solver can / can't do today
- [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) — which closure/scheme is enabled in which regime + the tunable knobs
- [`src/core/ocean/README.md`](src/core/ocean/README.md) — the dyn-core design contract (conventions, dispatch pipeline, slot map)
- [`docs/REFERENCE.md`](docs/REFERENCE.md) — physics + numerics + namelist reference
- [`FORTRAN_STYLE.md`](FORTRAN_STYLE.md) — coding style + the GPU (`do concurrent` + OpenACC) programming guide
- [`docs/codebase/INDEX.md`](docs/codebase/INDEX.md) — repo map + run lifecycle
- [`docs/ROADMAP_OCEAN.md`](docs/ROADMAP_OCEAN.md) — ocean dyn-core long-horizon roadmap
- **Per-module / per-procedure docs** — the source `!!` docstrings, rendered by FORD (built on PR to `main`)
