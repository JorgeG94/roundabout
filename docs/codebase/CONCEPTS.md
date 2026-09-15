# Concepts and Acronyms

Quick reference for the acronyms and short conceptual labels used
throughout the Roundabout codebase docs, source comments, and PR
discussion. Grouped by topic. If you find an acronym in the docs that
isn't here, please add it.

## Numerics — schemes and discretisation

| Acronym | Expansion | What it means in Roundabout |
|---------|-----------|--------------------------|
| **CWC** | Consistency With Continuity | The constraint that tracer-flux must agree with the same depth-integrated mass flux the barotropic continuity step used. Equivalently: `sum_k f_layer(k) = mass_flux_face`. Without CWC, uniform tracers don't stay uniform. |
| **PGF** | Pressure Gradient Force | The pressure-gradient term in the layered momentum equation. Selected by `&ocean_pgf_nml form=`: Montgomery, FV-lite, FV-Wright, gprime (2-layer) or FV-MOM6. (Older notes call the baroclinic part "BPG".) |
| **PCM / PLM / PPM** | Piecewise Constant / Linear / Parabolic Method | Increasing-order conservative remap reconstructions used by `rdb_remap_column`. PCM is donor-cell, PLM is minmod-limited, PPM is parabolic with monotonicity (Colella-Woodward 1984). |
| **CFL** | Courant-Friedrichs-Lewy | Stability number `u·dt/dx`. Sets the explicit timestep. |
| **SSP-RK** | Strong Stability Preserving Runge-Kutta | Time integrator. Roundabout uses SSP-RK2 (two-stage). Preserves monotonicity properties of the spatial discretisation. |
| **CG** | Conjugate Gradient | Iterative linear-system solver used for the non-hydrostatic Poisson pressure correction. |
| **PP81** | Pacanowski & Philander (1981) | Richardson-number-dependent vertical mixing scheme. `nu_t = nu_max / (1 + α·Ri)^p + nu_bg`, with p=2 for momentum and p=3 for tracers. |
| **EOS** | Equation Of State | Density formula. Roundabout uses the linear two-tracer EOS `rho = rho_0 + beta_S·(S - S_ref) - alpha_T·(T - T_ref)`. |
| **PGE** | Pressure Gradient Error | Spurious horizontal pressure gradient that arises in sigma coordinates over steep bathymetry, because constant-density surfaces don't align with sigma surfaces. Mitigated by the Jacobian BPG and by ZSTAR / ZSIGMA vertical coordinates. |

## Numerics — Roundabout-specific patterns

| Term | What it means |
|------|---------------|
| **Mode split** | Separating the depth-integrated barotropic step (cheap, fast wave speed) from the per-layer baroclinic distribution (slow). Each SSP-RK stage runs the barotropic update first, then distributes per-layer. |
| **Vanishing layer** | A layer with thickness ≤ `VANISHING_LAYER_TOL` (≈ `h_min`). Occurs in ZSTAR_FULL when the local depth `H < h_bed_ref`. Carries negligible mass; treated specially by the vanishing-layer-gated centred CWC advection (excluded from the per-layer flux split) and by per-layer kernels that gate on `VANISHING_LAYER_TOL` to avoid `1/h` blow-ups. |
| **Vanishing-layer-gated centred CWC** | The per-face per-layer mass-flux split used by the layered tracer advection. For each face, layers where either side is vanishing are excluded from the per-layer flux (`v_xL(k) = 0`); the residual is redistributed across active layers via a per-face symmetric `frac_face_active = h_face_active(k) / sum h_face_active`. Reduces to the original centred CWC for sigma / zsigma / zstar-lite. |
| **Conservative soft clamp** | Two-pass column-conservative tracer limiter in `ml_clamp_tracer`. Pass 1 clips violators to `[tr_min, tr_max]` and accumulates per-cell `excess_above` and `deficit_below`; pass 2 redistributes the net imbalance into spare capacity within the same column. Per-cell tracer total preserved exactly when the column-mean is in bounds. |
| **Outer-shim + flat-impl** | Pattern for kernels that need to read array-of-derived-types fields (e.g. `tracers(it)%hTr`) from inside a `do concurrent`. The outer shim dereferences once on the host and passes flat 3D arrays into the inner impl. NVHPC stdpar can't dereference array-of-derived-types indirection from device code; the shim works around that. |
| **`_impl` duplication** | NVHPC won't inline cross-module `acc routine seq` helpers across translation units. For hot kernels we *duplicate* the helper bodies as `*_impl` copies in the calling kernel module so the compiler sees the body in the same TU and can inline it. Trade-off: parallel implementation drift risk vs. measured 3-4× speedups on flux/extrapolate. |
| **ALE** | Arbitrary Lagrangian-Eulerian | A vertical-grid philosophy where layers can move with the flow during dynamics, with periodic remapping to a fixed target grid. MOM6 / NEMO use full ALE; Roundabout uses a "mode-split + remap" half-step that has some of the same benefits but rescales `dz_old` to fit `h_val`. The rescaling has known interactions with hard tracer clamps; Roundabout's conservative soft clamp absorbs them. |

## Computing — hardware and toolchain

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **GPU / CPU** | Graphics / Central Processing Unit | Roundabout targets GPUs as primary (via `do concurrent` + OpenACC). CPU multicore is supported via `-stdpar=multicore`. |
| **CUDA** | Compute Unified Device Architecture | NVIDIA's GPU programming stack. Used implicitly via NVHPC's stdpar / OpenACC backends. |
| **OpenACC** | Open Accelerators | Directive-based GPU programming. Roundabout uses bare directives (`!$acc enter data`, `!$acc parallel loop`) for data movement and reductions; `do concurrent` covers the rest. Inert as comments to non-OpenACC compilers. |
| **MPI** | Message Passing Interface | Multi-rank distributed memory communication. Used for multi-GPU / multi-node halo exchange and I/O server. Roundabout supports both `mpi_f08` (modern) and `mpi` (legacy) backends. |
| **NVHPC** | NVIDIA HPC SDK | NVIDIA's Fortran compiler suite (`nvfortran`). Required for GPU offload; Roundabout also supports gfortran and ifx for CPU builds. |
| **DGX** | Deep-learning GPU System | NVIDIA's multi-GPU server hardware. Most Roundabout development happens on a 4-V100 DGX. |
| **SoA / AoS** | Structure of Arrays / Array of Structures | Memory layout. Roundabout uses SoA (one allocatable per field, e.g. `h(:,:)`, `hu(:,:)`) for coalesced GPU access. |
| **TU** | Translation Unit | A single source file as seen by the compiler. Cross-TU inlining is the issue that motivates the `_impl` duplication pattern. |
| **IPA** | Interprocedural Analysis | Whole-program analysis that *would* enable cross-TU inlining (`-Mipa=fast,inline` in NVHPC). Slow and brittle; we don't use it — `_impl` duplication is the deliberate alternative. |
| **FP** | Floating Point | IEEE 754 floating-point arithmetic. Working precision is `wp = real64` (double). |
| **ULP** | Units in the Last Place | The smallest representable difference between two FP numbers. Used to describe rounding-level differences ("ULP-level drift"). |

## Domain — physics, geography, reference codes

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **NH** | Non-Hydrostatic | The full vertical-momentum equation including `dw/dt` (not just hydrostatic balance). Roundabout's NH extension uses a CG Poisson solver. Structured-only. |
| **PSU** | Practical Salinity Units | Salinity unit. Open-ocean is ~34-35 PSU; Roundabout uses `S_min = 0`, `S_max = 40` as default conservative bounds. |
| **SSH** | Sea Surface Height | The free-surface elevation `η` above a reference (z=0). `H = h_bed + η`. |
| **vcoord** | Vertical coordinate | Roundabout-specific shorthand. The choice of how to discretise the vertical (sigma, zsigma, zstar, zstar-sigma, zstar-full). |
| **MOM6** | Modular Ocean Model 6 (NOAA/GFDL) | A reference open-source ocean model. Roundabout borrows the full-z* design (`VCOORD_ZSTAR_FULL`) and the conservative tracer limiter pattern from MOM6. |
| **NEMO** | Nucleus for European Modelling of the Ocean | Another reference open-source ocean model. Similar tracer-limiter patterns. |
| **ROMS** | Regional Ocean Modeling System | Reference for the bottom-up `k=1=bed, k=nz=surface` indexing convention and the s-coordinate (Song & Haidvogel) stretching. |
| **SCHISM** | Semi-implicit Cross-scale Hydroscience Integrated System Model | Reference semi-implicit ocean model. Mentioned in capability discussions as the alternative regime to Roundabout's explicit time-stepping. |

## I/O — file formats and conventions

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **NetCDF** | Network Common Data Form | Self-describing binary format used for all Roundabout I/O (output, restart, bathymetry, forcing, gauges). Built on HDF5. |
| **HDF5** | Hierarchical Data Format 5 | The underlying storage layer of NetCDF-4. Has a known thread-safety issue under multi-threaded NVHPC builds — tests that touch NetCDF/HDF5 run with `OMP_NUM_THREADS=1` (set via the CTest `ENVIRONMENT` property; see the gotcha in `CLAUDE.md`). |
| **CF** | Climate and Forecast (metadata conventions) | Standard NetCDF metadata convention. Roundabout's gauge output is CF-1.8-compliant. |

## Documentation, build, and tooling

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **FORD** | FORtran Documenter | Tool that auto-generates HTML API docs from `!!`-prefixed Fortran comments. Roundabout uses `!!` after declarations for FORD discoverability. |
| **PR** | Pull Request | GitHub change-review unit. Most Roundabout features ship in one or more PRs. |
| **WIP** | Work In Progress | A branch / commit that isn't ready to merge yet. |
| **TODO** | "to do" | Inline marker for known follow-up work. Larger TODOs are tracked in the roadmap (`docs/ROADMAP_OCEAN.md`). |
| **CMake** | Cross-platform Make | Roundabout's build system. See the CMake-options table in `CLAUDE.md`. |
| **CI** | Continuous Integration | Automated build + test on every PR. Roundabout has GitHub Actions for the serial / GPU / MPI builds. |
| **MWE / MRE** | Minimal Working / Reproducing Example | Smallest possible code that demonstrates an issue — e.g. a throwaway `do concurrent` compiler-codegen reproducer used to isolate an NVHPC/flang bug. |

## Project glossary

| Term | What it means |
|------|---------------|
| **Roundabout** | The solver. Named for Canberra — the one big Australian city with no ocean — and its roundabouts, which look a lot like eddies. |
| **rdb_** | Source-file prefix for all Roundabout Fortran modules. |
| **`wp`** | Working precision — currently `real64` (`double` in C). All physical state and arithmetic are in `wp`. |
| **`hgrid_t`** | Horizontal-grid metadata type (cell counts, ghost width, dx/dy). |
| **`state_t`** | Top-level solution state (composition of barotropic, multilayer, NH, forcing, BC, decomp sub-objects). |
| **Mode-split** | (See Numerics section above.) |
| **Outer-shim + flat-impl** | (See Numerics section above.) |
| **`_impl` duplication** | (See Numerics section above.) |

---

Cross-references:
- Codebase map + run lifecycle → [INDEX.md](INDEX.md)
- Per-module / per-kernel / solver call-chain detail → the source `!!` docstrings, rendered by FORD (the old `KERNELS.md` / `SOLVER.md` / `DATA_FLOW.md` were retired in favour of this)
- Design contract → [`src/core/ocean/README.md`](../../src/core/ocean/README.md)
- GPU patterns and gotchas → [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md), [CLAUDE.md](../../CLAUDE.md) → "Gotchas"
- Capability matrix and known limitations → [CAPABILITIES_AND_LIMITATIONS.md](../CAPABILITIES_AND_LIMITATIONS.md)
