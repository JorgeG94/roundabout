# Roundabout Codebase Map

GPU-native 3D hydrostatic ocean solver (Arakawa C-grid). The canonical source is
Fortran `do concurrent` + OpenACC (portable across NVHPC GPU/
multicore, gfortran, ifx); an OpenMP-target variant for Intel/AMD GPUs is
generated from that same source by the acc→omp transformer (`tools/acc_to_omp.py`
/ `tools/dc_to_omp.py`, pushed to the `auto/*` branches).

> **Where to read first** depends on what you're working on:
>
> - **New to the codebase / acronyms**: [CONCEPTS.md](CONCEPTS.md)
> - **The dyn-core (Arakawa C-grid + continuity-PPM)**: [`src/core/ocean/README.md`](../../src/core/ocean/README.md) — slot map + design contract
> - **Which closure/scheme is enabled where + tunable knobs**: [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md)
> - **Capabilities + limitations**: [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md)
> - **Physics + namelist reference**: [`docs/REFERENCE.md`](../REFERENCE.md)
> - **Per-module / procedure detail**: the source `!!` docstrings, rendered by **FORD** (built on PR to `main`)

## Directory Structure

```
rdb/
├── app/
│   └── main.F90                        Thin entry — calls driver_run
├── src/
│   ├── core/                           Types, config, decomp, profiler + the
│   │                                   shared state types: rdb_grid (hgrid_t),
│   │                                   rdb_tracer, rdb_barotropic_state
│   │                                   (+ _workstate), rdb_multilayer_state,
│   │                                   rdb_state (register_default_tracers only)
│   │   ├── ocean/                      The C-grid dyn-core (operational, Tier-1)
│   │   │   ├── state/                   rdb_ocean_state (the god state), rdb_ocean_setup,
│   │   │   │                             rdb_ocean_metrics, rdb_ocean_porous, rdb_massless
│   │   │   ├── dynamics/split_rk2/      rdb_ocean_dyn (split-explicit RK2 — the step entry point)
│   │   │   ├── kernels/
│   │   │   │   ├── continuity_ppm/      rdb_continuity, rdb_ocean_min_thickness
│   │   │   │   ├── coriolis_adv/        rdb_coriolis_adv
│   │   │   │   ├── barotropic/          rdb_barotropic_substep, rdb_barotropic_coupling, rdb_bt_cont_type
│   │   │   │   └── vertical_advection/  rdb_ocean_vertical_advection
│   │   │   ├── vcoord/                  rdb_ocean_vcoord (ALE remap dispatch)
│   │   │   ├── diag/                    rdb_ocean_diag (+ _fills/_derived/_netcdf/_mask);
│   │   │   │                             rdb_ocean_budgets (conservation-check PRIMITIVES only)
│   │   │   ├── forcing/                 rdb_ocean_tides, rdb_ocean_tide_astro, rdb_ocean_p_surf,
│   │   │   │                             rdb_ocean_data_forcing
│   │   │   ├── boundary/                rdb_ocean_obc, rdb_ocean_sponge, rdb_ocean_periodic,
│   │   │   │                             rdb_ocean_fold(_apply), rdb_ocean_boundary_{types,data}
│   │   │   └── io/                      rdb_ocean_data_input, rdb_ocean_restart(_io), rdb_ocean_z_init
│   │   └── ice/                        Sea ice (SIS2 port)
│   │       ├── thermo/                  Winton column, enthalpy, optics, mass, frazil, snow
│   │       ├── itd/  transport/  dynamics/   category ITD, transport, C-grid EVP
│   │       └── state/  coupling/        rdb_ice_state, rdb_ice_init, rdb_ice_ocean_coupler
│   │
│   ├── ALE/                            rdb_vcoord + rdb_remap_column + rdb_ocean_remap
│   ├── tracer/                         rdb_recon_weno (PLM/WENO face-swept
│   │                                   primitives), rdb_ocean_ideal_age,
│   │                                   rdb_ocean_pseudo_salt
│   ├── equation_of_state/              rdb_eos (Wright / Roquet / linear)
│   ├── pressure_force/                 rdb_ocean_pressure_force (FV Boussinesq),
│   │                                   rdb_ocean_pgf_reconstruct
│   ├── parameterizations/
│   │   ├── vertical/                    rdb_ocean_vmix (PP81+KPP), rdb_ocean_vdiff,
│   │   │                                rdb_ocean_epbl, rdb_ocean_kappa_shear,
│   │   │                                rdb_ocean_tidal_mixing, rdb_ocean_bottom_drag,
│   │   │                                rdb_ocean_surface_{stress,flux}, rdb_ocean_geothermal
│   │   └── lateral/                     rdb_ocean_lateral_mix (Leith/Smag),
│   │                                    rdb_ocean_horizontal_viscosity, rdb_ocean_hdiff_tracer,
│   │                                    rdb_ocean_{meke,gm,redi,varmix,mle,isopycnal_slopes}
│   ├── framework/                      rdb_scratch_3d (scratch_3d_buffer_t), rdb_safe_math,
│   │                                    rdb_efp, rdb_mem_report
│   ├── driver/                         rdb_driver — full compute-rank lifecycle
│   ├── io/                             rdb_io_netcdf (+ output_rank_filename), rdb_banner,
│   │                                    rdb_bathymetry
│   └── comm/                           MPI facade — single source of truth for MPI
│       ├── mpi/                          Multi-rank backend (talks to pic_mpi_lib only — see no-mpi lint)
│       └── single/                       Single-process stubs
├── tests/                             Flat test_*.F90 — testdrive unit + analytical
│   └── mpi/                            MPI integration (halo, dyn-core)
├── benchmarks/                         bench_ocean — namelist-driven dyn-core throughput
├── validation_examples/ocean/          Canonical / analytical benchmark configs
├── tools/                              Lint scripts, post-processing, OpenMP↔OpenACC translators
├── docs/                               REFERENCE.md (physics + namelist), codebase/ map,
│                                       ROADMAP_OCEAN.md, …
├── cmake/                              Build options, toolchain guards, compiler flags,
│                                    dependency finders
├── CLAUDE.md                           Agent-facing project instructions
├── FORTRAN_STYLE.md                    Coding style guide
└── CMakeLists.txt                      Build system entry point
```

## Code-level reference → FORD

Per-module and per-procedure documentation comes from the `!!` docstrings in the
source, rendered by **FORD** (built on PR to `main`). We deliberately do **not**
hand-maintain a parallel set of markdowns describing what each module/kernel
does — that's the copy that rots (the previous `KERNELS.md` / `CORE.md` /
`API.md` / … did exactly that). The docstring next to the code stays current.

The surviving **curated / cross-cutting** docs — the synthesis FORD can't
generate — are:

| Document | Covers |
|----------|--------|
| [CONCEPTS.md](CONCEPTS.md) | Acronyms (PPM, PGF, CWC, PP81, KPP, ALE, …) + one-line concept definitions. Start here. |
| [`src/core/ocean/README.md`](../../src/core/ocean/README.md) | Design contract: god-state slot map + the four rules. |
| [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md) | Which closure/scheme is enabled + the tunable knobs. |
| [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md) | What the solver can / can't do today. |
| [`docs/REFERENCE.md`](../REFERENCE.md) | Physics, numerics, namelist keys. |
| [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) | Style + the GPU (`do concurrent` + OpenACC) programming guide. |
| [`docs/OPENMP_VARIANT_STATUS.md`](../OPENMP_VARIANT_STATUS.md) | Status of the generated OpenMP variant branches on ROCm — which build, and the flang bug that stops `dc-openmp-target`. |

## Run lifecycle (`app/main.F90` → `driver_run`)

```
comm_env_init                    MPI phase 1
read_config(namelist)
comm_env_setup_roles(...)        MPI phase 2 (split communicators)
decomp_init_from_config          subdomain extents
ocean_state%init_from_config     allocate
load_bathymetry (rank 0) → scatter → fill ghosts
ocean_state_seed_from_cfg        initial conditions
register_default_tracers         S + T (+ opt-in passive tracers)
configure_ocean_*                diag registry, forcing, metrics, tides,
                                 closures, PGF, BT split, BC, sponge,
                                 land mask, wave drag, porous, sea ice
ocean_halo_init → halo exchange  ghost fill on the host
ocean_state_enter_data           map state to the GPU (once)

  TIME LOOP:
    ocean_data_input_update_all  time-interpolated forcing (NetCDF builds)
    ocean_dyn_step_split         one split-explicit RK2 outer step
                                 (or ocean_dyn_step, unsplit)
    ice_* (frazil / EVP / thermo / transport)   when &ocean_ice_nml enable
    diag / restart               at their cadences

ocean_state_exit_data            map back from the GPU (once)
diag close / finalize
```

## Notes on the layout

The mental model:

- **The dyn-core**: `src/core/ocean/` — everything C-grid-specific lives here,
  with `src/core/ice/` alongside it for sea ice
- **Shared physics**: `src/{ALE, tracer, equation_of_state, pressure_force, parameterizations, framework}/`
- **Single source for MPI**: `src/comm/` — never `use mpi` / `use mpi_f08` anywhere
  else; enforced by the `no-mpi-in-rdb` pre-commit hook

The step is reached via `driver_run` → `driver_run_ocean` →
`ocean_dyn_step_split`. The god-state slot map and the "how to pick up a slot"
recipe live in `src/core/ocean/README.md`; remaining gaps (non-hydrostatic on
the C-grid, C-grid MPI halo) are tracked in `docs/ROADMAP_OCEAN.md`.

> The coastal A-grid (HLL/HLLC) path and the unstructured (triangular, KNP)
> backend used to live under `src/core/{coastal,unstructured,multilayer}/` and
> were split into their own repository — nothing in this tree references them.
