# The ocean dynamical core

This directory is the ocean solver (`sim_type = 'ocean'`, the only regime the repository builds). It is
a hydrostatic, layered ocean model on an Arakawa C-grid with a split-explicit time step. Layer
thickness is the prognostic variable: continuity is a transport equation solved with a piecewise-parabolic
(PPM) scheme, not a constraint enforced by a Poisson solve. Velocities, free surface, temperature, salinity
and any passive tracers advance together, and the vertical coordinate (sigma, z*, z-fixed, isopycnal, ...)
is restored by an ALE remap.

This page is a map for a first-time reader. It does not repeat what the code already says. Every module
and procedure carries a `!!` docstring that FORD renders, and the contracts that cross module boundaries
are written down once, in the docstring of the module that owns them (see the last section).

## The state

All ocean state lives in one object, `ocean_state_t` (`state/rdb_ocean_state.F90`), a composition of
"slots". A slot is a derived type that owns the arrays and workspace of one concern: the barotropic and
multilayer prognostics, the continuity scheme, the Coriolis term, the pressure force, each mixing or
viscosity closure, the sea-ice model, the diagnostics manager, and so on. `ocean_engine_t`
(`rdb_ocean_engine.F90`) owns an `ocean_state_t` plus the setup products that are not state, and runs the
lifecycle (`engine_setup`, `engine_step`, `engine_step_ice`, `engine_step_finalize`, `engine_teardown`).
The driver (`src/driver/rdb_driver.F90`), the C API and the benchmark all go through the engine.

Conventions you need before reading any kernel:

- Layers are bottom-up: `k = 1` is the bed and `k = nz` is the surface layer.
- `u_face_x(i,j,k)` is on the face between cell `(i,j)` and `(i+1,j)`; `v_face_y` is the same in `j`.
- Tracers are stored as content, `hTr = h*c` (thickness times concentration), never as concentration.
- Every data-parallel loop is a `do concurrent`. The GPU build has separate device memory, so every array a
  kernel touches must be mapped (see `FORTRAN_STYLE.md`).

## Directory map

| Path | What is there |
|---|---|
| `rdb_ocean_engine.F90` | Setup, step and teardown shared by the driver, the C API and the benchmark. |
| `rdb_ocean_halo_state.F90`, `rdb_ocean_ghost_poison.F90` | Outer-step ghost exchange of the prognostics; a NaN-ghost debug probe. |
| `state/` | `ocean_state_t` and everything that builds it: metrics, porous barriers, cavity geometry, bipolar grid map, setup, seeding, restart registry, balance, stability audit, status codes. |
| `state/setup/` | The `configure_ocean_*` routines (grid, barotropic split, vmix, lateral closures, boundaries, cavity). Order is load-bearing and is called from `engine_setup`. |
| `state/ocean_state/` | Submodules of `rdb_ocean_state`: initial-condition seeding, topography, restart registry. |
| `dynamics/split_rk2/` | The outer step: `rdb_ocean_dyn.F90` (`ocean_dyn_step_split`, `run_stage_split`; also the legacy unsplit `ocean_dyn_step`) and `rdb_ocean_bt_wide.F90` (the wide-halo barotropic clone). |
| `kernels/barotropic/` | The fast loop: forward-backward substep, barotropic-to-layer coupling, `DTBT`, flux-bounded continuity types. |
| `kernels/continuity_ppm/` | Continuity-PPM mass fluxes, layer-thickness update and windowed tracer advection. |
| `kernels/coriolis_adv/` | PV-conserving Coriolis and vector-invariant advection. |
| `kernels/vertical_advection/` | Vertical tracer and momentum advection. |
| `vcoord/` | The vertical-coordinate state and target-thickness generators (`VCOORD_*`). The remap itself is in `src/ALE/`. |
| `boundary/` | Open, tidal, sponge and periodic boundaries, the tripolar north fold, baroclinic OBC. |
| `forcing/` | Tides and self-attraction, surface-pressure loading, file and bulk-formula forcing, chlorophyll optics. |
| `diag/` | The diagnostics registry and its fills, NetCDF emit, console statistics, conservation budgets, probes. |
| `io/` | Time-varying NetCDF input, restart read/write, z-level T/S initialisation, lat-lon regridding. |

Physics that is shared with other parts of the tree sits outside this directory: the equation of state
(`src/equation_of_state/`), the pressure force (`src/pressure_force/`), the ALE remap and column
reconstruction (`src/ALE/`), tracer schemes (`src/tracer/`), vertical and lateral closures
(`src/parameterizations/`), sea ice (`src/core/ice/`), and the halo layer (`src/comm/`).

## Four design rules

1. **A handle, not a global.** `ocean_state` is passed as an argument. No module does `use ... only:
   the_one_state` or keeps it in a `save` variable. If you need state, take it as a dummy argument.
2. **Each slot owns its memory.** Allocation happens in the slot's `init`, device mapping in its
   `enter_data`, release in `destroy` / `exit_data`. Adding a slot means one `init` line in the parent, not
   edits to a central allocator. Default-off closures (EPBL, kappa-shear, GM, Redi, MEKE, ...) allocate
   nothing unless enabled; the gating rules are in `docs/codebase/CONTRACTS.md`.
3. **Read-mostly across slots.** A kernel reads other slots' fields but does not mutate them. The outer-step
   driver (`ocean_dyn`) is what advances the prognostics. Kernels stay diagnostic-agnostic: the diagnostics
   manager pulls from state through its registry, and a kernel never pushes to it.
4. **MPI only through `pic_mpi_lib`.** No ocean module says `use mpi`. Halo exchanges and reductions go
   through `src/comm/`, which has a single implementation that also runs on one rank. The `no-mpi-in-rdb`
   pre-commit hook enforces it.

## What one outer step does

`engine_step` is followed by `engine_step_ice` and `engine_step_finalize`; the order matters because the
sea-ice block writes surface-flux components that the finalize step assembles. The default scheme is
`pred_corr` (`&ocean_bt_nml split_scheme`); `ssp_rk2` differs in how the two stages are combined (item 5).

1. **Engine pre-step** (`engine_step`): update file-backed inputs and apply forcing (wind, heat, bulk formulae),
   refresh boundary data, the porous-barrier areas and the analytic sponge target, and refresh the
   barotropic `n_inner` if it is derived dynamically.
2. **Per-step refreshes** (top of `ocean_dyn_step_split`): equilibrium tide and self-attraction, the
   surface-pressure seam and `p_top`, mixed-layer-eddy transports, wave speed, isopycnal slopes, GM / VarMix /
   MEKE coefficients, Redi coefficients.
3. **Save** the layer state and tracer content (and the viscous bottom boundary layer, `vdiff_set_viscous_bbl`).
4. **Two stages** of `run_stage_split` (under `pred_corr` the second is the corrector, and its single
   full-step update is the step):
   1. fill ghosts and periodic or fold seams; derive the barotropic state `(η, u_bt, v_bt)` from the layers;
   2. evaluate the slow tendencies: equation of state, Coriolis and advection, pressure force, lateral
      viscosity, bottom, channel and top drag, surface stress;
   3. sum them into `F_slow` and take the thickness-weighted depth mean as the barotropic forcing;
   4. run the barotropic fast loop for `n_inner` substeps (forward-backward in η and `u_bt`), giving
      time-mean transports;
   5. run the continuity chain: PPM mass fluxes, layer-thickness update, tracer and vertical advection,
      tracer horizontal diffusion;
   6. apply the layer velocity tendencies, then the barotropic correction that makes the layer depth mean
      follow the fast loop;
   7. apply open-boundary baroclinic fills and sponges;
   8. surface heat and salt fluxes, shortwave penetration, geothermal heating, ideal-age source;
   9. implicit vertical mixing and friction (KPP or EPBL, kappa-shear, tidal and background mixing);
   10. halo exchange and seam reconciliation.
5. **Combine stages**: `ssp_rk2` averages the two stages; `pred_corr` skips the average.
6. **After the stages**: GM bolus transports, velocity and CFL truncation, drain of the windowed tracer
   advection (if `dt_tracer_advect_ratio > 1`), deferred thermodynamics, then the **ALE remap** of
   layers, tracers and face velocities onto the target coordinate (on remap steps).
7. **Clean-up**: ideal-age surface reset.
8. **Sea ice** (`engine_step_ice`, no-op if off): frazil, EVP dynamics and ice-ocean stress, and at thermo
   cadence ice transport and column thermodynamics.
9. **Finalize** (`engine_step_finalize`): ice-shelf melt, assembly of the surface-flux components into
   `Q_heat` / `Q_salt`, and one diagnostics-manager step.

Restart checkpoints, console statistics and diagnostic output are written by the driver after the step.

## Adding a slot or a physics knob

1. Add a namelist knob that defaults to off, so existing namelists and tests stay bit-identical
   (`docs/howto/add_namelist_knob.md`). Fail loud in `validate_config` for combinations you do not support.
2. Put the arrays in a slot type with `init`, `destroy`, `enter_data`, `exit_data` and `bytes`, and
   an `is_init` flag. If the closure is default-off, allocate only when enabled and include the slot in the
   memory report; latch `enable` before `init` (`docs/codebase/CONTRACTS.md`).
3. Write the kernel as a `do concurrent` over explicit-shape arrays, gated on the knob. Map every array it
   touches; there are no implicit host-device copies.
4. Join the checkpoint through `ocean_state_build_restart_registry` if the slot carries state across steps.
5. Add a unit test, preferably analytical, that covers the new path, plus a bit-identity test for the knob
   off.
6. Update the synthesis docs in the same change: `docs/CLOSURE_MATRIX.md`, `docs/CAPABILITIES_AND_LIMITATIONS.md`,
   and the matching guide under `docs/howto/`. Write the contract in the owning module's `!!` docstring.

## Where to read more

Contracts live in the module that owns them. `docs/codebase/CONTRACTS.md` has the full lookup table and
the rules that have no single owner (indexing, `is_init`, memory gating, one rho0, land state, real mass).

| Topic | Read |
|---|---|
| Vanished layers: invariant I1′, the pool, who may substitute what | header of `src/shared_module_utilities/rdb_vanished_layer.inc` |
| `p_top`, the three pressures, rules for a builder joining the seam | `p_top` in `src/core/rdb_multilayer_state.F90` |
| `eta_forcing` seam: tide + SAL + surface pressure | `src/core/ocean/forcing/rdb_ocean_p_surf.F90` |
| `F_slow` seam, the `bt_rem` accumulator, frozen-reference subtraction | `src/core/ocean/kernels/barotropic/rdb_barotropic_coupling.F90` |
| Cavity datum and the counted-once ice load | `state/rdb_ocean_cavity.F90`, `state/setup/rdb_ocean_setup_cavity.F90` |
| Basal melt, real freshwater mass | `src/parameterizations/vertical/rdb_ocean_cavity_flux.F90` |
| Surface-flux component ownership | `src/parameterizations/vertical/rdb_ocean_surface_flux.F90` |
| Boundary-layer friction velocity (`stress_mag` + `stress_shelf`) | `src/parameterizations/vertical/rdb_ocean_surface_stress.F90` |
| Reference densities | `src/pressure_force/rdb_ocean_pressure_force.F90`, `configure_ocean_reference_density` in `state/setup/rdb_ocean_setup_lateral.F90` |
| Land state | `ocean_state_seed_land_cells` in `state/ocean_state/rdb_ocean_state_seed.F90` |
| Which closures exist and their knobs | `docs/CLOSURE_MATRIX.md`, `docs/generated_nml_knobs.md` |
| What works and what does not | `docs/CAPABILITIES_AND_LIMITATIONS.md` |
| Physics and namelist reference | `docs/REFERENCE.md` |
| Style, `do concurrent` and OpenACC practice | `FORTRAN_STYLE.md` |
