# Roundabout — Capabilities and Limitations

A snapshot of what the solver can and cannot do today. Every "cannot" below is something the code does not do yet. It is not a promise that it never will. The roadmap (`ROADMAP_OCEAN.md`) lists the planned work, its scope tags and its validation criteria.

> **The code is the authority.** When this page and the source disagree, the source is right and this page is the bug.

*Last fully reconciled against the codebase: 2026-08-23 (coastal carve-out). The page was restructured on 2026-10-12. Only the statements listed in [section 17](#17-revision-notes) were re-checked against the source then. Everything else is as reconciled on 2026-08-23, plus the dated entries in the [archive](CAPABILITIES_ARCHIVE.md).*

**Scope note.** The `sim_type='coastal'` regime was split out into a separate repository. That regime was the A-grid HLL/HLLC path, the unstructured triangular (KNP central-upwind) backend, the semi-implicit (Casulli) free surface, the non-hydrostatic extension and the C/Python FFI. This tree is ocean-only. Nothing below describes those paths.

## Contents

- [1. What the model does today](#1-what-the-model-does-today)
- [2. Does the model fit my use case?](#2-does-the-model-fit-my-use-case)
- [3. Time stepping and the time-step limit](#3-time-stepping-and-the-time-step-limit)
- [4. Dynamical core](#4-dynamical-core)
- [5. Vertical coordinates](#5-vertical-coordinates)
- [6. Vertical mixing and eddy closures](#6-vertical-mixing-and-eddy-closures)
- [7. Horizontal grids and boundaries](#7-horizontal-grids-and-boundaries)
- [8. Surface forcing and tides](#8-surface-forcing-and-tides)
- [9. Sea ice](#9-sea-ice)
- [10. Ice-shelf cavities](#10-ice-shelf-cavities)
- [11. MPI and domain decomposition](#11-mpi-and-domain-decomposition)
- [12. GPU, compilers and performance](#12-gpu-compilers-and-performance)
- [13. Input, output, restart and initial conditions](#13-input-output-restart-and-initial-conditions)
- [14. Validation and conservation](#14-validation-and-conservation)
- [15. Open problems and the working envelope](#15-open-problems-and-the-working-envelope)
- [16. Data assimilation and machine learning](#16-data-assimilation-and-machine-learning)
- [17. Revision notes](#17-revision-notes)

## 1. What the model does today

- **One regime.** `sim_type='ocean'` is the only accepted value. It is a hydrostatic, Boussinesq ocean on a structured Arakawa C-grid.
- **Dynamics.** Layer thickness is a prognostic variable, moved by continuity-PPM transport. There is no Poisson solve. Momentum is vector-invariant, with a Sadourny PV-flux Coriolis term and a finite-volume pressure gradient.
- **Time stepping.** The step is split-explicit. The fast gravity-wave mode is sub-cycled, so the outer step is limited by advective and baroclinic CFL, not by `dx/√(gH)`. Typical steps are minutes (1200 s in the MOM6 double gyre).
- **Vertical coordinate.** ALE (advance Lagrangian, then remap conservatively). Ten families dispatch through one remap path.
- **Physics.** Wright and linear equations of state; PP81, KPP, EPBL, kappa-shear, tidal and convective vertical mixing; Leith, Smagorinsky, biharmonic, MEKE, GM, Redi and Fox-Kemper closures; bottom drag; a MOM6-style surface-flux set.
- **Grids and boundaries.** Cartesian, spherical, MOM6 supergrid and tripolar (with a north fold). Boundary types are wall, open (Flather), tidal, clamped, Chapman, sponge and periodic.
- **Optional components, all default off.** A port of the SIS2 sea-ice model, and an ice-shelf cavity with basal melt.
- **Parallelism.** One source for CPUs and NVIDIA GPUs (`do concurrent` plus OpenACC). MPI domain decomposition gives results that are bit-identical to the serial run, except the paths listed in [section 11](#not-refused-but-not-bit-identical-across-decompositions).
- **Maturity.** Tier-1 operational since 2026-05-24. The MOM6-reference double gyre runs stable to day 580 on one V100. The model is production-ready for regional hydrostatic configurations within the limits listed here.
- **Main deferrals.** Hollingsworth-Källén Coriolis correction, non-hydrostatic dynamics, the KPP bottom boundary layer, several closures on more than one MPI rank, and restarts across a change of rank count.

### How to read the tables

Each section starts with a table. The columns are feature, status, key knob and limits. Prose below the table covers the essential caveats. Longer notes, measurements, test names and history are in the [archive](CAPABILITIES_ARCHIVE.md). Each section ends with links to its archive entries.

| Status | Meaning |
|---|---|
| **Shipped** | Implemented and tested. Most new-physics knobs default to off, which keeps existing answers bit-identical. |
| **Experimental** | Implemented and tested, but with a known cost to the answer or a thin validation record. |
| **Refused** | The configure step stops the run with a message naming the knob. This is called "fail-loud". |
| **Not implemented** | The code does not do this today. |

### Words used on this page

| Term | Meaning |
|---|---|
| C-grid | Arakawa C staggering. `eta` and tracers sit at cell centres, `u` and `v` on cell faces, vorticity at corners. |
| ALE | Arbitrary Lagrangian-Eulerian. Layers move with the flow for a step, then are remapped onto a target grid. |
| PGF | Pressure-gradient force. |
| BT, barotropic | The depth-averaged, fast gravity-wave mode. It is sub-cycled inside each outer step. |
| Bit-identical | Equal in every bit of every owned value. "Byte-identical" means the same for output files. |
| Fail-loud | A configure-time check that refuses the run and says why, instead of running something unvalidated. |
| Filler, vanished layer | A layer thinner than `H_VANISHED` (1.5e-4 m). It exists only to keep the layer count fixed. |
| Closed faces | A static per-layer wall where a z-level layer is a filler on either side of a velocity face. |
| Cavity | An ice-shelf cavity: a static ice draft over the ocean. |
| rx0 | Slope stiffness of a terrain-following grid, `\|ΔH\|/(H_a+H_b)` across one face. |

More acronyms (CWC, PGF, PP81, ALE, EPBL, MEKE) and concept definitions (vanishing layer, conservative soft clamp, mode split) are in [`codebase/CONCEPTS.md`](codebase/CONCEPTS.md).

### Where to look next

| You want | Read |
|---|---|
| Which closure is enabled, and its knobs | [`CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md) |
| Physics and namelist reference | [`REFERENCE.md`](REFERENCE.md) |
| The planned and deferred work | `ROADMAP_OCEAN.md` |
| The design contract for the ocean state | [`../src/core/ocean/README.md`](../src/core/ocean/README.md) |
| Concepts and acronyms | [`codebase/CONCEPTS.md`](codebase/CONCEPTS.md) |

## 2. Does the model fit my use case?

| I want to... | Fit? | Notes |
|---|---|---|
| Run a regional or basin hydrostatic ocean at eddy-resolving scale | Yes | Tasman 2 km (KPP, ALE, sponges) costs about 34 s per simulated day on one V100. |
| Run a global 1-degree ocean on a tripolar grid | Yes | `validation_examples/ocean/global_1deg/` runs on 1, 2 and 4 ranks with identical state. |
| Run on NVIDIA GPUs | Yes | NVHPC is the primary target. GPU tests run locally, because CI has no GPU runner (see [section 12](#12-gpu-compilers-and-performance)). |
| Run on AMD or Intel GPUs | Not supported | Not a supported or CI-tested target. GPU offload goes through NVHPC (`-stdpar=gpu` plus OpenACC), and there is no HIP or SYCL backend. The code is also built with Intel and AMD toolchains to find compiler bugs, but that is not a GPU port. |
| Study an ice-shelf cavity with basal melt | Partly | Single rank only. The ISOMIP+ files are configuration deliverables, not validated results. See [section 10](#10-ice-shelf-cavities). |
| Include sea ice | Partly | A dynamical core and column model. No ridging, melt ponds or lateral melt. See [section 9](#9-sea-ice). |
| Run an ice-coupled global climate model | No | The multi-rank sea-ice fold is not wired, and the Hollingsworth-Källén Coriolis correction is not implemented. |
| Study coastal flooding, dam breaks, bores, estuaries | No | That work lives in the separate coastal repository. |
| Resolve non-hydrostatic or phase-resolving wave physics | No | The core is hydrostatic. There is no SWAN or WW3 coupling. |
| Use a semi-implicit (SCHISM or Casulli style) free surface | No | The gravity-wave mode is sub-cycled explicitly. There is no implicit `eta` solve. |
| Model biogeochemistry or sediment | No | No source terms and no FABM. Bathymetry is fixed. Passive tracers do work. |
| Assimilate data or plug in a machine-learning surrogate | No | Only relaxation (sponge, nudging) is available. |

### Standing physics gaps

| Gap | Status | What exists instead |
|---|---|---|
| GLS, k-ω and MY2.5 turbulence closures | Not implemented | PP81, KPP, EPBL, kappa-shear, tidal mixing, convective adjustment, double diffusion, two background schemes |
| TEOS-10 GSW equation of state | Not implemented | Wright (1997) and a Roquet et al. (2015) polynomial |
| Sediment transport, biogeochemistry, vegetation drag | Not implemented | Surface and bottom flux set only |
| Online bulk aerodynamic fluxes | Standalone kernel only | Large and Yeager (2004, 2009) is in `rdb_ocean_bulk_flux` but is wired to no forcing source or namelist. Forcing is pre-computed and read by `&ocean_dataovr_nml`. COARE and Zeng are not implemented. |
| Complete sea-ice model | Partial | See [section 9](#9-sea-ice). |

**Full notes in the archive:** [physics gaps and non-goals, in the original wording](CAPABILITIES_ARCHIVE.md#a1)

## 3. Time stepping and the time-step limit

The ocean path is **not** bound by the gravity-wave CFL. The barotropic (gravity-wave) mode is sub-cycled in a fast forward-backward loop. The outer step is limited by advective and baroclinic CFL instead, so it is minutes long. The substep count comes from `&ocean_bt_nml auto_n_inner`, which derives it once at configure from the gravity-wave CFL of each wet cell.

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Predictor-corrector outer step (MOM6 style) | **Shipped, default since 2026-09-14** | `&ocean_bt_nml split_scheme = "pred_corr"` | Refused with `eulerian_z`, wet/dry and `dt_tracer_advect_ratio > 1`. Slows the resting-state growth about 36 times but does not remove it. |
| Two-stage SSP-RK2 outer step | **Experimental** (fully supported and tested) | `split_scheme = "ssp_rk2"` | Makes energy out of a motionless stratified state. Needed for `eulerian_z`, wet/dry and windowed tracer advection. |
| Fixed outer step | Shipped | `&time_nml dt_fixed > 0` | Required. There is no adaptive-CFL helper on this path. |
| Barotropic sub-cycling | Shipped | `&ocean_bt_nml auto_n_inner`, `n_inner` | `n_inner` is set once at configure and is never re-derived. CFL truncation is counted, not controlled. |
| Backward-Euler vertical viscosity and diffusion | Shipped | `&ocean_vdiff_nml implicit_drag`, `implicit_stress` | Stable for any layer thickness. Bottom drag and wind stress can be folded into the same solve. |
| Viscous-remnant chain (MOM6 BBL glue plus barotropic damping) | Shipped, **on by default** | `&ocean_vdiff_nml hvel_mom6`, `bbl_glue`; `&ocean_bt_nml visc_rem_chain` | Refused with `substep_drag` and with `bt_halo > 0`. Setting all three knobs to `.false.` restores the old path. |
| Per-layer retro-correction of the baroclinic pressure | Shipped, opt-in | `&ocean_bt_nml correction_bc_pgf` | Refused unless the PGF is `fv_mom6`. |
| Barotropic correction weights | Uniform increment only | – | The h-weighted fold (`correction_h_weighted`) and `correction_visc_rem` are retired and refused. `accel_visc_rem` is retired and refused. |
| Implicit (semi-implicit) free surface | Not implemented | – | By design. The fast mode is sub-cycled explicitly. |

### Which scheme to use

Use the default `pred_corr` unless your case needs `eulerian_z`, wet/dry or windowed tracer advection. Those three need `ssp_rk2`, and `validate_config` tells you so. Some shipped namelists pin `ssp_rk2` for that reason.

`ssp_rk2` is labelled experimental because of what it costs, not because the code is shaky. Its two-stage average amplifies an internal gravity wave by `√(1+(ω·dt)⁴/4)` every step. So it makes energy out of a motionless stratified ocean. On `resting_stratified_channel.nml` the energy `En` reaches 2.992E-05 m²/s² by day 25 under `ssp_rk2`. Under `pred_corr` the same file holds 1.739E-09, which is 17 000 times less. A forced, energetic, viscous run sits far above this noise floor and never notices it. A quiescent, weakly damped or long spin-up run does notice it, and there the manufactured energy is the signal.

`pred_corr` is not perfect either. It slows the resting-state growth from a 2.5-day e-folding to an 83-day e-folding. It does not remove it. So `resting_stratified_channel` stays a scoped XFAIL on the settle gate.

Both schemes wrap the same nonlinear barotropic fast loop. The stability suite runs an `ssp_rk2` twin of every case whose namelist does not pin a scheme, so neither branch can rot.

### Other time-step facts

- The default barotropic damping is `&ocean_bt_nml bebt = 0.1` (MOM6 parity). The continuity renormalisation uses the donor-flip-consistent flux (`&ocean_continuity_nml renorm_consistent_flux = .true.`). Together they close the stepped-terrain barotropic grid-scale mode (FINDING B in [section 5](#5-vertical-coordinates)).
- Coriolis is a closed-form rotation, so it adds no time-step limit.
- This is a sub-cycled fast mode, not a semi-implicit free surface. The gravity-wave term is not treated implicitly in the SCHISM or Casulli sense.

**Full notes in the archive:**

- [measurements and history of the two schemes](CAPABILITIES_ARCHIVE.md#a2)
- [the split algorithm, the barotropic correction and the viscous-remnant chain (PR-1 to PR-3)](CAPABILITIES_ARCHIVE.md#a3)

## 4. Dynamical core

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Layer-thickness transport (continuity-PPM) | Shipped | `&ocean_continuity_nml` | No Poisson or FFT constraint. |
| Pressure gradient: Montgomery | **Shipped, default** | `&ocean_pgf_nml form = "mont"` | Valid over slopes and every vertical coordinate. Equals `fv_lite` when two columns have equal layer thicknesses. |
| Pressure gradient: FV-lite, FV-Wright, FV-MOM6 | Shipped | `form = "fv_lite"`, `"fv_wright"`, `"fv_mom6"` | Over a slope with unequal thicknesses no form is exact. That residual is the open sigma-coordinate PGF error. |
| Pressure gradient: reduced gravity | Shipped, 2 layers only | `form = "gprime"` | Refused unless `nz_layers = 2`. |
| Equation of state | Shipped | `&ocean_eos_nml eos` | Wright (1997) for production, a linear EOS for validation, and a Roquet et al. (2015) polynomial. No TEOS-10 GSW. |
| Freezing point | Shipped | `&ocean_eos_nml tfreeze_set` | Linear liquidus with two named coefficient sets: `"seaice"` (default) and `"isomip"`. A nonlinear liquidus is not implemented. |
| Coriolis | Shipped | `&ocean_coriolis_nml form`, `pv_adv_scheme` | `sadourny` (default), `sadourny_energy` and `sadourny_hk`. The full Hollingsworth-Källén refinement is deferred. `weno5` and `weno7` PV interpolation need `nghost` of 4 and 5. |
| Barotropic wave drag | Shipped, opt-in | `&ocean_bt_nml wave_drag` | `form = "file"` is refused. `"roughness_proxy"` is a placeholder for a true subgrid roughness map. |
| Bottom drag | Shipped | `&ocean_bdrag_nml form`, `hbbl` | `linear` or `quadratic` (default). Acts through the vertical viscosity by default. **No drag acts until a coefficient is set:** `cd`, `r` and `hbbl` all default to 0 (`form` defaults to `"quadratic"`, which with `cd = 0` is off). |
| Side-wall (channel) drag | Shipped, opt-in | `&ocean_bdrag_nml channel_drag` | The barotropic substep does not feel it during the fast loop. |
| Dynamic wetting and drying | Shipped, narrow envelope | `&ocean_wetdry_nml enable` | Sigma coordinate only, single rank, split solver and `ppm_limit_pos` required. |
| ALE remap | Shipped | `&vcoord_nml remap_method` | `ppm` (default), `ppm_h4`, `pqm`. Runs at the thermo cadence. |
| Vertical advection | Shipped | – | Eulerian. `w` is diagnosed from continuity. |
| Along-coordinate tracer Laplacian | Shipped, opt-in | `&ocean_hdiff_nml kappa_h` | Diffuses along model surfaces, not neutral surfaces. `&ocean_redi_nml` is the neutral path. |
| Passive tracers | Shipped | `&ocean_tracers_nml enable_ideal_age`, `enable_pseudo_salt` | No namelist-declared arbitrary dye, CFC or BGC package yet. |
| Configure-time stability audit | Shipped | – | Aborts or warns with named numbers. See the archive entry for this section. |
| Non-hydrostatic pressure | Not implemented | – | Hydrostatic only. |

### Notes on the choices

- **Pressure gradient.** Start with the default `mont`. Use `fv_mom6` with `reconstruct_for_pressure` over sloping topography. `correction_bc_pgf` and `p_top_in_bc` also need `fv_mom6`.
- **Bottom drag.** By default the drag acts on the layers through the vertical viscosity. This is the MOM6 BOTTOMDRAGLAW glue, and the barotropic fast loop feels the same friction. It sits on the first live layer above the bed (`k_bot`), which matters under `z_fixed`. Two limits remain. Under `zstar` and `zstar_full`, the bed-only mode (`hbbl = 0`) drags layer 1, which can be an inert filler, so set `hbbl > 0` there. And `implicit_drag` with `hbbl > 0` is refused.
- **Wet/dry.** The envelope is fail-loud. It is mutually exclusive with BT_cont, upstream-h, shortwave penetration, restoring and geothermal. One v1 gap remains: the surface-flux dry gate is binary, so a column held wet in the hysteresis band gets the full heat and salt flux over its thin top layer.
- **Remap.** PQM falls back to PPM silently when a column has fewer than 5 layers. Horizontal tracer advection is not yet on the coarser thermo cadence unless you use windowed tracer advection (`&ocean_vmix_nml dt_tracer_advect_ratio`, single rank).

**Full notes in the archive:**

- [pressure force, equation of state and Coriolis](CAPABILITIES_ARCHIVE.md#a4)
- [barotropic wave drag](CAPABILITIES_ARCHIVE.md#a5)
- [bottom drag and side-wall drag](CAPABILITIES_ARCHIVE.md#a6)
- [dynamic wetting and drying](CAPABILITIES_ARCHIVE.md#a7)
- [vertical advection, ALE remap, tracer diffusion, stability audit and passive tracers](CAPABILITIES_ARCHIVE.md#a8)

## 5. Vertical coordinates

Ten families dispatch through one ALE remap path. The enum of record is `src/core/rdb_constants.F90`. Choose with `&vcoord_nml vcoord_type`.

| Family | Status | `vcoord_type` | What to know |
|---|---|---|---|
| Lagrangian (isopycnal) | Shipped | `lagrangian` | The target is the live thickness, so the remap is a no-op. |
| Eulerian z | Shipped | `eulerian_z` | A stretched sigma with `eta` ignored. It is not a geopotential coordinate despite the name. Needs `ssp_rk2`. |
| Sigma | Shipped, the reference family | `sigma` | Terrain-following. Has the open sigma-coordinate PGF error over steep slopes. Required for wet/dry. |
| Sigma-to-z blend | **Refused** | `zsigma` | Its deep branch reads depths in metres, but the only writer of that table is dimensionless. |
| z\* (MOM6 z\*) | Shipped (new in 2026-10) | `zstar` | Dilates the `z_fixed` profile with the free surface. Refused under a cavity. Its rest-state envelope is not yet re-pinned. |
| z\* full | Shipped, slated for retirement | `zstar_full` | Terrain-following coarse zone. The sigma PGF error stays on its open faces. |
| z\*-sigma | Shipped | `zstar_sigma` | With the default uniform table it is numerically the same as sigma. 21 shipped namelists select it. |
| Fixed z | Shipped | `z_fixed` | Stretched nominal profiles via `z_fixed_profile`. The only z-like family allowed under an ice-shelf cavity. |
| Density (isopycnal) | Shipped, validation-grade | `rho` | Weakly stratified columns collapse. |
| HYCOM hybrid | Shipped, the production hybrid | `hycom` | z\* floor is in metres. On stepped topography it hits the open-staircase PGF defect (scoped XFAIL). |

### Rules to know before choosing

- **`z_fixed` over a stepped bed needs closed faces.** `&vcoord_nml zfixed_closed_faces` defaults to `.false.`. With it off, `z_fixed` over a bed whose columns fall in different nominal layers is **refused** at configure. The open-staircase pressure gradient is not a slow drift. It blew the 1° Southern Ocean up to about 250 m/s by step 3. A flat or step-free bed is still accepted.
- **What closed faces do.** Where a layer is a filler on either side of a velocity face, the face becomes a free-slip wall for that layer. The mask is static. It is built once at configure from the `eta = 0` target. It is accepted on `z_fixed`, `zstar` and `zstar_full`, and it composes with porous barriers.
- **Closed faces are refused with** wet/dry, `&ocean_bt_nml bt_halo > 0`, `&ocean_hvisc_nml stress_tensor` (and `kh_aniso`), `&ocean_bt_nml correction_bc_pgf`, and `pred_corr` with `&ocean_bc_nml mask_wall_velocity = .false.`.
- **Filler thickness.** `zstar_h_min` must stay at or below `H_VANISHED` on `zstar_full` and `z_fixed`. A larger value is refused.
- **`zsigma`** returns when its depth table is filled in metres.
- **Under an ice-shelf cavity** only `sigma` and `z_fixed` are accepted. See [section 10](#10-ice-shelf-cavities).

### Which family rests on which geometry

The rest-state matrix (`tests/regression/vcoord_matrix.py`) runs a motionless, stratified, f-plane ocean on eleven geometries. The table gives the largest slope stiffness `rx0` at which every viscous cell of the family passes every gate on both toolchains. Outside its envelope a family is not claimed to rest.

| Family | Viscous-leg envelope (v0.1.0) | Inviscid leg passes to | What caps it |
|---|---|---|---|
| `sigma` | rx0 up to 0.8 (top of the ladder) | 0.078 | Nothing measured. |
| `zstar` (MOM6 z\*, closed faces) | Not yet re-pinned | – | Re-measure with `vcoord_matrix_pin.py`. |
| `zstar_sigma`, `zstar_full` | rx0 up to 0.8 | 0.1 | Nothing measured. |
| `eulerian_z` (pins `ssp_rk2`) | rx0 up to 0.8 | Flat bed only | Nothing measured. |
| `lagrangian` | rx0 up to 0.8 | 0.8 | Nothing measured. |
| `z_fixed` (closed faces) | Flat bed only | Flat bed only | Salinity overshoot at the regrid where layers vanish. |
| `rho`, `hycom` | Flat bed only | Flat bed only | FINDING A (viscous leg) and the rest-state mode (inviscid leg). |

The inviscid leg grows a mode on the stiffness ladder. This is documented expected behaviour, because MOM6 grows the same mode.

**Full notes in the archive:**

- [the ten families and how each places its interfaces](CAPABILITIES_ARCHIVE.md#a9)
- [HYCOM's z* floor in metres](CAPABILITIES_ARCHIVE.md#a10)
- [each family under a displaced column top (ice lid)](CAPABILITIES_ARCHIVE.md#a11)
- [the rest-state matrix, the rx0 envelopes and FINDINGS A, B, C](CAPABILITIES_ARCHIVE.md#a12)
- [stretched z_fixed profiles, z_fixed under a shelf, and what it buys](CAPABILITIES_ARCHIVE.md#a13)
- [partial-step face closure (zfixed_closed_faces) on z_fixed, zstar_full and zstar](CAPABILITIES_ARCHIVE.md#a14)
- [the first-live-layer indices k_top and k_bot](CAPABILITIES_ARCHIVE.md#a15)

## 6. Vertical mixing and eddy closures

Every vertical closure adds into the interface diffusivities `kv` and `kt`. `ks` is derived from `kt` afterwards, and `vmix_assemble` applies the floors, ceilings and smoothing once. This table lists what each closure can and cannot do. The per-closure knob list and test names are in [`CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md).

### Vertical mixing

| Closure | Status | Key knob | Limits |
|---|---|---|---|
| PP81 interior (Richardson number) | Shipped | `&ocean_vmix_nml` | Uses one global reference pressure, so thermobaricity is neglected. |
| KPP boundary layer | **Shipped, on by default** | `&ocean_vmix_nml use_kpp` | Includes the `V_t²` term and the non-local transport. The bottom boundary layer is not implemented. The depth does not respond to a stabilising surface buoyancy flux. Excluded with EPBL. Refused with `z_fixed` under a cavity. |
| EPBL (Reichl and Hallberg 2018) | Shipped | `&ocean_epbl_nml enable` | Mutually exclusive with KPP. Includes Langmuir enhancement (`use_lt`) and a penetrating-shortwave TKE charge (`epbl_sw_ctke`, on by default). |
| Kappa-shear (JHL08) | Shipped | `&ocean_kappa_shear_nml enable` | Interior closure, additive with KPP and EPBL. No massless-layer merging yet. Optional vertex form (`at_vertex`). |
| Tidal mixing (St Laurent and Simmons) | Shipped | `&ocean_tidal_mixing_nml enable` | Bottom-intensified, additive. |
| Convective adjustment | Shipped | `&ocean_conv_nml enable` | Writes `kv` and `kt` only, never `ks`, so salt does not convect. Two v1 limits against MOM6. |
| Background mixing | Shipped | `&ocean_vmix_nml bkgnd_profile` **xor** `bkgnd_henyey` | The two schemes are mutually exclusive. Both default off. |
| Double diffusion | Shipped | `&ocean_ddiff_nml` | Makes `ks` differ from `kt`. |
| First-baroclinic wave speed and Rossby radius | Shipped, **diagnostic** | `&ocean_wavespeed_nml enable` | Feeds resolution-aware GM, Redi and MEKE scaling. Not called by the unsplit path. |

### Lateral closures and eddies

| Closure | Status | Key knob | Limits |
|---|---|---|---|
| Leith, Smagorinsky KH and AH, Leith-biharmonic, constant biharmonic, constant `nu_h` | Shipped | `&ocean_hvisc_nml` | Flow-aware closures give almost no viscosity in a fluid at rest. A quiescent run needs a constant floor. |
| MEKE, GM, Redi, VarMix, isopycnal slopes | Shipped | `&ocean_meke_nml`, `&ocean_gm_nml`, `&ocean_redi_nml`, `&ocean_varmix_nml`, `&ocean_slopes_nml` | Not detailed on this page. See [`CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md). |
| Fox-Kemper mixed-layer eddies | Shipped, opt-in | `&ocean_foxkemper_nml enable` | Needs EPBL to supply the mixed-layer depth. Three refinements are deferred. |
| Along-coordinate tracer Laplacian | Shipped, opt-in | `&ocean_hdiff_nml kappa_h` | See [section 4](#4-dynamical-core). |

### Notes

- **Where each closure sits.** KPP is the default boundary layer. Turn it off (`use_kpp = .false.`) when you enable EPBL.
- **Convective adjustment** raises `kt` to at least `kd_conv` (default 1.0 m²/s) where `N² < n2_thresh`, strictly below the active boundary layer. That large value is safe only because the vertical solves are backward-Euler.
- **Quiescent runs need a constant viscosity floor** over slopes. See [section 15](#15-open-problems-and-the-working-envelope).
- **Alpha and beta in KPP.** KPP's `B_0` and the double-diffusion density ratio use the scalar `alpha_T` and `beta_S` by default. Set `&ocean_vmix_nml buoyancy_coeffs = "eos"` to take them from the active EOS. The two are identical under a linear EOS.
- **Leith-biharmonic and constant biharmonic** closures are implemented and dispatched, along with Smagorinsky KH and AH. An older note that they were enum-only (`LMIX_SMAGORINSKY`, `LMIX_LEITH_BIHARM`, `LMIX_BIHARMONIC`) is stale.
- **The `alpha_T` default.** The historical default `alpha_T = 1.7e-4` is numerically the fractional coefficient, so it hides a missing `1/ρ₀`. It is not a physical configuration for a run that cares about surface buoyancy forcing.

**Full notes in the archive:**

- [lateral closures](CAPABILITIES_ARCHIVE.md#a16)
- [vertical mixing (PP81, KPP, background)](CAPABILITIES_ARCHIVE.md#a17)
- [convective adjustment](CAPABILITIES_ARCHIVE.md#a18)
- [EPBL, including Langmuir enhancement and penetrating shortwave](CAPABILITIES_ARCHIVE.md#a19)
- [kappa-shear](CAPABILITIES_ARCHIVE.md#a20)
- [first-baroclinic wave speed and Rossby radius](CAPABILITIES_ARCHIVE.md#a21)
- [Fox-Kemper mixed-layer eddy restratification](CAPABILITIES_ARCHIVE.md#a22)
- [boundary-layer scheme behaviour under ice (KPP and EPBL)](CAPABILITIES_ARCHIVE.md#a23)

## 7. Horizontal grids and boundaries

### Grids

Choose with `&ocean_grid_nml grid_config`. All kernels read full 2-D metric arrays, so one code path serves every grid.

| Grid | Status | Key knob | Limits |
|---|---|---|---|
| Cartesian | Shipped, default | `"cartesian"` | The bit-identity reference. |
| Spherical lon-lat sector | Shipped | `"spherical"` | – |
| MOM6 supergrid (mosaic) | Shipped | `"supergrid"` | The reader is windowed per rank. The grid rotation `angle_dx` is read, but nothing uses it yet. Vector forcing files must be pre-rotated. |
| Tripolar (Murray 1996 cap plus lon-lat) | Shipped | `"tripolar"`, `north = "tripolar_fold"` | Needs periodic west and east. Decomposes bit-identically. Output in the cap shows model-frame velocities. The analytic tripolar grid cannot be written out and read back as a supergrid. |
| Interior land masking | Shipped | – | Static free-slip walls by zeroing metrics. |
| Porous barriers (Adcroft 2013) | Shipped, opt-in | `&ocean_porous_nml` | Refused with `bt_halo > 0` and with wet/dry. |

### Boundary conditions

Choose per edge with `&ocean_bc_nml`.

| Type | Status | Key knob | Limits |
|---|---|---|---|
| Wall | Shipped, default | – | – |
| Open (Flather) | Shipped | `flather_form`, per-edge `*_ext_u/v` | Barotropic Flather plus a per-layer zero-gradient baroclinic anomaly. The `"full"` form supports an exterior velocity. |
| Tidal | Shipped | `obc_tidal_nodal` | Multi-constituent `eta`. |
| Clamped | Shipped | – | Dirichlet `eta`, `u`, `v` and per-tracer inflow values. |
| Chapman | Shipped | – | **Refused on more than one rank.** |
| Sponge | Shipped | `sponge_relax_tracers`, `&ocean_sponge_nml` | A closed wall at the outer face plus an interior band of momentum decay and optional tracer relaxation. |
| Periodic | Shipped | – | Needs `nghost >= 3` and paired edges. The seam is bit-exact. |
| Radiation and nudging | Shipped, opt-in | `radiation_scheme = "orlanski"`, `nudge_tau_in/out`, `res_lscale_in/out` | Per-layer implicit-upwind Orlanski radiation on the normal velocity, with a running-mean phase speed and an `rx_max` clip. The default `radiation_scheme = "anomaly"` is the Flather mean plus a zero-gradient baroclinic anomaly. `rx` and `u_prev` are not restart state, so a restart starts the scheme cold. |
| Nested | Reserved | – | Behaves as open until the parent-data orchestrator exists. |
| Inflow, discharge | Reserved | – | Tags only. |
| File-backed or tidal-table boundary data | Not implemented | – | Only the constant backend ships. The time-varying reader exists but does not yet deliver boundary segments. |

Other boundary facts:

- The sponge target can be an analytic profile that is independent of the initial condition (`&ocean_sponge_nml target_source = "linear_z"`). The default `"ic"` snapshot is taken on the seed layers. That is only right where the seed already sits on the running coordinate.
- Boundary-corner relative vorticity is zeroed at every non-periodic edge.
- Bathymetry from any source is periodic-wrapped and north-folded before any field is derived from it.
- **Orlanski radiation ships, as an opt-in.** Per-layer implicit-upwind Orlanski radiation is `RAD_ORLANSKI` in `rdb_ocean_obc_baroclinic.F90`, chosen with `radiation_scheme = "orlanski"`. An older "deferred" entry in the archive is superseded (C15 in [section 17](#17-revision-notes)). What is still missing: a per-cell adaptive Orlanski for the barotropic `eta` (the Chapman edge uses one edge-mean value), restart state for `rx` and `u_prev`, and file-backed boundary data.

**Full notes in the archive:**

- [boundary-condition dispatch](CAPABILITIES_ARCHIVE.md#a24)
- [horizontal grids, the tripolar north fold, cap poles and the supergrid reader](CAPABILITIES_ARCHIVE.md#a25)
- [interior-sponge targets](CAPABILITIES_ARCHIVE.md#a26)

## 8. Surface forcing and tides

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Wind stress | Shipped | `wind_config`, `&ocean_surface_*` | Analytic profiles or a file. MOM6 `DIRECT_STRESS` spreads it over `hmix_stress` metres. |
| Surface heat and salt flux (scalar) | Shipped | `&ocean_thermo_nml q_heat`, `q_salt` | Stamped into the top layer (`k = nz`). |
| Surface-flux component set | Shipped, opt-in | `&ocean_forcing_nml enable_components` | Mass fluxes carry enthalpy and salt bookkeeping but do **not** change column mass. `salt_flux` stays virtual. |
| Shortwave penetration | Shipped | `&ocean_thermo_nml sw_pen_frac`, `sw_source`, `sw_transmission` | Two-band transmission shared by the deposit and both boundary-layer schemes. `sw_source = "q_sw"` needs `enable_components`. Inert at `sw_pen_frac = 0`. |
| Surface buoyancy restoring | Shipped, opt-in | `&ocean_restore_nml` | – |
| Time-varying NetCDF forcing | Shipped, opt-in | `&ocean_dataovr_nml enable` | Six tags: `tau_x`, `tau_y`, `heat`, `salt`, `evap`, `lprec`. Files must already be on the model grid. One time mode is shared by all tags. |
| Atmospheric pressure loading (inverse barometer) | Shipped, opt-in | `&ocean_psurf_nml enable` | The version 1 fill is a uniform scalar, which is inert. Needs `enable_components` and the split solver. Excludes `bt_halo > 0`. |
| `p_surf` in the in-situ EOS pressure | Shipped, opt-in | `&ocean_psurf_nml in_eos` | The ice-shelf seam. See [section 10](#10-ice-shelf-cavities). |
| Equilibrium tide and scalar SAL | Shipped | `&ocean_tides_nml` | Scalar SAL is refused under a cavity. |
| Boundary-tide nodal correction | Shipped | `&ocean_bc_nml obc_tidal_nodal` | – |
| Geothermal bottom heat | Shipped, opt-in | `&ocean_geothermal_nml enable`, `q_geo` | Constant flux, heat only. |
| Bulk aerodynamic fluxes (Large and Yeager) | Standalone kernel | – | Not wired to a forcing source or a namelist. |
| In-model lat-lon regridding weights | Library only | – | `rdb_regrid_latlon`: bilinear and first-order conservative, onto spherical, supergrid or tripolar grids; decomposition bit-identical. Not wired to a reader or a namelist. No `patch` variant. |
| File-driven `p_surf` (reanalysis pressure) | Not implemented | – | `p_surf` is not one of the six tags. Adding it is a `register_tag` entry plus a destination slot. |
| Sea-ice mass loading into `p_surf` | Not implemented | – | The `eta_ib` seam is in place for it. |
| Real freshwater mass fluxes | Not implemented | – | Mass fluxes are stored, not applied. Meltwater under a shelf has an opt-in mass form. |

### Notes

- **The dataovr tags.** `tau_x` and `tau_y` feed the C-grid face stress. The stress file must carry `nx_phys+1` x-face values and `ny_phys+1` y-face values. `heat` and `salt` go to `Q_heat` and `Q_salt` directly, or to the `heat_added` and `salt_flux` components when `enable_components` is on. `evap` and `lprec` exist only in the component set.
- **No interpolation.** The reader does no horizontal interpolation, no vertical remap of a source z axis and no calendar. A per-tag time mode is the known follow-up.
- **Sea ice and wind.** The sea-ice atmospheric stress snapshot is taken once at configure, so it does not follow a time-varying wind.
- **Inverse barometer.** An atmospheric high depresses sea level about 1 cm per hPa. The term lands on the barotropic `eta_forcing` seam shared with the tide and SAL. Only the gradient of `p_surf` is physical, so a uniform load changes nothing. A genuine load needs `p_surf` wired up as a dataovr tag.
- **Tides.** Internal-tide drag on the barotropic mode ships separately, as `&ocean_bt_nml wave_drag` (see [section 4](#4-dynamical-core)).

**Full notes in the archive:**

- [wind, heat, salt, tide and pressure-loading entries in the original wording](CAPABILITIES_ARCHIVE.md#a27)
- [the time-varying NetCDF reader and file-backed surface forcing](CAPABILITIES_ARCHIVE.md#a28)
- [surface stress, heat and salt flux, the flux component set and geothermal heat](CAPABILITIES_ARCHIVE.md#a29)
- [what is not shipped in forcing, and the tides deferral note](CAPABILITIES_ARCHIVE.md#a30)

## 9. Sea ice

**This is a high-quality sea-ice dynamical core and column model. It is not yet a complete sea-ice model.** It is a port of SIS2 under `src/core/ice/`, enabled with `&ocean_ice_nml enable` (default off, so byte-identical). The physics that is there can be trusted. What is missing is most of the ice mass budget: no ridging, no melt ponds and no lateral melt.

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Winton (2000) two-layer column thermodynamics with enthalpy | Shipped | `&ocean_ice_nml enable`, `nk_ice = 2` | Stefan melt verified to under 1 %. |
| Multi-category ice thickness distribution (ITD) restore | Shipped | `ncat > 1` | Whole-category shift, not a Lipscomb (2001) remap. SIS2 does not ship one either. |
| Category ice and snow transport | Shipped | `&ocean_ice_nml transport` | `compress_ice` is not a ridging scheme. |
| C-grid EVP rheology | Shipped | `&ocean_ice_nml dynamics` | Refused with any open, tidal, sponge, clamped or Chapman edge. Runs across a tripolar fold on one rank only (A6, #164; needs `nghost >= 2`). |
| Frazil bank and uptake | Shipped | – | Checks only the top layer (`k = nz`). |
| Ice-to-ocean coupling | Shipped, partial | – | Salt is closed (virtual). Heat is closed, including transmitted shortwave. No freshwater or mass coupling. Ice applies no sea-surface loading. |
| Momentum conserved at fractional cover | Shipped, opt-in | `&ocean_ice_nml a_face_stress` | The default legacy form leaks `(1-a)(tau_a-fxoc)` per face. |
| Snowfall source | Shipped | `&ocean_ice_nml snowfall` | A uniform scalar rate, not a forcing file. Default 0. |
| Snow-ice flooding | Shipped, approximate | `&ocean_ice_nml snow_ice` | The converted ice takes the snow's enthalpy and zero salinity, so it is too cold and too fresh. |
| Ice initial condition | Shipped, analytic only | `&ocean_ice_ic_nml conc_config` | `"zero"` (default), `"uniform"` or `"latitudes"`. No file reader. |
| Ridging | **Not implemented** | – | `compress_ice` thickens ice in its own category, as SIS2 does with `DO_RIDGING=.false.`. |
| Melt ponds | **Not implemented** | – | Summer melt is biased low and the melt-albedo feedback is absent. |
| Lateral melt, floe-size effects | **Not implemented** | – | Marginal-ice-zone retreat is biased. |
| MPI decomposition | Shipped | – | Bit-identical to one rank. Refused on more than one rank under a tripolar north fold. |
| Restart | Shipped | – | Bitwise through the production engine path, on one rank and decomposed. |

### Why the missing pieces matter

- **No ridging.** Under convergence the ice thickens uniformly instead of building a thick ridged tail. The ITD is wrong in a convergent regime, and the EVP ice strength is biased with it. This bites hardest in convergent, near-shore and Antarctic-shear regimes. Porting SIS2's ridging option means vendoring Icepack. It is not a small follow-up.
- **No melt ponds.** The CSIM4 melting-albedo ramp is a crude stand-in.
- **Partial coupling.** The ocean's own background shortwave is not weighted by ice cover, and the ocean shortwave is a single broadband. The ice feels the full wind-stress snapshot, not a concentration-weighted drag law, unless `a_face_stress` is on.
- **Ice cannot sit under a tripolar fold on more than one rank.** On one rank, the category state and `dynamics` (EVP) both fold across the seam.

**Full notes in the archive:** [overview, line counts and the full limits list](CAPABILITIES_ARCHIVE.md#a31)

## 10. Ice-shelf cavities

An ice-shelf cavity is a **static, prescribed ice draft** over the ocean. The model treats the ice as a lid with a pressure load, a basal melt flux and an optional drag. Everything here is default off and byte-identical when off. All of it is **single rank** and runs only with the split solver.

**Status in one line.** The pieces ship and have analytical tests. The ISOMIP+ Ocean0/1/2 files in `validation_examples/ocean/isomip_plus/` are *configuration deliverables, not validated results*. They pass `validate_config` and run. Nothing is tuned to a melt rate. A quiescent sloping-lid cavity is not yet at rest.

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Static ice draft, analytic | Shipped | `&ocean_cavity_dyn_nml enable`, `draft_config = "flat"` or `"linear"` | A grounded column (less than `h_min_cavity` of water, default 10 m) becomes land. |
| Static ice draft from a file | Shipped | `draft_config = "file"`, `draft_sign` | Single rank. The file must already be on the model grid. Rank-3 `(x, y, t)` variable. |
| Barotropic datum with a draft | Shipped | – | `bt_eta` is the deviation from the loaded equilibrium. |
| Ice load in the PGF top boundary condition | Shipped | `&ocean_pgf_nml p_top_in_bc` | **Required** when the draft varies. FV-MOM6 only. |
| Surface pressure in the in-situ EOS | Shipped, opt-in | `&ocean_psurf_nml in_eos` | Ported for FV-Wright and EPBL only. Six other features are refused with it. |
| Basal melt (Holland and Jenkins 1999, three-equation) | Shipped | `&ocean_cavity_melt_nml enable` | Needs `tfreeze_set = "isomip"` and `enable_components`. Meltwater mass is virtual by default. |
| Meltwater mass | Shipped, opt-in | `freshwater = "mass"` | Refused with wet/dry and with `dt_tracer_advect_ratio > 1`. |
| Volume compensation | Shipped, opt-in | `volume_compensation = "uniform_open_ocean"` | Needs `freshwater = "mass"`. Default is `"none"`, the ISOMIP+ position. |
| Ice-shelf top drag | Shipped, opt-in | `&ocean_tdrag_nml enable` | Needs the cavity group. Shares one `C_d` with the melt. No diagnostic catalog entry yet. `cd` defaults to 0, so no drag acts until it is set (`form` defaults to `"quadratic"`). |
| Ice-cover mask on atmospheric forcing | Shipped | follows `&ocean_cavity_dyn_nml enable` | The cover is binary, so the calving front is a one-cell step. |
| KPP and EPBL under ice | Shipped | – | Both use `u_*² = (stress_mag + stress_shelf)/ρ₀`. KPP's depth ignores a stabilising `B_0`. |
| Sponge target independent of the initial condition | Shipped | `&ocean_sponge_nml target_source = "linear_z"` | Needed for ISOMIP+ Ocean1 and Ocean2. |
| ISOMIP+ analytic bedrock | Shipped | `&ocean_topo_nml topo_config = "isomip_plus"`, `x_origin` | – |
| Melt diagnostics | Shipped | `&ocean_diag_nml diags` | Values outside the cover are NaN, never zero. |
| Vertical coordinates | `sigma` and `z_fixed` only | `&vcoord_nml vcoord_type` | `zstar` is refused. Others are refused. `thickness_config` must not be `"uniform_z"`. |
| True isostatic load | **Not implemented** | `draft_source = "in_situ"` | Fails loud. |
| Fractional (area-blended) calving front | **Not implemented** | – | No producer for a fractional cover. |
| Multi-rank cavity | **Refused** | – | Grounding statistics are single-rank reductions. |

### Envelope: what the cavity refuses

Each of these is a fail-loud configure error.

- A PGF other than `fv_mom6`, or `gfs_scale` other than 1.
- A vertical coordinate other than `sigma` or `z_fixed`, or `thickness_config = "uniform_z"`.
- The unsplit driver (`n_inner = 0`), more than one rank, and `bt_halo > 0`.
- Wet/dry, porous barriers, sea ice and `&ocean_tides_nml use_sal`.
- `&ocean_dataovr_nml` file-driven forcing. The reader would restore the unmasked atmosphere under the shelf.
- A varying draft without `p_top_in_bc`.
- With `in_eos = .true.`: kappa-shear, tidal mixing, Redi, isopycnal slopes, `&ocean_pgf_nml reconstruct_for_pressure` and sea ice.

Under `z_fixed` plus a cavity, melt and top drag are accepted, and `kappa_h /= 0` is accepted only with `zfixed_closed_faces`. These are still refused: `use_kpp`, EPBL, ideal age, GM, Redi and slopes, `regrid_time_scale > 0`, kappa-shear, tidal mixing and `reconstruct_for_pressure`. A low `nu_h` (below 2 m²/s) only warns.

### Caveats to read before quoting a number

- **Meltwater is virtual by default.** It adds no volume and no direct buoyancy, only the salinity signature of dilution. For a cavity the meltwater buoyancy *is* the circulation. This is a first-order limit of the default, and configure warns on every virtual run. `freshwater = "mass"` removes it. In a closed box real mass raises sea level, so also set `volume_compensation`.
- **A resting loaded cavity is not exactly at rest.** Under sigma with a sloping draft, the leftover force is the sigma-coordinate PGF truncation. It grows as the cube of the draft slope. A grounding-line cell is decades worse. Measured at `max|u| = 7.08e-8 m/s` over 12 000 s on the small test case.
- **Resolved beds are far worse under sigma.** The ISOMIP+ channel wall gives `rx0 = 0.726`, which is 3.6 times the classical 0.2 bound. `ocean0_idealised_draft.nml` runs three days and then goes negative in thickness. No knob fixes this. Sloping-coordinate PGF corrections are a later phase. `ocean_stability_audit` warns at configure with the computed `rx0`.
- **`z_fixed` is the alternative, and it works with closed faces.** `ocean0_idealised_zfixed.nml` runs the full 30 days with melt and top drag on. The sigma twin blows up on day 15.9. Even so, `z_fixed` is not at rest on a sloping draft. See the archive entries for this section for the numbers, the energy plateau and the open budget question after day 15.
- **KPP does not shoal instantly under melt.** Its depth ignores a stabilising `B_0`. EPBL does respond.
- **Load approximation.** The load is `ρ₀·g·z_draft`, which leaves a depth-uniform `N²·z_draft·∇z_draft` residual in the raw PGF. The default split absorbs it as a small barotropic adjustment.

**Full notes in the archive:**

- [static draft geometry, the barotropic datum and grounding](CAPABILITIES_ARCHIVE.md#a32)
- [the ice load, the pressure seam and the FV-MOM6 top boundary condition](CAPABILITIES_ARCHIVE.md#a33)
- [basal melt, the ice-cover mask, freshwater mass and diagnostics](CAPABILITIES_ARCHIVE.md#a34)
- [ice-shelf top drag](CAPABILITIES_ARCHIVE.md#a35)
- [file-backed draft and the ISOMIP+ analytic bedrock](CAPABILITIES_ARCHIVE.md#a36)
- [z_fixed under a shelf, the version 0.1.0 envelope, the melt-off baseline and the sigma stiffness limit](CAPABILITIES_ARCHIVE.md#a37)
- [cavity features not yet shipped](CAPABILITIES_ARCHIVE.md#a38)

## 11. MPI and domain decomposition

**The contract: a decomposed run is the serial run, bit for bit.** Every owned value of every prognostic field matches the single-rank run on every supported decomposition. It is not a global integral that agrees to round-off. Any mismatch is a bug. The gate is `tests/mpi/test_ocean_decomp_bitid_mpi` (ctest at 1, 2 and 4 ranks).

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Rectangular domain decomposition | Shipped | `&mpi_nml px`, `py` | Left at `1 x 1` on several ranks, the engine picks the split. It is north-south under a tripolar fold. An explicit `px*py` that does not match the rank count is refused. |
| Bit-identity with the serial run | Shipped | – | Checked on sixteen configurations and every `px x py` split of 1, 2 and 4 ranks, under both split schemes. |
| Tripolar north fold | Shipped | `north = "tripolar_fold"` | East-west splits (`px > 1`) must be asked for explicitly. |
| Per-rank file readers | Shipped | – | Supergrid, bathymetry, z-level T/S initial condition and `&ocean_dataovr_nml` forcing are windowed per rank. |
| Console totals | Shipped | `&ocean_diag_nml reproducing_sums` (default on) | The status block prints the same digits on any rank count. |
| Sea ice | Shipped | – | Bit-identical, except under a tripolar north fold (refused). |
| Per-rank restart | Shipped | – | Resumes only on the same decomposition. |
| Restart on a different rank count | **Not implemented** | – | – |
| Parallel (collective) NetCDF | **Not implemented** | – | Output is per rank, merged offline by `tools/merge_output.py`. |
| CUDA-aware halo exchange | Shipped | `-DRDB_CUDA_AWARE_MPI=ON` | A multi-GPU run needs it, and `CUDA_VISIBLE_DEVICES` pinned per rank before `MPI_Init`. The driver refuses a host-staged multi-GPU build. |
| Barotropic march-in (wide halo) | Opt-in, **not bit-identical** | `&ocean_bt_nml bt_halo > 0` | See the list below. The default (`-1`, auto) resolves to 0. |

### Refused at configure on more than one rank

Each check is fail-loud on every rank and names the knob. The check keys on the actual rank count in `engine_setup`, so an auto-factored process grid cannot slip past it.

| Feature | Why |
|---|---|
| `&ocean_cavity_dyn_nml enable` (and the melt / top-drag paths that require it) | the grounding statistics are global reductions the configure does not take; `draft_config="file"` has no windowed reader |
| `&ocean_wetdry_nml enable` | the wet-mask / outflow-limiter halo exchange is not implemented |
| `&ocean_ice_nml enable` with a tripolar north fold | the single-rank fold ships (the ice category state and the EVP velocities are folded across the north seam), but the multi-rank ice fold is not wired: the fold-exchange reserve is not sized for the ice group and there is no decomposition bit-identity case for it (#143) |
| `&ocean_bc_nml` `'chapman'` edges | the edge-uniform eta target is a per-rank partial mean |
| `&ocean_vmix_nml dt_tracer_advect_ratio > 1` | the windowed drain's halo is not wired |
| tripolar fold with a north tile shorter than `nghost + 1` rows (`ny/py`), or tiles narrower than `nghost + 1` columns (`nx/px`, when `px > 1`) | the fold mirrors `nghost` rows below the fold line from the north tile itself; the width rule keeps every tile wider than its own ghost band |
| in-memory geometry injection (the API's staged bathymetry / supergrid arrays) | they describe the whole grid, not a tile (the file readers are windowed) |
| `nghost < 3` | the PPM stencil degrades to first order at a seam face (`ocean_halo_init`) |
| an explicit `&ocean_bt_nml bt_halo` wider than the smallest tile | the wide ghost ring must fit inside the tile |

### Not refused, but not bit-identical across decompositions

These are the only such paths found.

- **The barotropic march-in** (`&ocean_bt_nml bt_halo > 0`) differs from
  the serial run at round-off over variable bathymetry and grossly with
  open boundaries. It is therefore OPT-IN: the default (`bt_halo = -1`,
  AUTO) resolves to 0 on every rank count, and a multi-rank run that sets it
  explicitly logs a warning.
- **An EAST legacy `'sponge'` edge together with a `'clamped'` north edge**,
  split in x, drifts at round-off on the east boundary face.  A west sponge
  edge (alone, or with Orlanski-open and clamped edges) and every other
  combination tested are exact.
- **`eulerian_z` under `ssp_rk2` with EPBL + Fox-Kemper MLE** drifts at
  round-off (3e-11 relative after 24 steps); EPBL alone, the pair on any
  other coordinate, or under `pred_corr` is exact (row
  `decomp_eulerian_z_ssp_rk2`).

### Global 1-degree scaling

`validation_examples/ocean/global_1deg/` (360 x 320 x 50, dt = 1800 s) runs on 1, 2 and 4 ranks with north-south splits. The final state is bit-identical on all three. Wall time per simulated day, stepping only:

| Ranks | V100s, one per rank (nvfortran 26.5, HPC-X, CUDA-aware) | CPU cores, one per rank (gfortran 15.1, OpenMPI 5.0.5, Xeon E5-2698 v4) |
|---|---|---|
| 1 | 12.9 s / simulated day | 910 s / simulated day |
| 2 | 8.1 s / day (1.6x) | 492 s / day (1.85x) |
| 4 | 5.8 s / day (2.2x) | 252 s / day (3.6x) |

At 1 degree a V100 is under-filled (5.8 M cells, 1.4 M per GPU on 4); the
GPU scaling is halo-latency bound, the CPU scaling near-linear.

**Full notes in the archive:** [the bit-identity gate, the configuration list and the toolchains](CAPABILITIES_ARCHIVE.md#a39)

## 12. GPU, compilers and performance

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| One source for CPU and GPU (`do concurrent` plus OpenACC) | Shipped | `RDB_ENABLE_GPU` | Runs on NVHPC GPU and multicore, gfortran and ifx. An OpenMP-target build is also exercised. |
| NVHPC (GPU and multicore stdpar) | Shipped, primary target | – | Exercised in CI. |
| gfortran | Shipped | – | Exercised in CI for portability. |
| ifx | Buildable | – | Some OpenMP-target map-clause edge cases were seen in vendor tooling. They are not load-bearing for the canonical NVHPC build. |
| AMD and Intel GPUs | **Not implemented** | – | There is no HIP or SYCL backend. It would be a portability research project. |
| CUDA-aware MPI halo exchange | Shipped | `RDB_CUDA_AWARE_MPI` | See [section 11](#11-mpi-and-domain-decomposition). |
| Up-front memory report per slot | Shipped | `rdb_mem_report` | Reconciled against the measured device mapping, so a gated-off closure provably costs nothing. |
| Gated allocation for default-off closures | Shipped | – | EPBL, kappa-shear, tidal mixing, GM, Redi, MLE, VarMix, MEKE and isopycnal slopes. A plain run pays none of their multi-GB footprint. |
| Mixed-precision tracers (fp32 or bf16) | Not implemented | – | – |
| Dynamic load balancing for shifting wet fractions | Not implemented | – | – |
| Topology-aware (NVLink or NVSwitch) MPI tuning | Not implemented | – | – |
| Fault tolerance (ULFM, mid-run checkpoint-restart on rank failure) | Not implemented | – | – |
| GPU CI runner | Not implemented | – | GPU tests run locally today. |

### Headline numbers (single V100)

| Configuration | Wall time |
|---|---|
| Tasman 2 km (eddy-resolving, KPP, ALE, sponges) | about 34 s per simulated day |
| MOM6 double-gyre reference (44 x 40 x 2, dt = 1200 s) | about 16 s per 30 simulated days |

Multi-rank scaling of the global 1-degree case is in [section 11](#11-mpi-and-domain-decomposition).

## 13. Input, output, restart and initial conditions

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Per-rank NetCDF output (CF-1.8) | Shipped | `&ocean_diag_nml` | One file per rank. Merge offline with `tools/merge_output.py`. |
| Diagnostics registry with cadence dispatch | Shipped | `&ocean_diag_nml diags` | Per-diagnostic output vertical remap onto z, sigma, z* or density bins. Time operations: INSTANT, MEAN, MAX, MIN and INTEGRAL. |
| Deflate compression | Shipped, opt-in | – | 3 to 5 times smaller. |
| Single-precision diagnostic output | Shipped, opt-in | `&ocean_diag_nml output_precision = "single"` | Diagnostic stream only. Restarts and the console totals are always working precision. Regenerate any baseline that byte-compares diagnostic NetCDF. |
| Restart and warm start | Shipped, bit-exact | `cfg%restart_interval`, `cfg%restart_file` | Per rank. Resumes only on the same decomposition. Diagnostic MEAN, MAX and MIN windows reset on restart by design. |
| Analytic initial conditions | Shipped | `&ocean_ic_nml` | Uniform, linear `T(z)`, Eady front, geostrophic adjustment and per-layer `gprime` density. |
| Geopotential `T(z)`, `S(z)` overlay | Shipped, opt-in | `&ocean_zinit_nml enable`, `source = "file"` or `"linear"` | The file must be on the model horizontal grid. Needs `RDB_ENABLE_NETCDF=ON`, even for `source = "linear"`, which opens no file. |
| On-the-fly horizontal regrid with land flood-fill | Not implemented | – | Deferred to a Python preprocessor. |
| Conservative vertical remap of the initial condition | Not implemented | – | Version 1 uses point linear interpolation at layer centres. |
| MPI I/O server hand-off for the diagnostics manager | Not implemented | `&output_nml use_io_server` | Accepted, but it warns and is ignored. The emit is serial per-rank NetCDF. |
| Parallel collective NetCDF | Not implemented | – | – |
| Gauge or station time series | Not implemented | – | – |
| Lagrangian drifters | Not implemented | – | – |
| In-situ visualisation | Not implemented | – | – |

### Notes

- **Restarts are bit-exact.** A resume reproduces the writer's next steps bitwise, ghosts included. This holds on the production engine path as well as for hand-built states. It includes sea ice, wet/dry hysteresis and the `pred_corr` predictor's reused viscous tendency. A 4-rank 1/4-degree Southern Ocean resumed from day 730 matches the straight run bit for bit at day 732.
- **A missing required field on read is fatal.** Some fields are optional, so an older checkpoint resumes with a documented degradation.
- **Diagnostics that cannot run.** Requesting a canonical diagnostic whose feature gate is closed (for example `age` with `enable_ideal_age = .false.`) logs a warning that names the diagnostic and the likely gate. `name:off` stays a silent no-op.
- **`land_fill_t` and `land_fill_s` are inert on interior land cells.** The land-state contract (`h_layer = H_VANISHED`, `hTr = 0`) overrides them.

**Full notes in the archive:**

- [restart round-trips, engine path, sea ice and the flux-component state](CAPABILITIES_ARCHIVE.md#a40)
- [the diagnostics manager](CAPABILITIES_ARCHIVE.md#a41)
- [restart and checkpoint registry](CAPABILITIES_ARCHIVE.md#a42)
- [initial conditions and the geopotential T/S overlay](CAPABILITIES_ARCHIVE.md#a43)

## 14. Validation and conservation

### What is validated

| Area | Status | Where | Limits |
|---|---|---|---|
| Analytical and exact-solution tests | Shipped | `tests/` | Vertical diffusion against erfc, gravity-wave phase speed, geostrophic adjustment, lake-at-rest over a seamount, island-at-rest, Thacker moving shoreline (wet/dry), Eady growth rate, Nansen free drift and Stefan melt (sea ice). |
| Test suite | Shipped | `ctest -L ocean`, `ctest -L core` | Count the entries with `ctest -N -R rdb`, and the `ocean` and `core` subsets with `ctest -N -L ocean` and `ctest -N -L core`. Regime labels are mandatory and enforced by a pre-commit hook. |
| Multi-rank consistency | Shipped | `tests/mpi/` | C-grid halo tests and the bit-identity gate (see [section 11](#11-mpi-and-domain-decomposition)). |
| CPU and GPU bit comparison | Shipped | gcc plus nvhpc CI matrix | – |
| Golden-output regression harness | Shipped | `tests/regression/` | – |
| Canonical benchmark configurations | Shipped | `validation_examples/ocean/` | Each has a README with expected behaviour and failure-mode notes. |
| Observational comparison for a real domain | **Not implemented** | – | – |
| Performance regression CI | **Not implemented** | – | – |

`validation_examples/test_cases/tasman_validation/` carries the regional real-bathymetry case.

### Conservation contract

The periodic `[stats]` console block prints a `Mass`, `Salt` and `Heat` line, each with an `Error`. For a conserved quantity Q:

```
Error = (Q_total - Q_ref) + out_Q - src_Q
```

`Q_ref` is latched at the first status report. `out_Q` is the cumulative net outflux through open boundaries. `src_Q` is the cumulative surface source. On a closed, unforced, conservative run all three terms are zero, and `Error` is the numerical leak. With open boundaries or surface forcing, the budget terms keep `Error` at round-off instead of flagging the intended physics.

- **The Error floor.** With `&ocean_diag_nml reproducing_sums = .true.`, which is the default since v0.1.0, the totals are summed in Extended Fixed Point. They are bit-identical on every rank count. The residual floor drops to about 1e-54 relative. The knob changes diagnostic text only, never the trajectory. With `.false.`, the floor is summation noise of about 1e-14 to 1e-15 relative.
- **Cost.** About +0.5 % on the global 1-degree grid on one V100.
- **Instrumented paths.** Mass, salt and heat all close to a residual for the split and unsplit ocean paths. This includes windowed tracer advection (since 2026-09-12) and Redi with an open edge (since 2026-10-02).
- **Known limitations.** None at present. The `horiz_adv_budget_valid` gate in `ocean_budget_is_active` stays for the next un-instrumented transport path.

**Full notes in the archive:**

- [the full conservation contract, the summation-order floor and instrumented paths](CAPABILITIES_ARCHIVE.md#a44)
- [validation tests and canonical benchmarks in the original wording](CAPABILITIES_ARCHIVE.md#a45)

## 15. Open problems and the working envelope

### The envelope proven production-stable

This namelist set is the regime that reaches day 580 on the MOM6-reference double gyre. Outside it (lower `nu_h`, another drag form, another PGF on real bathymetry), each case needs its own validation.

```
&ocean_hvisc_nml      nu_h = 10000.0,  smag_ah = .true. /
&ocean_bdrag_nml      form = "linear", hbbl = 10.0, bg_vel = 0.1, r = 2.5e-5 /
&ocean_coriolis_nml   form = "sadourny" /
&ocean_pgf_nml        form = "gprime", maxvel = 6.0 /   ! 2-layer reduced gravity + clamp
                                                        ! (optional cfl_trunc = 0.5 adds an
                                                        !  advective-CFL velocity truncation
                                                        !  ahead of the absolute maxvel cap)
&ocean_bt_nml         auto_n_inner = .true. /
&vcoord_nml           vcoord_type = "sigma" /
&time_nml             dt_fixed = 1200.0 /
```

Ocean knobs live in per-concern sub-namelists named `&ocean_<group>_nml`, with the `ocean_` prefix dropped from each key. Run `tools/nml_split.py <file>` to migrate an old `&ocean_setup_nml` namelist.

### Open problems

| Problem | Where it bites | Workaround | Status |
|---|---|---|---|
| Nonlinear collapse of the Eady front | `validation_examples/ocean/eady/` past day about 65 | Keep to the shipped 60-day linear window. `dt = 300 s` survives to day 100 but with `En` above the front's available potential energy. | Open. |
| Resting-state growth | Quiescent stratified runs at low Laplacian viscosity | Use `pred_corr`, which slows it about 36 times. A constant `nu_h` (for example 100) hides it without changing the rate. | Attributed to the `ssp_rk2` outer split. `pred_corr` does not remove it. |
| Spurious 2Δy mode over a slope | Quiescent terrain-following runs on an f-plane, `ny >= 2` | Set `&vcoord_nml remap_boundary_extrap = .true.`. Failing that, use a constant `nu_h` of about 50 m²/s or `nu_4` of about 1e8 m⁴/s at 2 km. | Flow-aware closures (Smagorinsky, Leith) do not help, because they give no viscosity at rest. |
| Sloping-lid cavity not at rest | Sigma ice-shelf cavities | `z_fixed` with closed faces | Open. Needs the sloping-coordinate PGF corrections. See [section 10](#10-ice-shelf-cavities). |
| Open-staircase PGF defect | `z_fixed` without closed faces, `hycom` on stepped topography | For `z_fixed`: `zfixed_closed_faces = .true.`. None yet for `hycom`. | `z_fixed` is refused without the knob over a stepped bed. `hycom` is a scoped XFAIL. |
| KPP bottom boundary layer | KPP runs | – | Not implemented. The `V_t²` term and non-local transport ship. |
| Hollingsworth-Källén Coriolis correction | Production climate-scale runs | The PV form is production | Not implemented. It is a guard option only. |
| Non-hydrostatic pressure on the C-grid | – | – | Not implemented. |

### Why the slope mode matters

A sigma column at rest over a sloping bed, on an f-plane, with no constant lateral viscosity, grows a grid-scale 2Δy mode in `u`. It grows exponentially and does not saturate. The pressure-gradient truncation error is only the seed. The amplifier is the first-order boundary-cell closure of the remap: the top and bottom layers reconstruct as PCM, which injects a spurious diapycnal tracer flux there. The remap option `remap_boundary_extrap` removes the source. A constant viscosity damps it. The viscous requirement grows as `dy` shrinks, and it has not been measured as a function of resolution. A forced, energetic run sits decades above this floor and never notices it.

**Full notes in the archive:**

- [quiescent terrain-following runs over a slope (measurements and isolation)](CAPABILITIES_ARCHIVE.md#a46)
- [Eady-front collapse and resting-state growth](CAPABILITIES_ARCHIVE.md#a47)
- [Hollingsworth-Källén and non-hydrostatic items in the original wording](CAPABILITIES_ARCHIVE.md#a48)

## 16. Data assimilation and machine learning

| Feature | Status | Key knob | Limits |
|---|---|---|---|
| Relaxation toward a reference state | Shipped | `&ocean_sponge_nml`, per-edge asymmetric nudging | This is all the model offers for assimilation. It uses the namelist and restart surface. |
| In-model observation operators | Not implemented | – | – |
| Ensemble support | Not implemented | – | – |
| 4D-Var adjoint | Not implemented | – | – |
| In-process state-setter API | Not implemented | – | It left with the FFI carve-out. |
| Machine-learning surrogate hook | Not implemented | – | There is no in-process API surface for a surrogate to plug into, since the C/Python FFI was carved out with the coastal repository. That seam would have to be rebuilt first. |

### How this page relates to the roadmaps

The roadmap (`ROADMAP_OCEAN.md`) carries the full work checklist. It has project-scope tags (`[E]`, `[M]`, `[H]`, `[V]`) and a `Validate:` line for each planned item. This page summarises the *current* state. If the two disagree, the codebase is the authority and the roadmaps are the next-best source.

## 17. Revision notes

This page was restructured on 2026-10-12 so that people can read it. The [archive](CAPABILITIES_ARCHIVE.md) keeps the earlier wording. Two kinds of change were made to the facts. Statements that disagreed with the source were corrected, and pure duplicates were merged. Both are listed below for maintainers.

A review pass on the same day added rows C12 to C17 and the list of other changes. Rows C1 to C11 were already rewritten in place in the long-form text (now in the archive), so their earlier wording survives only in the "Earlier statement" column. Where earlier wording did survive in the archive (C5, C15, C17), it carries an inline *Superseded* marker.

**Full notes in the archive:**

- [Corrections made against the source on 2026-10-12](CAPABILITIES_ARCHIVE.md#a49)
- [Statements that look stale but were not changed](CAPABILITIES_ARCHIVE.md#a50)
- [Other changes made in the review pass (2026-10-12)](CAPABILITIES_ARCHIVE.md#a51)
- [Duplicates that were merged](CAPABILITIES_ARCHIVE.md#a52)
- [Where the old sections went](CAPABILITIES_ARCHIVE.md#a53)
