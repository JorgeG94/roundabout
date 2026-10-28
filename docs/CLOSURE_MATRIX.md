# Roundabout — Closure & Scheme Matrix

This page answers one question: **which closure or scheme is switched on in the
ocean model, how do I turn it on or off, and which test guards it?**

Use the tables as a quick lookup. When you need the full story (formulas,
measured numbers, MOM6 parameter mappings, history and `file:line` pointers),
follow the links at the end of each section to the
[archive](CLOSURE_MATRIX_ARCHIVE.md). Nothing was dropped when this page was
restructured; the long text moved into the archive.

Last reconciled with the code: **2026-10-13** (`roquet_rho` added; `roquet_spv` and `roquet_rho` answers changed by the PT→CT fix, see §6). The dated change log is in
[Appendix A](#appendix-a-change-history).

## How to read this page

- **One section per concern** (vertical mixing, lateral closures, pressure
  gradient, and so on). Each section starts with a short table, then links to the
  original long-form text in the archive.
- **Knob notation.** `&ocean_epbl_nml enable` means the namelist group
  `&ocean_epbl_nml` and the key `enable`. Ocean knobs live in
  `&ocean_<group>_nml`; the `ocean_` prefix is dropped from each *key* inside
  the group (for example `&ocean_bt_nml split_scheme`). Shared vertical-coordinate
  config lives in `&vcoord_nml`.
- **Regime.** This build ships the ocean regime only (`sim_type='ocean'`,
  Arakawa C-grid). The coastal A-grid and unstructured paths were split into a
  separate repository.
- **Defaults** in the tables were checked against `src/core/rdb_config.F90`.
  When this page and the code disagree, **the code is the authority and this
  page is the bug**.
- **Defaults are a convenience copy.** The authoritative, generated
  knob-to-default list is [`generated_nml_knobs.md`](generated_nml_knobs.md).
  The source of truth behind it is `src/core/rdb_config.F90`.

### Status words

| Word | Meaning |
|---|---|
| **default** | On unless you choose something else. |
| **opt-in** | Off until you set the knob. With the knob off, results are bit-identical to a build without the feature. |
| **available** | One of several choices for a selector (for example another equation of state). Not the default. |
| **refused** | A configure-time error. `validate_config` fails loud and `./rdb --validate-only case.nml` exits with code 3. |
| **experimental** | Supported and fully tested, but its answer has a known cost. The label describes the answer, not the code path. Only `ssp_rk2` carries it. |
| **retired** | The key still parses, but setting it fails loud and the message says why. |
| **reserved** | The enum value is fixed but not implemented. It is refused, never silently replaced by a default. |
| **validation-grade** | Works, but is not meant for production runs. |
| **not available (`—`)** | A documented absence. Each one has a row so it is not mistaken for an omission. |

### Jargon, defined once

| Term | Plain meaning |
|---|---|
| **bit-identical / byte-identical** | The output is exactly the same as before the feature existed. Most new knobs default off for this reason. |
| **fail-loud** | The model stops with a clear error instead of quietly doing something else. |
| **configure time** | Setup, before the first step. `./rdb --validate-only case.nml` runs only these checks. |
| **layer index** | Layers are numbered bottom-up. `k=1` is the bed and `k=nz` is the surface. |
| **`k_bot` / `k_top`** | `multilayer_state_t%k_bot` is the first layer that is really water, counting up from the bed. `k_top` is the first one counting down from the surface or ice base. They equal `1` and `nz` when no layer is a filler. |
| **filler / vanished layer** | A layer thinner than `H_VANISHED` (1.5e-4 m). It holds no water and takes no part in the physics. |
| **closed face (partial step)** | A velocity face that is a wall for one layer, because that layer is a filler on one side. Switched by `&vcoord_nml zfixed_closed_faces`. |
| **BT / barotropic** | The fast, depth-averaged part of the split-explicit time step. The slow part is the baroclinic (layered) step. |
| **thermo cadence** | Closures that run once per thermodynamic step instead of every dynamics stage. |
| **ALE remap** | The vertical re-gridding step. Layers move to their target thicknesses and tracers are redistributed conservatively. |
| **PGF** | Pressure-gradient force. |
| **PV** | Potential vorticity. |
| **En** | The energy diagnostic (m²/s²) printed on the console. |
| **XFAIL (scoped)** | A test known to fail for a documented physical reason. It is pinned so the failure cannot grow unnoticed. |
| **MOM6** | The reference ocean model Roundabout follows. "MOM6 parity" means matching its behaviour or defaults. |

More acronyms: [`codebase/CONCEPTS.md`](codebase/CONCEPTS.md). Prose and
limitations: [`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).
Physics and namelist reference: [`REFERENCE.md`](REFERENCE.md).

**Full notes in the archive:** [the original preface of this page](CLOSURE_MATRIX_ARCHIVE.md#a1)

## Contents

1. [Time stepping and the barotropic mode](#1-time-stepping-and-the-barotropic-mode)
2. [Vertical coordinates and ALE remap](#2-vertical-coordinates-and-ale-remap)
3. [Pressure gradient force](#3-pressure-gradient-force)
4. [Coriolis](#4-coriolis)
5. [Tracer advection and reconstruction](#5-tracer-advection-and-reconstruction)
6. [Equation of state](#6-equation-of-state)
7. [Vertical mixing](#7-vertical-mixing)
8. [Lateral closures and eddy parameterizations](#8-lateral-closures-and-eddy-parameterizations)
9. [Drag and surface stress](#9-drag-and-surface-stress)
10. [Subgrid topography: porous barriers and partial-step faces](#10-subgrid-topography-porous-barriers-and-partial-step-faces)
11. [Surface forcing, initial conditions and other fluxes](#11-surface-forcing-initial-conditions-and-other-fluxes)
12. [Sea ice](#12-sea-ice)
13. [Ice-shelf cavity](#13-ice-shelf-cavity)
14. [Boundaries](#14-boundaries)
15. [Diagnostics](#15-diagnostics)
16. [Tunable parameters](#16-tunable-parameters)
17. [Keeping this current](#17-keeping-this-current)
- [Appendix A: change history](#appendix-a-change-history)

---

## 1. Time stepping and the barotropic mode

The model uses a split-explicit step. A slow, layered (baroclinic) outer step
wraps a fast, depth-averaged (barotropic) loop. `&ocean_bt_nml split_scheme`
picks the outer step. Both schemes ship and both are tested; the stability suite
runs an `ssp_rk2` twin of every case whose namelist does not pin a scheme.

| Scheme or feature | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| `pred_corr` | **default** | `&ocean_bt_nml split_scheme = "pred_corr"` | MOM6 predictor-corrector. Neutrally stable to `ω·dt = 2`. Refused outside its v1 envelope (`eulerian_z`, wet/dry, `dt_tracer_advect_ratio > 1`); those cases pin `ssp_rk2`. | `test_ocean_pred_corr`, `test_ocean_cor_ref_seiche` |
| `ssp_rk2` | **experimental** (fully supported) | `split_scheme = "ssp_rk2"` | Two identical stages plus an SSP average. Widest envelope. Amplifies internal gravity waves by `√(1+(ω·dt)⁴/4)` per step. | `test_ocean_dyn_split`, `test_ocean_analytical` |
| Fast-loop wave damping (MOM6 `BEBT`) | **default** | `&ocean_bt_nml bebt` (`0.1`; `0.0` is neutral) | Damps barotropic gravity waves in the fast loop. | `test_ocean_barotropic_substep` |
| OM3-style barotropic step (MOM6 `DTBT`, `DTBT_RESET_PERIOD`) | opt-in; replaces `auto_n_inner`/`cfl_bt_safety` (refused together) | `&ocean_bt_nml dtbt` (unset), `dtbt_reset_period` (`-1`) | `1/dt_max² = ½(1+2·bebt)[(1/A_T)Σ g·D·dy/dx + Σ f²_corner]` over wet cells, global minimum; `dtbt < 0` is that fraction of `dt_max` (`0` means `-0.98`), `> 0` fixed seconds; `n_inner = ceiling(dt/step)`. Reset `0` re-derives `n_inner` every step from the live `Σh_layer`, `> 0` every that many seconds, `< 0` at setup only. | `test_ocean_bt_cfl` |
| BT forcing by depth-mean baroclinic PGF (MOM6 `BT_force`) | **default** | `bc_pgf_forcing` (`.true.`) | Keeps the baroclinic bottom-pressure gradient in the barotropic forcing. `.false.` is the legacy split. | `test_ocean_bt_baroclinic_forcing` |
| Per-layer bc-PGF correction (MOM6 `btstep_layer_accel`) | opt-in; **refused** with `zfixed_closed_faces` | `correction_bc_pgf` (`.false.`) | Adds a depth-mean-zero pressure response to each layer. Needs `&ocean_pgf_nml form = "fv_mom6"`. | `test_ocean_pgf_p_top_bc`, `test_ocean_zfixed_closed_faces` |
| Fast-loop Coriolis reference (MOM6 `ubt_Cor`) | **default** (follows the scheme) | none | Removes the Coriolis term already frozen into the slow forcing so it is not integrated twice. | `test_ocean_bt_cor_ref`, `test_ocean_cor_ref_seiche` |
| BT-correction fold weight (MOM6 `u_accel_bt`) | **default**: uniform, the only form | `correction_visc_rem`, `correction_h_weighted` are **retired** | Every layer gets the same barotropic increment. | `test_ocean_bt_correction_weight`, `test_ocean_visc_rem` |
| Viscous-remnant chain (MOM6 `av_rem` / `bt_rem`) | **default** on | `&ocean_bt_nml visc_rem_chain` (`.true.`) | One switch for the `visc_rem` producer plus `forcing_visc_rem`, `renorm_visc_rem` and `bt_rem_from_visc_rem`. Inert unless `bbl_glue` or `implicit_drag` is on. | `test_ocean_bt_rem_from_visc_rem`, `test_ocean_visc_rem` |
| Face-thickness mean for BT depth means | **default** `"hybrid"` | `&ocean_bt_nml frhat_scheme` (`"hybrid"` or `"arithmetic"`) | MOM6 `HVEL_SCHEME=HYBRID`. Active only with `zfixed_closed_faces`; otherwise it falls back to arithmetic. | `test_ocean_bt_frhat` |

**Why `ssp_rk2` is "experimental".** On `resting_stratified_channel.nml` (flat
bed, at rest, no energy source) `ssp_rk2` manufactures `En = 2.992E-05 m²/s²` by
day 25, where `pred_corr` holds `1.739E-09`. The label describes that cost. It
does not deprecate the scheme. Do not close the scoped XFAIL by raising
`en_rest_max`.

**Full notes in the archive:** [full split-scheme table and notes (original text)](CLOSURE_MATRIX_ARCHIVE.md#a2)

## 2. Vertical coordinates and ALE remap

`&vcoord_nml vcoord_type` picks how layers are placed. Ten families ship. The
ALE remap then moves the layers to their target thicknesses. The remap method is
chosen separately by `&vcoord_nml remap_method`.

| Coordinate | Status | `vcoord_type` value | What it does | Tests |
|---|---|---|---|---|
| Sigma | **default** | `"sigma"` | Terrain-following: layers share `(H+η)` in fixed fractions. | `test_ocean_vcoord_interface_depths` |
| Z-sigma hybrid | **refused** | `"zsigma"` | The deep branch reads a dimensionless table as metres, so the whole column collapses into the bed layer. | `test_ocean_vcoord_hygiene`, `test_ocean_vcoord_interface_depths` |
| Z-star (MOM6 z*) | available | `"zstar"` (alias `zstar_lite`) | The `z_fixed` profile stretched per column by `(H+η)/H`. Accepts `zfixed_closed_faces`. Refused under a cavity and with wet/dry. | `test_ocean_vcoord_interface_depths`, `test_ocean_zstar_closed_faces` |
| Z-star / sigma hybrid | available | `"zstar_sigma"` | Sigma in shallow water, z*-style in deep water. | — |
| Z-star full (per column) | available; **slated for retirement** in favour of `"zstar"` | `"zstar_full"` | Fine z-zone at the surface over a terrain-following coarse zone. | `test_ocean_zstar_full_closed_faces` |
| Eulerian-Z | available | `"eulerian_z"` | `H·dsig(k)`; `η` is ignored. Refused with the default `pred_corr` scheme, so pin `ssp_rk2`. | — |
| Z-fixed (gprime) | available | `"z_fixed"` (aliases `z_levels`, `gprime`) | Fixed-depth interfaces; bed-side layers vanish to fillers. The only z-like family legal under an ice-shelf cavity. | `test_ocean_vcoord_hygiene`, `test_ocean_zfixed_k_bot` |
| Lagrangian / isopycnal | available | `"lagrangian"` (alias `isopycnal`) | Target equals the live layers, so the remap is a no-op. | — |
| Isopycnal (rho) | **validation-grade** | `"rho"` | Interfaces sit on prescribed density surfaces. Weakly stratified columns collapse. | — |
| Hybrid z*/isopycnal (HYCOM) | available | `"hycom"` | Isopycnal interior with a z* near-surface floor in metres. The production isopycnal coordinate. | `test_ocean_vcoord_hycom` |

The `lagrangian` row is not in the old table; the alias list is from
`src/ALE/rdb_vcoord.F90:730`.

### ALE remap methods

All methods conserve column integrals to round-off and are monotone.

| Method | `remap_method` | Order | What it does |
|---|---|---|---|
| Piecewise constant | `"pcm"` | 0 | Donor cell. Diffusive. |
| Piecewise linear | `"plm"` | 1 | Minmod-limited. |
| Piecewise parabolic | `"ppm"` (**default**) | 2 | Colella & Woodward 1984. |
| PPM, non-uniform H4 edges | `"ppm_h4"` | 2 (4th-order edges) | Thickness-weighted edge estimate; less spurious mixing. |
| Piecewise quartic | `"pqm"` | 4 to 5 | White & Adcroft 2008. Falls back to PPM when `N<5`. |

| Remap add-on | Status | Knob (default) | What it does | Test |
|---|---|---|---|---|
| Boundary-cell closure | opt-in | `&vcoord_nml remap_boundary_extrap` (`.false.`) | Makes the two outermost cells linear-exact instead of dropping to PCM. | `test_remap_boundary_extrap` |
| Non-uniform-grid weights | opt-in | `&vcoord_nml remap_nonuniform_weights` (`.false.`) | Uses true layer thicknesses in the PLM and PPM stencils. Swaps PLM's limiter too. | `test_remap_nonuniform` |

**Full notes in the archive:**

- [full vertical-coordinate table (original cells, MOM6 mappings, cavity envelope)](CLOSURE_MATRIX_ARCHIVE.md#a3)
- [ALE remap reconstruction, boundary closure and non-uniform weights (formulas and measurements)](CLOSURE_MATRIX_ARCHIVE.md#a4)

## 3. Pressure gradient force

`&ocean_pgf_nml form` picks the pressure-gradient discretisation.

| Variant | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Montgomery potential | **default** | `form = "mont"` | Vertical recursion plus one horizontal difference. Valid on any bathymetry and vcoord; exact at rest in isopycnal columns. | `test_ocean_pgf_mont`, `test_ocean_pgf_fv` |
| FV-lite, FV-Wright | available | `form = "fv_lite"`, `"fv_wright"` | Finite-volume forms. Surface-relative. | `test_ocean_pgf_fv` |
| gprime | available | `form = "gprime"` | Two-layer reduced gravity; the interface term differences the interface HEIGHT `h_1 - b`, so a flat interface over a sloping bed is at rest. | `test_ocean_pgf_gprime` |
| FV-MOM6 | available | `form = "fv_mom6"` | MOM6 analytic finite-volume PGF with a pressure-stack surface BC. | `test_ocean_pgf_sigma_rest` |
| In-layer T/S reconstruction | opt-in; **refused** under `z_fixed` with a cavity | `reconstruct_for_pressure` (`.false.`), `recon_scheme` (`1`) | PLM or PPM T/S inside each layer. Exact at rest for linear EOS and linear T/S. FV-MOM6 only. | `test_ocean_pgf_reconstruct`, `test_ocean_pgf_sigma_rest`, `test_ocean_pgf_eos_fast`, `test_ocean_vcoord_hygiene` |
| In-situ-pressure density (MOM6 `RECONSTRUCT_FOR_PRESSURE=False`) | **default** (FV-MOM6, pressure-dependent EOS only) | `insitu_density` (`.true.`) | Evaluates density at in-situ pressure. `.false.` restores the potential-density integral. Refused with `&ocean_psurf_nml in_eos`. | `test_ocean_pgf_wright_analytic`, `test_ocean_pgf_eos_fast`, `test_ocean_pgf_insitu` |
| Top-of-column load in the BC | opt-in | `p_top_in_bc` (`.false.`) | Adds `ms%p_top` to the pressure-stack surface BC. Required for a cavity whose draft varies. | `test_ocean_pgf_p_top_bc` |
| Reference densities | **default** | `&ocean_ic_nml rho_0` (`1035`) | One configured `ρ₀` feeds the Boussinesq divisor and the FV-MOM6 baseline. | `test_ocean_pgf_rho_ref` |
| Grounded-layer gate | **default** (Lagrangian vcoord only) | `&ocean_isopycnal_nml pgf_skip_nonoverlap` (`.true.`) | Zeroes the face PGF where a grounded isopycnal layer does not overlap its neighbour. | `test_ocean_pgf_grounded`, `test_ocean_lagrangian_staircase` |

**Full notes in the archive:** [full PGF table (method descriptions, measurements, costs, tests)](CLOSURE_MATRIX_ARCHIVE.md#a5)

## 4. Coriolis

`&ocean_coriolis_nml form` picks the Coriolis-advection scheme.

| Scheme or add-on | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Sadourny enstrophy PV-flux | **default** | `form = "sadourny"` | Velocity form. | — |
| Sadourny energy-conserving | available | `form = "sadourny_energy"` | Transport form `q·vh` (MOM6 `SADOURNY75_ENERGY`). | — |
| Sadourny plus Hollingsworth-Källén guard | available | `form = "sadourny_hk"` | Arakawa-Hsu PV stencil. A pair-floored twin runs on partial-step z-level grids. | `test_ocean_coriolis_hk_vanished` |
| PV face interpolation | **default** `"centered"` | `pv_adv_scheme` (`"centered"`, or `"weno3"`, `"weno5"`, `"weno7"`) | Corner vorticity to face interpolation in the Sadourny path. WENO is upwind-biased (MOM6 `WENOVI*TH`). `weno5` and `weno7` need `nghost >= 3` and `4`. | `test_coriolis_multilayer` |
| Mass-consistent CorAdCalc | opt-in | `use_state_fluxes` (`.false.`) | The corrector reuses the predictor continuity's transports. Needs `sadourny_energy` and `pred_corr`. | `test_coriolis_multilayer` |
| `BOUND_CORIOLIS` clamp | opt-in | `bound_coriolis` (`.false.`) | Clamps the Coriolis term into the `(f+ζ)·v` range before the KE gradient. `sadourny_energy`, or `sadourny` (enstrophy) with `pv_adv_scheme = "centered"`; `sadourny_hk` and WENO PV fail loud. Inert on all-wet energy-form columns. | `test_coriolis_multilayer` (`bound_coriolis_enstrophy_*`) |
| PV corner thickness | **default** `"cell_mean"` | `corner_h` (`"mom6_area"` is available) | MOM6's area-weighted form. Algebraically identical above the `H_MIN_PV` floor. `sadourny_energy` only. | — |

The `pv_adv_scheme` and `use_state_fluxes` rows are not in the old table; they
come from `src/core/rdb_config.F90:207` and `:213`. The test links were checked
by `grep` in `tests/test_coriolis_multilayer.F90`.

**Full notes in the archive:** [full Coriolis table (HK pair floor, thin-layer behaviour, measurements)](CLOSURE_MATRIX_ARCHIVE.md#a6)

## 5. Tracer advection and reconstruction

Tracer transport has three separate reconstruction jobs: horizontal advection,
vertical (in-z) advection, and the conservative ALE vertical remap. They are
separate code paths even when they share a paper.

| Job or modifier | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Horizontal tracer advection | **default** `"ppm"` | `&ocean_vmix_nml tracer_advection_scheme` | Van Leer-slope PPM edges, upwind pick, no swept-volume correction; x then y every step; the two cells next to a wall drop to first-order donor. | `test_ocean_tracer_ppm_h3` |
| Horizontal tracer advection, H3 | opt-in; also drives the windowed drain (`dt_tracer_advect_ratio > 1`, requires `tracer_recon = "ppm"`) | `tracer_advection_scheme = "ppm_h3"` | Huynh (1997) 3rd-order bounded edges, Colella & Woodward (1984) limiter, swept-volume flux (OM3 `PPM:H3`). Width-3 stencil. In the drain it replaces the CW parabola / WENO face value and is fed per-pass volumes. | `test_ocean_tracer_ppm_h3`, `test_ocean_tracer_ppm_h3_drain` |
| Vertical (in-z) advection | **default**, no knob | hardcoded | First-order upwind. | — |
| ALE vertical remap | **default** `"ppm"` | `&vcoord_nml remap_method` | See [section 2](#2-vertical-coordinates-and-ale-remap). | — |
| Positivity limiter | opt-in | `&ocean_continuity_nml ppm_limit_pos` (`.false.`) | Modifier on the continuity-PPM face values. | — |
| Positive-definite continuity | opt-in; **refused** with wet/dry | `positive_definite` (`.false.`) | `2·h_lim` PPM edge floor plus a per-donor outflux limiter. Creates no mass. | — |
| Consistent barotropic renormalisation | **default** | `renorm_consistent_flux` (`.true.`) | MOM6 `zonal_flux_adjust` behaviour at donor flips. | `test_continuity_multilayer`, `test_ocean_step_bt_mode` |

Two things a reader will get wrong. First, the ocean ALE remap does honour
`remap_method`; the horizontal advection edge scheme is chosen by `tracer_advection_scheme`, not `remap_method`. Second, there are two
separate PPM implementations (ALE remap and continuity) on purpose; each fits its
own loop and conservation constraint.

The windowed tracer-advection path (`&ocean_vmix_nml dt_tracer_advect_ratio`) is
not tabulated here yet. See the TODO comment at the end of
[section 16](#16-tunable-parameters).

**Full notes in the archive:** [reconstruction primitives with file:line pointers, regime matrix, knob declarations, honesty notes](CLOSURE_MATRIX_ARCHIVE.md#a7)

## 6. Equation of state

| EOS | Status | Knob | Notes |
|---|---|---|---|
| Linear two-tracer | **default** | `&ocean_eos_nml eos = "linear"` | Uses the `&ocean_ic_nml` coefficients below. |
| Wright (1997) | available | `eos = "wright"` | Nonlinear. |
| Roquet et al. (2015) SpV | available | `eos = "roquet_spv"` | TEOS-10 class: 75-term specific-volume fit, ΔS = 24. Not with the `fv_wright` PGF. |
| Roquet et al. (2015) in-situ density | available | `eos = "roquet_rho"` | TEOS-10 class: 55-term density fit, ΔS = 32 (OM3's `ROQUET_RHO`). An independent fit from `roquet_spv`: they agree to ~1e-3 kg/m³, not bitwise. Not with the `fv_wright` PGF. The FV-MOM6 in-situ density PGF (`&ocean_pgf_nml insitu_density`, the default) has its own factored roquet_rho PCM + reconstruction kernels (Adcroft, Hallberg & Harrison 2008 analytic-FV form), so the combination VALIDATES and runs the in-situ path by default; `insitu_density = .false.` still opts into the potential-density integral explicitly. `test_ocean_eos_roquet_rho`, `test_ocean_pgf_roquet_rho`, `test_ocean_pgf_insitu :: thermobaric_face_matches_insitu_oracle_roquet_rho` / `rest_partial_steps_zero_pgf_roquet_rho`, `test_config :: roquet_rho_insitu_pgf_accepted` |

Both Roquet variants take the model's potential temperature and practical
salinity and convert them to CT and SR with the GSW-verified
`rdb_teos10_convert` (`roquet_model_ct`). **Answers changed 2026-10-13:** the
old conversion misplaced two polynomial terms. Over SP 0–42, θ −2 to 32 °C and
p 0/1000/5000 dbar the fix moves CT by up to 2.3e-2 °C and ρ by up to
4.5e-3 kg/m³. No regression golden selects either variant.

Linear-EOS reference state: `ρ = ρ_0 + β_S·(S−S_ref) − α_T·(T−T_ref)`.
`α_T` and `β_S` are **dimensional** (kg/m³ per unit). Multiply a fractional
coefficient by `ρ_0` first.

| Parameter | Knob | Default | Units | Test |
|---|---|---|---|---|
| `α_T` | `&ocean_ic_nml alpha_T` | `1.7e-4` | kg/m³ per °C | `test_ocean_linear_eos_knobs` |
| `β_S` | `&ocean_ic_nml beta_S` | `7.6e-4` | kg/m³ per PSU | `test_ocean_linear_eos_knobs` |
| `T_ref` | `&ocean_ic_nml T_ref` | `10.0` | °C | `test_ocean_linear_eos_knobs` |
| `S_ref` | `&ocean_ic_nml S_ref` | `35.0` | PSU | `test_ocean_linear_eos_knobs` |
| `ρ_0` | `&ocean_ic_nml rho_0` | `1035.0` | kg/m³ | `test_ocean_linear_eos_knobs` |

| Related choice | Status | Knob (default) | What it does | Test |
|---|---|---|---|---|
| Source of α/β for KPP `B_0` and double diffusion | **default** `"constant"` | `&ocean_vmix_nml buoyancy_coeffs` (`"constant"` or `"eos"`) | `"eos"` takes α and β from the active EOS. Byte-identical to `"constant"` under the linear EOS. | `test_ocean_eos_buoyancy` |
| Potential-density reference pressure | opt-in | `&ocean_eos_nml p_ref` (`0.0`) | Pressure that `ms%rho_layer` is referenced to. Horizontally uniform by design. | `test_ocean_eos_p_top` |
| T/S tracer convention | **default** `"pot_prac"` | `&ocean_eos_nml ts_convention` (`"pot_prac"` or `"cons_abs"`) | `"cons_abs"` makes the T/S tracers Conservative Temperature and Absolute Salinity (IOC/SCOR/IAPSO 2010): both Roquet et al. (2015) fits take T, S directly (no SR factor, no PT→CT, unit Jacobian); the heat capacity does not follow the convention — it is `&ocean_eos_nml c_p` (MOM6 `C_P`, default 3992 = ACCESS-OM3's value under CT/SA; TEOS-10 cp0 = 3991.86795711963, manual eq. 3.3.3, is opt-in). Fails loud unless `eos` is `roquet_rho`/`roquet_spv`, and with cavity melt. Under `cons_abs`: a file `&ocean_zinit_nml` holding θ/SP (`input_convention = "pot_prac"`, the default, as WOA) is converted on the host after interpolation, SA first (`SA = sr_from_sp(SP)`, Reference Salinity, no δSA) then `CT = ct_from_pt(SA, θ)`; `input_convention = "cons_abs"` reads CT/SA as is and is refused under a `pot_prac` model. `&ocean_bulk_flux_nml` evaluates its formulae and longwave emission at `θ = pt_from_ct(SA, CT)` of the top layer (heat content of precipitation keeps the model T with `c_p`). `&ocean_restore_nml` `restore_sst`/`restore_sss` are read in the model convention (OM3 restores toward a variable that is already Absolute Salinity). Derived diagnostics `temperature_pt` (θ, CMOR `thetao`) and `salinity_sp` (SP, CMOR `so`) are identity copies under `pot_prac`; under `cons_abs` `temperature`/`salinity` are CMOR `bigthetao`/`absso`. **Sea ice is allowed:** the coupler boundary converts once (ice side θ/SP: `sst_seam = pt_from_ct`, `ssurf_seam = sp_from_sr`, `tfw_seam = pt_from_ct(SA, CT_f)`; ocean side CT/SA: basal-flux available heat, frazil onset and heat, with `c_p` from the knob). `c_p` is bounded to 3000-4500. | `test_ocean_eos_roquet`, `test_ocean_freezing_point`, `test_ocean_zinit`, `test_ocean_diag`, `test_ocean_bulk_forcing`, `test_ocean_frazil`, `test_ocean_ice_coupling` |
| Liquidus form | **default** `"linear"` | `&ocean_eos_nml tfreeze_form` (`"linear"` or `"teos_poly"`), `tfreeze_sat_frac` | `"teos_poly"` is the TEOS-10 polynomial `CT_f(SA, p, sat_frac)` (manual section 3.33; McDougall & Barker 2011), about −1.92 °C at SA = 35.16504, p = 0. Needs `ts_convention = "cons_abs"`. `tfreeze_sat_frac` (default 0 = no dissolved air) is ACCESS-OM3 parity (`TFREEZE_FORM = "TEOS_POLY"`). | `test_ocean_freezing_point` |
| Freezing-point coefficient set | **default** `"seaice"` | `&ocean_eos_nml tfreeze_set` (`"seaice"` or `"isomip"`) | Named liquidus sets. They differ by about 0.031 °C at S = 34.5, enough to flip the sign of an ice-shelf melt rate. | `test_ocean_freezing_point` |
| Top-of-column pressure in the EOS | opt-in | `&ocean_psurf_nml in_eos` (`.false.`) | Evaluates in-situ density at `p_top + hydrostatic`. | `test_ocean_eos_p_top`, `test_ocean_bl_under_ice` |

`validate_config` **warns** on cavity melt with a nonlinear EOS and
`buoyancy_coeffs="constant"`. A known, out-of-scope KPP `B_0` units issue is
recorded in the archive.

**Full notes in the archive:** [where α and β come from, consumer table, measured α values, known KPP units issue](CLOSURE_MATRIX_ARCHIVE.md#a8)

## 7. Vertical mixing

Vertical mixing only runs in 3D multilayer mode (it is a no-op at `nz=1`). Every
interior or boundary-layer closure adds its diffusivity into `kv` (momentum) and
`kt` (temperature). The salt diffusivity `ks` is derived from `kt` afterwards.
One assembly step applies floors, ceilings, smoothing and an optional guard.
The solve is implicit (backward-Euler tridiagonal).

| Closure | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| PP81 (Pacanowski-Philander) | **default** | `&ocean_vmix_nml use_closure` (`.true.`), `interior_closure` (`"pp81"`) | Richardson-number interior mixing (Pacanowski & Philander 1981). `interior_closure = "none"` drops it: the interior `kv`/`kt` are seeded each stage with the background (scalar, Henyey or Bryan-Lewis; `prandtl·kt + kv_bg` for momentum when `bkgnd_prandtl > 0`) and the other contributors add on top, as on the PP81 path (ACCESS-OM3 runs with no PP81). `use_closure` stays the master switch for the whole pipeline. | `test_ocean_pp81`, `test_ocean_vmix_no_pp81` |
| KPP (Large et al. 1994) | **default** | `&ocean_vmix_nml use_kpp` (`.true.`) | Boundary-layer overlay on PP81 with non-local transport. Shortwave method `&ocean_thermo_nml kpp_sw_method` (`"mxl_sw"`). | `test_ocean_kpp*`, `test_ocean_buoyancy_flux`, `test_ocean_sw_bl_coupling` |
| EPBL (Reichl & Hallberg 2018) | opt-in; replaces KPP | `&ocean_epbl_nml enable` (`.false.`) | Energetics-based boundary layer with optional Langmuir term. Enabling it logs and disables `use_kpp`. | `test_ocean_epbl`, `test_ocean_sw_bl_coupling` |
| kappa-shear (Jackson, Hallberg & Legg 2008) | opt-in | `&ocean_kappa_shear_nml enable` (`.false.`) | Prognostic interior shear turbulence, added to the other closures. `at_vertex` selects the corner (MOM6 `VERTEX_SHEAR`) form. | `test_ocean_kappa_shear` |
| Tidal mixing (St Laurent / Simmons) | opt-in | `&ocean_tidal_mixing_nml enable` (`.false.`) | Bottom-intensified internal-tide diffusivity (St Laurent exponential profile). With `e_compute`, the energy input is `min(c_e·N_b, e_max)` (Jayne & St Laurent 2001); `tideamp_source`/`h2_source = "file"` read static 2-D `tideamp` and `h2` fields and fold the mask, `min_zbot` and roughness clamp into `c_e` at setup. OM3-parity refinements, all default to the legacy behaviour: `n2_bot_rough_avg` (near-bed `N²` averaged over the roughness height, Melet et al. 2013), `energy_limiter = "reinject"` (layer TKE capped at the `kd_max` equivalent, excess passed up; default `"clamp"` discards it), `bed_share = "lump"` (bed-layer power deposited in the first interior layer; default `"discard"`). `n2_layer_form = "mean_interfaces"` (layer `N²` = mean of the two bounding interface jumps over the thickness; default `"legacy"`, always used under Polzin). `profile = "polzin"` swaps the exponential for the Polzin (2009) algebraic profile in WKB-stretched height (`Φ = P·I_p·z_p*/(z_p*+z*)`, scale from `U`, `⟨h²⟩`, `N_b`, `N̄²`; Melet et al. 2013), through the same deposit/limiter sweep; constants `nu_polzin`, `nbotref_polzin`, `polzin_decay_scale_factor`, `polzin_scale_max_factor`, `polzin_min_decay_scale`; needs `e_compute`, default `"stlaurent"`. | `test_ocean_tidal_mixing`, `test_static_2d_reader`, `test_ocean_zfixed_k_bot` |
| Law-of-the-wall BBL mixing | opt-in | `&ocean_bbl_mixing_nml enable` (`.false.`) | Bottom-drag-driven tracer diffusivity. The BBL mean-KE loss `P_bbl` (per-face `u*` and BBL-mean along-face velocity, from the viscous BBL) times `effic` is spent along a column sweep from the bed: wall-law `Kd = κu*²z(D−z)/(u*D + \|f\|z(D−z))` (Killworth & Edwards 1999 rotation scale), power decaying upward with `1/λ = max(1/max_decay, \|f\|/u*)`, `Kd·Δz·max(N², Ω²)` bounded by the power left; unstratified interfaces fill to `kd_max`. Additive to `kt` only (the BBL viscosity stays `kv_bbl`); `use_omega`, `von_karman`, `kd_max` expose the OM3 values (`effic = 0.01`). Needs the per-face BBL (`bbl_glue` + bottom drag `hbbl > 0`), `c_d > 0`, thermodynamics; fails loud otherwise. Interface `N²` is the tidal kernel's `column_layer_n2`. Diagnostic `Kd_BBL`. | `test_ocean_bbl_mixing` |
| Penetrating shortwave (two-band, Paulson & Simpson 1977) | **default** (inert at `sw_pen_frac = 0`) | `&ocean_thermo_nml sw_pen_frac` (`0.0`), `sw_source` (`"net_heat"`), `sw_band_ratio`, `sw_zeta1`, `sw_zeta2` | Two-band exponential decay of a fraction of the surface heat. One transmission function is shared by the deposition kernel, the KPP `B_0` coupling (`kpp_sw_method`) and the EPBL SW TKE ledger (`epbl_sw_ctke`). | `test_ocean_sw_penetration`, `test_ocean_sw_bl_coupling` |
| Penetrating shortwave (chlorophyll, Manizza et al. 2005) | opt-in | `&ocean_thermo_nml sw_scheme` (`"default"`; `"manizza05"`), `sw_chl_file`, `sw_chl_var` (`"CHL"`), `sw_nbands` (`3`), `sw_blue_frac` (`0.5`), `sw_opacity_coef`, `sw_chl_power`, `sw_opacity_land` (`10.0`) | Replaces the two-band profile: the WHOLE shortwave is split 0.42 visible (blue/red by `sw_blue_frac`) and 0.58 near-infrared, and each band decays as `exp(-kappa_n z)` with a per-column `kappa_n = a_n + b_n chl^p_n` from a monthly-cyclic chlorophyll climatology. Deposition, KPP `B_0` and the EPBL ledger consume the same per-column band opacities. Refuses a `sw_pen_frac` other than 0 or 1, and a missing file or variable. The flux-absorption floor (`PEN_SW_FLUX_ABSORB`) is not applied. | `test_ocean_sw_chl` |
| Convective adjustment | opt-in | `&ocean_conv_nml enable` (`.false.`; `kd_conv` `1.0`, `n2_thresh` `0.0`) | Raises `kt` and `kv` where `N² < n2_thresh`, below the active boundary layer. Never writes `ks`. | `test_ocean_convection` |
| Double diffusion | opt-in | `&ocean_ddiff_nml enable` (`.false.`) | Salt fingering plus diffusive convection. Folded into the heat/salt split, so `ks` and `kt` differ. | `test_ocean_double_diffusion` |
| Background mixing | opt-in; the two forms are mutually exclusive | `&ocean_vmix_nml bkgnd_profile` **xor** `bkgnd_henyey` (both `.false.`) | Bryan-Lewis depth profile, or a Henyey latitude factor on the scalar floors (Henyey et al. 1986; Harrison & Hallberg 2008). Henyey is refused on a cartesian grid. `bkgnd_prandtl > 0` (default `-1`, off) ties the interior momentum floor to the tracer background in use, `kv ≥ prandtl·kt_background + kv_bg`, so the latitude factor and `bkgnd_kd_min` reach momentum (OM3: `bkgnd_prandtl = 1`, `pp81_nu_bg = 0`). | `test_ocean_bkgnd_mixing`, `test_ocean_vmix_no_pp81` | `test_ocean_bkgnd_mixing` |
| Assembly gate (`vmix_assemble`) | **default** (inert) | `kv_max`, `kd_max` (`huge`), `kd_smooth_iterations` (`0`), `vmix_guard` (`.false.`) | Single floors, ceilings, smoothing and guard stage. All defaults are bit-identical. | `test_ocean_vmix_assembly`, `test_ocean_vdiff` |

Not built: the two-equation k-ε closure and the Galperin/Canuto stability
functions (they left with the coastal path), and the GLS-generic family.

The `Double diffusion` row is not in the old table, and the old text said no
double-diffusion contributor was shipped. The code disagrees: see
`src/core/rdb_config.F90:2129` and `tests/CMakeLists.txt:192`.

**Full notes in the archive:** [full vertical-mixing table, diffusivity assembly (C11), MOM6 mappings](CLOSURE_MATRIX_ARCHIVE.md#a9)

## 8. Lateral closures and eddy parameterizations

This section covers horizontal viscosity, mixed-layer restratification, eddy
thickness diffusion (GM), neutral tracer diffusion (Redi) and the eddy-energy
budget (MEKE).

| Closure | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Constant ν_h floor | **default** | `&ocean_hvisc_nml nu_h` (`0.0`) | Scalar harmonic viscosity. | — |
| Smagorinsky Laplacian | opt-in | `lateral_closure = "smag"` (`"none"` default; `c_smag` `0.15`) | Strain-rate-scaled viscosity per face. | — (exercised by the multilayer tests) |
| Leith | opt-in | `lateral_closure = "leith"` (`c_leith` `1.0`) | Vorticity-gradient viscosity. | `test_ocean_leith` |
| Biharmonic, constant | opt-in | `nu_4` (`0.0`) | ∇⁴ with constant ν₄. CFL-capped. | `test_ocean_hvisc`, `test_ocean_hvisc_biharm_zfixed` |
| Smagorinsky biharmonic | opt-in | `smag_ah` (`.false.`) | Flow-aware biharmonic. | `test_ocean_smag_ah`, `test_ocean_hvisc_biharm_zfixed` |
| Leith biharmonic | opt-in | `lateral_closure = "leith_biharm"` (`c_leith_bi` `0.0`) | MOM6 `LEITH_AH`. | `test_ocean_hvisc_leith_biharm` |
| Resolution-scaled viscosity | opt-in | `resoln_scaled_visc` (`.false.`; needs `&ocean_varmix_nml enable`) | Hallberg 2013 scaling. | `test_ocean_hvisc_resoln` |
| Point Laplacian floors (stress path) | opt-in | `kh_bg_at_points` (`.false.`; needs `stress_tensor`) | `max(nu_h, kh_vel_scale·Δ, kh_sin_lat·\|sin φ\|^p)` evaluated AT the T centre and AT the Bu corner (local latitude, local harmonic `Δ`), `max`-combined into the `ah_t`/`ah_q` seed; not a face mean (the sin floor is convex). Statics `kh_bg_t`/`kh_bg_q` built at configure. | `test_ocean_hvisc_kh_points` |
| Latitude ν_h floor | opt-in | `kh_sin_lat` (`0.0`), `kh_pwr_of_sine` (`4.0`) | `max(floor, kh_sin_lat·\|sin φ\|^p)` at face-centre latitude, max-combined into `ah_face_*` under the closures (MOM6 `KH_SIN_LAT`). Refused on a cartesian grid. | `test_ocean_hvisc_om3_floors` |
| Biharmonic velocity-scale floor | opt-in | `ah_vel_scale` (`0.0`; needs `smag_ah`) | Background `max(nu_4_bg, U·(Δ²)^(3/2))`, `Δ²` the harmonic mean of the T-cell dx², dy² (MOM6 `AH_VEL_SCALE`). | `test_ocean_hvisc_om3_floors` |
| Coriolis bound on Smagorinsky biharmonic | opt-in | `bound_coriolis_biharm` (`.false.`; needs `smag_ah`), `bound_coriolis_vel` (`6.0`) | `A = \|S\|(C_b·Δ⁴ + Δ⁶·f_max·\|S\|/(5V²))`, `f_max` the largest `\|f\|` of the T cell's four corners (MOM6 `BOUND_CORIOLIS_BIHARM`; Griffies & Hallberg 2000). | `test_ocean_hvisc_om3_floors` |
| Live velocity-scale ν_h | opt-in | `kh_vel_scale_live` (`0.0`) | `Kh = U·dx·\|u\|`. | `test_ocean_hvisc_aniso` |
| Anisotropic ν_h | opt-in | `kh_aniso` (`0.0`) with `aniso_dir`; stress-tensor path only | Smith & McWilliams 2003. | `test_ocean_hvisc_aniso` |
| Stress-tensor Laplacian | opt-in | `stress_tensor` (`.false.`) | Thickness-weighted stress-divergence form, free-slip walls. Refused with `zfixed_closed_faces`. | — |
| Per-face harmonic clamp | opt-in | `bound_kh` (`.false.`; `bound_coef` `0.8`) | MOM6 `BOUND_KH`. Load-bearing on the split-explicit path. | `test_ocean_hvisc_kh_bound` |
| Fox-Kemper ML restratification | opt-in | `&ocean_foxkemper_nml enable` (`.false.`; needs EPBL) | Mixed-layer eddy overturning. | `test_ocean_foxkemper`, `test_ocean_mle_zfixed` |
| Gent-McWilliams thickness diffusion | opt-in | `&ocean_gm_nml enable` (`.false.`; needs `&ocean_slopes_nml enable`) | Eddy-induced bolus transport. `use_fgnv` (`.false.`; needs `&ocean_wavespeed_nml enable`) replaces the local streamfunction by the Ferrari, Griffies, Nurser & Vallis (2010) vertical boundary-value solution `c^2 psi'' - N^2 psi = -N^2 psi_unlim`, `psi = 0` at surface and bed, with `psi_unlim` amplified by `(1+fgnv_filter_scale)`, `c^2 = fgnv_filter_scale * cg1^2` (`fgnv_c_min` accepted for parity, not applied) and `N^2 >= (fgnv_strat_floor*Omega)^2`; the limiter chain runs unchanged on the solved psi. `gm_src` (MEKE's GM source) is the local-form PE release `rho0/4 Sum KH S^2 N^2 h` unless `gm_src_form = "work"` (`"local"` default; OM3 `MEKE_GM_SRC_ALT = False`): then it is the potential-energy release rate of the FINAL limited bolus transports (FGNV, EBT structure and limiters included), per face from the layer and interface terms `g[Sum_k uhD_k dRho_k ebar_k + Sum_K U_K (D_R - D_L)]` (`U_K` the transport above interface K, `D = (rho_below - rho_above) e_K`), averaged to cells as `0.5/A` of the four faces, computed after `uhD` is final and before it is applied; Boussinesq `g`, T/S from the slopes slot's thermo-cadence copies. | `test_ocean_gm`, `test_ocean_gm_zfixed`, `test_ocean_slopes_datum`, `test_ocean_gm_fgnv`, `test_ocean_gm_work` |
| GM vertical structure of KhTh (MOM6 `KHTH_USE_EBT_STRUCT`) | opt-in | `&ocean_gm_nml use_ebt_struct` (`.false.`; needs GM and `&ocean_wavespeed_nml enable` + `use_ebt`, refused otherwise) | At each interior interface the post-CFL-clamp 2-D face KhTh (base, Visbeck and MEKE together) is multiplied by the mean of the two columns' EBT structure in the layer above the interface (Hallberg 2013; Killworth & Blundell 2003), in the streamfunction and in the `gm_src` that MEKE reads. Structure 0 (degenerate column) gives KhTh 0, not floored. | `test_ocean_gm` |
| Redi neutral tracer diffusion | opt-in | `&ocean_redi_nml enable` (`.false.`), `check_diffusive_cfl` (`.false.`), `max_diffusive_cfl` (`-1`, off) | Along-isopycnal diffusion (Griffies et al. 1998; Shao et al. 2020 operator). Refused under a cavity. The two CFL knobs are the OM3 `CHECK_DIFFUSIVE_CFL` / `MAX_TR_DIFFUSION_CFL` (`True` / `2.0`): `max_diffusive_cfl > 0` caps each face coefficient so the cell CFL `2·Σ KhTr·dt·dy/dx / area` stays at or below it, and `check_diffusive_cfl` subcycles the flux `n_itt = ceil(max CFL)` times at `dt/n_itt`. Both off is bit-identical. Neutral diffusion replaces the along-coordinate Laplacian: OM3 parity is `&ocean_hdiff_nml kappa_h = 0`, and `kappa_h > 0` with Redi on warns at configure. | `test_ocean_redi`, `test_ocean_redi_zfixed` |
| MEKE | opt-in | `&ocean_meke_nml enable` (`.false.`; needs GM) | Prognostic eddy kinetic energy. | `test_ocean_meke` |
| MEKE backscatter | opt-in | `&ocean_meke_nml backscatter` (`.false.`) with `backscatter_visc_coeff_ku` | Negative-viscosity energy return. Needs a flow-aware closure and a biharmonic backstop. | `test_ocean_meke_backscatter` |
| Along-coordinate tracer Laplacian | opt-in | `&ocean_hdiff_nml kappa_h` (`0.0`) | Diffusion along the model coordinate, not neutral surfaces. | — |

`lateral_closure` is fail-loud: an unknown or unimplemented tag aborts at
configure. The knob spellings for the Smagorinsky Laplacian (`lateral_closure =
"smag"`) and the stress-tensor form (`stress_tensor`) were corrected from the old
table; the evidence is in the commit message.

**Full notes in the archive:**

- [viscosity notes (CFL cap, BOUND_KH, free-slip, biharmonic composition)](CLOSURE_MATRIX_ARCHIVE.md#a10)
- [full original closure table](CLOSURE_MATRIX_ARCHIVE.md#a11)
- [Fox-Kemper restratification (B5), including partial-step faces](CLOSURE_MATRIX_ARCHIVE.md#a12)
- [Gent-McWilliams (GM) and Redi, including partial-step faces and the MOM6 divergence note](CLOSURE_MATRIX_ARCHIVE.md#a13)
- [MEKE and MEKE backscatter](CLOSURE_MATRIX_ARCHIVE.md#a14)

## 9. Drag and surface stress

Three momentum sinks and sources live here: the wind stress, the bottom drag and
the ice-shelf top drag. The implicit folds put stress or drag into the
backward-Euler vertical-friction solve so thin layers stay stable.

### Wind stress

| Mode | Status | Knob (default) | What it does | Test |
|---|---|---|---|---|
| Top-layer stress | **default** | `&ocean_vmix_nml direct_stress` (`.false.`) | Stress lands in the surface layer. Generates Ekman and mixed-layer shear. | — |
| Distributed stress | opt-in | `direct_stress = .true.` (`hmix_stress` `20.0` m) | Spreads stress over the top `hmix_stress` metres. | — |
| Implicit fold | opt-in | `&ocean_vdiff_nml implicit_stress` (`.false.`) | Folds stress into the vdiff surface row. Incompatible with `direct_stress`. | `test_ocean_vdiff_implicit_stress_drag` |

### Bottom drag

`form` defaults to `"quadratic"`, but `cd`, `r` and `hbbl` all default to `0.0`
(`src/core/rdb_config.F90:1229` to `1237`). A coefficient of zero switches that
branch off, so no drag acts until you set `cd` or `r`.

| Form | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Quadratic | form **default** | `&ocean_bdrag_nml form = "quadratic"`, `cd` (`0.0`, typical `2.5e-3`) | Log-layer drag `-C_d·\|U\|·u/h`. | — |
| Linear (Rayleigh) | available | `form = "linear"`, `r` (`0.0`, 1/s) | `-r·u`. | — |
| Bed layer = first live layer | **default** | none | Every bottom-drag form sits on the first live layer (`kb_live`), not on a vanished sliver. | `test_ocean_zfixed_k_bot`, `test_ocean_bottom_drag_live_kbot` |
| HBBL-distributed | opt-in; required with `zstar_full` | `&ocean_bdrag_nml hbbl` (`0.0` m) | Spreads stress over the bottom `hbbl` metres. | — |
| Split-apply implicit | opt-in | `&ocean_bdrag_nml implicit` (`.false.`) | Backward-Euler drag. | — |
| Implicit fold into vdiff | opt-in | `&ocean_vdiff_nml implicit_drag` (`.false.`) | Drag on the vdiff bed-row diagonal. Exclusive with `implicit`. | `test_ocean_vdiff_implicit_stress_drag` |
| MOM6 `BOTTOMDRAGLAW` (BBL glue) | **default** (inert with no drag set) | `&ocean_vdiff_nml bbl_glue` (`.true.`; needs `hvel_mom6`) | Moves the drag into the vertical viscosity (MOM6 `set_viscous_BBL`). Replaces the explicit, split-apply and fold drag on the layers. | `test_ocean_vdiff_bbl` |
| MOM6 momentum face thickness | **default** | `&ocean_vdiff_nml hvel_mom6` (`.true.`) | MOM6 `hvel` and `h_shear` in the friction solve. | `test_ocean_vdiff_bbl` |

The `bbl_glue` and `hvel_mom6` defaults were `.false.` in the old tunable table.
The code says `.true.`: `src/core/rdb_config.F90:1601` and `:1638`.

### Ice-shelf top drag (`&ocean_tdrag_nml`)

Needs `&ocean_cavity_dyn_nml enable`. A face is under ice if **either** neighbouring
cell is, so the calving-front face is dragged. With `&ocean_cavity_melt_nml` on,
the melt and the momentum drag share one `C_d`; a disagreeing `cdrag_top` is
refused. `cd` defaults to `0.0` (`src/core/rdb_config.F90:1276`), so set it
(`2.5e-3` is the ISOMIP+ value) for a quadratic drag.

| Form | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Master switch | opt-in | `&ocean_tdrag_nml enable` (`.false.`) | Mirror of the bottom drag at `k=nz`. | `test_ocean_top_drag` |
| Quadratic | form **default** when enabled | `form = "quadratic"`, `cd` | ISOMIP+ log-layer drag. | `test_ocean_top_drag` |
| Linear | available | `form = "linear"`, `r` | Rayleigh drag with a closed-form spin-down. | `test_ocean_top_drag` |
| HTBL-distributed | opt-in | `htbl` (`0.0` m) | Spreads stress over the top `htbl` metres. | `test_ocean_top_drag` |
| Implicit in the drag kernel | opt-in | `implicit` (`.false.`) | Backward-Euler form of the drag. | `test_ocean_top_drag` |
| Implicit fold into vdiff | opt-in | `&ocean_vdiff_nml implicit_top_drag` (`.false.`) | Drag on the `k=nz` diagonal; wind is masked under cover. Exclusive with `implicit` and `htbl > 0`. | `test_ocean_top_drag` |
| Stress magnitude `stress_top` | **default** when enabled | none | Feeds `stress_shelf`, which KPP and EPBL read. | `test_ocean_bl_under_ice` |
| Under-ice `u_*` for KPP and EPBL | **default** when top drag or melt is on | `&ocean_tdrag_nml enable`, else `&ocean_cavity_melt_nml enable` | `u_* = √((stress_mag + stress_shelf)/ρ₀)`. | `test_ocean_bl_under_ice` |

**Full notes in the archive:**

- [wind-stress table (original cells)](CLOSURE_MATRIX_ARCHIVE.md#a15)
- [bottom-drag table (original cells, BBL glue formulas, kb_live rule)](CLOSURE_MATRIX_ARCHIVE.md#a16)
- [ice-shelf top drag (face rule, one-C_d rule, original cells)](CLOSURE_MATRIX_ARCHIVE.md#a17)

## 10. Subgrid topography: porous barriers and partial-step faces

Two independent, per-layer face gates. **Porous barriers** narrow a face
continuously for a sill the grid cannot resolve. **Partial-step closure** shuts a
face completely for a layer the coordinate has already declared dead. The
effective width is the product `dy_eff = dy_cu · por_face_area_u · open_u`.

| Scheme | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Porous barriers (Adcroft 2013) | opt-in | `&ocean_porous_nml enable` (`.false.`) | Per-layer open-area fraction from a three-parameter along-face fit; narrows continuity and transport-Coriolis fluxes, and the barotropic widths. Refused with `bt_halo > 0` and with wet/dry. Multi-rank bit-identical. | `test_ocean_porous` |
| Porous data source | **default** `"resolved"` | `source` (`"resolved"`; `"file"` is refused) | `"resolved"` samples the resolved bathymetry. It is a proxy that sees only variation along the face. | `test_ocean_porous` |
| Porous interface height | **default** `"max"` | `eta_interp` (`"max"`, `"min"`, `"arithmetic"`, `"harmonic"`) | `"max"` blocks the least. | — |
| Porous masking depth | **default** | `masking_depth` (`0.0` m) | Faces shallower than this stay fully open. | — |
| Partial-step z-level face closure (Adcroft, Hill & Marshall 1997; Losch 2008) | opt-in; **required** over a stepped `z_fixed` bed | `&vcoord_nml zfixed_closed_faces` (`.false.`) | A face where a layer is a filler on either side becomes a free-slip wall for that layer. Works with `z_fixed`, `zstar`, `zstar_full`. Single-rank. | `test_ocean_zfixed_closed_faces`, `test_ocean_zfixed_bt_seiche`, `test_ocean_zfixed_cor_ref`, `test_ocean_zstar_full_closed_faces`, `test_ocean_zstar_closed_faces`, `test_ocean_bt_upstream_zfixed`, `test_ocean_hvisc_biharm_zfixed`, `test_ocean_gm_zfixed`, `test_ocean_mle_zfixed`, `test_ocean_coriolis_hk_vanished` |

With `zfixed_closed_faces` off, `z_fixed` is **refused** wherever a wet face joins
two columns whose bed lies in different nominal layers. Flat or step-free beds
stay accepted. With the knob on, these stay **refused**: wet/dry, `bt_halo > 0`,
`stress_tensor` (and so `kh_aniso`), `correction_bc_pgf`, and `pred_corr` with
`&ocean_bc_nml mask_wall_velocity = .false.`.

Seven consumers change when the knob is on: continuity, the barotropic
renormaliser, transport Coriolis, the layer-velocity mask, `kappa_h`, the
harmonic and biharmonic viscosity (free-slip), and the momentum vertical
friction. GM, Redi and Fox-Kemper compose with it through each face's open
column.

**Full notes in the archive:**

- [porous barriers, proxy caveats, multi-rank behaviour](CLOSURE_MATRIX_ARCHIVE.md#a18)
- [partial-step closure (mask construction, every consumer, exclusions, measurements)](CLOSURE_MATRIX_ARCHIVE.md#a19)

## 11. Surface forcing, initial conditions and other fluxes

### Initial conditions (analytical seeds)

| Seed | Status | Knob | Test |
|---|---|---|---|
| Uniform T and S | **default** | `&tracer_nml initial_temperature`, `initial_salinity` | `test_ocean_salinity_ic` |
| Linear-in-layer T(z) (bed `k=1` to surface `k=nz`) | available | `&tracer_nml T_init_bottom` and `T_init_surface` (both non-zero) | `test_ocean_eady_ic` |
| Linear-in-layer S(z) | available | `&tracer_nml S_init_bottom` and `S_init_surface` (both non-zero; one alone fails loud) | `test_ocean_salinity_ic` |

### Boundary and surface fluxes

| Flux | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Surface heat and salt (top layer) | scalar | `&ocean_thermo_nml q_heat`, `q_salt` | Uniform scalar fluxes. | `test_ocean_surface_flux` |
| Flux component set | opt-in | `&ocean_forcing_nml enable_components` (`.false.`) | Per-component 2D fields (`q_sw`, `q_lw`, mass fluxes, `p_surf`, ...). `Q_heat` and `Q_salt` become derived views. | `test_ocean_surface_forcing_type` |
| Bulk air-sea fluxes off a prescribed atmosphere (JRA55-do) | opt-in, ocean + sea ice | `&ocean_atm_forcing_nml enable` + `&ocean_bulk_flux_nml enable` (`.false.`, required together) | Large & Yeager (2004, 2009) column kernel over the `atm_state_t` nine-field provider: wind rotated onto the grid axes (`angle_dx`), relative wind against the live surface current, writes `q_sens`/`q_lat`/`evap`/`q_sw`/`q_lw`/`lprec`/`fprec`/`p_surf_atm` as components, open-water-weighted (`ice%part_size(:,:,0)`) when sea ice is enabled; stress averaged onto the C-grid faces and seam-refreshed. Sea-ice seam (`ice_atm_forcing_from_atm`, OM3 C2c-3b): a SEB linearisation at `T_s = 0 degC` with Large & Yeager's constant over-ice transfer coefficients (`C_D = C_H = C_E = 1.63e-3`, no Monin-Obukhov iteration) feeds `atm_sf0`/`atm_dsfdt`; `atm_sw_dn = rsds`, `atm_fprec = prsn`; dispatched instead of `ice_atm_forcing_restoring` whenever atm forcing is on. | `test_ocean_bulk_forcing`, `test_ocean_bulk_flux`, `test_ocean_atm_state`, `test_ocean_atm_state_mpi` |
| Geothermal bottom heat | opt-in | `&ocean_geothermal_nml enable` (`.false.`), `q_geo` (`0.0` W/m²) | Heat into the lowest massive layer. | `test_ocean_geothermal` |
| Surface buoyancy restoring | opt-in | `&ocean_restore_nml enable_restore_temp`, `enable_restore_salt` (`.false.`) | Piston-velocity relaxation of top-layer T and S toward scalar targets. Non-conservative. | `test_ocean_restore` |
| File-based SSS restoring (OM3 `RESTORE_SALINITY`, `SALT_RESTORE_FILE`) | opt-in | `&ocean_restore_nml sss_file` (`""`), `sss_var`, `sss_time_cyclic` (`.true.`), `sss_cycle_period`, `max_delta_sss` (`999`), `sss_adjust_net_zero`, `sss_adjust_by_scaling` (`.false.`), `sss_as_sflux` (`.false.`, must be `.true.` with a file), `sss_convention` (`"model"`) | Salt-FLUX form (Griffies et al. 2009, CORE salinity restoring): `F = rho0 w dS` with `dS` clamped to `max_delta_sss`, in the `Q_salt` unit convention (g/kg x kg/m^2/s, i.e. 1e3 x the SI salt mass flux), deposited in the first live layer and mirrored into the salt budget contributor. `sss_adjust_net_zero` removes the global net: with `sss_adjust_by_scaling` the larger-magnitude side of `P = sum A max(F,0)` / `N = sum A min(F,0)` is scaled by `min(P,-N)/max(P,-N)` (zero contour unmoved), otherwise the wet-area mean is subtracted; the sums are EFP fixed point plus one packed allreduce, so they are decomposition independent. The target is blended off the shared reader (cyclic climatology or linear), seam-ghost exchanged, and under `cons_abs` with `sss_convention = "practical"` converted per cell with `sr_from_sp`. Replaces the scalar content-injection gate (no virtual-precipitation form exists; the salt-flux form is the only one). Cells with a non-finite target or a covered (cavity) column do not take part. Refused without `enable_restore_salt` and `piston_s > 0`, with `sss_as_sflux = .false.`, or when the clamp/net knobs are set without a file. | `test_ocean_sss_restore` |
| Tidal body forcing | opt-in | `&ocean_tides_nml enable` (`.false.`) | Equilibrium tide in the barotropic PGF. Ten constituents. Needs a lat/lon grid. | `test_ocean_tides_astronomy`, `test_ocean_tidal_forcing` |
| Barotropic linear wave drag | opt-in | `&ocean_bt_nml wave_drag` (`.false.`), `wave_drag_form` (`"uniform"`) | Static per-face piston velocity multiplied into `bt_rem`. Egbert & Ray 2001. | `test_ocean_wave_drag`, `test_ocean_bt_upstream_zfixed`, `test_ocean_bt_substep_drag` |
| Surface-pressure loading (inverse barometer) | opt-in | `&ocean_psurf_nml enable` (`.false.`) | `η_ib = −p_surf/(ρ₀·g_bt)` on the barotropic `eta_forcing` seam. Split solver only; needs `enable_components`. | `test_ocean_p_surf`, `test_ocean_p_surf_bitident` |
| Reference density of forcing terms | **default** | `&ocean_ic_nml rho_0` (`1035`) | One `ρ₀` feeds every surface source, the wind stress and the geothermal flux. | `test_ocean_forcing_rho_ref` |
| EOS in-situ top pressure | opt-in | `&ocean_psurf_nml in_eos` (`.false.`) | See [section 6](#6-equation-of-state). | `test_ocean_eos_p_top`, `test_ocean_bl_under_ice` |
| Potential-density reference pressure | **default** | `&ocean_eos_nml p_ref` (`0.0`) | See [section 6](#6-equation-of-state). | `test_ocean_eos_p_top` |
| Freezing-point set | **default** `"seaice"` | `&ocean_eos_nml tfreeze_set` | See [section 6](#6-equation-of-state). | `test_ocean_freezing_point` |

**Full notes in the archive:**

- [initial-condition table (original cells)](CLOSURE_MATRIX_ARCHIVE.md#a20)
- [boundary tracer flux table (original cells, seams, budgets, refusals)](CLOSURE_MATRIX_ARCHIVE.md#a21)

## 12. Sea ice

`&ocean_ice_nml enable` gates the whole subsystem. Off (the default) means the
slot is never initialised, mapped or stepped, so results are byte-identical. This
is a SIS2 port on the ocean C-grid (`src/core/ice/`). It is a high-quality
sea-ice dynamical core and column model, not yet a full sea-ice model. Rows
marked `—` are documented absences. The authoritative knob set is the
`ocean_ice_config_t` block of `src/core/rdb_config.F90`.

Envelope: the slot runs on the ocean decomposition, bit-identical to one rank.
The exception is a tripolar fold on more than one rank, which is refused.
`transport` runs across periodic edges and MPI seams. EVP `dynamics` allows
wall or periodic edges, or a tripolar fold on one rank (needs `nghost >= 2`),
and no open, tidal, sponge, clamped or Chapman edge.

| Capability | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Winton (2000) column thermodynamics | opt-in | `&ocean_ice_nml enable` (`.false.`), `nk_ice` (`2`) | Two-layer column with enthalpy. | `test_ocean_ice_column`, `test_ocean_ice_enthalpy` |
| Multi-category ice-thickness distribution (ITD) | opt-in | `enable` with `ncat > 1` (`ncat` `5`) | Thickness-space category restore. | `test_ocean_ice_itd` |
| Category transport and `compress_ice` | opt-in | `transport` (`.false.`; needs `ncat > 1`), `adv_substeps` (`1`) | Category-summed PPM transport. | `test_ocean_ice_transport` |
| C-grid EVP dynamics | opt-in | `dynamics` (`.false.`), `evp_sub_steps` (`432`), `del_sh_min_scale` (`2.0`), `tdamp` (`-0.2`) | Elastic-viscous-plastic momentum. Also `p0`, `c0`, `ec`, `cdw`, `rho_ocean`, `a_face_stress`, `cfl_trunc`, `cfl_trunc_dyn_its`, `project_ci`. | `test_ocean_ice_evp`, `test_ocean_ice_evp_fold` |
| Ridging and rafting | **not available** | — | `compress_ice` is area compaction (SIS2 `DO_RIDGING=.false.`), not a ridging scheme. | — |
| Snowfall source | opt-in | `snowfall` (`0.0`; needs `enable` and `&ocean_thermo_nml enable_thermodynamics`) | Accumulates `m_snow`. | `test_ocean_ice_snowfall` |
| Snow-ice flooding | opt-in | `snow_ice` (`.false.`; same prerequisites) | One Archimedes step, SIS2 `SN2IC`. Exchanges nothing with the ocean. | `test_ocean_ice_column` |
| Melt ponds | **not available** | — | `m_pond` is a dead local. | — |
| Lateral melt and floe size | **not available** | — | Surface and basal melt only. | — |
| CICE-named diagnostics | with `enable` | `&ocean_diag_nml diags` | Window masses from the column outputs (`congel`, `meltt`/`meltb` ice-only, `snoice`, `frazil`); rates in cm/day. `dvidtt` closes against the actual volume change. `meltl` is a structural zero. | `test_ocean_ice_driver_column` (`window_mass_closure_*`), `test_ocean_ice_diags` |
| Frazil | with `enable` (no knob) | none | Surface layer only. Virtual salt flux. | `test_ocean_frazil`, `test_ocean_ice_coupling` |
| Atmospheric forcing | v1 stub | `air_temp`, `restore_lambda`, `sw_down` (all `0.0`) | Scalar, uniform, slab-restoring. No bulk formulae, fields, files or time dependence. | `test_ocean_ice_driver_column` |
| Ice to ocean coupling | with `enable` | none | Salt closed (virtual). Heat closed, including transmitted shortwave. Momentum **not conserved** at fractional cover. Freshwater and mass **not coupled**. | `test_ocean_ice_coupling`, `test_ocean_ice_driver_column` |
| Ice initial condition | **default** `"zero"` | `&ocean_ice_ic_nml conc_config` (`"zero"`, `"uniform"`, `"latitudes"`) | Analytic seeding of a live pack before step one. File-backed ICs are deferred. | `test_ocean_ice_init` |

**Full notes in the archive:** [sea-ice table (original cells) and notes (SIS2 divergences, deferred items)](CLOSURE_MATRIX_ARCHIVE.md#a22)

## 13. Ice-shelf cavity

Three pieces make up the ice-shelf cavity capability. All are ocean-only and
default off. The top drag is in [section 9](#9-drag-and-surface-stress).

1. **Geometry and load** (`&ocean_cavity_dyn_nml`): a static ice draft, the
   barotropic datum that absorbs it, and the load on top of the water column.
2. **Basal melt** (`&ocean_cavity_melt_nml`): ice-ocean interface thermodynamics,
   delivered as a virtual salt flux by default.
3. **Top drag** (`&ocean_tdrag_nml`): momentum sink under the ice.

### Geometry and load

| Piece | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Static ice draft `z_draft` | opt-in | `&ocean_cavity_dyn_nml enable` (`.false.`), `draft_config` (`"none"`; `"flat"`, `"linear"`; `"file"` fails loud), `draft_source` (`"draft"`; `"thickness"`; `"in_situ"` fails loud) | Time-constant draft laid over the bed. | `test_ocean_cavity_draft` |
| Barotropic datum `bt_H_ref = b − z_draft` | with `enable` | none | `bt_eta` becomes the deviation from the loaded equilibrium. | `test_ocean_cavity_equivalence` |
| Grounding becomes land | with `enable` | `h_min_cavity` (`10` m), `grounded_max_frac` (`0.5`) | Columns thinner than `h_min_cavity` are land. | — |
| Isostatic load `p_ice_ref` | with `enable` | `trim_ic_for_p_surf` (`.false.`) | Flotation load `ρ_ref·g·z_draft`. | `test_ocean_cavity_load` |
| Load partition `ms%p_top` | with `enable` | `&ocean_pgf_nml p_top_in_bc`, `&ocean_psurf_nml in_eos` | `p_top = p_ice_ref + sf%p_surf`. | `test_ocean_cavity_load`, `test_ocean_cavity_equivalence` |
| `p_top_in_bc` for a varying draft | **required** (refused otherwise) | `&ocean_pgf_nml p_top_in_bc` | A uniform draft is exempt. | — |
| Rest residual under a sloping draft | measured, not a knob | none | Sigma-PGF truncation, cubic in the slope. | `test_ocean_cavity_load` |

**Envelope** (each fails loud at configure): `form = "fv_mom6"` with
`gfs_scale = 1`; `vcoord_type` in {`sigma`, `z_fixed`}; the split solver; one
rank; `bt_halo = 0`; and no wet/dry, porous barriers, sea ice or
`&ocean_tides_nml use_sal`. `&ocean_zinit_nml enable` composes.

**Full notes in the archive:** [geometry and load table (original cells, measurements) and envelope](CLOSURE_MATRIX_ARCHIVE.md#a23)

### Basal melt

Needs `&ocean_cavity_dyn_nml`, `&ocean_eos_nml tfreeze_set="isomip"` and
`&ocean_forcing_nml enable_components`. Each is a fail-loud configure
requirement. Kernel: `src/parameterizations/vertical/rdb_ocean_cavity_melt.F90`.
Coupling: `rdb_ocean_cavity_flux.F90`, called from `engine_step_finalize`.

| Capability | Status | Knob (default) | What it does | Tests |
|---|---|---|---|---|
| Three-equation interface solve | opt-in | `&ocean_cavity_melt_nml enable` (`.false.`) | Liquidus, heat and salt balance at the ice base. | `test_ocean_cavity_melt` |
| Two-equation variant | with `enable` | none | The `γ_S → ∞` limit. Matches the sea-ice basal flux law. | `test_ocean_cavity_melt` |
| Exchange law `const_gamma` | **default** | `exchange_law = "const_gamma"`, `gamma_t` (`2.2e-2`) | `γ = Γ·u*`. | `test_ocean_cavity_melt` |
| Exchange laws `hj99`, `yung25` | available | `exchange_law` | Implicit laws solved by bisection. | `test_ocean_cavity_melt` |
| Exchange laws `jenkins91`, `rosevear22`, `vt19`, `mk18`, `burchard22`, `jenkins21` | **reserved** (refused) | `exchange_law` | Never a silent fallback. | `test_ocean_cavity_melt` |
| Ice conduction | **default** `"insulating"` | `ice_conduction` (`"insulating"` or the H&J99 advective-diffusive form) | The purely diffusive option is reserved. | `test_ocean_cavity_melt` |
| Friction velocity `u*` | with `enable` | `cdrag_top` (`2.5e-3`), `ustar_min` (`1.0e-4`) | `u* = max(√(C_d(u²+v²+u_tide²)), u*_min)`. | `test_ocean_cavity_melt` |
| Far-field sampling | **default** `10` m | `far_field_depth` | Mean over metres below the ice base, never "layer nz". | `test_ocean_cavity_flux` |
| Interface pressure `ms%p_top` | with `enable` | none | Third consumer of `p_top`. | `test_ocean_cavity_flux` |
| Flux delivery | with `enable` | none | `heat_cavity` and `salt_cavity` fold into `Q_heat` and `Q_salt`. | `test_ocean_cavity_flux` |
| Virtual salt flux | **default** | `freshwater = "virtual"` | Dilution identity. Adds no volume. | `test_ocean_cavity_flux` |
| Real freshwater mass | opt-in | `freshwater = "mass"` | Adds `dh = m·dt/ρ₀` to the top layer. Refused with wet/dry and `dt_tracer_advect_ratio > 1`. | `test_ocean_cavity_freshwater` |
| Sea-level compensation | opt-in; needs `freshwater="mass"` | `volume_compensation` (`"none"`; `"uniform_open_ocean"`) | Removes the melt volume over uncovered cells. | `test_ocean_cavity_freshwater` |
| Budgets | with `enable` | none | Rides the existing heat and salt surface budgets; mass adds `ms%mass_src`. | `test_ocean_cavity_flux`, `test_ocean_cavity_freshwater`, `test_ocean_cavity_grounded_budget` |
| Status accounting | with `enable` | none | Non-finite input or state is fatal; other failures are counted and get a zero-melt safe state. | `test_ocean_cavity_flux` |
| Cover mask on atmospheric forcing | follows `&ocean_cavity_dyn_nml enable` | none | Atmospheric terms are scaled by `1 − cover_frac`. | `test_ocean_cavity_flux` |
| Diagnostics (13 entries) | with `enable` | `&ocean_cavity_melt_nml enable` (11), `&ocean_cavity_dyn_nml enable` (2) | NaN outside the cover. Fail-loud without the prerequisite knob. | `test_ocean_cavity_diags` |

**Full notes in the archive:** [basal-melt table (original cells, formulas, design rationale) and notes](CLOSURE_MATRIX_ARCHIVE.md#a24)

## 14. Boundaries

This page has no boundary-condition table. Open, wall, tidal, clamped, Chapman,
sponge and periodic edges are chosen per edge in `&ocean_bc_nml`, with Orlanski
radiation and asymmetric nudging. The types live in `rdb_ocean_boundary_types` and
the kernels in `src/core/ocean/boundary/`. Interior land masking and dynamic
wet/dry (`&ocean_wetdry_nml`) are separate. Surface and bottom *tracer* fluxes are
in [section 11](#11-surface-forcing-initial-conditions-and-other-fluxes).

## 15. Diagnostics

Diagnostics have no feedback on the dynamics.

| Diagnostic | Status | Knob (default) | What it does | Test |
|---|---|---|---|---|
| First-baroclinic wave speed and Rossby radius | opt-in | `&ocean_wavespeed_nml enable` (`.false.`) | Per-column rigid-lid eigensolve, once per outer step. Feeds the resolution scaling of GM, Redi and MEKE. Opt-in `use_ebt` (`.false.`) also computes the equivalent-barotropic speed `cg1_ebt` and the surface-normalised pressure-mode structure `ebt_struct` (surface Neumann, bed Dirichlet; PLM-remapped from the merged stack); GM consumes it with `&ocean_gm_nml use_ebt_struct`, and `mono_n2` stays unimplemented (refused with `use_ebt`). | `test_ocean_wave_speed` |
| Cavity diagnostics (13) | opt-in | see [section 13](#13-ice-shelf-cavity) | Melt rate, driving terms, exchange velocities, geometry. | `test_ocean_cavity_diags` |

**Full notes in the archive:** [wave-speed diagnostic (original cells)](CLOSURE_MATRIX_ARCHIVE.md#a25)

## 16. Tunable parameters

The curated set of key dials per scheme. Defaults are the namelist defaults
(`0` or `.false.` means off or inert). The exhaustive per-field list (units, CF
metadata) is in the `config_t` docstrings (FORD) and [`REFERENCE.md`](REFERENCE.md).
Groups covered: porous barriers, implicit stress and drag folds, wave drag, PP81
and KPP constants, tracer Laplacian, EPBL, kappa-shear, the vmix assembly gate,
geothermal heat and surface restoring.

**Porous barriers** — `&ocean_porous_nml` (ocean only)

| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch. Off ⇒ the open-area fields stay at their `(1,1,1)` placeholder and no porous kernel is launched (byte-identical). Fails loud with `&ocean_bt_nml bt_halo > 0` and with `&ocean_wetdry_nml enable`. |
| `source` | `"resolved"` | Along-face statistics source. `"resolved"` = wet-gated corner/midpoint samples of the RESOLVED bathymetry (a proxy). `"file"` = offline subgrid file — fails loud, deferred. |
| `eta_interp` | `"max"` | Interface height at the velocity point (MOM6 `PORBAR_ETA_INTERP`): `max` (the higher/shallower interface — the LEAST blocking, since `w` increases with height), `min` (the most blocking), `arithmetic`, `harmonic`. |
| `masking_depth` | `0.0` m | Faces whose mean along-face depth is SHALLOWER than this stay fully open (MOM6 `PORBAR_MASKING_DEPTH`, positive below the surface). |

**Implicit stress/drag fold** — `&ocean_vdiff_nml` (ocean path only; folds the
surface wind stress + bottom drag into the backward-Euler vertical-friction
tridiagonal as BCs instead of explicit pre-solve adds — thin-layer
(z*/ZSTAR_FULL pinch-out) CFL-robust). The folds (`implicit_stress`,
`implicit_drag`, `implicit_top_drag`) default off ⇒ bit-identical; `hvel_mom6` and
`bbl_glue` default ON (MOM6's face treatment).

| Knob | Default | Meaning |
|---|---|---|
| `implicit_stress` | `.false.` | Wind stress → vdiff surface (`k=nz`) RHS row (Neumann top-BC). Incompatible with `&ocean_vmix_nml direct_stress`. |
| `implicit_drag` | `.false.` | Bottom drag → vdiff bed-row diagonal (stress bottom-BC; the row is `k_bot_u/v`, `k=1` off `z_fixed`, and the rows below it are the identity). Mutually exclusive with `&ocean_bdrag_nml implicit`, and with HBBL (`hbbl>0` — the fold is one 2-D rate, a band needs a per-layer one) unless `bbl_glue` replaces it. |
| `hvel_mom6` | `.true.` | (DEFAULT ON since 2026-10.) MOM6 momentum face thickness `hvel` + arithmetic `h_shear` (`vertvisc_coef`); `.false.` = the historical arithmetic `h_u`. |
| `hvel_harmonic` | `.false.` | MOM6 `HARMONIC_VISC`: `.false.` = MOM6 default branch (arithmetic, harmonic near the bed for thin→thick flow, height above bed `max(zh, z_clear)` so face layers below the shallower bed of a step sit in the BBL); `.true.` = harmonic branch with the `hvel_upwind` blend. |
| `bbl_glue` | `.true.` | (DEFAULT ON since 2026-10; inert when no drag is configured.) MOM6 `BOTTOMDRAGLAW` — see the Bottom drag table. Requires `hvel_mom6`. |
| `hbbl_visc` | `10.0` m | botfn scale without the glue; the glue's HBBL when `&ocean_bdrag_nml hbbl = 0`. |
| `implicit_top_drag` | `.false.` | Ice-shelf top drag → vdiff surface (`k=nz`) diagonal, AND the wind-stress RHS on that row scaled by `(1 − cover)` so a covered face takes the drag and no wind. Requires `&ocean_tdrag_nml enable`; mutually exclusive with `&ocean_tdrag_nml implicit` and with `htbl>0`. Test: `test_ocean_top_drag`. |

**Barotropic linear wave drag** — `&ocean_bt_nml` (ocean path only; bulk
energy sink for the barotropic tide, Egbert & Ray 2001 / Jayne & St Laurent
2001). `wave_drag=.false.` (default) ⇒ bit-identical. `form="uniform"` needs
`wave_drag_r_uniform>0` (a warning fires if not). `form="roughness_proxy"` is
a documented PLACEHOLDER — resolved-bathymetry `⟨grad b⟩²` variance, not the
real subgrid `⟨h²⟩` (pending PR-14's file reader / PR-30). `form="file"` is
registered but fails loud at configure (not implemented — PR-14).

| Knob | Default | Meaning |
|---|---|---|
| `wave_drag` | `.false.` | Master switch. |
| `wave_drag_form` | `"uniform"` | `"uniform"` \| `"roughness_proxy"` \| `"file"` (fail-loud, unimplemented). |
| `wave_drag_scale` | `1.0` | Global tuning multiplier on `r_H`. |
| `wave_drag_r_uniform` | `0.0` | Piston velocity `r_H` (m/s) for `form="uniform"`. |
| `wave_drag_kappa` | `6.2832e-4` | Topographic wavenumber (1/m) for `form="roughness_proxy"` (matches `&ocean_tidal_mixing_nml kappa_itides`). |
| `wave_drag_n_bot` | `1.0e-3` | Reference bottom `N` (1/s) for `form="roughness_proxy"`. |
| `wave_drag_h2_max` | `2.5e4` | Ceiling on the `⟨h²⟩` proxy (m²) for `form="roughness_proxy"`. |

**PP81 + KPP constants** — `&ocean_vmix_nml` (routes to
`ocean_state%vmix`). `kv_bg`/`kt_bg`/`ks_bg` (the assembly floor,
`vmix_assemble`) are re-derived from `pp81_nu_bg`/`pp81_kappa_bg` at
configure time (`vmix_seed_backgrounds`) so the floor tracks a user-set
background rather than the type default.

| Knob | Default | Meaning |
|---|---|---|
| `interior_closure` | `"pp81"` | `"pp81"` or `"none"` (no shear-driven interior closure; `"large94"`/`"cvmix"` are reserved and fail loud). `"none"` needs `use_closure`. |
| `bkgnd_prandtl` | `-1` | Background Prandtl number; `> 0` makes the momentum floor `bkgnd_prandtl·kt_background + kv_bg` (`kv_bg` is `pp81_nu_bg`). Under `bkgnd_profile` the floor is `bkgnd_prandtl·kd_bg` and an unset value means 1. Only `-1` (exact) and `> 0` are legal; `(-1, 0]` is refused at configure (it used to mean a zero momentum floor under Bryan-Lewis). |
| `pp81_nu0` | `1e-2` | PP81 Richardson-dependent viscosity scale (m²/s). |
| `pp81_nu_bg` | `1e-4` | PP81 background viscosity (m²/s); also seeds `vmix%kv_bg`. |
| `pp81_kappa_bg` | `1e-5` | PP81 background diffusivity (m²/s); also seeds `vmix%kt_bg`/`ks_bg`. |
| `pp81_alpha` | `5.0` | Ri scaling in `ν = ν_bg + ν₀/(1+α·Ri)²` (paper value 5; implementations vary 4-10). |
| `shear2_floor` | `1e-10` | Floor on `\|∂u/∂z\|²+\|∂v/∂z\|²` in the PP81 Ri denominator (1/s²). |
| `kpp_ri_crit` | `0.3` | Critical bulk Richardson number for the KPP BL-depth sweep (LMD94 §3). |
| `kpp_cs_nonlocal` | `6.3` | Non-local (counter-gradient) transport coefficient `C_s` (LMD94 eq 20). |
| `kpp_c_vt2` | `1.8` | Unresolved-shear `V_t²` coefficient (LMD94 eq 23); `0` disables `V_t²`. |

**Along-coordinate tracer Laplacian (ocean path)** — `&ocean_hdiff_nml`
(`rdb_ocean_hdiff_tracer`). Diffusion along the model coordinate, NOT along
neutral surfaces (`&ocean_redi_nml` is the separate neutral/isopycnal path).
Default `kappa_h = 0.0` keeps the kernel's short-circuit intact ⇒
bit-identical. Configure-time guard aborts if
`kappa_h·dt_therm·(1/dx²+1/dy²) > 0.5` (Cartesian-metres check; skipped on
spherical/tripolar grids and when `dt_fixed` is adaptive).

| Knob | Default | Meaning |
|---|---|---|
| `kappa_h` | `0.0` | Horizontal tracer diffusivity (m²/s); `0` = no-op. |

**EPBL** — `&ocean_epbl_nml` (ocean path only; replaces the KPP overlay, PP81
interior + background continue underneath).  Full knob table + the MOM6 name
mapping: [`docs/generated_nml_knobs.md`](generated_nml_knobs.md); the
load-bearing dials:

| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (requires `use_closure` + thermodynamics; logs + disables `use_kpp`). |
| `mstar_scheme` | `"om4"` | Mechanical-TKE efficiency: `constant` / `om4` / `rh18`. |
| `mstar` | `1.2` | Constant-scheme mstar. |
| `nstar` | `0.2` | Convective-PE → TKE efficiency. |
| `tke_decay` | `2.5` | Ekman-depth / TKE-decay-scale ratio. |
| `mld_iteration` | `.true.` | Self-consistent MLD root-find (false position, `mld_tol = 1 m`). |
| `translay_scale` | `0.1` | Transition-layer mixing-length floor. |
| `prandtl` | `1.0` | `Kv = prandtl · Kd` into the momentum solve. |
| `combine` | `"add"` | Fold into the interior kv/kt additively (`"add"`), by `max` against the whole interior value (`"max"`), or `"max_shear"`: kappa-shear is merged first and ePBL adds only what exceeds it, `bg + max(ePBL, kappa-shear) + others`, kv likewise with `prandtl_turb`/`prandtl` (ACCESS-OM3, `EPBL_IS_ADDITIVE = False`). `"max_shear"` equals `"add"` when kappa-shear is off, is refused with `use_kpp`, and in the kappa-shear vertex form the cell-centred shear viscosity is zero so kv combines additively. Test `test_ocean_epbl_max_shear`. |
| `tke_diags` | `.false.` | Per-column TKE budget terms (close to round-off; asserted by `test_ocean_epbl`). |
| `use_lt` | `.false.` | Langmuir enhancement of mstar — LF17 wind-only statistical waves (one La per column from u* + MLD; no wave model). `lt_scheme = "rescale"\|"additive"`, Reichl & Li (2019) coefficients, `lt_lac1..5` stability modification. |

**kappa-shear** — `&ocean_kappa_shear_nml` (ocean path only; interior
closure, coexists with KPP/EPBL — its κ is ADDED to `kt` and
`prandtl_turb·κ` to `kv` every stage; the column solve runs at thermo
cadence).  Full knob table + the upstream-name mapping:
[`docs/generated_nml_knobs.md`](generated_nml_knobs.md); the
load-bearing dials:

| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (requires `use_closure` + thermodynamics; no exclusions). |
| `ri_crit` | `0.25` | Critical Richardson number for the shear source. |
| `shearmix_rate` | `0.089` | Source-rate coefficient (JHL08 calibration). |
| `kappa_0` | `1.0e-7` | Background κ floor / well-posedness smoothing (m²/s). |
| `tol_err` | `0.1` | Picard convergence tolerance (also scales the adaptive-dt bands). |
| `max_inner_it` / `max_substep_it` | `50` / `13` | Inner Picard / outer substep caps. |
| `prandtl_turb` | `1.0` | `Kv = prandtl_turb · Kd` on the momentum side of the merge. |
| `massless_merge` | `.false.` | D4: fold vanished (`< H_VANISHED`) layers onto the column's massive sub-grid before the solve (vs the blunt gather floor); identity columns bypass, so on-path is bit-identical on healthy envelopes. |
| `at_vertex` | `.false.` | MOM6 `VERTEX_SHEAR` (OM5-class production form): solve the columns at C-grid CORNERS from the native face velocities (no centre-average shear damping), then average corner Kd back to tracer points.  Kv is routed corner→face (`prandtl_turb·kd_corner` into `vdiff_apply_momentum`'s `kv_corner_source`, MOM6 `Kv_shear_Bu`) with the cell-centred kv merge suppressed; `tke_int` zeroed in vertex mode. |
| `vertex_geometric_mean` | `.false.` | Geometric (vs arithmetic) corner→centre Kd mean (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN`). |
| `vertex_geomean_kdmin` | `0.0` | Floor (m²/s) on each corner Kd BEFORE the geometric mean (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN_KDMIN`; inert unless geometric).  With 0 the geometric mean hard-zeros Kd at every shear-zone edge; OM5 configs use `1e-9`. |

**Vmix assembly gate** — `&ocean_vmix_nml` (ocean path only; the single
downstream floor/clip/smooth/guard stage `vmix_assemble`, run once per stage
after every closure has contributed into kv/kt and `vmix_split_kd_heat_salt`
has derived ks from kt).  All defaults are bit-identical to the pre-assembly
chain.

| Knob | Default | Meaning |
|---|---|---|
| `kv_max` | `huge` (off) | Ceiling on momentum viscosity kv (MOM6 `Kd_max` momentum). |
| `kd_max` | `huge` (off) | Ceiling on tracer diffusivity kt/ks (MOM6 `Kd_max`). |
| `kd_smooth_iterations` | `0` | 1-2-1 horizontal smoothing passes on kv/kt/ks (MOM6 `Kd_smooth`). |
| `vmix_guard` | `.false.` | Debug-gated negative/NaN diffusivity guard (error-stop, or status-returning testable path). |

**Geothermal bottom heat flux** — `&ocean_geothermal_nml` (ocean path only;
bed-side analogue of the surface heat flux, deposited into the lowest massive
layer)

| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (default off ⇒ bit-identical to no geothermal). |
| `q_geo` | `0.0` | Constant bottom heat flux (W/m², positive into the ocean; typical ~0.05–0.1). |

**Surface buoyancy restoring** — `&ocean_restore_nml` (ocean path only;
relaxes the surface top layer T/S toward scalar targets at rate
`λ = piston/h_top`; piston-velocity form, MOM6 `RESTOREBUOY`-inspired;
non-conservative — increment mirrored into the surface heat/salt budget)

| Knob | Default | Meaning |
|---|---|---|
| `enable_restore_temp` | `.false.` | Master switch for SST restoring (effective only when `piston_t /= 0`). |
| `enable_restore_salt` | `.false.` | Master switch for SSS restoring (effective only when `piston_s /= 0`). |
| `piston_t` | `0.0` | SST piston velocity (m/day, MOM6 `FLUXCONST_T`; converted to m/s at seed). |
| `piston_s` | `0.0` | SSS piston velocity (m/day, MOM6 `FLUXCONST_S`). |
| `restore_sst` | `0.0` | Scalar target SST (degC), in the model T convention (θ, or CT under `cons_abs`); never converted. |
| `restore_sss` | `0.0` | Scalar target SSS (PSU), in the model S convention (SP, or SA under `cons_abs`); never converted. |
| `sss_file` | `""` | File-based SSS target (OM3 `SALT_RESTORE_FILE`); blank keeps the scalar target. See the file-based restoring row of the forcing table. |
| `sss_var` | `"salt"` | Variable in `sss_file`. |
| `sss_time_cyclic` / `sss_cycle_period` | `.true.` / `31536000` s | Climatology wrap (reader `DATA_TIME_CYCLIC`) and its period. |
| `max_delta_sss` | `999` | Cap on `|S_target - S|` (g/kg); 999 = no clamp. |
| `sss_adjust_net_zero` / `sss_adjust_by_scaling` | `.false.` / `.false.` | Remove the global net salt flux; by scaling, else by subtracting the wet-area mean. |
| `sss_as_sflux` | `.false.` | Must be `.true.` with `sss_file` (salt-flux form is the only one). |
| `sss_convention` | `"model"` | `"practical"`: the file is SP, converted to SA under `cons_abs`. |

**Ocean** dials live in the per-concern `&ocean_<group>_nml` sub-namelists
(`vmix`, `hvisc`, `bdrag`, `pgf`, `coriolis`, …) — see the working envelope in
[`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).

<!-- TODO: windowed tracer advect (dt_tracer_advect_ratio) — add a row to a
     tracer-advection cadence section when the tracer-advection matrix branch
     (docs/tracer-advection-matrix) merges; the scheme lives in
     continuity_t (rdb_continuity.F90) and is gated by &ocean_vmix_nml. -->

## 17. Keeping this current

1. **Per-capability PR step.** The CLAUDE.md "adding a new capability" workflow
   gains a step: *namelist knob → kernel → unit test → **closure-matrix row** →
   pre-commit + fortitude → ctest → commit.* One row per scheme is cheap.
2. **Automated gate.** `tools/check_closure_matrix.py` parses the closure knobs
   from `src/core/rdb_config.F90` namelists and the ocean sub-namelists, asserts
   each appears somewhere in this file, and checks that every test named here
   exists in `tests/CMakeLists.txt`. Wired into pre-commit + CI so drift fails
   the build instead of rotting silently (this matrix was born because the prose
   in `CAPABILITIES_AND_LIMITATIONS.md` drifted — it still claims PP81 is the
   only coastal mixing scheme, three closures ago).

## Appendix A: change history

Dated reconciliation log, kept verbatim. Read it for the history of why a row
looks the way it does. It is not needed to use the tables.

- **2026-10-27** — `&ocean_epbl_nml combine = "max_shear"` (third value, no
  new knob): kappa-shear runs and merges BEFORE the ePBL merge on this branch
  only and ePBL adds `max(kd_epbl - kd_shear, 0)` (kv: `max(Pr_epbl kd_epbl -
  prandtl_turb kd_shear, 0)`), giving `bg + max(ePBL, kappa-shear)`. The
  `"add"`/`"max"` order is unchanged (bit-identical).
- **2026-10-27** — `&ocean_vmix_nml interior_closure` (`"pp81"` default |
  `"none"`) and `bkgnd_prandtl` (default `-1`, off). `"none"` runs no PP81
  kernel: the interior `kv`/`kt` are seeded each stage with the background
  (scalar, Henyey or Bryan-Lewis; additive with the contributors), with
  ePBL, kappa-shear, tides, convection and double diffusion adding on top
  (ACCESS-OM3 sets `USE_PP81 = False`). `bkgnd_prandtl > 0` ties the momentum
  floor to the tracer background in use. `bkgnd_prandtl` was previously a
  Bryan-Lewis-only knob defaulting to 1; the default is now the off sentinel
  and an unset value resolves to 1 under `bkgnd_profile`, so every existing
  configuration is bit-identical. KPP is not refused with `"none"`: its
  overlay only `max()`es into whatever is in `kv`/`kt`.
- **2026-10-27** — OM3 lateral-viscosity floors: `&ocean_hvisc_nml`
  `kh_sin_lat` / `kh_pwr_of_sine` (latitude Laplacian floor), `ah_vel_scale`
  (biharmonic velocity-scale floor) and `bound_coriolis_biharm` /
  `bound_coriolis_vel` (Coriolis bound on Smagorinsky_AH). All default off,
  bit-identical; the biharmonic knobs require `smag_ah`, `kh_sin_lat`
  requires a non-cartesian grid.
- **2026-10-27** — MEKE GM source from the bolus work: `&ocean_gm_nml
  gm_src_form` (`"local"` default, bit-identical | `"work"`; inverse of MOM6
  `MEKE_GM_SRC_ALT`, OM3 uses `"work"`). `gm_work_source` replaces `gm_src` by
  the PE-release rate of the final transports (sees FGNV, EBT structure and
  the limiters, which the local form does not). Test: `test_ocean_gm_work`.
- **2026-10-27** — GM PR-4: `&ocean_gm_nml use_fgnv`, `fgnv_filter_scale`,
  `fgnv_c_min`, `fgnv_strat_floor` (MOM6 `KHTH_USE_FGNV_STREAMFUNCTION`,
  `FGNV_FILTER_SCALE`, `FGNV_C_MIN`, `FGNV_STRAT_FLOOR`). Ferrari et al.
  (2010) streamfunction BVP per face column (symmetric tridiagonal, Thomas)
  on the (EBT-structured when `use_ebt_struct`) unlimited streamfunction,
  default off, bit-identical. Test: `test_ocean_gm_fgnv`.
- **2026-10-27** — GM PR-3: `&ocean_gm_nml use_ebt_struct` (default off,
  bit-identical) gives the thickness diffusivity the EBT vertical structure
  of the wave-speed slot; the "Deferred: EBT" item in the GM history is
  closed.  Refused without GM or without `use_ebt`.
- **2026-10-27** — File-based SSS restoring (OM3 `RESTORE_SALINITY`,
  `SALT_RESTORE_FILE`, `MAX_DELTA_SRESTORE`, `ADJUST_NET_SRESTORE_*`,
  `SRESTORE_AS_SFLUX`): `&ocean_restore_nml sss_file`, `sss_var`,
  `sss_time_cyclic`, `sss_cycle_period`, `max_delta_sss`,
  `sss_adjust_net_zero`, `sss_adjust_by_scaling`, `sss_as_sflux`,
  `sss_convention`. Salt-flux form with a per-cell clamp and a global
  zero-net adjustment (EFP sums); defaults off, bit-identical.
- **2026-10-27** — Redi diffusive-CFL control: `&ocean_redi_nml
  check_diffusive_cfl` and `max_diffusive_cfl` (OM3 `CHECK_DIFFUSIVE_CFL` /
  `MAX_TR_DIFFUSION_CFL`). Local face cap, global-max subcycle count,
  inter-pass tracer halo refresh. Defaults off, bit-identical. Configure now
  warns on `kappa_h > 0` with Redi on.
- **2026-10-26** — TEOS-10 PR 3: `&ocean_zinit_nml input_convention`
  (`pot_prac` default | `cons_abs`); file z-init, bulk flux and scalar
  restoring now run under `cons_abs` (their PR 2b refusals are lifted);
  derived diagnostics `temperature_pt` / `salinity_sp`. Defaults
  bit-identical. Still refused under `cons_abs`: sea ice (PR 4) and cavity
  melt.
- **2026-10-25** — TEOS-10 PR 4: `cons_abs` + `&ocean_ice_nml enable` is no
  longer refused. The ice-ocean boundary converts once (IOC/SCOR/IAPSO 2010;
  McDougall & Barker 2011): the ice stays θ/SP, the ocean heat budget stays
  CT/SA. Frazil and basal-flux heat use `eos%c_p` instead of the fixed 3992
  (same default, bit-identical). Frazil is evaluated at p = 0 only. No new
  knobs; `pot_prac` is unchanged.
- **2026-10-25** — `ts_convention = "cons_abs"`, `tfreeze_form = "teos_poly"`,
  `tfreeze_sat_frac` and `c_p` added (§6). OM3's `MOM_parameter_doc.all`
  (`USE_CONTEMP_ABSSAL = True`, `EQN_OF_STATE = "ROQUET_RHO"`,
  `TFREEZE_FORM = "TEOS_POLY"`, `C_P = 3992`) fixed the defaults: sat
  fraction 0 and `c_p = 3992` under both conventions. All default to the previous behaviour
  (bit-identical). Not yet covered under `cons_abs`: initial-condition
  conversion, PT/SP diagnostics, sea ice and cavity melt (refused);
  initial conditions, PT/SP diagnostics, bulk flux and restoring landed in PR 3.
- **2026-10-13** — `eos = "roquet_rho"` added (§6). **Answer change for
  `roquet_spv` and `roquet_rho`:** their PT→CT conversion carried two
  misplaced polynomial terms; it now uses the GSW-verified
  `rdb_teos10_convert`. Over SP 0–42, θ −2 to 32 °C, p 0/1000/5000 dbar:
  max |ΔCT| = 2.3e-2 °C, max |Δρ| = 4.5e-3 kg/m³ (both variants). Pinned by
  `test_ocean_eos_roquet :: roquet_model_ct_pins_gsw`.
- **2026-10-13** — `eos = "roquet_rho"` through the FV-MOM6 in-situ density
  PGF: `roquet_rho` now has its own factored in-situ PCM + reconstruction
  kernels (`roquet_rho_pcm_dpa_intz`/`_face`, `roquet_rho_recon_dpa_intz`/
  `_face` in `rdb_ocean_pressure_force.F90`, mirroring the Wright and
  Roquet-SpV paths; Adcroft, Hallberg & Harrison 2008 analytic-FV form),
  wired into the 3-way dispatch of `compute_fv_mom6_insitu_pcm_impl` and
  `compute_fv_mom6_reconstruct_impl`. The configure-time refusal of
  `eos="roquet_rho"` + `&ocean_pgf_nml form="fv_mom6"` with the default
  `insitu_density=.true.` is **lifted** — the combination now validates and
  runs the in-situ path by default, matching `roquet_spv`. Knob surface is
  unchanged (no new namelist key), and `roquet_rho` was refused before this
  change, so no existing accepted run re-baselines. Gated by
  `test_ocean_pgf_roquet_rho`, `test_ocean_pgf_insitu ::
  thermobaric_face_matches_insitu_oracle_roquet_rho` /
  `rest_partial_steps_zero_pgf_roquet_rho`, `test_config ::
  roquet_rho_insitu_pgf_accepted`.

**Full notes in the archive:** [the full dated change log (2026-08-04 to 2026-10-06)](CLOSURE_MATRIX_ARCHIVE.md#a26)
