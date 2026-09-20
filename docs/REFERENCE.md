# Roundabout — Physics & Numerics Reference

Canonical reference for the physics, numerical schemes, boundary
conditions, validation tests, and namelist configuration. **Codebase
structure, build system, GPU programming model, and kernel API live in
[`CLAUDE.md`](../CLAUDE.md) + [`docs/codebase/`](codebase/INDEX.md) —
this file is the science-and-config view.** Ocean dynamical-core design
contract: [`src/core/ocean/README.md`](../src/core/ocean/README.md) and
[`docs/ROADMAP_OCEAN.md`](ROADMAP_OCEAN.md). Per-closure knob tables:
[`docs/CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md) and the auto-generated
[`docs/generated_nml_knobs.md`](generated_nml_knobs.md).

> **Scope note.** The `sim_type='coastal'` regime — the A-grid HLL/HLLC
> path, the unstructured triangular (KNP central-upwind) backend, the
> semi-implicit (Casulli) free surface, the non-hydrostatic extension,
> and the C/Python FFI — was split out into a **separate repository**.
> This tree is ocean-only. Nothing below describes those paths.

*Last reconciled against the source: 2026-08-23 (coastal carve-out).*

## 1. Overview

Portable GPU-native ocean solver written in Fortran. One regime, one
backend:

- **`sim_type='ocean'`** — structured Arakawa C-grid hydrostatic
  dynamical core. Continuity-PPM layer transport, PV-conserving
  (Sadourny) Coriolis-advection, finite-volume pressure gradient,
  split-explicit RK2, Wright (1997) nonlinear EOS, ALE
  Lagrangian-then-remap vertical coordinate. Targets 1 km
  submesoscale-resolving regional-to-basin hydrostatic ocean.

`sim_type` is pinned: `"ocean"` is the default, the only value in the
config schema's enum, and `validate_config` aborts fail-loud on anything
else.

Horizontal grids are selected by `&ocean_grid_nml grid_config`:
`cartesian` (default), `spherical` lon-lat sector, `supergrid` (MOM6
mosaic reader), and `tripolar` (Murray 1996 bipolar Arctic cap closed by
a single-rank north fold).

## 2. Governing Equations

Layered hydrostatic Boussinesq primitive equations in a generalised
vertical coordinate, vector-invariant form (per layer `k`, bottom-up):

```
∂h_k/∂t + ∇·(h_k u_k) = 0

∂u_k/∂t + (f + ζ_k) ẑ × u_k + ∇(KE_k) =
    -(1/ρ_0) ∇p_k + ∇·(ν_h ∇u_k) + ∂/∂z(ν_v ∂u_k/∂z) + F_k

∂(h_k θ_k)/∂t + ∇·(h_k u_k θ_k) =
    ∇·(κ_h h_k ∇θ_k) + ∂/∂z(κ_v ∂θ_k/∂z) + Q_θ

p_k = p_surf + g Σ_{m>k} ρ_m h_m + ½ g ρ_k h_k      (hydrostatic)
ρ_k = ρ(S_k, T_k, p)                                (equation of state)
```

- `h_k` layer thickness — the prognostic continuity variable
- `u_k` layer horizontal velocity, stored on C-grid faces
- `η` free-surface elevation (`Σ_k h_k = H + η`)
- `ζ_k` relative vorticity at cell corners, `KE_k` kinetic energy at
  cell centres
- `f` Coriolis parameter; `θ` any registered tracer (salinity,
  temperature, passive)
- `F_k` surface wind stress / bottom drag / lateral-closure forcing

Continuity is a **transport equation**, not a constraint: the layer
thickness is prognostic, so there is no Poisson solve and no FFT
projection anywhere in the dynamical core.

## 3. Numerical Schemes

### Spatial discretisation

- **Grid**: Arakawa C-grid finite volume — `η`, `h`, `S`, `T`, `ρ` at
  cell centres; `u` on x-faces, `v` on y-faces; `ζ` and PV at corners.
  SoA `(nx, ny, nz)` for coalesced GPU access.
- **Continuity**: continuity-PPM (piecewise-parabolic thickness
  reconstruction with a positive-definite limiter) drives the per-layer
  mass fluxes; tracers ride the same fluxes (`&ocean_continuity_nml`).
- **Coriolis + advection**: Sadourny PV-flux vector-invariant form
  (`&ocean_coriolis_nml form=`) — enstrophy-conserving (default,
  velocity form), energy transport form (`sadourny_energy`), or the
  Hollingsworth-Källén guard (`sadourny_hk`). The orthogonal
  `pv_adv_scheme=` selects the corner-vorticity→face interpolation:
  `centered` (default, bit-identical 2-point average) or
  `weno3`/`weno5`/`weno7` (upwind-biased WENO-Z, MOM6 WENOVI{3,5,7}TH).
  `weno5`/`weno7` need `nghost ≥ 3`/`4` (fail-loud via
  `pv_adv_required_nghost`).
- **Pressure gradient**: finite-volume, `&ocean_pgf_nml form=` —
  Montgomery potential (`mont`, **default**), z-corrected FV
  (`fv_lite`), FV + Wright in-situ
  Picard re-evaluation (`fv_wright`), the MOM6 `PressureForce_FV_Bouss`
  port (`fv_mom6`), or the NK=2-only reduced-gravity `gprime`. Optional
  in-layer PLM/PPM T/S reconstruction (Boole quadrature).
- **Vertical**: ALE — advance Lagrangian, then conservatively remap onto
  the target coordinate (`&vcoord_nml remap_method`, PPM by default).

### Time integration

Split-explicit RK2. Each outer step:

1. **Slow tendencies** — Coriolis-advection, PGF, lateral viscosity,
   vertical mixing (backward-Euler tridiagonal), bottom drag, surface
   stress, tracer hdiff/vdiff, vertical advection.
2. **Barotropic fast loop** — `n_inner` forward-backward-Euler substeps
   on `(η, u_bt, v_bt)` plus the nonlinear `ζ + KE` terms.
   `&ocean_bt_nml auto_n_inner = .true.` derives `n_inner` from the
   gravity-wave CFL **once, at configure time** (`configure_ocean_bt_split`),
   writing the resolved value back into the config. It is not re-derived
   per step: CFL truncation is a counter (`dyn%ntrunc_total`), not a
   controller.
3. **BT correction** — distributes the barotropic `Δu` back into the
   layers; optionally h-weighted (`correction_h_weighted`) and further
   biased by the vdiff viscous remnant γ (`correction_visc_rem`).
4. **Stage 2** — under `ssp_rk2` a second, identical stage runs and the
   two stage outputs are averaged (Heun). Under the default `pred_corr`
   there is a single prognostic update: the predictor advances a
   provisional velocity to `pc_be * dt`, and the corrector takes the full
   step with tendencies evaluated on the barotropic-window means.

`&ocean_bt_nml split_scheme` selects the outer scheme: `pred_corr`
(default — the MOM6 predictor-corrector, which lifts the internal-wave
`dt` ceiling) or `ssp_rk2` (experimental — widest envelope, but it grows
internal gravity waves out of a stratified rest state; see
[`CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md)). The ALE remap and the thermodynamic kernels
are gated on the thermo cadence (`&ocean_vmix_nml dt_therm_ratio`, MOM6
`DT_THERM`); horizontal tracer advection has its own window
(`dt_tracer_advect_ratio`).

Because the gravity-wave mode is sub-cycled rather than resolved by the
outer step, the outer `Δt` is bounded by advective/baroclinic CFL, not by
`dx/√(gH)` — O(minutes) in practice (`dt_fixed = 1200 s` in the
MOM6-reference double gyre). The driver requires `dt_fixed > 0`; there is
no adaptive-CFL helper on this path.

### Design rationale

Explicit time-stepping with a sub-cycled fast mode maps to GPU
parallelism without global solves, outside the per-column tridiagonals.
`do concurrent` is portable across GPU vendors with a CPU fallback and no
separate kernel implementations. Continuity as a transport equation keeps
the vertical coordinate free to be Lagrangian, which is what makes the
ALE remap a drop-in rather than a rewrite.

## 4. Boundary Conditions

Per-edge tags via `&ocean_bc_nml` (`west` / `east` / `south` / `north`),
parsed onto the `OBC_*` taxonomy in `rdb_ocean_boundary_types`:

| Type | Code | Description |
|------|------|-------------|
| `OBC_WALL` | 1 | Closed wall (hard zero). Default for every edge. |
| `OBC_OPEN` | 2 | Flather radiation — gravity-wave outflow, η clamped to a reference (`data_eta_*`) |
| `OBC_TIDAL` | 3 | Prescribed multi-constituent η, composed into `data_eta_*` |
| `OBC_NESTED` | 4 | Reserved (not implemented); behaves like OPEN to the kernels |
| `OBC_INFLOW` | 5 | Reserved — prescribed normal velocity + tracer |
| `OBC_DISCHARGE` | 6 | Reserved — prescribed volume flux |
| `OBC_CLAMPED` | 7 | Hard Dirichlet on η + u + v + per-tracer values, sourced from `data_*` |
| `OBC_SPONGE` | 8 | Relaxation band — the BC kernel falls through to WALL at the outer face; the sponge kernel relaxes the interior band toward `data_*` targets |
| `OBC_CHAPMAN` | 9 | Orlanski radiation on η with implicit phase-speed estimation; persistent `eta_old_<edge>` state |
| `OBC_PERIODIC` | 10 | Ghost-wrap periodic. Requires `nghost ≥ 3` (PPM + biharmonic stencil) and must be paired (west ⟺ east; south ⟺ north). Cannot combine with `OBC_SPONGE` on the same edge. |
| `OBC_TRIPOLAR_FOLD` | 11 | Tripolar north-fold seam (Murray 1996). NORTH edge only; requires periodic west/east. |

Open edges carry a per-layer zero-gradient baroclinic anomaly on top of
the barotropic Flather mean, plus per-edge tracer inflow values and
asymmetric nudging. Boundary tides pick up the nodal correction when
`&ocean_bc_nml obc_tidal_nodal` is set — it shares the `&ocean_tides_nml`
reference epoch with the interior body tide, so the two stay
phase-consistent.

The sponge has two implementations: the legacy `bc`-band relaxation
(default) and the real per-cell sponge (`&ocean_sponge_nml enable`, with
per-edge `*_width` / `*_strength` overrides) which relaxes momentum
toward `u_ref`/`v_ref` and every registered tracer toward `ref_tracer`,
mirroring S/T into the sponge budget accumulators.

**Deferred**: per-layer Orlanski phase-speed radiation, file-backed
boundary-data backends (the polymorphic `update(t, bc)` call site is
wired in the driver; only the constant backend ships), and live nesting.

## 5. Layered Solver

The layered C-grid state (`multilayer_state_t`) carries `h_layer`,
`hu/hv_layer`, the tracer registry, `rho_layer`, and the budget
accumulators. `&nonhydrostatic_nml nz_layers` sets the layer count;
`use_multilayer` enables the coupled vertical stack.

### Vertical layer convention (load-bearing)

**Bottom-up (ROMS-style)**: `k=1` is the bed, `k=nz` is the surface
layer. Surface forcings land at `(:, :, nz)`; bed forcings at
`(:, :, 1)`. Any new kernel must follow this — PGF, `rho_ref`, per-layer
drag, vcoord output, surface tracer fluxes and vertical exchange all
assume it.

### Vertical coordinates

Dispatched via `vcoord_type` (namelist) → `parse_vcoord_type` →
`VCOORD_*` enum. Every family routes through the same ALE remap:

| Tag | Coord | Status | Notes |
|---|---|---|---|
| `sigma` | Pure sigma | Production | Terrain-following; layers rescaled by `dsig(k) * H`. Default. |
| `zsigma` | Smoothstep blend sigma→z | Production | Conservative remap. The blend depth/width are internal defaults on `ocean_vcoord_t`, not namelist knobs. |
| `zstar` | z*-lite | Production | Single global `z_ref` stretched per column. SSH-tracking. |
| `zstar_sigma` | Sigma + z*-lite hybrid | Production | Sigma shallow, z*-lite deep. Conserves by construction. |
| `zstar_full` | Full per-column z* | Production w/ caveat | Per-column `z_ref(0:nz)` from local bathymetry; surface layer anchored at `zstar_h_surf_target` regardless of H. Bed-side layers can vanish (`zstar_h_min`). Operators that divide by `h_layer` gate on `H_VANISHED = 1.5e-4 m`, and `zstar_h_min` is deliberately kept at or below that marker so the below-bed filler layers stay inert — `validate_config` warns on a larger value (see `rdb_vcoord :: vcoord_h_min_role`). **Caveat**: geometry with frequent wet/dry cycling still leaks salt across the cycle; prefer sigma or z*-lite there. |
| `eulerian_z` / `z` | Fixed z levels | Production | Static interface depths; the classic Eulerian z-coordinate. |
| `z_fixed` / `gprime` | Fixed interface stack | Production | Prescribed interface depths — the reduced-gravity setups. |
| `lagrangian` / `isopycnal` | No remap | Production | The layer *is* the coordinate; pair with `&ocean_isopycnal_nml` grounding controls. |
| `rho` | Isopycnal | Validation-grade | Places layer interfaces on prescribed potential-density surfaces `rho_target(0:nz)` (lightest→surface, densest→bed), referenced to `rho_ref_pressure` (Pa, default 2e7). PPM-reconstructs the column density profile (from T/S via the device EOS) and inverts (bracket sweep + fixed-8-iter Newton) for the depth where ρ equals each interior target; pre-compaction strips vanished source layers, MOM6 min-thickness inflation floors collapsed layers at `max(zstar_h_min, H_VANISHED)`. Conserves T·h / S·h to round-off. **Caveat**: weakly-stratified columns collapse — RHO alone is validation-grade; the `hycom` hybrid is the production follow-on. |
| `hycom` | Hybrid z*/isopycnal | Production isopycnal coord | The `rho` inversion + a z* surface-resolution floor (Bleck 2002 / MOM6 HYCOM1): after the density-space inversion, a top-down sweep enforces `z(k) ≥ Σ dsig·(H+η)` so the near-surface band keeps fixed z* resolution (no surface collapse) while the deep interior tracks isopycnals. Reuses `rho_target`/`rho_ref_pressure` + the existing `dsig` for the z* band; adds a bottom-up density monotonize before the inversion. Conserves to round-off; the z* floor is a *surface-side minimum* (deep interior may still collapse — by design). `only_improves` + interface-depth caps deferred. |

`&vcoord_nml thickness_config` selects the **initial** layer-thickness
profile, which is distinct from the running coordinate — see the note
under `&vcoord_nml` in §7.

### Tracer registry

Two prognostic tracers are registered by default:

- **Salinity** `S` (PSU), conserved as `hS = h_layer × S`
- **Temperature** `T` (°C), conserved as `hT = h_layer × T`

`multilayer_state_t%register_passive_tracer(grid, name, units,
long_name, idx)` grows the registry at setup (6-arg; `idx = 0` on
refusal; `registry_locked` after `enter_data` closes the `mem:separate`
device-map foot-gun). Every passive-transport kernel — advection, ALE
remap, vertical exchange, vertical/horizontal diffusion, Redi, sponge,
halo/fold, OBC ghost + reservoirs, restart — already loops the registry,
so a newly registered tracer rides them for free; only a named EOS /
surface-flux coupling and a diag `fill_<name>` need hand-wiring. Budget
attribution is `tracer_t%budget_id` (`NONE`/`HEAT`/`SALT`), never an
index comparison.

Shipped packages: **ideal age** (`&ocean_tracers_nml enable_ideal_age`)
and **pseudo-salt** (`enable_pseudo_salt`, Shao 2016 verification tracer
— seeded to S, given S's surface salt flux and non-local mirror; the
deviation measures the passive-vs-active transport-path error). Both
default off ⇒ bit-identical. Recipe:
[`docs/howto/add_passive_tracer.md`](howto/add_passive_tracer.md).

### Equation of state

`&ocean_eos_nml` selects the density kernel. The Wright (1997) nonlinear
rational fit is the production path; a Roquet et al. (2015) TEOS-10
polynomial and a linear form are available for validation. The linear
form is

```
ρ_k = ρ_0 + β_S × (S_k - S_ref) - α_T × (T_k - T_ref)
```

with `β_S` the haline contraction coefficient and `α_T` the thermal
expansion coefficient.

> **`α_T` and `β_S` are DIMENSIONAL — read this before writing any
> buoyancy expression.** This is the density-ANOMALY form, so
> `α_T = −∂ρ/∂T` in **kg m⁻³ per °C** and `β_S = +∂ρ/∂S` in **kg m⁻³ per
> PSU** (realistic seawater: `α_T ≈ 0.17–0.25`, `β_S ≈ 0.78`). They are
> **not** the fractional `(1/ρ)·∂ρ/∂T` coefficients (`≈ 2×10⁻⁴ K⁻¹`) that
> the KPP / mixed-layer literature also calls `α`. Every buoyancy,
> buoyancy flux, N², reduced gravity or thermal-wind expression built from
> them therefore carries an explicit `1/ρ₀`:
>
> | quantity | correct form |
> |---|---|
> | buoyancy | `b = −g·ρ′/ρ₀` |
> | Brunt–Väisälä | `N² = −(g/ρ₀)·∂ρ/∂z` |
> | reduced gravity | `g′ = (g/ρ₀)·Δρ` |
> | surface buoyancy flux | `B_0 = (g/ρ₀)·(α_T·F_T − β_S·F_S)` |
> | thermal wind | `∂U/∂z = −(g/(ρ₀·f))·∂ρ/∂y` |
>
> with the kinematic surface fluxes `F_T = Q_heat/(ρ₀·c_p)` [K m/s] and
> `F_S = Q_salt/ρ₀` [PSU m/s]. Dropping the `1/ρ₀` scales the result by
> ρ₀ ≈ 1035 — it is invisible only when `α_T` is *itself* given the
> fractional magnitude, which is exactly what the `1.7e-4` default below
> does. KPP's `B_0` shipped without the `1/ρ₀` until it was fixed; the
> regression gate is `tests/test_ocean_buoyancy_flux.F90`, which also
> pins KPP's `B_0` and EPBL's `b0` to the same number.


The whole reference state — all five of `α_T`, `β_S`, `T_ref`, `S_ref`,
`ρ_0` — is set from **`&ocean_ic_nml`**:

| Knob | Default | Units |
|---|---|---|
| `alpha_T` | `1.7e-4` | kg/m³ per °C |
| `beta_S` | `7.6e-4` | kg/m³ per PSU |
| `T_ref` | `10.0` | °C |
| `S_ref` | `35.0` | PSU |
| `rho_0` | `1035.0` | kg/m³ |

> **UNITS TRAP — `α_T`/`β_S` are DIMENSIONAL, not fractional.** The form
> above is the density-ANOMALY one, so the coefficients carry kg/m³ per
> unit T/S. Most protocols quote the FRACTIONAL coefficients of the
> equivalent `ρ = ρ_0·(1 − α·(T−T_ref) + β·(S−S_ref))`, with `α` in
> 1/°C and `β` in 1/PSU. **Convert by multiplying by `ρ_0`:**
> `alpha_T = ρ_0·α`, `beta_S = ρ_0·β`. For example ISOMIP+ (Asay-Davis
> et al. 2016: `α = 3.733e-5` 1/°C, `β = 7.843e-4` 1/PSU,
> `ρ_0 = 1027.51`) becomes
> `alpha_T = 3.8356948e-2`, `beta_S = 8.0587609e-1`. Feeding the
> fractional numbers straight in under-states the density response
> ~1000× — a plausible but far too weakly stratified run. Locked by
> `tests/test_ocean_linear_eos_knobs.F90`.

> **The `&tracer_nml` spellings are RETIRED.** `&tracer_nml alpha_T`,
> `beta_S`, `T_ref` and `S_ref` only ever reached
> `tracer_t%eos_coeff`/`eos_ref`, which no ocean kernel reads — setting
> them was silent. Moving any of them off its historical default is now
> a **fail-loud `validate_config` error** naming the `&ocean_ic_nml`
> replacement. Delete the key; it never did anything.

The shipped defaults are deliberately small (they suppress baroclinic
feedback from PPM round-off in tests that don't care about realistic
density gradients); density-driven cases must set `alpha_T`/`beta_S`
explicitly. Setting `α_T = β_S = 0` decouples the tracers from the
dynamics entirely — though the adiabatic reduced-gravity setup used by
the MOM6-reference double gyre gets there instead via
`&ocean_thermo_nml enable_thermodynamics = .false.` plus the gprime PGF,
which never calls the EOS at all.

> **`&ocean_ic_nml rho_0` (default `1035.0` kg/m³) is the single ρ₀ of
> record.** It is the linear-EOS reference density AND the Boussinesq divisor
> the pressure gradient uses (`du/dt = −(1/ρ₀)∂p/∂x`), plus the FV-MOM6
> anomaly baseline — `configure_ocean_pgf` copies `eos%rho0` into both PGF
> reference densities, so a run cannot end up with the EOS on one ρ₀ and the
> PGF on another. The same scalar reaches EPBL, kappa-shear, tidal mixing,
> wave speed, the `η_ib` surface-pressure seam, GM / MEKE / Redi / MLE and the
> isopycnal slopes, and — via `configure_ocean_reference_density`, the
> fan-out step that closes this list — the **surface heat/salt flux**
> (`dt/(ρ₀·cp)` and `dt/ρ₀` on every surface tracer source, sea-ice coupling
> included), the **surface wind stress** (`τ/(ρ₀·h_top)`, plain and
> DIRECT_STRESS), **KPP-vmix** (N², `u* = √(|τ|/ρ₀)`, the kinematic fluxes
> behind `B_0`, the PP81 and convective-adjustment N²) and the **geothermal**
> bed source (`dt·Q_geo/(ρ₀·cp)`). The implicit vertical-friction solve takes
> it as an argument from the wind-stress slot; omitting it with the implicit
> stress fold active now fails loud instead of falling back to a literal.
> Default `rho_0 = 1035` ⇒ every one of those is bit-identical to the old
> hard-coded value; a run that sets `rho_0 /= 1035` is answer-changing by
> `1035/ρ₀` on each of them, which is the bug being fixed, not a regression.
>
> Note `&nonhydrostatic_nml rho_0` is a DIFFERENT, dead knob (nothing reads
> it). Still NOT unified, deliberately: the `RHO_WATER` constant behind the
> console mass/heat/salt diagnostics and the OBC `mass_out` accounting,
> `&ocean_ice_nml rho_ocean` (the EVP ice-ocean drag reference), and the
> `ICE_RHO_*` constants. Those are separate physical constants with their own
> call sites, not copies of ρ₀.

#### Freezing point (liquidus) — `&ocean_eos_nml tfreeze_set`

`eos_freezing_point(eos, S, p)` evaluates the LINEAR liquidus

```
T_f = λ1 × S + λ2 + λ3 × p          (S in PSU, p in Pa, T_f in °C)
```

for **every** EOS variant — MOM6 keeps `TFREEZE_FORM = "LINEAR"` as its
default under any density branch, so this is parity, not a shortcut. The
three coefficients live on the EOS handle and are chosen as a **named
set**:

| `tfreeze_set` | λ1 (°C/PSU) | λ2 (°C) | λ3 (°C/Pa) | `T_f(34.5, 0)` | Source |
|---|---|---|---|---|---|
| `"seaice"` (**default**) | `-0.054` | `0` | `-7.53e-8` | −1.863 °C | SIS2/MOM6 sea-ice `T_Freeze` |
| `"isomip"` | `-0.0573` | `0.0832` | `-7.53e-8` | −1.89365 °C | ISOMIP+ — Asay-Davis et al. (2016) *GMD* **9**, 2471–2497, Table 4 p. 2483; consumed in their eq. (25) p. 2485 |

**Why it is selectable.** At S = 34.5 the two sets are **~0.031 °C**
apart — a few percent of a typical Antarctic thermal driving, and enough
to flip the **sign** of an ice-shelf basal melt rate over a 0.03 °C band
of ocean temperature. A sea-ice run wants the SIS2 number its column
model was ported and tested against; an ISOMIP+ cavity run is required by
its protocol to use the other. The difference is exactly
`(λ1ⁱˢᵒ − λ1ˢᵉᵃ)·S + λ2ⁱˢᵒ = −0.0033·S + 0.0832`, independent of `p`
(both sets share λ3).

Named sets only, on purpose: λ1/λ2/λ3 are a fitted triple, so there is no
free-form coefficient knob to mix λ1 from one paper with λ2 from another.
A **nonlinear** liquidus (MOM6 `MILLERO_78`, a TEOS-10
`t_freezing(SA, p)` polynomial) is a different functional *form* and would
arrive as its own selector at the documented seam in `eos_freezing_point`,
not as another member of this list.

An unrecognised value is a **fail-loud `validate_config` error**, not a
silent fallback — a mistyped liquidus has no run-time symptom. Default
`"seaice"` ⇒ bit-identical; the shipped consumers (frazil, frazil uptake,
basal flux) all call at `p = 0`, where the coefficients-on-the-handle
expression reproduces the previous one bit for bit. Locked by
`tests/test_ocean_freezing_point.F90`.

### Physics closures

The full per-closure detail — knobs, formulae, MOM6 parity notes,
per-path test names — lives in
[`docs/CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md); which
closure is active in which regime, and its tunable knobs, is tabulated in
[`docs/CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md) (drift-checked by
`tools/check_closure_matrix.py`). In brief:

- **Vertical mixing** — PP81 interior + KPP boundary overlay (default
  on), EPBL (`&ocean_epbl_nml`, mutually exclusive with KPP),
  kappa-shear (JHL08, `&ocean_kappa_shear_nml`), tidal mixing
  (`&ocean_tidal_mixing_nml`), convective adjustment (`&ocean_conv_nml`),
  double diffusion (`&ocean_ddiff_nml`), and two mutually exclusive
  background schemes (`&ocean_vmix_nml bkgnd_profile` Bryan-Lewis **xor**
  `bkgnd_henyey`). Every closure contributes into `kv`/`kt`, with
  `vmix_assemble` the single floors/ceilings/smoothing gate; `ks` is
  derived from `kt` by `vmix_split_kd_heat_salt` and equals it unless
  double diffusion is on.
- **Lateral closures** — Leith, Smagorinsky KH+AH, Leith-biharmonic,
  constant floors + per-cell CFL clamps, resolution-scaled viscosity,
  anisotropic + live velocity-scale (`&ocean_hvisc_nml`), MEKE
  backscatter (`&ocean_meke_nml`), GM (`&ocean_gm_nml`), Redi
  (`&ocean_redi_nml`), Fox-Kemper MLE (`&ocean_foxkemper_nml`), and the
  along-coordinate tracer Laplacian (`&ocean_hdiff_nml kappa_h`). The
  dispatcher is fail-loud; a biharmonic backstop with a NON-ZERO
  dissipation coefficient is mandatory under backscatter.
- **Bottom drag** (`&ocean_bdrag_nml`) — linear (Rayleigh) or quadratic
  (log-layer, default), each with an HBBL-distributed mode.
- **Implicit stress/drag fold** (`&ocean_vdiff_nml`) — folds wind stress
  and bottom drag into the backward-Euler vertical-friction solve,
  killing thin-layer CFL blow-up.
- **Surface forcing** — wind stress (2D `tau_x`/`tau_y`, with the MOM6
  `DIRECT_STRESS` distribution option), heat/salt flux and the
  MOM6-shaped component set (`&ocean_forcing_nml enable_components`),
  shortwave penetration (`&ocean_thermo_nml sw_pen_frac`), surface
  buoyancy restoring (`&ocean_restore_nml`), atmospheric pressure loading
  (`&ocean_psurf_nml`), geothermal bottom heat
  (`&ocean_geothermal_nml`), and time-varying NetCDF forcing
  (`&ocean_dataovr_nml`, off the shared `rdb_ocean_data_input` reader).
- **Tides** — equilibrium body forcing + scalar SAL
  (`&ocean_tides_nml`), the boundary-tide nodal correction, and
  barotropic linear wave drag (`&ocean_bt_nml wave_drag`).
- **Porous barriers** (`&ocean_porous_nml`, Adcroft 2013) — subgrid
  sill/strait blocking that narrows the continuity and barotropic face
  widths.
- **Interior land masking** (static free-slip walls via metric-zeroing)
  and **dynamic wet/dry** (`&ocean_wetdry_nml`; sigma / z*-lite,
  single-rank, positive-definite outflow limiter).
- **Sea ice** (`&ocean_ice_nml`, default off) — SIS2 port: Winton column
  thermodynamics, multi-category ITD, category transport, C-grid EVP
  dynamics, snowfall, Archimedes snow-ice flooding. Single-rank,
  wall/periodic edges, no tripolar fold.

Most new-physics knobs default off ⇒ bit-identical; KPP and the ALE remap
are the on-by-default exceptions. Each path is guarded by an analytical
test.

### Windowed tracer advection + the face-reconstruction ladder

`&ocean_vmix_nml dt_tracer_advect_ratio > 1` decouples the horizontal
tracer drain from the per-step dynamics (MOM6 `DT_TRACER_ADVECT`); the
drain runs a CW-PPM flux-form advection over the accumulated mass flux.
Single-rank.

`&ocean_vmix_nml tracer_recon` selects the FACE-VALUE reconstruction used
by that drain — the mass-flux split, the per-layer conservation shell and
the downstream clamp are unchanged, only the upwind face value swaps:

| `tracer_recon` | scheme | nghost | full-period L1 order | notes |
|---|---|---|---|---|
| `ppm` (default) | Colella-Woodward 1984 parabola | 2 | 2 | bit-identical to the pre-existing path |
| `weno5` | Jiang-Shu 1996 + Borges 2008 Z-weights, 3 quadratics | 3 | 3 | the recommended rung |
| `weno7` | Balsara-Shu 2000, 4 cubics, Z-weights | 4 | 4 | best front L1 |
| `weno9` | Balsara-Shu 2000, 5 quartics (simplified 3-term β), Z-weights | 5 | ~6 | documented-expensive |
| `pqm` | — | — | — | **not implemented, fail-loud** |

All rungs use the **CFL-swept-average** form (the face value integrates
the upwind reconstruction over the swept region, Godunov-style) — this is
single-stage stable, unlike the classical point-value WENO which is
weakly unstable in this integrator and needs a multi-stage RK. Z-weights
(Borges et al. 2008) are used throughout: better than Jiang-Shu weights
at coarse resolution, never worse.

**Boundary order-reduction ladder.** As a face approaches a wall / open
boundary and the wide stencil no longer fits, the rung degrades per face
`weno9 → weno7 → weno5 → plm → donor` (a stencil may reach at most one
ghost row, since deeper ghosts are zero-gradient copies that carry no
information). Guaranteed in-bounds and conservative; tested in
`test_ocean_tracer_weno`.

**Per-scheme `nghost` is fail-loud at configure** (weno5→3, weno7→4,
weno9→5): if `&grid_nml nghost` is smaller, `validate_config` aborts
naming both numbers. The default `nghost = 3` already satisfies
`ppm`/`weno5`; only `weno7`/`weno9` need an explicit extra column.

**GPU cost (V100, sm_70, the single WENO9-inclusive kernel).** 180
registers, **no spills** (16 KB stack frame, 0 spill stores/loads) ⇒ ~17%
theoretical occupancy — register-bound but not spilling (contrast the
255-register kappa-shear kernel). Cost rises with the rung's stencil
width; `weno5` is the accuracy-per-cost sweet spot. Default `ppm`
compiles the same kernel but never dispatches into it (bit-identical).

**Known-behaviour note — PPM overshoot.** The default PPM path is NOT
strictly monotone at sharp fronts: at a square-wave corner it overshoots
by ~1.8e-2 (≈3% of the jump), because the Colella-Woodward eq. 1.10
limiter without van-Leer-limited edge interpolation does not bound face
values by neighbour means. This is pre-existing, intended behaviour — the
downstream conservative tracer clamp is the safety net. The WENO rungs
overshoot two orders less and rely on the same clamp; no new limiter was
added.

> A separate `&tracer_nml tracer_recon` key survives from the coastal
> era. Its `weno*` rungs are rejected fail-loud on `sim_type='ocean'`, so
> on this tree that key is effectively pinned to `ppm`. Use the
> `&ocean_vmix_nml` knob above.

## 6. Validation Test Suite

Tests use the `test-drive` framework. PASS/FAIL against analytical
solutions or regression baselines. The suite is ~177 `rdb_*` ctest
entries over 180 test sources — 162 labelled `ocean`, 18 labelled `core`
(run a subset with `ctest -L ocean` / `-L core`). Regime labels are
mandatory per row in `tests/CMakeLists.txt` and enforced by the
`test-regime-labels` pre-commit hook.

| Test | What it tests | Notes |
|---|---|---|
| `test_ocean_analytical` | Exact-solution gates | vdiff vs erfc, gravity-wave phase speed, geostrophic-adjustment balance. |
| `test_ocean_validation` | Integrated invariants | Cooling-column heat budget, Ekman sign, conservation under no forcing. |
| `test_ocean_conservation_salt_heat` | Closed salt/heat budgets | Also the canonical `mem:separate` GPU device-mapping template. |
| `test_ocean_land_mask` | Interior land block stays exactly at rest | The land-mask trap. |
| `test_ocean_remap` / `test_ocean_remap_e2e` | ALE remap conservation | Drift ≤ 1e-15 over 10 remaps. |
| `test_vcoord_zstar_full` | z*-full column construction + layer protection | Surface anchoring, vanishing-layer floor. |
| `test_ocean_baroclinic_longrun` | Multi-day stability of the split-RK2 core | Bounded CFL, no drift. |
| `test_ocean_wetdry` / `test_ocean_wetdry_driver` | Dynamic wet/dry | Thacker moving shoreline; tracer conservation through a dry→wet→dry cycle. |
| `test_ocean_restart` | Restart round-trip | Bit-exact, including the wet/dry hysteresis registry. |
| `test_ocean_kpp` / `test_ocean_epbl` / `test_ocean_kappa_shear` | Boundary-layer + interior closures | Per-closure analytic anchors. |
| `test_ocean_tripolar` / `test_ocean_fold` / `test_ocean_bipolar` | Curvilinear grids + north fold | Vector sign flips, seam antisymmetry. |
| `test_halo_ocean_mpi` / `test_ocean_dyn_mpi` | C-grid halo exchange across ranks | MPI integration (`tests/mpi/`). |

Plus the per-closure, per-kernel and per-knob unit suites — PGF variants,
hvisc variants, vmix assembly, EOS, diag manager, sea ice, tides,
budgets, config schema, decomposition.

Canonical benchmark configs live under
[`validation_examples/ocean/`](../validation_examples/ocean):
`seamount/`, `geostrophic_adjustment/`, `eady/`, `eddy_test/`,
`double_gyre/`, `acc_channel/`, `neverworld2/`, `baroclinic_channel/`,
`island_at_rest/`, `flow_past_island/`, `double_drake/`, `sea_ice/`,
`tides/`, `sponge_demo/`, `ideal_age/`, `epbl_mld/`, and others. Each
carries a README with expected behaviour, analytical scales (where
applicable), and diagnostic interpretation of failure modes.

## 7. Runtime Configuration

All runtime parameters via Fortran namelist. Each namelist block is
optional — missing blocks fall back to defaults.

The tables below cover the **core** groups. The ~46 `&ocean_*_nml`
sub-namelists (one per physics concern — `coriolis`, `pgf`, `bt`,
`bdrag`, `hvisc`, `hdiff`, `vmix`, `vdiff`, `thermo`, `continuity`,
`topo`, `grid`, `ic`, `diag`, `bc`, `sponge`, `tides`, `epbl`,
`kappa_shear`, `meke`, `gm`, `redi`, `foxkemper`, `porous`, `wetdry`,
`ice`, …) are catalogued with defaults, units and descriptions in the
auto-generated
[`docs/generated_nml_knobs.md`](generated_nml_knobs.md), which is
regenerated from the live schema and is the authority for them.
`tools/nml_split.py` migrates a legacy monolithic `&ocean_setup_nml`.

### `&sim_nml`

| Parameter | Default | Description |
|---|---|---|
| `sim_type` | `"ocean"` | Simulation regime — `"ocean"` is the only value this build ships |

### `&grid_nml`

| Parameter | Default | Description |
|---|---|---|
| `nx`, `ny` | 200, 1 | Physical cells in x/y |
| `dx`, `dy` | 1.0, 1.0 | Cell size (m). **Derived** from `&ocean_grid_nml` when `axis_units = "degrees"` — do not pin them in that case. |
| `nghost` | 3 | Ghost cells per side (3 is the PPM-seam + wide-halo BT minimum) |

### `&time_nml`

| Parameter | Default | Description |
|---|---|---|
| `t_end` | 1.0 | Simulation end time (in `time_unit`) |
| `cfl` | 0.45 | CFL number for the adaptive timestep |
| `dt_max` | 1e10 | Maximum timestep (s) |
| `dt_fixed` | 0.0 | Fixed timestep. **Required `> 0`** — there is no adaptive-CFL helper on the ocean path. |
| `cfl_interval` | 1 | Recompute CFL every N steps |
| `time_unit` | `"s"` | Unit for the long-time fields: `s` / `min` / `hr` / `day` / `year` |

### `&physics_nml`

| Parameter | Default | Description |
|---|---|---|
| `manning_n` | 0.0 | Manning roughness coefficient |
| `wind_stress_x` | 0.0 | Wind stress x-component (Pa) |
| `wind_stress_y` | 0.0 | Wind stress y-component (Pa) |
| `coriolis_f` | 0.0 | Coriolis parameter f (1/s) |

### `&output_nml`

| Parameter | Default | Description |
|---|---|---|
| `output_to_file` | `.false.` | Enable file output |
| `output_dir` | `"./output"` | Output directory |
| `restart_interval` | 0.0 | Restart file interval (s; 0 = none) |
| `compress_output` | `.false.` | Enable deflate compression |
| `compress_level` | 1 | Deflate level (1–9) |
| `use_io_server` | `.false.` | Dedicated I/O rank per node. **Not supported on the ocean path** — the driver warns and the diag manager writes serial per-rank NetCDF. |
| `bathymetry_file` | `""` | NetCDF bathymetry file (empty = formula / flat) |
| `restart_file` | `""` | Warm start from restart file |

Snapshot output itself is driven by the diagnostics manager
(`&ocean_diag_nml`: `enabled`, `filename`, `dt_out`, the `diags`
selection knob, per-diag output vcoord remap, `output_precision`).
Per-rank files merge offline via `tools/merge_output.py`.

### `&boundary_nml`

| Parameter | Default | Description |
|---|---|---|
| `bc_west / east / south / north` | `"wall"` | BC type per edge (see §4; `&ocean_bc_nml` carries the per-edge ocean detail) |
| `n_tidal_constituents` | 0 | Number of tidal constituents (0 = legacy single) |
| `tidal_amp(10)` | 0.0 | Constituent amplitudes (m) |
| `tidal_phase(10)` | 0.0 | Constituent phases (rad) |
| `tidal_omega(10)` | 0.0 | Constituent angular frequencies (rad/s) |
| `inflow_salinity` | -1.0 | Inflow salinity (PSU); <0 → zero-gradient |
| `inflow_temperature` | -999.0 | Inflow temperature (°C); <0 → zero-gradient |
| `sponge_width` | 0 | Sponge layer width (cells) |
| `sponge_strength` | 0.0 | Sponge relaxation rate (1/s) |

### `&nonhydrostatic_nml`

The group name is historical; on this tree it carries the layer-stack
size and the PP81 / KPP constants.

| Parameter | Default | Description |
|---|---|---|
| `nz_layers` | 2 | Number of vertical layers |
| `use_multilayer` | `.false.` | Enable the coupled vertical layer stack |
| `rho_0` | 1000.0 | Reference density (kg/m³). **Dead on the ocean path** — nothing reads `cfg%rho_0`. The live reference density is `&ocean_ic_nml rho_0` (default 1035), which feeds the EOS *and* both PGF reference densities. |
| `kpp_ri_crit` | 0.3 | Critical bulk Richardson number for the KPP boundary-layer depth |
| `kpp_cs_nonlocal` | 6.3 | KPP non-local (counter-gradient) transport coefficient |
| `kpp_c_vt2` | 1.8 | KPP `V_t²` unresolved-turbulence coefficient (0 = off) |
| `hdiff_kappa` | 0.0 | Horizontal tracer diffusion (m²/s). The ocean along-coordinate equivalent is `&ocean_hdiff_nml kappa_h`. |

### `&tracer_nml`

Tracer initial conditions, linear-EOS coefficients, clamp bounds and
background vertical diffusivities. Running without this block is fine —
all fields have sensible defaults.

| Parameter | Default | Description |
|---|---|---|
| `initial_salinity` | 35.0 | Initial salinity (PSU, uniform IC) |
| `S_ref` | 0.0 | **RETIRED** — setting it fails loud; use `&ocean_ic_nml S_ref`. |
| `beta_S` | 0.78 | **RETIRED** — setting it fails loud; use `&ocean_ic_nml beta_S`. |
| `S_min` | 0.0 | Lower clamp (PSU) |
| `S_max` | 40.0 | Upper clamp (PSU) |
| `kappa_S_bg` | 1.0e-5 | Background vertical S diffusivity (m²/s) |
| `S_init_surface`, `S_init_bottom` | 0.0 | Surface (k=nz) / bed (k=1) salinity for the linear-in-layer stratified IC. Both must be non-zero (one alone fails loud); both zero ⇒ uniform `initial_salinity`. Stable polarity is the INVERSE of temperature's — salty water belongs at the bed, so `S_init_bottom > S_init_surface`. |
| `initial_temperature` | 15.0 | Initial potential temperature (°C, uniform IC) |
| `T_ref` | 15.0 | **RETIRED** — setting it fails loud; use `&ocean_ic_nml T_ref`. |
| `alpha_T` | 0.17 | **RETIRED** — setting it fails loud; use `&ocean_ic_nml alpha_T`. |
| `T_min` | -2.0 | Lower clamp (°C); seawater freezing |
| `T_max` | 40.0 | Upper clamp (°C) |
| `kappa_T_bg` | 1.0e-5 | Background vertical T diffusivity (m²/s) |
| `T_init_surface`, `T_init_bottom` | 0.0 | Surface (k=nz) / bed (k=1) temperature for the linear-in-layer stratified IC; both zero ⇒ uniform `initial_temperature` |
| `tracer_recon` | `"ppm"` | Legacy face-reconstruction selector; `weno*` is rejected fail-loud on `sim_type='ocean'`. Use `&ocean_vmix_nml tracer_recon`. |

The density kernel itself is selected by `&ocean_eos_nml`.

### `&vcoord_nml`

| Parameter | Default | Description |
|---|---|---|
| `vcoord_type` | `"sigma"` | `sigma`, `zsigma`, `zstar`, `zstar_sigma`, `zstar_full`, `z_fixed`, `rho`, `hycom`, `z`, `eulerian_z`, `isopycnal`, `lagrangian`, `gprime` |
| `thickness_config` | `"sigma"` | **Initial** layer thickness (distinct from `vcoord_type`, which selects the *running* coordinate). `sigma` = even split of the local depth, `h_layer = b/nz` (default ⇒ bit-identical). `uniform_z` = MOM6 `THICKNESS_CONFIG="uniform"` port: uniform z interfaces over the global `ocean_max_depth`, clipped bottom-up to the local bathymetry, sub-floor layers collapsed to `max(angstrom_h, 2·H_VANISHED)`. Exactly conservative (`Σh = b`). See the note below. |
| `remap_method` | `"ppm"` | `pcm`, `plm`, `ppm`, `ppm_h4`, `pqm` (conservative remap reconstruction; `pqm` = piecewise-quartic, 4th–5th order, W&A 2008, falls back to PPM for nz<5) |
| `zstar_h_surf_target` | 0.0 | Surface-anchored layer thickness for `zstar_full` (m; 0 = auto) |
| `zstar_h_min` | 1.0e-4 | Vanishing-layer thickness for `zstar_full` / `z_fixed` (m). Anti-zero armour for BELOW-BED filler layers, **not** a positivity floor: it belongs ≤ `H_VANISHED` (1.5e-4) so those layers keep reading as vanished. `validate_config` **refuses** ≤ 0 and **warns** above `H_VANISHED` (above it the filler goes dynamically live in EOS/PGF/remap-drain/vdiff while the coordinate still treats it as throwaway; it warns rather than aborts only because the Python worked example currently sits in that band). The `rho`/`hycom` regrid reads the same knob under the opposite contract — its collapsed layers carry tracer mass, so it floors at `max(zstar_h_min, 2·H_VANISHED)`. For a genuinely *live* minimum thickness use `&ocean_isopycnal_nml angstrom_h`. |
| `zstar_stretching` | `"log"` | `"log"` or `"uniform"` surface-concentration stretching |
| `zstar_n_surf` | 0 | Number of fine near-surface layers (`zstar_full`; 0 = auto) |
| `rho_ref_pressure` | 2.0e7 | Reference pressure (Pa) for the potential density that defines `rho` |
| `rho_target_light` | 1020.0 | Lightest (surface) target interface density (kg/m³) for the `rho` light→dense linspace fallback |
| `rho_target_dense` | 1030.0 | Densest (bed) target interface density (kg/m³) for the `rho` linspace fallback |
| `regrid_time_scale` | 0.0 | ALE regrid time-filter timescale (s). `>0` relaxes the coordinate toward the regrid target by `dt/(τ+dt)` each thermo step (MOM6 `REGRID_TIME_SCALE`), damping grid-motion shock / σ-PGE; `0` = jump to target (default, bit-identical). |
| `remap_vel_conserve_ke` | `.false.` | If true, rescale the remapped baroclinic velocity anomaly to conserve column KE (capped 1.25×; barotropic mean preserved) — MOM6 `REMAP_VEL_CONSERVE_KE`. Default off ⇒ momentum-only remap (bit-identical). |
| `rho_target(:)` | all-zero | Explicit monotone-increasing `nz+1` interface target densities (kg/m³); element 1 = surface. All-zero ⇒ build the `light→dense` linspace |

**When you need `thickness_config = "uniform_z"`.** Under a horizontally-uniform
density stack — the linear density-coordinate IC (`&ocean_ic_nml rho_lightest` /
`rho_range`, MOM6 `COORD_CONFIG="linear"`) with `enable_thermodynamics=.false.` —
the layer index *is* the density, so the layer interfaces *are* the isopycnals.
The default `"sigma"` seed then makes every isopycnal follow the bathymetry,
which stands the full `rho_range` contrast up across each shelf break at t=0.
On the MOM6 double-gyre spoon (`edge_depth=100`, `max_depth=2000`,
`rho_range=2`) that is a ~1.9 kg/m³ horizontal density jump — g′ ≈ 0.018 m/s²,
a ~1.3 m/s gravity current — over layers only `edge_depth/nz` thick. The rim
layers drain, the non-conservative `angstrom_h` floor starts injecting mass, and
the run grounds out. `"uniform_z"` seeds flat resting isopycnals instead, which
is the layered-model resting state MOM6 starts from. Pair it with
`vcoord_type="lagrangian"`. Mutually exclusive with `&ocean_wetdry_nml enable`
(fail-loud): that path pins its own emerged-column seed invariant.

### `&ocean_isopycnal_nml`

Grounding-stability controls for the Lagrangian (`vcoord_type="lagrangian"`)
vertical coordinate. All knobs default OFF ⇒ bit-identical.
Gates test `coord_type == VCOORD_LAGRANGIAN` specifically; other vcoords
(sigma, zstar, etc.) are unaffected.

| Parameter | Default | Description |
|---|---|---|
| `angstrom_h` | `0.0` | Minimum-thickness floor (m) on the Lagrangian continuity h-update — MOM6 GV%Angstrom_H analogue: `h_new = max(h - dt·div, angstrom_h)`. 0 = off (bit-identical). Recommended run value 1e-3..1e-2 m. **Non-conservative**: injects ≤ angstrom_h×cell-area per floored layer-cell; the mass `Error` diagnostic reflects the true divergence. R7 caveat: lifts h but not hTr, so Tr = hTr/h shifts on a floored layer — harmless for adiabatic isopycnal; thermo-on isopycnal correctness is out of scope for v1. |
| `reset_vanished_u` | `.false.` | Zero the face velocity of a layer vanished on **both** adjacent cells (h ≤ max(angstrom_h, H_VANISHED) on both sides). One-sided faces — a massive-to-thin grounding front — are unchanged so the re-wetting flux survives. |
| `cfl_ignore_vanished` | `.false.` | Exclude vanished layers (same both-sided criterion) from the console MaxCFL / CFL-panic / `apply_velocity_truncation`. A velocity spike in a vanished layer cannot falsely trigger an abort; the vanished face is zeroed (not CFL-clipped). |

### `&initial_condition_nml`

| Parameter | Default | Description |
|---|---|---|
| `h0` | 1.0 | Background depth (m) |

Ocean initial conditions proper are set by `&ocean_ic_nml`,
`&ocean_zinit_nml` and `&ocean_topo_nml`.

### `&mpi_nml`

| Parameter | Default | Description |
|---|---|---|
| `px`, `py` | 1, 1 | MPI process grid dimensions (0 = auto-factor) |

`&ocean_mpi_nml poison_ghosts` sentinel-NaNs the exchange-covered ghost
bands at each outer-step start, so consuming an unexchanged ghost becomes
a loud NaN instead of a plausible number.

### `&logging_nml`

| Parameter | Default | Description |
|---|---|---|
| `log_level` | `"info"` | `debug`, `verbose`, `info`, `performance`, `warning`, `error` |
| `status_interval` | 0.0 | Status print cadence (in `time_unit`; 0 = every 100 steps) |

### Example

Abridged from the canonical MOM6-reference double gyre
(`validation_examples/ocean/double_gyre/double_gyre_mom6.nml`), which runs
stable past day 580 on a single V100:

```fortran
&sim_nml
  sim_type = "ocean"
/

&grid_nml
  nx = 44
  ny = 40
  nghost = 3          ! dx/dy are DERIVED from &ocean_grid_nml degrees
/

&ocean_grid_nml
  axis_units = "degrees"
  len_lon    = 22.0
  len_lat    = 20.0
  rad_earth  = 6.378e6
/

&time_nml
  t_end     = 10.0
  time_unit = "day"
  dt_fixed  = 1200.0
/

&physics_nml
  coriolis_f = 9.4e-5
/

&nonhydrostatic_nml
  nz_layers = 2       ! the gprime PGF below is NK=2-only
/

&ocean_topo_nml
  topo_config    = "spoon"
  max_depth      = 2000.0
  edge_depth     = 100.0
  slope_scale    = 400000.0
  wind_config    = "2gyre"
  taux_magnitude = 0.1
  coriolis_beta  = 1.76e-11
/

&ocean_coriolis_nml
  form = "sadourny_energy"
/

&ocean_pgf_nml
  form        = "gprime"
  gprime_gfs  = 0.98
  gprime_gint = 0.0098
  maxvel      = 6.0
/

&ocean_bdrag_nml
  form   = "linear"
  r      = 2.5e-5
  hbbl   = 10.0
  bg_vel = 0.1
/

&ocean_hvisc_nml
  nu_h            = 10000.0
  smag_ah         = .true.
  lateral_closure = "smagorinsky"
  c_smag          = 0.15
/

&ocean_bt_nml
  auto_n_inner  = .true.
  cfl_bt_safety = 0.65
  bebt          = 0.2
/

&vcoord_nml
  vcoord_type         = "zstar_full"
  zstar_h_surf_target = 1000.0
  zstar_h_min         = 1.5e-4
/

&ocean_diag_nml
  enabled  = .true.
  filename = "double_gyre"
  dt_out   = 1.0      ! one frame per day (time_unit = "day")
/

&logging_nml
  log_level       = "info"
  status_interval = 1.0
/
```

## 8. Key Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Language | Fortran | Numerical performance, scientific community familiarity |
| GPU strategy | `do concurrent` + OpenACC | Single source for GPU, multicore, serial; no vendor lock-in |
| Reductions | `!$acc parallel loop reduction(...)` | Inert comment on non-OpenACC compilers; correct sequential fallback |
| GPU data | `!$acc enter/exit data`, `!$acc update` | Portable across OpenACC-capable compilers; inert elsewhere. `!$omp target` is used for the same effect where AMD/Intel portability matters |
| Compilers | NVHPC, gfortran, ifx, Cray ftn | Mixed-hardware HPC portability |
| Spatial scheme | Continuity-PPM + PV-conserving Coriolis on an Arakawa C-grid | MOM6-style; eddy-resolving-ready |
| Continuity | Transport equation, thickness prognostic | No Poisson constraint, no FFT projection; frees the vertical coordinate for ALE |
| Time scheme | Split-explicit RK2 with a sub-cycled barotropic fast loop | Outer `Δt` set by advective CFL, not gravity waves; no global solves |
| Vertical coordinate | ALE (Lagrangian-then-remap) | One remap path serves every `VCOORD_*` family |
| Data layout | SoA, `(nx, ny, nz)` | Coalesced GPU memory access on the first index |
| MPI | Routed through `pic_mpi_lib` only | Lint-enforced (`tools/no_mpi_in_rdb.sh`); never `use mpi` / `use mpi_f08` directly |
| Comm backend | One implementation in `src/comm/` | No stub twin; pic-mpi picks MPI vs its serial backend (`PIC_ENABLE_MPI`). Call sites must not reach point-to-point at one rank — see `src/comm/README.md` |
| I/O strategy | Per-rank files + offline merge (`tools/merge_output.py`) | Eliminates the MPI gather |
| Build system | CMake | Multi-configuration support, industry standard |
| Configuration | Fortran namelist, schema-validated | Zero dependencies, native to Fortran; the schema generates `docs/generated_nml_knobs.md` |
| Dependencies | `pic`, `pic-mpi`, `test-drive`, NetCDF-Fortran | Reuse tested types, MPI wrappers, testing |
