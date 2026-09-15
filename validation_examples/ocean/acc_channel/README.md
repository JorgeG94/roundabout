# Idealized ACC channel — the all-features ocean showcase

A re-entrant (periodic-x) **spherical** Southern-Ocean channel with a
seamount ridge, built to exercise the full ocean feature stack in a
single run. It is the closest thing in `validation_examples/` to an
end-to-end "everything on" integration case for the C-grid dyn-core.

## What it is

- **Spherical curvilinear channel** centred at ~45 S: `lon_west = 0`,
  24 deg lon x 14 deg lat, `lat_south = -52`, 0.5 deg resolution
  (48 x 28 cells), `nghost = 3` (periodic-x needs the deeper seam
  stencil). `nz = 16` layers over `H = 3000 m`.
- **Seamount ridge** mid-channel, a centred Gaussian bump rising to
  1500 m (50% of depth). (See the ridge caveat below — it is a radial
  bump, not a full meridional wall.)
- **Steady zonal westerly wind** (0.12 Pa) + **weak surface cooling**
  (-40 W/m^2), so the EPBL surface layer gets both wind and convective
  TKE. **Planetary Coriolis** (f = 2*Omega*sin(lat); southern
  hemisphere => f < 0; the closures take |f|).
- **Linear EOS** (`alpha_T = 0.2`, matching the `epbl_mld` thermal
  pattern), stratified linear T(z) IC (2 degC bed -> 8 degC surface),
  uniform S = 35.

## Features exercised (all alive in the 30-day output)

| Feature | Knob | Verified |
|---|---|---|
| Spherical curvilinear metrics | `&ocean_grid_nml grid_config="spherical"` | quiescent gate passes (1e-11 m/s) |
| Planetary Coriolis | `coriolis_scheme="planetary"` | zonal jet develops |
| Periodic-x OBC (re-entrant) | `&ocean_bc_nml west/east="periodic"` | seam invisible (ghost == wrap exactly) |
| Walls north/south | (default) | channel confined |
| EPBL surface boundary layer | `&ocean_epbl_nml enable=.true.` | MLD_EPBL 110-200 m, evolving |
| kappa-shear interior mixing | `&ocean_kappa_shear_nml enable=.true.` | Kd_KSHEAR up to 1.7e-2 m^2/s |
| Smagorinsky KH + AH | `&ocean_hvisc_nml lateral_closure="smagorinsky", smag_ah=.true.` | stable |
| z* vertical coord + ALE remap | `&vcoord_nml vcoord_type="zstar"` | mass conserved 1e-16 |
| Distributed linear bottom drag | `&ocean_bdrag_nml form="linear", hbbl=10` | stable |
| Diag manager | `&ocean_diag_nml` | SSH/T/S/u/v/KE + MLD_EPBL + Kd_EPBL + Kd_KSHEAR |

## How to run

Build once — load your NVHPC + NetCDF toolchain however your site does it
(`module load`, Spack: `environments/spack_env_*.yaml`, conda). Use **exactly
one** toolchain per shell; stacking a gfortran and an NVHPC environment puts
two incompatible NetCDF builds on the link line.

```bash
cmake --build build -j 40
```

Then, from this directory:

```bash
cd validation_examples/ocean/acc_channel
# (NVHPC + NetCDF toolchain loaded in this shell — and only that one)

# 1. Quiescent tier-1.5 trap FIRST (5 days, must stay at rest):
CUDA_VISIBLE_DEVICES=<gpu> ../../../build/rdb acc_channel_quiescent.nml

# 2. The forced showcase (30 days, ~17 s on a V100):
CUDA_VISIBLE_DEVICES=<gpu> ../../../build/rdb acc_channel.nml

# 3. The 3-panel movie (SSH+flow, MLD_EPBL, Kd_KSHEAR):
python3 animate_acc_channel.py
```

### One-knob Langmuir variant

Set `use_lt = .true.` in `&ocean_epbl_nml` to engage the LF17 wind-only
Langmuir-turbulence enhancement (COARE 3.5 + Phillips-spectrum Stokes,
Reichl & Li 2019 mstar boost) — no wave model needed. It deepens the
EPBL MLD modestly under the steady wind. Nothing else changes.

## Expected physics (30-day run, V100, 2026-06-11)

- **(a) Zonal jet**: mean surface u = +0.040 m/s (max 0.149), eastward
  through the channel; depth-mean u = +0.040 m/s. Crude depth-integrated
  zonal transport through a mid-section ~ O(100) Sv (over-read because
  the layer-thickness weighting here is approximate).
- **(b) Standing meander**: surface speed and v-structure differ
  upstream vs downstream of the ridge (mean surface speed 0.060 upstream
  vs 0.054 downstream; downstream surface |v| reaches 0.114 m/s vs 0.069
  over the crest) — the ridge steers the flow.
- **(c) MLD_EPBL**: 110-200 m, spatially varying and time-evolving
  (mean deepens ~180.4 -> 181.8 m over the run). It is alive but **near
  the first interface (187.5 m layers)** at this vertical resolution, so
  the MLD is effectively quantized to ~1 layer. For a crisp MLD-depth
  field, use thinner surface layers (see caveats).
- **(d) Kd_KSHEAR**: grows from 0 (day 1) to a max of 1.65e-2 m^2/s
  (day 30), concentrated at the **upper interface** in the wind-driven
  surface shear. The 187.5 m interior layers are too coarse to resolve a
  sheared ridge-wake interior, so kappa-shear fires near the surface
  rather than as a deep wake. Alive and evolving, but resolution-limited.
- **(e) Periodic seam invisible**: the periodic ghost column equals its
  wrapped partner exactly (max|delta| = 0.0 in u); physical SSH at
  column 1 vs column nx agree to within 6e-4 m (mean) / 2.5e-3 m
  (per-row max). No discontinuity across the x-edges.
- **(f) Stratification maintained**: per-layer mean T holds its profile;
  the surface-to-bed contrast goes 6.00 -> 5.83 degC (surface cooled
  from 8.0 to ~7.84 by the -40 W/m^2 flux + EPBL, exactly as intended).
  No column-integrated T blowup. Mass conserved to 1e-16, salt to 1e-13,
  heat budget closes to the imposed surface cooling (-0.17% over 30 d).

CFL stayed bounded < 0.012 throughout; `auto_n_inner` settled at
n_inner = 10 (c_ext ~ 171 m/s at H = 3000 m). The production envelope
(nu_h = 10000, Smag KH+AH, distributed linear drag, dt_fixed = 1200 s)
was used **as-is** from `double_gyre_mom6.nml` — **no envelope changes
were needed** for stability.

## Quiescent gate (the tier-1.5 trap)

`acc_channel_quiescent.nml` is the same grid / seamount bathymetry /
periodic-x BCs / closures with **zero wind, zero heat flux, and uniform
T,S** (no available potential energy). It stays at rest:
**max|u|,|v| = 1.0e-11 / 5.5e-12 m/s, En = 6e-25 J, MaxCFL = 0** over
5 days — i.e. roundoff. This proves the spherical metrics, planetary
Coriolis, periodic-x seam, and seamount bathymetry are mutually
consistent.

> **Finding — stratified resting state drifts over the ridge.** A
> *stratified* resting state (the `acc_channel.nml` T(z) IC with the
> forcing turned off) does **not** stay at rest: it drifts to ~0.16 m/s
> centred on the ridge crest, decaying smoothly to ~2e-4 at the
> periodic seam and the walls. This is **not** a spherical-metric or
> periodic-seam bug — it is the known **sigma/z* pressure-gradient
> residual over a slope under stratification (APE-over-slope)**:
>
> - It is reproduced **identically on a Cartesian-grid control** (same
>   seamount, same strat), so it is not spherical-specific.
> - sigma and z* give the same drift (at SSH = 0 they produce the same
>   layout over this bathymetry).
> - The `set_bathymetry_seamount` docstring itself specifies "uniform
>   T,S, no APE" for the zero-motion expectation.
>
> **Caveat added with the Montgomery PGF fix — this finding does not
> reproduce from this namelist's IC.** `T_init_surface`/`T_init_bottom`
> are applied **per layer index**, not per depth
> (`t_layer(k) = T_init_bottom + dT_dlayer*(k-1)`,
> `rdb_ocean_state.F90`), so every layer is isothermal along itself and
> the state is isopycnal **by construction** — there is no horizontal
> density gradient for any vcoord to mis-integrate. Measured on current
> HEAD, wind and `q_heat` zeroed, 5 days: `En = 0.000E+00` and
> `MaxCFL = 0.00000` for sigma AND z*, under all four of `fv_lite`,
> `mont`, `fv_wright` and `fv_mom6`. The ~0.16 m/s drift reported above
> must have come from a different stratification path; it is NOT
> evidence about any PGF form. Re-derive it from a depth-based T(z) IC
> before citing it again.

> The quiescent namelist therefore uses uniform T,S (the honest trap).
> The drift is a closure accuracy limitation of terrain-following
> coordinates over steep topography, tracked elsewhere in the ocean
> backlog — not something to tune away here.

## Caveats / follow-ups

- **Wind profile is constant (depth-uniform zonal stress)**, not a
  sinusoidal ACC jet. The ocean path's only meridionally-varying wind
  generator is the closed-basin `2gyre` cosine, which is not a channel
  westerly. A `tau_x(y)` jet (e.g. a single sin^2 hump centred on the
  channel) would be more ACC-like but needs a new `wind_config` kernel —
  deferred (no `src/` change in this pass).
- **The ridge is a radial Gaussian bump**, not a full meridional wall
  spanning the channel. `set_bathymetry_seamount` builds an isolated
  centred bump. A true cross-channel ridge (a function of lon only) would
  give a cleaner standing-meander signal but likewise needs a new topo
  generator. The radial bump still produces a clear up/downstream
  asymmetry.
- **Vertical resolution is coarse** (187.5 m layers): MLD_EPBL is
  quantized to ~1 layer and kappa-shear can only resolve the surface
  shear, not a deep ridge wake. Both features are demonstrably alive;
  for sharper MLD / interior-mixing structure, refine the surface layers
  (a z* layout with thin top layers, or more `nz_layers`).


## Eddy-resolving stress test

`acc_channel_eddy.nml`: the same channel at 0.1 deg / 50 layers
(240 x 140 x 50 = 1.7M cells, ~13x the showcase) — the single-GPU
compute stress test. For the full-fat 10M-cell version scale nx/ny x2
AND dx/dy /2 together (0.05 deg; scaling only nx/ny grows the domain
toward the poles and dies on the lon-lat CFL collapse — see the
design doc). Budget GPU memory for the 10M-cell case. Not yet
run-validated at scale; expect to tune dt/nu_h if the first day
misbehaves.

## "Kitchen sink" — exercise every new closure (profiling)

`acc_channel_kitchensink.nml`: the same channel on a deliberately
moderate grid (240 x 140 x 50 = ~1.7M cells) that turns on **every**
ocean-physics path added in the MOM6-parity push in ONE run — EPBL +
kappa-shear + **Fox-Kemper MLE** + **geothermal** + **PPM_H4 remap** +
**advective-CFL truncation** (and the B1 wave-speed slot, dormant until
its B2 consumer lands). A 1-day pass turns over in ~5 s on a V100, so
it's a fast `nsys`/profiler target for the new code, not a validation
case (the physics combos are stacked for coverage, not realism — use
the per-feature examples to validate any one closure). Confirmed stable:
mass/salt conserved to round-off, `ocean_ale_remap` (~14%) and
`ocean_foxkemper` (~1%) show in the profiler breakdown.

## WENO windowed-drain tracer advection

`acc_channel_weno5.nml`: the same channel physics, but the horizontal
tracer advection is decoupled from the per-step dynamics and run once per
window with a **WENO5-Z** face reconstruction. The key knobs
(`&ocean_vmix_nml`):

```
tracer_recon           = 'weno5'   ! WENO5-Z swept-average drain
dt_tracer_advect_ratio = 2         ! the drain fires only at ratio > 1
dt_therm_ratio         = 2         ! must be an integer multiple of the above
```

**WENO engages ONLY inside the windowed drain.** At the default
`dt_tracer_advect_ratio = 1` the tracer advection runs every dynamic step
via the per-step CW split and `tracer_recon` is never consulted — so
`ratio > 1` is **required** for WENO to take effect. Because the windowed
drain is thermo-cadence-tied, `dt_therm_ratio` must be an integer multiple
of `dt_tracer_advect_ratio` (fail-loud at configure), which is why this
variant also bumps `dt_therm_ratio` to 2. The base grid already carries
`nghost = 3` (periodic-x), which satisfies the weno5 halo requirement;
**weno7 needs `nghost >= 4`, weno9 `nghost >= 5`** — both fail-loud. This
is the SEPARATE ocean knob in `&ocean_vmix_nml`, independent of the coastal
`&tracer_nml tracer_recon`.

```bash
cd validation_examples/ocean/acc_channel
# (NVHPC + NetCDF toolchain loaded in this shell)
CUDA_VISIBLE_DEVICES=<gpu> ../../../build_gpu/rdb acc_channel_weno5.nml
```

**Reduced grid (load-bearing).** This variant runs at **120 x 70 x 20**,
not the canonical 480 x 280 x 75. The full-res showcase maps ~14 GB of
ocean state and grows the lazy kernel workspaces past 32 GB, so it **OOMs
on a single 32 GB V100 at the first kernel launch** — this is a pre-existing
case-sizing property, **not a WENO effect**: the stock PPM `acc_channel.nml`
OOMs identically at that resolution on a V100. 120 x 70 x 20 fits (~0.8 GB
mapped state, ~23 GB peak with workspaces) and turns over a 1-day run in
~3 s. Scale nx/ny/nz back up on an A100/H100 for the full-res showcase.

**Expected physics (1-day, V100, reduced grid):**
- Stable: `MaxCFL ~ 2.3e-4`, no NaN, no OOM.
- Conservation intact: mass error ~5e-15, salt error ~2e-14 (round-off);
  heat error ~-5.6e-5 is the imposed -40 W/m^2 surface cooling, not a leak.
- **No new tracer extrema**: T stays in [2.00, 8.00] degC (the IC bounds),
  S in [34.998, 35.002] — the WENO-Z weights are monotone, so the drain
  adds no over/undershoot.
- vs a same-resolution / same-`ratio` PPM control the tracer fields are
  nearly identical here: at this coarse grid with `smagorinsky` (c_smag =
  0.15) the lateral closure dominates the 1-day tracer evolution, so the
  reconstruction difference is small. For a front-contrast demonstration,
  soften the closure and lengthen the run (the coastal `river_plume` WENO
  example shows the sharper-front effect where diffusion is weaker).
