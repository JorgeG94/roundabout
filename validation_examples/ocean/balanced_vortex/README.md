# Isolated balanced vortex on a beta-plane (eddy-tilt experiment)

A single mesoscale warm- or cold-core eddy is started in geostrophic /
thermal-wind balance (`&ocean_ic_nml ic_config = "balanced_vortex"`), then
left to drift on a mid-latitude beta-plane. Question: do the upper and lower
parts of the vortex move at different speeds (a vertical TILT), and what does
that do to the SURFACE vorticity? Hydrostatic is the right model here (aspect
ratio ~ 0.01, Rossby number ~ 0.1).

Physics background, stated without attributing more than is certain:
beta drift of isolated vortices -- cyclones move north-west, anticyclones
south-west in the northern hemisphere (McWilliams & Flierl 1979); baroclinic
vortices on a beta-plane (Chassignet & Cushman-Roisin 1991); beta-induced tilt
and vertical alignment (Sutyrin & Flierl 1994); review in Carton (2001).

## Construction (see `docs/REFERENCE.md`, "Isolated balanced vortex")

`T = T_surf - dT_dz*d + dT*exp(-r^2/R^2)*exp(-d^2/Hv^2)` (`d` = depth), linear
EOS, uniform S. The hydrostatic potential `phi = p'/rho_0` comes from the
density anomaly integrated down from the surface; **the SSH is chosen so the
bottom pressure anomaly is exactly zero** (no deep flow, a surface-intensified
vortex; warm core -> SSH bump). Face velocities are the centred differences of
the corner-averaged layer potential divided by the LOCAL `f0 + beta*(y - y_ref)`
of each face; error vs the analytic speed is second order in `dx/R`.
`vortex_balance` = `"geostrophic"` (default), `"gradient"` (cyclostrophic
quadratic, fails loud with no real root) or `"none"` (same density + SSH at
rest -- the unbalanced start).

## Cases

| File | What |
|---|---|
| `vortex_beta.nml` | warm core (+2 K, anticyclone), beta = 2e-11, 40 days |
| `vortex_beta_cyclone.nml` | cold core (-2 K, cyclone) twin |
| `vortex_fplane.nml` | the warm vortex with beta = 0: must stay put (control) |

All: 400 km x 400 km closed box at 5 km (80x80), H = 4000 m flat, **40 layers**
(100 m; the Hv = 700 m core spans ~10 layers; sigma = z on a flat bed),
`pred_corr`, **`dt = 600 s`**, f0 = 1e-4, linear EOS (`alpha_T = 0.17`,
`rho_0 = 1035`), `dT_dz = 4e-3` (N = 2.5e-3 1/s), R = 50 km, Hv = 700 m,
`nu_h = 100`. Peak surface speed ~0.3 m/s, peak `zeta/f0 = -0.27` (warm),
first-baroclinic deformation radius N H/(pi f) ~ 32 km.

Two constraints that bit during development:

- **Static stability of a cold core.** The anomaly `-|dT| exp(-d^2/Hv^2)` warms
  with depth by up to `0.858 |dT|/Hv` (2.45e-3 K/m here), so `vortex_dT_dz` must
  exceed that or `N^2 < 0` in the core. The IC refuses it (message quotes the
  minimum). The first cyclone run used `dT_dz = 2e-3`, below the minimum: that
  statically unstable start is what produced a spurious "intensification"
  (below).
- **Internal-wave CFL.** With N = 2.5e-3 the first internal mode runs at
  c1 = N H/pi = 3.2 m/s; `dt = 1200 s` at 5 km blows up within ~40 steps,
  `dt = 600 s` is stable (tested at 20 and 40 layers). The examples were checked
  to run on CPU for a few steps only; the 40-layer, 40-day runs are for the GPU.

```bash
./build/rdb validation_examples/ocean/balanced_vortex/vortex_beta.nml
python3.10 validation_examples/ocean/balanced_vortex/analyze_vortex_tilt.py \
    out_vortex_beta/vortex_beta_rank_000000.nc \
    --dx 5000 --nghost 4 --depth 4000 --radius 50000 --tilt-depth 700
```

40 days at dt = 600 s is 5760 steps; the 20-layer version took ~15 min per
2900 steps on one CPU core (serial `do concurrent` build), the 40-layer files
roughly twice that per step: intended for the GPU. The script needs only numpy + netCDF4 (matplotlib optional, for `--png`).
It tracks the surface centre by the same-sign vorticity centroid, then reports,
per day, the surface and deep-layer centres (default: centroid of the
temperature anomaly against the layer's outer-rim mean, `--method zeta` uses
vorticity), the offset, `slope = offset/dz`, the tilt angle from the vertical
and its azimuth (direction of the UPPER part relative to the deep part, degrees
CCW from east), and the surface vorticity extremum and same-sign circulation.
Multi-rank output must be merged first with `tools/merge_output.py`.

## Measured (CPU, nvfortran)

Static IC checks (`tests/test_ocean_balanced_vortex.F90`): velocity error
vs analytic `0.99 (dx/R)^2` of the peak speed, falling 4.0x when dx is halved;
bottom pressure anomaly < 1e-3 of `g*eta`; f-plane 2 days (R = 40 km, dx = 10 km):
centre moves 0.000 cells, surface KE +1.2 %, anisotropy 0.006 -- whereas the same
density and SSH started from rest ends at 47 % of the balanced surface KE
(-53 %), so the balance test is not vacuous. A cold core below the minimum
`dT_dz` is refused.

**Earlier 40-day runs (20 layers, `dT_dz = 2e-3`, `dt = 1200 s`; NOT the shipped
configuration).** Warm core: drifts SW, centre (200, 200) -> (177.1, 166.3) km
(dx = -22.9, dy = -33.7); f-plane control stays at (200.0, 200.0) with
circulation within 2 %. The 100 m centre leads the 700 m centre by 0.1 km (day 1),
0.9 (day 4), 1.4 (day 8), 2.2 (day 16), 3.0 (day 28), 2.8 km (day 40) while the
azimuth of the upper part relative to the deep part rotates from east (352 deg)
to south (261 deg); the offset to 1500 m reaches 46 km, but that deep anomaly
(about 4 % of the surface anomaly's weight) is a westward-dispersing wake, not a
coherent tilted core. Surface `zeta/f0` of the warm core: -0.273 -> -0.219
(day 16) -> -0.295 (day 40). These were at 200 m vertical resolution, which is
too coarse to resolve a tilt of a 700 m core; rerun with the shipped 40-layer
files before quoting tilt numbers.

**The cyclone "intensification" was a statically unstable IC.** That earlier
cold-core run (same coarse configuration) gained 68 % domain energy
(`En` 2.214e-4 -> 3.717e-4 m^2/s^2 by day 20) and its peak `zeta/f0` tripled
(0.25 -> 0.73): `dT_dz = 2e-3` was below the 2.45e-3 K/m minimum, so the core
was convectively unstable. With `dT_dz = 4e-3` (20 layers, `dt = 600 s`) the
same cyclone has `En` 2.214e-4 (day 0), 1.984e-4 (5), 1.863e-4 (10), 1.771e-4
(15), 1.691e-4 (20): monotonic viscous decay (-24 %), no growth. Old for
comparison: 2.214e-4, 2.007e-4 (5), 2.377e-4 (10), 3.322e-4 (15), 3.717e-4 (20).

## Caveats

- Cartesian beta-plane, flat bed, linear EOS only (refused otherwise).
- The 400 km box and `nu_h = 100` are sized for a 40-day drift of ~40 km;
  the closed walls start to matter for longer runs.
- Single-rank analysis script; merge multi-rank output first.
