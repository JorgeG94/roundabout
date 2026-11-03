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

All: 400 km x 400 km closed box at 5 km (80x80), H = 4000 m flat, 20 layers
(sigma = z on a flat bed), `pred_corr`, `dt = 1200 s`, f0 = 1e-4, linear EOS
(`alpha_T = 0.17`, `rho_0 = 1035`), N = 1.8e-3 1/s (`dT_dz = 2e-3`), R = 50 km,
Hv = 700 m, `nu_h = 100`. Peak surface speed ~0.3 m/s, peak `zeta/f0 = -0.27`
(warm), first-baroclinic deformation radius N H/(pi f) ~ 23 km.

```bash
./build/rdb validation_examples/ocean/balanced_vortex/vortex_beta.nml
python3.10 validation_examples/ocean/balanced_vortex/analyze_vortex_tilt.py \
    out_vortex_beta/vortex_beta_rank_000000.nc \
    --dx 5000 --nghost 4 --depth 4000 --radius 50000 --tilt-depth 700
```

40 days is ~2900 steps, ~15 min per case on one CPU core (serial `do concurrent`
build). The script needs only numpy + netCDF4 (matplotlib optional, for `--png`).
It tracks the surface centre by the same-sign vorticity centroid, then reports,
per day, the surface and deep-layer centres (default: centroid of the
temperature anomaly against the layer's outer-rim mean, `--method zeta` uses
vorticity), the offset, `slope = offset/dz`, the tilt angle from the vertical
and its azimuth (direction of the UPPER part relative to the deep part, degrees
CCW from east), and the surface vorticity extremum and same-sign circulation.
Multi-rank output must be merged first with `tools/merge_output.py`.

## Measured (CPU, nvfortran, 2026-11)

Static IC checks (`tests/test_ocean_balanced_vortex.F90`): velocity error
vs analytic `0.99 (dx/R)^2` of the peak speed, falling 4.0x when dx is halved;
bottom pressure anomaly < 1e-3 of `g*eta`; f-plane 2 days (R = 40 km, dx = 10 km):
centre moves 0.000 cells, surface KE +1.2 %, anisotropy 0.006 -- whereas the same
density and SSH started from rest ends at 47 % of the balanced surface KE
(-53 %), so the balance test is not vacuous.

Beta drift of the 40-day runs (surface-vorticity centroid, box coordinates, km):

| Case | start | day 40 | displacement |
|---|---|---|---|
| warm / anticyclone | (200, 200) | (177.1, 166.3) | dx = -22.9, dy = -33.7 (SW) |
| cold / cyclone | (200, 200) | (186.1, 225.3) | dx = -13.9, dy = +25.3 (NW) |
| warm, f-plane | (200, 200) | (200.0, 200.0) | 0 (exactly symmetric) |

Tilt of the warm vortex (T-anomaly centres; upper = 100 m layer; azimuth of the
upper part relative to the deep part, deg CCW from east):

| day | off. vs 700 m (km) | slope | azimuth | off. vs 1500 m (km) |
|---|---|---|---|---|
| 1 | 0.11 | 0.18 | 352 | 3.4 |
| 4 | 0.93 | 1.55 | 339 | 21.6 |
| 8 | 1.42 | 2.37 | 319 | 36.8 |
| 16 | 2.24 | 3.73 | 303 | 46.6 |
| 28 | 3.04 | 5.06 | 282 | 48.6 |
| 40 | 2.75 | 4.59 | 261 | 46.1 |

So the vortex does tilt: the 100 m centre leads the 700 m centre by up to ~3 km
(0.6 cell) over 600 m of depth, the offset direction rotating from east to south as the surface vortex drifts SW. The offset to 1500 m is large
(46 km) but that is the deep temperature anomaly (about 4 % of the surface
anomaly's weight) being advected/dispersed westward as a Rossby-wave wake, not a coherent
vortex core; read it as the signature of the deep response rather than a rigid
tilt of one eddy. Layers below ~3000 m carry < 1 % of the surface
anomaly (their centroids are noise-dominated).

Surface vorticity, warm vortex: peak `zeta/f0` -0.273 (day 1) -> -0.219 (day 16)
-> -0.295 (day 40); same-sign circulation within +-3 % to day 24, then +12 % at
day 40. The f-plane control shows the viscous/wave baseline: peak -0.273 ->
-0.208 (day 20) -> -0.248, circulation within 2 %. The cold vortex INTENSIFIES:
peak `zeta/f0` 0.246 -> 0.73 (day 20) -> 0.51 (day 40), same-sign circulation
x2.6, domain `En` +68 % (2.2e-4 -> 3.7e-4 m^2/s^2 at day 20; the warm vortex's
`En` decays 7 %). Not diagnosed here: it is a cyclone-anticyclone asymmetry in a
vortex with a surface Rossby number of ~0.25 and R/Rd ~ 2 and may include
baroclinic instability; a gradient-wind start (`vortex_balance = "gradient"`)
or a higher-resolution run would separate those.

## Caveats

- Cartesian beta-plane, flat bed, linear EOS only (refused otherwise).
- The 400 km box and `nu_h = 100` are sized for a 40-day drift of ~40 km;
  the closed walls start to matter for longer runs.
- Single-rank analysis script; merge multi-rank output first.
