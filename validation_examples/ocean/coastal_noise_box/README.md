# Coastal-noise box

A minutes-scale reproducer for the surface-vorticity symptoms of the 1/4-degree
Southern Ocean run (`southern_ocean_025_wind.nml`): extra relative vorticity at
coastlines and in the north sponge band.  It keeps that run's physics stack
(z_fixed 50-level tanh grid, Wright EOS, FV-MOM6 PGF, SADOURNY75_ENSTRO,
Smagorinsky + SMAGORINSKY_AH 0.06, quadratic drag, dt 900 s, OBC_SPONGE north
edge) on a 40 x 20 degree, 1/4-degree re-entrant channel with a staircase south
coast, a continental slope, an island and a 5-degree north sponge.

```bash
python3 make_inputs.py                      # bathy.nc + ic.nc (stdlib only)
CUDA_VISIBLE_DEVICES=0 ../../../build_cc70/rdb coastal_noise_box.nml   # ~2 min / 30 days, 1 V100
python3 metrics.py .                        # needs ncdump on PATH
```

`make_inputs.py --profile wall` gives vertical staircase walls over a flat
4000 m bed; `--profile rough` gives a 20 m coast and +-25 % grid-scale bathymetry
roughness; `--res 0.125` builds the 1/8-degree twin (set `nx=320, ny=160,
dx=dy=0.125`).

## Metrics (`metrics.py`)

Everything is computed on the dynamics' own circulation-form corner vorticity of
the surface layer (`k = nz`), rebuilt from the final restart's face velocities
and masked free-slip at land corners, exactly as `coriolis_adv` forms `q_corner`.
The centre-averaged `vorticity_z` diag is NOT used for the 2dx measure: its
four-corner average annihilates a pure `(-1)^i` or `(-1)^(i+j)` corner mode.

| field | meaning |
|---|---|
| `rms_coast` | wet corners with a land T-cell within 2 cells |
| `rms_int` | corners at least 8 cells from land, outside the sponge rows |
| `rms_spg` | the 24 northern rows (band + 1 degree), not coastal |
| `ratio`, `spg/int` | the two headline ratios |
| `nyq_*` | `<(d2x^2 + d2y^2)/2> / <zeta^2>`, `d2 = (z- - 2z + z+)/4`: 1 = pure 2dx, 0.375 = white noise, ~0 = smooth |
| `prof` | rms by distance to land (cells) |

## What the box shows (2026-09-30)

Numbers are day 30 unless marked; "base" is the tree before the two fixes below.

| case | ratio | nyq_coast | rms_spg | spg/int |
|---|---|---|---|---|
| straight-walled closed box, domain-edge walls | 5.8 | 0.01 | – | – |
| same, closed by a periodic land strip | 4.7 | 0.01 | – | – |
| channel, base | 8.8 | 0.16 | 6.6e-6 | 13.6 |
| channel, no sponge | 11.9 | 0.07 | 7.2e-7 | 1.5 |
| channel, sponge `relax_tracers = .false.` | 8.6 | 0.16 | 2.5e-7 | 0.5 |
| channel, sponge-reference fix | 9.2 | 0.15 | 2.5e-7 | 0.5 |
| vertical staircase walls, base, day 120 | 1.6 | 0.17 | | |
| vertical staircase walls, free-slip viscosity fix, day 120 | 1.0 | 0.11 | | |

* **Straight walls carry no 2dx signal** (nyq 0.01), whether the wall is a
  domain edge or an interior land strip.  The coastal excess there is the
  resolved boundary current.
* **The sponge band was a defect.**  Under `z_fixed` + `&ocean_zinit_nml` the
  layers were seeded sigma-style, the `target_source = "ic"` reference was
  snapshotted on those layers, and only then did the step-1 regrid move the
  state onto the z_fixed grid.  The sponge therefore relaxed each layer toward
  the zinit value hundreds of metres deeper, and to a different depth in every
  column of a different `H`, so it forced the band with a grid-scale,
  bathymetry-following density pattern (band surface 5-6 degC too cold, domain
  mean T -0.55 degC in 5 days; the 1/4-degree run shows the same -0.22 degC
  step).  With the seed on the z_fixed target the band's rms falls 26x.
* **The coastal excess on a staircase is mostly a resolved coastal jet**, not a
  grid mode.  Its width is physical: one cell at 1/4 degree, cells 2-4 at
  1/8 degree (~15-25 km either way), and nyq_coast stays at 0.04-0.2.  The
  staircase adds a coherent vorticity spike of the opposite sign at every
  concave step corner (Adcroft & Marshall 1998).  The velocity-form viscosity
  used to treat the zero stored at a land face as a Dirichlet-0 wall (partial
  no-slip), and the biharmonic operator rang against it.  With the free-slip
  mask this excess is gone on vertical walls (ratio 1.6 to 1.0 at day 120); on
  the sloping shelf it is unchanged.
* **Grid-scale bathymetry roughness gives interior 2dx** (`--profile rough`:
  nyq_int 0.32, the same as the 1/4-degree run's open ocean, 0.28-0.32).  The
  Coriolis form (energy vs enstrophy), the viscosity fix and dropping the
  biharmonic operator each leave it unchanged.
