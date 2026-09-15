# Ocean tidal validation

Two runnable cases that exercise the ocean tidal forcing end-to-end (not just
the astronomy goldens): the **equilibrium body tide** (C1 / C2) and the
**OBC boundary tide** (C3). Both drive a single M2 constituent and are checked
by fitting a pure M2 sinusoid (period **12.4206 h**) to the SSH output — if the
fit explains essentially all the variance (R² → 1), the modelled signal *is* M2.

| Case | nml | Forcing | Grid | Validates |
|------|-----|---------|------|-----------|
| Body tide | `body_tide_basin.nml` | equilibrium potential (`&ocean_tides_nml`) | **spherical** closed basin | C1 (+C2 SAL, optional) |
| OBC tide | `seamount_tidal_obc.nml` | M2 prescribed at the west open boundary | cartesian seamount | C3 |

## Run

```bash
# from a scratch working dir (output lands in ./out_<name>/)
rdb validation_examples/ocean/tides/body_tide_basin.nml
rdb validation_examples/ocean/tides/seamount_tidal_obc.nml

# check the M2 signal
python validation_examples/ocean/tides/analyze_m2.py out_body_tide_basin/body_tide_basin_rank_000000.nc --spinup-days 3
python validation_examples/ocean/tides/analyze_m2.py out_seamount_tidal_obc/seamount_tidal_obc_rank_000000.nc
```

Both run serially on a CPU build in well under a minute (`RDB_ENABLE_GPU=OFF`,
`RDB_ENABLE_NETCDF=ON`).

## Expected results (single V100/CPU, reference)

**`body_tide_basin`** — 10-day run, drop 3-day spin-up:
- basin centre: M2 amp ≈ **0.8 cm** (the equilibrium tide), R² ≈ **0.99**
- wall (resonant) probe: M2 amp ≈ **28 cm**, R² ≈ **0.98**
- closed basin ⇒ mass/salt/heat conserve to round-off (~1e-15)

**`seamount_tidal_obc`** — 4-day run, drop 1.5-day spin-up:
- every probe R² ≈ **0.999–1.000** (an essentially pure M2 response)
- amp ≈ 40 cm (west) → 50 cm (seamount) → 64 cm (east); west→east **phase lag**
  ≈ −35° is the tide propagating across the domain
- open boundaries ⇒ the domain mass *breathes* with the tide: the reported
  ~1.2e-4 "mass error" is exactly the SSH signal (≈50 cm / 4000 m), **not** a
  conservation bug
- the interior west→east transect is the validation; the SE corner (Flather
  open-east × solid wall) shows an amplified but still-pure-M2 spike (R²≈1.0) —
  a known finite-difference artifact of that boundary junction, not a C3 defect

## Config notes (load-bearing)

- **The body tide is fail-loud on cartesian grids** — the astronomical potential
  needs real lat/lon, so `body_tide_basin` uses `&ocean_grid_nml
  grid_config="spherical"`. Its amplitude is physical/built-in; you only choose
  `constituents` (a whitespace/comma string, e.g. `"M2 S2 K1 O1"`).
- **The ocean C-grid reads `&ocean_bc_nml` only** — the coastal `&boundary_nml`
  block does *not* reach the ocean path.
- **OBC tidal constituents are keyed by angular frequency (rad/s), not name.**
  M2 = `1.4051890e-4`; amplitude in metres, phase in radians. With
  `obc_tidal_nodal=.true.` the C3 astronomical/nodal correction is folded in
  (then a parseable `&ocean_tides_nml ref_date` is required and the frequency
  must match the catalog within 1e-4 or it aborts).

## Variations

- **Exercise C2 SAL:** in `body_tide_basin.nml` set `use_sal=.true.`,
  `beta_sal=0.09` — the effective gravity becomes `−g(1−β)∇η` and the response
  amplitude shifts accordingly.
- **Nodal correction:** set `add_nodal=.true.` (body) / `obc_tidal_nodal=.true.`
  (OBC) to fold in the 18.6-yr nodal `f/u` factors.
- **Full resolution:** both nmls are shrunk for CPU turnaround; scale
  `nx/ny/nz_layers` up for a GPU run.
