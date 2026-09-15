# Double-gyre examples

Wind-driven double-gyre on the ocean C-grid dyn-core: the canonical
Munk/Stommel picture (Sverdrup interior, westward-intensified boundary
current, two counter-rotating gyres separated by a jet).

These namelists run through the production driver
(`app/main.F90 → driver_run → driver_run_ocean`), the stress-tested
path.  Run them with the main `rdb` binary:

```
./build/rdb validation_examples/ocean/double_gyre/double_gyre_mom6.nml
```

> For a **NetCDF-free** ocean throughput run (e.g. profiling, or building
> on Intel / AMD with `-DRDB_ENABLE_NETCDF=OFF`), use the standalone
> ocean benchmark `benchmarks/bench_ocean` — it reads the same namelists
> through the same setup path but skips the I/O.  See `benchmarks/`.

## `double_gyre_mom6.nml`

Reproduces NOAA-GFDL/MOM6-examples/ocean_only/double_gyre as closely as
the current Roundabout ocean stack allows: 44×40×2 spoon basin, gprime
reduced-gravity PGF, Sadourny+HK Coriolis, linear bottom drag, z*-full
vertical coord.  Stable to day 580 with the production envelope.  See
the header comments in the namelist for the full knob-by-knob mapping to
MOM6 and the known caveats (gprime vs Wright EOS, BOUND_CORIOLIS, HBBL
drag footprint, spherical→Cartesian projection).

Output: per-rank NetCDF (`out_double_gyre_mom6/double_gyre_rank_*.nc`)
with SSH, T, S, u, v, KE per model day via the diag manager.
`animate_double_gyre_3panel.py [NETCDF]` renders it.

## `double_gyre_linear_nk10.nml`

The same basin grown to **NK=10** using the linear density-range IC
(`&ocean_ic_nml rho_lightest = 1035.0, rho_range = 2.0`) — the MOM6
`COORD_CONFIG="linear"` + `DENSITY_RANGE` analogue.  Instead of a
hand-listed `layer_rho_init`, the layer densities self-populate
linearly, so you scale the column by editing `nz_layers` alone.  Because
`gprime` is a 2-layer reduced-gravity form, this NK>2 config uses the
z-corrected FV PGF (`form = "fv_lite"`) reading the per-layer
`rho_layer` the knob fills, on a distributed `zstar` coordinate.
Illustrative "it scales" demo, not a tuned validation target.

Quick smoke: set `t_end = 30.0` for a ~16 s spin-up run on a single V100.

### Comparison target

The MOM6 reference run lives at
`~/nci/projects/access-nri/cpu_MOM6/ocean_only/double_gyre/`:

- `prog__0001_006.nc` — instantaneous prognostic snapshots
- `ave_prog__0001_003.nc` — time-mean snapshots
- `ocean_geometry.nc` — grid + bathymetry (re-project onto our Cartesian grid)
