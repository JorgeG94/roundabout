# Neverworld2 idealized basin

A Roundabout port of **Neverworld2** (Marques et al. 2022, GMD, "Neverworld2:
an idealised model hierarchy to study ocean mesoscale eddies across
resolutions"; MOM6-inspired). A single Pangaea-style ocean basin on a
60° lon × 140° lat spherical sector (lon 0→60 E, lat 70 S→70 N) with a
**re-entrant (periodic-x) southern channel** — the Drake-Passage analog
that lets an ACC-like circumpolar jet form.

## Geometry + forcing

- **Bathymetry** (`topo_config="neverworld2"`,
  `set_bathymetry_neverworld2`): the fractional-depth formula in
  normalized `(x,y)∈[0,1]²` coordinates carves
  - the **great northern wall** + **Antarctica** (south wall) via
    `1.1·spike` terms, the southern wall breached by the Drake gap;
  - two meridional **continents** — South America at the x-edges (which
    straddle the periodic seam, so the *northern* seam sits inside land)
    and Africa mid-basin — that wall the northern wind-driven gyres;
  - the **Drake/Scotia ridge** system (cosbell bumps) setting the
    channel sill;
  - small-amplitude **roughness** (two cosines).

  `nl_continent_amp` scales the continent/ridge block (1.0 = full
  continents, **0.0 = aquaplanet with the southern channel only**, a
  good first bring-up); `nl_roughness_amp` scales the roughness. There
  is no upper depth clamp (matching MOM6 — only `D = max(D,0)`), so the
  roughness can lift the depth a few % above `max_depth`.

- **Wind** (`wind_config="neverworld2"`,
  `set_wind_stress_neverworld2`): the canonical 3-band zonal stress
  τ_x(lat) — southern westerlies / trades / polar-easterlies, peak
  `taux_magnitude` (0.2 Pa default), τ_y ≡ 0.

Both formulas are re-derived from the paper (MOM6-inspired) and verified
cell-by-cell against an independent reference in
`tests/test_ocean_neverworld2.F90`.

## What v1 covers (and what it doesn't)

This is the **geometry bring-up**: temperature-only linear EOS,
`zstar` + ALE remap, KPP boundary layer, the production lateral
(`nu_h=10000` + Smagorinsky) + linear-HBBL bottom-drag envelope. Forcing
is **wind-only** (`q_heat=0`) so the stratified IC spins up wind-driven
gyres + an ACC-like channel jet without spurious uniform convection.

Documented **fast-follows** (see `local_archive/specs/neverworld2_study.md`):
the latitude-banded **SST-restoring** buoyancy forcing (the paper's
smooth `T*(lat)` + piston velocity — the only genuinely new kernel), the
15-layer **exponential-thermocline IC**, finer resolutions (¼°, 1/16°),
and the `nw2_tracers` diagnostic ideal tracers. Roundabout runs `zstar` +
ALE rather than MOM6's isopycnal/Lagrangian layers — the right design
analog, not bit-for-bit.

## Run

```bash
cd validation_examples/ocean/neverworld2
# NVHPC + NetCDF toolchain loaded in this shell (module / spack / conda)
CUDA_VISIBLE_DEVICES=<gpu> ../../../build/rdb neverworld2.nml
```

1° resolution → 60 × 140 × 12 cells over H = 4000 m, `nghost = 3`
(periodic-x needs ≥3 for the PPM + biharmonic seam stencil). The
70 S/70 N extent stays well under the spherical |lat| ≤ 90 pole guard
even with ghosts. The shipped nml runs 2 days as a stability check; set
`t_end` higher for a spin-up.

Expected: no blowup; a coherent eastward jet through the southern
channel; wind-driven gyres in the enclosed northern sub-basins. The
aquaplanet variant (`nl_continent_amp=0`) is the simplest sanity run.
