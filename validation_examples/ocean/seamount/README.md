# Idealised seamount — quiescent-IC σ-coord bathymetry test

Tier-1.5 bridge between the analytical-test suite and real bathymetry
(Tasman 2 km). A Gaussian seamount over a closed flat-bottom basin
under **uniform ρ, zero velocity, no forcing**.

The answer is analytical and trivial: **zero motion forever**. Any
non-trivial flow that develops is a bug in BPG, Coriolis-adv,
continuity, or vmix interacting with a sloping h_layer field.

This is the configuration that caught three latent ocean bugs in
one May-21 session (PP81 wall discontinuity, FV-PGF z_centre sign,
ghost-cell bathymetry fill) — bugs that flat-bottom tests had missed
for months.

## Run

```bash
./build/rdb validation_examples/ocean/seamount/seamount.nml
```

100 × 100 km basin (50 × 50 cells @ 2 km), 15 σ-layers,
peak depth 200 m, basin depth 4000 m (20× depth ratio).
f-plane at f = 0. No wind, no surface flux, no sponges.

Single V100 wallclock: seconds.

## What to look for

- **max|u|, max|v| stay at machine zero** throughout the integration.
- **max|η| / max|h - H_init| stay at machine zero**.
- **CFL stays at 0** — no flow → no CFL.

Any drift > ~1e-10 m/s indicates an unbalanced operator. Typical
signatures:

- BPG sign or quadrature bug → growing |u| concentrated on the
  seamount flank.
- Ghost-row h_layer not filled → spurious face fluxes at the walls.
- Smag closure dividing by tiny h_layer → blowup near the peak.

## Variants in this dir

| File | What it changes |
|---|---|
| `seamount.nml` | Canonical setup (above). |
| `seamount_auto.nml` | Same geometry, auto-derived `n_inner`. |
| `seamount_flat.nml` | Flat-bottom control (sanity check that the dyn-core is quiet absent the seamount). |
| `seamount_bench_full.nml` | Production-style: linear T(z), f≠0, KPP+PP81, ALE remap on zstar_sigma, Laplacian ν_h, sponges, diag manager. Bridge to Tasman 2 km without the GEBCO real-bathy confound. |
