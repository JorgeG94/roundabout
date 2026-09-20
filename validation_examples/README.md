# Validation examples

Canonical / analytical-solution benchmarks for the Roundabout dyn-core.
Each subdir documents a setup that has a known reference behaviour
(analytical solution, published benchmark, or canonical regression
config) — runnable end-to-end with the in-tree binary and small
enough to spin up on a single GPU in minutes.

## Layout

```
validation_examples/
  ocean/      Tier-1+ canonical configs for the ocean dyn-core
    seamount/                  quiescent IC over a Gaussian seamount
    geostrophic_adjustment/    Rossby SSH-bump adjustment
    eady/                      baroclinic instability of a thermal-wind jet
    eddy_test/                 wind-driven double-gyre with stratification
    double_gyre/               MOM6 ocean_only/double_gyre reproduction
    ice_shelf_cavity/          resting cavities under a flat / sloping ice lid
    ... and the rest of ocean/ (acc_channel, baroclinic_channel, bc_inst,
        epbl_mld, ideal_age, island_at_rest, neverworld2, sea_ice, sponge_demo,
        tides, …) — each with its own README
  test_cases/ bathymetry-prep helpers for the Tasman regional validation
```

## Ocean benchmarks (in suggested run-up order)

| Subdir | What it tests | Reference |
| ------ | ------------- | --------- |
| `ocean/seamount/` | BPG over varying-h_layer bathymetry under quiescent uniform-ρ; should give zero motion forever. | Quiescent-IC validation; caught three latent ocean bugs in one May-21 session. |
| `ocean/geostrophic_adjustment/` | Coriolis coupling, IG-wave radiation at √(gH), 2π/f oscillation at the bump centre, energy conservation after the transient. | Rossby (1937); textbook problem with analytical scales. |
| `ocean/eady/` | Linear baroclinic-instability growth rate of a thermal-wind-balanced front; saturation to a turbulent eddy field. | Eady (1949); analytical σ ≈ 0.31 · f / √Ri. |
| `ocean/eddy_test/` | "Do we get eddies?" — stratified, wind-driven double-gyre on a β-plane. WBC instability + barotropic/baroclinic eddies. | Canonical MOM6 `ocean_only/double_gyre`. |
| `ocean/double_gyre/` | Closer MOM6 reproduction with `driver_run_ocean` + diag manager. Two entries (`double_gyre.nml`, `double_gyre_mom6.nml`). | MOM6 reference run at `~/nci/projects/access-nri/cpu_MOM6/ocean_only/double_gyre`. |
| `ocean/ice_shelf_cavity/` | A stratified ice-shelf cavity at rest must STAY at rest. Three cases differing by one ingredient each isolate the load bookkeeping (flat lid ⇒ exactly zero), the datum (uniform density ⇒ 1.3e-20) and the sigma-coordinate pressure-gradient truncation under a tilted ice base (the measured Phase-6 baseline). | Losch (2008); ISOMIP+ geometry/EOS from Asay-Davis et al. (2016) *GMD* **9**, 2471–2497; the correction target is Yung, Hallberg, Adcroft & Morrison (2026) *JAMES* **18**, e2025MS005645. |

## Running anything here

From the repo root:

```bash
./build/rdb validation_examples/ocean/<subdir>/<file>.nml
```

Each subdir has its own README with expected wallclock, what to
look for in the output, and the relevant analytical scales / tuning
knobs.

## Adding a new benchmark

Goes here if it has a published / textbook / analytical reference,
or is a canonical regression target. Real-bathymetry / case-study
configs go in [`examples/`](../examples/) instead.
