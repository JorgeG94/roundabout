# Baroclinic instability channel

The clean, fast, analytically‑anchored **eddy‑generation** test — the one the
suite was missing. A thermal‑wind front in a periodic‑x channel goes
baroclinically unstable and rolls up into a field of coherent mesoscale
eddies. Two configs:

| file | layers | grid | cells | runtime (V100) | use |
|---|---|---|---|---|---|
| `baroclinic_2layer.nml` | 2 | 192×96 @2 km | 37 k | ~3 min / 40 d | **fast regression** |
| `baroclinic_15layer.nml` | 15 | 384×192 @2 km | 1.1 M | ~21 min / 50 d | depth‑resolved demo |

## The physics

The Eady IC (`ic_config="eady"`) sets a meridional buoyancy gradient
`dT/dy < 0` balanced by a vertical shear `dU/dz` via thermal wind. The front
is **baroclinically unstable**: a tiny seed perturbation grows exponentially,
extracting the mean flow's available potential energy, then saturates by
overturning the isopycnals and pinching off into discrete vortices — the
textbook **linear growth → roll‑up → saturated turbulence** life cycle.

- **2 layers = the Phillips problem.** The instability lives on the *interior*
  PV jump between the two layers, so eddies fill the channel interior (unlike
  the continuous‑Eady mode, which is trapped at the top/bottom boundaries).
- Deformation radius `Ld = √(g'H_eff)/f ≈ 14 km` (2‑layer) / `NH/(πf) ≈ 18 km`
  (15‑layer), resolved at ~7–9 cells on the 2 km grid — fine enough that a
  **physical eddy is unambiguously distinguishable from grid noise** (the key
  reason this works where a marginal 5 km double‑gyre does not).

## Expected signature (the regression check)

Read it straight off the console stats — no NetCDF needed:

1. **`En`**: flat at the mean‑front baseline → **exponential growth** (onset
   ~day 18–22 for 2‑layer, ~16–20 for 15‑layer) → **saturation**
   (~0.09 for 2‑layer, ~0.42 for 15‑layer). *Unforced*, so it slowly **decays**
   after the peak as the front is consumed.
2. **Conservation**: Mass / Salt / Heat to ~1e‑13 relative throughout.
3. **`MaxCFL` bounded** (~0.3). A *monotonic growth to a CFL trip* would be a
   **numerical** instability, not physical eddies — that's the failure mode to
   catch.

If `En` grows and saturates with bounded CFL and exact conservation, the
dyn‑core's baroclinic instability is healthy.

## Why this example exists (lessons baked in)

- The flat‑bottom double‑gyre (`../eddy_test`) is only validated through its
  **laminar** spin‑up; at 5 km its "eddies" are **marginally resolved**
  (Ld ≈ 6 cells) and indistinguishable from grid noise. This channel resolves
  Ld well and uses the *interior* instability, so the test is clean.
- `bebt=0.2` (barotropic velocity projection) + the **2‑D gravity‑wave CFL**
  (`cfl_bt_safety` with the corrected `auto_n_inner`) keep the fast mode stable
  — without them the barotropic split can seed a spurious growth that *looks*
  like eddies. Don't drop them.
- Diag is kept **lean via coarse cadence** (12‑hourly): dumping every field ×
  every layer at a fine cadence costs ~46 % of wall time (measured). The
  dominant I/O lever is the frame count, so coarsen `dt_out`; the vorticity
  movie only needs the surface (top) layer, which the animation reads.

## Run + animate

```bash
# NVHPC + NetCDF toolchain loaded in this shell (module / spack / conda)
CUDA_VISIBLE_DEVICES=<gpu> ./build/rdb \
   validation_examples/ocean/baroclinic_channel/baroclinic_2layer.nml
python validation_examples/ocean/baroclinic_channel/animate_baroclinic.py \
   out_bci_2layer/bci_2layer_rank_000000.nc baroclinic.mp4
```

The animation renders surface relative vorticity `ζ/f` — coherent red/blue
vortices and filaments = a healthy baroclinic eddy field.

## Deferred

Surface buoyancy **restoring** (relax T toward `T*(y)` with a piston velocity)
would maintain the front and give a *statistically steady* eddy field instead
of the spin‑down — that's the Neverworld2 SST‑restore fast‑follow, not yet on
the tree.
