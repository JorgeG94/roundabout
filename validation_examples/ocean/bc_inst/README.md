# bc_inst — Two-Layer Baroclinic Instability


## The experiment

Two-layer, reduced-gravity, **isopycnal** shallow water on a spherical
re-entrant channel:

- lon ∈ [-193.75°, -171.25°] (span 22.5°), lat ∈ [53.625°, 64.875°], central
  latitude φ₀ = 59.25°, N×N (256 or 512), R = 6.371×10⁶ m, Ω = 7.2921×10⁻⁵.
- Periodic-x (re-entrant channel), wall-y (free-slip).
- Upper (surface) layer H₁ = 500 m, ρ₁ = 1025; lower (bed) H₂ = 1500 m,
  ρ₂ = 1027; g = 9.81, g' = 0.019141 m/s²; **R_d = 21.4 km** (~5 cells at
  256², ~10 at 512²).

The IC (`ic_config="baroclinic_jet"`) is a geostrophically balanced tanh jet
in the surface layer plus a front-localised sech² meander seed (3 zonal
wavelengths). Initial `v ≡ 0`, so **any growth of the console `En` (KE/mass)
or max|v| is the baroclinic instability**. Growth is exponential in the linear
phase (e-folding ~a few days for the tuned jet).


## Cases

| namelist | N | jet L | Δξ | dt | purpose |
|---|---|---|---|---|---|
| `bc_inst_tuned_512.nml` | 512 | 40 km | 200 m | 427 s | **headline perf case** (u_peak ~0.76 m/s, Ro ~0.15) |

## Run

```
# build with GPU on the NVHPC toolchain (single-GPU perf comparison)
module load nvhpc && cmake -B build -S . -DRDB_ENABLE_GPU=ON && cmake --build build -j
cd validation_examples/ocean/bc_inst
CUDA_VISIBLE_DEVICES=<gpu> ../../../build/rdb bc_inst_tuned_512.nml
```
