# Sea-ice analytic IC — pack under wind, from step 1

The first shipped namelist that turns sea ice on **without** growing it
from frazil first: `&ocean_ice_ic_nml conc_config = "uniform"` seeds a
2 m, 100%-covered pack over the whole (20×20 @ 10 km, flat, f-plane)
domain before the first step, pushed through the exact enthalpy
inversion so the seeded state is a provable fixed point of the ITD
restore on step 1 (see the namelist header). The same steady eastward
wind stress as `../sea_ice/polar_freezeup_dynamics.nml` (`tau_x = 0.05`
Pa) is applied from `t = 0`, so — unlike that case — there is a live pack
to push on immediately, with no frazil spin-up needed. Short run (2
days): the question is "does a pack that already exists respond to wind
from step 1," not long-horizon thermodynamics.

Formulations and citations are the same as `../sea_ice/README.md`
(Winton 2000 column thermodynamics; Hunke & Dukowicz 1997 EVP rheology;
Hibler 1979 ice strength; Nansen 1902 free-drift reference).

## Run

```bash
source nvhpc_env.sh
build_cc70/rdb validation_examples/ocean/sea_ice_pack/sea_ice_pack.nml
```

## Namelist fix applied for the v0.1.0 release candidate

**`sea_ice_pack.nml` shipped `dynamics = .true.` without
`transport = .true.`.** `./rdb --validate-only` warns on exactly this
combination: *"drift-only mode has no CFL backstop on u_ice (transport's
positivity abort is the only ice-velocity guard) and is unvalidated."*
`transport = .true.` costs nothing the example doesn't already satisfy
(`ncat = 5 > 1`, all edges wall, single rank) and is what actually lets
the EVP-computed velocity field move ice mass around — without it,
`dynamics` alone only computes `u_ice`/`v_ice`, never redistributes
`m_ice`/`part_size`. Added; the warning is gone and this README's
"mass response" measurement below depends on it.

## Results (GPU build, nvfortran 26.5, cc70, V100, single rank, 2026-10-06)

A **uniform, fully-covered** pack starts at maximum Hibler ice strength
everywhere (`pres_mice = (p0/rho_ice)*exp(-c0*max(1-ci,0))` saturates at
`ci = 1` from the first step — there is no concentration gradient
anywhere in the domain for it to converge against except the walls), so
the expected EVP behaviour is the opposite of the ice-free-start dynamics
case: **near-immobile, stress-limited drift from step 0**, not free drift
decaying into it.

Measured (re-run with `&ocean_diag_nml enabled = .true., diags = "ice_u
ice_v ice_speed"` to get the field, since the console summary line is
broken — see "Bugs found" below):

| t (day) | max `ice_speed` (m/s) | ice_thick range (m) |
|---|---|---|
| 0.25 | 1.47e-4 | 1.992 – 1.999 |
| 1.00 | 1.47e-4 | 1.994 – 2.002 |
| 2.00 | 1.46e-4 | 1.994 – 2.009 |

Max speed is **~900× below the Nansen free-drift speed** for the same
wind (0.1224 m/s, see `../sea_ice/README.md`) and **barely changes** over
the 2 days — a rigid, fully-jammed pack pressed against the (wall-only)
boundary, exactly the "stress-limited convergence, no grid-scale noise"
signature Hunke & Dukowicz (1997) EVP is designed to produce, and the
mirror image of the free-drift-then-jam curve measured in the
`polar_freezeup_dynamics.nml` case. Thickness stays within 0.5% of the
seeded 2.0 m everywhere (physical-index range of the re-run NetCDF file
— the ghost-diluted console/file mean reads 0.694444 of this, see below),
confirming the IC is close to an exact fixed point of both the ITD
restore and the EVP, as the namelist header claims, with only a small
residual thermodynamic adjustment (the IC's −4 °C / 4 PSU ice column is
not in exact thermal equilibrium with the `q_heat = −60 W/m²` / `air_temp
= −20 °C` forcing).

**Conservation**: Mass/Salt/Heat `Error` residuals stay at 1e-14–1e-22
relative through all 288 steps.

## Bugs found (localised, not patched here — see task scope)

Same two diagnostics-reporting defects as `../sea_ice/README.md`,
reproduced independently on this file's different grid size (useful
cross-check):

1. **Console `[stats] Ice : conc 0.0000 thick 0.0000` is wrong from
   `step 0`** — printed immediately after the run banner confirms the
   uniform 2 m / `conc = 1.0` IC was applied, i.e. before any physics
   could plausibly zero it. See `../sea_ice/README.md` for the
   localisation (`rdb_ocean_console_stats.F90`
   `compute_ice_totals`/`compute_ice_totals_efp`, fed from
   `rdb_driver.F90`).
2. **`[diag]`/NetCDF `ice_conc` mean is diluted by ghost padding**: this
   file's 20×20 physical / `nghost=2` grid (24×24 total) reads exactly
   `mean = 0.694444 = (20·20)/(24·24)` once the pack is at its true
   `conc = 1.0` everywhere — confirms the `(nx_phys·ny_phys)/
   (nx_total·ny_total)` dilution factor measured on the sibling
   40×30/`nghost=2` grid (0.802139) is grid-size-dependent, as the root
   cause (`fill_ice_conc_thick_impl` not NaN-sentinelling ghost cells)
   predicts. See `../sea_ice/README.md` for the full localisation.

## Known caveats (inherited, already documented)

- **No ridging** (`docs/CAPABILITIES_AND_LIMITATIONS.md`): irrelevant to
  this specific 2-day run (nothing converges enough to ridge) but shared
  by the whole `compress_ice` family.
- **No freshwater/mass coupling**: the 2 m ice column applies no dynamic
  surface loading on the ocean beneath it (documented, `docs/
  CLOSURE_MATRIX.md` "Sea ice").
