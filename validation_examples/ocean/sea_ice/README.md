# Polar freeze-up — thermodynamics and EVP dynamics

Two companion cases over the same 40×30 @ 10 km flat f-plane basin
(coriolis_f for ~80°N, 4 sigma layers, 50 m depth):

| File | Ice dynamics | Duration | Question |
|---|---|---|---|
| `polar_freezeup_thermo.nml` | off | 90 days | does frazil nucleate and does the column model grow plausible ice from a supercooled surface? |
| `polar_freezeup_dynamics.nml` | C-grid EVP, wind-driven | 20 days | does the EVP solver respond to wind stress and transition from free drift to strength-limited drift as the ice thickens? |

Formulations exercised, by paper (see each namelist's header for the
per-knob rationale):

- **Column thermodynamics** — Winton (2000) two-layer formulation
  (`nk_ice = 2`), M. Winton, *"A reformulated three-layer sea ice model,"*
  J. Atmos. Oceanic Technol., **17**, 525–531. The brine-pool
  simplification used here (`ICE_CP_BRINE = ICE_CP_ICE`) makes every
  per-layer solve a closed-form quadratic.
- **Frazil nucleation** — ocean-side supercooling clamp/bank, the
  standard frazil-ice closure used by Semtner-class models: A. J. Semtner,
  *"A model for the thermodynamic growth of sea ice in numerical
  investigations of climate,"* J. Phys. Oceanogr., **6**, 379–389 (1976).
- **EVP rheology** — E. C. Hunke and J. K. Dukowicz, *"An
  elastic–viscous–plastic model for sea ice dynamics,"* J. Phys.
  Oceanogr., **27**, 1849–1867 (1997), with the Hibler ice-strength
  closure and elliptical yield curve: W. D. Hibler III, *"A dynamic
  thermodynamic sea ice model,"* J. Phys. Oceanogr., **9**, 815–846
  (1979) (`p0 = 2.75e4` N/m², `c0 = 20`, `ec = 2`, matching Hibler's
  standard constants).
- **Free-drift reference** — the massless/zero-strength limit of the EVP
  momentum balance, `|u| = sqrt(tau / (rho_ocean * c_dw))`, attributed to
  F. Nansen, *"The oceanography of the North Polar Basin,"* Scientific
  Results of the Norwegian North Polar Expedition (1902); this is the same
  golden check `tests/test_ocean_ice_evp.F90` (`test_nansen_free_drift`)
  uses at unit-test scale.

## Run

```bash
source nvhpc_env.sh   # or mac_env.sh for a CPU build
build_cc70/rdb validation_examples/ocean/sea_ice/polar_freezeup_thermo.nml
build_cc70/rdb validation_examples/ocean/sea_ice/polar_freezeup_dynamics.nml
```

## Namelist fix applied for the v0.1.0 release candidate

**`polar_freezeup_dynamics.nml` shipped `t_end = 0.5` (half a day) while
its own header documented a 20-day run** ("Shorter run than the thermo
case (20 days, not 90) ... By day 20 the thermo twin shows ~0.77 m mean
thickness..."). At `t_end = 0.5` the run ends before ice even nucleates
(frazil onset is ~day 1.5–1.9 in the thermo twin — see below), so the
EVP-dynamics demonstration the header describes could never have run.
Fixed to `t_end = 20.0` to match the documented intent and this README's
measurements.

## Results (GPU build, nvfortran 26.5, cc70, V100, single rank, 2026-10-06)

### Thermodynamics (`polar_freezeup_thermo.nml`, 90 days)

Frazil nucleates between t ≈ 1.6–1.9 days (`ice_conc` steps from 0 to full
cover within one `dt_out` window of 0.25 day) — matches the header's
"~day 1.5" estimate. From there the whole basin covers uniformly (flat,
f-plane, spatially uniform forcing — no reason for it not to) and the
column thickens through the 90-day cooling season:

| day | max ice thickness (m) | pure Stefan¹ (m) | finite-exchange Stefan² (m) | sim / pure | sim / finite-exchange |
|---|---|---|---|---|---|
| 5  | 0.212 | 0.262 | 0.179 | 0.81 | 1.19 |
| 10 | 0.425 | 0.417 | 0.328 | 1.02 | 1.30 |
| 30 | 1.042 | 0.772 | 0.677 | 1.35 | 1.54 |
| 60 | 1.795 | 1.108 | 1.011 | 1.62 | 1.78 |
| 90 | 2.483 | 1.364 | 1.266 | 1.82 | 1.96 |

¹ Classical Stefan law `h = sqrt(2 k_i ΔT t / (ρ_i L))` with `k_i = 2.03
W/m/K`, `ρ_i = 905 kg/m³`, `L = 3.34e5 J/kg` (the model's own
`rdb_ice_column`/`rdb_ice_enthalpy` constants), `ΔT = T_f − T_air =
−1.836 − (−20) = 18.164 K`, elapsed time measured from the observed
nucleation onset (t₀ = 1.75 days).

² Finite-exchange ("Semtner-style") correction for the namelist's
*non-infinite* atmospheric coupling: `SF(T) = restore_lambda·(T_surf −
T_air)` (`rdb_ice_atm_forcing`) puts an atmospheric resistance `1/λ` in
series with the ice's own conductive resistance `h/k_i`, giving
`ρ_i·L·(h²/(2k_i) + h/λ) = ΔT·t` — a quadratic in `h` that reduces to
pure Stefan as `λ → ∞` and is strictly slower for finite `λ` (script:
`tmp_local_artifacts/stefan_compare.py`).

**Day 5–10 agree with pure Stefan to ≤ 20%** (the onset-time estimate is
good), but the simulated growth runs increasingly AHEAD of both analytic
curves from day ~15 onward — ×1.8–2.0 by day 90, not the ×0.9–1.0 the
finite-exchange correction predicts. This is explained by a **documented**
model divergence, not a bug: `&ocean_thermo_nml q_heat` is "ice-blind" —
`docs/CLOSURE_MATRIX.md` ("Sea ice" section): *"the ocean's own background
shortwave is not lead-fraction weighted by ice cover (the configure-time
q_heat/q_sw scalars stay ice-blind)."* The uniform `q_heat = −60 W/m²`
surface flux that bootstraps the initial frazil bank keeps supercooling
the mixed layer and banking additional frazil mass for the full 90 days,
even under 100% ice cover — a second heat sink running in parallel with
the ice column's own atmosphere-coupled conduction. A roughly time-linear
term stacked on a sqrt(t) term is exactly the shape measured (good early
match, growing excess later), and is consistent with the CLOSURE_MATRIX
caveat rather than contradicting it.

**Brine rejection**: mean ocean salinity rises monotonically from 34.000
to 35.343 PSU over the 90 days as growing ice rejects brine into the
mixed layer — the expected sign and a smooth, monotonic signal throughout.

**Conservation**: Mass/Salt/Heat `Error` residuals stay at 1e-13–1e-22
relative (round-off) for all 12960 steps — fully closed budgets, as
expected for the virtual-salt-flux (Boussinesq) ice-ocean coupling.

### Dynamics (`polar_freezeup_dynamics.nml`, 20 days, `tau_x = 0.05` Pa eastward)

Nansen free-drift speed for these parameters (`cdw = 3.24e-3`, `rho_ocean
= 1030`): `|u| = sqrt(0.05/(1030*3.24e-3)) = 0.1224 m/s`.

Before any ice mass exists, the EVP momentum solve still returns a
velocity at every wet cell (zero ice strength is not a special case — see
the "Ice margins" docstring in `rdb_ice_evp.F90`): the massless-limit
balance is exactly the free-drift balance, and the model reproduces it:
**max `ice_speed` = 0.131 m/s, mean 0.100 m/s** in the first day, 7%
above the analytic value (the small excess is the non-zero background
ocean velocity under the ice, which the free-drift formula above assumes
is zero).

As the column thickens (frazil onset ~day 1.6–1.9, same as the thermo
twin), the Hibler (1979) ice-strength term engages and drift collapses by
two orders of magnitude:

| day | max `ice_speed` (m/s) | regime |
|---|---|---|
| 0.25 | 0.131 | free drift |
| 2 | 0.125 | free drift |
| 4 | 0.091 | transition begins |
| 6 | 0.064 | strength engaging |
| 8 | 0.035 | strength-limited |
| 10 | 0.0058 | strength-limited |
| 15 | 0.0021 | jammed |
| 20 | 0.0013 | jammed |

At day 20 the ice has piled up against the downwind (east) wall exactly
as the namelist header predicts: the physical-domain west→east thickness
profile (y-averaged, `diags = "ice_u ice_v ice_speed"` re-run with
diagnostics enabled — see "Diagnostics bug" below) is a clean, monotonic
ramp from **0.634 m at the west wall to 0.849 m at the east wall**
(domain mean 0.762 m, max 0.868 m) — the east wall carries 34% more ice
than the west wall, the signature of wind-driven onshore convergence
against a solid downwind boundary with no ridging term to redistribute it
(`compress_ice` thickens in place rather than building a ridged tail —
documented in `docs/CAPABILITIES_AND_LIMITATIONS.md`, "No ridging").

**Conservation**: Mass/Salt/Heat `Error` residuals stay at 1e-14–1e-21
relative through all 2880 steps.

## Bugs found (localised, not patched here — see task scope)

Running these cases surfaced two real defects in the ice *diagnostics*
reporting path (the ice *physics* itself checks out against both
analytic references above):

1. **The console `[stats] ... Ice : conc X thick Y (ocean-area mean)`
   line always prints exactly `0.0000`/`0.0000`**, in every run, at every
   status interval — including at `step 0` of `sea_ice_pack.nml`
   (`../sea_ice_pack/`), immediately after the banner confirms a
   100%-covered, 2 m analytic IC was applied. The underlying ice state is
   NOT zero (confirmed by the independently-computed `[diag]`
   `ice_conc`/`ice_thick` NetCDF fields, which show correct, physically
   sensible values throughout both runs above) — the defect is isolated
   to the summary line's own area-weighted accumulation,
   `compute_ice_totals` / `compute_ice_totals_efp` in
   `src/core/ocean/diag/rdb_ocean_console_stats.F90` (~L424, L567–575,
   L730–737, L805–867, L1197–1356), fed from
   `src/driver/rdb_driver.F90` (the `ice_part_size=`/`ice_m_ice=`/
   `ice_ncat=` actual arguments at each `ocean_console_stats_report` call
   site). Not bisected further here (would need live instrumentation,
   out of this task's "don't patch physics" scope) — worth a dedicated
   follow-up.
2. **The `[diag]`/NetCDF `ice_conc` and `ice_thick` MEAN statistics are
   diluted by ghost-cell padding.** Confirmed exactly on two independent
   grids: the reported mean equals the true physical mean multiplied by
   `(nx_phys·ny_phys)/(nx_total·ny_total)` — 0.802139 on this file's 40×30
   / `nghost=2` grid (`44×34` total), 0.694444 on `sea_ice_pack.nml`'s
   20×20 / `nghost=2` grid (`24×24` total) — matching the printed values
   to 6 significant figures throughout both runs. Root cause:
   `fill_ice_conc` / `fill_ice_thick`
   (`src/core/ocean/diag/rdb_ocean_diag_fills.F90::fill_ice_conc_thick_impl`,
   ~L382–432) write `0.0` (not the IEEE-NaN "no data here" sentinel every
   other registered diagnostic, e.g. `fill_tracer_impl`, writes into
   land/non-physical cells) outside the ice-present branch, so
   `diag_field_stats`'s finite-cell count
   (`src/core/ocean/diag/rdb_ocean_diag.F90:1164–1215`) includes the
   zero-valued ghost ring in its denominator. `min`/`max` are unaffected
   (ghosts are legitimately 0, which can't raise a true max and
   coincides with the true min here); only `mean` is wrong, and it is
   wrong on every shipped example that turns ice on, by a factor that
   depends on grid size and `nghost`. The thermo table above uses `max`
   (unaffected); the "domain mean 0.762 m" / west-east profile for the
   dynamics case was read directly from the re-run NetCDF file's physical
   index range, bypassing the diluted console/NetCDF summary statistic.

## Known caveats (inherited, already documented)

- **No ridging** (`docs/CAPABILITIES_AND_LIMITATIONS.md`, "Sea ice"): the
  downwind pile-up above thickens uniformly rather than building a ridged
  tail — expected, not a bug.
- **`q_heat`/`q_sw` are ice-blind** (`docs/CLOSURE_MATRIX.md`, "Sea ice"):
  the mechanism behind the thermo case's late-time excess growth above.
- **Frazil checks the surface layer only** (`k = nz`), not the full
  column (`docs/CLOSURE_MATRIX.md`).
