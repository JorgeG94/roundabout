# Reading the sea-ice conservation ledger

For `&ocean_ice_nml scheme = "cice"` only. Default off; turn it on with
`&ocean_ice_nml ledger = .true.`. The module docstring of
`src/core/ice/coupling/rdb_ice_ledger_terms.F90` owns the definitions; this page
says how to read the output and how to add a term.

## The idea

The ice keeps three totals: mass `M = sum A (rhoi vice + rhos vsno)`, salt
`S = sum A rhoi vice Sref 1e-3` and enthalpy `E` (Icepack's reference: fresh liquid water at 0 degC, so ice enthalpy is negative). The step probes them, and
the flux accumulators `fresh`, `fsalt`, `fhocn`, `fswthru`, between its stages. A stage
is a NAMED TERM: `therm1`, `remap`, `newice`, `lateral`, `cleanup`, `dynamics`,
`transport`, `ridging`. For each stage and line:

```
water:  dM = (snow + evap + rain)_sources - dt sum A d(fresh)
salt:   dS = - dt sum A d(fsalt)
heat:   dE = (lw_in - snow_enth + atm_rad + atm_lat + frz_lat)_sources - dt sum A d(fhocn + fswthru)
```

All fluxes are Icepack's: `fresh`, `fsalt`, `fhocn` positive INTO the ocean. `evap > 0`
is condensation. `flat` is not a term (it cancels inside `fsurf - flat`). The residual of
a stage is `dstate - (sources - fluxes)`; its relative size is `|res|/scale`, `scale` being
the sum of the absolute values of every term of the line.

## The console line

At the status cadence, next to the ocean's block:

```
Ice ledger: 144 steps, |res|/scale water  1.60E-15 (lateral), salt  2.19E-15 (lateral), heat  1.01E-09 (therm1)
```

Each figure is the CUMULATIVE residual of that line, the largest over the stages,
relative to the cumulative scale; in brackets the stage that carries it. Expect mass and
salt at round-off (below `1e-11`; the stage named in brackets is just where the
round-off happens to be largest). Heat sits in `therm1` at about `1e-9` to `1e-8`: that
is the vertical solver's tolerance (`ICE_FERRMAX`, Bitz and Lipscomb 1999), never
round-off; everything else in the heat line is round-off too. A mass or salt residual above `1e-11`,
or a heat residual above the solver tolerance (or above `1e-11` outside `therm1`), in any
stage is a bug. The sums are EFP reproducing sums (Hallberg and Adcroft 2014), independent
of the summation order, so a restart reproduces the accumulators bitwise.

## What each stage books

| Stage | Moves the state | Books |
|---|---|---|
| `therm1` | growth, melt, snowfall, sublimation | the atmosphere exchange (`fswabs`, `flwout`, `fsens`, vapour), the ocean exchange, the snow and longwave sources |
| `remap` | nothing (ITD remap) | rain through the ice: a source and `fresh` of the same size |
| `newice` | frazil ice | the latent heat as a SINK of stored energy; `fresh`/`fsalt` with `update_ocn_f` |
| `lateral` | lateral melt | `fresh`, `fsalt`, `fhocn` |
| `cleanup` | rebin, zaps | the zap residual |
| `dynamics` | nothing | not probed: EVP changes velocities only; its effect is inside the `transport` probe |
| `transport` | advection | zero net on a closed domain |
| `ridging` | redistribution | the ridged snow the ocean receives |

**The new-ice sign.** New ice is the pure latent enthalpy `-rhoi Lfresh`; freezing a bank
of `E` joules lowers the stored energy by `E`. The `newice` heat term is negative and
equals minus the bank the ocean paid.

**What the ledger does not carry.** Water that leaves the ocean carries `cp T dM` against Icepack's
0 degC reference (`H_mass_conv`, about 2 % of the latent heat). The virtual-salt hand-off gives the ocean no
water, so its volume never changes and the term is zero by configuration; it comes back with the real-mass
coupling. `test_cice_closure` computes its size independently and reports it.

**Several ranks.** Every number above is the GLOBAL one on every rank (one EFP collective per ice step, in
`end_step`), bitwise the one-rank value on every `px x py` split; the console line is the same on every rank.
**Three named diagnostics.** `tab%clip_*` / `drop_*` (cumulative `ice_led_clip_*`): what the remap's guards drop
(booked as a sink of the `transport` stage, zero until a guard fires); `fout(:, cleanup)` (cumulative
`ice_led_zap_*`): the ITD clean-up zaps, mass and salt positive, heat negative; `tab%snow_reset` (cumulative
`ice_led_snow_reset`): the energy the thin-snow enthalpy reset adds inside `therm1`, REPORTED to attribute part of the
`therm1` heat residual, not subtracted.  `cum_zap` is the `cleanup` stage plus the post-transport clean-up (`transport`); the post-ridging clean-up's zaps are in `fout(:, ridging)` with the ridged snow lost to the ocean and are not accumulated.  With `advection = "upwind"` no clip record exists and the scheme's own sub-`puny` loss (at most `puny` times the thickness per cell) stays in the `transport` residual, unbooked.

**How big the thin-snow term is (reviewed estimate).** Per reset event the term adds at most about 70 J m-2 per kelvin
of snow cold (snow thinner than `hs_min = 1e-4 m`), about 700 J m-2 for snow at -10 degC. In the cell concerned that is
roughly 200 to 1300 times the solver's per-step bound, yet at most about 0.4 W m-2 spread over a step: physically negligible,
one-signed (it only ever adds energy), and it dominates the ledger's `therm1` heat verdict wherever cold snow is thinned below
`hs_min` (new ice diluting the snow of category 1, sublimation, transport or ridging carrying cold snow into thin-snow
cells, ice edges). `test_cice_closure` `snow_reset_reported` shows the term at 53 times the bound only because its forced snow
is very cold (`q = -1.5 rhos Lfresh`, about -79 degC); at -10 degC the same setup gives about 7 times. Whether the reset stays
is the maintainer's decision; the ledger only attributes it.

## Closing the coupled budget

`test_cice_closure` case `global_*` is the recipe. Heat: `dH_ocean + dE_ice` (from the ocean's
T/h arrays and the ice state) against the atmosphere's input to the ice (`flw`, `fsnow`, and the
ice's own `fswabs`, `flwout`, `fsens`, vapour), the open water's input rebuilt from the
BULK-FORCING outputs (not the assembler), the ice hand-over applied minus all the ice handed
over, the new-ice latent heat and the heat the ocean's own frazil clamp added. The residual is
tested against the vertical solver's tolerance as a BOUND; the ledger's own R is not
subtracted. Water: the ocean has no water budget under virtual salt; the ice mass change equals
the atmosphere's exchange minus the booked `fresh`, and the ocean's receipt of that `fresh`
is the salt dilution it shows in its own S/h arrays. Reported, not delivered: the rain enthalpy
on the ice-covered fraction, the precipitation that falls on the open water (enthalpy only in
the ocean today) and `H_mass_conv`.

## Adding a stage

1. Add its `ICE_LT_*` id (and name) in `rdb_ice_ledger_terms`; raise `ICE_LT_NSTAGE`.
2. Call `opts%ledger%probe(<id>, cice, flux)` after the stage in `ice_step_cice`.
3. If it has external sources, add them to `tl_end`.  Every value a probe keeps is an `efp_t` slot of the step vector (`ICE_LT_NV`, `ice_lt_ix_*`): a new reduced quantity needs a slot there, nothing else (the one collective carries it).
4. The restart tags (`ice_led_res_<line>_<stage>`) follow from the stage count; mind the
   256-entry restart registry.
5. Extend `test_cice_closure`.
