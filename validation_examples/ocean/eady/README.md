# Eady baroclinic-instability benchmark

Linear baroclinic instability of a thermal-wind-balanced front, with
analytical growth rate from Eady (1949):

    σ_max ≈ 0.31 · f / √Ri      ;     Ri = N² / (∂u/∂z)²

A meridional T gradient `dT/dy < 0` (cold-to-north) is balanced by a
vertical shear `dU/dz` via thermal wind. A small random T perturbation
seeds the instability; the fastest-growing channel mode emerges from the
noise and grows exponentially at the analytical rate.

**Periodic channel.** The domain is a true reentrant channel — periodic in
x via `&ocean_bc_nml west/east="periodic"` (`nghost = 3` required), walls
north/south. The zonal wavenumber is quantised by the 250 km box
(λ = 250 / kx km) and the gravest cross-channel structure is
sin(πy/L_y), so the mode that grows is kx = 2 (125 km), whose channel
growth rate is σ = 1.98e-6 1/s (the unbounded-domain maximum is 2.49e-6
at 157 km, which the box cannot hold).

**Re-baselined 2026-09-13.** Until then the file shipped a 4× weaker front
(`eady_dT_dy = -5e-6`) with `nu_h = 100` and *decayed* — it had never
demonstrated the instability it exists to test. It now ships the front
the analytical scales describe (`-2e-5`), `nu_h = 20`, and a 60-day
window. The full evidence (a nu_h × dT_dy sweep, the modal spectra, and
why 60 days) is in the `eady.nml` header.

## Run

```bash
./build/rdb validation_examples/ocean/eady/eady.nml
```

50 × 50 × 10 cells on a 250 × 250 km flat-bottom periodic-x channel,
60 days @ dt = 600 s. Single V100 wallclock ~2 min.

Reference run (2026-09-13, this nml, V100): the kx = 2 mode grows at
σ = 1.93e-6 1/s (console `max|v|` fit over days 30–60; 1.99e-6 on the
NetCDF modal amplitude over days 15–25) against 1.98e-6 from theory;
`max|v|` 1.3e-4 → 0.2 m/s (~1500×); grid-scale (< 4Δx) content of v
< 0.1 % throughout; MaxCFL 0.017 → 0.041, no velocity clamping; Mass,
Salt, Heat budgets closed to ~1e-14 relative. The two WENO variants
(`eady_weno5.nml`, `eady_weno7.nml`: same basin + the windowed
tracer-advect drain) measure 1.93e-6 and 1.91e-6; the gfortran CPU build
(different `random_number` realisation of the seed) 1.87e-6, ×1060.
All under the default `mont` PGF; `form = "fv_lite"` reproduces every
printed digit on this flat-bed z* setup.

## What to look for

- **Day 0–1** — IG-wave shedding from the unbalanced IC (no balanced η,
  see "Known caveats"). Smag_AH absorbs the worst of it.
- **Day 1–10** — `max|v|` sits at the random-noise floor (~1e-4 m/s)
  while the kx = 2 mode grows underneath it; kx ≥ 5 decay.
- **Day 10–60** — clean exponential growth of `max|v|` at σ ≈ 1.9–2.0e-6
  (e-folding 5.8 days); from day ~25 the 125 km mode holds ~50 % of
  `max|v|`. One-and-a-half to two wavelengths visible in v, SSH and T at
  z = −500 m.
- **`En`** (total KE) is *not* the growth signal: it is 4.3e-3 m²/s² of
  basic-state jet that drifts *down* ~6 % (wall drain) until the eddies
  reach finite amplitude at day ~50, then 6.9e-3 at day 60.
- **Conservation** — mass, S, T means flat to ~1e-14 relative.
- **CFL** ≤ 0.05 throughout.

## Why the run stops at day 60

Past day ~65 the front collapses violently: `max|u|` passes 0.5 m/s at
day ~68, 1 m/s at day ~72, reaches the 2 m/s `maxvel` clamp and the
NaN-catch fires (~day 75) — with `n_inner` 20 *and* 30, with the clamp
removed (`maxvel = 10`), and with the old `nu_h = 100` stack. Halving
`dt` to 300 s survives to day 100 but with `En` reaching 0.21 m²/s²,
above the ~0.17 m²/s² of available potential energy the front holds:
energy is being manufactured, not released. The nonlinear collapse of
this front at 5 km / 10 layers / dt = 600 s is outside the validated
envelope and is an open item; this is a **linear-growth** benchmark by
construction.

## Known caveats

- **Resting-state growth at low viscosity (open, 2026-09-13).** With
  `eady_dT_dy = 0` — a motionless, stably stratified channel seeded with
  the same ±0.5 mK white-noise T perturbation — and `nu_h = 0`
  (Smagorinsky + `nu_4 = 5e8` still on), `En` grows from 0 to 3e-5 m²/s²
  in 25 days (`max|v|` 2 cm/s, e-folding 2.5 days, broadband in x,
  symmetric high vertical mode). A motionless stratified fluid has no
  energy source; that growth is numerical. It is what the earlier
  "`nu_h = 10` recovers Eady growth for the weak front" measurement was
  actually seeing (wavenumber-independent growth, faster for a weaker
  front). `nu_h` lowers its level, not its rate (`nu_h = 20`: 230× lower
  `En` at day 25, still growing at ~2e-6; `nu_h = 100`: `En` 1.6e-11, i.e.
  the old file's viscosity was what held it down). With the jet present the v
  spectrum is the kx = 2 mode alone (kx ≥ 5 decay at `nu_h = 20`, < 0.1 %
  grid-scale content to day 60), which is what makes the benchmark usable;
  any growth here must be checked for wavenumber selectivity before it is
  called Eady.
- The IC is hydrostatic-only (no balanced η); the first ~24 h is
  IG-wave shedding rather than linear growth. The fit windows above
  skip it.
- The RNG realisation of the noise seed differs between compilers
  (gfortran vs nvfortran `random_number`), which moves the day the mode
  clears the noise by a few days; the growth rate does not change.

## Tuning

Knobs in the namelist:

| Knob | Effect |
|---|---|
| `eady_dT_dy` ↑ (less negative) | Weaker front, slower growth — below ~1e-5 the mode no longer clears the noise seed in a useful time |
| `eady_dT_dz` ↑ | Stronger N², larger Ri, slower growth |
| `eady_pert_amp` ↑ | Stronger perturbation seed; the mode clears the noise earlier |
| `nu_h` | 10–50 leaves the mode's rate alone (1.99 → 1.89e-6); ≥ 100 drains the jet through the side walls (−8 % En / 25 d) and delays the mode ~10 days; ≤ 10 lets the short waves creep (resting-state growth, above) |
| `c_smag` ↓ | Less Smag_KH damping; sharper eddies |
| `smag_bi_const` ↓ | Less Smag_AH damping; rougher fields |
