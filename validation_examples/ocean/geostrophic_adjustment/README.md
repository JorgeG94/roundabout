# Geostrophic adjustment (Rossby's classic problem)

Drop a Gaussian SSH bump on an f-plane and watch the system shed
inertia-gravity waves and settle (partially) into geostrophic balance.

Cheap diagnostic of:

- Coriolis coupling (oscillation period at the centre = 2π/f).
- Hydrostatic + barotropic gravity-wave physics (IG-wave front speed = √(gH)).
- Energy conservation (bounded KE + PE after the IG transient).
- PGF / continuity sign + scaling.

## Run

```bash
./build/rdb validation_examples/ocean/geostrophic_adjustment/geostrophic_adjustment.nml
```

1000 × 1000 km closed basin (100 × 100 cells @ 10 km), single
barotropic layer, depth H = 1000 m, f-plane at f = 10⁻⁴ s⁻¹.
200 km Gaussian bump, amp = 1 m.

Single V100 wallclock: seconds.

## Analytical scales

| Quantity | Value |
|---|---|
| Deformation radius `Rd = √(gH) / f` | 990 km |
| IG-wave speed `c = √(gH)` | 99 m/s |
| Oscillation period at centre `2π/f` | ~17.5 h |
| Bump-to-Rd ratio `L / Rd` | 0.20 (most energy radiates) |

## What to look for

- **Outgoing IG-wave fronts** propagate at c = √(gH). Time-snapshot
  η fields should show ring patterns expanding at this speed.
- **η at the bump centre oscillates at 2π/f**, damping toward the
  residual balanced value (small here, since L ≪ Rd).
- **Total energy** (KE + PE relative to background) is bounded —
  it doesn't drift up or down after the transient settles.

Tuning `ga_length_scale` to L ≫ Rd (e.g. 5000 km) switches the
behaviour: most energy is retained as a balanced anticyclonic ring
around the bump.

## Diagnostic interpretation

| Symptom | Likely cause |
|---|---|
| Bump amplitude grows in time | Wrong PGF sign, broken continuity |
| Oscillation period ≠ 2π/f | Wrong Coriolis coupling |
| Energy drifts after transient | Spurious sources/sinks in barotropic step or Coriolis-adv |
