# Tracer-advection cadence: `dt_tracer_advect_ratio` / `dt_therm_ratio`

How (and when) to advect tracers less often than the dynamics step — the
ocean path's equivalent of MOM6's `DT_TRACER_ADVECT` / `DT_THERM`.

## The two knobs (`&ocean_vmix_nml`)

| knob | default | meaning |
|---|---|---|
| `dt_tracer_advect_ratio` | `1` | horizontal tracer advection fires every `N` outer steps over the *accumulated* face transport (MOM6 `DT_TRACER_ADVECT`) |
| `dt_therm_ratio` | `1` | thermo step (ALE remap, EOS, surface fluxes) fires every `N` outer steps with effective `dt = N·dt_dyn` (MOM6 `DT_THERM`) |

Constraint (fail-loud at configure): `dt_therm_ratio` must be an integer
multiple of `dt_tracer_advect_ratio`, so every thermo/ALE step is also a
tracer-advect step (no accumulation window is ever cut by a remap).

## How tracer advection works at each setting

- **`dt_tracer_advect_ratio = 1` (default):** the fused continuity+tracer
  kernel advects tracers *in lockstep with the thickness, every RK2 stage*
  (twice per outer step), on the stage-instantaneous transport. Cheapest
  per step; tracer↔thickness consistency to RK2 truncation order.
- **`dt_tracer_advect_ratio = N > 1`:** during the dynamics the thickness
  advances while the tracer's CONCENTRATION `T = hTr/h_layer` is held fixed
  (`hTr` is re-weighted onto the new `h` each stage — MOM6's prognostic IS a
  concentration, so it gets this for free; freezing the CONTENT instead
  corrupts every `T = hTr/h` consumer, the EOS above all), and the layer mass
  transport is accumulated (`uhtr += ½·mass_flux·dt` per stage). Once per
  `N`-step window the accumulated *time-mean* transport advects every tracer
  in a single conservative swept-PPM **drain** (CFL sub-cycled),
  reconstructing `hprev = areaT·h_end + div(uhtr)` so transport telescopes
  *exactly* onto the thickness the dynamics produced — exact constancy,
  monotone, conservative (Adcroft & Hallberg 2006 / MOM6 `advect_tracer`).

  The drain and both halves of that concentration hold fill the
  `*_budget_horiz_adv` accumulators, so the console's CLOSED salt/heat budget
  stays valid at `ratio > 1` — including mid-window reports. (Before
  2026-09-12 nothing filled them, the console fell back to raw drift, and a
  run with a surface heat flux reported the flux itself as a ~5e-5 "leak".)

## Which to use

**Default is `1` deliberately — do not raise it globally.** A coarse cadence
trades accuracy for speed and is **not safe for every regime**: configs with
sharp density fronts (lock-exchange) or tight surface-flux budgets
(cooling + KPP) go *unstable* / break their budget closure when thermo and
advection only fire every few steps. This is exactly why MOM6 makes
`DT_THERM` a per-run tuned parameter, never a fixed `>1` default.

**For large, smooth production runs, set `ratio > 1` per config** — this is
where the windowed drain pays off, on both speed *and* consistency:

```
&ocean_vmix_nml
   dt_therm_ratio          = 4
   dt_tracer_advect_ratio  = 4   ! must divide dt_therm_ratio
/
```

Measured on the seamount bench (512×256×30, single V100): `ratio = 4` runs
**~24% faster** than the every-step default — the tracer drain fires 4× less
often (the swept-PPM scheme is heavier per call, so the win is the *cadence*,
not the scheme; at `ratio = 1` the heavier drain has no amortization and is
*slower* than the fused path, which is why `1` keeps the fast fused path).

Rule of thumb: raise the ratio until a validation run drifts beyond your
accuracy tolerance, then back off — per config, MOM6-style. Reduced-gravity
/ passive-tracer runs (where S/T do not feed the dynamics) can also just set
`do_horizontal_advection = .false.` on those tracers to skip the work
entirely.

## See also

- `continuity_tracer_drain` (`rdb_continuity.F90`) — the windowed swept-PPM
  drain. The sub-cycle budget is sized to the actual Courant
  (`ceiling(2·cfl_max)`, ≥1), since the limiter moves ≥ 0.5·CFL per pass.
- Known divergence from MOM6 (tracked): the drain falls back to PCM in the
  2 cells nearest a true wall, where MOM6 keeps full PPM (periodic seams are
  unaffected). See the `TODO(MOM6-fidelity)` in `drain_parabola_x`.
