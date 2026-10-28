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

## Ordering inside a thermo interval: `&ocean_vmix_nml diabatic_first`

Only matters at `dt_therm_ratio > 1`.  Both settings give the same answer at
`dt_therm_ratio = 1`.

| Work | `diabatic_first = .true.` (default) | `.false.` (ACCESS-OM3 `DIABATIC_FIRST = False`) |
|---|---|---|
| surface heat/salt fluxes, shortwave, restoring, geothermal, ideal-age source, vertical tracer diffusion + KPP nonlocal | first step of the interval (`1, N+1, ...`), inside the stage chain | last step (`N, 2N, ...`), AFTER that step's tracer-window drain (and for ratio 1 after the stage loop), BEFORE the remap: they act on the advected tracers |
| mixing coefficients (EPBL / KPP / kappa-shear / tidal / convective kt, ks, kv) | thermo step, stage chain | thermo step (N), stage chain, i.e. from the pre-drain state |
| lateral tracer diffusion (hdiff), Redi coefficients + apply | thermo step, stage chain | last step (`N`), after the window drain and a tracer-ghost refresh, once with the full tracer dt (Redi subcycles on that dt); the Redi coefficients are built from the post-advection T/S/h |
| vertical tracer advection | thermo step | unreachable (Eulerian-z is ssp_rk2-only) |
| ice thermo block, ice transport, flux assembly (engine, post-step) | after the last step of the interval, applied by the next interval's first step | after step `N-1`, applied by step `N` |
| Fox-Kemper MLE transports | thermo step | thermo step (reads the MLD of the previous EPBL call) |
| wave speed, isopycnal slopes, VarMix, MEKE (GM coefficient chain) | thermo step | first step of the interval (`1, N+1, ...`): thickness diffusion before the dynamics, OM3 `THICKNESSDIFFUSE_FIRST = True` |
| ALE remap | first step (no windowed advection) or last step (`dt_tracer_advect_ratio > 1`) | last step, after the drain that closes the last window |
| horizontal tracer window drain | steps `R, 2R, ...` | steps `R, 2R, ...` |

`diabatic_first = .false.` needs `split_scheme = "pred_corr"` and no
`&ocean_cavity_melt_nml` (both fail-loud at configure when `dt_therm_ratio > 1`).
Known limit: the diffusivities and the EPBL/KPP boundary-layer solve are
computed in the stage chain from the pre-drain state; only the tracer-changing
applies and the lateral tracer diffusion are deferred.  The tidal and
law-of-the-wall BBL `kd_int` are checkpointed (`tidal_kd_int`, `bbl_kd_int`),
since under `.false.` a resume always falls between their refreshes.  With
`&ocean_hdiff_nml kappa_h > 0` the single full-dt explicit step is checked at
configure (warning above half the diffusive-number limit, error above it).
