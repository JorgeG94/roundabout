# Flat-bottom eddy validation

Canonical "do we get eddies?" test for the Roundabout ocean dyn-core,
modelled after the MOM6 `ocean_only/double_gyre` example.

## Why this setup gets eddies

1. **Stratified initial condition.** Linear T(z) from 20°C at the
   surface to 5°C at the bed (over 4000 m depth) provides the
   available potential energy that baroclinic instability extracts.
   *Uniform T runs cannot produce baroclinic eddies* regardless of
   how the wind drives them — there's no stored energy to convert.
2. **Sinusoidal wind pattern** drives a subtropical anticyclone +
   subpolar cyclone (Sverdrup balance), with western intensification
   into a swift jet (~1 m/s at the WBC). The WBC + the lateral T
   gradient it generates is the unstable mean state.
3. **β-plane Coriolis** localises the WBC against the western wall
   (it doesn't form anywhere else on an f-plane).
4. **Flat bottom** removes the steep-bathy / sigma-BPG-error
   pathology that plagued the Tasman setup, so the integration
   stays stable long enough (~30 days) for eddies to spin up.

## Run

```bash
./build/rdb validation_examples/ocean/eddy_test/eddy_test.nml
```

Expected wallclock single V100: **~10-15 min for 30 days**
(300 × 300 × 15 cells, dt_outer = 300 s, 8640 outer steps × ~30 ms
each).

Output: `out_eddy_test/eddy_test_rank_000000.nc` — **180 frames per
variable** (every 4 hours over 30 days, set by `ocean_diag_dt_out =
14400.0`). Default variables: SSH, temperature, salinity, u, v, KE,
h_layer.

## What to look for

**Days 0-7** — spin-up. Wind drives Ekman in the surface layer +
Sverdrup transport interior. Western boundary current builds. KE
grows but stays smooth (no closed circulation patterns yet).

**Days 7-15** — WBC reaches mature speed (~0.8-1.2 m/s near the
western wall, in a band ~50 km wide). Still laminar.

**Days 15-25** — onset of barotropic / baroclinic instability of the
WBC. Wave-like meanders appear in the jet path. KE starts to
develop spatial structure beyond the wind forcing pattern.

**Days 25-30+** — eddy shedding. Closed circular features detach
from the WBC and propagate westward (anticyclones) or interact with
the gyre interior. **This is what "eddies" looks like.**

## Quick check

```python
import xarray as xr
ds = xr.open_dataset('out_eddy_test/eddy_test_rank_000000.nc',
                     decode_times=False)
ke_surf = ds['KE'].isel(z_KE=-1, time_KE=-1).values
print(f'Day 30 surface KE: peak {ke_surf.max():.3f}, mean {ke_surf.mean():.4f}')
# Eddies: peak KE >> mean KE (factor of ~10-50 in a developed eddy field)
# Pure laminar Ekman: peak/mean ~ 2-3
```

Or with matplotlib:
```python
import matplotlib.pyplot as plt
ssh_30 = ds['SSH'].isel(time_SSH=-1).values.squeeze()
plt.pcolormesh(ssh_30, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
plt.colorbar(label='SSH (m)')
plt.title('SSH at day 30')
# Eddies → multiple alternating high/low spots; gyre alone → a single
# basin-scale dipole.
```

## Tuning knobs (in `eddy_test.nml`)

| Knob | Effect |
|---|---|
| `ocean_nu_h` ↓ | Sharper WBC, eddies spin up faster (but watch for blowup) |
| `taux_magnitude` ↑ | Stronger WBC, faster instability onset |
| `T_init_surface - T_init_bottom` ↑ | Stronger stratification, faster baroclinic growth |
| `t_end` ↑ | More mature eddy field |
| `dx` ↓ (e.g. 2.5 km) | Resolves smaller eddies; ~4× wallclock |

## Known limitation: wind-driven 2-layer runs go runaway

A 2-layer cheat config (named `eddy_test_2layer.nml`, removed from
the tree as documented-failing) blew up at day 9-15 with both
vcoord variants attempted:

* `z_fixed` vcoord: interface anchored, can't tilt → BT mode
  absorbs all wind input → SSH ±12m → Inf.
* `isopycnal` vcoord: interface tilts freely → baroclinic mode
  runs away → u peaks 3+ m/s, h_layer ±30% → KE > 5 m²/s² → Inf.

The structural cause is that our explicit Smag biharmonic damping
is CFL-capped, while a wind-driven 2-layer setup at this domain
size pumps in energy faster than the cap allows us to dissipate.
TwoLayerSW.jl handles this via an AM3 corrector (implicit damping
in the time scheme) + Shapiro filter every step (continuous 2Δx
mode killing on `h`).  We have neither — our dissipation stack is
calibrated for dissipated regional ocean (real bathy + drag + KPP
+ thermo), not idealised energetic wind-driven 2-layer setups.

For a baroclinic-instability movie, **use `../eady/eady.nml`
instead** — IC is the unstable jet, no forcing, saturates at modest
amplitude, validated 25 days stable.  50×50 basin is small but
shows real eddies.

A Shapiro-filter port (the high-leverage low-effort item from the
session-5 agent comparison) would unlock the wind-driven 2-layer runs on this domain.  Not yet
prioritised.

## Status

First attempt at the eddy validation — calibrated against MOM6
defaults but not yet compared to a reference run.  Things to
check after each run:

1. Does the time loop survive 30 days? (`exit=0` from the binary)
2. Does KE develop spatial structure beyond the forcing pattern?
   (NetCDF SSH / vorticity slices should show closed circulations)
3. Does mass conservation stay at FP throughout? (the `Mass`
   column in `[stats]` lines stays constant)
4. Does surface relative vorticity develop ±0.1-0.5 Rossby-number
   patches by day 20? (the textbook eddy signature)
