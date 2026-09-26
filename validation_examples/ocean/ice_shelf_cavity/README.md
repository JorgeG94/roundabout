# Ice-shelf cavity at rest — the Phase-5 gate and the Phase-6 baseline

Three quiescent ice-shelf cavities. Each starts motionless, has **no energy
source of any kind** (no wind, no surface flux, no melt, no bottom drag, no
vertical mixing, no lateral viscosity), and must therefore stay motionless.
Every joule of kinetic energy any of them develops was manufactured by the
discretisation — and because the three differ by exactly one ingredient each,
the set says *which part of the discretisation made it*.

| file | lid | density | what it isolates |
|---|---|---|---|
| `cavity_flat_lid_rest.nml` | flat, 300 m | ISOMIP+ COLD stratification | the load bookkeeping + the datum. Every interface gap `Δe(K) = 0` ⇒ structurally bit-zero |
| `cavity_uniform_rho_rest.nml` | sloping + calving front | **uniform** | the load bookkeeping alone. `N² = 0` ⇒ every trapezoid error `G(K) = 0`, however steep the lid |
| `cavity_sloping_lid_rest.nml` | sloping + calving front | ISOMIP+ COLD stratification | **the headline**: the sigma-coordinate pressure-gradient truncation under a tilted ice base |

All three share the geometry (48 × 6 × 15 @ dx = dy = 2 km, flat bed 720 m),
the ISOMIP+ linear EOS, the f-plane at 75 °S, `VCOORD_SIGMA`,
`&ocean_pgf_nml form="fv_mom6"` + `p_top_in_bc=.true.`, and the analytic
geopotential initial condition. None pins `split_scheme`, so the stability
suite builds an `__ssp_rk2` twin of each automatically.

## How to run

```bash
./build/rdb validation_examples/ocean/ice_shelf_cavity/cavity_sloping_lid_rest.nml
```

~36 s for the 30-day run on one CPU core (gfortran Release). The three
namelists, and their three scheme twins, are rows in
`tests/regression/stability_manifest.py`:

```bash
python3 tests/regression/stability.py --tier 2 --build-dir build_gcc \
        --cases cavity_flat_lid_rest,cavity_uniform_rho_rest,cavity_sloping_lid_rest
```

## The initial condition is the subtle part

Under `VCOORD_SIGMA` the layer interfaces **tilt with the ice base**. The
shipped analytic stratifications — `&tracer_nml T_init_surface` /
`T_init_bottom` and their `S_init_*` twins — lay a profile out linearly in
**layer index**, i.e. *along the tilted coordinate*. Under a sloping lid that
tilts the isopycnals with the lid, which is a real available potential energy
and a real baroclinic adjustment: the run would measure physics, not
truncation, and would be wrong about both.

These cases therefore use `&ocean_zinit_nml source = "linear"`, which
evaluates `T(z)`/`S(z)` at each layer centre's **true geopotential depth**
(the column top sits at `z = −z_draft`, not at `z = 0`). Isopycnals are then
flat in `z` on every column, sloping lid included. That path — the draft
offset in `build_z_ctr` plus the analytic source — is what P5.3 added; before
it, `&ocean_cavity_dyn_nml` and `&ocean_zinit_nml` refused each other.

### …and the column top has to be trimmed to the load

The load `p_ice_ref = ρ₀·g·z_draft` is the displaced weight at the
**reference** density. The stratified ISOMIP+ COLD water it displaces is
lighter (its surface is 0.29 kg/m³ below `ρ₀`), by
`g·I(z_draft)`, `I = ∫_{−z_draft}^0 (ρ − ρ₀) dz = −0.2878·z_d + 4.197e-4·z_d²`
kg/m². At `η = 0` every level under the ice therefore sits `−g·I` off the open
ocean at the same `z` — a **depth-uniform** bottom-pressure gradient
`(g/ρ₀)(ρ(−z_draft) − ρ₀)·∇z_draft`, up to `1.8e-5 m/s²` on the sloping lid
(zero at `z_draft = 343 m`, where `ρ = ρ₀`), some 5000× the truncation below.
The legacy split discarded it with the depth mean; the MOM6 barotropic split
(`&ocean_bt_nml bc_pgf_forcing`, the default since 2026-09-25) hands it to
the barotropic mode, which adjusted to it on day 1 (En `1.061E-06`, over the
gate).

The sloping-lid case therefore sets `&ocean_cavity_dyn_nml
trim_ic_for_p_surf = .true.` — MOM6's `TRIM_IC_FOR_P_SURF` (`trim_for_ice`):
the **load is kept** (the ice mass is what is prescribed) and each loaded
column's initial top moves to the depth `s` where the water above weighs it,
`g·∫_{−s}^0 ρ dz = p_ice_ref`, an initial `η = z_draft − s` of −3.0 mm at the
front to −4.80 cm at `z_draft = 343 m`, with `T`/`S` evaluated at the trimmed
layer centres. The interface pressures then equal the open ocean's exactly
(a linear density is integrated exactly by the layer-midpoint stack), and the
depth-mean face force falls from `1.6e-5` to `3.9e-8 m/s²` on this geometry
(`test_ocean_cavity_load::trim_ic_balances_the_depth_mean_pfu`). The flat-lid
and uniform-ρ cases need no trim: a uniform draft has no gradient to balance,
and `T_ref`/`S_ref` on the uniform `T`/`S` make `ρ ≡ ρ₀` there.

## What the residual should be, from the algebra

Flat bed + sigma makes the vertical gap between two columns' `K`-th
interfaces `Δe(K) = σ_K·D`, with `D` the draft step across the face. The
FV-MOM6 Pass-3 numerator telescopes into a trapezoid error
`G(K) = −(Δe(K)³/12)·ρ₀N²`, and what survives the split solver replacing the
column mean with the barotropic solution is

```
PFu(k) − ⟨PFu⟩_h = N²·D³·(3σ_k² − 1) / (12·dx·H̄)
```

— **second order in `dx`, cubic in the draft slope, independent of `nz`**,
peaking at the ice base (`σ = 1`) at `a_peak = N²D³/(6·dx·H̄)`. Here
`N² = 8.01e-6 s⁻²`, `D = 13.8 m`, `dx = 2000 m`, `H̄ = 470 m`:

```
a_peak = 3.73e-9 m/s²    U = a_peak/|f| = 2.65e-5 m/s    En = ½U² = 3.5e-10
```

Under the MOM6 split the depth **mean** reaches the barotropic mode too, and
carries one more truncation term: the FV-MOM6 in-layer integral treats each
layer's density as uniform (short by `g·ρ_z·h³/12`), and sigma layers
`h = W/nz` differ across a face, so the mean carries a depth-uniform
`≈ N²·W·ΔW/(4·nz²·dx)` = 1–4e-8 m/s² — balanced by a sub-millimetre surface
tilt in this 2-D geometry, with a seiche of `|u| ~ a·L/(2c) ~ 2e-5 m/s`, inside
the plateau below.

## Measured (gfortran 15.1 Release, single rank, 2026-09-25)

`En` is the domain-mean specific kinetic energy off the daily `[stats]` line;
`|u|_rms = sqrt(2·En)`. Default MOM6 barotropic split; the sloping lid with
the trimmed IC.

| case | scheme | day 10 | day 30 | day 60 | `|u|_rms` at day 30 |
|---|---|---|---|---|---|
| flat lid | pred_corr | `0.000E+00` | `0.000E+00` | — | 0 |
| flat lid | ssp_rk2 | `0.000E+00` | `0.000E+00` | — | 0 |
| uniform ρ | pred_corr | `1.742E-21` | `2.017E-20` | — | 2.0e-10 m/s |
| uniform ρ | ssp_rk2 | `7.484E-21` | `1.303E-19` | — | 5.1e-10 m/s |
| **sloping lid** | **pred_corr** | **`1.877E-09`** | **`2.298E-08`** | `3.492E-06` | **2.14e-4 m/s** |
| sloping lid | ssp_rk2 | `1.875E-09` | `1.850E-04` | — | 1.9e-2 m/s |

Budget residuals stay at 8e-13 relative in every run; `MaxCFL ≤ 1e-3` at day
30; no CFL truncation and no velocity clamping anywhere.

Controls for the sloping lid (pred_corr, same tree): **untrimmed**, the
barotropic adjustment to the load shortfall reads `1.061E-06` (d1),
`1.020E-08` (d10), `1.316E-07` (d30); the **legacy split**
(`bc_pgf_forcing = .false.`, untrimmed), which discards the depth mean, reads
`1.626E-09` (d10), `3.043E-08` (d30), `3.894E-06` (d60) — the table this file
carried through 2026-09-24, within 30 % of the trimmed default at every daily
sample. The substitution table and the cross-toolchain check below were
measured on the legacy split.

**Cross-toolchain** (legacy split, 2026-09-20). The same three cases on the
GPU build (nvfortran 26.5, `-stdpar=gpu`, `cc70`, one Tesla V100-DGXS, 4320
steps in 74.6 s) at day 30:

| case | gfortran CPU | nvfortran V100 |
|---|---|---|
| flat lid | `0.000E+00` | `0.000E+00` |
| uniform ρ | `1.340E-20` | `3.168E-20` |
| **sloping lid** | **`3.043E-08`** | **`3.043E-08`** |

The sloping case agrees to every printed digit after 4320 steps on two
different compilers and two different execution models, which says the
measurement is the algorithm and not an FMA accident. The two machine-zero
cases differ at the ulp-accumulation level, as they must — `0` and `1e-20`
are the same statement.

### Read the sloping case in two parts

**(a) The plateau, days 1–20 — the measurement.** En sits at `0.9–1.9e-09`
(`|u|_rms = 4.2–6.2e-05 m/s`) against the derived `3.5e-10` / `2.6e-05 m/s`.
Two to three times, for a peak-acceleration estimate compared against a domain
rms, is as close as that formula can be asked to come. The two controls pin it
to the right term: the flat lid is **exactly** `0.000E+00` for 30 days (every
`Δe(K) = 0`) and the uniform-density twin is `2.0e-20` — twelve decades
below — so the residual is the `ρ₀N²Δe³` trapezoid error and not a
mis-cancelled 5.26 MPa ice load wearing its clothes.

**(b) The growth from day ~22 — the finding.** En leaves the plateau on a
~3.2-day e-folding and then **saturates**: `3.5e-06` by day 60 (2.6 mm/s rms)
with the rate visibly decaying, budgets still exact. It is bounded, not
runaway. Characterised by substitution (legacy split, 2026-09-20):

| substitution | result | conclusion |
|---|---|---|
| `dt` 600 → 300 s | day 30: `3.069E-08` vs `3.043E-08` | **not** an `(ω·dt)ⁿ` outer-split mode |
| ISOMIP+ Laplacian `nu_h = 6.0` (their Table 4) | peak `1.13E-09` (d2), min `9.6E-11` (d20), `9.65E-10` (d30), `1.24E-06` (d60) | only DELAYS it ~10 days |
| biharmonic `nu_4 = 1.0e7` (the dx⁴-scaled `resting_stratified_channel` value) | `3.58E-10` (d30), `6.69E-07` (d60) | also only delays — the mode is not purely grid-scale |
| flat lid | `0.000E+00` | the mode needs the sloping lid |
| uniform density | `1.34E-20` | the mode needs the stratification |

So the sloping-lid PGF truncation is the **source**, and no defensible
dissipation removes it — it only postpones it. `ssp_rk2` then amplifies the
same internal-wave field catastrophically from day 12 (`1.85E-04` at day 30,
8000× the default), which is what EXPERIMENTAL on that scheme means, stated as
a number; the same separation `validation_examples/ocean/eady/resting_stratified_channel.nml`
records for a lid-free channel.

### The gate, and what it deliberately does not do

`cavity_sloping_lid_rest` gates at `en_rest_max = REST_1MM_S` (`0.5e-06`
m²/s², 1 mm/s) — the bar the sloped resting family already uses, and the
velocity a resting sub-shelf cavity has no excuse to exceed. The measured
day-30 peak clears it by 22× in energy, 4.7× in velocity. At tier 1 (30
days) it carries a **scoped XFAIL on `energy:rest-settles` only**, because (b)
makes the final sample the peak at every horizon past day ~22 and short of
saturation; the 20-day tier-2 twin ends on the plateau (88 % of its day-10
peak) and passes both gates. Do not "fix" a failure here by turning the trim or
`bc_pgf_forcing` off — both put the load shortfall back.

It does **not** clear the tighter `REST_100UM_S` (`0.5e-08`) the Phase-5
design proposed: day 30 is 4.6× over it. That is recorded rather than legislated
away — restoring `REST_100UM_S` is the Phase-6 acceptance criterion. Do not
close either by widening the bar, by shortening the run, or by putting
viscosity into these files.

The flat-lid and uniform-density cases gate at `REST_1UM_S` (`0.5e-12`), the
"must be bit-zero" bar, with no XFAIL on either scheme.

## Against the literature — this is the Phase-6 baseline

The paper that fixes this class of error is

> Yung, C. K., Hallberg, R. W., Adcroft, A., & Morrison, A. K. (2026).
> *Assessment of a finite volume discretization of the horizontal pressure
> gradient force beneath sloping ice shelves.* **Journal of Advances in
> Modeling Earth Systems, 18**, e2025MS005645.
> https://doi.org/10.1029/2025MS005645

Their reported numbers, verbatim:

- Abstract — *"spurious velocities reduced to order 10⁻⁹ m s⁻¹ or smaller in
  quiet, linear stratification test cases."*
- §5.1.2, σ-coordinate icemount, linear stratification — *"inclusion of the
  nonlinear pressure reconstruction in the icemount test improves spurious
  currents from order 10⁻⁷ m s⁻¹ in the linear surface pressure construction
  (Figure 7c) to 10⁻¹² m s⁻¹ (Figure 7d)."*
- §5.2, 2-D idealised ice shelf — *"velocities in the idealized 2D ice shelf
  test cases have been reduced from order 10⁻² m s⁻¹ (not shown) to acceptably
  low velocities … order 10⁻⁹ m s⁻¹ or less."* (Figure 9b, the σ case, has a
  ±5 × 10⁻¹² m/s colourbar.)
- §6 — *"we achieve velocities of order 10⁻⁹ m s⁻¹ or smaller in all test
  cases (Figure 8)."*

Their three corrections are the **nonlinear surface pressure reconstruction**
(§3.1), the **interior reference interface** with its flattest-interface
fallback (§3.2), and **MWIPG**, mass weight in pressure gradient (§3.3.2).
**Roundabout implements none of them.** At Yung et al.'s own 10-day horizon
this build sits at `|u|_rms = 6.1e-05 m/s`.

Read that comparison with the configuration differences in front of you, both
of which cut against us:

- **Their runs are damped, ours is not.** The icemount tests carry a Laplacian
  viscosity of 1000 m² s⁻¹ plus top and bottom drag (`c_D = 0.002`,
  `U_bg = 0.05 m/s`, ~1-day spin-down), so a steady spurious force is
  *arrested* at `a/r` rather than integrated. These files carry no dissipation
  at all, by design, so a regression cannot hide behind viscosity. With
  ISOMIP+'s own `nu_h = 6.0 m² s⁻¹` this file read `3.73E-10`
  (`|u|_rms = 2.7e-05 m/s`) at day 10 on the legacy split — still above their
  *uncorrected* σ
  figure.
- **We quote a domain rms, they quote a domain maximum.** `sqrt(2·En)` is the
  rms over every wet face; their `max|u|` is the single worst cell. The rms is
  the smaller of the two, which makes the gap, if anything, wider than the
  table suggests.

The gap is the work, not the measurement. Phase 6 is where it closes.

## Known limits of these cases

- **Single rank.** `&ocean_cavity_dyn_nml` fences on `px*py > 1` in v1.
- **No steep calving cliff.** The front here is an 11 m step, chosen so the
  biggest `Δe` in the domain is the 13.8 m interior slope step. A 100 m
  ISOMIP+-style cliff is cubic-ally worse and would dominate the run; Yung
  et al. §3.4.2 record that MOM6-LAYER suppresses its calving criterion near
  the ice front "to minimise pressure gradient errors", and MPAS-Ocean caps
  the coordinate's Haney number at 5 for the same reason. Measuring the cliff
  regime honestly needs the corrections first.
- **No melt, no thermodynamics.** `&ocean_cavity_melt_nml` is off: a melt flux
  is an energy source and would destroy the "any motion is spurious" property.
- **Linear EOS, so the in-EOS load is inert.** Configure warns that
  `&ocean_psurf_nml in_eos` is off and the in-situ EOS therefore ignores the
  5.26 MPa load. With `eos = "linear"` there is no pressure dependence at all,
  so the warning is informational here; it matters for the nonlinear
  thermodynamic cases.
