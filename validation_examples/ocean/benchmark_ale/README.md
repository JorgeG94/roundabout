# benchmark_ALE — rdb's coupled-physics ALE integration target

Our analogue of MOM6's `ocean_only/benchmark_ALE`: a spherical, stratified,
wind-forced seamount basin with the **full ALE + coupled-physics stack on a
small cheap grid**. It is deliberately a *moving target* — the north-star
integration test we grow toward. Every capability we add (FNC1 stretched z*,
GM/Redi, MEKE, Henyey background, …) flips one more knob in `benchmark_ale.nml`
from a commented `TODO` to `ON`, and this README tracks that progress.

It is not a closed-form analytical regression (those live in `tests/`). Its job
is to **exercise the whole pipeline at once and surface integration bugs** — the
NaN-on-day-1, the silent conservation leak, the closure that destabilises only
in the presence of another closure — that no single-physics unit test can see.

## Sizes / tiers

| Tier | Grid | Run | Purpose |
|---|---|---|---|
| **this nml** | 120 × 60 × 75 | 5 day | MOM6-exact grid, 5-day spin-up; ~25 s on one V100, wind-driven gyre develops |
| MOM6 side-by-side | 120 × 60 × 75 | set `t_end=0.5` | matches MOM6 `DAYMAX` for a direct stats A/B; ~0.7 s |
| fast dev | 60 × 30 × 24 | shrink `nx/ny`, `nz_layers` | sub-second smoke iteration |

The committed nml is the MOM6-faithful geometry: a 100° × 50° spherical sector
from `SOUTHLAT=-45°`, `MAXIMUM_DEPTH=5500 m`, `DT=1800`/`DT_THERM=3600`
(`dt_therm_ratio=2`); MOM6's `DAYMAX` is 0.5 day but we run 5 days so the gyre
spins up. The only differences vs MOM6 are the physics gaps tracked below. Long-time fields (`t_end`, `dt_out`,
`status_interval`) are written in **days** via `&time_nml time_unit="day"`; the
CFL-bounded `dt_fixed` stays in seconds.

## How to run

```bash
# NVHPC + NetCDF toolchain loaded in this shell — exactly one per shell
cd build                            # GPU app build (TESTING=OFF)
OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=<gpu> \
   ./rdb ../validation_examples/ocean/benchmark_ale/benchmark_ale.nml
```

Watch the console status line (every 6 h sim time): SSH / KE / age should stay
finite and bounded. A NaN or runaway KE is a finding — capture the namelist
diff and the first bad step, not a "it broke."

## MOM6 knob → rdb status

✅ shipped & ON in the nml · 🟡 partial / lite path · ⬜ gap (commented TODO)

| MOM6 (benchmark_ALE) | rdb knob | Status |
|---|---|---|
| `OCEAN_GRID` spherical seamount | `&ocean_grid_nml grid_config="spherical"` + `topo_config="seamount"` | ✅ |
| `ROQUET_RHO` EOS | `&ocean_eos_nml eos="roquet_spv"` | ✅ |
| stratified T(z) IC | `T_init_surface/T_init_bottom` | ✅ |
| `WIND_CONFIG=gyres` | `&ocean_topo_nml wind_config="2gyre"` | ✅ |
| `SADOURNY75_ENSTRO` Coriolis | `&ocean_coriolis_nml form="sadourny"` | ✅ |
| FV PGF + `MASS_WEIGHT_IN_PRESSURE_GRADIENT` | `&ocean_pgf_nml form="fv_mom6" mass_weight=.true.` | ✅ |
| `SMAGORINSKY_AH` biharmonic | `&ocean_hvisc_nml smag_ah=.true.` | ✅ |
| `BEBT` barotropic blend | `&ocean_bt_nml bebt=0.2` | ✅ |
| `HBBL` + `CHANNEL_DRAG` bottom drag | `&ocean_bdrag_nml hbbl channel_drag=.true.` | ✅ |
| `DT_THERM` thermo cadence | `&ocean_vmix_nml dt_therm_ratio=2` | ✅ |
| `ENERGETICS_SFC_PBL` (EPBL) | `&ocean_epbl_nml enable=.true.` | ✅ |
| `USE_LA_LI2016` Langmuir | `&ocean_epbl_nml use_lt=.true.` | ✅ |
| `USE_JACKSON_PARAM` kappa-shear | `&ocean_kappa_shear_nml enable=.true.` | ✅ |
| `USE_IDEAL_AGE_TRACER` | `&ocean_tracers_nml enable_ideal_age=.true.` | ✅ |
| ALE conservative remap | always on (PPM); `vcoord_type="zstar_sigma"` | 🟡 lite z*; full GVC below |
| `RESTOREBUOY` SST/SSS restore | `q_heat`/`q_salt` forcing knobs | 🟡 partial |
| **`FNC1` stretched z* resolution** | uniform `dsig` only | ⬜ |
| **`vcoord_type="hycom"`** (isopycnal interior + z* floor) | shipped on `feat/ocean-vcoord-hycom` | ⬜ *(pending merge)* |
| **`remap_method="pqm"`** (piecewise-quartic) | shipped on `feat/ocean-remap-pqm` | ⬜ *(pending merge)* |
| **`regrid_time_scale`** grid time-filter | shipped on `feat/ocean-regrid-refine` | ⬜ *(pending merge)* |
| **`THICKNESSDIFFUSE` / GM** bolus transport | — | ⬜ |
| **`USE_STORED_SLOPES`** isopycnal slopes (Redi) | — | ⬜ |
| **`USE_MEKE`** eddy kinetic energy | — | ⬜ |
| **`USE_VARIABLE_MIXING` / `RESOLN_*`** Visbeck | — | ⬜ |
| **`HENYEY_IGW_BACKGROUND`** lat-dependent bg κ | — | ⬜ |
| **`USE_LOTW_BBL_DIFFUSIVITY`** law-of-the-wall BBL | — | ⬜ |
| **`FRAZIL` / `TFREEZE`** freezing | — | ⬜ |
| **`USE_BODNER23`** restratification | FK11 MLE (`use_mle`) is the analogue | 🟡 |
| **`TOPO_CONFIG="benchmark"`** exact topo | `"seamount"` substitute | 🟡 |

## Roadmap to the full benchmark

The remaining ⬜ blocks form two campaigns (see the ocean roadmap
`docs/ROADMAP_OCEAN.md` and `docs/CLOSURE_MATRIX.md`):

1. **GVC merge** — `hycom`, `pqm`, `regrid_time_scale`, `FNC1`: the vertical-
   coordinate fidelity. Three already shipped-pending-review; flip the commented
   knobs once they land on `main`.
2. **Mesoscale-eddy block** — isopycnal slopes → GM (`THICKNESSDIFFUSE`) → Redi →
   MEKE → variable/Visbeck mixing. A multi-capability campaign of its own; this is
   what makes a coarse benchmark_ALE *eddy-permitting* the MOM6 way.

When all rows are ✅ on the MOM6-exact 120×60×75 grid, rdb
reproduces the MOM6 benchmark_ALE coupled-physics envelope and the two can be
run side by side.
