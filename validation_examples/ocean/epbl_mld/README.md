# EPBL vs KPP — mixed-layer-depth comparison basin

The canonical `double_gyre_mom6.nml` is adiabatic (NK=2 reduced
gravity), so it cannot exercise a surface boundary layer.  This setup
is the thermally-active sibling built for comparing the surface
mixing schemes: flat 600 m basin, 30 × 20 m sigma layers, the same
2gyre wind / beta plane as the canonical case, plus uniform 100 W/m²
surface cooling and a linear T(z) stratification (10→18 °C,
N² ≈ 2.5e-5 s⁻²).

Three runs, identical except for the boundary-layer scheme:

| nml | scheme |
|---|---|
| `kpp_basin.nml` | PP81 interior + KPP overlay (current default) |
| `epbl_basin.nml` | PP81 interior + EPBL (`&ocean_epbl_nml enable`, OM4 mstar) |
| `epbl_lt_basin.nml` | as EPBL + LF17 Langmuir enhancement (`use_lt`) |

Expected physics (30 days): combined wind + convective deepening to
an MLD of roughly the encroachment scale `sqrt(2|B0|t)/N ≈ 97 m`,
spatially shaped by the gyres.  EPBL and KPP should agree on the
gross deepening; the differences (entrainment rate, transition-layer
sharpness, Langmuir deepening of the wind-driven component) are the
point of the comparison.

## Run

```bash
CUDA_VISIBLE_DEVICES=<gpu> ./build/rdb validation_examples/ocean/epbl_mld/kpp_basin.nml
CUDA_VISIBLE_DEVICES=<gpu> ./build/rdb validation_examples/ocean/epbl_mld/epbl_basin.nml
CUDA_VISIBLE_DEVICES=<gpu> ./build/rdb validation_examples/ocean/epbl_mld/epbl_lt_basin.nml
python3 validation_examples/ocean/epbl_mld/compare_epbl_kpp.py
```

(Output directories `out_{kpp,epbl,epbl_lt}_basin/` are created
wherever you launch from; pass `--base-dir` to the script if that is
not the repo root.)

The script computes a density-threshold MLD (Δσ = 0.03 kg/m³, de
Boyer Montégut convention) from the T output of **all** runs — the
apples-to-apples measure — and overlays EPBL's own `MLD_EPBL`
diagnostic, then writes `epbl_mld_comparison.png` and a summary
table.  The Δρ-criterion-vs-`MLD_EPBL` gap is itself informative:
the diagnostic is the *active mixing* depth, which leads the density
mixed layer during deepening.

## Knobs worth sweeping during tuning

- `&ocean_epbl_nml mstar_scheme` (`om4` vs `constant`), `mstar_cap`
- `tke_decay` (2.5 default — controls the wind-mixing depth under
  rotation), `translay_scale`
- `q_heat` in `&ocean_thermo_nml` (e.g. −400 for storm cooling, +100
  for restratification)
- `use_lt` + `lt_lac1..5` (the stability-modified La; zero them for
  the raw LF17 La)
