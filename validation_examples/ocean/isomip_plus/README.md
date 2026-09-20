# ISOMIP+ — Ocean0, Ocean1, Ocean2

Configurations for the ocean-only experiments of

> Asay-Davis, X. S., Cornford, S. L., Durand, G., Galton-Fenzi, B. K.,
> Gladstone, R. M., Gudmundsson, G. H., Hattermann, T., Holland, D. M.,
> Holland, D., Holland, P. R., Martin, D. F., Mathiot, P., Pattyn, F. and
> Seroussi, H. (2016), **"Experimental design for three interrelated
> marine ice sheet and ocean model intercomparison projects: MISMIP v. 3
> (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1)"**,
> *Geoscientific Model Development* **9**, 2471–2497,
> doi:[10.5194/gmd-9-2471-2016](https://doi.org/10.5194/gmd-9-2471-2016).

Every section, equation and table number below refers to that paper.

**These are configuration deliverables, not validated results.** They pass
`validate_config`, they run, and the parameter mapping is auditable
row-by-row. They have NOT been spun up, and nothing here has been tuned to
a melt rate — the shipped `gamma_t` is the protocol's *starting guess*
(Sect. 3.2.1), not a calibrated value. The blockers to a publishable
Ocean0 are listed at the bottom.

| File | Init | Far-field restoring | Draft | Duration |
|---|---|---|---|---|
| `ocean0.nml` | WARM | WARM | ISOMIP+ file | 1 yr |
| `ocean1.nml` | COLD | WARM | ISOMIP+ file | 20 yr |
| `ocean2.nml` | WARM | COLD | ISOMIP+ file (Ocean2 / retreated) | 20 yr |
| `ocean0_idealised_draft.nml` | WARM | WARM | **analytic** — runs today, no download | 1 yr |

---

## Where to get the geometry

The bedrock is analytic and ships in the model
(`&ocean_topo_nml topo_config = "isomip_plus"`, their Eqs. 1–4 + Table 1).
**The ice draft is not.** Sect. 3.1.1:

> "Because the ice-draft topography is derived from ice-sheet model
> results, it cannot be described by an analytic function. Instead, both
> the topography used for Ocean0–2 and the snapshots used to produce the
> dynamic topography for Ocean3–4 come from MISMIP+ BISICLES results, and
> are available in NetCDF format for download (Cornford and Asay-Davis,
> 2016)."

That citation is the paper's own data reference for the MISMIP+/ISOMIP+
geometry; follow it from the article's reference list
(doi:[10.5194/gmd-9-2471-2016](https://doi.org/10.5194/gmd-9-2471-2016)),
which is the authoritative pointer. Ocean0 and Ocean1 use the **initial
steady state** of MISMIP+ Ice1; Ocean2 uses the state at the **end of
Ice1r** (Sect. 3.2.3) — two different files.

**You must regrid it.** The protocol ships the geometry on a uniform 1 km
grid and says participants "are expected to interpolate the ice-sheet
topography to the ocean grid as part of whatever processing is required to
make the data ocean-model friendly" (Sect. 3.1.1). Roundabout's reader does
**no** horizontal interpolation — same rule as `bathymetry_file` and
`&ocean_zinit_nml file`. Apply the calving criterion (`H_calve = 100 m`,
Sect. 3.1.2) during that processing, as the paper instructs, because the
shipped file deliberately has it *not* applied.

The regridded file must satisfy the reader's contract:

* variable `iceDraft`, **rank 3**, declared `(time, y, x)` in CDL / C /
  Python (equivalently `(x, y, t)` in Fortran storage order) — this is
  already the ISOMIP+ file's own layout;
* a `time` coordinate variable; **record 1** is read and never re-read;
* horizontal extents `240 × 40`, matching `&grid_nml nx/ny`;
* values are the **elevation** `z_d ≤ 0`, which is why the namelists say
  `draft_sign = "elevation"`. Roundabout's `z_draft` is a DEPTH, positive
  down, so the reader negates. Getting this backwards is a fail-loud error
  (the non-negativity guard), never a silent 1000 m mistake.

---

## The protocol parameter map

"exact" = the protocol number reaches the kernel unchanged.
"approximated" = a defensible mapping onto a different-but-equivalent
knob, stated. "NOT AVAILABLE" = this build cannot express it.

### Domain and discretisation

| Protocol | Value | Roundabout knob | Status |
|---|---|---|---|
| `x0` (Table 3) | 320 km | `&ocean_topo_nml x_origin = 320000.0` | exact |
| `Lx`, `Ly` (Table 3) | 480 km × 80 km | `&grid_nml nx=240, ny=40, dx=dy=2000` | exact |
| `Δx = Δy` (Table 4) | 2 km | `&grid_nml dx/dy` | exact |
| Vertical layers (Sect. 3.1.5) | 36, distribution free | `&nonhydrostatic_nml nz_layers = 36`, `&vcoord_nml vcoord_type="sigma"` | exact (sigma is a legal "distribution at the modeller's discretion") |
| Bedrock (Eqs. 1–4, Table 1) | `B0 −150`, `B2 −728.8`, `B4 343.91`, `B6 −50.57`, `x̄ 300 km`, `d_c 500 m`, `f_c 4 km`, `w_c 24 km`, `z_b,deep −720 m` | `&ocean_topo_nml topo_config="isomip_plus"`, `max_depth = 720.0` | exact (named constants in `rdb_ocean_state`, gated by `test_ocean_topo_isomip_plus`) |
| Ice draft | BISICLES NetCDF | `&ocean_cavity_dyn_nml draft_config="file"` | exact **once regridded** — see above. `ocean0_idealised_draft.nml` substitutes an approximation. |
| `x_calve` (Table 1 / Sect. 3.1.1) | 640 km | implicit in the file; `draft_x1 = 320000.0` (model x) in the idealised variant | exact / approximated |
| `H_calve` (Sect. 3.1.2) | 100 m (draft above ≈ −90 m) | applied during the offline regrid | **offline** — not a model knob |
| Minimum ocean column (Sect. 3.1.5) | value left to the modeller | `&ocean_cavity_dyn_nml h_min_cavity = 20.0` | exact-by-construction. The protocol offers two options — modify the topography, or "remove the column from the ocean (i.e. mark it as 'land')". Roundabout takes the second: `b − z_draft < h_min_cavity` ⇒ the column goes through the ordinary `seed_wet_mask_impl` land path (metric-zeroed walls, finite held land state). There is deliberately **no thin film** under grounded ice. |
| Latitude / Coriolis (Sect. 3.1.1) | f-plane at 75° S | `&physics_nml coriolis_f = -1.409e-4`, `&ocean_topo_nml coriolis_beta = 0.0` | exact (`−2Ω sin 75° = −1.4089e-4`) |

### Forcing and boundaries

| Protocol | Value | Roundabout knob | Status |
|---|---|---|---|
| Surface forcing (Sect. 3.1.3) | none | `wind_stress_x/y = 0`, no `&ocean_restore_nml`, no `q_heat`/`q_salt`, `sw_pen_frac = 0` | exact (and the melt path *refuses* a non-zero wind fail-loud, so this cannot drift) |
| Restoring region (Table 3) | `x_r0 = 790 km`, `x_r1 = 800 km` | `&ocean_bc_nml east="sponge", sponge_width = 5` (= 10 km / 2 km) | exact |
| `γ0` (Table 3) | 10 day⁻¹ (τ = 0.1 day) | `&ocean_bc_nml sponge_strength = 1.1574074e-4` 1/s | exact |
| `γ(x)` shape (Eq. 20) | linear from 0 at `x_r0` to `γ0` at `x_r1` | `&ocean_sponge_nml ramp = "linear"` | exact at cell centres — `(band − d − 0.5)/band` IS Eq. (20) evaluated at the cell centre when the band spans `[x_r0, x_r1]` |
| Restoring tendency (Eqs. 18–19) | `∂T/∂t = −γ(x)(T − T_res(z))`, S twin | `&ocean_sponge_nml enable, target_source="linear_z", relax_tracers=.true., relax_uv=.false.` | exact — with one deliberate divergence: the discrete step is the EXACT exponential `φ ← φ·e^(−γΔt) + φ_ref(1−e^(−γΔt))` rather than a first-order Euler/backward-Euler step, which is the same ODE solved better |
| `T_res(z)`, `S_res(z)` (Eqs. 21–22) | linear from surface to `z_b,deep = −720 m` | `&ocean_sponge_nml lin_t_ref/lin_dt_dz/lin_s_ref/lin_ds_dz`, **re-evaluated on the live layer geometry every outer step** | exact. **Eq. (22) as printed has a typo** — `S_res(z) = S0 + (S_bot − T0)·z/z_b,deep` must be `(S_bot − S0)`. Signs: z is positive UP here, so WARM `dT/dz = −(1.0 − (−1.9))/720 = −4.0277778e-3` (WARM is thermally *unstable*, stabilised by salt) and `dS/dz = −(34.7 − 33.8)/720 = −1.25e-3`; COLD `dT/dz = 0`, `dS/dz = −(34.55 − 33.8)/720 = −1.0416667e-3`. |
| Lateral BCs (Sect. 3.1.4) | **no-slip** at all walls | `&ocean_bc_nml west/south/north = "wall"` | **NOT AVAILABLE** — Roundabout's static land mask is FREE-SLIP by construction (it zeroes the face metrics). The paper itself notes "free-slip or open boundary conditions may be more physically justifiable but no-slip boundary conditions are likely to be supported by the largest number of models", and asks that a different choice be noted. This is that note. |
| Melting/drag on vertical ice faces (Sect. 3.1.4) | none | none | exact (not implemented, which is what the protocol asks for) |
| Initial state (Sect. 3.1.4) | at rest, horizontally uniform T/S | `&ocean_zinit_nml source="linear"` (+ zero velocity seed) | exact. The zinit "linear" source evaluates `T(z)`/`S(z)` at each layer centre's TRUE GEOPOTENTIAL depth measured through the draft, so the isopycnals are flat in z under a tilted lid — `&tracer_nml T_init_surface` would lay the profile out along the tilted sigma coordinate and would NOT be a state of rest. |

### Physics

| Protocol | Value | Roundabout knob | Status |
|---|---|---|---|
| EOS (Eq. 23, Table 4) | `ρ = ρ_ref[1 − α_lin(T−T_ref) + β_lin(S−S_ref)]`, `ρ_ref 1027.51`, `T_ref −1`, `S_ref 34.2`, `α_lin 3.733e-5`, `β_lin 7.843e-4` | `&ocean_eos_nml eos="linear"`; `&ocean_ic_nml rho_0=1027.51, T_ref=-1.0, S_ref=34.2, alpha_T=3.8356948e-2, beta_S=8.0587609e-1` | exact. **Units trap:** Roundabout's `alpha_T`/`beta_S` are DIMENSIONAL (kg m⁻³ per unit) while the protocol's are FRACTIONAL (they multiply `ρ_ref` inside the bracket). `alpha_T = ρ_ref·α_lin`, `beta_S = ρ_ref·β_lin`. Feeding the fractional numbers in directly under-states the density response ~1000× and looks like a plausible weak run. |
| Liquidus (Eq. 25, Table 4) | `λ1 −0.0573`, `λ2 0.0832`, `λ3 −7.53e-8` | `&ocean_eos_nml tfreeze_set = "isomip"` | exact (named set; the default `"seaice"` SIS2 set is ~0.03 °C away, enough to flip the SIGN of melt over a 0.03 °C band) |
| Melt formulation (Sect. 3.1.8, Eqs. 24–27) | three-equation, constant `Γ_T`/`Γ_S` | `&ocean_cavity_melt_nml exchange_law="const_gamma"` | exact |
| `Γ_T` (Sect. 3.2.1) | starting guess 2.2e-2, **to be tuned** to ⟨m_w⟩ = 30 ± 2 m a⁻¹ | `gamma_t = 2.2e-2` | exact as the starting guess; **NOT tuned** (that is the Ocean0 calibration, Eq. 37, and needs a converged run) |
| `Γ_S = Γ_T/35` (Table 4) | derived | `gamma_s = -1.0` ⇒ derived | exact |
| `C_D,top`, `C_D,bot` (Table 4) | 2.5e-3 each | `&ocean_tdrag_nml cd`, `&ocean_bdrag_nml cd` (+ `&ocean_cavity_melt_nml cdrag_top`, which must agree — fail-loud) | exact |
| `u_tidal` (Table 4) | 0.01 m s⁻¹, melt `u*` only, never the drag | `&ocean_cavity_melt_nml u_tide = 1.0e-2` | exact (the kernel documents the same separation) |
| `κ_i` (Table 4) | 0 (perfectly insulating) | `ice_conduction = "insulating"` | exact |
| `c_w`, `L` (Table 4) | 3974, 3.34e5 | model constants | exact |
| `ρ_fw`, `ρ_sw` (Table 4) | 1000, 1028 | model constants / `rho_ice` | approximated — the load uses the Boussinesq-isostatic `ρ_0·g·z_draft`, which is what ISOMIP+ prescribes for models that take a draft rather than a pressure ("the pressure can be derived from the ice draft as `p_zd = −ρ_sw g z_d`", Sect. 3.1.1) |
| `ν_H` (Table 4) | 6.0 m² s⁻¹ harmonic | `&ocean_hvisc_nml nu_h = 6.0`, `nu_4 = 0.0`, no flow-aware closure | exact |
| `κ_H` (Table 4) | 1.0 m² s⁻¹ harmonic | `&ocean_hdiff_nml kappa_h = 1.0` | exact |
| `ν_stab`, `κ_stab` (Table 4) | 1e-3, 5e-5 m² s⁻¹, **constant** | `&ocean_vmix_nml use_closure=.true., use_kpp=.false., pp81_nu0=0.0, pp81_nu_bg=1.0e-3, pp81_kappa_bg=5.0e-5` | exact. PP81 with `nu0 = 0` degenerates to `kv ≡ nu_bg`, `kt ≡ ks ≡ kappa_bg` at every interface with no Richardson dependence left — i.e. exactly the protocol's constant-coefficient harmonic vertical mixing. KPP is off deliberately: a boundary-layer scheme is a different experiment (Sect. 3.1.6 allows it but asks that it be documented; this is the COM choice). |
| `ν_unstab`, `κ_unstab` (Table 4) | 0.1 m² s⁻¹ when stratification is unstable | `&ocean_conv_nml enable=.true., kd_conv=0.1, prandtl_conv=1.0, n2_thresh=0.0` | exact (`prandtl_conv = 1` ⇒ `Kv_conv = Kd_conv = 0.1`; the trigger is `N² < 0`) |
| Virtual vs volume melt flux (Eqs. 28–33) | either, documented | `&ocean_cavity_melt_nml freshwater = "virtual"` (default, what these files use) or `"mass"` | **BOTH AVAILABLE.** `"mass"` is the real volume flux: `dh = m·dt/ρ₀` on the top layer with `d(hS) = dh·s_ice`, `d(hT) = dh·T_b`. These namelists stay on `"virtual"` so the shipped configuration is unchanged; flip the one key for the volume form |
| Evaporative mass removal (Eqs. 34–36) | needed only for a VOLUME-flux model | `&ocean_cavity_melt_nml volume_compensation = "uniform_open_ocean"` | **AVAILABLE but NOT the protocol's choice for Ocean0-2.** Sect. 3.1.3 leaves the restored configurations uncompensated — their only seam is a restoring band that moves no volume — and allows compensation for the closed Ocean3/4. The knob is a uniform per-unit-area removal over the uncovered wet cells, each parcel carrying that cell's own `T`/`S`; it is NOT the protocol's evaporation field. Default `"none"` |

### Not prescribed by the protocol

| Choice | Value | Why |
|---|---|---|
| Outer time split | `pred_corr` (the default) | The protocol says nothing about the time split. `pred_corr` is the model default and has the better internal-wave behaviour. |
| `dt` | 300 s | Not prescribed ("time discretisation: e.g. ... explicit"). 300 s keeps MaxCFL ≈ 0.016 on the 1-day check with 36 sigma layers over a column that thins to ~25 m near the grounding line. |
| Barotropic substeps | `auto_n_inner = .true.` | Derived from the gravity-wave CFL at configure (28 substeps at 720 m). |
| PGF | `fv_mom6` + `p_top_in_bc` | Required by the cavity: the only form with an injectable top boundary condition, without which a varying draft is not represented. |
| `far_field_depth` | 10 m | Sect. 3.1.8 explicitly leaves the method to the modeller and asks it be documented: this samples a fixed 10 **metres** below the ice base, thickness-weighted with a partial last layer — never "layer nz", because that would make the melt rate a function of the vertical coordinate's cell thickness. |

---

## Running them

```bash
# needs no data:
./build/rdb validation_examples/ocean/isomip_plus/ocean0_idealised_draft.nml

# needs the regridded geometry next to the namelist (or an absolute
# draft_file path):
./build/rdb validation_examples/ocean/isomip_plus/ocean0.nml
```

`t_end` in each file is the **protocol** duration (1 yr for Ocean0, exactly
20 yr for Ocean1/Ocean2). For a smoke test, shorten it.

### The idealised draft

`ocean0_idealised_draft.nml` replaces the data file with a linear-in-x
draft anchored on the only two draft numbers the paper prints:

* the steady-state grounding line crosses the trough centreline at
  `x = 450 ± 10 km` (Sect. 2.1);
* the calving front is at `x = 640 km`, where ice thinner than
  `H_calve = 100 m` — "equivalent to an ice draft above ∼ −90 m"
  (Sect. 3.1.2) — has been removed.

A straight line through (450 km, 600 m) and (640 km, 90 m) gives
`d(draft)/dx = −2.6842e-3`, anchored at 975.79 m at `x = 310 km` of paper
x (10 km west of the model domain, so the ghost columns are iced rather
than stepping to open water at the wall). West of the grounding line the
formula is deeper than the bed, so those columns ground through
`h_min_cavity` and become land — which is what "grounded ice sheet" means
here. Measured at configure: the trough-centre grounding line lands at
`x = 449 km`, the draft at the calving front is exactly 90.0 m, and 3778
of 9600 interior columns (39.4 %) ground.

It is an **approximation of the shape, not of the ice sheet**: the real
draft varies in y, has a lateral shear margin, and is not linear in x. Use
it to exercise the configuration, never to produce an ISOMIP+ number.

### Measured 1-day checks (gfortran 15.1 Release, single rank, 2026-09-20)

`dt = 300 s`, 288 steps. `Salt`/`Temp` are the domain means from the
`[stats]` line.

| Case | En (day 1) | MaxCFL | Salt | Temp |
|---|---|---|---|---|
| `ocean0_idealised_draft` | 1.793E-06 | 0.0124 | 34.275 | −0.370 |
| `ocean0` (file draft) | 1.793E-06 | 0.0124 | 34.275 | −0.370 |
| `ocean1` (COLD → WARM) | 3.181E-06 | 0.0134 | 34.199 | −1.834 |
| `ocean2` (WARM → COLD) | 3.460E-06 | 0.0124 | 34.271 | −0.436 |

Ocean1 warms (−1.838 → −1.834 over the second half-day) and Ocean2 cools
(−0.431 → −0.436): the far field is pulling each toward a water mass it
did not start in, which is the property `target_source = "linear_z"`
exists for and which `target_source = "ic"` cannot express at all.

**Cross-toolchain.** The same three runs on the GPU build (nvfortran
26.5, `-stdpar=gpu`, `cc70`, one Tesla V100-DGXS) reproduce **every
printed digit** of the gfortran answers above — En, MaxCFL, Mass, Salt
and Temp, at day 0.5 and day 1, for `ocean0`, `ocean1` and
`ocean0_idealised_draft` alike — in 5.3–8.2 s instead of 149–181 s. That
covers the two new device-side code paths this configuration exercises:
the per-outer-step `linear_z` sponge-target refresh, and the cavity
draft read from NetCDF and sign-normalised on the host before the map.

`ocean0` and `ocean0_idealised_draft` agree to every printed digit because
the file used for that check was generated from the same analytic formula
— which is also an end-to-end confirmation that the NetCDF reader plus the
`draft_sign = "elevation"` negation reproduce the analytic setter exactly
on the interior.

---

## Blockers for a publishable Ocean0

Listed rather than worked around. The first is decisive on its own.

1. **CLOSED — the meltwater can now carry MASS.** `&ocean_cavity_melt_nml
   freshwater = "mass"` makes the meltwater a real Boussinesq volume
   source on the top layer (`dh = m·dt/ρ₀`, `d(hS) = dh·s_ice`,
   `d(hT) = dh·T_b`), with the virtual salt flux retained in `Q_salt`
   only as the `B_0` buoyancy forcing KPP/EPBL read and removed again
   from the tracer in the same stage.  Mass is a tracked budget source,
   so all three console residuals still sit at round-off at every step.
   **These namelists stay on the `"virtual"` default**, so they are
   unchanged; flip the one key for the volume form, and add
   `volume_compensation = "uniform_open_ocean"` if the sea-level rise of
   a closed box matters for the horizon you are running (ISOMIP+ itself
   leaves Ocean0-2 uncompensated).  What remains open is the TUNING
   (blocker 3): the mass form changes the answer, so `Γ_T` must be
   re-searched against it, not inherited from a virtual-salt run.
2. **Under-ice `u*` is not seen by the boundary-layer scheme.** The
   ice-shelf top drag publishes `stress_top`, but KPP/EPBL do not consume
   it. These namelists run with KPP off (the protocol prescribes constant
   vertical coefficients), so it is inert *here* — but a TYP-style
   configuration with a boundary-layer scheme would be missing the shear
   that drives sub-ice mixing.
3. **`Γ_T` is untuned.** Sect. 3.2.1 makes Ocean0 the calibration
   experiment: `Γ_T` must be searched until `⟨m_w⟩ = 30 ± 2 m a⁻¹` over
   `z_b < −300 m` and the final 6 months (Eq. 37). POP2x needed
   `Γ_T ≈ 0.11`, five times the starting guess shipped here. Nothing in
   this directory is tuned to a melt rate, deliberately.
4. **No-slip walls are NOT AVAILABLE** (Sect. 3.1.4). The static land mask
   is free-slip. The paper anticipates this and asks that it be noted.
5. **The sloping-lid resting-state growth mode is under investigation.**
   See `../ice_shelf_cavity/README.md`: a resting stratified cavity under a
   sloping lid holds the derived PGF-truncation plateau for ~20 days and
   then leaves it on a ~3.2-day e-folding before saturating. ISOMIP+'s own
   `ν_H = 6.0` only delays it ~10 days. It is bounded and diagnosed (it is
   entirely downstream of the sloping-lid PGF error), but it is a spurious
   energy source sitting under every run in this directory, and a 20-year
   Ocean1/Ocean2 integrates it for a long time.
6. **No-slip walls (4), the untuned `Γ_T` (3) and the sloping-lid mode (5)
   are the remaining blockers.** The grounded-cavity budget defect that used
   to be listed here is FIXED — see the next section — and so is the
   virtual-meltwater blocker (1).

These namelists are registered in
`tests/regression/stability_manifest.py` as
`isomip_plus_ocean0_idealised` (the idealised-draft file only; the three
file-backed ones need a download and stay out). It is a `forced` case run
for 1 simulated day at tier 1 and 100 steps at tier 2, asserting finite,
`conserve:{Mass,Salt,Heat}` at `1e-11` and the CFL guards — deliberately
NOT a melt-rate or a rest-state gate, because neither is meaningful until
the tuning blocker (3) and the sloping-lid mode (5) are closed. It is carried for the conservation assertion above
everything else: 39.4 % of its interior columns ground, the largest
grounded fraction anywhere in the corpus.

---

## Grounded columns and the budget (FIXED)

A column whose ice draft meets the bed (`b − z_draft < &ocean_cavity_dyn_nml
h_min_cavity`) is LAND, through the same `seed_wet_mask_impl` the
bathymetry uses. 3778 of this case's 9600 interior columns (39.4 %) are.

Until 2026-09-20 such a column was seeded wrong, and the console said so:

```
Mass : 1.277811397E+16  Error -5.005E-14
Salt : 4.379651448E+17  Error  6.088E-01  out 1.732E+06  src -4.348E+12
Heat :-4.723233226E+15  Error -7.994E-02  out 5.593E+06  src -1.062E+13
```

61 % of the initial salt content, appearing as a step change between step 0
and step 1 and then flat. **It was a real defect, not a baseline artefact.**
A grounded column's water column `b − z_draft` is NEGATIVE — −314 m on
average here — so its seeded `h_layer` was negative, and the land tracer
hold `val = hTr/max(h_old, H_VANISHED); hTr = val*H_VANISHED` is an exact
algebraic identity for `h_old ≤ H_VANISHED`: it left a FULL-COLUMN `hTr`
beside a thickness floored to `1.5e-4 m`. The budget latch integrated that;
the first ALE regrid discarded it (a land layer sits AT the vanish marker,
so the remap's `h > H_FLOOR` gate is false and it writes `hTr = 0`).

That also explains the shape of the numbers: the offset scaled with the
grounded columns' summed DEPTH DEFICIT `Σ(z_draft − b)`, not their area,
which is why 39.4 % and 14 % grounded gave 6.09E-01 and 1.87E-02 rather
than anything proportional to the fraction.

A land T-cell now holds `h_layer = H_VANISHED` and `hTr = 0` — the state
every vanished-gated operator already holds it at — and a grounded column's
barotropic datum is `0` rather than a negative water column. The same run
today:

```
Mass : 1.277811397E+16  Error -5.255E-14
Salt : 4.379651448E+17  Error -5.242E-14  out 1.732E+06  src -4.348E+12
Heat :-4.723233226E+15  Error  6.544E-14  out 5.593E+06  src -1.062E+13
```

The day-1 totals are unchanged to every printed digit: no wet cell was ever
affected (verified bitwise — `tests/test_ocean_cavity_grounded_budget.F90`
runs the same 64 cells as grounded ice and as ordinary island land and
finds every wet column identical). What moved is the step-0 latch, which
was the thing that was wrong.

The contract is in `src/core/ocean/README.md` ("the land-state contract"
and "the cavity datum contract"); the gate is
`tests/test_ocean_cavity_grounded_budget.F90`, which holds all three
relative residuals at `1e-12` for 20 steps of the full split solver, with
melt off and with melt on.
