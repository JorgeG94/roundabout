# Ocean dyn-core: what knob should I turn?

A symptom → fix guide for the `sim_type='ocean'` C-grid path, focused on
**real-bathymetry / regional** runs (steep shelf-breaks, coastlines, thin
shelf layers) — the regime that exposes instabilities flat-bottom and
idealized basins never hit.

> **The first thing to try, for almost any steep-topography blow-up:
> `&ocean_bt_nml bebt = 0.1`.** It is the single highest-leverage knob and
> the most common root cause (see #1). Set it before chasing anything else.

## Triage table

| Symptom | Most likely knob | Group |
|---|---|---|
| Real-coast/steep-shelf run blows up after days; **barotropic**, mass-conserved, **death-day ∝ dt** | **`bebt = 0.1`** | `&ocean_bt_nml` |
| NaN within hours at a **thin shallow-shelf bottom layer** | `implicit = .true.` | `&ocean_bdrag_nml` |
| Day-1/2 death at a **steep shelf-break**, density-driven | `form = "fv_mom6"` + `mass_weight = .true.` | `&ocean_pgf_nml` |
| Grid-scale noise hugging coasts; momentum diffusing *across* land | `stress_tensor = .true.` | `&ocean_hvisc_nml` |
| Under-damped flow in **narrow channels / over sloping shelves** | `channel_drag = .true.` + `cdrag_side` | `&ocean_bdrag_nml` |
| Spurious bottom-layer exchange at a steep wet-wet face at **marginal CFL** | `vol_cfl = .true.` | `&ocean_continuity_nml` |
| Barotropic NaN in laminar spin-up (deep, ~2 km) | `auto_n_inner = .true.` | `&ocean_bt_nml` |
| Nearshore KE not conserved (velocity-form Sadourny) | `form = "sadourny_energy"` (or `"sadourny_hk"`) | `&ocean_coriolis_nml` |
| Smagorinsky biharmonic NaN at small dx/dt | cap `nu_4_max` | `&ocean_hvisc_nml` |
| Real bathy: land "drowned" / continents not walling | regenerate bathy with land preserved (`b < 2 m`) | bathy prep |
| Isolated wet cells / 1-cell channels / NaN at the coast | remove degenerate 1-cell features (topology cleanup) | bathy prep |

---

## #1 — Undamped barotropic free surface (`bebt`)

**Symptom.** A real-coastline / steep-topography run is calm for days, then a
**free-surface (SSH) mode grows exponentially** and trips CFL → NaN. Tells:
the growth is **barotropic** (persists with uniform density — no
stratification needed), **mass-conserving** (it's an oscillation, not a
leak), **immune to bottom/side drag and to lateral viscosity** (those act on
layer momentum, not the barotropic substep), and the **death day scales with
the timestep** (halving dt ~roughly doubles the survival time — the signature
of an under-damped numerical mode, not a physical forcing).

**Cause.** The split-explicit barotropic substep is **forward-backward Euler**,
which is *neutral* (non-dissipative). With `bebt = 0` (the historical default)
there is **no damping of barotropic gravity modes**, so an under-resolved
mode at a steep `H` step (e.g. a 100 m shelf next to 5000 m, or against a
wall) grows unchecked. The lateral-viscosity closures (even `stress_tensor`)
act on the *layer* momentum and never reach the fast barotropic loop.

**Fix.** `&ocean_bt_nml bebt = 0.1` (MOM6's default; valid 0→1, "must be
greater than about 0.05 in practice" per MOM6) — a
backward-Euler velocity-projection bias (`ubt_trans = (1+bebt)·ubt^n −
bebt·ubt^{n-1}`) that damps the barotropic gravity modes. `bebt = 0` ⇒
bit-identical to the historical neutral scheme. We verified `bebt = 0.2`
(a touch more damping than MOM6's 0.1) gives the 30-day flat Tasmania run;
`0.1` matches MOM6 and is the recommended value.

**Evidence (Tasmania real-coast, 5 km, 20 σ-layers).** `bebt = 0` → NaN ~day 2;
with four MOM6 fidelity ports added but `bebt = 0` → still NaN ~day 5;
**`bebt = 0.2` alone → 30 days, `En` dead-flat, mass-exact.** The barotropic
damping was the whole cure; the ports only delayed the symptom.

> **Recommendation:** set `bebt = 0.1` (MOM6 default) for ANY steep-topography / real-coast
> run. (It is *not yet* the global default because flipping it perturbs the
> idealized configs validated at `bebt = 0`, e.g. double-gyre — a default
> change needs a revalidation pass. Until then: opt in per config.)

---

## #2 — Thin shelf bottom layer, explicit drag CFL (`implicit`)

Quadratic bottom drag `du/dt = −C_d·|U|·u/h` is **explicit** by default and
goes unstable when `C_d·|U|·dt/h > 1` — easy on a thin shelf bottom layer
(5 m layer, dt=300, C_d=2.5e-3 ⇒ |U|≈6.7 m/s). `&ocean_bdrag_nml implicit =
.true.` makes it backward-Euler (`u/(1+dt·C_d·|U|/h)`), unconditionally
stable for any `h` (matches MOM6's implicit bottom-BC drag). Default-off ⇒
bit-identical.

## #3 — Steep shelf-break baroclinic PGF (`mass_weight`)

The FV-PGF interface integral is a layer-midpoint average; at a steep
unequal-depth face it injects a spurious horizontal density gradient (the
classic σ/z* shelf-break PGE). `&ocean_pgf_nml form="fv_mom6", mass_weight=
.true.` ports MOM6's `hWght` mass-weighting (blend toward the thinner column),
cutting the shelf-break PGF error ~50×. `hWght=0` on aligned columns ⇒
bit-identical.

## #4 — Coastline grid-noise / momentum across land (`stress_tensor`)

The default lateral viscosity is a velocity Laplacian with global static
caps and **no land masking inside the stencil** (it diffuses momentum across
coastlines). `&ocean_hvisc_nml stress_tensor=.true.` switches to MOM6's
thickness-weighted stress-divergence with a **per-cell CFL viscosity limiter**
(`bound_coef`, lets you crank scale-selective viscosity on a non-uniform grid
without NaN) and **coast-masks** the strain/stress. Default-off ⇒ bit-identical.
The biharmonic add-on (`nu_4` / `smag_ah` / `leith_biharm`) composes with this
path (it previously did not — an early return silently disabled it whenever
`stress_tensor=.true.`); the harmonic part is momentum-conserving and
coast-masked, the velocity-form biharmonic part is not.

## #5 — Narrow channels / sloping shelves (`channel_drag`)

MOM6 applies lateral side-wall Rayleigh drag in partially-blocked cells to
*every* layer intersecting sloping bathymetry; we had only bed-layer drag.
`&ocean_bdrag_nml channel_drag=.true., cdrag_side=2.5e-3` adds it (implicit,
thin-layer stable). All-wet/flat ⇒ no-op.

## #6 — Marginal-CFL continuity at a steep wet-wet face (`vol_cfl`)

Continuity-PPM used the CFL-free PPM edge value for the face thickness (the
CFL→0 limit). `&ocean_continuity_nml vol_cfl=.true.` adds MOM6's swept-volume
`O(CFL)` term. Matters only at marginal CFL across a steep thickness jump;
negligible (and `cfl=0` bit-identical) otherwise.

## #7 — Barotropic CFL in spin-up (`auto_n_inner`)

Fixed `n_inner` can under-resolve the 2-D gravity-wave CFL in deep water →
barotropic NaN. `&ocean_bt_nml auto_n_inner=.true.` derives `n_inner` from
the 2-D `c·dt·√2/dx` bound each step.

## Bathymetry prep (real coastlines)

- **Land must survive as `b < 2 m`** (`LAND_DEPTH_THRESHOLD`) or the static
  land mask won't wall it. Stamp land at `b = 0` — clamp ocean depths first,
  *then* land (order matters).
- **Remove degenerate single-cell features** (isolated wet cells, 1-cell
  bays/channels) — a wet cell with ≥3 land neighbours has all its face
  metrics zeroed → no flux path → Inf. Apply a `≥3`-land-neighbour pass, and
  for aggressive cleanup a 1-cell-channel + largest-connected-ocean pass.
  Topology only — depths untouched.

---

## The method (how to find the knob)

When a real-bathy run blows up, **bisect before you theorize** — these four
tests localized the Tasmania residual to a single subsystem in an afternoon:

1. **Homogeneous run** (`T_init_surface = T_init_bottom`): identical death ⇒
   the mode is **barotropic / topographic**, not tracer/baroclinic. Different ⇒
   it's density-driven (EOS/PGF/tracer-advection).
2. **Submerge the land** (set land cells to deep ocean): unchanged ⇒ it's the
   **wet-wet bathy gradient**, not coast/land masking. Fixed ⇒ it's the coast.
3. **Halve dt**: death-day scales with dt ⇒ an **under-damped numerical mode**
   (look for missing dissipation); death-day unchanged ⇒ a fixed forcing/CFL.
4. **Toggle dampers** (drag, viscosity): no effect ⇒ the mode is *not* in that
   operator (e.g. drag-immune ⇒ it's in continuity/PGF/barotropic, not friction).

Then **localize**: dump the max-|U| (or max-|SSH|) cell + frame over the
calm→blowup window — is it bed vs surface, shelf vs deep, coast vs interior,
free-surface vs layer? That + the bisection names the subsystem; only then
read the kernel and pick the knob (or, if there's no knob, the port).

The Tasmania chain: homogeneous (→ barotropic) + submerge (→ wet-wet shelf) +
dt-halve (→ under-damped numerical) + drag/visc-immune (→ not friction) +
SSH-mode dump (→ free surface) + "is there barotropic viscosity?" (→ none) ⇒
**`bebt`**. One knob.
