# Cross-cutting contracts of the ocean dyn-core

Read this only when you touch one of the seams below. The per-module contracts are written once, in the
`!!` docstring of the module that owns them; the table at the top says where. What follows the table are
the few rules that have no single owner. If this file and the code disagree, the code wins and this file
is the bug.

## Where each contract lives

| Contract | Owner (read this) |
|---|---|
| Vanished-layer predicate, donor read and the consumer-substitution table | header of `src/shared_module_utilities/rdb_vanished_layer.inc` |
| Thin-layer constants `H_VANISHED` / `H_DIV_EPS` | `src/core/rdb_constants.F90` |
| Land state (`h = H_VANISHED`, `hTr = 0`) | `seed_land_tracer_hold_impl` and `ocean_state_seed_land_cells`, `src/core/ocean/state/ocean_state/rdb_ocean_state_seed.F90` |
| `stress_mag` vs `stress_shelf` (two upper-boundary stresses) | module header of `src/parameterizations/vertical/rdb_ocean_surface_stress.F90` |
| PGF reference densities `rho0` vs `rho_ref` | the two members in `src/pressure_force/rdb_ocean_pressure_force.F90` |
| One configured ρ₀ fanned out to forcing + KPP | `configure_ocean_reference_density`, `src/core/ocean/state/setup/rdb_ocean_setup_lateral.F90` |
| `eta_forcing` seam (tide + SAL + surface pressure) | module header of `src/core/ocean/forcing/rdb_ocean_p_surf.F90` |
| `p_top` seam (three pressures, rules for a joining builder) | comment on `p_top` in `src/core/rdb_multilayer_state.F90`; the third consumer (melt liquidus) in `configure_ocean_cavity_melt` |
| `F_slow` seam (which tendencies the fast loop sees) | `sum_slow_tendencies_into_F_slow`, `src/core/ocean/kernels/barotropic/rdb_barotropic_coupling.F90` |
| `bt_rem_u/v` multiplicative accumulator (one resetter) | `reset_bt_rem`, same file |
| Frozen-reference subtraction in the fast loop (PGF `η_PF`, Coriolis `cor_ref`) | `subtract_fast_cor_ref`, `set_cor_ref_velocity`, `set_fast_forcing_eta_pf`, same file |
| Cavity datum (`bt_H_ref = b − z_draft`, counted-once load) | module header of `src/core/ocean/state/rdb_ocean_cavity.F90`; assembly in `configure_ocean_cavity` |
| Surface-flux component ownership (one writer per component) | module header of `src/parameterizations/vertical/rdb_ocean_surface_flux.F90` |
| Real freshwater mass from basal melt | module header of `src/parameterizations/vertical/rdb_ocean_cavity_flux.F90`; summary below |
| Strided-hybrid column loops, scratch buffers | module header of `src/framework/rdb_scratch_3d.F90` |
| Joining the checkpoint | `ocean_state_build_restart_registry`, `src/core/ocean/state/ocean_state/rdb_ocean_state_restart.F90` |
| Diagnostics are pull-only | module header of `src/core/ocean/diag/rdb_ocean_diag.F90` |
| Order-invariant global sums | `src/framework/rdb_efp.F90` |
| Who reads / writes what in each slot | the `!!` docstring of each slot type, composed in `src/core/ocean/state/rdb_ocean_state.F90` |

## Conventions that apply everywhere

**Indexing.**
- Layers are bottom-up (ROMS style): `k = 1` is the bed, `k = nz` is the surface layer.
- `u_face_x(i,j,k)` sits on the face between cells `(i,j)` and `(i+1,j)`; its array runs from `i = 1` (west
  outer wall) to `i = nx+1` (east outer wall). `v_face_y` is the same in `j`.

**`is_init`.** Every slot type carries `logical :: is_init = .false.`. `init` sets it only after every
allocation and device attachment succeeded; `destroy` clears it before releasing anything. Always test
`slot%is_init`, never `allocated(slot%some_array)`: `allocated` sees only the host pointer, so it is true
for arrays that exist on the host but are not on the device yet. It is also the right guard against
double-init and use-after-destroy.

**Memory gating of default-off closures.** EPBL, kappa-shear, tidal mixing, GM, Redi, MLE, VarMix, MEKE and
the isopycnal slopes allocate only when their `enable` flag is set.
- Latch `enable` in `init_from_config` before `init`, so `if (this%slot%enable) call this%slot%init(...)`
  can read it.
- Gate the matching `enter_data` / `exit_data` calls in the parent walk the same way; those bodies carry no
  internal `allocated` guard.
- Add every array to the slot's `bytes()` through `arr_bytes`, so a gated-off slot counts 0 and the upfront
  memory report (reconciled against the measured device map) stays quiet.
- A gate *inside* a live slot (a scheme-variant buffer such as the FV_MOM6 reconstruction stack, or
  `continuity_t%windowed_advection`) follows the same shape, with two extra obligations. Prove per buffer,
  at the call site, that it is unreachable when not allocated: under `mem:separate` an unallocated array in
  a `do concurrent` does not reliably fault, it can read stale host data and give a plausible wrong
  answer; use a structural proof (a dispatch that returns before the buffer is named, or a consumer that
  already `error stop`s on the wrong variant). And default the gate to allocate-everything: direct
  `slot%init(...)` calls in tests set the discriminator *after* `init`, so only `init_from_config` turns
  the saving on. Where the discriminator is re-derived later, compare it to the latched value and
  `error stop` on a mismatch.

**One ρ₀.** `&ocean_ic_nml rho_0` lands on `eos%rho0`, the Boussinesq reference density of record.
Every slot that carries its own `rho0` copy must be assigned from it in its `configure_ocean_*` (or in
`configure_ocean_reference_density` for surface flux, surface stress, vmix and geothermal). A copy with a
type default that nothing assigns compiles, runs and conserves, and silently scales a force by
`1035/ρ₀`. `ocean_vmix_t%rho0` is read on the device, so it is correct only because the configure pass
runs before `ocean_state_enter_data`; anything that writes it later owes an `!$acc update device`.
`RHO_WATER` (console statistics, OBC mass flux, ice), `&ocean_ice_nml rho_ocean`, the melt law's own
constants and the `ICE_RHO_*` constants are deliberately *not* copies of ρ₀.

**Thin layers.** `H_VANISHED` (1.5e-4 m) is a skip/merge marker, never a positivity floor. `H_DIV_EPS`
(1e-20) is pure 1/0 armour. The runtime `angstrom_h` floor (`&ocean_isopycnal_nml`, Lagrangian only) lifts
`h` but not `hTr`, so it is non-conservative by a bounded amount per floored cell.

**Land.** A land T-cell holds `h = H_VANISHED` exactly and `hTr = 0` exactly, at `t = 0` and at every
step. The seed hands the budget latch the state every later step will have; any other seeded land content
is discarded by the first regrid, un-budgeted, and shows up as a step in the console error between
step 0 and step 1. Land is isolated from wet cells bit-for-bit
(`test_ocean_cavity_grounded_budget`). A grounded ice-shelf column goes through the wet mask, never a thin
film.

**Boundary-layer friction velocity.** KPP and EPBL both read `u*² = (stress_mag + stress_shelf)/ρ₀` and
nothing else. A new upper-boundary stress joins by the same test the `p_top` seam uses: if it is derived
from `tau`, it belongs in `stress_mag` and its writer owes `ocean_surface_stress_refresh_mag`; if it is a
separate momentum flux across the ocean top, it gets its own always-allocated companion and an inline
per-step fill, never a fold into `stress_mag`.

**One definition, many consumers.** Never hand-roll a vanished-layer test; use
`rdb_vl_*` from the `.inc`, and document a deliberate substitution in its consumer table.

**Kernels stay diagnostic-agnostic.** No `diag%send` in a kernel. New derived quantities go into the diag
registry as a fill.

**Global integrals.** A plain `reduction(+:)` plus allreduce drifts at round-off across decompositions.
Use the EFP seam (`rdb_efp`, `halo_allreduce_efp_list`) for any integral that must be rank-count
reproducible, and `efp_real_diff` for a difference against a latched reference.

## Surface flux: which paths change real column mass

The mass-flux components (`evap`, `lprec`, `fprec`, `vprec`, `lrunoff`, `frunoff`, `seaice_melt`) and
their `heat_content_*` companions are enthalpy and salt bookkeeping only. They do not change column mass.

| Path | Real column mass? |
|---|---|
| Ice-shelf basal melt, `&ocean_cavity_melt_nml freshwater = "mass"` | Yes. The only one. |
| Ice-shelf basal melt, `freshwater = "virtual"` (default) | No |
| Sea ice (`salt_flux`, `heat_added` through the coupler) | No. Virtual. |
| The seven atmospheric/river mass-flux components | No |

Under `"mass"` the melt is a Boussinesq volume source, `dh = m·dt/ρ₀`, added to `h_layer(:,:,nz)` by
`ocean_cavity_mass_step` in the stage, at thermo cadence, just after the surface-flux apply.
- Nothing is added to the barotropic forcing. `derive_bt_from_layers` rebuilds
  `bt_eta = Σh − bt_H_ref` at the top of every stage, so the volume appears in `bt_eta` one stage later.
  `bt_H_ref` stays static (it is the ice-base-to-bed datum, not a sea level).
- `ms%p_top` is prescribed and does not respond. In a closed domain the extra volume is a real sea-level
  rise, which is why `volume_compensation = "uniform_open_ocean"` exists.
- Budgets: `multilayer_state_t%mass_src` is the tracked source, accumulated with the same per-stage weight
  as `mass_out`, so the console residual `(M − M₀) + mass_out − mass_src` stays at round-off while the
  total grows. Salt and heat ride the existing `salt_budget_surface` / `heat_budget_surface`.

## Adding a seam: the short checklist

1. Decide which existing seam your quantity belongs to (table above) before adding a new one. The test is
   usually "what is it derived from, and who rebuilds it".
2. Write the contract in the `!!` docstring of the owning module, not here and not in a README.
3. Add a test that fails when the contract is broken (an analytical or bit-identity test). The existing
   gates are named in each owner's docstring.
4. If your change touches a row of the table above, update that owner in the same PR.
