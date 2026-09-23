"""Curated manifest of fast, physics-diverse OCEAN regression cases (P0).

This is the single source of truth for *which* bathymetry-free ocean namelist
inputs the regression runner exercises, and *how short* each one is run. It is
deliberately a hand-curated list (not a directory glob) so inclusion and the
per-case run-length / timeout knobs are explicit and reviewable.

Scope (P0): ocean only, run-clean / NaN gate only. Golden
comparison (P1) and gcov coverage (P2) are separate phases.

Each CASES entry is a dict:
    name        short unique id for the case (used for the scratch dir + report)
    nml         repo-relative path to the committed namelist
    n_steps     number of OUTER timesteps to run (the runner rewrites &time_nml
                so t_end = n_steps * dt_fixed seconds; keep this SMALL)
    timeout_s   hard per-case wallclock cap (catches hangs / runaway substeps)
    backends    set of backends this case is valid on -- {"cpu", "gpu"}
    tags        closures / schemes this case is meant to exercise (documentation
                + future coverage cross-check; not load-bearing for the P0 gate)
    tol         (optional, P1) per-case relative tolerance for the golden-summary
                compare (compare.py), overriding the loose global rtol. Present
                ONLY on cases whose measured CPU-vs-GPU spread genuinely exceeds
                the global bound for a well-understood, benign reason (documented
                inline) -- per the plan's "flag it, don't hide it under the
                global tol" rule.
    atol        (optional, P1) per-FIELD absolute tolerance for the compare,
                {field: value} -- `field` is a [diag] name (covers its min /
                max / mean) or "stats:<key>" (e.g. "stats:En"). There is NO
                global absolute floor: a field not named here is held to rtol
                alone, so a field that is PHYSICALLY ZERO and carries roundoff
                noise (a resting case's u / SSH / KE / En) must be named here,
                with the floor sized to that noise (not to the tolerance you
                would like). compare.py fails loud on a key that names no
                golden field, and reports every value that passed ONLY on its
                atol. Requires:
    atol_reason non-empty string: why each named field is legitimately ~0 and
                how its floor was sized.

Every listed nml is formula-bathymetry (no bathymetry_file/topo_file/dem_file/
input_file load) and small enough to finish O(seconds) at n_steps outer steps.
bench_scaling/* (perf benchmarks) are intentionally excluded.
"""

# All backends every case can currently run on. Kept as a module constant so a
# case that must be restricted (e.g. GPU-only, or CPU-only) overrides locally.
ALL_BACKENDS = {"cpu", "gpu"}

# Explicit absolute floors for the dynamic fields of a case whose exact
# answer is REST (see the `atol` key above). SI units: SSH m, u/v m/s, KE and
# En m^2/s^2. Sized from the measured noise of the resting seamount cases
# (SSH <= 1.05e-11 m, |u| <= 1.15e-11 m/s, KE <= 3.8e-22, En <= 2.1e-24):
# >= 100x above it, so a roundoff change (a new default, FMA, CPU vs GPU)
# passes, while a spin-up to 1 nm of SSH or 1 nm/s of current fails. KE/En
# use the velocity floor squared. Each user names only the fields its golden
# carries.
RESTING_NOISE_ATOL = {"SSH": 1e-9, "u": 1e-9, "v": 1e-9, "KE": 1e-18,
                      "stats:En": 1e-18}
RESTING_NOISE_REASON = (
    "exact answer is rest: SSH/u/v/KE/En are roundoff noise (SSH ~1e-11 m, "
    "u ~1e-13..1e-11 m/s), whose relative drift is meaningless; floors are "
    ">=100x the measured noise so a spin-up still fails")


CASES = [
    {
        # Canonical MOM6 double-gyre reference (NK=2): split-explicit RK2,
        # gprime reduced-gravity PGF, Sadourny PV Coriolis, linear bottom
        # drag, z* vcoord. The core dyn-core smoke.
        "name": "double_gyre_mom6",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_mom6.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pgf_gprime", "coriolis_sadourny", "bdrag_linear",
                 "vcoord_zstar", "pred_corr"],
    },
    {
        # Multi-layer (NK=10) double-gyre: FV-lite PGF reading per-layer rho,
        # Sadourny-energy Coriolis, Smagorinsky lateral viscosity, linear
        # density-range IC. Exercises the multi-layer FV pressure path.
        "name": "double_gyre_linear_nk10",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_linear_nk10.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pgf_fv_lite", "coriolis_sadourny_energy", "hvisc_smagorinsky",
                 "ic_linear_rho", "vcoord_zstar"],
    },
    {
        # PR-15 file-backed surface forcing: the double-gyre reference driven
        # from a time-varying NetCDF instead of the analytic 2gyre formula.
        # The wind file is GENERATED into the scratch dir by "setup" (it is
        # derived data, not committed), which is why this is the first case
        # with an input-file dependency -- see the module docstring's
        # formula-bathymetry note; the exception is the forcing file, not the
        # bathymetry, which is still formula "spoon".
        #
        # --steady writes the analytic 2gyre profile at every record, so this
        # case's answer must track double_gyre_mom6: it is a same-suite A/B on
        # the forcing PATH with the forcing VALUES held fixed.
        "name": "double_gyre_dataovr",
        "nml": "validation_examples/ocean/data_forcing/double_gyre_dataovr.nml",
        "setup": ["{python}", "{repo}/tools/make_forcing_nc.py", "wind.nc",
                  "--nx", "44", "--ny", "40", "--nt", "2", "--steady"],
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["dataovr_file_forcing", "data_input_cyclic", "pgf_gprime",
                 "coriolis_sadourny", "vcoord_zstar", "pred_corr"],
    },
    {
        # Stratified seamount: seamount topography, nonlinear EOS, quiescent
        # start -> spurious-pressure-gradient / z*-full remap stress test.
        "name": "seamount",
        "nml": "validation_examples/ocean/seamount/seamount.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["seamount_topo", "eos_nonlinear", "pgf", "ale_remap"],
        # Uniform T/S, f = 0, no forcing: the exact answer is rest, and every
        # dynamic field is roundoff noise -- SSH ~1e-11 m, u/v ~1e-13 m/s,
        # KE ~1e-25, En ~1e-27 (gfortran, 10 steps). Relative drift of noise is
        # meaningless: bebt 0 -> 0.1 alone moves SSH:max 1.05e-11 -> 3.2e-12
        # (-70 %). The floors sit >= 100x above the largest noise measured
        # (SSH 1.05e-11, u 1.15e-11) so a noise change passes, but a spin-up
        # past 1 nm / 1 nm/s fails. KE/En floors are the same velocity squared.
        "atol": RESTING_NOISE_ATOL,
        "atol_reason": RESTING_NOISE_REASON,
    },
    {
        # Barotropic geostrophic adjustment (NK=1): fast-mode / free-surface
        # + Coriolis geostrophic balance. Single-layer barotropic path.
        "name": "geostrophic_adjustment",
        "nml": "validation_examples/ocean/geostrophic_adjustment/geostrophic_adjustment.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["barotropic", "coriolis", "free_surface", "nk1"],
    },
    {
        # Island at rest: interior land masking (free-slip walls via
        # metric-zeroing) with a quiescent start -> should stay motionless.
        "name": "island_at_rest",
        "nml": "validation_examples/ocean/island_at_rest/island_at_rest.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["land_mask", "island", "quiescent", "wall_bc"],
        # Rest: u/v/KE/En are exactly 0.0 on gfortran. The floor keeps a
        # roundoff-level nonzero (another compiler, FMA) from failing a
        # bit-exact 0 while a spin-up past 1 nm/s still fails. SSH is NOT
        # noise here (max 6e-4 m) and stays on rtol.
        "atol": {k: RESTING_NOISE_ATOL[k]
                 for k in ("u", "v", "KE", "stats:En")},
        "atol_reason": RESTING_NOISE_REASON,
    },
    {
        # Coriolis + coastline wall: rotation against a solid boundary, wall BC
        # + Coriolis interaction.
        "name": "coriolis_coast",
        "nml": "validation_examples/ocean/coriolis_coast/coriolis_coast.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["coriolis", "wall_bc", "coastline_geometry"],
    },
    {
        # EPBL energetics boundary-layer mixing (Reichl-Hallberg) on a small
        # basin -> the EPBL vertical-mixing path.
        "name": "epbl_basin",
        "nml": "validation_examples/ocean/epbl_mld/epbl_basin.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["vmix_epbl", "mld", "surface_flux"],
    },
    {
        # KPP boundary-layer mixing on the same basin -> the KPP vmix path
        # (mutually exclusive with EPBL; both worth exercising).
        "name": "kpp_basin",
        "nml": "validation_examples/ocean/epbl_mld/kpp_basin.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["vmix_kpp", "mld", "surface_flux"],
    },
    {
        # Ideal-age passive tracer: the tracer registry + passive-transport
        # path (advection / remap / vertical exchange over a registered tracer).
        "name": "ideal_age_demo",
        "nml": "validation_examples/ocean/ideal_age/ideal_age_demo.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["passive_tracer", "ideal_age", "tracer_registry"],
    },
    {
        # Equilibrium body-tide forcing + scalar SAL in a basin -> the tidal
        # body-forcing path off the astronomical generator.
        "name": "body_tide_basin",
        "nml": "validation_examples/ocean/tides/body_tide_basin.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["tides_body", "sal", "astro_forcing"],
    },
    {
        # Mesoscale eddy spin-down (NK=15). Larger grid (300x300) -> kept to a
        # small step count to stay fast; exercises baroclinic eddy dynamics.
        "name": "eddy_test",
        "nml": "validation_examples/ocean/eddy_test/eddy_test.nml",
        "n_steps": 6,
        "timeout_s": 240,
        "backends": ALL_BACKENDS,
        "tags": ["baroclinic_eddy", "multilayer", "hvisc"],
    },
    {
        # Baroclinic-instability channel (NK=2): the classic re-entrant
        # channel baroclinic-instability setup.
        "name": "baroclinic_2layer",
        "nml": "validation_examples/ocean/baroclinic_channel/baroclinic_2layer.nml",
        "n_steps": 10,
        "timeout_s": 180,
        "backends": ALL_BACKENDS,
        "tags": ["baroclinic_instability", "channel", "periodic_bc"],
        # P1 per-case tol: early-transient baroclinic instability. At 10 steps
        # the along-channel jet u (~0.08 m/s) agrees CPU<->GPU tightly, but the
        # cross-channel v (~9e-4) and free-surface SSH (~4e-5) are the linear
        # growth of the perturbation -- 2-3 orders smaller, sitting near the
        # noise floor -- and their min/max EXTREMA diverge up to ~15% on
        # CPU(gfortran) vs GPU(nvfortran) from FMA / reduction ordering. The
        # integrated En/Mass/Salt/Temp all agree within the global 1e-3, so this
        # is benign transient-extrema divergence, not a physics regression. 0.25
        # covers the measured 0.146 with margin (GPU is deterministic run-to-run)
        # while still catching a real break (which moves u/En/Mass or blows the
        # extrema by orders).
        "tol": 0.25,
    },
    {
        # Eady problem: linear baroclinic (Eady) growth in a small box ->
        # stratified shear / thermal-wind balance.
        "name": "eady",
        "nml": "validation_examples/ocean/eady/eady.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["eady", "baroclinic", "thermal_wind", "multilayer"],
        # P1 per-case tol: same benign mechanism as baroclinic_2layer. The Eady
        # jet u (~0.145 m/s since the 2026-09-13 re-baseline to dT_dy = -2e-5)
        # agrees CPU<->GPU tightly; the tiny cross-channel v (~1e-4) and SSH
        # extrema are noise-floor quantities at 10 steps and diverge at the
        # 10% level while En/Mass/Salt/Temp agree within 1e-3. 0.25 covers it
        # with margin (deterministic run-to-run).
        "tol": 0.25,
    },
    {
        # Sponge-boundary relaxation demo -> the SPONGE BC nudging kernel.
        "name": "sponge_real_demo",
        "nml": "validation_examples/ocean/sponge_demo/sponge_real_demo.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["sponge_bc", "nudging", "open_boundary"],
    },
    {
        # ALE benchmark (NK=75): the PPM ALE remap path exercised across many
        # layers (MOM6 benchmark analogue).
        "name": "benchmark_ale",
        "nml": "validation_examples/ocean/benchmark_ale/benchmark_ale.nml",
        "n_steps": 8,
        "timeout_s": 240,
        "backends": ALL_BACKENDS,
        "tags": ["ale_remap_ppm", "many_layers", "vcoord"],
    },
    {
        # Flow past an island: wake formation around an interior land mass ->
        # land masking under a driven flow.
        "name": "flow_past_island",
        "nml": "validation_examples/ocean/flow_past_island/flow_past_island.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["land_mask", "island_wake", "driven_flow"],
    },
    {
        # Polar freeze-up (dynamics): the sea-ice EVP dynamics + coupling path
        # -> exercises the ice model, not just the ocean dyn-core.
        "name": "polar_freezeup_dynamics",
        "nml": "validation_examples/ocean/sea_ice/polar_freezeup_dynamics.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["sea_ice_evp", "ice_ocean_coupling", "ice_dynamics"],
    },

    # ---- P4 gap-filling cases (seamount / spoon formula topo) -------------
    # Added to exercise closures the P0-P2 corpus left at 0-10% coverage.
    # Each is built on FORMULA bathymetry (seamount/Lagrangian) with a
    # stratified, eddying state so the target closures do real work, and is
    # kept small + short to stay well within the runtime budget.
    {
        # Mesoscale lateral-closure stack on a stratified seamount: isopycnal
        # slopes + Gent-McWilliams thickness diffusion + Redi neutral tracer
        # diffusion + MEKE prognostic eddy energy, all wired on together
        # (MEKE requires GM requires slopes). One case lights up all four
        # mesoscale files that the base corpus never touched.
        "name": "seamount_gm_redi_meke",
        "nml": "validation_examples/ocean/seamount/seamount_gm_redi_meke.nml",
        "n_steps": 8,
        "timeout_s": 180,
        "backends": ALL_BACKENDS,
        "tags": ["gm", "redi", "meke", "isopycnal_slopes", "mesoscale",
                 "seamount_topo", "eos_nonlinear"],
    },
    {
        # St-Laurent/Simmons internal-tide bottom-intensified diapycnal
        # mixing on a stratified seamount (e_uniform > 0 -> non-zero Kd every
        # step). Exercises rdb_ocean_tidal_mixing (interior-closure path).
        "name": "seamount_tidal_mixing",
        "nml": "validation_examples/ocean/seamount/seamount_tidal_mixing.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["tidal_mixing", "vmix_closure", "seamount_topo"],
    },
    {
        # FV_MOM6 pressure-gradient with in-layer PLM T/S reconstruction of
        # the density integral (reconstruct_for_pressure). Exercises
        # rdb_ocean_pgf_reconstruct over the seamount's vertical density gradient.
        "name": "seamount_pgf_reconstruct",
        "nml": "validation_examples/ocean/seamount/seamount_pgf_reconstruct.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pgf_reconstruct", "pgf_fv_mom6", "seamount_topo"],
    },
    {
        # Lagrangian (isopycnal) seamount with the conservative minimum-
        # thickness borrow. The uniform-z seed grounds bottom layers on the
        # bump; conservative_floor borrows the sub-floor deficit from surplus
        # layers instead of injecting mass. Exercises rdb_ocean_min_thickness.
        "name": "seamount_conservative_floor",
        "nml": "validation_examples/ocean/seamount/seamount_conservative_floor.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["min_thickness", "conservative_floor", "lagrangian",
                 "isopycnal", "seamount_topo"],
        # Resting seamount again: En ~8.4e-27 is roundoff (its golden carries
        # no [diag] fields, only En and Mass). Mass stays on rtol.
        "atol": {"stats:En": RESTING_NOISE_ATOL["stats:En"]},
        "atol_reason": RESTING_NOISE_REASON,
    },
    {
        # Double diffusion (salt fingering) on a stratified seamount: warm+
        # salty over cool+fresh (R_rho ~= 1.30, inside (1, strat_param_max=
        # 2.55)) fires the fingering branch at every interior interface, so
        # ks diverges from kt and drives salinity through vdiff's two-source
        # solve. Exercises rdb_ocean_vmix vmix_split_ddiff_impl in a real
        # integration (the unit test covers the exact kernel math).
        "name": "seamount_ddiff",
        "nml": "validation_examples/ocean/seamount/seamount_ddiff.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["double_diffusion", "salt_fingering", "vmix_closure",
                 "seamount_topo"],
    },
    {
        # The same fingering column with split_scheme="pred_corr" PINNED:
        # double diffusion lives in the shared vmix_apply_in_stage, applied
        # in the pred_corr corrector.  Since the 2026-09-14 default flip the
        # sibling above runs the same scheme, so what this pair now proves is
        # that the pin and the default agree -- and the pin is what holds
        # this golden still if the default ever moves again.
        "name": "seamount_ddiff_pred_corr",
        "nml": "validation_examples/ocean/seamount/seamount_ddiff_pred_corr.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["double_diffusion", "salt_fingering", "pred_corr",
                 "vmix_closure", "seamount_topo"],
    },

    # ---- P4 gap-filling cases, batch 2 -----------------------------------
    # Second targeted batch: raise the three most under-covered closure files
    # that batch-1 left low -- MEKE (backscatter + length-scale + BBL off),
    # PGF reconstruct (only the PLM branch covered), and the baroclinic OBC
    # kernels (a closed basin can never reach them). Each is FORMULA-bathy
    # (seamount), stratified + eddying so the target closures do real work,
    # and short enough to stay well within the runtime budget.
    {
        # MEKE backscatter + mixing-length weights + BBL drag, wired on top of
        # the GM/Redi/slopes chain. Enables the harmonic momentum backscatter
        # (negative viscosity, meke_ku_closure + meke_backscatter_apply), the
        # deformation/Rhines/grid mixing-length arms, and the resolved-|u_bed|
        # BBL drag -- the MEKE code batch-1's seamount_gm_redi_meke left cold.
        # Backscatter is fail-loud: needs a flow-aware harmonic closure whose
        # ah_face the hvisc consumes (lateral_closure="smagorinsky") AND a
        # non-zero biharmonic backstop (smag_ah + default smag_bi_const=0.06).
        "name": "seamount_meke_backscatter",
        "nml": "validation_examples/ocean/seamount/seamount_meke_backscatter.nml",
        "n_steps": 8,
        "timeout_s": 180,
        "backends": ALL_BACKENDS,
        "tags": ["meke", "meke_backscatter", "meke_length_scales", "bbl_drag",
                 "hvisc_smagorinsky", "smag_ah", "gm", "redi", "seamount_topo"],
    },
    {
        # FV_MOM6 PGF with in-layer PPM (recon_scheme=2) T/S reconstruction --
        # the branch batch-1's seamount_pgf_reconstruct (PLM, recon_scheme=1)
        # left cold. Covers ppm_edges_column (implicit-h4 edges + CW limiter,
        # the module's largest uncovered block) and the parabolic branch of
        # the Boole density-integral quadrature (boole_dpa_intz_layer).
        "name": "seamount_pgf_ppm",
        "nml": "validation_examples/ocean/seamount/seamount_pgf_ppm.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pgf_reconstruct", "pgf_ppm", "pgf_fv_mom6", "seamount_topo"],
    },
    {
        # Baroclinic open-boundary kernels on an M2-tidal-west / Flather-open-
        # east stratified seamount: the Orlanski per-layer radiation branch
        # (radiation_scheme="orlanski") with Marchesiello nudging + the
        # implicit open-edge tracer reservoirs (res_lscale_out/in > 0). These
        # rdb_ocean_obc_baroclinic paths (~6% before) are unreachable in a
        # closed basin -- an OPEN boundary + NK>1 stratification is required.
        "name": "seamount_obc_baroclinic",
        "nml": "validation_examples/ocean/seamount/seamount_obc_baroclinic.nml",
        "n_steps": 10,
        "timeout_s": 180,
        "backends": ALL_BACKENDS,
        "tags": ["obc_baroclinic", "orlanski", "obc_reservoir", "obc_nudging",
                 "open_boundary", "tidal_bc", "seamount_topo"],
    },

    # ---- WENO tracer-reconstruction coverage -----------------------------
    # The base corpus runs every tracer with the default `ppm` reconstruction,
    # leaving src/tracer/structured/rdb_recon_weno.F90 (the shipped
    # weno5/7/9 + WENO-Z swept-average face helpers) completely uncovered. The
    # ocean windowed tracer-advect drain reuses those face helpers, so turning
    # it on with a WENO face reconstruction lights the module up. Both cases
    # are built on the Eady front (eady.nml) — its IC gives temperature a REAL
    # horizontal gradient (dT/dy) from step 0, so WENO reconstructs a non-flat
    # field and does genuine work (NOT a quiescent/uniform-tracer flat recon).
    {
        # WENO5-Z windowed drain: dt_tracer_advect_ratio=2 turns the windowed
        # tracer-advect drain ON, tracer_recon="weno5" swaps its face
        # reconstruction from the CW parabola to the WENO5-Z swept-average
        # ladder. Covers weno5_face_swept + plm_face_swept (near-boundary rung
        # degradation) + recon_rung_for_face + the ml_advect kernel path.
        # nghost=3 already satisfies weno5. dt_therm_ratio=2 keeps the ALE
        # remap cadence an integer multiple of the advect window (required).
        "name": "eady_weno5",
        "nml": "validation_examples/ocean/eady/eady_weno5.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["tracer_weno", "weno5", "windowed_tracer_advect",
                 "eady", "baroclinic"],
        # P1 per-case tol: inherits the parent eady's benign CPU<->GPU
        # divergence — the tiny cross-channel v (~1.5e-4) and SSH (~5e-6)
        # EXTREMA drift ~10% at the noise floor (FMA / reduction order) while
        # En/Mass/Salt/Temp agree within the global 1e-3. Measured GPU drift
        # 0.098 (v:min); 0.25 covers it with margin, same as eady/baroclinic_2layer.
        "tol": 0.25,
    },
    {
        # WENO7-Z windowed drain: same base + knobs but tracer_recon="weno7"
        # with nghost bumped to 4 (weno7 needs the wider stencil). Reaches the
        # weno7 cubic-candidate branch (weno7_face_swept) that weno5 cannot,
        # covering the higher-order stencil the weno5 case leaves cold.
        "name": "eady_weno7",
        "nml": "validation_examples/ocean/eady/eady_weno7.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["tracer_weno", "weno7", "windowed_tracer_advect",
                 "eady", "baroclinic"],
        # P1 per-case tol: same benign eady mechanism (noise-floor v/SSH extrema
        # diverge CPU<->GPU; integrated quantities agree). weno7's sharper
        # reconstruction lifts the drift a touch — measured 0.153 (v:max) —
        # still well under 0.25.
        "tol": 0.25,
    },

    # ---- MOM6 predictor-corrector split scheme -----------------------------
    # These two PIN `&ocean_bt_nml split_scheme="pred_corr"` explicitly.
    # They were written when the predictor-corrector was the NON-default
    # scheme and existed to keep its production branches in
    # rdb_ocean_dyn.F90 (predictor stage, corrector on the time-mean
    # u_av/h_av, restore_state, the is_pc conditionals) out of the cold. The
    # default flipped on 2026-09-14, so the whole base corpus now exercises
    # those branches and the two cases are no longer load-bearing for THAT —
    # but they stay, for two reasons: they pin the scheme, which is what
    # keeps them (and their goldens) invariant across any future default
    # move; and each is a distinct configuration in its own right (a
    # NON-Lagrangian ALE vcoord under the pc loop, #395). The branches that
    # are cold in the base corpus now are the ssp_rk2 ones — covered by
    # `eady_weno5` / `eady_weno7`, which pin ssp_rk2 because the pred_corr
    # v1 envelope refuses their windowed tracer-advect configuration.
    {
        # Quiescent stratified-free seamount (zstar_sigma, NK=15) under the
        # MOM6 predictor-corrector: a validated rest-preservation guard — the
        # pc outer loop holds the rest state to machine precision (En~2.5e-21,
        # MaxCFL 0) on a non-Lagrangian ALE coord, exactly as ssp_rk2 does.
        "name": "seamount_pred_corr",
        "nml": "validation_examples/ocean/seamount/seamount_pred_corr.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pred_corr", "predictor_corrector", "split_scheme",
                 "seamount_topo", "vcoord_zstar_sigma"],
        # Rest, same as `seamount`. Its committed golden (u ~1.1e-11 m/s, En
        # ~2.1e-24) is already NOT what origin/main produces (u 3.8e-13, En
        # 3.7e-27, byte-identical to `seamount` now that pred_corr is the
        # default): a 97 % noise move that the old global atol passed without
        # a word. Under this floor it passes BY NAME; refreshing the golden is
        # the maintainer's call.
        "atol": RESTING_NOISE_ATOL,
        "atol_reason": RESTING_NOISE_REASON,
    },
    {
        # Active wind-driven double-gyre (fv_lite PGF, NK=10, zstar) under the
        # MOM6 predictor-corrector: the corrector runs on REAL (non-trivial)
        # tendencies, so its golden is a meaningful drift target for the
        # corrector-on-time-means path — unlike the quiescent seamount guard.
        "name": "double_gyre_pred_corr",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_pred_corr.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["pred_corr", "predictor_corrector", "split_scheme",
                 "pgf_fv_lite", "coriolis_sadourny_energy", "vcoord_zstar"],
    },
    {
        # WENO3 PV reconstruction (F1) on the Sadourny ENSTROPHY Coriolis path:
        # the wind-driven western boundary current is a sharp PV front, so the
        # weno3 branch of coriolis_adv_compute_tendencies_sadourny does real
        # work every stage. Exercises rdb_coriolis_adv weno3_recon in a live
        # integration (the unit tests cover the exact reconstruction math).
        "name": "double_gyre_weno3",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_weno3.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["weno_pv", "coriolis_sadourny", "pv_adv_weno3",
                 "pgf_fv_lite", "vcoord_zstar"],
    },
    {
        # The identical WENO3-on-enstrophy gyre under split_scheme="pred_corr":
        # the Coriolis-adv kernel runs every stage of both integrators, so this
        # confirms the weno3 branch fires clean under the predictor-corrector.
        "name": "double_gyre_weno3_pred_corr",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_weno3_pred_corr.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["weno_pv", "coriolis_sadourny", "pv_adv_weno3", "pred_corr",
                 "split_scheme", "vcoord_zstar"],
    },
    {
        # weno5 PV reconstruction (radius-3 stencil, nghost=4): exercises the
        # weno5_recon path of the Coriolis kernel on the wind-driven gyre.
        "name": "double_gyre_weno5",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_weno5.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["weno_pv", "coriolis_sadourny", "pv_adv_weno5", "nghost4",
                 "vcoord_zstar"],
    },
    {
        # weno7 PV reconstruction (radius-4 stencil, nghost=4): exercises the
        # weno7_recon path (Balsara-Shu smoothness) on the wind-driven gyre.
        "name": "double_gyre_weno7",
        "nml": "validation_examples/ocean/double_gyre/double_gyre_weno7.nml",
        "n_steps": 10,
        "timeout_s": 120,
        "backends": ALL_BACKENDS,
        "tags": ["weno_pv", "coriolis_sadourny", "pv_adv_weno7", "nghost4",
                 "vcoord_zstar"],
    },
]


# ---------------------------------------------------------------------------
# Two-tier stability suite (stability.py)
# ---------------------------------------------------------------------------
# `CASES` above is the SHORT-RUN golden-drift corpus: 6-10 outer steps per case
# (1.7 simulated hours at dt=600) compared against a committed golden. That has
# two structural blind spots -- it stops before most physics exists, and a
# golden captured from a broken run passes forever.
#
# `STABILITY_CASES` is the answer to both: every one of the 68 tracked ocean
# namelists, run long enough that its physics manifests, asserted against
# PHYSICS (finite / closed budgets / energy must not grow in an unforced case /
# bounded CFL / the case's own stated claims) rather than against a stored
# answer. It lives in its own module because it is large and independently
# reviewable; see tests/regression/stability_manifest.py and README.md.
from stability_manifest import STABILITY_CASES  # noqa: E402,F401
