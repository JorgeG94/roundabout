"""Two-tier STABILITY manifest -- every tracked ocean namelist, both tiers.

This is the source of truth for `stability.py` and is imported into
`manifest.py` as `STABILITY_CASES`, alongside the older short-run `CASES`
list that `compare.py` / `coverage.py` still use.

Why a second list rather than more keys on the first
====================================================
`CASES` answers "did 10 steps of this namelist drift from its golden". This
list answers a different question -- "does this namelist produce PHYSICS that
holds up over a run long enough for its physics to exist" -- and it covers all
68 tracked ocean namelists, not the 35 the golden suite curates. Keeping them
separate means neither gate has to compromise for the other.

Reading an entry
================
    name            unique id (also the scratch dir and report row)
    nml             repo-relative path to the COMMITTED namelist -- never a
                    generated or edited copy. A twin is built by OVERRIDING
                    keys at run time; the shipped file is never modified.
    setup           optional argv run in the scratch dir first (generated
                    input files); "{python}" / "{repo}" are substituted.
    physics         the assertion block (see below)
    tier1           full-scale spec: n_steps, timeout_s, samples, and
                    optional `overrides` (`t1_overrides`) -- the same run-time
                    key overrides a tier-2 twin uses, for a tier-1 leg that is
                    a documented VARIANT of a committed namelist (melt off, a
                    viscosity twin) rather than the file as shipped.
    tier2           downscaled-twin spec: the same, plus `overrides` (namelist
                    keys the twin changes) and `dimensionless` (the inputs
                    downscale.py checks). `skip` with a `reason` marks a case
                    that CANNOT be honestly downscaled to CI size.
    known_failure   this case fails TODAY, for a known and documented reason.
                    Reported as XFAIL, with the reason printed, so CI red is
                    informative rather than flaky -- and so a reader can tell
                    a known defect from a regression their PR just introduced.

The `physics` block
===================
    regime          rest       -- starts motionless, no forcing: must STAY
                                  motionless. Any energy is spurious.
                    adiabatic  -- energetic IC, nothing driving it: energy may
                                  decay, MUST NOT GROW.
                    baroclinic -- an instability case: energy SHOULD grow.
                                  Gated on boundedness, plus the analytical
                                  growth rate at tier 1 where it is
                                  measurable.
                    forced     -- wind / heat / tide driven: spin-up growth is
                                  expected; gated on a physical ceiling and,
                                  at tier 1, on saturating.
    en_rest_max     rest only: the peak mean specific KE a "motionless" run may
                    reach, m2/s2. Stated as a VELOCITY bar: 5e-7 is 1 mm/s of
                    spurious current, 5e-13 is 1 um/s (a case that should be
                    bit-zero). Chosen from physics, never from what the run
                    happens to do.
    en_growth_max   adiabatic only: peak En as a multiple of its early-run
                    reference. 1.5 allows a transient adjustment; the
                    2026-09-11 thermo-cadence instability ran to 30-60x.
    en_max          baroclinic/forced: absolute ceiling, m2/s2. 0.5 = 1 m/s rms
                    over the whole domain, which no case here should approach.
    cfl_max         advective CFL ceiling.
    budget_tol      {"Mass"|"Salt"|"Heat": relative tolerance} on the MODEL'S
                    OWN closed-budget residual. These are printed every step;
                    they sit at 1e-14/1e-15 in a healthy run.
    claims          the namelist header's own testable statements. NEVER
                    weakened to pass, NEVER satisfied by editing the namelist:
                    when a header and the code disagree, that IS the finding.
                    Kinds: mass_rel, en_ratio, field_ratio (a [diag] field's
                    extrema amplification from `from_day` to the end), and
                    the day-anchored en_at_day (En(day) < max), en_day_ratio
                    (En(num_day)/En(den_day) < max) and en_rate_decel (the
                    log-rate over the `late` window below the `early` one).
                    A
                    claim that needs the full-length run names its tiers
                    (`"tiers": [1]`) and is SKIPped, not failed, elsewhere.

Stdlib only.
"""

# --- shared defaults --------------------------------------------------------
# The model's own budget residuals sit at 1e-14..1e-15 in a healthy run, so
# these bounds are ~100x roundoff: loose enough not to false-trip on reduction
# order, tight enough that a real leak (which shows up at 1e-6 and worse) is
# caught on the first step it happens.
BUDGET_EXACT = {"Mass": 1e-11, "Salt": 1e-11, "Heat": 1e-11}
# Cases with an open boundary, a sponge, or a surface flux exchange mass/salt/
# heat with the outside on purpose. The model still reports a closed residual
# (it accounts for `out` and `src`), so the same bound applies -- but the
# arithmetic has more terms, hence one extra decade.
BUDGET_OPEN = {"Mass": 1e-9, "Salt": 1e-9, "Heat": 1e-9}

# Velocity bars for a "must stay at rest" case, expressed as mean specific KE
# (En = 0.5 * u_rms^2), chosen from what the configuration can PHYSICALLY
# generate -- never from what a run happens to produce.
#
#   REST_SEAMOUNT  a terrain-following / hybrid coordinate over a tall
#                  seamount has a known, quantified spurious pressure
#                  gradient. The seamount test case's conventional acceptance
#                  bar is ~1 cm/s of spurious velocity (Beckmann & Haidvogel
#                  1993 and the sigma-coordinate PGF-error literature that
#                  follows it). MEASURED on this tree for reference, so any
#                  drift toward the bar is visible in review rather than
#                  hidden by it: the 96x96x20 zstar_sigma family sits at
#                  1.08 mm/s (9x under the bar) and `seamount` /
#                  `seamount_auto` / `seamount_pred_corr` at well under that;
#                  `seamount_flat` is 1e-23, i.e. exactly zero.
#   REST_1MM_S     a sloped case with no excuse for more than a mm/s.
#   REST_1UM_S     a case with NO slope and no forcing: nothing can generate
#                  a pressure-gradient error, so it must be bit-zero.
REST_SEAMOUNT = 0.5e-4  # 1 cm/s -- the seamount test's conventional bar
REST_1MM_S = 0.5e-6     # 1 mm/s of spurious current
REST_100UM_S = 0.5e-8   # 0.1 mm/s -- a resting case with an internal-wave
                        # field: a T-noise seed excites real IGWs, so exact
                        # zero is not attainable, but nothing may FEED them
REST_1UM_S = 0.5e-12    # 1 um/s -- a case that should be bit-zero

EN_CEIL = 0.5           # 1 m/s domain-rms: no case here should approach it
EN_CEIL_HI = 2.0        # energetic configurations (neverworld2, ACC)


def _case(name, nml, regime, t1_steps, t2_steps, **kw):
    """Build one entry, filling the parts that are the same for every case."""
    phys = {
        "regime": regime,
        "cfl_max": kw.pop("cfl_max", 0.9),
        "budget_tol": kw.pop("budget_tol", BUDGET_EXACT),
    }
    if regime == "rest":
        phys["en_rest_max"] = kw.pop("en_rest_max", REST_1MM_S)
        # OPT-IN, per row. `rest_sigma_max` arms the fitted GROWTH-RATE gate
        # and `matrix` arms the tracer-bound / thickness / counter gates that
        # the vertical-coordinate matrix carries (`assert_matrix`). They are
        # not switched on corpus-wide because several existing rest cases
        # have scoped `known_failure`s that would not cover a NEW assertion,
        # and an uncovered failure on a known-failing case turns the row red
        # for the wrong reason. Extending them is a separate change with its
        # own measurement.
        for k in ("rest_sigma_max", "matrix_gates", "tracer_slack",
                  "h_min_floor", "rest_trend_skip_tiers",
                  "rest_trend_skip_reason"):
            if k in kw:
                phys[k] = kw.pop(k)
    elif regime == "adiabatic":
        phys["en_growth_max"] = kw.pop("en_growth_max", 1.5)
    else:
        phys["en_max"] = kw.pop("en_max", EN_CEIL)
        if regime == "forced":
            phys["late_growth_max"] = kw.pop("late_growth_max", 4.0)
        if regime == "baroclinic":
            # An unforced instability case may grow at its PHYSICAL rate and
            # no faster. Where the case has an analytical band
            # (`sigma_expected`) the bar is derived from it; otherwise this
            # default applies. 3e-6 1/s is a 3.9-day e-folding:
            #   * the fastest PHYSICAL windowed rate measured anywhere in this
            #     corpus is 3.5e-7 (bc_inst_tuned_512, a 33-day e-folding), so
            #     the bar sits ~9x above real balanced growth;
            #   * the 2026-09-11 thermo-cadence instability ran at ~1.0e-5
            #     (27-hour e-folding), 3.3x ABOVE the bar.
            # The gap between balanced baroclinic growth and a grid-scale
            # numerical mode is one to two decades, which is what makes a
            # single bar workable here.
            phys["sigma_fast_max"] = kw.pop("sigma_fast_max", 3.0e-6)
    for k in ("sigma_expected", "claims"):
        if k in kw:
            phys[k] = kw.pop(k)

    entry = {
        "name": name,
        "nml": nml,
        "physics": phys,
        "tier1": {"n_steps": t1_steps,
                  "timeout_s": kw.pop("t1_timeout", 900),
                  "samples": kw.pop("t1_samples", 80)},
    }
    # [diag] emissions across the run (default 4: the [diag] line only carries
    # field extrema, and NetCDF is forced off, so more is only console). A
    # case whose growth-rate gate reads a FIELD's extrema (`sigma_expected`
    # with a "field") needs a dense series here -- eady reads max|v|.
    if "t1_diag_samples" in kw:
        entry["tier1"]["diag_samples"] = kw.pop("t1_diag_samples")
    if "t1_overrides" in kw:
        entry["tier1"]["overrides"] = kw.pop("t1_overrides")
    if kw.pop("t2_skip", False):
        entry["tier2"] = {"skip": True, "reason": kw.pop("t2_reason", "")}
    else:
        entry["tier2"] = {
            "n_steps": t2_steps,
            "timeout_s": kw.pop("t2_timeout", 180),
            "samples": kw.pop("t2_samples", 40),
            "overrides": kw.pop("t2_overrides", {}),
            "dimensionless": kw.pop("t2_dimensionless", {}),
        }
        if "t2_physics" in kw:
            entry["tier2"]["physics_override"] = kw.pop("t2_physics")
    for k in ("setup", "known_failure", "tags", "note", "matrix"):
        if k in kw:
            entry[k] = kw.pop(k)
    if kw:
        raise ValueError("{}: unknown keys {}".format(name, sorted(kw)))
    return entry


# ---------------------------------------------------------------------------
# Shared downscale recipes
# ---------------------------------------------------------------------------
# A periodic channel can be SHORTENED (fewer cells, SAME dx) without touching
# a single dimensionless number: the resolution, the deformation radius in
# cells, the Munk width in cells, the viscous CFL and the stratification are
# all unchanged; only the number of wavelengths the box holds goes down. That
# is the cleanest downscale available and is preferred wherever the geometry
# allows it. Rule R7 (>= 8 deformation radii across the box) is what stops it
# going too far.
#
# A CLOSED basin cannot be shortened that way -- the basin IS the case -- so
# there the only downscale is to coarsen (nx/2, dx*2), which moves Rd/dx and
# the Munk width in cells by the same factor of 2. For an eddy-resolving case
# that is fatal (rule R4), which is why several of them are tier-1 only.

SEAMOUNT_96_TO_48 = {
    # 96x96 @ 4 km -> 48x48 @ 8 km: the SAME 384 km domain and the same
    # seamount, at half the resolution. Legitimate here because the case is a
    # QUIESCENT spurious-pressure-gradient test -- there is no boundary
    # current to under-resolve (R1 n/a) and no eddy field to lose (R4 n/a) --
    # and because the repo already ships `seamount_conservative_floor` at
    # exactly 48x48 @ 8 km, so the coarse geometry is the author's own.
    "t2_overrides": {"grid_nml": {"nx": 48, "ny": 48, "dx": 8000.0, "dy": 8000.0}},
    "t2_dimensionless": {
        "dx": 8000.0, "dy": 8000.0, "dt": 300.0, "nu_h": 500.0,
        "nghost": 3, "has_western_boundary": False, "eddying": False,
    },
}

ACC_HALF_RES = {
    # 120x70 @ 0.2 deg -> 60x35 @ 0.4 deg. A re-entrant channel at -52 S:
    # periodic east-west, so there is no western boundary layer to resolve
    # (R1 n/a). Not declared eddying -- as shipped these ACC cases have no
    # meridional density gradient to draw on, so they produce no eddies at
    # ANY resolution (rule R5; see the acc_channel note below).
    "t2_overrides": {"grid_nml": {"nx": 60, "ny": 35, "dx": 0.4, "dy": 0.4}},
    "t2_dimensionless": {
        "dx": 0.4 * 111320.0, "dy": 0.4 * 111320.0, "dt": 600.0,
        "nu_h": 500.0, "ah_max": 1.0e4, "lat_deg": -52.0,
        "nghost": 3, "has_western_boundary": False, "eddying": False,
    },
}

_BENCH_REASON = (
    "performance benchmark: the problem SIZE is the thing being measured, so "
    "a downscaled twin is not the same test. Gated at tier 1 only.")


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------
V = "validation_examples/ocean/"

STABILITY_CASES = [

    # ===================== double gyre (wind-forced, closed basin) =========
    # 44x40, dt=1200. The canonical MOM6 reference and the dyn-core smoke.
    # Cheap enough that both tiers run the full shipped 10-day integration.
    _case("double_gyre_mom6", V + "double_gyre/double_gyre_mom6.nml",
          "forced", 720, 400,
          tags=["pgf_gprime", "coriolis_sadourny", "vcoord_zstar_full"]),
    _case("double_gyre_linear_nk10", V + "double_gyre/double_gyre_linear_nk10.nml",
          "forced", 720, 300, tags=["pgf_fv_lite", "nk10"]),
    _case("double_gyre_pred_corr", V + "double_gyre/double_gyre_pred_corr.nml",
          "forced", 720, 300, tags=["pred_corr"]),
    _case("double_gyre_weno3", V + "double_gyre/double_gyre_weno3.nml",
          "forced", 720, 300, tags=["pv_adv_weno3"]),
    _case("double_gyre_weno3_pred_corr", V + "double_gyre/double_gyre_weno3_pred_corr.nml",
          "forced", 720, 300, tags=["pv_adv_weno3", "pred_corr"]),
    _case("double_gyre_weno5", V + "double_gyre/double_gyre_weno5.nml",
          "forced", 720, 300, tags=["pv_adv_weno5"],
          t2_dimensionless={"nghost": 4, "schemes": ["weno5"]}),
    _case("double_gyre_weno7", V + "double_gyre/double_gyre_weno7.nml",
          "forced", 720, 300, tags=["pv_adv_weno7"],
          t2_dimensionless={"nghost": 4, "schemes": ["weno7"]}),
    _case("double_gyre_dataovr", V + "data_forcing/double_gyre_dataovr.nml",
          "forced", 720, 300,
          setup=["{python}", "{repo}/tools/make_forcing_nc.py", "wind.nc",
                 "--nx", "44", "--ny", "40", "--nt", "2", "--steady"],
          tags=["dataovr_file_forcing"]),

    # ===================== seamount family (QUIESCENT) =====================
    # Every one of these starts motionless over topography with no forcing.
    # The only thing that can move them is a discretisation error -- most
    # often the spurious pressure gradient a terrain-following coordinate
    # generates over a slope. So the gate is a VELOCITY BAR, not a golden:
    # a resting ocean must not develop a current above 1 mm/s.
    _case("seamount", V + "seamount/seamount.nml", "rest", 8640, 200,
          en_rest_max=REST_SEAMOUNT, tags=["seamount_topo", "eos_nonlinear"]),
    _case("seamount_auto", V + "seamount/seamount_auto.nml", "rest", 8640, 200,
          en_rest_max=REST_SEAMOUNT, tags=["auto_n_inner"]),
    _case("seamount_flat", V + "seamount/seamount_flat.nml", "rest", 8640, 200,
          # Flat bottom AND no forcing: with no slope there is no spurious
          # PGF to generate, so this one must be bit-zero, not merely small.
          en_rest_max=REST_1UM_S, tags=["flat_bottom", "rest_exact"]),
    _case("seamount_pred_corr", V + "seamount/seamount_pred_corr.nml", "rest", 4000, 200,
          en_rest_max=REST_SEAMOUNT, tags=["pred_corr", "vcoord_zstar_sigma"]),
    # Was a known failure: peak En 2.24e-3 => 6.7 cm/s of spurious current on a
    # motionless ocean, 45x over the seamount bar and still growing at the end
    # of the run. Root cause was the two-point PGF Jacobian on GROUNDED
    # isopycnal layers -- a layer squeezed onto the angstrom_h floor on the
    # shallow side of a face while massive on the deep side has its two layer
    # centres hundreds of metres apart in z, so the
    # `Delta p_centre + g*rho_layer*Delta z_centre` cancellation fails and
    # leaves g*(rho_layer - rho_ambient)*dz/dx of acceleration AT REST.
    # `&ocean_isopycnal_nml pgf_skip_nonoverlap` (default ON, VCOORD_LAGRANGIAN
    # only) zeroes the face PGF where the layers do not overlap in z. The case
    # now holds En ~ 1e-26 -- machine zero -- for 30 simulated days, so it gets
    # the REST_1UM_S "must be bit-zero" bar rather than the 1 cm/s seamount
    # tolerance: with the grounded-layer gradient gone there is nothing left
    # for this configuration to generate.
    _case("seamount_conservative_floor", V + "seamount/seamount_conservative_floor.nml",
          "rest", 576, 200, en_rest_max=REST_1UM_S,
          tags=["lagrangian", "bound_kh", "conservative_floor",
                "pgf_skip_nonoverlap", "rest_exact"]),
    _case("seamount_pgf_ppm", V + "seamount/seamount_pgf_ppm.nml", "rest", 576, 130,
          en_rest_max=REST_SEAMOUNT, tags=["pgf_ppm"], **SEAMOUNT_96_TO_48),
    _case("seamount_pgf_reconstruct", V + "seamount/seamount_pgf_reconstruct.nml",
          "rest", 576, 130, en_rest_max=REST_SEAMOUNT, tags=["pgf_reconstruct"],
          **SEAMOUNT_96_TO_48),
    _case("seamount_ddiff", V + "seamount/seamount_ddiff.nml", "rest", 576, 130,
          en_rest_max=REST_SEAMOUNT, tags=["double_diffusion"], **SEAMOUNT_96_TO_48),
    _case("seamount_ddiff_pred_corr", V + "seamount/seamount_ddiff_pred_corr.nml",
          "rest", 576, 130, en_rest_max=REST_SEAMOUNT,
          tags=["double_diffusion", "pred_corr"], **SEAMOUNT_96_TO_48),
    _case("seamount_gm_redi_meke", V + "seamount/seamount_gm_redi_meke.nml",
          "rest", 576, 130, en_rest_max=REST_SEAMOUNT,
          tags=["gm", "redi", "meke"], **SEAMOUNT_96_TO_48),
    _case("seamount_meke_backscatter", V + "seamount/seamount_meke_backscatter.nml",
          "rest", 576, 130, en_rest_max=REST_SEAMOUNT,
          tags=["meke_backscatter"], **SEAMOUNT_96_TO_48),
    _case("seamount_tidal_mixing", V + "seamount/seamount_tidal_mixing.nml",
          "rest", 576, 130, en_rest_max=REST_SEAMOUNT,
          tags=["tidal_mixing"], **SEAMOUNT_96_TO_48),

    # Seamount variants that ARE driven (open / tidal boundaries), so "rest"
    # would be the wrong assertion for them.
    _case("seamount_obc_baroclinic", V + "seamount/seamount_obc_baroclinic.nml",
          "forced", 1152, 200, budget_tol=BUDGET_OPEN,
          tags=["obc", "flather", "tidal_bc"]),
    _case("seamount_tidal_obc", V + "tides/seamount_tidal_obc.nml",
          "forced", 1152, 160, budget_tol=BUDGET_OPEN,
          tags=["obc", "tides", "nodal"]),
    _case("seamount_bench_full", V + "seamount/seamount_bench_full.nml",
          "rest", 300, 0, en_rest_max=REST_SEAMOUNT,
          t1_timeout=900, t2_skip=True, t2_reason=_BENCH_REASON,
          tags=["benchmark"]),

    # ===================== eady / baroclinic instability ===================
    # These exist to GROW. Asserting "energy must not grow" here would be
    # wrong -- and asserting a golden would enshrine whatever they currently
    # do, which is the trap `eady` fell into before 2026-09-13 (its golden WAS
    # the over-damped answer: dT_dy = -5e-6 + nu_h = 100 decayed at
    # sigma = -2.8e-7). Re-baselined that day to the front the header always
    # described (dT_dy = -2e-5, nu_h = 20, 60 days): the 125 km channel mode
    # now grows at 1.93e-6 1/s against 1.98e-6 from linear theory.
    #
    # The growth gate reads max|v| from the [diag] line, NOT En: En is the
    # TOTAL kinetic energy and is 4.3e-3 m2/s2 of basic-state jet that drifts
    # DOWN ~6% (wall drain) while the mode climbs three decades underneath;
    # it only turns up at day ~50. v is zero in the basic state, so its
    # extrema ARE the perturbation. The fit starts at day 30, by which time
    # the 125 km mode holds ~50% of max|v| (measured on the NetCDF modal
    # amplitudes); earlier the extrema are still the decaying noise seed.
    # `t1_diag_samples` makes the [diag] series dense enough to fit.
    _case("eady", V + "eady/eady.nml", "baroclinic", 8640, 240,
          en_max=EN_CEIL, t1_diag_samples=120, t1_timeout=1200,
          # Calibrated against the 2026-09-11 thermo-cadence instability
          # (sigma ~ 1.0e-5): 6e-6 catches it with 1.7x margin and leaves the
          # physical mode (1.9e-6) 3x of headroom. The default derivation
          # (5x the band's upper edge) would land at 1.25e-5 and miss it.
          sigma_fast_max=6.0e-6,
          sigma_expected={
              "min": 1.5e-6, "max": 2.5e-6, "field": "v", "fit_from_day": 30.0,
              "basis": "Eady (1949) linear theory for the kx = 2 (125 km) "
                       "mode of this 250 km channel at dT_dy = -2e-5 K/m: "
                       "Ri = 155, mu = Rd*sqrt(k^2 + (pi/L_y)^2) = 2.08, "
                       "sigma = 1.98e-6 1/s (unbounded-domain maximum "
                       "0.31 f/sqrt(Ri) = 2.49e-6)",
              "meaning":
                  "The case exists to demonstrate baroclinic instability at "
                  "the analytical rate. Measured 2026-09-13 on a V100: "
                  "sigma = 1.93e-6 from the console max|v| over days 30-60 "
                  "(1.99e-6 on the NetCDF kx = 2 modal amplitude, days "
                  "15-25). Too LOW means over-damping (nu_h back up, Smag_AH, "
                  "PPM truncation) or a fit window still inside the noise "
                  "floor; too HIGH, or growth at the same rate at every "
                  "zonal wavenumber, means a spurious source (the "
                  "resting-state growth documented in the nml header, an "
                  "EOS sign error, Coriolis-adv KE non-conservation).",
              "ref": "validation_examples/ocean/eady/eady.nml header; "
                     "tmp_local_artifacts/eady_hunt/p3"},
          claims=[{
              "kind": "field_ratio", "name": "header-1500x", "field": "v",
              "from_day": 1.0, "min": 300.0, "max": 1.0e4, "tiers": [1],
              "text": "the header advertises max|v| ~1500x over 60 days "
                      "(1.3e-4 -> 0.2 m/s)",
              "meaning":
                  "Measured 1480x (nvfortran, V100). The RNG realisation of "
                  "the noise seed differs between compilers, which moves the "
                  "day the mode clears the noise by a few days and the ratio "
                  "by a factor of ~2-3 either way; 300x is ~5 e-folds short "
                  "of the measurement and still 100x more than the old "
                  "over-damped file could ever show. Above 1e4 the run has "
                  "left the linear window (the header's day-60 limit).",
              "ref": "validation_examples/ocean/eady/eady.nml header"}],
          tags=["eady", "baroclinic"]),
    _case("eady_weno5", V + "eady/eady_weno5.nml", "baroclinic", 8640, 240,
          en_max=EN_CEIL, t1_diag_samples=120, t1_timeout=1200,
          sigma_fast_max=6.0e-6,
          t2_dimensionless={"nghost": 3, "schemes": ["weno5", "periodic"]},
          # Same basin, same re-baseline, same gates as `eady`: measured
          # 1.93e-6 / 1485x on this file (2026-09-13), i.e. the windowed WENO5
          # drain neither damps nor excites the mode.
          sigma_expected={
              "min": 1.5e-6, "max": 2.5e-6, "field": "v", "fit_from_day": 30.0,
              "basis": "as eady: kx = 2 channel mode, sigma = 1.98e-6 1/s",
              "meaning": "as eady; a rate that differs from plain eady's is "
                         "the windowed tracer-advect path (dt_tracer_advect_"
                         "ratio = 2, dt_therm_ratio = 2, weno5) acting on the "
                         "mode, which it must not.",
              "ref": "validation_examples/ocean/eady/eady.nml header"},
          claims=[{
              "kind": "field_ratio", "name": "header-1500x", "field": "v",
              "from_day": 1.0, "min": 300.0, "max": 1.0e4, "tiers": [1],
              "text": "max|v| ~1500x over 60 days, as eady",
              "meaning": "as eady.",
              "ref": "validation_examples/ocean/eady/eady.nml header"}],
          note=
              "reads as a weno5 reconstruction test but changes THREE knobs at once: "
              "tracer_recon=weno5 AND dt_tracer_advect_ratio=2 AND dt_therm_ratio=2. "
              "Bisected on 2026-09-11: the WENO reconstruction is innocent (clean "
              "for 25 days on its own); the thermo-cadence lag drove a grid-scale "
              "internal-wave instability that grew En 30x+ with every budget exact. "
              "This branch carries the fixes (08c7d3fe recompute the EOS every "
              "dynamics step, 19fee596 hold tracer CONCENTRATION across an advect "
              "window), which is what makes this case a regression guard for those "
              "two commits -- the energy:no-fast-growth bar is the one that sees it.",
          tags=["tracer_weno5", "windowed_advect", "dt_therm_ratio"]),
    _case("eady_weno7", V + "eady/eady_weno7.nml", "baroclinic", 8640, 240,
          en_max=EN_CEIL, t1_diag_samples=120, t1_timeout=1200,
          sigma_fast_max=6.0e-6,
          t2_dimensionless={"nghost": 4, "schemes": ["weno7", "periodic"]},
          # Measured 1.91e-6 / 1060x on this file (2026-09-13).
          sigma_expected={
              "min": 1.5e-6, "max": 2.5e-6, "field": "v", "fit_from_day": 30.0,
              "basis": "as eady: kx = 2 channel mode, sigma = 1.98e-6 1/s",
              "meaning": "as eady_weno5, for the weno7 rung.",
              "ref": "validation_examples/ocean/eady/eady.nml header"},
          claims=[{
              "kind": "field_ratio", "name": "header-1500x", "field": "v",
              "from_day": 1.0, "min": 300.0, "max": 1.0e4, "tiers": [1],
              "text": "max|v| ~1000-1500x over 60 days, as eady",
              "meaning": "as eady.",
              "ref": "validation_examples/ocean/eady/eady.nml header"}],
          note=
              "reads as a weno7 reconstruction test but changes THREE knobs at once: "
              "tracer_recon=weno7 AND dt_tracer_advect_ratio=2 AND dt_therm_ratio=2. "
              "Bisected on 2026-09-11: the WENO reconstruction is innocent (clean "
              "for 25 days on its own); the thermo-cadence lag drove a grid-scale "
              "internal-wave instability that grew En 30x+ with every budget exact. "
              "This branch carries the fixes (08c7d3fe recompute the EOS every "
              "dynamics step, 19fee596 hold tracer CONCENTRATION across an advect "
              "window), which is what makes this case a regression guard for those "
              "two commits -- the energy:no-fast-growth bar is the one that sees it.",
          tags=["tracer_weno7", "windowed_advect", "dt_therm_ratio"]),
    # ---- the OUTER-SPLIT rest-preservation gate --------------------------
    # `eady.nml` with the front off and nu_h = 0: a motionless, stably
    # stratified, flat-bottomed periodic channel seeded with +-0.5 mK of
    # white noise. No forcing, no slope, no front -- so it has NO energy
    # source and every joule of En it develops is manufactured by the
    # numerics. Nothing else in this corpus isolates the OUTER time-split
    # that way: the seamount family tests the pressure gradient over
    # topography, `island_at_rest` and `seamount_flat` are bit-zero (no
    # internal-wave field to excite at all), and every baroclinic case has a
    # real instability whose growth masks a numerical one.
    #
    # It is the case that separates the two outer schemes by four decades
    # (measured day 25, both compilers):
    #
    #     ssp_rk2      En 2.992E-05 (nvfortran) / 2.72e-05 (gfortran)
    #     pred_corr    En 1.739E-09 (nvfortran) / 1.76e-09 (gfortran)
    #
    # -- which is why `pred_corr` is the DEFAULT since 2026-09-14, and why the
    # `__ssp_rk2` twin the scheme axis builds is an XFAIL rather than a
    # deleted case: ssp_rk2 is still shipped and still supported, and this is
    # what "experimental" on it means, stated as a number.
    _case("resting_stratified_channel",
          V + "eady/resting_stratified_channel.nml", "rest", 3600, 1440,
          en_rest_max=REST_100UM_S, t1_timeout=900,
          # 0.1 mm/s. The pred_corr answer is 59 um/s and reproduces to 1%
          # across gfortran and nvfortran despite a DIFFERENT RNG
          # realisation of the noise seed, so the 2.9x margin is real
          # margin and not compiler luck. Set from the velocity a resting
          # ocean may not develop, not from the run: it sits 5700x BELOW
          # the ssp_rk2 answer, so it is the defect it has to catch that
          # fixes the decade, not the measurement.
          # The tier-2 twin is the FULL grid at the FULL dt, shortened to 10
          # days (1440 steps, ~60 s). Not a downscale -- 50x50x10 is already
          # CI-sized -- but it MUST stay long enough for the defect to be
          # visible: ssp_rk2 crosses the 0.5e-08 bar at day 7.5 and reaches
          # 4.8e-08 (10x over) by day 10, while pred_corr is flat at 1.7e-09
          # from day 1. At the 300 steps a naive twin would use (2 days) BOTH
          # schemes read ~1e-09 and the case gates nothing -- which is what it
          # did on the first run, XPASSing its own XFAIL. That matters MORE
          # now that ssp_rk2 is the twin rather than the base: the run length
          # is what makes the twin's XFAIL real.
          t2_overrides={"time_nml": {"dt_fixed": 600.0}},
          t2_dimensionless={"dx": 5000.0, "dy": 5000.0, "dt": 600.0,
                            "nu_h": 0.0, "nghost": 3, "schemes": ["periodic"],
                            "has_western_boundary": False, "eddying": False},
          # Under the DEFAULT (pred_corr) this case PASSES the magnitude
          # gate -- En 1.739E-09 against the 0.5e-08 bar -- and XFAILs only
          # the SETTLE gate, because the residual growth is slowed by ~36x
          # rather than removed (83-day e-folding). The `__ssp_rk2` twin the
          # scheme axis builds is the one that XFAILs BOTH: 2.992E-05 by day
          # 25, 6000x over the bar, still climbing on a 2.5-day e-folding.
          # That four-decade separation is what makes this file worth
          # shipping. Fixing either is a scheme change, not a tuning change
          # -- do NOT close it by raising en_rest_max or by putting nu_h back
          # into the file (nu_h = 100 hides it at En 1.6e-11 without touching
          # its rate).
          known_failure={
              "assertions": ["energy:rest-settles"],
              "reason":
                  "OPEN, and scoped: pred_corr (the DEFAULT) does not "
                  "ELIMINATE the resting-state growth, it slows it by ~36x. "
                  "Measured over 120 days on a V100, En goes 1.21e-09 -> "
                  "4.82e-09 -- an 83-day e-folding with no sign of "
                  "saturating, against ssp_rk2's 2.5-day one -- so the "
                  "SETTLE gate cannot pass at any horizon. The MAGNITUDE "
                  "gate (energy:rest) is deliberately NOT excused here: it "
                  "is what separates the schemes (1.739E-09 vs 2.992E-05) "
                  "and pred_corr passes it with margin. Do not close this by "
                  "raising en_rest_max or en_rest_trend_floor -- the fix is "
                  "whatever still feeds the internal-wave field.",
              "ref": "validation_examples/ocean/eady/"
                     "resting_stratified_channel.nml",
          },
          note="the regression test for the class of bug the 2026-09-13 "
               "split-scheme flip fixed: a spurious internal-gravity-wave "
               "instability of the OUTER time-split. Substitution exonerated "
               "the Coriolis form, the ALE remap and the PGF form (all within "
               "0.1% of the ssp_rk2 baseline); removing the stratification "
               "dropped En 119x; removing the seed made the run bit-zero.",
          tags=["rest", "outer_split", "internal_waves", "vcoord_zstar"]),

    _case("baroclinic_2layer", V + "baroclinic_channel/baroclinic_2layer.nml",
          "baroclinic", 17280, 300, en_max=EN_CEIL,
          # Re-entrant channel: SHORTEN the box, keep dx -- every dimensionless
          # number is preserved exactly (see the recipe note above).
          t2_overrides={"grid_nml": {"nx": 96, "ny": 48}},
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 200.0,
                            "nu_h": 100.0, "nghost": 3,
                            "schemes": ["periodic"],
                            "has_western_boundary": False, "eddying": False},
          tags=["baroclinic_channel", "nk2"]),
    _case("baroclinic_15layer", V + "baroclinic_channel/baroclinic_15layer.nml",
          # The FULL shipped 50-day integration, and it has to be: this case's
          # instability does not appear until day ~30. En decays for 28 days
          # (1.72e-2 -> 1.62e-2), turns at day 30 and then runs away to 3.44e-1
          # by day 50 -- a 20x amplification, MaxCFL 0.018 -> 0.207. Run for
          # 14 days it reports "DECAYS EVERYWHERE" and looks exactly like
          # eady's genuine over-damping. This is the clearest case in the
          # corpus for the rule that tier-1 length is set by WHEN THE PHYSICS
          # APPEARS, per case, and never by a global step count. ~1000 s on a
          # V100, which is why it is nightly and not CI.
          "baroclinic", 36000, 90, en_max=EN_CEIL_HI, t1_timeout=2400,
          t2_overrides={"grid_nml": {"nx": 128, "ny": 64}},
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 120.0,
                            "nu_h": 80.0, "nghost": 3,
                            "schemes": ["periodic"],
                            "has_western_boundary": False, "eddying": False},
          tags=["baroclinic_channel", "nk15"]),
    _case("bc_inst_tuned_512", V + "bc_inst/bc_inst_tuned_512.nml",
          "baroclinic", 10000, 200, en_max=EN_CEIL,
          t2_overrides={"grid_nml": {"nx": 128, "ny": 128}},
          t2_dimensionless={"dx": 0.0439453125 * 111320.0,
                            "dy": 0.0439453125 * 111320.0, "dt": 200.0,
                            "nu_h": 50.0, "lat_deg": 53.625, "nghost": 3,
                            "schemes": ["periodic"],
                            "has_western_boundary": False, "eddying": False},
          tags=["baroclinic_instability", "lagrangian"]),

    # ===================== ACC channel =====================================
    # Re-entrant channel at -52 S. NOTE (rule R5): as shipped these have a
    # surface heat flux but no meridional density gradient to convert into
    # eddy kinetic energy, so they are NOT eddy-resolving tests however long
    # they run -- they are gated as forced, bounded flows, which is what they
    # actually are.
    _case("acc_channel", V + "acc_channel/acc_channel.nml", "forced", 288, 120,
          en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          # nu_h was raised to 22000 on 2026-09-11 to satisfy the Munk
          # criterion, and ah_max with it (R2) -- raising one without the
          # other would have been inert.
          t2_dimensionless={"dx": 0.5 * 111320.0, "dy": 0.5 * 111320.0,
                            "dt": 300.0, "nu_h": 22000.0, "ah_max": 2.2e4,
                            "lat_deg": -52.0, "nghost": 3,
                            "schemes": ["periodic"],
                            "has_western_boundary": False},
          tags=["acc", "channel", "munk_sized"]),
    _case("acc_channel_quiescent", V + "acc_channel/acc_channel_quiescent.nml",
          "rest", 360, 200, en_rest_max=REST_1UM_S, budget_tol=BUDGET_OPEN,
          note="q_heat = 0 and no other driver: a zero-forcing channel must "
               "stay bit-zero. Weak by construction (see the aquaplanet "
               "finding: a quiescent uniform-density run stays at zero "
               "whether the numerics are right or wrong) -- kept because it "
               "is nearly free and catches a spurious source.",
          tags=["acc", "quiescent"]),
    _case("acc_channel_weno5", V + "acc_channel/acc_channel_weno5.nml",
          "forced", 288, 140, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          t2_overrides=dict(ACC_HALF_RES["t2_overrides"]),
          t2_dimensionless=dict(ACC_HALF_RES["t2_dimensionless"],
                                dt=300.0, schemes=["weno5", "periodic"]),
          tags=["acc", "tracer_weno5", "windowed_advect"]),
    _case("acc_channel_kitchensink", V + "acc_channel/acc_channel_kitchensink.nml",
          "forced", 144, 100, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          **{k: dict(v) for k, v in ACC_HALF_RES.items()},
          tags=["acc", "kitchensink", "all_closures"]),
    _case("acc_channel_sw_penetration", V + "acc_channel/acc_channel_sw_penetration.nml",
          "forced", 144, 100, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          **{k: dict(v) for k, v in ACC_HALF_RES.items()},
          tags=["acc", "shortwave_penetration"]),
    _case("acc_channel_eddy", V + "acc_channel/acc_channel_eddy.nml",
          "forced", 144, 0, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          t1_timeout=900, t2_skip=True,
          t2_reason="eddy-resolving at 0.1 deg: the only affordable downscale "
                    "is to coarsen, which halves Rd/dx and violates rule R4 "
                    "(>= 4 cells per deformation radius). Below ~2 cells no "
                    "eddies form at ANY run length, so the twin would not be "
                    "the same test. Tier 1 only.",
          tags=["acc", "eddying"]),
    _case("acc_channel_kitchensink_xl", V + "acc_channel/acc_channel_kitchensink_xl.nml",
          "forced", 144, 0, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          t1_timeout=900, t2_skip=True,
          t2_reason="this IS the XL twin of acc_channel_kitchensink, which is "
                    "already in tier 2 at half the resolution. Downscaling it "
                    "would duplicate that case exactly.",
          tags=["acc", "kitchensink", "xl"]),

    # ===================== boundary-layer mixing ===========================
    _case("epbl_basin", V + "epbl_mld/epbl_basin.nml", "forced", 2160, 300,
          budget_tol=BUDGET_OPEN, tags=["epbl", "mld"]),
    _case("epbl_lt_basin", V + "epbl_mld/epbl_lt_basin.nml", "forced", 2160, 300,
          budget_tol=BUDGET_OPEN, tags=["epbl", "langmuir"]),
    _case("kpp_basin", V + "epbl_mld/kpp_basin.nml", "forced", 2160, 300,
          budget_tol=BUDGET_OPEN, tags=["kpp", "mld"]),

    # ============ ice-shelf cavity at rest (the sloping-lid PGF) ===========
    # Three quiescent cavities that differ by ONE ingredient each, so the set
    # says WHICH part of the discretisation manufactured any energy. All
    # three: flat bed, no wind, no surface flux, no melt, no bottom drag, no
    # vertical mixing and NO lateral viscosity, so nothing damps the answer
    # into the floor and nothing can hide a pressure-gradient regression.
    # Full derivation, the substitution table and the Yung et al. (2026)
    # comparison live in validation_examples/ocean/ice_shelf_cavity/README.md.
    _case("cavity_flat_lid_rest",
          V + "ice_shelf_cavity/cavity_flat_lid_rest.nml",
          "rest", 4320, 1440, en_rest_max=REST_1UM_S,
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 600.0,
                            "nu_h": 0.0, "nghost": 2,
                            "has_western_boundary": True, "eddying": False},
          note="a FLAT ice lid over a flat bed makes every interface gap "
               "De(K) = 0, so every FV-MOM6 trapezoid error G(K) is "
               "structurally ABSENT and the only correct answer is exactly "
               "zero motion. It is the control for three things at once: the "
               "load cancellation pa(nz+1) = rho_ref*g*eta_geo + "
               "rho_ref*g*z_draft (the same product twice, opposite signs), "
               "the datum bt_H_ref = b - z_draft (which makes bt_eta = 0 at "
               "rest), and the geopotential zinit overlay. MEASURED "
               "En = 0.000E+00 at every daily sample out to 30 days, under "
               "BOTH outer schemes -- hence the bit-zero bar.",
          tags=["cavity", "ice_shelf", "rest_exact", "vcoord_sigma",
                "pgf_fv_mom6", "zinit_linear"]),
    _case("cavity_flat_lid_rest_zfixed",
          V + "ice_shelf_cavity/cavity_flat_lid_rest_zfixed.nml",
          "rest", 4320, 1440, en_rest_max=REST_1UM_S,
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 600.0,
                            "nu_h": 0.0, "nghost": 2,
                            "has_western_boundary": True, "eddying": False},
          note="the P6.2 gate: the flat-lid cavity on the QUASI-GEOPOTENTIAL "
               "coordinate (vcoord_type='z_fixed'), the first z-like family "
               "taught about the ice base. Under a UNIFORM draft every column "
               "vanishes the SAME layers to the inert filler and cuts its "
               "first live layer at the SAME depth, so every interface offset "
               "De(K) is identically zero -- fillers included -- and the only "
               "correct answer is again exactly zero motion. It separates 'the "
               "rigid-top target builder is RIGHT' from 'the rigid-top target "
               "builder happens to be small', which the sigma twin cannot: "
               "under sigma nothing vanishes at all. MEASURED En = 0.000E+00 "
               "at every daily sample to 30 days, mass and salt residuals "
               "-1.8E-12 / -2.0E-12 relative at day 30 and EXACTLY 0.000E+00 "
               "at step 0 -- the fillers are seeded at hTr = 0, so the first "
               "regrid's drain has nothing to discard and no budget step "
               "appears. h_min_cavity = 96 m = 2*h_nominal is the ISOMIP+ "
               "minimum-column rule (Asay-Davis et al. 2016 3.1.5) that "
               "z_fixed x cavity fails loud on; nothing is grounded at 40 or "
               "at 96 here, so the grounding line does not move between the "
               "legs. The SLOPING twin is deliberately NOT tracked: its "
               "outcome is toolchain-dependent (NaN during day 18 on "
               "gfortran; nvfortran/GPU saturates at En = 3.845E-04, four "
               "decades over the sigma leg) because the ice-base staircase "
               "needs Yung et al. (2026) 3.3.1/3.3.2/3.2, none of which are "
               "in this build. See "
               "validation_examples/ocean/ice_shelf_cavity/"
               "cavity_sloping_lid_rest_zfixed.nml for that measurement.",
          tags=["cavity", "ice_shelf", "rest_exact", "vcoord_z_fixed",
                "pgf_fv_mom6", "zinit_linear"]),
    _case("cavity_uniform_rho_rest",
          V + "ice_shelf_cavity/cavity_uniform_rho_rest.nml",
          "rest", 4320, 1440, en_rest_max=REST_1UM_S,
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 600.0,
                            "nu_h": 0.0, "nghost": 2,
                            "has_western_boundary": True, "eddying": False},
          note="the sloping lid and the calving front of the headline case, "
               "with the stratification removed (lin_ds_dz = 0). G(K) = "
               "-(De(K)^3/12)*rho_0*N^2, so at N^2 = 0 the truncation is not "
               "small but ABSENT however steep the lid -- whatever is left is "
               "the LOAD bookkeeping plus rounding. T_ref/S_ref sit ON the "
               "uniform T/S so rho == rho_0 and the rho_0*g*z_draft load is "
               "the exact displaced weight: under the MOM6 barotropic split "
               "(bc_pgf_forcing) a mismatch is a real bottom-pressure "
               "gradient (4 mm/s with rho - rho_0 = -0.29 kg/m3). MEASURED "
               "En = 2.017E-20 (pred_corr) / 1.303E-19 (ssp_rk2) at day 30: "
               "machine zero, "
               "and twelve decades under the stratified twin at the same day. "
               "That is what says the stratified case's residual really is "
               "the rho_0*N^2*De^3 term and not a mis-cancelled 5.26 MPa ice "
               "load pretending to be one.",
          tags=["cavity", "ice_shelf", "rest_exact", "uniform_density",
                "pgf_fv_mom6"]),
    _case("cavity_sloping_lid_rest",
          V + "ice_shelf_cavity/cavity_sloping_lid_rest.nml",
          "rest", 4320, 2880, en_rest_max=REST_1MM_S, t1_timeout=900,
          # THE BAR. REST_1MM_S is the sloped-rest family's existing bar and
          # the velocity a resting sub-shelf cavity has no excuse to exceed
          # (real sub-shelf flows are cm/s, so 1 mm/s of spurious current is
          # already a serious contaminant). Measured day-30 peak 2.298E-08 =>
          # 2.14e-04 m/s: 22x under it in energy, 4.7x in velocity (gfortran
          # 15.1, 2026-09-25, trimmed IC, MOM6 split).
          #
          # THE BALANCED IC. The namelist sets &ocean_cavity_dyn_nml
          # trim_ic_for_p_surf (MOM6 TRIM_IC_FOR_P_SURF): the rho_ref*g*z_draft
          # load is lighter than the stratified water it displaces by
          # g*int(rho - rho_ref), a depth-uniform force up to 1.8e-5 m/s^2 that
          # the MOM6 split (bc_pgf_forcing) hands the barotropic mode. Without
          # the trim the day-1 adjustment reads En 1.061E-06 -- over this bar;
          # with it, 1.117E-09. Do NOT turn the trim off, or bc_pgf_forcing,
          # to "fix" a failure here: both put the load shortfall back.
          #
          # It does NOT clear the tighter REST_100UM_S (0.5e-08) the Phase-5
          # design proposed -- day 30 is 4.6x over -- and that is recorded here
          # rather than legislated away: restoring REST_100UM_S is the
          # Phase-6 acceptance criterion. Do NOT close the gap by widening
          # this bar, by shortening the run, or by putting viscosity back
          # into the namelist.
          #
          # The tier-2 twin is the FULL grid at the FULL dt, shortened to 20
          # days (2880 steps, ~24 s). Not a downscale -- 48x6x15 is already
          # CI-sized -- but the LENGTH is load-bearing: ssp_rk2 tracks
          # pred_corr to within 2% for 11 days and only then leaves
          # (6.06E-08 at day 13, 4.24E-05 at day 20, while pred_corr is still
          # at 1.65E-09). At the 10 days a naive twin would use, BOTH schemes
          # read ~1.9E-09 and the twin's energy:rest XFAIL would XPASS
          # against a case that gates nothing -- the same trap
          # `resting_stratified_channel` documents. At 20 days the BASE case
          # ends on its plateau (final 1.645E-09 = 88% of the 1.877E-09
          # day-10 peak) and SETTLES, so its known_failure is tier-1 only.
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 600.0,
                            "nu_h": 0.0, "nghost": 2,
                            "has_western_boundary": True, "eddying": False},
          known_failure={
              "tiers": [1],
              "assertions": ["energy:rest-settles"],
              "reason":
                  "OPEN, scoped, and the POINT of the case: this build "
                  "carries none of the sloping-surface PGF corrections, so "
                  "the sigma truncation under a tilted ice base does not "
                  "merely hold a static spurious current -- it seeds a mode "
                  "that leaves the day-1..20 plateau (En 0.9-1.9E-09, the "
                  "derived N^2*D^3/(6*dx*Hbar) scale) on a ~3.2-day "
                  "e-folding. It is BOUNDED, not runaway: 2.298E-08 at day "
                  "30, 3.492E-06 at day 60 with the rate visibly decaying, "
                  "budgets exact to 8e-13 throughout, no CFL truncation and "
                  "no clamping. So the final sample is the peak at every "
                  "horizon past day ~22 and short of saturation, and "
                  "energy:rest-settles cannot pass at tier 1 (the 20-day "
                  "tier-2 twin ends on the plateau and settles). "
                  "Characterised by substitution (legacy split, whose curve "
                  "the trimmed default tracks to 30%): dt 600 -> 300 "
                  "reproduces the curve (3.069E-08 vs 3.043E-08), so it is "
                  "NOT an (omega*dt)^n outer-split mode; ISOMIP+'s own "
                  "nu_h = 6.0 and a dx^4-scaled biharmonic nu_4 = 1e7 each "
                  "only DELAY it ~10 days; and the flat-lid "
                  "(En = 0.000E+00) and uniform-density (1.34E-20) siblings "
                  "are bit-zero, so the sloping-lid PGF error is the whole "
                  "source. The MAGNITUDE gate energy:rest is deliberately "
                  "NOT excused -- it is what separates the schemes "
                  "(2.298E-08 vs ssp_rk2's 1.850E-04) and pred_corr passes "
                  "it with margin. The fix is Yung et al. (2026)'s three "
                  "corrections, not a tolerance.",
              "ref": "validation_examples/ocean/ice_shelf_cavity/README.md",
          },
          note="THE Phase-5 physics gate and the Phase-6 baseline: a "
               "stratified cavity under a linearly sloping ice shelf with a "
               "calving front, initialised at rest with the isopycnals flat "
               "in GEOPOTENTIAL z (&ocean_zinit_nml source='linear'; a "
               "layer-index profile would tilt them with the sigma "
               "coordinate and give the run real APE to convert). The load "
               "partition -- datum bt_H_ref = b - z_draft, "
               "p_top = rho_ref*g*z_draft in the pa(nz+1) BC -- plus the "
               "MOM6 trimmed IC (trim_ic_for_p_surf: the column top starts "
               "where the stratified water above it weighs the load) makes "
               "the barotropic state at rest, so what is measured is the "
               "FV-MOM6 truncation alone. Yung, Hallberg, Adcroft & Morrison "
               "(2026), JAMES 18, e2025MS005645 report order 1e-9 m/s for "
               "their CORRECTED algorithm (abstract; sigma icemount 1e-7 -> "
               "1e-12, their SS5.1.2) on damped configurations; this build "
               "implements none of their corrections and reads "
               "|u|_rms = 6.1e-05 m/s at their 10-day horizon with no "
               "dissipation at all.",
          tags=["cavity", "ice_shelf", "sloping_lid", "calving_front",
                "vcoord_sigma", "pgf_fv_mom6", "zinit_linear", "rest"]),

    _case("isomip_plus_ice_free_zfixed",
          V + "isomip_plus/ocean0_ice_free_zfixed.nml",
          "rest", 576, 144, en_rest_max=1.0e-05, t1_timeout=1200,
          t2_timeout=300, t2_skip=False,
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 300.0,
                            "nu_h": 6.0, "nghost": 2,
                            "has_western_boundary": True, "eddying": False},
          note="THE PARTIAL-STEP GATE, and the cheapest honest statement of "
               "the z-level staircase defect: the ISOMIP+ bed and trough on "
               "vcoord_type='z_fixed' with NO ice shelf, NO cavity code and "
               "NO melt, at rest, with &vcoord_nml zfixed_closed_faces on. "
               "A face whose layer is an inert `zstar_h_min` filler on either "
               "side is a z-LEVEL WALL for that layer (Adcroft, Hill & "
               "Marshall 1997; Losch 2008 2.1) -- no normal velocity, no mass "
               "or tracer flux, free-slip -- and leaving it OPEN lets the FV "
               "pressure gradient integrate across a staircase step of up to "
               "h_nominal = 20 m. MEASURED (gfortran 15.1 Release, "
               "pred_corr, 2 days): knob ON En 1.006E-08 at day 2 and "
               "SATURATED (1.000E-08 at day 1.75), MaxCFL <= 0.0008, Mass "
               "Error 4.2E-14; knob OFF 5.629E-04 at day 2 and still "
               "climbing, MaxCFL to 0.086 -- 56000x. With the face mask but "
               "WITHOUT the open-column barotropic weighting the same leg "
               "sat at 2.700E-06 and had not settled, which is what the "
               "barotropic-consistency slice bought. Note the throwaway "
               "spike reported NaN by day 0.25 "
               "for its own cavity-off variant of ocean0_idealised_draft; "
               "THIS namelist does not go non-finite in 2 days, so the bar "
               "here is the ENERGY, not survival. nu_h = 6 and kappa_h = 1 "
               "are the protocol values and are deliberately NOT zeroed: "
               "they are what exercises the free-slip closure in the "
               "harmonic velocity-Laplacian and the per-layer mask on the "
               "lateral tracer flux. split_scheme is not pinned, so the "
               "ssp_rk2 twin is built from the same file. conserve:Salt/Heat "
               "were scoped XFAIL here (a one-time -1.524E-06 / -2.030E-06 "
               "step at the first remap: the IC relamp draining the "
               "fillers' tracer content) until the vanished-layer content "
               "rule I1' landed; MEASURED with it (gfortran 15.1 Release, "
               "tier 2, 0.5 d): worst |Error| Salt 1.663E-13, Heat "
               "3.396E-13, Mass 1.804E-14 (ssp_rk2 twin 2.108E-13 / "
               "3.836E-13 / 2.043E-14), peak En 1.025E-08.",
          tags=["isomip_plus", "vcoord_z_fixed", "closed_faces",
                "partial_steps", "rest", "pgf_fv_mom6", "zinit_linear"]),
    _case("isomip_plus_ocean0_idealised",
          V + "isomip_plus/ocean0_idealised_draft.nml",
          "forced", 1000, 100, t1_timeout=1200, t2_timeout=300,
          # THE GROUNDED-CAVITY CONSERVATION GATE, at case scale.
          #
          # 3778 of 9600 interior columns (39.4 %) GROUND -- the ice draft
          # meets the bed and they become land through `seed_wet_mask_impl`
          # -- which is the largest grounded fraction anywhere in this
          # corpus and the configuration that exposed the land-state defect
          # (a grounded column's seeded `h_layer` is NEGATIVE, so the old
          # land tracer hold was an algebraic identity and left a
          # full-column `hTr` beside a floored `h`).  Measured at day 1
          # pre-fix: Salt Error 6.088E-01, Heat -7.994E-02 -- a step change
          # at step 1, flat thereafter -- against a Mass that closed at
          # -5.0E-14.  Post-fix all three read ~5E-14.
          #
          # So this row is carried for `conserve:*` (BUDGET_EXACT, 1e-11)
          # above everything else: `en_max` and `cfl:no-runaway` are the
          # ordinary finite/stability guards, and the case is not a
          # rest-state gate at all (melt and the northern sponge drive it).
          # A budget that steps is the defect, and only a case with real
          # grounded ice can see it.
          #
          # Cost: 0.63 s/step on one core. Tier 1 was the protocol's own 1
          # simulated day (288 steps at dt = 300) until 2026-09-20, and that
          # length is exactly why the day-3.2 blow-up below shipped unseen:
          # the case is CLEAN at every one of those 288 steps. Tier 1 is now
          # 1000 steps (day 3.47), which is past the failure; the run aborts
          # at step 910, so it costs ~575 s, not 630. Tier 2 stays 100 steps
          # -- it is the cheap smoke test and lengthening it buys nothing the
          # tier-1 row does not already say. Not downscaled -- 240x40x36 is
          # the COM resolution the grounding line is resolved at, and halving
          # it moves the grounding line rather than the numerics.
          known_failure={
              "tiers": [1],
              # Scoped to `completed` ALONE, and measured: at 1000 steps the
              # run aborts at step 910 and EVERY other assertion still
              # passes -- conserve:{Mass,Salt,Heat} at round-off, energy and
              # CFL inside their bars over the 3 days before the abort. So
              # this marker excuses the abort and nothing else, and a budget
              # or energy regression still turns the row red.
              "assertions": ["completed"],
              "reason":
                  "OPEN, and the sigma pressure-gradient truncation over the "
                  "RESOLVED ISOMIP+ BED -- the same term as "
                  "cavity_sloping_lid_rest, 3958x larger. The ISOMIP+ "
                  "channel-wall term (their Eq. 4, d_c = 500 m, f_c = 4 km) "
                  "drops the bed 122 m across one 2 km cell, so at the trough "
                  "sidewall a 23.1 m water column sits beside a 145.5 m one: "
                  "stiffness rx0 = |dH|/(Ha+Hb) = 0.726, 3.6x the classical "
                  "Beckmann-Haidvogel (1993) sigma bound of 0.2, with 168 "
                  "wet-wet faces over it. Cubed, a_peak = N^2 De^3/(6 dx Hbar) "
                  "= 1.48E-05 m/s^2 against the sloping lid's 3.73E-09, i.e. "
                  "U = a/|f| = 10.5 cm/s of spurious geostrophic flow on one "
                  "row of cells. En e-folds on ~1 day, hits &ocean_pgf_nml "
                  "maxvel at step 907, drives bt_eta to -36.8 m in a <= 23 m "
                  "column at step 909, and mints a NEGATIVE h_layer "
                  "(-0.710 m) at step 910; the melt driver's non-finite guard "
                  "is the detector, not the site. Substituted and INERT: melt "
                  "(off merely moves it to day 3.6), the freshwater form, "
                  "both implicit drag folds, the sponge, the calving front, "
                  "convective adjustment, Gamma_T, n_inner (4x), and the "
                  "outer split. Uniform density holds En = 8.5E-22 -- machine "
                  "zero -- for 7 days on the identical geometry; nz = "
                  "12/18/36/72 all fail within 0.2 day of each other (the "
                  "truncation is nz-independent, as derived); a FLAT lid over "
                  "the same bed still fails, so the tilted boundary that "
                  "matters is the bed, not the ice base. Every palliative only "
                  "postpones (h_min_cavity 25/30/50 m -> day 7.6/11.2/19.8; "
                  "nu_h 20/60 -> day 14.2/25.2), and the only two changes that "
                  "reach 30 days -- nu_h = 600 (100x Table 4) and dt = 75 s -- "
                  "do it by letting the mode SATURATE at En 1.4E-04 / 9.4E-04, "
                  "which is a quiet wrong answer rather than a fix. The "
                  "namelist is therefore NOT tuned to pass this row. "
                  "conserve:{Mass,Salt,Heat} are deliberately NOT excused: "
                  "they are what this case is carried for, they are evaluated "
                  "over the 3 days before the abort, and they still hold at "
                  "round-off. The fix is Phase 6 (the sloping-coordinate PGF "
                  "corrections), not a tolerance and not a shorter run.",
              "ref": "validation_examples/ocean/isomip_plus/README.md",
          },
          t2_dimensionless={"dx": 2000.0, "dy": 2000.0, "dt": 300.0,
                            "nu_h": 6.0, "nghost": 2,
                            "has_western_boundary": False, "eddying": False},
          note="ISOMIP+ Ocean0 with the idealised linear draft -- the "
               "grounded-cavity conservation gate.  See "
               "validation_examples/ocean/isomip_plus/README.md and "
               "tests/test_ocean_cavity_grounded_budget.F90, which pins the "
               "same statement at unit scale and to 1e-12.",
          tags=["cavity", "ice_shelf", "grounded", "basal_melt", "sponge",
                "vcoord_sigma", "pgf_fv_mom6", "zinit_linear", "conserve"]),

    # ============ ISOMIP+ Ocean0 melt OFF under z_fixed: gate E6 ============
    # The v0.1.0 stability gate for the z_fixed x cavity envelope
    # (python_prototypes design/cavity_rest_growth_diagnosis.md section Q).
    # The shipped `ocean0_idealised_zfixed.nml` with basal melt and top drag
    # OFF, carried to 180 days -- the length is load-bearing: the numerical
    # regime (a one-row calving-front partial-cell jet) does not appear until
    # day ~105 and saturates by ~165, so a 30- or 90-day leg sees only the
    # physical regime and gates nothing about it.  Two regimes, both bounded:
    #   1. d0-~100, PHYSICAL: the kappa_v = 5E-05 diffusive boundary anomaly
    #      against the insulating SLOPED boundaries (trough sidewalls, ice
    #      base) drives an along-slope geostrophic current, En ~ t^2 -> t
    #      (Phillips 1970 / Wunsch 1970).  Not a mode; protocol physics.
    #   2. d~105-165, NUMERICAL: one row of 7.3 m partial top cells at the
    #      calving front grows a one-cell jet (Re_dx ~ 7) on a 15-18-day
    #      e-folding and SATURATES at ~3-4E-06 m2/s2 (~2 cm/s); nu_h >= 30
    #      removes it outright.
    # Each gate below is section Q.6 fix 1's, and each fails on nu_h = 0 and
    # passes on the protocol.  The day-anchored claims read the TWICE-DAILY
    # [stats] series (360 samples over 51840 steps), so they are exact days.
    # Tier 1 only (V100, ~18 min per leg): 240x40x36 at 2 km is the COM
    # resolution the grounding line and the calving-front cut are resolved
    # at, and 180 days is ~50x a tier-2 budget.
    _case("isomip_plus_ocean0_zfixed_meltoff",
          V + "isomip_plus/ocean0_idealised_zfixed.nml",
          "rest", 51840, 0, en_rest_max=1.0e-05, budget_tol=BUDGET_OPEN,
          t1_timeout=3000, t1_samples=360,
          t1_overrides={"ocean_cavity_melt_nml": {"enable": False},
                        "ocean_tdrag_nml": {"enable": False},
                        # the shipped melt diagnostics refuse to register
                        # with melt off (a plane of missing values).
                        "ocean_diag_nml": {"diags": '""'}},
          t2_skip=True,
          t2_reason="the numerical regime needs ~105 days to appear and "
                    "~150 to saturate, and 240x40x36 @ 2 km is the "
                    "resolution the calving-front partial cells are cut at; "
                    "coarsening moves the cut, shortening removes the "
                    "regime. Tier 1 only.",
          claims=[
              {"kind": "en_at_day", "name": "en30-regime1", "day": 30.0,
               "max": 1.0e-07,
               "text": "section Q.6 gate (a): the 30-day energy is the "
                       "diffusive boundary current alone (5.632E-08 measured)",
               "meaning": "Above 1E-07 at day 30 the inviscid mode (nu_h = 0: "
                          "1.44E-06) or a new amplifier is back.",
               "ref": "validation_examples/ocean/isomip_plus/"
                      "ocean0_idealised_zfixed.nml header"},
              {"kind": "en_rate_decel", "name": "decelerates-by-d30",
               "early": [15.0, 20.0], "late": [25.0, 30.0],
               "text": "section Q.6 gate (a): regime 1 is a power law, so its "
                       "5-day log-rate FALLS (0.066 -> 0.051 /day measured)",
               "meaning": "An accelerating rate at day 30 is an exponential "
                          "mode (nu_h = 0: 0.142 -> 0.220 /day), not the "
                          "t^2 -> t boundary current."},
              {"kind": "en_at_day", "name": "en180-bounded", "day": 180.0,
               "max": 1.0e-05,
               "text": "section Q.6 gate (c): the calving-front jet saturates "
                       "near 3.5E-06 (3.147E-06 measured at day 180)"},
              {"kind": "en_day_ratio", "name": "saturated-d150-180",
               "num_day": 180.0, "den_day": 150.0, "max": 2.0,
               "text": "section Q.6 gate (c): regime 2 SATURATES over "
                       "d150-180 (1.39 measured; growth stops by ~d165)",
               "meaning": "A ratio >= 2 over the last 30 days is a 17-day "
                          "e-folding still running: the jet did not "
                          "saturate."},
          ],
          # No known_failure: energy:rest-settles PASSES here (final En
          # 3.147E-06 = 83% of the 3.804E-06 peak at d173.5, V100
          # 2026-09-24, bebt = 0.1 defaults). Regime 2 OSCILLATES about its
          # saturation level, so this gate reads the phase of that
          # oscillation at day 180 -- it missed by a hair (95.3%) on the
          # pre-bebt tree. If it trips again, localise before re-marking.
          note="ISOMIP+ Ocean0 idealised, z_fixed + zfixed_closed_faces, melt "
               "and top drag OFF, 180 days: gate E6 of v0.1.0 (bounded and "
               "explained). See the namelist header and "
               "docs/CAPABILITIES_AND_LIMITATIONS.md.",
          tags=["isomip_plus", "cavity", "ice_shelf", "vcoord_z_fixed",
                "closed_faces", "partial_steps", "sponge", "pgf_fv_mom6",
                "zinit_linear", "rest", "long_run"]),
    _case("isomip_plus_ocean0_zfixed_meltoff_nu30",
          V + "isomip_plus/ocean0_idealised_zfixed.nml",
          "rest", 51840, 0, en_rest_max=REST_1MM_S, budget_tol=BUDGET_OPEN,
          t1_timeout=3000, t1_samples=360,
          t1_overrides={"ocean_cavity_melt_nml": {"enable": False},
                        "ocean_tdrag_nml": {"enable": False},
                        "ocean_hvisc_nml": {"nu_h": 30.0},
                        "ocean_diag_nml": {"diags": '""'}},
          t2_skip=True,
          t2_reason="twin of isomip_plus_ocean0_zfixed_meltoff; the same "
                    "180-day horizon and COM grid. Tier 1 only.",
          # THE BAR is section Q.6 gate (c)'s twin: En(180 d) < 5E-07, which
          # is REST_1MM_S exactly. At nu_h = 30 the calving-front jet never
          # exists (viscous decay ~nu/dx^2 = 0.65 /day beats its ~0.2 /day
          # generation) and the run stays on the regime-1 power law
          # (2.12E-07 measured at day 180, En ~ t^0.9).
          known_failure={
              "tiers": [1],
              "assertions": ["energy:rest-settles"],
              "reason":
                  "BY DESIGN: with the calving-front jet gone the run sits "
                  "on the regime-1 power law to the end (En ~ t^0.9, "
                  "1.160E-07 d90 -> 2.121E-07 d180, measured on a V100 "
                  "2026-09-24), so its final sample IS its peak. That is "
                  "the kappa_v boundary current (Phillips 1970 / Wunsch "
                  "1970), protocol physics, not a mode. energy:rest at "
                  "5E-07 is the gate this row exists for.",
              "ref": "docs/CAPABILITIES_AND_LIMITATIONS.md",
          },
          note="the nu_h = 30 twin of isomip_plus_ocean0_zfixed_meltoff: the "
               "calving-front partial-cell jet is ABSENT, so peak En stays "
               "under 1 mm/s rms for 180 days.",
          tags=["isomip_plus", "cavity", "ice_shelf", "vcoord_z_fixed",
                "closed_faces", "partial_steps", "sponge", "pgf_fv_mom6",
                "zinit_linear", "rest", "long_run"]),

    # ===================== geometry / masking ==============================
    _case("island_at_rest", V + "island_at_rest/island_at_rest.nml",
          "rest", 2880, 300, en_rest_max=REST_1UM_S,
          note="interior land masking with a motionless start: the free-slip "
               "metric-zeroed walls must not leak a wall-normal velocity.",
          tags=["land_mask", "rest_exact"]),
    _case("coriolis_coast", V + "coriolis_coast/coriolis_coast.nml",
          "forced", 2880, 300, tags=["coriolis", "wall_bc"],
          ),
    _case("flow_past_island", V + "flow_past_island/flow_past_island.nml",
          "forced", 1440, 300, tags=["island", "wake"]),
    _case("double_drake", V + "double_drake/double_drake.nml",
          "forced", 1440, 250, tags=["double_drake", "spherical"]),
    _case("geostrophic_adjustment",
          V + "geostrophic_adjustment/geostrophic_adjustment.nml",
          "adiabatic", 4320, 400, en_growth_max=1.5,
          note="NK=1 barotropic adjustment from an initial SSH bump: the "
               "front radiates gravity waves and settles into geostrophic "
               "balance. Nothing adds energy, so En must not grow.",
          tags=["barotropic", "nk1", "geostrophic"]),

    # ===================== tracers / sponges / tides =======================
    _case("ideal_age_demo", V + "ideal_age/ideal_age_demo.nml",
          "forced", 864, 250, tags=["ideal_age", "passive_tracer"]),
    _case("sponge_real_demo", V + "sponge_demo/sponge_real_demo.nml",
          "forced", 864, 250, budget_tol=BUDGET_OPEN,
          tags=["sponge", "nudging"]),
    _case("body_tide_basin", V + "tides/body_tide_basin.nml",
          "forced", 2880, 400, tags=["tides", "body_forcing", "sal"]),
    _case("neverworld2", V + "neverworld2/neverworld2.nml",
          "forced", 192, 150, en_max=EN_CEIL_HI, budget_tol=BUDGET_OPEN,
          tags=["neverworld2", "idealised_global"]),

    # ===================== eddying gyres ===================================
    _case("eddy_test", V + "eddy_test/eddy_test.nml", "forced", 400, 0,
          t1_timeout=900, t2_skip=True,
          t2_reason="eddy-resolving closed basin at 5 km. A closed basin "
                    "cannot be shortened (the basin IS the case), so the only "
                    "downscale is to coarsen to 10 km, which halves Rd/dx and "
                    "violates rule R4. Tier 1 only.",
          tags=["eddying", "gyre"]),
    _case("eddy_test_quick", V + "eddy_test/eddy_test_quick.nml", "forced", 100, 0,
          t1_timeout=900, t2_skip=True,
          t2_reason="1024x1024x15 closed basin -- 15.7M cells. Too large for a "
                    "hosted runner even at 20 steps, and coarsening breaks "
                    "rule R4 as for eddy_test. Tier 1 only.",
          tags=["eddying", "gyre", "large"]),

    # ===================== PGF-form coverage ===============================
    # `mont` and `fv_wright` were selected by ZERO shipped namelists, which is
    # exactly how each came to carry an undetected sloping-bathymetry defect
    # (`mont` NaNs on every sloping case; `fv_wright` NaNs at step ~24 on the
    # stratified zstar_sigma seamount). These two cases exercise them on a FLAT
    # BED, where both are valid, so a regression is caught.
    #
    # Cross-form invariant: on a flat bed mont / fv_wright / fv_lite must agree.
    # Measured on this config, all three give En = 1.220E-05.
    _case("pgf_mont", V + "pgf_forms/pgf_mont.nml", "forced", 720, 60,
          t1_timeout=1200,
          tags=["pgf", "pgf_form_coverage", "flat_bed", "mont"]),

    _case("pgf_fv_wright", V + "pgf_forms/pgf_fv_wright.nml", "forced", 720, 60,
          t1_timeout=1200,
          tags=["pgf", "pgf_form_coverage", "flat_bed", "fv_wright"]),

    # ===================== ALE benchmark ===================================
    _case("benchmark_ale", V + "benchmark_ale/benchmark_ale.nml",
          "forced", 240, 50, budget_tol=BUDGET_OPEN,
          t2_overrides={"grid_nml": {"nx": 60, "ny": 30,
                                     "dx": 1.66666, "dy": 1.66666}},
          t2_dimensionless={"dx": 1.66666 * 111320.0, "dy": 1.66666 * 111320.0,
                            "dt": 1800.0, "nu_h": 0.0, "lat_deg": -45.0,
                            "nghost": 3, "has_western_boundary": False,
                            "eddying": False},
          tags=["ale_remap", "zstar_sigma"]),

    # ===================== sea ice =========================================
    _case("polar_freezeup_thermo", V + "sea_ice/polar_freezeup_thermo.nml",
          "rest", 2000, 300, en_rest_max=REST_1UM_S, budget_tol=BUDGET_OPEN,
          note="no wind: the ocean stays motionless while the ice thermo "
               "runs. Any current is spurious.",
          tags=["sea_ice", "thermo"]),
    _case("polar_freezeup_dynamics", V + "sea_ice/polar_freezeup_dynamics.nml",
          "forced", 72, 72, budget_tol=BUDGET_OPEN,
          tags=["sea_ice", "evp", "dynamics"]),
    _case("sea_ice_pack", V + "sea_ice_pack/sea_ice_pack.nml",
          "forced", 288, 288, budget_tol=BUDGET_OPEN,
          tags=["sea_ice", "itd", "transport"]),

    # ===================== performance benchmarks (tier 1 only) ============
    # These live in bench_scaling/ and exist to measure throughput at size.
    # They still must not NaN, so they are gated at tier 1; a downscaled twin
    # of a scaling benchmark measures nothing.
    _case("bench_benchmark_ale_big", V + "bench_scaling/benchmark_ale_big.nml",
          "forced", 192, 0, budget_tol=BUDGET_OPEN, t1_timeout=900,
          t2_skip=True, t2_reason=_BENCH_REASON, tags=["benchmark"]),
    _case("bench_double_gyre_big", V + "bench_scaling/double_gyre_big.nml",
          "forced", 150, 0, t1_timeout=900, t2_skip=True,
          t2_reason=_BENCH_REASON, tags=["benchmark"]),
    _case("bench_dg50_bt0", V + "bench_scaling/dg50_bt0.nml", "forced", 144, 0,
          t1_timeout=900, t2_skip=True, t2_reason=_BENCH_REASON,
          tags=["benchmark", "bt0"]),
    _case("bench_dg50_bt8", V + "bench_scaling/dg50_bt8.nml", "forced", 144, 0,
          t1_timeout=900, t2_skip=True, t2_reason=_BENCH_REASON,
          tags=["benchmark", "bt_halo"]),
    _case("bench_dg50_bt8_fv", V + "bench_scaling/dg50_bt8_fv.nml", "forced", 144, 0,
          t1_timeout=900, t2_skip=True, t2_reason=_BENCH_REASON,
          tags=["benchmark", "bt_halo", "pgf_fv"]),
    _case("bench_seamount_bt0", V + "bench_scaling/seamount_bt0.nml", "rest", 300, 0,
          en_rest_max=REST_SEAMOUNT, t1_timeout=900, t2_skip=True,
          t2_reason=_BENCH_REASON, tags=["benchmark", "bt0"]),
    _case("bench_seamount_bt8", V + "bench_scaling/seamount_bt8.nml", "rest", 300, 0,
          en_rest_max=REST_SEAMOUNT, t1_timeout=900, t2_skip=True,
          t2_reason=_BENCH_REASON, tags=["benchmark", "bt_halo"]),
    _case("bench_seamount_extreme_bt8", V + "bench_scaling/seamount_extreme_bt8.nml",
          "rest", 150, 0, en_rest_max=REST_SEAMOUNT, t1_timeout=900, t2_skip=True,
          t2_reason=_BENCH_REASON, tags=["benchmark", "extreme"]),
]


# ---------------------------------------------------------------------------
# The OUTER SPLIT-SCHEME axis
# ---------------------------------------------------------------------------
# `&ocean_bt_nml split_scheme` selects the outer split-explicit integrator.
# Every case above runs the DEFAULT, `"pred_corr"` -- the MOM6
# predictor-corrector. `"ssp_rk2"`, the two-stage SSP average, is the
# EXPERIMENTAL alternative: still shipped, still supported, and the only
# scheme wired through `eulerian_z`, dynamic wet/dry and
# `dt_tracer_advect_ratio > 1`, which is why six namelists pin it. It is
# labelled experimental for ONE measured reason -- it manufactures internal
# gravity waves out of a stratified REST state, four decades above the
# default (see `resting_stratified_channel`) -- and it must stay under test
# rather than rot into a cold branch of the dispatcher.
#
# This axis is how. Each entry is re-run with `split_scheme` forced to
# `ssp_rk2` under the name `<case>__ssp_rk2`, against the SAME physics
# assertions -- a scheme change must not change what a case is ALLOWED to do,
# so anywhere the two schemes disagree, the disagreement is the finding.
#
# THE AXIS INVERTED ON 2026-09-14, when `pred_corr` became the default. Before
# that the base corpus ran `ssp_rk2` and the twins were `<case>__mom6_pc`;
# now the base corpus runs `pred_corr` and the twins are `<case>__ssp_rk2`.
# The case COUNT is unchanged -- the same namelists are excluded, for the
# same reasons -- but every twin row is renamed, and the XFAIL that records
# the resting-state defect moved with the scheme it belongs to.
#
# WHY A FULL AXIS AND NOT A CURATED SUBSET. Measured on this tree, 4 CPU
# workers, gfortran: the whole tier-2 corpus is 2 min 43 s. Doubling it costs
# ~2.5 minutes of CI, which is cheaper than the reasoning required to defend
# any particular subset -- and a subset argued from "these are the cases where
# the schemes could differ" is exactly the argument that would have missed
# BOTH of the cases that actually separate them: a quiescent stratified
# channel and a 48^2 idealised coastline. Tier 1 is a different budget
# (full-scale, minutes per case, GPU), so there the axis IS curated, by
# `SCHEME_AXIS_TIER1`.
#
# EXCLUDED, and why -- a case is skipped iff it cannot legally run both:
#   * the six `dt_tracer_advect_ratio = 2` namelists pin `ssp_rk2` in the
#     file, because the pred_corr v1 envelope refuses that configuration
#     fail-loud (the `nz = 1` namelist used to pin for the same reason and
#     no longer does -- the single-layer vertical-friction tridiagonal was
#     repaired and the refusal lifted, 2026-09-14);
#   * the four `*_pred_corr` namelists pin `pred_corr` -- they were written
#     when that scheme was the non-default one and their whole subject is
#     it; they now pin the default explicitly, which keeps them out of the
#     axis and keeps the twin count where it was;
#   * the tier-1-only scaling benchmarks measure throughput at size; a second
#     scheme does not change what they measure.
# The skip is DERIVED from the namelist, not listed by hand, so a file that
# gains or loses a pin cannot silently fall out of the axis.

# ---------------------------------------------------------------------------
# The VERTICAL-COORDINATE REST MATRIX
# ---------------------------------------------------------------------------
# Every case above runs its namelist under whatever vertical coordinate that
# namelist chose. None of them runs the SAME problem under every coordinate,
# which is the question "which coordinate is trustworthy on which geometry?"
# and the question that hid five separate defects on the cavity branch (see
# `vcoord_matrix.__doc__`).
#
# `vcoord_matrix.build_matrix` emits one namelist per (geometry x family x
# stratification x EOS) cell into `tmp_local_artifacts/vcoord_matrix/` at
# import time and returns the rows. The namelists are generated rather than
# committed: dozens of near-identical files rot, and a reader cannot tell
# which cell of the matrix a given file is. They are left on disk, so any
# cell is reproducible by hand.
import sys as _sys
import os as _os

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import vcoord_matrix as _vcm  # noqa: E402


class _Bars(object):
    """The bar constants, handed to the generator so it cannot drift from
    the values the rest of this file uses."""
    REST_SEAMOUNT = REST_SEAMOUNT
    REST_1MM_S = REST_1MM_S
    REST_100UM_S = REST_100UM_S
    REST_1UM_S = REST_1UM_S


STABILITY_CASES += _vcm.build_matrix(_case, _Bars)


SCHEME_AXIS_SCHEME = "ssp_rk2"

# Cases whose tier-1 physics gates were calibrated under ONE scheme and whose
# measured value differs under the other. Each entry REPLACES keys in the
# twin's physics block; it never widens the base case's.
SCHEME_AXIS_PHYSICS = {}

# Cases that FAIL under the second scheme for a known, documented reason.
# Reported as XFAIL so CI red stays informative -- and so that when someone
# fixes the underlying defect the case turns XPASS and says so. `assertions`
# scopes the marker so one documented defect does not excuse a case from every
# other gate it has.
SCHEME_AXIS_KNOWN_FAILURE = {

    # `coriolis_coast` lived here as the BLOCKER entry until 2026-09-14.
    # It is gone because the defect is FIXED, not because the gate was
    # relaxed: the fast-loop Coriolis reference `subtract_fast_cor_ref`
    # removes was being evaluated on the stage-entry `u^n` while the slow
    # Coriolis inside `F_bt` was evaluated on `u_av`, so under pred_corr the
    # uncancelled `f x (v_av_bar - v^n_bar)` forced every barotropic substep
    # and pumped the basin's gravest Poincare seiche. `set_cor_ref_velocity`
    # now builds the reference from the same velocity the slow tendency used
    # (MOM6 `ubt_Cor`). Measured on a V100, 20 simulated days, same binary,
    # same file: before, En 6.98e-07 (d8) -> 1.23e-05 (d10) -> 2.78e-04 (d12)
    # -> 1.57e-01 (d16) -> NaN (d18); after, En stays in 2.0e-07..1.1e-06 for
    # the full 20 days and exits 0 -- the same band ssp_rk2 holds. The
    # permanent guard is `tests/test_ocean_cor_ref_seiche.F90`, an unforced
    # rotating closed basin that asserts KE+PE does not grow.
    "resting_stratified_channel": {
        "assertions": ["energy:rest", "energy:rest-settles"],
        "reason":
            "THE reason ssp_rk2 is labelled EXPERIMENTAL, and this twin is "
            "its permanent regression test. The SSP two-stage average "
            "amplifies a gravity wave by sqrt(1 + (omega dt)^4 / 4) per "
            "step, so a motionless stratified channel seeded with +-0.5 mK "
            "of noise develops 7.7 mm/s of current in 25 days out of NO "
            "energy source: En 2.992E-05 (nvfortran) / 2.72e-05 (gfortran) "
            "against the 0.5e-08 bar -- 6000x over -- still climbing on a "
            "2.5-day e-folding when the run ends. The base case, on the "
            "default pred_corr, holds 1.739E-09 on the identical file: "
            "17000x less. The cause is the OUTER split and nothing else -- "
            "substituting the Coriolis form, the ALE remap and the PGF form "
            "each moved the answer < 0.1%, so all three are EXONERATED; "
            "removing the lateral viscosity RAISED it. Do NOT close this by "
            "raising en_rest_max or by putting nu_h back into the file "
            "(nu_h = 100 hides it at En 1.6e-11 without touching its rate). "
            "ssp_rk2 stays supported and stays on this axis.",
        "ref": "validation_examples/ocean/eady/"
               "resting_stratified_channel.nml",
    },
    "cavity_sloping_lid_rest": {
        "assertions": ["energy:rest", "energy:rest-settles"],
        "reason":
            "The same ssp_rk2 defect as `resting_stratified_channel`, seen "
            "on an internal-wave field this case GENERATES for itself rather "
            "than one a noise seed puts there: the sloping-lid sigma "
            "truncation keeps forcing the cavity, and the SSP two-stage "
            "average amplifies the resulting waves by "
            "sqrt(1 + (omega dt)^4 / 4) per step. The two schemes are "
            "indistinguishable for 10 days (En 1.875E-09 vs 1.877E-09, "
            "trimmed IC), then ssp_rk2 leaves: 6.06E-08 (d13), 5.76E-06 "
            "(d15), 4.24E-05 (d20), 1.850E-04 (d30) -- 370x over the 0.5e-06 "
            "bar and 8000x the default pred_corr's 2.298E-08 on the "
            "identical file. The base case keeps GATING energy:rest for "
            "exactly that reason and is excused only on energy:rest-settles "
            "(tier 1). Do NOT close this by "
            "raising en_rest_max or by adding viscosity to the namelist "
            "(nu_h = 6 and nu_4 = 1e7 each delay the growth ~10 days and "
            "change nothing about its rate). ssp_rk2 stays supported and "
            "stays on this axis.",
        "ref": "validation_examples/ocean/ice_shelf_cavity/README.md",
    },
}

# Tier 1 is minutes per case on a GPU, so its axis is curated rather than
# complete. The selection covers, deliberately:
#   * every vcoord family that reaches the outer loop -- zstar_full
#     (double_gyre_mom6), zstar (eady, double_gyre_weno3), zstar_sigma
#     (benchmark_ale), sigma (seamount), lagrangian (bc_inst_tuned_512);
#   * both PGF-sensitive cases (pgf_mont, pgf_fv_wright), which is where a
#     change in WHEN the pressure gradient is evaluated would show;
#   * the baroclinic cases, where the schemes must agree on a PHYSICAL
#     growth rate (eady, baroclinic_2layer, bc_inst_tuned_512);
#   * the quiescent cases, where they must agree on zero
#     (resting_stratified_channel, seamount, island_at_rest);
#   * one surface-flux case and one open-boundary case (kpp_basin,
#     seamount_obc_baroclinic), because the salt/heat budget ledger is
#     weighted PER SCHEME (`ocean_budget_stage_weight`) and a wrong weight
#     there is invisible in every case with no source and no boundary;
#   * the LAND-MASKED cases (coriolis_coast, island_at_rest,
#     flow_past_island). This is where pred_corr's land-mask NaN was found;
#     it is fixed (see SCHEME_AXIS_KNOWN_FAILURE), and coriolis_coast stays
#     on the axis as the regression guard for it.
SCHEME_AXIS_TIER1 = [
    "resting_stratified_channel", "eady", "double_gyre_mom6",
    "double_gyre_weno3", "baroclinic_2layer", "bc_inst_tuned_512",
    "seamount", "island_at_rest", "benchmark_ale", "pgf_mont",
    "pgf_fv_wright", "kpp_basin", "seamount_obc_baroclinic",
    "coriolis_coast", "flow_past_island",
]


def _nml_pins_scheme(path):
    """True iff the COMMITTED namelist sets `split_scheme` itself.

    Read from the file rather than listed here: a namelist that gains or
    loses the pin must fall in or out of the axis on its own.
    """
    import os
    full = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, path)
    try:
        with open(full) as fh:
            for line in fh:
                bare = line.split("!", 1)[0].strip().lower()
                if bare.startswith("split_scheme"):
                    return True
    except OSError:
        return False
    return False


def _scheme_twin(case, scheme):
    """Build the `<name>__<scheme>` twin of one base case."""
    import copy
    twin = copy.deepcopy(case)
    twin["name"] = "{}__{}".format(case["name"], scheme)
    twin["split_scheme"] = scheme
    twin["tags"] = list(case.get("tags", [])) + ["split_scheme", scheme]
    for tier in ("tier1", "tier2"):
        spec = twin.get(tier)
        if not spec or spec.get("skip"):
            continue
        spec.setdefault("overrides", {}).setdefault(
            "ocean_bt_nml", {})["split_scheme"] = '"{}"'.format(scheme)
    # tier 1 is curated; drop the twin's tier-1 spec when not selected.
    if case["name"] not in SCHEME_AXIS_TIER1:
        twin["tier1"] = {"skip": True,
                         "reason": "tier-1 scheme axis is curated for wall "
                                   "time; see SCHEME_AXIS_TIER1"}
    over = SCHEME_AXIS_PHYSICS.get(case["name"])
    if over:
        twin["physics"] = dict(twin.get("physics") or {})
        twin["physics"].update(over)
    kf = SCHEME_AXIS_KNOWN_FAILURE.get(case["name"])
    if kf:
        twin["known_failure"] = dict(kf)
    elif case.get("matrix"):
        # A vcoord-matrix row's twin follows its own rule:
        #   * a REFUSAL row asserts that `validate_config` rejects the
        #     configuration, which is a property of the CONFIGURATION and not
        #     of the outer time-split, so the twin inherits the refusal;
        #   * a RUNNABLE cell's twin is MEASURED in its own right --
        #     `vcoord_matrix_pin.py` pins the `__ssp_rk2` runs separately
        #     (`MEASURED_TWIN`) -- so it carries its own record, or none when
        #     it was measured passing. The two schemes genuinely differ on
        #     this matrix (ssp_rk2 amplifies the inviscid mode, and is the
        #     scheme that RESTS the viscous leg's stepped cells), so
        #     inheriting the base marker would either excuse a twin-only
        #     failure or report a permanent XPASS.
        tk = case["matrix"].get("twin_known_failure", "inherit")
        if case["matrix"].get("expect") == "run" and tk != "inherit":
            if tk:
                twin["known_failure"] = dict(tk)
            else:
                twin.pop("known_failure", None)
        elif case.get("known_failure"):
            twin["known_failure"] = dict(case["known_failure"])
        else:
            twin.pop("known_failure", None)
    else:
        twin.pop("known_failure", None)
    return twin


def _build_scheme_axis():
    out = []
    for case in STABILITY_CASES:
        if _nml_pins_scheme(case["nml"]):
            continue
        t1_skip = case.get("tier1", {}).get("skip")
        t2_skip = case.get("tier2", {}).get("skip")
        if t2_skip and case["name"] not in SCHEME_AXIS_TIER1:
            # tier-1-only scaling benchmark, not selected for the tier-1
            # axis: the twin would have nothing to run.
            continue
        if t1_skip and t2_skip:
            continue
        out.append(_scheme_twin(case, SCHEME_AXIS_SCHEME))
    return out


# ---------------------------------------------------------------------------
# Remap precondition guard -- ON in every case whose coordinate REMAPS
# ---------------------------------------------------------------------------
# `&vcoord_nml remap_check_preconditions` asserts, once per ALE remap, the two
# things the overlap sweep has always ASSUMED: non-negative source and target
# thicknesses, and matching column totals. Outside them the sweep silently
# creates or deletes tracer mass -- no NaN, no bounds hit, no budget entry --
# so a stability run is exactly where the check earns its cost (one device
# reduction per thermo step). It is switched on here, by OVERRIDE, never by
# editing the shipped namelists: it is a diagnostic, not physics, and default
# off is what users get.
#
# Which families: `ocean_apply_ale_remap_step` returns before touching
# anything for `eulerian_z` and `lagrangian` (the former holds h at H*dsig by
# vertical-advection cancellation, the latter's target IS the current h), so
# the pair the check judges (`remap_h_old`, `target_h`) is never written
# there and the check would assert on stale scratch. Every other family --
# sigma (the namelist default), zsigma, zstar, zstar_sigma, zstar_full,
# z_fixed, rho, hycom -- relayers every thermo step and gets the guard.
#
# DERIVED from the namelist (plus any tier override of `vcoord_type`), not
# listed by hand, so a file that changes coordinate falls in or out on its
# own.
#
# `&vcoord_nml check_vanished_content` (the vanished-layer content assertion)
# is NOT switched on: it lands with `fix/remap-vanished-layer-content`, which
# is not on this branch's base. Add it alongside this one when it merges.
REMAP_CHECK_KNOB = "remap_check_preconditions"
_NO_REMAP_VCOORDS = ("lagrangian", "isopycnal", "eulerian_z", "z")
_REMAP_VCOORDS = ("sigma", "zsigma", "z-sigma", "z_sigma", "zstar", "z-star",
                  "z_star", "zstar_lite", "zstar_full", "z-star-full",
                  "z_star_full", "zstarfull", "zstar_sigma", "z-star-sigma",
                  "z_star_sigma", "zstarsigma", "z_fixed", "z_levels",
                  "gprime", "rho", "isopycnic", "rho_target", "hycom",
                  "hybrid")


def _nml_vcoord_type(path):
    """`&vcoord_nml vcoord_type` of the COMMITTED namelist, lower-cased.

    Absent => "sigma", the `config_t` default (`rdb_config`).
    """
    import os
    full = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, path)
    group = None
    with open(full) as fh:
        for line in fh:
            bare = line.split("!", 1)[0].strip()
            if bare.startswith("&"):
                group = bare[1:].split()[0].lower() if len(bare) > 1 else None
                continue
            if bare.startswith("/"):
                group = None
                continue
            if group == "vcoord_nml" and "=" in bare:
                key, val = bare.split("=", 1)
                if key.strip().lower() == "vcoord_type":
                    return val.strip().rstrip(",").rstrip("/").strip() \
                              .strip("'\"").strip().lower()
    return "sigma"


def _case_remaps(case, spec):
    over = (spec.get("overrides") or {}).get("vcoord_nml", {})
    vc = over.get("vcoord_type")
    if vc is None:
        vc = _nml_vcoord_type(case["nml"])
    vc = str(vc).strip("'\"").strip().lower()
    if vc in _NO_REMAP_VCOORDS:
        return False
    if vc in _REMAP_VCOORDS:
        return True
    raise ValueError("{}: unknown vcoord_type {!r} -- teach "
                     "_REMAP_VCOORDS / _NO_REMAP_VCOORDS about it".format(
                         case["name"], vc))


def _add_remap_checks():
    for case in STABILITY_CASES:
        for tier in ("tier1", "tier2"):
            spec = case.get(tier)
            if not spec or spec.get("skip"):
                continue
            if not _case_remaps(case, spec):
                continue
            spec.setdefault("overrides", {}).setdefault(
                "vcoord_nml", {})[REMAP_CHECK_KNOB] = True


_add_remap_checks()
STABILITY_CASES += _build_scheme_axis()


def _self_check():
    names = [c["name"] for c in STABILITY_CASES]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise ValueError("duplicate stability case names: {}".format(sorted(dup)))
    base = {c["name"] for c in STABILITY_CASES if not c.get("split_scheme")}
    for n in SCHEME_AXIS_TIER1:
        if n not in base:
            raise ValueError("SCHEME_AXIS_TIER1 names a case that does not "
                             "exist: {}".format(n))
    for tbl, label in ((SCHEME_AXIS_PHYSICS, "SCHEME_AXIS_PHYSICS"),
                       (SCHEME_AXIS_KNOWN_FAILURE, "SCHEME_AXIS_KNOWN_FAILURE")):
        for n in tbl:
            if n not in base:
                raise ValueError("{} names a case that does not exist: "
                                 "{}".format(label, n))


_self_check()
