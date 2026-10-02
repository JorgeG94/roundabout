"""What the compatibility matrix EXPECTS to fail, and why.

`compat_matrix.py` runs a pairwise covering array of configurations.  Every
refusal and every runtime failure it meets must be explained by a row here,
or the cell FAILS.  Two classes of row:

* **PHYSICAL** -- an exclusion that is the design, forever (KPP xor EPBL, the
  Henyey latitude factor on a grid with no latitude, ...).
* **KNOWN_GAP** -- a combination that SHOULD work and does not yet, with the
  reason, an owner, and the v0.1.0 tracker item
  (`python_prototypes/design/v010_blockers.md`) that will close it.

A KNOWN_GAP row is an XFAIL: when the gap is fixed the cell passes, the row
turns XPASS, and an XPASS FAILS the suite until the row is deleted.  The list
can only shrink.

Rows match on FEATURES, which are read off the cell's merged namelist (with
the model's own defaults), never off the axis value names -- so a row still
matches when a value is renamed or a second axis value turns the same knob
on.  A `refused` row must also name the refusal: its `message` regex has to
match the line the model logged, so a cell refused for a DIFFERENT reason is
still a failure.  A `runtime` row lists the outcomes it expects (CRASH,
NONFINITE, BUDGET).

Adding a row: run the matrix, read the FAIL, decide whether it is physics or a
gap, and give it the narrowest `when` that explains it.  Never widen a row to
swallow a failure you have not read.
"""

import re

V010 = "python_prototypes/design/v010_blockers.md"


def _g(nml, grp, key, default=None):
    return nml.get(grp, {}).get(key, default)


def _vtype(nml):
    return str(_g(nml, "vcoord_nml", "vcoord_type", "sigma"))


def _edges(nml):
    return [str(_g(nml, "ocean_bc_nml", e, "wall")) for e in ("west", "east", "south", "north")]


# name -> (meaning, predicate over the merged namelist).  Model defaults are
# spelled out where the knob may be absent (e.g. KPP is ON unless disabled).
FEATURES = {
    # vertical coordinate
    "vc_sigma": ("sigma", lambda n: _vtype(n) == "sigma"),
    "vc_zstar": ("z*-lite", lambda n: _vtype(n) == "zstar"),
    "vc_zstar_full": ("per-column z*", lambda n: _vtype(n) == "zstar_full"),
    "vc_zstar_sigma": ("sigma shallow / z* deep", lambda n: _vtype(n) == "zstar_sigma"),
    "vc_z_fixed": ("fixed z levels", lambda n: _vtype(n) == "z_fixed"),
    "closed_faces": ("z_fixed partial-step closed faces",
                     lambda n: _vtype(n) == "z_fixed" and bool(_g(n, "vcoord_nml", "zfixed_closed_faces", False))),
    "open_steps": ("z_fixed WITHOUT closed faces",
                   lambda n: _vtype(n) == "z_fixed" and not _g(n, "vcoord_nml", "zfixed_closed_faces", False)),
    "vc_hycom": ("hybrid z*/isopycnal", lambda n: _vtype(n) == "hycom"),
    "vc_rho": ("isopycnal targets", lambda n: _vtype(n) == "rho"),
    "vc_eulerian_z": ("eulerian z (H*dsig)", lambda n: _vtype(n) == "eulerian_z"),
    "vc_lagrangian": ("pure Lagrangian", lambda n: _vtype(n) == "lagrangian"),
    "vc_zsigma": ("smoothstep sigma->z", lambda n: _vtype(n) == "zsigma"),
    # outer split
    "pred_corr": ("predictor-corrector split",
                  lambda n: _g(n, "ocean_bt_nml", "split_scheme", "pred_corr") == "pred_corr"),
    "ssp_rk2": ("SSP-RK2 split", lambda n: _g(n, "ocean_bt_nml", "split_scheme", "pred_corr") == "ssp_rk2"),
    # vertical mixing
    "kpp": ("KPP boundary layer", lambda n: bool(_g(n, "ocean_vmix_nml", "use_closure", True))
            and bool(_g(n, "ocean_vmix_nml", "use_kpp", True))),
    "epbl": ("energetic PBL", lambda n: bool(_g(n, "ocean_epbl_nml", "enable", False))),
    "kappa_shear": ("JHL08 shear mixing", lambda n: bool(_g(n, "ocean_kappa_shear_nml", "enable", False))),
    "kappa_shear_vertex": ("kappa-shear at corners",
                           lambda n: bool(_g(n, "ocean_kappa_shear_nml", "enable", False))
                           and bool(_g(n, "ocean_kappa_shear_nml", "at_vertex", False))),
    "tidal_mixing": ("St-Laurent tidal mixing", lambda n: bool(_g(n, "ocean_tidal_mixing_nml", "enable", False))),
    "conv": ("convective adjustment", lambda n: bool(_g(n, "ocean_conv_nml", "enable", False))),
    "ddiff": ("double diffusion", lambda n: bool(_g(n, "ocean_ddiff_nml", "enable", False))),
    "bryan_lewis": ("Bryan-Lewis background", lambda n: bool(_g(n, "ocean_vmix_nml", "bkgnd_profile", False))),
    "henyey": ("Henyey background", lambda n: bool(_g(n, "ocean_vmix_nml", "bkgnd_henyey", False))),
    # lateral momentum
    "smag": ("Smagorinsky harmonic", lambda n: _g(n, "ocean_hvisc_nml", "lateral_closure", "none") == "smagorinsky"),
    "smag_ah": ("flow-aware biharmonic (Smagorinsky AH)", lambda n: bool(_g(n, "ocean_hvisc_nml", "smag_ah", False))),
    "leith": ("Leith harmonic", lambda n: _g(n, "ocean_hvisc_nml", "lateral_closure", "none") == "leith"),
    "leith_biharm": ("Leith biharmonic", lambda n: _g(n, "ocean_hvisc_nml", "lateral_closure", "none") == "leith_biharm"),
    "nu_4": ("constant biharmonic", lambda n: float(_g(n, "ocean_hvisc_nml", "nu_4", 0.0)) > 0.0),
    "stress_tensor": ("MOM6 stress-tensor operator", lambda n: bool(_g(n, "ocean_hvisc_nml", "stress_tensor", False))),
    "kh_aniso": ("anisotropic viscosity", lambda n: float(_g(n, "ocean_hvisc_nml", "kh_aniso", 0.0)) > 0.0),
    "meke_backscatter": ("MEKE backscatter", lambda n: bool(_g(n, "ocean_meke_nml", "backscatter", False))),
    # eddy parameterisation
    "slopes": ("isopycnal slopes", lambda n: bool(_g(n, "ocean_slopes_nml", "enable", False))),
    "gm": ("Gent-McWilliams", lambda n: bool(_g(n, "ocean_gm_nml", "enable", False))),
    "redi": ("Redi", lambda n: bool(_g(n, "ocean_redi_nml", "enable", False))),
    "meke": ("MEKE", lambda n: bool(_g(n, "ocean_meke_nml", "enable", False))),
    "varmix": ("VarMix", lambda n: bool(_g(n, "ocean_varmix_nml", "enable", False))),
    "mle": ("Fox-Kemper MLE", lambda n: bool(_g(n, "ocean_foxkemper_nml", "enable", False))),
    # tracers
    "ideal_age": ("ideal age", lambda n: bool(_g(n, "ocean_tracers_nml", "enable_ideal_age", False))),
    "pseudo_salt": ("pseudo-salt", lambda n: bool(_g(n, "ocean_tracers_nml", "enable_pseudo_salt", False))),
    # PGF / EOS
    "pgf_mont": ("Montgomery PGF", lambda n: _g(n, "ocean_pgf_nml", "form", "mont") == "mont"),
    "pgf_fv_mom6": ("FV-MOM6 PGF", lambda n: _g(n, "ocean_pgf_nml", "form", "mont") == "fv_mom6"),
    "pgf_recon": ("in-layer T/S reconstruction for the PGF",
                  lambda n: bool(_g(n, "ocean_pgf_nml", "reconstruct_for_pressure", False))),
    "eos_linear": ("linear EOS", lambda n: _g(n, "ocean_eos_nml", "eos", "linear") == "linear"),
    "eos_wright": ("Wright EOS", lambda n: _g(n, "ocean_eos_nml", "eos", "linear") == "wright"),
    "eos_roquet": ("Roquet SpV EOS", lambda n: _g(n, "ocean_eos_nml", "eos", "linear") == "roquet_spv"),
    # Coriolis
    "cor_sadourny": ("Sadourny enstrophy", lambda n: _g(n, "ocean_coriolis_nml", "form", "sadourny") == "sadourny"),
    "cor_energy": ("Sadourny energy", lambda n: _g(n, "ocean_coriolis_nml", "form", "sadourny") == "sadourny_energy"),
    "cor_hk": ("Hollingsworth-Kallen", lambda n: _g(n, "ocean_coriolis_nml", "form", "sadourny") == "sadourny_hk"),
    "pv_weno": ("any WENO PV interpolation",
                lambda n: str(_g(n, "ocean_coriolis_nml", "pv_adv_scheme", "centered")).startswith("weno")),
    # barotropic options
    "bt_bc_pgf": ("BT corrector baroclinic-PGF retro-correction",
                  lambda n: bool(_g(n, "ocean_bt_nml", "correction_bc_pgf", False))),
    "substep_drag": ("BT substep drag", lambda n: bool(_g(n, "ocean_bt_nml", "substep_drag", False))),
    "wave_drag": ("BT linear wave drag", lambda n: bool(_g(n, "ocean_bt_nml", "wave_drag", False))),
    "h_weighted": ("h-weighted BT corrector", lambda n: bool(_g(n, "ocean_bt_nml", "correction_h_weighted", False))),
    "visc_rem": ("visc_rem BT corrector", lambda n: bool(_g(n, "ocean_bt_nml", "correction_visc_rem", False))),
    "implicit_drag": ("implicit bottom-drag fold", lambda n: bool(_g(n, "ocean_vdiff_nml", "implicit_drag", False))),
    # geometry / grid / forcing
    "periodic_x": ("re-entrant in x", lambda n: _edges(n)[0] == "periodic"),
    "sponge": ("sponge edge", lambda n: "sponge" in _edges(n)),
    "obc_open": ("Flather open edge", lambda n: "open" in _edges(n)),
    "walls_only": ("closed basin", lambda n: all(e == "wall" for e in _edges(n))),
    "tripolar": ("tripolar fold", lambda n: _edges(n)[3] == "tripolar_fold"),
    "cavity": ("ice-shelf cavity", lambda n: bool(_g(n, "ocean_cavity_dyn_nml", "enable", False))),
    "cartesian": ("Cartesian grid", lambda n: _g(n, "ocean_grid_nml", "grid_config", "cartesian") == "cartesian"),
    "spherical": ("spherical sector", lambda n: _g(n, "ocean_grid_nml", "grid_config", "cartesian") == "spherical"),
    "sw_pen": ("penetrating shortwave", lambda n: float(_g(n, "ocean_thermo_nml", "sw_pen_frac", 0.0)) > 0.0),
    "cooling": ("surface cooling", lambda n: float(_g(n, "ocean_thermo_nml", "q_heat", 0.0)) < 0.0),
}


def features(nml):
    return {name: bool(fn(nml)) for name, (_, fn) in FEATURES.items()}


class Row(object):
    """One expected failure.  See the module docstring."""

    def __init__(self, rid, cls, kind, when, reason, message=None, expect=(),
                 unless=(), owner="-", link="-", scope="cell", witness=None):
        assert cls in ("PHYSICAL", "KNOWN_GAP"), cls
        assert kind in ("refused", "runtime"), kind
        assert scope in ("cell", "any"), scope
        assert kind == "runtime" or scope == "cell", rid + ": a refusal is deterministic"
        assert (scope == "any") == (witness is not None), \
            rid + ": a scope='any' row (and only one) pins the witness cell that fails"
        assert message, rid + ": a row must name the refusal / failure text it explains"
        assert kind != "runtime" or expect, rid + ": a runtime row must list its outcomes"
        assert cls != "KNOWN_GAP" or (owner != "-" and link != "-"), rid + ": a gap needs owner + link"
        self.rid, self.cls, self.kind = rid, cls, kind
        self.when, self.unless = tuple(when), tuple(unless)
        self.reason, self.owner, self.link = reason, owner, link
        self.message = re.compile(message) if message else None
        self.expect = frozenset(expect)
        # "cell": every matching cell must fail (a passing one is an XPASS).
        # "any":  the gap bites in SOME combinations only.  The row pins one
        #         `witness` cell ({axis: value} over compat_matrix.BASE_CELL)
        #         that is evaluated on every run and MUST fail with the
        #         row's signature -- when it passes, the row is an XPASS.
        #         Other matching cells may pass or fail (XFAIL) freely.
        self.scope = scope
        self.witness = dict(witness) if witness else None

    def matches(self, feats):
        return all(feats[f] for f in self.when) and not any(feats[f] for f in self.unless)

    def explains(self, msg):
        return bool(self.message and self.message.search(msg))


def _gap(rid, kind, when, reason, item, message=None, expect=(), unless=(), owner="orchestrator",
         scope="cell", witness=None):
    return Row(rid, "KNOWN_GAP", kind, when, reason, message=message, expect=expect,
               unless=unless, owner=owner, link="{}: {}".format(V010, item), scope=scope,
               witness=witness)


def _phys(rid, when, reason, message, unless=()):
    return Row(rid, "PHYSICAL", "refused", when, reason, message=message, unless=unless)


# The remap precondition guard (`remap_check_preconditions`, ON in every
# cell) stopping the very first step: the regrid handed it a bad column.
_PRECOND_STEP1 = r"remap preconditions at step 1;"

ROWS = [
    # ===================================================================
    # PHYSICAL -- design exclusions, expected forever.
    # ===================================================================
    _phys("pv_weno_needs_sadourny", ("pv_weno",),
          "The WENO PV face interpolation is MOM6's WENOVI family, which is itself an "
          "enstrophy-form Coriolis scheme; it has no energy-form or HK counterpart.",
          r"pv_adv_scheme='weno\d' is only wired into form='sadourny'", unless=("cor_sadourny",)),
    _phys("henyey_needs_latitude", ("henyey", "cartesian"),
          "The Henyey background is a latitude factor; a Cartesian grid has no latitude "
          "(geolatT = 0 would make every column equatorial).",
          r"bkgnd_henyey=\.true\. requires a non-cartesian"),

    # ===================================================================
    # KNOWN_GAP -- refusals.
    # ===================================================================
    _gap("zsigma_units", "refused", ("vc_zsigma",),
         "z_ref_global is filled dimensionless (k/nz) but the zsigma deep branch reads metres; "
         "the whole column would collapse into the bed layer.",
         "Other coordinates, 'zsigma refused (z_ref_global units defect)'",
         message=r"vcoord_type = 'zsigma' is refused on the ocean path"),
    _gap("pred_corr_eulerian_z", "refused", ("vc_eulerian_z", "pred_corr"),
         "pred_corr's v1 envelope: the legacy eulerian_z per-stage vertical-advection + "
         "h-rescale path is not wired into the predictor-corrector.",
         "NOT TRACKED (CLAUDE.md 'pred_corr v1 envelope')",
         message=r"split_scheme='pred_corr' requires an ALE vertical coordinate"),
    _gap("mle_needs_epbl", "refused", ("mle",),
         "Fox-Kemper MLE reads epbl%mld only; MOM6 takes the mixed-layer depth from any "
         "boundary-layer scheme (KPP included).",
         "NOT TRACKED (B5 reads epbl%mld)",
         message=r"ocean_foxkemper_nml: enable=\.true\. requires ocean_epbl_nml enable=\.true\.",
         unless=("epbl",)),
] + [
    _gap("closed_faces_" + f, "refused", ("closed_faces", f),
         "GM / Redi / MLE fold their face fluxes from a 2-D wet gate after the per-layer "
         "closed-face mask: the transports would leak through a closed face.",
         "item {} (feat/{}-zfixed-closed-faces)".format(item, f),
         message=r"zfixed_closed_faces does not yet compose with GM / Redi / MLE")
    for f, item in (("gm", 1), ("redi", 2), ("mle", 3))
] + [
    _gap("closed_faces_nu_4", "refused", ("closed_faces", "nu_4"),
         "The free-slip closure of a closed face exists for the harmonic velocity-Laplacian "
         "kernels only.",
         "item 4 (fix/biharmonic-zfixed-closed-faces lifts the nu_4 refusal)",
         message=r"zfixed_closed_faces does not yet compose with the BIHARMONIC viscosity"),
    _gap("closed_faces_stress_tensor", "refused", ("closed_faces", "stress_tensor"),
         "stress_tensor's tension/shear use 2-D wet masks (kh_aniso rides on it).",
         "item 10 (queued)",
         message=r"zfixed_closed_faces does not yet compose with the BIHARMONIC viscosity"),
    _gap("closed_faces_bc_pgf", "refused", ("closed_faces", "bt_bc_pgf"),
         "compute_pbce / compute_gtot_faces / the bc-PGF block weight by the FULL column.",
         "item 9 (queued)",
         message=r"zfixed_closed_faces does not yet compose with &ocean_bt_nml correction_bc_pgf"),
    _gap("closed_faces_substep_drag", "refused", ("closed_faces", "substep_drag"),
         "compute_bt_rem damps from the FULL-column face depth.",
         "item 5 (fix/bt-upstream-h-face-closed-faces ports it)",
         message=r"zfixed_closed_faces does not yet compose with &ocean_bt_nml substep_drag"),
    _gap("closed_faces_wave_drag", "refused", ("closed_faces", "wave_drag"),
         "compute_bt_rem_wave_drag damps from the FULL-column face depth.",
         "item 5 (fix/bt-upstream-h-face-closed-faces ports it)",
         message=r"zfixed_closed_faces does not yet compose with &ocean_bt_nml wave_drag"),

    # ----- under an ice-shelf cavity (single-rank row) -------------------
    _gap("cavity_vcoord", "refused", ("cavity",),
         "Under a cavity v1 accepts sigma / zstar (they rescale the live column) and z_fixed "
         "(taught the ice base) only; the others are unvalidated or draft-following.",
         "'Cavity-only refusals' (+ design/phase6_zlike_coordinates_under_ice.md)",
         message=r"ocean_cavity_dyn_nml enable=\.true\. accepts vcoord_type='sigma', 'zstar'",
         unless=("vc_sigma", "vc_zstar", "vc_z_fixed")),
] + [
    _gap("cavity_zfixed_" + f, "refused", ("cavity", "vc_z_fixed", f),
         "Under a cavity the z_fixed top layers inside the draft are fillers, and this closure "
         "still closes its surface row on k = nz (not the first live layer k_top).",
         "'Cavity-only refusals' ({})".format(note),
         message=r"{} is refused with vcoord_type='z_fixed' under a cavity".format(msg))
    for f, msg, note in (
        ("kpp", r"use_kpp=\.true\.", "follow-up P6.4"),
        ("epbl", r"&ocean_epbl_nml enable=\.true\.", "follow-up P6.4"),
        ("kappa_shear", r"&ocean_kappa_shear_nml enable=\.true\.", "P6.3/P6.4"),
        ("tidal_mixing", r"&ocean_tidal_mixing_nml enable=\.true\.", "P6.3/P6.4"),
        ("ideal_age", r"enable_ideal_age=\.true\.", "follow-up P6.3"),
        ("gm", r"&ocean_gm_nml enable=\.true\.", "coordinate study"),
        ("slopes", r"&ocean_redi_nml / &ocean_slopes_nml enable=\.true\.", "coordinate study"),
        ("pgf_recon", r"reconstruct_for_pressure=\.true\.", "partial top cell reads fillers"),
    )
] + [
    Row("cavity_needs_fv_mom6", "PHYSICAL", "refused", ("cavity", "pgf_mont"),
        "The ice load enters through the FV pressure-stack top boundary condition "
        "pa(nz+1); Montgomery hard-zeroes M(nz), so it has nowhere to put it.",
        message=r"(p_top_in_bc=\.true\.|ocean_cavity_dyn_nml enable=\.true\.) requires "
                r"(&ocean_pgf_nml )?form='fv_mom6'"),

    # ===================================================================
    # KNOWN_GAP -- runtime failures of accepted configurations.
    # ===================================================================
    _gap("land_column_regrid_rho_hycom", "runtime", ("vc_hycom",),
         "The rho/hycom regrid writes a negative thickness in a LAND column (10 x 1.5e-4 m "
         "minus 9 inflated layers = -1.2e-3 m); the remap precondition guard stops step 1.",
         "item C3 (fix/rho-regrid-land-negative-thickness)", expect=("CRASH",),
         message=_PRECOND_STEP1),
    _gap("land_column_regrid_rho", "runtime", ("vc_rho",),
         "As land_column_regrid_rho_hycom, on the pure isopycnal coordinate.",
         "item C3 (fix/rho-regrid-land-negative-thickness)", expect=("CRASH",),
         message=_PRECOND_STEP1),
    _gap("land_column_target_zstar_full", "runtime", ("vc_zstar_full",),
         "zstar_full builds a LAND column's target as nz x zstar_h_min (1.0e-3 m) against a "
         "column of nz x H_VANISHED (1.5e-3 m): a 1/3 column-total mismatch the remap "
         "precondition guard stops at step 1.  Same land-column class as C3; new site.",
         "item C3 (NEW site: the zstar_full target builder)", expect=("CRASH",),
         message=_PRECOND_STEP1),
    _gap("zfixed_open_steps", "runtime", ("open_steps",),
         "z_fixed WITHOUT closed faces takes the full staircase PGF at every step face: "
         "En ~30x the closed-face run in 24 steps and, in some combinations, a negative "
         "thickness the remap guard stops.  Decided: refuse it (bed steps) now.",
         "item 11 (fix/zfixed-require-closed-faces)", expect=("CRASH", "NONFINITE"),
         message=r"remap preconditions at step|nan-catch|I1' tripwire", scope="any",
         # measured 2026-10-02: a negative thickness at step 20
         witness={"vcoord": "z_fixed_open", "vmix_extra": "conv", "vmix_bg": "bryan_lewis",
                  "lateral": "stress_tensor", "eddy": "gm_varmix_resscaled",
                  "tracers": "pseudo_salt", "pgf": "fv_mom6_ppm", "eos": "linear",
                  "coriolis": "sadourny", "pv_adv": "weno7", "bt": "substep_drag",
                  "grid": "spherical"}),
    _gap("bc_pgf_needs_fv_mom6", "runtime", ("bt_bc_pgf", "pgf_mont"),
         "correction_bc_pgf reads pgf%e_face, which only the fv_mom6 PGF fills -- but the "
         "combination is ACCEPTED at configure and `error stop`s inside step 1 "
         "(compute_pbce).  The fix is a validate_config refusal.",
         "NOT TRACKED (found by this matrix, 2026-10-02)", expect=("CRASH",),
         message=r"compute_pbce: requires ocean_pgf_form = 'fv_mom6'"),
    _gap("redi_obc_salt_budget", "runtime", ("redi", "obc_open"),
         "Redi with a Flather open edge: the model's own salt budget misses ~4e-5 of the "
         "salt content in 24 steps (every such cell; Redi on walls / a periodic channel "
         "closes to 1e-15).  Either the neutral flux through the open face is not in the "
         "budget's boundary term or it is a real leak.",
         "NOT TRACKED (found by this matrix, 2026-10-02)", expect=("BUDGET",),
         message=r"^Salt residual"),
]


# ---------------------------------------------------------------------------
# Who tests the classifier (run by `compat_matrix.py self-test`).  Each takes
# the compat_matrix module and returns (ok, what).
# ---------------------------------------------------------------------------
def _t_unexplained(cm):
    cell = dict(cm.BASE_CELL)
    nml = cm.merged_namelist(cell)
    val = {"status": "refused", "rc": 3, "stage": "engine_setup",
           "messages": ["some refusal nobody wrote a row for"]}
    cls, _, _ = cm.classify(cell, nml, val, None)
    return cls == "FAIL", "an unexplained refusal is a FAIL"


def _t_validate_crash(cm):
    cell = dict(cm.BASE_CELL)
    val = {"status": "crashed", "rc": -11, "stage": None, "messages": ["Segmentation fault"]}
    cls, _, _ = cm.classify(cell, cm.merged_namelist(cell), val, None)
    return cls == "FAIL", "a crashing --validate-only is a FAIL"


def _t_runtime_unexpected(cm):
    cell = dict(cm.BASE_CELL)
    val = {"status": "accepted", "rc": 0, "stage": None, "messages": []}
    run = {"outcome": "CRASH", "detail": "x"}
    cls, _, _ = cm.classify(cell, cm.merged_namelist(cell), val, run)
    return cls == "FAIL", "an unexplained crash of an accepted cell is a FAIL"


def _t_rows_contract(cm):
    """Synthetic rows: an explained refusal is REFUSED_GAP; the same row on
    an accepted cell is XPASS; a runtime row turns a crash into XFAIL and a
    pass into XPASS; a refusal with the WRONG message stays a FAIL."""
    global ROWS
    saved = ROWS
    try:
        ROWS = [_gap("t_ref", "refused", ("closed_faces",), "t", 0, message=r"^synthetic refusal"),
                _gap("t_run", "runtime", ("closed_faces",), "t", 0, expect=("CRASH",),
                     message=r"^synthetic crash")]
        cell = dict(cm.BASE_CELL)          # z_fixed + closed faces
        nml = cm.merged_namelist(cell)
        ref = {"status": "refused", "rc": 3, "stage": "engine_setup",
               "messages": ["synthetic refusal: closed faces"]}
        wrong = dict(ref, messages=["a different reason"])
        acc = {"status": "accepted", "rc": 0, "stage": None, "messages": []}
        got = [cm.classify(cell, nml, ref, None)[0],
               cm.classify(cell, nml, wrong, None)[0],
               cm.classify(cell, nml, acc, {"outcome": "PASS", "detail": ""})[0]]
        ROWS = [ROWS[1]]
        got += [cm.classify(cell, nml, acc, {"outcome": "CRASH", "detail": "synthetic crash"})[0],
                cm.classify(cell, nml, acc, {"outcome": "PASS", "detail": ""})[0],
                cm.classify(cell, nml, acc, {"outcome": "BUDGET", "detail": "synthetic crash"})[0],
                cm.classify(cell, nml, acc, {"outcome": "CRASH", "detail": "another crash"})[0]]
        want = ["REFUSED_GAP", "FAIL", "XPASS", "XFAIL", "XPASS", "FAIL", "FAIL"]
        return got == want, ("row contract (refused / wrong message / xpass / xfail / xpass / "
                             "wrong outcome / wrong signature) {}".format(got))
    finally:
        ROWS = saved


def _t_witnesses(cm):
    """Every scope='any' row pins a real cell that its own predicate matches,
    and the row-level XPASS fires when that witness passes."""
    bad = []
    for r in ROWS:
        if r.scope != "any":
            continue
        cell = cm.witness_cell(r)
        if any(v not in cm.VALUE_NAMES.get(a, ()) for a, v in r.witness.items()):
            bad.append(r.rid + ": unknown axis/value")
            continue
        if not r.matches(features(cm.merged_namelist(cell))):
            bad.append(r.rid + ": witness does not match the row")
            continue
        rec = {"axes": cell, "class": "PASS", "rows": []}
        if (r.rid, "PASS") not in cm.row_xpasses([rec]):
            bad.append(r.rid + ": a passing witness is not a row XPASS")
        rec = {"axes": cell, "class": "XFAIL", "rows": [r.rid]}
        if any(rid == r.rid for rid, _ in cm.row_xpasses([rec])):
            bad.append(r.rid + ": a failing witness is a row XPASS")
    return not bad, "scope='any' witnesses are real, matching, and XPASS when they pass {}".format(bad)


SELF_TESTS = [("unexplained", _t_unexplained), ("validate_crash", _t_validate_crash),
              ("runtime_unexpected", _t_runtime_unexpected), ("rows_contract", _t_rows_contract),
              ("witnesses", _t_witnesses)]
