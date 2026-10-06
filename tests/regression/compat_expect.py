"""What the compatibility matrix EXPECTS to fail, and why.

`compat_matrix.py` runs a pairwise covering array of configurations.  Every
refusal and every runtime failure it meets must be explained by a row here,
or the cell FAILS.  Two classes of row:

* **PHYSICAL** -- an exclusion that is the design, forever (KPP xor EPBL, the
  Henyey latitude factor on a grid with no latitude, ...).
* **KNOWN_GAP** -- a combination that SHOULD work and does not yet, with the
  reason, an owner, and the v0.1.0 tracker item
  (`python_prototypes/design/v010_blockers.md`) that will close it.

Three kinds of row: `refused` (a configure-time refusal, `rdb
--validate-only`), `runtime` (an accepted cell that fails a check: CRASH,
NONFINITE, BUDGET, ENERGY, DECOMP, RESTART, CROSS_BACKEND -- the outcome of
each check, see `compat_matrix.OUTCOME_CHECK`), and `multirank` (a cell
that runs on one rank but whose DECOMP leg the model refuses at configure: a
single-rank feature).  A runtime row is judged only when a check it expects
ran, so a DECOMP gap neither explains nor XPASSes on a run without the MPI
leg.

A KNOWN_GAP row is an XFAIL: when the gap is fixed the cell passes, the row
turns XPASS, and an XPASS FAILS the suite until the row is deleted.  The list
can only shrink.

Rows match on FEATURES, which are read off the cell's merged namelist (with
the model's own defaults), never off the axis value names -- so a row still
matches when a value is renamed or a second axis value turns the same knob
on.  A `refused` row must also name the refusal: its `message` regex has to
match the line the model logged, so a cell refused for a DIFFERENT reason is
still a failure.  A `runtime` row lists the outcomes it expects (any of
`OUTCOMES`).

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
    "mle_mld_filter": ("MLE running-mean MLD filter",
                       lambda n: bool(_g(n, "ocean_foxkemper_nml", "enable", False))
                       and float(_g(n, "ocean_foxkemper_nml", "mld_decay_time", 0.0)) > 0.0),
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
    "visc_rem": ("visc_rem chain (D1 follow-up: was correction_visc_rem, "
                 "retired 2026-10)", lambda n: bool(_g(n, "ocean_bt_nml", "visc_rem_chain", False))),
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


OUTCOMES = frozenset(("CRASH", "NONFINITE", "BUDGET", "ENERGY", "DECOMP", "RESTART",
                      "CROSS_BACKEND"))


class Row(object):
    """One expected failure.  See the module docstring."""

    def __init__(self, rid, cls, kind, when, reason, message=None, expect=(),
                 unless=(), owner="-", link="-", scope="cell", witness=None, backend=None):
        assert cls in ("PHYSICAL", "KNOWN_GAP"), cls
        assert kind in ("refused", "runtime", "multirank"), kind
        assert set(expect) <= OUTCOMES, "{}: unknown outcome(s) {}".format(rid, set(expect) - OUTCOMES)
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
        # A backend-specific gap ("gpu" matches every gpu-* --backend): the
        # row is invisible to a run on any other backend.
        self.backend = backend

    def on(self, backend):
        return self.backend is None or str(backend).startswith(self.backend)

    def matches(self, feats):
        return all(feats[f] for f in self.when) and not any(feats[f] for f in self.unless)

    def explains(self, msg):
        return bool(self.message and self.message.search(msg))


def _gap(rid, kind, when, reason, item, message=None, expect=(), unless=(), owner="orchestrator",
         scope="cell", witness=None, backend=None):
    return Row(rid, "KNOWN_GAP", kind, when, reason, message=message, expect=expect,
               unless=unless, owner=owner, link="{}: {}".format(V010, item), scope=scope,
               witness=witness, backend=backend)


def _phys(rid, when, reason, message, unless=()):
    return Row(rid, "PHYSICAL", "refused", when, reason, message=message, unless=unless)


def _single_rank(rid, when, reason, item, message, unless=()):
    """A feature the model runs on ONE rank only and refuses, at configure,
    on more (the DECOMP leg's refusal).  Should work decomposed: a gap."""
    return _gap(rid, "multirank", when, reason, item, message=message, unless=unless)


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
    _phys("bc_pgf_needs_fv_mom6", ("bt_bc_pgf", "pgf_mont"),
          "correction_bc_pgf's compute_pbce builds the per-layer pressure response from the "
          "FV_MOM6 interface-height stack (pgf%e_face), which no other PGF form fills; the "
          "surface-relative forms carry no free-surface term (fix/bc-pgf-needs-fv-mom6).",
          r"correction_bc_pgf=\.true\. requires &ocean_pgf_nml form='fv_mom6'"),
    _phys("zfixed_open_steps", ("open_steps",),
          "z_fixed WITHOUT closed faces over a stepped bed pairs a live layer with a bed "
          "filler at every step face, and the staircase pressure gradient drives flow from "
          "rest; refused at configure (fix/zfixed-require-closed-faces).",
          r"zfixed_closed_faces = \.false\. is refused with vcoord_type='z_fixed' over a "
          r"stepped bed"),
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
    _gap("closed_faces_stress_tensor", "refused", ("closed_faces", "stress_tensor"),
         "stress_tensor's tension/shear use 2-D wet masks (kh_aniso rides on it).",
         "item 10 (queued)",
         message=r"zfixed_closed_faces does not yet compose with &ocean_hvisc_nml stress_tensor"),
    _gap("closed_faces_bc_pgf", "refused", ("closed_faces", "bt_bc_pgf"),
         "compute_pbce / compute_gtot_faces / the bc-PGF block weight by the FULL column.",
         "item 9 (queued)",
         message=r"zfixed_closed_faces does not yet compose with &ocean_bt_nml correction_bc_pgf"),

    # ----- under an ice-shelf cavity (single-rank row) -------------------
    _gap("cavity_vcoord", "refused", ("cavity",),
         "Under a cavity v1 accepts sigma (it rescales the live column) and z_fixed (taught "
         "the ice base) only; the others are unvalidated or draft-following -- zstar too, "
         "since it became MOM6 z* (a geopotential z_fixed stack dilated per column).",
         "'Cavity-only refusals' (+ design/phase6_zlike_coordinates_under_ice.md)",
         message=r"ocean_cavity_dyn_nml enable=\.true\. accepts vcoord_type='sigma' or 'z_fixed'",
         unless=("vc_sigma", "vc_z_fixed")),
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
    _gap("land_column_target_zstar_full", "runtime", ("vc_zstar_full",),
         "zstar_full builds a LAND column's target as nz x zstar_h_min (1.0e-3 m) against a "
         "column of nz x H_VANISHED (1.5e-3 m): a 1/3 column-total mismatch the remap "
         "precondition guard stops at step 1.  Same land-column class as C3; new site.",
         "item C3 (NEW site: the zstar_full target builder)", expect=("CRASH",),
         message=_PRECOND_STEP1),
    _gap("rho_runtime_crash", "runtime", ("vc_rho",),
         "Pure isopycnal (rho) stops within 2-12 steps on the matrix domains: the remap "
         "precondition guard, or the console CFL panic.  Before the land-column fix (item C3) "
         "the same cells stopped at step 1; with it they run further and die later "
         "(CLAUDE.md: rho is validation-grade alone, weakly stratified columns collapse).  "
         "Under the visc_rem chain (bt=visc_rem) the barotropic friction keeps some of these "
         "cells alive the full 24 steps, at ~190x the channel PASS energy instead of stopping "
         "(c042: chain off stops on the remap guard; chain on reaches En 1.5; the same cell "
         "under ssp_rk2 runs at 1.7e-2): the same rho collapse, a later symptom.",
         "NOT TRACKED (found by this matrix, 2026-10-05)", expect=("CRASH", "ENERGY"),
         message=r"remap preconditions at step \d+|console stats: CFL > panic threshold"
                 r"|x the \w+ PASS-population reference",
         scope="any",
         # c045 of the 2026-10-05 train run (as generated, not minimised)
         witness={"vcoord": "rho", "split": "pred_corr", "vmix_bl": "epbl",
                  "vmix_extra": "ddiff", "vmix_bg": "bryan_lewis", "lateral": "nu_4",
                  "eddy": "mle", "tracers": "pseudo_salt", "pgf": "fv_mom6_plm", "eos": "linear",
                  "coriolis": "sadourny", "pv_adv": "weno7", "bt": "correction_bc_pgf",
                  "geometry": "closed", "grid": "spherical", "forcing": "cool"}),
    _gap("zstar_open_steps_stress_tensor", "runtime", ("vc_zstar", "stress_tensor"),
         "zstar's open stepped bed (closed faces off: a live layer faces a 1e-4 m filler) "
         "with the MOM6 stress-tensor viscosity drives a layer negative and stops on the remap "
         "guard at step 2-3.  The corner shear stress is weighted by the ARITHMETIC 4-cell "
         "mean h_q (`hvisc_compute_stress`, Phase 2) while the divergence divides by the face "
         "thickness, so on a filler face beside a live corner the explicit viscous step is "
         "amplified by h_q/h_u ~ 1e5.  MOM6 forms hq as the harmonic-type mean of the four "
         "face thicknesses (MOM_hor_visc.F90 `hq = 2*h2uq*h2vq/(...)`), small whenever one "
         "face is vanished; substituting it runs this witness clean (En 6.46e-3 against the "
         "sigma twin's 6.73e-3).  The same operator defect is the vcoord matrix's FINDING A "
         "(thin density-space layers driven negative); the port changes every stress_tensor "
         "answer, so it is its own PR.",
         "NOT TRACKED (found by this matrix, 2026-10-04; vcoord matrix FINDING A)",
         expect=("CRASH",), message=r"remap preconditions at step \d+", scope="any",
         # minimised 2026-10-04 from c011 (greedy; every other axis at base)
         witness={"vcoord": "zstar", "lateral": "stress_tensor"}),
    _gap("zstar_open_steps_stress_tensor_gm", "runtime", ("vc_zstar", "stress_tensor", "gm"),
         "The zstar_open_steps_stress_tensor defect (arithmetic corner h_q against the face "
         "divisor on a filler face) with GM on: the cell runs 24 steps with GM off and dies "
         "with GM on.  Folded into continuity, GM's bolus drove a layer negative at step 2-3 "
         "(remap guard); as its own operator on the current h it lasts to step 11-19 and then "
         "trips the remap guard (sequential only) or the CFL panic (with MOM6 bottom-blocking). "
         "Lifts with the h_q port.",
         "NOT TRACKED (found by this matrix, 2026-10-05; GM sequential operator)",
         expect=("CRASH",),
         message=r"console stats: CFL > panic threshold|remap preconditions at step \d+",
         scope="any",
         # c010 of the 2026-10-05 train run (as generated, not minimised)
         witness={"vcoord": "zstar", "split": "ssp_rk2", "vmix_bl": "pp81",
                  "vmix_extra": "tidal", "vmix_bg": "scalar", "lateral": "kh_aniso",
                  "eddy": "gm_meke", "tracers": "ideal_age", "pgf": "mont", "eos": "linear",
                  "coriolis": "sadourny_hk", "pv_adv": "centered", "bt": "visc_rem",
                  "geometry": "closed", "grid": "cartesian", "forcing": "warm_sw"}),
    _gap("rho_gm_energy", "runtime", ("vc_rho", "gm"),
         "Pure isopycnal (rho) with GM: cells that stopped within a few steps under "
         "rho_runtime_crash (GM folded into continuity) now run all 24 steps with GM as its own "
         "operator on the current h, but at 13-200x the PASS-population energy -- the same "
         "validation-grade collapse of weakly stratified rho columns, reached later.",
         "NOT TRACKED (found by this matrix, 2026-10-05; GM sequential operator)",
         expect=("ENERGY",), message=r"x the \w+ PASS-population reference", scope="any",
         # c044 of the 2026-10-05 train run (as generated, not minimised)
         witness={"vcoord": "rho", "split": "pred_corr", "vmix_bl": "epbl",
                  "vmix_extra": "conv", "vmix_bg": "bryan_lewis", "lateral": "leith_biharm",
                  "eddy": "gm_meke", "tracers": "pseudo_salt", "pgf": "fv_mom6_plm",
                  "eos": "linear", "coriolis": "sadourny", "pv_adv": "weno5",
                  "bt": "substep_drag", "geometry": "obc", "grid": "spherical",
                  "forcing": "cool"}),
    _gap("hycom_runtime_crash", "runtime", ("vc_hycom",),
         "hycom with an open boundary stops at step 4 (remap precondition guard + nan-catch) "
         "now that the land-column crash (item C3) no longer stops it at step 1.",
         "NOT TRACKED (found by this matrix, 2026-10-05)", expect=("CRASH",),
         message=r"remap preconditions at step \d+", scope="any",
         # c002 of the 2026-10-05 train run (as generated, not minimised)
         witness={"vcoord": "hycom", "split": "pred_corr", "vmix_bl": "epbl",
                  "vmix_extra": "none", "vmix_bg": "henyey", "lateral": "kh_aniso",
                  "eddy": "mle", "tracers": "ts", "pgf": "fv_mom6", "eos": "roquet",
                  "coriolis": "sadourny_hk", "pv_adv": "centered", "bt": "substep_drag",
                  "geometry": "obc", "grid": "spherical", "forcing": "warm_sw"}),

    # ===================================================================
    # KNOWN_GAP -- the legs (phase 3, 2026-10-04): RESTART, DECOMP.
    # ===================================================================
    # `restart_visc_rem` CLOSED by PR-2 (bt-rem-from-av-rem, 2026-10-05):
    # the root cause (`visc_rem_precompute` reading the PREVIOUS stage's
    # `vmix%kv`, which the restart registry did not carry) is fixed by
    # registering `vmix_kv` (`ocean_state_build_restart_registry`,
    # `src/core/ocean/state/rdb_ocean_state.F90`) -- `test_engine_bit_exact_
    # visc_rem` (`tests/test_ocean_restart_engine.F90`) now runs the FULL
    # round trip (not resume-point-only) and is bit-exact.  Row deleted so
    # it does not XPASS.
    _gap("restart_mle_mld_filter", "runtime", ("mle_mld_filter",),
         "Fox-Kemper MLE with mld_decay_time > 0 does not resume bit-exact: its running-mean "
         "`mld_filtered` is persistent state but not in the restart registry, so a resume "
         "re-seeds it from the instantaneous MLD and the restratification flux changes.",
         "NOT TRACKED (found by this matrix, 2026-10-04)", expect=("RESTART",),
         message=r"differ after a warm restart"),
    _gap("decomp_eulerian_z_ssp_rk2", "runtime", ("vc_eulerian_z", "ssp_rk2", "epbl", "mle"),
         "eulerian_z under ssp_rk2 (its legacy per-stage vertical-advection + h-rescale path) "
         "with EPBL + Fox-Kemper MLE is not decomposition-invariant: last-bit differences in "
         "every owned cell on 2x2 and 4x1 (max|diff|/max|field| 3e-11 after 24 steps); "
         "EPBL alone, and the same pair on z_fixed / sigma or under pred_corr, are bitwise.  "
         "(Its visc_rem half was the visc_rem-weighted BT fold leaving non-image ghost "
         "velocities; fixed by the post-fold face refresh in `run_stage_split`.)",
         "NOT TRACKED (found by this matrix, 2026-10-04)", expect=("DECOMP",),
         message=r"\dx\d: \d+ field mismatch", scope="any",
         # re-minimised 2026-10-05 (greedy, every other axis at base)
         witness={"vcoord": "eulerian_z", "split": "ssp_rk2", "vmix_bl": "epbl", "eddy": "mle"}),

    # ===================================================================
    # KNOWN_GAP -- the GPU leg (nvfortran cc70) only.
    # ===================================================================
    _gap("gpu_eulerian_z_epbl_mle_drift", "runtime", ("vc_eulerian_z", "ssp_rk2", "epbl", "mle"),
         "eulerian_z + ssp_rk2 + EPBL + MLE ends 2e-7 (field norm) away from gfortran in 24 "
         "steps, 3000x the PASS population's spread: the same combination is not "
         "decomposition-invariant either (decomp_eulerian_z_ssp_rk2), so an order-sensitive "
         "operation amplified through the EPBL / MLE thresholds.",
         "NOT TRACKED (found by this matrix, 2026-10-04)", expect=("CROSS_BACKEND",),
         message=r"norm .* > band", scope="any", backend="gpu",
         # c049 of the 2026-10-04 array, as generated (not minimised: GPU-only)
         witness={"vcoord": "eulerian_z", "split": "ssp_rk2", "vmix_bl": "epbl",
                  "vmix_extra": "kappa_shear", "vmix_bg": "bryan_lewis",
                  "lateral": "meke_backscatter", "eddy": "mle", "tracers": "pseudo_salt",
                  "pgf": "mont", "eos": "linear", "coriolis": "sadourny", "pv_adv": "weno5",
                  "bt": "wave_drag", "geometry": "obc", "grid": "spherical"}),

    # ===================================================================
    # KNOWN_GAP -- single-rank features: the DECOMP leg's configure-time
    # refusal on more than one rank.
    # ===================================================================
    _single_rank("cavity_single_rank", ("cavity",),
                 "The ice-shelf cavity is single-rank in v1: the grounding statistics are "
                 "global reductions its configure does not take.",
                 "'Single-rank-only features to lift: ... cavity'",
                 message=r"&ocean_cavity_dyn_nml enable=\.true\. is single-rank in v1"),
]


# ---------------------------------------------------------------------------
# The per-PR slice: axis value -> the source paths that implement it.
#
# `compat_matrix.py run --diff-base REF` (or `--changed-files`) runs only the
# cells holding a value whose paths the diff touches.  A changed Fortran
# source NO value claims is shared code (the dynamical core, the continuity,
# the engine, the configuration) and runs the full matrix -- so a path
# missing here costs time, never coverage.  Globs are fnmatch patterns over
# repo-relative paths.  A value shared by all of an axis (the remap driver
# under every coordinate) touches every cell, which is the honest answer.
# ---------------------------------------------------------------------------
_P = "src/parameterizations/"
_VCOORD_ALL = ["src/core/ocean/vcoord/*", "src/ALE/*",
               "src/shared_module_utilities/rdb_vanished_layer.inc"]
_SPLIT = ["src/core/ocean/dynamics/split_rk2/*", "src/core/ocean/kernels/barotropic/*",
          "src/core/rdb_barotropic_workstate.F90"]
_VMIX = [_P + "vertical/rdb_ocean_vmix.F90", _P + "vertical/rdb_ocean_vdiff.F90"]
_HVISC = [_P + "lateral/rdb_ocean_horizontal_viscosity.F90"]
_PGF = ["src/pressure_force/rdb_ocean_pressure_force.F90"]
_EOS = ["src/equation_of_state/*", "src/core/ocean/state/rdb_ocean_eos_compute.F90"]
_COR = ["src/core/ocean/kernels/coriolis_adv/*",
        "src/shared_module_utilities/rdb_rel_vort_corner.inc"]
_SLOPES = [_P + "lateral/rdb_ocean_isopycnal_slopes.F90", _P + "lateral/rdb_ocean_lateral_mix.F90"]
_GMP = [_P + "lateral/rdb_ocean_gm.F90"] + _SLOPES
_MEKEP = [_P + "lateral/rdb_ocean_meke.F90"]
_VARMIXP = [_P + "lateral/rdb_ocean_varmix.F90", _P + "vertical/rdb_ocean_wave_speed.F90"]
_REDIP = [_P + "lateral/rdb_ocean_redi.F90"] + _SLOPES

VALUE_PATHS = {}
for _v in ("sigma", "zstar", "zstar_full", "zstar_sigma", "z_fixed_cf", "z_fixed_open", "hycom",
           "rho", "eulerian_z", "lagrangian", "zsigma"):
    VALUE_PATHS[("vcoord", _v)] = list(_VCOORD_ALL)
VALUE_PATHS[("vcoord", "z_fixed_cf")] += ["src/core/ocean/state/rdb_ocean_metrics.F90"]
for _v in ("pred_corr", "ssp_rk2"):
    VALUE_PATHS[("split", _v)] = list(_SPLIT)
VALUE_PATHS.update({
    ("vmix_bl", "kpp"): list(_VMIX),
    ("vmix_bl", "pp81"): list(_VMIX),
    ("vmix_bl", "epbl"): _VMIX + [_P + "vertical/rdb_ocean_epbl.F90"],
    ("vmix_extra", "kappa_shear"): [_P + "vertical/rdb_ocean_kappa_shear.F90"],
    ("vmix_extra", "kappa_shear_vertex"): [_P + "vertical/rdb_ocean_kappa_shear.F90"],
    ("vmix_extra", "tidal"): [_P + "vertical/rdb_ocean_tidal_mixing.F90"],
    ("vmix_extra", "conv"): list(_VMIX),
    ("vmix_extra", "ddiff"): list(_VMIX),
    ("vmix_bg", "scalar"): list(_VMIX),
    ("vmix_bg", "bryan_lewis"): list(_VMIX),
    ("vmix_bg", "henyey"): list(_VMIX),
    ("lateral", "meke_backscatter"): _HVISC + _MEKEP + _GMP,
    ("eddy", "gm"): list(_GMP),
    ("eddy", "gm_meke"): _GMP + _MEKEP + _VARMIXP,
    ("eddy", "redi"): list(_REDIP),
    ("eddy", "gm_redi_meke"): _GMP + _REDIP + _MEKEP + _VARMIXP,
    ("eddy", "gm_varmix_resscaled"): _GMP + _VARMIXP,
    ("eddy", "mle"): [_P + "lateral/rdb_ocean_mle.F90", _P + "vertical/rdb_ocean_epbl.F90"],
    ("tracers", "ts"): ["src/core/rdb_tracer.F90"],
    ("tracers", "ideal_age"): ["src/core/rdb_tracer.F90", "src/tracer/rdb_ocean_ideal_age.F90"],
    ("tracers", "pseudo_salt"): ["src/core/rdb_tracer.F90", "src/tracer/rdb_ocean_pseudo_salt.F90"],
    ("pgf", "mont"): list(_PGF),
    ("pgf", "fv_mom6"): list(_PGF),
    ("pgf", "fv_mom6_plm"): _PGF + ["src/pressure_force/rdb_ocean_pgf_reconstruct.F90"],
    ("pgf", "fv_mom6_ppm"): _PGF + ["src/pressure_force/rdb_ocean_pgf_reconstruct.F90"],
    ("eos", "wright"): list(_EOS),
    ("eos", "roquet"): _EOS + ["src/shared_module_utilities/rdb_roquet_spv.inc"],
    ("eos", "linear"): list(_EOS),
    ("geometry", "channel"): ["src/core/ocean/boundary/rdb_ocean_periodic.F90",
                              "src/core/ocean/boundary/rdb_ocean_sponge.F90"],
    ("geometry", "obc"): ["src/core/ocean/boundary/rdb_ocean_obc*.F90",
                          "src/core/ocean/boundary/rdb_ocean_boundary_data.F90"],
    ("geometry", "tripolar"): ["src/core/ocean/boundary/rdb_ocean_fold*.F90",
                               "src/core/ocean/state/rdb_ocean_bipolar.F90"],
    ("geometry", "cavity"): ["src/core/ocean/state/rdb_ocean_cavity.F90",
                             _P + "vertical/rdb_ocean_cavity_*.F90",
                             _P + "vertical/rdb_ocean_top_drag.F90"],
    ("forcing", "cool"): [_P + "vertical/rdb_ocean_surface_flux.F90"],
    ("forcing", "warm_sw"): [_P + "vertical/rdb_ocean_surface_flux.F90"],
})
for _v in ("const_nu_h", "smagorinsky", "leith", "leith_biharm", "nu_4", "stress_tensor", "kh_aniso"):
    VALUE_PATHS[("lateral", _v)] = list(_HVISC)
for _v in ("sadourny", "sadourny_energy", "sadourny_hk"):
    VALUE_PATHS[("coriolis", _v)] = list(_COR)
for _v in ("centered", "weno3", "weno5", "weno7"):
    VALUE_PATHS[("pv_adv", _v)] = list(_COR)
for _v in ("default", "correction_bc_pgf", "substep_drag", "wave_drag", "visc_rem"):
    VALUE_PATHS[("bt", _v)] = list(_SPLIT)
VALUE_PATHS[("bt", "visc_rem")] += [_P + "vertical/rdb_ocean_vdiff.F90"]
VALUE_PATHS[("bt", "substep_drag")] += [_P + "vertical/rdb_ocean_bottom_drag.F90"]

# Paths that are the matrix itself: any change runs everything.
ALWAYS_FULL = ["tests/regression/compat_*.py", "app/main.F90", "src/driver/*"]


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
        saved = cm.BACKEND[0]
        try:
            # judged on the row's own backend; invisible on any other
            cm.BACKEND[0] = (r.backend + "-selftest") if r.backend else saved
            rec = {"axes": cell, "class": "PASS", "rows": []}
            if (r.rid, "PASS") not in cm.row_xpasses([rec]):
                bad.append(r.rid + ": a passing witness is not a row XPASS")
            rec = {"axes": cell, "class": "XFAIL", "rows": [r.rid]}
            if any(rid == r.rid for rid, _ in cm.row_xpasses([rec])):
                bad.append(r.rid + ": a failing witness is a row XPASS")
            if r.backend:
                cm.BACKEND[0] = "other-backend"
                rec = {"axes": cell, "class": "PASS", "rows": []}
                if any(rid == r.rid for rid, _ in cm.row_xpasses([rec])):
                    bad.append(r.rid + ": a backend row judged on another backend")
        finally:
            cm.BACKEND[0] = saved
    return not bad, "scope='any' witnesses are real, matching, and XPASS when they pass {}".format(bad)


def _t_legs_contract(cm):
    """A leg row (RESTART) is judged only when its leg ran; a multirank row
    turns a skipped DECOMP into a PASS and a DECOMP that ran into an XPASS."""
    global ROWS
    saved = ROWS
    try:
        ROWS = [_gap("t_rst", "runtime", ("closed_faces",), "t", 0, expect=("RESTART",),
                     message=r"^synthetic restart"),
                _single_rank("t_mr", ("closed_faces",), "t", 0, message=r"^synthetic single")]
        cell = dict(cm.BASE_CELL)
        nml = cm.merged_namelist(cell)
        acc = {"status": "accepted", "rc": 0, "stage": None, "messages": []}
        ok = {"outcome": "PASS", "detail": ""}
        legs = cm.PHASE_A + ("restart", "decomp")
        got = [cm.classify(cell, nml, acc, ok)[0],                       # restart not run
               cm.classify(cell, nml, acc, {"outcome": "RESTART", "detail": "synthetic restart"},
                           ran=legs, decomp={"status": "skipped", "rows": ["t_mr"]})[0],
               cm.classify(cell, nml, acc, ok, ran=legs,
                           decomp={"status": "skipped", "rows": ["t_mr"]})[0],
               cm.classify(cell, nml, acc, ok, ran=cm.PHASE_A + ("decomp",),
                           decomp={"status": "skipped", "rows": ["t_mr"]})[0],
               cm.classify(cell, nml, acc, ok, ran=cm.PHASE_A + ("decomp",),
                           decomp={"status": "accepted", "rows": []})[0]]
        want = ["PASS", "XFAIL", "XPASS", "PASS", "XPASS"]
        return got == want, ("leg rows (not run / xfail / xpass) and multirank rows "
                             "(skipped / decomposes) {}".format(got))
    finally:
        ROWS = saved


SELF_TESTS = [("unexplained", _t_unexplained), ("validate_crash", _t_validate_crash),
              ("legs_contract", _t_legs_contract),
              ("runtime_unexpected", _t_runtime_unexpected), ("rows_contract", _t_rows_contract),
              ("witnesses", _t_witnesses)]
