"""THE REST-STATE MATRIX -- the same problem under EVERY vertical coordinate.

Why this exists
===============
`stability_manifest.py` runs each shipped namelist under whatever vertical
coordinate that namelist happens to choose.  Nothing in the tree runs ONE
problem under ALL of them, so nothing in the tree can answer the question a
user actually has: *which coordinate is trustworthy on which geometry?*

That gap is not theoretical.  In two days of work on the cavity branch it
hid, simultaneously:

  * `VCOORD_ZSIGMA` collapsing its whole column into the bed layer (the
    deep branch reads `z_ref_global` as METRES; its only writer fills it
    with the dimensionless `k/nz`) -- with `sum target_h = H + eta` still
    exact, which is why no conservation test caught it;
  * `VCOORD_ZSTAR_SIGMA` being numerically indistinguishable from
    `VCOORD_SIGMA` in all 21 shipped namelists that select it;
  * `VCOORD_ZSTAR_FULL` resolving the WRONG HALF of the column under a lid;
  * a rest-state growth mode over ANY sloping boundary that no shipped case
    saw, because every sloping case carries viscosity or runs short;
  * sigma failing outright on steep bathymetry (`rx0 = 0.73`).

Every one of those is invisible to a per-namelist suite and obvious in a
family x geometry table.  Producing that table is what this module is for --
and the table itself, recorded in `tests/regression/README.md`, IS the
product.

How it works
============
ONE canonical template (`vcoord_templates/rest_matrix.nml.in`) plus one
parameter dict per PROBLEM and one per FAMILY.  `build_matrix()` emits a
namelist per cell into `tmp_local_artifacts/vcoord_matrix/` at import time
and returns manifest rows built through `stability_manifest._case`.  The
emitted files are plain text and are left on disk after a run, so any cell
can be reproduced by hand with `./rdb <that file>`.

Dozens of near-identical namelists are deliberately NOT checked in: they
rot, a reader cannot tell which cell of the matrix a given file is, and a
change to the shared problem then has to be applied by hand N times.

What a row asserts
==================
The problem is AT REST with no energy source of any kind (see the template
header), so:

  * `energy:rest`            -- a velocity bar, per geometry (below);
  * `energy:rest-growth-rate`-- the fitted EXPONENTIAL RATE of `En`, which
                                is the assertion the existing suite was
                                missing.  A spurious pressure gradient that
                                EQUILIBRATES is a tolerable discretisation
                                error; one that is still exponential when
                                the clock runs out is an instability wearing
                                a small number, and on the cavity case a
                                10 900x smaller seed bought only 143 days
                                before it breached the same bar (see
                                `validation_examples/ocean/
                                ice_shelf_cavity/README.md`).  A LEVEL bar
                                cannot tell those apart; a RATE bar can.
  * `conserve:{Mass,Salt,Heat}`, `finite`, `cfl:bounded`, `cfl:no-runaway`
                              -- the existing gates, unchanged;
  * `tracer:no-new-extrema`   -- at rest, with no source and no mixing, the
                                T and S extrema may not leave their initial
                                range.  Spurious diapycnal mixing shows here
                                before it shows anywhere else.
  * `thickness:positive`      -- the minimum `h_layer` over the run stays
                                strictly positive (read off the `h_layer`
                                derived diagnostic's console `[diag]` line).
  * `counters:no-truncation`  -- a rest case that CFL-truncates or clamps is
                                not resting.

The bars, and where they come from
==================================
Never from what a run happens to do.  Three, all pre-existing in
`stability_manifest.py`:

  * FLAT geometry -> `REST_1UM_S` (1 um/s).  With no slope there is no
    spurious pressure gradient to generate and the only correct answer is
    machine zero.  Same bar as `seamount_flat` and `cavity_flat_lid_rest`,
    both measured bit-zero.
  * SLOPED geometry -> `REST_1MM_S` (1 mm/s).  The velocity a resting ocean
    has no excuse to exceed; real sub-shelf and shelf-break flows are cm/s,
    so 1 mm/s of spurious current is already a serious contaminant.  Same
    bar as `cavity_sloping_lid_rest`.
  * SEAMOUNT geometry -> `REST_SEAMOUNT` (1 cm/s).  The seamount test's own
    conventional acceptance bar in the terrain-following literature
    (Beckmann & Haidvogel 1993, J. Phys. Oceanogr. 23, 1736-1753; Haney
    1991, J. Phys. Oceanogr. 21, 610-619; Mellor, Ezer & Oey 1994).

and one new one:

  * `REST_SIGMA_MAX` -- the GROWTH-RATE bar, 0.05 per day of `En`, i.e. a
    20-day e-folding.  Derived, not tuned: the sloping-boundary sigma
    instability this matrix exists to find runs at `sigma_En = 0.333/day`
    (a 3.0-day e-folding) on the measured cavity case and never below
    0.094/day anywhere on its slope ladder, so the bar sits 1.9x under the
    slowest measured instance of the defect; the residual creep the DEFAULT
    outer split leaves behind is 0.012/day (an 83-day e-folding), so the bar
    sits 4x above what a healthy-but-imperfect run does.  The gap between
    those two is one decade, which is what makes a single bar workable.
    Expressed in 1/s of AMPLITUDE (`En ~ exp(2 sigma t)`) to match
    `stability.py`'s convention.

The truncation estimate, stated in the case docstrings
======================================================
For a terrain-following family the PLATEAU the spurious velocity settles at
is not arbitrary -- it is the second-order pressure-gradient truncation,

    a_peak = N^2 * de^3 / (6 * dx * Hbar),        de = |H_a - H_b| per face

(derived in `validation_examples/ocean/ice_shelf_cavity/README.md` and
verified there against a four-decade slope ladder to within 33%), which a
non-rotating run balances against nothing and a rotating one balances
geostrophically at `U = a_peak / |f|`, i.e.

    En_plateau ~ 0.5 * (a_peak / |f|)^2.

`problem_truncation_estimate()` evaluates that for each geometry so the row
carries the number it OUGHT to sit at beside the number it does.  It is an
estimate and is labelled as one: it assumes a single dominant face step and
a geostrophic balance, and the measured cavity case tracks it to within a
factor of ~2.

NOT DONE, and not stubbed
=========================
  * The NONLINEAR (exponential thermocline) stratification axis.
    `&ocean_zinit_nml source="linear"` is AFFINE in z by construction and
    `source="file"` aborts at `validate_config` (PR-23b), so there is no way
    to lay a non-affine rest state today without either a new analytic
    profile or the file reader.  The axis is absent rather than faked.

Stdlib only.
"""

import math
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
TEMPLATE = os.path.join(THIS_DIR, "vcoord_templates", "rest_matrix.nml.in")
OUT_REL = os.path.join("tmp_local_artifacts", "vcoord_matrix")
OUT_DIR = os.path.join(REPO_ROOT, OUT_REL)

# --- the shared problem geometry -------------------------------------------
NX, NY, NZ = 48, 6, 15
DX = DY = 2000.0
NGHOST = 2
DT = 600.0
CORIOLIS_F = -1.409e-4          # f-plane at 75 S (ISOMIP+ 3.1.1)
MAX_DEPTH = 1000.0
CAVITY_DEPTH = 720.0            # ISOMIP+ z_b,deep

# --- the stratification, laid FLAT IN GEOPOTENTIAL z -----------------------
# S falls linearly from 33.8 PSU at z = 0 to 34.55 PSU at z = -1000 m, with T
# uniform, so the whole density signal is in S and the linear EOS makes
# N^2 = g * beta_S * |dS/dz| / rho_0 = 5.77e-06 s^-2 exactly.  ISOMIP+
# Table 4 water masses (Asay-Davis et al. 2016); the same profile the shipped
# cavity rest cases use, rescaled from their 720 m column to 1000 m.
LIN_T_REF = -1.9
LIN_S_REF = 33.8
LIN_DS_DZ = -7.5e-4             # PSU/m, z positive UP
G = 9.81
RHO_0 = 1027.51
BETA_S = 8.0587609e-1           # kg/m^3/PSU (ISOMIP+ linear EOS)
N2_LINEAR = G * BETA_S * abs(LIN_DS_DZ) / RHO_0

STRATIFICATIONS = {
    # id      : (lin_dt_dz, lin_ds_dz, N^2, note)
    "linear": (0.0, LIN_DS_DZ, N2_LINEAR,
               "uniform N^2 = {:.3e} s^-2 from a linear S(z)".format(N2_LINEAR)),
    "unstrat": (0.0, 0.0, 0.0,
                "N^2 = 0. The sigma truncation G(K) goes as rho_0*N^2*de^3, "
                "so at N^2 = 0 it is not small but ABSENT however steep the "
                "geometry -- whatever a family produces here is load "
                "bookkeeping plus rounding, and the control that says the "
                "stratified answer really is the N^2 term."),
}

# --- the coordinate families ------------------------------------------------
# `status` is what the matrix EXPECTS of the family on a non-cavity problem:
#   "run"     -- runs, and is gated
#   "refused" -- `validate_config` refuses it by design; the row asserts the
#                REFUSAL, so an accidental un-refusal turns the row XPASS
FAMILIES = [
    # (id, vcoord_type, status, pins_scheme, extra vcoord_nml lines, note)
    ("lagrangian", "lagrangian", "run", None, "",
     "pure Lagrangian/isopycnal: target_h IS the live h_layer and the ALE "
     "remap is a no-op (early return), so this is the DATUM-FREE control leg "
     "-- the one family a displaced boundary cannot mis-place."),
    ("eulerian_z", "eulerian_z", "run", "ssp_rk2", "",
     "`H*dsig` with eta DROPPED -- a stretched sigma, not a geopotential "
     "coordinate despite the name. It pins split_scheme='ssp_rk2' because "
     "the pred_corr v1 envelope refuses eulerian_z fail-loud, which also "
     "keeps this row out of the automatic scheme axis."),
    ("sigma", "sigma", "run", None, "",
     "terrain-following. The reference leg: every sigma PGF-error result in "
     "the literature (Haney 1991; Beckmann & Haidvogel 1993; Mellor, Ezer & "
     "Oey 1994) is a statement about this row."),
    ("zstar", "zstar", "run", None, "",
     "z*-lite. On the ocean path it SHARES the sigma branch "
     "(`case (VCOORD_SIGMA, VCOORD_ZSTAR)`), so a row that differs from the "
     "sigma row is a defect in one of them, not a coordinate difference."),
    ("zstar_sigma", "zstar_sigma", "run", None, "",
     "sigma in shallow water, z*-lite in deep -- consuming `z_ref_global` "
     "FRACTIONALLY, so it is immune to the ZSIGMA units defect. With the "
     "uniform default table its deep branch is numerically "
     "indistinguishable from sigma, which is precisely what this row is "
     "here to make visible: 21 shipped namelists select it believing "
     "otherwise."),
    ("zstar_full", "zstar_full", "run", None,
     "   zstar_h_surf_target = 20.0\n   zstar_n_surf = 3",
     "per-column z_ref table from the local bathymetry, with a genuinely "
     "NON-UNIFORM stack (a 3 x 20 m fine surface band). The uniform default "
     "would make this row a duplicate of sigma and prove nothing."),
    ("z_fixed", "z_fixed", "run", None, "",
     "quasi-geopotential: fixed-z interfaces with bed-side layers vanishing "
     "to the inert filler. `z_fixed_h_ref` is taken from "
     "`&ocean_topo_nml max_depth`, so the nominal spacing is max_depth/nz."),
    ("rho", "rho", "run", None,
     "   rho_target_light = 1027.20\n   rho_target_dense = 1027.90\n"
     "   rho_ref_pressure = 0.0",
     "isopycnal: interfaces placed on prescribed potential densities by "
     "inverting a PPM reconstruction of the column density. The target band "
     "brackets the ISOMIP+ profile's actual range (1027.22..1027.83 kg/m^3 "
     "under the linear EOS), so the coordinate has something to resolve. "
     "Validation-grade: weakly-stratified columns collapse, which is why "
     "the `unstrat` cell of this row is expected to be degenerate."),
    ("hycom", "hycom", "run", None,
     "   rho_target_light = 1027.20\n   rho_target_dense = 1027.90\n"
     "   rho_ref_pressure = 0.0",
     "hybrid z*/isopycnal (Bleck 2002, MOM6 coord_hycom): the same density "
     "inversion plus a bottom-up monotonize before it and a z* nominal-floor "
     "sweep after. The production GVC coordinate."),
    ("zsigma", "zsigma", "refused", None, "",
     "REFUSED at configure on EVERY path: the deep branch reads "
     "`z_ref_global` as absolute depths in METRES while its only writer "
     "fills it with the dimensionless `k/nz`, so every z-level interval is "
     "`1/nz` metres and the whole column collapses into the bed layer -- "
     "with `sum target_h = H + eta` still exact, which is why no "
     "conservation test caught it. This row asserts the REFUSAL: if it ever "
     "turns XPASS, either the units defect is fixed (delete the row) or the "
     "refusal was dropped without fixing it (a regression)."),
]

# Families the CAVITY (`&ocean_cavity_dyn_nml enable`) accepts.  Everything
# else is refused fail-loud by `validate_config`, each with its own reason --
# the two that rescale the live column and so follow the ice base for free,
# plus the one that has been TAUGHT the ice base.
CAVITY_ACCEPTED = ("sigma", "zstar", "z_fixed")

# --- the geometries ---------------------------------------------------------
# `bar` names the `en_rest_max` constant; `de` is the worst per-face depth
# step the geometry presents, used for the truncation estimate.
PROBLEMS = {
    "flat": {
        "class": "flat", "bar": "REST_1UM_S", "de": 0.0, "hbar": MAX_DEPTH,
        "topo": "flat", "tier2": True,
        "doc": "flat bed, no lid. Every per-face interface offset de(K) is "
               "identically ZERO, so every terrain-following truncation term "
               "is structurally ABSENT and the only correct answer is exactly "
               "no motion. The control for the whole matrix: a family that "
               "moves HERE has a defect that has nothing to do with "
               "topography.",
    },
    "slope": {
        "class": "slope", "bar": "REST_1MM_S",
        "de": (MAX_DEPTH - 500.0) / (NX - 1), "hbar": 750.0,
        "topo": "file", "bathy": ("slope", {"deep": MAX_DEPTH, "shallow": 500.0}),
        "tier2": True,
        "doc": "a GENTLE constant-gradient bed, 1000 m -> 500 m across 48 "
               "cells. Every wet-wet face carries the same depth step "
               "(10.64 m) and so the same interface offset, which is what "
               "makes this the clean geometry for the de^3 truncation "
               "scaling -- a seamount convolves the scaling with the "
               "curvature of its own slope.",
    },
    "seamount_gentle": {
        "class": "seamount", "bar": "REST_SEAMOUNT", "de": 30.0, "hbar": 850.0,
        "topo": "seamount", "edge_depth": 300.0, "slope_scale": 40000.0,
        "tier2": False,
        "doc": "Gaussian seamount, e-folding scale 40 km = 20 cells. The "
               "Beckmann & Haidvogel (1993) test geometry at its GENTLE "
               "setting: they smooth their seamount to rx0 <= 0.2 and this "
               "sits under that.",
    },
    "seamount_steep": {
        "class": "seamount", "bar": "REST_SEAMOUNT", "de": 75.0, "hbar": 800.0,
        "topo": "seamount", "edge_depth": 300.0, "slope_scale": 15000.0,
        "tier2": False,
        "doc": "the same Gaussian seamount at 15 km = 7.5 cells, i.e. past "
               "the Beckmann & Haidvogel stiffness bound. Paired with the "
               "gentle twin this is a two-point de^3 check on every family "
               "at once.",
    },
}

# The rx0 LADDER -- a two-level shelf/trough whose terrain-following
# stiffness is dialled exactly.  `rdb_ocean_stability_audit` bounds
# rx0 = |dH|/(H_a + H_b) at 0.2 (Beckmann & Haidvogel 1993 section 2c) and
# WARNS above it; the ladder walks past the bound so the matrix can say where
# each family stops being trustworthy rather than asserting the bound.
# The motivating failure sits at 0.73 (ISOMIP+ trough sidewall), between the
# 0.6 and 0.8 rungs.
RX0_HBAR = 750.0
RX0_RUNGS = (0.1, 0.2, 0.4, 0.6, 0.8)
for _r in RX0_RUNGS:
    PROBLEMS["rx0_{:03d}".format(int(round(_r * 100)))] = {
        "class": "rx0", "bar": "REST_1MM_S",
        "de": 2.0 * RX0_HBAR * _r, "hbar": RX0_HBAR,
        "topo": "file", "bathy": ("rx0_ladder", {"hbar": RX0_HBAR, "rx0": _r}),
        # `max_depth` must be the DEEPEST column the geometry contains, not
        # the mean: `z_fixed_h_ref` is taken from it, so a nominal stack cut
        # for 750 m would leave a 1350 m column unable to reach its own bed.
        "max_depth": RX0_HBAR * (1.0 + _r),
        "tier2": (abs(_r - 0.6) < 1e-9),
        "rx0": _r,
        "doc": "rx0 ladder rung {:.1f}: a single wet-wet face joining a "
               "{:.0f} m column to a {:.0f} m one, so rx0 = |dH|/(H_a+H_b) "
               "= {:.1f} EXACTLY at that face and 0 everywhere else. The "
               "Beckmann & Haidvogel (1993) sigma bound is 0.2; the ISOMIP+ "
               "trough sidewall that takes the shipped case non-finite at "
               "day 3.2 is 0.73.".format(
                   _r, RX0_HBAR * (1 + _r), RX0_HBAR * (1 - _r), _r),
    }

# The CAVITY geometries -- problems 5 and 6.  A rigid ice lid instead of a
# free surface, which is the SAME question with the tilted boundary at the
# other end of the column (proved to four significant figures by the mirror
# experiment in the cavity diagnosis: an ice-free domain with the geometry
# turned upside down reproduces the cavity case at every sample).
PROBLEMS["lid_flat"] = {
    "class": "cavity", "bar": "REST_1UM_S", "de": 0.0, "hbar": 500.0,
    "topo": "flat", "max_depth": CAVITY_DEPTH, "tier2": False,
    "cavity": {"draft_config": "flat", "draft_depth": 220.0},
    "doc": "a FLAT ice lid over a flat bed: every interface gap is zero at "
           "BOTH boundaries, so like `flat` the only correct answer is "
           "exactly no motion -- but it additionally exercises the ice-load "
           "cancellation (pa(nz+1) = rho_ref*g*eta_geo + rho_ref*g*z_draft, "
           "the same product twice with opposite signs) and the "
           "bt_H_ref = b - z_draft datum.",
}
PROBLEMS["lid_slope"] = {
    "class": "cavity", "bar": "REST_1MM_S", "de": 13.8, "hbar": 500.0,
    "topo": "flat", "max_depth": CAVITY_DEPTH, "tier2": False,
    "cavity": {"draft_config": "linear", "draft_depth": 570.0,
               "draft_slope": -6.9e-3, "draft_x0": -10000.0,
               "draft_x1": 72000.0},
    "doc": "a linearly SLOPING ice lid with a calving front. The headline "
           "cavity geometry, and the one the Phase-6 PGF corrections (Yung, "
           "Hallberg, Adcroft & Morrison 2026, JAMES 18, e2025MS005645) "
           "exist for. Carried here so every coordinate family is measured "
           "on it, not only the sigma the shipped case pins.",
}



# ---------------------------------------------------------------------------
# THE BASELINE TABLE -- measured, then pinned
# ---------------------------------------------------------------------------
# Every cell below FAILS today.  Each entry is
#     "problem/family[/strat][/eos]": (reason-key, [assertions], measured En)
# and is MEASURED, never chosen: the table was produced by running the whole
# matrix at 3.33 simulated days on gfortran 15.1 Release, single rank, MPI
# off (`tmp_local_artifacts/vcm_measure`), and the numbers below are what
# came out.  Nothing here was tuned to pass; nothing was excused wholesale --
# each marker names the ASSERTIONS it covers, so a cell that is XFAIL on
# `energy:rest-settles` still gates `conserve:*`, `finite` and the level bar.
#
# WHEN A CELL FLIPS.  A fix that removes one of these turns the row XPASS
# and the suite says so.  The two dyn-core fixes that live on other branches
# and are NOT in this tree are named in the reasons below:
#   * `origin/fix/sigma-pgf-rest-state`      -- the exact FV pressure-gradient
#     under `reconstruct_for_pressure`, which is what the `slope` and
#     `seamount_*` plateaus are measuring;
#   * `origin/fix/ale-remap-rest-amplifier`  -- `&vcoord_nml
#     remap_boundary_extrap`, the first-order boundary cell in the remap,
#     which is the leading suspect for the `z_fixed` budget leak.
# Both must be re-measured against this table when they land.
XFAIL_REASONS = {
    "blowup":
        "BLOWS UP. The run goes non-finite and aborts inside the first "
        "simulated day, from a state at rest with no energy source and no "
        "dissipation of any kind. This is not a tolerance question: the "
        "terrain-following pressure-gradient truncation over a single-face "
        "depth step of this size produces a spurious acceleration the "
        "barotropic mode cannot absorb, and NO shipped coordinate family "
        "survives it. Beckmann & Haidvogel (1993) bound the stiffness "
        "rx0 = |dH|/(H_a+H_b) at 0.2 for exactly this reason; the "
        "`rdb_ocean_stability_audit` warning fires at configure with the "
        "number. The fix is to smooth the topography or to implement a "
        "pressure-gradient form that is exact on a step -- not a bar.",
    "rx0_plateau":
        "SURVIVES, BUT DOES NOT REST. The run completes and its budgets "
        "close at round-off, and it develops centimetres per second of "
        "spurious current out of nothing and is still at its maximum when "
        "the clock stops. A forced, viscous configuration would never "
        "notice; a quiescent or long spin-up one measures this instead of "
        "the physics. Same mechanism as the `blowup` rungs, one or two "
        "decades weaker.",
    "z_fixed_leak":
        "LEAKS SALT AND HEAT, and makes new tracer extrema, on ANY geometry "
        "with vanishing layers -- a slope, a seamount or an ice base. The "
        "budget residual is 1e-6 RELATIVE against a 1e-11 bar, i.e. five "
        "decades over, and it is a step at the first regrid rather than a "
        "drift, which points at the ALE drain of the inert filler layers "
        "rather than at transport. `z_fixed` is the one family whose "
        "shallow columns carry a full stack of fillers, and it is the only "
        "family in the matrix that leaks. Suspect: the first-order boundary "
        "cell in the remap reconstruction "
        "(`origin/fix/ale-remap-rest-amplifier`, `&vcoord_nml "
        "remap_boundary_extrap`), which this tree does not carry -- but "
        "that is a hypothesis, not a measurement, and this marker must be "
        "re-measured when that branch lands.",
    "density_coord":
        "A DENSITY-SPACE COORDINATE ON A GEOMETRIC REST STATE. `rho` and "
        "`hycom` place their interfaces on prescribed potential densities, "
        "so on a linearly stratified column over a slope the regrid moves "
        "every interface every step and the remap error it pays is three to "
        "four decades above the geometric families on the identical "
        "problem. `rho` is documented validation-grade for exactly this "
        "reason; `hycom` is the production GVC coordinate and its number "
        "here is the one to watch. Neither is refused, so this is the "
        "envelope statement rather than a fence.",
    "zstar_full":
        "z*-FULL's PER-COLUMN TABLE. The family builds a separate `z_ref` "
        "stack per column from the LOCAL bed, so two neighbouring columns "
        "of different depth get different interface depths and the "
        "resulting offset is larger than sigma's on the same geometry -- "
        "measured 3.2x (slope) to 46x (seamount_gentle) above the sigma "
        "cell. It is the only geometric family whose interface offset does "
        "not shrink with the bathymetric gradient.",
    "seamount_settle":
        "DOES NOT EQUILIBRATE WITHIN THE RUN. The LEVEL is inside the "
        "seamount bar (Beckmann & Haidvogel 1993's conventional 1 cm/s) and "
        "the budgets close, but the spurious energy is still at its maximum "
        "when the run ends, so nothing has bounded it yet. Scoped to the "
        "SETTLE gate alone for that reason: the magnitude gate is live and "
        "passes. Whether this is a slow approach to a plateau or an "
        "instability is what the tier-1 30-day horizon and the "
        "`energy:rest-growth-rate` fit are there to answer.",
    "cavity_settle":
        "THE SLOPING-LID SIGMA TRUNCATION, at matrix scale -- the same "
        "defect `validation_examples/ocean/ice_shelf_cavity/"
        "cavity_sloping_lid_rest.nml` carries as a scoped XFAIL, and for "
        "the same reason: this build implements none of the "
        "sloping-surface pressure-gradient corrections of Yung, Hallberg, "
        "Adcroft & Morrison (2026), JAMES 18, e2025MS005645, so the "
        "spurious energy leaves its plateau on a ~3-day e-folding and the "
        "final sample is the peak at every horizon short of saturation. "
        "Scoped to the SETTLE gate: the MAGNITUDE gate passes with margin "
        "and is deliberately left live, because it is what separates the "
        "two outer split schemes.",
    "slope_plateau":
        "THE TERRAIN-FOLLOWING PLATEAU, AND WHETHER IT IS ONE. Over a "
        "constant-gradient bed or a Gaussian seamount the spurious energy "
        "sits at the second-order pressure-gradient truncation "
        "a_peak = N^2 de^3/(6 dx Hbar) -- which is a tolerable "
        "discretisation error IF it equilibrates. At the tier-2 horizon "
        "(3.33 days) several of these cells are still at their maximum when "
        "the clock stops, and at the tier-1 horizon (30 days) the fitted "
        "growth rate says which of them are actually still growing. That is "
        "the whole reason the RATE gate exists: a level bar reads a "
        "slow instability and a settled truncation identically. The fix is "
        "a pressure-gradient form that is exact on a sloping coordinate "
        "surface -- `origin/fix/sigma-pgf-rest-state` carries the exact FV "
        "form under `reconstruct_for_pressure` and is NOT in this tree, so "
        "these cells must be re-measured when it lands.",
    "generic":
        "Measured failing on this cell; see the baseline table in "
        "tests/regression/README.md.",
}

MEASURED_XFAIL = {
    "lid_slope/sigma":
        ("cavity_settle", ['energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 1.988e-09 (final 1.988e-09); tier2 (3.33 d) peak En 1.02e-09 (final 1.004e-09)"),
    "lid_slope/z_fixed":
        ("blowup", ['completed', 'conserve:Heat', 'conserve:Salt', 'energy:rest', 'finite', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 3.478e-05 (final 2.909e-05)"),
    "lid_slope/zstar":
        ("cavity_settle", ['energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 1.988e-09 (final 1.988e-09); tier2 (3.33 d) peak En 1.02e-09 (final 1.004e-09)"),
    "rx0_010/eulerian_z":
        ("blowup", ['completed', 'energy:rest', 'energy:rest-settles', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 1.954e-06 (final 1.954e-06)"),
    "rx0_010/hycom":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_010/lagrangian":
        ("blowup", ['completed', 'energy:rest', 'energy:rest-settles', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 0.0001066 (final 0.0001066)"),
    "rx0_010/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_010/sigma":
        ("rx0_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.000902 (final 0.000902); tier2 (3.33 d) peak En 1.445e-06 (final 1.445e-06)"),
    "rx0_010/z_fixed":
        ("blowup", ['completed', 'conserve:Heat', 'conserve:Salt', 'energy:rest', 'energy:rest-settles', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 0.0007616 (final 0.0007616)"),
    "rx0_010/zstar":
        ("rx0_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.000902 (final 0.000902); tier2 (3.33 d) peak En 1.445e-06 (final 1.445e-06)"),
    "rx0_010/zstar_full":
        ("rx0_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.0007728 (final 0.0007728); tier2 (3.33 d) peak En 2.261e-06 (final 2.261e-06)"),
    "rx0_010/zstar_sigma":
        ("rx0_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.000902 (final 0.000902); tier2 (3.33 d) peak En 1.445e-06 (final 1.445e-06)"),
    "rx0_020/eulerian_z":
        ("blowup", ['completed', 'energy:rest', 'energy:rest-settles', 'finite', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 0.0001767 (final 0.0001767)"),
    "rx0_020/hycom":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_020/lagrangian":
        ("rx0_plateau", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_020/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_020/sigma":
        ("rx0_plateau", ['energy:rest', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.006354 (final 0.006222); tier2 (3.33 d) peak En 0.0001335 (final 0.0001335)"),
    "rx0_020/z_fixed":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_020/zstar":
        ("rx0_plateau", ['energy:rest', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.006354 (final 0.006222); tier2 (3.33 d) peak En 0.0001335 (final 0.0001335)"),
    "rx0_020/zstar_full":
        ("rx0_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.005049 (final 0.004983); tier2 (3.33 d) peak En 0.0001344 (final 0.0001344)"),
    "rx0_020/zstar_sigma":
        ("rx0_plateau", ['energy:rest', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.006358 (final 0.006249); tier2 (3.33 d) peak En 0.0001335 (final 0.0001335)"),
    "rx0_040/eulerian_z":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/hycom":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/lagrangian":
        ("rx0_plateau", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/z_fixed":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/zstar":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/zstar_full":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_040/zstar_sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/eulerian_z":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/hycom":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/lagrangian":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/z_fixed":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/zstar":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/zstar_full":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_060/zstar_sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/eulerian_z":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/hycom":
        ("blowup", ['completed'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/lagrangian":
        ("blowup", ['completed'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/sigma/unstrat":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/z_fixed":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/zstar":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/zstar_full":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "rx0_080/zstar_sigma":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "seamount_gentle/eulerian_z":
        ("slope_plateau", ['energy:rest', 'energy:rest-growth-rate', 'tracer:no-new-extrema'], [1], "tier1 (30 d) peak En 0.0005058 (final 0.0004234); tier2 (3.33 d) peak En 9.924e-09 (final 8.983e-09)"),
    "seamount_gentle/hycom":
        ("density_coord", ['energy:rest', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 7.601e-05 (final 7.563e-05); tier2 (3.33 d) peak En 2.374e-05 (final 2.361e-05)"),
    "seamount_gentle/lagrangian":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1], "tier1 (30 d) peak En 8.143e-08 (final 8.143e-08); tier2 (3.33 d) peak En 1.025e-08 (final 8.508e-09)"),
    "seamount_gentle/rho":
        ("blowup", ['completed', 'energy:rest-settles', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 1.808e-05 (final 1.808e-05)"),
    "seamount_gentle/sigma":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 4.911e-06 (final 4.911e-06); tier2 (3.33 d) peak En 1.048e-08 (final 1.006e-08)"),
    "seamount_gentle/z_fixed":
        ("blowup", ['completed', 'conserve:Heat', 'conserve:Salt', 'energy:rest', 'energy:rest-settles', 'finite', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 7.392e-05 (final 7.255e-05)"),
    "seamount_gentle/zstar":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 4.911e-06 (final 4.911e-06); tier2 (3.33 d) peak En 1.048e-08 (final 1.006e-08)"),
    "seamount_gentle/zstar_full":
        ("zstar_full", ['energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 2.481e-05 (final 2.481e-05); tier2 (3.33 d) peak En 4.877e-07 (final 4.877e-07)"),
    "seamount_gentle/zstar_sigma":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 4.911e-06 (final 4.911e-06); tier2 (3.33 d) peak En 1.048e-08 (final 1.006e-08)"),
    "seamount_steep/eulerian_z":
        ("blowup", ['completed', 'energy:rest-settles', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 1.599e-06 (final 1.599e-06)"),
    "seamount_steep/hycom":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "seamount_steep/lagrangian":
        ("slope_plateau", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "seamount_steep/rho":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "seamount_steep/sigma":
        ("slope_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.0001186 (final 0.000116); tier2 (3.33 d) peak En 1.084e-06 (final 1.084e-06)"),
    "seamount_steep/z_fixed":
        ("blowup", ['completed', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En aborted"),
    "seamount_steep/zstar":
        ("slope_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.0001186 (final 0.000116); tier2 (3.33 d) peak En 1.084e-06 (final 1.084e-06)"),
    "seamount_steep/zstar_full":
        ("zstar_full", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.000268 (final 0.000268); tier2 (3.33 d) peak En 2.402e-06 (final 2.402e-06)"),
    "seamount_steep/zstar_sigma":
        ("slope_plateau", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 0.0001186 (final 0.000116); tier2 (3.33 d) peak En 1.084e-06 (final 1.084e-06)"),
    "slope/eulerian_z":
        ("slope_plateau", ['energy:rest', 'tracer:no-new-extrema'], [1], "tier1 (30 d) peak En 0.001131 (final 0.001028); tier2 (3.33 d) peak En 8.391e-11 (final 6.29e-11)"),
    "slope/hycom":
        ("density_coord", ['energy:rest', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 1.721e-05 (final 1.689e-05); tier2 (3.33 d) peak En 1.094e-05 (final 1.094e-05)"),
    "slope/hycom/wright":
        ("blowup", ['completed', 'energy:rest'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 1.49e-05 (final 1.35e-05)"),
    "slope/rho":
        ("blowup", ['completed', 'energy:rest', 'finite'], [1, 2], "tier1 (30 d) peak En aborted; tier2 (3.33 d) peak En 1.298e-05 (final 1.203e-05)"),
    "slope/sigma":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1], "tier1 (30 d) peak En 2.995e-09 (final 2.995e-09); tier2 (3.33 d) peak En 8.85e-11 (final 6.584e-11)"),
    "slope/sigma/wright":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1], "tier1 (30 d) peak En 3.246e-09 (final 3.246e-09); tier2 (3.33 d) peak En 8.982e-11 (final 6.652e-11)"),
    "slope/z_fixed":
        ("slope_plateau", ['conserve:Heat', 'conserve:Salt', 'energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En 3.23e-05 (final 3.23e-05); tier2 (3.33 d) peak En 5.109e-06 (final 5.109e-06)"),
    "slope/z_fixed/wright":
        ("z_fixed_leak", ['conserve:Heat', 'conserve:Salt', 'energy:rest', 'energy:rest-settles', 'tracer:no-new-extrema'], [1, 2], "tier1 (30 d) peak En 2.936e-05 (final 2.936e-05); tier2 (3.33 d) peak En 4.755e-06 (final 4.755e-06)"),
    "slope/zstar":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1], "tier1 (30 d) peak En 2.995e-09 (final 2.995e-09); tier2 (3.33 d) peak En 8.85e-11 (final 6.584e-11)"),
    "slope/zstar_full":
        ("zstar_full", ['energy:rest', 'energy:rest-growth-rate', 'energy:rest-settles'], [1, 2], "tier1 (30 d) peak En 4.501e-06 (final 4.501e-06); tier2 (3.33 d) peak En 2.836e-07 (final 2.836e-07)"),
    "slope/zstar_sigma":
        ("slope_plateau", ['energy:rest-growth-rate', 'energy:rest-settles'], [1], "tier1 (30 d) peak En 2.995e-09 (final 2.995e-09); tier2 (3.33 d) peak En 8.85e-11 (final 6.584e-11)"),
}

def problem_truncation_estimate(prob, n2):
    """Derived plateau estimate for a terrain-following family, m2/s2.

    `a_peak = N^2 * de^3 / (6 * dx * Hbar)` is the second-order
    pressure-gradient truncation a sigma column carries across a face whose
    two interface positions differ by `de` (derived in
    `validation_examples/ocean/ice_shelf_cavity/README.md`; verified there
    against a four-decade slope ladder to within 33%).  A rotating run
    balances it geostrophically at `U = a_peak / |f|`, so

        En_plateau ~ 0.5 * U^2.

    Returns `(a_peak, en_plateau)`, both 0 when `de` or `N^2` is 0 -- which
    is a STATEMENT, not a degenerate case: at `N^2 = 0` the term is absent
    however steep the geometry, and at `de = 0` there is nothing to truncate.
    """
    de = prob["de"]
    if de <= 0.0 or n2 <= 0.0:
        return 0.0, 0.0
    a_peak = n2 * de ** 3 / (6.0 * DX * prob["hbar"])
    u = a_peak / abs(CORIOLIS_F)
    return a_peak, 0.5 * u * u


# ---------------------------------------------------------------------------
# Namelist emission
# ---------------------------------------------------------------------------
def _fmt(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, float):
        return repr(v)
    return str(v)


def render(name, prob, fam, strat, eos, remap="ppm"):
    """Render one matrix cell's namelist text from the canonical template."""
    fam_id, vcoord, status, pin, extra, _note = fam
    dt_dz, ds_dz, _n2, _sn = STRATIFICATIONS[strat]
    max_depth = prob.get("max_depth", MAX_DEPTH)
    cav = prob.get("cavity")

    bathy_block = ""
    if prob["topo"] == "file":
        # `bathymetry_file` lives in `&output_nml`; the path is RELATIVE so it
        # resolves against the scratch dir the setup hook writes it into,
        # which is also the binary's cwd.
        bathy_block = '   bathymetry_file = "bathy.nc"'

    cavity_block = ""
    if cav:
        lines = ["&ocean_cavity_dyn_nml", "   enable = .true."]
        for k, v in sorted(cav.items()):
            lines.append("   {} = {}".format(
                k, '"{}"'.format(v) if isinstance(v, str) else _fmt(v)))
        # 96 m = 2*h_nominal (720/15), the ISOMIP+ minimum-column rule
        # (Asay-Davis et al. 2016 3.1.5) that `z_fixed` x cavity fails
        # loud on. Applied to EVERY family so the grounding line is in
        # the same place on every leg -- nothing is grounded at either
        # 40 or 96 in this geometry, so it changes no other row.
        lines.append("   h_min_cavity = 96.0")
        lines.append("/")
        cavity_block = "\n".join(lines)

    # The cavity path needs the ONE pressure-gradient form with an injectable
    # top boundary condition; everywhere else the matrix uses the default so
    # the coordinate is the only thing that changes across a row.
    pgf_form = "fv_mom6" if cav else "mont"
    p_top = ".true." if (cav and cav.get("draft_config") == "linear") else ".false."

    bt_extra = ""
    if pin:
        bt_extra = '   split_scheme = "{}"\n'.format(pin)

    subs = {
        "NAME": name, "NX": NX, "NY": NY, "NZ": NZ, "DX": DX, "DY": DY,
        "NGHOST": NGHOST, "DT": DT, "CORIOLIS_F": CORIOLIS_F,
        "VCOORD": vcoord, "REMAP": remap, "VCOORD_EXTRA": extra,
        "TOPO_CONFIG": prob["topo"], "MAX_DEPTH": max_depth,
        "EDGE_DEPTH": prob.get("edge_depth", max_depth),
        "SLOPE_SCALE": prob.get("slope_scale", 4.0e5),
        "BATHY_FILE": bathy_block, "CAVITY": cavity_block,
        "PGF_FORM": pgf_form, "P_TOP_IN_BC": p_top, "EOS": eos,
        "LIN_T_REF": LIN_T_REF, "LIN_DT_DZ": dt_dz,
        "LIN_S_REF": LIN_S_REF, "LIN_DS_DZ": ds_dz,
        "BT_EXTRA": bt_extra,
    }
    text = open(TEMPLATE).read()
    for k, v in subs.items():
        text = text.replace("{{" + k + "}}", _fmt(v) if not isinstance(v, str) else v)
    if "{{" in text:
        raise ValueError("{}: unsubstituted token in the template".format(name))
    return text


def _bathy_setup(prob):
    """The `setup` argv that writes this geometry's bathymetry NetCDF.

    Run by `stability.run_case` in the case's scratch dir, which is also the
    binary's cwd -- so the namelist's relative `bathymetry_file = "bathy.nc"`
    resolves to exactly this file.
    """
    kind, kw = prob["bathy"]
    argv = ["{python}", "{repo}/tools/make_bathy_nc.py", "bathy.nc",
            "--nx", str(NX), "--ny", str(NY), "--profile", kind]
    for k, v in sorted(kw.items()):
        argv += ["--" + k, repr(v)]
    return argv


# ---------------------------------------------------------------------------
# Manifest rows
# ---------------------------------------------------------------------------
# Tier-1 (GPU, nightly) runs 30 simulated days; tier-2 (CPU, the CI gate)
# runs 3.3 days on a thin slice of the geometries.  The tier-2 length is NOT
# enough to fit a growth rate -- the sloping-boundary mode does not leave its
# plateau until day 25 at this resolution -- so the rate assertion is TIER 1
# ONLY and says so rather than being evaluated on noise.
T1_STEPS = 4320          # 30 days at dt = 600
T2_STEPS = 480           # 3.33 days
REST_SIGMA_MAX = 0.05 / 2.0 / 86400.0
    #: Amplitude growth-rate bar, 1/s.  0.05 per day of `En` = a 20-day
    #: e-folding; `En ~ exp(2 sigma t)` so the amplitude rate is half that.
    #: See the module header for the derivation.


def build_matrix(_case, bars, out_dir=None):
    """Emit every cell's namelist and return its manifest rows.

    `_case` is `stability_manifest._case` and `bars` the module's bar
    constants, both passed in rather than imported, so this module has NO
    import cycle with the manifest that consumes it.
    """
    out = out_dir or OUT_DIR
    os.makedirs(out, exist_ok=True)
    rows = []
    for pid in sorted(PROBLEMS):
        prob = PROBLEMS[pid]
        is_cavity = prob["class"] == "cavity"
        for fam in FAMILIES:
            fam_id, _vc, status, pin, _extra, fam_note = fam
            for strat in ("linear",):
                for eos in ("linear",):
                    rows.append(_one(_case, bars, out, pid, prob, fam, strat,
                                     eos, is_cavity, status, pin, fam_note))
        # The N^2 = 0 CONTROL, on the two geometries where it says the most:
        # the steepest ladder rung and the sloping lid.  It is a control, not
        # a coordinate test, so it runs on sigma alone -- the family whose
        # truncation term the control is isolating.
        if pid in ("rx0_080", "lid_slope"):
            fam = [f for f in FAMILIES if f[0] == "sigma"][0]
            rows.append(_one(_case, bars, out, pid, prob, fam, "unstrat",
                             "linear", is_cavity, "run", None, fam[5]))
        # The nonlinear-EOS leg, on the clean slope geometry: Wright (1997)
        # is the production EOS and its thermobaricity is exactly what a
        # linear-EOS rest state cannot see.
        if pid == "slope":
            for fam in FAMILIES:
                if fam[0] in ("sigma", "z_fixed", "hycom"):
                    rows.append(_one(_case, bars, out, pid, prob, fam,
                                     "linear", "wright", is_cavity, fam[2],
                                     fam[3], fam[5]))
    return rows


def _one(_case, bars, out, pid, prob, fam, strat, eos, is_cavity, status,
         pin, fam_note):
    """Build one matrix cell: write its namelist, return its manifest row."""
    fam_id = fam[0]
    parts = ["vcm", pid, fam_id]
    if strat != "linear":
        parts.append(strat)
    if eos != "linear":
        parts.append(eos)
    name = "_".join(parts)

    refused = (status == "refused") or (is_cavity and fam_id not in CAVITY_ACCEPTED)
    text = render(name, prob, fam, strat, eos)
    path = os.path.join(out, name + ".nml")
    with open(path, "w") as fh:
        fh.write(_header(name, pid, prob, fam, strat, eos, refused, fam_note))
        fh.write(text)

    _dt_dz, _ds_dz, n2, strat_note = STRATIFICATIONS[strat]
    a_peak, en_est = problem_truncation_estimate(prob, n2)
    bar_name = prob["bar"] if n2 > 0.0 else "REST_1UM_S"
    kw = {
        "en_rest_max": getattr(bars, bar_name),
        "tags": ["vcoord_matrix", "rest", "vcoord_" + fam_id, prob["class"]],
        "note": _row_note(pid, prob, fam, strat, eos, a_peak, en_est,
                          strat_note, bar_name, refused),
        "t1_timeout": 900, "t2_timeout": 300,
        "t1_samples": 60, "t2_samples": 30,
        "matrix_gates": True,
        "matrix": {"problem": pid, "family": fam_id,
                   "stratification": strat, "eos": eos,
                   "expect": "refused" if refused else "run",
                   "class": prob["class"], "rx0": prob.get("rx0"),
                   "de": prob["de"], "a_peak": a_peak,
                   "en_estimate": en_est},
    }
    if prob["topo"] == "file":
        kw["setup"] = _bathy_setup(prob)
    if not refused:
        # The RATE gate. Tier 1 only: a 3.3-day twin cannot fit an
        # exponential whose e-folding the bar puts at 20 days.
        kw["rest_sigma_max"] = REST_SIGMA_MAX
    else:
        kw["known_failure"] = {
            "assertions": ["completed"],
            "tiers": [1, 2],
            "reason": _refusal_reason(pid, prob, fam, is_cavity),
            "ref": "src/core/rdb_config.F90 :: validate_config",
        }
    # A cell that is MEASURED failing today carries its own scoped marker,
    # naming the assertions it covers and the number behind them. See
    # MEASURED_XFAIL -- the numbers are measurements, not choices.
    mkey = "/".join([pid, fam_id]
                    + ([strat] if strat != "linear" else [])
                    + ([eos] if eos != "linear" else []))
    if not refused and mkey in MEASURED_XFAIL:
        rkey, assertions, tiers, measured = MEASURED_XFAIL[mkey]
        kw["known_failure"] = {
            "assertions": list(assertions),
            "tiers": list(tiers),
            "reason": "{}  MEASURED on this cell -- {} (tier 1 = 30 "
                      "simulated days on a V100/nvfortran, tier 2 = 3.33 "
                      "days on gfortran 15.1 Release, single rank, "
                      "pred_corr). The marker is scoped to the TIERS the "
                      "defect is visible at as well as to the assertions it "
                      "covers: several of these only appear at the 30-day "
                      "horizon, and marking them known-failing at tier 2 "
                      "would report a permanent XPASS there. Do NOT close "
                      "this by widening en_rest_max, by shortening the run, "
                      "or by putting viscosity into the template -- a "
                      "viscosity that removes the spurious energy removes "
                      "the measurement with it (nu_h = 50 m2/s does exactly "
                      "that on the sibling cavity case).".format(
                          XFAIL_REASONS[rkey], measured),
            "ref": "tests/regression/README.md (the baseline table)",
        }
    entry = _case(name, os.path.join(OUT_REL, name + ".nml"), "rest",
                  T1_STEPS, T2_STEPS if prob.get("tier2") else 0, **kw)
    if not prob.get("tier2"):
        entry["tier2"] = {"skip": True, "reason":
                          "the tier-2 slice is one problem per geometry "
                          "class (flat / slope / rx0); the rest of the "
                          "matrix is tier-1 nightly. See "
                          "tests/regression/README.md."}
    return entry


def _refusal_reason(pid, prob, fam, is_cavity):
    fam_id = fam[0]
    if is_cavity and fam_id not in CAVITY_ACCEPTED:
        return ("REFUSED BY DESIGN, and this row asserts the refusal. "
                "`&ocean_cavity_dyn_nml enable=.true.` accepts "
                "vcoord_type='sigma', 'zstar' or 'z_fixed' only -- the two "
                "families that rescale the live column and so follow the ice "
                "base for free, plus the one that has been TAUGHT it "
                "(z_fixed reads vcoord%z_top). '{}' is refused; "
                "`validate_config` prints the per-family reason. If this row "
                "ever turns XPASS the envelope widened: either the family was "
                "validated under a draft (delete the row and add the real "
                "one) or the fence was dropped (a regression).".format(fam_id))
    return fam[5]


def _row_note(pid, prob, fam, strat, eos, a_peak, en_est, strat_note,
              bar_name, refused):
    bits = ["{} x {} ({}, {} EOS).".format(pid, fam[0], strat, eos),
            prob["doc"], fam[5], strat_note]
    if refused:
        bits.append("REFUSAL ROW: the assertion is that the run does NOT "
                    "start.")
    elif a_peak > 0.0:
        bits.append(
            "DERIVED truncation estimate for a terrain-following family on "
            "this geometry: a_peak = N^2*de^3/(6*dx*Hbar) = {:.3e} m/s^2 "
            "with de = {:.4g} m, which a rotating run balances "
            "geostrophically at U = a_peak/|f| = {:.3e} m/s, i.e. "
            "En ~ {:.3e} m2/s2. That is an ESTIMATE (one dominant face step, "
            "geostrophic balance); the measured cavity case tracks it to "
            "within a factor of ~2. The bar is {} and is independent of it."
            .format(a_peak, prob["de"], a_peak / abs(CORIOLIS_F), en_est,
                    bar_name))
    else:
        bits.append("No truncation term exists on this cell (de = 0 or "
                    "N^2 = 0), so the correct answer is machine zero and the "
                    "bar is {}.".format(bar_name))
    return " ".join(bits)


def _header(name, pid, prob, fam, strat, eos, refused, fam_note):
    return (
        "! GENERATED by tests/regression/vcoord_matrix.py -- DO NOT EDIT.\n"
        "! Matrix cell: problem={}  family={}  stratification={}  eos={}\n"
        "! Expectation: {}\n"
        "! Geometry: {}\n"
        "! Family:   {}\n"
        "! Reproduce by hand:  ./rdb tmp_local_artifacts/vcoord_matrix/{}.nml\n"
        "! Regenerate:         python3 tests/regression/stability.py --self-test\n"
        "!\n".format(
            pid, fam[0], strat, eos,
            "REFUSED at configure" if refused else "runs, and is gated",
            " ".join(prob["doc"].split()),
            " ".join(fam_note.split()), name))
