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

Two legs, one fixed configuration
=================================
Every cell runs the configuration v0.1.0 recommends over sloping topography
(exact FV PGF, linear-exact remap boundary cells + non-uniform weights, the
fail-loud remap precondition guard, closed z_fixed staircase faces).  The
INVISCID leg (`vcm_*`) is the hard probe; its terrain-following growth is
documented expected behaviour shared with MOM6.  The VISCOUS leg (`vcmv_*`)
carries MOM6's shipped seamount closure translated onto this grid (see LEGS)
and is the PASS/FAIL gate: a viscous cell inside its family's rx0 ENVELOPE
may not carry a marker.  The numbers -- markers and envelopes -- are pinned
from sweeps on both toolchains by `vcoord_matrix_pin.py` into
`vcoord_matrix_measured.py`; the reasons are policy and live here.

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

# --- the two LEGS: the same problem, with and without dissipation ----------
# INVISCID is the hard probe (nothing damps anything, so whatever grows is the
# discretisation's own); VISCOUS is the PASS/FAIL gate, carrying the closure
# MOM6 ships with its own seamount rest test (`ocean_only/seamount`,
# dev/gfdl d74a11f9c, `MOM_parameter_doc.all`):
#
#   LAPLACIAN = True, KH = 1000 m2/s, KH_VEL_SCALE = 0.003 m/s, BOUND_KH
#   BOTTOMDRAGLAW = LINEAR_DRAG = True, CDRAG = 0.002, DRAG_BG_VEL = 0.05 m/s,
#   HBBL = 10 m;  BIHARMONIC = False;  on a 5 km grid at dt = 900 s.
#
# THE TRANSLATION, and why each number is what it is:
#
#   * VISCOSITY.  Both codes' Laplacian coefficient is the nu of du/dt =
#     nu * del^2 u: MOM6's `diffu = (1/h) div(h Kh strain)` (MOM_hor_visc.F90,
#     tension dudx-dvdy + shear dvdx+dudy) and roundabout's scalar path
#     (`hvisc_compute_scalar_impl`, 5-point velocity Laplacian) both reduce to
#     Kh * del^2 u on a uniform grid with uniform h -- the cross terms of the
#     stress form cancel -- so KH maps onto `nu_h` one-for-one in DEFINITION.
#     The OPERATORS differ where h varies between neighbours: MOM6
#     thickness-weights and conserves momentum, the scalar path does neither.
#     roundabout's `stress_tensor=.true.` IS MOM6's operator (tension/shear
#     stress times h, divided by h_u + h_neglect, per-cell BOUND_KH clamp,
#     coast-masked), so the viscous leg selects it.  It is refused under
#     `zfixed_closed_faces`, so the z_fixed cells carry the scalar operator
#     (with the closed-face free-slip masks).  Measured, and the reason this
#     matters: the SCALAR operator at nu_h >= 40 m2/s goes explosively
#     unstable on `rx0_060 x sigma` under pred_corr at dt = 600 s, where the
#     stress-divergence operator at nu_h = 160 m2/s rests -- see the FINDING
#     in tests/regression/README.md.
#     The MAGNITUDE is translated, not copied: what a Laplacian does to a
#     grid-scale mode is damp it at nu * k_grid^2 ~ nu / dx^2, and the mode
#     this matrix exists for is 2-3 dx wide.  MOM6's shipped Kh/dx^2 =
#     1000 / 5000^2 = 4.0e-5 1/s; the same grid-scale damping on this 2 km
#     grid is nu_h = 4.0e-5 * 2000^2 = 160 m2/s.  Copying 1000 m2/s instead
#     would be 6.25x MOM6's grid-scale damping -- and MOM6 itself would not
#     run it here: its BOUND_KH ceiling 0.1 * dx^2 / dt is 667 m2/s at
#     (2 km, 600 s).  KH_VEL_SCALE * dx = 6 m2/s is below either, so it
#     does not bind (roundabout's `kh_vel_scale` only seeds `ah_bg` for the
#     flow-aware closures and would be inert on this path anyway).
#   * DRAG.  MOM6 LINEAR_DRAG is a bottom STRESS tau/rho_0 = CDRAG *
#     DRAG_BG_VEL * u_bbl = 1.0e-4 m/s * u_bbl, grid-independent, carried over
#     the bottom HBBL.  roundabout's distributed linear form applies
#     du_k/dt = -r * u_k * (h_in_bbl_k / h_k), whose column integral is
#     r * HBBL * u_bbl, so r = 1.0e-4 / HBBL = 1.0e-5 1/s with hbbl = 10 m
#     reproduces the stress exactly.
#   * NOT carried: MOM6's KV = 1e-4 m2/s vertical viscosity. Its damping of
#     the first baroclinic mode, KV * (pi/H)^2 ~ 1e-9 1/s, is four decades
#     under the growth rates this matrix measures; `use_closure=.false.`
#     keeps KD = 0 too, which the tracer-extrema gate needs.
MOM6_SEAMOUNT_DX = 5000.0
MOM6_KH = 1000.0
MOM6_KH_VEL_SCALE = 0.003
MOM6_CDRAG = 0.002
MOM6_DRAG_BG_VEL = 0.05
MOM6_HBBL = 10.0
VISC_NU_H = round(max(MOM6_KH * (DX / MOM6_SEAMOUNT_DX) ** 2,
                      MOM6_KH_VEL_SCALE * DX), 6)
VISC_BDRAG_R = MOM6_CDRAG * MOM6_DRAG_BG_VEL / MOM6_HBBL

LEGS = {
    "inviscid": {"prefix": "vcm", "nu_h": 0.0, "bdrag_form": "quadratic",
                 "bdrag_r": 0.0, "bdrag_hbbl": 0.0, "stress_tensor": False},
    "viscous": {"prefix": "vcmv", "nu_h": VISC_NU_H, "bdrag_form": "linear",
                "bdrag_r": VISC_BDRAG_R, "bdrag_hbbl": MOM6_HBBL,
                "stress_tensor": True},
}
# `stress_tensor` is refused under `&vcoord_nml zfixed_closed_faces`, so the
# z_fixed cells of the viscous leg carry the scalar velocity Laplacian.
SCALAR_LAPLACIAN_FAMILIES = ("z_fixed",)

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
    ("z_fixed", "z_fixed", "run", None, "   zfixed_closed_faces = .true.",
     "quasi-geopotential: fixed-z interfaces with bed-side layers vanishing "
     "to the inert filler. `z_fixed_h_ref` is taken from "
     "`&ocean_topo_nml max_depth`, so the nominal spacing is max_depth/nz. "
     "Runs with `zfixed_closed_faces` (the staircase faces closed, the "
     "barotropic mode open), which is how v0.1.0 recommends z_fixed be "
     "run and without which the staircase PGF residual is not bounded."),
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
        "tier2": False, "tier2_viscous": True,
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
        "tier2": (abs(_r - 0.6) < 1e-9), "tier2_viscous": False,
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
# THE BASELINE TABLE -- measured, then pinned (by a tool, never by hand)
# ---------------------------------------------------------------------------
# The NUMBERS live in `vcoord_matrix_measured.py`, generated by
# `vcoord_matrix_pin.py` from `stability.py --out` JSONs of a tier-1 and a
# tier-2 sweep on EACH toolchain (see that module's header for provenance).
# What lives HERE is policy: which documented reason a failing cell carries.
# Each record names the ASSERTIONS it covers (the union over toolchains and
# tiers -- the tolerance band), so a cell that is XFAIL on
# `energy:rest-growth-rate` still gates `conserve:*`, `finite` and the level
# bar, and a cell that fails a NEW assertion on any toolchain turns FAIL.
#
# The MOM6 evidence the inviscid reasons cite is recorded in
# tests/regression/README.md ("The MOM6 baseline"): MOM6 dev/gfdl d74a11f9c,
# ocean_only/seamount, sigma, at rest, f = 1e-4, every dissipation off.
MOM6_EVIDENCE = (
    "MOM6 SHOWS THE SAME MODE: its own seamount at rest in sigma "
    "coordinates, inviscid, f = 1e-4 s^-1, grows kinetic energy "
    "exponentially out of round-off with an e-folding of 0.81-0.85 days "
    "(R^2 > 0.9998 over five decades), independent of rx0 over 0.10-0.76, "
    "identically zero on a flat bed and absent at f = 0; nu_h = 10 m2/s does "
    "not stop it there, and MOM6's shipped closure (KH = 1000 m2/s + "
    "KH_VEL_SCALE + linear drag) does (MOM6 dev/gfdl d74a11f9c; "
    "tests/regression/README.md, 'The MOM6 baseline'). It is a property of "
    "the C-grid/ALE method on a sloping boundary, not a roundabout defect.")

XFAIL_REASONS = {
    "inviscid_mode":
        "DOCUMENTED EXPECTED BEHAVIOUR (inviscid leg). With no dissipation of "
        "any kind, a terrain-following column over a slope amplifies "
        "round-off in the pressure difference at the steepest face (1.1e-16 "
        "m/s^2, below one ulp of the hydrostatic pressure -- the FV PGF is "
        "exact everywhere else) into a baroclinic, grid-scale (2-3 dx), "
        "rotation-dependent, stratification-dependent, dt-INDEPENDENT mode "
        "(design/vcoord_ale_audit.md, 2026-09-21 forensics Q2). It grows "
        "exponentially for the whole run and, on the steep rungs, reaches "
        "the CFL wall, where continuity writes a negative layer and the "
        "remap guard stops the run. " + MOM6_EVIDENCE,
    "unstrat_control":
        "THE N^2 = 0 CONTROL (inviscid leg). With no stratification the "
        "baroclinic mode is absent; what survives on the steepest rung is a "
        "slower BAROTROPIC residual growth (measured: KE_bc/KE_bt 1.7e-05, "
        "forensics Q2b) that still reaches the CFL wall inside 30 days at "
        "rx0 = 0.8. Carried as a control, not as a coordinate verdict.",
    "z_fixed_leak":
        "z_fixed LEAKS SALT AND HEAT wherever layers vanish -- a step at the "
        "first regrid, 1e-7 relative against the 1e-11 bar, scaling with "
        "the number of filler layers, not with the energy. "
        "`remap_boundary_extrap` bought 71x of it and it is still four "
        "decades over (audit 2026-09-21, Result 4). The fix is "
        "`origin/fix/remap-vanished-layer-content` (vanished-layer contract "
        "I1: a sub-threshold filler may not carry content), not yet on main. "
        "MEASURED with it (viscous leg, gfortran, 30 d): the salt residual "
        "falls 4.4e-07 -> 1.7e-13 (slope), 1.2e-06 -> 8.3e-14 "
        "(seamount_gentle), 7.3e-07 -> 1.2e-14 (rx0_060), i.e. to round-off; "
        "what remains is `tracer:no-new-extrema` (1-3 mPSU of salinity "
        "overshoot from the regrid) and, under the sloping lid, the "
        "saturated staircase residual (1.15 cm/s). Re-measure when it "
        "lands. The flat geometries, which have no fillers, close at "
        "round-off.",
    "z_fixed_staircase":
        "z_fixed's STAIRCASE PGF RESIDUAL: a forced, bounded truncation that "
        "saturates (static in level, dt-invariant; audit forensics Q2e), "
        "not an instability -- but it sits over the rest bar.",
    "density_coord":
        "A DENSITY-SPACE COORDINATE ON A GEOMETRIC REST STATE (inviscid). "
        "`rho`/`hycom` hand the remap a column whose layers track "
        "isopycnals; on a slope the rest-state mode drives a layer negative "
        "within days and the guard stops the run (the 1.00008 column-total "
        "mismatch the density target builder then reports is downstream of "
        "that negative layer -- audit forensics Q1.4). `rho` is "
        "validation-grade by documentation.",
    "eulerian_z":
        "`eulerian_z` is a stretched sigma with the free surface dropped "
        "(and pins ssp_rk2, which pred_corr refuses): on a slope it carries "
        "the same rest-state mode, amplified by the ssp_rk2 two-stage "
        "average, and six decades above sigma's level.",
    "lagrangian":
        "`lagrangian` inherits the sigma-shaped initial thickness and then "
        "never regrids: over a step the layers deform until one collapses "
        "and the run goes non-finite. Expected of a pure isopycnal "
        "coordinate carrying a geometric rest state.",
    "finding_visc_pred_corr":
        "FINDING B (viscous leg) -- NOT expected, NOT hidden. With the "
        "MOM6-comparable Laplacian viscosity (nu_h = 160 m2/s, either "
        "operator) a terrain-following column over a step that REMAINS "
        "FINITE without viscosity goes explosively unstable under the "
        "default pred_corr at dt = 600 s: a round-off creep at the step "
        "face's top layer turns, after 8-15 days, into a domain-wide "
        "BAROTROPIC grid-scale (2 dx) mode e-folding in under an hour, and "
        "continuity drives a layer negative. First-failure localisation "
        "(tests/regression/README.md, 'FINDING B'). LOCALISED: not the "
        "viscosity -- a 1e-6 m/s barotropic seed blows the cell up in 2-3 "
        "days inviscid, at f = 0, at dt = 300 s and on ONE layer. The "
        "uhbt renormalisation in continuity has no root when a layer's "
        "upwind donor flips across the step's thickness jump, returns a "
        "wrong-sign transport, and the layer eta leaves the barotropic "
        "eta_end, pumping the (under pred_corr undamped) 2 dx barotropic "
        "mode. Fixed by &ocean_continuity_nml renorm_consistent_flux "
        "(MOM6 zonal_flux_adjust parity), not yet in this template: adding "
        "it moves pinned markers on both toolchains.",
    "finding_stress_density":
        "FINDING A (viscous leg) -- NOT expected, NOT hidden. The "
        "thickness-weighted stress-divergence viscosity "
        "(`stress_tensor=.true.`, MOM6's operator) drives a thin layer of a "
        "density-space column negative within the first 4-8 outer steps "
        "(mm to dm of negative h, every sloping geometry, both split "
        "schemes); the scalar velocity Laplacian and no viscosity at all "
        "are clean. The one MOM6 thin-layer safeguard the port lacks is "
        "`hrat_min = min(1, h_min/h)` scaling the BOUND_KH ceiling "
        "(MOM_hor_visc.F90) -- a hypothesis, not a measurement.",
    "z_fixed_viscous_abort":
        "The viscous z_fixed cell additionally ABORTS on the remap guard "
        "(the inviscid twin completes). z_fixed carries the scalar "
        "Laplacian (stress_tensor is refused under closed faces), so the "
        "pred_corr x Laplacian-viscosity instability of FINDING B is the "
        "suspect; it was not localised on z_fixed.",
    "finding_gpu_wright_density":
        "FINDING C (GPU only) -- NOT expected, NOT hidden. On nvfortran 26.5 "
        "/ V100 EVERY rho or hycom run with eos = 'wright' dies at the "
        "first regrid with CUDA_ERROR_ILLEGAL_ADDRESS inside "
        "`ocean_vcoord_compute_target_h_rho_impl` (the column do concurrent, "
        "rdb_ocean_vcoord.F90) -- on a FLAT bed at rest too; the linear EOS "
        "is clean on the GPU, and gfortran with -fcheck=all finds no bounds "
        "violation on the same run. The Wright coefficients are parameters, "
        "so the point EOS call is not the suspect; the kernel's `associate` "
        "over `this%` components around the do concurrent (the NVHPC "
        "mapping hazard CLAUDE.md records) is -- a hypothesis.",
    "outside_envelope":
        "OUTSIDE THE FAMILY'S DOCUMENTED rx0 ENVELOPE (viscous leg). The "
        "envelope is the largest geometry rx0 at which EVERY viscous cell "
        "of the family passes on every toolchain measured (ENVELOPES, "
        "docs/CAPABILITIES_AND_LIMITATIONS.md); past it the family is not "
        "claimed to rest.",
}

try:
    import vcoord_matrix_measured as _measured
    ENVELOPES = dict(_measured.ENVELOPES)
    MEASURED = _measured.MEASURED
    MEASURED_TWIN = _measured.MEASURED_TWIN
    MEASURED_PROVENANCE = (_measured.__doc__ or "").split("Provenance:", 1)[-1].strip()
except ImportError:          # bootstrapping a first measurement
    ENVELOPES, MEASURED_PROVENANCE = {}, "(none)"
    MEASURED = MEASURED_TWIN = {"inviscid": {}, "viscous": {}}

TERRAIN_FOLLOWING = ("sigma", "zstar", "zstar_sigma", "zstar_full")


def reason_for(leg, key, rec):
    """Which documented reasons a measured failing cell carries (policy).

    Returns a list of `XFAIL_REASONS` keys, primary first: a cell can carry
    two independent defects (a z_fixed leak AND an abort), and its marker
    must name both rather than let one excuse the other in prose.
    """
    parts = key.split("/")
    fam = parts[1]
    a = set(rec["assertions"])
    out = []
    if "unstrat" in parts:
        return ["unstrat_control"]
    if fam == "z_fixed":
        out.append("z_fixed_leak" if a & {"conserve:Salt", "conserve:Heat"}
                   else "z_fixed_staircase")
        if "completed" in a and leg == "viscous":
            out.append("z_fixed_viscous_abort")
        return out
    if leg == "inviscid":
        if fam in TERRAIN_FOLLOWING:
            out.append("inviscid_mode")
        elif fam in ("rho", "hycom"):
            out.append("density_coord")
        else:
            out.append(fam if fam in XFAIL_REASONS else "outside_envelope")
    elif fam in ("rho", "hycom"):
        out.append("finding_stress_density")
    elif fam in TERRAIN_FOLLOWING and "completed" in a:
        out.append("finding_visc_pred_corr")
    elif fam == "lagrangian":
        out.append("lagrangian")
    else:
        out.append("outside_envelope")
    if fam in ("rho", "hycom") and "wright" in parts:
        out.append("finding_gpu_wright_density")
    return out


def _marker(leg, key, rec):
    rkeys = reason_for(leg, key, rec)
    rkey = rkeys[0]
    nums = "; ".join("{}: {}".format(k, v)
                     for k, v in sorted(rec["measured"].items()))
    kf = {"assertions": list(rec["assertions"]),
          "tiers": list(rec["tiers"]),
          "reason": "{}  MEASURED ({}): {}.{}".format(
              "  ALSO: ".join(XFAIL_REASONS[k] for k in rkeys),
              MEASURED_PROVENANCE, nums,
              "  TOOLCHAIN-DEPENDENT: passes outright on at least one "
              "toolchain, so a pass reports PASS, not XPASS."
              if rec["toolchain_dependent"] else ""),
          "ref": "tests/regression/README.md (the baseline table)",
          "reason_key": rkey, "reason_keys": rkeys}
    if rec["toolchain_dependent"]:
        kf["toolchain_dependent"] = True
    return kf


def _marker_to_known_failure(leg, mkey, marker):
    return _marker(leg, mkey, marker)


MARKERS = MEASURED


def geometry_rx0(prob):
    """The worst wet-wet Beckmann-Haidvogel stiffness rx0 the geometry has.

    `rx0 = |H_a - H_b| / (H_a + H_b)` over every interior face.  Exact for
    the ladder (its whole point) and for the Gaussian seamounts (evaluated
    from the same formula `set_bathymetry_seamount` fills); `de / (2 Hbar)`
    for the constant-gradient slope and the lids, whose step per face is
    uniform.  This is what places EVERY problem -- not only the ladder
    rungs -- on the envelope axis.
    """
    if "rx0" in prob:
        return prob["rx0"]
    if prob["topo"] == "seamount":
        hmax = prob.get("max_depth", MAX_DEPTH)
        dep = hmax - prob["edge_depth"]
        lw = prob["slope_scale"]
        xc, yc = 0.5 * NX * DX, 0.5 * NY * DY

        def d(i, j):
            x = (i - 0.5) * DX
            y = (j - 0.5) * DY
            return hmax - dep * math.exp(-((x - xc) ** 2 + (y - yc) ** 2) / lw ** 2)
        worst = 0.0
        for j in range(1, NY + 1):
            for i in range(1, NX + 1):
                for (a, b) in ((i + 1, j), (i, j + 1)):
                    if a > NX or b > NY:
                        continue
                    h1, h2 = d(i, j), d(a, b)
                    worst = max(worst, abs(h1 - h2) / (h1 + h2))
        return worst
    if prob["de"] <= 0.0:
        return 0.0
    return prob["de"] / (2.0 * prob["hbar"])


LEG_NOTES = {
    "inviscid": "INVISCID leg (nu_h = 0, no drag): the hard probe. Nothing "
                "damps anything, so whatever grows is the discretisation's "
                "own; its sigma-family growth is documented EXPECTED "
                "behaviour shared with MOM6.",
    "viscous": "VISCOUS leg: MOM6's shipped seamount closure translated onto "
               "this grid (nu_h = {:g} m2/s, linear drag r = {:g} 1/s over "
               "hbbl = {:g} m). This leg is the PASS/FAIL gate.".format(
                   VISC_NU_H, VISC_BDRAG_R, MOM6_HBBL),
}

def in_envelope(fam_id, prob):
    """True/False when the family has a documented rx0 envelope; else None."""
    lim = ENVELOPES.get(fam_id)
    if lim is None:
        return None
    return geometry_rx0(prob) <= lim + 1e-12


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


def render(name, prob, fam, strat, eos, remap="ppm", leg="inviscid"):
    """Render one matrix cell's namelist text from the canonical template."""
    fam_id, vcoord, status, pin, extra, _note = fam
    lg = LEGS[leg]
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

    # Every cell runs the FIXED configuration's exact FV pressure gradient
    # (`form="fv_mom6"` + `reconstruct_for_pressure`, set in the template);
    # the sloping lid additionally injects the ice load as the stack's top
    # boundary condition.
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
        "P_TOP_IN_BC": p_top, "EOS": eos, "LEG": leg.upper(),
        "NU_H": lg["nu_h"], "BDRAG_FORM": lg["bdrag_form"],
        "BDRAG_R": lg["bdrag_r"], "BDRAG_HBBL": lg["bdrag_hbbl"],
        "STRESS_TENSOR": bool(lg["stress_tensor"]
                              and fam_id not in SCALAR_LAPLACIAN_FAMILIES),
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
    for leg in ("inviscid", "viscous"):
        rows += _build_leg(_case, bars, out, leg)
    return rows


def _build_leg(_case, bars, out, leg):
    """Every cell of ONE leg.

    The VISCOUS leg carries the runnable cells only: a REFUSAL is a property
    of the configuration, not of the dissipation, so asserting it twice says
    nothing new; and the N^2 = 0 control exists to isolate the inviscid
    truncation term, which is what the inviscid leg measures.
    """
    rows = []
    for pid in sorted(PROBLEMS):
        prob = PROBLEMS[pid]
        is_cavity = prob["class"] == "cavity"
        for fam in FAMILIES:
            fam_id, _vc, status, pin, _extra, fam_note = fam
            refused = (status == "refused") or (
                is_cavity and fam_id not in CAVITY_ACCEPTED)
            if refused and leg != "inviscid":
                continue
            rows.append(_one(_case, bars, out, pid, prob, fam, "linear",
                             "linear", is_cavity, status, pin, fam_note, leg))
        # The N^2 = 0 CONTROL, on the two geometries where it says the most:
        # the steepest ladder rung and the sloping lid.  It is a control, not
        # a coordinate test, so it runs on sigma alone -- the family whose
        # truncation term the control is isolating.
        if pid in ("rx0_080", "lid_slope") and leg == "inviscid":
            fam = [f for f in FAMILIES if f[0] == "sigma"][0]
            rows.append(_one(_case, bars, out, pid, prob, fam, "unstrat",
                             "linear", is_cavity, "run", None, fam[5], leg))
        # The nonlinear-EOS leg, on the clean slope geometry: Wright (1997)
        # is the production EOS and its thermobaricity is exactly what a
        # linear-EOS rest state cannot see.
        if pid == "slope":
            for fam in FAMILIES:
                if fam[0] in ("sigma", "z_fixed", "hycom"):
                    rows.append(_one(_case, bars, out, pid, prob, fam,
                                     "linear", "wright", is_cavity, fam[2],
                                     fam[3], fam[5], leg))
    return rows


def cell_key(pid, fam_id, strat="linear", eos="linear"):
    """The marker key of one cell within a leg: `problem/family[/strat][/eos]`."""
    return "/".join([pid, fam_id]
                    + ([strat] if strat != "linear" else [])
                    + ([eos] if eos != "linear" else []))


def cell_name(leg, pid, fam_id, strat="linear", eos="linear"):
    """The manifest / scratch-dir name of one cell."""
    parts = [LEGS[leg]["prefix"], pid, fam_id]
    if strat != "linear":
        parts.append(strat)
    if eos != "linear":
        parts.append(eos)
    return "_".join(parts)


def _one(_case, bars, out, pid, prob, fam, strat, eos, is_cavity, status,
         pin, fam_note, leg="inviscid"):
    """Build one matrix cell: write its namelist, return its manifest row."""
    fam_id = fam[0]
    name = cell_name(leg, pid, fam_id, strat, eos)

    refused = (status == "refused") or (is_cavity and fam_id not in CAVITY_ACCEPTED)
    text = render(name, prob, fam, strat, eos, leg=leg)
    path = os.path.join(out, name + ".nml")
    with open(path, "w") as fh:
        fh.write(_header(name, pid, prob, fam, strat, eos, refused, fam_note,
                         leg))
        fh.write(text)

    _dt_dz, _ds_dz, n2, strat_note = STRATIFICATIONS[strat]
    a_peak, en_est = problem_truncation_estimate(prob, n2)
    bar_name = prob["bar"] if n2 > 0.0 else "REST_1UM_S"
    kw = {
        "en_rest_max": getattr(bars, bar_name),
        "tags": ["vcoord_matrix", "vcoord_matrix_" + leg, "rest",
                 "vcoord_" + fam_id, prob["class"]],
        "note": _row_note(pid, prob, fam, strat, eos, a_peak, en_est,
                          strat_note, bar_name, refused, leg),
        "t1_timeout": 1800, "t2_timeout": 300,
        "t1_samples": 60, "t2_samples": 30,
        "matrix_gates": True,
        "matrix": {"problem": pid, "family": fam_id, "leg": leg,
                   "stratification": strat, "eos": eos,
                   "expect": "refused" if refused else "run",
                   "class": prob["class"], "rx0": prob.get("rx0"),
                   "de": prob["de"], "a_peak": a_peak,
                   "rx0_geometry": geometry_rx0(prob),
                   "en_estimate": en_est,
                   "envelope": in_envelope(fam_id, prob)},
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
    # MARKERS -- the numbers are measurements, not choices.
    mkey = cell_key(pid, fam_id, strat, eos)
    marker = MARKERS[leg].get(mkey)
    if not refused and marker:
        kw["known_failure"] = _marker_to_known_failure(leg, mkey, marker)
    if not refused:
        # The `__ssp_rk2` twin carries ITS OWN measured record (or none: it
        # was measured passing), read by stability_manifest._scheme_twin.
        twin = MEASURED_TWIN[leg].get(mkey)
        kw["matrix"]["twin_known_failure"] = (
            _marker(leg, mkey, twin) if twin else None)
    t2 = prob.get("tier2_" + leg, prob.get("tier2"))
    entry = _case(name, os.path.join(OUT_REL, name + ".nml"), "rest",
                  T1_STEPS, T2_STEPS if t2 else 0, **kw)
    if not t2:
        entry["tier2"] = {"skip": True, "reason":
                          "the tier-2 slice is a few problems per leg (see "
                          "TIER2_SLICE in vcoord_matrix.py); the rest of the "
                          "matrix is tier-1, local. See "
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
              bar_name, refused, leg="inviscid"):
    bits = ["{} x {} ({}, {} EOS), {} leg.".format(pid, fam[0], strat, eos,
                                                  leg.upper()),
            prob["doc"], fam[5], strat_note, LEG_NOTES[leg]]
    env = in_envelope(fam[0], prob)
    if not refused and env is not None:
        bits.append("This cell is {} the {} family's documented rx0 envelope "
                    "(ENVELOPES in vcoord_matrix.py; "
                    "docs/CAPABILITIES_AND_LIMITATIONS.md).".format(
                        "INSIDE" if env else "OUTSIDE", fam[0]))
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


def _header(name, pid, prob, fam, strat, eos, refused, fam_note,
            leg="inviscid"):
    return (
        "! GENERATED by tests/regression/vcoord_matrix.py -- DO NOT EDIT.\n"
        "! Matrix cell: leg={}  problem={}  family={}  stratification={}  "
        "eos={}\n"
        "! Expectation: {}\n"
        "! Geometry: {}\n"
        "! Family:   {}\n"
        "! Reproduce by hand:  ./rdb tmp_local_artifacts/vcoord_matrix/{}.nml\n"
        "! Regenerate:         python3 tests/regression/stability.py --self-test\n"
        "!\n".format(
            leg, pid, fam[0], strat, eos,
            "REFUSED at configure" if refused else "runs, and is gated",
            " ".join(prob["doc"].split()),
            " ".join(fam_note.split()), name))
