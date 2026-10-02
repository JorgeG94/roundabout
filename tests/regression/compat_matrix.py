#!/usr/bin/env python3
"""Pairwise compatibility matrix: does every configuration the model accepts
also RUN, on a small realistic domain, with every pair of features met?

Why this exists
---------------
Every composition bug found in the autumn of 2026 was two features meeting for
the first time: GM x z_fixed fillers (NaN in 3 steps), closed faces x
GM/Redi/MLE, sponge x periodic seam, carried tendency x restart registry...
None needed exotic physics; each needed one specific PAIR of features that no
hand-written test held.  Per-feature tests cover features, not pairs.  This
suite enumerates pairs (design: python_prototypes/design/compat_matrix_plan.md).

What it does
------------
1. **The domain** (`merged_namelist`).  ONE synthetic family, sized for
   seconds: 24 x 16 x 10, dx = 20 km, f-plane + beta, a staircase shelf (with
   an along-shelf offset so steps face both ways) + an island (so the z-like
   coordinates carry closed faces and filler layers), stratified T/S with a
   front whose position tilts with depth (so slopes / GM / Redi are
   non-trivial), wind + a surface heat flux, 24 steps of 900 s.  The
   bathymetry and the z-level T/S are written as classic NetCDF by the
   stdlib writer below (`write_netcdf_classic`) -- nothing to install.
2. **The covering array** (`ipog`).  A seeded, deterministic t = 2 IPOG
   generator over the axes in `AXES`: every PAIR of axis values appears in at
   least one cell.  The exclusion rules are NOT re-encoded here (they would
   drift from the model): forbidden tuples are DISCOVERED by running
   `rdb --validate-only` on candidate cells and classifying each refusal
   against `compat_expect.py`.  An expected refusal -- and a deterministic
   runtime gap -- forbids the tuple that caused it, the array is regenerated
   around it to a fixed point (`build_matrix`), and the cell that produced the
   tuple stays in the report as its witness.
3. **The checks** (`evaluate_cell`), in order, stop at the first failure:
     1  REFUSED  -- expected (PHYSICAL / KNOWN_GAP) or FAIL;
     2  CRASH / NONFINITE -- non-zero exit, a nan-catch, the vanished-content
        tripwire (`check_vanished_content`, ON in every cell), the remap
        precondition guard, a non-finite console scalar; a crashing cell is
        re-run with `&ocean_debug_nml chksum` over the last steps to name the
        phase that MINTED the first non-finite;
     3  BUDGET -- the model's own closed mass/salt/heat residuals.
   Checks 4-6 (MPI decomposition, restart, GPU cross-backend band) slot into
   the same record (`CHECKS`); this phase runs 1-3 on CPU.
4. **The verdict**.  A KNOWN_GAP row that does not fail where it says it
   will is an XPASS, and an XPASS FAILS the suite until the row is deleted --
   so the gap list can only shrink.

Usage
-----
    # the cell list (no model)
    python3 tests/regression/compat_matrix.py list
    # the whole pairwise set, CPU (gfortran toolchain loaded in this shell)
    python3 tests/regression/compat_matrix.py run --build-dir build_gfortran \\
            --jobs 4 --out tmp_local_artifacts/compat/last.json
    # every vertical coordinate x split x edge variant x grid, base closures
    python3 tests/regression/compat_matrix.py domain --build-dir build_gfortran
    # one cell by hand: emit its namelist + inputs and run rdb on it
    python3 tests/regression/compat_matrix.py emit --cell c017 --out DIR
    # who tests the test (no model; also in ctest)
    python3 tests/regression/compat_matrix.py self-test
    # the --validate-only contract (needs the binary; in ctest)
    python3 tests/regression/compat_matrix.py validate-smoke --binary build_gfortran/rdb

Stdlib only -- never `pip install` anything for this.
"""

import argparse
import hashlib
import itertools
import json
import math
import os
import random
import re
import shutil
import struct
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, "tools"))

import compat_expect  # noqa: E402

DEFAULT_SCRATCH = os.path.join(REPO_ROOT, "tmp_local_artifacts", "compat_matrix")
SEED = 20261002

# ===========================================================================
# 1. The domain
# ===========================================================================
NX, NY, NZ = 24, 16, 10
DX = 20.0e3
DT = 900.0
N_STEPS = 24
MAX_DEPTH = 2000.0
F0 = 1.0e-4
BETA = 2.0e-11

# The staircase: depth by (shelf-relative) row, south to north.  Each step is
# a distinct z_fixed partial column at nz = 10 / dz = 200 m, and rx0 stays
# <= 0.5 so the terrain-following coordinates are stressed but not doomed
# (the configure audit WARNS above rx0 = 0.2; it never refuses).
STAIRCASE = (100.0, 100.0, 250.0, 250.0, 500.0, 800.0, 1200.0, 1600.0)
SHELF_OFFSET_I = range(8, 16)   # the shelf is one row wider here: x-facing steps
ISLAND = (range(10, 13), range(9, 11))   # (i, j), 0-based, land

# The tilted front: T, S on z-levels (positive-down metres).
Z_SRC = (0.0, 10.0, 25.0, 50.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 1500.0,
         2000.0, 2500.0)
Y_FRONT = 0.5 * NY * DX      # front centre at the surface (m)
FRONT_TILT = 40.0            # metres of northward displacement per metre depth
FRONT_WIDTH = 50.0e3


def bathymetry():
    """Depth (m, positive down, 0 = land) as rows [j][i], j = 0 the south."""
    b = []
    for j in range(NY):
        row = []
        for i in range(NX):
            js = j - (1 if i in SHELF_OFFSET_I else 0)
            d = STAIRCASE[js] if 0 <= js < len(STAIRCASE) else (
                STAIRCASE[0] if js < 0 else MAX_DEPTH)
            if i in ISLAND[0] and j in ISLAND[1]:
                d = 0.0
            row.append(d)
        b.append(row)
    return b


def ts_profile(j, z):
    """(T degC, S PSU) at row j and depth z (m, positive down)."""
    y = (j + 0.5) * DX
    xi = (y - Y_FRONT - FRONT_TILT * z) / FRONT_WIDTH
    t = 3.0 + 13.0 * math.exp(-z / 500.0) + 1.5 * math.tanh(xi) * math.exp(-z / 1000.0)
    s = 34.6 + 0.4 * (1.0 - math.exp(-z / 400.0)) - 0.1 * math.tanh(xi) * math.exp(-z / 800.0)
    return t, s


# ---------------------------------------------------------------------------
# A classic-format (CDF-1) NetCDF writer: enough for a few double variables.
# The format is a fixed big-endian header followed by the data; see the
# NetCDF "Classic Format Specification".  Stdlib only, by repo rule.
# ---------------------------------------------------------------------------
_NC_DIMENSION, _NC_VARIABLE, _NC_DOUBLE = 0x0A, 0x0B, 6


def _nc_name(name):
    raw = name.encode("ascii")
    return struct.pack(">i", len(raw)) + raw + b"\0" * ((4 - len(raw) % 4) % 4)


def write_netcdf_classic(path, dims, variables):
    """Write `variables` [(name, (dimname, ...), flat C-order values)] over
    `dims` [(name, length)] as a CDF-1 file of doubles (no attributes)."""
    dim_index = {name: k for k, (name, _) in enumerate(dims)}
    head = b"CDF\x01" + struct.pack(">i", 0)
    head += struct.pack(">ii", _NC_DIMENSION, len(dims))
    for name, length in dims:
        head += _nc_name(name) + struct.pack(">i", length)
    head += struct.pack(">ii", 0, 0)                      # no global attributes
    var_hdr_len = 8
    for name, vdims, _ in variables:
        var_hdr_len += len(_nc_name(name)) + 4 + 4 * len(vdims) + 8 + 4 + 4 + 4
    offset = len(head) + var_hdr_len
    var_hdr, blobs = b"", []
    for name, vdims, values in variables:
        n = 1
        for d in vdims:
            n *= dims[dim_index[d]][1]
        if len(values) != n:
            raise ValueError("variable {} has {} values, dims want {}".format(name, len(values), n))
        blob = struct.pack(">{}d".format(n), *values)
        var_hdr += _nc_name(name) + struct.pack(">i", len(vdims))
        var_hdr += b"".join(struct.pack(">i", dim_index[d]) for d in vdims)
        var_hdr += struct.pack(">ii", 0, 0)                # no variable attributes
        var_hdr += struct.pack(">iii", _NC_DOUBLE, len(blob), offset)
        offset += len(blob)
        blobs.append(blob)
    with open(path, "wb") as fh:
        fh.write(head + struct.pack(">ii", _NC_VARIABLE, len(variables)) + var_hdr)
        for blob in blobs:
            fh.write(blob)


BATHY_FILE = "compat_bathy.nc"
TS_FILE = "compat_ts.nc"


def write_domain_inputs(directory):
    """Write the bathymetry and the z-level T/S the namelists point at."""
    b = bathymetry()
    write_netcdf_classic(os.path.join(directory, BATHY_FILE), [("y", NY), ("x", NX)],
                         [("b", ("y", "x"), [b[j][i] for j in range(NY) for i in range(NX)])])
    temp, salt = [], []
    for z in Z_SRC:
        for j in range(NY):
            t, s = ts_profile(j, z)
            temp.extend([t] * NX)
            salt.extend([s] * NX)
    write_netcdf_classic(os.path.join(directory, TS_FILE),
                         [("z", len(Z_SRC)), ("y", NY), ("x", NX)],
                         [("z_src", ("z",), list(Z_SRC)),
                          ("temp", ("z", "y", "x"), temp),
                          ("salt", ("z", "y", "x"), salt)])


# ===========================================================================
# 2. The base namelist and the axes
# ===========================================================================
# Base: what every cell carries before the axis overlays.  The vertical-
# coordinate safety configuration v0.1.0 recommends is ON everywhere (the
# remap guards and the vanished-content tripwire), so a cell that breaks an
# invariant fails loud instead of drifting.
BASE = {
    "sim_nml": {"sim_type": "ocean"},
    "grid_nml": {"nx": NX, "ny": NY, "dx": DX, "dy": DX, "nghost": 4},
    "time_nml": {"t_end": N_STEPS * DT, "dt_fixed": DT, "time_unit": "s"},
    "physics_nml": {"coriolis_f": F0, "wind_stress_x": 0.1, "wind_stress_y": 0.0},
    "nonhydrostatic_nml": {"nz_layers": NZ},
    "vcoord_nml": {
        "remap_method": "ppm", "zstar_h_min": 1.0e-4,
        "remap_boundary_extrap": True, "remap_nonuniform_weights": True,
        "remap_check_preconditions": True, "check_vanished_content": True,
        # rho / hycom targets: potential density at the surface, spanning
        # the IC's 1025.5 .. 1028.4 kg/m3 under every EOS on the axis.
        "rho_ref_pressure": 0.0, "rho_target_light": 1025.2, "rho_target_dense": 1028.6,
    },
    "ocean_topo_nml": {"topo_config": "file", "max_depth": MAX_DEPTH,
                       "wind_config": "constant", "coriolis_beta": BETA,
                       "coriolis_y_ref": 0.5 * NY * DX},
    "output_nml": {"output_dir": ".", "output_to_file": False,
                   "bathymetry_file": BATHY_FILE},
    "ocean_zinit_nml": {"enable": True, "source": "file", "file": TS_FILE},
    # Ghost columns are not covered by the z-level file: give them the same
    # stable profile, linear in layer index.
    "tracer_nml": {"initial_temperature": 9.0, "initial_salinity": 34.8,
                   "T_init_surface": 16.0, "T_init_bottom": 3.0,
                   "S_init_surface": 34.6, "S_init_bottom": 35.0},
    "ocean_ic_nml": {"alpha_T": 0.2, "beta_S": 0.77, "T_ref": 10.0, "S_ref": 35.0,
                     "rho_0": 1027.0},
    "ocean_bdrag_nml": {"form": "quadratic", "cd": 3.0e-3, "hbbl": 10.0, "bg_vel": 0.05},
    "ocean_bt_nml": {"auto_n_inner": True},
    "ocean_diag_nml": {"enabled": False},
    "logging_nml": {"log_level": "info", "status_interval": DT},
}

# Reusable closure blocks: a prerequisite two axis values share is spelled
# ONCE, so their overlays merge without a conflict.
_SLOPES = {"ocean_slopes_nml": {"enable": True}}
_GM = {"ocean_gm_nml": {"enable": True, "khth": 500.0}}
_REDI = {"ocean_redi_nml": {"enable": True, "khtr": 500.0}}
_MEKE = {"ocean_meke_nml": {"enable": True, "gmcoeff": 0.15, "khcoeff": 1.0,
                            "damping": 1.0e-6}}
_VARMIX = {"ocean_wavespeed_nml": {"enable": True},
           "ocean_varmix_nml": {"enable": True, "use_visbeck": True,
                                "khth_slope_cff": 0.1, "khtr_slope_cff": 0.1,
                                "visbeck_l_scale": 3.0e4, "khth_max": 2000.0,
                                "khtr_max": 2000.0}}
_SMAG = {"ocean_hvisc_nml": {"nu_h": 200.0, "lateral_closure": "smagorinsky",
                             "c_smag": 0.15, "smag_ah": True}}


def _merge(*blocks):
    out = {}
    for blk in blocks:
        for grp, kv in blk.items():
            out.setdefault(grp, {}).update(kv)
    return out


def _vc(vtype, **extra):
    return {"vcoord_nml": dict({"vcoord_type": vtype}, **extra)}


# The axes.  ORDER IS PART OF THE SEED: append new axes at the end and new
# values at the end of an axis, or every cell id moves.  Each value is
# (name, overlay).  `compat_expect.FEATURES` reads the merged namelist, never
# these names, so a row keeps matching when a value is renamed.
AXES = [
    ("vcoord", [
        ("sigma", _vc("sigma")),
        ("zstar", _vc("zstar")),
        ("zstar_full", _vc("zstar_full")),
        ("zstar_sigma", _vc("zstar_sigma")),
        ("z_fixed_cf", _vc("z_fixed", zfixed_closed_faces=True)),
        ("z_fixed_open", _vc("z_fixed", zfixed_closed_faces=False)),
        ("hycom", _vc("hycom")),
        ("rho", _vc("rho")),
        ("eulerian_z", _vc("eulerian_z")),
        ("lagrangian", _vc("lagrangian")),
        ("zsigma", _vc("zsigma")),
    ]),
    ("split", [
        ("pred_corr", {"ocean_bt_nml": {"split_scheme": "pred_corr"}}),
        ("ssp_rk2", {"ocean_bt_nml": {"split_scheme": "ssp_rk2"}}),
    ]),
    ("vmix_bl", [
        ("kpp", {"ocean_vmix_nml": {"use_closure": True, "use_kpp": True}}),
        ("epbl", {"ocean_vmix_nml": {"use_closure": True, "use_kpp": False},
                  "ocean_epbl_nml": {"enable": True, "mstar_scheme": "om4",
                                     "mld_use_prev_guess": True}}),
        ("pp81", {"ocean_vmix_nml": {"use_closure": True, "use_kpp": False}}),
    ]),
    ("vmix_extra", [
        ("none", {}),
        ("kappa_shear", {"ocean_kappa_shear_nml": {"enable": True}}),
        ("kappa_shear_vertex", {"ocean_kappa_shear_nml": {"enable": True, "at_vertex": True}}),
        ("tidal", {"ocean_tidal_mixing_nml": {"enable": True, "e_uniform": 1.0e-3}}),
        ("conv", {"ocean_conv_nml": {"enable": True, "kd_conv": 0.1}}),
        ("ddiff", {"ocean_ddiff_nml": {"enable": True}}),
    ]),
    ("vmix_bg", [
        ("scalar", {}),
        ("bryan_lewis", {"ocean_vmix_nml": {"bkgnd_profile": True}}),
        ("henyey", {"ocean_vmix_nml": {"bkgnd_henyey": True}}),
    ]),
    ("lateral", [
        ("const_nu_h", {"ocean_hvisc_nml": {"nu_h": 1000.0}}),
        ("smagorinsky", _SMAG),
        ("leith", {"ocean_hvisc_nml": {"nu_h": 200.0, "lateral_closure": "leith"}}),
        ("leith_biharm", {"ocean_hvisc_nml": {"nu_h": 200.0, "lateral_closure": "leith_biharm",
                                              "c_leith_bi": 1.0}}),
        ("nu_4", {"ocean_hvisc_nml": {"nu_h": 200.0, "nu_4": 1.0e11}}),
        ("stress_tensor", {"ocean_hvisc_nml": {"nu_h": 1000.0, "stress_tensor": True}}),
        ("kh_aniso", {"ocean_hvisc_nml": {"nu_h": 1000.0, "stress_tensor": True,
                                          "kh_aniso": 500.0}}),
        # MEKE backscatter needs GM (MEKE sources from it) and a non-zero
        # biharmonic backstop; both ride with the value.
        ("meke_backscatter", _merge(_SMAG, _SLOPES, _GM, _MEKE, {
            "ocean_meke_nml": {"backscatter": True, "backscatter_visc_coeff_ku": 1.0,
                               "alpha_deform": 1.0, "alpha_grid": 1.0}})),
    ]),
    ("eddy", [
        ("none", {}),
        ("gm", _merge(_SLOPES, _GM)),
        ("gm_meke", _merge(_SLOPES, _GM, _MEKE, _VARMIX)),
        ("redi", _merge(_SLOPES, _REDI)),
        ("gm_redi_meke", _merge(_SLOPES, _GM, _REDI, _MEKE, _VARMIX)),
        ("gm_varmix_resscaled", _merge(_SLOPES, _GM, _VARMIX, {
            "ocean_varmix_nml": {"resoln_scaled_khth": True, "resoln_scaled_khtr": True}})),
        ("mle", {"ocean_foxkemper_nml": {"enable": True, "use_mom_mixrate": True,
                                         "mld_decay_time": 86400.0}}),
    ]),
    ("tracers", [
        ("ts", {}),
        ("ideal_age", {"ocean_tracers_nml": {"enable_ideal_age": True}}),
        ("pseudo_salt", {"ocean_tracers_nml": {"enable_pseudo_salt": True}}),
    ]),
    ("pgf", [
        ("mont", {"ocean_pgf_nml": {"form": "mont"}}),
        ("fv_mom6", {"ocean_pgf_nml": {"form": "fv_mom6"}}),
        ("fv_mom6_plm", {"ocean_pgf_nml": {"form": "fv_mom6", "reconstruct_for_pressure": True,
                                           "recon_scheme": 1}}),
        ("fv_mom6_ppm", {"ocean_pgf_nml": {"form": "fv_mom6", "reconstruct_for_pressure": True,
                                           "recon_scheme": 2}}),
    ]),
    ("eos", [
        ("wright", {"ocean_eos_nml": {"eos": "wright"}}),
        ("roquet", {"ocean_eos_nml": {"eos": "roquet_spv"}}),
        # (TEOS-10 is refused on every backend today -- "not yet
        # device-callable" -- and is not on the plan's v1 axis.)
        ("linear", {"ocean_eos_nml": {"eos": "linear"}}),
    ]),
    ("coriolis", [
        ("sadourny", {"ocean_coriolis_nml": {"form": "sadourny"}}),
        ("sadourny_energy", {"ocean_coriolis_nml": {"form": "sadourny_energy"}}),
        ("sadourny_hk", {"ocean_coriolis_nml": {"form": "sadourny_hk"}}),
    ]),
    ("pv_adv", [
        ("centered", {"ocean_coriolis_nml": {"pv_adv_scheme": "centered"}}),
        ("weno3", {"ocean_coriolis_nml": {"pv_adv_scheme": "weno3"}}),
        ("weno5", {"ocean_coriolis_nml": {"pv_adv_scheme": "weno5"}}),
        ("weno7", {"ocean_coriolis_nml": {"pv_adv_scheme": "weno7"}}),
    ]),
    ("bt", [
        ("default", {}),
        ("correction_bc_pgf", {"ocean_bt_nml": {"correction_bc_pgf": True}}),
        ("substep_drag", {"ocean_bt_nml": {"substep_drag": True}}),
        ("wave_drag", {"ocean_bt_nml": {"wave_drag": True, "wave_drag_r_uniform": 1.0e-3}}),
        # The visc_rem family needs the implicit drag fold, which refuses an
        # HBBL-distributed drag: the value brings bed-only drag with it.
        ("visc_rem", {"ocean_bt_nml": {"correction_h_weighted": True,
                                       "correction_visc_rem": True},
                      "ocean_vdiff_nml": {"implicit_drag": True},
                      "ocean_bdrag_nml": {"hbbl": 0.0}}),
    ]),
    ("geometry", [
        # Walls all round; the staircase shelf on the south; the island.
        ("closed", {}),
        # Re-entrant in x, a sponge band against the north wall.
        ("channel", {"ocean_bc_nml": {"west": "periodic", "east": "periodic",
                                      "north": "sponge", "sponge_width": 3,
                                      "sponge_strength": 1.0e-4},
                     "ocean_sponge_nml": {"enable": True, "damp_source": "band",
                                          "target_source": "ic"}}),
        # A Flather open edge in the east.
        ("obc", {"ocean_bc_nml": {"east": "open"}}),
        # SINGLE-RANK rows (the MPI legs of a later phase skip them).
        # The analytic tripolar generator: a 15 x 1 degree ring from 59 N,
        # re-entrant in x, a bipolar cap above 70 N closed by the north fold
        # (px = 1) -- the geometry `test_ocean_tripolar_fold_mpi` proves.
        # (A coarse cap -- 7-degree rows from 40 S, join at 60 N -- loses a
        # column to a negative thickness at step 1 on every coordinate; the
        # cells next to the cap poles degenerate.  Not a matrix question.)
        # It IS a grid, so it excludes the `grid` axis's own generators.
        ("tripolar", {"grid_nml": {"dx": 15.0, "dy": 1.0},
                      "ocean_grid_nml": {"grid_config": "tripolar", "lon_west": 0.0,
                                         "lat_south": 59.0, "phi_join": 70.0,
                                         "lon_pole": 0.0, "coriolis_scheme": "planetary"},
                      "ocean_bc_nml": {"west": "periodic", "east": "periodic",
                                       "north": "tripolar_fold"}}),
        # A flat 200 m ice-shelf draft over the deep northern third
        # (open water south of y = 220 km), walls all round.
        ("cavity", {"ocean_cavity_dyn_nml": {"enable": True, "draft_config": "linear",
                                             "draft_depth": 200.0, "draft_slope": 0.0,
                                             "draft_x0": -1.0e5, "draft_y0": 11.0 * DX,
                                             # z_fixed's 2 x h_nominal rule (ISOMIP+):
                                             # nothing under this 1800 m cavity grounds.
                                             "h_min_cavity": 2.0 * MAX_DEPTH / NZ},
                    "ocean_pgf_nml": {"p_top_in_bc": True}}),
    ]),
    ("grid", [
        ("cartesian", {}),
        # ~20 km cells on a 40-43.5 N sector; f from the latitude.
        ("spherical", {"grid_nml": {"dx": 0.25, "dy": 0.18},
                       "ocean_grid_nml": {"grid_config": "spherical", "lon_west": 0.0,
                                          "lat_south": 40.0, "coriolis_scheme": "planetary"}}),
    ]),
    ("forcing", [
        # Wind + cooling: a convecting, wind-mixed surface layer.
        ("cool", {"ocean_thermo_nml": {"enable_thermodynamics": True, "q_heat": -60.0}}),
        # Wind + heating with penetrating shortwave.
        ("warm_sw", {"ocean_thermo_nml": {"enable_thermodynamics": True, "q_heat": 60.0,
                                          "sw_pen_frac": 0.4}}),
    ]),
]
AXIS_NAMES = [a for a, _ in AXES]
AXIS_INDEX = {a: k for k, a in enumerate(AXIS_NAMES)}
VALUE_NAMES = {a: [v for v, _ in vals] for a, vals in AXES}
OVERLAY = {(a, v): ov for a, vals in AXES for v, ov in vals}

# The base closures every vertical coordinate is first proved on (`domain`).
BASE_CELL = {"vcoord": "z_fixed_cf", "split": "pred_corr", "vmix_bl": "kpp",
             "vmix_extra": "none", "vmix_bg": "scalar", "lateral": "smagorinsky",
             "eddy": "none", "tracers": "ts", "pgf": "fv_mom6", "eos": "wright",
             "coriolis": "sadourny_energy", "pv_adv": "centered", "bt": "default",
             "geometry": "closed", "grid": "cartesian", "forcing": "cool"}


class BuilderConflict(Exception):
    """Two axis values set the same knob to different values: a table bug."""


def merged_namelist(cell):
    """The full {group: {key: value}} of a cell (base + every overlay)."""
    nml = {g: dict(kv) for g, kv in BASE.items()}
    owner = {}
    for axis in AXIS_NAMES:
        ov = OVERLAY[(axis, cell[axis])]
        for grp, kv in ov.items():
            tgt = nml.setdefault(grp, {})
            for k, v in kv.items():
                prev = owner.get((grp, k))
                if prev is not None and tgt.get(k) != v:
                    raise BuilderConflict("{}%{}: {}={} vs {}={}".format(
                        grp, k, prev, tgt.get(k), axis, v))
                tgt[k] = v
                owner[(grp, k)] = axis
    return nml


def _fortran_value(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, str):
        return '"{}"'.format(v)
    if isinstance(v, (list, tuple)):
        return ", ".join(_fortran_value(x) for x in v)
    raise TypeError("cannot render {!r}".format(v))


def render_namelist(nml, header=""):
    out = []
    for line in header.splitlines():
        out.append("! " + line if line else "!")
    for grp in sorted(nml):
        out.append("&" + grp)
        for k in sorted(nml[grp]):
            out.append("   {} = {}".format(k, _fortran_value(nml[grp][k])))
        out.append("/")
        out.append("")
    return "\n".join(out)


def cell_key(cell):
    return ";".join("{}={}".format(a, cell[a]) for a in AXIS_NAMES)


def cell_hash(cell):
    return hashlib.sha1(cell_key(cell).encode()).hexdigest()[:10]


def cell_header(cell, cid, extra=""):
    lines = ["GENERATED by tests/regression/compat_matrix.py -- do not edit a copy.",
             "cell {} ({})".format(cid, cell_hash(cell))]
    lines += ["  {:<11s} {}".format(a, cell[a]) for a in AXIS_NAMES]
    if extra:
        lines += extra.splitlines()
    return "\n".join(lines)


# ===========================================================================
# 3. The covering array: seeded IPOG with forbidden tuples
# ===========================================================================
def _violates(test, forbidden):
    """`test` (list of value indices or None) contains a forbidden tuple."""
    for tup in forbidden:
        for a, v in tup:
            if test[a] != v:
                break
        else:
            return True
    return False


def ipog(sizes, forbidden=(), seed=SEED):
    """A t = 2 covering array over axes of `sizes` values (In-Parameter-Order-
    General, Lei et al. 2007).  `forbidden` is a collection of tuples of
    (axis, value) pairs no row may contain in full.  Deterministic for a
    given (sizes, forbidden, seed).  Returns (rows, uncoverable) -- rows are
    lists of value indices; `uncoverable` lists the pairs no valid row holds.
    """
    forbidden = [tuple(sorted(t)) for t in sorted(set(tuple(sorted(t)) for t in forbidden))]
    n = len(sizes)
    rng = random.Random(seed)
    pref = []                       # a seeded, fixed value-preference order per axis
    for s in sizes:
        order = list(range(s))
        rng.shuffle(order)
        pref.append(order)
    order = sorted(range(n), key=lambda a: (-sizes[a], a))

    def ok(t):
        return not _violates(t, forbidden)

    def pair_ok(a, va, b, vb):
        t = [None] * n
        t[a], t[b] = va, vb
        return ok(t)

    a0, a1 = order[0], order[1]
    rows = []
    for va in pref[a0]:
        for vb in pref[a1]:
            t = [None] * n
            t[a0], t[a1] = va, vb
            if ok(t):
                rows.append(t)
    done = [a0, a1]
    uncoverable = []
    for k in order[2:]:
        need = set()
        for a in done:
            for va in range(sizes[a]):
                for vk in range(sizes[k]):
                    if pair_ok(a, va, k, vk):
                        need.add((a, va, vk))
        # Horizontal growth: give each row the value of k covering most.
        for t in rows:
            best, best_gain = None, -1
            for vk in pref[k]:
                t[k] = vk
                if not ok(t):
                    continue
                gain = sum(1 for a in done if t[a] is not None and (a, t[a], vk) in need)
                if gain > best_gain:
                    best, best_gain = vk, gain
            t[k] = best
            if best is not None:
                for a in done:
                    if t[a] is not None:
                        need.discard((a, t[a], best))
        # Vertical growth: fit each still-uncovered pair into a row with
        # don't-cares, else open a new row.
        for (a, va, vk) in sorted(need, key=lambda p: (p[0], pref[p[0]].index(p[1]),
                                                       pref[k].index(p[2]))):
            placed = False
            for t in rows:
                if t[a] not in (None, va) or t[k] not in (None, vk):
                    continue
                old = (t[a], t[k])
                t[a], t[k] = va, vk
                if ok(t):
                    placed = True
                    break
                t[a], t[k] = old
            if not placed:
                t = [None] * n
                t[a], t[k] = va, vk
                rows.append(t)
        done.append(k)

    # Fill the don't-cares with valid values (backtracking, seeded order).
    def fill(t, idx):
        while idx < n and t[idx] is not None:
            idx += 1
        if idx == n:
            return True
        for v in pref[idx]:
            t[idx] = v
            if ok(t) and fill(t, idx + 1):
                return True
        t[idx] = None
        return False

    full, seen = [], set()
    for t in rows:
        t = list(t)
        if not fill(t, 0):
            continue
        key = tuple(t)
        if key not in seen:
            seen.add(key)
            full.append(t)
    # A row whose don't-cares could not be completed (a higher-arity
    # forbidden tuple closed every option) dropped its pairs: re-seat each
    # lost pair in a fresh row of its own.
    covered = {(a, t[a], b, t[b]) for t in full for a, b in itertools.combinations(range(n), 2)}
    for a, b in itertools.combinations(range(n), 2):
        for va in pref[a]:
            for vb in pref[b]:
                if (a, va, b, vb) in covered or not pair_ok(a, va, b, vb):
                    continue
                t = [None] * n
                t[a], t[b] = va, vb
                if fill(t, 0) and tuple(t) not in seen:
                    seen.add(tuple(t))
                    full.append(t)
                    covered |= {(x, t[x], y, t[y]) for x, y in itertools.combinations(range(n), 2)}
    # Report every valid pair no row ended up holding (higher-arity
    # constraints can make a pairwise-valid pair unreachable).
    covered = set()
    for t in full:
        for a, b in itertools.combinations(range(n), 2):
            covered.add((a, t[a], b, t[b]))
    for a, b in itertools.combinations(range(n), 2):
        for va in range(sizes[a]):
            for vb in range(sizes[b]):
                if pair_ok(a, va, b, vb) and (a, va, b, vb) not in covered:
                    uncoverable.append((a, va, b, vb))
    return full, uncoverable


def builder_conflicts():
    """Value pairs whose overlays set the same knob differently (the tripolar
    geometry IS a grid, so it cannot meet `grid=spherical`).  A property of
    these tables, not a model rule -- excluded structurally."""
    out = []
    for (a, b) in itertools.combinations(AXIS_NAMES, 2):
        for va in VALUE_NAMES[a]:
            for vb in VALUE_NAMES[b]:
                probe = dict(BASE_CELL, **{a: va, b: vb})
                try:
                    merged_namelist(probe)
                except BuilderConflict:
                    out.append(((a, va), (b, vb)))
    return out


def generate_cells(forbidden=(), seed=SEED):
    """IPOG over AXES with `forbidden` tuples of (axis_name, value_name),
    plus the structural `builder_conflicts`."""
    sizes = [len(VALUE_NAMES[a]) for a in AXIS_NAMES]
    fidx = []
    for tup in list(forbidden) + builder_conflicts():
        fidx.append(tuple((AXIS_INDEX[a], VALUE_NAMES[a].index(v)) for a, v in tup))
    rows, unc = ipog(sizes, fidx, seed)
    cells = [{AXIS_NAMES[a]: VALUE_NAMES[AXIS_NAMES[a]][r[a]] for a in range(len(r))}
             for r in rows]
    uncoverable = [((AXIS_NAMES[a], VALUE_NAMES[AXIS_NAMES[a]][va]),
                    (AXIS_NAMES[b], VALUE_NAMES[AXIS_NAMES[b]][vb])) for a, va, b, vb in unc]
    return cells, uncoverable


# ===========================================================================
# 4. Running the model
# ===========================================================================
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?|[-+]?NaN|[-+]?Infinity|[-+]?Inf"
# The console formats, as `stability.py` parses them (kept in step with it).
_STATS_RE = re.compile(
    r"\[stats\]\s+Day\s+(?P<day>" + _NUM + r")\s+step\s+(?P<step>\d+)"
    r"\s+En\s+(?P<en>" + _NUM + r")"
    r"(?:\s+MaxCFL\s+(?P<cfl>" + _NUM + r"))?"
    r"(?:\s+Mass\s+(?P<mass>" + _NUM + r"))?"
    r"(?:\s+Salt\s+(?P<salt>" + _NUM + r"))?"
    r"(?:\s+Temp\s+(?P<temp>" + _NUM + r"))?")
_BUDGET_RE = re.compile(
    r"^\s+(?P<what>Mass|Salt|Heat)\s+:\s+(?P<total>" + _NUM + r")"
    r"\s+Error\s+(?P<err>" + _NUM + r")")
_TOTAL_STEPS_RE = re.compile(r"^\s*Total steps:\s*(\d+)", re.MULTILINE)
_NAN_CATCH_RE = re.compile(r"\[nan-catch\] (?:outer step (\d+)|stage \d+ step (\d+))")
_I1_RE = re.compile(r"\[I1'\] vanished-layer invariant violated at outer step (\d+)")
_PRECOND_RE = re.compile(r"ALE remap preconditions violated at outer step\s+(\d+)")
_CRASH_RE = re.compile(r"error stop|segmentation fault|floating point exception|"
                       r"program received signal|\baborted\b", re.IGNORECASE)
_NOISE_RE = re.compile(r"^\s*$|^\[gpu-bind\]|^VALIDATE-ONLY:|^ERROR STOP 3\s*$|"
                       r"^Configuration validation failed|^Note: The following "
                       r"floating-point|^IEEE_")

# Budget tolerance: every cell carries a surface flux (and two of the three
# geometries an open or relaxed edge), so the model's OWN residuals are held
# to the open/forced band `stability_manifest.BUDGET_OPEN` uses.
BUDGET_TOL = {"Mass": 1.0e-9, "Salt": 1.0e-9, "Heat": 1.0e-9}

# The checks a record carries, in the order the plan runs them.  Phase 2
# fills 1-3; 4-6 stay "not_run" until their legs land (the record format is
# fixed now so they slot in without a schema change).
CHECKS = ("validate", "run", "budget", "decomp", "restart", "cross_backend")


def _f(tok):
    t = tok.strip().replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        low = t.lower().lstrip("+-")
        if low.startswith("nan"):
            return float("nan")
        if low.startswith("inf"):
            return float("-inf") if t.startswith("-") else float("inf")
        raise


def _finite(x):
    return x is not None and not math.isnan(x) and not math.isinf(x)


def find_binary(build_dir=None, binary=None):
    if binary:
        return os.path.abspath(binary)
    cands = [build_dir] if build_dir else ["build_gfortran", "build_gcc", "build"]
    for b in cands:
        p = os.path.join(b if os.path.isabs(b) else os.path.join(REPO_ROOT, b), "rdb")
        if os.path.isfile(p):
            return p
    raise SystemExit("rdb binary not found (pass --build-dir or --binary)")


def prepare_cell_dir(cell, cid, root, nml_patch=None, header_extra=""):
    """Write the cell's namelist + inputs into root/cid; return the nml path."""
    d = os.path.join(root, cid)
    os.makedirs(d, exist_ok=True)
    write_domain_inputs(d)
    nml = merged_namelist(cell)
    for grp, kv in (nml_patch or {}).items():
        nml.setdefault(grp, {}).update(kv)
    path = os.path.join(d, "cell.nml")
    with open(path, "w") as fh:
        fh.write(render_namelist(nml, cell_header(cell, cid, header_extra)))
    return path


def _run(binary, args, cwd, timeout):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    t0 = time.time()
    try:
        p = subprocess.run([binary] + args, cwd=cwd, env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, timeout=timeout)
        rc, out, err = p.returncode, p.stdout.decode("utf-8", "replace"), \
            p.stderr.decode("utf-8", "replace")
    except subprocess.TimeoutExpired as exc:
        rc = "timeout"
        out = (exc.stdout or b"").decode("utf-8", "replace")
        err = (exc.stderr or b"").decode("utf-8", "replace")
    return rc, out, err, time.time() - t0


def validate_cell(binary, cell, cid, root, timeout=60):
    """`rdb --validate-only` on the cell, at log_level "error": stdout is
    then exactly the refusal reasons.  Returns a dict."""
    path = prepare_cell_dir(cell, cid, root,
                            nml_patch={"logging_nml": {"log_level": "error"}})
    rc, out, err, wall = _run(binary, ["--validate-only", os.path.basename(path)],
                              os.path.dirname(path), timeout)
    lines = [ln.strip() for ln in (out + "\n" + err).splitlines()]
    msgs = [ln for ln in lines if not _NOISE_RE.search(ln)
            and not ln.startswith("#") and "Error termination" not in ln
            and not ln.startswith("at ")]
    stage = None
    for ln in lines:
        m = re.match(r"VALIDATE-ONLY: REFUSED by (\S+):", ln)
        if m:
            stage = m.group(1)
    if rc == 0:
        status = "accepted"
    elif rc == 3 or (isinstance(rc, int) and rc in (1, 2) and re.search(r"ERROR STOP", err)
                     and not re.search(r"signal|segmentation", err, re.I)):
        status = "refused"
    else:
        status = "crashed"
    return {"status": status, "rc": rc, "stage": stage, "messages": msgs, "wall_s": round(wall, 3)}


def parse_run(text):
    stats, budget = [], {}
    for line in text.splitlines():
        m = _STATS_RE.search(line)
        if m:
            stats.append({k: (_f(m.group(g)) if m.group(g) is not None else None)
                          for k, g in (("day", "day"), ("En", "en"), ("MaxCFL", "cfl"),
                                       ("Mass", "mass"), ("Salt", "salt"), ("Temp", "temp"))})
            stats[-1]["step"] = int(m.group("step"))
            continue
        m = _BUDGET_RE.match(line)
        if m:
            budget.setdefault(m.group("what"), []).append(_f(m.group("err")))
    total = None
    for m in _TOTAL_STEPS_RE.finditer(text):
        total = int(m.group(1))
    nan_steps = [int(a or b) for a, b in _NAN_CATCH_RE.findall(text)]
    i1 = _I1_RE.search(text)
    pre = _PRECOND_RE.search(text)
    crash = _CRASH_RE.search(text)
    stop = re.search(r"ERROR STOP\s*(.*)", text)
    return {"stats": stats, "budget": budget, "total_steps": total,
            "error_stop": stop.group(1).strip()[:200] if stop else None,
            "nan_catch_steps": nan_steps,
            "i1_step": int(i1.group(1)) if i1 else None,
            "precond_step": int(pre.group(1)) if pre else None,
            "crash_marker": crash.group(0) if crash else None}


def classify_run(rc, p):
    """Check 2 then check 3.  Returns (outcome, detail, first_bad_step)."""
    last = p["stats"][-1]["step"] if p["stats"] else 0
    bad_scalar = None
    for s in p["stats"]:
        for k in ("En", "MaxCFL", "Mass", "Salt", "Temp"):
            if s.get(k) is not None and not _finite(s[k]):
                bad_scalar = (k, s["step"])
                break
        if bad_scalar:
            break
    steps = [x for x in ([p["i1_step"], p["precond_step"]] + p["nan_catch_steps"]
                         + ([bad_scalar[1]] if bad_scalar else [])) if x]
    first = min(steps) if steps else None
    if rc != 0 or p["total_steps"] != N_STEPS:
        why = []
        if p["i1_step"]:
            why.append("I1' tripwire at step {}".format(p["i1_step"]))
        if p["precond_step"]:
            why.append("remap preconditions at step {}".format(p["precond_step"]))
        if p["nan_catch_steps"]:
            why.append("nan-catch from step {}".format(min(p["nan_catch_steps"])))
        if bad_scalar:
            why.append("{} non-finite at step {}".format(*bad_scalar))
        if p["error_stop"] is not None:
            why.append("ERROR STOP {}".format(p["error_stop"]))
        elif p["crash_marker"]:
            why.append("'{}'".format(p["crash_marker"]))
        why.append("exit {} after {} of {} steps".format(rc, p["total_steps"] or last, N_STEPS))
        return "CRASH", "; ".join(why), first or (last + 1)
    if p["nan_catch_steps"] or bad_scalar:
        why = []
        if p["nan_catch_steps"]:
            why.append("nan-catch zeroed non-finite faces from step {}".format(
                min(p["nan_catch_steps"])))
        if bad_scalar:
            why.append("{} non-finite at step {}".format(*bad_scalar))
        return "NONFINITE", "; ".join(why), first
    worst = {}
    for what, tol in BUDGET_TOL.items():
        errs = [e for e in p["budget"].get(what, []) if e is not None]
        if not errs:
            continue
        w = max(errs, key=lambda e: abs(e) if _finite(e) else float("inf"))
        worst[what] = w
        if not _finite(w) or abs(w) > tol:
            return "BUDGET", "{} residual {:.3e} > {:.0e}".format(what, w, tol), None
    if not p["budget"]:
        return "BUDGET", "no budget lines printed", None
    return "PASS", "En {:.3e}  MaxCFL {:.3g}  worst budget {}".format(
        p["stats"][-1]["En"] if p["stats"] else float("nan"),
        max((s["MaxCFL"] or 0.0) for s in p["stats"]) if p["stats"] else float("nan"),
        "  ".join("{} {:.1e}".format(k, v) for k, v in sorted(worst.items()))), None


def chksum_attribution(binary, cell, cid, root, step, timeout):
    """Re-run a crashed cell with the chksum probe over its last few steps and
    name the (step, stage, phase, field) that MINTED the first non-finite."""
    import read_chksum   # tools/read_chksum.py
    lo = max(1, step - 2)
    path = prepare_cell_dir(cell, cid + "_chk", root, nml_patch={
        "ocean_debug_nml": {"chksum": True, "chksum_start_step": lo,
                            "chksum_end_step": step}},
        header_extra="chksum re-run, steps {}..{}".format(lo, step))
    rc, out, err, _ = _run(binary, [os.path.basename(path)], os.path.dirname(path), timeout)
    log = os.path.join(os.path.dirname(path), "run.log")
    with open(log, "w") as fh:
        fh.write(out + err)
    for rec in read_chksum.parse_log(log):
        if rec.get("kind") == "chksum" and rec.get("nonfin", 0) > 0:
            return "first non-finite minted at step {} s{} phase {} field {} (nonfin {})".format(
                rec["step"], rec["stage"], rec["phase"], rec["field"], rec["nonfin"])
    return "chksum window {}..{}: no non-finite row (corruption is finite or later)".format(lo, step)


def run_metrics(p):
    """The per-cell numbers a later leg compares across backends / rank
    counts / a restart: the final console scalars, the CFL peak and the
    worst budget residuals (None where the run did not print them)."""
    last = p["stats"][-1] if p["stats"] else {}
    worst = {}
    for what, errs in p["budget"].items():
        errs = [e for e in errs if e is not None]
        if errs:
            worst[what] = max(errs, key=lambda e: abs(e) if _finite(e) else float("inf"))
    cfl = [s["MaxCFL"] for s in p["stats"] if s.get("MaxCFL") is not None]

    def _j(x):   # JSON has no NaN/Inf: keep them as strings
        return x if x is None or _finite(x) else repr(x)
    return {"step": last.get("step"), "En": _j(last.get("En")), "Mass": _j(last.get("Mass")),
            "Salt": _j(last.get("Salt")), "Temp": _j(last.get("Temp")),
            "MaxCFL_max": _j(max(cfl)) if cfl else None,
            "budget_worst": {k: _j(v) for k, v in sorted(worst.items())}}


def run_cell(binary, cell, cid, root, timeout=300, attribute=True):
    path = prepare_cell_dir(cell, cid, root)
    rc, out, err, wall = _run(binary, [os.path.basename(path)], os.path.dirname(path), timeout)
    with open(os.path.join(os.path.dirname(path), "run.log"), "w") as fh:
        fh.write(out + err)
    p = parse_run(out + "\n" + err)
    outcome, detail, first = classify_run(rc if rc != "timeout" else -1, p)
    if rc == "timeout":
        outcome, detail = "CRASH", "timeout after {} s".format(timeout)
    res = {"outcome": outcome, "detail": detail, "rc": rc, "wall_s": round(wall, 2),
           "steps": p["total_steps"], "metrics": run_metrics(p)}
    if outcome in ("CRASH", "NONFINITE") and attribute and first:
        try:
            res["attribution"] = chksum_attribution(binary, cell, cid, root, first, timeout)
        except Exception as exc:   # attribution is a hint, never a verdict
            res["attribution"] = "attribution failed: {}".format(exc)
    return res


# ===========================================================================
# 5. Classification against compat_expect
# ===========================================================================
def classify(cell, nml, val, run):
    """The record's verdict.  Returns (cls, rows, note).

    cls: PASS | REFUSED_PHYSICAL | REFUSED_GAP | XFAIL | XPASS | FAIL
    """
    feats = compat_expect.features(nml)
    refuse_rows = [r for r in compat_expect.ROWS if r.kind == "refused" and r.matches(feats)]
    if val["status"] == "crashed":
        return "FAIL", [], "--validate-only crashed (rc {}): {}".format(
            val["rc"], " | ".join(val["messages"][-3:]))
    if val["status"] == "refused":
        unexplained, used = [], []
        for msg in val["messages"]:
            # The NARROWEST explaining row: it names the tuple to forbid.
            hit = sorted((r for r in refuse_rows if r.explains(msg)),
                         key=lambda r: (len(r.when) + len(r.unless), r.rid))
            if hit:
                if hit[0] not in used:
                    used.append(hit[0])
            else:
                unexplained.append(msg)
        if unexplained or not used:
            return "FAIL", used, "UNEXPECTED refusal: " + " | ".join(
                m[:240] for m in (unexplained or val["messages"]))
        cls = "REFUSED_PHYSICAL" if all(r.cls == "PHYSICAL" for r in used) else "REFUSED_GAP"
        return cls, used, "expected refusal ({})".format(", ".join(r.rid for r in used))
    # accepted
    if refuse_rows:
        return "XPASS", refuse_rows, "accepted, but row(s) {} say it is refused: delete them".format(
            ", ".join(r.rid for r in refuse_rows))
    run_rows = [r for r in compat_expect.ROWS if r.kind == "runtime" and r.matches(feats)]
    if run is None:
        return "FAIL", [], "accepted but not run"
    if run["outcome"] == "PASS":
        strict = [r for r in run_rows if r.scope == "cell"]
        if strict:
            return "XPASS", strict, "passes, but row(s) {} say it fails: delete them".format(
                ", ".join(r.rid for r in strict))
        return "PASS", [], run["detail"]
    # A row explains a failure only if the OUTCOME and its SIGNATURE (the
    # crash text / budget line) both match -- a different crash on a cell a
    # row covers is still a FAIL.
    hit = [r for r in run_rows if run["outcome"] in r.expect and r.explains(run["detail"])]
    if hit:
        return "XFAIL", hit, "{}: {}".format(run["outcome"], run["detail"])
    return "FAIL", [], "{}: {}".format(run["outcome"], run["detail"])


def witness_cell(row):
    """The full cell a `scope="any"` row pins (its axes over BASE_CELL)."""
    return dict(BASE_CELL, **row.witness)


def row_xpasses(records):
    """Row-level XPASS for `scope="any"` rows: the gap bites only in some
    combinations, so the row pins a witness cell that must fail with its
    signature on every run; a witness that passes (or fails for another
    reason) means the row is stale.  Returns [(rid, witness class)]."""
    out = []
    for row in compat_expect.ROWS:
        if row.scope != "any":
            continue
        key = cell_key(witness_cell(row))
        rec = [r for r in records if cell_key(r["axes"]) == key]
        if not rec or not (rec[0]["class"] == "XFAIL" and row.rid in rec[0]["rows"]):
            out.append((row.rid, rec[0]["class"] if rec else "not evaluated"))
    return out


def forbidden_tuple(cell, nml, rows):
    """The (axis, value) tuple an expected refusal forbids.

    For each feature the explaining rows require, the RESPONSIBLE axes are
    those where some other value of that axis alone turns the feature off.
    A feature two axes both supply (GM from `eddy` and from
    `lateral=meke_backscatter`) has no single responsible axis; then every
    axis whose own overlay carries it is taken.  Nothing found (an emergent
    feature) forbids the whole cell, which is always sound."""
    need, unless = set(), set()
    for r in rows:
        need |= set(r.when)
        unless |= set(r.unless)
    ref = compat_expect.features(BASE)
    tup = set()
    # An `unless` feature is part of the cause by its ABSENCE: the axes
    # whose other values would turn it on belong to the tuple too.
    for f in sorted(unless):
        for a in AXIS_NAMES:
            for w in VALUE_NAMES[a]:
                if w == cell[a]:
                    continue
                try:
                    alt = compat_expect.features(merged_namelist(dict(cell, **{a: w})))
                except BuilderConflict:
                    continue
                if alt[f]:
                    tup.add((a, cell[a]))
                    break
    for f in sorted(need):
        resp = []
        for a in AXIS_NAMES:
            for w in VALUE_NAMES[a]:
                if w == cell[a]:
                    continue
                try:
                    alt = compat_expect.features(merged_namelist(dict(cell, **{a: w})))
                except BuilderConflict:
                    continue
                if not alt[f]:
                    resp.append(a)
                    break
        if not resp:
            for a in AXIS_NAMES:
                solo = {g: dict(kv) for g, kv in BASE.items()}
                for g, kv in OVERLAY[(a, cell[a])].items():
                    solo.setdefault(g, {}).update(kv)
                if compat_expect.features(solo)[f] and not ref[f]:
                    resp.append(a)
        if not resp:
            return tuple(sorted((a, cell[a]) for a in AXIS_NAMES))
        tup |= {(a, cell[a]) for a in resp}
    return tuple(sorted(tup))


# ===========================================================================
# 6. The driver
# ===========================================================================
_PRINT_LOCK = threading.Lock()


def _log(msg):
    with _PRINT_LOCK:
        print(msg, flush=True)


def _pmap(fn, items, jobs):
    if jobs <= 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(fn, items))


def evaluate_cell(binary, cell, root, timeout, attribute):
    """Checks 1-3 on one cell -> (val, run, cls, rows, note).  Scratch lives
    in root/<hash>; a refused cell is never run."""
    cid = "h" + cell_hash(cell)
    val = validate_cell(binary, cell, cid, os.path.join(root, "validate"))
    run = None
    if val["status"] == "accepted":
        run = run_cell(binary, cell, cid, os.path.join(root, "run"), timeout=timeout,
                       attribute=attribute)
    cls, rows, note = classify(cell, merged_namelist(cell), val, run)
    return val, run, cls, rows, note


def _forbids(cls, rows):
    """Does this verdict take its cause out of the covering array?  Expected
    refusals do; so does a DETERMINISTIC runtime gap (every matching cell
    fails), because every other pair packed into such a cell is untested.
    A `scope="any"` gap bites only in some combinations: it stays in."""
    if cls in ("REFUSED_PHYSICAL", "REFUSED_GAP"):
        return True
    return cls == "XFAIL" and all(r.scope == "cell" for r in rows)


def build_matrix(binary, root, jobs, seed=SEED, timeout=300, attribute=True, max_iter=40):
    """The fixed point the plan prescribes: generate -> evaluate -> forbid the
    tuple behind every expected refusal and every deterministic XFAIL ->
    regenerate, until no new tuple appears.  The pairs a forbidden tuple
    makes unreachable are reported as `uncoverable`.

    Returns (final cells, results by cell key, forbidden tuples, iterations,
    uncoverable pairs, witness cells -- one per forbidden tuple, the cell
    that produced it, kept in the report so the gap is still proved)."""
    results, forbidden, witness = {}, set(), {}
    it = 0
    while True:
        it += 1
        cells, unc = generate_cells(sorted(forbidden), seed)
        todo = [c for c in cells if cell_key(c) not in results]

        def _e(c):
            return cell_key(c), evaluate_cell(binary, c, root, timeout, attribute)

        for k, res in _pmap(_e, todo, jobs):
            results[k] = res
        new = 0
        for c in cells:
            val, run, cls, rows, note = results[cell_key(c)]
            if not _forbids(cls, rows):
                continue
            # Each explaining row is an independent, sufficient cause (one
            # refusal message each), so each forbids its OWN minimal tuple.
            for row in rows:
                tup = forbidden_tuple(c, merged_namelist(c), [row])
                if tup not in forbidden:
                    forbidden.add(tup)
                    witness[tup] = c
                    new += 1
        _log("  iteration {}: {} cells, {} newly evaluated, {} forbidden tuples (+{})".format(
            it, len(cells), len(todo), len(forbidden), new))
        if not new or it >= max_iter:
            break
    return cells, results, sorted(forbidden), it, unc, [witness[t] for t in sorted(witness)]


def make_record(cid, cell, val, run, cls, rows, note, backend, role):
    nml = merged_namelist(cell)
    checks = {c: {"status": "not_run"} for c in CHECKS}
    checks["validate"] = {"status": val["status"], "stage": val["stage"],
                          "messages": val["messages"][:6], "wall_s": val["wall_s"]}
    matching = []
    if run is not None:
        checks["run"] = {"status": run["outcome"] if run["outcome"] in ("CRASH", "NONFINITE")
                         else "PASS", "detail": run["detail"], "steps": run["steps"],
                         "wall_s": run["wall_s"]}
        checks["run"]["metrics"] = run["metrics"]
        if "attribution" in run:
            checks["run"]["attribution"] = run["attribution"]
        if run["outcome"] in ("PASS", "BUDGET"):
            checks["budget"] = {"status": run["outcome"], "detail": run["detail"]}
        feats = compat_expect.features(nml)
        matching = [r.rid for r in compat_expect.ROWS if r.kind == "runtime" and r.matches(feats)]
    return {"id": cid, "hash": cell_hash(cell), "role": role, "axes": dict(cell),
            "backend": backend, "ranks": "1", "class": cls, "rows": [r.rid for r in rows],
            "note": note, "matching_rows": matching, "checks": checks}


def cmd_run(args):
    binary = find_binary(args.build_dir, args.binary)
    root = os.path.abspath(args.scratch_root)
    if os.path.isdir(root):
        shutil.rmtree(root)
    os.makedirs(root, exist_ok=True)
    t0 = time.time()
    _log("compat matrix: binary {}\n  axes: {}".format(binary, ", ".join(
        "{}({})".format(a, len(VALUE_NAMES[a])) for a in AXIS_NAMES)))
    cells, results, forbidden, iters, unc, witnesses = build_matrix(
        binary, root, args.jobs, args.seed, timeout=args.timeout,
        attribute=not args.no_attribution)
    wall = time.time() - t0
    pinned = [witness_cell(r) for r in compat_expect.ROWS if r.scope == "any"]
    todo = [c for c in pinned if cell_key(c) not in results]
    for k, res in _pmap(lambda c: (cell_key(c), evaluate_cell(
            binary, c, root, args.timeout, not args.no_attribution)), todo, args.jobs):
        results[k] = res
    wall = time.time() - t0
    records, seen = [], set()
    for role, group, prefix in (("cover", cells, "c"), ("witness", witnesses, "w"),
                                ("pinned", pinned, "p")):
        for k, c in enumerate(group):
            if cell_key(c) in seen:
                continue
            seen.add(cell_key(c))
            val, run, cls, rows, note = results[cell_key(c)]
            records.append(make_record("{}{:03d}".format(prefix, k), c, val, run, cls, rows,
                                       note, args.backend, role))
    summary = summarise(records, forbidden, unc, iters, results, wall)
    summary["row_xpass"] = row_xpasses(records)
    report = {"schema": 1, "seed": args.seed, "binary": binary, "backend": args.backend,
              "axes": {a: VALUE_NAMES[a] for a in AXIS_NAMES},
              "forbidden": [list(map(list, t)) for t in forbidden],
              "uncoverable": [list(map(list, u)) for u in unc],
              "summary": summary, "cells": records}
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=1, sort_keys=True)
    print_report(report, args.previous)
    bad = (summary["counts"].get("FAIL", 0) + summary["counts"].get("XPASS", 0)
           + len(summary["row_xpass"]))
    if not args.keep_scratch and bad == 0:
        shutil.rmtree(root, ignore_errors=True)
    return 1 if bad else 0


def summarise(records, forbidden, unc, iters, results, wall):
    counts = {}
    for r in records:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    v_sum = sum(res[0]["wall_s"] for res in results.values())
    r_sum = sum(res[1]["wall_s"] for res in results.values() if res[1] is not None)
    return {"counts": counts, "cells": len(records),
            "cover_cells": sum(1 for r in records if r["role"] == "cover"),
            "witness_cells": sum(1 for r in records if r["role"] == "witness"),
            "run_cells": sum(1 for r in records if r["checks"]["run"]["status"] != "not_run"),
            "evaluated_cells": len(results),
            "forbidden_tuples": len(forbidden), "uncoverable_pairs": len(unc),
            "iterations": iters, "wall_total_s": round(wall, 1),
            "cpu_validate_s_serial_sum": round(v_sum, 1),
            "cpu_run_s_serial_sum": round(r_sum, 1)}


def print_report(report, previous=None):
    s = report["summary"]
    print("\n=== compat matrix: {} covering cells + {} witnesses ({} run); {} forbidden "
          "tuples after {} iterations ({} cells evaluated)".format(
              s["cover_cells"], s["witness_cells"], s["run_cells"], s["forbidden_tuples"],
              s["iterations"], s["evaluated_cells"]))
    print("    wall {} s; serial CPU sums: validate {} s, run {} s".format(
        s["wall_total_s"], s["cpu_validate_s_serial_sum"], s["cpu_run_s_serial_sum"]))
    for k in ("PASS", "REFUSED_PHYSICAL", "REFUSED_GAP", "XFAIL", "XPASS", "FAIL"):
        print("    {:17s} {}".format(k, s["counts"].get(k, 0)))
    for rid, cls in s.get("row_xpass", []):
        print("    ROW XPASS         {}: its pinned witness is {}, not the row's XFAIL -- "
              "delete or re-pin the row".format(rid, cls))
    if s["uncoverable_pairs"]:
        print("    {} pair(s) are reachable only through refused / XFAIL cells "
              "(listed in the JSON report)".format(s["uncoverable_pairs"]))
    for cls in ("FAIL", "XPASS"):
        rows = [r for r in report["cells"] if r["class"] == cls]
        if rows:
            print("\n--- {} ---".format(cls))
        for r in rows:
            print("  {} {}".format(r["id"], " ".join("{}={}".format(a, r["axes"][a])
                                                    for a in AXIS_NAMES)))
            print("      {}".format(r["note"]))
            att = r["checks"]["run"].get("attribution")
            if att:
                print("      chksum: {}".format(att))
    used = {}
    for r in report["cells"]:
        for rid in r["rows"]:
            used[rid] = used.get(rid, 0) + 1
    print("\n--- KNOWN_GAP rows (owner) ---")
    for row in compat_expect.ROWS:
        if row.cls != "KNOWN_GAP":
            continue
        print("  {:28s} {:3d} cell(s)  owner {:12s} {}".format(
            row.rid, used.get(row.rid, 0), row.owner, row.link))
    unused = [row.rid for row in compat_expect.ROWS if row.rid not in used]
    if unused:
        print("  (no cell exercised: {})".format(", ".join(unused)))
    if previous and os.path.isfile(previous):
        with open(previous) as fh:
            prev = {r["hash"]: r["class"] for r in json.load(fh)["cells"]}
        diff = [(r["id"], prev[r["hash"]], r["class"]) for r in report["cells"]
                if r["hash"] in prev and prev[r["hash"]] != r["class"]]
        print("\n--- diff vs {} ---".format(previous))
        for d in diff:
            print("  {}: {} -> {}".format(*d))
        if not diff:
            print("  no class changes on common cells")


def cmd_list(args):
    cells, unc = generate_cells((), args.seed)
    print("{} cells over {} axes (no forbidden tuples; `run` regenerates around the "
          "refusals it finds)".format(len(cells), len(AXIS_NAMES)))
    for a in AXIS_NAMES:
        print("  {:11s} {:2d}  {}".format(a, len(VALUE_NAMES[a]), ", ".join(VALUE_NAMES[a])))
    if args.verbose:
        for k, c in enumerate(cells):
            print("c{:03d} {}".format(k, " ".join(c[a] for a in AXIS_NAMES)))
    return 0


def cmd_emit(args):
    root = os.path.abspath(args.out)
    if args.cell == "base":
        cell = dict(BASE_CELL)
    else:
        cells, _ = generate_cells((), args.seed)
        cell = cells[int(args.cell.lstrip("c"))]
    for kv in args.set or []:
        a, v = kv.split("=", 1)
        if v not in VALUE_NAMES[a]:
            raise SystemExit("{}: no value {} (have {})".format(a, v, VALUE_NAMES[a]))
        cell[a] = v
    path = prepare_cell_dir(cell, args.cell, root)
    print(path)
    return 0


def cmd_domain(args):
    """Phase-1 viability: every vertical coordinate x both split schemes on
    every edge variant and grid, the base closures otherwise; validate + run.
    (Under the cavity the base runs PP81 without KPP: KPP is refused on a
    covered z_fixed column, and this sweep is about the coordinate.)"""
    binary = find_binary(args.build_dir, args.binary)
    root = os.path.abspath(args.scratch_root)
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root, exist_ok=True)
    cases = []
    for g in VALUE_NAMES["geometry"]:
        for gr in VALUE_NAMES["grid"]:
            for sp in VALUE_NAMES["split"]:
                for v in VALUE_NAMES["vcoord"]:
                    c = dict(BASE_CELL, vcoord=v, geometry=g, grid=gr, split=sp)
                    if g == "cavity":
                        c["vmix_bl"] = "pp81"
                    try:
                        merged_namelist(c)
                    except BuilderConflict:
                        continue
                    cases.append(("d_{}_{}_{}_{}".format(v, sp, g, gr), c))

    def _go(item):
        cid, c = item
        val = validate_cell(binary, c, cid, os.path.join(root, "validate"))
        run = None
        if val["status"] == "accepted":
            run = run_cell(binary, c, cid, os.path.join(root, "run"), timeout=args.timeout,
                           attribute=not args.no_attribution)
        cls, rows, note = classify(c, merged_namelist(c), val, run)
        return cid, c, val, run, cls, note

    t0 = time.time()
    res = _pmap(_go, cases, args.jobs)
    print("{:48s} {:9s} {:16s} {}".format("case", "validate", "class", "detail"))
    for cid, c, val, run, cls, note in res:
        print("{:48s} {:9s} {:16s} {}".format(cid, val["status"], cls,
                                              (run["outcome"] + ": " + run["detail"]
                                               if run else note)[:160]))
        if run and run.get("attribution"):
            print("{:48s} {:9s} {:16s} chksum: {}".format("", "", "", run["attribution"]))
    # The coordinate summary: what each vertical coordinate does today.
    print("\nvcoord        " + "  ".join("{:>9s}".format(g[:9]) for g in VALUE_NAMES["geometry"]))
    for v in VALUE_NAMES["vcoord"]:
        row = []
        for g in VALUE_NAMES["geometry"]:
            outs = sorted({(run["outcome"] if run else val["status"].upper())
                           for cid, c, val, run, cls, note in res
                           if c["vcoord"] == v and c["geometry"] == g})
            row.append("{:>9s}".format("/".join(o[:4] for o in outs)))
        print("{:13s} {}".format(v, "  ".join(row)))
    print("wall {:.1f} s".format(time.time() - t0))
    return 0


# ===========================================================================
# 7. Who tests the test
# ===========================================================================
def _check(cond, what, fails):
    print("  {} {}".format("ok  " if cond else "FAIL", what))
    if not cond:
        fails.append(what)


def cmd_self_test(args):
    fails = []
    print("compat_matrix self-test (no model)")
    # IPOG covers every pair, deterministically, around forbidden tuples.
    sizes = [5, 4, 3, 3, 2, 2]
    rows, unc = ipog(sizes, [], 7)
    cov = {(a, t[a], b, t[b]) for t in rows for a, b in itertools.combinations(range(6), 2)}
    allp = {(a, va, b, vb) for a, b in itertools.combinations(range(6), 2)
            for va in range(sizes[a]) for vb in range(sizes[b])}
    _check(cov == allp and not unc, "ipog covers all {} pairs in {} rows".format(len(allp), len(rows)),
           fails)
    _check(rows == ipog(sizes, [], 7)[0], "ipog is deterministic for a seed", fails)
    forb = [((0, 1), (1, 2)), ((2, 0),)]
    rows, unc = ipog(sizes, forb, 7)
    _check(not any(_violates(t, forb) for t in rows), "ipog never emits a forbidden tuple", fails)
    cov = {(a, t[a], b, t[b]) for t in rows for a, b in itertools.combinations(range(6), 2)}
    want = {p for p in allp if not (p[0] == 0 and p[1] == 1 and p[2] == 1 and p[3] == 2)
            and not (p[0] == 2 and p[1] == 0) and not (p[2] == 2 and p[3] == 0)}
    _check(want <= cov and not unc, "ipog covers every pair a forbidden tuple leaves legal", fails)
    # The real axes build without conflicts, and the base cell is in range.
    cells, _ = generate_cells((), SEED)
    conflicts = []
    for c in cells:
        try:
            merged_namelist(c)
        except BuilderConflict as exc:
            conflicts.append(str(exc))
    _check(not conflicts, "{} generated cells merge without knob conflicts {}".format(
        len(cells), conflicts[:2]), fails)
    _check(cells == generate_cells((), SEED)[0], "the real cell list is seed-stable", fails)
    # An expected refusal forbids exactly the tuple that causes it.
    rows = {r.rid: r for r in compat_expect.ROWS}
    probes = [
        (dict(BASE_CELL, eddy="gm"), "closed_faces_gm",
         (("eddy", "gm"), ("vcoord", "z_fixed_cf"))),
        (dict(BASE_CELL, eddy="mle"), "mle_needs_epbl",       # `unless` joins the tuple
         (("eddy", "mle"), ("vmix_bl", "kpp"))),
        (dict(BASE_CELL, pv_adv="weno5"), "pv_weno_needs_sadourny",
         (("coriolis", "sadourny_energy"), ("pv_adv", "weno5"))),
        (dict(BASE_CELL, eddy="gm", lateral="meke_backscatter"), "closed_faces_gm",
         (("eddy", "gm"), ("lateral", "meke_backscatter"), ("vcoord", "z_fixed_cf"))),
    ]
    for cell, rid, want in probes:
        got = forbidden_tuple(cell, merged_namelist(cell), [rows[rid]])
        _check(got == want, "forbidden tuple for {} -> {}".format(rid, got), fails)
    # The NetCDF writer emits a well-formed CDF-1 header.
    d = os.path.join(os.path.abspath(args.scratch_root), "selftest")
    os.makedirs(d, exist_ok=True)
    write_domain_inputs(d)
    with open(os.path.join(d, BATHY_FILE), "rb") as fh:
        raw = fh.read()
    # magic+numrecs 8, dim list 8 + 2 x 12, no gatts 8, var list 8 + one
    # 40-byte var entry ("b", 2 dims) = 96 bytes of header, then the doubles.
    _check(raw[:4] == b"CDF\x01" and len(raw) == 96 + NX * NY * 8,
           "classic NetCDF bathymetry: magic + size ({} bytes)".format(len(raw)), fails)
    b = bathymetry()
    _check(sum(1 for r in b for x in r if x == 0.0) == 6 and min(min(r) for r in b if min(r) > 0)
           == STAIRCASE[0], "the domain has the island and the shelf", fails)
    # Every T/S column is statically stable (dT/dz dominates the front).
    unstable = 0
    for j in range(NY):
        for z0, z1 in zip(Z_SRC[:-1], Z_SRC[1:]):
            t0, s0 = ts_profile(j, z0)
            t1, s1 = ts_profile(j, z1)
            if (0.77 * (s1 - s0) - 0.2 * (t1 - t0)) <= 0.0:
                unstable += 1
    _check(unstable == 0, "the tilted-front IC is statically stable (linear EOS)", fails)
    # Classification: an expected refusal, an unexplained one, an XPASS.
    for name, fn in compat_expect.SELF_TESTS:
        ok, what = fn(sys.modules[__name__])
        _check(ok, "compat_expect: " + what, fails)
    rid_dupes = [r.rid for r in compat_expect.ROWS
                 if sum(1 for x in compat_expect.ROWS if x.rid == r.rid) > 1]
    _check(not rid_dupes, "row ids are unique {}".format(sorted(set(rid_dupes))), fails)
    bad_feats = [(r.rid, f) for r in compat_expect.ROWS for f in r.when + r.unless
                 if f not in compat_expect.FEATURES]
    _check(not bad_feats, "every row names known features {}".format(bad_feats[:3]), fails)
    shutil.rmtree(d, ignore_errors=True)
    print("{} failure(s)".format(len(fails)))
    return 1 if fails else 0


def cmd_validate_smoke(args):
    """The `--validate-only` contract, end to end (ctest `rdb_validate_only`):
    accepted -> exit 0 and NO output written; a validate_config refusal and
    an engine_setup (configure-stage) refusal -> exit 3 with the reason."""
    binary = find_binary(args.build_dir, args.binary)
    root = os.path.abspath(args.scratch_root)
    shutil.rmtree(root, ignore_errors=True)
    fails = []
    print("rdb --validate-only contract ({})".format(binary))
    cell = dict(BASE_CELL)
    # Diagnostics ON with an output directory: a real run would create the
    # directory, the per-rank NetCDF stream and the parameter-doc dumps.
    d = os.path.dirname(prepare_cell_dir(cell, "accept", root, nml_patch={
        "ocean_diag_nml": {"enabled": True, "filename": "smoke", "dt_out": DT},
        "output_nml": {"output_dir": "out"}}))
    before = set(os.listdir(d))
    t0 = time.time()
    rc, out, err, _ = _run(binary, ["--validate-only", "cell.nml"], d, 60)
    _check(rc == 0 and "VALIDATE-ONLY: ACCEPTED" in out,
           "the base cell is accepted (rc {})".format(rc), fails)
    _check(set(os.listdir(d)) == before, "an accepted validation writes no file ({})".format(
        sorted(set(os.listdir(d)) - before)), fails)
    _check("[stats]" not in out and "Total steps" not in out, "and takes no step", fails)
    print("    ({:.2f} s)".format(time.time() - t0))
    # Positive control: the SAME namelist run for one step does write --
    # so the no-file assertion above has teeth.
    d1 = os.path.dirname(prepare_cell_dir(cell, "control", root, nml_patch={
        "ocean_diag_nml": {"enabled": True, "filename": "smoke", "dt_out": DT},
        "output_nml": {"output_dir": "out"}, "time_nml": {"t_end": DT}}))
    rc, out, err, _ = _run(binary, ["cell.nml"], d1, 120)
    _check(rc == 0 and os.path.isdir(os.path.join(d1, "out")) and any(
        n.endswith(".nc") for n in os.listdir(os.path.join(d1, "out"))),
        "control: the same namelist WITHOUT --validate-only writes out/*.nc (rc {})".format(rc),
        fails)
    # validate_config refusal: a WENO PV scheme on the energy form.
    v = validate_cell(binary, dict(cell, coriolis="sadourny_energy", pv_adv="weno5"),
                      "refuse_validate", root)
    _check(v["status"] == "refused" and v["rc"] == 3 and v["stage"] == "validate_config"
           and any("pv_adv_scheme" in m for m in v["messages"]),
           "validate_config refusal -> rc 3, reason logged ({} / {})".format(v["rc"], v["messages"][:1]),
           fails)
    # engine_setup refusal: the configure-time viscous-CFL audit (a
    # PHYSICAL bound, so this contract does not move when a gap closes).
    d2 = os.path.dirname(prepare_cell_dir(cell, "refuse_engine", root, nml_patch={
        "ocean_hvisc_nml": {"nu_h": 1.0e8}, "logging_nml": {"log_level": "error"}}))
    rc, out, err, _ = _run(binary, ["--validate-only", "cell.nml"], d2, 60)
    _check(rc == 3 and "REFUSED by engine_setup" in out and "nu_h" in out,
           "engine_setup (configure-stage) refusal -> rc 3, reason logged (rc {})".format(rc), fails)
    _check(not any(n.endswith(".nc") and n not in (BATHY_FILE, TS_FILE) for n in os.listdir(d2)),
           "a refused validation writes no output file", fails)
    print("{} failure(s)".format(len(fails)))
    if not fails:
        shutil.rmtree(root, ignore_errors=True)
    return 1 if fails else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    def common(p, scratch):
        p.add_argument("--build-dir", default=None)
        p.add_argument("--binary", default=None)
        p.add_argument("--scratch-root", default=os.path.join(DEFAULT_SCRATCH, scratch))
        p.add_argument("--jobs", type=int, default=4)
        p.add_argument("--timeout", type=int, default=300)
        p.add_argument("--seed", type=int, default=SEED)
        p.add_argument("--no-attribution", action="store_true",
                       help="skip the chksum re-run of crashing cells")

    p = sub.add_parser("run", help="the full pairwise matrix, checks 1-3")
    common(p, "run")
    p.add_argument("--out", default=None, help="JSON report path")
    p.add_argument("--previous", default=None, help="last report, for a class diff")
    p.add_argument("--backend", default="cpu-gfortran")
    p.add_argument("--keep-scratch", action="store_true")
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("domain", help="every vcoord x geometry x grid on the base closures")
    common(p, "domain")
    p.set_defaults(fn=cmd_domain)
    p = sub.add_parser("list", help="the axes and the unconstrained cell list")
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(fn=cmd_list)
    p = sub.add_parser("emit", help="write one cell's namelist + inputs")
    p.add_argument("--cell", default="base", help="'base' or cNNN of the unconstrained list")
    p.add_argument("--set", action="append", help="axis=value override (repeatable)")
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=SEED)
    p.set_defaults(fn=cmd_emit)
    p = sub.add_parser("self-test", help="no model: generator, builder, classifier")
    p.add_argument("--scratch-root", default=DEFAULT_SCRATCH)
    p.set_defaults(fn=cmd_self_test)
    p = sub.add_parser("validate-smoke", help="the rdb --validate-only contract")
    common(p, "validate_smoke")
    p.set_defaults(fn=cmd_validate_smoke)
    args = ap.parse_args(argv)
    if not getattr(args, "fn", None):
        ap.print_help()
        return 2
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
