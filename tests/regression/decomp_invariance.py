#!/usr/bin/env python3
"""Decomposition-invariance gate: same case, 1 rank vs N, same answer.

**Why this exists.** Two real bugs shipped because they are INVISIBLE at a
single rank *by construction*:

  * `configure_ocean_land_mask` passed the global edge tag without AND-ing
    the per-rank `has_<edge>`, so every rank masked its own MPI seams as
    land. At np=1 `has_west` is `.true.`, so the buggy and fixed
    expressions are literally identical -- no single-rank test can differ.
  * `seed_baroclinic_jet_ic` built its coordinate maps from LOCAL extents.
    At np=1 local extents ARE the global extents, so again identical.

No amount of unit testing catches that class. The only gate that can is
running the SAME case on different rank counts and demanding the same
answer, which is what this does.

It also defends against a mistake I made by hand: an A/B measured against
a stale binary. This always builds nothing and runs one binary, so the two
legs cannot disagree about what code they are testing.

**Coverage is the point.** The sweep deliberately walks every `ic_config`
and every formula `topo_config`, because the seeders are exactly where
local-vs-global index bugs live. Adding a new IC or bathymetry formula
means adding a row here.

Usage:
    python3 tests/regression/decomp_invariance.py --build-dir build_mpi
    python3 tests/regression/decomp_invariance.py --ranks 4 --cases eady,spoon

Requires an MPI build (`-DRDB_ENABLE_MPI=ON`). Stdlib only -- this box
has no numpy.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toml_util  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# name -> (TOML input, n_steps).  One row per ic_config / topo_config so the
# sweep covers every seeder.  `n_steps` is small on purpose: a
# decomposition bug shows up in ONE step (both bugs above did), so a long
# run only buys chaotic divergence that muddies the comparison.
CASES = [
    # --- ic_config coverage ---
    ("ic_eady",        "validation_examples/ocean/eady/eady.toml", 3),
    ("ic_geoadj",      "validation_examples/ocean/geostrophic_adjustment/"
                       "geostrophic_adjustment.toml", 3),
    ("ic_bcjet",       "validation_examples/ocean/bc_inst/bc_inst_tuned_512.toml", 3),
    # --- topo_config coverage (ic_config = "" default) ---
    ("topo_spoon",     "validation_examples/ocean/double_gyre/double_gyre_mom6.toml", 3),
    ("topo_seamount",  "validation_examples/ocean/seamount/seamount_pred_corr.toml", 3),
    ("topo_island",    "validation_examples/ocean/island_at_rest/island_at_rest.toml", 3),
    ("topo_ddrake",    "validation_examples/ocean/double_drake/double_drake.toml", 3),
    ("topo_nw2",       "validation_examples/ocean/neverworld2/neverworld2.toml", 3),
    ("topo_flat",      "validation_examples/ocean/epbl_mld/epbl_basin.toml", 3),
]

STATS = re.compile(
    r"\[stats\].*?En\s+([0-9.eE+-]+)\s+MaxCFL\s+([0-9.eE+-]+)\s+Mass\s+([0-9.eE+-]+)")


def patch(text, n_steps, px, py, outdir):
    """Rewrite the TOML input for a short, output-isolated, decomposed run.

    Two settings are forced, on BOTH legs so the comparison stays valid:

    * `bt_halo = 0` -- the auto value scales with the domain and trips
      `ng_wide > min(nx_phys, ny_phys)` on small tiles.
    * `nghost >= 3` -- `ocean_halo_init` requires it for any decomposed
      run, and several shipped nmls use 2.

    Both are legitimate CONFIGURATION limits rather than the property
    under test; letting either fire would silently skip ranks.
    """
    dt = toml_util.get_float(text, "time_nml", "dt_fixed")
    over = {"ocean_bt_nml": {"bt_halo": 0},
            "mpi_nml": {"px": px, "py": py}}
    if dt is not None:
        over["time_nml"] = {"t_end": float("%.10g" % (n_steps * dt)),
                            "time_unit": "second"}
    ng = toml_util.get_float(text, "grid_nml", "nghost")
    if ng is not None:
        over["grid_nml"] = {"nghost": max(int(ng), 3)}
    if toml_util.get(text, "output_nml", "output_dir") is not None:
        over["output_nml"] = {"output_dir": outdir}
    if toml_util.get(text, "ocean_diag_nml", "enabled") is not None:
        over["ocean_diag_nml"] = {"enabled": False}   # no NetCDF churn
    return toml_util.patch(text, over)


def run(binary, nml_text, workdir, nranks):
    path = os.path.join(workdir, "case.toml")
    with open(path, "w") as fh:
        fh.write(nml_text)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        p = subprocess.run(["mpirun", "--oversubscribe", "-np", str(nranks),
                            binary, path],
                           cwd=workdir, env=env, capture_output=True,
                           text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return None, "timeout"
    hits = STATS.findall(p.stdout)
    if not hits:
        tail = (p.stdout + p.stderr).strip().splitlines()[-1:] or ["no output"]
        return None, tail[0][:110]
    return hits[-1], None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default="build_mpi")
    ap.add_argument("--ranks", type=int, default=4,
                    help="rank count for the decomposed leg (default 4)")
    ap.add_argument("--layout", default=None,
                    help="PXxPY for the decomposed leg (default: PX=ranks, PY=1)")
    ap.add_argument("--cases", default=None, help="comma-separated subset")
    ap.add_argument("--keep", action="store_true", help="keep scratch dirs")
    args = ap.parse_args()

    binary = os.path.join(REPO, args.build_dir, "rdb")
    if not os.path.isfile(binary):
        sys.exit("no rdb binary at %s -- configure with "
                 "-DRDB_ENABLE_MPI=ON" % binary)
    if args.layout:
        px, py = (int(v) for v in args.layout.lower().split("x"))
    else:
        px, py = args.ranks, 1
    if px * py != args.ranks:
        sys.exit("layout %dx%d != --ranks %d" % (px, py, args.ranks))

    cases = CASES
    if args.cases:
        want = set(args.cases.split(","))
        cases = [c for c in CASES if c[0] in want]

    scratch = tempfile.mkdtemp(prefix="decomp_inv_")
    print("1 rank vs %d ranks (%dx%d), binary %s\n" % (args.ranks, px, py, binary))
    print("  %-16s %-34s %-34s %s" % ("case", "np=1  En / MaxCFL",
                                      "x-split En / MaxCFL", "verdict"))
    print("  " + "-" * 104)
    failures = []
    for name, nml, steps in cases:
        src = os.path.join(REPO, nml)
        if not os.path.isfile(src):
            print("  %-16s SKIP (missing %s)" % (name, nml)); continue
        text = open(src).read()
        # BOTH an x-split and a y-split.  A single layout has a blind
        # spot: `seed_eady_ic` varies only in y, so an x-only
        # decomposition cannot expose its local-index bug (and did not).
        legs = [("ref", 1, 1, 1), ("x", px, py, args.ranks),
                ("y", py, px, args.ranks)]
        res, bad = {}, None
        for tag, a, b, n in legs:
            wd = os.path.join(scratch, "%s_%s" % (name, tag))
            os.makedirs(wd, exist_ok=True)
            res[tag], err = run(binary, patch(text, steps, a, b, wd), wd, n)
            if err:
                bad = "%s leg: %s" % (tag, err)
                break
        if bad:
            print("  %-16s ERROR  %s" % (name, bad)); failures.append(name); continue
        r = res["ref"]
        diffs = [t for t in ("x", "y") if res[t] != r]
        print("  %-16s %-34s %-34s %s"
              % (name, "%s / %s" % (r[0], r[1]),
                 "%s / %s" % (res["x"][0], res["x"][1]),
                 "ok" if not diffs else "*** DIFFERS (%s) ***" % ",".join(diffs)))
        if diffs:
            failures.append(name)
            print("  %-16s   np=1   : En %s  MaxCFL %s  Mass %s"
                  % ("", r[0], r[1], r[2]))
            for t in diffs:
                print("  %-16s   %s-split: En %s  MaxCFL %s  Mass %s"
                      % ("", t, res[t][0], res[t][1], res[t][2]))
    if not args.keep:
        shutil.rmtree(scratch, ignore_errors=True)
    print()
    if failures:
        print("FAILED (%d/%d): %s" % (len(failures), len(cases), ", ".join(failures)))
        print("A difference here is a decomposition bug, not roundoff: these runs "
              "are only %d steps and MaxCFL is a max reduction (order-invariant)."
              % cases[0][2])
        return 1
    print("all %d cases decomposition-invariant" % len(cases))
    return 0


if __name__ == "__main__":
    sys.exit(main())
