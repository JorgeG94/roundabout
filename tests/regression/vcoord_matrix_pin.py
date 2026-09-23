#!/usr/bin/env python3
"""Pin the vertical-coordinate matrix from MEASURED runs -- never by hand.

Reads `stability.py --out` JSON files (one per toolchain per tier), and
writes `vcoord_matrix_measured.py`: the per-family rx0 ENVELOPES of the
viscous leg and one MEASURED record per failing cell.  `vcoord_matrix.py`
turns each record into a scoped `known_failure`; the REASON text is policy
and lives there, the NUMBERS are data and live in the generated module.

Rules, stated once so a reader can check them against the output:

  * A record's `assertions` are the UNION of the failing assertions over
    every toolchain and every tier measured -- the tolerance band.  The bars
    themselves are physics bars and are never read from a run; only WHICH of
    them a cell fails is pinned, so a cell whose number drifts inside its
    band stays XFAIL on every toolchain, and a cell that fails a NEW
    assertion on any toolchain turns FAIL.
  * `tiers` lists the tiers at which the cell failed on ANY toolchain.
  * `toolchain_dependent` is set when, at some tier, one toolchain passed
    the cell outright and another did not: a pass is then one of the
    measured outcomes and reports PASS rather than XPASS.
  * `measured` quotes every toolchain's own number at every tier.
  * ENVELOPES (viscous leg, base scheme only): per family, the largest
    geometry rx0 such that EVERY viscous cell of that family at or below it
    passes on EVERY toolchain measured.  None when even rx0 = 0 fails.

Usage (after a tier-1 and a tier-2 sweep on each toolchain):

    python3 tests/regression/vcoord_matrix_pin.py \\
        --t1 gfortran=t1_gcc.json --t1 nvfortran=t1_gpu.json \\
        --t2 gfortran=t2_gcc.json --t2 nvfortran=t2_gpu.json \\
        --provenance "gfortran 15.1 CPU; nvfortran 25.5 V100 cc70; ..."

Stdlib only.
"""

import argparse
import json
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
OUT = os.path.join(THIS_DIR, "vcoord_matrix_measured.py")


def _load(spec):
    label, path = spec.split("=", 1)
    data = json.load(open(path))
    return label, {c["case"]: c for c in data["cases"]}


def _failed(case):
    return sorted(a["name"] for a in case["assertions"]
                  if not a["ok"] and not a["skipped"])


def _obs(case, name):
    for a in case["assertions"]:
        if a["name"] == name:
            return a["observed"]
    return ""


def _summary(case, dt):
    f = _failed(case)
    if not f:
        return "passes"
    if "completed" in f:
        m = re.search(r"at step (\d+)", _obs(case, "completed"))
        s = "aborts"
        if m:
            s += " day {:.1f}".format(int(m.group(1)) * dt / 86400.0)
        pc = _obs(case, "remap:preconditions")
        if pc and pc != "no violation":
            m = re.search(r"most negative h (\S+) m", pc)
            s += " on the remap guard (h {} m)".format(m.group(1) if m else "?")
        elif "finite" in f:
            s += " non-finite"
        return s
    m = re.search(r"peak En = (\S+)", _obs(case, "energy:rest"))
    s = "peak En {}".format(m.group(1).rstrip(",")) if m else "completes"
    if "energy:rest-growth-rate" in f:
        m = re.search(r"En e-folding ([\d.]+) days", _obs(case, "energy:rest-growth-rate"))
        if m:
            s += ", En e-folding {} d".format(m.group(1))
    if "conserve:Salt" in f:
        m = re.search(r"= (\S+)", _obs(case, "conserve:Salt"))
        if m:
            s += ", salt residual {}".format(m.group(1))
    return s


def _record(name, runs, dt):
    """The pinned record of ONE case name across toolchains and tiers."""
    rec = {"assertions": set(), "tiers": set(), "measured": {},
           "toolchain_dependent": False}
    for tier in (1, 2):
        outcomes = []
        for label, table in runs[tier]:
            run = table.get(name)
            if run is None:
                continue
            f = _failed(run)
            outcomes.append(bool(f))
            if f:
                rec["assertions"].update(f)
                rec["tiers"].add(tier)
            rec["measured"].setdefault(label, []).append(
                "tier {}: {}".format(tier, _summary(run, dt)))
        if outcomes and any(outcomes) and not all(outcomes):
            rec["toolchain_dependent"] = True
    if not rec["assertions"]:
        return None
    return {"assertions": sorted(rec["assertions"]),
            "tiers": sorted(rec["tiers"]),
            "toolchain_dependent": rec["toolchain_dependent"],
            "measured": {k: "; ".join(v) for k, v in sorted(rec["measured"].items())}}


def _write_table(fh, name, table):
    fh.write("{} = {{\n".format(name))
    for leg in ("inviscid", "viscous"):
        fh.write("    {!r}: {{\n".format(leg))
        for key, rec in sorted(table[leg].items()):
            fh.write("        {!r}: {{\n".format(key))
            fh.write("            \"assertions\": {!r},\n".format(rec["assertions"]))
            fh.write("            \"tiers\": {!r},\n".format(rec["tiers"]))
            fh.write("            \"toolchain_dependent\": {!r},\n".format(
                rec["toolchain_dependent"]))
            fh.write("            \"measured\": {\n")
            for lab, txt in sorted(rec["measured"].items()):
                fh.write("                {!r}: {!r},\n".format(lab, txt))
            fh.write("            },\n        },\n")
        fh.write("    },\n")
    fh.write("}\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--t1", action="append", default=[], metavar="LABEL=JSON")
    p.add_argument("--t2", action="append", default=[], metavar="LABEL=JSON")
    p.add_argument("--provenance", default="")
    p.add_argument("--out", default=OUT)
    args = p.parse_args(argv)

    import vcoord_matrix as vcm
    import manifest

    cases = {c["name"]: c for c in manifest.STABILITY_CASES if c.get("matrix")}
    runs = {1: [_load(s) for s in args.t1], 2: [_load(s) for s in args.t2]}

    measured = {"inviscid": {}, "viscous": {}}
    twins = {"inviscid": {}, "viscous": {}}
    for name, c in sorted(cases.items()):
        mx = c["matrix"]
        if mx["expect"] == "refused" or c.get("split_scheme"):
            continue
        key = vcm.cell_key(mx["problem"], mx["family"], mx["stratification"],
                           mx["eos"])
        for suffix, table_out in (("", measured), ("__ssp_rk2", twins)):
            rec = _record(name + suffix, runs, vcm.DT)
            if rec:
                table_out[mx["leg"]][key] = rec

    # ---- envelopes: viscous leg, base scheme, tier 1, every toolchain ----
    fam_cells = {}
    for name, c in cases.items():
        mx = c["matrix"]
        if mx["leg"] != "viscous" or c.get("split_scheme"):
            continue
        ok = True
        seen = False
        for _label, table in runs[1]:
            run = table.get(name)
            if run is None:
                continue
            seen = True
            ok = ok and not _failed(run)
        if seen:
            fam_cells.setdefault(mx["family"], []).append(
                (mx["rx0_geometry"], ok))
    envelopes = {}
    for fam, cells in sorted(fam_cells.items()):
        lim = None
        for rx0, ok in sorted(cells):
            if not ok:
                break
            lim = rx0
        envelopes[fam] = None if lim is None else round(lim, 4)

    with open(args.out, "w") as fh:
        fh.write('"""GENERATED by tests/regression/vcoord_matrix_pin.py -- DO NOT EDIT.\n\n'
                 "The MEASURED vertical-coordinate matrix: per-family viscous-leg rx0\n"
                 "ENVELOPES and one record per failing cell (union of failing\n"
                 "assertions over toolchains and tiers, every toolchain's number\n"
                 "quoted).  Re-measure and regenerate; never hand-edit a number.\n\n"
                 "Provenance: {}\n\"\"\"\n\n".format(args.provenance or "(not given)"))
        fh.write("ENVELOPES = {\n")
        for fam, lim in sorted(envelopes.items()):
            fh.write("    {!r}: {!r},\n".format(fam, lim))
        fh.write("}\n\n# Base cells (the default pred_corr).\n")
        _write_table(fh, "MEASURED", measured)
        fh.write("\n# The `__ssp_rk2` twins (tier 2 only), pinned separately so a\n"
                 "# twin-only failure does not loosen the base cell's marker.\n")
        _write_table(fh, "MEASURED_TWIN", twins)
    print("wrote {} ({} inviscid, {} viscous, {} twin records)".format(
        args.out, len(measured["inviscid"]), len(measured["viscous"]),
        len(twins["inviscid"]) + len(twins["viscous"])))
    for fam, lim in sorted(envelopes.items()):
        print("  envelope {:12s} rx0 <= {}".format(fam, lim))
    return 0


if __name__ == "__main__":
    sys.exit(main())
