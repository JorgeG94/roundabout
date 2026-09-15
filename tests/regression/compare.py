#!/usr/bin/env python3
"""P1 golden-summary regression compare (high-tolerance, NOT bitwise).

Layered on top of the P0 runner (run_regression.py): reuses its case-execution
machinery (temp-nml patch + output isolation + run-clean/NaN gate) and adds a
*golden field-summary* compare so the suite catches numerical DRIFT, not just
crashes.

What a golden is
----------------
A tiny per-case final-state SUMMARY -- NOT full fields. For each case we
capture, at the final step, the `min / max / mean` of every prognostic the run
prints as a `[diag]` line (SSH/eta, u, v, T, S, KE, any active tracer, ice_*)
plus the scalar diagnostics on the `[stats]` line (En, total Mass, and Salt /
Temp when thermodynamics are reported). These come straight from the model's
own console output -- we parse the LAST occurrence of each field.

Guaranteeing a final-state print
--------------------------------
The committed nmls print diags/stats on an hours-or-days cadence, so a short
regression run would otherwise only emit at t=0. run_regression.patch_namelist
(force_emit=True, the default) rewrites `&logging_nml status_interval` and
`&ocean_diag_nml dt_out` down to half an outer step (in seconds, after the
time_unit normalization) so BOTH fire every step -- the final step therefore
always emits a usable summary. This module relies on that.

Tolerance model (why not a plain relative delta)
------------------------------------------------
CPU (gfortran) and GPU (nvfortran) differ by more than bitwise -- FMA / reduction
order / transcendental libraries. A plain relative delta also explodes on fields
that are PHYSICALLY zero but carry roundoff noise (a quiescent seamount's
velocity min/max sit at ~1e-12 and CPU vs GPU differ by whole factors there,
which is meaningless). So each value is compared numpy-`isclose` style:

    allowance = atol + rtol * scale
    pass      = |golden - run| <= allowance

where `scale` is the FIELD's magnitude (for a diag field: max |min|,|max| over
both golden and run; for a stats scalar: max |golden|,|run|). A roundoff-zero
field has scale ~1e-12, so `allowance` collapses to `atol` and the noise is
absorbed; an active field is bounded by `rtol * scale`; a regression that
spins up 0.1 m/s where the golden was ~0 lifts the scale to 0.1 and fails
loudly. `rtol` is the tunable golden tolerance (default LOOSE, per-case
overridable via a manifest `tol` key); `atol` is a small SI floor.

Modes
-----
--update-golden : run the suite, (re)write tests/regression/golden/<case>.json
                  from the parsed summaries. Generate these from a trusted CPU
                  build and commit them.
--compare       : (default) run the suite, parse the same summaries, compare
                  each value to the committed golden within tolerance. A case
                  FAILS if any value drifts past its allowance OR the P0
                  run-clean/NaN gate trips OR a golden field is missing from
                  the run. Reports per-case PASS/FAIL with the worst-drifting
                  field + its relative delta.

Stdlib only (the repo forbids `pip install`).
"""

import argparse
import json
import math
import os
import queue
import re
import sys
import threading
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
GOLDEN_DIR = os.path.join(THIS_DIR, "golden")
sys.path.insert(0, THIS_DIR)
import manifest  # noqa: E402
import run_regression as rr  # noqa: E402  (reuse the P0 case-execution machinery)

# Default tolerances -- set from the MEASURED CPU(gfortran)-vs-GPU(nvfortran)
# spread over this corpus, not guessed (see README + the commit message):
#   * The worst REAL-field (non-roundoff) CPU<->GPU relative drift across the 22
#     well-behaved cases is 1.4e-6 (ideal_age_demo SSH:min). rtol=1e-3 is ~700x
#     that -- loose enough to absorb FMA / reduction-order / transcendental-lib
#     divergence (and GPU run-to-run wobble; measured here = zero, deterministic)
#     yet tight enough that a real physics regression (fields move by
#     percent-to-orders) fails loudly.
#   * atol is an absolute floor that absorbs PHYSICALLY-ZERO fields carrying
#     roundoff noise: a quiescent seamount's velocity min/max sit at ~8e-12 and
#     CPU vs GPU differ ~6% *relatively* there (|abs diff| ~5e-13), which is
#     meaningless. atol=1e-7 clears that noise by >1e5 while masking only
#     sub-micro-SI signal (negligible for a physics regression suite).
#   * Two early-transient baroclinic-instability cases (baroclinic_2layer, eady)
#     whose tiny cross-channel v / free-surface SSH EXTREMA diverge 9-15% at the
#     noise floor -- while their integrated En/Mass/Salt/Temp agree within rtol
#     -- carry a documented per-case `tol` in the manifest instead of inflating
#     this global bound (the plan's "flag it, don't hide it" rule).
DEFAULT_RTOL = 1e-3
DEFAULT_ATOL = 1e-7

# A number in Fortran console output: ES/E/D exponent forms, plain decimals,
# and the IEEE spellings (which the P0 gate catches first, but parse anyway).
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?|[-+]?NaN|[-+]?Inf(?:inity)?"

# [diag] t=<time> <name> [<units>]  min=<v>  max=<v>  mean=<v>
_DIAG_RE = re.compile(
    r"\[diag\]\s+t=\s*(?P<t>" + _NUM + r")\s+"
    r"(?P<name>.+?)\s+\[.*?\]\s+"
    r"min=\s*(?P<min>" + _NUM + r")\s+"
    r"max=\s*(?P<max>" + _NUM + r")\s+"
    r"mean=\s*(?P<mean>" + _NUM + r")"
)
# [stats] ... En <v> ... Mass <v> [ Salt <v>  Temp <v> ]
_STATS_LINE_RE = re.compile(r"\[stats\]\s")
_EN_RE = re.compile(r"\bEn\s+(" + _NUM + r")")
_MASS_RE = re.compile(r"\bMass\s+(" + _NUM + r")")
_SALT_RE = re.compile(r"\bSalt\s+(" + _NUM + r")")
_TEMP_RE = re.compile(r"\bTemp\s+(" + _NUM + r")")


def _to_float(tok):
    """Parse a Fortran numeric token to float (handles d-exponents / IEEE)."""
    t = tok.strip().replace("D", "E").replace("d", "e")
    return float(t)


# ---------------------------------------------------------------------------
# Summary parsing
# ---------------------------------------------------------------------------
def parse_summary(text):
    """Parse the final-state golden summary out of a run's combined stdout.

    Returns {"diag": {name: {"min","max","mean"}}, "stats": {En,Mass[,Salt,Temp]}}.
    For each diag FIELD we keep the LAST occurrence (final step). For stats we
    parse the LAST [stats] line. Missing sections yield empty dicts.
    """
    diag = {}
    for m in _DIAG_RE.finditer(text):
        name = m.group("name").strip()
        # Last occurrence wins (dict overwrite as we scan top->bottom).
        diag[name] = {
            "min": _to_float(m.group("min")),
            "max": _to_float(m.group("max")),
            "mean": _to_float(m.group("mean")),
        }

    stats = {}
    last_stats = None
    for line in text.splitlines():
        if _STATS_LINE_RE.search(line):
            last_stats = line
    if last_stats is not None:
        for key, rx in (("En", _EN_RE), ("Mass", _MASS_RE),
                        ("Salt", _SALT_RE), ("Temp", _TEMP_RE)):
            mm = rx.search(last_stats)
            if mm:
                stats[key] = _to_float(mm.group(1))

    return {"diag": diag, "stats": stats}


# ---------------------------------------------------------------------------
# Golden JSON I/O
# ---------------------------------------------------------------------------
def golden_path(case_name):
    return os.path.join(GOLDEN_DIR, case_name + ".json")


def write_golden(case, summary, backend):
    os.makedirs(GOLDEN_DIR, exist_ok=True)
    payload = {
        "case": case["name"],
        "backend_generated": backend,
        "n_steps": case["n_steps"],
        "summary": summary,
    }
    with open(golden_path(case["name"]), "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")


def load_golden(case_name):
    path = golden_path(case_name)
    if not os.path.isfile(path):
        return None
    with open(path, "r") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def _isclose(g, r, rtol, atol, scale):
    """numpy-style closeness on one value against a field-magnitude `scale`.

    Returns (passed, abs_diff, allowance, rel_drift). `rel_drift` is |g-r| over
    the field scale (with a tiny floor) -- a scale-free "how far it moved"
    number for reporting; `passed` uses the atol+rtol*scale allowance so
    roundoff-zero fields do not false-fail.
    """
    if math.isnan(g) or math.isnan(r) or math.isinf(g) or math.isinf(r):
        return (False, float("inf"), 0.0, float("inf"))
    abs_diff = abs(g - r)
    allowance = atol + rtol * scale
    rel_drift = abs_diff / max(scale, 1e-30)
    return (abs_diff <= allowance, abs_diff, allowance, rel_drift)


def compare_summary(golden_summary, run_summary, rtol, atol):
    """Compare a run summary to a golden summary.

    Returns (passed, worst) where `worst` describes the largest-drift value:
        {"field","stat","golden","run","abs_diff","allowance","rel_drift","ok"}
    plus a `problems` list of every value that exceeded its allowance (incl.
    golden fields entirely absent from the run).
    """
    problems = []
    worst = None

    def consider(field, stat, g, r, scale):
        nonlocal worst
        ok, abs_diff, allowance, rel_drift = _isclose(g, r, rtol, atol, scale)
        rec = {
            "field": field, "stat": stat, "golden": g, "run": r,
            "abs_diff": abs_diff, "allowance": allowance,
            "rel_drift": rel_drift, "ok": ok,
        }
        # "worst" ranks by how far past the allowance we are (drift/allowance),
        # so the reported field is the closest to (or furthest past) failing.
        margin = abs_diff / allowance if allowance > 0 else float("inf")
        if worst is None or margin > worst["_margin"]:
            worst = dict(rec, _margin=margin)
        if not ok:
            problems.append(rec)

    gdiag = golden_summary.get("diag", {})
    rdiag = run_summary.get("diag", {})
    for name, gvals in gdiag.items():
        if name not in rdiag:
            problems.append({
                "field": name, "stat": "*", "golden": None, "run": None,
                "abs_diff": float("inf"), "allowance": 0.0,
                "rel_drift": float("inf"), "ok": False,
                "note": "field absent from run",
            })
            worst = worst or {}
            worst = {"field": name, "stat": "*", "golden": None, "run": None,
                     "abs_diff": float("inf"), "allowance": 0.0,
                     "rel_drift": float("inf"), "ok": False, "_margin": float("inf")}
            continue
        rvals = rdiag[name]
        # Field magnitude: the largest |extreme| across golden AND run, so a
        # spurious spin-up lifts the scale and the near-zero mean can't hide it.
        scale = max(abs(gvals["min"]), abs(gvals["max"]),
                    abs(rvals["min"]), abs(rvals["max"]))
        for stat in ("min", "max", "mean"):
            consider(name, stat, gvals[stat], rvals[stat], scale)

    gstats = golden_summary.get("stats", {})
    rstats = run_summary.get("stats", {})
    for key, gval in gstats.items():
        if key not in rstats:
            problems.append({
                "field": "stats", "stat": key, "golden": gval, "run": None,
                "abs_diff": float("inf"), "allowance": 0.0,
                "rel_drift": float("inf"), "ok": False,
                "note": "stat absent from run",
            })
            continue
        rval = rstats[key]
        scale = max(abs(gval), abs(rval))
        consider("stats", key, gval, rval, scale)

    if worst is not None:
        worst.pop("_margin", None)
    passed = len(problems) == 0
    return passed, worst, problems


# ---------------------------------------------------------------------------
# Case execution (reuses run_regression.run_case for the model run + P0 gate)
# ---------------------------------------------------------------------------
def _run_and_parse(case, binary, scratch_root, gpu_id, keep_scratch):
    """Run one case (P0 gate) and attach its parsed summary. Never raises."""
    res = rr.run_case(case, binary, scratch_root, gpu_id=gpu_id,
                      keep_scratch=keep_scratch, keep_output=True)
    output = res.pop("output", "")  # keep it out of the results JSON (large)
    res["summary"] = parse_summary(output)
    return res


def execute_cases(cases, binary, scratch_root, backend, gpus, keep_scratch):
    """Run all cases, returning {name: result-with-summary}. CPU serial;
    GPU farms one case per device via a work queue (same rule as the P0 GPU
    path -- a GPU is never shared)."""
    out = {}
    if backend == "cpu":
        for case in cases:
            out[case["name"]] = _run_and_parse(
                case, binary, scratch_root, None, keep_scratch)
        return out

    work = queue.Queue()
    for case in cases:
        work.put(case)
    lock = threading.Lock()

    def worker(gpu_id):
        while True:
            try:
                case = work.get_nowait()
            except queue.Empty:
                return
            try:
                res = _run_and_parse(case, binary, scratch_root, gpu_id,
                                     keep_scratch)
                with lock:
                    out[case["name"]] = res
            finally:
                work.task_done()

    threads = [threading.Thread(target=worker, args=(g,), daemon=True)
               for g in gpus]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return out


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _fmt(x):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isinf(x):
            return "inf"
        return "{:.4g}".format(x)
    return str(x)


def print_compare_report(rows):
    print("\n{:<28} {:<6} {:>8}  {:<22} {:>10}  {}".format(
        "case", "result", "wall_s", "worst field", "rel drift", "note"))
    print("-" * 96)
    for r in rows:
        worst = r.get("worst") or {}
        wf = "{}:{}".format(worst.get("field", "-"), worst.get("stat", "-"))
        print("{:<28} {:<6} {:>8}  {:<22} {:>10}  {}".format(
            r["case"],
            "PASS" if r["passed"] else "FAIL",
            "{:.2f}".format(r.get("wall_s", 0.0)),
            wf,
            _fmt(worst.get("rel_drift")),
            r.get("reason", ""),
        ))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="P1 golden-summary regression compare (high-tolerance).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--compare", action="store_true", default=True,
                      help="(default) run + compare to committed goldens.")
    mode.add_argument("--update-golden", action="store_true",
                      help="run + (re)write golden/<case>.json from the runs.")
    p.add_argument("--backend", choices=["cpu", "gpu"], default="cpu")
    p.add_argument("--build-dir", default="build_gcc",
                   help="Build dir with the 'rdb' binary (rel to repo root "
                        "unless absolute).")
    p.add_argument("--binary", default=None, help="Explicit rdb path.")
    p.add_argument("--gpus", default=None,
                   help="Comma-separated GPU ids for --backend gpu (e.g. "
                        "'0,1,2,3'). Default: auto-detect.")
    p.add_argument("--cases", default=None,
                   help="Comma-separated case-name filter (default: all).")
    p.add_argument("--rtol", type=float, default=DEFAULT_RTOL,
                   help="Relative tolerance (per-field magnitude). Manifest "
                        "'tol' key overrides per case.")
    p.add_argument("--atol", type=float, default=DEFAULT_ATOL,
                   help="Absolute floor (absorbs roundoff-zero fields).")
    p.add_argument("--out", default=None,
                   help="Write machine-readable JSON results here.")
    p.add_argument("--budget-min", type=float, default=rr.DEFAULT_BUDGET_MIN,
                   help="Per-backend wallclock budget (minutes); exceeding it "
                        "is a loud failure.")
    p.add_argument("--keep-scratch", action="store_true",
                   help="Keep per-case scratch dirs even on success.")
    p.add_argument("--scratch-root", default=None,
                   help="Root for scratch dirs (default: "
                        "<repo>/tmp_local_artifacts/regression).")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    update = args.update_golden

    build_dir = args.build_dir
    if not os.path.isabs(build_dir):
        build_dir = os.path.join(REPO_ROOT, build_dir)
    try:
        binary = rr.locate_binary(build_dir, args.binary)
    except FileNotFoundError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    scratch_root = args.scratch_root or os.path.join(
        REPO_ROOT, "tmp_local_artifacts", "regression")
    os.makedirs(scratch_root, exist_ok=True)

    cases = list(manifest.CASES)
    if args.cases:
        wanted = {c.strip() for c in args.cases.split(",") if c.strip()}
        cases = [c for c in cases if c["name"] in wanted]
        missing = wanted - {c["name"] for c in cases}
        if missing:
            print("ERROR: unknown case(s): {}".format(", ".join(sorted(missing))),
                  file=sys.stderr)
            return 2
    cases = [c for c in cases if args.backend in c.get("backends", {args.backend})]
    if not cases:
        print("ERROR: no cases for backend '{}'".format(args.backend),
              file=sys.stderr)
        return 2

    gpus = ["cpu"]
    if args.backend == "gpu":
        gpus = rr.detect_gpus(args.gpus)

    print("rdb binary : {}".format(binary))
    print("mode          : {}".format("UPDATE-GOLDEN" if update else "COMPARE"))
    print("backend       : {}".format(args.backend))
    print("cases         : {}".format(len(cases)))
    if args.backend == "gpu":
        print("gpus          : {}".format(",".join(gpus)))
    if not update:
        print("tolerance     : rtol={:g}  atol={:g}".format(args.rtol, args.atol))

    t0 = time.time()
    runs = execute_cases(cases, binary, scratch_root, args.backend,
                         gpus if args.backend == "gpu" else None,
                         args.keep_scratch)
    total_wall = time.time() - t0
    budget_blown = total_wall > args.budget_min * 60.0

    rows = []
    n_pass = 0
    for case in cases:
        name = case["name"]
        res = runs.get(name, {})
        row = {
            "case": name,
            "backend": args.backend,
            "wall_s": res.get("wall_s", 0.0),
            "gate_passed": res.get("passed", False),
            "gate_reason": res.get("reason", "did not run"),
            "passed": False,
            "reason": "",
            "worst": None,
        }
        summary = res.get("summary", {"diag": {}, "stats": {}})

        # The P0 run-clean/NaN gate must hold first.
        if not res.get("passed", False):
            row["reason"] = "gate: {}".format(res.get("reason", "did not run"))
            rows.append(row)
            continue

        if update:
            write_golden(case, summary, args.backend)
            row["passed"] = True
            row["reason"] = "golden written"
            n_pass += 1
            rows.append(row)
            continue

        golden = load_golden(name)
        if golden is None:
            row["reason"] = "NO GOLDEN (run --update-golden first)"
            rows.append(row)
            continue

        rtol = case.get("tol", args.rtol)
        ok, worst, problems = compare_summary(
            golden.get("summary", {}), summary, rtol, args.atol)
        row["worst"] = worst
        row["problems"] = problems
        row["rtol"] = rtol
        if ok:
            row["passed"] = True
            row["reason"] = "ok"
            n_pass += 1
        else:
            row["passed"] = False
            row["reason"] = "{} value(s) drifted".format(len(problems))
        rows.append(row)

    print_compare_report(rows)

    print("\n" + "=" * 96)
    verb = "written" if update else "passed"
    print("SUMMARY: {}/{} {}  |  wall {:.1f}s ({:.1f} min)  |  budget {:.0f} min".format(
        n_pass, len(cases), verb, total_wall, total_wall / 60.0, args.budget_min))
    if budget_blown:
        print("         !! BUDGET BLOWN -- failing loud.")

    # Detail the failing values so a drift is actionable.
    if not update:
        for r in rows:
            if not r["passed"] and r.get("problems"):
                print("\n  {} drifted values:".format(r["case"]))
                for pr in r["problems"][:12]:
                    print("    {:<14} {:<5} golden={:<14} run={:<14} "
                          "|d|={:<12} allow={:<12}".format(
                              pr["field"], pr["stat"],
                              _fmt(pr["golden"]), _fmt(pr["run"]),
                              _fmt(pr["abs_diff"]), _fmt(pr["allowance"])))

    if args.out:
        meta = {
            "mode": "update-golden" if update else "compare",
            "backend": args.backend, "binary": binary,
            "rtol": args.rtol, "atol": args.atol,
            "n_cases": len(cases), "n_pass": n_pass,
            "total_wall_s": round(total_wall, 2),
            "budget_min": args.budget_min, "budget_blown": budget_blown,
        }
        # Drop bulky per-value problem lists' redundancy but keep them.
        with open(args.out, "w") as fh:
            json.dump({"meta": meta, "results": rows}, fh, indent=2,
                      sort_keys=True)
        print("         JSON: {}".format(os.path.abspath(args.out)))

    failed = (n_pass != len(cases)) or budget_blown
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
