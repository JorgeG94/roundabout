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

Tolerance model
---------------
CPU (gfortran) and GPU (nvfortran) differ by more than bitwise -- FMA / reduction
order / transcendental libraries -- so each value is compared numpy-`isclose`
style against ITS FIELD'S tolerance:

    allowance = atol_field + rtol * scale
    pass      = |golden - run| <= allowance

where `scale` is the FIELD's magnitude (for a diag field: max |min|,|max| over
both golden and run; for a stats scalar: max |golden|,|run|). `rtol` is the
golden tolerance (default 1e-3, per-case overridable via a manifest `tol` key).

`atol_field` is ZERO unless the case's manifest entry declares it EXPLICITLY:
`"atol": {"SSH": 1e-9, "stats:En": 1e-18, ...}` plus a non-empty
`"atol_reason"`. There is deliberately NO global absolute floor. The old global
`atol = 1e-7` (SI units, applied to every field of every case) is how
`seamount` / `seamount_pred_corr` reported PASS while their SSH / u / KE / En
moved by 70-100 %: each whole field sat at 1e-11 or below, so ANY change up to
1e-7 -- a 10^4-fold growth of the rest-state noise included -- was inside the
floor. It also hid drifts on small-but-physical fields (ideal_age_demo SSH,
~8e-7 m, moved 9 %). A physically-zero field carrying roundoff noise (a resting
case's u / SSH / KE / En) must name its floor in the manifest, sized to the
field's own noise, so the pass is visible and reviewable. The verdict is
re-derived from every checked value (`_verdict_invariant`) and can never be
PASS while one of them exceeds its allowance. `--self-test` pins all of this
without running a model.

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
#   * Two early-transient baroclinic-instability cases (baroclinic_2layer, eady)
#     whose tiny cross-channel v / free-surface SSH EXTREMA diverge 9-15% at the
#     noise floor -- while their integrated En/Mass/Salt/Temp agree within rtol
#     -- carry a documented per-case `tol` in the manifest instead of inflating
#     this global bound (the plan's "flag it, don't hide it" rule).
#   * There is NO default absolute floor. A PHYSICALLY-ZERO field carrying
#     roundoff noise (a quiescent seamount's velocity extrema at ~1e-13, which
#     differ ~6% relatively CPU vs GPU) gets a per-field `atol` declared in its
#     manifest entry with an `atol_reason` -- see the module docstring for how
#     the old global atol=1e-7 passed 70-100% drifts without saying so.
DEFAULT_RTOL = 1e-3

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
def stats_key(key):
    """Manifest `atol` key for a `[stats]` scalar (a diag field uses its name)."""
    return "stats:" + key


def case_atol(case, golden_summary=None):
    """Validate and return a case's EXPLICIT per-field absolute tolerances.

    The manifest spells them `"atol": {field: value}` -- `field` is a `[diag]`
    field name (the floor applies to its min / max / mean) or `stats:<key>` --
    together with a non-empty `"atol_reason"`. Raises ValueError on a malformed
    entry or a missing reason, and -- when `golden_summary` is given -- on a
    key naming no field of the golden (a typo would otherwise check nothing).
    A case without `atol` gets `{}`: NO absolute floor on any field.
    """
    atol = case.get("atol")
    if atol is None:
        if case.get("atol_reason"):
            raise ValueError("{}: 'atol_reason' without 'atol'".format(
                case["name"]))
        return {}
    if not isinstance(atol, dict) or not atol:
        raise ValueError("{}: 'atol' must be a non-empty dict "
                         "{{field: value}}".format(case["name"]))
    reason = case.get("atol_reason", "")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("{}: 'atol' needs a non-empty 'atol_reason' saying "
                         "why each field is legitimately near zero".format(
                             case["name"]))
    for key, val in atol.items():
        if (isinstance(val, bool) or not isinstance(val, (int, float))
                or not math.isfinite(val) or val <= 0.0):
            raise ValueError("{}: atol[{!r}] = {!r} must be a finite number "
                             "> 0".format(case["name"], key, val))
    if golden_summary is not None:
        known = set(golden_summary.get("diag", {}))
        known |= {stats_key(k) for k in golden_summary.get("stats", {})}
        unknown = sorted(set(atol) - known)
        if unknown:
            raise ValueError("{}: 'atol' names field(s) absent from the "
                             "golden: {}".format(case["name"],
                                                 ", ".join(unknown)))
    return {k: float(v) for k, v in atol.items()}


def _isclose(g, r, rtol, atol, scale):
    """numpy-style closeness on one value against a field-magnitude `scale`.

    Returns (passed, abs_diff, allowance, rel_drift). `rel_drift` is |g-r| over
    the field scale (with a tiny floor) -- a scale-free "how far it moved"
    number for reporting; `passed` uses the atol+rtol*scale allowance, where
    `atol` is the field's EXPLICIT manifest floor (0.0 when none is declared).
    """
    if math.isnan(g) or math.isnan(r) or math.isinf(g) or math.isinf(r):
        return (False, float("inf"), 0.0, float("inf"))
    abs_diff = abs(g - r)
    allowance = atol + rtol * scale
    rel_drift = abs_diff / max(scale, 1e-30)
    return (abs_diff <= allowance, abs_diff, allowance, rel_drift)


def _verdict_invariant(passed, checked):
    """Refuse a PASS verdict while any checked value is outside its allowance.

    Re-derived from the raw golden/run numbers, not from the `ok` flags, so a
    future edit to the pass logic cannot re-open the silent-pass hole (the old
    global atol) without this assertion firing.
    """
    if not passed:
        return passed
    for rec in checked:
        g, r = rec["golden"], rec["run"]
        bad = (g is None or r is None or not math.isfinite(g)
               or not math.isfinite(r) or abs(g - r) > rec["allowance"])
        if bad:
            raise AssertionError(
                "compare verdict PASS with {}:{} outside its tolerance "
                "(golden={} run={} allowance={})".format(
                    rec["field"], rec["stat"], g, r, rec["allowance"]))
    return passed


def compare_summary(golden_summary, run_summary, rtol, atol=None):
    """Compare a run summary to a golden summary.

    `atol` is the case's explicit per-field floor map (see `case_atol`); a
    field it does not name gets NO absolute floor.

    Returns (passed, worst, problems, checked). `worst` is the value with the
    largest drift/allowance margin:
        {"field","stat","golden","run","abs_diff","allowance","rel_drift",
         "atol","ok"}
    `problems` lists every value that exceeded its allowance (incl. golden
    fields entirely absent from the run); `checked` lists every value compared.
    """
    atol = atol or {}
    problems = []
    checked = []
    worst = None

    def rank(rec, margin):
        nonlocal worst
        if worst is None or margin > worst["_margin"]:
            worst = dict(rec, _margin=margin)

    def consider(field, stat, g, r, scale, field_atol):
        ok, abs_diff, allowance, rel_drift = _isclose(g, r, rtol, field_atol,
                                                      scale)
        rec = {
            "field": field, "stat": stat, "golden": g, "run": r,
            "abs_diff": abs_diff, "allowance": allowance,
            "rel_drift": rel_drift, "atol": field_atol, "ok": ok,
        }
        checked.append(rec)
        # "worst" ranks by how far past the allowance we are (drift/allowance),
        # so the reported field is the closest to (or furthest past) failing.
        if allowance > 0:
            margin = abs_diff / allowance
        else:
            margin = 0.0 if abs_diff == 0 else float("inf")
        rank(rec, margin)
        if not ok:
            problems.append(rec)

    def absent(field, stat, g, note):
        rec = {"field": field, "stat": stat, "golden": g, "run": None,
               "abs_diff": float("inf"), "allowance": 0.0,
               "rel_drift": float("inf"), "atol": 0.0, "ok": False,
               "note": note}
        checked.append(rec)
        problems.append(rec)
        rank(rec, float("inf"))

    gdiag = golden_summary.get("diag", {})
    rdiag = run_summary.get("diag", {})
    for name, gvals in gdiag.items():
        if name not in rdiag:
            absent(name, "*", None, "field absent from run")
            continue
        rvals = rdiag[name]
        # Field magnitude: the largest |extreme| across golden AND run, so a
        # spurious spin-up lifts the scale and the near-zero mean can't hide it.
        scale = max(abs(gvals["min"]), abs(gvals["max"]),
                    abs(rvals["min"]), abs(rvals["max"]))
        for stat in ("min", "max", "mean"):
            consider(name, stat, gvals[stat], rvals[stat], scale,
                     atol.get(name, 0.0))

    gstats = golden_summary.get("stats", {})
    rstats = run_summary.get("stats", {})
    for key, gval in gstats.items():
        if key not in rstats:
            absent("stats", key, gval, "stat absent from run")
            continue
        rval = rstats[key]
        scale = max(abs(gval), abs(rval))
        consider("stats", key, gval, rval, scale,
                 atol.get(stats_key(key), 0.0))

    if worst is not None:
        worst.pop("_margin", None)
    passed = _verdict_invariant(len(problems) == 0, checked)
    return passed, worst, problems, checked



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


def _margin(worst):
    """|golden - run| / allowance of a record: <= 1 is inside its tolerance."""
    if not worst or worst.get("abs_diff") is None:
        return None
    allowance = worst.get("allowance") or 0.0
    if allowance > 0:
        return worst["abs_diff"] / allowance
    return 0.0 if worst["abs_diff"] == 0 else float("inf")


def print_compare_report(rows):
    # `|d|/allow` is the verdict quantity (a PASS row is always <= 1);
    # `rel drift` is |d|/scale, reported so a large relative move that an
    # EXPLICIT manifest atol absorbs is still visible next to its PASS.
    print("\n{:<28} {:<6} {:>8}  {:<22} {:>10} {:>10}  {}".format(
        "case", "result", "wall_s", "worst field", "|d|/allow", "rel drift",
        "note"))
    print("-" * 108)
    for r in rows:
        worst = r.get("worst") or {}
        wf = "{}:{}".format(worst.get("field", "-"), worst.get("stat", "-"))
        print("{:<28} {:<6} {:>8}  {:<22} {:>10} {:>10}  {}".format(
            r["case"],
            "PASS" if r["passed"] else "FAIL",
            "{:.2f}".format(r.get("wall_s", 0.0)),
            wf,
            _fmt(_margin(worst)),
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
    p.add_argument("--self-test", action="store_true",
                   help="verify the compare verdict itself against synthetic "
                        "and measured summaries (incl. the seamount silent "
                        "pass), and validate every manifest 'atol' against "
                        "its golden, then exit. Runs no model.")
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


def self_test():
    """Who tests the test: pin the compare VERDICT without running a model.

    Every summary pair below is MEASURED (gfortran CPU, 10 steps) unless
    marked synthetic. Case (1) is the silent pass this harness shipped with:
    under the old global atol=1e-7 it returned PASS; it must now FAIL.
    """
    fails = []

    def check(label, ok, detail=""):
        print("  {:<64} {}".format(label, "ok" if ok else "FAILED " + detail))
        if not ok:
            fails.append(label)

    def diag(**fields):
        return {k: {"min": v[0], "max": v[1], "mean": v[2]}
                for k, v in fields.items()}

    # seamount golden (committed) vs the integration-branch run (bebt = 0.1).
    sm_gold = {"diag": diag(SSH=(-8.18545e-12, 1.04592e-11, -8.61618e-14),
                            u=(-3.80073e-13, 3.80073e-13, -5.03544e-32),
                            KE=(0.0, 2.88911e-25, 1.16445e-26)),
               "stats": {"En": 3.673e-27, "Mass": 3.37493e+16}}
    sm_run = {"diag": diag(SSH=(-3.63798e-12, 3.18323e-12, -2.22227e-13),
                           u=(-1.96474e-13, 1.96474e-13, -5.0e-32),
                           KE=(0.0, 2.39667e-26, 9.38137e-28)),
              "stats": {"En": 3.484e-28, "Mass": 3.37493e+16}}
    sm_atol = {"SSH": 1e-9, "u": 1e-9, "KE": 1e-18, "stats:En": 1e-18}

    # (1) The silent pass. The OLD rule (a global floor of 1e-7 on every
    #     value) passed it; the NEW rule, with no manifest atol, fails it.
    def old_rule(g, r, old_atol=1e-7):
        _, _, _, chk = compare_summary(g, r, DEFAULT_RTOL)
        return all(c["abs_diff"] <= old_atol + c["allowance"] for c in chk)
    check("(1a) old global atol=1e-7 PASSED the seamount 70% SSH drift",
          old_rule(sm_gold, sm_run))
    ok, worst, probs, _ = compare_summary(sm_gold, sm_run, DEFAULT_RTOL)
    check("(1b) ...and with no manifest atol it now FAILS",
          not ok and worst["field"] in ("SSH", "u", "KE", "stats"))
    check("(1c) ...naming SSH:max among the drifted values",
          any(p["field"] == "SSH" and p["stat"] == "max" for p in probs))

    # (2) The same pair with an EXPLICIT per-field floor passes -- visibly.
    ok, _, probs, chk = compare_summary(sm_gold, sm_run, DEFAULT_RTOL,
                                        sm_atol)
    check("(2a) explicit manifest atol (SSH/u 1e-9, KE/En 1e-18) -> PASS",
          ok and not probs)
    check("(2b) ...and the values it absorbed are identifiable",
          any(c["atol"] > 0 and c["abs_diff"] > c["allowance"] - c["atol"]
              for c in chk))
    # ...but a 1000x spin-up of the noise (SSH 1e-8 m) still fails it.
    spun = json.loads(json.dumps(sm_run))
    spun["diag"]["SSH"]["max"] = 1.0e-8
    ok, _, _, _ = compare_summary(sm_gold, spun, DEFAULT_RTOL, sm_atol)
    check("(2c) a 1000x SSH spin-up (1e-8 m) FAILS even with the atol",
          not ok)
    # ...and a floor on SSH does not leak onto an undeclared field.
    ok, _, probs, _ = compare_summary(sm_gold, sm_run, DEFAULT_RTOL,
                                      {"SSH": 1e-9})
    check("(2d) atol on SSH only -> u/KE/En still FAIL",
          not ok and {p["field"] for p in probs} >= {"u", "KE", "stats"})

    # (3) A small-but-physical field: ideal_age_demo SSH:max 8.3409e-07 ->
    #     7.5544e-07 (9.4%) was under the old floor; it must fail now.
    g = {"diag": diag(SSH=(-3.67178e-07, 8.34093e-07, 1.0e-9)), "stats": {}}
    r = {"diag": diag(SSH=(-3.49592e-07, 7.55444e-07, 1.0e-9)), "stats": {}}
    check("(3) ideal_age_demo SSH 9.4% drift at 8e-7 m -> FAIL",
          not compare_summary(g, r, DEFAULT_RTOL)[0])

    # (4) rtol on an active field (synthetic): 5e-4 passes, 2e-3 fails.
    g = {"diag": diag(u=(-0.1, 0.1, 0.0)), "stats": {"En": 1.0e-3}}
    for rel, want in ((5e-4, True), (2e-3, False)):
        r = {"diag": diag(u=(-0.1, 0.1 * (1 + rel), 0.0)),
             "stats": {"En": 1.0e-3}}
        check("(4) active-field drift {:g} -> {}".format(
            rel, "PASS" if want else "FAIL"),
            compare_summary(g, r, DEFAULT_RTOL)[0] is want)

    # (5) exact zeros agree (island_at_rest shape); NaN / absent always fail.
    z = {"diag": diag(u=(0.0, 0.0, 0.0)), "stats": {"En": 0.0}}
    check("(5a) exact 0.0 == 0.0 with no atol -> PASS",
          compare_summary(z, z, DEFAULT_RTOL)[0])
    nz = {"diag": diag(u=(0.0, 1e-20, 0.0)), "stats": {"En": 0.0}}
    check("(5b) 0.0 vs 1e-20 with no atol -> FAIL (no hidden floor)",
          not compare_summary(z, nz, DEFAULT_RTOL)[0])
    nan = {"diag": diag(u=(0.0, float("nan"), 0.0)), "stats": {"En": 0.0}}
    check("(5c) NaN in the run -> FAIL even under a huge atol",
          not compare_summary(z, nan, DEFAULT_RTOL, {"u": 1e30})[0])
    check("(5d) golden field absent from the run -> FAIL",
          not compare_summary(z, {"diag": {}, "stats": {"En": 0.0}},
                              DEFAULT_RTOL)[0])
    check("(5e) golden stat absent from the run -> FAIL",
          not compare_summary(z, {"diag": z["diag"], "stats": {}},
                              DEFAULT_RTOL)[0])

    # (6) The verdict invariant refuses a PASS over an out-of-tolerance value.
    try:
        _verdict_invariant(True, [{"field": "SSH", "stat": "max",
                                   "golden": 1.0, "run": 2.0,
                                   "allowance": 0.5}])
        raised = False
    except AssertionError:
        raised = True
    check("(6) _verdict_invariant raises on PASS with |d| > allowance",
          raised)

    # (7) Manifest atol validation.
    bad = [({"name": "x", "atol": {"SSH": 1e-9}}, "no atol_reason"),
           ({"name": "x", "atol": {"SSH": 0.0}, "atol_reason": "r"},
            "zero atol"),
           ({"name": "x", "atol": {"SSH": float("nan")}, "atol_reason": "r"},
            "NaN atol"),
           ({"name": "x", "atol": {}, "atol_reason": "r"}, "empty atol"),
           ({"name": "x", "atol_reason": "r"}, "reason without atol"),
           ({"name": "x", "atol": {"SHH": 1e-9}, "atol_reason": "r"},
            "atol key naming no golden field")]
    for case, what in bad:
        try:
            case_atol(case, sm_gold)
            refused = False
        except ValueError:
            refused = True
        check("(7) case_atol refuses: {}".format(what), refused)

    # (8) The REAL manifest: no global floor survives, and every declared
    #     atol is well-formed and names fields of its committed golden.
    check("(8a) no DEFAULT_ATOL global floor in this module",
          "DEFAULT_ATOL" not in globals())
    for case in manifest.CASES:
        if "atol" not in case and "atol_reason" not in case:
            continue
        golden = load_golden(case["name"])
        try:
            case_atol(case, golden.get("summary", {}) if golden else None)
            valid = golden is not None
        except ValueError as exc:
            valid = False
            print("    {}".format(exc))
        check("(8b) manifest atol valid vs golden: {}".format(case["name"]),
              valid)

    print("\ncompare.py self-test: {}".format(
        "ALL PASS" if not fails else "{} FAILED".format(len(fails))))
    return 1 if fails else 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.self_test:
        print("compare.py self-test -- no model is run\n")
        return self_test()
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
        print("tolerance     : rtol={:g} (manifest 'tol' per case); NO global "
              "atol -- only a case's explicit manifest 'atol'".format(args.rtol))

    t0 = time.time()
    runs = execute_cases(cases, binary, scratch_root, args.backend,
                         gpus if args.backend == "gpu" else None,
                         args.keep_scratch)
    total_wall = time.time() - t0
    budget_blown = total_wall > args.budget_min * 60.0

    rows = []
    n_pass = 0
    skipped_on_update = []
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
            if update:
                # A gate-failing case gets NO new golden -- its old file stays
                # on disk.  Say so by name: a bare "FAIL gate:" row in a table
                # of 30 "golden written" rows is how five goldens silently
                # survived two default changes (see README, "Skipped on
                # update").
                row["reason"] += "  [golden NOT written]"
                skipped_on_update.append(name)
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
        try:
            atol = case_atol(case, golden.get("summary", {}))
        except ValueError as exc:
            row["reason"] = "manifest: {}".format(exc)
            rows.append(row)
            continue
        ok, worst, problems, checked = compare_summary(
            golden.get("summary", {}), summary, rtol, atol)
        row["worst"] = worst
        row["problems"] = problems
        row["rtol"] = rtol
        row["atol"] = atol
        if ok:
            row["passed"] = True
            # A value that is inside its allowance ONLY because of the case's
            # explicit manifest atol is named, never folded into a bare "ok".
            # (allowance - atol) is the relative part, rtol * scale.
            on_atol = sorted("{}:{}".format(c["field"], c["stat"])
                             for c in checked
                             if c["atol"] > 0
                             and c["abs_diff"] > c["allowance"] - c["atol"])
            row["passed_on_atol"] = on_atol
            row["reason"] = ("ok" if not on_atol else
                             "ok; {} value(s) inside manifest atol only"
                             .format(len(on_atol)))
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
    if skipped_on_update:
        print("         !! {} golden(s) NOT written -- the case failed the "
              "run-clean gate, so its".format(len(skipped_on_update)))
        print("            OLD golden is still on disk and is now STALE "
              "relative to this update:")
        for nm in skipped_on_update:
            print("              - {}".format(nm))
        print("            Fix the gate failure and re-run --update-golden "
              "for these cases. Exiting non-zero.")

    # Detail the failing values so a drift is actionable, and name every value
    # that passed only on a manifest atol so such a pass is never silent.
    if not update:
        for r in rows:
            if r["passed"] and r.get("passed_on_atol"):
                print("\n  {} passed on its manifest atol only: {}".format(
                    r["case"], ", ".join(r["passed_on_atol"])))
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
            "rtol": args.rtol, "atol": "per-case manifest only",
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
