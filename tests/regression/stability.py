#!/usr/bin/env python3
"""Two-tier STABILITY-CHECKED validation suite for the shipped ocean namelists.

What this adds that the golden-drift suite cannot do
====================================================
`compare.py` runs every case for 6-10 outer steps and compares the final state
to a committed golden. That has two structural blind spots, and on 2026-09-11
both of them shipped real defects:

  1. **10 steps is 1.7 simulated hours.** The bugs found that day manifest at
     hour 3, hour 7, step 131 and "days". The suite stopped before the physics
     existed.
  2. **A golden enshrines whatever the run did.** A golden captured from a
     broken run passes forever. `eady` passed its golden for months BECAUSE
     the golden was the over-damped answer (fixed 2026-09-13).

This module asserts on the model's own PHYSICS instead, from its own console
time series -- no golden, nothing to enshrine:

  * **finite**       -- no NaN/Inf in any conserved scalar or field extremum.
  * **conservation** -- the model's own closed-budget residuals (Mass / Salt /
                        Heat `Error`) stay at roundoff.
  * **stability**    -- in an unforced or adiabatic case KINETIC ENERGY MUST
                        NOT GROW. This is the gate that catches what NaN and
                        budget checks miss: on 2026-09-11 a thermo-cadence
                        instability grew `En` 30-60x while Mass, Salt and Temp
                        stayed exact to every printed digit and nothing
                        crashed.
  * **bounded CFL**  -- a monotone climb to a trip is numerical instability,
                        not physics.
  * **the case's own stated expectation**, where its header makes a testable
    one. Several of those turn out to be FALSE; that is reported as a case
    failure, never patched away by editing the namelist.

The two tiers
=============
TIER 1 -- full scale, GPU (NVHPC), LOCAL / nightly.
    Each case runs long enough that its physics actually manifests: a
    per-case duration, not a global step count. Tens of seconds to minutes
    per case. This is the tier that can see a baroclinic growth rate or a
    multi-day saturation. Not wired into GitHub Actions: it needs a GPU.

TIER 2 -- downscaled twins, gfortran CPU, the CI GATE.
    Sized for a standard GitHub-hosted runner (no NVIDIA hardware, 2-4 cores,
    ~7 GB RAM). Each twin is seconds. Twins are built by preserving the
    DIMENSIONLESS numbers, never by naive grid shrinking -- see
    `downscale.py`, whose rules are CHECKED for every twin before it runs.
    Tier 2 keeps the finite / conservation / energy-growth / CFL gates; it
    drops only the assertions that genuinely need a long integration (a
    baroclinic growth rate cannot be measured in 40 simulated hours), and
    says so rather than weakening them.

Usage
-----
    python3 tests/regression/stability.py --tier 2 --build-dir build_gcc
    python3 tests/regression/stability.py --tier 1 --backend gpu \
            --build-dir build --gpus 0,1,2,3
    python3 tests/regression/stability.py --tier 2 --cases eady,acc_channel -v

Stdlib only (the repo forbids `pip install`).
"""

import argparse
import json
import math
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
sys.path.insert(0, THIS_DIR)
import downscale  # noqa: E402
import manifest  # noqa: E402
import stability_manifest as _stab_manifest  # noqa: E402

# The scheme a BASE case runs, i.e. the `&ocean_bt_nml split_scheme` default in
# `rdb_config.F90`.  Printed in the run header so a result file can never be
# read as "the bare name means whatever is default today": the axis suffix
# names the OTHER scheme, and these two must be different.
_DEFAULT_SPLIT_SCHEME = "pred_corr"
import run_regression as rr  # noqa: E402  (binary discovery + bad-marker list)

# ---------------------------------------------------------------------------
# Console parsing
# ---------------------------------------------------------------------------
# A Fortran numeric token, including the IEEE spellings we must be able to READ
# (so a NaN becomes a parsed float('nan') that an assertion fails on, rather
# than an unparsed line the suite silently skips).
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?|[-+]?NaN|[-+]?Infinity|[-+]?Inf"

_STATS_RE = re.compile(
    r"\[stats\]\s+Day\s+(?P<day>" + _NUM + r")\s+step\s+(?P<step>\d+)"
    r"\s+En\s+(?P<en>" + _NUM + r")"
    r"(?:\s+MaxCFL\s+(?P<cfl>" + _NUM + r"))?"
    r"(?:\s+Mass\s+(?P<mass>" + _NUM + r"))?"
    r"(?:\s+Salt\s+(?P<salt>" + _NUM + r"))?"
    r"(?:\s+Temp\s+(?P<temp>" + _NUM + r"))?"
)
# "    Mass :  6.468750000E+16  Error -3.215E-15  out  3.034E-08  src  0.000E+00"
_BUDGET_RE = re.compile(
    r"^\s+(?P<what>Mass|Salt|Heat|En)\s+:\s+(?P<total>" + _NUM + r")"
    r"(?:\s+Error\s+(?P<err>" + _NUM + r"))?"
)
# "[diag] t=  600.00 temperature [degC]  min= ..  max= ..  mean= .."
_DIAG_RE = re.compile(
    r"\[diag\]\s+t=\s*(?P<t>" + _NUM + r")\s+(?P<name>.+?)\s+\[.*?\]\s+"
    r"min=\s*(?P<min>" + _NUM + r")\s+max=\s*(?P<max>" + _NUM + r")"
    r"\s+mean=\s*(?P<mean>" + _NUM + r")"
)

# Crash markers. Deliberately NARROWER than run_regression.BAD_MARKERS in two
# ways, both learned the hard way:
#
#   * no bare `\bnan\b`. That scan cannot tell a genuine blow-up from the
#     MISSING-DATA SENTINEL `fill_tracer_impl` writes into vanished layers by
#     design (see DIAG_MEAN_NAN_NOTE), and it fails the flagship
#     `double_gyre_mom6` on sight. This suite parses the NUMBERS and decides
#     per field instead.
#   * no bare `\babort\b`. `sea_ice_pack` prints a configuration warning
#     containing the word ("transport's positivity abort is the only
#     ice-velocity guard") and was reported as having crashed while exiting 0.
#
# What is left only matches text a healthy run cannot produce. A non-zero exit
# code is checked separately and is the primary signal.
_CRASH_MARKERS = tuple(re.compile(p, re.IGNORECASE) for p in (
    r"error stop", r"\baborted\b", r"\baborting\b", r"program abort",
    r"segmentation fault", r"floating point exception", r"backtrace",
))
# The solver's own in-flight repair mechanism. Its presence is not a crash --
# the run continues -- but it means non-finite values WERE produced and zeroed,
# which is a finding in its own right, so it is counted and reported.
_NAN_CATCH_RE = re.compile(r"\[nan-catch\]", re.IGNORECASE)
# "CFL truncations:  3 this report (17 total)" -- the solver's own count of
# faces whose velocity it had to truncate to hold the barotropic CFL. On a
# REST case any non-zero value means the run was not resting; it is a
# counter, not a controller (`dyn%ntrunc_total`), so nothing else reports it.
_NTRUNC_RE = re.compile(r"CFL truncations:\s+\d+\s+this report\s+\((\d+)\s+total\)")
# The configure-time terrain-following stiffness audit, which prints the
# number the rx0 ladder is built to dial:
#   "sigma stiffness rx0 = 0.726 at ... exceeds 0.200"
_RX0_RE = re.compile(r"rx0\s*=\s*(" + _NUM + r")")
# The model's authoritative end-of-run step count, independent of the [stats]
# emission cadence (which lands on a stride and rarely on the very last step).
_TOTAL_STEPS_RE = re.compile(r"^\s*Total steps:\s*(\d+)", re.MULTILINE)


def _f(tok):
    """Parse a Fortran numeric token to float; NaN/Inf spellings included."""
    if tok is None:
        return None
    t = tok.strip().replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        low = t.lower().lstrip("+-")
        if low.startswith("nan"):
            return float("nan")
        if low.startswith("inf"):
            return float("-inf") if t.strip().startswith("-") else float("inf")
        raise


def _finite(x):
    return x is not None and not math.isnan(x) and not math.isinf(x)


def parse_series(text):
    """Parse a run's stdout into time series.

    Returns:
        {"stats":   [{day, step, En, MaxCFL, Mass, Salt, Temp}, ...],
         "budget":  {"Mass": [err,...], "Salt": [...], "Heat": [...]},
         "diag":    {field: [{t, min, max, mean}, ...]},
         "crash":   first crash marker matched, or None}
    """
    stats, diag, budget = [], {}, {}
    for line in text.splitlines():
        m = _STATS_RE.search(line)
        if m:
            stats.append({
                "day": _f(m.group("day")), "step": int(m.group("step")),
                "En": _f(m.group("en")), "MaxCFL": _f(m.group("cfl")),
                "Mass": _f(m.group("mass")), "Salt": _f(m.group("salt")),
                "Temp": _f(m.group("temp")),
            })
            continue
        m = _BUDGET_RE.match(line)
        if m and m.group("err") is not None:
            budget.setdefault(m.group("what"), []).append(_f(m.group("err")))
            continue
        m = _DIAG_RE.search(line)
        if m:
            diag.setdefault(m.group("name").strip(), []).append({
                "t": _f(m.group("t")), "min": _f(m.group("min")),
                "max": _f(m.group("max")), "mean": _f(m.group("mean")),
            })

    crash = None
    for rx in _CRASH_MARKERS:
        mm = rx.search(text)
        if mm:
            crash = mm.group(0)
            break
    mt = None
    for mm in _TOTAL_STEPS_RE.finditer(text):
        mt = int(mm.group(1))
    ntrunc = 0
    for mm in _NTRUNC_RE.finditer(text):
        ntrunc = max(ntrunc, int(mm.group(1)))
    rx0 = None
    mm = _RX0_RE.search(text)
    if mm:
        rx0 = _f(mm.group(1))
    return {"stats": stats, "budget": budget, "diag": diag, "crash": crash,
            "total_steps": mt, "ntrunc": ntrunc, "rx0": rx0,
            "nan_catch": len(_NAN_CATCH_RE.findall(text))}


# ---------------------------------------------------------------------------
# Known, deliberate artifacts that must NOT be read as failures
# ---------------------------------------------------------------------------
# `fill_tracer_impl` (src/core/ocean/diag/rdb_ocean_diag_fills.F90) writes an
# IEEE quiet NaN into any cell whose layer has vanished (h <= H_VANISHED) --
# deliberately, so a pinched-out or below-bottom cell reads as MISSING rather
# than as a plausible 0 degC / 0 PSU. Land, halo and vanished bed cells are all
# such cells, so EVERY tracer-concentration diagnostic carries NaN, and the
# whole-array `mean` printed on the `[diag]` line is therefore NaN for any case
# with land, a halo or a vanishing coordinate.
#
# Two consequences this suite has to handle honestly:
#   * The printed `min`/`max` are NOT NaN, because `minval`/`maxval` (and the
#     GPU `reduce(min:)/reduce(max:)`) are NaN-blind -- they launder it. So the
#     extrema stay meaningful and ARE gated; only the mean is not.
#   * `run_regression.py`'s bare-text `\bnan\b` scan cannot make that
#     distinction and fails the flagship `double_gyre_mom6` case on sight.
# The mean being unusable is itself a reportable defect (the reduction should
# skip the sentinel), tracked as DIAG-MEAN-NAN in the suite report -- but it is
# a diagnostic-reporting defect, not a solver blow-up, so it is surfaced as an
# ARTIFACT line rather than silently ignored OR conflated with a real NaN.
DIAG_MEAN_NAN_NOTE = (
    "tracer-concentration diag `mean` is NaN because vanished/land/halo cells "
    "carry the deliberate missing-data sentinel from fill_tracer_impl and the "
    "mean reduction does not skip it (min/max are NaN-blind so they survive). "
    "Reported as an artifact, not a blow-up; the extrema ARE gated."
)


# ---------------------------------------------------------------------------
# Assertions
# ---------------------------------------------------------------------------
class Verdict(object):
    """One assertion's verdict, with everything a reader needs to act on it."""

    __slots__ = ("name", "ok", "expected", "observed", "meaning", "ref", "skipped")

    def __init__(self, name, ok, expected, observed, meaning="", ref="",
                 skipped=False):
        self.name = name
        self.ok = ok
        self.expected = expected
        self.observed = observed
        self.meaning = meaning
        self.ref = ref
        self.skipped = skipped

    def render(self, case, tier):
        head = "{} {} [tier{}] {}".format(
            "SKIP" if self.skipped else ("ok  " if self.ok else "FAIL"),
            case, tier, self.name)
        if self.ok and not self.skipped:
            return head + "   " + self.observed
        body = [head,
                "     expected : " + self.expected,
                "     observed : " + self.observed]
        if self.meaning:
            for i, chunk in enumerate(_wrap(self.meaning, 84)):
                body.append(("     meaning  : " if i == 0 else "                ")
                            + chunk)
        if self.ref:
            body.append("     see      : " + self.ref)
        return "\n".join(body)


def _wrap(text, width):
    words, line, out = text.split(), "", []
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w) if line else w
    if line:
        out.append(line)
    return out


def _fmt(x):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isnan(x):
            return "NaN"
        if math.isinf(x):
            return "Inf"
        return "{:.4g}".format(x)
    return str(x)


def assert_ran(series, want_steps, res):
    """The run must have reached the requested step count without aborting."""
    steps = series["stats"]
    last = steps[-1]["step"] if steps else -1
    if res.get("timed_out"):
        return Verdict(
            "completed", False,
            "reach step {} within {} s".format(want_steps, res["timeout_s"]),
            "TIMED OUT at step {} after {:.1f} s".format(last, res["wall_s"]),
            "The case did not finish inside its wallclock cap. Either it hangs, "
            "or its per-step cost is far above what the manifest budgeted.")
    if series["crash"]:
        return Verdict(
            "completed", False,
            "reach step {} and exit 0".format(want_steps),
            "ABORTED at step {} ('{}', rc={})".format(
                last, series["crash"], res.get("returncode")),
            "The solver aborted. The last [stats] line above the abort is where "
            "to start.")
    if res.get("returncode") not in (0, None):
        return Verdict(
            "completed", False, "exit 0",
            "exit rc={} at step {}".format(res.get("returncode"), last))
    if not steps:
        return Verdict("completed", False,
                       "at least one [stats] line",
                       "no [stats] output parsed -- the run produced nothing to "
                       "assert on",
                       "Usually a namelist rejected at configure time; read the "
                       "scratch run.log.")
    # The model prints its own authoritative "Total steps: N" at the end. Use
    # that, NOT the last [stats] line: [stats] fires on a cadence stride and so
    # normally lands a few steps short of the end. Reading the cadence as an
    # early exit reported nine healthy cases as failures on the first sweep.
    done = series.get("total_steps")
    if done is None:
        done = last
    if done < want_steps:
        return Verdict(
            "completed", False, "reach step {}".format(want_steps),
            "stopped at step {} (last [stats] sample at step {})".format(
                done, last),
            "The run ended early without an abort marker.")
    return Verdict("completed", True, "reach step {}".format(want_steps),
                   "ran {} steps ({:.2f} simulated days) in {:.1f} s".format(
                       done, steps[-1]["day"] or 0.0, res["wall_s"]))


def assert_finite(series):
    """No NaN/Inf in any conserved scalar, or in any field EXTREMUM.

    Field `mean` is excluded on purpose -- see DIAG_MEAN_NAN_NOTE. Returns
    (verdict, artifacts) where `artifacts` lists the deliberate-sentinel means
    so they are reported rather than hidden.
    """
    bad = []
    for s in series["stats"]:
        for k in ("En", "MaxCFL", "Mass", "Salt", "Temp"):
            v = s.get(k)
            if v is not None and not _finite(v):
                bad.append("{} = {} at step {}".format(k, _fmt(v), s["step"]))
                break
        if bad:
            break
    artifacts = []
    for name, recs in sorted(series["diag"].items()):
        for r in recs:
            for stat in ("min", "max"):
                if not _finite(r[stat]):
                    bad.append("diag {}:{} = {} at t={}".format(
                        name, stat, _fmt(r[stat]), _fmt(r["t"])))
                    break
            if not _finite(r["mean"]) and name not in artifacts:
                artifacts.append(name)
    if bad:
        return Verdict(
            "finite", False,
            "every conserved scalar and every field extremum stays finite",
            "; ".join(bad[:3]),
            "A NaN or Inf reached a prognostic. This is a blow-up, not a "
            "missing-data sentinel -- the sentinel only ever appears in a "
            "tracer-concentration MEAN."), artifacts
    return Verdict("finite", True,
                   "no NaN/Inf in En/MaxCFL/Mass/Salt/Temp or any field extremum",
                   "clean over {} samples".format(len(series["stats"]))), artifacts


def assert_conservation(series, budget_tol):
    """The model's own closed-budget residuals stay at roundoff.

    `budget_tol` maps a budget name (Mass/Salt/Heat) to its relative tolerance;
    a budget absent from the run is skipped, never assumed good.
    """
    out = []
    for what, tol in sorted(budget_tol.items()):
        errs = series["budget"].get(what)
        if not errs:
            out.append(Verdict(
                "conserve:" + what, True,
                "|Error| <= {:.1e}".format(tol),
                "not reported by this case (no {} budget)".format(what),
                skipped=True))
            continue
        worst = max(errs, key=lambda e: (0.0 if e is None or math.isnan(e)
                                         else abs(e)))
        worst = worst if worst is not None else 0.0
        ok = _finite(worst) and abs(worst) <= tol
        out.append(Verdict(
            "conserve:" + what, ok,
            "the model's own closed {} budget residual |Error| <= {:.1e} at "
            "every step".format(what, tol),
            "worst |Error| = {} over {} steps".format(_fmt(abs(worst) if _finite(worst) else worst), len(errs)),
            "The {} budget is not closing. Unlike an energy blow-up this is a "
            "LEAK: a kernel is creating or destroying {} rather than "
            "transporting it.".format(what, what.lower())))
    return out


def _early_reference(en, frac=0.10, min_n=3):
    """Reference energy for a growth ratio: the max over the first `frac` of
    the run (never the single t=0 sample, which is often exactly 0 for a
    from-rest case and would make every ratio infinite)."""
    n = max(min_n, int(math.ceil(len(en) * frac)))
    head = [e for e in en[:n] if _finite(e)]
    return max(head) if head else None


def _max_window_rate(days, en, window_frac=0.2, min_pts=4):
    """The FASTEST growth rate of log(En) over any window of the run, 1/s.

    Why a sliding window and not one fit over the whole run: the 2026-09-11
    thermo-cadence instability took `eady` from En 2.68e-4 to 8.59e-3 in TWO
    DAYS and then sat near that level for the remaining 23. A least-squares
    fit over the whole 25 days returns sigma = 8.6e-7 -- squarely inside the
    band the physical Eady mode occupies, so a whole-run fit CANNOT tell the
    instability from the physics it is impersonating. The two-day burst runs
    at sigma ~ 1e-5, twenty times outside it, and a windowed maximum sees that.

    Returns the amplitude rate (En ~ exp(2*sigma*t)), or None.
    """
    pts = [(d * 86400.0, math.log(e)) for d, e in zip(days, en)
           if _finite(d) and _finite(e) and e > 0.0]
    if len(pts) < min_pts:
        return None
    w = max(min_pts, int(math.ceil(len(pts) * window_frac)))
    best = None
    for i in range(0, len(pts) - w + 1):
        seg = pts[i:i + w]
        dt_span = seg[-1][0] - seg[0][0]
        if dt_span <= 0.0:
            continue
        r = (seg[-1][1] - seg[0][1]) / dt_span
        if best is None or r > best:
            best = r
    return None if best is None else best / 2.0


def _efold_rate(days, en):
    """Least-squares growth rate sigma of log(En) vs time, 1/s.

    En is a quadratic energy, so the AMPLITUDE growth rate is sigma/2. Returns
    None when there is not enough positive signal to fit.
    """
    pts = [(d * 86400.0, math.log(e)) for d, e in zip(days, en)
           if _finite(d) and _finite(e) and e > 0.0]
    if len(pts) < 4:
        return None
    n = float(len(pts))
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    if den == 0.0:
        return None
    return (n * sxy - sx * sy) / den


def fit_log_rate(days, en, from_frac=0.4):
    """Least-squares exponential rate of `En` over the TAIL of a run.

    Returns `(sigma_amp, r2, n, (day_first, day_last))` where `sigma_amp` is
    the AMPLITUDE rate in 1/s (`En ~ exp(2 sigma t)`), or
    `(None, None, n, None)` when there is not enough positive signal.

    Why the tail and not the whole run.  A rest case's `En` has two regimes:
    an initial adjustment (the spurious pressure gradient spinning up its
    geostrophic balance, over a few inertial periods) and then either a
    PLATEAU or an EXPONENTIAL.  Only the second is the thing being measured,
    and a whole-run fit mixes the two -- on the measured cavity case a
    whole-run fit over days 0-60 reads 0.13/day against 0.333/day fitted
    over the clean window, i.e. it understates the instability by 2.5x.
    `from_frac = 0.4` drops the leading 40% of samples.

    Why R^2 is returned and reported rather than gated.  A PLATEAU has a
    near-zero slope and an essentially meaningless R^2; an INSTABILITY has a
    slope and an R^2 near 1.  Reporting both lets a reader tell "flat, so the
    rate is noise" from "exponential, and here is the rate" without the suite
    having to pick a threshold for a quantity whose distribution under the
    null is not known.  The GATE is on sigma alone.
    """
    pts = [(d * 86400.0, math.log(e)) for d, e in zip(days, en)
           if _finite(d) and _finite(e) and e > 0.0]
    if len(pts) < 6:
        return None, None, len(pts), None
    pts = pts[int(len(pts) * from_frac):]
    if len(pts) < 4:
        return None, None, len(pts), None
    n = float(len(pts))
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    if den == 0.0:
        return None, None, len(pts), None
    slope = (n * sxy - sx * sy) / den
    inter = (sy - slope * sx) / n
    ybar = sy / n
    ss_tot = sum((p[1] - ybar) ** 2 for p in pts)
    ss_res = sum((p[1] - (slope * p[0] + inter)) ** 2 for p in pts)
    r2 = (1.0 - ss_res / ss_tot) if ss_tot > 0.0 else 0.0
    return slope / 2.0, r2, len(pts), (pts[0][0] / 86400.0, pts[-1][0] / 86400.0)


def assert_matrix(series, res, phys, tier):
    """The extra gates the VERTICAL-COORDINATE MATRIX rows carry.

    Opt-in per row (`phys["matrix_gates"]`), deliberately: every one of these is a
    correct thing to ask of ANY rest case, but switching them on for the
    whole corpus at once would turn existing rows red on assertions their
    scoped `known_failure`s do not cover, which is exactly the failure mode
    `_status`'s assertion scoping exists to prevent.  Extending them to the
    rest of the corpus is a separate change with its own measurement.

    Four assertions, each answering something `energy:rest` cannot:

      * `energy:rest-growth-rate` -- the RATE, not the level.  A spurious
        pressure gradient that EQUILIBRATES is a tolerable discretisation
        error; one that is still exponential at the end of the run is an
        instability that happens to be small yet.  Measured on the cavity
        case: a 10 900x smaller seed bought 143 days and the run still
        breached the same level bar, because the RATE was unchanged.
        TIER 1 ONLY -- a 3.3-day twin cannot fit an exponential whose bar is
        a 20-day e-folding.
      * `tracer:no-new-extrema` -- at rest, with no surface flux and no
        mixing closure, T and S may not leave their initial range.  Spurious
        diapycnal mixing from the coordinate's own regrid shows up here
        before it shows up in the energy.
      * `thickness:positive` -- the minimum `h_layer` over the run stays
        strictly positive.  A coordinate that drives a layer to zero divides
        by it downstream.
      * `counters:no-truncation` -- a rest case that CFL-truncates is not
        resting, whatever its energy says.
    """
    out = []
    sig_bar = phys.get("rest_sigma_max")
    if sig_bar:
        if tier != 1:
            out.append(Verdict(
                "energy:rest-growth-rate", True,
                "fitted amplitude rate <= {:.3g} 1/s".format(sig_bar),
                "not measurable at tier-2 length (the bar is a 20-day "
                "e-folding) -- TIER 1 ONLY", skipped=True))
        else:
            days = [s["day"] for s in series["stats"]]
            en = [s["En"] for s in series["stats"]]
            sig, r2, npts, span = fit_log_rate(days, en)
            ok = sig is None or sig <= sig_bar
            out.append(Verdict(
                "energy:rest-growth-rate", ok,
                "a resting case has NO energy source, so its spurious energy "
                "must equilibrate rather than grow: fitted amplitude rate "
                "<= {:.3g} 1/s (En e-folding >= {:.1f} days)".format(
                    sig_bar, 1.0 / (2.0 * sig_bar) / 86400.0),
                "sigma = {} 1/s ({}), R2 = {}, {} samples over days {}".format(
                    _fmt(sig),
                    "En e-folding {:.2f} days".format(
                        1.0 / (2.0 * sig) / 86400.0)
                    if sig and sig > 0 else "not growing",
                    _fmt(r2), npts,
                    "{:.1f}-{:.1f}".format(*span) if span else "?"),
                "The spurious motion is still growing EXPONENTIALLY when the "
                "run ends. A level bar cannot see this: weakening the source "
                "only postpones it. On the measured cavity case a 10 900x "
                "smaller seed bought 143 days and breached the same bar.",
                "validation_examples/ocean/ice_shelf_cavity/README.md"))

    if not phys.get("matrix_gates"):
        return out

    # --- tracer bounds -----------------------------------------------------
    slack = phys.get("tracer_slack", 1.0e-6)
    worst, worst_name = 0.0, None
    seen = False
    for name in ("temperature", "salinity"):
        rows = [r for r in (series["diag"].get(name) or [])
                if _finite(r["min"]) and _finite(r["max"])]
        if len(rows) < 2:
            continue
        seen = True
        lo0, hi0 = rows[0]["min"], rows[0]["max"]
        for r in rows[1:]:
            d = max(lo0 - r["min"], r["max"] - hi0, 0.0)
            if d > worst:
                worst, worst_name = d, name
    if seen:
        out.append(Verdict(
            "tracer:no-new-extrema", worst <= slack,
            "at rest, with no surface flux and no mixing closure, T and S "
            "may not leave their INITIAL range by more than {:.1g}".format(slack),
            "worst overshoot = {} ({})".format(
                _fmt(worst), worst_name or "none"),
            "A new tracer extremum in a motionless, unforced, unmixed ocean "
            "was created by the coordinate's own regrid -- an unlimited "
            "reconstruction, or a vanished layer the remap drained and "
            "refilled. It is spurious diapycnal mixing, and it shows here "
            "before it shows in the energy."))

    # --- thickness positivity ----------------------------------------------
    rows = [r for r in (series["diag"].get("h_layer") or []) if _finite(r["min"])]
    if rows:
        h_min = min(r["min"] for r in rows)
        floor = phys.get("h_min_floor", 0.0)
        out.append(Verdict(
            "thickness:positive", h_min > floor,
            "the minimum layer thickness over the run stays strictly above "
            "{:.3g} m".format(floor),
            "min h_layer = {} m over {} samples".format(_fmt(h_min), len(rows)),
            "A layer reached zero or negative thickness. Every operator that "
            "divides by h_layer is armoured at H_DIV_EPS = 1e-20, so this "
            "does not crash -- it produces a plausible wrong answer."))

    # --- the solver's own counters -----------------------------------------
    nt = series.get("ntrunc", 0) or 0
    out.append(Verdict(
        "counters:no-truncation", nt == 0,
        "a case that starts at rest and has no forcing must never hit the "
        "barotropic CFL truncation",
        "{} CFL truncation(s)".format(nt),
        "The solver truncated a face velocity to hold the CFL. In a resting, "
        "unforced ocean that velocity is entirely spurious, and the "
        "truncation is a clamp -- it bounds the symptom and hides the size "
        "of the error from every other assertion here."))
    return out


def _diag_amplitude_rate(series, field, fit_from_day=0.0):
    """Amplitude growth rate sigma (1/s) of a [diag] FIELD's extrema.

    Why not `En`: `En` is the TOTAL kinetic energy, and in a thermal-wind
    case it is dominated by the basic-state jet (eady: 4.3e-3 m2/s2 from a
    +-0.16 m/s shear) until the eddies reach finite amplitude -- which for a
    correctly growing Eady mode is ~day 50 of a 60-day run. An En fit sees
    the wall-drain of the jet (-6% over 40 days) and calls the case "decaying"
    while the unstable mode is going up three decades underneath. The
    cross-channel velocity `v` is zero in the basic state, so its extrema are a
    pure perturbation measure and grow at the mode's own rate once the mode
    has out-grown the random-noise seed. `fit_from_day` excludes the phase
    where max|v| is still the decaying noise floor (eady: the 125 km mode
    holds ~50% of max|v| from day ~25 on, measured on the NetCDF modal
    amplitudes; before that the fit reads the noise, not the mode).

    Amplitude a(t) = max(|min|, |max|) at each [diag] emission; least-squares
    fit of log a(t). Returns (sigma, (day_first, day_last), n_points,
    amplification) or (None, None, 0, None) when there is not enough signal.
    """
    rows = series.get("diag", {}).get(field) or []
    pts = []
    for r in rows:
        t, lo, hi = r.get("t"), r.get("min"), r.get("max")
        if not (_finite(t) and _finite(lo) and _finite(hi)):
            continue
        a = max(abs(lo), abs(hi))
        if a > 0.0 and t / 86400.0 >= fit_from_day:
            pts.append((t, math.log(a)))
    if len(pts) < 4:
        return None, None, len(pts), None
    n = float(len(pts))
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    if den == 0.0:
        return None, None, len(pts), None
    sigma = (n * sxy - sx * sy) / den
    span = (pts[0][0] / 86400.0, pts[-1][0] / 86400.0)
    return sigma, span, len(pts), math.exp(pts[-1][1] - pts[0][1])


def assert_energy(series, phys, tier):
    """The stability gate. Regime-specific, because "must not grow" is only a
    correct statement for a case with nothing driving it."""
    en = [s["En"] for s in series["stats"]]
    days = [s["day"] for s in series["stats"]]
    if not en:
        return [Verdict("energy", False, "an En time series", "none parsed")]
    regime = phys.get("regime", "forced")
    out = []
    ref = _early_reference(en)
    final = en[-1]
    peak = max((e for e in en if _finite(e)), default=float("nan"))

    ref_note = (
        "An unforced run has no energy source. Energy that appears anyway was "
        "MANUFACTURED by the discretisation. Budgets do not see this: on "
        "2026-09-11 En grew 30-60x while Mass, Salt and Temp stayed exact to "
        "every printed digit and nothing crashed.")
    ref_link = "tmp_local_artifacts/eady_hunt/FINDINGS_windowed_advect.md"

    if regime == "rest":
        # Starts at rest and must stay there. Any energy at all is spurious --
        # a seamount's spurious pressure gradient, a wall-normal velocity leak,
        # an ALE remap that is not idempotent on a resting column. TWO
        # assertions, because they say different things:
        #
        #   (a) a MAGNITUDE bar. A terrain-following coordinate over a tall
        #       seamount has a known, quantified spurious pressure gradient;
        #       the seamount test's conventional acceptance bar is ~1 cm/s
        #       (Beckmann & Haidvogel 1993 and the sigma-coordinate error
        #       literature that follows it), which is what `en_rest_max`
        #       encodes for a sloped case. A case with NO slope has nothing to
        #       generate and gets a 1 um/s bar instead.
        #   (b) a TREND bar, which is the sharper of the two: a spurious
        #       pressure gradient that EQUILIBRATES is a tolerable
        #       discretisation error; one that is still at its maximum when the
        #       run ends has not equilibrated and is an instability wearing a
        #       small number. `seamount_pgf_ppm` peaks at 5.8e-7 and decays to
        #       2.8e-7 (equilibrated); `seamount_conservative_floor` ends AT
        #       its peak, 4000x higher.
        cap = phys.get("en_rest_max", 1e-8)
        ok = _finite(peak) and peak <= cap
        out.append(Verdict(
            "energy:rest", ok,
            "a quiescent case must stay at rest: peak En <= {:.2g} m2/s2 "
            "(|u|_rms <= {:.2g} m/s)".format(cap, math.sqrt(2.0 * cap)),
            "peak En = {} (final {}), i.e. |u|_rms ~ {} m/s".format(
                _fmt(peak), _fmt(final),
                _fmt(math.sqrt(2.0 * peak) if _finite(peak) and peak > 0 else 0.0)),
            "Motion appeared in a run with no forcing and a motionless initial "
            "state. The usual sources are a spurious pressure gradient over "
            "topography, a wall-normal velocity leak, or an ALE remap that is "
            "not idempotent on a resting column.",
            "src/core/ocean/README.md"))
        # Only meaningful where there is spurious motion to settle. `seamount`
        # and friends sit at En ~ 1e-20 -- machine zero, |u| ~ 1e-10 m/s --
        # where "peak == final" just means the number never changed. Asserting
        # a trend on roundoff is noise, so the trend gate starts at 1e-12
        # (|u|_rms ~ 1.4e-6 m/s), still a thousand times below the magnitude
        # bar it complements.
        trend_floor = phys.get("en_rest_trend_floor", 1e-12)
        if _finite(peak) and _finite(final) and peak > trend_floor:
            # "Still at the peak when the clock ran out" == not equilibrated.
            settled = final <= 0.95 * peak
            out.append(Verdict(
                "energy:rest-settles", settled,
                "spurious motion must EQUILIBRATE: final En below 95% of the "
                "run's peak",
                "peak En = {} -> final {} ({:.0f}% of peak)".format(
                    _fmt(peak), _fmt(final), 100.0 * final / peak),
                "The spurious velocity is still at its maximum at the end of "
                "the run, so nothing bounded it -- it is growing, not "
                "saturating. A tolerable discretisation error settles; an "
                "instability does not. Re-run longer to see which this is.",
                "src/core/ocean/README.md"))

    elif regime in ("adiabatic", "unforced"):
        # An IC with energy, nothing adding more: En may decay, must not grow.
        gmax = phys.get("en_growth_max", 1.5)
        ratio = (final / ref) if (ref and _finite(final) and ref > 0) else None
        pratio = (peak / ref) if (ref and _finite(peak) and ref > 0) else None
        ok = pratio is not None and pratio <= gmax
        out.append(Verdict(
            "energy:no-growth", ok,
            "an unforced/adiabatic run's kinetic energy MUST NOT GROW: "
            "peak En <= {:.2f} x its early-run reference".format(gmax),
            "En {} (early ref) -> peak {} -> final {}  =  {}x peak, {}x final"
            .format(_fmt(ref), _fmt(peak), _fmt(final),
                    _fmt(pratio), _fmt(ratio)),
            ref_note, ref_link))

    elif regime == "baroclinic":
        # Energy SHOULD grow -- that is the case's purpose. Two separate
        # assertions: the growth must be bounded (not a numerical blow-up),
        # and, at tier 1 where the integration is long enough to measure it,
        # the rate must match the analytical expectation the case is FOR.
        cap = phys.get("en_max", 1.0)
        ok = _finite(peak) and peak <= cap
        out.append(Verdict(
            "energy:bounded", ok,
            "a baroclinic-instability case may grow, but must stay physical: "
            "peak En <= {:.3g} m2/s2".format(cap),
            "peak En = {} (|u|_rms ~ {} m/s), final {}".format(
                _fmt(peak), _fmt(math.sqrt(2.0 * peak)
                                 if _finite(peak) and peak > 0 else 0.0),
                _fmt(final)),
            "Unbounded growth in a case that is SUPPOSED to grow is the hard "
            "case: the instability is real but the discretisation has stopped "
            "saturating it. Check MaxCFL and the grid-scale variance."))
        # The sharp gate, at BOTH tiers: an unforced case may grow at its
        # PHYSICAL rate and no faster. A numerical instability in an adiabatic
        # run is not subtle in rate terms -- it is one to two orders of
        # magnitude above the mode it corrupts -- but it IS subtle in
        # amplitude, which is why an absolute En ceiling misses it entirely
        # (eady's instability topped out at 1.1e-2, well under any sane
        # ceiling) and why a whole-run growth fit misses it too (see
        # _max_window_rate).
        # The analogue of energy:responds for an instability case: it must
        # actually DESTABILISE. A baroclinic case that decays is not a passing
        # test, it is a case that has never demonstrated the physics it exists
        # for -- which is exactly what `eady` turned out to be doing until
        # 2026-09-13 (a 4x-too-weak front under nu_h = 100 decayed at
        # sigma = -2.8e-7). Gated at tier 1 only: a growth rate cannot be
        # measured over 40 simulated hours.
        if phys.get("expect_growth", True) and tier == 1:
            wr_any = _max_window_rate(days, en)
            ok = wr_any is not None and wr_any > 0.0
            out.append(Verdict(
                "energy:destabilises", ok,
                "a baroclinic-instability case must GROW somewhere in the run: "
                "fastest windowed amplitude rate > 0",
                "fastest windowed sigma = {} 1/s ({})".format(
                    _fmt(wr_any),
                    "grows" if (wr_any or 0) > 0 else "DECAYS EVERYWHERE"),
                "The case exists to demonstrate an instability and the energy "
                "only ever went down. Either the front is too weak for its "
                "mode to clear the noise seed in the run (eady's defect until "
                "2026-09-13: dT_dy = -5e-6 gave an 18.6-day e-folding under a "
                "wall-draining nu_h = 100), the initial perturbation is too "
                "small, or the run is far too short for the e-folding time. "
                "Lengthening the run is the first thing to rule out -- and "
                "check the [diag] extrema of a perturbation field, because "
                "En can be all basic-state jet.",
                "validation_examples/ocean/eady/eady.nml header"))
        fast = phys.get("sigma_fast_max")
        if fast is None and phys.get("sigma_expected"):
            # 5x the physical band's upper edge. Calibrated against the
            # 2026-09-11 measurement, not against a passing run: eady's
            # thermo-cadence instability ran at sigma ~ 1.0e-5 (27-hour
            # e-folding) against a physical band topping out at 1.2e-6, so a
            # 5x bar (6e-6) catches it with 1.7x margin while leaving the
            # physical mode 10x of headroom underneath.
            fast = 5.0 * phys["sigma_expected"]["max"]
        if fast:
            wr = _max_window_rate(days, en)
            # A windowed exponential fit needs something to fit. When the run
            # never actually grew -- peak En within GROWTH_FLOOR of its early
            # reference -- what the sliding window measures is the case's
            # inertia-gravity OSCILLATION, and calling that a growth rate is
            # an extrapolation, not a measurement.
            #
            # Measured, and why this is not a loosened bar. `bc_inst_tuned_512`
            # at tier 2 runs 200 steps = 11 SIMULATED HOURS. Its En sits at
            # 1.37e-2 and wobbles +-3% for the whole run -- a DYNAMIC RANGE of
            # 1.064, no trend, under BOTH outer split schemes. The windowed fit
            # turns that 3% wobble into "sigma = 3.19e-6, e-folding 3.6 days"
            # and fails the generic 3.0e-6 bar, while the ssp_rk2 twin's
            # slightly smaller wobble reads 2.80e-6 and passes. The verdict was
            # reporting the phase of an oscillation.
            #
            # The discriminator is the run's dynamic range max(En)/min(En), NOT
            # peak-over-early-reference: `_early_reference` is a MAX over the
            # first 10% of the run, and the 2026-09-11 burst happens INSIDE
            # that window (En 2.68e-4 -> 8.59e-3 by day 2 of 25), so a
            # peak/reference ratio reads 1.29 for the very defect this gate
            # exists to catch. Its dynamic range is 41x. The floor therefore
            # sits ~27x below the defect and 1.4x above the oscillation, and
            # `self_test` replays both series to keep it there.
            GROWTH_FLOOR = 1.5
            pos = [e for e in en if _finite(e) and e > 0.0]
            rng = (max(pos) / min(pos)) if pos else None
            grew = rng is None or rng >= GROWTH_FLOOR
            ok = wr is None or wr <= fast or not grew
            if not grew:
                out.append(Verdict(
                    "energy:no-fast-growth", True,
                    "a windowed growth rate, once the run has actually grown",
                    "not asserted: En spans only {:.4g}x over the whole run "
                    "({} to {}), below the {}x floor, so the fastest window "
                    "({} 1/s) is measuring an oscillation, not growth"
                    .format(rng, _fmt(min(pos)), _fmt(max(pos)),
                            GROWTH_FLOOR, _fmt(wr)),
                    skipped=True))
            else:
                out.append(Verdict(
                    "energy:no-fast-growth", ok,
                    "an unforced case may grow at its PHYSICAL rate and no "
                    "faster: fastest windowed amplitude rate <= {:.3g} 1/s "
                    "(e-folding >= {:.2f} days)".format(
                        fast, 1.0 / fast / 86400.0),
                    "fastest windowed sigma = {} 1/s ({})".format(
                        _fmt(wr),
                        "e-folding {:.2f} days".format(1.0 / wr / 86400.0)
                        if wr and wr > 0 else "not growing"),
                    "Energy is growing far faster than the physical mode this "
                    "case exists to demonstrate. This is the signature the "
                    "2026-09-11 thermo-cadence instability left: En up 32x in "
                    "two days, budgets exact to every digit, no NaN, exit "
                    "code 0. An absolute energy ceiling does not see it and "
                    "neither does a whole-run growth fit.",
                    "tmp_local_artifacts/eady_hunt/FINDINGS_windowed_advect.md"))
        sig = phys.get("sigma_expected")
        if sig and tier == 1:
            lo, hi = sig["min"], sig["max"]
            field = sig.get("field")
            if field:
                # Perturbation-amplitude fit on a [diag] field's extrema (see
                # _diag_amplitude_rate for why En cannot be used here).
                amp, span, npts, growth = _diag_amplitude_rate(
                    series, field, sig.get("fit_from_day", 0.0))
                how = "max|{}| over days [{}, {}], {} samples, x{}".format(
                    field,
                    "{:.1f}".format(span[0]) if span else "?",
                    "{:.1f}".format(span[1]) if span else "?",
                    npts, _fmt(growth) if growth is not None else "?")
            else:
                rate = _efold_rate(days, en)
                # En ~ exp(2*sigma*t) for an amplitude growth rate sigma.
                amp = (rate / 2.0) if rate is not None else None
                how = "whole-run fit of log(En)/2"
            ok = amp is not None and lo <= amp <= hi
            out.append(Verdict(
                "energy:growth-rate", ok,
                "amplitude growth rate sigma in [{:.3g}, {:.3g}] 1/s -- the "
                "analytical expectation this case exists to demonstrate ({})"
                .format(lo, hi, sig.get("basis", "theory")),
                "measured sigma = {} 1/s ({}; {})".format(
                    _fmt(amp),
                    "DECAYS" if (amp is not None and amp < 0) else "grows",
                    how),
                sig.get("meaning", ""), sig.get("ref", "")))
        elif sig and tier == 2:
            out.append(Verdict(
                "energy:growth-rate", True,
                "amplitude growth rate in [{:.3g}, {:.3g}] 1/s".format(
                    sig["min"], sig["max"]),
                "not measurable at tier-2 length (needs several e-folding "
                "times; this twin runs {:.2f} days) -- TIER 1 ONLY".format(
                    days[-1] if days and days[-1] else 0.0),
                skipped=True))

    else:  # "forced"
        # Driven by wind / heat flux: growth is expected during spin-up, so the
        # honest gate is (a) an absolute physical ceiling and (b) at tier 1,
        # where the run is long enough to have saturated, no runaway.
        cap = phys.get("en_max", 1.0)
        ok = _finite(peak) and peak <= cap
        out.append(Verdict(
            "energy:bounded", ok,
            "a forced case spins up, but must stay physical: peak En <= "
            "{:.3g} m2/s2 (|u|_rms <= {:.2f} m/s)".format(
                cap, math.sqrt(2.0 * cap)),
            "peak En = {} (|u|_rms ~ {} m/s), final {}".format(
                _fmt(peak), _fmt(math.sqrt(2.0 * peak)
                                 if _finite(peak) and peak > 0 else 0.0),
                _fmt(final)),
            "The forced solution left the physical envelope. At this point a "
            "NaN is only a matter of time.", ref_link))
        # A forced case must actually RESPOND. This looks trivial and is not:
        # `coriolis_coast` ships wind_config="constant", taux_magnitude=0.05
        # and produces En = 0.000E+00 at EVERY step of a 10-day integration --
        # the wind never reaches the ocean. No NaN gate, no golden and no
        # conservation check can see a forcing path that is silently a no-op,
        # because doing nothing is perfectly conservative and perfectly finite.
        if phys.get("expect_response", True):
            floor = phys.get("en_response_min", 1e-12)
            ok = _finite(peak) and peak > floor
            out.append(Verdict(
                "energy:responds", ok,
                "a FORCED case must actually move: peak En > {:.1g} m2/s2"
                .format(floor),
                "peak En = {} over {} samples ({:.2f} simulated days)".format(
                    _fmt(peak), len(en),
                    days[-1] if days and days[-1] else 0.0),
                "The case declares a forcing and the ocean did not respond to "
                "it. A forcing path that is silently a no-op passes every NaN, "
                "budget and golden check ever written -- doing nothing is "
                "exactly conservative. This assertion is the only thing that "
                "sees it."))
        lg = phys.get("late_growth_max")
        if lg and tier == 1 and len(en) >= 8 and _finite(peak) and peak > 0.0:
            half = len(en) // 2
            a = _early_reference(en[half:], frac=0.25)
            b = max((e for e in en[half:] if _finite(e)), default=None)
            ratio = (b / a) if (a and b and a > 0) else None
            ok = ratio is not None and ratio <= lg
            out.append(Verdict(
                "energy:saturating", ok,
                "over the SECOND half of the run a forced case must be "
                "saturating: En grows by <= {:.2f}x".format(lg),
                "second-half En {} -> {} = {}x".format(
                    _fmt(a), _fmt(b), _fmt(ratio)),
                "Still climbing steeply at the end of a spin-up-length "
                "integration. Either the case is under-run (raise its tier-1 "
                "duration and re-check) or energy is being manufactured.",
                ref_link))
    return out


def assert_cfl(series, phys):
    """Bounded CFL, and no monotone climb (the signature of a trip)."""
    cfl = [s["MaxCFL"] for s in series["stats"] if s["MaxCFL"] is not None]
    if not cfl:
        return [Verdict("cfl", True, "a MaxCFL series",
                        "not reported by this case", skipped=True)]
    cap = phys.get("cfl_max", 0.9)
    peak = max((c for c in cfl if _finite(c)), default=float("nan"))
    out = [Verdict(
        "cfl:bounded", _finite(peak) and peak <= cap,
        "MaxCFL <= {:.2f} at every step".format(cap),
        "peak MaxCFL = {} at step {}".format(
            _fmt(peak),
            next((s["step"] for s in series["stats"]
                  if s["MaxCFL"] == peak), "?")),
        "The advective CFL exceeded its bound. Past this the scheme is not "
        "merely inaccurate, it is unstable.")]
    # Monotone climb over the final quarter: physics oscillates, a trip does
    # not. But "climbing" alone is not enough to call a runaway -- a quiescent
    # seamount's spurious mode is still spinning up at step 300 and climbs
    # from MaxCFL 6e-4 to 2e-3, which is nowhere near a trip and never will be.
    # So the assertion is PROJECTIVE: take the growth factor the final quarter
    # actually achieved and continue it for four more quarters (one further
    # run length). If MaxCFL would pass its ceiling in that time, the run was
    # on its way to blowing up and merely ran out of clock; if it would not,
    # a monotone climb is spin-up and is reported, not failed.
    # NOT applied to a `baroclinic` case. There the declared physics IS
    # exponential growth, so "MaxCFL is growing exponentially" carries no
    # information: `baroclinic_15layer` ends its shipped 50-day run with the
    # eddy field still spinning up (MaxCFL 0.061 -> 0.207 over the last
    # quarter), which projects to 27 and means nothing. The energy assertions
    # -- `no-fast-growth` in particular -- are the right instrument for that
    # regime, and `cfl:bounded` above still applies.
    tail = [c for c in cfl[-max(6, len(cfl) // 4):] if _finite(c)]
    if phys.get("regime") != "baroclinic" and len(tail) >= 6 and tail[0] > 0.0:
        climbing = all(b >= a for a, b in zip(tail, tail[1:]))
        factor = tail[-1] / tail[0]
        projected = tail[-1] * (factor ** 4) if factor > 1.0 else tail[-1]
        # ...and only once MaxCFL is actually within reach of its ceiling.
        # Extrapolating an exponential from 5% of the cap is not evidence of
        # anything: `baroclinic_2layer` climbs 0.018 -> 0.049 over the last
        # quarter of its shipped 40-day run because its eddy field is still
        # SPINNING UP, and projecting that growth four more quarters forward
        # "predicts" 2.76 -- a number the case has no intention of reaching.
        # A genuine trip passes 10% of the ceiling long before it NaNs.
        near = tail[-1] > 0.1 * cap
        runaway = climbing and factor > 1.5 and projected > cap and near
        out.append(Verdict(
            "cfl:no-runaway", not runaway,
            "a monotone MaxCFL climb that is already past 10% of the "
            "ceiling must not project to a trip: continuing the final "
            "quarter's growth for one more run length must keep MaxCFL "
            "<= {:.2f}".format(cap),
            "final-quarter MaxCFL {} -> {} ({:.2f}x, {}); projects to {} after "
            "one more run length; now at {:.1f}% of the ceiling".format(
                _fmt(tail[0]), _fmt(tail[-1]), factor,
                "monotone" if climbing else "not monotone", _fmt(projected),
                100.0 * tail[-1] / cap if cap else 0.0),
            "A monotonic climb to a trip is numerical instability, not "
            "physics: real flow wobbles, a runaway does not. The run may well "
            "have exited 0 -- it just had not blown up YET."))
    return out


def assert_claims(series, phys, tier=1):
    """The case's OWN stated expectation, where its header makes a testable one.

    Never weakened to pass and never satisfied by editing the namelist: when a
    header and the code disagree, that IS the finding. A claim may name the
    tiers it is measurable at (`"tiers": [1]`): a "grows 1500x over 60 days"
    figure cannot be checked on a 40-hour twin, and reporting that as a
    failure would teach the reader to ignore claim failures. It is reported
    as SKIP ... TIER 1 ONLY instead, like energy:growth-rate.
    """
    out = []
    for claim in phys.get("claims", []) or []:
        kind = claim["kind"]
        if tier not in claim.get("tiers", (1, 2)):
            out.append(Verdict(
                "claim:" + claim["name"], True,
                claim["text"],
                "not measurable at tier-{} length -- TIER {} ONLY".format(
                    tier, "/".join(str(t) for t in claim["tiers"])),
                skipped=True))
            continue
        if kind == "mass_rel":
            mass = [s["Mass"] for s in series["stats"] if s["Mass"] is not None]
            if len(mass) < 2 or not mass[0]:
                continue
            drift = abs(mass[-1] - mass[0]) / abs(mass[0])
            out.append(Verdict(
                "claim:" + claim["name"], drift <= claim["tol"],
                "{} (the namelist header's own claim)".format(claim["text"]),
                "relative mass drift = {}".format(_fmt(drift)),
                claim.get("meaning", ""), claim.get("ref", "")))
        elif kind == "field_ratio":
            # Amplification of a [diag] field's extrema between `from_day`
            # and the end of the run -- the header's own "grows N x" figure,
            # for a perturbation field that starts at zero in the basic state
            # (eady's cross-channel `v`), where En says nothing.
            field = claim["field"]
            rows = [r for r in (series["diag"].get(field) or [])
                    if _finite(r["t"]) and _finite(r["min"]) and _finite(r["max"])]
            start = [r for r in rows if r["t"] / 86400.0 >= claim.get("from_day", 0.0)]
            ratio = None
            if len(start) >= 2:
                a0 = max(abs(start[0]["min"]), abs(start[0]["max"]))
                a1 = max(abs(start[-1]["min"]), abs(start[-1]["max"]))
                ratio = (a1 / a0) if a0 > 0.0 else None
            lo, hi = claim["min"], claim["max"]
            out.append(Verdict(
                "claim:" + claim["name"],
                ratio is not None and lo <= ratio <= hi,
                "{} => max|{}| amplification in [{:.4g}, {:.4g}] from day {:g} "
                "to the end of the run".format(
                    claim["text"], field, lo, hi, claim.get("from_day", 0.0)),
                "measured amplification = {} ({} -> {})".format(
                    _fmt(ratio),
                    _fmt(max(abs(start[0]["min"]), abs(start[0]["max"]))) if len(start) >= 2 else "?",
                    _fmt(max(abs(start[-1]["min"]), abs(start[-1]["max"]))) if len(start) >= 2 else "?"),
                claim.get("meaning", ""), claim.get("ref", "")))
        elif kind == "en_ratio":
            en = [s["En"] for s in series["stats"]]
            ref = _early_reference(en)
            peak = max((e for e in en if _finite(e)), default=None)
            ratio = (peak / ref) if (ref and peak and ref > 0) else None
            lo, hi = claim["min"], claim["max"]
            out.append(Verdict(
                "claim:" + claim["name"],
                ratio is not None and lo <= ratio <= hi,
                "{} => En amplification in [{:.4g}, {:.4g}]".format(
                    claim["text"], lo, hi),
                "measured amplification = {}".format(_fmt(ratio)),
                claim.get("meaning", ""), claim.get("ref", "")))
    return out


# ---------------------------------------------------------------------------
# Namelist patching (time window + arbitrary group overrides)
# ---------------------------------------------------------------------------
def _fmt_value(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, int):
        return str(v)
    return str(v)


def patch_namelist(src_text, spec, scratch_dir):
    """Return `src_text` patched for one stability run.

    * `&<group> key = value` for every entry of `spec["overrides"]`, injecting
      the key (or the whole group) when the source omits it.
    * `&time_nml`: `time_unit="second"`, `dt_fixed` from the override if given,
      `t_end = n_steps * dt`.
    * `&output_nml output_dir` -> scratch.
    * `&logging_nml status_interval` / `&ocean_diag_nml dt_out` -> a cadence
      that yields ~`spec["samples"]` [stats]/[diag] emissions across the run.
      A TIME SERIES is the whole point here -- unlike compare.py, which only
      needs the final state -- but emitting every step on a multi-thousand-step
      tier-1 run is both slow and unreadable, so the cadence is derived from a
      target sample count instead.
    """
    over = {g.lower(): dict(kv) for g, kv in (spec.get("overrides") or {}).items()}

    lines = src_text.splitlines()

    def scan(group, key):
        g = None
        for line in lines:
            s = line.split("!", 1)[0].strip()
            if s.startswith("&"):
                g = s[1:].split()[0].lower()
                continue
            if s.startswith("/"):
                g = None
                continue
            if g == group and "=" in s and s.split("=", 1)[0].strip().lower() == key:
                return s.split("=", 1)[1].strip().rstrip(",").rstrip("/").strip()
        return None

    dt = over.get("time_nml", {}).get("dt_fixed")
    if dt is None:
        raw = scan("time_nml", "dt_fixed")
        if raw is None:
            raise ValueError("no dt_fixed in &time_nml")
        dt = float(raw.replace("d", "e").replace("D", "e"))
    dt = float(dt)
    n_steps = int(spec["n_steps"])
    t_end = n_steps * dt
    samples = max(2, int(spec.get("samples", 60)))
    # Emit on a STRIDE of whole outer steps that DIVIDES n_steps, so the last
    # sample lands exactly on the final step. (A stride that does not divide it
    # leaves the series ending several steps early, which the completion
    # assertion would then have to guess about.)
    stride = max(1, int(round(n_steps / float(samples))))
    while stride > 1 and n_steps % stride:
        stride -= 1
    cadence = stride * dt

    time_over = over.setdefault("time_nml", {})
    time_over["t_end"] = t_end
    time_over["time_unit"] = '"second"'
    time_over["dt_fixed"] = dt
    out_nml = over.setdefault("output_nml", {})
    out_nml["output_dir"] = '"{}"'.format(scratch_dir)
    # This suite asserts on the model's CONSOLE time series and never reads a
    # NetCDF file, so writing them is pure cost -- and at tier 1 (tens of
    # samples across grids up to 600x600x50) it is enough cost to fill a disk,
    # which is exactly what happened on the first full sweep. Force it off.
    out_nml["output_to_file"] = False
    over.setdefault("logging_nml", {})["status_interval"] = cadence
    # The DIAG manager gets a MUCH coarser cadence than [stats]. Two reasons:
    #   * the assertions that need a dense time series (energy growth, CFL
    #     runaway, budget residuals) all read the [stats] line, which is free;
    #     the [diag] line only supplies field EXTREMA, for which a handful of
    #     samples is plenty.
    #   * `&ocean_diag_nml enabled` gates the console line and the NetCDF write
    #     together -- there is no console-only mode -- and at tier 1 a frame of
    #     a 600x600x50 case is ~0.9 GB. Eighty frames of that filled the disk
    #     twice on the first sweeps. Four frames is bounded and sufficient.
    diag_samples = max(1, int(spec.get("diag_samples", 4)))
    diag_stride = max(stride, (n_steps // diag_samples) or 1)
    if any(l.split("!", 1)[0].strip().lower().startswith("&ocean_diag_nml")
           for l in lines):
        over.setdefault("ocean_diag_nml", {})["dt_out"] = diag_stride * dt

    out, group, seen = [], None, {}
    for line in lines:
        s = line.split("!", 1)[0].strip()
        if s.startswith("&"):
            group = s[1:].split()[0].lower()
            seen[group] = set()
            out.append(line)
            continue
        if s.startswith("/"):
            for k, v in sorted(over.get(group, {}).items()):
                if k not in seen.get(group, ()):
                    out.append("   {} = {}".format(k, _fmt_value(v)))
            out.append(line)
            group = None
            continue
        if group and "=" in s:
            key = s.split("=", 1)[0].strip().lower()
            if key in over.get(group, {}):
                seen[group].add(key)
                out.append("   {} = {}".format(
                    key, _fmt_value(over[group][key])))
                continue
        out.append(line)

    # Any group the source omits entirely gets appended.
    present = {l.split("!", 1)[0].strip()[1:].split()[0].lower()
               for l in lines if l.split("!", 1)[0].strip().startswith("&")}
    for g, kv in sorted(over.items()):
        if g in present:
            continue
        out.append("")
        out.append("&" + g)
        for k, v in sorted(kv.items()):
            out.append("   {} = {}".format(k, _fmt_value(v)))
        out.append("/")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Case execution
# ---------------------------------------------------------------------------
def run_case(case, tier, binary, scratch_root, gpu_id=None, keep=False):
    spec = case["tier{}".format(tier)]
    name = case["name"]
    scratch = os.path.join(scratch_root, "tier{}".format(tier), name)
    shutil.rmtree(scratch, ignore_errors=True)
    os.makedirs(scratch, exist_ok=True)

    res = {"case": name, "tier": tier, "nml": case["nml"],
           "n_steps": spec["n_steps"], "timeout_s": spec["timeout_s"],
           "wall_s": 0.0, "returncode": None, "timed_out": False}

    try:
        src = open(os.path.join(REPO_ROOT, case["nml"])).read()
        patched = patch_namelist(src, spec, scratch)
    except (OSError, ValueError) as exc:
        res["fatal"] = "namelist patch failed: {}".format(exc)
        return res, {"stats": [], "budget": {}, "diag": {}, "crash": None}

    nml_path = os.path.join(scratch, name + ".nml")
    open(nml_path, "w").write(patched)

    setup = case.get("setup")
    if setup:
        argv = [a.replace("{python}", sys.executable).replace("{repo}", REPO_ROOT)
                for a in setup]
        sp = subprocess.run(argv, cwd=scratch, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=300)
        if sp.returncode != 0:
            res["fatal"] = "setup exited {}".format(sp.returncode)
            return res, {"stats": [], "budget": {}, "diag": {}, "crash": None}

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"      # NetCDF/HDF5 here is not thread-safe
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    t0 = time.time()
    try:
        proc = subprocess.run([binary, nml_path], cwd=scratch, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=spec["timeout_s"], check=False)
        out = proc.stdout.decode("utf-8", "replace")
        res["returncode"] = proc.returncode
    except subprocess.TimeoutExpired as exc:
        out = (exc.stdout or b"").decode("utf-8", "replace")
        res["timed_out"] = True
    res["wall_s"] = round(time.time() - t0, 2)

    open(os.path.join(scratch, "run.log"), "w").write(out)
    series = parse_series(out)
    res["scratch"] = scratch
    res["nan_catch"] = series.get("nan_catch", 0)
    if not keep:
        # ALWAYS drop the NetCDF, pass or fail: the suite asserts on the
        # console series and never reads a file, and leaving them behind on a
        # failing tier-1 case is what filled the disk. run.log is kept.
        for f in os.listdir(scratch):
            if f.endswith(".nc"):
                try:
                    os.remove(os.path.join(scratch, f))
                except OSError:
                    pass
    return res, series


def evaluate(case, tier, res, series):
    """Run every assertion for one case and return (verdicts, artifacts)."""
    spec = case["tier{}".format(tier)]
    phys = dict(case.get("physics") or {})
    phys.update(spec.get("physics_override") or {})

    if res.get("fatal"):
        return [Verdict("completed", False, "a runnable namelist",
                        res["fatal"])], []

    verdicts = [assert_ran(series, spec["n_steps"], res)]
    fin, artifacts = assert_finite(series)
    verdicts.append(fin)
    if verdicts[0].ok:
        verdicts += assert_conservation(series, phys.get("budget_tol", {}))
        verdicts += assert_energy(series, phys, tier)
        verdicts += assert_cfl(series, phys)
        verdicts += assert_claims(series, phys, tier)
        verdicts += assert_matrix(series, res, phys, tier)
    return verdicts, artifacts


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def _force_scheme(case, scheme):
    """Return a copy of `case` with `split_scheme` forced on both tiers.

    A case that already pins a scheme -- in its namelist or through the
    manifest scheme axis -- is returned UNTOUCHED: forcing those would
    silently retarget a case whose whole point is the scheme it names.
    """
    import copy
    if case.get("split_scheme"):
        return case
    src = open(os.path.join(REPO_ROOT, case["nml"])).read()
    for line in src.splitlines():
        if line.split("!", 1)[0].strip().lower().startswith("split_scheme"):
            return case
    out = copy.deepcopy(case)
    for t in ("tier1", "tier2"):
        spec = out.get(t)
        if not spec or spec.get("skip"):
            continue
        spec.setdefault("overrides", {}).setdefault(
            "ocean_bt_nml", {})["split_scheme"] = '"{}"'.format(scheme)
    out["split_scheme"] = scheme
    return out


def select_cases(tier, only=None, tags=None):
    """Cases runnable at `tier`, optionally filtered by NAME or by TAG.

    The tag filter exists because the vertical-coordinate matrix is 115
    cases whose names are generated: naming them all on a `--cases` line is
    not something a human can type, and running the whole tier-1 corpus to
    re-measure one table is an hour of GPU for nothing.
    """
    out = []
    for c in manifest.STABILITY_CASES:
        spec = c.get("tier{}".format(tier))
        if not spec or spec.get("skip"):
            continue
        if only and c["name"] not in only:
            continue
        if tags and not (tags & set(c.get("tags") or ())):
            continue
        out.append(c)
    return out


def run_all(cases, tier, binary, scratch_root, gpus, keep, jobs):
    results = {}
    lock = threading.Lock()
    work = queue.Queue()
    for c in cases:
        work.put(c)

    def worker(gpu_id):
        while True:
            try:
                c = work.get_nowait()
            except queue.Empty:
                return
            try:
                res, series = run_case(c, tier, binary, scratch_root,
                                       gpu_id=gpu_id, keep=keep)
                verdicts, artifacts = evaluate(c, tier, res, series)
                with lock:
                    results[c["name"]] = (c, res, verdicts, artifacts)
                    _print_case(c, tier, res, verdicts, artifacts)
            finally:
                work.task_done()

    if gpus:
        slots = list(gpus)
    else:
        slots = [None] * max(1, jobs)
    threads = [threading.Thread(target=worker, args=(s,), daemon=True)
               for s in slots]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return results


def _status(case, verdicts, tier=None):
    """PASS / FAIL / XFAIL / XPASS for one case.

    A `known_failure` may name the tiers it applies to (`"tiers": [1]`) --
    several known defects are only DETECTABLE at tier 1, and marking them
    known-failing everywhere would report a permanent XPASS at tier 2, which
    reads as "the defect is fixed, drop the marker" and is wrong.

    It may also name the ASSERTIONS it covers (`"assertions": [...]`). Without
    that, one documented defect excuses a case from EVERY gate it has, which
    is how a real regression hides behind an old XFAIL. `resting_stratified_
    channel` is the case that forced this: its spurious energy is genuinely
    still creeping at day 120 under BOTH outer schemes (`energy:rest-settles`
    can never pass), while the MAGNITUDE bar `energy:rest` is the whole point
    of the case and separates the two schemes by four decades. Scoping the
    marker keeps the second gate live.
    """
    failed = [v for v in verdicts if not v.ok and not v.skipped]
    known = case.get("known_failure")
    if known and tier is not None and tier not in known.get("tiers", (1, 2)):
        known = None
    covered = known.get("assertions") if known else None
    if failed and known:
        if covered is not None:
            uncovered = [v for v in failed if v.name not in covered]
            if uncovered:
                return "FAIL", uncovered
        return "XFAIL", failed
    if failed:
        return "FAIL", failed
    if known:
        return "XPASS", []
    return "PASS", []


_PRINT_LOCK = threading.Lock()


def _print_case(case, tier, res, verdicts, artifacts):
    status, failed = _status(case, verdicts, tier)
    with _PRINT_LOCK:
        print("\n{} {:<30} {:>7.1f}s  {}".format(
            {"PASS": "[ PASS]", "FAIL": "[*FAIL]", "XFAIL": "[XFAIL]",
             "XPASS": "[XPASS]"}[status], case["name"], res["wall_s"],
            case["nml"]))
        for v in verdicts:
            if (v.ok or v.skipped) and not _VERBOSE:
                continue
            print(v.render(case["name"], tier))
        if res.get("nan_catch"):
            print("     NAN-CATCH: the solver's in-flight repair fired {} time(s). "
                  "The run continued, but non-finite values WERE produced and "
                  "zeroed -- treat as a defect even when every other assertion "
                  "passes.".format(res["nan_catch"]))
        if artifacts and _VERBOSE:
            print("     ARTIFACT : diag mean NaN on [{}] -- {}".format(
                ", ".join(artifacts), DIAG_MEAN_NAN_NOTE))
        if status == "XFAIL":
            kf = case["known_failure"]
            print("     KNOWN    : {}".format(kf["reason"]))
            print("     see      : {}".format(kf.get("ref", "-")))
        if status == "XPASS":
            print("     XPASS    : this case is marked known-failing "
                  "({}) but every assertion passed. If the defect is fixed, "
                  "drop `known_failure` from the manifest."
                  .format(case["known_failure"]["reason"]))


_VERBOSE = False


# ---------------------------------------------------------------------------
# Self-test: does the gate actually catch what it was built for?
# ---------------------------------------------------------------------------
def _series(day_en, cfl=None, budget=None):
    """Build a minimal parsed-series dict from (day, En) pairs."""
    stats = [{"day": d, "step": i, "En": e,
              "MaxCFL": (cfl[i] if cfl else 0.01),
              "Mass": 1.0e16, "Salt": 35.0, "Temp": 5.0}
             for i, (d, e) in enumerate(day_en)]
    return {"stats": stats, "diag": {},
            "budget": budget or {"Mass": [1e-15] * len(stats)},
            "crash": None, "total_steps": len(stats) - 1, "nan_catch": 0}


def self_test():
    """Assert the gate classifies the 2026-09-11 signatures correctly.

    Every series here is taken from a MEASURED run, not invented, so this is a
    regression test on the gate itself: if someone loosens a bound, one of
    these flips and says so. Who tests the test.
    """
    import math as _m
    fails = []

    def check(label, ok, detail=""):
        print("  {:<52} {}".format(label, "ok" if ok else "FAILED " + detail))
        if not ok:
            fails.append(label)

    def named(vs, name):
        return next((v for v in vs if v.name == name), None)

    # (1) The thermo-cadence instability on eady (FINDINGS_windowed_advect.md):
    #     En 2.676e-4 -> 8.585e-3 by day 2, then flat to 1.105e-2 at day 25.
    #     Budgets exact, no NaN, exit 0. The gate MUST fail it.
    broken = [(0.0, 2.676e-4), (1.0, 1.5e-3), (2.0, 8.585e-3)] + \
             [(2.0 + d, 8.6e-3 + d * 9.8e-5) for d in range(1, 24)]
    phys = {"regime": "baroclinic", "en_max": 0.5, "cfl_max": 0.9,
            "sigma_fast_max": 3.0e-6}
    v = assert_energy(_series(broken), phys, 1)
    check("2026-09-11 instability -> energy:no-fast-growth FAILS",
          named(v, "energy:no-fast-growth") is not None
          and not named(v, "energy:no-fast-growth").ok)
    check("...and energy:bounded does NOT see it (En peaks at 1.1e-2)",
          named(v, "energy:bounded").ok)
    check("...and the whole-run growth fit does NOT see it either",
          abs(_efold_rate([d for d, _ in broken],
                          [e for _, e in broken]) / 2.0) < 3.0e-6)

    # (2) The control from the same table: 2.676e-4 -> 2.454e-4 over 25 days.
    good = [(float(d), 2.676e-4 - d * 8.9e-7) for d in range(26)]
    v = assert_energy(_series(good), phys, 1)
    check("the healthy control PASSES energy:no-fast-growth",
          named(v, "energy:no-fast-growth").ok)

    # (2b) The growth FLOOR: a bounded oscillation is not a growth rate.
    #      Measured series -- bc_inst_tuned_512's tier-2 twin, 11 simulated
    #      hours, En wobbling +-3% around 1.4e-2 with no trend. The windowed
    #      fit reads 3.19e-6 (over the 3.0e-6 bar); the verdict must SKIP.
    import random as _rnd
    _rnd.seed(7)
    wobble = [(i * 0.0116, 1.40e-2 * (1.0 + 0.03 * math.sin(i * 1.1)))
              for i in range(41)]
    v = assert_energy(_series(wobble), phys, 2)
    nfg = named(v, "energy:no-fast-growth")
    check("a +-3% oscillation -> energy:no-fast-growth SKIPS (not FAIL)",
          nfg is not None and nfg.skipped and nfg.ok)
    check("...and the raw windowed fit on it WOULD have tripped the bar",
          (_max_window_rate([d for d, _ in wobble], [e for _, e in wobble])
           or 0.0) > 3.0e-6)
    check("...and the floor does NOT reach the 2026-09-11 burst "
          "(range 41x vs 1.06x)",
          max(e for _, e in broken) / min(e for _, e in broken) > 1.5
          and max(e for _, e in wobble) / min(e for _, e in wobble) < 1.5)
    _ = _rnd

    # (3) A forcing path that is silently a no-op (coriolis_coast, measured).
    v = assert_energy(_series([(float(d), 0.0) for d in range(11)]),
                      {"regime": "forced", "en_max": 0.5, "cfl_max": 0.9}, 1)
    check("a forced case stuck at En=0 -> energy:responds FAILS",
          not named(v, "energy:responds").ok)

    # (4) A quiescent case that never equilibrates (seamount_conservative_floor).
    stuck = [(float(d), 2.242e-3) for d in range(10)]
    v = assert_energy(_series(stuck),
                      {"regime": "rest", "en_rest_max": 0.5e-4,
                       "cfl_max": 0.9}, 1)
    check("a rest case at 6.7 cm/s -> energy:rest FAILS",
          not named(v, "energy:rest").ok)
    check("...and energy:rest-settles FAILS (ends at its peak)",
          not named(v, "energy:rest-settles").ok)

    # (5) A rest case at machine zero must NOT trip the trend gate (seamount).
    v = assert_energy(_series([(float(d), 1.778e-20) for d in range(10)]),
                      {"regime": "rest", "en_rest_max": 0.5e-4,
                       "cfl_max": 0.9}, 1)
    check("a rest case at En=1.8e-20 passes (roundoff is not a trend)",
          all(x.ok for x in v))

    # (6) A budget leak.
    v = assert_conservation(_series([(0.0, 1e-4)],
                                    budget={"Heat": [1e-15, 5.6e-5]}),
                            {"Heat": 1e-9})
    check("a Heat residual of 5.6e-5 -> conserve:Heat FAILS", not v[0].ok)

    # (7) A CFL runaway vs a spin-up climb, same monotone shape.
    trip = [0.05 * (1.6 ** i) for i in range(12)]
    v = assert_cfl(_series([(float(i), 1e-3) for i in range(12)], cfl=trip),
                   {"regime": "forced", "cfl_max": 0.9})
    check("a CFL climb projecting past the ceiling -> cfl:no-runaway FAILS",
          not named(v, "cfl:no-runaway").ok)
    spin = [6.1e-4 * (1.35 ** i) for i in range(12)]     # seamount_bench_full
    v = assert_cfl(_series([(float(i), 1e-3) for i in range(12)], cfl=spin),
                   {"regime": "forced", "cfl_max": 0.9})
    check("a quiescent spin-up climb at 0.2% of the ceiling PASSES",
          named(v, "cfl:no-runaway").ok)

    # (7b) The REST GROWTH-RATE fit, against the two MEASURED series it has
    #      to separate. Both are from
    #      validation_examples/ocean/ice_shelf_cavity/README.md and the
    #      cavity diagnosis: the sloping-lid sigma instability (En 1.030E-09
    #      at day 1 rising to 3.894E-06 at day 60, sigma_En = 0.333/day over
    #      the clean window) and the pred_corr residual creep on the resting
    #      channel (En 1.21e-09 -> 4.82e-09 over 120 days, an 83-day
    #      e-folding). The bar must fail the first and pass the second.
    bar = 0.05 / 2.0 / 86400.0
    inst = [(1.0, 1.030e-9), (5.0, 1.346e-9), (10.0, 1.626e-9),
            (20.0, 2.070e-9), (25.0, 5.275e-9), (30.0, 3.043e-8),
            (35.0, 1.901e-7), (40.0, 7.133e-7), (50.0, 2.770e-6),
            (60.0, 3.894e-6)]
    sig_i, r2_i, _n, _sp = fit_log_rate([d for d, _ in inst],
                                        [e for _, e in inst])
    check("the sloping-lid instability -> rest-growth-rate FAILS "
          "(sigma_En {:.3f}/day, R2 {:.3f})".format(
              (sig_i or 0) * 2 * 86400.0, r2_i or 0),
          sig_i is not None and sig_i > bar)
    creep = [(float(d), 1.21e-9 * math.exp(d / 83.0)) for d in range(0, 121, 5)]
    sig_c, r2_c, _n, _sp = fit_log_rate([d for d, _ in creep],
                                        [e for _, e in creep])
    check("the pred_corr 83-day creep PASSES it "
          "(sigma_En {:.4f}/day, R2 {:.3f})".format(
              (sig_c or 0) * 2 * 86400.0, r2_c or 0),
          sig_c is not None and sig_c <= bar)
    flat = [(float(d), 1.7e-9 * (1.0 + 0.02 * math.sin(d))) for d in range(60)]
    sig_f, _r2, _n, _sp = fit_log_rate([d for d, _ in flat],
                                       [e for _, e in flat])
    check("a PLATEAU reads as not growing", abs(sig_f or 0.0) <= bar)

    # (8) Every manifest twin satisfies the downscaling rules.
    bad = [c["name"] for c in manifest.STABILITY_CASES
           if c.get("tier2", {}).get("dimensionless")
           and not downscale.twin_ok(
               downscale.check_twin(c["tier2"]["dimensionless"]))]
    check("every tier-2 twin satisfies downscale.py's rules", not bad, str(bad))

    # (9) Manifest coverage.
    base = [c for c in manifest.STABILITY_CASES if not c.get("split_scheme")]
    twins = [c for c in manifest.STABILITY_CASES if c.get("split_scheme")]
    n2 = sum(1 for c in manifest.STABILITY_CASES
             if not c["tier2"].get("skip"))
    # 71 shipped namelists + the 115 generated vertical-coordinate matrix
    # cells. The tripwire is on the TOTAL so that a matrix cell silently
    # disappearing (a template token renamed, a family dropped from the
    # FAMILIES list) fails here rather than reporting a smaller green table.
    vcm = [c for c in base if c.get("matrix")]
    check("manifest covers 71 shipped namelists + the vcoord matrix "
          "({} base cases, {} of them matrix cells, {} at tier 2)"
          .format(len(base), len(vcm), n2),
          len(base) - len(vcm) == 71 and len(vcm) == 115)
    # Every matrix cell must name its coordinates in the table. Without this
    # a row whose `matrix` block lost a key would simply vanish from the
    # printed family x problem table with nothing said.
    check("every vcoord-matrix row carries problem/family/expect",
          all(set(("problem", "family", "expect", "stratification", "eos"))
              <= set(c["matrix"]) for c in vcm))
    # The REFUSAL rows are the ones that catch an envelope widening silently.
    refused = [c for c in vcm if c["matrix"]["expect"] == "refused"]
    check("the matrix asserts every by-design REFUSAL ({} rows)"
          .format(len(refused)),
          len(refused) >= 11
          and all(c["known_failure"]["assertions"] == ["completed"]
                  for c in refused))
    # The rate gate is what this matrix adds over the existing rest gates;
    # a refactor that dropped it would leave a table of level bars only.
    rated = [c for c in vcm if c["physics"].get("rest_sigma_max")]
    check("the matrix arms the fitted GROWTH-RATE gate on every runnable "
          "cell ({} rows)".format(len(rated)),
          len(rated) == len(vcm) - len(refused))
    # The scheme axis must EXIST -- a refactor that quietly stops building it
    # would leave the non-default scheme an untested branch of the dispatcher
    # while every report still said PASS.
    check("the outer split-scheme axis is built ({} {} twins)"
          .format(len(twins), _stab_manifest.SCHEME_AXIS_SCHEME),
          len(twins) >= 40)
    check("...every twin forces split_scheme on the tier it runs",
          all(t[tier]["overrides"]["ocean_bt_nml"]["split_scheme"]
              == '"{}"'.format(_stab_manifest.SCHEME_AXIS_SCHEME)
              for t in twins for tier in ("tier1", "tier2")
              if not t[tier].get("skip")))
    # A scoped known_failure must not silently cover an assertion that the
    # case does not have -- that is how a marker outlives its defect.
    scoped = [c for c in manifest.STABILITY_CASES
              if (c.get("known_failure") or {}).get("assertions")]
    check("scoped known_failures name assertions, not whole cases ({})"
          .format(len(scoped)),
          all(isinstance(c["known_failure"]["assertions"], list)
              and c["known_failure"]["assertions"] for c in scoped))
    # The resting-channel case exists to SEPARATE the two outer schemes on the
    # MAGNITUDE bar: ssp_rk2 reaches En 4.8e-08 at the tier-2 horizon
    # (2.992E-05 at tier 1) against a 0.5e-08 bar, pred_corr holds 1.4e-09.
    # Since the 2026-09-14 default flip the BASE case runs pred_corr, so it is
    # the base case that must keep gating `energy:rest`: if its known_failure
    # ever grew to cover the magnitude bar, a regression of the default back
    # to ssp_rk2 behaviour would report XFAIL and the case would gate nothing.
    # It is allowed to excuse ONLY the trend gate, which neither scheme can
    # pass (both are still creeping). The `__ssp_rk2` twin, by contrast, is
    # EXPECTED to miss both -- that is what the label on ssp_rk2 means.
    rc = next((c for c in manifest.STABILITY_CASES
               if c["name"] == "resting_stratified_channel"), None)
    check("the resting-channel base case (pred_corr) still GATES energy:rest",
          rc is not None
          and "energy:rest" not in rc["known_failure"]["assertions"]
          and "energy:rest-settles" in rc["known_failure"]["assertions"])
    tw = next((c for c in manifest.STABILITY_CASES
               if c["name"] == "resting_stratified_channel__"
               + _stab_manifest.SCHEME_AXIS_SCHEME), None)
    check("the resting-channel {} twin records the defect as XFAIL"
          .format(_stab_manifest.SCHEME_AXIS_SCHEME),
          tw is not None
          and "energy:rest" in tw["known_failure"]["assertions"]
          and "energy:rest-settles" in tw["known_failure"]["assertions"])
    check("...and the bar sits BETWEEN the two measured regimes",
          _stab_manifest.REST_100UM_S > 2.0e-9
          and _stab_manifest.REST_100UM_S < 2.0e-8)
    # The suite's idea of "what a base case runs" must be the model's.  A
    # default flip that forgot this line would mislabel every result file and
    # would silently make the axis run the default TWICE.
    import re as _re
    _cfg = open(os.path.join(REPO_ROOT, "src", "core", "rdb_config.F90")).read()
    _m = _re.search(r'character\(len=16\) :: split_scheme = "(\w+)"', _cfg)
    check("_DEFAULT_SPLIT_SCHEME tracks rdb_config.F90 ({})"
          .format(_m.group(1) if _m else "?"),
          _m is not None and _m.group(1) == _DEFAULT_SPLIT_SCHEME
          and _DEFAULT_SPLIT_SCHEME != _stab_manifest.SCHEME_AXIS_SCHEME)
    check("...and no twin is built from a namelist that pins the scheme",
          not any(_stab_manifest._nml_pins_scheme(t["nml"])
                  for t in twins))
    # The remap precondition guard rides every run whose coordinate remaps,
    # and no run whose coordinate does not (there it would judge stale
    # scratch).  Derived per tier from the namelist + overrides.
    _knob = _stab_manifest.REMAP_CHECK_KNOB
    _bad = []
    for c in manifest.STABILITY_CASES:
        for t in ("tier1", "tier2"):
            s = c.get(t)
            if not s or s.get("skip"):
                continue
            on = bool((s.get("overrides") or {}).get("vcoord_nml", {})
                      .get(_knob, False))
            if on != _stab_manifest._case_remaps(c, s):
                _bad.append("{}:{}".format(c["name"], t))
    check("{} is ON exactly where the coordinate remaps".format(_knob),
          not _bad)

    print("\nself-test: {}".format(
        "ALL PASS" if not fails else "{} FAILED".format(len(fails))))
    _ = _m
    return 1 if fails else 0


def main(argv=None):
    global _VERBOSE
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tier", type=int, choices=(1, 2), default=2,
                   help="1 = full-scale GPU (local/nightly), 2 = downscaled "
                        "CPU twins (the CI gate).")
    p.add_argument("--backend", choices=("cpu", "gpu"), default=None,
                   help="default: cpu for tier 2, gpu for tier 1.")
    p.add_argument("--build-dir", default=None,
                   help="default: build_gcc for cpu, build for gpu.")
    p.add_argument("--binary", default=None)
    p.add_argument("--gpus", default=None, help="e.g. 0,1,2,3")
    p.add_argument("--jobs", type=int, default=1,
                   help="parallel CPU workers (GPU uses one worker per device).")
    p.add_argument("--cases", default=None, help="comma-separated name filter.")
    p.add_argument("--tags", default=None,
                   help="comma-separated TAG filter, e.g. `vcoord_matrix` to "
                        "run only the vertical-coordinate rest matrix. "
                        "Composes with --cases (both must match).")
    p.add_argument("--out", default=None, help="write JSON results here.")
    p.add_argument("--keep", action="store_true", help="keep NetCDF output.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="print every assertion, not only failures.")
    p.add_argument("--self-test", action="store_true",
                   help="verify the ASSERTIONS themselves against the measured "
                        "2026-09-11 signatures, then exit. Runs no model.")
    p.add_argument("--skip-rule-check", action="store_true",
                   help="do not validate the tier-2 downscale rules first.")
    p.add_argument("--split-scheme", default=None,
                   choices=("ssp_rk2", "pred_corr"),
                   help="force &ocean_bt_nml split_scheme on every case that "
                        "does not already pin one, for an A/B sweep of the "
                        "outer time-splitting. NOT a substitute for the "
                        "permanent scheme axis in the manifest -- it is the "
                        "tool for measuring a candidate default before "
                        "committing to it.")
    args = p.parse_args(argv)
    _VERBOSE = args.verbose
    if args.self_test:
        print("stability.py self-test -- no model is run\n")
        return self_test()

    print("split scheme  : base cases run the &ocean_bt_nml DEFAULT "
          "(\"{}\"); `<case>__{}` twins force the other."
          .format(_DEFAULT_SPLIT_SCHEME, _stab_manifest.SCHEME_AXIS_SCHEME))
    backend = args.backend or ("cpu" if args.tier == 2 else "gpu")
    build_dir = args.build_dir or ("build_gcc" if backend == "cpu" else "build")
    if not os.path.isabs(build_dir):
        build_dir = os.path.join(REPO_ROOT, build_dir)
    try:
        binary = rr.locate_binary(build_dir, args.binary)
    except FileNotFoundError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    only = {c.strip() for c in args.cases.split(",")} if args.cases else None
    tags = {t.strip() for t in args.tags.split(",")} if args.tags else None
    cases = select_cases(args.tier, only, tags)
    if args.split_scheme:
        cases = [_force_scheme(c, args.split_scheme) for c in cases]
    gpus = ([g.strip() for g in args.gpus.split(",")]
            if (backend == "gpu" and args.gpus) else
            (rr.detect_gpus(None) if backend == "gpu" else None))

    # A forced-scheme A/B sweep gets its OWN scratch root: two sweeps sharing
    # one would have the second overwrite the run logs of the first, which is
    # exactly the comparison the sweep exists to make.
    scratch_root = os.path.join(
        REPO_ROOT, "tmp_local_artifacts",
        "stability" + ("_" + args.split_scheme if args.split_scheme else ""))
    os.makedirs(scratch_root, exist_ok=True)

    print("=" * 92)
    print("Roundabout ocean stability suite -- TIER {} ({})".format(
        args.tier,
        "full-scale, GPU, local/nightly" if args.tier == 1
        else "downscaled twins, CPU, the CI gate"))
    print("  binary : {}".format(binary))
    print("  cases  : {} of {} in the manifest".format(
        len(cases), len(manifest.STABILITY_CASES)))
    if gpus:
        print("  gpus   : {}".format(",".join(gpus)))
    print("=" * 92)

    # Tier-2 twins must satisfy the dimensionless-number rules BEFORE they run:
    # a violating twin is not the same test, however green it looks.
    rule_fail = 0
    if args.tier == 2 and not args.skip_rule_check:
        for c in cases:
            spec = (c["tier2"].get("dimensionless") or {})
            if not spec:
                continue
            res = downscale.check_twin(spec)
            if not downscale.twin_ok(res):
                rule_fail += 1
                print("\nDOWNSCALE RULE VIOLATION -- this twin is not the same "
                      "test as its tier-1 parent:")
                print(downscale.format_report(c["name"], res))
        if rule_fail:
            print("\n{} twin(s) violate the downscaling rules (downscale.py). "
                  "Fix the manifest; do not run.".format(rule_fail))
            return 3

    t0 = time.time()
    results = run_all(cases, args.tier, binary, scratch_root, gpus,
                      args.keep, args.jobs)
    wall = time.time() - t0

    order = [c["name"] for c in cases]
    rows = []
    for n in order:
        if n not in results:
            continue
        c, res, verdicts, artifacts = results[n]
        status, failed = _status(c, verdicts, args.tier)
        rows.append((n, status, res["wall_s"],
                     ";".join(v.name for v in failed), artifacts))

    print("\n" + "=" * 92)
    print("TIER {} SUMMARY   ({} cases, {:.1f} s wall)".format(
        args.tier, len(rows), wall))
    print("-" * 92)
    print("{:<32} {:<7} {:>8}  {}".format("case", "status", "wall_s",
                                          "failed assertions"))
    for n, status, w, f, _a in rows:
        print("{:<32} {:<7} {:>8.1f}  {}".format(n, status, w, f))
    counts = {}
    for _n, s, _w, _f, _a in rows:
        counts[s] = counts.get(s, 0) + 1
    print("-" * 92)
    print("  " + "  ".join("{} {}".format(counts.get(k, 0), k)
                           for k in ("PASS", "FAIL", "XFAIL", "XPASS")))
    art = sorted({a for _n, _s, _w, _f, arts in rows for a in arts})
    if art:
        print("  ARTIFACT (not a failure): diag mean NaN on [{}]".format(
            ", ".join(art)))
        print("           " + "\n           ".join(_wrap(DIAG_MEAN_NAN_NOTE, 78)))
    print("=" * 92)

    matrix_rows = _matrix_summary(results, args.tier)
    if matrix_rows:
        _print_matrix_table(matrix_rows, args.tier)

    if args.out:
        payload = []
        for n in order:
            if n not in results:
                continue
            c, res, verdicts, artifacts = results[n]
            status, _f = _status(c, verdicts, args.tier)
            payload.append({
                "case": n, "tier": args.tier, "status": status,
                "wall_s": res["wall_s"], "nml": c["nml"],
                "artifacts": artifacts,
                "assertions": [{"name": v.name, "ok": v.ok,
                                "skipped": v.skipped,
                                "expected": v.expected,
                                "observed": v.observed} for v in verdicts],
            })
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"tier": args.tier, "wall_s": wall,
                       "cases": payload, "vcoord_matrix": matrix_rows}, fh,
                      indent=1)
        print("JSON: {}".format(args.out))
        if matrix_rows:
            mpath = os.path.splitext(args.out)[0] + "_vcoord_matrix.json"
            with open(mpath, "w") as fh:
                json.dump({"tier": args.tier, "rows": matrix_rows}, fh, indent=1)
            print("JSON: {}  (the family x problem baseline table)".format(mpath))

    return 1 if counts.get("FAIL") else 0


# ---------------------------------------------------------------------------
# The vertical-coordinate matrix summary
# ---------------------------------------------------------------------------
def _matrix_summary(results, tier):
    """One machine-readable record per matrix cell that ran.

    Written next to the ordinary report so drift in the family x problem
    table can be tracked over time -- the table IS the product of the
    matrix, and a table that only ever appears in a terminal cannot be
    diffed against last week's.
    """
    rows = []
    for name, (case, res, verdicts, _art) in sorted(results.items()):
        mx = case.get("matrix")
        if not mx:
            continue
        status, failed = _status(case, verdicts, tier)
        by = {v.name: v for v in verdicts}

        def observed(key):
            v = by.get(key)
            return v.observed if v else None

        rows.append({
            "case": name,
            "problem": mx["problem"], "family": mx["family"],
            "stratification": mx["stratification"], "eos": mx["eos"],
            "scheme": case.get("split_scheme") or _DEFAULT_SPLIT_SCHEME,
            "expect": mx["expect"],
            "status": status,
            "failed": [v.name for v in failed],
            "wall_s": res.get("wall_s"),
            "rx0": res.get("rx0"),
            "ntrunc": res.get("ntrunc"),
            "en_rest": observed("energy:rest"),
            "growth_rate": observed("energy:rest-growth-rate"),
            "tracer": observed("tracer:no-new-extrema"),
            "thickness": observed("thickness:positive"),
            "truncation_estimate": mx.get("en_estimate"),
        })
    return rows


def _print_matrix_table(rows, tier):
    """The family x problem BASELINE TABLE, printed at the end of a run.

    This is what a user reads to decide which coordinate to trust on which
    geometry. One cell per (problem, family): the verdict, and the key
    number behind it.
    """
    problems = sorted({r["problem"] for r in rows})
    families = sorted({r["family"] for r in rows})
    cell = {}
    for r in rows:
        key = (r["problem"], r["family"])
        mark = {"PASS": "ok", "FAIL": "FAIL", "XFAIL": "xfail",
                "XPASS": "XPASS"}[r["status"]]
        if r["expect"] == "refused":
            mark = "refused" if r["status"] == "XFAIL" else "NOT-REFUSED"
        # A twin lands in the same cell; keep the worse verdict.
        prev = cell.get(key)
        rank = {"ok": 0, "xfail": 1, "refused": 0, "FAIL": 3, "XPASS": 2,
                "NOT-REFUSED": 3}
        if prev is None or rank[mark] > rank[prev]:
            cell[key] = mark
    w = max(12, max(len(f) for f in families) + 1)
    print("\n" + "=" * 92)
    print("VERTICAL-COORDINATE MATRIX -- tier {}   (family x problem)".format(tier))
    print("  ok = every assertion passed   xfail = documented defect   "
          "refused = rejected at configure, as designed")
    print("  FAIL = a live defect          XPASS/NOT-REFUSED = the envelope "
          "moved; re-measure")
    print("-" * 92)
    head = "{:<18}".format("problem") + "".join(
        "{:<{w}}".format(f[:w - 1], w=w) for f in families)
    print(head)
    for p in problems:
        line = "{:<18}".format(p) + "".join(
            "{:<{w}}".format(cell.get((p, f), "-"), w=w) for f in families)
        print(line)
    print("=" * 92)


if __name__ == "__main__":
    sys.exit(main())
