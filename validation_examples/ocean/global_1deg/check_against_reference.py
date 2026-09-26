#!/usr/bin/env python3
"""Check a global 1-degree run against the committed reference — standard library only.

    python3 check_against_reference.py RUN.log [--reference reference_daily.csv]
                                               [--days N] [--quick]

``RUN.log`` is the console output of a run of ``global_1deg_unforced.nml``:

* the ``rdb`` executable prints, every simulated day, the
  ``[stats] Day ... En ... MaxCFL ...`` line and the ``Mass`` / ``Salt`` /
  ``Heat`` budget lines — all of it is checked;
* ``run_global_1deg.py`` (the Python interface) prints its own
  ``[py] day D En mass_drift ...`` line.  The Python API exposes the
  kinetic energy and the total mass but not MaxCFL or the salt / heat
  totals, so a Python-driven log is checked on En and the mass residual
  only (the rest shows as n/a).  The two runs are bit-identical, so the
  executable's check covers the rest.

Every day the log carries is compared with the same day of
``reference_daily.csv``; ``--days N`` additionally requires the log to reach
day N (a run that died early FAILS).  ``--quick`` is ``--days 10``: the
ten-day slice a new machine can run in about two minutes.

What is checked, and why these tolerances
-----------------------------------------

The reference was measured with one toolchain (nvfortran 26.5 on a V100,
see the CSV header).  Another compiler, another GPU, or the CPU build does
not produce the same bits: FMA contraction, reduction order and libm
differ at round-off, and round-off grows.  So the two kinds of number are
checked differently.

**Energy (En, and MaxCFL): a relative band that widens with time.**
En is a global integral of a smooth, large-scale adjustment (the
barotropic mode answering the WOA pressure field within days, a plateau
through the first hundred days, then the closures spinning it down), so it stays close across toolchains long after
individual features have decorrelated — but it does decorrelate, so the
band widens:

====================  ===========  ==========================================
days                   En band      reasoning
====================  ===========  ==========================================
1 - 10  (quick mode)   0.5 %        deterministic adjustment (En peaks on
                                    day 3); the console prints En to 4
                                    digits (rounding alone is up to 0.2 %),
                                    cross-toolchain round-off is still far
                                    below that
11 - 30                2 %          spin-up; the one-cell straits and
                                    shelves that carry the fastest water
                                    start to decorrelate
31 - 90                5 %          En plateau (secondary maximum on day
                                    63); eddying part of the flow
                                    decorrelated
91 - 365               10 %         slow spin-down of a decorrelated flow; the
                                    global integral is still pinned by the
                                    initial state and the closures
====================  ===========  ==========================================

On the reference toolchain the run is deterministic and matches every
printed digit (the previous reference year re-ran all 365 days identically),
so the band is slack there; it exists for the others.  On that previous
reference, gfortran 15.1 on the CPU (serial) also matched En and MaxCFL to
every printed digit over the 10 quick days (not yet repeated on this one).
The later bands are not measured across toolchains (a CPU year is days of
wall time); they are set from how the flow evolves, and are the numbers to
revisit when a second toolchain's year is on record.
MaxCFL is a pointwise maximum, far more sensitive to where one fast cell
sits than the integral, so its band is twice En's (and it must stay
below 0.5 every day, well above the reference maximum of 0.243).

**Budgets (Mass, Salt, Heat Error): round-off in THIS run — not equality
with the reference.**  The Error is the relative closure residual of a
closed domain; it is pure accumulated round-off, and round-off is
toolchain-specific, so comparing its digits with the reference would be
meaningless.  Instead each series must

1. stay below a round-off envelope ``|Error(d)| <= floor + rate * d``.
   Two round-off sources add.  (a) Per-step round-off accumulates:
   ``rate`` = 1e-13 per day for mass and 1e-14 per day for salt and heat —
   with 48 steps a day ~2e-15 and ~2e-16 per step, a few machine epsilons;
   the reference sits at 1.8e-14 (mass) and 4e-16 (salt, heat) per day, 5x
   and 25x inside.  (b) The global totals are themselves double sums over
   5.8 M cells, and HOW they are summed is toolchain-specific: the GPU's
   tree reductions are nearly exact, a serial CPU loop is not (its error
   grows like eps * sqrt(N)).  Measured with gfortran 15.1 (serial CPU,
   --quick): mass 1.30e-12 - 1.46e-12, salt and heat jittering within
   +-2.1e-13 over days 1 - 10 — a flat measurement floor, not a drift,
   while En and MaxCFL match the reference to every printed digit.
   ``floor`` = 1e-11 covers it 7x.
   A real leak — a missing boundary flux, a non-conservative remap — is
   1e-10 or more within a day and fails on day 1.
2. grow no faster than linearly: the mean daily increment over the second
   half of the run may not exceed 4x the mean daily increment over the
   first half, plus 1e-12 per day for the day-to-day jitter of the
   summation floor.  Round-off accumulates at a steady rate; an
   accelerating residual is a leak that feeds on the flow.

Every value must also be finite.  Exit status 0 = PASS, 1 = FAIL.
"""

import argparse
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REFERENCE = os.path.join(HERE, "reference_daily.csv")
QUICK_DAYS = 10

# (last day, relative En band); the MaxCFL band is CFL_FACTOR times this.
EN_BANDS = [(10, 0.005), (30, 0.02), (90, 0.05), (10 ** 9, 0.10)]
CFL_FACTOR = 2.0
CFL_CEILING = 0.5
# Round-off envelope of the relative residual: BUDGET_FLOOR + BUDGET_RATE * day.
BUDGET_RATE = {"mass": 1.0e-13, "salt": 1.0e-14, "heat": 1.0e-14}
BUDGET_FLOOR = 1.0e-11
LINEARITY_FACTOR = 4.0
LINEARITY_FLOOR = 1.0e-12
NAN = float("nan")

_NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?|NaN|nan|Infinity|-?Inf)"
STATS_RE = re.compile(r"\[stats\]\s+Day\s+" + _NUM + r"\s+step\s+(\d+)\s+En\s+" + _NUM
                      + r"\s+MaxCFL\s+" + _NUM)
BUDGET_RE = re.compile(r"\b(Mass|Salt|Heat)\s*:\s*\S+\s+Error\s+" + _NUM)
# run_global_1deg.py: "[py] day  D  En  (mass - mass0)/mass0  max_surface_speed  wall ..."
PY_RE = re.compile(r"\[py\] day\s+(\d+)\s+" + _NUM + r"\s+" + _NUM)


def num(s):
    s = s.replace("D", "E").replace("d", "e")
    try:
        return float(s)
    except ValueError:
        return NAN


def en_band(day):
    for last, band in EN_BANDS:
        if day <= last:
            return band
    return EN_BANDS[-1][1]


def parse_log(path):
    """{day: {"En", "MaxCFL", "mass", "salt", "heat"}} for whole days > 0.
    Quantities a log does not carry are absent (the Python driver's log)."""
    days, cur = {}, None
    with open(path, errors="replace") as f:
        for line in f:
            m = PY_RE.search(line)
            if m:
                days.setdefault(int(m.group(1)), {"En": num(m.group(2)),
                                                   "mass": num(m.group(3))})
                continue
            m = STATS_RE.search(line)
            if m:
                d = num(m.group(1))
                if d > 0.0 and abs(d - round(d)) < 1e-6:
                    cur = int(round(d))
                    days[cur] = {"En": num(m.group(3)), "MaxCFL": num(m.group(4))}
                else:
                    cur = None
                continue
            m = BUDGET_RE.search(line)
            if m and cur is not None:
                days[cur].setdefault(m.group(1).lower(), num(m.group(2)))
    return days


def read_reference(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            rows[int(r["day"])] = {k: float(v) for k, v in r.items() if k != "day"}
    return rows


def rel(a, b):
    return abs(a - b) / abs(b) if b else abs(a - b)


def check(run, ref, want_days):
    """-> (failures, table rows, budget quantities present in the log)."""
    failures = []
    days = sorted(d for d in run if d in ref)
    if not days:
        return ["no daily [stats] / [py] line of the log matches a reference day"], [], []
    last = days[-1]
    if want_days and last < want_days:
        failures.append(f"log ends at day {last}, expected day {want_days}")
    missing = [d for d in range(1, last + 1) if d not in run]
    if missing:
        failures.append(f"days missing from the log: {missing[:10]}"
                        f"{' ...' if len(missing) > 10 else ''}")
    # The executable's log carries all three budgets, the Python driver's only
    # the mass; whatever day 1 carries is then required every day.
    quantities = [q for q in ("mass", "salt", "heat") if q in run[days[0]]]
    if "mass" not in quantities:
        failures.append("the log carries no mass budget")

    table = []
    for d in days:
        r, x = run[d], ref[d]
        band = en_band(d)
        cfl = r.get("MaxCFL")
        row = {"day": d, "En": r["En"], "En_ref": x["En"], "dEn": rel(r["En"], x["En"]),
               "band": band, "MaxCFL": NAN if cfl is None else cfl,
               "MaxCFL_ref": x["MaxCFL"], "flags": []}
        for q in ("mass", "salt", "heat"):
            row[q] = r.get(q, NAN)
            if q not in quantities:
                continue
            if not math.isfinite(row[q]):
                row["flags"].append(f"{q} Error missing or non-finite")
            elif abs(row[q]) > BUDGET_FLOOR + BUDGET_RATE[q] * d:
                row["flags"].append(f"|{q} Error| {abs(row[q]):.2e} > "
                                    f"{BUDGET_FLOOR + BUDGET_RATE[q] * d:.2e}")
        if not math.isfinite(r["En"]):
            row["flags"].append("non-finite En")
        elif row["dEn"] > band:
            row["flags"].append(f"En off by {100 * row['dEn']:.2f}% > {100 * band:g}%")
        if cfl is not None:
            dcfl = rel(cfl, x["MaxCFL"])
            if not math.isfinite(cfl):
                row["flags"].append("non-finite MaxCFL")
            elif dcfl > CFL_FACTOR * band:
                row["flags"].append(f"MaxCFL off by {100 * dcfl:.2f}% > "
                                    f"{100 * CFL_FACTOR * band:g}%")
            elif cfl >= CFL_CEILING:
                row["flags"].append(f"MaxCFL {cfl} >= {CFL_CEILING}")
        failures += [f"day {d}: {msg}" for msg in row["flags"]]
        table.append(row)

    # Linearity: no acceleration of the round-off residual.
    if len(days) >= 4:
        for q in quantities:
            series = [0.0] + [run.get(d, {}).get(q, NAN) for d in range(1, last + 1)]
            if not all(math.isfinite(v) for v in series):
                continue
            inc = [abs(series[i] - series[i - 1]) for i in range(1, len(series))]
            h = len(inc) // 2
            early, late = sum(inc[:h]) / h, sum(inc[h:]) / (len(inc) - h)
            if late > LINEARITY_FACTOR * early + LINEARITY_FLOOR:
                failures.append(f"{q} Error accelerates: mean daily increment {late:.2e} in the "
                                f"second half vs {early:.2e} in the first (> {LINEARITY_FACTOR:g}x)")
    return failures, table, quantities


def fmt(v, spec, width):
    return f"{v:{width}{spec}}" if math.isfinite(v) else f"{'n/a':>{width}}"


def print_table(table, every):
    print(f"{'day':>4} {'En':>10} {'En ref':>10} {'dEn':>7} {'band':>5} "
          f"{'MaxCFL':>8} {'ref':>8} {'mass Err':>10} {'salt Err':>10} {'heat Err':>10}  status")
    n = len(table)
    for i, r in enumerate(table):
        if r["flags"] or i < 10 or i == n - 1 or r["day"] % every == 0:
            print(f"{r['day']:4d} {r['En']:10.3e} {r['En_ref']:10.3e} {100 * r['dEn']:6.2f}% "
                  f"{100 * r['band']:4g}% {fmt(r['MaxCFL'], '.5f', 8)} {r['MaxCFL_ref']:8.5f} "
                  f"{fmt(r['mass'], '.2e', 10)} {fmt(r['salt'], '.2e', 10)} "
                  f"{fmt(r['heat'], '.2e', 10)}  {'FAIL' if r['flags'] else 'ok'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("log", help="console log of the run")
    ap.add_argument("--reference", default=DEFAULT_REFERENCE)
    ap.add_argument("--days", type=int, default=0,
                    help="the log must reach this day (0: whatever it reaches)")
    ap.add_argument("--quick", action="store_true", help=f"= --days {QUICK_DAYS}")
    ap.add_argument("--every", type=int, default=30, help="table row stride after day 10")
    a = ap.parse_args()
    want = QUICK_DAYS if a.quick else a.days

    run = parse_log(a.log)
    ref = read_reference(a.reference)
    if want and want > max(ref):
        sys.exit(f"--days {want}: the reference only reaches day {max(ref)}")
    run = {d: v for d, v in run.items() if not want or d <= want}
    failures, table, quantities = check(run, ref, want)
    print(f"log:       {os.path.abspath(a.log)}")
    print(f"reference: {os.path.abspath(a.reference)}")
    if table:
        print_table(table, a.every)
    if failures:
        print(f"\nFAIL ({len(failures)} problem{'s' if len(failures) > 1 else ''}):")
        for msg in failures[:40]:
            print(f"  - {msg}")
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
        return 1
    worst = max(table, key=lambda r: r["dEn"] / r["band"])
    absent = [q for q in ("MaxCFL",) if not any(math.isfinite(r["MaxCFL"]) for r in table)]
    absent += [q for q in ("salt", "heat") if q not in quantities]
    print(f"\nPASS: {len(table)} days (1-{table[-1]['day']}); largest En deviation "
          f"{100 * worst['dEn']:.3f}% on day {worst['day']} (band {100 * worst['band']:g}%); "
          f"{'/'.join(quantities)} Error within the round-off envelope and linear"
          + (f"; not in this log: {', '.join(absent)}" if absent else "") + ".")
    return 0


if __name__ == "__main__":
    sys.exit(main())
