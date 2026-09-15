#!/usr/bin/env python3
"""Parse and diff the CHKSUM / HOTFACE forensic rows emitted by the ocean
`&ocean_debug_nml chksum` probe (`src/core/ocean/diag/rdb_ocean_chksum.F90`).

The probe writes fixed-format greppable rows straight to stdout (bypassing the
logger by design) at every split-RK2 tendency seam:

    CHKSUM  <step> s<stage> <phase> <field> <sum> <min> <max> <nonfin> <bits>
    HOTFACE <step> s<stage> <phase> u|v (i,j,k) <val> <hL> <hR> <rem> <hdn> <hup>

Two uses:

  * Attribution — the FIRST seam whose `nonfin > 0` MINTED the corruption
    (later seams only inherit it).  Grep or `--first-nonfin`.
  * Decomposition regression — the `bits` column is a POPCNT reduction of the
    IEEE-754 bit pattern, summed as int64.  Integer addition is exactly
    associative + commutative, so `bits` is IDENTICAL across rank counts and
    loop orders (unlike the FP `sum`, which only reproduces within one binary).
    `--diff LOG_A LOG_B` aligns rows by (step, stage, phase, field) and reports
    the FIRST diverging row: a `bits` mismatch is the authoritative
    "operator X diverged" signal; `sum` is compared only within a tiny FP
    tolerance as a secondary hint.

Robustness: Fortran `es20.12` can emit 3-digit exponents (`1.234E+100`,
`-5.6E-105`) and, on some compilers, drop the `E` entirely (`1.234-105`).
Both forms are handled.

Standard library only (hard repo rule — no pip, no third-party imports).
"""
import argparse
import math
import re
import sys


# A Fortran real: optional sign, mantissa, optional exponent.  The exponent
# may be `E+NN` / `E-NNN`, a bare `+NN`, or — on compilers that overflow the
# `Ew.d` exponent field — a sign glued to the mantissa with no `E`
# (`1.234-105`).  Capture the whole token, normalise in `_to_float`.
_FLOAT_TOKEN = r"[+-]?\d+\.\d+(?:[EeDd][+-]?\d+|[+-]\d{2,3})?"
_NANINF = r"(?:[+-]?(?:NaN|Inf(?:inity)?))"
_NUM = rf"(?:{_FLOAT_TOKEN}|{_NANINF})"

# CHKSUM row.  `bits` is optional so pre-PR-A logs (no bits column) still parse.
_CHKSUM_RE = re.compile(
    r"^\s*CHKSUM\s+(?P<step>\d+)\s+s(?P<stage>\d+)\s+"
    r"(?P<phase>\S+)\s+(?P<field>\S+)\s+"
    rf"(?P<sum>{_NUM})\s+(?P<min>{_NUM})\s+(?P<max>{_NUM})\s+"
    r"(?P<nonfin>\d+)(?:\s+(?P<bits>\d+))?\s*$"
)

_HOTFACE_RE = re.compile(
    r"^\s*HOTFACE\s+(?P<step>\d+)\s+s(?P<stage>\d+)\s+"
    r"(?P<phase>\S+)\s+(?P<comp>[uv])\s+"
    r"\(\s*(?P<i>\d+)\s*,\s*(?P<j>\d+)\s*,\s*(?P<k>\d+)\s*\)\s+"
    r"(?P<rest>.*\S)\s*$"
)


def _to_float(token):
    """Parse a Fortran real token, repairing the missing-`E` exponent form."""
    t = token.strip().replace("D", "E").replace("d", "E")
    low = t.lower()
    if "nan" in low:
        return math.nan
    if "inf" in low:
        return -math.inf if t.lstrip().startswith("-") else math.inf
    if "e" not in low:
        # Possibly `1.234-105` (exponent sign glued to mantissa with no E).
        m = re.match(r"^([+-]?\d+\.\d+)([+-]\d{2,3})$", t)
        if m:
            t = m.group(1) + "E" + m.group(2)
    return float(t)


def parse_log(path):
    """Yield parsed records from a run log.

    Each record is a dict with a `kind` ('chksum'|'hotface') and a `key`
    tuple (step, stage, phase, field) for alignment.
    """
    records = []
    with open(path, "r", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if "CHKSUM" in line:
                m = _CHKSUM_RE.match(line)
                if not m:
                    continue
                g = m.groupdict()
                records.append(
                    {
                        "kind": "chksum",
                        "lineno": lineno,
                        "key": (int(g["step"]), int(g["stage"]), g["phase"], g["field"]),
                        "sum": _to_float(g["sum"]),
                        "min": _to_float(g["min"]),
                        "max": _to_float(g["max"]),
                        "nonfin": int(g["nonfin"]),
                        "bits": int(g["bits"]) if g["bits"] is not None else None,
                    }
                )
            elif "HOTFACE" in line:
                m = _HOTFACE_RE.match(line)
                if not m:
                    continue
                g = m.groupdict()
                records.append(
                    {
                        "kind": "hotface",
                        "lineno": lineno,
                        "key": (int(g["step"]), int(g["stage"]), g["phase"], g["comp"]),
                        "ijk": (int(g["i"]), int(g["j"]), int(g["k"])),
                        "vals": [_to_float(tok) for tok in g["rest"].split()],
                    }
                )
    return records


def _fmt_key(key):
    step, stage, phase, field = key
    return f"step={step} s{stage} {phase}/{field}"


def cmd_dump(args):
    recs = parse_log(args.log)
    for r in recs:
        if args.first_nonfin and not (r["kind"] == "chksum" and r["nonfin"] > 0):
            continue
        if r["kind"] == "chksum":
            print(
                f"{_fmt_key(r['key'])}  sum={r['sum']:.12e}  "
                f"min={r['min']:.6e}  max={r['max']:.6e}  "
                f"nonfin={r['nonfin']}  bits={r['bits']}"
            )
            if args.first_nonfin:
                print("  ^ FIRST non-finite-minting seam", file=sys.stderr)
                return 0
        else:
            print(f"{_fmt_key(r['key'])}  ijk={r['ijk']}  vals={r['vals']}")
    if args.first_nonfin:
        print("no non-finite seam found", file=sys.stderr)
    return 0


def cmd_diff(args):
    recs_a = [r for r in parse_log(args.log_a) if r["kind"] == "chksum"]
    recs_b = [r for r in parse_log(args.log_b) if r["kind"] == "chksum"]
    by_key_b = {}
    for r in recs_b:
        by_key_b.setdefault(r["key"], r)

    checked = 0
    for ra in recs_a:
        rb = by_key_b.get(ra["key"])
        if rb is None:
            continue
        checked += 1
        # bits is the authoritative signal — exact int compare.
        if ra["bits"] is not None and rb["bits"] is not None and ra["bits"] != rb["bits"]:
            print(f"DIVERGE (bits) at {_fmt_key(ra['key'])}")
            print(f"  A: bits={ra['bits']}  sum={ra['sum']:.12e}")
            print(f"  B: bits={rb['bits']}  sum={rb['sum']:.12e}")
            print("  -> a bits mismatch names the first operator that diverged "
                  "(decomposition-invariant).")
            return 1
        # Secondary hint when bits is absent (old logs) or matches: FP sum.
        sa, sb = ra["sum"], rb["sum"]
        finite = math.isfinite(sa) and math.isfinite(sb)
        if not finite and (math.isnan(sa) != math.isnan(sb)
                           or math.isinf(sa) != math.isinf(sb)):
            print(f"DIVERGE (sum non-finite) at {_fmt_key(ra['key'])}")
            print(f"  A: sum={sa}   B: sum={sb}")
            return 1
        if finite:
            scale = max(1.0, abs(sa), abs(sb))
            if abs(sa - sb) > args.tol * scale:
                print(f"DIVERGE (sum) at {_fmt_key(ra['key'])}")
                print(f"  A: sum={sa:.12e}   B: sum={sb:.12e}   "
                      f"|delta|={abs(sa - sb):.3e} > tol*{scale:.3e}")
                return 1

    if checked == 0:
        print("no common (step,stage,phase,field) CHKSUM keys between the two logs",
              file=sys.stderr)
        return 2
    print(f"no divergence over {checked} aligned CHKSUM rows")
    return 0


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode")

    d = sub.add_parser("dump", help="parse and print CHKSUM/HOTFACE rows from one log")
    d.add_argument("log")
    d.add_argument("--first-nonfin", action="store_true",
                   help="print only the FIRST seam with nonfin > 0 (the minter)")
    d.set_defaults(func=cmd_dump)

    df = sub.add_parser("diff", help="align two logs by key, report the first diverging row")
    df.add_argument("log_a")
    df.add_argument("log_b")
    df.add_argument("--tol", type=float, default=1e-12,
                    help="relative FP tolerance for the secondary `sum` compare "
                         "(default 1e-12); `bits` is always compared exactly")
    df.set_defaults(func=cmd_diff)
    return p


def main(argv=None):
    parser = build_parser()
    # Convenience: `--diff A B` is an alias for the `diff` subcommand so the
    # brief's `read_chksum.py --diff LOG_A LOG_B` invocation works verbatim.
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "--diff":
        argv = ["diff"] + argv[1:]
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
