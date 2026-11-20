#!/usr/bin/env python3
"""Prove that a TOML input resolves to the same configuration as its namelist.

    tools/check_nml_toml_equiv.py --build-dir build_gfortran OLD.nml NEW.toml [...]

For each OLD/NEW pair: run `rdb_cfg_dump` on both and require BOTH outputs
(the doc dump and the 18-digit `.exact` dump) to be byte-identical, then run
`rdb --validate-only` on both from NEW's directory (where the case's relative
input files resolve) and require equal exit codes.  Prints one row per pair;
exit 0 iff every pair matches.  Scratch goes under `tmp_local_artifacts/`.
"""
import argparse
import filecmp
import os
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def run(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=300)
    return p.returncode, p.stdout + p.stderr


def check_pair(build, old, new, scratch):
    dump = os.path.join(build, "rdb_cfg_dump")
    rdb = os.path.join(build, "rdb")
    tag = new.replace(os.sep, "_")
    a, b = os.path.join(scratch, tag + ".old"), os.path.join(scratch, tag + ".new")
    rc_a, _ = run([dump, old, a])
    rc_b, _ = run([dump, new, b])
    if rc_a != rc_b:
        return False, "cfg_dump rc {} vs {}".format(rc_a, rc_b)
    same = True
    if rc_a == 0:
        same = filecmp.cmp(a, b, shallow=False) and \
            filecmp.cmp(a + ".exact", b + ".exact", shallow=False)
    cwd = os.path.dirname(os.path.abspath(new))
    v_a = run([rdb, "--validate-only", os.path.abspath(old)], cwd=cwd)[0]
    v_b = run([rdb, "--validate-only", os.path.abspath(new)], cwd=cwd)[0]
    ok = same and v_a == v_b
    return ok, "dump rc={} identical={} validate rc {} vs {}".format(
        rc_a, same, v_a, v_b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True)
    ap.add_argument("pairs", nargs="+", help="OLD NEW [OLD NEW ...]")
    args = ap.parse_args()
    if len(args.pairs) % 2:
        sys.exit("pairs must be OLD NEW ...")
    scratch = os.path.join(REPO, "tmp_local_artifacts", "nml_toml_equiv")
    os.makedirs(scratch, exist_ok=True)
    build = os.path.abspath(args.build_dir)
    bad = 0
    for old, new in zip(args.pairs[::2], args.pairs[1::2]):
        ok, msg = check_pair(build, os.path.abspath(old), new, scratch)
        print("{}  {}  {}".format("OK  " if ok else "FAIL", new, msg), flush=True)
        bad += not ok
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
