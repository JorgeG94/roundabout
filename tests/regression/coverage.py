#!/usr/bin/env python3
"""P2 gcov coverage mode -- measure what fraction of the OCEAN physics/closure
source the regression corpus exercises, and surface the biggest UNCOVERED
closure files as the P4 gap list.

What it does (end to end)
-------------------------
1. **Build a coverage binary.** Configure + build `<build-dir>` (default
   `build_cov`) with the existing `RDB_ENABLE_COVERAGE=ON` option (gfortran
   `-O0 -g --coverage`), GPU + MPI off. If the binary already exists it is
   reused (pass --rebuild to force). Load your gfortran / NetCDF toolchain
   first (module, Spack -- see environments/ --, conda) so it is on PATH; this
   script inherits that environment. Do not stack a second toolchain in the
   same shell -- two NetCDF builds on one link line fail in confusing ways.
   (Run via `run_all.py` on a configured dev box with RDB_ON_DEV=1 and it
   sources the local site env script for you; plain invocation never does.)

2. **Run the ocean manifest corpus** against the coverage binary, SERIALLY
   (coverage instrumentation is slow; that is fine -- it has its own looser
   budget). This reuses `run_regression.run_case` verbatim: the same temp-nml
   run-length patch + output isolation to a scratch dir. The `.gcda` counters
   accumulate across every case (counters are zeroed once up front).

3. **Capture coverage** two ways:
   - Always: a **stdlib-only gcov fallback** -- run `gcov` on every
     instrumented object's `.gcno`/`.gcda` under the build dir and parse the
     per-line counts out of the emitted `.gcov` files (exact, version-stable).
   - If `lcov` + `genhtml` are on PATH: additionally render the HTML report
     (like the CMake `coverage` target) and parse lcov's per-file numbers as a
     cross-check. lcov is NOT required -- the env is flaky, so the gcov numbers
     are authoritative.

4. **Report**, focused on the physics/closures:
   - overall line coverage (whole instrumented tree + ocean-closure subset),
   - a per-file table for the OCEAN CLOSURE sources, sorted by coverage ascending,
   - an explicit GAPS list: closure files at 0% / very low coverage -- the
     physics the corpus does NOT exercise (the P4 case-addition targets),
   - a machine-readable `coverage_summary.json` plus a human table to stdout.

Ocean-closure source set
------------------------
Files under `src/core/ocean/**`, the ocean parameterizations
(`src/parameterizations/{vertical,lateral}/structured/rdb_ocean_*`), and the
ocean-specific shared-physics files (`src/{equation_of_state,pressure_force,
ALE,tracer}/**` whose basename contains `ocean`).

Stdlib only (the repo forbids `pip install`); no third-party imports.

Run `python3 coverage.py --help` for the full flag list.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# Paths / reuse the P0 runner machinery (per-case temp-nml patch + isolation).
# ---------------------------------------------------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
sys.path.insert(0, THIS_DIR)
import manifest          # noqa: E402
import run_regression    # noqa: E402  (reuse run_case / locate_binary verbatim)

DEFAULT_BUILD_DIR = "build_cov"
DEFAULT_BUDGET_MIN = 60.0   # looser than the correctness pass; instrumented = slow


# ---------------------------------------------------------------------------
# Ocean-closure source classification
# ---------------------------------------------------------------------------
def is_ocean_closure(relpath):
    """True if `relpath` (repo-relative, forward slashes) is an ocean physics /
    closure source file we want in the per-closure report."""
    p = relpath.replace(os.sep, "/")
    if not (p.endswith(".F90") or p.endswith(".f90")):
        return False
    base = os.path.basename(p).lower()
    # (a) everything under the ocean dyn-core tree
    if p.startswith("src/core/ocean/"):
        return True
    # (b) ocean parameterizations (vertical + lateral, structured)
    if (p.startswith("src/parameterizations/vertical/structured/")
            or p.startswith("src/parameterizations/lateral/structured/")):
        return base.startswith("rdb_ocean_")
    # (c) ocean-specific shared physics (eos / pgf / ALE / tracer)
    for root in ("src/equation_of_state/", "src/pressure_force/",
                 "src/ALE/", "src/tracer/"):
        if p.startswith(root) and "ocean" in base:
            return True
    return False


# ---------------------------------------------------------------------------
# Coverage build
# ---------------------------------------------------------------------------
def ensure_coverage_binary(build_dir, rebuild, jobs, log):
    """Configure + build the coverage `rdb` binary if needed. Return its path.

    Reuses an existing binary unless --rebuild. Streams build output to `log`.
    """
    build_dir = os.path.abspath(build_dir)
    binary = os.path.join(build_dir, "rdb")
    if os.path.isfile(binary) and not rebuild:
        log("coverage binary present: {} (reuse; --rebuild to force)"
            .format(binary))
        return binary

    log("configuring coverage build in {} ...".format(build_dir))
    cfg = subprocess.run(
        ["cmake", "-B", build_dir, "-S", REPO_ROOT,
         "-DRDB_ENABLE_COVERAGE=ON",
         "-DRDB_ENABLE_GPU=OFF",
         "-DRDB_ENABLE_MPI=OFF"],
        cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    if cfg.returncode != 0:
        sys.stdout.write(cfg.stdout.decode("utf-8", "replace"))
        raise RuntimeError("cmake configure failed (rc={})".format(cfg.returncode))

    log("building coverage `rdb` target (-j {}) ...".format(jobs))
    bld = subprocess.run(
        ["cmake", "--build", build_dir, "--target", "rdb", "-j", str(jobs)],
        cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    if bld.returncode != 0:
        sys.stdout.write(bld.stdout.decode("utf-8", "replace"))
        raise RuntimeError("cmake build failed (rc={})".format(bld.returncode))
    if not os.path.isfile(binary):
        raise RuntimeError("build succeeded but {} is missing".format(binary))
    log("coverage binary built: {}".format(binary))
    return binary


# ---------------------------------------------------------------------------
# Counter management
# ---------------------------------------------------------------------------
def zero_gcda(build_dir, log):
    """Delete every `.gcda` under build_dir so the corpus starts from zero
    counters (a fresh, reproducible measurement)."""
    n = 0
    for root, _dirs, files in os.walk(build_dir):
        for f in files:
            if f.endswith(".gcda"):
                os.remove(os.path.join(root, f))
                n += 1
    log("zeroed {} stale .gcda counter file(s)".format(n))


# ---------------------------------------------------------------------------
# Corpus run (reuses run_regression.run_case)
# ---------------------------------------------------------------------------
def run_corpus(cases, binary, scratch_root, budget_s, log):
    """Run every case serially against the coverage binary so `.gcda` counters
    accumulate. Return (results, total_wall_s, budget_blown)."""
    results = []
    start = time.time()
    blown = False
    for case in cases:
        res = run_regression.run_case(
            case, binary, scratch_root, gpu_id=None, keep_scratch=False)
        results.append(res)
        log("  {:<26} {:<5} {:>7.2f}s  {}".format(
            res["case"], "PASS" if res["passed"] else "FAIL",
            res["wall_s"], res["reason"]))
        if time.time() - start > budget_s:
            log("!! COVERAGE BUDGET EXCEEDED ({:.1f} min) after '{}'"
                .format(budget_s / 60.0, case["name"]))
            blown = True
            break
    return results, time.time() - start, blown


# ---------------------------------------------------------------------------
# gcov capture (stdlib fallback, always run)
# ---------------------------------------------------------------------------
def _parse_gcov_file(path):
    """Parse one .gcov file -> (source_abspath, executed_lines, total_exec_lines).

    Line format is `<count>:<lineno>:<src>`; count is `-` (non-executable),
    `#####`/`=====` (executable, never run) or an integer (times executed).
    """
    source = None
    executed = 0
    total = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            # Metadata rows: `<space>-:    0:Source:<path>` etc.
            if ":" not in line:
                continue
            count_field = line.split(":", 1)[0].strip()
            # tag rows have lineno 0
            parts = line.split(":", 2)
            if len(parts) >= 3 and parts[1].strip() == "0":
                if parts[2].startswith("Source:"):
                    source = parts[2][len("Source:"):].strip()
                continue
            if count_field == "-":
                continue  # non-executable line
            total += 1
            if count_field in ("#####", "====="):
                continue  # executable but never run
            # otherwise an integer (possibly with a '*' actual-block marker)
            digits = count_field.rstrip("*")
            if digits.isdigit() and int(digits) > 0:
                executed += 1
    return source, executed, total


def capture_gcov(build_dir, gcov_bin, gcov_scratch, log):
    """Run gcov over every instrumented object under build_dir and return a dict
    {repo_relative_source: {"executed": e, "total": t, "pct": p}}.

    Exact per-line counts are read from the emitted .gcov files (version-stable);
    an object with no .gcda (its TU never executed) reports 0% off the .gcno.
    """
    if os.path.isdir(gcov_scratch):
        shutil.rmtree(gcov_scratch, ignore_errors=True)
    os.makedirs(gcov_scratch, exist_ok=True)

    gcnos = sorted(glob.glob(os.path.join(build_dir, "**", "*.gcno"),
                             recursive=True))
    log("running gcov over {} instrumented object(s) ...".format(len(gcnos)))

    coverage = {}
    for gcno in gcnos:
        objdir = os.path.dirname(gcno)
        # Clear any stale .gcov from the previous iteration.
        for old in glob.glob(os.path.join(gcov_scratch, "*.gcov")):
            os.remove(old)
        subprocess.run(
            [gcov_bin, "-o", objdir, gcno],
            cwd=gcov_scratch, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, check=False,
        )
        for gcov_file in glob.glob(os.path.join(gcov_scratch, "*.gcov")):
            source, executed, total = _parse_gcov_file(gcov_file)
            if source is None or total == 0:
                continue
            src_abs = os.path.abspath(source)
            try:
                rel = os.path.relpath(src_abs, REPO_ROOT)
            except ValueError:
                rel = src_abs
            if rel.startswith(".."):
                continue  # outside the repo (system / dependency headers)
            rel = rel.replace(os.sep, "/")
            # Merge (a .gcov for the same source may appear from >1 object).
            prev = coverage.get(rel)
            if prev is None or executed > prev["executed"]:
                coverage[rel] = {
                    "executed": executed,
                    "total": total,
                    "pct": round(100.0 * executed / total, 2),
                }
    shutil.rmtree(gcov_scratch, ignore_errors=True)
    return coverage


# ---------------------------------------------------------------------------
# lcov capture (optional -- HTML report + cross-check numbers)
# ---------------------------------------------------------------------------
def lcov_available():
    return bool(shutil.which("lcov") and shutil.which("genhtml"))


def capture_lcov(build_dir, html_dir, log):
    """If lcov + genhtml are present, capture a .info and render HTML. Return the
    lcov overall (hit, total) as a cross-check, or None if unavailable/failed."""
    if not lcov_available():
        log("lcov/genhtml not on PATH -- skipping HTML report (gcov fallback "
            "is authoritative).")
        return None
    info = os.path.join(build_dir, "coverage.info")
    filtered = os.path.join(build_dir, "coverage_filtered.info")
    # lcov 2.x wants explicit ignore flags for the mismatch/unused warnings that
    # abort an otherwise-fine capture; harmless on 1.x-that-knows-them, so probe.
    ver = subprocess.run(["lcov", "--version"], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, check=False)
    vtxt = ver.stdout.decode("utf-8", "replace")
    major = 1
    for tok in vtxt.replace(".", " ").split():
        if tok.isdigit():
            major = int(tok)
            break
    cap_ignore = (["--ignore-errors", "mismatch,gcov,source,unused"]
                  if major >= 2 else ["--ignore-errors", "gcov,source"])
    rm_ignore = ["--ignore-errors", "unused"] if major >= 2 else []
    try:
        rc = subprocess.run(
            ["lcov", "--directory", build_dir, "--capture",
             "--output-file", info] + cap_ignore,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        if rc.returncode != 0:
            log("lcov capture failed (rc={}); using gcov fallback only."
                .format(rc.returncode))
            return None
        subprocess.run(
            ["lcov", "--remove", info, "/usr/*", "/opt/*", "*/build*/*",
             "*/tests/*", "--output-file", filtered] + rm_ignore,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        subprocess.run(
            ["genhtml", filtered, "--output-directory", html_dir,
             "--ignore-errors", "source"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        log("lcov HTML report: {}".format(os.path.join(html_dir, "index.html")))
    except OSError as exc:
        log("lcov invocation error: {}".format(exc))
        return None
    # Parse the filtered .info for an overall cross-check.
    hit = total = 0
    src = filtered if os.path.isfile(filtered) else info
    try:
        with open(src, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("LH:"):
                    hit += int(line[3:].strip())
                elif line.startswith("LF:"):
                    total += int(line[3:].strip())
    except OSError:
        return None
    return (hit, total)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def aggregate(cov, predicate=None):
    """Sum executed/total over cov entries whose relpath satisfies predicate."""
    ex = tot = 0
    for rel, d in cov.items():
        if predicate is None or predicate(rel):
            ex += d["executed"]
            tot += d["total"]
    pct = round(100.0 * ex / tot, 2) if tot else 0.0
    return ex, tot, pct


def build_report(cov, gap_threshold):
    """Return (ocean_rows, gaps) sorted ascending by coverage.

    ocean_rows: [(rel, pct, executed, total)] for every ocean-closure source.
    gaps:       the subset at/under gap_threshold percent (physics not exercised).
    """
    rows = []
    for rel, d in cov.items():
        if is_ocean_closure(rel):
            rows.append((rel, d["pct"], d["executed"], d["total"]))
    rows.sort(key=lambda r: (r[1], -r[3], r[0]))  # pct asc, then larger files
    gaps = [r for r in rows if r[1] <= gap_threshold]
    return rows, gaps


def print_table(rows, title):
    print("\n" + title)
    print("-" * 74)
    print("{:<52} {:>7} {:>12}".format("file", "cov%", "exec/total"))
    print("-" * 74)
    for rel, pct, ex, tot in rows:
        short = rel
        if len(short) > 52:
            short = "..." + short[-49:]
        print("{:<52} {:>6.1f}% {:>12}".format(
            short, pct, "{}/{}".format(ex, tot)))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="P2 gcov coverage mode: measure ocean-closure coverage of "
                    "the regression corpus + emit the per-closure gap list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--build-dir", default=DEFAULT_BUILD_DIR,
                   help="Coverage build dir (relative to repo root unless "
                        "absolute). Reused if it already holds a rdb binary.")
    p.add_argument("--rebuild", action="store_true",
                   help="Force a fresh coverage configure + build even if the "
                        "binary exists.")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 4,
                   help="Parallel build jobs for the coverage build.")
    p.add_argument("--gcov", default="gcov",
                   help="gcov executable (must match the gfortran that built "
                        "the objects -- load that toolchain first).")
    p.add_argument("--cases", default=None,
                   help="Comma-separated case-name filter (default: all ocean "
                        "cases in the manifest).")
    p.add_argument("--no-lcov", action="store_true",
                   help="Skip the lcov HTML report even if lcov is on PATH "
                        "(the gcov fallback still runs).")
    p.add_argument("--gap-threshold", type=float, default=10.0,
                   help="A closure file at/under this coverage%% is a GAP.")
    p.add_argument("--budget-min", type=float, default=DEFAULT_BUDGET_MIN,
                   help="Corpus-run wallclock budget in minutes (looser than "
                        "the correctness pass; instrumented runs are slow).")
    p.add_argument("--skip-run", action="store_true",
                   help="Do NOT run the corpus; capture coverage from whatever "
                        ".gcda counters already exist (debugging).")
    p.add_argument("--out", default=None,
                   help="Path for the JSON summary (default: "
                        "<build-dir>/coverage_summary.json).")
    p.add_argument("--scratch-root", default=None,
                   help="Root for per-case scratch dirs (default: "
                        "<repo>/tmp_local_artifacts/coverage).")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    def log(msg):
        print("[coverage] {}".format(msg), flush=True)

    t_start = time.time()

    build_dir = args.build_dir
    if not os.path.isabs(build_dir):
        build_dir = os.path.join(REPO_ROOT, build_dir)

    # --- 1. coverage binary ------------------------------------------------
    try:
        binary = ensure_coverage_binary(build_dir, args.rebuild, args.jobs, log)
    except RuntimeError as exc:
        log("ERROR: {}".format(exc))
        return 2

    scratch_root = args.scratch_root or os.path.join(
        REPO_ROOT, "tmp_local_artifacts", "coverage")
    os.makedirs(scratch_root, exist_ok=True)

    cases = list(manifest.CASES)
    if args.cases:
        wanted = {c.strip() for c in args.cases.split(",") if c.strip()}
        cases = [c for c in cases if c["name"] in wanted]
    cases = [c for c in cases if "cpu" in c.get("backends", {"cpu"})]

    # --- 2. run corpus (accumulate .gcda) ----------------------------------
    corpus_results = []
    corpus_wall = 0.0
    blown = False
    if args.skip_run:
        log("--skip-run: capturing from existing .gcda counters.")
    else:
        zero_gcda(build_dir, log)
        log("running {} ocean case(s) SERIALLY against the coverage binary "
            "...".format(len(cases)))
        corpus_results, corpus_wall, blown = run_corpus(
            cases, binary, scratch_root, args.budget_min * 60.0, log)
        n_pass = sum(1 for r in corpus_results if r["passed"])
        log("corpus done: {}/{} clean, {:.1f}s ({:.1f} min) wall".format(
            n_pass, len(corpus_results), corpus_wall, corpus_wall / 60.0))

    # --- 3a. gcov fallback (always) ----------------------------------------
    gcov_scratch = os.path.join(scratch_root, "_gcov")
    cov = capture_gcov(build_dir, args.gcov, gcov_scratch, log)
    if not cov:
        log("ERROR: gcov produced no coverage data -- was the corpus run?")
        return 2

    # --- 3b. lcov HTML + cross-check (optional) ----------------------------
    lcov_xcheck = None
    used_lcov = False
    if not args.no_lcov:
        html_dir = os.path.join(build_dir, "coverage_report")
        lcov_xcheck = capture_lcov(build_dir, html_dir, log)
        used_lcov = lcov_xcheck is not None

    # --- 4. report ---------------------------------------------------------
    ex_all, tot_all, pct_all = aggregate(cov)
    ex_oc, tot_oc, pct_oc = aggregate(cov, is_ocean_closure)
    ocean_rows, gaps = build_report(cov, args.gap_threshold)

    print("\n" + "=" * 74)
    print("OCEAN-CLOSURE COVERAGE  (regression corpus, {} cases)".format(
        len(corpus_results) if not args.skip_run else "existing .gcda"))
    print("=" * 74)
    print("overall (all instrumented src) : {:6.2f}%  ({}/{} lines)".format(
        pct_all, ex_all, tot_all))
    print("overall (ocean closures only)  : {:6.2f}%  ({}/{} lines)".format(
        pct_oc, ex_oc, tot_oc))
    if lcov_xcheck:
        lh, lf = lcov_xcheck
        lpct = round(100.0 * lh / lf, 2) if lf else 0.0
        print("lcov cross-check (filtered)    : {:6.2f}%  ({}/{} lines)".format(
            lpct, lh, lf))

    print_table(ocean_rows, "PER-CLOSURE COVERAGE  (ascending -- worst first)")

    print("\nGAPS -- ocean-closure files at/under {:.0f}% (P4 case targets)"
          .format(args.gap_threshold))
    print("-" * 74)
    if gaps:
        for rel, pct, ex, tot in gaps:
            print("  {:6.1f}%  {}  ({} exec lines)".format(pct, rel, tot))
    else:
        print("  (none -- every ocean-closure file is above the threshold)")
    print("=" * 74)

    total_wall = time.time() - t_start
    log("total coverage-run wallclock: {:.1f}s ({:.1f} min)".format(
        total_wall, total_wall / 60.0))

    # --- machine-readable summary -----------------------------------------
    out_path = args.out or os.path.join(build_dir, "coverage_summary.json")
    summary = {
        "meta": {
            "phase": "P2-gcov-coverage",
            "binary": binary,
            "build_dir": build_dir,
            "gcov": args.gcov,
            "lcov_used": used_lcov,
            "gap_threshold_pct": args.gap_threshold,
            "corpus_wall_s": round(corpus_wall, 2),
            "total_wall_s": round(total_wall, 2),
            "budget_blown": blown,
            "n_cases": len(corpus_results),
            "n_clean": sum(1 for r in corpus_results if r["passed"]),
        },
        "overall_all": {"executed": ex_all, "total": tot_all, "pct": pct_all},
        "overall_ocean_closures": {
            "executed": ex_oc, "total": tot_oc, "pct": pct_oc},
        "lcov_crosscheck": (
            {"hit": lcov_xcheck[0], "total": lcov_xcheck[1]}
            if lcov_xcheck else None),
        "per_closure": [
            {"file": rel, "pct": pct, "executed": ex, "total": tot}
            for (rel, pct, ex, tot) in ocean_rows
        ],
        "gaps": [
            {"file": rel, "pct": pct, "total": tot}
            for (rel, pct, ex, tot) in gaps
        ],
        "corpus_results": [
            {"case": r["case"], "passed": r["passed"],
             "wall_s": r["wall_s"], "reason": r["reason"]}
            for r in corpus_results
        ],
    }
    with open(out_path, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    log("JSON summary: {}".format(os.path.abspath(out_path)))

    # Exit nonzero only on a hard failure (budget blown or no coverage data),
    # NOT on low coverage -- low coverage is the report's whole point.
    return 1 if blown else 0


if __name__ == "__main__":
    sys.exit(main())
