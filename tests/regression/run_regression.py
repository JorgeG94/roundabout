#!/usr/bin/env python3
"""P0 ocean regression runner -- run-clean / NaN gate (no golden compare yet).

Runs each curated ocean input case (TOML) (see manifest.py) for a short number of
outer timesteps and asserts it RUNS CLEAN: exit code 0 and no non-finite /
crash markers in the captured output. This is the P0 bar. Golden-summary
comparison is P1; gcov coverage is P2.

Stdlib only (the repo forbids `pip install`); no third-party imports.

Per-case execution
------------------
1. Copy the committed TOML input into a per-case scratch dir.
2. Patch `[time]` so the run is `n_steps` outer steps:
     t_end     = n_steps * dt_fixed   (seconds)
     time_unit = "second"             (normalized; so t_end is read in seconds)
   The committed input is never mutated.
3. Redirect `[output] output_dir` to the scratch dir AND run the process with
   cwd=<scratch>, so no run ever writes into the repo tree.
4. Run `rdb <temp.toml>` with the per-case timeout; capture rc + stdout+stderr.
5. PASS iff rc == 0 AND no bad markers (nan / inf / error stop / not finite /
   abort / panic / segmentation) appear in the output (case-insensitive).

Backends
--------
--backend cpu : run cases serially (the gfortran CPU build).
--backend gpu : farm cases across visible GPUs, ONE case per GPU at a time (a
                GPU is never shared), by exporting CUDA_VISIBLE_DEVICES=<id> per
                case. Up to (num GPUs) cases run concurrently via a work queue.

Budget
------
The runner tracks cumulative wallclock and FAILS LOUD if a backend exceeds the
budget (default 30 min) -- a silent 45-min suite is a bug. Per-case timeouts
catch individual hangs.

Reporting
---------
A fixed-width table to stdout plus a machine-readable JSON results file
(--out). Exit code is nonzero if any case failed (or the budget was blown).

Run `python3 run_regression.py --help` for the full flag list.
"""

import argparse
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time

# ---------------------------------------------------------------------------
# Paths / manifest import (stdlib only; add this file's dir to sys.path so the
# runner works regardless of the caller's cwd).
# ---------------------------------------------------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
# repo root is two levels up: <repo>/tests/regression/run_regression.py
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
sys.path.insert(0, THIS_DIR)
import manifest  # noqa: E402  (local, after sys.path tweak)
import toml_util  # noqa: E402

# Output markers that signal a broken / non-finite / crashed run, as
# case-insensitive regexes scanned against combined stdout+stderr.
#
# nan / inf / infinity are WORD-BOUNDARY matched -- a bare substring "inf"
# false-matches benign banner text ("informed physical vibes ...") and "nan"
# false-matches "nanometre" etc. `\b(inf|infinity)\b` still catches the
# Fortran IEEE prints "Inf", "+Inf", "-Inf", "Infinity" while ignoring
# "informed"; `\bnan\b` catches "NaN" without matching longer words.
BAD_MARKERS = tuple(re.compile(pat, re.IGNORECASE) for pat in (
    r"\bnan\b",
    r"\binf(inity)?\b",
    r"error stop",
    r"not finite",
    r"\babort(ed|ing)?\b",
    r"\bpanic\b",
    r"segmentation",
    r"floating point exception",
    r"backtrace",
))

DEFAULT_BUDGET_MIN = 30.0


# ---------------------------------------------------------------------------
# Input patching
# ---------------------------------------------------------------------------
def patch_input(src_text, n_steps, scratch_dir, force_emit=True):
    """Return a patched copy of the TOML input `src_text` for a short,
    output-isolated run.

    - `[time]`: dt_fixed is read, t_end is set to n_steps*dt_fixed seconds
      and time_unit is forced to "second".
    - `[output]`: output_dir is redirected to `scratch_dir`.
    - When `force_emit` (default True): the console-status cadence
      (`[logging] status_interval`) and the ocean diag-manager cadence
      (`[ocean.diag] dt_out`) are forced to a value strictly below one
      outer step (`0.5*dt_fixed` seconds, once time_unit is normalized to
      seconds) so BOTH fire every step -- guaranteeing a final-state
      `[stats]` line and a final-state `[diag]` summary for every field.
      Without this a short run only emits at t=0 (the committed cadences are
      hours/days), leaving the P1 golden-summary compare (compare.py) with no
      usable final state. Harmless to the P0 run-clean gate (it only adds log
      lines / more NaN sampling). If the source omits `status_interval` /
      `dt_out` inside a table that IS present, the key is added to that
      table; tables that are entirely absent are left untouched.

    Raises ValueError if dt_fixed cannot be found (we cannot bound the run).
    The edits go through `toml_util.patch` (line-level, stdlib only).
    """
    dt_fixed = toml_util.get_float(src_text, "time_nml", "dt_fixed")
    if dt_fixed is None:
        raise ValueError("could not find dt_fixed in [time]")

    over = {"time_nml": {"t_end": n_steps * dt_fixed, "time_unit": "second"}}
    if toml_util.has_group(src_text, "output_nml"):
        over["output_nml"] = {"output_dir": scratch_dir}
    if force_emit:
        # A cadence strictly below one outer step => the diag/status probes
        # fire on every step, so the last step always emits a final summary.
        emit_cadence = 0.5 * dt_fixed
        if toml_util.has_group(src_text, "logging_nml"):
            over["logging_nml"] = {"status_interval": emit_cadence}
        if toml_util.has_group(src_text, "ocean_diag_nml"):
            over["ocean_diag_nml"] = {"dt_out": emit_cadence}
    return toml_util.patch(src_text, over)


# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------
def locate_binary(build_dir, explicit):
    """Return an absolute path to the rdb executable.

    `explicit` (from --binary) wins. Otherwise look at <build_dir>/rdb, then
    fall back to a recursive search under build_dir.
    """
    if explicit:
        path = os.path.abspath(explicit)
        if not os.path.isfile(path):
            raise FileNotFoundError("--binary not found: {}".format(path))
        return path

    build_dir = os.path.abspath(build_dir)
    direct = os.path.join(build_dir, "rdb")
    if os.path.isfile(direct):
        return direct
    for root, _dirs, files in os.walk(build_dir):
        if "rdb" in files:
            cand = os.path.join(root, "rdb")
            if os.access(cand, os.X_OK):
                return cand
    raise FileNotFoundError(
        "could not find 'rdb' under {} -- build the CPU app first "
        "(see tests/regression/README.md)".format(build_dir)
    )


# ---------------------------------------------------------------------------
# GPU discovery
# ---------------------------------------------------------------------------
def detect_gpus(explicit):
    """Return a list of GPU id strings.

    `explicit` (from --gpus, e.g. "0,1,2,3") wins. Otherwise query
    `nvidia-smi -L`; tolerate its absence by falling back to ["0"].
    """
    if explicit:
        return [tok.strip() for tok in explicit.split(",") if tok.strip()]
    try:
        out = subprocess.run(
            ["nvidia-smi", "-L"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            timeout=30,
            check=False,
        ).stdout.decode("utf-8", "replace")
    except (OSError, subprocess.SubprocessError):
        return ["0"]
    ids = []
    for line in out.splitlines():
        line = line.strip()
        if line.lower().startswith("gpu ") and ":" in line:
            # "GPU 0: NVIDIA ..." -> "0"
            ids.append(line.split(":", 1)[0].split()[1])
    return ids or ["0"]


# ---------------------------------------------------------------------------
# Single-case execution
# ---------------------------------------------------------------------------
def scan_output(text):
    """Return the first bad marker matched in `text`, or None.

    Markers are word-boundary-aware regexes (see BAD_MARKERS) so benign banner
    text like "informed physical vibes" does not false-trip the "inf" gate.
    """
    for marker in BAD_MARKERS:
        m = marker.search(text)
        if m:
            return m.group(0)
    return None


def run_case(case, binary, scratch_root, gpu_id=None, keep_scratch=False,
             keep_output=False):
    """Run one case end-to-end. Return a result dict.

    gpu_id: if not None, CUDA_VISIBLE_DEVICES is set to it for the child.
    keep_output: if True, the FULL combined stdout/stderr is attached to the
        result under the "output" key (compare.py parses the [diag]/[stats]
        summary lines from it). Left off by default so the P0 results JSON
        stays small -- only the 25-line "stdout_tail" is kept there.
    """
    name = case["name"]
    scratch = os.path.join(scratch_root, name)
    # Fresh scratch dir per run.
    if os.path.isdir(scratch):
        shutil.rmtree(scratch, ignore_errors=True)
    os.makedirs(scratch, exist_ok=True)

    result = {
        "case": name,
        "backend": "gpu" if gpu_id is not None else "cpu",
        "gpu_id": gpu_id,
        "passed": False,
        "wall_s": 0.0,
        "returncode": None,
        "reason": "",
        "nml": case["nml"],
        "n_steps": case["n_steps"],
        "stdout_tail": "",
    }

    src_path = os.path.join(REPO_ROOT, case["nml"])
    try:
        with open(src_path, "r") as fh:
            src_text = fh.read()
    except OSError as exc:
        result["reason"] = "cannot read input: {}".format(exc)
        return result

    try:
        patched = patch_input(src_text, case["n_steps"], scratch)
    except ValueError as exc:
        result["reason"] = "input patch failed: {}".format(exc)
        return result

    temp_input = os.path.join(scratch, name + ".toml")
    with open(temp_input, "w") as fh:
        fh.write(patched)

    # Optional per-case setup command, run INSIDE the scratch dir before the
    # model. This exists for cases whose input names an input file the repo
    # deliberately does not commit -- e.g. the PR-15 file-forcing case, whose
    # wind NetCDF is generated data. The model runs with cwd=scratch, so a
    # bare filename in the input resolves to whatever setup dropped there.
    # argv is taken verbatim except for "{python}", which becomes the
    # interpreter running this script, and "{repo}", the repo root.
    setup = case.get("setup")
    if setup:
        argv = [a.replace("{python}", sys.executable).replace("{repo}", REPO_ROOT)
                for a in setup]
        try:
            sp = subprocess.run(argv, cwd=scratch, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        except (OSError, subprocess.TimeoutExpired) as exc:
            result["reason"] = "setup failed: {}".format(exc)
            return result
        if sp.returncode != 0:
            tail = sp.stdout.decode("utf-8", "replace").strip().splitlines()[-5:]
            result["reason"] = "setup exited {}: {}".format(
                sp.returncode, " | ".join(tail))
            return result

    env = dict(os.environ)
    # NetCDF/HDF5 in this tree is not thread-safe under stdpar-multicore; keep
    # the child single-threaded (matches the ctest ENVIRONMENT convention).
    env["OMP_NUM_THREADS"] = "1"
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    t0 = time.time()
    timed_out = False
    try:
        proc = subprocess.run(
            [binary, temp_input],
            cwd=scratch,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=case["timeout_s"],
            check=False,
        )
        rc = proc.returncode
        output = proc.stdout.decode("utf-8", "replace")
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        rc = None
        output = (exc.stdout or b"").decode("utf-8", "replace")
    result["wall_s"] = round(time.time() - t0, 2)
    result["returncode"] = rc

    # Keep a tail of the output for the report / debugging.
    tail_lines = output.strip().splitlines()[-25:]
    result["stdout_tail"] = "\n".join(tail_lines)
    if keep_output:
        result["output"] = output

    # ---- Gate ----
    if timed_out:
        result["reason"] = "TIMEOUT after {}s".format(case["timeout_s"])
    elif rc != 0:
        marker = scan_output(output)
        result["reason"] = "exit rc={}{}".format(
            rc, " ({})".format(marker) if marker else ""
        )
    else:
        marker = scan_output(output)
        if marker:
            result["reason"] = "bad marker: '{}'".format(marker)
        else:
            result["passed"] = True
            result["reason"] = "ok"

    # Clean scratch on success unless asked to keep it; keep on failure.
    if result["passed"] and not keep_scratch:
        shutil.rmtree(scratch, ignore_errors=True)

    return result


# ---------------------------------------------------------------------------
# Backend drivers
# ---------------------------------------------------------------------------
def run_cpu(cases, binary, scratch_root, budget_s, keep_scratch):
    """Run cases serially on CPU. Returns (results, budget_blown)."""
    results = []
    start = time.time()
    for case in cases:
        res = run_case(case, binary, scratch_root, gpu_id=None,
                       keep_scratch=keep_scratch)
        results.append(res)
        _print_row(res)
        if time.time() - start > budget_s:
            print("\n!! BUDGET EXCEEDED ({:.1f} min) after case '{}' -- "
                  "stopping.".format(budget_s / 60.0, case["name"]))
            return results, True
    return results, False


def run_gpu(cases, binary, scratch_root, gpus, budget_s, keep_scratch):
    """Farm cases across GPUs, one case per GPU at a time. Returns (results, blown).

    A work queue hands cases to a pool of worker threads -- exactly one worker
    per GPU id, so a given GPU never runs two cases at once. Each worker pins
    its case with CUDA_VISIBLE_DEVICES=<its gpu id>.
    """
    work = queue.Queue()
    for case in cases:
        work.put(case)

    results = []
    results_lock = threading.Lock()
    start = time.time()
    blown = threading.Event()

    def worker(gpu_id):
        while True:
            if blown.is_set():
                return
            try:
                case = work.get_nowait()
            except queue.Empty:
                return
            try:
                res = run_case(case, binary, scratch_root, gpu_id=gpu_id,
                               keep_scratch=keep_scratch)
                with results_lock:
                    results.append(res)
                    _print_row(res)
                    if time.time() - start > budget_s:
                        print("\n!! BUDGET EXCEEDED ({:.1f} min) on GPU {} -- "
                              "stopping.".format(budget_s / 60.0, gpu_id))
                        blown.set()
            finally:
                work.task_done()

    threads = [threading.Thread(target=worker, args=(g,), daemon=True)
               for g in gpus]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    return results, blown.is_set()


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
_HEADER_PRINTED = False


def _print_header():
    global _HEADER_PRINTED
    if _HEADER_PRINTED:
        return
    print("{:<26} {:<7} {:<6} {:>8}  {}".format(
        "case", "backend", "result", "wall_s", "reason"))
    print("-" * 78)
    _HEADER_PRINTED = True


def _print_row(res):
    _print_header()
    print("{:<26} {:<7} {:<6} {:>8}  {}".format(
        res["case"],
        res["backend"],
        "PASS" if res["passed"] else "FAIL",
        "{:.2f}".format(res["wall_s"]),
        res["reason"],
    ))


def write_json(path, results, meta):
    payload = {"meta": meta, "results": results}
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="P0 ocean regression runner (run-clean / NaN gate).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--backend", choices=["cpu", "gpu"], default="cpu",
                   help="cpu: serial; gpu: one case per GPU, farmed across "
                        "visible devices.")
    p.add_argument("--build-dir", default="build_gcc",
                   help="Build directory containing the 'rdb' binary "
                        "(relative to repo root unless absolute).")
    p.add_argument("--binary", default=None,
                   help="Explicit path to the rdb binary (overrides "
                        "--build-dir discovery).")
    p.add_argument("--out", default="results.json",
                   help="Path to write the machine-readable JSON results.")
    p.add_argument("--gpus", default=None,
                   help="Comma-separated GPU ids for --backend gpu (e.g. "
                        "'0,1,2,3'). Default: auto-detect via nvidia-smi -L.")
    p.add_argument("--cases", default=None,
                   help="Comma-separated case-name filter (default: all in the "
                        "manifest).")
    p.add_argument("--budget-min", type=float, default=DEFAULT_BUDGET_MIN,
                   help="Per-backend wallclock budget in minutes; exceeding it "
                        "is a loud failure.")
    p.add_argument("--keep-scratch", action="store_true",
                   help="Keep per-case scratch dirs even on success (debugging).")
    p.add_argument("--scratch-root", default=None,
                   help="Root for per-case scratch dirs (default: "
                        "<repo>/tmp_local_artifacts/regression).")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # Resolve build dir relative to repo root when not absolute.
    build_dir = args.build_dir
    if not os.path.isabs(build_dir):
        build_dir = os.path.join(REPO_ROOT, build_dir)

    try:
        binary = locate_binary(build_dir, args.binary)
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
    # Restrict to cases valid on the chosen backend.
    cases = [c for c in cases if args.backend in c.get("backends", {args.backend})]
    if not cases:
        print("ERROR: no cases to run for backend '{}'".format(args.backend),
              file=sys.stderr)
        return 2

    budget_s = args.budget_min * 60.0

    print("rdb binary : {}".format(binary))
    print("backend       : {}".format(args.backend))
    print("cases         : {}".format(len(cases)))
    print("budget        : {:.0f} min".format(args.budget_min))

    t0 = time.time()
    if args.backend == "cpu":
        print("mode          : serial CPU\n")
        results, blown = run_cpu(cases, binary, scratch_root, budget_s,
                                 args.keep_scratch)
    else:
        gpus = detect_gpus(args.gpus)
        print("gpus          : {}\n".format(",".join(gpus)))
        results, blown = run_gpu(cases, binary, scratch_root, gpus, budget_s,
                                 args.keep_scratch)
    total_wall = time.time() - t0

    n_pass = sum(1 for r in results if r["passed"])
    n_total = len(results)
    n_ran = n_total  # results only holds cases that actually ran

    print("\n" + "=" * 78)
    print("SUMMARY: {}/{} passed  |  total wall {:.1f}s ({:.1f} min)  |  "
          "budget {:.0f} min".format(
              n_pass, len(cases), total_wall, total_wall / 60.0,
              args.budget_min))
    if n_ran < len(cases):
        print("         {} case(s) did not run (budget stop).".format(
            len(cases) - n_ran))
    if blown:
        print("         !! BUDGET BLOWN -- failing loud.")

    meta = {
        "backend": args.backend,
        "binary": binary,
        "n_cases": len(cases),
        "n_ran": n_ran,
        "n_pass": n_pass,
        "total_wall_s": round(total_wall, 2),
        "budget_min": args.budget_min,
        "budget_blown": blown,
        "phase": "P0-run-clean-nan-gate",
    }
    write_json(args.out, results, meta)
    print("         JSON: {}".format(os.path.abspath(args.out)))

    failed = (n_pass != len(cases)) or blown
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
