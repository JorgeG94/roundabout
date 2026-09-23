#!/usr/bin/env python3
"""Single entry point for the ocean physics-coverage regression suite.

Invokes the whole flow and returns a CI-ready exit code (0 = all good):

  1. Golden-summary regression compare (run-clean + NaN gate + drift-vs-golden)
     on the chosen backend(s) -- CPU (gfortran) and/or GPU (nvfortran).  This
     subsumes the P0 run-clean gate: `compare.py` runs every case and checks
     both liveness (exit 0, no NaN/Inf) and drift against the committed goldens.
  2. gcov physics-coverage measurement (gfortran --coverage build) -- on by
     default; skip with --no-coverage.  Coverage is a measurement, so it only
     fails the suite on an infrastructure error (budget blown / no data), never
     on "coverage is low".

This is deliberately a STANDALONE script, not a make/CTest target -- CI can call
it later.  Each stage runs in its own subshell, so `--backend both` works from a
single invocation without ever stacking a CPU and a GPU toolchain in one shell
(two NetCDF builds on one link line fail in confusing ways).

Toolchain: load gfortran / NVHPC + NetCDF however your site does it (module,
Spack -- see environments/ --, conda) before running; each stage inherits the
environment it is launched in.  On a configured dev box that keeps untracked
site-specific env scripts at the repo root, export RDB_ON_DEV=1 (and optionally
--cpu-env / --gpu-env) to have each stage's subshell source them; a silent
no-op otherwise.

The application binaries are assumed already built: `build_gcc/rdb` (CPU),
`build_gpu/rdb` (GPU).  The coverage stage builds its own `build_cov`.

Usage:
  python3 tests/regression/run_all.py                     # CPU compare + coverage
  python3 tests/regression/run_all.py --backend both      # CPU + GPU + coverage
  python3 tests/regression/run_all.py --backend gpu --gpus 0,1,2,3 --no-coverage
  python3 tests/regression/run_all.py --update-golden      # regenerate goldens (CPU)
"""
import argparse
import os
import pathlib
import shlex
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent  # tests/regression


# Root-only markers: `tests/` and other subdirectories carry their own
# CMakeLists.txt, so that file is NOT a usable marker.
ROOT_MARKERS = ("CMakePresets.json", "fpm.toml", ".git")


def find_root(start):
    """Walk up to the repo/worktree root (the dir holding CMakePresets.json)."""
    p = start
    for _ in range(8):
        if any((p / m).exists() for m in ROOT_MARKERS):
            return p
        p = p.parent
    return start.parents[1]


ROOT = find_root(HERE)


def dev_env_script(name):
    """Return a local dev toolchain script to source for a stage, or None.

    Roundabout ships no machine-specific env scripts, so the normal path is
    None: the stage runs in the environment it inherits (toolchain already
    loaded via module / Spack / conda).  On a configured dev box that keeps
    untracked site-specific scripts at the repo root, export RDB_ON_DEV=1 and
    the named script is sourced instead.  Missing script => still None, so
    this never breaks a clean checkout.
    """
    if not name or os.environ.get("RDB_ON_DEV", "0") != "1":
        return None
    path = pathlib.Path(name)
    if not path.is_absolute():
        path = ROOT / path
    return path if path.is_file() else None


def run_stage(title, env_script, script, args):
    """Run one suite stage in a fresh subshell.

    `env_script` is an optional path sourced first (see `dev_env_script`).

    Returns the stage's exit code (0 = pass).
    """
    inner = (
        "python3 " + shlex.quote("tests/regression/" + script) + " "
        + " ".join(shlex.quote(a) for a in args)
    )
    cmd = "cd " + shlex.quote(str(ROOT)) + " && " + inner
    if env_script is not None:
        cmd = "source " + shlex.quote(str(env_script)) + " && " + cmd
    print("\n" + "=" * 72)
    print(">>> " + title)
    print("    env : " + (str(env_script) if env_script is not None
                          else "(inherited)"))
    print("    run : " + script + " " + " ".join(args))
    print("=" * 72, flush=True)
    rc = subprocess.run(["bash", "-c", cmd]).returncode
    print("<<< " + title + ": " + ("PASS" if rc == 0 else "FAIL (rc=%d)" % rc), flush=True)
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Run the full ocean physics-coverage regression suite "
        "(CI-ready exit code; standalone, not a make/CTest target).")
    ap.add_argument("--backend", choices=["cpu", "gpu", "both"], default="cpu",
                    help="which backend(s) to run the regression compare on.")
    ap.add_argument("--cpu-build-dir", default="build_gcc",
                    help="gfortran app build dir (holds rdb).")
    ap.add_argument("--gpu-build-dir", default="build_gpu",
                    help="nvfortran app build dir (holds rdb).")
    ap.add_argument("--coverage-build-dir", default="build_cov",
                    help="gcov --coverage build dir (coverage.py builds it).")
    ap.add_argument("--gpus", default=None,
                    help="comma GPU ids for the gpu farm, e.g. 0,1,2,3.")
    ap.add_argument("--no-coverage", action="store_true",
                    help="skip the gcov coverage stage.")
    ap.add_argument("--update-golden", action="store_true",
                    help="regenerate the goldens from the CPU build instead of comparing "
                    "(implies --backend cpu, skips coverage).")
    ap.add_argument("--cpu-env", default="gcc_env.sh",
                    help="local dev toolchain script for the CPU/coverage "
                         "stages; sourced ONLY when RDB_ON_DEV=1 and it "
                         "exists. Otherwise the inherited environment is used.")
    ap.add_argument("--gpu-env", default="nvhpc_env.sh",
                    help="local dev toolchain script for the GPU stage; "
                         "sourced ONLY when RDB_ON_DEV=1 and it exists. "
                         "Otherwise the inherited environment is used.")
    ap.add_argument("--rtol", default=None, help="override compare relative tolerance.")
    args = ap.parse_args(argv)

    cpu_env = dev_env_script(args.cpu_env)
    gpu_env = dev_env_script(args.gpu_env)

    stages = {}

    if args.update_golden:
        cargs = ["--update-golden", "--backend", "cpu", "--build-dir", args.cpu_build_dir]
        stages["golden-update-cpu"] = run_stage(
            "Golden update (cpu)", cpu_env, "compare.py", cargs)
    else:
        backends = ["cpu", "gpu"] if args.backend == "both" else [args.backend]
        for b in backends:
            env = gpu_env if b == "gpu" else cpu_env
            bd = args.gpu_build_dir if b == "gpu" else args.cpu_build_dir
            cargs = ["--compare", "--backend", b, "--build-dir", bd]
            if b == "gpu" and args.gpus:
                cargs += ["--gpus", args.gpus]
            if args.rtol:
                cargs += ["--rtol", args.rtol]
            stages["compare-" + b] = run_stage(
                "Regression compare (" + b + ")", env, "compare.py", cargs)

        if not args.no_coverage:
            stages["coverage"] = run_stage(
                "Coverage (gcov)", cpu_env, "coverage.py",
                ["--build-dir", args.coverage_build_dir])

    print("\n" + "=" * 72)
    print("=== REGRESSION SUITE SUMMARY ===")
    ok = True
    for name, rc in stages.items():
        print("  %-24s %s" % (name, "PASS" if rc == 0 else "FAIL (rc=%d)" % rc))
        ok = ok and rc == 0
    print("  %-24s %s" % ("OVERALL", "PASS" if ok else "FAIL"))
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
