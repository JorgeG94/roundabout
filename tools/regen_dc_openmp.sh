#!/usr/bin/env bash
# regen_dc_openmp.sh — regenerate the dc-openmp variant tree from main.
#
# ONE variant ships: `dc-openmp-target`.
#   !$acc data/reduction directives → !$omp target equivalents, but the
#   `do concurrent` loops are KEPT verbatim. The compiler maps them to the
#   device itself (NVHPC -stdpar=gpu -mp=gpu, LLVM flang
#   -fdo-concurrent-to-openmp=device, ifx -fopenmp-target-do-concurrent) or to
#   the host (gfortran -fopenmp, plain `do concurrent`). Produces the tree
#   that lives on branch  auto/dc-openmp.
#
# There is NO full-OpenMP variant. `do concurrent` is never rewritten to an
# `!$omp` worksharing construct — the `dc_to_omp.py` rewriter and the
# `auto/openmp` / `auto/openmp-cpu` branches it fed were retired; see
# git history if you need the old tool.
#
# Runs in the repo root and modifies the working tree IN PLACE:
#   1. Lints `do concurrent` loops (dc_audit.py --strict).
#   2. Translates !$acc directives to !$omp via acc_to_omp.py --write
#      (--target gpu: offload directives, not the host `parallel do` form).
#   3. Adds explicit map() clauses for OPTIONAL array dummies read inside
#      the emitted target regions (omp_optional_map.py).
#   4. Applies every overlay patch under patches/openmp/*/*.patch.
#   5. Flips the CMake default backend to openmp (cmake/options.cmake).
#
# Use a clean checkout of main to call this script; do NOT call from a branch
# that already has uncommitted source changes.
#
# Usage:
#   bash tools/regen_dc_openmp.sh                 # dc-openmp-target (default)
#   bash tools/regen_dc_openmp.sh dc-openmp-target
#
# An optional second argument sets the worker-process count for the per-file
# passes (0 = one per available core; default 1 = serial):
#   bash tools/regen_dc_openmp.sh dc-openmp-target 0
#   RDB_REGEN_JOBS=32 bash tools/regen_dc_openmp.sh
#
# Exits non-zero if any stage fails. Intended to be called from CI and
# (occasionally) by developers wanting to test the variant locally.

set -euo pipefail

variant="${1:-dc-openmp-target}"
if [[ "$variant" != "dc-openmp-target" ]]; then
  echo "error: unknown variant '$variant'" >&2
  echo "       only 'dc-openmp-target' ships -- do concurrent is never" >&2
  echo "       rewritten, so there is no 'openmp' / 'openmp-cpu' variant" >&2
  echo "       any more." >&2
  exit 2
fi
# Worker processes for the per-file passes. Second positional arg, or
# RDB_REGEN_JOBS. 0 = one per available core (respects SLURM binding);
# 1 (default) = serial, matching the historical behaviour.
jobs="${2:-${RDB_REGEN_JOBS:-1}}"
omp_target=gpu          # acc_to_omp --target (gpu offload vs cpu host); the
                         # dc-openmp-target variant always emits the offload
                         # form -- do concurrent, not the translator, decides
                         # what runs where.

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

echo "==> variant: $variant  (translator --target $omp_target, -j $jobs)"

echo "==> [1/5] dc_audit --strict"
python tools/dc_audit.py --strict -j "$jobs" src/ app/ benchmarks/ tests/

echo "==> [2/5] acc_to_omp --write"
# acc_to_omp.py walks the directories given on the command line. Pass src/,
# app/, benchmarks/ and tests/ so directives in app/main.F90, the bench
# drivers (bench_ocean's data + update self) and the unit tests are all
# translated — otherwise those keep !$acc, which is inert on an OpenMP build.
python tools/acc_to_omp.py --target "$omp_target" -j "$jobs" --write src/ app/ benchmarks/ tests/

echo "==> [3/5] omp_optional_map --write (guard absent OPTIONAL dummies)"
# LLVM Flang maps an OPTIONAL explicit-shape array dummy read inside a target
# region unconditionally: bounds from the declaration, base address NULL when
# the argument is absent, so libomptarget tries hsa_amd_memory_lock(0x0, size)
# and the run aborts. Naming it in an explicit map() clause takes the
# presence-guarded path.
python tools/omp_optional_map.py -j "$jobs" --write src/ app/ benchmarks/ tests/

echo "==> [4/5] apply overlay patches"
shopt -s nullglob
patches=(patches/openmp/*/*.patch)
shopt -u nullglob
if [[ ${#patches[@]} -eq 0 ]]; then
  echo "    (no patches found under patches/openmp/)"
else
  for p in "${patches[@]}"; do
    echo "    applying $p"
    git apply --whitespace=nowarn "$p"
  done
fi

echo "==> [5/5] flip CMake default backend openacc -> openmp"
# In-place sed; preserves indentation. cmake/options.cmake declares the cache
# variable with its default on its own quoted line, so the substitution is
# precise. macOS sed differs from GNU sed in -i semantics; pipe-through-temp
# is portable.
if grep -q '^    "openacc"$' cmake/options.cmake; then
  tmp=$(mktemp)
  awk '
    /set\(RDB_PARALLEL_BACKEND$/ { in_block = 1 }
    in_block && /^    "openacc"$/ { sub("openacc", "openmp"); in_block = 0 }
    { print }
  ' cmake/options.cmake > "$tmp"
  mv "$tmp" cmake/options.cmake
  echo "    done"
else
  echo "    SKIP — cmake/options.cmake default backend already not 'openacc'"
fi

echo "==> regen complete ($variant)"
echo ""
echo "Next steps (manual, for local validation):"
echo "  cmake -B build_omp -S . -DRDB_PARALLEL_BACKEND=openmp -DRDB_ENABLE_GPU=OFF \\"
echo "        -DRDB_ENABLE_THREADS=ON"
echo "  cmake --build build_omp -j"
echo "  cd build_omp && ctest -R rdb --output-on-failure"
