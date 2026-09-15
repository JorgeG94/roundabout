#!/usr/bin/env bash
# regen_dc_openmp.sh — regenerate an OpenMP-target variant tree from main.
#
# Two variants, selected by the first argument:
#
#   dc-openmp-target  (default)
#       !$acc data/reduction directives → !$omp target equivalents, but the
#       `do concurrent` loops are KEPT verbatim. The compiler maps them to the
#       device itself (NVHPC -stdpar=gpu, LLVM flang). Produces the tree that
#       lives on  auto/dc-openmp.
#       NOTE: this variant does NOT currently compile on ROCm/amdflang — a
#       crash in flang's DoConcurrentConversion pass, not a translator bug.
#       See docs/OPENMP_VARIANT_STATUS.md.
#
#   openmp
#       Everything dc-openmp-target does, PLUS dc_to_omp.py rewrites every
#       `do concurrent` as `!$omp target teams distribute parallel do`
#       (collapsing the header, mapping local()→private(), and stripping
#       `pure` from any procedure that ends up holding a target region — the
#       whole-program cascade). No `do concurrent` survives. Produces the tree
#       that lives on  auto/openmp.
#
#   openmp-cpu
#       Like openmp, but both translators run with --target cpu: compute loops
#       become host `!$omp parallel do` (NOT `target teams distribute`, which is
#       flaky on some host fallbacks, notably GNU), and every device-data
#       directive becomes an inert no-op comment. Build it with
#       -DRDB_ENABLE_GPU=OFF. Produces the tree that lives on  auto/openmp-cpu.
#
# Runs in the repo root and modifies the working tree IN PLACE:
#   1. Lints `do concurrent` loops (dc_audit.py --strict).
#   2. Translates !$acc directives to !$omp via acc_to_omp.py --write.
#   3. (openmp variant only) Rewrites do-concurrent loops via dc_to_omp.py.
#   4. Adds explicit map() clauses for OPTIONAL array dummies read inside
#      the emitted target regions (omp_optional_map.py).
#   5. Applies every overlay patch under patches/openmp/*/*.patch.
#   6. Flips the CMake default backend to openmp.
#
# Use a clean checkout of main to call this script; do NOT call from a branch
# that already has uncommitted source changes.
#
# Usage:
#   bash tools/regen_dc_openmp.sh                 # dc-openmp-target (default)
#   bash tools/regen_dc_openmp.sh dc-openmp-target
#   bash tools/regen_dc_openmp.sh openmp
#   bash tools/regen_dc_openmp.sh openmp-cpu
#
# An optional second argument sets the worker-process count for the per-file
# passes (0 = one per available core; default 1 = serial):
#   bash tools/regen_dc_openmp.sh openmp 0
#   RDB_REGEN_JOBS=32 bash tools/regen_dc_openmp.sh openmp
#
# Exits non-zero if any stage fails. Intended to be called from CI and
# (occasionally) by developers wanting to test a variant locally.

set -euo pipefail

variant="${1:-dc-openmp-target}"
# Worker processes for the per-file passes. Second positional arg, or
# RDB_REGEN_JOBS. 0 = one per available core (respects SLURM binding);
# 1 (default) = serial, matching the historical behaviour.
jobs="${2:-${RDB_REGEN_JOBS:-1}}"
omp_target=gpu          # acc_to_omp / dc_to_omp --target (gpu offload vs cpu host)
case "$variant" in
  dc-openmp-target) run_dc_to_omp=0 ;;
  openmp)           run_dc_to_omp=1 ;;
  openmp-cpu)       run_dc_to_omp=1; omp_target=cpu ;;
  *)
    echo "error: unknown variant '$variant'" >&2
    echo "       expected 'dc-openmp-target', 'openmp' or 'openmp-cpu'" >&2
    exit 2
    ;;
esac

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

echo "==> variant: $variant  (translator --target $omp_target, -j $jobs)"

echo "==> [1/6] dc_audit --strict"
python tools/dc_audit.py --strict -j "$jobs" src/ app/ benchmarks/ tests/

echo "==> [2/6] acc_to_omp --write"
# acc_to_omp.py walks the directories given on the command line. Pass src/,
# app/, benchmarks/ and tests/ so directives in app/main.F90, the bench
# drivers (bench_ocean's data + update self) and the unit tests are all
# translated — otherwise those keep !$acc, which is inert on an OpenMP build.
python tools/acc_to_omp.py --target "$omp_target" -j "$jobs" --write src/ app/ benchmarks/ tests/

if [[ "$run_dc_to_omp" -eq 1 ]]; then
  echo "==> [3/6] dc_to_omp --write (full OpenMP: rewrite do concurrent loops)"
  # Runs AFTER acc_to_omp: the pure-strip cascade also sees the worksharing
  # regions acc_to_omp emitted from !$acc parallel loop (gpu: !$omp target;
  # cpu: !$omp parallel do), so pure is stripped consistently across both
  # directive sources.  Same dir set as above so the whole-program pure cascade
  # spans benchmarks/ and tests/ too.
  python tools/dc_to_omp.py --target "$omp_target" -j "$jobs" --write src/ app/ benchmarks/ tests/
else
  echo "==> [3/6] dc_to_omp SKIPPED (dc-openmp-target keeps do concurrent)"
fi

if [[ "$omp_target" == "gpu" ]]; then
  echo "==> [4/6] omp_optional_map --write (guard absent OPTIONAL dummies)"
  # Runs AFTER both translators so it covers target regions from EITHER source
  # (dc_to_omp's plain loops and acc_to_omp's reduction loops). LLVM Flang maps
  # an OPTIONAL explicit-shape array dummy read inside a target region
  # unconditionally: bounds from the declaration, base address NULL when the
  # argument is absent, so libomptarget tries hsa_amd_memory_lock(0x0, size)
  # and the run aborts. Naming it in an explicit map() clause takes the
  # presence-guarded path.
  python tools/omp_optional_map.py -j "$jobs" --write src/ app/ benchmarks/ tests/
else
  echo "==> [4/6] omp_optional_map SKIPPED (cpu target emits no map clauses)"
fi

echo "==> [5/6] apply overlay patches"
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

echo "==> [6/6] flip CMake default backend openacc -> openmp"
# In-place sed; preserves indentation. The CMakeLists option block uses double
# quotes around the default value so the substitution is precise.
if grep -q '^    "openacc"$' CMakeLists.txt; then
  # macOS sed differs from GNU sed in -i semantics; pipe-through-temp is portable
  tmp=$(mktemp)
  awk '
    /set\(RDB_PARALLEL_BACKEND$/ { in_block = 1 }
    in_block && /^    "openacc"$/ { sub("openacc", "openmp"); in_block = 0 }
    { print }
  ' CMakeLists.txt > "$tmp"
  mv "$tmp" CMakeLists.txt
  echo "    done"
else
  echo "    SKIP — CMakeLists.txt default backend already not 'openacc'"
fi

echo "==> regen complete ($variant)"
echo ""
echo "Next steps (manual, for local validation):"
if [[ "$omp_target" == "cpu" ]]; then
  echo "  cmake -B build_omp -S . -DRDB_PARALLEL_BACKEND=openmp \\"
  echo "        -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_THREADS=ON"
else
  echo "  cmake -B build_omp -S . -DRDB_PARALLEL_BACKEND=openmp -DRDB_ENABLE_GPU=OFF"
fi
echo "  cmake --build build_omp -j"
