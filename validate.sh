#!/usr/bin/env bash
# -----------------------------------------------------------------------
# validate.sh — Build and test Roundabout across all supported backends
#
# Backends:
#   gcc-serial          gfortran, no acceleration, no MPI
#   gcc-mpi             gfortran, no acceleration, MPI (4 ranks)
#   nvhpc-multicore     nvfortran -stdpar=multicore -mp=multicore, no MPI
#   nvhpc-gpu           nvfortran -stdpar=gpu -mp=gpu, no MPI
#   nvhpc-gpu-mpi       nvfortran -stdpar=gpu -mp=gpu, MPI (4 ranks)
#
# Usage:
#   ./validate.sh                  # run all backends
#   ./validate.sh gcc-serial       # run one backend
#   ./validate.sh gcc-serial nvhpc-gpu  # run selected backends
#   JOBS=8 ./validate.sh           # override parallelism
#
# Toolchain: load gfortran / NVHPC + NetCDF however your site does it
# (module load, Spack -- see environments/ --, conda) before running.
# See "Environment loaders" below.
# -----------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_ROOT="${SCRIPT_DIR}/build_validate"
RESULTS_FILE="${BUILD_ROOT}/results.txt"
JOBS="${JOBS:-$(nproc)}"

# -----------------------------------------------------------------------
# Environment loaders
#
# Roundabout ships no machine-specific module scripts: load your toolchain
# however your site does it (`module load`, Spack -- see environments/ --,
# conda) BEFORE running this.  Each backend runs in its own subshell, so a
# CPU and a GPU toolchain never end up stacked in one shell -- which is the
# point: two NetCDF builds on one link line fail in confusing ways.
#
# Local dev convenience: on a configured dev box that keeps untracked
# site-specific `gcc_env.sh` / `nvhpc_env.sh` module scripts at the repo
# root, export RDB_ON_DEV=1 and they are sourced per backend.  A silent
# no-op in a clean checkout.
# -----------------------------------------------------------------------
load_dev_env() {
    # load_dev_env <script-name> -- source a local dev env script, if enabled.
    local script="${SCRIPT_DIR}/$1"
    if [[ "${RDB_ON_DEV:-0}" == "1" && -f "${script}" ]]; then
        # shellcheck disable=SC1090
        source "${script}"
    fi
}

load_gcc_env() {
    load_dev_env gcc_env.sh
}

load_nvhpc_env() {
    load_dev_env nvhpc_env.sh
}

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------
log()  { printf "\033[1;34m==>\033[0m %s\n" "$*"; }
ok()   { printf "\033[1;32m  PASS\033[0m %s\n" "$1"; }
fail() { printf "\033[1;31m  FAIL\033[0m %s\n" "$1"; }
skip() { printf "\033[1;33m  SKIP\033[0m %s\n" "$1"; }

record() {
    # record <backend> <result>
    echo "$1=$2" >> "${RESULTS_FILE}"
}

has_gpu() {
    nvidia-smi --query-gpu=name --format=csv,noheader &>/dev/null
}

run_backend() {
    local name="$1"
    shift
    local cmake_args=("$@")
    local build_dir="${BUILD_ROOT}/${name}"

    # Always start each backend from a clean build dir.  Incremental
    # rebuilds across `validate.sh` invocations have repeatedly failed in
    # subtle ways -- in particular, NVHPC's pgacclnk on RDB_BUILD_SHARED
    # links can pick up the previous-run libcore_rdb.so as both input
    # and output:
    #     /usr/bin/ld: input file 'libcore_rdb.so' is the same as output file
    # The validate flow's contract is "reproducible from scratch", so
    # wiping the dir is the right invariant.  Set NO_CLEAN=1 to skip.
    if [[ "${NO_CLEAN:-0}" != "1" ]]; then
        rm -rf "${build_dir}"
    fi

    log "${name}: configuring"
    if ! cmake -B "${build_dir}" -S "${SCRIPT_DIR}" \
         -DCMAKE_BUILD_TYPE=Release \
         "${cmake_args[@]}" 2>&1 | tail -5; then
        fail "${name}: configure failed"
        record "${name}" "FAIL (configure)"
        return 1
    fi

    log "${name}: building (${JOBS} jobs)"
    if ! cmake --build "${build_dir}" -j"${JOBS}" 2>&1 | tail -5; then
        fail "${name}: build failed"
        record "${name}" "FAIL (build)"
        return 1
    fi

    log "${name}: running tests"
    local ctest_output
    ctest_output=$(cd "${build_dir}" && ctest -R "rdb" --output-on-failure 2>&1) || true
    echo "${ctest_output}" | tail -15

    local summary
    summary=$(echo "${ctest_output}" | grep "tests passed" || echo "")
    if echo "${summary}" | grep -q "100% tests passed"; then
        ok "${name}: ${summary}"
        record "${name}" "PASS"
    else
        fail "${name}: ${summary:-tests failed}"
        record "${name}" "FAIL (tests)"
        return 1
    fi
}

# -----------------------------------------------------------------------
# Backend definitions
# -----------------------------------------------------------------------
backend_gcc_serial() {
    load_gcc_env
    run_backend "gcc-serial" \
        -DCMAKE_Fortran_COMPILER=gfortran \
        -DRDB_ENABLE_GPU=OFF \
        -DRDB_ENABLE_MPI=OFF
}

backend_gcc_mpi() {
    load_gcc_env
    run_backend "gcc-mpi" \
        -DCMAKE_Fortran_COMPILER=mpifort \
        -DRDB_ENABLE_GPU=OFF \
        -DRDB_ENABLE_MPI=ON
}

backend_nvhpc_serial() {
    load_nvhpc_env
    run_backend "nvhpc-serial" \
        -DCMAKE_Fortran_COMPILER=nvfortran \
        -DRDB_ENABLE_GPU=OFF \
        -DRDB_ENABLE_THREADS=OFF \
        -DRDB_ENABLE_MPI=OFF
}

backend_nvhpc_multicore() {
    load_nvhpc_env
    run_backend "nvhpc-multicore" \
        -DCMAKE_Fortran_COMPILER=nvfortran \
        -DRDB_ENABLE_GPU=OFF \
        -DRDB_ENABLE_THREADS=ON \
        -DRDB_ENABLE_MPI=OFF
}

backend_nvhpc_gpu() {
    if ! has_gpu; then
        skip "nvhpc-gpu: no GPU detected"
        record "nvhpc-gpu" "SKIP"
        return 0
    fi
    load_nvhpc_env
    run_backend "nvhpc-gpu" \
        -DCMAKE_Fortran_COMPILER=nvfortran \
        -DRDB_ENABLE_GPU=ON \
        -DRDB_ENABLE_MPI=OFF \
        -DRDB_GPU_ARCH=cc70
}

backend_nvhpc_gpu_mpi() {
    if ! has_gpu; then
        skip "nvhpc-gpu-mpi: no GPU detected"
        record "nvhpc-gpu-mpi" "SKIP"
        return 0
    fi
    load_nvhpc_env
    run_backend "nvhpc-gpu-mpi" \
        -DCMAKE_Fortran_COMPILER=mpifort \
        -DRDB_ENABLE_GPU=ON \
        -DRDB_ENABLE_MPI=ON \
        -DRDB_GPU_ARCH=cc70
}

# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
ALL_BACKENDS=(
    gcc-serial
    gcc-mpi
    nvhpc-serial
    nvhpc-multicore
    nvhpc-gpu
    nvhpc-gpu-mpi
)

# Parse arguments: if none given, run all
if [[ $# -gt 0 ]]; then
    BACKENDS=("$@")
else
    BACKENDS=("${ALL_BACKENDS[@]}")
fi

log "Roundabout validation — $(date -Iseconds)"
log "Backends: ${BACKENDS[*]}"
log "Build root: ${BUILD_ROOT}"
log "Jobs: ${JOBS}"
echo ""

mkdir -p "${BUILD_ROOT}"
: > "${RESULTS_FILE}"

for backend in "${BACKENDS[@]}"; do
    func_name="backend_${backend//-/_}"
    if declare -f "${func_name}" >/dev/null 2>&1; then
        echo ""
        log "========== ${backend} =========="
        # Run in subshell so environment changes don't leak between backends
        ( "${func_name}" ) || true
    else
        fail "Unknown backend: ${backend}"
        record "${backend}" "UNKNOWN"
    fi
done

# -----------------------------------------------------------------------
# Summary — read results back from file
# -----------------------------------------------------------------------
declare -A RESULTS
PASS=0
FAIL=0
SKIP=0

while IFS='=' read -r key value; do
    [[ -z "${key}" ]] && continue
    RESULTS["${key}"]="${value}"
    case "${value}" in
        PASS*)  ((PASS++)) || true ;;
        FAIL*)  ((FAIL++)) || true ;;
        SKIP*)  ((SKIP++)) || true ;;
    esac
done < "${RESULTS_FILE}"

echo ""
echo ""
log "==================== SUMMARY ===================="
printf "  %-25s %s\n" "Backend" "Result"
printf "  %-25s %s\n" "-------" "------"
for backend in "${BACKENDS[@]}"; do
    result="${RESULTS[${backend}]:-UNKNOWN}"
    case "${result}" in
        PASS*)  printf "  \033[1;32m%-25s %s\033[0m\n" "${backend}" "${result}" ;;
        FAIL*)  printf "  \033[1;31m%-25s %s\033[0m\n" "${backend}" "${result}" ;;
        SKIP*)  printf "  \033[1;33m%-25s %s\033[0m\n" "${backend}" "${result}" ;;
        *)      printf "  %-25s %s\n" "${backend}" "${result}" ;;
    esac
done
echo ""
log "Pass: ${PASS}  Fail: ${FAIL}  Skip: ${SKIP}"
log "================================================="

if [[ ${FAIL} -gt 0 ]]; then
    exit 1
fi
