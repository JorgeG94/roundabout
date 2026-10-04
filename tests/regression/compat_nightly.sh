#!/usr/bin/env bash
# The compatibility matrix, every leg, on one machine -- the nightly (and,
# with --t3, the weekly) run that has a GPU.  GitHub's hosted runners have no
# GPU, so .github/workflows/compat-matrix.yml runs the CPU + MPI legs only;
# this script is the GPU leg's home and the local twin of that workflow.
#
#   tests/regression/compat_nightly.sh \
#       --cpu-env  ./gcc_env.sh   --cpu-build build_gfortran \
#       --mpi-build build_gfortran_mpi \
#       --gpu-env  ./nvhpc_env.sh --gpu-build build_cc70 \
#       [--out-dir tmp_local_artifacts/compat_nightly] [--jobs 4] [--gpu-jobs 4]
#       [--gpu-device N] [--t3] [--previous DIR]
#
# The three builds must already exist (CLAUDE.md "Build"; the MPI one with
# -DRDB_ENABLE_MPI=ON, the GPU one with -DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc70).
# Each env file is sourced in its OWN subshell: the gfortran- and the
# nvfortran-built netcdf-fortran share a soname, so the two toolchains must
# never share one LD_LIBRARY_PATH.  `nccopy` (netcdf-c) must be on the CPU
# env's PATH -- the restart and decomposition legs read checkpoints with it.
#
# Legs:  1. CPU  -- checks 1-3 + ENERGY + RESTART on the serial build, DECOMP
#                   (2x2, 4x1) on the MPI build             -> cpu.json
#        2. GPU  -- the SAME cells (`--cells-from cpu.json`: no second fixed
#                   point), checks 1-3 + ENERGY + CROSS_BACKEND against
#                   cpu.json                                -> gpu.json
# Budget (README "Cost"): the GPU leg of the pairwise set fits well inside one
# V100-hour; the t = 3 slice is the weekly job.
# Exit status: non-zero if either leg has a FAIL or an XPASS.
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cpu_env="" cpu_build="build_gfortran" mpi_build="" gpu_env="" gpu_build=""
out_dir="$repo/tmp_local_artifacts/compat_nightly" jobs=4 gpu_jobs=4 gpu_device="" design=pairwise
previous=""
while [ $# -gt 0 ]; do
  case "$1" in
    --cpu-env) cpu_env="$2"; shift 2 ;;
    --cpu-build) cpu_build="$2"; shift 2 ;;
    --mpi-build) mpi_build="$2"; shift 2 ;;
    --gpu-env) gpu_env="$2"; shift 2 ;;
    --gpu-build) gpu_build="$2"; shift 2 ;;
    --out-dir) out_dir="$2"; shift 2 ;;
    --jobs) jobs="$2"; shift 2 ;;
    --gpu-jobs) gpu_jobs="$2"; shift 2 ;;
    --gpu-device) gpu_device="$2"; shift 2 ;;
    --t3) design=t3; shift ;;
    --previous) previous="$2"; shift 2 ;;
    -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$out_dir"
matrix="python3 $here/compat_matrix.py run --design $design"
prev() { [ -n "$previous" ] && [ -f "$previous/$1" ] && echo "--previous $previous/$1"; }

src() { [ -n "$1" ] && echo ". '$1' >/dev/null 2>&1;"; }

status=0
mpi_args=""
[ -n "$mpi_build" ] && mpi_args="--mpi-build-dir $mpi_build"
echo "=== CPU leg ($design): $cpu_build${mpi_build:+ + MPI $mpi_build}"
t0=$(date +%s)
bash -c "$(src "$cpu_env") cd '$repo' && $matrix --build-dir $cpu_build $mpi_args \
    --jobs $jobs --backend cpu-gfortran --out '$out_dir/cpu.json' $(prev cpu.json)" \
    > "$out_dir/cpu.log" 2>&1 || status=1
echo "    $(( $(date +%s) - t0 )) s, exit $status -- $out_dir/cpu.log"

if [ -n "$gpu_build" ]; then
  # GPU leg: one device (the least loaded unless pinned), several cells at a
  # time -- each is a ~1 s, few-hundred-MB process.
  if [ -z "$gpu_device" ] && command -v nvidia-smi >/dev/null; then
    gpu_device=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
                 | sort -t, -k2 -n | head -1 | cut -d, -f1)
  fi
  echo "=== GPU leg ($design): $gpu_build on device ${gpu_device:-default}"
  t0=$(date +%s)
  gstat=0
  bash -c "$(src "$gpu_env") cd '$repo' && \
      ${gpu_device:+CUDA_VISIBLE_DEVICES=$gpu_device} $matrix --build-dir $gpu_build \
      --jobs $gpu_jobs --backend gpu-nvfortran-cc70 --legs energy,cross_backend \
      --reference '$out_dir/cpu.json' --cells-from '$out_dir/cpu.json' --scratch-root '$repo/tmp_local_artifacts/compat_matrix/gpu' \
      --out '$out_dir/gpu.json' $(prev gpu.json)" > "$out_dir/gpu.log" 2>&1 || gstat=1
  echo "    $(( $(date +%s) - t0 )) s, exit $gstat -- $out_dir/gpu.log"
  [ $gstat -ne 0 ] && status=1
fi
exit $status
