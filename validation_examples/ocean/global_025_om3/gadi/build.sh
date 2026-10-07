#!/bin/bash
# Configure + build roundabout for the global 0.25 degree OM3 2-node / 8-GPU
# (H200) scaling test. Run on a Gadi login or compute node with the NVHPC
# module stack loaded -- module names below are PLAUSIBLE Gadi names, not
# verified from here (no Gadi access); run `module avail nvhpc netcdf hdf5`
# first and fix names before trusting this script. Mirrors
# ../../southern_ocean_025/README.md Sec. 1/5 (same HPC-X/openmpi module
# collision gotcha -- CLAUDE.md "Multi-GPU one node").
set -euo pipefail

cd "$(dirname "$0")/../../../.."   # repo root

# check: confirm these module names/versions with `module avail` on Gadi --
# copied from southern_ocean_025/README.md Sec. 1, which used nvhpc 25.5.
module unload openmpi/5.0.5 2>/dev/null || true
module load nvhpc
module load misc/nvhpc-build/25.5/netcdf-c misc/nvhpc-build/25.5/netcdf-fortran
module load hdf5

BUILD_DIR=build_nvhpc_cc90_mpi

cmake -B "$BUILD_DIR" -S . -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_C_COMPILER=nvc \
    -DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc90 \
    -DRDB_ENABLE_MPI=ON -DRDB_CUDA_AWARE_MPI=ON

cmake --build "$BUILD_DIR" -j 8

echo "Built: $BUILD_DIR/rdb"
echo "check: confirm 'GPU offload:' and MPI lines in the cmake configure"
echo "       log above say ON as expected before submitting run_2node.pbs."
