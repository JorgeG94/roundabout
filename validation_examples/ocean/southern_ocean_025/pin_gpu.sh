#!/bin/bash
# One GPU per MPI rank -- CUDA_VISIBLE_DEVICES must be set BEFORE MPI_Init,
# which under CUDA-aware MPI (RDB_CUDA_AWARE_MPI=ON) creates a CUDA primary
# context immediately, before the app's own acc_set_device_num call runs.
# See CLAUDE.md, "Multi-GPU one node".
export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK
exec "$@"
