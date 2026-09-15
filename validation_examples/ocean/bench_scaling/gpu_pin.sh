#!/bin/bash
# One GPU per local rank, set BEFORE MPI_Init (UCX/hpcx primary-context
# gotcha: without this every rank also lands a context on device 0).
export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
exec "$@"
