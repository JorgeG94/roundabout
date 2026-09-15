#!/usr/bin/env bash
# Pre-commit hook: ban `use mpi` / `use mpi_f08` anywhere in Roundabout.
# All MPI calls go through the `pic_mpi_lib` wrapper from the pic
# dependency — including the `src/comm/` backend.  That wrapper is
# the single layer that knows about MPI; everything above it (comm
# facade, solvers, kernels) stays portable across single-process and
# multi-rank builds.
#
# See src/core/ocean/README.md (Design contract rule 4) for the
# rationale; the rule applies to the coastal path too.
#
# Receives candidate files as args from pre-commit.  Can also be run
# manually: `bash tools/no_mpi_in_rdb.sh <files...>`.
set -euo pipefail

if [[ $# -eq 0 ]]; then
    exit 0
fi

violations=$(grep -EnH '^[[:space:]]*use[[:space:]]+(mpi|mpi_f08)\b' "$@" || true)

if [[ -n "$violations" ]]; then
    echo "$violations"
    echo
    echo "ERROR: direct 'use mpi' / 'use mpi_f08' is banned everywhere in Roundabout."
    echo "All MPI calls go through 'pic_mpi_lib' — including inside src/comm/."
    echo "See src/core/ocean/README.md (Design contract rule 4)."
    exit 1
fi

exit 0
