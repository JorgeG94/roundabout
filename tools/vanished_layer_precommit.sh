#!/usr/bin/env bash
# Pre-commit wrapper around the vanished-layer lint.
#
# The lint is diff-aware (only NEW offenders fail), which needs a base ref.
# Locally that is origin/main; some checkouts don't have it (shallow CI
# clones, fresh worktrees).  In that case resolve the first available base
# ref.  If none is available, skip cleanly.
#
# Mirror of tools/dc_assumed_shape_precommit.sh.
set -euo pipefail

base_ref=""
for cand in origin/main main; do
    if git rev-parse --verify --quiet "${cand}^{commit}" >/dev/null 2>&1; then
        base_ref="${cand}"
        break
    fi
done

if [ -z "${base_ref}" ]; then
    echo "vanished-layer: no base ref (origin/main) in this checkout — skipping."
    exit 0
fi

exec python3 tools/vanished_layer_lint.py src/ app/ --diff "${base_ref}"
