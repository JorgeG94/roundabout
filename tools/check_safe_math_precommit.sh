#!/usr/bin/env bash
# Pre-commit wrapper for tools/check_safe_math.py (diff-aware).
# Base: $SAFE_MATH_BASE, else the first of origin/beta, origin/main, main that
# exists.  No base (shallow clone, fresh worktree): skip cleanly.
set -euo pipefail
base_ref="${SAFE_MATH_BASE:-}"
if [ -z "${base_ref}" ]; then
    for cand in origin/beta origin/main main; do
        if git rev-parse --verify --quiet "${cand}^{commit}" >/dev/null 2>&1; then
            base_ref="${cand}"
            break
        fi
    done
fi
if [ -z "${base_ref}" ]; then
    echo "safe-math-intrinsics: no base ref in this checkout -- skipping."
    exit 0
fi
exec python3 tools/check_safe_math.py --diff "${base_ref}"
