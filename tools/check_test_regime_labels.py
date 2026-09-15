#!/usr/bin/env python3
"""Pre-commit guard: every rdb ctest entry must carry a regime LABEL.

Every test registered in tests/CMakeLists.txt must be attributable to a regime
suite (`ctest -L ocean` / `-L core`).  (The coastal regime was carved out into
its own repository; `coastal` is no longer a legal label — tests/CMakeLists.txt
is the authority and fails configure on it.)
Three things are enforced, all in tests/CMakeLists.txt:

  1. Every RDB_TESTS row has exactly 5 `|`-separated fields and the 5th
     (regime) is one of: ocean / core / ocean+core.
  2. Every test named in MPI_TESTS has a matching
     `set(MPI_TEST_REGIME_<name> <regime>)` with a valid single regime.
  3. Every literally-named `add_test(NAME rdb_...)` (the extra MPI
     rank-count legs) appears in some `set_tests_properties(... LABELS ...)`.

The pic / test-drive dependency self-tests are registered by their own
(vendored) build systems and are deliberately NOT covered here — they are not
rdb tests.  CMake enforces 1-2 fail-loud at configure time as well; this
hook catches it before a configure is ever run.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS_CMAKE = ROOT / "tests" / "CMakeLists.txt"

VALID_REGIMES = {"ocean", "core", "ocean+core"}


def main() -> int:
    text = TESTS_CMAKE.read_text(encoding="utf-8", errors="replace")
    errors: list[str] = []

    # ------------------------------------------------------------------ 1.
    # RDB_TESTS rows: every quoted |-row must have 5 fields, valid regime.
    rows = [m for m in re.findall(r'"(test_[^"]*)"', text) if "|" in m]
    if not rows:
        errors.append("no RDB_TESTS rows found — parser broken or list renamed")
    for row in rows:
        fields = row.split("|")
        if len(fields) != 5:
            errors.append(
                f"RDB_TESTS row '{row}' has {len(fields)} fields; expected 5 "
                "(<test-name>|<module-name>|<collect-fn>|<suite-display>|<regime>)"
            )
        elif fields[4] not in VALID_REGIMES:
            errors.append(
                f"RDB_TESTS row '{fields[0]}': invalid regime '{fields[4]}' "
                f"(must be one of: {', '.join(sorted(VALID_REGIMES))})"
            )

    # ------------------------------------------------------------------ 2.
    # MPI_TESTS entries need a valid MPI_TEST_REGIME_<name>.
    mpi_block = re.search(r"set\(MPI_TESTS\s+([^)]*)\)", text)
    mpi_tests = mpi_block.group(1).split() if mpi_block else []
    regime_sets = dict(
        re.findall(r"set\(MPI_TEST_REGIME_(\w+)\s+(\S+)\)", text)
    )
    for t in mpi_tests:
        regime = regime_sets.get(t)
        if regime is None:
            errors.append(
                f"MPI test '{t}' has no set(MPI_TEST_REGIME_{t} <regime>)"
            )
        elif regime not in VALID_REGIMES - {"ocean+core"}:
            errors.append(
                f"MPI test '{t}': invalid regime '{regime}' "
                "(must be ocean or core)"
            )

    # ------------------------------------------------------------------ 3.
    # Literally-named add_test entries (extra MPI rank-count legs) must be
    # labelled via an explicit set_tests_properties(... LABELS ...).
    labelled: set[str] = set()
    for props in re.findall(r"set_tests_properties\(([^)]*)\)", text):
        if "LABELS" not in props:
            continue
        labelled.update(re.findall(r"\brdb_\w+", props))
    for name in re.findall(r"add_test\(\s*NAME\s+(rdb_\w+)", text):
        # names containing ${...} are handled by the loops above; the regex
        # only matches fully literal names.
        if name not in labelled:
            errors.append(
                f"ctest entry '{name}' is registered with add_test but never "
                "given a regime LABEL via set_tests_properties(... LABELS ...)"
            )

    if errors:
        print("check_test_regime_labels: FAIL")
        for e in errors:
            print(f"  - {e}")
        print(
            "\nEvery rdb ctest entry must carry a regime label so the "
            "ocean / core suites stay separable (ctest -L <regime>)."
        )
        return 1
    print(
        f"check_test_regime_labels: OK "
        f"({len(rows)} RDB_TESTS rows, {len(mpi_tests)} MPI tests)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
