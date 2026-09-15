"""P4 anti-drift gate: `python/rdb/_config_generated.py` must be an
exact rendering of the CURRENTLY BUILT `nml_schema_t`.

Runs the real pipeline (`rdb_nml_json` -> `tools/gen_python_config.py`)
into a temp file and diffs it against the checked-in file. This is the
whole point of the P4 phase: without this gate, a knob added, renamed, or
retired in Fortran would silently desynchronise from the Python layer,
which is worse than `docs/generated_nml_knobs.md` going stale (that is
merely a stale doc; a stale generated CONFIG LAYER accepts a retired
knob and silently does nothing -- exactly the dead-knob problem this
phase's `dead_on_ocean_path` annotation exists to prevent, but for knobs
Python doesn't even know changed).

Needs the `rdb_nml_json` binary from the SAME build tree as
`RDB_LIB` (both land in the build tree's top-level directory, no
`RUNTIME_OUTPUT_DIRECTORY` override) -- skips (not fails) if it cannot be
found, e.g. a build tree that only shipped the installed `.so`.
"""

import filecmp
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
GEN_SCRIPT = ROOT / "tools" / "gen_python_config.py"
CHECKED_IN = ROOT / "python" / "rdb" / "_config_generated.py"


def _find_nml_json_binary():
    import rdb._ffi as _ffi

    lib_path = None
    try:
        lib_path = Path(_ffi.find_library())
    except _ffi.LibraryNotFoundError:
        pass
    candidates = []
    if lib_path is not None:
        candidates.append(lib_path.parent / "rdb_nml_json")
    for cand in sorted(ROOT.glob("build*/rdb_nml_json")):
        candidates.append(cand)
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def test_config_generated_matches_live_schema(tmp_path):
    binary = _find_nml_json_binary()
    if binary is None:
        pytest.skip(
            "rdb_nml_json not found next to librdb_core.so -- build "
            "it with `cmake --build <build_shared_dir> --target "
            "rdb_nml_json` to enable this drift check")

    schema_json = tmp_path / "schema.json"
    regen = tmp_path / "_config_generated.py"

    result = subprocess.run(
        [str(binary), str(schema_json)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert result.returncode == 0, (
        f"rdb_nml_json failed: {result.stdout}\n{result.stderr}")

    result = subprocess.run(
        [sys.executable, str(GEN_SCRIPT), str(schema_json), str(regen)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert result.returncode == 0, (
        f"gen_python_config.py failed: {result.stdout}\n{result.stderr}")

    if not filecmp.cmp(str(regen), str(CHECKED_IN), shallow=False):
        diff = subprocess.run(
            ["diff", "-u", str(CHECKED_IN), str(regen)],
            capture_output=True, text=True).stdout
        pytest.fail(
            "python/rdb/_config_generated.py has drifted from the live "
            "nml_schema_t. Regenerate with:\n\n"
            f"    {binary} tmp_local_artifacts/schema.json\n"
            f"    python3 {GEN_SCRIPT.relative_to(ROOT)} "
            "tmp_local_artifacts/schema.json "
            "python/rdb/_config_generated.py\n\n"
            f"Diff (checked-in vs regenerated):\n{diff}")
