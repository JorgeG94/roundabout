"""P4 rule: the generated config layer never imports numpy (P3/P4 rule --
numpy exists only in the conda dev env, `pip install` is forbidden, and
`import rdb` must work from a plain checkout with no Fortran build).

Two checks: a static source scan (always runs, catches an accidental
`import numpy` at review time) and, when a numpy-free system interpreter
can be found, an actual subprocess import under it -- this is what
directly demonstrates "imports on bare system python3 (3.12.3, no
numpy)" rather than merely inferring it from this suite's own
interpreter, which DOES have numpy in some dev environments.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "python" / "rdb"

_CONFIG_MODULES = [
    "_knob.py", "_config.py", "_config_generated.py",
]


@pytest.mark.parametrize("filename", _CONFIG_MODULES)
def test_no_numpy_import_in_source(filename):
    text = (PKG / filename).read_text(encoding="utf-8")
    for line in text.splitlines():
        stripped = line.strip()
        assert not stripped.startswith("import numpy"), (
            f"{filename}: {stripped!r} -- numpy is forbidden in the P4 "
            f"config layer (see this package's design constraint D6.6)")
        assert "from numpy" not in stripped, (
            f"{filename}: {stripped!r} -- numpy is forbidden in the P4 "
            f"config layer (see this package's design constraint D6.6)")


def _bare_system_python():
    """A python3 interpreter that is NOT the one running this test suite
    (which may have numpy) and does NOT have numpy importable -- as close
    as this environment gets to the deployment target this rule protects
    ("bare system python3, no numpy")."""
    candidates = ["/usr/bin/python3", "python3"]
    for cand in candidates:
        if cand == sys.executable:
            continue
        if importlib.util.find_spec is None:
            continue
        probe = subprocess.run(
            [cand, "-c", "import numpy"], capture_output=True)
        has_numpy = probe.returncode == 0
        version = subprocess.run(
            [cand, "-c", "import sys; print(sys.executable)"],
            capture_output=True, text=True)
        if version.returncode == 0 and not has_numpy:
            return cand
    return None


def test_config_layer_imports_on_bare_system_python():
    interp = _bare_system_python()
    if interp is None:
        pytest.skip(
            "no numpy-free system python3 found distinct from the test "
            "runner's own interpreter -- cannot demonstrate the bare-"
            "system-python import independently of this suite's env")
    result = subprocess.run(
        [interp, "-c",
         "import sys; sys.path.insert(0, %r); import rdb; "
         "cfg = rdb.Config(); "
         "assert cfg.to_namelist() == ''; "
         "print('OK', rdb.Config)" % str(ROOT / "python")],
        capture_output=True, text=True)
    assert result.returncode == 0, (
        f"rdb failed to import under {interp}:\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}")
    assert "OK" in result.stdout
