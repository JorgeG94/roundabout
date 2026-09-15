"""End-to-end pytest coverage for the two "resting state must stay at
rest" validation examples ported to the Python API:
``validation_examples/ocean/acc_channel/acc_channel_quiescent.py`` and
``validation_examples/ocean/seamount/seamount.py``.

These are cheap (small grids, short integration here) and strong
correctness traps: a genuinely at-rest configuration (uniform T/S,
zero forcing) that drifts is a real dyn-core bug, not noise -- see each
script's own module docstring for the physics. The full-duration
thresholds documented there (5 days / 30 days, max|u|,|v| ~1e-10 to
1e-11 m/s) are the ones actually asserted by each script's own
``main()``; this test runs a much SHORTER window (a handful of steps,
via ``model.run(steps=...)`` rather than the scripts' own
``until=``) so the suite stays fast, while exercising the exact same
``build_model()`` construction + assertion logic end-to-end (create,
step, inspect fields, close). The scripts' own byte-identity-vs-.nml
config check is a separate, non-pytest gate (see the port report /
``tmp_local_artifacts/compare_config.py``) -- this file only checks
behaviour: does the Python-built model actually create and stay quiet.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_ACC_DIR = _REPO_ROOT / "validation_examples" / "ocean" / "acc_channel"
_SEAMOUNT_DIR = _REPO_ROOT / "validation_examples" / "ocean" / "seamount"


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def acc_channel_quiescent_module():
    return _load_module(
        _ACC_DIR / "acc_channel_quiescent.py", "_val_acc_channel_quiescent")


@pytest.fixture(scope="module")
def seamount_module():
    return _load_module(_SEAMOUNT_DIR / "seamount.py", "_val_seamount")


def _max_abs(field):
    """Max |value| over a whole 3D Field via ``.copy()`` -- stdlib only,
    no numpy dependency (``Field.__getitem__``'s partial-slicing path
    needs numpy; ``.copy()``/full-ellipsis does not -- see D6.6)."""
    data = field.copy()
    return max(abs(v) for row in data for col in row for v in col)


def _max_uv(model):
    return max(_max_abs(model.u), _max_abs(model.v))


def test_acc_channel_quiescent_stays_at_rest(acc_channel_quiescent_module):
    """Spherical re-entrant channel, seamount ridge, EPBL + kappa-shear
    ON, zero wind, zero heat flux, uniform T/S. A short window here
    (20 steps of the model's own dt=1200s, ~6.7 hours simulated) --
    the script's own ``main()`` asserts the full 5-day threshold."""
    model = acc_channel_quiescent_module.build_model()
    with model:
        model.run(steps=20)
        max_speed = _max_uv(model)
    assert max_speed < 1.0e-6, (
        f"acc_channel_quiescent drifted: max|u|,|v| = {max_speed:.3e} m/s "
        f"after 20 steps (expect ~roundoff for a genuinely at-rest IC)")


def test_seamount_stays_at_rest(seamount_module):
    """Closed-basin Gaussian seamount, uniform rho, f=0, zero forcing,
    zero viscosity/diffusivity -- the cleanest possible BPG/Coriolis-
    adv/continuity trap. Short window here (20 steps of dt=300s, 100
    minutes simulated) -- the script's own ``main()`` asserts the full
    30-day threshold."""
    model = seamount_module.build_model()
    with model:
        model.run(steps=20)
        max_speed = _max_uv(model)
    assert max_speed < 1.0e-6, (
        f"seamount drifted: max|u|,|v| = {max_speed:.3e} m/s after 20 "
        f"steps (expect ~roundoff -- uniform rho + sloping bathymetry "
        f"should cancel exactly in the BPG)")


def test_acc_channel_quiescent_matches_nml_config(acc_channel_quiescent_module):
    """The config-identity gate, as a pytest (in addition to the
    stand-alone tmp_local_artifacts/compare_config.py harness used
    during porting): every knob acc_channel_quiescent.nml sets
    explicitly must be reproduced (explicitly, or via an equal Fortran
    default) by the curated build_model()."""
    from rdb._config import ALL_GROUPS
    from rdb._nml_parse import config_from_namelist_text, _strip_comments

    model = acc_channel_quiescent_module.build_model()
    py_explicit = dict(
        ((g, k), v) for g, k, v in model.config.explicit_knobs())

    nml_path = _ACC_DIR / "acc_channel_quiescent.nml"
    # See compare_config.py's docstring: _strip_comments works around a
    # real bug in parse_namelist_text's group-boundary scanner tripped
    # by an English contraction apostrophe inside a comment.
    nml_text = _strip_comments(nml_path.read_text())
    nml_cfg = config_from_namelist_text(nml_text)
    name_to_attr = {cls._nml_name: attr for attr, cls in ALL_GROUPS.items()}

    for g, k, want in nml_cfg.explicit_knobs():
        attr = name_to_attr[g]
        if (g, k) in py_explicit:
            got = py_explicit[(g, k)]
        else:
            got = getattr(ALL_GROUPS[attr], k).default
        assert got == want, f"&{g}_nml {k}: nml={want!r} python={got!r}"


def test_seamount_matches_nml_config(seamount_module):
    from rdb._config import ALL_GROUPS
    from rdb._nml_parse import config_from_namelist_text, _strip_comments

    model = seamount_module.build_model()
    py_explicit = dict(
        ((g, k), v) for g, k, v in model.config.explicit_knobs())

    nml_path = _SEAMOUNT_DIR / "seamount.nml"
    nml_text = _strip_comments(nml_path.read_text())
    nml_cfg = config_from_namelist_text(nml_text)
    name_to_attr = {cls._nml_name: attr for attr, cls in ALL_GROUPS.items()}

    for g, k, want in nml_cfg.explicit_knobs():
        attr = name_to_attr[g]
        if (g, k) in py_explicit:
            got = py_explicit[(g, k)]
        else:
            got = getattr(ALL_GROUPS[attr], k).default
        # Documented, verified-equivalent divergence (see the port
        # report): RectilinearGrid always routes through
        # &ocean_grid_nml len_lon/len_lat (letting the Fortran derive
        # dx/dy), never through &grid_nml dx/dy directly -- both this
        # .nml's dx=2000.0/dy=2000.0 and the curated grid's
        # len_lon=len_lat=100000.0 describe the SAME 50-cell/100 km
        # grid (100000.0 / 50 == 2000.0).
        if (g, k) in (("grid", "dx"), ("grid", "dy")):
            continue
        assert got == want, f"&{g}_nml {k}: nml={want!r} python={got!r}"
