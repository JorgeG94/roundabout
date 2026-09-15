"""pytest fixtures for the rdb Python package tests.

Adds ``python/`` to sys.path (so ``import rdb`` works from a checkout
with no ``pip install``) and skips the whole session with a clear message
if the shared library cannot be located/loaded -- this suite exercises the
real ``.so``, not a mock.
"""

import os
import sys

import pytest

_PY_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PY_DIR not in sys.path:
    sys.path.insert(0, _PY_DIR)

try:
    import rdb
except Exception as exc:  # pragma: no cover - import-time environment issue
    rdb = None
    _import_error = exc
else:
    _import_error = None


def pytest_configure(config):
    if rdb is None:
        pytest.exit(
            f"rdb package failed to import: {_import_error!r}. "
            f"Set RDB_LIB to librdb_core.so, or build it with "
            f"`cmake -B build_shared -S . -DRDB_BUILD_SHARED=ON && "
            f"cmake --build build_shared`.",
            returncode=0,
        )


@pytest.fixture
def ocean_nml():
    """Minimal valid ocean namelist -- 8x6x3, matches the Fortran ctest
    reference (tests/test_ocean_api_p2.F90's `ocean_nml()`), so a failure
    here is directly comparable against the Fortran-side gate."""
    return (
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /\n"
        "&nonhydrostatic_nml nz_layers = 3 /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        '&ocean_topo_nml topo_config = "flat", max_depth = 200.0 /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        "&ocean_diag_nml enabled = .false. /\n"
        "&output_nml output_to_file = .false. /\n"
    )


@pytest.fixture
def bad_nml():
    """A malformed namelist -- unknown key, strict schema parse failure."""
    return "&grid_nml nx = 8, ny = 6, totally_bogus_key = 1 /\n"


@pytest.fixture(autouse=True)
def _no_leaked_model():
    """Guard against a test leaking a live handle into the next test (the
    single-live-handle guard would otherwise turn an unrelated later test
    into a confusing AlreadyExistsError)."""
    yield
    import rdb._model as _m
    if _m._live_model_ref is not None:
        m = _m._live_model_ref()
        if m is not None:
            m.close()
