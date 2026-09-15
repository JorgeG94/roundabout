"""Roundabout ocean C-grid dyn-core -- Python runtime driver (P3).

Stdlib only (no numpy; see ``_ffi.py``'s module docstring). A handle-based
driver over the ocean C ABI (``include/rdb_ocean.h``): build a config
(namelist string, for now), create a solver, step it, read and write state
arrays, destroy it.

    import rdb
    with rdb.Model(namelist_text) as model:
        model.step(10)
        print(model.time, model.total_mass)
        sst = model.temperature[:, :, -1]   # surface layer, degC

See ``docs/ocean_python_api_plan.md`` and
``tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`` for the
full design. The typed generated config layer (P4) and curated physics
objects (P5) are NOT part of this phase.
"""

from ._config import Config
from ._errors import (AlreadyExistsError, BadShapeError,
                       BathymetrySignError, ConfigConflictError,
                       ConfigParseError, ConfigValidationError,
                       InvalidHandleError, NoBoundaryLayerWarning,
                       NotInitialisedError, NotPendingError,
                       RdbClosedError, RdbConfigError, RdbError,
                       RdbHandleError, RdbIOError,
                       RdbInitialConditionError, RdbReadOnlyError,
                       RdbRestartError, RdbSetupError,
                       RdbUnsupportedError, RestartDecompMismatch,
                       RestartGridMismatch, RestartSchemaMismatch,
                       TracerNotFoundError)
from ._ffi import (BATHY_DEPTH_POSITIVE_DOWN, BATHY_HEIGHT_POSITIVE_UP,
                    LibraryNotFoundError, available_diagnostics,
                    find_library, get_lib)
from ._field import DerivedField, DiagnosticField, Field, Loc
from ._knob import RdbDeadKnobWarning
from ._model import Model
from .diagnostics import Diagnostics, Duration, Restart

__all__ = [
    "Model", "Field", "DerivedField", "DiagnosticField", "Loc",
    "Config", "RdbDeadKnobWarning",
    "Diagnostics", "Duration", "Restart", "available_diagnostics",
    "find_library", "get_lib", "LibraryNotFoundError",
    "BATHY_DEPTH_POSITIVE_DOWN", "BATHY_HEIGHT_POSITIVE_UP",
    "RdbError", "RdbConfigError", "ConfigParseError",
    "ConfigValidationError", "RdbSetupError",
    "RdbInitialConditionError", "RdbIOError", "RdbHandleError",
    "InvalidHandleError", "AlreadyExistsError", "NotInitialisedError",
    "TracerNotFoundError", "BadShapeError", "NotPendingError",
    "RdbRestartError", "RestartSchemaMismatch", "RestartDecompMismatch",
    "RestartGridMismatch", "RdbClosedError", "RdbReadOnlyError",
    "RdbUnsupportedError", "BathymetrySignError", "ConfigConflictError",
    "NoBoundaryLayerWarning",
]


def working_precision() -> int:
    """Bytes per working-precision float in the loaded build (4 or 8)."""
    return get_lib().wp_bytes
