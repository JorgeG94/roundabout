"""D4 — the exception hierarchy and the error-ring reader.

The problem this exists to solve: every ocean C ABI entry point returns a
coarse STAGE code (``rdb_ocean_status``). The SPECIFIC reason a namelist was
rejected, a tracer name wasn't found, or a bathymetry array normalised to
zero wet cells lives in the Fortran error ring
(``rdb_ocean_last_error``), not in the integer. A bare
``RuntimeError: error 3`` is exactly the "awful" this package exists to
avoid, so ``check()`` below is the ONLY place a status code is allowed to
become an exception, and it always drains the ring first.

See ``tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`` D4
and its sketch ``sketch/04_errors.py`` for the full design rationale (this
module follows both, adapted to drop numpy and to cover every status code
the current ``include/rdb_ocean.h`` defines, including the P2/P2.5
additions the sketch predates: NOT_FOUND(13), BAD_SHAPE(14),
BATHYMETRY_SIGN(15), NOT_PENDING(16)).
"""

from __future__ import annotations

import ctypes

from . import _ffi


class RdbError(Exception):
    """Base for everything this package raises.

    Carries the machine-readable parts (``code``, ``stage``) alongside the
    human message so callers can branch without string matching. ``trace``
    holds the full error-ring dump, deepest-first; ``str(exc)`` is always
    just the deepest (most specific) message.
    """

    code = None
    stage = None

    def __init__(self, message, *, code=None, trace=()):
        super().__init__(message)
        if code is not None:
            self.code = code
        self.trace = tuple(trace)


# --- config / setup, mapped from rdb_ocean_status ----------------------

class RdbConfigError(RdbError):
    stage = "config"


class ConfigParseError(RdbConfigError):
    """OCEAN_STATUS_ERR_CONFIG_PARSE (1). Unknown group/key, type/range/
    enum violation. Carries the schema parser's own `path:line:` prefix
    and did-you-mean suggestion through the ring, verbatim."""
    code = 1


class ConfigValidationError(RdbConfigError):
    """OCEAN_STATUS_ERR_CONFIG_VALIDATE (2). A cross-knob semantic check
    failed (``validate_config``)."""
    code = 2


class RdbSetupError(RdbError):
    """OCEAN_STATUS_ERR_SETUP (3). The `configure_ocean_*` setup chain (or
    a direct callee: metrics/EOS/BC validation) rejected the config."""
    code = 3
    stage = "setup"


class RdbInitialConditionError(RdbError):
    """OCEAN_STATUS_ERR_IC_SEED (4)."""
    code = 4
    stage = "ic"


class RdbIOError(RdbError, OSError):
    """OCEAN_STATUS_ERR_IO (5). A setup-time file read failed. Also an
    OSError so `except OSError` catches a missing forcing file the way an
    outer handler would expect."""
    code = 5
    stage = "io"


# --- handle lifecycle ----------------------------------------------------

class RdbHandleError(RdbError):
    stage = "handle"


class InvalidHandleError(RdbHandleError):
    """OCEAN_STATUS_ERR_BAD_HANDLE (10): null, stale, or garbage handle."""
    code = 10


class AlreadyExistsError(RdbHandleError):
    """OCEAN_STATUS_ERR_ALREADY_EXISTS (11) — one live ocean handle.

    Raised by :mod:`rdb._model` BEFORE the C call whenever possible, so
    the message can name the model that is still open — information the
    Fortran side cannot have.
    """
    code = 11


class NotInitialisedError(RdbHandleError):
    """OCEAN_STATUS_ERR_NOT_INITIALISED (12): a step/query/getter landed on
    a handle that never finished create() (still pending, or bad)."""
    code = 12


class TracerNotFoundError(RdbHandleError):
    """OCEAN_STATUS_ERR_NOT_FOUND (13): no tracer with that name."""
    code = 13


class BadShapeError(RdbHandleError):
    """OCEAN_STATUS_ERR_BAD_SHAPE (14): a setter/stage array's extents do
    not match the live state's physical interior (or the grid, once it
    exists)."""
    code = 14


class NotPendingError(RdbHandleError):
    """OCEAN_STATUS_ERR_NOT_PENDING (16): a stage_*()/create_finalize()
    call landed on a handle that is not currently in the pending window."""
    code = 16


# --- restart (rdb_ocean_status.F90, numbered from 20) -------------------

class RdbRestartError(RdbError):
    stage = "restart"


class RestartSchemaMismatch(RdbRestartError):
    code = 20


class RestartDecompMismatch(RdbRestartError):
    """The restart was written by a different rank decomposition."""
    code = 21


class RestartGridMismatch(RdbRestartError):
    """Grid / vcoord / tracer-registry mismatch against the restart."""
    code = 22


# --- Python-side only ------------------------------------------------

class RdbClosedError(RdbError):
    """A Field or Model was used after close(). Never reaches C."""
    stage = "python"


class RdbReadOnlyError(RdbError):
    stage = "python"


class RdbUnsupportedError(RdbError):
    """A real concept with no Roundabout analogue, or a request this package
    deliberately does not support (partial-slice writes to derived fields
    without numpy, coordinates at an unsupported stagger point, ...)."""
    stage = "python"


class BathymetrySignError(RdbUnsupportedError):
    """OCEAN_STATUS_ERR_BATHYMETRY_SIGN (15). A staged bathymetry array
    normalised to (near-)zero wet cells under the given sign convention --
    the single most dangerous silent failure this API guards against
    (D6.2): passing GEBCO/ETOPO's negative-down convention through
    unflipped makes every cell read as land, giving a clean, crash-free,
    entirely wrong quiescent run."""
    code = 15


class ConfigConflictError(RdbConfigError):
    """Two writers set the same knob to different values. Python-side only
    (P5 curated-object composition); not reachable from P3."""
    stage = "python"


class NoBoundaryLayerWarning(UserWarning):
    """`rdb._compose` warns with this when `closures=[...]` leaves
    NEITHER `KPP` nor `EPBL` active -- e.g. `KPP(enabled=False)` with no
    `EPBL(...)` in the same list. PP81 interior mixing + background
    floors + convective adjustment still run, but nothing represents
    Ekman / wind-driven surface boundary-layer entrainment: a real,
    forced run without either scheme gets an unrealistically shallow
    mixed layer and the wrong near-surface velocity structure.
    Legitimate for a closed-basin, no-wind "pure dyn-core" test (see
    `validation_examples/ocean/seamount/seamount.py`) -- not a mistake
    by construction, hence a warning rather than a raise."""


_BY_CODE = {
    1: ConfigParseError,
    2: ConfigValidationError,
    3: RdbSetupError,
    4: RdbInitialConditionError,
    5: RdbIOError,
    10: InvalidHandleError,
    11: AlreadyExistsError,
    12: NotInitialisedError,
    13: TracerNotFoundError,
    14: BadShapeError,
    15: BathymetrySignError,
    16: NotPendingError,
    20: RestartSchemaMismatch,
    21: RestartDecompMismatch,
    22: RestartGridMismatch,
}

def _read_ring(cap: int = 512, slots: int = 16):
    """Drain rdb_ocean_last_error(i, buf, cap) into a deepest-first list
    of decoded strings. Returns an empty list if the ring has nothing (a
    path that predates the `fail()` conversion, or a handle-lifecycle code
    that never logs)."""
    lib = _ffi.get_lib()
    out = []
    buf = ctypes.create_string_buffer(cap)
    for i in range(slots):
        n = lib.rdb_ocean_last_error(i, buf, cap)
        if n <= 0:
            break
        out.append(buf.raw[:min(n, cap)].split(b"\x00", 1)[0].decode(
            "utf-8", "replace"))
    return out


def check(status: int, context: str = ""):
    """Raise the right exception, with the right (specific) message, or
    return None on success. The ONLY place a status code becomes an
    exception -- one place to fix if the code<->class map ever changes.
    """
    if status == _ffi.OK:
        return
    trace = _read_ring()
    if trace:
        message = trace[0]
    else:
        message = (
            f"{context or 'rdb'} failed with status {status} and left "
            "no diagnostic in the error ring (this path may predate the "
            "fail() conversion -- please report it with the config that "
            "produced it)"
        )
    cls = _BY_CODE.get(status, RdbError)
    raise cls(message, code=status, trace=trace)
