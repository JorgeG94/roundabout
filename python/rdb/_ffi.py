"""ctypes bindings to librdb_core.so's ocean C ABI (include/rdb_ocean.h).

Stdlib only — no numpy, no third-party packages (project rule: never `pip
install`; system python3 has no numpy at all). Arrays cross this boundary as
raw ``ctypes`` pointers; ``_field.py`` turns them into Python-usable views
(including an ``__array_interface__`` dict that numpy can consume without
this module ever importing it).

Library discovery order (first hit wins), per the P3 spec:

  1. ``RDB_LIB`` environment variable — explicit path to the .so.
  2. The build tree: ``build_shared/librdb_core.so`` next to the repo
     root (also probes any other ``build*/`` directory, so a differently
     named shared build is still found).
  3. The normal dynamic-linker search path (``ctypes.CDLL("librdb_core.so")``
     with no directory component — respects ``LD_LIBRARY_PATH`` etc).

Every prototype in ``include/rdb_ocean.h`` is declared here so a caller
gets a `TypeError`/`ArgumentError` from ctypes on a mismatched call rather
than silent memory corruption.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

# ---------------------------------------------------------------------
# Status codes -- must stay byte-for-byte in sync with
# include/rdb_ocean.h and src/core/ocean/state/rdb_ocean_status.F90.
# ---------------------------------------------------------------------
OK = 0
ERR_CONFIG_PARSE = 1
ERR_CONFIG_VALIDATE = 2
ERR_SETUP = 3
ERR_IC_SEED = 4
ERR_IO = 5
ERR_BAD_HANDLE = 10
ERR_ALREADY_EXISTS = 11
ERR_NOT_INITIALISED = 12
ERR_NOT_FOUND = 13
ERR_BAD_SHAPE = 14
ERR_BATHYMETRY_SIGN = 15
ERR_NOT_PENDING = 16
ERR_RESTART_SCHEMA = 20
ERR_RESTART_DECOMP = 21
ERR_RESTART_GRID = 22

# Bathymetry sign conventions (rdb_ocean_stage_bathymetry).
BATHY_DEPTH_POSITIVE_DOWN = 1
BATHY_HEIGHT_POSITIVE_UP = 2

_LIB_NAME = "librdb_core.so"


class LibraryNotFoundError(FileNotFoundError):
    """Raised when no usable librdb_core.so can be located."""


def _candidate_paths():
    env = os.environ.get("RDB_LIB")
    if env:
        yield Path(env)

    # python/rdb/_ffi.py -> repo root is two parents up.
    repo_root = Path(__file__).resolve().parents[2]
    yield repo_root / "build_shared" / _LIB_NAME
    for cand in sorted(repo_root.glob(f"build*/{_LIB_NAME}")):
        yield cand
    for cand in sorted(repo_root.glob(f"build*/*/{_LIB_NAME}")):
        yield cand


def _loader_search():
    """Last-resort lookup via the normal dynamic-linker search path
    (LD_LIBRARY_PATH / ldconfig cache / rpath). Imported lazily so module
    import itself never touches the filesystem beyond the candidates above.
    """
    import ctypes.util

    return ctypes.util.find_library("rdb_core")


def find_library() -> str:
    """Locate librdb_core.so. Raises LibraryNotFoundError with a clear
    account of every location that was tried."""
    tried = []
    for cand in _candidate_paths():
        tried.append(str(cand))
        if cand.is_file():
            return str(cand)

    found = _loader_search()
    if found:
        return found

    raise LibraryNotFoundError(
        "Could not find librdb_core.so. Looked for:\n"
        + "\n".join(f"  - {t}" for t in tried)
        + "\n  - the normal dynamic-linker search path (LD_LIBRARY_PATH)\n"
        "Set the RDB_LIB environment variable to the absolute path of "
        "librdb_core.so, or build it with:\n"
        "  cmake -B build_shared -S . -DRDB_BUILD_SHARED=ON && "
        "cmake --build build_shared"
    )


class _Lib:
    """Lazily-loaded, module-singleton wrapper around the CDLL handle."""

    def __init__(self):
        self._cdll = None
        self._path = None
        self.wp_bytes = None
        self.wp_ctype = None
        self.wp_typestr = None

    def ensure_loaded(self):
        if self._cdll is not None:
            return
        path = find_library()
        self._cdll = ctypes.CDLL(path)
        self._path = path
        _declare_prototypes(self._cdll)
        self.wp_bytes = self._cdll.rdb_working_precision()
        if self.wp_bytes == 4:
            self.wp_ctype = ctypes.c_float
            self.wp_typestr = "<f4"
        elif self.wp_bytes == 8:
            self.wp_ctype = ctypes.c_double
            self.wp_typestr = "<f8"
        else:
            raise RuntimeError(
                f"rdb_working_precision() returned {self.wp_bytes}, "
                "expected 4 or 8"
            )

    @property
    def path(self):
        self.ensure_loaded()
        return self._path

    def __getattr__(self, name):
        self.ensure_loaded()
        return getattr(self._cdll, name)


_lib = _Lib()


def get_lib() -> _Lib:
    """Return the process-wide library singleton, loading it on first use."""
    _lib.ensure_loaded()
    return _lib


def _read_catalog(size_fn, name_fn, cap=64):
    """Shared tail for the four P7 catalog getters: call `size_fn()`,
    then `name_fn(i, buf, cap)` for each index, decoding the returned
    name. `size_fn`/`name_fn` are bound ctypes functions (module-level,
    no handle)."""
    names = []
    buf = ctypes.create_string_buffer(cap)
    n = size_fn()
    for i in range(n):
        length = name_fn(i, buf, cap)
        names.append(buf.raw[:length].decode("utf-8"))
    return names


def available_diagnostics():
    """Every diagnostic name reachable via `&ocean_diag_nml diags` — the
    CANONICAL catalog (``SSH``, ``temperature``, ``u``, ...) plus the
    DERIVED catalog (``vorticity_z``, ``ke_total``, ``mld_density``,
    ...). Static: no live :class:`~rdb.Model` is required, so this
    can be called before ``create()`` — the biggest remaining
    discoverability gap the P7 design set out to close (a user
    previously had no way to learn what they could ask for without
    reading source).

    Some canonical names are conditionally registered by another
    namelist group (e.g. ``temperature``/``salinity`` need
    ``&ocean_thermo_nml enable_thermodynamics``) — this is "what CAN be
    requested" on SOME configuration, not "what THIS run has
    registered". For the live, actually-registered set on an open
    model use ``model.diagnostics.selected``.
    """
    lib = get_lib()
    names = _read_catalog(lib.rdb_ocean_canonical_catalog_size,
                          lib.rdb_ocean_canonical_catalog_name)
    names += _read_catalog(lib.rdb_ocean_derived_catalog_size,
                           lib.rdb_ocean_derived_catalog_name)
    return names


# ---------------------------------------------------------------------
# Prototype declarations -- one block per include/rdb_ocean.h section.
# ---------------------------------------------------------------------

HANDLE = ctypes.c_void_p
HANDLE_OUT = ctypes.POINTER(ctypes.c_void_p)
c_int = ctypes.c_int
c_double = ctypes.c_double
c_char_p = ctypes.c_char_p
c_char = ctypes.c_char
p_int = ctypes.POINTER(c_int)
p_double = ctypes.POINTER(c_double)
c_void_p = ctypes.c_void_p
p_void_p = ctypes.POINTER(c_void_p)


def _declare_prototypes(lib):
    # ---- Lifecycle ----
    lib.rdb_ocean_create_from_string.argtypes = [c_char_p, c_int, HANDLE_OUT]
    lib.rdb_ocean_create_from_string.restype = c_int

    lib.rdb_ocean_step.argtypes = [HANDLE, c_int]
    lib.rdb_ocean_step.restype = c_int

    lib.rdb_ocean_destroy.argtypes = [HANDLE_OUT]
    lib.rdb_ocean_destroy.restype = c_int

    # ---- P2.5: pre-create geometry injection ----
    lib.rdb_ocean_create_pending.argtypes = [c_char_p, c_int, HANDLE_OUT]
    lib.rdb_ocean_create_pending.restype = c_int

    lib.rdb_ocean_stage_bathymetry.argtypes = [
        HANDLE, p_double, c_int, c_int, c_int,
    ]
    lib.rdb_ocean_stage_bathymetry.restype = c_int

    lib.rdb_ocean_stage_metrics.argtypes = [
        HANDLE, p_double, p_double, p_double, p_double, p_double,
        c_int, c_int, c_int, c_int,
    ]
    lib.rdb_ocean_stage_metrics.restype = c_int

    lib.rdb_ocean_stage_topology.argtypes = [HANDLE, c_int, c_int]
    lib.rdb_ocean_stage_topology.restype = c_int

    lib.rdb_ocean_create_finalize.argtypes = [HANDLE_OUT]
    lib.rdb_ocean_create_finalize.restype = c_int

    lib.rdb_ocean_required_halo.argtypes = [
        c_char_p, c_int, c_char_p, c_int, c_int, c_int, c_int, c_int,
    ]
    lib.rdb_ocean_required_halo.restype = c_int

    # ---- Queries ----
    lib.rdb_ocean_get_time.argtypes = [HANDLE, p_double]
    lib.rdb_ocean_get_time.restype = c_int

    lib.rdb_ocean_get_step_count.argtypes = [HANDLE, p_int]
    lib.rdb_ocean_get_step_count.restype = c_int

    lib.rdb_ocean_get_grid_info.argtypes = [HANDLE, p_int, p_int, p_int, p_int]
    lib.rdb_ocean_get_grid_info.restype = c_int

    lib.rdb_ocean_get_total_mass.argtypes = [HANDLE, p_double]
    lib.rdb_ocean_get_total_mass.restype = c_int

    lib.rdb_working_precision.argtypes = []
    lib.rdb_working_precision.restype = c_int

    # ---- Error ring ----
    lib.rdb_ocean_last_error.argtypes = [c_int, ctypes.c_char_p, c_int]
    lib.rdb_ocean_last_error.restype = c_int

    lib.rdb_flush_logs.argtypes = []
    lib.rdb_flush_logs.restype = None

    # ---- D<->H refresh ----
    lib.rdb_ocean_refresh_host.argtypes = [HANDLE]
    lib.rdb_ocean_refresh_host.restype = c_int

    # ---- Getters: (ptr, extents..., generation) ----
    _getter_3d = [HANDLE, p_void_p, p_int, p_int, p_int, p_int]
    _getter_2d = [HANDLE, p_void_p, p_int, p_int, p_int]

    for name in (
        "rdb_ocean_get_h_layer_ptr",
        "rdb_ocean_get_u_face_x_layer_ptr",
        "rdb_ocean_get_v_face_y_layer_ptr",
        "rdb_ocean_get_hu_ptr",
        "rdb_ocean_get_hv_ptr",
        "rdb_ocean_get_w_interface_ptr",
        "rdb_ocean_get_rho_layer_ptr",
        "rdb_ocean_get_kv_ptr",
        "rdb_ocean_get_kt_ptr",
        "rdb_ocean_get_ks_ptr",
    ):
        f = getattr(lib, name)
        f.argtypes = _getter_3d
        f.restype = c_int

    for name in (
        "rdb_ocean_get_b_ptr",
        "rdb_ocean_get_bt_eta_ptr",
        "rdb_ocean_get_tau_x_ptr",
        "rdb_ocean_get_tau_y_ptr",
        "rdb_ocean_get_q_heat_ptr",
        "rdb_ocean_get_q_salt_ptr",
        "rdb_ocean_get_wet_t_ptr",
    ):
        f = getattr(lib, name)
        f.argtypes = _getter_2d
        f.restype = c_int

    # ---- Tracers by name ----
    lib.rdb_ocean_get_tracer_count.argtypes = [HANDLE, p_int]
    lib.rdb_ocean_get_tracer_count.restype = c_int

    lib.rdb_ocean_list_tracers.argtypes = [HANDLE, c_int, c_char_p, c_int]
    lib.rdb_ocean_list_tracers.restype = c_int

    lib.rdb_ocean_get_tracer_ptr.argtypes = [
        HANDLE, c_char_p, c_int, p_void_p, p_int, p_int, p_int, p_int,
    ]
    lib.rdb_ocean_get_tracer_ptr.restype = c_int

    # ---- Narrow setters ----
    lib.rdb_ocean_set_h.argtypes = [HANDLE, p_double, c_int, c_int, c_int]
    lib.rdb_ocean_set_h.restype = c_int

    lib.rdb_ocean_set_u.argtypes = [HANDLE, p_double, c_int, c_int, c_int]
    lib.rdb_ocean_set_u.restype = c_int

    lib.rdb_ocean_set_v.argtypes = [HANDLE, p_double, c_int, c_int, c_int]
    lib.rdb_ocean_set_v.restype = c_int

    lib.rdb_ocean_set_bathymetry.argtypes = [HANDLE, p_double, c_int, c_int]
    lib.rdb_ocean_set_bathymetry.restype = c_int

    lib.rdb_ocean_set_wind.argtypes = [
        HANDLE, p_double, p_double, c_int, c_int,
    ]
    lib.rdb_ocean_set_wind.restype = c_int

    lib.rdb_ocean_set_heat_flux.argtypes = [HANDLE, p_double, c_int, c_int]
    lib.rdb_ocean_set_heat_flux.restype = c_int

    lib.rdb_ocean_set_salt_flux.argtypes = [HANDLE, p_double, c_int, c_int]
    lib.rdb_ocean_set_salt_flux.restype = c_int

    lib.rdb_ocean_set_tracer.argtypes = [
        HANDLE, c_char_p, c_int, p_double, c_int, c_int, c_int,
    ]
    lib.rdb_ocean_set_tracer.restype = c_int

    # ---- Scalar diagnostics ----
    lib.rdb_ocean_get_kinetic_energy.argtypes = [HANDLE, p_double]
    lib.rdb_ocean_get_kinetic_energy.restype = c_int

    # ---- P7: diagnostics discoverability. The two catalog-size/-name
    # pairs take NO handle -- they are static build-time information,
    # reachable before create(). ----
    lib.rdb_ocean_derived_catalog_size.argtypes = []
    lib.rdb_ocean_derived_catalog_size.restype = c_int

    lib.rdb_ocean_derived_catalog_name.argtypes = [c_int, c_char_p, c_int]
    lib.rdb_ocean_derived_catalog_name.restype = c_int

    lib.rdb_ocean_canonical_catalog_size.argtypes = []
    lib.rdb_ocean_canonical_catalog_size.restype = c_int

    lib.rdb_ocean_canonical_catalog_name.argtypes = [c_int, c_char_p, c_int]
    lib.rdb_ocean_canonical_catalog_name.restype = c_int

    lib.rdb_ocean_get_diag_count.argtypes = [HANDLE, p_int]
    lib.rdb_ocean_get_diag_count.restype = c_int

    lib.rdb_ocean_list_diags.argtypes = [HANDLE, c_int, c_char_p, c_int]
    lib.rdb_ocean_list_diags.restype = c_int

    # ---- P7: in-memory diagnostic access ----
    lib.rdb_ocean_get_diagnostic_ptr.argtypes = [
        HANDLE, c_char_p, c_int, p_void_p, p_int, p_int, p_int, p_int,
    ]
    lib.rdb_ocean_get_diagnostic_ptr.restype = c_int


# ---------------------------------------------------------------------
# Pure-Python array plumbing (no numpy). Fortran arrays are column-major:
# for shape (n0, n1, n2, ...), index 0 is fastest-varying.
# ---------------------------------------------------------------------


def f_strides(shape, itemsize):
    """Column-major (Fortran) byte strides for `shape`."""
    strides = []
    acc = itemsize
    for n in shape:
        strides.append(acc)
        acc *= n
    return tuple(strides)


def array_interface(address: int, shape, itemsize: int, typestr: str,
                     readonly: bool = True) -> dict:
    """Build a numpy-compatible ``__array_interface__`` dict for a raw
    address, Fortran-strided. Never imports numpy -- this is pure data."""
    strides = f_strides(shape, itemsize)
    return {
        "shape": tuple(shape),
        "typestr": typestr,
        "data": (address, readonly),
        "strides": strides,
        "version": 3,
    }


def read_element(address: int, ctype) -> float:
    """Read one scalar of `ctype` at raw address `address`."""
    return ctype.from_address(address).value


def read_flat(address: int, n: int, ctype) -> list:
    """Read `n` contiguous elements of `ctype` starting at `address` into a
    flat Python list."""
    arr_type = ctype * n
    return list(arr_type.from_address(address))


def strided_offset(index, strides):
    """Byte offset for a full integer index tuple, given Fortran strides."""
    return sum(i * s for i, s in zip(index, strides))


def nest_flat(flat, shape):
    """Reshape a flat, Fortran-ordered list into nested Python lists whose
    indexing matches numpy's: nested[i][j][k] == flat[i + j*n0 + k*n0*n1]."""
    if len(shape) == 1:
        return list(flat)
    strides = f_strides(shape, 1)
    return [_nest_axis(flat, shape, strides, (i,)) for i in range(shape[0])]


def _nest_axis(flat, shape, strides, prefix):
    depth = len(prefix)
    if depth == len(shape):
        idx = sum(p * s for p, s in zip(prefix, strides))
        return flat[idx]
    n = shape[depth]
    return [_nest_axis(flat, shape, strides, prefix + (i,)) for i in range(n)]


def flatten_nested(value, shape):
    """Flatten a nested list/tuple (numpy-style indexing: value[i][j][k])
    into a flat, Fortran-ordered list of floats matching `shape`.

    Also accepts anything exposing ``__array_interface__`` (e.g. a numpy
    array) without importing numpy: its own buffer is read directly via
    ctypes, respecting ITS strides, in the same iteration order.
    """
    if hasattr(value, "__array_interface__"):
        info = value.__array_interface__
        vshape = info["shape"]
        if tuple(vshape) != tuple(shape):
            raise ValueError(f"expected shape {shape}, got {vshape}")
        address, _ro = info["data"]
        typestr = info["typestr"]
        itemsize = int(typestr[2:])
        ctype = ctypes.c_double if itemsize == 8 else ctypes.c_float
        strides = info.get("strides")
        if strides is None:
            strides = f_strides(vshape, itemsize)

        # Iterate in Fortran order (first axis fastest) to match our own
        # flatten convention.
        def walk_fortran():
            total = 1
            for n in vshape:
                total *= n
            idx = [0] * len(vshape)
            for _ in range(total):
                off = sum(i * s for i, s in zip(idx, strides))
                flat.append(ctype.from_address(address + off).value)
                for d in range(len(idx)):
                    idx[d] += 1
                    if idx[d] < vshape[d]:
                        break
                    idx[d] = 0

        flat = []
        walk_fortran()
        return flat

    # Plain nested list/tuple: value[i][j][k]... numpy-style indexing.
    flat = [0.0] * _prod(shape)
    strides = f_strides(shape, 1)

    def rec(v, prefix):
        if len(prefix) == len(shape):
            idx = sum(p * s for p, s in zip(prefix, strides))
            flat[idx] = float(v)
            return
        if len(v) != shape[len(prefix)]:
            raise ValueError(
                f"expected extent {shape[len(prefix)]} at axis "
                f"{len(prefix)}, got {len(v)}"
            )
        for i, sub in enumerate(v):
            rec(sub, prefix + (i,))

    rec(value, ())
    return flat


def _prod(shape):
    p = 1
    for n in shape:
        p *= n
    return p


def make_double_array(flat):
    n = len(flat)
    return (c_double * n)(*flat)
