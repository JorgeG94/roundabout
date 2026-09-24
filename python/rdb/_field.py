"""D3 — field/state access: interior-by-default, lazy-sync, read-only views.

Implements the five decisions from
``tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`` D3
(and its sketch, ``sketch/03_field_semantics.py``), adapted to the
stdlib-only P3 constraint (no numpy anywhere in this module):

  1. INTERIOR BY DEFAULT. ``model.h`` is ``(nx, ny, nz)``, never the
     ghost-inclusive full extent. Ghosts are reachable only through the
     explicitly named ``.with_halo``.
  2. CONCENTRATION, NOT hTr. ``model.temperature`` is degC -- a DERIVED
     copy (the store is h*Tr), dividing only where ``h > H_VANISHED``
     and writing NaN elsewhere. The raw store is available as
     ``field.hTr``.
  3. STAGGERING IS AN ATTRIBUTE. Every field carries ``.location``, shown
     in ``repr()``.
  4. LAZY SYNC BY GENERATION, SILENT RE-SYNC. A field remembers the outer
     step count at issue and re-checks ON ACCESS -- never at the step
     boundary, and never raises for staleness. It raises only for a
     CLOSED model (the one case a stale pointer would actually dangle).
  5. WRITES ARE CALLS, NOT VIEW MUTATIONS. Raw views are returned
     read-only; ``field[...] = x`` routes to the field's narrow Fortran
     setter.
"""

from __future__ import annotations

import enum
import math

from . import _ffi
from ._errors import (RdbClosedError, RdbReadOnlyError,
                       RdbUnsupportedError, check)

# rdb_constants.F90 -- the dynamic-vanish threshold. NOT H_DIV_EPS
# (1e-20), which is pure 1/0 armour and would let a ~0 m layer produce an
# absurd "temperature".
H_VANISHED = 1.5e-4


class Loc(enum.Enum):
    CENTER = "Center"
    FACE = "Face"


class Field:
    """A live window onto one Roundabout state array.

    Holds a raw host pointer, the FULL (ghost-inclusive) Fortran extents,
    the outer step count the pointer's contents were last known good at
    ("generation"), a stagger location, and (optionally) a narrow setter.
    Never a numpy array; ``__array_interface__`` lets numpy consume it
    zero-copy when numpy happens to be present, without this module ever
    importing it.
    """

    __slots__ = ("_model", "_name", "_ptr", "_full_shape", "_gen",
                 "location", "units", "long_name", "_ghost", "_setter",
                 "_interior", "_epoch")

    def __init__(self, model, name, ptr, full_shape, generation, location,
                 units="", long_name="", setter=None, interior=True):
        self._model = model
        self._name = name
        self._ptr = ptr
        self._full_shape = tuple(full_shape)
        self._gen = generation
        self._epoch = model._epoch
        self.location = location
        self.units = units
        self.long_name = long_name
        self._ghost = model.nghost
        self._setter = setter
        self._interior = interior

    # ---- lifetime -----------------------------------------------------
    def _alive(self):
        if self._model._epoch != self._epoch:
            raise RdbClosedError(
                f"field {self._name!r} belongs to a model that has been "
                f"closed; re-create the model or take a .copy() before "
                f"leaving the `with` block")

    def _fresh(self):
        """Bring the host buffer up to date if the solver has moved on.
        Silent re-sync; never raises for staleness (D3.4) -- see module
        docstring point 4."""
        self._alive()
        # ALWAYS ask the library: `rdb_ocean_refresh_host` is gated on the
        # handle's own `host_is_current` flag (a no-op flag check when
        # nothing moved), which is the only authority on whether the HOST
        # copy is current.  Comparing the step count against this Field's
        # own `_gen` is not: a Field constructed AFTER a step is stamped
        # with the current step count by its getter, so the old
        # `gen_now != self._gen` test skipped the refresh and every such
        # Field read stale host memory on the GPU build (the initial
        # state, if nothing else had refreshed it) -- silently, on
        # `-gpu=mem:separate` only.
        check(self._model._lib.rdb_ocean_refresh_host(
            self._model._handle), "rdb_ocean_refresh_host")
        self._gen = self._model._step_count()

    # ---- shape ---------------------------------------------------------
    @property
    def shape(self):
        if not self._interior:
            return self._full_shape
        g = self._ghost
        out = []
        for i, n in enumerate(self._full_shape):
            out.append(n - 2 * g if i < 2 else n)
        return tuple(out)

    @property
    def ndim(self):
        return len(self._full_shape)

    def _elem_offset(self):
        """Element offset (not bytes) of this view's origin into the full
        buffer -- ghost*stride on the first two (horizontal) axes only;
        there is no vertical halo."""
        if not self._interior:
            return 0
        g = self._ghost
        strides = _elem_strides(self._full_shape)
        off = 0
        for i in range(min(2, len(strides))):
            off += g * strides[i]
        return off

    # ---- numpy handshake (works without numpy ever being imported) ----
    @property
    def __array_interface__(self):
        """Zero-copy, F-ordered, READ-ONLY.

        Read-only is deliberate (decision 5): a writable view would let
        ``arr[...] = x`` mutate the mapped HOST array with no device
        push, silently discarded on the next kernel launch under
        ``mem:separate``.
        """
        self._fresh()
        itemsize = self._model._wp_bytes
        typestr = self._model._wp_typestr
        strides = _ffi.f_strides(self._full_shape, itemsize)
        offset_bytes = self._elem_offset() * itemsize
        return {
            "shape": self.shape,
            "typestr": typestr,
            "data": (self._ptr + offset_bytes, True),
            "strides": strides,
            "version": 3,
        }

    @property
    def with_halo(self) -> "Field":
        """The same array including ghost cells. Named, not the default --
        a ghost-inclusive default is what forces the ``[2:-2, 2:-2]``
        idiom this API exists to retire."""
        return Field(self._model, self._name, self._ptr, self._full_shape,
                     self._gen, self.location, self.units, self.long_name,
                     setter=None, interior=False)

    # ---- reads ----------------------------------------------------------
    def _flat(self):
        """Interior (or with_halo) contents as a flat Fortran-ordered
        Python list of floats. Pure ctypes; no numpy."""
        self._fresh()
        itemsize = self._model._wp_bytes
        ctype = self._model._wp_ctype
        base_elem = self._elem_offset()
        address = self._ptr + base_elem * itemsize
        shape = self.shape
        full = self._full_shape
        if shape == full:
            return _ffi.read_flat(address, _ffi._prod(shape), ctype)
        # Strided read: walk the view's own shape against the FULL array's
        # strides (element units), since the view may be an interior crop
        # of a larger buffer. The base offset (ghost skip) must be added
        # here too -- forgetting it silently reads ghost cells instead of
        # the interior (caught by test_derived_field_differs_from_raw_htr
        # after a write: it read stale ghost data instead of the just-
        # written interior value).
        full_strides = _elem_strides(full)
        out = []
        idx = [0] * len(shape)
        total = _ffi._prod(shape)
        for _ in range(total):
            off = base_elem + sum(i * s for i, s in zip(idx, full_strides))
            out.append(_ffi.read_element(self._ptr + off * itemsize, ctype))
            for d in range(len(idx)):
                idx[d] += 1
                if idx[d] < shape[d]:
                    break
                idx[d] = 0
        return out

    def __getitem__(self, idx):
        shape = self.shape
        if idx is Ellipsis or idx == slice(None):
            return _ffi.nest_flat(self._flat(), shape)
        if isinstance(idx, tuple) and len(idx) == len(shape) and all(
                isinstance(i, int) for i in idx):
            self._fresh()
            itemsize = self._model._wp_bytes
            ctype = self._model._wp_ctype
            full_strides = _elem_strides(self._full_shape)
            base = self._elem_offset()
            norm = tuple(
                (i if i >= 0 else i + n) for i, n in zip(idx, shape))
            off = base + sum(i * s for i, s in zip(norm, full_strides))
            return _ffi.read_element(self._ptr + off * itemsize, ctype)
        # General slicing (a mix of int/slice, or fewer indices than
        # dims) has no pure-stdlib implementation here -- fall back to
        # numpy's OWN slicing over the SAME zero-copy
        # __array_interface__ handshake `np.asarray(field)` already
        # gives it (D6.6: "install numpy and use np.asarray(field) for
        # general slicing" -- this just does that automatically when
        # numpy happens to be importable, rather than making every
        # caller spell it out).
        try:
            import numpy as _np
        except ImportError:
            raise RdbUnsupportedError(
                f"{self._name}: partial slicing without numpy is not "
                f"supported in this build (P3 is stdlib-only). Use "
                f"{self._name}[...] or {self._name}.copy() for the full "
                f"array, a full integer index for a scalar, or install "
                f"numpy and use np.asarray({self._name}) for general "
                f"slicing.") from None
        return _np.asarray(self)[idx]

    def copy(self):
        """An owned snapshot that does NOT track the solver -- the
        documented way to keep data across steps."""
        return _ffi.nest_flat(self._flat(), self.shape)

    def nodes(self):
        """Coordinates at this field's stagger point (deferred: no grid
        metrics accessor is exposed by the C ABI yet, so this always
        raises in P3)."""
        raise RdbUnsupportedError(
            f"{self._name}.nodes(): grid coordinate accessors are not "
            f"part of the P3 C ABI surface (rdb_ocean.h has no "
            f"geolon/geolat getter yet)."
        )

    # ---- writes (decision 5) -------------------------------------------
    def __setitem__(self, idx, value):
        self._alive()
        if self._setter is None:
            raise RdbReadOnlyError(
                f"{self._name!r} is read-only. Writable fields: "
                f"{', '.join(sorted(self._model.writable_fields))}. "
                f"(ssh/eta are diagnostic -- sum_k h - b -- set `h` or "
                f"`bathymetry` instead.)")
        full = (idx is Ellipsis) or (idx == slice(None))
        shape = self.shape
        if full:
            flat = _ffi.flatten_nested(value, shape)
        else:
            nested = self.copy()
            _assign_into(nested, idx, value)
            flat = _ffi.flatten_nested(nested, shape)
        arr = _ffi.make_double_array(flat)
        status = self._setter(self._model._handle, arr, *shape)
        check(status, f"{self._name} setter")
        # The setter left the DEVICE canonical; poison the generation so
        # the next read refreshes rather than pulling eagerly now.
        self._gen = -1

    def __repr__(self):
        locs = ", ".join(loc.value for loc in self.location)
        scope = "interior" if self._interior else "with_halo"
        return (f"Field({self._name}, ({locs}), shape={self.shape}, "
                f"unit={self.units!r}, {scope}, generation={self._gen})")


class DerivedField(Field):
    """A field COMPUTED from the store (hTr / h), so it is copy-only.

    THE GUARD: concentration = hTr / h only where ``h > H_VANISHED``
    (1.5e-4 m); elsewhere NaN -- not 0.0 (0 degC / 0 PSU are both legal
    ocean values, so 0.0 would read a vanished layer as ice-point
    freshwater rather than as missing; see D3.2).
    """

    __slots__ = ("_h_field",)

    def __init__(self, model, name, htr_field, h_field, units="",
                 long_name="", setter=None):
        self._model = model
        self._name = name
        self._ptr = htr_field._ptr
        self._full_shape = htr_field._full_shape
        self._gen = htr_field._gen
        self._epoch = model._epoch
        self.location = htr_field.location
        self.units = units
        self.long_name = long_name
        self._ghost = model.nghost
        self._setter = setter
        self._interior = True
        self._h_field = h_field

    @property
    def hTr(self):
        """The raw h*Tr store -- a zero-copy Field view."""
        return Field(self._model, self._name + ".hTr", self._ptr,
                     self._full_shape, self._gen, self.location,
                     units=f"m*({self.units})" if self.units else "",
                     long_name=self.long_name, setter=None, interior=True)

    @property
    def __array_interface__(self):
        raise RdbUnsupportedError(
            f"{self._name} is derived (hTr / h) and has no contiguous "
            f"buffer of its own. Use {self._name}.hTr for the raw store, "
            f"or {self._name}.copy() / {self._name}[...] for the "
            f"concentration.")

    def copy(self):
        htr = Field(self._model, self._name, self._ptr, self._full_shape,
                    self._gen, self.location, interior=True)
        htr_nested = htr.copy()
        h_nested = self._h_field.copy()
        return _divide_guarded(htr_nested, h_nested)

    def __getitem__(self, idx):
        nested = self.copy()
        if idx is Ellipsis or idx == slice(None):
            return nested
        if isinstance(idx, tuple) and len(idx) == self.ndim and all(
                isinstance(i, int) for i in idx):
            v = nested
            for i in idx:
                v = v[i]
            return v
        # Same numpy fallback as Field.__getitem__, but over the COPY
        # (`nested`) -- a DerivedField has no contiguous buffer of its
        # own to hand numpy zero-copy (see __array_interface__ above),
        # so this is np.asarray(a Python nested list), not a view.
        try:
            import numpy as _np
        except ImportError:
            raise RdbUnsupportedError(
                f"{self._name}: partial slicing without numpy is not "
                f"supported in this build. Use {self._name}[...] or "
                f"{self._name}.copy(), or install numpy.") from None
        return _np.asarray(nested)[idx]


def _divide_guarded(htr_nested, h_nested):
    if isinstance(htr_nested, list):
        return [_divide_guarded(a, b) for a, b in zip(htr_nested, h_nested)]
    return htr_nested / h_nested if h_nested > H_VANISHED else math.nan


def _assign_into(nested, idx, value):
    """Best-effort read-modify-write into a nested-list snapshot for a
    small set of common index shapes (full integer tuple, or a single
    leading-axis index/slice). Anything fancier: use numpy."""
    if isinstance(idx, tuple) and all(isinstance(i, int) for i in idx):
        v = nested
        for i in idx[:-1]:
            v = v[i]
        v[idx[-1]] = value
        return
    if isinstance(idx, int):
        nested[idx] = value
        return
    raise RdbUnsupportedError(
        "partial-slice writes beyond a full integer index are not "
        "supported without numpy in this build; write the whole field "
        "(field[...] = x) or install numpy."
    )


class DiagnosticField(Field):
    """A read-only view onto the diag manager's OWN ``output_buffer`` for
    one registered diagnostic (P7) — no NetCDF round-trip.

    Unlike :class:`Field`, the generation is the diagnostic's own
    per-fire counter, not the outer step count, and needs no separate
    refresh call: the Fortran side (`ocean_diag_step`) pulls
    ``output_buffer`` host-ward UNCONDITIONALLY at every cadence fire,
    synchronously inside ``rdb_ocean_step`` — so by the time
    :meth:`Model.step` returns, the data is already host-current for
    anything that fired. ``_fresh`` is overridden to skip the
    generation-vs-step-count comparison (meaningless here — a
    diagnostic on a multi-hour cadence should not look "stale" between
    fires) and the ``rdb_ocean_refresh_host`` call (which does not
    even touch ``output_buffer``).

    Read-only: there is no Fortran setter for a diagnostic output
    buffer (`register_default_diags`/derived fills own the write path),
    so ``field[...] = x`` raises :class:`~rdb._errors.RdbReadOnlyError`
    via the inherited ``__setitem__`` (``setter=None``).
    """

    def _fresh(self):
        self._alive()


def _elem_strides(shape):
    return _ffi.f_strides(shape, 1)
