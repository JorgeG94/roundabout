"""Array setters: numpy-protocol inputs (no numpy import) and layout.

Regression: per the ``__array_interface__`` protocol ``strides is None``
means C-contiguous.  ``flatten_nested`` once treated it as Fortran order, so
a default (C-order) numpy array was silently scrambled on write.
"""

import array
import itertools

import pytest

import rdb
from rdb import _ffi


class FakeArray:
    """Minimal ``__array_interface__`` object over an ``array.array``."""

    def __init__(self, shape, order="C", code="d", typestr="<f8", fn=None):
        n = 1
        for s in shape:
            n *= s
        isz = array.array(code).itemsize
        strides = (_ffi.c_strides(shape, isz) if order == "C"
                   else _ffi.f_strides(shape, isz))
        buf = array.array(code, [0] * n)
        for idx in itertools.product(*[range(s) for s in shape]):
            buf[sum(i * s for i, s in zip(idx, strides)) // isz] = fn(idx)
        self._buf = buf
        self._shape = tuple(shape)
        self._strides = strides if order == "F" else None  # C => None
        self._typestr = typestr

    @property
    def __array_interface__(self):
        return {"shape": self._shape, "typestr": self._typestr,
                "data": (self._buf.buffer_info()[0], False),
                "strides": self._strides, "version": 3}


def _val(idx):
    return 1.0 + sum((i + 1) * 10.0 ** k for k, i in enumerate(idx[::-1]))


def _check(field, shape):
    got = field[...]
    for idx in itertools.product(*[range(s) for s in shape]):
        g = got
        for i in idx:
            g = g[i]
        assert g == pytest.approx(_val(idx), rel=1e-12), idx  # T is hTr/h


@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("name", ["h", "temperature"])
def test_3d_setter_layout(ocean_nml, name, order):
    with rdb.Model(ocean_nml) as m:
        f = getattr(m, name)
        arr = FakeArray(f.shape, order, fn=_val)
        f[...] = arr
        _check(getattr(m, name), f.shape)


@pytest.mark.parametrize("order", ["C", "F"])
def test_2d_setter_layout(ocean_nml, order):
    with rdb.Model(ocean_nml) as m:
        f = m.bathymetry
        arr = FakeArray(f.shape, order, fn=_val)
        f[...] = arr
        _check(m.bathymetry, f.shape)


def test_float32_c_order(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        shape = m.h.shape
        arr = FakeArray(shape, "C", code="f", typestr="<f4", fn=_val)
        m.h[...] = arr
        _check(m.h, shape)


def test_nested_list_still_works(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        nx, ny, nz = m.h.shape
        m.h[...] = [[[_val((i, j, k)) for k in range(nz)] for j in range(ny)]
                    for i in range(nx)]
        _check(m.h, (nx, ny, nz))


def test_f_contiguous_is_zero_copy():
    shape = (4, 3, 2)
    arr = FakeArray(shape, "F", fn=_val)
    ptr = _ffi.as_double_buffer(arr, shape)
    addr = ctypes_addr(ptr)
    assert addr == arr._buf.buffer_info()[0]
    assert ptr._keepalive is arr


def ctypes_addr(ptr):
    import ctypes
    return ctypes.cast(ptr, ctypes.c_void_p).value


def test_c_order_flatten_is_fortran_ordered():
    shape = (3, 4)
    arr = FakeArray(shape, "C", fn=_val)
    flat = _ffi.flatten_nested(arr, shape)
    for i in range(3):
        for j in range(4):
            assert flat[i + 3 * j] == _val((i, j))


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        _ffi.as_double_buffer(FakeArray((2, 3), "C", fn=_val), (3, 2))


class BroadcastArray:
    """``__array_interface__`` view with a ZERO stride on one axis, as
    ``np.broadcast_to`` produces: every index along `axis` reads the same
    element of a smaller C-ordered buffer."""

    def __init__(self, shape, axis, fn):
        base_shape = tuple(1 if d == axis else s for d, s in enumerate(shape))
        n = 1
        for s in base_shape:
            n *= s
        buf = array.array("d", [0.0] * n)
        cs = _ffi.c_strides(base_shape, 8)
        for idx in itertools.product(*[range(s) for s in base_shape]):
            buf[sum(i * s for i, s in zip(idx, cs)) // 8] = fn(idx)
        self._buf = buf
        self._shape = tuple(shape)
        self._strides = tuple(0 if d == axis else s for d, s in enumerate(cs))

    @property
    def __array_interface__(self):
        return {"shape": self._shape, "typestr": "<f8",
                "data": (self._buf.buffer_info()[0], False),
                "strides": self._strides, "version": 3}


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_zero_stride_broadcast_flattens(axis):
    shape = (4, 3, 5)
    fn = lambda idx: _val(tuple(0 if d == axis else i for d, i in enumerate(idx)))
    flat = _ffi.flatten_nested(BroadcastArray(shape, axis, fn), shape)
    fs = _ffi.f_strides(shape, 1)
    for idx in itertools.product(*[range(s) for s in shape]):
        assert flat[sum(i * s for i, s in zip(idx, fs))] == fn(idx), idx


def test_zero_stride_broadcast_through_setter(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        shape = m.h.shape
        fn = lambda idx: _val((0,) + tuple(idx[1:]))
        m.h[...] = BroadcastArray(shape, 0, fn)
        got = m.h[...]
        for idx in itertools.product(*[range(s) for s in shape]):
            g = got
            for i in idx:
                g = g[i]
            assert g == fn(idx), idx
