"""D3: interior-by-default shape, read-only views, tracer concentration
fields, and the lazy-sync generation counter -- silent re-sync on access,
no explicit sync_host() call ever required."""

import rdb


def test_field_read_returns_interior_shape_not_ghost_shape(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        info = m.grid_info
        h = m.h
        assert h.shape == (info["nx"], info["ny"], info["nz"])
        # The with_halo view is strictly larger on the horizontal axes.
        full = h.with_halo
        g = info["nghost"]
        assert full.shape[0] == info["nx"] + 2 * g
        assert full.shape[1] == info["ny"] + 2 * g


def test_view_is_read_only(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        h = m.h
        info = h.__array_interface__
        assert info["data"][1] is True  # (address, readonly) -- readonly


def test_tracer_field_is_the_stored_concentration(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        T = m.temperature
        t_val = T[0, 0, 0]
        # temperature is the namelist's initial_temperature=12.0 (the store
        # IS the concentration, not h*T, so h plays no part).
        assert abs(t_val - 12.0) < 1e-6
        # a zero-copy raw field: it exposes a buffer
        assert T.__array_interface__["data"][1] is True  # read-only view


def test_field_read_after_step_silently_resyncs(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        h = m.h  # Field object taken BEFORE stepping
        before = h[0, 0, 0]
        m.step(5)
        # No explicit sync_host()/refresh call -- just read again.
        after = h[0, 0, 0]
        assert h._gen == m.step_count  # generation silently caught up
        # (values may or may not have changed physically; the point is
        # the read succeeded and re-synced without raising or hanging.)
        assert isinstance(before, float) and isinstance(after, float)


def test_writable_field_write_then_read_back(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        h = m.h
        nx, ny, nz = h.shape
        new_h = [[[42.0 for _ in range(nz)] for _ in range(ny)]
                 for _ in range(nx)]
        h[...] = new_h
        assert m.h[0, 0, 0] == 42.0
        assert m.h[nx - 1, ny - 1, nz - 1] == 42.0


def test_full_copy_reads_interior_not_ghosts_after_write(ocean_nml):
    """Regression: Field._flat()'s strided-read branch (taken whenever
    nghost > 0, i.e. shape != with_halo shape) once forgot to add the
    interior offset, so .copy()/[...] silently read GHOST cells -- still
    holding the pre-write value -- instead of the just-written interior.
    A single-element read (m.h[0, 0, 0]) took a different code path and
    was unaffected, which is exactly why this needs its own full-array
    assertion rather than trusting the scalar-index tests above."""
    with rdb.Model(ocean_nml) as m:
        assert m.grid_info["nghost"] > 0  # otherwise this test proves nothing
        T = m.temperature
        nx, ny, nz = T.shape
        new_t = [[[20.0 for _ in range(nz)] for _ in range(ny)]
                 for _ in range(nx)]
        T[...] = new_t
        scalar = m.temperature[0, 0, 0]
        full = m.temperature[...]
        assert scalar == 20.0
        assert full[0][0][0] == 20.0
        assert full[0][0][0] == scalar


def test_readonly_field_setitem_raises(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        rho = m.rho
        try:
            rho[...] = rho.copy()
            assert False, "expected RdbReadOnlyError"
        except rdb.RdbReadOnlyError:
            pass


def test_closed_model_field_access_raises(ocean_nml):
    m = rdb.Model(ocean_nml)
    h = m.h
    m.close()
    try:
        _ = h[0, 0, 0]
        assert False, "expected RdbClosedError"
    except rdb.RdbClosedError:
        pass


def test_field_taken_after_step_reads_current_state(ocean_nml):
    """Regression (GPU `mem:separate` only): a Field constructed AFTER
    stepping is stamped by its getter with the current step count, so a
    `_fresh()` that compared only against the Field's own stamp never
    refreshed the host copy and read the INITIAL state. A constant wind
    spins up the surface layer from rest; a Field taken only after the
    steps must see that motion."""
    nml = ocean_nml + "&physics_nml wind_stress_x = 0.1 /\n"
    with rdb.Model(nml) as m:
        m.step(3)
        u = m.u.copy()  # first Field touch happens after the steps
        umax = max(abs(v) for plane in u for col in plane for v in col)
        assert umax > 0.0
