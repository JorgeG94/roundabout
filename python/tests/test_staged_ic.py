"""Staged initial state (stage_h / stage_u / stage_v / stage_tracer) in the
pending-create window: the data lands after the built-in seed and before every
setup-time snapshot.  The sponge-reference and pseudo-salt cases need engine
internals and live in tests/test_ocean_staged_ic.F90; this file covers what is
observable through the public surface."""

import itertools

import pytest

import rdb

NX, NY, NZ = 8, 6, 3
EPS = 2.220446049250313e-16


def _nml(wind=False, depth=200.0):
    return (
        '&sim_nml sim_type = "ocean" /\n'
        f"&grid_nml nx = {NX}, ny = {NY}, dx = 2000.0, dy = 2000.0 /\n"
        f"&nonhydrostatic_nml nz_layers = {NZ} /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        f'&ocean_topo_nml topo_config = "flat", max_depth = {depth} /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        + ("&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.05 /\n" if wind else "")
        + "&ocean_diag_nml enabled = .false. /\n"
        "&output_nml output_to_file = .false. /\n"
    )


def _field(nx, ny, nz, fn):
    return [[[fn(i, j, k) for k in range(nz)] for j in range(ny)]
            for i in range(nx)]


def _h(depth=200.0):
    # weights 1:2:3 bottom-up, columns sum to `depth`
    return _field(NX, NY, NZ, lambda i, j, k: depth * (k + 1) / 6.0)


def _temp():
    return _field(NX, NY, NZ,
                  lambda i, j, k: 8.0 + 0.1 * i + 0.03 * j + 1.7 * k)


def test_staged_h_then_T_gives_staged_concentration():
    T = _temp()
    h = _h()
    h[2][3] = [h[2][3][0] * 1.2, h[2][3][1] * 0.9, h[2][3][2] * 0.95]
    s = sum(h[2][3])
    h[2][3] = [x * 200.0 / s for x in h[2][3]]
    m = rdb.Model.pending(_nml())
    try:
        m.stage_tracer("temperature", T)    # T before h: order must not matter
        m.stage_h(h)
        m.finalize()
        got = m.temperature[...]
        hh = m.h[...]
        for i, j, k in itertools.product(range(NX), range(NY), range(NZ)):
            assert abs(got[i][j][k] - T[i][j][k]) <= 4 * EPS * abs(T[i][j][k])
            assert hh[i][j][k] == h[i][j][k]
        # unstaged salinity keeps its concentration across the staged h
        sal = m.salinity[...]
        assert abs(sal[2][3][1] - 35.0) <= 8 * EPS * 35.0
    finally:
        m.close()


def test_staged_u_zeroed_on_land_faces():
    depth = [[200.0] * NY for _ in range(NX)]
    depth[3][2] = 0.0                       # one land cell
    u = _field(NX + 1, NY, NZ, lambda i, j, k: 0.5)
    v = _field(NX, NY + 1, NZ, lambda i, j, k: 0.25)
    m = rdb.Model.pending(_nml())
    try:
        m.stage_bathymetry(depth, rdb.BATHY_DEPTH_POSITIVE_DOWN)
        m.stage_u(u)
        m.stage_v(v)
        m.finalize()
        uu = m.u[...]
        vv = m.v[...]
        for k in range(NZ):
            assert uu[3][2][k] == 0.0 and uu[4][2][k] == 0.0
            assert vv[3][2][k] == 0.0 and vv[3][3][k] == 0.0
        assert uu[5][4][0] == 0.5           # open interior face is untouched
        assert vv[5][4][0] == 0.25
    finally:
        m.close()


def test_pending_without_staging_is_bit_identical_to_one_shot():
    nml = _nml(wind=True)
    with rdb.Model(nml) as a:
        a.step(20)
        ref = (a.h[...], a.u[...], a.v[...], a.temperature[...],
               a.salinity[...])
    m = rdb.Model.pending(nml)
    try:
        m.finalize()
        m.step(20)
        got = (m.h[...], m.u[...], m.v[...], m.temperature[...],
               m.salinity[...])
    finally:
        m.close()
    assert got == ref


def test_staged_h_to_uniform_equals_default_seed():
    """Staging exactly what the seed would have produced changes nothing."""
    nml = _nml(wind=True)
    with rdb.Model(nml) as a:
        a.step(5)
        ref = (a.h[...], a.u[...], a.temperature[...])
    m = rdb.Model.pending(nml)
    try:
        m.stage_h(_field(NX, NY, NZ, lambda i, j, k: 200.0 / NZ))
        m.stage_tracer("temperature", _field(NX, NY, NZ, lambda i, j, k: 12.0))
        m.finalize()
        m.step(5)
        got = (m.h[...], m.u[...], m.temperature[...])
    finally:
        m.close()
    for r, g in zip(ref, got):
        for i, j, k in itertools.product(range(len(r)), range(len(r[0])),
                                         range(len(r[0][0]))):
            assert abs(g[i][j][k] - r[i][j][k]) <= 1e-12 * (1 + abs(r[i][j][k]))


def test_wrong_shape_raises_bad_shape():
    m = rdb.Model.pending(_nml())
    try:
        m.stage_u(_field(NX, NY, NZ, lambda i, j, k: 0.0))   # needs nx+1
        with pytest.raises(rdb.BadShapeError):
            m.finalize()
    finally:
        m.close()


def test_unknown_tracer_raises_not_found():
    m = rdb.Model.pending(_nml())
    try:
        m.stage_tracer("no_such_tracer", _field(NX, NY, NZ, lambda i, j, k: 1.0))
        with pytest.raises(rdb.TracerNotFoundError):
            m.finalize()
    finally:
        m.close()


def test_nonpositive_h_raises_ic_error():
    h = _h()
    h[1][1][1] = 0.0
    m = rdb.Model.pending(_nml())
    try:
        m.stage_h(h)
        with pytest.raises(rdb.RdbInitialConditionError):
            m.finalize()
    finally:
        m.close()


def test_stage_on_non_pending_model_raises(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        with pytest.raises(rdb.NotInitialisedError):
            m.stage_h(_h())
