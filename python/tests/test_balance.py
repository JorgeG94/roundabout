"""Balanced velocity from the model's own PGF: ``Model.balance`` /
``Model.stage_balance``.

The physical test (B3) seeds an elliptical Gaussian warm eddy through the
staged tracer path on an f-plane and compares a balanced start with a twin
started from rest.  The balanced eddy must stay near its initial kinetic
energy and sea surface; the rest twin must visibly adjust (otherwise the
bands prove nothing).  The residency check (B6) is the live-vs-staged
equality: on the GPU build ``balance()`` runs on device-resident state, so a
missed device sync shows up as a mismatch against the staged path, which
balances at finalize.
"""

import math

import pytest

import rdb

NX, NY, NZ = 24, 24, 4
DX = 10000.0
DEPTH = 1000.0
F0 = 1.0e-4
DT = 900.0
T_BG = 10.0


def _nml(dt=DT, t_end=None):
    t_end = t_end if t_end is not None else 100 * dt
    return (
        '&sim_nml sim_type = "ocean" /\n'
        f"&time_nml t_end = {t_end}, dt_fixed = {dt} /\n"
        f"&nonhydrostatic_nml nz_layers = {NZ} /\n"
        f"&tracer_nml initial_temperature = {T_BG}, initial_salinity = 35.0 /\n"
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&ocean_diag_nml enabled = .false. /\n"
        f"&grid_nml nx = {NX}, ny = {NY}, nghost = 3, dx = {DX}, dy = {DX} /\n"
        f"&physics_nml coriolis_f = {F0} /\n"
        "&vcoord_nml vcoord_type = 'sigma' /\n"
        f"&ocean_topo_nml topo_config = 'flat', max_depth = {DEPTH}, "
        "coriolis_beta = 0.0, coriolis_y_ref = 0.0 /\n"
        "&ocean_bc_nml west = 'periodic', east = 'periodic', "
        "south = 'periodic', north = 'periodic' /\n"
        "&ocean_eos_nml eos = 'linear' /\n"
        "&ocean_thermo_nml enable_thermodynamics = .true. /\n"
        "&ocean_ic_nml alpha_T = 0.17, rho_0 = 1035.0 /\n"
        "&output_nml output_to_file = .false. /\n"
    )


def _eddy(tilt=0.0, amp=1.0, ax=5.0, ay=3.0):
    """Elliptical Gaussian warm anomaly, surface-intensified (k = 0 bed);
    ``tilt`` shifts the centre by that many cells per layer."""
    cx, cy = NX / 2.0, NY / 2.0

    def fn(i, j, k):
        x = i + 0.5 - cx - tilt * (k - (NZ - 1) / 2.0)
        y = j + 0.5 - cy
        return T_BG + amp * (k + 1) / NZ * math.exp(
            -(x / ax) ** 2 - (y / ay) ** 2)
    return [[[fn(i, j, k) for k in range(NZ)] for j in range(NY)]
            for i in range(NX)]


def _staged(temp, balance=None, nml=None, **bal_kw):
    m = rdb.Model.pending(nml or _nml())
    m.stage_tracer("temperature", temp)
    if balance:
        m.stage_balance(balance, **bal_kw)
    m.finalize()
    return m


def _eta_stats(m, ref=None):
    e = m.eta[...]
    flat = [e[i][j] for i in range(NX) for j in range(NY)]
    if ref is None:
        return max(abs(v) for v in flat), flat
    return max(abs(a - b) for a, b in zip(flat, ref)), flat


def _run(m, n):
    ke0 = m.kinetic_energy
    emax0, eta0 = _eta_stats(m)
    m.step(n)
    ke1 = m.kinetic_energy
    de, _ = _eta_stats(m, eta0)
    return ke0, ke1, emax0, de


@pytest.mark.parametrize("tilt", [0.0, 1.5])
def test_balanced_eddy_is_nearly_steady_and_rest_twin_is_not(tilt):
    temp = _eddy(tilt=tilt)
    n = 70                     # ~ one inertial period (2 pi / f = 62.8 ks)
    b = _staged(temp, "geostrophic")
    try:
        bke0, bke1, bemax, bde = _run(b, n)
    finally:
        b.close()
    r = _staged(temp)
    try:
        rke0, rke1, remax, rde = _run(r, n)
    finally:
        r.close()
    print(f"\ntilt={tilt}: balanced KE {bke0:.4e}->{bke1:.4e} "
          f"(ratio {bke1 / bke0:.3f}) eta_max {bemax:.4e} d_eta/eta_max "
          f"{bde / bemax:.3f}; rest KE {rke0:.4e}->{rke1:.4e} "
          f"(of balanced {rke1 / bke0:.3f}) d_eta/balanced eta_max "
          f"{rde / bemax:.3f}")
    # Measured (24x24x4, dx 10 km, f 1e-4, dt 900 s, 70 steps): balanced KE
    # ratio 1.025 / 1.029, balanced d_eta/eta_max 0.043 / 0.077; rest twin
    # d_eta/balanced eta_max 0.50 / 0.53 and KE still 2-3 % of the balanced one.
    assert rke0 == 0.0 and bke0 > 0.0
    assert 0.9 < bke1 / bke0 < 1.1
    assert bde < 0.15 * bemax
    assert rde > 0.3 * bemax and rde > 4.0 * bde
    assert rke1 < 0.1 * bke0


def _uvh(m):
    return m.u[...], m.v[...], m.h[...], m.eta[...]


def test_live_balance_equals_staged_balance_bit_for_bit():
    """B6 / residency: the live ``balance()`` (device-resident state on the GPU
    build) must reproduce the staged path, which balances at finalize --
    including after stepping, so a stale device copy cannot hide."""
    temp = _eddy(tilt=1.5)
    a = _staged(temp, "geostrophic")
    try:
        ua = _uvh(a)
        a.step(3)
        ua3 = _uvh(a)
    finally:
        a.close()
    b = rdb.Model(_nml())
    try:
        b.temperature[...] = temp
        b.balance("geostrophic", "bottom_rest")
        ub = _uvh(b)
        b.step(3)
        ub3 = _uvh(b)
    finally:
        b.close()
    assert max(abs(x) for plane in ua[0] for row in plane for x in row) > 0.0
    assert ua == ub
    assert ua3 == ub3


def test_keep_h_leaves_the_column_thickness_alone():
    temp = _eddy()
    m = rdb.Model(_nml())
    try:
        m.temperature[...] = temp
        h0 = m.h[...]
        m.balance("geostrophic", "keep_h")
        assert m.h[...] == h0
        assert max(abs(x) for plane in m.u[...] for row in plane
                   for x in row) > 0.0
    finally:
        m.close()


def test_bottom_rest_relevels_the_sea_surface():
    m = rdb.Model(_nml())
    try:
        m.temperature[...] = _eddy()
        m.balance()
        assert max(abs(e) for row in m.eta[...] for e in row) > 1.0e-2
    finally:
        m.close()


def test_gradient_mode_is_not_implemented():
    m = rdb.Model(_nml())
    try:
        with pytest.raises(rdb.RdbInitialConditionError):
            m.balance("gradient")
    finally:
        m.close()
    p = rdb.Model.pending(_nml())
    try:
        with pytest.raises(rdb.RdbInitialConditionError):
            p.stage_balance("gradient")
        p.finalize()          # the refused request staged nothing
    finally:
        p.close()


def test_balance_refused_after_the_first_step():
    with rdb.Model(_nml()) as m:
        m.step(1)
        with pytest.raises(rdb.RdbSetupError):
            m.balance()


def test_bad_arguments_and_wrong_window():
    with rdb.Model(_nml()) as m:
        with pytest.raises(ValueError):
            m.balance("geostrophic", "nonsense")
        with pytest.raises(ValueError):
            m.balance("nonsense")
        with pytest.raises(rdb.NotInitialisedError):
            m.stage_balance()
    p = rdb.Model.pending(_nml())
    try:
        with pytest.raises(ValueError):
            p.stage_balance("geostrophic", "nonsense")
    finally:
        p.close()


def test_equatorial_f_is_refused_and_velocity_untouched():
    nml = _nml().replace(f"coriolis_f = {F0}", "coriolis_f = 0.0")
    with rdb.Model(nml) as m:
        m.temperature[...] = _eddy()
        u0 = m.u[...]
        with pytest.raises(rdb.RdbInitialConditionError):
            m.balance()
        assert m.u[...] == u0


def test_balance_overrides_staged_velocity():
    p = rdb.Model.pending(_nml())
    p.stage_tracer("temperature", _eddy())
    p.stage_u([[[0.5] * NZ for _ in range(NY)] for _ in range(NX + 1)])
    p.stage_balance()
    p.finalize()
    try:
        peak = max(abs(x) for plane in p.u[...] for row in plane for x in row)
        assert 0.0 < peak < 0.5
    finally:
        p.close()


def _nml_no_thermo(layer_rho=None):
    nml = _nml().replace("enable_thermodynamics = .true.",
                         "enable_thermodynamics = .false.")
    if layer_rho:
        nml = nml.replace(
            "&ocean_ic_nml alpha_T = 0.17, rho_0 = 1035.0 /",
            "&ocean_ic_nml alpha_T = 0.17, rho_0 = 1035.0, layer_rho_init = "
            + ", ".join(str(v) for v in layer_rho) + " /")
    return nml


def test_thermo_off_refuses_balance_of_written_temperature():
    """T/S cannot change rho with thermodynamics off: refuse, don't return a
    silent zero velocity."""
    with rdb.Model(_nml_no_thermo()) as m:
        m.temperature[...] = _eddy()
        with pytest.raises(rdb.RdbInitialConditionError):
            m.balance()
    p = rdb.Model.pending(_nml_no_thermo())
    try:
        p.stage_tracer("temperature", _eddy())
        with pytest.raises(rdb.RdbInitialConditionError):
            p.stage_balance()
        p.finalize()
    finally:
        p.close()
    p = rdb.Model.pending(_nml_no_thermo())
    try:
        p.stage_balance()
        with pytest.raises(rdb.RdbInitialConditionError):
            p.stage_tracer("salinity", _eddy())
        p.finalize()
    finally:
        p.close()


def test_thermo_off_with_prescribed_layer_densities_is_allowed():
    rho = [1034.0, 1034.3, 1034.6, 1034.9]
    with rdb.Model(_nml_no_thermo(rho)) as m:
        m.balance()
    p = rdb.Model.pending(_nml_no_thermo(rho))
    try:
        p.stage_balance()
        p.finalize()
    finally:
        p.close()
