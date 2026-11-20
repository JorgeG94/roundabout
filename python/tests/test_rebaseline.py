"""``Model.rebaseline()``: re-deriving the setup-time snapshots after
post-create sets.  The sponge reference and the pseudo-salt seed are not
readable through the public surface, so the discriminating check is the
10-step bit-identity against the staged build (the sponge relaxes toward the
reference, so a wrong or un-pushed reference changes T/S within a few steps);
the direct sponge/pseudo-salt comparison lives in
tests/test_ocean_rebaseline.F90.  On the GPU build this is the device check:
every compared field is read back through the API getters."""

import pytest

import rdb

NX, NY, NZ = 10, 8, 3
DEPTH = 300.0
LAND = (3, 2)           # 0-based island cell


def _toml(sponge=True, closed=False):
    bc = ("[ocean.bc]\nwest = \"wall\"\neast = \"wall\"\nsouth = \"wall\"\n"
          "north = \"wall\"\n" if closed else
          "[ocean.bc]\nwest = \"periodic\"\neast = \"periodic\"\n"
          "south = \"wall\"\nnorth = \"sponge\"\nsponge_width = 3\n"
          "sponge_strength = 2.0e-3\nsponge_relax_tracers = true\n")
    sp = ("[ocean.sponge]\nenable = true\ndamp_source = \"band\"\n"
          "target_source = \"ic\"\nrelax_uv = true\nrelax_tracers = true\n"
          if sponge and not closed else "")
    return (
        f"""[sim]
sim_type = "ocean"
[time]
t_end = 3600.0
dt_fixed = 300.0
[nonhydrostatic]
nz_layers = {NZ}
[tracer]
initial_temperature = 12.0
initial_salinity = 35.0
[ocean.bt]
auto_n_inner = true
[ocean.diag]
enabled = false
[grid]
nx = {NX}
ny = {NY}
nghost = 3
dx = 20000.0
dy = 20000.0
[physics]
coriolis_f = 1.0e-4
[vcoord]
vcoord_type = "sigma"
[ocean.topo]
topo_config = "flat"
max_depth = {DEPTH}
"""
        + bc + sp +
        """[ocean.tracers]
enable_pseudo_salt = true
[output]
output_to_file = false
"""
    )


def _field(nx, ny, nz, fn):
    return [[[fn(i, j, k) for k in range(nz)] for j in range(ny)]
            for i in range(nx)]


def _depth():
    d = [[DEPTH] * NY for _ in range(NX)]
    d[LAND[0]][LAND[1]] = 0.0
    return d


def _ic():
    h = _field(NX, NY, NZ, lambda i, j, k: DEPTH * (k + 1) / 6.0
               * (1.0 + 1.0e-3 * (i + 1)))
    return dict(
        h=h,
        T=_field(NX, NY, NZ, lambda i, j, k: 8.0 + 0.1 * i + 0.03 * j + 1.7 * k),
        S=_field(NX, NY, NZ, lambda i, j, k: 33.0 + 0.05 * i + 0.01 * j + 0.4 * k),
        u=_field(NX + 1, NY, NZ, lambda i, j, k: 0.01 + 0.01 * i - 0.002 * j + 0.01 * k),
        v=_field(NX, NY + 1, NZ, lambda i, j, k: 0.004 * j - 0.001 * i + 0.02 * k),
    )


def _staged(cfg_text, ic):
    m = rdb.Model.pending(cfg_text)
    m.stage_bathymetry(_depth(), rdb.BATHY_DEPTH_POSITIVE_DOWN)
    m.stage_h(ic["h"])
    m.stage_tracer("temperature", ic["T"])
    m.stage_tracer("salinity", ic["S"])
    m.stage_u(ic["u"])
    m.stage_v(ic["v"])
    m.finalize()
    return m


def _post_create(cfg_text, ic):
    m = rdb.Model.pending(cfg_text)
    m.stage_bathymetry(_depth(), rdb.BATHY_DEPTH_POSITIVE_DOWN)
    m.finalize()
    m.h[...] = ic["h"]                 # h first: the tracer setters read it
    m.temperature[...] = ic["T"]
    m.salinity[...] = ic["S"]
    m.u[...] = ic["u"]
    m.v[...] = ic["v"]
    return m


def _read(m):
    out = {name: getattr(m, name)[...] for name in ("h", "u", "v")}
    out["T"] = m.temperature[...]
    out["S"] = m.salinity[...]
    out["P"] = m.tracer("pseudo_salt")[...]
    out["eta"] = m.eta[...]
    rho = m.rho[...]
    out["rho"] = [[rho[i][j] for j in range(NY)] for i in range(NX)]
    out["rho"][LAND[0]][LAND[1]] = None      # EOS on a land column: not compared
    return out


def test_post_create_plus_rebaseline_matches_staged_for_ten_steps():
    cfg_text = _toml()
    ic = _ic()
    a = _staged(cfg_text, ic)
    try:
        a0 = _read(a)
        a.step(10)
        a10 = _read(a)
    finally:
        a.close()
    b = _post_create(cfg_text, ic)
    try:
        b.rebaseline()
        b0 = _read(b)
        # eta is current straight after the rebaseline: sum(h) - H
        for i in range(NX):
            for j in range(NY):
                if (i, j) == LAND:
                    continue
                want = DEPTH * 1.0e-3 * (i + 1)
                assert abs(b0["eta"][i][j] - want) < 1e-9
        b.step(10)
        b10 = _read(b)
    finally:
        b.close()
    for key in ("h", "u", "v", "T", "S", "P"):
        assert a0[key] == b0[key], f"t=0 {key}"
        assert a10[key] == b10[key], f"t=10 {key}"
    assert a0["rho"] == b0["rho"]
    assert a10["rho"] == b10["rho"]
    assert a10["eta"] == b10["eta"]


def test_without_rebaseline_the_sponge_target_is_stale():
    """Anti-vacuity: the 10-step comparison above must NOT hold without the
    rebaseline (the sponge still relaxes toward the built-in uniform IC)."""
    cfg_text = _toml()
    ic = _ic()
    a = _staged(cfg_text, ic)
    try:
        a.step(10)
        a10 = _read(a)
    finally:
        a.close()
    b = _post_create(cfg_text, ic)
    try:
        b.step(10)
        b10 = _read(b)
    finally:
        b.close()
    assert a10["T"] != b10["T"]


def test_rebaseline_after_step_raises_setup_error():
    with rdb.Model(_toml()) as m:
        m.step(1)
        with pytest.raises(rdb.RdbSetupError):
            m.rebaseline()


def test_rebaseline_keeps_user_set_pseudo_salt():
    cfg_text = _toml()
    ic = _ic()
    m = _post_create(cfg_text, ic)
    try:
        m.tracer("pseudo_salt")[...] = _field(
            NX, NY, NZ, lambda i, j, k: 20.0 + k)
        m.rebaseline()
        p = m.tracer("pseudo_salt")[...]
        s = m.salinity[...]
        assert abs(p[0][0][1] - 21.0) < 1e-12
        assert p[0][0][1] != s[0][0][1]
        # the flag was consumed: a second set of S re-seeds pseudo-salt
        m.salinity[...] = ic["S"]
        m.rebaseline()
        p = m.tracer("pseudo_salt")[...]
        s = m.salinity[...]
        assert p[2][2][1] == s[2][2][1]
    finally:
        m.close()


def test_salt_and_heat_totals_close_after_rebaseline_in_a_closed_basin():
    """A9: after the rebaseline (land mask, ghosts, EOS) one step of a closed,
    unforced basin conserves sum(h*S) and sum(h*T) to roundoff."""
    cfg_text = _toml(closed=True)
    ic = _ic()
    m = _post_create(cfg_text, ic)
    try:
        m.rebaseline()

        def totals():
            h = m.h[...]
            t = m.temperature[...]
            s = m.salinity[...]
            tot = [0.0, 0.0, 0.0]
            for i in range(NX):
                for j in range(NY):
                    if (i, j) == LAND:      # land column: h = 0, no content
                        continue
                    for k in range(NZ):
                        tot[0] += h[i][j][k]
                        tot[1] += h[i][j][k] * t[i][j][k]
                        tot[2] += h[i][j][k] * s[i][j][k]
            return tot

        before = totals()
        m.step(1)
        after = totals()
        for b, a in zip(before, after):
            assert abs(a - b) <= 1e-12 * abs(b)
    finally:
        m.close()
