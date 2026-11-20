"""P6 -- `Model.from_toml`: the migration path. Reads a config FILE
in Python (relative paths inside it resolved against its own directory,
no `os.chdir()`), and the resulting `.config` is the SAME TYPED `Config`
a hand-written script would build -- so a user migrates one knob at a
time. A conflict between the file and a curated override raises rather
than silently letting one win.
"""

from __future__ import annotations

import os

import pytest

import rdb
from rdb._knob import MISSING


@pytest.fixture
def toml_path(tmp_path):
    text = (
        """[sim]
sim_type = "ocean"
[grid]
nx = 8
ny = 6
dx = 2000.0
dy = 2000.0
[nonhydrostatic]
nz_layers = 3
[time]
t_end = 86400.0
dt_fixed = 300.0
[ocean.topo]
topo_config = "flat"
max_depth = 200.0
[ocean.bt]
auto_n_inner = true
[tracer]
initial_salinity = 35.0
initial_temperature = 12.0
[ocean.diag]
enabled = false
[output]
output_to_file = false
[ocean.hvisc]
nu_h = 5000.0
smag_ah = true
bound_kh = true
"""
        # bound_kh=true: nu_h=5000 at this fixture's 2 km cells gives a
        # viscous CFL of 0.375 against the configure-time stability
        # audit's 0.125 bound (rdb_ocean_stability_audit.F90) -- bound_kh
        # is the runtime per-cell clamp that makes that configuration
        # safe (see FINDINGS.md's real global-aquaplanet failure, fixed
        # the same way), keeping this fixture's nu_h=5000.0 intact for
        # the assertion in test_from_toml_curated_override_without_conflict_composes.
    )
    p = tmp_path / "case.toml"
    p.write_text(text)
    return p


def test_from_toml_yields_typed_equivalent_config(toml_path):
    model = rdb.Model.from_toml(str(toml_path), defer=True)
    cfg = model.config
    assert cfg.grid.nx == 8
    assert cfg.grid.ny == 6
    assert cfg.grid.dx == 2000.0
    assert cfg.nonhydrostatic.nz_layers == 3
    assert cfg.time.dt_fixed == 300.0
    assert cfg.ocean_topo.topo_config == "flat"
    assert cfg.ocean_topo.max_depth == 200.0
    assert cfg.ocean_hvisc.nu_h == 5000.0
    assert cfg.ocean_hvisc.smag_ah is True
    # A knob the file never set stays MISSING, same as a hand-built Config.
    assert cfg.ocean_hvisc.c_smag is MISSING
    model.close()


def test_from_toml_matches_hand_built_config_on_to_toml(toml_path):
    from rdb import Config

    model = rdb.Model.from_toml(str(toml_path), defer=True)

    hand = Config()
    hand.sim.sim_type = "ocean"
    hand.grid.nx = 8
    hand.grid.ny = 6
    hand.grid.dx = 2000.0
    hand.grid.dy = 2000.0
    hand.nonhydrostatic.nz_layers = 3
    hand.time.t_end = 86400.0
    hand.time.dt_fixed = 300.0
    hand.ocean_topo.topo_config = "flat"
    hand.ocean_topo.max_depth = 200.0
    hand.ocean_bt.auto_n_inner = True
    hand.tracer.initial_salinity = 35.0
    hand.tracer.initial_temperature = 12.0
    hand.ocean_diag.enabled = False
    hand.output.output_to_file = False
    hand.ocean_hvisc.nu_h = 5000.0
    hand.ocean_hvisc.smag_ah = True
    hand.ocean_hvisc.bound_kh = True

    # Same explicit-knob SET (order doesn't matter -- compare as sets).
    got = set(model.config.explicit_knobs())
    want = set(hand.explicit_knobs())
    assert got == want
    model.close()


def test_from_toml_and_model_direct_create_agree(toml_path):
    text = toml_path.read_text()
    m1 = rdb.Model(text)
    grid_info_1 = m1.grid_info
    m1.close()
    with rdb.Model.from_toml(str(toml_path)) as m2:
        assert grid_info_1 == m2.grid_info


def test_from_toml_curated_override_conflict_raises(toml_path):
    from rdb.closures import ConstantViscosity

    with pytest.raises(rdb.ConfigConflictError, match="nu_h"):
        rdb.Model.from_toml(
            str(toml_path), defer=True,
            closures=[ConstantViscosity(nu_h=1.0)])


def test_from_toml_curated_override_without_conflict_composes(toml_path):
    from rdb.forcing import SurfaceHeatFlux

    model = rdb.Model.from_toml(
        str(toml_path), defer=True, forcing=[SurfaceHeatFlux(q=-30.0)])
    assert model.config.ocean_thermo.q_heat == -30.0
    assert model.config.ocean_hvisc.nu_h == 5000.0  # file's own knob intact
    model.close()


def test_from_toml_resolves_relative_paths_against_file_directory(
        tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    bathy_file = sub / "bathy.nc"
    bathy_file.write_bytes(b"not a real netcdf, just needs to exist")

    cfg = tmp_path / "case.toml"
    cfg.write_text(
        """[sim]
sim_type = "ocean"
[grid]
nx = 4
ny = 4
dx = 1000.0
dy = 1000.0
[ocean.topo]
topo_config = "file"
[output]
bathymetry_file = "sub/bathy.nc"
"""
    )
    # cwd is NOT tmp_path -- the point of the test.
    old_cwd = os.getcwd()
    try:
        os.chdir("/")
        model = rdb.Model.from_toml(str(cfg), defer=True)
    finally:
        os.chdir(old_cwd)
    assert model.config.output.bathymetry_file == str(bathy_file.resolve())
    model.close()


def test_from_toml_unknown_table_raises(tmp_path):
    path = tmp_path / "bogus.toml"
    path.write_text("[totally_bogus_group]\nx = 1\n")
    with pytest.raises(rdb.ConfigParseError,
                       match=r"unknown table '\[totally_bogus_group\]'"):
        rdb.Model.from_toml(str(path), defer=True)


def test_from_toml_unknown_key_raises(tmp_path):
    path = tmp_path / "bogus_key.toml"
    path.write_text("[grid]\nnx = 8\ntotally_bogus_key = 1\n")
    with pytest.raises(rdb.ConfigParseError,
                       match="unknown key 'totally_bogus_key' in .grid."):
        rdb.Model.from_toml(str(path), defer=True)


def test_from_toml_keys_are_exact_case(tmp_path):
    # The TOML front end is exact-case: `NX` is not the knob `nx`.
    path = tmp_path / "case.toml"
    path.write_text("[grid]\nNX = 8\n")
    with pytest.raises(rdb.ConfigParseError, match="unknown key 'NX'"):
        rdb.Model.from_toml(str(path), defer=True)


def test_from_toml_accepts_an_integer_for_a_real_knob(tmp_path):
    path = tmp_path / "int_real.toml"
    path.write_text("[grid]\nnx = 8\ndx = 2000\n")
    model = rdb.Model.from_toml(str(path), defer=True)
    assert model.config.grid.dx == 2000.0
    model.close()


def test_from_toml_ignores_comments_naming_tables(tmp_path):
    # A header comment that names a table and quotes a value is prose: it
    # must neither open a phantom table nor set a knob.
    p = tmp_path / "commented.toml"
    p.write_text(
        '# The IC comes from [ocean.zinit] source = "file" (see README).\n'
        "# Apostrophes too: the model's grid.\n"
        '[sim]\nsim_type = "ocean"  # trailing [grid] nx = 99\n'
        "[grid]\nnx = 8\nny = 6\n")
    model = rdb.Model.from_toml(str(p), defer=True)
    assert model.config.grid.nx == 8
    assert model.config.ocean_zinit.source is MISSING
    model.close()


def test_to_toml_round_trips_through_from_toml(tmp_path):
    from rdb import Config

    hand = Config()
    hand.sim.sim_type = "ocean"
    hand.grid.nx = 8
    hand.grid.dx = 2000.5
    hand.ocean_bt.auto_n_inner = True
    hand.ocean_hvisc.nu_h = 5000.0
    hand.output.output_dir = 'a "quoted" dir'
    path = tmp_path / "round.toml"
    hand.write_toml(str(path))
    model = rdb.Model.from_toml(str(path), defer=True)
    assert set(model.config.explicit_knobs()) == set(hand.explicit_knobs())
    assert model.config.to_toml() == hand.to_toml()
    model.close()
