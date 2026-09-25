"""P6 -- `Model.from_namelist`: the migration path. Reads a namelist FILE
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
def nml_path(tmp_path):
    text = (
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /\n"
        "&nonhydrostatic_nml nz_layers = 3 /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        '&ocean_topo_nml topo_config = "flat", max_depth = 200.0 /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        "&ocean_diag_nml enabled = .false. /\n"
        "&output_nml output_to_file = .false. /\n"
        # bound_kh=.true.: nu_h=5000 at this fixture's 2 km cells gives a
        # viscous CFL of 0.375 against the configure-time stability
        # audit's 0.125 bound (rdb_ocean_stability_audit.F90) -- bound_kh
        # is the runtime per-cell clamp that makes that configuration
        # safe (see FINDINGS.md's real global-aquaplanet failure, fixed
        # the same way), keeping this fixture's nu_h=5000.0 intact for
        # the assertion in test_from_namelist_curated_override_without_conflict_composes.
        "&ocean_hvisc_nml nu_h = 5000.0, smag_ah = .true., bound_kh = .true. /\n"
    )
    p = tmp_path / "case.nml"
    p.write_text(text)
    return p


def test_from_namelist_yields_typed_equivalent_config(nml_path):
    model = rdb.Model.from_namelist(str(nml_path), defer=True)
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


def test_from_namelist_matches_hand_built_config_on_to_namelist(nml_path):
    from rdb import Config

    model = rdb.Model.from_namelist(str(nml_path), defer=True)

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


def test_from_namelist_and_model_direct_create_agree(nml_path):
    text = nml_path.read_text()
    m1 = rdb.Model(text)
    grid_info_1 = m1.grid_info
    m1.close()
    with rdb.Model.from_namelist(str(nml_path)) as m2:
        assert grid_info_1 == m2.grid_info


def test_from_namelist_curated_override_conflict_raises(nml_path):
    from rdb.closures import ConstantViscosity

    with pytest.raises(rdb.ConfigConflictError, match="nu_h"):
        rdb.Model.from_namelist(
            str(nml_path), defer=True,
            closures=[ConstantViscosity(nu_h=1.0)])


def test_from_namelist_curated_override_without_conflict_composes(nml_path):
    from rdb.forcing import SurfaceHeatFlux

    model = rdb.Model.from_namelist(
        str(nml_path), defer=True, forcing=[SurfaceHeatFlux(q=-30.0)])
    assert model.config.ocean_thermo.q_heat == -30.0
    assert model.config.ocean_hvisc.nu_h == 5000.0  # file's own knob intact
    model.close()


def test_from_namelist_resolves_relative_paths_against_file_directory(
        tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    bathy_file = sub / "bathy.nc"
    bathy_file.write_bytes(b"not a real netcdf, just needs to exist")

    nml = tmp_path / "case.nml"
    nml.write_text(
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 4, ny = 4, dx = 1000.0, dy = 1000.0 /\n"
        '&ocean_topo_nml topo_config = "file" /\n'
        '&output_nml bathymetry_file = "sub/bathy.nc" /\n'
    )
    # cwd is NOT tmp_path -- the point of the test.
    old_cwd = os.getcwd()
    try:
        os.chdir("/")
        model = rdb.Model.from_namelist(str(nml), defer=True)
    finally:
        os.chdir(old_cwd)
    assert model.config.output.bathymetry_file == str(bathy_file.resolve())
    model.close()


def test_from_namelist_unknown_group_raises():
    p_text = "&totally_bogus_group_nml x = 1 /\n"
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".nml", delete=False) as f:
        f.write(p_text)
        path = f.name
    try:
        with pytest.raises(rdb.ConfigParseError, match="unknown namelist group"):
            rdb.Model.from_namelist(path, defer=True)
    finally:
        os.unlink(path)


def test_from_namelist_ignores_group_names_inside_comments(tmp_path):
    # A header comment that names a group and quotes a value (as
    # validation_examples/ocean/global_1deg/global_1deg_unforced.nml does)
    # is prose: it must neither open a phantom group nor set a knob.
    p = tmp_path / "commented.nml"
    p.write_text(
        '! The IC comes from &ocean_zinit_nml source = "file" (see README).\n'
        "! Apostrophes too: the model's grid.\n"
        '&sim_nml sim_type = "ocean" /  ! trailing &grid_nml nx = 99\n'
        "&grid_nml nx = 8, ny = 6 /\n")
    model = rdb.Model.from_namelist(str(p), defer=True)
    assert model.config.grid.nx == 8
    assert model.config.ocean_zinit.source is MISSING
    model.close()
