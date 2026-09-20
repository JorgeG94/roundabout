"""P4 -- the generated typed config layer (`rdb.Config`).

Round-trip through the REAL C ABI (matches this suite's convention: it
exercises the real `.so`, not a mock, per `conftest.py`), the
bit-identity/absence property, validation-on-assign, and the
`dead_on_ocean_path` warning.
"""

import warnings

import pytest

import rdb
from rdb import Config
from rdb._config import ALL_GROUPS, N_GROUPS, N_KNOBS
from rdb._errors import ConfigParseError
from rdb._knob import RdbDeadKnobWarning


def test_group_and_knob_counts():
    # 60 groups / 657 knobs total, i.e. exactly what `rdb_nml_json` +
    # `tools/gen_python_config.py` report for the live schema (the last
    # bump: &ocean_vdiff_nml implicit_top_drag, the vdiff k=nz fold for
    # &ocean_tdrag_nml's 8 ice-shelf top-drag knobs).
    assert N_GROUPS == 60
    assert N_KNOBS == 657
    assert "ocean_bc" in ALL_GROUPS


def test_absence_is_the_default():
    """A config with nothing set serialises to nothing -- the
    bit-identity property `read_config_from_string` relies on."""
    cfg = Config()
    assert cfg.to_namelist() == ""
    assert list(cfg.explicit_knobs()) == []


def test_class_level_introspection():
    assert Config.ocean_epbl.mstar.default == 1.2
    assert Config.vcoord.zstar_stretching.dead_on_ocean_path


def test_to_namelist_emits_only_explicit_knobs():
    cfg = Config()
    cfg.ocean_hvisc.nu_h = 5000.0
    cfg.ocean_hvisc.smag_ah = True
    text = cfg.to_namelist()
    assert "&ocean_hvisc_nml" in text
    assert "nu_h = 5000.0" in text
    assert "smag_ah = .true." in text
    # Only ONE group emitted (nothing else was ever touched).
    assert text.count("&") == 1
    assert list(cfg.explicit_knobs()) == [
        ("ocean_hvisc", "nu_h", 5000.0),
        ("ocean_hvisc", "smag_ah", True),
    ]


def test_validation_range_error_is_typed():
    cfg = Config()
    with pytest.raises(ConfigParseError):
        cfg.grid.nx = 0  # below vmin=1


def test_validation_enum_error_is_typed():
    cfg = Config()
    with pytest.raises(ConfigParseError):
        cfg.ocean_coriolis.form = "not_a_real_form"


def test_validation_type_error_is_typed():
    cfg = Config()
    with pytest.raises(ConfigParseError):
        cfg.grid.nx = 5.5  # int knob, got a float


def test_enum_stores_canonical_spelling_case_insensitive():
    cfg = Config()
    cfg.vcoord.vcoord_type = "ZSTAR_FULL"
    assert cfg.vcoord.vcoord_type == "zstar_full"


def test_dead_knob_warns_on_assignment():
    cfg = Config()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cfg.vcoord.zstar_stretching = "log"
    assert any(issubclass(w.category, RdbDeadKnobWarning) for w in caught)


def test_live_knob_does_not_warn():
    cfg = Config()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cfg.ocean_hvisc.nu_h = 1000.0
    assert not any(issubclass(w.category, RdbDeadKnobWarning)
                    for w in caught)


def test_model_accepts_config_and_matches_string_path(ocean_nml):
    """A Config built to say the exact same thing as `ocean_nml` (the
    fixture also driving the Fortran ctest reference) is accepted and
    produces the same grid."""
    cfg = Config()
    cfg.sim.sim_type = "ocean"
    cfg.grid.nx = 8
    cfg.grid.ny = 6
    cfg.grid.dx = 2000.0
    cfg.grid.dy = 2000.0
    cfg.nonhydrostatic.nz_layers = 3
    cfg.time.t_end = 86400.0
    cfg.time.dt_fixed = 300.0
    cfg.ocean_topo.topo_config = "flat"
    cfg.ocean_bt.auto_n_inner = True
    cfg.tracer.initial_salinity = 35.0
    cfg.tracer.initial_temperature = 12.0
    cfg.ocean_diag.enabled = False

    with rdb.Model(cfg) as m_cfg:
        grid_from_config = m_cfg.grid_info
        assert m_cfg.config is cfg

    with rdb.Model(ocean_nml) as m_str:
        grid_from_string = m_str.grid_info

    assert grid_from_config == grid_from_string


def test_model_pending_accepts_config(ocean_nml):
    cfg = Config()
    cfg.sim.sim_type = "ocean"
    cfg.grid.nx = 8
    cfg.grid.ny = 6
    cfg.grid.dx = 2000.0
    cfg.grid.dy = 2000.0
    cfg.nonhydrostatic.nz_layers = 3
    cfg.time.t_end = 86400.0
    cfg.time.dt_fixed = 300.0
    cfg.ocean_topo.topo_config = "flat"
    cfg.ocean_bt.auto_n_inner = True
    cfg.tracer.initial_salinity = 35.0
    cfg.tracer.initial_temperature = 12.0
    cfg.ocean_diag.enabled = False

    m = rdb.Model.pending(cfg)
    try:
        m.finalize()
        assert m.grid_info["nx"] == 8
    finally:
        m.close()


def test_model_config_property_on_string_build_is_empty_config(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        cfg = m.config
        assert isinstance(cfg, Config)
        assert cfg.to_namelist() == ""


def test_model_rejects_unrecognised_config_type():
    with pytest.raises(TypeError):
        rdb.Model(12345)


def test_ocean_bc_stub_reachable_and_round_trips():
    cfg = Config()
    cfg.ocean_bc.west = "open"
    cfg.ocean_bc.mask_wall_velocity = False
    text = cfg.to_namelist()
    assert "&ocean_bc_nml" in text
    assert 'west = "open"' in text
    assert "mask_wall_velocity = .false." in text
