"""P1-C -- tests for the curated routes added to close the escape-hatch
gaps measured by porting the four validation_examples/ocean scripts
(acc_channel[.py]/acc_channel_quiescent.py/seamount.py/
double_gyre_mom6.py) to the Python API.

Covers: KPP(enabled=False) (the boundary-layer off-switch + the
NoBoundaryLayerWarning it can trigger), EPBL() turning KPP off itself,
PacanowskiPhilander(enabled=False), Smagorinsky/Leith's nu_h/
kh_vel_scale/nu_4 constant-floor knobs (and their biharmonic_coeff()
composition), NearSurfaceViscosity/HarmonicFaceThickness/
ThermoSubcycling, Diagnostics' output_dir/status_interval/log_level,
Duration, RectilinearGrid's radius=, StratifiedInitialCondition, and
Thermodynamics.
"""

from __future__ import annotations

import warnings

import pytest

import rdb
from rdb import _compose
from rdb.closures import (
    EPBL,
    KPP,
    HarmonicFaceThickness,
    Leith,
    NearSurfaceViscosity,
    PacanowskiPhilander,
    Smagorinsky,
    ThermoSubcycling,
)
from rdb.diagnostics import Diagnostics, Duration
from rdb.forcing import StratifiedInitialCondition, Thermodynamics
from rdb.grids import RectilinearGrid


def _cart_grid():
    return RectilinearGrid(size=(4, 4, 2), extent=(4000.0, 4000.0))


# ----------------------------------------------------------------------
# KPP(enabled=False) / EPBL() disabling KPP itself
# ----------------------------------------------------------------------


def test_kpp_enabled_false_turns_off_use_kpp():
    cfg = rdb.Config()
    KPP(enabled=False).apply(cfg)
    assert cfg.ocean_vmix.use_kpp is False
    # every other keyword is ignored -- nothing else should be set
    assert "kpp_ri_crit" not in cfg.ocean_vmix._explicit_names()


def test_kpp_default_still_enabled():
    cfg = rdb.Config()
    KPP().apply(cfg)
    assert cfg.ocean_vmix.use_kpp is True
    assert cfg.ocean_vmix.kpp_ri_crit == 0.3


def test_epbl_turns_off_kpp_itself():
    """EPBL() alone (no KPP object at all) must leave the composed
    config with use_kpp = False -- selecting EPBL structurally implies
    deselecting KPP (KPP defaults ON in the schema)."""
    composed = _compose.compose(grid=_cart_grid(), closures=[EPBL()])
    assert composed.config.ocean_vmix.use_kpp is False
    assert composed.config.ocean_epbl.enable is True


def test_kpp_enabled_false_with_epbl_does_not_raise():
    """KPP(enabled=False) alongside EPBL(...) is legal (if redundant --
    EPBL already turns use_kpp off itself)."""
    composed = _compose.compose(
        grid=_cart_grid(), closures=[EPBL(), KPP(enabled=False)])
    assert composed.config.ocean_vmix.use_kpp is False
    assert composed.config.ocean_epbl.enable is True


def test_kpp_true_with_epbl_still_raises():
    with pytest.raises(rdb.ConfigConflictError, match="mutually exclusive"):
        _compose.compose(grid=_cart_grid(), closures=[KPP(), EPBL()])


def test_kpp_enabled_false_without_epbl_warns_no_boundary_layer():
    with pytest.warns(rdb.NoBoundaryLayerWarning, match="no surface boundary-layer"):
        _compose.compose(grid=_cart_grid(), closures=[KPP(enabled=False)])


def test_neither_kpp_nor_epbl_mentioned_does_not_warn():
    """Omitting BOTH objects entirely leaves the schema's own KPP-on-
    by-default in effect -- NOT the "neither" case, so no warning."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        _compose.compose(grid=_cart_grid(), closures=[])
    assert not any(issubclass(w.category, rdb.NoBoundaryLayerWarning)
                   for w in record)


# ----------------------------------------------------------------------
# PacanowskiPhilander(enabled=False)
# ----------------------------------------------------------------------


def test_pacanowski_philander_enabled_false_turns_off_use_closure():
    cfg = rdb.Config()
    PacanowskiPhilander(enabled=False).apply(cfg)
    assert cfg.ocean_vmix.use_closure is False


def test_pacanowski_philander_default_enabled():
    cfg = rdb.Config()
    PacanowskiPhilander().apply(cfg)
    assert cfg.ocean_vmix.use_closure is True
    assert cfg.ocean_vmix.pp81_nu0 == 1e-2


# ----------------------------------------------------------------------
# Smagorinsky / Leith constant-floor knobs
# ----------------------------------------------------------------------


def test_smagorinsky_constant_nu_h_floor_composes():
    """CLAUDE.md's own production envelope: Smagorinsky + a constant
    nu_h floor, in the SAME &ocean_hvisc_nml group."""
    cfg = rdb.Config()
    Smagorinsky(C=0.15, nu_h=10000.0, kh_vel_scale=3.0e-3).apply(cfg)
    assert cfg.ocean_hvisc.lateral_closure == "smagorinsky"
    assert cfg.ocean_hvisc.nu_h == 10000.0
    assert cfg.ocean_hvisc.kh_vel_scale == 3.0e-3
    text = cfg.to_namelist()
    assert 'lateral_closure = "smagorinsky"' in text
    assert "nu_h = 10000.0" in text


def test_smagorinsky_biharmonic_coeff_reflects_nu_4_floor():
    """The MEKE-backscatter compose-time check inspects
    biharmonic_coeff() -- a non-zero nu_4 constant floor is a real
    dissipation backstop even with biharmonic=False."""
    s = Smagorinsky(biharmonic=False, nu_4=1e9)
    assert s.biharmonic_coeff() == 1e9

    s2 = Smagorinsky(biharmonic=True, C_biharmonic=0.06, nu_4=0.0)
    assert s2.biharmonic_coeff() == pytest.approx(0.06)

    s3 = Smagorinsky(biharmonic=True, C_biharmonic=0.06, nu_4=1e9)
    assert s3.biharmonic_coeff() == 1e9


def test_smagorinsky_nu_4_floor_satisfies_meke_backscatter():
    from rdb.closures import MEKE, GentMcWilliams

    composed = _compose.compose(
        grid=_cart_grid(),
        closures=[GentMcWilliams(kappa=100.0),
                  MEKE(backscatter=True, backscatter_ku=1.0),
                  Smagorinsky(biharmonic=False, nu_4=1e9)])
    assert composed.config.ocean_meke.backscatter is True


def test_leith_constant_nu_h_floor_composes():
    cfg = rdb.Config()
    Leith(C=1.0, nu_h=500.0, nu_4=1e8).apply(cfg)
    assert cfg.ocean_hvisc.nu_h == 500.0
    assert cfg.ocean_hvisc.nu_4 == 1e8


# ----------------------------------------------------------------------
# NearSurfaceViscosity / HarmonicFaceThickness / ThermoSubcycling
# ----------------------------------------------------------------------


def test_near_surface_viscosity_writes_kv_ml_invz2_and_hmix_fixed():
    cfg = rdb.Config()
    NearSurfaceViscosity(kv=0.01, hmix=20.0).apply(cfg)
    assert cfg.ocean_vmix.kv_ml_invz2 == 0.01
    assert cfg.ocean_vmix.hmix_fixed == 20.0


def test_harmonic_face_thickness_writes_harmonic_visc():
    cfg = rdb.Config()
    HarmonicFaceThickness().apply(cfg)
    assert cfg.ocean_vmix.harmonic_visc is True


def test_thermo_subcycling_writes_dt_therm_ratio():
    cfg = rdb.Config()
    ThermoSubcycling(ratio=2).apply(cfg)
    assert cfg.ocean_vmix.dt_therm_ratio == 2


# ----------------------------------------------------------------------
# Diagnostics: output_dir / status_interval / log_level
# ----------------------------------------------------------------------


def test_diagnostics_output_dir_status_interval_log_level_round_trip():
    cfg = rdb.Config()
    Diagnostics(enabled=True, filename="foo", every=1.0,
                output_dir="./out_foo", status_interval=3600.0,
                log_level="warning").apply(cfg)
    assert cfg.output.output_dir == "./out_foo"
    assert cfg.logging.status_interval == 3600.0
    assert cfg.logging.log_level == "warning"
    text = cfg.to_namelist()
    assert 'output_dir     = "./out_foo"' in text or "./out_foo" in text
    assert "status_interval = 3600.0" in text
    assert 'log_level = "warning"' in text


def test_diagnostics_output_dir_status_interval_applied_even_when_disabled():
    """output_dir/status_interval/log_level are console/output plumbing,
    not diagnostic selection -- they must survive enabled=False."""
    cfg = rdb.Config()
    Diagnostics(enabled=False, output_dir="./out_bar",
                status_interval=60.0, log_level="debug").apply(cfg)
    assert cfg.ocean_diag.enabled is False
    assert cfg.output.output_dir == "./out_bar"
    assert cfg.logging.status_interval == 60.0
    assert cfg.logging.log_level == "debug"


def test_diagnostics_invalid_log_level_raises():
    with pytest.raises(ValueError, match="log_level"):
        Diagnostics(log_level="nonsense")


# ----------------------------------------------------------------------
# Duration
# ----------------------------------------------------------------------


def test_duration_writes_t_end_and_time_unit():
    cfg = rdb.Config()
    Duration(t_end=10.0, time_unit="day").apply(cfg)
    assert cfg.time.t_end == 10.0
    assert cfg.time.time_unit == "day"


def test_duration_invalid_time_unit_raises():
    with pytest.raises(ValueError, match="time_unit"):
        Duration(t_end=1.0, time_unit="fortnight")


def test_model_duration_kwarg_composes():
    composed = _compose.compose(
        grid=_cart_grid(), duration=Duration(t_end=5.0, time_unit="hr"))
    assert composed.config.time.t_end == 5.0
    assert composed.config.time.time_unit == "hr"


# ----------------------------------------------------------------------
# RectilinearGrid radius=
# ----------------------------------------------------------------------


def test_rectilinear_grid_radius_none_leaves_rad_earth_unset():
    grid = RectilinearGrid(size=(4, 4, 2), extent=(4000.0, 4000.0))
    cfg = rdb.Config()
    grid.apply(cfg)
    assert "rad_earth" not in cfg.ocean_grid._explicit_names()
    text = cfg.to_namelist()
    assert "rad_earth" not in text


def test_rectilinear_grid_degrees_with_radius_writes_rad_earth():
    grid = RectilinearGrid(size=(44, 40, 2), extent=(22.0, 20.0),
                            axis_units="degrees", radius=6.378e6)
    cfg = rdb.Config()
    grid.apply(cfg)
    assert cfg.ocean_grid.axis_units == "degrees"
    assert cfg.ocean_grid.len_lon == 22.0
    assert cfg.ocean_grid.rad_earth == 6.378e6


# ----------------------------------------------------------------------
# StratifiedInitialCondition
# ----------------------------------------------------------------------


def test_stratified_initial_condition_writes_both():
    cfg = rdb.Config()
    StratifiedInitialCondition(T_surface=8.0, T_bottom=2.0).apply(cfg)
    assert cfg.tracer.T_init_surface == 8.0
    assert cfg.tracer.T_init_bottom == 2.0


def test_stratified_initial_condition_requires_both_together():
    with pytest.raises(ValueError, match="TOGETHER"):
        StratifiedInitialCondition(T_surface=8.0)
    with pytest.raises(ValueError, match="TOGETHER"):
        StratifiedInitialCondition(T_bottom=2.0)


def test_stratified_initial_condition_none_is_a_no_op():
    cfg = rdb.Config()
    StratifiedInitialCondition().apply(cfg)
    text = cfg.to_namelist()
    assert "T_init_surface" not in text
    assert "T_init_bottom" not in text


# ----------------------------------------------------------------------
# Thermodynamics
# ----------------------------------------------------------------------


def test_thermodynamics_enabled_false_writes_flag():
    cfg = rdb.Config()
    Thermodynamics(enabled=False).apply(cfg)
    assert cfg.ocean_thermo.enable_thermodynamics is False


def test_thermodynamics_default_enabled():
    cfg = rdb.Config()
    Thermodynamics().apply(cfg)
    assert cfg.ocean_thermo.enable_thermodynamics is True
