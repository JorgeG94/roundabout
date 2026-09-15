"""P5 D2.7/D2.8 -- the compose-time checks a curated `closures=[...]` list
buys over 612 flat knobs. Every check here raises before `create()` is
ever called (text-only composition via `rdb._compose.compose`), naming
the actual problem -- the thing 06_python_surface_design.md calls out as
"the single clearest thing the curated layer buys".
"""

from __future__ import annotations

import pytest

import rdb
from rdb import _compose
from rdb.closures import (
    KPP,
    ConstantViscosity,
    EPBL,
    GentMcWilliams,
    HenyeyBackground,
    BryanLewisBackground,
    ImplicitVerticalFriction,
    LinearDrag,
    MEKE,
    Smagorinsky,
)
from rdb.grids import RectilinearGrid, LatitudeLongitudeGrid, Topology


def _cart_grid():
    return RectilinearGrid(size=(4, 4, 2), extent=(4000.0, 4000.0))


def _spherical_grid():
    return LatitudeLongitudeGrid(size=(4, 4, 2), longitude=(0.0, 4.0),
                                  latitude=(0.0, 4.0))


def test_kpp_xor_epbl_raises():
    with pytest.raises(rdb.ConfigConflictError, match="mutually exclusive"):
        _compose.compose(grid=_cart_grid(), closures=[KPP(), EPBL()])


def test_bryan_lewis_xor_henyey_raises():
    with pytest.raises(rdb.ConfigConflictError, match="mutually exclusive"):
        _compose.compose(
            grid=_spherical_grid(),
            closures=[BryanLewisBackground(), HenyeyBackground()])


def test_henyey_on_cartesian_grid_raises():
    with pytest.raises(rdb.RdbUnsupportedError, match="non-cartesian"):
        _compose.compose(grid=_cart_grid(),
                          closures=[HenyeyBackground()])


def test_meke_backscatter_without_biharmonic_raises():
    with pytest.raises(rdb.ConfigConflictError, match="NON-ZERO"):
        _compose.compose(
            grid=_cart_grid(),
            closures=[GentMcWilliams(kappa=100.0),
                      MEKE(backscatter=True),
                      Smagorinsky(biharmonic=False)])


def test_meke_backscatter_with_zero_biharmonic_coeff_raises():
    with pytest.raises(rdb.ConfigConflictError, match="NON-ZERO"):
        _compose.compose(
            grid=_cart_grid(),
            closures=[GentMcWilliams(kappa=100.0),
                      MEKE(backscatter=True),
                      Smagorinsky(biharmonic=True, C_biharmonic=0.0)])


def test_meke_backscatter_with_nonzero_biharmonic_composes():
    composed = _compose.compose(
        grid=_cart_grid(),
        closures=[GentMcWilliams(kappa=100.0),
                  MEKE(backscatter=True, backscatter_ku=1.0),
                  Smagorinsky(biharmonic=True, C_biharmonic=0.06)])
    text = composed.config.to_namelist()
    assert "backscatter = .true." in text


def test_meke_without_gm_raises():
    with pytest.raises(rdb.ConfigConflictError, match="GentMcWilliams"):
        _compose.compose(grid=_cart_grid(), closures=[MEKE()])


def test_implicit_vertical_friction_bbl_glue_missing_all_three():
    with pytest.raises(rdb.ConfigConflictError,
                        match="harmonic_thickness.*drag.*LinearDrag|"
                              "THREE prerequisites"):
        _compose.compose(
            grid=_cart_grid(),
            closures=[ImplicitVerticalFriction(bbl_glue=True)])


def test_implicit_vertical_friction_bbl_glue_missing_linear_drag_only():
    with pytest.raises(rdb.ConfigConflictError, match="LinearDrag"):
        _compose.compose(
            grid=_cart_grid(),
            closures=[ImplicitVerticalFriction(
                bbl_glue=True, harmonic_thickness=True, drag=True)])


def test_implicit_vertical_friction_bbl_glue_satisfied():
    composed = _compose.compose(
        grid=_cart_grid(),
        closures=[ImplicitVerticalFriction(
            bbl_glue=True, harmonic_thickness=True, drag=True),
            LinearDrag(r=1e-4)])
    text = composed.config.to_namelist()
    assert "bbl_glue = .true." in text
    assert 'form = "linear"' in text


def test_constant_viscosity_warns_above_biharmonic_cfl_cap():
    with pytest.warns(UserWarning, match="biharmonic"):
        ConstantViscosity(nu_4=5e8)


def test_sadourny_energy_only_knobs_raise_off_energy_form():
    from rdb.coriolis import Sadourny
    with pytest.raises(rdb.RdbUnsupportedError, match="energy-form-only"):
        Sadourny(form="enstrophy", bound_coriolis=True)


def test_grid_topology_reaches_ocean_bc_without_a_boundaries_object():
    """A grid's PERIODIC topology must reach &ocean_bc_nml even when the
    user passes no Boundaries() object.

    Regression for the worst silent-wrongness bug of the P5 surface: a grid
    built with topology=(PERIODIC, BOUNDED) derived nghost correctly and set
    .periodic_x, but nothing wrote west/east -- so the edges stayed at their
    "wall" default and a re-entrant channel became a CLOSED BASIN that ran,
    conserved mass, and was wrong. TripolarGrid escaped it only because it
    also sets .forced_edges.
    """
    from rdb.grids import RectilinearGrid, Topology
    from rdb._compose import compose

    grid = RectilinearGrid(size=(32, 16, 4), extent=(3.2e5, 1.6e5, 4.0e3),
                           topology=(Topology.PERIODIC, Topology.BOUNDED))
    composed = compose(grid=grid)
    cfg = getattr(composed, "cfg", None) or composed.config

    assert cfg.ocean_bc.west == "periodic"
    assert cfg.ocean_bc.east == "periodic"
    # y is BOUNDED -- must NOT be silently made periodic
    text = cfg.to_namelist()
    assert 'south = "periodic"' not in text
    assert 'north = "periodic"' not in text


def test_grid_topology_periodic_y_reaches_ocean_bc():
    """The y-axis half of the same rule."""
    from rdb.grids import RectilinearGrid, Topology
    from rdb._compose import compose

    grid = RectilinearGrid(size=(32, 16, 4), extent=(3.2e5, 1.6e5, 4.0e3),
                           topology=(Topology.BOUNDED, Topology.PERIODIC))
    composed = compose(grid=grid)
    cfg = getattr(composed, "cfg", None) or composed.config

    assert cfg.ocean_bc.south == "periodic"
    assert cfg.ocean_bc.north == "periodic"
