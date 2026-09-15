"""P5/P6 -- `Model(...)` construction-path rules: curated kwargs and a raw
namelist/Config are mutually exclusive entry points (one Model, one
recipe), the single-live-handle guard still holds for curated
construction, and `defer=True` requires an explicit `create()`.
"""

from __future__ import annotations

import pytest

import rdb
from rdb.barotropic import BarotropicSolver
from rdb.grids import FlatBottom, RectilinearGrid


def _grid():
    return RectilinearGrid(size=(4, 4, 2), extent=(4000.0, 4000.0))


def _make(**kwargs):
    model = rdb.Model(grid=_grid(), bathymetry=FlatBottom(100.0),
                          timestep=60.0,
                          barotropic=BarotropicSolver(auto=True), **kwargs)
    # No curated Diagnostics()/Restart() object exists yet (deferred,
    # see the phase report) and the Fortran defaults leave diag/output
    # ON with a relative "./output/..." path; a failed NetCDF write
    # there is an uncaught `ERROR STOP` (P0 review F5), not a typed
    # exception -- kills the whole interpreter. Disable both explicitly
    # so this test suite never depends on the cwd having a writable
    # ./output/ directory.
    model.config.ocean_diag.enabled = False
    model.config.output.output_to_file = False
    return model


def test_curated_and_raw_namelist_together_raises_type_error():
    with pytest.raises(TypeError, match="not both"):
        rdb.Model("&sim_nml sim_type = \"ocean\" /\n", grid=_grid())


def test_curated_construction_defers_create_until_with_block():
    model = _make()
    assert repr(model) == "Model(uncreated, grid=?)"
    with pytest.raises(rdb.NotInitialisedError):
        model.step(1)
    with model:
        assert model.grid_info["nx"] == 4
        model.step(1)
    assert model._closed


def test_defer_true_requires_explicit_create_even_inside_with():
    model = _make(defer=True)
    with model:
        with pytest.raises(rdb.NotInitialisedError):
            model.step(1)
        model.create()
        model.step(1)
        assert model.step_count == 1
    assert model._closed


def test_create_twice_raises():
    model = _make()
    with model:
        with pytest.raises(rdb.NotInitialisedError):
            model.create()


def test_second_curated_model_while_first_open_raises_already_exists():
    model1 = _make()
    with model1:
        with pytest.raises(rdb.AlreadyExistsError):
            _make().create()


def test_config_escape_hatch_and_curated_object_compose_in_one_script():
    """A curated object AND a raw model.config.<group>.<knob>
    assignment compose in one script, and to_namelist() shows both
    (D2.0 rule 1)."""
    from rdb.closures import HorizontalDiffusivity

    model = _make(defer=True, closures=[HorizontalDiffusivity(kappa_h=25.0)])
    model.config.ocean_bdrag.cd = 0.0025  # raw escape-hatch assignment
    text = model.config.to_namelist()
    assert "kappa_h = 25.0" in text          # from the curated object
    assert "cd = 0.0025" in text             # from the raw assignment
    model.close()
