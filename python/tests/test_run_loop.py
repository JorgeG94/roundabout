"""P6 -- the run loop: `Model.run(until=/steps=, callback=, every=)`, and
the D3.4 "no sync() call anywhere" contract. There is no `sync()` method
on `Model` at all (grep of the class confirms it -- see
`test_no_sync_method_exists`); every Field re-checks its generation on
access instead.
"""

from __future__ import annotations

import rdb
from rdb.barotropic import BarotropicSolver
from rdb.grids import FlatBottom, RectilinearGrid


def _small_model():
    grid = RectilinearGrid(size=(6, 5, 2), extent=(6000.0, 5000.0))
    model = rdb.Model(grid=grid, bathymetry=FlatBottom(100.0),
                          timestep=120.0,
                          barotropic=BarotropicSolver(auto=True))
    # See test_curated_model_construction.py's _make(): no curated
    # Diagnostics()/Restart() object exists yet, and the Fortran
    # defaults leave diag/output on with a relative path whose failure
    # is an uncaught ERROR STOP, not a typed exception.
    model.config.ocean_diag.enabled = False
    model.config.output.output_to_file = False
    return model


def test_no_sync_method_exists():
    assert not hasattr(rdb.Model, "sync")
    assert not hasattr(rdb.Model, "sync_host")


def test_run_with_steps_and_callback():
    model = _small_model()
    model.config.time.t_end = 1.0e6
    calls = []
    with model:
        model.run(steps=6, callback=lambda m: calls.append(m.step_count),
                   every=2)
    assert calls == [2, 4, 6]


def test_run_with_until_and_callback_fires_on_time_boundaries():
    model = _small_model()
    model.config.time.t_end = 1.0e6
    times = []
    with model:
        model.run(until=6 * 120.0,
                   callback=lambda m: times.append(m.time), every=240.0)
    # dt=120s: boundaries at 240, 480, 720 s.
    assert times == [240.0, 480.0, 720.0]


def test_run_requires_exactly_one_of_until_or_steps():
    model = _small_model()
    model.config.time.t_end = 1.0e6
    with model:
        try:
            model.run()
            assert False, "expected ValueError"
        except ValueError:
            pass
        try:
            model.run(until=100.0, steps=3)
            assert False, "expected ValueError"
        except ValueError:
            pass


def test_field_reads_update_across_steps_without_any_sync_call():
    """The callback reads state after every step with NO manual sync --
    the field's generation counter silently refreshes on access (D3.4).
    Each read carries the CURRENT `generation` (== step_count at read
    time); asserting that sequence is strictly increasing across the
    run is the observable proxy for "no stale reads", since nothing in
    this test (or anywhere in run()) ever calls a sync method."""
    model = _small_model()
    model.config.time.t_end = 1.0e6
    generations = []
    with model:
        def cb(m):
            generations.append(m.h._gen)
            assert m.kinetic_energy >= 0.0

        model.run(steps=4, callback=cb, every=1)
    assert generations == [1, 2, 3, 4]
