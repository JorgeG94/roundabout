"""P7 -- first-class diagnostics.

Covers all five pieces of the phase:
  1. create() with a missing output directory does not abort the
     interpreter (the headline F1 fix) -- and DOES abort typed-error-only
     when the directory cannot be created at all.
  2. available_diagnostics() / model.diagnostics.available / .selected.
  3. Diagnostics(enabled=False) actually suppresses output.
  4. model.diagnostic(name) matches the value the NetCDF path would write
     (compared against sum(h) - b, the same formula the diag manager's
     own fill_ssh uses).
  5. The NaN sentinel change is exercised at the Fortran unit level
     (tests/test_ocean_diag.F90::fill_temperature_vanished_layer_is_nan);
     this file only asserts a NORMAL diagnostic read is not NaN, so a
     regression in the other direction (everything reads NaN) is also
     caught here.
"""

from __future__ import annotations

import math

import pytest

import rdb


def _diag_nml(output_dir, *, filename="p7test", enabled=True):
    on = ".true." if enabled else ".false."
    return (
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /\n"
        "&nonhydrostatic_nml nz_layers = 3 /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        '&ocean_topo_nml topo_config = "flat", max_depth = 200.0 /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        f'&ocean_diag_nml enabled = {on}, dt_out = 300.0, '
        f'filename = "{filename}" /\n'
        f'&output_nml output_to_file = {on}, output_dir = "{output_dir}" /\n'
    )


# ---- 1. the headline crash fix -----------------------------------------

def test_missing_output_dir_is_created_not_a_crash(tmp_path):
    """A single missing path component (the common case: a fresh
    checkout with no ./output/) must be created transparently, and
    create() must succeed -- not abort the interpreter, not raise."""
    out_dir = tmp_path / "output"
    assert not out_dir.exists()
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        m.step(1)
    assert out_dir.is_dir()
    assert (out_dir / "p7test_rank_000000.nc").is_file()


def test_deeply_missing_output_dir_raises_typed_error_not_a_crash(tmp_path):
    """A path whose GRANDPARENT is also missing cannot be created by a
    single mkdir -- this must surface as a typed RdbIOError (proving
    the process is still alive to raise it), never an interpreter
    abort."""
    out_dir = tmp_path / "does" / "not" / "exist" / "output"
    with pytest.raises(rdb.RdbIOError):
        rdb.Model(_diag_nml(str(out_dir)))
    # The interpreter is provably still alive: a second, valid Model
    # can still be constructed in the same process.
    good_dir = tmp_path / "output2"
    with rdb.Model(_diag_nml(str(good_dir))) as m:
        m.step(1)


# ---- 2. discoverability -------------------------------------------------

def test_available_diagnostics_is_nonempty_with_known_names():
    names = rdb.available_diagnostics()
    assert names, "available_diagnostics() must not be empty"
    for known in ("SSH", "temperature", "salinity", "u", "v", "KE"):
        assert known in names
    # Derived-catalog names are also reachable.
    assert "vorticity_z" in names
    assert "ke_total" in names


def test_diagnostics_available_matches_module_level(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        assert m.diagnostics.available == rdb.available_diagnostics()


def test_diagnostics_selected_reflects_live_registration(tmp_path):
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        selected = m.diagnostics.selected
        for name in ("SSH", "temperature", "salinity", "u", "v", "KE"):
            assert name in selected
        # Never-requested derived diagnostics are NOT auto-registered.
        assert "vorticity_z" not in selected


# ---- 3. the curated Diagnostics() off switch -----------------------------

def test_diagnostics_enabled_false_suppresses_output(tmp_path):
    out_dir = tmp_path / "output"
    cfg = rdb.Config()
    cfg.grid.nx = 8
    cfg.grid.ny = 6
    cfg.grid.dx = 2000.0
    cfg.grid.dy = 2000.0
    cfg.nonhydrostatic.nz_layers = 3
    cfg.time.t_end = 86400.0
    cfg.time.dt_fixed = 300.0
    cfg.ocean_topo.topo_config = "flat"
    cfg.ocean_topo.max_depth = 200.0
    cfg.ocean_bt.auto_n_inner = True
    cfg.tracer.initial_salinity = 35.0
    cfg.tracer.initial_temperature = 12.0
    cfg.output.output_dir = str(out_dir)
    rdb.Diagnostics(enabled=False).apply(cfg)

    with rdb.Model(cfg) as m:
        m.step(2)
        assert m.diagnostics.selected == []
    # No output directory should even have been created.
    assert not out_dir.exists()


def test_diagnostics_bad_token_raises_before_create():
    """A malformed &ocean_diag_nml diags token is a Python ValueError at
    construction time -- parse_one_spec_token's unknown-attribute branch
    (rdb_ocean_diag.F90) is one of the paths P0.1 did not convert, so
    reaching it via the raw namelist would still error stop today."""
    try:
        rdb.Diagnostics(names=["temperature:not_a_real_attribute"])
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "temperature:not_a_real_attribute" in str(exc)


def test_diagnostics_names_and_cadence_reach_the_config():
    cfg = rdb.Config()
    rdb.Diagnostics(names=["KE:off", "temperature:6h:mean"],
                       every=1800.0, filename="mydiag",
                       precision="single").apply(cfg)
    text = cfg.to_namelist()
    assert "diags = " in text
    assert "KE:off" in text
    assert "temperature:6h:mean" in text
    assert "dt_out = 1800.0" in text
    assert 'filename = "mydiag"' in text
    assert 'output_precision = "single"' in text


# ---- 4. in-memory access matches the diagnostic-stream formula ----------

def test_diagnostic_ssh_matches_sum_h_minus_b(tmp_path):
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        m.step(1)
        ssh_diag = m.diagnostic("SSH")
        nx, ny, nz = ssh_diag.shape
        assert nz == 1
        ssh_formula = m.ssh  # sum_k h - b, the same formula fill_ssh uses
        diag_vals = ssh_diag[...]
        for i in range(nx):
            for j in range(ny):
                assert abs(diag_vals[i][j][0] - ssh_formula[i][j]) < 1e-9


def test_diagnostic_temperature_matches_tracer_concentration(tmp_path):
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        m.step(1)
        t_diag = m.diagnostic("temperature")
        t_field = m.temperature
        nx, ny, nz = t_diag.shape
        diag_vals = t_diag[...]
        field_vals = t_field[...]
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    assert abs(diag_vals[i][j][k] - field_vals[i][j][k]) < 1e-9


def test_diagnostic_not_registered_raises_typed_error(tmp_path):
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        try:
            m.diagnostic("does_not_exist")
            assert False, "expected a not-found error"
        except rdb.TracerNotFoundError as exc:
            assert "does_not_exist" in str(exc)


def test_diagnostic_field_is_read_only(tmp_path):
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        m.step(1)
        ssh_diag = m.diagnostic("SSH")
        try:
            ssh_diag[...] = 1.0
            assert False, "expected RdbReadOnlyError"
        except rdb.RdbReadOnlyError:
            pass


# ---- 5. NaN sentinel (see also the Fortran unit test) --------------------

def test_diagnostic_normal_cell_is_not_nan(tmp_path):
    """The counterpart of the Fortran vanished-layer regression test: an
    ordinary, non-vanished cell must NOT read NaN -- guards against a
    fix that is too aggressive (everything reads NaN)."""
    out_dir = tmp_path / "output"
    with rdb.Model(_diag_nml(str(out_dir))) as m:
        m.step(1)
        t_diag = m.diagnostic("temperature")
        for val in [v for row in t_diag[...] for col in row for v in col]:
            assert not math.isnan(val)
