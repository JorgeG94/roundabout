"""P2.5 two-phase create: pending -> stage_bathymetry -> finalize, and the
D6.2 bathymetry-sign guard -- the single most dangerous argument in this
API (see 06_python_surface_design.md D6.2)."""

import rdb


def _inject_nml(max_depth=999.0):
    return (
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /\n"
        "&nonhydrostatic_nml nz_layers = 3 /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        f'&ocean_topo_nml topo_config = "flat", max_depth = {max_depth} /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        "&ocean_diag_nml enabled = .false. /\n"
        "&output_nml output_to_file = .false. /\n"
    )


def test_pending_stage_finalize_roundtrip():
    depth = [[1234.0] * 6 for _ in range(8)]
    m = rdb.Model.pending(_inject_nml())
    try:
        m.stage_bathymetry(depth, rdb.BATHY_DEPTH_POSITIVE_DOWN)
        m.finalize()
        b = m.bathymetry[...]
        assert b[0][0] == 1234.0
        assert b[7][5] == 1234.0
    finally:
        m.close()


def test_all_negative_gebco_style_bathymetry_raises_sign_error():
    """A GEBCO-style array (negative-down) staged under the WRONG
    convention (claiming it is already positive-down) normalises to zero
    wet cells and must raise BathymetrySignError -- not silently produce
    an all-land quiescent run (D6.2)."""
    gebco = [[-4000.0] * 6 for _ in range(8)]
    m = rdb.Model.pending(_inject_nml())
    try:
        try:
            m.stage_bathymetry(gebco, rdb.BATHY_DEPTH_POSITIVE_DOWN)
            m.finalize()
            assert False, "expected BathymetrySignError"
        except rdb.BathymetrySignError as exc:
            assert exc.code == 15
            msg = str(exc)
            assert "wet cells" in msg or "convention" in msg
    finally:
        m.close()


def test_stage_on_non_pending_model_raises_not_pending(ocean_nml):
    import rdb

    with rdb.Model(ocean_nml) as m:  # one-shot: never entered the pending window
        try:
            m.stage_bathymetry([[1.0] * 6 for _ in range(8)],
                               rdb.BATHY_DEPTH_POSITIVE_DOWN)
            assert False, "expected NotInitialisedError"
        except rdb.NotInitialisedError:
            pass


def test_gebco_style_bathymetry_with_correct_convention_succeeds():
    gebco = [[-4000.0] * 6 for _ in range(8)]
    m = rdb.Model.pending(_inject_nml())
    try:
        m.stage_bathymetry(gebco, rdb.BATHY_HEIGHT_POSITIVE_UP)
        m.finalize()
        b = m.bathymetry[...]
        assert b[0][0] == 4000.0  # sign-normalised to positive-down
    finally:
        m.close()
