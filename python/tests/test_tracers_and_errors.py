"""Tracer-by-name access and a broader sweep of the error-ring mapping."""

import rdb


def test_tracer_names_include_salinity_and_temperature(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        names = m.tracer_names
        assert "salinity" in names
        assert "temperature" in names


def test_tracer_not_found_raises_typed_exception(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        try:
            m.tracer("definitely_not_a_tracer")
            assert False, "expected TracerNotFoundError"
        except rdb.TracerNotFoundError as exc:
            assert exc.code == 13
            assert "definitely_not_a_tracer" in str(exc)


def test_kinetic_energy_and_total_mass_are_finite(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        m.step(2)
        ke = m.kinetic_energy
        mass = m.total_mass
        assert ke == ke  # not NaN
        assert mass > 0.0


def test_numpy_array_interface_interop(ocean_nml):
    np = __import__("pytest").importorskip("numpy")
    with rdb.Model(ocean_nml) as m:
        h = m.h
        arr = np.asarray(h)
        assert arr.shape == h.shape
        assert arr.dtype == np.dtype(f"<f{m._wp_bytes}")
        # Read-only: numpy must refuse a write into the mapped view.
        try:
            arr[0, 0, 0] = 1.0
            assert False, "expected numpy to refuse write to readonly array"
        except (ValueError, RuntimeError):
            pass


def test_total_mass_is_area_weighted_on_a_spherical_grid():
    """On a curvilinear grid `grid%dx*grid%dy` is not a cell area (on a
    spherical sector dx/dy are DEGREES), so `total_mass` must weight each
    column by its own `areaT`.  Flat 200 m lat-lon sector at rest:
    mass = rho0 * H * sum_j R^2 cos(lat_j) dlon dlat, rho0 = 1025 (the
    API's documented reference density)."""
    import math
    nml = (
        '&sim_nml sim_type = "ocean" /\n'
        "&grid_nml nx = 8, ny = 6, dx = 1.0, dy = 1.0, nghost = 3 /\n"
        '&ocean_grid_nml grid_config = "spherical", lon_west = 0.0, '
        'lat_south = 10.0, rad_earth = 6371000.0 /\n'
        "&nonhydrostatic_nml nz_layers = 3 /\n"
        "&time_nml t_end = 86400.0, dt_fixed = 300.0 /\n"
        '&ocean_topo_nml topo_config = "flat", max_depth = 200.0 /\n'
        "&ocean_bt_nml auto_n_inner = .true. /\n"
        "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /\n"
        "&ocean_diag_nml enabled = .false. /\n"
        "&output_nml output_to_file = .false. /\n"
    )
    r, d = 6371000.0, math.radians(1.0)
    area = sum(8 * r * math.cos(math.radians(10.0 + j + 0.5)) * d * r * d
               for j in range(6))
    with rdb.Model(nml) as m:
        mass = m.total_mass
        assert abs(mass - 1025.0 * 200.0 * area) <= 1e-10 * mass
        assert m.kinetic_energy == 0.0
