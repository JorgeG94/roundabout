"""Tracer-by-name access and a broader sweep of the error-ring mapping."""

import rdb
import pytest


def test_tracer_names_include_salinity_and_temperature(ocean_toml):
    with rdb.Model(ocean_toml) as m:
        names = m.tracer_names
        assert "salinity" in names
        assert "temperature" in names


def test_tracer_not_found_raises_typed_exception(ocean_toml):
    with rdb.Model(ocean_toml) as m:
        try:
            m.tracer("definitely_not_a_tracer")
            assert False, "expected TracerNotFoundError"
        except rdb.TracerNotFoundError as exc:
            assert exc.code == 13
            assert "definitely_not_a_tracer" in str(exc)


def test_kinetic_energy_and_total_mass_are_finite(ocean_toml):
    with rdb.Model(ocean_toml) as m:
        m.step(2)
        ke = m.kinetic_energy
        mass = m.total_mass
        assert ke == ke  # not NaN
        assert mass > 0.0


def test_numpy_array_interface_interop(ocean_toml):
    np = __import__("pytest").importorskip("numpy")
    with rdb.Model(ocean_toml) as m:
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
    cfg_text = (
        """[sim]
sim_type = "ocean"
[grid]
nx = 8
ny = 6
dx = 1.0
dy = 1.0
nghost = 3
[ocean.grid]
grid_config = "spherical"
lon_west = 0.0
lat_south = 10.0
rad_earth = 6371000.0
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
"""
    )
    r, d = 6371000.0, math.radians(1.0)
    area = sum(8 * r * math.cos(math.radians(10.0 + j + 0.5)) * d * r * d
               for j in range(6))
    with rdb.Model(cfg_text) as m:
        mass = m.total_mass
        assert abs(mass - 1025.0 * 200.0 * area) <= 1e-10 * mass
        assert m.kinetic_energy == 0.0


def test_toml_parse_errors_carry_the_front_end_line_through_the_error_ring():
    # The strict TOML front end names `<source>:<line>`; the C API's error
    # ring must hand that text to Python unchanged.
    with pytest.raises(rdb.ConfigParseError,
                       match=r"<config-string>:4: unknown key "
                             r"'totally_bogus_key' in table \[grid\]"):
        rdb.Model("[grid]\nnx = 8\n\ntotally_bogus_key = 1\n")
    with pytest.raises(rdb.ConfigParseError, match="TOML syntax error"):
        rdb.Model("[grid\nnx = 8\n")
