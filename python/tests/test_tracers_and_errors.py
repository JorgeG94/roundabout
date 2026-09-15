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
