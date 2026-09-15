"""Library discovery, working_precision, and the create/step/destroy
lifecycle, including the P0/P2/P2.5 hard-constraints from
docs/ocean_python_api_plan.md S2."""

import rdb


def test_library_loads_and_precision_resolves():
    path = rdb.find_library()
    assert ".so" in path  # may be a versioned SONAME, e.g. libfoo.so.0.1.0
    wp = rdb.working_precision()
    assert wp in (4, 8)


def test_create_step_destroy(ocean_nml):
    m = rdb.Model(ocean_nml)
    try:
        assert m.time == 0.0
        assert m.step_count == 0
        info = m.grid_info
        assert info == {"nx": 8, "ny": 6, "nz": 3, "nghost": info["nghost"]}
        m.step(3)
        assert m.step_count == 3
        assert m.time > 0.0
    finally:
        m.close()


def test_context_manager_closes_deterministically(ocean_nml):
    with rdb.Model(ocean_nml) as m:
        m.step(1)
        assert m.step_count == 1
    assert m._closed


def test_destroy_twice_is_clean_noop(ocean_nml):
    m = rdb.Model(ocean_nml)
    m.close()
    m.close()  # must not raise, must not touch freed memory
    m.close()  # a third time, for good measure


def test_second_live_model_raises_already_exists(ocean_nml):
    m1 = rdb.Model(ocean_nml)
    try:
        try:
            rdb.Model(ocean_nml)
            assert False, "expected AlreadyExistsError"
        except rdb.AlreadyExistsError as exc:
            assert exc.code == 11
            assert "already open" in str(exc) or "already exists" in str(exc)
    finally:
        m1.close()

    # Guard lifts after close(): a fresh model can now be created.
    m2 = rdb.Model(ocean_nml)
    m2.close()


def test_bad_config_raises_typed_exception_with_fortran_message(bad_nml):
    try:
        rdb.Model(bad_nml)
        assert False, "expected ConfigParseError"
    except rdb.ConfigParseError as exc:
        assert exc.code == 1
        msg = str(exc)
        # Not a bare "error 1" -- the actual Fortran schema-parser sentence,
        # naming the offending key.
        assert "totally_bogus_key" in msg
        assert msg != "error 1"


def test_step_on_closed_model_raises(ocean_nml):
    m = rdb.Model(ocean_nml)
    m.close()
    try:
        m.step(1)
        assert False, "expected NotInitialisedError"
    except rdb.NotInitialisedError:
        pass
