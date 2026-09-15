"""P5 D2.12 -- tracers. `&ocean_tracers_nml` + the registry.

Ordering (D6.7): tracer registration must precede `enter_data`, which
closes `registry_locked` (`rdb_multilayer_state.F90:533`). So tracers are
a `Model()` constructor argument, never a post-create method.
"""

from __future__ import annotations

from ._errors import RdbUnsupportedError


class IdealAge:
    """``IdealAge(young_value=0.0, surface_growth_rate=0.0)``
    -> `&ocean_tracers_nml enable_ideal_age` + the two knobs."""

    def __init__(self, young_value=0.0, surface_growth_rate=0.0):
        self.young_value = young_value
        self.surface_growth_rate = surface_growth_rate

    def apply(self, config):
        config.ocean_tracers.enable_ideal_age = True
        config.ocean_tracers.ideal_age_young_val = float(self.young_value)
        config.ocean_tracers.ideal_age_sfc_growth_rate = float(
            self.surface_growth_rate)


class PseudoSalt:
    """``PseudoSalt()`` -> `enable_pseudo_salt`. Shao (2016) verification
    tracer: seeded to S, given S's surface salt flux and KPP nonlocal
    mirror; its deviation from S measures the passive-vs-active
    transport-path error."""

    def apply(self, config):
        config.ocean_tracers.enable_pseudo_salt = True


class TracerBounds:
    """``TracerBounds(salinity=(0.0, 40.0), temperature=(-2.0, 40.0))``
    -> `&tracer_nml S_min/S_max/T_min/T_max` (coastal-legacy group name,
    but these two clamps are read on the ocean path too)."""

    def __init__(self, salinity=None, temperature=None):
        self.salinity, self.temperature = salinity, temperature

    def apply(self, config):
        if self.salinity is not None:
            config.tracer.S_min = float(self.salinity[0])
            config.tracer.S_max = float(self.salinity[1])
        if self.temperature is not None:
            config.tracer.T_min = float(self.temperature[0])
            config.tracer.T_max = float(self.temperature[1])


class PassiveTracer:
    """``PassiveTracer(name, units="", long_name="")``

    -> `multilayer_state_t%register_passive_tracer` at setup.

    **NOT FUNCTIONAL in this phase.** Registering an arbitrary named
    passive tracer needs a pre-create Fortran entry point
    (`register_passive_tracer` reached before `enter_data` closes
    `registry_locked`), and no such entry point is exposed on
    `include/rdb_ocean.h` today -- the current C ABI (P2/P2.5) only
    reaches tracers that are ALREADY registered (`enable_ideal_age` /
    `enable_pseudo_salt`, both namelist-driven booleans; use
    :class:`IdealAge` / :class:`PseudoSalt` for those). This is a real
    Fortran-side gap, not something P5 (pure Python, no Fortran changes)
    can paper over, so this raises rather than silently doing nothing.
    """

    def __init__(self, name, units="", long_name=""):
        self.name, self.units, self.long_name = name, units, long_name
        raise RdbUnsupportedError(
            f"PassiveTracer({name!r}) has no reachable entry point: "
            f"arbitrary passive-tracer registration needs a pre-create "
            f"call to register_passive_tracer "
            f"(rdb_multilayer_state.F90:533's registry_locked closes at "
            f"enter_data) and include/rdb_ocean.h exposes no such "
            f"entry point yet (P2/P2.5 gap). Use IdealAge() or "
            f"PseudoSalt() -- both are namelist-driven and work today -- "
            f"or extend the C ABI with a pre-create tracer-registration "
            f"call.")

    def apply(self, config):  # pragma: no cover - unreachable, __init__ raises
        raise RdbUnsupportedError("PassiveTracer never reaches apply()")
