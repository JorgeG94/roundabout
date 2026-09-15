"""P5 D2.5 -- equation of state. `&ocean_eos_nml` + `&ocean_ic_nml`."""

from __future__ import annotations

from ._errors import RdbUnsupportedError


class LinearEquationOfState:
    """``LinearEquationOfState(thermal_expansion=1.7e-4,
    reference_density=1035.0)``

    -> `&ocean_eos_nml eos="linear"`, `&ocean_ic_nml alpha_T, rho_0`.

    A REAL TRAP THIS OBJECT CLOSES: there are two `alpha_T` knobs with the
    same documented units and different values -- `&tracer_nml alpha_T`
    (coastal-legacy, feeds `tr_T%eos_coeff`) and `&ocean_ic_nml alpha_T`
    (the ocean path, feeds `eos%alpha_T`). This object writes the ocean
    one only.

    `haline_contraction` is NOT offered: on the ocean path there is
    nowhere to put it (`beta_S` exists only on the coastal-legacy
    `&tracer_nml` path). Offering a keyword that does nothing would be
    exactly the "plausible-looking wrong answer" this design exists to
    prevent (memory `hk-validation`).
    """

    def __init__(self, thermal_expansion=1.7e-4, reference_density=1035.0):
        self.thermal_expansion = thermal_expansion
        self.reference_density = reference_density

    def apply(self, config):
        config.ocean_eos.eos = "linear"
        config.ocean_ic.alpha_T = float(self.thermal_expansion)
        config.ocean_ic.rho_0 = float(self.reference_density)


class Wright1997:
    """``Wright1997()`` -> `eos="wright"`. The production nonlinear EOS."""

    def apply(self, config):
        config.ocean_eos.eos = "wright"


class RoquetSPV:
    """``RoquetSPV()`` -> `eos="roquet_spv"`."""

    def apply(self, config):
        config.ocean_eos.eos = "roquet_spv"


class TEOS10:
    """``TEOS10()`` -> `eos="teos10"`.

    RAISES AT CONSTRUCTION, carrying the Fortran's own reason ("TEOS-10
    not yet device-callable", `eos_validate`). Constructing an object
    guaranteed to abort `create()` twenty lines later is worse than an
    immediate, specific refusal.
    """

    def __init__(self):
        raise RdbUnsupportedError(
            "TEOS10 is not offered: the Fortran's own eos_validate refuses "
            "it ('TEOS-10 not yet device-callable'). Constructing this "
            "object would be guaranteed to abort create() later with a "
            "less specific message -- use Wright1997() (production) or "
            "RoquetSPV() instead.")
