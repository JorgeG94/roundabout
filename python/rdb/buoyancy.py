"""P5 D2.5 -- equation of state. `&ocean_eos_nml` + `&ocean_ic_nml`."""

from __future__ import annotations

from ._errors import RdbUnsupportedError


class LinearEquationOfState:
    """``LinearEquationOfState(thermal_expansion=1.7e-4,
    reference_density=1035.0, haline_contraction=None,
    reference_temperature=None, reference_salinity=None)``

    -> `&ocean_eos_nml eos="linear"`, `&ocean_ic_nml alpha_T, beta_S,
    T_ref, S_ref, rho_0` -- the whole linear reference state:

        rho = rho_0 + beta_S*(S - S_ref) - alpha_T*(T - T_ref)

    A REAL TRAP THIS OBJECT CLOSES: there used to be two `alpha_T` knobs
    with the same documented units and different values -- `&tracer_nml
    alpha_T` (coastal-legacy, fed `tr_T%eos_coeff`, read by nothing on
    the ocean path) and `&ocean_ic_nml alpha_T` (the ocean path, feeds
    `eos%alpha_T`). This object writes the ocean one; the legacy quartet
    (`alpha_T`, `beta_S`, `T_ref`, `S_ref` on `&tracer_nml`) is now
    RETIRED and fails loud if you set it.

    A SECOND TRAP: `thermal_expansion`/`haline_contraction` are
    **DIMENSIONAL** (kg/m^3 per degC / per PSU), because the form above
    is the density-ANOMALY one. Protocols usually quote the FRACTIONAL
    coefficients of ``rho = rho_0*(1 - alpha*dT + beta*dS)`` in 1/degC
    and 1/PSU -- multiply those by `reference_density` before passing
    them in. ISOMIP+ (Asay-Davis et al. 2016), for instance::

        LinearEquationOfState(
            thermal_expansion=1027.51 * 3.733e-5,   # 3.8357e-2
            haline_contraction=1027.51 * 7.843e-4,  # 8.0588e-1
            reference_temperature=-1.0,
            reference_salinity=34.2,
            reference_density=1027.51)

    Passing the fractional numbers straight through under-states the
    density response ~1000x -- a plausible-looking, far too weakly
    stratified run.

    The three salinity/reference keywords default to ``None`` = "leave
    the Fortran's own default alone" rather than to a repeated literal,
    so they are written to the namelist only when you actually ask for
    them -- an existing call site that names neither is serialised
    exactly as before.
    """

    def __init__(self, thermal_expansion=1.7e-4, reference_density=1035.0,
                 haline_contraction=None, reference_temperature=None,
                 reference_salinity=None):
        self.thermal_expansion = thermal_expansion
        self.reference_density = reference_density
        self.haline_contraction = haline_contraction
        self.reference_temperature = reference_temperature
        self.reference_salinity = reference_salinity

    def apply(self, config):
        config.ocean_eos.eos = "linear"
        config.ocean_ic.alpha_T = float(self.thermal_expansion)
        config.ocean_ic.rho_0 = float(self.reference_density)
        if self.haline_contraction is not None:
            config.ocean_ic.beta_S = float(self.haline_contraction)
        if self.reference_temperature is not None:
            config.ocean_ic.T_ref = float(self.reference_temperature)
        if self.reference_salinity is not None:
            config.ocean_ic.S_ref = float(self.reference_salinity)


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
