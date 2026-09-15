"""P5 D2.4 -- rotation (the PARAMETER) and momentum/PV discretisation (the
SCHEME): two objects on purpose, because the namelist splits one concept
across `&physics_nml` / `&ocean_topo_nml` / `&ocean_grid_nml` (the
parameter) and `&ocean_coriolis_nml` (the scheme) in a way that reads
backwards. See 06_python_surface_design.md D2.4.
"""

from __future__ import annotations

from ._errors import RdbUnsupportedError


class FPlane:
    """``FPlane(f)`` -> `&physics_nml coriolis_f`, `&ocean_topo_nml
    coriolis_beta=0`."""

    def __init__(self, f):
        self.f = f

    def apply(self, config):
        config.physics.coriolis_f = float(self.f)
        config.ocean_topo.coriolis_beta = 0.0


class BetaPlane:
    """``BetaPlane(f0, beta, y_ref=0.0)``

    -> `&physics_nml coriolis_f=f0`, `&ocean_topo_nml
    coriolis_beta=beta, coriolis_y_ref=y_ref`, `&ocean_grid_nml
    coriolis_scheme="beta_plane"`.
    """

    def __init__(self, f0, beta, y_ref=0.0):
        self.f0, self.beta, self.y_ref = f0, beta, y_ref

    def apply(self, config):
        config.physics.coriolis_f = float(self.f0)
        config.ocean_topo.coriolis_beta = float(self.beta)
        config.ocean_topo.coriolis_y_ref = float(self.y_ref)
        config.ocean_grid.coriolis_scheme = "beta_plane"


class SphericalCoriolis:
    """``SphericalCoriolis(rotation_rate=7.2921e-5)``

    -> `&ocean_grid_nml coriolis_scheme="planetary", omega=rotation_rate`.
    f is filled from `geolat` by the metrics; there is no f knob to set.
    Requires a non-cartesian grid (checked by the composer, which is the
    only place that has both this object and the grid object in view).
    """

    def __init__(self, rotation_rate=7.2921e-5):
        self.rotation_rate = rotation_rate

    def apply(self, config):
        config.ocean_grid.coriolis_scheme = "planetary"
        config.ocean_grid.omega = float(self.rotation_rate)


_FORM_MAP = {"enstrophy": "sadourny", "energy": "sadourny_energy",
             "hk": "sadourny_hk"}


class WENO:
    """``WENO(order=5)`` -> the `"weno{order}"` `pv_adv_scheme` spelling.
    Orders 3, 5, 7 only (the PV-advection path)."""

    def __init__(self, order=5):
        if order not in (3, 5, 7):
            raise ValueError(f"WENO order must be 3, 5 or 7 (pv-advection "
                              f"path), got {order}")
        self.order = order

    @property
    def scheme(self):
        return f"weno{self.order}"


class Sadourny:
    """The PV/Coriolis-advection discretisation. -> `&ocean_coriolis_nml`.

    ``Sadourny(form="enstrophy", pv_advection="centered",
    use_state_fluxes=False, bound_coriolis=False, corner_h="cell_mean")``

    ``form``: ``"enstrophy"`` (default, velocity form) ->
    ``form="sadourny"``; ``"energy"`` (transport form) ->
    ``"sadourny_energy"``; ``"hk"`` (Hollingsworth-Kallen) ->
    ``"sadourny_hk"``.

    ``pv_advection``: ``"centered"`` (default, bit-identical) or a
    :class:`WENO` instance -> ``pv_adv_scheme = "centered" | "weno3" |
    "weno5" | "weno7"``.

    ``bound_coriolis``/``corner_h`` are ``energy``-form-only: passing a
    non-default value with ``form != "energy"`` raises rather than being
    silently inert.
    """

    def __init__(self, form="enstrophy", pv_advection="centered",
                 use_state_fluxes=False, bound_coriolis=False,
                 corner_h="cell_mean"):
        if form not in _FORM_MAP:
            raise ValueError(f"form must be one of {tuple(_FORM_MAP)}, "
                              f"got {form!r}")
        self.form = form
        self.pv_advection = pv_advection
        self.use_state_fluxes = use_state_fluxes
        self.bound_coriolis = bound_coriolis
        self.corner_h = corner_h
        if form != "energy" and (bound_coriolis is not False
                                  or corner_h != "cell_mean"):
            raise RdbUnsupportedError(
                "bound_coriolis / corner_h are energy-form-only "
                "(&ocean_coriolis_nml form='sadourny_energy'); passing a "
                "non-default value with form='enstrophy' or form='hk' "
                "would be silently inert in the Fortran, so this raises "
                "instead of accepting it.")

    @property
    def pv_adv_scheme(self):
        if isinstance(self.pv_advection, WENO):
            return self.pv_advection.scheme
        if self.pv_advection == "centered":
            return "centered"
        raise ValueError(f"pv_advection must be 'centered' or a WENO(...) "
                          f"instance, got {self.pv_advection!r}")

    def apply(self, config):
        config.ocean_coriolis.form = _FORM_MAP[self.form]
        config.ocean_coriolis.pv_adv_scheme = self.pv_adv_scheme
        config.ocean_coriolis.use_state_fluxes = bool(self.use_state_fluxes)
        if self.form == "energy":
            config.ocean_coriolis.bound_coriolis = bool(self.bound_coriolis)
            config.ocean_coriolis.corner_h = self.corner_h
