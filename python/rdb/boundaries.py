"""P5 D2.11 -- boundaries. `&ocean_bc_nml` (`rdb._config_generated.OceanBc`,
schema-generated since P4.5) plus `&ocean_sponge_nml` (`SpongeMap`).

`Periodic()` is NOT here -- it is a grid topology (`rdb.grids.Topology`).
`Boundaries.apply` raises if an edge the grid's topology has made periodic
(or, for `TripolarGrid`, structurally fold-owned) is also given a condition
here, mirroring Oceananigans' "Cannot set $side $bc in a Periodic
direction!". There are no `bottom=`/`top=` slots: `&ocean_bc_nml` has no z
edges (the vertical boundaries are forcing, not `&ocean_bc_nml`).

Several `&ocean_bc_nml` knobs are SHARED across whichever edges use them
(`radiation_scheme`/`orlanski_rx_max`/`orlanski_gamma`, `nudge_tau_in/out`,
`sponge_width`/`sponge_strength`, `flather_form`, `obc_tidal_nodal`) rather
than per-edge -- a Fortran design already in place, not something this
layer invents. `Boundaries.apply` raises `ConfigConflictError` if two edges
ask for different values of the same shared knob (a real silent
last-writer-wins bug the raw namelist has no way to catch).
"""

from __future__ import annotations

from ._errors import ConfigConflictError, RdbUnsupportedError

_EDGES = ("west", "east", "south", "north")


class _SharedTracker:
    """Collects (group, key) -> (value, first_edge) across edges; raises on
    a differing second write."""

    def __init__(self):
        self._seen = {}

    def set(self, config, group_attr, key, value, edge):
        group = getattr(config, group_attr)
        prior = self._seen.get((group_attr, key))
        if prior is not None and prior[0] != value:
            raise ConfigConflictError(
                f"boundaries: edge {edge!r} sets {group_attr}.{key} = "
                f"{value!r}, but edge {prior[1]!r} already set it to "
                f"{prior[0]!r} -- this is a SHARED knob (one value for "
                f"whichever edges use it), not per-edge. Make the two "
                f"edges agree.")
        self._seen[(group_attr, key)] = (value, edge)
        setattr(group, key, value)


class Wall:
    """-> `"wall"`."""

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "wall")


class Flather:
    """``Flather(external_normal_velocity=0.0, form="legacy",
    radiation=None, nudging=None)``

    -> `"open"` + `flather_form` + `<edge>_ext_u`/`<edge>_ext_v`.

    ``form`` is `"legacy"` (v1 Flather, u_ext=0, no interior velocity --
    the Fortran default) or `"full"` (half-characteristic form, Flather
    1976) -> the SHARED `flather_form` knob. (The design sketch showed a
    default of `"anomaly"`, which is not a legal `flather_form` member --
    that string belongs to `radiation_scheme`, i.e. the `radiation=`
    argument below; this implementation uses the Fortran's own
    `flather_form` default, `"legacy"`, to avoid shipping an invalid
    default.)
    """

    def __init__(self, external_normal_velocity=0.0, form="legacy",
                 radiation=None, nudging=None):
        if form not in ("legacy", "full"):
            raise ValueError(f"form must be 'legacy' or 'full', got "
                              f"{form!r}")
        self.external_normal_velocity = external_normal_velocity
        self.form = form
        self.radiation = radiation
        self.nudging = nudging

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "open")
        shared.set(config, "ocean_bc", "flather_form", self.form, edge)
        vel_key = f"{edge}_ext_u" if edge in ("west", "east") else \
            f"{edge}_ext_v"
        setattr(config.ocean_bc, vel_key,
                float(self.external_normal_velocity))
        if self.radiation is not None:
            self.radiation._apply_shared(config, edge, shared)
        if self.nudging is not None:
            self.nudging._apply_shared(config, edge, shared)


class Clamped:
    """``Clamped(eta=None, velocity=None)`` -> `"clamped"` +
    `<edge>_clamped_eta` / `_clamped_u` / `_clamped_v`."""

    def __init__(self, eta=None, velocity=None):
        self.eta, self.velocity = eta, velocity

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "clamped")
        if self.eta is not None:
            setattr(config.ocean_bc, f"{edge}_clamped_eta", float(self.eta))
        if self.velocity is not None:
            vel_key = f"{edge}_clamped_u" if edge in ("west", "east") else \
                f"{edge}_clamped_v"
            setattr(config.ocean_bc, vel_key, float(self.velocity))


class Chapman:
    """-> `"chapman"`."""

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "chapman")


class Sponge:
    """``Sponge(width, strength, relax_tracers=True)`` -> `"sponge"` +
    the SHARED `sponge_width`/`sponge_strength`/`sponge_relax_tracers`."""

    def __init__(self, width, strength, relax_tracers=True):
        self.width, self.strength, self.relax_tracers = (
            width, strength, relax_tracers)

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "sponge")
        shared.set(config, "ocean_bc", "sponge_width", int(self.width), edge)
        shared.set(config, "ocean_bc", "sponge_strength",
                    float(self.strength), edge)
        shared.set(config, "ocean_bc", "sponge_relax_tracers",
                    bool(self.relax_tracers), edge)


class TidalBoundary:
    """``TidalBoundary(constituents=[(amp, phase, omega), ...],
    nodal=False)``

    -> `"tidal"` + `<edge>_n_tidal`/`_tidal_amp`/`_tidal_phase`/
    `_tidal_omega` + the SHARED `obc_tidal_nodal`. At most 8
    constituents per edge (`OBC_MAX_TIDAL_CFG`).
    """

    _MAX = 8

    def __init__(self, constituents, nodal=False):
        if len(constituents) > self._MAX:
            raise ValueError(
                f"at most {self._MAX} tidal constituents per edge, got "
                f"{len(constituents)}")
        self.constituents, self.nodal = constituents, nodal

    def _apply_edge(self, config, edge, shared):
        setattr(config.ocean_bc, edge, "tidal")
        n = len(self.constituents)
        setattr(config.ocean_bc, f"{edge}_n_tidal", n)
        pad = self._MAX - n
        amp = [float(c[0]) for c in self.constituents] + [0.0] * pad
        phase = [float(c[1]) for c in self.constituents] + [0.0] * pad
        omega = [float(c[2]) for c in self.constituents] + [0.0] * pad
        setattr(config.ocean_bc, f"{edge}_tidal_amp", amp)
        setattr(config.ocean_bc, f"{edge}_tidal_phase", phase)
        setattr(config.ocean_bc, f"{edge}_tidal_omega", omega)
        shared.set(config, "ocean_bc", "obc_tidal_nodal", bool(self.nodal),
                   edge)


class Radiation:
    """``Radiation(scheme="orlanski", rx_max=10.0, gamma=1.0)`` MODIFIER,
    passed as `Flather(radiation=...)`. -> the SHARED
    `radiation_scheme`/`orlanski_rx_max`/`orlanski_gamma`."""

    def __init__(self, scheme="orlanski", rx_max=10.0, gamma=1.0):
        if scheme not in ("anomaly", "orlanski"):
            raise ValueError(f"scheme must be 'anomaly' or 'orlanski', "
                              f"got {scheme!r}")
        self.scheme, self.rx_max, self.gamma = scheme, rx_max, gamma

    def _apply_shared(self, config, edge, shared):
        shared.set(config, "ocean_bc", "radiation_scheme", self.scheme, edge)
        shared.set(config, "ocean_bc", "orlanski_rx_max",
                    float(self.rx_max), edge)
        shared.set(config, "ocean_bc", "orlanski_gamma", float(self.gamma),
                   edge)


class Nudging:
    """``Nudging(tau_in=0.0, tau_out=0.0)`` MODIFIER (asymmetric),
    passed as `Flather(nudging=...)`. -> the SHARED `nudge_tau_in`/
    `nudge_tau_out`."""

    def __init__(self, tau_in=0.0, tau_out=0.0):
        self.tau_in, self.tau_out = tau_in, tau_out

    def _apply_shared(self, config, edge, shared):
        shared.set(config, "ocean_bc", "nudge_tau_in", float(self.tau_in),
                   edge)
        shared.set(config, "ocean_bc", "nudge_tau_out", float(self.tau_out),
                   edge)


class Boundaries:
    """``Boundaries(west=None, east=None, south=None, north=None,
    sponge_width=None, sponge_strength=None, mask_wall_velocity=None)``

    An edge left `None` keeps the Fortran default (`"wall"`) unwritten --
    consistent with "absence is the default" (D2.0 rule 2). An edge whose
    grid topology is PERIODIC (or, on `TripolarGrid`, the structurally
    fold-owned `north`) must not be given a condition; doing so raises.

    `mask_wall_velocity` defaults to the Fortran's own default (`False`,
    i.e. off) even though on is the physically correct behaviour (memory
    `wall-vel-mask-bug`) -- flipping the default re-baselines
    `double_gyre` and is a solver PR, not an API decision.
    """

    def __init__(self, west=None, east=None, south=None, north=None,
                 sponge_width=None, sponge_strength=None,
                 mask_wall_velocity=None):
        self.west, self.east, self.south, self.north = west, east, south, \
            north
        self.sponge_width = sponge_width
        self.sponge_strength = sponge_strength
        self.mask_wall_velocity = mask_wall_velocity

    def apply(self, config, *, periodic_x=False, periodic_y=False,
              forced_edges=None):
        forced_edges = forced_edges or {}
        shared = _SharedTracker()
        edges = {"west": self.west, "east": self.east,
                 "south": self.south, "north": self.north}
        for edge, cond in edges.items():
            locked = ((edge in ("west", "east") and periodic_x)
                      or (edge in ("south", "north") and periodic_y)
                      or edge in forced_edges)
            if locked:
                if cond is not None:
                    raise RdbUnsupportedError(
                        f"Boundaries({edge}=...) was given, but the grid "
                        f"topology makes {edge!r} PERIODIC (or, on a "
                        f"TripolarGrid, structurally fold-owned) -- a "
                        f"boundary condition cannot be set on a periodic "
                        f"edge. Remove {edge}= or change the grid "
                        f"topology (Topology.PERIODIC lives on the grid, "
                        f"not on Boundaries -- see D2.11).")
                continue
            if cond is None:
                continue
            cond._apply_edge(config, edge, shared)
        if self.sponge_width is not None:
            shared.set(config, "ocean_bc", "sponge_width",
                       int(self.sponge_width), "Boundaries(sponge_width=)")
        if self.sponge_strength is not None:
            shared.set(config, "ocean_bc", "sponge_strength",
                       float(self.sponge_strength),
                       "Boundaries(sponge_strength=)")
        if self.mask_wall_velocity is not None:
            config.ocean_bc.mask_wall_velocity = bool(
                self.mask_wall_velocity)


class SpongeMap:
    """Map-driven sponge. -> `&ocean_sponge_nml enable=.true.`

    ``SpongeMap(damp_source="band", target_source="ic", relax_uv=True,
    relax_tracers=True, widths=None)``

    ``widths``: optional ``{"west": (width, strength), ...}`` -- any
    subset of edges. ``relax_h`` is NOT IMPLEMENTED in the Fortran v1 and
    is not offered.
    """

    def __init__(self, damp_source="band", target_source="ic",
                 relax_uv=True, relax_tracers=True, widths=None):
        self.damp_source, self.target_source = damp_source, target_source
        self.relax_uv, self.relax_tracers = relax_uv, relax_tracers
        self.widths = widths or {}

    def apply(self, config):
        config.ocean_sponge.enable = True
        config.ocean_sponge.damp_source = self.damp_source
        config.ocean_sponge.target_source = self.target_source
        config.ocean_sponge.relax_uv = bool(self.relax_uv)
        config.ocean_sponge.relax_tracers = bool(self.relax_tracers)
        for edge, (width, strength) in self.widths.items():
            if edge not in _EDGES:
                raise ValueError(f"unknown edge {edge!r}")
            setattr(config.ocean_sponge, f"{edge}_width", int(width))
            setattr(config.ocean_sponge, f"{edge}_strength",
                    float(strength))
