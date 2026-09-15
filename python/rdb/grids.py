"""P5 D2.1/D2.2/D2.3 -- topology, grids, and the immersed bottom (bathymetry).

Every class here is a WRITER over the P4 generated `Config` (D2.0 rule 1):
its only effect is `.apply(config)` setting named knobs, or -- for the one
array-shaped object, `Bathymetry` -- carrying a normalised interior array
plus enough metadata for `rdb._compose` to stage it through the P2.5
pending/finalize window. See
`tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md` D2.1-D2.3
and D6.2/D6.4/D6.5 for the full rationale; this module states the mapping
and the guards, not the essay.

Stdlib only (no numpy at import time; `Bathymetry` accepts a numpy array at
runtime through `_ffi.flatten_nested`'s `__array_interface__` duck-typing
without ever importing numpy itself).
"""

from __future__ import annotations

import enum
import warnings
from pathlib import Path

from . import _ffi
from ._errors import (BathymetrySignError, RdbIOError,
                       RdbUnsupportedError)

#: Roundabout's own LAND_DEPTH_THRESHOLD (rdb_constants.F90:71) -- the Fortran
#: constant this class's `land_threshold` default mirrors.
LAND_DEPTH_THRESHOLD = 2.0


def _reject_z(z):
    if z is not None:
        raise RdbUnsupportedError(
            "z= is not supported: Roundabout has no explicit vertical-grid "
            "input in any form (no interface list) -- only nz_layers + "
            "vcoord_type + scalar family parameters. Use nz= (on the grid's "
            "`size`) and vcoord= (on Model(...)) instead. See D6.4.")


class Topology(enum.Enum):
    """Per-dimension domain topology (mirrors Oceananigans Periodic/Bounded).

    PERIODIC wraps BOTH edges by construction (the namelist's per-edge
    'periodic' spelling collapses to this single per-dimension flag,
    D2.1). BOUNDED has two physical ends; what happens there is a
    `boundaries=` condition, set separately on `Model(...)`.

    There is no FLAT -- Oceananigans' Flat is compile-time dispatch to
    zero/identity operators; Roundabout's ny=1 still runs live differencing.
    The string "flat" (any case) in a topology tuple raises
    RdbUnsupportedError (D6.5).
    """
    PERIODIC = "periodic"
    BOUNDED = "bounded"


def _check_topology(topology):
    if len(topology) != 2:
        raise ValueError("topology must be a 2-tuple (x, y)")
    out = []
    for t in topology:
        if isinstance(t, str) and t.strip().lower() == "flat":
            raise RdbUnsupportedError(
                "Topology.FLAT has no analogue: Oceananigans' Flat is "
                "compile-time dispatch to zero/identity operators; "
                "Roundabout's ny=1 (or nx=1) still runs live differencing "
                "kernels. Use a size-1 dimension with BOUNDED, not Flat. "
                "See D6.5.")
        if isinstance(t, Topology):
            out.append(t)
        else:
            raise TypeError(f"topology entries must be Topology members, "
                             f"got {t!r}")
    return tuple(out)


def required_halo(*, periodic=False, pv_advection="centered",
                   tracer_recon="ppm", north_fold=False, decomposed=False,
                   kappa_shear_at_vertex=False) -> int:
    """Minimum `nghost` for this grid + numerics -- max() over every rule,
    queried directly from the Fortran (`rdb_ocean_required_halo`, a
    stateless ABI call, no handle needed) so this NEVER re-implements the
    six scattered minimum-nghost rules in Python, where they would drift
    (D2.2). `halo=None` on a grid derives through this; `halo=` given and
    short raises, naming the binding rule.
    """
    lib = _ffi.get_lib()
    pv = pv_advection.encode("utf-8")
    tr = tracer_recon.encode("utf-8")
    return lib.rdb_ocean_required_halo(
        pv, len(pv), tr, len(tr), int(bool(periodic)), int(bool(north_fold)),
        int(bool(decomposed)), int(bool(kappa_shear_at_vertex)))


class _GridBase:
    """Shared grid-object contract the composer (`rdb._compose`) reads:
    `.periodic_x` / `.periodic_y` (bool), `.requested_halo` (user's `halo=`
    or None), `.tripolar_fold` (bool), `.forced_edges` (dict of edge name
    -> `&ocean_bc_nml` value the grid itself owns, e.g. tripolar's fold),
    `.nx` / `.ny` / `.nz`.
    """
    tripolar_fold = False
    forced_edges: dict = {}

    def apply(self, config):  # pragma: no cover - overridden
        raise NotImplementedError


class RectilinearGrid(_GridBase):
    """Uniform Cartesian grid. -> `&ocean_grid_nml grid_config="cartesian"`.

    ``RectilinearGrid(*, size, extent=None, x=None, y=None, halo=None,
    topology=(Topology.BOUNDED, Topology.BOUNDED), z=None,
    axis_units="meters", radius=None)``

    ``size=(nx, ny, nz)`` -- nz is a LAYER COUNT (-> `&nonhydrostatic_nml
    nz_layers`), not an interface list (D6.4). Give exactly one of
    ``extent=(Lx, Ly[, Lz])`` or ``x=(lo, hi), y=(lo, hi)`` -- in
    METRES under the default ``axis_units="meters"``, in DEGREES under
    ``axis_units="degrees"`` (see below; a non-zero ``x=``/``y=`` origin
    is recorded for `repr()` only -- Roundabout's cartesian metric has no
    origin). Both route to ``&ocean_grid_nml len_lon=Lx, len_lat=Ly``,
    from which the Fortran itself derives ``dx``/``dy`` -- this object
    never computes them, so there is one source of truth for the
    derivation. ``Lz`` (from a 3-element ``extent``), if given, is
    recorded as a flat-bottom reference depth for `Model` to cross-check
    against an explicit ``bathymetry=``.

    ``axis_units`` mirrors MOM6 ``AXIS_UNITS``: ``"meters"`` (default),
    ``"degrees"`` or ``"km"``. With ``"degrees"`` the extent is read as
    degrees and ``dx = rad_earth * len_lon * pi/180 / nx`` (arc length, no
    cos(lat)) -- this is the MOM6-parity mode the double_gyre reference
    uses. ``radius=`` (Earth radius, metres) is the matching
    ``&ocean_grid_nml rad_earth`` -- meaningless under ``"meters"``/
    ``"km"`` (nothing reads it) so it is left unset (Fortran default)
    unless given explicitly -- give it under ``"degrees"`` rather than
    relying on the Fortran's own default, which silently changes the
    derived ``dx``/``dy`` if it is not the radius you had in mind.

    Default topology is ``(BOUNDED, BOUNDED)``, deliberately NOT
    Oceananigans' periodic-by-default: a doubly-periodic ocean basin is
    rarely intended, and PERIODIC additionally forces ``nghost >= 3``.
    """

    def __init__(self, *, size, extent=None, x=None, y=None, halo=None,
                 topology=(Topology.BOUNDED, Topology.BOUNDED), z=None,
                 axis_units="meters", radius=None):
        _reject_z(z)
        if axis_units not in ("meters", "degrees", "km"):
            raise ValueError(
                f"axis_units must be one of 'meters', 'degrees', 'km' "
                f"(mirroring MOM6 AXIS_UNITS); got {axis_units!r}")
        self.axis_units = axis_units
        self.radius = radius
        if len(size) != 3:
            raise ValueError("size must be (nx, ny, nz)")
        self.nx, self.ny, self.nz = size
        self.topology = _check_topology(topology)
        self.periodic_x = self.topology[0] is Topology.PERIODIC
        self.periodic_y = self.topology[1] is Topology.PERIODIC
        self.requested_halo = halo
        self.origin = None
        self.flat_bottom_depth = None

        if extent is not None and (x is not None or y is not None):
            raise ValueError("give either extent= or x=/y=, not both")
        if extent is not None:
            if len(extent) == 3:
                Lx, Ly, Lz = extent
                self.flat_bottom_depth = Lz
            else:
                Lx, Ly = extent
        elif x is not None and y is not None:
            Lx = x[1] - x[0]
            Ly = y[1] - y[0]
            self.origin = (x[0], y[0])
        else:
            raise ValueError("give extent=(Lx, Ly[, Lz]) or x=(lo, hi), "
                              "y=(lo, hi)")
        self.len_lon, self.len_lat = Lx, Ly

    def apply(self, config):
        config.ocean_grid.grid_config = "cartesian"
        config.ocean_grid.axis_units = self.axis_units
        config.ocean_grid.len_lon = float(self.len_lon)
        config.ocean_grid.len_lat = float(self.len_lat)
        if self.radius is not None:
            config.ocean_grid.rad_earth = float(self.radius)
        config.grid.nx = self.nx
        config.grid.ny = self.ny
        config.nonhydrostatic.nz_layers = self.nz


class LatitudeLongitudeGrid(_GridBase):
    """Lon-lat sector. -> `grid_config="spherical"`.

    ``LatitudeLongitudeGrid(*, size, longitude, latitude, radius=6.378e6,
    halo=None, topology=None, z=None)``

    ``longitude``/``latitude`` are ``(lo, hi)`` DEGREES, increasing.
    Roundabout stores an SW corner + uniform spacing, so ``dx``/``dy`` are
    computed here in DEGREES (`rdb_ocean_setup.F90:144` -- the single most
    confusing thing in the raw namelist) and written directly to
    ``&grid_nml``.

    ``topology=None`` derives Oceananigans-style: x is PERIODIC iff the
    longitude extent is exactly 360 degrees, else BOUNDED; y is always
    BOUNDED (PERIODIC in latitude raises).
    """

    def __init__(self, *, size, longitude, latitude, radius=6.378e6,
                 halo=None, topology=None, z=None):
        _reject_z(z)
        if len(size) != 3:
            raise ValueError("size must be (nx, ny, nz)")
        self.nx, self.ny, self.nz = size
        self.longitude, self.latitude = tuple(longitude), tuple(latitude)
        if not (-90.0 <= self.latitude[0] <= 90.0
                and -90.0 <= self.latitude[1] <= 90.0):
            raise ValueError("latitude must be within [-90, 90]")
        self.radius = radius
        self.requested_halo = halo

        if topology is None:
            x_periodic = abs((self.longitude[1] - self.longitude[0])
                              - 360.0) < 1e-9
            self.topology = (Topology.PERIODIC if x_periodic
                              else Topology.BOUNDED, Topology.BOUNDED)
        else:
            self.topology = _check_topology(topology)
            if self.topology[1] is Topology.PERIODIC:
                raise RdbUnsupportedError(
                    "PERIODIC latitude is not representable on "
                    "LatitudeLongitudeGrid -- the pole is a physical "
                    "boundary, not a wrap. Use TripolarGrid for a global "
                    "grid with an Arctic cap.")
        self.periodic_x = self.topology[0] is Topology.PERIODIC
        self.periodic_y = False

        self.dx = (self.longitude[1] - self.longitude[0]) / self.nx
        self.dy = (self.latitude[1] - self.latitude[0]) / self.ny

    def apply(self, config):
        config.ocean_grid.grid_config = "spherical"
        config.ocean_grid.lon_west = float(self.longitude[0])
        config.ocean_grid.lat_south = float(self.latitude[0])
        config.ocean_grid.rad_earth = float(self.radius)
        config.grid.nx = self.nx
        config.grid.ny = self.ny
        config.grid.dx = float(self.dx)
        config.grid.dy = float(self.dy)
        config.nonhydrostatic.nz_layers = self.nz


class TripolarGrid(_GridBase):
    """Murray (1996) tripolar global grid. -> `grid_config="tripolar"`.

    ``TripolarGrid(*, size, latitude_south=-80.0, phi_join=65.0,
    lon_pole=100.0, radius=6.378e6, halo=None, z=None)``

    Deliberately NO `topology` argument: the topology is
    ``(PERIODIC, north-folded)`` by construction, exactly as
    Oceananigans hardcodes ``(Periodic, fold, TZ)``. `halo` defaults to 3
    (the fold needs `nghost >= 3`). `phi_join` must exceed
    `latitude_south`. Single-rank only.
    """

    def __init__(self, *, size, latitude_south=-80.0, phi_join=65.0,
                 lon_pole=100.0, radius=6.378e6, halo=None, z=None):
        _reject_z(z)
        if len(size) != 3:
            raise ValueError("size must be (nx, ny, nz)")
        if phi_join <= latitude_south:
            raise ValueError(
                f"phi_join ({phi_join}) must exceed latitude_south "
                f"({latitude_south})")
        self.nx, self.ny, self.nz = size
        self.latitude_south = latitude_south
        self.phi_join = phi_join
        self.lon_pole = lon_pole
        self.radius = radius
        self.requested_halo = 3 if halo is None else halo
        self.periodic_x = True
        self.periodic_y = False
        self.tripolar_fold = True
        self.forced_edges = {"west": "periodic", "east": "periodic",
                              "north": "tripolar_fold"}
        self.dx = 360.0 / self.nx
        self.dy = (90.0 - self.latitude_south) / self.ny

    def apply(self, config):
        config.ocean_grid.grid_config = "tripolar"
        config.ocean_grid.lat_south = float(self.latitude_south)
        config.ocean_grid.phi_join = float(self.phi_join)
        config.ocean_grid.lon_pole = float(self.lon_pole)
        config.ocean_grid.rad_earth = float(self.radius)
        config.grid.nx = self.nx
        config.grid.ny = self.ny
        config.grid.dx = float(self.dx)
        config.grid.dy = float(self.dy)
        config.nonhydrostatic.nz_layers = self.nz


class SupergridGrid(_GridBase):
    """MOM6 mosaic-backed grid. -> `grid_config="supergrid"`.

    ``SupergridGrid(*, size, file, halo=None, z=None)``

    The in-memory array form (`x, y, dx, dy, area` arrays, the ONLY route
    to a stretched horizontal mesh) needs
    `metrics_assemble_from_supergrid_arrays`
    (`rdb_ocean_metrics.F90:1061`) exported past `private` -- a one-line
    Fortran change out of scope for this (pure-Python) phase, so it raises
    here. Only `file=` works today.
    """

    def __init__(self, *, size, file=None, x=None, y=None, dx=None, dy=None,
                 area=None, halo=None, z=None):
        _reject_z(z)
        if any(v is not None for v in (x, y, dx, dy, area)):
            raise RdbUnsupportedError(
                "SupergridGrid's in-memory array form needs "
                "metrics_assemble_from_supergrid_arrays "
                "(rdb_ocean_metrics.F90:1061) exported past `private` -- a "
                "one-line Fortran change not made in this (pure-Python) "
                "phase. Use file= (a MOM6 mosaic NetCDF) instead.")
        if file is None:
            raise ValueError("SupergridGrid requires file=")
        p = Path(file)
        if not p.is_file():
            raise RdbIOError(f"supergrid file not found: {file}")
        if len(size) != 3:
            raise ValueError("size must be (nx, ny, nz)")
        self.nx, self.ny, self.nz = size
        self.file = str(p)
        self.requested_halo = halo
        self.periodic_x = False
        self.periodic_y = False

    def apply(self, config):
        config.ocean_grid.grid_config = "supergrid"
        config.ocean_grid.supergrid_file = self.file
        config.grid.nx = self.nx
        config.grid.ny = self.ny
        config.nonhydrostatic.nz_layers = self.nz


# ======================================================================
# Bathymetry -- the immersed bottom
# ======================================================================

_CONVENTIONS = ("depth_positive_down", "height_positive_up")


class Bathymetry:
    """Array bathymetry -- the Oceananigans `GridFittedBottom` analogue.

    ``Bathymetry(depth, *, convention, land_threshold=2.0,
    minimum_depth=None)``

    ``depth``: interior-sized ``(nx, ny)`` array / nested lists / anything
    with ``__array_interface__`` (numpy, without this module ever
    importing numpy). Ghost rows are filled by the library at stage time,
    never by the caller.

    ``convention`` is REQUIRED, no default -- the single most dangerous
    argument in this API (D6.2): ``"depth_positive_down"`` is Roundabout's own
    ``%barotropic%b`` (a 4000 m deep cell is +4000.0);
    ``"height_positive_up"`` is GEBCO/ETOPO/Oceananigans (the same cell is
    -4000.0). Getting this backwards makes every cell read as land -- a
    clean, crash-free, entirely wrong quiescent run.

    Normalises to positive-down HERE, in Python, then checks the wet
    fraction of the NORMALISED array: zero wet cells raises
    `BathymetrySignError` (naming the median and the other convention);
    under 1% warns. This is an EARLY, Python-side instance of the same
    check `rdb_ocean_stage_bathymetry` performs in Fortran -- both run;
    neither alone is trusted as the only guard.
    """

    def __init__(self, depth, *, convention, land_threshold=2.0,
                 minimum_depth=None):
        if convention not in _CONVENTIONS:
            raise ValueError(
                f"convention must be one of {_CONVENTIONS!r} (no default "
                f"exists -- see D6.2), got {convention!r}")
        self.convention = convention
        self.land_threshold = land_threshold
        self.minimum_depth = minimum_depth

        if hasattr(depth, "__array_interface__"):
            shape = tuple(depth.__array_interface__["shape"])
        else:
            shape = (len(depth), len(depth[0]))
        self.nx, self.ny = shape

        flat = _ffi.flatten_nested(depth, shape)
        if convention == "height_positive_up":
            flat = [-v for v in flat]

        wet = [v > land_threshold for v in flat]
        n_wet = sum(wet)
        n = len(flat)
        if minimum_depth is not None:
            flat = [max(v, minimum_depth) if w else v
                    for v, w in zip(flat, wet)]

        wet_fraction = n_wet / n if n else 0.0
        if wet_fraction == 0.0:
            sorted_vals = sorted(flat)
            median = sorted_vals[n // 2] if n else float("nan")
            other = ("height_positive_up" if convention ==
                     "depth_positive_down" else "depth_positive_down")
            raise BathymetrySignError(
                f"Bathymetry: normalised to 0 wet cells under "
                f"convention={convention!r} (median depth {median!r} m, "
                f"land_threshold={land_threshold} m). This is almost "
                f"always the sign convention read backwards -- try "
                f"convention={other!r}. See D6.2.")
        if wet_fraction < 0.01:
            warnings.warn(
                f"Bathymetry: only {wet_fraction * 100:.3f}% wet cells "
                f"under convention={convention!r} -- double-check the "
                f"sign convention.", stacklevel=2)

        wet_vals = [v for v, w in zip(flat, wet) if w]
        self._min = min(wet_vals)
        self._max = max(wet_vals)
        self.wet_fraction = wet_fraction
        self.depth = _ffi.nest_flat(flat, shape)  # interior, positive-down

    def __repr__(self):
        return (f"Bathymetry({self.nx}x{self.ny}, depth "
                f"{self._min:.1f}-{self._max:.1f} m, "
                f"{self.wet_fraction * 100:.1f}% wet, from "
                f"{self.convention})")


class FlatBottom:
    """``FlatBottom(depth)`` -> `topo_config="flat"`, `max_depth=depth`
    (POSITIVE)."""

    def __init__(self, depth):
        self.depth = depth

    def apply(self, config):
        config.ocean_topo.topo_config = "flat"
        config.ocean_topo.max_depth = float(self.depth)


class Spoon:
    """``Spoon(max_depth=2000.0, edge_depth=100.0, slope_scale=4.0e5)``
    -> `topo_config="spoon"`.

    `slope_scale` is in METRES; the Fortran converts to grid units at
    dispatch (`topo_length_to_grid_units`), so this object never needs to
    know the grid's units (the spherical-seamount bug this closes).
    """

    def __init__(self, max_depth=2000.0, edge_depth=100.0,
                 slope_scale=4.0e5):
        self.max_depth, self.edge_depth, self.slope_scale = (
            max_depth, edge_depth, slope_scale)

    def apply(self, config):
        config.ocean_topo.topo_config = "spoon"
        config.ocean_topo.max_depth = float(self.max_depth)
        config.ocean_topo.edge_depth = float(self.edge_depth)
        config.ocean_topo.slope_scale = float(self.slope_scale)


class Seamount:
    """``Seamount(max_depth=2000.0, edge_depth=100.0, slope_scale=4.0e5)``
    -> `topo_config="seamount"`. Same metres-vs-grid-units note as
    `Spoon`."""

    def __init__(self, max_depth=2000.0, edge_depth=100.0,
                 slope_scale=4.0e5):
        self.max_depth, self.edge_depth, self.slope_scale = (
            max_depth, edge_depth, slope_scale)

    def apply(self, config):
        config.ocean_topo.topo_config = "seamount"
        config.ocean_topo.max_depth = float(self.max_depth)
        config.ocean_topo.edge_depth = float(self.edge_depth)
        config.ocean_topo.slope_scale = float(self.slope_scale)


class Island:
    """-> `topo_config="island"`."""

    def apply(self, config):
        config.ocean_topo.topo_config = "island"


class DoubleDrake:
    """-> `topo_config="double_drake"`."""

    def apply(self, config):
        config.ocean_topo.topo_config = "double_drake"


class Neverworld2:
    """``Neverworld2(continent_amp=1.0, roughness_amp=0.05,
    min_depth=500.0)`` -> `topo_config="neverworld2"` + `nl_*`."""

    def __init__(self, continent_amp=1.0, roughness_amp=0.05,
                 min_depth=500.0):
        self.continent_amp, self.roughness_amp, self.min_depth = (
            continent_amp, roughness_amp, min_depth)

    def apply(self, config):
        config.ocean_topo.topo_config = "neverworld2"
        config.ocean_topo.nl_continent_amp = float(self.continent_amp)
        config.ocean_topo.nl_roughness_amp = float(self.roughness_amp)
        config.ocean_topo.nl_min_depth = float(self.min_depth)


class BathymetryFile:
    """``BathymetryFile(path)`` -> `topo_config="file"` +
    `&output_nml bathymetry_file`.

    The knob lives in `&output_nml`, not `&ocean_topo_nml`, despite
    `topo_config="file"` being its only consumer -- cosmetic in a
    namelist, confusing in a typed surface; this object hides the
    seam. Pre-flights the file's existence (the Fortran reader still
    `error stop`s on a missing/malformed file).
    """

    def __init__(self, path):
        p = Path(path)
        if not p.is_file():
            raise RdbIOError(f"bathymetry file not found: {path}")
        self.path = str(p)

    def apply(self, config):
        config.ocean_topo.topo_config = "file"
        config.output.bathymetry_file = self.path
