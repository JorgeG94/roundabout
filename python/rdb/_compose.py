"""P5/P6 -- the composer. Turns `Model(grid=..., closures=[...], ...)`'s
curated keyword arguments into ONE `Config` plus the small amount of
metadata `Model` needs to decide whether the P2.5 pending/finalize window
is required (an array bathymetry, or any periodic/fold edge).

This is the ONE place cross-object compose-time checks live (D2.7/D2.8):
KPP xor EPBL, Bryan-Lewis xor Henyey, Henyey needs a non-cartesian grid,
MEKE backscatter needs a non-zero biharmonic backstop,
`ImplicitVerticalFriction(bbl_glue=True)`'s three-prerequisite chain
(two knobs on itself, one on a `LinearDrag` object elsewhere in the same
`closures=[...]` list) -- plus a few smaller ones (SphericalCoriolis /
EquilibriumTide / HenyeyBackground need a non-cartesian grid;
AtmosphericPressure needs SurfaceFluxComponents; BaroclinicJetIC needs
nz=2; a grid's flat-bottom `extent=(.., Lz)` conflicts with an explicit
`bathymetry=`). None of these re-implement `validate_config`'s cross-knob
rules (D5.5) -- every one of them is ALSO enforced independently by the
Fortran; this layer just names the problem one call earlier, with both
objects in view.
"""

from __future__ import annotations

import warnings

from . import closures as C
from . import coriolis as CO
from . import forcing as F
from . import grids as G
from ._config import Config
from ._errors import (ConfigConflictError, NoBoundaryLayerWarning,
                       RdbUnsupportedError)
from ._ffi import BATHY_DEPTH_POSITIVE_DOWN
from ._knob import MISSING


def _is_cartesian(grid):
    return isinstance(grid, G.RectilinearGrid)


class Composed:
    """Result of :func:`compose`: the assembled `Config` plus what
    `Model._create_from_composed` needs to decide one-shot vs.
    pending/stage/finalize."""

    __slots__ = ("config", "needs_pending", "bathymetry_stage",
                 "periodic_x", "periodic_y", "flux_components_enabled")

    def __init__(self, config, needs_pending, bathymetry_stage,
                 periodic_x, periodic_y, flux_components_enabled):
        self.config = config
        self.needs_pending = needs_pending
        self.bathymetry_stage = bathymetry_stage
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        self.flux_components_enabled = flux_components_enabled


def _compose_closures(closures, config, grid):
    # KPP(enabled=False) is the curated OFF-switch (mirrors
    # Diagnostics(enabled=False)) -- it must NOT trip the KPP-xor-EPBL
    # check below, since "KPP(enabled=False)" alongside an EPBL() is a
    # legal (if redundant -- EPBL.apply() already turns use_kpp off
    # itself) way to compose them explicitly.
    kpp_objs = [c for c in closures if isinstance(c, C.KPP)]
    kpp = [c for c in kpp_objs if c.enabled]
    epbl = [c for c in closures if isinstance(c, C.EPBL)]
    if kpp and epbl:
        raise ConfigConflictError(
            "KPP and EPBL are mutually exclusive boundary-layer closures "
            "(&ocean_vmix_nml use_kpp vs &ocean_epbl_nml enable) -- "
            "remove one of KPP(...) / EPBL(...) from closures=[...]")
    if kpp_objs and not kpp and not epbl:
        # KPP(enabled=False) given explicitly, and no EPBL(...) either --
        # NEITHER boundary-layer scheme is active. See
        # NoBoundaryLayerWarning's own docstring for why this warns
        # rather than raises (legitimate for a resting/no-wind
        # dyn-core-only test).
        warnings.warn(
            "closures=[...] has KPP(enabled=False) and no EPBL(...) -- "
            "no surface boundary-layer scheme is active at all (PP81 "
            "interior mixing + background floors + convective "
            "adjustment still run). Fine for a resting/no-wind "
            "dyn-core-only test; a wind-forced run wants EPBL(...) or "
            "KPP() (omit KPP(enabled=False) to keep the schema's own "
            "KPP-on-by-default).", NoBoundaryLayerWarning, stacklevel=3)

    bl = [c for c in closures if isinstance(c, C.BryanLewisBackground)]
    hy = [c for c in closures if isinstance(c, C.HenyeyBackground)]
    if bl and hy:
        raise ConfigConflictError(
            "BryanLewisBackground and HenyeyBackground are mutually "
            "exclusive background-diffusivity profiles (MOM6's "
            "one-background-scheme rule, &ocean_vmix_nml bkgnd_profile "
            "vs bkgnd_henyey) -- remove one from closures=[...]")
    if hy and grid is not None and _is_cartesian(grid):
        raise RdbUnsupportedError(
            "HenyeyBackground requires a non-cartesian grid: on a "
            "cartesian grid geolatT == 0 makes every column equatorial, "
            "and the Fortran fails loud. Use a LatitudeLongitudeGrid, "
            "TripolarGrid, or SupergridGrid instead.")

    meke = [c for c in closures if isinstance(c, C.MEKE)]
    if meke:
        gm = [c for c in closures if isinstance(c, C.GentMcWilliams)]
        if not gm:
            raise ConfigConflictError(
                "MEKE requires GentMcWilliams(...) in the same "
                "closures=[...] list (the Fortran requires "
                "&ocean_gm_nml enable=.true.)")
        if any(m.backscatter for m in meke):
            lateral = [c for c in closures if hasattr(c, "biharmonic_coeff")]
            coeff = max((c.biharmonic_coeff() for c in lateral), default=0.0)
            if coeff <= 0.0:
                raise ConfigConflictError(
                    "MEKE(backscatter=True) requires a biharmonic "
                    "backstop with a NON-ZERO dissipation coefficient in "
                    "the same closures=[...] list -- "
                    "Smagorinsky(biharmonic=True, C_biharmonic=...>0), "
                    "Leith(biharmonic=True, C_biharmonic=...>0), or "
                    "ConstantViscosity(nu_4=...>0). This is not invented "
                    "for you: picking a dissipation scale is guessing at "
                    "physics.")

    ivf = [c for c in closures if isinstance(c, C.ImplicitVerticalFriction)]
    for f in ivf:
        if not f.bbl_glue:
            continue
        missing = []
        if not f.harmonic_thickness:
            missing.append("harmonic_thickness=True (-> hvel_mom6)")
        if not f.drag:
            missing.append("drag=True (-> implicit_drag)")
        if not any(isinstance(c, C.LinearDrag) for c in closures):
            missing.append("a LinearDrag(...) bottom-drag object in the "
                            "same closures=[...] list (bbl_glue requires "
                            "a LINEAR bottom drag)")
        if missing:
            raise ConfigConflictError(
                "ImplicitVerticalFriction(bbl_glue=True) needs THREE "
                "prerequisites (harmonic_thickness + drag + a linear "
                "bottom drag); missing: " + "; ".join(missing))

    for c in closures:
        c.apply(config)


def compose(*, grid=None, bathymetry=None, timestep=None, coriolis=None,
            momentum_advection=None, buoyancy=None, pressure=None,
            vcoord=None, closures=(), boundaries=None, forcing=(),
            tracers=(), initial_condition=None, barotropic=None,
            continuity=None, diagnostics=None, restart=None,
            duration=None, config=None):
    cfg = config if config is not None else Config()

    periodic_x = periodic_y = False
    forced_edges = {}
    if grid is not None:
        grid.apply(cfg)
        periodic_x = getattr(grid, "periodic_x", False)
        periodic_y = getattr(grid, "periodic_y", False)
        forced_edges = getattr(grid, "forced_edges", {}) or {}
        for edge, value in forced_edges.items():
            setattr(cfg.ocean_bc, edge, value)

    if timestep is not None:
        cfg.time.dt_fixed = float(timestep)

    if coriolis is not None:
        if (isinstance(coriolis, CO.SphericalCoriolis) and grid is not None
                and _is_cartesian(grid)):
            raise RdbUnsupportedError(
                "SphericalCoriolis requires a non-cartesian grid: f is "
                "filled from geolat, which a RectilinearGrid does not "
                "have. Use BetaPlane/FPlane on a cartesian grid instead.")
        coriolis.apply(cfg)

    if momentum_advection is not None:
        momentum_advection.apply(cfg)

    if buoyancy is not None:
        buoyancy.apply(cfg)

    if pressure is not None:
        pressure.apply(cfg)

    if vcoord is not None:
        vcoord.apply(cfg)

    _compose_closures(list(closures), cfg, grid)

    if boundaries is not None:
        boundaries.apply(cfg, periodic_x=periodic_x, periodic_y=periodic_y,
                          forced_edges=forced_edges)
    else:
        # A grid's topology MUST reach &ocean_bc_nml even when the user
        # supplies no Boundaries() object.  Without this, a grid built with
        # topology=(PERIODIC, BOUNDED) derives nghost correctly, sets
        # .periodic_x -- and then silently leaves west/east at their "wall"
        # default, turning a re-entrant channel into a CLOSED BASIN that runs,
        # conserves mass, and is wrong.  TripolarGrid escaped this only
        # because it also sets .forced_edges.
        if periodic_x:
            cfg.ocean_bc.west = "periodic"
            cfg.ocean_bc.east = "periodic"
        if periodic_y:
            cfg.ocean_bc.south = "periodic"
            cfg.ocean_bc.north = "periodic"

    flux_components_enabled = any(
        isinstance(f, F.SurfaceFluxComponents) for f in forcing)
    has_atmp = any(isinstance(f, F.AtmosphericPressure) for f in forcing)
    if has_atmp and not flux_components_enabled:
        raise ConfigConflictError(
            "AtmosphericPressure requires SurfaceFluxComponents() in the "
            "same forcing=[...] list (the Fortran only reads p_surf_const "
            "when the component-assembly path is enabled)")
    for f in forcing:
        if (isinstance(f, F.EquilibriumTide) and grid is not None
                and _is_cartesian(grid)):
            raise RdbUnsupportedError(
                "EquilibriumTide requires a non-cartesian grid (the "
                "equilibrium-tide generator needs geolat/geolon)")
        f.apply(cfg)

    ic_list = ([] if initial_condition is None else
               (list(initial_condition)
                if isinstance(initial_condition, (list, tuple))
                else [initial_condition]))
    for ic in ic_list:
        if (isinstance(ic, F.BaroclinicJetIC) and grid is not None
                and grid.nz != 2):
            raise ValueError(
                f"BaroclinicJetIC requires nz=2, got nz={grid.nz}")
        ic.apply(cfg)

    for t in tracers:
        t.apply(cfg)

    if barotropic is not None:
        barotropic.apply(cfg)
    if continuity is not None:
        continuity.apply(cfg)
    if diagnostics is not None:
        diagnostics.apply(cfg)
    if restart is not None:
        restart.apply(cfg)
    if duration is not None:
        duration.apply(cfg)

    bathymetry_stage = None
    if isinstance(bathymetry, G.Bathymetry):
        if grid is not None and (bathymetry.nx, bathymetry.ny) != (
                grid.nx, grid.ny):
            raise ValueError(
                f"bathymetry shape ({bathymetry.nx}x{bathymetry.ny}) does "
                f"not match grid size ({grid.nx}x{grid.ny})")
        bathymetry_stage = (bathymetry.depth, BATHY_DEPTH_POSITIVE_DOWN)
    elif bathymetry is not None:
        flat_z = getattr(grid, "flat_bottom_depth", None)
        if flat_z is not None:
            raise ConfigConflictError(
                "grid=...(extent=(Lx, Ly, Lz)) gave a flat-bottom "
                "reference depth AND an explicit bathymetry= was also "
                "given -- these are mutually exclusive (D2.1); drop one.")
        bathymetry.apply(cfg)
    else:
        flat_z = getattr(grid, "flat_bottom_depth", None)
        if flat_z is not None:
            G.FlatBottom(flat_z).apply(cfg)

    if grid is not None:
        pv_scheme = "centered"
        if momentum_advection is not None and hasattr(
                momentum_advection, "pv_adv_scheme"):
            pv_scheme = momentum_advection.pv_adv_scheme
        recon = cfg.ocean_vmix.tracer_recon
        tracer_recon = "ppm" if recon is MISSING else recon
        kappa_vertex = any(
            isinstance(c, C.KappaShear) and c.at_vertex for c in closures)
        periodic_any = periodic_x or periodic_y
        tripolar_fold = getattr(grid, "tripolar_fold", False)
        minimum = G.required_halo(
            periodic=periodic_any, pv_advection=pv_scheme,
            tracer_recon=tracer_recon, north_fold=tripolar_fold,
            decomposed=False, kappa_shear_at_vertex=kappa_vertex)
        requested = getattr(grid, "requested_halo", None)
        if requested is not None:
            if requested < minimum:
                raise RdbUnsupportedError(
                    f"halo={requested} is short for this grid + numerics "
                    f"(periodic={periodic_any}, pv_advection={pv_scheme!r}, "
                    f"tracer_recon={tracer_recon!r}, "
                    f"tripolar_fold={tripolar_fold}, "
                    f"kappa_shear_at_vertex={kappa_vertex}); needs at "
                    f"least {minimum}. Roundabout cannot inflate-and-rebuild "
                    f"(D2.2) -- raise halo= or drop the binding option.")
            nghost = requested
        else:
            nghost = minimum
        cfg.grid.nghost = nghost

    needs_pending = bool(
        periodic_x or periodic_y or forced_edges
        or bathymetry_stage is not None)
    return Composed(cfg, needs_pending, bathymetry_stage, periodic_x,
                     periodic_y, flux_components_enabled)
