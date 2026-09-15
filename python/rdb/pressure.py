"""P5 D2.6 -- pressure gradient force. `&ocean_pgf_nml`."""

from __future__ import annotations

_RECON = {"plm": 1, "ppm": 2}


class Montgomery:
    """-> `form="mont"` (the namelist default)."""

    def apply(self, config):
        config.ocean_pgf.form = "mont"


class FiniteVolumeLite:
    """-> `form="fv_lite"`."""

    def apply(self, config):
        config.ocean_pgf.form = "fv_lite"


class FiniteVolumeWright:
    """-> `form="fv_wright"`."""

    def apply(self, config):
        config.ocean_pgf.form = "fv_wright"


class ReducedGravity:
    """``ReducedGravity(g_free_surface=9.81, g_internal=0.0098)``
    -> `form="gprime"`, `gprime_gfs`, `gprime_gint`. Two-layer only."""

    def __init__(self, g_free_surface=9.81, g_internal=0.0098):
        self.g_free_surface, self.g_internal = g_free_surface, g_internal

    def apply(self, config):
        config.ocean_pgf.form = "gprime"
        config.ocean_pgf.gprime_gfs = float(self.g_free_surface)
        config.ocean_pgf.gprime_gint = float(self.g_internal)


class FiniteVolumeMOM6:
    """``FiniteVolumeMOM6(gfs_scale=1.0, mass_weight=False,
    reconstruct_for_pressure=False, recon="plm")``

    -> `form="fv_mom6"`, `gfs_scale`, `mass_weight`,
    `reconstruct_for_pressure`, `recon_scheme` (1=PLM, 2=PPM). `recon` is
    spelled `"plm"`/`"ppm"` and mapped to the integer -- a bare 1/2 in a
    Python script is unreadable.
    """

    def __init__(self, gfs_scale=1.0, mass_weight=False,
                 reconstruct_for_pressure=False, recon="plm"):
        if recon not in _RECON:
            raise ValueError(f"recon must be 'plm' or 'ppm', got {recon!r}")
        self.gfs_scale = gfs_scale
        self.mass_weight = mass_weight
        self.reconstruct_for_pressure = reconstruct_for_pressure
        self.recon = recon

    def apply(self, config):
        config.ocean_pgf.form = "fv_mom6"
        config.ocean_pgf.gfs_scale = float(self.gfs_scale)
        config.ocean_pgf.mass_weight = bool(self.mass_weight)
        config.ocean_pgf.reconstruct_for_pressure = bool(
            self.reconstruct_for_pressure)
        config.ocean_pgf.recon_scheme = _RECON[self.recon]


class VelocityTruncation:
    """``VelocityTruncation(maxvel=0.0, cfl_trunc=0.0)`` -> same group.

    A numerical safety net, not a PGF variant -- kept as a separate
    object because CLAUDE.md's `if/else` NaN-laundering rule makes it
    load-bearing: an unguarded clamp turns corruption into plausible
    extreme values under `-fast`.
    """

    def __init__(self, maxvel=0.0, cfl_trunc=0.0):
        self.maxvel, self.cfl_trunc = maxvel, cfl_trunc

    def apply(self, config):
        config.ocean_pgf.maxvel = float(self.maxvel)
        config.ocean_pgf.cfl_trunc = float(self.cfl_trunc)
