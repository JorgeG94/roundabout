"""P5 D2.10 -- vertical coordinate. `&vcoord_nml`. NOT spelled `z=` -- see
D6.4: Roundabout has no explicit vertical-grid input in any form.
"""

from __future__ import annotations


class Sigma:
    """``Sigma(remap="ppm")`` -> `vcoord_type="sigma"`. No remap in
    practice (terrain-following, no vertical remap step)."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "sigma"
        config.vcoord.remap_method = self.remap


class ZSigma:
    """``ZSigma(remap="ppm")`` -> `"zsigma"`. Smoothstep sigma->z
    blend.

    REFUSED by `validate_config` today. The deep branch reads
    `z_ref_global` as absolute depths in metres, but the only writer of
    that table is the dimensionless `k/nz` init, so every z-level
    interval is `1/nz` metres and the whole column collapses into the
    bed layer (with `sum(target_h) == H + eta` still exact, which is why
    it looked healthy). Kept on the schema because it returns once the
    table is filled in metres; use `ZStarSigma` for a sigma/z* blend
    meanwhile."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "zsigma"
        config.vcoord.remap_method = self.remap


class ZStar:
    """``ZStar(remap="ppm")`` -> `"zstar"`. z*-lite: one global
    `z_ref`."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "zstar"
        config.vcoord.remap_method = self.remap


class ZStarSigma:
    """-> `"zstar_sigma"`. Sigma in shallow water, z*-lite in deep."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "zstar_sigma"
        config.vcoord.remap_method = self.remap


class ZStarFull:
    """``ZStarFull(h_surf_target=0.0, h_min=1e-4, remap="ppm")``

    -> `"zstar_full"` + `zstar_h_surf_target`, `zstar_h_min`. Per-column
    `z_ref` from local bathymetry; surface layer anchored at
    `h_surf_target` regardless of `H`.

    NOT offered: `stretching` / `n_surf`. `&vcoord_nml zstar_stretching`
    (default `"log"`) and `zstar_n_surf` are registered on the schema but
    never copied onto `ocean_state%vcoord` -- the ocean path always runs
    `STRETCH_UNIFORM` / `n_surf=0` regardless (D5.6). Offering them here
    would expose two dead knobs as live; use
    `model.config.vcoord.zstar_stretching` directly if you want the
    `RdbDeadKnobWarning` documenting exactly this.

    CAVEAT: intertidal domains leak 1-2% salt/cycle under this
    coordinate; prefer `Sigma`/`ZStar` there.

    CAVEAT: `h_min` must stay at or below `H_VANISHED = 1.5e-4` m. On
    this family it is the anti-zero thickness of filler layers that are
    MEANT to read as vanished downstream, so a larger value promotes them
    to dynamically live (EOS / PGF / remap-drain / vdiff) while the
    coordinate still treats them as throwaway. `validate_config` refuses
    it. For a genuinely live minimum layer thickness use
    `&ocean_isopycnal_nml angstrom_h`; the `Isopycnal` / `Hycom`
    coordinates carry the opposite (keep-alive) contract on the same knob.
    """

    def __init__(self, h_surf_target=0.0, h_min=1e-4, remap="ppm"):
        self.h_surf_target, self.h_min, self.remap = (
            h_surf_target, h_min, remap)

    def apply(self, config):
        config.vcoord.vcoord_type = "zstar_full"
        config.vcoord.remap_method = self.remap
        config.vcoord.zstar_h_surf_target = float(self.h_surf_target)
        config.vcoord.zstar_h_min = float(self.h_min)


class Isopycnal:
    """``Isopycnal(rho_light=1020.0, rho_dense=1030.0, ref_pressure=2e7,
    remap="ppm")`` -> `"rho"` + `rho_target_*` + `rho_ref_pressure`."""

    def __init__(self, rho_light=1020.0, rho_dense=1030.0,
                 ref_pressure=2e7, remap="ppm"):
        self.rho_light, self.rho_dense = rho_light, rho_dense
        self.ref_pressure, self.remap = ref_pressure, remap

    def apply(self, config):
        config.vcoord.vcoord_type = "rho"
        config.vcoord.remap_method = self.remap
        config.vcoord.rho_target_light = float(self.rho_light)
        config.vcoord.rho_target_dense = float(self.rho_dense)
        config.vcoord.rho_ref_pressure = float(self.ref_pressure)


class Hycom:
    """-> `"hycom"`. Hybrid."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "hycom"
        config.vcoord.remap_method = self.remap


class Lagrangian:
    """-> `"lagrangian"`. Pairs with `&ocean_isopycnal_nml` grounding
    knobs, which are generated-layer only."""

    def __init__(self, remap="ppm"):
        self.remap = remap

    def apply(self, config):
        config.vcoord.vcoord_type = "lagrangian"
        config.vcoord.remap_method = self.remap


class ALE:
    """``ALE(remap="ppm", regrid_time_scale=0.0, conserve_ke=False,
    thickness_config="sigma")``

    -> `remap_method`, `regrid_time_scale`, `remap_vel_conserve_ke`,
    `thickness_config`. A MODIFIER on any of the vcoord objects above --
    apply it AFTER the vcoord object so its `remap_method` wins if both
    set one.
    """

    def __init__(self, remap="ppm", regrid_time_scale=0.0,
                 conserve_ke=False, thickness_config="sigma"):
        self.remap = remap
        self.regrid_time_scale = regrid_time_scale
        self.conserve_ke = conserve_ke
        self.thickness_config = thickness_config

    def apply(self, config):
        config.vcoord.remap_method = self.remap
        config.vcoord.regrid_time_scale = float(self.regrid_time_scale)
        config.vcoord.remap_vel_conserve_ke = bool(self.conserve_ke)
        config.vcoord.thickness_config = self.thickness_config
