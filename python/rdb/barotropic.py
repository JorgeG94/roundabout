"""P5 D2.9 -- barotropic solver, continuity, porous barriers, wet/dry.

The `&ocean_bt_nml` MOM6-parity chain (`use_cont_type`,
`cont_corr_bounds`, `upstream_h_face`, `correction_h_weighted`,
`correction_visc_rem`, `renorm_visc_rem`, `forcing_visc_rem`,
`correction_bc_pgf`, `substep_zeta_ke`) is DELIBERATELY generated-layer
only: it has documented prerequisite CHAINS that belong in the Fortran's
own `cross_check` callback, not re-encoded in Python where they would
drift (D2.9). Use `model.config.ocean_bt.<knob>` for those.
"""

from __future__ import annotations


class WaveDrag:
    """``WaveDrag(form="uniform", r=0.0, scale=1.0)`` -> `&ocean_bt_nml
    wave_drag=.true.` + `wave_drag_form`/`r_uniform`/`scale`.
    `form="file"` raises (reserved for PR-14, unused today)."""

    def __init__(self, form="uniform", r=0.0, scale=1.0):
        if form == "file":
            raise ValueError(
                "WaveDrag(form='file') is reserved for PR-14 and unused "
                "today")
        if form not in ("uniform", "roughness_proxy"):
            raise ValueError(f"form must be 'uniform' or "
                              f"'roughness_proxy', got {form!r}")
        self.form, self.r, self.scale = form, r, scale

    def _apply(self, config):
        config.ocean_bt.wave_drag = True
        config.ocean_bt.wave_drag_form = self.form
        config.ocean_bt.wave_drag_r_uniform = float(self.r)
        config.ocean_bt.wave_drag_scale = float(self.scale)


class BarotropicSolver:
    """``BarotropicSolver(n_inner=0, auto=False, cfl_safety=0.65,
    scheme="ssp_rk2", pc_be=0.6, bebt=0.1, halo=None, wave_drag=None,
    substep_drag=False)``

    -> `n_inner`, `auto_n_inner`, `cfl_bt_safety`, `split_scheme`,
    `pc_be`, `bebt`, `bt_halo`, `substep_drag`.

    `auto=False` is the DEFAULT, matching `&ocean_bt_nml auto_n_inner =
    .false.` (rdb_config.F90:499). An earlier version defaulted to
    `auto=True`, so `BarotropicSolver(n_inner=60)` silently ignored the
    `n_inner` it was handed and derived one from the gravity-wave CFL
    instead.
    """

    def __init__(self, n_inner=0, auto=False, cfl_safety=0.65,
                 scheme="ssp_rk2", pc_be=0.6, bebt=0.1, halo=None,
                 wave_drag=None, substep_drag=False):
        self.n_inner, self.auto, self.cfl_safety = n_inner, auto, cfl_safety
        self.scheme, self.pc_be, self.bebt = scheme, pc_be, bebt
        self.halo, self.wave_drag = halo, wave_drag
        self.substep_drag = substep_drag

    def apply(self, config):
        config.ocean_bt.n_inner = int(self.n_inner)
        config.ocean_bt.auto_n_inner = bool(self.auto)
        config.ocean_bt.cfl_bt_safety = float(self.cfl_safety)
        config.ocean_bt.split_scheme = self.scheme
        config.ocean_bt.pc_be = float(self.pc_be)
        config.ocean_bt.bebt = float(self.bebt)
        if self.halo is not None:
            config.ocean_bt.bt_halo = int(self.halo)
        config.ocean_bt.substep_drag = bool(self.substep_drag)
        if self.wave_drag is not None:
            self.wave_drag._apply(config)


class Continuity:
    """``Continuity(h_min=1e-6, positive_definite=False, vol_cfl=False,
    ppm_limit_pos=False)`` -> `&ocean_continuity_nml`."""

    def __init__(self, h_min=1e-6, positive_definite=False, vol_cfl=False,
                 ppm_limit_pos=False):
        self.h_min, self.positive_definite = h_min, positive_definite
        self.vol_cfl, self.ppm_limit_pos = vol_cfl, ppm_limit_pos

    def apply(self, config):
        config.ocean_continuity.h_min = float(self.h_min)
        config.ocean_continuity.positive_definite = bool(
            self.positive_definite)
        config.ocean_continuity.vol_cfl = bool(self.vol_cfl)
        config.ocean_continuity.ppm_limit_pos = bool(self.ppm_limit_pos)


class PorousBarriers:
    """``PorousBarriers(source="resolved", eta_interp="max",
    masking_depth=0.0)`` -> `&ocean_porous_nml enable=.true.`

    `source="file"` raises (deferred in the Fortran). Incompatible with
    `BarotropicSolver(halo>0)` and with `WetDry` -- checked at compose
    time, where both objects are visible.
    """

    def __init__(self, source="resolved", eta_interp="max",
                 masking_depth=0.0):
        if source == "file":
            raise ValueError(
                "PorousBarriers(source='file') is deferred in the "
                "Fortran and not offered")
        self.source, self.eta_interp = source, eta_interp
        self.masking_depth = masking_depth

    def apply(self, config):
        config.ocean_porous.enable = True
        config.ocean_porous.source = self.source
        config.ocean_porous.eta_interp = self.eta_interp
        config.ocean_porous.masking_depth = float(self.masking_depth)


class WetDry:
    """``WetDry(dry_depth=0.05, rewet_depth=0.1, land_margin=5.0)`` ->
    `&ocean_wetdry_nml enable=.true.` Single-rank."""

    def __init__(self, dry_depth=0.05, rewet_depth=0.1, land_margin=5.0):
        self.dry_depth, self.rewet_depth = dry_depth, rewet_depth
        self.land_margin = land_margin

    def apply(self, config):
        config.ocean_wetdry.enable = True
        config.ocean_wetdry.dry_depth = float(self.dry_depth)
        config.ocean_wetdry.rewet_depth = float(self.rewet_depth)
        config.ocean_wetdry.land_margin = float(self.land_margin)
