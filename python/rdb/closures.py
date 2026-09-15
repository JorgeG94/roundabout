"""P5 D2.7-D2.9 -- vertical closures, lateral closures, and bottom drag: one
composed ``closures=[...]`` list on `Model(...)`.

Every object here writes only its own knobs (D2.0 rule 1); the CROSS-object
compose-time checks (KPP xor EPBL, Bryan-Lewis xor Henyey, Henyey needs a
non-cartesian grid, MEKE backscatter needs a non-zero biharmonic backstop,
`ImplicitVerticalFriction(bbl_glue=True)`'s three-prerequisite chain) live in
`rdb._compose._compose_closures`, because they need two (or three)
objects -- and, for Henyey, the grid -- in view at once, which no single
object's `apply()` can see. See 06_python_surface_design.md D2.7/D2.8.
"""

from __future__ import annotations

# ======================================================================
# Vertical
# ======================================================================


class PacanowskiPhilander:
    """PP81 interior closure. -> `&ocean_vmix_nml use_closure=.true.`

    ``PacanowskiPhilander(nu0=1e-2, nu_bg=1e-3, kappa_bg=1e-4, alpha=5.0,
    shear2_floor=1e-9, enabled=True)``

    `use_closure` is ON by default in the namelist, so OMITTING this
    object does not turn PP81 off -- pass `PacanowskiPhilander(enabled=
    False)` for that (mirrors `Diagnostics(enabled=False)` /
    `KPP(enabled=False)`): sets `use_closure = .false.` and every other
    keyword is ignored (nothing to configure).
    """

    def __init__(self, nu0=1e-2, nu_bg=1e-3, kappa_bg=1e-4, alpha=5.0,
                 shear2_floor=1e-9, enabled=True):
        self.nu0, self.nu_bg, self.kappa_bg = nu0, nu_bg, kappa_bg
        self.alpha, self.shear2_floor = alpha, shear2_floor
        self.enabled = bool(enabled)

    def apply(self, config):
        config.ocean_vmix.use_closure = self.enabled
        if not self.enabled:
            return
        config.ocean_vmix.pp81_nu0 = float(self.nu0)
        config.ocean_vmix.pp81_nu_bg = float(self.nu_bg)
        config.ocean_vmix.pp81_kappa_bg = float(self.kappa_bg)
        config.ocean_vmix.pp81_alpha = float(self.alpha)
        config.ocean_vmix.shear2_floor = float(self.shear2_floor)


class KPP:
    """KPP boundary-layer overlay. -> `&ocean_vmix_nml use_kpp=.true.`

    ``KPP(ri_crit=0.3, cs_nonlocal=6.3, c_vt2=1.8,
    shortwave_method="mxl_sw", enabled=True)``

    Mutually exclusive with `EPBL` (checked at compose time -- and
    `EPBL` already turns `use_kpp` OFF itself, since KPP defaults ON in
    the schema; see `EPBL`'s own docstring). `KPP(enabled=False)` is an
    EXPLICIT off-switch for the boundary-layer overlay entirely (mirrors
    `Diagnostics(enabled=False)`): sets `use_kpp = .false.` and every
    other keyword is ignored. Used WITHOUT an `EPBL(...)` in the same
    `closures=[...]` list, this leaves NO surface boundary-layer scheme
    at all (PP81 interior + background floors + convective adjustment
    still run) -- the composer warns (`NoBoundaryLayerWarning`) unless
    that is clearly deliberate (e.g. a closed-basin resting-state
    dyn-core test with no wind).

    Requires an interior closure (`use_closure`, on by default).
    """

    def __init__(self, ri_crit=0.3, cs_nonlocal=6.3, c_vt2=1.8,
                 shortwave_method="mxl_sw", enabled=True):
        self.ri_crit, self.cs_nonlocal, self.c_vt2 = (
            ri_crit, cs_nonlocal, c_vt2)
        self.shortwave_method = shortwave_method
        self.enabled = bool(enabled)

    def apply(self, config):
        config.ocean_vmix.use_kpp = self.enabled
        if not self.enabled:
            return
        config.ocean_vmix.kpp_ri_crit = float(self.ri_crit)
        config.ocean_vmix.kpp_cs_nonlocal = float(self.cs_nonlocal)
        config.ocean_vmix.kpp_c_vt2 = float(self.c_vt2)
        config.ocean_thermo.kpp_sw_method = self.shortwave_method


class Langmuir:
    """``Langmuir(scheme="rescale", coef=0.447, exponent=-1.33,
    max_enhance=5.0, sl_fraction=0.04)`` -> `&ocean_epbl_nml lt_*`.
    Passed as `EPBL(langmuir=...)`."""

    def __init__(self, scheme="rescale", coef=0.447, exponent=-1.33,
                 max_enhance=5.0, sl_fraction=0.04):
        self.scheme, self.coef, self.exponent = scheme, coef, exponent
        self.max_enhance, self.sl_fraction = max_enhance, sl_fraction

    def _apply(self, config):
        config.ocean_epbl.use_lt = True
        config.ocean_epbl.lt_scheme = self.scheme
        config.ocean_epbl.lt_enhance_coef = float(self.coef)
        config.ocean_epbl.lt_enhance_exp = float(self.exponent)
        config.ocean_epbl.lt_max_enhance = float(self.max_enhance)
        config.ocean_epbl.la_frac_hbl = float(self.sl_fraction)


class EPBL:
    """Reichl-Hallberg energetics PBL. -> `&ocean_epbl_nml enable=.true.`

    ``EPBL(mstar_scheme="om4", mstar=1.2, nstar=0.2, prandtl=1.0,
    combine="add", mld_iteration=True, mld_use_prev_guess=False,
    langmuir=None)``

    Mutually exclusive with `KPP` (checked at compose time) -- and this
    object ITSELF sets `&ocean_vmix_nml use_kpp = .false.` (KPP defaults
    ON in the schema, so selecting EPBL structurally implies deselecting
    KPP; this is the same "stated requirement" pattern `GentMcWilliams`
    uses for `&ocean_slopes_nml`, not a new rule). Nothing further is
    needed to turn KPP off when composing `EPBL(...)`.

    46 knobs live in this group; this constructor curates the ~8 a user
    tunes -- the RH18 fit coefficients, the mixing-length shape, the
    MLD root-find controls beyond `mld_iteration`/`mld_use_prev_guess`,
    and the five Langmuir `lt_lac*` coefficients are
    `model.config.ocean_epbl.*` only.
    """

    def __init__(self, mstar_scheme="om4", mstar=1.2, nstar=0.2,
                 prandtl=1.0, combine="add", mld_iteration=True,
                 mld_use_prev_guess=False, langmuir=None):
        self.mstar_scheme, self.mstar, self.nstar = (
            mstar_scheme, mstar, nstar)
        self.prandtl, self.combine = prandtl, combine
        self.mld_iteration = mld_iteration
        self.mld_use_prev_guess = mld_use_prev_guess
        self.langmuir = langmuir

    def apply(self, config):
        config.ocean_epbl.enable = True
        config.ocean_vmix.use_kpp = False
        config.ocean_epbl.mstar_scheme = self.mstar_scheme
        config.ocean_epbl.mstar = float(self.mstar)
        config.ocean_epbl.nstar = float(self.nstar)
        config.ocean_epbl.prandtl = float(self.prandtl)
        config.ocean_epbl.combine = self.combine
        config.ocean_epbl.mld_iteration = bool(self.mld_iteration)
        config.ocean_epbl.mld_use_prev_guess = bool(self.mld_use_prev_guess)
        if self.langmuir is not None:
            self.langmuir._apply(config)


class KappaShear:
    """JHL08 prognostic interior shear turbulence. ->
    `&ocean_kappa_shear_nml`

    ``KappaShear(ri_crit=0.25, prandtl=1.0, at_vertex=False,
    vertex_geometric_mean=False, vertex_kdmin=0.0)``

    `at_vertex=True` is the MOM6 VERTEX_SHEAR/OM5 corner solve; it raises
    `required_halo()` to 2 and suppresses the cell-centred `kv` merge.
    The 12 iteration-control knobs go to `model.config.ocean_kappa_shear`.
    """

    def __init__(self, ri_crit=0.25, prandtl=1.0, at_vertex=False,
                 vertex_geometric_mean=False, vertex_kdmin=0.0):
        self.ri_crit, self.prandtl = ri_crit, prandtl
        self.at_vertex = at_vertex
        self.vertex_geometric_mean = vertex_geometric_mean
        self.vertex_kdmin = vertex_kdmin

    def apply(self, config):
        config.ocean_kappa_shear.enable = True
        config.ocean_kappa_shear.ri_crit = float(self.ri_crit)
        config.ocean_kappa_shear.prandtl_turb = float(self.prandtl)
        config.ocean_kappa_shear.at_vertex = bool(self.at_vertex)
        config.ocean_kappa_shear.vertex_geometric_mean = bool(
            self.vertex_geometric_mean)
        config.ocean_kappa_shear.vertex_geomean_kdmin = float(
            self.vertex_kdmin)


class TidalMixing:
    """St-Laurent/Simmons. -> `&ocean_tidal_mixing_nml enable=.true.`

    ``TidalMixing(gamma=0.3333, mu=0.2, decay_scale=500.0, energy=0.0,
    kd_max=1e-2, prandtl=1.0)``

    `energy` -> `e_uniform`. The v1.1 state-dependent `E` recompute is
    generated-layer only.
    """

    def __init__(self, gamma=0.3333, mu=0.2, decay_scale=500.0, energy=0.0,
                 kd_max=1e-2, prandtl=1.0):
        self.gamma, self.mu, self.decay_scale = gamma, mu, decay_scale
        self.energy, self.kd_max, self.prandtl = energy, kd_max, prandtl

    def apply(self, config):
        config.ocean_tidal_mixing.enable = True
        config.ocean_tidal_mixing.gamma = float(self.gamma)
        config.ocean_tidal_mixing.mu = float(self.mu)
        config.ocean_tidal_mixing.zeta = float(self.decay_scale)
        config.ocean_tidal_mixing.e_uniform = float(self.energy)
        config.ocean_tidal_mixing.kd_max = float(self.kd_max)
        config.ocean_tidal_mixing.prandtl_tidal = float(self.prandtl)


class ConvectiveAdjustment:
    """-> `&ocean_conv_nml enable=.true.`

    ``ConvectiveAdjustment(kd=1.0, prandtl=1.0, n2_threshold=0.0)``
    Writes `kv`/`kt` only, never `ks` (a `max()` contributor)."""

    def __init__(self, kd=1.0, prandtl=1.0, n2_threshold=0.0):
        self.kd, self.prandtl, self.n2_threshold = kd, prandtl, n2_threshold

    def apply(self, config):
        config.ocean_conv.enable = True
        config.ocean_conv.kd_conv = float(self.kd)
        config.ocean_conv.prandtl_conv = float(self.prandtl)
        config.ocean_conv.n2_thresh = float(self.n2_threshold)


class DoubleDiffusion:
    """Salt fingering + diffusive convection. -> `&ocean_ddiff_nml`

    ``DoubleDiffusion(kappa_s=1e-4, strat_param_max=2.55, form="mc76")``
    `form="k90"` -> `use_k90=.true.`. Folds the asymmetric `Kd_extra`
    into the `kt`/`ks` split -- the ONLY way `ks` stops being identical
    to `kt`.
    """

    def __init__(self, kappa_s=1e-4, strat_param_max=2.55, form="mc76"):
        if form not in ("mc76", "k90"):
            raise ValueError(f"form must be 'mc76' or 'k90', got {form!r}")
        self.kappa_s, self.strat_param_max, self.form = (
            kappa_s, strat_param_max, form)

    def apply(self, config):
        config.ocean_ddiff.enable = True
        config.ocean_ddiff.kappa_ddiff_s = float(self.kappa_s)
        config.ocean_ddiff.strat_param_max = float(self.strat_param_max)
        config.ocean_ddiff.use_k90 = (self.form == "k90")


class BryanLewisBackground:
    """Depth-varying background Kd. -> `&ocean_vmix_nml
    bkgnd_profile=.true.`

    ``BryanLewisBackground(kd_surface=1e-5, kd_deep=1.3e-4, z0=2500.0,
    delta=222.0, prandtl=1.0)``

    Mutually exclusive with `HenyeyBackground` (MOM6's one-background
    rule, checked at compose time).
    """

    def __init__(self, kd_surface=1e-5, kd_deep=1.3e-4, z0=2500.0,
                 delta=222.0, prandtl=1.0):
        self.kd_surface, self.kd_deep = kd_surface, kd_deep
        self.z0, self.delta, self.prandtl = z0, delta, prandtl

    def apply(self, config):
        config.ocean_vmix.bkgnd_profile = True
        config.ocean_vmix.bkgnd_kd_sfc = float(self.kd_surface)
        config.ocean_vmix.bkgnd_kd_deep = float(self.kd_deep)
        config.ocean_vmix.bkgnd_z0 = float(self.z0)
        config.ocean_vmix.bkgnd_delta = float(self.delta)
        config.ocean_vmix.bkgnd_prandtl = float(self.prandtl)


class HenyeyBackground:
    """Latitude-scaled background Kd. -> `bkgnd_henyey=.true.`

    ``HenyeyBackground(kd_min=-1.0, n0_2omega=20.0, max_lat=95.0)``

    Requires a non-cartesian grid: on a cartesian grid `geolatT == 0`
    makes every column equatorial, and the Fortran fails loud. The
    composer checks the grid type here and says so first.
    """

    def __init__(self, kd_min=-1.0, n0_2omega=20.0, max_lat=95.0):
        self.kd_min, self.n0_2omega, self.max_lat = (
            kd_min, n0_2omega, max_lat)

    def apply(self, config):
        config.ocean_vmix.bkgnd_henyey = True
        config.ocean_vmix.bkgnd_kd_min = float(self.kd_min)
        config.ocean_vmix.bkgnd_henyey_n0_2omega = float(self.n0_2omega)
        config.ocean_vmix.bkgnd_henyey_max_lat = float(self.max_lat)


class ConstantBackground:
    """``ConstantBackground(kv=1e-3, kt=1e-4)`` -> `pp81_nu_bg` /
    `pp81_kappa_bg`, the scalar floor `BryanLewis`/`Henyey` replace."""

    def __init__(self, kv=1e-3, kt=1e-4):
        self.kv, self.kt = kv, kt

    def apply(self, config):
        config.ocean_vmix.pp81_nu_bg = float(self.kv)
        config.ocean_vmix.pp81_kappa_bg = float(self.kt)


class NearSurfaceViscosity:
    """MOM6 `KV_ML_INVZ2` near-surface viscosity band. ->
    `&ocean_vmix_nml kv_ml_invz2`/`hmix_fixed`.

    ``NearSurfaceViscosity(kv=0.0, hmix=20.0)``

    Adds `kv_extra(z) = kv * (hmix/z)^2` over the top `hmix` metres to
    `kv` ONLY -- never `kt`/`ks` (`vmix_add_kv_ml_invz2`,
    `rdb_ocean_vmix.F90:1121`). This is a MOMENTUM-side surface
    viscosity floor, not a boundary-layer scheme: it COMPOSES with
    `KPP(...)`/`EPBL(...)`/neither, rather than selecting between them
    -- do not confuse it with the `KPP`/`EPBL` toggle.
    """

    def __init__(self, kv=0.0, hmix=20.0):
        self.kv, self.hmix = kv, hmix

    def apply(self, config):
        config.ocean_vmix.kv_ml_invz2 = float(self.kv)
        config.ocean_vmix.hmix_fixed = float(self.hmix)


class HarmonicFaceThickness:
    """MOM6 `HARMONIC_VISC` vdiff-assembly parity. -> `&ocean_vmix_nml
    harmonic_visc = .true.`

    ``HarmonicFaceThickness()`` -- harmonic-mean (rather than
    arithmetic-mean) face thickness when building the vertical-friction
    matrix. A pure discretisation choice, independent of which (if any)
    boundary-layer scheme is active; presence in `closures=[...]` is
    the enable switch, matching `Island()`/`DoubleDrake()`'s pattern.
    """

    def apply(self, config):
        config.ocean_vmix.harmonic_visc = True


class ThermoSubcycling:
    """MOM6 `DT_THERM` sub-cycling ratio. -> `&ocean_vmix_nml
    dt_therm_ratio`.

    ``ThermoSubcycling(ratio=1)`` -- the thermo step (EOS / vertical
    mixing / ALE remap) runs at `ratio * dt_dyn` instead of every
    dynamics step. `ratio=1` (the default) is a no-op; this object
    exists for `ratio > 1`.
    """

    def __init__(self, ratio=1):
        self.ratio = ratio

    def apply(self, config):
        config.ocean_vmix.dt_therm_ratio = int(self.ratio)


class ImplicitVerticalFriction:
    """-> `&ocean_vdiff_nml`

    ``ImplicitVerticalFriction(stress=False, drag=False,
    harmonic_thickness=False, bbl_glue=False, hbbl=10.0,
    bbl_piston=3e-4)``

    `bbl_glue` has THREE prerequisites: `harmonic_thickness` (->
    `hvel_mom6`) + `drag` (-> `implicit_drag`) + a LINEAR bottom drag
    (lives on a different object -- `LinearDrag` in the same
    `closures=` list). The composer checks all three and names whichever
    is missing.
    """

    def __init__(self, stress=False, drag=False, harmonic_thickness=False,
                 bbl_glue=False, hbbl=10.0, bbl_piston=3e-4):
        self.stress, self.drag = stress, drag
        self.harmonic_thickness, self.bbl_glue = (
            harmonic_thickness, bbl_glue)
        self.hbbl, self.bbl_piston = hbbl, bbl_piston

    def apply(self, config):
        config.ocean_vdiff.implicit_stress = bool(self.stress)
        config.ocean_vdiff.implicit_drag = bool(self.drag)
        config.ocean_vdiff.hvel_mom6 = bool(self.harmonic_thickness)
        config.ocean_vdiff.bbl_glue = bool(self.bbl_glue)
        config.ocean_vdiff.hbbl_visc = float(self.hbbl)
        config.ocean_vdiff.bbl_piston = float(self.bbl_piston)


class DiffusivityLimits:
    """``DiffusivityLimits(kv_max=1e6, kd_max=1e6, smoothing=0)`` ->
    `&ocean_vmix_nml kv_max`, `kd_max`, `kd_smooth_iterations`. The
    single floors/ceilings/smoothing gate over `kv`/`kt`/`ks`."""

    def __init__(self, kv_max=1e6, kd_max=1e6, smoothing=0):
        self.kv_max, self.kd_max, self.smoothing = kv_max, kd_max, smoothing

    def apply(self, config):
        config.ocean_vmix.kv_max = float(self.kv_max)
        config.ocean_vmix.kd_max = float(self.kd_max)
        config.ocean_vmix.kd_smooth_iterations = int(self.smoothing)


# ======================================================================
# Lateral
# ======================================================================


class Anisotropy:
    """``Anisotropy(kh, direction=(1.0, 0.0))`` MODIFIER, passed as
    `Smagorinsky(anisotropy=...)`. Requires `form="stress_tensor"`."""

    def __init__(self, kh, direction=(1.0, 0.0)):
        self.kh, self.direction = kh, direction


class Smagorinsky:
    """-> `&ocean_hvisc_nml lateral_closure="smagorinsky"`

    ``Smagorinsky(C=0.15, biharmonic=False, C_biharmonic=0.06,
    form="laplacian", free_slip=True, bound_kh=False, bound_coef=0.8,
    resolution_scaled=False, anisotropy=None, ah_max=1e4, nu_4_max=1e12,
    nu_h=0.0, kh_vel_scale=0.0, nu_4=0.0)``

    `nu_h`/`kh_vel_scale`/`nu_4` are the CONSTANT-floor knobs that
    compose ALONGSIDE the flow-aware Smagorinsky term -- real configs
    pair them (`nu_h=10000` + Smagorinsky is the production envelope,
    CLAUDE.md's double-gyre reference). They are not exclusive with
    `lateral_closure="smagorinsky"`; use `ConstantViscosity` instead if
    you want ONLY a constant floor with no flow-aware term at all.

    `form="stress_tensor"` -> `stress_tensor=.true.` (MOM6
    thickness-weighted stress-div operator); it COMPOSES with the
    biharmonic add-on, matching MOM6.

    `free_slip=True` is the DEFAULT, matching the Fortran schema's own
    `&ocean_hvisc_nml no_slip = .false.` (rdb_config.F90:864). An earlier
    version defaulted to `free_slip=False` and claimed that matched the
    Fortran -- it did not; it wrote `no_slip = .true.`, silently inverting
    the wall boundary condition relative to a namelist that omits the knob.
    Pass `free_slip=False` explicitly for no-slip walls.

    `anisotropy=Anisotropy(...)` requires `form="stress_tensor"`; passing
    it otherwise raises.
    """

    def __init__(self, C=0.15, biharmonic=False, C_biharmonic=0.06,
                 form="laplacian", free_slip=True, bound_kh=False,
                 bound_coef=0.8, resolution_scaled=False, anisotropy=None,
                 ah_max=1e4, nu_4_max=1e12, nu_h=0.0, kh_vel_scale=0.0,
                 nu_4=0.0):
        if form not in ("laplacian", "stress_tensor"):
            raise ValueError(f"form must be 'laplacian' or "
                              f"'stress_tensor', got {form!r}")
        if anisotropy is not None and form != "stress_tensor":
            raise ValueError(
                "anisotropy= requires form='stress_tensor'")
        self.C, self.biharmonic, self.C_biharmonic = C, biharmonic, \
            C_biharmonic
        self.form, self.free_slip = form, free_slip
        self.bound_kh, self.bound_coef = bound_kh, bound_coef
        self.resolution_scaled, self.anisotropy = (
            resolution_scaled, anisotropy)
        self.ah_max, self.nu_4_max = ah_max, nu_4_max
        self.nu_h, self.kh_vel_scale, self.nu_4 = nu_h, kh_vel_scale, nu_4

    def biharmonic_coeff(self):
        """The dissipation coefficient the MEKE-backscatter compose-time
        check inspects -- the larger of the flow-aware biharmonic term
        (0.0 unless `biharmonic=True`) and the constant `nu_4` floor,
        since either alone is a real, non-zero dissipation backstop."""
        flow_aware = float(self.C_biharmonic) if self.biharmonic else 0.0
        return max(flow_aware, float(self.nu_4))

    def apply(self, config):
        config.ocean_hvisc.lateral_closure = "smagorinsky"
        config.ocean_hvisc.c_smag = float(self.C)
        config.ocean_hvisc.smag_ah = bool(self.biharmonic)
        config.ocean_hvisc.smag_bi_const = float(self.C_biharmonic)
        config.ocean_hvisc.stress_tensor = (self.form == "stress_tensor")
        config.ocean_hvisc.no_slip = not self.free_slip
        config.ocean_hvisc.bound_kh = bool(self.bound_kh)
        config.ocean_hvisc.bound_coef = float(self.bound_coef)
        config.ocean_hvisc.resoln_scaled_visc = bool(self.resolution_scaled)
        config.ocean_hvisc.ah_max = float(self.ah_max)
        config.ocean_hvisc.nu_4_max = float(self.nu_4_max)
        config.ocean_hvisc.nu_h = float(self.nu_h)
        config.ocean_hvisc.kh_vel_scale = float(self.kh_vel_scale)
        config.ocean_hvisc.nu_4 = float(self.nu_4)
        if self.anisotropy is not None:
            config.ocean_hvisc.kh_aniso = float(self.anisotropy.kh)
            config.ocean_hvisc.aniso_dir = list(self.anisotropy.direction)


class Leith:
    """-> `lateral_closure="leith"`. Same keyword set as `Smagorinsky`
    (``C`` here scales `c_leith`/`c_leith_bi`), including the
    `nu_h`/`kh_vel_scale`/`nu_4` constant-floor knobs that compose
    alongside the flow-aware term."""

    def __init__(self, C=1.0, biharmonic=False, C_biharmonic=0.0,
                 free_slip=True, bound_kh=False, bound_coef=0.8,
                 ah_max=1e4, nu_4_max=1e12, nu_h=0.0, kh_vel_scale=0.0,
                 nu_4=0.0):
        self.C, self.biharmonic, self.C_biharmonic = C, biharmonic, \
            C_biharmonic
        self.free_slip, self.bound_kh, self.bound_coef = (
            free_slip, bound_kh, bound_coef)
        self.ah_max, self.nu_4_max = ah_max, nu_4_max
        self.nu_h, self.kh_vel_scale, self.nu_4 = nu_h, kh_vel_scale, nu_4

    def biharmonic_coeff(self):
        flow_aware = float(self.C_biharmonic) if self.biharmonic else 0.0
        return max(flow_aware, float(self.nu_4))

    def apply(self, config):
        config.ocean_hvisc.lateral_closure = (
            "leith_biharm" if self.biharmonic else "leith")
        config.ocean_hvisc.c_leith = float(self.C)
        config.ocean_hvisc.c_leith_bi = float(self.C_biharmonic)
        config.ocean_hvisc.no_slip = not self.free_slip
        config.ocean_hvisc.bound_kh = bool(self.bound_kh)
        config.ocean_hvisc.bound_coef = float(self.bound_coef)
        config.ocean_hvisc.ah_max = float(self.ah_max)
        config.ocean_hvisc.nu_4_max = float(self.nu_4_max)
        config.ocean_hvisc.nu_h = float(self.nu_h)
        config.ocean_hvisc.kh_vel_scale = float(self.kh_vel_scale)
        config.ocean_hvisc.nu_4 = float(self.nu_4)


class ConstantViscosity:
    """``ConstantViscosity(nu_h=0.0, nu_4=0.0)`` -> `lateral_closure="none"`
    + `nu_h`/`nu_4`.

    Warns above `nu_4 ~ 3e8`: the forward-Euler biharmonic CFL caps it
    there, so MOM6's 1e9-1e10 is not reachable explicitly (memory
    `biharmonic-cfl`) -- better a warning at construction than a blow-up
    on day 2.
    """

    def __init__(self, nu_h=0.0, nu_4=0.0):
        if nu_4 > 3e8:
            import warnings
            warnings.warn(
                f"ConstantViscosity(nu_4={nu_4!r}): the forward-Euler "
                f"biharmonic CFL caps nu_4 around 3e8 -- this run is "
                f"likely to blow up (memory: biharmonic-cfl). MOM6's "
                f"1e9-1e10 range is not reachable with an explicit "
                f"biharmonic term.", stacklevel=2)
        self.nu_h, self.nu_4 = nu_h, nu_4

    def apply(self, config):
        config.ocean_hvisc.lateral_closure = "none"
        config.ocean_hvisc.nu_h = float(self.nu_h)
        config.ocean_hvisc.nu_4 = float(self.nu_4)

    def biharmonic_coeff(self):
        return float(self.nu_4)


class HorizontalDiffusivity:
    """``HorizontalDiffusivity(kappa_h=0.0)`` -> `&ocean_hdiff_nml
    kappa_h`. ALONG-COORDINATE, not neutral -- use `Redi` for neutral
    diffusion."""

    def __init__(self, kappa_h=0.0):
        self.kappa_h = kappa_h

    def apply(self, config):
        config.ocean_hdiff.kappa_h = float(self.kappa_h)


class MEKE:
    """Prognostic mesoscale EKE. -> `&ocean_meke_nml enable=.true.`

    ``MEKE(gm_coeff=0.0, fr_coeff=0.0, damping=0.0, kh_coeff=1.0,
    backscatter=False, backscatter_ku=0.0, cdrag=2.5e-3)``

    Requires `GentMcWilliams` in the same `closures=` list.
    `backscatter=True` additionally requires a biharmonic backstop with a
    NON-ZERO coefficient -- the composer checks the lateral object in the
    list and raises if it carries none; it never invents one (picking a
    dissipation scale for the user is guessing at physics).
    """

    def __init__(self, gm_coeff=0.0, fr_coeff=0.0, damping=0.0,
                 kh_coeff=1.0, backscatter=False, backscatter_ku=0.0,
                 cdrag=2.5e-3):
        self.gm_coeff, self.fr_coeff, self.damping = (
            gm_coeff, fr_coeff, damping)
        self.kh_coeff = kh_coeff
        self.backscatter, self.backscatter_ku = backscatter, backscatter_ku
        self.cdrag = cdrag

    def apply(self, config):
        config.ocean_meke.enable = True
        config.ocean_meke.gmcoeff = float(self.gm_coeff)
        config.ocean_meke.frcoeff = float(self.fr_coeff)
        config.ocean_meke.damping = float(self.damping)
        config.ocean_meke.khcoeff = float(self.kh_coeff)
        config.ocean_meke.backscatter = bool(self.backscatter)
        config.ocean_meke.backscatter_visc_coeff_ku = float(
            self.backscatter_ku)
        config.ocean_meke.cdrag = float(self.cdrag)


class GentMcWilliams:
    """-> `&ocean_gm_nml enable=.true.`

    ``GentMcWilliams(kappa=0.0, max_cfl=0.1, slope_max=0.01)``
    Implicitly enables `&ocean_slopes_nml` (its stated requirement); the
    composer adds that automatically.
    """

    def __init__(self, kappa=0.0, max_cfl=0.1, slope_max=0.01):
        self.kappa, self.max_cfl, self.slope_max = kappa, max_cfl, slope_max

    def apply(self, config):
        config.ocean_gm.enable = True
        config.ocean_gm.khth = float(self.kappa)
        config.ocean_gm.khth_max_cfl = float(self.max_cfl)
        config.ocean_gm.khth_slope_max = float(self.slope_max)
        config.ocean_slopes.enable = True


class Redi:
    """``Redi(kappa=0.0)`` -> `&ocean_redi_nml enable=.true.`, `khtr`.
    Also implies slopes. `continuous=False` is rejected by the Fortran,
    so it is not an argument."""

    def __init__(self, kappa=0.0):
        self.kappa = kappa

    def apply(self, config):
        config.ocean_redi.enable = True
        config.ocean_redi.khtr = float(self.kappa)
        config.ocean_redi.continuous = True
        config.ocean_slopes.enable = True


class IsopycnalSlopes:
    """``IsopycnalSlopes(kd_smooth=1e-6, min_dz_for_n2=1.0)`` ->
    `&ocean_slopes_nml`. Rarely constructed directly (`GentMcWilliams`/
    `Redi`/`VariableMixing` imply it); offered so the smoothing knobs
    are reachable curated."""

    def __init__(self, kd_smooth=1e-6, min_dz_for_n2=1.0):
        self.kd_smooth, self.min_dz_for_n2 = kd_smooth, min_dz_for_n2

    def apply(self, config):
        config.ocean_slopes.enable = True
        config.ocean_slopes.kd_smooth = float(self.kd_smooth)
        config.ocean_slopes.min_dz_for_n2 = float(self.min_dz_for_n2)


class VariableMixing:
    """Spatially-varying GM/Redi coefficients. -> `&ocean_varmix_nml`

    ``VariableMixing(visbeck=False, resolution_scaled_khth=False,
    resolution_scaled_khtr=False, khth=0.0, khtr=0.0, res_fn_power=2,
    res_scale_coef=1.0)``

    Requires slopes + `WaveSpeed` (both implied automatically by the
    composer).
    """

    def __init__(self, visbeck=False, resolution_scaled_khth=False,
                 resolution_scaled_khtr=False, khth=0.0, khtr=0.0,
                 res_fn_power=2, res_scale_coef=1.0):
        self.visbeck = visbeck
        self.resolution_scaled_khth = resolution_scaled_khth
        self.resolution_scaled_khtr = resolution_scaled_khtr
        self.khth, self.khtr = khth, khtr
        self.res_fn_power, self.res_scale_coef = (
            res_fn_power, res_scale_coef)

    def apply(self, config):
        config.ocean_varmix.enable = True
        config.ocean_varmix.use_visbeck = bool(self.visbeck)
        config.ocean_varmix.resoln_scaled_khth = bool(
            self.resolution_scaled_khth)
        config.ocean_varmix.resoln_scaled_khtr = bool(
            self.resolution_scaled_khtr)
        config.ocean_varmix.khth = float(self.khth)
        config.ocean_varmix.khtr = float(self.khtr)
        config.ocean_varmix.kh_res_fn_power = int(self.res_fn_power)
        config.ocean_varmix.kh_res_scale_coef = float(self.res_scale_coef)
        config.ocean_slopes.enable = True
        config.ocean_wavespeed.enable = True


class Bodner:
    """``Bodner(cr=..., mstar=..., nstar=...)`` MODIFIER, passed as
    `MixedLayerRestratification(bodner=...)`. Bodner (2023) TTW/
    frontogenesis-arrest refinement -- OVERRIDES `ce`/`mom_mixrate`."""

    def __init__(self, cr=1.0, mstar=0.75, nstar=1.0):
        self.cr, self.mstar, self.nstar = cr, mstar, nstar


class MixedLayerRestratification:
    """Fox-Kemper MLE. -> `&ocean_foxkemper_nml enable=.true.`

    ``MixedLayerRestratification(ce=0.0625, mom_mixrate=False,
    mld_decay_time=0.0, bodner=None)``

    `bodner=Bodner(...)` OVERRIDES `ce`/`mom_mixrate` (the Fortran's own
    rule); passing both `ce` (non-default) and `bodner` raises rather
    than silently picking one.
    """

    def __init__(self, ce=0.0625, mom_mixrate=False, mld_decay_time=0.0,
                 bodner=None):
        if bodner is not None and ce != 0.0625:
            raise ValueError(
                "bodner= overrides ce/mom_mixrate in the Fortran; passing "
                "a non-default ce= together with bodner= would silently "
                "be ignored, so this raises instead")
        self.ce, self.mom_mixrate = ce, mom_mixrate
        self.mld_decay_time, self.bodner = mld_decay_time, bodner

    def apply(self, config):
        config.ocean_foxkemper.enable = True
        config.ocean_foxkemper.mld_decay_time = float(self.mld_decay_time)
        if self.bodner is not None:
            config.ocean_foxkemper.use_bodner = True
            config.ocean_foxkemper.cr = float(self.bodner.cr)
            config.ocean_foxkemper.bodner_mstar = float(self.bodner.mstar)
            config.ocean_foxkemper.bodner_nstar = float(self.bodner.nstar)
        else:
            config.ocean_foxkemper.ce = float(self.ce)
            config.ocean_foxkemper.use_mom_mixrate = bool(self.mom_mixrate)


class WaveSpeed:
    """``WaveSpeed(cadence=1)`` -> `&ocean_wavespeed_nml enable=.true.`,
    `n_wavespeed`. Diagnostic; required by `VariableMixing` (implied
    automatically). `mono_n2`/`use_ebt` are DEFERRED in the Fortran and
    are not offered."""

    def __init__(self, cadence=1):
        self.cadence = cadence

    def apply(self, config):
        config.ocean_wavespeed.enable = True
        config.ocean_wavespeed.n_wavespeed = int(self.cadence)


# ======================================================================
# Bottom drag   ->  &ocean_bdrag_nml
# ======================================================================


class QuadraticDrag:
    """``QuadraticDrag(cd=0.0, hbbl=0.0, bg_vel=0.0, bbl_thick_min=0.0,
    bed_factor=1.0, implicit=False)`` -> `form="quadratic"`."""

    def __init__(self, cd=0.0, hbbl=0.0, bg_vel=0.0, bbl_thick_min=0.0,
                 bed_factor=1.0, implicit=False):
        self.cd, self.hbbl, self.bg_vel = cd, hbbl, bg_vel
        self.bbl_thick_min, self.bed_factor = bbl_thick_min, bed_factor
        self.implicit = implicit

    def apply(self, config):
        config.ocean_bdrag.form = "quadratic"
        config.ocean_bdrag.cd = float(self.cd)
        config.ocean_bdrag.hbbl = float(self.hbbl)
        config.ocean_bdrag.bg_vel = float(self.bg_vel)
        config.ocean_bdrag.bbl_thick_min = float(self.bbl_thick_min)
        config.ocean_bdrag.bed_factor = float(self.bed_factor)
        config.ocean_bdrag.implicit = bool(self.implicit)


class LinearDrag:
    """``LinearDrag(r=0.0, hbbl=0.0, bg_vel=0.0, bbl_thick_min=0.0,
    bed_factor=1.0, implicit=False)`` -> `form="linear"`. Required by
    `ImplicitVerticalFriction(bbl_glue=True)`."""

    def __init__(self, r=0.0, hbbl=0.0, bg_vel=0.0, bbl_thick_min=0.0,
                 bed_factor=1.0, implicit=False):
        self.r, self.hbbl, self.bg_vel = r, hbbl, bg_vel
        self.bbl_thick_min, self.bed_factor = bbl_thick_min, bed_factor
        self.implicit = implicit

    def apply(self, config):
        config.ocean_bdrag.form = "linear"
        config.ocean_bdrag.r = float(self.r)
        config.ocean_bdrag.hbbl = float(self.hbbl)
        config.ocean_bdrag.bg_vel = float(self.bg_vel)
        config.ocean_bdrag.bbl_thick_min = float(self.bbl_thick_min)
        config.ocean_bdrag.bed_factor = float(self.bed_factor)
        config.ocean_bdrag.implicit = bool(self.implicit)


class ChannelDrag:
    """``ChannelDrag(cdrag_side=0.0)`` -> `channel_drag=.true.` +
    `cdrag_side`. A MODIFIER, composed alongside a `QuadraticDrag`/
    `LinearDrag` in the same `closures=` list."""

    def __init__(self, cdrag_side=0.0):
        self.cdrag_side = cdrag_side

    def apply(self, config):
        config.ocean_bdrag.channel_drag = True
        config.ocean_bdrag.cdrag_side = float(self.cdrag_side)
