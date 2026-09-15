"""P5 D2.13 -- surface forcing and initial-condition overlays."""

from __future__ import annotations

from pathlib import Path

from ._errors import RdbIOError

# ======================================================================
# Wind
# ======================================================================


class ConstantWind:
    """``ConstantWind(tau_x=0.0, tau_y=0.0)`` -> `&physics_nml
    wind_stress_x/y` with `&ocean_topo_nml wind_config="constant"`."""

    def __init__(self, tau_x=0.0, tau_y=0.0):
        self.tau_x, self.tau_y = tau_x, tau_y

    def apply(self, config):
        config.ocean_topo.wind_config = "constant"
        config.physics.wind_stress_x = float(self.tau_x)
        config.physics.wind_stress_y = float(self.tau_y)


class DoubleGyreWind:
    """``DoubleGyreWind(magnitude=0.1)`` -> `wind_config="2gyre"`,
    `taux_magnitude`."""

    def __init__(self, magnitude=0.1):
        self.magnitude = magnitude

    def apply(self, config):
        config.ocean_topo.wind_config = "2gyre"
        config.ocean_topo.taux_magnitude = float(self.magnitude)


class Neverworld2Wind:
    """-> `wind_config="neverworld2"`."""

    def apply(self, config):
        config.ocean_topo.wind_config = "neverworld2"


# ======================================================================
# Heat / salt / shortwave
# ======================================================================


class Thermodynamics:
    """``Thermodynamics(enabled=True)`` -> `&ocean_thermo_nml
    enable_thermodynamics`.

    `enabled=False` is the MOM6 `ENABLE_THERMODYNAMICS = False` analogue
    -- Roundabout's multilayer kernel always advects T/S regardless (there
    is no way to skip that), but this skips the EOS call entirely, so
    density stays at `rho_0` and T/S drift cannot reach the dynamics
    (see `double_gyre_mom6.py`'s own note on why it ALSO zeroes the
    coastal-legacy `tracer.alpha_T`/`beta_S` belt-and-braces -- those
    are dead on the ocean path either way, per the P4 dead-knob sweep).
    """

    def __init__(self, enabled=True):
        self.enabled = enabled

    def apply(self, config):
        config.ocean_thermo.enable_thermodynamics = bool(self.enabled)


class SurfaceHeatFlux:
    """``SurfaceHeatFlux(q=0.0)`` -> `&ocean_thermo_nml q_heat` (W/m^2,
    positive down)."""

    def __init__(self, q=0.0):
        self.q = q

    def apply(self, config):
        config.ocean_thermo.q_heat = float(self.q)


class SurfaceSaltFlux:
    """``SurfaceSaltFlux(q=0.0)`` -> `q_salt` (kg/m^2/s, positive
    salinifies)."""

    def __init__(self, q=0.0):
        self.q = q

    def apply(self, config):
        config.ocean_thermo.q_salt = float(self.q)


class ShortwavePenetration:
    """``ShortwavePenetration(fraction=0.0, band_ratio=0.58, zeta1=0.35,
    zeta2=23.0, source="net_heat", epbl_tke_ledger=True)``

    -> `sw_pen_frac`, `sw_band_ratio`, `sw_zeta1`, `sw_zeta2`,
    `sw_source`, `epbl_sw_ctke`. All inert at `fraction=0`.
    """

    def __init__(self, fraction=0.0, band_ratio=0.58, zeta1=0.35,
                 zeta2=23.0, source="net_heat", epbl_tke_ledger=True):
        self.fraction, self.band_ratio = fraction, band_ratio
        self.zeta1, self.zeta2 = zeta1, zeta2
        self.source, self.epbl_tke_ledger = source, epbl_tke_ledger

    def apply(self, config):
        config.ocean_thermo.sw_pen_frac = float(self.fraction)
        config.ocean_thermo.sw_band_ratio = float(self.band_ratio)
        config.ocean_thermo.sw_zeta1 = float(self.zeta1)
        config.ocean_thermo.sw_zeta2 = float(self.zeta2)
        config.ocean_thermo.sw_source = self.source
        config.ocean_thermo.epbl_sw_ctke = bool(self.epbl_tke_ledger)


class BuoyancyRestoring:
    """``BuoyancyRestoring(sst=None, sss=None, piston_t=0.0,
    piston_s=0.0)`` -> `&ocean_restore_nml`.

    Passing `sst`/`sss` without a matching non-zero piston raises: the
    Fortran would silently restore with zero piston velocity (a no-op
    that looks configured).
    """

    def __init__(self, sst=None, sss=None, piston_t=0.0, piston_s=0.0):
        if sst is not None and piston_t == 0.0:
            raise ValueError(
                "BuoyancyRestoring(sst=...) needs a non-zero piston_t -- "
                "the Fortran would silently restore at zero rate")
        if sss is not None and piston_s == 0.0:
            raise ValueError(
                "BuoyancyRestoring(sss=...) needs a non-zero piston_s -- "
                "the Fortran would silently restore at zero rate")
        self.sst, self.sss = sst, sss
        self.piston_t, self.piston_s = piston_t, piston_s

    def apply(self, config):
        if self.sst is not None:
            config.ocean_restore.enable_restore_temp = True
            config.ocean_restore.restore_sst = float(self.sst)
            config.ocean_restore.piston_t = float(self.piston_t)
        if self.sss is not None:
            config.ocean_restore.enable_restore_salt = True
            config.ocean_restore.restore_sss = float(self.sss)
            config.ocean_restore.piston_s = float(self.piston_s)


class Geothermal:
    """``Geothermal(q=0.0)`` -> `&ocean_geothermal_nml enable` +
    `q_geo`."""

    def __init__(self, q=0.0):
        self.q = q

    def apply(self, config):
        config.ocean_geothermal.enable = True
        config.ocean_geothermal.q_geo = float(self.q)


class AtmosphericPressure:
    """``AtmosphericPressure(p_const=0.0)`` -> `&ocean_psurf_nml enable`
    + `p_surf_const`. Requires `SurfaceFluxComponents()` in the same
    `forcing=` list (checked at compose time)."""

    def __init__(self, p_const=0.0):
        self.p_const = p_const

    def apply(self, config):
        config.ocean_psurf.enable = True
        config.ocean_psurf.p_surf_const = float(self.p_const)


class SurfaceFluxComponents:
    """``SurfaceFluxComponents()`` -> `&ocean_forcing_nml
    enable_components`.

    Turning this on makes `Q_heat`/`Q_salt` DERIVED: a direct
    `model.Q_heat[:] = ...` write would be overwritten at the next
    thermo step. `Model` records this knob and the `Q_heat`/`Q_salt`
    setters then refuse with that reason, rather than writing into the
    void.
    """

    def apply(self, config):
        config.ocean_forcing.enable_components = True


class FileVar:
    """One `&ocean_dataovr_nml` file/variable pair.
    ``FileVar(path, var, scale=1.0, add=0.0)``."""

    def __init__(self, path, var, scale=1.0, add=0.0):
        self.path, self.var, self.scale, self.add = path, var, scale, add


def _as_filevar(value):
    if value is None:
        return None
    if isinstance(value, FileVar):
        return value
    path, var = value
    return FileVar(path, var)


class FileForcing:
    """``FileForcing(tau_x=None, tau_y=None, heat=None, evap=None,
    lprec=None, salt=None, time_mode="linear", cycle_period=0.0,
    t_offset=0.0, clamp=False)`` -> `&ocean_dataovr_nml`.

    Each argument is a `(path, var)` pair or a `FileVar(path, var,
    scale=1.0, add=0.0)`. This is the ERA5 keystone and also the
    single largest remaining `error stop` surface in the Fortran (20
    NetCDF-backed sites, P0 review F5) -- this object pre-flights path
    existence before `create()` ever runs.
    """

    _FIELDS = ("tau_x", "tau_y", "heat", "evap", "lprec", "salt")

    def __init__(self, tau_x=None, tau_y=None, heat=None, evap=None,
                 lprec=None, salt=None, time_mode="linear",
                 cycle_period=0.0, t_offset=0.0, clamp=False):
        _args = {"tau_x": tau_x, "tau_y": tau_y, "heat": heat,
                  "evap": evap, "lprec": lprec, "salt": salt}
        self.vars = {name: _as_filevar(_args[name])
                     for name in self._FIELDS}
        for name, fv in self.vars.items():
            if fv is not None and not Path(fv.path).is_file():
                raise RdbIOError(
                    f"FileForcing: {name} file not found: {fv.path}")
        self.time_mode = time_mode
        self.cycle_period, self.t_offset, self.clamp = (
            cycle_period, t_offset, clamp)

    def apply(self, config):
        config.ocean_dataovr.enable = True
        config.ocean_dataovr.time_mode = self.time_mode
        config.ocean_dataovr.cycle_period = float(self.cycle_period)
        config.ocean_dataovr.t_offset = float(self.t_offset)
        config.ocean_dataovr.oor_clamp = bool(self.clamp)
        for name, fv in self.vars.items():
            if fv is None:
                continue
            setattr(config.ocean_dataovr, f"{name}_file", str(fv.path))
            setattr(config.ocean_dataovr, f"{name}_var", fv.var)
            setattr(config.ocean_dataovr, f"{name}_scale", float(fv.scale))
            setattr(config.ocean_dataovr, f"{name}_add", float(fv.add))


# ======================================================================
# Tides
# ======================================================================


class ScalarSAL:
    """``ScalarSAL(beta=0.09)`` MODIFIER, passed as
    `EquilibriumTide(sal=...)`."""

    def __init__(self, beta=0.09):
        self.beta = beta


class EquilibriumTide:
    """``EquilibriumTide(constituents="M2 S2 N2 K2 K1 O1 P1 Q1",
    ref_date="1900-01-01", nodal=False, sal=None)`` -> `&ocean_tides_nml`.

    `sal=ScalarSAL(beta=...)` -> `use_sal` + `beta_sal`. Requires a
    non-cartesian grid (checked at compose time).
    """

    def __init__(self, constituents="M2 S2 N2 K2 K1 O1 P1 Q1",
                 ref_date="1900-01-01", nodal=False, sal=None):
        self.constituents, self.ref_date = constituents, ref_date
        self.nodal, self.sal = nodal, sal

    def apply(self, config):
        config.ocean_tides.enable = True
        config.ocean_tides.constituents = self.constituents
        config.ocean_tides.ref_date = self.ref_date
        config.ocean_tides.add_nodal = bool(self.nodal)
        if self.sal is not None:
            config.ocean_tides.use_sal = True
            config.ocean_tides.beta_sal = float(self.sal.beta)


# ======================================================================
# Initial-condition overlays
# ======================================================================


class EadyIC:
    """``EadyIC(dT_dy=-2e-5, dT_dz=1e-2, T_ref=10.0, perturbation=1e-3,
    seed=12345)`` -> `ic_config="eady"` + `eady_*`."""

    def __init__(self, dT_dy=-2e-5, dT_dz=1e-2, T_ref=10.0,
                 perturbation=1e-3, seed=12345):
        self.dT_dy, self.dT_dz, self.T_ref = dT_dy, dT_dz, T_ref
        self.perturbation, self.seed = perturbation, seed

    def apply(self, config):
        config.ocean_ic.ic_config = "eady"
        config.ocean_ic.eady_dT_dy = float(self.dT_dy)
        config.ocean_ic.eady_dT_dz = float(self.dT_dz)
        config.ocean_ic.eady_T_ref = float(self.T_ref)
        config.ocean_ic.eady_pert_amp = float(self.perturbation)
        config.ocean_ic.eady_pert_seed = int(self.seed)


class GeostrophicAdjustmentIC:
    """``GeostrophicAdjustmentIC(amplitude=1.0, length_scale=5e4,
    centre=None)`` -> `ic_config="geostrophic_adjustment"` + `ga_*`."""

    def __init__(self, amplitude=1.0, length_scale=5e4, centre=None):
        self.amplitude, self.length_scale, self.centre = (
            amplitude, length_scale, centre)

    def apply(self, config):
        config.ocean_ic.ic_config = "geostrophic_adjustment"
        config.ocean_ic.ga_eta_amp = float(self.amplitude)
        config.ocean_ic.ga_length_scale = float(self.length_scale)
        if self.centre is not None:
            config.ocean_ic.ga_x_center = float(self.centre[0])
            config.ocean_ic.ga_y_center = float(self.centre[1])


class BaroclinicJetIC:
    """``BaroclinicJetIC(half_width=4e4, interface_amp=200.0,
    meander_fraction=0.2, wavenumber=3, upper_thickness=500.0)`` ->
    `ic_config="baroclinic_jet"`. `nz` must be 2 (checked at compose
    time -- the Fortran rejects otherwise)."""

    def __init__(self, half_width=4e4, interface_amp=200.0,
                 meander_fraction=0.2, wavenumber=3,
                 upper_thickness=500.0):
        self.half_width, self.interface_amp = half_width, interface_amp
        self.meander_fraction, self.wavenumber = (
            meander_fraction, wavenumber)
        self.upper_thickness = upper_thickness

    def apply(self, config):
        config.ocean_ic.ic_config = "baroclinic_jet"
        config.ocean_ic.jet_half_width = float(self.half_width)
        config.ocean_ic.interface_amp = float(self.interface_amp)
        config.ocean_ic.pert_amp_frac = float(self.meander_fraction)
        config.ocean_ic.pert_nx = int(self.wavenumber)
        config.ocean_ic.upper_layer_rest = float(self.upper_thickness)


class LayerDensities:
    """``LayerDensities([...])`` -> `&ocean_ic_nml layer_rho_init`.

    `k=1` BED first, `k=nz` surface last -- the same bottom-up order as
    every 3-D field. This object does NOT reverse.
    """

    def __init__(self, densities):
        self.densities = list(densities)

    def apply(self, config):
        config.ocean_ic.layer_rho_init = [float(v) for v in self.densities]


class DensityRange:
    """``DensityRange(lightest=1020.0, contrast=2.0)`` -> `rho_lightest`,
    `rho_range`."""

    def __init__(self, lightest=1020.0, contrast=2.0):
        self.lightest, self.contrast = lightest, contrast

    def apply(self, config):
        config.ocean_ic.rho_lightest = float(self.lightest)
        config.ocean_ic.rho_range = float(self.contrast)


class ZLevelInitialCondition:
    """``ZLevelInitialCondition(file, t_var=None, s_var=None,
    z_var=None, land_fill=(10.0, 35.0))`` -> `&ocean_zinit_nml`.
    Requires `RDB_ENABLE_NETCDF`."""

    def __init__(self, file, t_var=None, s_var=None, z_var=None,
                 land_fill=(10.0, 35.0)):
        p = Path(file)
        if not p.is_file():
            raise RdbIOError(f"ZLevelInitialCondition: file not found: "
                                 f"{file}")
        self.file = str(p)
        self.t_var, self.s_var, self.z_var = t_var, s_var, z_var
        self.land_fill = land_fill

    def apply(self, config):
        config.ocean_zinit.enable = True
        config.ocean_zinit.file = self.file
        if self.t_var is not None:
            config.ocean_zinit.t_var = self.t_var
        if self.s_var is not None:
            config.ocean_zinit.s_var = self.s_var
        if self.z_var is not None:
            config.ocean_zinit.z_var = self.z_var
        config.ocean_zinit.land_fill_t = float(self.land_fill[0])
        config.ocean_zinit.land_fill_s = float(self.land_fill[1])


class UniformInitialCondition:
    """``UniformInitialCondition(temperature=15.0, salinity=35.0)`` ->
    `&tracer_nml initial_temperature`/`initial_salinity`."""

    def __init__(self, temperature=15.0, salinity=35.0):
        self.temperature, self.salinity = temperature, salinity

    def apply(self, config):
        config.tracer.initial_temperature = float(self.temperature)
        config.tracer.initial_salinity = float(self.salinity)


class StratifiedInitialCondition:
    """``StratifiedInitialCondition(T_surface=None, T_bottom=None)`` ->
    `&tracer_nml T_init_surface`/`T_init_bottom` -- linear T(z), bed
    (k=1) at `T_bottom`, surface (k=nz) at `T_surface`.

    The Fortran only builds this linear profile when BOTH are non-zero
    AND `nz > 1` (`rdb_ocean_state.F90:1204`'s `stratify` condition);
    otherwise it falls back to the UNIFORM `tracer.initial_temperature`
    (`UniformInitialCondition`). Passing `T_surface == T_bottom` (both
    non-zero) is the documented way to get a genuinely uniform T while
    still routing through the stratified code path (see
    `acc_channel_quiescent.py`).

    There is NO curated route for the analogous SALINITY profile
    (`S_init_surface`/`S_init_bottom`) -- unlike their temperature
    siblings, the P4 dead-knob sweep found these are accepted and
    validated but READ NOWHERE on the ocean path (2026-09-10). Use
    `UniformInitialCondition(salinity=...)` for a uniform S; there is
    no ocean-path stratified-salinity IC today.
    """

    def __init__(self, T_surface=None, T_bottom=None):
        if (T_surface is None) != (T_bottom is None):
            raise ValueError(
                "StratifiedInitialCondition needs T_surface= and "
                "T_bottom= TOGETHER -- the Fortran's stratify condition "
                "requires both non-zero; passing only one would "
                "silently fall back to the uniform initial_temperature")
        self.T_surface, self.T_bottom = T_surface, T_bottom

    def apply(self, config):
        if self.T_surface is not None:
            config.tracer.T_init_surface = float(self.T_surface)
            config.tracer.T_init_bottom = float(self.T_bottom)
