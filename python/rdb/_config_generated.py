"""Generated typed config layer -- DO NOT EDIT.

Emitted by `tools/gen_python_config.py` from the live `nml_schema_t`
(`src/core/rdb_config.F90` -> `build_rdb_schema`, walked via
`schema%render_json`, `src/core/rdb_nml_schema.F90`). One class per
namelist group; one descriptor attribute per knob, carrying the schema's
own doc/units/default/bounds/enum-set (see `python/rdb/_knob.py`).

ABSENCE IS THE DEFAULT, NOT THE FORTRAN VALUE: every knob starts unset,
and only explicitly-assigned knobs serialise via `to_namelist()`
(`python/rdb/_config.py`). This is what gives bit-identity to existing
namelists -- a config of all defaults serialises to nothing.

Regenerate with:

    ./build_shared/rdb_nml_json tmp_local_artifacts/schema.json
    python3 tools/gen_python_config.py tmp_local_artifacts/schema.json \
            python/rdb/_config_generated.py

`python/tests/test_config_drift.py` is the anti-drift gate: it runs this
pipeline into a temp file and diffs it against this checked-in file,
failing with the regeneration command above if they differ.
"""

from __future__ import annotations

from ._knob import Bool, Enum, Group, Int, Real, RealArray, Str

class Sim(Group):
    """`&sim_nml` -- Simulation-regime selector."""

    _nml_name = 'sim'

    sim_type = Enum(
        'sim_type',
        doc='Simulation regime (ocean is the only regime this build ships)',
        units='',
        required=False,
        default='ocean',
        allowed=('ocean',),
    )


class Grid(Group):
    """`&grid_nml` -- Structured-grid geometry."""

    _nml_name = 'grid'

    nx = Int(
        'nx',
        doc='Number of physical cells in x',
        units='',
        required=False,
        default=200,
        has_min=True,
        vmin=1,
    )

    ny = Int(
        'ny',
        doc='Number of physical cells in y',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )

    dx = Real('dx', doc='Cell size in x', units='m', required=False, default=1.0)

    dy = Real('dy', doc='Cell size in y', units='m', required=False, default=1.0)

    nghost = Int(
        'nghost',
        doc='Ghost cells on each side',
        units='',
        required=False,
        default=3,
        has_min=True,
        vmin=1,
    )


class Time(Group):
    """`&time_nml` -- Time-integration controls."""

    _nml_name = 'time'

    t_end = Real(
        't_end',
        doc='Simulation end time (in time_unit)',
        units='',
        required=False,
        default=1.0,
    )

    cfl = Real(
        'cfl',
        doc='CFL number for adaptive timestep',
        units='',
        required=False,
        default=0.45,
    )

    dt_max = Real(
        'dt_max',
        doc='Maximum allowable timestep',
        units='s',
        required=False,
        default=10000000000.0,
        dead_on_ocean_path="accepted and type/range-validated but read nowhere in src/ -- the ocean path's adaptive-timestep ceiling comes from the barotropic gravity-wave CFL (auto_n_inner), not this knob (found by the P4 dead-knob sweep, 2026-09-10).",
    )

    dt_fixed = Real(
        'dt_fixed',
        doc='Fixed timestep (0 = adaptive CFL)',
        units='s',
        required=False,
        default=0.0,
    )

    cfl_interval = Int(
        'cfl_interval',
        doc='Recompute CFL timestep every N steps',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
        dead_on_ocean_path='accepted and range-validated but read nowhere in src/ -- there is no adaptive-CFL recompute cadence on the ocean path to interval-gate (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    time_unit = Str(
        'time_unit',
        doc='Unit for the long-time fields: s/min/hr/day/year',
        units='',
        required=False,
        default='s',
        max_len=8,
    )


class Mpi(Group):
    """`&mpi_nml` -- MPI domain decomposition."""

    _nml_name = 'mpi'

    px = Int(
        'px',
        doc='MPI processes in x',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )

    py = Int(
        'py',
        doc='MPI processes in y',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )


class Logging(Group):
    """`&logging_nml` -- Logger verbosity + status cadence."""

    _nml_name = 'logging'

    log_level = Enum(
        'log_level',
        doc='Log verbosity',
        units='',
        required=False,
        default='info',
        allowed=('debug', 'verbose', 'info', 'performance', 'warning', 'error'),
    )

    status_interval = Real(
        'status_interval',
        doc='Status-print cadence (in time_unit; 0 = every 100 steps)',
        units='',
        required=False,
        default=0.0,
    )


class Vcoord(Group):
    """`&vcoord_nml` -- Vertical-coordinate + ALE-remap controls."""

    _nml_name = 'vcoord'

    vcoord_type = Enum(
        'vcoord_type',
        doc='Vertical coordinate type',
        units='',
        required=False,
        default='sigma',
        allowed=('sigma', 'zsigma', 'zstar', 'zstar_full', 'zstar_sigma', 'z', 'eulerian_z', 'isopycnal', 'lagrangian', 'gprime', 'z_fixed', 'rho', 'hycom'),
    )

    thickness_config = Enum(
        'thickness_config',
        doc='Initial layer-thickness profile (ocean path)',
        units='',
        required=False,
        default='sigma',
        allowed=('sigma', 'uniform_z'),
    )

    remap_method = Enum(
        'remap_method',
        doc='Vertical remapping method',
        units='',
        required=False,
        default='ppm',
        allowed=('pcm', 'plm', 'ppm', 'ppm_h4', 'pqm'),
    )

    zstar_h_surf_target = Real(
        'zstar_h_surf_target',
        doc='z*-full: target surface-layer thickness (0 = auto)',
        units='m',
        required=False,
        default=0.0,
    )

    zstar_stretching = Enum(
        'zstar_stretching',
        doc='z*-full: surface-concentration stretching',
        units='',
        required=False,
        default='log',
        allowed=('log', 'uniform'),
        dead_on_ocean_path='registered but never copied onto ocean_state%vcoord -- the ocean path always runs STRETCH_UNIFORM regardless of this setting (rdb_ocean_vcoord.F90:179). See docs/ocean_python_api_plan.md P4.5.',
    )

    zstar_h_min = Real(
        'zstar_h_min',
        doc='z*-full: vanishing-layer floor',
        units='m',
        required=False,
        default=0.0001,
    )

    zstar_n_surf = Int(
        'zstar_n_surf',
        doc='z*-full: number of fine near-surface layers (0 = auto)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
        dead_on_ocean_path='registered but never copied onto ocean_state%vcoord -- the ocean path always runs n_surf=0 regardless of this setting (rdb_ocean_vcoord.F90:182). See docs/ocean_python_api_plan.md P4.5.',
    )

    rho_ref_pressure = Real(
        'rho_ref_pressure',
        doc='rho-coord: reference pressure for potential density',
        units='Pa',
        required=False,
        default=20000000.0,
    )

    rho_target_light = Real(
        'rho_target_light',
        doc='rho-coord: lightest (surface) target density',
        units='kg/m^3',
        required=False,
        default=1020.0,
    )

    rho_target_dense = Real(
        'rho_target_dense',
        doc='rho-coord: densest (bed) target density',
        units='kg/m^3',
        required=False,
        default=1030.0,
    )

    regrid_time_scale = Real(
        'regrid_time_scale',
        doc='ALE regrid grid time-filter timescale (0 = jump to target)',
        units='s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    remap_vel_conserve_ke = Bool(
        'remap_vel_conserve_ke',
        doc='ALE velocity remap: KE-conserving baroclinic-anomaly rescale',
        units='',
        required=False,
        default=False,
    )

    zfixed_closed_faces = Bool(
        'zfixed_closed_faces',
        doc='z_fixed partial steps: close every face whose layer is an inert filler on either side (z-level wall, free-slip)',
        units='',
        required=False,
        default=False,
    )


class Physics(Group):
    """`&physics_nml` -- Barotropic physics: bottom drag, wind stress, Coriolis."""

    _nml_name = 'physics'

    manning_n = Real(
        'manning_n',
        doc='Manning roughness coefficient',
        units='',
        required=False,
        default=0.0,
        dead_on_ocean_path="coastal-legacy A-grid bottom-drag knob; the A-grid god state that read it is gone (see rdb_state.F90's module docstring). The ocean path's drag is &ocean_bdrag_nml.",
    )

    wind_stress_x = Real(
        'wind_stress_x',
        doc='Surface wind stress in x',
        units='Pa',
        required=False,
        default=0.0,
    )

    wind_stress_y = Real(
        'wind_stress_y',
        doc='Surface wind stress in y',
        units='Pa',
        required=False,
        default=0.0,
    )

    coriolis_f = Real(
        'coriolis_f',
        doc='Coriolis parameter f',
        units='1/s',
        required=False,
        default=0.0,
    )


class Output(Group):
    """`&output_nml` -- File output, I/O server, forcing/restart/gauge files."""

    _nml_name = 'output'

    output_to_file = Bool(
        'output_to_file',
        doc='Enable file output (NetCDF snapshots)',
        units='',
        required=False,
        default=False,
        dead_on_ocean_path='accepted and validated but read nowhere in src/ -- NetCDF output is gated by RDB_ENABLE_NETCDF at build time and by the diag/restart registries at run time, not by this flag (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    output_dir = Str(
        'output_dir',
        doc='Directory for output files',
        units='',
        required=False,
        default='./output',
        max_len=256,
    )

    restart_interval = Real(
        'restart_interval',
        doc='Time between restart writes (0 = none)',
        units='s',
        required=False,
        default=0.0,
    )

    compress_output = Bool(
        'compress_output',
        doc='Deflate-compress NetCDF output',
        units='',
        required=False,
        default=False,
    )

    compress_level = Int(
        'compress_level',
        doc='Deflate level (1=fast, 9=max)',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
        has_max=True,
        vmax=9,
    )

    use_io_server = Bool(
        'use_io_server',
        doc='Dedicate one MPI rank per node as I/O server',
        units='',
        required=False,
        default=False,
    )

    bathymetry_file = Str(
        'bathymetry_file',
        doc='NetCDF bathymetry file (empty = flat)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    restart_file = Str(
        'restart_file',
        doc='Restart file for warm start (empty = cold)',
        units='',
        required=False,
        default='',
        max_len=256,
    )


class Boundary(Group):
    """`&boundary_nml` -- Boundary conditions, tidal forcing, inflow/discharge/sponge/nesting."""

    _nml_name = 'boundary'

    bc_west = Enum(
        'bc_west',
        doc='West boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman'),
        dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc and nothing else. The ocean path's edge selector is &ocean_bc_nml west/east/south/north.",
    )

    bc_east = Enum(
        'bc_east',
        doc='East boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman'),
        dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc and nothing else. The ocean path's edge selector is &ocean_bc_nml west/east/south/north.",
    )

    bc_south = Enum(
        'bc_south',
        doc='South boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman'),
        dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc and nothing else. The ocean path's edge selector is &ocean_bc_nml west/east/south/north.",
    )

    bc_north = Enum(
        'bc_north',
        doc='North boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman'),
        dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc and nothing else. The ocean path's edge selector is &ocean_bc_nml west/east/south/north.",
    )

    n_tidal_constituents = Int(
        'n_tidal_constituents',
        doc='Active tidal constituents (0 = legacy single)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    tidal_amp = RealArray(
        'tidal_amp',
        doc='Constituent amplitudes',
        units='m',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=10,
    )

    tidal_phase = RealArray(
        'tidal_phase',
        doc='Constituent phases',
        units='rad',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=10,
    )

    tidal_omega = RealArray(
        'tidal_omega',
        doc='Constituent angular frequencies',
        units='rad/s',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=10,
    )

    inflow_salinity = Real(
        'inflow_salinity',
        doc='Inflow salinity (<0 = zero-gradient)',
        units='PSU',
        required=False,
        default=-1.0,
        dead_on_ocean_path='stored onto tracer_t%tr_inflow via register_default_tracers but tr_inflow is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    inflow_temperature = Real(
        'inflow_temperature',
        doc='Inflow temperature (<0 = zero-gradient)',
        units='degC',
        required=False,
        default=-999.0,
        dead_on_ocean_path='stored onto tracer_t%tr_inflow via register_default_tracers but tr_inflow is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    sponge_width = Int(
        'sponge_width',
        doc='Sponge layer width in cells',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    sponge_strength = Real(
        'sponge_strength',
        doc='Sponge relaxation rate',
        units='1/s',
        required=False,
        default=0.0,
    )


class Tracer(Group):
    """`&tracer_nml` -- Salinity + temperature IC/EOS/bounds, sediment transport."""

    _nml_name = 'tracer'

    initial_salinity = Real(
        'initial_salinity',
        doc='Initial salinity (uniform IC)',
        units='PSU',
        required=False,
        default=35.0,
    )

    S_ref = Real(
        'S_ref',
        doc='EOS reference salinity',
        units='PSU',
        required=False,
        default=0.0,
        dead_on_ocean_path='RETIRED -- stored onto tracer_t%eos_ref via register_default_tracers but eos_ref is read nowhere on the ocean path. The live ocean-path spelling is &ocean_ic_nml S_ref; setting THIS one to anything other than its default is a fail-loud configure error (validate_config).',
    )

    beta_S = Real(
        'beta_S',
        doc='Haline contraction coefficient',
        units='kg/m^3/PSU',
        required=False,
        default=0.78,
        dead_on_ocean_path='RETIRED -- stored onto tracer_t%eos_coeff via register_default_tracers but eos_coeff is read nowhere on the ocean path. The live ocean-path spelling is &ocean_ic_nml beta_S; setting THIS one to anything other than its default is a fail-loud configure error (validate_config).',
    )

    S_min = Real(
        'S_min',
        doc='Lower physical bound for salinity',
        units='PSU',
        required=False,
        default=0.0,
        dead_on_ocean_path='stored onto tracer_t%tr_min via register_default_tracers but tr_min is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    S_max = Real(
        'S_max',
        doc='Upper physical bound for salinity',
        units='PSU',
        required=False,
        default=40.0,
        dead_on_ocean_path='stored onto tracer_t%tr_max via register_default_tracers but tr_max is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    kappa_S_bg = Real(
        'kappa_S_bg',
        doc='Background vertical salinity diffusivity',
        units='m^2/s',
        required=False,
        default=1e-05,
        dead_on_ocean_path="stored onto tracer_t%kappa_bg via register_default_tracers but kappa_bg is read nowhere on the ocean path -- the ocean path's background salt diffusivity is &ocean_vmix_nml ks_bg (found by the P4 dead-knob sweep, 2026-09-10).",
    )

    S_init_surface = Real(
        'S_init_surface',
        doc='Initial surface salinity at k=nz (linear-in-layer stratified IC; needs S_init_bottom non-zero too)',
        units='PSU',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    S_init_bottom = Real(
        'S_init_bottom',
        doc='Initial bed salinity at k=1 (linear-in-layer stratified IC; needs S_init_surface non-zero too)',
        units='PSU',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    initial_temperature = Real(
        'initial_temperature',
        doc='Initial temperature (uniform IC)',
        units='degC',
        required=False,
        default=15.0,
    )

    T_ref = Real(
        'T_ref',
        doc='EOS reference temperature',
        units='degC',
        required=False,
        default=15.0,
        dead_on_ocean_path='RETIRED -- stored onto tracer_t%eos_ref via register_default_tracers but eos_ref is read nowhere on the ocean path. The live ocean-path spelling is &ocean_ic_nml T_ref; setting THIS one to anything other than its default is a fail-loud configure error (validate_config).',
    )

    alpha_T = Real(
        'alpha_T',
        doc='Thermal expansion coefficient',
        units='kg/m^3/degC',
        required=False,
        default=0.17,
        dead_on_ocean_path='RETIRED -- the coastal-legacy alpha_T (D2.5): stored onto tracer_t%eos_coeff via register_default_tracers but eos_coeff is read nowhere on the ocean path. The live ocean-path spelling is &ocean_ic_nml alpha_T; setting THIS one to anything other than its default is a fail-loud configure error (validate_config).',
    )

    T_min = Real(
        'T_min',
        doc='Lower physical bound for temperature',
        units='degC',
        required=False,
        default=-2.0,
        dead_on_ocean_path='stored onto tracer_t%tr_min via register_default_tracers but tr_min is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    T_max = Real(
        'T_max',
        doc='Upper physical bound for temperature',
        units='degC',
        required=False,
        default=40.0,
        dead_on_ocean_path='stored onto tracer_t%tr_max via register_default_tracers but tr_max is read nowhere on the ocean path (found by the P4 dead-knob sweep, 2026-09-10).',
    )

    kappa_T_bg = Real(
        'kappa_T_bg',
        doc='Background vertical temperature diffusivity',
        units='m^2/s',
        required=False,
        default=1e-05,
        dead_on_ocean_path="stored onto tracer_t%kappa_bg via register_default_tracers but kappa_bg is read nowhere on the ocean path -- the ocean path's background heat diffusivity is &ocean_vmix_nml kt_bg (found by the P4 dead-knob sweep, 2026-09-10).",
    )

    T_init_surface = Real(
        'T_init_surface',
        doc='Initial surface temperature (stratified IC)',
        units='degC',
        required=False,
        default=0.0,
    )

    T_init_bottom = Real(
        'T_init_bottom',
        doc='Initial bed temperature (stratified IC)',
        units='degC',
        required=False,
        default=0.0,
    )


class InitialCondition(Group):
    """`&initial_condition_nml` -- Coastal initial-condition selector + parameters."""

    _nml_name = 'initial_condition'

    h0 = Real(
        'h0',
        doc='Background depth (gaussian_hump IC)',
        units='m',
        required=False,
        default=1.0,
    )


class Nonhydrostatic(Group):
    """`&nonhydrostatic_nml` -- Non-hydrostatic / multilayer: CG-Poisson, PP81/KPP/k-eps vmix, Smagorinsky, BPG, mode split."""

    _nml_name = 'nonhydrostatic'

    nz_layers = Int(
        'nz_layers',
        doc='Number of sigma layers',
        units='',
        required=False,
        default=2,
        has_min=True,
        vmin=1,
    )

    use_multilayer = Bool(
        'use_multilayer',
        doc='Enable coupled hydrostatic vertical layers',
        units='',
        required=False,
        default=False,
    )

    rho_0 = Real(
        'rho_0',
        doc='Reference density for EOS',
        units='kg/m^3',
        required=False,
        default=1000.0,
    )

    kpp_ri_crit = Real(
        'kpp_ri_crit',
        doc='Critical bulk Richardson number for KPP BL-depth',
        units='',
        required=False,
        default=0.3,
    )

    kpp_cs_nonlocal = Real(
        'kpp_cs_nonlocal',
        doc='KPP non-local (counter-gradient) transport coefficient',
        units='',
        required=False,
        default=6.3,
    )

    kpp_c_vt2 = Real(
        'kpp_c_vt2',
        doc='KPP V_t^2 unresolved-turbulence coefficient (0 = off)',
        units='',
        required=False,
        default=1.8,
    )

    hdiff_kappa = Real(
        'hdiff_kappa',
        doc="Horizontal tracer diffusion coefficient (COASTAL path only; ocean's along-coordinate equivalent is &ocean_hdiff_nml kappa_h)",
        units='m^2/s',
        required=False,
        default=0.0,
        dead_on_ocean_path="coastal-legacy, per this key's own doc string: stored onto tracer_t%hdiff_kappa via register_default_tracers but hdiff_kappa is read nowhere on the ocean path. Use &ocean_hdiff_nml kappa_h.",
    )


class OceanKappaShear(Group):
    """`&ocean_kappa_shear_nml` -- JHL08 shear-driven interior turbulence (kappa-shear)."""

    _nml_name = 'ocean_kappa_shear'

    enable = Bool(
        'enable',
        doc='Master switch (requires use_closure + thermodynamics)',
        units='',
        required=False,
        default=False,
    )

    ri_crit = Real(
        'ri_crit',
        doc='Critical Richardson number',
        units='nondim',
        required=False,
        default=0.25,
    )

    shearmix_rate = Real(
        'shearmix_rate',
        doc='Shear source-rate coefficient',
        units='',
        required=False,
        default=0.089,
    )

    fri_curvature = Real(
        'fri_curvature',
        doc='Ri-function curvature in the shear source',
        units='',
        required=False,
        default=-0.97,
    )

    c_n = Real(
        'c_n',
        doc='TKE decay-rate coefficient vs stratification N',
        units='',
        required=False,
        default=0.24,
    )

    c_s = Real(
        'c_s',
        doc='TKE decay-rate coefficient vs shear S',
        units='',
        required=False,
        default=0.14,
    )

    lambda_ = Real(
        'lambda',
        doc='Buoyancy mixing-length-scale coefficient',
        units='',
        required=False,
        default=0.82,
    )

    lz_rescale = Real(
        'lz_rescale',
        doc='Boundary-distance length-scale rescale factor',
        units='',
        required=False,
        default=1.0,
    )

    kappa_0 = Real(
        'kappa_0',
        doc='Background diffusivity (pre-step kappa)',
        units='m^2/s',
        required=False,
        default=1e-07,
    )

    kappa_seed = Real(
        'kappa_seed',
        doc='Iteration seed diffusivity',
        units='m^2/s',
        required=False,
        default=1.0,
    )

    kappa_trunc = Real(
        'kappa_trunc',
        doc='Diffusivity truncated to 0 below this',
        units='m^2/s',
        required=False,
        default=1e-09,
    )

    tke_bg = Real(
        'tke_bg',
        doc='Background TKE (Q denominator floor)',
        units='m^2/s^2',
        required=False,
        default=0.0,
    )

    tol_err = Real(
        'tol_err',
        doc='Picard convergence tolerance',
        units='',
        required=False,
        default=0.1,
    )

    max_inner_it = Int(
        'max_inner_it',
        doc='Inner Picard iteration cap',
        units='',
        required=False,
        default=50,
        has_min=True,
        vmin=1,
    )

    max_substep_it = Int(
        'max_substep_it',
        doc='Outer adaptive-substep iteration cap',
        units='',
        required=False,
        default=13,
        has_min=True,
        vmin=1,
    )

    src_max_chg = Real(
        'src_max_chg',
        doc='Adaptive-dt source-change tolerance band',
        units='',
        required=False,
        default=10.0,
    )

    prandtl_turb = Real(
        'prandtl_turb',
        doc='Kv = prandtl_turb * Kd into the momentum solve',
        units='',
        required=False,
        default=1.0,
    )

    vel_underflow = Real(
        'vel_underflow',
        doc='Velocity snap-to-zero magnitude',
        units='m/s',
        required=False,
        default=0.0,
    )

    massless_merge = Bool(
        'massless_merge',
        doc='Merge vanished (<H_VANISHED) layers onto the massive sub-grid before the column solve (default off; identity columns bypass)',
        units='',
        required=False,
        default=False,
    )

    at_vertex = Bool(
        'at_vertex',
        doc='Solve the JHL08 columns at C-grid corners (vorticity points) from the native face velocities, averaging corner Kd back to tracer points (MOM6 VERTEX_SHEAR; v1 = Kd only)',
        units='',
        required=False,
        default=False,
    )

    vertex_geometric_mean = Bool(
        'vertex_geometric_mean',
        doc='Geometric (vs arithmetic) mean in the corner->centre Kd average (MOM6 VERTEX_SHEAR_GEOMETRIC_MEAN)',
        units='',
        required=False,
        default=False,
    )

    vertex_geomean_kdmin = Real(
        'vertex_geomean_kdmin',
        doc='Floor applied to each corner Kd before the geometric mean (inert unless vertex_geometric_mean; OM5 configs use 1e-9)',
        units='m^2/s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )


class OceanSlopes(Group):
    """`&ocean_slopes_nml` -- Isopycnal (neutral) slope diagnostics (Griffies 1998)."""

    _nml_name = 'ocean_slopes'

    enable = Bool(
        'enable',
        doc='Master switch (diagnostic; default off ⇒ no-op)',
        units='',
        required=False,
        default=False,
    )

    kd_smooth = Real(
        'kd_smooth',
        doc='Vert-fill smoothing diffusivity (× dt fills massless layers)',
        units='m^2/s',
        required=False,
        default=1e-06,
    )

    min_dz_for_n2 = Real(
        'min_dz_for_n2',
        doc='Minimum layer thickness floored in the N²/drdz denominator',
        units='m',
        required=False,
        default=1.0,
    )


class OceanGm(Group):
    """`&ocean_gm_nml` -- Gent-McWilliams thickness diffusion (eddy bolus transport)."""

    _nml_name = 'ocean_gm'

    enable = Bool(
        'enable',
        doc='Master switch (requires &ocean_slopes_nml enable; default off ⇒ no-op)',
        units='',
        required=False,
        default=False,
    )

    khth = Real(
        'khth',
        doc='Thickness diffusivity KhTh (constant-fill; production 1e2-1e3)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khth_max_cfl = Real(
        'khth_max_cfl',
        doc='Fraction of the diffusive CFL the face KH may use',
        units='nondim',
        required=False,
        default=0.1,
    )

    khth_slope_max = Real(
        'khth_slope_max',
        doc='Slope magnitude above which the safe-streamfunction blend takes over',
        units='nondim',
        required=False,
        default=0.01,
    )


class OceanRedi(Group):
    """`&ocean_redi_nml` -- Redi continuous neutral (along-isopycnal) tracer diffusion."""

    _nml_name = 'ocean_redi'

    enable = Bool(
        'enable',
        doc='Master switch (default off ⇒ no-op; augments hdiff_tracer)',
        units='',
        required=False,
        default=False,
    )

    continuous = Bool(
        'continuous',
        doc='Continuous variant (discontinuous deferred R4; .false. rejected)',
        units='',
        required=False,
        default=True,
    )

    khtr = Real(
        'khtr',
        doc='Redi neutral diffusivity KhTr (production 1e2-1e3)',
        units='m^2/s',
        required=False,
        default=0.0,
    )


class OceanVarmix(Group):
    """`&ocean_varmix_nml` -- Spatially-varying GM/Redi lateral-diffusivity coefficients."""

    _nml_name = 'ocean_varmix'

    enable = Bool(
        'enable',
        doc='Master switch (requires slopes + wavespeed; default off ⇒ GM uses const khth)',
        units='',
        required=False,
        default=False,
    )

    use_visbeck = Bool(
        'use_visbeck',
        doc='Add the Visbeck/Eady khth_slope_cff·L²·SN baroclinicity term',
        units='',
        required=False,
        default=False,
    )

    resoln_scaled_khth = Bool(
        'resoln_scaled_khth',
        doc='Scale the assembled KhTh by the resolution function',
        units='',
        required=False,
        default=False,
    )

    resoln_scaled_khtr = Bool(
        'resoln_scaled_khtr',
        doc='Scale the assembled KhTr by the resolution function',
        units='',
        required=False,
        default=False,
    )

    gill_equatorial_ld = Bool(
        'gill_equatorial_ld',
        doc='Gill (1982) equatorial-Ld convention (factor 2 in beta_dx2)',
        units='',
        required=False,
        default=True,
    )

    interpolate_res_fn = Bool(
        'interpolate_res_fn',
        doc='Interpolate centre Res_fn to faces (else interpolate cg1, MOM6 default)',
        units='',
        required=False,
        default=False,
    )

    kh_res_fn_power = Int(
        'kh_res_fn_power',
        doc='Resolution-function power p (even)',
        units='',
        required=False,
        default=2,
    )

    kh_res_scale_coef = Real(
        'kh_res_scale_coef',
        doc='Resolution-function alpha ((alpha·cg1)^p denominator coef)',
        units='nondim',
        required=False,
        default=1.0,
    )

    khth = Real(
        'khth',
        doc='Background thickness diffusivity KhTh',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khtr = Real(
        'khtr',
        doc='Background tracer diffusivity KhTr (future Redi)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khth_slope_cff = Real(
        'khth_slope_cff',
        doc='Visbeck coefficient for the KhTh chain',
        units='nondim',
        required=False,
        default=0.0,
    )

    khtr_slope_cff = Real(
        'khtr_slope_cff',
        doc='Visbeck coefficient for the KhTr chain',
        units='nondim',
        required=False,
        default=0.0,
    )

    khth_min = Real(
        'khth_min',
        doc='Lower clamp on KhTh',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khth_max = Real(
        'khth_max',
        doc='Upper clamp on KhTh (<= 0 ⇒ no cap)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khtr_min = Real(
        'khtr_min',
        doc='Lower clamp on KhTr',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    khtr_max = Real(
        'khtr_max',
        doc='Upper clamp on KhTr (<= 0 ⇒ no cap)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    visbeck_l_scale = Real(
        'visbeck_l_scale',
        doc='Visbeck length scale L (m); if < 0, |L|²·areaCu',
        units='m',
        required=False,
        default=0.0,
    )

    visbeck_max_slope = Real(
        'visbeck_max_slope',
        doc='S² limiter scale (<= 0 ⇒ no limit)',
        units='nondim',
        required=False,
        default=0.0,
    )


class OceanMeke(Group):
    """`&ocean_meke_nml` -- Prognostic mesoscale eddy kinetic energy (GM<->eddy loop)."""

    _nml_name = 'ocean_meke'

    enable = Bool(
        'enable',
        doc='Master switch (requires &ocean_gm_nml enable; default off ⇒ no-op)',
        units='',
        required=False,
        default=False,
    )

    gmcoeff = Real(
        'gmcoeff',
        doc='PE->MEKE conversion efficiency (< 0 ⇒ off)',
        units='nondim',
        required=False,
        default=-1.0,
    )

    frcoeff = Real(
        'frcoeff',
        doc='Frictional mean->eddy conversion (< 0 ⇒ off; >=0 sources hvisc KE dissipation)',
        units='nondim',
        required=False,
        default=-1.0,
    )

    bgsrc = Real(
        'bgsrc',
        doc='Background energy source',
        units='m^2/s^3',
        required=False,
        default=0.0,
    )

    damping = Real(
        'damping',
        doc='Linear MEKE dissipation rate',
        units='1/s',
        required=False,
        default=0.0,
    )

    kh = Real(
        'kh',
        doc='Background lateral diffusion of MEKE (< 0 ⇒ off)',
        units='m^2/s',
        required=False,
        default=-1.0,
    )

    k4 = Real(
        'k4',
        doc='Background biharmonic diffusion of MEKE (< 0 ⇒ off)',
        units='m^4/s',
        required=False,
        default=-1.0,
    )

    khcoeff = Real(
        'khcoeff',
        doc='MEKE->Kh scaling (<= 0 ⇒ closure off)',
        units='nondim',
        required=False,
        default=1.0,
    )

    cd_scale = Real(
        'cd_scale',
        doc='Bottom/column eddy-velocity ratio',
        units='nondim',
        required=False,
        default=0.0,
    )

    cb = Real(
        'cb',
        doc='Coefficient in gamma_bot (bottomFac2)',
        units='nondim',
        required=False,
        default=25.0,
    )

    ct = Real(
        'ct',
        doc='Coefficient in gamma_bt (barotrFac2)',
        units='nondim',
        required=False,
        default=50.0,
    )

    min_gamma2 = Real(
        'min_gamma2',
        doc='Floor on gamma_b^2/gamma_t^2',
        units='nondim',
        required=False,
        default=0.0001,
    )

    uscale = Real(
        'uscale',
        doc='Background eddy velocity scale for bottom drag',
        units='m/s',
        required=False,
        default=0.0,
    )

    dtscale = Real(
        'dtscale',
        doc='Time-stepping acceleration factor',
        units='nondim',
        required=False,
        default=1.0,
    )

    khth_fac = Real(
        'khth_fac',
        doc='Geom-mean kh -> VarMix KhTh factor (0 ⇒ inert)',
        units='nondim',
        required=False,
        default=0.0,
    )

    khtr_fac = Real(
        'khtr_fac',
        doc='Geom-mean kh -> VarMix KhTr factor (0 ⇒ inert)',
        units='nondim',
        required=False,
        default=0.0,
    )

    backscatter = Bool(
        'backscatter',
        doc='Enable MEKE -> momentum harmonic backscatter (negative viscosity); default off ⇒ bit-identical',
        units='',
        required=False,
        default=False,
    )

    backscatter_visc_coeff_ku = Real(
        'backscatter_visc_coeff_ku',
        doc='MEKE_VISCOSITY_COEFF_KU: harmonic backscatter efficiency Ku=coeff*sqrt(2*gt2*E)*Lmix (0 ⇒ inert)',
        units='nondim',
        required=False,
        default=0.0,
    )

    khmeke_fac = Real(
        'khmeke_fac',
        doc='meke%kh -> MEKE self-diffusivity factor',
        units='nondim',
        required=False,
        default=0.0,
    )

    advection_factor = Real(
        'advection_factor',
        doc='Barotropic-transport advection scaling (0 ⇒ off)',
        units='nondim',
        required=False,
        default=0.0,
    )

    cdrag = Real(
        'cdrag',
        doc='Bottom drag coefficient for MEKE',
        units='nondim',
        required=False,
        default=0.0025,
    )

    use_bbl_drag = Bool(
        'use_bbl_drag',
        doc='Add resolved |u_bed|^2 to the MEKE bottom-drag rate (default off ⇒ bit-identical)',
        units='',
        required=False,
        default=False,
    )

    alpha_deform = Real(
        'alpha_deform',
        doc='Weight on deformation length scale',
        units='nondim',
        required=False,
        default=0.0,
    )

    alpha_rhines = Real(
        'alpha_rhines',
        doc='Weight on Rhines length scale (v1 default 0 ⇒ inert)',
        units='nondim',
        required=False,
        default=0.0,
    )

    alpha_eady = Real(
        'alpha_eady',
        doc='Weight on Eady length scale (needs VarMix SN)',
        units='nondim',
        required=False,
        default=0.0,
    )

    alpha_frict = Real(
        'alpha_frict',
        doc='Weight on frictional-arrest length scale',
        units='nondim',
        required=False,
        default=0.0,
    )

    alpha_grid = Real(
        'alpha_grid',
        doc='Weight on grid length scale',
        units='nondim',
        required=False,
        default=0.0,
    )


class OceanTidalMixing(Group):
    """`&ocean_tidal_mixing_nml` -- St-Laurent/Simmons internal-tide interior diapycnal mixing."""

    _nml_name = 'ocean_tidal_mixing'

    enable = Bool(
        'enable',
        doc='Master switch (requires use_closure + thermodynamics)',
        units='',
        required=False,
        default=False,
    )

    gamma = Real(
        'gamma',
        doc='Local-dissipation fraction q (GAMMA_ITIDES)',
        units='nondim',
        required=False,
        default=0.3333,
    )

    mu = Real(
        'mu',
        doc='Mixing efficiency Gamma_mix (MU_ITIDES)',
        units='nondim',
        required=False,
        default=0.2,
    )

    zeta = Real(
        'zeta',
        doc='Bottom decay scale (INT_TIDE_DECAY_SCALE)',
        units='m',
        required=False,
        default=500.0,
    )

    kd_max = Real(
        'kd_max',
        doc='Per-layer physical Kd cap (<0 => no cap)',
        units='m^2/s',
        required=False,
        default=0.01,
    )

    prandtl_tidal = Real(
        'prandtl_tidal',
        doc='Kv = prandtl_tidal * Kd',
        units='',
        required=False,
        default=1.0,
    )

    min_zbot = Real(
        'min_zbot',
        doc='Mask off where column depth H < min_zbot',
        units='m',
        required=False,
        default=0.0,
    )

    e_uniform = Real(
        'e_uniform',
        doc='Uniform bottom internal-tide energy input E',
        units='W m-2',
        required=False,
        default=0.0,
    )

    e_compute = Bool(
        'e_compute',
        doc='State-dependent E = min(TKE_coef*N_bot, e_max) (v1.1)',
        units='',
        required=False,
        default=False,
    )

    kappa_itides = Real(
        'kappa_itides',
        doc='Topographic wavenumber (v1.1 E recompute)',
        units='m^-1',
        required=False,
        default=0.00062832,
    )

    kappa_h2 = Real(
        'kappa_h2',
        doc='KAPPA_H2_FACTOR (v1.1 E recompute)',
        units='',
        required=False,
        default=1.0,
    )

    utide = Real(
        'utide',
        doc='RMS barotropic tidal velocity (v1.1 E recompute)',
        units='m/s',
        required=False,
        default=0.0,
    )

    h2_rough = Real(
        'h2_rough',
        doc='Sub-grid topographic roughness variance <h^2>',
        units='m^2',
        required=False,
        default=0.0,
    )

    frac_rough = Real(
        'frac_rough',
        doc='Roughness clamp <h^2> <= (frac_rough*H)^2',
        units='',
        required=False,
        default=0.1,
    )

    e_max = Real(
        'e_max',
        doc='TKE_itide_max cap on E',
        units='W m-2',
        required=False,
        default=1000.0,
    )


class OceanConv(Group):
    """`&ocean_conv_nml` -- Brunt-Vaisala-triggered convective adjustment (interior closure contributor)."""

    _nml_name = 'ocean_conv'

    enable = Bool(
        'enable',
        doc='Master switch (requires use_closure + thermodynamics)',
        units='',
        required=False,
        default=False,
    )

    kd_conv = Real(
        'kd_conv',
        doc='Convective tracer diffusivity (KD_CONV)',
        units='m^2/s',
        required=False,
        default=1.0,
    )

    prandtl_conv = Real(
        'prandtl_conv',
        doc='Kv_conv = prandtl_conv * kd_conv',
        units='',
        required=False,
        default=1.0,
    )

    n2_thresh = Real(
        'n2_thresh',
        doc='Trigger threshold on N^2 (BV_SQR_CONV)',
        units='s^-2',
        required=False,
        default=0.0,
    )


class OceanPorous(Group):
    """`&ocean_porous_nml` -- Porous barriers: subgrid sill/strait blocking of the C-grid face widths."""

    _nml_name = 'ocean_porous'

    enable = Bool(
        'enable',
        doc='Master switch (ocean multilayer path only; fails loud with bt_halo > 0 or wetdry enable)',
        units='',
        required=False,
        default=False,
    )

    source = Enum(
        'source',
        doc="Along-face bathymetry-statistics source ('resolved' is a wet-gated proxy; 'file' is deferred)",
        units='',
        required=False,
        default='resolved',
        allowed=('resolved', 'file'),
    )

    eta_interp = Enum(
        'eta_interp',
        doc="Interface height at the velocity point (PORBAR_ETA_INTERP); 'max' is the LEAST blocking, 'min' the most",
        units='',
        required=False,
        default='max',
        allowed=('max', 'min', 'arithmetic', 'harmonic'),
    )

    masking_depth = Real(
        'masking_depth',
        doc='Faces shallower than this stay fully open (PORBAR_MASKING_DEPTH, positive below the surface)',
        units='m',
        required=False,
        default=0.0,
    )


class OceanDdiff(Group):
    """`&ocean_ddiff_nml` -- Double diffusion: salt fingering + diffusive convection (folded into the heat/salt split)."""

    _nml_name = 'ocean_ddiff'

    enable = Bool(
        'enable',
        doc='Master switch (requires thermodynamics)',
        units='',
        required=False,
        default=False,
    )

    strat_param_max = Real(
        'strat_param_max',
        doc='R_rho salt-fingering cutoff (STRAT_PARAM_MAX)',
        units='',
        required=False,
        default=2.55,
    )

    kappa_ddiff_s = Real(
        'kappa_ddiff_s',
        doc='Leading salt-fingering diffusivity K_f (KAPPA_DDIFF_S)',
        units='m^2/s',
        required=False,
        default=0.0001,
    )

    ddiff_exp1 = Real(
        'ddiff_exp1',
        doc='Inner fingering exponent (DDIFF_EXP1)',
        units='',
        required=False,
        default=1.0,
    )

    ddiff_exp2 = Real(
        'ddiff_exp2',
        doc='Outer fingering exponent (DDIFF_EXP2)',
        units='',
        required=False,
        default=3.0,
    )

    param1 = Real(
        'param1',
        doc='MC76 convection exterior coeff (KAPPA_DDIFF_PARAM1)',
        units='',
        required=False,
        default=0.909,
    )

    param2 = Real(
        'param2',
        doc='MC76 convection middle coeff (KAPPA_DDIFF_PARAM2)',
        units='',
        required=False,
        default=4.6,
    )

    param3 = Real(
        'param3',
        doc='MC76 convection interior coeff (KAPPA_DDIFF_PARAM3)',
        units='',
        required=False,
        default=-0.54,
    )

    mol_diff = Real(
        'mol_diff',
        doc='Molecular diffusivity scaling convection (MOL_DIFF)',
        units='m^2/s',
        required=False,
        default=1.5e-06,
    )

    use_k90 = Bool(
        'use_k90',
        doc='Convection form: MC76 (default) vs Kelley-90',
        units='',
        required=False,
        default=False,
    )


class OceanTides(Group):
    """`&ocean_tides_nml` -- Equilibrium (astronomical) body-force tide (C1) + scalar SAL (C2)."""

    _nml_name = 'ocean_tides'

    enable = Bool(
        'enable',
        doc='Master switch (requires a non-cartesian grid)',
        units='',
        required=False,
        default=False,
    )

    use_sal = Bool(
        'use_sal',
        doc='Apply scalar self-attraction & loading (C2)',
        units='',
        required=False,
        default=False,
    )

    beta_sal = Real(
        'beta_sal',
        doc='Scalar SAL factor beta (~0.085-0.12)',
        units='',
        required=False,
        default=0.0,
    )

    add_nodal = Bool(
        'add_nodal',
        doc='Apply the 18.6-yr nodal f/u corrections',
        units='',
        required=False,
        default=False,
    )

    constituents = Str(
        'constituents',
        doc="Active constituent list (e.g. 'M2 S2 K1 O1')",
        units='',
        required=False,
        default='M2 S2 N2 K2 K1 O1 P1 Q1',
        max_len=64,
    )

    ref_date = Str(
        'ref_date',
        doc='Astronomical reference date YYYY-MM-DD',
        units='',
        required=False,
        default='1900-01-01',
        max_len=16,
    )

    nodal_ref_date = Str(
        'nodal_ref_date',
        doc="Nodal reference date ('' => ref_date)",
        units='',
        required=False,
        default='',
        max_len=16,
    )


class OceanPsurf(Group):
    """`&ocean_psurf_nml` -- Atmospheric surface-pressure loading / inverse barometer (PR-17)."""

    _nml_name = 'ocean_psurf'

    enable = Bool(
        'enable',
        doc='Master switch (split-solver only; requires &ocean_forcing_nml enable_components=.true.)',
        units='',
        required=False,
        default=False,
    )

    in_eos = Bool(
        'in_eos',
        doc='Also feed the surface load to the equation of state as the top-of-column pressure p_top (requires enable=.true.; refused with the unported pressure builders)',
        units='',
        required=False,
        default=False,
    )

    p_surf_const = Real(
        'p_surf_const',
        doc='Uniform atmospheric surface pressure (Pa) seeded into p_surf_atm (uniform => provably inert)',
        units='',
        required=False,
        default=0.0,
    )


class OceanCavityDyn(Group):
    """`&ocean_cavity_dyn_nml` -- Static ice-shelf cavity geometry: prescribed draft + barotropic datum bt_H_ref = b - z_draft."""

    _nml_name = 'ocean_cavity_dyn'

    enable = Bool(
        'enable',
        doc='Master switch (single-rank, split solver, fv_mom6 PGF, sigma/zstar only)',
        units='',
        required=False,
        default=False,
    )

    draft_config = Enum(
        'draft_config',
        doc="Draft source: analytic shape, or 'file' (static 2-D NetCDF on the model grid, single rank)",
        units='',
        required=False,
        default='none',
        allowed=('none', 'flat', 'linear', 'file'),
    )

    draft_source = Enum(
        'draft_source',
        doc="Whether the formula gives the ice-base DEPTH or an ice THICKNESS ('in_situ' isostasy is deferred)",
        units='',
        required=False,
        default='draft',
        allowed=('draft', 'thickness', 'in_situ'),
    )

    draft_depth = Real(
        'draft_depth',
        doc="Draft amplitude (ice thickness under draft_source='thickness')",
        units='m',
        required=False,
        default=0.0,
    )

    draft_slope = Real(
        'draft_slope',
        doc="d(draft)/dx for draft_config='linear' (dimensionless; converted to grid units)",
        units='',
        required=False,
        default=0.0,
    )

    draft_x0 = Real(
        'draft_x0',
        doc="Western edge of the shelf box, and the anchor of the 'linear' profile (+/-1e30 => no limit)",
        units='m',
        required=False,
        default=-1e+30,
    )

    draft_x1 = Real(
        'draft_x1',
        doc='Eastern edge of the shelf box = the calving front (+/-1e30 => no limit)',
        units='m',
        required=False,
        default=1e+30,
    )

    draft_y0 = Real(
        'draft_y0',
        doc='Southern edge of the shelf box (+/-1e30 => no limit)',
        units='m',
        required=False,
        default=-1e+30,
    )

    draft_y1 = Real(
        'draft_y1',
        doc='Northern edge of the shelf box (+/-1e30 => no limit)',
        units='m',
        required=False,
        default=1e+30,
    )

    draft_file = Str(
        'draft_file',
        doc="draft_config='file': NetCDF path (variable must be (x,y,t) Fortran order, on the model grid; record 1 read)",
        units='',
        required=False,
        default='',
        max_len=256,
    )

    draft_var = Str(
        'draft_var',
        doc="draft_config='file': 2-D variable name (ISOMIP+ ships 'iceDraft')",
        units='',
        required=False,
        default='iceDraft',
        max_len=64,
    )

    draft_sign = Enum(
        'draft_sign',
        doc="draft_config='file': sign convention of the file values (ISOMIP+ iceDraft is an ELEVATION)",
        units='',
        required=False,
        default='depth',
        allowed=('depth', 'positive_down', 'elevation', 'positive_up'),
    )

    h_min_cavity = Real(
        'h_min_cavity',
        doc='Grounding cutoff: b - z_draft below this is LAND (never a thin film under grounded ice)',
        units='m',
        required=False,
        default=10.0,
    )

    grounded_max_frac = Real(
        'grounded_max_frac',
        doc='Fail loud if more than this fraction of the interior columns ground',
        units='',
        required=False,
        default=0.5,
    )

    rho_ice = Real(
        'rho_ice',
        doc="Ice density, consulted only by draft_source='thickness'",
        units='kg/m^3',
        required=False,
        default=918.0,
    )


class OceanCavityMelt(Group):
    """`&ocean_cavity_melt_nml` -- Ice-shelf basal-melt thermodynamics: the three-equation interface, its exchange law, and the far-field sampling depth."""

    _nml_name = 'ocean_cavity_melt'

    enable = Bool(
        'enable',
        doc="Master switch (requires &ocean_cavity_dyn_nml, tfreeze_set='isomip' and the surface-flux component set)",
        units='',
        required=False,
        default=False,
    )

    exchange_law = Enum(
        'exchange_law',
        doc='Turbulent exchange law; laws other than const_gamma/hj99/yung25 are RESERVED and refused at configure',
        units='',
        required=False,
        default='const_gamma',
        allowed=('const_gamma', 'hj99', 'yung25', 'jenkins91', 'rosevear22', 'vt19', 'mk18', 'burchard22', 'jenkins21'),
    )

    gamma_t = Real(
        'gamma_t',
        doc='Dimensionless heat-transfer coefficient Gamma_T (ISOMIP+ starting guess; tune per coordinate)',
        units='',
        required=False,
        default=0.022,
        has_min=True,
        vmin=0.0,
    )

    gamma_s = Real(
        'gamma_s',
        doc='Dimensionless salt-transfer coefficient Gamma_S (negative = unset = gamma_t/35)',
        units='',
        required=False,
        default=-1.0,
    )

    cdrag_top = Real(
        'cdrag_top',
        doc='Top drag coefficient for the MELT friction velocity (no momentum drag yet - that is Phase 4)',
        units='',
        required=False,
        default=0.0025,
        has_min=True,
        vmin=0.0,
    )

    u_tide = Real(
        'u_tide',
        doc='RMS tidal velocity in the melt u* only, never the drag',
        units='m/s',
        required=False,
        default=0.01,
        has_min=True,
        vmin=0.0,
    )

    ustar_min = Real(
        'ustar_min',
        doc='Friction-velocity floor (Yung et al. 2025 eq. 14)',
        units='m/s',
        required=False,
        default=0.0001,
        has_min=True,
        vmin=0.0,
    )

    ice_conduction = Enum(
        'ice_conduction',
        doc="Ice-side conduction; 'diffusive' is RESERVED and refused (it changes the melt/freeze branch logic)",
        units='',
        required=False,
        default='insulating',
        allowed=('insulating', 'adv_diff', 'diffusive'),
    )

    t_ice = Real(
        't_ice',
        doc="Ice interior temperature; read by ice_conduction='adv_diff' only",
        units='degC',
        required=False,
        default=-25.0,
    )

    s_ice = Real(
        's_ice',
        doc='Ice salinity; must stay strictly below the far-field salinity',
        units='g/kg',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    far_field_depth = Real(
        'far_field_depth',
        doc='Thickness below the ice base the far-field T/S/u are averaged over (METRES, not layers)',
        units='m',
        required=False,
        default=10.0,
        has_min=True,
        vmin=0.0,
    )

    freshwater = Enum(
        'freshwater',
        doc="Meltwater delivery: 'virtual' (default, fixed column mass) or 'mass' (real Boussinesq volume on the top layer)",
        units='',
        required=False,
        default='virtual',
        allowed=('virtual', 'mass'),
    )

    volume_compensation = Enum(
        'volume_compensation',
        doc="Sea-level compensation for freshwater='mass': 'none' (default) or 'uniform_open_ocean' (remove the melt volume again over uncovered wet cells)",
        units='',
        required=False,
        default='none',
        allowed=('none', 'uniform_open_ocean'),
    )


class OceanEpbl(Group):
    """`&ocean_epbl_nml` -- RH18 energetics-based planetary boundary layer."""

    _nml_name = 'ocean_epbl'

    enable = Bool(
        'enable',
        doc='Master switch (replaces the KPP overlay; requires use_closure)',
        units='',
        required=False,
        default=False,
    )

    mstar_scheme = Enum(
        'mstar_scheme',
        doc='Surface TKE mstar scheme',
        units='',
        required=False,
        default='om4',
        allowed=('constant', 'om4', 'rh18'),
    )

    mstar = Real('mstar', doc='Constant-scheme mstar', units='', required=False, default=1.2)

    mstar_cap = Real(
        'mstar_cap',
        doc='Cap for OM4/RH18 mstar; off when < 0',
        units='',
        required=False,
        default=-1.0,
    )

    mstar_coef1 = Real(
        'mstar_coef1',
        doc='OM4 stabilizing coefficient',
        units='',
        required=False,
        default=0.3,
    )

    c_ek = Real('c_ek', doc='OM4 Ekman coefficient', units='', required=False, default=0.085)

    mstar_conv_adj = Real(
        'mstar_conv_adj',
        doc='Convective mstar reduction in [0,1]',
        units='',
        required=False,
        default=0.0,
    )

    rh18_cn1 = Real(
        'rh18_cn1',
        doc='RH18 mstar fit coefficient cn1',
        units='',
        required=False,
        default=0.275,
    )

    rh18_cn2 = Real(
        'rh18_cn2',
        doc='RH18 mstar fit coefficient cn2',
        units='',
        required=False,
        default=8.0,
    )

    rh18_cn3 = Real(
        'rh18_cn3',
        doc='RH18 mstar fit coefficient cn3',
        units='',
        required=False,
        default=-5.0,
    )

    rh18_cs1 = Real(
        'rh18_cs1',
        doc='RH18 mstar fit coefficient cs1',
        units='',
        required=False,
        default=0.2,
    )

    rh18_cs2 = Real(
        'rh18_cs2',
        doc='RH18 mstar fit coefficient cs2',
        units='',
        required=False,
        default=0.4,
    )

    nstar = Real(
        'nstar',
        doc='Convective PE -> TKE efficiency',
        units='',
        required=False,
        default=0.2,
    )

    tke_decay = Real(
        'tke_decay',
        doc='Ekman-depth / TKE-decay-scale ratio',
        units='',
        required=False,
        default=2.5,
    )

    wstar_ustar_coef = Real(
        'wstar_ustar_coef',
        doc='Convective weight in the velocity scale',
        units='',
        required=False,
        default=1.0,
    )

    vel_scale_scheme = Enum(
        'vel_scale_scheme',
        doc='Velocity-scale scheme',
        units='',
        required=False,
        default='cube_root',
        allowed=('cube_root', 'rh18'),
    )

    vstar_scale_fac = Real(
        'vstar_scale_fac',
        doc='Overall vstar multiplier',
        units='',
        required=False,
        default=1.0,
    )

    vstar_surf_fac = Real(
        'vstar_surf_fac',
        doc='RH18 mechanical surface vstar factor',
        units='',
        required=False,
        default=1.2,
    )

    von_karman = Real(
        'von_karman',
        doc='von Karman kappa in Kd = vstar*kappa*mixlen',
        units='',
        required=False,
        default=0.41,
    )

    ekman_scale_coef = Real(
        'ekman_scale_coef',
        doc='Rotational mixing-length rolloff',
        units='',
        required=False,
        default=1.0,
    )

    min_mix_len = Real(
        'min_mix_len',
        doc='Mixing-length floor',
        units='m',
        required=False,
        default=0.0,
    )

    mixlen_exponent = Real(
        'mixlen_exponent',
        doc='Shape-function exponent',
        units='',
        required=False,
        default=2.0,
    )

    translay_scale = Real(
        'translay_scale',
        doc='Transition-layer shape floor (in [0,1) when iterating)',
        units='',
        required=False,
        default=0.1,
    )

    mld_iteration = Bool(
        'mld_iteration',
        doc='Self-consistent MLD root-find',
        units='',
        required=False,
        default=True,
    )

    mld_tol = Real(
        'mld_tol',
        doc='MLD convergence tolerance',
        units='m',
        required=False,
        default=1.0,
    )

    mld_max_its = Int(
        'mld_max_its',
        doc='Max MLD iterations',
        units='',
        required=False,
        default=20,
        has_min=True,
        vmin=1,
    )

    mld_bisection = Bool(
        'mld_bisection',
        doc='Bisection instead of false position',
        units='',
        required=False,
        default=False,
    )

    mld_use_prev_guess = Bool(
        'mld_use_prev_guess',
        doc="Seed from the previous step's MLD",
        units='',
        required=False,
        default=False,
    )

    omega = Real(
        'omega',
        doc='Earth rotation rate',
        units='1/s',
        required=False,
        default=7.2921e-05,
    )

    omega_frac = Real(
        'omega_frac',
        doc='Blend |f| with 2*Omega',
        units='',
        required=False,
        default=0.0,
    )

    prandtl = Real(
        'prandtl',
        doc='Kv = prandtl*Kd into the momentum solve',
        units='',
        required=False,
        default=1.0,
    )

    combine = Enum(
        'combine',
        doc='Combine vs interior closure kv/kt',
        units='',
        required=False,
        default='add',
        allowed=('add', 'max'),
    )

    tke_diags = Bool(
        'tke_diags',
        doc='Compute per-column TKE budget diagnostics',
        units='',
        required=False,
        default=False,
    )

    use_lt = Bool(
        'use_lt',
        doc='Langmuir-turbulence enhancement (LF17 wind-only)',
        units='',
        required=False,
        default=False,
    )

    lt_scheme = Enum(
        'lt_scheme',
        doc='Langmuir enhancement scheme',
        units='',
        required=False,
        default='rescale',
        allowed=('rescale', 'additive'),
    )

    lt_enhance_coef = Real(
        'lt_enhance_coef',
        doc='Langmuir enhancement coefficient',
        units='',
        required=False,
        default=0.447,
    )

    lt_enhance_exp = Real(
        'lt_enhance_exp',
        doc='Langmuir-number exponent',
        units='',
        required=False,
        default=-1.33,
    )

    lt_max_enhance = Real(
        'lt_max_enhance',
        doc='Cap on the multiplicative enhancement',
        units='',
        required=False,
        default=5.0,
    )

    la_frac_hbl = Real(
        'la_frac_hbl',
        doc='Stokes SL-average depth fraction',
        units='',
        required=False,
        default=0.04,
    )

    lt_lac1 = Real(
        'lt_lac1',
        doc='Stability-modified La coefficient 1',
        units='',
        required=False,
        default=-0.87,
    )

    lt_lac2 = Real(
        'lt_lac2',
        doc='Stability-modified La coefficient 2',
        units='',
        required=False,
        default=0.0,
    )

    lt_lac3 = Real(
        'lt_lac3',
        doc='Stability-modified La coefficient 3',
        units='',
        required=False,
        default=0.0,
    )

    lt_lac4 = Real(
        'lt_lac4',
        doc='Stability-modified La coefficient 4',
        units='',
        required=False,
        default=0.95,
    )

    lt_lac5 = Real(
        'lt_lac5',
        doc='Stability-modified La coefficient 5',
        units='',
        required=False,
        default=0.95,
    )


class OceanWavespeed(Group):
    """`&ocean_wavespeed_nml` -- First-baroclinic wave speed + Rossby deformation radius (diagnostic)."""

    _nml_name = 'ocean_wavespeed'

    enable = Bool(
        'enable',
        doc='Master switch (diagnostic; default off)',
        units='',
        required=False,
        default=False,
    )

    mono_n2 = Real(
        'mono_n2',
        doc='DEFERRED N2-monotonising depth (EBT path); < 0 = off',
        units='',
        required=False,
        default=-1.0,
    )

    use_ebt = Bool(
        'use_ebt',
        doc='DEFERRED equivalent-barotropic variant',
        units='',
        required=False,
        default=False,
    )

    n_wavespeed = Int(
        'n_wavespeed',
        doc='Recompute cadence (every N steps)',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )


class OceanFoxkemper(Group):
    """`&ocean_foxkemper_nml` -- Fox-Kemper mixed-layer-eddy restratification (B5)."""

    _nml_name = 'ocean_foxkemper'

    enable = Bool(
        'enable',
        doc='Master switch (off => bit-identity)',
        units='',
        required=False,
        default=False,
    )

    ce = Real(
        'ce',
        doc='FK08 coefficient Ce (0.06-0.08)',
        units='',
        required=False,
        default=0.0625,
    )

    f_floor = Real(
        'f_floor',
        doc='|f| regularisation floor',
        units='1/s',
        required=False,
        default=1e-05,
    )

    mld_decay_time = Real(
        'mld_decay_time',
        doc='Running-mean MLD filter time-scale (0 = off, instantaneous MLD)',
        units='s',
        required=False,
        default=0.0,
    )

    tail_dh = Real(
        'tail_dh',
        doc='mu cubic-tail extension (0 = exact mu)',
        units='',
        required=False,
        default=0.0,
    )

    use_mom_mixrate = Bool(
        'use_mom_mixrate',
        doc='FK11 momentum-mixrate timescale (PRODUCTION-RECOMMENDED) vs bare Ce/|f|',
        units='',
        required=False,
        default=False,
    )

    resolution_taper = Bool(
        'resolution_taper',
        doc='B2 res_fn hook (hard error if on without B2)',
        units='',
        required=False,
        default=False,
    )

    use_bodner = Bool(
        'use_bodner',
        doc='Bodner 2023 frontogenesis-arrest MLE (overrides ce/mixrate)',
        units='',
        required=False,
        default=False,
    )

    cr = Real(
        'cr',
        doc='Bodner 2023 efficiency coefficient Cr (0 = off)',
        units='',
        required=False,
        default=0.0,
    )

    bodner_mstar = Real(
        'bodner_mstar',
        doc="Bodner mechanical (u*) weight in w'u'",
        units='',
        required=False,
        default=0.5,
    )

    bodner_nstar = Real(
        'bodner_nstar',
        doc="Bodner convective (w*) weight in w'u'",
        units='',
        required=False,
        default=0.066,
    )

    min_wstar2 = Real(
        'min_wstar2',
        doc="Floor on w'u' (1/0 armour)",
        units='m^2/s^2',
        required=False,
        default=1e-24,
    )


class OceanGrid(Group):
    """`&ocean_grid_nml` -- Horizontal-grid generator + geometry (ocean curvilinear stream)."""

    _nml_name = 'ocean_grid'

    grid_config = Enum(
        'grid_config',
        doc='Horizontal grid generator',
        units='',
        required=False,
        default='cartesian',
        allowed=('cartesian', 'spherical', 'supergrid', 'tripolar'),
    )

    lon_west = Real(
        'lon_west',
        doc='West edge of the domain (spherical)',
        units='degrees_east',
        required=False,
        default=0.0,
    )

    lat_south = Real(
        'lat_south',
        doc='South edge of the domain (spherical)',
        units='degrees_north',
        required=False,
        default=0.0,
        has_min=True,
        vmin=-90.0,
        has_max=True,
        vmax=90.0,
    )

    rad_earth = Real(
        'rad_earth',
        doc='Earth radius for the spherical metric',
        units='m',
        required=False,
        default=6378000.0,
        has_min=True,
        vmin=0.0,
    )

    supergrid_file = Str(
        'supergrid_file',
        doc='Path to the MOM6 supergrid NetCDF (supergrid grid_config)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    coriolis_scheme = Enum(
        'coriolis_scheme',
        doc='Coriolis source for the metrics f-fill',
        units='',
        required=False,
        default='beta_plane',
        allowed=('beta_plane', 'planetary'),
    )

    omega = Real(
        'omega',
        doc='Planetary rotation rate (planetary scheme)',
        units='rad/s',
        required=False,
        default=7.2921e-05,
    )

    phi_join = Real(
        'phi_join',
        doc='Join latitude for the tripolar bipolar cap',
        units='degrees_north',
        required=False,
        default=65.0,
        has_min=True,
        vmin=-90.0,
        has_max=True,
        vmax=90.0,
    )

    lon_pole = Real(
        'lon_pole',
        doc='Longitude of the first tripolar cap pole (partner +180)',
        units='degrees_east',
        required=False,
        default=100.0,
    )

    axis_units = Enum(
        'axis_units',
        doc='Units of the Cartesian domain extent (MOM6 AXIS_UNITS)',
        units='',
        required=False,
        default='meters',
        allowed=('meters', 'degrees', 'km'),
    )

    len_lon = Real(
        'len_lon',
        doc='Total x-extent of the Cartesian domain in axis_units (MOM6 LENLON; derives dx when > 0)',
        units='',
        required=False,
        default=0.0,
    )

    len_lat = Real(
        'len_lat',
        doc='Total y-extent of the Cartesian domain in axis_units (MOM6 LENLAT; derives dy when > 0)',
        units='',
        required=False,
        default=0.0,
    )


class OceanCoriolis(Group):
    """`&ocean_coriolis_nml` -- Coriolis-advection scheme selector."""

    _nml_name = 'ocean_coriolis'

    form = Enum(
        'form',
        doc='Coriolis-advection variant',
        units='',
        required=False,
        default='sadourny',
        allowed=('sadourny', 'sadourny_hk', 'sadourny_energy'),
    )

    pv_adv_scheme = Enum(
        'pv_adv_scheme',
        doc='PV face interpolation (Sadourny path): centered (default) or weno3/weno5/weno7 (WENO-Z); weno5/weno7 need nghost>=3/4',
        units='',
        required=False,
        default='centered',
        allowed=('centered', 'weno3', 'weno5', 'weno7'),
    )

    use_state_fluxes = Bool(
        'use_state_fluxes',
        doc="mom6-corrector CorAdv consumes continuity's renormalised mass fluxes (MOM6 mass-consistent uh/vh)",
        units='',
        required=False,
        default=False,
    )

    bound_coriolis = Bool(
        'bound_coriolis',
        doc='clamp the energy-scheme Coriolis accel to the (f+zeta)*v velocity-form range (MOM6 BOUND_CORIOLIS; sadourny_energy only)',
        units='',
        required=False,
        default=False,
    )

    corner_h = Enum(
        'corner_h',
        doc='PV corner-thickness construction (energy scheme)',
        units='',
        required=False,
        default='cell_mean',
        allowed=('cell_mean', 'mom6_area'),
    )


class OceanThermo(Group):
    """`&ocean_thermo_nml` -- Thermodynamics master switch + scalar surface fluxes."""

    _nml_name = 'ocean_thermo'

    sw_source = Enum(
        'sw_source',
        doc='Shortwave irradiance source: net_heat (legacy) or q_sw',
        units='',
        required=False,
        default='net_heat',
        allowed=('net_heat', 'q_sw'),
    )

    kpp_sw_method = Enum(
        'kpp_sw_method',
        doc='KPP shortwave-in-BL method: all_sw | mxl_sw | lv1_sw',
        units='',
        required=False,
        default='mxl_sw',
        allowed=('all_sw', 'mxl_sw', 'lv1_sw'),
    )

    epbl_sw_ctke = Bool(
        'epbl_sw_ctke',
        doc='Charge the EPBL TKE ledger for penetrating shortwave',
        units='',
        required=False,
        default=True,
    )

    enable_thermodynamics = Bool(
        'enable_thermodynamics',
        doc='Run EOS + tracer advection + vertical mixing',
        units='',
        required=False,
        default=True,
    )

    q_heat = Real(
        'q_heat',
        doc='Net surface heat flux (positive down)',
        units='W/m^2',
        required=False,
        default=0.0,
    )

    q_salt = Real(
        'q_salt',
        doc='Net surface salt flux (positive salinifies)',
        units='kg/m^2/s',
        required=False,
        default=0.0,
    )

    sw_pen_frac = Real(
        'sw_pen_frac',
        doc='Penetrating fraction of q_heat (0 = off, all at surface)',
        units='',
        required=False,
        default=0.0,
    )

    sw_band_ratio = Real(
        'sw_band_ratio',
        doc='Two-band shortwave band-1 weight R (Jerlov type I)',
        units='',
        required=False,
        default=0.58,
    )

    sw_zeta1 = Real(
        'sw_zeta1',
        doc='Shortwave band-1 e-folding depth',
        units='m',
        required=False,
        default=0.35,
    )

    sw_zeta2 = Real(
        'sw_zeta2',
        doc='Shortwave band-2 e-folding depth',
        units='m',
        required=False,
        default=23.0,
    )


class OceanForcing(Group):
    """`&ocean_forcing_nml` -- Surface-flux component-set gate (PR-12): allocate the q_sw/evap/heat_content_*/... set and derive Q_heat/Q_salt from it every thermo step."""

    _nml_name = 'ocean_forcing'

    enable_components = Bool(
        'enable_components',
        doc='Allocate the surface-flux component set and run the assembler (default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )


class OceanIce(Group):
    """`&ocean_ice_nml` -- Sea-ice model (SIS2 port) master switch + category/layer counts."""

    _nml_name = 'ocean_ice'

    enable = Bool(
        'enable',
        doc='Master switch for the sea-ice slot (default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )

    ncat = Int(
        'ncat',
        doc='Number of ice thickness categories (cat 0 = open water)',
        units='',
        required=False,
        default=5,
        has_min=True,
        vmin=1,
    )

    hlim = RealArray(
        'hlim',
        doc='ITD category lower thickness edges; unset => SIS2 default table',
        units='m',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=16,
    )

    nk_ice = Int(
        'nk_ice',
        doc='Vertical ice layers per category (2 = Winton two-layer)',
        units='',
        required=False,
        default=2,
        has_min=True,
        vmin=1,
    )

    air_temp = Real(
        'air_temp',
        doc='Prescribed slab-atmosphere air temperature for the v1 restoring filler',
        units='degC',
        required=False,
        default=0.0,
    )

    restore_lambda = Real(
        'restore_lambda',
        doc='Surface-flux restoring coefficient dSF/dT (0 = passive column)',
        units='W/m^2/K',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    sw_down = Real(
        'sw_down',
        doc='Downwelling shortwave into the ice top for the v1 filler',
        units='W/m^2',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    snowfall = Real(
        'snowfall',
        doc='Uniform frozen-precipitation rate onto the ice top for the v1 filler',
        units='kg/m^2/s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    snow_ice = Bool(
        'snow_ice',
        doc='Enable Archimedes snow-ice flooding conversion (SIS2 SN2IC; default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )

    transport = Bool(
        'transport',
        doc='Enable horizontal category ice/snow transport (default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )

    adv_substeps = Int(
        'adv_substeps',
        doc='Advective sub-iterations per transport call (SIS2 NSTEPS_ADV)',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )

    roll_factor = Real(
        'roll_factor',
        doc='Thin-ice rolling floor factor (SIS2 SEA_ICE_ROLL_FACTOR); 0 disables rolling',
        units='',
        required=False,
        default=1.0,
        has_min=True,
        vmin=0.0,
    )

    dynamics = Bool(
        'dynamics',
        doc='Enable C-grid EVP ice dynamics (default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )

    a_face_stress = Bool(
        'a_face_stress',
        doc='Weight EVP wind stress + ice-ocean drag by face ice concentration (momentum-conserving; default off => byte-identical)',
        units='',
        required=False,
        default=False,
    )

    p0 = Real(
        'p0',
        doc='Ice-strength pressure constant (SIS2 ICE_STRENGTH_PSTAR)',
        units='Pa',
        required=False,
        default=27500.0,
        has_min=True,
        vmin=0.0,
    )

    c0 = Real(
        'c0',
        doc='Ice-strength exponent constant (SIS2 ICE_STRENGTH_CSTAR)',
        units='',
        required=False,
        default=20.0,
        has_min=True,
        vmin=0.0,
    )

    ec = Real(
        'ec',
        doc='Yield-curve axis ratio (SIS2 ICE_YIELD_ELLIPTICITY); 0 = cavitating fluid',
        units='',
        required=False,
        default=2.0,
        has_min=True,
        vmin=0.0,
    )

    cdw = Real(
        'cdw',
        doc='Ice-ocean drag coefficient (SIS2 ICE_CDRAG_WATER)',
        units='',
        required=False,
        default=0.00324,
        has_min=True,
        vmin=0.0,
    )

    rho_ocean = Real(
        'rho_ocean',
        doc='Ice-drag reference density (SIS2 RHO_OCEAN)',
        units='kg/m^3',
        required=False,
        default=1030.0,
        has_min=True,
        vmin=0.0,
    )

    evp_sub_steps = Int(
        'evp_sub_steps',
        doc='EVP subcycles per slow step (SIS2 NSTEPS_DYN)',
        units='',
        required=False,
        default=432,
        has_min=True,
        vmin=1,
    )

    del_sh_min_scale = Real(
        'del_sh_min_scale',
        doc='Viscosity-floor scale (SIS2 ICE_DEL_SH_MIN_SCALE)',
        units='',
        required=False,
        default=2.0,
        has_min=True,
        vmin=0.0,
    )

    tdamp = Real(
        'tdamp',
        doc='Elastic damping timescale rule (SIS2 ICE_TDAMP_ELASTIC): >0 seconds, ==0 auto (0.2*dt_slow), <0 fraction of dt_slow',
        units='',
        required=False,
        default=-0.2,
    )

    cfl_trunc = Real(
        'cfl_trunc',
        doc='Transport-CFL ceiling on the final ice velocity (SIS2 CFL_TRUNCATE, default 0.5 there); 0 disables',
        units='',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    cfl_trunc_dyn_its = Bool(
        'cfl_trunc_dyn_its',
        doc='Also clip the ice velocity every EVP subcycle (SIS2 CFL_TRUNC_DYN_ITS)',
        units='',
        required=False,
        default=False,
    )

    project_ci = Bool(
        'project_ci',
        doc='Project ice concentration forward within the EVP subcycle loop and recompute the ice strength (SIS2 PROJECT_ICE_CONCENTRATION)',
        units='',
        required=False,
        default=False,
    )


class OceanIceIc(Group):
    """`&ocean_ice_ic_nml` -- Sea-ice analytic initial-condition path (PR 24)."""

    _nml_name = 'ocean_ice_ic'

    conc_config = Enum(
        'conc_config',
        doc='Sea-ice initial concentration configuration',
        units='',
        required=False,
        default='zero',
        allowed=('zero', 'uniform', 'latitudes'),
    )

    conc = Real(
        'conc',
        doc="Uniform-mode concentration ('uniform' only)",
        units='',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
        has_max=True,
        vmax=1.0,
    )

    h_ice = Real(
        'h_ice',
        doc='Ice thickness where seeded (SIS2 ICE_INIT_MASS as a thickness)',
        units='m',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    h_snow = Real(
        'h_snow',
        doc='Snow thickness where seeded (SIS2 SNOW_INIT_MASS as a thickness)',
        units='m',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    t_ice = Real(
        't_ice',
        doc='Ice/snow temperature fed through the exact enthalpy inversion (SIS2 ICE_TEMPERATURE_IC)',
        units='degC',
        required=False,
        default=-4.0,
    )

    s_ice = Real(
        's_ice',
        doc='Ice bulk salinity (SIS2 ICE_SALINITY_IC)',
        units='PSU',
        required=False,
        default=4.0,
        has_min=True,
        vmin=0.0,
    )

    arctic_edge = Real(
        'arctic_edge',
        doc="'latitudes' Arctic ice edge (SIS2 ARCTIC_ICE_EDGE_IC)",
        units='degrees_north',
        required=False,
        default=91.0,
        has_min=True,
        vmin=-91.0,
        has_max=True,
        vmax=91.0,
    )

    antarctic_edge = Real(
        'antarctic_edge',
        doc="'latitudes' Antarctic ice edge (SIS2 ANTARCTIC_ICE_EDGE_IC)",
        units='degrees_north',
        required=False,
        default=-91.0,
        has_min=True,
        vmin=-91.0,
        has_max=True,
        vmax=91.0,
    )


class OceanRestore(Group):
    """`&ocean_restore_nml` -- Surface buoyancy restoring (MOM6 RESTOREBUOY): piston-velocity relaxation of top-layer T / S toward targets."""

    _nml_name = 'ocean_restore'

    enable_restore_temp = Bool(
        'enable_restore_temp',
        doc='Master switch for SST restoring',
        units='',
        required=False,
        default=False,
    )

    enable_restore_salt = Bool(
        'enable_restore_salt',
        doc='Master switch for SSS restoring',
        units='',
        required=False,
        default=False,
    )

    piston_t = Real(
        'piston_t',
        doc='SST piston velocity (MOM6 FLUXCONST_T)',
        units='m/day',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    piston_s = Real(
        'piston_s',
        doc='SSS piston velocity (MOM6 FLUXCONST_S)',
        units='m/day',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    restore_sst = Real(
        'restore_sst',
        doc='Scalar target SST',
        units='degC',
        required=False,
        default=0.0,
    )

    restore_sss = Real(
        'restore_sss',
        doc='Scalar target SSS',
        units='PSU',
        required=False,
        default=0.0,
    )


class OceanGeothermal(Group):
    """`&ocean_geothermal_nml` -- Geothermal bottom heat flux (bed-side analogue of the surface heat flux)."""

    _nml_name = 'ocean_geothermal'

    enable = Bool(
        'enable',
        doc='Apply a constant geothermal bottom heat flux to the bed layer',
        units='',
        required=False,
        default=False,
    )

    q_geo = Real(
        'q_geo',
        doc='Constant bottom heat flux (positive into the ocean from below)',
        units='W/m^2',
        required=False,
        default=0.0,
    )


class OceanSponge(Group):
    """`&ocean_sponge_nml` -- Map-driven sponge: per-cell Idamp + 3-D reference state (PR-23)."""

    _nml_name = 'ocean_sponge'

    enable = Bool(
        'enable',
        doc='Master switch (default off; legacy band kernels run when off)',
        units='',
        required=False,
        default=False,
    )

    damp_source = Enum(
        'damp_source',
        doc='How idamp_h/u/v are filled',
        units='',
        required=False,
        default='band',
        allowed=('band', 'file'),
    )

    target_source = Enum(
        'target_source',
        doc='Reference-state source',
        units='',
        required=False,
        default='ic',
        allowed=('ic', 'linear_z', 'file'),
    )

    ramp = Enum(
        'ramp',
        doc="Band ramp shape from the sponge wall inward ('linear' is ISOMIP+ Eq. 20)",
        units='',
        required=False,
        default='cosine',
        allowed=('cosine', 'linear'),
    )

    lin_t_ref = Real(
        'lin_t_ref',
        doc="target_source='linear_z': T at the z = 0 datum",
        units='degC',
        required=False,
        default=0.0,
    )

    lin_dt_dz = Real(
        'lin_dt_dz',
        doc="target_source='linear_z': dT/dz, z positive UP (stable => > 0)",
        units='degC/m',
        required=False,
        default=0.0,
    )

    lin_s_ref = Real(
        'lin_s_ref',
        doc="target_source='linear_z': S at the z = 0 datum",
        units='PSU',
        required=False,
        default=35.0,
    )

    lin_ds_dz = Real(
        'lin_ds_dz',
        doc="target_source='linear_z': dS/dz, z positive UP (stable => < 0)",
        units='PSU/m',
        required=False,
        default=0.0,
    )

    relax_uv = Bool(
        'relax_uv',
        doc='Relax u/v toward u_ref/v_ref',
        units='',
        required=False,
        default=True,
    )

    relax_tracers = Bool(
        'relax_tracers',
        doc='Relax every registered tracer toward its 3-D reference field',
        units='',
        required=False,
        default=True,
    )

    relax_h = Bool(
        'relax_h',
        doc="Interior-interface thickness damping (NOT IMPLEMENTED in v1 — deferred to PR-23b, requires vcoord_type='lagrangian')",
        units='',
        required=False,
        default=False,
    )

    west_width = Int(
        'west_width',
        doc='Per-edge band-width override, cells (<0 => inherit &ocean_bc_nml sponge_width)',
        units='',
        required=False,
        default=-1,
    )

    east_width = Int('east_width', doc='as west_width', units='', required=False, default=-1)

    south_width = Int('south_width', doc='as west_width', units='', required=False, default=-1)

    north_width = Int('north_width', doc='as west_width', units='', required=False, default=-1)

    west_strength = Real(
        'west_strength',
        doc='Per-edge peak relaxation-rate override (<0 => inherit &ocean_bc_nml sponge_strength)',
        units='1/s',
        required=False,
        default=-1.0,
    )

    east_strength = Real(
        'east_strength',
        doc='as west_strength',
        units='1/s',
        required=False,
        default=-1.0,
    )

    south_strength = Real(
        'south_strength',
        doc='as west_strength',
        units='1/s',
        required=False,
        default=-1.0,
    )

    north_strength = Real(
        'north_strength',
        doc='as west_strength',
        units='1/s',
        required=False,
        default=-1.0,
    )


class OceanTracers(Group):
    """`&ocean_tracers_nml` -- Prognostic-tracer registry switches."""

    _nml_name = 'ocean_tracers'

    enable_ideal_age = Bool(
        'enable_ideal_age',
        doc='Register a passive ideal-age tracer',
        units='',
        required=False,
        default=False,
    )

    ideal_age_young_val = Real(
        'ideal_age_young_val',
        doc="Surface-band ideal-age Dirichlet value (0 = today's hard-coded reset)",
        units='s',
        required=False,
        default=0.0,
    )

    ideal_age_sfc_growth_rate = Real(
        'ideal_age_sfc_growth_rate',
        doc='Exponential growth rate of the ideal-age surface value (0 = constant)',
        units='1/s',
        required=False,
        default=0.0,
    )

    enable_pseudo_salt = Bool(
        'enable_pseudo_salt',
        doc="Register the pseudo-salt verification tracer (diagnostic; seeded to S, given salinity's boundary fluxes)",
        units='',
        required=False,
        default=False,
    )


class OceanBt(Group):
    """`&ocean_bt_nml` -- Split-explicit barotropic substep controls."""

    _nml_name = 'ocean_bt'

    n_inner = Int(
        'n_inner',
        doc='Barotropic substeps per outer step (0 = unsplit)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    auto_n_inner = Bool(
        'auto_n_inner',
        doc='Derive n_inner from the gravity-wave CFL at setup',
        units='',
        required=False,
        default=False,
    )

    cfl_bt_safety = Real(
        'cfl_bt_safety',
        doc='Safety fraction on the BT CFL when auto_n_inner',
        units='',
        required=False,
        default=0.65,
    )

    bebt = Real(
        'bebt',
        doc='Forward-velocity-projection weight (MOM6 BEBT)',
        units='',
        required=False,
        default=0.0,
    )

    use_cont_type = Bool(
        'use_cont_type',
        doc='Use the BT_cont flux-bounded closure (MOM6 USE_BT_CONT_TYPE)',
        units='',
        required=False,
        default=False,
    )

    cont_corr_bounds = Bool(
        'cont_corr_bounds',
        doc='Use BT_cont flux limits for the eta-correction bound',
        units='',
        required=False,
        default=False,
    )

    upstream_h_face = Bool(
        'upstream_h_face',
        doc='Use per-face upstream-PPM column-sum thickness in the BT chain',
        units='',
        required=False,
        default=False,
    )

    correction_h_weighted = Bool(
        'correction_h_weighted',
        doc='h-weight the post-substep BT corrector (MOM6 frhatu)',
        units='',
        required=False,
        default=False,
    )

    correction_visc_rem = Bool(
        'correction_visc_rem',
        doc='h*visc_rem joint corrector weight (requires correction_h_weighted; visc_rem is produced by vdiff and is inert, =1, without ocean_vdiff_nml implicit_drag)',
        units='',
        required=False,
        default=False,
    )

    split_scheme = Enum(
        'split_scheme',
        doc='Outer split-explicit time scheme: pred_corr (DEFAULT; MOM6 predictor-corrector, slow tendencies on the u_av/h_av step time-means, forward-backward gravity-wave pairing; lifts the internal-wave dt ceiling) or ssp_rk2 (EXPERIMENTAL; two-stage SSP average, widest envelope — the only scheme wired through eulerian_z, wet/dry and dt_tracer_advect_ratio>1 — but it spuriously grows internal gravity waves out of a stratified REST state, En 2.992E-05 vs 1.739E-09 at day 25 on resting_stratified_channel.nml; a (omega*dt)^4 noise floor, so forced viscous runs sit decades above it and quiescent or long spin-up runs do not)',
        units='',
        required=False,
        default='pred_corr',
        allowed=('pred_corr', 'ssp_rk2'),
    )

    pc_be = Real(
        'pc_be',
        doc='pred_corr predictor fraction BE (MOM6 BE, 0.6 reference)',
        units='',
        required=False,
        default=0.6,
        has_min=True,
        vmin=0.0,
        has_max=True,
        vmax=1.0,
    )

    renorm_visc_rem = Bool(
        'renorm_visc_rem',
        doc='gamma-weighted continuity transport-matching inversion (MOM6 u_cor = u + du*visc_rem; requires correction_visc_rem)',
        units='',
        required=False,
        default=False,
    )

    forcing_visc_rem = Bool(
        'forcing_visc_rem',
        doc='MOM6 wt_u parity: h*visc_rem-weight the BT forcing depth-mean so friction-damped (glued) layers do not force the fast loop (requires correction_visc_rem)',
        units='',
        required=False,
        default=False,
    )

    correction_bc_pgf = Bool(
        'correction_bc_pgf',
        doc='Per-layer baroclinic-PGF retro-correction for the eta change',
        units='',
        required=False,
        default=False,
    )

    substep_drag = Bool(
        'substep_drag',
        doc='Apply a per-substep BT velocity damping factor',
        units='',
        required=False,
        default=False,
    )

    substep_zeta_ke = Bool(
        'substep_zeta_ke',
        doc='Integrate live zeta_bt + KE-gradient in the BT fast loop (.false. = MOM6 parity: planetary Coriolis only, zeta/KE frozen in the slow forcing)',
        units='',
        required=False,
        default=True,
    )

    wave_drag = Bool(
        'wave_drag',
        doc='Barotropic linear wave drag master switch (MOM6 BT_LINEAR_WAVE_DRAG)',
        units='',
        required=False,
        default=False,
    )

    wave_drag_form = Enum(
        'wave_drag_form',
        doc='Wave-drag r_H filler',
        units='',
        required=False,
        default='uniform',
        allowed=('uniform', 'roughness_proxy', 'file'),
    )

    wave_drag_scale = Real(
        'wave_drag_scale',
        doc='Global tuning multiplier on r_H (MOM6 BT_WAVE_DRAG_SCALE)',
        units='',
        required=False,
        default=1.0,
    )

    wave_drag_r_uniform = Real(
        'wave_drag_r_uniform',
        doc="Piston velocity r_H for wave_drag_form='uniform'",
        units='m/s',
        required=False,
        default=0.0,
    )

    wave_drag_kappa = Real(
        'wave_drag_kappa',
        doc="Topographic wavenumber for wave_drag_form='roughness_proxy'",
        units='1/m',
        required=False,
        default=0.00062832,
    )

    wave_drag_n_bot = Real(
        'wave_drag_n_bot',
        doc="Reference bottom N for wave_drag_form='roughness_proxy'",
        units='1/s',
        required=False,
        default=0.001,
    )

    wave_drag_h2_max = Real(
        'wave_drag_h2_max',
        doc="Ceiling on <h^2> proxy for wave_drag_form='roughness_proxy'",
        units='m^2',
        required=False,
        default=25000.0,
    )

    wave_drag_file = Str(
        'wave_drag_file',
        doc='Reserved for PR-14 (MOM6 BT_WAVE_DRAG_FILE); unused today',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    wave_drag_var = Str(
        'wave_drag_var',
        doc='Reserved for PR-14 (MOM6 BT_WAVE_DRAG_VAR); unused today',
        units='',
        required=False,
        default='rH',
        max_len=64,
    )

    bt_halo = Int(
        'bt_halo',
        doc='Wide-halo BT march-in width (-1 = auto: 8 under multi-rank if compatible, else 0; 0 = explicit off, bit-identical)',
        units='',
        required=False,
        default=-1,
        has_min=True,
        vmin=-1,
    )


class OceanDebug(Group):
    """`&ocean_debug_nml` -- Forensic probes for the ocean dyn-core (all heavy, all default off)."""

    _nml_name = 'ocean_debug'

    budget = Bool(
        'budget',
        doc='Per-stage BT power-budget probe (BUDGET rows)',
        units='',
        required=False,
        default=False,
    )

    ke_attr = Bool(
        'ke_attr',
        doc='Per-segment layer-KE attribution meter (KE_ATTR rows)',
        units='',
        required=False,
        default=False,
    )

    ke_attr_start_step = Int(
        'ke_attr_start_step',
        doc='First outer step the KE meter samples (0 = from start)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    ke_attr_end_step = Int(
        'ke_attr_end_step',
        doc='Last outer step the KE meter samples (0 = unbounded)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    chksum = Bool(
        'chksum',
        doc='MOM6-style per-phase field checksums (CHKSUM rows) + HOTFACE argmax-face anatomy at the tendency seams',
        units='',
        required=False,
        default=False,
    )

    chksum_start_step = Int(
        'chksum_start_step',
        doc='First outer step the chksum probe samples (0 = from start)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    chksum_end_step = Int(
        'chksum_end_step',
        doc='Last outer step the chksum probe samples (0 = unbounded)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    chksum_interior = Bool(
        'chksum_interior',
        doc='Reduce over physical cells only, making `bits` a valid 1-rank-vs-N-rank gate',
        units='',
        required=False,
        default=False,
    )


class OceanMpi(Group):
    """`&ocean_mpi_nml` -- Multi-rank MPI debug/tuning controls."""

    _nml_name = 'ocean_mpi'

    poison_ghosts = Bool(
        'poison_ghosts',
        doc='Sentinel-NaN the exchange-covered ghost bands at each outer-step start; unexchanged-ghost consumption becomes a loud NaN. Default off = bit-identical.',
        units='',
        required=False,
        default=False,
    )


class OceanWetdry(Group):
    """`&ocean_wetdry_nml` -- Dynamic wetting/drying for the split-explicit barotropic substep."""

    _nml_name = 'ocean_wetdry'

    enable = Bool(
        'enable',
        doc='Dynamic wet/dry: upwind BT face thickness + positive-definite outflow limiter + bed-blocking momentum gate',
        units='',
        required=False,
        default=False,
    )

    dry_depth = Real(
        'dry_depth',
        doc='Total-depth dry threshold',
        units='m',
        required=False,
        default=0.05,
        has_min=True,
        vmin=0.0,
    )

    rewet_depth = Real(
        'rewet_depth',
        doc='Hysteresis re-wet threshold (must exceed dry_depth)',
        units='m',
        required=False,
        default=0.1,
        has_min=True,
        vmin=0.0,
    )

    land_margin = Real(
        'land_margin',
        doc='Static-land headroom above rest MSL (intertidal cutoff)',
        units='m',
        required=False,
        default=5.0,
        has_min=True,
        vmin=0.0,
    )


class OceanPgf(Group):
    """`&ocean_pgf_nml` -- Pressure-gradient-force kernel selector + knobs."""

    _nml_name = 'ocean_pgf'

    form = Enum(
        'form',
        doc='PGF kernel variant',
        units='',
        required=False,
        default='mont',
        allowed=('mont', 'montgomery', 'fv_lite', 'fv_wright', 'gprime', 'fv_mom6'),
    )

    gprime_gfs = Real(
        'gprime_gfs',
        doc='Free-surface gravity for the gprime PGF',
        units='m/s^2',
        required=False,
        default=9.81,
    )

    gprime_gint = Real(
        'gprime_gint',
        doc='Internal-interface reduced gravity (gprime PGF)',
        units='m/s^2',
        required=False,
        default=0.0098,
    )

    gfs_scale = Real(
        'gfs_scale',
        doc='Free-surface gravity scaling (FV_MOM6, MOM6 GFS_scale)',
        units='',
        required=False,
        default=1.0,
    )

    maxvel = Real(
        'maxvel',
        doc='Velocity-truncation clamp (0 = disabled)',
        units='m/s',
        required=False,
        default=0.0,
    )

    cfl_trunc = Real(
        'cfl_trunc',
        doc='Advective-CFL velocity truncation threshold (0 = disabled)',
        units='',
        required=False,
        default=0.0,
    )

    mass_weight = Bool(
        'mass_weight',
        doc='FV_MOM6 shelf-break hWght mass-weighting at unequal-depth faces',
        units='',
        required=False,
        default=False,
    )

    reconstruct_for_pressure = Bool(
        'reconstruct_for_pressure',
        doc='FV_MOM6 in-layer PLM/PPM T/S reconstruction for the density integral',
        units='',
        required=False,
        default=False,
    )

    recon_scheme = Int(
        'recon_scheme',
        doc='In-layer reconstruction scheme: 1=PLM, 2=PPM',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
        has_max=True,
        vmax=2,
    )

    p_top_in_bc = Bool(
        'p_top_in_bc',
        doc='FV_MOM6: add the top-of-column load ms%p_top to the pressure-stack surface boundary condition pa(nz+1) = rho_ref*g*eta + p_top',
        units='',
        required=False,
        default=False,
    )


class OceanEos(Group):
    """`&ocean_eos_nml` -- Equation-of-state variant selector, liquidus set + reference pressure."""

    _nml_name = 'ocean_eos'

    eos = Enum(
        'eos',
        doc='Equation-of-state variant',
        units='',
        required=False,
        default='linear',
        allowed=('linear', 'wright', 'roquet_spv', 'teos10'),
    )

    tfreeze_set = Enum(
        'tfreeze_set',
        doc="Named liquidus coefficient set for eos_freezing_point (T_f = l1*S + l2 + l3*p): 'seaice' = SIS2/MOM6 (-0.054, 0, -7.53e-8), 'isomip' = ISOMIP+ (-0.0573, 0.0832, -7.53e-8)",
        units='',
        required=False,
        default='seaice',
        allowed=('seaice', 'isomip'),
    )

    p_ref = Real(
        'p_ref',
        doc='Reference pressure for the potential density ms%rho_layer (horizontally uniform by design; 0 => surface density)',
        units='Pa',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )


class OceanBdrag(Group):
    """`&ocean_bdrag_nml` -- Bottom-drag selector + coefficients."""

    _nml_name = 'ocean_bdrag'

    form = Enum(
        'form',
        doc='Bottom-drag variant',
        units='',
        required=False,
        default='quadratic',
        allowed=('quadratic', 'linear'),
    )

    cd = Real(
        'cd',
        doc='Quadratic drag coefficient (0 disables)',
        units='',
        required=False,
        default=0.0,
    )

    r = Real(
        'r',
        doc='Linear Rayleigh coefficient (0 disables)',
        units='1/s',
        required=False,
        default=0.0,
    )

    hbbl = Real(
        'hbbl',
        doc='BBL thickness for distributed drag (0 = bed-only)',
        units='m',
        required=False,
        default=0.0,
    )

    bg_vel = Real(
        'bg_vel',
        doc='Background velocity floor for distributed drag',
        units='m/s',
        required=False,
        default=0.0,
    )

    bbl_thick_min = Real(
        'bbl_thick_min',
        doc='Minimum effective BBL thickness',
        units='m',
        required=False,
        default=0.0,
    )

    bed_factor = Real(
        'bed_factor',
        doc='Multiplier on the bed-layer drag tendency only',
        units='',
        required=False,
        default=1.0,
    )

    channel_drag = Bool(
        'channel_drag',
        doc='Enable per-layer lateral side-wall (channel) Rayleigh drag',
        units='',
        required=False,
        default=False,
    )

    cdrag_side = Real(
        'cdrag_side',
        doc='Side-wall drag coefficient (0 disables channel drag)',
        units='',
        required=False,
        default=0.0,
    )

    implicit = Bool(
        'implicit',
        doc='Backward-Euler implicit drag (stable for thin bottom layers)',
        units='',
        required=False,
        default=False,
    )


class OceanTdrag(Group):
    """`&ocean_tdrag_nml` -- Ice-shelf top-drag selector + coefficients (mirror of &ocean_bdrag_nml)."""

    _nml_name = 'ocean_tdrag'

    enable = Bool(
        'enable',
        doc='Enable the ice-shelf top drag (requires &ocean_cavity_dyn_nml enable)',
        units='',
        required=False,
        default=False,
    )

    form = Enum(
        'form',
        doc='Top-drag variant',
        units='',
        required=False,
        default='quadratic',
        allowed=('quadratic', 'linear'),
    )

    cd = Real(
        'cd',
        doc='Quadratic top-drag coefficient (0 disables); must equal &ocean_cavity_melt_nml cdrag_top when melt is on',
        units='',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    r = Real(
        'r',
        doc='Linear Rayleigh top-drag coefficient (0 disables)',
        units='1/s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    htbl = Real(
        'htbl',
        doc='Top-boundary-layer thickness for distributed drag (0 = layer-nz only)',
        units='m',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    bg_vel = Real(
        'bg_vel',
        doc='Background velocity floor in the quadratic top-drag speed',
        units='m/s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    tbl_thick_min = Real(
        'tbl_thick_min',
        doc='Minimum effective TBL thickness (0 = fall back to h_min)',
        units='m',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    implicit = Bool(
        'implicit',
        doc='Backward-Euler top drag inside the drag kernel (stable for thin top layers)',
        units='',
        required=False,
        default=False,
    )


class OceanHdiff(Group):
    """`&ocean_hdiff_nml` -- Along-coordinate (not neutral) constant-coefficient tracer diffusion."""

    _nml_name = 'ocean_hdiff'

    kappa_h = Real(
        'kappa_h',
        doc="Horizontal tracer diffusivity, along the model coordinate (0 disables; ocean path only — coastal's equivalent knob is top-level hdiff_kappa)",
        units='m^2/s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )


class OceanHvisc(Group):
    """`&ocean_hvisc_nml` -- Lateral-viscosity closure + coefficients."""

    _nml_name = 'ocean_hvisc'

    lateral_closure = Enum(
        'lateral_closure',
        doc='Lateral-mixing closure tag',
        units='',
        required=False,
        default='none',
        allowed=('none', 'leith', 'smagorinsky', 'smag', 'biharmonic', 'leith_biharm'),
    )

    c_smag = Real(
        'c_smag',
        doc='Smagorinsky coefficient',
        units='',
        required=False,
        default=0.15,
    )

    c_leith = Real('c_leith', doc='Leith coefficient', units='', required=False, default=1.0)

    kh_vel_scale = Real(
        'kh_vel_scale',
        doc='Velocity scale for the resolution viscosity floor',
        units='m/s',
        required=False,
        default=0.0,
    )

    ah_bg = Real(
        'ah_bg',
        doc='Background harmonic viscosity (<0 = derive)',
        units='m^2/s',
        required=False,
        default=-1.0,
    )

    ah_max = Real(
        'ah_max',
        doc='Cap on per-face harmonic viscosity',
        units='m^2/s',
        required=False,
        default=10000.0,
    )

    smag_ah = Bool(
        'smag_ah',
        doc='Flow-aware biharmonic viscosity (MOM6 SMAGORINSKY_AH)',
        units='',
        required=False,
        default=False,
    )

    smag_bi_const = Real(
        'smag_bi_const',
        doc='Biharmonic Smagorinsky constant (MOM6 SMAG_BI_CONST)',
        units='',
        required=False,
        default=0.06,
    )

    c_leith_bi = Real(
        'c_leith_bi',
        doc='Biharmonic Leith constant (MOM6 LEITH_BI_CONST)',
        units='',
        required=False,
        default=0.0,
    )

    nu_4_bg = Real(
        'nu_4_bg',
        doc='Background biharmonic viscosity floor',
        units='m^4/s',
        required=False,
        default=0.0,
    )

    nu_4_max = Real(
        'nu_4_max',
        doc='Cap on per-face biharmonic viscosity',
        units='m^4/s',
        required=False,
        default=1000000000000.0,
    )

    nu_h = Real(
        'nu_h',
        doc='Constant horizontal eddy viscosity',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    nu_4 = Real(
        'nu_4',
        doc='Constant biharmonic eddy viscosity',
        units='m^4/s',
        required=False,
        default=0.0,
    )

    no_slip = Bool(
        'no_slip',
        doc='Coastal lateral BC: .false.=free-slip (×wet_q), .true.=no-slip (×(2-wet_q))',
        units='',
        required=False,
        default=False,
    )

    stress_tensor = Bool(
        'stress_tensor',
        doc='MOM6 thickness-weighted stress-div operator + per-cell CFL + coast-mask',
        units='',
        required=False,
        default=False,
    )

    bound_coef = Real(
        'bound_coef',
        doc='Per-cell viscosity-CFL safety coefficient (MOM6 HORVISC_BOUND_COEF)',
        units='',
        required=False,
        default=0.8,
    )

    bound_kh = Bool(
        'bound_kh',
        doc='Per-face harmonic viscosity CFL clamp on the velocity-Laplacian paths (MOM6 BOUND_KH)',
        units='',
        required=False,
        default=False,
    )

    resoln_scaled_visc = Bool(
        'resoln_scaled_visc',
        doc='Scale dynamic LAPLACIAN viscosity (Leith/Smag_KH A_h, not biharmonic) by VarMix Res_fn (MOM6 RESOLN_SCALED_KH)',
        units='',
        required=False,
        default=False,
    )

    kh_vel_scale_live = Real(
        'kh_vel_scale_live',
        doc='Live velocity-scale viscosity Kh=U*dx*|u| (0=off; distinct from kh_vel_scale floor)',
        units='m/s',
        required=False,
        default=0.0,
    )

    kh_aniso = Real(
        'kh_aniso',
        doc='Anisotropic Laplacian viscosity magnitude (MOM6 KH_ANISO; stress_tensor path only)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    aniso_mode = Int(
        'aniso_mode',
        doc='Anisotropy direction mode (0=grid-relative aniso_dir; MOM6 ANISOTROPIC_MODE)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    aniso_dir = RealArray(
        'aniso_dir',
        doc='Anisotropy direction (n1,n2) in grid i,j components (MOM6 ANISO_GRID_DIR)',
        units='',
        required=False,
        default=(1.0, 0.0),
        size=2,
    )


class OceanVmix(Group):
    """`&ocean_vmix_nml` -- Vertical-mixing module switches + knobs."""

    _nml_name = 'ocean_vmix'

    use_closure = Bool(
        'use_closure',
        doc='Master switch for the interior closure (PP81)',
        units='',
        required=False,
        default=True,
    )

    use_kpp = Bool(
        'use_kpp',
        doc='KPP surface-boundary-layer overlay (needs use_closure)',
        units='',
        required=False,
        default=True,
    )

    direct_stress = Bool(
        'direct_stress',
        doc='Spread wind stress over hmix_stress (MOM6 DIRECT_STRESS)',
        units='',
        required=False,
        default=False,
    )

    hmix_stress = Real(
        'hmix_stress',
        doc='Surface-slab thickness for direct_stress',
        units='m',
        required=False,
        default=20.0,
    )

    kv_ml_invz2 = Real(
        'kv_ml_invz2',
        doc='Near-surface 1/(z hmix)^2 viscosity (MOM6 KV_ML_INVZ2)',
        units='m^2/s',
        required=False,
        default=0.0,
    )

    hmix_fixed = Real(
        'hmix_fixed',
        doc='Mixed-layer thickness for the KV_ML_INVZ2 profile',
        units='m',
        required=False,
        default=20.0,
    )

    harmonic_visc = Bool(
        'harmonic_visc',
        doc='Harmonic-mean face thickness in vdiff (MOM6 HARMONIC_VISC)',
        units='',
        required=False,
        default=False,
    )

    dt_therm_ratio = Int(
        'dt_therm_ratio',
        doc='Thermo step runs at ratio*dt_dyn (MOM6 DT_THERM)',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )

    dt_tracer_advect_ratio = Int(
        'dt_tracer_advect_ratio',
        doc='Horizontal tracer advect runs at ratio*dt_dyn (MOM6 DT_TRACER_ADVECT)',
        units='',
        required=False,
        default=1,
        has_min=True,
        vmin=1,
    )

    tracer_recon = Enum(
        'tracer_recon',
        doc='Windowed tracer-advect drain face reconstruction (ocean): ppm|weno5|weno7|weno9',
        units='',
        required=False,
        default='ppm',
        allowed=('ppm', 'weno5', 'weno7', 'weno9'),
    )

    kv_max = Real(
        'kv_max',
        doc='Assembly ceiling on kv; huge=off (MOM6 Kd_max momentum)',
        units='m^2/s',
        required=False,
        default=1.7976931348623157e+308,
    )

    kd_max = Real(
        'kd_max',
        doc='Assembly ceiling on kt/ks; huge=off (MOM6 Kd_max)',
        units='m^2/s',
        required=False,
        default=1.7976931348623157e+308,
    )

    kd_smooth_iterations = Int(
        'kd_smooth_iterations',
        doc='1-2-1 smoothing passes on kv/kt (MOM6 Kd_smooth)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    vmix_guard = Bool(
        'vmix_guard',
        doc='Debug-gated negative/NaN diffusivity guard (assembly)',
        units='',
        required=False,
        default=False,
    )

    bkgnd_profile = Bool(
        'bkgnd_profile',
        doc='Bryan-Lewis depth-varying background diffusivity',
        units='',
        required=False,
        default=False,
    )

    bkgnd_kd_sfc = Real(
        'bkgnd_kd_sfc',
        doc='Bryan-Lewis surface-asymptote background Kd',
        units='m^2/s',
        required=False,
        default=1e-05,
    )

    bkgnd_kd_deep = Real(
        'bkgnd_kd_deep',
        doc='Bryan-Lewis deep-asymptote background Kd',
        units='m^2/s',
        required=False,
        default=0.00013,
    )

    bkgnd_z0 = Real(
        'bkgnd_z0',
        doc='Bryan-Lewis transition-centre depth',
        units='m',
        required=False,
        default=2500.0,
    )

    bkgnd_delta = Real(
        'bkgnd_delta',
        doc='Bryan-Lewis transition half-width',
        units='m',
        required=False,
        default=222.0,
        has_min=True,
        vmin=1e-06,
    )

    bkgnd_prandtl = Real(
        'bkgnd_prandtl',
        doc='Background Prandtl number Kv_bg=prandtl*Kd_bg',
        units='',
        required=False,
        default=1.0,
        has_min=True,
        vmin=0.0,
    )

    bkgnd_henyey = Bool(
        'bkgnd_henyey',
        doc='Henyey IGW latitude factor on the scalar background (excludes bkgnd_profile; needs a non-cartesian grid)',
        units='',
        required=False,
        default=False,
    )

    bkgnd_kd_min = Real(
        'bkgnd_kd_min',
        doc='Minimum background Kd under the Henyey scaling (negative = 0.01*kt_bg)',
        units='m^2/s',
        required=False,
        default=-1.0,
    )

    bkgnd_henyey_n0_2omega = Real(
        'bkgnd_henyey_n0_2omega',
        doc='Henyey N0/(2*Omega) reference stratification ratio',
        units='',
        required=False,
        default=20.0,
        has_min=True,
        vmin=1.0,
    )

    bkgnd_henyey_max_lat = Real(
        'bkgnd_henyey_max_lat',
        doc='Latitude poleward of which the Henyey factor floors',
        units='degN',
        required=False,
        default=95.0,
        has_min=True,
        vmin=0.0,
    )

    pp81_nu0 = Real(
        'pp81_nu0',
        doc='PP81 Richardson-dependent viscosity scale',
        units='m^2/s',
        required=False,
        default=0.01,
        has_min=True,
        vmin=0.0,
    )

    pp81_nu_bg = Real(
        'pp81_nu_bg',
        doc='PP81 background viscosity (also seeds vmix%kv_bg)',
        units='m^2/s',
        required=False,
        default=0.0001,
        has_min=True,
        vmin=0.0,
    )

    pp81_kappa_bg = Real(
        'pp81_kappa_bg',
        doc='PP81 background diffusivity (also seeds vmix%kt_bg/ks_bg)',
        units='m^2/s',
        required=False,
        default=1e-05,
        has_min=True,
        vmin=0.0,
    )

    pp81_alpha = Real(
        'pp81_alpha',
        doc='PP81 Richardson-number scaling coefficient (paper value 5)',
        units='',
        required=False,
        default=5.0,
        has_min=True,
        vmin=1e-12,
    )

    shear2_floor = Real(
        'shear2_floor',
        doc='Floor on |du/dz|^2+|dv/dz|^2 in the PP81 Ri denominator',
        units='1/s^2',
        required=False,
        default=1e-10,
        has_min=True,
        vmin=1e-30,
    )

    kpp_ri_crit = Real(
        'kpp_ri_crit',
        doc='Critical bulk Richardson number for KPP BL-depth',
        units='',
        required=False,
        default=0.3,
        has_min=True,
        vmin=1e-12,
    )

    kpp_cs_nonlocal = Real(
        'kpp_cs_nonlocal',
        doc='KPP non-local (counter-gradient) transport coefficient C_s',
        units='',
        required=False,
        default=6.3,
        has_min=True,
        vmin=0.0,
    )

    kpp_c_vt2 = Real(
        'kpp_c_vt2',
        doc='KPP unresolved-turbulence V_t^2 coefficient (0 disables V_t^2)',
        units='',
        required=False,
        default=1.8,
        has_min=True,
        vmin=0.0,
    )

    buoyancy_coeffs = Enum(
        'buoyancy_coeffs',
        doc='Source of alpha/beta for KPP B_0 + double diffusion: constant (scalar &ocean_ic_nml pair) | eos (active EOS derivatives)',
        units='',
        required=False,
        default='constant',
        allowed=('constant', 'eos'),
    )


class OceanVdiff(Group):
    """`&ocean_vdiff_nml` -- Backward-Euler vertical-friction solver knobs."""

    _nml_name = 'ocean_vdiff'

    implicit_stress = Bool(
        'implicit_stress',
        doc='Fold wind stress into the vdiff surface (k=nz) RHS',
        units='',
        required=False,
        default=False,
    )

    accel_visc_rem = Bool(
        'accel_visc_rem',
        doc='MOM6 parity: attenuate the slow explicit accelerations by the per-layer viscous remnant (u = u0 + visc_rem*(u-u0) after the applies); requires ocean_bt_nml correction_visc_rem',
        units='',
        required=False,
        default=False,
    )

    implicit_drag = Bool(
        'implicit_drag',
        doc='Fold bottom drag into the vdiff bed (k=1) diagonal',
        units='',
        required=False,
        default=False,
    )

    implicit_top_drag = Bool(
        'implicit_top_drag',
        doc='Fold the ice-shelf top drag into the vdiff surface (k=nz) diagonal, masking the wind RHS under cover',
        units='',
        required=False,
        default=False,
    )

    hvel_mom6 = Bool(
        'hvel_mom6',
        doc='MOM6 HARMONIC_VISC parity: harmonic momentum face thickness with the near-bed upwind blend, and arithmetic h_shear. Suppresses grounded-sliver momentum as MOM6 does',
        units='',
        required=False,
        default=False,
    )

    hbbl_visc = Real(
        'hbbl_visc',
        doc='Bottom-layer scale for the hvel_mom6 botfn blend (MOM6 HBBL)',
        units='',
        required=False,
        default=10.0,
    )

    bbl_glue = Bool(
        'bbl_glue',
        doc='MOM6 bottomdraglaw coupling parity: kv_bbl botfn glue at near-bed interfaces + piston bed drag. Absorbs the spurious grounded-layer PGF as MOM6 does (PGF_BUG.md par.9). Requires hvel_mom6 + implicit_drag + linear bottom drag',
        units='',
        required=False,
        default=False,
    )

    bbl_piston = Real(
        'bbl_piston',
        doc='BBL drag piston velocity u* for bbl_glue (MOM6 CDRAG*DRAG_BG_VEL)',
        units='m/s',
        required=False,
        default=0.0003,
        has_min=True,
        vmin=0.0,
    )

    hvel_upwind = Bool(
        'hvel_upwind',
        doc='Near-bed upwind blend in the hvel_mom6 thickness build (.false. = pure harmonic; the u-sign blend flip-flops on roundoff at rest and collapses the BBL glue, PGF_BUG.md par.9.8)',
        units='',
        required=False,
        default=True,
    )


class OceanContinuity(Group):
    """`&ocean_continuity_nml` -- Continuity-PPM positivity controls."""

    _nml_name = 'ocean_continuity'

    h_min = Real(
        'h_min',
        doc='PPM positivity-limiter thickness floor',
        units='m',
        required=False,
        default=1e-06,
    )

    ppm_limit_pos = Bool(
        'ppm_limit_pos',
        doc='Positivity-preserving PPM face limiter (MOM6 PPM_limit_pos)',
        units='',
        required=False,
        default=False,
    )

    vol_cfl = Bool(
        'vol_cfl',
        doc='Swept-volume continuity-PPM face thickness (MOM6 vol_CFL)',
        units='',
        required=False,
        default=False,
    )

    positive_definite = Bool(
        'positive_definite',
        doc='Positive-definite split continuity (h>=h_lim, zero mass created; reconstruction floor + per-donor outflux limiter)',
        units='',
        required=False,
        default=False,
    )


class OceanIsopycnal(Group):
    """`&ocean_isopycnal_nml` -- Lagrangian grounding-stability controls."""

    _nml_name = 'ocean_isopycnal'

    angstrom_h = Real(
        'angstrom_h',
        doc='Minimum-thickness floor on the Lagrangian continuity h-update (MOM6 Angstrom_H analogue)',
        units='m',
        required=False,
        default=0.0,
    )

    reset_vanished_u = Bool(
        'reset_vanished_u',
        doc='Zero face velocity when layer vanished on both adjacent cells',
        units='',
        required=False,
        default=False,
    )

    cfl_ignore_vanished = Bool(
        'cfl_ignore_vanished',
        doc='Exclude vanished layers from MaxCFL / panic / CFL truncation',
        units='',
        required=False,
        default=False,
    )

    pgf_skip_nonoverlap = Bool(
        'pgf_skip_nonoverlap',
        doc="Zero the face PGF where a grounded layer's z-extents do not overlap across the face (VCOORD_LAGRANGIAN only, mont / fv_lite / fv_wright / fv_mom6; default ON — kills the spurious at-rest grounded-layer pressure gradient)",
        units='',
        required=False,
        default=True,
    )

    conservative_floor = Bool(
        'conservative_floor',
        doc='Conservative min-thickness borrow (replaces the injecting angstrom_h floor)',
        units='',
        required=False,
        default=False,
    )

    check_h_positive = Bool(
        'check_h_positive',
        doc='DEBUG: abort on the first negative h_layer, naming the stage + (i,j,k)',
        units='',
        required=False,
        default=False,
    )


class OceanTopo(Group):
    """`&ocean_topo_nml` -- Basin geometry + surface forcing + Coriolis tilt."""

    _nml_name = 'ocean_topo'

    topo_config = Enum(
        'topo_config',
        doc='Bathymetry profile selector',
        units='',
        required=False,
        default='flat',
        allowed=('flat', 'spoon', 'seamount', 'neverworld2', 'island', 'double_drake', 'isomip_plus', 'file'),
    )

    max_depth = Real(
        'max_depth',
        doc='Basin maximum depth',
        units='m',
        required=False,
        default=2000.0,
    )

    edge_depth = Real(
        'edge_depth',
        doc='Spoon-bathymetry edge depth',
        units='m',
        required=False,
        default=100.0,
    )

    slope_scale = Real(
        'slope_scale',
        doc='Spoon-bathymetry exponential decay scale',
        units='m',
        required=False,
        default=400000.0,
    )

    nl_continent_amp = Real(
        'nl_continent_amp',
        doc='Neverworld2 continent amplitude (1=full continents, 0=aquaplanet+channel)',
        units='nondim',
        required=False,
        default=1.0,
        has_min=True,
        vmin=0.0,
    )

    nl_roughness_amp = Real(
        'nl_roughness_amp',
        doc='Neverworld2 bathymetry roughness amplitude',
        units='nondim',
        required=False,
        default=0.05,
        has_min=True,
        vmin=0.0,
    )

    nl_min_depth = Real(
        'nl_min_depth',
        doc='Neverworld2 minimum-depth floor (MOM6 MINIMUM_DEPTH analogue)',
        units='m',
        required=False,
        default=500.0,
        has_min=True,
        vmin=0.0,
    )

    wind_config = Enum(
        'wind_config',
        doc='Surface wind-stress dispatch',
        units='',
        required=False,
        default='constant',
        allowed=('constant', '2gyre', 'neverworld2'),
    )

    taux_magnitude = Real(
        'taux_magnitude',
        doc='Peak zonal wind stress for wind_config=2gyre/neverworld2',
        units='Pa',
        required=False,
        default=0.1,
    )

    coriolis_beta = Real(
        'coriolis_beta',
        doc='Meridional gradient of f (0 = f-plane)',
        units='1/(s m)',
        required=False,
        default=0.0,
    )

    coriolis_y_ref = Real(
        'coriolis_y_ref',
        doc='Reference y where f = coriolis_f under beta-plane',
        units='m',
        required=False,
        default=0.0,
    )

    x_origin = Real(
        'x_origin',
        doc="Absolute x of the domain west edge, for topo_config='isomip_plus' (ISOMIP+ ocean box starts at the MISMIP+ x = 320 km)",
        units='m',
        required=False,
        default=0.0,
    )


class OceanIc(Group):
    """`&ocean_ic_nml` -- Initial-condition overlay + EOS reference state."""

    _nml_name = 'ocean_ic'

    ic_config = Enum(
        'ic_config',
        doc='IC overlay tag',
        units='',
        required=False,
        default='',
        allowed=('', 'eady', 'geostrophic_adjustment', 'baroclinic_jet'),
    )

    alpha_T = Real(
        'alpha_T',
        doc='Linear-EOS thermal-expansion coefficient (DIMENSIONAL: multiply a fractional 1/degC coefficient by rho_0)',
        units='kg/m^3/degC',
        required=False,
        default=0.00017,
    )

    beta_S = Real(
        'beta_S',
        doc='Linear-EOS haline contraction coefficient (DIMENSIONAL: multiply a fractional 1/PSU coefficient by rho_0)',
        units='kg/m^3/PSU',
        required=False,
        default=0.00076,
        has_min=True,
        vmin=0.0,
    )

    T_ref = Real(
        'T_ref',
        doc='Linear-EOS reference temperature',
        units='degC',
        required=False,
        default=10.0,
        has_min=True,
        vmin=-273.15,
    )

    S_ref = Real(
        'S_ref',
        doc='Linear-EOS reference salinity',
        units='PSU',
        required=False,
        default=35.0,
        has_min=True,
        vmin=0.0,
    )

    rho_0 = Real(
        'rho_0',
        doc='Reference density for the linear EOS / Boussinesq PGF',
        units='kg/m^3',
        required=False,
        default=1035.0,
    )

    layer_rho_init = RealArray(
        'layer_rho_init',
        doc='Per-layer initial density (k=1 bed -> k=nz surface; -1 = EOS init)',
        units='kg/m^3',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=64,
    )

    rho_lightest = Real(
        'rho_lightest',
        doc='Linear density-range IC: surface (lightest) layer density; -1 = off',
        units='kg/m^3',
        required=False,
        default=-1.0,
    )

    rho_range = Real(
        'rho_range',
        doc='Linear density-range IC: top-to-bottom density contrast (MOM6 DENSITY_RANGE)',
        units='kg/m^3',
        required=False,
        default=2.0,
    )

    eady_dT_dy = Real(
        'eady_dT_dy',
        doc='Eady IC meridional T gradient',
        units='degC/m',
        required=False,
        default=-2e-05,
    )

    eady_dT_dz = Real(
        'eady_dT_dz',
        doc='Eady IC vertical stratification',
        units='degC/m',
        required=False,
        default=0.01,
    )

    eady_T_ref = Real(
        'eady_T_ref',
        doc='Eady IC reference temperature at z=0,y=y_mid',
        units='degC',
        required=False,
        default=10.0,
    )

    eady_pert_amp = Real(
        'eady_pert_amp',
        doc='Eady IC symmetry-breaking perturbation amplitude',
        units='degC',
        required=False,
        default=0.001,
    )

    eady_pert_seed = Int(
        'eady_pert_seed',
        doc='RNG seed for the Eady IC perturbation',
        units='',
        required=False,
        default=12345,
    )

    ga_eta_amp = Real(
        'ga_eta_amp',
        doc='Geostrophic-adjustment IC: SSH-bump amplitude',
        units='m',
        required=False,
        default=1.0,
    )

    ga_length_scale = Real(
        'ga_length_scale',
        doc='Geostrophic-adjustment IC: SSH-bump e-folding scale',
        units='m',
        required=False,
        default=50000.0,
    )

    ga_x_center = Real(
        'ga_x_center',
        doc='Geostrophic-adjustment IC: bump x-centre (<0 = auto)',
        units='m',
        required=False,
        default=-1.0,
    )

    ga_y_center = Real(
        'ga_y_center',
        doc='Geostrophic-adjustment IC: bump y-centre (<0 = auto)',
        units='m',
        required=False,
        default=-1.0,
    )

    jet_half_width = Real(
        'jet_half_width',
        doc='Baroclinic-jet IC: tanh jet half-width L',
        units='m',
        required=False,
        default=40000.0,
    )

    interface_amp = Real(
        'interface_amp',
        doc='Baroclinic-jet IC: interface displacement amplitude',
        units='m',
        required=False,
        default=200.0,
    )

    pert_amp_frac = Real(
        'pert_amp_frac',
        doc='Baroclinic-jet IC: meander amplitude as a fraction of interface_amp',
        units='',
        required=False,
        default=0.2,
    )

    pert_nx = Int(
        'pert_nx',
        doc='Baroclinic-jet IC: zonal perturbation wavenumber (integer)',
        units='',
        required=False,
        default=3,
    )

    upper_layer_rest = Real(
        'upper_layer_rest',
        doc='Baroclinic-jet IC: upper (surface) layer rest thickness H1',
        units='m',
        required=False,
        default=500.0,
    )


class OceanZinit(Group):
    """`&ocean_zinit_nml` -- Z-level T/S initial-condition overlay (A2)."""

    _nml_name = 'ocean_zinit'

    enable = Bool(
        'enable',
        doc='Master switch (default off; requires RDB_ENABLE_NETCDF=ON)',
        units='',
        required=False,
        default=False,
    )

    source = Enum(
        'source',
        doc='Where the T(z)/S(z) profile comes from: a pre-regridded NetCDF, or the analytic affine lin_* profile (no file)',
        units='',
        required=False,
        default='file',
        allowed=('file', 'linear'),
    )

    file = Str(
        'file',
        doc='Path to the model-grid T/S NetCDF',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    t_var = Str(
        't_var',
        doc='Temperature variable-name override (blank tries temp/T/temperature)',
        units='',
        required=False,
        default='',
        max_len=32,
    )

    s_var = Str(
        's_var',
        doc='Salinity variable-name override (blank tries salt/S/salinity)',
        units='',
        required=False,
        default='',
        max_len=32,
    )

    z_var = Str(
        'z_var',
        doc='Source-axis variable-name override (blank tries z_src/z/depth/lev)',
        units='',
        required=False,
        default='',
        max_len=32,
    )

    land_fill_t = Real(
        'land_fill_t',
        doc='Fallback temperature for dry columns',
        units='degC',
        required=False,
        default=10.0,
    )

    land_fill_s = Real(
        'land_fill_s',
        doc='Fallback salinity for dry columns',
        units='PSU',
        required=False,
        default=35.0,
    )

    lin_t_ref = Real(
        'lin_t_ref',
        doc="source='linear': temperature at the z = 0 datum",
        units='degC',
        required=False,
        default=0.0,
    )

    lin_dt_dz = Real(
        'lin_dt_dz',
        doc="source='linear': dT/dz, z positive UP (stable > 0)",
        units='degC/m',
        required=False,
        default=0.0,
    )

    lin_s_ref = Real(
        'lin_s_ref',
        doc="source='linear': salinity at the z = 0 datum",
        units='PSU',
        required=False,
        default=35.0,
    )

    lin_ds_dz = Real(
        'lin_ds_dz',
        doc="source='linear': dS/dz, z positive UP (stable < 0)",
        units='PSU/m',
        required=False,
        default=0.0,
    )


class OceanData(Group):
    """`&ocean_data_nml` -- Shared time-varying NetCDF input reader (PR-14): registry sizing only."""

    _nml_name = 'ocean_data'

    max_fields = Int(
        'max_fields',
        doc="Size of the reader's field registry",
        units='',
        required=False,
        default=16,
        has_min=True,
        vmin=1,
    )

    verbose = Bool(
        'verbose',
        doc='Log every bracket advance (field, records, weight, time)',
        units='',
        required=False,
        default=False,
    )


class OceanDataovr(Group):
    """`&ocean_dataovr_nml` -- File-backed surface forcing (PR-15): wind stress, heat, freshwater."""

    _nml_name = 'ocean_dataovr'

    enable = Bool(
        'enable',
        doc='Master switch (default off; requires RDB_ENABLE_NETCDF=ON)',
        units='',
        required=False,
        default=False,
    )

    time_mode = Enum(
        'time_mode',
        doc='Shared time-axis mode for every active tag',
        units='',
        required=False,
        default='linear',
        allowed=('linear', 'cyclic', 'static'),
    )

    cycle_period = Real(
        'cycle_period',
        doc='Climatology period (required > 0 when time_mode=cyclic)',
        units='s',
        required=False,
        default=0.0,
        has_min=True,
        vmin=0.0,
    )

    t_offset = Real(
        't_offset',
        doc='Added to the model time before the file-axis lookup',
        units='s',
        required=False,
        default=0.0,
    )

    oor_clamp = Bool(
        'oor_clamp',
        doc='Clamp an out-of-range query to the end record instead of aborting',
        units='',
        required=False,
        default=False,
    )

    tau_x_file = Str(
        'tau_x_file',
        doc='Path to the zonal wind stress NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    tau_x_var = Str(
        'tau_x_var',
        doc='Variable name inside tau_x_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    tau_x_scale = Real(
        'tau_x_scale',
        doc='Multiplier applied to the zonal wind stress slab at read',
        units='Pa',
        required=False,
        default=1.0,
    )

    tau_x_add = Real(
        'tau_x_add',
        doc='Offset added to the zonal wind stress slab after scale',
        units='Pa',
        required=False,
        default=0.0,
    )

    tau_y_file = Str(
        'tau_y_file',
        doc='Path to the meridional wind stress NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    tau_y_var = Str(
        'tau_y_var',
        doc='Variable name inside tau_y_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    tau_y_scale = Real(
        'tau_y_scale',
        doc='Multiplier applied to the meridional wind stress slab at read',
        units='Pa',
        required=False,
        default=1.0,
    )

    tau_y_add = Real(
        'tau_y_add',
        doc='Offset added to the meridional wind stress slab after scale',
        units='Pa',
        required=False,
        default=0.0,
    )

    heat_file = Str(
        'heat_file',
        doc='Path to the net surface heat flux (positive down) NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    heat_var = Str(
        'heat_var',
        doc='Variable name inside heat_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    heat_scale = Real(
        'heat_scale',
        doc='Multiplier applied to the net surface heat flux (positive down) slab at read',
        units='W m-2',
        required=False,
        default=1.0,
    )

    heat_add = Real(
        'heat_add',
        doc='Offset added to the net surface heat flux (positive down) slab after scale',
        units='W m-2',
        required=False,
        default=0.0,
    )

    evap_file = Str(
        'evap_file',
        doc='Path to the evaporative mass flux (<= 0) NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    evap_var = Str(
        'evap_var',
        doc='Variable name inside evap_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    evap_scale = Real(
        'evap_scale',
        doc='Multiplier applied to the evaporative mass flux (<= 0) slab at read',
        units='kg m-2 s-1',
        required=False,
        default=1.0,
    )

    evap_add = Real(
        'evap_add',
        doc='Offset added to the evaporative mass flux (<= 0) slab after scale',
        units='kg m-2 s-1',
        required=False,
        default=0.0,
    )

    lprec_file = Str(
        'lprec_file',
        doc='Path to the liquid precipitation (>= 0) NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    lprec_var = Str(
        'lprec_var',
        doc='Variable name inside lprec_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    lprec_scale = Real(
        'lprec_scale',
        doc='Multiplier applied to the liquid precipitation (>= 0) slab at read',
        units='kg m-2 s-1',
        required=False,
        default=1.0,
    )

    lprec_add = Real(
        'lprec_add',
        doc='Offset added to the liquid precipitation (>= 0) slab after scale',
        units='kg m-2 s-1',
        required=False,
        default=0.0,
    )

    salt_file = Str(
        'salt_file',
        doc='Path to the surface salt flux (positive salinifies) NetCDF (blank ⇒ tag not file-driven)',
        units='',
        required=False,
        default='',
        max_len=256,
    )

    salt_var = Str(
        'salt_var',
        doc='Variable name inside salt_file (required when it is set)',
        units='',
        required=False,
        default='',
        max_len=64,
    )

    salt_scale = Real(
        'salt_scale',
        doc='Multiplier applied to the surface salt flux (positive salinifies) slab at read',
        units='kg m-2 s-1',
        required=False,
        default=1.0,
    )

    salt_add = Real(
        'salt_add',
        doc='Offset added to the surface salt flux (positive salinifies) slab after scale',
        units='kg m-2 s-1',
        required=False,
        default=0.0,
    )


class OceanDiag(Group):
    """`&ocean_diag_nml` -- Ocean diag-manager output controls."""

    _nml_name = 'ocean_diag'

    enabled = Bool(
        'enabled',
        doc='Enable the per-step diag-manager hook',
        units='',
        required=False,
        default=True,
    )

    filename = Str(
        'filename',
        doc='Output file basename (per-rank suffix appended)',
        units='',
        required=False,
        default='ocean_diag',
        max_len=256,
    )

    dt_out = Real(
        'dt_out',
        doc='Diag-fire cadence (in time_unit)',
        units='',
        required=False,
        default=3600.0,
    )

    vgrid = Enum(
        'vgrid',
        doc='Default output vertical grid for layer-shaped diags',
        units='',
        required=False,
        default='layer',
        allowed=('layer', 'z_fixed', 'sigma', 'zstar', 'density'),
    )

    z_levels = RealArray(
        'z_levels',
        doc='Output z-levels when vgrid=z_fixed (positive down)',
        units='m',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=64,
    )

    n_z_levels = Int(
        'n_z_levels',
        doc='Number of z_levels entries used (0 = none)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    sigma_levels = RealArray(
        'sigma_levels',
        doc='Output sigma fractions (0..1) for sigma output; empty = auto-uniform',
        units='',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=64,
    )

    n_sigma_levels = Int(
        'n_sigma_levels',
        doc='Number of sigma_levels entries used (0 = auto)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    zstar_levels = RealArray(
        'zstar_levels',
        doc='Output z* reference depths for zstar output; empty = auto-uniform',
        units='m',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=64,
    )

    n_zstar_levels = Int(
        'n_zstar_levels',
        doc='Number of zstar_levels entries used (0 = auto)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    rho_levels = RealArray(
        'rho_levels',
        doc='Output potential-density bin edges when vgrid=density; strictly increasing, light->dense',
        units='kg/m^3',
        required=False,
        default=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
        size=64,
    )

    n_rho_levels = Int(
        'n_rho_levels',
        doc='Number of rho_levels entries used (0 = none)',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
    )

    mask_vanished_layers = Bool(
        'mask_vanished_layers',
        doc='Mask below-bottom/pinched remap cells to missing_value (default off)',
        units='',
        required=False,
        default=False,
    )

    reproducing_sums = Bool(
        'reproducing_sums',
        doc='Order-invariant EFP console totals: rank-count reproducible, exact drift residual (default off = byte-identical)',
        units='',
        required=False,
        default=False,
    )

    output_precision = Enum(
        'output_precision',
        doc="Element width of the diag NetCDF data vars; 'single' halves the bytes written (restarts/gauges/console totals stay double regardless)",
        units='',
        required=False,
        default='double',
        allowed=('double', 'single'),
    )

    diag_remap_scheme = Enum(
        'diag_remap_scheme',
        doc='Reconstruction for the conservative diagnostic vertical remap',
        units='',
        required=False,
        default='ppm',
        allowed=('pcm', 'plm', 'ppm', 'ppm_h4', 'pqm'),
    )

    diags = Str(
        'diags',
        doc="Unified diag list modifying the default set: name[:off][:cadence][:op] (e.g. 'vorticity_z:1d  KE:off  temperature:6h:mean'); default empty = canonical set",
        units='',
        required=False,
        default='',
        max_len=512,
    )


class OceanBc(Group):
    """`&ocean_bc_nml` -- Open-boundary condition config: per-edge BC type, clamped/inflow Dirichlet values, sponge band, Orlanski radiation, full Flather, per-edge tidal constituents, tidal-OBC nodal correction, wall-velocity masking."""

    _nml_name = 'ocean_bc'

    west = Enum(
        'west',
        doc='West edge boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman', 'periodic', 'tripolar_fold'),
    )

    east = Enum(
        'east',
        doc='East edge boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman', 'periodic', 'tripolar_fold'),
    )

    south = Enum(
        'south',
        doc='South edge boundary type',
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman', 'periodic', 'tripolar_fold'),
    )

    north = Enum(
        'north',
        doc="North edge boundary type ('tripolar_fold' is meaningful here only — Murray bipolar cap; requires grid_config='tripolar' + periodic west/east)",
        units='',
        required=False,
        default='wall',
        allowed=('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge', 'clamped', 'sponge', 'chapman', 'periodic', 'tripolar_fold'),
    )

    west_clamped_eta = Real(
        'west_clamped_eta',
        doc='Clamped SSH, west edge',
        units='m',
        required=False,
        default=0.0,
    )

    east_clamped_eta = Real(
        'east_clamped_eta',
        doc='Clamped SSH, east edge',
        units='m',
        required=False,
        default=0.0,
    )

    south_clamped_eta = Real(
        'south_clamped_eta',
        doc='Clamped SSH, south edge',
        units='m',
        required=False,
        default=0.0,
    )

    north_clamped_eta = Real(
        'north_clamped_eta',
        doc='Clamped SSH, north edge',
        units='m',
        required=False,
        default=0.0,
    )

    west_clamped_u = Real(
        'west_clamped_u',
        doc='Clamped normal velocity, west edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    east_clamped_u = Real(
        'east_clamped_u',
        doc='Clamped normal velocity, east edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    south_clamped_v = Real(
        'south_clamped_v',
        doc='Clamped normal velocity, south edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    north_clamped_v = Real(
        'north_clamped_v',
        doc='Clamped normal velocity, north edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    west_inflow_S = Real(
        'west_inflow_S',
        doc='Inflow salinity, west edge',
        units='PSU',
        required=False,
        default=35.0,
    )

    west_inflow_T = Real(
        'west_inflow_T',
        doc='Inflow temperature, west edge',
        units='degC',
        required=False,
        default=10.0,
    )

    east_inflow_S = Real(
        'east_inflow_S',
        doc='Inflow salinity, east edge',
        units='PSU',
        required=False,
        default=35.0,
    )

    east_inflow_T = Real(
        'east_inflow_T',
        doc='Inflow temperature, east edge',
        units='degC',
        required=False,
        default=10.0,
    )

    south_inflow_S = Real(
        'south_inflow_S',
        doc='Inflow salinity, south edge',
        units='PSU',
        required=False,
        default=35.0,
    )

    south_inflow_T = Real(
        'south_inflow_T',
        doc='Inflow temperature, south edge',
        units='degC',
        required=False,
        default=10.0,
    )

    north_inflow_S = Real(
        'north_inflow_S',
        doc='Inflow salinity, north edge',
        units='PSU',
        required=False,
        default=35.0,
    )

    north_inflow_T = Real(
        'north_inflow_T',
        doc='Inflow temperature, north edge',
        units='degC',
        required=False,
        default=10.0,
    )

    sponge_width = Int(
        'sponge_width',
        doc='Sponge band width in cells (0 disables)',
        units='',
        required=False,
        default=0,
    )

    sponge_strength = Real(
        'sponge_strength',
        doc='Peak sponge relaxation rate at the outer face of the band',
        units='1/s',
        required=False,
        default=0.0,
    )

    sponge_relax_tracers = Bool(
        'sponge_relax_tracers',
        doc='Extend the sponge relaxation to h_layer + tracers',
        units='',
        required=False,
        default=False,
    )

    res_lscale_out = Real(
        'res_lscale_out',
        doc='Outflow reservoir length scale (0 = disabled)',
        units='m',
        required=False,
        default=0.0,
    )

    res_lscale_in = Real(
        'res_lscale_in',
        doc='Inflow reservoir length scale (0 = instantaneous inflow)',
        units='m',
        required=False,
        default=0.0,
    )

    radiation_scheme = Enum(
        'radiation_scheme',
        doc="'anomaly' = v1 BT-mean + zero-gradient anomaly (default); 'orlanski' = per-layer implicit-upwind radiation (Orlanski 1976)",
        units='',
        required=False,
        default='anomaly',
        allowed=('anomaly', 'orlanski'),
    )

    orlanski_rx_max = Real(
        'orlanski_rx_max',
        doc='Clamp on the Orlanski nondimensional phase speed',
        units='grid cells/step',
        required=False,
        default=10.0,
    )

    orlanski_gamma = Real(
        'orlanski_gamma',
        doc='Running-mean weight (0=full running mean, 1=instant)',
        units='',
        required=False,
        default=1.0,
    )

    nudge_tau_in = Real(
        'nudge_tau_in',
        doc='Inflow nudging timescale (0 = off)',
        units='s',
        required=False,
        default=0.0,
    )

    nudge_tau_out = Real(
        'nudge_tau_out',
        doc='Outflow nudging timescale (0 = off)',
        units='s',
        required=False,
        default=0.0,
    )

    flather_form = Enum(
        'flather_form',
        doc="'legacy' = v1 Flather (u_ext=0, no interior vel, default); 'full' = half-characteristic form (Flather 1976)",
        units='',
        required=False,
        default='legacy',
        allowed=('legacy', 'full'),
    )

    west_ext_u = Real(
        'west_ext_u',
        doc='Exterior barotropic u, west edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    east_ext_u = Real(
        'east_ext_u',
        doc='Exterior barotropic u, east edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    south_ext_v = Real(
        'south_ext_v',
        doc='Exterior barotropic v, south edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    north_ext_v = Real(
        'north_ext_v',
        doc='Exterior barotropic v, north edge',
        units='m/s',
        required=False,
        default=0.0,
    )

    west_n_tidal = Int(
        'west_n_tidal',
        doc='Active tidal constituents, west edge',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
        has_max=True,
        vmax=8,
    )

    west_tidal_amp = RealArray(
        'west_tidal_amp',
        doc='Constituent amplitudes, west edge',
        units='m',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    west_tidal_phase = RealArray(
        'west_tidal_phase',
        doc='Constituent phases, west edge',
        units='rad',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    west_tidal_omega = RealArray(
        'west_tidal_omega',
        doc='Constituent angular frequencies, west edge',
        units='rad/s',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    east_n_tidal = Int(
        'east_n_tidal',
        doc='Active tidal constituents, east edge',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
        has_max=True,
        vmax=8,
    )

    east_tidal_amp = RealArray(
        'east_tidal_amp',
        doc='Constituent amplitudes, east edge',
        units='m',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    east_tidal_phase = RealArray(
        'east_tidal_phase',
        doc='Constituent phases, east edge',
        units='rad',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    east_tidal_omega = RealArray(
        'east_tidal_omega',
        doc='Constituent angular frequencies, east edge',
        units='rad/s',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    south_n_tidal = Int(
        'south_n_tidal',
        doc='Active tidal constituents, south edge',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
        has_max=True,
        vmax=8,
    )

    south_tidal_amp = RealArray(
        'south_tidal_amp',
        doc='Constituent amplitudes, south edge',
        units='m',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    south_tidal_phase = RealArray(
        'south_tidal_phase',
        doc='Constituent phases, south edge',
        units='rad',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    south_tidal_omega = RealArray(
        'south_tidal_omega',
        doc='Constituent angular frequencies, south edge',
        units='rad/s',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    north_n_tidal = Int(
        'north_n_tidal',
        doc='Active tidal constituents, north edge',
        units='',
        required=False,
        default=0,
        has_min=True,
        vmin=0,
        has_max=True,
        vmax=8,
    )

    north_tidal_amp = RealArray(
        'north_tidal_amp',
        doc='Constituent amplitudes, north edge',
        units='m',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    north_tidal_phase = RealArray(
        'north_tidal_phase',
        doc='Constituent phases, north edge',
        units='rad',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    north_tidal_omega = RealArray(
        'north_tidal_omega',
        doc='Constituent angular frequencies, north edge',
        units='rad/s',
        required=False,
        default=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        size=8,
    )

    obc_tidal_nodal = Bool(
        'obc_tidal_nodal',
        doc='Apply the 18.6-yr nodal factor + equilibrium/nodal phase to the open-boundary tidal elevation forcing (shares the &ocean_tides_nml astro generator); when true *_tidal_phase becomes a Greenwich phase lag',
        units='',
        required=False,
        default=False,
    )

    mask_wall_velocity = Bool(
        'mask_wall_velocity',
        doc='Zero the T-cell wet-mask in the ghost cells beyond every solid WALL edge at setup, so wall-normal C-grid face velocities are masked rather than left as garbage (default ON; set .false. to reproduce a pre-fix legacy closed-basin baseline)',
        units='',
        required=False,
        default=True,
    )


#: schema group name -> generated Group subclass, in the order build_rdb_schema
#: registers them. Combined with OceanBc (python/rdb/_config_bc.py)
#: by python/rdb/_config.py to build the full Config object.
GENERATED_GROUPS = {
    'sim': Sim,
    'grid': Grid,
    'time': Time,
    'mpi': Mpi,
    'logging': Logging,
    'vcoord': Vcoord,
    'physics': Physics,
    'output': Output,
    'boundary': Boundary,
    'tracer': Tracer,
    'initial_condition': InitialCondition,
    'nonhydrostatic': Nonhydrostatic,
    'ocean_kappa_shear': OceanKappaShear,
    'ocean_slopes': OceanSlopes,
    'ocean_gm': OceanGm,
    'ocean_redi': OceanRedi,
    'ocean_varmix': OceanVarmix,
    'ocean_meke': OceanMeke,
    'ocean_tidal_mixing': OceanTidalMixing,
    'ocean_conv': OceanConv,
    'ocean_porous': OceanPorous,
    'ocean_ddiff': OceanDdiff,
    'ocean_tides': OceanTides,
    'ocean_psurf': OceanPsurf,
    'ocean_cavity_dyn': OceanCavityDyn,
    'ocean_cavity_melt': OceanCavityMelt,
    'ocean_epbl': OceanEpbl,
    'ocean_wavespeed': OceanWavespeed,
    'ocean_foxkemper': OceanFoxkemper,
    'ocean_grid': OceanGrid,
    'ocean_coriolis': OceanCoriolis,
    'ocean_thermo': OceanThermo,
    'ocean_forcing': OceanForcing,
    'ocean_ice': OceanIce,
    'ocean_ice_ic': OceanIceIc,
    'ocean_restore': OceanRestore,
    'ocean_geothermal': OceanGeothermal,
    'ocean_sponge': OceanSponge,
    'ocean_tracers': OceanTracers,
    'ocean_bt': OceanBt,
    'ocean_debug': OceanDebug,
    'ocean_mpi': OceanMpi,
    'ocean_wetdry': OceanWetdry,
    'ocean_pgf': OceanPgf,
    'ocean_eos': OceanEos,
    'ocean_bdrag': OceanBdrag,
    'ocean_tdrag': OceanTdrag,
    'ocean_hdiff': OceanHdiff,
    'ocean_hvisc': OceanHvisc,
    'ocean_vmix': OceanVmix,
    'ocean_vdiff': OceanVdiff,
    'ocean_continuity': OceanContinuity,
    'ocean_isopycnal': OceanIsopycnal,
    'ocean_topo': OceanTopo,
    'ocean_ic': OceanIc,
    'ocean_zinit': OceanZinit,
    'ocean_data': OceanData,
    'ocean_dataovr': OceanDataovr,
    'ocean_diag': OceanDiag,
    'ocean_bc': OceanBc,
}

N_GROUPS = 60
N_KNOBS = 670
