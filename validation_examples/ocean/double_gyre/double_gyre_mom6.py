"""Roundabout ocean-path double-gyre -- MOM6 reference reproduction,
Python-API port.

Adapted from ``double_gyre_mom6.nml``'s own header comment:

Mirrors NOAA-GFDL/MOM6-examples/ocean_only/double_gyre exactly up to the
limits of the current ocean stack:

    NIGLOBAL = 44, NJGLOBAL = 40
    22 deg lon x 20 deg lat at mid-lat ~40N -> Cartesian dx/dy
    NK = 2, MAX_DEPTH = 2000 m, spoon bathymetry (D_edge=100, slope=400 km)
    TAUX_MAGNITUDE = 0.1 Pa, 2gyre wind profile
    DT = 1200 s, DAYMAX = 10 days
    Beta-plane: f_0 = 9.4e-5 (40N), beta = 1.76e-11

Known caveats vs MOM6 (see ``README.md``):

  - MOM6 uses gprime / adiabatic reduced gravity, NK=2, with
    ENABLE_THERMODYNAMICS = False. Roundabout's multilayer kernel always
    runs T/S advection, so we instead set alpha_T = beta_S = 0 below --
    the linear EOS returns rho_0 regardless of T/S so any
    tracer-advection drift cannot reach the dynamics. SSH amplitudes
    will not bit-match MOM6, but the mass + energy budget is decoupled
    from thermodynamics in the same spirit.
  - MOM6 BOUND_CORIOLIS = True (Sadourny + HK correction). Roundabout ships
    the Arakawa-Hsu / HK kernel -- engaged below via
    ocean_coriolis_form = "sadourny_energy" (matches MOM6's
    CORIOLIS_SCHEME = "SADOURNY75_ENERGY" default for THIS config).
  - MOM6 LINEAR_DRAG with HBBL=10m distributes the linear drag across
    the bottom HBBL metres. Roundabout's bottom drag applies to the bed
    layer only -- same Rayleigh form but a different BBL footprint.

THE FORMERLY MOST SIGNIFICANT ESCAPE HATCH IN THIS PORT, NOW CLOSED:
MOM6's "degrees-per-cell, arc-length-derived Cartesian" grid mode --
``&ocean_grid_nml axis_units = "degrees"`` on an otherwise
``grid_config = "cartesian"`` grid, with ``dx``/``dy`` DERIVED from
``len_lon``/``len_lat``/``rad_earth`` as arc lengths (no ``cos(lat)``
factor) -- used to have NO curated route at all: ``RectilinearGrid``
already accepted ``axis_units="degrees"`` (it forwards whatever
``extent=`` it is given straight to ``len_lon``/``len_lat``
regardless of unit), but had no way to reach ``&ocean_grid_nml
rad_earth``, which the "degrees" derivation needs. ``RectilinearGrid``
now takes ``radius=`` for exactly this; see ``build_model()`` below.

Note ``extent=(22.0, 20.0)`` on the ``RectilinearGrid`` below is in
DEGREES here, not metres -- ``RectilinearGrid``'s docstring is explicit
that ``extent=``'s units follow ``axis_units=``.

Run:
    cd validation_examples/ocean/double_gyre
    # first load your NVHPC + NetCDF toolchain (module load / spack / conda;
    # see environments/) -- one toolchain per shell, never stack two
    CUDA_VISIBLE_DEVICES=<gpu> python3 double_gyre_mom6.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make `rdb` importable when run directly from any cwd -- see
# acc_channel.py's matching bootstrap for the full rationale.
_REPO_ROOT = Path(__file__).resolve().parents[3]
_PY_DIR = str(_REPO_ROOT / "python")
if _PY_DIR not in sys.path:
    sys.path.insert(0, _PY_DIR)

import rdb
from rdb.barotropic import BarotropicSolver
from rdb.closures import (
    HarmonicFaceThickness,
    LinearDrag,
    NearSurfaceViscosity,
    Smagorinsky,
    ThermoSubcycling,
)
from rdb.coriolis import BetaPlane, Sadourny
from rdb.diagnostics import Diagnostics, Duration
from rdb.forcing import (
    DoubleGyreWind,
    LayerDensities,
    Thermodynamics,
    UniformInitialCondition,
)
from rdb.grids import RectilinearGrid, Spoon, Topology
from rdb.pressure import ReducedGravity, VelocityTruncation
from rdb.vcoord import ZStarFull


def build_model():
    grid = RectilinearGrid(
        size=(44, 40, 2),                  # &grid_nml nx/ny, &nonhydro-
                                            # static_nml nz_layers
        extent=(22.0, 20.0),                # DEGREES here -- see the
                                             # module docstring: extent=
                                             # follows axis_units=
        axis_units="degrees",               # MOM6-parity arc-length dx/dy
        radius=6.378e6,                     # &ocean_grid_nml rad_earth
        topology=(Topology.BOUNDED, Topology.BOUNDED),
        halo=3,   # matches the .nml's explicit &grid_nml nghost=3 (pinned
                  # explicitly, matching seamount.py's convention, rather
                  # than relying on the derivation landing on 3 too).
    )

    model = rdb.Model(
        grid=grid,
        timestep=1200.0,   # &time_nml dt_fixed

        coriolis=BetaPlane(f0=9.4e-5, beta=1.76e-11, y_ref=1113200.0),
        momentum_advection=Sadourny(form="energy"),  # -> sadourny_energy
        # No buoyancy= object here, DELIBERATELY: the .nml never touches
        # &ocean_eos_nml (default "linear" already matches) or
        # &ocean_ic_nml alpha_T/rho_0 (stays at the Fortran's own
        # defaults) -- and, per the .nml's own comment, ocean_ic.alpha_T
        # is unused anyway once enable_thermodynamics=.false. below
        # skips the EOS call entirely. The "set alpha_T = beta_S = 0"
        # belt-and-braces the .nml describes is the COASTAL-LEGACY
        # &tracer_nml pair (dead_on_ocean_path in the schema), not this
        # one -- see the tracer.alpha_T/beta_S escape-hatch lines below.
        # An earlier draft of this script called
        # LinearEquationOfState(thermal_expansion=0.0, ...) here, which
        # wrote a NON-default &ocean_ic_nml alpha_T=0.0 the .nml never
        # sets -- caught by the config-identity gate (see the port
        # report); removed rather than papered over.
        pressure=ReducedGravity(g_free_surface=0.98, g_internal=0.0098),

        vcoord=ZStarFull(h_surf_target=1000.0, h_min=1.5e-4),

        closures=[
            Smagorinsky(C=0.15, biharmonic=True, C_biharmonic=0.06,
                        free_slip=True,  # &ocean_hvisc_nml no_slip is
                                          # False by Fortran default and
                                          # the .nml never touches it --
                                          # see acc_channel.py's fuller
                                          # note on this same gap.
                        nu_h=10000.0, kh_vel_scale=3.0e-3),
            LinearDrag(r=2.5e-5, hbbl=10.0, bg_vel=0.1,
                       bbl_thick_min=0.1),
            # MOM6 KV_ML_INVZ2 near-surface momentum viscosity band --
            # composes alongside whatever (if any) boundary-layer scheme
            # is active; this config has neither KPP nor EPBL.
            NearSurfaceViscosity(kv=0.01, hmix=20.0),
            # MOM6 HARMONIC_VISC vdiff-assembly parity.
            HarmonicFaceThickness(),
            # MOM6 DT_THERM=2*DT sub-cycling ratio.
            ThermoSubcycling(ratio=2),
        ],

        forcing=[
            DoubleGyreWind(magnitude=0.1),
            # MOM6 ENABLE_THERMODYNAMICS = False analogue -- skips the
            # EOS call entirely (see the module docstring's caveat on
            # why the reduced-gravity form is used instead of a full
            # EOS-driven baroclinic run).
            Thermodynamics(enabled=False),
        ],

        initial_condition=[
            UniformInitialCondition(temperature=15.0, salinity=35.0),
            LayerDensities([1036.0, 1035.0]),  # k=1 (bed) first, k=nz last
        ],

        barotropic=BarotropicSolver(
            auto=True, cfl_safety=0.65, bebt=0.2),

        diagnostics=Diagnostics(
            enabled=True, filename="double_gyre", every=1.0,
            output_dir="./out_double_gyre_mom6", status_interval=1.0),

        duration=Duration(t_end=10.0, time_unit="day"),
    )

    cfg = model.config

    # Bathymetry: Spoon() IS a curated writer.
    Spoon(max_depth=2000.0, edge_depth=100.0, slope_scale=400000.0).apply(cfg)

    # PGF safety clamp: no maxvel param on ReducedGravity (a separate
    # object per D2.6).
    VelocityTruncation(maxvel=6.0).apply(cfg)

    # NOTE: &ocean_vmix_nml direct_stress/hmix_stress and &ocean_hvisc_nml
    # ah_max are NOT set here -- the .nml's own values (direct_stress=
    # False, hmix_stress=20.0, ah_max=1e4) are already the Fortran's OWN
    # schema defaults (and Smagorinsky()'s own ah_max default), so
    # setting them explicitly would be a behaviourally-inert no-op; see
    # the port report for the full accounting.

    # Belt-and-braces tracer drift guard (see .nml comment): both knobs
    # are dead_on_ocean_path per the generated schema's own docstring
    # (tracer.alpha_T/beta_S feed a coastal-legacy path read nowhere on
    # the ocean path, per the P4 dead-knob sweep) -- reproduced verbatim
    # anyway for configuration fidelity, not because it does anything.
    # DELIBERATELY left on the escape hatch: giving a dead knob a nice
    # curated route would misrepresent it as meaningful.
    cfg.tracer.alpha_T = 0.0
    cfg.tracer.beta_S = 0.0

    # &output_nml output_to_file: cosmetic-only dead knob, see
    # acc_channel.py's matching note (Diagnostics(enabled=True) sets it
    # True; the .nml sets it False independently; kept for text
    # fidelity, not because it changes anything).
    cfg.output.output_to_file = False

    return model


def main():
    model = build_model()
    print(model.config.to_namelist()[:400], "...\n")

    with model:
        ssh_min_series = []

        def daily(m):
            eta = m.ssh  # nested list, (nx, ny)
            ssh_min = min(min(row) for row in eta)
            ssh_min_series.append(ssh_min)
            print(f"day {m.time / 86400.0:5.2f}  SSH_min {ssh_min:+.4f} m  "
                  f"KE {m.kinetic_energy:.4e} J")

        model.run(until=10 * 86400.0, callback=daily, every=86400.0)

        print(f"\nday-10 SSH minimum: {ssh_min_series[-1]:+.4f} m "
              f"(MOM6 reference is roughly -0.23 m at 30+ days; see "
              f"README.md's MOM6-parity caveats -- 10 days here is a "
              f"quick smoke run, not the validation horizon)")


if __name__ == "__main__":
    main()
