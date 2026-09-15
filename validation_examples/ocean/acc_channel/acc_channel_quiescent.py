"""Quiescent ACC channel -- the SPHERICAL tier-1.5 trap, Python-API port.

Adapted from ``acc_channel_quiescent.nml``'s own header comment:

SAME grid / seamount bathymetry / periodic-x BCs / closures (EPBL +
kappa-shear ON) as ``acc_channel.py``, but ZERO wind, ZERO surface heat
flux, and UNIFORM T,S (no available potential energy). This is the
honest tier-1.5 trap for a sigma/zstar dyn-core: with no APE the resting
state MUST stay identically at rest -- max|u|,|v| < ~1e-10 m/s over 5
days. Any drift here is a real finding (metric/PGF inconsistency on the
spherical grid at the ridge, or a periodic-seam metric issue), NOT a
tuning target.

VERIFIED (2026-06-11, V100, Fortran/.nml path): max|u|,|v| = 1.0e-11 /
5.5e-12 m/s, En = 6e-25 J, MaxCFL = 0 -- passes at roundoff. The
spherical metrics + planetary Coriolis + periodic-x seam + seamount
bathymetry are mutually consistent.

NOTE: a STRATIFIED resting state over the seamount (the acc_channel.py
T(z) IC with zero forcing) does NOT stay at rest -- it drifts to
~0.16 m/s over the ridge crest. That is the well-known sigma/zstar PGF
residual over a slope under stratification (APE-over-slope), not a
metric/seam bug -- see ``acc_channel.py``'s own module docstring and
``README.md``. The seamount setter's own docstring specifies "uniform
T,S, no APE" for the zero-motion expectation; this script honours that
by using a UNIFORM temperature profile (T_init_surface == T_init_bottom).

Unlike ``acc_channel.py`` this grid (48x28x16) is small enough to run
as-is -- there is no SMOKE_GRID/CANONICAL_GRID split here.

CURATED-API GAP FROM AN EARLIER PORT, NOW CLOSED: an earlier version of
this script found ``LatitudeLongitudeGrid(topology=(Topology.PERIODIC,
...))`` alone did NOT make the live model periodic-x, and worked around
it with an explicit ``cfg.ocean_bc.west = cfg.ocean_bc.east =
"periodic"`` escape hatch. That gap was fixed in
``rdb._compose.compose`` (commit ``36abbfbf``, "grid topology must
reach &ocean_bc_nml without a Boundaries object") -- topology=(PERIODIC,
...) alone is now sufficient; the escape hatch below has been removed
(verified directly: the solver's own startup banner now prints ``OBC:
west=periodic east=periodic ...`` with no manual override).

Run:
    cd validation_examples/ocean/acc_channel
    # first load your NVHPC + NetCDF toolchain (module load / spack / conda;
    # see environments/) -- one toolchain per shell, never stack two
    CUDA_VISIBLE_DEVICES=<gpu> python3 acc_channel_quiescent.py
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
from rdb.buoyancy import LinearEquationOfState
from rdb.closures import (
    EPBL,
    HarmonicFaceThickness,
    KappaShear,
    LinearDrag,
    NearSurfaceViscosity,
    Smagorinsky,
)
from rdb.coriolis import Sadourny, SphericalCoriolis
from rdb.diagnostics import Diagnostics, Duration
from rdb.forcing import (
    ConstantWind,
    StratifiedInitialCondition,
    SurfaceHeatFlux,
    UniformInitialCondition,
)
from rdb.grids import LatitudeLongitudeGrid, Seamount, Topology
from rdb.pressure import FiniteVolumeLite, VelocityTruncation
from rdb.vcoord import ZStar

GRID = dict(nx=48, ny=28, nz=16)  # matches the .nml's &grid_nml exactly


def build_model():
    grid = LatitudeLongitudeGrid(
        size=(GRID["nx"], GRID["ny"], GRID["nz"]),
        longitude=(0.0, 0.0 + GRID["nx"] * 0.5),
        latitude=(-52.0, -52.0 + GRID["ny"] * 0.5),
        radius=6371000.0,
        topology=(Topology.PERIODIC, Topology.BOUNDED),
        halo=None,  # derives to 3 (periodic-x PPM/biharmonic seam)
    )

    model = rdb.Model(
        grid=grid,
        timestep=1200.0,   # &time_nml dt_fixed

        coriolis=SphericalCoriolis(rotation_rate=7.2921e-5),
        momentum_advection=Sadourny(form="enstrophy"),
        buoyancy=LinearEquationOfState(
            thermal_expansion=0.2, reference_density=1035.0),
        pressure=FiniteVolumeLite(),

        vcoord=ZStar(),

        closures=[
            EPBL(mstar_scheme="om4", mld_use_prev_guess=True),
            # EPBL() above already sets use_kpp = .false. itself.
            KappaShear(),
            Smagorinsky(C=0.15, biharmonic=True, C_biharmonic=0.06,
                        bound_kh=False, nu_h=10000.0, kh_vel_scale=3.0e-3),
            LinearDrag(r=2.5e-5, hbbl=10.0, bg_vel=0.1,
                       bbl_thick_min=0.1),
            NearSurfaceViscosity(kv=0.01, hmix=20.0),
            HarmonicFaceThickness(),
        ],

        # ZERO forcing -- the whole point of the trap.
        forcing=[
            ConstantWind(tau_x=0.0, tau_y=0.0),
            SurfaceHeatFlux(q=0.0),
        ],

        # UNIFORM T,S -- no APE. initial_salinity=35.0 already matches
        # the Fortran default so it is a no-op text-wise either way.
        # T_surface == T_bottom (both 8.0) routes through the SAME
        # stratified code path as acc_channel.py but with zero gradient
        # -- i.e. genuinely uniform T (see module docstring).
        initial_condition=[
            UniformInitialCondition(temperature=8.0, salinity=35.0),
            StratifiedInitialCondition(T_surface=8.0, T_bottom=8.0),
        ],

        barotropic=BarotropicSolver(auto=True, cfl_safety=0.65, bebt=0.2),

        diagnostics=Diagnostics(
            enabled=True, filename="acc_channel_quiescent", every=1.0,
            output_dir="./out_acc_channel_quiescent", status_interval=1.0),

        duration=Duration(t_end=5.0, time_unit="day"),
    )

    # ---- escape hatch (same reasons as acc_channel.py; see its module
    # docstring / comments for the fuller explanation) ----
    cfg = model.config

    # &physics_nml coriolis_f: SphericalCoriolis fills f from geolat, so
    # this legacy value is provably inert under coriolis_scheme=
    # "planetary" -- see acc_channel.py's matching comment for why this
    # deliberately has no curated route (NOT the same gap as FPlane,
    # which already writes this knob correctly).
    cfg.physics.coriolis_f = 1.0e-4

    Seamount(max_depth=3000.0, edge_depth=1500.0, slope_scale=3.0).apply(cfg)
    VelocityTruncation(maxvel=2.0).apply(cfg)

    # &output_nml output_to_file: Diagnostics(enabled=True) above sets
    # this True (it also controls &ocean_diag_nml enabled, its real
    # gate); the .nml independently sets it False. Purely COSMETIC --
    # the P4 dead-knob sweep found output_to_file is validated but read
    # NOWHERE on the ocean path (RdbDeadKnobWarning fires on this
    # line) -- kept only for namelist-text fidelity with the .nml.
    cfg.output.output_to_file = False

    return model


def _max_abs(field):
    """Max |value| over a whole 3D Field, via ``.copy()`` (a plain
    nested Python list -- no numpy dependency; ``Field.__getitem__``
    only supports numpy-backed PARTIAL slicing, see D6.6 / the port
    report's numpy-availability note)."""
    data = field.copy()
    return max(abs(v) for row in data for col in row for v in col)


def main():
    model = build_model()
    print("acc_channel_quiescent: grid=", GRID)

    with model:
        max_speed = 0.0

        def daily(m):
            nonlocal max_speed
            day_max = max(_max_abs(m.u), _max_abs(m.v))
            max_speed = max(max_speed, day_max)
            print(f"day {m.time / 86400.0:5.2f}  max|u,v| this frame "
                  f"{day_max:.3e} m/s  KE {m.kinetic_energy:.3e} J")

        model.run(until=5 * 86400.0, callback=daily, every=86400.0)

        print(f"\nmax|u|,|v| over the run: {max_speed:.3e} m/s "
              f"(expect ~1e-10 or smaller -- a resting ocean must stay "
              f"at rest)")
        assert max_speed < 1.0e-6, (
            f"quiescent trap FAILED: max speed {max_speed:.3e} m/s is "
            f"well above roundoff")


if __name__ == "__main__":
    main()
