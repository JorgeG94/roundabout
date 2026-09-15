"""Idealized ACC -- the all-features ocean showcase, Python-API port.

Adapted from ``acc_channel.nml``'s own header comment (the physics is
identical; only the front door changed):

A RE-ENTRANT (periodic-x) SPHERICAL channel in the Southern Ocean (~45 S)
with a seamount ridge mid-channel, exercising in ONE run:

  * spherical curvilinear metrics + PLANETARY Coriolis
    (f = 2*Omega*sin(lat), southern hemisphere => f < 0; the closures
    take |f|),
  * periodic-x OBC (re-entrant channel) + walls north/south,
  * EPBL surface boundary layer (RH18 energetics-based PBL),
  * kappa-shear (JHL08) interior mixing,
  * the production envelope: linear EOS (alpha_T=0.2, matching the
    epbl_mld thermal pattern), Smagorinsky KH+AH, zstar vcoord,
    distributed linear bottom drag, the diag manager.

Forcing: steady zonal (westerly) wind stress 0.12 Pa + weak surface
cooling (-40 W/m^2) so EPBL has BOTH wind and convective TKE.  The wind
is depth-uniform "constant" (there is no kernel for a meridional
tau_x(y) ramp on the ocean path beyond the closed-basin 2gyre cosine --
see ``README.md``'s caveats section; do not invent one here).

GRID SIZE -- READ THIS FIRST.  ``acc_channel.nml``'s own header comment
says "24 deg lon x 14 deg lat ... 0.5 deg resolution => 48 x 28 cells
... nz=16 over H=3000 m", but the file's actual ``&grid_nml``/
``&nonhydrostatic_nml`` values are ``nx=480, ny=280, nz_layers=75`` --
TEN TIMES the described lon/lat extent and ~4.7x the described layer
count (longitude 0-240 deg, latitude -52 to +88 deg). This is a
pre-existing mismatch in the checked-in ``.nml`` between its own prose
and its own numbers (not something introduced by this port -- the
``.nml`` is treated read-only per the porting brief). ``README.md``
separately claims this size "OOMs on a single 32 GB V100 at the first
kernel launch" (~14 GB mapped state + workspaces past 32 GB) -- that
was NOT reproduced while porting this script: a 30-outer-step run at
the literal 480x280x75 size peaked at ~22 GB / 32 GB device memory on
a V100 (see the port report's byte-identity gate, which ran BOTH this
script's emitted config and the stock .nml at this exact size to
completion). The OOM may be real for a much longer run (more lazy
workspace growth over more steps) or a different build/GPU than the
one used here; ``--full`` is still opt-in rather than default, both
because 480x280x75 is an unusually large "example" to run by default
and because the discrepancy with the .nml's own header prose (above)
makes 48x28x16 the more defensible default for a quick local smoke
run. ``CANONICAL_GRID`` below reproduces the ``.nml``'s ACTUAL numbers
byte-for-byte for the configuration-identity gate regardless of which
default ``main()`` picks.

Run:
    cd validation_examples/ocean/acc_channel
    # first load your NVHPC + NetCDF toolchain (module load / spack / conda;
    # see environments/) -- one toolchain per shell, never stack two
    CUDA_VISIBLE_DEVICES=<gpu> python3 acc_channel.py

See ``acc_channel_quiescent.py`` for the resting-state tier-1.5 trap
(same grid/bathy/BCs, zero forcing -- must stay at rest).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make `rdb` importable when this script is run directly (`python3
# acc_channel.py`) from any cwd, without requiring the caller to set
# PYTHONPATH by hand -- mirrors python/tests/conftest.py's own
# sys.path bootstrap for the same "no pip install" package.
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

#: The .nml's ACTUAL &grid_nml/&nonhydrostatic_nml numbers (not the
#: header prose's 48x28x16 -- see module docstring).
CANONICAL_GRID = dict(nx=480, ny=280, nz=75)
#: A grid small enough to actually run on one GPU -- the size the
#: .nml's OWN header text describes, not a value this port invented.
SMOKE_GRID = dict(nx=48, ny=28, nz=16)


def build_model(*, nx, ny, nz):
    """Compose the curated `rdb.Model` for the ACC channel at grid
    size (nx, ny, nz). Everything else (bathymetry ridge, forcing,
    closures, BCs) is fixed -- only the grid resolution is parametrised
    so this can run at SMOKE_GRID locally and CANONICAL_GRID for the
    configuration-identity gate.
    """
    grid = LatitudeLongitudeGrid(
        size=(nx, ny, nz),
        longitude=(0.0, 0.0 + nx * 0.5),   # dx = 0.5 deg/cell (&grid_nml)
        latitude=(-52.0, -52.0 + ny * 0.5),  # dy = 0.5 deg/cell
        radius=6371000.0,                  # &ocean_grid_nml rad_earth
        # Re-entrant channel: periodic in x (west/east), walls north/south.
        # NOTE this is a SECTOR re-entrant channel, not a full 360-degree
        # wrap -- Topology.PERIODIC is what forces the ghost-wrap seam
        # regardless of the actual longitude extent (see D2.1).
        topology=(Topology.PERIODIC, Topology.BOUNDED),
        halo=None,   # derived via required_halo() -- must land on 3
                     # (periodic-x needs nghost>=3 for the PPM +
                     # biharmonic seam stencil); matches &grid_nml
                     # nghost=3 explicitly.
    )

    model = rdb.Model(
        grid=grid,
        timestep=300.0,                    # &time_nml dt_fixed

        coriolis=SphericalCoriolis(rotation_rate=7.2921e-5),
        momentum_advection=Sadourny(form="enstrophy"),  # -> form="sadourny"
        buoyancy=LinearEquationOfState(
            thermal_expansion=0.2, reference_density=1035.0),
        pressure=FiniteVolumeLite(),

        vcoord=ZStar(),

        closures=[
            EPBL(mstar_scheme="om4", mld_use_prev_guess=True),
            # EPBL() above ALREADY sets &ocean_vmix_nml use_kpp = .false.
            # itself (KPP defaults ON in the schema; EPBL/KPP are
            # mutually exclusive) -- no separate off-switch needed.
            KappaShear(),                  # ri_crit/prandtl at their own
                                            # (Fortran-matching) defaults
            Smagorinsky(C=0.15, biharmonic=True, C_biharmonic=0.06,
                        bound_kh=False, nu_h=10000.0, kh_vel_scale=3.0e-3),
            LinearDrag(r=2.5e-5, hbbl=10.0, bg_vel=0.1,
                       bbl_thick_min=0.1),
            # MOM6 KV_ML_INVZ2 near-surface momentum viscosity band --
            # composes alongside EPBL, does not select a BL scheme.
            NearSurfaceViscosity(kv=0.01, hmix=20.0),
            # MOM6 HARMONIC_VISC vdiff-assembly parity.
            HarmonicFaceThickness(),
        ],

        forcing=[
            ConstantWind(tau_x=0.12, tau_y=0.0),
            SurfaceHeatFlux(q=-40.0),
        ],

        initial_condition=[
            UniformInitialCondition(temperature=8.0, salinity=35.0),
            # Stratified linear T(z): 2 degC bed (k=1) -> 8 degC surface
            # (k=nz); the UNIFORM initial_temperature above is only a
            # fallback the Fortran does not read once this is set.
            StratifiedInitialCondition(T_surface=8.0, T_bottom=2.0),
        ],

        barotropic=BarotropicSolver(auto=True, cfl_safety=0.65, bebt=0.2),

        diagnostics=Diagnostics(
            enabled=True, filename="acc_channel", every=1.0,
            output_dir="./out_acc_channel", status_interval=1.0),

        # &time_nml t_end/time_unit: Model.run(until=...) drives the
        # outer loop directly in raw seconds and never reads either
        # knob, but t_end still gates create() (see Duration's own
        # docstring) -- pin it explicitly rather than relying on the
        # Fortran's default t_end=1.0 SECOND, which the Diagnostics
        # every=1.0-DAY cadence above would otherwise exceed.
        duration=Duration(t_end=1.0, time_unit="day"),
    )

    # ------------------------------------------------------------------
    # Escape hatch: knobs with no curated writer at all, or where a
    # curated object's own constructor default would silently pick a
    # DIFFERENT value than the .nml sets.
    # ------------------------------------------------------------------
    cfg = model.config

    # &physics_nml coriolis_f: unused by the metric f-fill under
    # coriolis_scheme="planetary" (SphericalCoriolis fills f from
    # geolat), but the .nml keeps it non-zero "for any legacy guards".
    # DELIBERATELY left on the escape hatch, not given a curated route:
    # FPlane/BetaPlane already write &physics_nml.coriolis_f correctly
    # (this is NOT the same gap) -- SphericalCoriolis genuinely has
    # nothing to do with this value, and inventing a
    # SphericalCoriolis(legacy_coriolis_f=...) parameter would encode a
    # provably-inert number as if it meant something (the same judgment
    # call `ConstantWind` makes for `taux_magnitude` in seamount.py).
    cfg.physics.coriolis_f = 1.0e-4

    # Seamount ridge: Seamount() IS a curated writer, but slope_scale is
    # documented as METRES converted to grid units by the Fortran
    # (topo_length_to_grid_units) -- the .nml passes the raw number 3.0
    # unchanged (see module docstring: this looks like it was intended
    # as "~3 deg" in prose but is actually 3.0 METRES at the Fortran
    # dispatch). Reproduced verbatim, not "fixed", per the porting brief.
    Seamount(max_depth=3000.0, edge_depth=1500.0, slope_scale=3.0).apply(cfg)

    # FiniteVolumeLite() carries no maxvel knob (that lives on the
    # separate VelocityTruncation object -- a numerical safety net, not
    # a PGF variant, D2.6):
    VelocityTruncation(maxvel=2.0).apply(cfg)

    # &output_nml output_to_file: cosmetic-only, see
    # acc_channel_quiescent.py's matching comment -- a dead knob (P4
    # dead-knob sweep) kept for namelist-text fidelity with the .nml.
    cfg.output.output_to_file = False

    return model


def _surface_mean(field):
    """Mean over the surface layer (k = nz-1, the last index -- k=nz is
    the surface in the bottom-up layer convention) of a 3D Field, via
    ``.copy()`` (a plain nested Python list, no numpy dependency --
    ``Field.__getitem__`` only supports numpy-backed PARTIAL slicing,
    see D6.6 / the port report's numpy-availability note)."""
    data = field.copy()   # nested [i][j][k], full interior shape
    surface = [[col[-1] for col in row] for row in data]
    flat = [v for row in surface for v in row]
    return sum(flat) / len(flat)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full", action="store_true",
        help="use the .nml's literal grid size (480x280x75) -- matches "
             "the .nml exactly; off by default uses the SMOKE_GRID size "
             "the .nml's own header prose describes (48x28x16) so this "
             "runs fast even on a modest GPU")
    parser.add_argument("--days", type=float, default=1.0)
    args = parser.parse_args()

    grid_size = CANONICAL_GRID if args.full else SMOKE_GRID
    model = build_model(**grid_size)

    print(f"acc_channel: grid={grid_size}, "
          f"{'CANONICAL (matches .nml)' if args.full else 'SMOKE (reduced for local run)'}")
    print(model.config.to_namelist()[:400], "...\n")

    with model:
        u_series = []

        def daily(m):
            mean_u = _surface_mean(m.u)
            u_series.append(mean_u)
            print(f"day {m.time / 86400.0:6.2f}  mean surface u "
                  f"{mean_u:+.4f} m/s  KE {m.kinetic_energy:.4e} J")

        model.run(until=args.days * 86400.0, callback=daily, every=86400.0)

        print("\nFinal surface-u snapshot (zonal jet indicator):",
              f"mean={u_series[-1]:+.4f} m/s" if u_series else "n/a")


if __name__ == "__main__":
    main()
