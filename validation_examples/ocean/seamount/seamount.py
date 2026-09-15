"""Idealised seamount -- Tier-1.5 sigma-coord bathymetry bridge test,
Python-API port.

Adapted from ``seamount.nml``'s own header comment:

Goal: validate the ocean dyn-core operators (BPG, Coriolis-adv,
continuity, vmix) over a varying-h_layer bathymetry WITHOUT the
confounds that a real-bathymetry regional run carries (50x depth ratio,
IC stratification feeding APE, closed-basin baroclinic resonance,
sponge edge effects).

Setup:
  * 100 km x 100 km closed basin (50 x 50 cells @ 2 km).
  * Single Gaussian seamount centred at (50 km, 50 km): peak depth
    200 m, basin depth 4000 m, half-width 25 km. Maximum slope
    ~150 m/km -- comparable to a real shelf break.
  * 15 sigma-layers.
  * UNIFORM T = 15 degC, S = 35 PSU -- no stratification, no APE.
    Density is horizontally and vertically uniform; the BPG over
    uniform rho + sloping bathy should cancel exactly.
  * f-plane Coriolis at f = 0 -- no Coriolis activity from residual
    flow. Eliminates baroclinic-instability resonance entirely.
  * No wind, no surface flux, no sponges. Closed walls.

Expected outcome: IDENTICALLY ZERO MOTION for 30 days. En = 0,
max|u| = 0, SSH = 0 +/- machine roundoff. Any non-zero u, v, eta is a
sigma-coord numerical artefact in the BPG, Coriolis-adv, or continuity
kernels -- this is the cleanest possible test of those operators'
behaviour under varying h_layer. (This is the run picked as "the most
canonical resting-state seamount case" among the variants in this
directory -- the README explicitly labels ``seamount.nml`` "Canonical
setup".)

Run:
    cd validation_examples/ocean/seamount
    # first load your NVHPC + NetCDF toolchain (module load / spack / conda;
    # see environments/) -- one toolchain per shell, never stack two
    CUDA_VISIBLE_DEVICES=<gpu> python3 seamount.py
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
from rdb.closures import KPP, ConstantViscosity, PacanowskiPhilander
from rdb.coriolis import FPlane, Sadourny
from rdb.diagnostics import Diagnostics, Duration
from rdb.forcing import ConstantWind, UniformInitialCondition
from rdb.grids import RectilinearGrid, Seamount, Topology
from rdb.vcoord import ZStarSigma

GRID = dict(nx=50, ny=50, nz=15)


def build_model():
    grid = RectilinearGrid(
        size=(GRID["nx"], GRID["ny"], GRID["nz"]),
        extent=(100000.0, 100000.0),   # 100 km x 100 km, metres
        topology=(Topology.BOUNDED, Topology.BOUNDED),  # closed walls
        halo=2,   # matches the .nml's explicit &grid_nml nghost = 2
                  # (the derived minimum for this non-periodic, centered
                  # PV, PPM-tracer, no-fold, no-vertex-kappa-shear case
                  # would also land on 2 -- pinned explicitly here for
                  # 1:1 fidelity with the .nml rather than relying on
                  # the derivation matching by construction).
    )

    model = rdb.Model(
        grid=grid,
        timestep=300.0,           # &time_nml dt_fixed, 5 min outer

        coriolis=FPlane(f=0.0),   # f-plane, f = 0 (no Coriolis feedback)
        momentum_advection=Sadourny(form="enstrophy"),  # -> "sadourny"

        # No buoyancy=/pressure= object: the .nml never touches
        # &ocean_eos_nml (default "linear") or &ocean_pgf_nml form
        # (default "mont").
        vcoord=ZStarSigma(),   # "Production default -- same as Tasman"
                                # per the .nml's own comment; NOT the
                                # plain "sigma" default.

        closures=[
            # "No viscosity -- pure dyn-core test": ConstantViscosity()
            # is the curated route to an explicit nu_h=0.0/nu_4=0.0
            # floor (lateral_closure="none" is already the Fortran
            # default, so this line is behaviourally a no-op, kept for
            # the same reason the .nml states it explicitly).
            ConstantViscosity(nu_h=0.0, nu_4=0.0),
            # "Pure dyn-core test": no interior closure, no boundary
            # layer at all. PacanowskiPhilander(enabled=False) ->
            # use_closure=.false.; KPP(enabled=False) -> use_kpp=
            # .false. (KPP defaults ON in the schema). Composing
            # KPP(enabled=False) with no EPBL(...) is EXACTLY the
            # deliberate "no boundary-layer scheme" case
            # NoBoundaryLayerWarning documents -- expected and fine
            # here (closed basin, zero wind, uniform rho).
            PacanowskiPhilander(enabled=False),
            KPP(enabled=False),
        ],

        # wind_config="constant" with zero stress -- ConstantWind(0, 0)
        # (the .nml ALSO sets ocean_topo.taux_magnitude=0.0, a 2gyre-only
        # knob that is inert under wind_config="constant"; matches its
        # own default 0.1->left unset would diverge textually so it is
        # NOT reproduced here -- see the port report's fallback list).
        forcing=[ConstantWind(tau_x=0.0, tau_y=0.0)],

        initial_condition=UniformInitialCondition(
            temperature=15.0, salinity=35.0),

        # This .nml pins a FIXED n_inner=60 and never touches
        # auto_n_inner, so auto=False (now also the curated default,
        # matching the Fortran) is stated explicitly for the record.
        barotropic=BarotropicSolver(n_inner=60, auto=False),

        diagnostics=Diagnostics(
            enabled=True, filename="seamount", every=21600.0,  # 6-hourly
            output_dir="./out_seamount", status_interval=86400.0,
            log_level="info"),  # log_level matches the Fortran default;
                                 # kept for 1:1 textual fidelity with the
                                 # .nml, which also states it explicitly.

        duration=Duration(t_end=2592000.0),   # 30 days, time_unit="s"
                                               # (the default -- the
                                               # .nml never touches it).
    )

    cfg = model.config

    Seamount(max_depth=4000.0, edge_depth=200.0,
             slope_scale=25000.0).apply(cfg)
    # taux_magnitude is a "2gyre"-wind_config-only knob, inert here
    # under wind_config="constant" -- but the .nml sets it (to 0.0,
    # non-default 0.1) explicitly anyway, and ConstantWind() has no
    # route to it (it only writes &physics_nml wind_stress_x/y):
    cfg.ocean_topo.taux_magnitude = 0.0

    cfg.output.output_to_file = False   # cosmetic-only dead knob, see
                                         # acc_channel.py's matching note

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
    print(model.config.to_namelist()[:400], "...\n")

    with model:
        max_speed = 0.0

        def daily6h(m):
            nonlocal max_speed
            frame_max = max(_max_abs(m.u), _max_abs(m.v))
            max_speed = max(max_speed, frame_max)
            print(f"t={m.time / 3600.0:7.1f} h  max|u,v| this frame "
                  f"{frame_max:.3e} m/s  KE {m.kinetic_energy:.3e} J")

        model.run(until=2 * 86400.0, callback=daily6h, every=21600.0)

        print(f"\nmax|u|,|v| over the run: {max_speed:.3e} m/s "
              f"(expect machine zero -- uniform rho + sloping bathymetry "
              f"should produce identically zero motion)")
        assert max_speed < 1.0e-6, (
            f"seamount quiescent trap FAILED: max speed {max_speed:.3e} "
            f"m/s is well above roundoff")


if __name__ == "__main__":
    main()
