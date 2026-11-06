"""Isolated balanced vortex on a beta-plane -- Python-API port of the
three namelists in this directory.

    python3 vortex_beta.py             # vortex_beta.nml (warm core, anticyclone)
    python3 vortex_beta.py cyclone     # vortex_beta_cyclone.nml (cold core)
    python3 vortex_beta.py fplane      # vortex_fplane.nml (beta = 0 control)

A single mesoscale eddy starts in geostrophic / thermal-wind balance and
drifts on beta: anticyclones SW, cyclones NW (McWilliams & Flierl 1979).
Physics, construction and the measured 40-day drift + tilt: README.md.
Analyse the output with analyze_vortex_tilt.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_PY_DIR = str(_REPO_ROOT / "python")
if _PY_DIR not in sys.path:
    sys.path.insert(0, _PY_DIR)

import rdb
from rdb.barotropic import BarotropicSolver
from rdb.closures import ConstantViscosity
from rdb.coriolis import BetaPlane
from rdb.diagnostics import Diagnostics, Duration
from rdb.forcing import Thermodynamics
from rdb.grids import FlatBottom, RectilinearGrid, Topology
from rdb.vcoord import Sigma

DAY = 86400.0

CASES = {
    #  name       (filename,              vortex_dT, beta)
    "warm":    ("vortex_beta",          2.0,  2.0e-11),
    "cyclone": ("vortex_beta_cyclone", -2.0,  2.0e-11),
    "fplane":  ("vortex_fplane",        2.0,  0.0),
}


def build_model(case="warm"):
    filename, vortex_dT, beta = CASES[case]

    grid = RectilinearGrid(
        size=(80, 80, 40),                 # 400 km x 400 km at 5 km, 40 x 100 m layers
        extent=(400.0e3, 400.0e3),
        topology=(Topology.BOUNDED, Topology.BOUNDED),
        halo=4,
    )

    model = rdb.Model(
        grid=grid,
        timestep=600.0,                    # internal-wave CFL: c1 = N H / pi = 3.2 m/s
        coriolis=BetaPlane(f0=1.0e-4, beta=beta, y_ref=200.0e3),
        vcoord=Sigma(),                    # flat bed: sigma = z levels
        closures=[ConstantViscosity(nu_h=100.0)],
        forcing=[Thermodynamics(enabled=True)],
        barotropic=BarotropicSolver(auto=True, scheme="pred_corr"),
        diagnostics=Diagnostics(
            ["vorticity_z:1d", "SSH:1d", "temperature:1d",
             "u:off", "v:off", "KE:off"],
            every=DAY, filename=filename,
            output_dir=f"./out_{filename}", status_interval=1.0),
        duration=Duration(t_end=40.0, time_unit="day"),
    )

    cfg = model.config
    FlatBottom(depth=4000.0).apply(cfg)

    # Linear EOS (alpha_T is dimensional: 1.7e-4 1/K x rho_0), no vertical mixing
    cfg.ocean_eos.eos = "linear"
    cfg.ocean_ic.alpha_T = 0.17
    cfg.ocean_ic.rho_0 = 1035.0
    cfg.ocean_vmix.use_closure = False
    cfg.ocean_vmix.use_kpp = False
    cfg.tracer.initial_salinity = 35.0

    # The eddy: Gaussian T anomaly, R = 50 km, Hv = 700 m, N = 2.5e-3 1/s
    ic = cfg.ocean_ic
    ic.ic_config = "balanced_vortex"
    ic.vortex_radius = 50.0e3
    ic.vortex_hv = 700.0
    ic.vortex_dT = vortex_dT               # + warm core / anticyclone, - cold core / cyclone
    ic.vortex_T_surf = 20.0
    ic.vortex_dT_dz = 4.0e-3               # a cold core needs > 0.858 |dT| / Hv
    ic.vortex_balance = "geostrophic"
    return model


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "warm"
    with build_model(case) as model:
        def daily(m):
            print(f"day {m.time / DAY:5.1f}  KE {m.kinetic_energy:.4e} J")
        model.run(until=40 * DAY, callback=daily, every=DAY)


if __name__ == "__main__":
    main()
