"""P5/P6 acceptance test -- a shrunk port of
``tmp_local_artifacts/python_ffi_scope/sketch/01_regional_run.py``, the
document's own stated bar for "not awful": *"if a user cannot write
approximately this, the API has failed."*

Grid shrunk from 600x500x50 to 10x8x3 and the run from 30 days to 5
minutes so this runs in well under a second; everything else keeps the
SAME shape as the sketch: a LatitudeLongitudeGrid, GEBCO-convention array
bathymetry, Wright-1997 EOS, FV-MOM6 PGF, a composed closure list, wind +
heat + shortwave forcing, an ideal-age tracer, ICs written through the
fields, a run() with a callback, and post-run inspection.

Three DELIBERATE deviations from the sketch, each with a reason (see the
inline comments at the point of deviation, and the phase report for the
full account):

  1. ``boundaries=`` is dropped from the live model. Composing a
     `Boundaries(...)` object into `to_namelist()` text works correctly
     (asserted separately, `test_boundaries_composes_correct_text`
     below) but ACTUALLY CREATING a model with any `&ocean_bc_nml`
     content fails on the NVHPC/GPU build specifically -- a pre-existing
     Fortran bug in `read_ocean_bc_nml`'s in-memory (``lines(:)``)
     namelist read, independent of this phase (confirmed: the identical
     content creates successfully through the FILE path on gfortran).
  2. ``ImplicitVerticalFriction(drag=True)`` is dropped to
     ``drag=False``: combined with ``QuadraticDrag(hbbl=10.0,
     implicit=True)`` as the sketch has it, the Fortran's own
     `validate_config` rejects the combination ("ocean_vdiff
     implicit_drag is mutually exclusive with ocean_bdrag implicit").
     This is the sketch's OWN closure list conflicting with the
     Fortran's cross-knob rule -- not a P5 bug.
  3. The escape hatch gains two lines the sketch doesn't have:
     ``ocean_bt.auto_n_inner = True`` and ``time.t_end`` bumped up --
     neither is reachable from any curated constructor argument
     (`auto_n_inner` needs `barotropic=BarotropicSolver(auto=True)`,
     which this test also demonstrates; `t_end` has no curated argument
     at all), and both default to values that make `create()` refuse a
     namelist-driven run outright.

The dropped `PassiveTracer("dye", ...)` call is exercised separately
in `test_passive_tracer_unsupported` -- see its docstring.
"""

from __future__ import annotations

import numpy as np
import pytest

import rdb
from rdb.barotropic import BarotropicSolver
from rdb.boundaries import Boundaries, Flather, Radiation, Sponge, Wall
from rdb.buoyancy import Wright1997
from rdb.closures import (
    KPP,
    BryanLewisBackground,
    ConvectiveAdjustment,
    HorizontalDiffusivity,
    ImplicitVerticalFriction,
    KappaShear,
    PacanowskiPhilander,
    QuadraticDrag,
    Smagorinsky,
)
from rdb.coriolis import Sadourny, SphericalCoriolis, WENO
from rdb.forcing import ConstantWind, ShortwavePenetration, SurfaceHeatFlux
from rdb.grids import Bathymetry, LatitudeLongitudeGrid, Topology
from rdb.pressure import FiniteVolumeMOM6
from rdb.tracers import IdealAge, PassiveTracer
from rdb.vcoord import ZStarFull


def _build_grid():
    return LatitudeLongitudeGrid(
        size=(10, 8, 3),
        longitude=(140.0, 141.0),
        latitude=(-10.0, -9.2),
        radius=6.378e6,
        topology=(Topology.BOUNDED, Topology.BOUNDED),
        halo=None,
    )


def _build_bathymetry():
    depth = -500.0 * np.ones((10, 8))  # GEBCO-style: negative height
    bathymetry = Bathymetry(
        depth,
        convention="height_positive_up",  # required keyword, no default
        land_threshold=2.0,
        minimum_depth=10.0,
    )
    assert bathymetry.wet_fraction == 1.0
    assert "10x8" in repr(bathymetry)
    return bathymetry


def test_worked_example_runs():
    grid = _build_grid()
    bathymetry = _build_bathymetry()

    model = rdb.Model(
        grid=grid,
        bathymetry=bathymetry,
        timestep=60.0,
        coriolis=SphericalCoriolis(),
        momentum_advection=Sadourny(form="enstrophy",
                                     pv_advection=WENO(order=5)),
        buoyancy=Wright1997(),
        pressure=FiniteVolumeMOM6(gfs_scale=0.98, mass_weight=True),
        vcoord=ZStarFull(h_surf_target=2.0, h_min=1.0e-3, remap="ppm"),
        closures=[
            PacanowskiPhilander(nu0=1e-2, nu_bg=1e-4, kappa_bg=1e-5),
            KPP(ri_crit=0.3, shortwave_method="mxl_sw"),
            KappaShear(ri_crit=0.25, at_vertex=True),
            ConvectiveAdjustment(kd=1.0),
            BryanLewisBackground(kd_surface=1e-5, kd_deep=1.3e-4),
            # deviation 2: drag=False (drag=True conflicts with
            # QuadraticDrag(implicit=True) below per the Fortran's own
            # validate_config -- see module docstring).
            ImplicitVerticalFriction(stress=True, drag=False),
            Smagorinsky(C=0.15, biharmonic=True, C_biharmonic=0.06,
                        form="stress_tensor", bound_kh=True),
            HorizontalDiffusivity(kappa_h=25.0),
            QuadraticDrag(cd=2.5e-3, hbbl=10.0, bg_vel=0.1, implicit=True),
        ],
        # deviation 1: boundaries= dropped -- see module docstring;
        # exercised compose-only in test_boundaries_composes_correct_text.
        forcing=[
            ConstantWind(tau_x=0.05, tau_y=0.0),
            SurfaceHeatFlux(q=-20.0),
            ShortwavePenetration(fraction=0.45, band_ratio=0.58,
                                  zeta1=0.35, zeta2=23.0),
        ],
        tracers=[IdealAge()],
        # deviation 3a: needed to get a non-zero barotropic substep count
        # at all -- the sketch never sets barotropic=.
        barotropic=BarotropicSolver(auto=True),
    )

    # ---- the escape hatch: same knobs the curated objects above write,
    # reachable and composing in the SAME script (D2.0 rule 1).
    model.config.ocean_epbl.mstar = 1.2               # inert (EPBL off)
    model.config.ocean_debug.chksum = True
    model.config.ocean_debug.chksum_start_step = 100
    model.config.ocean_bt.split_scheme = "pred_corr"
    model.config.ocean_bt.pc_be = 0.6
    # deviation 3b: no curated t_end=/duration= argument exists at all.
    model.config.time.t_end = 1.0e6
    # No curated Diagnostics()/Restart() object exists yet -- the Fortran
    # defaults leave diag/output on with a relative "./output/..." path
    # whose failure is an uncaught ERROR STOP (P0 review F5), not a typed
    # exception. Disable both so this test never depends on the cwd
    # having a writable ./output/ directory.
    model.config.ocean_diag.enabled = False
    model.config.output.output_to_file = False

    text = model.config.to_namelist()
    assert "&ocean_hvisc_nml" in text
    assert "mstar = 1.2" in text
    assert "split_scheme = \"pred_corr\"" in text

    with model:
        model.temperature[:] = np.full(model.temperature.shape, 18.0)
        model.salinity[:] = np.full(model.salinity.shape, 34.7)

        sst_series = []

        def daily(m):
            sst = m.temperature[:, :, -1]         # surface layer
            sst_series.append(float(np.nanmean(sst)))

        model.run(until=5 * 60.0, callback=daily, every=60.0)

        assert model.step_count == 5
        assert len(sst_series) == 5
        assert all(np.isfinite(s) for s in sst_series)

        u = np.asarray(model.u)
        h = np.asarray(model.h)
        assert u.shape == (11, 8, 3)
        assert h.shape == (10, 8, 3)

        eta = model.ssh
        assert len(eta) == 10 and len(eta[0]) == 8

        snapshot = model.temperature.copy()
        assert len(snapshot) == 10


def test_z_kwarg_raises_on_model_and_grid():
    """z= is accepted ONLY to raise a specific, actionable error -- never
    silently ignored via **kwargs (D6.4)."""
    with pytest.raises(rdb.RdbUnsupportedError, match="vcoord="):
        rdb.Model(z=[0, -10, -20])
    with pytest.raises(rdb.RdbUnsupportedError, match="vcoord="):
        LatitudeLongitudeGrid(size=(4, 4, 2), longitude=(0, 1),
                               latitude=(0, 1), z=[0, -10])


def test_flat_topology_raises():
    with pytest.raises(rdb.RdbUnsupportedError, match="Flat"):
        LatitudeLongitudeGrid(size=(4, 4, 2), longitude=(0, 1),
                               latitude=(0, 1), topology=("flat", "flat"))


def test_bathymetry_convention_is_required():
    with pytest.raises(TypeError):
        Bathymetry(-500.0 * np.ones((4, 4)))  # convention= missing


def test_bathymetry_wrong_convention_raises_sign_error():
    gebco = -500.0 * np.ones((4, 4))  # negative-down (height convention)
    with pytest.raises(rdb.BathymetrySignError, match="wet cells"):
        Bathymetry(gebco, convention="depth_positive_down")


def test_boundaries_composes_correct_text():
    """`Boundaries(...)` with the SAME shape as the sketch composes the
    correct namelist text (D2.0 rule 1) even though actually CREATING a
    model with it is blocked by the NVHPC-build ocean_bc bug (module
    docstring). Verifies the writer half of the contract independently
    of that Fortran gap."""
    from rdb._config import Config

    boundaries = Boundaries(
        west=Wall(),
        east=Flather(radiation=Radiation("orlanski", rx_max=0.7)),
        south=Flather(radiation=Radiation("orlanski", rx_max=0.7)),
        north=Sponge(width=2, strength=1.0 / 86400.0),
    )
    cfg = Config()
    boundaries.apply(cfg, periodic_x=False, periodic_y=False)
    text = cfg.to_namelist()
    assert '&ocean_bc_nml' in text
    assert 'west = "wall"' in text
    assert 'east = "open"' in text
    assert 'south = "open"' in text
    assert 'north = "sponge"' in text
    assert 'radiation_scheme = "orlanski"' in text
    assert 'sponge_width = 2' in text


def test_passive_tracer_unsupported():
    """The sketch's ``PassiveTracer("dye", units="1", long_name="passive
    dye")`` call has no reachable Fortran entry point in this phase (no
    Fortran changes) -- see rdb.tracers.PassiveTracer's docstring."""
    with pytest.raises(rdb.RdbUnsupportedError, match="register_passive_tracer"):
        PassiveTracer("dye", units="1", long_name="passive dye")
