#!/usr/bin/env python3
"""Run the global 1-degree case in-process through roundabout's Python API
and make its movie on the fly — standard library only.

    export RDB_DATA_DIR=/path/to/data          # holds OM_1deg/ (tools/fetch_om1deg.py)
    export RDB_LIB=/path/to/build/librdb_core.so
    export CUDA_VISIBLE_DEVICES=0              # GPU build: one device
    python3 validation_examples/ocean/global_1deg/run_global_1deg.py --days 365

What it shows:

1. The configuration is built with the typed ``rdb.Config`` API, knob by
   knob — the same knobs ``global_1deg_unforced.nml`` sets.  Before the run
   the script checks that against the namelist file (``--no-check`` skips
   it), so the two cannot drift.
2. ``rdb.Model(config)`` creates the solver; ``model.step(n)`` advances it.
3. Every simulated day the diagnostic manager's own in-memory buffers
   (``model.diagnostic("SSH")``, ``"u"``, ``"v"`` — the top-10 m daily means
   the namelist asks for) are read straight out of the running model, with
   no NetCDF round trip, and a movie frame is rendered from them
   (``global_movie.py``: stdlib PPM writer + ffmpeg).
4. A daily log line: kinetic energy per unit mass, total mass drift, and the
   extreme of the daily-mean surface current (``--max-velocity`` adds the
   3-D max |u|, |v| of the instantaneous state and where it sits).

The ``rdb`` executable reading ``global_1deg_unforced.nml`` is the reference
run (its console carries the mass / salt / heat budget residuals); this
script drives the SAME configuration through the C ABI, and its diagnostic
file is bit-identical to the executable's.
"""

import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
sys.path.insert(0, HERE)

import rdb  # noqa: E402
from global_movie import FrameRenderer, encode, model_tpoints, speed_from_uv  # noqa: E402

NAMELIST = os.path.join(HERE, "global_1deg_unforced.nml")
RHO0_API = 1025.0  # rdb_ocean_get_total_mass's reference density (RHO0_DIAG)


def build_config(data_dir, days, out_dir):
    """The global 1-degree unforced case, knob by knob (the Config API)."""
    cfg = rdb.Config()
    cfg.sim.sim_type = "ocean"

    # Grid: MOM6 OM_1deg mosaic supergrid, periodic in x, tripolar fold.
    cfg.grid.nx, cfg.grid.ny, cfg.grid.nghost = 360, 320, 3
    cfg.ocean_grid.grid_config = "supergrid"
    cfg.ocean_grid.supergrid_file = os.path.join(data_dir, "ocean_hgrid.nc")
    cfg.ocean_grid.coriolis_scheme = "planetary"
    cfg.ocean_grid.rad_earth = 6.371e6
    cfg.ocean_bc.west = cfg.ocean_bc.east = "periodic"
    cfg.ocean_bc.south = "wall"
    cfg.ocean_bc.north = "tripolar_fold"

    # Time: MOM6 OM_1deg DT; t_end bounds nothing here (we step ourselves)
    # but the diag cadence is validated against it.
    cfg.time.t_end = float(days)
    cfg.time.time_unit = "day"
    cfg.time.dt_fixed = 1800.0
    cfg.time.cfl_interval = 1

    # Unforced: no wind, no surface fluxes, no ice.
    cfg.physics.wind_stress_x = 0.0
    cfg.physics.wind_stress_y = 0.0
    cfg.ocean_topo.wind_config = "constant"
    cfg.ocean_topo.taux_magnitude = 0.0

    # Vertical: 50 stretched z levels (2 m at the top), partial-step faces closed.
    cfg.nonhydrostatic.nz_layers = 50
    cfg.vcoord.vcoord_type = "z_fixed"
    cfg.vcoord.zfixed_closed_faces = True
    cfg.vcoord.z_fixed_profile = "tanh"
    cfg.vcoord.z_fixed_dz_top = 2.0
    cfg.vcoord.z_fixed_tanh_center = 0.5
    cfg.vcoord.z_fixed_tanh_width = 0.25
    cfg.vcoord.check_vanished_content = True

    # Bathymetry (OM_1deg topog.nc, MINIMUM_DEPTH 9.5 m) + WOA13 January T/S.
    cfg.ocean_topo.topo_config = "file"
    cfg.ocean_topo.max_depth = 6500.0
    cfg.output.bathymetry_file = os.path.join(data_dir, "bathy_om1deg.nc")
    cfg.output.output_dir = out_dir
    cfg.output.output_to_file = False
    cfg.ocean_zinit.enable = True
    cfg.ocean_zinit.source = "file"
    cfg.ocean_zinit.file = os.path.join(data_dir, "ic_woa13_jan.nc")

    # Dynamics and closures (OM_1deg-like).
    cfg.ocean_pgf.form = "fv_mom6"
    cfg.ocean_pgf.maxvel = 6.0
    cfg.ocean_coriolis.form = "sadourny_energy"
    cfg.ocean_hvisc.nu_h = 2000.0
    cfg.ocean_hvisc.lateral_closure = "smagorinsky"
    cfg.ocean_hvisc.c_smag = 0.15
    cfg.ocean_hvisc.smag_ah = True
    cfg.ocean_hvisc.smag_bi_const = 0.06
    cfg.ocean_hvisc.bound_kh = True
    cfg.ocean_bdrag.form = "quadratic"
    cfg.ocean_bdrag.cd = 3.0e-3
    cfg.ocean_bdrag.hbbl = 10.0
    cfg.ocean_bdrag.bg_vel = 0.1
    cfg.ocean_bt.auto_n_inner = True
    cfg.ocean_eos.eos = "wright"

    # Daily top-10 m means of u, v, T, S plus SSH — the movie's source.
    cfg.ocean_diag.enabled = True
    cfg.ocean_diag.filename = "global_1deg"
    cfg.ocean_diag.dt_out = 1.0
    cfg.ocean_diag.vgrid = "z_fixed"
    cfg.ocean_diag.z_levels = [10.0]
    cfg.ocean_diag.n_z_levels = 1
    cfg.ocean_diag.output_precision = "single"
    cfg.ocean_diag.diags = "KE:off"
    cfg.logging.status_interval = 1.0
    return cfg


def check_against_namelist(cfg):
    """Every knob of the shipped namelist, and nothing else, is set here
    with the same value (file paths and run length excepted)."""
    ref = rdb.Model.from_namelist(NAMELIST, defer=True)
    try:
        mine = {(g, k): v for g, k, v in cfg.explicit_knobs()}
        theirs = {(g, k): v for g, k, v in ref.config.explicit_knobs()}
    finally:
        ref.close()
    skip = {("ocean_grid", "supergrid_file"), ("output", "bathymetry_file"),
            ("ocean_zinit", "file"), ("output", "output_dir"), ("time", "t_end")}
    diff = sorted(k for k in set(mine) | set(theirs)
                  if k not in skip and mine.get(k) != theirs.get(k))
    if diff:
        raise SystemExit("Config API and global_1deg_unforced.nml disagree on: " +
                         ", ".join(f"&{g}_nml {k} ({mine.get((g, k))!r} vs "
                                   f"{theirs.get((g, k))!r})" for g, k in diff))
    print(f"config: {len(mine)} knobs, identical to {os.path.basename(NAMELIST)}")


def max_face_velocity(field, ng):
    """(max |value|, i, j, k) over a face-velocity Field, interior indices
    (0-based, k bottom-up).  Reads the whole ghosted array once through the
    stdlib fast path (`_flat`, one ctypes block copy) rather than nesting
    5.9 M values into lists."""
    f = field.with_halo
    n0, n1, _ = f.shape
    flat = f._flat()
    best, at = -1.0, 0
    for idx, val in enumerate(flat):
        a = -val if val < 0.0 else val
        if a > best:
            best, at = a, idx
    k, r = divmod(at, n0 * n1)
    j, i = divmod(r, n0)
    return best, i - ng, j - ng, k


def flat_jmajor(nested, nx, ny):
    """DiagnosticField[...] nesting [i][j][k] -> flat j-major list (k = 0)."""
    return [nested[i][j][0] for j in range(ny) for i in range(nx)]


def main():
    root = os.environ.get("RDB_DATA_DIR")
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--data", default=os.path.join(root, "OM_1deg") if root else None,
                    help="OM_1deg data directory (default $RDB_DATA_DIR/OM_1deg)")
    ap.add_argument("--out", default="global_1deg_py", help="output directory")
    ap.add_argument("--frame-every", type=int, default=1, help="days between frames")
    ap.add_argument("--no-movie", action="store_true")
    ap.add_argument("--no-check", action="store_true",
                    help="skip the Config-vs-namelist comparison")
    ap.add_argument("--max-velocity", action="store_true",
                    help="also log the daily 3-D max |u|, |v| and where (~2 s/day)")
    a = ap.parse_args()
    if not a.data:
        ap.error("pass --data DIR or set RDB_DATA_DIR (see tools/fetch_om1deg.py)")
    data = os.path.abspath(a.data)
    out = os.path.abspath(a.out)
    frames = os.path.join(out, "frames")
    os.makedirs(frames, exist_ok=True)

    cfg = build_config(data, a.days, out)
    if not a.no_check:
        check_against_namelist(cfg)
    renderer = None if a.no_movie else FrameRenderer(
        os.path.join(data, "ocean_hgrid.nc"), os.path.join(data, "bathy_om1deg.nc"))

    log = open(os.path.join(out, "daily.txt"), "w")
    log.write("# day  En[m2/s2]  mass_drift_rel  max_surface_speed[m/s]  wall[s]"
              + ("  max|u| i j k lon lat  max|v| i j k lon lat" if a.max_velocity else "") + "\n")
    if a.max_velocity:
        _, _, lon, lat = model_tpoints(os.path.join(data, "ocean_hgrid.nc"))
    t_start = time.time()
    with rdb.Model(cfg) as model:
        info = model.grid_info
        nx, ny = info["nx"], info["ny"]
        spd = int(round(86400.0 / cfg.time.dt_fixed))
        mass0 = model.total_mass
        nframe = 0
        for day in range(1, a.days + 1):
            t0 = time.time()
            model.step(spd)
            ke, mass = model.kinetic_energy, model.total_mass
            en = ke / (mass / RHO0_API)
            u = flat_jmajor(model.diagnostic("u")[...], nx, ny)
            v = flat_jmajor(model.diagnostic("v")[...], nx, ny)
            speed = speed_from_uv(u, v)
            smax = max((s for s, w in zip(speed, renderer.wet if renderer else speed) if w),
                       default=0.0)
            if not all(math.isfinite(x) for x in (ke, mass, smax)):
                raise SystemExit(f"day {day}: non-finite state (KE {ke}, mass {mass})")
            if renderer and (day - 1) % a.frame_every == 0:
                ssh = flat_jmajor(model.diagnostic("SSH")[...], nx, ny)
                nframe += 1
                renderer.render(day, speed, ssh).write_ppm(
                    os.path.join(frames, f"frame_{nframe:04d}.ppm"))
            line = (f"{day:5d}  {en:.6e}  {(mass - mass0) / mass0:+.3e}  {smax:.4f}  "
                    f"{time.time() - t0:.1f}")
            if a.max_velocity:
                for fld in (model.u, model.v):
                    vmax, i, j, k = max_face_velocity(fld, info["nghost"])
                    ic, jc = min(max(i, 0), nx - 1), min(max(j, 0), ny - 1)
                    line += (f"  {vmax:.3f} {i} {j} {k} {lon[jc][ic]:.1f} {lat[jc][ic]:.1f}")
            log.write(line + "\n")
            log.flush()
            print(f"[py] day {line}", flush=True)
    print(f"{a.days} days in {time.time() - t_start:.0f} s")
    if renderer and nframe:
        mp4, gif = encode(frames, os.path.join(out, "global_1deg"))
        print(f"movie: {mp4}\n       {gif}")


if __name__ == "__main__":
    main()
