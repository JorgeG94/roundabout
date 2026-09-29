"""South-polar-stereographic movie for the Southern Ocean 1-degree run
(stdlib only).

Renders daily frames of surface relative vorticity, normalised by the local
Coriolis parameter (zeta/|f|, a diverging colour map, symmetric limits from
the data's own robust percentiles) on a south-polar stereographic
projection centred on Antarctica and extending out to the domain's open
north edge (``--lat-edge``, default -30 deg, matching
``southern_ocean_1deg_wind.nml``'s sponge edge), with a small surface-speed
inset in the top-left corner, a day label, and -- when the diagnostic file
carries ``transport_x`` -- a Drake Passage transport time strip along the
bottom (the same idea as ``python_prototypes/global_1deg/make_movie.py``'s
``add_drake_strip``, folded in here so the whole movie is one reproducible
command). Reuses ``global_1deg/global_movie.py``'s ``Canvas`` /
``make_lut`` / bitmap font / colour anchor tables / ``encode`` (ffmpeg) and
``tools/om1deg_prepare_inputs.py``'s stdlib NetCDF-3 reader -- no numpy, no
netCDF4, nothing installed.

    python3 southern_movie.py \\
        RUN/output/southern_ocean_1deg_wind_rank_000000.nc OUT_DIR \\
        --data /home/jorge/nci/cdx/data/OM_1deg_southern

Needs ``nccopy`` (NetCDF-C) on PATH to flatten the NetCDF-4 diagnostic file,
and ``ffmpeg`` to encode. ``--fps60`` adds a 60 fps ``minterpolate``
motion-smoothed MP4 on top of the plain 24 fps one (the technique in
``python_prototypes/global_1deg/make_movie.py``).

**Known simplification**: the tripolar->image remap here is NEAREST-CELL
(one flood-filled owner per pixel, the same algorithm
``global_movie.LatLonMap`` uses), not the bilinear-in-index-space remap
``global_movie.BilinearLatLonMap`` does for the equirectangular movie --
the tangent-plane inversion there is written for a lon-lat rectangle, and a
bilinear polar version needs its own pixel -> (lat, lon) inverse. At 1 deg
resolution the visible difference is a slightly blockier coastline, not
missing structure; a bilinear polar remap is a reasonable follow-up.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "global_1deg")))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools")))
from global_movie import (Canvas, make_lut, colorbar, encode,  # noqa: E402
                          INFERNO, BALANCE, LAND, BG, INK, speed_from_uv)
from om1deg_prepare_inputs import NC3, model_tpoints  # noqa: E402

OMEGA = 7.2921159e-5  # Earth rotation rate, rad/s


class PolarSouthMap:
    """Nearest-model-cell lookup for a south-polar stereographic image.

    Standard polar stereographic, true at the pole: a point at colatitude
    ``c`` (angular distance from the south pole, degrees) sits at pixel
    radius ``rho = R*tan(c/2)`` -- so ``lat_edge`` maps exactly to
    ``radius_px``. Same "drop each model cell in its pixel, then
    breadth-first flood fill" approach as ``global_movie.LatLonMap``, just
    without the periodic x-wrap (a polar image is not periodic in pixel
    space) and gated on ``lat <= lat_edge`` (poleward of the domain edge
    only -- there is no data north of it to draw).
    """

    def __init__(self, lon, lat, wet, size=900, lat_edge=-30.0, margin=40):
        ny, nx = len(lon), len(lon[0])
        self.w = self.h = size
        cx = cy = size / 2.0
        radius_px = size / 2.0 - margin
        colat_edge = math.radians(90.0 + lat_edge)
        r_scale = radius_px / math.tan(colat_edge / 2.0)
        self.cx, self.cy, self.r_scale = cx, cy, r_scale
        self.lat_edge = lat_edge
        owner = [-1] * (size * size)
        q = []
        for j in range(ny):
            for i in range(nx):
                la = lat[j][i]
                if la > lat_edge:
                    continue
                px, py = self._project(la, lon[j][i])
                x, y = int(px), int(py)
                if 0 <= x < size and 0 <= y < size:
                    p = y * size + x
                    if owner[p] < 0:
                        owner[p] = j * nx + i
                        q.append(p)
        head = 0
        while head < len(q):
            p = q[head]
            head += 1
            y, x = divmod(p, size)
            for yy, xx in ((y, x + 1), (y, x - 1), (y + 1, x), (y - 1, x)):
                if 0 <= xx < size and 0 <= yy < size:
                    n = yy * size + xx
                    if owner[n] < 0:
                        owner[n] = owner[p]
                        q.append(n)
        owner = [o if (o >= 0 and wet[o] > 0) else -1 for o in owner]
        # The flood fill above has no radius cutoff, so pixels in the four
        # corners of the square canvas (well outside the projected data
        # disk -- nothing was ever seeded there) inherit whichever seeded
        # pixel's wavefront reaches them first, via the shortest path along
        # the boundary ring: a blocky, wedge-shaped LAND-coloured pattern
        # with no data behind it.  Mask anything past the disk radius back
        # to a separate "outside" sentinel (-2) so `draw` paints it as
        # background instead of land.
        r2 = (radius_px + 1.0) ** 2
        self.outside = [((x - cx) ** 2 + (y - cy) ** 2) > r2
                        for y in range(size) for x in range(size)]
        self.owner = owner

    def _project(self, lat_deg, lon_deg):
        colat = math.radians(90.0 + lat_deg)
        rho = self.r_scale * math.tan(colat / 2.0)
        theta = math.radians(lon_deg)
        return (self.cx + rho * math.sin(theta), self.cy - rho * math.cos(theta))

    def draw(self, canvas, x0, y0, values, lut, vmin, vmax):
        n = len(lut) - 1
        scale = n / (vmax - vmin)
        land = bytes(LAND)
        bg = bytes(BG)
        row = bytearray(3 * self.w)
        for y in range(self.h):
            base = y * self.w
            for x in range(self.w):
                p = base + x
                if self.outside[p]:
                    row[3 * x:3 * x + 3] = bg
                    continue
                o = self.owner[p]
                if o < 0:
                    row[3 * x:3 * x + 3] = land
                else:
                    v = values[o]
                    s = int((v - vmin) * scale) if v == v else 0
                    row[3 * x:3 * x + 3] = lut[0 if s < 0 else (n if s > n else s)]
            off = 3 * ((y0 + y) * canvas.w + x0)
            canvas.px[off:off + 3 * self.w] = row


def percentile(sorted_vals, p):
    if not sorted_vals:
        return 1.0
    k = min(max(int(p * (len(sorted_vals) - 1)), 0), len(sorted_vals) - 1)
    return sorted_vals[k]


def add_drake_strip(frame, day, series, strip=140):
    """New canvas = `frame` + a strip with the Drake transport up to `day`."""
    c = Canvas(frame.w, frame.h + strip)
    c.px[:len(frame.px)] = frame.px
    y0 = frame.h
    x0, x1 = 90, frame.w - 20
    ya, yb = y0 + 30, y0 + strip - 20
    ndays = max(series) if series else 1
    vals = list(series.values()) + [0.0]
    lo, hi = min(vals), max(vals)
    pad = 0.08 * (hi - lo or 1.0)
    lo, hi = lo - pad, hi + pad

    def px(d):
        return x0 + round((d - 1) * (x1 - x0) / max(ndays - 1, 1))

    def py(v):
        return yb - round((v - lo) * (yb - ya) / (hi - lo))

    grid = (70, 70, 80)
    c.rect(x0, ya, x1 + 1, ya + 1, grid)
    c.rect(x0, yb, x1 + 1, yb + 1, grid)
    z = py(0.0)
    for x in range(x0, x1, 4):
        c.rect(x, z, x + 2, z + 1, grid)
    pts = [(px(d), py(v)) for d, v in sorted(series.items()) if d <= day]
    for (xa, ya_), (xb, yb_) in zip(pts, pts[1:]):
        n = max(abs(xb - xa), abs(yb_ - ya_), 1)
        for k in range(n + 1):
            x = xa + (xb - xa) * k // n
            y = ya_ + (yb_ - ya_) * k // n
            c.rect(x, y, x + 3, y + 3, (255, 200, 60))
    d = min(max(int(round(day)), 1), ndays)
    cur = series.get(d, series.get(float(d), 0.0))
    lab = f"DRAKE PASSAGE TRANSPORT (SV, 67.5W)   DAY {d:4d}: {cur:+6.1f}"
    c.text(12, y0 + 8, lab)
    return c


class SouthernFrameRenderer:
    HEADER = 48
    INSET = 220

    def __init__(self, hgrid_path, bathy_path, size=900, lat_edge=-30.0,
                 title="ROUNDABOUT  SOUTHERN OCEAN 1 DEG  JRA55-DO WIND"):
        nx, ny, lon, lat = model_tpoints(hgrid_path)
        depth = NC3(bathy_path).read("depth")
        self.nx, self.ny = nx, ny
        self.wet = [1 if d > 0.0 else 0 for d in depth]
        self.map = PolarSouthMap(lon, lat, self.wet, size=size, lat_edge=lat_edge)
        self.speed_map = PolarSouthMap(lon, lat, self.wet, size=self.INSET + 20, lat_edge=lat_edge)
        self.title = title
        self.w = self.map.w
        self.h = self.HEADER + self.map.h + 40
        self.h += self.h % 2
        self.vort_lut = make_lut(BALANCE)
        self.speed_lut = make_lut(INFERNO)
        # Coriolis parameter per T-cell, |f| = 2*Omega*|sin(lat)|.
        self.f_abs = [2.0 * OMEGA * abs(math.sin(math.radians(lat[j][i])))
                     for j in range(ny) for i in range(nx)]

    def render(self, day, zeta, speed, vmin, vmax, speed_max):
        c = Canvas(self.w, self.h)
        c.text(12, 12, self.title, scale=2)
        label = f"DAY {day:6.1f}"
        c.text(self.w - Canvas.text_width(label, 2) - 12, 12, label, scale=2)
        ratio = [z / f if f > 1e-12 else float("nan") for z, f in zip(zeta, self.f_abs)]
        y = self.HEADER
        self.map.draw(c, 0, y, ratio, self.vort_lut, vmin, vmax)
        # Speed inset, top-left corner, boxed.
        ix0, iy0 = 16, y + 16
        self.speed_map.draw(c, ix0, iy0, speed, self.speed_lut, 0.0, speed_max)
        box = (255, 255, 255)
        c.rect(ix0 - 2, iy0 - 2, ix0 + self.speed_map.w + 2, iy0 - 1, box)
        c.rect(ix0 - 2, iy0 + self.speed_map.h + 1, ix0 + self.speed_map.w + 2,
               iy0 + self.speed_map.h + 3, box)
        c.rect(ix0 - 2, iy0 - 2, ix0 - 1, iy0 + self.speed_map.h + 2, box)
        c.rect(ix0 + self.speed_map.w + 1, iy0 - 2, ix0 + self.speed_map.w + 3,
               iy0 + self.speed_map.h + 2, box)
        c.text(ix0, iy0 + self.speed_map.h + 6, "SPEED (M/S, TOP 10 M)", scale=1)
        y += self.map.h
        c.text(12, y + 8, "SURFACE RELATIVE VORTICITY / |f|  (DIMENSIONLESS)")
        colorbar(c, self.w - 24 - 420, y + 4, 420, 12, self.vort_lut, vmin, vmax,
                [vmin, vmin / 2, 0.0, vmax / 2, vmax], lambda t: f"{t:+.2f}")
        return c


def flatten(diag, out, names):
    subprocess.run([shutil.which("nccopy") or "nccopy", "-u", "-k", "64-bit-offset",
                    "-V", ",".join(names), diag, out], check=True)


def frames_from_diag(diag_nc, out_dir, data_dir, every=1, size=900, lat_edge=-30.0,
                     speed_max=1.0, title=None, drake_lon=-67.5, drake_lat=(-70.0, -52.0)):
    os.makedirs(out_dir, exist_ok=True)
    flat = os.path.join(out_dir, "diag_nc3.nc")
    names = ["time_SSH", "SSH", "time_u", "u", "time_v", "v",
             "time_vorticity_z", "vorticity_z", "time_transport_x", "transport_x"]
    flatten(diag_nc, flat, names)
    nc = NC3(flat)
    kw = {"title": title} if title else {}
    r = SouthernFrameRenderer(os.path.join(data_dir, "ocean_hgrid.nc"),
                              os.path.join(data_dir, "bathy_om1deg.nc"),
                              size=size, lat_edge=lat_edge, **kw)
    nt, nyg, nxg = nc.shape("vorticity_z")
    ng = (nxg - r.nx) // 2
    t_ssh = nc.read("time_SSH")
    zeta_all, u_all, v_all = nc.read("vorticity_z"), nc.read("u"), nc.read("v")
    has_transport = "transport_x" in nc.vars

    def interior(a, t):
        base = t * nyg * nxg
        return [a[base + (j + ng) * nxg + i + ng] for j in range(r.ny) for i in range(r.nx)]

    # ---- robust colour limits: sample every 5th day's |zeta/f| ----
    samples = []
    for t in range(0, nt, max(1, nt // 100 or 1)):
        z = interior(zeta_all, t)
        for zz, f in zip(z, r.f_abs):
            if zz == zz and f > 1e-12:
                samples.append(abs(zz / f))
    samples.sort()
    vmax = max(percentile(samples, 0.98), 1e-3)
    vmin = -vmax

    # ---- Drake Passage transport series, if the diag carries it ----
    drake_series = {}
    if has_transport:
        gxy = NC3(os.path.join(data_dir, "ocean_hgrid.nc"))
        nyp, nxp = gxy.shape("x")
        X, Y, DY = gxy.read("x"), gxy.read("y"), gxy.read("dy")
        nxm = (nxp - 1) // 2
        lon_col = [X[(2 * 0 + 1) * nxp + 2 * i + 1] for i in range(nxm)]
        i_d = min(range(nxm), key=lambda i: abs(((lon_col[i] + 180) % 360 - 180)
                                                 - ((drake_lon + 180) % 360 - 180)))
        lat_col = [Y[(2 * j + 1) * nxp + 2 * i_d + 1] for j in range(r.ny)]
        dyt_col = [DY[(2 * j) * nxp + 2 * i_d + 1] + DY[(2 * j + 1) * nxp + 2 * i_d + 1]
                  for j in range(r.ny)]
        sec = [j for j in range(r.ny) if drake_lat[0] <= lat_col[j] <= drake_lat[1]
              and r.wet[j * r.nx + i_d]]
        tx_all = nc.read("transport_x")

        def drake(t):
            base = t * nyg * nxg
            return sum(tx_all[base + (j + ng) * nxg + i_d + ng] * dyt_col[j] for j in sec) / 1e6

        for t in range(nt):
            drake_series[float(t + 1)] = drake(t)

    n = 0
    for t in range(0, min(nt, nc.shape("u")[0]), every):
        zeta = interior(zeta_all, t)
        spd = speed_from_uv(interior(u_all, t), interior(v_all, t))
        day = t_ssh[t] / 86400.0
        frame = r.render(day, zeta, spd, vmin, vmax, speed_max)
        if drake_series:
            frame = add_drake_strip(frame, day, drake_series)
        n += 1
        frame.write_ppm(os.path.join(out_dir, f"frame_{n:04d}.ppm"))
    os.remove(flat)
    return n, drake_series


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("diag_nc")
    ap.add_argument("out_dir")
    ap.add_argument("--data", required=True, help="OM_1deg_southern data directory")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--stem", default="southern_ocean_1deg_wind")
    ap.add_argument("--size", type=int, default=900, help="square canvas side, px")
    ap.add_argument("--lat-edge", type=float, default=-30.0)
    ap.add_argument("--speed-max", type=float, default=1.0)
    ap.add_argument("--title", default=None)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--fps60", action="store_true", help="also render a 60 fps minterpolate MP4")
    ap.add_argument("--keep-frames", action="store_true")
    a = ap.parse_args()
    nfr, drake = frames_from_diag(a.diag_nc, a.out_dir, a.data, a.every, a.size, a.lat_edge,
                                  a.speed_max, a.title)
    print(f"{nfr} frames in {a.out_dir}; Drake series: {len(drake)} days")
    mp4, gif = encode(a.out_dir, os.path.join(a.out_dir, a.stem), fps=a.fps)
    print("encoded:", mp4, gif)
    if a.fps60:
        vf = "minterpolate=fps=60:mi_mode=mci:mc_mode=aobmc:me_mode=bidir:vsbmc=1"
        name = os.path.join(a.out_dir, a.stem + "_smooth60.mp4")
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i",
                        os.path.join(a.out_dir, "frame_%04d.ppm"), "-vf", vf, "-c:v", "libx264",
                        "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart", name],
                      check=True)
        print("encoded:", name)
    last = os.path.join(a.out_dir, f"frame_{nfr:04d}.ppm")
    png = os.path.join(a.out_dir, a.stem + "_last_frame.png")
    if shutil.which("convert"):
        subprocess.run(["convert", last, png], check=True)
        print("last frame:", png)
    if not a.keep_frames:
        for k in range(1, nfr + 1):
            p = os.path.join(a.out_dir, f"frame_{k:04d}.ppm")
            if os.path.exists(p) and p != last:
                os.remove(p)
