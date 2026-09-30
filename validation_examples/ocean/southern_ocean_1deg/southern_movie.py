"""South-polar-stereographic movie for the Southern Ocean 1-degree run
(stdlib only).

Renders daily frames on a south-polar stereographic projection centred on
Antarctica and extending out to the domain's open north edge (``--lat-edge``,
default -30 deg, matching ``southern_ocean_1deg_wind.nml``'s sponge edge).
Two panels are drawn, each independently selectable from
``{vorticity, speed, ssh}`` via ``--main``/``--inset`` (default
``main=speed``, ``inset=ssh``; the vorticity view from earlier versions of
this script is still available on either panel):

* ``vorticity`` -- surface relative vorticity normalised by the local
  Coriolis parameter (zeta/|f|, diverging colour map, symmetric limits from
  the run's own robust percentile of |zeta/f|).
* ``speed`` -- top-10 m current speed |u| (sequential colour map, limits
  ``[0, robust percentile of |u|]``).
* ``ssh`` -- sea-surface-height anomaly about the AREA-WEIGHTED domain mean,
  recomputed per frame (diverging colour map, symmetric limits from the
  run's own robust percentile of the anomaly).

All colour limits are sampled once across the run and held FIXED for every
frame. The main panel is the full-size disc; the inset is a proper-size
second disc placed in its own column to the left of the main disc (never
overlapping it), each with its own colorbar stating the field and its
units. A day label, the title, and -- when the diagnostic file carries
``transport_x`` -- a Drake Passage transport time strip along the bottom
(the same idea as ``python_prototypes/global_1deg/make_movie.py``'s
``add_drake_strip``, folded in here so the whole movie is one reproducible
command) round out the frame. Reuses ``global_1deg/global_movie.py``'s
``Canvas`` / ``make_lut`` / bitmap font / colour anchor tables / ``encode``
(ffmpeg) and ``tools/om1deg_prepare_inputs.py``'s stdlib NetCDF-3 reader --
no numpy, no netCDF4, nothing installed.

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
missing structure; a bilinear polar remap is a reasonable follow-up. This
grid is genuinely periodic in its i-index (checked directly: the longitude
step across the ``i = nx-1 -> i = 0`` seam is the same 1 degree as
everywhere else in the row) and the flood fill is correctly NON-periodic in
PIXEL space -- a polar disc image has no left/right wrap, unlike the
equirectangular map's pixel columns. What *did* need a fix (see
``PolarSouthMap``'s ``pole_void`` mask): this is a regional cut that stops
at the domain's south WALL (``lat ~ -77.8``, well short of the geographic
pole), so the area strictly poleward of that wall has no model data at
all, and the naive flood fill filled that gap with whichever boundary cell
happened to win the pixel-space race -- occasionally an open-ocean cell,
cutting a false "open ocean" wedge through what should be (and, at the
model's own south wall, IS) solid Antarctica.
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

# ----------------------------------------------------------------------------
# Field registry -- what each of --main/--inset can show
# ----------------------------------------------------------------------------
FIELD_META = {
    "vorticity": dict(label="SURFACE RELATIVE VORTICITY / |F|", units="DIMENSIONLESS",
                      cmap=BALANCE, diverging=True, pctl=0.98),
    "speed": dict(label="SURFACE SPEED", units="M/S, TOP 10 M",
                 cmap=INFERNO, diverging=False, pctl=0.995),
    "ssh": dict(label="SEA SURFACE HEIGHT ANOMALY", units="M, AREA-WEIGHTED MEAN REMOVED",
               cmap=BALANCE, diverging=True, pctl=0.995),
}


def vorticity_ratio(zeta, f_abs):
    # |zeta| > 1e10 is the diag's land sentinel (DIAG_MISSING_VALUE = 1e20).
    # The stdlib NC3 reader does not apply _FillValue, and files written
    # before vorticity_z advertised it on every path carry it unflagged, so
    # treat it as missing here rather than letting it set the colour range.
    return [z / f if f > 1e-12 and abs(z) < 1e10 else float("nan")
            for z, f in zip(zeta, f_abs)]


def area_weighted_mean(values, area, wet):
    num = den = 0.0
    for v, a, w in zip(values, area, wet):
        if w and v == v:
            num += v * a
            den += a
    return num / den if den > 0.0 else 0.0


def model_tarea(hgrid_path):
    """T-cell area (supergrid units), flat j-major list (ny, nx).

    MOM6 supergrid convention: ``area`` is defined on the same doubled grid
    as ``x``/``y`` but one row/column shorter (cell quadrants, not nodes) --
    a T-cell's area is the sum of its 4 quadrants.
    """
    nc = NC3(hgrid_path)
    nyp, nxp = nc.shape("x")
    nx, ny = (nxp - 1) // 2, (nyp - 1) // 2
    naxp = nxp - 1
    a = nc.read("area")
    area = []
    for j in range(ny):
        j0 = 2 * j
        for i in range(nx):
            i0 = 2 * i
            area.append(a[j0 * naxp + i0] + a[j0 * naxp + i0 + 1] +
                       a[(j0 + 1) * naxp + i0] + a[(j0 + 1) * naxp + i0 + 1])
    return area


def _fmt_for(scale, signed):
    """Tick formatter with enough decimals to be legible at this scale."""
    decimals = 2
    if scale < 0.1:
        decimals = 3
    if scale < 0.01:
        decimals = 4
    spec = f"{{:+.{decimals}f}}" if signed else f"{{:.{decimals}f}}"
    return lambda t: spec.format(t)


def panel_ticks_fmt(field, vmin, vmax):
    meta = FIELD_META[field]
    if meta["diverging"]:
        ticks = [vmin, vmin / 2, 0.0, vmax / 2, vmax]
    else:
        ticks = [vmin, vmax / 4, vmax / 2, 3 * vmax / 4, vmax]
    fmt = _fmt_for(max(abs(vmin), abs(vmax)), meta["diverging"])
    return ticks, fmt


class PolarSouthMap:
    """Nearest-model-cell lookup for a south-polar stereographic image.

    Standard polar stereographic, true at the pole: a point at colatitude
    ``c`` (angular distance from the south pole, degrees) sits at pixel
    radius ``rho = R*tan(c/2)`` -- so ``lat_edge`` maps exactly to
    ``radius_px``. Same "drop each model cell in its pixel, then
    breadth-first flood fill" approach as ``global_movie.LatLonMap``, just
    without the periodic x-wrap (a polar image is not periodic in pixel
    space -- there is no left/right seam to cross) and gated on
    ``lat <= lat_edge`` (poleward of the domain edge only -- there is no
    data north of it to draw).
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
        min_lat = None
        for j in range(ny):
            for i in range(nx):
                la = lat[j][i]
                if la > lat_edge:
                    continue
                if min_lat is None or la < min_lat:
                    min_lat = la
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
        # --- the pole-void seam fix -------------------------------------
        # This is a REGIONAL cut that stops at the domain's south wall
        # (``min_lat``, e.g. -77.8 deg for OM_1deg_southern) -- well short
        # of the geographic pole -- so nothing is ever seeded inside the
        # radius that row projects to (`rho_min` below): there is no model
        # data there, at any longitude. The flood fill still assigns those
        # pixels an owner (whichever boundary-ring cell's BFS wavefront
        # reaches them first over the shortest PIXEL path), and near the
        # projection singularity (rho -> 0 as colatitude -> 0) that is NOT
        # the same as shortest physical distance: a single open-ocean
        # boundary cell can win a wide fan of pixels reaching all the way
        # to the disc centre, cutting a false open-water wedge through
        # what should be solid Antarctica -- the seam this class used to
        # render. Force that whole no-data disc to LAND explicitly: it is
        # the physically correct call (south of the wall IS the
        # continent), and it is what >90% of that disc already resolves to
        # on its own -- this just closes the gap consistently instead of
        # leaving it to flood-fill happenstance.
        if min_lat is None:
            rho_min = 0.0
        else:
            colat_min = math.radians(90.0 + min_lat)
            rho_min = r_scale * math.tan(colat_min / 2.0)
        r2min = rho_min * rho_min
        self.pole_void = [((x - cx) ** 2 + (y - cy) ** 2) < r2min
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
                elif self.pole_void[p] or self.owner[p] < 0:
                    row[3 * x:3 * x + 3] = land
                else:
                    v = values[self.owner[p]]
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
    """Main disc (``main_field``) + a proper-size, non-overlapping inset
    disc (``inset_field``) in its own column to the left, each with its own
    colorbar naming the field and its units."""

    HEADER = 48
    INSET_SIZE = 300
    LEFT_MARGIN = 48
    GAP = 40
    RIGHT_MARGIN = 44
    MAIN_BAR_H = 64
    INSET_BAR_H = 46

    def __init__(self, hgrid_path, bathy_path, main_field="speed", inset_field="ssh",
                 size=900, lat_edge=-30.0,
                 title="ROUNDABOUT  SOUTHERN OCEAN 1 DEG  JRA55-DO WIND"):
        nx, ny, lon, lat = model_tpoints(hgrid_path)
        depth = NC3(bathy_path).read("depth")
        self.nx, self.ny = nx, ny
        self.wet = [1 if d > 0.0 else 0 for d in depth]
        self.area = model_tarea(hgrid_path)
        self.map = PolarSouthMap(lon, lat, self.wet, size=size, lat_edge=lat_edge)
        self.inset_map = PolarSouthMap(lon, lat, self.wet, size=self.INSET_SIZE,
                                       lat_edge=lat_edge)
        self.main_field, self.inset_field = main_field, inset_field
        self.title = title
        self.lut = {k: make_lut(v["cmap"]) for k, v in FIELD_META.items()}
        self.f_abs = [2.0 * OMEGA * abs(math.sin(math.radians(lat[j][i])))
                     for j in range(ny) for i in range(nx)]
        self.x_inset = self.LEFT_MARGIN
        self.x_main = self.x_inset + self.INSET_SIZE + 4 + self.GAP
        self.w = self.x_main + self.map.w + self.RIGHT_MARGIN
        inset_col_h = self.HEADER + self.inset_map.h + self.INSET_BAR_H
        main_col_h = self.HEADER + self.map.h + self.MAIN_BAR_H
        self.h = max(inset_col_h, main_col_h)
        self.h += self.h % 2
        self.ranges = {}

    def set_ranges(self, ranges):
        """`ranges`: {field_name: (vmin, vmax)}, fixed for the whole run."""
        self.ranges = ranges

    def render(self, day, values):
        """`values`: {field_name: flat j-major list} for at least
        `self.main_field` and `self.inset_field`."""
        c = Canvas(self.w, self.h)
        c.text(12, 12, self.title, scale=2)
        label = f"DAY {day:6.1f}"
        c.text(self.w - Canvas.text_width(label, 2) - 12, 12, label, scale=2)
        y0 = self.HEADER

        # Main panel (right of the inset column, full requested size).
        mmeta = FIELD_META[self.main_field]
        mvmin, mvmax = self.ranges[self.main_field]
        self.map.draw(c, self.x_main, y0, values[self.main_field],
                      self.lut[self.main_field], mvmin, mvmax)
        y_bot = y0 + self.map.h
        c.text(self.x_main, y_bot + 8, f"{mmeta['label']}  ({mmeta['units']})")
        mticks, mfmt = panel_ticks_fmt(self.main_field, mvmin, mvmax)
        colorbar(c, self.x_main + self.map.w - 420, y_bot + 28, 420, 12,
                self.lut[self.main_field], mvmin, mvmax, mticks, mfmt)

        # Inset panel: its own column, own box, own colorbar -- never
        # overlaps the main disc.
        imeta = FIELD_META[self.inset_field]
        ivmin, ivmax = self.ranges[self.inset_field]
        ix0, iy0 = self.x_inset, y0
        self.inset_map.draw(c, ix0, iy0, values[self.inset_field],
                            self.lut[self.inset_field], ivmin, ivmax)
        box = (255, 255, 255)
        iw, ih = self.inset_map.w, self.inset_map.h
        c.rect(ix0 - 2, iy0 - 2, ix0 + iw + 2, iy0 - 1, box)
        c.rect(ix0 - 2, iy0 + ih + 1, ix0 + iw + 2, iy0 + ih + 3, box)
        c.rect(ix0 - 2, iy0 - 2, ix0 - 1, iy0 + ih + 2, box)
        c.rect(ix0 + iw + 1, iy0 - 2, ix0 + iw + 3, iy0 + ih + 2, box)
        c.text(ix0, iy0 + ih + 8, imeta["label"], scale=1)
        c.text(ix0, iy0 + ih + 18, f"({imeta['units']})", scale=1)
        iticks, ifmt = panel_ticks_fmt(self.inset_field, ivmin, ivmax)
        colorbar(c, ix0, iy0 + ih + 32, iw, 10, self.lut[self.inset_field],
                ivmin, ivmax, iticks, ifmt)
        return c


def flatten(diag, out, names):
    subprocess.run([shutil.which("nccopy") or "nccopy", "-u", "-k", "64-bit-offset",
                    "-V", ",".join(names), diag, out], check=True)


def frames_from_diag(diag_nc, out_dir, data_dir, every=1, size=900, lat_edge=-30.0,
                     main_field="speed", inset_field="ssh", speed_max=None,
                     ssh_max=None, title=None, drake_lon=-67.5, drake_lat=(-70.0, -52.0),
                     bathy_name="bathy_om1deg.nc"):
    os.makedirs(out_dir, exist_ok=True)
    flat = os.path.join(out_dir, "diag_nc3.nc")
    names = ["time_SSH", "SSH", "time_u", "u", "time_v", "v",
             "time_vorticity_z", "vorticity_z", "time_transport_x", "transport_x"]
    flatten(diag_nc, flat, names)
    nc = NC3(flat)
    kw = {"title": title} if title else {}
    r = SouthernFrameRenderer(os.path.join(data_dir, "ocean_hgrid.nc"),
                              os.path.join(data_dir, bathy_name),
                              main_field=main_field, inset_field=inset_field,
                              size=size, lat_edge=lat_edge, **kw)
    nt, nyg, nxg = nc.shape("vorticity_z")
    ng = (nxg - r.nx) // 2
    t_ssh = nc.read("time_SSH")
    zeta_all, u_all, v_all, ssh_all = (nc.read("vorticity_z"), nc.read("u"),
                                       nc.read("v"), nc.read("SSH"))
    has_transport = "transport_x" in nc.vars

    def interior(a, t):
        base = t * nyg * nxg
        return [a[base + (j + ng) * nxg + i + ng] for j in range(r.ny) for i in range(r.nx)]

    # ---- robust, fixed-for-the-run colour limits ----
    # Sample ~100 frames spread over the run for each of the three fields
    # (cheap either way at this resolution), independent of which two are
    # actually selected for main/inset.
    vort_s, speed_s, ssh_s = [], [], []
    step = max(1, nt // 100 or 1)
    for t in range(0, nt, step):
        zeta = interior(zeta_all, t)
        vort_s.extend(abs(x) for x in vorticity_ratio(zeta, r.f_abs) if x == x)
        u, v = interior(u_all, t), interior(v_all, t)
        speed_s.extend(x for x in speed_from_uv(u, v) if x == x)
        ssh = interior(ssh_all, t)
        mean = area_weighted_mean(ssh, r.area, r.wet)
        ssh_s.extend(abs(s - mean) for s in ssh if s == s)
    vort_s.sort()
    speed_s.sort()
    ssh_s.sort()
    vort_max = max(percentile(vort_s, FIELD_META["vorticity"]["pctl"]), 1e-3)
    smax = speed_max if speed_max is not None else \
        max(percentile(speed_s, FIELD_META["speed"]["pctl"]), 1e-3)
    hmax = ssh_max if ssh_max is not None else \
        max(percentile(ssh_s, FIELD_META["ssh"]["pctl"]), 1e-3)
    r.set_ranges({"vorticity": (-vort_max, vort_max), "speed": (0.0, smax),
                 "ssh": (-hmax, hmax)})

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
        u, v = interior(u_all, t), interior(v_all, t)
        ssh = interior(ssh_all, t)
        mean = area_weighted_mean(ssh, r.area, r.wet)
        values = {
            "vorticity": vorticity_ratio(zeta, r.f_abs),
            "speed": speed_from_uv(u, v),
            "ssh": [s - mean if s == s else s for s in ssh],
        }
        day = t_ssh[t] / 86400.0
        frame = r.render(day, values)
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
    ap.add_argument("--data", required=True, help="OM_1deg_southern (or other cut) data directory")
    ap.add_argument("--bathy", default="bathy_om1deg.nc",
                    help="bathymetry file name inside --data (classic NetCDF-3, a 'depth(y,x)' "
                    "variable is all this script reads) -- e.g. 'ocean_topog.nc' for a raw "
                    "MOM6-style topog file instead of the preprocessed OM_1deg one")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--stem", default="southern_ocean_1deg_wind")
    ap.add_argument("--size", type=int, default=900, help="square canvas side, px")
    ap.add_argument("--lat-edge", type=float, default=-30.0)
    ap.add_argument("--main", choices=("vorticity", "speed", "ssh"), default="speed",
                    help="main-disc field (default: speed)")
    ap.add_argument("--inset", choices=("vorticity", "speed", "ssh"), default="ssh",
                    help="inset-disc field (default: ssh)")
    ap.add_argument("--speed-max", type=float, default=None,
                    help="speed colour-scale top, m/s (default: run's own "
                         f"{FIELD_META['speed']['pctl'] * 100:g}th percentile)")
    ap.add_argument("--ssh-max", type=float, default=None,
                    help="SSH-anomaly colour-scale +/- limit, m (default: run's own "
                         f"{FIELD_META['ssh']['pctl'] * 100:g}th percentile)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--fps60", action="store_true", help="also render a 60 fps minterpolate MP4")
    ap.add_argument("--keep-frames", action="store_true")
    a = ap.parse_args()
    nfr, drake = frames_from_diag(a.diag_nc, a.out_dir, a.data, a.every, a.size, a.lat_edge,
                                  a.main, a.inset, a.speed_max, a.ssh_max, a.title,
                                  bathy_name=a.bathy)
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
