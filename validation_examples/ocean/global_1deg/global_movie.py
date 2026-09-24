"""Movie frames for the global 1-degree case — standard library only.

Renders a two-panel frame (surface speed over sea-surface height) of the
tripolar model grid on a regular longitude-latitude image, writes it as a
binary PPM, and encodes a frame directory to MP4 (H.264, yuv420p) and GIF
with ``ffmpeg``.  No numpy / matplotlib: the projection is a nearest-cell
lookup table built once, the colour maps are anchor tables, and the labels
use a built-in 5x7 bitmap font.

Used two ways:

* in-process by ``run_global_1deg.py`` (fields pulled straight out of the
  running model through the ``rdb`` Python interface), and
* offline on the ``rdb`` executable's diagnostic file
  (``python3 global_movie.py DIAG.nc OUTDIR``; the NetCDF-4 file is first
  flattened to NetCDF-3 with ``nccopy -u -k 64-bit-offset``).
"""

import collections
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools")))
from om1deg_prepare_inputs import NC3, model_tpoints  # noqa: E402

# ----------------------------------------------------------------------------
# Colour maps (anchor tables, linear in RGB between anchors)
# ----------------------------------------------------------------------------
INFERNO = [(0, 0, 4), (40, 11, 84), (101, 21, 110), (159, 42, 99), (212, 72, 66),
           (245, 125, 21), (250, 193, 39), (252, 255, 164)]
BALANCE = [(24, 28, 67), (12, 94, 190), (117, 155, 207), (205, 218, 232),
           (241, 236, 235), (231, 197, 184), (208, 126, 104), (164, 38, 41),
           (60, 9, 18)]
LAND = (150, 150, 150)
BG = (255, 255, 255)
INK = (30, 30, 30)


def make_lut(anchors, n=256):
    lut = []
    for s in range(n):
        x = s / (n - 1) * (len(anchors) - 1)
        a = min(int(x), len(anchors) - 2)
        f = x - a
        lut.append(bytes(int(round(anchors[a][c] + f * (anchors[a + 1][c] - anchors[a][c])))
                         for c in range(3)))
    return lut


# ----------------------------------------------------------------------------
# 5x7 bitmap font (the characters the labels need)
# ----------------------------------------------------------------------------
_FONT = {
    "A": "01110 10001 10001 11111 10001 10001 10001",
    "B": "11110 10001 10001 11110 10001 10001 11110",
    "C": "01110 10001 10000 10000 10000 10001 01110",
    "D": "11110 10001 10001 10001 10001 10001 11110",
    "E": "11111 10000 10000 11110 10000 10000 11111",
    "F": "11111 10000 10000 11110 10000 10000 10000",
    "G": "01110 10001 10000 10111 10001 10001 01111",
    "H": "10001 10001 10001 11111 10001 10001 10001",
    "I": "01110 00100 00100 00100 00100 00100 01110",
    "J": "00111 00010 00010 00010 00010 10010 01100",
    "K": "10001 10010 10100 11000 10100 10010 10001",
    "L": "10000 10000 10000 10000 10000 10000 11111",
    "M": "10001 11011 10101 10101 10001 10001 10001",
    "N": "10001 11001 10101 10011 10001 10001 10001",
    "O": "01110 10001 10001 10001 10001 10001 01110",
    "P": "11110 10001 10001 11110 10000 10000 10000",
    "Q": "01110 10001 10001 10001 10101 10010 01101",
    "R": "11110 10001 10001 11110 10100 10010 10001",
    "S": "01111 10000 10000 01110 00001 00001 11110",
    "T": "11111 00100 00100 00100 00100 00100 00100",
    "U": "10001 10001 10001 10001 10001 10001 01110",
    "V": "10001 10001 10001 10001 10001 01010 00100",
    "W": "10001 10001 10001 10101 10101 10101 01010",
    "X": "10001 10001 01010 00100 01010 10001 10001",
    "Y": "10001 10001 01010 00100 00100 00100 00100",
    "Z": "11111 00001 00010 00100 01000 10000 11111",
    "0": "01110 10001 10011 10101 11001 10001 01110",
    "1": "00100 01100 00100 00100 00100 00100 01110",
    "2": "01110 10001 00001 00010 00100 01000 11111",
    "3": "11110 00001 00001 01110 00001 00001 11110",
    "4": "00010 00110 01010 10010 11111 00010 00010",
    "5": "11111 10000 11110 00001 00001 10001 01110",
    "6": "00110 01000 10000 11110 10001 10001 01110",
    "7": "11111 00001 00010 00100 01000 01000 01000",
    "8": "01110 10001 10001 01110 10001 10001 01110",
    "9": "01110 10001 10001 01111 00001 00010 01100",
    " ": "00000 00000 00000 00000 00000 00000 00000",
    ".": "00000 00000 00000 00000 00000 01100 01100",
    ",": "00000 00000 00000 00000 01100 00100 01000",
    "-": "00000 00000 00000 11111 00000 00000 00000",
    "+": "00000 00100 00100 11111 00100 00100 00000",
    "(": "00010 00100 01000 01000 01000 00100 00010",
    ")": "01000 00100 00010 00010 00010 00100 01000",
    "/": "00001 00010 00010 00100 01000 01000 10000",
    ":": "00000 01100 01100 00000 01100 01100 00000",
    "|": "00100 00100 00100 00100 00100 00100 00100",
}


class Canvas:
    """RGB byte canvas with a PPM writer."""

    def __init__(self, width, height, fill=BG):
        self.w, self.h = width, height
        self.px = bytearray(bytes(fill) * (width * height))

    def set(self, x, y, rgb):
        if 0 <= x < self.w and 0 <= y < self.h:
            o = 3 * (y * self.w + x)
            self.px[o:o + 3] = rgb

    def rect(self, x0, y0, x1, y1, rgb):
        for y in range(max(y0, 0), min(y1, self.h)):
            o = 3 * (y * self.w + max(x0, 0))
            n = min(x1, self.w) - max(x0, 0)
            if n > 0:
                self.px[o:o + 3 * n] = bytes(rgb) * n

    def text(self, x, y, s, scale=2, rgb=INK):
        for ch in s.upper():
            rows = _FONT.get(ch, _FONT[" "]).split()
            for r, row in enumerate(rows):
                for c, bit in enumerate(row):
                    if bit == "1":
                        self.rect(x + c * scale, y + r * scale,
                                  x + (c + 1) * scale, y + (r + 1) * scale, rgb)
            x += 6 * scale
        return x

    @staticmethod
    def text_width(s, scale=2):
        return 6 * scale * len(s)

    def write_ppm(self, path):
        with open(path, "wb") as f:
            f.write(b"P6\n%d %d\n255\n" % (self.w, self.h))
            f.write(self.px)


# ----------------------------------------------------------------------------
# Tripolar grid -> regular lon-lat image
# ----------------------------------------------------------------------------
class LatLonMap:
    """Nearest-model-cell lookup for every pixel of a lon-lat panel.

    ``lon``/``lat`` are the model T-point coordinates ([j][i] nested lists,
    ``model_tpoints``); ``wet`` a flat j-major list of 0/1.  Each model cell
    is dropped into the pixel holding its centre, then a breadth-first flood
    (periodic in longitude) gives every other pixel its nearest binned cell
    — the tripolar cap and the ordinary lon-lat band alike, with no special
    case for the fold.  Pixels whose cell is land draw grey.
    """

    def __init__(self, lon, lat, wet, width=1080, lon0=20.0, lat_s=-80.0, lat_n=90.0):
        ny, nx = len(lon), len(lon[0])
        self.w = width
        self.ppd = width / 360.0
        self.h = int(round((lat_n - lat_s) * self.ppd))
        self.lon0, self.lat_n = lon0, lat_n
        owner = [-1] * (self.w * self.h)
        q = collections.deque()
        for j in range(ny):
            for i in range(nx):
                x = int(((lon[j][i] - lon0) % 360.0) * self.ppd)
                y = int((lat_n - lat[j][i]) * self.ppd)
                if 0 <= y < self.h:
                    p = y * self.w + min(x, self.w - 1)
                    if owner[p] < 0:
                        owner[p] = j * nx + i
                        q.append(p)
        while q:
            p = q.popleft()
            y, x = divmod(p, self.w)
            for yy, xx in ((y, (x + 1) % self.w), (y, (x - 1) % self.w), (y + 1, x), (y - 1, x)):
                if 0 <= yy < self.h:
                    n = yy * self.w + xx
                    if owner[n] < 0:
                        owner[n] = owner[p]
                        q.append(n)
        self.owner = [o if wet[o] > 0 else -1 for o in owner]

    def draw(self, canvas, x0, y0, values, lut, vmin, vmax, gamma=1.0):
        """Blit one field (flat j-major, model shape) at (x0, y0).

        `gamma` < 1 stretches the low end: colour index ~ ((v-vmin)/range)**gamma.
        """
        n = len(lut) - 1
        if gamma != 1.0:
            values = [((max(v - vmin, 0.0) / (vmax - vmin)) ** gamma) * (vmax - vmin) + vmin
                      if v == v else v for v in values]
        scale = n / (vmax - vmin)
        land = bytes(LAND)
        row = bytearray(3 * self.w)
        for y in range(self.h):
            base = y * self.w
            for x in range(self.w):
                o = self.owner[base + x]
                if o < 0:
                    row[3 * x:3 * x + 3] = land
                else:
                    v = values[o]
                    s = int((v - vmin) * scale) if v == v else 0
                    row[3 * x:3 * x + 3] = lut[0 if s < 0 else (n if s > n else s)]
            off = 3 * ((y0 + y) * canvas.w + x0)
            canvas.px[off:off + 3 * self.w] = row


def colorbar(canvas, x0, y0, width, height, lut, vmin, vmax, ticks, fmt, gamma=1.0):
    n = len(lut) - 1
    for x in range(width):
        canvas.rect(x0 + x, y0, x0 + x + 1, y0 + height, lut[int(x / (width - 1) * n)])
    for t in ticks:
        x = x0 + int(((t - vmin) / (vmax - vmin)) ** gamma * (width - 1))
        canvas.rect(x, y0 + height, x + 1, y0 + height + 4, INK)
        s = fmt(t)
        canvas.text(x - Canvas.text_width(s) // 2, y0 + height + 7, s)


# ----------------------------------------------------------------------------
# The frame
# ----------------------------------------------------------------------------
SPEED_MAX = 0.5   # m/s
SPEED_GAMMA = 0.5  # square-root colour stretch: the unforced ocean is slow
SSH_MAX = 0.3     # m, symmetric about the global mean


class FrameRenderer:
    """Two stacked panels: surface speed (top 10 m) and SSH anomaly."""

    HEADER = 44
    BAR = 44

    def __init__(self, hgrid_path, bathy_path, width=1080):
        nx, ny, lon, lat = model_tpoints(hgrid_path)
        depth = NC3(bathy_path).read("depth")
        self.nx, self.ny = nx, ny
        self.wet = [1 if d > 0.0 else 0 for d in depth]
        self.map = LatLonMap(lon, lat, self.wet, width=width)
        self.w = self.map.w
        self.h = self.HEADER + 2 * (self.map.h + self.BAR)
        self.h += self.h % 2
        self.speed_lut = make_lut(INFERNO)
        self.ssh_lut = make_lut(BALANCE)

    def render(self, day, speed, ssh, title="ROUNDABOUT  GLOBAL 1 DEG TRIPOLAR  UNFORCED"):
        """`speed`, `ssh`: flat j-major lists of nx*ny values (land ignored)."""
        c = Canvas(self.w, self.h)
        c.text(12, 12, title, scale=3)
        label = f"DAY {day:5.1f}"
        c.text(self.w - Canvas.text_width(label, 3) - 12, 12, label, scale=3)
        wet_ssh = [s for s, m in zip(ssh, self.wet) if m and s == s]
        mean = sum(wet_ssh) / max(len(wet_ssh), 1)
        anom = [s - mean for s in ssh]
        y = self.HEADER
        self.map.draw(c, 0, y, speed, self.speed_lut, 0.0, SPEED_MAX, SPEED_GAMMA)
        y += self.map.h
        c.text(12, y + 10, "SURFACE SPEED  (M/S, TOP 10 M)")
        colorbar(c, self.w - 36 - 420, y + 6, 420, 12, self.speed_lut, 0.0, SPEED_MAX,
                 [0.0, 0.02, 0.1, 0.2, 0.3, 0.5], lambda t: f"{t:g}", SPEED_GAMMA)
        y += self.BAR
        self.map.draw(c, 0, y, anom, self.ssh_lut, -SSH_MAX, SSH_MAX)
        y += self.map.h
        c.text(12, y + 10, "SEA SURFACE HEIGHT  (M, GLOBAL MEAN REMOVED)")
        colorbar(c, self.w - 36 - 420, y + 6, 420, 12, self.ssh_lut, -SSH_MAX, SSH_MAX,
                 [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3], lambda t: f"{t:+.1f}" if t else "0")
        return c


def speed_from_uv(u, v):
    return [math.sqrt(a * a + b * b) for a, b in zip(u, v)]


# ----------------------------------------------------------------------------
# Encoding
# ----------------------------------------------------------------------------
def encode(frame_dir, out_stem, fps=24, gif_width=540, gif_fps=12):
    """frame_%04d.ppm in `frame_dir` -> `out_stem`.mp4 (H.264) + .gif."""
    ffmpeg = shutil.which("ffmpeg") or "/usr/bin/ffmpeg"
    pattern = os.path.join(frame_dir, "frame_%04d.ppm")
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pattern,
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-movflags", "+faststart", out_stem + ".mp4"], check=True)
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pattern,
                    "-vf", f"fps={gif_fps},scale={gif_width}:-1:flags=lanczos,split[a][b];"
                           "[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=bayer",
                    out_stem + ".gif"], check=True)
    return out_stem + ".mp4", out_stem + ".gif"


# ----------------------------------------------------------------------------
# Offline: frames from the rdb executable's diagnostic file
# ----------------------------------------------------------------------------
def frames_from_diag(diag_nc, out_dir, data_dir, every=1):
    """Render every `every`-th daily frame of `diag_nc` into `out_dir`."""
    os.makedirs(out_dir, exist_ok=True)
    flat = os.path.join(out_dir, "diag_nc3.nc")
    subprocess.run([shutil.which("nccopy") or "nccopy", "-u", "-k", "64-bit-offset", "-V",
                    "time_SSH,SSH,time_u,u,time_v,v", diag_nc, flat], check=True)
    nc = NC3(flat)
    r = FrameRenderer(os.path.join(data_dir, "ocean_hgrid.nc"),
                      os.path.join(data_dir, "bathy_om1deg.nc"))
    nt, nyg, nxg = nc.shape("SSH")
    ng = (nxg - r.nx) // 2
    t_ssh = nc.read("time_SSH")
    ssh_all, u_all, v_all = nc.read("SSH"), nc.read("u"), nc.read("v")

    def interior(a, t):
        base = t * nyg * nxg
        return [a[base + (j + ng) * nxg + i + ng] for j in range(r.ny) for i in range(r.nx)]

    n = 0
    for t in range(0, min(nt, nc.shape("u")[0]), every):
        spd = speed_from_uv(interior(u_all, t), interior(v_all, t))
        frame = r.render(t_ssh[t] / 86400.0, spd, interior(ssh_all, t))
        n += 1
        frame.write_ppm(os.path.join(out_dir, f"frame_{n:04d}.ppm"))
    os.remove(flat)
    return n


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Movie from a global_1deg diagnostic file")
    ap.add_argument("diag_nc")
    ap.add_argument("out_dir")
    ap.add_argument("--data", default=os.path.join(os.environ.get("RDB_DATA_DIR", ""), "OM_1deg"))
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--stem", default="global_1deg")
    a = ap.parse_args()
    nfr = frames_from_diag(a.diag_nc, a.out_dir, a.data, a.every)
    print(f"{nfr} frames in {a.out_dir}")
    print("encoded:", *encode(a.out_dir, os.path.join(a.out_dir, a.stem)))
