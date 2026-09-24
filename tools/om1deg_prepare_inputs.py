#!/usr/bin/env python3
"""OM_1deg + WOA13 -> roundabout model-grid inputs (standard library only).

Writes, next to the fetched data (``tools/fetch_om1deg.py``):

``bathy_om1deg.nc``
    ``depth(y, x)`` (m, positive down) from MOM6's ``topog.nc`` with the
    OM_1deg ``MOM_input`` limits applied the way MOM6 applies them
    (``limit_topography`` with ``MASKING_DEPTH = 0``): every wet cell
    (``depth > 0``) is raised to ``MINIMUM_DEPTH`` (9.5 m) and capped at
    ``MAXIMUM_DEPTH`` (6500 m); ``depth <= 0`` stays land (0).

``ic_woa13_jan.nc``
    WOA13 decav January potential temperature (``ptemp_an``) and salinity
    (``s_an``) on the model's 360 x 320 tracer points at WOA's own 102
    standard depths — the input of ``&ocean_zinit_nml source = "file"``,
    which interpolates linearly in depth onto the model layers at run time.
    Horizontal: nearest neighbour (the WOA 1-degree cell containing the
    model T-point, from the supergrid's cell-centre nodes).  Missing data
    (land, and depths below WOA's own bottom) are filled LEVEL BY LEVEL
    from the nearest WOA cell with data at that same depth (breadth-first
    flood fill on the periodic WOA grid), so every model column, wet or
    dry, gets a complete profile and no surface water is carried to depth.

NetCDF-3 (classic / 64-bit offset) is read and written directly — every
input file of this configuration is classic, so no NetCDF library, numpy or
ncgen is needed.

    python3 tools/om1deg_prepare_inputs.py [--data DIR]

``--data`` defaults to ``$RDB_DATA_DIR/OM_1deg``.
"""

import argparse
import array
import collections
import os
import struct
import sys

MINIMUM_DEPTH = 9.5  # OM_1deg MOM_input
MAXIMUM_DEPTH = 6500.0  # OM_1deg MOM_input
WOA_FILL = 9.0e36  # WOA13 _FillValue is 9.96921e+36

# ----------------------------------------------------------------------------
# Minimal NetCDF-3 reader/writer
# ----------------------------------------------------------------------------
NC_DIMENSION, NC_VARIABLE, NC_ATTRIBUTE = 10, 11, 12
_TYPES = {1: ("b", 1), 2: ("c", 1), 3: ("h", 2), 4: ("i", 4), 5: ("f", 4), 6: ("d", 8)}


class NC3:
    """Read-only view of a classic / 64-bit-offset NetCDF file."""

    def __init__(self, path):
        self.f = open(path, "rb")
        magic = self.f.read(4)
        if magic[:3] != b"CDF" or magic[3] not in (1, 2):
            raise ValueError(f"{path}: not a NetCDF-3 classic/64-bit-offset file")
        self.offset_size = 8 if magic[3] == 2 else 4
        self.numrecs = self._int()
        self.dims = []  # [(name, length)], length 0 = record dim
        self._expect_list(NC_DIMENSION, self._read_dim, self.dims)
        self._skip_atts()
        self.vars = {}
        vlist = []
        self._expect_list(NC_VARIABLE, self._read_var, vlist)
        rec = [v for v in vlist if v["is_rec"]]
        self.recsize = sum(v["vsize"] for v in rec)
        if len(rec) == 1:  # single record variable: no padding between records
            v = rec[0]
            self.recsize = v["nelem"] * _TYPES[v["type"]][1]
        for v in vlist:
            self.vars[v["name"]] = v

    def _int(self):
        return struct.unpack(">i", self.f.read(4))[0]

    def _name(self):
        n = self._int()
        s = self.f.read(n).decode()
        self.f.read((4 - n % 4) % 4)
        return s

    def _expect_list(self, tag, reader, out):
        t, n = self._int(), self._int()
        if t == 0 and n == 0:
            return
        if t != tag:
            raise ValueError(f"bad NetCDF header tag {t} (expected {tag})")
        for _ in range(n):
            out.append(reader())

    def _read_dim(self):
        return (self._name(), self._int())

    def _skip_atts(self):
        t, n = self._int(), self._int()
        if t == 0 and n == 0:
            return
        for _ in range(n):
            self._name()
            typ, nel = self._int(), self._int()
            nb = nel * _TYPES[typ][1]
            self.f.read(nb + (4 - nb % 4) % 4)

    def _read_var(self):
        name = self._name()
        nd = self._int()
        dimids = [self._int() for _ in range(nd)]
        self._skip_atts()
        typ = self._int()
        vsize = self._int()
        begin = struct.unpack(">q" if self.offset_size == 8 else ">i",
                              self.f.read(self.offset_size))[0]
        shape = [self.dims[d][1] for d in dimids]
        is_rec = nd > 0 and self.dims[dimids[0]][1] == 0
        per_rec = shape[1:] if is_rec else shape
        nelem = 1
        for s in per_rec:
            nelem *= s
        return {"name": name, "shape": shape, "type": typ, "vsize": vsize,
                "begin": begin, "is_rec": is_rec, "nelem": nelem}

    def read(self, name, record=0):
        """Whole variable (or one record of a record variable), flat, C order."""
        v = self.vars[name]
        code, size = _TYPES[v["type"]]
        off = v["begin"] + (record * self.recsize if v["is_rec"] else 0)
        self.f.seek(off)
        a = array.array(code)
        a.frombytes(self.f.read(v["nelem"] * size))
        if sys.byteorder == "little":
            a.byteswap()
        return a

    def shape(self, name):
        return self.vars[name]["shape"]


def write_nc3(path, dims, variables):
    """Write a 64-bit-offset NetCDF file.

    dims: [(name, length)]; variables: [(name, [dim names], 'f'|'d',
    {attr: str}, flat C-order data)].  No record dimension.
    """
    dim_index = {d[0]: i for i, d in enumerate(dims)}

    def name_bytes(s):
        b = s.encode()
        return struct.pack(">i", len(b)) + b + b"\0" * ((4 - len(b) % 4) % 4)

    def atts_bytes(atts):
        if not atts:
            return struct.pack(">ii", 0, 0)
        out = struct.pack(">ii", NC_ATTRIBUTE, len(atts))
        for k, val in atts.items():
            b = val.encode()
            out += name_bytes(k) + struct.pack(">ii", 2, len(b)) + b + \
                b"\0" * ((4 - len(b) % 4) % 4)
        return out

    def header(begins):
        h = b"CDF\x02" + struct.pack(">i", 0)
        h += struct.pack(">ii", NC_DIMENSION, len(dims))
        for n, ln in dims:
            h += name_bytes(n) + struct.pack(">i", ln)
        h += struct.pack(">ii", 0, 0)  # no global attributes
        h += struct.pack(">ii", NC_VARIABLE, len(variables))
        for (n, dn, code, atts, data), b in zip(variables, begins):
            h += name_bytes(n) + struct.pack(">i", len(dn))
            h += b"".join(struct.pack(">i", dim_index[d]) for d in dn)
            h += atts_bytes(atts)
            typ = 5 if code == "f" else 6
            nb = len(data) * (4 if code == "f" else 8)
            vsize = nb + (4 - nb % 4) % 4
            h += struct.pack(">iiq", typ, min(vsize, 2**31 - 1), b)
        return h

    sizes = []
    for n, dn, code, atts, data in variables:
        nelem = 1
        for d in dn:
            nelem *= dims[dim_index[d]][1]
        if nelem != len(data):
            raise ValueError(f"{n}: {len(data)} values for shape {nelem}")
        nb = nelem * (4 if code == "f" else 8)
        sizes.append(nb + (4 - nb % 4) % 4)
    hlen = len(header([0] * len(variables)))
    begins, pos = [], hlen
    for s in sizes:
        begins.append(pos)
        pos += s
    tmp = path + ".part"
    with open(tmp, "wb") as f:
        f.write(header(begins))
        for (n, dn, code, atts, data), s in zip(variables, sizes):
            a = array.array(code, data)
            if sys.byteorder == "little":
                a.byteswap()
            b = a.tobytes()
            f.write(b + b"\0" * (s - len(b)))
    os.replace(tmp, path)


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
def model_tpoints(hgrid):
    """(lon, lat) of the model tracer points, j-major lists (ny, nx)."""
    nc = NC3(hgrid)
    nyp, nxp = nc.shape("x")
    x, y = nc.read("x"), nc.read("y")
    nx, ny = (nxp - 1) // 2, (nyp - 1) // 2
    lon = [[x[(2 * j + 1) * nxp + 2 * i + 1] for i in range(nx)] for j in range(ny)]
    lat = [[y[(2 * j + 1) * nxp + 2 * i + 1] for i in range(nx)] for j in range(ny)]
    return nx, ny, lon, lat


def make_bathymetry(topog, out, nx, ny):
    nc = NC3(topog)
    if nc.shape("depth") != [ny, nx]:
        raise SystemExit(f"topog depth shape {nc.shape('depth')} != {[ny, nx]}")
    d = nc.read("depth")
    out_d = []
    nwet = nraise = ncap = 0
    for v in d:
        if v > 0.0:
            nwet += 1
            if v < MINIMUM_DEPTH:
                nraise += 1
                v = MINIMUM_DEPTH
            if v > MAXIMUM_DEPTH:
                ncap += 1
                v = MAXIMUM_DEPTH
        else:
            v = 0.0
        out_d.append(float(v))
    write_nc3(out, [("y", ny), ("x", nx)],
              [("depth", ["y", "x"], "d",
                {"units": "m", "positive": "down",
                 "source": "MOM6 OM_1deg topog.nc; MINIMUM_DEPTH 9.5 / "
                           "MAXIMUM_DEPTH 6500 (MASKING_DEPTH 0)"}, out_d)])
    print(f"bathymetry: {out}\n  wet {nwet}, raised to {MINIMUM_DEPTH} m: {nraise}, "
          f"capped at {MAXIMUM_DEPTH} m: {ncap}")


def woa_filled_levels(path, var):
    """January `var` with every missing value filled, as one flat lat-major
    list of 64800 values per WOA depth level.

    Filled LEVEL BY LEVEL (MOM6's ``fill_miss_2d`` idea): at each depth a
    missing cell takes the value of the nearest cell that has data AT THAT
    SAME DEPTH (breadth-first flood on the periodic WOA grid).  Extending a
    column downward from its own deepest value instead is wrong wherever
    WOA's coastline/bathymetry is shallower than the model's: a WOA coastal
    column with data only at the surface would carry 29 degC water to the
    model's 3500 m bed (it did — and drove 6 m/s bed jets within 4 steps).
    A level with no data anywhere copies the level above.
    """
    nc = NC3(path)
    nt, nz, nlat, nlon = nc.shape(var)
    depth = list(nc.read("depth"))
    a = nc.read(var, record=0)  # (depth, lat, lon), January
    ncol = nlat * nlon
    levels = []
    for k in range(nz):
        lev = [a[k * ncol + c] for c in range(ncol)]
        have = [abs(v) < WOA_FILL for v in lev]
        if not any(have):
            levels.append(list(levels[-1]))
            continue
        q = collections.deque(c for c in range(ncol) if have[c])
        while q:
            c = q.popleft()
            jj, ii = divmod(c, nlon)
            for dj, di in ((0, 1), (0, -1), (1, 0), (-1, 0)):
                j2 = jj + dj
                if j2 < 0 or j2 >= nlat:
                    continue
                n = j2 * nlon + (ii + di) % nlon
                if not have[n]:
                    have[n] = True
                    lev[n] = lev[c]
                    q.append(n)
        levels.append(lev)
    lat0 = nc.read("lat")[0]
    lon0 = nc.read("lon")[0]
    return depth, nlat, nlon, lat0, lon0, levels


def make_ic(ptemp, salt, out, nx, ny, lon, lat):
    fields = {}
    for key, path, var in (("temp", ptemp, "ptemp_an"), ("salt", salt, "s_an")):
        depth, nlat, nlon, lat0, lon0, levels = woa_filled_levels(path, var)
        nz = len(depth)
        flat = [0.0] * (nz * ny * nx)
        for j in range(ny):
            for i in range(nx):
                jw = int((lat[j][i] - (lat0 - 0.5)) // 1.0)
                jw = min(max(jw, 0), nlat - 1)
                iw = int((lon[j][i] - (lon0 - 0.5)) // 1.0) % nlon
                c = jw * nlon + iw
                base = j * nx + i
                for k in range(nz):
                    flat[k * ny * nx + base] = levels[k][c]
        fields[key] = flat
        print(f"  {var}: min {min(flat):.3f} max {max(flat):.3f} on {nz} levels")
    write_nc3(out, [("z", nz), ("y", ny), ("x", nx)], [
        ("z_src", ["z"], "d", {"units": "m", "positive": "down"}, depth),
        ("temp", ["z", "y", "x"], "f",
         {"units": "degC", "long_name": "WOA13 decav January potential temperature"},
         fields["temp"]),
        ("salt", ["z", "y", "x"], "f",
         {"units": "psu", "long_name": "WOA13 decav January salinity"}, fields["salt"]),
    ])
    print(f"initial condition: {out}")


def main():
    root = os.environ.get("RDB_DATA_DIR")
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--data", default=None,
                    help="data directory (default: $RDB_DATA_DIR/OM_1deg)")
    args = ap.parse_args()
    if args.data is None and not root:
        ap.error("pass --data DIR or set RDB_DATA_DIR (see tools/fetch_om1deg.py)")
    data = os.path.abspath(args.data or os.path.join(root, "OM_1deg"))
    p = lambda n: os.path.join(data, n)  # noqa: E731
    nx, ny, lon, lat = model_tpoints(p("ocean_hgrid.nc"))
    print(f"model grid {nx} x {ny} from {p('ocean_hgrid.nc')}")
    make_bathymetry(p("topog.nc"), p("bathy_om1deg.nc"), nx, ny)
    make_ic(p("woa13_decav_ptemp_monthly_fulldepth_01.nc"),
            p("woa13_decav_s_monthly_fulldepth_01.nc"),
            p("ic_woa13_jan.nc"), nx, ny, lon, lat)


if __name__ == "__main__":
    main()
