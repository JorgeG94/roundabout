#!/usr/bin/env python3
"""Generate a synthetic time-varying surface-forcing NetCDF for
`&ocean_dataovr_nml` (PR-15).

Writes NetCDF-3 *classic* (CDF-1) by hand from the standard library only.
This box has no `netCDF4`, no `scipy` and no `numpy`, and the project
forbids installing packages -- so the format is emitted byte-for-byte
rather than through a library. CDF-1 is small and fully specified, which
makes this a ~150-line writer instead of a dependency.

Two deliberate simplifications:

  * **No unlimited dimension.** `time` is a FIXED dimension, so every
    variable lives in the non-record data section and the layout is a
    simple concatenation. The reader only needs `field(x, y, time)` plus
    a `time` coordinate variable; nothing requires `time` to be
    unlimited.
  * **Everything is NC_DOUBLE**, matching the solver's working precision,
    so no scale/offset games are needed to round-trip exactly.

Dimension order: netCDF stores dimension ids in C order (first varies
SLOWEST), while the netCDF-Fortran API takes them in Fortran order. A
variable that Fortran must see as `taux(x, y, time)` is therefore
written here with the dimid list `[time, y, x]` and data laid out with
x varying fastest. Get this backwards and the reader aborts with a
horizontal-dim mismatch (which is exactly what that check is for).

Usage:
    make_forcing_nc.py OUT.nc --nx 44 --ny 40 [--nt 4] [--period 8640000]
                              [--tau0 0.1]
"""

import argparse
import math
import struct

# --- CDF-1 tags -----------------------------------------------------------
NC_DIMENSION = 0x0A
NC_VARIABLE = 0x0B
NC_ATTRIBUTE = 0x0C
NC_CHAR = 2
NC_DOUBLE = 6
ABSENT = struct.pack(">II", 0, 0)


def _pad4(n):
    """Bytes of zero padding needed to reach a 4-byte boundary."""
    return (4 - (n % 4)) % 4


def _name(s):
    """CDF-1 `name`: length prefix, the bytes, then padding to 4."""
    b = s.encode("ascii")
    return struct.pack(">I", len(b)) + b + b"\0" * _pad4(len(b))


def _text_attr(name, value):
    b = value.encode("ascii")
    return (_name(name) + struct.pack(">II", NC_CHAR, len(b))
            + b + b"\0" * _pad4(len(b)))


def _doubles(values):
    return struct.pack(">%dd" % len(values), *values)


def build_cdf1(dims, variables):
    """Serialise a CDF-1 file.

    dims      : list of (name, length), in declaration order
    variables : list of dicts {name, dimids (C order), attrs, data}
                where `data` is a flat list of floats, laid out with the
                LAST dimid varying fastest.

    Offsets are resolved in two passes: the header cannot be written
    until every `begin` is known, and every `begin` depends on the total
    header size. Pass 1 emits the header with placeholder offsets purely
    to measure it; pass 2 rewrites it with the real ones.
    """
    dim_list = (struct.pack(">II", NC_DIMENSION, len(dims))
                + b"".join(_name(n) + struct.pack(">I", ln) for n, ln in dims))

    def var_entries(offsets):
        out = struct.pack(">II", NC_VARIABLE, len(variables))
        for var, begin in zip(variables, offsets):
            out += _name(var["name"])
            out += struct.pack(">I", len(var["dimids"]))
            out += b"".join(struct.pack(">I", d) for d in var["dimids"])
            attrs = var.get("attrs", {})
            if attrs:
                out += struct.pack(">II", NC_ATTRIBUTE, len(attrs))
                for k, v in attrs.items():
                    out += _text_attr(k, v)
            else:
                out += ABSENT
            nbytes = len(var["data"]) * 8
            out += struct.pack(">III", NC_DOUBLE, nbytes, begin)
        return out

    placeholder = [0] * len(variables)
    header_len = (4 + 4 + len(dim_list) + len(ABSENT)
                  + len(var_entries(placeholder)))

    offsets, cursor = [], header_len
    for var in variables:
        offsets.append(cursor)
        nbytes = len(var["data"]) * 8
        cursor += nbytes + _pad4(nbytes)

    header = (b"CDF\x01" + struct.pack(">I", 0)      # magic, numrecs = 0
              + dim_list + ABSENT                    # no global attributes
              + var_entries(offsets))
    assert len(header) == header_len, "header size drifted between passes"

    body = b""
    for var in variables:
        blob = _doubles(var["data"])
        body += blob + b"\0" * _pad4(len(blob))
    return header + body


def double_gyre_wind(nx, ny, nt, tau0, period, steady):
    """The MOM6 double-gyre zonal wind, optionally with a seasonal cycle.

    The spatial profile is a byte-for-byte match of what
    `wind_config="2gyre"` builds internally
    (`ocean_surfstress_set_2gyre`):

        tau_x(j) = tau0 * (1 - cos(2*pi * (j_phys - 0.5) / ny))

    -- note `(1 - cos)`, which is non-negative and peaks at 2*tau0, NOT
    a bare `-cos`. Matching it exactly is what lets `--steady` be a
    checkable claim rather than a vibe: a run driven from a steady file
    must reproduce the analytic run.

    `steady=True` writes that profile unchanged at every record, so the
    time interpolation is exercised but returns the same field for any
    weight. `steady=False` modulates it by `sin(2*pi*t/period)`, which
    makes the interpolation visible in the data (and reverses the wind
    in the second half of the cycle).

    Returns (t_axis, taux_flat, tauy_flat), the tau arrays flattened in
    C order [time][y][x] -- x fastest.
    """
    t_axis = [period * k / float(nt) for k in range(nt)]
    taux, tauy = [], []
    for k in range(nt):
        season = 1.0 if steady else math.sin(2.0 * math.pi * t_axis[k] / period)
        # taux is (time, y, nx+1): cell-row profile, one extra face in x.
        for j in range(ny):
            y_rel = (j + 0.5) / float(ny)
            prof = tau0 * (1.0 - math.cos(2.0 * math.pi * y_rel)) * season
            for _ in range(nx + 1):
                taux.append(prof)
        # tauy is (time, ny+1, nx): identically zero in this profile.
        for _ in range(ny + 1):
            for _ in range(nx):
                tauy.append(0.0)
    return t_axis, taux, tauy


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="output .nc path")
    ap.add_argument("--nx", type=int, required=True,
                    help="physical cells in x (must match the model grid)")
    ap.add_argument("--ny", type=int, required=True,
                    help="physical cells in y (must match the model grid)")
    ap.add_argument("--nt", type=int, default=4, help="time records")
    ap.add_argument("--period", type=float, default=8640000.0,
                    help="seasonal period in seconds (default 100 days)")
    ap.add_argument("--tau0", type=float, default=0.1,
                    help="wind stress magnitude, Pa (profile peaks at 2*tau0)")
    ap.add_argument("--steady", action="store_true",
                    help="write the analytic 2gyre profile at every record "
                         "(no seasonal modulation) -- the mode whose run must "
                         "reproduce an analytic wind_config=\"2gyre\" run")
    args = ap.parse_args()

    t_axis, taux, tauy = double_gyre_wind(args.nx, args.ny, args.nt,
                                          args.tau0, args.period, args.steady)

    # C-grid staggering: an x-face array spans nx+1 faces, a y-face array
    # ny+1. The reader registers tau_x with nx_extra=1 and tau_y with
    # ny_extra=1, so the file must carry those extra rows/columns -- the
    # trailing face is owned by the local rank and no halo exchange can
    # supply it. Hence two x dims and two y dims in one file.
    # dimids index into this list; C order puts time first.
    dims = [("time", args.nt), ("y", args.ny), ("x", args.nx),
            ("xf", args.nx + 1), ("yf", args.ny + 1)]
    variables = [
        {"name": "time", "dimids": [0],
         "attrs": {"units": "seconds"}, "data": t_axis},
        {"name": "taux", "dimids": [0, 1, 3],       # (time, y, xf)
         "attrs": {"units": "Pa"}, "data": taux},
        {"name": "tauy", "dimids": [0, 4, 2],       # (time, yf, x)
         "attrs": {"units": "Pa"}, "data": tauy},
    ]

    with open(args.out, "wb") as fh:
        fh.write(build_cdf1(dims, variables))

    print("wrote %s: %dx%d, %d records, tau0 %.6g Pa, %s"
          % (args.out, args.nx, args.ny, args.nt, args.tau0,
             "steady (analytic 2gyre)" if args.steady
             else "seasonal, period %.6g s" % args.period))


if __name__ == "__main__":
    main()
