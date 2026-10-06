"""The compatibility matrix's legs after the first run: checks 4-6 and the
energy bound, plus the readers and comparators they share.

`compat_matrix.py` decides WHICH cells exist (the covering array, fixed by
checks 1-3 alone, so the cell list is the same on every backend and every
leg) and runs each accepted cell once.  This module holds what happens to a
cell that PASSES those three checks:

* ENERGY -- the run's own kinetic-energy series against a bound derived
  from the PASS population (`energy_check`, constants below).
* RESTART (check 5) -- 24 steps straight vs 12 + checkpoint + warm restart
  through the production engine path + 12, every restart-registry field of
  the two final checkpoints BITWISE, ghosts included (`compare_full`).
* DECOMP (check 4) -- the 1-rank run vs 2x2 and 4x1 MPI decompositions, the
  OWNED window of every restart-registry field BITWISE, with the exemptions
  of `tests/mpi/test_ocean_decomp_bitid_mpi.F90` read from that file
  (`decomp_exempt_tags`) so the list lives in ONE place
  (`compare_decomp`).
* CROSS_BACKEND (check 6) -- a GPU (nvfortran) run against the CPU
  (gfortran) record of the same cell: the L2 and max norms of every
  restart-registry field of the final checkpoint (owned cells, full
  precision) and the per-step console En / MaxCFL series, each inside a
  relative band that widens with the step count, and the GPU budgets inside
  the round-off envelope (`cross_backend_check`, constants below).

The restart files are NetCDF-4 (HDF5).  The standard library cannot read
HDF5, so `read_nc` converts with `nccopy -k '64-bit offset'` (shipped with
netcdf-c, next to `ncdump`) and parses the classic format itself; bitwise
means the raw IEEE bytes are compared, never decoded values.

Stdlib only -- never `pip install` anything for this.
"""

import fnmatch
import math
import os
import re
import shutil
import struct
import subprocess

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))
DECOMP_TEST = os.path.join(REPO_ROOT, "tests", "mpi", "test_ocean_decomp_bitid_mpi.F90")

# ===========================================================================
# Band and bound constants (chosen from the data; README "Checks 4-6")
# ===========================================================================
# Modelled on validation_examples/ocean/global_1deg/check_against_reference.py
# (a relative En / MaxCFL band that widens with time, budgets in the
# round-off envelope), with the numbers from the data (README, "The
# cross-backend band"; gfortran 15.1 vs nvfortran 26.5 cc70, V100,
# 2026-10-04, 63 cells passing checks 1-3 on both):
#
# CROSS_BACKEND, the real test -- the L2 and max norm of EVERY restart-
# registry field of the final checkpoint, owned cells, at full precision:
# |gpu - cpu| / |cpu| <= FIELD_BAND_FLOOR + FIELD_BAND_RATE * n  (1.06e-9 at
# n = 24).  Measured: median 1e-13, PASS population max 5.7e-11
# (kappa-shear kd), the z_fixed open-step XFAIL cells 1.1e-10 -- 18x margin
# over the PASS maximum.  The out-of-band cells sat at 5.9e-9 .. 4e26.
FIELD_BAND_FLOOR = 1.0e-10
FIELD_BAND_RATE = 4.0e-11
# ... and the console En / MaxCFL every step:  EN_BAND_FLOOR + EN_BAND_RATE
# * n (MaxCFL: CFL_FACTOR times that).  The console prints 4 significant
# digits, so this is a coarse guard (print rounding alone is up to 1e-3
# between two values); measured: every step of every passing cell agreed to
# the last printed digit.
EN_BAND_FLOOR = 1.0e-3
EN_BAND_RATE = 1.0e-4
CFL_FACTOR = 2.0
# ... and every budget residual of the GPU run must sit inside the round-off
# envelope  |Error| <= BUDGET_RO_FLOOR + BUDGET_RO_RATE * n  (3.4e-13 at
# n = 24) -- a level the CPU run need not match digit for digit (summation
# order differs) but a leak exceeds by orders of magnitude.  Measured: GPU
# PASS-cell worst Mass 1.7e-15, Heat 2.1e-15, Salt 7.9e-16 (CPU: 1.6e-15,
# 1.8e-15, 7.9e-16).
BUDGET_RO_FLOOR = 1.0e-13
BUDGET_RO_RATE = 1.0e-14

# ENERGY -- every cell starts from rest under the same 0.1 Pa wind, so the
# kinetic energy it holds after 24 steps is set by the wind and the
# geometry, and the closures move it only by a factor of about two: the
# PASS population of each geometry sits within 0.5-1.6x of its median
# (README, "The energy bound").  A cell holding ENERGY_RATIO_MAX times its
# geometry's reference is being driven by something other than the wind --
# the z_fixed open-staircase PGF sits at 5-61x.  EN_REF is that median
# (m2/s2, En at step 24, gfortran, measured 2026-10-04, z_fixed_open cells
# excluded); the tripolar ring is 15 x 1 degree at 59-70 N, a different flow.
# "cliff" (2026-10-05, with MOM6's BBL glue the default): the median of the
# `domain` sweep's PASS population on the cliff geometry (every coordinate x
# both splits x both grids, base closures) -- 7.8e-3 against that sweep's
# closed-geometry 6.5e-3; the pairwise PASS population there is only 4 cells.
EN_REF = {"closed": 6.9e-3, "channel": 7.9e-3, "obc": 7.6e-3, "tripolar": 1.7e-4,
          "cavity": 9.5e-3, "cliff": 7.8e-3}
ENERGY_RATIO_MAX = 2.5
# ... and its growth over the last third of the run: from rest a wind
# spin-up grows En like t^2 at most (En(24)/En(16) <= 2.25) and the PASS
# population tops out at 1.8 as the closures engage; an instability that
# sets in late accelerates past it.
GROWTH_WINDOW = (16, 24)
GROWTH_MAX = 3.0


# ===========================================================================
# NetCDF: the classic format, read raw (bitwise comparison)
# ===========================================================================
# Scalar bookkeeping variables in a checkpoint: never compared as fields.
_META_VARS = ("time", "n_steps", "outer_step_count")
_NC_SIZES = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}
_HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"


class NcVar(object):
    __slots__ = ("name", "dims", "shape", "nc_type", "raw")

    def __init__(self, name, dims, shape, nc_type, raw):
        self.name, self.dims, self.shape, self.nc_type, self.raw = name, dims, shape, nc_type, raw

    @property
    def itemsize(self):
        return _NC_SIZES[self.nc_type]

    def value(self, flat):
        """Decoded value at a flat C-order index (for reports only)."""
        sz = self.itemsize
        fmt = {1: ">b", 2: ">c", 3: ">h", 4: ">i", 5: ">f", 6: ">d"}[self.nc_type]
        return struct.unpack(fmt, self.raw[flat * sz:(flat + 1) * sz])[0]


def nccopy_binary():
    exe = shutil.which("nccopy")
    if not exe:
        raise RuntimeError("nccopy not found on PATH: the restart / decomposition legs read "
                           "NetCDF-4 checkpoints through it (it ships with netcdf-c)")
    return exe


def to_classic(path, out):
    """`nccopy -k '64-bit offset'` path -> out (no-op copy for a classic file)."""
    with open(path, "rb") as fh:
        head = fh.read(8)
    if head[:3] == b"CDF":
        shutil.copyfile(path, out)
        return out
    if head != _HDF5_MAGIC:
        raise ValueError("{}: neither classic NetCDF nor HDF5".format(path))
    p = subprocess.run([nccopy_binary(), "-k", "64-bit offset", path, out],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError("nccopy {} failed: {}".format(path, p.stderr.decode("utf-8", "replace")))
    return out


def parse_classic(data):
    """Parse a CDF-1 / CDF-2 file held in `data` -> (gatts, {name: NcVar})."""
    if data[:3] != b"CDF" or data[3] not in (1, 2):
        raise ValueError("not a CDF-1/CDF-2 file")
    off_fmt = ">i" if data[3] == 1 else ">q"
    pos = [4]

    def i32():
        v = struct.unpack(">i", data[pos[0]:pos[0] + 4])[0]
        pos[0] += 4
        return v

    def name():
        n = i32()
        s = data[pos[0]:pos[0] + n].decode("ascii")
        pos[0] += n + (4 - n % 4) % 4
        return s

    def values(t, n):
        sz = _NC_SIZES[t] * n
        raw = data[pos[0]:pos[0] + sz]
        pos[0] += sz + (4 - sz % 4) % 4
        if t == 2:
            return raw.decode("ascii", "replace").rstrip("\0")
        fmt = {1: "b", 3: "h", 4: "i", 5: "f", 6: "d"}[t]
        vals = struct.unpack(">{}{}".format(n, fmt), raw)
        return vals[0] if n == 1 else list(vals)

    def atts():
        tag, n = i32(), i32()
        out = {}
        for _ in range(n if tag else 0):
            k = name()
            t = i32()
            out[k] = values(t, i32())
        return out

    i32()                                   # numrecs
    tag, n = i32(), i32()
    dims = []
    for _ in range(n if tag else 0):
        dims.append((name(), i32()))
    gatts = atts()
    tag, n = i32(), i32()
    variables = {}
    for _ in range(n if tag else 0):
        vname = name()
        ids = [i32() for _ in range(i32())]
        atts()
        t = i32()
        i32()                               # vsize
        begin = struct.unpack(off_fmt, data[pos[0]:pos[0] + struct.calcsize(off_fmt)])[0]
        pos[0] += struct.calcsize(off_fmt)
        shape = tuple(dims[d][1] for d in ids)
        if any(s == 0 for s in shape):
            raise ValueError("record variable {}: not supported".format(vname))
        count = 1
        for s in shape:
            count *= s
        variables[vname] = NcVar(vname, tuple(dims[d][0] for d in ids), shape, t,
                                 data[begin:begin + count * _NC_SIZES[t]])
    return gatts, variables


def read_nc(path, scratch):
    """(global attributes, {name: NcVar}) of a NetCDF-4 or classic file."""
    os.makedirs(scratch, exist_ok=True)
    out = os.path.join(scratch, os.path.basename(path) + ".cdf2")
    to_classic(path, out)
    with open(out, "rb") as fh:
        data = fh.read()
    os.remove(out)
    return parse_classic(data)


# ===========================================================================
# Check 5: RESTART -- two checkpoints of the same decomposition, bitwise
# ===========================================================================
def _first_diff(a, b, sz):
    """Flat index of the first differing item of two equal-length buffers,
    and the number of differing items."""
    n = len(a) // sz
    first, count = None, 0
    step = 512 * sz
    for blk in range(0, len(a), step):
        if a[blk:blk + step] == b[blk:blk + step]:
            continue
        for k in range(blk // sz, min(n, (blk + step) // sz)):
            if a[k * sz:(k + 1) * sz] != b[k * sz:(k + 1) * sz]:
                count += 1
                if first is None:
                    first = k
    return first, count


def _unflat(shape, flat):
    """C-order flat index -> Fortran (i, j[, k]) 1-based (the file stores
    a Fortran (x, y[, z]) array, so the C dims are reversed)."""
    idx = []
    for s in reversed(shape):
        idx.append(flat % s + 1)
        flat //= s
    return tuple(idx)


def compare_full(path_a, path_b, scratch):
    """Every variable of two checkpoints, bitwise (full local arrays: owned
    cells AND ghosts, plus time / n_steps).  -> list of difference strings."""
    _, va = read_nc(path_a, scratch)
    _, vb = read_nc(path_b, scratch)
    out = []
    if sorted(va) != sorted(vb):
        out.append("variable sets differ: only in A {}, only in B {}".format(
            sorted(set(va) - set(vb)), sorted(set(vb) - set(va))))
    for name in sorted(set(va) & set(vb)):
        a, b = va[name], vb[name]
        if a.shape != b.shape:
            out.append("{}: shape {} vs {}".format(name, a.shape, b.shape))
            continue
        if a.raw == b.raw:
            continue
        first, count = _first_diff(a.raw, b.raw, a.itemsize)
        out.append("{}: {} value(s) differ, first at (i,j,k)={} {!r} vs {!r}".format(
            name, count, _unflat(a.shape, first), a.value(first), b.value(first)))
    return out


def field_norms(path, scratch):
    """{field: [L2 norm, max |x|]} of every array field of a checkpoint, at
    full precision -- what the CROSS_BACKEND leg compares (the console
    prints En / MaxCFL to 4 digits only, too coarse to see two toolchains
    differ in 24 steps)."""
    ga, va = read_nc(path, scratch)
    ng, nxl, nyl = int(ga["nghost"]), int(ga["nx_local"]), int(ga["ny_local"])
    out = {}
    for name, v in sorted(va.items()):
        if not v.shape or v.nc_type != 6 or name in _META_VARS:
            continue
        allv = struct.unpack(">{}d".format(len(v.raw) // 8), v.raw)
        # The OWNED window (a staggered field owns both edge faces), as the
        # DECOMP leg compares it: ghost cells are not answers -- a ghost
        # nothing reads may legitimately hold anything on a device build.
        sh = ((1,) + v.shape) if len(v.shape) == 2 else v.shape
        nk, nya, nxa = sh
        ex, ey = nxa - (nxl + 2 * ng), nya - (nyl + 2 * ng)
        vals = [allv[(k * nya + j) * nxa + i] for k in range(nk)
                for j in range(ng, ng + nyl + ey) for i in range(ng, ng + nxl + ex)]
        fin = [x for x in vals if _finite(x)]
        if len(fin) != len(vals):
            out[name] = ["nan", "nan"]
            continue
        out[name] = [math.sqrt(math.fsum(x * x for x in fin)), max(abs(x) for x in fin)]
    return out


# ===========================================================================
# Check 4: DECOMP -- 1 rank vs a px x py decomposition, owned cells bitwise
# ===========================================================================
_EXEMPT_RE = re.compile(r"function\s+carried_tendency\s*\(\s*tag\s*\)(?P<body>.*?)"
                        r"end\s+function\s+carried_tendency", re.S | re.I)


def decomp_exempt_tags(path=DECOMP_TEST):
    """The restart-registry tags `test_ocean_decomp_bitid_mpi` does not
    compare across decompositions, read from ITS `carried_tendency` -- the
    one place the list (and the reason for it) lives."""
    with open(path) as fh:
        m = _EXEMPT_RE.search(fh.read())
    if not m:
        raise RuntimeError("{}: function carried_tendency not found -- the decomposition "
                           "exemptions moved; update compat_legs.decomp_exempt_tags".format(path))
    code = "\n".join(ln.split("!", 1)[0] for ln in m.group("body").splitlines())
    tags = re.findall(r"trim\s*\(\s*tag\s*\)\s*==\s*[\"']([^\"']+)[\"']", code)
    if not tags:
        raise RuntimeError("{}: carried_tendency names no tag".format(path))
    return tuple(tags)


def compare_decomp(ref_path, rank_paths, scratch, exempt):
    """Each rank's OWNED window of every array field against the matching
    window of the 1-rank checkpoint, bitwise; a staggered face array owns
    both of its tile's edge faces, so a seam face is checked on both ranks
    (the `test_ocean_decomp_bitid_mpi::compare` rule).  Rank-0 scalars and
    the time / step metadata are not compared (rank-local / global).
    -> (difference strings, number of fields compared, the worst
    normwise relative difference: max |diff| / max |field|)."""
    _, ref = read_nc(ref_path, scratch)
    out, nfields, worst = [], 0, 0.0
    for rank, path in enumerate(rank_paths):
        ga, dec = read_nc(path, scratch)
        ng = int(ga["nghost"])
        io, jo = int(ga["i_start"]) - 1, int(ga["j_start"]) - 1
        nxl, nyl = int(ga["nx_local"]), int(ga["ny_local"])
        names = [n for n in dec if dec[n].shape and n not in _META_VARS]
        if sorted(names) != sorted(n for n in ref if ref[n].shape and n not in _META_VARS):
            out.append("rank {}: field set differs from the 1-rank checkpoint".format(rank))
            continue
        for name in sorted(names):
            if name in exempt:
                continue
            a, b = dec[name], ref[name]
            sz = a.itemsize
            if len(a.shape) == 2:
                ash, bsh = (1,) + a.shape, (1,) + b.shape
            else:
                ash, bsh = a.shape, b.shape
            nk, nya, nxa = ash
            _, nyb, nxb = bsh
            ex, ey = nxa - (nxl + 2 * ng), nya - (nyl + 2 * ng)
            i0, i1 = ng, ng + nxl + ex          # 0-based, half-open
            j0, j1 = ng, ng + nyl + ey
            count, first, dmax = 0, None, 0.0
            scale = None
            for k in range(nk):
                for j in range(j0, j1):
                    ra = (k * nya + j) * nxa
                    rb = (k * nyb + j + jo) * nxb + io
                    sa = a.raw[(ra + i0) * sz:(ra + i1) * sz]
                    sb = b.raw[(rb + i0) * sz:(rb + i1) * sz]
                    if sa == sb:
                        continue
                    for i in range(i0, i1):
                        if sa[(i - i0) * sz:(i - i0 + 1) * sz] != sb[(i - i0) * sz:(i - i0 + 1) * sz]:
                            count += 1
                            va_, vb_ = a.value(ra + i), b.value(rb + i)
                            if first is None:
                                first = (i + 1, j + 1, k + 1, va_, vb_)
                            dmax = max(dmax, abs(va_ - vb_)) if _finite(va_ - vb_) \
                                else float("inf")
            if rank == 0:
                nfields += 1
            if count:
                if scale is None:
                    vals = struct.unpack(">{}d".format(len(b.raw) // 8), b.raw)
                    scale = max([abs(x) for x in vals if _finite(x)] or [0.0])
                rel = dmax / scale if scale > 0.0 else float("inf")
                worst = max(worst, rel)
                out.append("rank {} {}: {} owned value(s) differ (max |diff| / max |field| "
                           "{:.1e}), first at local (i,j,k)=({},{},{}) {!r} vs {!r}".format(
                               rank, name, count, rel, *first))
    return out, nfields, worst


# ===========================================================================
# ENERGY
# ===========================================================================
def energy_check(series, geometry):
    """`series` = [(step, En, MaxCFL), ...] of a run that passed checks 1-3,
    on the `geometry` axis value.  -> (ok, detail)."""
    en = {s: e for s, e, _ in series if e is not None}
    if not en:
        return False, "no En series"
    if geometry not in EN_REF:
        raise RuntimeError("compat_legs.EN_REF has no energy reference for geometry {!r}: "
                           "measure the PASS-population median En(24) for it and add it to "
                           "EN_REF".format(geometry))
    ref = EN_REF[geometry]
    last = max(en)
    e_last = en[last]
    ratio = e_last / ref
    lo, hi = GROWTH_WINDOW
    growth = (en[hi] / en[lo]) if (lo in en and hi in en and en[lo] > 0.0) else None
    detail = "En({}) {:.3e} = {:.2f} x EN_REF[{}]".format(last, e_last, ratio, geometry)
    if growth is not None:
        detail += ", En({})/En({}) {:.2f}".format(hi, lo, growth)
    if ratio > ENERGY_RATIO_MAX:
        return False, "En({}) {:.3e} is {:.1f} x the {} PASS-population reference {:.1e} " \
                      "(bound {:.1f}x)".format(last, e_last, ratio, geometry, ref, ENERGY_RATIO_MAX)
    if growth is not None and growth > GROWTH_MAX:
        return False, "En grew {:.2f}x over steps {}..{} (bound {:.1f}x)".format(
            growth, lo, hi, GROWTH_MAX)
    return True, detail


# ===========================================================================
# Check 6: CROSS_BACKEND
# ===========================================================================
def en_band(step):
    return EN_BAND_FLOOR + EN_BAND_RATE * step


def field_band(step):
    return FIELD_BAND_FLOOR + FIELD_BAND_RATE * step


def _rel(a, b):
    return abs(a - b) / abs(b) if b else abs(a - b)


def cross_backend_check(here, ref):
    """`here` / `ref`: the `run.metrics` of the same cell on this backend and
    on the reference backend.  -> (ok, detail, worst {'En': (rel, band, step), ...})."""
    sh = {s: (e, c) for s, e, c in here.get("series", [])}
    sr = {s: (e, c) for s, e, c in ref.get("series", [])}
    common = sorted(set(sh) & set(sr))
    if not common:
        return False, "no common [stats] step with the reference", {}
    worst = {"En": (0.0, 0.0, 0), "MaxCFL": (0.0, 0.0, 0)}
    bad = []
    for s in common:
        for k, (a, b), fac in (("En", (sh[s][0], sr[s][0]), 1.0),
                               ("MaxCFL", (sh[s][1], sr[s][1]), CFL_FACTOR)):
            if a is None or b is None:
                continue
            if not (_finite(a) and _finite(b)):
                bad.append("{} non-finite at step {}".format(k, s))
                continue
            r, band = _rel(a, b), fac * en_band(s)
            if r / band > worst[k][0] / max(worst[k][1], 1e-300):
                worst[k] = (r, band, s)
            if r > band and not bad:
                bad.append("{} at step {}: {:.4e} vs reference {:.4e} (rel {:.2e} > band "
                           "{:.2e})".format(k, s, a, b, r, band))
    nsteps = max(common)
    # The field norms of the final checkpoint, at full precision.
    fh, fr = here.get("field_norms"), ref.get("field_norms")
    worst["field"] = (0.0, field_band(nsteps), "")
    if fh is None or fr is None:
        bad.append("no checkpoint field norms on {} (run both legs with a checkpoint)".format(
            "this backend" if fh is None else "the reference backend"))
    else:
        for name in sorted(set(fh) & set(fr)):
            for k, label in ((0, "L2"), (1, "max")):
                a, b = fh[name][k], fr[name][k]
                if not (_finite(a) and _finite(b)):
                    bad.append("{} {} non-finite ({} vs {})".format(name, label, a, b))
                    continue
                if b == 0.0 and a == 0.0:
                    continue
                r = _rel(a, b)
                if r > worst["field"][0]:
                    worst["field"] = (r, field_band(nsteps), "{} {}".format(name, label))
                if r > field_band(nsteps):
                    bad.append("{} {} norm {:.16e} vs reference {:.16e} (rel {:.2e} > band "
                               "{:.2e})".format(name, label, a, b, r, field_band(nsteps)))
        if sorted(fh) != sorted(fr):
            bad.append("checkpoint field sets differ")
    for what, err in sorted((here.get("budget_worst") or {}).items()):
        lim = BUDGET_RO_FLOOR + BUDGET_RO_RATE * nsteps
        if not _finite(err) or abs(err) > lim:
            bad.append("{} residual {} outside the round-off envelope {:.1e}".format(what, err, lim))
    if len(common) < len(sr):
        bad.append("ran {} of the reference's {} steps".format(len(common), len(sr)))
    detail = ("worst field-norm rel {:.2e} ({}; band {:.2e}); console En rel {:.2e} (band "
              "{:.2e}), MaxCFL rel {:.2e}").format(
        worst["field"][0], worst["field"][2] or "-", worst["field"][1], worst["En"][0],
        worst["En"][1], worst["MaxCFL"][0])
    return (not bad), ("; ".join(bad) if bad else detail), worst


def _finite(x):
    return isinstance(x, (int, float)) and not math.isnan(x) and not math.isinf(x)


# ===========================================================================
# The per-PR slice: which axis values does a diff touch?
# ===========================================================================
def touched_values(changed, value_paths, always_full):
    """`changed`: repo-relative paths a diff touches.  `value_paths`:
    {(axis, value): [glob, ...]}.  -> ("full", reason) | ("slice", set of
    (axis, value)) | ("none", reason).

    A changed Fortran source no value claims is SHARED code (the dynamical
    core, the remap driver, the config) and can break any cell: full run.
    So does any path matching `always_full` (the matrix itself)."""
    hits, unclaimed = set(), []
    for f in changed:
        if any(fnmatch.fnmatch(f, g) for g in always_full):
            return "full", "{} is the matrix itself".format(f)
        mine = {av for av, globs in value_paths.items() if any(fnmatch.fnmatch(f, g) for g in globs)}
        if mine:
            hits |= mine
        elif re.match(r"(src|app)/.*\.(F90|f90|inc)$", f) or f in ("CMakeLists.txt",) \
                or f.startswith("cmake/"):
            unclaimed.append(f)
    if unclaimed:
        return "full", "shared source touched: {}{}".format(
            ", ".join(unclaimed[:4]), " ..." if len(unclaimed) > 4 else "")
    if not hits:
        return "none", "no model source touched"
    return "slice", hits


# ===========================================================================
# Who tests the legs (no model; run by `compat_matrix.py self-test`)
# ===========================================================================
def self_tests(cm, scratch):
    """-> [(ok, what)].  `cm` is the compat_matrix module (its classic
    NetCDF writer builds the synthetic checkpoints)."""
    out = []
    os.makedirs(scratch, exist_ok=True)
    # --- the exemptions are read from the MPI test, not copied ---
    try:
        tags = decomp_exempt_tags()
        out.append(("hvisc_du_visc" in tags and "hvisc_dv_visc" in tags,
                    "decomp exemptions read from test_ocean_decomp_bitid_mpi {}".format(tags)))
    except RuntimeError as exc:
        out.append((False, str(exc)))
    # --- a synthetic 6 x 4 domain, ng = 1, split 2 x 1 ---
    nx, ny, ng = 6, 4, 1

    def field(nxa, nya, base):          # C-order (y, x), full extent incl. ghosts
        return [base + 100.0 * j + i for j in range(nya) for i in range(nxa)]

    gx, gy = nx + 2 * ng, ny + 2 * ng
    ref = {"c": (gx, gy, field(gx, gy, 0.0)), "u": (gx + 1, gy, field(gx + 1, gy, 0.5)),
           "hvisc_du_visc": (gx + 1, gy, field(gx + 1, gy, 0.25))}

    def write(path, flds, i_start, nxl):
        dims, variables = [], []
        for n, (nxa, nya, vals) in sorted(flds.items()):
            dims += [("x_" + n, nxa), ("y_" + n, nya)]
            variables.append((n, ("y_" + n, "x_" + n), vals))
        cm.write_netcdf_classic(path, dims, variables, gatts=[
            ("i_start", i_start), ("j_start", 1), ("nx_local", nxl), ("ny_local", ny),
            ("nghost", ng)])

    write(os.path.join(scratch, "ref.nc"), ref, 1, nx)

    def tile(i_start, nxl, poke=None):
        flds = {}
        for n, (nxa, nya, vals) in ref.items():
            w = nxl + 2 * ng + (nxa - gx)
            loc = [vals[j * nxa + (i_start - 1) + i] for j in range(nya) for i in range(w)]
            if poke and poke[0] == n:
                i, j = poke[1]
                loc[j * w + i] += 1.0
            flds[n] = (w, nya, loc)
        return flds

    def run(pokes):
        paths = []
        for r, (i0, nxl) in enumerate(((1, 3), (4, 3))):
            p = os.path.join(scratch, "rank{}.nc".format(r))
            write(p, tile(i0, nxl, pokes.get(r)), i0, nxl)
            paths.append(p)
        return compare_decomp(os.path.join(scratch, "ref.nc"), paths, scratch,
                              ("hvisc_du_visc",))[:2]

    d, nf = run({})
    out.append((d == [] and nf == 2, "decomp: identical tiles compare clean ({} fields)".format(nf)))
    d, _ = run({1: ("c", (0, 2))})               # rank 1's WEST GHOST column
    out.append((d == [], "decomp: a stale ghost is not an owned mismatch {}".format(d)))
    d, _ = run({1: ("c", (2, 2))})               # owned: local (3, 3) 1-based
    out.append((len(d) == 1 and "rank 1 c: 1 owned" in d[0] and "(3,3,1)" in d[0],
                "decomp: an owned mismatch is found and located {}".format(d)))
    d, _ = run({0: ("u", (4, 1))})               # rank 0's EAST seam face (owned)
    out.append((len(d) == 1 and "rank 0 u" in d[0],
                "decomp: a seam face is checked on the rank west of it {}".format(d)))
    d, _ = run({0: ("hvisc_du_visc", (2, 2))})
    out.append((d == [], "decomp: the exempt carried tendency is skipped {}".format(d)))
    # --- RESTART: full arrays, ghosts included ---
    write(os.path.join(scratch, "a.nc"), ref, 1, nx)
    poked = dict(ref)
    vals = list(ref["c"][2])
    vals[0] += 1.0                               # a corner GHOST
    poked["c"] = (gx, gy, vals)
    write(os.path.join(scratch, "b.nc"), poked, 1, nx)
    d0 = compare_full(os.path.join(scratch, "a.nc"), os.path.join(scratch, "a.nc"), scratch)
    d1 = compare_full(os.path.join(scratch, "a.nc"), os.path.join(scratch, "b.nc"), scratch)
    out.append((d0 == [] and len(d1) == 1 and "(1, 1)" in d1[0],
                "restart: ghosts are compared too {}".format(d1)))
    # --- ENERGY ---
    ref = EN_REF["closed"]
    spin_up = [(s, ref * (s / 24.0) ** 2, 0.01) for s in range(1, 25)]
    hot = [(s, 2.0 * ENERGY_RATIO_MAX * ref * (s / 24.0) ** 2, 0.01) for s in range(1, 25)]
    burst = [(s, 0.2 * ref * (1.0 if s <= 16 else 4.0 ** ((s - 16) / 8.0)), 0.01)
             for s in range(1, 25)]
    out.append((energy_check(spin_up, "closed")[0] and not energy_check(hot, "closed")[0]
                and not energy_check(burst, "closed")[0]
                and not energy_check(spin_up, "tripolar")[0],
                "energy: a wind spin-up passes; 2x the bound, a 4x late burst and a basin's "
                "energy on the tripolar ring fail"))
    # --- CROSS_BACKEND ---
    fb = field_band(24)
    ref_m = {"series": [[s, 1.0e-3 * s, 0.01 * s] for s in range(1, 25)], "budget_worst": {},
             "field_norms": {"h": [100.0, 2.0], "u": [3.0, 0.1]}}
    near = {"series": [[s, 1.0e-3 * s * (1 + 0.5 * en_band(s)), 0.01 * s] for s in range(1, 25)],
            "budget_worst": {"Mass": 1e-15},
            "field_norms": {"h": [100.0 * (1 + 0.5 * fb), 2.0], "u": [3.0, 0.1]}}
    far = dict(near, series=[[s, 1.0e-3 * s * (1 + 2.0 * en_band(s)), 0.01 * s]
                             for s in range(1, 25)])
    drift = dict(near, field_norms={"h": [100.0, 2.0], "u": [3.0 * (1 + 2.0 * fb), 0.1]})
    leak = dict(near, budget_worst={"Salt": 1e-10})
    short = dict(near, series=near["series"][:12])
    nonorm = {k: v for k, v in near.items() if k != "field_norms"}
    got = [cross_backend_check(x, ref_m)[0] for x in (near, far, drift, leak, short, nonorm)]
    out.append((got == [True, False, False, False, False, False],
                "cross-backend: in band / En out / field norm out / leak / short run / no "
                "norms {}".format(got)))
    # --- the per-PR slice ---
    vp = {("eddy", "gm"): ["src/x/gm.F90"], ("geometry", "obc"): ["src/b/obc*.F90"]}
    t = [touched_values(["src/x/gm.F90", "docs/a.md"], vp, ["tests/compat_*.py"]),
         touched_values(["src/b/obc_baroclinic.F90"], vp, []),
         touched_values(["src/core/dyn.F90"], vp, []),
         touched_values(["docs/a.md"], vp, []),
         touched_values(["tests/compat_x.py"], vp, ["tests/compat_*.py"])]
    want = [("slice", {("eddy", "gm")}), ("slice", {("geometry", "obc")})]
    out.append((t[:2] == want and t[2][0] == "full" and t[3][0] == "none" and t[4][0] == "full",
                "slice: claimed path -> its value; shared source -> full; docs -> none"))
    shutil.rmtree(scratch, ignore_errors=True)
    return out
