#!/usr/bin/env python3
"""Fetch the MOM6 OM_1deg grid/topography + WOA13 January T/S (standard library only).

Reproduces the input set of `validation_examples/ocean/global_1deg/`:

* ``OM_1deg.tgz`` (~9 MB) from GFDL's public MOM6-testing ftp tree —
  ``ocean_hgrid.nc`` (the 720x640 supergrid), ``topog.nc``, ``ocean_mask.nc``,
  ``vgrid_75_2m.nc``, ... — extracted in full.
* ``obs.woa13.tgz`` (~1.06 GB) from the same tree — STREAMED; only the two
  WOA13 decav monthly full-depth files are written
  (``woa13_decav_ptemp_monthly_fulldepth_01.nc``,
  ``woa13_decav_s_monthly_fulldepth_01.nc``); the tarball itself is never
  stored.

Nothing is committed to the repository: data land in ``--dest DIR``, or in
``$RDB_DATA_DIR/OM_1deg`` when ``--dest`` is omitted (one of the two is
required). Moving the data later is just ``mv`` + a new ``RDB_DATA_DIR``.

Idempotent: a file already present with the expected SHA-256 is skipped, and an
archive is only downloaded when one of its wanted members is missing or wrong.
The transferred byte count of each archive is checked against its published
size, and every extracted file against its recorded SHA-256.

    python3 tools/fetch_om1deg.py [--dest DIR] [--skip-woa]

``--jra-wind [--year 1958]`` instead fetches one year of JRA55-do v1.4.0
near-surface wind (``uas``, ``vas``; 3-hourly point values on the TL319
640x320 grid) for ``global_1deg_wind.nml``: ``reanalysis.tgz`` (~12.4 GB, the
OM_1deg forcing set) is STREAMED from the same GFDL tree and only the two wind
files are written (~3.4 GB), to ``--jra-dest DIR`` or
``$RDB_DATA_DIR/JRA55do``. GFDL's tarball holds the year 1958 only (the
``padded`` files: the 2920 records of 1958 plus 1959-01-01 00:00); any other
year has to come from input4MIPs/ESGF (see the global_1deg README).

    python3 tools/fetch_om1deg.py --jra-wind [--year 1958] [--jra-dest DIR]
"""

import argparse
import hashlib
import os
import sys
import tarfile
import urllib.request

BASE = "ftp://ftp.gfdl.noaa.gov/perm/Alistair.Adcroft/MOM6-testing/"

# archive -> (published byte size, {member basename: sha256})
ARCHIVES = {
    "OM_1deg.tgz": (9337533, {
        "ocean_hgrid.nc": "247c01a410e88760ca724edba4447aec2adb17c06fc769beaa0742236c342666",
        "topog.nc": "14336d911a69be668c18b8fb0a2d587853c7e58982f4f9854c33f8d0cc7cbd00",
        "ocean_mask.nc": "0bf9e5a456207e97712ee88339471348e13adcf6b16f6835f7ea66256a8710da",
        "vgrid_75_2m.nc": "b4e027c57e0f18177a2a928112e89741235cb239b707b521fa909b6c436754af",
    }),
    "obs.woa13.tgz": (1058331974, {
        "woa13_decav_ptemp_monthly_fulldepth_01.nc":
            "bcc2472907c9501c8c3101bfc8840f7540f27dc7b160dec7f3ff14d11d52614a",
        "woa13_decav_s_monthly_fulldepth_01.nc":
            "d32c66aedf44ea0220b2491cb5ed77801ea5fa31a05851eb4115221b98a60fdb",
    }),
}
# OM_1deg.tgz: extract every regular .nc member (grid, masks, vgrids, ...).
EXTRACT_ALL = {"OM_1deg.tgz"}

# JRA55-do v1.4.0 near-surface wind, per year: archive -> (size, members).
# reanalysis.tgz carries all nine OM_1deg state/flux variables for 1958
# (``reanalysis/JRA55-do/v1.4.0/padded/*.padded.nc``); only uas/vas are kept.
JRA_TAG = "input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-4-0_gr"
JRA_ARCHIVES = {
    1958: {"reanalysis.tgz": (12382889092, {
        f"uas_{JRA_TAG}_195801010000-195812312100.padded.nc":
            "e8e52ae9457bcb57e56609a2e5b27248975711db85e71f3e27545a2111f85aeb",
        f"vas_{JRA_TAG}_195801010000-195812312100.padded.nc":
            "ecbb30e1ec3baaa83d82286efe3350a09293083f9e43599671939c349707e892",
    })},
}


def sha256(path, bufsize=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def member_ok(dest, name, want):
    path = os.path.join(dest, name)
    if not os.path.isfile(path):
        return False
    return sha256(path) == want


class CountingReader:
    """File-like wrapper that counts bytes and prints coarse progress."""

    def __init__(self, raw, total, label):
        self.raw, self.total, self.label = raw, total, label
        self.n, self._next = 0, 0

    def read(self, size=-1):
        b = self.raw.read(size)
        self.n += len(b)
        if self.n >= self._next:
            pct = 100.0 * self.n / self.total if self.total else 0.0
            sys.stderr.write(f"\r  {self.label}: {self.n / 1e6:9.1f} MB ({pct:5.1f} %)")
            sys.stderr.flush()
            self._next = self.n + 20_000_000
        return b


def fetch_archive(dest, archive, size, members):
    url = BASE + archive
    print(f"fetching {url} (streamed, {size / 1e6:.1f} MB)")
    extract_all = archive in EXTRACT_ALL
    got = []
    with urllib.request.urlopen(url, timeout=120) as resp:
        reader = CountingReader(resp, size, archive)
        with tarfile.open(fileobj=reader, mode="r|gz") as tar:
            for ti in tar:
                base = os.path.basename(ti.name)
                if not ti.isfile():
                    continue
                if base not in members and not (extract_all and base.endswith(".nc")):
                    print(f"  (skipping member {ti.name}, {ti.size} bytes)", file=sys.stderr)
                    continue
                src = tar.extractfile(ti)
                tmp = os.path.join(dest, base + ".part")
                with open(tmp, "wb") as out:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        out.write(b)
                os.replace(tmp, os.path.join(dest, base))
                got.append(base)
            # Drain the rest of the stream so the byte count is complete.
            while reader.read(1 << 20):
                pass
    sys.stderr.write("\n")
    if reader.n != size:
        raise SystemExit(f"{archive}: transferred {reader.n} bytes, expected {size}")
    missing = [m for m in members if m not in got]
    if missing:
        raise SystemExit(f"{archive}: members not found in archive: {missing}")


def main():
    root = os.environ.get("RDB_DATA_DIR")
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--dest", default=None,
                    help="output directory (default: $RDB_DATA_DIR/OM_1deg)")
    ap.add_argument("--skip-woa", action="store_true", help="skip the 1 GB WOA13 archive")
    ap.add_argument("--jra-wind", action="store_true",
                    help="fetch one year of JRA55-do uas/vas instead of the OM_1deg set")
    ap.add_argument("--year", type=int, default=1958, help="JRA55-do year (default 1958)")
    ap.add_argument("--jra-dest", default=None,
                    help="JRA output directory (default: $RDB_DATA_DIR/JRA55do)")
    args = ap.parse_args()
    if args.jra_wind:
        fetch_jra_wind(args, root, ap)
        return
    if args.dest is None and not root:
        ap.error("no destination: pass --dest DIR, or set RDB_DATA_DIR (data then go to "
                 "$RDB_DATA_DIR/OM_1deg). About 4 GB of free space is needed; the data "
                 "are re-downloadable, so a scratch/RAID partition is fine.")
    dest = os.path.abspath(args.dest or os.path.join(root, "OM_1deg"))
    os.makedirs(dest, exist_ok=True)

    fetch_all(dest, ARCHIVES, skip={"obs.woa13.tgz"} if args.skip_woa else set())


def fetch_all(dest, archives, skip=()):
    """Fetch every archive whose wanted members are not all present + verified."""
    for archive, (size, members) in archives.items():
        if archive in skip:
            continue
        if all(member_ok(dest, m, h) for m, h in members.items()):
            print(f"{archive}: all members present and verified — skipped")
            continue
        fetch_archive(dest, archive, size, members)
        for m, h in members.items():
            path = os.path.join(dest, m)
            digest = sha256(path)
            if digest != h:
                raise SystemExit(f"{m}: sha256 {digest} != expected {h}")
            print(f"  {m}: {os.path.getsize(path)} bytes sha256 {digest}")

    print(f"\ndata directory: {dest}")
    for name in sorted(os.listdir(dest)):
        p = os.path.join(dest, name)
        if os.path.isfile(p):
            print(f"  {os.path.getsize(p):>12d}  {name}")


def fetch_jra_wind(args, root, ap):
    """One year of JRA55-do v1.4.0 uas/vas (see the module docstring)."""
    if args.year not in JRA_ARCHIVES:
        raise SystemExit(
            f"--year {args.year}: GFDL's reanalysis.tgz holds only "
            f"{sorted(JRA_ARCHIVES)}. For another year download "
            f"uas/vas_{JRA_TAG}_<year>01010000-<year>12312100.nc from input4MIPs "
            "(ESGF, project input4MIPs, source MRI-JRA55-do-1-4-0) into the JRA "
            "directory by hand; tools/om1deg_prepare_wind.py takes --year.")
    if args.jra_dest is None and not root:
        ap.error("no destination: pass --jra-dest DIR, or set RDB_DATA_DIR (data then "
                 "go to $RDB_DATA_DIR/JRA55do). About 3.5 GB of free space is needed.")
    dest = os.path.abspath(args.jra_dest or os.path.join(root, "JRA55do"))
    os.makedirs(dest, exist_ok=True)
    fetch_all(dest, JRA_ARCHIVES[args.year])


if __name__ == "__main__":
    main()
