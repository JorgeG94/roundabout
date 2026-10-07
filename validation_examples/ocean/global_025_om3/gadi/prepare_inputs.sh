#!/bin/bash
# Build the roundabout INPUT/ directory for the global 0.25 degree OM3
# 2-node / 8-GPU scaling test, from the real ACCESS-OM3 25 km JRA55-do IAF
# config inputs on /g/data (source: release-MC_25km_jra_iaf,
# github.com/ACCESS-NRI/access-om3-configs -- see
# tmp_local_artifacts/global025/OM3_25KM_INPUTS.md for every path + its
# citation and every UNKNOWN). Run on a Gadi login or copyq node (needs
# /g/data read access + enough local disk for the derived files -- NOT a
# compute-node job).
#
# Usage:
#   RDB_DATA_DIR=/scratch/<project>/<user>/rdb_data ./prepare_inputs.sh
#
# Writes to $RDB_DATA_DIR/OM3_025/ (created if absent). Every /g/data
# source path is checked BEFORE anything is copied; any missing path fails
# loud, listing every missing path at once (not just the first), so a
# stale module/path assumption surfaces immediately instead of partway
# through a multi-hour prep.
set -euo pipefail

if [[ -z "${RDB_DATA_DIR:-}" ]]; then
    echo "ERROR: RDB_DATA_DIR is not set. Example:" >&2
    echo "  RDB_DATA_DIR=/scratch/<project>/<user>/rdb_data $0" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="$RDB_DATA_DIR/OM3_025"
mkdir -p "$OUT"

# ---- STEP 1: every /g/data source path (OM3_25KM_INPUTS.md citations) ----
HGRID_SRC=/g/data/vk83/configurations/inputs/access-om3/share/grids/global.25km/2026.06.11/ocean_hgrid.nc
TOPOG_SRC=/g/data/vk83/configurations/inputs/access-om3/share/grids/global.25km/2026.06.11/topog.nc
VGRID_SRC=/g/data/vk83/configurations/inputs/access-om3/mom/grids/vertical/global.25km/2026.03.16/ocean_vgrid.nc
IC_SRC=/g/data/vk83/configurations/inputs/access-om3/mom/initial_conditions/global.25km/2026.06.11/woa23_ts_01_mom.nc
# check: this is a DIRECTORY of per-variable/per-year input4MIPs files, not
# one file -- confirm the uas/vas filenames inside it against what
# tools/om1deg_wind_regrid.f90's glob expects (README below, Sec. "Wind").
JRA_ATMOS_DIR=/g/data/qv56/replicas/input4MIPs/CMIP6Plus/OMIP/MRI/MRI-JRA55-do-1-6-0/atmos

missing=0
for p in "$HGRID_SRC" "$TOPOG_SRC" "$VGRID_SRC" "$IC_SRC" "$JRA_ATMOS_DIR"; do
    if [[ ! -e "$p" ]]; then
        echo "MISSING: $p" >&2
        missing=1
    fi
done
if [[ "$missing" -ne 0 ]]; then
    echo "ERROR: one or more OM3 25km input paths above do not exist on this" >&2
    echo "system. These come from release-MC_25km_jra_iaf's config.yaml /" >&2
    echo "MOM_input as fetched on 2026-10-07 -- the branch may have moved to a" >&2
    echo "newer input vintage (directory date stamps change) since then. Check" >&2
    echo "https://github.com/ACCESS-NRI/access-om3-configs/tree/release-MC_25km_jra_iaf" >&2
    echo "for the current paths and edit this script's SRC variables." >&2
    exit 1
fi

echo "All /g/data source paths present. Copying into $OUT ..."
cp -n "$HGRID_SRC" "$OUT/ocean_hgrid.nc"
cp -n "$TOPOG_SRC" "$OUT/topog.nc"
cp -n "$VGRID_SRC" "$OUT/ocean_vgrid.nc"     # kept for reference -- see nml header, not read directly
cp -n "$IC_SRC" "$OUT/woa23_ts_01_mom.nc"

# ---- STEP 2: bathymetry -- MOM6 limits (MINIMUM_DEPTH/MASKING_DEPTH are
# UNKNOWN for this OM3 config, see OM3_25KM_INPUTS.md -- reusing OM4_025's
# values per the task brief's "closest known analogue"; override via env
# if the maintainer finds the real values) ----
MIN_DEPTH="${OM3_MIN_DEPTH:-9.5}"
MAX_DEPTH="${OM3_MAX_DEPTH:-6000.0}"       # release-MC_25km_jra_iaf MOM_input, confirmed
MASKING_DEPTH="${OM3_MASKING_DEPTH:-0.0}"

module use /g/data/xp65/public/modules 2>/dev/null || true
module load conda/analysis3 2>/dev/null || true
python3 "$REPO_ROOT/tools/om_prepare_bathy.py" \
    --topog "$OUT/topog.nc" --hgrid "$OUT/ocean_hgrid.nc" \
    --out "$OUT/topog_limited.nc" \
    --min-depth "$MIN_DEPTH" --max-depth "$MAX_DEPTH" --masking-depth "$MASKING_DEPTH"

# ---- STEP 3: T/S IC -- vertical fill-gap repair + layout conversion.
# woa23_ts_01_mom.nc is ALREADY on the OM3 model grid (MOM6
# TEMP_SALT_Z_INIT_FILE convention), so no horizontal interpolation is
# needed -- same situation as southern_ocean_025's WOA05 file. Variable
# names inside it are UNKNOWN from here; override via env after
# `ncdump -h woa23_ts_01_mom.nc` on Gadi if these defaults are wrong. ----
TEMP_VAR="${OM3_TEMP_VAR:-ptemp}"
SALT_VAR="${OM3_SALT_VAR:-salt}"
LEVEL_VAR="${OM3_LEVEL_VAR:-level}"

# check: confirm these module names/versions with `module avail` on Gadi.
module load nvhpc 2>/dev/null || true
module load misc/nvhpc-build/25.5/netcdf-c misc/nvhpc-build/25.5/netcdf-fortran 2>/dev/null || true
module load hdf5 2>/dev/null || true

nvfortran -O2 $(nf-config --fflags) -o "$OUT/om_zclim_to_zinit" \
    "$REPO_ROOT/tools/om_zclim_to_zinit.f90" $(nf-config --flibs)
"$OUT/om_zclim_to_zinit" "$OUT/woa23_ts_01_mom.nc" "$OUT/ic_om3_zinit.nc" \
    --temp-var "$TEMP_VAR" --salt-var "$SALT_VAR" --level-var "$LEVEL_VAR"

# ---- STEP 4: wind stress -- JRA55-do 1958 daily-mean, regridded onto OM3's
# own supergrid. check: the qv56 input4MIPs atmos directory's uas/vas
# filenames may not match the *.padded.nc glob
# tools/om1deg_wind_regrid.f90 expects (that padding was done for the
# OM_1deg/OM4_025 recipes' own JRA55do copy, see
# ../../southern_ocean_025/README.md Sec. 2); if this step fails on a glob
# miss, either re-run tools/fetch_om1deg.py --jra-wind's padding step
# against the qv56 files first, or point --jra at a directory that already
# has them. ----
nvfortran -O2 $(nf-config --fflags) -o "$OUT/om1deg_subset" \
    "$REPO_ROOT/tools/om1deg_subset.f90" $(nf-config --flibs) 2>/dev/null || true
python3 "$REPO_ROOT/tools/om1deg_prepare_wind.py" \
    --hgrid "$OUT/ocean_hgrid.nc" --jra "$JRA_ATMOS_DIR" \
    --out "$OUT/wind_jra55do_1958_24h.nc" \
    --year 1958 --cd ly04 --avg-hours 24 --fc nvfortran \
    --build-dir "$OUT"

echo "Done. INPUT files are in $OUT"
echo "  ln -s $OUT INPUT   # from the run directory"
