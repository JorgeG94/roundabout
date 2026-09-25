#!/usr/bin/env bash
# Reproduce the global 1-degree unforced run in one command, and check it
# against the committed reference (reference_daily.csv).
#
#   validation_examples/ocean/global_1deg/reproduce.sh --build-dir BUILD [options]
#
# Steps: fetch the OM_1deg + WOA13 inputs (idempotent) -> prepare the model-grid
# bathymetry and initial condition (skipped when up to date) -> set up the run
# directory (INPUT symlink + namelist copy with t_end = --days) -> run ->
# movie -> check_against_reference.py.  Standard-library Python only; nothing
# is installed.  Run with --help for the options.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NML=global_1deg_unforced.nml
QUICK_DAYS=10

usage() {
    cat <<EOF
Usage: $(basename "$0") --build-dir DIR [options]

  --data-dir DIR   data root holding OM_1deg/ (default: \$RDB_DATA_DIR; ~4 GB,
                   re-downloadable, so any scratch disk will do)
  --build-dir DIR  roundabout build tree: the rdb executable and, for
                   --python, librdb_core.so (configure -DRDB_BUILD_SHARED=ON)
                   (default: \$RDB_BUILD_DIR)
  --run-dir DIR    where the run happens (default: ./global_1deg_run)
  --days N         simulated days (default: 365)
  --quick          $QUICK_DAYS days, checked against the reference's first
                   $QUICK_DAYS days (about two minutes on one V100)
  --python         drive the run through the Python interface
                   (run_global_1deg.py) instead of the rdb executable; the
                   movie is rendered in-process as the run goes
  --movie          require the movie (fail early if a tool is missing)
  --no-movie       skip the movie
  --gpu N          export CUDA_VISIBLE_DEVICES=N (one device per run)
  --skip-fetch     do not run tools/fetch_om1deg.py (data already in place)
  --no-check       do not compare against the reference
  -h, --help       this text

By default the rdb executable runs the namelist and the movie is rendered
afterwards from its diagnostic file (global_movie.py; needs nccopy from
NetCDF-C and ffmpeg — without them the movie is skipped with a note).
EOF
}

die() { echo "reproduce.sh: $*" >&2; exit 2; }

data_root="${RDB_DATA_DIR:-}"
build_dir="${RDB_BUILD_DIR:-}"
run_dir="./global_1deg_run"
days=365
quick=0
use_python=0
movie=auto
skip_fetch=0
do_check=1
while [ $# -gt 0 ]; do
    case "$1" in
        --data-dir) data_root="${2:?--data-dir needs a value}"; shift 2 ;;
        --build-dir) build_dir="${2:?--build-dir needs a value}"; shift 2 ;;
        --run-dir) run_dir="${2:?--run-dir needs a value}"; shift 2 ;;
        --days) days="${2:?--days needs a value}"; shift 2 ;;
        --quick) quick=1; shift ;;
        --python) use_python=1; shift ;;
        --movie) movie=yes; shift ;;
        --no-movie) movie=no; shift ;;
        --gpu) export CUDA_VISIBLE_DEVICES="${2:?--gpu needs a device number}"; shift 2 ;;
        --skip-fetch) skip_fetch=1; shift ;;
        --no-check) do_check=0; shift ;;
        -h | --help) usage; exit 0 ;;
        *) usage >&2; die "unknown option: $1" ;;
    esac
done
[ "$quick" -eq 1 ] && days=$QUICK_DAYS
case "$days" in '' | *[!0-9]* | 0) die "--days must be a positive integer (got '$days')" ;; esac

# --- prerequisites ----------------------------------------------------------
[ -n "$data_root" ] || die "no data directory: pass --data-dir DIR or set RDB_DATA_DIR (the inputs, ~4 GB, are fetched into DIR/OM_1deg)"
[ -n "$build_dir" ] || die "no build directory: pass --build-dir DIR (a roundabout build tree with the rdb executable), or set RDB_BUILD_DIR"
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null || die "$PY not found (standard-library Python 3 is all that is needed)"
mkdir -p "$data_root"
data_root="$(cd "$data_root" && pwd)"
build_dir="$(cd "$build_dir" 2>/dev/null && pwd)" || die "build directory not found: $build_dir"
data="$data_root/OM_1deg"
if [ "$use_python" -eq 1 ]; then
    lib="$build_dir/librdb_core.so"
    [ -f "$lib" ] || die "--python needs $lib (configure the build with -DRDB_BUILD_SHARED=ON)"
else
    exe="$build_dir/rdb"
    [ -x "$exe" ] || die "no rdb executable in $build_dir"
fi
if [ "$movie" != no ]; then
    need="ffmpeg"
    [ "$use_python" -eq 0 ] && need="nccopy ffmpeg"
    for tool in $need; do
        if ! command -v "$tool" >/dev/null; then
            [ "$movie" = yes ] && die "--movie needs $tool on PATH"
            echo "note: $tool not found — the movie is skipped (the run and the check are not affected)"
            movie=no
        fi
    done
fi

# --- 1. inputs --------------------------------------------------------------
echo "== 1/5 inputs: $data"
if [ "$skip_fetch" -eq 0 ]; then
    "$PY" "$REPO/tools/fetch_om1deg.py" --dest "$data"
else
    echo "fetch skipped (--skip-fetch)"
fi
raw="topog.nc ocean_hgrid.nc woa13_decav_ptemp_monthly_fulldepth_01.nc woa13_decav_s_monthly_fulldepth_01.nc"
for f in $raw; do
    [ -f "$data/$f" ] || die "$data/$f is missing (run without --skip-fetch)"
done

echo "== 2/5 prepare model-grid inputs"
stale=0
for out in bathy_om1deg.nc ic_woa13_jan.nc; do
    [ -f "$data/$out" ] || { stale=1; continue; }
    for f in $raw; do
        [ "$data/$out" -nt "$data/$f" ] || stale=1
    done
done
if [ "$stale" -eq 1 ]; then
    "$PY" "$REPO/tools/om1deg_prepare_inputs.py" --data "$data"
else
    echo "bathy_om1deg.nc and ic_woa13_jan.nc are newer than their inputs — skipped"
fi

# --- 3. run directory -------------------------------------------------------
echo "== 3/5 run directory"
mkdir -p "$run_dir"
run_dir="$(cd "$run_dir" && pwd)"
if [ -L "$run_dir/INPUT" ]; then
    rm "$run_dir/INPUT"
elif [ -e "$run_dir/INPUT" ]; then
    die "$run_dir/INPUT exists and is not a symlink; refusing to replace it"
fi
ln -s "$data" "$run_dir/INPUT"
sed "s/^\([[:space:]]*t_end[[:space:]]*=[[:space:]]*\)[0-9.]*/\1${days}.0/" \
    "$HERE/$NML" >"$run_dir/$NML"
grep -q "t_end *= *${days}\.0" "$run_dir/$NML" || die "could not set t_end in $run_dir/$NML"
rm -f "$run_dir"/output/global_1deg_rank_*.nc
log="$run_dir/run.log"
echo "run directory $run_dir: INPUT -> $data, $NML (t_end = $days days)"

# --- 4. run -----------------------------------------------------------------
echo "== 4/5 run ($days days, $([ "$use_python" -eq 1 ] && echo "Python interface" || echo "rdb executable"), CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset})"
t0=$(date +%s)
if [ "$use_python" -eq 1 ]; then
    pyargs=(--days "$days" --data "$data" --out "$run_dir/output")
    [ "$movie" = no ] && pyargs+=(--no-movie)
    (cd "$run_dir" && RDB_LIB="$lib" "$PY" "$HERE/run_global_1deg.py" "${pyargs[@]}") 2>&1 | tee "$log"
else
    (cd "$run_dir" && "$exe" "$NML") 2>&1 | tee "$log"
fi
echo "run: $(($(date +%s) - t0)) s wall"

# --- 5. movie + check -------------------------------------------------------
echo "== 5/5 movie and check"
diag="$run_dir/output/global_1deg_rank_000000.nc"
if [ "$movie" != no ] && [ "$use_python" -eq 0 ]; then
    "$PY" "$HERE/global_movie.py" "$diag" "$run_dir/movie" --data "$data"
fi
status=0
if [ "$do_check" -eq 1 ]; then
    if [ "$days" -le 365 ]; then
        "$PY" "$HERE/check_against_reference.py" "$log" --days "$days" | tee "$run_dir/check.txt" || status=1
    else
        echo "check skipped: the reference covers 365 days"
    fi
fi

echo
echo "data:        $data"
echo "run dir:     $run_dir"
echo "console log: $log"
[ -f "$diag" ] && echo "diag file:   $diag"
if [ "$movie" != no ]; then
    if [ "$use_python" -eq 1 ]; then
        echo "movie:       $run_dir/output/global_1deg.mp4 (+ .gif, frames/)"
    else
        echo "movie:       $run_dir/movie/global_1deg.mp4 (+ .gif, frames)"
    fi
fi
[ "$do_check" -eq 1 ] && echo "check:       $run_dir/check.txt"
exit "$status"
