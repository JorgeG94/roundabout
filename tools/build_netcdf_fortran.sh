#!/usr/bin/env bash
# -----------------------------------------------------------------------
# build_netcdf_fortran.sh — build netcdf-fortran for ONE Fortran compiler
#
# Roundabout's only compiler-coupled dependency is netcdf-fortran. Fortran
# `.mod` files and module symbol manglings are compiler-specific, so a
# netcdf-fortran built by gfortran cannot be consumed by nvfortran, ifx or
# flang -- you need one per compiler.
#
# Everything BELOW it -- netcdf-c, HDF5, zlib -- is C, and a C library is
# reusable across all of these compilers on Linux x86-64. So you do NOT
# need a per-compiler netcdf-c: take a prebuilt one from your distro,
# module tree, or conda, and build only this Fortran wrapper on top. That
# is minutes of work instead of an HDF5-from-source stack.
#
# Usage:
#   tools/build_netcdf_fortran.sh --fc nvfortran
#   tools/build_netcdf_fortran.sh --fc ifx --prefix ~/opt/netcdff-ifx
#   tools/build_netcdf_fortran.sh --fc gfortran --netcdf-c /usr
#
# netcdf-c is located in this order:
#   1. --netcdf-c <prefix>
#   2. $NETCDFC_DIR
#   3. `nc-config --prefix`, if nc-config is on PATH
#
# Note that Roundabout needs NEITHER parallel netcdf NOR HDF5's Fortran
# bindings: its I/O is per-rank serial `nf90_create` with an offline merge
# (tools/merge_output.py), and it imports `netcdf` and nothing else. A
# serial netcdf-c against a serial `hdf5 ~fortran` is enough.
# -----------------------------------------------------------------------
set -euo pipefail

NCF_VERSION="4.6.2"
NCF_SHA256="df26b99d9003c93a8bc287b58172bf1c279676f8c10d6dd0daf8bc7204877096"
NCF_URL_BASE="https://downloads.unidata.ucar.edu/netcdf-fortran"

FC="${FC:-}"
PREFIX=""
NETCDF_C=""
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
KEEP_BUILD=0

die() { printf '\nerror: %s\n' "$*" >&2; exit 1; }
note() { printf '%s\n' "$*"; }

usage() {
   sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'
   exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
   case "$1" in
      --fc)        FC="${2:-}"; shift 2 ;;
      --prefix)    PREFIX="${2:-}"; shift 2 ;;
      --netcdf-c)  NETCDF_C="${2:-}"; shift 2 ;;
      --version)   NCF_VERSION="${2:-}"; NCF_SHA256=""; shift 2 ;;
      --jobs|-j)   JOBS="${2:-}"; shift 2 ;;
      --keep-build) KEEP_BUILD=1; shift ;;
      -h|--help)   usage 0 ;;
      *)           printf 'unknown argument: %s\n' "$1" >&2; usage 1 ;;
   esac
done

# --- the Fortran compiler ----------------------------------------------
[[ -n "$FC" ]] || die "no Fortran compiler given. Pass --fc <compiler> (or set FC).
Roundabout supports gfortran (>= 15), nvfortran, ifx, and flang-new."
command -v "$FC" >/dev/null 2>&1 || die "Fortran compiler '$FC' is not on PATH.
Load your toolchain first (module load / source setvars.sh / conda activate)."
FC_PATH="$(command -v "$FC")"

# A stable tag for the default prefix, so two compilers never collide.
FC_TAG="$(basename "$FC")"

# --- netcdf-c ----------------------------------------------------------
if [[ -z "$NETCDF_C" && -n "${NETCDFC_DIR:-}" ]]; then
   NETCDF_C="$NETCDFC_DIR"
   NETCDF_C_SRC="\$NETCDFC_DIR"
elif [[ -z "$NETCDF_C" ]] && command -v nc-config >/dev/null 2>&1; then
   NETCDF_C="$(nc-config --prefix)"
   NETCDF_C_SRC="nc-config on PATH"
else
   NETCDF_C_SRC="--netcdf-c"
fi

[[ -n "$NETCDF_C" ]] || die "could not find netcdf-c.

Give it explicitly with --netcdf-c <prefix>, or install one -- ANY C build
will do, it does not have to match your Fortran compiler:

  Debian/Ubuntu   sudo apt install libnetcdf-dev
  Fedora/RHEL     sudo dnf install netcdf-devel
  conda           conda install -c conda-forge libnetcdf
  HPC site        module load netcdf-c   (name varies)"

[[ -d "$NETCDF_C" ]] || die "netcdf-c prefix '$NETCDF_C' does not exist (from $NETCDF_C_SRC)."
[[ -f "$NETCDF_C/include/netcdf.h" ]] || die \
   "'$NETCDF_C' has no include/netcdf.h -- that is not a netcdf-c prefix (from $NETCDF_C_SRC)."

NC_CONFIG="$NETCDF_C/bin/nc-config"
if [[ -x "$NC_CONFIG" ]]; then
   NC_VERSION="$("$NC_CONFIG" --version 2>/dev/null || echo unknown)"
else
   NC_VERSION="unknown (no nc-config in that prefix)"
fi

# --- prefix ------------------------------------------------------------
[[ -n "$PREFIX" ]] || PREFIX="$HOME/opt/netcdf-fortran-${NCF_VERSION}-${FC_TAG}"

note "-------------------------------------------------------------------"
note "  netcdf-fortran   ${NCF_VERSION}"
note "  Fortran compiler ${FC_PATH}"
note "  netcdf-c         ${NETCDF_C}"
note "                   ${NC_VERSION}  (found via ${NETCDF_C_SRC})"
note "  install prefix   ${PREFIX}"
note "  parallel jobs    ${JOBS}"
note "-------------------------------------------------------------------"

for tool in cmake curl tar; do
   command -v "$tool" >/dev/null 2>&1 || die "'$tool' is required but not on PATH."
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/netcdff-build-XXXXXX")"
cleanup() {
   if [[ $KEEP_BUILD -eq 1 ]]; then
      note "build tree kept at ${WORK}"
   else
      rm -rf "$WORK"
   fi
}
trap cleanup EXIT

TARBALL="$WORK/netcdf-fortran-${NCF_VERSION}.tar.gz"
URL="${NCF_URL_BASE}/${NCF_VERSION}/netcdf-fortran-${NCF_VERSION}.tar.gz"

note ""
note "==> downloading ${URL}"
curl -fsSL -o "$TARBALL" "$URL" || die "download failed: $URL"

if [[ -n "$NCF_SHA256" ]]; then
   note "==> verifying sha256"
   GOT="$(sha256sum "$TARBALL" | awk '{print $1}')"
   [[ "$GOT" == "$NCF_SHA256" ]] || die "checksum mismatch for $URL
  expected  $NCF_SHA256
  got       $GOT"
else
   note "==> skipping checksum (pin overridden by --version)"
fi

note "==> unpacking"
tar xzf "$TARBALL" -C "$WORK"
SRC="$WORK/netcdf-fortran-${NCF_VERSION}"
[[ -d "$SRC" ]] || die "unpacked tree not found at $SRC"

note "==> configuring"
# netcdf-fortran finds netcdf-c through find_package(netCDF) first and falls
# back to FIND_LIBRARY; CMAKE_PREFIX_PATH satisfies both. The generated
# netCDF-FortranConfig.cmake then records THIS netcdf-c, which is what
# Roundabout's find_package(netCDF-Fortran) will resolve against.
cmake -S "$SRC" -B "$WORK/build" \
   -DCMAKE_Fortran_COMPILER="$FC_PATH" \
   -DCMAKE_INSTALL_PREFIX="$PREFIX" \
   -DCMAKE_PREFIX_PATH="$NETCDF_C" \
   -DCMAKE_BUILD_TYPE=Release \
   -DBUILD_SHARED_LIBS=ON \
   -DENABLE_TESTS=OFF \
   -DBUILD_EXAMPLES=OFF \
   -DCMAKE_INSTALL_RPATH="$PREFIX/lib;$NETCDF_C/lib" \
   -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
   > "$WORK/configure.log" 2>&1 \
   || { tail -40 "$WORK/configure.log"; die "configure failed; full log: $WORK/configure.log
(re-run with --keep-build to keep it)"; }

note "==> building (-j${JOBS})"
cmake --build "$WORK/build" -j "$JOBS" > "$WORK/build.log" 2>&1 \
   || { tail -40 "$WORK/build.log"; die "build failed; full log: $WORK/build.log
(re-run with --keep-build to keep it)"; }

note "==> installing"
cmake --install "$WORK/build" > "$WORK/install.log" 2>&1 \
   || { tail -40 "$WORK/install.log"; die "install failed; full log: $WORK/install.log"; }

[[ -f "$PREFIX/include/netcdf.mod" ]] || die "install finished but $PREFIX/include/netcdf.mod is missing."

# --- prove it actually works before claiming success --------------------
note "==> smoke test (compile, link, write a file, read it back)"
cat > "$WORK/smoke.f90" <<'EOF'
program smoke
   use netcdf
   implicit none
   integer :: ncid, dimid, varid, ierr
   real :: wrote(4) = [1.0, 2.0, 3.0, 4.0], readback(4)
   ierr = nf90_create("smoke.nc", NF90_CLOBBER, ncid); if (ierr /= nf90_noerr) stop 1
   ierr = nf90_def_dim(ncid, "x", 4, dimid);            if (ierr /= nf90_noerr) stop 2
   ierr = nf90_def_var(ncid, "v", NF90_FLOAT, [dimid], varid)
   ierr = nf90_enddef(ncid)
   ierr = nf90_put_var(ncid, varid, wrote)
   ierr = nf90_close(ncid);                             if (ierr /= nf90_noerr) stop 3
   ierr = nf90_open("smoke.nc", NF90_NOWRITE, ncid);    if (ierr /= nf90_noerr) stop 4
   ierr = nf90_inq_varid(ncid, "v", varid)
   ierr = nf90_get_var(ncid, varid, readback)
   ierr = nf90_close(ncid)
   if (any(abs(readback - wrote) > 0.0)) stop 5
   print *, "netcdf-fortran smoke test OK, library ", trim(nf90_inq_libvers())
end program
EOF
(
   cd "$WORK"
   "$FC_PATH" -o smoke smoke.f90 \
      -I"$PREFIX/include" -L"$PREFIX/lib" -lnetcdff \
      -L"$NETCDF_C/lib" -lnetcdf > smoke_compile.log 2>&1 \
      || { tail -20 smoke_compile.log; die "smoke test failed to compile/link."; }
   LD_LIBRARY_PATH="$PREFIX/lib:$NETCDF_C/lib:${LD_LIBRARY_PATH:-}" ./smoke \
      || die "smoke test compiled but failed at runtime."
)

note ""
note "-------------------------------------------------------------------"
note "netcdf-fortran ${NCF_VERSION} installed for ${FC_TAG}"
note ""
note "Configure Roundabout against it with:"
note ""
note "  cmake -B build_${FC_TAG} -S . \\"
note "        -DCMAKE_Fortran_COMPILER=${FC_PATH} \\"
note "        -DCMAKE_PREFIX_PATH='${PREFIX};${NETCDF_C}'"
note ""
note "At runtime the libraries must be findable:"
note ""
note "  export LD_LIBRARY_PATH=${PREFIX}/lib:${NETCDF_C}/lib:\$LD_LIBRARY_PATH"
note ""
note "(The install carries an RPATH for both, so this is usually already"
note "handled -- set it if you relocate either prefix.)"
note "-------------------------------------------------------------------"
