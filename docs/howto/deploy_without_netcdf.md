# Deploy on a machine with HDF5 but no netcdf-c

Roundabout's I/O needs **netcdf-fortran built by your Fortran compiler**, on top
of a netcdf-c with NetCDF-4/HDF5 support. Most machines already have a netcdf-c
(distro, module tree, conda), and then
[`tools/build_netcdf_fortran.sh`](../../tools/build_netcdf_fortran.sh) builds
only the Fortran wrapper on top of it — see the README's *Getting
NetCDF-Fortran*.

This page is for the other case: a compiler, (maybe) MPI and an **HDF5**, but
**no netcdf-c you can use** — for example NERSC Perlmutter with `cray-hdf5` and
a `PrgEnv-*`. [`tools/build_netcdf.sh`](../../tools/build_netcdf.sh) builds
both netcdf-c and netcdf-fortran into one prefix, against the HDF5 you already
have. It never builds HDF5.

## What it does

1. Fetches pinned Unidata release tarballs and checks their SHA-256
   (`--download-only` / `--offline` split the fetch from the build).
2. Finds HDF5: `--hdf5 <prefix>`, else `$HDF5_DIR`, `$HDF5_ROOT` (set by
   `cray-hdf5`), then `h5cc -show` (`h5pcc` with `--mpi`), then
   `pkg-config hdf5`. If none of those finds it, the script stops and says so.
3. Builds **netcdf-c** with its autotools `configure`. Remote access (DAP,
   byte-range, S3), NCZarr, filter plugins, libxml2, examples and tests are
   off, so libcurl and libxml2 are not needed. NetCDF-4/HDF5 stays on, and
   szip is linked only if your HDF5 was built with it.
4. Builds **netcdf-fortran** with CMake, using your Fortran compiler.
5. Smoke-tests the pair: a Fortran program writes a deflated NetCDF-4 file
   and reads it back, finding the libraries through RPATH only.
6. Writes `<prefix>/roundabout-netcdf.env` and prints the `cmake` line for
   Roundabout.

| | Pinned version | SHA-256 (Unidata tarball) |
|---|---|---|
| netcdf-c | 4.9.3 | `a474149844e6144566673facf097fea253dc843c37bc0a7d3de047dc8adda5dd` |
| netcdf-fortran | 4.6.2 | `df26b99d9003c93a8bc287b58172bf1c279676f8c10d6dd0daf8bc7204877096` |

netcdf-c 4.10.1 and netcdf-fortran 4.6.4, the newest releases, are also in
the checksum table (`--netcdf-c-version 4.10.1`). Any other version needs `--netcdf-c-sha256` /
`--netcdf-fortran-sha256`. Without one, the script builds it anyway, warns,
and prints the hash it computed.

**Why autotools for netcdf-c?** `configure` accepts the Cray wrappers (`cc`)
and NVHPC's `nvc` as a plain `CC=`, and takes HDF5 from `CPPFLAGS`/`LDFLAGS`
whatever the install layout (Debian's `include/hdf5/serial`, Cray's
`$HDF5_DIR`). It does not rely on CMake's `FindHDF5` finding the right one.
netcdf-fortran stays on CMake, the same build CI uses for every compiler.

**Serial is enough.** Roundabout writes per-rank serial files with an offline
merge, so it needs no parallel NetCDF, even in an MPI build. `--mpi` exists
only for machines whose HDF5 is *parallel*: netcdf-c then has to be built with
the MPI compiler, and Roundabout has to be an MPI build
(`-DRDB_ENABLE_MPI=ON`). The script refuses a parallel HDF5 without `--mpi`,
and `--mpi` with a serial one.

## Generic Linux

```bash
module load gcc hdf5                     # or: apt install libhdf5-dev
tools/build_netcdf.sh --fc gfortran --prefix ~/opt/netcdf-gnu -j 8

source ~/opt/netcdf-gnu/roundabout-netcdf.env
cmake -B build_gnu -S . -DCMAKE_Fortran_COMPILER=gfortran \
      -DCMAKE_PREFIX_PATH=~/opt/netcdf-gnu
cmake --build build_gnu -j 8
```

The C compiler follows the Fortran one unless you pass `--cc` or set `$CC`
(`gfortran` → `gcc`, `nvfortran` → `nvc`, `ifx` → `icx`, `flang` → `clang`,
`ftn` → `cc`). Configure has to print

```
-- NetCDF-Fortran found via CMake config: <prefix>/lib/cmake/netCDF-Fortran
```

If `<dir>` is somewhere else, CMake found another netcdf-fortran first:
unload that module, or put the prefix first in `CMAKE_PREFIX_PATH`.

The installed libraries carry RPATHs to each other and to HDF5, so `rdb`
runs without `LD_LIBRARY_PATH`. The env file sets `LD_LIBRARY_PATH` anyway,
for a moved prefix. For NVHPC, the libraries target the build host's CPU
unless you pass `--cflags -tp=<cpu> --fflags -tp=<cpu>` (see
`RDB_EXTRA_FORTRAN_FLAGS` in the README). The script warns about this.

## NERSC Perlmutter

Perlmutter's login nodes and its CPU and GPU compute nodes are all AMD EPYC
Milan (zen3). A library built on a login node therefore runs on the compute
nodes. `-tp=zen3` makes that explicit for NVHPC.

Install into `/global/common/software/<project>`, which NERSC recommends for
software because it loads fastest at job start. `$SCRATCH` works too, but it
is purged.

### Fetch the sources (login node)

```bash
cd $HOME/roundabout
tools/build_netcdf.sh --download-only --download-dir $HOME/netcdf-src
```

Build on the login node, or on a compute node with
`--offline $HOME/netcdf-src` if that node has no outbound network.

### PrgEnv-gnu + cray-hdf5 (CPU build)

Roundabout needs gfortran 15 or newer. Load a `gcc-native` / `gcc` module
that provides it, or use PrgEnv-nvidia below.

```bash
module load PrgEnv-gnu cray-hdf5 cmake     # cray-hdf5 sets $HDF5_DIR
P=/global/common/software/<project>/netcdf-gnu
tools/build_netcdf.sh --offline $HOME/netcdf-src \
      --cc cc --fc ftn --prefix $P -j 16

source $P/roundabout-netcdf.env
cmake -B build_gnu -S . \
      -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_C_COMPILER=cc \
      -DCMAKE_PREFIX_PATH=$P
cmake --build build_gnu -j 16
```

### PrgEnv-nvidia + cray-hdf5 (A100 GPU build)

```bash
module load PrgEnv-nvidia cray-hdf5 cmake cudatoolkit craype-accel-nvidia80
P=/global/common/software/<project>/netcdf-nvidia
tools/build_netcdf.sh --offline $HOME/netcdf-src \
      --cc cc --fc ftn --prefix $P -j 16 \
      --cflags -tp=zen3 --fflags -tp=zen3

source $P/roundabout-netcdf.env
cmake -B build_cc80 -S . \
      -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_C_COMPILER=cc \
      -DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc80 \
      -DRDB_EXTRA_FORTRAN_FLAGS=-tp=zen3 \
      -DCMAKE_PREFIX_PATH=$P
cmake --build build_cc80 -j 16
```

For multi-GPU runs, add `-DRDB_ENABLE_MPI=ON`, which uses cray-mpich through
the `ftn` wrapper. NetCDF stays serial. Read the README's MPI section about
pinning `CUDA_VISIBLE_DEVICES` before `MPI_Init`.

### Cray notes

- The script detects a Cray PE from `$CRAYPE_VERSION`. With no `--fc`, it
  defaults to `ftn`/`cc`. `CRAYPE_LINK_TYPE=static` turns on `--static`.
- `--static` builds `.a` libraries only. Then pass netcdf-c's dependency list
  to the link, as the script prints:
  `-DCMAKE_Fortran_STANDARD_LIBRARIES="$(<prefix>/bin/nc-config --libs --static)"`.
  The Cray wrappers add HDF5 themselves while `cray-hdf5` is loaded.
- Use `cray-hdf5`, not `cray-hdf5-parallel`, unless you mean `--mpi`.

## Options

`tools/build_netcdf.sh --help` has the full list. The main ones:

| Option | Meaning |
|---|---|
| `--fc`, `--cc` | compilers (default `$FC`/`$CC`; `ftn`/`cc` on a Cray PE) |
| `--hdf5 <prefix>` | the HDF5 to use (default `$HDF5_DIR`, `$HDF5_ROOT`, `h5cc`, pkg-config) |
| `--prefix`, `--jobs` | install location, make parallelism |
| `--mpi` | parallel HDF5: `mpicc`/`mpifort` (`cc`/`ftn` on Cray), parallel netcdf-c |
| `--static` | static libraries only |
| `--cflags`, `--fflags` | extra flags, e.g. `-tp=zen3` for NVHPC |
| `--download-only`, `--download-dir`, `--offline` | split the fetch from the build |
| `--check` | also run the upstream netcdf-c/-fortran test suites |
| `--dry-run` | find everything and print every command, but build nothing |
| `--fortran-only --netcdf-c <prefix>` | the `build_netcdf_fortran.sh` mode (wrapper only) |

## What has been verified, and where

**Verified on the dev box** (x86-64 Ubuntu 24.04, V100, 2026-09-24):

| Toolchain | HDF5 | Result |
|---|---|---|
| gfortran 15.1 + gcc 15.1 | Debian HDF5 1.10.10, found via `h5cc -show` (split `include/hdf5/serial` layout) | netcdf-c 4.9.3 + netcdf-fortran 4.6.2 built and smoke-tested; Roundabout found the bootstrapped config; `ctest -R rdb` 226/226 |
| nvfortran + nvc 26.5 | module HDF5 1.14.2, found via `$HDF5_DIR` | built and smoke-tested; Roundabout GPU build (`-DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc70`) found the bootstrapped config; `ctest -R rdb` 226/226 on a V100 |
| gfortran 15.1, `--check` | Debian HDF5 1.10.10 | upstream suites: netcdf-c 186/186, netcdf-fortran 38/38 |
| gfortran 15.1, `--static` | Debian HDF5 1.10.10 | static build smoke-tested; a CMake `find_package(netCDF-Fortran)` consumer linked with the printed `nc-config --libs --static` line |
| gfortran 15.1, netcdf-c 4.10.1 + netcdf-fortran 4.6.4 | Debian HDF5 1.10.10 | built and smoke-tested |
| `--download-only`, then `--offline` | — | the offline build ran with `curl`/`wget` replaced by stubs that fail if called |
| `build_netcdf_fortran.sh` (the `--fortran-only` wrapper) | — | gfortran netcdf-fortran on the nvc-built netcdf-c, smoke-tested |

No NetCDF module was loaded for any of these runs, and `rdb` resolved every
NetCDF and HDF5 library through RPATH with an empty environment (`ldd`).

**Not verified:**

- **Anything on Perlmutter.** That covers the module names and versions
  above, the `cc`/`ftn` wrappers, the `cray-hdf5` `$HDF5_DIR` layout,
  `CRAYPE_LINK_TYPE=static`, and the `-tp=zen3` builds. The recipe follows
  NERSC's documented environment, but nobody has run it there.
- **`--mpi`.** The dev box has no parallel HDF5. Only the refusals are tested
  (a serial HDF5 with `--mpi`).
- **ifx / flang.** They pair with `icx` / `clang`, but those builds have not
  been run.
