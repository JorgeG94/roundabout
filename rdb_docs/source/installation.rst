.. _installation:

------------
Installation
------------

.. contents::
   :local:


Obtaining Roundabout
====================

Clone the git repository:

.. code-block:: bash

   git clone https://github.com/JorgeG94/roundabout.git
   cd roundabout


What you need
=============

**A Fortran compiler.** Roundabout is written against Fortran 2018, and it
leans hard on ``do concurrent`` *with locality specifiers* — ``local(...)``,
``local_init(...)``, ``shared(...)``. Those clauses are how a kernel declares
its per-iteration private state so that one source compiles to both a CPU
loop and a GPU kernel. They are load-bearing, not stylistic, and they
constrain which compilers work:

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Compiler
     - Minimum
     - Notes
   * - NVIDIA ``nvfortran`` (NVHPC)
     - —
     - The only GPU-offload toolchain. Also the best-tested CPU multicore
       path.
   * - GNU ``gfortran``
     - **15.0**
     - CPU only. GCC implements the ``do concurrent`` locality clauses from
       GCC 15; CMake refuses to configure below that.
   * - Intel ``ifx``
     - —
     - CPU, and GPU via the OpenMP-target backend.

.. warning::

   The gfortran floor is enforced at configure time, not discovered at
   compile time:

   .. code-block:: text

      gfortran 13.2.0 is too old to build Roundabout.

        Required: GNU Fortran >= 15.0
        Found:    13.2.0

      WHY: the ocean dynamical core is written with F2018 `do concurrent`
      locality specifiers -- `local(...)`, `local_init(...)`, `shared(...)`
      -- which are how the kernels declare per-iteration private state so
      they parallelise on CPU and GPU from the SAME source. GCC implements
      those clauses only from GCC 15. Older gfortran accepts `do concurrent`
      itself but rejects the locality clauses, so the build would fail with
      a wall of syntax errors inside the kernels instead of this message.

   This is deliberate. A GCC 14 build would otherwise fail deep inside the
   kernels with hundreds of syntax errors and no indication of the cause.

**CMake 3.25 or newer**, and the **Unix Makefiles** generator. Ninja is not
supported: the analytical-test helper library races under Ninja and fails
with ``Unable to open MODULE file rdb_constants.mod``. The shipped CMake
presets pin the generator for you.

**NetCDF-Fortran**, unless you configure with ``-DRDB_ENABLE_NETCDF=OFF``.
This is the one dependency you must supply yourself. CMake looks for a
CMake config package first — honouring the ``NETCDFF_DIR`` and
``NETCDF_FORTRAN_DIR`` environment variables as hints — and falls back to
pkg-config's ``netcdf-fortran``.

**Three Fortran libraries are fetched for you** by CMake's ``FetchContent``
at configure time; you do not install them:

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Library
     - Source
     - Provides
   * - ``pic``
     - ``github.com/JorgeG94/pic``
     - Base types, string helpers, timers, the logger
   * - ``pic-mpi``
     - ``github.com/JorgeG94/pic-mpi``
     - The MPI wrapper layer. Always a dependency — ``RDB_ENABLE_MPI=OFF``
       selects its *serial backend* rather than removing it.
   * - ``test-drive``
     - ``github.com/fortran-lang/test-drive``
     - The unit-test framework

.. note::

   All three are pinned — ``pic`` and ``pic-mpi`` to release tags,
   ``test-drive`` to a commit — in ``cmake/Findpic.cmake``,
   ``cmake/Findpic-mpi.cmake`` and ``cmake/Findtest-drive.cmake``. A clone
   builds the same dependency versions CI does. To develop against a local
   checkout instead, use the stock FetchContent override, e.g.
   ``-DFETCHCONTENT_SOURCE_DIR_PIC=/path/to/local/pic``.


Building
========

The short version:

.. code-block:: bash

   cmake -B build -S .
   cmake --build build -j
   cd build && ctest --output-on-failure

That gives you a **CPU, single-rank, double-precision, NetCDF-enabled**
build — GPU offload is opt-in.

Configuring prints a summary block; check it before you trust a build:

.. code-block:: text

   -- === Roundabout Build Configuration ===
   -- Compiler:       NVHPC (/opt/nvidia/hpc_sdk/.../nvfortran)
   -- Build type:     Release
   -- Release flags:  -O3 -fast
   -- GPU offload:    ON
   -- Double prec:    ON
   -- NetCDF I/O:     ON
   -- Tests:          ON
   -- Benchmarks:     OFF
   -- =================================

The ``GPU offload:`` line is the one that catches the classic mistake of
running a bare ``cmake -B build -S .`` in a tree you thought was a GPU
build. If it says ``OFF``, it is off.

Build type defaults to ``Release`` when you do not set one.


The three toolchains
====================

NVHPC — GPU
-----------

The production configuration.

.. code-block:: bash

   cmake -B build_cc70 -S . \
         -DCMAKE_Fortran_COMPILER=nvfortran \
         -DRDB_ENABLE_GPU=ON \
         -DRDB_GPU_ARCH=cc70
   cmake --build build_cc70 -j

``RDB_GPU_ARCH`` must be an NVIDIA ``ccNN`` string on this path —
``cc70`` (Volta, the default), ``cc80`` (Ampere), ``cc90`` (Hopper).
Passing anything that does not match ``ccNN`` is a configure-time fatal
error, not a silent fallback.

The compiler flags this selects are

.. code-block:: text

   -stdpar=gpu -acc=gpu -gpu=<arch>,mem:separate

``mem:separate`` is the single most consequential flag in the project.
There is **no** managed or unified memory: the device has its own address
space, and nothing moves between host and device unless the code says so.
Every array a kernel touches must have been explicitly mapped. This shapes
how kernels and tests are written; see :ref:`extending`.

NVHPC — CPU multicore
---------------------

The same source, the same ``do concurrent`` loops, threaded across CPU
cores:

.. code-block:: bash

   cmake -B build_multicore -S . \
         -DCMAKE_Fortran_COMPILER=nvfortran \
         -DRDB_ENABLE_THREADS=ON
   cmake --build build_multicore -j

This selects ``-stdpar=multicore -acc=multicore``. It is useful for
debugging kernels without a GPU, but note the caveat in
:ref:`extending`: a green multicore run proves nothing about device data
motion.

gfortran — CPU
--------------

The everyday development build, and the one most contributors will use:

.. code-block:: bash

   cmake -B build_serial -S . \
         -DCMAKE_Fortran_COMPILER=gfortran
   cmake --build build_serial -j

Requesting ``RDB_ENABLE_GPU=ON`` with gfortran is a fatal configure error —
gfortran has no ``do concurrent`` GPU offload path, and the build refuses
to quietly give you a CPU binary instead.

ifx — CPU and Intel GPU
-----------------------

.. code-block:: bash

   cmake -B build_ifx -S . -DCMAKE_Fortran_COMPILER=ifx

``-qopenmp`` is added unconditionally on this compiler and is **not
optional**: without it, ifx at ``-O2`` and above silently miscompiles
multi-index ``do concurrent`` bodies that write ``intent(out)`` scalars
through pure-subroutine calls — the writes become zero, with no error. The
compiler's own warning says as much (*"Locality information is ignored
without one of these command line qualifiers"*).

For Intel GPU offload you must also switch the parallel backend:

.. code-block:: bash

   cmake -B build_ifx_gpu -S . \
         -DCMAKE_Fortran_COMPILER=ifx \
         -DRDB_ENABLE_GPU=ON \
         -DRDB_PARALLEL_BACKEND=openmp

which adds ``-fopenmp-targets=spir64 -fopenmp-target-do-concurrent``.
Asking for ``RDB_ENABLE_GPU=ON`` with the default OpenACC backend on ifx is
a fatal error naming the fix.


CMake options
=============

The options you are most likely to touch:

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Option
     - Default
     - What it does
   * - ``RDB_ENABLE_GPU``
     - ``OFF``
     - GPU offload (``do concurrent`` + OpenACC). Opt-in. Fatal error on a
       compiler with no GPU path — never a silent CPU downgrade.
   * - ``RDB_ENABLE_THREADS``
     - ``OFF``
     - CPU multicore threading.
   * - ``RDB_ENABLE_DOUBLE``
     - ``ON``
     - Working precision ``real64``. All shipped science runs at this
       setting.
   * - ``RDB_ENABLE_MPI``
     - ``OFF``
     - Link a real MPI and build multi-rank. ``OFF`` still uses pic-mpi,
       via its serial backend.
   * - ``RDB_CUDA_AWARE_MPI``
     - ``OFF``
     - GPU-direct halo exchange. See the multi-GPU warning below.
   * - ``RDB_GPU_ARCH``
     - ``cc70``
     - NVIDIA ``ccNN`` under NVHPC. (Under LLVM flang it defaults to
       ``gfx90a`` and takes AMD ``gfxNNN`` strings instead.)
   * - ``RDB_ENABLE_NETCDF``
     - ``ON``
     - The whole NetCDF-backed I/O subsystem — driver, output, restart,
       forcing. Turning it off also drops the ``rdb`` executable; kernels
       and benchmarks still build. Use it for portability testing where
       NetCDF is inconvenient.
   * - ``RDB_ENABLE_TESTING``
     - ``ON``
     - Build the test suite.
   * - ``RDB_ENABLE_NVTX``
     - ``OFF``
     - NVTX annotations for Nsight Systems.
   * - ``RDB_VERBOSE_COMPILE``
     - ``OFF``
     - ``-Minfo=all`` / ``-fopt-info-loop``. Useful when you want to know
       whether a loop actually offloaded.
   * - ``RDB_ENABLE_COVERAGE``
     - ``OFF``
     - gcov instrumentation, gfortran only. Adds a ``coverage`` build
       target that renders ``build/coverage_report/index.html``.
   * - ``RDB_NZ_STACK_MAX``
     - ``128``
     - Maximum vertical layers for per-column stack workspaces. Lower it if
       a device link hits stack-frame limits.
   * - ``RDB_EXTRA_FORTRAN_FLAGS``
     - *(empty)*
     - Extra Fortran flags, appended after every flag the project sets and
       applied to the fetched dependencies too, so the whole build agrees.
       The case that motivates it is CPU targeting: **nvfortran targets the
       build host's CPU by default** — nothing in
       ``cmake/compiler_flags.cmake`` sets ``-tp``, and neither ``-fast`` nor
       the build type changes that. Build on a login node, run on a compute
       node with a different CPU, and you get ``SIGILL``. Pin it with
       ``-DRDB_EXTRA_FORTRAN_FLAGS=-tp=haswell`` on NVHPC, or
       ``-march=x86-64-v3`` on GNU / ifx / flang. The value prints on the
       configure summary as ``Extra flags:``.

The full list, including the advanced knobs, is at the top of
``CMakeLists.txt``.


What gets built
===============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Artifact
     - Notes
   * - ``build/rdb``
     - The solver. Built only when ``RDB_ENABLE_NETCDF=ON``.
   * - ``build/libcore_rdb.a``
     - The static core library.
   * - ``build/librdb_core.so``
     - The shared core, built when tests link shared (the default) or
       ``RDB_BUILD_SHARED=ON``.
   * - ``build/rdb_nml_doc``, ``build/rdb_nml_json``
     - Namelist-schema dumpers. ``rdb_nml_doc`` is what regenerates
       ``docs/generated_nml_knobs.md``.
   * - ``build/modules/``
     - Fortran ``.mod`` files.
   * - ``build/tests/``
     - Test binaries.


Environment setup
=================

Roundabout ships no machine-specific environment scripts. Load a toolchain
however your site does it — ``module load``, Spack, conda, or a hand-rolled
script — and then configure and build. What the build needs on ``PATH``:

* a Fortran compiler: gfortran **15+**, NVHPC (``nvfortran``) for the GPU
  build, ifx, or flang;
* CMake 3.25+;
* netcdf-fortran and its netcdf-c / hdf5;
* MPI only if you configure ``RDB_ENABLE_MPI=ON``.

A module-based site typically looks like:

.. code-block:: bash

   module load gcc netcdf-c netcdf-fortran hdf5    # names are site-specific
   cmake -B build -S . && cmake --build build

.. warning::

   Set up **exactly one** toolchain per shell. Stacking a gfortran environment
   and an NVHPC one puts two incompatible NetCDF builds on the link line and
   the NetCDF tests fail in confusing ways. Use a fresh shell to switch.

.. tip::

   The project's own tooling (``validate.sh``,
   ``tests/regression/run_all.py``) assumes the toolchain is already loaded.
   If you keep a setup script, point ``RDB_ENV_SCRIPT`` at it — or
   ``RDB_CPU_ENV_SCRIPT`` / ``RDB_GPU_ENV_SCRIPT`` when the CPU and GPU
   toolchains need different scripts — and the tooling will source it in each
   stage's own subshell. Unset, they run in the environment they inherit.

.. _getting-netcdf-fortran:

Getting netcdf-fortran
======================

This is the one dependency you have to supply, and it is worth two minutes
of explanation because the obvious approach — "install the whole NetCDF
stack with my package manager" — is both the slowest option and the one
that most often fails.

**Only one library in the tree is compiler-coupled.** Fortran ``.mod``
files and module symbol manglings are compiler-specific, so a
netcdf-fortran built by gfortran cannot be consumed by nvfortran, ifx or
flang. You need one netcdf-fortran per compiler.

**Everything below it is C, and one C build serves every compiler.**
netcdf-c, HDF5 and zlib are C libraries following the platform ABI. You do
*not* need a per-compiler netcdf-c. This is verified, not assumed: the full
187-test suite passes with an ``nvfortran`` solver linked against a
**gfortran-built** netcdf-c.

So the fast path is to take a prebuilt netcdf-c from wherever is convenient
and build only the Fortran wrapper on top:

.. code-block:: bash

   # any netcdf-c will do -- distro, module tree, conda, Spack
   sudo apt install libnetcdf-dev          # or dnf install netcdf-devel

   module load nvhpc                       # whatever compiler you want
   tools/build_netcdf_fortran.sh --fc nvfortran

That takes about 30 seconds. The script downloads a pinned, checksummed
netcdf-fortran release, builds it against the netcdf-c it finds, runs a
smoke test that writes and reads a file, and prints the exact ``cmake``
line to configure Roundabout against the result.

It locates netcdf-c in this order: ``--netcdf-c <prefix>``, then
``$NETCDFC_DIR``, then ``nc-config --prefix`` from ``PATH``. Pass
``--prefix`` to choose where it installs; the default embeds the compiler
name so two compilers never collide.

Two things Roundabout specifically does **not** need, which is most of what
makes a NetCDF stack slow to build:

* **Parallel NetCDF.** The I/O is per-rank serial ``nf90_create`` with an
  offline merge (``tools/merge_output.py``). A serial netcdf-c is enough
  even for an MPI build, and MPI-enabled HDF5 is the most fragile part of
  that stack.
* **HDF5's Fortran bindings.** The source imports ``netcdf`` and nothing
  else, so ``hdf5 +fortran`` is pure build time.

Other routes, if the script does not suit:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Route
     - Notes
   * - Site modules
     - Nothing to build. Most HPC sites ship netcdf-fortran per compiler.
   * - conda-forge
     - ``conda install -c conda-forge netcdf-fortran``. One minute, no
       build — but **gfortran builds only**, so it is no help for
       nvfortran, ifx or flang.
   * - ``environments/spack.yaml``
     - Builds the *C* layer (zlib, HDF5, netcdf-c) for machines with no
       usable system one. Deliberately pins no compiler and builds no
       netcdf-fortran: that half is compiler-agnostic, so you build it once
       and reuse it, and a stock build gets binary-cache hits that
       ``%nvhpc`` never would. Put the script on top.
   * - ``-DRDB_ENABLE_NETCDF=OFF``
     - No NetCDF at all. Kernels, unit tests and benchmarks build; the
       ``rdb`` executable does not. This is the portability-testing
       configuration, not a way to run simulations.

.. note::

   Prefer a system or module netcdf-c over conda's when you have the
   choice. Conda's drags its own libcurl and OpenSSL along with RPATHs
   pointing into the conda prefix, which gets messy in a mixed link.

.. warning::

   An **MPI** build has a *second* compiler-coupled dependency: pic-mpi
   uses ``mpi_f08``, so the MPI library's Fortran bindings must match your
   compiler too. A serial build has only netcdf-fortran to worry about.

   Distro packages will not do here. Ubuntu's ``libopenmpi-dev`` is built
   against the distro gfortran (13 on 24.04), below this project's floor of
   15, and the mismatch is fatal::

      Fatal Error: Cannot read module file 'mpi_f08.mod' opened at (1),
      because it was created by a different version of GNU Fortran

   Three things that do work:

   * **Your site's MPI module**, built per compiler — how every HPC site
     already ships it, and the right answer on a cluster.
   * **conda-forge**, for gfortran: its ``openmpi`` is a gfortran-15 build,
     so ``gfortran`` and ``openmpi`` resolved in one solve match. Nothing in
     conda's metadata encodes ``.mod`` compatibility, so check the pair
     rather than assume it — the ``gha3mi/setup-fortran-conda`` action does
     exactly that, compiling each binding and running a two-rank job before
     it reports success, which is what the MPI CI leg relies on.
   * **NVHPC**, which bundles an HPC-X OpenMPI whose Fortran modules are
     nvfortran-built, so for the GPU toolchain the MPI comes with the
     compiler and there is nothing to match.

There is also ``environments/roundabout_env_3.13.yml``, a conda environment for
the Python post-processing and test tooling. It deliberately **excludes**
``netcdf4`` and ``mpi4py``: loading conda's HDF5 into the same process as
``librdb_core.so`` produces ABI mismatches and MPI symbol clashes.


Running the tests
=================

.. code-block:: bash

   cd build
   ctest --output-on-failure

Every Roundabout test is named with an ``rdb_`` prefix, so

.. code-block:: bash

   ctest -R rdb --output-on-failure

runs the project's own suite and skips the dependency self-tests that
``pic`` and ``test-drive`` register. Tests also carry regime labels, so
``ctest -L ocean`` and ``ctest -L core`` select subsets.

.. warning::

   **Never pass ``-j N`` to ctest on a GPU build.** Every worker shares the
   one GPU, and parallel execution produces spurious failures and hangs.
   Building with ``-j`` is fine; it is only ``ctest -j`` that breaks.


Multi-GPU on one node
=====================

If you build with ``RDB_CUDA_AWARE_MPI=ON`` and run several ranks on one
node, **each rank must see only its own GPU, and that must be arranged
before ``MPI_Init``**.

The reason: UCX creates a CUDA primary context during ``MPI_Init``, which
runs *before* the solver calls ``acc_set_device_num``. If every GPU is
visible to every rank, every rank lands a context on device 0 in addition
to its intended device — ``nvidia-smi`` shows ``{0}, {0,1}, {0,2}, {0,3}``
— while the solver's own binding diagnostic still looks perfectly correct.

The fix is to set ``CUDA_VISIBLE_DEVICES`` from the launcher's local-rank
variable before the first CUDA call:

.. code-block:: bash

   export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK

With that pinning, the solver's ``mod(node_rank, n_devices)`` clamp binds
the single visible device and a plain ``mpirun -np N ./rdb run.nml``
works. Host-staged MPI (``RDB_CUDA_AWARE_MPI=OFF``) never touches CUDA
during ``MPI_Init`` and does not hit this.
