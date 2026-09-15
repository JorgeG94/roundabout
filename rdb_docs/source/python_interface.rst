.. _python_interface:

----------------
Python interface
----------------

.. contents::
   :local:

Roundabout can be driven from Python. Fortran still does the arithmetic and
still owns the device; Python builds the configuration, creates the solver,
advances it, and reads and writes its state arrays — in-process, through
``ctypes``, against the same ``librdb_core`` the executable links.

This is a **runtime driver**, not a namelist generator. ``model.h`` is not a
copy of a file that was written earlier; it is a window onto the live
``h_layer`` array of a solver that is halfway through a run. The reason to
reach for it is *decisions in between steps* — a namelist states one run from
beginning to end, a script can look at the mixed-layer depth at day 30 and
decide what to do about it. Diagnosing a run interactively, sweeping a
parameter without writing twenty ``.nml`` files, checking an analytical
solution cell by cell, driving the model from a notebook: those are what this
is for. A production run that is fully described before it starts is still
better expressed as a namelist and handed to the ``rdb`` binary.

The package is ``python/rdb`` in the repository. It has **no dependencies at
all** — stdlib plus ``ctypes``, no numpy, nothing to install. numpy is used
if it happens to be importable, and is never required.


The whole interface, at a glance
================================

One object does the work, and everything else describes it:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Object
     - What it is
   * - ``rdb.Model``
     - A live solver handle. Created from a namelist string, a
       :class:`~rdb.Config`, a namelist file, or a set of curated physics
       objects. Steps, exposes state, and must be closed.
   * - ``rdb.Field``
     - A window onto one state array — layer thickness, a velocity
       component, a surface flux. Knows its shape, its units and where it
       sits on the C-grid.
   * - ``rdb.DerivedField``
     - A field that is *computed* rather than stored: tracer concentration,
       which the solver keeps as ``h*Tr``.
   * - ``rdb.DiagnosticField``
     - A read-only view onto the diagnostic manager's own output buffer, for
       a diagnostic that is registered on this run.
   * - ``rdb.Loc``
     - ``Loc.CENTER`` / ``Loc.FACE`` — the stagger point of each axis of a
       field.
   * - ``rdb.Config``
     - Every namelist group and every knob, typed and validated on
       assignment. Generated from the solver's own schema.
   * - ``rdb.Diagnostics``, ``rdb.Duration``, ``rdb.Restart``
     - Curated writers over ``&ocean_diag_nml``, ``&time_nml`` and the
       restart half of ``&output_nml``.
   * - ``rdb.RdbError``
     - The base of a typed exception hierarchy. Every failure the Fortran can
       return arrives as a specific subclass carrying the Fortran's own
       sentence.

.. code-block:: python

   import rdb

   with rdb.Model(open("case.nml").read()) as model:
       model.step(10)
       print(model.time, model.total_mass)
       print(model.temperature[0, 0, -1])     # surface layer, degC

That is the whole shape of it: create, step, read, and let the context
manager close the handle.


Building and importing
======================

The interface loads ``librdb_core.so``. It is built by default — the test
suite links against it, so an ordinary configure produces one:

.. code-block:: bash

   cmake -B build -S . && cmake --build build
   ls build/librdb_core.so

If you want the shared core explicitly, and nothing else:

.. code-block:: bash

   cmake -B build_shared -S . -DRDB_BUILD_SHARED=ON
   cmake --build build_shared

There is no ``install()`` rule and nothing to ``pip install``. Point
``PYTHONPATH`` at ``python/`` and import:

.. code-block:: bash

   export PYTHONPATH=$PWD/python
   python3 -c "import rdb; print(rdb.find_library())"

The toolchain environment that built the library has to be loaded in the shell
that imports it, or ``dlopen`` will not find NetCDF and HDF5 — load it exactly
as you did for the build (``module load``, Spack, conda; see
:doc:`installation`), and load only that one.

Finding the library
-------------------

:func:`rdb.find_library` runs a fixed search, first hit wins:

#. ``$RDB_LIB`` — an explicit path to the ``.so``.
#. ``<repo root>/build_shared/librdb_core.so``.
#. Any other ``build*/librdb_core.so`` or ``build*/*/librdb_core.so`` under
   the repository root, sorted.
#. The normal dynamic-linker search path (``LD_LIBRARY_PATH``, the ldconfig
   cache, rpath) via ``ctypes.util.find_library("rdb_core")``.

The repository root is derived from the package's own location
(``python/rdb/_ffi.py``, two directories up), so an in-tree checkout with any
``build*`` directory usually needs no configuration at all. A build somewhere
else is named outright:

.. code-block:: bash

   export RDB_LIB=/scratch/me/roundabout/build_gpu/librdb_core.so

Failure is explicit. ``LibraryNotFoundError`` lists every path that was
tried and tells you the CMake line that would produce one — it never falls
back to a stale library it happened to find.

:func:`rdb.get_lib` returns the process-wide singleton that wraps the
``CDLL`` handle. The library is loaded lazily, on first use, so importing
``rdb`` touches nothing. Ordinary scripts never need ``get_lib``; it exists
for ``.path``, and for the rare case of forcing the load early to find out
whether the build is usable at all.

Working precision
-----------------

The library reports the precision it was compiled with:

.. code-block:: python

   >>> rdb.working_precision()
   8

Bytes per float: ``8`` for the default ``RDB_ENABLE_DOUBLE=ON`` build, ``4``
for a single-precision one. Everything downstream follows it — the
``typestr`` in a field's ``__array_interface__``, the dtype numpy sees, the
size of a read. It is worth checking once in any script that compares
against a saved reference, because a single-precision build will not
reproduce a double-precision one and nothing else will tell you.

numpy is optional
-----------------

The package imports no third-party module anywhere, including in the array
plumbing. A :class:`~rdb.Field` exposes ``__array_interface__``, which numpy
consumes zero-copy when it is present:

.. code-block:: python

   import numpy as np
   arr = np.asarray(model.h)        # a view, not a copy

Without numpy, ``field[...]`` and ``field.copy()`` return nested Python
lists, and a full integer index returns a float. What is *not* available
without numpy is general slicing — ``field[:, :, -1]`` raises
:class:`~rdb.RdbUnsupportedError` with that advice in the message. When numpy
is importable the same expression works, because the field quietly hands
numpy its own buffer and lets numpy do the slicing.


A first run
===========

The shortest complete thing. A namelist string, ten steps, one number out:

.. code-block:: python

   import rdb

   nml = """
   &sim_nml            sim_type = "ocean" /
   &grid_nml           nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /
   &nonhydrostatic_nml nz_layers = 3 /
   &time_nml           t_end = 86400.0, dt_fixed = 300.0 /
   &ocean_topo_nml     topo_config = "flat", max_depth = 200.0 /
   &ocean_bt_nml       auto_n_inner = .true. /
   &tracer_nml         initial_salinity = 35.0, initial_temperature = 12.0 /
   &ocean_diag_nml     enabled = .false. /
   """

   with rdb.Model(nml) as model:
       print(model.grid_info)
       model.step(10)
       print(model.time, model.step_count)
       print(model.h.shape, model.temperature[0, 0, -1])

.. code-block:: text

   {'nx': 8, 'ny': 6, 'nz': 3, 'nghost': 3}
   3000.0 10
   (8, 6, 3) 12.0

Three things in that output are worth naming now.

``nghost`` is **derived, not asked for**. The namelist above never mentions a
halo width; the solver worked out that this configuration needs three ghost
cells and said so.

``model.h.shape`` is ``(8, 6, 3)`` — the **interior**, not the
ghost-inclusive ``(14, 12, 3)`` the Fortran array actually is. Ghosts are
reachable, by name, and never by default.

``model.temperature[0, 0, -1]`` is degrees Celsius at the **surface**. The
layer stack is bottom-up: ``k=0`` in Python (``k=1`` in Fortran) is the bed
and the last index is the surface layer. This is the single convention most
worth keeping in your head; see :ref:`vertical_coordinates`.

**The handle must be closed.** One live ocean solver exists per process — the
halo module, the profiler and the diagnostic fills all hold module-scope
state, so a second live handle is refused rather than allowed to corrupt the
first. ``with rdb.Model(...)`` closes on the way out, including when the body
raises. Without the ``with``, call ``model.close()`` yourself; it is
idempotent, so closing twice is a clean no-op, and ``__del__`` calls it as a
backstop. Trying to create a second model while one is open raises
:class:`~rdb.AlreadyExistsError`.


Configuring a model
===================

There are two ways to say what a run is, and they compose in one script.

The **curated objects** are named physics: ``RectilinearGrid``,
``Smagorinsky``, ``EPBL``, ``Flather``. They are the readable surface, and
they check things that cross object boundaries.

The **generated config** is every namelist knob the schema carries — 57
groups, 613 keys — typed and validated, reachable by name. It is the escape
hatch, and nothing is missing from it.

The curated objects
-------------------

Each curated object is a *writer*: its whole effect is to set named knobs on
a :class:`~rdb.Config`. Passing them to ``Model(...)`` composes them into one
configuration:

.. code-block:: python

   import rdb
   from rdb.barotropic import BarotropicSolver
   from rdb.closures import ConstantViscosity, KPP, PacanowskiPhilander
   from rdb.coriolis import BetaPlane
   from rdb.diagnostics import Diagnostics, Duration
   from rdb.forcing import ConstantWind, UniformInitialCondition
   from rdb.grids import FlatBottom, RectilinearGrid, Topology
   from rdb.vcoord import ZStarSigma

   model = rdb.Model(
       grid=RectilinearGrid(size=(20, 16, 4), extent=(200e3, 160e3),
                            topology=(Topology.PERIODIC, Topology.BOUNDED)),
       bathymetry=FlatBottom(1000.0),
       timestep=300.0,
       coriolis=BetaPlane(f0=-1.0e-4, beta=1.5e-11),
       vcoord=ZStarSigma(),
       closures=[PacanowskiPhilander(), KPP(), ConstantViscosity(nu_h=200.0)],
       forcing=[ConstantWind(tau_x=0.05)],
       initial_condition=UniformInitialCondition(temperature=15.0,
                                                 salinity=35.0),
       barotropic=BarotropicSolver(auto=True),
       diagnostics=Diagnostics(enabled=False),
       duration=Duration(t_end=10.0, time_unit="day"),
   )

   with model:
       model.run(steps=6, every=3,
                 callback=lambda m: print(m.step_count, m.kinetic_energy))

The keyword names are fixed and exhaustive — ``grid``, ``bathymetry``,
``timestep``, ``coriolis``, ``momentum_advection``, ``buoyancy``,
``pressure``, ``vcoord``, ``closures``, ``boundaries``, ``forcing``,
``tracers``, ``initial_condition``, ``barotropic``, ``continuity``,
``diagnostics``, ``restart``, ``duration`` — so a misspelling is a
``TypeError`` at the call rather than a setting that silently did nothing.

The modules, and what each covers:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Module
     - Contents
   * - ``rdb.grids``
     - ``RectilinearGrid``, ``LatitudeLongitudeGrid``, ``TripolarGrid``,
       ``SupergridGrid``; ``Topology``; ``required_halo``; and the bottom —
       ``Bathymetry``, ``FlatBottom``, ``Spoon``, ``Seamount``, ``Island``,
       ``DoubleDrake``, ``Neverworld2``, ``BathymetryFile``.
   * - ``rdb.closures``
     - Vertical: ``PacanowskiPhilander``, ``KPP``, ``EPBL``, ``Langmuir``,
       ``KappaShear``, ``TidalMixing``, ``ConvectiveAdjustment``,
       ``DoubleDiffusion``, ``BryanLewisBackground``, ``HenyeyBackground``,
       ``ConstantBackground``, ``ImplicitVerticalFriction``,
       ``DiffusivityLimits``, ``ThermoSubcycling``. Lateral:
       ``Smagorinsky``, ``Leith``, ``ConstantViscosity``,
       ``HorizontalDiffusivity``, ``Anisotropy``, ``MEKE``,
       ``GentMcWilliams``, ``Redi``, ``VariableMixing``,
       ``MixedLayerRestratification``, ``Bodner``. Bottom drag:
       ``QuadraticDrag``, ``LinearDrag``, ``ChannelDrag``.
   * - ``rdb.forcing``
     - ``ConstantWind``, ``DoubleGyreWind``, ``Neverworld2Wind``,
       ``Thermodynamics``, ``SurfaceHeatFlux``, ``SurfaceSaltFlux``,
       ``ShortwavePenetration``, ``BuoyancyRestoring``, ``Geothermal``,
       ``AtmosphericPressure``, ``SurfaceFluxComponents``, ``FileForcing``,
       ``ScalarSAL``, ``EquilibriumTide``, and the initial-condition
       overlays (``UniformInitialCondition``,
       ``StratifiedInitialCondition``, ``ZLevelInitialCondition``,
       ``LayerDensities``, ``EadyIC``, ``BaroclinicJetIC``, …).
   * - ``rdb.boundaries``
     - ``Wall``, ``Flather``, ``Clamped``, ``Chapman``, ``Sponge``,
       ``TidalBoundary``, ``Radiation``, ``Nudging``, and ``Boundaries``,
       which holds one per edge.
   * - ``rdb.vcoord``
     - ``Sigma``, ``ZSigma``, ``ZStar``, ``ZStarSigma``, ``ZStarFull``,
       ``Isopycnal``, ``Hycom``, ``Lagrangian``, ``ALE``.
   * - ``rdb.pressure``
     - ``Montgomery``, ``FiniteVolumeLite``, ``FiniteVolumeWright``,
       ``FiniteVolumeMOM6``, ``ReducedGravity``, ``VelocityTruncation``.
   * - ``rdb.coriolis``
     - ``FPlane``, ``BetaPlane``, ``SphericalCoriolis`` (the rotation), and
       ``Sadourny`` + ``WENO`` (the discretisation).
   * - ``rdb.buoyancy``
     - ``LinearEquationOfState``, ``Wright1997``, ``RoquetSPV``, ``TEOS10``.
   * - ``rdb.tracers``
     - ``IdealAge``, ``PseudoSalt``, ``TracerBounds``.
   * - ``rdb.barotropic``
     - ``BarotropicSolver``, ``Continuity``, ``WaveDrag``,
       ``PorousBarriers``, ``WetDry``.
   * - ``rdb.diagnostics``
     - ``Diagnostics``, ``Duration``, ``Restart``.

What the curated layer adds over the raw knobs is the checks that need two
objects in view at once, which no namelist can express and no single knob can
catch. KPP and EPBL are mutually exclusive, so asking for both raises
:class:`~rdb.ConfigConflictError`. ``HenyeyBackground`` on a
``RectilinearGrid`` raises, because a Cartesian grid has no latitude and
every column would read as equatorial. ``MEKE(backscatter=True)`` without a
non-zero biharmonic backstop somewhere in the same ``closures=[...]`` raises,
naming the three objects that could provide one — it does not invent a
dissipation scale, because picking one is guessing at physics.

None of these replace the Fortran's own ``validate_config``. Every one of
them is *also* enforced there; the Python check simply names the problem one
call earlier, with both objects in view.

The generated config
--------------------

:class:`rdb.Config` is the other half. Every namelist group is an attribute,
every knob is an attribute of that:

.. code-block:: python

   cfg = rdb.Config()
   cfg.sim.sim_type = "ocean"
   cfg.grid.nx = 8
   cfg.grid.ny = 6
   cfg.grid.dx = 2000.0
   cfg.grid.dy = 2000.0
   cfg.nonhydrostatic.nz_layers = 3
   cfg.time.dt_fixed = 300.0
   cfg.ocean_topo.topo_config = "flat"
   cfg.ocean_topo.max_depth = 200.0
   cfg.ocean_bt.auto_n_inner = True
   cfg.ocean_diag.enabled = False

   with rdb.Model(cfg) as model:
       ...

``python/rdb/_config_generated.py`` is **generated**, not written. It is
produced by ``tools/gen_python_config.py`` from the solver's live
``nml_schema_t`` — the same schema object that parses a namelist file, that
the strict parser checks against, and that ``rdb_nml_doc`` dumps to
``docs/generated_nml_knobs.md``. That is why this page does not list the
knobs: any list here would be a second copy with its own decay rate. Regenerate
it when the schema changes, and read the generated file (or
``docs/generated_nml_knobs.md``) for what a knob means.

Validation on assignment mirrors the Fortran parser arm for arm — range
checks on reals and ints, case-insensitive membership for enums stored back
in the canonical spelling, a length cap on strings:

.. code-block:: python

   >>> cfg.grid.nx = 0
   rdb._errors.ConfigParseError: key 'nx' = 0 below min 1

   >>> cfg.ocean_coriolis.form = "not_a_form"
   rdb._errors.ConfigParseError: key 'form': 'not_a_form' not in allowed set
   {'sadourny', 'sadourny_hk', 'sadourny_energy'}

   >>> cfg.vcoord.vcoord_type = "ZSTAR_FULL"
   >>> cfg.vcoord.vcoord_type
   'zstar_full'

What it deliberately does *not* re-implement is the cross-knob semantics —
``validate_config`` and the schema's per-group ``cross_check`` callbacks stay
in Fortran, and surface as :class:`~rdb.ConfigValidationError` when the
composed namelist reaches the library. Two copies of a cross-knob rule is two
rules, and one of them will be wrong.

The knob metadata is readable from the class, without building anything:

.. code-block:: python

   >>> rdb.Config.ocean_epbl.mstar.default
   1.2
   >>> rdb.Config.ocean_epbl.mstar.doc
   'Constant-scheme mstar'

Absence is the default
----------------------

A fresh ``Config`` has every knob **unset**, and ``to_namelist()`` emits only
what was explicitly assigned:

.. code-block:: python

   >>> rdb.Config().to_namelist()
   ''
   >>> cfg = rdb.Config()
   >>> cfg.ocean_hvisc.nu_h = 5000.0
   >>> cfg.ocean_hvisc.smag_ah = True
   >>> print(cfg.to_namelist())
   &ocean_hvisc_nml
     nu_h = 5000.0
     smag_ah = .true.
   /

This is what makes a Python-built run bit-identical to the equivalent
namelist: an empty buffer leaves the Fortran's own pristine defaults
untouched, so nothing is silently re-stated at a Python-side idea of the
default. ``explicit_knobs()`` yields ``(group, key, value)`` for what was
actually set, and ``write_namelist(path)`` puts the same text on disk — which
is the honest way to record what a script ran.

Whatever the entry point, **that string is the only thing that crosses the
ABI**. A curated model, a hand-built ``Config`` and a namelist file all end
up in ``read_config_from_string``, through the identical parse, the identical
cross-checks and the identical ``validate_config``. There is no second,
mirrored configuration path to fall out of sync.

The escape hatch
----------------

A curated ``Model(...)`` defers creation until the first ``with`` block (or
an explicit ``model.create()``), precisely so the composed config can still
be edited:

.. code-block:: python

   model = rdb.Model(grid=..., closures=[...], ...)
   model.config.ocean_hvisc.bound_kh = True        # not a curated argument
   model.config.ocean_debug.chksum = True
   print(model.config.to_namelist())               # exactly what will be sent

   with model:                                     # create() happens here
       model.run(steps=100)

``defer=True`` requires the explicit ``create()`` even inside ``with``, for
scripts that want the creation point named.

Curated keywords and a raw namelist/``Config`` are **mutually exclusive
entry points** — one model, one recipe. Passing both raises ``TypeError``
rather than picking a winner. The escape hatch above is not a mix of the two;
it is editing the config the curated objects composed.

After creation ``model.config`` is a record, not a control panel: there is no
reconfigure-in-place entry point, so assigning to it later changes nothing in
the running solver. On a model built from a raw namelist *string*,
``model.config`` is an empty ``Config`` — the string is never parsed back
into typed knobs. ``Model.from_namelist`` is the path that does parse, below.

Dead knobs
----------

Some knobs are registered on the schema — accepted, type-checked,
round-tripped through ``to_namelist()`` — and read nowhere on the ocean path.
They were swept out once and flagged in the generator, so assigning one warns:

.. code-block:: python

   >>> cfg.vcoord.zstar_stretching = "log"
   RdbDeadKnobWarning: 'zstar_stretching' (z*-full: surface-concentration
   stretching) has no effect on the ocean path: ...

The assignment still *happens* — the knob is set, and it will appear in the
namelist text — because the warning is about effect, not legality. The best
known example is ``&output_nml output_to_file``, which several shipped
namelists still set and which does nothing; ``&ocean_diag_nml enabled`` is
the real switch. Promote the warning to an error while developing a script if
you want the stricter contract:

.. code-block:: python

   import warnings
   warnings.simplefilter("error", rdb.RdbDeadKnobWarning)

Starting from an existing namelist
----------------------------------

``Model.from_namelist`` is the migration path. It reads the file in Python,
resolves relative paths inside it against *the file's own directory* (so a
script's working directory stops mattering), and hands back a model whose
``.config`` is the same typed object a hand-written script would have built:

.. code-block:: python

   >>> model = rdb.Model.from_namelist("validation/case.nml", defer=True)
   >>> model.config
   Config(14 knob(s) set)
   >>> model.config.ocean_topo.max_depth
   200.0
   >>> model.config.ocean_hvisc.c_smag           # never set in the file
   MISSING

   >>> model.config.ocean_hvisc.nu_h = 50.0      # change one knob
   >>> model.create()
   >>> model.step(1)

A knob the file never mentioned reads back as ``MISSING``, exactly as on a
hand-built config — the object records what the file *said*, not what the
Fortran will default it to. That is what lets a user move one knob at a time
from a namelist into curated objects without ever holding two descriptions of
the same run.

Curated keyword arguments can be passed alongside, and a disagreement is
refused rather than resolved:

.. code-block:: python

   >>> rdb.Model.from_namelist("case.nml", closures=[ConstantViscosity(nu_h=1.0)])
   rdb._errors.ConfigConflictError: ... nu_h ...


Grids, geometry and topology
============================

The grid objects follow Oceananigans: size and extent on the grid, topology
on the grid, the bottom as an immersed boundary rather than a namelist
setting elsewhere. The spellings are Roundabout's, and where a concept has no
analogue it raises rather than approximating one.

.. code-block:: python

   from rdb.grids import (LatitudeLongitudeGrid, RectilinearGrid,
                          Topology, TripolarGrid)

   RectilinearGrid(size=(600, 500, 50), extent=(1200e3, 1000e3),
                   topology=(Topology.PERIODIC, Topology.BOUNDED))

   LatitudeLongitudeGrid(size=(360, 200, 50),
                         longitude=(0.0, 360.0), latitude=(-70.0, 70.0))

   TripolarGrid(size=(360, 240, 50), latitude_south=-80.0, phi_join=65.0)

``size`` is always ``(nx, ny, nz)`` and ``nz`` is a **layer count**, not a
list of interfaces. There is no ``z=`` argument anywhere, because Roundabout
has no explicit vertical-grid input in any form: the vertical is
``nz_layers`` plus a ``vcoord=`` family plus that family's scalar
parameters. Passing ``z=`` to a grid or to ``Model`` raises
:class:`~rdb.RdbUnsupportedError` saying so, rather than silently swallowing
it through ``**kwargs``.

``RectilinearGrid`` takes either ``extent=(Lx, Ly)`` or ``x=(lo, hi),
y=(lo, hi)``, and the object never computes ``dx``/``dy`` itself — both route
to ``len_lon``/``len_lat`` and the Fortran does the derivation, so there is
one source of truth for it. ``axis_units`` mirrors MOM6's ``AXIS_UNITS``:
under ``"degrees"`` the extent is read as degrees and ``dx`` becomes an arc
length, which is the mode the MOM6 double-gyre reference uses.

``SupergridGrid`` reads a MOM6 mosaic NetCDF and is the route to a stretched
horizontal mesh. Only ``file=`` works: the in-memory array form needs a
Fortran symbol that is still ``private``, and asking for it raises
:class:`~rdb.RdbUnsupportedError` saying which one.

``LatitudeLongitudeGrid`` derives topology when none is given: x is periodic
exactly when the longitude extent is 360 degrees, y is always bounded.
Periodic latitude raises — a pole is a physical boundary, not a wrap, and the
object points at ``TripolarGrid`` instead. ``TripolarGrid`` has no
``topology`` argument at all: it is periodic in x and north-folded by
construction, and it forces the matching ``&ocean_bc_nml`` edges itself.

Topology lives on the grid
--------------------------

``Topology.PERIODIC`` wraps both edges of a dimension; ``Topology.BOUNDED``
has two physical ends, and what happens at them is a ``boundaries=``
condition set separately. The default is ``(BOUNDED, BOUNDED)`` —
deliberately *not* Oceananigans' periodic-by-default, because a doubly
periodic ocean basin is rarely what was meant.

There is no ``Topology.FLAT``. Oceananigans' ``Flat`` is compile-time
dispatch to zero operators; a Roundabout run with ``ny=1`` still executes
live differencing kernels, so the concept does not carry over. The string
``"flat"`` in a topology tuple raises rather than being quietly accepted.

A grid's topology reaches ``&ocean_bc_nml`` whether or not a ``Boundaries``
object was given. This matters more than it sounds: without it, a grid built
periodic in x would derive its halo correctly, set its periodic flag — and
leave west and east at the default wall, turning a re-entrant channel into a
closed basin that runs, conserves mass, and is wrong.

The halo is derived, not guessed
--------------------------------

``halo=None`` (the default) asks the library what this grid and these
numerics need:

.. code-block:: python

   >>> from rdb.grids import required_halo
   >>> required_halo()
   2
   >>> required_halo(pv_advection="weno5")
   3
   >>> required_halo(periodic=True, north_fold=True)
   3

:func:`rdb.grids.required_halo` is a stateless call into Fortran — no handle,
no model — so it is the *same* max-over-every-rule the solver applies, never
a Python re-derivation that could drift from it. Give ``halo=`` explicitly
and it is checked: a value below the minimum raises, naming which option
binds. Roundabout cannot inflate a halo and rebuild, so this has to be right
before creation rather than fixed afterwards.

Bathymetry and the sign convention
----------------------------------

Array bathymetry is the ``GridFittedBottom`` analogue, and its ``convention``
argument is **required with no default**:

.. code-block:: python

   from rdb.grids import Bathymetry

   b = Bathymetry(depth,                            # (nx, ny), interior
                  convention="height_positive_up",  # GEBCO / ETOPO
                  land_threshold=2.0,
                  minimum_depth=10.0)

.. warning::

   ``"depth_positive_down"`` is Roundabout's own convention — a cell 4000 m
   deep is ``+4000``. ``"height_positive_up"`` is GEBCO's, ETOPO's and
   Oceananigans' — the same cell is ``-4000``. Getting this backwards makes
   every cell read as land, and the result is a clean, crash-free, entirely
   wrong quiescent run. There is no default because there is no safe one.

Which is why the array is checked, twice, rather than documented once. The
Python side normalises to positive-down and then looks at the wet fraction of
the normalised array; zero wet cells raises
:class:`~rdb.BathymetrySignError`, naming the median depth and the other
convention:

.. code-block:: python

   >>> Bathymetry([[-4000.0] * 6 for _ in range(8)],
   ...            convention="depth_positive_down")
   rdb._errors.BathymetrySignError: Bathymetry: normalised to 0 wet cells
   under convention='depth_positive_down' (median depth -4000.0 m,
   land_threshold=2.0 m). This is almost always the sign convention read
   backwards -- try convention='height_positive_up'.

Under one percent wet warns instead of raising. The Fortran stage call runs
the same check independently; neither side is trusted as the only guard.

Get it right and the object tells you what it holds:

.. code-block:: python

   >>> Bathymetry([[-4000.0] * 6 for _ in range(8)],
   ...            convention="height_positive_up")
   Bathymetry(8x6, depth 4000.0-4000.0 m, 100.0% wet, from height_positive_up)

``depth`` is interior-sized: ghost rows are filled by the library at stage
time, never by the caller. It accepts nested lists, or anything exposing
``__array_interface__`` — a numpy array crosses without this package ever
importing numpy.

The formula bathymetries are the other route, and take no array at all:
``FlatBottom(1000.0)``, ``Seamount(...)``, ``Spoon(...)``, ``Island()``,
``DoubleDrake()``, ``Neverworld2(...)``, ``BathymetryFile("topog.nc")``.
Length scales on ``Spoon`` and ``Seamount`` are in **metres**; the Fortran
converts to grid units at dispatch, so the same number means the same thing
on a Cartesian and a spherical grid.

Two-phase creation
------------------

Array geometry has to be injected before setup runs, which is what the
pending window is for. A curated ``Model(bathymetry=Bathymetry(...))`` — or
any periodic edge — uses it automatically. It is also directly available:

.. code-block:: python

   model = rdb.Model.pending(nml)          # config parsed and validated;
                                           # nothing allocated, no device map
   model.stage_bathymetry(depth, rdb.BATHY_DEPTH_POSITIVE_DOWN)
   model.stage_topology(periodic_x=True, periodic_y=False)
   model.finalize()                        # setup runs here

   print(model.bathymetry[0, 0])

``stage_bathymetry``, ``stage_metrics`` and ``stage_topology`` may be called
in any subset and any order between ``pending()`` and ``finalize()``. Calling
one outside that window raises :class:`~rdb.NotInitialisedError`; a staged
array whose extents do not match the configured interior raises
:class:`~rdb.BadShapeError`. If ``finalize()`` fails, the whole handle is
destroyed — there is no half-created model to inspect afterwards.


Fields
======

A :class:`~rdb.Field` is a live window onto one state array. It holds a host
pointer, the full Fortran extents, the step count the contents were last known
good at, a stagger location, units, and — if the array is writable — the
narrow Fortran setter that writes it.

What is on the model
--------------------

.. list-table::
   :header-rows: 1
   :widths: 16 20 12 52

   * - Attribute
     - Interior shape
     - Writable
     - What it is
   * - ``h``
     - ``(nx, ny, nz)``
     - yes
     - Layer thickness, m — the prognostic continuity variable.
   * - ``u``
     - ``(nx+1, ny, nz)``
     - yes
     - x-velocity on the west face, m/s.
   * - ``v``
     - ``(nx, ny+1, nz)``
     - yes
     - y-velocity on the south face, m/s.
   * - ``hu``, ``hv``
     - as ``u``, ``v``
     - no
     - Layer transports, m²/s.
   * - ``w``
     - ``(nx, ny, nz+1)``
     - no
     - Vertical velocity at layer interfaces, m/s.
   * - ``rho``
     - ``(nx, ny, nz)``
     - no
     - In-situ density, kg/m³.
   * - ``kv``, ``kt``, ``ks``
     - ``(nx, ny, nz+1)``
     - no
     - Vertical viscosity, heat diffusivity, salt (and passive-tracer)
       diffusivity, m²/s, at interfaces.
   * - ``bathymetry``
     - ``(nx, ny)``
     - yes
     - Depth, positive down, m.
   * - ``eta``
     - ``(nx, ny)``
     - no
     - The raw barotropic free-surface diagnostic.
   * - ``tau_x``, ``tau_y``
     - ``(nx+1, ny)``, ``(nx, ny+1)``
     - no
     - Wind stress, N/m². Written through ``model.set_wind(tau_x, tau_y)``.
   * - ``Q_heat``, ``Q_salt``
     - ``(nx, ny)``
     - yes
     - Net surface heat flux (W/m²) and salt flux.
   * - ``wet_mask``
     - ``(nx, ny)``
     - no
     - 1 wet, 0 land, at T points.
   * - ``temperature``, ``salinity``
     - ``(nx, ny, nz)``
     - yes
     - Concentrations, degC and PSU. See `Concentration, not the raw store`_.
   * - ``tracer(name)``
     - ``(nx, ny, nz)``
     - yes
     - Any registered tracer, **by name** — never by index.
   * - ``ssh``
     - ``(nx, ny)`` list
     - —
     - Sea-surface height as ``sum_k h - b``, the same definition the
       diagnostic manager uses. This is the one to prefer for "is the free
       surface where I expect".

A field carries its own description:

.. code-block:: python

   >>> model.h
   Field(h, (Center, Center, Center), shape=(8, 6, 3), unit='m', interior,
         generation=10)
   >>> model.u.location
   (<Loc.FACE: 'Face'>, <Loc.CENTER: 'Center'>, <Loc.CENTER: 'Center'>)

``Loc.FACE`` on the first axis is what makes ``u`` one wider than ``h`` in
that direction. Staggering is an attribute rather than something the caller
is expected to remember from the shape.

Interior by default
-------------------

``model.h.shape`` is the physical interior. The ghost-inclusive array is
reachable under a name:

.. code-block:: python

   >>> model.grid_info
   {'nx': 8, 'ny': 6, 'nz': 3, 'nghost': 3}
   >>> model.h.shape
   (8, 6, 3)
   >>> model.h.with_halo.shape
   (14, 12, 3)

There is no vertical halo, so only the horizontal axes grow. Making the
interior the default is what retires the ``[2:-2, 2:-2]`` idiom that
ghost-inclusive arrays force on every caller — and, more to the point, the
off-by-one that idiom produces when ``nghost`` turns out to be 3.

Concentration, not the raw store
--------------------------------

Tracers are stored as ``h*Tr``. ``model.temperature`` is degrees Celsius:

.. code-block:: python

   >>> T = model.temperature
   >>> T[0, 0, 0]                 # degC
   12.0
   >>> T.hTr[0, 0, 0]             # m*degC, the raw store
   799.9999999999999

:class:`~rdb.DerivedField` is computed, not mapped, so it is copy-only: it
has no contiguous buffer of its own, ``__array_interface__`` raises saying
so, and ``.hTr`` is the zero-copy view of the underlying store.

The division is guarded. Concentration is computed only where
``h > 1.5e-4 m`` (the solver's own dynamic-vanish threshold) and is **NaN**
elsewhere — not ``0.0``, because 0 degC and 0 PSU are both perfectly legal
ocean values and zero would read a vanished layer as ice-point freshwater
rather than as missing. On a vanishing-layer vertical coordinate, use
``math.isnan`` (or ``np.nanmean``) rather than assuming every cell has a
value.

Lazy synchronisation
--------------------

This is the design centrepiece, and it is worth understanding rather than
working around.

The GPU build is compiled with separate device memory. Between
``enter_data`` and ``exit_data``, the host copies of the state arrays are
**stale by default** — the solver runs entirely on the device and never
copies anything back unless something asks. An API that guaranteed
host-current data at every step boundary would force a device-to-host
transfer on every step, which is exactly the cost the design exists to avoid.

So Roundabout does not synchronise at step boundaries. It synchronises **when
you read**, and only then.

Every getter returns a pointer, the extents, and a **generation** — the
handle's outer step count at the moment the field was issued. ``step()``
merely invalidates: it clears the library's host-current flag and moves the
step count on. Nothing is copied. Then, on the next access to any field, the
field compares its own generation against the current step count; if they
differ it calls the library's refresh (an ``!$acc update self`` over the leaf
component arrays) and records the new generation. If they match, the read is
a pointer dereference and costs nothing.

.. code-block:: python

   with rdb.Model(nml) as model:
       h = model.h                    # issued at generation 0
       print(h._gen)                  # 0

       model.step(5)                  # no copy happens here
       print(h._gen)                  # still 0 -- nothing has been read

       value = h[0, 0, 0]             # THIS triggers the device -> host copy
       print(h._gen)                  # 5

Run that and the printed generations are ``0``, ``0``, ``5``. The field held
a stale pointer across five steps without anyone caring, and paid for exactly
one transfer, at the point where the answer was actually wanted.

The consequences worth stating plainly:

* **A field object stays valid across steps.** Taking ``h = model.h`` once
  outside a loop and reading it inside is correct, and is the intended usage.
  It is not a snapshot.
* **There is no ``sync()``.** Not on ``Model``, not on ``Field`` — the method
  does not exist, and the test suite asserts that it does not. If you find
  yourself wanting one, you want ``.copy()``.
* **Staleness never raises.** A field that has fallen behind re-syncs
  silently. The one case that *does* raise is a field belonging to a model
  that has been **closed** — there the pointer would genuinely dangle, so
  it raises :class:`~rdb.RdbClosedError` telling you to take a ``.copy()``
  before leaving the ``with`` block.
* **Refreshes are per-handle, not per-field.** The first read after a step
  brings the whole host side up to date, so reading ten fields in one
  callback costs one transfer, not ten.
* **A write poisons the generation rather than pulling eagerly.** After
  ``field[...] = x`` the device is canonical; the next read re-syncs.

``.copy()`` is the way to keep data across steps. It returns an owned
snapshot — nested Python lists — that does not track the solver:

.. code-block:: python

   before = model.temperature.copy()
   model.run(steps=100)
   after = model.temperature.copy()      # `before` is still the old field

Writes are calls, not view mutations
------------------------------------

Raw views are handed out **read-only**, on purpose. A writable mapped view
would let ``arr[...] = x`` mutate host memory that the next kernel launch
silently discards — a write that appears to work and does nothing. numpy is
told the buffer is read-only and refuses the assignment.

Writing goes through the field's own narrow Fortran setter instead:

.. code-block:: python

   nx, ny, nz = model.temperature.shape
   model.temperature[...] = [[[18.0] * nz for _ in range(ny)]
                             for _ in range(nx)]

   model.set_wind(tau_x, tau_y)          # the two components together

Each setter pushes only what changed, and does the bookkeeping that write
implies — flushing the tracer window before a tracer write, recomputing
density after one, re-deriving the barotropic reference depth after a
bathymetry write. A broad push would rewind live device state to a stale host
snapshot, which is why there is no "write everything" call.

Writing a read-only field raises with the list of the ones that are not:

.. code-block:: python

   >>> model.rho[...] = model.rho.copy()
   rdb._errors.RdbReadOnlyError: 'rho' is read-only. Writable fields:
   Q_heat, Q_salt, bathymetry, h, salinity, temperature, u, v. (ssh/eta are
   diagnostic -- sum_k h - b -- set `h` or `bathymetry` instead.)

Reading without numpy
---------------------

Three access patterns work with no dependencies:

.. code-block:: python

   model.h[...]          # nested lists, numpy-style indexing: a[i][j][k]
   model.h.copy()        # the same thing, named for what it is
   model.h[3, 2, 0]      # a full integer index -> one float

Anything else — ``model.h[:, :, -1]``, a mix of ints and slices, fewer
indices than dimensions — is handed to numpy over the same zero-copy
``__array_interface__``, when numpy is importable. When it is not, it raises
:class:`~rdb.RdbUnsupportedError` and says which of the three forms above to
use instead. The nested lists index the way numpy does, so
``full[i][j][k]`` and ``arr[i, j, k]`` agree.


Stepping
========

``step(n)`` advances ``n`` fixed-dt outer steps and does nothing else:

.. code-block:: python

   model.step()          # one
   model.step(100)       # a hundred, in one call into Fortran

``run()`` is the loop with a callback in it. Exactly one of ``until=``
(simulated seconds) or ``steps=`` (outer steps):

.. code-block:: python

   def report(m):
       print(f"t={m.time/3600:7.1f} h  KE={m.kinetic_energy:.3e} "
             f"mass={m.total_mass:.6e}")

   model.run(until=10 * 86400.0, callback=report, every=6 * 3600.0)
   model.run(steps=1000, callback=report, every=100)

``every`` is read in the same units as the bound: seconds under ``until=``,
outer steps under ``steps=``. Omit it and the callback fires every step.

``run()`` steps one outer step at a time, which is obviously correct for a
callback-driven loop and slightly wasteful without one — for a long unattended
stretch, ``step(n)`` is one call instead of *n*.

There is no synchronisation anywhere in that loop. The callback reads fields;
each read re-syncs itself if it needs to. A callback that reads nothing costs
nothing.

Scalars are properties on the model and always host-side:

.. code-block:: python

   model.time             # simulated seconds
   model.step_count       # outer steps taken
   model.total_mass
   model.kinetic_energy
   model.grid_info        # {'nx', 'ny', 'nz', 'nghost'}
   model.tracer_names     # ['salinity', 'temperature', ...]


Diagnostics
===========

The diagnostic machinery is the same registry and cadence dispatch a
namelist-driven run uses — the Python layer adds discovery and in-memory
reads on top of it.

What can be asked for
---------------------

:func:`rdb.available_diagnostics` is **static**: it needs no model, and can
be called before anything is created.

.. code-block:: python

   >>> rdb.available_diagnostics()
   ['SSH', 'temperature', 'salinity', 'age', 'u', 'v', 'KE', 'MLD_EPBL',
    'Kd_EPBL', 'Kd_KSHEAR', 'ice_conc', 'ice_thick', 'pseudo_salt',
    'pseudo_salt_diff', 'h_layer', 'rho_layer', 'vorticity_z', 'ke_total',
    'transport_x', 'transport_y', 'mld_density', 'ice_speed', 'ice_u',
    'ice_v']

That is the canonical catalog followed by the derived one, read out of the
library rather than listed here. It answers "what *can* be requested on some
configuration" — several names are conditionally registered (``temperature``
and ``salinity`` need thermodynamics enabled, the ``ice_*`` names need the
sea-ice model). For what *this* run actually registered:

.. code-block:: python

   >>> model.diagnostics.selected
   ['SSH', 'temperature', 'salinity', 'u', 'v', 'KE']
   >>> model.diagnostics.available == rdb.available_diagnostics()
   True

A name in ``available`` but not in ``selected`` was either gated off by
another namelist group or simply never asked for. Derived diagnostics are
never auto-registered.

Configuring the stream
----------------------

:class:`rdb.Diagnostics` writes ``&ocean_diag_nml``, plus the pieces of
``&output_nml`` and ``&logging_nml`` that go with it:

.. code-block:: python

   from rdb.diagnostics import Diagnostics, Duration

   rdb.Model(
       ...,
       diagnostics=Diagnostics(
           names=["temperature:6h:mean", "vorticity_z:z:1d", "KE:off"],
           every=3600.0,
           vgrid="z_fixed", levels=[0.0, 50.0, 100.0, 500.0, 2000.0],
           filename="tasman", output_dir="./out",
           precision="single", deflate=4,
           status_interval=86400.0),
       duration=Duration(t_end=30.0, time_unit="day"),
   )

``names=`` is the unified token-list selection: a name, optionally followed
by ``:``-separated attributes — ``off``, a time operation
(``instant``/``mean``/``max``/``min``/``integral``), an output coordinate
(``layer``/``z``/``zstar``/``sigma``/``density``), or a cadence like ``6h``
or ``1d``. Each token is validated against the same grammar the Fortran
parser enforces, so a typo is a ``ValueError`` at construction:

.. code-block:: python

   >>> Diagnostics(names=["temperature:not_a_real_attribute"])
   ValueError: unknown diagnostic attribute 'not_a_real_attribute' in spec
   token 'temperature:not_a_real_attribute'; valid: off |
   instant|mean|max|min|integral | layer|z|zstar|sigma|density|rho |
   <int>s/m/h/d

This matters because that particular Fortran branch still ``error stop``\ s,
which from Python means the interpreter dies. Catching it in the constructor
is the difference between a traceback and a lost session.

``Diagnostics(enabled=False)`` is the off switch, and the reason the object
exists: with no ``diagnostics=`` at all the Fortran defaults leave output on,
writing to a relative ``./output/`` path. ``output_dir``, ``status_interval``
and ``log_level`` are applied even when disabled, because they are console
and output plumbing rather than diagnostic selection.

.. note::

   ``Diagnostics(...).apply()`` also sets ``&output_nml output_to_file`` for
   parity with the hand-written escape hatch, which emits a
   ``RdbDeadKnobWarning`` — that knob is dead on the ocean path. The warning
   is accurate and harmless; ``&ocean_diag_nml enabled`` is the gate that
   actually does the work.

Reading a diagnostic out of memory
----------------------------------

``model.diagnostic(name)`` returns a :class:`~rdb.DiagnosticField` — a
read-only view onto the diagnostic manager's own output buffer. No NetCDF
file is written, and none is read back:

.. code-block:: python

   >>> ssh = model.diagnostic("SSH")
   >>> ssh.shape
   (8, 6, 1)
   >>> ssh[0, 0, 0]
   -2.842170943040401e-14

This is most of the point of driving the model in process. Otherwise you
write a file and read it back to obtain what memory already held.

The shape is ``(nx, ny, nz)`` in the usual interior convention: ``nz`` is the
layer count for a layered diagnostic on the default LAYER grid, ``1`` for a
2D one, or the remapped level count when that diagnostic has its own output
vertical grid.

A ``DiagnosticField`` does not track the outer step count — its generation is
the diagnostic's own per-fire counter, and it needs no refresh, because the
Fortran pulls the output buffer host-ward unconditionally at every cadence
fire, inside ``step()``. By the time ``step()`` returns, anything that fired
is already current. Reading a diagnostic that has never fired gives its
allocation-time value, which is zero.

Asking for a name that is not registered on *this* run raises
:class:`~rdb.TracerNotFoundError`, even if the name appears in
``available``; check ``diagnostics.selected``.


Time and restarts
=================

:class:`rdb.Duration` reaches ``&time_nml``:

.. code-block:: python

   from rdb.diagnostics import Duration
   Duration(t_end=30.0, time_unit="day")

``t_end`` does **not** bound ``run()`` — ``until=`` and ``steps=`` do that
directly, in raw simulated seconds. But it does gate creation, and it
defaults to **1.0 second**, which almost any real diagnostic cadence exceeds:
the Fortran cross-check refuses a ``dt_out`` larger than ``t_end`` once both
are scaled through ``time_unit``. Diagnostics on plus no ``duration=`` was a
reliable first-run stumble, which is what this object is for.

``time_unit`` scales ``t_end``, ``status_interval`` and ``dt_out`` together.
It does not scale ``dt_fixed``, which is always seconds — the same trap the
namelist has, carried over unchanged.

:class:`rdb.Restart` sets the restart half of ``&output_nml`` without
dropping to the escape hatch:

.. code-block:: python

   from rdb.diagnostics import Restart
   Restart(every=86400.0, file="tasman_restart.nc")

Restart is the thinnest part of this interface, and it is worth being plain
about why.

Writing a restart *from a live handle* has no entry point at all.
``model.write_restart()`` raises :class:`~rdb.RdbUnsupportedError` naming the
gap rather than quietly doing nothing. The same is true of
``model.write_diagnostics()``: the diagnostic cadence already runs inside
``step()``, but the C ABI exposes no way to force an out-of-cadence write, so
set ``every=`` small enough that the cadence covers what you need.

**Warm-starting from a restart file is not wired through this path.** The C
ABI's create routine deliberately does not pass ``restart_file`` through to
the setup chain — its own comment says so, in as many words, and calls it a
future phase's concern. ``Restart(file=...)`` writes the knob, and the knob
is carried in the namelist text faithfully, but a model created from Python
starts from its initial condition. Restarting a run is a job for the ``rdb``
executable today.


Error handling
==============

Every C ABI entry point returns a coarse stage code. The *specific* reason —
which key was unknown, which tracer was not found, which file could not be
read — lives in the Fortran's own error ring. The Python layer drains the
ring first and only then turns the code into an exception, so what you catch
carries the Fortran's own sentence:

.. code-block:: python

   >>> rdb.Model("&grid_nml nx = 8, ny = 6, totally_bogus_key = 1 /\n")
   rdb._errors.ConfigParseError: <config-string>:1: unknown key
   'totally_bogus_key' in group 'grid'

Not ``RuntimeError: error 1``. Every exception carries ``.code`` (the status
integer), ``.stage`` (a coarse category) and ``.trace`` (the whole ring dump,
deepest first); ``str(exc)`` is always the deepest, most specific message.

The hierarchy
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 8 62

   * - Exception
     - Code
     - Raised when
   * - ``RdbError``
     - —
     - Base of everything this package raises.
   * - ``RdbConfigError``
     - —
     - Base for configuration failures.
   * - ``ConfigParseError``
     - 1
     - Unknown group or key, wrong type, out of range, not in an enum —
       from the Fortran parser *or* from assignment to a
       :class:`~rdb.Config` knob.
   * - ``ConfigValidationError``
     - 2
     - A cross-knob semantic check failed (``validate_config``).
   * - ``RdbSetupError``
     - 3
     - The setup chain rejected the configuration — metrics, EOS, boundary
       validation.
   * - ``RdbInitialConditionError``
     - 4
     - Seeding the initial condition failed.
   * - ``RdbIOError``
     - 5
     - A setup-time file read or write failed. Also an ``OSError``, so an
       outer ``except OSError`` catches a missing forcing file.
   * - ``RdbHandleError``
     - —
     - Base for handle-lifecycle failures.
   * - ``InvalidHandleError``
     - 10
     - Null, stale or garbage handle.
   * - ``AlreadyExistsError``
     - 11
     - A second live ocean model. Raised Python-side where possible, so the
       message can name the model still open.
   * - ``NotInitialisedError``
     - 12
     - A step, query or getter on a handle that is closed, deferred, or was
       never created.
   * - ``TracerNotFoundError``
     - 13
     - No tracer — or no registered diagnostic — with that name.
   * - ``BadShapeError``
     - 14
     - A setter or staged array's extents do not match the interior.
   * - ``BathymetrySignError``
     - 15
     - A bathymetry array normalised to (near-)zero wet cells.
   * - ``NotPendingError``
     - 16
     - A ``stage_*`` or ``finalize`` call outside the pending window.
   * - ``RdbRestartError``
     - —
     - Base for restart failures.
   * - ``RestartSchemaMismatch``
     - 20
     - The restart file's schema does not match.
   * - ``RestartDecompMismatch``
     - 21
     - The restart was written by a different rank decomposition.
   * - ``RestartGridMismatch``
     - 22
     - Grid, vertical coordinate or tracer registry disagrees with the
       restart.
   * - ``RdbClosedError``
     - —
     - A field or model used after ``close()``. Python-side; never reaches C.
   * - ``RdbReadOnlyError``
     - —
     - A write to a field that has no setter.
   * - ``RdbUnsupportedError``
     - —
     - A real concept with no Roundabout analogue, or a deliberate refusal —
       ``z=``, ``Topology.FLAT``, partial slicing without numpy,
       ``write_restart()``.
   * - ``ConfigConflictError``
     - —
     - Two curated objects set the same knob to different values, or an
       object's prerequisite is missing. Python-side.
   * - ``LibraryNotFoundError``
     - —
     - No usable ``librdb_core.so``. A ``FileNotFoundError``, and listed
       alongside these because it is what an import-time failure looks like.

Catching a specific one
-----------------------

The point of the hierarchy is branching without matching on strings:

.. code-block:: python

   import rdb

   try:
       model = rdb.Model(candidate_config)
   except rdb.ConfigParseError as exc:
       print(f"typo in the config (code {exc.code}): {exc}")
       raise
   except rdb.ConfigValidationError as exc:
       print(f"the knobs disagree with each other: {exc}")
       raise
   except rdb.RdbIOError as exc:
       print(f"a file this run needs is missing: {exc}")
       raise

Or coarsely, by stage:

.. code-block:: python

   except rdb.RdbConfigError:     # parse + validate
   except rdb.RdbHandleError:     # anything about the handle's lifecycle
   except rdb.RdbError:           # everything

A parameter sweep is the case that earns it — a configuration that fails to
validate should skip to the next point, while a missing file should stop the
whole sweep:

.. code-block:: python

   results = {}
   for nu_h in (10.0, 100.0, 1000.0, 10000.0):
       cfg = base_config()
       cfg.ocean_hvisc.nu_h = nu_h
       try:
           with rdb.Model(cfg) as model:
               model.run(until=5 * 86400.0)
               results[nu_h] = model.kinetic_energy
       except rdb.ConfigValidationError as exc:
           print(f"nu_h={nu_h} rejected: {exc}")
           continue

Warnings
--------

Two warning types, neither of which stops anything:

``RdbDeadKnobWarning``
   A knob was assigned that the schema accepts and nothing reads. See `Dead
   knobs`_.

``NoBoundaryLayerWarning``
   ``closures=[...]`` left neither KPP nor EPBL active — typically
   ``KPP(enabled=False)`` with no ``EPBL(...)`` beside it. Interior mixing,
   background floors and convective adjustment still run, but nothing
   represents wind-driven surface boundary-layer entrainment, so a forced run
   gets an unrealistically shallow mixed layer. It warns rather than raises
   because it is exactly right for a closed-basin, no-wind dynamical-core
   test.


What it cannot do
=================

**One live model per process.** The halo module, the profiler, the
ice-ocean coupler and the diagnostic fills all hold module-scope state, so
multi-instance is not a missing feature but an impossible one. A second
``Model`` while one is open raises :class:`~rdb.AlreadyExistsError`. Ensemble
work means one process per member.

**No MPI.** Multi-rank through this path is explicitly out of scope: there is
no communicator argument, no rank query, and no decomposition control.
``comm_env_finalize`` calls ``MPI_Finalize`` and is irreversible, so it is
deliberately not exposed. A multi-rank run is a job for the ``rdb``
executable and a namelist.

**The GPU build works.** This was verified rather than assumed: the NVHPC
GPU ``.so`` does ``dlopen`` from stock CPython, ``acc_init`` creates a CUDA
context inside the loaded library and sees the devices, and a ``do
concurrent`` device kernel with Fortran I/O in it runs under a CPython main.
No Fortran main is needed. The lazy-sync design above is what makes that
worth doing — state stays device-resident, and a transfer happens only where
a script actually reads.

**No restart or out-of-cadence diagnostic write from a live handle.** The C
ABI has no entry point for either. Both raise
:class:`~rdb.RdbUnsupportedError` naming the gap. Restart and diagnostic
output are namelist-driven: set them before ``create()``.

**No grid coordinates.** ``field.nodes()`` raises — the C ABI exposes no
longitude/latitude accessor yet. Derive coordinates from ``grid_info`` and
the grid object's own parameters.

**No passive-tracer registration from Python.** The registry is locked by the
device map, so a tracer must be registered before ``enter_data``, which is
why tracers are a constructor argument rather than a method. The shipped
packages — ``IdealAge``, ``PseudoSalt`` — are reachable;
``PassiveTracer(...)`` for an arbitrary new name raises, since it needs a
Fortran-side entry point that does not exist yet.

**No vertical grid input, in any form.** Not a limitation of this layer —
Roundabout itself has none. See `Grids, geometry and topology`_.

**Configuration is one-way.** There is no reconfigure-in-place: editing
``model.config`` after ``create()`` changes a record, not a solver. Close the
model and build another.

**Errors in the Fortran are mostly, but not entirely, returnable.** The setup
chain was converted to return status codes, which is what makes the typed
hierarchy above possible. A few paths — notably some branches of the
diagnostic spec parser — still ``error stop``, and from inside CPython that
takes the interpreter down without a traceback. Where such a path was
reachable from a plausible mistake, the Python layer duplicates the *grammar*
check in front of it (see `Configuring the stream`_) so the mistake raises
instead. If a script dies with no Python traceback at all, that is the
signature.
