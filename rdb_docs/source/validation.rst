.. _validation:

----------
Validation
----------

.. contents::
   :local:

Roundabout is checked at three levels: unit and analytical tests through
CTest, a golden-comparison regression suite, and a two-tier **stability
suite** that asserts on the model's own console output over physically
meaningful integrations.


Unit and analytical tests
=========================

Built with the project (``RDB_ENABLE_TESTING``, default ``ON``) and run
through CTest:

.. code-block:: bash

   cd build
   ctest --output-on-failure

Every Roundabout test is named with an ``rdb_`` prefix, so

.. code-block:: bash

   ctest -R rdb --output-on-failure

runs the project's own suite — around 177 tests, roughly 35 seconds —
without the dependency self-tests that ``pic`` and ``test-drive`` also
register. Tests carry mandatory regime labels, so ``ctest -L ocean`` and
``ctest -L core`` select subsets.

.. warning::

   **Never pass ``-j N`` to ctest on a GPU build.** Every worker shares the
   one GPU; parallel execution produces spurious failures and hangs.
   Building with ``-j`` is fine.

A report of "N/N passing" should say which scope it ran — the ``rdb``
subset and the full suite are different numbers.

The tests that earn their keep
------------------------------

The unit suite covers the usual ground, but the **analytical** tests are
the ones that have repeatedly caught real bugs, because they compare
against a closed-form answer rather than against yesterday's output:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Test
     - What it pins down
   * - ``test_ocean_analytical``
     - Vertical diffusion against the complementary error function;
       gravity-wave phase speed; geostrophic-adjustment balance.
   * - ``test_ocean_validation``
     - Cooling-column heat budget, Ekman-transport sign, conservation
       under no forcing.
   * - ``test_ocean_conservation_salt_heat``
     - Closed salt and heat budgets. Doubles as the canonical template for
       GPU device mapping.
   * - ``test_ocean_land_mask``
     - An interior land block stays **exactly** at rest.
   * - ``test_ocean_remap`` / ``test_ocean_remap_e2e``
     - ALE remap conservation — drift :math:`\le 10^{-15}` over ten
       remaps.
   * - ``test_ocean_restart``
     - Bit-exact restart round-trip, including the wet/dry hysteresis
       registry.
   * - ``test_ocean_tripolar`` / ``test_ocean_fold``
     - Vector sign flips and seam antisymmetry across the north fold.

The pattern worth copying when you add physics: find a configuration whose
answer you know independently, and assert on that. An analytical anchor
cannot be satisfied by a plausible-looking wrong answer, and it does not
have to be re-baselined when something legitimate changes.


The stability suite
===================

Unit tests are short by nature, and a golden baseline records whatever the
model did on the day it was captured. Neither catches a defect that takes
hours of simulated time to appear, and a golden captured from a
misconfigured run will keep passing forever.

The stability suite closes both gaps. It runs shipped namelists for long
enough that their physics actually manifests, and it asserts on the model's
**own console time series** — no stored baseline at all. Five families of
assertion:

#. **Finite** — nothing is NaN or infinite.
#. **Conservation** — the model's own closed-budget ``Error`` term stays at
   round-off.
#. **Energy** — kinetic energy behaves as the case's regime requires.
#. **CFL** — the Courant number stays bounded.
#. **The case's own claim** — the testable statements in the namelist
   header.

Regimes
-------

Each case declares what kind of thing it is, which determines what the
energy assertion means:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Regime
     - Energy must
   * - ``rest``
     - Stay at zero. Any energy at all is spurious.
   * - ``adiabatic``
     - Not grow. Decay is permitted.
   * - ``baroclinic``
     - Grow — and at tier 1, at the analytically predicted rate.
   * - ``forced``
     - Spin up, stay under a physical ceiling, and at tier 1 saturate.

The bars are chosen from physics, not from what a run happens to do. The
at-rest seamount bar, for instance, is a mean kinetic energy corresponding
to about 1 cm/s of spurious current — the Beckmann & Haidvogel (1993)
figure for sigma-coordinate pressure-gradient error — and cases that should
be exactly zero are gated near machine precision instead.

The two tiers
-------------

**Tier 1** runs each case at full scale on a GPU, for a per-case duration
long enough that its physics appears — tens of seconds to minutes per case.
This is the tier that can measure a baroclinic growth rate or watch a
multi-day saturation. It needs a GPU and is not part of the automated
gate.

**Tier 2** runs downscaled twins on a CPU with gfortran, sized for an
ordinary CI runner. Each twin is seconds. It keeps the finite,
conservation, energy and CFL assertions, and drops only the ones that
genuinely need a long integration — a growth rate cannot be measured in
forty simulated hours — saying so rather than weakening them.

The important part is **how** a twin is built. Not by shrinking the grid
and hoping: by preserving the **dimensionless numbers** that make the case
what it is. A twin must still resolve its boundary current, still retain an
eddy field, still span enough deformation radii across the box. Those rules
are checked *before* any twin runs, and a violation aborts the whole run
without executing the model:

.. code-block:: text

   DOWNSCALE RULE VIOLATION -- this twin is not the same test as its
   tier-1 parent:

A grid-shrunk twin that no longer resolves what its parent resolved is not
a cheaper version of the test — it is a different test that happens to
share a name.

The scheme axis
---------------

``&ocean_bt_nml split_scheme`` selects the outer time integration
(:ref:`discretisation`), and results are not bit-identical between the two
values. Both are therefore kept under test.

Every eligible case is run a second time with the non-default scheme
forced, under the name ``<case>__ssp_rk2``, against the **same** physics
assertions. The reasoning is that a change of integrator must not change
what a case is *allowed* to do — so anywhere the two disagree, the
disagreement is itself the finding.

Cases that pin a scheme in their own namelist are excluded from the axis,
and that exclusion is **derived from the file** rather than maintained as a
hand-written list. Tier 2 runs the full axis; tier 1 runs a curated subset
chosen to cover each vertical-coordinate family, both pressure-gradient-
sensitive cases, the surface-flux and open-boundary cases, and the
land-masked cases.

Running it
----------

Standard-library Python only — the project does not install packages.

.. code-block:: bash

   # The CI gate: downscaled CPU twins
   python3 tests/regression/stability.py --tier 2 --build-dir build_gcc

   # Full scale on a GPU farm
   python3 tests/regression/stability.py --tier 1 --backend gpu \
           --build-dir build --gpus 0,1,2,3

   # A couple of cases, verbosely
   python3 tests/regression/stability.py --tier 2 --cases eady,acc_channel -v

Useful flags: ``--cases`` to filter, ``--jobs`` for CPU parallelism,
``--gpus`` for one worker per device, ``--out FILE.json`` for machine-
readable results, ``--keep`` to retain the NetCDF (it is deleted by default,
pass or fail), and ``--self-test`` to verify the assertions themselves
against recorded signatures without running the model at all.

Reading the output
------------------

One line per case:

.. code-block:: text

   [ PASS] eady                              123.4s  validation_examples/ocean/eady/eady.nml

with ``[*FAIL]``, ``[XFAIL]`` and ``[XPASS]`` as the other states. Only
failing assertions print unless you pass ``-v``.

Two subtleties are worth knowing.

A ``NAN-CATCH:`` line means the solver's in-flight NaN repair fired. Treat
that as a defect **even when every other assertion passes** — the model
papered over something.

An expected-failure marker is scoped to specific tiers and specific
assertions. Any failing assertion *not* covered by the marker is a real
FAIL, not an XFAIL. Without that scoping, one documented issue would excuse
a case from every gate it has, which is how a genuine regression hides
behind an old marker.

Exit codes: ``0`` no failures, ``1`` at least one failure, ``2`` binary not
found, ``3`` a downscale rule violation (nothing ran).

Artefacts land in ``tmp_local_artifacts/stability/tier<N>/<case>/``: the
patched namelist actually used — the committed file is never modified — and
``run.log``, the full output of the run.


The golden-comparison suite
===========================

A separate, older corpus of 28 bathymetry-free cases, each a handful of
outer steps, comparing field summaries against committed baselines. It
covers the run-clean gate (exit zero, no NaN) and drift against the
goldens.

.. code-block:: bash

   python3 tests/regression/run_all.py                    # CPU compare + coverage
   python3 tests/regression/run_all.py --backend both     # CPU + GPU + coverage
   python3 tests/regression/run_all.py --update-golden    # regenerate baselines

``run_all.py`` orchestrates the golden comparison and a gcov coverage
measurement; each stage runs in its own subshell sourcing the right
toolchain environment, so ``--backend both`` works from one invocation. It
does **not** run the stability suite — invoke that separately.

Coverage is treated as a measurement, not a gate: it fails the suite only
on an infrastructure error, never on "coverage is low".


Before you commit
=================

The project's standing order, in order:

#. ``pre-commit run --all`` — whitespace, trailing newlines, end-of-file.
#. ``fortitude check`` — Fortran static analysis.
#. ``ctest`` — the full suite, **without** ``-j``, green on the GPU build.
#. Re-read the hand-maintained synthesis documents. They rot silently, and
   a stale cell in a table is a bug. Where a document and the code
   disagree, **the code is the authority and the document is the bug**.

When adding a capability, the shape is: a namelist knob defaulting **off**
so existing configurations stay bit-identical; the kernel behind that knob;
an analytical test covering the new path; then the three checks above. One
capability per commit, one per pull request.
