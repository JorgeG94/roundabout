.. _getting_started:

---------------
Getting started
---------------

.. contents::
   :local:

This page walks through one complete run: choosing a case, starting it,
reading what the model tells you while it works, and finding the output. It
assumes you have a working build — if not, see :ref:`installation`.


Running the model
=================

The solver takes exactly one argument, the namelist file:

.. code-block:: bash

   ./build/rdb validation_examples/ocean/double_gyre/double_gyre_mom6.nml

That is the entire command-line interface. There are no flags — no
``--help``, no ``--version``. Everything the model does is set in the
namelist. Give it no argument, or a path that does not exist, and it prints

.. code-block:: text

   FATAL: no input file provided.
   Usage: rdb <input_file.nml>

and exits with status 2.

Output lands in a directory the namelist chooses (``&output_nml
output_dir``, default ``./output``), **relative to the directory you run
from**. So it is usually easier to run from the case directory:

.. code-block:: bash

   cd validation_examples/ocean/double_gyre
   ../../../build/rdb double_gyre_mom6.nml

On a GPU build, pin the device explicitly:

.. code-block:: bash

   CUDA_VISIBLE_DEVICES=0 ../../../build/rdb double_gyre_mom6.nml

Multi-rank runs go through ``mpirun`` in the usual way, with the per-rank
device pinning described in :ref:`installation`:

.. code-block:: bash

   mpirun -np 4 ./build/rdb case.nml


Choosing a first case
=====================

``validation_examples/ocean/`` holds the canonical configurations. Three are
good starting points, for different reasons.

**``double_gyre/double_gyre_mom6.nml``** — the MOM6-reference double gyre.
44 × 40 cells, two layers, :math:`\Delta t = 1200\ \mathrm{s}`, ten
simulated days. Tiny, fast, and it exercises the full driver path: the
split-explicit dynamical core, the barotropic sub-cycle, bathymetry, wind
forcing and the diagnostics manager. This is the one to run first.

**``seamount/seamount.nml``** — an at-rest test. A 100 km closed basin with
a Gaussian seamount, uniform density, :math:`f = 0`, no wind and no forcing
of any kind. The correct answer is **no motion at all**, for thirty
simulated days. Any velocity it produces is a numerical artefact of the
vertical coordinate or the pressure gradient. About 30 s on a single V100.
Tests of this shape are unusually informative: there is no tuning that can
make a wrong answer look plausible.

**``eady/eady.nml``** — a 250 km reentrant channel with an Eady front,
50 × 50 × 10, run for sixty days. It has an analytical answer — linear
theory predicts a growth rate of :math:`1.98\times10^{-6}\ \mathrm{s^{-1}}`
for the 125 km mode in this channel, and the model reproduces it to about
3 %. Roughly two minutes on a single V100.

Every shipped namelist carries a header comment block describing the setup,
the expected outcome, and how to interpret a failure. Read it before you
run the case; several of them are the best documentation in the repository.


Reading the console
===================

The model announces itself with a banner, then a summary of the
configuration it is about to run:

.. code-block:: text

     Regime:   ocean (C-grid + continuity-PPM + split-RK2)
     Grid:     44 x 40
     Spacing:  dx = 55659.0 m, dy = 55660.0 m
     Layers:   2
     Duration: 864000.0 s
     Fixed dt: 1200.0 s
     Diag:     double_gyre every 86400.0 s

followed by a memory budget per state slot and the setup timings. Check the
grid spacing here. Several namelists *derive* ``dx`` and ``dy`` from a
domain extent in degrees rather than stating them, and a configuration
error at this point produces a basin of the wrong size that runs perfectly
happily.

The status line
---------------

At a cadence set by ``&logging_nml status_interval`` — or every 100 steps
if that is left at its default of zero — the model prints a ``[stats]``
line:

.. code-block:: text

   [stats] Day    0.000  step        0  En  0.000E+00  MaxCFL  0.00000  Mass  9.09864E+18  Salt   35.000  Temp    5.000

.. list-table::
   :header-rows: 1
   :widths: 14 14 72

   * - Field
     - Units
     - What it is
   * - ``Day``
     - days
     - Simulated time. Not wallclock.
   * - ``step``
     - —
     - Outer steps taken.
   * - ``En``
     - :math:`\mathrm{m^2\,s^{-2}}`
     - Kinetic energy **per unit mass** — so it reads as a velocity
       squared. :math:`10^{-4}` is a domain-rms of about 1 cm/s.
   * - ``MaxCFL``
     - —
     - The largest advective Courant number in the domain.
   * - ``Mass``
     - kg
     - Total mass, :math:`\sum h \, A \, \rho_0`.
   * - ``Salt``
     - PSU
     - Volume-mean salinity.
   * - ``Temp``
     - :math:`{}^\circ\mathrm{C}`
     - Volume-mean temperature.

``Salt`` and ``Temp`` appear only when thermodynamics is enabled
(``&ocean_thermo_nml enable_thermodynamics``); the double gyre above turns
it off, so its lines carry only the first five fields.

**How to read these.** ``En`` is the health indicator. What counts as
healthy depends on what the case is: in a forced run it should spin up and
then level off; in an unforced adiabatic run it may decay but must not
grow; in an at-rest case like the seamount it should stay at zero. ``En``
climbing steadily in a case with no energy source is the signal that
something is wrong.

``MaxCFL`` should sit comfortably below 1 — typically well below 0.1 in the
shipped cases, since the outer step is limited by advection and not by
gravity waves. If it approaches 0.9 the run aborts deliberately rather than
producing nonsense.

``Mass``, ``Salt`` and ``Temp`` are **conservation checks**, not physics.
Under closed boundaries and no surface flux they should not move at all.
Watching a digit change in ``Mass`` is how you catch a leak.

The budget block
----------------

Below each status line, one indented line per conserved quantity:

.. code-block:: text

       Mass :  9.098635585E+18  Error  0.000E+00
       Salt :  3.184522455E+20  Error  0.000E+00  out  0.000E+00  src  0.000E+00
       Heat :  4.549317792E+19  Error  0.000E+00  out  0.000E+00  src  0.000E+00
       En   :  0.000000000E+00  Growth         —

This is the model's **own closed budget**: it tracks everything that
entered through a source and everything that left through a boundary, and
``Error`` is the residual that is unaccounted for. That is a far stronger
statement than "the total did not change much" — a run with inflow and
outflow can have a wildly varying total and still close its budget exactly.
``Error`` should sit near round-off, of order :math:`10^{-11}` relative or
better for a closed domain.

The ``En`` line here reports **absolute** kinetic energy in joules, and its
growth factor relative to the start — unlike the ``En`` column on the
compact line above, which is per unit mass. The two are different numbers
with the same name.

If the sea-ice module is active, an ``Ice`` line joins the block with mean
concentration and thickness.

The progress table
------------------

Finally, a bare numeric row per status interval:

.. code-block:: text

        Step         t (s)         dt (s)     Wall (s)  Remaining (s)
     -------  ------------  ------------  -----------  -------------
          72      86400.00     1200.0000         0.63           5.65

This is wallclock, and it is where you find out whether the run will finish
today.


Finding the output
==================

Diagnostics are written **per rank**, with no gather:

.. code-block:: text

   <output_dir>/<filename>_rank_NNNNNN.nc

where ``<filename>`` is ``&ocean_diag_nml filename`` (default
``ocean_diag``). The double gyre above, with ``output_dir =
"./out_double_gyre_mom6"`` and ``filename = "double_gyre"``, produces

.. code-block:: text

   out_double_gyre_mom6/double_gyre_rank_000000.nc

Write cadence is ``&ocean_diag_nml dt_out``. A single-rank run produces a
single file you can open directly:

.. code-block:: bash

   ncdump -h out_double_gyre_mom6/double_gyre_rank_000000.nc

For a multi-rank run, ``tools/merge_output.py`` stitches the per-rank files
into one offline.

Rank 0 also writes two parameter dumps into the same directory —
``rdb_parameter_doc.all`` and ``rdb_parameter_doc.short`` — recording every
namelist value the run actually used. When you come back to output in six
months and cannot remember how it was configured, these are the answer.

Restarts, if ``&output_nml restart_interval`` is set, follow the same
pattern with a fixed prefix:

.. code-block:: text

   <output_dir>/restart_rank_NNNNNN.nc

They are overwritten in place rather than versioned by time. Point
``&output_nml restart_file`` at the file (or its directory) to warm-start
from one.

.. note::

   ``&output_nml output_to_file`` is a **dead knob** on this path. Several
   shipped namelists still set it; it parses and validates and then does
   nothing. To turn diagnostics off, use
   ``&ocean_diag_nml enabled = .false.``.


Where to go next
================

* :ref:`running_simulations` — how a namelist is organised, and walkthroughs
  of real cases.
* :ref:`theory` — the equations being solved.
* :ref:`validation` — how to check that a build is behaving.
