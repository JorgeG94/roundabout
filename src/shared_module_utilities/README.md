# `src/shared_module_utilities/` — shared include bodies

One place for small helper **bodies** that several modules need *inside* their
`do concurrent` kernels.

## What belongs here

A file here holds `pure` (and, where it is called from a device kernel,
`!$acc routine seq`) **procedure bodies** — no `module` statement, no `use`,
no `contains`. A consuming module `#include`s it in its own `contains`
section and so gets a module-LOCAL, `private` copy:

```fortran
module rdb_something
   use rdb_constants, only: wp, H_VANISHED, NZ_STACK_MAX
   implicit none
   private
   public :: something_public
contains

   subroutine something_public(...)
      ...
   end subroutine something_public

#include "rdb_vanished_layer.inc"

end module rdb_something
```

Admission test — all four:

1. **It is a rule, not a routine.** Several modules must agree on it, and two
   of them disagreeing is a bug rather than a style difference (the
   vanished-layer rule; a shared limiter; a shared predicate).
2. **It is small and leaf.** No state, no logging, no allocation, no I/O, no
   calls out except to intrinsics. `pure` always.
3. **It is called from inside a kernel**, where a module boundary is a
   performance question rather than a free abstraction.
4. **One concern per file**, named for the concern.

If it does not meet all four it is an ordinary module procedure. Put it in a
module.

## Why an include and not a module

CLAUDE.md used to prescribe hand-duplicating a helper body into each calling
kernel module as a `*_impl` copy, because NVHPC's device codegen was not
inlining `pure !$acc routine seq` helpers across a module boundary. Hand
duplication gives you N copies to keep in step, and they drift.

An include gives every consumer its own module-local copy of ONE source text.
The compiler sees a same-module leaf procedure — the case every toolchain
inlines most readily — and the tree has a single definition to edit.

**Measured** (nvfortran 26.5, `-O2 -stdpar=gpu -gpu=cc70,mem:separate`, one
V100; reproducer kept out of the tree, see the commit that added this
directory): a per-cell accessor written three ways — raw inline expression, a
module-local `#include`d helper, and a cross-module `!$acc routine seq` helper
— produces device PTX with **zero `call` instructions in all three**; all three
kernels are within two lines of each other in PTX length and carry exactly one
inlined `div.rn.f64`. On a 600×600×50 bandwidth-bound sweep the three are
indistinguishable (0.0179 / 0.0173 / 0.0172 s for 30 launches).

So on **this** compiler version the cross-module form inlines too, and the
include's justification is **one source of truth**, not a speed-up. The
include also never depends on that staying true, which hand duplication was
invented to work around.

One observable difference worth knowing: passing an array element to a
function (either helper form) makes nvfortran report the enclosing array as
`implicit copy(a(:,:,:))` rather than the sliced `copyin(a(:nx,:ny,:nz))` it
reports for the raw expression. Under `mem:separate` with the array already
device-resident this is "if not already present" and costs nothing — but a
kernel whose arrays are *not* mapped would pay for it, which is one more
reason the data-residency rules in CLAUDE.md are not optional.

## Rules

- **`#include` (cpp), not Fortran `include`.** Sources are `.F90`, so the
  preprocessor runs on every toolchain the project builds with, and the
  include path is the ordinary `-I` one. The directory is on the include path
  of every target (library, app, tests, benchmarks) — see
  `RDB_SHARED_INCLUDE_DIR` in the top-level `CMakeLists.txt`.
- **Distinctive, prefixed names.** Every procedure in a file here is prefixed
  (`rdb_vl_*` for the vanished-layer rule) and every local variable inside it
  carries the same suffix (`k_vl`, `orphan_vl`), because the body lands in a
  module that has its own names. Consumers keep the included names `private`
  — which the repo's `private`-by-default module header already does.
- **The consuming module supplies the `use`.** An include cannot carry one;
  the file's header comment says what it needs (for `rdb_vanished_layer.inc`:
  `wp`, `H_VANISHED`, and `NZ_STACK_MAX` for the column merge).
- **The file is the authority.** Do not copy a body out of here into a module
  "for speed" — measure first, and if you must, say so where you do it.

## Tooling

`.inc` does not match the `\.(f|F|f90|F90)$` filter the `fprettify` and
`fortitude` pre-commit hooks use, so a file here is not auto-formatted or
linted. Match the house style by hand; it is a small file by construction.
The bodies ARE compiled (and so type-checked, and FORD-documented) through
every module that includes them.

## Contents

| File | Concern |
|---|---|
| `rdb_vanished_layer.inc` | The vanished-layer rule — `rdb_vl_is_live`, `rdb_vl_conc`, `rdb_vl_merge_content`. Contract: `src/core/ocean/README.md`, "The vanished-layer content rule". |
