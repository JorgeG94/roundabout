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

So on **this** compiler version the cross-module form inlines too — for a
SMALL helper — and the include's justification there is **one source of
truth**, not a speed-up. The include also never depends on that staying true,
which hand duplication was invented to work around.

**It does not hold for a large helper.** `rdb_eos`'s `roquet_spv_ts_coeffs`
(the (T, S) half of the Roquet SpV polynomial: two sqrt, a degree-7 PT→CT
polynomial, ~50 terms) called from the FV-MOM6 PGF kernels stayed a real
`call rdb_eos_roquet_spv_ts_coeffs_` in the PTX (nvfortran 26.5, cc70), its
four results through a 48-byte device stack frame, 130+ registers in the face
kernels. Included module-locally (`rdb_roquet_spv.inc`) it inlines: no stack,
86–88 registers, and the global 1° `ocean_pgf` under Roquet fell 3.64 → 2.59 s
with bit-identical `[stats]`. Check the PTX (`-gpu=keep`, grep `call`) rather
than assuming either way.

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
  carries the same suffix (`k_vl`, `q_sum_vl`), because the body lands in a
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
| `rdb_rel_vort_corner.inc` | The C-grid corner relative vorticity ζ (circulation form, with the C1 free-/no-slip factor) — `rdb_rvc_zeta_corner`. One consumer today, the `vorticity_z` derived diag; `rdb_coriolis_adv`'s four `q_corner` passes still carry inline twins of it, marked in the source, pending migration. |
| `rdb_roquet_spv.inc` | The VALUE of the Roquet et al. (2015) SpV polynomial (ΔS = 24), split into its (T, S) part (`rdb_roq_ts_coeffs`) and its pressure part (`rdb_roq_spv_p`) — `rdb_eos`'s `roquet_spv_ts_coeffs` / `roquet_spv_value` and the FV-MOM6 PGF's per-EOS density integrals. Also carries `rdb_roq_ct_of`, the shared PT->CT helper `rdb_roquet_rho.inc` calls too (one of three copies of the TEOS-10 polynomial tree-wide — see the file header). The consumer imports the coefficient table from `rdb_eos` (see the file header). |
| `rdb_roquet_rho.inc` | The VALUE of the Roquet et al. (2015) in-situ-DENSITY polynomial (ΔS = 32, independent fit from the SpV one above — same paper, not bit-identical), split into its (T, S) part (`rdb_roq_rho_ts_coeffs`, calling `rdb_roquet_spv.inc`'s `rdb_roq_ct_of`) and its pressure part (`rdb_roq_rho_p`, no reciprocal) — `rdb_eos`'s `roquet_rho_ts_coeffs` / `roquet_rho_value` and the FV-MOM6 PGF's `roquet_rho_pcm_dpa_intz` / `_face`. Must be `#include`d AFTER `rdb_roquet_spv.inc` in the same module `contains` (needs `rdb_roq_ct_of` in scope). |
| `rdb_vanished_layer.inc` | The vanished-layer predicate and donor read (tracers are concentrations; a filler keeps its last value) — `rdb_vl_is_live`, `rdb_vl_column_conc`. Contract: the header of the file itself. |
| `rdb_bed_live_face.inc` | The "is this layer live at this face" rule (#178) used to walk a static bed index `k_bot_u/v` up to the first row still live after a free-surface drawdown — `rdb_blf_is_live`. Consumers: `rdb_ocean_vdiff`'s bed-BC row placement, `rdb_ocean_bottom_drag`'s explicit/implicit-fold drag and HBBL band walk — they must agree on `kb_live` for the same face or the drag law samples one row while the solver glues the BC to another. |
| `rdb_tracer_conc.inc` | Tracer source deposits into a concentration: `rdb_conc_add_source` (`c + src/(h + h_neglect)`), `rdb_conc_accepted_source` (the part a layer takes) and `rdb_conc_add_source_mass` (source plus a thickness change). A layer at or below `H_VANISHED` refuses a deposit and keeps `c`; callers book `rdb_conc_accepted_source`. Consumers: the surface flux, SSS restore, geothermal, cavity-flux and KPP non-local deposits. |
| `rdb_safe_math_consts.inc` + `rdb_safe_math_{core,pow,trig,atan,hyp,erfc}.inc` | Position-independent transcendental functions (`rdb_sm_*`), written only in IEEE `+ - * /`, `sqrt` and compile-time tables, so a value does not depend on vector position, tile, math library or optimisation level; special values propagate (NaN in, NaN out; `exp(+Inf) = +Inf`, `log(0) = -Inf`; `pow` with a NaN argument is NaN, deliberately unlike C99); every result with |r| < `tiny` is flushed to exactly zero (the sign of a flushed zero is not guaranteed: ifx returns +0 where gfortran keeps the sign) (`rdb_sm_flush`, so flushing and non-flushing compilers agree); subnormal arguments: `log` finite unless the compiler reads them as 0 (then `-Inf`), `cbrt` reads them as 0. **Include set:** the consts file ONCE in the module's SPECIFICATION part; then, in `contains`, only the body files needed, in this order: `_core` (two-sum/two-prod, `exp`, `log`; required by `_pow`, `_hyp`, `_erfc`), `_pow` (`pow`, `cbrt`), `_trig` (`sin`, `cos`, `tan`; needs `_core` for the flush helper), `_atan` (`atan`, `atan2`; needs `_core`), `_hyp` (`tanh`, `cosh`, `acosh`), `_erfc` (`erfc`). Each file's header lists what it provides and requires. The consumer supplies `use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite`, `iso_fortran_env`'s `int64` and `wp` (real64: a compile-time check stops the build otherwise). Infinity and NaN are `SM_INF` / `SM_NAN` constants, not `ieee_value()`: on the device nvfortran materialises that call's class argument as a derived-type global and two classes in one module fail in NVVM. Public face: module `rdb_safe_math` (`src/framework/`), elemental `safe_*` wrappers for host / setup code and the golden-bit test. **Contraction policy:** the library adds no per-file flag and follows the build's global FMA setting. Golden bits (`tests/test_safe_math_golden.inc`) are a contract for the contraction-OFF CPU builds only (gfortran `-ffp-contract=off`, nvfortran `-Mnofma`, ifx `-no-fma`); a GPU or other contracting build differs by 1-2 ulp on some arguments and its test asserts position independence on the device, the ulp bounds and the special values instead. Multi-term sums are written one step per statement on purpose: nvfortran `-O3 -fast` re-associated the final four-term sum of `log` when written as one expression (`-Mnofprelaxed` restored it; the split form needs no flag); an explicit `-Mfprelaxed=order` still breaks the golden table. Fitted tables: `tools/gen_safe_math_coeffs.py`; golden table: `tools/gen_safe_math_golden.F90`. **Verified:** the golden table reproduces bit for bit on gfortran 15, nvfortran 26.5 CPU (both via CMake) and ifx 2022.1 (standalone driver only); flang-new, LFortran (a `transfer` in a parameter) and an ifx CMake build are untested. **Admission-rule exception:** this set is large (seven files, ~900 lines) and its constants file is included in the SPECIFICATION part, not `contains`; it is admitted because the rule it carries (one bitwise-reproducible implementation of each function) is the whole point and the pieces are small leaf `pure` routines called inside kernels. Adopters so far (module-local includes): `rdb_barotropic_coupling` (consts + core + pow, for `bt_rem`), `rdb_atm_state` (core + trig, per-cell solar zenith), `rdb_ocean_bulk_forcing` (core + trig, grid-angle rotation), `rdb_ocean_bulk_flux` (core + atan, saturation humidity and stability functions), `rdb_ocean_surface_flux` (core + pow, shortwave transmission and Manizza opacity), `rdb_ocean_epbl` (core + pow + erfc), `rdb_ocean_vmix` (core + pow + trig + atan + hyp), `rdb_ocean_cavity_melt` (core + pow), `rdb_ocean_tidal_mixing` and `rdb_ocean_bbl_mixing` (core), `rdb_ocean_meke` (core + pow), `rdb_ocean_mle` (core + pow). Host/setup users of the `rdb_safe_math` facade: `rdb_barotropic_substep` (tidal OBC target), `rdb_ocean_tides`, `rdb_ocean_tide_astro`, `rdb_atm_state` (Spencer solar position), `rdb_ocean_surface_stress` (analytic wind profiles), `rdb_ocean_lateral_mix` (latitude floor), `rdb_ocean_ideal_age` (growth `exp`). |
