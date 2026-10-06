# Roundabout Fortran Style Guide

The conventions for all Fortran in this repository: `src/`, `app/`, `tests/`
and `benchmarks/`. Dependencies (`pic`, `pic-mpi`, `test-drive`) follow their
own.

This guide is the *how*. [`CLAUDE.md`](CLAUDE.md) is the companion: its
**Gotchas** section records the traps that have already cost real debugging
time (GPU data motion, compiler bugs, the vanished-layer rule), and it wins
where the two disagree. Much of what follows is enforced by the pre-commit
hooks and CI (see [Enforcement](#enforcement)); the rest is checked in review,
including by the `pr-review` skill (`.claude/skills/pr-review/SKILL.md`).

## Naming Conventions

### Files

- Source files: `rdb_<name>.F90`. Every file is `.F90` (preprocessed), so
  `#ifdef`, `#include` and the CMake-defined macros work everywhere.
- Tests: `tests/test_<name>.F90`; MPI tests in `tests/mpi/test_<name>_mpi.F90`.
- Shared helper bodies: `src/shared_module_utilities/rdb_<concern>.inc` (see
  [Sharing a helper body](#sharing-a-helper-body-across-kernel-modules-include-not-a-hand-copied-_impl)).
- Free-form source only; no `.f` / `.F` files.

### Modules

- The module name matches the file name: `rdb_<name>`.
- One module per file (submodules excepted).

```fortran
! Good
module rdb_ocean_sponge

! Bad
module OceanSponge
module ocean_sponge        ! missing rdb_ prefix
```

### Derived Types

- `snake_case` with a `_t` suffix.

```fortran
! Good
type :: multilayer_state_t
type :: ocean_metrics_t
type :: sponge_band_t

! Bad
type :: MultilayerState
type :: ocean_metrics      ! missing _t suffix
```

### Variables and Procedures

- `snake_case`, descriptive. Single letters only for loop indices
  (`i`, `j`, `k`, with `k = 1` the bed) or a very local, obvious quantity.
- Grid sizes follow the existing vocabulary: `nx`, `ny`, `nz`, `ng` (ghost
  width), `nx_total = nx + 2*ng`.

```fortran
! Good
integer :: n_tracers
real(wp) :: layer_thickness
subroutine compute_bottom_drag(ms, metrics, cfg)

! Bad
integer :: nTracers        ! camelCase
real(wp) :: lt             ! too short to read
```

### Constants

- `UPPER_SNAKE_CASE` for every `parameter`.

```fortran
! Good
integer, parameter :: MAX_REMAP_ITERATIONS = 100
real(wp), parameter :: SPONGE_MIN_TIMESCALE = 3600.0_wp   ! s

! Bad
integer, parameter :: maxIter = 100
```

## Required Practices

### Use Statements

- **Always use `only`** (fortitude enforces it), and keep the list to what the
  file actually uses.
- Intrinsic modules need the `intrinsic` attribute (fortitude rule C122).

```fortran
! Good
use rdb_constants, only: wp, GRAVITY, H_VANISHED
use rdb_multilayer_state, only: multilayer_state_t
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

! Bad
use rdb_constants                 ! fails the linter
use iso_fortran_env, only: int64  ! missing `intrinsic`
```

### Implicit None

- `implicit none` in every module and program.

### Intent Declarations

- **Every dummy argument declares its intent.** No exceptions.

```fortran
! Good
subroutine apply_bottom_drag(ms, metrics, dt)
   type(multilayer_state_t), intent(inout) :: ms
   type(ocean_metrics_t), intent(in) :: metrics
   real(wp), intent(in) :: dt

! Bad
subroutine apply_bottom_drag(ms, metrics, dt)
   type(multilayer_state_t) :: ms      ! missing intent
```

### Declaration Order

- Declare integer dimensions **before** the explicit-shape arrays that use
  them. ifx warns otherwise (#8586); gfortran and nvfortran accept it silently.
  Enforced by the `decl-order` hook.

```fortran
! Good
integer, intent(in) :: nx, ny, nz
real(wp), intent(inout) :: h_layer(nx, ny, nz)

! Bad - the bound references a dummy declared later
real(wp), intent(inout) :: h_layer(nx, ny, nz)
integer, intent(in) :: nx, ny, nz
```

### Private by Default

- Modules are `private` by default, with an explicit `public ::` list.
- Make something public only when another module needs it.

```fortran
module rdb_ocean_example
   use rdb_constants, only: wp
   implicit none
   private

   public :: example_state_t
   public :: example_apply

   type :: example_state_t
      ! ...
   end type example_state_t

   type :: example_workspace_t   ! stays private
      ! ...
   end type example_workspace_t
end module rdb_ocean_example
```

### Pure by Default

- **Make a new procedure `pure`** unless it genuinely needs a side effect:
  I/O, logging, MPI, mutating module state, or an `allocate` +
  `!$acc enter data`. The bar is "can this be pure?", not "should it be?".
- Every procedure called from inside a `do concurrent` must be pure.
- Use `elemental` for scalar operations that should also work on arrays.
- A procedure that is impure only because it lazily allocates scratch wants
  that scratch lifted onto a state slot, allocated at init and freed at
  destroy (see [Persistent kernel workspaces](#what-remains-genuinely-painful)).

```fortran
! Good - pure, scalar, callable from a do concurrent body
pure function face_thickness(h_left, h_right) result(h_face)
   real(wp), intent(in) :: h_left, h_right
   real(wp) :: h_face
   h_face = 0.5_wp*(h_left + h_right)
end function face_thickness
```

### Limit Procedure Arguments

- A **public** procedure takes **six arguments or fewer**. Group related ones
  into a derived type: roundabout passes `cfg`, `grid`, `metrics` and the
  state objects (`ms`, `bt_work`, ...) rather than their components.
- **Exceptions:** private `*_impl` kernels that receive explicit-shape arrays
  (the shim pattern below needs flat arrays, not structs), and simple
  utilities.

```fortran
! Bad - a public entry point taking the pieces
subroutine sponge_apply(h_layer, hTr, u, v, nx, ny, nz, ng, tau, width, dt)

! Good - the public entry point takes the objects ...
subroutine sponge_apply(sponge, ms, grid, dt)
   type(ocean_sponge_t), intent(in) :: sponge
   type(multilayer_state_t), intent(inout) :: ms
   type(hgrid_t), intent(in) :: grid
   real(wp), intent(in) :: dt
   ! ... and hands flat explicit-shape arrays to a private *_impl kernel
```

## Kinds, Constants and Units

### Kind Parameters

- Working precision is **`wp` from `rdb_constants`** (currently `real64`), with
  `_wp` literals. `RDB_ENABLE_DOUBLE` selects it; never hard-code a kind.
- Integer kinds that must be exact (`int64` for the EFP bins, file offsets)
  come from `iso_fortran_env`.

```fortran
! Good
use rdb_constants, only: wp
real(wp) :: dt
dt = 1200.0_wp

! Bad
real(8) :: dt              ! non-portable
real*8 :: dt               ! obsolete syntax
double precision :: dt     ! ignores RDB_ENABLE_DOUBLE
dt = 1200.0d0              ! literal kind
```

### Physical and Numerical Constants

- **Physical constants live in `rdb_constants` and nowhere else**: `GRAVITY`,
  `RHO_WATER`, and the rest. Never re-declare or truncate one locally — a
  hard-coded `9.81` against `GRAVITY = 9.80665` once produced a spurious
  barotropic force wherever the free surface sloped.
- Thin-layer constants have documented roles: `H_VANISHED` (1.5e-4 m) for
  dynamic vanishing (skip or merge, don't clamp) and `H_DIV_EPS` (1e-20) as
  pure divide-by-zero armour. Pick the right one.

### No Magic Numbers

- Name any non-obvious literal (thresholds, tolerances, iteration caps) as a
  `parameter`, with its unit in a comment. `0.5_wp` in a two-point average,
  or `2.0_wp` in a formula, needs no name.

```fortran
! Bad
if (h_layer(i, j, k) < 0.01_wp) cycle
tol = 1.0e-10_wp

! Good
real(wp), parameter :: MIN_MIXING_THICKNESS = 0.01_wp      ! m
real(wp), parameter :: REMAP_CONSERVATION_TOL = 1.0e-10_wp ! relative
```

### Units

- **SI throughout**: m, s, kg, Pa, kg m⁻³, °C (potential or conservative
  temperature, as the EOS defines), psu (or g/kg for absolute salinity).
- Write the unit in the `!!` docstring of every physical field, and in a
  comment wherever a quantity is not in SI (e.g. a lat/lon in degrees).
- **Formula-bathymetry length scales are in grid units**: metres on Cartesian
  grids, degrees on spherical and curvilinear ones. Convert a metres knob
  with `topo_length_to_grid_units`.

```fortran
real(wp) :: eta            !! sea-surface height (m)
real(wp) :: tau_x          !! zonal wind stress (Pa)
real(wp) :: lat_deg        ! degrees north, not SI
```

### The Vertical Layer Convention

- The layer stack is **bottom-up**: `k = 1` is the bed and `k = nz` is the
  surface. Surface forcing lands at `(:, :, nz)`, bed forcing at `(:, :, 1)`.
  Every layered kernel follows this; don't flip it.

## Error Handling

- **Fail loud, never silently fall back.** An unsupported combination is
  refused at configure time with a message that names the knob, not
  "corrected" behind the user's back. A `select case` over a scheme ends in
  `case default` + an error, not a fallback.
- **Status codes, not exceptions.** Procedures on the solver-creation path
  (config parse and validate, the `configure_ocean_*` chain, IC seeding) take
  an optional `integer, intent(out) :: ierr` and return one of the named
  `OCEAN_STATUS_*` codes from `rdb_ocean_status`. A library cannot abort its
  host process, so when `ierr` is present they return the code; when it is
  absent they `error stop`.
- Log the reason with `global_logger%error` before returning, and check the
  code at the call site.

```fortran
subroutine configure_example(cfg, ocean_state, ierr)
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   type(config_t), intent(in) :: cfg
   type(ocean_state_t), intent(inout) :: ocean_state
   integer, intent(out), optional :: ierr

   if (present(ierr)) ierr = OCEAN_STATUS_OK
   if (cfg%ocean%example%width < 1) then
      call global_logger%error("&ocean_example_nml width must be >= 1")
      if (present(ierr)) then
         ierr = OCEAN_STATUS_ERR_SETUP
         return
      end if
      error stop "configure_example: invalid &ocean_example_nml width"
   end if
end subroutine configure_example
```

## Forbidden Practices

### No GOTO Statements

- `goto` is forbidden under any circumstance; fortitude catches it and
  `continue` labels. Use structured control flow: `do`, `if`,
  `select case`, `exit`, `cycle`, named `block ... end block`.

```fortran
! Bad - bail with goto
   if (bad) goto 999
   ! ...
999 continue
   return

! Good - early return
   if (bad) then
      call global_logger%error("...")
      return
   end if
```

**In `test-drive` tests** (the most common place this pattern appeared), wrap
the checks in a named `block` and `exit` it, so the cleanup after the block
still runs:

```fortran
! Bad - labelled cleanup target
   call check(error, c1, "m1")
   if (allocated(error)) goto 99
   call check(error, c2, "m2")
   if (allocated(error)) goto 99
99 call cor%destroy()
   call ms%destroy()

! Good - named block, structured exit, single cleanup
   checks: block
      call check(error, c1, "m1")
      if (allocated(error)) exit checks
      call check(error, c2, "m2")
      if (allocated(error)) exit checks
   end block checks
   call cor%destroy()
   call ms%destroy()
```

For **phased cleanup** (each resource released in reverse order of
acquisition), nest one block per acquire/release pair. Each release sits right
after its `end block`:

```fortran
   call nc_open_read(fname, ncid)
   ncid_scope: block
      call check(error, ...)
      if (allocated(error)) exit ncid_scope
      allocate (read_buf(...))
      buf_scope: block
         call check(error, ...)
         if (allocated(error)) exit buf_scope
      end block buf_scope
      deallocate (read_buf)
   end block ncid_scope
   call nc_close(ncid)
```

### No Arithmetic IF, COMMON, EQUIVALENCE or `external`

- Use `if`/`else if` or `select case`; modules and derived types; proper type
  conversions; and `use` of a module procedure, respectively.

### No Implicit Module State

- No module variable that silently carries state between calls.
- **The one sanctioned exception is a persistent kernel workspace:** a
  module-level `allocatable` scratch array, allocated once by a
  `*_workspace_ensure`, mapped to the device, and released by a `*_cleanup`
  in the exit-data path. That is the required pattern for per-step GPU
  kernels (see below), not a violation of this rule. Document it as such.

### No Naked Output

- **Never** `print *`, `write(*,*)` or `write(6,*)`. Use `pic_logger`'s
  `global_logger` with a real level. Watch for debug writes left behind.

```fortran
use pic_logger, only: global_logger
use pic_strings, only: to_string

call global_logger%info("ALE remap: " // to_string(n_remapped) // " columns")
call global_logger%warning("sponge width exceeds the tile; clipped")
call global_logger%error("&ocean_bc_nml: CHAPMAN edges are single-rank")
```

### No Emojis in Fortran

- None in `.F90` / `.inc` files: comments, strings and docstrings included.
  Python tooling may use them.

### No Assumed-Size Arrays

- Use explicit-shape (kernels) or assumed-shape (host code) dummies, never
  `arr(*)`. Inside `do concurrent` kernels, assumed-shape is out too (see
  below).

### No Scalar Assignment to a Whole Deferred-Length Character Array

- `character(len=:), allocatable :: lines(:)` followed by `lines = ""`
  re-allocates `lines` with `len=0` (F2018 10.2.1.3). ifx does this correctly;
  gfortran and nvfortran don't, which is how a reader silently discarded
  every in-memory namelist on ifx. Blank-fill through a section: `lines(:) = ""`.

## Recommended Practices

### Use pic and pic-mpi

| Library | Use for | Instead of |
|---|---|---|
| `pic_logger` | logging (`global_logger`) | `print *` / `write(*,*)` (**forbidden**) |
| `pic_timer` | timing | manual `cpu_time` |
| `pic_strings` | string helpers (`to_string`, ...) | hand-rolled conversions |
| `pic_mpi_lib` | **all** MPI (`comm_t`, `allreduce`, `isend`, ...) | `use mpi_f08` / `use mpi` (**forbidden**, `no-mpi-in-rdb` hook) |
| `testdrive` | tests (`new_unittest`, `check`) | ad-hoc test drivers |

- `pic_mpi_lib` abstracts `mpi_f08` versus legacy `mpi` (`USE_LEGACY_MPI`),
  which some systems only provide one of. Comm code lives in `src/comm/`; a
  kernel never talks to MPI directly.
- There is no BLAS or LAPACK in this code. Don't introduce one for a small
  dense operation a loop does fine.

### Array Operations

- Prefer whole-array syntax on the host when it is clear. On the device path,
  write the loop as a `do concurrent` (a whole-array assignment on a
  device-resident array may run on the host).

### Block Constructs

- Use `block` to scope temporaries to the branch or loop that needs them.

```fortran
do it = 1, ms%n_tracers
   block
      real(wp) :: total_before
      total_before = column_content(ms, it)
      call remap_tracer(ms, it)
      call check_conservation(ms, it, total_before)
   end block
end do
```

### Associate — Host Code Only

- `associate` is fine for readability in host code.
- **Never wrap a `do concurrent` kernel in `associate` over a derived-type
  component.** Under `-qopenmp`, ifx 2025/2026 reads such a name as zero
  inside the loop body (it once deleted the whole Coriolis term), and NVHPC
  has a sibling failure with `!$omp target`. Spell the component out, or pass
  it to an explicit-shape `*_impl` dummy.

```fortran
! Bad - ifx reads f_corner as 0 inside the loop
associate (f => metrics%f_corner)
   do concurrent (j = 2:ny, i = 2:nx)
      q(i, j) = f(i, j)*...
   end do
end associate

! Good
do concurrent (j = 2:ny, i = 2:nx)
   q(i, j) = metrics%f_corner(i, j)*...
end do
```

### Documentation

- `!!` is the FORD documentation marker, placed **below** the declaration it
  documents. FORD renders the docs from these on every PR to `main`; don't
  write a parallel markdown description of a module or kernel (that copy rots).
- Document every public type, procedure, type component and namelist field:
  what it is and its unit, not what the code obviously does.
- Use `!` for ordinary comments inside a routine body.

```fortran
type :: ocean_sponge_t
   !! North/south/east/west relaxation bands toward a target T/S/u/v state.

   real(wp), allocatable :: rate(:, :)
      !! Relaxation rate per cell (s^-1); zero outside the band.
   integer :: width = 0
      !! Band width in cells, counted inward from the edge.
end type ocean_sponge_t
```

- When a change adds, removes or renames a closure, scheme, knob or extension
  seam, update the hand-maintained synthesis docs in the same PR:
  `docs/CLOSURE_MATRIX.md`, `docs/CAPABILITIES_AND_LIMITATIONS.md`,
  `src/core/ocean/README.md` and the relevant `docs/howto/` guide.

### Citing Algorithms

- Cite the **paper**, not another model's source file: `!! Wright (1997)
  nonlinear EOS`, not a comment naming the file it was ported from. Source
  paths rot, read as "we copied this", and bury the provenance a reader
  needs.
- Reading a non-GPL reference implementation (MOM6, ROMS, ...) to understand
  an algorithm is fine; implement it from the maths and cite the paper it
  cites. **GPL code is off-limits** (e.g. GOTM): don't read, port or
  line-by-line reference it. A closure whose only reference implementation is
  GPL is written clean-room from the published equations. The interface (what
  a scheme consumes and returns) is fine to mirror.

### Memory Management

- Every type with allocatable components has a `destroy`, and a type that maps
  arrays to the device has matching `enter_data` / `exit_data` (see below).
- Every array a type owns is counted in its `bytes()` (the
  `check-bytes-accounting` hook checks this).

### Prefer Allocatable Over Pointer

- `allocatable` unless you need pointer semantics (aliasing, linked
  structures). Allocatables deallocate automatically and optimise better.

### Allocatable Character Strings

- `character(len=:), allocatable` for strings of unknown length, not a fixed
  `character(len=256)`.

### Avoid Deep Nesting

- Three or four levels at most. Use early `return`, `cycle` and `exit`
  (outside `do concurrent`), or extract a helper.

```fortran
! Bad
do k = 1, nz
   if (active) then
      if (h_layer(i, j, k) > H_VANISHED) then
         ! work buried here
      end if
   end if
end do

! Good
if (.not. active) return
do k = 1, nz
   if (h_layer(i, j, k) <= H_VANISHED) cycle
   ! work at a readable depth
end do
```

### Clamps Must Not Launder NaN

- Under nvfortran's relaxed floating point (`-fast`, no `-Kieee`), an
  `if/else` clamp is lowered to a NaN-blind min/max, so a NaN comes out as the
  bound. Any clamp that can see non-finite data needs an `ieee_is_finite`
  guard. Comparisons with NaN are false, so `if (abs(u) > thresh)` silently
  skips a NaN too.
- Fixed-point (EFP) reductions propagate non-finite input through their
  `poison` counter; never convert a NaN to an integer bin.

## Do Concurrent + OpenACC (Roundabout's GPU Path)

Roundabout's GPU offloading is **`do concurrent` + OpenACC**, compiled with
NVHPC `-stdpar=gpu -acc=gpu` and `-gpu=...,mem:separate`. This gives us:

- data-parallel loops in portable Fortran (`do concurrent`): the same source
  runs serially, on multicore (`-stdpar=multicore`) and on the GPU;
- OpenACC directives (`!$acc enter data`, `!$acc update`,
  `!$acc parallel loop reduction(...)`) for what `do concurrent` can't
  express: explicit device data management and reductions. They are inert
  comments on compilers built without `-acc`, so they are written directly,
  with no `#ifdef`.

NVHPC has matured a lot, and much older folklore no longer applies. This
section records what is true of the current (26.x) compiler.

### Device data residency (`mem:separate`)

There is **no managed or unified memory** (`-gpu=managed` is banned). A
kernel gets no implicit host↔device copies, so every array it touches must
already be on the device, or it silently reads stale memory. State arrays stay
device-resident; only diagnostic output and user callbacks move data back.

- Map state once, at the start of a long-lived scope (`enter_data`), and
  keep it resident.
- `create`-mapped arrays carry no host values: push host-set inputs with
  `!$acc update device(...)` after the map, and read results back with
  `!$acc update self(...)`.
- **Setup-time host editors run before the map.** Masking a metrics slot after
  it has been mapped leaves the device copy unmasked while every host check
  passes. Pass the mask to the constructor
  (`make_cartesian_metrics(metrics, grid, wet_mask=...)`) so the order is
  unskippable.
- **Never `update` or `copyin` a whole derived type with allocatable
  components**: it overwrites the host descriptors with device addresses.
  Update the component arrays.
- Tests must map everything they hand to a kernel, state **and** scratch
  companions. A test that only runs on the CPU build passes while being
  GPU-broken; verify device data motion on the GPU build. `CLAUDE.md`'s
  *Writing GPU tests* gotcha has the full checklist.

### Required pattern for every `*_enter_data` / `*_exit_data`

Every derived type whose components are referenced from `do concurrent` (or
`!$acc parallel loop`) bodies attaches the parent struct with
`!$acc enter data copyin(this)` **before** any component, and reverses the
order on exit (components first, parent last).

```fortran
subroutine foo_enter_data(this)
   class(foo_t), intent(inout) :: this
   !$acc enter data copyin(this)              ! PARENT FIRST
   !$acc enter data copyin(this%arr1, this%arr2)
   call this%sub_struct%enter_data()          ! recursive: same rule
end subroutine foo_enter_data

subroutine foo_exit_data(this)
   class(foo_t), intent(inout) :: this
   call this%sub_struct%exit_data()
   !$acc exit data delete(this%arr1, this%arr2)
   !$acc exit data delete(this)               ! PARENT LAST
end subroutine foo_exit_data
```

**Why:** when a `do concurrent` body references `foo%arr1(i, j)`, NVHPC needs
the parent's descriptor on the device to resolve the component. Without
`copyin(this)` it falls back to page-faulting the descriptor on every launch;
`nsys` shows memcpys eating 90%+ of kernel time. Adding the missing lines to
about 13 ocean `enter_data` routines turned a 60 s smoke run into 5 s.
`ocean_bc_state_t`'s `enter_data` in
`src/core/ocean/boundary/rdb_ocean_boundary_types.F90` is a current example.

Omitting a component, even one no kernel reads, leaves an invalid pointer in
the parent's device descriptor; the first kernel to touch anything on it fails
cryptically.

### What works well

- **`do concurrent` over large `k, j, i` nests with `local(...)`.** Compiles
  cleanly and coalesces well.
- **`!$acc update self/device`** for in-place host↔device sync, with no
  exit/enter round trip.
- **Reductions with `!$acc parallel loop collapse(n) reduction(...)`**, e.g.
  the console sums in `src/core/ocean/diag/rdb_ocean_console_stats.F90`.

```fortran
! Good - independent iterations, scratch via local()
do concurrent (j = ng + 1:ny - ng, i = ng + 1:nx - ng) local(flux_e, flux_w)
   flux_e = u(i + 1, j)*0.5_wp*(h(i, j) + h(i + 1, j))
   flux_w = u(i, j)*0.5_wp*(h(i - 1, j) + h(i, j))
   dhdt(i, j) = -(flux_e - flux_w)*idx(i, j)
end do
```

### What remains genuinely painful

- **Indirect access through a derived-type-array component inside a kernel.**
  Writing `ms%tracers(it)%hTr(i, j, k)` directly inside a `do concurrent` or
  `!$acc parallel loop` hits `CUDA_ERROR_ILLEGAL_ADDRESS`: two levels of
  descriptor indirection aren't reachable on the device. **Fix: outer shim +
  inner impl.** The outer routine dereferences the registry on the host and
  passes flat arrays to a private `*_impl` with explicit-shape dummies; the
  call collapses the indirection. Every tracer-registry loop does this, e.g.
  the sponge's `relax_map_tracer_impl` and `relax_band_x_impl`
  (`src/core/ocean/boundary/rdb_ocean_sponge.F90`) and
  `ocean_halo_exchange_ml_state`.

  ```fortran
  ! Bad - the kernel chases the registry indirection on the device
  do concurrent (k = 1:nz, j = 1:ny, i = 1:nx)
     ms%tracers(it)%hTr(i, j, k) = ms%tracers(it)%hTr(i, j, k)*(1.0_wp - rate(i, j))
  end do

  ! Good - the shim hands a flat array to the kernel
  call relax_tracer_impl(ms%tracers(it)%hTr, rate, nx, ny, nz)

  pure subroutine relax_tracer_impl(hTr, rate, nx, ny, nz)
     integer, intent(in) :: nx, ny, nz
     real(wp), intent(inout) :: hTr(nx, ny, nz)
     real(wp), intent(in) :: rate(nx, ny)
     ! ...
  end subroutine relax_tracer_impl
  ```

- **Explicit-shape dummies in kernels.** With an assumed-shape dummy
  (`arr(:, :, :)`) referenced inside a `do concurrent`, NVHPC walks the array
  descriptor on every launch: the vmix incident produced 1.4 M per-launch
  memcpys in 52 s. Pass the integer dimensions and declare `arr(nx, ny, nz)`,
  as `kappa_shear_merge_into_kv_kt` and `epbl_merge_into_kv_kt` do.
  Cadence-bounded paths (diagnostic fills once per output frame, init-only
  routines, registry shims receiving arrays of varying size) may waive with
  `! assumed-shape-ok: <reason>` on or just above the declaration.
  Enforced diff-aware by the `dc-assumed-shape` hook.

- **Persistent kernel workspaces.** Never `allocate` scratch inside a
  per-step kernel on `-stdpar=gpu`: it stalls the device on entry and exit.
  Use a module-level `allocatable`, lazily allocated by `*_workspace_ensure`
  and released by `*_cleanup` in the exit-data path.

- **A host-gated call that passes a state array still costs, even when not
  taken.** Handing e.g. `ms%mass_flux_x_layer` as an `intent(inout)` actual
  to another subroutine makes nvfortran treat the array as escaping and
  pessimises every `do concurrent` in the calling routine: +4.8 % total solver
  time for a call behind a knob that was off. Write the guarded pass inline
  as a `do concurrent` in the same routine; a same-module helper does not help.

- **Default-initialised derived types in device code.** An
  `!$acc routine seq` procedure that builds a type with default initialisers
  (function result, `intent(out)` dummy or local) can segfault nvfortran's
  `fort2` with no diagnostic, depending on file and module name length. Pass
  scalars instead; `intent(in)` bundles are fine.

- **Device stack size.** The GPU thread stack is small, so per-thread
  automatic arrays are bounded by `NZ_STACK_MAX` (`rdb_constants`, set by the
  `RDB_NZ_STACK_MAX` CMake option). Size column scratch with it, and let
  `validate_config` refuse a larger `nz`. The LFortran build pins a
  module-local copy under `#ifdef LFORTRAN_PASSING`; a new importer needs one
  too.

- **Cross-iteration read/write of the same array (Jacobi vs Gauss–Seidel).**
  A `do concurrent` that reads a neighbour's element and writes its own element
  of the **same** array gives different answers per backend: the GPU reads the
  pre-update image (Jacobi), the serial CPU path sees earlier writes
  (Gauss–Seidel). The result is a small, accumulating CPU↔GPU divergence that
  only a bit-compare catches. **Fix:** read from a snapshot and write to the
  live array, as the barotropic substep does with its `bt_eta_new` scratch.

- **`transfer()` inside `do concurrent`** is miscompiled by nvfortran 26.5
  (most elements silently wrong); keep such loops on `!$acc parallel loop`.
  Locals that shadow an intrinsic name inside a `do concurrent` hit a codegen
  bug on NVHPC and ifx. Both are hook-enforced (`dc-transfer`,
  `dc-intrinsic-shadow`).

- **Multi-GPU device selection.** On a multi-GPU node,
  `omp_set_default_device(gpu_rank)` sets only the OpenMP device; the OpenACC
  device stays 0 unless `!$acc set device_num(gpu_rank)` is called too.
  `rdb_comm_env::comm_env_setup_roles` does both. With CUDA-aware MPI, also pin
  `CUDA_VISIBLE_DEVICES` per rank before `MPI_Init` (see `CLAUDE.md`, *MPI*).

- **OpenMP target offload is not our path.** We compile `-stdpar=gpu
  -acc=gpu`, not `-mp=gpu`. Don't mix `!$omp target` directives with
  `do concurrent` in new code.

### Other NVHPC / gfortran codegen quirks

- **Put the contiguous (leading) array index LAST** in the `do concurrent`
  header: `do concurrent (k = ..., j = ..., i = ...)` for `arr(i, j, k)`. On
  NVHPC `-stdpar=gpu` this is order-immune (stdpar auto-collapses the nest and
  maps the contiguous index to the warp lane; `-Minfo` prints
  `auto-collapsed-innermost`; all permutations measured identical). A compiler
  that lowers the nest to a literal collapsed loop is not: a non-contiguous
  innermost index strides global memory, historically 16–23× slower on a
  `do concurrent`-to-`!$omp` rewrite. So `(k, j, i)` is free on stdpar and may
  not be free elsewhere.
- **`-gpu=managed` is banned.** A host-allocated scratch written in a
  `do concurrent` and read back with `sum()` silently returns 0 under managed
  memory. Use `!$acc parallel loop reduction(...)` for reductions over fresh
  scratch.
- **Don't split a `present(optional)` test into two `do concurrent` loops.**
  NVHPC compiled the no-optional branch of a split remap about 25× slower.
  Gate the optional inside one loop.
- **Keep a loop-invariant `select case` inside one `do concurrent`.** One loop
  per case doubles the launches (about 1.8× slower); a uniform branch is
  essentially free.
- **Don't hand loop-swap or hoist invariants on bandwidth-bound kernels.**
  NVHPC already does it, and the extra scratch arrays add register spills (a
  MUSCL rewrite came out about 20 % slower). Micro-optimise only when a profile
  shows the kernel is compute-bound.
- **CUDA-graph capture.** A bare `do concurrent` can't be captured (stdpar
  inserts a sync per loop); wrapping the body in `!$acc kernels async(q)`
  makes it capturable, about 5× on launch-bound sequences.
- **gfortran 15.1 `local(x)` reassignment.** gfortran corrupts a `local(x)`
  scalar reassigned across several `if/else` branches (NVHPC is fine). Use one
  local per write site (`h_face_e/w/n/s`, not one reused `h_face`).
- **Workspaces mapped `create` / `delete` can't be unit-tested by host
  pre-fill**: the host buffer isn't the device buffer. Assert by
  run-twice-and-compare instead.
- **An `ieee_get_flag` assertion only sees this thread's FPU.** A device
  kernel raises no host flag, so such a check is vacuous on the GPU build;
  skip it there (`KERNEL_IN_HOST_FPENV`), never let it "pass".
- **No GPU kernel behind overridable (polymorphic) dispatch.** A
  `do concurrent`, or a routine that launches one or maps state, invoked as a
  type-bound procedure on a `class(...)` actual goes through a runtime vtable
  that NVHPC can't inline into device code. At best that is a missed inline;
  at worst the `class(...)` box is mapped to the device and the run crashes
  (seen on AMD). **Rule:** a TBP that contains or launches a kernel is
  `non_overridable`, or it is a thin polymorphic wrapper that immediately
  delegates (via `select type`) to a non-polymorphic `*_impl` that does the
  device work, the ocean-state pattern in `src/core/ocean/README.md`.

### Sharing a helper body across kernel modules: `#include`, not a hand-copied `*_impl`

When several kernel modules must agree on a small `pure` (+ `!$acc routine
seq`) leaf helper called from inside a `do concurrent`, put the **body** in
`src/shared_module_utilities/<concern>.inc` and `#include` it in each
consumer's `contains` section:

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

Every consumer gets a module-local, `private` copy, which every toolchain
inlines readily, and the tree has one source text to edit. The include path
(`RDB_SHARED_INCLUDE_DIR`) is on every target. The admission test (a rule, not
a routine; small and leaf; called inside a kernel; one concern per file), the
naming rules (prefixed procedure names and suffixed locals, since the body
lands in a module with its own names) and the tooling caveats are in
`src/shared_module_utilities/README.md`.

This replaces the older advice to hand-copy a helper body into each kernel
module as a `*_impl`. Measured on nvfortran 26.5 (`-O2 -stdpar=gpu
-gpu=cc70,mem:separate`, V100), a per-cell accessor written as a raw inline
expression, as an included module-local helper and as a **cross-module**
`!$acc routine seq` helper produced device PTX with zero `call` instructions in
all three, and timed the same. So on this compiler the cross-module form
inlines too (older releases failed with `NVFORTRAN-S-1061` or wrong results),
and the include's value is a single source of truth rather than speed, with
the bonus that it does not depend on that staying true. Existing
hand-duplicated copies are historical; don't add more, and don't keep a second
inline copy of a formula an include already holds.

Passing an array **element** to a helper makes nvfortran report the enclosing
array as a whole-array `implicit copy`. That is free when the array is already
device-resident, and one more reason the residency rules are not optional.

### Restrictions that apply to any `do concurrent`

- No `exit`, `cycle`, `return` or `goto` in the body (a language rule).
- No side effects: I/O, logging and impure calls are undefined behaviour.
- Reductions use `!$acc parallel loop reduction(...)`, not the loop alone.

### Rule of thumb

- **Default to `do concurrent` + `local(...)`** for data-parallel kernels.
- **Drop to `!$acc parallel loop`** for reductions, explicit collapse, or the
  constructs `do concurrent` mis-compiles (above).
- **Never rely on `do concurrent` for a reduction or an ordered operation.**
  F2023's `do concurrent (...) reduce(+:s)` is legal, but we don't use it:
  `!$acc parallel loop reduction(...)` is an inert comment to compilers
  without `-acc`, so the same source is a correct serial reduction on
  gfortran and ifx and a GPU reduction on NVHPC, while `reduce` support is new
  and uneven across the three. A plain `do concurrent` with `s = s + ...` and
  no clause is a data race.

## Compiler Portability

CI builds and tests every PR on gcc 15, Intel 2025 and ifx + Intel MPI,
nvfortran (+ HPC-X), gfortran + OpenMPI, **flang-new and LFortran with NetCDF
off**. Code has to compile on all of them.

- **NetCDF-only code** sits under `#ifndef RDB_NO_NETCDF`: the source, the
  test cases, and their registration in `tests/tester.F90` and
  `tests/CMakeLists.txt`. A test that needs NetCDF as a whole goes in
  `NETCDF_REQUIRED_TESTS`.
- **LFortran:** a new importer of `NZ_STACK_MAX` needs the module-local
  `#ifdef LFORTRAN_PASSING` parameter its siblings carry, and the file count
  in `.github/workflows/portability.yml` follows.
- Declaration order, `use, intrinsic` and the deferred-length character
  assignment above are portability rules too: each one is a compiler that
  accepted the code silently while another rejected or miscompiled it.

## Submodules for Large Modules

- Use submodules to separate a module's interface from its implementation,
  so an implementation change doesn't recompile every user of the module.
  Worth it for the largest modules.

```fortran
! rdb_heavy_module.F90 - interface only
module rdb_heavy_module
   use rdb_constants, only: wp
   implicit none
   private
   public :: compute_expensive_thing

   interface
      module subroutine compute_expensive_thing(input, output)
         real(wp), intent(in) :: input(:)
         real(wp), intent(out) :: output(:)
      end subroutine compute_expensive_thing
   end interface
end module rdb_heavy_module

! rdb_heavy_module_impl.F90 - implementation
submodule (rdb_heavy_module) rdb_heavy_module_impl
contains
   module subroutine compute_expensive_thing(input, output)
      real(wp), intent(in) :: input(:)
      real(wp), intent(out) :: output(:)
      ! changes here don't recompile modules that use rdb_heavy_module
   end subroutine compute_expensive_thing
end submodule rdb_heavy_module_impl
```

## File Structure Template

```fortran
module rdb_ocean_example
   !! One-line description of the module.
   !!
   !! What it computes, the paper it implements, and where it sits in the
   !! step (which slot of ocean_state_t owns it).
   use rdb_constants, only: wp, H_VANISHED
   use rdb_multilayer_state, only: multilayer_state_t
   use pic_logger, only: global_logger
   implicit none
   private

   public :: example_t

   type :: example_t
      !! Per-run state of the example closure.
      real(wp), allocatable :: rate(:, :)
         !! Relaxation rate per cell (s^-1).
   contains
      procedure :: init => example_init
      procedure :: enter_data => example_enter_data
      procedure :: exit_data => example_exit_data
      procedure :: destroy => example_destroy
   end type example_t

contains

   subroutine example_init(this, nx, ny)
      !! Allocate and zero the rate field.
      class(example_t), intent(inout) :: this
      integer, intent(in) :: nx, ny
      allocate (this%rate(nx, ny), source=0.0_wp)
   end subroutine example_init

   subroutine example_enter_data(this)
      !! Map to the device: parent first, then components.
      class(example_t), intent(inout) :: this
      !$acc enter data copyin(this)
      !$acc enter data copyin(this%rate)
   end subroutine example_enter_data

   subroutine example_exit_data(this)
      !! Unmap: components first, parent last.
      class(example_t), intent(inout) :: this
      !$acc exit data delete(this%rate)
      !$acc exit data delete(this)
   end subroutine example_exit_data

   subroutine example_destroy(this)
      !! Free the owned arrays.
      class(example_t), intent(inout) :: this
      if (allocated(this%rate)) deallocate (this%rate)
   end subroutine example_destroy

end module rdb_ocean_example
```

## Enforcement

- **pre-commit** (`.pre-commit-config.yaml`): `fprettify` (layout),
  `fortitude` (language hygiene; `fortitude.toml`, line length 178; skips
  `tests/`, `benchmarks/` and `tools/`), and the project hooks:
  `no-mpi-in-rdb`, `dc-assumed-shape`, `dc-intrinsic-shadow`, `dc-transfer`,
  `decl-order`, `vanished-layer`, `openmp-portability`,
  `check-closure-matrix`, `test-regime-labels`, `check-bytes-accounting`.
  Run `pre-commit run --all` before committing.
- **CI**: the compiler matrix above, the MPI legs, the tier-2 stability sweep
  and the pre-commit hooks run on every PR.
- **Review**: everything the hooks can't check, using this guide and
  `CLAUDE.md`. The `pr-review` skill (`.claude/skills/pr-review/SKILL.md`)
  encodes the checklist; comment `/claude-review` on a PR to run it.
