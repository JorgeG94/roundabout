!! Ocean diagnostics manager.
module rdb_ocean_diag
   !! Owns the diagnostics pipeline for the ocean dynamical core.
   !! Built around the following design points (Phase 6 work):
   !!
   !!   1. **Pull model, not push.**  Kernels write to `state` and
   !!      stay diag-agnostic; this manager iterates a registry each
   !!      output step, fetches the source array (by pointer / by
   !!      enum tag), remaps onto the configured output vertical
   !!      grid, accumulates time means, and hands the ready buffer
   !!      to the I/O server.  Kernels never call `diag%send(...)`.
   !!
   !!   2. **Compute-side remap.**  Layer→z, layer→isopycnal, and
   !!      time-mean accumulation all happen on the compute rank
   !!      before hand-off.  The I/O server then only has to do
   !!      compression + netcdf write.  Keeps the I/O server lean
   !!      (no physics-aware code) at the cost of carrying the
   !!      remap target buffers on every compute rank — which is
   !!      cheap since they're shape `(nx, ny, nz_out)` not
   !!      `(nx_global, ny_global, nz_out)`.
   !!
   !!   3. **I/O server interaction.**  When a `diag_var_t`'s
   !!      cadence fires, this manager calls `io_server_send(buf,
   !!      meta)` on the per-rank send queue.  Backpressure (server
   !!      slower than compute) shows up as a hung send and is
   !!      logged.  Hand-off contract: the manager owns the buffer
   !!      lifetime; the I/O server gets a non-owning view.
   !!
   !! Phase 0e status: empty scaffold.  Components are declared so
   !! Phase 6 can fill in `register / step / send_ready` without
   !! restructuring the god state.  The bound procedures init/destroy
   !! are no-ops on the empty registry.
   !!
   !! Collaborator hand-off: this slot is independent of the
   !! dynamical core — once `ocean_diag_register` lands, every other
   !! kernel can be diagnosed without modification.  Good first task
   !! for a new contributor with FMS / diag_manager experience.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp
#else
   use rdb_constants, only: NZ_STACK_MAX, wp
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_diag_mask, only: diag_mask_t, diag_mask_destroy
   use pic_logger, only: logger => global_logger
   use, intrinsic :: iso_fortran_env, only: int64, real32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   implicit none
   private
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ocean_diag_t, diag_var_t, diag_fill_proc, diag_remap_proc
   public :: diag_emit_proc, ocean_diag_nc_stream_t
   public :: diag_spec_t, parse_diag_spec
   public :: diag_reduce_stats, diag_field_stats
      !! Exposed for `test_ocean_diag_reduce`: the emit-path statistics
      !! reduction is the only device kernel in this module, so it is unit
      !! tested directly (numerics + `mem:separate` residency) rather than
      !! only through `step`.

   abstract interface
      subroutine diag_fill_proc(state_handle, buf)
         !! Fill routine signature.  `state_handle` is polymorphic so
         !! this module stays decoupled from `ocean_state_t` (which
         !! itself composes `ocean_diag_t`).  Implementations live in
         !! `rdb_ocean_diag_fills` and do a `select type` cast to the
         !! concrete ocean state inside the body.  `buf` is the
         !! pre-allocated output buffer the manager owns.
         import :: wp
         implicit none
         class(*), intent(in) :: state_handle
         real(wp), intent(inout) :: buf(:, :, :)
      end subroutine diag_fill_proc

      subroutine diag_emit_proc(diag, ivar, t)
         !! Optional post-fire emit hook.  Called by `step` once per
         !! cadence-fire AFTER the log line; used by the NetCDF
         !! writer to append a slice for the var that just fired.
         !! Kept abstract so the diag module stays free of NetCDF /
         !! I/O-server deps — `rdb_ocean_diag_netcdf` provides the
         !! concrete implementation.
         import :: wp
         implicit none
         class(*), intent(inout) :: diag
         integer, intent(in) :: ivar
         real(wp), intent(in) :: t
      end subroutine diag_emit_proc

      subroutine diag_remap_proc(state_handle, z_out, layer_buf, output_buf, is_extensive)
         !! Vertical-remap routine signature.  Called by the manager
         !! after `fill` when a var's `output_vgrid` is not LAYER.
         !! `layer_buf` is the manager-owned native-grid scratch the
         !! fill wrote into; `output_buf` is the remap target on the
         !! configured output vgrid (z-levels / isopycnals / 2D).
         !! State-side data needed for the remap (e.g. `h_layer`,
         !! bathy) is fetched via `select type` on `state_handle`.
         !!
         !! Conservation contract:  `is_extensive` (forwarded by the
         !! manager from `diag_var_t%is_extensive`) selects the operator.
         !! `.false.` — INTENSIVE: the target value is the
         !! thickness-weighted average of the overlapping source layers.
         !! `.true.` — EXTENSIVE: the column integral is conservatively
         !! REDISTRIBUTED (Σ output == Σ layer_buf), for thickness-
         !! integrated quantities (hTr, h·KE, transports).
         import :: wp
         implicit none
         class(*), intent(in) :: state_handle
         real(wp), intent(in) :: z_out(:)
         real(wp), intent(in) :: layer_buf(:, :, :)
         real(wp), intent(inout) :: output_buf(:, :, :)
         logical, intent(in) :: is_extensive
      end subroutine diag_remap_proc
   end interface

   ! ---- Source-array location tags ----
   ! Identifies where in `ocean_state` the source data lives, so the
   ! remap step knows which stagger / shape to consume.
   integer, parameter, public :: DIAG_LOC_CENTER = 1  !! cell centre
   integer, parameter, public :: DIAG_LOC_FACE_X = 2  !! east face
   integer, parameter, public :: DIAG_LOC_FACE_Y = 3  !! north face
   integer, parameter, public :: DIAG_LOC_CORNER = 4  !! south-west corner

   ! ---- Vertical-grid tags ----
   ! Source vgrid: what the model stores.  Output vgrid: what we
   ! remap to before handing the buffer to the I/O server.
   integer, parameter, public :: DIAG_VGRID_LAYER = 1
      !! Model layers (k=1 bed, k=nz surface).
   integer, parameter, public :: DIAG_VGRID_Z_FIXED = 2
      !! Fixed z-levels (e.g., 50 standard depths).
   integer, parameter, public :: DIAG_VGRID_DENSITY = 3
      !! Isopycnal bins for watermass analysis.
   integer, parameter, public :: DIAG_VGRID_SURFACE = 4
      !! 2D surface slice (no remap needed).
   integer, parameter, public :: DIAG_VGRID_BOTTOM = 5
      !! 2D bed slice.
   integer, parameter, public :: DIAG_VGRID_ZSTAR = 6
      !! Fixed z*-levels: SSH-tracking stretched-depth output grid
      !! (each column's target depths scale with col_h / resting depth).
   integer, parameter, public :: DIAG_VGRID_SIGMA = 7
      !! Fixed sigma-levels: terrain-following fractional-depth output
      !! grid (target interface depth = sigma fraction * col_h).

   real(wp), parameter, public :: DIAG_MISSING_VALUE = 1.0e20_wp
      !! Sentinel written into remapped output cells that overlap no water
      !! (below-bottom / pinched-out) when vanished-target masking is on,
      !! and advertised as the NetCDF `_FillValue` / `missing_value`.

   ! ---- Time-operator tags ----
   integer, parameter, public :: DIAG_OP_INSTANT = 1  !! snapshot at cadence
   integer, parameter, public :: DIAG_OP_MEAN = 2  !! dt-weighted time mean over cadence
   integer, parameter, public :: DIAG_OP_MAX = 3
   integer, parameter, public :: DIAG_OP_MIN = 4
   integer, parameter, public :: DIAG_OP_INTEGRAL = 5
      !! Cumulative time integral over the cadence window —
      !! emits Σ(sample · dt) without dividing.  Used for budget
      !! closure: time-integrated fluxes through a surface or
      !! through an OBC face are the conserved quantity, not
      !! their per-step rate.

   ! ---- Unified diagnostic-selection spec ----
   integer, parameter, public :: DIAG_OP_UNSET = -1
      !! `diag_spec_t%time_op` sentinel: attribute not given => keep the
      !! diagnostic's canonical default time-operator.
   integer, parameter, public :: DIAG_COORD_UNSET = 0
      !! `diag_spec_t%coord` sentinel: no `:coord` attribute given => keep
      !! the diagnostic's default output vgrid (the global `&ocean_diag_nml
      !! vgrid`).  Set values are `DIAG_VGRID_*` (layer / z / zstar / sigma).

   type :: diag_spec_t
      !! One parsed entry of the `&ocean_diag_nml diags` selection list.
      !! Produced by `parse_diag_spec`; consumed by `register_default_diags`
      !! (canonical-default overrides / skips) and the derived-diagnostic
      !! orchestrator (`apply_diag_selection`).  Unset attributes carry
      !! sentinels so the consumer falls back to the canonical default.
      character(len=64) :: name = ""        !! diagnostic name (catalog or canonical)
      logical  :: off = .false.             !! `:off` => skip this diagnostic
      integer  :: time_op = DIAG_OP_UNSET   !! `:instant/mean/max/min` override
      real(wp) :: dt_out = -1.0_wp          !! `:<n>s/m/h/d` cadence override (<0 = unset)
      integer  :: coord = DIAG_COORD_UNSET  !! `:layer/z/zstar/sigma` output-vgrid override
   end type diag_spec_t

   type :: ocean_diag_nc_stream_t
      !! NetCDF output stream state.  Carried inline on `ocean_diag_t`
      !! so the diag module stays NetCDF-free (no `use netcdf` here);
      !! the actual file operations live in `rdb_ocean_diag_netcdf`
      !! which manipulates these fields.
      logical :: is_open = .false.
      integer :: ncid = -1
      character(len=256) :: filename = ""
      integer :: x_dimid = -1
      integer :: y_dimid = -1
      integer :: xtype = -1
         !! Resolved NetCDF element type for the DATA variables of this
         !! stream (`NC_WP` / `NC_R4` — see `rdb_io_netcdf`).  `-1` = unset;
         !! `open_stream` resolves it from `&ocean_diag_nml
         !! output_precision` ("double", the default, => `NC_WP` =>
         !! byte-identical output).  Held as a bare integer so this module
         !! stays NetCDF-free; only `rdb_ocean_diag_netcdf` interprets it.
         !!
         !! Scope note: this governs the DIAGNOSTIC stream only.  Restarts,
         !! gauges, coastal output, console conservation totals and
         !! checksums are unconditionally `NC_WP`/`wp` — there is
         !! deliberately no knob that can make a restart lossy.
      real(real32), allocatable :: stage(:, :, :)
         !! Host staging buffer for the fp64 -> fp32 conversion done
         !! immediately before `nf90_put_var`.  Allocated by `open_stream`
         !! ONLY when `xtype == NC_R4`, sized to the largest registered
         !! variable's `output_buffer`; unallocated (zero cost) on the
         !! default double-precision path.
   end type ocean_diag_nc_stream_t

   type :: diag_var_t
      !! Per-variable diagnostic record.  Owned by ocean_diag_t.
      character(len=64)  :: name = ""
         !! Short netcdf variable name.
      character(len=128) :: long_name = ""
         !! CF-compliant long_name attribute.
      character(len=32)  :: units = ""
      character(len=64)  :: standard_name = ""

      ! ---- Source binding ----
      ! Procedure pointer to the per-var fill routine.  Set at
      ! `register` time by the slot that owns the source data; the
      ! manager invokes it on cadence-fire to populate the buffer.
      ! `nopass` because the registry holds the bind, not the var.
      procedure(diag_fill_proc), pointer, nopass :: fill => null()
      procedure(diag_remap_proc), pointer, nopass :: remap => null()
         !! Optional vertical-remap routine.  Set at register time
         !! when `output_vgrid /= DIAG_VGRID_LAYER`; null otherwise.
      integer :: source_loc = DIAG_LOC_CENTER
      integer :: source_vgrid = DIAG_VGRID_LAYER

      ! ---- Output binding ----
      integer  :: output_vgrid = DIAG_VGRID_LAYER
      integer  :: time_op = DIAG_OP_MEAN
      real(wp) :: dt_out = 3600.0_wp
         !! Output cadence (s).  Manager keeps a per-var counter and
         !! fires when the accumulated dt passes the threshold.

      ! ---- Buffers ----
      ! `output_buffer` is the post-remap buffer that gets shipped
      ! to the I/O server.  `layer_buffer` is the manager-owned
      ! native-grid scratch the `fill` routine writes into when a
      ! remap step is required; for `output_vgrid == LAYER` the
      ! fill writes directly to `output_buffer` and `layer_buffer`
      ! stays unallocated.  `accumulator` is the running-sum buffer
      ! used when `time_op /= DIAG_OP_INSTANT`.
      real(wp), allocatable :: output_buffer(:, :, :)
      real(wp), allocatable :: layer_buffer(:, :, :)
      real(wp), allocatable :: accumulator(:, :, :)
      integer :: n_accum = 0
         !! Number of contributions in the current accumulator window.
      real(wp) :: dt_accum = 0.0_wp
         !! Wall time accumulated into this var's window (s).

      ! ---- Region restriction ----
      type(diag_mask_t), allocatable :: mask
         !! Optional region mask.  When allocated, `fold_sample`
         !! multiplies each sample by `mask%weight(i, j)` before
         !! folding into the accumulator — cells outside the region
         !! contribute zero.  Output buffer shape is unchanged (full
         !! domain with zeros outside the mask); scalar-aggregating
         !! reductions land with the Phase D budget plumbing or as
         !! a Phase C v2 follow-on.

      ! ---- Vertical quantity kind ----
      logical :: is_extensive = .false.
         !! `.false.` (default): field is INTENSIVE — per-unit-thickness
         !! quantity like temperature, salinity, velocity, density.
         !! When remapped to a non-LAYER vgrid, the column value is
         !! weight-averaged across overlapping source layers.
         !!
         !! `.true.`: field is EXTENSIVE — already thickness-integrated,
         !! e.g. `hTr` (tracer · m), KE per layer (h · 0.5 · |u|²),
         !! transport per layer (h · u).  When remapped, the column
         !! values must be conservatively REDISTRIBUTED across the
         !! target layers (sum preserved), not averaged.
         !!
         !! Phase B v1 status: the existing `remap_layer_to_z` is
         !! intensive-only.  Extensive remap (and the corresponding
         !! split inside `diag_remap_proc`) lands with isopycnal /
         !! density-bin remap in Phase E.  Flag is here now so calling
         !! code can declare intent and the upgrade is non-breaking.

      ! ---- Enable gate ----
      logical :: enabled = .true.
         !! When `.false.` the dispatcher skips this var entirely — no
         !! `fill_*` kernel, no accumulator fold, no emit — so an
         !! unrequested diagnostic costs zero GPU work.  Lets a run turn
         !! off individual default diagnostics it does not want.  Default
         !! `.true.` ⇒ every registered var runs (bit-identical).

      ! ---- Vanished-target masking ----
      logical :: has_missing = .false.
         !! When `.true.` the remap fills target cells that overlap no water
         !! (below-bottom / pinched-out in a shallow column) with
         !! `DIAG_MISSING_VALUE` instead of 0, and the NetCDF writer tags the
         !! variable with a `_FillValue` / `missing_value` attribute.  Set at
         !! register time for non-LAYER diagnostics when masking is enabled
         !! (`&ocean_diag_nml mask_vanished_layers`).  Default `.false.` =>
         !! below-bottom cells read 0 (bit-identical to the legacy writer).

      ! ---- I/O server binding ----
      integer :: stream_id = 0
         !! ID of the output stream this var ships to (Phase 6 wires
         !! the I/O server stream table).

      ! ---- Per-var NetCDF state ----
      ! Populated by `rdb_ocean_diag_netcdf` when the stream is
      ! opened; consumed by the NetCDF emit hook on each fire.
      integer :: nc_varid = -1
      integer :: nc_time_dimid = -1
      integer :: nc_time_varid = -1
      integer :: nc_z_dimid = -1
         !! -1 for 2D vars; set for 3D vars when the stream is opened.
      integer :: nc_time_index = 0
         !! Number of slices already written for this var (1-indexed
         !! position of the NEXT write).

      ! ---- In-memory access generation (P7) ----
      integer :: fire_count = 0
         !! Number of times this var's `output_buffer` has been
         !! refreshed (filled/folded/finalised + pulled host-ward) since
         !! registration. Incremented unconditionally in `ocean_diag_step`
         !! at the SAME point as the per-fire `update self` — independent
         !! of whether a NetCDF stream is open (`nc_time_index` stays 0
         !! with diagnostics disabled or output suppressed, so it cannot
         !! serve this role). This is the `generation` the C ABI's
         !! `rdb_ocean_get_diagnostic_ptr` returns: a Python `Field`
         !! re-checks it on access and knows `output_buffer` is already
         !! host-current the instant it changes (no separate refresh call
         !! needed — the pull above already happened synchronously).
   contains
      procedure, non_overridable :: bytes => diag_var_bytes
   end type diag_var_t

   type :: ocean_diag_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      logical :: on_device = .false.
         !! `.true.` between `enter_data` and `exit_data`.  The emit-time
         !! statistics reduction reads `output_buffer` where it actually
         !! lives: on the device via `do concurrent ... reduce` when mapped,
         !! on the host otherwise (unit tests that skip `enter_data`).  A
         !! device read of an unmapped buffer under `-gpu=mem:separate`
         !! would return garbage silently, so this is not optional.

      ! ---- Registry ----
      integer :: nvars = 0
         !! Live count of registered diagnostics.
      integer :: nvars_max = 0
         !! Capacity of `vars(:)`.  Grown via reallocation on
         !! `register` overflow.
      type(diag_var_t), allocatable :: vars(:)

      ! ---- Output vertical grids ----
      ! Configured at init from the namelist.  The manager looks up
      ! `vars(i)%output_vgrid` and dispatches the remap onto the
      ! matching array.  Phase 6 wires the remap kernels.
      integer  :: nz_out = 0
      real(wp), allocatable :: z_out(:)
         !! Output z-levels (m, positive up; surface at index nz_out).
      integer  :: n_rho_out = 0
      real(wp), allocatable :: rho_out(:)
         !! Output isopycnal bin edges (kg/m^3).
      integer  :: n_sigma_out = 0
      real(wp), allocatable :: sigma_out(:)
         !! Output sigma levels: cumulative fractions (0..1, shallow->deep).
      integer  :: n_zstar_out = 0
      real(wp), allocatable :: zstar_out(:)
         !! Output z* reference interface depths (m, positive-down; deepest =
         !! reference total depth H_ref).  Per-column grid stretched by
         !! col_h / H_ref (SSH-tracking).

      ! ---- I/O-server send buffer ----
      ! Pre-allocated, reused across diag pushes.  Sized to the
      ! largest registered variable's `output_buffer`.  Sending uses
      ! a zero-copy `c_loc` view into this buffer.
      real(wp), allocatable :: send_buf(:, :, :)

      ! ---- Triggering scalars ----
      logical :: enabled = .true.
         !! Master switch.  Phase 6 reads from namelist.
      real(wp) :: dt_last_eval = 0.0_wp
         !! Wall time at last evaluation pass (s).

      ! ---- NetCDF stream ----
      ! Inline state (no `use netcdf` here); the writer module
      ! `rdb_ocean_diag_netcdf` opens / writes / closes against this.
      type(ocean_diag_nc_stream_t) :: nc_stream
      procedure(diag_emit_proc), pointer, nopass :: emit_post_fire => null()
         !! Optional post-fire hook the manager calls after the log
         !! line.  `rdb_ocean_diag_netcdf::open_stream` binds this to
         !! the NetCDF writer; unbound = log-only behaviour.

      type(hgrid_t) :: grid
         !! Cached grid (scalar-only struct, cheap to copy).  Derived
         !! diagnostic fills read `grid%nghost` from here (and the
         !! per-cell spacing from `state%metrics` directly — design D5);
         !! they only see `state_handle` so the grid has to be reachable
         !! through the state composition.
   contains
      procedure, non_overridable :: init => ocean_diag_init
      procedure, non_overridable :: destroy => ocean_diag_destroy
      procedure, non_overridable :: register => ocean_diag_register
      procedure, non_overridable :: disable => ocean_diag_disable
      procedure, non_overridable :: is_registered => ocean_diag_is_registered
      procedure, non_overridable :: step => ocean_diag_step
      procedure, non_overridable :: set_output_z_levels => ocean_diag_set_output_z_levels
      procedure, non_overridable :: set_output_density_levels => ocean_diag_set_output_density_levels
      procedure, non_overridable :: set_output_sigma_levels => ocean_diag_set_output_sigma_levels
      procedure, non_overridable :: set_output_zstar_levels => ocean_diag_set_output_zstar_levels
      procedure, non_overridable :: enter_data => ocean_diag_enter_data
      procedure, non_overridable :: exit_data => ocean_diag_exit_data
      procedure, non_overridable :: bytes => ocean_diag_bytes
   end type ocean_diag_t

   integer, parameter :: INITIAL_CAPACITY = 16

contains

   function parse_diag_spec(spec) result(specs)
      !! Parse the `&ocean_diag_nml diags` selection string into structured
      !! entries.  Entries are whitespace/comma-separated; within an entry,
      !! colon-separated attributes are self-identifying (order-free):
      !!
      !!   * `off`                        -> skip the diagnostic
      !!   * `instant`/`mean`/`max`/`min` -> time-operator override
      !!   * `<int><unit>` (s/m/h/d)      -> cadence override
      !!
      !! The first colon field is always the diagnostic name (case-sensitive,
      !! matched against registered / catalog names).  Attribute matching is
      !! case-insensitive.  An empty/blank string yields a zero-length array
      !! (the bit-identical default).  An unrecognised attribute fails loud.
      character(len=*), intent(in) :: spec
      type(diag_spec_t), allocatable :: specs(:)
      character(len=len(spec)) :: buf
      integer :: i, j, n, ntok, k

      ! Normalise commas to spaces so both separators work.
      buf = spec
      do i = 1, len(buf)
         if (buf(i:i) == ",") buf(i:i) = " "
      end do
      n = len_trim(buf)

      ! Pass 1: count whitespace-delimited tokens.
      ntok = 0
      i = 1
      do while (i <= n)
         if (buf(i:i) == " ") then
            i = i + 1
            cycle
         end if
         ntok = ntok + 1
         do while (i <= n)
            if (buf(i:i) == " ") exit
            i = i + 1
         end do
      end do

      allocate (specs(ntok))
      if (ntok == 0) return

      ! Pass 2: parse each token into name + attributes.
      k = 0
      i = 1
      do while (i <= n)
         if (buf(i:i) == " ") then
            i = i + 1
            cycle
         end if
         j = i
         do while (j <= n)
            if (buf(j:j) == " ") exit
            j = j + 1
         end do
         k = k + 1
         call parse_one_spec_token(buf(i:j - 1), specs(k))
         i = j + 1
      end do
   end function parse_diag_spec

   subroutine parse_one_spec_token(tok, s)
      !! Parse a single `name[:attr]...` token into a `diag_spec_t`.
      !! Fails loud on an unrecognised attribute.
      character(len=*), intent(in) :: tok
      type(diag_spec_t), intent(out) :: s
      integer :: p, q, m
      real(wp) :: secs
      logical :: ok
      character(len=:), allocatable :: attr

      m = len_trim(tok)
      ! First colon field is the name.
      p = index(tok(1:m), ":")
      if (p == 0) then
         s%name = tok(1:m)
         return
      end if
      s%name = tok(1:p - 1)

      ! Remaining colon-separated fields are attributes.
      p = p + 1
      do while (p <= m)
         q = index(tok(p:m), ":")
         if (q == 0) then
            attr = lower_ascii(tok(p:m))
            p = m + 1
         else
            attr = lower_ascii(tok(p:p + q - 2))
            p = p + q
         end if
         if (len_trim(attr) == 0) cycle

         select case (attr)
         case ("off")
            s%off = .true.
         case ("instant")
            s%time_op = DIAG_OP_INSTANT
         case ("mean")
            s%time_op = DIAG_OP_MEAN
         case ("max")
            s%time_op = DIAG_OP_MAX
         case ("min")
            s%time_op = DIAG_OP_MIN
         case ("integral")
            s%time_op = DIAG_OP_INTEGRAL
         case ("layer")
            s%coord = DIAG_VGRID_LAYER
         case ("z")
            s%coord = DIAG_VGRID_Z_FIXED
         case ("zstar", "z*")
            s%coord = DIAG_VGRID_ZSTAR
         case ("sigma")
            s%coord = DIAG_VGRID_SIGMA
         case ("density", "rho")
            s%coord = DIAG_VGRID_DENSITY
         case default
            call parse_cadence_attr(attr, secs, ok)
            if (ok) then
               s%dt_out = secs
            else
               call logger%error("unknown diag attribute '"//attr// &
                                 "' in spec token '"//trim(tok)// &
                                 "'; valid: off | instant|mean|max|min|integral | "// &
                                 "layer|z|zstar|sigma|density | <int>s/m/h/d")
               error stop "parse_diag_spec: unknown attribute"
            end if
         end select
      end do
   end subroutine parse_one_spec_token

   pure subroutine parse_cadence_attr(s, secs, ok)
      !! Parse a cadence attribute `<int><unit>` (unit s/m/h/d) to seconds.
      !! `ok=.false.` if `s` is not a well-formed positive cadence.
      character(len=*), intent(in) :: s
      real(wp), intent(out) :: secs
      logical, intent(out) :: ok
      integer :: n, ios, val
      character :: unit
      character(len=128) :: ferr

      secs = -1.0_wp
      ok = .false.
      n = len_trim(s)
      if (n < 2) return
      unit = s(n:n)
      read (s(1:n - 1), *, iostat=ios, iomsg=ferr) val
      if (ios /= 0 .or. val <= 0) return
      select case (unit)
      case ("s")
         secs = real(val, wp)
      case ("m")
         secs = real(val, wp)*60.0_wp
      case ("h")
         secs = real(val, wp)*3600.0_wp
      case ("d")
         secs = real(val, wp)*86400.0_wp
      case default
         return
      end select
      ok = .true.
   end subroutine parse_cadence_attr

   pure function lower_ascii(s) result(out)
      !! ASCII lowercase a string (attribute matching is case-insensitive).
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            out(i:i) = achar(c + 32)
         else
            out(i:i) = s(i:i)
         end if
      end do
   end function lower_ascii

   subroutine ocean_diag_init(this, grid)
      !! Allocate an empty registry sized at `INITIAL_CAPACITY` slots.
      !! Subsequent `register` calls grow the array via doubling.
      class(ocean_diag_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      this%grid = grid
      this%nvars = 0
      this%nvars_max = INITIAL_CAPACITY
      allocate (this%vars(this%nvars_max))
      this%is_init = .true.
   end subroutine ocean_diag_init

   subroutine ocean_diag_destroy(this)
      class(ocean_diag_t), intent(inout) :: this
      integer :: i
      this%is_init = .false.
      if (allocated(this%vars)) then
         do i = 1, this%nvars
            if (allocated(this%vars(i)%output_buffer)) deallocate (this%vars(i)%output_buffer)
            if (allocated(this%vars(i)%layer_buffer)) deallocate (this%vars(i)%layer_buffer)
            if (allocated(this%vars(i)%accumulator)) deallocate (this%vars(i)%accumulator)
            if (allocated(this%vars(i)%mask)) then
               call diag_mask_destroy(this%vars(i)%mask)
               deallocate (this%vars(i)%mask)
            end if
            nullify (this%vars(i)%fill)
            nullify (this%vars(i)%remap)
         end do
         deallocate (this%vars)
      end if
      this%nvars = 0
      this%nvars_max = 0
      if (allocated(this%z_out)) deallocate (this%z_out)
      if (allocated(this%rho_out)) deallocate (this%rho_out)
      if (allocated(this%sigma_out)) deallocate (this%sigma_out)
      if (allocated(this%zstar_out)) deallocate (this%zstar_out)
      if (allocated(this%send_buf)) deallocate (this%send_buf)
      ! Host-only fp32 staging buffer for the single-precision diag stream
      ! (never device-mapped, so no exit_data pairing).  Normally released
      ! by `close_stream`; freed here too so a destroy without a close
      ! does not leak.
      if (allocated(this%nc_stream%stage)) deallocate (this%nc_stream%stage)
   end subroutine ocean_diag_destroy

   subroutine ocean_diag_enter_data(this)
      !! Attach each registered var's per-buffer allocatables to the device.
      !! Called by `ocean_state_enter_data` after `register_default_diags`
      !! has populated `vars(:)` — buffers are sized at register time, so
      !! the descriptors here are valid.
      !!
      !! Per-var bookkeeping (dt_accum, n_accum, time_op, ...) stays on
      !! host; the manager step loop is host code that dispatches device
      !! kernels via the outer-shim + flat-impl pattern.  Only the bulk
      !! buffers (output_buffer, accumulator, layer_buffer, mask weight)
      !! need device residency.
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! map base is the heap object, not a polymorphic stack box (AMD
      !! libomptarget cross-slot-overlap fix).
      class(ocean_diag_t), intent(inout) :: this
      select type (this)
      type is (ocean_diag_t)
         call ocean_diag_enter_data_impl(this)
      end select
   end subroutine ocean_diag_enter_data

   subroutine ocean_diag_enter_data_impl(this)
      type(ocean_diag_t), intent(inout) :: this
      integer :: i
      if (.not. this%is_init) return
      if (allocated(this%vars)) then
         do i = 1, this%nvars
            associate (v => this%vars(i))
               if (allocated(v%output_buffer)) then
                  !$acc enter data copyin(v%output_buffer)
               end if
               if (allocated(v%accumulator)) then
                  !$acc enter data copyin(v%accumulator)
               end if
               if (allocated(v%layer_buffer)) then
                  !$acc enter data copyin(v%layer_buffer)
               end if
               if (allocated(v%mask)) then
                  !$acc enter data copyin(v%mask)
                  associate (m => v%mask)
                     if (allocated(m%weight)) then
                        !$acc enter data copyin(m%weight)
                     end if
                  end associate
               end if
            end associate
         end do
      end if
      if (allocated(this%z_out)) then
         !$acc enter data copyin(this%z_out)
      end if
      if (allocated(this%rho_out)) then
         !$acc enter data copyin(this%rho_out)
      end if
      if (allocated(this%sigma_out)) then
         !$acc enter data copyin(this%sigma_out)
      end if
      if (allocated(this%zstar_out)) then
         !$acc enter data copyin(this%zstar_out)
      end if
      if (allocated(this%send_buf)) then
         !$acc enter data copyin(this%send_buf)
      end if
      this%on_device = .true.
   end subroutine ocean_diag_enter_data_impl

   subroutine ocean_diag_exit_data(this)
      !! Detach in reverse order of enter_data.  Idempotency-safe via
      !! `is_init` gate — repeated calls without intervening enter_data
      !! become no-ops once the manager is destroyed.
      class(ocean_diag_t), intent(inout) :: this
      select type (this)
      type is (ocean_diag_t)
         call ocean_diag_exit_data_impl(this)
      end select
   end subroutine ocean_diag_exit_data

   subroutine ocean_diag_exit_data_impl(this)
      type(ocean_diag_t), intent(inout) :: this
      integer :: i
      if (.not. this%is_init) return
      if (allocated(this%send_buf)) then
         !$acc exit data delete(this%send_buf)
      end if
      if (allocated(this%zstar_out)) then
         !$acc exit data delete(this%zstar_out)
      end if
      if (allocated(this%sigma_out)) then
         !$acc exit data delete(this%sigma_out)
      end if
      if (allocated(this%rho_out)) then
         !$acc exit data delete(this%rho_out)
      end if
      if (allocated(this%z_out)) then
         !$acc exit data delete(this%z_out)
      end if
      if (allocated(this%vars)) then
         do i = this%nvars, 1, -1
            associate (v => this%vars(i))
               if (allocated(v%mask)) then
                  associate (m => v%mask)
                     if (allocated(m%weight)) then
                        !$acc exit data delete(m%weight)
                     end if
                  end associate
                  !$acc exit data delete(v%mask)
               end if
               if (allocated(v%layer_buffer)) then
                  !$acc exit data delete(v%layer_buffer)
               end if
               if (allocated(v%accumulator)) then
                  !$acc exit data delete(v%accumulator)
               end if
               if (allocated(v%output_buffer)) then
                  !$acc exit data delete(v%output_buffer)
               end if
            end associate
         end do
      end if
      this%on_device = .false.
   end subroutine ocean_diag_exit_data_impl

   subroutine ocean_diag_register(this, name, units, fill, n1, n2, n3, &
                                  long_name, standard_name, time_op, dt_out, &
                                  output_vgrid, remap, mask, is_extensive, has_missing)
      !! Register a new diagnostic variable.  Grows the registry via
      !! capacity doubling on overflow.  Buffer allocation depends on
      !! `output_vgrid`:
      !!   * `LAYER` (default): one `output_buffer(n1, n2, n3)` —
      !!     `fill` writes directly into it.
      !!   * `Z_FIXED` (or any non-LAYER target): two buffers —
      !!     `layer_buffer(n1, n2, n3)` for the fill, plus
      !!     `output_buffer(n1, n2, this%nz_out)` for the remapped
      !!     result.  Caller must have configured `nz_out` via
      !!     `set_output_z_levels` first, and bind a `remap` proc.
      !! Caller binds `fill` to a routine that knows how to populate
      !! the layer-native buffer from the state handle.
      !!
      !! `mask` (optional): restricts accumulation to the region
      !! where `mask%weight > 0`.  See `diag_mask_t` builders in
      !! `rdb_ocean_diag_mask`.
      class(ocean_diag_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: units
      procedure(diag_fill_proc) :: fill
      integer, intent(in) :: n1, n2, n3
      character(len=*), intent(in), optional :: long_name, standard_name
      integer, intent(in), optional :: time_op
      real(wp), intent(in), optional :: dt_out
      integer, intent(in), optional :: output_vgrid
      procedure(diag_remap_proc), optional :: remap
      type(diag_mask_t), intent(in), optional :: mask
      logical, intent(in), optional :: is_extensive
      logical, intent(in), optional :: has_missing
      type(diag_var_t), allocatable :: tmp(:)
      integer :: i, ovgrid, nzout

      if (.not. this%is_init) return

      if (this%nvars == this%nvars_max) then
         allocate (tmp(2*this%nvars_max))
         do i = 1, this%nvars
            tmp(i) = this%vars(i)
         end do
         call move_alloc(tmp, this%vars)
         this%nvars_max = size(this%vars)
      end if

      ovgrid = DIAG_VGRID_LAYER
      if (present(output_vgrid)) ovgrid = output_vgrid

      this%nvars = this%nvars + 1
      associate (v => this%vars(this%nvars))
         v%name = name
         v%units = units
         v%fill => fill
         v%output_vgrid = ovgrid
         if (present(remap)) v%remap => remap
         if (present(long_name)) v%long_name = long_name
         if (present(standard_name)) v%standard_name = standard_name
         if (present(time_op)) v%time_op = time_op
         if (present(dt_out)) v%dt_out = dt_out
         if (present(has_missing)) v%has_missing = has_missing
         if (present(is_extensive)) v%is_extensive = is_extensive

         if (ovgrid == DIAG_VGRID_LAYER) then
            allocate (v%output_buffer(n1, n2, n3), source=0.0_wp)
         else
            allocate (v%layer_buffer(n1, n2, n3), source=0.0_wp)
            select case (ovgrid)
            case (DIAG_VGRID_DENSITY)
               nzout = this%n_rho_out
            case (DIAG_VGRID_SIGMA)
               nzout = this%n_sigma_out
            case (DIAG_VGRID_ZSTAR)
               nzout = this%n_zstar_out
            case default
               nzout = this%nz_out
            end select
            if (nzout <= 0) nzout = n3
            allocate (v%output_buffer(n1, n2, nzout), source=0.0_wp)
         end if

         if (v%time_op /= DIAG_OP_INSTANT) then
            allocate (v%accumulator(size(v%output_buffer, 1), &
                                    size(v%output_buffer, 2), &
                                    size(v%output_buffer, 3)))
            call reset_accumulator(v)
         end if
         v%n_accum = 0
         v%dt_accum = 0.0_wp

         if (present(mask)) then
            allocate (v%mask, source=mask)
         end if
      end associate
   end subroutine ocean_diag_register

   subroutine ocean_diag_disable(this, name)
      !! Turn OFF the registered diagnostic `name` so the dispatcher skips
      !! it entirely (no fill, no fold, no emit).  Fail-loud if `name`
      !! matches no registered var — a typo must not silently leave a
      !! diagnostic running.  Lists the registered names on abort.
      class(ocean_diag_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: i
      character(len=1024) :: avail

      do i = 1, this%nvars
         if (trim(this%vars(i)%name) == trim(name)) then
            this%vars(i)%enabled = .false.
            return
         end if
      end do

      avail = ""
      do i = 1, this%nvars
         avail = trim(avail)//" "//trim(this%vars(i)%name)
      end do
      call logger%error("ocean_diag%disable: no registered diagnostic named '"// &
                        trim(name)//"' — registered:"//trim(avail))
      error stop "ocean_diag%disable: unknown diagnostic name"
   end subroutine ocean_diag_disable

   pure function ocean_diag_is_registered(this, name) result(yes)
      !! `.true.` iff a diagnostic named `name` is registered (enabled or
      !! not).  Registration state only — says nothing about whether it
      !! will actually fire (see `enabled`); a `disable`d diagnostic is
      !! still registered and this returns `.true.` for it.
      class(ocean_diag_t), intent(in) :: this
      character(len=*), intent(in) :: name
      logical :: yes
      integer :: i
      yes = .false.
      do i = 1, this%nvars
         if (trim(this%vars(i)%name) == trim(name)) then
            yes = .true.
            return
         end if
      end do
   end function ocean_diag_is_registered

   pure subroutine reset_accumulator(v)
      !! Zero (MEAN), -huge (MAX), or +huge (MIN) the accumulator —
      !! so the first accumulate step seeds correctly.  Runs on device
      !! via the flat-impl shim because `v%accumulator` is reached
      !! through the `vars(:)` array-of-derived-types indirection that
      !! NVHPC can't follow inside a `do concurrent`.
      type(diag_var_t), intent(inout) :: v
      real(wp) :: seed
      if (.not. allocated(v%accumulator)) return
      select case (v%time_op)
      case (DIAG_OP_MAX)
         seed = -huge(0.0_wp)
      case (DIAG_OP_MIN)
         seed = huge(0.0_wp)
      case default
         seed = 0.0_wp
      end select
      call fill_buffer_impl(v%accumulator, seed)
   end subroutine reset_accumulator

   pure subroutine fill_buffer_impl(buf, val)
      !! Device-side scalar fill.  Used by reset_accumulator and (when
      !! `idx_temperature` / `idx_salinity` is unset) the manager could
      !! also seed `output_buffer` via this; today only the accumulator
      !! reset goes through here.
      ! assumed-shape-ok: diag accumulator reset — fires once per output frame.
      real(wp), intent(inout) :: buf(:, :, :)
      real(wp), intent(in)    :: val
      integer :: i, j, k, nx, ny, nz
      nx = size(buf, 1)
      ny = size(buf, 2)
      nz = size(buf, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = val
      end do
   end subroutine fill_buffer_impl

   subroutine ocean_diag_set_output_z_levels(this, z, ierr)
      !! Configure the fixed-z output grid.  Stored copy; the original
      !! `z` array is not retained.  Must be called BEFORE any
      !! `register` with `output_vgrid == DIAG_VGRID_Z_FIXED` so the
      !! manager knows the output buffer shape.
      class(ocean_diag_t), intent(inout) :: this
      real(wp), intent(in) :: z(:)
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`) when `size(z) >
         !! NZ_STACK_MAX`, when present; absent behaves as today
         !! (`error stop`). (F5 residual, P2.4)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. this%is_init) return
      ! The conservative remap pads to max(nz, nz_out) in NZ_STACK_MAX-sized
      ! per-column stack buffers; more output levels than that would overrun
      ! them (silent device illegal-address).  Fail loud instead.
      if (size(z) > NZ_STACK_MAX) then
         call fail("set_output_z_levels: number of z-levels exceeds "// &
                   "NZ_STACK_MAX; raise NZ_STACK_MAX or use fewer levels", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (allocated(this%z_out)) deallocate (this%z_out)
      this%nz_out = size(z)
      allocate (this%z_out(this%nz_out), source=z)
   end subroutine ocean_diag_set_output_z_levels

   subroutine ocean_diag_set_output_density_levels(this, rho, ierr)
      !! Configure the isopycnal (DENSITY) output grid — monotone-increasing
      !! target potential densities (kg/m³).  Stored copy; must be called
      !! BEFORE any `register` with `output_vgrid == DIAG_VGRID_DENSITY` so
      !! the manager knows the output buffer shape (one cell per target).
      class(ocean_diag_t), intent(inout) :: this
      real(wp), intent(in) :: rho(:)
      integer, intent(out), optional :: ierr
      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. this%is_init) return
      if (size(rho) > NZ_STACK_MAX) then
         call fail("set_output_density_levels: number of density bins "// &
                   "exceeds NZ_STACK_MAX; raise NZ_STACK_MAX or use fewer bins", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (allocated(this%rho_out)) deallocate (this%rho_out)
      this%n_rho_out = size(rho)
      allocate (this%rho_out(this%n_rho_out), source=rho)
   end subroutine ocean_diag_set_output_density_levels

   subroutine ocean_diag_set_output_sigma_levels(this, sigma, ierr)
      !! Configure the terrain-following (SIGMA) output grid — cumulative
      !! sigma fractions (0..1, monotone shallow->deep).  Stored copy; must
      !! be called BEFORE any `register` with `output_vgrid ==
      !! DIAG_VGRID_SIGMA` so the manager knows the output buffer shape.
      class(ocean_diag_t), intent(inout) :: this
      real(wp), intent(in) :: sigma(:)
      integer, intent(out), optional :: ierr
      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. this%is_init) return
      if (size(sigma) > NZ_STACK_MAX) then
         call fail("set_output_sigma_levels: number of sigma levels "// &
                   "exceeds NZ_STACK_MAX; raise NZ_STACK_MAX or use fewer levels", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (allocated(this%sigma_out)) deallocate (this%sigma_out)
      this%n_sigma_out = size(sigma)
      allocate (this%sigma_out(this%n_sigma_out), source=sigma)
   end subroutine ocean_diag_set_output_sigma_levels

   subroutine ocean_diag_set_output_zstar_levels(this, zstar, ierr)
      !! Configure the SSH-tracking (ZSTAR) output grid — reference interface
      !! depths (m, positive-down, monotone shallow->deep; deepest = H_ref).
      !! Stored copy; must be called BEFORE any `register` with
      !! `output_vgrid == DIAG_VGRID_ZSTAR` so the manager knows the buffer
      !! shape.
      class(ocean_diag_t), intent(inout) :: this
      real(wp), intent(in) :: zstar(:)
      integer, intent(out), optional :: ierr
      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. this%is_init) return
      if (size(zstar) > NZ_STACK_MAX) then
         call fail("set_output_zstar_levels: number of z* levels "// &
                   "exceeds NZ_STACK_MAX; raise NZ_STACK_MAX or use fewer levels", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (allocated(this%zstar_out)) deallocate (this%zstar_out)
      this%n_zstar_out = size(zstar)
      allocate (this%zstar_out(this%n_zstar_out), source=zstar)
   end subroutine ocean_diag_set_output_zstar_levels

   subroutine ocean_diag_step(this, state_handle, dt, t)
      !! Advance every registered variable.  Behaviour by `time_op`:
      !!
      !!   * `INSTANT`: skip until `dt_accum >= dt_out`, then fill once
      !!     (snapshot) and emit.
      !!   * `MEAN` / `INTEGRAL` / `MAX` / `MIN`: fill EVERY step, fold
      !!     the sample into the accumulator.  When `dt_accum >= dt_out`,
      !!     finalise the accumulator into `output_buffer` (see
      !!     `finalise_accumulator`), emit, reset.
      !!
      !! Fill writes into `output_buffer` for LAYER vgrid or into
      !! `layer_buffer` (then remap to `output_buffer`) for Z_FIXED.
      !! The post-remap value is the unit of accumulation — so the
      !! mean of T at fixed z-levels is computed in z-space, not
      !! layer-space (matters when h_layer drifts).
      !!
      !! Time stamping (CF-1.8 convention):
      !!   * INSTANT, MAX, MIN -> emit at `t` (end of window).
      !!   * MEAN, INTEGRAL    -> emit at the centre of the window,
      !!     `t - dt_accum / 2`.  The centre lines up with the
      !!     temporal centroid of the time-weighted average so
      !!     downstream tools (CF-aware analysis) plot it correctly.
      !!
      !! The emitted `[diag]` line's min / max / mean are taken over the
      !! FINITE cells only (`diag_field_stats`) — land columns and
      !! vanished layers carry the NaN missing-data sentinel — and the
      !! line gains a `missing=<excluded>/<total>` suffix whenever any
      !! cell was excluded.
      class(ocean_diag_t), intent(inout) :: this
      class(*), intent(in) :: state_handle
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: t
      integer :: i, n_valid, n_total
      real(wp) :: vmin, vmax, vmean, t_emit
      logical :: fire
      character(len=256) :: line
      character(len=64) :: missing_tag

      if (.not. this%is_init) return
      if (.not. this%enabled) return

      do i = 1, this%nvars
         associate (v => this%vars(i))
            v%dt_accum = v%dt_accum + dt
            fire = (v%dt_accum + 1.0e-9_wp >= v%dt_out)

            if (.not. associated(v%fill)) cycle
            if (.not. v%enabled) cycle

            select case (v%time_op)
            case (DIAG_OP_INSTANT)
               if (.not. fire) cycle
               call fill_and_remap(v, state_handle, this%z_out, this%rho_out, &
                                   this%sigma_out, this%zstar_out)
            case (DIAG_OP_MEAN, DIAG_OP_INTEGRAL, DIAG_OP_MAX, DIAG_OP_MIN)
               call fill_and_remap(v, state_handle, this%z_out, this%rho_out, &
                                   this%sigma_out, this%zstar_out)
               call fold_sample(v, dt)
               if (.not. fire) cycle
               call finalise_accumulator(v)
            case default
               error stop "ocean_diag%step: unknown v%time_op"
            end select

            select case (v%time_op)
            case (DIAG_OP_MEAN, DIAG_OP_INTEGRAL)
               t_emit = t - 0.5_wp*v%dt_accum
            case default
               t_emit = t
            end select

            ! Single device→host pull at cadence fire — the fills, fold,
            ! and finalise above all ran on device, so `output_buffer`
            ! lives on the device.  Pull it once here so both the log-line
            ! reductions and the NetCDF write below see fresh host data.
            !
            ! `if_present` keeps lightweight unit tests (which skip
            ! `ocean_state_enter_data` and run on host) working — when
            ! the buffer hasn't been attached, the directive is a no-op
            ! and the host buffer is already current from the host-side
            ! fills.
            ! Statistics first, where the data already is.  On a GPU build
            ! the buffer is device-resident, so one fused device reduction
            ! returns a handful of scalars instead of dragging the whole
            ! array through three host passes.  Without device residency
            ! (unit tests that skip `enter_data`) the host path is the only
            ! correct one -- a device read of an unmapped buffer under
            ! `-gpu=mem:separate` returns garbage without failing.  Both
            ! paths skip the NaN missing-data sentinel; see
            ! `diag_field_stats`.
            call diag_field_stats(v%output_buffer, this%on_device, &
                                  vmin, vmax, vmean, n_valid, n_total)

            ! The NetCDF write below needs the array itself on the host.
            !$acc update self(v%output_buffer) if_present
            v%fire_count = v%fire_count + 1
            write (line, "(A,F12.2,A,A,A,A,A,ES13.5,A,ES13.5,A,ES13.5)") &
               "[diag] t=", t_emit, " ", trim(v%name), " [", trim(v%units), &
               "]  min=", vmin, "  max=", vmax, "  mean=", vmean
            ! Only appended when cells were actually excluded, so a run with
            ! no masked cells emits the byte-identical line it always did.
            if (n_valid < n_total) then
               write (missing_tag, "(A,I0,A,I0)") "  missing=", n_total - n_valid, &
                  "/", n_total
               line = trim(line)//trim(missing_tag)
            end if
            call logger%info(trim(line))
            if (associated(this%emit_post_fire)) then
               call this%emit_post_fire(this, i, t_emit)
            end if
            v%dt_accum = 0.0_wp
            if (v%time_op /= DIAG_OP_INSTANT) then
               call reset_accumulator(v)
               v%n_accum = 0
            end if
         end associate
      end do
      this%dt_last_eval = t
   end subroutine ocean_diag_step

   pure subroutine diag_field_stats(buf, on_device, vmin, vmax, vmean, n_valid, n_total)
      !! The `[diag]` console line's min / max / mean for one diagnostic
      !! buffer, **over the finite cells only**.
      !!
      !! A diagnostic buffer legitimately carries IEEE NaN as its "no water
      !! here" sentinel: `fill_tracer_impl` writes one into every land
      !! column and every dynamically vanished layer (ZSTAR_FULL bed layers
      !! pinched out below `zstar_h_min`), because 0 degC / 0 PSU are legal
      !! ocean values and must not be confused with missing data.  Those
      !! cells are not data, so none of the three statistics may see them
      !! and the mean divides by `n_valid`, not by the array size.
      !!
      !! What this replaced, and why it was wrong: a plain
      !! `minval`/`maxval`/`sum` over the whole buffer.  Comparisons with
      !! NaN are FALSE, so `minval`/`maxval` silently skipped the sentinel
      !! cells while `sum` propagated them — emitting a self-contradictory
      !! `min= 1.5E+01  max= 1.5E+01  mean= NaN` for a run whose state was
      !! entirely healthy, and tripping the regression suite's NaN gate on
      !! every masked configuration (island / coastline / vanishing-layer
      !! cases).  Leaning on NaN-false comparisons is not portable either:
      !! nvfortran's relaxed-FP default may lower an unguarded `min`/`max`
      !! to a NaN-blind select (see CLAUDE.md's clamp-laundering gotcha).
      !!
      !! `n_valid == n_total` (the overwhelmingly common case — no land, no
      !! vanished layer) takes the unmasked intrinsics, so the emitted
      !! numbers stay BIT-IDENTICAL to the pre-fix behaviour on every
      !! all-finite field.  `n_valid == 0` reports `DIAG_MISSING_VALUE` for
      !! all three rather than the reduction's untouched `+huge`/`-huge`
      !! seeds or a zero that reads as a legal value.
      real(wp), intent(in) :: buf(:, :, :)  ! assumed-shape-ok: diag emit — cadence-bounded
      logical, intent(in) :: on_device
         !! `.true.` when `buf` is device-resident, so the fused device
         !! reduction is the correct (and only correct) reader — a host
         !! read of a mapped buffer under `-gpu=mem:separate` is stale.
      real(wp), intent(out) :: vmin, vmax, vmean
      integer, intent(out) :: n_valid, n_total
      real(wp) :: vsum

      n_total = size(buf)
      if (on_device) then
         call diag_reduce_stats(buf, size(buf, 1), size(buf, 2), size(buf, 3), &
                                vmin, vmax, vsum, n_valid)
      else
         n_valid = count(ieee_is_finite(buf))
         if (n_valid == n_total) then
            vmin = minval(buf)
            vmax = maxval(buf)
            vsum = sum(buf)
         else
            vmin = minval(buf, mask=ieee_is_finite(buf))
            vmax = maxval(buf, mask=ieee_is_finite(buf))
            vsum = sum(buf, mask=ieee_is_finite(buf))
         end if
      end if

      if (n_valid > 0) then
         vmean = vsum/real(n_valid, wp)
      else
         vmin = DIAG_MISSING_VALUE
         vmax = DIAG_MISSING_VALUE
         vmean = DIAG_MISSING_VALUE
      end if
   end subroutine diag_field_stats

   pure subroutine diag_reduce_stats(buf, n1, n2, n3, vmin, vmax, vsum, n_valid)
      !! Whole-array min / max / sum of a diagnostic buffer in ONE pass,
      !! over the FINITE cells only.
      !!
      !! Replaces three separate host passes (`minval`/`maxval`/`sum`) with
      !! a single `do concurrent ... reduce`, so on a GPU build this runs
      !! where the buffer already lives and only a few scalars come back.
      !! Explicit-shape dummies (never assumed-shape) so NVHPC does not walk
      !! a descriptor per launch; index order is `(k, j, i)` with the
      !! contiguous index last.
      !!
      !! **Missing data.**  A diagnostic buffer legitimately carries IEEE
      !! NaN as the "no water here" sentinel — `fill_tracer_impl` writes it
      !! into every land column and every dynamically vanished layer (see
      !! its docstring, and `test_fill_vanished_nan`).  Those cells are not
      !! data and must not enter the statistics, so every cell is
      !! `ieee_is_finite`-guarded and `n_valid` counts the cells that did
      !! contribute — the caller divides the sum by THAT, not by the array
      !! size.  Relying on "comparisons with NaN are false" to make
      !! `min`/`max` skip them is not enough and not portable: it leaves
      !! `sum` poisoned (the whole field's mean becomes NaN next to a
      !! perfectly finite min/max) and nvfortran's relaxed-FP default is
      !! free to lower an unguarded `min`/`max` to a NaN-blind select.
      !!
      !! For an all-finite buffer the guard is always taken, so the
      !! reduction order — and therefore the result — is bit-identical to
      !! the unguarded form.  `n_valid == 0` (nothing finite anywhere)
      !! leaves the `+huge` / `-huge` / `0` seeds untouched; the caller
      !! substitutes the missing-value sentinel.
      !!
      !! The caller MUST only invoke this when the buffer is device-resident
      !! on a GPU build — see `ocean_diag_t%on_device`.
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: buf(n1, n2, n3)
      real(wp), intent(out) :: vmin, vmax, vsum
      integer, intent(out), optional :: n_valid
         !! Number of finite cells folded in. Optional so the pre-existing
         !! three-scalar call sites keep working unchanged.
      integer :: i, j, k, nv
      vmin = huge(1.0_wp)
      vmax = -huge(1.0_wp)
      vsum = 0.0_wp
      nv = 0
      do concurrent(k=1:n3, j=1:n2, i=1:n1) reduce(min:vmin) reduce(max:vmax) &
         reduce(+:vsum) reduce(+:nv)
         if (ieee_is_finite(buf(i, j, k))) then
            vmin = min(vmin, buf(i, j, k))
            vmax = max(vmax, buf(i, j, k))
            vsum = vsum + buf(i, j, k)
            nv = nv + 1
         end if
      end do
      if (present(n_valid)) n_valid = nv
   end subroutine diag_reduce_stats

   subroutine fill_and_remap(v, state_handle, z_out, rho_out, sigma_out, zstar_out)
      !! Invoke the var's `fill` (and `remap` if non-LAYER) so
      !! `output_buffer` holds the current sample.  The accumulation
      !! / log emission step consumes `output_buffer` after this.
      !!
      !! The single target-levels array the `remap` proc receives is
      !! selected by the var's `output_vgrid`: Z_FIXED gets the z-level
      !! depths (`z_out`), DENSITY the isopycnal targets (`rho_out`), SIGMA
      !! the sigma fractions (`sigma_out`), ZSTAR the reference depths
      !! (`zstar_out`); other non-LAYER vgrids fall back to `z_out`.
      type(diag_var_t), intent(inout) :: v
      class(*), intent(in) :: state_handle
      ! Allocatable so an unconfigured (unallocated) target array is a legal
      ! actual argument — the select case below dereferences only the array
      ! matching the var's output_vgrid, which is allocated when that vgrid
      ! is in use.  (Non-allocatable dummies would make passing an
      ! unallocated `this%*_out` undefined behaviour, F2018 15.5.2.4.)
      real(wp), intent(in), allocatable :: z_out(:)
      real(wp), intent(in), allocatable :: rho_out(:)
      real(wp), intent(in), allocatable :: sigma_out(:)
      real(wp), intent(in), allocatable :: zstar_out(:)
      if (v%output_vgrid == DIAG_VGRID_LAYER) then
         call v%fill(state_handle, v%output_buffer)
      else
         call v%fill(state_handle, v%layer_buffer)
         if (associated(v%remap)) then
            select case (v%output_vgrid)
            case (DIAG_VGRID_DENSITY)
               call v%remap(state_handle, rho_out, v%layer_buffer, v%output_buffer, &
                            v%is_extensive)
            case (DIAG_VGRID_SIGMA)
               call v%remap(state_handle, sigma_out, v%layer_buffer, v%output_buffer, &
                            v%is_extensive)
            case (DIAG_VGRID_ZSTAR)
               call v%remap(state_handle, zstar_out, v%layer_buffer, v%output_buffer, &
                            v%is_extensive)
            case default
               call v%remap(state_handle, z_out, v%layer_buffer, v%output_buffer, &
                            v%is_extensive)
            end select
         else
            ! Identity passthrough when caller registered a non-LAYER
            ! output vgrid but no remap proc — device-side copy because
            ! both buffers are device-resident under the GPU-resident
            ! pipeline.
            call finalise_copy_impl(v%output_buffer, v%layer_buffer)
         end if
      end if
   end subroutine fill_and_remap

   subroutine fold_sample(v, dt)
      !! Combine the current `output_buffer` sample into `accumulator`
      !! per `time_op`.
      !!
      !! MEAN and INTEGRAL fold dt-weighted: `accumulator += sample · dt`.
      !! With fixed `dt` this matches a count-weighted sum exactly, so
      !! callers running uniform timesteps see no change.  With variable
      !! `dt` it produces the true `(1/T) ∫ f dt` (MEAN) or `∫ f dt`
      !! (INTEGRAL), which a count-weighted scheme would not.
      !!
      !! MAX / MIN are dt-independent — the running extremum doesn't
      !! care about sample weight.
      !!
      !! When `v%mask` is allocated, the sample is multiplied by the
      !! mask weight (broadcast across z) before folding.  Cells with
      !! weight = 0 contribute nothing to MEAN / INTEGRAL.  For MAX /
      !! MIN, masked cells are taken as `-huge` / `+huge` respectively
      !! so they never win — masked output stays at the seed value.
      type(diag_var_t), intent(inout) :: v
      real(wp), intent(in) :: dt
      if (.not. allocated(v%accumulator)) return
      if (allocated(v%mask)) then
         call fold_sample_masked_impl(v%accumulator, v%output_buffer, &
                                      v%mask%weight, v%mask%nx, v%mask%ny, &
                                      dt, v%time_op)
      else
         call fold_sample_unmasked_impl(v%accumulator, v%output_buffer, &
                                        dt, v%time_op)
      end if
      v%n_accum = v%n_accum + 1
   end subroutine fold_sample

   pure subroutine fold_sample_unmasked_impl(accum, out, dt, time_op)
      !! Whole-buffer fold without a region mask.  One `do concurrent`
      !! per op so the compiler can specialise — case-inside-loop blocks
      !! NVHPC device codegen.
      ! assumed-shape-ok: diag fold — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: accum(:, :, :)
      real(wp), intent(in)    :: out(:, :, :)  ! assumed-shape-ok: diag fold — cadence-bounded
      real(wp), intent(in)    :: dt
      integer, intent(in)    :: time_op
      integer :: i, j, k, nx, ny, nz
      nx = size(accum, 1)
      ny = size(accum, 2)
      nz = size(accum, 3)
      select case (time_op)
      case (DIAG_OP_MEAN, DIAG_OP_INTEGRAL)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            accum(i, j, k) = accum(i, j, k) + out(i, j, k)*dt
         end do
      case (DIAG_OP_MAX)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            accum(i, j, k) = max(accum(i, j, k), out(i, j, k))
         end do
      case (DIAG_OP_MIN)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            accum(i, j, k) = min(accum(i, j, k), out(i, j, k))
         end do
      case default
         ! INSTANT and unknown ops don't accumulate.  The caller in
         ! `fold_sample` already returns early when `accumulator` is
         ! unallocated (which is the INSTANT path), so reaching here
         ! with an unknown op is a defensive no-op.
      end select
   end subroutine fold_sample_unmasked_impl

   pure subroutine fold_sample_masked_impl(accum, out, weight, mask_nx, mask_ny, &
                                           dt, time_op)
      !! Masked fold — multiplies sample by `weight(i, j)` before folding.
      !! MAX/MIN treat `weight == 0` as "don't update" so masked-out cells
      !! keep their seed value.  Loop bounds clip to whichever extent is
      !! smaller (accumulator vs mask) so an undersized mask doesn't OOB.
      ! assumed-shape-ok: diag fold — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: accum(:, :, :)
      real(wp), intent(in)    :: out(:, :, :)  ! assumed-shape-ok: diag fold — cadence-bounded
      real(wp), intent(in)    :: weight(:, :)  ! assumed-shape-ok: diag fold — cadence-bounded
      integer, intent(in)    :: mask_nx, mask_ny
      real(wp), intent(in)    :: dt
      integer, intent(in)    :: time_op
      integer :: i, j, k, nx, ny, nz
      real(wp) :: w
      nx = min(size(accum, 1), mask_nx)
      ny = min(size(accum, 2), mask_ny)
      nz = size(accum, 3)
      select case (time_op)
      case (DIAG_OP_MEAN, DIAG_OP_INTEGRAL)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            accum(i, j, k) = accum(i, j, k) + &
                             weight(i, j)*out(i, j, k)*dt
         end do
      case (DIAG_OP_MAX)
         do concurrent(k=1:nz, j=1:ny, i=1:nx) &
            local(w)
            w = weight(i, j)
            if (w > 0.0_wp) then
               accum(i, j, k) = max(accum(i, j, k), out(i, j, k))
            end if
         end do
      case (DIAG_OP_MIN)
         do concurrent(k=1:nz, j=1:ny, i=1:nx) &
            local(w)
            w = weight(i, j)
            if (w > 0.0_wp) then
               accum(i, j, k) = min(accum(i, j, k), out(i, j, k))
            end if
         end do
      case default
         ! INSTANT and unknown ops don't accumulate (caller-gated).
      end select
   end subroutine fold_sample_masked_impl

   subroutine finalise_accumulator(v)
      !! Copy the accumulator into `output_buffer`.  Called at
      !! cadence-fire BEFORE log / NetCDF emission consumes
      !! `output_buffer`.
      !!
      !! MEAN     -> `accumulator / dt_accum`  (dt-weighted time mean)
      !! INTEGRAL -> `accumulator`             (raw cumulative integral)
      !! MAX/MIN  -> `accumulator`             (running extremum)
      type(diag_var_t), intent(inout) :: v
      real(wp), parameter :: TINY_DT = 1.0e-12_wp
      if (.not. allocated(v%accumulator)) return
      if (v%n_accum == 0) return
      select case (v%time_op)
      case (DIAG_OP_MEAN)
         if (v%dt_accum < TINY_DT) return
         call finalise_scale_impl(v%output_buffer, v%accumulator, &
                                  1.0_wp/v%dt_accum)
      case (DIAG_OP_INTEGRAL, DIAG_OP_MAX, DIAG_OP_MIN)
         call finalise_copy_impl(v%output_buffer, v%accumulator)
      case default
         return
      end select
   end subroutine finalise_accumulator

   pure subroutine finalise_scale_impl(out, accum, scale_factor)
      !! Device-side `out = accum * scale_factor`.  Used for MEAN
      !! finalise where `scale_factor = 1/dt_accum`.
      !!
      !! Local var named `scale_factor` (not `scale`) to dodge the
      !! NVHPC 26.3 intrinsic-shadow bug
      !! (`feedback_nvhpc_local_intrinsic_shadow.md`).
      ! assumed-shape-ok: diag finalise — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: out(:, :, :)
      real(wp), intent(in)    :: accum(:, :, :)  ! assumed-shape-ok: diag finalise — cadence-bounded
      real(wp), intent(in)    :: scale_factor
      integer :: i, j, k, nx, ny, nz
      nx = size(out, 1)
      ny = size(out, 2)
      nz = size(out, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         out(i, j, k) = accum(i, j, k)*scale_factor
      end do
   end subroutine finalise_scale_impl

   pure subroutine finalise_copy_impl(out, accum)
      ! assumed-shape-ok: diag finalise — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: out(:, :, :)
      real(wp), intent(in)    :: accum(:, :, :)  ! assumed-shape-ok: diag finalise — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = size(out, 1)
      ny = size(out, 2)
      nz = size(out, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         out(i, j, k) = accum(i, j, k)
      end do
   end subroutine finalise_copy_impl

   pure function diag_var_bytes(this) result(nbytes)
      !! Counted allocatable footprint of ONE registered diagnostic
      !! (0 for every buffer that is unallocated).
      !!
      !! All four terms are device-mapped by `ocean_diag_enter_data_impl`
      !! and, until this function existed, none of them appeared in any
      !! `bytes()` total — the single largest hole in the startup estimate
      !! (~3.2 GB for the default catalog at 1000x800x50, up to ~9.9 GB
      !! with the optional diags).  Conditionality rides `arr_bytes`:
      !! `layer_buffer` is allocated only when `output_vgrid /= LAYER`,
      !! `accumulator` only when `time_op /= INSTANT`, `mask` only when a
      !! region was supplied at register time.
      class(diag_var_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%output_buffer) &
               + arr_bytes(this%layer_buffer) &
               + arr_bytes(this%accumulator)
      if (allocated(this%mask)) nbytes = nbytes + this%mask%bytes()
   end function diag_var_bytes

   pure function ocean_diag_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the diagnostics slot: the flat
      !! remap-level arrays, the host NetCDF staging buffer, AND the
      !! per-variable registry buffers (0 when unallocated).  One
      !! arr_bytes term per array — add a term here when a new allocatable
      !! joins the type; new `diag_var_t` buffers go in `diag_var_bytes`.
      class(ocean_diag_t), intent(in) :: this
      integer(int64) :: nbytes
      integer :: i
      nbytes = arr_bytes(this%z_out) &
               + arr_bytes(this%rho_out) &
               + arr_bytes(this%sigma_out) &
               + arr_bytes(this%zstar_out) &
               + arr_bytes(this%send_buf)
      ! `nc_stream%stage` is real32, which `arr_bytes` (real(wp) only) does
      ! not cover — counted inline at 4 bytes/element.  Unallocated on the
      ! default double-precision diag stream.
      if (allocated(this%nc_stream%stage)) then
         nbytes = nbytes + 4_int64*int(size(this%nc_stream%stage), int64)
      end if
      ! Registry: `nvars` live entries out of `nvars_max` capacity — the
      ! tail slots past `nvars` are default-initialised and unallocated,
      ! but bound the loop by `nvars` anyway so the term tracks what is
      ! actually registered (and mapped).
      if (allocated(this%vars)) then
         do i = 1, min(this%nvars, size(this%vars))
            nbytes = nbytes + this%vars(i)%bytes()
         end do
      end if
   end function ocean_diag_bytes

end module rdb_ocean_diag
