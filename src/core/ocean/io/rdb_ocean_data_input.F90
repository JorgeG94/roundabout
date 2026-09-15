!! The shared time-varying NetCDF input reader (PR-14).
module rdb_ocean_data_input
   !! Generic, decomposition-aware, device-resident reader for
   !! time-varying NetCDF fields `field(x, y[, z], t)`.  Owns the file
   !! handles, the time axis, the bracket bookkeeping, the H->D motion
   !! (only when the bracket advances), and the per-step linear-in-time
   !! blend kernel.  Owns NO field names and NO state-slot mapping: a
   !! consumer registers its own `(file, variable, destination)` triple
   !! via `ocean_data_input_register_2d/_3d` (or the edge-segment
   !! variants for OBC files) and gets back an opaque `id`.
   !!
   !! **The binding seam (see `PLAN_PR14_netcdf_input_reader.md` §13):**
   !!
   !!   * The reader never touches state.  It blends on-device into the
   !!     caller's own, caller-MAPPED, explicit-shape, WHOLE (never a
   !!     section) array, at a registration-time offset
   !!     `(dest_i0, dest_j0)`.  Passing a section into an explicit-shape
   !!     dummy would trigger a host copy-in/out and silently detach
   !!     from the device map — never do that.
   !!   * The consumer owns the halo exchange (physical cells only are
   !!     filled) and the vertical (k) orientation — a source file's z
   !!     axis is read AS STORED; nothing here flips it.  A 3-D field
   !!     read from a surface-first source file is handed back
   !!     bottom-last; if that's wrong for your consumer, flip it there.
   !!   * `&ocean_data_nml` carries no per-field entries — only
   !!     `max_fields`/`verbose`.  Each consumer defines its own file/
   !!     variable keys in its own namelist group.
   !!   * Out-of-range time is an ABORT by default (`DATA_OOR_ERROR`);
   !!     `DATA_OOR_CLAMP` is opt-in per field.
   !!   * No in-core horizontal interpolation (files are pre-regridded
   !!     to the model's global horizontal extent), no vertical
   !!     remapping of a source z axis, no calendar.  Time interpolation
   !!     is linear between the two bracketing records; `DATA_TIME_CYCLIC`
   !!     wraps a climatology through an explicit period; `DATA_TIME_STATIC`
   !!     reads record 1 once and never touches the file again.
   !!   * Time interpolation is a PURE function of `t` — no reader-side
   !!     restart state.  A resumed run reproduces every field
   !!     bit-for-bit from `t_current` alone.
   !!   * `ocean_data_input_update_all(this, t)` refreshes every
   !!     registered field's bracket (+ pushes any new slab to the
   !!     device) and is the driver's one-line per-step hook — a no-op
   !!     when `nfields == 0` (every shipped namelist today).  It does
   !!     NOT write into any consumer array (the reader stores the
   !!     registration-time `dest` shape as metadata only, never a
   !!     pointer).  Each consumer then calls `update_2d/_3d` with its
   !!     OWN array at its own point in the step, and that call is
   !!     checked against the `t` `update_all` most recently refreshed
   !!     with — calling `update_2d/_3d` before `update_all` has run
   !!     this step's `t` is a fail-loud ordering bug, not a silent
   !!     stale read.
   !!   * `ocean_data_input_fill_static_host` is a one-shot HOST-side
   !!     fill for `DATA_TIME_STATIC` fields, legal at SETUP time before
   !!     the destination is device-mapped (PR-24 ice IC, PR-30 tidal
   !!     maps).  Touches no device memory; fails loud on a non-static
   !!     field.
   !!
   !! Storage order: files must be FORTRAN-ordered `(x, y[, z], t)` —
   !! v1 detects a horizontal-dimension mismatch and aborts rather than
   !! silently transposing (a transposed field is exactly the class of
   !! bug the analytical tests are built to catch); regrid/reorder the
   !! file if this trips.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: ocean_data_config_t
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, nc_get_var_1d, &
                            nc_get_var_slab_3d, nc_get_att_text
   use netcdf, only: nf90_noerr, nf90_inq_varid, nf90_inquire_variable, &
                     nf90_inquire_dimension
   use rdb_mem_report, only: arr_bytes
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_ascii, only: to_lower
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP, OCEAN_STATUS_ERR_IO
   use rdb_error_ring, only: fail
   implicit none
   private

   public :: ocean_data_input_t
   public :: data_input_field_t

   ! --- time modes ---
   public :: DATA_TIME_LINEAR, DATA_TIME_CYCLIC, DATA_TIME_STATIC
   ! --- out-of-range policy ---
   public :: DATA_OOR_ERROR, DATA_OOR_CLAMP
   ! --- edge tags (segment registration, PR-22) ---
   public :: DATA_EDGE_WEST, DATA_EDGE_EAST, DATA_EDGE_SOUTH, DATA_EDGE_NORTH

   public :: ocean_data_input_register_2d
   public :: ocean_data_input_register_3d
   public :: ocean_data_input_register_segment_2d
   public :: ocean_data_input_register_segment_3d
   public :: ocean_data_input_update_2d
   public :: ocean_data_input_update_3d
   public :: ocean_data_input_fill_static_host
   public :: ocean_data_input_fill_static_host_3d
   public :: ocean_data_input_update_all

   public :: data_input_locate
   public :: data_input_time_scale_from_units
   public :: data_input_time_mode_is_implemented
   public :: data_input_time_mode_from_string
   public :: data_input_dims_ok

   integer, parameter :: DATA_TIME_LINEAR = 1
   integer, parameter :: DATA_TIME_CYCLIC = 2
   integer, parameter :: DATA_TIME_STATIC = 3

   integer, parameter :: DATA_OOR_ERROR = 1
   integer, parameter :: DATA_OOR_CLAMP = 2

   integer, parameter :: DATA_EDGE_WEST = 1
   integer, parameter :: DATA_EDGE_EAST = 2
   integer, parameter :: DATA_EDGE_SOUTH = 3
   integer, parameter :: DATA_EDGE_NORTH = 4

   type :: data_input_field_t
      !! One registered `(file, variable, destination-shape)` triple.
      !! Array-of-derived-type element — NEVER dereference this from
      !! inside a `do concurrent`/device-kernel body (outer-shim +
      !! flat-impl rule); every place this type is touched below is a
      !! host-side registry walk that hands plain explicit-shape arrays
      !! to a `pure` flat-impl kernel.
      logical :: active = .false.
      integer :: ncid = -1
      integer :: varid = -1
      logical :: is_3d = .false.
      integer :: time_mode = DATA_TIME_LINEAR
      integer :: oor = DATA_OOR_ERROR
      logical :: oor_warned = .false.

      ! Source-slab extents (this rank's own subdomain slab, source-k
      ! order, source-z count for a 3-D field, else 1).
      integer :: nx = 0, ny = 0, nz = 1
      ! File-side 1-based start offsets for THIS rank's horizontal slab.
      integer :: i0 = 1, j0 = 1
      ! Where the slab lands in the consumer's `dest` array.
      integer :: dest_i0 = 1, dest_j0 = 1
      ! Shape metadata of the consumer's `dest` — validated, not stored
      ! as a pointer (see module docstring).
      integer :: dest_n1 = 0, dest_n2 = 0, dest_n3 = 1

      real(wp) :: scale = 1.0_wp
      real(wp) :: add_offset = 0.0_wp
      real(wp) :: t_offset = 0.0_wp
      real(wp) :: cycle_period = 0.0_wp

      integer :: nt = 0
      real(wp), allocatable :: t_axis(:)
         !! File time axis, converted to seconds (scale + t_offset NOT
         !! applied here — t_offset is applied to the QUERY time, not
         !! the axis; see `data_input_refresh_brackets`).

      real(wp), allocatable :: f0(:, :, :), f1(:, :, :)
         !! Bracketing records, shape `(nx, ny, nz)`.  Device-mapped by
         !! `ocean_data_input_t%enter_data`.
      integer :: rec0 = -1, rec1 = -1
         !! Currently-loaded record indices (1-based file record
         !! numbers); -1 = nothing loaded yet.
      real(wp) :: w = 0.0_wp
         !! Current blend weight: `f = (1-w)*f0 + w*f1`.
      real(wp) :: t_last = -huge(1.0_wp)
         !! Model time `update_all` most recently refreshed this field
         !! at.  `update_2d/_3d` asserts the caller's `t` matches (fail
         !! loud on an ordering bug — see module docstring).  Ignored
         !! for `DATA_TIME_STATIC` (never refreshed after registration).
      integer :: nreads = 0
         !! Slab-read counter (test hook — T4 asserts exactly one read
         !! for a STATIC field's whole run).
   end type data_input_field_t

   type :: ocean_data_input_t
      logical :: is_init = .false.
      logical :: verbose = .false.
      integer :: nfields = 0
      integer :: nfields_max = 0
      type(data_input_field_t), allocatable :: fields(:)
   contains
      procedure, non_overridable :: init => ocean_data_input_init
      procedure, non_overridable :: destroy => ocean_data_input_destroy
      procedure, non_overridable :: enter_data => ocean_data_input_enter_data
      procedure, non_overridable :: exit_data => ocean_data_input_exit_data
      procedure, non_overridable :: bytes => ocean_data_input_bytes
   end type ocean_data_input_t

   ! Module-level host read-scratch (registration-time-only allocation
   ! per CLAUDE.md rule 5 — no per-step allocate on `-stdpar=gpu`; this
   ! buffer never touches the device, so it carries no `!$acc` directives
   ! at all).  Copy of the `remap_workspace_ensure` idiom
   ! (`rdb_ml_dynamics.F90`).
   real(wp), allocatable :: data_input_ws(:, :, :)

contains

   ! ======================================================================
   ! Pure, independently-testable helpers
   ! ======================================================================

   pure subroutine data_input_locate(t_axis, nt, mode, cycle_period, t, n0, n1, w, out_of_range)
      !! Bracket search + blend weight for a query time `t` (already in
      !! file-time units/offset — the caller applies `t_offset`/`t_scale`
      !! before calling).  `nt == 1` is degenerate: always returns
      !! `n0 = n1 = 1`, `w = 0`, never out of range.
      integer, intent(in) :: nt, mode
      real(wp), intent(in) :: t_axis(nt)
      real(wp), intent(in) :: cycle_period, t
      integer, intent(out) :: n0, n1
      real(wp), intent(out) :: w
      logical, intent(out) :: out_of_range

      integer :: k
      real(wp) :: t_eff, gap

      out_of_range = .false.

      if (nt <= 1) then
         n0 = 1
         n1 = 1
         w = 0.0_wp
         return
      end if

      if (mode == DATA_TIME_CYCLIC) then
         t_eff = t_axis(1) + modulo(t - t_axis(1), cycle_period)
         if (t_eff >= t_axis(nt)) then
            ! The seam bracket (nt, 1). Gap is NOT the mean record
            ! spacing — it is the wrap distance from the last record to
            ! the first record of the NEXT cycle.
            n0 = nt
            n1 = 1
            gap = (t_axis(1) + cycle_period) - t_axis(nt)
            if (gap > 0.0_wp) then
               w = (t_eff - t_axis(nt))/gap
            else
               w = 0.0_wp
            end if
         else
            n0 = 1
            do k = 1, nt - 1
               if (t_eff >= t_axis(k) .and. t_eff < t_axis(k + 1)) then
                  n0 = k
                  exit
               end if
            end do
            n1 = n0 + 1
            w = (t_eff - t_axis(n0))/(t_axis(n1) - t_axis(n0))
         end if
      else
         ! DATA_TIME_LINEAR (STATIC never reaches here — the caller
         ! returns before calling this for a static field).
         if (t < t_axis(1)) then
            n0 = 1
            n1 = 1
            w = 0.0_wp
            out_of_range = .true.
         else if (t > t_axis(nt)) then
            n0 = nt
            n1 = nt
            w = 0.0_wp
            out_of_range = .true.
         else if (t >= t_axis(nt)) then
            n0 = nt
            n1 = nt
            w = 0.0_wp
         else
            n0 = 1
            do k = 1, nt - 1
               if (t >= t_axis(k) .and. t < t_axis(k + 1)) then
                  n0 = k
                  exit
               end if
            end do
            n1 = n0 + 1
            w = (t - t_axis(n0))/(t_axis(n1) - t_axis(n0))
         end if
      end if
   end subroutine data_input_locate

   pure subroutine data_input_time_scale_from_units(units, scale, ok)
      !! CF `units` attribute ("seconds since ...", "hours since ...", ...)
      !! -> a multiplier converting the raw file time axis to seconds.
      !! Only the leading unit word matters (no calendar, no reference
      !! date — see module docstring).  `ok = .false.` for an
      !! unrecognised leading word; caller decides the fallback.
      character(len=*), intent(in) :: units
      real(wp), intent(out) :: scale
      logical, intent(out) :: ok

      character(len=:), allocatable :: u

      u = to_lower(adjustl(units))
      ok = .true.
      if (index(u, "second") == 1 .or. index(u, "sec") == 1) then
         scale = 1.0_wp
      else if (index(u, "minute") == 1 .or. index(u, "min") == 1) then
         scale = 60.0_wp
      else if (index(u, "hour") == 1) then
         scale = 3600.0_wp
      else if (index(u, "day") == 1) then
         scale = 86400.0_wp
      else
         scale = 1.0_wp
         ok = .false.
      end if
   end subroutine data_input_time_scale_from_units

   pure logical function data_input_time_mode_is_implemented(tag) result(ok)
      !! `.true.` iff `tag` (case-insensitive) names a shipped time mode.
      !! Drives the fail-loud dispatch in `data_input_time_mode_from_string`
      !! — house idiom, see `lateral_closure_is_implemented`.
      character(len=*), intent(in) :: tag
      select case (to_lower(trim(tag)))
      case ("linear", "cyclic", "static")
         ok = .true.
      case default
         ok = .false.
      end select
   end function data_input_time_mode_is_implemented

   function data_input_time_mode_from_string(tag) result(mode)
      !! Translate a namelist/registration-time string into a
      !! `DATA_TIME_*` code.  Fail-loud on anything not covered by
      !! `data_input_time_mode_is_implemented` — an unrecognised tag
      !! must never silently fall back to a default mode.
      character(len=*), intent(in) :: tag
      integer :: mode
      select case (to_lower(trim(tag)))
      case ("linear")
         mode = DATA_TIME_LINEAR
      case ("cyclic")
         mode = DATA_TIME_CYCLIC
      case ("static")
         mode = DATA_TIME_STATIC
      case default
         call logger%error("ocean_data_input: unimplemented time_mode '"//trim(tag)// &
                           "' (implemented: linear, cyclic, static)")
         error stop "ocean_data_input: unimplemented time_mode"
      end select
   end function data_input_time_mode_from_string

   pure logical function dims_geometry_ok(nd, expect_nd, d_horiz1, d_horiz2, &
                                          expect1, expect2) result(ok)
      !! Low-level geometry check used inline by `register_common`:
      !! `.true.` iff the variable has the expected rank and its two
      !! validated-length dims are each at least as long as required
      !! (decomposition-safe: a subdomain slab only needs `offset+extent`
      !! to fit, not to equal the file's global length).  Never aborts.
      integer, intent(in) :: nd, expect_nd, d_horiz1, d_horiz2, expect1, expect2
      ok = (nd == expect_nd) .and. (d_horiz1 >= expect1) .and. (d_horiz2 >= expect2)
   end function dims_geometry_ok

   function data_input_dims_ok(filename, var, nx_phys, ny_phys, is_3d, nz_src) result(ok)
      !! Self-contained, non-erroring dimension validator — the
      !! single-rank (`i_offset_global = j_offset_global = 0`) testable
      !! twin of `register_common`'s inline checks; mirrors
      !! `zinit_dims_ok`.  Opens `filename`, inspects `var`'s rank and
      !! lengths, and returns `.false.` on any mismatch, missing
      !! variable, or I/O error — NEVER aborts, so test code can probe
      !! the false branch without a subprocess/death-test harness (the
      !! house convention — see `zinit_dims_ok`).
      use netcdf, only: nf90_open, nf90_close, nf90_nowrite
      character(len=*), intent(in) :: filename, var
      integer, intent(in) :: nx_phys, ny_phys
      logical, intent(in) :: is_3d
      integer, intent(in), optional :: nz_src
      logical :: ok

      integer :: ncid, varid, ierr, var_ndims, expect_nd
      integer :: var_dimids(4)
      integer :: d1, d2, d3

      ok = .false.
      ierr = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (ierr /= nf90_noerr) return

      expect_nd = merge(4, 3, is_3d)
      ierr = nf90_inq_varid(ncid, trim(var), varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_inquire_variable(ncid, varid, ndims=var_ndims, &
                                      dimids=var_dimids(1:expect_nd))
      end if
      if (ierr == nf90_noerr .and. var_ndims /= expect_nd) ierr = -1
      if (ierr == nf90_noerr) ierr = nf90_inquire_dimension(ncid, var_dimids(1), len=d1)
      if (ierr == nf90_noerr) ierr = nf90_inquire_dimension(ncid, var_dimids(2), len=d2)
      if (is_3d .and. ierr == nf90_noerr) then
         ierr = nf90_inquire_dimension(ncid, var_dimids(3), len=d3)
      end if

      if (ierr == nf90_noerr) then
         ok = (d1 >= nx_phys) .and. (d2 >= ny_phys)
         if (is_3d .and. present(nz_src)) ok = ok .and. (d3 == nz_src)
      end if

      ierr = nf90_close(ncid)
   end function data_input_dims_ok

   ! ======================================================================
   ! Blend kernels — pure, explicit-shape, integer dims declared first
   ! (decl-order), a single `do concurrent` each.  Called ONLY from the
   ! outer-shim `update_2d/_3d` below with plain arrays pulled out of the
   ! registry — never pass `this%fields(id)%f0` etc. into a kernel body.
   ! ======================================================================

   pure subroutine data_input_blend_2d_impl(nxs, nys, dn1, dn2, i0, j0, f0, f1, w, dest)
      integer, intent(in) :: nxs, nys, dn1, dn2, i0, j0
      real(wp), intent(in) :: f0(nxs, nys), f1(nxs, nys)
      real(wp), intent(in) :: w
      real(wp), intent(inout) :: dest(dn1, dn2)
      integer :: i, j
      do concurrent(j=1:nys, i=1:nxs)
         dest(i0 + i - 1, j0 + j - 1) = (1.0_wp - w)*f0(i, j) + w*f1(i, j)
      end do
   end subroutine data_input_blend_2d_impl

   pure subroutine data_input_blend_3d_impl(nxs, nys, nzs, dn1, dn2, dn3, i0, j0, f0, f1, w, dest)
      integer, intent(in) :: nxs, nys, nzs, dn1, dn2, dn3, i0, j0
      real(wp), intent(in) :: f0(nxs, nys, nzs), f1(nxs, nys, nzs)
      real(wp), intent(in) :: w
      real(wp), intent(inout) :: dest(dn1, dn2, dn3)
      integer :: i, j, k
      do concurrent(k=1:nzs, j=1:nys, i=1:nxs)
         dest(i0 + i - 1, j0 + j - 1, k) = (1.0_wp - w)*f0(i, j, k) + w*f1(i, j, k)
      end do
   end subroutine data_input_blend_3d_impl

   ! ======================================================================
   ! Lifecycle
   ! ======================================================================

   subroutine ocean_data_input_init(this, cfg)
      !! Allocate the field registry.  `cfg` optional so `ocean_state_init`
      !! (which many tests call directly with no `config_t` in hand) can
      !! construct a default-sized (16 slots, quiet) reader identically to
      !! the scaffold it replaces; `ocean_state_init_from_config` re-calls
      !! this with the real `cfg%ocean%data` once `cfg` is available —
      !! safe because nothing is registered between the two calls.
      class(ocean_data_input_t), intent(inout) :: this
      type(ocean_data_config_t), intent(in), optional :: cfg
      integer :: max_fields
      logical :: verbose

      max_fields = 16
      verbose = .false.
      if (present(cfg)) then
         max_fields = cfg%max_fields
         verbose = cfg%verbose
      end if

      call this%destroy()
      allocate (this%fields(max(max_fields, 1)))
      this%nfields = 0
      this%nfields_max = max_fields
      this%verbose = verbose
      this%is_init = .true.
   end subroutine ocean_data_input_init

   subroutine ocean_data_input_destroy(this)
      class(ocean_data_input_t), intent(inout) :: this
      integer :: i
      if (allocated(this%fields)) then
         do i = 1, this%nfields
            if (this%fields(i)%ncid >= 0) call nc_close(this%fields(i)%ncid)
         end do
         deallocate (this%fields)
      end if
      this%nfields = 0
      this%nfields_max = 0
      this%is_init = .false.
      call data_input_workspace_cleanup()
   end subroutine ocean_data_input_destroy

   subroutine ocean_data_input_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl (the
      !! AMD libomptarget cross-slot-overlap fix; see
      !! `rdb_ocean_surface_stress.F90`).
      class(ocean_data_input_t), intent(inout) :: this
      select type (this)
      type is (ocean_data_input_t)
         call ocean_data_input_enter_data_impl(this)
      end select
   end subroutine ocean_data_input_enter_data

   subroutine ocean_data_input_enter_data_impl(this)
      type(ocean_data_input_t), intent(inout) :: this
      integer :: i
      if (.not. this%is_init) return
      do i = 1, this%nfields
         if (.not. allocated(this%fields(i)%f0)) cycle
         !$acc enter data copyin(this%fields(i)%f0, this%fields(i)%f1)
         ! Rule (2): create-mapped arrays do NOT carry the host values
         ! read at registration time — push them explicitly.
         !$acc update device(this%fields(i)%f0, this%fields(i)%f1)
      end do
   end subroutine ocean_data_input_enter_data_impl

   subroutine ocean_data_input_exit_data(this)
      class(ocean_data_input_t), intent(inout) :: this
      select type (this)
      type is (ocean_data_input_t)
         call ocean_data_input_exit_data_impl(this)
      end select
   end subroutine ocean_data_input_exit_data

   subroutine ocean_data_input_exit_data_impl(this)
      type(ocean_data_input_t), intent(inout) :: this
      integer :: i
      if (.not. this%is_init) return
      do i = 1, this%nfields
         if (.not. allocated(this%fields(i)%f0)) cycle
         !$acc exit data delete(this%fields(i)%f0, this%fields(i)%f1)
      end do
   end subroutine ocean_data_input_exit_data_impl

   pure function ocean_data_input_bytes(this) result(nbytes)
      !! Counted allocatable footprint (0 when unallocated).  Summed over
      !! every registered field's `f0`/`f1` — the `t_axis` is small
      !! (`nt` reals) and intentionally excluded, matching the house
      !! convention of counting device-resident footprint only.
      class(ocean_data_input_t), intent(in) :: this
      integer(int64) :: nbytes
      integer :: i
      nbytes = 0_int64
      if (.not. allocated(this%fields)) return
      do i = 1, this%nfields
         nbytes = nbytes + arr_bytes(this%fields(i)%f0) + arr_bytes(this%fields(i)%f1)
      end do
   end function ocean_data_input_bytes

   ! ======================================================================
   ! Registration
   ! ======================================================================

   subroutine ocean_data_input_register_2d(this, file, var, grid, dest_n1, dest_n2, &
                                           dest_i0, dest_j0, id, time_mode, cycle_period, &
                                           t_offset, t_scale, scale, add_offset, oor, &
                                           nx_extra, ny_extra, ierr)
      !! Register a 2-D time-varying field `f(x, y, t)`.  Opens the file
      !! now and keeps the handle for the run.  See the module docstring
      !! for the destination-offset contract.
      !!
      !! `nx_extra`/`ny_extra` (default 0) widen this rank's slab by that
      !! many values in x/y.  They exist for C-grid FACE fields: an
      !! x-face array spans `nx_phys + 1` faces, not `nx_phys`, and the
      !! trailing face is OWNED by this rank (the west/south-owns-the-seam
      !! rule), so no exchange can supply it — it has to come from the
      !! file.  A cell-centred field leaves both at 0.
      class(ocean_data_input_t), intent(inout) :: this
      character(len=*), intent(in) :: file, var
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: dest_n1, dest_n2, dest_i0, dest_j0
      integer, intent(out) :: id
      character(len=*), intent(in), optional :: time_mode
      real(wp), intent(in), optional :: cycle_period, t_offset, t_scale, scale, add_offset
      integer, intent(in), optional :: oor
      integer, intent(in), optional :: nx_extra, ny_extra
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`/`OCEAN_STATUS_ERR_IO`) on a
         !! registry/file/dimension failure when present; absent behaves
         !! as today (`error stop`).

      integer :: nxe, nye

      nxe = 0
      if (present(nx_extra)) nxe = nx_extra
      nye = 0
      if (present(ny_extra)) nye = ny_extra

      call register_common(this, file, var, .false., &
                           grid%i_offset_global + 1, grid%j_offset_global + 1, &
                           grid%nx_phys + nxe, grid%ny_phys + nye, 1, &
                           dest_i0, dest_j0, dest_n1, dest_n2, 1, &
                           time_mode, cycle_period, t_offset, t_scale, scale, add_offset, oor, id, &
                           ierr=ierr)
   end subroutine ocean_data_input_register_2d

   subroutine ocean_data_input_register_3d(this, file, var, grid, nz_src, dest_n1, dest_n2, &
                                           dest_n3, dest_i0, dest_j0, id, time_mode, &
                                           cycle_period, t_offset, t_scale, scale, add_offset, oor, &
                                           ierr)
      !! Register a 3-D time-varying field `f(x, y, z, t)`.  `nz_src` is
      !! the source z-level count (the consumer's own concern — the
      !! reader validates it against the file's z dim and does NOT flip
      !! or remap k; see the module docstring).
      class(ocean_data_input_t), intent(inout) :: this
      character(len=*), intent(in) :: file, var
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_src, dest_n1, dest_n2, dest_n3, dest_i0, dest_j0
      integer, intent(out) :: id
      character(len=*), intent(in), optional :: time_mode
      real(wp), intent(in), optional :: cycle_period, t_offset, t_scale, scale, add_offset
      integer, intent(in), optional :: oor
      integer, intent(out), optional :: ierr

      call register_common(this, file, var, .true., &
                           grid%i_offset_global + 1, grid%j_offset_global + 1, &
                           grid%nx_phys, grid%ny_phys, nz_src, &
                           dest_i0, dest_j0, dest_n1, dest_n2, dest_n3, &
                           time_mode, cycle_period, t_offset, t_scale, scale, add_offset, oor, id, &
                           ierr=ierr)
   end subroutine ocean_data_input_register_3d

   subroutine ocean_data_input_register_segment_2d(this, file, var, grid, edge, &
                                                   dest_n1, dest_n2, id, time_mode, &
                                                   cycle_period, t_offset, t_scale, &
                                                   scale, add_offset, oor, ierr)
      !! Register a 2-D OBC-segment field.  The degenerate horizontal
      !! axis (x for west/east, y for south/north) reads as `start=1,
      !! count=1`; the along-edge axis slices from the global index
      !! range exactly as `register_2d` does.  `dest_i0`/`dest_j0` are
      !! implied by `edge` (degenerate axis -> index 1; along-edge axis
      !! -> `grid%nghost + 1`) — not arguments (see module docstring).
      class(ocean_data_input_t), intent(inout) :: this
      character(len=*), intent(in) :: file, var
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: edge, dest_n1, dest_n2
      integer, intent(out) :: id
      character(len=*), intent(in), optional :: time_mode
      real(wp), intent(in), optional :: cycle_period, t_offset, t_scale, scale, add_offset
      integer, intent(in), optional :: oor
      integer, intent(out), optional :: ierr

      integer :: i0, j0, nx, ny, di0, dj0

      call segment_geometry(grid, edge, i0, j0, nx, ny, di0, dj0, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_common(this, file, var, .false., i0, j0, nx, ny, 1, &
                           di0, dj0, dest_n1, dest_n2, 1, &
                           time_mode, cycle_period, t_offset, t_scale, scale, add_offset, oor, id, &
                           edge=edge, ierr=ierr)
   end subroutine ocean_data_input_register_segment_2d

   subroutine ocean_data_input_register_segment_3d(this, file, var, grid, edge, nz_src, &
                                                   dest_n1, dest_n2, dest_n3, id, time_mode, &
                                                   cycle_period, t_offset, t_scale, &
                                                   scale, add_offset, oor, ierr)
      !! 3-D twin of `register_segment_2d` (PR-22 OBC segment files).
      class(ocean_data_input_t), intent(inout) :: this
      character(len=*), intent(in) :: file, var
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: edge, nz_src, dest_n1, dest_n2, dest_n3
      integer, intent(out) :: id
      character(len=*), intent(in), optional :: time_mode
      real(wp), intent(in), optional :: cycle_period, t_offset, t_scale, scale, add_offset
      integer, intent(in), optional :: oor
      integer, intent(out), optional :: ierr

      integer :: i0, j0, nx, ny, di0, dj0

      call segment_geometry(grid, edge, i0, j0, nx, ny, di0, dj0, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_common(this, file, var, .true., i0, j0, nx, ny, nz_src, &
                           di0, dj0, dest_n1, dest_n2, dest_n3, &
                           time_mode, cycle_period, t_offset, t_scale, scale, add_offset, oor, id, &
                           edge=edge, ierr=ierr)
   end subroutine ocean_data_input_register_segment_3d

   subroutine segment_geometry(grid, edge, i0, j0, nx, ny, dest_i0, dest_j0, ierr)
      !! Degenerate-axis + along-edge slab geometry for an OBC-segment
      !! registration.  `i0`/`j0`/`nx`/`ny` are FILE-side (start, count);
      !! `dest_i0`/`dest_j0` are where the slab lands in the consumer's
      !! staging array.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: edge
      integer, intent(out) :: i0, j0, nx, ny, dest_i0, dest_j0
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`) on an unrecognised `edge`
         !! tag when present; absent behaves as today (`error stop`).

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      select case (edge)
      case (DATA_EDGE_WEST, DATA_EDGE_EAST)
         i0 = 1
         nx = 1
         j0 = grid%j_offset_global + 1
         ny = grid%ny_phys
         dest_i0 = 1
         dest_j0 = grid%nghost + 1
      case (DATA_EDGE_SOUTH, DATA_EDGE_NORTH)
         i0 = grid%i_offset_global + 1
         nx = grid%nx_phys
         j0 = 1
         ny = 1
         dest_i0 = grid%nghost + 1
         dest_j0 = 1
      case default
         call fail("ocean_data_input: unknown edge tag "//to_string(edge)// &
                   " (expected DATA_EDGE_WEST/EAST/SOUTH/NORTH)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end select
   end subroutine segment_geometry

   subroutine register_common(this, file, var, is_3d, i0, j0, nx, ny, nz, &
                              dest_i0, dest_j0, dest_n1, dest_n2, dest_n3, &
                              time_mode, cycle_period, t_offset, t_scale, scale, add_offset, &
                              oor, id, edge, ierr)
      !! Shared registration body for register_2d/_3d/_segment_2d/_segment_3d.
      !! `(i0, j0, nx, ny, nz)` are FILE-side start/count (already resolved
      !! by the caller — plain global-offset slicing for the base
      !! variants, degenerate-axis geometry for the segment variants).
      class(ocean_data_input_t), intent(inout) :: this
      character(len=*), intent(in) :: file, var
      logical, intent(in) :: is_3d
      integer, intent(in) :: i0, j0, nx, ny, nz
      integer, intent(in) :: dest_i0, dest_j0, dest_n1, dest_n2, dest_n3
      character(len=*), intent(in), optional :: time_mode
      real(wp), intent(in), optional :: cycle_period, t_offset, t_scale, scale, add_offset
      integer, intent(in), optional :: oor
      integer, intent(out) :: id
      integer, intent(in), optional :: edge
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`/`OCEAN_STATUS_ERR_IO`) on a
         !! registry/file/dimension failure when present; absent behaves
         !! as today (`error stop`).

      integer :: ncid, varid, tvarid, tdimid, status
      integer :: var_ndims
      integer :: var_dimids(4)
      integer :: expect_nd, d1, d2, dz, dt, nt
      character(len=:), allocatable :: units_str
      logical :: units_ok
      real(wp) :: t_scale_eff
      integer :: k
      integer :: local_ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      if (this%nfields >= this%nfields_max) then
         call fail("ocean_data_input: registry full (max_fields = "// &
                   to_string(this%nfields_max)//"); raise &ocean_data_nml max_fields", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      call nc_open_read(trim(file), ncid, ierr=local_ierr)
      if (.not. reg_io_ok(local_ierr, ierr)) return
      call nc_check(nf90_inq_varid(ncid, trim(var), varid), &
                    "ocean_data_input: finding variable '"//trim(var)//"' in "//trim(file), local_ierr)
      if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_variable(ncid, varid, ndims=var_ndims, &
                                          dimids=var_dimids(1:merge(4, 3, is_3d))), &
                    "ocean_data_input: querying variable '"//trim(var)//"'", local_ierr)
      if (.not. reg_io_ok(local_ierr, ierr, ncid)) return

      expect_nd = merge(4, 3, is_3d)
      if (var_ndims /= expect_nd) then
         call nc_close(ncid)
         call fail("ocean_data_input: '"//trim(var)//"' in "//trim(file)// &
                   " has "//to_string(var_ndims)//" dims; expected "// &
                   to_string(expect_nd)//" (x,y[,z],t) — Fortran storage order required, "// &
                   "no in-core transpose (regrid/reorder the file)", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      call nc_check(nf90_inquire_dimension(ncid, var_dimids(1), len=d1), &
                    "ocean_data_input: dim 1", local_ierr)
      if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(2), len=d2), &
                    "ocean_data_input: dim 2", local_ierr)
      if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
      if (is_3d) then
         call nc_check(nf90_inquire_dimension(ncid, var_dimids(3), len=dz), &
                       "ocean_data_input: dim 3 (z)", local_ierr)
         if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
         call nc_check(nf90_inquire_dimension(ncid, var_dimids(4), len=dt), &
                       "ocean_data_input: dim 4 (t)", local_ierr)
         if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
      else
         dz = 1
         call nc_check(nf90_inquire_dimension(ncid, var_dimids(3), len=dt), &
                       "ocean_data_input: dim 3 (t)", local_ierr)
         if (.not. reg_io_ok(local_ierr, ierr, ncid)) return
      end if

      if (present(edge)) then
         ! Degenerate-axis fail-loud check (§5.7.1): the horizontal
         ! extent this registration claims is 1 (a segment file) must
         ! actually BE 1 in the file.
         if (nx == 1 .and. d1 /= 1) then
            call nc_close(ncid)
            call fail("ocean_data_input: segment '"//trim(var)//"' edge "// &
                      to_string(edge)//" expects a degenerate x axis (extent 1); "// &
                      "file has extent "//to_string(d1), ierr, OCEAN_STATUS_ERR_IO)
            return
         end if
         if (ny == 1 .and. d2 /= 1) then
            call nc_close(ncid)
            call fail("ocean_data_input: segment '"//trim(var)//"' edge "// &
                      to_string(edge)//" expects a degenerate y axis (extent 1); "// &
                      "file has extent "//to_string(d2), ierr, OCEAN_STATUS_ERR_IO)
            return
         end if
      end if

      if (.not. dims_geometry_ok(var_ndims, expect_nd, d1, d2, i0 + nx - 1, j0 + ny - 1)) then
         call nc_close(ncid)
         call fail("ocean_data_input: '"//trim(var)//"' in "//trim(file)// &
                   " horizontal dims too small for this rank's slab: file ("// &
                   to_string(d1)//","//to_string(d2)//"), need offset+extent ("// &
                   to_string(i0 + nx - 1)//","//to_string(j0 + ny - 1)//")", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if
      if (is_3d .and. dz /= nz) then
         call nc_close(ncid)
         call fail("ocean_data_input: '"//trim(var)//"' in "//trim(file)// &
                   " has "//to_string(dz)//" z-levels; caller expects nz_src = "// &
                   to_string(nz), ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      nt = dt
      if (nt < 1) then
         call nc_close(ncid)
         call fail("ocean_data_input: '"//trim(var)//"' in "//trim(file)//" has no time records", &
                   ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      ! Time axis: the time DIMENSION and its coordinate VARIABLE share
      ! the CF name by convention; try that, then fall back to "time".
      tdimid = var_dimids(expect_nd)
      call resolve_time_var(ncid, tdimid, tvarid, status)
      if (status /= nf90_noerr) then
         call nc_close(ncid)
         call fail("ocean_data_input: no time coordinate variable found for '"// &
                   trim(var)//"' in "//trim(file), ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      id = this%nfields + 1
      this%nfields = id
      associate (fld => this%fields(id))
         fld%active = .true.
         fld%ncid = ncid
         fld%varid = varid
         fld%is_3d = is_3d
         fld%nx = nx
         fld%ny = ny
         fld%nz = merge(nz, 1, is_3d)
         fld%i0 = i0
         fld%j0 = j0
         fld%dest_i0 = dest_i0
         fld%dest_j0 = dest_j0
         fld%dest_n1 = dest_n1
         fld%dest_n2 = dest_n2
         fld%dest_n3 = dest_n3
         fld%nt = nt

         fld%scale = 1.0_wp
         if (present(scale)) fld%scale = scale
         fld%add_offset = 0.0_wp
         if (present(add_offset)) fld%add_offset = add_offset
         fld%t_offset = 0.0_wp
         if (present(t_offset)) fld%t_offset = t_offset
         fld%cycle_period = 0.0_wp
         if (present(cycle_period)) fld%cycle_period = cycle_period
         fld%oor = DATA_OOR_ERROR
         if (present(oor)) fld%oor = oor

         fld%time_mode = DATA_TIME_LINEAR
         if (present(time_mode)) fld%time_mode = data_input_time_mode_from_string(time_mode)
         if (fld%time_mode == DATA_TIME_CYCLIC .and. fld%cycle_period <= 0.0_wp) then
            fld%active = .false.
            this%nfields = this%nfields - 1
            call nc_close(ncid)
            call fail("ocean_data_input: '"//trim(var)//"' time_mode='cyclic' requires "// &
                      "cycle_period > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         allocate (fld%t_axis(nt))
         call nc_get_var_1d(ncid, tvarid, fld%t_axis)

         t_scale_eff = 1.0_wp
         if (present(t_scale)) then
            t_scale_eff = t_scale
         else
            call nc_get_att_text(ncid, tvarid, "units", units_str, units_ok)
            if (units_ok) then
               call data_input_time_scale_from_units(trim(units_str), t_scale_eff, units_ok)
            end if
            ! units_ok = .false. (absent/unrecognised attribute) silently
            ! keeps the seconds-default — the common case for a
            ! model-native file with no CF units attribute.
         end if
         fld%t_axis = fld%t_axis*t_scale_eff

         do k = 2, nt
            if (fld%t_axis(k) <= fld%t_axis(k - 1)) then
               call logger%error("ocean_data_input: '"//trim(var)//"' time axis not "// &
                                 "monotonically increasing at record "//to_string(k)//" ("// &
                                 to_string(fld%t_axis(k))//" <= "//to_string(fld%t_axis(k - 1))//")")
               error stop "ocean_data_input: non-monotonic time axis"
            end if
         end do

         allocate (fld%f0(fld%nx, fld%ny, fld%nz))
         allocate (fld%f1(fld%nx, fld%ny, fld%nz))

         if (fld%time_mode == DATA_TIME_STATIC) then
            call data_input_read_slab_impl(fld, 1)
            fld%f1 = fld%f0
            fld%rec0 = 1
            fld%rec1 = 1
            fld%w = 0.0_wp
         end if
      end associate
   end subroutine register_common

   function reg_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` reader calls in `register_common` into the caller's `ierr`
      !! contract: `.true.` on success; on failure, `.false.` with
      !! `ierr = OCEAN_STATUS_ERR_IO` when `ierr` is present (closing
      !! `ncid` first, when given, so a mid-registration failure does not
      !! leak the file handle), or `error stop`s with a generic message
      !! when `ierr` is absent — the SPECIFIC reason was already logged by
      !! `nc_check`/`fail` before this returns, so the legacy (no `ierr`)
      !! log output is unchanged; only the raw `error stop` text is
      !! generic (same idiom as `rdb_bathymetry::bathy_io_ok`, P0.1 F1/F2).
      integer, intent(in) :: local_ierr
      integer, intent(out), optional :: ierr
      integer, intent(in), optional :: ncid
      logical :: ok

      integer :: discard_ierr

      ok = (local_ierr == 0)
      if (ok) return

      if (present(ierr)) then
         if (present(ncid)) call nc_close(ncid, ierr=discard_ierr)
         ierr = OCEAN_STATUS_ERR_IO
         return
      end if

      error stop "ocean_data_input: NetCDF registration failure"
   end function reg_io_ok

   subroutine resolve_time_var(ncid, tdimid, tvarid, status)
      !! Look up the time coordinate variable by the CF convention (dim
      !! name == var name); fall back to a variable literally named
      !! "time". `status` is `nf90_noerr` iff found.
      integer, intent(in) :: ncid, tdimid
      integer, intent(out) :: tvarid, status
      character(len=64) :: dname

      status = nf90_inquire_dimension(ncid, tdimid, name=dname)
      if (status == nf90_noerr) status = nf90_inq_varid(ncid, trim(dname), tvarid)
      if (status /= nf90_noerr) status = nf90_inq_varid(ncid, "time", tvarid)
   end subroutine resolve_time_var

   ! ======================================================================
   ! Per-step path
   ! ======================================================================

   subroutine ocean_data_input_update_all(this, t)
      !! Driver hook: refresh every registered field's bracket for model
      !! time `t` (reading + pushing a new slab to the device only when
      !! the bracket actually advances).  No-op when `nfields == 0` —
      !! every shipped namelist today.  Does NOT write into any
      !! consumer array (see module docstring) — call `update_2d/_3d`
      !! afterwards for that.
      class(ocean_data_input_t), intent(inout) :: this
      real(wp), intent(in) :: t
      integer :: i
      if (this%nfields == 0) return
      do i = 1, this%nfields
         call data_input_refresh_brackets(this, i, t)
      end do
   end subroutine ocean_data_input_update_all

   subroutine data_input_refresh_brackets(this, id, t)
      !! Host-side registry walk (outer shim) for field `id`.
      class(ocean_data_input_t), intent(inout) :: this
      integer, intent(in) :: id
      real(wp), intent(in) :: t
      integer :: n0, n1
      real(wp) :: w, t_query
      logical :: oor

      associate (fld => this%fields(id))
         if (fld%time_mode == DATA_TIME_STATIC) then
            fld%t_last = t
            return
         end if

         t_query = t + fld%t_offset
         call data_input_locate(fld%t_axis, fld%nt, fld%time_mode, fld%cycle_period, &
                                t_query, n0, n1, w, oor)

         if (oor) then
            if (fld%oor == DATA_OOR_ERROR) then
               call logger%error("ocean_data_input: query time "//to_string(t_query)// &
                                 "s is outside the file time axis ["// &
                                 to_string(fld%t_axis(1))//", "// &
                                 to_string(fld%t_axis(fld%nt))//"]s (field id "// &
                                 to_string(id)//"); pass oor=DATA_OOR_CLAMP to clamp instead")
               error stop "ocean_data_input: query time out of range"
            end if
            if (.not. fld%oor_warned) then
               call logger%warning("ocean_data_input: field id "//to_string(id)// &
                                   "clamped at the time-axis boundary (query time "// &
                                   to_string(t_query)//"s out of range)")
               fld%oor_warned = .true.
            end if
         end if

         ! `!$acc update device(...)` is only correct HERE (the runtime
         ! per-step path, always called after `ocean_state_enter_data`
         ! has mapped `fld%f0`/`fld%f1`) — never inside
         ! `data_input_read_slab_impl` itself, which is also called from
         ! `register_common` at SETUP time for a DATA_TIME_STATIC field,
         ! before anything is device-mapped (CLAUDE.md gotcha: an
         ! `update device` on an unmapped array is undefined).
         if (n0 /= fld%rec0) then
            call data_input_read_slab_impl(fld, n0, into_f1=.false.)
            !$acc update device(fld%f0)
         end if
         if (n1 /= fld%rec1) then
            if (n1 == n0) then
               fld%f1 = fld%f0
               !$acc update device(fld%f1)
            else
               call data_input_read_slab_impl(fld, n1, into_f1=.true.)
               !$acc update device(fld%f1)
            end if
         end if
         if (n0 /= fld%rec0 .or. n1 /= fld%rec1) then
            if (this%verbose) then
               call logger%info("ocean_data_input: field id "//to_string(id)// &
                                " bracket advanced to records ("//to_string(n0)//", "// &
                                to_string(n1)//"), w = "//to_string(w)// &
                                ", t = "//to_string(t))
            end if
         end if
         fld%rec0 = n0
         fld%rec1 = n1
         fld%w = w
         fld%t_last = t
      end associate
   end subroutine data_input_refresh_brackets

   subroutine data_input_read_slab_impl(fld, rec, into_f1)
      !! Host-ONLY NetCDF slab read for file record `rec` -> `fld%f0`
      !! (default) or `fld%f1` (`into_f1 = .true.`), applying
      !! `scale`/`add_offset` once at read.  Uses the module-level host
      !! workspace (registration/bracket-advance-time only — never a
      !! per-step allocation).  Deliberately carries NO `!$acc` directive:
      !! called both before `enter_data` (STATIC field, at registration)
      !! and after it (LINEAR/CYCLIC bracket advance) — the caller pushes
      !! to the device itself, only when that is actually correct
      !! (`data_input_refresh_brackets`).
      type(data_input_field_t), intent(inout) :: fld
      integer, intent(in) :: rec
      logical, intent(in), optional :: into_f1
      integer :: start4(4), count4(4)
      logical :: to_f1

      to_f1 = .false.
      if (present(into_f1)) to_f1 = into_f1

      call data_input_workspace_ensure(fld%nx, fld%ny, fld%nz)

      if (fld%is_3d) then
         start4 = [fld%i0, fld%j0, 1, rec]
         count4 = [fld%nx, fld%ny, fld%nz, 1]
         call nc_get_var_slab_3d(fld%ncid, fld%varid, start4, count4, data_input_ws)
      else
         start4(1:3) = [fld%i0, fld%j0, rec]
         count4(1:3) = [fld%nx, fld%ny, 1]
         call nc_get_var_slab_3d(fld%ncid, fld%varid, start4(1:3), count4(1:3), data_input_ws)
      end if

      data_input_ws = fld%scale*data_input_ws + fld%add_offset

      if (to_f1) then
         fld%f1 = data_input_ws(1:fld%nx, 1:fld%ny, 1:fld%nz)
      else
         fld%f0 = data_input_ws(1:fld%nx, 1:fld%ny, 1:fld%nz)
      end if
      fld%nreads = fld%nreads + 1
   end subroutine data_input_read_slab_impl

   subroutine data_input_workspace_ensure(nx, ny, nz)
      integer, intent(in) :: nx, ny, nz
      logical :: need_alloc
      need_alloc = .false.
      if (.not. allocated(data_input_ws)) then
         need_alloc = .true.
      else if (size(data_input_ws, 1) /= nx .or. size(data_input_ws, 2) /= ny .or. &
               size(data_input_ws, 3) /= nz) then
         deallocate (data_input_ws)
         need_alloc = .true.
      end if
      if (need_alloc) allocate (data_input_ws(nx, ny, nz))
   end subroutine data_input_workspace_ensure

   subroutine data_input_workspace_cleanup()
      if (allocated(data_input_ws)) deallocate (data_input_ws)
   end subroutine data_input_workspace_cleanup

   subroutine ocean_data_input_update_2d(this, id, t, n1, n2, dest)
      !! Blend field `id`'s current bracket into the caller's WHOLE,
      !! device-mapped `dest(n1, n2)` array at the registration-time
      !! offset.  `t` must equal the value `update_all` most recently
      !! refreshed this field with (fail-loud ordering check — a
      !! consumer calling this before `update_all` has run for the
      !! current step is a real bug, not a silent stale read).  Exempt
      !! for `DATA_TIME_STATIC` fields, whose bracket never changes.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id, n1, n2
      real(wp), intent(in) :: t
      real(wp), intent(inout) :: dest(n1, n2)

      call check_registered(this, id, is_3d=.false., n1=n1, n2=n2)
      call check_fresh(this, id, t)
      call data_input_blend_2d_impl(this%fields(id)%nx, this%fields(id)%ny, n1, n2, &
                                    this%fields(id)%dest_i0, this%fields(id)%dest_j0, &
                                    this%fields(id)%f0, this%fields(id)%f1, &
                                    this%fields(id)%w, dest)
   end subroutine ocean_data_input_update_2d

   subroutine ocean_data_input_update_3d(this, id, t, n1, n2, n3, dest)
      !! 3-D twin of `update_2d`.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id, n1, n2, n3
      real(wp), intent(in) :: t
      real(wp), intent(inout) :: dest(n1, n2, n3)

      call check_registered(this, id, is_3d=.true., n1=n1, n2=n2, n3=n3)
      call check_fresh(this, id, t)
      call data_input_blend_3d_impl(this%fields(id)%nx, this%fields(id)%ny, this%fields(id)%nz, &
                                    n1, n2, n3, &
                                    this%fields(id)%dest_i0, this%fields(id)%dest_j0, &
                                    this%fields(id)%f0, this%fields(id)%f1, &
                                    this%fields(id)%w, dest)
   end subroutine ocean_data_input_update_3d

   subroutine ocean_data_input_fill_static_host(this, id, n1, n2, dest)
      !! One-shot HOST-side fill of a `DATA_TIME_STATIC` field.  Copies
      !! the already-read record-1 slab into `dest` at the
      !! registration-time offsets.  Touches NO device memory —
      !! callable at setup, before `dest` is mapped.  Plain host `do`
      !! loops (NOT `do concurrent` — on `-stdpar=gpu` a bare `do
      !! concurrent` is unconditionally offloaded regardless of whether
      !! `dest` is mapped, which is exactly the silent-stale-write bug
      !! this routine exists to avoid).  Fails loud for any non-static
      !! field.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id, n1, n2
      real(wp), intent(inout) :: dest(n1, n2)
      integer :: i, j

      call check_registered(this, id, is_3d=.false., n1=n1, n2=n2)
      if (this%fields(id)%time_mode /= DATA_TIME_STATIC) then
         call logger%error("ocean_data_input: fill_static_host called on field id "// &
                           to_string(id)//", which is not DATA_TIME_STATIC")
         error stop "ocean_data_input: fill_static_host on a non-static field"
      end if
      do j = 1, this%fields(id)%ny
         do i = 1, this%fields(id)%nx
            dest(this%fields(id)%dest_i0 + i - 1, this%fields(id)%dest_j0 + j - 1) = &
               this%fields(id)%f0(i, j, 1)
         end do
      end do
   end subroutine ocean_data_input_fill_static_host

   subroutine ocean_data_input_fill_static_host_3d(this, id, n1, n2, n3, dest)
      !! 3-D twin of `fill_static_host`.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id, n1, n2, n3
      real(wp), intent(inout) :: dest(n1, n2, n3)
      integer :: i, j, k

      call check_registered(this, id, is_3d=.true., n1=n1, n2=n2, n3=n3)
      if (this%fields(id)%time_mode /= DATA_TIME_STATIC) then
         call logger%error("ocean_data_input: fill_static_host_3d called on field id "// &
                           to_string(id)//", which is not DATA_TIME_STATIC")
         error stop "ocean_data_input: fill_static_host_3d on a non-static field"
      end if
      do k = 1, this%fields(id)%nz
         do j = 1, this%fields(id)%ny
            do i = 1, this%fields(id)%nx
               dest(this%fields(id)%dest_i0 + i - 1, this%fields(id)%dest_j0 + j - 1, k) = &
                  this%fields(id)%f0(i, j, k)
            end do
         end do
      end do
   end subroutine ocean_data_input_fill_static_host_3d

   subroutine check_registered(this, id, is_3d, n1, n2, n3)
      !! Fail-loud guard shared by every per-step/fill accessor: `id`
      !! must name an active field of the right rank, and the caller's
      !! `dest` shape must match what was declared at registration.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id, n1, n2
      logical, intent(in) :: is_3d
      integer, intent(in), optional :: n3
      integer :: dn3

      if (id < 1 .or. id > this%nfields) then
         call logger%error("ocean_data_input: invalid field id "//to_string(id))
         error stop "ocean_data_input: invalid field id"
      end if
      if (.not. this%fields(id)%active) then
         call logger%error("ocean_data_input: field id "//to_string(id)//" is not active")
         error stop "ocean_data_input: field id not active"
      end if
      if (this%fields(id)%is_3d .neqv. is_3d) then
         call logger%error("ocean_data_input: field id "//to_string(id)// &
                           " rank mismatch (registered "// &
                           merge("3-D", "2-D", this%fields(id)%is_3d)//", called as "// &
                           merge("3-D", "2-D", is_3d)//")")
         error stop "ocean_data_input: field rank mismatch"
      end if
      dn3 = 1
      if (present(n3)) dn3 = n3
      if (n1 /= this%fields(id)%dest_n1 .or. n2 /= this%fields(id)%dest_n2 .or. &
          dn3 /= this%fields(id)%dest_n3) then
         call logger%error("ocean_data_input: field id "//to_string(id)// &
                           " dest shape ("//to_string(n1)//","//to_string(n2)//","// &
                           to_string(dn3)//") does not match the registered shape ("// &
                           to_string(this%fields(id)%dest_n1)//","// &
                           to_string(this%fields(id)%dest_n2)//","// &
                           to_string(this%fields(id)%dest_n3)//")")
         error stop "ocean_data_input: dest shape mismatch"
      end if
   end subroutine check_registered

   subroutine check_fresh(this, id, t)
      !! Ordering guard: `update_2d/_3d` must be called with the same
      !! `t` `update_all` most recently refreshed this field's bracket
      !! with (skipped for STATIC — its bracket never changes, so
      !! "freshness" is meaningless).  Catches a consumer calling
      !! `update_2d/_3d` before `update_all` has run this step, which
      !! would otherwise silently blend a stale bracket.
      class(ocean_data_input_t), intent(in) :: this
      integer, intent(in) :: id
      real(wp), intent(in) :: t
      if (this%fields(id)%time_mode == DATA_TIME_STATIC) return
      if (t /= this%fields(id)%t_last) then
         call logger%error("ocean_data_input: field id "//to_string(id)// &
                           " queried at t = "//to_string(t)// &
                           "s but its bracket was last refreshed at t = "// &
                           to_string(this%fields(id)%t_last)// &
                           "s — call ocean_data_input_update_all(this, t) first")
         error stop "ocean_data_input: stale bracket (update_all not called for this t)"
      end if
   end subroutine check_fresh

end module rdb_ocean_data_input
