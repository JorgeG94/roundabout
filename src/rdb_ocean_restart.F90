!! Ocean restart / checkpoint manager + per-slot field registry.
module rdb_ocean_restart
   !! Restart cadence + per-slot checkpoint registry for the ocean
   !! dyn-core.  Each prognostic-owning slot registers its arrays at
   !! birth; the manager walks the registry to do the I/O, so new
   !! state-carrying slots join by registering rather than editing here.
   !!
   !! Contract (MPI-native):
   !!   1. FULL local arrays — every registered array is written/read at
   !!      its TOTAL local extent (interior + ghosts).  Physical wall
   !!      ghosts are genuine owned boundary state (nothing rebuilds them),
   !!      so an interior-only write breaks the bit-exact gate; periodic/
   !!      fold/halo seam ghosts are redundant but harmless to save.
   !!   2. Decomposition + grid metadata (px, py, dims, i/j_start, nz,
   !!      nghost, vcoord, tracer set) live in the file; resume requires
   !!      the SAME decomposition and error-stops on mismatch.  Cross-rank
   !!      redistribution is an offline tool, never in-Fortran.
   !!   3. Global scalars (time, step, outer_step_count) are rank-0's.
   !!   4. A registered field is REQUIRED by default (missing on read is
   !!      FATAL); mark genuinely optional entries `optional=.true.`.
   !!
   !! Device sync: `device_mapped` arrays (default) are pulled with
   !! `!$acc update self` before the host NetCDF write; host-only state
   !! (`device_mapped=.false.`) is skipped so `update self` is never
   !! issued on an unmapped array.  The read path runs BEFORE
   !! `ocean_state_enter_data`, writing the host interior so the
   !! subsequent enter_data carries values up to the device.
   use rdb_constants, only: wp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_restart_t
   public :: restart_registry_t
   public :: restart_entry_t
   public :: RESTART_SCHEMA_VERSION

   integer, parameter :: RESTART_SCHEMA_VERSION = 1
      !! Bumped on any non-back-compatible on-disk field-set/layout
      !! change.  Written as a global attribute, validated on read.

   integer, parameter :: MAX_RESTART_ENTRIES = 256
      !! Fixed cap on registered fields; bump if a slot family exceeds it.

   type :: restart_entry_t
      !! One registered checkpoint field.  Holds a pointer to the slot's
      !! host array (or host scalar) plus the interior extents needed to
      !! slice owned cells.  Exactly one of `p0 / p2 / p3` is associated
      !! (per `rank`: 0 = host scalar, 2/3 = array).
      character(len=64) :: tag = ""
         !! NetCDF variable name (unique within the file).
      integer :: rank = 0
         !! 0 (host scalar), 2, or 3.
      real(wp), pointer :: p0 => null()
         !! Host scalar (rank-0 persistent state, e.g. Chapman eta_old).
      real(wp), pointer :: p2(:, :) => null()
         !! Host array for rank-2 fields (incl. ghosts).
      real(wp), pointer :: p3(:, :, :) => null()
         !! Host array for rank-3 fields (incl. ghosts).
      integer :: ng = 0
         !! Ghost width (interior starts at ng+1 in x and y).
      integer :: nx_phys = 0, ny_phys = 0
         !! Owned-cell extents in x, y.
      integer :: nk = 0
         !! Third-dim extent for rank-3 fields (no vertical ghosts —
         !! layers/interfaces are all owned).
      logical :: optional = .false.
         !! .true. => read path warns-and-seeds if absent.  Default
         !! .false. => a missing field is FATAL.
      logical :: device_mapped = .true.
         !! .true. => write path pulls host-ward via `!$acc update self`
         !! first.  Host-only state sets .false. so `update self` is never
         !! issued on an unmapped array (crashes on GPU).
   end type restart_entry_t

   type :: restart_registry_t
      !! Append-only registry of checkpoint fields.  Slots call
      !! `register_2d` / `register_3d` from their (or the state's) init.
      integer :: n = 0
      type(restart_entry_t) :: entries(MAX_RESTART_ENTRIES)
   contains
      procedure, non_overridable :: register_scalar => registry_register_scalar
      procedure, non_overridable :: register_2d => registry_register_2d
      procedure, non_overridable :: register_3d => registry_register_3d
      procedure, non_overridable :: clear => registry_clear
   end type restart_registry_t

   type :: ocean_restart_t
      !! Restart-manager lifecycle marker.  The registry is built per-run
      !! by `ocean_state_build_restart_registry` (it holds live pointers
      !! into the freshly allocated slots), not stored here.  Cadence is
      !! driven by the driver off `cfg%restart_interval`/`cfg%output_dir`.
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
   contains
      procedure, non_overridable :: init => ocean_restart_init
      procedure, non_overridable :: destroy => ocean_restart_destroy
   end type ocean_restart_t

contains

   subroutine ocean_restart_init(this)
      class(ocean_restart_t), intent(inout) :: this
      this%is_init = .true.
   end subroutine ocean_restart_init

   subroutine ocean_restart_destroy(this)
      class(ocean_restart_t), intent(inout) :: this
      this%is_init = .false.
   end subroutine ocean_restart_destroy

   subroutine registry_clear(this)
      class(restart_registry_t), intent(inout) :: this
      integer :: i
      do i = 1, this%n
         this%entries(i)%p0 => null()
         this%entries(i)%p2 => null()
         this%entries(i)%p3 => null()
         this%entries(i)%tag = ""
         this%entries(i)%rank = 0
         this%entries(i)%optional = .false.
         this%entries(i)%device_mapped = .true.
      end do
      this%n = 0
   end subroutine registry_clear

   subroutine registry_register_scalar(this, tag, scal, optional)
      !! Register a host scalar (rank-0 persistent state, e.g. the Chapman
      !! `eta_old_chapman_*` corner values).  Host-only — never device-
      !! mapped, so the write path skips `update self` for it.
      class(restart_registry_t), intent(inout) :: this
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: scal
      logical, intent(in), optional :: optional

      if (this%n >= MAX_RESTART_ENTRIES) then
         call logger%error("restart registry full ("//to_string(MAX_RESTART_ENTRIES)// &
                           " entries); bump MAX_RESTART_ENTRIES")
         error stop "restart registry overflow"
      end if
      this%n = this%n + 1
      this%entries(this%n)%tag = tag
      this%entries(this%n)%rank = 0
      this%entries(this%n)%p0 => scal
      this%entries(this%n)%device_mapped = .false.
      if (present(optional)) this%entries(this%n)%optional = optional
   end subroutine registry_register_scalar

   subroutine registry_register_2d(this, tag, arr, ng, nx_phys, ny_phys, optional, device_mapped)
      !! Register a rank-2 owned field.  `arr` is the full (ghosted)
      !! host array; the owned slice is `(ng+1:ng+nx_phys, ng+1:ng+ny_phys)`.
      class(restart_registry_t), intent(inout) :: this
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: arr(:, :)
      integer, intent(in) :: ng, nx_phys, ny_phys
      logical, intent(in), optional :: optional, device_mapped

      if (this%n >= MAX_RESTART_ENTRIES) then
         call logger%error("restart registry full ("//to_string(MAX_RESTART_ENTRIES)// &
                           " entries); bump MAX_RESTART_ENTRIES")
         error stop "restart registry overflow"
      end if
      this%n = this%n + 1
      this%entries(this%n)%tag = tag
      this%entries(this%n)%rank = 2
      this%entries(this%n)%p2 => arr
      this%entries(this%n)%ng = ng
      this%entries(this%n)%nx_phys = nx_phys
      this%entries(this%n)%ny_phys = ny_phys
      this%entries(this%n)%nk = 0
      if (present(optional)) this%entries(this%n)%optional = optional
      if (present(device_mapped)) this%entries(this%n)%device_mapped = device_mapped
   end subroutine registry_register_2d

   subroutine registry_register_3d(this, tag, arr, ng, nx_phys, ny_phys, optional, device_mapped)
      !! Register a rank-3 owned field.  The vertical extent is taken
      !! from `size(arr,3)` (layers or interfaces — both fully owned).
      class(restart_registry_t), intent(inout) :: this
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: arr(:, :, :)
      integer, intent(in) :: ng, nx_phys, ny_phys
      logical, intent(in), optional :: optional, device_mapped

      if (this%n >= MAX_RESTART_ENTRIES) then
         call logger%error("restart registry full ("//to_string(MAX_RESTART_ENTRIES)// &
                           " entries); bump MAX_RESTART_ENTRIES")
         error stop "restart registry overflow"
      end if
      this%n = this%n + 1
      this%entries(this%n)%tag = tag
      this%entries(this%n)%rank = 3
      this%entries(this%n)%p3 => arr
      this%entries(this%n)%ng = ng
      this%entries(this%n)%nx_phys = nx_phys
      this%entries(this%n)%ny_phys = ny_phys
      this%entries(this%n)%nk = size(arr, 3)
      if (present(optional)) this%entries(this%n)%optional = optional
      if (present(device_mapped)) this%entries(this%n)%device_mapped = device_mapped
   end subroutine registry_register_3d

end module rdb_ocean_restart
