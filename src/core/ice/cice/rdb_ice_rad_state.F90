!! Radiation state slot of the CICE-port sea ice (roundabout-original).
module rdb_ice_rad_state
   !! `ice_rad_state_t` holds the per-step OUTPUTS of the CCSM3 shortwave driver
   !! (`rdb_ice_radiation`; Briegleb et al. 2004, NCAR/TN-463+STR), laid out as
   !! structure-of-arrays with the horizontal indices FIRST (MASTER_PLAN C.1):
   !!
   !! | array | shape | unit | role |
   !! |---|---|---|---|
   !! | `alvn, alin` | (nx,ny,ncat) | 1 | category albedo, visible / near-IR (direct == diffuse here: ONE value per band; the four coupler names are aliases at the seam) |
   !! | `fswsfcn, fswintn, fswthrun` | (nx,ny,ncat) | W m-2 of category area | absorbed at the surface / in the interior / transmitted to the ocean |
   !! | `Iswabsn` | (nx,ny,ncat,nilyr) | W m-2 | absorbed in ice layer k (`k = 1` TOP) |
   !! | `alv_ai, ali_ai` | (nx,ny) | 1 | partial sums `sum_n aicen*albedo_n` over the ice categories: not divided by the ice area, no open water (see below) |
   !! | `errn` | int32 (nx,ny,ncat) | | per-category error word, merged into the cell word by the aggregate sweep |
   !!
   !! **Aggregate.**  A grid-cell albedo is `alv_ai + aice0*(ocean albedo)`; an ice-only mean is `alv_ai/aice`.
   !!
   !! There is no snow-layer absorption array: `Sswabs` is identically zero under CCSM3 and the
   !! wiring PR passes a zero snow absorption to the thermodynamics.
   !!
   !! **Not restart state.**  With `sw_lag = .false.` (the only mode, MASTER_PLAN
   !! C.6-9) every array is RECOMPUTED from the current ice state and the current
   !! shortwave at the start of each step, before therm1 reads it; there is no
   !! `prep_radiation` rescale and no previous-step value is ever consumed.  So
   !! nothing here is registered with the restart registry, and the slot is not
   !! composed into the live ice state yet (the wiring PR does that).
   !!
   !! **Device lifecycle** (`mem:separate`): `enter_data` `copyin`s the COMPONENT
   !! arrays (host zeros: defined device contents, never `create`d garbage);
   !! `exit_data` deletes them (per-step scratch; inspect with `!$acc update self`
   !! first).  The aggregate type is never mapped.
   !!
   !! **Setup.**  `init` then `configure(ahmax, snowpatch)` (once), then `enter_data`.
   !!
   !! **Ice core.**  The slot takes explicit extents (not an `hgrid_t`) so that it
   !! lives in `rdb_ice_core_objs`; its `bytes()` therefore counts with
   !! `storage_size` rather than the main library's `arr_bytes`.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX
   use rdb_ice_shortwave, only: ice_albedo_fhtan
   implicit none
   private

   public :: ice_rad_state_t
   public :: ICE_RAD_INIT_OK, ICE_RAD_INIT_BAD_EXTENT

   integer, parameter :: ICE_RAD_INIT_OK = 0
      !! `init` allocated the slot.
   integer, parameter :: ICE_RAD_INIT_BAD_EXTENT = 1
      !! An extent is non-positive or beyond the `rdb_ice_dims` bounds; nothing allocated.

   type :: ice_rad_state_t
      !! The radiation slot (see the module docstring).
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data`.
      integer :: nx = 0
         !! Total extent in x including ghosts.
      integer :: ny = 0
         !! Total extent in y including ghosts.
      integer :: ncat = 0
         !! Thickness categories.
      integer :: nilyr = 0
         !! Ice layers per category.
      logical :: configured = .false.
         !! True after `configure`.
      real(wp) :: fhtan = 0.0_wp
         !! `atan(4 ahmax)`, the thin-ice ramp normaliser: per configuration, computed once by `configure`.
      real(wp) :: snowpatch = 0.0_wp
         !! Snow-cover length scale (m), set by `configure`.
      real(wp), allocatable :: alvn(:, :, :)
         !! Visible albedo per category (nx,ny,ncat); direct == diffuse.
      real(wp), allocatable :: alin(:, :, :)
         !! Near-IR albedo per category.
      real(wp), allocatable :: fswsfcn(:, :, :)
         !! Shortwave absorbed at the surface [W m-2 of category area].
      real(wp), allocatable :: fswintn(:, :, :)
         !! Shortwave absorbed in the ice interior [W m-2].
      real(wp), allocatable :: fswthrun(:, :, :)
         !! Shortwave transmitted to the ocean [W m-2]; all visible under CCSM3.
      real(wp), allocatable :: Iswabsn(:, :, :, :)
         !! Absorbed in ice layer k (nx,ny,ncat,nilyr), `k = 1` TOP [W m-2].
      real(wp), allocatable :: alv_ai(:, :)
         !! Visible partial sum `sum_n aicen*alvn` (nx,ny); no division by `aice`, no open water.
      real(wp), allocatable :: ali_ai(:, :)
         !! Near-IR partial sum `sum_n aicen*alin`.
      integer(int32), allocatable :: errn(:, :, :)
         !! Per-category error word (nx,ny,ncat), bits of `rdb_ice_errors`.
   contains
      procedure, non_overridable :: init => rad_state_init
      procedure, non_overridable :: configure => rad_state_configure
      procedure, non_overridable :: destroy => rad_state_destroy
      procedure, non_overridable :: enter_data => rad_state_enter_data
      procedure, non_overridable :: exit_data => rad_state_exit_data
      procedure, non_overridable :: bytes => rad_state_bytes
   end type ice_rad_state_t

contains

   subroutine rad_state_init(this, nx, ny, ncat, nilyr, ierr)
      !! Allocate and zero the slot.  `ierr` is `ICE_RAD_INIT_*`; on a non-zero value
      !! nothing is allocated.
      class(ice_rad_state_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, nilyr
      integer, intent(out) :: ierr

      ierr = ICE_RAD_INIT_OK
      if (nx < 1 .or. ny < 1 .or. ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 &
          .or. nilyr > ICE_NILYR_MAX) then
         ierr = ICE_RAD_INIT_BAD_EXTENT
         return
      end if
      if (this%is_init) call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%nilyr = nilyr
      this%configured = .false.
      allocate (this%alvn(nx, ny, ncat), this%alin(nx, ny, ncat), source=0.0_wp)
      allocate (this%fswsfcn(nx, ny, ncat), this%fswintn(nx, ny, ncat), this%fswthrun(nx, ny, ncat), &
                source=0.0_wp)
      allocate (this%Iswabsn(nx, ny, ncat, nilyr), source=0.0_wp)
      allocate (this%alv_ai(nx, ny), this%ali_ai(nx, ny), source=0.0_wp)
      allocate (this%errn(nx, ny, ncat), source=0_int32)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine rad_state_init

   subroutine rad_state_configure(this, ahmax, snowpatch)
      !! Evaluate the per-configuration constants ONCE (`atan(4 ahmax)`) and keep `snowpatch`.
      class(ice_rad_state_t), intent(inout) :: this
      real(wp), intent(in) :: ahmax
         !! Thickness above which the ice albedo is constant (m).
      real(wp), intent(in) :: snowpatch
         !! Snow-cover length scale (m).
      this%fhtan = ice_albedo_fhtan(ahmax)
      this%snowpatch = snowpatch
      this%configured = .true.
   end subroutine rad_state_configure

   subroutine rad_state_destroy(this)
      !! Reverse of `init` (unmap first).
      class(ice_rad_state_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      this%configured = .false.
      if (allocated(this%alvn)) deallocate (this%alvn)
      if (allocated(this%alin)) deallocate (this%alin)
      if (allocated(this%fswsfcn)) deallocate (this%fswsfcn)
      if (allocated(this%fswintn)) deallocate (this%fswintn)
      if (allocated(this%fswthrun)) deallocate (this%fswthrun)
      if (allocated(this%Iswabsn)) deallocate (this%Iswabsn)
      if (allocated(this%alv_ai)) deallocate (this%alv_ai)
      if (allocated(this%ali_ai)) deallocate (this%ali_ai)
      if (allocated(this%errn)) deallocate (this%errn)
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%nilyr = 0
   end subroutine rad_state_destroy

   subroutine rad_state_enter_data(this)
      !! Map the slot to the device (poly TBP delegating to a `type(...)` impl).
      class(ice_rad_state_t), intent(inout) :: this
      select type (this)
      type is (ice_rad_state_t)
         call rad_state_enter_data_impl(this)
      end select
   end subroutine rad_state_enter_data

   subroutine rad_state_enter_data_impl(this)
      !! `copyin` every component array (host zeros).  COMPONENT arrays only.
      type(ice_rad_state_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%alvn, this%alin)
      !$acc enter data copyin(this%fswsfcn, this%fswintn, this%fswthrun)
      !$acc enter data copyin(this%Iswabsn)
      !$acc enter data copyin(this%alv_ai, this%ali_ai, this%errn)
      this%mapped = .true.
   end subroutine rad_state_enter_data_impl

   subroutine rad_state_exit_data(this)
      !! Unmap the slot (poly TBP delegating to a `type(...)` impl).
      class(ice_rad_state_t), intent(inout) :: this
      select type (this)
      type is (ice_rad_state_t)
         call rad_state_exit_data_impl(this)
      end select
   end subroutine rad_state_exit_data

   subroutine rad_state_exit_data_impl(this)
      !! `delete` every component array (per-step scratch, nothing to copy out).
      type(ice_rad_state_t), intent(inout) :: this
      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%alvn, this%alin)
      !$acc exit data delete(this%fswsfcn, this%fswintn, this%fswthrun)
      !$acc exit data delete(this%Iswabsn)
      !$acc exit data delete(this%alv_ai, this%ali_ai, this%errn)
      this%mapped = .false.
   end subroutine rad_state_exit_data_impl

   pure function rad_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the slot (one term per array).
      class(ice_rad_state_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = rb(size(this%alvn, kind=int64), storage_size(this%alvn)) &
               + rb(size(this%alin, kind=int64), storage_size(this%alin)) &
               + rb(size(this%fswsfcn, kind=int64), storage_size(this%fswsfcn)) &
               + rb(size(this%fswintn, kind=int64), storage_size(this%fswintn)) &
               + rb(size(this%fswthrun, kind=int64), storage_size(this%fswthrun)) &
               + rb(size(this%Iswabsn, kind=int64), storage_size(this%Iswabsn)) &
               + rb(size(this%alv_ai, kind=int64), storage_size(this%alv_ai)) &
               + rb(size(this%ali_ai, kind=int64), storage_size(this%ali_ai)) &
               + rb(size(this%errn, kind=int64), storage_size(this%errn))
   end function rad_state_bytes

   pure function rb(n, bits) result(b)
      !! Bytes of `n` elements of `bits` bits each.
      integer(int64), intent(in) :: n
      integer, intent(in) :: bits
      integer(int64) :: b
      b = n*int(bits, int64)/8_int64
   end function rb

end module rdb_ice_rad_state
