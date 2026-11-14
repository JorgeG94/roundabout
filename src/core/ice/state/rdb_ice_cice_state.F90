!! State slot of the CICE-port sea ice (Icepack layout; roundabout-original).
module rdb_ice_cice_state
   !! `ice_cice_state_t` holds the prognostic state of the CICE-port ice
   !! (`&ocean_ice_nml scheme = "cice"`), laid out as named structure-of-arrays
   !! with the HORIZONTAL indices FIRST and the ocean's ghost extents
   !! (`nx_total, ny_total`), exactly the table of MASTER_PLAN C.1:
   !!
   !! | array | shape | unit | role |
   !! |---|---|---|---|
   !! | `aicen, vicen, vsnon` | (nx,ny,ncat) | 1, m, m | category area fraction; ice and snow volume per CELL area |
   !! | `Tsfcn` | (nx,ny,ncat) | degC | surface temperature (area-weighted tracer, persistent) |
   !! | `qicen, sicen` | (nx,ny,ncat,nilyr) | J/m3 of ice, ppt | layer enthalpy and salinity, `k = 1` at the TOP |
   !! | `qsnon` | (nx,ny,ncat,nslyr) | J/m3 of snow | snow enthalpy |
   !! | `aice0` | (nx,ny) | 1 | open water, `1 - sum(aicen)` (rebuilt by `ice_update_state`) |
   !! | `aice, vice, vsno` | (nx,ny) | 1, m, m | cell aggregates |
   !! | `aice_init, aicen_init, vicen_init, vsnon_init, Tsfcn_init` | | | start-of-step copies (scratch; written by the step driver) |
   !! | `hin_max(0:ncat)` | | m | category bounds |
   !! | `ice_err(nx,ny)` | int32 | | per-cell error word (`rdb_ice_errors`) |
   !!
   !! **Layer order.**  The ice column is TOP-DOWN (`k = 1` under the surface,
   !! `k = nilyr` at the ice base): the documented exception to the ocean's
   !! bottom-up rule (see `rdb_ice_dims` and CLAUDE.md).  `sicen` is stored
   !! although BL99 never changes it (a prescribed profile), so that every
   !! tracer routine, the ITD and the restart are index-identical to the
   !! reference.
   !!
   !! **Coexistence with the SIS2 ice.**  The slot is composed into
   !! `ocean_sea_ice_t%cice` but is never initialised unless
   !! `scheme = "cice"` is selected (which cannot yet be done at run time), so
   !! under the default `scheme = "sis2"` nothing here is allocated, mapped or
   !! registered.
   !!
   !! **Device lifecycle** (`mem:separate`): `enter_data` maps the COMPONENT
   !! arrays (never the aggregate type, whose allocatable descriptors would be
   !! overwritten by an aggregate update) and the parameter bundle `p` ONCE
   !! (all scalars).  Kernels reach `p` ONLY by receiving a sub-object as a
   !! dummy argument (`st%p%c`), never by writing `this%p%c%rhoi` inside a
   !! `do concurrent` (that would implicitly copy the whole parent, the
   !! forbidden aggregate pattern); and the parent aggregate must never be
   !! mapped later (its `p` sub-object is already present: "partially present").
   !! Persistent state, the aggregates and `ice_err` are `copyin` (host
   !! initialised or restart-seeded first); the start-of-step copies are
   !! `copyin` (host zeros until `ice_save_init` writes them).  `exit_data`
   !! copies the persistent state and aggregates back and deletes the rest.
   !!
   !! **Aggregates and warm starts (for the driver PRs).**  `aice0, aice, vice,
   !! vsno` are NOT restart state: after a warm-start read they arrive as the
   !! `init` values (0, 0, 0, `aice0 = 1`), so the driver must call
   !! `ice_update_state(.true.)` after `enter_data` before anything reads them.
   !! Outside that, `ice_update_state` is called only at the two sites the plan
   !! allows (MASTER_PLAN C.6-5): `aice0` is STATE through transport (D10 writes
   !! it), so rebuilding it elsewhere would overwrite the transport result.  A
   !! ghost aggregate holds the aggregate of its own ghost state only when
   !! `include_ghosts = .true.`; otherwise ghosts need a halo exchange before
   !! anything reads them.
   !!
   !! **Restart.**  `register_restart` adds the persistent arrays under the v2
   !! tags (`ice_aicen`, `ice_vicen`, `ice_vsnon`, `ice_Tsfcn`, and the rank-4
   !! arrays as rank-3 slices `ice_qicen_k<n>`, `ice_sicen_k<n>`,
   !! `ice_qsnon_k<n>`), all optional, plus the version scalar
   !! `ice_state_version` (= `ICE_STATE_VERSION`, 2: Icepack layout; the SIS2
   !! state is version 1 and carries `ice_part_size`).  Full local arrays
   !! (ghosts included) are written, like every other registry entry.
   !!
   !! Every entry is `optional`, so a SIS2 file (no `ice_state_version`) read
   !! under `scheme = "cice"` would silently cold-start: the PR that makes the
   !! scheme runnable MUST refuse or convert when `ice_state_version` is not
   !! found (`reg%entry_found`).  This slot adds 1 + 4 + 2*nilyr + nslyr entries
   !! (20 at nilyr = 7, 22 at the maximum nilyr = 8, version scalar
   !! included) against `MAX_RESTART_ENTRIES = 256` in `rdb_ocean_restart`.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_restart, only: restart_registry_t
   use pic_strings, only: to_string
   use rdb_ice_params, only: ice_params_t
   use rdb_ice_itd_params, only: ice_itd_bounds, ITD_BOUNDS_OK
   use rdb_ice_aggregate, only: ice_aggregate_state_impl
   use rdb_ice_seed, only: ice_seed_uniform_impl, ICE_SEED_OK, ICE_SEED_BAD_EXTENT
   implicit none
   private

   public :: ice_cice_state_t
   public :: ice_update_state, ice_state_seed_uniform
   public :: ICE_STATE_VERSION
   public :: ICE_STATE_INIT_OK, ICE_STATE_INIT_BAD_PARAMS, ICE_STATE_INIT_BAD_BOUNDS

   real(wp), parameter :: ICE_STATE_VERSION = 2.0_wp
      !! Value of the `ice_state_version` restart scalar (1 = SIS2 state).

   integer, parameter :: ICE_STATE_INIT_OK = 0
      !! `init` allocated the slot.
   integer, parameter :: ICE_STATE_INIT_BAD_PARAMS = 1
      !! `p%bad` is set (unknown constants set or a structure outside the bounds).
   integer, parameter :: ICE_STATE_INIT_BAD_BOUNDS = 2
      !! The category bounds for `p%itd` could not be built (`ITD_BOUNDS_*`).

   type :: ice_cice_state_t
      !! The CICE-port ice state slot (see the module docstring).
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data` (the arrays are on the device).
      integer :: nx = 0
         !! Total extent in x including ghosts (`nx_total`).
      integer :: ny = 0
         !! Total extent in y including ghosts (`ny_total`).
      integer :: nghost = 0
         !! Ghost width; the owned cells are `nghost+1 : nghost+nx_phys` (and y).
      integer :: nx_phys = 0
         !! Owned cells in x.
      integer :: ny_phys = 0
         !! Owned cells in y.
      integer :: i_offset_global = 0
         !! Global index offset of this tile in x (for `ice_err_reduce`).
      integer :: j_offset_global = 0
         !! Global index offset in y.
      integer :: nx_global = 0
         !! Global physical extent in x (for `ice_err_reduce`).
      integer :: ncat = 0
         !! Thickness categories.
      integer :: nilyr = 0
         !! Ice layers per category.
      integer :: nslyr = 0
         !! Snow layers per category.
      type(ice_params_t) :: p
         !! Parameter bundle the kernels read; mapped once by `enter_data`.
      real(wp) :: state_version = ICE_STATE_VERSION
         !! Host scalar written to the restart file as `ice_state_version`.

      real(wp), allocatable :: aicen(:, :, :)
         !! Category area fraction (nx,ny,ncat); no category 0.  Restart.
      real(wp), allocatable :: vicen(:, :, :)
         !! Ice volume per cell area [m] (nx,ny,ncat).  Restart.
      real(wp), allocatable :: vsnon(:, :, :)
         !! Snow volume per cell area [m] (nx,ny,ncat).  Restart.
      real(wp), allocatable :: Tsfcn(:, :, :)
         !! Surface temperature [degC] (nx,ny,ncat).  Restart.
      real(wp), allocatable :: qicen(:, :, :, :)
         !! Ice layer enthalpy [J/m3 of ice] (nx,ny,ncat,nilyr), `k = 1` TOP.  Restart.
      real(wp), allocatable :: sicen(:, :, :, :)
         !! Ice layer salinity [ppt] (nx,ny,ncat,nilyr), `k = 1` TOP.  Restart.
      real(wp), allocatable :: qsnon(:, :, :, :)
         !! Snow layer enthalpy [J/m3 of snow] (nx,ny,ncat,nslyr).  Restart.
      real(wp), allocatable :: aice0(:, :)
         !! Open-water fraction (nx,ny); rebuilt by `ice_update_state`.
      real(wp), allocatable :: aice(:, :)
         !! Total ice area fraction (nx,ny).
      real(wp), allocatable :: vice(:, :)
         !! Total ice volume per cell area [m] (nx,ny).
      real(wp), allocatable :: vsno(:, :)
         !! Total snow volume per cell area [m] (nx,ny).
      real(wp), allocatable :: aice_init(:, :)
         !! Start-of-step `aice` (nx,ny); scratch.
      real(wp), allocatable :: aicen_init(:, :, :)
         !! Start-of-step `aicen` (nx,ny,ncat); scratch.
      real(wp), allocatable :: vicen_init(:, :, :)
         !! Start-of-step `vicen` (nx,ny,ncat); scratch.
      real(wp), allocatable :: vsnon_init(:, :, :)
         !! Start-of-step `vsnon` (nx,ny,ncat); scratch.
      real(wp), allocatable :: Tsfcn_init(:, :, :)
         !! Start-of-step `Tsfcn` (nx,ny,ncat); scratch.
      real(wp), allocatable :: hin_max(:)
         !! Category bounds [m], `0:ncat`.
      integer(int32), allocatable :: ice_err(:, :)
         !! Per-cell error word (nx,ny), bits of `rdb_ice_errors`.
   contains
      procedure, non_overridable :: init => cice_state_init
      procedure, non_overridable :: destroy => cice_state_destroy
      procedure, non_overridable :: enter_data => cice_state_enter_data
      procedure, non_overridable :: exit_data => cice_state_exit_data
      procedure, non_overridable :: bytes => cice_state_bytes
      procedure, non_overridable :: register_restart => cice_state_register_restart
   end type ice_cice_state_t

contains

   subroutine cice_state_init(this, grid, p, ierr)
      !! Allocate and zero the slot for `grid` and the parameter bundle `p`
      !! (its `ncat`, `nilyr`, `nslyr`, category-bound knobs).  State arrays
      !! are zero, `aice0 = 1` (open water everywhere) and `hin_max` is built
      !! from `p%itd`; a zero state is NOT a thermodynamically valid ice state
      !! (zero enthalpy is warmer than any ice) -- seed it with
      !! `ice_state_seed_uniform` or read a restart before stepping.  `ierr` is
      !! `ICE_STATE_INIT_*`; on a non-zero value nothing is allocated.
      class(ice_cice_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      type(ice_params_t), intent(in) :: p
      integer, intent(out) :: ierr

      integer :: nx, ny, nc, nl, ns, status
      real(wp), allocatable :: hb(:)

      ierr = ICE_STATE_INIT_OK
      if (p%bad) then
         ierr = ICE_STATE_INIT_BAD_PARAMS
         return
      end if
      allocate (hb(0:p%ncat))
      call ice_itd_bounds(p, hb, status)
      if (status /= ITD_BOUNDS_OK) then
         ierr = ICE_STATE_INIT_BAD_BOUNDS
         return
      end if
      if (this%is_init) call this%destroy()

      nx = grid%nx_total
      ny = grid%ny_total
      nc = p%ncat
      nl = p%nilyr
      ns = p%nslyr
      this%nx = nx
      this%ny = ny
      this%nghost = grid%nghost
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%i_offset_global = grid%i_offset_global
      this%j_offset_global = grid%j_offset_global
      this%nx_global = grid%nx_global
      this%ncat = nc
      this%nilyr = nl
      this%nslyr = ns
      this%p = p
      this%state_version = ICE_STATE_VERSION

      allocate (this%aicen(nx, ny, nc), this%vicen(nx, ny, nc), this%vsnon(nx, ny, nc), source=0.0_wp)
      allocate (this%Tsfcn(nx, ny, nc), source=0.0_wp)
      allocate (this%qicen(nx, ny, nc, nl), this%sicen(nx, ny, nc, nl), source=0.0_wp)
      allocate (this%qsnon(nx, ny, nc, ns), source=0.0_wp)
      allocate (this%aice0(nx, ny), source=1.0_wp)
      allocate (this%aice(nx, ny), this%vice(nx, ny), this%vsno(nx, ny), source=0.0_wp)
      allocate (this%aice_init(nx, ny), source=0.0_wp)
      allocate (this%aicen_init(nx, ny, nc), this%vicen_init(nx, ny, nc), &
                this%vsnon_init(nx, ny, nc), this%Tsfcn_init(nx, ny, nc), source=0.0_wp)
      allocate (this%hin_max(0:nc), source=0.0_wp)
      this%hin_max(0:nc) = hb(0:nc)
      allocate (this%ice_err(nx, ny), source=0_int32)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine cice_state_init

   subroutine cice_state_destroy(this)
      !! Reverse of `init` (unmap first).  Clears `is_init` before releasing.
      class(ice_cice_state_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%aicen)) deallocate (this%aicen)
      if (allocated(this%vicen)) deallocate (this%vicen)
      if (allocated(this%vsnon)) deallocate (this%vsnon)
      if (allocated(this%Tsfcn)) deallocate (this%Tsfcn)
      if (allocated(this%qicen)) deallocate (this%qicen)
      if (allocated(this%sicen)) deallocate (this%sicen)
      if (allocated(this%qsnon)) deallocate (this%qsnon)
      if (allocated(this%aice0)) deallocate (this%aice0)
      if (allocated(this%aice)) deallocate (this%aice)
      if (allocated(this%vice)) deallocate (this%vice)
      if (allocated(this%vsno)) deallocate (this%vsno)
      if (allocated(this%aice_init)) deallocate (this%aice_init)
      if (allocated(this%aicen_init)) deallocate (this%aicen_init)
      if (allocated(this%vicen_init)) deallocate (this%vicen_init)
      if (allocated(this%vsnon_init)) deallocate (this%vsnon_init)
      if (allocated(this%Tsfcn_init)) deallocate (this%Tsfcn_init)
      if (allocated(this%hin_max)) deallocate (this%hin_max)
      if (allocated(this%ice_err)) deallocate (this%ice_err)
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%nilyr = 0
      this%nslyr = 0
   end subroutine cice_state_destroy

   subroutine cice_state_enter_data(this)
      !! Map the slot to the device (poly TBP delegating to a `type(...)` impl).
      class(ice_cice_state_t), intent(inout) :: this
      select type (this)
      type is (ice_cice_state_t)
         call cice_state_enter_data_impl(this)
      end select
   end subroutine cice_state_enter_data

   subroutine cice_state_enter_data_impl(this)
      !! `copyin` the persistent state, the aggregates, `hin_max`, `ice_err` and
      !! the parameter bundle; `copyin` the start-of-step copies too (host zeros, so
      !! a reader before the first `ice_save_init` sees zeros).  COMPONENT
      !! arrays only.  Restart reads and seeding land on the host first.
      type(ice_cice_state_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%p)
      !$acc enter data copyin(this%aicen, this%vicen, this%vsnon, this%Tsfcn)
      !$acc enter data copyin(this%qicen, this%sicen, this%qsnon)
      !$acc enter data copyin(this%aice0, this%aice, this%vice, this%vsno)
      !$acc enter data copyin(this%hin_max, this%ice_err)
      !$acc enter data copyin(this%aice_init, this%aicen_init, this%vicen_init)
      !$acc enter data copyin(this%vsnon_init, this%Tsfcn_init)
      this%mapped = .true.
   end subroutine cice_state_enter_data_impl

   subroutine cice_state_exit_data(this)
      !! Unmap the slot (poly TBP delegating to a `type(...)` impl).
      class(ice_cice_state_t), intent(inout) :: this
      select type (this)
      type is (ice_cice_state_t)
         call cice_state_exit_data_impl(this)
      end select
   end subroutine cice_state_exit_data

   subroutine cice_state_exit_data_impl(this)
      !! `copyout` the persistent state and the aggregates (post-run host
      !! inspection and the restart write read the host copy); `delete` the rest.
      type(ice_cice_state_t), intent(inout) :: this
      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data copyout(this%aicen, this%vicen, this%vsnon, this%Tsfcn)
      !$acc exit data copyout(this%qicen, this%sicen, this%qsnon)
      !$acc exit data copyout(this%aice0, this%aice, this%vice, this%vsno)
      !$acc exit data delete(this%hin_max, this%ice_err)
      !$acc exit data delete(this%aice_init, this%aicen_init, this%vicen_init)
      !$acc exit data delete(this%vsnon_init, this%Tsfcn_init)
      !$acc exit data delete(this%p)
      this%mapped = .false.
   end subroutine cice_state_exit_data_impl

   pure function cice_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the slot (one `arr_bytes` term per array).
      class(ice_cice_state_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = arr_bytes(this%aicen) &
               + arr_bytes(this%vicen) &
               + arr_bytes(this%vsnon) &
               + arr_bytes(this%Tsfcn) &
               + arr_bytes(this%qicen) &
               + arr_bytes(this%sicen) &
               + arr_bytes(this%qsnon) &
               + arr_bytes(this%aice0) &
               + arr_bytes(this%aice) &
               + arr_bytes(this%vice) &
               + arr_bytes(this%vsno) &
               + arr_bytes(this%aice_init) &
               + arr_bytes(this%aicen_init) &
               + arr_bytes(this%vicen_init) &
               + arr_bytes(this%vsnon_init) &
               + arr_bytes(this%Tsfcn_init) &
               + arr_bytes(this%hin_max) &
               + arr_bytes(this%ice_err)
   end function cice_state_bytes

   subroutine ice_update_state(this, include_ghosts)
      !! Rebuild `aice, vice, vsno, aice0` from the category state: ONE
      !! `do concurrent` sweep (`ice_aggregate_state_impl`, in place over the SoA
      !! arrays).  Over the OWNED cells, or the whole ghosted array when
      !! `include_ghosts` (then a ghost holds the aggregate of its own ghost
      !! state).  Thin host shim: unpacks the component arrays into the
      !! explicit-shape kernel.  On an offload build the slot must be mapped.
      type(ice_cice_state_t), intent(inout) :: this
      logical, intent(in) :: include_ghosts
      integer :: is, ie, js, je

      if (.not. this%is_init) return
      if (include_ghosts) then
         is = 1
         ie = this%nx
         js = 1
         je = this%ny
      else
         is = this%nghost + 1
         ie = this%nghost + this%nx_phys
         js = this%nghost + 1
         je = this%nghost + this%ny_phys
      end if
      call ice_aggregate_state_impl(this%nx, this%ny, this%ncat, is, ie, js, je, &
                                    this%aicen, this%vicen, this%vsnon, &
                                    this%aice, this%vice, this%vsno, this%aice0, this%ice_err)
   end subroutine ice_update_state

   subroutine ice_state_seed_uniform(this, aice_in, hi, hs, Tsfc, Tf, ierr)
      !! Seed every cell (ghosts included) with one uniform slab: concentration
      !! `aice_in`, mean ice thickness `hi` [m] and snow depth `hs` [m] over the
      !! ice, surface temperature `Tsfc` [degC], ocean freezing point `Tf` [degC]
      !! at the ice base.  The slab goes into the category whose bounds contain
      !! `hi`; profiles are documented in `rdb_ice_seed`.  `ierr` is `ICE_SEED_*`.
      !!
      !! HOST code, by design: it runs once at setup, evaluates a few column
      !! profiles and broadcasts them, so a kernel would add a launch and a
      !! device-presence requirement for nothing.  Called before `enter_data` the
      !! copyin carries the values up; called on a mapped slot it follows the host
      !! write with `!$acc update device` of the COMPONENT arrays.  The aggregates
      !! are set here on the host from the closed forms (`aice = aice_in`, ...), so
      !! no device kernel is needed to make the slot consistent.  Not `pure`: the
      !! device update is a directive.
      type(ice_cice_state_t), intent(inout) :: this
      real(wp), intent(in) :: aice_in, hi, hs, Tsfc, Tf
      integer, intent(out) :: ierr

      ierr = ICE_SEED_OK
      if (.not. this%is_init) then
         ierr = ICE_SEED_BAD_EXTENT   ! no slot
         return
      end if
      call ice_seed_uniform_impl(this%nx, this%ny, this%ncat, this%nilyr, this%nslyr, this%p, &
                                 this%hin_max, aice_in, hi, hs, Tsfc, Tf, &
                                 this%aicen, this%vicen, this%vsnon, this%Tsfcn, &
                                 this%qicen, this%sicen, this%qsnon, ierr)
      if (ierr /= ICE_SEED_OK) return
      this%aice = aice_in
      this%vice = aice_in*hi
      this%vsno = aice_in*hs
      this%aice0 = max(1.0_wp - aice_in, 0.0_wp)
      this%ice_err = 0_int32
      if (this%mapped) then
         !$acc update device(this%aicen, this%vicen, this%vsnon, this%Tsfcn)
         !$acc update device(this%qicen, this%sicen, this%qsnon)
         !$acc update device(this%aice0, this%aice, this%vice, this%vsno, this%ice_err)
      end if
   end subroutine ice_state_seed_uniform

   subroutine cice_state_register_restart(this, reg)
      !! Register the persistent arrays and the version scalar with the restart
      !! registry under the v2 tags (module docstring).  Every entry is
      !! `optional`.  Full local arrays: `ng = 0`, extents `(nx, ny)`.
      class(ice_cice_state_t), target, intent(in) :: this
      type(restart_registry_t), intent(inout) :: reg
      integer :: k

      if (.not. this%is_init) return
      call reg%register_scalar("ice_state_version", this%state_version, optional=.true.)
      call reg%register_3d("ice_aicen", this%aicen, 0, this%nx, this%ny, optional=.true.)
      call reg%register_3d("ice_vicen", this%vicen, 0, this%nx, this%ny, optional=.true.)
      call reg%register_3d("ice_vsnon", this%vsnon, 0, this%nx, this%ny, optional=.true.)
      call reg%register_3d("ice_Tsfcn", this%Tsfcn, 0, this%nx, this%ny, optional=.true.)
      do k = 1, this%nilyr
         call reg%register_3d("ice_qicen_k"//to_string(k), this%qicen(:, :, :, k), 0, &
                              this%nx, this%ny, optional=.true.)
         call reg%register_3d("ice_sicen_k"//to_string(k), this%sicen(:, :, :, k), 0, &
                              this%nx, this%ny, optional=.true.)
      end do
      do k = 1, this%nslyr
         call reg%register_3d("ice_qsnon_k"//to_string(k), this%qsnon(:, :, :, k), 0, &
                              this%nx, this%ny, optional=.true.)
      end do
   end subroutine cice_state_register_restart

end module rdb_ice_cice_state
