!! Ocean-ice flux accumulators and start-of-step snapshots (roundabout-original).
module rdb_ice_flux_state
   !! `ice_flux_t` is the slot holding the ocean-ice flux contract of MASTER_PLAN
   !! C.4: cell-level accumulators `(nx, ny)` that are ZEROED once per ice step by
   !! `ice_flux_zero` and ADDED INTO by every producer (no producer ever sets: the
   !! SIS2 overwrite-order trap is gone).  Every flux is per unit CELL area and
   !! per second; signs are ICEPACK's, taken from the upstream column physics:
   !! positive into the ocean for `fresh`, `fsalt` (salinifies) and `fhocn`;
   !! positive DOWN for `fswthru` and `fswabs`; for the atmosphere-side fields
   !! `fsens` and `flat` are positive DOWNWARD (into the surface), `flwout` is
   !! NEGATIVE (it is `-emis*sigma*Tsfc^4`, net of the reflected incoming
   !! longwave once merged), and `evap > 0` is CONDENSATION (mass gain; sublimation
   !! is `evap < 0`).  `fresh` already contains `evap`
   !! (`freshn = evapn - (rhoi*dhi + rhos*dhs)/dt`).  The docstring on each
   !! component states its unit and producers.  The merged fields are valid on the
   !! OWNED cells only (the merge never writes ghosts): consumers loop owned cells
   !! or halo-refresh first.  The four-band split of the shortwave through the ice
   !! is owned by the shortwave and coupling PRs, not by this slot.  The set is
   !! exactly the `flux%` row of C.1, less `strocnT_x/y`, which are RESTART state
   !! owned by the dynamics PR (one-step lag), and `strair_x_T/y_T`, owned by the
   !! atmosphere PR; both are separate C.1 rows, not accumulators.
   !!
   !! **Cell producers versus category producers (GPU audit amendment to C.4).**
   !! The vertical thermodynamics, shortwave and atmosphere kernels run
   !! category-parallel (`do concurrent (n, j, i)`), and different `n` must not add
   !! into one cell (a race; atomics make the sum order non-reproducible).  So the
   !! slot ALSO holds the per-category outputs `(nx, ny, ncat)`:
   !! `freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn` and the
   !! thickness bookkeeping `melttn, meltbn, congeln, snoicen`.  The category
   !! kernels WRITE them (set, per category, once per step); the one merge kernel
   !! `ice_flux_merge_categories` then ADDS `sum_{n=1..ncat} aicen_init*xn` into the
   !! twelve matching cell fields in the fixed order `n = 1..ncat` (reproducible).
   !! Cell-coupled producers (lateral melt, new ice, ridging, the cleanup zaps, the
   !! rain pass-through, and the shortwave aggregate for `fswthru`/`fswthru_vis`,
   !! which is a sum over per-category values held in the radiation slot) add
   !! directly into the cell fields.  Which later PR writes which field:
   !!
   !! | fields | written by |
   !! |---|---|
   !! | `freshn, fsaltn, fhocnn, melttn, meltbn, congeln, snoicen` | therm1 / thickness changes (D3.4) |
   !! | `evapn, fsensn, flatn, flwoutn, fswabsn` | atmosphere boundary layer + BL99 (D8.1, D3.4) |
   !! | `fresh, fsalt, fhocn` (cell adds) | new ice (D4.2), lateral melt (D4.1), zaps (D5.3), ridging (D6), rain pass-through (therm2) |
   !! | `fswthru, fswthru_vis` | shortwave aggregate (D7) |
   !! | `frazil_vol, meltl` | new ice, lateral melt |
   !! | `dardg1dt, dardg2dt, dvirdgdt, opening` | ridging driver (D6.3) |
   !!
   !! **Lifecycle.**  `init(grid, ncat)` allocates every array at the ghosted
   !! extent, zeroed.  `enter_data` maps the COMPONENT arrays (never the aggregate
   !! type) with `copyin`, so the device holds defined zeros from the start;
   !! `exit_data` deletes them (the fluxes are per-step scratch, not state, and are
   !! not restart-carried).  `bytes()` counts every array.
   !!
   !! **Shims.**  `ice_flux_zero(flux)` zeroes everything (cell fields,
   !! per-category outputs, ghosts included) in ONE launch.
   !! `ice_flux_merge_categories(flux, cice, flw, ierr)` merges over the owned
   !! cells (`flw(nx,ny)`: incoming longwave, for the `flwout` merge); it first
   !! checks that both slots agree on `nx, ny, ncat, nghost` and, on a mismatch,
   !! returns `ICE_FLUX_ERR_SHAPE` without running the kernel.
   !! `ice_save_init(cice)` snapshots `aicen, vicen, vsnon, Tsfcn, aice` into their
   !! `*_init` twins (ONE launch, ghosts included).  The `*_init` arrays are
   !! `copyin`-mapped by `ice_cice_state_t` (host zeros), so before the first
   !! `ice_save_init` they read as ZERO, never as undefined device memory (a merge
   !! then adds nothing).  The ice step driver still calls `ice_flux_zero` then
   !! `ice_save_init` as its first two actions, before any reader (therm1 merge,
   !! linear ITD, atmosphere BL, EVP).  All three kernels are flat `do concurrent` bodies in the
   !! ice core (`rdb_ice_flux`); the shims here only unpack the components.
   !!
   !! **Coexistence.**  Nothing composes this slot into `ocean_sea_ice_t` yet (a
   !! later PR adds the guarded `flux` component, as `cice` was); under the default
   !! `scheme = "sis2"` nothing here is allocated, mapped or counted.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux, only: ice_flux_zero_impl, ice_flux_merge_impl, ice_save_init_impl
   implicit none
   private

   public :: ice_flux_t, ice_flux_zero, ice_flux_merge_categories, ice_save_init
   public :: ICE_FLUX_INIT_OK, ICE_FLUX_INIT_BAD_NCAT, ICE_FLUX_OK, ICE_FLUX_ERR_SHAPE, ice_slots_agree

   integer, parameter :: ICE_FLUX_INIT_OK = 0
      !! `init` allocated the slot.
   integer, parameter :: ICE_FLUX_INIT_BAD_NCAT = 1
      !! `ncat` outside `1:ICE_NCAT_MAX`; nothing allocated.

   integer, parameter :: ICE_FLUX_OK = 0
      !! `ice_flux_merge_categories` ran.
   integer, parameter :: ICE_FLUX_ERR_SHAPE = 1
      !! The two slots disagree (or an argument has the wrong shape); nothing ran.

   type :: ice_flux_t
      !! The flux accumulator slot (see the module docstring).
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data`.
      integer :: nx = 0
         !! Total extent in x including ghosts.
      integer :: ny = 0
         !! Total extent in y including ghosts.
      integer :: nghost = 0
         !! Ghost width; the owned cells are `nghost+1 : nghost+nx_phys` (and y).
      integer :: nx_phys = 0
         !! Owned cells in x.
      integer :: ny_phys = 0
         !! Owned cells in y.
      integer :: ncat = 0
         !! Thickness categories (extent of the per-category outputs).
      real(wp), allocatable :: fresh(:, :)
         !! [kg m-2 s-1] + into the ocean.  Contains `evap` (sublimation cancels), excludes new
         !! snow.  Producers: therm1 merge, rain pass-through, new ice, lateral melt, zaps, ridging
         !! snow.
      real(wp), allocatable :: fsalt(:, :)
         !! [kg m-2 s-1] + salinifies the ocean.  Producers: therm1 merge, new ice, lateral melt,
         !! zaps.
      real(wp), allocatable :: fhocn(:, :)
         !! [W m-2] + into the ocean (never frazil).  Producers: therm1 merge, lateral melt, zaps,
         !! ridged snow.
      real(wp), allocatable :: fswthru(:, :)
         !! [W m-2] + down, through the ice.  Producers: shortwave aggregate (sum aicen*fswthrun).
      real(wp), allocatable :: fswthru_vis(:, :)
         !! [W m-2] + down; the visible part (all of it under CCSM3).  Producers: shortwave
         !! aggregate.
      real(wp), allocatable :: evap(:, :)
         !! [kg m-2 s-1] + CONDENSATION (mass gain of the ice/snow); `< 0` is sublimation.  Per cell
         !! area.  Producers: therm1 merge.
      real(wp), allocatable :: fsens(:, :)
         !! [W m-2] sensible heat, + DOWNWARD (into the surface).  Producers: therm1 merge.
      real(wp), allocatable :: flat(:, :)
         !! [W m-2] latent heat, + DOWNWARD.  Leaves with the vapour mass: not in the ice heat
         !! budget.  Producers: therm1 merge.
      real(wp), allocatable :: flwout(:, :)
         !! [W m-2] outgoing longwave, NEGATIVE (`-emis*sigma*Tsfc^4`); the merge subtracts the
         !! reflected
         !! incoming part `(1-emis)*flw`.  Producers: therm1 merge.
      real(wp), allocatable :: fswabs(:, :)
         !! [W m-2] + shortwave absorbed by snow, ice AND the part passed to the ocean (`= fswsfc +
         !! fswint + fswthru`).  Producers: therm1 merge.
      real(wp), allocatable :: frazil_vol(:, :)
         !! [m per step] new-ice volume per cell area (diagnostic).  Producers: new ice.
      real(wp), allocatable :: meltt(:, :)
         !! [m per step] top ice melt per cell area (diagnostic).  Producers: therm1 merge.
      real(wp), allocatable :: meltb(:, :)
         !! [m per step] basal ice melt per cell area (diagnostic).  Producers: therm1 merge.
      real(wp), allocatable :: meltl(:, :)
         !! [m per step] lateral ice melt per cell area (diagnostic).  Producers: lateral melt.
      real(wp), allocatable :: congel(:, :)
         !! [m per step] basal congelation growth per cell area (diagnostic).  Producers: therm1
         !! merge.
      real(wp), allocatable :: snoice(:, :)
         !! [m per step] snow-ice formation per cell area (diagnostic).  Producers: therm1 merge.
      real(wp), allocatable :: dardg1dt(:, :)
         !! [1/s] area loss to ridging, ridged part (diagnostic).  Producers: ridging.
      real(wp), allocatable :: dardg2dt(:, :)
         !! [1/s] area gain by ridging (diagnostic).  Producers: ridging.
      real(wp), allocatable :: dvirdgdt(:, :)
         !! [m/s] ice volume ridged (diagnostic).  Producers: ridging.
      real(wp), allocatable :: opening(:, :)
         !! [1/s] lead-opening rate (diagnostic).  Producers: ridging.
      real(wp), allocatable :: freshn(:, :, :)
         !! Per-category twin of `fresh` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: fsaltn(:, :, :)
         !! Per-category twin of `fsalt` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: fhocnn(:, :, :)
         !! Per-category twin of `fhocn` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: evapn(:, :, :)
         !! Per-category twin of `evap` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: fsensn(:, :, :)
         !! Per-category twin of `fsens` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: flatn(:, :, :)
         !! Per-category twin of `flat` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: flwoutn(:, :, :)
         !! Per-category twin of `flwout` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: fswabsn(:, :, :)
         !! Per-category twin of `fswabs` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: melttn(:, :, :)
         !! Per-category twin of `meltt` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: meltbn(:, :, :)
         !! Per-category twin of `meltb` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: congeln(:, :, :)
         !! Per-category twin of `congel` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
      real(wp), allocatable :: snoicen(:, :, :)
         !! Per-category twin of `snoice` (nx,ny,ncat), per unit category area; written by the
         !! category-parallel kernels, merged by `ice_flux_merge_categories`.
   contains
      procedure, non_overridable :: init => flux_init
      procedure, non_overridable :: destroy => flux_destroy
      procedure, non_overridable :: enter_data => flux_enter_data
      procedure, non_overridable :: exit_data => flux_exit_data
      procedure, non_overridable :: bytes => flux_bytes
   end type ice_flux_t

contains

   subroutine flux_init(this, grid, ncat, ierr)
      !! Allocate and zero the slot for `grid` and `ncat` categories.  `ierr` is
      !! `ICE_FLUX_INIT_*`; on a non-zero value nothing is allocated.
      class(ice_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: ncat
      integer, intent(out) :: ierr

      integer :: nx, ny, nc

      ierr = ICE_FLUX_INIT_OK
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX) then
         ierr = ICE_FLUX_INIT_BAD_NCAT
         return
      end if
      if (this%is_init) call this%destroy()
      nx = grid%nx_total
      ny = grid%ny_total
      nc = ncat
      this%nx = nx
      this%ny = ny
      this%nghost = grid%nghost
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%ncat = nc
      allocate (this%fresh(nx, ny), source=0.0_wp)
      allocate (this%fsalt(nx, ny), source=0.0_wp)
      allocate (this%fhocn(nx, ny), source=0.0_wp)
      allocate (this%fswthru(nx, ny), source=0.0_wp)
      allocate (this%fswthru_vis(nx, ny), source=0.0_wp)
      allocate (this%evap(nx, ny), source=0.0_wp)
      allocate (this%fsens(nx, ny), source=0.0_wp)
      allocate (this%flat(nx, ny), source=0.0_wp)
      allocate (this%flwout(nx, ny), source=0.0_wp)
      allocate (this%fswabs(nx, ny), source=0.0_wp)
      allocate (this%frazil_vol(nx, ny), source=0.0_wp)
      allocate (this%meltt(nx, ny), source=0.0_wp)
      allocate (this%meltb(nx, ny), source=0.0_wp)
      allocate (this%meltl(nx, ny), source=0.0_wp)
      allocate (this%congel(nx, ny), source=0.0_wp)
      allocate (this%snoice(nx, ny), source=0.0_wp)
      allocate (this%dardg1dt(nx, ny), source=0.0_wp)
      allocate (this%dardg2dt(nx, ny), source=0.0_wp)
      allocate (this%dvirdgdt(nx, ny), source=0.0_wp)
      allocate (this%opening(nx, ny), source=0.0_wp)
      allocate (this%freshn(nx, ny, nc), source=0.0_wp)
      allocate (this%fsaltn(nx, ny, nc), source=0.0_wp)
      allocate (this%fhocnn(nx, ny, nc), source=0.0_wp)
      allocate (this%evapn(nx, ny, nc), source=0.0_wp)
      allocate (this%fsensn(nx, ny, nc), source=0.0_wp)
      allocate (this%flatn(nx, ny, nc), source=0.0_wp)
      allocate (this%flwoutn(nx, ny, nc), source=0.0_wp)
      allocate (this%fswabsn(nx, ny, nc), source=0.0_wp)
      allocate (this%melttn(nx, ny, nc), source=0.0_wp)
      allocate (this%meltbn(nx, ny, nc), source=0.0_wp)
      allocate (this%congeln(nx, ny, nc), source=0.0_wp)
      allocate (this%snoicen(nx, ny, nc), source=0.0_wp)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine flux_init

   subroutine flux_destroy(this)
      !! Reverse of `init` (unmap first).
      class(ice_flux_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%fresh)) deallocate (this%fresh)
      if (allocated(this%fsalt)) deallocate (this%fsalt)
      if (allocated(this%fhocn)) deallocate (this%fhocn)
      if (allocated(this%fswthru)) deallocate (this%fswthru)
      if (allocated(this%fswthru_vis)) deallocate (this%fswthru_vis)
      if (allocated(this%evap)) deallocate (this%evap)
      if (allocated(this%fsens)) deallocate (this%fsens)
      if (allocated(this%flat)) deallocate (this%flat)
      if (allocated(this%flwout)) deallocate (this%flwout)
      if (allocated(this%fswabs)) deallocate (this%fswabs)
      if (allocated(this%frazil_vol)) deallocate (this%frazil_vol)
      if (allocated(this%meltt)) deallocate (this%meltt)
      if (allocated(this%meltb)) deallocate (this%meltb)
      if (allocated(this%meltl)) deallocate (this%meltl)
      if (allocated(this%congel)) deallocate (this%congel)
      if (allocated(this%snoice)) deallocate (this%snoice)
      if (allocated(this%dardg1dt)) deallocate (this%dardg1dt)
      if (allocated(this%dardg2dt)) deallocate (this%dardg2dt)
      if (allocated(this%dvirdgdt)) deallocate (this%dvirdgdt)
      if (allocated(this%opening)) deallocate (this%opening)
      if (allocated(this%freshn)) deallocate (this%freshn)
      if (allocated(this%fsaltn)) deallocate (this%fsaltn)
      if (allocated(this%fhocnn)) deallocate (this%fhocnn)
      if (allocated(this%evapn)) deallocate (this%evapn)
      if (allocated(this%fsensn)) deallocate (this%fsensn)
      if (allocated(this%flatn)) deallocate (this%flatn)
      if (allocated(this%flwoutn)) deallocate (this%flwoutn)
      if (allocated(this%fswabsn)) deallocate (this%fswabsn)
      if (allocated(this%melttn)) deallocate (this%melttn)
      if (allocated(this%meltbn)) deallocate (this%meltbn)
      if (allocated(this%congeln)) deallocate (this%congeln)
      if (allocated(this%snoicen)) deallocate (this%snoicen)
      this%nx = 0
      this%ny = 0
      this%nghost = 0
      this%nx_phys = 0
      this%ny_phys = 0
      this%ncat = 0
   end subroutine flux_destroy

   subroutine flux_enter_data(this)
      !! Map the slot to the device (poly TBP delegating to a `type(...)` impl).
      class(ice_flux_t), intent(inout) :: this
      select type (this)
      type is (ice_flux_t)
         call flux_enter_data_impl(this)
      end select
   end subroutine flux_enter_data

   subroutine flux_enter_data_impl(this)
      !! `copyin` every component array (host zeros).  COMPONENT arrays only.
      type(ice_flux_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%fresh, this%fsalt, this%fhocn, this%fswthru, this%fswthru_vis)
      !$acc enter data copyin(this%evap, this%fsens, this%flat, this%flwout, this%fswabs)
      !$acc enter data copyin(this%frazil_vol, this%meltt, this%meltb, this%meltl, this%congel)
      !$acc enter data copyin(this%snoice, this%dardg1dt, this%dardg2dt, this%dvirdgdt, this%opening)
      !$acc enter data copyin(this%freshn, this%fsaltn, this%fhocnn, this%evapn, this%fsensn)
      !$acc enter data copyin(this%flatn, this%flwoutn, this%fswabsn, this%melttn, this%meltbn)
      !$acc enter data copyin(this%congeln, this%snoicen)
      this%mapped = .true.
   end subroutine flux_enter_data_impl

   subroutine flux_exit_data(this)
      !! Unmap the slot (poly TBP delegating to a `type(...)` impl).
      class(ice_flux_t), intent(inout) :: this
      select type (this)
      type is (ice_flux_t)
         call flux_exit_data_impl(this)
      end select
   end subroutine flux_exit_data

   subroutine flux_exit_data_impl(this)
      !! `delete` every component array (per-step scratch, nothing to copy out;
      !! inspect with `!$acc update self(...)` before this).
      type(ice_flux_t), intent(inout) :: this
      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%fresh, this%fsalt, this%fhocn, this%fswthru, this%fswthru_vis)
      !$acc exit data delete(this%evap, this%fsens, this%flat, this%flwout, this%fswabs)
      !$acc exit data delete(this%frazil_vol, this%meltt, this%meltb, this%meltl, this%congel)
      !$acc exit data delete(this%snoice, this%dardg1dt, this%dardg2dt, this%dvirdgdt, this%opening)
      !$acc exit data delete(this%freshn, this%fsaltn, this%fhocnn, this%evapn, this%fsensn)
      !$acc exit data delete(this%flatn, this%flwoutn, this%fswabsn, this%melttn, this%meltbn)
      !$acc exit data delete(this%congeln, this%snoicen)
      this%mapped = .false.
   end subroutine flux_exit_data_impl

   pure function flux_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the slot (one `arr_bytes` term per array).
      class(ice_flux_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = 0_int64 &
               + arr_bytes(this%fresh) &
               + arr_bytes(this%fsalt) &
               + arr_bytes(this%fhocn) &
               + arr_bytes(this%fswthru) &
               + arr_bytes(this%fswthru_vis) &
               + arr_bytes(this%evap) &
               + arr_bytes(this%fsens) &
               + arr_bytes(this%flat) &
               + arr_bytes(this%flwout) &
               + arr_bytes(this%fswabs) &
               + arr_bytes(this%frazil_vol) &
               + arr_bytes(this%meltt) &
               + arr_bytes(this%meltb) &
               + arr_bytes(this%meltl) &
               + arr_bytes(this%congel) &
               + arr_bytes(this%snoice) &
               + arr_bytes(this%dardg1dt) &
               + arr_bytes(this%dardg2dt) &
               + arr_bytes(this%dvirdgdt) &
               + arr_bytes(this%opening) &
               + arr_bytes(this%freshn) &
               + arr_bytes(this%fsaltn) &
               + arr_bytes(this%fhocnn) &
               + arr_bytes(this%evapn) &
               + arr_bytes(this%fsensn) &
               + arr_bytes(this%flatn) &
               + arr_bytes(this%flwoutn) &
               + arr_bytes(this%fswabsn) &
               + arr_bytes(this%melttn) &
               + arr_bytes(this%meltbn) &
               + arr_bytes(this%congeln) &
               + arr_bytes(this%snoicen) &
               + 0_int64
   end function flux_bytes

   subroutine ice_flux_zero(flux)
      !! Zero every accumulator and every per-category output, ghosts included, in
      !! ONE launch.  Call once per ice step before any producer.  On an offload
      !! build the slot must be mapped.
      type(ice_flux_t), intent(inout) :: flux

      if (.not. flux%is_init) return
      call ice_flux_zero_impl(flux%nx, flux%ny, flux%ncat, &
                              flux%fresh, flux%fsalt, flux%fhocn, flux%fswthru, &
                              flux%fswthru_vis, flux%evap, flux%fsens, flux%flat, &
                              flux%flwout, flux%fswabs, flux%frazil_vol, flux%meltt, &
                              flux%meltb, flux%meltl, flux%congel, flux%snoice, &
                              flux%dardg1dt, flux%dardg2dt, flux%dvirdgdt, flux%opening, &
                              flux%freshn, flux%fsaltn, flux%fhocnn, flux%evapn, &
                              flux%fsensn, flux%flatn, flux%flwoutn, flux%fswabsn, &
                              flux%melttn, flux%meltbn, flux%congeln, flux%snoicen)
   end subroutine ice_flux_zero

   subroutine ice_flux_merge_categories(flux, cice, flw, ierr)
      !! Fold the per-category outputs into the cell accumulators over the OWNED
      !! cells: `x += sum_{n=1..ncat} cice%aicen_init(n) * xn(n)` for
      !! `fresh, fsalt, fhocn, evap, fsens, flat, fswabs, meltt, meltb, congel,
      !! snoice`, and `flwout += sum_n aicen_init*(flwoutn - (1 - emis)*flw)` with
      !! `emis = cice%p%thermo%emissivity` and `flw(nx,ny)` the incoming longwave
      !! [W m-2].  ADDITIVE (pre-existing cell values are kept), category order
      !! fixed, ONE launch, weights are the START-OF-STEP areas (a prior
      !! `ice_save_init`; before one they are the zeros of the host).  A non-finite
      !! contribution propagates into the cell and sets `ICE_ERR_NONFINITE` in
      !! `cice%ice_err`.  Both slots must be mapped on an offload build.
      !! `ierr = ICE_FLUX_ERR_SHAPE` (and nothing is run) if the two slots disagree
      !! on `nx, ny, ncat, nghost`, `flw` is not `(nx, ny)`, or either is not
      !! initialised.  Valid on owned cells only.
      type(ice_flux_t), intent(inout) :: flux
      type(ice_cice_state_t), intent(inout) :: cice
      real(wp), intent(in) :: flw(:, :)
      integer, intent(out) :: ierr
      integer :: is, ie, js, je

      ierr = ICE_FLUX_OK
      if (.not. ice_slots_agree(flux, cice) .or. size(flw, 1) /= flux%nx .or. size(flw, 2) /= flux%ny) then
         ierr = ICE_FLUX_ERR_SHAPE
         return
      end if
      is = flux%nghost + 1
      ie = flux%nghost + flux%nx_phys
      js = flux%nghost + 1
      je = flux%nghost + flux%ny_phys
      call ice_flux_merge_impl(flux%nx, flux%ny, flux%ncat, is, ie, js, je, &
                               cice%p%thermo%emissivity, cice%p%c%puny, flw, cice%aicen_init, &
                               flux%freshn, flux%fsaltn, flux%fhocnn, &
                               flux%evapn, flux%fsensn, flux%flatn, &
                               flux%flwoutn, flux%fswabsn, flux%melttn, &
                               flux%meltbn, flux%congeln, flux%snoicen, &
                               flux%fresh, flux%fsalt, flux%fhocn, &
                               flux%evap, flux%fsens, flux%flat, &
                               flux%flwout, flux%fswabs, flux%meltt, &
                               flux%meltb, flux%congel, flux%snoice, cice%ice_err)
   end subroutine ice_flux_merge_categories

   pure function ice_slots_agree(flux, cice) result(ok)
      !! True when both slots are initialised and agree on `nx, ny, ncat, nghost`
      !! and the owned extents.
      type(ice_flux_t), intent(in) :: flux
      type(ice_cice_state_t), intent(in) :: cice
      logical :: ok

      ok = flux%is_init .and. cice%is_init
      if (.not. ok) return
      ok = flux%nx == cice%nx .and. flux%ny == cice%ny .and. flux%ncat == cice%ncat .and. &
           flux%nghost == cice%nghost .and. flux%nx_phys == cice%nx_phys .and. flux%ny_phys == cice%ny_phys
   end function ice_slots_agree

   subroutine ice_save_init(cice)
      !! Snapshot `aicen, vicen, vsnon, Tsfcn, aice` into `aicen_init, vicen_init,
      !! vsnon_init, Tsfcn_init, aice_init` (ghosts included, ONE launch).  First
      !! action of an ice step, before any reader of the `*_init` arrays.  On an
      !! offload build the slot must be mapped.
      type(ice_cice_state_t), intent(inout) :: cice

      if (.not. cice%is_init) return
      call ice_save_init_impl(cice%nx, cice%ny, cice%ncat, &
                              cice%aicen, cice%vicen, cice%vsnon, cice%Tsfcn, cice%aice, &
                              cice%aicen_init, cice%vicen_init, cice%vsnon_init, cice%Tsfcn_init, &
                              cice%aice_init)
   end subroutine ice_save_init

end module rdb_ice_flux_state
