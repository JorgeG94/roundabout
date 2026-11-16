!! The slots of one `scheme = "cice"` run, owned together (roundabout-original).
module rdb_ice_cice_run
   !! `ice_cice_run_t` bundles everything the CICE-port ice step needs besides the
   !! category state itself (`ocean_sea_ice_t%cice`) and the C-grid dynamics slot
   !! (`ocean_sea_ice_t%dync`): the flux slot, the radiation slot, the step work slot
   !! (inputs, scratch and the persistent centre ice-ocean stress), the window-mean
   !! flux averager, the ocean imports, the transport workspace, the static grid
   !! rotation, and the uniform atmosphere used when no prescribed atmospheric state
   !! is configured.  It is composed into `ocean_state_t` (`icec`), initialised ONLY
   !! under `scheme = "cice"` (so under `"sis2"` nothing is allocated, mapped,
   !! registered or counted), mapped once by `enter_data` and released by `exit_data`.
   !!
   !! **Restart.**  `register_restart` adds the work slot's persistent centre stress
   !! (`ice_strocnTx/y`; the ONE owner of that field, the dynamics slot's own copy of
   !! `strocnT` is a per-step scratch the engine copies from) and the averager's mid-window
   !! state (`ice_avg_*`).  The category state registers itself
   !! (`ice_cice_state_t%register_restart`), and the dynamics slot its stresses/masks.
   !!
   !! **Uniform atmosphere.**  `uni_*` are the nine fields of `atm_state_t` as constant
   !! arrays (`&ocean_ice_atm_nml uni_*`, `&ocean_ice_nml sw_down`, `snowfall`), used by
   !! the `atm_prepare` hook when `&ocean_atm_forcing_nml enable = .false.`.  The heights
   !! of the forcing are `&ocean_bulk_flux_nml z_wind` and the fixed 10 m of the scalars.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_safe_math, only: safe_cos, safe_sin
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_step_cice, only: ice_step_work_t, ICE_STEP_INIT_OK
   use rdb_ice_flux_avg, only: ice_flux_avg_t
   use rdb_ice_import_ocean, only: ice_ocean_in_t
   use rdb_ice_transport_cice, only: transport_cice_ws_t, transport_result_t
   use rdb_ice_remap_flux, only: remap_ledger_t
   use rdb_ice_params, only: ICE_ADVECT_REMAP
   use rdb_ice_ledger_terms, only: ice_term_ledger_t
   implicit none
   private

   public :: ice_cice_run_t
   public :: ICE_LCUM_N, ICE_LCUM_NCLIP, ICE_LCUM_AREA, ICE_LCUM_VICE, ICE_LCUM_VSNO, ICE_LCUM_EICE, ICE_LCUM_ESNO

   integer, parameter :: ICE_LCUM_NCLIP = 1, ICE_LCUM_AREA = 2, ICE_LCUM_VICE = 3, ICE_LCUM_VSNO = 4
   integer, parameter :: ICE_LCUM_EICE = 5, ICE_LCUM_ESNO = 6, ICE_LCUM_N = 6
      !! Slots of `ice_cice_run_t%ledger_cum`.

   type :: ice_cice_run_t
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      real(wp) :: hfreeze = -1.0_wp
         !! `&ocean_ice_nml hfreeze` [m]: depth of the melt-potential provider.
      real(wp) :: z_wind = 10.0_wp
         !! Height of the wind and state forcing [m] (`&ocean_bulk_flux_nml z_wind`).
      logical :: resume_cold = .false.
         !! The restart read found NO ice state (an ocean-only checkpoint): the ice cold-starts from
         !! `&ocean_ice_ic_nml` (the engine applies it) and says so in the log.
      integer :: n_steps = 0
         !! Ice steps taken (host counter, diagnostics only).
      type(ice_flux_t) :: flux
      type(ice_rad_state_t) :: rad
      type(ice_step_work_t) :: work
      type(ice_flux_avg_t) :: avg
      type(ice_ocean_in_t) :: oin
      type(transport_cice_ws_t) :: tws
      type(transport_result_t) :: tres
         !! The result of the LAST transport call (host scalars, not restart state): `tres%ledger` is
         !! what the remap update's guards created or dropped that step (clipped cells and area,
         !! dropped volume and enthalpy content), for the conservation ledger to book.
      real(wp) :: ledger_cum(ICE_LCUM_N) = 0.0_wp
         !! The same, summed over the whole run on THIS RANK's owned cells with plain sums (a diagnostic: with the ledger on
         !! the GLOBAL, bitwise record is `ledger%cum_clip`, restart scalars `ice_led_clip_*`) (RESTART state: optional scalars `ice_remap_*`, so a resumed run
         !! continues the total; a checkpoint without them restarts it at 0).  Entries `ICE_LCUM_NCLIP`
         !! (clipped cells, as a real), `_AREA`, `_VICE`, `_VSNO`, `_EICE`, `_ESNO`: the fields of `remap_ledger_t`.
      type(ice_term_ledger_t) :: ledger
         !! The named-term conservation ledger (`&ocean_ice_nml ledger`; inert otherwise).
      logical :: skip_thermo = .false.
         !! DIAGNOSTIC: skip the thermodynamic stages (the dynamics-only invariance test).  Never set in a run.
      real(wp), allocatable :: cosA(:, :), sinA(:, :)
         !! Cosine and sine of the grid angle at T (static, set by `finalize`).
      real(wp), allocatable :: u_tas(:, :), u_huss(:, :), u_uas(:, :), u_vas(:, :), u_psl(:, :)
      real(wp), allocatable :: u_rlds(:, :), u_prra(:, :), u_prsn(:, :), u_rsds(:, :)
         !! The uniform atmosphere (see the module docstring).
   contains
      procedure, non_overridable :: init => run_init
      procedure, non_overridable :: finalize => run_finalize
      procedure, non_overridable :: enter_data => run_enter_data
      procedure, non_overridable :: exit_data => run_exit_data
      procedure, non_overridable :: destroy => run_destroy
      procedure, non_overridable :: bytes => run_bytes
      procedure, non_overridable :: register_restart => run_register_restart
   end type ice_cice_run_t

contains

   subroutine run_init(this, grid, cfg, cice, ierr)
      !! Allocate every slot for `grid` and the structure of `cice`; host, once.  `ierr /= 0`
      !! names the slot that refused (1 flux, 2 radiation, 3 work, 4 transport).
      class(ice_cice_run_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      type(ice_cice_state_t), intent(in) :: cice
      integer, intent(out) :: ierr
      integer :: e, nx, ny

      ierr = 0
      call this%destroy()
      nx = grid%nx_total
      ny = grid%ny_total
      this%nx = nx
      this%ny = ny
      this%hfreeze = cfg%ocean%ice%hfreeze
      this%z_wind = cfg%ocean%bulk_flux%z_wind
      call this%flux%init(grid, cice%ncat, e)
      if (e /= 0) then
         ierr = 1
         return
      end if
      call this%rad%init(nx, ny, cice%ncat, cice%nilyr, e)
      if (e /= 0) then
         ierr = 2
         return
      end if
      call this%rad%configure(cice%p%sw%ahmax, cice%p%c%snowpatch)
      call this%work%init(cice, e)
      if (e /= ICE_STEP_INIT_OK) then
         ierr = 3
         return
      end if
      call this%avg%init(grid)
      call this%ledger%init(cice, cfg%ocean%ice%ledger)
      call this%oin%init(grid)
      call this%tws%init(nx, ny, cice%ncat, cice%nilyr, cice%nslyr, .false., e, &
                         remap=(cice%p%transport%advection == ICE_ADVECT_REMAP), &
                         remap_flux_bytes_max=int(cfg%ocean%ice_dyn%transport_buffer_mib, int64)*1048576_int64, &
                         ledger=cfg%ocean%ice%ledger)
      if (e /= 0) then
         ierr = 4
         return
      end if
      allocate (this%cosA(nx, ny), source=1.0_wp)
      allocate (this%sinA(nx, ny), source=0.0_wp)
      allocate (this%u_tas(nx, ny), source=cfg%ocean%ice_atm%uni_tas)
      allocate (this%u_huss(nx, ny), source=cfg%ocean%ice_atm%uni_huss)
      allocate (this%u_uas(nx, ny), source=cfg%ocean%ice_atm%uni_uas)
      allocate (this%u_vas(nx, ny), source=cfg%ocean%ice_atm%uni_vas)
      allocate (this%u_psl(nx, ny), source=cfg%ocean%ice_atm%uni_psl)
      allocate (this%u_rlds(nx, ny), source=cfg%ocean%ice_atm%uni_rlds)
      allocate (this%u_prra(nx, ny), source=0.0_wp)
      allocate (this%u_prsn(nx, ny), source=cfg%ocean%ice%snowfall)
      allocate (this%u_rsds(nx, ny), source=cfg%ocean%ice%sw_down)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine run_init

   subroutine run_finalize(this, wet_T, angle_dx, tarea)
      !! Setup values that need the final metrics (host loops, BEFORE `enter_data`): the wet
      !! mask of the work slot (`wet_T > 0.5`) and `cos/sin` of the grid angle.
      class(ice_cice_run_t), intent(inout) :: this
      real(wp), intent(in) :: wet_T(:, :), angle_dx(:, :), tarea(:, :)
         !! `tarea`: the T-cell areas of the C-grid metrics (the ledger's weights).
      integer :: i, j

      if (.not. this%is_init) return
      call this%ledger%set_area(tarea)
      call this%work%set_wet(wet_T > 0.5_wp)
      do j = 1, this%ny
         do i = 1, this%nx
            this%cosA(i, j) = safe_cos(angle_dx(i, j))
            this%sinA(i, j) = safe_sin(angle_dx(i, j))
         end do
      end do
   end subroutine run_finalize

   subroutine run_enter_data(this)
      class(ice_cice_run_t), intent(inout) :: this

      if (.not. this%is_init .or. this%mapped) return
      call this%flux%enter_data()
      call this%rad%enter_data()
      call this%work%enter_data()
      call this%avg%enter_data()
      call this%oin%enter_data()
      call this%tws%enter_data()
      call this%ledger%enter_data()
      !$acc enter data copyin(this%cosA, this%sinA)
      !$acc enter data copyin(this%u_tas, this%u_huss, this%u_uas, this%u_vas, this%u_psl)
      !$acc enter data copyin(this%u_rlds, this%u_prra, this%u_prsn, this%u_rsds)
      this%mapped = .true.
   end subroutine run_enter_data

   subroutine run_exit_data(this)
      class(ice_cice_run_t), intent(inout) :: this

      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%u_tas, this%u_huss, this%u_uas, this%u_vas, this%u_psl)
      !$acc exit data delete(this%u_rlds, this%u_prra, this%u_prsn, this%u_rsds)
      !$acc exit data delete(this%cosA, this%sinA)
      call this%ledger%exit_data()
      call this%tws%exit_data()
      call this%oin%exit_data()
      call this%avg%exit_data()
      call this%work%exit_data()
      call this%rad%exit_data()
      call this%flux%exit_data()
      this%mapped = .false.
   end subroutine run_exit_data

   subroutine run_destroy(this)
      class(ice_cice_run_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      call this%tws%destroy()
      call this%ledger%destroy()
      call this%oin%destroy()
      call this%avg%destroy()
      call this%work%destroy()
      call this%rad%destroy()
      call this%flux%destroy()
      if (allocated(this%cosA)) deallocate (this%cosA, this%sinA)
      if (allocated(this%u_tas)) deallocate (this%u_tas, this%u_huss, this%u_uas, this%u_vas, this%u_psl)
      if (allocated(this%u_rlds)) deallocate (this%u_rlds, this%u_prra, this%u_prsn, this%u_rsds)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
   end subroutine run_destroy

   pure function run_bytes(this) result(nbytes)
      class(ice_cice_run_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = this%flux%bytes() + this%rad%bytes() + this%work%bytes() + this%avg%bytes() &
               + this%oin%bytes() + this%tws%bytes() + this%ledger%bytes() &
               + arr_bytes(this%cosA) + arr_bytes(this%sinA) &
               + arr_bytes(this%u_tas) + arr_bytes(this%u_huss) + arr_bytes(this%u_uas) &
               + arr_bytes(this%u_vas) + arr_bytes(this%u_psl) + arr_bytes(this%u_rlds) &
               + arr_bytes(this%u_prra) + arr_bytes(this%u_prsn) + arr_bytes(this%u_rsds)
   end function run_bytes

   subroutine run_register_restart(this, reg)
      !! The persistent centre stress and the mid-window averager (all `optional`; the caller
      !! checks the window counters of a v2 file after the read).
      class(ice_cice_run_t), intent(in), target :: this
      type(restart_registry_t), intent(inout) :: reg

      if (.not. this%is_init) return
      call this%work%register_restart(reg)
      call this%avg%register_restart(reg)
      call this%ledger%register_restart(reg)
      call reg%register_scalar("ice_remap_n_clip", this%ledger_cum(ICE_LCUM_NCLIP), optional=.true.)
      call reg%register_scalar("ice_remap_clip_area", this%ledger_cum(ICE_LCUM_AREA), optional=.true.)
      call reg%register_scalar("ice_remap_drop_vice", this%ledger_cum(ICE_LCUM_VICE), optional=.true.)
      call reg%register_scalar("ice_remap_drop_vsno", this%ledger_cum(ICE_LCUM_VSNO), optional=.true.)
      call reg%register_scalar("ice_remap_drop_eice", this%ledger_cum(ICE_LCUM_EICE), optional=.true.)
      call reg%register_scalar("ice_remap_drop_esno", this%ledger_cum(ICE_LCUM_ESNO), optional=.true.)
   end subroutine run_register_restart

end module rdb_ice_cice_run
