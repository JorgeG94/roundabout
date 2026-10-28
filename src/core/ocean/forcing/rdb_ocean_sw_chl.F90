!! Chlorophyll-dependent penetrating-shortwave optics (Manizza et al. 2005).
module rdb_ocean_sw_chl
   !! Consumer of the shared NetCDF reader that turns a monthly surface
   !! chlorophyll climatology into the per-column, per-band shortwave
   !! opacities `ocean_surface_flux_t%sw_kappa` (Manizza et al. 2005,
   !! GRL 32, L05603):
   !!
   !!   kappa_n = a_n + b_n * chl**p_n         [m-1], uniform in depth,
   !!
   !! land columns carry `sw_opacity_land`.  The band split
   !! (`sw_band_fractions`) and the N-band transmission
   !! (`sw_transmission_nb`) live in `rdb_ocean_surface_flux`, so the
   !! deposition kernel, the KPP `B_0` coupling and the EPBL SW ledger all
   !! read one definition.
   !!
   !! **Two-call contract**, as `rdb_ocean_data_forcing`: `sw_chl_configure`
   !! runs once in the configure phase BEFORE `ocean_state_enter_data`
   !! (registration opens the file and allocates the reader brackets);
   !! `sw_chl_apply` runs once per outer step right after
   !! `ocean_data_input_update_all(t)`.  The chlorophyll destination
   !! (`sf%sw_chl_field`) and the opacity array are mapped by the
   !! surface-flux slot, so this slot adds no `enter_data` term.
   !!
   !! **Time axis.**  The climatology is read `DATA_TIME_CYCLIC` with a
   !! 365-day period; model `t = 0` is taken as the start of the file's
   !! year (day 0), so the first record (mid-January) is reached at about
   !! day 14.  **Floor.**  Chlorophyll is floored at `SW_CHL_FLOOR` before
   !! the power so land / fill values (negative) and zeros are defined.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: ocean_thermo_config_t, MAX_SW_BANDS
   use rdb_ocean_data_input, only: ocean_data_input_t, &
                                   ocean_data_input_register_2d, &
                                   ocean_data_input_update_2d, DATA_OOR_ERROR
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, sw_opacity_manizza
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_sw_chl_t
   public :: sw_chl_configure
   public :: sw_chl_apply
   public :: sw_chl_fill_kappa
   public :: SW_CHL_FLOOR, SW_CHL_YEAR_SECONDS

   real(wp), parameter :: SW_CHL_FLOOR = 1.0e-4_wp
      !! Lower bound on chlorophyll (mg m-3) before `chl**p`.
   real(wp), parameter :: SW_CHL_YEAR_SECONDS = 365.0_wp*86400.0_wp
      !! Cycle period of the monthly climatology (noleap year).

   type :: ocean_sw_chl_t
      !! Registration bookkeeping + the (host) opacity coefficients.
      logical :: active = .false.
      integer :: id_chl = 0
      integer :: nb = 0
      real(wp) :: a(MAX_SW_BANDS) = 0.0_wp
      real(wp) :: b(MAX_SW_BANDS) = 0.0_wp
      real(wp) :: p(MAX_SW_BANDS) = 0.0_wp
      real(wp) :: land = 10.0_wp
   end type ocean_sw_chl_t

contains

   subroutine sw_chl_configure(cfg, reader, grid, sf, this, ierr)
      !! Seed the N-band state on `sf` and register the chlorophyll field.
      !! No-op unless `cfg%sw_scheme == "manizza05"`.  Fails loud (via
      !! `ierr` when present) on a missing file / variable, a bad
      !! `sw_pen_frac`, or a band beyond 3 without its own coefficients.
      type(ocean_thermo_config_t), intent(in) :: cfg
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_sw_chl_t), intent(out) :: this
      integer, intent(out), optional :: ierr

      integer :: n, nx, ny, ng, ierr_l

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (trim(cfg%sw_scheme) /= "manizza05") return

      if (len_trim(cfg%sw_chl_file) == 0) then
         call fail("&ocean_thermo_nml sw_scheme='manizza05' requires sw_chl_file", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%sw_pen_frac /= 0.0_wp .and. cfg%sw_pen_frac /= 1.0_wp) then
         call fail("&ocean_thermo_nml sw_scheme='manizza05' band-absorbs the WHOLE "// &
                   "shortwave: sw_pen_frac must be 0 or 1, got "// &
                   to_string(cfg%sw_pen_frac), ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      this%nb = cfg%sw_nbands
      do n = 1, this%nb
         this%a(n) = cfg%sw_opacity_coef(2*n - 1)
         this%b(n) = cfg%sw_opacity_coef(2*n)
         this%p(n) = cfg%sw_chl_power(n)
         if (this%a(n) <= 0.0_wp) then
            call fail("&ocean_thermo_nml sw_opacity_coef: band "//to_string(n)// &
                      " needs a positive clear-water opacity a_n", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end do
      this%land = cfg%sw_opacity_land

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost
      call ocean_data_input_register_2d(reader, trim(cfg%sw_chl_file), trim(cfg%sw_chl_var), &
                                        grid, nx, ny, ng + 1, ng + 1, this%id_chl, &
                                        time_mode="cyclic", cycle_period=SW_CHL_YEAR_SECONDS, &
                                        oor=DATA_OOR_ERROR, ierr=ierr_l)
      if (ierr_l /= OCEAN_STATUS_OK) then
         if (present(ierr)) then
            ierr = ierr_l
            return
         end if
         error stop "sw_chl_configure: chlorophyll registration failed"
      end if

      call sf%set_sw_chl(this%nb, cfg%sw_blue_frac)
      sf%sw_chl_field = SW_CHL_FLOOR
      this%active = .true.
      call logger%info("sw_chl: Manizza (2005) "//to_string(this%nb)//"-band shortwave, chl <- '"// &
                       trim(cfg%sw_chl_var)//"' in "//trim(cfg%sw_chl_file))
   end subroutine sw_chl_configure

   subroutine sw_chl_apply(this, reader, grid, wet_mask, sf, t)
      !! Blend the chlorophyll bracket for time `t` into
      !! `sf%sw_chl_field` and refresh every band's opacity.  Call right
      !! after `ocean_data_input_update_all(t)`.
      type(ocean_sw_chl_t), intent(in) :: this
      type(ocean_data_input_t), intent(in) :: reader
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: wet_mask(grid%nx_total, grid%ny_total)
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: t

      if (.not. this%active) return
      call ocean_data_input_update_2d(reader, this%id_chl, t, grid%nx_total, grid%ny_total, &
                                      sf%sw_chl_field)
      call sw_chl_fill_kappa(this, grid%nx_total, grid%ny_total, wet_mask, &
                             sf%sw_chl_field, sf%sw_kappa)
   end subroutine sw_chl_apply

   subroutine sw_chl_fill_kappa(this, nx, ny, wet_mask, chl, kappa)
      !! Host shim: one device kernel per band over the whole 3-D `kappa`
      !! (never a section -- a section into an explicit-shape dummy would
      !! detach from the device map).
      type(ocean_sw_chl_t), intent(in) :: this
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: chl(nx, ny)
      real(wp), intent(inout) :: kappa(nx, ny, this%nb)
      integer :: n
      do n = 1, this%nb
         call sw_kappa_band_impl(nx, ny, this%nb, n, wet_mask, chl, this%a(n), this%b(n), &
                                 this%p(n), this%land, SW_CHL_FLOOR, kappa)
      end do
   end subroutine sw_chl_fill_kappa

   pure subroutine sw_kappa_band_impl(nx, ny, nb, n, wet_mask, chl, a, b, p, land, chl_floor, kappa)
      !! `kappa(:,:,n) = a + b*max(chl, floor)**p` on wet cells,
      !! `land` elsewhere.
      integer, intent(in) :: nx, ny, nb, n
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: chl(nx, ny)
      real(wp), intent(in) :: a, b, p, land, chl_floor
      real(wp), intent(inout) :: kappa(nx, ny, nb)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         if (wet_mask(i, j) > 0.5_wp) then
            kappa(i, j, n) = sw_opacity_manizza(chl(i, j), a, b, p, chl_floor)
         else
            kappa(i, j, n) = land
         end if
      end do
   end subroutine sw_kappa_band_impl

end module rdb_ocean_sw_chl
