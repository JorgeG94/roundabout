!! JRA55-do-style prescribed atmospheric-state provider (OM3 wave-3, C2c-2b).
module rdb_atm_state
   !! Builds the nine-field atmospheric state `tas, huss, uas, vas, psl,
   !! rlds, rsds, prra, prsn` (Tsujino et al. 2018, JRA55-do) on the model
   !! T-grid: time-blend each field's current bracket on the shared PR-14
   !! reader into a reusable longitude-circle x latitude-band WINDOW
   !! buffer, then regrid that window onto the physical interior with the
   !! in-model lat-lon weights (`rdb_regrid_latlon`, Jones 1999).  No MPI
   !! communication in the provider itself — every rank reads its own
   !! hyperslab (the reader's `win_i0/win_j0/win_nx/win_ny` override, C2c-1)
   !! and regrids independently.
   !!
   !! **`rsds` -- held interval + diurnal cosine-zenith factor
   !! (`sw_interp = "coszen"`, the namelist default, C2c-2b).**  The JRA55-
   !! do `rsds` file carries 3-hour interval MEANS; registering it with
   !! `hold = .true.` and `t_offset = sw_offset_s` (default -5400 s, half
   !! the 3-hour bucket) makes the shared reader return that bucket's
   !! value PIECEWISE-CONSTANT (`rsds_rec`, weight forced to 0 -- no
   !! blend toward the next record) rather than linearly interpolated.
   !! The diurnal shape within the bucket is restored by
   !! `rsds = rsds_rec * cosz / cosz_mean` (0 where `cosz_mean < 1e-6`,
   !! i.e. the whole bucket is night): `cosz` is the INSTANTANEOUS cosine
   !! of the solar zenith angle at the query time `t`; `cosz_mean` is a
   !! 12-point midpoint quadrature of `max(0, cos Z)` over the bucket
   !! `[t_rec, t_rec + 3h)` that `rsds_rec` itself represents (`t_rec =
   !! floor(t / 3h) * 3h`, the JRA55-do bucket boundary) -- so
   !! `cosz_mean * 3h` is (to the quadrature's own accuracy) the bucket's
   !! true clear-sky-insolation-weighted mean, matching the file's
   !! interval-mean convention and preserving the daily energy total
   !! exactly in the limit where the model's own step cadence samples the
   !! same 12 nodes (`test_ocean_atm_state`'s energy-identity gate).  The
   !! solar position (declination + equation of time) comes from the
   !! Spencer (1971) Fourier fit -- a `day-of-year` + `hour-of-day`
   !! function only, computed HOST-SIDE via `rdb_calendar` (calendar
   !! arithmetic is not device-safe) and handed to the per-cell device
   !! kernel as a handful of scalars/short arrays; the kernel itself does
   !! only elementwise trigonometry on `geolonT`/`geolatT`, never touches
   !! the calendar.  `orb_year` is the reserved Berger (1978) orbital-
   !! parameter follow-up -- `validate_config` refuses any value other
   !! than its default, since this chunk's declination uses the FIXED
   !! Spencer coefficients.  `sw_interp = "linear"` disables all of the
   !! above (`hold = .false.`, `t_offset = 0`, plain blend, same answer as
   !! every other flux) -- `rsds` is then `rsds_rec` copied straight
   !! through.
   !!
   !! There is NO bulk-flux/ice-seam CONSUMER here (that is C2c-3) -- the
   !! engine only configures this slot and refreshes it every step;
   !! nothing reads `tas`/`huss`/... yet.
   !!
   !! **Two reusable window buffers, not one per field.**  `win_bl`/`win_cv`
   !! are each sized to ONE weight set's source footprint
   !! (`regrid_window`) and reused by every field that uses that weight
   !! set: `tas, huss, uas, vas, psl` always use the bilinear set (`wb`,
   !! `&ocean_atm_forcing_nml states_regrid` has no other legal value);
   !! `rlds, rsds, prra, prsn` use `wb` or the conservative set (`wc`)
   !! together, selected once at configure by `fluxes_regrid`.  Both
   !! weight sets are always built (bit-identical either way the knob is
   !! set) -- only the per-field BINDING differs.  Each field's call
   !! sequence is blend-into-window, regrid-into-a-physical-size scratch
   !! (`tgt_bl`/`tgt_cv`, a persistent workspace -- never a per-step
   !! local-allocate, CLAUDE.md "Persistent kernel workspaces"), then
   !! scatter that scratch into the field's own ghosted array at the
   !! `nghost` offset -- so the window/scratch buffers are safely reused
   !! read-after-write within one field's update, before the next field
   !! overwrites them.
   !!
   !! **Ghosts are untouched.**  Every field array is `(nx_total,
   !! ny_total)` but only the physical interior is written here -- no
   !! consumer in this chunk reads a ghost cell, and extrapolating into
   !! one is the data-forcing anti-pattern `rdb_ocean_data_forcing`'s
   !! docstring already warns against.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, TWO_PI, DEG2RAD
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_data_input, only: ocean_data_input_t, &
                                   ocean_data_input_register_2d_filelist, &
                                   ocean_data_input_update_2d
   use rdb_regrid_latlon, only: regrid_weights_t, regrid_bnds_to_edges, &
                                regrid_weights_bilinear, regrid_weights_conservative, &
                                regrid_apply_2d_impl
   use rdb_calendar, only: date_t, parse_date, seconds_to_date, gregorian_day_number
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP, OCEAN_STATUS_ERR_IO
   use rdb_error_ring, only: fail
   use rdb_io_netcdf, only: nc_open_read, nc_close, nc_get_varid, nc_get_var_1d, &
                            nc_get_var_2d, nc_get_dim_len
   use pic_strings, only: to_string
   implicit none
   private

   public :: atm_state_t
   public :: atm_state_configure
   public :: atm_state_update

   integer, parameter :: IDX_TAS = 1, IDX_HUSS = 2, IDX_UAS = 3, IDX_VAS = 4, IDX_PSL = 5
   integer, parameter :: IDX_RLDS = 6, IDX_PRRA = 7, IDX_PRSN = 8, IDX_RSDS = 9
   integer, parameter :: N_ATM_FIELDS = 9

   real(wp), parameter :: ATM_SW_BUCKET_S = 10800.0_wp
      !! JRA55-do shortwave accumulation bucket width (3 h; Tsujino et al.
      !! 2018) -- the `coszen` kernel's quadrature window, fixed by the
      !! dataset convention rather than derived from the reader's own
      !! bracket (see the module docstring's `rsds` paragraph).
   integer, parameter :: ATM_COSZ_NQUAD = 12
      !! Midpoint-quadrature node count for `cosz_mean`.

   type :: atm_state_t
      !! Nine-field atmospheric state, device-resident, physical cells
      !! only.  Owns its own regrid weights + window/scratch buffers --
      !! unlike `ocean_data_forcing_t`, which writes into OTHER slots'
      !! arrays and therefore owns no device-mapped state of its own.
      logical :: enable = .false.
         !! `.true.` once `atm_state_configure` has registered every
         !! field (mirrors `&ocean_atm_forcing_nml enable`).  `.false.`
         !! makes `atm_state_update` (and every lifecycle method below)
         !! an immediate no-op -- the default, bit-identical state.
      logical :: on_device = .false.
      logical :: sw_coszen = .false.
         !! `.true.` iff `&ocean_atm_forcing_nml sw_interp = "coszen"`
         !! (the default) -- `rsds` is registered `hold = .true.` with
         !! `t_offset = sw_offset_s` and `atm_state_update` runs the
         !! diurnal cosine-zenith kernel every step.  `.false.`
         !! (`sw_interp = "linear"`) reproduces the plain interval-mean
         !! blend every other field uses.
      type(date_t) :: epoch
         !! `&time_nml start_date`, parsed once at configure -- the
         !! calendar epoch `atm_state_update`'s Spencer (1971) solar
         !! position needs every step (`t` is seconds since this).
      logical :: fluxes_conservative = .false.
         !! `.true.` => `rlds/rsds/prra/prsn` regrid through `wc`
         !! (`fluxes_regrid = "conservative"`, the namelist default);
         !! `.false.` => through `wb` (`"bilinear"`).  `tas/huss/uas/
         !! vas/psl` always use `wb` (`states_regrid` has no other
         !! legal value).
      integer :: nlon = 0, nlat = 0
         !! Source grid size, read from the first `filelist_tas` file.
      integer :: nx_total = 0, ny_total = 0, nx_phys = 0, ny_phys = 0, nghost = 0
         !! Target-tile geometry, cached from `grid` at configure time.
      integer :: id(N_ATM_FIELDS) = 0
         !! Reader registration ids, indexed by the `IDX_*` constants.
      type(regrid_weights_t) :: wb
         !! Bilinear weights (states; fluxes when `fluxes_regrid =
         !! "bilinear"`).
      type(regrid_weights_t) :: wc
         !! Conservative weights (fluxes when `fluxes_regrid =
         !! "conservative"`, the default).
      real(wp), allocatable :: win_bl(:, :)
         !! `(nlon, wb%nsy)` -- the bilinear source window, reused by
         !! every field bound to `wb`.
      real(wp), allocatable :: win_cv(:, :)
         !! `(nlon, wc%nsy)` -- the conservative source window.
      real(wp), allocatable :: tgt_bl(:, :), tgt_cv(:, :)
         !! `(nx_phys, ny_phys)` persistent regrid-output scratch, one
         !! per weight set -- never locally allocated per step.
      real(wp), allocatable :: tas(:, :), huss(:, :), uas(:, :), vas(:, :), psl(:, :)
         !! Instantaneous states, blended (no hold).
      real(wp), allocatable :: rlds(:, :), prra(:, :), prsn(:, :), rsds(:, :)
         !! Interval-mean fluxes, blended. `rsds` is the FINAL consumer-
         !! facing shortwave field -- see `rsds_rec`/`cosz`/`cosz_mean`
         !! and the module docstring's `rsds` paragraph.
      real(wp), allocatable :: rsds_rec(:, :)
         !! The held (or, under `sw_interp = "linear"`, blended) bucket
         !! value straight off the reader -- `rsds = rsds_rec *
         !! cosz / cosz_mean` under `coszen`, `rsds = rsds_rec` under
         !! `linear`.
      real(wp), allocatable :: cosz(:, :)
         !! Instantaneous cosine of the solar zenith angle at the
         !! current query time `t`, clamped to `>= 0` (night is not
         !! negative insolation). `coszen` only -- stays 0 under
         !! `linear`.
      real(wp), allocatable :: cosz_mean(:, :)
         !! 12-point midpoint quadrature of `max(0, cos Z)` over the
         !! `rsds_rec` bucket `[t_rec, t_rec+3h)`. `coszen` only.
   contains
      procedure, non_overridable :: enter_data => atm_state_enter_data
      procedure, non_overridable :: exit_data => atm_state_exit_data
      procedure, non_overridable :: bytes => atm_state_bytes
      procedure, non_overridable :: destroy => atm_state_destroy
   end type atm_state_t

contains

   ! ======================================================================
   ! Lifecycle: device mapping, footprint, teardown.
   ! ======================================================================

   subroutine atm_state_enter_data(this)
      class(atm_state_t), intent(inout) :: this
      if (.not. this%enable) return
      !$omp target enter data map(to: this%tas, this%huss, this%uas, this%vas, this%psl)
      !$omp target enter data map(to: this%rlds, this%prra, this%prsn, this%rsds)
      !$omp target enter data map(to: this%rsds_rec, this%cosz, this%cosz_mean)
      !$omp target enter data map(alloc: this%win_bl, this%win_cv, this%tgt_bl, this%tgt_cv)
      call this%wb%enter_data()
      call this%wc%enter_data()
      this%on_device = .true.
   end subroutine atm_state_enter_data

   subroutine atm_state_exit_data(this)
      class(atm_state_t), intent(inout) :: this
      if (.not. this%enable) return
      call this%wb%exit_data()
      call this%wc%exit_data()
      !$omp target exit data map(delete: this%win_bl, this%win_cv, this%tgt_bl, this%tgt_cv)
      !$omp target exit data map(delete: this%tas, this%huss, this%uas, this%vas, this%psl)
      !$omp target exit data map(delete: this%rlds, this%prra, this%prsn, this%rsds)
      !$omp target exit data map(delete: this%rsds_rec, this%cosz, this%cosz_mean)
      this%on_device = .false.
   end subroutine atm_state_exit_data

   pure function atm_state_bytes(this) result(nbytes)
      !! Device-resident footprint: the nine field arrays, both window
      !! buffers, both regrid-scratch buffers, and both weight sets'
      !! own `bytes()`.
      class(atm_state_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (.not. this%enable) return
      nbytes = arr2_bytes(this%tas) + arr2_bytes(this%huss) + arr2_bytes(this%uas) + &
               arr2_bytes(this%vas) + arr2_bytes(this%psl) + arr2_bytes(this%rlds) + &
               arr2_bytes(this%prra) + arr2_bytes(this%prsn) + arr2_bytes(this%rsds) + &
               arr2_bytes(this%rsds_rec) + arr2_bytes(this%cosz) + arr2_bytes(this%cosz_mean) + &
               arr2_bytes(this%win_bl) + arr2_bytes(this%win_cv) + &
               arr2_bytes(this%tgt_bl) + arr2_bytes(this%tgt_cv) + &
               this%wb%bytes() + this%wc%bytes()
   end function atm_state_bytes

   pure function arr2_bytes(a) result(nbytes)
      real(wp), allocatable, intent(in) :: a(:, :)
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (allocated(a)) nbytes = int(size(a), int64)*8_int64
   end function arr2_bytes

   subroutine atm_state_destroy(this)
      class(atm_state_t), intent(inout) :: this
      if (this%on_device) call this%exit_data()
      call this%wb%destroy()
      call this%wc%destroy()
      if (allocated(this%win_bl)) deallocate (this%win_bl)
      if (allocated(this%win_cv)) deallocate (this%win_cv)
      if (allocated(this%tgt_bl)) deallocate (this%tgt_bl)
      if (allocated(this%tgt_cv)) deallocate (this%tgt_cv)
      if (allocated(this%tas)) deallocate (this%tas)
      if (allocated(this%huss)) deallocate (this%huss)
      if (allocated(this%uas)) deallocate (this%uas)
      if (allocated(this%vas)) deallocate (this%vas)
      if (allocated(this%psl)) deallocate (this%psl)
      if (allocated(this%rlds)) deallocate (this%rlds)
      if (allocated(this%prra)) deallocate (this%prra)
      if (allocated(this%prsn)) deallocate (this%prsn)
      if (allocated(this%rsds)) deallocate (this%rsds)
      if (allocated(this%rsds_rec)) deallocate (this%rsds_rec)
      if (allocated(this%cosz)) deallocate (this%cosz)
      if (allocated(this%cosz_mean)) deallocate (this%cosz_mean)
      this%enable = .false.
      this%sw_coszen = .false.
      this%nlon = 0
      this%nlat = 0
      this%nx_total = 0
      this%ny_total = 0
      this%nx_phys = 0
      this%ny_phys = 0
      this%nghost = 0
      this%id = 0
   end subroutine atm_state_destroy

   ! ======================================================================
   ! Configure (once, before ocean_state_enter_data)
   ! ======================================================================

   subroutine atm_state_configure(cfg, reader, grid, metrics, this, ierr)
      !! A no-op when `&ocean_atm_forcing_nml enable = .false.` (the
      !! default) -- nothing registers, `this%enable` stays `.false.`,
      !! `atm_state_update` is then a no-op too: bit-identical.
      !!
      !! Fails loud (never silently falls back) on: an unrecognised
      !! `sw_interp`/`fluxes_regrid`/`states_regrid` value, a malformed
      !! `&time_nml start_date`, a missing or empty filelist, or a
      !! `rdb_regrid_latlon` weight-build failure. `validate_config` has
      !! already refused a cartesian grid, a blank `start_date`, a
      !! missing ninth filelist, a non-default `orb_year`, and
      !! `states_regrid = "patch"` before this ever runs.
      type(config_t), intent(in) :: cfg
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(atm_state_t), intent(out) :: this
      integer, intent(out), optional :: ierr

      character(len=:), allocatable :: grid_files(:)
      character(len=16) :: sw_interp_l, states_regrid_l, fluxes_regrid_l
      integer :: ncid, varid, nlon, nlat, local_ierr
      integer :: win_j0_flux, nsy_flux
      real(wp), allocatable :: lon_c(:), lat_c(:), lon_bnds(:, :), lat_bnds(:, :)
      real(wp), allocatable :: lon_e(:), lat_e(:)
      real(wp) :: sw_offset_eff
      logical :: epoch_ok

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      this%enable = .false.
      if (.not. cfg%ocean%atm_forcing%enable) return

      sw_interp_l = trim(cfg%ocean%atm_forcing%sw_interp)
      select case (sw_interp_l)
      case ("coszen")
         this%sw_coszen = .true.
      case ("linear")
         this%sw_coszen = .false.
      case default
         call fail("atm_state_configure: unknown &ocean_atm_forcing_nml sw_interp '"// &
                   trim(sw_interp_l)//"' (expected 'coszen' or 'linear')", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end select

      call parse_date(trim(cfg%start_date), this%epoch, epoch_ok)
      if (.not. epoch_ok) then
         call fail("atm_state_configure: bad &time_nml start_date '"//trim(cfg%start_date)// &
                   "' (want 'YYYY-MM-DD' or 'YYYY-MM-DD hh:mm:ss') -- the coszen kernel's "// &
                   "Spencer (1971) solar position needs a calendar epoch", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      states_regrid_l = trim(cfg%ocean%atm_forcing%states_regrid)
      if (states_regrid_l /= "bilinear") then
         call fail("atm_state_configure: &ocean_atm_forcing_nml states_regrid='"// &
                   trim(states_regrid_l)//"' is not implemented (only 'bilinear' is)", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      fluxes_regrid_l = trim(cfg%ocean%atm_forcing%fluxes_regrid)
      select case (fluxes_regrid_l)
      case ("conservative")
         this%fluxes_conservative = .true.
      case ("bilinear")
         this%fluxes_conservative = .false.
      case default
         call fail("atm_state_configure: unknown &ocean_atm_forcing_nml fluxes_regrid '"// &
                   trim(fluxes_regrid_l)//"' (expected 'conservative' or 'bilinear')", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end select

      ! --- Source grid definition: lon/lat/lon_bnds/lat_bnds of the
      ! first file in the tas filelist (every field shares one source
      ! grid -- JRA55-do ships every variable on the same lat-lon mesh). ---
      call read_filelist(trim(cfg%ocean%atm_forcing%filelist_tas), grid_files, ierr=ierr)
      if (setup_failed(ierr)) return

      call nc_open_read(trim(grid_files(1)), ncid, ierr=ierr)
      if (setup_failed(ierr)) return

      call nc_get_dim_len(ncid, "lon", nlon, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if
      call nc_get_dim_len(ncid, "lat", nlat, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if

      allocate (lon_c(nlon), lat_c(nlat), lon_bnds(2, nlon), lat_bnds(2, nlat))
      allocate (lon_e(nlon + 1), lat_e(nlat + 1))

      call nc_get_varid(ncid, "lon", varid, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if
      call nc_get_var_1d(ncid, varid, lon_c, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if

      call nc_get_varid(ncid, "lat", varid, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if
      call nc_get_var_1d(ncid, varid, lat_c, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if

      call nc_get_varid(ncid, "lon_bnds", varid, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if
      call nc_get_var_2d(ncid, varid, lon_bnds, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if

      call nc_get_varid(ncid, "lat_bnds", varid, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if
      call nc_get_var_2d(ncid, varid, lat_bnds, ierr=ierr)
      if (setup_failed(ierr)) then
         call nc_close(ncid)
         return
      end if

      call nc_close(ncid, ierr=ierr)
      if (setup_failed(ierr)) return

      call regrid_bnds_to_edges(nlon, lon_bnds, lon_e)
      call regrid_bnds_to_edges(nlat, lat_bnds, lat_e)

      this%nlon = nlon
      this%nlat = nlat
      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%nghost = grid%nghost

      ! --- Weights: always build both sets (bit-identical either way
      ! fluxes_regrid is set -- only the per-field BINDING below differs). ---
      call regrid_weights_bilinear(this%wb, nlon, nlat, lon_c, lat_c, &
                                   this%nx_phys, this%ny_phys, this%nghost, &
                                   this%nx_total, this%ny_total, &
                                   metrics%geolonT, metrics%geolatT, local_ierr)
      if (local_ierr /= 0) then
         call fail("atm_state_configure: regrid_weights_bilinear failed (ierr="// &
                   to_string(local_ierr)//")", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      call regrid_weights_conservative(this%wc, nlon, nlat, lon_e, lat_e, &
                                       this%nx_phys, this%ny_phys, this%nghost, &
                                       this%nx_total, this%ny_total, &
                                       metrics%geolonBu, metrics%geolatBu, &
                                       metrics%geolonT, metrics%geolatT, local_ierr)
      if (local_ierr /= 0) then
         call fail("atm_state_configure: regrid_weights_conservative failed (ierr="// &
                   to_string(local_ierr)//")", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      allocate (this%win_bl(nlon, this%wb%nsy), source=0.0_wp)
      allocate (this%win_cv(nlon, this%wc%nsy), source=0.0_wp)
      allocate (this%tgt_bl(this%nx_phys, this%ny_phys), source=0.0_wp)
      allocate (this%tgt_cv(this%nx_phys, this%ny_phys), source=0.0_wp)
      allocate (this%tas(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%huss(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%uas(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%vas(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%psl(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%rlds(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%prra(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%prsn(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%rsds(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%rsds_rec(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%cosz(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%cosz_mean(this%nx_total, this%ny_total), source=0.0_wp)

      ! --- Registration: states always on wb; fluxes on wb or wc. ---
      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_tas), &
                          trim(cfg%ocean%atm_forcing%var_tas), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "limit", &
                          this%wb%j_lo, nlon, this%wb%nsy, this%id(IDX_TAS), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_huss), &
                          trim(cfg%ocean%atm_forcing%var_huss), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "limit", &
                          this%wb%j_lo, nlon, this%wb%nsy, this%id(IDX_HUSS), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_uas), &
                          trim(cfg%ocean%atm_forcing%var_uas), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "limit", &
                          this%wb%j_lo, nlon, this%wb%nsy, this%id(IDX_UAS), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_vas), &
                          trim(cfg%ocean%atm_forcing%var_vas), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "limit", &
                          this%wb%j_lo, nlon, this%wb%nsy, this%id(IDX_VAS), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_psl), &
                          trim(cfg%ocean%atm_forcing%var_psl), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "limit", &
                          this%wb%j_lo, nlon, this%wb%nsy, this%id(IDX_PSL), ierr=ierr)
      if (setup_failed(ierr)) return

      if (this%fluxes_conservative) then
         win_j0_flux = this%wc%j_lo
         nsy_flux = this%wc%nsy
      else
         win_j0_flux = this%wb%j_lo
         nsy_flux = this%wb%nsy
      end if

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_rlds), &
                          trim(cfg%ocean%atm_forcing%var_rlds), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                          win_j0_flux, nlon, nsy_flux, this%id(IDX_RLDS), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_prra), &
                          trim(cfg%ocean%atm_forcing%var_prra), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                          win_j0_flux, nlon, nsy_flux, this%id(IDX_PRRA), ierr=ierr)
      if (setup_failed(ierr)) return

      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_prsn), &
                          trim(cfg%ocean%atm_forcing%var_prsn), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                          win_j0_flux, nlon, nsy_flux, this%id(IDX_PRSN), ierr=ierr)
      if (setup_failed(ierr)) return

      sw_offset_eff = 0.0_wp
      if (this%sw_coszen) sw_offset_eff = cfg%ocean%atm_forcing%sw_offset_s
      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_rsds), &
                          trim(cfg%ocean%atm_forcing%var_rsds), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                          win_j0_flux, nlon, nsy_flux, this%id(IDX_RSDS), &
                          hold=this%sw_coszen, t_offset=sw_offset_eff, ierr=ierr)
      if (setup_failed(ierr)) return

      this%enable = .true.
   end subroutine atm_state_configure

   subroutine register_field(reader, grid, filelist_path, var, start_date, time_mode, &
                             ryf_start_date, ryf_period_days, oor_mode, win_j0, nlon, nsy, &
                             id, hold, t_offset, ierr)
      !! Read `filelist_path` into a path array, then register it on the
      !! shared reader with the given weight set's window
      !! `(win_i0=1, win_j0, win_nx=nlon, win_ny=nsy)` -- destination is
      !! the FULL window buffer (`dest_i0=dest_j0=1`). `hold`/`t_offset`
      !! (both optional, forwarded as-is -- absent stays absent) are the
      !! `rsds` coszen knobs; every other field leaves them off.
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      character(len=*), intent(in) :: filelist_path, var, start_date, time_mode
      character(len=*), intent(in) :: ryf_start_date, oor_mode
      integer, intent(in) :: win_j0, nlon, nsy
      real(wp), intent(in) :: ryf_period_days
      integer, intent(out) :: id
      logical, intent(in), optional :: hold
      real(wp), intent(in), optional :: t_offset
      integer, intent(out), optional :: ierr

      character(len=:), allocatable :: files(:)

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      id = 0

      call read_filelist(filelist_path, files, ierr=ierr)
      if (setup_failed(ierr)) return

      call ocean_data_input_register_2d_filelist(reader, files, var, grid, &
                                                 nlon, nsy, 1, 1, start_date, id, &
                                                 oor_mode=oor_mode, time_mode=time_mode, &
                                                 ryf_start_date=ryf_start_date, &
                                                 ryf_period_days=ryf_period_days, &
                                                 win_i0=1, win_j0=win_j0, win_nx=nlon, &
                                                 win_ny=nsy, hold=hold, t_offset=t_offset, &
                                                 ierr=ierr)
   end subroutine register_field

   subroutine read_filelist(path, files, ierr)
      !! Host-only: split a `tools/jra_filelist.py`-style one-path-per-
      !! line text file into a character array.  Blank lines are
      !! skipped; an empty or unreadable list fails loud.
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: files(:)
      integer, intent(out), optional :: ierr

      integer :: u, ios, n, i
      character(len=1024) :: line
      character(len=256) :: imsg
      logical :: exists

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      inquire (file=trim(path), exist=exists)
      if (.not. exists) then
         call fail("atm_state: filelist '"//trim(path)//"' not found", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      open (newunit=u, file=trim(path), status="old", action="read", iostat=ios, iomsg=imsg)
      if (ios /= 0) then
         call fail("atm_state: could not open filelist '"//trim(path)//"' ("//trim(imsg)//")", &
                   ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      n = 0
      do
         read (u, "(A)", iostat=ios, iomsg=imsg) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) n = n + 1
      end do

      if (n < 1) then
         close (u)
         call fail("atm_state: filelist '"//trim(path)//"' has no entries", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      allocate (character(len=1024) :: files(n))
      rewind (u)
      i = 0
      do
         read (u, "(A)", iostat=ios, iomsg=imsg) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) then
            i = i + 1
            files(i) = trim(line)
         end if
      end do
      close (u)
   end subroutine read_filelist

   pure function setup_failed(ierr) result(bad)
      !! `.true.` iff `ierr` is present and non-OK -- the "shall I return
      !! early" test after every ierr-threaded call above (copy of
      !! `rdb_ocean_engine`'s helper of the same name; when `ierr` is
      !! absent the callee already `error stop`ped on failure, so this
      !! is always `.false.` here).
      integer, intent(in), optional :: ierr
      logical :: bad
      bad = .false.
      if (present(ierr)) bad = (ierr /= OCEAN_STATUS_OK)
   end function setup_failed

   ! ======================================================================
   ! Per-step update (immediately after ocean_data_input_update_all)
   ! ======================================================================

   subroutine atm_state_update(this, reader, metrics, t)
      !! Blend + regrid every one of the nine fields' current bracket,
      !! then (`sw_interp = "coszen"`) restore `rsds`'s diurnal shape
      !! from the held bucket mean -- see the module docstring's `rsds`
      !! paragraph. `t` must be the same model time
      !! `ocean_data_input_update_all` was just called with -- the
      !! reader enforces this and aborts otherwise (see
      !! `ocean_data_input_update_2d`'s fail-loud freshness check).
      !! `metrics` supplies `geolonT`/`geolatT` for the coszen kernel
      !! (unused, and not read, when `sw_interp = "linear"`).
      class(atm_state_t), intent(inout) :: this
      type(ocean_data_input_t), intent(in) :: reader
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: t

      if (.not. this%enable) return

      call atm_apply_field(reader, this%id(IDX_TAS), t, &
                           this%nx_phys, this%ny_phys, this%nghost, &
                           this%nx_total, this%ny_total, &
                           this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                           this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                           this%win_bl, this%tgt_bl, this%tas)
      call atm_apply_field(reader, this%id(IDX_HUSS), t, &
                           this%nx_phys, this%ny_phys, this%nghost, &
                           this%nx_total, this%ny_total, &
                           this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                           this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                           this%win_bl, this%tgt_bl, this%huss)
      call atm_apply_field(reader, this%id(IDX_UAS), t, &
                           this%nx_phys, this%ny_phys, this%nghost, &
                           this%nx_total, this%ny_total, &
                           this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                           this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                           this%win_bl, this%tgt_bl, this%uas)
      call atm_apply_field(reader, this%id(IDX_VAS), t, &
                           this%nx_phys, this%ny_phys, this%nghost, &
                           this%nx_total, this%ny_total, &
                           this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                           this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                           this%win_bl, this%tgt_bl, this%vas)
      call atm_apply_field(reader, this%id(IDX_PSL), t, &
                           this%nx_phys, this%ny_phys, this%nghost, &
                           this%nx_total, this%ny_total, &
                           this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                           this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                           this%win_bl, this%tgt_bl, this%psl)

      if (this%fluxes_conservative) then
         call atm_apply_field(reader, this%id(IDX_RLDS), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wc%nmax, this%wc%nsx, this%wc%nsy, &
                              this%wc%idx_i, this%wc%idx_j, this%wc%w, &
                              this%win_cv, this%tgt_cv, this%rlds)
         call atm_apply_field(reader, this%id(IDX_PRRA), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wc%nmax, this%wc%nsx, this%wc%nsy, &
                              this%wc%idx_i, this%wc%idx_j, this%wc%w, &
                              this%win_cv, this%tgt_cv, this%prra)
         call atm_apply_field(reader, this%id(IDX_PRSN), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wc%nmax, this%wc%nsx, this%wc%nsy, &
                              this%wc%idx_i, this%wc%idx_j, this%wc%w, &
                              this%win_cv, this%tgt_cv, this%prsn)
         call atm_apply_field(reader, this%id(IDX_RSDS), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wc%nmax, this%wc%nsx, this%wc%nsy, &
                              this%wc%idx_i, this%wc%idx_j, this%wc%w, &
                              this%win_cv, this%tgt_cv, this%rsds_rec)
      else
         call atm_apply_field(reader, this%id(IDX_RLDS), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                              this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                              this%win_bl, this%tgt_bl, this%rlds)
         call atm_apply_field(reader, this%id(IDX_PRRA), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                              this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                              this%win_bl, this%tgt_bl, this%prra)
         call atm_apply_field(reader, this%id(IDX_PRSN), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                              this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                              this%win_bl, this%tgt_bl, this%prsn)
         call atm_apply_field(reader, this%id(IDX_RSDS), t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                              this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                              this%win_bl, this%tgt_bl, this%rsds_rec)
      end if

      if (this%sw_coszen) then
         call atm_sw_coszen_refresh(this, metrics, t)
      else
         call atm_state_copy_impl(this%nx_total, this%ny_total, this%rsds_rec, this%rsds)
      end if
   end subroutine atm_state_update

   subroutine atm_sw_coszen_refresh(this, metrics, t)
      !! Host-side Spencer (1971) solar-position scalars (calendar
      !! arithmetic, not device-safe) for the instantaneous query time
      !! `t` and `ATM_COSZ_NQUAD` midpoint-quadrature nodes spanning the
      !! `rsds_rec` bucket `[t_rec, t_rec+3h)`, then ONE device kernel
      !! (`atm_coszen_apply_impl`) doing the per-cell trigonometry.
      class(atm_state_t), intent(inout) :: this
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: t

      real(wp) :: t_rec, t_node, dt_node
      real(wp) :: hour_frac_inst, decl_inst, eqtime_inst
      real(wp) :: hour_frac_q(ATM_COSZ_NQUAD), decl_q(ATM_COSZ_NQUAD), &
                  eqtime_q(ATM_COSZ_NQUAD)
      integer :: k

      call atm_spencer_sun(this%epoch, t, hour_frac_inst, decl_inst, eqtime_inst)

      t_rec = floor(t/ATM_SW_BUCKET_S)*ATM_SW_BUCKET_S
      dt_node = ATM_SW_BUCKET_S/real(ATM_COSZ_NQUAD, wp)
      do k = 1, ATM_COSZ_NQUAD
         t_node = t_rec + (real(k, wp) - 0.5_wp)*dt_node
         call atm_spencer_sun(this%epoch, t_node, hour_frac_q(k), decl_q(k), eqtime_q(k))
      end do

      call atm_coszen_apply_impl(this%nx_total, this%ny_total, this%nghost, &
                                 this%nx_phys, this%ny_phys, ATM_COSZ_NQUAD, &
                                 hour_frac_q, decl_q, eqtime_q, &
                                 hour_frac_inst, decl_inst, eqtime_inst, &
                                 metrics%geolonT, metrics%geolatT, &
                                 this%rsds_rec, this%cosz, this%cosz_mean, this%rsds)
   end subroutine atm_sw_coszen_refresh

   pure subroutine atm_spencer_sun(epoch, t, hour_frac, decl, eqtime_min)
      !! Spencer (1971) Fourier fit: UTC hour-of-day `hour_frac`, solar
      !! declination `decl` (radians) and equation of time `eqtime_min`
      !! (minutes) at model time `t` seconds after `epoch`. HOST-ONLY --
      !! `rdb_calendar`'s date arithmetic is not `!$acc routine`-safe
      !! (module docstring); called from `atm_sw_coszen_refresh`, never
      !! from inside a `do concurrent`.
      type(date_t), intent(in) :: epoch
      real(wp), intent(in) :: t
      real(wp), intent(out) :: hour_frac, decl, eqtime_min

      type(date_t) :: d
      integer :: doy
      real(wp) :: gamma

      call seconds_to_date(t, epoch, d)
      doy = int(gregorian_day_number(d%y, d%m, d%d) - gregorian_day_number(d%y, 1, 1)) + 1
      hour_frac = real(d%h, wp) + real(d%mi, wp)/60.0_wp + d%s/3600.0_wp

      gamma = TWO_PI*(real(doy - 1, wp) + (hour_frac - 12.0_wp)/24.0_wp)/365.0_wp

      decl = 0.006918_wp - 0.399912_wp*cos(gamma) + 0.070257_wp*sin(gamma) &
             - 0.006758_wp*cos(2.0_wp*gamma) + 0.000907_wp*sin(2.0_wp*gamma) &
             - 0.002697_wp*cos(3.0_wp*gamma) + 0.001480_wp*sin(3.0_wp*gamma)

      eqtime_min = 229.18_wp*(0.000075_wp + 0.001868_wp*cos(gamma) - 0.032077_wp*sin(gamma) &
                              - 0.014615_wp*cos(2.0_wp*gamma) - 0.040849_wp*sin(2.0_wp*gamma))
   end subroutine atm_spencer_sun

   pure subroutine atm_coszen_apply_impl(nx_total, ny_total, nghost, nx_phys, ny_phys, nq, &
                                         hour_frac_q, decl_q, eqtime_q, &
                                         hour_frac_inst, decl_inst, eqtime_inst, &
                                         geolonT, geolatT, rsds_rec, cosz, cosz_mean, rsds)
      !! The ONE device kernel: per physical T-cell, the instantaneous
      !! `cos Z` (clamped >= 0) at `(hour_frac_inst, decl_inst,
      !! eqtime_inst)`, the `nq`-point midpoint-quadrature mean of
      !! `max(0, cos Z)` over the `(hour_frac_q, decl_q, eqtime_q)`
      !! nodes, and the coszen-restored `rsds`. Explicit-shape
      !! throughout; no cross-module calls -- every quantity but
      !! `geolonT`/`geolatT` is a host-precomputed scalar/short array
      !! (CLAUDE.md `do-concurrent` dummy + device-safety rules).
      integer, intent(in) :: nx_total, ny_total, nghost, nx_phys, ny_phys, nq
      real(wp), intent(in) :: hour_frac_q(nq), decl_q(nq), eqtime_q(nq)
      real(wp), intent(in) :: hour_frac_inst, decl_inst, eqtime_inst
      real(wp), intent(in) :: geolonT(nx_total, ny_total), geolatT(nx_total, ny_total)
      real(wp), intent(in) :: rsds_rec(nx_total, ny_total)
      real(wp), intent(inout) :: cosz(nx_total, ny_total), cosz_mean(nx_total, ny_total)
      real(wp), intent(inout) :: rsds(nx_total, ny_total)

      integer :: i, j, ii, jj, k
      real(wp) :: lat_r, lon_deg, tst_min, ha, cz, acc, mean_val

      do concurrent(j=1:ny_phys, i=1:nx_phys)
         ii = nghost + i
         jj = nghost + j
         lat_r = geolatT(ii, jj)*DEG2RAD
         lon_deg = geolonT(ii, jj)

         tst_min = hour_frac_inst*60.0_wp + eqtime_inst + 4.0_wp*lon_deg
         ha = (tst_min/4.0_wp - 180.0_wp)*DEG2RAD
         cz = sin(lat_r)*sin(decl_inst) + cos(lat_r)*cos(decl_inst)*cos(ha)
         cosz(ii, jj) = max(0.0_wp, cz)

         acc = 0.0_wp
         do k = 1, nq
            tst_min = hour_frac_q(k)*60.0_wp + eqtime_q(k) + 4.0_wp*lon_deg
            ha = (tst_min/4.0_wp - 180.0_wp)*DEG2RAD
            cz = sin(lat_r)*sin(decl_q(k)) + cos(lat_r)*cos(decl_q(k))*cos(ha)
            acc = acc + max(0.0_wp, cz)
         end do
         mean_val = acc/real(nq, wp)
         cosz_mean(ii, jj) = mean_val

         if (mean_val < 1.0e-6_wp) then
            rsds(ii, jj) = 0.0_wp
         else
            rsds(ii, jj) = rsds_rec(ii, jj)*cosz(ii, jj)/mean_val
         end if
      end do
   end subroutine atm_coszen_apply_impl

   pure subroutine atm_state_copy_impl(n1, n2, src, dst)
      !! `sw_interp = "linear"` path: `rsds = rsds_rec`, unchanged.
      !! Explicit-shape `do concurrent` copy (never a bare whole-array
      !! assignment -- CLAUDE.md `do-concurrent` dummy rule).
      integer, intent(in) :: n1, n2
      real(wp), intent(in) :: src(n1, n2)
      real(wp), intent(inout) :: dst(n1, n2)
      integer :: i, j
      do concurrent(j=1:n2, i=1:n1)
         dst(i, j) = src(i, j)
      end do
   end subroutine atm_state_copy_impl

   subroutine atm_apply_field(reader, id, t, nx_phys, ny_phys, nghost, nx_total, ny_total, &
                              nmax, nsx, nsy, idx_i, idx_j, w, win, tgt, dst)
      !! One field's blend-into-window -> regrid-into-scratch ->
      !! scatter-into-physical-interior sequence.  Explicit-shape
      !! throughout (CLAUDE.md `do-concurrent` dummy rule) -- `idx_i`/
      !! `idx_j`/`w` are the weight set's own arrays, passed straight
      !! through with no slicing.
      type(ocean_data_input_t), intent(in) :: reader
      integer, intent(in) :: id, nx_phys, ny_phys, nghost, nx_total, ny_total
      integer, intent(in) :: nmax, nsx, nsy
      real(wp), intent(in) :: t
      integer, intent(in) :: idx_i(nx_phys, ny_phys, nmax), idx_j(nx_phys, ny_phys, nmax)
      real(wp), intent(in) :: w(nx_phys, ny_phys, nmax)
      real(wp), intent(inout) :: win(nsx, nsy)
      real(wp), intent(inout) :: tgt(nx_phys, ny_phys)
      real(wp), intent(inout) :: dst(nx_total, ny_total)

      if (id <= 0) return
      call ocean_data_input_update_2d(reader, id, t, nsx, nsy, win)
      call regrid_apply_2d_impl(nx_phys, ny_phys, nmax, nsx, nsy, idx_i, idx_j, w, win, tgt)
      call atm_state_scatter_impl(nx_total, ny_total, nghost, nx_phys, ny_phys, tgt, dst)
   end subroutine atm_apply_field

   pure subroutine atm_state_scatter_impl(n1, n2, nghost, nx, ny, src, dst)
      !! Copy a physical-size regrid result into the ghosted field array
      !! at the `nghost` offset -- ghosts untouched, same offset-write
      !! idiom as `data_input_blend_2d_impl`.
      integer, intent(in) :: n1, n2, nghost, nx, ny
      real(wp), intent(in) :: src(nx, ny)
      real(wp), intent(inout) :: dst(n1, n2)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         dst(nghost + i, nghost + j) = src(i, j)
      end do
   end subroutine atm_state_scatter_impl

end module rdb_atm_state
