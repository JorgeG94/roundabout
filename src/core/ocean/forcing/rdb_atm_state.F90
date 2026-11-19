!! JRA55-do-style prescribed atmospheric-state provider (OM3 wave-3, C2c-2b).
module rdb_atm_state
   !! **Reproducibility.**  The solar-position sin/cos are `safe_sin`/`safe_cos`
   !! (host, `rdb_safe_math`) and `rdb_sm_sin`/`rdb_sm_cos`/`rdb_sm_atan2` (the per-cell
   !! zenith kernel), so `cosz` does not depend on a cell's vector position or tile width.
   !! Every phase is bounded (day-of-year angle in [0, 2 pi), hour angle within a
   !! few hundred degrees), so no time-growing argument reaches the library.
   !!
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
   !! do `rsds` file carries 3-hour interval MEANS stamped at the interval
   !! CENTRE (01:30, 04:30, ...), record `k` meaning the mean over the UTC
   !! bucket `[3k h, 3k h + 3 h)`.  Registering it with `hold = .true.` and
   !! `t_offset = sw_offset_s` (default `+5400` s, half the bucket) makes
   !! the shared reader return, for model time `t`, the record whose stamp
   !! is the last one `<= t + 5400 s`: exactly the record of the bucket
   !! that CONTAINS `t`, PIECEWISE-CONSTANT (`rsds_rec`, blend weight
   !! forced to 0).  (The offset is added to the QUERY time; a negative
   !! one selects the PREVIOUS bucket's record -- it did until 2026-11,
   !! and divided it by the current bucket's cosine mean, which blows up
   !! at every sunset/sunrise bucket boundary.)
   !!
   !! The diurnal shape within the bucket is restored by
   !! `rsds = rsds_rec * cosz / cosz_mean`, where, with the Sun's
   !! declination and the equation of time frozen at their value at the
   !! bucket midpoint (Spencer 1971 Fourier fit, host-side `rdb_calendar`
   !! arithmetic) and `omega` the hour angle, `cos Z(omega) = a + b cos
   !! omega` with `a = sin(lat) sin(dec)` and `b = cos(lat) cos(dec)`:
   !!   * `cosz_mean` is the EXACT mean of `max(0, cos Z)` over the UTC
   !!     bucket `[t_rec, t_rec + 3 h)` that `rsds_rec` represents;
   !!   * `cosz` is the EXACT mean of `max(0, cos Z)` over the model step
   !!     `[t, t + dt)` that the flux is applied for (not the instantaneous
   !!     value at `t`: a step that straddles sunset must carry the
   !!     step's mean insolation, not the value at its first instant).
   !! Both are differences of ONE closed-form antiderivative of
   !! `max(0, a + b cos omega)` (clipped at the sunrise/sunset hour angle
   !! `acos(-a/b)`; polar day and polar night are the unclipped / empty
   !! cases), so they telescope: for steps that tile the bucket
   !! (`dt` divides 3 h, `t` on the `dt` grid) `Sum_steps rsds * dt =
   !! rsds_rec * 3 h` to round-off in EVERY cell, however thin the sunlit
   !! sliver; and since each step's integral cannot exceed the bucket's,
   !! `rsds <= rsds_rec * 3 h / dt` (6 x the record at `dt = 1800`),
   !! exactly 0 when the Sun is down for the whole step, and 0 in a
   !! bucket whose mean cosine is below `1e-9` (a sliver of about 1 s;
   !! the record there is ~1e-6 W/m2 and the 0/0 would otherwise
   !! be noise).
   !!
   !! **The window is the interval the flux is APPLIED over.**  With
   !! `dt_therm_ratio = N > 1` the surface flux sampled on a thermo step is
   !! applied over `N*dt`, not `dt`: `diabatic_first` samples at the interval
   !! start and applies over `[t, t + N dt)`; `diabatic_first = .false.`
   !! samples on the interval's last step and applies over `[t - (N-1) dt,
   !! t + dt)`.  `atm_sw_window` is the ONE place that says so (the engine
   !! and the tests both call it) and `atm_state_update` integrates over
   !! exactly that window.  A window needs neither to divide the bucket nor to
   !! start on any grid of the UTC clock (`start_date` time of day included).
   !! One that lies inside a single bucket uses the formula above unchanged.
   !! One that STRADDLES a bucket edge (7200 s windows against 10800 s buckets:
   !! every third window) is split at the edge into pieces `p`, piece `p` lying
   !! in bucket `B_p` with record `rec_p`, and
   !!   `rsds = (bucket / wlen) * Sum_p rec_p * I_p / I_{B_p}`,
   !! `I_p` and `I_{B_p}` the integrals of `max(0, cos Z)` over the piece and
   !! over its whole bucket (each bucket's declination / equation of time frozen
   !! at ITS midpoint, a bucket below the mean-cosine floor contributing 0).
   !! The pieces of one bucket, taken from the windows that overlap it,
   !! telescope to `I_B`, so `Sum_windows rsds * wlen = Sum_buckets rec * bucket`
   !! to round-off in every cell for ANY tiling of time by windows.  Each piece's
   !! share `I_p / I_B` is capped at 1 (round-off guard; it cannot exceed it).
   !! The record of the neighbouring bucket comes from two further registrations
   !! of the same files (`id_rsds_prev`/`id_rsds_next`, held, query shifted by
   !! `-/+ bucket`, so the repeat-year wrap and the per-rank hyperslab windows
   !! come for free); they are refreshed by the reader with every other field,
   !! but regridded (into `rsds_alt`) only on steps whose window straddles an
   !! edge.  They are pure functions of `t`: nothing extra is checkpointed.
   !! A window LONGER than the bucket is refused (configure time when `dt_fixed`
   !! is known, `atm_sw_coszen_refresh` otherwise): it would need records more
   !! than one bucket away.
   !! **Bucket length** is not a constant: `atm_sw_check_axis` derives it at
   !! configure from the `rsds` record spacing (must divide 86400 s, the
   !! records must be uniform and stamped at the interval CENTRE) and
   !! refuses a `sw_offset_s` that is not half of it (a namelist still
   !! carrying `-5400` would silently reproduce the old one-bucket-late
   !! record).  The reader is given `sw_offset_s + 0.5 s` so the record
   !! switch never rides a float-equality tie on the `dt` grid.
   !! The solar-position scalars are computed HOST-SIDE (calendar
   !! arithmetic is not device-safe); the per-cell device kernel does only
   !! elementwise trigonometry on `geolonT`/`geolatT`.  Time of day is UTC
   !! (JRA55-do's convention), taken from `start_date` + `t`, so a
   !! `start_date` that is not on a 3 h boundary is handled.  `orb_year`
   !! is the reserved Berger (1978) orbital-parameter follow-up --
   !! `validate_config` refuses any value other than its default, since
   !! this chunk's declination uses the FIXED Spencer coefficients.
   !! `sw_interp = "linear"` disables all of the above (`hold = .false.`,
   !! `t_offset = 0`, plain blend, same answer as every other flux) --
   !! `rsds` is then `rsds_rec` copied straight through.
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
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use rdb_constants, only: wp, TWO_PI, DEG2RAD
   use rdb_safe_math, only: safe_sin, safe_cos
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
   public :: atm_sw_window

   integer, parameter :: IDX_TAS = 1, IDX_HUSS = 2, IDX_UAS = 3, IDX_VAS = 4, IDX_PSL = 5
   integer, parameter :: IDX_RLDS = 6, IDX_PRRA = 7, IDX_PRSN = 8, IDX_RSDS = 9
   integer, parameter :: N_ATM_FIELDS = 9

   real(wp), parameter :: ATM_SW_TIE_S = 0.5_wp
      !! Added to `sw_offset_s` when the `rsds` record is registered: the
      !! record switches at `t + offset = stamp`, and with `offset = bucket/2`
      !! that is a float-equality tie on every bucket edge.  Any effective
      !! offset in `[bucket/2, bucket/2 + dt)` selects the same record for
      !! every `t` on the `dt` grid, so half a second moves the switch off
      !! the tie without changing a single selection.
   real(wp), parameter :: ATM_SW_TOL_S = 1.0e-3_wp
      !! Tolerance (s) for "on the grid" / "at the bucket edge" decisions about a
      !! shortwave window (the calendar carries whole seconds; this is far below
      !! any real offset).
   real(wp), parameter :: ATM_EARTH_OMEGA = TWO_PI/86400.0_wp
      !! Hour-angle rate (rad/s): 15 degrees per hour.
   real(wp), parameter :: ATM_COSZ_MEAN_MIN = 1.0e-9_wp
      !! A bucket whose exact mean `max(0, cos Z)` is below this (a sunlit
      !! sliver of about 1 s, carrying ~1e-6 W/m2 of record) gets
      !! `rsds = 0`: the ratio of two such round-off-sized integrals is
      !! noise, not a flux.  Bounds the ratio's denominator only -- the
      !! ratio itself is bounded by `3 h / dt` by construction.

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
      real(wp) :: sw_bucket = 10800.0_wp
         !! Length (s) of the `rsds` averaging interval, derived at configure
         !! from the file's record spacing (`atm_sw_check_axis`; JRA55-do: 3 h).
         !! `coszen` only.
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
      integer :: id_rsds_prev = 0, id_rsds_next = 0
         !! `coszen` only: the SAME `rsds` files registered a second and third
         !! time, held, with the query shifted by `-bucket` / `+bucket`, so at
         !! model time `t` they return the records of the buckets before / after
         !! the one holding `t`.  Used only when a window straddles a bucket
         !! edge (`rsds_alt`).
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
      real(wp), allocatable :: rsds_alt(:, :)
         !! `coszen` only: the regridded record of the NEIGHBOURING bucket
         !! (previous or next, whichever the straddling window reaches into),
         !! refilled from the reader on every step whose window straddles a
         !! bucket edge and never read otherwise.  A pure function of `t`:
         !! nothing to checkpoint, a resumed run rebuilds it on its first step.
      real(wp), allocatable :: cosz(:, :)
         !! Mean of `max(0, cos Z)` (cosine of the solar zenith angle,
         !! clamped to `>= 0`: night is not negative insolation) over the
         !! model step `[t, t + dt)` the flux is applied for. `coszen`
         !! only -- stays 0 under `linear`.
      real(wp), allocatable :: cosz_mean(:, :)
         !! Exact mean of `max(0, cos Z)` over the `rsds_rec` bucket
         !! `[t_rec, t_rec+3h)` (UTC 3-hour bucket containing `t`); for a
         !! straddling window, over the bucket holding the window START.
         !! `coszen` only.
   contains
      procedure, non_overridable :: enter_data => atm_state_enter_data
      procedure, non_overridable :: exit_data => atm_state_exit_data
      procedure, non_overridable :: bytes => atm_state_bytes
      procedure, non_overridable :: destroy => atm_state_destroy
   end type atm_state_t

   ! Position-independent sin/cos for the per-cell zenith kernel
   ! (`rdb_sm_*`): constants once here, bodies at the end of `contains`.
#include "rdb_safe_math_consts.inc"

contains

   ! ======================================================================
   ! Lifecycle: device mapping, footprint, teardown.
   ! ======================================================================

   subroutine atm_state_enter_data(this)
      class(atm_state_t), intent(inout) :: this
      if (.not. this%enable) return
      !$acc enter data copyin(this%tas, this%huss, this%uas, this%vas, this%psl)
      !$acc enter data copyin(this%rlds, this%prra, this%prsn, this%rsds)
      !$acc enter data copyin(this%rsds_rec, this%cosz, this%cosz_mean)
      !$acc enter data copyin(this%rsds_alt)
      !$acc enter data create(this%win_bl, this%win_cv, this%tgt_bl, this%tgt_cv)
      call this%wb%enter_data()
      call this%wc%enter_data()
      this%on_device = .true.
   end subroutine atm_state_enter_data

   subroutine atm_state_exit_data(this)
      class(atm_state_t), intent(inout) :: this
      if (.not. this%enable) return
      call this%wb%exit_data()
      call this%wc%exit_data()
      !$acc exit data delete(this%win_bl, this%win_cv, this%tgt_bl, this%tgt_cv)
      !$acc exit data delete(this%tas, this%huss, this%uas, this%vas, this%psl)
      !$acc exit data delete(this%rlds, this%prra, this%prsn, this%rsds)
      !$acc exit data delete(this%rsds_rec, this%cosz, this%cosz_mean)
      !$acc exit data delete(this%rsds_alt)
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
               arr2_bytes(this%rsds_alt) + &
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
      if (allocated(this%rsds_alt)) deallocate (this%rsds_alt)
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
      this%id_rsds_prev = 0
      this%id_rsds_next = 0
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
      real(wp) :: sw_offset_eff, wlen_cfg
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
      allocate (this%rsds_alt(this%nx_total, this%ny_total), source=0.0_wp)
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
      if (this%sw_coszen) sw_offset_eff = cfg%ocean%atm_forcing%sw_offset_s + ATM_SW_TIE_S
      call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_rsds), &
                          trim(cfg%ocean%atm_forcing%var_rsds), trim(cfg%start_date), &
                          trim(cfg%ocean%atm_forcing%time_mode), &
                          trim(cfg%ocean%atm_forcing%ryf_start_date), &
                          cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                          win_j0_flux, nlon, nsy_flux, this%id(IDX_RSDS), &
                          hold=this%sw_coszen, t_offset=sw_offset_eff, ierr=ierr)
      if (setup_failed(ierr)) return
      if (this%sw_coszen) then
         call atm_sw_check_axis(this, reader, cfg%ocean%atm_forcing%sw_offset_s, ierr=ierr)
         if (setup_failed(ierr)) return
         ! The thermo window the flux is applied over, `dt * dt_therm_ratio`
         ! (`validate_config` guarantees `dt_fixed > 0`).  It may straddle one
         ! bucket edge (the neighbour record is read for it) but not span a
         ! whole bucket: a window longer than the bucket would need records
         ! two buckets away.  Refused HERE so `--validate-only` reports it.
         ! A config that never set `dt_fixed` (unit tests that drive the provider
         ! directly) skips this check; `atm_sw_coszen_refresh` repeats it per step.
         wlen_cfg = cfg%dt_fixed*real(max(cfg%ocean%vmix%dt_therm_ratio, 1), wp)
         if (cfg%dt_fixed > 0.0_wp .and. wlen_cfg > this%sw_bucket + ATM_SW_TOL_S) then
            call fail("atm_state_configure: sw_interp = 'coszen' needs dt * dt_therm_ratio <= "// &
                      "the rsds interval ("//to_string(this%sw_bucket)//" s); the shortwave "// &
                      "window is "//to_string(wlen_cfg)//" s. Lower dt_therm_ratio or use "// &
                      "sw_interp = 'linear'", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! Neighbour-bucket views of the same record series (see `id_rsds_prev`).
         call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_rsds), &
                             trim(cfg%ocean%atm_forcing%var_rsds), trim(cfg%start_date), &
                             trim(cfg%ocean%atm_forcing%time_mode), &
                             trim(cfg%ocean%atm_forcing%ryf_start_date), &
                             cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                             win_j0_flux, nlon, nsy_flux, this%id_rsds_prev, &
                             hold=.true., t_offset=sw_offset_eff - this%sw_bucket, &
                             quiet_clamp=.true., ierr=ierr)
         if (setup_failed(ierr)) return
         call register_field(reader, grid, trim(cfg%ocean%atm_forcing%filelist_rsds), &
                             trim(cfg%ocean%atm_forcing%var_rsds), trim(cfg%start_date), &
                             trim(cfg%ocean%atm_forcing%time_mode), &
                             trim(cfg%ocean%atm_forcing%ryf_start_date), &
                             cfg%ocean%atm_forcing%ryf_period_days, "extend", &
                             win_j0_flux, nlon, nsy_flux, this%id_rsds_next, &
                             hold=.true., t_offset=sw_offset_eff + this%sw_bucket, &
                             quiet_clamp=.true., ierr=ierr)
         if (setup_failed(ierr)) return
      end if

      this%enable = .true.
   end subroutine atm_state_configure

   subroutine atm_sw_check_axis(this, reader, sw_offset_s, ierr)
      !! `coszen` only: derive the shortwave averaging interval (bucket) from
      !! the `rsds` time axis (first/last record may be padding) and fail loud unless the file is what the
      !! kernel assumes -- uniform records, spacing dividing one day, stamped
      !! at the interval CENTRE on the UTC clock -- and `sw_offset_s` is half
      !! the spacing.  Sets `this%sw_bucket`.
      type(atm_state_t), intent(inout) :: this
      type(ocean_data_input_t), intent(in) :: reader
      real(wp), intent(in) :: sw_offset_s
      integer, intent(out), optional :: ierr

      real(wp) :: spacing, sod0, phase, r
      integer :: k, nt

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      nt = reader%fields(this%id(IDX_RSDS))%nt
      if (nt < 3) then
         call fail("atm_state_configure: sw_interp = 'coszen' needs >= 3 rsds records to "// &
                   "infer the averaging interval", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! The first and last record may be padding (the padded JRA55-do files add a
      ! 00:00 record in front and one after the year end): the interval is read
      ! off the INTERIOR spacing and every interior record must follow it.
      spacing = reader%fields(this%id(IDX_RSDS))%t_axis(3) - reader%fields(this%id(IDX_RSDS))%t_axis(2)
      do k = 3, nt - 1
         if (abs(reader%fields(this%id(IDX_RSDS))%t_axis(k) - &
                 reader%fields(this%id(IDX_RSDS))%t_axis(k - 1) - spacing) > 1.0_wp) then
            call fail("atm_state_configure: sw_interp = 'coszen' needs uniformly spaced rsds "// &
                      "records (non-uniform spacing at record "//to_string(k)// &
                      "); use sw_interp = 'linear' for this file", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end do
      r = modulo(86400.0_wp, spacing)
      if (spacing < 60.0_wp .or. min(r, spacing - r) > 1.0e-3_wp) then
         call fail("atm_state_configure: sw_interp = 'coszen' supports rsds record spacings that "// &
                   "divide 24 h (1 h, 3 h, 6 h ...); this file has "//to_string(spacing)//" s", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      sod0 = real(this%epoch%h, wp)*3600.0_wp + real(this%epoch%mi, wp)*60.0_wp + this%epoch%s
      phase = modulo(reader%fields(this%id(IDX_RSDS))%t_axis(2) + sod0, spacing)
      if (abs(phase - 0.5_wp*spacing) > 1.0_wp) then
         call fail("atm_state_configure: sw_interp = 'coszen' needs the rsds records stamped at "// &
                   "the interval CENTRE (JRA55-do: 01:30, 04:30, ...); second record (the first interior one) is at "// &
                   to_string(phase)//" s into its "//to_string(spacing)//" s interval", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (abs(sw_offset_s - 0.5_wp*spacing) > 1.0_wp) then
         call fail("atm_state_configure: &ocean_atm_forcing_nml sw_offset_s = "// &
                   to_string(sw_offset_s)//" is inconsistent with the rsds file: records every "// &
                   to_string(spacing)//" s stamped at the interval centre need sw_offset_s = +"// &
                   to_string(0.5_wp*spacing)//" (a negative value selects the PREVIOUS interval's record)", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      this%sw_bucket = spacing
   end subroutine atm_sw_check_axis

   pure subroutine atm_sw_window(dt, therm_ratio, diabatic_first, thermo_step, shift, length)
      !! The interval `[t + shift, t + shift + length)` over which the
      !! shortwave sampled at model time `t` is APPLIED.  A step that is not
      !! a thermo step (or `therm_ratio <= 1`): `[t, t + dt)`.  Thermo step
      !! at `therm_ratio = N > 1`: the surface flux is applied over
      !! `N*dt` -- `[t, t + N dt)` under `diabatic_first` (sampled at the
      !! interval start), `[t - (N-1) dt, t + dt)` otherwise (sampled on the
      !! interval's last step).  The engine and the tests both call this.
      real(wp), intent(in), value :: dt
      integer, intent(in), value :: therm_ratio
      logical, intent(in), value :: diabatic_first, thermo_step
      real(wp), intent(out) :: shift, length

      shift = 0.0_wp
      length = dt
      if (thermo_step .and. therm_ratio > 1) then
         length = real(therm_ratio, wp)*dt
         if (.not. diabatic_first) shift = -real(therm_ratio - 1, wp)*dt
      end if
   end subroutine atm_sw_window

   subroutine register_field(reader, grid, filelist_path, var, start_date, time_mode, &
                             ryf_start_date, ryf_period_days, oor_mode, win_j0, nlon, nsy, &
                             id, hold, t_offset, quiet_clamp, ierr)
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
      logical, intent(in), optional :: hold, quiet_clamp
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
                                                 quiet_clamp=quiet_clamp, ierr=ierr)
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

   subroutine atm_state_update(this, reader, metrics, t, dt, win_shift, win_len)
      !! Blend + regrid every one of the nine fields' current bracket,
      !! then (`sw_interp = "coszen"`) restore `rsds`'s diurnal shape
      !! from the held bucket mean -- see the module docstring's `rsds`
      !! paragraph. `t` must be the same model time
      !! `ocean_data_input_update_all` was just called with -- the
      !! reader enforces this and aborts otherwise (see
      !! `ocean_data_input_update_2d`'s fail-loud freshness check).
      !! `dt` is the model step; the coszen window is `[t, t + dt)` unless the
      !! caller passes `win_shift`/`win_len` (from `atm_sw_window`: the
      !! interval a thermo-cadence flux is really applied over).  Unused
      !! under `sw_interp = "linear"`.
      !! `metrics` supplies `geolonT`/`geolatT` for the coszen kernel
      !! (unused, and not read, when `sw_interp = "linear"`).
      class(atm_state_t), intent(inout) :: this
      type(ocean_data_input_t), intent(in) :: reader
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: t, dt
      real(wp), intent(in), optional :: win_shift, win_len
      real(wp) :: shift, wlen

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
         shift = 0.0_wp
         wlen = dt
         if (present(win_shift)) shift = win_shift
         if (present(win_len)) wlen = win_len
         call atm_sw_coszen_refresh(this, reader, metrics, t, shift, wlen)
      else
         call atm_state_copy_impl(this%nx_total, this%ny_total, this%rsds_rec, this%rsds)
      end if
   end subroutine atm_state_update

   subroutine atm_sw_coszen_refresh(this, reader, metrics, t, shift, wlen)
      !! Host-side Spencer (1971) solar-position scalars (calendar
      !! arithmetic, not device-safe) for the UTC bucket(s) holding the
      !! application window `[t + shift, t + shift + wlen)`, then ONE device
      !! kernel doing the per-cell closed-form integrals.
      !!
      !! A window inside one bucket takes the original single-bucket path
      !! (`atm_coszen_apply_impl`) with the original arithmetic, bit for bit.
      !! A window across a bucket edge (`atm_coszen_straddle_impl`) is split
      !! at the edge: the record of the bucket holding the window START and
      !! that of the bucket holding its END are paired with the integral of
      !! `max(0, cos Z)` over their own piece.  One of the two is the record
      !! the reader selected for `t` (`rsds_rec`); the other is read through
      !! the previous- or next-bucket view (`id_rsds_prev/next`) into
      !! `rsds_alt`.  Fails loud only on a window longer than the bucket (also
      !! refused at configure when `dt_fixed` is known).
      class(atm_state_t), intent(inout) :: this
      type(ocean_data_input_t), intent(in) :: reader
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: t, shift, wlen

      type(date_t) :: d
      real(wp), parameter :: TOL_S = ATM_SW_TOL_S
      real(wp) :: sod, sod_a, sod_a0, r, rb, t_off, t_rec, t_rec_e, t_cur, sod_edge
      real(wp) :: hour_frac_rec, hour_frac_mid, hour_frac_e
      real(wp) :: decl_mid, eqtime_mid, decl_e, eqtime_e, decl_unused, eqtime_unused
      logical :: cur_is_start

      if (.not. (wlen > 0.0_wp)) then
         call fail("atm_sw_coszen_refresh: non-positive shortwave window length "// &
                   to_string(wlen)//" s", code=OCEAN_STATUS_ERR_SETUP)
      end if
      if (wlen > this%sw_bucket + TOL_S) then
         call fail("atm_sw_coszen_refresh: sw_interp = 'coszen' needs the step over which the "// &
                   "flux is applied ("//to_string(wlen)//" s = dt * dt_therm_ratio) to be no "// &
                   "longer than the "//to_string(this%sw_bucket)//" s rsds interval; lower "// &
                   "dt_therm_ratio or use sw_interp = 'linear'", code=OCEAN_STATUS_ERR_SETUP)
      end if

      ! UTC time of day of the window start; the buckets are fixed to the
      ! UTC clock, not to the model epoch.
      call seconds_to_date(t, this%epoch, d)
      sod = real(d%h, wp)*3600.0_wp + real(d%mi, wp)*60.0_wp + d%s
      sod_a0 = sod + shift
      ! Snap a start that sits (within a millisecond) on the window-length grid
      ! of the clock, else on a bucket edge, so a start a hair below an edge
      ! cannot fall into the previous bucket.  Any other start is used as is:
      ! there is no grid requirement.
      sod_a = sod_a0
      r = modulo(sod_a0, wlen)
      rb = modulo(sod_a0, this%sw_bucket)
      if (min(r, wlen - r) <= TOL_S) then
         sod_a = anint(sod_a0/wlen)*wlen
      else if (min(rb, this%sw_bucket - rb) <= TOL_S) then
         sod_a = anint(sod_a0/this%sw_bucket)*this%sw_bucket
      end if
      t_off = modulo(sod_a, this%sw_bucket)
      t_rec = t + shift + (sod_a - sod_a0) - t_off

      if (t_off + wlen <= this%sw_bucket + TOL_S) then
         ! ---- one bucket: the original path, unchanged ----
         ! Hour of day at the bucket start; declination + equation of time
         ! frozen at the bucket midpoint (so every window's integral and the
         ! bucket's come from ONE function and telescope exactly).
         call atm_spencer_sun(this%epoch, t_rec, hour_frac_rec, decl_unused, eqtime_unused)
         call atm_spencer_sun(this%epoch, t_rec + 0.5_wp*this%sw_bucket, hour_frac_mid, &
                              decl_mid, eqtime_mid)

         call atm_coszen_apply_impl(this%nx_total, this%ny_total, this%nghost, &
                                    this%nx_phys, this%ny_phys, &
                                    hour_frac_rec, safe_sin(decl_mid), safe_cos(decl_mid), &
                                    eqtime_mid, t_off, wlen, this%sw_bucket, &
                                    metrics%geolonT, metrics%geolatT, &
                                    this%rsds_rec, this%cosz, this%cosz_mean, this%rsds)
         return
      end if

      ! ---- straddle: [t_off, bucket) of the start bucket + [0, t_off+wlen-bucket)
      ! of the next one.  Which of the two is the reader's record for `t`?
      ! (The reader switches records `ATM_SW_TIE_S` before the edge.)
      sod_edge = floor((sod + ATM_SW_TIE_S)/this%sw_bucket)*this%sw_bucket
      t_cur = t - (sod - sod_edge)
      t_rec_e = t_rec + this%sw_bucket
      if (abs(t_cur - t_rec) < 0.5_wp*this%sw_bucket) then
         cur_is_start = .true.
         call atm_apply_rsds_alt(this, reader, this%id_rsds_next, t)
      else if (abs(t_cur - t_rec_e) < 0.5_wp*this%sw_bucket) then
         cur_is_start = .false.
         call atm_apply_rsds_alt(this, reader, this%id_rsds_prev, t)
      else
         call fail("atm_sw_coszen_refresh: the straddling shortwave window does not touch the "// &
                   "bucket of the record selected for t (internal inconsistency)", &
                   code=OCEAN_STATUS_ERR_SETUP)
      end if

      call atm_spencer_sun(this%epoch, t_rec, hour_frac_rec, decl_unused, eqtime_unused)
      call atm_spencer_sun(this%epoch, t_rec + 0.5_wp*this%sw_bucket, hour_frac_mid, &
                           decl_mid, eqtime_mid)
      call atm_spencer_sun(this%epoch, t_rec_e, hour_frac_e, decl_unused, eqtime_unused)
      call atm_spencer_sun(this%epoch, t_rec_e + 0.5_wp*this%sw_bucket, hour_frac_mid, &
                           decl_e, eqtime_e)

      if (cur_is_start) then
         call atm_coszen_straddle_impl(this%nx_total, this%ny_total, this%nghost, &
                                       this%nx_phys, this%ny_phys, &
                                       hour_frac_rec, safe_sin(decl_mid), safe_cos(decl_mid), &
                                       eqtime_mid, hour_frac_e, safe_sin(decl_e), &
                                       safe_cos(decl_e), eqtime_e, t_off, wlen, this%sw_bucket, &
                                       metrics%geolonT, metrics%geolatT, &
                                       this%rsds_rec, this%rsds_alt, &
                                       this%cosz, this%cosz_mean, this%rsds)
      else
         call atm_coszen_straddle_impl(this%nx_total, this%ny_total, this%nghost, &
                                       this%nx_phys, this%ny_phys, &
                                       hour_frac_rec, safe_sin(decl_mid), safe_cos(decl_mid), &
                                       eqtime_mid, hour_frac_e, safe_sin(decl_e), &
                                       safe_cos(decl_e), eqtime_e, t_off, wlen, this%sw_bucket, &
                                       metrics%geolonT, metrics%geolatT, &
                                       this%rsds_alt, this%rsds_rec, &
                                       this%cosz, this%cosz_mean, this%rsds)
      end if
   end subroutine atm_sw_coszen_refresh

   subroutine atm_apply_rsds_alt(this, reader, id, t)
      !! Blend + regrid the neighbour-bucket `rsds` view `id` (previous or
      !! next record) into `rsds_alt`, through the same weight set as `rsds`.
      class(atm_state_t), intent(inout) :: this
      type(ocean_data_input_t), intent(in) :: reader
      integer, intent(in) :: id
      real(wp), intent(in) :: t

      if (this%fluxes_conservative) then
         call atm_apply_field(reader, id, t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wc%nmax, this%wc%nsx, this%wc%nsy, &
                              this%wc%idx_i, this%wc%idx_j, this%wc%w, &
                              this%win_cv, this%tgt_cv, this%rsds_alt)
      else
         call atm_apply_field(reader, id, t, &
                              this%nx_phys, this%ny_phys, this%nghost, &
                              this%nx_total, this%ny_total, &
                              this%wb%nmax, this%wb%nsx, this%wb%nsy, &
                              this%wb%idx_i, this%wb%idx_j, this%wb%w, &
                              this%win_bl, this%tgt_bl, this%rsds_alt)
      end if
   end subroutine atm_apply_rsds_alt

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

      decl = 0.006918_wp - 0.399912_wp*safe_cos(gamma) + 0.070257_wp*safe_sin(gamma) &
             - 0.006758_wp*safe_cos(2.0_wp*gamma) + 0.000907_wp*safe_sin(2.0_wp*gamma) &
             - 0.002697_wp*safe_cos(3.0_wp*gamma) + 0.001480_wp*safe_sin(3.0_wp*gamma)

      eqtime_min = 229.18_wp*(0.000075_wp + 0.001868_wp*safe_cos(gamma) - 0.032077_wp*safe_sin(gamma) &
                              - 0.014615_wp*safe_cos(2.0_wp*gamma) - 0.040849_wp*safe_sin(2.0_wp*gamma))
   end subroutine atm_spencer_sun

   pure subroutine atm_coszen_apply_impl(nx_total, ny_total, nghost, nx_phys, ny_phys, &
                                         hour_frac_rec, sdec, cdec, eqtime, t_off, dt, bucket, &
                                         geolonT, geolatT, rsds_rec, cosz, cosz_mean, rsds)
      !! The ONE device kernel.  Per physical T-cell, with `cos Z(omega) =
      !! a + b cos omega` (`a = sin(lat) sin(dec)`, `b = cos(lat) cos(dec)`,
      !! `omega` the hour angle, advancing `ATM_EARTH_OMEGA` rad/s): the
      !! exact mean of `max(0, cos Z)` over the step `[t, t+dt)` (`cosz`)
      !! and over the whole `bucket` (`cosz_mean`) -- both differences of
      !! `atm_cosz_cumulative` -- and `rsds = rsds_rec * cosz /
      !! cosz_mean` (see the module docstring for the conservation and
      !! bound this gives).  `t_off` = seconds from the bucket start to
      !! `t`; `hour_frac_rec` = UTC hour of day at the bucket start;
      !! `sdec`/`cdec`/`eqtime` frozen at the bucket midpoint.
      !! Explicit-shape throughout; no cross-module calls -- every quantity
      !! but `geolonT`/`geolatT` is a host-precomputed scalar.
      integer, intent(in), value :: nx_total, ny_total, nghost, nx_phys, ny_phys
      real(wp), intent(in), value :: hour_frac_rec, sdec, cdec, eqtime, t_off, dt, bucket
      real(wp), intent(in) :: geolonT(nx_total, ny_total), geolatT(nx_total, ny_total)
      real(wp), intent(in) :: rsds_rec(nx_total, ny_total)
      real(wp), intent(inout) :: cosz(nx_total, ny_total), cosz_mean(nx_total, ny_total)
      real(wp), intent(inout) :: rsds(nx_total, ny_total)

      integer :: i, j, ii, jj
      real(wp) :: lat_r, slat, clat, lon_deg, a, b, w_bucket0, w_step0, ratio
      real(wp) :: step_mean, bucket_mean, ratio_max

      ratio_max = bucket/max(dt, tiny(1.0_wp))
      do concurrent(j=1:ny_phys, i=1:nx_phys)
         ii = nghost + i
         jj = nghost + j
         lat_r = geolatT(ii, jj)*DEG2RAD
         lon_deg = geolonT(ii, jj)
         slat = rdb_sm_sin(lat_r)
         clat = rdb_sm_cos(lat_r)
         a = slat*sdec
         b = clat*cdec

         ! Hour angle at the bucket start, and at the step start.
         w_bucket0 = (15.0_wp*hour_frac_rec + 0.25_wp*eqtime + lon_deg - 180.0_wp)*DEG2RAD
         w_step0 = w_bucket0 + ATM_EARTH_OMEGA*t_off

         bucket_mean = (atm_cosz_cumulative(w_bucket0 + ATM_EARTH_OMEGA*bucket, a, b) &
                        - atm_cosz_cumulative(w_bucket0, a, b))/(ATM_EARTH_OMEGA*bucket)
         step_mean = (atm_cosz_cumulative(w_step0 + ATM_EARTH_OMEGA*dt, a, b) &
                      - atm_cosz_cumulative(w_step0, a, b))/(ATM_EARTH_OMEGA*dt)
         bucket_mean = max(0.0_wp, bucket_mean)
         step_mean = max(0.0_wp, step_mean)
         cosz(ii, jj) = step_mean
         cosz_mean(ii, jj) = bucket_mean

         if (bucket_mean < ATM_COSZ_MEAN_MIN) then
            rsds(ii, jj) = 0.0_wp
         else
            ratio = min(step_mean/bucket_mean, ratio_max)
            rsds(ii, jj) = rsds_rec(ii, jj)*ratio
         end if
      end do
   end subroutine atm_coszen_apply_impl

   pure subroutine atm_coszen_straddle_impl(nx_total, ny_total, nghost, nx_phys, ny_phys, &
                                            hf_s, sdec_s, cdec_s, eq_s, hf_e, sdec_e, cdec_e, &
                                            eq_e, t_off, wlen, bucket, geolonT, geolatT, &
                                            rec_s, rec_e, cosz, cosz_mean, rsds)
      !! The window `[t_off, t_off + wlen)` seconds after the START of bucket S
      !! runs past the end of S (`t_off + wlen > bucket`) into the next bucket
      !! E.  Piece 1 is `[t_off, bucket)` of S (record `rec_s`, Sun frozen at S's
      !! midpoint: `hf_s`/`sdec_s`/`cdec_s`/`eq_s`), piece 2 is `[0, t_off + wlen
      !! - bucket)` of E (`rec_e`, `hf_e`/...).  With `I_p` the integral of
      !! `max(0, cos Z)` over the piece and `I_B` over its whole bucket (the
      !! same `atm_cosz_cumulative` differences as the one-bucket kernel),
      !!   rsds = (bucket / wlen) * Sum_p rec_p * I_p / I_B_p .
      !! Each bucket then hands out exactly its own `rec * bucket` between the
      !! windows that overlap it (the pieces telescope), so energy is conserved
      !! per cell for ANY tiling of time by windows.  A bucket whose mean cosine
      !! is below `ATM_COSZ_MEAN_MIN` contributes 0.  The bound kept from the
      !! one-bucket kernel: each piece's share `I_p / I_B` is capped at 1 (the
      !! piece cannot carry more than its bucket; exact arithmetic never needs
      !! it, it guards round-off), i.e. `rsds <= (rec_s + rec_e) * bucket /
      !! wlen`.  `cosz` = the window mean of `max(0, cos Z)` (sum of both
      !! pieces over `wlen`), `cosz_mean` = the mean over bucket S.
      !! Explicit-shape, no cross-module calls; every scalar is host-computed.
      integer, intent(in), value :: nx_total, ny_total, nghost, nx_phys, ny_phys
      real(wp), intent(in), value :: hf_s, sdec_s, cdec_s, eq_s, hf_e, sdec_e, cdec_e, eq_e
      real(wp), intent(in), value :: t_off, wlen, bucket
      real(wp), intent(in) :: geolonT(nx_total, ny_total), geolatT(nx_total, ny_total)
      real(wp), intent(in) :: rec_s(nx_total, ny_total), rec_e(nx_total, ny_total)
      real(wp), intent(inout) :: cosz(nx_total, ny_total), cosz_mean(nx_total, ny_total)
      real(wp), intent(inout) :: rsds(nx_total, ny_total)

      integer :: i, j, ii, jj
      real(wp) :: lat_r, slat, clat, lon_deg, a_s, b_s, a_e, b_e, w0_s, w0_e, len2
      real(wp) :: int_s, int_e, tot_s, tot_e, mean_s, mean_e, f_s, f_e, acc

      len2 = t_off + wlen - bucket
      do concurrent(j=1:ny_phys, i=1:nx_phys)
         ii = nghost + i
         jj = nghost + j
         lat_r = geolatT(ii, jj)*DEG2RAD
         lon_deg = geolonT(ii, jj)
         slat = rdb_sm_sin(lat_r)
         clat = rdb_sm_cos(lat_r)
         a_s = slat*sdec_s
         b_s = clat*cdec_s
         a_e = slat*sdec_e
         b_e = clat*cdec_e
         w0_s = (15.0_wp*hf_s + 0.25_wp*eq_s + lon_deg - 180.0_wp)*DEG2RAD
         w0_e = (15.0_wp*hf_e + 0.25_wp*eq_e + lon_deg - 180.0_wp)*DEG2RAD

         tot_s = atm_cosz_cumulative(w0_s + ATM_EARTH_OMEGA*bucket, a_s, b_s) &
                 - atm_cosz_cumulative(w0_s, a_s, b_s)
         int_s = atm_cosz_cumulative(w0_s + ATM_EARTH_OMEGA*bucket, a_s, b_s) &
                 - atm_cosz_cumulative(w0_s + ATM_EARTH_OMEGA*t_off, a_s, b_s)
         tot_e = atm_cosz_cumulative(w0_e + ATM_EARTH_OMEGA*bucket, a_e, b_e) &
                 - atm_cosz_cumulative(w0_e, a_e, b_e)
         int_e = atm_cosz_cumulative(w0_e + ATM_EARTH_OMEGA*len2, a_e, b_e) &
                 - atm_cosz_cumulative(w0_e, a_e, b_e)
         tot_s = max(0.0_wp, tot_s)
         tot_e = max(0.0_wp, tot_e)
         int_s = max(0.0_wp, int_s)
         int_e = max(0.0_wp, int_e)
         mean_s = tot_s/(ATM_EARTH_OMEGA*bucket)
         mean_e = tot_e/(ATM_EARTH_OMEGA*bucket)
         cosz(ii, jj) = (int_s + int_e)/(ATM_EARTH_OMEGA*wlen)
         cosz_mean(ii, jj) = mean_s

         f_s = 0.0_wp
         f_e = 0.0_wp
         if (mean_s >= ATM_COSZ_MEAN_MIN) f_s = min(int_s/tot_s, 1.0_wp)
         if (mean_e >= ATM_COSZ_MEAN_MIN) f_e = min(int_e/tot_e, 1.0_wp)
         acc = rec_s(ii, jj)*f_s + rec_e(ii, jj)*f_e
         rsds(ii, jj) = acc*(bucket/wlen)
      end do
   end subroutine atm_coszen_straddle_impl

   pure function atm_cosz_cumulative(x, a, b) result(c)
      !! Antiderivative `C(x) = Integral_0^x max(0, a + b cos w) dw`
      !! of the clamped cosine of the zenith angle, `b >= 0`, up to an
      !! additive constant that cancels in every difference.  Three cases:
      !! `a >= b` (polar day: never clamped, `C = a x + b sin x`);
      !! `a + b <= 0` (polar night: `C = 0`); else the Sun is up for
      !! `|w mod 2 pi| < w0`, `cos w0 = -a/b`, and `C` adds the per-day
      !! integral `P = 2 (a w0 + b sin w0)` once per elapsed day.
      !! Device-callable (`routine seq`, same-module).
      !$acc routine seq
      real(wp), intent(in), value :: x, a, b
      real(wp) :: c

      real(wp) :: s0, w0, per, y, yc, n

      if (a >= b) then
         c = a*x + b*rdb_sm_sin(x)
      else if (a + b <= 0.0_wp) then
         c = 0.0_wp
      else
         w0 = rdb_sm_atan2(sqrt((b - a)*(b + a))/b, -a/b)
         s0 = rdb_sm_sin(w0)
         per = 2.0_wp*(a*w0 + b*s0)
         n = floor((x + 0.5_wp*TWO_PI)/TWO_PI)
         y = x - n*TWO_PI
         yc = max(-w0, min(w0, y))
         c = n*per + a*(yc + w0) + b*(rdb_sm_sin(yc) + s0)
      end if
   end function atm_cosz_cumulative

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

#include "rdb_safe_math_core.inc"

#include "rdb_safe_math_trig.inc"

#include "rdb_safe_math_atan.inc"

end module rdb_atm_state
