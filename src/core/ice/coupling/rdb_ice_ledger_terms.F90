!! The named-term conservation ledger of the CICE-port sea ice (roundabout-original).
module rdb_ice_ledger_terms
   !! `ice_term_ledger_t` books, for ONE ice step, every change of the three conserved
   !! quantities of the ice (ice + snow mass, salt, heat) against the sources and sinks the
   !! step booked, BY STAGE.  The sequencer (`ice_step_cice`) probes the state and the flux
   !! accumulators between the stages; a stage is a named TERM of the ledger.  Everything is
   !! summed with the order-invariant EFP reproducing sums (Hallberg and Adcroft, 2014): a sum does
   !! not depend on the summation order, the thread count or the device reduction tree, which is
   !! what lets the accumulators be restart state (a restart reproduces them bitwise).  The
   !! fluxes and the state being summed are the model's own, so two backends still differ by
   !! what the model differs by.
   !!
   !! **Conserved totals** (area-weighted over the owned cells, cell areas `A`): the mass
   !! `M = sum A (rhoi vice + rhos vsno)` [kg], the salt `S = sum A rhoi vice Sref 1e-3` [kg]
   !! (the flux convention: the ice is charged the reference salinity, not its stored
   !! profile) and the enthalpy `E = sum A sum_n (sum_k qice vicen/nilyr + sum_k qsno
   !! vsnon/nslyr)` [J].  Icepack's reference is FRESH LIQUID WATER AT 0 degC (enthalpy zero), so
   !! ice and snow enthalpy are negative.  Water that leaves the ocean at its own temperature
   !! carries `cp T dM` relative to that reference (`H_mass_conv`, about 2 % of the latent heat
   !! at the freezing point): it is booked NOWHERE, and is invisible only because the ocean
   !! volume never changes (the virtual-salt hand-off gives the ocean no water): zero by
   !! configuration, to reappear with the real-mass coupling.  The closure test computes its
   !! magnitude independently and reports it.
   !!
   !! **Stages = terms** (`ICE_LT_*`, in step order):
   !! | stage | what moves state | what it books |
   !! |---|---|---|
   !! | `therm1` | surface/basal growth and melt, snowfall, sublimation | atmosphere exchange (`fswabs`, `flwout`, `fsens`, vapour), ocean exchange, snow and longwave sources |
   !! | `remap` | rain pass-through, ITD remap (state-neutral) | rain: external source and `fresh` |
   !! | `newice` | frazil / new ice in open water | the latent heat (a SINK of stored energy), `fresh`/`fsalt` when `update_ocn_f` |
   !! | `lateral` | lateral melt | `fresh`, `fsalt`, `fhocn` |
   !! | `cleanup` | rebin, zaps of vanishing categories, aggregates | the zap residual in `fresh`, `fsalt`, `fhocn` |
   !! | `dynamics` | nothing (EVP changes velocities, not state) | NOT probed: its (null) effect is inside the `transport` probe, which spans both stages |
   !! | `transport` | advection (zero net on a closed domain) | the zaps of the post-transport cleanup, if any |
   !! | `ridging` | redistribution between categories | the ridged-snow fraction lost to the ocean |
   !!
   !! **The lines** (Icepack signs; `fresh`, `fsalt`, `fhocn` positive INTO the ocean,
   !! `fswthru` positive down; `evap > 0` condensation):
   !! - water:  `dM = [snow + rain + evap]_ext - dt sum A dfresh`;  `snow = dt sum A fsnow aice_init`
   !!   (falls on the ice-covered fraction of the START of the step), `rain = dt sum A frain aice_init`
   !!   (passes through: it enters `fresh` in the same stage, so its net effect on the ice is zero);
   !! - salt:   `dS = -dt sum A dfsalt`;
   !! - heat:   `dE = [lw_in - snow_enth + atm_rad + atm_lat + frz_lat]_ext - dt sum A (dfhocn + dfswthru)`
   !!   with `lw_in = dt sum A aice_init flw` (the merge stores the longwave net of the reflected
   !!   share), `snow_enth = Lfresh dt sum A aice_init fsnow` (falling snow carries `-rhos Lfresh`),
   !!   `atm_rad = dt sum A (fswabs + flwout + fsens)`, `atm_lat = Lvap dt sum A evap` (the vapour
   !!   side of the energy check; `flat` is NOT a term: it cancels inside `fsurf - flat`).
   !! **The sign of the new-ice energy (settled, D11.5).**  New ice is made of the pure latent
   !! enthalpy `qi0new = -rhoi Lfresh` (Bitz and Lipscomb, 1999), so freezing a bank of `E_frz =
   !! frzmlt dt` of heat stores `-E_frz`: the stored energy FALLS by the latent heat the ocean
   !! kept (the frazil clamp warmed the water back to freezing and banked that heat; the new ice
   !! carries the latent enthalpy against the 0 degC reference above).  The new-ice term is therefore
   !! `frz_lat = -rhoi Lfresh dt sum A frazil_vol`, a NEGATIVE source of the ice; the earlier
   !! docstring of `rdb_ice_ledger` that said `+ E_frz` had the sign wrong (the D11.3 smoke test
   !! was right).  Ocean plus ice conserve because the ocean's bank is `+E_frz` of ocean heat.
   !!
   !! **Residuals.**  Per stage and line `res = dstate - (sources - ocean-ward fluxes)`, its
   !! scale is the sum over the stages of every absolute term of the line (`scale(line)`), and the
   !! relative residual of a stage is `|res|/scale`.  There are no dropped-flux terms: this port
   !! conserves, so a non-zero residual is a bug (heat carries the vertical solver's tolerance,
   !! `ICE_FERRMAX`, in `therm1` only).  Mass and salt residuals are round-off.  ONE known extra in the `therm1`
   !! heat residual: the vertical-profile initialisation resets snow thinner than `hs_min` to the melting enthalpy
   !! before the thermodynamics' own energy check forms its initial energy (upstream behaviour, ported verbatim), an
   !! energy creation that check cannot see; `test_cice_closure` computes it independently (96 % of the `therm1` heat
   !! residual in a supercooled-ocean run, 1.8e12 J in 96 steps).
   !!
   !! **Cost and cadence.**  Opt-in (`&ocean_ice_nml ledger`).  One probe is one state-fill launch
   !! plus 6 host-synchronised EFP reductions (three totals, `fresh`, `fsalt`, `fhocn`); the
   !! `therm1` probe adds a fill and four reductions (`fswthru`, `evap`, atmosphere), the `newice`
   !! probe one reduction (`frazil_vol`); `begin` is 3 + 10 reductions.  Seven probes per step
   !! (the dynamics stage is not probed: EVP changes no total) and `therm2` runs as its four
   !! sweeps instead of three.  Off: the sequencer takes none of this and the run is bit-identical.
   !!
   !! **Restart.**  The cumulative residuals, the cumulative scales, the worst relative residual
   !! per line and the step counter are host scalars registered `optional` (`register_restart`):
   !! a checkpoint resumes them bitwise (they are sums of EFP-exact numbers in a fixed order).
   !! A file without them (ledger off when written) restarts the accumulators at zero.
   !!
   !! **Single rank.**  The sums are local; `begin` refuses a ledger on more than one rank
   !! (`valid = .false.`) rather than under-count.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_mem_report, only: arr_bytes
   use rdb_efp, only: efp_t, efp_carry, efp_to_real
   use rdb_ocean_console_stats, only: efp_decompose_impl, efp_summands_guard
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t
   implicit none
   private

   public :: ice_term_ledger_t, ice_ledger_tab_t
   public :: ice_lt_stage_name, ice_lt_line_name
   public :: ICE_LT_THERM1, ICE_LT_REMAP, ICE_LT_NEWICE, ICE_LT_LATERAL, ICE_LT_CLEANUP
   public :: ICE_LT_DYNAMICS, ICE_LT_TRANSPORT, ICE_LT_RIDGING, ICE_LT_NSTAGE
   public :: ICE_LT_WATER, ICE_LT_SALT, ICE_LT_HEAT

   integer, parameter :: ICE_LT_THERM1 = 1
   integer, parameter :: ICE_LT_REMAP = 2
   integer, parameter :: ICE_LT_NEWICE = 3
   integer, parameter :: ICE_LT_LATERAL = 4
   integer, parameter :: ICE_LT_CLEANUP = 5
   integer, parameter :: ICE_LT_DYNAMICS = 6
   integer, parameter :: ICE_LT_TRANSPORT = 7
   integer, parameter :: ICE_LT_RIDGING = 8
   integer, parameter :: ICE_LT_NSTAGE = 8
   integer, parameter :: ICE_LT_WATER = 1
   integer, parameter :: ICE_LT_SALT = 2
   integer, parameter :: ICE_LT_HEAT = 3

   type :: ice_ledger_tab_t
      !! One step's table (see the module docstring for the signs).
      real(wp) :: dstate(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Change of `M` [kg], `S` [kg], `E` [J] across each stage.
      real(wp) :: fout(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Ocean-ward `dt sum A d(fresh, fsalt, fhocn)` booked in each stage [kg, kg, J].
      real(wp) :: fsw(ICE_LT_NSTAGE) = 0.0_wp
         !! `dt sum A d(fswthru)` [J] booked in each stage.
      real(wp) :: snow = 0.0_wp, rain = 0.0_wp, evap = 0.0_wp
         !! Water sources [kg]: snowfall on the ice, rain through the ice, condensation.
      real(wp) :: lw_in = 0.0_wp, snow_enth = 0.0_wp, atm_rad = 0.0_wp, atm_lat = 0.0_wp, frz_lat = 0.0_wp
         !! Heat sources [J]: longwave in, enthalpy of falling snow (subtracted), shortwave
         !! absorbed + longwave out + sensible, vapour latent, new-ice latent (negative).
      real(wp) :: res(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Residual of each line in each stage.
      real(wp) :: scale(3) = 0.0_wp
         !! Sum of the absolute values of every term of each line.
   end type ice_ledger_tab_t

   type :: ice_term_ledger_t
      !! The ledger slot (module docstring).  Initialised and mapped ONLY under the new
      !! scheme with `&ocean_ice_nml ledger`.
      logical :: on = .false.
         !! The sequencer probes this ledger.
      logical :: is_init = .false.
      logical :: mapped = .false.
      logical :: is_open = .false.
      logical :: valid = .false.
         !! The last `end_step` produced residuals.
      logical :: finite = .true.
      integer :: nx = 0, ny = 0, ncat = 0, nilyr = 0, nslyr = 0
      integer :: is = 0, ie = 0, js = 0, je = 0
      integer :: nranks = 1
      real(wp) :: dt = 0.0_wp
      real(wp), allocatable :: area(:, :)
         !! Cell areas [m2] `(nx, ny)` (device copy owned here).
      real(wp), allocatable :: scr(:, :), scrs(:, :), scre(:, :)
         !! Per-cell scratch of the quantities being summed (mass, salt, energy).
      real(wp) :: m_prev = 0.0_wp, s_prev = 0.0_wp, e_prev = 0.0_wp
      real(wp) :: f_prev(4) = 0.0_wp
      type(ice_ledger_tab_t) :: tab
         !! The last closed step.
      real(wp) :: cum_res(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Cumulative residual per line and stage (restart state).
      real(wp) :: cum_scale(3) = 0.0_wp
         !! Cumulative scale per line (restart state).
      real(wp) :: worst(3) = 0.0_wp
         !! Worst per-step `max_stage |res|/scale` per line (restart state).
      real(wp) :: nsteps = 0.0_wp
         !! Closed steps (real so it registers as a restart scalar).
   contains
      procedure, non_overridable :: init => tl_init
      procedure, non_overridable :: set_area => tl_set_area
      procedure, non_overridable :: enter_data => tl_enter_data
      procedure, non_overridable :: exit_data => tl_exit_data
      procedure, non_overridable :: destroy => tl_destroy
      procedure, non_overridable :: bytes => tl_bytes
      procedure, non_overridable :: begin => tl_begin
      procedure, non_overridable :: probe => tl_probe
      procedure, non_overridable :: end_step => tl_end
      procedure, non_overridable :: totals => tl_totals
      procedure, non_overridable :: reset => tl_reset
      procedure, non_overridable :: register_restart => tl_register_restart
      procedure, non_overridable :: cum_rel => tl_cum_rel
      procedure, non_overridable :: step_rel => tl_step_rel
      procedure, non_overridable :: console_line => tl_console_line
   end type ice_term_ledger_t

contains

   pure function ice_lt_stage_name(s) result(name)
      integer, intent(in) :: s
      character(len=9) :: name

      select case (s)
      case (ICE_LT_THERM1)
         name = "therm1"
      case (ICE_LT_REMAP)
         name = "remap"
      case (ICE_LT_NEWICE)
         name = "newice"
      case (ICE_LT_LATERAL)
         name = "lateral"
      case (ICE_LT_CLEANUP)
         name = "cleanup"
      case (ICE_LT_DYNAMICS)
         name = "dynamics"
      case (ICE_LT_TRANSPORT)
         name = "transport"
      case (ICE_LT_RIDGING)
         name = "ridging"
      case default
         name = "?"
      end select
   end function ice_lt_stage_name

   pure function ice_lt_line_name(q) result(name)
      integer, intent(in) :: q
      character(len=5) :: name

      select case (q)
      case (ICE_LT_WATER)
         name = "water"
      case (ICE_LT_SALT)
         name = "salt"
      case (ICE_LT_HEAT)
         name = "heat"
      case default
         name = "?"
      end select
   end function ice_lt_line_name

   ! ------------------------------------------------------------------
   ! lifecycle
   ! ------------------------------------------------------------------

   subroutine tl_init(this, cice, on)
      !! Allocate for the structure of `cice`; host, once.  `on = .false.` leaves the slot
      !! inert (nothing allocated, nothing mapped).
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      logical, intent(in) :: on

      call this%destroy()
      this%on = on
      if (.not. on) return
      this%nx = cice%nx
      this%ny = cice%ny
      this%ncat = cice%ncat
      this%nilyr = cice%nilyr
      this%nslyr = cice%nslyr
      this%is = cice%nghost + 1
      this%ie = cice%nghost + cice%nx_phys
      this%js = cice%nghost + 1
      this%je = cice%nghost + cice%ny_phys
      allocate (this%area(this%nx, this%ny), source=0.0_wp)
      allocate (this%scr(this%nx, this%ny), source=0.0_wp)
      allocate (this%scrs(this%nx, this%ny), source=0.0_wp)
      allocate (this%scre(this%nx, this%ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine tl_init

   pure function tl_bytes(this) result(nbytes)
      !! Counted allocatable footprint (`area` and the three scratch fields).
      class(ice_term_ledger_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = arr_bytes(this%area) + arr_bytes(this%scr) + arr_bytes(this%scrs) + arr_bytes(this%scre)
   end function tl_bytes

   subroutine tl_set_area(this, tarea)
      !! Host, BEFORE `enter_data`: the T-cell areas.
      class(ice_term_ledger_t), intent(inout) :: this
      real(wp), intent(in) :: tarea(:, :)

      if (.not. this%is_init) return
      this%area = tarea(1:this%nx, 1:this%ny)
   end subroutine tl_set_area

   subroutine tl_enter_data(this)
      class(ice_term_ledger_t), intent(inout) :: this

      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%area, this%scr, this%scrs, this%scre)
      this%mapped = .true.
   end subroutine tl_enter_data

   subroutine tl_exit_data(this)
      class(ice_term_ledger_t), intent(inout) :: this

      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%area, this%scr, this%scrs, this%scre)
      this%mapped = .false.
   end subroutine tl_exit_data

   subroutine tl_destroy(this)
      class(ice_term_ledger_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      if (allocated(this%area)) deallocate (this%area)
      if (allocated(this%scr)) deallocate (this%scr)
      if (allocated(this%scrs)) deallocate (this%scrs)
      if (allocated(this%scre)) deallocate (this%scre)
      this%is_init = .false.
      this%on = .false.
      this%is_open = .false.
      this%valid = .false.
      this%nx = 0
      this%ny = 0
   end subroutine tl_destroy

   subroutine tl_reset(this)
      !! Clear the accumulators (a fresh run).
      class(ice_term_ledger_t), intent(inout) :: this

      this%cum_res = 0.0_wp
      this%cum_scale = 0.0_wp
      this%worst = 0.0_wp
      this%nsteps = 0.0_wp
      this%valid = .false.
      this%tab = ice_ledger_tab_t()
   end subroutine tl_reset

   subroutine tl_register_restart(this, reg)
      !! Register the accumulators (all `optional`).  No-op when the ledger is off.
      class(ice_term_ledger_t), intent(in), target :: this
      type(restart_registry_t), intent(inout) :: reg
      integer :: q, s
      character(len=64) :: tag

      if (.not. this%is_init) return
      do s = 1, ICE_LT_NSTAGE
         do q = 1, 3
            write (tag, "(a,a,a,a)") "ice_led_res_", trim(ice_lt_line_name(q)), "_", trim(ice_lt_stage_name(s))
            call reg%register_scalar(trim(tag), this%cum_res(q, s), optional=.true.)
         end do
      end do
      do q = 1, 3
         write (tag, "(a,a)") "ice_led_scale_", trim(ice_lt_line_name(q))
         call reg%register_scalar(trim(tag), this%cum_scale(q), optional=.true.)
         write (tag, "(a,a)") "ice_led_worst_", trim(ice_lt_line_name(q))
         call reg%register_scalar(trim(tag), this%worst(q), optional=.true.)
      end do
      call reg%register_scalar("ice_led_nsteps", this%nsteps, optional=.true.)
   end subroutine tl_register_restart

   ! ------------------------------------------------------------------
   ! EFP sums (device)
   ! ------------------------------------------------------------------

   subroutine efp_sum_prod(nx, ny, is, ie, js, je, a, x, total)
      !! `sum a(i,j) x(i,j)` over the owned cells with the order-invariant EFP bins: bitwise the
      !! same on any backend and thread count.  One reduction of seven integer accumulators (the
      !! proven footprint).  NaN or Inf in a summand returns NaN.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: a(nx, ny), x(nx, ny)
      real(wp), intent(out) :: total
      type(efp_t) :: acc
      integer :: i, j
      integer(int64) :: e1, e2, e3, e4, e5, e6, ep, d1, d2, d3, d4, d5, d6, dp
      real(real64) :: v

      call efp_summands_guard(ie - is + 1, je - js + 1, "ice ledger EFP sum")
      e1 = 0_int64
      e2 = 0_int64
      e3 = 0_int64
      e4 = 0_int64
      e5 = 0_int64
      e6 = 0_int64
      ep = 0_int64
      ! acc-kept: EFP reproducing-sum path (maintainer rule: sums that must be bitwise reproducible keep the directive form)
      !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6,ep) &
      !$acc&    private(v, d1, d2, d3, d4, d5, d6, dp) present(a, x)
      do j = js, je
         do i = is, ie
            v = real(a(i, j)*x(i, j), real64)
            call efp_decompose_impl(v, d1, d2, d3, d4, d5, d6, dp)
            e1 = e1 + d1
            e2 = e2 + d2
            e3 = e3 + d3
            e4 = e4 + d4
            e5 = e5 + d5
            e6 = e6 + d6
            ep = ep + dp
         end do
      end do
      acc%v = [e1, e2, e3, e4, e5, e6]
      acc%poison = ep
      call efp_carry(acc%v)
      total = real(efp_to_real(acc), wp)
   end subroutine efp_sum_prod

   subroutine fill_state(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, rhoi, rhos, sref, &
                         vicen, vsnon, qicen, qsnon, scr_m, scr_s, scr_e)
      !! Per-cell `M`, `S`, `E` from the category arrays, ONE launch.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, is, ie, js, je
      real(wp), intent(in) :: rhoi, rhos, sref
      real(wp), intent(in) :: vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: scr_m(nx, ny), scr_s(nx, ny), scr_e(nx, ny)
      integer :: i, j, n, k
      real(wp) :: vi, vs, e

      !$acc data present(vicen, vsnon, qicen, qsnon, scr_m, scr_s, scr_e)
      do concurrent(j=js:je, i=is:ie) local(n, k, vi, vs, e)
         vi = 0.0_wp
         vs = 0.0_wp
         e = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            vi = vi + vicen(i, j, n)
            vs = vs + vsnon(i, j, n)
            do k = 1, nilyr
               e = e + qicen(i, j, n, k)*vicen(i, j, n)/real(nilyr, wp)
            end do
            do k = 1, nslyr
               e = e + qsnon(i, j, n, k)*vsnon(i, j, n)/real(nslyr, wp)
            end do
         end do
         scr_m(i, j) = rhoi*vi + rhos*vs
         scr_s(i, j) = rhoi*vi*sref*1.0e-3_wp
         scr_e(i, j) = e
      end do
      !$acc end data
   end subroutine fill_state

   subroutine fill_atm(nx, ny, is, ie, js, je, Lvap, fswabs, flwout, fsens, evap, scr_rad, scr_lat)
      !! Per-cell `fswabs + flwout + fsens` and `Lvap evap`, one launch.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: Lvap
      real(wp), intent(in) :: fswabs(nx, ny), flwout(nx, ny), fsens(nx, ny), evap(nx, ny)
      real(wp), intent(inout) :: scr_rad(nx, ny), scr_lat(nx, ny)
      integer :: i, j

      !$acc data present(fswabs, flwout, fsens, evap, scr_rad, scr_lat)
      do concurrent(j=js:je, i=is:ie)
         scr_rad(i, j) = (fswabs(i, j) + flwout(i, j)) + fsens(i, j)
         scr_lat(i, j) = Lvap*evap(i, j)
      end do
      !$acc end data
   end subroutine fill_atm

   subroutine fill_prod(nx, ny, is, ie, js, je, x, y, scr)
      !! `scr = x*y` over the owned cells.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: x(nx, ny), y(nx, ny)
      real(wp), intent(inout) :: scr(nx, ny)
      integer :: i, j

      !$acc data present(x, y, scr)
      do concurrent(j=js:je, i=is:ie)
         scr(i, j) = x(i, j)*y(i, j)
      end do
      !$acc end data
   end subroutine fill_prod

   subroutine state_totals(this, cice, m, s, e)
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(out) :: m, s, e

      call fill_state(this%nx, this%ny, this%ncat, this%nilyr, this%nslyr, this%is, this%ie, this%js, this%je, &
                      cice%p%c%rhoi, cice%p%c%rhos, cice%p%thermo%ice_ref_salinity, &
                      cice%vicen, cice%vsnon, cice%qicen, cice%qsnon, this%scr, this%scrs, this%scre)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, m)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scrs, s)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scre, e)
   end subroutine state_totals

   subroutine flux_totals(this, flux, f, with_sw)
      !! `sum A (fresh, fsalt, fhocn)` [flux units x m2] and, when `with_sw`, `sum A fswthru` (only the
      !! `therm1` stage writes it: the other probes carry the previous value).
      class(ice_term_ledger_t), intent(in) :: this
      type(ice_flux_t), intent(in) :: flux
      real(wp), intent(inout) :: f(4)
      logical, intent(in) :: with_sw

      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fresh, f(1))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fsalt, f(2))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fhocn, f(3))
      if (with_sw) then
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fswthru, f(4))
      end if
   end subroutine flux_totals

   ! ------------------------------------------------------------------
   ! the step
   ! ------------------------------------------------------------------

   subroutine tl_totals(this, cice, m, s, e)
      !! The conserved totals `M` [kg], `S` [kg], `E` [J] of the ice state now (EFP sums).  Host call;
      !! the state and the ledger must be device-mapped.  Used by the tests and the console.
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(out) :: m, s, e

      m = 0.0_wp
      s = 0.0_wp
      e = 0.0_wp
      if (.not. this%is_init .or. .not. this%mapped) return
      call state_totals(this, cice, m, s, e)
   end subroutine tl_totals

   subroutine tl_begin(this, cice, flux, dt, fsnow, frain, flw, nranks, with_sources)
      !! Open a step of length `dt` [s]: snapshot the totals (the fluxes are zero) and form the
      !! external sources of the step from this step's atmosphere inputs and the start-of-step
      !! concentration `aice_init` (call it after `ice_save_init`).  `fsnow`, `frain` [kg m-2 s-1]
      !! and `flw` [W m-2] are per cell area.
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_flux_t), intent(in) :: flux
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: fsnow(this%nx, this%ny), frain(this%nx, this%ny), flw(this%nx, this%ny)
      integer, intent(in) :: nranks
      logical, intent(in) :: with_sources
         !! `.false.` when the step skips the thermodynamics (no atmosphere exchange: the sources stay zero).
      real(wp) :: v

      this%valid = .false.
      this%is_open = .false.
      this%nranks = nranks
      if (.not. this%is_init .or. .not. this%mapped .or. nranks > 1) return
      this%dt = dt
      this%tab = ice_ledger_tab_t()
      call state_totals(this, cice, this%m_prev, this%s_prev, this%e_prev)
      this%f_prev = 0.0_wp
      call flux_totals(this, flux, this%f_prev, .true.)
      if (.not. with_sources) then
         this%is_open = .true.
         return
      end if
      call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, fsnow, cice%aice_init, this%scr)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, v)
      this%tab%snow = dt*v
      this%tab%snow_enth = cice%p%c%Lfresh*dt*v
      call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, frain, cice%aice_init, this%scr)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, v)
      this%tab%rain = dt*v
      call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, flw, cice%aice_init, this%scr)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, v)
      this%tab%lw_in = dt*v
      this%is_open = .true.
   end subroutine tl_begin

   subroutine tl_probe(this, stage, cice, flux)
      !! Close stage `stage` (`ICE_LT_*`): the change of the totals and of the flux accumulators
      !! since the previous probe are booked to it.  After `therm1` the atmosphere terms are
      !! read; after `newice` the new-ice volume.
      class(ice_term_ledger_t), intent(inout) :: this
      integer, intent(in) :: stage
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_flux_t), intent(in) :: flux
      real(wp) :: m, s, e, v
      real(wp) :: f(4)

      if (.not. this%is_open) return
      call state_totals(this, cice, m, s, e)
      f = this%f_prev
      call flux_totals(this, flux, f, stage == ICE_LT_THERM1)
      this%tab%dstate(1, stage) = m - this%m_prev
      this%tab%dstate(2, stage) = s - this%s_prev
      this%tab%dstate(3, stage) = e - this%e_prev
      this%tab%fout(1, stage) = this%dt*(f(1) - this%f_prev(1))
      this%tab%fout(2, stage) = this%dt*(f(2) - this%f_prev(2))
      this%tab%fout(3, stage) = this%dt*(f(3) - this%f_prev(3))
      this%tab%fsw(stage) = this%dt*(f(4) - this%f_prev(4))
      this%m_prev = m
      this%s_prev = s
      this%e_prev = e
      this%f_prev = f
      if (stage == ICE_LT_THERM1) then
         call fill_atm(this%nx, this%ny, this%is, this%ie, this%js, this%je, cice%p%c%Lvap, flux%fswabs, &
                       flux%flwout, flux%fsens, flux%evap, this%scr, this%scrs)
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, v)
         this%tab%atm_rad = this%dt*v
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scrs, v)
         this%tab%atm_lat = this%dt*v
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%evap, v)
         this%tab%evap = this%dt*v
      else if (stage == ICE_LT_NEWICE) then
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%frazil_vol, v)
         this%tab%frz_lat = -cice%p%c%rhoi*cice%p%c%Lfresh*v
      end if
   end subroutine tl_probe

   subroutine tl_end(this)
      !! Close the step: evaluate the residuals of every stage and line, add them to the
      !! cumulative accumulators.
      class(ice_term_ledger_t), intent(inout) :: this
      integer :: s, q
      real(wp) :: r
      real(wp) :: srcw(ICE_LT_NSTAGE), srch(ICE_LT_NSTAGE), worst_step(3)
      logical :: fin

      this%valid = .false.
      if (.not. this%is_open) return
      this%is_open = .false.
      associate (t => this%tab)
         srcw = 0.0_wp
         srch = 0.0_wp
         srcw(ICE_LT_THERM1) = t%snow + t%evap
         srcw(ICE_LT_REMAP) = t%rain
         srch(ICE_LT_THERM1) = t%lw_in - t%snow_enth + t%atm_rad + t%atm_lat
         srch(ICE_LT_NEWICE) = t%frz_lat
         t%scale = 0.0_wp
         do s = 1, ICE_LT_NSTAGE
            t%res(ICE_LT_WATER, s) = t%dstate(1, s) - (srcw(s) - t%fout(1, s))
            t%res(ICE_LT_SALT, s) = t%dstate(2, s) + t%fout(2, s)
            t%res(ICE_LT_HEAT, s) = t%dstate(3, s) - (srch(s) - t%fout(3, s) - t%fsw(s))
            t%scale(ICE_LT_WATER) = t%scale(ICE_LT_WATER) + abs(t%dstate(1, s)) + abs(srcw(s)) + abs(t%fout(1, s))
            t%scale(ICE_LT_SALT) = t%scale(ICE_LT_SALT) + abs(t%dstate(2, s)) + abs(t%fout(2, s))
            t%scale(ICE_LT_HEAT) = t%scale(ICE_LT_HEAT) + abs(t%dstate(3, s)) + abs(srch(s)) + abs(t%fout(3, s)) &
                                   + abs(t%fsw(s))
         end do
         fin = all(ieee_is_finite(t%res)) .and. all(ieee_is_finite(t%scale))
         this%finite = fin
         if (.not. fin) then
            ! STICKY: a non-finite step poisons the accumulators (and so the console) for good
            t%res = ieee_value(1.0_wp, ieee_quiet_nan)
            this%cum_res = ieee_value(1.0_wp, ieee_quiet_nan)
            this%cum_scale = ieee_value(1.0_wp, ieee_quiet_nan)
            this%worst = ieee_value(1.0_wp, ieee_quiet_nan)
            this%nsteps = this%nsteps + 1.0_wp
            this%valid = .true.
            return
         end if
         worst_step = 0.0_wp
         do q = 1, 3
            do s = 1, ICE_LT_NSTAGE
               r = abs(t%res(q, s))/max(t%scale(q), tiny(1.0_wp))
               worst_step(q) = max(worst_step(q), r)
               this%cum_res(q, s) = this%cum_res(q, s) + t%res(q, s)
            end do
            this%cum_scale(q) = this%cum_scale(q) + t%scale(q)
            this%worst(q) = max(this%worst(q), worst_step(q))
         end do
      end associate
      this%nsteps = this%nsteps + 1.0_wp
      this%valid = .true.
   end subroutine tl_end

   pure function tl_step_rel(this, q, s) result(r)
      !! `|res|/scale` of line `q` in stage `s` of the last closed step (stage 0: summed over the
      !! stages, the whole-step residual); NaN if not valid.
      class(ice_term_ledger_t), intent(in) :: this
      integer, intent(in) :: q, s
      real(wp) :: r

      r = ieee_value(1.0_wp, ieee_quiet_nan)
      if (.not. this%valid .or. .not. this%finite) return
      if (q < 1 .or. q > 3 .or. s < 0 .or. s > ICE_LT_NSTAGE) return
      if (s == 0) then
         r = abs(sum(this%tab%res(q, :)))/max(this%tab%scale(q), tiny(1.0_wp))
      else
         r = abs(this%tab%res(q, s))/max(this%tab%scale(q), tiny(1.0_wp))
      end if
   end function tl_step_rel

   pure function tl_cum_rel(this, q) result(r)
      !! Cumulative residual of line `q` relative to the cumulative scale: the largest over the
      !! stages of `|cum_res|/cum_scale`.  NaN as soon as any accumulator is non-finite.
      class(ice_term_ledger_t), intent(in) :: this
      integer, intent(in) :: q
      real(wp) :: r

      if (.not. all(ieee_is_finite(this%cum_res(q, :))) .or. .not. ieee_is_finite(this%cum_scale(q))) then
         r = ieee_value(1.0_wp, ieee_quiet_nan)
      else
         r = maxval(abs(this%cum_res(q, :)))/max(this%cum_scale(q), tiny(1.0_wp))
      end if
   end function tl_cum_rel

   function tl_console_line(this) result(line)
      !! The console line: cumulative relative residual per line, the stage that carries the
      !! largest one, and the worst single-step value.
      class(ice_term_ledger_t), intent(in) :: this
      character(len=200) :: line
      character(len=24) :: w(3), st(3)
      integer :: q, s, sw
      real(wp) :: c

      do q = 1, 3
         sw = 1
         c = -1.0_wp
         do s = 1, ICE_LT_NSTAGE
            if (abs(this%cum_res(q, s)) > c) then
               c = abs(this%cum_res(q, s))
               sw = s
            end if
         end do
         write (w(q), "(es9.2)") this%cum_rel(q)
         write (st(q), "(a)") trim(ice_lt_stage_name(sw))
      end do
      write (line, "(a,i0,13a)") "Ice ledger: ", nint(this%nsteps), " steps, |res|/scale water ", &
         trim(w(1)), " (", trim(st(1)), "), salt ", trim(w(2)), " (", trim(st(2)), "), heat ", trim(w(3)), " (", &
         trim(st(3)), ")"
   end function tl_console_line

end module rdb_ice_ledger_terms
