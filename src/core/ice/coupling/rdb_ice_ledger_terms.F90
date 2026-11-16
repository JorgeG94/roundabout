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
   !! | `transport` | advection (zero net on a closed domain) | the zaps of the post-transport cleanup, if any; the remap's CLIP RECORD (below) |
   !! | `ridging` | redistribution between categories | the ridged-snow fraction lost to the ocean |
   !!
   !! **The transport clip record (named term of `transport`).**  The remap update clips a new area in `[-puny, 0)` to 0
   !! and leaves the tracer products of a cell whose new area (or parent) is 0 behind: that content is DROPPED from the
   !! state with no flux to the ocean.  The remap writes it per cell (`remap_flux_workspace_t%clip`: clipped count, area
   !! created, ice and snow volume, ice and snow enthalpy; the snow one in the packed frame, shifted by `rhos Lfresh`),
   !! `book_clip` reduces the planes with the EFP sums, and `end_step` books it as a SOURCE of the `transport` stage
   !! (`clip_w = -(rhoi vice + rhos vsno)`, `clip_s = -rhoi vice Sref 1e-3`, `clip_h = -(eice + esno - rhos Lfresh vsno)`,
   !! all NEGATIVE: sinks of the state), so the stage closes to round-off whatever the guards removed.  It was zero in
   !! every run so far; the term is what makes that a measured statement.  (The clip area is reported, not a line of the
   !! ledger.)  With `advection = "upwind"` no record exists and none is booked: the first-order scheme's own
   !! sub-`puny` drop (a new area in `(-puny, 0]` unpacked as 0 with its volume, at most `puny` times the thickness per
   !! cell, round-off of the scheme) stays in the `transport` residual, where it has been below the closure tolerance in every run.
   !!
   !! **The zaps (named term of `cleanup`, and of the post-transport and post-ridging clean-ups).**  The ITD clean-up's
   !! zaps (a category with `|aicen| <= puny`, the round-off excess `1 < aice < 1 + puny`, cells with `aice` or the ice
   !! mass below the minimum, snow outside its temperature range) hand the removed ice, snow, salt and enthalpy to the
   !! ocean through `fresh`, `fsalt`, `fhocn` with the flux slot's signs (mass and salt positive, heat negative): in the
   !! table they are `fout(:, cleanup)` (and `fout(:, transport)`), and `cum_zap` accumulates them.
   !!
   !! **The thin-snow reset (REPORTED, not booked).**  The vertical-profile initialisation resets snow thinner than
   !! `hs_min` to the melting enthalpy `-rhos Lfresh` BEFORE the thermodynamics' own energy check forms its initial
   !! energy (upstream behaviour, awaiting a maintainer decision; the physics is unchanged).  The energy that adds,
   !! `sum A (-rhos Lfresh - q_k) vsno/nslyr` over the categories with `aicen > puny` and a snow layer thickness not above
   !! `hs_min/nslyr`, from the START-of-step state, is `tab%snow_reset` (cumulative `cum_snow_reset`, restart
   !! `ice_led_snow_reset`): a diagnostic that ATTRIBUTES part of the `therm1` heat residual; the closure bound is
   !! unchanged and the term is not subtracted.
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
   !! **Cost and cadence.**  Opt-in (`&ocean_ice_nml ledger`).  A probe is one state-fill launch plus 6 host-synchronised
   !! LOCAL EFP reductions (three totals, `fresh`, `fsalt`, `fhocn`); the `therm1` probe adds a fill and four more
   !! (`fswthru`, `evap`, atmosphere), the `newice` probe one (`frazil_vol`); `begin` is 3 + 5 (the state, the three
   !! sources, the thin-snow reset and the heat bound, one more fill each for the products); the transport clip record
   !! adds one launch inside the remap and six reductions.  Seven probes per step (the dynamics stage is not probed: EVP
   !! changes no total) and `therm2` runs as its four sweeps instead of three.  Off: the sequencer takes none of this
   !! and the run is bit-identical (the ledger is an OBSERVER: it reads the state and the flux accumulators and writes
   !! only its own arrays).
   !!
   !! **Several ranks.**  Every probe sums over this rank's owned cells only and keeps the LOCAL EFP value (the
   !! integer bins and the poison counter, `lv`); nothing is communicated until `end_step`, which combines the whole
   !! step's vector (`NV` values: every total, flux sum, source, the clip record, the thin-snow reset and the bound) in
   !! ONE collective (`seam_sum_efp`, bin-wise over the ranks, poison counters included) and only then forms the
   !! differences, residuals and scales on the host.  That is ONE collective per ice step, a no-op on one rank.  The
   !! combined bins are exact integers, so every number the ledger reports is the GLOBAL value and is the same bits
   !! on every rank and on every `px x py` split, and the same as one rank's (`test_ocean_decomp_bitid_mpi`); a
   !! non-finite summand on any rank makes every rank read NaN identically (the poison counter rides the same
   !! `MPI_SUM`).  The ledger must be called on every rank at the same places (the sequencer's decisions are
   !! rank-uniform, so a stopped step stops everywhere before `end_step`).
   !!
   !! **Restart.**  The cumulative residuals, the cumulative scales, the worst relative residual per line, the
   !! cumulative transport clip record, the cumulative thin-snow reset and zap terms, and the step counter are host scalars
   !! registered `optional` (`register_restart`).  They are built from the GLOBAL step values (after `end_step`'s
   !! reduction), so EVERY rank holds the same numbers: a decomposed run's per-rank checkpoint resumes bitwise (each
   !! rank restores the global value) and its cumulative report equals the one-rank report.  The alternative
   !! (per-rank partials reduced at report time) would need a collective at every report and would make the file
   !! depend on the decomposition.  A file without them (ledger off when written) restarts the accumulators at zero.
   !! The older `ice_remap_*` scalars of `rdb_ice_cice_run` stay rank-local plain sums (the transport result's own
   !! record, ledger on or off); the ledger's global clip record is `ice_led_clip_*`.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_mem_report, only: arr_bytes
   use rdb_efp, only: efp_t, efp_carry, efp_to_real
   use rdb_ocean_console_stats, only: efp_decompose_impl, efp_summands_guard
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t
   use rdb_ice_bl99, only: ICE_FERRMAX
   use rdb_ice_transport_seam, only: transport_seam_t, seam_sum_efp
   implicit none
   private

   public :: ice_term_ledger_t, ice_ledger_tab_t
   public :: ice_lt_stage_name, ice_lt_line_name
   public :: ICE_LT_THERM1, ICE_LT_REMAP, ICE_LT_NEWICE, ICE_LT_LATERAL, ICE_LT_CLEANUP
   public :: ICE_LT_DYNAMICS, ICE_LT_TRANSPORT, ICE_LT_RIDGING, ICE_LT_NSTAGE
   public :: ICE_LT_WATER, ICE_LT_SALT, ICE_LT_HEAT
   public :: ICE_LT_NV, ICE_LT_NCLIP, ice_lt_ix_state, ice_lt_ix_flux, ice_lt_ix_clip
   public :: ice_lt_test_local

   integer, save :: ice_lt_test_local = 0
      !! TEST HOOK, never set by the model: the index (`ice_lt_ix_*`) of ONE value of the step's vector that `end_step`
      !! leaves RANK-LOCAL (its cross-rank sum is discarded).  The decomposition test sets it to show that a ledger
      !! term left unreduced is detected.  0: every value is reduced.

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

   integer, parameter :: ICE_LT_NCLIP = 6
      !! Values of the transport clip record: clipped cells, area created, ice and snow volume dropped, ice and snow
      !! enthalpy dropped (`remap_ledger_t`, per cell, in the planes of `remap_flux_workspace_t%clip`).
   integer, parameter :: IX_FLUX0 = 3*(ICE_LT_NSTAGE + 1)
   integer, parameter :: IX_SW = IX_FLUX0 + 3*ICE_LT_NSTAGE + 1
   integer, parameter :: IX_ATM = IX_SW
   integer, parameter :: IX_FRZ = IX_ATM + 4
   integer, parameter :: IX_SRC = IX_FRZ
   integer, parameter :: IX_CLIP = IX_SRC + 4
   integer, parameter :: IX_SNR = IX_CLIP + ICE_LT_NCLIP
   integer, parameter :: IX_BND = IX_SNR + 1
   integer, parameter :: ICE_LT_NV = IX_BND
      !! Number of EFP values a step reduces (ONE collective): 3 state totals at each of the `begin` and the 8
      !! stages, 3 flux sums per stage, `fswthru`, 3 atmosphere sums, the new-ice volume, 3 sources, the 6 clip
      !! values, the thin-snow reset and the heat bound.

   type :: ice_ledger_tab_t
      !! One step's table (see the module docstring for the signs).
      real(wp) :: dstate(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Change of `M` [kg], `S` [kg], `E` [J] across each stage.
      real(wp) :: fout(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Ocean-ward `dt sum A d(fresh, fsalt, fhocn)` booked in each stage [kg, kg, J].  In the `cleanup` stage (and in
         !! the post-transport and post-ridging clean-ups) this IS the zap term: the ice, snow, salt and enthalpy the
         !! ITD clean-up removed (small areas, the round-off excess `aice > 1`, near-empty cells, snow outside its
         !! temperature range) and handed to the ocean, with the flux slot's signs (mass and salt positive, the heat
         !! NEGATIVE because the removed enthalpy is).
      real(wp) :: fsw(ICE_LT_NSTAGE) = 0.0_wp
         !! `dt sum A d(fswthru)` [J] booked in each stage.
      real(wp) :: snow = 0.0_wp, rain = 0.0_wp, evap = 0.0_wp
         !! Water sources [kg]: snowfall on the ice, rain through the ice, condensation.
      real(wp) :: lw_in = 0.0_wp, snow_enth = 0.0_wp, atm_rad = 0.0_wp, atm_lat = 0.0_wp, frz_lat = 0.0_wp
         !! Heat sources [J]: longwave in, enthalpy of falling snow (subtracted), shortwave
         !! absorbed + longwave out + sensible, vapour latent, new-ice latent (negative).
      real(wp) :: n_clip = 0.0_wp, clip_area = 0.0_wp, drop_vice = 0.0_wp, drop_vsno = 0.0_wp
      real(wp) :: drop_eice = 0.0_wp, drop_esno = 0.0_wp
         !! The remap transport's clip record of the step (global): cells whose new area was clipped to 0, the area
         !! that created [m2], and the ice and snow volume [m3] and enthalpy [J, the snow one in the packed frame]
         !! left behind where the new area (or its parent) is 0 and so dropped from the state.
      real(wp) :: clip_w = 0.0_wp, clip_s = 0.0_wp, clip_h = 0.0_wp
         !! The same clip record as a SOURCE of the transport stage in the three lines [kg, kg, J] (negative: a
         !! sink of the state): `-(rhoi vice + rhos vsno)`, `-rhoi vice Sref 1e-3`, `-(eice + esno - rhos Lfresh vsno)`.
      real(wp) :: snow_reset = 0.0_wp
         !! DIAGNOSTIC, not a booked term [J]: the energy the vertical-profile initialisation of the
         !! thermodynamics ADDS to the state by resetting the snow thinner than `hs_min` to the melting enthalpy
         !! `-rhos Lfresh` before its own energy check forms its initial energy (a known unbooked upstream behaviour:
         !! `sum A (-rhos Lfresh - q_k) vsno/nslyr` over the layers of the thinner-than-`hs_min` snow of categories with
         !! `aicen > puny`, from the START-of-step state).  It is part of the `therm1` heat residual and the closure
         !! bound does not subtract it.
      real(wp) :: heat_bound = 0.0_wp
         !! The solver's accepted energy error of the step, `0.9 ICE_FERRMAX dt sum A aice_init` [J].
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
      logical :: with_sources = .false.
         !! The open step books atmosphere sources (the thermodynamics ran).
      logical :: probed(0:ICE_LT_NSTAGE) = .false.
         !! Which probes of the open step ran (0: `begin`).
      logical :: have_clip = .false.
         !! `book_clip` ran in the open step.
      integer :: nx = 0, ny = 0, ncat = 0, nilyr = 0, nslyr = 0
      integer :: is = 0, ie = 0, js = 0, je = 0
      real(wp) :: dt = 0.0_wp
      real(wp) :: rhoi = 0.0_wp, rhos = 0.0_wp, lfresh = 0.0_wp, sref = 0.0_wp
         !! Constants of the open step (so `end_step` needs no state).
      real(wp), allocatable :: area(:, :)
         !! Cell areas [m2] `(nx, ny)` (device copy owned here).
      real(wp), allocatable :: scr(:, :), scrs(:, :), scre(:, :)
         !! Per-cell scratch of the quantities being summed (mass, salt, energy).
      type(efp_t) :: lv(ICE_LT_NV)
         !! THIS RANK's EFP values of the open step (`ice_lt_ix_*`); combined in ONE collective by `end_step`.
      type(ice_ledger_tab_t) :: tab
         !! The last closed step.
      real(wp) :: cum_res(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Cumulative residual per line and stage (restart state).
      real(wp) :: cum_scale(3) = 0.0_wp
         !! Cumulative scale per line (restart state).
      real(wp) :: worst(3) = 0.0_wp
         !! Worst per-step `max_stage |res|/scale` per line (restart state).
      real(wp) :: cum_clip(ICE_LT_NCLIP) = 0.0_wp
         !! Cumulative GLOBAL transport clip record (restart state): clipped cells, area, vice, vsno, eice, esno.
      real(wp) :: cum_snow_reset = 0.0_wp
         !! Cumulative thin-snow reset energy [J] (restart state).
      real(wp) :: cum_zap(3) = 0.0_wp
         !! Cumulative ocean-ward flux of the clean-up zaps (`cleanup` + `transport` stages) in the three lines (restart state).
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
      procedure, non_overridable :: book_clip => tl_book_clip
      procedure, non_overridable :: end_step => tl_end
      procedure, non_overridable :: totals => tl_totals
      procedure, non_overridable :: reset => tl_reset
      procedure, non_overridable :: register_restart => tl_register_restart
      procedure, non_overridable :: cum_rel => tl_cum_rel
      procedure, non_overridable :: step_rel => tl_step_rel
      procedure, non_overridable :: console_line => tl_console_line
   end type ice_term_ledger_t

contains

   pure function ice_lt_ix_state(q, s) result(ix)
      !! Index in `lv` of the total `q` (`ICE_LT_WATER..HEAT`) after stage `s` (0: at `begin`).
      integer, intent(in) :: q, s
      integer :: ix

      ix = q + 3*s
   end function ice_lt_ix_state

   pure function ice_lt_ix_flux(q, s) result(ix)
      !! Index in `lv` of `sum A (fresh, fsalt, fhocn)(q)` after stage `s` (1..`ICE_LT_NSTAGE`).
      integer, intent(in) :: q, s
      integer :: ix

      ix = IX_FLUX0 + q + 3*(s - 1)
   end function ice_lt_ix_flux

   pure function ice_lt_ix_clip(k) result(ix)
      !! Index in `lv` of clip value `k` (1..`ICE_LT_NCLIP`).
      integer, intent(in) :: k
      integer :: ix

      ix = IX_CLIP + k - 1
   end function ice_lt_ix_clip

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
      this%cum_clip = 0.0_wp
      this%cum_snow_reset = 0.0_wp
      this%cum_zap = 0.0_wp
      this%nsteps = 0.0_wp
      this%valid = .false.
      this%tab = ice_ledger_tab_t()
   end subroutine tl_reset

   subroutine tl_register_restart(this, reg)
      !! Register the accumulators (all `optional`).  No-op when the ledger is off.  They are GLOBAL values (module
      !! docstring), the same on every rank.
      class(ice_term_ledger_t), intent(in), target :: this
      type(restart_registry_t), intent(inout) :: reg
      integer :: q, s
      character(len=64) :: tag
      character(len=7), parameter :: CLIPNAME(ICE_LT_NCLIP) = &
                                     [character(len=7) :: "n", "area", "vice", "vsno", "eice", "esno"]

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
         write (tag, "(a,a)") "ice_led_zap_", trim(ice_lt_line_name(q))
         call reg%register_scalar(trim(tag), this%cum_zap(q), optional=.true.)
      end do
      do q = 1, ICE_LT_NCLIP
         write (tag, "(a,a)") "ice_led_clip_", trim(CLIPNAME(q))
         call reg%register_scalar(trim(tag), this%cum_clip(q), optional=.true.)
      end do
      call reg%register_scalar("ice_led_snow_reset", this%cum_snow_reset, optional=.true.)
      call reg%register_scalar("ice_led_nsteps", this%nsteps, optional=.true.)
   end subroutine tl_register_restart

   ! ------------------------------------------------------------------
   ! EFP sums (device)
   ! ------------------------------------------------------------------

   subroutine efp_sum_prod(nx, ny, is, ie, js, je, a, x, acc)
      !! `sum a(i,j) x(i,j)` over the owned cells of THIS rank with the order-invariant EFP bins: bitwise the
      !! same on any backend and thread count.  One reduction of seven integer accumulators (the proven footprint).
      !! The result is the LOCAL value (carried bins and the poison counter; NaN or Inf in a summand poisons it),
      !! to be combined over the ranks by `end_step`.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: a(nx, ny), x(nx, ny)
      type(efp_t), intent(out) :: acc
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
   end subroutine efp_sum_prod

   subroutine efp_sum_plane(nx, ny, nq, q, is, ie, js, je, x, acc)
      !! `sum x(i,j,q)` over the owned cells of this rank, as `efp_sum_prod` (same footprint, same local EFP result).
      integer, intent(in) :: nx, ny, nq, q, is, ie, js, je
      real(wp), intent(in) :: x(nx, ny, nq)
      type(efp_t), intent(out) :: acc
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
      !$acc&    private(v, d1, d2, d3, d4, d5, d6, dp) present(x)
      do j = js, je
         do i = is, ie
            v = real(x(i, j, q), real64)
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
   end subroutine efp_sum_plane

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

   subroutine fill_snow_reset(nx, ny, ncat, nslyr, is, ie, js, je, puny, hs_min, rhos, Lfresh, aicen, vsnon, qsnon, scr)
      !! Per-cell energy the vertical-profile initialisation adds by resetting the snow thinner than `hs_min` to the
      !! melting enthalpy `-rhos Lfresh` (module docstring): for every category with `aicen > puny`, if the layer
      !! thickness `vsnon/aicen/nslyr` is not above `hs_min/nslyr`, `sum_k (-rhos Lfresh - q_k) vsnon/nslyr`.  ONE launch.
      integer, intent(in) :: nx, ny, ncat, nslyr, is, ie, js, je
      real(wp), intent(in) :: puny, hs_min, rhos, Lfresh
      real(wp), intent(in) :: aicen(nx, ny, ncat), vsnon(nx, ny, ncat), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: scr(nx, ny)
      integer :: i, j, n, k
      real(wp) :: e, rns, hslyr

      rns = real(nslyr, wp)
      !$acc data present(aicen, vsnon, qsnon, scr)
      do concurrent(j=js:je, i=is:ie) local(n, k, e, hslyr)
         e = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            if (aicen(i, j, n) <= puny) cycle
            hslyr = (vsnon(i, j, n)/aicen(i, j, n))/rns
            if (hslyr > hs_min/rns) cycle
            do k = 1, nslyr
               e = e + ((-rhos*Lfresh) - qsnon(i, j, n, k))*vsnon(i, j, n)/rns
            end do
         end do
         scr(i, j) = e
      end do
      !$acc end data
   end subroutine fill_snow_reset

   subroutine state_totals(this, cice, ixm)
      !! Local EFP totals of `M`, `S`, `E` of the state into `lv(ixm:ixm+2)`.
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      integer, intent(in) :: ixm

      call fill_state(this%nx, this%ny, this%ncat, this%nilyr, this%nslyr, this%is, this%ie, this%js, this%je, &
                      cice%p%c%rhoi, cice%p%c%rhos, cice%p%thermo%ice_ref_salinity, &
                      cice%vicen, cice%vsnon, cice%qicen, cice%qsnon, this%scr, this%scrs, this%scre)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(ixm))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scrs, this%lv(ixm + 1))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scre, this%lv(ixm + 2))
   end subroutine state_totals

   ! ------------------------------------------------------------------
   ! the step
   ! ------------------------------------------------------------------

   subroutine tl_totals(this, cice, m, s, e, seam)
      !! The conserved totals `M` [kg], `S` [kg], `E` [J] of the ice state now (EFP sums, GLOBAL when `seam` is given
      !! and has several ranks: one collective, every rank calls it).  Host call; the state and the ledger must be
      !! device-mapped.  Used by the tests and the console.
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(out) :: m, s, e
      type(transport_seam_t), intent(in), optional :: seam
      type(efp_t) :: t(3), keep(3)

      m = 0.0_wp
      s = 0.0_wp
      e = 0.0_wp
      if (.not. this%is_init .or. .not. this%mapped) return
      keep = this%lv(1:3)      ! (the slots of the `begin` totals: an open step is not disturbed)
      call state_totals(this, cice, 1)
      t = this%lv(1:3)
      this%lv(1:3) = keep
      if (present(seam)) call seam_sum_efp(seam, 3, t)
      m = real(efp_to_real(t(1)), wp)
      s = real(efp_to_real(t(2)), wp)
      e = real(efp_to_real(t(3)), wp)
   end subroutine tl_totals

   subroutine tl_begin(this, cice, dt, fsnow, frain, flw, with_sources)
      !! Open a step of length `dt` [s]: snapshot this rank's totals (the fluxes are zero: the sequencer zeroes them
      !! just before) and form the external sources of the step from this step's atmosphere inputs and the
      !! start-of-step concentration `aice_init` (call it after `ice_save_init`).  `fsnow`, `frain` [kg m-2 s-1] and
      !! `flw` [W m-2] are per cell area.  No collective: everything stays local until `end_step`.
      class(ice_term_ledger_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: fsnow(this%nx, this%ny), frain(this%nx, this%ny), flw(this%nx, this%ny)
      logical, intent(in) :: with_sources
         !! `.false.` when the step skips the thermodynamics (no atmosphere exchange: the sources stay zero).

      this%valid = .false.
      this%is_open = .false.
      if (.not. this%is_init .or. .not. this%mapped) return
      this%dt = dt
      this%tab = ice_ledger_tab_t()
      this%lv = efp_t()
      this%probed = .false.
      this%have_clip = .false.
      this%with_sources = with_sources
      this%rhoi = cice%p%c%rhoi
      this%rhos = cice%p%c%rhos
      this%lfresh = cice%p%c%Lfresh
      this%sref = cice%p%thermo%ice_ref_salinity
      call state_totals(this, cice, ice_lt_ix_state(1, 0))
      this%probed(0) = .true.
      ! the start-of-step state of the thin-snow reset and the solver bound (both stay local until `end_step`)
      call fill_snow_reset(this%nx, this%ny, this%ncat, this%nslyr, this%is, this%ie, this%js, this%je, cice%p%c%puny, &
                           cice%p%c%hs_min, cice%p%c%rhos, cice%p%c%Lfresh, cice%aicen, cice%vsnon, cice%qsnon, this%scr)
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(IX_SNR))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, cice%aice_init, this%lv(IX_BND))
      if (with_sources) then
         call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, fsnow, cice%aice_init, this%scr)
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(IX_SRC + 1))
         call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, frain, cice%aice_init, this%scr)
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(IX_SRC + 2))
         call fill_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, flw, cice%aice_init, this%scr)
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(IX_SRC + 3))
      end if
      this%is_open = .true.
   end subroutine tl_begin

   subroutine tl_probe(this, stage, cice, flux)
      !! Close stage `stage` (`ICE_LT_*`): keep this rank's totals and flux accumulators after it (the differences
      !! with the previous probe are formed by `end_step`, on the global values).  After `therm1` the atmosphere
      !! terms are read; after `newice` the new-ice volume.  No collective.
      class(ice_term_ledger_t), intent(inout) :: this
      integer, intent(in) :: stage
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_flux_t), intent(in) :: flux

      if (.not. this%is_open) return
      call state_totals(this, cice, ice_lt_ix_state(1, stage))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fresh, &
                        this%lv(ice_lt_ix_flux(1, stage)))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fsalt, &
                        this%lv(ice_lt_ix_flux(2, stage)))
      call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fhocn, &
                        this%lv(ice_lt_ix_flux(3, stage)))
      this%probed(stage) = .true.
      if (stage == ICE_LT_THERM1) then
         ! only the `therm1` stage writes `fswthru`: the other probes carry the previous value
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%fswthru, this%lv(IX_SW))
         call fill_atm(this%nx, this%ny, this%is, this%ie, this%js, this%je, cice%p%c%Lvap, flux%fswabs, &
                       flux%flwout, flux%fsens, flux%evap, this%scr, this%scrs)
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scr, this%lv(IX_ATM + 1))
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, this%scrs, this%lv(IX_ATM + 2))
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%evap, this%lv(IX_ATM + 3))
      else if (stage == ICE_LT_NEWICE) then
         call efp_sum_prod(this%nx, this%ny, this%is, this%ie, this%js, this%je, this%area, flux%frazil_vol, this%lv(IX_FRZ))
      end if
   end subroutine tl_probe

   subroutine tl_book_clip(this, clip)
      !! Keep this rank's transport clip record of the step: the per-cell planes of the remap
      !! (`remap_flux_workspace_t%clip`: clipped count, area created, ice and snow volume and enthalpy dropped),
      !! summed over the owned cells with the EFP sums (six reductions, no collective).  Call it after the
      !! transport and before the `transport` probe; `end_step` combines it with everything else and books it as a
      !! source of the `transport` stage (module docstring).  A step whose transport does not run the remap books none.
      class(ice_term_ledger_t), intent(inout) :: this
      real(wp), intent(in) :: clip(this%nx, this%ny, ICE_LT_NCLIP)
      integer :: k

      if (.not. this%is_open) return
      do k = 1, ICE_LT_NCLIP
         call efp_sum_plane(this%nx, this%ny, ICE_LT_NCLIP, k, this%is, this%ie, this%js, this%je, clip, &
                            this%lv(ice_lt_ix_clip(k)))
      end do
      this%have_clip = .true.
   end subroutine tl_book_clip

   subroutine tl_end(this, seam)
      !! Close the step: combine the step's local EFP values over the ranks in ONE collective (`seam`, default one
      !! rank), form the differences of the totals between probes, the booked sources and the residuals of every
      !! stage and line on the global values, add them to the cumulative accumulators.
      class(ice_term_ledger_t), intent(inout) :: this
      type(transport_seam_t), intent(in), optional :: seam
      integer :: s, q, k
      real(wp) :: r, dt
      real(wp) :: prev(3), fprev(4), cur(3), f(4)
      real(wp) :: srcw(ICE_LT_NSTAGE), srch(ICE_LT_NSTAGE), srcs(ICE_LT_NSTAGE), worst_step(3), g(ICE_LT_NV)
      logical :: fin
      type(efp_t) :: keep

      this%valid = .false.
      if (.not. this%is_open) return
      this%is_open = .false.
      if (ice_lt_test_local > 0) keep = this%lv(ice_lt_test_local)
      if (present(seam)) call seam_sum_efp(seam, ICE_LT_NV, this%lv)
      if (ice_lt_test_local > 0) this%lv(ice_lt_test_local) = keep
      do k = 1, ICE_LT_NV
         g(k) = real(efp_to_real(this%lv(k)), wp)
      end do
      dt = this%dt
      associate (t => this%tab)
         ! the changes of the totals and of the flux accumulators between the probes that ran
         do q = 1, 3
            prev(q) = g(ice_lt_ix_state(q, 0))
         end do
         fprev = 0.0_wp
         do s = 1, ICE_LT_NSTAGE
            if (.not. this%probed(s)) cycle
            do q = 1, 3
               cur(q) = g(ice_lt_ix_state(q, s))
               t%dstate(q, s) = cur(q) - prev(q)
               f(q) = g(ice_lt_ix_flux(q, s))
               t%fout(q, s) = dt*(f(q) - fprev(q))
            end do
            f(4) = fprev(4)
            if (s == ICE_LT_THERM1) f(4) = g(IX_SW)
            t%fsw(s) = dt*(f(4) - fprev(4))
            prev = cur
            fprev = f
         end do
         if (this%probed(ICE_LT_THERM1)) then
            t%atm_rad = dt*g(IX_ATM + 1)
            t%atm_lat = dt*g(IX_ATM + 2)
            t%evap = dt*g(IX_ATM + 3)
         end if
         if (this%probed(ICE_LT_NEWICE)) t%frz_lat = -this%rhoi*this%lfresh*g(IX_FRZ)
         if (this%with_sources) then
            t%snow = dt*g(IX_SRC + 1)
            t%snow_enth = this%lfresh*dt*g(IX_SRC + 1)
            t%rain = dt*g(IX_SRC + 2)
            t%lw_in = dt*g(IX_SRC + 3)
         end if
         t%snow_reset = g(IX_SNR)
         t%heat_bound = 0.9_wp*ICE_FERRMAX*dt*g(IX_BND)
         srcw = 0.0_wp
         srch = 0.0_wp
         srcs = 0.0_wp
         if (this%have_clip) then
            t%n_clip = g(ice_lt_ix_clip(1))
            t%clip_area = g(ice_lt_ix_clip(2))
            t%drop_vice = g(ice_lt_ix_clip(3))
            t%drop_vsno = g(ice_lt_ix_clip(4))
            t%drop_eice = g(ice_lt_ix_clip(5))
            t%drop_esno = g(ice_lt_ix_clip(6))
            t%clip_w = -(this%rhoi*t%drop_vice + this%rhos*t%drop_vsno)
            t%clip_s = -(this%rhoi*t%drop_vice*this%sref*1.0e-3_wp)
            t%clip_h = -(t%drop_eice + t%drop_esno - this%rhos*this%lfresh*t%drop_vsno)
            srcw(ICE_LT_TRANSPORT) = t%clip_w
            srcs(ICE_LT_TRANSPORT) = t%clip_s
            srch(ICE_LT_TRANSPORT) = t%clip_h
         end if
         srcw(ICE_LT_THERM1) = t%snow + t%evap
         srcw(ICE_LT_REMAP) = t%rain
         srch(ICE_LT_THERM1) = t%lw_in - t%snow_enth + t%atm_rad + t%atm_lat
         srch(ICE_LT_NEWICE) = t%frz_lat
         t%scale = 0.0_wp
         do s = 1, ICE_LT_NSTAGE
            t%res(ICE_LT_WATER, s) = t%dstate(1, s) - (srcw(s) - t%fout(1, s))
            t%res(ICE_LT_SALT, s) = (t%dstate(2, s) + t%fout(2, s)) - srcs(s)
            t%res(ICE_LT_HEAT, s) = t%dstate(3, s) - (srch(s) - t%fout(3, s) - t%fsw(s))
            t%scale(ICE_LT_WATER) = t%scale(ICE_LT_WATER) + abs(t%dstate(1, s)) + abs(srcw(s)) + abs(t%fout(1, s))
            t%scale(ICE_LT_SALT) = t%scale(ICE_LT_SALT) + abs(t%dstate(2, s)) + abs(t%fout(2, s)) + abs(srcs(s))
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
            this%cum_clip = ieee_value(1.0_wp, ieee_quiet_nan)
            this%cum_zap = ieee_value(1.0_wp, ieee_quiet_nan)
            this%cum_snow_reset = ieee_value(1.0_wp, ieee_quiet_nan)
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
            this%cum_zap(q) = this%cum_zap(q) + (t%fout(q, ICE_LT_CLEANUP) + t%fout(q, ICE_LT_TRANSPORT))
         end do
         this%cum_clip = this%cum_clip + [t%n_clip, t%clip_area, t%drop_vice, t%drop_vsno, t%drop_eice, t%drop_esno]
         this%cum_snow_reset = this%cum_snow_reset + t%snow_reset
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
