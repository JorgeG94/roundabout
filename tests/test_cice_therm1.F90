!! The vertical thermodynamics step of the CICE-port ice (D3.4): the Stefan growth law
!! through the whole glue, the coupled 1000-step column ledger (water, salt, heat), a
!! melt season with shortwave, empty categories left untouched bitwise, error
!! isolation (one flagged category does not touch another category or cell), the merge
!! against the hand sum, goldens from the real upstream step, and a device sweep of every
!! output of every cell and category.
!!
!! **Stefan (A7).**  Fresh ice (`l_brine` off, zero salinity, Maykut-Untersteiner
!! conductivity `kice`) under a fixed surface temperature `Ts` (the BL99 test surface
!! `ICE_SFC_LINEAR` with a stiff derivative), freezing point 0 degC, no snow, no ocean
!! heat.  From the steady linear profile the first step is exactly
!! `dhi = k |Ts| dt/(h0 rhoi Lfresh)`; the long-time law is
!! `h^2 = h0^2 + 2 k |Ts| t/(rhoi Lfresh)` (Stefan 1891; the zero-layer limit of
!! Maykut and Untersteiner 1971), whose error against the four-layer scheme (the
!! sensible heat of the new ice, a concave profile) is measured and pinned below.
!!
!! **Ledger.**  Per step and per cell, with `M = rhoi vice + rhos vsno`, `S = rhoi vice
!! Sref 1e-3`, `E = sum_n (sum_k qicen vicen/nilyr + qsnon vsnon/nslyr)`:
!! water `dM = dt (fsnow aice_init + evap - fresh)`, salt `dS = -dt fsalt`, heat `dE = dt
!! (fswabs - fswthru + flwout + fsens - fhocn + aice_init flw - Lfresh fsnow aice_init)`
!! (`flat` leaves with the vapour; `fswthru = sum aicen_init fswthrun` is the radiation
!! driver's).  Brine (BL99) conserves energy to the check bound `0.9 ferrmax = 9e-4 W/m2`
!! per category, not to round-off; the residual is reported relative to the sum of the
!! absolute values of every term.
!!
!! **Goldens.**  `icepack_step_therm1` of upstream (Icepack 9a1f3e4, gfortran -O0
!! -ffp-contract=off, both constant sets) on 24 cells x 5 categories (nilyr 7 and 4,
!! highfreq on and off, `sw_redist` on and off, growth, melt, snow, snow-ice,
!! sublimation, condensation, empty and tiny categories): the post-step state and
!! `fsensn, flatn, melttn, meltbn, congeln, snoicen` per category and the 14 merged
!! cell fluxes, the stress and the basal heat.  Inputs and outputs are literals.  What
!! differs from upstream is the reproducible transcendentals of the atmosphere layer and
!! the surface balance (our `rdb_sm_exp/log/atan` against libm); the measured worst
!! relative differences are in the test log line and in the PR report.
module test_cice_therm1
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_is_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_THERMO_CONSERV
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_CONDUCT_MU71, ice_q_from_T, ice_qsno_from_T, ice_salinity_profile
   use rdb_ice_basal_heat, only: ice_basal_heat_step, ice_frzmlt_bottom_lateral_cell
   use rdb_ice_flux, only: ice_flux_merge_impl
   use rdb_ice_thermo_shared, only: ICE_SFC_BULK, ICE_SFC_LINEAR
   use rdb_ice_therm1, only: ice_step_therm1, ice_therm1_impl, ice_therm1_cell, ice_therm1_merge_impl
   use rdb_ice_bl99, only: ICE_NITERMAX, ICE_SW_FRAC, ICE_SW_DTEMP
   use rdb_ice_thickness, only: ICE_FERRMAX
   implicit none
   private

   public :: collect_cice_therm1_tests

   real(wp), parameter :: TWOPI = 6.283185307179586_wp

#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: DRTOL = 1.0e-12_wp
      !! Device-vs-host (and cell-vs-column) relative tolerance, GPU build (FMA on, host code too).
   real(wp), parameter :: SUMTOL = 1.0e-13_wp
      !! Merge against the hand sum, GPU build.
   real(wp), parameter :: GRTOL = 1.0e-12_wp
      !! Golden relative tolerance, GPU build.
#else
   real(wp), parameter :: DRTOL = 0.0_wp
   real(wp), parameter :: SUMTOL = 0.0_wp
   real(wp), parameter :: GRTOL = 4.0e-14_wp
#endif

   type :: world_t
      !! One test world: the arrays of `ice_step_therm1` at the ghost-less extent `(nx, ny, ncat, nilyr)`.
      integer :: nx = 0, ny = 0, ncat = 0, nilyr = 0, nslyr = 1, is = 1, ie = 0, js = 1, je = 0
      type(ice_params_t) :: p
      logical :: sw_redist = .false.
      integer :: sfc_kind = ICE_SFC_BULK
      logical :: manual = .false.   !! run the four launches by hand (needed for `merge_puny`)
      real(wp) :: sf_0 = 0.0_wp, dsf_dT = 0.0_wp, dt = 900.0_wp, zlvl = 10.0_wp, zlvs = 10.0_wp
      real(wp) :: merge_puny = -1.0_wp   !! puny of the manual-path merge (the constants' value unless a test sets it)
      real(wp), allocatable :: uatm(:, :), vatm(:, :), wind(:, :), potT(:, :), Qa(:, :), rhoa(:, :), &
                               flw(:, :), fsnow(:, :), uvel(:, :), vvel(:, :), Tf(:, :), sst(:, :), frzmlt(:, :), &
                               tau_x(:, :), tau_y(:, :), fresh(:, :), fsalt(:, :), fhocn(:, :), evap(:, :), fsens(:, :), &
                               flat(:, :), flwout(:, :), fswabs(:, :), meltt(:, :), meltb(:, :), congel(:, :), snoice(:, :), &
                               strair_x_T(:, :), strair_y_T(:, :), fbot(:, :), rside(:, :), wlat(:, :)
      integer(int32), allocatable :: ice_err(:, :)
      real(wp), allocatable :: aicen_init(:, :, :), aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), &
                               Tsfcn(:, :, :), fswsfcn(:, :, :), fswintn(:, :, :), fswthrun(:, :, :), freshn(:, :, :), &
                               fsaltn(:, :, :), fhocnn(:, :, :), evapn(:, :, :), fsensn(:, :, :), flatn(:, :, :), &
                               flwoutn(:, :, :), fswabsn(:, :, :), melttn(:, :, :), meltbn(:, :, :), congeln(:, :, :), &
                               snoicen(:, :, :), strairxn(:, :, :), strairyn(:, :, :)
      integer(int32), allocatable :: ierr_n(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), Iswabsn(:, :, :, :), &
                               qsnon(:, :, :, :), Sswabsn(:, :, :, :)
   end type world_t

contains

   subroutine collect_cice_therm1_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("stefan_first_step_closed_form", test_stefan_first), &
                  new_unittest("stefan_sqrt_t_law", test_stefan_law), &
                  new_unittest("coupled_column_ledger_1000_steps", test_ledger), &
                  new_unittest("melt_season_with_shortwave", test_melt), &
                  new_unittest("empty_category_untouched_bitwise", test_empty), &
                  new_unittest("error_isolation_cell_and_category", test_errors), &
                  new_unittest("merge_equals_hand_sum_and_is_additive", test_merge), &
                  new_unittest("reflected_longwave_only_for_active_categories", test_flwout_gate), &
                  new_unittest("golden_vs_upstream_step_therm1", test_golden), &
                  new_unittest("device_sweep_every_output", test_device) &
                  ]
   end subroutine collect_cice_therm1_tests

   ! ------------------------------------------------------------------
   ! the world
   ! ------------------------------------------------------------------

   subroutine world_alloc(w, nx, ny, ncat, nilyr, set)
      type(world_t), intent(inout) :: w
      integer, intent(in) :: nx, ny, ncat, nilyr, set

      w%nx = nx
      w%ny = ny
      w%ncat = ncat
      w%nilyr = nilyr
      w%nslyr = 1
      w%is = 1
      w%ie = nx
      w%js = 1
      w%je = ny
      w%p = ice_params_default(set, ncat, nilyr, 1)
      w%merge_puny = w%p%c%puny
      if (allocated(w%uatm)) deallocate (w%uatm)
      if (allocated(w%vatm)) deallocate (w%vatm)
      if (allocated(w%wind)) deallocate (w%wind)
      if (allocated(w%potT)) deallocate (w%potT)
      if (allocated(w%Qa)) deallocate (w%Qa)
      if (allocated(w%rhoa)) deallocate (w%rhoa)
      if (allocated(w%flw)) deallocate (w%flw)
      if (allocated(w%fsnow)) deallocate (w%fsnow)
      if (allocated(w%uvel)) deallocate (w%uvel)
      if (allocated(w%vvel)) deallocate (w%vvel)
      if (allocated(w%Tf)) deallocate (w%Tf)
      if (allocated(w%sst)) deallocate (w%sst)
      if (allocated(w%frzmlt)) deallocate (w%frzmlt)
      if (allocated(w%tau_x)) deallocate (w%tau_x)
      if (allocated(w%tau_y)) deallocate (w%tau_y)
      if (allocated(w%fresh)) deallocate (w%fresh)
      if (allocated(w%fsalt)) deallocate (w%fsalt)
      if (allocated(w%fhocn)) deallocate (w%fhocn)
      if (allocated(w%evap)) deallocate (w%evap)
      if (allocated(w%fsens)) deallocate (w%fsens)
      if (allocated(w%flat)) deallocate (w%flat)
      if (allocated(w%flwout)) deallocate (w%flwout)
      if (allocated(w%fswabs)) deallocate (w%fswabs)
      if (allocated(w%meltt)) deallocate (w%meltt)
      if (allocated(w%meltb)) deallocate (w%meltb)
      if (allocated(w%congel)) deallocate (w%congel)
      if (allocated(w%snoice)) deallocate (w%snoice)
      if (allocated(w%strair_x_T)) deallocate (w%strair_x_T)
      if (allocated(w%strair_y_T)) deallocate (w%strair_y_T)
      if (allocated(w%fbot)) deallocate (w%fbot)
      if (allocated(w%rside)) deallocate (w%rside)
      if (allocated(w%wlat)) deallocate (w%wlat)
      if (allocated(w%ice_err)) deallocate (w%ice_err)
      if (allocated(w%aicen_init)) deallocate (w%aicen_init)
      if (allocated(w%aicen)) deallocate (w%aicen)
      if (allocated(w%vicen)) deallocate (w%vicen)
      if (allocated(w%vsnon)) deallocate (w%vsnon)
      if (allocated(w%Tsfcn)) deallocate (w%Tsfcn)
      if (allocated(w%fswsfcn)) deallocate (w%fswsfcn)
      if (allocated(w%fswintn)) deallocate (w%fswintn)
      if (allocated(w%fswthrun)) deallocate (w%fswthrun)
      if (allocated(w%freshn)) deallocate (w%freshn)
      if (allocated(w%fsaltn)) deallocate (w%fsaltn)
      if (allocated(w%fhocnn)) deallocate (w%fhocnn)
      if (allocated(w%evapn)) deallocate (w%evapn)
      if (allocated(w%fsensn)) deallocate (w%fsensn)
      if (allocated(w%flatn)) deallocate (w%flatn)
      if (allocated(w%flwoutn)) deallocate (w%flwoutn)
      if (allocated(w%fswabsn)) deallocate (w%fswabsn)
      if (allocated(w%melttn)) deallocate (w%melttn)
      if (allocated(w%meltbn)) deallocate (w%meltbn)
      if (allocated(w%congeln)) deallocate (w%congeln)
      if (allocated(w%snoicen)) deallocate (w%snoicen)
      if (allocated(w%strairxn)) deallocate (w%strairxn)
      if (allocated(w%strairyn)) deallocate (w%strairyn)
      if (allocated(w%ierr_n)) deallocate (w%ierr_n)
      if (allocated(w%qicen)) deallocate (w%qicen)
      if (allocated(w%sicen)) deallocate (w%sicen)
      if (allocated(w%Iswabsn)) deallocate (w%Iswabsn)
      if (allocated(w%qsnon)) deallocate (w%qsnon)
      if (allocated(w%Sswabsn)) deallocate (w%Sswabsn)
      allocate (w%uatm(nx, ny))
      w%uatm = 0.0_wp
      allocate (w%vatm(nx, ny))
      w%vatm = 0.0_wp
      allocate (w%wind(nx, ny))
      w%wind = 0.0_wp
      allocate (w%potT(nx, ny))
      w%potT = 0.0_wp
      allocate (w%Qa(nx, ny))
      w%Qa = 0.0_wp
      allocate (w%rhoa(nx, ny))
      w%rhoa = 0.0_wp
      allocate (w%flw(nx, ny))
      w%flw = 0.0_wp
      allocate (w%fsnow(nx, ny))
      w%fsnow = 0.0_wp
      allocate (w%uvel(nx, ny))
      w%uvel = 0.0_wp
      allocate (w%vvel(nx, ny))
      w%vvel = 0.0_wp
      allocate (w%Tf(nx, ny))
      w%Tf = 0.0_wp
      allocate (w%sst(nx, ny))
      w%sst = 0.0_wp
      allocate (w%frzmlt(nx, ny))
      w%frzmlt = 0.0_wp
      allocate (w%tau_x(nx, ny))
      w%tau_x = 0.0_wp
      allocate (w%tau_y(nx, ny))
      w%tau_y = 0.0_wp
      allocate (w%fresh(nx, ny))
      w%fresh = 0.0_wp
      allocate (w%fsalt(nx, ny))
      w%fsalt = 0.0_wp
      allocate (w%fhocn(nx, ny))
      w%fhocn = 0.0_wp
      allocate (w%evap(nx, ny))
      w%evap = 0.0_wp
      allocate (w%fsens(nx, ny))
      w%fsens = 0.0_wp
      allocate (w%flat(nx, ny))
      w%flat = 0.0_wp
      allocate (w%flwout(nx, ny))
      w%flwout = 0.0_wp
      allocate (w%fswabs(nx, ny))
      w%fswabs = 0.0_wp
      allocate (w%meltt(nx, ny))
      w%meltt = 0.0_wp
      allocate (w%meltb(nx, ny))
      w%meltb = 0.0_wp
      allocate (w%congel(nx, ny))
      w%congel = 0.0_wp
      allocate (w%snoice(nx, ny))
      w%snoice = 0.0_wp
      allocate (w%strair_x_T(nx, ny))
      w%strair_x_T = 0.0_wp
      allocate (w%strair_y_T(nx, ny))
      w%strair_y_T = 0.0_wp
      allocate (w%fbot(nx, ny))
      w%fbot = 0.0_wp
      allocate (w%rside(nx, ny))
      w%rside = 0.0_wp
      allocate (w%wlat(nx, ny))
      w%wlat = 0.0_wp
      allocate (w%ice_err(nx, ny))
      w%ice_err = 0_int32
      allocate (w%aicen_init(nx, ny, ncat))
      w%aicen_init = 0.0_wp
      allocate (w%aicen(nx, ny, ncat))
      w%aicen = 0.0_wp
      allocate (w%vicen(nx, ny, ncat))
      w%vicen = 0.0_wp
      allocate (w%vsnon(nx, ny, ncat))
      w%vsnon = 0.0_wp
      allocate (w%Tsfcn(nx, ny, ncat))
      w%Tsfcn = 0.0_wp
      allocate (w%fswsfcn(nx, ny, ncat))
      w%fswsfcn = 0.0_wp
      allocate (w%fswintn(nx, ny, ncat))
      w%fswintn = 0.0_wp
      allocate (w%fswthrun(nx, ny, ncat))
      w%fswthrun = 0.0_wp
      allocate (w%freshn(nx, ny, ncat))
      w%freshn = 0.0_wp
      allocate (w%fsaltn(nx, ny, ncat))
      w%fsaltn = 0.0_wp
      allocate (w%fhocnn(nx, ny, ncat))
      w%fhocnn = 0.0_wp
      allocate (w%evapn(nx, ny, ncat))
      w%evapn = 0.0_wp
      allocate (w%fsensn(nx, ny, ncat))
      w%fsensn = 0.0_wp
      allocate (w%flatn(nx, ny, ncat))
      w%flatn = 0.0_wp
      allocate (w%flwoutn(nx, ny, ncat))
      w%flwoutn = 0.0_wp
      allocate (w%fswabsn(nx, ny, ncat))
      w%fswabsn = 0.0_wp
      allocate (w%melttn(nx, ny, ncat))
      w%melttn = 0.0_wp
      allocate (w%meltbn(nx, ny, ncat))
      w%meltbn = 0.0_wp
      allocate (w%congeln(nx, ny, ncat))
      w%congeln = 0.0_wp
      allocate (w%snoicen(nx, ny, ncat))
      w%snoicen = 0.0_wp
      allocate (w%strairxn(nx, ny, ncat))
      w%strairxn = 0.0_wp
      allocate (w%strairyn(nx, ny, ncat))
      w%strairyn = 0.0_wp
      allocate (w%ierr_n(nx, ny, ncat))
      w%ierr_n = 0_int32
      allocate (w%qicen(nx, ny, ncat, nilyr))
      w%qicen = 0.0_wp
      allocate (w%sicen(nx, ny, ncat, nilyr))
      w%sicen = 0.0_wp
      allocate (w%Iswabsn(nx, ny, ncat, nilyr))
      w%Iswabsn = 0.0_wp
      allocate (w%qsnon(nx, ny, ncat, 1))
      w%qsnon = 0.0_wp
      allocate (w%Sswabsn(nx, ny, ncat, 1))
      w%Sswabsn = 0.0_wp
      call world_forcing_default(w)
   end subroutine world_alloc

   subroutine world_forcing_default(w)
      type(world_t), intent(inout) :: w
      w%potT = 263.0_wp
      w%Qa = 6.0e-4_wp
      w%rhoa = 1.3_wp
      w%wind = 6.0_wp
      w%uatm = 6.0_wp
      w%vatm = 0.0_wp
      w%flw = 200.0_wp
      w%Tf = -1.8_wp
      w%sst = -1.8_wp
      w%tau_x = 0.05_wp
      w%tau_y = 0.0_wp
   end subroutine world_forcing_default

   subroutine world_zero_fluxes(w)
      !! The per-step zeroing (`ice_flux_zero`) of the cell accumulators.
      type(world_t), intent(inout) :: w
      w%fresh = 0.0_wp; w%fsalt = 0.0_wp; w%fhocn = 0.0_wp; w%evap = 0.0_wp; w%fsens = 0.0_wp; w%flat = 0.0_wp
      w%flwout = 0.0_wp; w%fswabs = 0.0_wp; w%meltt = 0.0_wp; w%meltb = 0.0_wp; w%congel = 0.0_wp; w%snoice = 0.0_wp
   end subroutine world_zero_fluxes

   subroutine world_map(w)
      type(world_t), intent(inout) :: w
      !$acc enter data copyin(w%p, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, &
      !$acc&  w%vvel, w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, &
      !$acc&  w%flat, w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%strair_x_T, w%strair_y_T, &
      !$acc&  w%fbot, w%rside, w%wlat, w%ice_err, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, &
      !$acc&  w%fswsfcn, w%fswintn, w%fswthrun, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, &
      !$acc&  w%flwoutn, w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, &
      !$acc&  w%ierr_n, w%qicen, w%sicen, w%Iswabsn, w%qsnon, w%Sswabsn)
   end subroutine world_map

   subroutine world_unmap(w)
      type(world_t), intent(inout) :: w
      !$acc exit data delete(w%p, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, &
      !$acc&  w%vvel, w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, &
      !$acc&  w%flat, w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%strair_x_T, w%strair_y_T, &
      !$acc&  w%fbot, w%rside, w%wlat, w%ice_err, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, &
      !$acc&  w%fswsfcn, w%fswintn, w%fswthrun, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, &
      !$acc&  w%flwoutn, w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, &
      !$acc&  w%ierr_n, w%qicen, w%sicen, w%Iswabsn, w%qsnon, w%Sswabsn)
   end subroutine world_unmap

   subroutine world_step(w)
      !! One ice step on the (mapped) world: push, run the step, pull.  Production path `ice_step_therm1`; with
      !! the BL99 test hook (`sfc_kind /= ICE_SFC_BULK`) the same four launches by hand.
      type(world_t), intent(inout) :: w
      integer :: istat
      !$acc update device(w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel, &
      !$acc&  w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, &
      !$acc&  w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%strair_x_T, w%strair_y_T, w%fbot, &
      !$acc&  w%rside, w%wlat, w%ice_err, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%fswsfcn, &
      !$acc&  w%fswintn, w%fswthrun, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, &
      !$acc&  w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, w%ierr_n, &
      !$acc&  w%qicen, w%sicen, w%Iswabsn, w%qsnon, w%Sswabsn)
      istat = 0
      if (w%sfc_kind == ICE_SFC_BULK .and. .not. w%manual) then
         call ice_step_therm1(w%nx, w%ny, w%ncat, w%nilyr, w%nslyr, w%is, w%ie, w%js, w%je, w%dt, w%p, w%sw_redist, &
                              w%zlvl, w%zlvs, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, &
                              w%vvel, w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, &
                              w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                              w%fswsfcn, w%fswintn, w%fswthrun, w%Iswabsn, w%Sswabsn, &
                              w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, w%meltt, &
                              w%meltb, w%congel, w%snoice, &
                              w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, &
                              w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, &
                              w%strair_x_T, w%strair_y_T, w%fbot, w%rside, w%wlat, w%ierr_n, w%ice_err, istat)
      else
         call ice_basal_heat_step(w%nx, w%ny, w%ncat, w%nilyr, w%nslyr, w%is, w%ie, w%js, w%je, w%dt, w%p, &
                                  w%aicen, w%vicen, w%vsnon, w%qicen, w%qsnon, w%frzmlt, w%sst, w%Tf, w%tau_x, w%tau_y, &
                                  w%fbot, w%rside, w%wlat, w%ice_err)
         call ice_therm1_impl(w%nx, w%ny, w%ncat, w%nilyr, w%nslyr, w%is, w%ie, w%js, w%je, w%dt, w%p%c, w%p%atm, &
                              w%p%thermo%conduct, w%p%thermo%l_brine, w%p%thermo%emissivity, &
                              w%p%thermo%ice_ref_salinity, w%sw_redist, ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, &
                              w%sfc_kind, w%sf_0, w%dsf_dT, w%zlvl, w%zlvs, &
                              w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel, w%Tf, &
                              w%fbot, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                              w%fswsfcn, w%fswintn, w%fswthrun, w%Iswabsn, w%Sswabsn, &
                              w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, &
                              w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, w%ice_err, w%ierr_n)
         call ice_flux_merge_impl(w%nx, w%ny, w%ncat, w%is, w%ie, w%js, w%je, w%p%thermo%emissivity, w%merge_puny, w%flw, &
                                  w%aicen_init, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, &
                                  w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, &
                                  w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, &
                                  w%meltt, w%meltb, w%congel, w%snoice, w%ice_err)
         call ice_therm1_merge_impl(w%nx, w%ny, w%ncat, w%is, w%ie, w%js, w%je, w%aicen_init, w%strairxn, &
                                    w%strairyn, w%ierr_n, w%strair_x_T, w%strair_y_T, w%ice_err)
      end if
      !$acc update self(w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel, w%Tf, &
      !$acc&  w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, &
      !$acc&  w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%strair_x_T, w%strair_y_T, w%fbot, &
      !$acc&  w%rside, w%wlat, w%ice_err, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%fswsfcn, &
      !$acc&  w%fswintn, w%fswthrun, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, &
      !$acc&  w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, w%ierr_n, &
      !$acc&  w%qicen, w%sicen, w%Iswabsn, w%qsnon, w%Sswabsn)
   end subroutine world_step

   subroutine set_column(w, i, j, n, a, hin, hsn, Tsfc, Ttop, Tbase, fresh_ice)
      !! One category: area `a`, ice `hin`, snow `hsn` [m], surface `Tsfc`, a linear temperature profile from `Ttop` (top
      !! layer centre) to `Tbase` (base layer centre) [degC], the BL99 salinity profile (zero for `fresh_ice`).
      type(world_t), intent(inout) :: w
      integer, intent(in) :: i, j, n
      real(wp), intent(in) :: a, hin, hsn, Tsfc, Ttop, Tbase
      logical, intent(in) :: fresh_ice
      integer :: k
      real(wp) :: S(w%nilyr + 1), T, fr

      call ice_salinity_profile(w%nilyr, w%p%c%saltmax, w%p%c%min_salin, S)
      w%aicen(i, j, n) = a
      w%aicen_init(i, j, n) = a
      w%vicen(i, j, n) = a*hin
      w%vsnon(i, j, n) = a*hsn
      w%Tsfcn(i, j, n) = Tsfc
      do k = 1, w%nilyr
         fr = merge(0.0_wp, real(k - 1, wp)/real(w%nilyr - 1, wp), w%nilyr == 1)
         T = Ttop + (Tbase - Ttop)*fr
         if (fresh_ice) then
            w%sicen(i, j, n, k) = 0.0_wp
         else
            w%sicen(i, j, n, k) = S(k)
         end if
         w%qicen(i, j, n, k) = ice_q_from_T(T, w%sicen(i, j, n, k), w%p%c, w%p%thermo%l_brine)
      end do
      w%qsnon(i, j, n, 1) = ice_qsno_from_T(min(0.5_wp*Tsfc, -0.01_wp), w%p%c)
   end subroutine set_column

   function near(a, b, rtol, atol) result(ok)
      real(wp), intent(in) :: a, b, rtol, atol
      logical :: ok
      ok = .false.
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) return
      ok = abs(a - b) <= atol + rtol*max(abs(a), abs(b))
   end function near

   function uniform(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function uniform

   ! ------------------------------------------------------------------
   ! Stefan
   ! ------------------------------------------------------------------

   subroutine stefan_world(w, h0, Ts, dt)
      !! One fresh-ice category on the steady linear profile for the fixed surface temperature `Ts`.
      type(world_t), intent(inout) :: w
      real(wp), intent(in) :: h0, Ts, dt
      integer :: k
      real(wp) :: z, T

      call world_alloc(w, 1, 1, 1, 4, ICE_CONSTANTS_ICEPACK)
      w%p%thermo%conduct = ICE_CONDUCT_MU71
      w%p%thermo%l_brine = .false.
      w%p%c%min_salin = 0.0_wp   ! fresh ice: zero salinity must pass the profile's range check
      w%dt = dt
      w%sfc_kind = ICE_SFC_LINEAR
      w%dsf_dT = 1.0e6_wp
      w%sf_0 = -w%dsf_dT*Ts
      w%Tf = 0.0_wp
      w%sst = 0.0_wp
      w%tau_x = 0.0_wp
      w%aicen(1, 1, 1) = 1.0_wp
      w%aicen_init(1, 1, 1) = 1.0_wp
      w%vicen(1, 1, 1) = h0
      w%Tsfcn(1, 1, 1) = Ts
      do k = 1, w%nilyr
         z = (real(k, wp) - 0.5_wp)*h0/real(w%nilyr, wp)
         T = Ts*(1.0_wp - z/h0)
         w%sicen(1, 1, 1, k) = 0.0_wp
         w%qicen(1, 1, 1, k) = -w%p%c%rhoi*(w%p%c%Lfresh - w%p%c%cp_ice*T)
      end do
      w%qsnon(1, 1, 1, 1) = -w%p%c%rhos*w%p%c%Lfresh
   end subroutine stefan_world

   subroutine test_stefan_first(error)
      !! First step from the steady profile: `dhi = k |Ts| dt/(h0 rhoi Lfresh)` (the base heat flux is `k |Ts|/h0` exactly
      !! for a linear profile and a fixed surface), `congel = dhi`, `fsalt = -rhoi dhi Sref 1e-3/dt`, `fresh = -rhoi
      !! dhi/dt`, `fhocn = 0`, and the merge with the area 1 weight.  The surface sits `k|Ts|/(h0 dsf_dT)` from `Ts`.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp), parameter :: H0 = 1.0_wp, TS = -20.0_wp, DT = 3600.0_wp
      real(wp) :: dhi, dhi_got, k

      call stefan_world(w, H0, TS, DT)
      call world_map(w)
      call world_step(w)
      call world_unmap(w)
      k = w%p%c%kice
      dhi = k*abs(TS)*DT/(H0*w%p%c%rhoi*w%p%c%Lfresh)
      dhi_got = w%vicen(1, 1, 1) - H0
      if (w%ice_err(1, 1) /= 0_int32) then
         write (output_unit, '(a,i0,a,i0)') "ice_err = ", w%ice_err(1, 1), " cat err ", w%ierr_n(1, 1, 1)
         call test_failed(error, "no error bit")
         return
      end if
      call check(error, near(dhi_got, dhi, 2.0e-5_wp, 0.0_wp), "first-step growth is the closed form")
      if (allocated(error)) return
      call check(error, near(w%congeln(1, 1, 1), dhi_got, 1.0e-12_wp, 0.0_wp), "congelation is the thickness gain")
      if (allocated(error)) return
      call check(error, near(w%fsaltn(1, 1, 1), -w%p%c%rhoi*dhi_got*w%p%thermo%ice_ref_salinity*0.001_wp/DT, 1.0e-12_wp, 0.0_wp), &
                 "salt flux books the total thickness change")
      if (allocated(error)) return
      call check(error, near(w%freshn(1, 1, 1), -w%p%c%rhoi*dhi_got/DT, 1.0e-12_wp, 0.0_wp), "fresh water books the ice mass")
      if (allocated(error)) return
      call check(error, abs(w%fhocnn(1, 1, 1)) < 1.0e-9_wp, "no heat to the ocean when only growing")
      if (allocated(error)) return
      call check(error, near(w%fsalt(1, 1), w%fsaltn(1, 1, 1), 1.0e-14_wp, 0.0_wp), "cell flux is the weight-1 category flux")
   end subroutine test_stefan_first

   subroutine test_stefan_law(error)
      !! 800 steps of 3600 s from `h0 = 0.3 m`: `h^2 - h0^2` against `2 k |Ts| t/(rhoi Lfresh)` and the thickness of the
      !! quasi-steady profile.  The measured relative error of `h^2 - h0^2` (the four-layer scheme with the sensible heat
      !! of the new ice) is bracketed by `-0.045 < rel < -0.02` (the sensible heat of the new ice, Stefan number
      !! cp|Ts|/L = 0.126, slows the growth: about 4 % below the idealised law; a scheme without heat capacity gives 0,
      !! so a symmetric band around zero is not enough); measured -0.0328; the energy ledger of the run closes in the
      !! sense of `test_ledger`.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp), parameter :: H0 = 0.3_wp, TS = -20.0_wp, DT = 3600.0_wp
      integer, parameter :: NSTEP = 800
      integer :: s
      real(wp) :: h, hlaw, rel

      call stefan_world(w, H0, TS, DT)
      call world_map(w)
      do s = 1, NSTEP
         w%aicen_init = w%aicen
         call world_zero_fluxes(w)
         call world_step(w)
         if (w%ice_err(1, 1) /= 0_int32) then
            call world_unmap(w)
            call test_failed(error, "error bit during the Stefan run")
            return
         end if
      end do
      call world_unmap(w)
      h = w%vicen(1, 1, 1)
      hlaw = H0**2 + 2.0_wp*w%p%c%kice*abs(TS)*real(NSTEP, wp)*DT/(w%p%c%rhoi*w%p%c%Lfresh)
      rel = (h*h - H0*H0)/(hlaw - H0*H0) - 1.0_wp
      write (output_unit, '(a,f8.4,a,f8.4,a,es10.3)') "Stefan: h(end) = ", h, " m (law ", sqrt(hlaw), "), relative error of h^2 - h0^2 = ", rel
      call check(error, ieee_is_finite(rel), "finite")
      if (allocated(error)) return
      call check(error, rel > -0.045_wp .and. rel < -0.02_wp, "thickness follows the Stefan law on the side the sensible heat predicts (Stefan number 0.126: about 4 % below)")
      if (allocated(error)) return
      call check(error, abs(w%Tsfcn(1, 1, 1) - TS) < 1.0e-3_wp, "the surface is held at Ts")
   end subroutine test_stefan_law

   ! ------------------------------------------------------------------
   ! the ledger
   ! ------------------------------------------------------------------

   subroutine ledger_world(w)
      !! Five categories of unequal thickness, snow on the thick ones, a cold start.
      type(world_t), intent(inout) :: w
      real(wp), parameter :: A(5) = [0.05_wp, 0.15_wp, 0.30_wp, 0.25_wp, 0.10_wp]
      real(wp), parameter :: H(5) = [0.12_wp, 0.45_wp, 1.0_wp, 1.9_wp, 3.1_wp]
      real(wp), parameter :: HS(5) = [0.0_wp, 0.03_wp, 0.08_wp, 0.15_wp, 0.2_wp]
      integer :: n

      call world_alloc(w, 1, 1, 5, 4, ICE_CONSTANTS_ICEPACK)
      w%dt = 3600.0_wp
      w%sw_redist = .true.
      do n = 1, 5
         call set_column(w, 1, 1, n, A(n), H(n), HS(n), -9.0_wp, -7.0_wp, -2.0_wp, .false.)
      end do
   end subroutine ledger_world

   subroutine ledger_totals(w, M, S, E)
      type(world_t), intent(in) :: w
      real(wp), intent(out) :: M, S, E
      integer :: n, k
      real(wp) :: vi, vs, ei

      vi = 0.0_wp
      vs = 0.0_wp
      E = 0.0_wp
      do n = 1, w%ncat
         vi = vi + w%vicen(1, 1, n)
         vs = vs + w%vsnon(1, 1, n)
         ei = 0.0_wp
         do k = 1, w%nilyr
            ei = ei + w%qicen(1, 1, n, k)*w%vicen(1, 1, n)/real(w%nilyr, wp)
         end do
         E = E + ei + w%qsnon(1, 1, n, 1)*w%vsnon(1, 1, n)
      end do
      M = w%p%c%rhoi*vi + w%p%c%rhos*vs
      S = w%p%c%rhoi*vi*w%p%thermo%ice_ref_salinity*0.001_wp
   end subroutine ledger_totals

   subroutine ledger_force(w, s, nstep, kind, SW)
      !! Seasonal forcing of step `s` of `nstep` (kind 1: a full cycle, winter - summer - winter; kind 2: a warm sunny melt
      !! season).  Sets the air, the snowfall, the ocean and the shortwave inputs of every category.
      type(world_t), intent(inout) :: w
      integer, intent(in) :: s, nstep, kind
      real(wp), intent(out) :: SW
      real(wp) :: ph, Tair, qs
      integer :: n, k
      real(wp) :: wk(w%nilyr), wsum

      ph = TWOPI*real(s - 1, wp)/real(nstep, wp)
      if (kind == 1) then
         Tair = 266.0_wp - 18.0_wp*cos(ph)
         SW = 260.0_wp*(0.5_wp - 0.5_wp*cos(ph))**2
      else
         Tair = 277.0_wp + 2.0_wp*cos(ph)
         SW = 260.0_wp
      end if
      w%potT = Tair
      w%rhoa = 1.3_wp
      qs = 11637800.0_wp*exp(-5897.8_wp/Tair)/w%rhoa(1, 1)
      w%Qa = 0.8_wp*qs
      w%flw = 0.8_wp*5.67e-8_wp*Tair**4
      w%wind = 6.0_wp
      w%uatm = 6.0_wp
      w%vatm = 0.0_wp
      w%fsnow = merge(1.5e-4_wp, 0.0_wp, Tair < 270.0_wp .and. mod(s, 50) < 5)
      w%frzmlt = merge(-30.0_wp*max(0.0_wp, -cos(ph)), 0.0_wp, kind == 1)
      if (kind == 2) w%frzmlt = -20.0_wp
      w%sst = w%Tf + merge(0.2_wp, 0.0_wp, SW > 100.0_wp)
      w%tau_x = 0.03_wp
      wsum = 0.0_wp
      do k = 1, w%nilyr
         wk(k) = 0.5_wp**k
         wsum = wsum + wk(k)
      end do
      do n = 1, w%ncat
         w%fswsfcn(1, 1, n) = 0.25_wp*SW
         w%fswthrun(1, 1, n) = 0.03_wp*SW*real(n, wp)/5.0_wp
         do k = 1, w%nilyr
            w%Iswabsn(1, 1, n, k) = 0.08_wp*SW*wk(k)/wsum
         end do
         w%Sswabsn(1, 1, n, 1) = 0.0_wp
         w%fswintn(1, 1, n) = 0.08_wp*SW
      end do
   end subroutine ledger_force

   subroutine run_ledger(error, kind, nstep, w, wres, sres, hres, hratio, nmeltt, nsnow)
      !! `nstep` steps of the ledger world; returns the worst relative residual of the water, salt and heat lines
      !! (residual over the sum of the absolute values of every term of the line) and the number of steps with surface melt.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kind, nstep
      type(world_t), intent(inout) :: w
      real(wp), intent(out) :: wres, sres, hres, hratio
      integer, intent(out) :: nmeltt, nsnow
      integer :: s, n
      real(wp) :: SW, M0, S0, E0, M1, S1, E1, dt, a0, fsw_thru, lw, th, sc, fsn
      real(wp) :: res_w, res_s, res_h, term(8)

      wres = 0.0_wp
      sres = 0.0_wp
      hres = 0.0_wp
      hratio = 0.0_wp
      nmeltt = 0
      nsnow = 0
      dt = w%dt
      call world_map(w)
      do s = 1, nstep
         call ledger_force(w, s, nstep, kind, SW)
         w%aicen_init = w%aicen
         call world_zero_fluxes(w)
         call ledger_totals(w, M0, S0, E0)
         call world_step(w)
         if (w%ice_err(1, 1) /= 0_int32) then
            call world_unmap(w)
            call test_failed(error, "error bit raised in the ledger run")
            return
         end if
         call ledger_totals(w, M1, S1, E1)
         a0 = 0.0_wp
         fsw_thru = 0.0_wp
         do n = 1, w%ncat
            if (w%aicen_init(1, 1, n) > w%p%c%puny) then
               a0 = a0 + w%aicen_init(1, 1, n)
               fsw_thru = fsw_thru + w%aicen_init(1, 1, n)*w%fswthrun(1, 1, n)
            end if
         end do
         fsn = w%fsnow(1, 1)
         ! water
         term = [dt*fsn*a0, dt*w%evap(1, 1), -dt*w%fresh(1, 1), M1, M0, 0.0_wp, 0.0_wp, 0.0_wp]
         res_w = (M1 - M0) - (dt*fsn*a0 + dt*(w%evap(1, 1) - w%fresh(1, 1)))
         sc = sum(abs(term))
         wres = max(wres, abs(res_w)/sc)
         if (.not. ieee_is_finite(res_w)) wres = huge(1.0_wp)
         ! salt
         res_s = (S1 - S0) + dt*w%fsalt(1, 1)
         sc = abs(S1) + abs(S0) + abs(dt*w%fsalt(1, 1))
         sres = max(sres, abs(res_s)/max(sc, tiny(1.0_wp)))
         if (.not. ieee_is_finite(res_s)) sres = huge(1.0_wp)
         ! heat
         lw = dt*(w%fswabs(1, 1) - fsw_thru + w%flwout(1, 1) + w%fsens(1, 1) - w%fhocn(1, 1) + w%p%c%Lvap*w%evap(1, 1))
         th = dt*(a0*w%flw(1, 1) - w%p%c%Lfresh*fsn*a0)
         res_h = (E1 - E0) - (lw + th)
         sc = abs(E1) + abs(E0) + dt*(abs(w%fswabs(1, 1)) + abs(fsw_thru) + abs(w%flwout(1, 1)) + abs(w%fsens(1, 1)) &
                                      + abs(w%fhocn(1, 1)) + a0*w%flw(1, 1) + w%p%c%Lfresh*fsn*a0 + w%p%c%Lvap*abs(w%evap(1, 1)))
         hres = max(hres, abs(res_h)/sc)
         ! the BL99 solver accepts a pass at ferr <= 0.9*ferrmax W/m2 per ice area and the check at 1.1*ferrmax, so
         ! the cell residual over a step is bounded by 1.1*ferrmax*dt*sum(aicen_init) (plus round-off of the terms)
         hratio = max(hratio, abs(res_h)/(1.1_wp*ICE_FERRMAX*dt*a0 + 1.0e-15_wp*sc))
         if (.not. ieee_is_finite(res_h)) then
            hres = huge(1.0_wp)
            hratio = huge(1.0_wp)
         end if
         if (w%meltt(1, 1) > 0.0_wp) nmeltt = nmeltt + 1
         if (sum(w%vsnon(1, 1, :)) > 0.0_wp) nsnow = nsnow + 1
      end do
      call world_unmap(w)
   end subroutine run_ledger

   subroutine test_ledger(error)
      !! 1000 steps of a seasonal cycle (cold, then a sunny melt season with ocean heat, then cold): water, salt and heat
      !! close every step.  Water and salt to round-off; heat to the BL99 check bound (reported relative to the throughput).
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp) :: wres, sres, hres, hratio
      integer :: nm, ns

      call ledger_world(w)
      call run_ledger(error, 1, 1000, w, wres, sres, hres, hratio, nm, ns)
      if (allocated(error)) return
      write (output_unit, '(a,3es11.3,a,es10.3,a,i0,a,i0)') "ledger 1000 steps, worst relative residual water/salt/heat: ", wres, &
         sres, hres, "; worst heat residual / BL99 bound ", hratio, "; steps with top melt ", nm, ", steps with snow ", ns
      call check(error, nm > 50, "the season exercises top melt")
      if (allocated(error)) return
      call check(error, ns > 50, "the season exercises snow")
      if (allocated(error)) return
      call check(error, wres < 1.0e-14_wp, "water closes every step")
      if (allocated(error)) return
      call check(error, sres < 1.0e-14_wp, "salt closes every step")
      if (allocated(error)) return
      call check(error, hratio <= 1.0_wp, "heat closes every step within the BL99 solver bound")
   end subroutine test_ledger

   subroutine test_melt(error)
      !! A warm sunny spell (4 degC air, 260 W/m2 shortwave, ocean heat, `sw_redist` on): the snow melts first and the ice
      !! melts from the top and the base, the surface sits at the melting point, and the ledger closes.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp) :: wres, sres, hres, hratio, snow0, ice0
      integer :: nm, ns

      call ledger_world(w)
      call set_column(w, 1, 1, 3, 0.30_wp, 1.0_wp, 0.15_wp, -0.2_wp, -0.3_wp, -1.0_wp, .false.)
      w%Tsfcn(1, 1, :) = -0.2_wp
      snow0 = sum(w%vsnon)
      ice0 = sum(w%vicen)
      call run_ledger(error, 2, 300, w, wres, sres, hres, hratio, nm, ns)
      if (allocated(error)) return
      write (output_unit, '(a,3es11.3,a,i0,a,f6.3,a,f6.3)') "melt season 300 steps, worst relative residual water/salt/heat: ", &
         wres, sres, hres, "; melting steps ", nm, "; snow volume ", sum(w%vsnon), " m (from ", snow0
      call check(error, nm > 100, "the surface melts")
      if (allocated(error)) return
      call check(error, sum(w%vsnon) < 0.2_wp*snow0, "the snow is mostly gone")
      if (allocated(error)) return
      call check(error, sum(w%vicen) < ice0, "the ice thins")
      if (allocated(error)) return
      call check(error, wres < 1.0e-14_wp .and. sres < 1.0e-14_wp .and. hratio <= 1.0_wp, "the ledger closes through the melt")
   end subroutine test_melt

   ! ------------------------------------------------------------------
   ! categories and cells that must not interact
   ! ------------------------------------------------------------------

   subroutine fill_grid_world(w, nx, ny, ncat, nilyr, seed)
      !! A small multi-cell, multi-category world where every cell and category differs.
      type(world_t), intent(inout) :: w
      integer, intent(in) :: nx, ny, ncat, nilyr
      integer(int64), intent(in) :: seed
      integer :: i, j
      integer(int64) :: st

      call world_alloc(w, nx, ny, ncat, nilyr, ICE_CONSTANTS_ICEPACK)
      w%sw_redist = .true.
      st = seed
      do j = 1, ny
         do i = 1, nx
            call random_cell(w, i, j, st)
         end do
      end do
      w%aicen_init = w%aicen
   end subroutine fill_grid_world

   subroutine random_cell(w, i, j, st)
      !! Forcing and `ncat` categories of cell `(i, j)` from the stream `st`.
      type(world_t), intent(inout) :: w
      integer, intent(in) :: i, j
      integer(int64), intent(inout) :: st
      integer :: n, k
      real(wp) :: u, rsum, r(w%ncat), hin, hsn, Tsfc, Sk, Tm, T, Tair, SW, wk, ws
      logical :: summer

      summer = uniform(st) < 0.4_wp
      Tair = merge(271.5_wp + 13.0_wp*uniform(st), 245.0_wp + 30.0_wp*uniform(st), summer)
      w%potT(i, j) = Tair + (2.0_wp*uniform(st) - 1.0_wp)
      w%Qa(i, j) = 5.0e-4_wp + 3.5e-3_wp*uniform(st)
      w%rhoa(i, j) = 1.15_wp + 0.2_wp*uniform(st)
      u = 0.5_wp + 14.0_wp*uniform(st)
      w%uatm(i, j) = u*(2.0_wp*uniform(st) - 1.0_wp)
      w%vatm(i, j) = u*(2.0_wp*uniform(st) - 1.0_wp)
      w%wind(i, j) = sqrt(w%uatm(i, j)**2 + w%vatm(i, j)**2)
      w%flw(i, j) = 150.0_wp + 170.0_wp*uniform(st)
      w%fsnow(i, j) = merge(3.0e-4_wp*uniform(st), 0.0_wp, uniform(st) > 0.35_wp)
      w%uvel(i, j) = 0.4_wp*(2.0_wp*uniform(st) - 1.0_wp)
      w%vvel(i, j) = 0.4_wp*(2.0_wp*uniform(st) - 1.0_wp)
      w%Tf(i, j) = -1.9_wp + 0.3_wp*uniform(st)
      w%sst(i, j) = w%Tf(i, j) + 0.4_wp*uniform(st)
      w%frzmlt(i, j) = merge(-150.0_wp*uniform(st), 30.0_wp*uniform(st), uniform(st) < 0.5_wp)
      w%tau_x(i, j) = 0.3_wp*(2.0_wp*uniform(st) - 1.0_wp)
      w%tau_y(i, j) = 0.3_wp*(2.0_wp*uniform(st) - 1.0_wp)
      rsum = 0.0_wp
      do n = 1, w%ncat
         u = uniform(st)
         if (u < 0.2_wp) then
            r(n) = 0.0_wp
         else
            r(n) = 0.2_wp + 0.8_wp*uniform(st)
            rsum = rsum + r(n)
         end if
      end do
      if (rsum == 0.0_wp) then
         r(1) = 1.0_wp
         rsum = 1.0_wp
      end if
      do n = 1, w%ncat
         w%aicen(i, j, n) = 0.92_wp*r(n)/rsum
         u = uniform(st)
         hin = 0.03_wp + 3.5_wp*u*u
         hsn = merge(0.45_wp*uniform(st), 0.0_wp, uniform(st) > 0.35_wp)
         w%vicen(i, j, n) = w%aicen(i, j, n)*hin
         w%vsnon(i, j, n) = w%aicen(i, j, n)*hsn
         Tsfc = -0.05_wp - 34.0_wp*uniform(st)
         w%Tsfcn(i, j, n) = Tsfc
         do k = 1, w%nilyr
            Sk = 0.3_wp + 5.5_wp*uniform(st)
            w%sicen(i, j, n, k) = Sk
            Tm = -0.054_wp*Sk
            T = Tm - 1.0e-3_wp + (min(Tsfc, Tm) - Tm)*0.9_wp*uniform(st)
            w%qicen(i, j, n, k) = ice_q_from_T(T, Sk, w%p%c, w%p%thermo%l_brine)
         end do
         T = min(Tsfc*(0.2_wp + 0.8_wp*uniform(st)), -0.01_wp)
         w%qsnon(i, j, n, 1) = ice_qsno_from_T(T, w%p%c)
         SW = merge(120.0_wp + 140.0_wp*uniform(st), 90.0_wp*uniform(st), summer)
         w%fswsfcn(i, j, n) = SW
         w%fswthrun(i, j, n) = merge(25.0_wp, 15.0_wp, summer)*uniform(st)
         ws = 0.0_wp
         do k = 1, w%nilyr
            wk = 10.0_wp*uniform(st)*0.5_wp**k
            w%Iswabsn(i, j, n, k) = wk
            ws = ws + wk
         end do
         w%Sswabsn(i, j, n, 1) = merge(5.0_wp*uniform(st), 0.0_wp, hsn > 1.0e-3_wp)
         w%fswintn(i, j, n) = ws + w%Sswabsn(i, j, n, 1)
      end do
   end subroutine random_cell

   subroutine test_empty(error)
      !! A category with `aicen_init <= puny` (exactly 0, or half of `puny`) holds garbage state: it must come out BITWISE
      !! unchanged, with zero per-category outputs and a clean error word, and add nothing to the cell.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, ref
      integer :: n, k

      call fill_grid_world(w, 3, 2, 3, 4, 777_int64)
      ! categories 2 and 3 of every cell carry garbage state and (almost) no area at the start of the step
      do n = 2, 3
         w%aicen_init(:, :, n) = merge(0.0_wp, 0.5_wp*w%p%c%puny, n == 2)
         w%aicen(:, :, n) = w%aicen_init(:, :, n)
         w%vicen(:, :, n) = 0.7_wp
         w%vsnon(:, :, n) = 0.3_wp
         w%Tsfcn(:, :, n) = -3.0_wp
         w%qicen(:, :, n, :) = -2.5e8_wp
         w%qsnon(:, :, n, 1) = -1.0e8_wp
      end do
      ref = w
      call world_map(w)
      call world_step(w)
      call world_unmap(w)
      do n = 2, 3
         call check(error, all(transfer(w%aicen(:, :, n), [0_int64]) == transfer(ref%aicen(:, :, n), [0_int64])), "aicen untouched")
         if (allocated(error)) return
         call check(error, all(transfer(w%vicen(:, :, n), [0_int64]) == transfer(ref%vicen(:, :, n), [0_int64])), "vicen untouched")
         if (allocated(error)) return
         call check(error, all(transfer(w%vsnon(:, :, n), [0_int64]) == transfer(ref%vsnon(:, :, n), [0_int64])), "vsnon untouched")
         if (allocated(error)) return
         call check(error, all(transfer(w%Tsfcn(:, :, n), [0_int64]) == transfer(ref%Tsfcn(:, :, n), [0_int64])), "Tsfcn untouched")
         if (allocated(error)) return
         do k = 1, 4
            call check(error, all(transfer(w%qicen(:, :, n, k), [0_int64]) == transfer(ref%qicen(:, :, n, k), [0_int64])), &
                       "qicen untouched")
            if (allocated(error)) return
         end do
         call check(error, all(transfer(w%qsnon(:, :, n, 1), [0_int64]) == transfer(ref%qsnon(:, :, n, 1), [0_int64])), &
                    "qsnon untouched")
         if (allocated(error)) return
         call check(error, all(w%ierr_n(:, :, n) == 0_int32), "clean per-category error word")
         if (allocated(error)) return
         call check(error, all(w%freshn(:, :, n) == 0.0_wp) .and. all(w%fhocnn(:, :, n) == 0.0_wp) .and. &
                    all(w%strairxn(:, :, n) == 0.0_wp) .and. all(w%melttn(:, :, n) == 0.0_wp), "outputs are zero")
         if (allocated(error)) return
      end do
      call check(error, all(w%ice_err == 0_int32), "no cell error")
      if (allocated(error)) return
      ! the active category (1) is the same as in a world without the other two categories
      call check(error, all(w%vicen(:, :, 1) /= ref%vicen(:, :, 1) .or. ref%aicen(:, :, 1) == 0.0_wp), "category 1 stepped")
   end subroutine test_empty

   subroutine test_errors(error)
      !! NaN in the enthalpy of ONE category of ONE cell flags that cell only; the other categories of the cell, the other
      !! cells, and the flagged category's state are as in a clean run (the state of the flagged category is untouched).  A
      !! non-finite forcing flags its own cell; a NaN `aicen_init` is flagged, never read as "no ice".  The random world is
      !! chosen (first seed) so that the CLEAN run raises no error of its own.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, clean, base
      integer :: i, j, n, ib, jb, nb, k, iseed
      integer(int32) :: bits
      logical :: same, ok

      ib = 2
      jb = 2
      ok = .false.
      do iseed = 1, 60
         call fill_grid_world(base, 4, 3, 3, 4, 9000_int64 + iseed)
         w = base
         call world_map(w)
         call world_step(w)
         call world_unmap(w)
         if (all(w%ice_err == 0_int32) .and. maxval(base%aicen_init(ib, jb, :)) > 0.1_wp) then
            ok = .true.
            exit
         end if
      end do
      call check(error, ok, "found a random world whose clean run raises no error")
      if (allocated(error)) return
      nb = maxloc(base%aicen_init(ib, jb, :), dim=1)
      ! the basal heat of a cell reads the enthalpy of all its categories while `frzmlt < 0`: with a melt potential the NaN
      ! would refuse the whole cell (fbot = 0), so this cell has none
      base%frzmlt(ib, jb) = 10.0_wp
      w = base
      clean = base
      call world_map(clean)
      call world_step(clean)
      call world_unmap(clean)
      call check(error, all(clean%ice_err == 0_int32), "the clean run is error-free")
      if (allocated(error)) return
      w%qicen(ib, jb, nb, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
      call world_map(w)
      call world_step(w)
      call world_unmap(w)
      do j = 1, 3
         do i = 1, 4
            bits = w%ice_err(i, j)
            if (i == ib .and. j == jb) then
               call check(error, iand(bits, ICE_ERR_NONFINITE) /= 0_int32, "the flagged cell carries the bit")
               if (allocated(error)) return
               call check(error, w%ierr_n(i, j, nb) /= 0_int32, "the flagged category carries the bit")
               if (allocated(error)) return
               call check(error, w%vicen(i, j, nb) == base%vicen(i, j, nb) .and. w%aicen(i, j, nb) == base%aicen(i, j, nb), &
                          "flagged category not committed")
               if (allocated(error)) return
               call check(error, w%freshn(i, j, nb) == 0.0_wp .and. w%fhocnn(i, j, nb) == 0.0_wp, "flagged category outputs zero")
               if (allocated(error)) return
               do n = 1, 3
                  if (n == nb) cycle
                  call check(error, w%ierr_n(i, j, n) == 0_int32, "other categories of the cell are clean")
                  if (allocated(error)) return
                  same = w%vicen(i, j, n) == clean%vicen(i, j, n) .and. w%Tsfcn(i, j, n) == clean%Tsfcn(i, j, n) .and. &
                         w%freshn(i, j, n) == clean%freshn(i, j, n) .and. w%qicen(i, j, n, 2) == clean%qicen(i, j, n, 2)
                  call check(error, same, "other categories of the flagged cell are as in the clean run")
                  if (allocated(error)) return
               end do
            else
               call check(error, bits == 0_int32, "no other cell is flagged")
               if (allocated(error)) return
               same = .true.
               do n = 1, 3
                  same = same .and. w%vicen(i, j, n) == clean%vicen(i, j, n) .and. w%freshn(i, j, n) == clean%freshn(i, j, n) &
                         .and. w%fsaltn(i, j, n) == clean%fsaltn(i, j, n)
                  do k = 1, 4
                     same = same .and. w%qicen(i, j, n, k) == clean%qicen(i, j, n, k)
                  end do
               end do
               same = same .and. w%fresh(i, j) == clean%fresh(i, j) .and. w%strair_x_T(i, j) == clean%strair_x_T(i, j)
               call check(error, same, "other cells are bitwise as in the clean run")
               if (allocated(error)) return
            end if
         end do
      end do

      ! a non-finite forcing and a NaN start-of-step area flag their own cells only
      w = base
      w%potT(3, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      w%aicen_init(1, 3, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call world_map(w)
      call world_step(w)
      call world_unmap(w)
      do j = 1, 3
         do i = 1, 4
            if (i == 3 .and. j == 1) then
               call check(error, iand(w%ice_err(i, j), ICE_ERR_NONFINITE) /= 0_int32, "non-finite forcing is flagged")
            else if (i == 1 .and. j == 3) then
               call check(error, iand(w%ice_err(i, j), ICE_ERR_NONFINITE) /= 0_int32, "NaN area is flagged")
            else
               call check(error, w%ice_err(i, j) == 0_int32, "unaffected cells are clean")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_errors

   subroutine test_flwout_gate(error)
      !! The reflected incoming longwave `(1 - emis)*flw` of the `flwout` merge belongs to the categories the thermodynamics
      !! ran (upstream merges only those, `aicen_init > puny`): with the merge's `puny` raised to 0.3, the categories of area
      !! 0.2 and 0.1 add nothing to `flwout` (a skipped category's own `flwoutn` is zero), the one of area 0.5 adds
      !! `a*(flwoutn - (1 - emis)*flw)`.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp), parameter :: A(3) = [0.5_wp, 0.2_wp, 0.1_wp]
      real(wp) :: want, wrong, emis
      integer :: n, iseed
      logical :: ok

      ok = .false.
      do iseed = 1, 400
         call fill_grid_world(w, 1, 1, 3, 4, 31000_int64 + iseed)
         do n = 1, 3
            ! the random columns keep their thickness and snow depth; only the areas are set
            if (w%aicen(1, 1, n) > 0.0_wp) then
               w%vicen(1, 1, n) = w%vicen(1, 1, n)/w%aicen(1, 1, n)*A(n)
               w%vsnon(1, 1, n) = w%vsnon(1, 1, n)/w%aicen(1, 1, n)*A(n)
            else
               w%vicen(1, 1, n) = 0.5_wp*A(n)
               w%vsnon(1, 1, n) = 0.0_wp
            end if
            w%aicen(1, 1, n) = A(n)
         end do
         w%aicen_init = w%aicen
         w%manual = .true.
         w%frzmlt = 0.0_wp
         w%merge_puny = 0.3_wp
         call world_map(w)
         call world_step(w)
         call world_unmap(w)
         if (w%ice_err(1, 1) == 0_int32) then
            ok = .true.
            exit
         end if
      end do
      call check(error, ok, "a clean random column was found")
      if (allocated(error)) return
      emis = w%p%thermo%emissivity
      want = 0.0_wp
      wrong = 0.0_wp
      do n = 1, 3
         if (A(n) > 0.3_wp) want = want + A(n)*(w%flwoutn(1, 1, n) - (1.0_wp - emis)*w%flw(1, 1))
         wrong = wrong + A(n)*(w%flwoutn(1, 1, n) - (1.0_wp - emis)*w%flw(1, 1))
      end do
      call check(error, near(w%flwout(1, 1), want, SUMTOL + 1.0e-15_wp, 0.0_wp), "reflected longwave only for active categories")
      if (allocated(error)) return
      call check(error, abs(wrong - want) > 1.0_wp, "the gate is visible (the ungated sum differs)")
   end subroutine test_flwout_gate

   subroutine test_merge(error)
      !! The cell fluxes are the `aicen_init`-weighted sums of the category values (hand-summed in category order) ADDED to
      !! what the cell held; the flwout merge subtracts the reflected longwave; the stress is SET (not added to what was
      !! there); only the owned range is written.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      integer :: i, j, n
      real(wp) :: s(14), pre(12), a, emis, x(14)

      call fill_grid_world(w, 5, 4, 5, 7, 99_int64)
      w%is = 2
      w%ie = 4
      w%js = 2
      w%je = 3
      w%fresh = 0.37_wp; w%fsalt = 0.011_wp; w%fhocn = 3.5_wp; w%evap = 1.0e-5_wp; w%fsens = -7.0_wp; w%flat = 2.5_wp
      w%flwout = -300.0_wp; w%fswabs = 12.0_wp; w%meltt = 0.002_wp; w%meltb = 0.003_wp; w%congel = 0.004_wp; w%snoice = 0.005_wp
      w%strair_x_T = 9.0_wp
      w%strair_y_T = 9.0_wp
      pre = [0.37_wp, 0.011_wp, 3.5_wp, 1.0e-5_wp, -7.0_wp, 2.5_wp, -300.0_wp, 12.0_wp, 0.002_wp, 0.003_wp, 0.004_wp, 0.005_wp]
      emis = w%p%thermo%emissivity
      call world_map(w)
      call world_step(w)
      call world_unmap(w)
      do j = 1, 4
         do i = 1, 5
            if (i < 2 .or. i > 4 .or. j < 2 .or. j > 3) then
               call check(error, w%fresh(i, j) == 0.37_wp .and. w%strair_x_T(i, j) == 9.0_wp .and. w%flwout(i, j) == -300.0_wp, &
                          "outside the owned range nothing is written")
               if (allocated(error)) return
               cycle
            end if
            s = 0.0_wp
            do n = 1, 5
               a = w%aicen_init(i, j, n)
               s(1) = s(1) + a*w%freshn(i, j, n)
               s(2) = s(2) + a*w%fsaltn(i, j, n)
               s(3) = s(3) + a*w%fhocnn(i, j, n)
               s(4) = s(4) + a*w%evapn(i, j, n)
               s(5) = s(5) + a*w%fsensn(i, j, n)
               s(6) = s(6) + a*w%flatn(i, j, n)
               s(7) = s(7) + a*(w%flwoutn(i, j, n) - (1.0_wp - emis)*w%flw(i, j))
               s(8) = s(8) + a*w%fswabsn(i, j, n)
               s(9) = s(9) + a*w%melttn(i, j, n)
               s(10) = s(10) + a*w%meltbn(i, j, n)
               s(11) = s(11) + a*w%congeln(i, j, n)
               s(12) = s(12) + a*w%snoicen(i, j, n)
               s(13) = s(13) + a*w%strairxn(i, j, n)
               s(14) = s(14) + a*w%strairyn(i, j, n)
            end do
            x = [w%fresh(i, j), w%fsalt(i, j), w%fhocn(i, j), w%evap(i, j), w%fsens(i, j), w%flat(i, j), w%flwout(i, j), &
                 w%fswabs(i, j), w%meltt(i, j), w%meltb(i, j), w%congel(i, j), w%snoice(i, j), w%strair_x_T(i, j), &
                 w%strair_y_T(i, j)]
            do n = 1, 12
               call check(error, near(x(n), pre(n) + s(n), SUMTOL, 0.0_wp), "cell flux = pre-existing + weighted category sum")
               if (allocated(error)) return
            end do
          call check(error, near(x(13), s(13), SUMTOL, 0.0_wp) .and. near(x(14), s(14), SUMTOL, 0.0_wp), "stress is set to the sum")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_merge

   ! ------------------------------------------------------------------
   ! goldens
   ! ------------------------------------------------------------------

   subroutine load_goldens(set, GA, GB, GC)
      integer, intent(in) :: set
      real(wp), intent(out) :: GA(273, 12), GB(273, 6), GC(213, 6)
#include "test_cice_therm1_golden.inc"
   end subroutine load_goldens

   subroutine load_row(w, v, nil)
      !! Fill the single-cell world from an input vector of the oracle table (layout in the oracle header).
      type(world_t), intent(inout) :: w
      real(wp), intent(in) :: v(:)
      integer, intent(in) :: nil
      integer :: n, k, p0

      w%dt = v(1)
      w%potT(1, 1) = v(2)
      w%uatm(1, 1) = v(3)
      w%vatm(1, 1) = v(4)
      w%wind(1, 1) = v(5)
      w%Qa(1, 1) = v(6)
      w%rhoa(1, 1) = v(7)
      w%flw(1, 1) = v(8)
      w%fsnow(1, 1) = v(9)
      w%uvel(1, 1) = v(10)
      w%vvel(1, 1) = v(11)
      w%Tf(1, 1) = v(12)
      w%sst(1, 1) = v(13)
      w%frzmlt(1, 1) = v(14)
      w%tau_x(1, 1) = v(15)
      w%tau_y(1, 1) = v(16)
      p0 = 16
      do n = 1, 5
         w%aicen(1, 1, n) = v(p0 + 1)
         w%aicen_init(1, 1, n) = v(p0 + 1)
         w%vicen(1, 1, n) = v(p0 + 2)
         w%vsnon(1, 1, n) = v(p0 + 3)
         w%Tsfcn(1, 1, n) = v(p0 + 4)
         p0 = p0 + 4
         do k = 1, nil
            w%qicen(1, 1, n, k) = v(p0 + k)
            w%sicen(1, 1, n, k) = v(p0 + nil + k)
         end do
         p0 = p0 + 2*nil
         w%qsnon(1, 1, n, 1) = v(p0 + 1)
         w%fswsfcn(1, 1, n) = v(p0 + 2)
         w%fswintn(1, 1, n) = v(p0 + 3)
         w%fswthrun(1, 1, n) = v(p0 + 4)
         p0 = p0 + 4
         do k = 1, nil
            w%Iswabsn(1, 1, n, k) = v(p0 + k)
         end do
         p0 = p0 + nil
         w%Sswabsn(1, 1, n, 1) = v(p0 + 1)
         p0 = p0 + 1
      end do
   end subroutine load_row

   subroutine collect_row(w, nil, got)
      !! The outputs in the order of the oracle table.
      type(world_t), intent(in) :: w
      integer, intent(in) :: nil
      real(wp), intent(out) :: got(:)
      integer :: n, k, q

      q = 0
      do n = 1, 5
         got(q + 1) = w%aicen(1, 1, n)
         got(q + 2) = w%vicen(1, 1, n)
         got(q + 3) = w%vsnon(1, 1, n)
         got(q + 4) = w%Tsfcn(1, 1, n)
         q = q + 4
         do k = 1, nil
            got(q + k) = w%qicen(1, 1, n, k)
         end do
         q = q + nil
         got(q + 1) = w%qsnon(1, 1, n, 1)
         got(q + 2) = w%melttn(1, 1, n)
         got(q + 3) = w%meltbn(1, 1, n)
         got(q + 4) = w%congeln(1, 1, n)
         got(q + 5) = w%snoicen(1, 1, n)
         got(q + 6) = w%fsensn(1, 1, n)
         got(q + 7) = w%flatn(1, 1, n)
         q = q + 7
      end do
      got(q + 1:q + 17) = [w%fresh(1, 1), w%fsalt(1, 1), w%fhocn(1, 1), w%evap(1, 1), w%fsens(1, 1), w%flat(1, 1), &
                           w%flwout(1, 1), w%fswabs(1, 1), w%meltt(1, 1), w%meltb(1, 1), w%congel(1, 1), w%snoice(1, 1), &
                           w%strair_x_T(1, 1), w%strair_y_T(1, 1), w%fbot(1, 1), w%rside(1, 1), w%wlat(1, 1)]
   end subroutine collect_row

   subroutine golden_group(error, G, ncell, nil, hf, swr, set, tag, worst)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: G(:, :)
      integer, intent(in) :: ncell, nil, set
      logical, intent(in) :: hf, swr
      character(len=*), intent(in) :: tag
      real(wp), intent(inout) :: worst(4)
      type(world_t) :: w
      integer :: r, q, nin, nout, n, cls
      real(wp), allocatable :: got(:)
      real(wp) :: want, rel
      character(len=200) :: msg

      nin = 16 + 5*(9 + 3*nil)
      nout = 5*(11 + nil) + 17
      allocate (got(nout))
      do r = 1, ncell
         call world_alloc(w, 1, 1, 5, nil, set)
         w%p%atm%highfreq = hf
         w%sw_redist = swr
         call load_row(w, G(1:nin, r), nil)
         call world_map(w)
         call world_step(w)
         call world_unmap(w)
         if (w%ice_err(1, 1) /= 0_int32) then
            write (msg, '(a,a,i0,a,i0)') tag, " row ", r, ": error word ", w%ice_err(1, 1)
            call test_failed(error, trim(msg))
            return
         end if
         call collect_row(w, nil, got)
         do q = 1, nout
            want = G(nin + q, r)
            if (.not. ieee_is_finite(got(q))) then
               write (msg, '(a,a,i0,a,i0)') tag, " row ", r, ": non-finite output ", q
               call test_failed(error, trim(msg))
               return
            end if
            ! output class: 1 category state, 2 category flux, 3 cell flux/stress, 4 basal heat
            if (q <= 5*(11 + nil)) then
               n = mod(q - 1, 11 + nil) + 1
               cls = merge(1, 2, n <= 4 + nil + 1)
            else
               cls = merge(3, 4, q - 5*(11 + nil) <= 14)
            end if
            ! scale of the output's group: a flux that cancels to ~0 is compared against its terms' size
            rel = abs(got(q) - want)/max(abs(want), 1.0e-30_wp)
            if (abs(want) < 1.0e-14_wp .and. abs(got(q)) < 1.0e-14_wp) rel = 0.0_wp
            if (rel > worst(cls)) worst(cls) = rel
            if (.not. near(got(q), want, GRTOL, golden_atol(q, nil))) then
               write (msg, '(a,a,i0,a,i0,a,es24.16,a,es24.16)') tag, " row ", r, " output ", q, ": got ", got(q), " want ", want
               call test_failed(error, trim(msg))
               return
            end if
         end do
      end do
   end subroutine golden_group

   function golden_atol(q, nil) result(a)
      !! Absolute floor of the golden comparison for output `q` (table order): an output that cancels to round-off (a
      !! `1e-19` melt residue, a heat flux that is zero up to its terms) is compared to the unit of its quantity times 1e-12,
      !! never to itself; everything else is purely relative.
      integer, intent(in) :: q, nil
      real(wp) :: a
      integer :: n, c0
      a = 0.0_wp
      if (q <= 5*(11 + nil)) then
         n = mod(q - 1, 11 + nil) + 1
         if (n <= 3) then
            a = 1.0e-15_wp                  ! area, volumes [m]
         else if (n == 4) then
            a = 1.0e-12_wp                  ! surface temperature [degC]
         else if (n <= 4 + nil + 1) then
            a = 1.0e-3_wp                   ! enthalpy [J/m3]
         else if (n <= 4 + nil + 5) then
            a = 1.0e-15_wp                  ! melt / growth accumulators [m]
         else
            a = 1.0e-9_wp                   ! sensible, latent [W/m2]
         end if
      else
         c0 = q - 5*(11 + nil)
         select case (c0)
         case (1, 2, 4)
            a = 1.0e-17_wp                  ! fresh, fsalt, evap [kg/m2/s]
         case (3, 5, 6, 7, 8, 15)
            a = 1.0e-9_wp                   ! heat fluxes [W/m2]
         case (9:12, 16, 17)
            a = 1.0e-15_wp                  ! accumulators [m], rside, wlat
         case default
            a = 1.0e-12_wp                  ! stress [N/m2]
         end select
      end if
   end function golden_atol

   subroutine test_golden(error)
      !! Both constant sets, three groups (nilyr 7 / highfreq off / sw_redist off; nilyr 7 / on / on; nilyr 4 / off / on).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: GA(273, 12), GB(273, 6), GC(213, 6)
      real(wp) :: worst(4)
      integer :: set

      worst = 0.0_wp
      do set = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         call load_goldens(set, GA, GB, GC)
         call golden_group(error, GA, 12, 7, .false., .false., set, "group A", worst)
         if (allocated(error)) return
         call golden_group(error, GB, 6, 7, .true., .true., set, "group B", worst)
         if (allocated(error)) return
         call golden_group(error, GC, 6, 4, .false., .true., set, "group C", worst)
         if (allocated(error)) return
      end do
    write (output_unit, '(a,4es10.3)') "golden worst relative difference [category state, category flux, cell flux, basal]: ", worst
   end subroutine test_golden

   ! ------------------------------------------------------------------
   ! the device sweep
   ! ------------------------------------------------------------------

   subroutine cmp_all(error, a, b, label)
      !! Compare every array of two worlds (every state, forcing-derived output and error word), NaN-proof.
      type(error_type), allocatable, intent(out) :: error
      type(world_t), intent(in) :: a, b
      character(len=*), intent(in) :: label
      call cmp2(error, label//" uatm", a%uatm, b%uatm, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" vatm", a%vatm, b%vatm, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" wind", a%wind, b%wind, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" potT", a%potT, b%potT, .true.)
      if (allocated(error)) return
      call cmp2(error, label//" Qa", a%Qa, b%Qa, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" rhoa", a%rhoa, b%rhoa, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" flw", a%flw, b%flw, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fsnow", a%fsnow, b%fsnow, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" uvel", a%uvel, b%uvel, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" vvel", a%vvel, b%vvel, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" Tf", a%Tf, b%Tf, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" sst", a%sst, b%sst, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" frzmlt", a%frzmlt, b%frzmlt, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" tau_x", a%tau_x, b%tau_x, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" tau_y", a%tau_y, b%tau_y, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fresh", a%fresh, b%fresh, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fsalt", a%fsalt, b%fsalt, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fhocn", a%fhocn, b%fhocn, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" evap", a%evap, b%evap, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fsens", a%fsens, b%fsens, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" flat", a%flat, b%flat, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" flwout", a%flwout, b%flwout, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fswabs", a%fswabs, b%fswabs, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" meltt", a%meltt, b%meltt, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" meltb", a%meltb, b%meltb, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" congel", a%congel, b%congel, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" snoice", a%snoice, b%snoice, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" strair_x_T", a%strair_x_T, b%strair_x_T, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" strair_y_T", a%strair_y_T, b%strair_y_T, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" fbot", a%fbot, b%fbot, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" rside", a%rside, b%rside, .false.)
      if (allocated(error)) return
      call cmp2(error, label//" wlat", a%wlat, b%wlat, .false.)
      if (allocated(error)) return
      call cmp_i2(error, label//" ice_err", a%ice_err, b%ice_err)
      if (allocated(error)) return
      call cmp3(error, label//" aicen_init", a%aicen_init, b%aicen_init, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" aicen", a%aicen, b%aicen, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" vicen", a%vicen, b%vicen, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" vsnon", a%vsnon, b%vsnon, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" Tsfcn", a%Tsfcn, b%Tsfcn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fswsfcn", a%fswsfcn, b%fswsfcn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fswintn", a%fswintn, b%fswintn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fswthrun", a%fswthrun, b%fswthrun, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" freshn", a%freshn, b%freshn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fsaltn", a%fsaltn, b%fsaltn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fhocnn", a%fhocnn, b%fhocnn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" evapn", a%evapn, b%evapn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fsensn", a%fsensn, b%fsensn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" flatn", a%flatn, b%flatn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" flwoutn", a%flwoutn, b%flwoutn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" fswabsn", a%fswabsn, b%fswabsn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" melttn", a%melttn, b%melttn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" meltbn", a%meltbn, b%meltbn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" congeln", a%congeln, b%congeln, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" snoicen", a%snoicen, b%snoicen, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" strairxn", a%strairxn, b%strairxn, .false.)
      if (allocated(error)) return
      call cmp3(error, label//" strairyn", a%strairyn, b%strairyn, .false.)
      if (allocated(error)) return
      call cmp_i3(error, label//" ierr_n", a%ierr_n, b%ierr_n)
      if (allocated(error)) return
      call cmp4(error, label//" qicen", a%qicen, b%qicen, .true.)
      if (allocated(error)) return
      call cmp4(error, label//" sicen", a%sicen, b%sicen, .false.)
      if (allocated(error)) return
      call cmp4(error, label//" Iswabsn", a%Iswabsn, b%Iswabsn, .false.)
      if (allocated(error)) return
      call cmp4(error, label//" qsnon", a%qsnon, b%qsnon, .false.)
      if (allocated(error)) return
      call cmp4(error, label//" Sswabsn", a%Sswabsn, b%Sswabsn, .false.)
      if (allocated(error)) return
   end subroutine cmp_all

   function ieee_finite_mask2(y) result(m)
      real(wp), intent(in) :: y(:, :)
      logical :: m(size(y, 1), size(y, 2))
      m = abs(y) < huge(1.0_wp)
   end function ieee_finite_mask2

   function ieee_finite_mask3(y) result(m)
      real(wp), intent(in) :: y(:, :, :)
      logical :: m(size(y, 1), size(y, 2), size(y, 3))
      m = abs(y) < huge(1.0_wp)
   end function ieee_finite_mask3

   subroutine cmp2(error, label, x, y, allow_nan)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: x(:, :), y(:, :)
      logical, intent(in) :: allow_nan
      integer :: i, j
      real(wp) :: atol
      character(len=200) :: msg
      atol = DRTOL*maxval(abs(y(:, :)), mask=ieee_finite_mask2(y))   ! scale-aware floor on the GPU build, 0 on CPU builds
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (allow_nan .and. ieee_is_nan(x(i, j)) .and. ieee_is_nan(y(i, j))) cycle
            if (.not. (ieee_is_finite(x(i, j)) .and. ieee_is_finite(y(i, j)))) then
               write (msg, '(a,a,2(1x,i0))') "non-finite ", label, i, j
               call test_failed(error, trim(msg))
               return
            end if
            if (.not. near(x(i, j), y(i, j), DRTOL, atol)) then
               write (msg, '(a,a,2(1x,i0),2es24.16)') "differs ", label, i, j, x(i, j), y(i, j)
               call test_failed(error, trim(msg))
               return
            end if
         end do
      end do
   end subroutine cmp2

   subroutine cmp3(error, label, x, y, allow_nan)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: x(:, :, :), y(:, :, :)
      logical, intent(in) :: allow_nan
      integer :: i, j, n
      real(wp) :: atol
      character(len=200) :: msg
      atol = DRTOL*maxval(abs(y(:, :, :)), mask=ieee_finite_mask3(y))
      do n = 1, size(x, 3)
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               if (allow_nan .and. ieee_is_nan(x(i, j, n)) .and. ieee_is_nan(y(i, j, n))) cycle
               if (.not. (ieee_is_finite(x(i, j, n)) .and. ieee_is_finite(y(i, j, n)))) then
                  write (msg, '(a,a,3(1x,i0))') "non-finite ", label, i, j, n
                  call test_failed(error, trim(msg))
                  return
               end if
               if (.not. near(x(i, j, n), y(i, j, n), DRTOL, atol)) then
                  write (msg, '(a,a,3(1x,i0),2es24.16)') "differs ", label, i, j, n, x(i, j, n), y(i, j, n)
                  call test_failed(error, trim(msg))
                  return
               end if
            end do
         end do
      end do
   end subroutine cmp3

   subroutine cmp4(error, label, x, y, allow_nan)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: x(:, :, :, :), y(:, :, :, :)
      logical, intent(in) :: allow_nan
      integer :: k
      do k = 1, size(x, 4)
         call cmp3(error, label, x(:, :, :, k), y(:, :, :, k), allow_nan)
         if (allocated(error)) return
      end do
   end subroutine cmp4

   subroutine cmp_i2(error, label, x, y)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer(int32), intent(in) :: x(:, :), y(:, :)
      character(len=200) :: msg
      integer :: i, j
      if (any((x /= 0_int32) .neqv. (y /= 0_int32)) .or. any(iand(x, ICE_ERR_NONFINITE) /= iand(y, ICE_ERR_NONFINITE))) then
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               if ((x(i, j) /= 0_int32) .neqv. (y(i, j) /= 0_int32)) then
                  write (output_unit, '(a,2i5,2i12)') "first differing error word (i, j, device, host): ", i, j, x(i, j), y(i, j)
                  exit
               end if
            end do
            if (any((x(:, j) /= 0_int32) .neqv. (y(:, j) /= 0_int32))) exit
         end do
         write (msg, '(a,a,i0)') "error words differ ", label, count(x /= y)
         call test_failed(error, trim(msg))
      end if
   end subroutine cmp_i2

   subroutine cmp_i3(error, label, x, y)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer(int32), intent(in) :: x(:, :, :), y(:, :, :)
      character(len=200) :: msg
      if (any((x /= 0_int32) .neqv. (y /= 0_int32)) .or. any(iand(x, ICE_ERR_NONFINITE) /= iand(y, ICE_ERR_NONFINITE))) then
         write (msg, '(a,a,i0)') "error words differ ", label, count(x /= y)
         call test_failed(error, trim(msg))
      end if
   end subroutine cmp_i3

   subroutine make_smooth(w, i0, j0)
      !! Every cell gets the forcing and the categories of cell `(i0, j0)`.
      type(world_t), intent(inout) :: w
      integer, intent(in) :: i0, j0
      integer :: n, k

      w%uatm = w%uatm(i0, j0)
      w%vatm = w%vatm(i0, j0)
      w%wind = w%wind(i0, j0)
      w%potT = w%potT(i0, j0)
      w%Qa = w%Qa(i0, j0)
      w%rhoa = w%rhoa(i0, j0)
      w%flw = w%flw(i0, j0)
      w%fsnow = w%fsnow(i0, j0)
      w%uvel = w%uvel(i0, j0)
      w%vvel = w%vvel(i0, j0)
      w%Tf = w%Tf(i0, j0)
      w%sst = w%sst(i0, j0)
      w%frzmlt = w%frzmlt(i0, j0)
      w%tau_x = w%tau_x(i0, j0)
      w%tau_y = w%tau_y(i0, j0)
      do n = 1, w%ncat
         w%aicen_init(:, :, n) = w%aicen_init(i0, j0, n)
         w%aicen(:, :, n) = w%aicen(i0, j0, n)
         w%vicen(:, :, n) = w%vicen(i0, j0, n)
         w%vsnon(:, :, n) = w%vsnon(i0, j0, n)
         w%Tsfcn(:, :, n) = w%Tsfcn(i0, j0, n)
         w%fswsfcn(:, :, n) = w%fswsfcn(i0, j0, n)
         w%fswintn(:, :, n) = w%fswintn(i0, j0, n)
         w%fswthrun(:, :, n) = w%fswthrun(i0, j0, n)
         w%qsnon(:, :, n, 1) = w%qsnon(i0, j0, n, 1)
         w%Sswabsn(:, :, n, 1) = w%Sswabsn(i0, j0, n, 1)
         do k = 1, w%nilyr
            w%qicen(:, :, n, k) = w%qicen(i0, j0, n, k)
            w%sicen(:, :, n, k) = w%sicen(i0, j0, n, k)
            w%Iswabsn(:, :, n, k) = w%Iswabsn(i0, j0, n, k)
         end do
      end do
      w%ice_err = 0_int32
   end subroutine make_smooth

   subroutine report_word_differences(a, b)
      !! Print the per-category error words of the first cells whose words differ between the two worlds (diagnosis aid).
      type(world_t), intent(in) :: a, b
      integer :: i, j, nshown
      nshown = 0
      do j = a%js, a%je
         do i = a%is, a%ie
            if (a%ice_err(i, j) /= b%ice_err(i, j) .and. nshown < 4) then
               nshown = nshown + 1
               write (output_unit, '(a,2i5,a,5i6,a,5i6,a,2es12.4)') "error words differ at ", i, j, ": device ", &
                  a%ierr_n(i, j, :), " host ", b%ierr_n(i, j, :), " potT, Qa ", a%potT(i, j), a%Qa(i, j)
            end if
         end do
      end do
   end subroutine report_word_differences

   subroutine test_device(error)
      !! 512 x 512 x 5 x 7: the production step on the device against plain host loops calling the same cell kernels, over an
      !! INTERIOR range; every output of every cell and category is compared and the untouched ring is asserted.  The error
      !! words are compared as `nonzero` and by the NON-FINITE bit: on the offload build a NaN may run further down the solver
      !! before it is caught (the relaxed-FP comparison hazard), so which additional FATAL bits follow it is not portable.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 512, NY = 512, NC = 5, NIL = 7
      type(world_t) :: w, h
      integer :: i, j, n
      integer(int64) :: st, t0, t1, rate
      real(wp) :: tdev, thost, ncell
      integer(int64) :: nbad
      integer :: rep, nexact

      call world_alloc(w, NX, NY, NC, NIL, ICE_CONSTANTS_ICEPACK)
      w%sw_redist = .true.
      w%is = 3
      w%ie = NX - 2
      w%js = 3
      w%je = NY - 2
      st = 20261113_int64
      do j = 1, NY
         do i = 1, NX
            call random_cell(w, i, j, st)
            ! one cell in 1500 and one category in 3000 get a NaN: error words must differ between cells too
            if (mod(i*7 + j*13, 1500) == 0) w%potT(i, j) = ieee_value(1.0_wp, ieee_quiet_nan)
            if (mod(i*11 + j*5, 3000) == 0) w%qicen(i, j, 2, 2) = ieee_value(1.0_wp, ieee_quiet_nan)
         end do
      end do
      w%aicen_init = w%aicen
      ! sentinels in every output so that an unwritten element is detected, and non-zero cell accumulators
      w%freshn = -777.0_wp; w%fsaltn = -777.0_wp; w%fhocnn = -777.0_wp; w%evapn = -777.0_wp; w%fsensn = -777.0_wp
      w%flatn = -777.0_wp; w%flwoutn = -777.0_wp; w%fswabsn = -777.0_wp; w%melttn = -777.0_wp; w%meltbn = -777.0_wp
      w%congeln = -777.0_wp; w%snoicen = -777.0_wp; w%strairxn = -555.0_wp; w%strairyn = -555.0_wp; w%ierr_n = 12345_int32
      w%strair_x_T = -444.0_wp; w%strair_y_T = -444.0_wp; w%fbot = -333.0_wp; w%rside = -333.0_wp; w%wlat = -333.0_wp
      do j = 1, NY
         do i = 1, NX
            w%fresh(i, j) = 1.0e-3_wp*real(1 + mod(i + 3*j, 17), wp)
         end do
      end do
      w%fsalt = 1.0e-5_wp; w%fhocn = 0.5_wp; w%evap = 1.0e-6_wp; w%fsens = -2.0_wp; w%flat = -1.0_wp; w%flwout = -250.0_wp
      w%fswabs = 3.0_wp; w%meltt = 1.0e-3_wp; w%meltb = 2.0e-3_wp; w%congel = 3.0e-3_wp; w%snoice = 4.0e-3_wp
      h = w

      ! host reference: plain loops over the same cell kernels
      call system_clock(t0, rate)
      call host_step(h)
      call system_clock(t1)
      thost = real(t1 - t0, wp)/real(rate, wp)

      call world_map(w)
      call system_clock(t0)
      call world_step(w)
      call system_clock(t1)
      tdev = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,f8.3,a,f8.3,a)') "therm1 step ", NX, "x", NY, "x", NC, "x", NIL, &
         ": first step incl. the transfers of every array ", tdev, " s; host loops ", thost, " s"
      nbad = count(h%ierr_n(3:NX - 2, 3:NY - 2, :) /= 0_int32)
      write (output_unit, '(a,i0,a,i0)') "error words nonzero: ", nbad, " of ", int(NX - 4, int64)*(NY - 4)*NC
      call check(error, nbad > 10 .and. nbad < 4000, "the state exercises the error words")
      call report_word_differences(w, h)
      if (.not. allocated(error)) call cmp_all(error, w, h, "device vs host")
      if (.not. allocated(error)) then
         ! the ring: untouched (still the sentinels), in every output
         call check(error, all(w%freshn(1:2, :, :) == -777.0_wp) .and. all(w%ierr_n(:, NY - 1:NY, :) == 12345_int32) .and. &
                    all(w%strair_x_T(:, 1:2) == -444.0_wp) .and. all(w%fbot(NX - 1:NX, :) == -333.0_wp), "the ring is untouched")
      end if
      ! the injected-NaN cells: the words are compared EXACTLY (every category and the cell), device against host
      if (.not. allocated(error)) then
         nexact = 0
         do j = 3, NY - 2
            do i = 3, NX - 2
               if (mod(i*7 + j*13, 1500) == 0 .or. mod(i*11 + j*5, 3000) == 0) then
                  nexact = nexact + 1
                  if (w%ice_err(i, j) /= h%ice_err(i, j) .or. any(w%ierr_n(i, j, :) /= h%ierr_n(i, j, :))) then
                     write (output_unit, '(a,2i5,a,i0,a,i0)') "injected-NaN cell ", i, j, ": device word ", w%ice_err(i, j), &
                        " host word ", h%ice_err(i, j)
                     call test_failed(error, "error words of an injected-NaN cell differ")
                     exit
                  end if
               end if
            end do
            if (allocated(error)) exit
         end do
       if (.not. allocated(error)) write (output_unit, '(a,i0,a)') "error words compared exactly on ", nexact, " injected-NaN cells"
      end if
      ! throughput with the data resident: the category launch alone and the whole four-launch step (state keeps evolving)
      if (.not. allocated(error)) then
         ncell = real(NX - 4, wp)*real(NY - 4, wp)*NC
         call system_clock(t0)
         do rep = 1, 5
            call ice_therm1_impl(w%nx, w%ny, w%ncat, w%nilyr, w%nslyr, w%is, w%ie, w%js, w%je, w%dt, w%p%c, w%p%atm, &
                                 w%p%thermo%conduct, w%p%thermo%l_brine, w%p%thermo%emissivity, &
                                 w%p%thermo%ice_ref_salinity, w%sw_redist, ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, &
                                 ICE_SFC_BULK, 0.0_wp, 0.0_wp, w%zlvl, w%zlvs, &
                                 w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel, w%Tf, &
                                 w%fbot, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                                 w%fswsfcn, w%fswintn, w%fswthrun, w%Iswabsn, w%Sswabsn, &
                                 w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, &
                                 w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, w%ice_err, w%ierr_n)
         end do
         !$acc wait
         call system_clock(t1)
         tdev = real(t1 - t0, wp)/real(rate, wp)/5.0_wp
         write (output_unit, '(a,f8.4,a,es10.3,a,f8.3,a,es10.3,a)') "category launch alone (data resident): ", tdev, " s (", &
            ncell/max(tdev, 1.0e-9_wp), " cell-cat/s); host loops ", thost, " s (", ncell/max(thost, 1.0e-9_wp), " cell-cat/s)"
      end if
      ! the same on a smooth, spatially coherent state (every cell the column of one clean cell): no divergence between
      ! neighbouring threads, so the upper end of what the launch can do
      if (.not. allocated(error)) then
         call make_smooth(w, NX/2, NY/2)
         call world_step(w)
         call system_clock(t0)
         do rep = 1, 5
            call ice_therm1_impl(w%nx, w%ny, w%ncat, w%nilyr, w%nslyr, w%is, w%ie, w%js, w%je, w%dt, w%p%c, w%p%atm, &
                                 w%p%thermo%conduct, w%p%thermo%l_brine, w%p%thermo%emissivity, &
                                 w%p%thermo%ice_ref_salinity, w%sw_redist, ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, &
                                 ICE_SFC_BULK, 0.0_wp, 0.0_wp, w%zlvl, w%zlvs, &
                                 w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel, w%Tf, &
                                 w%fbot, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                                 w%fswsfcn, w%fswintn, w%fswthrun, w%Iswabsn, w%Sswabsn, &
                                 w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, &
                                 w%melttn, w%meltbn, w%congeln, w%snoicen, w%strairxn, w%strairyn, w%ice_err, w%ierr_n)
         end do
         !$acc wait
         call system_clock(t1)
         tdev = real(t1 - t0, wp)/real(rate, wp)/5.0_wp
         write (output_unit, '(a,f8.4,a,es10.3,a,i0)') "category launch alone, smooth state: ", tdev, " s (", &
            ncell/max(tdev, 1.0e-9_wp), " cell-cat/s); error words nonzero after the first smooth step: ", &
            count(w%ierr_n(3:NX - 2, 3:NY - 2, :) /= 0_int32)
      end if
      call world_unmap(w)
   contains

      subroutine host_step(x)
         type(world_t), intent(inout) :: x
         integer :: ii, jj, nn, k
         real(wp) :: a, s(14)

         do jj = x%js, x%je
            do ii = x%is, x%ie
               call ice_frzmlt_bottom_lateral_cell(ii, jj, x%nx, x%ny, x%ncat, x%nilyr, x%nslyr, x%dt, x%p%thermo%ustar_min, &
                                                   x%p%thermo%floediam, x%p%thermo%floeshape, x%p%c%rhow, x%p%c%cp_ocn, &
                                                   x%p%c%pi, x%p%c%puny, x%aicen, x%vicen, x%vsnon, x%qicen, x%qsnon, &
                                                   x%frzmlt, x%sst, x%Tf, x%tau_x, x%tau_y, x%fbot, x%rside, x%wlat, x%ice_err)
            end do
         end do
         do nn = 1, x%ncat
            do jj = x%js, x%je
               do ii = x%is, x%ie
                  call ice_therm1_cell(ii, jj, nn, x%nx, x%ny, x%ncat, x%nilyr, x%nslyr, x%dt, x%p%c, x%p%atm, &
                                       x%p%thermo%conduct, x%p%thermo%l_brine, x%p%thermo%emissivity, &
                                       x%p%thermo%ice_ref_salinity, x%sw_redist, ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, &
                                       ICE_SFC_BULK, 0.0_wp, 0.0_wp, x%zlvl, x%zlvs, &
                                       x%uatm, x%vatm, x%wind, x%potT, x%Qa, x%rhoa, x%flw, x%fsnow, x%uvel, x%vvel, x%Tf, &
                                       x%fbot, x%aicen_init, x%aicen, x%vicen, x%vsnon, x%Tsfcn, x%qicen, x%sicen, x%qsnon, &
                                       x%fswsfcn, x%fswintn, x%fswthrun, x%Iswabsn, x%Sswabsn, &
                                       x%freshn, x%fsaltn, x%fhocnn, x%evapn, x%fsensn, x%flatn, x%flwoutn, x%fswabsn, &
                                       x%melttn, x%meltbn, x%congeln, x%snoicen, x%strairxn, x%strairyn, x%ice_err, x%ierr_n)
               end do
            end do
         end do
         do jj = x%js, x%je
            do ii = x%is, x%ie
               s = 0.0_wp
               do nn = 1, x%ncat
                  a = x%aicen_init(ii, jj, nn)
                  s(1) = s(1) + a*x%freshn(ii, jj, nn)
                  s(2) = s(2) + a*x%fsaltn(ii, jj, nn)
                  s(3) = s(3) + a*x%fhocnn(ii, jj, nn)
                  s(4) = s(4) + a*x%evapn(ii, jj, nn)
                  s(5) = s(5) + a*x%fsensn(ii, jj, nn)
                  s(6) = s(6) + a*x%flatn(ii, jj, nn)
                  s(7) = s(7) + a*(x%flwoutn(ii, jj, nn) - (1.0_wp - x%p%thermo%emissivity)*x%flw(ii, jj))
                  s(8) = s(8) + a*x%fswabsn(ii, jj, nn)
                  s(9) = s(9) + a*x%melttn(ii, jj, nn)
                  s(10) = s(10) + a*x%meltbn(ii, jj, nn)
                  s(11) = s(11) + a*x%congeln(ii, jj, nn)
                  s(12) = s(12) + a*x%snoicen(ii, jj, nn)
                  s(13) = s(13) + a*x%strairxn(ii, jj, nn)
                  s(14) = s(14) + a*x%strairyn(ii, jj, nn)
                  x%ice_err(ii, jj) = ior(x%ice_err(ii, jj), x%ierr_n(ii, jj, nn))
               end do
               x%fresh(ii, jj) = x%fresh(ii, jj) + s(1)
               x%fsalt(ii, jj) = x%fsalt(ii, jj) + s(2)
               x%fhocn(ii, jj) = x%fhocn(ii, jj) + s(3)
               x%evap(ii, jj) = x%evap(ii, jj) + s(4)
               x%fsens(ii, jj) = x%fsens(ii, jj) + s(5)
               x%flat(ii, jj) = x%flat(ii, jj) + s(6)
               x%flwout(ii, jj) = x%flwout(ii, jj) + s(7)
               x%fswabs(ii, jj) = x%fswabs(ii, jj) + s(8)
               x%meltt(ii, jj) = x%meltt(ii, jj) + s(9)
               x%meltb(ii, jj) = x%meltb(ii, jj) + s(10)
               x%congel(ii, jj) = x%congel(ii, jj) + s(11)
               x%snoice(ii, jj) = x%snoice(ii, jj) + s(12)
               x%strair_x_T(ii, jj) = s(13)
               x%strair_y_T(ii, jj) = s(14)
               k = 0
            end do
         end do
      end subroutine host_step

   end subroutine test_device

end module test_cice_therm1
