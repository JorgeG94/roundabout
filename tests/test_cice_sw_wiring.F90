!! The shortwave wiring (`rdb_ice_sw_wiring`, D7.4): radiation -> vertical thermodynamics, with a ZERO snow
!! absorption, and the `fswthru` / `fswthru_vis` producer.
!!
!! * `heat_closes_300_steps`     -- one 5-category column, 300 hourly steps of daytime shortwave with `sw_redist`, heat
!!                                  ledger `dE = dt (fswabs - fswthru + flwout + fsens - fhocn + Lvap evap) + H_ext` BY HAND
!!                                  (the shipped ledger type needs the ocean-side state slots; this test links the ice core
!!                                  only), residual against the BL99 solver bound `1.1 ferrmax dt sum(aicen_init)` exactly as
!!                                  `test_cice_therm1` asserts it; per step `fswabs - fswthru` against the hand sum over
!!                                  categories and layers of what the column absorbed (`fswsfcn + sum_k Iswabsn`);
!! * `zero_shortwave_bitwise`    -- `sw_down = 0` gives results bitwise equal to `ice_step_therm1` called with zero shortwave
!!                                  inputs, and `fswthru = fswthru_vis = 0`;
!! * `redist_conserves_absorbed` -- `sw_redist` on vs off: the absorbed total of every category and the cell `fswabs` agree to
!!                                  round-off, and the redistribution did move shortwave;
!! * `producer_hand_sum`         -- `fswthru` is the hand sum, additive on a pre-set value, `0 <= fswthru_vis <= fswthru`, and
!!                                  equal (all transmitted light is visible);
!! * `radiation_error_stops`     -- a refused radiation call stops the step (state, fluxes bitwise untouched, `istat = -1`); a
!!                                  non-finite `sw_down` flags and freezes that cell only;
!! * `device_vs_host_every_output`-- every output of every cell and category of a random world against the host cell kernels
!!                                  (interior), ring untouched; the snow absorption array is still zero on the device;
!! * `throughput`                -- 512 x 512 x 5 x 7 cell-categories per second of the whole sequence.
module test_cice_sw_wiring
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_bl99, only: ICE_NITERMAX, ICE_SW_FRAC, ICE_SW_DTEMP
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_sw_params_t, ICE_CONSTANTS_ICEPACK, ice_q_from_T, &
                             ice_qsno_from_T
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_radiation, only: ice_radiation_cell, ICE_RAD_OK, ICE_RAD_ERR_NOT_CONFIGURED
   use rdb_ice_sw_wiring, only: ice_step_sw_therm1, ice_sw_wiring_cleanup, ice_sw_wiring_zero_snow_ok
   use rdb_ice_therm1, only: ice_step_therm1, ice_therm1_cell
   use rdb_ice_thermo_shared, only: ICE_SFC_BULK
   use rdb_ice_thickness, only: ICE_FERRMAX
   implicit none
   private

   public :: collect_cice_sw_wiring_tests

   real(wp), parameter :: TWOPI = 6.283185307179586_wp, SNOWPATCH = 0.02_wp, SENT = -777.0_wp
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: DRTOL = 1.0e-12_wp
      !! device vs host, scale-aware (FMA on the GPU build, host code too)
   real(wp), parameter :: STOL = 1.0e-12_wp
      !! hand-sum / conservation tolerance relative to the incoming flux
#else
   real(wp), parameter :: DRTOL = 0.0_wp
      !! host: the same statements in the same order, measured worst difference 0
   real(wp), parameter :: STOL = 1.0e-13_wp
      !! measured worst hand-sum / conservation residual 2e-15 of the flux scale (about 50x)
#endif

   type :: world_t
      !! The arguments of `ice_step_sw_therm1` on a ghost-ringed `(nx, ny)` grid; `is:ie, js:je` the interior.
      integer :: nx = 0, ny = 0, ncat = 0, nil = 0, nsl = 1, is = 1, ie = 0, js = 1, je = 0
      real(wp) :: dt = 3600.0_wp
      logical :: redist = .true.
      type(ice_params_t) :: p
      type(ice_sw_params_t) :: swp
      type(ice_rad_state_t) :: rad
      real(wp), allocatable :: sw(:, :), uatm(:, :), vatm(:, :), wind(:, :), potT(:, :), Qa(:, :), rhoa(:, :), flw(:, :), &
                               fsnow(:, :), uvel(:, :), vvel(:, :), Tf(:, :), sst(:, :), frzmlt(:, :), tau_x(:, :), tau_y(:, :), &
                               fresh(:, :), fsalt(:, :), fhocn(:, :), evap(:, :), fsens(:, :), flat(:, :), flwout(:, :), &
                               fswabs(:, :), meltt(:, :), meltb(:, :), congel(:, :), snoice(:, :), sx(:, :), sy(:, :), &
                               fbot(:, :), rside(:, :), wlat(:, :), fswthru(:, :), fswthru_vis(:, :)
      integer(int32), allocatable :: ice_err(:, :), ierr_n(:, :, :)
      real(wp), allocatable :: aicen_init(:, :, :), aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :), &
                               freshn(:, :, :), fsaltn(:, :, :), fhocnn(:, :, :), evapn(:, :, :), fsensn(:, :, :), &
                               flatn(:, :, :), flwoutn(:, :, :), fswabsn(:, :, :), melttn(:, :, :), meltbn(:, :, :), &
                               congeln(:, :, :), snoicen(:, :, :), sxn(:, :, :), syn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
   end type world_t

contains

   subroutine collect_cice_sw_wiring_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("heat_closes_300_steps", test_closure), &
                  new_unittest("zero_shortwave_bitwise", test_zero_sw), &
                  new_unittest("redist_conserves_absorbed", test_redist), &
                  new_unittest("producer_hand_sum", test_producer), &
                  new_unittest("radiation_error_stops", test_errors), &
                  new_unittest("device_vs_host_every_output", test_device_vs_host), &
                  new_unittest("throughput", test_throughput) &
                  ]
   end subroutine collect_cice_sw_wiring_tests

   ! ------------------------------------------------------------------
   ! the world
   ! ------------------------------------------------------------------

   function uniform(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function uniform

   subroutine world_alloc(w, nx, ny, ncat, nil)
      type(world_t), intent(inout) :: w
      integer, intent(in) :: nx, ny, ncat, nil
      integer :: ierr

      w%nx = nx; w%ny = ny; w%ncat = ncat; w%nil = nil; w%nsl = 1
      w%is = merge(1, 2, nx == 1); w%ie = merge(1, nx - 1, nx == 1)
      w%js = merge(1, 2, ny == 1); w%je = merge(1, ny - 1, ny == 1)
      w%p = ice_params_default(ICE_CONSTANTS_ICEPACK, ncat, nil, 1)
      w%swp = ice_sw_params_t(albicev=0.86_wp, albicei=0.44_wp, albsnowv=0.98_wp, albsnowi=0.70_wp, ahmax=0.1_wp, &
                              sw_redist=.true., sw_frac=0.9_wp, sw_dtemp=0.02_wp, sw_vis_frac=0.52_wp)
      call w%rad%init(nx, ny, ncat, nil, ierr)
      call w%rad%configure(w%swp%ahmax, SNOWPATCH)
      allocate (w%sw(nx, ny), w%uatm(nx, ny), w%vatm(nx, ny), w%wind(nx, ny), w%potT(nx, ny), w%Qa(nx, ny), w%rhoa(nx, ny), &
                w%flw(nx, ny), w%fsnow(nx, ny), w%uvel(nx, ny), w%vvel(nx, ny), w%Tf(nx, ny), w%sst(nx, ny), w%frzmlt(nx, ny), &
                w%tau_x(nx, ny), w%tau_y(nx, ny), w%fresh(nx, ny), w%fsalt(nx, ny), w%fhocn(nx, ny), w%evap(nx, ny), &
                w%fsens(nx, ny), w%flat(nx, ny), w%flwout(nx, ny), w%fswabs(nx, ny), w%meltt(nx, ny), w%meltb(nx, ny), &
                w%congel(nx, ny), w%snoice(nx, ny), w%sx(nx, ny), w%sy(nx, ny), w%fbot(nx, ny), w%rside(nx, ny), &
                w%wlat(nx, ny), w%fswthru(nx, ny), w%fswthru_vis(nx, ny), w%ice_err(nx, ny), w%ierr_n(nx, ny, ncat), &
                w%aicen_init(nx, ny, ncat), w%aicen(nx, ny, ncat), w%vicen(nx, ny, ncat), w%vsnon(nx, ny, ncat), &
                w%Tsfcn(nx, ny, ncat), w%freshn(nx, ny, ncat), w%fsaltn(nx, ny, ncat), w%fhocnn(nx, ny, ncat), &
                w%evapn(nx, ny, ncat), w%fsensn(nx, ny, ncat), w%flatn(nx, ny, ncat), w%flwoutn(nx, ny, ncat), &
                w%fswabsn(nx, ny, ncat), w%melttn(nx, ny, ncat), w%meltbn(nx, ny, ncat), w%congeln(nx, ny, ncat), &
                w%snoicen(nx, ny, ncat), w%sxn(nx, ny, ncat), w%syn(nx, ny, ncat), w%qicen(nx, ny, ncat, nil), &
                w%sicen(nx, ny, ncat, nil), w%qsnon(nx, ny, ncat, 1))
      w%sw = 0.0_wp
      w%uatm = 6.0_wp; w%vatm = 0.0_wp; w%wind = 6.0_wp; w%potT = 263.0_wp; w%Qa = 6.0e-4_wp; w%rhoa = 1.3_wp
      w%flw = 200.0_wp; w%fsnow = 0.0_wp; w%uvel = 0.0_wp; w%vvel = 0.0_wp; w%Tf = -1.8_wp; w%sst = -1.8_wp
      w%frzmlt = 0.0_wp; w%tau_x = 0.03_wp; w%tau_y = 0.0_wp
      call world_zero_fluxes(w)
      w%sx = 0.0_wp; w%sy = 0.0_wp; w%fbot = 0.0_wp; w%rside = 0.0_wp; w%wlat = 0.0_wp
      w%ice_err = 0_int32; w%ierr_n = 0_int32
      w%freshn = SENT; w%fsaltn = SENT; w%fhocnn = SENT; w%evapn = SENT; w%fsensn = SENT; w%flatn = SENT
      w%flwoutn = SENT; w%fswabsn = SENT; w%melttn = SENT; w%meltbn = SENT; w%congeln = SENT; w%snoicen = SENT
      w%sxn = SENT; w%syn = SENT
      w%rad%fswsfcn = SENT; w%rad%fswintn = SENT; w%rad%fswthrun = SENT; w%rad%Iswabsn = SENT
      w%rad%alvn = SENT; w%rad%alin = SENT; w%rad%alv_ai = SENT; w%rad%ali_ai = SENT; w%rad%errn = 0_int32
   end subroutine world_alloc

   subroutine world_zero_fluxes(w)
      !! The per-step zeroing (`ice_flux_zero`) of the cell accumulators.
      type(world_t), intent(inout) :: w
      w%fresh = 0.0_wp; w%fsalt = 0.0_wp; w%fhocn = 0.0_wp; w%evap = 0.0_wp; w%fsens = 0.0_wp; w%flat = 0.0_wp
      w%flwout = 0.0_wp; w%fswabs = 0.0_wp; w%meltt = 0.0_wp; w%meltb = 0.0_wp; w%congel = 0.0_wp; w%snoice = 0.0_wp
      w%fswthru = 0.0_wp; w%fswthru_vis = 0.0_wp
   end subroutine world_zero_fluxes

   subroutine world_map(w)
      type(world_t), intent(inout) :: w
      !$acc enter data copyin(w%p, w%sw, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel) &
      !$acc& copyin(w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat) &
      !$acc& copyin(w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%sx, w%sy, w%fbot, w%rside, w%wlat) &
      !$acc& copyin(w%fswthru, w%fswthru_vis, w%ice_err, w%ierr_n, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn) &
      !$acc& copyin(w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, w%melttn) &
      !$acc& copyin(w%meltbn, w%congeln, w%snoicen, w%sxn, w%syn, w%qicen, w%sicen, w%qsnon)
      call w%rad%enter_data()
   end subroutine world_map

   subroutine world_unmap(w)
      type(world_t), intent(inout) :: w
      call w%rad%exit_data()
      !$acc exit data delete(w%p, w%sw, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel) &
      !$acc& delete(w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat) &
      !$acc& delete(w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%sx, w%sy, w%fbot, w%rside, w%wlat) &
      !$acc& delete(w%fswthru, w%fswthru_vis, w%ice_err, w%ierr_n, w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn) &
      !$acc& delete(w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, w%fswabsn, w%melttn) &
      !$acc& delete(w%meltbn, w%congeln, w%snoicen, w%sxn, w%syn, w%qicen, w%sicen, w%qsnon)
   end subroutine world_unmap

   subroutine world_push(w)
      !! Host inputs (and the zeroed accumulators) to the device.
      type(world_t), intent(inout) :: w
      !$acc update device(w%sw, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, w%vvel) &
      !$acc& device(w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat) &
      !$acc& device(w%flwout, w%fswabs, w%meltt, w%meltb, w%congel, w%snoice, w%fswthru, w%fswthru_vis, w%ice_err) &
      !$acc& device(w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon)
   end subroutine world_push

   subroutine world_pull(w)
      type(world_t), intent(inout) :: w
      !$acc update self(w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, w%meltt, w%meltb) &
      !$acc& self(w%congel, w%snoice, w%sx, w%sy, w%fbot, w%rside, w%wlat, w%fswthru, w%fswthru_vis, w%ice_err, w%ierr_n) &
      !$acc& self(w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%qsnon, w%freshn, w%fsaltn, w%fhocnn, w%evapn) &
      !$acc& self(w%fsensn, w%flatn, w%flwoutn, w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%sxn, w%syn)
      !$acc update self(w%rad%fswsfcn, w%rad%fswintn, w%rad%fswthrun, w%rad%Iswabsn, w%rad%alvn, w%rad%alin, w%rad%errn)
   end subroutine world_pull

   subroutine world_step(w, ierr, istat)
      !! One ice step on the MAPPED world: zero the accumulators, push, `ice_step_sw_therm1`, pull.
      type(world_t), intent(inout) :: w
      integer, intent(out) :: ierr, istat
      call world_zero_fluxes(w)
      call world_push(w)
      call ice_step_sw_therm1(w%rad, w%swp, w%nx, w%ny, w%ncat, w%nil, w%nsl, w%is, w%ie, w%js, w%je, w%dt, w%p, w%redist, &
                              10.0_wp, 10.0_wp, w%sw, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, w%fsnow, w%uvel, &
                              w%vvel, w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, &
                              w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                              w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, w%meltt, w%meltb, &
                              w%congel, w%snoice, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, w%flwoutn, &
                              w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%sxn, w%syn, w%sx, w%sy, &
                              w%fbot, w%rside, w%wlat, w%ierr_n, w%ice_err, w%fswthru, w%fswthru_vis, ierr, istat)
      call world_pull(w)
   end subroutine world_step

   subroutine world_random(w, nx, ny, ncat, nil, seed)
      !! Every cell and category differs: summer and winter, thin and thick, snow, near-melting layers, unlit cells.
      type(world_t), intent(inout) :: w
      integer, intent(in) :: nx, ny, ncat, nil
      integer(int64), intent(in) :: seed
      integer :: i, j, n, k
      integer(int64) :: st
      real(wp) :: u, rsum, r(ncat), hin, hsn, Tsfc, Sk, Tm, T, Tair
      logical :: summer

      call world_alloc(w, nx, ny, ncat, nil)
      st = seed
      do j = 1, ny
         do i = 1, nx
            summer = uniform(st) < 0.5_wp
            Tair = merge(271.5_wp + 4.0_wp*uniform(st), 245.0_wp + 30.0_wp*uniform(st), summer)
            w%potT(i, j) = Tair
            w%Qa(i, j) = 5.0e-4_wp + 3.0e-3_wp*uniform(st)
            w%rhoa(i, j) = 1.15_wp + 0.2_wp*uniform(st)
            u = 0.5_wp + 12.0_wp*uniform(st)
            w%uatm(i, j) = u*(2.0_wp*uniform(st) - 1.0_wp)
            w%vatm(i, j) = u*(2.0_wp*uniform(st) - 1.0_wp)
            w%wind(i, j) = sqrt(w%uatm(i, j)**2 + w%vatm(i, j)**2)
            w%flw(i, j) = 150.0_wp + 170.0_wp*uniform(st)
            w%fsnow(i, j) = merge(3.0e-4_wp*uniform(st), 0.0_wp, uniform(st) > 0.5_wp)
            w%Tf(i, j) = -1.9_wp + 0.3_wp*uniform(st)
            w%sst(i, j) = w%Tf(i, j) + 0.4_wp*uniform(st)
            w%frzmlt(i, j) = merge(-100.0_wp*uniform(st), 20.0_wp*uniform(st), uniform(st) < 0.5_wp)
            w%tau_x(i, j) = 0.2_wp*(2.0_wp*uniform(st) - 1.0_wp)
            w%tau_y(i, j) = 0.2_wp*(2.0_wp*uniform(st) - 1.0_wp)
            w%sw(i, j) = merge(0.0_wp, merge(150.0_wp + 200.0_wp*uniform(st), 80.0_wp*uniform(st), summer), uniform(st) < 0.15_wp)
            rsum = 0.0_wp
            do n = 1, ncat
               if (uniform(st) < 0.2_wp) then
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
            do n = 1, ncat
               w%aicen(i, j, n) = 0.92_wp*r(n)/rsum
               u = uniform(st)
               hin = 0.03_wp + 3.5_wp*u*u
               hsn = merge(0.3_wp*uniform(st), 0.0_wp, uniform(st) > 0.5_wp)
               w%vicen(i, j, n) = w%aicen(i, j, n)*hin
               w%vsnon(i, j, n) = w%aicen(i, j, n)*hsn
               Tsfc = merge(-0.05_wp - 0.5_wp*uniform(st), -0.05_wp - 30.0_wp*uniform(st), summer .and. uniform(st) < 0.6_wp)
               w%Tsfcn(i, j, n) = Tsfc
               do k = 1, nil
                  Sk = 0.3_wp + 5.5_wp*uniform(st)
                  w%sicen(i, j, n, k) = Sk
                  Tm = -0.054_wp*Sk
                  ! the upper layers of a summer category sit close to melting so the redistribution has work to do
                  T = Tm - 1.0e-3_wp - merge(0.05_wp, 3.0_wp, summer .and. k <= 2)*uniform(st)
                  w%qicen(i, j, n, k) = ice_q_from_T(T, Sk, w%p%c, w%p%thermo%l_brine)
               end do
               T = min(Tsfc*(0.2_wp + 0.8_wp*uniform(st)), -0.01_wp)
               w%qsnon(i, j, n, 1) = ice_qsno_from_T(T, w%p%c)
            end do
         end do
      end do
      w%aicen_init = w%aicen
   end subroutine world_random

   function worst_rel(n, a, b) result(worst)
      !! Largest `|a - b|` over the larger of the array's own scale and tiny; -1 if anything is non-finite.
      integer, intent(in) :: n
      real(wp), intent(in) :: a(*), b(*)
      real(wp) :: worst, sc
      integer :: m

      worst = 0.0_wp
      sc = 0.0_wp
      do m = 1, n
         if (.not. (ieee_is_finite(a(m)) .and. ieee_is_finite(b(m)))) then
            worst = -1.0_wp
            return
         end if
         sc = max(sc, abs(b(m)))
      end do
      do m = 1, n
         worst = max(worst, abs(a(m) - b(m)))
      end do
      worst = worst/max(sc, tiny(1.0_wp))
   end function worst_rel

   ! ------------------------------------------------------------------
   ! closure
   ! ------------------------------------------------------------------

   subroutine test_closure(error)
      !! 300 hourly steps of daytime shortwave on a melting multi-category column; the heat line by hand.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      real(wp), parameter :: A(5) = [0.05_wp, 0.15_wp, 0.30_wp, 0.25_wp, 0.10_wp]
      real(wp), parameter :: H(5) = [0.12_wp, 0.45_wp, 1.0_wp, 1.9_wp, 3.1_wp]
      real(wp), parameter :: HS(5) = [0.0_wp, 0.03_wp, 0.08_wp, 0.15_wp, 0.2_wp]
      integer, parameter :: NSTEP = 300, NIL = 4
      integer :: s, n, k, ierr, istat, nmelt, nthru
      real(wp) :: ph, Sp, Tz, E0, E1, a0, lw, th, sc, resid, ratio, worst_ratio, absorbed, dev, worst_dev
      real(wp) :: worst_h, thr_max
      real(wp), parameter :: LVAP = 2.501e6_wp

      call world_alloc(w, 1, 1, 5, NIL)
      do n = 1, 5
         w%aicen(1, 1, n) = A(n)
         w%vicen(1, 1, n) = A(n)*H(n)
         w%vsnon(1, 1, n) = A(n)*HS(n)
         w%Tsfcn(1, 1, n) = -3.0_wp
         do k = 1, NIL
            Sp = 4.0_wp - 1.0_wp*real(k - 1, wp)
            w%sicen(1, 1, n, k) = Sp
            Tz = -2.5_wp + 1.5_wp*real(k - 1, wp)/real(NIL - 1, wp)
            w%qicen(1, 1, n, k) = ice_q_from_T(Tz, Sp, w%p%c, w%p%thermo%l_brine)
         end do
         w%qsnon(1, 1, n, 1) = ice_qsno_from_T(-2.0_wp, w%p%c)
      end do
      w%aicen_init = w%aicen
      call world_map(w)
      nmelt = 0
      nthru = 0
      worst_ratio = 0.0_wp
      worst_h = 0.0_wp
      worst_dev = 0.0_wp
      thr_max = 0.0_wp
      do s = 1, NSTEP
         ph = TWOPI*real(s - 1, wp)/24.0_wp
         w%sw = 320.0_wp*max(0.0_wp, sin(ph))
         w%potT = 272.0_wp + 1.5_wp*sin(ph)
         w%flw = 0.85_wp*5.67e-8_wp*w%potT(1, 1)**4
         w%Qa = 0.8_wp*11637800.0_wp*exp(-5897.8_wp/w%potT(1, 1))/1.3_wp
         w%fsnow = merge(1.0e-4_wp, 0.0_wp, mod(s, 40) < 3)
         w%frzmlt = -15.0_wp
         w%sst = w%Tf + 0.1_wp
         w%aicen_init = w%aicen
         E0 = enthalpy(w)
         call world_step(w, ierr, istat)
         call check(error, ierr, ICE_RAD_OK, "radiation accepted")
         if (allocated(error)) exit
         call check(error, istat, 0, "therm1 accepted")
         if (allocated(error)) exit
         if (w%ice_err(1, 1) /= 0_int32) then
            call test_failed(error, "error word raised")
            exit
         end if
         E1 = enthalpy(w)
         a0 = 0.0_wp
         absorbed = 0.0_wp
         do n = 1, 5
            if (w%aicen_init(1, 1, n) > w%p%c%puny) then
               a0 = a0 + w%aicen_init(1, 1, n)
               absorbed = absorbed + w%aicen_init(1, 1, n)*(w%rad%fswsfcn(1, 1, n) + sum(w%rad%Iswabsn(1, 1, n, :)))
            end if
         end do
         ! what the column absorbed, summed by hand over categories and layers, against fswabs - fswthru
         dev = abs((w%fswabs(1, 1) - w%fswthru(1, 1)) - absorbed)/max(w%sw(1, 1), 1.0_wp)
         worst_dev = max(worst_dev, dev)
         if (.not. ieee_is_finite(dev)) worst_dev = huge(1.0_wp)
         lw = w%dt*(w%fswabs(1, 1) - w%fswthru(1, 1) + w%flwout(1, 1) + w%fsens(1, 1) - w%fhocn(1, 1) + LVAP*w%evap(1, 1))
         th = w%dt*(a0*w%flw(1, 1) - w%p%c%Lfresh*w%fsnow(1, 1)*a0)
         resid = (E1 - E0) - (lw + th)
         sc = abs(E1) + abs(E0) + w%dt*(abs(w%fswabs(1, 1)) + abs(w%fswthru(1, 1)) + abs(w%flwout(1, 1)) + abs(w%fsens(1, 1)) &
                                        + abs(w%fhocn(1, 1)) + a0*w%flw(1, 1) + w%p%c%Lfresh*w%fsnow(1, 1)*a0 &
                                        + LVAP*abs(w%evap(1, 1)))
         ratio = abs(resid)/(1.1_wp*ICE_FERRMAX*w%dt*a0 + 1.0e-15_wp*sc)
         if (.not. ieee_is_finite(ratio)) ratio = huge(1.0_wp)
         worst_ratio = max(worst_ratio, ratio)
         worst_h = max(worst_h, abs(resid)/sc)
         thr_max = max(thr_max, w%fswthru(1, 1))
         if (w%meltt(1, 1) > 0.0_wp) nmelt = nmelt + 1
         if (w%fswthru(1, 1) > 0.0_wp) nthru = nthru + 1
      end do
      call world_unmap(w)
      if (allocated(error)) return
   write (output_unit, '(a,es10.3,a,es10.3,a,es10.3,a,i0,a,i0,a,f7.2)') "sw wiring, 300 steps: worst heat residual / BL99 bound ", &
         worst_ratio, "; relative ", worst_h, "; worst |fswabs - fswthru - hand sum|/sw ", worst_dev, "; melting steps ", nmelt, &
         ", steps with fswthru ", nthru, "; max fswthru ", thr_max
      call check(error, nmelt > 20 .and. nthru > 100, "the run melts and transmits light")
      if (allocated(error)) return
      call check(error, worst_ratio <= 1.0_wp, "heat closes within the BL99 solver bound with shortwave on")
      if (allocated(error)) return
      call check(error, worst_dev <= STOL, "fswabs - fswthru is what the column absorbed")
   end subroutine test_closure

   function enthalpy(w) result(E)
      !! `E = sum_n (sum_k qicen vicen/nilyr + qsnon vsnon/nslyr)` of the single cell (1, 1) [J/m2].
      type(world_t), intent(in) :: w
      real(wp) :: E
      integer :: n, k

      E = 0.0_wp
      do n = 1, w%ncat
         do k = 1, w%nil
            E = E + w%qicen(1, 1, n, k)*w%vicen(1, 1, n)/real(w%nil, wp)
         end do
         E = E + w%qsnon(1, 1, n, 1)*w%vsnon(1, 1, n)
      end do
   end function enthalpy

   ! ------------------------------------------------------------------
   ! zero shortwave
   ! ------------------------------------------------------------------

   subroutine test_zero_sw(error)
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, r
      real(wp), allocatable :: z3(:, :, :), z3b(:, :, :), z3c(:, :, :), z4(:, :, :, :), s4(:, :, :, :)
      integer :: ierr, istat, ist2
      logical :: same

      call world_random(w, 6, 5, 4, 7, 4242_int64)
      w%sw = 0.0_wp
      r = w
      ! reference: ice_step_therm1 on zero shortwave inputs
      allocate (z3(w%nx, w%ny, w%ncat), z3b(w%nx, w%ny, w%ncat), z3c(w%nx, w%ny, w%ncat), source=0.0_wp)
      allocate (z4(w%nx, w%ny, w%ncat, w%nil), source=0.0_wp)
      allocate (s4(w%nx, w%ny, w%ncat, 1), source=0.0_wp)
      call world_map(w)
      call world_map(r)
      !$acc enter data copyin(z3, z3b, z3c, z4, s4)
      call world_step(w, ierr, istat)
      call world_zero_fluxes(r)
      call world_push(r)
      call ice_step_therm1(r%nx, r%ny, r%ncat, r%nil, 1, r%is, r%ie, r%js, r%je, r%dt, r%p, r%redist, 10.0_wp, 10.0_wp, &
                           r%uatm, r%vatm, r%wind, r%potT, r%Qa, r%rhoa, r%flw, r%fsnow, r%uvel, r%vvel, r%Tf, r%sst, &
                           r%frzmlt, r%tau_x, r%tau_y, r%aicen_init, r%aicen, r%vicen, r%vsnon, r%Tsfcn, r%qicen, r%sicen, &
                           r%qsnon, z3, z3b, z3c, z4, s4, r%fresh, r%fsalt, r%fhocn, r%evap, r%fsens, r%flat, r%flwout, r%fswabs, &
                           r%meltt, r%meltb, r%congel, r%snoice, r%freshn, r%fsaltn, r%fhocnn, r%evapn, r%fsensn, r%flatn, &
                           r%flwoutn, r%fswabsn, r%melttn, r%meltbn, r%congeln, r%snoicen, r%sxn, r%syn, r%sx, r%sy, r%fbot, &
                           r%rside, r%wlat, r%ierr_n, r%ice_err, ist2)
      call world_pull(r)
      !$acc exit data delete(z3, z3b, z3c, z4, s4)
      call world_unmap(w)
      call world_unmap(r)
      call check(error, ierr, ICE_RAD_OK, "radiation accepted")
      if (allocated(error)) return
      call check(error, istat + ist2, 0, "therm1 accepted both ways")
      if (allocated(error)) return
      same = all(bits3(w%aicen) == bits3(r%aicen)) .and. all(bits3(w%vicen) == bits3(r%vicen)) .and. &
             all(bits3(w%vsnon) == bits3(r%vsnon)) .and. all(bits3(w%Tsfcn) == bits3(r%Tsfcn)) .and. &
             all(bits3(w%fhocnn) == bits3(r%fhocnn)) .and. all(bits3(w%freshn) == bits3(r%freshn)) .and. &
             all(bits3(w%fswabsn) == bits3(r%fswabsn)) .and. all(bits3(w%flwoutn) == bits3(r%flwoutn))
      call check(error, same, "per-category state and fluxes bitwise equal to therm1 with zero shortwave inputs")
      if (allocated(error)) return
      same = all(bits2(w%fswabs) == bits2(r%fswabs)) .and. all(bits2(w%fhocn) == bits2(r%fhocn)) .and. &
             all(bits2(w%fresh) == bits2(r%fresh)) .and. all(bits2(w%fsens) == bits2(r%fsens)) .and. &
             all(bits2(w%meltt) == bits2(r%meltt))
      call check(error, same, "cell fluxes bitwise equal")
      if (allocated(error)) return
      same = all(abs(reshape(w%qicen, [size(w%qicen)]) - reshape(r%qicen, [size(r%qicen)])) == 0.0_wp)
      call check(error, same, "layer enthalpies bitwise equal")
      if (allocated(error)) return
      call check(error, all(w%fswthru == 0.0_wp) .and. all(w%fswthru_vis == 0.0_wp), "no transmitted shortwave")
      if (allocated(error)) return
      call check(error, all(w%rad%fswsfcn(w%is:w%ie, w%js:w%je, :) == 0.0_wp) .and. &
                 all(w%rad%Iswabsn(w%is:w%ie, w%js:w%je, :, :) == 0.0_wp), "radiation absorbs nothing")
   end subroutine test_zero_sw

   pure function bits2(a) result(b)
      real(wp), intent(in) :: a(:, :)
      integer(int64) :: b(size(a))
      b = transfer(a, b)
   end function bits2

   pure function bits3(a) result(b)
      real(wp), intent(in) :: a(:, :, :)
      integer(int64) :: b(size(a))
      b = transfer(a, b)
   end function bits3

   ! ------------------------------------------------------------------
   ! redistribution
   ! ------------------------------------------------------------------

   subroutine test_redist(error)
      !! `sw_redist` moves absorbed shortwave from the layers to the surface and must not change the total.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: on, off
      integer :: ierr, istat, i, j, n, nmoved, ncheck
      real(wp) :: tot_on, tot_off, worst, dfsw, worst_fswabs

      call world_random(on, 8, 6, 5, 7, 99_int64)
      off = on
      off%redist = .false.
      call world_map(on)
      call world_step(on, ierr, istat)
      call world_unmap(on)
      call world_map(off)
      call world_step(off, ierr, istat)
      call world_unmap(off)
      nmoved = 0
      ncheck = 0
      worst = 0.0_wp
      worst_fswabs = 0.0_wp
      do j = on%js, on%je
         do i = on%is, on%ie
            do n = 1, on%ncat
               if (on%aicen_init(i, j, n) <= on%p%c%puny .or. on%ierr_n(i, j, n) /= 0_int32) cycle
               ncheck = ncheck + 1
               tot_on = on%rad%fswsfcn(i, j, n) + sum(on%rad%Iswabsn(i, j, n, :))
               tot_off = off%rad%fswsfcn(i, j, n) + sum(off%rad%Iswabsn(i, j, n, :))
               if (.not. (ieee_is_finite(tot_on) .and. ieee_is_finite(tot_off))) then
                  call test_failed(error, "non-finite absorbed total")
                  return
               end if
               worst = max(worst, abs(tot_on - tot_off)/max(off%sw(i, j), 1.0_wp))
               dfsw = abs(on%rad%fswsfcn(i, j, n) - off%rad%fswsfcn(i, j, n))
               if (dfsw > 1.0e-6_wp) nmoved = nmoved + 1
               ! the per-category absorbed term of the step is fswsfcn + fswintn + fswthrun in both
               if (.not. (ieee_is_finite(on%fswabsn(i, j, n)) .and. ieee_is_finite(off%fswabsn(i, j, n)))) then
                  call test_failed(error, "non-finite fswabsn")
                  return
               end if
               worst_fswabs = max(worst_fswabs, abs(on%fswabsn(i, j, n) - off%fswabsn(i, j, n))/max(off%sw(i, j), 1.0_wp))
            end do
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,es10.3,a,es10.3)') "sw_redist on vs off: ", nmoved, " of ", ncheck, &
         " categories had shortwave moved; worst absorbed-total difference / sw ", worst, "; worst fswabsn difference / sw ", &
         worst_fswabs
      call check(error, nmoved >= 3, "the redistribution moved shortwave in this world")
      if (allocated(error)) return
      call check(error, worst <= STOL, "per-category absorbed total (surface + layers) unchanged by sw_redist")
      if (allocated(error)) return
      call check(error, worst_fswabs <= STOL, "fswabsn unchanged by sw_redist")
   end subroutine test_redist

   ! ------------------------------------------------------------------
   ! producer
   ! ------------------------------------------------------------------

   subroutine test_producer(error)
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, v
      integer :: ierr, istat, i, j, n, ncells, nthru
      real(wp) :: s, worst, pre
      logical :: ok_order, ok_eq

      call world_random(w, 8, 6, 5, 7, 31337_int64)
      call world_map(w)
      call world_step(w, ierr, istat)
      call world_unmap(w)
      call check(error, ierr + istat, 0, "step accepted")
      if (allocated(error)) return
      worst = 0.0_wp
      ok_order = .true.
      ok_eq = .true.
      ncells = 0
      nthru = 0
      do j = w%js, w%je
         do i = w%is, w%ie
            s = 0.0_wp
            do n = 1, w%ncat
               if (w%aicen_init(i, j, n) > w%p%c%puny .and. w%ierr_n(i, j, n) == 0_int32) &
                  s = s + w%aicen_init(i, j, n)*w%rad%fswthrun(i, j, n)
            end do
            if (.not. (ieee_is_finite(s) .and. ieee_is_finite(w%fswthru(i, j)) .and. ieee_is_finite(w%fswthru_vis(i, j)))) then
               call test_failed(error, "non-finite transmitted shortwave")
               return
            end if
            worst = max(worst, abs(w%fswthru(i, j) - s)/max(w%sw(i, j), 1.0_wp))
            ok_order = ok_order .and. w%fswthru_vis(i, j) <= w%fswthru(i, j) .and. w%fswthru(i, j) >= 0.0_wp &
                       .and. w%fswthru_vis(i, j) >= 0.0_wp
            ok_eq = ok_eq .and. w%fswthru_vis(i, j) == w%fswthru(i, j)
            ncells = ncells + 1
            if (s > 0.0_wp) nthru = nthru + 1
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,es10.3)') "producer: ", nthru, " of ", ncells, " cells transmit light; worst |fswthru - hand sum|/sw ", &
         worst
      call check(error, nthru >= ncells/3, "a third of the cells transmit")
      if (allocated(error)) return
      call check(error, worst <= STOL, "fswthru is the aicen_init-weighted hand sum of fswthrun")
      if (allocated(error)) return
      call check(error, ok_order, "0 <= fswthru_vis <= fswthru")
      if (allocated(error)) return
      call check(error, ok_eq, "fswthru_vis == fswthru (every transmitted watt is visible under CCSM3)")
      if (allocated(error)) return
      ! the ring is outside the owned range: untouched (zero accumulators)
      call check(error, w%fswthru(1, 1) == 0.0_wp .and. w%fswthru(w%nx, w%ny) == 0.0_wp, "ring untouched")
      if (allocated(error)) return
      ! additive: a second step on a pre-set accumulator keeps it (the world step zeroes first, so go through the shim directly)
      pre = 7.0_wp
      call world_random(v, 4, 4, 3, 4, 5_int64)
      call world_map(v)
      call world_zero_fluxes(v)
      v%fswthru = pre
      v%fswthru_vis = pre
      call world_push(v)
      call ice_step_sw_therm1(v%rad, v%swp, v%nx, v%ny, v%ncat, v%nil, v%nsl, v%is, v%ie, v%js, v%je, v%dt, v%p, v%redist, &
                              10.0_wp, 10.0_wp, v%sw, v%uatm, v%vatm, v%wind, v%potT, v%Qa, v%rhoa, v%flw, v%fsnow, v%uvel, &
                              v%vvel, v%Tf, v%sst, v%frzmlt, v%tau_x, v%tau_y, &
                              v%aicen_init, v%aicen, v%vicen, v%vsnon, v%Tsfcn, v%qicen, v%sicen, v%qsnon, &
                              v%fresh, v%fsalt, v%fhocn, v%evap, v%fsens, v%flat, v%flwout, v%fswabs, v%meltt, v%meltb, &
                              v%congel, v%snoice, v%freshn, v%fsaltn, v%fhocnn, v%evapn, v%fsensn, v%flatn, v%flwoutn, &
                              v%fswabsn, v%melttn, v%meltbn, v%congeln, v%snoicen, v%sxn, v%syn, v%sx, v%sy, &
                              v%fbot, v%rside, v%wlat, v%ierr_n, v%ice_err, v%fswthru, v%fswthru_vis, ierr, istat)
      call world_pull(v)
      call world_unmap(v)
      call check(error, ierr + istat, 0, "second step accepted")
      if (allocated(error)) return
      call check(error, all(v%fswthru(v%is:v%ie, v%js:v%je) >= pre) .and. all(v%fswthru(v%is:v%ie, v%js:v%je) < pre + 400.0_wp), &
                 "the producer ADDS to the accumulator")
      if (allocated(error)) return
      call check(error, v%fswthru(1, 1) == pre .and. v%fswthru_vis(1, 1) == pre, "ring still untouched (still the pre-set value)")
   end subroutine test_producer

   ! ------------------------------------------------------------------
   ! errors
   ! ------------------------------------------------------------------

   subroutine test_errors(error)
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, clean, ref
      integer :: ierr, istat, i, j, ib, jb
      logical :: same

      ib = 3
      jb = 3
      call world_random(w, 6, 5, 4, 4, 2718_int64)
      w%sw(ib, jb) = 250.0_wp   ! the flagged cell is lit and iced
      w%aicen(ib, jb, 1) = 0.4_wp
      w%vicen(ib, jb, 1) = 0.4_wp
      w%aicen_init = w%aicen
      clean = w
      ref = w
      call world_map(clean)
      call world_step(clean, ierr, istat)
      call world_unmap(clean)
      call check(error, ierr + istat, 0, "clean run accepted")
      if (allocated(error)) return
      ! (1) a refused radiation call stops the step: nothing runs, nothing is touched
      w%rad%configured = .false.
      call world_map(w)
      call world_step(w, ierr, istat)
      call world_unmap(w)
      call check(error, ierr, ICE_RAD_ERR_NOT_CONFIGURED, "unconfigured slot refused")
      if (allocated(error)) return
      call check(error, istat, -1, "therm1 not run")
      if (allocated(error)) return
      same = all(bits3(w%aicen) == bits3(ref%aicen)) .and. all(bits3(w%vicen) == bits3(ref%vicen)) .and. &
             all(bits3(w%Tsfcn) == bits3(ref%Tsfcn)) .and. all(w%fswthru == 0.0_wp) .and. all(w%fswabs == 0.0_wp) .and. &
             all(w%fresh == 0.0_wp)
      call check(error, same, "state and fluxes untouched after a refused radiation call")
      if (allocated(error)) return
      ! (2) a non-finite sw_down flags and freezes its own cell only
      w = ref
      w%sw(ib, jb) = ieee_value(1.0_wp, ieee_quiet_nan)
      call world_map(w)
      call world_step(w, ierr, istat)
      call world_unmap(w)
      call check(error, ierr + istat, 0, "step accepted (the flag is per cell)")
      if (allocated(error)) return
      call check(error, iand(w%ice_err(ib, jb), ICE_ERR_NONFINITE) /= 0_int32, "the corrupt cell carries the error bit")
      if (allocated(error)) return
      same = w%vicen(ib, jb, 1) == ref%vicen(ib, jb, 1) .and. w%aicen(ib, jb, 1) == ref%aicen(ib, jb, 1) .and. &
             w%Tsfcn(ib, jb, 1) == ref%Tsfcn(ib, jb, 1) .and. w%fswthru(ib, jb) == 0.0_wp .and. w%fswabs(ib, jb) == 0.0_wp
      call check(error, same, "the corrupt cell is not committed and transmits nothing")
      if (allocated(error)) return
      do j = w%js, w%je
         do i = w%is, w%ie
            if (i == ib .and. j == jb) cycle
            same = w%ice_err(i, j) == clean%ice_err(i, j) .and. w%fswthru(i, j) == clean%fswthru(i, j) .and. &
                   w%fswabs(i, j) == clean%fswabs(i, j) .and. all(w%vicen(i, j, :) == clean%vicen(i, j, :)) .and. &
                   all(w%Tsfcn(i, j, :) == clean%Tsfcn(i, j, :))
            call check(error, same, "the other cells are as in the clean run")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_errors

   ! ------------------------------------------------------------------
   ! device vs host
   ! ------------------------------------------------------------------

   subroutine test_device_vs_host(error)
      !! The shim against the host cell kernels (`ice_radiation_cell`, `ice_therm1_cell`) on a random world; the basal heat
      !! `fbot` (an independently tested launch) is taken from the shim's own output.  Whole arrays are compared, ring
      !! included, so an untouched ring is part of the check.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w, r
      real(wp), allocatable :: s4(:, :, :, :)
      real(wp) :: worst, wa(40)
      character(len=12) :: names(40)
      integer :: ierr, istat, i, j, n, na, nbad, u
      real(wp) :: fthru_ref, fabs_ref

      call world_random(w, 9, 7, 5, 7, 20261113_int64)
      r = w
      call world_map(w)
      call world_step(w, ierr, istat)
      call check(error, ierr + istat, 0, "step accepted")
      if (.not. allocated(error)) then
         call check(error, ice_sw_wiring_zero_snow_ok(), "the snow absorption is still exactly zero on the device")
      end if
      call world_unmap(w)
      if (allocated(error)) return
      ! ---- host reference
      allocate (s4(r%nx, r%ny, r%ncat, 1), source=0.0_wp)
      r%fbot = w%fbot
      r%rside = w%rside
      r%wlat = w%wlat
      do n = 1, r%ncat
         do j = r%js, r%je
            do i = r%is, r%ie
               call ice_radiation_cell(i, j, n, r%nx, r%ny, r%ncat, r%nil, r%swp%albicev, r%swp%albicei, r%swp%albsnowv, &
                                       r%swp%albsnowi, r%rad%fhtan, r%rad%snowpatch, r%swp%sw_vis_frac, r%aicen, r%vicen, &
                                       r%vsnon, r%Tsfcn, r%sw, r%rad%alvn, r%rad%alin, r%rad%fswsfcn, r%rad%fswintn, &
                                       r%rad%fswthrun, r%rad%Iswabsn, r%rad%errn)
            end do
         end do
      end do
      do j = r%js, r%je
         do i = r%is, r%ie
            r%ice_err(i, j) = ior(r%ice_err(i, j), ior(ior(r%rad%errn(i, j, 1), r%rad%errn(i, j, 2)), &
                                                       ior(r%rad%errn(i, j, 3), ior(r%rad%errn(i, j, 4), r%rad%errn(i, j, 5)))))
         end do
      end do
      do n = 1, r%ncat
         do j = r%js, r%je
            do i = r%is, r%ie
               call ice_therm1_cell(i, j, n, r%nx, r%ny, r%ncat, r%nil, 1, r%dt, r%p%c, r%p%atm, r%p%thermo%conduct, &
                                    r%p%thermo%l_brine, r%p%thermo%emissivity, r%p%thermo%ice_ref_salinity, r%redist, &
                                    ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, ICE_SFC_BULK, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, &
                                    r%uatm, r%vatm, r%wind, r%potT, r%Qa, r%rhoa, r%flw, r%fsnow, r%uvel, r%vvel, r%Tf, r%fbot, &
                                    r%aicen_init, r%aicen, r%vicen, r%vsnon, r%Tsfcn, r%qicen, r%sicen, r%qsnon, &
                                    r%rad%fswsfcn, r%rad%fswintn, r%rad%fswthrun, r%rad%Iswabsn, s4, &
                                    r%freshn, r%fsaltn, r%fhocnn, r%evapn, r%fsensn, r%flatn, r%flwoutn, r%fswabsn, &
                                    r%melttn, r%meltbn, r%congeln, r%snoicen, r%sxn, r%syn, r%ice_err, r%ierr_n)
            end do
         end do
      end do
      do j = r%js, r%je
         do i = r%is, r%ie
            fthru_ref = 0.0_wp
            fabs_ref = 0.0_wp
            do n = 1, r%ncat
               if (r%aicen_init(i, j, n) > r%p%c%puny) then
                  fabs_ref = fabs_ref + r%aicen_init(i, j, n)*r%fswabsn(i, j, n)
                  if (r%fswabsn(i, j, n) > 0.0_wp) fthru_ref = fthru_ref + r%aicen_init(i, j, n)*r%rad%fswthrun(i, j, n)
               end if
            end do
            r%fswabs(i, j) = fabs_ref
            r%fswthru(i, j) = fthru_ref
            r%fswthru_vis(i, j) = fthru_ref
         end do
      end do
      ! ---- compare every output
      na = 0
      call add(w%aicen, r%aicen, "aicen", w%nx*w%ny*w%ncat)
      call add(w%vicen, r%vicen, "vicen", w%nx*w%ny*w%ncat)
      call add(w%vsnon, r%vsnon, "vsnon", w%nx*w%ny*w%ncat)
      call add(w%Tsfcn, r%Tsfcn, "Tsfcn", w%nx*w%ny*w%ncat)
      call add(w%qicen, r%qicen, "qicen", w%nx*w%ny*w%ncat*w%nil)
      call add(w%qsnon, r%qsnon, "qsnon", w%nx*w%ny*w%ncat)
      call add(w%freshn, r%freshn, "freshn", w%nx*w%ny*w%ncat)
      call add(w%fsaltn, r%fsaltn, "fsaltn", w%nx*w%ny*w%ncat)
      call add(w%fhocnn, r%fhocnn, "fhocnn", w%nx*w%ny*w%ncat)
      call add(w%evapn, r%evapn, "evapn", w%nx*w%ny*w%ncat)
      call add(w%fsensn, r%fsensn, "fsensn", w%nx*w%ny*w%ncat)
      call add(w%flatn, r%flatn, "flatn", w%nx*w%ny*w%ncat)
      call add(w%flwoutn, r%flwoutn, "flwoutn", w%nx*w%ny*w%ncat)
      call add(w%fswabsn, r%fswabsn, "fswabsn", w%nx*w%ny*w%ncat)
      call add(w%melttn, r%melttn, "melttn", w%nx*w%ny*w%ncat)
      call add(w%meltbn, r%meltbn, "meltbn", w%nx*w%ny*w%ncat)
      call add(w%congeln, r%congeln, "congeln", w%nx*w%ny*w%ncat)
      call add(w%snoicen, r%snoicen, "snoicen", w%nx*w%ny*w%ncat)
      call add(w%sxn, r%sxn, "strairxn", w%nx*w%ny*w%ncat)
      call add(w%rad%alvn, r%rad%alvn, "alvn", w%nx*w%ny*w%ncat)
      call add(w%rad%alin, r%rad%alin, "alin", w%nx*w%ny*w%ncat)
      call add(w%rad%fswsfcn, r%rad%fswsfcn, "fswsfcn", w%nx*w%ny*w%ncat)
      call add(w%rad%fswintn, r%rad%fswintn, "fswintn", w%nx*w%ny*w%ncat)
      call add(w%rad%fswthrun, r%rad%fswthrun, "fswthrun", w%nx*w%ny*w%ncat)
      call add(w%rad%Iswabsn, r%rad%Iswabsn, "Iswabsn", w%nx*w%ny*w%ncat*w%nil)
      call add(w%fswthru, r%fswthru, "fswthru", w%nx*w%ny)
      call add(w%fswthru_vis, r%fswthru_vis, "fswthru_vis", w%nx*w%ny)
      call add(w%fswabs, r%fswabs, "fswabs", w%nx*w%ny)
      write (output_unit, '(a)', advance='no') "device vs host, worst relative difference per output:"
      worst = 0.0_wp
      nbad = 0
      do u = 1, na
         write (output_unit, '(1x,a,1x,es8.1)', advance='no') trim(names(u)), wa(u)
         if (wa(u) < 0.0_wp .or. wa(u) > DRTOL) nbad = nbad + 1
         worst = max(worst, wa(u))
      end do
      write (output_unit, '(a)') ""
      call check(error, nbad, 0, "every output of every cell and category matches the host cell kernels")
      if (allocated(error)) return
      call check(error, any(w%rad%fswthrun(w%is:w%ie, w%js:w%je, :) > 0.0_wp) .and. any(w%fswthru > 0.0_wp), &
                 "the world transmits light")
   contains
      subroutine add(a, b, nm, m)
         real(wp), intent(in) :: a(*), b(*)
         character(len=*), intent(in) :: nm
         integer, intent(in) :: m
         na = na + 1
         names(na) = nm
         wa(na) = worst_rel(m, a, b)
      end subroutine add
   end subroutine test_device_vs_host

   ! ------------------------------------------------------------------
   ! throughput
   ! ------------------------------------------------------------------

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: w
      integer, parameter :: NX = 512, NY = 512, NC = 5, NL = 7, NREP = 3
      integer :: ierr, istat, rep, i, j, n, k
      integer(int64) :: t0, t1, rate
      real(wp) :: secs, Sk, T

      call world_alloc(w, NX, NY, NC, NL)
      w%sw = 200.0_wp
      w%potT = 271.5_wp
      w%frzmlt = -10.0_wp
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               w%aicen(i, j, n) = 0.18_wp + 0.02_wp*real(mod(i + j, 3), wp)
               w%vicen(i, j, n) = w%aicen(i, j, n)*(0.2_wp + 0.5_wp*real(n, wp))
               w%vsnon(i, j, n) = w%aicen(i, j, n)*merge(0.0_wp, 0.05_wp*real(n, wp), n < 2)
               w%Tsfcn(i, j, n) = -1.0_wp
               do k = 1, NL
                  Sk = 3.0_wp
                  w%sicen(i, j, n, k) = Sk
                  T = -0.5_wp - 0.1_wp*real(k, wp)
                  w%qicen(i, j, n, k) = ice_q_from_T(T, Sk, w%p%c, w%p%thermo%l_brine)
               end do
               w%qsnon(i, j, n, 1) = ice_qsno_from_T(-1.0_wp, w%p%c)
            end do
         end do
      end do
      w%aicen_init = w%aicen
      call world_map(w)
      call world_step(w, ierr, istat)     ! warm-up: workspace, first launch
      call check(error, ierr + istat, 0, "warm-up step accepted")
      if (.not. allocated(error)) then
         call system_clock(t0, rate)
         do rep = 1, NREP
            call zero_fluxes_impl(w%nx, w%ny, w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, w%meltt, &
                                  w%meltb, w%congel, w%snoice, w%fswthru, w%fswthru_vis)
            call ice_step_sw_therm1(w%rad, w%swp, w%nx, w%ny, w%ncat, w%nil, w%nsl, w%is, w%ie, w%js, w%je, w%dt, w%p, &
                                    w%redist, 10.0_wp, 10.0_wp, w%sw, w%uatm, w%vatm, w%wind, w%potT, w%Qa, w%rhoa, w%flw, &
                                    w%fsnow, w%uvel, w%vvel, w%Tf, w%sst, w%frzmlt, w%tau_x, w%tau_y, &
                                    w%aicen_init, w%aicen, w%vicen, w%vsnon, w%Tsfcn, w%qicen, w%sicen, w%qsnon, &
                                    w%fresh, w%fsalt, w%fhocn, w%evap, w%fsens, w%flat, w%flwout, w%fswabs, w%meltt, &
                                    w%meltb, w%congel, w%snoice, w%freshn, w%fsaltn, w%fhocnn, w%evapn, w%fsensn, w%flatn, &
                                    w%flwoutn, w%fswabsn, w%melttn, w%meltbn, w%congeln, w%snoicen, w%sxn, w%syn, w%sx, &
                                    w%sy, w%fbot, w%rside, w%wlat, w%ierr_n, w%ice_err, w%fswthru, w%fswthru_vis, ierr, istat)
         end do
         !$acc update self(w%fswthru, w%ice_err)
         call system_clock(t1)
         secs = real(t1 - t0, wp)/real(rate, wp)
         write (output_unit, '(a,f8.1,a,es10.3,a)') "sw wiring 512x512x5x7: ", 1.0e3_wp*secs/real(NREP, wp), " ms/step, ", &
            real(NX*NY*NC, wp)*real(NREP, wp)/secs, " cell-categories/s (radiation + therm1 + producer)"
         call check(error, ierr + istat, 0, "timed steps accepted")
         if (.not. allocated(error)) then
           call check(error, ieee_is_finite(sum(w%fswthru(2:NX - 1, 2:NY - 1))) .and. sum(w%fswthru(2:NX - 1, 2:NY - 1)) > 0.0_wp, &
                       "finite, positive transmitted light")
         end if
      end if
      call world_unmap(w)
      call ice_sw_wiring_cleanup()
   end subroutine test_throughput

   subroutine zero_fluxes_impl(nx, ny, fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, &
                               snoice, fswthru, fswthru_vis)
      !! Zero the cell accumulators on the device (a launch, not a transfer); explicit-shape arrays.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), evap(nx, ny), fsens(nx, ny), flat(nx, ny)
      real(wp), intent(inout) :: flwout(nx, ny), fswabs(nx, ny), meltt(nx, ny), meltb(nx, ny), congel(nx, ny), snoice(nx, ny)
      real(wp), intent(inout) :: fswthru(nx, ny), fswthru_vis(nx, ny)
      integer :: i, j
      !$acc data present(fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, fswthru, fswthru_vis)
      do concurrent(j=1:ny, i=1:nx)
         fresh(i, j) = 0.0_wp; fsalt(i, j) = 0.0_wp; fhocn(i, j) = 0.0_wp; evap(i, j) = 0.0_wp
         fsens(i, j) = 0.0_wp; flat(i, j) = 0.0_wp; flwout(i, j) = 0.0_wp; fswabs(i, j) = 0.0_wp
         meltt(i, j) = 0.0_wp; meltb(i, j) = 0.0_wp; congel(i, j) = 0.0_wp; snoice(i, j) = 0.0_wp
         fswthru(i, j) = 0.0_wp; fswthru_vis(i, j) = 0.0_wp
      end do
      !$acc end data
   end subroutine zero_fluxes_impl

end module test_cice_sw_wiring
