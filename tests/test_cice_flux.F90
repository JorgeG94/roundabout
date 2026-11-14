!! `ice_flux_t` accumulators, `ice_flux_zero`, `ice_flux_merge_categories`,
!! `ice_save_init` and the ledger host mirror: lifecycle, device residency, every
!! array on every cell (ghosts included), brute-force host references, the
!! conservation ledger closing on a synthetic step with known fluxes, and the
!! non-finite paths.
!!
!! Device tests map the slots with `enter_data`, drive kernels that read and
!! write every array and compare CELL BY CELL with a host reference on states
!! whose cells all differ (a conserved total cannot see a wrong-cell read or
!! write).  Component arrays only are pulled back with `!$acc update self`.  On
!! host builds the directives are inert and the same code checks the kernels
!! against the references: EXACT on a CPU build (the merge sum order is fixed),
!! a stated tolerance on the offload build (FMA contraction is on there).
!! Set `RDB_CICE_BENCH=1` to also run the throughput case (512 x 512 x 5).
module test_cice_flux
   use, intrinsic :: iso_fortran_env, only: int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_grid, only: hgrid_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t, ice_flux_zero, ice_flux_merge_categories, ice_save_init, &
                                 ICE_FLUX_INIT_BAD_NCAT, ICE_FLUX_ERR_SHAPE
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   implicit none
   private

   public :: collect_cice_flux_tests

   integer, parameter :: NXP = 7, NYP = 5, NG = 2, NCAT = 5, NILYR = 7
   integer, parameter :: NX = NXP + 2*NG, NY = NYP + 2*NG
   integer, parameter :: NCELL = 20, NCATF = 12
      !! Cell accumulators and per-category outputs in `ice_flux_t`.
   integer, parameter :: NMERGE = NCATF
   integer, parameter :: MERGED(NMERGE) = [1, 2, 3, 6, 7, 8, 9, 10, 12, 13, 15, 16]
      !! Index in the cell list of the field each per-category output merges into.
   integer, parameter :: NKIND = 8
   real(wp), parameter :: DINJ(NKIND) = [1.0e-3_wp, 1.0e-5_wp, 1.0_wp, 0.0_wp, 1.0e-3_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      !! Injected imbalance (flux units) per `synthetic_step` kind.
   integer, parameter :: LINE_OF(NKIND) = [1, 2, 3, 0, 1, 3, 3, 0]
      !! Ledger line (1 water, 2 salt, 3 heat) each kind lands on; 0 = none.
   integer, parameter :: SGN_OF(NKIND) = [1, 1, 1, 0, -1, 1, -1, 0]
      !! Sign of the injected residual on that line.
   integer, parameter :: NB = 512
   real(wp), save :: FLW_BENCH(NB, NB)
   integer, parameter :: IS = NG + 1, IE = NG + NXP, JS = NG + 1, JE = NG + NYP
      !! Owned range.
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: TOL_REL = 1.0e-12_wp
      !! Offload build: FMA contraction on host and device (bounded by the measured
      !! worst deviation, reported with `RDB_CICE_BENCH=1`).
#else
   real(wp), parameter :: TOL_REL = 0.0_wp
      !! CPU build: exact (fixed summation order, no contraction).
#endif

contains

   subroutine collect_cice_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_destroy_bytes", test_init_destroy_bytes), &
                  new_unittest("init_refuses_bad_ncat", test_init_refuses), &
                  new_unittest("device_canary_every_array", test_device_canary), &
                  new_unittest("zero_clears_every_array_every_cell", test_zero), &
                  new_unittest("merge_random_cell_by_cell_additive", test_merge), &
                  new_unittest("merge_nonfinite_flags_only_its_cell", test_merge_nan), &
                  new_unittest("save_init_copies_exactly", test_save_init), &
                  new_unittest("merge_refuses_mismatched_slots", test_shape_mismatch), &
                  new_unittest("ledger_closes_icepack_sign_column", test_ledger_icepack), &
                  new_unittest("ledger_refuses_multirank_and_mismatch", test_ledger_refusals), &
                  new_unittest("ledger_closes_synthetic_steps", test_ledger_closes), &
                  new_unittest("ledger_reports_injected_imbalance", test_ledger_imbalance), &
                  new_unittest("ledger_nonfinite_and_unopened", test_ledger_nonfinite), &
                  new_unittest("throughput_512x512x5", test_throughput) &
                  ]
   end subroutine collect_cice_flux_tests

   ! ------------------------------------------------------------------
   subroutine make_slots(st, fx, ierr, nc, nxp_in, nyp_in)
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_flux_t), intent(inout) :: fx
      integer, intent(out) :: ierr
      integer, intent(in), optional :: nc, nxp_in, nyp_in
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      integer :: ncat_l, nxp_l, nyp_l

      ncat_l = NCAT
      nxp_l = NXP
      nyp_l = NYP
      if (present(nc)) ncat_l = nc
      if (present(nxp_in)) nxp_l = nxp_in
      if (present(nyp_in)) nyp_l = nyp_in
      call g%init(nxp_l, nyp_l, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, ncat_l, NILYR, 1)
      call st%init(g, p, ierr)
      if (ierr /= 0) return
      call fx%init(g, ncat_l, ierr)
   end subroutine make_slots

   function pt2(fx, k) result(p)
      !! Pointer to the k-th cell accumulator.
      type(ice_flux_t), target, intent(in) :: fx
      integer, intent(in) :: k
      real(wp), pointer :: p(:, :)
      p => null()
      select case (k)
      case (1)
         p => fx%fresh
      case (2)
         p => fx%fsalt
      case (3)
         p => fx%fhocn
      case (4)
         p => fx%fswthru
      case (5)
         p => fx%fswthru_vis
      case (6)
         p => fx%evap
      case (7)
         p => fx%fsens
      case (8)
         p => fx%flat
      case (9)
         p => fx%flwout
      case (10)
         p => fx%fswabs
      case (11)
         p => fx%frazil_vol
      case (12)
         p => fx%meltt
      case (13)
         p => fx%meltb
      case (14)
         p => fx%meltl
      case (15)
         p => fx%congel
      case (16)
         p => fx%snoice
      case (17)
         p => fx%dardg1dt
      case (18)
         p => fx%dardg2dt
      case (19)
         p => fx%dvirdgdt
      case (20)
         p => fx%opening
      case default
         continue
      end select
   end function pt2

   function pt3(fx, k) result(q)
      !! Pointer to the k-th per-category output.
      type(ice_flux_t), target, intent(in) :: fx
      integer, intent(in) :: k
      real(wp), pointer :: q(:, :, :)
      q => null()
      select case (k)
      case (1)
         q => fx%freshn
      case (2)
         q => fx%fsaltn
      case (3)
         q => fx%fhocnn
      case (4)
         q => fx%evapn
      case (5)
         q => fx%fsensn
      case (6)
         q => fx%flatn
      case (7)
         q => fx%flwoutn
      case (8)
         q => fx%fswabsn
      case (9)
         q => fx%melttn
      case (10)
         q => fx%meltbn
      case (11)
         q => fx%congeln
      case (12)
         q => fx%snoicen
      case default
         continue
      end select
   end function pt3

   subroutine pull_flux(fx)
      !! Component-array D->H pull (no-op on host builds).
      type(ice_flux_t), intent(inout) :: fx
      !$acc update self(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis)
      !$acc update self(fx%evap, fx%fsens, fx%flat, fx%flwout, fx%fswabs)
      !$acc update self(fx%frazil_vol, fx%meltt, fx%meltb, fx%meltl, fx%congel)
      !$acc update self(fx%snoice, fx%dardg1dt, fx%dardg2dt, fx%dvirdgdt, fx%opening)
      !$acc update self(fx%freshn, fx%fsaltn, fx%fhocnn, fx%evapn, fx%fsensn)
      !$acc update self(fx%flatn, fx%flwoutn, fx%fswabsn, fx%melttn, fx%meltbn)
      !$acc update self(fx%congeln, fx%snoicen)
   end subroutine pull_flux

   subroutine push_flux(fx)
      !! Component-array H->D push (no-op on host builds).
      type(ice_flux_t), intent(inout) :: fx
      !$acc update device(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis)
      !$acc update device(fx%evap, fx%fsens, fx%flat, fx%flwout, fx%fswabs)
      !$acc update device(fx%frazil_vol, fx%meltt, fx%meltb, fx%meltl, fx%congel)
      !$acc update device(fx%snoice, fx%dardg1dt, fx%dardg2dt, fx%dvirdgdt, fx%opening)
      !$acc update device(fx%freshn, fx%fsaltn, fx%fhocnn, fx%evapn, fx%fsensn)
      !$acc update device(fx%flatn, fx%flwoutn, fx%fswabsn, fx%melttn, fx%meltbn)
      !$acc update device(fx%congeln, fx%snoicen)
   end subroutine push_flux

   subroutine fill_flux_random(fx, seed, scale)
      !! Every array of the slot gets a per-cell different value (ghosts too).
      type(ice_flux_t), target, intent(inout) :: fx
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: scale
      integer(int64) :: s
      integer :: k, i, j, n
      real(wp), pointer :: p(:, :), q(:, :, :)

      s = seed
      do k = 1, NCELL
         p => pt2(fx, k)
         do j = 1, NY
            do i = 1, NX
               s = lcg_next(s)
               p(i, j) = scale*(lcg_unit(s) - 0.5_wp)*real(k, wp)
            end do
         end do
      end do
      do k = 1, NCATF
         q => pt3(fx, k)
         do n = 1, fx%ncat
            do j = 1, NY
               do i = 1, NX
                  s = lcg_next(s)
                  q(i, j, n) = scale*(lcg_unit(s) - 0.5_wp)*real(k, wp)
               end do
            end do
         end do
      end do
   end subroutine fill_flux_random

   subroutine same(error, got, want, what, rel)
      !! Finite-checked comparison: `|got-want| <= rel*max(|want|,1)`.
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: got, want, rel
      character(len=*), intent(in) :: what
      if (allocated(error)) return
      if (.not. ieee_is_finite(got)) then
         call test_failed(error, what//": non-finite value")
         return
      end if
      if (abs(got - want) > rel*max(abs(want), 1.0_wp)) call test_failed(error, what//": mismatch")
   end subroutine same

   ! ------------------------------------------------------------------
   subroutine test_init_destroy_bytes(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      integer :: ierr
      integer(int64) :: want

      call check(error, fx%bytes() == 0_int64, "bytes() of an uninitialised slot is 0")
      if (allocated(error)) return
      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init ok")
      if (allocated(error)) return
      call check(error, fx%is_init .and. .not. fx%mapped, "is_init, not mapped")
      call check(error, size(fx%fresh, 1) == NX .and. size(fx%fresh, 2) == NY, "cell shape is the ghosted extent")
      call check(error, size(fx%freshn, 1) == NX .and. size(fx%freshn, 2) == NY .and. size(fx%freshn, 3) == NCAT, &
                 "per-category shape (nx,ny,ncat)")
      call check(error, all(pt2(fx, 1) == 0.0_wp) .and. all(fx%snoicen == 0.0_wp), "inited to zero")
      want = int(NX*NY, int64)*8_int64*(int(NCELL, int64) + int(NCATF, int64)*int(NCAT, int64))
      call check(error, fx%bytes() == want, "bytes() counts every array (closed form)")
      call fx%destroy()
      call check(error,.not. fx%is_init .and. .not. allocated(fx%fresh) .and. fx%bytes() == 0_int64, &
                 "destroy releases everything")
   end subroutine test_init_destroy_bytes

   subroutine test_init_refuses(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_flux_t) :: fx
      type(hgrid_t) :: g
      integer :: ierr

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call fx%init(g, 0, ierr)
      call check(error, ierr, ICE_FLUX_INIT_BAD_NCAT, "ncat = 0 refused")
      call fx%init(g, 9, ierr)
      if (.not. allocated(error)) call check(error, ierr, ICE_FLUX_INIT_BAD_NCAT, "ncat beyond the bound refused")
      if (.not. allocated(error)) call check(error,.not. fx%is_init .and. .not. allocated(fx%fresh), "nothing allocated")
   end subroutine test_init_refuses

   ! ------------------------------------------------------------------
   subroutine canary_cell(nx_, ny_, fresh, fsalt, fhocn, fswthru, fswthru_vis, evap, &
                          fsens, flat, flwout, fswabs, frazil_vol, meltt, &
                          meltb, meltl, congel, snoice, dardg1dt, dardg2dt, &
                          dvirdgdt, opening)
      integer, intent(in) :: nx_, ny_
      real(wp), intent(inout) :: fresh(nx_, ny_)
      real(wp), intent(inout) :: fsalt(nx_, ny_)
      real(wp), intent(inout) :: fhocn(nx_, ny_)
      real(wp), intent(inout) :: fswthru(nx_, ny_)
      real(wp), intent(inout) :: fswthru_vis(nx_, ny_)
      real(wp), intent(inout) :: evap(nx_, ny_)
      real(wp), intent(inout) :: fsens(nx_, ny_)
      real(wp), intent(inout) :: flat(nx_, ny_)
      real(wp), intent(inout) :: flwout(nx_, ny_)
      real(wp), intent(inout) :: fswabs(nx_, ny_)
      real(wp), intent(inout) :: frazil_vol(nx_, ny_)
      real(wp), intent(inout) :: meltt(nx_, ny_)
      real(wp), intent(inout) :: meltb(nx_, ny_)
      real(wp), intent(inout) :: meltl(nx_, ny_)
      real(wp), intent(inout) :: congel(nx_, ny_)
      real(wp), intent(inout) :: snoice(nx_, ny_)
      real(wp), intent(inout) :: dardg1dt(nx_, ny_)
      real(wp), intent(inout) :: dardg2dt(nx_, ny_)
      real(wp), intent(inout) :: dvirdgdt(nx_, ny_)
      real(wp), intent(inout) :: opening(nx_, ny_)
      integer :: i, j

      !$acc data present(fresh, fsalt, fhocn, fswthru, fswthru_vis, evap) &
      !$acc& present(fsens, flat, flwout, fswabs, frazil_vol, meltt) &
      !$acc& present(meltb, meltl, congel, snoice, dardg1dt, dardg2dt) &
      !$acc& present(dvirdgdt, opening)
      do concurrent(j=1:ny_, i=1:nx_)
         fresh(i, j) = 2.0_wp*fresh(i, j) + 1.0_wp
         fsalt(i, j) = 2.0_wp*fsalt(i, j) + 2.0_wp
         fhocn(i, j) = 2.0_wp*fhocn(i, j) + 3.0_wp
         fswthru(i, j) = 2.0_wp*fswthru(i, j) + 4.0_wp
         fswthru_vis(i, j) = 2.0_wp*fswthru_vis(i, j) + 5.0_wp
         evap(i, j) = 2.0_wp*evap(i, j) + 6.0_wp
         fsens(i, j) = 2.0_wp*fsens(i, j) + 7.0_wp
         flat(i, j) = 2.0_wp*flat(i, j) + 8.0_wp
         flwout(i, j) = 2.0_wp*flwout(i, j) + 9.0_wp
         fswabs(i, j) = 2.0_wp*fswabs(i, j) + 10.0_wp
         frazil_vol(i, j) = 2.0_wp*frazil_vol(i, j) + 11.0_wp
         meltt(i, j) = 2.0_wp*meltt(i, j) + 12.0_wp
         meltb(i, j) = 2.0_wp*meltb(i, j) + 13.0_wp
         meltl(i, j) = 2.0_wp*meltl(i, j) + 14.0_wp
         congel(i, j) = 2.0_wp*congel(i, j) + 15.0_wp
         snoice(i, j) = 2.0_wp*snoice(i, j) + 16.0_wp
         dardg1dt(i, j) = 2.0_wp*dardg1dt(i, j) + 17.0_wp
         dardg2dt(i, j) = 2.0_wp*dardg2dt(i, j) + 18.0_wp
         dvirdgdt(i, j) = 2.0_wp*dvirdgdt(i, j) + 19.0_wp
         opening(i, j) = 2.0_wp*opening(i, j) + 20.0_wp
      end do
      !$acc end data
   end subroutine canary_cell

   subroutine canary_cat(nx_, ny_, nc, freshn, fsaltn, fhocnn, evapn, fsensn, flatn, &
                         flwoutn, fswabsn, melttn, meltbn, congeln, snoicen)
      integer, intent(in) :: nx_, ny_, nc
      real(wp), intent(inout) :: freshn(nx_, ny_, nc)
      real(wp), intent(inout) :: fsaltn(nx_, ny_, nc)
      real(wp), intent(inout) :: fhocnn(nx_, ny_, nc)
      real(wp), intent(inout) :: evapn(nx_, ny_, nc)
      real(wp), intent(inout) :: fsensn(nx_, ny_, nc)
      real(wp), intent(inout) :: flatn(nx_, ny_, nc)
      real(wp), intent(inout) :: flwoutn(nx_, ny_, nc)
      real(wp), intent(inout) :: fswabsn(nx_, ny_, nc)
      real(wp), intent(inout) :: melttn(nx_, ny_, nc)
      real(wp), intent(inout) :: meltbn(nx_, ny_, nc)
      real(wp), intent(inout) :: congeln(nx_, ny_, nc)
      real(wp), intent(inout) :: snoicen(nx_, ny_, nc)
      integer :: i, j, n

      !$acc data present(freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, melttn, meltbn, congeln, snoicen)
      do concurrent(n=1:nc, j=1:ny_, i=1:nx_)
         freshn(i, j, n) = 2.0_wp*freshn(i, j, n) + 1.0_wp
         fsaltn(i, j, n) = 2.0_wp*fsaltn(i, j, n) + 2.0_wp
         fhocnn(i, j, n) = 2.0_wp*fhocnn(i, j, n) + 3.0_wp
         evapn(i, j, n) = 2.0_wp*evapn(i, j, n) + 4.0_wp
         fsensn(i, j, n) = 2.0_wp*fsensn(i, j, n) + 5.0_wp
         flatn(i, j, n) = 2.0_wp*flatn(i, j, n) + 6.0_wp
         flwoutn(i, j, n) = 2.0_wp*flwoutn(i, j, n) + 7.0_wp
         fswabsn(i, j, n) = 2.0_wp*fswabsn(i, j, n) + 8.0_wp
         melttn(i, j, n) = 2.0_wp*melttn(i, j, n) + 9.0_wp
         meltbn(i, j, n) = 2.0_wp*meltbn(i, j, n) + 10.0_wp
         congeln(i, j, n) = 2.0_wp*congeln(i, j, n) + 11.0_wp
         snoicen(i, j, n) = 2.0_wp*snoicen(i, j, n) + 12.0_wp
      end do
      !$acc end data
   end subroutine canary_cat

   subroutine test_device_canary(error)
      !! Every array is read and written on the device, per-cell-different values; the
      !! result comes back through component `update self`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      real(wp), allocatable :: c0(:, :, :), k0(:, :, :, :)
      real(wp), pointer :: p(:, :), q(:, :, :)
      integer :: ierr, i, j, n, k

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_flux_random(fx, 20261113_int64, 4.0_wp)
      allocate (c0(NX, NY, NCELL), k0(NX, NY, NCAT, NCATF))
      do k = 1, NCELL
         c0(:, :, k) = pt2(fx, k)
      end do
      do k = 1, NCATF
         k0(:, :, :, k) = pt3(fx, k)
      end do
      call fx%enter_data()
      call check(error, fx%mapped, "enter_data maps")
      call canary_cell(NX, NY, fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis, fx%evap, &
                       fx%fsens, fx%flat, fx%flwout, fx%fswabs, fx%frazil_vol, fx%meltt, &
                       fx%meltb, fx%meltl, fx%congel, fx%snoice, fx%dardg1dt, fx%dardg2dt, &
                       fx%dvirdgdt, fx%opening)
      call canary_cat(NX, NY, NCAT, fx%freshn, fx%fsaltn, fx%fhocnn, fx%evapn, fx%fsensn, fx%flatn, &
                      fx%flwoutn, fx%fswabsn, fx%melttn, fx%meltbn, fx%congeln, fx%snoicen)
      call pull_flux(fx)
      do k = 1, NCELL
         p => pt2(fx, k)
         do j = 1, NY
            do i = 1, NX
               call same(error, p(i, j), 2.0_wp*c0(i, j, k) + real(k, wp), "cell accumulator canary", 1.0e-14_wp)
            end do
         end do
      end do
      do k = 1, NCATF
         q => pt3(fx, k)
         do n = 1, NCAT
            do j = 1, NY
               do i = 1, NX
                  call same(error, q(i, j, n), 2.0_wp*k0(i, j, n, k) + real(k, wp), "per-category canary", 1.0e-14_wp)
               end do
            end do
         end do
      end do
      call fx%exit_data()
      call check(error,.not. fx%mapped, "exit_data unmaps")
      call fx%destroy()
   end subroutine test_device_canary

   ! ------------------------------------------------------------------
   subroutine test_zero(error)
      !! Garbage in every array (including NaN and Inf on some cells), then ONE
      !! `ice_flux_zero`: every array, every cell, ghosts included, is exactly 0.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      real(wp), pointer :: p(:, :), q(:, :, :)
      integer :: ierr, i, j, n, k
      real(wp) :: keep_aicen

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_flux_random(fx, 7_int64, 100.0_wp)
      fx%fresh(3, 4) = ieee_value(1.0_wp, ieee_quiet_nan)
      fx%freshn(1, 1, 2) = ieee_value(1.0_wp, ieee_positive_inf)
      fx%opening(NX, NY) = ieee_value(1.0_wp, ieee_quiet_nan)
      fx%snoicen(NX, NY, NCAT) = ieee_value(1.0_wp, ieee_quiet_nan)
      st%aicen = 0.25_wp
      call st%enter_data()
      call fx%enter_data()
      call push_flux(fx)
      call ice_flux_zero(fx)
      call pull_flux(fx)
      do k = 1, NCELL
         p => pt2(fx, k)
         do j = 1, NY
            do i = 1, NX
               if (.not. (ieee_is_finite(p(i, j)) .and. p(i, j) == 0.0_wp)) then
                  call test_failed(error, "cell accumulator not exactly zero")
                  return
               end if
            end do
         end do
      end do
      do k = 1, NCATF
         q => pt3(fx, k)
         do n = 1, NCAT
            do j = 1, NY
               do i = 1, NX
                  if (.not. (ieee_is_finite(q(i, j, n)) .and. q(i, j, n) == 0.0_wp)) then
                     call test_failed(error, "per-category output not exactly zero")
                     return
                  end if
               end do
            end do
         end do
      end do
      keep_aicen = 0.25_wp
      !$acc update self(st%aicen)
      call check(error, all(st%aicen == keep_aicen), "ice_flux_zero touches no state array")
      call fx%exit_data()
      call st%exit_data()
      call fx%destroy()
      call st%destroy()
   end subroutine test_zero

   ! ------------------------------------------------------------------
   subroutine fill_state_random(st, seed, empty_frac)
      !! Per-cell different category state in every cell (ghosts too), some empty.
      type(ice_cice_state_t), intent(inout) :: st
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: empty_frac
      integer(int64) :: s
      integer :: i, j, n, k
      real(wp) :: u

      s = seed
      do n = 1, st%ncat
         do j = 1, st%ny
            do i = 1, st%nx
               s = lcg_next(s)
               u = lcg_unit(s)
               st%aicen(i, j, n) = 0.0_wp
               st%vicen(i, j, n) = 0.0_wp
               st%vsnon(i, j, n) = 0.0_wp
               st%Tsfcn(i, j, n) = -1.0_wp - 20.0_wp*u
               if (u > empty_frac) then
                  st%aicen(i, j, n) = 0.2_wp*u
                  st%vicen(i, j, n) = 3.0_wp*u*u
                  st%vsnon(i, j, n) = 0.4_wp*(1.0_wp - u)
               end if
               do k = 1, st%nilyr
                  s = lcg_next(s)
                  st%qicen(i, j, n, k) = -(1.0e8_wp + 2.0e8_wp*lcg_unit(s))
               end do
               s = lcg_next(s)
               st%qsnon(i, j, n, 1) = -(1.0e8_wp + 1.0e8_wp*lcg_unit(s))
            end do
         end do
      end do
      do j = 1, st%ny
         do i = 1, st%nx
            s = lcg_next(s)
            st%aice(i, j) = 0.1_wp + 0.8_wp*lcg_unit(s)
         end do
      end do
   end subroutine fill_state_random

   subroutine push_state(st)
      type(ice_cice_state_t), intent(inout) :: st
      !$acc update device(st%aicen, st%vicen, st%vsnon, st%Tsfcn)
      !$acc update device(st%qicen, st%sicen, st%qsnon)
      !$acc update device(st%aice0, st%aice, st%vice, st%vsno, st%ice_err)
   end subroutine push_state

   subroutine test_merge(error)
      !! The merge against a brute-force host sum in fixed category order, on a
      !! random per-cell-different state (some categories empty): every owned
      !! cell is the host value (exact on CPU), pre-existing cell values kept
      !! (additive), ghost cells untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      real(wp), allocatable :: pre(:, :, :), a_saved(:, :, :)
      real(wp), pointer :: p(:, :), q(:, :, :)
      real(wp) :: flw(NX, NY)
      integer :: ierr, i, j, n, k, ierr_m
      real(wp) :: sum_k, want, worst, emis

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_state_random(st, 99_int64, 0.3_wp)
      call fill_flux_random(fx, 1234_int64, 3.0_wp)
      allocate (pre(NX, NY, NCATF))
      do k = 1, NCATF
         pre(:, :, k) = pt2(fx, MERGED(k))
      end do
      do j = 1, NY
         do i = 1, NX
            flw(i, j) = 100.0_wp + 7.0_wp*real(i, wp) + 3.0_wp*real(j*j, wp)
         end do
      end do
      emis = st%p%thermo%emissivity
      st%aicen_init = -3.0_wp   ! host scratch; the device copy is what the kernel reads
      a_saved = st%aicen
      call st%enter_data()
      call fx%enter_data()
      !$acc enter data copyin(flw)
      call ice_save_init(st)
      ! change the LIVE aicen after the snapshot: a merge weighted by it would fail below
      st%aicen = 0.37_wp*a_saved + 0.01_wp
      !$acc update device(st%aicen)
      call ice_flux_merge_categories(fx, st, flw, ierr_m)
      call check(error, ierr_m, 0, "merge accepted matching slots")
      if (allocated(error)) return
      call pull_flux(fx)
      !$acc update self(st%ice_err)
      worst = 0.0_wp
      do k = 1, NCATF
         p => pt2(fx, MERGED(k))
         q => pt3(fx, k)
         do j = 1, NY
            do i = 1, NX
               want = pre(i, j, k)
               if (i >= IS .and. i <= IE .and. j >= JS .and. j <= JE) then
                  sum_k = 0.0_wp
                  do n = 1, NCAT
                     if (MERGED(k) == 9) then
                        sum_k = sum_k + a_saved(i, j, n)*(q(i, j, n) - (1.0_wp - emis)*flw(i, j))
                     else
                        sum_k = sum_k + a_saved(i, j, n)*q(i, j, n)
                     end if
                  end do
                  want = pre(i, j, k) + sum_k
               end if
               if (.not. ieee_is_finite(p(i, j))) then
                  call test_failed(error, "merged cell value not finite")
                  return
               end if
               worst = max(worst, abs(p(i, j) - want)/max(abs(want), 1.0_wp))
               call same(error, p(i, j), want, "merge vs brute force", TOL_REL)
               if (allocated(error)) return
            end do
         end do
      end do
      call report("merge worst relative deviation", worst)
      call check(error, all(st%ice_err == 0_int32), "finite input sets no error bit")
      !$acc exit data delete(flw)
      call fx%exit_data()
      call st%exit_data()
   end subroutine test_merge

   subroutine test_merge_nan(error)
      !! A NaN in one category output of one cell propagates into that cell's merged
      !! field (no clamp launders it) and sets `ICE_ERR_NONFINITE` there only.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      real(wp) :: flw(NX, NY)
      integer :: ierr, i, j, ierr_m
      integer, parameter :: IB = IS + 2, JB = JS + 1

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_state_random(st, 5_int64, 0.0_wp)
      call fill_flux_random(fx, 55_int64, 1.0_wp)
      fx%fsaltn(IB, JB, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
      st%aicen(IB, JB, 3) = 0.0_wp            ! zero weight: NaN * 0 must still poison
      st%aicen(IS, JS, 1) = 0.1_wp
      fx%melttn(IS, JS, 1) = ieee_value(1.0_wp, ieee_positive_inf)
      flw = 200.0_wp
      call st%enter_data()
      call fx%enter_data()
      !$acc enter data copyin(flw)
      call ice_save_init(st)
      call ice_flux_merge_categories(fx, st, flw, ierr_m)
      !$acc exit data delete(flw)
      call pull_flux(fx)
      !$acc update self(st%ice_err)
      call check(error,.not. ieee_is_finite(fx%fsalt(IB, JB)), "NaN reaches the cell accumulator")
      call check(error,.not. ieee_is_finite(fx%meltt(IS, JS)), "Inf reaches the cell accumulator")
      if (.not. allocated(error)) call check(error, iand(st%ice_err(IB, JB), ICE_ERR_NONFINITE) /= 0, "bit set (NaN cell)")
      if (.not. allocated(error)) call check(error, iand(st%ice_err(IS, JS), ICE_ERR_NONFINITE) /= 0, "bit set (Inf cell)")
      if (allocated(error)) return
      do j = JS, JE
         do i = IS, IE
            if ((i == IB .and. j == JB) .or. (i == IS .and. j == JS)) cycle
            if (st%ice_err(i, j) /= 0_int32) then
               call test_failed(error, "bit leaked to a finite cell")
               return
            end if
            if (.not. ieee_is_finite(fx%fsalt(i, j) + fx%meltt(i, j))) then
               call test_failed(error, "a finite cell became non-finite")
               return
            end if
         end do
      end do
      call fx%exit_data()
      call st%exit_data()
   end subroutine test_merge_nan

   subroutine test_save_init(error)
      !! The snapshots are exact copies (every cell, ghosts included) and the
      !! sources are untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      real(wp), allocatable :: a0(:, :, :), v0(:, :, :), s0(:, :, :), t0(:, :, :), ai0(:, :)
      integer :: ierr

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_state_random(st, 31_int64, 0.2_wp)
      a0 = st%aicen
      v0 = st%vicen
      s0 = st%vsnon
      t0 = st%Tsfcn
      ai0 = st%aice
      st%aicen_init = -7.0_wp   ! stale host values must not matter
      call st%enter_data()
      call ice_save_init(st)
      !$acc update self(st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%aice)
      !$acc update self(st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init, st%aice_init)
      call check(error, all(st%aicen_init == a0) .and. all(st%vicen_init == v0) .and. all(st%vsnon_init == s0) &
                 .and. all(st%Tsfcn_init == t0) .and. all(st%aice_init == ai0), "snapshots are exact copies, ghosts too")
      if (.not. allocated(error)) call check(error, all(st%aicen == a0) .and. all(st%vicen == v0) .and. &
                                             all(st%vsnon == s0) .and. all(st%Tsfcn == t0) .and. all(st%aice == ai0), &
                                             "sources unchanged")
      call st%exit_data()
   end subroutine test_save_init

   ! ------------------------------------------------------------------
   subroutine report(what, val)
      character(len=*), intent(in) :: what
      real(wp), intent(in) :: val
      character(len=8) :: e
      integer :: stat
      call get_environment_variable("RDB_CICE_BENCH", e, status=stat)
      if (stat == 0) write (error_unit, "(a,a,es10.3)") "  [cice_flux] ", what//": ", val
   end subroutine report

   subroutine synthetic_step(st, fx, led, area, seed, dt, kind, error)
      !! One ice step the TEST applies: it perturbs the category state and sets the
      !! cell flux accumulators so that the C.5 identities hold exactly (to
      !! round-off) on the OWNED cells, with the external sources it also hands to
      !! the ledger.  Ghost cells get unrelated garbage.  `kind` injects an
      !! imbalance: 0 none, 1 water (`fresh += d`), 2 salt (`fsalt += d`),
      !! 3 heat (`fhocn += d`), 4 NaN in `fresh`, 5 `evap += d`, 6 `fswthru += d`,
      !! 7 `fsens += d` (an atmosphere term), 8 `flat += d` (must change NOTHING).  `d = DINJ(kind)` on cell (IS+1, JS+1).
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_flux_t), intent(inout) :: fx
      type(ice_ledger_mirror_t), intent(inout) :: led
      real(wp), intent(in) :: area(NX, NY)
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: dt
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(inout) :: error
      integer(int64) :: s
      integer :: i, j, n, k
      real(wp) :: u, m0, m1, s0, s1, e0, e1, we, he, w_ext, h_ext, atm, rhoi, rhos, sref, vi0, vs0, vi1, vs1
      real(wp) :: flw_c

      rhoi = st%p%c%rhoi
      rhos = st%p%c%rhos
      sref = st%p%thermo%ice_ref_salinity
      call led%open(st, area, dt, 1)
      s = seed
      w_ext = 0.0_wp
      h_ext = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            vi0 = 0.0_wp
            vs0 = 0.0_wp
            e0 = 0.0_wp
            do n = 1, NCAT
               vi0 = vi0 + st%vicen(i, j, n)
               vs0 = vs0 + st%vsnon(i, j, n)
               do k = 1, NILYR
                  e0 = e0 + st%qicen(i, j, n, k)*st%vicen(i, j, n)/real(NILYR, wp)
               end do
               e0 = e0 + st%qsnon(i, j, n, 1)*st%vsnon(i, j, n)
            end do
            ! perturb the state (all categories), then the same sums after
            do n = 1, NCAT
               s = lcg_next(s)
               u = lcg_unit(s)
               st%vicen(i, j, n) = st%vicen(i, j, n)*(0.9_wp + 0.2_wp*u)
               s = lcg_next(s)
               st%vsnon(i, j, n) = st%vsnon(i, j, n)*(0.85_wp + 0.3_wp*lcg_unit(s))
               do k = 1, NILYR
                  s = lcg_next(s)
                  st%qicen(i, j, n, k) = st%qicen(i, j, n, k)*(0.95_wp + 0.1_wp*lcg_unit(s))
               end do
               s = lcg_next(s)
               st%qsnon(i, j, n, 1) = st%qsnon(i, j, n, 1)*(0.95_wp + 0.1_wp*lcg_unit(s))
            end do
            vi1 = 0.0_wp
            vs1 = 0.0_wp
            e1 = 0.0_wp
            do n = 1, NCAT
               vi1 = vi1 + st%vicen(i, j, n)
               vs1 = vs1 + st%vsnon(i, j, n)
               do k = 1, NILYR
                  e1 = e1 + st%qicen(i, j, n, k)*st%vicen(i, j, n)/real(NILYR, wp)
               end do
               e1 = e1 + st%qsnon(i, j, n, 1)*st%vsnon(i, j, n)
            end do
            m0 = rhoi*vi0 + rhos*vs0
            m1 = rhoi*vi1 + rhos*vs1
            s0 = rhoi*vi0*sref*1.0e-3_wp
            s1 = rhoi*vi1*sref*1.0e-3_wp
            ! fluxes: random atmosphere terms, closing identities for water, salt, ocean heat
            s = lcg_next(s)
            fx%evap(i, j) = 2.0e-5_wp*(lcg_unit(s) - 0.5_wp)   ! + condensation, - sublimation
            s = lcg_next(s)
            fx%fswabs(i, j) = 50.0_wp*lcg_unit(s)
            s = lcg_next(s)
            fx%flwout(i, j) = -300.0_wp*lcg_unit(s)   ! negative (Icepack)
            s = lcg_next(s)
            fx%fsens(i, j) = 20.0_wp*(lcg_unit(s) - 0.5_wp)
            s = lcg_next(s)
            fx%flat(i, j) = 30.0_wp*(lcg_unit(s) - 0.5_wp)
            s = lcg_next(s)
            fx%fswthru(i, j) = 10.0_wp*lcg_unit(s)
            s = lcg_next(s)
            we = 1.0e-4_wp*dt*lcg_unit(s)
            s = lcg_next(s)
            he = 1.0e4_wp*(lcg_unit(s) - 0.3_wp)
            s = lcg_next(s)
            flw_c = 100.0_wp + 200.0_wp*lcg_unit(s)
            he = he + dt*st%aice(i, j)*flw_c   ! incoming longwave on the ice fraction
            ! flat is not in the heat line; the vapour leaves with its own enthalpy: + Lvap*evap (upstream's final energy)
            atm = ((fx%fswabs(i, j) + fx%flwout(i, j)) + fx%fsens(i, j)) + st%p%c%Lvap*fx%evap(i, j)
            fx%fsalt(i, j) = -(s1 - s0)/dt
            fx%fresh(i, j) = (we - (m1 - m0))/dt + fx%evap(i, j)
            fx%fhocn(i, j) = atm - fx%fswthru(i, j) + (he - (e1 - e0))/dt
            if (i >= IS .and. i <= IE .and. j >= JS .and. j <= JE) then
               w_ext = w_ext + area(i, j)*we
               h_ext = h_ext + area(i, j)*he
               if (i == IS + 1 .and. j == JS + 1) then
                  select case (kind)
                  case (1)
                     fx%fresh(i, j) = fx%fresh(i, j) + DINJ(1)
                  case (2)
                     fx%fsalt(i, j) = fx%fsalt(i, j) + DINJ(2)
                  case (3)
                     fx%fhocn(i, j) = fx%fhocn(i, j) + DINJ(3)
                  case (4)
                     fx%fresh(i, j) = ieee_value(1.0_wp, ieee_quiet_nan)
                  case (5)
                     fx%evap(i, j) = fx%evap(i, j) + DINJ(5)
                  case (6)
                     fx%fswthru(i, j) = fx%fswthru(i, j) + DINJ(6)
                  case (7)
                     fx%fsens(i, j) = fx%fsens(i, j) + DINJ(7)
                  case (8)
                     fx%flat(i, j) = fx%flat(i, j) + DINJ(8)
                  case default
                     continue
                  end select
               end if
            else
               ! ghosts: unrelated garbage (the ledger must not read them)
               fx%fresh(i, j) = 1.0e3_wp
               fx%fsalt(i, j) = -1.0e3_wp
               fx%fhocn(i, j) = 1.0e6_wp
               st%vicen(i, j, 1) = st%vicen(i, j, 1) + 5.0_wp
            end if
         end do
      end do
      call push_state(st)
      call push_flux(fx)
      call led%close(st, fx, area, w_ext, h_ext)
      if (.not. led%valid) call test_failed(error, "ledger close produced no result")
   end subroutine synthetic_step

   subroutine setup_ledger_case(st, fx, area, error)
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_flux_t), intent(inout) :: fx
      real(wp), intent(inout) :: area(NX, NY)
      type(error_type), allocatable, intent(inout) :: error
      integer :: ierr, i, j

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_state_random(st, 4242_int64, 0.0_wp)
      do j = 1, NY
         do i = 1, NX
            area(i, j) = 1.0e6_wp*(1.0_wp + 0.01_wp*real(i, wp) + 0.003_wp*real(j*j, wp))   ! cells differ
         end do
      end do
      call st%enter_data()
      call fx%enter_data()
      !$acc enter data copyin(area)
   end subroutine setup_ledger_case

   subroutine finish_ledger_case(st, fx, area)
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_flux_t), intent(inout) :: fx
      real(wp), intent(inout) :: area(NX, NY)
      !$acc exit data delete(area)
      call fx%exit_data()
      call st%exit_data()
   end subroutine finish_ledger_case

   subroutine test_ledger_closes(error)
      !! Two synthetic steps with known fluxes and the matching state change: every
      !! C.5 line closes to round-off; the cumulative interface sums equal the host
      !! sums of `dt*A*flux` over both steps.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      type(ice_ledger_mirror_t) :: led
      real(wp) :: area(NX, NY)
      real(wp), parameter :: DT = 900.0_wp
      real(wp), parameter :: TOL = 1.0e-11_wp
         !! Measured worst 3.1e-13 (CPU): the line scale is built from NET totals, and
         !! cell changes of both signs cancel in them, so round-off of the per-cell
         !! values (1e-16 of the state) shows up at ~1e-13 of the scale; 1e-11 leaves
         !! a factor 30.  An injected imbalance is ~1e-6 of the scale.
      real(wp) :: cf, cs, ch, cw, worst
      integer :: step, i, j

      call setup_ledger_case(st, fx, area, error)
      if (allocated(error)) return
      cf = 0.0_wp
      cs = 0.0_wp
      ch = 0.0_wp
      cw = 0.0_wp
      do step = 1, 2
         call synthetic_step(st, fx, led, area, 600_int64 + step, DT, 0, error)
         if (allocated(error)) return
         call check(error, led%finite, "totals finite")
         worst = max(led%rel_residual(ICE_LEDGER_WATER), led%rel_residual(ICE_LEDGER_SALT), &
                     led%rel_residual(ICE_LEDGER_HEAT))
         call report("ledger worst relative residual", worst)
         call check(error, led%closed(TOL), "all three C.5 lines close")
         if (allocated(error)) return
         do j = JS, JE
            do i = IS, IE
               cf = cf + DT*area(i, j)*fx%fresh(i, j)
               cs = cs + DT*area(i, j)*fx%fsalt(i, j)
               ch = ch + DT*area(i, j)*fx%fhocn(i, j)
               cw = cw + DT*area(i, j)*fx%fswthru(i, j)
            end do
         end do
      end do
      call check(error, led%nsteps, 2, "two closed steps")
      call same(error, led%cum_fresh, cf, "cumulative fresh", 1.0e-12_wp)
      call same(error, led%cum_fsalt, cs, "cumulative fsalt", 1.0e-12_wp)
      call same(error, led%cum_fhocn, ch, "cumulative fhocn", 1.0e-12_wp)
      call same(error, led%cum_fswthru, cw, "cumulative fswthru", 1.0e-12_wp)
      call finish_ledger_case(st, fx, area)
   end subroutine test_ledger_closes

   subroutine test_ledger_imbalance(error)
      !! A deliberate imbalance on ONE owned cell appears as exactly the injected
      !! residual (signed) on its own line and only there, for the ocean terms and
      !! for `evap`, `fswthru` and the atmosphere term `fsens`, and `closed`
      !! fails; an injected `flat` changes nothing (it is not in the heat line).
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      type(ice_ledger_mirror_t) :: led
      real(wp) :: area(NX, NY)
      real(wp), parameter :: DT = 900.0_wp, TOL = 1.0e-11_wp
      real(wp) :: inj, want_res
      integer :: kind, line

      call setup_ledger_case(st, fx, area, error)
      if (allocated(error)) return
      do kind = 1, 8
         if (kind == 4) cycle
         inj = DT*area(IS + 1, JS + 1)*DINJ(kind)
         call synthetic_step(st, fx, led, area, 700_int64 + kind, DT, kind, error)
         if (allocated(error)) return
         if (LINE_OF(kind) == 0) then
            call check(error, led%closed(TOL), "an injected flat leaves every line closed")
         else
            call check(error,.not. led%closed(TOL), "closure fails with an injected imbalance")
         end if
         if (allocated(error)) return
         do line = 1, 3
            if (line == LINE_OF(kind) .or. (kind == 5 .and. line == 3)) then
               ! an injected `evap` (kind 5) is also vapour enthalpy: it lands on the heat line as -Lvap*d*A*dt
               want_res = real(SGN_OF(kind), wp)*inj
               if (kind == 5 .and. line == 3) want_res = -st%p%c%Lvap*inj
               if (.not. ieee_is_finite(led%res(line))) then
                  call test_failed(error, "injected line non-finite")
               else if (abs(led%res(line) - want_res) > 1.0e-9_wp*led%scale(line)) then
                  call test_failed(error, "residual is not the signed injected d*A*dt")
               end if
            else
               if (.not. ieee_is_finite(led%rel_residual(line))) then
                  call test_failed(error, "other line non-finite")
               else if (led%rel_residual(line) > TOL) then
                  call test_failed(error, "imbalance leaked to another line")
               end if
            end if
            if (allocated(error)) return
         end do
      end do
      call finish_ledger_case(st, fx, area)
   end subroutine test_ledger_imbalance

   subroutine test_shape_mismatch(error)
      !! The merge refuses slots that disagree on `ncat` (and a wrong-shape `flw`):
      !! `ICE_FLUX_ERR_SHAPE`, no kernel run, nothing written.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx, fx_bad
      type(hgrid_t) :: g
      real(wp) :: flw(NX, NY), flw_small(NX - 1, NY)
      integer :: ierr

      call make_slots(st, fx, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call fx_bad%init(g, NCAT - 1, ierr)
      call check(error, ierr, 0, "init (other ncat)")
      if (allocated(error)) return
      fx_bad%fresh = 1.0_wp
      flw = 0.0_wp
      flw_small = 0.0_wp
      call ice_flux_merge_categories(fx_bad, st, flw, ierr)
      call check(error, ierr, ICE_FLUX_ERR_SHAPE, "ncat mismatch refused")
      if (.not. allocated(error)) call check(error, all(fx_bad%fresh == 1.0_wp), "nothing written on refusal")
      if (.not. allocated(error)) then
         call ice_flux_merge_categories(fx, st, flw_small, ierr)
         call check(error, ierr, ICE_FLUX_ERR_SHAPE, "wrong-shape flw refused")
      end if
   end subroutine test_shape_mismatch

   subroutine test_ledger_refusals(error)
      !! More than one rank, or a flux slot that disagrees with the state: the
      !! ledger is invalid (no result, NaN residuals, never closed).
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx, fx_bad
      type(ice_ledger_mirror_t) :: led
      type(hgrid_t) :: g
      real(wp) :: area(NX, NY)
      integer :: ierr

      call setup_ledger_case(st, fx, area, error)
      if (allocated(error)) return
      call led%open(st, area, 900.0_wp, 2)
      call led%close(st, fx, area)
      call check(error,.not. led%valid .and. .not. led%closed(huge(1.0_wp)) .and. &
                 .not. ieee_is_finite(led%rel_residual(ICE_LEDGER_SALT)), "two ranks: ledger refused")
      if (allocated(error)) return
      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call fx_bad%init(g, NCAT - 1, ierr)
      call led%open(st, area, 900.0_wp, 1)
      call led%close(st, fx_bad, area)
      call check(error,.not. led%valid .and. .not. led%closed(huge(1.0_wp)), "mismatched flux slot: ledger refused")
      call finish_ledger_case(st, fx, area)
   end subroutine test_ledger_refusals

   subroutine test_ledger_icepack(error)
      !! Hand-built Icepack-sign columns, closed with the UPSTREAM expressions and
      !! not with the ledger's own equation.  Per category and cell:
      !!   mass:   d(rhoi*hi + rhos*hs) = (fsnow + evapn - freshn)*dt   (snow takes it),
      !!           with pure sublimation (`evapn < 0`, `freshn = 0`) in a third of the cells;
      !!   energy: dE_state = einp + Lvap*evapn*dt, einp = (fsurfn - flatn + fswint - fhocnn - fsnow*Lfresh)*dt,
      !!           (upstream's final energy carries the vapour at -Lvap*evapn*dt),
      !!           fsurfn = fswsfc + emis*flw + flwoutn + fsensn + flatn   (flwoutn < 0,
      !!           fsensn, flatn positive downward), fswabsn = fswsfc + fswint + fswthrun.
      !! The cell fluxes come from the real merge (`flw` reflected part included);
      !! the caller-side sources are `W_ext = dt*A*aice*fsnow` and
      !! `H_ext = dt*A*aice*(flw - Lfresh*fsnow)`.  ICE_NSLYR_MAX = 1, so the
      !! energy sum with more than one snow layer cannot be exercised yet.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t), target :: fx
      type(ice_ledger_mirror_t) :: led
      real(wp) :: area(NX, NY), flw(NX, NY)
      real(wp), parameter :: DT = 900.0_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
         !! Measured worst 2.6e-11: the snow enthalpy is rebuilt as (E + dE)/v, so the
         !! round-off of E (1e8 J/m3 scale) is amplified by E/dE ~ 1e3 against the
         !! line scale.  A wrong sign or a stray `flat` term is ~1e-1 of the scale.
      integer(int64) :: sd
      integer :: i, j, n, ierr, ierr_m
      real(wp) :: emis, lf, rhos, fsnow, ai, w_ext, h_ext, u
      real(wp) :: fswsfc, fswint, fswthrun, flwoutn, fsensn, flatn, fhocnn, evapn, freshn, fsurfn, einp, dm
      real(wp) :: fswthru_cell, aicen_n, vs0, e_n
      real(wp) :: a_save(NX, NY, NCAT)

      call setup_ledger_case(st, fx, area, error)
      if (allocated(error)) return
      emis = st%p%thermo%emissivity
      lf = st%p%c%Lfresh
      rhos = st%p%c%rhos
      sd = 2718_int64
      do j = 1, NY
         do i = 1, NX
            ai = 0.0_wp
            do n = 1, NCAT
               ai = ai + st%aicen(i, j, n)
               st%vsnon(i, j, n) = max(st%vsnon(i, j, n), 0.05_wp)
            end do
            st%aice(i, j) = ai
            sd = lcg_next(sd)
            flw(i, j) = 150.0_wp + 150.0_wp*lcg_unit(sd)
         end do
      end do
      a_save = st%aicen
      call push_state(st)
      !$acc enter data copyin(flw)
      call ice_flux_zero(fx)
      call ice_save_init(st)
      call led%open(st, area, DT, 1)
      w_ext = 0.0_wp
      h_ext = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            sd = lcg_next(sd)
            fsnow = 2.0e-5_wp*lcg_unit(sd)
            ai = st%aice(i, j)
            fswthru_cell = 0.0_wp
            do n = 1, NCAT
               aicen_n = a_save(i, j, n)
               vs0 = st%vsnon(i, j, n)
               sd = lcg_next(sd)
               fswsfc = 80.0_wp*lcg_unit(sd)
               sd = lcg_next(sd)
               fswint = 20.0_wp*lcg_unit(sd)
               sd = lcg_next(sd)
               fswthrun = 10.0_wp*lcg_unit(sd)
               sd = lcg_next(sd)
               flwoutn = -(250.0_wp + 150.0_wp*lcg_unit(sd))
               sd = lcg_next(sd)
               fsensn = 30.0_wp*(lcg_unit(sd) - 0.5_wp)
               sd = lcg_next(sd)
               flatn = 40.0_wp*(lcg_unit(sd) - 0.5_wp)
               sd = lcg_next(sd)
               fhocnn = 60.0_wp*lcg_unit(sd)
               sd = lcg_next(sd)
               u = lcg_unit(sd)
               sd = lcg_next(sd)
               evapn = -1.0e-5_wp*lcg_unit(sd)
               if (u < 0.33_wp) then
                  freshn = 0.0_wp                 ! pure sublimation: fresh = 0, mass is lost
               else
                  sd = lcg_next(sd)
                  freshn = 1.0e-5_wp*lcg_unit(sd)
               end if
               fsurfn = fswsfc + emis*flw(i, j) + flwoutn + fsensn + flatn
               einp = (fsurfn - flatn + fswint - fhocnn - fsnow*lf)*DT
               dm = (fsnow + evapn - freshn)*DT          ! per ice area [kg m-2]
               e_n = st%qsnon(i, j, n, 1)*vs0
               st%vsnon(i, j, n) = vs0 + aicen_n*dm/rhos
               ! upstream: the energy check closes ice + vapour, efinal = -evap*Lvap + sum(h*q), so the STORED
               ! enthalpy changes by einp + Lvap*evap*dt (evap = evapn*dt in kg m-2)
               st%qsnon(i, j, n, 1) = (e_n + aicen_n*(einp + st%p%c%Lvap*evapn*DT))/st%vsnon(i, j, n)
               fx%freshn(i, j, n) = freshn
               fx%evapn(i, j, n) = evapn
               fx%fhocnn(i, j, n) = fhocnn
               fx%fsensn(i, j, n) = fsensn
               fx%flatn(i, j, n) = flatn
               fx%flwoutn(i, j, n) = flwoutn
               fx%fswabsn(i, j, n) = fswsfc + fswint + fswthrun
               fswthru_cell = fswthru_cell + aicen_n*fswthrun
            end do
            fx%fswthru(i, j) = fswthru_cell
            if (i >= IS .and. i <= IE .and. j >= JS .and. j <= JE) then
               w_ext = w_ext + DT*area(i, j)*ai*fsnow
               h_ext = h_ext + DT*area(i, j)*ai*(flw(i, j) - lf*fsnow)
            end if
         end do
      end do
      call push_state(st)
      call push_flux(fx)
      call ice_flux_merge_categories(fx, st, flw, ierr_m)
      call check(error, ierr_m, 0, "merge ok")
      if (allocated(error)) return
      call led%close(st, fx, area, w_ext, h_ext)
      call report("icepack-sign column worst relative residual", &
                  max(led%rel_residual(ICE_LEDGER_WATER), led%rel_residual(ICE_LEDGER_HEAT)))
      call check(error, led%valid .and. led%closed(TOL), "water and heat close with upstream sign expressions")
      ierr = 0
      !$acc exit data delete(flw)
      call finish_ledger_case(st, fx, area)
   end subroutine test_ledger_icepack

   subroutine test_ledger_nonfinite(error)
      !! A NaN flux makes the ledger non-finite and not closed; a `close` without
      !! an `open` yields no result.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      type(ice_ledger_mirror_t) :: led
      real(wp) :: area(NX, NY)

      call setup_ledger_case(st, fx, area, error)
      if (allocated(error)) return
      call led%close(st, fx, area)
      call check(error,.not. led%valid .and. .not. led%closed(1.0_wp), "close without open is not a result")
      if (allocated(error)) return
      call synthetic_step(st, fx, led, area, 800_int64, 900.0_wp, 4, error)
      if (allocated(error)) return
      call check(error,.not. led%finite, "NaN flux is reported")
      if (.not. allocated(error)) call check(error,.not. led%closed(huge(1.0_wp)), "never closed on NaN (even at huge tol)")
      if (.not. allocated(error)) call check(error,.not. ieee_is_finite(led%rel_residual(ICE_LEDGER_WATER)), "residual is NaN")
      call finish_ledger_case(st, fx, area)
   end subroutine test_ledger_nonfinite

   ! ------------------------------------------------------------------
   subroutine test_throughput(error)
      !! Opt-in (`RDB_CICE_BENCH` set): cells/s of the three kernels on a
      !! 512 x 512 x 5 grid, printed to stderr.  Otherwise a no-op.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      character(len=8) :: e
      integer :: stat, ierr, r
      integer, parameter :: REPS = 20
      integer(int64) :: t0, t1, rate
      real(wp) :: t_zero, t_merge, t_save, cells

      call get_environment_variable("RDB_CICE_BENCH", e, status=stat)
      if (stat /= 0) return
      call make_slots(st, fx, ierr, nxp_in=508, nyp_in=508)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_state_random(st, 3_int64, 0.2_wp)
      FLW_BENCH = 250.0_wp
      !$acc enter data copyin(FLW_BENCH)
      call st%enter_data()
      call fx%enter_data()
      cells = real(st%nx, wp)*real(st%ny, wp)*real(REPS, wp)
      call ice_save_init(st)
      call ice_flux_zero(fx)
      call ice_flux_merge_categories(fx, st, FLW_BENCH, ierr)   ! warm up
      !$acc update self(fx%fresh)
      call system_clock(t0, rate)
      do r = 1, REPS
         call ice_flux_zero(fx)
      end do
      !$acc update self(fx%fresh)
      call system_clock(t1)
      t_zero = real(t1 - t0, wp)/real(rate, wp)
      call system_clock(t0)
      do r = 1, REPS
         call ice_flux_merge_categories(fx, st, FLW_BENCH, ierr)
      end do
      !$acc update self(fx%fresh)
      call system_clock(t1)
      t_merge = real(t1 - t0, wp)/real(rate, wp)
      call system_clock(t0)
      do r = 1, REPS
         call ice_save_init(st)
      end do
      !$acc update self(st%aicen_init)
      call system_clock(t1)
      t_save = real(t1 - t0, wp)/real(rate, wp)
      write (error_unit, "(a,i0,a,i0,a,i0)") "  [cice_flux] throughput on ", st%nx, " x ", st%ny, " x ncat ", NCAT
      write (error_unit, "(a,es10.3,a)") "  [cice_flux] ice_flux_zero             ", cells/t_zero, " cells/s"
      write (error_unit, "(a,es10.3,a)") "  [cice_flux] ice_flux_merge_categories ", cells/t_merge, " cells/s"
      write (error_unit, "(a,es10.3,a)") "  [cice_flux] ice_save_init             ", cells/t_save, " cells/s"
      !$acc exit data delete(FLW_BENCH)
      call fx%exit_data()
      call st%exit_data()
   end subroutine test_throughput

end module test_cice_flux
