!! Tests of the therm2 driver (`rdb_ice_therm2`, D4.3): rain -> aggregate -> linear remap -> new ice -> lateral melt
!! -> clean-up in one per-cell launch.
!!
!! * ORDER: the driver (`ice_step_therm2`, staged launches) equals the five stages called by hand in upstream's order (`ice_linear_itd_step`,
!!   `ice_new_ice_step`, `ice_lateral_melt_step`, `ice_cleanup_itd_step` + the rain term) bit for bit on CPU builds, on a random
!!   multi-cell state that holds growth, melt, remap fallbacks, new ice, lateral melt, rain, ice-free cells and zaps; and it
!!   DIFFERS from three wrong orders on a state built to tell them apart (lateral melt before the new ice, new ice before the
!!   remap, clean-up before the lateral melt);
!! * T12 (no double count): a supercooled ocean column's heat deficit equals the enthalpy of the new ice, `fhocn` untouched,
!!   the bank spent only where `frzmlt > 0` (the ocean-side accumulator lives in the ocean core, which an ice-only test does
!!   not link: the bank is built by hand from the column);
!! * T13: the virtual-salt adapter `1000*fsalt - S*fresh` equals `m*(S - S_ice)/dt`;
!! * L4: a cell holding only new ice is melted laterally as new ice (volume, enthalpy, fluxes in closed form);
!! * the `frzmlt` sign split (freezing-only, melting-only and neutral cells), a cell with a fatal bit untouched bitwise, NaN
!!   input flags only its own cell, a remap-fallback cell ends inside its bins;
!! * exact closure of mass, salt and energy over the whole step on 100 000 random cells, against the flux increments and the
!!   frazil input;
!! * the Icepack ORACLE: the real upstream `icepack_step_therm2`, 44 cells per constant set (literals in
!!   `test_cice_therm2_golden.inc`), bitwise on CPU builds and scale-aware on the offload build;
!! * a device sweep of the host shim over an INTERIOR range against the host cell kernel (every output of every cell, the ring
!!   untouched), and the throughput on a random and on a spatially coherent ice-free state.
module test_cice_therm2
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_REMAP_FALLBACK, ICE_ERR_NEG_THICK, ICE_ERR_AICE_RANGE, ICE_ERR_FATAL_MASK
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_ntrcr, ice_trcr_get, ice_trcr_set, ice_qsno_from_T, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_itd_shift, only: itd_aggregate
   use rdb_ice_itd_linear, only: ice_linear_itd_step
   use rdb_ice_new_ice, only: ice_new_ice_step
   use rdb_ice_lateral_melt, only: ice_lateral_melt_step
   use rdb_ice_itd_cleanup, only: ice_cleanup_itd_step
   use rdb_ice_therm2, only: ice_step_therm2, ice_therm2_remap_cell, ice_therm2_newice_cell, ice_therm2_lateral_cell, &
                             ice_therm2_cleanup_cell
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_therm2_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TF = -1.8_wp
   integer, parameter :: NCAT = 5
   integer, parameter :: NCO = 44
   real(wp), parameter :: HIN(0:NCAT) = [0.0_wp, 0.6445072168194258_wp, 1.3914334975763036_wp, &
                                         2.4701793819598885_wp, 4.567287918850491_wp, 999.9_wp]
   real(wp), parameter :: SCL(8) = [0.3_wp, 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 0.6_wp, 3.0_wp, 8.0_wp]
#ifdef RDB_GPU_OFFLOAD
   logical, parameter :: GPU_BUILD = .true.
   real(wp), parameter :: TOL_GOLD = 1.0e-12_wp
   real(wp), parameter :: TOL_EXT = 7.0e-12_wp
      !! Device vs host area and volumes against the cell's totals (measured worst 6.9e-13).
   real(wp), parameter :: TOL_INT = 5.0e-11_wp
   real(wp), parameter :: TOL_FLX = 3.0e-15_wp
      !! Device vs host flux increments, relative to their gross magnitude.
      !! Device vs host intensive tracers of healthy categories, relative (measured worst 4.7e-12).
   real(wp), parameter :: TOL_ORD = 1.0e-13_wp
      !! Fused vs staged (both on the device), relative.
#else
   logical, parameter :: GPU_BUILD = .false.
   real(wp), parameter :: TOL_GOLD = 0.0_wp   ! bitwise on CPU builds
   real(wp), parameter :: TOL_EXT = 0.0_wp
   real(wp), parameter :: TOL_INT = 0.0_wp
   real(wp), parameter :: TOL_FLX = 0.0_wp
   real(wp), parameter :: TOL_ORD = 0.0_wp
#endif
   real(wp), parameter :: CLOSE_TOL = 1.0e-12_wp
      !! Round-off closure of the conserved sums (relative to the sum of magnitudes).

   integer, save :: GNIL(NCO)
   real(wp), save :: GII(126, NCO), GIO(102, NCO), GCI(126, NCO), GCO(102, NCO)

   type :: st_t
      !! A rectangular patch of cells: the whole state, the snapshots, the cell fields and the fluxes.
      integer :: nx = 0, ny = 0, nilyr = 7
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), T(:, :, :)
      real(wp), allocatable :: q(:, :, :, :), si(:, :, :, :), qs(:, :, :, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :), sc(:, :, :), qsc(:, :, :)
      real(wp), allocatable :: fr(:, :), fs(:, :), fh(:, :), ml(:, :), fv(:, :)
      real(wp), allocatable :: Tf(:, :), frz(:, :), frain(:, :), rside(:, :), bank(:, :)
      logical, allocatable :: wet(:, :)
      integer(int32), allocatable :: ierr(:, :)
   end type st_t

contains

   subroutine collect_cice_therm2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fused_equals_staged_upstream_order", test_order_equal), &
                  new_unittest("wrong_orders_differ", test_order_differs), &
                  new_unittest("t12_ice_side", test_t12), &
                  new_unittest("t13_virtual_salt_adapter", test_t13), &
                  new_unittest("l4_lateral_melt_after_new_ice", test_l4), &
                  new_unittest("frzmlt_sign_split_cells", test_sign_split), &
                  new_unittest("rside_zero_is_skipping_the_stage", test_rside_zero), &
                  new_unittest("aggregates_are_consistent_after_step", test_aggregates), &
                  new_unittest("fatal_bit_cell_untouched", test_fatal_untouched), &
                  new_unittest("midway_fatal_skips_cleanup", test_midway_fatal), &
                  new_unittest("nan_input_flags_only_its_cell", test_nan), &
                  new_unittest("remap_fallback_ends_inside_bins", test_fallback_bins), &
                  new_unittest("conservation_100k_cells", test_conserve), &
                  new_unittest("oracle_therm2_icepack", test_oracle_icepack), &
                  new_unittest("oracle_therm2_cesm", test_oracle_cesm), &
                  new_unittest("device_sweep_and_throughput", test_device) &
                  ]
   end subroutine collect_cice_therm2_tests

   subroutine fill_golden()
#include "test_cice_therm2_golden.inc"
   end subroutine fill_golden

   ! ------------------------------------------------------------------
   ! state container, device mapping, comparison
   ! ------------------------------------------------------------------

   subroutine st_alloc(st, nx, ny, nilyr)
      type(st_t), intent(out) :: st
      integer, intent(in) :: nx, ny, nilyr

      st%nx = nx; st%ny = ny; st%nilyr = nilyr
      allocate (st%ai0(nx, ny, NCAT), st%vi0(nx, ny, NCAT), st%a(nx, ny, NCAT), st%v(nx, ny, NCAT), st%s(nx, ny, NCAT), &
                st%T(nx, ny, NCAT), st%q(nx, ny, NCAT, nilyr), st%si(nx, ny, NCAT, nilyr), st%qs(nx, ny, NCAT, 1), &
                st%aice(nx, ny), st%vice(nx, ny), st%vsno(nx, ny), st%aice0(nx, ny), st%Tc(nx, ny), st%qc(nx, ny, nilyr), &
                st%sc(nx, ny, nilyr), st%qsc(nx, ny, 1), st%fr(nx, ny), st%fs(nx, ny), st%fh(nx, ny), st%ml(nx, ny), &
                st%fv(nx, ny), st%Tf(nx, ny), st%frz(nx, ny), st%frain(nx, ny), st%rside(nx, ny), st%bank(nx, ny), &
                st%wet(nx, ny), st%ierr(nx, ny))
      st%ai0 = 0.0_wp; st%vi0 = 0.0_wp; st%a = 0.0_wp; st%v = 0.0_wp; st%s = 0.0_wp; st%T = TF
      st%q = 0.0_wp; st%si = 0.0_wp; st%qs = 0.0_wp
      st%aice = 0.0_wp; st%vice = -777.0_wp; st%vsno = -777.0_wp; st%aice0 = -777.0_wp; st%Tc = -777.0_wp
      st%qc = -777.0_wp; st%sc = -777.0_wp; st%qsc = -777.0_wp
      st%fr = 0.0_wp; st%fs = 0.0_wp; st%fh = 0.0_wp; st%ml = 0.0_wp; st%fv = 0.0_wp
      st%Tf = TF; st%frz = 0.0_wp; st%frain = 0.0_wp; st%rside = 0.0_wp; st%bank = 0.0_wp
      st%wet = .true.
      st%ierr = 0_int32
   end subroutine st_alloc

   subroutine st_enter(st)
      type(st_t), intent(inout) :: st
      !$acc enter data copyin(st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, st%ml, st%fv, st%Tf, st%frz, st%frain, st%rside, st%bank, st%wet, st%ierr)
   end subroutine st_enter

   subroutine st_exit(st)
      type(st_t), intent(inout) :: st
      !$acc update self(st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, st%ml, st%fv, st%Tf, st%frz, st%frain, st%rside, st%bank, st%wet, st%ierr)
      !$acc exit data delete(st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, st%ml, st%fv, st%Tf, st%frz, st%frain, st%rside, st%bank, st%wet, st%ierr)
   end subroutine st_exit

   subroutine st_set_cell(st, i, j, ai, vi, a, v, s, tr)
      !! Load cell `(i, j)` from the packed layout `tr(ICE_NTRCR_MAX, NCAT)` (Tsfc | qice | qsno | sice).
      type(st_t), intent(inout) :: st
      integer, intent(in) :: i, j
      real(wp), intent(in) :: ai(NCAT), vi(NCAT), a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, it

      do n = 1, NCAT
         st%ai0(i, j, n) = ai(n); st%vi0(i, j, n) = vi(n)
         st%a(i, j, n) = a(n); st%v(i, j, n) = v(n); st%s(i, j, n) = s(n)
         do it = 1, ice_ntrcr(st%nilyr, 1)
            call ice_trcr_set(it, i, j, n, st%nx, st%ny, NCAT, st%nilyr, 1, tr(it, n), st%T, st%q, st%si, st%qs)
         end do
      end do
   end subroutine st_set_cell

   subroutine st_get_cell(st, i, j, a, v, s, tr)
      type(st_t), intent(in) :: st
      integer, intent(in) :: i, j
      real(wp), intent(out) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, it

      tr = 0.0_wp
      do n = 1, NCAT
         a(n) = st%a(i, j, n); v(n) = st%v(i, j, n); s(n) = st%s(i, j, n)
         do it = 1, ice_ntrcr(st%nilyr, 1)
            tr(it, n) = ice_trcr_get(it, i, j, n, st%nx, st%ny, NCAT, st%nilyr, 1, st%T, st%q, st%si, st%qs)
         end do
      end do
   end subroutine st_get_cell

   subroutine rain_stage(nx, ny, i_lo, i_hi, j_lo, j_hi, wet, frain, aice, fresh)
      !! Stage 1 by hand: `fresh += frain*aice`.
      integer, intent(in) :: nx, ny, i_lo, i_hi, j_lo, j_hi
      logical, intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: frain(nx, ny), aice(nx, ny)
      real(wp), intent(inout) :: fresh(nx, ny)
      integer :: i, j

      !$acc data present(wet, frain, aice, fresh)
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         if (wet(i, j)) fresh(i, j) = fresh(i, j) + frain(i, j)*aice(i, j)
      end do
      !$acc end data
   end subroutine rain_stage

   subroutine run_fused(st, p, dt, i_lo, i_hi, j_lo, j_hi, sal_in, fuse_bc)
      !! The driver: map, one launch, bring back.  `sal_in`: the new-ice salinity profile (default: the BL99 one).
      type(st_t), intent(inout) :: st
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in), optional :: sal_in(:)
      logical, intent(in), optional :: fuse_bc
      real(wp) :: hin_l(0:NCAT), sal(st%nilyr)

      hin_l = HIN
      sal = salinz_of(p, st%nilyr)
      if (present(sal_in)) sal = sal_in(1:st%nilyr)
      call st_enter(st)
      !$acc enter data copyin(hin_l, sal)
      call ice_step_therm2(st%nx, st%ny, NCAT, st%nilyr, 1, i_lo, i_hi, j_lo, j_hi, dt, p, hin_l, sal, st%wet, st%Tf, st%frz, &
                           st%frain, st%rside, st%bank, st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, &
                           st%aice, st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, st%ml, &
                           st%fv, st%ierr, fuse_bc)
      !$acc exit data delete(hin_l, sal)
      call st_exit(st)
   end subroutine run_fused

   subroutine run_staged(st, p, dt, i_lo, i_hi, j_lo, j_hi, order)
      !! The stages by hand through the owning modules' shims, in the given order (1 rain, 2 remap, 3 new ice,
      !! 4 lateral melt, 5 clean-up).  Upstream's order is `[1, 2, 3, 4, 5]`.  Maps and unmaps the state.
      type(st_t), intent(inout) :: st
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi, order(:)
      real(wp) :: hin_l(0:NCAT), sal(st%nilyr)

      hin_l = HIN
      sal = salinz_of(p, st%nilyr)
      call st_enter(st)
      !$acc enter data copyin(hin_l, sal)
      call staged_stages(st, p, dt, i_lo, i_hi, j_lo, j_hi, order, hin_l, sal)
      !$acc exit data delete(hin_l, sal)
      call st_exit(st)
   end subroutine run_staged

   subroutine staged_stages(st, p, dt, i_lo, i_hi, j_lo, j_hi, order, hin_l, sal)
      !! The stage launches of `run_staged` on an already mapped state.
      type(st_t), intent(inout) :: st
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi, order(:)
      real(wp), intent(in) :: hin_l(0:NCAT), sal(st%nilyr)
      integer :: k

      do k = 1, size(order)
         select case (order(k))
         case (1)
            call rain_stage(st%nx, st%ny, i_lo, i_hi, j_lo, j_hi, st%wet, st%frain, st%aice, st%fr)
         case (2)
            call ice_linear_itd_step(st%nx, st%ny, NCAT, st%nilyr, 1, i_lo, i_hi, j_lo, j_hi, p, hin_l, st%Tf, st%ai0, st%vi0, &
                                     st%a, st%v, st%s, st%T, st%q, st%si, st%qs, .false., st%aice, st%aice0, st%ierr)
         case (3)
            call ice_new_ice_step(st%nx, st%ny, NCAT, st%nilyr, i_lo, i_hi, j_lo, j_hi, dt, p, hin_l, sal, st%wet, st%frz, &
                                  st%bank, st%Tf, st%a, st%v, st%T, st%q, st%si, st%aice0, st%fr, st%fs, st%fv, st%ierr)
         case (4)
            call ice_lateral_melt_step(st%nx, st%ny, NCAT, st%nilyr, 1, i_lo, i_hi, j_lo, j_hi, dt, p, st%rside, st%a, st%v, &
                                       st%s, st%q, st%qs, st%fr, st%fs, st%fh, st%ml, st%ierr)
         case (5)
            call ice_cleanup_itd_step(st%nx, st%ny, NCAT, st%nilyr, 1, i_lo, i_hi, j_lo, j_hi, p, hin_l, dt, st%Tf, &
                                      .true., .true., .true., st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, st%vice, &
                                      st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, st%ierr)
         end select
      end do
   end subroutine staged_stages

   subroutine run_host_cells(st, p, dt, i_lo, i_hi, j_lo, j_hi)
      !! The host reference: the cell kernel called from a serial host loop on the (unmapped) host arrays.
      type(st_t), intent(inout) :: st
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi
      real(wp) :: sal(st%nilyr)
      integer :: i, j

      sal = salinz_of(p, st%nilyr)
      do j = j_lo, j_hi
         do i = i_lo, i_hi
            if (.not. st%wet(i, j)) cycle
            call ice_therm2_remap_cell(i, j, st%nx, st%ny, NCAT, st%nilyr, 1, p%itd%kitd == 1, p%c%rhos, p%c%Lfresh, &
                                       p%thermo%hi_min, p%c%puny, HIN, st%Tf(i, j), st%frz(i, j), st%frain(i, j), &
                                       st%rside(i, j), st%aice, st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, &
                                       st%qs, st%fr, st%ierr)
            call ice_therm2_newice_cell(i, j, st%nx, st%ny, NCAT, st%nilyr, dt, p%c%rhoi, p%c%Lfresh, p%c%puny, &
                                        p%thermo%ice_ref_salinity, p%thermo%hfrazilmin, p%thermo%update_ocn_f, HIN, sal, &
                                        st%Tf(i, j), st%frz(i, j), st%bank, st%a, st%v, st%T, st%q, st%si, st%aice0, &
                                        st%fr, st%fs, st%fv, st%ierr)
            call ice_therm2_lateral_cell(i, j, st%nx, st%ny, NCAT, st%nilyr, 1, dt, p%c%rhoi, p%c%rhos, &
                                         p%thermo%ice_ref_salinity, p%c%puny, st%rside(i, j), st%a, st%v, st%s, st%q, &
                                         st%qs, st%fr, st%fs, st%fh, st%ml, st%ierr)
            call ice_therm2_cleanup_cell(i, j, st%nx, st%ny, NCAT, st%nilyr, 1, HIN, dt, st%Tf(i, j), p%c%rhoi, p%c%rhos, &
                                         p%c%Lfresh, p%c%cp_ice, p%c%hs_min, p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, &
                                         p%thermo%ice_ref_salinity, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, &
                                         st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, &
                                         st%ierr)
         end do
      end do
   end subroutine run_host_cells

   pure subroutine note(worst, x, nbad)
      !! NaN-proof running maximum: a non-finite `x` is counted in `nbad`, never hidden by `max`.
      real(wp), intent(inout) :: worst
      real(wp), intent(in) :: x
      integer, intent(inout) :: nbad

      if (.not. ieee_is_finite(x)) then
         nbad = nbad + 1
      else if (x > worst) then
         worst = x
      end if
   end subroutine note

   function salinz_of(p, nilyr) result(sal)
      !! The BL99 salinity profile of new ice (`saltmax*(1 - cos(pi*z**(nsal/(msal + z))))/2`, midpoint of each layer).
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: nilyr
      real(wp) :: sal(nilyr)
      real(wp) :: zn, ex
      integer :: k

      do k = 1, nilyr
         zn = (real(k, wp) - 0.5_wp)/real(nilyr, wp)
         ex = 0.407_wp/(0.573_wp + zn)
         sal(k) = 0.5_wp*p%c%saltmax*(1.0_wp - cos(3.14159265358979323846_wp*zn**ex))
      end do
   end function salinz_of

   pure logical function same(a, b, tol)
      !! Bitwise (`tol = 0`; NaN matches NaN) or relative.
      real(wp), intent(in) :: a, b, tol

      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same = ieee_is_nan(a) .and. ieee_is_nan(b)
      else if (tol == 0.0_wp) then
         same = (a == b)
      else
         same = a == b .or. abs(a - b) <= tol*max(abs(a), abs(b))
      end if
   end function same

   subroutine cmp_flat(n, x, y, tol, nbad, worst)
      integer, intent(in) :: n
      real(wp), intent(in) :: x(n), y(n), tol
      integer, intent(inout) :: nbad
      real(wp), intent(inout) :: worst
      integer :: i
      real(wp) :: d

      do i = 1, n
         if (.not. same(x(i), y(i), tol)) nbad = nbad + 1
         if (ieee_is_finite(x(i)) .and. ieee_is_finite(y(i)) .and. x(i) /= y(i)) then
            d = abs(x(i) - y(i))/max(abs(x(i)), abs(y(i)))
            if (d > worst) worst = d
         end if
      end do
   end subroutine cmp_flat

   subroutine st_compare(x, y, tol, nbad, worst)
      !! Every array of two states (all cells), `nbad` elements outside `tol`, `worst` relative difference of the finite pairs.
      type(st_t), intent(in) :: x, y
      real(wp), intent(in) :: tol
      integer, intent(out) :: nbad
      real(wp), intent(out) :: worst

      nbad = 0; worst = 0.0_wp
      call cmp_flat(size(x%a), x%a, y%a, tol, nbad, worst); call cmp_flat(size(x%v), x%v, y%v, tol, nbad, worst)
      call cmp_flat(size(x%s), x%s, y%s, tol, nbad, worst); call cmp_flat(size(x%T), x%T, y%T, tol, nbad, worst)
      call cmp_flat(size(x%q), x%q, y%q, tol, nbad, worst); call cmp_flat(size(x%si), x%si, y%si, tol, nbad, worst)
      call cmp_flat(size(x%qs), x%qs, y%qs, tol, nbad, worst)
      call cmp_flat(size(x%aice), x%aice, y%aice, tol, nbad, worst); call cmp_flat(size(x%vice), x%vice, y%vice, tol, nbad, worst)
      call cmp_flat(size(x%vsno), x%vsno, y%vsno, tol, nbad, worst)
      call cmp_flat(size(x%aice0), x%aice0, y%aice0, tol, nbad, worst); call cmp_flat(size(x%Tc), x%Tc, y%Tc, tol, nbad, worst)
      call cmp_flat(size(x%qc), x%qc, y%qc, tol, nbad, worst); call cmp_flat(size(x%sc), x%sc, y%sc, tol, nbad, worst)
      call cmp_flat(size(x%qsc), x%qsc, y%qsc, tol, nbad, worst)
      call cmp_flat(size(x%fr), x%fr, y%fr, tol, nbad, worst); call cmp_flat(size(x%fs), x%fs, y%fs, tol, nbad, worst)
      call cmp_flat(size(x%fh), x%fh, y%fh, tol, nbad, worst); call cmp_flat(size(x%ml), x%ml, y%ml, tol, nbad, worst)
      call cmp_flat(size(x%fv), x%fv, y%fv, tol, nbad, worst); call cmp_flat(size(x%bank), x%bank, y%bank, tol, nbad, worst)
      nbad = nbad + count(x%ierr /= y%ierr)
   end subroutine st_compare

   function cell_identical(x, y, i, j) result(ok)
      !! Every array of cell `(i, j)` of two states is bitwise equal.
      type(st_t), intent(in) :: x, y
      integer, intent(in) :: i, j
      logical :: ok

      ok = all(same_v(x%a(i, j, :), y%a(i, j, :))) .and. all(same_v(x%v(i, j, :), y%v(i, j, :))) .and. &
           all(same_v(x%s(i, j, :), y%s(i, j, :))) .and. all(same_v(x%T(i, j, :), y%T(i, j, :))) .and. &
           all(same_v(reshape(x%q(i, j, :, :), [NCAT*x%nilyr]), reshape(y%q(i, j, :, :), [NCAT*x%nilyr]))) .and. &
           all(same_v(reshape(x%si(i, j, :, :), [NCAT*x%nilyr]), reshape(y%si(i, j, :, :), [NCAT*x%nilyr]))) .and. &
           all(same_v(x%qs(i, j, :, 1), y%qs(i, j, :, 1))) .and. &
           all(same_v([x%aice(i, j), x%vice(i, j), x%vsno(i, j), x%aice0(i, j), x%Tc(i, j), x%fr(i, j), x%fs(i, j), &
                       x%fh(i, j), x%ml(i, j), x%fv(i, j), x%bank(i, j)], &
                      [y%aice(i, j), y%vice(i, j), y%vsno(i, j), y%aice0(i, j), y%Tc(i, j), y%fr(i, j), y%fs(i, j), &
                       y%fh(i, j), y%ml(i, j), y%fv(i, j), y%bank(i, j)])) .and. &
           all(same_v(x%qc(i, j, :), y%qc(i, j, :))) .and. all(same_v(x%sc(i, j, :), y%sc(i, j, :))) .and. &
           x%ierr(i, j) == y%ierr(i, j)
   end function cell_identical

   pure function same_v(x, y) result(r)
      real(wp), intent(in) :: x(:), y(:)
      logical :: r(size(x))
      integer :: k

      do k = 1, size(x)
         r(k) = same(x(k), y(k), 0.0_wp)
      end do
   end function same_v

   ! ------------------------------------------------------------------
   ! random states
   ! ------------------------------------------------------------------

   subroutine fill_state(st, seed, mode)
      !! `mode = 1`: wild random cells (growth and melt up to several metres per step, out-of-bin thicknesses, snow some
      !! warm, ice-free cells, new ice with and without a bank, rain, lateral melt, a few dry cells, a few cells with a stale
      !! bank and `frzmlt <= 0`).  `mode = 2`: the realistic case, ~85 % ice free in 32 x 32 blocks, ice cells clean and
      !! nearly stationary, 5 % of the ice-free cells freezing.
      type(st_t), intent(inout) :: st
      integer(int64), intent(in) :: seed
      integer, intent(in) :: mode
      integer(int64) :: stt
      integer :: i, j, n, k
      real(wp) :: u1, u2, u3, u4, lo, hi, h0, dh, sc, tot
      logical :: ice

      do j = 1, st%ny
         do i = 1, st%nx
            stt = lcg_next(seed + int(i, int64)*7919_int64 + int(j, int64)*104729_int64*65537_int64)
            stt = lcg_next(stt); u1 = lcg_unit(stt)
            sc = SCL(1 + int(u1*8.0_wp))
            ice = .true.
            if (mode == 1) then
               stt = lcg_next(stt); ice = lcg_unit(stt) > 0.06_wp                       ! 6 % ice-free cells
            end if
            if (mode == 2) then
               stt = lcg_next(seed + int((i - 1)/32, int64)*7919_int64 + int((j - 1)/32, int64)*104729_int64*65537_int64)
               stt = lcg_next(stt); ice = lcg_unit(stt) < 0.15_wp
               stt = lcg_next(seed + int(i, int64)*31_int64 + int(j, int64)*104729_int64*131_int64)
               sc = 0.05_wp
            end if
            tot = 0.0_wp
            do n = 1, NCAT
               stt = lcg_next(stt); u1 = lcg_unit(stt)
               stt = lcg_next(stt); u2 = lcg_unit(stt)
               stt = lcg_next(stt); u3 = lcg_unit(stt)
               stt = lcg_next(stt); u4 = lcg_unit(stt)
               st%ai0(i, j, n) = 0.0_wp; st%vi0(i, j, n) = 0.0_wp; st%a(i, j, n) = 0.0_wp; st%v(i, j, n) = 0.0_wp
               st%s(i, j, n) = 0.0_wp
               st%T(i, j, n) = TF; st%q(i, j, n, :) = 0.0_wp; st%si(i, j, n, :) = 0.0_wp; st%qs(i, j, n, 1) = 0.0_wp
               if (ice .and. u1 < merge(0.75_wp, 0.8_wp, mode == 1)) then
                  st%ai0(i, j, n) = 0.02_wp + 0.17_wp*u2
                  lo = HIN(n - 1); hi = min(HIN(n), 8.0_wp)
                  h0 = lo + (hi - lo)*(0.05_wp + 0.9_wp*u3)
                  if (mode == 1 .and. u4 > 0.97_wp) h0 = 0.05_wp + 6.0_wp*u3      ! out of its bin
                  dh = sc*0.2_wp*(u4 - 0.4_wp)
                  dh = max(dh, -0.9_wp*h0)
                  if (h0 + dh < 0.002_wp) dh = 0.002_wp - h0
                  st%vi0(i, j, n) = st%ai0(i, j, n)*h0
                  st%a(i, j, n) = st%ai0(i, j, n)
                  st%v(i, j, n) = st%a(i, j, n)*(h0 + dh)
                  tot = tot + st%a(i, j, n)
                  stt = lcg_next(stt)
                  if (lcg_unit(stt) > 0.3_wp) then
                     stt = lcg_next(stt); st%s(i, j, n) = st%a(i, j, n)*0.4_wp*lcg_unit(stt)
                  end if
                  stt = lcg_next(stt); st%T(i, j, n) = -30.0_wp*lcg_unit(stt)
                  do k = 1, st%nilyr
                     stt = lcg_next(stt); st%q(i, j, n, k) = -4.0e8_wp + 3.0e8_wp*lcg_unit(stt)
                     stt = lcg_next(stt); st%si(i, j, n, k) = 0.1_wp + 5.9_wp*lcg_unit(stt)
                  end do
                  stt = lcg_next(stt)
                  if (mode == 1) then
                     st%qs(i, j, n, 1) = -2.0e8_wp + 0.95e8_wp*lcg_unit(stt)         ! some snow above the melting point
                  else
                     st%qs(i, j, n, 1) = -1.45e8_wp + 0.25e8_wp*lcg_unit(stt)
                  end if
               end if
            end do
            if (tot > 0.97_wp) then
               st%ai0(i, j, :) = st%ai0(i, j, :)*(0.97_wp/tot); st%a(i, j, :) = st%a(i, j, :)*(0.97_wp/tot)
               st%vi0(i, j, :) = st%vi0(i, j, :)*(0.97_wp/tot); st%v(i, j, :) = st%v(i, j, :)*(0.97_wp/tot)
               st%s(i, j, :) = st%s(i, j, :)*(0.97_wp/tot)
            end if
            st%aice(i, j) = sum(st%ai0(i, j, :))
            ! cell-level inputs
            st%Tf(i, j) = TF
            st%frz(i, j) = 0.0_wp; st%bank(i, j) = 0.0_wp; st%rside(i, j) = 0.0_wp; st%frain(i, j) = 0.0_wp
            st%wet(i, j) = .true.
            stt = lcg_next(stt); u1 = lcg_unit(stt)
            stt = lcg_next(stt); u2 = lcg_unit(stt)
            if (mode == 1) then
               if (u1 < 0.25_wp) then
                  st%bank(i, j) = 10.0_wp**(3.0_wp + 4.0_wp*u2); st%frz(i, j) = st%bank(i, j)/2400.0_wp
               else if (u1 < 0.40_wp) then
                  st%frz(i, j) = -10.0_wp - 90.0_wp*u2; st%bank(i, j) = 1.0e3_wp + 1.0e5_wp*u2    ! stale bank: retained
               end if
               stt = lcg_next(stt); u3 = lcg_unit(stt)
               stt = lcg_next(stt); u4 = lcg_unit(stt)
               if (u3 < 0.3_wp) st%rside(i, j) = 0.2_wp*u4
               if (u3 > 0.99_wp) st%rside(i, j) = 1.0_wp
               stt = lcg_next(stt); u3 = lcg_unit(stt)
               stt = lcg_next(stt); u4 = lcg_unit(stt)
               if (u3 < 0.4_wp) st%frain(i, j) = 5.0e-5_wp*u4
               stt = lcg_next(stt); u3 = lcg_unit(stt)
               st%Tf(i, j) = -1.6_wp - 0.3_wp*u3
               stt = lcg_next(stt); u3 = lcg_unit(stt)
               if (u3 < 0.03_wp) st%wet(i, j) = .false.
            else
               if (.not. ice .and. u1 < 0.05_wp) then
                  st%bank(i, j) = 10.0_wp**(3.0_wp + 3.0_wp*u2); st%frz(i, j) = st%bank(i, j)/2400.0_wp
               end if
            end if
         end do
      end do
   end subroutine fill_state

   ! ------------------------------------------------------------------
   ! helpers of the closed-form tests
   ! ------------------------------------------------------------------

   function params_of(nilyr, updf) result(p)
      integer, intent(in) :: nilyr
      logical, intent(in) :: updf
      type(ice_params_t) :: p

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, nilyr, 1)
      p%thermo%update_ocn_f = updf
   end function params_of

   subroutine set_simple_cell(st, p, i, j, a, h, sd, Tsn)
      !! A clean, stationary column: areas `a`, thicknesses `h`, snow depths `sd`, snow of temperature `Tsn`, Tsfc -5,
      !! qice(k) = -(2e8 + 1e7 k), sice(k) = 3 + 0.2 k; the snapshots equal the state.
      type(st_t), intent(inout) :: st
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: i, j
      real(wp), intent(in) :: a(NCAT), h(NCAT), sd, Tsn
      integer :: n, k

      do n = 1, NCAT
         st%ai0(i, j, n) = a(n); st%a(i, j, n) = a(n)
         st%vi0(i, j, n) = a(n)*h(n); st%v(i, j, n) = a(n)*h(n)
         st%s(i, j, n) = a(n)*sd
         st%T(i, j, n) = -5.0_wp
         do k = 1, st%nilyr
            st%q(i, j, n, k) = -(2.0e8_wp + 1.0e7_wp*real(k, wp))
            st%si(i, j, n, k) = 3.0_wp + 0.2_wp*real(k, wp)
         end do
         st%qs(i, j, n, 1) = ice_qsno_from_T(Tsn, p%c)
      end do
      st%aice(i, j) = sum(a)
   end subroutine set_simple_cell

   function cell_energy(st, i, j) result(e)
      !! `sum_n sum_k qice*vice/nilyr + qsno*vsno` of cell `(i, j)` [J/m2].
      type(st_t), intent(in) :: st
      integer, intent(in) :: i, j
      real(wp) :: e
      integer :: n, k

      e = 0.0_wp
      do n = 1, NCAT
         do k = 1, st%nilyr
            e = e + st%q(i, j, n, k)*st%v(i, j, n)/real(st%nilyr, wp)
         end do
         e = e + st%qs(i, j, n, 1)*st%s(i, j, n)
      end do
   end function cell_energy

   ! ------------------------------------------------------------------
   ! order
   ! ------------------------------------------------------------------

   subroutine test_order_equal(error)
      !! The fused kernel against the stages called by hand in upstream's order: random wild cells, every output.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: s1, s2, s3
      type(ice_params_t) :: p
      integer :: nbad, nbad3, i, j, nfb, nnew, nlat, nempty, nrain
      real(wp) :: worst, worst3

      p = params_of(7, .true.)
      call st_alloc(s1, 64, 48, 7)
      call fill_state(s1, 424242_int64, 1)
      s1%wet = .true.
      s2 = s1
      s3 = s1
      call run_fused(s1, p, 2400.0_wp, 1, 64, 1, 48)
      call run_fused(s3, p, 2400.0_wp, 1, 64, 1, 48, fuse_bc=.false.)
      call run_staged(s2, p, 2400.0_wp, 1, 64, 1, 48, [1, 2, 3, 4, 5])
      call st_compare(s1, s2, TOL_ORD, nbad, worst)
      call st_compare(s3, s2, TOL_ORD, nbad3, worst3)
      nbad = nbad + nbad3
      nfb = 0; nnew = 0; nlat = 0; nempty = 0; nrain = 0
      do j = 1, 48
         do i = 1, 64
            if (iand(s1%ierr(i, j), ICE_ERR_REMAP_FALLBACK) /= 0_int32) nfb = nfb + 1
            if (s1%fv(i, j) > 0.0_wp) nnew = nnew + 1
            if (s1%ml(i, j) > 0.0_wp) nlat = nlat + 1
            if (s1%frain(i, j) > 0.0_wp) nrain = nrain + 1
            if (all(s1%ai0(i, j, :) == 0.0_wp)) nempty = nempty + 1
         end do
      end do
    write (output_unit, '(a,i0,a,es10.2,a,5i6)') "order test: elements differing from the staged run ", nbad, ", worst relative ", &
         worst, "; cells with fallback, new ice, lateral melt, rain, ice-free: ", nfb, nnew, nlat, nrain, nempty
      call check(error, nbad, 0, "the driver equals the five stages called in upstream's order")
      if (allocated(error)) return
      call check(error, nfb > 20 .and. nnew > 200 .and. nlat > 200 .and. nrain > 500 .and. nempty > 20, &
                 "the state exercises every stage")
   end subroutine test_order_equal

   subroutine test_order_differs(error)
      !! A cell built to tell the orders apart: category 1 grows across its bound, there is open water, a freezing potential
      !! and a lateral-melt fraction.  Upstream's order reproduces the driver; four wrong orders do not (the last pins that the rain term uses the ENTRY aggregate).
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: base, ref, w
      type(ice_params_t) :: p
      integer :: nbad, k
      real(wp) :: worst, dt, vtarget
      integer, parameter :: ORD(5, 5) = reshape([1, 2, 3, 4, 5, &       ! upstream
                                                 1, 2, 4, 3, 5, &       ! lateral melt before the new ice
                                                 1, 3, 2, 4, 5, &       ! new ice before the remap
                                                 1, 2, 3, 5, 4, &       ! clean-up before the lateral melt
                                                 2, 3, 4, 5, 1], [5, 5])    ! rain after everything (stale-aggregate bug)
      character(len=24), parameter :: NAME(5) = ["upstream order        ", "lateral before new ice", &
                                                 "new ice before remap  ", "cleanup before lateral", &
                                                 "rain last             "]

      dt = 3600.0_wp
      p = params_of(7, .true.)
      call st_alloc(base, 3, 1, 7)
      vtarget = 0.1_wp
      call set_simple_cell(base, p, 1, 1, [0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], [0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      base%v(1, 1, 1) = 0.3_wp*0.9_wp                                   ! category 1: 0.3 m -> 0.9 m, across its bound
      base%frz(1, 1) = vtarget*p%c%rhoi*p%c%Lfresh/dt; base%bank(1, 1) = base%frz(1, 1)*dt
      base%rside(1, 1) = 0.1_wp; base%frain(1, 1) = 1.0e-5_wp
      ! cells 2 and 3: the same with growth in two categories and a thick bank
      call set_simple_cell(base, p, 2, 1, [0.2_wp, 0.3_wp, 0.1_wp, 0.0_wp, 0.0_wp], [0.4_wp, 1.0_wp, 2.0_wp, 0.0_wp, 0.0_wp], &
                           0.1_wp, -10.0_wp)
      base%v(2, 1, 1) = 0.2_wp*0.8_wp; base%v(2, 1, 2) = 0.3_wp*1.5_wp
      base%frz(2, 1) = 0.2_wp*p%c%rhoi*p%c%Lfresh/dt; base%bank(2, 1) = base%frz(2, 1)*dt; base%rside(2, 1) = 0.05_wp
      call set_simple_cell(base, p, 3, 1, [0.1_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp], [0.2_wp, 0.9_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      base%v(3, 1, 2) = 0.1_wp*1.4_wp; base%frz(3, 1) = 0.03_wp*p%c%rhoi*p%c%Lfresh/dt; base%rside(3, 1) = 0.2_wp
      ref = base
      call run_fused(ref, p, dt, 1, 3, 1, 1)
      do k = 1, 5
         w = base
         call run_staged(w, p, dt, 1, 3, 1, 1, ORD(:, k))
         call st_compare(ref, w, 0.0_wp, nbad, worst)
         write (output_unit, '(a,a,a,i0,a,es10.2)') "order ", NAME(k), ": elements differing from the driver ", nbad, &
            ", worst relative ", worst
         if (k == 1) then
            call st_compare(ref, w, TOL_ORD, nbad, worst)
            call check(error, nbad, 0, "upstream's order equals the driver")
         else
            call check(error, nbad >= 1 .and. worst > 1.0e-4_wp, "a wrong order is distinguishable: "//trim(NAME(k)))
         end if
         if (allocated(error)) return
      end do
   end subroutine test_order_differs

   ! ------------------------------------------------------------------
   ! T12, T13, L4
   ! ------------------------------------------------------------------

   subroutine test_t12(error)
      !! T12, ICE SIDE ONLY.  A bank is formed by hand from a supercooled ocean column (8 layers, heat deficit
      !! `rho0*cp*sum(h*dT)`); the ice turns it into new ice.  What this proves: the enthalpy of the new ice equals the
      !! deficit (`deficit + dE_ice = 0`), `fhocn` is bitwise untouched by the new ice, the bank is reset exactly where
      !! `frzmlt > 0` and kept where `frzmlt < 0` or the cell is dry.  What it does NOT prove: the ocean side (accumulator ->
      !! freezing potential -> ice) lives in the ocean core, which an ice-only test does not link; the real no-double-count
      !! test belongs to the coupling PRs.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, st0
      type(ice_params_t) :: p
      real(wp) :: dt, deficit, rho0, cpo, dTk, e0(4), e1(4), fh0, v_ice, err_e, err_f, err_m
      integer :: k, c, nbadf
      real(wp), parameter :: HK = 2.0_wp

      dt = 3600.0_wp
      p = params_of(7, .true.)
      rho0 = 1026.0_wp
      cpo = p%c%cp_ocn
      deficit = 0.0_wp
      do k = 1, 8
         dTk = 0.05_wp*real(k, wp)                    ! K below the freezing point
         deficit = deficit + rho0*cpo*HK*dTk          ! J/m2: the heat the clamp took out of the ice's future
      end do
      call st_alloc(st, 4, 1, 7)
      call set_simple_cell(st, p, 1, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      call set_simple_cell(st, p, 2, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      call set_simple_cell(st, p, 3, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      st%fh = 7.25_wp
      ! cell 1 and 4: freezing (4 over open water only); cell 2: a melt potential (frzmlt < 0) with a stale bank; cell 3 dry
      st%bank(1, 1) = deficit; st%frz(1, 1) = deficit/dt
      st%bank(4, 1) = deficit; st%frz(4, 1) = deficit/dt
      st%bank(2, 1) = 5.0e4_wp; st%frz(2, 1) = -30.0_wp
      st%bank(3, 1) = deficit; st%frz(3, 1) = deficit/dt; st%wet(3, 1) = .false.
      st0 = st
      do c = 1, 4
         e0(c) = cell_energy(st0, c, 1)
      end do
      call run_fused(st, p, dt, 1, 4, 1, 1)
      do c = 1, 4
         e1(c) = cell_energy(st, c, 1)
      end do
      fh0 = 7.25_wp
      err_e = 0.0_wp; err_f = 0.0_wp; err_m = 0.0_wp; nbadf = 0
      do c = 1, 4, 3
         call note(err_e, abs(deficit + (e1(c) - e0(c)))/deficit, nbadf)
         call note(err_f, abs(st%fv(c, 1)*p%c%rhoi*p%c%Lfresh - deficit)/deficit, nbadf)
         v_ice = sum(st%v(c, 1, :)) - sum(st0%v(c, 1, :))
         call note(err_m, abs(p%c%rhoi*v_ice + dt*(st%fr(c, 1) - 0.0_wp))/(p%c%rhoi*v_ice), nbadf)
         call check(error, st%fh(c, 1) == fh0, "new ice never writes fhocn (bitwise)")
         if (allocated(error)) return
         call check(error, st%bank(c, 1) == 0.0_wp, "the bank is spent where frzmlt > 0")
         if (allocated(error)) return
      end do
      write (output_unit, '(a,3es10.2)') "T12: worst |deficit + dE_ice|/deficit, |frazil_vol*rhoi*Lf - deficit|/deficit, mass: ", &
         err_e, err_f, err_m
      call check(error, nbadf, 0, "T12 ledger is finite")
      if (allocated(error)) return
      call check(error, err_e < 1.0e-13_wp, "bank deficit + ice enthalpy = 0")
      if (allocated(error)) return
      call check(error, err_f < 1.0e-14_wp, "frazil volume x latent content = the bank")
      if (allocated(error)) return
      call check(error, err_m < 1.0e-13_wp, "mass booked to fresh = the ice formed")
      if (allocated(error)) return
      call check(error, st%bank(2, 1) == 5.0e4_wp .and. st%fv(2, 1) == 0.0_wp, "frzmlt < 0: bank retained, no new ice")
      if (allocated(error)) return
      call check(error, abs(e1(2) - e0(2)) <= 1.0e-13_wp*abs(e0(2)), "frzmlt < 0 and no growth: ice energy unchanged")
      if (allocated(error)) return
      call check(error, cell_identical(st, st0, 3, 1), "dry cell untouched, bank retained")
   end subroutine test_t12

   subroutine test_t13(error)
      !! The virtual-salt adapter `Q = 1000*fsalt - S*fresh` against the legacy `m*(S - S_ice)/dt` form with `S_ice` = the
      !! reference salinity of the flux convention: new ice (m > 0) and lateral melt (m < 0).
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, st0
      type(ice_params_t) :: p
      real(wp) :: dt, sss, m, adapter, legacy, worst
      integer :: c, nbadf

      dt = 3600.0_wp
      sss = 34.0_wp
      p = params_of(7, .true.)
      call st_alloc(st, 2, 1, 7)
      call set_simple_cell(st, p, 1, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      call set_simple_cell(st, p, 2, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.0_wp, -10.0_wp)
      st%frz(1, 1) = 0.1_wp*p%c%rhoi*p%c%Lfresh/dt; st%bank(1, 1) = st%frz(1, 1)*dt     ! new ice
      st%rside(2, 1) = 0.1_wp                                                           ! lateral melt (no snow)
      st0 = st
      call run_fused(st, p, dt, 1, 2, 1, 1)
      worst = 0.0_wp; nbadf = 0
      do c = 1, 2
         m = p%c%rhoi*(sum(st%v(c, 1, :)) - sum(st0%v(c, 1, :)))        ! signed ice mass change [kg/m2]
         adapter = 1000.0_wp*st%fs(c, 1) - sss*st%fr(c, 1)
         legacy = m*(sss - p%thermo%ice_ref_salinity)/dt
         call note(worst, abs(adapter - legacy)/abs(legacy), nbadf)
         call check(error, ieee_is_finite(adapter) .and. abs(adapter - legacy) <= 1.0e-13_wp*abs(legacy), &
                    "adapter equals the legacy form")
         if (allocated(error)) return
      end do
      write (output_unit, '(a,es10.2)') "T13: worst relative |adapter - legacy|: ", worst
      call check(error, (1000.0_wp*st%fs(1, 1) - sss*st%fr(1, 1)) > 0.0_wp .and. &
                 (1000.0_wp*st%fs(2, 1) - sss*st%fr(2, 1)) < 0.0_wp, "freezing salinifies, melting freshens")
   end subroutine test_t13

   subroutine test_l4(error)
      !! Order after the new ice: a cell holding ONLY new ice is melted laterally as new ice.  Closed forms.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(ice_params_t) :: p
      real(wp) :: dt, vv, r, rl, expect(6), got(6), worst
      integer :: k, nbadf

      dt = 3600.0_wp; vv = 0.1_wp; r = 0.1_wp
      p = params_of(7, .true.)
      rl = p%c%rhoi*p%c%Lfresh
      call st_alloc(st, 1, 1, 7)
      st%frz(1, 1) = vv*rl/dt; st%bank(1, 1) = st%frz(1, 1)*dt; st%rside(1, 1) = r
      call run_fused(st, p, dt, 1, 1, 1, 1)
      ! open water 1, hi0new = V (above hfrazilmin 0.05): the new ice covers the cell; the melt takes the fraction r of it
      got = [st%a(1, 1, 1), st%v(1, 1, 1), st%fh(1, 1), st%fr(1, 1), st%ml(1, 1), st%fs(1, 1)]
      expect = [1.0_wp - r, vv*(1.0_wp - r), -rl*vv*r/dt, -p%c%rhoi*vv*(1.0_wp - r)/dt, vv*r, &
                p%thermo%ice_ref_salinity*0.001_wp*(-p%c%rhoi*vv*(1.0_wp - r)/dt)]
      worst = 0.0_wp; nbadf = 0
      do k = 1, 6
         call note(worst, abs(got(k) - expect(k))/abs(expect(k)), nbadf)
      end do
      write (output_unit, '(a,es10.2)') "L4: worst relative difference from the closed forms: ", worst
      call check(error, nbadf, 0, "L4 closed forms are finite")
      if (allocated(error)) return
      call check(error, worst < 1.0e-13_wp, "area, volume, fhocn, fresh, meltl, fsalt of new ice melted laterally")
      if (allocated(error)) return
      worst = 0.0_wp
      do k = 1, 7
         call note(worst, abs(st%q(1, 1, 1, k) + rl)/rl, nbadf)
      end do
      call check(error, nbadf, 0, "L4 layer ledger is finite")
      if (allocated(error)) return
      call check(error, worst < 1.0e-14_wp, "the layers of the melted new ice keep the enthalpy -rhoi*Lfresh")
      if (allocated(error)) return
      call check(error, abs(st%T(1, 1, 1) - TF) < 1.0e-14_wp .and. abs(st%fv(1, 1) - vv) < 1.0e-14_wp, &
                 "surface temperature = Tf, frazil volume = V before the melt")
   end subroutine test_l4

   subroutine test_sign_split(error)
      !! Freezing-only, melting-only and neutral cells.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, st0
      type(ice_params_t) :: p
      real(wp) :: dt, rl, vv
      integer :: nbad
      real(wp) :: worst

      dt = 3600.0_wp; vv = 0.05_wp
      p = params_of(7, .true.)
      rl = p%c%rhoi*p%c%Lfresh
      call st_alloc(st, 3, 1, 7)
      call set_simple_cell(st, p, 1, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.1_wp, -10.0_wp)
      call set_simple_cell(st, p, 2, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.1_wp, -10.0_wp)
      call set_simple_cell(st, p, 3, 1, [0.2_wp, 0.2_wp, 0.1_wp, 0.1_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp], &
                           0.1_wp, -10.0_wp)
      st%fh = 7.25_wp; st%fr = 0.125_wp; st%fs = -0.5_wp
      st%frz(1, 1) = vv*rl/dt; st%bank(1, 1) = st%frz(1, 1)*dt                      ! freezing only
      st%frz(2, 1) = -60.0_wp; st%bank(2, 1) = 2.0e4_wp; st%rside(2, 1) = 0.08_wp   ! melting only
      st0 = st
      call run_fused(st, p, dt, 1, 3, 1, 1)
      ! freezing only: ice volume up by V, fhocn and the melt diagnostics untouched, bank spent
      call check(error, abs((sum(st%v(1, 1, :)) - sum(st0%v(1, 1, :))) - vv) < 1.0e-13_wp*vv, "freezing: ice volume rises by V")
      if (allocated(error)) return
      call check(error, st%fh(1, 1) == 7.25_wp .and. st%ml(1, 1) == 0.0_wp .and. st%bank(1, 1) == 0.0_wp, &
                 "freezing: fhocn and meltl untouched, bank spent")
      if (allocated(error)) return
      ! melting only: no new ice, the stale bank retained, heat taken from the ocean (fhocn falls), mass to the ocean
      call check(error, st%fv(2, 1) == 0.0_wp .and. st%bank(2, 1) == 2.0e4_wp, "melting: no new ice, bank retained")
      if (allocated(error)) return
      call check(error, st%fh(2, 1) < 7.25_wp .and. st%fr(2, 1) > 0.125_wp .and. st%ml(2, 1) > 0.0_wp, &
                 "melting: fhocn falls, fresh and meltl rise")
      if (allocated(error)) return
      ! neutral: nothing is added to any flux and the clean stationary column stays put
      call check(error, st%fh(3, 1) == 7.25_wp .and. st%fr(3, 1) == 0.125_wp .and. st%fs(3, 1) == -0.5_wp .and. &
                 st%ml(3, 1) == 0.0_wp .and. st%fv(3, 1) == 0.0_wp, "neutral: no flux increment")
      if (allocated(error)) return
      st%wet = .true.
      call st_compare(st, st0, 1.0e-14_wp, nbad, worst)
      write (output_unit, '(a,es10.2)') "sign split: neutral cell worst relative state change ", &
         maxval([abs(st%v(3, 1, :) - st0%v(3, 1, :)), abs(st%a(3, 1, :) - st0%a(3, 1, :))])
      call check(error, all(abs(st%v(3, 1, :) - st0%v(3, 1, :)) <= 1.0e-14_wp*max(st0%v(3, 1, :), 1.0e-30_wp)) .and. &
                 all(abs(st%a(3, 1, :) - st0%a(3, 1, :)) <= 1.0e-14_wp), "neutral: the clean state is unchanged")
   end subroutine test_sign_split

   subroutine test_rside_zero(error)
      !! `rside = 0` everywhere: the fused kernel equals the staged run WITHOUT the lateral-melt stage, and the melt
      !! diagnostic is exactly untouched.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: s1, s2
      type(ice_params_t) :: p
      integer :: nbad
      real(wp) :: worst

      p = params_of(7, .true.)
      call st_alloc(s1, 40, 30, 7)
      call fill_state(s1, 777777_int64, 1)
      s1%wet = .true.; s1%rside = 0.0_wp; s1%ml = 0.5_wp
      s2 = s1
      call run_fused(s1, p, 2400.0_wp, 1, 40, 1, 30)
      call run_staged(s2, p, 2400.0_wp, 1, 40, 1, 30, [1, 2, 3, 5])
      call st_compare(s1, s2, TOL_ORD, nbad, worst)
      call check(error, nbad, 0, "rside = 0: the lateral stage is a no-op")
      if (allocated(error)) return
      call check(error, all(s1%ml == 0.5_wp), "rside = 0: meltl bitwise untouched")
   end subroutine test_rside_zero

   subroutine test_aggregates(error)
      !! The aggregates left by the driver are the exact aggregates of the state it left (the lateral melt does not refresh
      !! `aice0`; the clean-up's last action does).
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(ice_params_t) :: p
      integer :: i, j, n, nbad, nbadf
      real(wp) :: sa, sv, ss, worst

      p = params_of(7, .true.)
      call st_alloc(st, 50, 40, 7)
      call fill_state(st, 135791_int64, 1)
      st%wet = .true.
      call run_fused(st, p, 2400.0_wp, 1, 50, 1, 40)
      nbad = 0; worst = 0.0_wp; nbadf = 0
      do j = 1, 40
         do i = 1, 50
            if (iand(st%ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) then
               nbad = nbad + 1
               cycle
            end if
            sa = 0.0_wp; sv = 0.0_wp; ss = 0.0_wp
            do n = 1, NCAT
               sa = sa + st%a(i, j, n); sv = sv + st%v(i, j, n); ss = ss + st%s(i, j, n)
            end do
            call note(worst, abs(st%aice(i, j) - sa), nbadf); call note(worst, abs(st%vice(i, j) - sv), nbadf)
            call note(worst, abs(st%vsno(i, j) - ss), nbadf); call note(worst, abs(st%aice0(i, j) - max(1.0_wp - sa, 0.0_wp)), nbadf)
         end do
      end do
      write (output_unit, '(a,es10.2,a,i0)') "aggregates: worst absolute difference from a recomputation ", worst, &
         ", fatal cells ", nbad
      call check(error, nbad, 0, "no fatal bit on the random state")
      if (allocated(error)) return
      call check(error, nbadf, 0, "the aggregates are finite")
      if (allocated(error)) return
      call check(error, worst <= 1.0e-14_wp, "aice, vice, vsno, aice0 are the aggregates of the final state")
   end subroutine test_aggregates

   subroutine test_fatal_untouched(error)
      !! A cell that already carries a FATAL bit is skipped by every stage: bitwise untouched (state, fluxes, bank).
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, st0
      type(ice_params_t) :: p
      integer :: i, nchanged

      p = params_of(7, .true.)
      call st_alloc(st, 8, 1, 7)
      call fill_state(st, 99991_int64, 1)
      st%wet = .true.
      do i = 1, 8
         st%frz(i, 1) = 150.0_wp; st%bank(i, 1) = 150.0_wp*2400.0_wp; st%rside(i, 1) = 0.05_wp; st%frain(i, 1) = 2.0e-5_wp
      end do
      st%ierr(3, 1) = ICE_ERR_NEG_THICK
      st%ierr(6, 1) = ior(ICE_ERR_NEG_THICK, ICE_ERR_REMAP_FALLBACK)
      st0 = st
      call run_fused(st, p, 2400.0_wp, 1, 8, 1, 1)
      call check(error, cell_identical(st, st0, 3, 1) .and. cell_identical(st, st0, 6, 1), "fatal-bit cells are untouched bitwise")
      if (allocated(error)) return
      nchanged = 0
      do i = 1, 8
         if (i /= 3 .and. i /= 6 .and. .not. cell_identical(st, st0, i, 1)) nchanged = nchanged + 1
      end do
      call check(error, nchanged >= 4, "the other cells are processed")
   end subroutine test_fatal_untouched

   subroutine test_nan(error)
      !! A NaN in one cell's input flags that cell only: a NaN scalar leaves the cell untouched, a NaN in the state is
      !! flagged; the healthy cells equal a run without any NaN.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, clean, st0
      type(ice_params_t) :: p
      real(wp) :: qn
      integer :: i

      p = params_of(7, .true.)
      qn = ieee_value(1.0_wp, ieee_quiet_nan)
      call st_alloc(st, 8, 1, 7)
      call fill_state(st, 31337_int64, 1)
      st%wet = .true.
      do i = 1, 8
         st%frz(i, 1) = 100.0_wp; st%bank(i, 1) = 100.0_wp*2400.0_wp; st%rside(i, 1) = 0.03_wp; st%frain(i, 1) = 1.0e-5_wp
         st%a(i, 1, 1) = 0.2_wp; st%ai0(i, 1, 1) = 0.2_wp; st%v(i, 1, 1) = 0.1_wp; st%vi0(i, 1, 1) = 0.08_wp
         st%T(i, 1, 1) = -5.0_wp; st%q(i, 1, 1, :) = -2.5e8_wp; st%si(i, 1, 1, :) = 3.0_wp
         st%aice(i, 1) = sum(st%ai0(i, 1, :))
      end do
      clean = st
      call run_fused(clean, p, 2400.0_wp, 1, 8, 1, 1)
      st%frain(2, 1) = qn
      st%rside(3, 1) = qn
      st%frz(4, 1) = qn
      st%a(5, 1, 1) = qn
      st%q(6, 1, 1, 2) = qn
      st%aice(7, 1) = qn
      st0 = st
      st0%ierr(2, 1) = ICE_ERR_NONFINITE; st0%ierr(3, 1) = ICE_ERR_NONFINITE      ! the only change: the flag
      st0%ierr(4, 1) = ICE_ERR_NONFINITE; st0%ierr(7, 1) = ICE_ERR_NONFINITE
      call run_fused(st, p, 2400.0_wp, 1, 8, 1, 1)
      do i = 1, 8
         if (i == 1 .or. i == 8) then
            call check(error, cell_identical(st, clean, i, 1), "healthy cells equal the run without NaN")
         else
            call check(error, iand(st%ierr(i, 1), ICE_ERR_NONFINITE) /= 0_int32, "the NaN cell is flagged")
            if (allocated(error)) return
            if (i == 2 .or. i == 3 .or. i == 4 .or. i == 7) &
               call check(error, cell_identical(st, st0, i, 1), "a NaN scalar leaves the cell untouched")
         end if
         if (allocated(error)) return
      end do
    call check(error, iand(st%ierr(1, 1), ICE_ERR_FATAL_MASK) == 0_int32 .and. iand(st%ierr(8, 1), ICE_ERR_FATAL_MASK) == 0_int32, &
                 "no flag on the healthy cells")
   end subroutine test_nan

   subroutine test_fallback_bins(error)
      !! Remap-fallback cells (growth or melt too fast for the remap) end the step inside their bins.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(ice_params_t) :: p
      integer :: i, j, n, nfb, nviol, ncat_checked
      real(wp) :: h

      p = params_of(7, .true.)
      call st_alloc(st, 200, 100, 7)
      call fill_state(st, 2468013_int64, 1)
      st%wet = .true.
      call run_fused(st, p, 2400.0_wp, 1, 200, 1, 100)
      nfb = 0; nviol = 0; ncat_checked = 0
      do j = 1, 100
         do i = 1, 200
            if (iand(st%ierr(i, j), ICE_ERR_REMAP_FALLBACK) /= 0_int32) nfb = nfb + 1
            do n = 1, NCAT
               if (st%a(i, j, n) > PUNY) then
                  ncat_checked = ncat_checked + 1
                  h = st%v(i, j, n)/st%a(i, j, n)
                  if (h < HIN(n - 1)*(1.0_wp - 1.0e-12_wp) - 1.0e-14_wp .or. h > HIN(n)*(1.0_wp + 1.0e-12_wp)) nviol = nviol + 1
               end if
            end do
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0)') "fallback: cells with a remap fallback ", nfb, ", populated categories ", &
         ncat_checked, ", outside their bin ", nviol
      call check(error, nfb > 100, "the state contains remap fallbacks")
      if (allocated(error)) return
      call check(error, nviol, 0, "every populated category is inside its bin after the step")
   end subroutine test_fallback_bins

   ! ------------------------------------------------------------------
   ! conservation
   ! ------------------------------------------------------------------

   subroutine test_conserve(error)
      !! Mass (ice + snow), salt and energy over the WHOLE step on 100 000 random cells against the flux increments and the
      !! frazil input: `dM + dt*(fresh_inc - rain) = 0`, `dS + dt*fsalt_inc = 0`, `dE + dt*fhocn_inc - qi0new*V_frazil = 0`
      !! (the new ice brings its enthalpy `-rhoi*Lfresh` per volume from the ocean's bank, not through `fhocn`), and the
      !! frazil volume against the bank.  Worst relative residuals (to the sum of the magnitudes) are printed.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 500, ny = 200
      type(st_t) :: st, st0
      type(ice_params_t) :: p
      integer :: i, j, n, k, nfatal, nbank, nbadf
      real(wp) :: dt, rhoi, rhos, sref, m0, m1, s0, s1, e0, e1, dfr, rm, rs, re, sm, ss, se, w(4), qi0, rbank, nrain

      dt = 2400.0_wp
      p = params_of(7, .true.)
      rhoi = p%c%rhoi; rhos = p%c%rhos; sref = p%thermo%ice_ref_salinity; qi0 = -rhoi*p%c%Lfresh
      call st_alloc(st, nx, ny, 7)
      call fill_state(st, 8675309_int64, 1)
      st%wet = .true.
      st0 = st
      call run_fused(st, p, dt, 1, nx, 1, ny)
      w = 0.0_wp; nfatal = 0; nbank = 0; nbadf = 0
      do j = 1, ny
         do i = 1, nx
            if (iand(st%ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) then
               nfatal = nfatal + 1
               cycle
            end if
            m0 = 0.0_wp; m1 = 0.0_wp; s0 = 0.0_wp; s1 = 0.0_wp
            do n = 1, NCAT
               m0 = m0 + rhoi*st0%v(i, j, n) + rhos*st0%s(i, j, n); m1 = m1 + rhoi*st%v(i, j, n) + rhos*st%s(i, j, n)
               s0 = s0 + rhoi*st0%v(i, j, n)*sref*0.001_wp; s1 = s1 + rhoi*st%v(i, j, n)*sref*0.001_wp
            end do
            e0 = cell_energy(st0, i, j); e1 = cell_energy(st, i, j)
            dfr = st%fr(i, j) - st0%frain(i, j)*st0%aice(i, j)
            sm = abs(m0) + abs(m1) + dt*abs(dfr) + rhoi*st%fv(i, j)       ! the gross terms: new ice melted again cancels in the flux
            ss = abs(s0) + abs(s1) + dt*abs(st%fs(i, j)) + rhoi*st%fv(i, j)*sref*0.001_wp
            se = abs(e0) + abs(e1) + dt*abs(st%fh(i, j)) + abs(qi0*st%fv(i, j))
            rm = abs((m1 - m0) + dt*dfr)/max(sm, tiny(1.0_wp))
            rs = abs((s1 - s0) + dt*st%fs(i, j))/max(ss, tiny(1.0_wp))
            re = abs((e1 - e0) + dt*st%fh(i, j) - qi0*st%fv(i, j))/max(se, tiny(1.0_wp))
            call note(w(1), rm, nbadf); call note(w(2), rs, nbadf); call note(w(3), re, nbadf)
            if (st0%frz(i, j) > 0.0_wp) then
               nbank = nbank + 1
               rbank = abs(st%fv(i, j)*rhoi*p%c%Lfresh - st0%bank(i, j))/st0%bank(i, j)
               call note(w(4), rbank, nbadf)
               if (st%bank(i, j) /= 0.0_wp) nfatal = nfatal + 1000000
            else
               if (st%bank(i, j) /= st0%bank(i, j) .or. st%fv(i, j) /= st0%fv(i, j)) nfatal = nfatal + 1000000
            end if
         end do
      end do
      nrain = real(count(st0%frain > 0.0_wp), wp)
      write (output_unit, '(a,i0,a,i0,a)') "conservation: ", nx*ny, " cells, ", nbank, " freezing cells"
      write (output_unit, '(a,4es10.2)') "  worst relative residual: mass, salt, energy, frazil-vs-bank: ", w
      call check(error, nbadf, 0, "the ledgers are finite")
      if (allocated(error)) return
      call check(error, nfatal, 0, "no fatal bit; the bank is spent exactly where frzmlt > 0")
      if (allocated(error)) return
      call check(error, w(1) < CLOSE_TOL, "mass closes")
      if (allocated(error)) return
      call check(error, w(2) < CLOSE_TOL, "salt closes")
      if (allocated(error)) return
      call check(error, w(3) < CLOSE_TOL, "energy closes")
      if (allocated(error)) return
      call check(error, w(4) < 2.0e-15_wp*8.0_wp, "the frazil volume is the bank")
   end subroutine test_conserve

   ! ------------------------------------------------------------------
   ! Icepack oracle
   ! ------------------------------------------------------------------

   subroutine oracle_run(error, kset, gin, gold)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kset
      real(wp), intent(in) :: gin(126, NCO), gold(102, NCO)
      type(st_t) :: st
      type(ice_params_t) :: p
      integer :: c, nilyr, ntr, nout, bad, it, nfb, nnew, nlat, nrain, nzap, nempty, nup, ndown, nrebin, n
      real(wp) :: dt, sal(7), tr(ICE_NTRCR_MAX, NCAT), ai(NCAT), vi(NCAT), a(NCAT), v(NCAT), s(NCAT), got(102), worst, d, aic
      real(wp) :: vic, vso, a0c, trm(ICE_NTRCR_MAX), ex(ICE_NTRCR_MAX, NCAT), ao(NCAT), vo(NCAT), so(NCAT)
      integer(int32) :: ie2
      character(len=200) :: msg

      call fill_golden()
      worst = 0.0_wp
      nfb = 0; nnew = 0; nlat = 0; nrain = 0; nzap = 0; nempty = 0; nup = 0; ndown = 0; nrebin = 0
      do c = 1, NCO
         nilyr = GNIL(c)
         ntr = ice_ntrcr(nilyr, 1)
         p = ice_params_default(kset, NCAT, nilyr, 1)
         p%thermo%hfrazilmin = gin(3, c)
         p%thermo%update_ocn_f = gin(4, c) > 0.5_wp
         p%thermo%ice_ref_salinity = gin(5, c)
         dt = gin(1, c)
         sal = gin(15:21, c)
         ai = gin(22:26, c); vi = gin(27:31, c); a = gin(32:36, c); v = gin(37:41, c); s = gin(42:46, c)
         tr = 0.0_wp
         tr(1:ntr, :) = reshape(gin(47:46 + ntr*NCAT, c), [ntr, NCAT])
         call st_alloc(st, 1, 1, nilyr)
         call st_set_cell(st, 1, 1, ai, vi, a, v, s, tr)
         st%Tf(1, 1) = gin(2, c); st%frz(1, 1) = gin(6, c); st%frain(1, 1) = gin(7, c); st%rside(1, 1) = gin(8, c)
         st%aice(1, 1) = gin(9, c); st%fr(1, 1) = gin(10, c); st%fs(1, 1) = gin(11, c); st%fh(1, 1) = gin(12, c)
         st%ml(1, 1) = gin(13, c); st%fv(1, 1) = gin(14, c)
         if (gin(6, c) > 0.0_wp) st%bank(1, 1) = gin(6, c)*dt
         call run_fused(st, p, dt, 1, 1, 1, 1, sal)
         call check(error, iand(st%ierr(1, 1), ICE_ERR_FATAL_MASK), 0, "oracle case: no fatal bit")
         if (allocated(error)) return
         if (iand(st%ierr(1, 1), ICE_ERR_REMAP_FALLBACK) /= 0_int32) nfb = nfb + 1
         call st_get_cell(st, 1, 1, a, v, s, tr)
         got = 0.0_wp
         got(1:5) = a; got(6:10) = v; got(11:15) = s
         got(16:15 + ntr*NCAT) = reshape(tr(1:ntr, :), [ntr*NCAT])
         nout = 15 + ntr*NCAT
         got(nout + 3) = st%fr(1, 1); got(nout + 4) = st%fs(1, 1); got(nout + 5) = st%fh(1, 1)
         got(nout + 6) = st%ml(1, 1); got(nout + 7) = st%fv(1, 1)
         do it = 1, nout + 7
            if (it == nout + 1 .or. it == nout + 2) cycle      ! upstream's intermediate aggregates
            if (gold(it, c) /= 0.0_wp .and. ieee_is_finite(got(it))) worst = max(worst, abs(got(it) - gold(it, c))/abs(gold(it, c)))
            if (.not. (ieee_is_finite(got(it)) .and. same_gold(got(it), gold(it, c)))) then
               bad = it
               write (msg, '(a,i0,a,i0,a,es24.16,a,es24.16)') "oracle case ", c, " output ", bad, " got ", got(it), &
                  " expected ", gold(it, c)
               call test_failed(error, trim(msg))
               return
            end if
         end do
         ! the cell aggregates equal the aggregates of the oracle's output state
         ao = gold(1:5, c); vo = gold(6:10, c); so = gold(11:15, c)
         ex = 0.0_wp
         ex(1:ntr, :) = reshape(gold(16:15 + ntr*NCAT, c), [ntr, NCAT])
         call itd_aggregate(NCAT, ICE_NTRCR_MAX, nilyr, 1, ao, vo, so, ex, gin(2, c), PUNY, aic, vic, vso, a0c, trm, ie2)
         call check(error, same_gold(st%aice(1, 1), aic) .and. same_gold(st%vice(1, 1), vic) .and. same_gold(st%vsno(1, 1), vso), &
                    "oracle case: aice, vice, vsno are the aggregates of the oracle state")
         if (allocated(error)) return
         call check(error, abs(st%aice0(1, 1) - a0c) <= max(TOL_GOLD, 1.0e-15_wp), "oracle case: aice0")
         if (allocated(error)) return
         ! coverage
         if (gin(6, c) > 0.0_wp) nnew = nnew + 1
         if (gin(8, c) > 0.0_wp) nlat = nlat + 1
         if (gin(7, c) > 0.0_wp) nrain = nrain + 1
         if (all(gin(22:26, c) == 0.0_wp)) nempty = nempty + 1
         do n = 1, NCAT
            if (gin(32 + n - 1, c) > 0.0_wp .and. gold(n, c) == 0.0_wp .and. gin(8, c) < 1.0_wp) nzap = nzap + 1
         end do
         if (gin(22, c) > 0.0_wp .and. gold(1, c) < gin(22, c) .and. any(gold(2:5, c) > gin(23:26, c))) nup = nup + 1
         if (any(gold(1:5, c) < gin(32:36, c)*(1.0_wp - 0.001_wp)) .and. gin(8, c) == 0.0_wp) ndown = ndown + 1
         nrebin = nrebin
      end do
      write (output_unit, '(a,es10.2)') "oracle worst relative difference (every output of every case) ", worst
      write (output_unit, '(a,8i4)') "oracle coverage (new ice, lateral melt, rain, ice-free, fallback, area moved up, area lost or moved down, small cats zapped): ", &
         nnew, nlat, nrain, nempty, nfb, nup, ndown, nzap
      call check(error, nnew >= 12 .and. nlat >= 8 .and. nrain >= 6 .and. nempty >= 2 .and. nfb >= 2 .and. nup >= 3 .and. &
                 ndown >= 3, "oracle set covers every stage")
   end subroutine oracle_run

   pure logical function same_gold(x, y)
      real(wp), intent(in) :: x, y

      if (TOL_GOLD == 0.0_wp) then
         same_gold = x == y
      else
         same_gold = abs(x - y) <= TOL_GOLD*max(abs(x), abs(y), 1.0_wp)
      end if
   end function same_gold

   subroutine test_oracle_icepack(error)
      type(error_type), allocatable, intent(out) :: error

      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_ICEPACK, GII, GIO)
   end subroutine test_oracle_icepack

   subroutine test_oracle_cesm(error)
      type(error_type), allocatable, intent(out) :: error

      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_CESM, GCI, GCO)
   end subroutine test_oracle_cesm

   ! ------------------------------------------------------------------
   ! device sweep and throughput
   ! ------------------------------------------------------------------

   subroutine test_device(error)
      !! The host shim over an INTERIOR range of a random 64 x 48 state (with dry cells) against the host cell kernels run
      !! serially.  Exact: error words, the zero pattern of every category, the new-ice decisions (a cell formed new ice;
      !! the final open water exceeds `puny`).  Category area, ice volume and snow volume: against `tol x` the CELL's aice,
      !! vice, vsno.  Flux increments: against their gross magnitude (sum of the absolute terms).  Intensive tracers only
      !! where the category holds `> 1e-3` of the cell's ice area.  The ring and the dry cells must be untouched.  The worst
      !! of each class and where it is are printed.  Then the throughput.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 64, ny = 48
      type(st_t) :: dev, ref, orig
      type(ice_params_t) :: p
      integer :: i, j, n, it, nring, nflip, nbadf
      real(wp) :: wv(4), dt, rhoi, rhos, sref, lf, sumv0, sums0, eabs, fv, g_fr, g_fs, g_fh, g_ml, fl
      character(len=16) :: wn(4)
      integer :: wi(4), wj(4)
      logical :: inside

      dt = 2400.0_wp
      p = params_of(7, .true.)
      rhoi = p%c%rhoi; rhos = p%c%rhos; sref = p%thermo%ice_ref_salinity; lf = p%c%Lfresh
      call st_alloc(dev, nx, ny, 7)
      call fill_state(dev, 515151_int64, 1)
      orig = dev
      ref = dev
      call run_fused(dev, p, dt, 2, nx - 1, 2, ny - 1)
      call run_host_cells(ref, p, dt, 2, nx - 1, 2, ny - 1)
      nring = 0; nflip = 0; nbadf = 0
      wv = 0.0_wp; wi = 0; wj = 0; wn = ''
      do j = 1, ny
         do i = 1, nx
            inside = i >= 2 .and. i <= nx - 1 .and. j >= 2 .and. j <= ny - 1
            if (.not. inside .or. .not. orig%wet(i, j)) then
               if (.not. cell_identical(dev, orig, i, j)) nring = nring + 1
               cycle
            end if
            if (dev%ierr(i, j) /= ref%ierr(i, j)) then
               nflip = nflip + 1
               write (output_unit, '(a,2i4,2i9)') "   FLIP error word, cell ", i, j, ref%ierr(i, j), dev%ierr(i, j)
            end if
            if (iand(ref%ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) cycle
  if (((ref%aice0(i, j) > PUNY) .neqv. (dev%aice0(i, j) > PUNY)) .or. ((ref%fv(i, j) > 0.0_wp) .neqv. (dev%fv(i, j) > 0.0_wp))) then
               nflip = nflip + 1
               write (output_unit, '(a,2i4)') "   FLIP new-ice / open-water decision, cell ", i, j
            end if
            ! gross magnitudes from the input state
            sumv0 = sum(orig%v(i, j, :)); sums0 = sum(orig%s(i, j, :)); fv = ref%fv(i, j)
            eabs = 0.0_wp
            do n = 1, NCAT
               eabs = eabs + sum(abs(orig%q(i, j, n, :)))*orig%v(i, j, n)/7.0_wp + abs(orig%qs(i, j, n, 1))*orig%s(i, j, n)
            end do
            g_fr = (rhoi*sumv0 + rhos*sums0 + rhoi*fv)/dt + orig%frain(i, j)*orig%aice(i, j)
            g_fs = rhoi*(sumv0 + fv)*sref*0.001_wp/dt
            g_fh = (eabs + rhoi*lf*fv)/dt
            g_ml = sumv0 + fv
            do n = 1, NCAT
               if (((ref%a(i, j, n) == 0.0_wp) .neqv. (dev%a(i, j, n) == 0.0_wp)) .or. &
                   ((ref%v(i, j, n) == 0.0_wp) .neqv. (dev%v(i, j, n) == 0.0_wp)) .or. &
                   ((ref%s(i, j, n) == 0.0_wp) .neqv. (dev%s(i, j, n) == 0.0_wp))) then
                  nflip = nflip + 1
                  if (nflip <= 5) write (output_unit, '(a,3i4)') "   FLIP zero pattern, cell/cat ", i, j, n
               end if
               call dif(1, "area", ref%a(i, j, n), dev%a(i, j, n), ref%aice(i, j), i, j)
               call dif(1, "ice volume", ref%v(i, j, n), dev%v(i, j, n), ref%vice(i, j), i, j)
               call dif(1, "snow volume", ref%s(i, j, n), dev%s(i, j, n), ref%vsno(i, j), i, j)
               if (ref%a(i, j, n) > 1.0e-3_wp*ref%aice(i, j)) then
                  do it = 1, ice_ntrcr(7, 1)
                     fl = 0.0_wp
                     if (it == 1) fl = 1.0_wp
                     call dif(2, "tracer", ice_trcr_get(it, i, j, n, nx, ny, NCAT, 7, 1, ref%T, ref%q, ref%si, ref%qs), &
                              ice_trcr_get(it, i, j, n, nx, ny, NCAT, 7, 1, dev%T, dev%q, dev%si, dev%qs), &
                              max(abs(ice_trcr_get(it, i, j, n, nx, ny, NCAT, 7, 1, ref%T, ref%q, ref%si, ref%qs)), fl), i, j)
                  end do
               end if
            end do
            call dif(4, "fresh", ref%fr(i, j), dev%fr(i, j), g_fr, i, j)
            call dif(4, "fsalt", ref%fs(i, j), dev%fs(i, j), g_fs, i, j)
            call dif(4, "fhocn", ref%fh(i, j), dev%fh(i, j), g_fh, i, j)
            call dif(4, "meltl", ref%ml(i, j), dev%ml(i, j), g_ml, i, j)
            call dif(4, "frazil", ref%fv(i, j), dev%fv(i, j), max(abs(fv), tiny(1.0_wp)), i, j)
            call dif(1, "aice", ref%aice(i, j), dev%aice(i, j), ref%aice(i, j), i, j)
            call dif(1, "vice", ref%vice(i, j), dev%vice(i, j), ref%vice(i, j), i, j)
            call dif(1, "vsno", ref%vsno(i, j), dev%vsno(i, j), ref%vsno(i, j), i, j)
            call dif(3, "aice0", ref%aice0(i, j), dev%aice0(i, j), 1.0_wp, i, j)
            call check(error, ref%bank(i, j) == dev%bank(i, j), "the bank is spent identically")
            if (allocated(error)) return
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0)') "therm2 device sweep: decision flips ", nflip, ", ring/dry violations ", nring, &
         ", non-finite differences ", nbadf
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst area/volume difference / cell aice,vice,vsno ", wv(1), " at ", trim(wn(1)), wi(1), wj(1)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst tracer relative difference (a_n > 1e-3 aice) ", wv(2), " at ", trim(wn(2)), wi(2), wj(2)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst aice0 absolute difference ", wv(3), " at ", trim(wn(3)), wi(3), wj(3)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst flux increment difference / gross magnitude ", wv(4), " at ", trim(wn(4)), wi(4), wj(4)
      call check(error, nring, 0, "the ring and the dry cells are untouched (all outputs)")
      if (allocated(error)) return
      call check(error, nflip + nbadf, 0, "every error word, zero pattern and decision equals the host reference")
      if (allocated(error)) return
      call check(error, wv(1) <= TOL_EXT, "area and volumes agree against the cell's totals")
      if (allocated(error)) return
      call check(error, wv(2) <= TOL_INT, "intensive tracers of healthy categories agree relatively")
      if (allocated(error)) return
      call check(error, wv(3) <= 1.0e-15_wp, "aice0 agrees to round-off of 1")
      if (allocated(error)) return
      call check(error, wv(4) <= TOL_FLX, "flux increments agree against their gross magnitude")
      if (allocated(error)) return
      call throughput(error, p, 1)
      if (allocated(error)) return
      call throughput(error, p, 2)
   contains
      subroutine dif(cls, name, rv, gv, scale_f, ci, cj)
         integer, intent(in) :: cls, ci, cj
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: rv, gv, scale_f
         real(wp) :: dd

         if (.not. (ieee_is_finite(rv) .and. ieee_is_finite(gv) .and. ieee_is_finite(scale_f))) then
            nbadf = nbadf + 1
            return
         end if
         if (rv == gv) return
         dd = abs(rv - gv)/max(abs(scale_f), tiny(1.0_wp))
         if (dd > wv(cls)) then
            wv(cls) = dd; wn(cls) = name; wi(cls) = ci; wj(cls) = cj
         end if
      end subroutine dif
   end subroutine test_device

   subroutine restore_state(nx, ny, nilyr, a, v, s, T, q, si, qs, aice, bank, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, aice_s, bank_s, &
                            fr, fs, fh, ml, fv)
      !! Reset the sweep's inputs and zero its flux accumulators on the device.
      integer, intent(in) :: nx, ny, nilyr
      real(wp), intent(inout) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), T(nx, ny, NCAT)
      real(wp), intent(inout) :: q(nx, ny, NCAT, nilyr), si(nx, ny, NCAT, nilyr), qs(nx, ny, NCAT, 1), aice(nx, ny), bank(nx, ny)
      real(wp), intent(in) :: a_s(nx, ny, NCAT), v_s(nx, ny, NCAT), s_s(nx, ny, NCAT), T_s(nx, ny, NCAT)
  real(wp), intent(in) :: q_s(nx, ny, NCAT, nilyr), si_s(nx, ny, NCAT, nilyr), qs_s(nx, ny, NCAT, 1), aice_s(nx, ny), bank_s(nx, ny)
      real(wp), intent(inout) :: fr(nx, ny), fs(nx, ny), fh(nx, ny), ml(nx, ny), fv(nx, ny)

      !$acc kernels present(a, v, s, T, q, si, qs, aice, bank, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, aice_s, bank_s, fr, fs, fh, ml, fv)
      a = a_s; v = v_s; s = s_s; T = T_s; q = q_s; si = si_s; qs = qs_s; aice = aice_s; bank = bank_s
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp; ml = 0.0_wp; fv = 0.0_wp
      !$acc end kernels
   end subroutine restore_state

   subroutine throughput(error, p, kind)
      !! 512 x 512 cells x 5 categories x 7 layers, repeated sweeps on re-seeded data, for three launch structures: the
      !! shim's four launches, the shim with the new ice and lateral melt fused (three), and the five stage shims of the
      !! owning modules.  `kind = 1`: the wild random state; `kind = 2`: the realistic one (~85 % ice free in coherent
      !! blocks, clean nearly stationary ice).
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: kind
      integer, parameter :: nx = 512, ny = 512, nilyr = 7, nrep = 5
      type(st_t) :: st, sv
      real(wp) :: hin_l(0:NCAT), sal(nilyr), secs(3)
      integer(int64) :: t0, t1, rate
      integer :: k, v, nfatal
      character(len=24), parameter :: LBL(3) = [character(len=24) :: "4 launches (shim)", "3 launches (B+C fused)", "5 launches (stage shims)"]

      call st_alloc(st, nx, ny, nilyr)
      call fill_state(st, merge(616161_int64, 818181_int64, kind == 1), kind)
      st%wet = .true.
      sv = st
      hin_l = HIN
      sal = salinz_of(p, nilyr)
      call st_enter(st)
      call st_enter(sv)
      !$acc enter data copyin(hin_l, sal)
      secs = 0.0_wp
      call system_clock(count_rate=rate)
      do v = 1, 3
         do k = 0, nrep
            call restore_state(nx, ny, nilyr, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, st%aice, st%bank, sv%a, sv%v, sv%s, &
                               sv%T, sv%q, sv%si, sv%qs, sv%aice, sv%bank, st%fr, st%fs, st%fh, st%ml, st%fv)
            call system_clock(t0)
            if (v == 3) then
               call staged_stages(st, p, 2400.0_wp, 1, nx, 1, ny, [1, 2, 3, 4, 5], hin_l, sal)
            else
               call ice_step_therm2(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, 2400.0_wp, p, hin_l, sal, st%wet, st%Tf, st%frz, &
                                    st%frain, st%rside, st%bank, st%ai0, st%vi0, st%a, st%v, st%s, st%T, st%q, st%si, st%qs, &
                                    st%aice, st%vice, st%vsno, st%aice0, st%Tc, st%qc, st%sc, st%qsc, st%fr, st%fs, st%fh, &
                                    st%ml, st%fv, st%ierr, v == 2)
            end if
            !$acc wait
            call system_clock(t1)
            if (k > 0) secs(v) = secs(v) + real(t1 - t0, wp)/real(rate, wp)
         end do
      end do
      !$acc exit data delete(hin_l, sal)
      call st_exit(st)
      call st_exit(sv)
      nfatal = count(iand(st%ierr, ICE_ERR_FATAL_MASK) /= 0_int32)
      do v = 1, 3
         write (output_unit, '(a,i0,a,a,f10.3,a,es11.3,a)') "therm2 sweep 512x512x5x7, state ", kind, ": ", LBL(v), &
            secs(v)/real(nrep, wp)*1.0e3_wp, " ms/sweep, ", real(nx*ny, wp)*real(nrep, wp)/secs(v), " cells/s"
      end do
      call check(error, nfatal, 0, "no fatal bit in the throughput state")
   end subroutine throughput

   subroutine test_midway_fatal(error)
      !! A finite `rside = 1.5` passes the finite check of stage A but makes the lateral melt raise its range bit: the
      !! clean-up must then NOT run on that cell (its aggregates stay untouched), while a healthy neighbour is cleaned.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(ice_params_t) :: p
      integer :: c

      p = params_of(7, .true.)
      call st_alloc(st, 2, 1, 7)
      do c = 1, 2
         call set_simple_cell(st, p, c, 1, [0.2_wp, 0.3_wp, 0.1_wp, 0.0_wp, 0.0_wp], [0.3_wp, 1.0_wp, 1.9_wp, 0.0_wp, 0.0_wp], &
                              0.1_wp, -10.0_wp)
      end do
      st%rside(1, 1) = 1.5_wp; st%rside(2, 1) = 0.05_wp
      call run_fused(st, p, 3600.0_wp, 1, 2, 1, 1)
      call check(error, iand(st%ierr(1, 1), ICE_ERR_AICE_RANGE) /= 0_int32, "the lateral melt raised its range bit")
      if (allocated(error)) return
      call check(error, st%vice(1, 1) == -777.0_wp .and. st%aice0(1, 1) == -777.0_wp .and. st%vsno(1, 1) == -777.0_wp .and. &
                 st%Tc(1, 1) == -777.0_wp, "the clean-up did not run: the aggregates are untouched")
      if (allocated(error)) return
      call check(error, st%fr(1, 1) == 0.0_wp .and. st%fh(1, 1) == 0.0_wp, "the cell got no flux")
      if (allocated(error)) return
      call check(error, iand(st%ierr(2, 1), ICE_ERR_FATAL_MASK) == 0_int32 .and. st%vice(2, 1) /= -777.0_wp, &
                 "the healthy neighbour is cleaned")
   end subroutine test_midway_fatal

end module test_cice_therm2
