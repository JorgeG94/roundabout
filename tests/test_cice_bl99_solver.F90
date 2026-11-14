module test_cice_bl99_solver
   !! Tests of the BL99 temperature solve `ice_temperature_changes`.  Expected
   !! values are closed forms or independent evaluations written here, never the
   !! port's own output; the comparison with the upstream routine is in
   !! `test_cice_golden_bl99_solver`.
   !!
   !! A4  energy closure.  Over one step the change of the column enthalpy
   !!     (recomputed here from the new temperatures with the BL99 enthalpy
   !!     formulas) equals the net flux into the column times `dt`:
   !!     `fcondtopn - fcondbot + fswint`.  Fresh ice with the linear test
   !!     surface closes to round-off (measured, tolerance below); briny ice
   !!     (the production case) closes only to the iteration bound
   !!     `0.9*ferrmax = 9e-4 W/m2`, because the scheme accepts a pass with that
   !!     much energy error (condition 5), not because of round-off.  1, 2, 4, 7
   !!     and 8 ice layers, with and without snow, cold and melting surfaces.
   !! A5  melting surface (`Tsf = 0` exactly, surplus flux, layers capped at
   !!     their melting point), the branch switches in both directions, the
   !!     flux corrections (the corrected parts sum to the corrected total),
   !!     `sw_redist` (closed form, conservation, no-op when off, bounded by
   !!     `fswint`), the error bits (`ICE_ERR_TSFC_NOCONV` from a capped
   !!     iteration, `ICE_ERR_NONFINITE` with every in/out argument bitwise
   !!     untouched), and the steady linear profile (A3 through the solver).
   !! A7 (flux half): `fcondbot = k*dT/h` for a column in steady conduction.
   !! Device: the solver called from a `do concurrent` over mapped arrays.
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_TSFC_NOCONV, ICE_ERR_TIN_RANGE
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71
   use rdb_ice_thermo_shared, only: ICE_SFC_BULK, ICE_SFC_LINEAR
   use rdb_ice_bl99, only: ICE_NITERMAX, ICE_FERRMAX, ICE_SW_FRAC, ICE_SW_DTEMP, ice_temperature_changes, &
                           ice_temperature_changes_cell
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_bl99_solver_tests

   real(wp), parameter :: TOL_CLOSE_FRESH = 1.0e-9_wp
      !! Fresh-ice closure [W/m2]: round-off of enthalpies of order 1e9 J/m3
      !! times a 1 m layer divided by dt (measured worst below).
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: TOL_DEVICE = 1.0e-11_wp
      !! GPU build (FMA on for host code too): measured worst 1.7e-13 (V100, nvfortran 26.5, 512 x 512 grid) against the serial host sweep
#else
   real(wp), parameter :: TOL_DEVICE = 0.0_wp
#endif
   real(wp), parameter :: ERRMAX_BRINE = 0.9_wp*ICE_FERRMAX
      !! The brine closure bound: the iteration accepts `ferr <= 0.9*ferrmax`.

   type :: col_t
      !! One column and the atmosphere above it (host test data only).
      integer :: nil = 4
      real(wp) :: hilyr = 0.0_wp, hslyr = 0.0_wp
      real(wp) :: S(8) = 0.0_wp, T(8) = 0.0_wp, q(8) = 0.0_wp, Isw(8) = 0.0_wp
      real(wp) :: Tsn = 0.0_wp, qsn = 0.0_wp, Ssw = 0.0_wp
      real(wp) :: fswsfc = 0.0_wp, fswint = 0.0_wp, Tsf = 0.0_wp, Tbot = -1.8_wp
      real(wp) :: dt = 3600.0_wp, rhoa = 1.3_wp, flw = 200.0_wp, potT = 258.15_wp, Qa = 5.0e-4_wp
      real(wp) :: shcoef = 8.0_wp, lhcoef = 2.0e4_wp, emis = 0.985_wp, sf0 = 0.0_wp, dsf = 0.0_wp
      integer :: sfc_kind = ICE_SFC_BULK, conduct = ICE_CONDUCT_BUBBLY, nitermax = ICE_NITERMAX
      logical :: brine = .true., swr = .false.
      real(wp) :: sw_frac = ICE_SW_FRAC, sw_dtemp = ICE_SW_DTEMP
      ! results
      real(wp) :: fsens = 0.0_wp, flat = 0.0_wp, flwout = 0.0_wp, fsurf = 0.0_wp, fct = 0.0_wp, fcb = 0.0_wp
      real(wp) :: einex = 0.0_wp, diag(3) = 0.0_wp, einit = 0.0_wp
      integer(int32) :: ierr = 0_int32
   end type col_t

   type :: grid_t
      !! Grid-shaped state of the device test (host container; only its
      !! component arrays are ever mapped).
      integer :: nx = 0, ny = 0, nc = 0, nl = 0
      real(wp), allocatable :: potT(:, :), flw(:, :), Qa(:, :)
      real(wp), allocatable :: hilyr(:, :, :), hslyr(:, :, :), einit(:, :, :), fsw(:, :, :), fsi(:, :, :), Tsf(:, :, :)
      real(wp), allocatable :: fsens(:, :, :), flat(:, :, :), flwout(:, :, :), fsurf(:, :, :), fct(:, :, :), fcb(:, :, :)
      real(wp), allocatable :: einex(:, :, :)
      integer(int32), allocatable :: niter(:, :, :)
      real(wp), allocatable :: S(:, :, :, :), Isw(:, :, :, :), q(:, :, :, :), T(:, :, :, :)
      real(wp), allocatable :: Ssw(:, :, :, :), qsn(:, :, :, :), Tsn(:, :, :, :)
      integer(int32), allocatable :: ierr(:, :, :)
   end type grid_t

contains

   subroutine collect_cice_bl99_solver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("energy_closure_fresh_round_off", test_closure_fresh), &
                  new_unittest("energy_closure_brine_bound", test_closure_brine), &
                  new_unittest("steady_profile_ice_only", test_steady_ice), &
                  new_unittest("steady_profile_snow_ice", test_steady_snow), &
                  new_unittest("conductive_flux_half_of_stefan", test_stefan_flux), &
                  new_unittest("melting_surface_pins_zero", test_melting), &
                  new_unittest("branch_switch_cold_to_melting", test_switch_to_melt), &
                  new_unittest("branch_switch_melting_to_cold", test_switch_to_cold), &
                  new_unittest("flux_corrections_sum_to_total", test_flux_corrections), &
                  new_unittest("sw_redist_closed_form", test_swredist), &
                  new_unittest("sw_redist_off_is_noop", test_swredist_off), &
                  new_unittest("sw_redist_bounded_by_fswint", test_swredist_bound), &
                  new_unittest("no_convergence_bit", test_noconv), &
                  new_unittest("nonfinite_flags_and_leaves_arguments", test_nonfinite), &
                  new_unittest("device_matches_host", test_device) &
                  ]
   end subroutine collect_cice_bl99_solver_tests

   ! ------------------------------------------------------------------
   ! helpers: independent enthalpies, column set-up, the call
   ! ------------------------------------------------------------------

   pure function q_ice(T, S, c, brine) result(q)
      !! BL99 ice enthalpy [J/m3], written out independently here.
      real(wp), intent(in) :: T, S
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: brine
      real(wp) :: q, Tm
      if (brine) then
         Tm = -c%depressT*S
         q = -c%rhoi*(c%cp_ice*(Tm - T) + c%Lfresh*(1.0_wp - Tm/T) - c%cp_ocn*Tm)
      else
         q = c%rhoi*(c%cp_ice*T - c%Lfresh)
      end if
   end function q_ice

   pure function q_snow(T, c) result(q)
      real(wp), intent(in) :: T
      type(ice_constants_t), intent(in) :: c
      real(wp) :: q
      q = -c%rhos*(c%Lfresh - c%cp_ice*T)
   end function q_snow

   pure function col_energy(nil, hilyr, hslyr, S, T, Tsn, c, brine) result(e)
      !! Column enthalpy per area [J/m2].
      integer, intent(in) :: nil
      real(wp), intent(in) :: hilyr, hslyr, S(8), T(8), Tsn
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: brine
      real(wp) :: e
      integer :: k
      e = hslyr*q_snow(Tsn, c)
      do k = 1, nil
         e = e + hilyr*q_ice(T(k), S(k), c, brine)
      end do
   end function col_energy

   function rnd(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function rnd

   subroutine set_q(col, c)
      !! Fill the enthalpies and the initial energy from `T`, `S`, `Tsn`.
      type(col_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      integer :: k
      do k = 1, col%nil
         col%q(k) = q_ice(col%T(k), col%S(k), c, col%brine)
      end do
      col%qsn = q_snow(col%Tsn, c)
      col%einit = col_energy(col%nil, col%hilyr, col%hslyr, col%S, col%T, col%Tsn, c, col%brine)
   end subroutine set_q

   subroutine run_col(col, c)
      !! Call the solver on `col` (in/out arguments updated in place).
      type(col_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp) :: Ssw(1), zqs(1), zTs(1)
      Ssw(1) = col%Ssw
      zqs(1) = col%qsn
      zTs(1) = col%Tsn
      col%ierr = 0_int32
      call ice_temperature_changes(col%nil, 1, col%dt, col%rhoa, col%flw, col%potT, col%Qa, col%shcoef, col%lhcoef, &
                                   col%emis, col%sfc_kind, col%sf0, col%dsf, c, col%conduct, col%brine, &
                                   col%swr, col%sw_frac, col%sw_dtemp, col%nitermax, &
                                   col%hilyr, col%hslyr, col%S, col%Tbot, col%einit, &
                                   col%fswsfc, col%fswint, Ssw, col%Isw, col%q, col%T, zqs, zTs, col%Tsf, &
                                   col%fsens, col%flat, col%flwout, col%fsurf, col%fct, col%fcb, col%einex, col%diag, &
                                   col%ierr)
      col%Ssw = Ssw(1)
      col%qsn = zqs(1)
      col%Tsn = zTs(1)
   end subroutine run_col

   subroutine random_col(st, nil, snow, brine, col, c)
      !! A random sensible column: ice 0.2-3.2 m, optional snow 0.02-0.4 m,
      !! temperatures linear between a skin value and the freezing point and
      !! below the layer melting points, modest absorbed shortwave.
      integer(int64), intent(inout) :: st
      integer, intent(in) :: nil
      logical, intent(in) :: snow, brine
      type(col_t), intent(out) :: col
      type(ice_constants_t), intent(in) :: c
      integer :: k
      real(wp) :: hin, skin, zn, Tm, w

      col%nil = nil
      col%brine = brine
      hin = 0.2_wp + 3.0_wp*rnd(st)
      col%hilyr = hin/real(nil, wp)
      col%hslyr = 0.0_wp
      if (snow) col%hslyr = 0.02_wp + 0.38_wp*rnd(st)
      skin = -3.0_wp - 30.0_wp*rnd(st)
      do k = 1, nil
         zn = (real(k, wp) - 0.5_wp)/real(nil, wp)
         col%S(k) = merge(0.1_wp + 4.5_wp*rnd(st), 0.0_wp, brine)
         Tm = merge(-c%depressT*col%S(k), 0.0_wp, brine)
         col%T(k) = min(skin + (col%Tbot - skin)*zn, Tm - 0.05_wp)
      end do
      col%Tsn = merge(0.5_wp*skin, 0.0_wp, snow)
      col%Tsf = skin
      w = 5.0_wp + 40.0_wp*rnd(st)
      col%fswint = 0.0_wp
      do k = 1, nil
         col%Isw(k) = w*exp(-0.5_wp*real(k - 1, wp))/real(nil, wp)
         col%fswint = col%fswint + col%Isw(k)
      end do
      col%Ssw = 0.0_wp
      col%fswsfc = 10.0_wp*rnd(st)
      col%dt = merge(3600.0_wp, 1800.0_wp, rnd(st) > 0.5_wp)
      col%potT = 273.15_wp + skin - 5.0_wp*rnd(st)
      col%flw = 140.0_wp + 80.0_wp*rnd(st)
      col%Qa = 4.0e-4_wp + 6.0e-4_wp*rnd(st)
      col%conduct = merge(ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71, rnd(st) > 0.5_wp)
      call set_q(col, c)
   end subroutine random_col

   function net_flux_in(col) result(f)
      !! Net flux into the column (top conduction minus bottom conduction plus
      !! the interior shortwave) [W/m2].
      type(col_t), intent(in) :: col
      real(wp) :: f
      f = col%fct - col%fcb + col%fswint
   end function net_flux_in

   function closure_error(col, c) result(r)
      !! |d(column enthalpy)/dt - net flux in| [W/m2] for a solved column
      !! (`col%einit` is the start energy).
      type(col_t), intent(in) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp) :: r, eend
      eend = col_energy(col%nil, col%hilyr, col%hslyr, col%S, col%T, col%Tsn, c, col%brine)
      r = abs((eend - col%einit)/col%dt - net_flux_in(col))
   end function closure_error

   ! ------------------------------------------------------------------
   ! A4: energy closure
   ! ------------------------------------------------------------------

   subroutine closure_sweep(error, label, brine, linear, tol, worst)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: label
      logical, intent(in) :: brine, linear
      real(wp), intent(in) :: tol
      real(wp), intent(out) :: worst
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer, parameter :: NILS(5) = [1, 2, 4, 7, 8]
      integer(int64) :: st
      integer :: i, j, nil, nsnow, nmelt
      real(wp) :: r
      character(len=160) :: msg

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      st = 987654321_int64
      worst = 0.0_wp
      nmelt = 0
      do i = 1, 5
         nil = NILS(i)
         do j = 1, 24
            nsnow = mod(j, 2)
            call random_col(st, nil, nsnow == 1, brine, col, c)
            if (linear) then
               ! linear surface: heating or cooling of a few tens of W/m2, so both branches occur
               col%sfc_kind = ICE_SFC_LINEAR
               col%dsf = 8.0_wp
               col%sf0 = -60.0_wp + 120.0_wp*rnd(st) - 8.0_wp*col%Tsf
               if (mod(j, 4) == 0) col%sf0 = -200.0_wp*rnd(st)
            else
               if (mod(j, 6) == 0) then
                  col%potT = 274.0_wp + 4.0_wp*rnd(st)      ! warm air: melting surface
                  col%flw = 300.0_wp + 40.0_wp*rnd(st)
                  col%Qa = 4.0e-3_wp
                  col%fswsfc = 50.0_wp + 100.0_wp*rnd(st)
                  col%Tsf = 0.0_wp
               end if
            end if
            call run_col(col, c)
            write (msg, '(a,a,i0,a,i0,a,i0)') label, ": nilyr ", nil, " case ", j, " error bits ", col%ierr
            if (col%ierr /= 0_int32) then
               call test_failed(error, trim(msg))
               return
            end if
            if (.not. (ieee_is_finite(col%Tsf) .and. ieee_is_finite(col%fct) .and. ieee_is_finite(col%fcb) &
                       .and. all(ieee_is_finite(col%T(1:nil))))) then
               call test_failed(error, trim(msg)//" non-finite output")
               return
            end if
            if (col%Tsf == 0.0_wp) nmelt = nmelt + 1
            r = closure_error(col, c)
            if (.not. ieee_is_finite(r)) then
               call test_failed(error, trim(msg)//" non-finite closure")
               return
            end if
            if (r > worst) worst = r
            if (.not. (r <= tol)) then
               write (msg, '(a,a,i0,a,i0,a,es10.3,a,es10.3)') label, ": nilyr ", nil, " case ", j, &
                  " closure error ", r, " > ", tol
               call test_failed(error, trim(msg))
               return
            end if
         end do
      end do
      call check(error, nmelt > 0, label//": the sweep contains melting-surface columns")
   end subroutine closure_sweep

   subroutine test_closure_fresh(error)
      !! Fresh ice, linear surface: the converged solution satisfies the
      !! discrete energy equation, so the closure is round-off (measured
      !! worst 1.4e-10 W/m2 on gfortran, enthalpies of order 1e9 J/m3).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: worst
      call closure_sweep(error, "fresh/linear", .false., .true., TOL_CLOSE_FRESH, worst)
   end subroutine test_closure_fresh

   subroutine test_closure_brine(error)
      !! Briny ice, bulk surface (the production case): the iteration
      !! accepts a pass whose energy error is below `0.9*ferrmax`, and does
      !! not do better than round-off.  The assertion is that bound, stated
      !! here and not a round-off claim.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: worst
      call closure_sweep(error, "brine/bulk", .true., .false., ERRMAX_BRINE, worst)
   end subroutine test_closure_brine

   ! ------------------------------------------------------------------
   ! A3 through the solver, A7 flux half
   ! ------------------------------------------------------------------

   subroutine steady_march(col, c, nstep)
      !! Step the column to a steady state with the linear surface: `dt` far
      !! beyond the diffusion time (backward Euler is unconditionally stable).
      type(col_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nstep
      integer :: n
      do n = 1, nstep
         call set_q(col, c)
         call run_col(col, c)
         if (col%ierr /= 0_int32) return
      end do
   end subroutine steady_march

   subroutine test_steady_ice(error)
      !! Ice only, fresh, constant `k = kice` (MU71, S = 0), linear surface
      !! `fsurf = lam*(Ta - Tsf)`: the steady state is the straight line with
      !! `Tsf = (lam*Ta + K*Tbot)/(lam + K)`, `K = kice/hi`, and the
      !! conductive flux `K*(Tsf - Tbot)` at the top and the bottom.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer, parameter :: NILS(4) = [1, 2, 4, 8]
      integer :: i, k
      real(wp) :: lam, Ta, K_, Tsf_e, hi, zn, worst

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      lam = 12.0_wp
      Ta = -25.0_wp
      worst = 0.0_wp
      do i = 1, 4
         col = col_t()
         col%nil = NILS(i)
         col%brine = .false.
         col%conduct = ICE_CONDUCT_MU71
         col%sfc_kind = ICE_SFC_LINEAR
         hi = 1.5_wp
         col%hilyr = hi/real(col%nil, wp)
         col%dt = 1.0e6_wp
         col%dsf = lam
         col%sf0 = -lam*Ta
         col%T(1:col%nil) = -5.0_wp
         col%Tsf = -5.0_wp
         call steady_march(col, c, 30)
         call check(error, col%ierr == 0_int32, "steady ice: no error bits")
         if (allocated(error)) return
         K_ = c%kice/hi
         Tsf_e = (lam*Ta + K_*col%Tbot)/(lam + K_)
         worst = max(worst, abs(col%Tsf - Tsf_e))
         do k = 1, col%nil
            zn = (real(k, wp) - 0.5_wp)/real(col%nil, wp)
            if (.not. (abs(col%T(k) - (Tsf_e + (col%Tbot - Tsf_e)*zn)) <= 1.0e-9_wp)) then
               call test_failed(error, "steady ice: profile is not the straight line")
               return
            end if
         end do
         call check(error, abs(col%Tsf - Tsf_e) <= 1.0e-9_wp, "steady ice: surface temperature")
         if (allocated(error)) return
         call check(error, abs(col%fct - K_*(Tsf_e - col%Tbot)) <= 1.0e-8_wp, "steady ice: top conductive flux")
         if (allocated(error)) return
         call check(error, abs(col%fcb - K_*(Tsf_e - col%Tbot)) <= 1.0e-8_wp, "steady ice: bottom conductive flux")
         if (allocated(error)) return
      end do
   end subroutine test_steady_ice

   subroutine test_steady_snow(error)
      !! Snow over ice, fresh: `K = 1/(hs/ksno + hi/kice)`, the interface
      !! temperature `T_i = Tbot + (Tsf - Tbot)*(hi/kice)/(hs/ksno + hi/kice)`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      real(wp) :: lam, Ta, K_, Tsf_e, hi, hs, Ti_e, Tsn_e

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      lam = 12.0_wp
      Ta = -30.0_wp
      hi = 2.0_wp
      hs = 0.3_wp
      col = col_t()
      col%nil = 4
      col%brine = .false.
      col%conduct = ICE_CONDUCT_MU71
      col%sfc_kind = ICE_SFC_LINEAR
      col%hilyr = hi/4.0_wp
      col%hslyr = hs
      col%dt = 1.0e6_wp
      col%dsf = lam
      col%sf0 = -lam*Ta
      col%T(1:4) = -6.0_wp
      col%Tsn = -8.0_wp
      col%Tsf = -10.0_wp
      call steady_march(col, c, 30)
      call check(error, col%ierr == 0_int32, "steady snow: no error bits")
      if (allocated(error)) return
      K_ = 1.0_wp/(hs/c%ksno + hi/c%kice)
      Tsf_e = (lam*Ta + K_*col%Tbot)/(lam + K_)
      call check(error, abs(col%Tsf - Tsf_e) <= 1.0e-9_wp, "steady snow: surface temperature")
      if (allocated(error)) return
      call check(error, abs(col%fct - K_*(Tsf_e - col%Tbot)) <= 1.0e-8_wp, "steady snow: top flux")
      if (allocated(error)) return
      call check(error, abs(col%fcb - K_*(Tsf_e - col%Tbot)) <= 1.0e-8_wp, "steady snow: bottom flux")
      if (allocated(error)) return
      ! the snow layer sits at the middle of the snow: temperature at hs/2 below the surface
      Tsn_e = Tsf_e - K_*(Tsf_e - col%Tbot)*(0.5_wp*hs/c%ksno)
      Ti_e = Tsf_e - K_*(Tsf_e - col%Tbot)*(hs/c%ksno)
      call check(error, abs(col%Tsn - Tsn_e) <= 1.0e-8_wp, "steady snow: snow temperature")
      if (allocated(error)) return
      ! top ice layer midpoint is hi/8 below the interface
      call check(error, abs(col%T(1) - (Ti_e - K_*(Tsf_e - col%Tbot)*(hi/8.0_wp/c%kice))) <= 1.0e-8_wp, &
                 "steady snow: top ice layer temperature")
   end subroutine test_steady_snow

   subroutine test_stefan_flux(error)
      !! A7 conductive half: a column whose surface is held at the (cold) air
      !! temperature by a very large exchange coefficient conducts
      !! `fcondbot = kice*(Tsf - Tbot)/h`, the flux that sets the Stefan
      !! growth rate `dh/dt = fcondbot/(rhoi*Lfresh)` (Stefan 1891; Maykut and
      !! Untersteiner 1971), for a range of thicknesses and air temperatures.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer :: i, j
      real(wp) :: hs_(4), Ta_(3), lam, K_, dhdt, dhdt_e, Tsf_e

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      hs_ = [0.1_wp, 0.5_wp, 1.5_wp, 3.0_wp]
      Ta_ = [-8.0_wp, -20.0_wp, -35.0_wp]
      lam = 1.0e7_wp
      do i = 1, 4
         do j = 1, 3
            col = col_t()
            col%nil = 4
            col%brine = .false.
            col%conduct = ICE_CONDUCT_MU71
            col%sfc_kind = ICE_SFC_LINEAR
            col%hilyr = hs_(i)/4.0_wp
            col%dt = 1.0e6_wp
            col%dsf = lam
            col%sf0 = -lam*Ta_(j)
            col%T(1:4) = Ta_(j)
            col%Tsf = Ta_(j)
            call steady_march(col, c, 40)
            call check(error, col%ierr == 0_int32, "stefan flux: no error bits")
            if (allocated(error)) return
            K_ = c%kice/hs_(i)
            call check(error, abs(col%fcb - K_*(Ta_(j) - col%Tbot)) <= 1.0e-5_wp*abs(col%fcb), &
                       "stefan flux: fcondbot = k dT / h (surface held within 1e-5 K of the air)")
            if (allocated(error)) return
            ! exact steady state of the finite exchange coefficient
            Tsf_e = (lam*Ta_(j) + K_*col%Tbot)/(lam + K_)
            call check(error, abs(col%fcb - K_*(Tsf_e - col%Tbot)) <= 1.0e-8_wp*abs(col%fcb), &
                       "stefan flux: exact steady flux")
            if (allocated(error)) return
            dhdt = col%fcb/(c%rhoi*c%Lfresh)
            dhdt_e = c%kice*(Ta_(j) - col%Tbot)/(hs_(i)*c%rhoi*c%Lfresh)
            call check(error, abs(dhdt - dhdt_e) <= 1.0e-5_wp*abs(dhdt_e), "stefan flux: growth rate")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_stefan_flux

   ! ------------------------------------------------------------------
   ! A5: melting surface and branch switching
   ! ------------------------------------------------------------------

   function warm_col(c, nil, snow) result(col)
      !! A column heated at the surface: linear surface with `100 W/m2` at 0 degC.
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nil
      logical, intent(in) :: snow
      type(col_t) :: col
      integer :: k
      col%nil = nil
      col%brine = .true.
      col%sfc_kind = ICE_SFC_LINEAR
      col%hilyr = 1.2_wp/real(nil, wp)
      col%hslyr = merge(0.1_wp, 0.0_wp, snow)
      col%dsf = 10.0_wp
      col%sf0 = -100.0_wp
      do k = 1, nil
         col%S(k) = 3.0_wp
         col%T(k) = -4.0_wp + 2.0_wp*real(k - 1, wp)/real(nil, wp)
      end do
      col%Tsn = merge(-5.0_wp, 0.0_wp, snow)
      col%Tsf = -2.0_wp
      col%dt = 3600.0_wp
      call set_q(col, c)
   end function warm_col

   subroutine test_melting(error)
      !! Heating through the surface: `Tsf = 0` exactly, the surface flux
      !! exceeds the conduction (the surplus melts, booked by the caller), the
      !! layers stay at or below their melting points, and the closure holds.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer :: i, k
      integer, parameter :: NILS(3) = [1, 4, 7]

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do i = 1, 3
         col = warm_col(c, NILS(i), mod(i, 2) == 1)
         call run_col(col, c)
         call check(error, col%ierr == 0_int32, "melting: no error bits")
         if (allocated(error)) return
         call check(error, col%Tsf == 0.0_wp, "melting: Tsf is exactly 0")
         if (allocated(error)) return
         call check(error, col%fsurf > col%fct, "melting: surface flux exceeds the conduction")
         if (allocated(error)) return
         ! the linear surface: fsurf = -(sf0 + dsf*Tsf) = 100 at Tsf = 0
         call check(error, abs(col%fsurf - 100.0_wp) <= 1.0e-12_wp, "melting: fsurf = 100 W/m2")
         if (allocated(error)) return
         do k = 1, NILS(i)
            if (.not. (col%T(k) <= -c%depressT*col%S(k) + 1.0e-12_wp)) then
               call test_failed(error, "melting: layer above its melting point")
               return
            end if
         end do
         call check(error, col%einex == 0.0_wp, "melting: einex_sfc_flux is 0 for calc_Tsfc")
         if (allocated(error)) return
         call check(error, closure_error(col, c) <= ERRMAX_BRINE, "melting: closure within the bound")
         if (allocated(error)) return
      end do
   end subroutine test_melting

   subroutine test_switch_to_melt(error)
      !! Start cold (`Tsf = -2`) under strong heating: the first pass
      !! solves the cold branch, finds `Tsf > puny`, resets it to 0 and
      !! switches to the melting branch, so more than one pass is needed.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      col = warm_col(c, 4, .false.)
      col%sf0 = -300.0_wp
      col%Tsf = -2.0_wp
      call run_col(col, c)
      call check(error, col%ierr == 0_int32, "switch to melting: no error bits")
      if (allocated(error)) return
      call check(error, col%Tsf == 0.0_wp, "switch to melting: ends at Tsf = 0")
      if (allocated(error)) return
      call check(error, col%fsurf > col%fct, "switch to melting: surplus surface flux")
   end subroutine test_switch_to_melt

   subroutine test_switch_to_cold(error)
      !! Start at `Tsf = 0` under a surface that cools (`fsurf < fcondtopn`):
      !! the pass resets `Tsf = -puny`, takes the cold branch and the
      !! converged state balances `fsurfn = fcondtopn` (linear surface: exact
      !! to the iteration bound) with `Tsf < 0`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      col = warm_col(c, 4, .true.)
      col%sf0 = 60.0_wp       ! fsurf = -60 at 0 degC: a net loss
      col%Tsf = 0.0_wp
      call run_col(col, c)
      call check(error, col%ierr == 0_int32, "switch to cold: no error bits")
      if (allocated(error)) return
      call check(error, col%Tsf < 0.0_wp, "switch to cold: ends below 0")
      if (allocated(error)) return
      call check(error, abs(col%fsurf - col%fct) <= 1.0e-9_wp, "switch to cold: surface flux balances the conduction")
   end subroutine test_switch_to_cold

   subroutine test_flux_corrections(error)
      !! Bulk surface: after convergence the corrected parts (outgoing
      !! longwave, sensible, latent, each advanced by `dTsf_prev` times its
      !! derivative) add to the corrected total: `fswsfc + emis*flw + flwout +
      !! fsens + flat = fsurfn` to round-off, and `fsurfn = fcondtopn` on the
      !! cold branch to the iteration bound.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer(int64) :: st
      integer :: j
      real(wp) :: sum_parts

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      st = 4242_int64
      do j = 1, 20
         call random_col(st, 4, mod(j, 2) == 0, .true., col, c)
         call run_col(col, c)
         call check(error, col%ierr == 0_int32, "flux corrections: no error bits")
         if (allocated(error)) return
         sum_parts = col%fswsfc + col%emis*col%flw + col%flwout + col%fsens + col%flat
         call check(error, abs(sum_parts - col%fsurf) <= 1.0e-10_wp*max(1.0_wp, abs(col%fsurf)), &
                    "flux corrections: parts add to the corrected total")
         if (allocated(error)) return
         if (col%Tsf < 0.0_wp) then
            call check(error, abs(col%fsurf - col%fct) <= 1.0e-2_wp, "flux corrections: cold surface balance")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_flux_corrections

   ! ------------------------------------------------------------------
   ! sw_redist
   ! ------------------------------------------------------------------

   function swr_col(c, brine) result(col)
      !! Top ice layer 0.01 K below the redistribution band edge region: layer
      !! 1 inside the band (all its absorbed shortwave moves up), layer 2 below
      !! the band with more shortwave than `sw_frac` of its warm-up energy,
      !! layer 3 cold enough that nothing moves; warm snow inside the cap.
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: brine
      type(col_t) :: col
      integer :: k
      col%nil = 3
      col%brine = brine
      col%sfc_kind = ICE_SFC_LINEAR
      col%hilyr = 1.0e-3_wp
      col%hslyr = 0.05_wp
      col%dsf = 10.0_wp
      col%sf0 = 10.0_wp
      col%dt = 3600.0_wp
      col%swr = .true.
      do k = 1, 3
         col%S(k) = merge(2.0_wp, 0.0_wp, brine)
      end do
      col%T(1) = merge(-0.108_wp + 0.005_wp, -0.005_wp, brine)   ! inside the band (within sw_dtemp of melting)
      col%T(2) = merge(-0.108_wp - 0.05_wp, -0.05_wp, brine)     ! just below the band
      col%T(3) = -8.0_wp
      col%Tsn = -0.05_wp
      col%Isw(1) = 30.0_wp
      col%Isw(2) = 50.0_wp
      col%Isw(3) = 1.0_wp
      col%Ssw = 4.0_wp
      col%fswint = col%Isw(1) + col%Isw(2) + col%Isw(3) + col%Ssw
      col%fswsfc = 5.0_wp
      col%Tsf = -2.0_wp
      call set_q(col, c)
   end function swr_col

   subroutine test_swredist(error)
      !! Closed forms.  Layer in the band: all of its shortwave moves to the
      !! surface.  Layer below the band: kept up to `sw_frac*(Tmlt - T)*ci/
      !! (dt/(rhoi*hilyr))` with `ci = cp_ice - Lfresh*Tmlt/T^2` (brine) or
      !! `cp_ice`, the rest moves.  Snow: `-sw_frac*T/etas`.  A cold layer
      !! keeps everything.  `fswsfc + fswint` is conserved to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      integer :: ib
      logical :: brine
      real(wp) :: tot0, Tm, ci, cap2, etas, caps, f1, moved
      real(wp) :: Isw0(3), Ssw0, fsw0, fsi0, T2

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do ib = 1, 2
         brine = ib == 1
         col = swr_col(c, brine)
         Isw0 = col%Isw(1:3)
         Ssw0 = col%Ssw
         fsw0 = col%fswsfc
         fsi0 = col%fswint
         tot0 = fsw0 + fsi0
         Tm = merge(-c%depressT*2.0_wp, 0.0_wp, brine)
         T2 = col%T(2)
         ! closed forms from the START temperatures
         if (brine) then
            ci = c%cp_ice - c%Lfresh*Tm/(T2*T2)
            cap2 = ICE_SW_FRAC*(Tm - T2)*ci/(col%dt/(c%rhoi*col%hilyr))
         else
            cap2 = ICE_SW_FRAC*(-T2)*c%cp_ice/(col%dt/(c%rhoi*col%hilyr))
         end if
         etas = col%dt/(c%rhos*c%cp_ice*col%hslyr)
         caps = -ICE_SW_FRAC*col%Tsn/etas
         call run_col(col, c)
         call check(error, col%ierr == 0_int32, "sw_redist: no error bits")
         if (allocated(error)) return
         call check(error, col%Isw(1) == 0.0_wp, "sw_redist: layer in the band keeps nothing")
         if (allocated(error)) return
         call check(error, abs(col%Isw(2) - cap2) <= 1.0e-12_wp*cap2, "sw_redist: layer below the band keeps the cap")
         if (allocated(error)) return
         call check(error, col%Isw(3) == Isw0(3), "sw_redist: cold layer unchanged")
         if (allocated(error)) return
         call check(error, abs(col%Ssw - caps) <= 1.0e-12_wp*caps, "sw_redist: snow keeps its cap")
         if (allocated(error)) return
         moved = (Isw0(1) - col%Isw(1)) + (Isw0(2) - col%Isw(2)) + (Ssw0 - col%Ssw)
         f1 = col%fswsfc - fsw0
         call check(error, abs(f1 - moved) <= 1.0e-9_wp, "sw_redist: the surface receives what the layers lose")
         if (allocated(error)) return
         call check(error, abs((col%fswsfc + col%fswint) - tot0) <= 1.0e-9_wp, "sw_redist: fswsfc + fswint conserved")
         if (allocated(error)) return
         call check(error, abs(col%fswint - (col%Isw(1) + col%Isw(2) + col%Isw(3) + col%Ssw)) <= 1.0e-9_wp, &
                    "sw_redist: fswint equals what the layers now absorb")
         if (allocated(error)) return
      end do
   end subroutine test_swredist

   subroutine test_swredist_off(error)
      !! With the switch off the shortwave bookkeeping is untouched, bitwise.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col
      real(wp) :: Isw0(3), Ssw0, fsw0, fsi0

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      col = swr_col(c, .true.)
      col%swr = .false.
      Isw0 = col%Isw(1:3)
      Ssw0 = col%Ssw
      fsw0 = col%fswsfc
      fsi0 = col%fswint
      call run_col(col, c)
      call check(error, col%ierr == 0_int32, "sw_redist off: no error bits")
      if (allocated(error)) return
      call check(error, all(col%Isw(1:3) == Isw0) .and. col%Ssw == Ssw0 .and. col%fswsfc == fsw0 &
                 .and. col%fswint == fsi0, "sw_redist off: nothing moved")
   end subroutine test_swredist_off

   subroutine test_swredist_bound(error)
      !! The amount moved is bounded by the interior shortwave left
      !! (`dswabs = min(excess, fswint)`): with a deliberately small `fswint`
      !! the surface cannot receive more than it holds, and `fswint` never
      !! goes negative (upstream drops the part it cannot move; kept).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      col = swr_col(c, .true.)
      col%fswint = 10.0_wp
      ! the bookkeeping is inconsistent on purpose (the layers absorb more than fswint, and
      ! upstream drops the part it cannot move), so the heat balance cannot close: one pass only
      col%nitermax = 1
      call run_col(col, c)
      call check(error, col%ierr == ICE_ERR_TSFC_NOCONV, "sw_redist bound: only the capped-iteration bit")
      if (allocated(error)) return
      call check(error, col%fswint >= 0.0_wp, "sw_redist bound: fswint stays non-negative")
      if (allocated(error)) return
      call check(error, abs(col%fswsfc - 15.0_wp) <= 1.0e-12_wp, "sw_redist bound: the surface got exactly fswint")
   end subroutine test_swredist_bound

   ! ------------------------------------------------------------------
   ! error bits
   ! ------------------------------------------------------------------

   subroutine test_noconv(error)
      !! A capped iteration raises `ICE_ERR_TSFC_NOCONV` (and only it), the
      !! outputs are the last iterate (finite), `diag(1)` is the cap; the
      !! same column converges with the full cap.  Earlier bits are kept.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(col_t) :: col, ref
      integer(int64) :: st

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      st = 31_int64
      call random_col(st, 4, .true., .true., col, c)
      ref = col
      call run_col(ref, c)
      call check(error, ref%ierr == 0_int32 .and. ref%diag(1) > 1.0_wp, "noconv: reference converges in more than one pass")
      if (allocated(error)) return
      col%nitermax = 1
      call run_col(col, c)
      call check(error, col%ierr == ICE_ERR_TSFC_NOCONV, "noconv: exactly the NOCONV bit")
      if (allocated(error)) return
      call check(error, nint(col%diag(1)) == 1, "noconv: one pass done")
      if (allocated(error)) return
      call check(error, ieee_is_finite(col%Tsf) .and. all(ieee_is_finite(col%T(1:4))) .and. ieee_is_finite(col%fct), &
                 "noconv: outputs are the finite last iterate")
      if (allocated(error)) return
      ! the cap of the iteration count is a real input: 2 passes still short of the converged count
      col = ref
      call random_col(st, 4, .true., .true., col, c)
      col%nitermax = ref%nitermax
      ! earlier bits are kept (OR-ed into)
      call random_col(st, 4, .true., .true., col, c)
      col%nitermax = 1
      col%ierr = ICE_ERR_TIN_RANGE
      call ice_call_keep_bits(col, c)
      call check(error, col%ierr == ior(ICE_ERR_TIN_RANGE, ICE_ERR_TSFC_NOCONV), "noconv: earlier bits kept")
   end subroutine test_noconv

   subroutine ice_call_keep_bits(col, c)
      !! `run_col` resets `ierr`; call the solver once with the bits preset.
      type(col_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp) :: Ssw(1), zqs(1), zTs(1)
      Ssw(1) = col%Ssw
      zqs(1) = col%qsn
      zTs(1) = col%Tsn
      call ice_temperature_changes(col%nil, 1, col%dt, col%rhoa, col%flw, col%potT, col%Qa, col%shcoef, col%lhcoef, &
                                   col%emis, col%sfc_kind, col%sf0, col%dsf, c, col%conduct, col%brine, &
                                   col%swr, col%sw_frac, col%sw_dtemp, col%nitermax, &
                                   col%hilyr, col%hslyr, col%S, col%Tbot, col%einit, &
                                   col%fswsfc, col%fswint, Ssw, col%Isw, col%q, col%T, zqs, zTs, col%Tsf, &
                                   col%fsens, col%flat, col%flwout, col%fsurf, col%fct, col%fcb, col%einex, col%diag, &
                                   col%ierr)
   end subroutine ice_call_keep_bits

   function same_bits(a, b) result(ok)
      !! Bitwise equality of two columns' in/out data (no NaN involved).
      type(col_t), intent(in) :: a, b
      logical :: ok
      ok = all(a%T == b%T) .and. all(a%q == b%q) .and. all(a%Isw == b%Isw) .and. a%Tsn == b%Tsn .and. a%qsn == b%qsn &
           .and. a%Ssw == b%Ssw .and. a%fswsfc == b%fswsfc .and. a%fswint == b%fswint .and. a%Tsf == b%Tsf &
           .and. a%fsens == b%fsens .and. a%flat == b%flat .and. a%flwout == b%flwout .and. a%fsurf == b%fsurf &
           .and. a%fct == b%fct .and. a%fcb == b%fcb .and. a%einex == b%einex .and. all(a%diag == b%diag)
   end function same_bits

   subroutine test_nonfinite(error)
      !! A NaN or Inf in any input, a bad layer count or step, or a surface
      !! temperature at or below absolute zero raises `ICE_ERR_NONFINITE`.
      !! Failures of the inputs are caught before the first write (arguments
      !! bitwise untouched; no NaN is laundered through `min`/comparisons into
      !! a plausible converged state); one found mid-iteration (Tsf below absolute zero, 0/0) only promises the bit.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c, cbad
      type(col_t) :: col, ref
      integer(int64) :: st
      integer :: i
      real(wp) :: nan, inf
      character(len=64) :: what

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      inf = ieee_value(1.0_wp, ieee_positive_inf)
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do i = 1, 24
         st = 77_int64
         call random_col(st, 4, .true., .true., col, c)
         ref = col
         cbad = c
         select case (i)
         case (1); col%rhoa = nan; what = "rhoa"
         case (2); col%flw = nan; what = "flw"
         case (3); col%potT = nan; what = "potT"
         case (4); col%Qa = inf; what = "Qa"
         case (5); col%shcoef = nan; what = "shcoef"
         case (6); col%lhcoef = nan; what = "lhcoef"
         case (7); col%hilyr = nan; what = "hilyr"
         case (8); col%hslyr = nan; what = "hslyr"
         case (9); col%Tbot = nan; what = "Tbot"
         case (10); col%einit = nan; what = "einit"
         case (11); col%fswsfc = nan; what = "fswsfc"
         case (12); col%fswint = nan; what = "fswint"
         case (13); col%Tsf = nan; what = "Tsf"
         case (14); col%T(2) = nan; what = "zTin"
         case (15); col%sf0 = nan; what = "sf_0"
         case (16); col%S(1) = nan; what = "zSin"
         case (17); col%Isw(2) = nan; what = "Iswabs"
         case (18); col%Tsn = nan; what = "zTsn"
         case (19); col%dt = 0.0_wp; what = "dt = 0"
         case (20); col%hilyr = 0.0_wp; what = "hilyr = 0"
         case (21); col%nitermax = 0; what = "nitermax = 0"
         case (22); col%Tsf = -300.0_wp; what = "Tsf below absolute zero"
         case (23); col%rhoa = -inf; what = "rhoa = -Inf"
         case default
            ! zero start temperature with zero salinity: 0/0 in the heat capacity, caught inside the iteration
            col%S(1) = 0.0_wp
            col%T(1) = 0.0_wp
            what = "0/0 heat capacity"
         end select
         ref = col
         col%ierr = 0_int32
         call run_col_c(col, cbad)
         call check(error, iand(col%ierr, ICE_ERR_NONFINITE) /= 0_int32, "nonfinite flagged: "//trim(what))
         if (allocated(error)) return
         ! input failures are caught before any write; a failure found mid-iteration (the last case)
         ! only promises the bit
         if (i < 22 .or. i == 23) then
            call check(error, same_bits_in(col, ref), "nonfinite leaves arguments untouched: "//trim(what))
            if (allocated(error)) return
         end if
      end do
      ! nilyr out of range
      st = 77_int64
      call random_col(st, 4, .true., .true., col, c)
      col%nil = 9
      ref = col
      call run_col(col, c)
      call check(error, iand(col%ierr, ICE_ERR_NONFINITE) /= 0_int32, "nonfinite flagged: nilyr = 9")
      if (allocated(error)) return
      call check(error, same_bits_in(col, ref), "nilyr = 9 leaves arguments untouched")
   end subroutine test_nonfinite

   subroutine run_col_c(col, c)
      !! `run_col` with a possibly corrupted constants bundle.
      type(col_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      call run_col(col, c)
   end subroutine run_col_c

   function same_bits_in(a, b) result(ok)
      !! The in/out data of two columns agree bitwise; NaNs are compared by
      !! their bit pattern.
      type(col_t), intent(in) :: a, b
      logical :: ok
      ok = all(transfer(a%T, 1_int64, 8) == transfer(b%T, 1_int64, 8)) .and. same_real_bits(a%q, b%q) &
           .and. same_real_bits(a%Isw, b%Isw) .and. same_real_bits([a%Tsn, a%qsn, a%Ssw, a%fswsfc, a%fswint, a%Tsf, &
                                                                    a%fsens, a%flat, a%flwout, a%fsurf, a%fct, a%fcb, &
                                                                    a%einex, a%diag(1), a%diag(2), a%diag(3)], &
                                                                   [b%Tsn, b%qsn, b%Ssw, b%fswsfc, b%fswint, b%Tsf, &
                                                                    b%fsens, b%flat, b%flwout, b%fsurf, b%fct, b%fcb, &
                                                                    b%einex, b%diag(1), b%diag(2), b%diag(3)])
   end function same_bits_in

   function same_real_bits(x, y) result(ok)
      real(wp), intent(in) :: x(:), y(:)
      logical :: ok
      integer :: i
      ok = size(x) == size(y)
      if (.not. ok) return
      do i = 1, size(x)
         if (transfer(x(i), 1_int64) /= transfer(y(i), 1_int64)) then
            ok = .false.
            return
         end if
      end do
   end function same_real_bits

   ! ------------------------------------------------------------------
   ! device
   ! ------------------------------------------------------------------

   subroutine init_grid(g, c)
      !! Deterministic columns on the whole grid: ice 0.15-3.15 m over `nl`
      !! layers, bare or snow covered, skin temperatures from -3 to -28 degC,
      !! the atmosphere from -28 to +3 degC at `(i, j)` (so cold and melting
      !! surfaces both occur), every category different.
      type(grid_t), intent(inout) :: g
      type(ice_constants_t), intent(in) :: c
      integer :: i, j, n, k
      real(wp) :: h1, h2, h3, hin, skin, zn, Tm, ph, w, e

      do j = 1, g%ny
         do i = 1, g%nx
            ph = real(mod(i*7 + j*13, 97), wp)/97.0_wp
            g%potT(i, j) = 273.15_wp - 28.0_wp + 31.0_wp*ph
            g%flw(i, j) = 150.0_wp + 190.0_wp*ph
            g%Qa(i, j) = 5.0e-4_wp + 3.5e-3_wp*ph
         end do
      end do
      do n = 1, g%nc
         do j = 1, g%ny
            do i = 1, g%nx
               h1 = real(mod(i*7 + j*13 + n*29, 97), wp)/97.0_wp
               h2 = real(mod(i*11 + j*3 + n*17, 89), wp)/89.0_wp
               h3 = real(mod(i*5 + j*19 + n*7, 83), wp)/83.0_wp
               hin = 0.15_wp + 3.0_wp*h1
               g%hilyr(i, j, n) = hin/real(g%nl, wp)
               g%hslyr(i, j, n) = merge(0.0_wp, 0.02_wp + 0.3_wp*h2, h2 < 0.3_wp)
               skin = -3.0_wp - 25.0_wp*h3
               e = g%hslyr(i, j, n)*q_snow(0.5_wp*skin*merge(0.0_wp, 1.0_wp, h2 < 0.3_wp), c)
               w = 0.0_wp
               do k = 1, g%nl
                  zn = (real(k, wp) - 0.5_wp)/real(g%nl, wp)
                  g%S(i, j, n, k) = 0.5_wp + 2.5_wp*zn
                  Tm = -c%depressT*g%S(i, j, n, k)
                  g%T(i, j, n, k) = min(skin + (-1.8_wp - skin)*zn, Tm - 0.05_wp)
                  g%q(i, j, n, k) = q_ice(g%T(i, j, n, k), g%S(i, j, n, k), c, .true.)
                  g%Isw(i, j, n, k) = 10.0_wp*0.6_wp**(k - 1)*(0.5_wp + h1)
                  w = w + g%Isw(i, j, n, k)
                  e = e + g%hilyr(i, j, n)*g%q(i, j, n, k)
               end do
               g%Tsn(i, j, n, 1) = merge(0.0_wp, 0.5_wp*skin, h2 < 0.3_wp)
               g%qsn(i, j, n, 1) = q_snow(g%Tsn(i, j, n, 1), c)
               g%Ssw(i, j, n, 1) = 0.0_wp
               g%einit(i, j, n) = e
               g%fsw(i, j, n) = 10.0_wp
               g%fsi(i, j, n) = w
               g%Tsf(i, j, n) = skin
            end do
         end do
      end do
      g%fsens = 0.0_wp; g%flat = 0.0_wp; g%flwout = 0.0_wp; g%fsurf = 0.0_wp
      g%fct = 0.0_wp; g%fcb = 0.0_wp; g%einex = 0.0_wp; g%niter = 0_int32
      g%ierr = 0_int32
   end subroutine init_grid

   subroutine sweep_dc(nx, ny, nc, nl, c, potT, flw, Qa, hilyr, hslyr, einit, S, fsw, fsi, Ssw, Isw, q, T, qsn, Tsn, Tsf, &
                       fsens, flat, flwout, fsurf, fct, fcb, einex, niter, ierr)
      !! Flat implementation of the driver: one `do concurrent` over
      !! (category, j, i), every array explicit-shape and device-present.
      integer, intent(in) :: nx, ny, nc, nl
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: potT(nx, ny), flw(nx, ny), Qa(nx, ny)
      real(wp), intent(in) :: hilyr(nx, ny, nc), hslyr(nx, ny, nc), einit(nx, ny, nc), S(nx, ny, nc, nl)
      real(wp), intent(inout) :: fsw(nx, ny, nc), fsi(nx, ny, nc), Tsf(nx, ny, nc)
      real(wp), intent(inout) :: Ssw(nx, ny, nc, 1), qsn(nx, ny, nc, 1), Tsn(nx, ny, nc, 1)
      real(wp), intent(inout) :: Isw(nx, ny, nc, nl), q(nx, ny, nc, nl), T(nx, ny, nc, nl)
      real(wp), intent(inout) :: fsens(nx, ny, nc), flat(nx, ny, nc), flwout(nx, ny, nc), fsurf(nx, ny, nc)
      real(wp), intent(inout) :: fct(nx, ny, nc), fcb(nx, ny, nc), einex(nx, ny, nc)
      integer(int32), intent(inout) :: niter(nx, ny, nc), ierr(nx, ny, nc)
      integer :: i, j, n

      ! The loop indices reach the kernel by value (`integer, value` dummies): by reference,
      ! nvfortran 26.5 -stdpar=gpu read the wrong cells.  Both (n, j, i) and (i, j, n) are correct.
      do concurrent(n=1:nc, j=1:ny, i=1:nx)
         call ice_temperature_changes_cell(i, j, n, nx, ny, nc, nl, 1, 3600.0_wp, 1.3_wp, flw(i, j), potT(i, j), Qa(i, j), &
                                           8.0_wp, 2.0e4_wp, 0.985_wp, ICE_SFC_BULK, 0.0_wp, 0.0_wp, c, ICE_CONDUCT_BUBBLY, &
                                           .true., .true., ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, -1.8_wp, &
                                           hilyr, hslyr, S, einit, fsw, fsi, Ssw, Isw, q, T, qsn, Tsn, Tsf, &
                                           fsens, flat, flwout, fsurf, fct, fcb, einex, niter, ierr)
      end do
   end subroutine sweep_dc

   subroutine sweep_host(nx, ny, nc, nl, c, potT, flw, Qa, hilyr, hslyr, einit, S, fsw, fsi, Ssw, Isw, q, T, qsn, Tsn, Tsf, &
                         fsens, flat, flwout, fsurf, fct, fcb, einex, niter, ierr)
      !! The same sweep as plain serial loops (the reference, and the CPU
      !! timing on a GPU build).
      integer, intent(in) :: nx, ny, nc, nl
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: potT(nx, ny), flw(nx, ny), Qa(nx, ny)
      real(wp), intent(in) :: hilyr(nx, ny, nc), hslyr(nx, ny, nc), einit(nx, ny, nc), S(nx, ny, nc, nl)
      real(wp), intent(inout) :: fsw(nx, ny, nc), fsi(nx, ny, nc), Tsf(nx, ny, nc)
      real(wp), intent(inout) :: Ssw(nx, ny, nc, 1), qsn(nx, ny, nc, 1), Tsn(nx, ny, nc, 1)
      real(wp), intent(inout) :: Isw(nx, ny, nc, nl), q(nx, ny, nc, nl), T(nx, ny, nc, nl)
      real(wp), intent(inout) :: fsens(nx, ny, nc), flat(nx, ny, nc), flwout(nx, ny, nc), fsurf(nx, ny, nc)
      real(wp), intent(inout) :: fct(nx, ny, nc), fcb(nx, ny, nc), einex(nx, ny, nc)
      integer(int32), intent(inout) :: niter(nx, ny, nc), ierr(nx, ny, nc)
      integer :: i, j, n

      do n = 1, nc
         do j = 1, ny
            do i = 1, nx
               call ice_temperature_changes_cell(i, j, n, nx, ny, nc, nl, 1, 3600.0_wp, 1.3_wp, flw(i, j), potT(i, j), &
                                                 Qa(i, j), 8.0_wp, 2.0e4_wp, 0.985_wp, ICE_SFC_BULK, 0.0_wp, 0.0_wp, c, &
                                                 ICE_CONDUCT_BUBBLY, .true., .true., ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, &
                                                 -1.8_wp, hilyr, hslyr, S, einit, fsw, fsi, Ssw, Isw, q, T, qsn, Tsn, Tsf, &
                                                 fsens, flat, flwout, fsurf, fct, fcb, einex, niter, ierr)
            end do
         end do
      end do
   end subroutine sweep_host

   subroutine run_host(g, c)
      type(grid_t), intent(inout) :: g
      type(ice_constants_t), intent(in) :: c
      call sweep_host(g%nx, g%ny, g%nc, g%nl, c, g%potT, g%flw, g%Qa, g%hilyr, g%hslyr, g%einit, g%S, g%fsw, g%fsi, g%Ssw, &
                      g%Isw, g%q, g%T, g%qsn, g%Tsn, g%Tsf, g%fsens, g%flat, g%flwout, g%fsurf, g%fct, g%fcb, g%einex, &
                      g%niter, g%ierr)
   end subroutine run_host

   subroutine run_dc(g, c)
      type(grid_t), intent(inout) :: g
      type(ice_constants_t), intent(in) :: c
      call sweep_dc(g%nx, g%ny, g%nc, g%nl, c, g%potT, g%flw, g%Qa, g%hilyr, g%hslyr, g%einit, g%S, g%fsw, g%fsi, g%Ssw, &
                    g%Isw, g%q, g%T, g%qsn, g%Tsn, g%Tsf, g%fsens, g%flat, g%flwout, g%fsurf, g%fct, g%fcb, g%einex, &
                    g%niter, g%ierr)
   end subroutine run_dc

   subroutine test_device(error)
      !! The solver driven from `do concurrent (n, j, i)` over mapped
      !! `(nx, ny, ncat, nlyr)` arrays (512 x 512 x 5 categories x 7 layers =
      !! 1,310,720 columns) against the same sweep as serial host loops, on
      !! every column.  The pass count of the iteration may differ between the
      !! builds only by what fused multiply-add moves; the converged answers
      !! agree to the tolerance below (measured worst printed by the test).
      !! Reports columns per second of both sweeps.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(grid_t) :: gd, gh
      integer, parameter :: NX = 512, NY = 512, NC = 5, NL = 7
      integer(int64) :: t0, t1, t2, rate
      real(wp) :: worst, d, tdev, thost
      integer :: i, j, n, nerr, nmelt, nbad
      character(len=200) :: msg

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call alloc_grid(gd, NX, NY, NC, NL)
      call alloc_grid(gh, NX, NY, NC, NL)
      call init_grid(gd, c)
      call init_grid(gh, c)

      call system_clock(count_rate=rate)
      call system_clock(t0)
      call run_host(gh, c)
      call system_clock(t1)
      thost = real(t1 - t0, wp)/real(rate, wp)

      !$acc enter data copyin(c, gd%potT, gd%flw, gd%Qa, gd%hilyr, gd%hslyr, gd%einit, gd%S, gd%fsw, gd%fsi, gd%Ssw, gd%Isw, &
      !$acc& gd%q, gd%T, gd%qsn, gd%Tsn, gd%Tsf, gd%fsens, gd%flat, gd%flwout, gd%fsurf, gd%fct, gd%fcb, gd%einex, gd%niter, &
      !$acc& gd%ierr)
      call system_clock(t1)
      call run_dc(gd, c)
      call system_clock(t2)
      tdev = real(t2 - t1, wp)/real(rate, wp)
      !$acc update self(gd%T, gd%q, gd%Isw, gd%Tsn, gd%qsn, gd%Ssw, gd%Tsf, gd%fsw, gd%fsi, gd%fsens, gd%flat, gd%flwout, &
      !$acc& gd%fsurf, gd%fct, gd%fcb, gd%einex, gd%niter, gd%ierr)
      !$acc exit data delete(c, gd%potT, gd%flw, gd%Qa, gd%hilyr, gd%hslyr, gd%einit, gd%S, gd%fsw, gd%fsi, gd%Ssw, gd%Isw, &
      !$acc& gd%q, gd%T, gd%qsn, gd%Tsn, gd%Tsf, gd%fsens, gd%flat, gd%flwout, gd%fsurf, gd%fct, gd%fcb, gd%einex, gd%niter, &
      !$acc& gd%ierr)

write (output_unit, '(a,i0,a,f8.3,a,f10.0,a,f8.3,a,f10.0,a)') "  bl99 solver sweep, ", NX*NY*NC, " columns: do concurrent sweep ", &
         tdev, " s (", real(NX*NY*NC, wp)/tdev, " columns/s); serial host loop ", thost, " s (", real(NX*NY*NC, wp)/thost, &
         " columns/s)"

      nerr = count(gh%ierr /= 0_int32)
      nmelt = count(gh%Tsf == 0.0_wp)
      call check(error, nerr == 0, "device grid: no error bits on the host reference")
      if (allocated(error)) return
      call check(error, all(gd%ierr == gh%ierr), "device grid: error bits equal the host's")
      if (allocated(error)) return
      worst = 0.0_wp
      nbad = 0
      call check(error, all(gd%niter == gh%niter), "device grid: pass counts equal the host's, every column")
      if (allocated(error)) return
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               d = max(cell_diff4(gd%T, gh%T, i, j, n, NL, 1.0_wp), cell_diff4(gd%q, gh%q, i, j, n, NL, 1.0e8_wp), &
                       cell_diff4(gd%Isw, gh%Isw, i, j, n, NL, 1.0_wp), cell_diff4(gd%Tsn, gh%Tsn, i, j, n, 1, 1.0_wp), &
                       cell_diff4(gd%qsn, gh%qsn, i, j, n, 1, 1.0e7_wp), cell_diff4(gd%Ssw, gh%Ssw, i, j, n, 1, 1.0_wp), &
                       cell_diff3(gd%Tsf, gh%Tsf, i, j, n), cell_diff3(gd%fsw, gh%fsw, i, j, n), &
                       cell_diff3(gd%fsi, gh%fsi, i, j, n), cell_diff3(gd%fsens, gh%fsens, i, j, n), &
                       cell_diff3(gd%flat, gh%flat, i, j, n), cell_diff3(gd%flwout, gh%flwout, i, j, n), &
                       cell_diff3(gd%fsurf, gh%fsurf, i, j, n), cell_diff3(gd%fct, gh%fct, i, j, n), &
                       cell_diff3(gd%fcb, gh%fcb, i, j, n), cell_diff3(gd%einex, gh%einex, i, j, n))
               if (.not. ieee_is_finite(d)) then
                  call test_failed(error, "device grid: non-finite difference")
                  return
               end if
               if (d > worst) worst = d
               if (d > TOL_DEVICE) nbad = nbad + 1
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,i0,a,i0)') "  bl99 solver grid: worst difference ", worst, ", columns above tolerance ", nbad, &
         ", melting columns ", nmelt
      write (output_unit, '(a)') trim(msg)
      call check(error, nmelt > 0 .and. nmelt < NX*NY*NC, "device grid: both cold and melting columns")
      if (allocated(error)) return
      call check(error, worst <= TOL_DEVICE, "device sweep == host sweep")
   end subroutine test_device

   pure function cell_diff3(a, b, i, j, n) result(d)
      !! Difference of one cell of two 3-D fields, normalised by max(|b|, 1).
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      integer, intent(in) :: i, j, n
      real(wp) :: d
      d = abs(a(i, j, n) - b(i, j, n))/max(abs(b(i, j, n)), 1.0_wp)
   end function cell_diff3

   pure function cell_diff4(a, b, i, j, n, nk, scale) result(d)
      !! Worst difference over the layers of one cell of two 4-D fields, normalised by max(|b|, scale).
      real(wp), intent(in) :: a(:, :, :, :), b(:, :, :, :), scale
      integer, intent(in) :: i, j, n, nk
      real(wp) :: d
      integer :: k
      d = 0.0_wp
      do k = 1, nk
         d = max(d, abs(a(i, j, n, k) - b(i, j, n, k))/max(abs(b(i, j, n, k)), scale))
      end do
   end function cell_diff4

   subroutine alloc_grid(g, nx, ny, nc, nl)
      type(grid_t), intent(out) :: g
      integer, intent(in) :: nx, ny, nc, nl
      g%nx = nx; g%ny = ny; g%nc = nc; g%nl = nl
      allocate (g%potT(nx, ny), g%flw(nx, ny), g%Qa(nx, ny))
      allocate (g%hilyr(nx, ny, nc), g%hslyr(nx, ny, nc), g%einit(nx, ny, nc), g%fsw(nx, ny, nc), g%fsi(nx, ny, nc), &
                g%Tsf(nx, ny, nc), g%fsens(nx, ny, nc), g%flat(nx, ny, nc), g%flwout(nx, ny, nc), g%fsurf(nx, ny, nc), &
                g%fct(nx, ny, nc), g%fcb(nx, ny, nc), g%einex(nx, ny, nc), g%niter(nx, ny, nc), g%ierr(nx, ny, nc))
      allocate (g%S(nx, ny, nc, nl), g%Isw(nx, ny, nc, nl), g%q(nx, ny, nc, nl), g%T(nx, ny, nc, nl))
      allocate (g%Ssw(nx, ny, nc, 1), g%qsn(nx, ny, nc, 1), g%Tsn(nx, ny, nc, 1))
   end subroutine alloc_grid

end module test_cice_bl99_solver
