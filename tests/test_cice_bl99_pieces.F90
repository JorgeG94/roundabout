module test_cice_bl99_pieces
   !! Tests of the BL99 building blocks (`rdb_ice_bl99`, `rdb_ice_vprofile`):
   !! conductivity closed forms for both options (fresh ice, the melting
   !! point, the bubbly minimum), the Thomas solver against an independent
   !! dense elimination on random diagonally dominant systems and hand cases,
   !! the matrix-row properties the scheme guarantees (row sums, scaled
   !! symmetry, diagonal dominance), A3 steady conduction (the rows marched
   !! to steady state reproduce the analytic profile and the uniform flux for
   !! constant and layer-varying conductivity, melting and cold surface, with
   !! and without snow), the vertical-profile closed forms and repairs,
   !! NaN-proof error flagging with the columns left untouched, and a
   !! `do concurrent` device test.  Expected values are closed forms written
   !! here, never the port's own output.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_NEG_THICK, ICE_ERR_TSN_RANGE, ICE_ERR_SIN_LOW, ICE_ERR_TIN_RANGE
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71, ice_q_from_T
   use rdb_ice_bl99, only: ICE_NMAT_MAX, ICE_KIMIN, ICE_BETAK, ice_layer_conductivity, ice_conductivity, &
                           ice_matrix_calc_Tsfc, ice_tridiag_solver, ice_conductivity_core, ice_matrix_calc_Tsfc_core, &
                           ice_tridiag_solver_core
   use rdb_ice_vprofile, only: ICE_TMIN, ice_init_vertical_profile
   use lcg_deterministic, only: lcg_next, lcg_unit
   use test_cice_util, only: ice_random_columns
   implicit none
   private

   public :: collect_cice_bl99_pieces_tests

contains

   subroutine collect_cice_bl99_pieces_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("conductivity_closed_forms", test_cond_closed), &
                  new_unittest("conductivity_limits", test_cond_limits), &
                  new_unittest("interface_conductances", test_kh_closed), &
                  new_unittest("tridiag_hand_cases", test_tri_hand), &
                  new_unittest("tridiag_vs_dense_random", test_tri_dense), &
                  new_unittest("matrix_row_properties", test_matrix_props), &
                  new_unittest("steady_conduction_constant_k", test_a3_const), &
                  new_unittest("steady_conduction_snow_ice", test_a3_snow), &
                  new_unittest("steady_conduction_layered_k", test_a3_layered), &
                  new_unittest("steady_conduction_cold_surface", test_a3_cold), &
                  new_unittest("vprofile_closed_forms", test_vp_closed), &
                  new_unittest("vprofile_repairs", test_vp_repairs), &
                  new_unittest("vprofile_range_errors", test_vp_range), &
                  new_unittest("nan_in_flags_and_leaves_columns", test_nan), &
                  new_unittest("core_tridiag_nonfinite_flags", test_core_tridiag_nonfinite), &
                  new_unittest("core_chain_nonfinite_flags", test_core_chain_nonfinite), &
                  new_unittest("wrappers_equal_cores_bitwise", test_wrappers_equal_cores), &
                  new_unittest("device_chain_matches_host", test_device) &
                  ]
   end subroutine collect_cice_bl99_pieces_tests

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   pure function q_fwd(T, S, c) result(q)
      !! BL99 ice enthalpy (brine form), written out independently here.
      real(wp), intent(in) :: T, S
      type(ice_constants_t), intent(in) :: c
      real(wp) :: q, Tm
      Tm = -c%depressT*S
      q = -c%rhoi*(c%cp_ice*(Tm - T) + c%Lfresh*(1.0_wp - Tm/T) - c%cp_ocn*Tm)
   end function q_fwd

   pure function qs_fwd(T, c) result(q)
      real(wp), intent(in) :: T
      type(ice_constants_t), intent(in) :: c
      real(wp) :: q
      q = -c%rhos*(c%Lfresh - c%cp_ice*T)
   end function qs_fwd

   function rnd(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function rnd

   subroutine dense_solve(n, a_in, b, x)
      !! Gaussian elimination with partial pivoting (independent reference).
      integer, intent(in) :: n
      real(wp), intent(in) :: a_in(n, n), b(n)
      real(wp), intent(out) :: x(n)
      real(wp) :: a(n, n), r(n), f, tmp(n), tr
      integer :: i, j, p
      a = a_in
      r = b
      do j = 1, n
         p = j
         do i = j + 1, n
            if (abs(a(i, j)) > abs(a(p, j))) p = i
         end do
         if (p /= j) then
            tmp = a(j, :); a(j, :) = a(p, :); a(p, :) = tmp
            tr = r(j); r(j) = r(p); r(p) = tr
         end if
         do i = j + 1, n
            f = a(i, j)/a(j, j)
            a(i, :) = a(i, :) - f*a(j, :)
            r(i) = r(i) - f*r(j)
         end do
      end do
      do i = n, 1, -1
         x(i) = r(i)
         do j = i + 1, n
            x(i) = x(i) - a(i, j)*x(j)
         end do
         x(i) = x(i)/a(i, i)
      end do
   end subroutine dense_solve

   ! ------------------------------------------------------------------
   ! conductivity
   ! ------------------------------------------------------------------

   subroutine test_cond_closed(error)
      !! Hand evaluations of both forms at (T, S) points.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: k

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ! bubbly: (2.11 - 0.011 T + 0.09 S/T)*rhoi/917; rhoi = 917 so the factor is 1
      k = ice_layer_conductivity(-10.0_wp, 0.0_wp, c, ICE_CONDUCT_BUBBLY)
      call check(error, abs(k - 2.22_wp) <= 1.0e-14_wp, "bubbly fresh -10C = 2.11+0.11")
      if (allocated(error)) return
      k = ice_layer_conductivity(-2.0_wp, 3.2_wp, c, ICE_CONDUCT_BUBBLY)
      call check(error, abs(k - (2.11_wp + 0.022_wp - 0.144_wp)) <= 1.0e-14_wp, "bubbly saline -2C")
      if (allocated(error)) return
      c%rhoi = 900.0_wp
      k = ice_layer_conductivity(-10.0_wp, 0.0_wp, c, ICE_CONDUCT_BUBBLY)
      call check(error, abs(k - 2.22_wp*900.0_wp/917.0_wp) <= 1.0e-14_wp, "bubbly scales with rhoi/917")
      if (allocated(error)) return
      ! MU71: kice + 0.13 S/T
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      k = ice_layer_conductivity(-10.0_wp, 3.2_wp, c, ICE_CONDUCT_MU71)
      call check(error, abs(k - (2.03_wp - 0.0416_wp)) <= 1.0e-14_wp, "MU71 saline -10C")
      if (allocated(error)) return
      k = ice_layer_conductivity(-0.5_wp, 1.0_wp, c, ICE_CONDUCT_MU71)
      call check(error, abs(k - (2.03_wp - 0.26_wp)) <= 1.0e-14_wp, "MU71 -0.5C S=1")
      if (allocated(error)) return
      call check(error, ICE_KIMIN == 0.10_wp .and. ICE_BETAK == 0.13_wp, "constants of record")
   end subroutine test_cond_closed

   subroutine test_cond_limits(error)
      !! Fresh ice, the melting point and the bubbly minimum.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: T, k
      integer :: i

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ! fresh ice: MU71 is exactly kice at any temperature
      do i = 1, 6
         T = -real(i*i, wp)
         k = ice_layer_conductivity(T, 0.0_wp, c, ICE_CONDUCT_MU71)
         call check(error, k == c%kice, "MU71 fresh ice is kice")
         if (allocated(error)) return
      end do
      ! fresh bubbly is linear in T: 2.11 - 0.011 T
      k = ice_layer_conductivity(-30.0_wp, 0.0_wp, c, ICE_CONDUCT_BUBBLY)
      call check(error, abs(k - 2.44_wp) <= 1.0e-14_wp, "bubbly fresh -30C")
      if (allocated(error)) return
      ! at (and above) the melting point the 1/T term is guarded by min(-puny, T)
      ! and the conductivity falls to the minimum for saline ice
      do i = 1, 3
         select case (i)
         case (1); T = -1.0e-12_wp
         case (2); T = 0.0_wp
         case default; T = 2.0_wp
         end select
         k = ice_layer_conductivity(T, 3.2_wp, c, ICE_CONDUCT_BUBBLY)
         call check(error, k == ICE_KIMIN, "bubbly saline at the melting point -> kimin")
         if (allocated(error)) return
         k = ice_layer_conductivity(T, 3.2_wp, c, ICE_CONDUCT_MU71)
         call check(error, k == ICE_KIMIN, "MU71 saline at the melting point -> kimin")
         if (allocated(error)) return
      end do
      ! S=3.2 at -0.06 degC: 2.11 + 0.00066 - 0.09*3.2/0.06 is negative, so the floor is active
      k = ice_layer_conductivity(-0.06_wp, 3.2_wp, c, ICE_CONDUCT_BUBBLY)
      call check(error, k == ICE_KIMIN, "S=3.2 at -0.06C: the formula is negative, floor active")
   end subroutine test_cond_limits

   subroutine test_kh_closed(error)
      !! Interface conductances: harmonic means, ends, snow and no snow.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      real(wp) :: kh(ICE_NMAT_MAX), T(4), S(4), hi, hs, k(4), exp1, exp2
      integer :: j

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      T = [-20.0_wp, -12.0_wp, -6.0_wp, -3.0_wp]
      S = [0.1_wp, 1.0_wp, 2.0_wp, 3.2_wp]
      hi = 0.8_wp
      hs = 0.15_wp
      do j = 1, 4
         k(j) = max(2.11_wp - 0.011_wp*T(j) + 0.09_wp*S(j)/T(j), ICE_KIMIN)
      end do
      ierr = 0_int32
      kh = -1.0_wp
      call ice_conductivity(.true., 4, 1, hi, hs, T, S, c, ICE_CONDUCT_BUBBLY, kh, ierr)
      call check(error, ierr == 0_int32, "no error")
      if (allocated(error)) return
      call check(error, abs(kh(1) - 2.0_wp*c%ksno/hs) <= 1.0e-13_wp, "snow top 2 ksno/hs")
      if (allocated(error)) return
      exp1 = 2.0_wp*c%ksno*k(1)/(c%ksno*hi + k(1)*hs)
      call check(error, abs(kh(2) - exp1) <= 1.0e-13_wp, "snow/ice interface")
      if (allocated(error)) return
      do j = 2, 4
         exp2 = 2.0_wp*k(j - 1)*k(j)/((k(j - 1) + k(j))*hi)
         call check(error, abs(kh(j + 1) - exp2) <= 1.0e-13_wp, "interior ice interface")
         if (allocated(error)) return
      end do
      call check(error, abs(kh(6) - 2.0_wp*k(4)/hi) <= 1.0e-13_wp, "ice base 2 k/h")
      if (allocated(error)) return
      ! no snow
      kh = -1.0_wp
      call ice_conductivity(.false., 4, 1, hi, 0.0_wp, T, S, c, ICE_CONDUCT_BUBBLY, kh, ierr)
      call check(error, kh(1) == 0.0_wp .and. abs(kh(2) - 2.0_wp*k(1)/hi) <= 1.0e-13_wp, "no snow: top of ice 2 k/h")
      if (allocated(error)) return
      ! uniform column: interior kh = k/h
      T = -10.0_wp
      S = 0.0_wp
      call ice_conductivity(.false., 4, 1, hi, 0.0_wp, T, S, c, ICE_CONDUCT_MU71, kh, ierr)
      call check(error, abs(kh(3) - c%kice/hi) <= 1.0e-14_wp .and. abs(kh(2) - 2.0_wp*c%kice/hi) <= 1.0e-14_wp, &
                 "uniform column harmonic mean = k/h")
   end subroutine test_kh_closed

   ! ------------------------------------------------------------------
   ! tridiagonal solver
   ! ------------------------------------------------------------------

   subroutine test_tri_hand(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x(ICE_NMAT_MAX)
      integer(int32) :: ierr

      ierr = 0_int32
      ! 3x3: [2 -1 0; -1 2 -1; 0 -1 2] x = [1 0 1] -> x = [1 1 1]
      call ice_tridiag_solver(3, [0.0_wp, -1.0_wp, -1.0_wp], [2.0_wp, 2.0_wp, 2.0_wp], &
                              [-1.0_wp, -1.0_wp, 0.0_wp], [1.0_wp, 0.0_wp, 1.0_wp], x, ierr)
      call check(error, ierr == 0_int32 .and. all(abs(x(1:3) - 1.0_wp) <= 1.0e-15_wp), "Laplacian 3x3")
      if (allocated(error)) return
      ! 2x2: [4 1; 2 3] x = [6 8] -> x = [1 2]
      call ice_tridiag_solver(2, [0.0_wp, 2.0_wp], [4.0_wp, 3.0_wp], [1.0_wp, 0.0_wp], [6.0_wp, 8.0_wp], x, ierr)
      call check(error, abs(x(1) - 1.0_wp) <= 1.0e-15_wp .and. abs(x(2) - 2.0_wp) <= 1.0e-15_wp, "2x2")
      if (allocated(error)) return
      ! 1x1
      call ice_tridiag_solver(1, [0.0_wp], [4.0_wp], [0.0_wp], [10.0_wp], x, ierr)
      call check(error, x(1) == 2.5_wp, "1x1")
      if (allocated(error)) return
      ! diagonal system (the dummy rows of the matrix builder): x = rhs/diag
      call ice_tridiag_solver(3, [0.0_wp, 0.0_wp, 0.0_wp], [1.0_wp, 2.0_wp, 4.0_wp], [0.0_wp, 0.0_wp, 0.0_wp], &
                              [3.0_wp, 3.0_wp, 3.0_wp], x, ierr)
      call check(error, x(1) == 3.0_wp .and. x(2) == 1.5_wp .and. x(3) == 0.75_wp .and. ierr == 0_int32, "diagonal system")
   end subroutine test_tri_hand

   subroutine test_tri_dense(error)
      !! Random diagonally dominant systems, nmat = 1..ICE_NMAT_MAX, 100 each.
      !! Measured worst relative difference vs dense elimination is stated in
      !! the tolerance (measured worst 6.4e-16 relative to max|x|, tolerance
      !! 1e-14; the row residual measured 2.1e-14, tolerance 3e-13).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TOL = 1.0e-14_wp
      integer :: n, rep, k
      integer(int64) :: st
      real(wp) :: sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX), x(ICE_NMAT_MAX)
      real(wp) :: a(ICE_NMAT_MAX, ICE_NMAT_MAX), xd(ICE_NMAT_MAX), worst, d, res, rmax
      integer(int32) :: ierr

      st = 424242_int64
      worst = 0.0_wp
      rmax = 0.0_wp
      ierr = 0_int32
      do n = 1, ICE_NMAT_MAX
         do rep = 1, 100
            a = 0.0_wp
            do k = 1, n
               sb(k) = -3.0_wp*rnd(st)
               sp(k) = -3.0_wp*rnd(st)
               rh(k) = 60.0_wp*rnd(st) - 30.0_wp
            end do
            sb(1) = 0.0_wp
            sp(n) = 0.0_wp
            do k = 1, n
               dg(k) = 1.0_wp - sb(k) - sp(k) + 2.0_wp*rnd(st)
               a(k, k) = dg(k)
               if (k > 1) a(k, k - 1) = sb(k)
               if (k < n) a(k, k + 1) = sp(k)
            end do
            x = 0.0_wp
            call ice_tridiag_solver(n, sb, dg, sp, rh, x, ierr)
            call dense_solve(n, a(1:n, 1:n), rh(1:n), xd(1:n))
            d = maxval(abs(x(1:n) - xd(1:n)))/max(maxval(abs(xd(1:n))), 1.0_wp)
            if (.not. (d <= worst)) worst = d
            ! residual of the returned solution against the rows
            res = 0.0_wp
            do k = 1, n
               d = dg(k)*x(k) - rh(k)
               if (k > 1) d = d + sb(k)*x(k - 1)
               if (k < n) d = d + sp(k)*x(k + 1)
               if (.not. (abs(d) <= res)) res = abs(d)
            end do
            if (.not. (res <= rmax)) rmax = res
         end do
      end do
      call check(error, ierr == 0_int32, "no error bits")
      if (allocated(error)) return
      call check(error, worst <= TOL, "Thomas == dense elimination")
      if (allocated(error)) return
      call check(error, rmax <= 3.0e-13_wp, "row residual of the Thomas solution")
   end subroutine test_tri_dense

   ! ------------------------------------------------------------------
   ! matrix rows
   ! ------------------------------------------------------------------

   subroutine test_matrix_props(error)
      !! For every layer count 1..8, snow/no snow and cold/melting surface:
      !! row sums, scaled symmetry of the layer couplings (the scheme is the
      !! symmetric conduction operator divided by the layer capacities) and
      !! strict diagonal dominance of the layer rows.
      type(error_type), allocatable, intent(out) :: error
      integer :: nil, is, ic, k, nmat, rep, r
      integer(int64) :: st
      logical :: l_snow, l_cold
      real(wp) :: Tin(8), Tsn(1), kh(ICE_NMAT_MAX), Ssw(1), Isw(8), etai(8), etas(1)
      real(wp) :: sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX)
      real(wp) :: Tsf, fs, dfs, e1, e2, wsum, wsym, wdom
      integer(int32) :: ierr

      st = 99_int64
      wsum = 0.0_wp
      wsym = 0.0_wp
      wdom = huge(1.0_wp)
      ierr = 0_int32
      do nil = 1, 8
         nmat = nil + 2
         do is = 0, 1
            do ic = 0, 1
               l_snow = is == 1
               l_cold = ic == 1
               do rep = 1, 5
                  do k = 1, nil
                     Tin(k) = -30.0_wp*rnd(st)
                     Isw(k) = 5.0_wp*rnd(st)
                     etai(k) = 0.001_wp + 0.05_wp*rnd(st)
                  end do
                  Tsn(1) = -20.0_wp*rnd(st)
                  Ssw(1) = 5.0_wp*rnd(st)
                  etas(1) = 0.005_wp + 0.3_wp*rnd(st)
                  do k = 1, nmat
                     kh(k) = 0.5_wp + 40.0_wp*rnd(st)
                  end do
                  if (.not. l_snow) kh(1) = 0.0_wp
                  Tsf = -10.0_wp*rnd(st)
                  fs = 100.0_wp*rnd(st) - 50.0_wp
                  dfs = -(3.0_wp + 12.0_wp*rnd(st))
                  call ice_matrix_calc_Tsfc(l_snow, l_cold, nil, 1, Tsf, -1.8_wp, fs, dfs, Tin, Tsn, kh, Ssw, Isw, &
                                            etai, etas, sb, dg, sp, rh, ierr)
                  ! row sums: rows are 1 + eta*(kh(k)+kh(k+1)) - eta*kh(k) - eta*kh(k+1) = 1
                  ! except where a neighbour is the pinned Tsf (top) or Tbot (bottom)
                  do r = 3, nmat
                     k = r - 2
                     e1 = dg(r) + sb(r) + sp(r) - 1.0_wp
                     if (r == 3 .and. .not. (l_snow .or. l_cold)) then
                        e2 = etai(1)*kh(2)         ! pinned top neighbour (melting, no snow)
                     else
                        e2 = 0.0_wp
                     end if
                     if (r == nmat) e2 = e2 + etai(nil)*kh(nmat)       ! Tbot neighbour
                     if (.not. (abs(e1 - e2) <= wsum)) wsum = abs(e1 - e2)
                     ! strict dominance of the layer rows: diag > |sb| + |sp|
                     if (.not. (dg(r) - abs(sb(r)) - abs(sp(r)) >= 0.0_wp)) wdom = -1.0_wp
                     wdom = min(wdom, dg(r) - abs(sb(r)) - abs(sp(r)))
                  end do
                  if (l_snow) then
                     e1 = dg(2) + sb(2) + sp(2) - 1.0_wp
                     e2 = merge(0.0_wp, etas(1)*kh(1), l_cold)
                     if (.not. (abs(e1 - e2) <= wsum)) wsum = abs(e1 - e2)
                  end if
                  if (l_cold) then
                     r = merge(1, 2, l_snow)
                     e1 = dg(r) + sb(r) + sp(r) - dfs
                     if (.not. (abs(e1) <= wsum)) wsum = abs(e1)
                  end if
                  ! scaled symmetry: the coupling kh between rows r and r+1 appears as
                  ! sp(r)/s_r and sb(r+1)/s_(r+1) with s = eta for a layer row and s = -1 for
                  ! the cold surface row, which must agree
                  if (l_snow .and. l_cold) then                  ! surface -> snow
                     e1 = sp(1)
                     e2 = -sb(2)/etas(1)
                     if (.not. (abs(e1 - e2) <= wsym)) wsym = abs(e1 - e2)
                  end if
                  if (l_snow) then                               ! snow -> top ice
                     e1 = sp(2)/etas(1)
                     e2 = sb(3)/etai(1)
                     if (.not. (abs(e1 - e2) <= wsym)) wsym = abs(e1 - e2)
                  else if (l_cold) then                          ! surface (row 2) -> top ice
                     e1 = sp(2)
                     e2 = -sb(3)/etai(1)
                     if (.not. (abs(e1 - e2) <= wsym)) wsym = abs(e1 - e2)
                  end if
                  do r = 3, nmat - 1                             ! ice -> ice
                     e1 = sp(r)/etai(r - 2)
                     e2 = sb(r + 1)/etai(r - 1)
                     if (.not. (abs(e1 - e2) <= wsym)) wsym = abs(e1 - e2)
                  end do
                  ! the symmetric part itself is the conductance: sp(r)/eta = -kh(r+1)
                  do r = 3, nmat - 1
                     e1 = sp(r)/etai(r - 2) + kh(r)
                     if (.not. (abs(e1) <= wsym)) wsym = abs(e1)
                  end do
                  if (.not. (wdom >= 0.0_wp)) wdom = -1.0_wp
               end do
            end do
         end do
      end do
      call check(error, ierr == 0_int32, "no error bits")
      if (allocated(error)) return
      call check(error, wsum <= 2.0e-14_wp, "row sums")
      if (allocated(error)) return
      call check(error, wsym <= 1.0e-13_wp, "scaled symmetry of the couplings")
      if (allocated(error)) return
      call check(error, wdom >= 1.0_wp - 1.0e-14_wp, "layer rows dominant: diag - |sb| - |sp| >= 1")
   end subroutine test_matrix_props

   ! ------------------------------------------------------------------
   ! A3: steady conduction
   ! ------------------------------------------------------------------

   subroutine march_to_steady(nil, l_snow, l_cold, Tsf, Tbot, lam, Ta, kh, etai0, etas0, Tin, Tsn, nstep, ok, ierr)
      !! March the backward-Euler rows (no absorbed shortwave, linear
      !! surface `fsurf = lam*(Ta - Tsf)` for a cold surface, pinned `Tsf`
      !! otherwise) until no temperature changes by more than 1e-14 K.
      integer, intent(in) :: nil
      logical, intent(in) :: l_snow, l_cold
      real(wp), intent(inout) :: Tsf
      real(wp), intent(in) :: Tbot, lam, Ta, kh(nil + 2), etai0, etas0
      real(wp), intent(inout) :: Tin(nil), Tsn
      integer, intent(out) :: nstep
      logical, intent(out) :: ok
      integer(int32), intent(inout) :: ierr
      real(wp) :: sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX), x(ICE_NMAT_MAX)
      real(wp) :: etai(nil), Isw(nil), Ssw(1), Tsn1(1), Tnew(nil), dmax, d
      integer :: nmat

      nmat = nil + 2
      etai = etai0
      Isw = 0.0_wp
      Ssw = 0.0_wp
      ok = .false.
      do nstep = 1, 200000
         Tsn1(1) = Tsn
         call ice_matrix_calc_Tsfc(l_snow, l_cold, nil, 1, Tsf, Tbot, lam*(Ta - Tsf), -lam, Tin, Tsn1, kh, Ssw, Isw, &
                                   etai, [etas0], sb, dg, sp, rh, ierr)
         call ice_tridiag_solver(nmat, sb, dg, sp, rh, x, ierr)
         Tnew = x(3:nmat)
         dmax = maxval(abs(Tnew - Tin))
         if (l_snow) then
            d = abs(x(2) - Tsn)
            if (.not. (d <= dmax)) dmax = d
            Tsn = x(2)
         end if
         if (l_cold) then
            d = abs(x(merge(1, 2, l_snow)) - Tsf)
            if (.not. (d <= dmax)) dmax = d
            Tsf = x(merge(1, 2, l_snow))
         end if
         Tin = Tnew
         if (ierr /= 0_int32) return
         if (dmax <= 1.0e-14_wp) then
            ok = .true.
            return
         end if
      end do
   end subroutine march_to_steady

   subroutine test_a3_const(error)
      !! Fixed top and bottom temperatures, constant conductivity (fresh ice,
      !! MU71): the steady state is the straight line, the flux is uniform
      !! and equal to k*DeltaT/H.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer :: nil, j, nstep
      real(wp) :: kh(ICE_NMAT_MAX), Tz(8), S(8), Tin(8), Tsn, Tsf, Tbot, hi, hl, F, worst, d, ferr
      integer(int32) :: ierr
      logical :: ok

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      Tsf = -20.0_wp
      Tbot = -1.8_wp
      hi = 1.6_wp
      worst = 0.0_wp
      ferr = 0.0_wp
      ierr = 0_int32
      do nil = 1, 8
         hl = hi/real(nil, wp)
         Tz(1:nil) = -5.0_wp
         S(1:nil) = 0.0_wp
         call ice_conductivity(.false., nil, 1, hl, 0.0_wp, Tz(1:nil), S(1:nil), c, ICE_CONDUCT_MU71, kh, ierr)
         Tin(1:nil) = -8.0_wp       ! start away from the answer
         Tsn = 0.0_wp
         call march_to_steady(nil, .false., .false., Tsf, Tbot, 0.0_wp, 0.0_wp, kh, 8.0_wp, 1.0_wp, Tin(1:nil), Tsn, &
                              nstep, ok, ierr)
         call check(error, ok .and. ierr == 0_int32, "marching converged")
         if (allocated(error)) return
         F = c%kice*(Tsf - Tbot)/hi
         do j = 1, nil
            d = abs(Tin(j) - (Tsf + (Tbot - Tsf)*(real(j, wp) - 0.5_wp)/real(nil, wp)))
            if (.not. (d <= worst)) worst = d
         end do
         ! flux through every interface: top 2k/h half cell, interior k/h, base 2k/h
         d = abs(kh(2)*(Tsf - Tin(1)) - F)/abs(F)
         if (.not. (d <= ferr)) ferr = d
         do j = 2, nil
            d = abs(kh(j + 1)*(Tin(j - 1) - Tin(j)) - F)/abs(F)
            if (.not. (d <= ferr)) ferr = d
         end do
         d = abs(kh(nil + 2)*(Tin(nil) - Tbot) - F)/abs(F)
         if (.not. (d <= ferr)) ferr = d
      end do
      call check(error, ferr <= 1.0e-13_wp, "uniform flux = k dT/H through every interface")
      if (allocated(error)) return
      call check(error, worst <= 1.0e-13_wp, "steady profile linear and flux uniform = k dT/H")
   end subroutine test_a3_const

   subroutine test_a3_snow(error)
      !! Piecewise conductivity: snow (ksno) over fresh ice (kice), pinned
      !! top and bottom.  F = DeltaT/(hs/ksno + hi/kice); layer centres follow
      !! from the cumulative thermal resistance.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: kh(ICE_NMAT_MAX), Tz(NIL), S(NIL), Tin(NIL), Tsn, Tsf, Tbot, hi, hs, hl, F, R, worst, d, Fi
      integer(int32) :: ierr
      integer :: j, nstep
      logical :: ok

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      Tsf = -25.0_wp
      Tbot = -1.8_wp
      hi = 1.2_wp
      hs = 0.3_wp
      hl = hi/real(NIL, wp)
      Tz = -5.0_wp
      S = 0.0_wp
      ierr = 0_int32
      call ice_conductivity(.true., NIL, 1, hl, hs, Tz, S, c, ICE_CONDUCT_MU71, kh, ierr)
      Tin = -10.0_wp
      Tsn = -15.0_wp
      call march_to_steady(NIL, .true., .false., Tsf, Tbot, 0.0_wp, 0.0_wp, kh, 8.0_wp, 5.0_wp, Tin, Tsn, nstep, ok, ierr)
      call check(error, ok .and. ierr == 0_int32, "marching converged")
      if (allocated(error)) return
      R = hs/c%ksno + hi/c%kice
      F = (Tsf - Tbot)/R
      worst = abs(Tsn - (Tsf - F*hs/(2.0_wp*c%ksno)))
      do j = 1, NIL
         d = abs(Tin(j) - (Tsf - F*(hs/c%ksno + (real(j, wp) - 0.5_wp)*hl/c%kice)))
         if (.not. (d <= worst)) worst = d
      end do
      call check(error, worst <= 1.0e-13_wp, "piecewise-linear steady profile")
      if (allocated(error)) return
      Fi = kh(1)*(Tsf - Tsn)
      call check(error, abs(Fi - F) <= 1.0e-13_wp*abs(F), "flux at snow top")
      if (allocated(error)) return
      Fi = kh(2)*(Tsn - Tin(1))
      call check(error, abs(Fi - F) <= 1.0e-13_wp*abs(F), "flux at the snow/ice interface")
      if (allocated(error)) return
      Fi = kh(NIL + 2)*(Tin(NIL) - Tbot)
      call check(error, abs(Fi - F) <= 1.0e-13_wp*abs(F), "flux at the ice base")
   end subroutine test_a3_snow

   subroutine test_a3_layered(error)
      !! Layer-varying (bubbly, T and S dependent) conductivity frozen at a
      !! prescribed profile, as the solver holds it over an iteration: the
      !! steady flux is DeltaT/sum(h/k_j) and the centres follow the
      !! cumulative resistance.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: kh(ICE_NMAT_MAX), Tz(NIL), S(NIL), Tin(NIL), Tsn, Tsf, Tbot, hi, hl, F, kj(NIL), cum, worst, d, Fi
      integer(int32) :: ierr
      integer :: j, nstep
      logical :: ok

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      Tsf = -18.0_wp
      Tbot = -1.8_wp
      hi = 2.0_wp
      hl = hi/real(NIL, wp)
      Tz = [-15.0_wp, -9.0_wp, -4.0_wp, -2.0_wp]
      S = [0.1_wp, 1.5_wp, 2.5_wp, 3.2_wp]
      ierr = 0_int32
      call ice_conductivity(.false., NIL, 1, hl, 0.0_wp, Tz, S, c, ICE_CONDUCT_BUBBLY, kh, ierr)
      do j = 1, NIL
         kj(j) = max(2.11_wp - 0.011_wp*Tz(j) + 0.09_wp*S(j)/Tz(j), ICE_KIMIN)
      end do
      call check(error, kj(1) > kj(4), "conductivity actually varies across the column")
      if (allocated(error)) return
      Tin = -5.0_wp
      Tsn = 0.0_wp
      call march_to_steady(NIL, .false., .false., Tsf, Tbot, 0.0_wp, 0.0_wp, kh, 6.0_wp, 1.0_wp, Tin, Tsn, nstep, ok, ierr)
      call check(error, ok .and. ierr == 0_int32, "marching converged")
      if (allocated(error)) return
      F = (Tsf - Tbot)/sum(hl/kj)
      worst = 0.0_wp
      cum = 0.0_wp
      do j = 1, NIL
         d = abs(Tin(j) - (Tsf - F*(cum + 0.5_wp*hl/kj(j))))
         if (.not. (d <= worst)) worst = d
         cum = cum + hl/kj(j)
      end do
      call check(error, worst <= 1.0e-13_wp, "centres follow the cumulative resistance")
      if (allocated(error)) return
      do j = 2, NIL
         Fi = kh(j + 1)*(Tin(j - 1) - Tin(j))
         call check(error, abs(Fi - F) <= 1.0e-13_wp*abs(F), "uniform flux through an interior interface")
         if (allocated(error)) return
      end do
      Fi = kh(2)*(Tsf - Tin(1))
      call check(error, abs(Fi - F) <= 1.0e-13_wp*abs(F), "flux at the top")
   end subroutine test_a3_layered

   subroutine test_a3_cold(error)
      !! Cold-surface row with a linear surface balance `lam*(Ta - Tsf)`:
      !! the Newton step is exact, so the steady state is
      !! F = (Ta - Tbot)/(1/lam + R), Tsf = Ta - F/lam; with and without snow.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: kh(ICE_NMAT_MAX), Tz(NIL), S(NIL), Tin(NIL), Tsn, Tsf, Tbot, hi, hs, hl, F, R, worst, Ta, lam
      integer(int32) :: ierr
      integer :: j, nstep, ic
      logical :: ok, l_snow

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      Ta = -22.0_wp
      lam = 12.0_wp
      Tbot = -1.8_wp
      hi = 1.0_wp
      hs = 0.2_wp
      hl = hi/real(NIL, wp)
      Tz = -5.0_wp
      S = 0.0_wp
      ierr = 0_int32
      do ic = 0, 1
         l_snow = ic == 1
         call ice_conductivity(l_snow, NIL, 1, hl, merge(hs, 0.0_wp, l_snow), Tz, S, c, ICE_CONDUCT_MU71, kh, ierr)
         Tin = -9.0_wp
         Tsn = -9.0_wp
         Tsf = -9.0_wp
         call march_to_steady(NIL, l_snow, .true., Tsf, Tbot, lam, Ta, kh, 6.0_wp, 3.0_wp, Tin, Tsn, nstep, ok, ierr)
         call check(error, ok .and. ierr == 0_int32, "marching converged")
         if (allocated(error)) return
         R = hi/c%kice
         if (l_snow) R = R + hs/c%ksno
         F = (Ta - Tbot)/(1.0_wp/lam + R)
         worst = abs(Tsf - (Ta - F/lam))
         if (l_snow) then
            if (.not. (abs(Tsn - (Tsf - F*hs/(2.0_wp*c%ksno))) <= worst)) &
               worst = abs(Tsn - (Tsf - F*hs/(2.0_wp*c%ksno)))
         end if
         do j = 1, NIL
            R = (real(j, wp) - 0.5_wp)*hl/c%kice
            if (l_snow) R = R + hs/c%ksno
            if (.not. (abs(Tin(j) - (Tsf - F*R)) <= worst)) worst = abs(Tin(j) - (Tsf - F*R))
         end do
         call check(error, worst <= 1.0e-13_wp, "cold-surface steady state")
         if (allocated(error)) return
         ! the surface balance closes: lam*(Ta - Tsf) == conductive flux
         call check(error, abs(lam*(Ta - Tsf) - F) <= 1.0e-13_wp*abs(F), "surface balance closes at the solution")
         if (allocated(error)) return
      end do
   end subroutine test_a3_cold

   ! ------------------------------------------------------------------
   ! vertical profile
   ! ------------------------------------------------------------------

   subroutine test_vp_closed(error)
      !! Thicknesses, temperatures, melting points and the initial energy of
      !! hand-built columns, both constant sets.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer :: iset, k
      integer, parameter :: NIL = 4
      real(wp) :: aicen, vicen, vsnon, zq(NIL), zS(NIL), zT(NIL), zqs(1), zTs(1), Tm(NIL)
      real(wp) :: Texp(NIL), Tsexp, hin, hl, hsn, hsl, einit, eexp, worst
      integer(int32) :: ierr

      do iset = 0, 1
         c = ice_constants(merge(ICE_CONSTANTS_CESM, ICE_CONSTANTS_ICEPACK, iset == 1))
         aicen = 0.5_wp
         vicen = 1.0_wp
         vsnon = 0.1_wp
         zS = [0.2_wp, 1.5_wp, 2.7_wp, 3.2_wp]
         Texp = [-15.0_wp, -9.0_wp, -4.0_wp, -2.0_wp]
         Tsexp = -12.0_wp
         eexp = 0.0_wp
         do k = 1, NIL
            zq(k) = q_fwd(Texp(k), zS(k), c)
         end do
         zqs(1) = qs_fwd(Tsexp, c)
         ierr = 0_int32
         call ice_init_vertical_profile(aicen, vicen, vsnon, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                        hin, hl, hsn, hsl, einit, ierr)
         call check(error, ierr == 0_int32, "no error")
         if (allocated(error)) return
         call check(error, hin == 2.0_wp .and. abs(hl - 0.5_wp) <= 1.0e-15_wp .and. abs(hsn - 0.2_wp) <= 1.0e-15_wp &
                    .and. abs(hsl - 0.2_wp) <= 1.0e-15_wp, "thicknesses")
         if (allocated(error)) return
         worst = 0.0_wp
         do k = 1, NIL
            if (.not. (abs(zT(k) - Texp(k)) <= worst)) worst = abs(zT(k) - Texp(k))
            if (.not. (abs(Tm(k) + c%depressT*zS(k)) <= worst)) worst = abs(Tm(k) + c%depressT*zS(k))
            eexp = eexp + hl*q_fwd(Texp(k), zS(k), c)
         end do
         if (.not. (abs(zTs(1) - Tsexp) <= worst)) worst = abs(zTs(1) - Tsexp)
         call check(error, worst <= 2.0e-13_wp, "temperatures round-trip and Tmlts = -depressT*S")
         if (allocated(error)) return
         eexp = eexp + hsl*qs_fwd(Tsexp, c)
         call check(error, abs(einit - eexp) <= 1.0e-12_wp*abs(eexp), "einit = sum h q")
         if (allocated(error)) return
         ! no snow: the thin-snow reset makes the snow layer the melting reference
         zqs(1) = 5.0e6_wp
         call ice_init_vertical_profile(aicen, vicen, 0.0_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                        hin, hl, hsn, hsl, einit, ierr)
         call check(error, zqs(1) == -c%rhos*c%Lfresh .and. zTs(1) == 0.0_wp .and. hsn == 0.0_wp .and. hsl == 0.0_wp, &
                    "bare ice: snow reset to the reference enthalpy")
         if (allocated(error)) return
      end do
   end subroutine test_vp_closed

   subroutine test_vp_repairs(error)
      !! The round-off repairs upstream applies: a layer at/above its melting
      !! point is put `puny` below it; positive fresh-ice and snow
      !! temperatures within the allowance are set to 0 with the reference
      !! enthalpy; fresh ice uses the linear enthalpy.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: zq(NIL), zS(NIL), zT(NIL), zqs(1), zTs(1), Tm(NIL), hin, hl, hsn, hsl, einit
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      zS = [1.0_wp, 2.0_wp, 3.0_wp, 3.2_wp]
      zq = [q_fwd(-8.0_wp, 1.0_wp, c), q_fwd(-6.0_wp, 2.0_wp, c), q_fwd(-4.0_wp, 3.0_wp, c), q_fwd(-2.0_wp, 3.2_wp, c)]
      ! layer 3 exactly at its melting point: q = rhoi*cp_ocn*Tmlt
      zq(3) = c%rhoi*c%cp_ocn*(-3.0_wp*c%depressT)
      zqs(1) = qs_fwd(0.5e-12_wp, c)      ! snow T = +5e-13, inside the allowance
      ierr = 0_int32
      call ice_init_vertical_profile(0.7_wp, 0.7_wp, 0.07_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == 0_int32, "no error")
      if (allocated(error)) return
      call check(error, abs(zT(3) - (-3.0_wp*c%depressT - c%puny)) <= 1.0e-15_wp, "layer at Tmlt put puny below")
      if (allocated(error)) return
      call check(error, abs(zq(3) - (-3.0_wp*c%depressT*c%rhoi*c%cp_ocn)) <= 1.0e-14_wp*abs(zq(3)), "its enthalpy rewritten")
      if (allocated(error)) return
      call check(error, zTs(1) == 0.0_wp .and. zqs(1) == -c%rhos*c%Lfresh, "positive snow round-off -> 0 and reference q")
      if (allocated(error)) return
      ! fresh ice (no brine pockets): linear enthalpy, T = (Lfresh + q/rhoi)/cp_ice
      zS = 0.1_wp
      zq = c%rhoi*(c%cp_ice*[-20.0_wp, -10.0_wp, -5.0_wp, 0.5e-13_wp] - c%Lfresh)
      zqs(1) = qs_fwd(-3.0_wp, c)
      call ice_init_vertical_profile(0.7_wp, 0.7_wp, 0.07_wp, NIL, 1, c, .false., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == 0_int32, "fresh: no error")
      if (allocated(error)) return
      call check(error, abs(zT(1) + 20.0_wp) <= 1.0e-12_wp .and. abs(zT(3) + 5.0_wp) <= 1.0e-12_wp, "fresh-ice temperatures")
      if (allocated(error)) return
      ! last layer: T = +5e-14 is within the (rhos-based) allowance, set to 0, then (S=0.1) put puny
      ! below the melting point -0.0054
      call check(error, abs(zT(4) - (-0.1_wp*c%depressT - c%puny)) <= 1.0e-15_wp, "fresh round-off layer repaired")
   end subroutine test_vp_repairs

   subroutine test_vp_range(error)
      !! Each documented fatal bit, and that the layer arrays stay untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: zq(NIL), zS(NIL), zT(NIL), zqs(1), zTs(1), Tm(NIL), hin, hl, hsn, hsl, einit
      real(wp) :: zq0(NIL), zqs0(1)
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      zS = [1.0_wp, 2.0_wp, 3.0_wp, 3.2_wp]
      zq = [q_fwd(-8.0_wp, 1.0_wp, c), q_fwd(-6.0_wp, 2.0_wp, c), q_fwd(-4.0_wp, 3.0_wp, c), q_fwd(-2.0_wp, 3.2_wp, c)]
      zqs = qs_fwd(-5.0_wp, c)
      zT = 7.0_wp
      zTs = 7.0_wp
      Tm = 7.0_wp
      zq0 = zq
      zqs0 = zqs

      ! snow too cold
      zqs = qs_fwd(-150.0_wp, c)
      zqs0 = zqs
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, iand(ierr, ICE_ERR_TSN_RANGE) /= 0_int32 .and. ierr == ICE_ERR_TSN_RANGE, "snow below Tmin")
      if (allocated(error)) return
      call check(error, all(zT == 7.0_wp) .and. all(zq == zq0) .and. zqs(1) == zqs0(1) .and. zTs(1) == 7.0_wp, "arrays untouched")
      if (allocated(error)) return
      ! snow too warm: T = +0.5 degC
      zqs = qs_fwd(0.5_wp, c)
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ICE_ERR_TSN_RANGE, "snow above Tmax")
      if (allocated(error)) return
      call check(error, hin == 0.0_wp .and. hl == 0.0_wp .and. hsn == 0.0_wp .and. hsl == 0.0_wp .and. einit == 0.0_wp, &
                 "range error zeroes every scalar output")
      if (allocated(error)) return
      ! salinity below the threshold
      zqs = qs_fwd(-5.0_wp, c)
      zS(2) = 0.05_wp
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ICE_ERR_SIN_LOW, "salinity below min_salin")
      if (allocated(error)) return
      zS(2) = 2.0_wp
      ! ice too cold
      zq(3) = q_fwd(-120.0_wp, 3.0_wp, c)
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ICE_ERR_TIN_RANGE, "ice below Tmin")
      if (allocated(error)) return
      call check(error, hin == 0.0_wp .and. hl == 0.0_wp .and. hsn == 0.0_wp .and. hsl == 0.0_wp .and. einit == 0.0_wp, &
                 "ice range error zeroes every scalar output")
      if (allocated(error)) return
      ! fresh ice well above 0 degC
      zq = c%rhoi*(c%cp_ice*[-20.0_wp, -10.0_wp, 1.0_wp, -4.0_wp] - c%Lfresh)
      zS = 0.1_wp
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .false., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ICE_ERR_TIN_RANGE, "fresh ice above the allowance")
      if (allocated(error)) return
      ! earlier bits are preserved (OR-ed into)
      ierr = ICE_ERR_NONFINITE
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, 0.1_wp, NIL, 1, c, .false., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ior(ICE_ERR_NONFINITE, ICE_ERR_TIN_RANGE), "bits OR-ed, none cleared")
      if (allocated(error)) return
      call check(error, ICE_TMIN == -100.0_wp, "Tmin of record")
      if (allocated(error)) return
      ! round-off negative snow volume is "no snow", not an error (upstream: thin-snow reset)
      zq = [q_fwd(-8.0_wp, 1.0_wp, c), q_fwd(-6.0_wp, 2.0_wp, c), q_fwd(-4.0_wp, 3.0_wp, c), q_fwd(-2.0_wp, 3.2_wp, c)]
      zS = [1.0_wp, 2.0_wp, 3.0_wp, 3.2_wp]
      zqs = 3.0e6_wp
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 1.0_wp, -1.0e-20_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == 0_int32 .and. zqs(1) == -c%rhos*c%Lfresh .and. zTs(1) == 0.0_wp, "negative snow volume = no snow")
      if (allocated(error)) return
      ! non-positive ice volume: NEG_THICK, scalars zero
      ierr = 0_int32
      call ice_init_vertical_profile(0.5_wp, 0.0_wp, 0.0_wp, NIL, 1, c, .true., zS, zq, zqs, zT, zTs, Tm, &
                                     hin, hl, hsn, hsl, einit, ierr)
      call check(error, ierr == ICE_ERR_NEG_THICK .and. hin == 0.0_wp .and. einit == 0.0_wp, "vicen = 0 flagged NEG_THICK")
   end subroutine test_vp_range

   ! ------------------------------------------------------------------
   ! unchecked cores: non-finite input must still end in the error bit
   ! ------------------------------------------------------------------

   subroutine test_core_tridiag_nonfinite(error)
      !! `ice_tridiag_solver_core` has no entry scan: NaN, +Inf and -Inf placed in
      !! `diag(k)`, `sbdiag(k)`, `spdiag(k-1)` and `rhs(k)` (every k > 1), with
      !! and without a coupled super-diagonal, must each raise `ICE_ERR_NONFINITE`
      !! (through the pivot or the solution test).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NM = 5
      real(wp) :: bad(3), sb(NM), dg(NM), sp(NM), rh(NM), x(NM), spbase
      integer(int32) :: ierr
      integer :: ib, k, which, isp
      character(len=64) :: msg

      bad = [ieee_value(1.0_wp, ieee_quiet_nan), ieee_value(1.0_wp, ieee_positive_inf), &
             ieee_value(1.0_wp, ieee_negative_inf)]
      do isp = 1, 2
         spbase = merge(0.0_wp, -0.1_wp, isp == 1)
         do ib = 1, 3
            do k = 2, NM
               do which = 1, 4
                  sb = -0.1_wp; dg = 2.0_wp; sp = spbase; rh = 1.0_wp
                  if (which == 1) dg(k) = bad(ib)
                  if (which == 2) sb(k) = bad(ib)
                  if (which == 3) sp(k - 1) = bad(ib)
                  if (which == 4) rh(k) = bad(ib)
                  x = 0.0_wp
                  ierr = 0_int32
                  call ice_tridiag_solver_core(NM, sb, dg, sp, rh, x, ierr)
                  write (msg, '(a,i0,a,i0,a,i0,a,i0)') "core tridiag flags: bad ", ib, " k ", k, " array ", which, " sp ", isp
                  call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, trim(msg))
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_core_tridiag_nonfinite

   subroutine test_core_chain_nonfinite(error)
      !! The unchecked matrix + Thomas chain flags a NaN `kh(k)` or `etai(k)` and
      !! an Inf `kh(1)` on a snow column, for the cold and the melting surface.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NIL = 4, NM = NIL + 2
      real(wp) :: kh(NM), etai(NIL), Tin(NIL), Isw(NIL), Ts1(1), etas(1), Ssw(1)
      real(wp) :: sb(NM), dg(NM), sp(NM), rh(NM), x(NM), nan, inf
      integer(int32) :: ierr
      integer :: icase, icold
      logical :: l_cold
      character(len=64) :: msg

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      inf = ieee_value(1.0_wp, ieee_positive_inf)
      Tin = -5.0_wp
      Isw = 0.0_wp
      Ts1 = -5.0_wp
      Ssw = 0.0_wp
      etas = 0.1_wp
      do icold = 1, 2
         l_cold = icold == 1
         ! cases 1..NM: NaN kh(k); NM+1..NM+NIL: NaN etai(k); last: Inf kh(1)
         do icase = 1, NM + NIL + 1
            kh = 3.0_wp
            etai = 0.02_wp
            if (icase <= NM) kh(icase) = nan
            if (icase > NM .and. icase <= NM + NIL) etai(icase - NM) = nan
            if (icase == NM + NIL + 1) kh(1) = inf
            sb = 0.0_wp; dg = 1.0_wp; sp = 0.0_wp; rh = 0.0_wp; x = 0.0_wp
            ierr = 0_int32
            call ice_matrix_calc_Tsfc_core(.true., l_cold, NIL, 1, -0.5_wp, -1.8_wp, 1.0_wp, -10.0_wp, Tin, Ts1, kh, Ssw, Isw, &
                                           etai, etas, sb, dg, sp, rh)
            call ice_tridiag_solver_core(NM, sb, dg, sp, rh, x, ierr)
            write (msg, '(a,i0,a,i0)') "core chain flags: case ", icase, " cold ", icold
            call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, trim(msg))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_core_chain_nonfinite

   subroutine test_wrappers_equal_cores(error)
      !! Each public checked routine equals its core bitwise on finite input.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4, NM = NIL + 2
      real(wp) :: T(NIL), S(NIL), khw(NM), khc(NM), Tin(NIL), Isw(NIL), etai(NIL), Ts1(1), Ssw(1), etas(1)
      real(wp) :: sbw(NM), dgw(NM), spw(NM), rhw(NM), sbc(NM), dgc(NM), spc(NM), rhc(NM), xw(NM), xc(NM)
      integer(int32) :: ierr
      integer :: icold, isnow
      logical :: l_cold, l_snow

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      T = [-9.0_wp, -6.5_wp, -4.0_wp, -2.5_wp]
      S = [1.0_wp, 2.0_wp, 3.0_wp, 3.5_wp]
      Tin = T
      Isw = [3.0_wp, 2.0_wp, 1.0_wp, 0.5_wp]
      etai = [0.021_wp, 0.019_wp, 0.023_wp, 0.02_wp]
      Ts1 = -7.0_wp
      Ssw = 0.7_wp
      etas = 0.13_wp
      do isnow = 1, 2
         l_snow = isnow == 1
         do icold = 1, 2
            l_cold = icold == 1
            khw = -1.0_wp
            khc = -1.0_wp
            ierr = 0_int32
            call ice_conductivity(l_snow, NIL, 1, 0.31_wp, 0.07_wp, T, S, c, ICE_CONDUCT_BUBBLY, khw, ierr)
            call ice_conductivity_core(l_snow, NIL, 1, 0.31_wp, 0.07_wp, T, S, c, ICE_CONDUCT_BUBBLY, khc)
            call check(error, ierr == 0_int32 .and. same_bits(khw, khc), "conductivity wrapper == core")
            if (allocated(error)) return
            sbw = 9.0_wp; dgw = 9.0_wp; spw = 9.0_wp; rhw = 9.0_wp
            sbc = 9.0_wp; dgc = 9.0_wp; spc = 9.0_wp; rhc = 9.0_wp
            call ice_matrix_calc_Tsfc(l_snow, l_cold, NIL, 1, -0.4_wp, -1.8_wp, 2.0_wp, -9.0_wp, Tin, Ts1, khw, Ssw, Isw, &
                                      etai, etas, sbw, dgw, spw, rhw, ierr)
            call ice_matrix_calc_Tsfc_core(l_snow, l_cold, NIL, 1, -0.4_wp, -1.8_wp, 2.0_wp, -9.0_wp, Tin, Ts1, khc, Ssw, Isw, &
                                           etai, etas, sbc, dgc, spc, rhc)
            call check(error, ierr == 0_int32 .and. same_bits(sbw, sbc) .and. same_bits(dgw, dgc) .and. same_bits(spw, spc) &
                       .and. same_bits(rhw, rhc), "matrix wrapper == core")
            if (allocated(error)) return
            xw = 0.0_wp
            xc = 0.0_wp
            call ice_tridiag_solver(NM, sbw, dgw, spw, rhw, xw, ierr)
            call ice_tridiag_solver_core(NM, sbc, dgc, spc, rhc, xc, ierr)
            call check(error, ierr == 0_int32 .and. same_bits(xw, xc), "tridiag wrapper == core")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_wrappers_equal_cores

   function same_bits(a, b) result(ok)
      real(wp), intent(in) :: a(:), b(:)
      logical :: ok
      integer :: i
      ok = size(a) == size(b)
      if (.not. ok) return
      do i = 1, size(a)
         if (transfer(a(i), 1_int64) /= transfer(b(i), 1_int64)) then
            ok = .false.
            return
         end if
      end do
   end function same_bits

   ! ------------------------------------------------------------------
   ! NaN-in, device
   ! ------------------------------------------------------------------

   subroutine test_nan(error)
      !! A NaN in any input flags NONFINITE and leaves the outputs untouched
      !! (no NaN laundering through `min`/comparisons).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NIL = 4
      real(wp) :: nan, kh(ICE_NMAT_MAX), kh0(ICE_NMAT_MAX), T(NIL), S(NIL)
      real(wp) :: sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX), sb0(ICE_NMAT_MAX)
      real(wp) :: x(ICE_NMAT_MAX), x0(ICE_NMAT_MAX), Isw(NIL), etai(NIL), Tin(NIL), Ts1(1), zq(NIL), zT(NIL), zqs(1), zTs(1)
      real(wp) :: Tm(NIL), hin, hl, hsn, hsl, einit, ai
      integer(int32) :: ierr
      integer :: i

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      T = -5.0_wp
      S = 1.0_wp
      kh0 = 3.0_wp
      kh = kh0
      ! conductivity: NaN temperature / salinity / thickness, zero thickness
      do i = 1, 4
         T = -5.0_wp
         S = 1.0_wp
         hl = 0.3_wp
         if (i == 1) T(2) = nan
         if (i == 2) S(4) = nan
         if (i == 3) hl = nan
         if (i == 4) hl = 0.0_wp
         ierr = 0_int32
         call ice_conductivity(.false., NIL, 1, hl, 0.0_wp, T, S, c, ICE_CONDUCT_BUBBLY, kh, ierr)
         call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "conductivity flags")
         if (allocated(error)) return
         call check(error, all(kh == kh0), "conductivity leaves kh untouched")
         if (allocated(error)) return
      end do
      ! matrix rows
      Tin = -5.0_wp
      Isw = 0.0_wp
      etai = 0.02_wp
      Ts1 = -5.0_wp
      sb = 1.0_wp; dg = 2.0_wp; sp = 3.0_wp; rh = 4.0_wp
      sb0 = sb
      do i = 1, 4
         kh = 3.0_wp
         if (i == 1) kh(3) = nan
         ierr = 0_int32
         call ice_matrix_calc_Tsfc(.true., .true., NIL, 1, merge(nan, -3.0_wp, i == 2), -1.8_wp, merge(nan, 1.0_wp, i == 3), &
                                   -10.0_wp, merge([nan, 1.0_wp, 1.0_wp, 1.0_wp], Tin, i == 4), Ts1, kh, [0.0_wp], Isw, etai, &
                                   [0.1_wp], sb, dg, sp, rh, ierr)
         call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "matrix flags")
         if (allocated(error)) return
         call check(error, all(sb == 1.0_wp) .and. all(dg == 2.0_wp) .and. all(sp == 3.0_wp) .and. all(rh == 4.0_wp), &
                    "matrix leaves rows untouched")
         if (allocated(error)) return
      end do
      ! tridiagonal solver: NaN rhs, NaN diag, zero pivot
      x0 = 9.0_wp
      do i = 1, 3
         sb = 0.0_wp; dg = 2.0_wp; sp = 0.0_wp; rh = 1.0_wp
         if (i == 1) rh(2) = nan
         if (i == 2) dg(3) = nan
         if (i == 3) dg(1) = 0.0_wp
         x = x0
         ierr = 0_int32
         call ice_tridiag_solver(4, sb, dg, sp, rh, x, ierr)
         call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "tridiag flags")
         if (allocated(error)) return
         call check(error, all(x == x0), "tridiag leaves x untouched")
         if (allocated(error)) return
      end do
      ! vertical profile
      zq = q_fwd(-5.0_wp, 1.0_wp, c)
      zqs = qs_fwd(-5.0_wp, c)
      do i = 1, 5
         S = 1.0_wp
         zq = q_fwd(-5.0_wp, 1.0_wp, c)
         zqs = qs_fwd(-5.0_wp, c)
         zT = 7.0_wp
         zTs = 7.0_wp
         if (i == 1) zq(2) = nan
         if (i == 2) S(3) = nan
         if (i == 3) zqs(1) = nan
         ai = 0.5_wp
         if (i == 4) ai = nan
         if (i == 5) ai = 0.0_wp
         ierr = 0_int32
         call ice_init_vertical_profile(ai, 1.0_wp, 0.1_wp, NIL, 1, c, .true., S, zq, zqs, zT, zTs, Tm, &
                                        hin, hl, hsn, hsl, einit, ierr)
         if (i == 5) then
            call check(error, ierr == ICE_ERR_NEG_THICK, "aicen = 0 flagged NEG_THICK")
         else
            call check(error, ierr == ICE_ERR_NONFINITE, "vprofile flags NONFINITE")
         end if
         if (allocated(error)) return
         call check(error, all(zT == 7.0_wp) .and. zTs(1) == 7.0_wp .and. einit == 0.0_wp .and. hin == 0.0_wp &
                    .and. hl == 0.0_wp .and. hsn == 0.0_wp .and. hsl == 0.0_wp, "vprofile leaves columns untouched, scalars zero")
         if (allocated(error)) return
      end do
      ! layer count out of range
      ierr = 0_int32
      call ice_tridiag_solver(ICE_NMAT_MAX + 1, sb, dg, sp, rh, x, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "nmat out of range flagged")
   end subroutine test_nan

   pure subroutine chain(i, hi, Tsfc, Tl, Sl, c, xo, ierr)
      !! One column through the whole BL99 setup: profile -> conductivity ->
      !! matrix -> solve.  Used by the host reference and the device kernel.
      !$acc routine seq
      integer, intent(in) :: i
      real(wp), intent(in) :: hi, Tsfc
      real(wp), intent(in) :: Tl(4), Sl(4)
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(out) :: xo(6)
      integer(int32), intent(out) :: ierr
      integer :: k
      real(wp) :: zq(4), zT(4), Tm(4), zqs(1), zTs(1), hin, hl, hsn, hsl, einit
      real(wp) :: kh(6), etai(4), Isw(4), sb(6), dg(6), sp(6), rh(6), ci, ai, as, Ssw(1), etas(1)
      logical :: l_snow

      ierr = 0_int32
      xo = 0.0_wp
      do k = 1, 4
         zq(k) = ice_q_from_T(Tl(k), Sl(k), c, .true.)
      end do
      ai = 0.4_wp + 0.5_wp*real(mod(i, 7), wp)/7.0_wp
      l_snow = mod(i, 3) /= 0
      as = merge(0.12_wp, 0.0_wp, l_snow)
      zqs(1) = -c%rhos*(c%Lfresh - c%cp_ice*(-7.0_wp))
      call ice_init_vertical_profile(ai, ai*hi, ai*as, 4, 1, c, .true., Sl, zq, zqs, zT, zTs, Tm, hin, hl, hsn, hsl, einit, ierr)
      if (ierr /= 0_int32) return
      call ice_conductivity(hsl > c%hs_min, 4, 1, hl, hsl, zT, Sl, c, ICE_CONDUCT_BUBBLY, kh, ierr)
      do k = 1, 4
         ci = c%cp_ice - c%Lfresh*Tm(k)/(zT(k)*zT(k))
         etai(k) = 600.0_wp/(c%rhoi*ci*hl)
         Isw(k) = 0.5_wp
      end do
      Ssw(1) = 1.0_wp
      etas(1) = 600.0_wp/(c%rhos*c%cp_ice*max(hsl, c%hs_min))
      call ice_matrix_calc_Tsfc(hsl > c%hs_min, .true., 4, 1, Tsfc, -1.8_wp, 20.0_wp - Tsfc, -9.0_wp, zT, zTs, kh, Ssw, Isw, &
                                etai, etas, sb, dg, sp, rh, ierr)
      call ice_tridiag_solver(6, sb, dg, sp, rh, xo, ierr)
   end subroutine chain

   subroutine test_device(error)
      !! The whole chain called from one `do concurrent` over mapped arrays
      !! (every array device-resident under `mem:separate`, the constants
      !! bundle copied in) matches the host evaluation (measured worst 3.6e-16 on the V100, 0 on the host; tolerance 5e-15: FMA is on
      !! for the GPU build).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: N = 4096
      real(wp) :: hi(N), hs(N), Tsfc(N), S(4, N), T(4, N), xd(6, N), xh(6, N), worst, d
      integer(int32) :: ierr(N), ierr_h(N), ierr_l
      integer :: i

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call ice_random_columns(N, 4, 31337_int64, hi, hs, Tsfc, S, T)
      do i = 1, N
         call chain(i, hi(i), Tsfc(i), T(:, i), S(:, i), c, xh(:, i), ierr_h(i))
      end do
      xd = 0.0_wp
      ierr = 0_int32

      !$acc enter data copyin(c, hi, Tsfc, S, T, ierr) create(xd)
      do concurrent(i=1:N) local(ierr_l)
         call chain(i, hi(i), Tsfc(i), T(:, i), S(:, i), c, xd(:, i), ierr_l)
         ierr(i) = ierr_l
      end do
      !$acc update self(xd, ierr)
      !$acc exit data delete(c, hi, Tsfc, S, T, ierr, xd)

      call check(error, all(ierr == 0_int32) .and. all(ierr_h == 0_int32), "no error bits")
      if (allocated(error)) return
      worst = 0.0_wp
      do i = 1, N
         d = maxval(abs(xd(:, i) - xh(:, i))/max(abs(xh(:, i)), 1.0_wp))
         if (.not. (d <= worst)) worst = d
      end do
      call check(error, worst <= 5.0e-15_wp, "device chain == host chain")
   end subroutine test_device

end module test_cice_bl99_pieces
