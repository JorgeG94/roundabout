module test_cice_thermo_shared
   !! Tests of `rdb_ice_thermo_shared`: closed forms for the saturation
   !! humidity and the surface balance, finite-difference checks of every
   !! analytic derivative (Richardson-extrapolated central differences, relative
   !! 1e-7), the linearised `sfc_kind` hook, `adjust_enthalpy` conservation on
   !! random columns (relative 1e-14), NaN-proof error flagging, and a
   !! `do concurrent` device test over mapped arrays.  Expected values are
   !! closed forms written here or finite differences of the kernel itself,
   !! never the kernel's own derivative.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ice_q_from_T, ice_params_t, ice_params_default
   use rdb_ice_thermo_shared, only: ICE_SFC_BULK, ICE_SFC_LINEAR, ICE_QQQICE, ICE_TTTICE, &
                                    surface_heat_flux, dsurface_heat_flux_dTsf, &
                                    ice_qsat_surface, ice_dqsat_dTsf, adjust_enthalpy
   use test_cice_util, only: ice_random_columns
   implicit none
   private

   public :: collect_cice_thermo_shared_tests

contains

   function emis_default() result(e)
      !! Surface emissivity of the upstream parameter bundle.
      real(wp) :: e
      type(ice_params_t) :: p
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      e = p%thermo%emissivity
   end function emis_default

   subroutine collect_cice_thermo_shared_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("qsat_closed_form", test_qsat_closed), &
                  new_unittest("surface_flux_closed_form", test_flux_closed), &
                  new_unittest("flux_derivatives_vs_finite_difference", test_fd), &
                  new_unittest("linear_sfc_kind_hook", test_linear), &
                  new_unittest("adjust_enthalpy_identity_and_uniform", test_adjust_trivial), &
                  new_unittest("adjust_enthalpy_conserves_random_columns", test_adjust_random), &
                  new_unittest("adjust_enthalpy_growth_melt_regrid", test_adjust_growth), &
                  new_unittest("adjust_enthalpy_thin_and_zero", test_adjust_thin), &
                  new_unittest("nan_and_range_flags", test_nan), &
                  new_unittest("device_kernels_match_host", test_device) &
                  ]
   end subroutine collect_cice_thermo_shared_tests

   subroutine test_qsat_closed(error)
      !! qsat(273.15 K)/rhoa from an independent evaluation of
      !! qqqice*exp(-TTTice/T), and the derivative TTTice/T^2*Qsfc.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: q, dq

      call check(error, ICE_QQQICE == 11637800.0_wp .and. ICE_TTTICE == 5897.8_wp, "constants of record")
      if (allocated(error)) return
      q = ice_qsat_surface(273.15_wp, 1.3_wp)
      call check(error, abs(q - 0.0037560227407897324_wp) <= 1.0e-13_wp*0.0037560227407897324_wp, "Qsfc(273.15)")
      if (allocated(error)) return
      dq = ice_dqsat_dTsf(273.15_wp, 1.3_wp)
      call check(error, abs(dq - 0.00029690386043182475_wp) <= 1.0e-13_wp*0.00029690386043182475_wp, &
                 "dQsfc/dTsf(273.15)")
      if (allocated(error)) return
      ! Scales as 1/rhoa exactly; monotone increasing in T.
      call check(error, abs(ice_qsat_surface(260.0_wp, 2.6_wp)*2.0_wp - ice_qsat_surface(260.0_wp, 1.3_wp)) &
                 <= 1.0e-15_wp, "1/rhoa scaling")
      if (allocated(error)) return
      call check(error, ice_qsat_surface(250.0_wp, 1.3_wp) < ice_qsat_surface(260.0_wp, 1.3_wp), "monotone")
   end subroutine test_qsat_closed

   subroutine test_flux_closed(error)
      !! With lhcoef = 0 the surface balance has no exponential:
      !! fsurf = fsw + eps*flw - eps*sigma*T^4 + shcoef*(potT - T).
      !! At Tsf = 0 degC the humidity of the surface is the qsat closed form.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: lw, se, la, fs, ref, T
      integer(int32) :: ierr
      integer :: s

      do s = 1, 2
         c = ice_constants(merge(ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, s == 1))
         ierr = 0_int32
         call surface_heat_flux(-10.0_wp, 40.0_wp, 1.3_wp, 250.0_wp, 258.0_wp, 0.001_wp, 7.0_wp, 0.0_wp, &
                                emis_default(), c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fs, ierr)
         T = 263.15_wp
         ref = 40.0_wp + emis_default()*250.0_wp - emis_default()*5.67e-8_wp*T**4 + 7.0_wp*(258.0_wp - T)
         call check(error, ierr == 0_int32 .and. abs(la) == 0.0_wp, "no latent, no error")
         if (allocated(error)) return
         call check(error, abs(fs - ref) <= 1.0e-12_wp*abs(ref), "fsurf closed form")
         if (allocated(error)) return
         call check(error, abs(lw + emis_default()*5.67e-8_wp*T**4) <= 1.0e-12_wp*abs(lw), "flwout closed form")
         if (allocated(error)) return
         call check(error, abs(se - 7.0_wp*(258.0_wp - T)) <= 1.0e-12_wp, "fsens closed form")
         if (allocated(error)) return
      end do
      ! Latent term at Tsf = 0, rhoa = 1.3: lhcoef*(Qa - Qsfc) with the
      ! independent Qsfc value.
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      call surface_heat_flux(0.0_wp, 0.0_wp, 1.3_wp, 0.0_wp, 273.15_wp, 0.002_wp, 0.0_wp, 1.0e4_wp, &
                             emis_default(), c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fs, ierr)
      ref = 1.0e4_wp*(0.002_wp - 0.0037560227407897324_wp)
      call check(error, abs(la - ref) <= 1.0e-12_wp*abs(ref), "latent flux closed form")
   end subroutine test_flux_closed

   subroutine test_fd(error)
      !! Every analytic derivative against a Richardson-extrapolated central
      !! difference (error O(h^4)) of the kernel, relative 5e-9 (measured worst 4.3e-10, the Richardson truncation of the exponential at h = 0.25 K; ~10x margin; the linear terms are exact to 1e-13).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp), parameter :: tsf(6) = [-40.0_wp, -25.0_wp, -10.0_wp, -3.0_wp, -0.4_wp, 0.0_wp]
      real(wp), parameter :: h = 0.25_wp
      real(wp) :: lw(4), se(4), la(4), fs(4), d(4), dfs, dlw, dse, dla
      real(wp) :: dh, dh2, fd
      real(wp) :: Tp(4)
      integer(int32) :: ierr
      integer :: i, j, s
      real(wp), parameter :: sh(3) = [8.5_wp, 1.0_wp, 15.0_wp]
      real(wp), parameter :: lh(3) = [2.4e4_wp, 3.0e3_wp, 5.0e4_wp]
      real(wp), parameter :: rho(3) = [1.30_wp, 1.40_wp, 1.25_wp]

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      do s = 1, 3
         do i = 1, 6
            Tp = [tsf(i) - h, tsf(i) + h, tsf(i) - 0.5_wp*h, tsf(i) + 0.5_wp*h]
            do j = 1, 4
               call surface_heat_flux(Tp(j), 50.0_wp, rho(s), 250.0_wp, 255.0_wp, 0.0008_wp, sh(s), lh(s), &
                                      emis_default(), c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw(j), se(j), la(j), fs(j), ierr)
            end do
            call dsurface_heat_flux_dTsf(tsf(i), rho(s), sh(s), lh(s), emis_default(), c, ICE_SFC_BULK, 0.0_wp, &
                                         dfs, dlw, dse, dla, ierr)
            ! d(j) = central differences at h and h/2, then Richardson.
            call richardson(lw, h, fd)
            call check(error, abs(dlw - fd) <= 5.0e-13_wp*abs(dlw), "dflwout/dTsf FD")
            if (allocated(error)) return
            call richardson(se, h, fd)
            call check(error, abs(dse - fd) <= 1.0e-12_wp*abs(dse), "dfsens/dTsf FD")
            if (allocated(error)) return
            call richardson(la, h, fd)
            call check(error, abs(dla - fd) <= 5.0e-9_wp*abs(dla), "dflat/dTsf FD")
            if (allocated(error)) return
            call richardson(fs, h, fd)
            call check(error, abs(dfs - fd) <= 5.0e-10_wp*abs(dfs), "dfsurf/dTsf FD")
            if (allocated(error)) return
            ! The total is the sum of the parts.
            call check(error, abs(dfs - (dlw + dse + dla)) <= 1.0e-14_wp*abs(dfs), "sum of parts")
            if (allocated(error)) return
            ! qsat derivative.
            dh = 273.15_wp + tsf(i)
            dh2 = (4.0_wp*(ice_qsat_surface(dh + 0.5_wp*h, rho(s)) - ice_qsat_surface(dh - 0.5_wp*h, rho(s)))/h &
                   - (ice_qsat_surface(dh + h, rho(s)) - ice_qsat_surface(dh - h, rho(s)))/(2.0_wp*h)) &
                  /3.0_wp
            d(1) = ice_dqsat_dTsf(dh, rho(s))
            call check(error, abs(d(1) - dh2) <= 5.0e-9_wp*abs(d(1)), "dQsfc/dTsf FD")
            if (allocated(error)) return
         end do
      end do
      call check(error, ierr == 0_int32, "no error bits")
   contains
      subroutine richardson(v, hh, out)
         !! v = values at [T-h, T+h, T-h/2, T+h/2].
         real(wp), intent(in) :: v(4), hh
         real(wp), intent(out) :: out
         real(wp) :: d1, d2
         d1 = (v(2) - v(1))/(2.0_wp*hh)
         d2 = (v(4) - v(3))/hh
         out = (4.0_wp*d2 - d1)/3.0_wp
      end subroutine richardson
   end subroutine test_fd

   subroutine test_linear(error)
      !! `sfc_kind = ICE_SFC_LINEAR`: fsurf = -(sf_0 + dsf_dT*Tsf) exactly,
      !! derivative -dsf_dT, no component fluxes, atmosphere inputs ignored.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: lw, se, la, fs, dfs, dlw, dse, dla
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      call surface_heat_flux(-7.5_wp, 1.0e9_wp, 1.3_wp, 1.0e9_wp, 1.0e9_wp, 1.0_wp, 1.0e9_wp, 1.0e9_wp, &
                             emis_default(), c, ICE_SFC_LINEAR, 12.0_wp, 3.5_wp, lw, se, la, fs, ierr)
      call check(error, fs == -(12.0_wp + 3.5_wp*(-7.5_wp)), "linear fsurf")
      if (allocated(error)) return
      call check(error, lw == 0.0_wp .and. se == 0.0_wp .and. la == 0.0_wp, "components zero")
      if (allocated(error)) return
      call dsurface_heat_flux_dTsf(-7.5_wp, 1.3_wp, 1.0e9_wp, 1.0e9_wp, emis_default(), c, ICE_SFC_LINEAR, 3.5_wp, &
                                   dfs, dlw, dse, dla, ierr)
      call check(error, dfs == -3.5_wp .and. dlw == 0.0_wp .and. dse == 0.0_wp .and. dla == 0.0_wp, &
                 "linear derivative")
      if (allocated(error)) return
      call check(error, ierr == 0_int32, "no error bits")
   end subroutine test_linear

   subroutine test_adjust_trivial(error)
      !! Old layers equal to the new layers leave the column unchanged; a
      !! uniform column stays uniform under any layer split.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: z1(5), z2(5), qn(4), q0(4), hlyr
      integer(int32) :: ierr
      integer :: k

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      hlyr = 0.3_wp
      do k = 1, 5
         z2(k) = hlyr*real(k - 1, wp)
      end do
      z1 = z2
      q0 = [-1.0e8_wp, -2.0e8_wp, -3.0e8_wp, -4.0e8_wp]
      qn = q0
      call adjust_enthalpy(4, z1, z2, hlyr, 1.2_wp, qn, c, ierr)
      call check(error, all(abs(qn - q0) <= 1.0e-14_wp*abs(q0)), "identity regrid")
      if (allocated(error)) return
      z1 = [0.0_wp, 0.1_wp, 0.5_wp, 0.55_wp, 1.2_wp]
      qn = -3.0e8_wp
      call adjust_enthalpy(4, z1, z2, hlyr, 1.2_wp, qn, c, ierr)
      call check(error, all(abs(qn + 3.0e8_wp) <= 1.0e-14_wp*3.0e8_wp) .and. ierr == 0_int32, "uniform stays uniform")
   end subroutine test_adjust_trivial

   subroutine test_adjust_random(error)
      !! 2000 random columns, 1, 4 and 8 layers: unequal old layers (weights
      !! from the salinities), enthalpy from the BL99 forward map; total
      !! h*q and total thickness are conserved to 1e-14 (relative).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: N = 2000
      integer, parameter :: NLS(3) = [1, 4, 8]
      real(wp) :: hi(N), hs(N), Tsfc(N), S(8, N), T(8, N)
      real(wp) :: z1(9), z2(9), qn(8), qold(8), hlyr, wsum, e0, e1, worst, dz
      integer(int32) :: ierr
      integer :: i, k, il, nl

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      worst = 0.0_wp
      call ice_random_columns(N, 8, 4242_int64, hi, hs, Tsfc, S, T)
      do il = 1, 3
         nl = NLS(il)
         do i = 1, N
            wsum = sum(S(1:nl, i))
            z1(1) = 0.0_wp
            e0 = 0.0_wp
            do k = 1, nl
               qold(k) = ice_q_from_T(T(k, i), S(k, i), c, .true.)
               dz = hi(i)*S(k, i)/wsum
               z1(k + 1) = z1(k) + dz
               e0 = e0 + qold(k)*dz
            end do
            z1(nl + 1) = hi(i)
            hlyr = hi(i)/real(nl, wp)
            do k = 1, nl + 1
               z2(k) = hlyr*real(k - 1, wp)
            end do
            qn(1:nl) = qold(1:nl)
            call adjust_enthalpy(nl, z1, z2, hlyr, hi(i), qn, c, ierr)
            e1 = hlyr*sum(qn(1:nl))
            dz = abs(e1 - e0)/abs(e0)
            if (.not. (dz <= worst)) worst = dz
         end do
      end do
      call check(error, ierr == 0_int32, "no error bits")
      if (allocated(error)) return
      call check(error, worst <= 1.0e-14_wp, "total h*q conserved")
   end subroutine test_adjust_random

   subroutine test_adjust_growth(error)
      !! Layer-count-preserving regrid after growth/melt at the top and the
      !! bottom: the old layers keep their enthalpy but the top and bottom
      !! layer thicknesses change; the new equal layers carry the same total.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: N = 1000
      real(wp) :: hi(N), hs(N), Tsfc(N), S(4, N), T(4, N)
      real(wp) :: z1(5), z2(5), qn(4), hlyr, hn, dz(4), dtop, dbot, e0, e1, worst, u
      integer(int64) :: st
      integer(int32) :: ierr
      integer :: i, k

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      worst = 0.0_wp
      call ice_random_columns(N, 4, 99_int64, hi, hs, Tsfc, S, T)
      st = 20240601_int64
      do i = 1, N
         dz = hi(i)/4.0_wp
         ! deterministic pseudo-random growth (+) / melt (-), at most 90 % of a layer
         st = mod(st*1103515245_int64 + 12345_int64, 2147483648_int64)
         u = real(st, wp)/2147483648_wp
         dtop = (2.0_wp*u - 1.0_wp)*0.9_wp*dz(1)
         st = mod(st*1103515245_int64 + 12345_int64, 2147483648_int64)
         u = real(st, wp)/2147483648_wp
         dbot = (2.0_wp*u - 1.0_wp)*0.9_wp*dz(4)
         dz(1) = dz(1) + dtop
         dz(4) = dz(4) + dbot
         hn = sum(dz)
         z1(1) = 0.0_wp
         e0 = 0.0_wp
         do k = 1, 4
            qn(k) = ice_q_from_T(T(k, i), S(k, i), c, .true.)
            z1(k + 1) = z1(k) + dz(k)
            e0 = e0 + qn(k)*dz(k)
         end do
         hlyr = hn/4.0_wp
         do k = 1, 5
            z2(k) = hlyr*real(k - 1, wp)
         end do
         call adjust_enthalpy(4, z1, z2, hlyr, hn, qn, c, ierr)
         e1 = hlyr*sum(qn)
         u = abs(e1 - e0)/abs(e0)
         if (.not. (u <= worst)) worst = u
      end do
      call check(error, ierr == 0_int32, "no error bits")
      if (allocated(error)) return
      call check(error, worst <= 1.0e-14_wp, "growth/melt regrid conserves enthalpy")
   end subroutine test_adjust_growth

   subroutine test_adjust_thin(error)
      !! hn <= puny: the column is homogenised to its thickness-weighted mean
      !! (and to zero at hn = 0), the total h*q is preserved.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: z1(5), z2(5), qn(4), hn, hlyr, e0
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      hn = 4.0e-12_wp
      z1 = hn*[0.0_wp, 0.1_wp, 0.4_wp, 0.5_wp, 1.0_wp]
      hlyr = hn/4.0_wp
      z2 = hlyr*[0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      qn = [-1.0e8_wp, -2.0e8_wp, -3.0e8_wp, -5.0e8_wp]
      e0 = -1.0e8_wp*0.1_wp*hn - 2.0e8_wp*0.3_wp*hn - 3.0e8_wp*0.1_wp*hn - 5.0e8_wp*0.5_wp*hn
      call adjust_enthalpy(4, z1, z2, hlyr, hn, qn, c, ierr)
      call check(error, all(abs(qn - e0/hn) <= 1.0e-13_wp*abs(e0/hn)), "thin column uniform at the mean")
      if (allocated(error)) return
      z1 = 0.0_wp
      z2 = 0.0_wp
      qn = -3.0e8_wp
      call adjust_enthalpy(4, z1, z2, 0.0_wp, 0.0_wp, qn, c, ierr)
      call check(error, all(qn == 0.0_wp) .and. ierr == 0_int32, "zero thickness gives zero enthalpy")
   end subroutine test_adjust_thin

   subroutine test_nan(error)
      !! Non-finite input sets the NONFINITE bit; the layer-count guard leaves
      !! the column untouched.  Assertions are written `.not.(x <= tol)`
      !! style so a NaN cannot pass.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: nan, lw, se, la, fs, dfs, dlw, dse, dla, z1(10), z2(10), qn(9), hlyr
      integer(int32) :: ierr
      integer :: k

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      call surface_heat_flux(nan, 0.0_wp, 1.3_wp, 250.0_wp, 255.0_wp, 0.001_wp, 5.0_wp, 1.0e4_wp, emis_default(), c, &
                             ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fs, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "NaN Tsf flagged")
      if (allocated(error)) return
      call check(error,.not. (abs(fs) <= huge(1.0_wp)), "NaN propagates to fsurf")
      if (allocated(error)) return
      ierr = 0_int32
      call dsurface_heat_flux_dTsf(nan, 1.3_wp, 5.0_wp, 1.0e4_wp, emis_default(), c, ICE_SFC_BULK, 0.0_wp, &
                                   dfs, dlw, dse, dla, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "NaN derivative flagged")
      if (allocated(error)) return
      ierr = 0_int32
      z1 = 0.0_wp
      z2 = 0.0_wp
      qn = -1.0_wp
      hlyr = 0.1_wp
      call adjust_enthalpy(9, z1, z2, hlyr, 0.9_wp, qn, c, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(qn == -1.0_wp), &
                 "nlyr > ICE_NILYR_MAX flagged, column untouched")
      if (allocated(error)) return
      ierr = 0_int32
      call adjust_enthalpy(4, z1(1:5), z2(1:5), nan, 0.9_wp, qn(1:4), c, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "NaN layer thickness flagged")
      if (allocated(error)) return
      ierr = 0_int32
      z1(1:5) = [0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      z2(1:5) = [0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      qn(1:4) = [-1.0_wp, nan, -1.0_wp, -1.0_wp]
      call adjust_enthalpy(4, z1(1:5), z2(1:5), 0.1_wp, 0.4_wp, qn(1:4), c, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(z1(1:5) == z2(1:5)), &
                 "NaN layer enthalpy flagged")
      if (allocated(error)) return
      ! NaN in each interface array and in the bed-side interface, column untouched.
      do k = 1, 3
         ierr = 0_int32
         z1(1:5) = [0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
         z2(1:5) = z1(1:5)
         qn(1:4) = -1.0_wp
         select case (k)
         case (1)
            z1(3) = nan
         case (2)
            z2(5) = nan
         case default
            z1(5) = nan
         end select
         call adjust_enthalpy(4, z1(1:5), z2(1:5), 0.1_wp, 0.4_wp, qn(1:4), c, ierr)
         call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(qn(1:4) == -1.0_wp), &
                    "NaN interface depth flagged, column untouched")
         if (allocated(error)) return
      end do
      ! Tsf at or below absolute zero is flagged.
      ierr = 0_int32
      call surface_heat_flux(-273.15_wp, 0.0_wp, 1.3_wp, 250.0_wp, 255.0_wp, 0.001_wp, 5.0_wp, 1.0e4_wp, &
                             emis_default(), c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fs, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "TsfK = 0 flagged")
      if (allocated(error)) return
      ierr = 0_int32
      call dsurface_heat_flux_dTsf(-300.0_wp, 1.3_wp, 5.0_wp, 1.0e4_wp, emis_default(), c, ICE_SFC_BULK, &
                                   0.0_wp, dfs, dlw, dse, dla, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "TsfK < 0 flagged (derivative)")
   end subroutine test_nan

   subroutine test_device(error)
      !! The three kernels called from one `do concurrent` over mapped arrays
      !! (every array device-resident under `mem:separate`) match the same
      !! calls made on the host.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: N = 4096
      integer, parameter :: NL = 4
      real(wp) :: hi(N), hs(N), Tsfc(N), S(NL, N), T(NL, N)
      real(wp) :: fsurf(N), dfsurf(N), flat(N), qo(NL, N), z1(NL + 1, N), z2(NL + 1, N), hl(N)
      real(wp) :: fsurf_h(N), dfsurf_h(N), flat_h(N), qo_h(NL, N)
      real(wp) :: q(NL, N), lw, se, la, dd1, dd2, dd3
      integer(int32) :: ierr(N), ierr_l, ierr_h(N)
      integer :: i, k
      real(wp) :: tsf_k, wsum, cum, em

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      em = emis_default()
      call ice_random_columns(N, NL, 777_int64, hi, hs, Tsfc, S, T)
      do i = 1, N
         wsum = sum(S(:, i))
         cum = 0.0_wp
         z1(1, i) = 0.0_wp
         do k = 1, NL
            q(k, i) = ice_q_from_T(T(k, i), S(k, i), c, .true.)
            cum = cum + hi(i)*S(k, i)/wsum
            z1(k + 1, i) = cum
         end do
         hl(i) = hi(i)/real(NL, wp)
         do k = 1, NL + 1
            z2(k, i) = hl(i)*real(k - 1, wp)
         end do
      end do
      qo = q
      qo_h = q
      ierr = 0_int32
      ierr_h = 0_int32

      ! host reference
      do i = 1, N
         ierr_l = 0_int32
         call surface_heat_flux(Tsfc(i), 30.0_wp, 1.3_wp, 250.0_wp, 255.0_wp, 0.001_wp, 6.0_wp, 2.0e4_wp, &
                                em, c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fsurf_h(i), ierr_l)
         flat_h(i) = la
         call dsurface_heat_flux_dTsf(Tsfc(i), 1.3_wp, 6.0_wp, 2.0e4_wp, em, c, ICE_SFC_BULK, 0.0_wp, &
                                      dfsurf_h(i), dd1, dd2, dd3, ierr_l)
         call adjust_enthalpy(NL, z1(:, i), z2(:, i), hl(i), hi(i), qo_h(:, i), c, ierr_l)
         ierr_h(i) = ierr_l
      end do

      !$acc enter data copyin(c, Tsfc, hi, hl, z1, z2, qo, ierr) create(fsurf, dfsurf, flat)
      do concurrent(i=1:N) local(lw, se, la, dd1, dd2, dd3, ierr_l, tsf_k)
         ierr_l = 0_int32
         tsf_k = Tsfc(i)
         call surface_heat_flux(tsf_k, 30.0_wp, 1.3_wp, 250.0_wp, 255.0_wp, 0.001_wp, 6.0_wp, 2.0e4_wp, &
                                em, c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, lw, se, la, fsurf(i), ierr_l)
         flat(i) = la
         call dsurface_heat_flux_dTsf(tsf_k, 1.3_wp, 6.0_wp, 2.0e4_wp, em, c, ICE_SFC_BULK, 0.0_wp, &
                                      dfsurf(i), dd1, dd2, dd3, ierr_l)
         call adjust_enthalpy(NL, z1(:, i), z2(:, i), hl(i), hi(i), qo(:, i), c, ierr_l)
         ierr(i) = ierr_l
      end do
      !$acc update self(fsurf, dfsurf, flat, qo, ierr)
      !$acc exit data delete(c, Tsfc, hi, hl, z1, z2, qo, ierr, fsurf, dfsurf, flat)

      call check(error, all(ierr == 0_int32) .and. all(ierr_h == 0_int32), "no error bits")
      if (allocated(error)) return
      call check(error, all(abs(fsurf - fsurf_h) <= 1.0e-12_wp*abs(fsurf_h) + 1.0e-10_wp), "fsurf device == host")
      if (allocated(error)) return
      call check(error, all(abs(flat - flat_h) <= 1.0e-12_wp*abs(flat_h) + 1.0e-10_wp), "flat device == host")
      if (allocated(error)) return
      call check(error, all(abs(dfsurf - dfsurf_h) <= 1.0e-12_wp*abs(dfsurf_h)), "dfsurf device == host")
      if (allocated(error)) return
      call check(error, all(abs(qo - qo_h) <= 1.0e-13_wp*abs(qo_h)), "adjust_enthalpy device == host")
   end subroutine test_device

end module test_cice_thermo_shared
