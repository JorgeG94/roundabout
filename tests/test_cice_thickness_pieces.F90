module test_cice_thickness_pieces
   !! Tests of `rdb_ice_thickness`: the freeboard worked case (the interface
   !! lands on the waterline, `wk1 = 0`), mass / enthalpy conservation of the
   !! snow-ice conversion on random columns for both constant sets together with
   !! the fresh-water and salt convention (zero fresh, salt booked from the
   !! total ice-thickness change), the bitwise no-op at positive freeboard,
   !! relayering conservation (relative 1e-14) for nilyr in {1, 2, 4, 7, 8},
   !! the state update and its vanishing cases, NaN flagging with the column
   !! untouched, and a `do concurrent` device test over mapped arrays.
   !! Expected values are closed forms written here.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite, ieee_is_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ice_params_t, ice_params_default
   use rdb_ice_thickness, only: ice_freeboard, ice_relayer, ice_update_state_vthermo
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_thickness_pieces_tests

contains

   subroutine collect_cice_thickness_pieces_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("freeboard_worked_case", test_worked), &
                  new_unittest("freeboard_conserves_mass_enthalpy_salt", test_fb_conserve), &
                  new_unittest("freeboard_noop_when_freeboard_positive", test_fb_noop), &
                  new_unittest("relayer_conserves_enthalpy", test_relayer), &
                  new_unittest("relayer_uniform_thin_and_empty", test_relayer_edge), &
                  new_unittest("update_state_roundtrip_and_vanishing", test_update), &
                  new_unittest("nan_flags_leave_column_untouched", test_nan), &
                  new_unittest("device_kernels_match_host", test_device) &
                  ]
   end subroutine collect_cice_thickness_pieces_tests

   function urand(st) result(u)
      !! Next uniform number in (0,1) from the deterministic LCG state `st`.
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function urand

   subroutine random_column(st, n, hin, hsn, dzi, zqin, zqsn, c)
      !! Unequal ice layers summing to `hin` in [0.05, 3] m, snow `hsn` in
      !! [0, 1.2] m (a third bare), BL99-magnitude enthalpies.
      integer(int64), intent(inout) :: st
      integer, intent(in) :: n
      real(wp), intent(out) :: hin, hsn, dzi(n), zqin(n), zqsn
      type(ice_constants_t), intent(in) :: c
      real(wp) :: w(n), sw, Tsn
      integer :: k

      hin = 0.05_wp + 2.95_wp*urand(st)
      hsn = 0.0_wp
      if (urand(st) > 0.33_wp) hsn = 1.2_wp*urand(st)
      sw = 0.0_wp
      do k = 1, n
         w(k) = 0.2_wp + 1.6_wp*urand(st)
         sw = sw + w(k)
      end do
      do k = 1, n
         dzi(k) = hin*w(k)/sw
         zqin(k) = -(1.5e8_wp + 1.5e8_wp*urand(st))
      end do
      Tsn = -20.0_wp*urand(st)
      zqsn = -c%rhos*(c%Lfresh - c%cp_ice*Tsn)
   end subroutine random_column

   subroutine test_worked(error)
      !! Plan worked case (rhow 1026, rhoi 917, rhos 330): hin = hsn = 0.5 m,
      !! four equal layers.  `wk1 = 0.5 - 0.5*109/330`, `dhsn = wk1*917/1026`,
      !! `dhin = dhsn*330/917`; afterwards the waterline residual is zero.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: hin, hsn, dzi(4), dzs(1), zqin(4), zqsn(1), snoice, dsnow
      real(wp) :: wk1, dhsn, dhin, wk_after
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      hin = 0.5_wp
      hsn = 0.5_wp
      dzi = 0.125_wp
      dzs = 0.5_wp
      zqin = [-2.7e8_wp, -1.8e8_wp, -2.6e8_wp, -1.9e8_wp]
      zqsn = -1.1e8_wp
      snoice = 0.0_wp
      dsnow = 0.0_wp
      ierr = 0_int32
      wk1 = 0.5_wp - 0.5_wp*109.0_wp/330.0_wp
      dhsn = wk1*917.0_wp/1026.0_wp
      dhin = dhsn*330.0_wp/917.0_wp
      call ice_freeboard(4, 1, c, hin, hsn, zqin, zqsn, dzi, dzs, snoice, dsnow, ierr)
      call check(error, ierr == 0_int32, "no error bits")
      if (allocated(error)) return
      call check(error, abs(wk1 - 0.33485_wp) < 1.0e-5_wp .and. abs(dhsn - 0.29927_wp) < 1.0e-5_wp &
                 .and. abs(dhin - 0.10770_wp) < 1.0e-5_wp, "worked-case intermediates match the plan")
      if (allocated(error)) return
      call check(error, abs(hin - 0.6077_wp) < 1.0e-4_wp .and. abs(hsn - 0.2007_wp) < 1.0e-4_wp, "final thicknesses")
      if (allocated(error)) return
      call check(error, abs(hin - (0.5_wp + dhin)) <= 1.0e-15_wp .and. abs(hsn - (0.5_wp - dhsn)) <= 1.0e-15_wp, &
                 "thicknesses from the closed form")
      if (allocated(error)) return
      wk_after = hsn - hin*(1026.0_wp - 917.0_wp)/330.0_wp
      call check(error, abs(wk_after) <= 1.0e-15_wp, "waterline residual wk1 = 0 afterwards")
      if (allocated(error)) return
      call check(error, abs(snoice - dhin) <= 1.0e-15_wp .and. abs(dsnow + dhsn) <= 1.0e-15_wp, "snoice, dsnow")
      if (allocated(error)) return
      call check(error, abs(dzi(1) - (0.125_wp + dhin)) <= 1.0e-15_wp .and. all(dzi(2:4) == 0.125_wp), &
                 "only the top ice layer grew")
      if (allocated(error)) return
      call check(error, abs(zqin(1) - (0.125_wp*(-2.7e8_wp) + dhsn*(-1.1e8_wp))/(0.125_wp + dhin)) &
                 <= 1.0e-6_wp .and. all(zqin(2:4) == [-1.8e8_wp, -2.6e8_wp, -1.9e8_wp]), &
                 "top-layer enthalpy is the thickness-weighted mix with the snow")
   end subroutine test_worked

   subroutine test_fb_conserve(error)
      !! 3000 random columns per constant set.  Snow mass lost equals ice mass
      !! gained (relative 1e-14), the column enthalpy sum(dz*q) is unchanged
      !! (relative 1e-14), the waterline residual vanishes, `hin = sum(dzi)`,
      !! and the ocean booking from the total thickness change is zero fresh
      !! water and `fsalt = -rhoi*dhi*Sref*1e-3/dt` with the salt ledger
      !! `fsalt*dt + d(rhoi*hin*Sref*1e-3) = 0` closed.  The salt-ledger line
      !! is algebraic (it restates the booking formula against the thickness
      !! change): it documents the M-D3-3 convention for the next PR and tests
      !! nothing in the kernel; the fresh-water line does check that the
      !! conversion changes no total mass.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(ice_params_t) :: p
      integer, parameter :: NC = 3000
      integer, parameter :: NLS(5) = [1, 2, 4, 7, 8]
      integer :: iset, ic, n, nconv, k
      integer(int32) :: ierr
      integer(int64) :: st
      real(wp) :: hin, hsn, dzi(8), dzs(1), zqin(8), zqsn(1), snoice, dsnow, hin0, hsn0
      real(wp) :: e0, e1, dt, fresh, fsalt, sref, w_mass, w_e, w_wk, w_salt, w_fresh, w_sum, x
      integer :: ierr_any

      dt = 3600.0_wp
      ierr_any = 0
      nconv = 0
      w_mass = 0.0_wp; w_e = 0.0_wp; w_wk = 0.0_wp; w_salt = 0.0_wp; w_fresh = 0.0_wp; w_sum = 0.0_wp
      do iset = 1, 2
         if (iset == 1) c = ice_constants(ICE_CONSTANTS_ICEPACK)
         if (iset == 2) c = ice_constants(ICE_CONSTANTS_CESM)
         p = ice_params_default(merge(ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, iset == 1), 5, 4, 1)
         sref = p%thermo%ice_ref_salinity
         st = 90210_int64 + iset
         do ic = 1, NC
            n = NLS(1 + mod(ic, 5))
            call random_column(st, n, hin, hsn, dzi(1:n), zqin(1:n), zqsn(1), c)
            hin0 = hin
            hsn0 = hsn
            dzs(1) = hsn
            e0 = sum(dzi(1:n)*zqin(1:n)) + dzs(1)*zqsn(1)
            snoice = 0.0_wp
            dsnow = 0.0_wp
            ierr = 0_int32
            call ice_freeboard(n, 1, c, hin, hsn, zqin(1:n), zqsn, dzi(1:n), dzs, snoice, dsnow, ierr)
            if (ierr /= 0_int32) ierr_any = ierr_any + 1
            if (.not. (ieee_is_finite(hin) .and. ieee_is_finite(hsn) .and. ieee_is_finite(zqin(1)))) then
               call check(error, .false., "non-finite result on finite input")
               return
            end if
            if (snoice <= 0.0_wp) cycle
            nconv = nconv + 1
            e1 = sum(dzi(1:n)*zqin(1:n)) + dzs(1)*zqsn(1)
            x = abs(c%rhos*(hsn0 - hsn) - c%rhoi*(hin - hin0))/(c%rhoi*hin)
            if (.not. ieee_is_finite(x)) then
               call check(error, .false., "mass residual not finite")
               return
            end if
            w_mass = max(w_mass, x)
            x = abs(e1 - e0)/abs(e0)
            if (.not. ieee_is_finite(x)) then
               call check(error, .false., "energy residual not finite")
               return
            end if
            w_e = max(w_e, x)
            x = abs(hsn - hin*(c%rhow - c%rhoi)/c%rhos)
            if (.not. ieee_is_finite(x)) then
               call check(error, .false., "waterline residual not finite")
               return
            end if
            w_wk = max(w_wk, x)
            w_sum = max(w_sum, abs(hin - sum(dzi(1:n)))/hin)
            ! ocean booking from the total thickness change (caller's convention)
            fresh = -(c%rhoi*(hin - hin0) + c%rhos*(hsn - hsn0))/dt
            fsalt = -c%rhoi*(hin - hin0)*sref*1.0e-3_wp/dt
            w_fresh = max(w_fresh, abs(fresh)*dt/(c%rhoi*hin))
            w_salt = max(w_salt, abs(fsalt*dt + (c%rhoi*hin*sref*1.0e-3_wp - c%rhoi*hin0*sref*1.0e-3_wp)) &
                         /(c%rhoi*hin*sref*1.0e-3_wp))
            do k = 1, n
               if (dzi(k) <= 0.0_wp) then
                  call check(error, .false., "non-positive layer")
                  return
               end if
            end do
         end do
      end do
      call check(error, ierr_any == 0, "no error bits on finite input")
      if (allocated(error)) return
      call check(error, nconv > 2000, "enough converting columns in the sample")
      if (allocated(error)) return
      call check(error, w_mass <= 1.0e-14_wp, "snow mass lost == ice mass gained")
      if (allocated(error)) return
      call check(error, w_e <= 1.0e-14_wp, "column enthalpy conserved")
      if (allocated(error)) return
      call check(error, w_wk <= 1.0e-14_wp, "interface on the waterline (m)")
      if (allocated(error)) return
      call check(error, w_sum <= 1.0e-15_wp, "hin == sum(dzi)")
      if (allocated(error)) return
      call check(error, w_fresh <= 1.0e-14_wp, "zero fresh-water booking")
      if (allocated(error)) return
      call check(error, w_salt <= 1.0e-14_wp, "salt ledger closes")
   end subroutine test_fb_conserve

   subroutine test_fb_noop(error)
      !! Positive freeboard, no snow, snow below `puny`, and a tolerance-edge
      !! case leave every output bitwise unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: hin, hsn, dzi(4), dzs(1), zqin(4), zqsn(1), snoice, dsnow
      real(wp) :: hin0, hsn0, dzi0(4), dzs0(1), zqin0(4), snoice0, dsnow0
      integer(int32) :: ierr
      integer :: ic
      real(wp), parameter :: HS(4) = [0.1_wp, 0.0_wp, 1.0e-14_wp, 0.4_wp]
      real(wp), parameter :: HI(4) = [2.0_wp, 1.0_wp, 1.0e-9_wp, 3.0_wp]

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do ic = 1, 4
         hin = HI(ic); hsn = HS(ic)
         dzi = hin/4.0_wp
         dzs = hsn
         zqin = [-2.7e8_wp, -1.8e8_wp, -2.6e8_wp, -1.9e8_wp]
         zqsn = -1.1e8_wp
         snoice = 0.013_wp
         dsnow = -0.002_wp
         hin0 = hin; hsn0 = hsn; dzi0 = dzi; dzs0 = dzs; zqin0 = zqin; snoice0 = snoice; dsnow0 = dsnow
         ierr = 0_int32
         call ice_freeboard(4, 1, c, hin, hsn, zqin, zqsn, dzi, dzs, snoice, dsnow, ierr)
         call check(error, hin == hin0 .and. hsn == hsn0 .and. all(dzi == dzi0) .and. all(dzs == dzs0) &
                    .and. all(zqin == zqin0) .and. snoice == snoice0 .and. dsnow == dsnow0 .and. ierr == 0_int32, &
                    "freeboard is a bitwise no-op")
         if (allocated(error)) return
      end do
   end subroutine test_fb_noop

   subroutine test_relayer(error)
      !! 4000 random columns, nilyr in {1, 2, 4, 7, 8}, both constant sets:
      !! sum(hilyr*q) is conserved to 1e-14 relative, hilyr = hin/nilyr and the
      !! layers come out equal; a uniform input is a fixed point (1e-14).
      !! nslyr = 1: the snow enthalpy is untouched, `hslyr = hsn`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NC = 4000
      integer, parameter :: NLS(5) = [1, 2, 4, 7, 8]
      integer :: iset, ic, n, k, nbad
      integer(int32) :: ierr
      integer(int64) :: st
      real(wp) :: hin, hsn, dzi(8), dzs(1), zqin(8), zqsn(1), zq0(8), hilyr, hslyr, e0, e1, x, w_e, w_u, hin0, zsn0
      real(wp) :: dzeq(8), qeq(8)

      w_e = 0.0_wp
      w_u = 0.0_wp
      nbad = 0
      do iset = 1, 2
         c = ice_constants(merge(ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, iset == 1))
         st = 31337_int64 + iset
         do ic = 1, NC
            n = NLS(1 + mod(ic, 5))
            call random_column(st, n, hin, hsn, dzi(1:n), zqin(1:n), zqsn(1), c)
            hin0 = hin
            zq0(1:n) = zqin(1:n)
            zsn0 = zqsn(1)
            dzs(1) = hsn
            e0 = sum(dzi(1:n)*zqin(1:n))
            ierr = 0_int32
            call ice_relayer(n, 1, c, dzi(1:n), dzs, hin, hsn, zqin(1:n), zqsn, hilyr, hslyr, ierr)
            if (ierr /= 0_int32) nbad = nbad + 1
            e1 = hilyr*sum(zqin(1:n))
            x = abs(e1 - e0)/abs(e0)
            if (.not. ieee_is_finite(x)) then
               call check(error, .false., "relayer residual not finite")
               return
            end if
            w_e = max(w_e, x)
            if (abs(hilyr - hin0/real(n, wp)) > 1.0e-16_wp .or. hin /= hin0 &
                .or. any(abs(dzi(1:n) - hilyr) > 0.0_wp) .or. hslyr /= hsn .or. zqsn(1) /= zsn0) then
               call check(error, .false., "layer thicknesses / snow pass-through")
               return
            end if
            ! uniform input: fixed point
            dzeq(1:n) = hin0/real(n, wp)
            do k = 1, n
               qeq(k) = zq0(1)
            end do
            call ice_relayer(n, 1, c, dzeq(1:n), dzs, hin, hsn, qeq(1:n), zqsn, hilyr, hslyr, ierr)
            x = maxval(abs(qeq(1:n) - zq0(1)))/abs(zq0(1))
            if (.not. ieee_is_finite(x)) then
               call check(error, .false., "fixed-point residual not finite")
               return
            end if
            w_u = max(w_u, x)
         end do
      end do
      call check(error, nbad == 0, "no error bits")
      if (allocated(error)) return
      call check(error, w_e <= 1.0e-14_wp, "sum(h*q) conserved to 1e-14")
      if (allocated(error)) return
      call check(error, w_u <= 1.0e-14_wp, "uniform column is a fixed point")
   end subroutine test_relayer

   subroutine test_relayer_edge(error)
      !! `hin <= puny`: homogenised to the thickness-weighted mean; `hin = 0`:
      !! zero enthalpy and `hilyr = 0`; a negative remnant is clamped to 0.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: hin, hsn, dzi(4), dzs(1), zqin(4), zqsn(1), hilyr, hslyr, mean
      integer(int32) :: ierr

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ierr = 0_int32
      ! thin: 1e-12 m total, unequal layers
      dzi = [4.0e-13_wp, 3.0e-13_wp, 2.0e-13_wp, 1.0e-13_wp]
      hin = 1.0e-12_wp
      zqin = [-1.0e8_wp, -2.0e8_wp, -3.0e8_wp, -4.0e8_wp]
      mean = (4.0e-13_wp*(-1.0e8_wp) + 3.0e-13_wp*(-2.0e8_wp) + 2.0e-13_wp*(-3.0e8_wp) + 1.0e-13_wp*(-4.0e8_wp))/1.0e-12_wp
      hsn = 0.0_wp
      dzs = 0.0_wp
      zqsn = -1.0e8_wp
      call ice_relayer(4, 1, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)
      call check(error, all(abs(zqin - mean) <= 1.0e-12_wp*abs(mean)) .and. abs(hilyr - 2.5e-13_wp) <= 1.0e-27_wp, &
                 "thin column homogenised to the weighted mean")
      if (allocated(error)) return
      ! empty
      dzi = 0.0_wp
      hin = 0.0_wp
      zqin = -1.0e8_wp
      call ice_relayer(4, 1, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)
      call check(error, all(zqin == 0.0_wp) .and. hilyr == 0.0_wp .and. hin == 0.0_wp .and. hslyr == 0.0_wp, &
                 "empty column: zero enthalpy and thickness")
      if (allocated(error)) return
      ! negative remnant
      dzi = 0.0_wp
      dzi(1) = -1.0e-20_wp
      hin = -1.0e-20_wp
      hsn = -1.0e-20_wp
      zqin = -1.0e8_wp
      call ice_relayer(4, 1, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)
      call check(error, hin == 0.0_wp .and. hsn == 0.0_wp .and. hilyr == 0.0_wp .and. all(zqin == 0.0_wp) &
                 .and. ierr == 0_int32, "negative remnants clamp to exactly zero")
   end subroutine test_relayer_edge

   subroutine test_update(error)
      !! Round trip `vicen/aicen = hin`, `vsnon/aicen = hsn`; ice gone: area,
      !! volumes and every enthalpy are zero and `Tsf = Tf`; snow gone alone:
      !! `vsnon = 0` with the snow enthalpy and `Tsf` left alone (upstream).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: zqin(4), zqsn(1), Tsf, aicen, vicen, vsnon
      integer(int32) :: ierr

      ierr = 0_int32
      zqin = [-2.0e8_wp, -2.1e8_wp, -2.2e8_wp, -2.3e8_wp]
      zqsn = -1.1e8_wp
      Tsf = -7.0_wp
      aicen = 0.6_wp
      vicen = 99.0_wp
      vsnon = 99.0_wp
      call ice_update_state_vthermo(4, 1, -1.8_wp, 1.5_wp, 0.25_wp, Tsf, zqin, zqsn, aicen, vicen, vsnon, ierr)
      call check(error, abs(vicen/aicen - 1.5_wp) <= 2.0e-16_wp*1.5_wp .and. abs(vsnon/aicen - 0.25_wp) <= 1.0e-16_wp &
                 .and. aicen == 0.6_wp .and. Tsf == -7.0_wp .and. zqin(1) == -2.0e8_wp .and. ierr == 0_int32, &
                 "round trip, state untouched")
      if (allocated(error)) return
      ! snow only vanishes
      call ice_update_state_vthermo(4, 1, -1.8_wp, 1.5_wp, 0.0_wp, Tsf, zqin, zqsn, aicen, vicen, vsnon, ierr)
      call check(error, vsnon == 0.0_wp .and. zqsn(1) == -1.1e8_wp .and. Tsf == -7.0_wp .and. vicen == 0.6_wp*1.5_wp, &
                 "snow alone vanishes: snow enthalpy untouched")
      if (allocated(error)) return
      ! ice vanishes
      call ice_update_state_vthermo(4, 1, -1.8_wp, 0.0_wp, 0.3_wp, Tsf, zqin, zqsn, aicen, vicen, vsnon, ierr)
      call check(error, aicen == 0.0_wp .and. vicen == 0.0_wp .and. vsnon == 0.0_wp .and. Tsf == -1.8_wp &
                 .and. all(zqin == 0.0_wp) .and. zqsn(1) == 0.0_wp .and. ierr == 0_int32, &
                 "ice vanishes: everything zero, Tsf = Tf")
   end subroutine test_update

   logical function same(a, b)
      !! Bitwise-equal, treating a NaN pair as equal.
      real(wp), intent(in) :: a, b
      same = (ieee_is_nan(a) .and. ieee_is_nan(b)) .or. (a == b)
   end function same

   subroutine test_nan(error)
      !! A NaN in any input of each routine sets `ICE_ERR_NONFINITE` (OR-ed
      !! into an existing word) and leaves every argument bitwise untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: nan
      real(wp) :: hin, hsn, dzi(4), dzs(1), zqin(4), zqsn(1), snoice, dsnow, hilyr, hslyr
      real(wp) :: hin0, hsn0, dzi0(4), dzs0(1), zqin0(4), zqsn0(1), snoice0, dsnow0
      real(wp) :: Tsf, aicen, vicen, vsnon, Tf, Tsf0, a0, vi0, vs0
      integer(int32) :: ierr
      integer :: which, nbad

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      nbad = 0
      do which = 1, 7
         hin = 0.5_wp; hsn = 0.5_wp
         dzi = 0.125_wp; dzs = 0.5_wp
         zqin = [-2.7e8_wp, -1.8e8_wp, -2.6e8_wp, -1.9e8_wp]
         zqsn = -1.1e8_wp
         snoice = 0.01_wp; dsnow = -0.01_wp
         select case (which)
         case (1); hin = nan
         case (2); hsn = nan
         case (3); zqin(1) = nan
         case (4); zqsn(1) = nan
         case (5); dzi(1) = nan
         case (6); dzs(1) = nan
         case (7); snoice = nan
         end select
         hin0 = hin; hsn0 = hsn; dzi0 = dzi; dzs0 = dzs; zqin0 = zqin; zqsn0 = zqsn; snoice0 = snoice; dsnow0 = dsnow
         ierr = 4_int32
         call ice_freeboard(4, 1, c, hin, hsn, zqin, zqsn, dzi, dzs, snoice, dsnow, ierr)
         if (iand(ierr, ICE_ERR_NONFINITE) == 0_int32 .or. iand(ierr, 4_int32) == 0_int32) nbad = nbad + 1
         if (.not. (same(hin, hin0) .and. same(hsn, hsn0) .and. same(snoice, snoice0) .and. same(dsnow, dsnow0) &
                    .and. same(dzi(1), dzi0(1)) .and. same(dzi(2), dzi0(2)) .and. same(dzs(1), dzs0(1)) &
                    .and. same(zqin(1), zqin0(1)) .and. same(zqin(2), zqin0(2)) .and. same(zqin(4), zqin0(4)) &
                    .and. same(zqsn(1), zqsn0(1)))) nbad = nbad + 100
      end do
      call check(error, nbad == 0, "freeboard: NaN flagged, untouched")
      if (allocated(error)) return

      nbad = 0
      do which = 1, 4
         dzi = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
         hin = 1.0_wp; hsn = 0.2_wp
         dzs = 0.2_wp
         zqin = [-2.7e8_wp, -1.8e8_wp, -2.6e8_wp, -1.9e8_wp]
         zqsn = -1.1e8_wp
         select case (which)
         case (1); hin = nan
         case (2); dzi(3) = nan
         case (3); zqin(4) = nan
         case (4); zqsn(1) = nan
         end select
         hin0 = hin; hsn0 = hsn; dzi0 = dzi; dzs0 = dzs; zqin0 = zqin; zqsn0 = zqsn
         ierr = 0_int32
         call ice_relayer(4, 1, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)
         if (iand(ierr, ICE_ERR_NONFINITE) == 0_int32) nbad = nbad + 1
         if (.not. (same(hin, hin0) .and. same(hsn, hsn0) .and. same(dzi(3), dzi0(3)) .and. same(dzi(1), dzi0(1)) &
                    .and. same(zqin(4), zqin0(4)) .and. same(zqin(1), zqin0(1)) .and. same(zqsn(1), zqsn0(1)))) nbad = nbad + 100
      end do
      call check(error, nbad == 0, "relayer: NaN flagged, untouched")
      if (allocated(error)) return

      nbad = 0
      do which = 1, 6
         zqin = [-2.7e8_wp, -1.8e8_wp, -2.6e8_wp, -1.9e8_wp]
         zqsn = -1.1e8_wp
         Tsf = -5.0_wp; aicen = 0.5_wp; vicen = 0.7_wp; vsnon = 0.2_wp; Tf = -1.8_wp; hin = 1.0_wp; hsn = 0.2_wp
         select case (which)
         case (1); hin = nan
         case (2); hsn = nan
         case (3); aicen = nan
         case (4); Tf = nan
         case (5); zqin(2) = nan
         case (6); Tsf = nan
         end select
         Tsf0 = Tsf; a0 = aicen; vi0 = vicen; vs0 = vsnon; zqin0 = zqin; zqsn0 = zqsn
         ierr = 0_int32
         call ice_update_state_vthermo(4, 1, Tf, hin, hsn, Tsf, zqin, zqsn, aicen, vicen, vsnon, ierr)
         if (iand(ierr, ICE_ERR_NONFINITE) == 0_int32) nbad = nbad + 1
         if (.not. (same(Tsf, Tsf0) .and. same(aicen, a0) .and. same(vicen, vi0) .and. same(vsnon, vs0) &
                    .and. same(zqin(2), zqin0(2)) .and. same(zqsn(1), zqsn0(1)))) nbad = nbad + 100
      end do
      call check(error, nbad == 0, "update_state: NaN flagged, untouched")
   end subroutine test_nan

   subroutine test_device(error)
      !! The three routines called per (i, j, n) from one
      !! `do concurrent (n, j, i)` over mapped `(nx, ny, ncat[, nlyr])` arrays of
      !! 256 x 256 x 5 x 7 (every array and the constants bundle
      !! device-resident under `mem:separate`), compared with the same kernel
      !! run serially on the host (tolerance, not bitwise: FMA is on in the GPU
      !! build; 1e-12 relative plus a 1e-14 m absolute floor, 1e-6 J/m3 on the
      !! enthalpies).  The timings are informational (one line on stdout).
      !!
      !! The loop indices reach the kernel as `value` dummies: passing the
      !! `do concurrent` indices by reference to the `routine seq` kernel made
      !! nvfortran 26.5 read wrong cells (silently, no conversion in 99.9 % of
      !! the cells).
      !!
      !! Layout note: the state's layer index is LAST, so the kernel gathers
      !! the two layered columns (`dzi`, `zqin`) into per-thread locals once and
      !! scatters them back; the per-call local footprint is reported in the
      !! header of `cell_kernel`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      integer, parameter :: NX = 256, NY = 256, NC = 5, NL = 7
      real(wp), allocatable :: hin(:, :, :), hsn(:, :, :), snoice(:, :, :), dsnow(:, :, :), aicen(:, :, :)
      real(wp), allocatable :: Tsf(:, :, :), vicen(:, :, :), vsnon(:, :, :), dzi(:, :, :, :), zqin(:, :, :, :)
      real(wp), allocatable :: zqsn(:, :, :)
      real(wp), allocatable :: hin_h(:, :, :), hsn_h(:, :, :), snoice_h(:, :, :), dsnow_h(:, :, :), aicen_h(:, :, :)
      real(wp), allocatable :: Tsf_h(:, :, :), vicen_h(:, :, :), vsnon_h(:, :, :), dzi_h(:, :, :, :), zqin_h(:, :, :, :)
      real(wp), allocatable :: zqsn_h(:, :, :)
      integer(int32), allocatable :: ierr(:, :, :), ierr_h(:, :, :)
      integer(int64) :: st, t0, t1, t2, t3, rate
      integer :: i, j, n, nconv, nmelt
      logical :: bad

      allocate (hin(NX, NY, NC), hsn(NX, NY, NC), snoice(NX, NY, NC), dsnow(NX, NY, NC), aicen(NX, NY, NC))
      allocate (Tsf(NX, NY, NC), vicen(NX, NY, NC), vsnon(NX, NY, NC), dzi(NX, NY, NC, NL), zqin(NX, NY, NC, NL))
      allocate (zqsn(NX, NY, NC), ierr(NX, NY, NC))
      allocate (hin_h, hsn_h, snoice_h, dsnow_h, aicen_h, Tsf_h, vicen_h, vsnon_h, zqsn_h, source=hin)
      allocate (dzi_h, zqin_h, source=dzi)
      allocate (ierr_h, source=ierr)

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      st = 4711_int64
      nmelt = 0
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               call random_column(st, NL, hin(i, j, n), hsn(i, j, n), dzi(i, j, n, :), zqin(i, j, n, :), zqsn(i, j, n), c)
               aicen(i, j, n) = 0.05_wp + 0.95_wp*urand(st)
               Tsf(i, j, n) = -10.0_wp*urand(st)
               if (mod(i + 3*j + n, 97) == 0) then   ! melt-out columns
                  hin(i, j, n) = 0.0_wp
                  dzi(i, j, n, :) = 0.0_wp
                  nmelt = nmelt + 1
               end if
            end do
         end do
      end do
      snoice = 0.0_wp; dsnow = 0.0_wp; vicen = 0.0_wp; vsnon = 0.0_wp; ierr = 0_int32
      hin_h = hin; hsn_h = hsn; snoice_h = snoice; dsnow_h = dsnow; aicen_h = aicen; Tsf_h = Tsf
      vicen_h = vicen; vsnon_h = vsnon; zqsn_h = zqsn; dzi_h = dzi; zqin_h = zqin; ierr_h = ierr

      call system_clock(t0, rate)
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               call cell_kernel(c, NX, NY, NC, NL, i, j, n, hin_h, hsn_h, dzi_h, zqin_h, zqsn_h, snoice_h, dsnow_h, &
                                Tsf_h, aicen_h, vicen_h, vsnon_h, ierr_h)
            end do
         end do
      end do
      call system_clock(t1)

      !$acc enter data copyin(c, hin, hsn, dzi, zqin, zqsn, snoice, dsnow, Tsf, aicen, vicen, vsnon, ierr)
      call system_clock(t2)
      do concurrent(n=1:NC, j=1:NY, i=1:NX)
         call cell_kernel(c, NX, NY, NC, NL, i, j, n, hin, hsn, dzi, zqin, zqsn, snoice, dsnow, Tsf, aicen, &
                          vicen, vsnon, ierr)
      end do
      call system_clock(t3)
      !$acc update self(hin, hsn, dzi, zqin, zqsn, snoice, dsnow, Tsf, aicen, vicen, vsnon, ierr)
      !$acc exit data delete(c, hin, hsn, dzi, zqin, zqsn, snoice, dsnow, Tsf, aicen, vicen, vsnon, ierr)
      write (*, '(a,es10.3,a,es10.3,a)') "  thickness kernel throughput [cell-categories/s]: host serial ", &
         real(NX*NY*NC, wp)/(real(t1 - t0, wp)/real(rate, wp)), "  do concurrent ", &
         real(NX*NY*NC, wp)/(real(max(t3 - t2, 1_int64), wp)/real(rate, wp)), ""

      call check(error, all(ierr == 0_int32) .and. all(ierr_h == 0_int32), "no error bits")
      if (allocated(error)) return
      bad = .false.
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               if (.not. (ieee_is_finite(hin(i, j, n)) .and. ieee_is_finite(zqin(i, j, n, 1)) &
                          .and. ieee_is_finite(zqin(i, j, n, NL)) .and. ieee_is_finite(vicen(i, j, n)) &
                          .and. ieee_is_finite(vsnon(i, j, n)) .and. ieee_is_finite(snoice(i, j, n)))) bad = .true.
            end do
         end do
      end do
      call check(error,.not. bad, "device results finite")
      if (allocated(error)) return
      call check(error, all(abs(hin - hin_h) <= 1.0e-12_wp*abs(hin_h) + 1.0e-14_wp) .and. &
                 all(abs(hsn - hsn_h) <= 1.0e-12_wp*abs(hsn_h) + 1.0e-14_wp), "thicknesses device == host")
      if (allocated(error)) return
      call check(error, all(abs(zqin - zqin_h) <= 1.0e-12_wp*abs(zqin_h) + 1.0e-6_wp), "zqin device == host")
      if (allocated(error)) return
      call check(error, all(abs(vicen - vicen_h) <= 1.0e-12_wp*abs(vicen_h) + 1.0e-14_wp) .and. &
                 all(abs(vsnon - vsnon_h) <= 1.0e-12_wp*abs(vsnon_h) + 1.0e-14_wp), "volumes device == host")
      if (allocated(error)) return
      call check(error, all(abs(snoice - snoice_h) <= 1.0e-12_wp*abs(snoice_h) + 1.0e-14_wp), "snoice device == host")
      if (allocated(error)) return
      nconv = count(snoice_h > 0.0_wp)
      call check(error, nconv > 10000 .and. nmelt > 1000, "sample has conversions and melt-out columns")
   end subroutine test_device

   pure subroutine cell_kernel(c, nx, ny, nc, nl, i, j, n, hin, hsn, dzi, zqin, zqsn, snoice, dsnow, Tsf, aicen, &
                               vicen, vsnon, ierr)
      !! Freeboard, relayer and state update of the column (i, j, n) of
      !! `(nx, ny, ncat[, nlyr])` state arrays, in place.  The layer index is
      !! the last one, so `dzi` and `zqin` are gathered once into per-thread
      !! locals and scattered back: 2 x 7 x 8 B = 112 B of columns, plus the
      !! callee locals (relayer: zi1, zi2 9 x 8 B each, one more 7 x 8 B copy of
      !! the enthalpy; `adjust_enthalpy`: 7 x 8 B), about 350 B per call.
      !$acc routine seq
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nx, ny, nc, nl
      integer, value :: i, j, n
      real(wp), intent(inout) :: hin(nx, ny, nc), hsn(nx, ny, nc), dzi(nx, ny, nc, nl), zqin(nx, ny, nc, nl)
      real(wp), intent(inout) :: zqsn(nx, ny, nc), snoice(nx, ny, nc), dsnow(nx, ny, nc), Tsf(nx, ny, nc)
      real(wp), intent(inout) :: aicen(nx, ny, nc), vicen(nx, ny, nc), vsnon(nx, ny, nc)
      integer(int32), intent(inout) :: ierr(nx, ny, nc)
      integer :: k
      real(wp) :: dz(ICE_NILYR_MAX), zq(ICE_NILYR_MAX), dzs(1), zqs(1), hil, hsl, hi, hs, ts, ai, vi, vs, sno, dsn
      integer(int32) :: ie

      ie = ierr(i, j, n)
      do k = 1, nl
         dz(k) = dzi(i, j, n, k)
         zq(k) = zqin(i, j, n, k)
      end do
      hi = hin(i, j, n)
      hs = hsn(i, j, n)
      dzs(1) = hs
      zqs(1) = zqsn(i, j, n)
      ts = Tsf(i, j, n)
      ai = aicen(i, j, n)
      sno = snoice(i, j, n)
      dsn = dsnow(i, j, n)
      call ice_freeboard(nl, 1, c, hi, hs, zq, zqs, dz, dzs, sno, dsn, ie)
      call ice_relayer(nl, 1, c, dz, dzs, hi, hs, zq, zqs, hil, hsl, ie)
      vi = 0.0_wp
      vs = 0.0_wp
      call ice_update_state_vthermo(nl, 1, -1.8_wp, hi, hs, ts, zq, zqs, ai, vi, vs, ie)
      do k = 1, nl
         dzi(i, j, n, k) = dz(k)
         zqin(i, j, n, k) = zq(k)
      end do
      snoice(i, j, n) = sno
      dsnow(i, j, n) = dsn
      hin(i, j, n) = hi
      hsn(i, j, n) = hs
      zqsn(i, j, n) = zqs(1)
      Tsf(i, j, n) = ts
      aicen(i, j, n) = ai
      vicen(i, j, n) = vi
      vsnon(i, j, n) = vs
      ierr(i, j, n) = ie
   end subroutine cell_kernel

end module test_cice_thickness_pieces
