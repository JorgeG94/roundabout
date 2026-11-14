module test_cice_atm_inputs
   !! Tests of `rdb_ice_atm_inputs` (D8.2): the closed forms of `potT`, `rhoa`, `wind` for two
   !! atmospheres (hand numbers); the kelvin / degC handling; the rotation covariance of the
   !! wind stress through the REAL boundary-layer kernel (wind rotated: stress rotates, scalar
   !! fluxes unchanged; wind and grid angle rotated together: nothing changes); the
   !! high-frequency T velocity (uniform face velocity, masked faces excluded, no ice gives
   !! zero) and its effect on the boundary layer (air moving with the ice behaves as calm
   !! air with the 0.5 m/s floor, `highfreq` off ignores the ice velocity bitwise); the
   !! T-to-face stress map (hand numbers, uniform field); land and non-finite behaviour; the
   !! builder feeding `ice_step_therm1` against hand-built inputs; a device sweep against the
   !! host on every output of every cell (the ring untouched) and the throughput.
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, &
                             ice_q_from_T, ice_qsno_from_T, ice_salinity_profile
   use rdb_ice_atmo, only: ice_atm_boundary_ice
   use rdb_ice_therm1, only: ice_step_therm1
   use rdb_ice_atm_inputs, only: ice_atm_prepare_cell_inputs, ice_atm_cell_inputs, ice_atm_heights, &
                                 ice_atm_velocity_T, ice_atm_stress_to_faces, ice_atm_grid_rotation, ICE_ATM_BENIGN_POTT, &
                                 ICE_ATM_BENIGN_RHOA
   implicit none
   private

   public :: collect_cice_atm_inputs_tests

#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: DRTOL = 1.0e-12_wp
      !! Device-vs-host (and hand-built-vs-builder) relative tolerance, GPU build (FMA on, host code too).
#else
   real(wp), parameter :: DRTOL = 0.0_wp
#endif
   real(wp), parameter :: SENT = -777.0_wp

contains

   subroutine collect_cice_atm_inputs_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("potT_rhoa_wind_closed_forms", test_closed_forms), &
                  new_unittest("kelvin_and_degC_sources", test_units), &
                  new_unittest("rotation_covariance_of_stress", test_rotation), &
                  new_unittest("highfreq_velocity_T_average", test_velocity), &
                  new_unittest("highfreq_boundary_layer_effect", test_highfreq_bl), &
                  new_unittest("stress_to_faces_map", test_faces), &
                  new_unittest("land_bad_cells_and_heights", test_land_bad), &
                  new_unittest("builder_into_therm1_matches_hand_inputs", test_therm1), &
                  new_unittest("device_matches_host_every_cell", test_device), &
                  new_unittest("throughput_512x512", test_throughput) &
                  ]
   end subroutine collect_cice_atm_inputs_tests

   logical function near(x, y, rtol, atol)
      !! NaN-proof closeness: a non-finite `x` is never near.
      real(wp), intent(in) :: x, y, rtol, atol
      near = .false.
      if (.not. ieee_is_finite(x)) return
      near = abs(x - y) <= atol + rtol*abs(y)
   end function near

   function uni(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function uni

   subroutine d_rotation(nx, ny, is, ie, js, je, ang, cosA, sinA)
      !! `ice_atm_grid_rotation` with every array mapped for the call.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: ang(nx, ny)
      real(wp), intent(inout) :: cosA(nx, ny), sinA(nx, ny)
      !$acc enter data copyin(ang, cosA, sinA)
      call ice_atm_grid_rotation(nx, ny, is, ie, js, je, ang, cosA, sinA)
      !$acc update self(cosA, sinA)
      !$acc exit data delete(ang, cosA, sinA)
   end subroutine d_rotation

   subroutine d_prepare(nx, ny, is, ie, js, je, off, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, &
                        potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ice_err)
      !! The setup rotation followed by `ice_atm_prepare_cell_inputs`, every array mapped for the calls (GPU
      !! build: `mem:separate`).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: off
      real(wp), dimension(nx, ny), intent(in) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn
      real(wp), dimension(nx, ny), intent(inout) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      integer(int32), intent(inout) :: ice_err(nx, ny)
      real(wp) :: cosA(nx, ny), sinA(nx, ny)
      cosA = 1.0_wp
      sinA = 0.0_wp
      call d_rotation(nx, ny, is, ie, js, je, ang, cosA, sinA)
      !$acc enter data copyin(hm, cosA, sinA, tas, huss, uas, vas, psl, rlds, prra, prsn, potT, rhoa, wind, uatm, vatm, &
      !$acc&  Qa, flw, fsn, frn, ice_err)
      call ice_atm_prepare_cell_inputs(nx, ny, is, ie, js, je, off, hm, cosA, sinA, tas, huss, uas, vas, psl, rlds, prra, &
                                       prsn, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ice_err)
      !$acc update self(potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ice_err)
      !$acc exit data delete(hm, cosA, sinA, tas, huss, uas, vas, psl, rlds, prra, prsn, potT, rhoa, wind, uatm, vatm, &
      !$acc&  Qa, flw, fsn, frn, ice_err)
   end subroutine d_prepare

   subroutine d_velocity(nx, ny, is, ie, js, je, amin, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ice_err)
      !! `ice_atm_velocity_T` with every array mapped for the call.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: amin
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), dimension(nx, ny), intent(in) :: ea, na, ua, uvm, hm, ai
      real(wp), dimension(nx, ny), intent(inout) :: uT, vT
      integer(int32), intent(inout) :: ice_err(nx, ny)
      !$acc enter data copyin(u, v, ea, na, ua, uvm, hm, ai, uT, vT, ice_err)
      call ice_atm_velocity_T(nx, ny, is, ie, js, je, amin, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ice_err)
      !$acc update self(uT, vT, ice_err)
      !$acc exit data delete(u, v, ea, na, ua, uvm, hm, ai, uT, vT, ice_err)
   end subroutine d_velocity

   subroutine d_faces(nx, ny, is, ie, js, je, sxT, syT, ta, ea, na, sxE, syN)
      !! `ice_atm_stress_to_faces` with every array mapped for the call.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), dimension(nx, ny), intent(in) :: sxT, syT, ta, ea, na
      real(wp), dimension(nx, ny), intent(inout) :: sxE, syN
      !$acc enter data copyin(sxT, syT, ta, ea, na, sxE, syN)
      call ice_atm_stress_to_faces(nx, ny, is, ie, js, je, sxT, syT, ta, ea, na, sxE, syN)
      !$acc update self(sxE, syN)
      !$acc exit data delete(sxT, syT, ta, ea, na, sxE, syN)
   end subroutine d_faces

   subroutine cin(tas, huss, uas, vas, psl, rlds, prra, prsn, angle, off, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
      !! The scalar kernel with the grid angle given as an angle (libm `cos`, `sin`, independent of the kernel's).
      real(wp), intent(in) :: tas, huss, uas, vas, psl, rlds, prra, prsn, angle, off
      real(wp), intent(out) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      logical, intent(out) :: ok
      call ice_atm_cell_inputs(tas, huss, uas, vas, psl, rlds, prra, prsn, cos(angle), sin(angle), off, &
                               potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
   end subroutine cin

   subroutine one(tas, huss, uas, vas, psl, angle, off, potT, rhoa, wind, uatm, vatm, ok)
      !! The scalar kernel with the other fields at neutral values.
      real(wp), intent(in) :: tas, huss, uas, vas, psl, angle, off
      real(wp), intent(out) :: potT, rhoa, wind, uatm, vatm
      logical, intent(out) :: ok
      real(wp) :: Qa, flw, fsn, frn
      call cin(tas, huss, uas, vas, psl, 250.0_wp, 1.0e-5_wp, 2.0e-5_wp, angle, off, &
               potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
   end subroutine one

   subroutine test_closed_forms(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      logical :: ok

      ! Atmosphere A: psl = 1000 hPa, so potT = tas exactly (a unit ratio); rhoa from the ideal-gas law
      call cin(270.0_wp, 0.002_wp, 3.0_wp, -4.0_wp, 1.0e5_wp, 230.0_wp, 1.0e-5_wp, 2.0e-5_wp, &
               0.0_wp, 0.0_wp, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
      call check(error, ok, "A ok")
      if (allocated(error)) return
      call check(error, potT == 270.0_wp, "A potT = tas at 1000 hPa")
      if (allocated(error)) return
      call check(error, near(rhoa, 1.2886613401655285_wp, 1.0e-14_wp, 0.0_wp), "A rhoa hand number")
      if (allocated(error)) return
      call check(error, wind == 5.0_wp .and. uatm == 3.0_wp .and. vatm == -4.0_wp, "A wind 3-4-5 at zero angle")
      if (allocated(error)) return
      call check(error, Qa == 0.002_wp .and. flw == 230.0_wp .and. fsn == 2.0e-5_wp .and. frn == 1.0e-5_wp, &
                 "A passthrough Qa, flw, fsnow, frain")
      if (allocated(error)) return
      ! Atmosphere B: 980 hPa, 250 K, 0.5 g/kg
      call cin(250.0_wp, 0.0005_wp, 0.0_wp, 0.0_wp, 98000.0_wp, 200.0_wp, 0.0_wp, 0.0_wp, &
               0.0_wp, 0.0_wp, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
      call check(error, ok .and. near(potT, 251.44867474559197_wp, 1.0e-14_wp, 0.0_wp), "B potT hand number")
      if (allocated(error)) return
      call check(error, near(rhoa, 1.3651626786783937_wp, 1.0e-14_wp, 0.0_wp), "B rhoa hand number")
      if (allocated(error)) return
      call check(error, wind == 0.0_wp, "B calm")
      if (allocated(error)) return
      ! ABSOLUTE rotation sign: angle = pi/2, uas = 1, vas = 0 gives uatm = 0, vatm = -1 (u = cosA uas + sinA vas,
      ! v = -sinA uas + cosA vas); the covariance test cannot see a flipped sign
      call cin(270.0_wp, 0.002_wp, 1.0_wp, 0.0_wp, 1.0e5_wp, 230.0_wp, 0.0_wp, 0.0_wp, 1.5707963267948966_wp, 0.0_wp, &
               potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ok)
      call check(error, ok .and. abs(uatm) < 1.0e-15_wp .and. near(vatm, -1.0_wp, 1.0e-15_wp, 0.0_wp), &
                 "angle pi/2: uas = 1 gives (uatm, vatm) = (0, -1)")
   end subroutine test_closed_forms

   subroutine test_units(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: p1, r1, w1, u1, v1, p2, r2, w2, u2, v2
      logical :: ok1, ok2

      call one(270.0_wp, 0.002_wp, 3.0_wp, 4.0_wp, 99000.0_wp, 0.3_wp, 0.0_wp, p1, r1, w1, u1, v1, ok1)
      call one(270.0_wp - 273.15_wp, 0.002_wp, 3.0_wp, 4.0_wp, 99000.0_wp, 0.3_wp, 273.15_wp, p2, r2, w2, u2, v2, ok2)
      call check(error, ok1 .and. ok2, "ok")
      if (allocated(error)) return
      call check(error, near(p2, p1, 1.0e-13_wp, 0.0_wp) .and. near(r2, r1, 1.0e-13_wp, 0.0_wp), &
                 "degC with the offset equals the kelvin source")
      if (allocated(error)) return
      call check(error, w2 == w1 .and. u2 == u1 .and. v2 == v1, "the wind does not see the offset")
      if (allocated(error)) return
      ! degC without the offset is a negative kelvin: refused, never silently used
      call one(-3.15_wp, 0.002_wp, 3.0_wp, 4.0_wp, 99000.0_wp, 0.3_wp, 0.0_wp, p2, r2, w2, u2, v2, ok2)
      call check(error,.not. ok2 .and. p2 == ICE_ATM_BENIGN_POTT .and. r2 == ICE_ATM_BENIGN_RHOA, "degC passed as K refused")
   end subroutine test_units

   subroutine bl(p, potT, rhoa, Qa, uatm, vatm, uvel, vvel, Tsf, sh, lh, sx, sy)
      !! The real boundary-layer kernel on one surface (heights 10 m).
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: potT, rhoa, Qa, uatm, vatm, uvel, vvel, Tsf
      real(wp), intent(out) :: sh, lh, sx, sy
      real(wp) :: dt, dq
      integer(int32) :: ierr
      ierr = 0_int32
      call ice_atm_boundary_ice(Tsf, potT, uatm, vatm, uvel, vvel, sqrt(uatm*uatm + vatm*vatm), 10.0_wp, 10.0_wp, Qa, &
                                rhoa, p%c, p%atm, sh, lh, sx, sy, dt, dq, ierr)
   end subroutine bl

   subroutine test_rotation(error)
      !! Driver test 6: rotating the earth-relative wind by `d` with the grid angle fixed rotates the model-axes
      !! wind and the stress by `d` and leaves the scalar fluxes alone; rotating wind and angle together changes
      !! nothing.  With `highfreq` the ice velocity rotates too.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer :: m, k
      real(wp), parameter :: ANG(4) = [0.0_wp, 0.4_wp, -1.1_wp, 2.5_wp]
      real(wp), parameter :: DEL(3) = [0.3_wp, 1.7_wp, -2.2_wp]
      real(wp) :: uas, vas, a, d, cd, sd, potT, rhoa, wind, ua0, va0, ua1, va1, ua2, va2, Qa, flw, fsn, frn
      real(wp) :: sh0, lh0, sx0, sy0, sh1, lh1, sx1, sy1, sh2, lh2, sx2, sy2, ui, vi, ui1, vi1, ex, ey
      logical :: ok

      do m = 0, 1
         p = ice_params_default(ICE_CONSTANTS_ICEPACK, 1, 4, 1)
         p%atm%highfreq = m == 1
         do k = 1, 12
            uas = 3.0_wp + 0.7_wp*real(k, wp)
            vas = -2.0_wp + 0.9_wp*real(k, wp)
            a = ANG(mod(k, 4) + 1)
            d = DEL(mod(k, 3) + 1)
            cd = cos(d)
            sd = sin(d)
            ui = 0.12_wp*real(k, wp)
            vi = -0.05_wp*real(k, wp)
            ui1 = cd*ui - sd*vi
            vi1 = sd*ui + cd*vi
            ! base
            call cin(262.0_wp, 6.0e-4_wp, uas, vas, 1.01e5_wp, 200.0_wp, 0.0_wp, 0.0_wp, a, 0.0_wp, &
                     potT, rhoa, wind, ua0, va0, Qa, flw, fsn, frn, ok)
            call bl(p, potT, rhoa, Qa, ua0, va0, ui, vi, -8.0_wp, sh0, lh0, sx0, sy0)
            ! wind rotated, angle fixed
            call cin(262.0_wp, 6.0e-4_wp, cd*uas - sd*vas, sd*uas + cd*vas, 1.01e5_wp, 200.0_wp, &
                     0.0_wp, 0.0_wp, a, 0.0_wp, potT, rhoa, wind, ua1, va1, Qa, flw, fsn, frn, ok)
            call bl(p, potT, rhoa, Qa, ua1, va1, ui1, vi1, -8.0_wp, sh1, lh1, sx1, sy1)
            ex = cd*sx0 - sd*sy0
            ey = sd*sx0 + cd*sy0
            call check(error, near(sx1, ex, 1.0e-12_wp, 1.0e-15_wp) .and. near(sy1, ey, 1.0e-12_wp, 1.0e-15_wp), &
                       "stress rotates by the wind rotation")
            if (allocated(error)) return
            call check(error, near(sh1, sh0, 1.0e-12_wp, 0.0_wp) .and. near(lh1, lh0, 1.0e-12_wp, 0.0_wp), &
                       "scalar fluxes unchanged")
            if (allocated(error)) return
            call check(error, near(hypot(ua1, va1), hypot(ua0, va0), 1.0e-13_wp, 0.0_wp), "model wind speed unchanged")
            if (allocated(error)) return
            ! wind and angle rotated together: the model-axes wind is the same
            call cin(262.0_wp, 6.0e-4_wp, cd*uas - sd*vas, sd*uas + cd*vas, 1.01e5_wp, 200.0_wp, &
                     0.0_wp, 0.0_wp, a + d, 0.0_wp, potT, rhoa, wind, ua2, va2, Qa, flw, fsn, frn, ok)
            call bl(p, potT, rhoa, Qa, ua2, va2, ui, vi, -8.0_wp, sh2, lh2, sx2, sy2)
            call check(error, near(ua2, ua0, 1.0e-12_wp, 1.0e-13_wp) .and. near(va2, va0, 1.0e-12_wp, 1.0e-13_wp), &
                       "grid-angle covariance of the model-axes wind")
            if (allocated(error)) return
            call check(error, near(sx2, sx0, 1.0e-11_wp, 1.0e-14_wp) .and. near(sy2, sy0, 1.0e-11_wp, 1.0e-14_wp) &
                       .and. near(sh2, sh0, 1.0e-11_wp, 0.0_wp), "stress and fluxes invariant under the joint rotation")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_rotation

   subroutine grid_world(nx, ny, hm, earea, narea, uarea, uvm, tarea)
      !! A uniform grid for the helpers: every cell wet, face areas 3 (east) and 2 (north), corner area 1, T area 4.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: hm(nx, ny), earea(nx, ny), narea(nx, ny), uarea(nx, ny), uvm(nx, ny), tarea(nx, ny)
      hm = 1.0_wp
      uvm = 1.0_wp
      earea = 3.0_wp
      narea = 2.0_wp
      uarea = 1.0_wp
      tarea = 4.0_wp
   end subroutine grid_world

   subroutine test_velocity(error)
      !! The reference C-grid T velocity: faces -> corners (masked two-face mean) -> T (corner-area weights, unmasked).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 7
      real(wp) :: u(NX + 1, NY), v(NX, NY + 1), uT(NX, NY), vT(NX, NY), hm(NX, NY), ea(NX, NY), na(NX, NY)
      real(wp) :: ua(NX, NY), uvm(NX, NY), ta(NX, NY), ai(NX, NY), nan
      integer(int32) :: ie(NX, NY)
      integer :: k, dj, di
      real(wp), parameter :: WX(-1:1) = [1.0_wp, 2.0_wp, 1.0_wp]

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      call grid_world(NX, NY, hm, ea, na, ua, uvm, ta)
      ai = 0.8_wp
      ie = 0_int32
      ! uniform field is its own average (exactly: power-of-two weights), ring untouched
      u = 0.3_wp
      v = -0.2_wp
      uT = SENT
      vT = SENT
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, all(uT(2:NX - 1, 2:NY - 1) == 0.3_wp) .and. all(vT(2:NX - 1, 2:NY - 1) == -0.2_wp), &
                 "uniform face velocity is its own T average")
      if (allocated(error)) return
      call check(error, all(uT(1, :) == SENT) .and. all(vT(:, NY) == SENT), "ring untouched")
      if (allocated(error)) return
      ! stencil weights by hand on the uniform grid: cell (4,4); the east column u(5, 3:5), the west column u(4, 3:5)
      ! carry (1/8)(1, 2, 1) each; nothing else enters
      do dj = -1, 1
         do k = 0, 1
            u = 0.0_wp
            v = 0.0_wp
            u(4 + k, 4 + dj) = 8.0_wp
            call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
            call check(error, near(uT(4, 4), WX(dj), 1.0e-15_wp, 0.0_wp) .and. vT(4, 4) == 0.0_wp, &
                       "east/west column face weight (1,2,1)/8")
            if (allocated(error)) return
         end do
      end do
      u = 0.0_wp
      u(6, 4) = 8.0_wp          ! the next column east of the stencil
      u(3, 4) = 8.0_wp
      u(5, 6) = 8.0_wp          ! two rows up
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, uT(4, 4) == 0.0_wp, "faces outside the six-face stencil do not enter")
      if (allocated(error)) return
      ! y: south row v(3:5, 4) and north row v(3:5, 5) of cell (4,4)
      do di = -1, 1
         do k = 0, 1
            u = 0.0_wp
            v = 0.0_wp
            v(4 + di, 4 + k) = 8.0_wp
            call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
            call check(error, near(vT(4, 4), WX(di), 1.0e-15_wp, 0.0_wp) .and. uT(4, 4) == 0.0_wp, &
                       "south/north row face weight (1,2,1)/8")
            if (allocated(error)) return
         end do
      end do
      ! coast: the two west corners of cell (4,4) masked (closed faces carry 0): the mean shrinks to one half
      u = 0.3_wp
      v = 0.0_wp
      uvm(3, 3) = 0.0_wp
      uvm(3, 4) = 0.0_wp
      u(4, 3:5) = 0.0_wp
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, near(uT(4, 4), 0.5_wp*0.3_wp, 1.0e-15_wp, 0.0_wp), "next to a coast the mean shrinks (masked corners keep weight)")
      if (allocated(error)) return
      ! corner areas weight the mean: the NE corner twice the others
      uvm = 1.0_wp
      ua = 1.0_wp
      ua(4, 4) = 2.0_wp
      u = 0.0_wp
      u(5, 3:5) = 0.6_wp         ! east column: NE corner U(4,4) (area 2) and SE corner U(4,3) (area 1)
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, near(uT(4, 4), 0.6_wp*3.0_wp/5.0_wp, 1.0e-15_wp, 0.0_wp), "corner-area weights")
      if (allocated(error)) return
      ! a non-finite face flags the cells that read it and returns zero there
      ua = 1.0_wp
      u = 0.3_wp
      u(5, 4) = nan
      ie = 0_int32
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, uT(4, 4) == 0.0_wp .and. ie(4, 4) == ICE_ERR_NONFINITE .and. all(ieee_is_finite(uT)), &
                 "NaN face: flagged, zero, no NaN out")
      if (allocated(error)) return
      ! gate addition: no ice and land give exactly zero
      u = 0.3_wp
      ie = 0_int32
      ai(4, 4) = 0.0_wp
      hm(5, 4) = 0.0_wp
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      call check(error, uT(4, 4) == 0.0_wp .and. vT(4, 4) == 0.0_wp .and. uT(5, 4) == 0.0_wp .and. vT(5, 4) == 0.0_wp &
                 .and. all(ie == 0_int32), "no ice and land cells are exactly zero")
   end subroutine test_velocity

   subroutine test_highfreq_bl(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: NX = 5, NY = 5
      real(wp) :: u(NX + 1, NY), v(NX, NY + 1), uT(NX, NY), vT(NX, NY), hm(NX, NY), ea(NX, NY), na(NX, NY)
      real(wp) :: ua(NX, NY), uvm(NX, NY), ta(NX, NY), ai(NX, NY)
      integer(int32) :: ie(NX, NY)
      real(wp) :: sh, lh, sx, sy, sh0, lh0, sx0, sy0, sh1, lh1, sx1, sy1, potT, rhoa, Qa
      real(wp) :: dt, dq
      integer(int32) :: ierr

      call grid_world(NX, NY, hm, ea, na, ua, uvm, ta)
      ai = 0.9_wp
      u = 0.45_wp
      v = -0.35_wp
      ie = 0_int32
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      potT = 258.0_wp
      rhoa = 1.32_wp
      Qa = 5.0e-4_wp
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 1, 4, 1)
      ! highfreq ON: air moving exactly with the ice -> zero relative wind: zero stress, coefficients of calm air
      p%atm%highfreq = .true.
      call bl(p, potT, rhoa, Qa, uT(3, 3), vT(3, 3), uT(3, 3), vT(3, 3), -9.0_wp, sh, lh, sx, sy)
      call bl(p, potT, rhoa, Qa, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -9.0_wp, sh0, lh0, sx0, sy0)
      call check(error, sx == 0.0_wp .and. sy == 0.0_wp, "no relative wind, no stress")
      if (allocated(error)) return
      call check(error, sh == sh0 .and. lh == lh0 .and. sh > 1.0_wp, "coefficients of the 0.5 m/s minimum speed (Galilean)")
      if (allocated(error)) return
      ! highfreq OFF: the ice velocity is ignored bitwise
      p%atm%highfreq = .false.
      call bl(p, potT, rhoa, Qa, 6.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, -9.0_wp, sh0, lh0, sx0, sy0)
      ierr = 0_int32
      call ice_atm_boundary_ice(-9.0_wp, potT, 6.0_wp, -2.0_wp, uT(3, 3), vT(3, 3), sqrt(40.0_wp), 10.0_wp, 10.0_wp, Qa, &
                                rhoa, p%c, p%atm, sh1, lh1, sx1, sy1, dt, dq, ierr)
      call check(error, sh1 == sh0 .and. lh1 == lh0 .and. sx1 == sx0 .and. sy1 == sy0, "highfreq off ignores the ice velocity")
      if (allocated(error)) return
      ! highfreq ON with a moving ice sheet: weaker stress than the ice at rest
      p%atm%highfreq = .true.
      call bl(p, potT, rhoa, Qa, 6.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, -9.0_wp, sh0, lh0, sx0, sy0)
      call bl(p, potT, rhoa, Qa, 6.0_wp, -2.0_wp, uT(3, 3), vT(3, 3), -9.0_wp, sh1, lh1, sx1, sy1)
      call check(error, abs(sx1) < abs(sx0) .and. sx1 > 0.0_wp, "ice moving with the wind reduces the stress")
   end subroutine test_highfreq_bl

   subroutine test_faces(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 4
      real(wp) :: sxT(NX, NY), syT(NX, NY), ta(NX, NY), ea(NX, NY), na(NX, NY), sxE(NX, NY), syN(NX, NY)
      integer :: i, j

      ta = 4.0_wp
      ea = 4.0_wp
      na = 4.0_wp
      sxT = 0.0_wp
      syT = 0.0_wp
      sxE = SENT
      syN = SENT
      ! hand numbers: T areas 2 and 6, face area 3, stress 0.2 and 0.4 -> 0.5 (0.4 + 2.4)/3
      ta(2, 2) = 2.0_wp
      ta(3, 2) = 6.0_wp
      ea(2, 2) = 3.0_wp
      sxT(2, 2) = 0.2_wp
      sxT(3, 2) = 0.4_wp
      ta(2, 3) = 2.0_wp
      na(2, 2) = 3.0_wp
      syT(2, 2) = -0.3_wp
      syT(2, 3) = 0.1_wp
      ta(2, 3) = 6.0_wp
      call d_faces(NX, NY, 2, 3, 2, 2, sxT, syT, ta, ea, na, sxE, syN)
      call check(error, near(sxE(2, 2), 0.5_wp*(0.2_wp*2.0_wp + 0.4_wp*6.0_wp)/3.0_wp, 4.0e-16_wp, 0.0_wp), "E hand number")
      if (allocated(error)) return
      call check(error, near(syN(2, 2), 0.5_wp*(-0.3_wp*2.0_wp + 0.1_wp*6.0_wp)/3.0_wp, 1.0e-15_wp, 1.0e-17_wp), &
                 "N hand number")
      if (allocated(error)) return
      call check(error, sxE(1, 1) == SENT .and. syN(4, 4) == SENT, "outside the owned range untouched")
      if (allocated(error)) return
      ! uniform field on a uniform grid is preserved (face area = T area)
      sxT = 0.07_wp
      syT = -0.03_wp
      ta = 4.0_wp
      ea = 4.0_wp
      na = 4.0_wp
      call d_faces(NX, NY, 1, NX - 1, 1, NY - 1, sxT, syT, ta, ea, na, sxE, syN)
      do j = 1, NY - 1
         do i = 1, NX - 1
            if (.not. (near(sxE(i, j), 0.07_wp, 4.0e-16_wp, 0.0_wp) .and. near(syN(i, j), -0.03_wp, 4.0e-16_wp, 0.0_wp))) then
               call test_failed(error, "uniform stress not preserved")
               return
            end if
         end do
      end do
      ! zero face area gives zero (never a division)
      ea(2, 2) = 0.0_wp
      call d_faces(NX, NY, 1, NX - 1, 1, NY - 1, sxT, syT, ta, ea, na, sxE, syN)
      call check(error, sxE(2, 2) == 0.0_wp, "zero face area")
      if (allocated(error)) return
      call faces_identity(error)
   end subroutine test_faces

   subroutine faces_identity(error)
      !! The force identity `sum_E A_E S_E = sum_i A_i S_i` (every face area positive, so no face is dropped):
      !! on a periodic row (faces 1..N, T N+1 is the wrap of T 1) and on a closed row (land cells with zero stress
      !! on both ends, faces 1..N+1).  And the land-neighbour face gets half the cell's stress (area ratio 1),
      !! with no renormalisation.
      type(error_type), allocatable, intent(inout) :: error
      integer, parameter :: N = 9
      real(wp) :: sxT(N + 2, 1), syT(N + 2, 1), ta(N + 2, 1), ea(N + 2, 1), na(N + 2, 1), sE(N + 2, 1), sN(N + 2, 1)
      real(wp) :: lhs, rhs
      integer :: i
      integer(int64) :: st

      st = 31_int64
      syT = 0.0_wp
      na = 1.0_wp
      do i = 1, N + 2
         ta(i, 1) = 2.0_wp + uni(st)
         ea(i, 1) = 1.0_wp + 3.0_wp*uni(st)
         sxT(i, 1) = uni(st) - 0.5_wp
      end do
      ! periodic: cells 1..N, ghost N+1 = copy of cell 1
      sxT(N + 1, 1) = sxT(1, 1)
      ta(N + 1, 1) = ta(1, 1)
      sE = SENT
      call d_faces(N + 2, 1, 1, N, 1, 1, sxT, syT, ta, ea, na, sE, sN)
      lhs = sum(ea(1:N, 1)*sE(1:N, 1))
      rhs = sum(ta(1:N, 1)*sxT(1:N, 1))
      call check(error, near(lhs, rhs, 1.0e-13_wp, 1.0e-14_wp), "force identity on a periodic row")
      if (allocated(error)) return
      ! closed: land cells 1 and N+2 (zero stress), ice 2..N+1, faces 1..N+1
      sxT(1, 1) = 0.0_wp
      sxT(N + 2, 1) = 0.0_wp
      sE = SENT
      call d_faces(N + 2, 1, 1, N + 1, 1, 1, sxT, syT, ta, ea, na, sE, sN)
      lhs = sum(ea(1:N + 1, 1)*sE(1:N + 1, 1))
      rhs = sum(ta(2:N + 1, 1)*sxT(2:N + 1, 1))
      call check(error, near(lhs, rhs, 1.0e-13_wp, 1.0e-14_wp), "force identity on a closed row")
      if (allocated(error)) return
      ea(1, 1) = ta(2, 1)
      call d_faces(N + 2, 1, 1, N + 1, 1, 1, sxT, syT, ta, ea, na, sE, sN)
      call check(error, near(sE(1, 1), 0.5_wp*sxT(2, 1), 1.0e-15_wp, 0.0_wp), "face next to land: half the cell stress")
   end subroutine faces_identity

   subroutine test_land_bad(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 3
      real(wp), dimension(NX, NY) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn
      real(wp), dimension(NX, NY) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      integer(int32) :: ie(NX, NY)
      real(wp) :: nan, zl, zs
      logical :: ok

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      hm = 1.0_wp
      ang = 0.0_wp
      tas = 265.0_wp
      huss = 1.0e-3_wp
      uas = 4.0_wp
      vas = 3.0_wp
      psl = 1.0e5_wp
      rlds = 210.0_wp
      prra = 0.0_wp
      prsn = 1.0e-5_wp
      hm(1, 1) = 0.0_wp
      tas(1, 1) = nan            ! land with NaN: harmless
      psl(1, 1) = nan
      huss(2, 1) = nan           ! wet NaN
      tas(3, 1) = 0.0_wp         ! wet, non-positive temperature
      psl(4, 1) = 0.0_wp         ! wet, non-positive pressure
      rlds(1, 2) = nan           ! wet NaN in a pass-through field
      ang(2, 2) = nan            ! wet, non-finite grid angle
      ie = 0_int32
      ie(4, 3) = 2_int32         ! a prior bit stays (OR-ed, never cleared)
      potT = SENT
      call d_prepare(NX, NY, 1, NX, 1, NY, 0.0_wp, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, &
                     potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ie)
      call check(error, ie(1, 1) == 0_int32, "land never flagged")
      if (allocated(error)) return
      call check(error, ie(2, 1) == ICE_ERR_NONFINITE .and. ie(3, 1) == ICE_ERR_NONFINITE .and. ie(4, 1) == ICE_ERR_NONFINITE &
                 .and. ie(1, 2) == ICE_ERR_NONFINITE .and. ie(2, 2) == ICE_ERR_NONFINITE, "every bad wet cell flagged")
      if (allocated(error)) return
      call check(error, count(ie /= 0_int32) == 6 .and. ie(4, 3) == 2_int32, "only those, prior bits kept")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(potT)) .and. all(ieee_is_finite(rhoa)) .and. all(ieee_is_finite(wind)) &
                 .and. all(ieee_is_finite(uatm)) .and. all(ieee_is_finite(flw)) .and. all(ieee_is_finite(Qa)) &
                 .and. all(ieee_is_finite(fsn)) .and. all(ieee_is_finite(frn)) .and. all(ieee_is_finite(vatm)), &
                 "every output defined (no NaN) incl. land and bad cells")
      if (allocated(error)) return
      call check(error, potT(1, 1) == ICE_ATM_BENIGN_POTT .and. rhoa(2, 1) == ICE_ATM_BENIGN_RHOA .and. wind(1, 1) == 0.0_wp &
                 .and. fsn(1, 2) == 0.0_wp, "benign values on land and bad cells")
      if (allocated(error)) return
      call check(error, wind(3, 3) == 5.0_wp .and. potT(3, 3) == 265.0_wp .and. fsn(3, 3) == 1.0e-5_wp, "good cell")
      if (allocated(error)) return
      call ice_atm_heights(10.0_wp, zl, zs, ok)
      call check(error, ok .and. zl == 10.0_wp .and. zs == 10.0_wp, "heights")
      if (allocated(error)) return
      call ice_atm_heights(2.0_wp, zl, zs, ok)
      call check(error, ok .and. zl == 2.0_wp .and. zs == 10.0_wp, "zlvl follows z_wind, zlvs fixed")
      if (allocated(error)) return
      call ice_atm_heights(nan, zl, zs, ok)
      call check(error,.not. ok .and. zl == 0.0_wp, "NaN height refused")
      if (allocated(error)) return
      call ice_atm_heights(-1.0_wp, zl, zs, ok)
      call check(error,.not. ok, "negative height refused")
   end subroutine test_land_bad

   ! ------------------------------------------------------------------
   ! builder -> ice_step_therm1
   ! ------------------------------------------------------------------

   subroutine run_therm1(nx, ny, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, &
                         fresh, fsens, flat, flwout, fhocn, sx, sy, aicen, vicen, Tsfcn, ice_err, istat)
      !! One `ice_step_therm1` on a small three-category world with the given cell atmosphere (every other
      !! input fixed); returns the merged fluxes, the stress, the state and the cell error word.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: uatm(nx, ny), vatm(nx, ny), wind(nx, ny), potT(nx, ny), Qa(nx, ny), rhoa(nx, ny)
      real(wp), intent(in) :: flw(nx, ny), fsnow(nx, ny), uvel(nx, ny), vvel(nx, ny)
      real(wp), intent(out) :: fresh(nx, ny), fsens(nx, ny), flat(nx, ny), flwout(nx, ny), fhocn(nx, ny), sx(nx, ny), sy(nx, ny)
      real(wp), intent(out) :: aicen(nx, ny, 3), vicen(nx, ny, 3), Tsfcn(nx, ny, 3)
      integer(int32), intent(out) :: ice_err(nx, ny)
      integer, intent(out) :: istat
      integer, parameter :: NC = 3, NL = 4, NS = 1
      type(ice_params_t) :: p
      real(wp), dimension(nx, ny) :: Tf, sst, frzmlt, tau_x, tau_y, fsalt, evap, fswabs, meltt, meltb, congel, snoice
      real(wp), dimension(nx, ny) :: fbot, rside, wlat
      real(wp), dimension(nx, ny, NC) :: aicen_init, vsnon, fswsfcn, fswintn, fswthrun, freshn, fsaltn, fhocnn, evapn
      real(wp), dimension(nx, ny, NC) :: fsensn, flatn, flwoutn, fswabsn, melttn, meltbn, congeln, snoicen, sxn, syn
      integer(int32) :: ierr_n(nx, ny, NC)
      real(wp) :: qicen(nx, ny, NC, NL), sicen(nx, ny, NC, NL), Iswabsn(nx, ny, NC, NL), qsnon(nx, ny, NC, NS)
      real(wp) :: Sswabsn(nx, ny, NC, NS), S(NL + 1)
      integer :: i, j, k, n
      real(wp), parameter :: AR(NC) = [0.3_wp, 0.4_wp, 0.0_wp], HI(NC) = [0.5_wp, 1.5_wp, 0.0_wp]
      real(wp), parameter :: HS(NC) = [0.0_wp, 0.1_wp, 0.0_wp]

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NC, NL, NS)
      p%atm%highfreq = .true.
      call ice_salinity_profile(NL, p%c%saltmax, p%c%min_salin, S)
      Tf = -1.8_wp
      sst = -1.8_wp
      frzmlt = 0.0_wp
      tau_x = 0.02_wp
      tau_y = 0.0_wp
      fresh = 0.0_wp; fsalt = 0.0_wp; fhocn = 0.0_wp; evap = 0.0_wp; fsens = 0.0_wp; flat = 0.0_wp
      flwout = 0.0_wp; fswabs = 0.0_wp; meltt = 0.0_wp; meltb = 0.0_wp; congel = 0.0_wp; snoice = 0.0_wp
      sx = 0.0_wp; sy = 0.0_wp; fbot = 0.0_wp; rside = 0.0_wp; wlat = 0.0_wp
      ice_err = 0_int32
      ierr_n = 0_int32
      vsnon = 0.0_wp; fswsfcn = 0.0_wp; fswintn = 0.0_wp; fswthrun = 0.0_wp; freshn = 0.0_wp; fsaltn = 0.0_wp
      fhocnn = 0.0_wp; evapn = 0.0_wp; fsensn = 0.0_wp; flatn = 0.0_wp; flwoutn = 0.0_wp; fswabsn = 0.0_wp
      melttn = 0.0_wp; meltbn = 0.0_wp; congeln = 0.0_wp; snoicen = 0.0_wp; sxn = 0.0_wp; syn = 0.0_wp
      qicen = 0.0_wp; sicen = 0.0_wp; Iswabsn = 0.0_wp; qsnon = 0.0_wp; Sswabsn = 0.0_wp
      do n = 1, NC
         do j = 1, ny
            do i = 1, nx
               aicen(i, j, n) = AR(n)
               aicen_init(i, j, n) = AR(n)
               vicen(i, j, n) = AR(n)*HI(n)
               vsnon(i, j, n) = AR(n)*HS(n)
               Tsfcn(i, j, n) = merge(-10.0_wp, 0.0_wp, AR(n) > 0.0_wp)
               do k = 1, NL
                  sicen(i, j, n, k) = S(k)
                  if (AR(n) > 0.0_wp) then
                     qicen(i, j, n, k) = ice_q_from_T(-6.0_wp + 4.0_wp*real(k - 1, wp)/real(NL - 1, wp), S(k), p%c, &
                                                      p%thermo%l_brine)
                  end if
               end do
               if (AR(n) > 0.0_wp) qsnon(i, j, n, 1) = ice_qsno_from_T(-5.0_wp, p%c)
            end do
         end do
      end do

      !$acc enter data copyin(p, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
      !$acc&  fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, sx, sy, fbot, rside, &
      !$acc&  wlat, ice_err, aicen_init, aicen, vicen, vsnon, Tsfcn, fswsfcn, fswintn, fswthrun, freshn, fsaltn, fhocnn, &
      !$acc&  evapn, fsensn, flatn, flwoutn, fswabsn, melttn, meltbn, congeln, snoicen, sxn, syn, ierr_n, qicen, sicen, &
      !$acc&  Iswabsn, qsnon, Sswabsn)
      call ice_step_therm1(nx, ny, NC, NL, NS, 1, nx, 1, ny, 900.0_wp, p, .false., 10.0_wp, 10.0_wp, &
                           uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
                           aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                           fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                           fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, &
                           freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                           melttn, meltbn, congeln, snoicen, sxn, syn, sx, sy, fbot, rside, wlat, ierr_n, ice_err, istat)
      !$acc update self(fresh, fsens, flat, flwout, fhocn, sx, sy, aicen, vicen, Tsfcn, ice_err)
      !$acc exit data delete(p, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
      !$acc&  fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, sx, sy, fbot, rside, &
      !$acc&  wlat, ice_err, aicen_init, aicen, vicen, vsnon, Tsfcn, fswsfcn, fswintn, fswthrun, freshn, fsaltn, fhocnn, &
      !$acc&  evapn, fsensn, flatn, flwoutn, fswabsn, melttn, meltbn, congeln, snoicen, sxn, syn, ierr_n, qicen, sicen, &
      !$acc&  Iswabsn, qsnon, Sswabsn)
   end subroutine run_therm1

   subroutine test_therm1(error)
      !! Builder -> `ice_step_therm1` on a 3-category world against the same step with the atmosphere typed in as
      !! LITERALS: two atmospheres (non-unit pressure ratio, non-zero grid angle, columns 1-3 and 4-6) whose
      !! `potT, rhoa, uatm, vatm, wind` were computed independently (Python `math`, 17 digits).  The builder agrees
      !! with the literals to the kernel's accuracy (`rdb_sm_pow` against libm: a few ulp), so the thermodynamics
      !! outputs agree to `TH_RTOL`, the measured spread times ten.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 6, NY = 5
      real(wp), parameter :: TH_RTOL = 1.0e-11_wp
      real(wp), parameter :: LP(2) = [257.3330600324957_wp, 267.25270685323864_wp]
      real(wp), parameter :: LR(2) = [1.3661243519972421_wp, 1.2889133089278122_wp]
      real(wp), parameter :: LU(2) = [3.127739698515642_wp, -4.815229607298926_wp]
      real(wp), parameter :: LV(2) = [-5.168872640946891_wp, -1.3466862399949844_wp]
      real(wp), parameter :: LW(2) = [6.041522986797286_wp, 5.0_wp]
      real(wp), dimension(NX, NY) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, ea, na, ua, uvm, ta, ai
      real(wp), dimension(NX, NY) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, uT, vT
      real(wp), dimension(NX, NY) :: potT_l, rhoa_l, wind_l, uatm_l, vatm_l
      real(wp) :: u(NX + 1, NY), v(NX, NY + 1)
      real(wp), dimension(NX, NY) :: fr1, fs1, fl1, fo1, fh1, sx1, sy1, fr2, fs2, fl2, fo2, fh2, sx2, sy2
      real(wp) :: a1(NX, NY, 3), v1(NX, NY, 3), t1(NX, NY, 3), a2(NX, NY, 3), v2(NX, NY, 3), t2(NX, NY, 3)
      integer(int32) :: e0(NX, NY), e1(NX, NY), e2(NX, NY)
      integer :: i, j, c, st1, st2
      integer(int64) :: st

      st = 4242_int64
      call grid_world(NX, NY, hm, ea, na, ua, uvm, ta)
      ai = 0.7_wp
      do j = 1, NY
         do i = 1, NX
            c = merge(1, 2, i <= 3)
            tas(i, j) = merge(258.3_wp, 266.1_wp, c == 1)
            huss(i, j) = merge(4.2e-4_wp, 7.5e-4_wp, c == 1)
            psl(i, j) = merge(1.0132e5_wp, 9.85e4_wp, c == 1)
            uas(i, j) = merge(5.5_wp, -3.0_wp, c == 1)
            vas(i, j) = merge(-2.5_wp, 4.0_wp, c == 1)
            ang(i, j) = merge(0.6_wp, -1.2_wp, c == 1)
            rlds(i, j) = merge(190.0_wp, 215.0_wp, c == 1)
            prra(i, j) = 0.0_wp
            prsn(i, j) = merge(1.0e-5_wp, 0.0_wp, c == 1)
            potT_l(i, j) = LP(c)
            rhoa_l(i, j) = LR(c)
            uatm_l(i, j) = LU(c)
            vatm_l(i, j) = LV(c)
            wind_l(i, j) = LW(c)
         end do
      end do
      do j = 1, NY
         do i = 1, NX + 1
            u(i, j) = 0.3_wp*uni(st) - 0.15_wp
         end do
      end do
      do j = 1, NY + 1
         do i = 1, NX
            v(i, j) = 0.3_wp*uni(st) - 0.15_wp
         end do
      end do
      e0 = 0_int32
      call d_prepare(NX, NY, 1, NX, 1, NY, 0.0_wp, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, &
                     potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, e0)
      uT = 0.0_wp
      vT = 0.0_wp
      call d_velocity(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, e0)
      call check(error, all(e0 == 0_int32), "no error word from the builders")
      if (allocated(error)) return
      do j = 1, NY
         do i = 1, NX
            if (.not. (near(potT(i, j), potT_l(i, j), 1.0e-14_wp, 0.0_wp) .and. near(rhoa(i, j), rhoa_l(i, j), 1.0e-14_wp, 0.0_wp) &
                       .and. near(uatm(i, j), uatm_l(i, j), 1.0e-14_wp, 1.0e-15_wp) &
                       .and. near(vatm(i, j), vatm_l(i, j), 1.0e-14_wp, 1.0e-15_wp) &
                       .and. near(wind(i, j), wind_l(i, j), 1.0e-14_wp, 0.0_wp))) then
               call test_failed(error, "builder differs from the independently computed literals")
               return
            end if
         end do
      end do
      call check(error, all(potT /= tas) .and. all(abs(vatm) > 1.0_wp), "the case is non-trivial: pressure ratio, rotation")
      if (allocated(error)) return
      call run_therm1(NX, NY, uatm, vatm, wind, potT, Qa, rhoa, flw, fsn, uT, vT, &
                      fr1, fs1, fl1, fo1, fh1, sx1, sy1, a1, v1, t1, e1, st1)
      call run_therm1(NX, NY, uatm_l, vatm_l, wind_l, potT_l, huss, rhoa_l, rlds, prsn, uT, vT, &
                      fr2, fs2, fl2, fo2, fh2, sx2, sy2, a2, v2, t2, e2, st2)
      call check(error, st1 == 0 .and. st2 == 0 .and. all(e1 == 0_int32) .and. all(e2 == 0_int32), "therm1 runs clean")
      if (allocated(error)) return
      if (.not. (all(ieee_is_finite(fs1)) .and. all(ieee_is_finite(sx1)) .and. all(ieee_is_finite(fo1)))) then
         call test_failed(error, "non-finite flux")
         return
      end if
      call check(error, all(sx1(2:NX - 1, 2:NY - 1) /= 0.0_wp), "the stress is produced")
      if (allocated(error)) return
      call cmp_all(error, fr1, fr2, "fresh", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, fs1, fs2, "fsens", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, fl1, fl2, "flat", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, fo1, fo2, "flwout", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, fh1, fh2, "fhocn", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, sx1, sx2, "strair_x_T", TH_RTOL)
      if (allocated(error)) return
      call cmp_all(error, sy1, sy2, "strair_y_T", TH_RTOL)
      if (allocated(error)) return
      call cmp1(error, reshape(v1, [NX*NY*3]), reshape(v2, [NX*NY*3]), "vicen", TH_RTOL)
      if (allocated(error)) return
      call cmp1(error, reshape(t1, [NX*NY*3]), reshape(t2, [NX*NY*3]), "Tsfcn", TH_RTOL)
   end subroutine test_therm1

   subroutine cmp_all(error, a, b, what, rtol)
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: a(:, :), b(:, :), rtol
      character(len=*), intent(in) :: what
      call cmp1(error, reshape(a, [size(a)]), reshape(b, [size(b)]), what, rtol)
   end subroutine cmp_all

   subroutine cmp1(error, a, b, what, rtol)
      !! NaN-proof agreement to `rtol` (relative, with a tiny absolute floor for exact zeros).
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: a(:), b(:), rtol
      character(len=*), intent(in) :: what
      integer :: k
      do k = 1, size(a)
         if (.not. ieee_is_finite(a(k))) then
            call test_failed(error, what//": non-finite")
            return
         end if
         if (.not. (a(k) == b(k) .or. near(a(k), b(k), rtol, 1.0e-300_wp))) then
            call test_failed(error, what//": differs from the literal-input step")
            return
         end if
      end do
   end subroutine cmp1

   ! ------------------------------------------------------------------
   ! device sweep and throughput
   ! ------------------------------------------------------------------

   subroutine fill_atm(nx, ny, st, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn)
      integer, intent(in) :: nx, ny
      integer(int64), intent(inout) :: st
      real(wp), dimension(nx, ny), intent(out) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn
      integer :: i, j
      do j = 1, ny
         do i = 1, nx
            hm(i, j) = merge(0.0_wp, 1.0_wp, uni(st) < 0.15_wp)
            ang(i, j) = 6.0_wp*uni(st) - 3.0_wp
            tas(i, j) = 235.0_wp + 55.0_wp*uni(st)
            huss(i, j) = 3.0e-3_wp*uni(st)
            uas(i, j) = -15.0_wp + 30.0_wp*uni(st)
            vas(i, j) = -15.0_wp + 30.0_wp*uni(st)
            psl(i, j) = 9.5e4_wp + 8.0e3_wp*uni(st)
            rlds(i, j) = 150.0_wp + 150.0_wp*uni(st)
            prra(i, j) = 1.0e-4_wp*uni(st)
            prsn(i, j) = 1.0e-4_wp*uni(st)
         end do
      end do
   end subroutine fill_atm

   subroutine ref_cell(i, j, hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, prsn, nx, ny, potT, rhoa, wind, uatm, vatm, Qa, &
                       flw, fsn, frn, ie)
      !! Host reference of one cell of the builder: the scalar kernel, or the benign values on land.
      integer, intent(in) :: i, j, nx, ny
      real(wp), dimension(nx, ny), intent(in) :: hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, prsn
      real(wp), dimension(nx, ny), intent(inout) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      integer(int32), intent(inout) :: ie(nx, ny)
      logical :: ok
      if (hm(i, j) > 0.5_wp) then
         call ice_atm_cell_inputs(tas(i, j), huss(i, j), uas(i, j), vas(i, j), psl(i, j), rlds(i, j), prra(i, j), &
                                  prsn(i, j), cA(i, j), sA(i, j), 0.0_wp, potT(i, j), rhoa(i, j), wind(i, j), uatm(i, j), &
                                  vatm(i, j), Qa(i, j), flw(i, j), fsn(i, j), frn(i, j), ok)
         if (.not. ok) ie(i, j) = ior(ie(i, j), ICE_ERR_NONFINITE)
      else
         potT(i, j) = ICE_ATM_BENIGN_POTT
         rhoa(i, j) = ICE_ATM_BENIGN_RHOA
         wind(i, j) = 0.0_wp
         uatm(i, j) = 0.0_wp
         vatm(i, j) = 0.0_wp
         Qa(i, j) = 0.0_wp
         flw(i, j) = 0.0_wp
         fsn(i, j) = 0.0_wp
         frn(i, j) = 0.0_wp
      end if
   end subroutine ref_cell

   function corner_mean(a1, w1, a2, w2) result(s)
      !! The two-face corner mean of the dynamics (zero weights give zero).
      real(wp), intent(in) :: a1, w1, a2, w2
      real(wp) :: s
      s = 0.0_wp
      if (w1 + w2 /= 0.0_wp) s = (a1*w1 + a2*w2)/(w1 + w2)
   end function corner_mean

   subroutine test_device(error)
      !! 24 x 18 cells, owned 3..22 x 3..16: every output of every cell against the scalar kernel run on the
      !! host (NaN, non-positive and land cells included), the ring untouched, the velocity and stress helpers.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 24, NY = 18, I0 = 3, I1 = 22, J0 = 3, J1 = 16
      real(wp), dimension(NX, NY) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, ea, na, ua, uvm, ta, ai, cA, sA
      real(wp), dimension(NX, NY) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      real(wp), dimension(NX, NY) :: hp, hr, hw, hu, hv, hq, hf, hs, hn
      real(wp), dimension(NX, NY) :: uT, vT, uTh, vTh, sxE, syN, sxEh, syNh, sxT, syT
      real(wp) :: u(NX + 1, NY), v(NX, NY + 1)
      integer(int32), dimension(NX, NY) :: ie, ieh, ie2
      integer(int64) :: st
      integer :: i, j
      real(wp) :: nan, wt
      real(wp) :: cu(4), cv(4)
      integer :: q, ci(4), cj(4)

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      st = 99_int64
      call fill_atm(NX, NY, st, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn)
      huss(5, 5) = nan
      tas(7, 6) = -3.0_wp
      psl(9, 7) = 0.0_wp
      tas(1, 1) = nan                      ! ring, never read
      hm(11, 8) = 0.0_wp
      tas(11, 8) = nan                     ! land with NaN
      cA = 1.0_wp
      sA = 0.0_wp
      call d_rotation(NX, NY, I0, I1, J0, J1, ang, cA, sA)
      call check(error, all(abs(cA(I0:I1, J0:J1) - cos(ang(I0:I1, J0:J1))) < 4.0e-16_wp*4.0_wp) &
                 .and. all(abs(sA(I0:I1, J0:J1) - sin(ang(I0:I1, J0:J1))) < 4.0e-16_wp*4.0_wp), "setup rotation vs libm")
      if (allocated(error)) return
      ieh = 0_int32
      hp = SENT; hr = SENT; hw = SENT; hu = SENT; hv = SENT; hq = SENT; hf = SENT; hs = SENT; hn = SENT
      do j = J0, J1
         do i = I0, I1
            call ref_cell(i, j, hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, prsn, NX, NY, hp, hr, hw, hu, hv, hq, hf, &
                          hs, hn, ieh)
         end do
      end do
      uvm = 1.0_wp
      do j = 1, NY
         do i = 1, NX
            ea(i, j) = 2.0_wp + uni(st)
            na(i, j) = 2.0_wp + uni(st)
            ta(i, j) = 3.0_wp + uni(st)
            ai(i, j) = uni(st)
            ua(i, j) = 0.5_wp + uni(st)
            if (uni(st) < 0.15_wp) uvm(i, j) = 0.0_wp
            sxT(i, j) = uni(st) - 0.5_wp
            syT(i, j) = uni(st) - 0.5_wp
         end do
      end do
      do j = 1, NY
         do i = 1, NX + 1
            u(i, j) = uni(st) - 0.5_wp
         end do
      end do
      do j = 1, NY + 1
         do i = 1, NX
            v(i, j) = uni(st) - 0.5_wp
         end do
      end do
      uTh = SENT
      vTh = SENT
      sxEh = SENT
      syNh = SENT
      do j = J0, J1
         do i = I0, I1
            uTh(i, j) = 0.0_wp
            vTh(i, j) = 0.0_wp
            if (hm(i, j) > 0.5_wp .and. ai(i, j) > 1.0e-11_wp) then
               ci = [i, i - 1, i, i - 1]
               cj = [j, j, j - 1, j - 1]
               wt = 0.0_wp
               do q = 1, 4
                  cu(q) = corner_mean(u(ci(q) + 1, cj(q)), ea(ci(q), cj(q)), u(ci(q) + 1, cj(q) + 1), ea(ci(q), cj(q) + 1)) &
                          *uvm(ci(q), cj(q))
                  cv(q) = corner_mean(v(ci(q), cj(q) + 1), na(ci(q), cj(q)), v(ci(q) + 1, cj(q) + 1), na(ci(q) + 1, cj(q))) &
                          *uvm(ci(q), cj(q))
                  wt = wt + ua(ci(q), cj(q))
               end do
               uTh(i, j) = (cu(1)*ua(i, j) + cu(2)*ua(i - 1, j) + cu(3)*ua(i, j - 1) + cu(4)*ua(i - 1, j - 1))/wt
               vTh(i, j) = (cv(1)*ua(i, j) + cv(2)*ua(i - 1, j) + cv(3)*ua(i, j - 1) + cv(4)*ua(i - 1, j - 1))/wt
            end if
            sxEh(i, j) = 0.5_wp*(sxT(i, j)*ta(i, j) + sxT(i + 1, j)*ta(i + 1, j))/ea(i, j)
            syNh(i, j) = 0.5_wp*(syT(i, j)*ta(i, j) + syT(i, j + 1)*ta(i, j + 1))/na(i, j)
         end do
      end do
      ie = 0_int32
      ie2 = 0_int32
      potT = SENT; rhoa = SENT; wind = SENT; uatm = SENT; vatm = SENT; Qa = SENT; flw = SENT; fsn = SENT; frn = SENT
      uT = SENT
      vT = SENT
      sxE = SENT
      syN = SENT

      call d_prepare(NX, NY, I0, I1, J0, J1, 0.0_wp, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn, &
                     potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ie)
      call d_velocity(NX, NY, I0, I1, J0, J1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie2)
      call d_faces(NX, NY, I0, I1, J0, J1, sxT, syT, ta, ea, na, sxE, syN)

      call check(error, all(ie == ieh) .and. all(ie2 == 0_int32) .and. count(ie /= 0_int32) == 3, &
                 "error words match the host, 3 bad cells")
      if (allocated(error)) return
      do j = 1, NY
         do i = 1, NX
            if (i >= I0 .and. i <= I1 .and. j >= J0 .and. j <= J1) then
               if (.not. (dcmp(potT(i, j), hp(i, j)) .and. dcmp(rhoa(i, j), hr(i, j)) .and. dcmp(wind(i, j), hw(i, j)) &
                          .and. dcmp(uatm(i, j), hu(i, j)) .and. dcmp(vatm(i, j), hv(i, j)) .and. dcmp(Qa(i, j), hq(i, j)) &
                          .and. dcmp(flw(i, j), hf(i, j)) .and. dcmp(fsn(i, j), hs(i, j)) .and. dcmp(frn(i, j), hn(i, j)) &
                          .and. dcmp(uT(i, j), uTh(i, j)) .and. dcmp(vT(i, j), vTh(i, j)) &
                          .and. dcmp(sxE(i, j), sxEh(i, j)) .and. dcmp(syN(i, j), syNh(i, j)))) then
                  call test_failed(error, "device output differs from the host reference (owned cell)")
                  return
               end if
            else
               if (.not. (potT(i, j) == SENT .and. rhoa(i, j) == SENT .and. wind(i, j) == SENT .and. uatm(i, j) == SENT &
                          .and. vatm(i, j) == SENT .and. Qa(i, j) == SENT .and. flw(i, j) == SENT .and. fsn(i, j) == SENT &
                          .and. frn(i, j) == SENT .and. uT(i, j) == SENT .and. vT(i, j) == SENT .and. sxE(i, j) == SENT &
                          .and. syN(i, j) == SENT)) then
                  call test_failed(error, "ring written")
                  return
               end if
            end if
         end do
      end do
   end subroutine test_device

   logical function dcmp(a, b)
      !! Device vs host: exact on CPU builds; scale-aware on the GPU build; a non-finite never matches.
      real(wp), intent(in) :: a, b
      dcmp = .false.
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) return
      dcmp = a == b .or. near(a, b, 1.0e-12_wp*merge(1.0_wp, 0.0_wp, DRTOL > 0.0_wp), 1.0e-300_wp)
   end function dcmp

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 512, NY = 512, NREP = 20
      real(wp), allocatable, dimension(:, :) :: hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn
      real(wp), allocatable, dimension(:, :) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn
      real(wp), allocatable, dimension(:, :) :: ea, na, ua, uvm, ta, ai, uT, vT, sx, sy, sxE, syN, cA, sA
      real(wp), allocatable :: u(:, :), v(:, :)
      integer(int32), allocatable :: ie(:, :)
      integer(int64) :: st, t0, t1, rate
      integer :: k
      real(wp) :: sec

      allocate (hm(NX, NY), ang(NX, NY), tas(NX, NY), huss(NX, NY), uas(NX, NY), vas(NX, NY), psl(NX, NY), rlds(NX, NY))
      allocate (prra(NX, NY), prsn(NX, NY), potT(NX, NY), rhoa(NX, NY), wind(NX, NY), uatm(NX, NY), vatm(NX, NY))
      allocate (Qa(NX, NY), flw(NX, NY), fsn(NX, NY), frn(NX, NY), ea(NX, NY), na(NX, NY), ua(NX, NY), uvm(NX, NY))
      allocate (cA(NX, NY), sA(NX, NY))
      cA = 1.0_wp
      sA = 0.0_wp
      allocate (ta(NX, NY), ai(NX, NY), uT(NX, NY), vT(NX, NY), sx(NX, NY), sy(NX, NY), sxE(NX, NY), syN(NX, NY))
      allocate (u(NX + 1, NY), v(NX, NY + 1), ie(NX, NY))
      st = 7_int64
      call grid_world(NX, NY, hm, ea, na, ua, uvm, ta)
      call fill_atm(NX, NY, st, hm, ang, tas, huss, uas, vas, psl, rlds, prra, prsn)
      ai = 0.8_wp
      u = 0.2_wp
      v = -0.1_wp
      sx = 0.05_wp
      sy = -0.02_wp
      ie = 0_int32
      potT = 0.0_wp; rhoa = 0.0_wp; wind = 0.0_wp; uatm = 0.0_wp; vatm = 0.0_wp; Qa = 0.0_wp; flw = 0.0_wp
      fsn = 0.0_wp; frn = 0.0_wp; uT = 0.0_wp; vT = 0.0_wp; sxE = 0.0_wp; syN = 0.0_wp
      call d_rotation(NX, NY, 2, NX - 1, 2, NY - 1, ang, cA, sA)
      !$acc enter data copyin(hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, prsn, ie, potT, rhoa, wind, uatm, vatm, Qa, &
      !$acc&  flw, fsn, frn, ea, na, ua, uvm, ta, ai, u, v, uT, vT, sx, sy, sxE, syN)
      call ice_atm_prepare_cell_inputs(NX, NY, 2, NX - 1, 2, NY - 1, 0.0_wp, hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, &
                                       prsn, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ie)
      call system_clock(t0, rate)
      do k = 1, NREP
         call ice_atm_prepare_cell_inputs(NX, NY, 2, NX - 1, 2, NY - 1, 0.0_wp, hm, cA, sA, tas, huss, uas, vas, psl, rlds, &
                                          prra, prsn, potT, rhoa, wind, uatm, vatm, Qa, flw, fsn, frn, ie)
      end do
      call system_clock(t1)
      sec = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,es10.3,a)') "  ice_atm_prepare_cell_inputs: ", real(NX - 2, wp)*real(NY - 2, wp)*NREP/sec, &
         " cells/s (512x512)"
      call system_clock(t0)
      do k = 1, NREP
         call ice_atm_velocity_T(NX, NY, 2, NX - 1, 2, NY - 1, 1.0e-11_wp, u, v, ea, na, ua, uvm, hm, ai, uT, vT, ie)
      end do
      call system_clock(t1)
      sec = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,es10.3,a)') "  ice_atm_velocity_T:          ", real(NX - 2, wp)*real(NY - 2, wp)*NREP/sec, &
         " cells/s (512x512)"
      call system_clock(t0)
      do k = 1, NREP
         call ice_atm_stress_to_faces(NX, NY, 2, NX - 1, 2, NY - 1, sx, sy, ta, ea, na, sxE, syN)
      end do
      call system_clock(t1)
      sec = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,es10.3,a)') "  ice_atm_stress_to_faces:     ", real(NX - 2, wp)*real(NY - 2, wp)*NREP/sec, &
         " cells/s (512x512)"
      !$acc update self(potT, rhoa, wind, uT, sxE, ie)
      !$acc exit data delete(hm, cA, sA, tas, huss, uas, vas, psl, rlds, prra, prsn, ie, potT, rhoa, wind, uatm, vatm, Qa, &
      !$acc&  flw, fsn, frn, ea, na, ua, uvm, ta, ai, u, v, uT, vT, sx, sy, sxE, syN)
      call check(error, all(ieee_is_finite(potT)) .and. all(ieee_is_finite(rhoa)) .and. all(ieee_is_finite(sxE)) &
                 .and. all(ieee_is_finite(uT)), "finite results")
   end subroutine test_throughput

end module test_cice_atm_inputs
