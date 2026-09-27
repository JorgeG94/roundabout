!! Analytic tests for the closed-form Wright (1997) in-situ density
!! integrals of the FV_MOM6 constant-by-layer PGF (`wright_pcm_dpa_intz`,
!! `wright_pcm_dpa_face` — MOM6 `int_density_dz_wright`).
!!
!! 1. layer_integral_matches_quadrature — for a constant-T/S layer the
!!    closed form `dpa = g*int (rho - rho_ref) dz` and its first moment
!!    from the layer top, `intz_dpa`, must equal a composite 5-point
!!    Gauss-Legendre quadrature (64 panels, exact to degree 9 per panel)
!!    of the generic `eos_density_point` at `p = -g*rho0*z`, to round-off.
!!    The sweep covers cold/fresh to warm/salty water, the surface to
!!    6000 m (p ~ 6.1e7 Pa), 1 mm to 6000 m thick layers, an interface
!!    above the datum, and both the anomaly (`rho_ref = rho0`) and the
!!    full-density (`rho_ref = 0`) integrand.
!! 2. face_integral_matches_quadrature — the cross-face Boole rule with
!!    the analytic vertical integral at every lateral sub-column must
!!    equal the same lateral Boole rule over Gauss-Legendre vertical
!!    integrals, with and without non-trivial near-bottom mass weights.
!! 3. boole_path_agrees — the generic 5-point Boole path the analytic one
!!    replaces (`boole_dpa_intz_layer` / `boole_dpa_face_pcm`) agrees with
!!    it to its own (tiny) truncation error: in pressure the Wright density
!!    varies on the scale `p + p0 + lambda/alpha0 ~ 8e8 Pa`, so a layer's
!!    vertical Boole rule is accurate to `(g*rho0*dz/8e8)^6`.
module test_ocean_pgf_wright_analytic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_WRIGHT_97
   use rdb_ocean_pressure_force, only: wright_pcm_dpa_intz, wright_pcm_dpa_face
   use rdb_ocean_pgf_reconstruct, only: boole_dpa_intz_layer, boole_dpa_face_pcm
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_ocean_pgf_wright_analytic_tests

   real(wp), parameter :: RHO0 = 1035.0_wp
   integer, parameter :: N_T = 5, N_S = 3, N_LAY = 10
   real(wp), parameter :: T_SET(N_T) = [-1.9_wp, 2.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      !! Temperatures (degC): freezing to tropical surface.
   real(wp), parameter :: S_SET(N_S) = [32.0_wp, 34.7_wp, 37.0_wp]
      !! Salinities (PSU).
   real(wp), parameter :: E_TOP_SET(N_LAY) = [0.0_wp, 0.0_wp, 1.0_wp, -50.0_wp, &
                                              -500.0_wp, -1000.0_wp, -4000.0_wp, &
                                              -5800.0_wp, 0.0_wp, -5999.0_wp]
      !! Layer-top heights (m): surface, above the datum, deep.
   real(wp), parameter :: DZ_SET(N_LAY) = [1.0e-3_wp, 2.0_wp, 5.0_wp, 10.0_wp, &
                                           100.0_wp, 1000.0_wp, 500.0_wp, &
                                           200.0_wp, 6000.0_wp, 1.0e-3_wp]
      !! Layer thicknesses (m): 1 mm to the whole 6000 m column.
   real(wp), parameter :: TOL_LAYER = 1.0e-13_wp
      !! Round-off bound on the layer integrals, relative to the full
      !! hydrostatic scale `g*rho0*dz` (and `0.5*g*rho0*dz^2`).
   real(wp), parameter :: TOL_ANOM = 1.0e-11_wp
      !! Round-off bound relative to the ANOMALY scale
      !! `g*dz*max(|rho_mean - rho_ref|, 1 kg/m^3)`.
   real(wp), parameter :: TOL_BOOLE = 1.0e-9_wp
      !! Bound on the retired 5-point Boole path's truncation error,
      !! relative to `g*rho0*dz`.

contains

   subroutine collect_ocean_pgf_wright_analytic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("layer_integral_matches_quadrature", test_layer), &
                  new_unittest("face_integral_matches_quadrature", test_face), &
                  new_unittest("boole_path_agrees", test_boole_agrees) &
                  ]
   end subroutine collect_ocean_pgf_wright_analytic_tests

   subroutine make_wright(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%rho0 = RHO0
      eos%is_init = .true.
   end subroutine make_wright

   pure subroutine gl_layer(eos, t, s, e_top, dz, rho_ref, dpa, intz_dpa)
      !! Composite 5-point Gauss-Legendre reference for one constant-T/S
      !! layer: `dpa = g*int_0^dz rho'(s) ds`, `intz_dpa = g*int_0^dz
      !! rho'(s)*(dz - s) ds`, `s` the depth below the layer top and
      !! `rho' = EOS(T, S, g*rho0*(s - e_top)) - rho_ref`.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: t, s, e_top, dz, rho_ref
      real(wp), intent(out) :: dpa, intz_dpa
      integer, parameter :: NPAN = 64
      real(wp), parameter :: X5(5) = [-0.9061798459386640_wp, -0.5384693101056831_wp, &
                                      0.0_wp, 0.5384693101056831_wp, 0.9061798459386640_wp]
      real(wp), parameter :: W5(5) = [0.2369268850561891_wp, 0.4786286704993665_wp, &
                                      0.5688888888888889_wp, 0.4786286704993665_wp, &
                                      0.2369268850561891_wp]
      real(wp) :: hp, sc, sd, r, acc0, acc1
      integer :: n, m
      hp = dz/real(NPAN, wp)
      acc0 = 0.0_wp
      acc1 = 0.0_wp
      do n = 1, NPAN
         sc = (real(n, wp) - 0.5_wp)*hp
         do m = 1, 5
            sd = sc + 0.5_wp*hp*X5(m)
            r = eos_density_point(eos, t, s, GRAVITY*RHO0*(sd - e_top)) - rho_ref
            acc0 = acc0 + W5(m)*r
            acc1 = acc1 + W5(m)*r*(dz - sd)
         end do
      end do
      dpa = GRAVITY*0.5_wp*hp*acc0
      intz_dpa = GRAVITY*0.5_wp*hp*acc1
   end subroutine gl_layer

   subroutine test_layer(error)
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: rho_ref, dpa, intz, dpa_q, intz_q, e_full, e_anom, anom
      real(wp) :: max_full, max_anom
      integer :: it, is, il, ir
      character(len=200) :: msg

      call make_wright(eos)
      max_full = 0.0_wp
      max_anom = 0.0_wp
      do ir = 1, 2
         rho_ref = merge(RHO0, 0.0_wp, ir == 1)
         do il = 1, N_LAY
            do is = 1, N_S
               do it = 1, N_T
                  call wright_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                           RHO0, rho_ref, dpa, intz)
                  call gl_layer(eos, T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                rho_ref, dpa_q, intz_q)
                  anom = max(abs(dpa_q)/(GRAVITY*DZ_SET(il)), 1.0_wp)
                  e_full = max(abs(dpa - dpa_q)/(GRAVITY*RHO0*DZ_SET(il)), &
                               abs(intz - intz_q)/(0.5_wp*GRAVITY*RHO0*DZ_SET(il)**2))
                  e_anom = max(abs(dpa - dpa_q)/(GRAVITY*anom*DZ_SET(il)), &
                               abs(intz - intz_q)/(0.5_wp*GRAVITY*anom*DZ_SET(il)**2))
                  max_full = max(max_full, e_full)
                  max_anom = max(max_anom, e_anom)
                  if (e_full > TOL_LAYER .or. e_anom > TOL_ANOM) then
                     write (msg, '(a,2f7.2,2es11.3,a,2es10.2)') &
                        "layer T,S,e_top,dz=", T_SET(it), S_SET(is), E_TOP_SET(il), &
                        DZ_SET(il), " err full/anom=", e_full, e_anom
                     call check(error, .false., trim(msg))
                     return
                  end if
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3)') "wright analytic layer: max rel err (full scale) ", &
         max_full, ", (anomaly scale) ", max_anom
      call global_logger%info(trim(msg))
   end subroutine test_layer

   pure subroutine face_reference(eos, e_top_l, e_top_r, dz_l, dz_r, t_l, t_r, s_l, s_r, &
                                  hwt_ll, hwt_lr, hwt_rr, hwt_rl, rho_ref, &
                                  dpa_l, dpa_r, dpa_face)
      !! Lateral 5-point Boole over Gauss-Legendre vertical integrals at the
      !! same sub-columns `wright_pcm_dpa_face` samples.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: e_top_l, e_top_r, dz_l, dz_r, t_l, t_r, s_l, s_r
      real(wp), intent(in) :: hwt_ll, hwt_lr, hwt_rr, hwt_rl, rho_ref
      real(wp), intent(out) :: dpa_l, dpa_r, dpa_face
      real(wp), parameter :: BW(5) = [7.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp]
      real(wp) :: wl, wr, wtl, wtr, d, dm, acc
      integer :: m
      call gl_layer(eos, t_l, s_l, e_top_l, dz_l, rho_ref, dpa_l, d)
      call gl_layer(eos, t_r, s_r, e_top_r, dz_r, rho_ref, dpa_r, d)
      acc = BW(1)*dpa_l + BW(5)*dpa_r
      do m = 2, 4
         wr = 0.25_wp*real(m - 1, wp)
         wl = 1.0_wp - wr
         wtl = wl*hwt_ll + wr*hwt_rl
         wtr = wl*hwt_lr + wr*hwt_rr
         call gl_layer(eos, wtl*t_l + wtr*t_r, wtl*s_l + wtr*s_r, &
                       wl*e_top_l + wr*e_top_r, wl*dz_l + wr*dz_r, rho_ref, dm, d)
         acc = acc + BW(m)*dm
      end do
      dpa_face = acc/90.0_wp
   end subroutine face_reference

   subroutine test_face(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCASE = 5
      ! Column pairs: (e_top_l, e_top_r, dz_l, dz_r, t_l, t_r, s_l, s_r)
      real(wp), parameter :: C(8, NCASE) = reshape([ &
                                                   0.0_wp, 0.5_wp, 2.0_wp, 2.5_wp, 25.0_wp, 18.0_wp, 35.5_wp, 34.0_wp, &
                                                   -200.0_wp, -260.0_wp, 150.0_wp, 90.0_wp, 12.0_wp, 9.0_wp, 35.0_wp, 34.6_wp, &
                                                   -3000.0_wp, -3000.0_wp, 800.0_wp, 30.0_wp, 1.5_wp, 2.5_wp, 34.7_wp, 34.9_wp, &
                                                   -5500.0_wp, -5700.0_wp, 500.0_wp, 300.0_wp, -0.5_wp, 1.0_wp, 34.65_wp, 34.7_wp, &
                                                   -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, 4.0_wp, 1.0_wp, 34.9_wp, 34.7_wp], &
                                                   [8, NCASE])
      real(wp) :: hw(4, 2), eos_dummy_ref, dpa_l, dpa_r, ref, ana, scale, err, max_err
      type(eos_t) :: eos
      integer :: n, iw, ir
      character(len=200) :: msg

      call make_wright(eos)
      ! (hwt_ll, hwt_lr, hwt_rr, hwt_rl): plain linear, then mass-weighted.
      hw(:, 1) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
      hw(:, 2) = [0.8_wp, 0.2_wp, 0.65_wp, 0.35_wp]
      max_err = 0.0_wp
      do ir = 1, 2
         eos_dummy_ref = merge(RHO0, 0.0_wp, ir == 1)
         do iw = 1, 2
            do n = 1, NCASE
               call face_reference(eos, C(1, n), C(2, n), C(3, n), C(4, n), &
                                   C(5, n), C(6, n), C(7, n), C(8, n), &
                                   hw(1, iw), hw(2, iw), hw(3, iw), hw(4, iw), &
                                   eos_dummy_ref, dpa_l, dpa_r, ref)
               call wright_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                        C(5, n), C(6, n), C(7, n), C(8, n), &
                                        dpa_l, dpa_r, hw(1, iw), hw(2, iw), hw(3, iw), &
                                        hw(4, iw), RHO0, eos_dummy_ref, ana)
               scale = GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n))
               err = abs(ana - ref)/scale
               max_err = max(max_err, err)
               if (err > TOL_LAYER) then
                  write (msg, '(a,i0,a,i0,a,i0,a,es10.3)') "face case ", n, " weights ", iw, &
                     " rho_ref ", ir, " rel err ", err
                  call check(error, .false., trim(msg))
                  return
               end if
            end do
         end do
      end do
      write (msg, '(a,es10.3)') "wright analytic face: max rel err (full scale) ", max_err
      call global_logger%info(trim(msg))
   end subroutine test_face

   subroutine test_boole_agrees(error)
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: dpa, intz, dpa_b, intz_b, err, max_err, fb, fa
      integer :: it, is, il
      character(len=200) :: msg

      call make_wright(eos)
      max_err = 0.0_wp
      do il = 1, N_LAY
         do is = 1, N_S
            do it = 1, N_T
               call wright_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                        RHO0, RHO0, dpa, intz)
               call boole_dpa_intz_layer(eos, RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                         T_SET(it), T_SET(it), T_SET(it), &
                                         S_SET(is), S_SET(is), S_SET(is), &
                                         .false., dpa_b, intz_b)
               err = max(abs(dpa - dpa_b)/(GRAVITY*RHO0*DZ_SET(il)), &
                         abs(intz - intz_b)/(0.5_wp*GRAVITY*RHO0*DZ_SET(il)**2))
               max_err = max(max_err, err)
            end do
         end do
      end do
      ! One strongly-tilted, thick face through the old and new face rules.
      call wright_pcm_dpa_intz(4.0_wp, 34.9_wp, -10.0_wp, 5990.0_wp, RHO0, RHO0, dpa, intz)
      call wright_pcm_dpa_intz(1.0_wp, 34.7_wp, -4000.0_wp, 2000.0_wp, RHO0, RHO0, dpa_b, intz_b)
      call wright_pcm_dpa_face(-10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, 4.0_wp, 1.0_wp, &
                               34.9_wp, 34.7_wp, dpa, dpa_b, 0.8_wp, 0.2_wp, 0.65_wp, 0.35_wp, &
                               RHO0, RHO0, fa)
      call boole_dpa_face_pcm(eos, RHO0, RHO0, -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, &
                              4.0_wp, 1.0_wp, 34.9_wp, 34.7_wp, dpa, dpa_b, &
                              0.8_wp, 0.2_wp, 0.65_wp, 0.35_wp, fb)
      max_err = max(max_err, abs(fa - fb)/(GRAVITY*RHO0*0.5_wp*(5990.0_wp + 2000.0_wp)))
      write (msg, '(a,es10.3)') "wright analytic vs 5-point Boole: max rel diff ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_BOOLE, trim(msg))
   end subroutine test_boole_agrees

end module test_ocean_pgf_wright_analytic
