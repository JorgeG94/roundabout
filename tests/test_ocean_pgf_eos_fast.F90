!! Tests for the fast per-EOS density integrals of the FV_MOM6 PGF and the
!! value-only / fused EOS entry points.
!!
!! 1. roquet_layer_matches_quadrature -- the factored Roquet SpV layer
!!    integral `roquet_pcm_dpa_intz` (5-point Boole of `1/SV`, the (T, S)
!!    part of the EOS evaluated once) against a composite 5-point
!!    Gauss-Legendre quadrature (64 panels) of the generic
!!    `eos_density_point`, surface to 6000 m (p ~ 6.1e7 Pa), 1 mm to
!!    6000 m layers.  `1/SV(p)` is rational in depth, so there is no closed
!!    form; the bound is the Boole rule's own truncation error.
!! 2. roquet_matches_generic_boole -- the factored rule and its face twin
!!    `roquet_pcm_dpa_face` against the generic Boole path they replace
!!    (`boole_dpa_intz_layer` / `boole_dpa_face_pcm` through the `eos_t`
!!    handle), to round-off: the same five densities, weights and order.
!! 3. value_and_fused_entry_points -- `roquet_spv_value` is the `sv` of
!!    `roquet_spv_point`, and `eos_density_specvol_derivs` returns what
!!    `eos_density_point` + `eos_specvol_derivs` return, for every variant.
!! 4. wright_recon_twin_matches_generic -- the reconstruct-for-pressure
!!    Wright twins (`boole_dpa_intz_layer_wright` / `boole_dpa_face_wright`,
!!    density inline, no handle) against the generic rule on PLM and PPM
!!    profiles, to round-off.
!! 5. face_end_points_are_the_columns -- `boole_dpa_face` with the columns'
!!    own `dpa` as its `w = 0, 1` end points equals the full 5-sub-column
!!    rule it replaced (whose end sub-columns ARE the two columns).
module test_ocean_pgf_eos_fast
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      eos_density_specvol_derivs, roquet_spv_value, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV
   use rdb_ocean_pressure_force, only: roquet_pcm_dpa_intz, roquet_pcm_dpa_face
   use rdb_ocean_pgf_reconstruct, only: boole_dpa_intz_layer, boole_dpa_face, &
                                        boole_dpa_face_pcm, boole_dpa_intz_layer_wright, &
                                        boole_dpa_face_wright
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_ocean_pgf_eos_fast_tests

   real(wp), parameter :: RHO0 = 1035.0_wp
   integer, parameter :: N_T = 5, N_S = 3, N_LAY = 10
   real(wp), parameter :: T_SET(N_T) = [-1.9_wp, 2.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      !! Potential temperatures (degC): freezing to tropical surface.
   real(wp), parameter :: S_SET(N_S) = [32.0_wp, 34.7_wp, 37.0_wp]
      !! Practical salinities (PSU).
   real(wp), parameter :: E_TOP_SET(N_LAY) = [0.0_wp, 0.0_wp, 1.0_wp, -50.0_wp, &
                                              -500.0_wp, -1000.0_wp, -4000.0_wp, &
                                              -5800.0_wp, 0.0_wp, -5999.0_wp]
      !! Layer-top heights (m): surface, above the datum, deep.
   real(wp), parameter :: DZ_SET(N_LAY) = [1.0e-3_wp, 2.0_wp, 5.0_wp, 10.0_wp, &
                                           100.0_wp, 1000.0_wp, 500.0_wp, &
                                           200.0_wp, 6000.0_wp, 1.0e-3_wp]
      !! Layer thicknesses (m): 1 mm to the whole 6000 m column.
   real(wp), parameter :: TOL_ROUND = 1.0e-13_wp
      !! Round-off bound relative to the full hydrostatic scale `g*rho0*dz`.
   real(wp), parameter :: TOL_BOOLE_ROQ = 2.0e-9_wp
      !! Bound on the 5-point Boole truncation error of the Roquet layer
      !! integral, relative to `g*rho0*dz`.  Measured (gfortran 15.1):
      !! 8.0e-10 for the whole 6000 m column as ONE layer, 1.1e-13 for every
      !! layer up to 1000 m.
   real(wp), parameter :: TOL_ULP = 4.0_wp*epsilon(1.0_wp)
      !! Relative bound for "the same expression, compiled twice" -- equal
      !! but for an FMA contraction choice.

contains

   subroutine collect_ocean_pgf_eos_fast_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("roquet_layer_matches_quadrature", test_roquet_quadrature), &
                  new_unittest("roquet_matches_generic_boole", test_roquet_generic), &
                  new_unittest("value_and_fused_entry_points", test_entry_points), &
                  new_unittest("wright_recon_twin_matches_generic", test_wright_twin), &
                  new_unittest("face_end_points_are_the_columns", test_face_end_points) &
                  ]
   end subroutine collect_ocean_pgf_eos_fast_tests

   subroutine make_eos(eos, variant)
      type(eos_t), intent(out) :: eos
      integer, intent(in) :: variant
      eos%variant = variant
      eos%rho0 = RHO0
      eos%alpha_T = 0.17_wp
      eos%beta_S = 0.78_wp
      eos%is_init = .true.
   end subroutine make_eos

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

   subroutine test_roquet_quadrature(error)
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: rho_ref, dpa, intz, dpa_q, intz_q, err, max_err, max_err_thin
      integer :: it, is, il, ir
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_SPV)
      max_err = 0.0_wp
      max_err_thin = 0.0_wp
      do ir = 1, 2
         rho_ref = merge(RHO0, 0.0_wp, ir == 1)
         do il = 1, N_LAY
            do is = 1, N_S
               do it = 1, N_T
                  call roquet_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                           RHO0, rho_ref, dpa, intz)
                  call gl_layer(eos, T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                rho_ref, dpa_q, intz_q)
                  err = max(abs(dpa - dpa_q)/(GRAVITY*RHO0*DZ_SET(il)), &
                            abs(intz - intz_q)/(0.5_wp*GRAVITY*RHO0*DZ_SET(il)**2))
                  max_err = max(max_err, err)
                  if (DZ_SET(il) <= 1000.0_wp) max_err_thin = max(max_err_thin, err)
                  if (err > TOL_BOOLE_ROQ) then
                     write (msg, '(a,2f7.2,2es11.3,a,es10.3)') "roquet layer T,S,e_top,dz=", &
                        T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), " rel err ", err
                     call check(error, .false., trim(msg))
                     return
                  end if
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3)') "roquet layer vs 64-panel GL: max rel err ", max_err, &
         " (layers <= 1000 m: ", max_err_thin
      call global_logger%info(trim(msg)//")")
   end subroutine test_roquet_quadrature

   subroutine test_roquet_generic(error)
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
      type(eos_t) :: eos
      real(wp) :: hw(4, 2), dpa, intz, dpa_b, intz_b, dpa_l, dpa_r, fa, fb, err, max_err
      integer :: it, is, il, n, iw
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_SPV)
      max_err = 0.0_wp
      do il = 1, N_LAY
         do is = 1, N_S
            do it = 1, N_T
               call roquet_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
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
      hw(:, 1) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
      hw(:, 2) = [0.8_wp, 0.2_wp, 0.65_wp, 0.35_wp]
      do iw = 1, 2
         do n = 1, NCASE
            call roquet_pcm_dpa_intz(C(5, n), C(7, n), C(1, n), C(3, n), RHO0, RHO0, dpa_l, intz)
            call roquet_pcm_dpa_intz(C(6, n), C(8, n), C(2, n), C(4, n), RHO0, RHO0, dpa_r, intz)
            call roquet_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                     C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                     hw(1, iw), hw(2, iw), hw(3, iw), hw(4, iw), &
                                     RHO0, RHO0, fa)
            call boole_dpa_face_pcm(eos, RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                    C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                    hw(1, iw), hw(2, iw), hw(3, iw), hw(4, iw), fb)
            max_err = max(max_err, abs(fa - fb)/(GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n))))
         end do
      end do
      write (msg, '(a,es10.3)') "roquet factored vs generic Boole: max rel diff ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ROUND, trim(msg))
   end subroutine test_roquet_generic

   subroutine test_entry_points(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: VARS(3) = [EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                                       EOS_VARIANT_ROQUET_SPV]
      real(wp), parameter :: P_SET(4) = [0.0_wp, 1.0e6_wp, 2.0e7_wp, 6.0e7_wp]
      type(eos_t) :: eos
      real(wp) :: sv_v, sv_p, d1, d2, rho_a, rho_b, dt_a, ds_a, dt_b, ds_b, err, max_err
      integer :: it, is, ip, iv
      character(len=200) :: msg

      max_err = 0.0_wp
      do ip = 1, size(P_SET)
         do is = 1, N_S
            do it = 1, N_T
               sv_v = roquet_spv_value(T_SET(it), S_SET(is), P_SET(ip))
               call make_eos(eos, EOS_VARIANT_ROQUET_SPV)
               rho_b = eos_density_point(eos, T_SET(it), S_SET(is), P_SET(ip))
               call eos_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), d1, d2)
               sv_p = 1.0_wp/rho_b
               max_err = max(max_err, abs(sv_v - sv_p)/sv_p)
               do iv = 1, size(VARS)
                  call make_eos(eos, VARS(iv))
                  rho_a = eos_density_point(eos, T_SET(it), S_SET(is), P_SET(ip))
                  call eos_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), dt_a, ds_a)
                  call eos_density_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), &
                                                  rho_b, dt_b, ds_b)
                  err = max(abs(rho_a - rho_b)/abs(rho_a), abs(dt_a - dt_b)/abs(dt_a), &
                            abs(ds_a - ds_b)/abs(ds_a))
                  max_err = max(max_err, err)
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3)') "value-only / fused EOS vs separate calls: max rel diff ", &
         max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ULP, trim(msg))
   end subroutine test_entry_points

   subroutine column_edges(t_t, t_b, t_m, s_t, s_b, s_m)
      !! A PLM/PPM-like sub-layer profile: edges off the mean.
      real(wp), intent(in) :: t_m, s_m
      real(wp), intent(out) :: t_t, t_b, s_t, s_b
      t_t = t_m + 0.7_wp
      t_b = t_m - 0.4_wp
      s_t = s_m - 0.05_wp
      s_b = s_m + 0.08_wp
   end subroutine column_edges

   subroutine test_wright_twin(error)
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: t_t, t_b, s_t, s_b, t_t2, t_b2, s_t2, s_b2
      real(wp) :: dpa, intz, dpa_b, intz_b, dpa_l, dpa_r, fa, fb, err, max_err
      integer :: it, is, il, ip
      logical :: parabolic
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_WRIGHT_97)
      max_err = 0.0_wp
      do ip = 1, 2
         parabolic = (ip == 2)
         do il = 1, N_LAY
            do is = 1, N_S
               do it = 1, N_T
                  call column_edges(t_t, t_b, T_SET(it), s_t, s_b, S_SET(is))
                  call boole_dpa_intz_layer_wright(RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                                   t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                                   parabolic, dpa, intz)
                  call boole_dpa_intz_layer(eos, RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                            t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                            parabolic, dpa_b, intz_b)
                  err = max(abs(dpa - dpa_b)/(GRAVITY*RHO0*DZ_SET(il)), &
                            abs(intz - intz_b)/(0.5_wp*GRAVITY*RHO0*DZ_SET(il)**2))
                  max_err = max(max_err, err)
               end do
            end do
         end do
         ! A tilted, thick face.
         call column_edges(t_t, t_b, 4.0_wp, s_t, s_b, 34.9_wp)
         call column_edges(t_t2, t_b2, 1.0_wp, s_t2, s_b2, 34.7_wp)
         call boole_dpa_intz_layer(eos, RHO0, RHO0, -10.0_wp, 5990.0_wp, t_t, t_b, 4.0_wp, &
                                   s_t, s_b, 34.9_wp, parabolic, dpa_l, intz)
         call boole_dpa_intz_layer(eos, RHO0, RHO0, -4000.0_wp, 2000.0_wp, t_t2, t_b2, 1.0_wp, &
                                   s_t2, s_b2, 34.7_wp, parabolic, dpa_r, intz)
         call boole_dpa_face_wright(RHO0, RHO0, -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, &
                                    t_t, t_b, 4.0_wp, t_t2, t_b2, 1.0_wp, &
                                    s_t, s_b, 34.9_wp, s_t2, s_b2, 34.7_wp, &
                                    dpa_l, dpa_r, parabolic, fa)
         call boole_dpa_face(eos, RHO0, RHO0, -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, &
                             t_t, t_b, 4.0_wp, t_t2, t_b2, 1.0_wp, &
                             s_t, s_b, 34.9_wp, s_t2, s_b2, 34.7_wp, &
                             dpa_l, dpa_r, parabolic, fb)
         max_err = max(max_err, abs(fa - fb)/(GRAVITY*RHO0*0.5_wp*(5990.0_wp + 2000.0_wp)))
      end do
      write (msg, '(a,es10.3)') "wright recon twin vs generic Boole: max rel diff ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ROUND, trim(msg))
   end subroutine test_wright_twin

   subroutine test_face_end_points(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: BW(5) = [7.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp]
      type(eos_t) :: eos
      real(wp) :: t_t, t_b, s_t, s_b, t_t2, t_b2, s_t2, s_b2
      real(wp) :: dpa_l, dpa_r, dpa_m, intz, full, reuse, wl, wr, acc, err, max_err
      integer :: iv, ip, m
      logical :: parabolic
      character(len=200) :: msg

      max_err = 0.0_wp
      do iv = 1, 2
         call make_eos(eos, merge(EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV, iv == 1))
         do ip = 1, 2
            parabolic = (ip == 2)
            call column_edges(t_t, t_b, 12.0_wp, s_t, s_b, 35.0_wp)
            call column_edges(t_t2, t_b2, 9.0_wp, s_t2, s_b2, 34.6_wp)
            ! The retired rule: all five sub-columns integrated.
            acc = 0.0_wp
            do m = 1, 5
               wr = 0.25_wp*real(m - 1, wp)
               wl = 1.0_wp - wr
               call boole_dpa_intz_layer(eos, RHO0, RHO0, wl*(-200.0_wp) + wr*(-260.0_wp), &
                                         wl*150.0_wp + wr*90.0_wp, &
                                         wl*t_t + wr*t_t2, wl*t_b + wr*t_b2, &
                                         wl*12.0_wp + wr*9.0_wp, &
                                         wl*s_t + wr*s_t2, wl*s_b + wr*s_b2, &
                                         wl*35.0_wp + wr*34.6_wp, parabolic, dpa_m, intz)
               acc = acc + BW(m)*dpa_m
            end do
            full = acc/90.0_wp
            ! The shipped rule: end points are the columns' own integrals.
            call boole_dpa_intz_layer(eos, RHO0, RHO0, -200.0_wp, 150.0_wp, t_t, t_b, 12.0_wp, &
                                      s_t, s_b, 35.0_wp, parabolic, dpa_l, intz)
            call boole_dpa_intz_layer(eos, RHO0, RHO0, -260.0_wp, 90.0_wp, t_t2, t_b2, 9.0_wp, &
                                      s_t2, s_b2, 34.6_wp, parabolic, dpa_r, intz)
            call boole_dpa_face(eos, RHO0, RHO0, -200.0_wp, -260.0_wp, 150.0_wp, 90.0_wp, &
                                t_t, t_b, 12.0_wp, t_t2, t_b2, 9.0_wp, &
                                s_t, s_b, 35.0_wp, s_t2, s_b2, 34.6_wp, &
                                dpa_l, dpa_r, parabolic, reuse)
            err = abs(full - reuse)/(GRAVITY*RHO0*0.5_wp*(150.0_wp + 90.0_wp))
            max_err = max(max_err, err)
         end do
      end do
      write (msg, '(a,es10.3)') "face end-point reuse vs 5-sub-column rule: max rel diff ", &
         max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ROUND, trim(msg))
   end subroutine test_face_end_points

end module test_ocean_pgf_eos_fast
