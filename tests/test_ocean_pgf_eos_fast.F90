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
!! 6. recon_twins_match_generic_and_quadrature -- the reconstruct-for-
!!    pressure per-EOS twins, Wright (`boole_dpa_intz_layer_wright` /
!!    `boole_dpa_face_wright`) and Roquet (`roquet_recon_dpa_intz` /
!!    `roquet_recon_dpa_face`, the SpV value from `rdb_roquet_spv.inc`), on
!!    PLM and PPM profiles: against the generic Boole rule (round-off) and
!!    against Gauss-Legendre of the generic EOS along the same profile, 64
!!    panels in the vertical and 16 across a face (the Boole truncation),
!!    surface to 6000 m (p to 6.1e7 Pa), 1 mm to 6000 m layers.
!! 7. pcm_face_matches_quadrature -- the PCM cross-face rules
!!    (`roquet_pcm_dpa_face`, `wright_pcm_dpa_face`) against the same 2-D
!!    Gauss-Legendre reference, with and without mass weighting.
!! 8. edges_layer_match_column_reference -- the per-layer PLM / PPM edge
!!    helpers of the reconstruct kernel's 3-D Pass 0 (`plm_edges_layer`,
!!    `ppm_edges_layer`) against the retired per-column routines (kept here,
!!    verbatim, as the oracle), 1- to 12-layer columns with vanished layers.
module test_ocean_pgf_eos_fast
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      eos_density_specvol_derivs, roquet_spv_value, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV
   use rdb_ocean_pressure_force, only: roquet_pcm_dpa_intz, roquet_pcm_dpa_face, &
                                       wright_pcm_dpa_face, &
                                       boole_dpa_intz_layer_wright, boole_dpa_face_wright, &
                                       roquet_recon_dpa_intz, roquet_recon_dpa_face, &
                                       plm_edges_layer, ppm_edges_layer
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_ocean_pgf_reconstruct, only: boole_dpa_intz_layer, boole_dpa_face, &
                                        boole_dpa_face_pcm
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
   real(wp), parameter :: TOL_ULP_CANCEL = 16.0_wp*epsilon(1.0_wp)
      !! The same, for d(SV)/dT, a sum of polynomial terms that CANCEL near
      !! the temperature of maximum density: nvfortran -fast reassociates
      !! the two compilations differently and measured 9 ulp (2.07e-15) of
      !! d(SV)/dT at T = -1.9 degC, where it is 2.6e-8, against 1.1e-16 on
      !! gfortran.  Measured against DSVDT_SCALE below, not its own value.
   real(wp), parameter :: DSVDT_SCALE = 1.0e-7_wp
      !! Physical scale of d(SV)/dT (m3/kg/K): alpha ~ 1e-4 /K times
      !! SV ~ 9.7e-4 m3/kg.  A relative error measured against d(SV)/dT
      !! itself blows up where it passes through zero.

   real(wp), parameter :: TOL_BOOLE_PROF = 1.0e-8_wp
      !! Bound on the in-layer Boole truncation on a PLM / PPM sub-layer
      !! profile, relative to `g*rho0*dz`: the profile's T/S variation
      !! (1.1 K, 0.13 PSU across the test layers) makes the integrand a
      !! high-degree polynomial in depth even where the pressure is not.
      !! Measured (gfortran 15.1): 4.5e-9, the 6000 m column as ONE layer.
   real(wp), parameter :: TOL_BOOLE_THIN = 1.0e-9_wp
      !! The same for layers up to 1000 m.  Measured 1.8e-10.
   real(wp), parameter :: TOL_BOOLE_FACE = 1.0e-10_wp
      !! Bound on the cross-face Boole truncation, relative to
      !! `g*rho0*mean(dz)`, over the tilted test faces up to a
      !! 5990 m / 2000 m pair.  Measured 2.3e-11 (reconstructed profiles),
      !! 9.8e-12 (PCM).

contains

   subroutine collect_ocean_pgf_eos_fast_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("roquet_layer_matches_quadrature", test_roquet_quadrature), &
                  new_unittest("roquet_matches_generic_boole", test_roquet_generic), &
                  new_unittest("value_and_fused_entry_points", test_entry_points), &
                  new_unittest("wright_recon_twin_matches_generic", test_wright_twin), &
                  new_unittest("face_end_points_are_the_columns", test_face_end_points), &
                  new_unittest("recon_twins_match_generic_and_quadrature", test_recon_twins), &
                  new_unittest("pcm_face_matches_quadrature", test_pcm_face_quadrature), &
                  new_unittest("edges_layer_match_column_reference", test_edges_layer) &
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
      real(wp) :: max_err_dt
      integer :: it, is, ip, iv
      character(len=200) :: msg

      max_err = 0.0_wp
      max_err_dt = 0.0_wp
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
                  err = max(abs(rho_a - rho_b)/abs(rho_a), abs(ds_a - ds_b)/abs(ds_a))
                  max_err = max(max_err, err)
                  max_err_dt = max(max_err_dt, &
                                   abs(dt_a - dt_b)/max(abs(dt_a), DSVDT_SCALE))
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3)') "value-only / fused EOS vs separate calls: max rel diff ", &
         max_err, "; d(SV)/dT vs its scale ", max_err_dt
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ULP .and. max_err_dt <= TOL_ULP_CANCEL, trim(msg))
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

   pure subroutine gl_layer_profile(eos, t_t, t_b, t_m, s_t, s_b, s_m, parabolic, &
                                    e_top, dz, rho_ref, dpa, intz_dpa)
      !! `gl_layer` for a layer whose T/S follow the reconstruct-for-pressure
      !! sub-layer profile: `q(wt) = wt*q_t + (1-wt)*q_b + q6*wt*(1-wt)`,
      !! `wt = 1 - s/dz` the fraction of the way UP from the bottom edge,
      !! `q6 = 3*(2*q_m - (q_t + q_b))` when `parabolic` (else 0).
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: t_t, t_b, t_m, s_t, s_b, s_m
      logical, intent(in) :: parabolic
      real(wp), intent(in) :: e_top, dz, rho_ref
      real(wp), intent(out) :: dpa, intz_dpa
      integer, parameter :: NPAN = 64
      real(wp), parameter :: X5(5) = [-0.9061798459386640_wp, -0.5384693101056831_wp, &
                                      0.0_wp, 0.5384693101056831_wp, 0.9061798459386640_wp]
      real(wp), parameter :: W5(5) = [0.2369268850561891_wp, 0.4786286704993665_wp, &
                                      0.5688888888888889_wp, 0.4786286704993665_wp, &
                                      0.2369268850561891_wp]
      real(wp) :: hp, sc, sd, r, acc0, acc1, t6, s6, wt, tq, sq
      integer :: n, m
      t6 = 0.0_wp
      s6 = 0.0_wp
      if (parabolic) then
         t6 = 3.0_wp*(2.0_wp*t_m - (t_t + t_b))
         s6 = 3.0_wp*(2.0_wp*s_m - (s_t + s_b))
      end if
      hp = dz/real(NPAN, wp)
      acc0 = 0.0_wp
      acc1 = 0.0_wp
      do n = 1, NPAN
         sc = (real(n, wp) - 0.5_wp)*hp
         do m = 1, 5
            sd = sc + 0.5_wp*hp*X5(m)
            wt = 1.0_wp - sd/dz
            tq = wt*t_t + (1.0_wp - wt)*t_b + t6*wt*(1.0_wp - wt)
            sq = wt*s_t + (1.0_wp - wt)*s_b + s6*wt*(1.0_wp - wt)
            r = eos_density_point(eos, tq, sq, GRAVITY*RHO0*(sd - e_top)) - rho_ref
            acc0 = acc0 + W5(m)*r
            acc1 = acc1 + W5(m)*r*(dz - sd)
         end do
      end do
      dpa = GRAVITY*0.5_wp*hp*acc0
      intz_dpa = GRAVITY*0.5_wp*hp*acc1
   end subroutine gl_layer_profile

   pure function gl_face(eos, e_l, e_r, dz_l, dz_r, tt, sv, hw, parabolic) result(dpa_face)
      !! Along-face mean of the layer `dpa` by composite Gauss-Legendre in
      !! BOTH directions: 16 panels x 5 points across the face (fraction
      !! `w` from the left column), each sub-column integrated in the
      !! vertical by `gl_layer_profile` (64 panels).  The sub-column at `w`
      !! interpolates the interface height and thickness linearly and the
      !! T/S edge + mean triples with the mass-weighted fractions
      !! `wL = (1-w)*hw(1) + w*hw(4)`, `wR = (1-w)*hw(2) + w*hw(3)`
      !! (`hw = [1, 0, 1, 0]`: plain linear interpolation) -- the
      !! sub-column definition of the Boole face rules.  `tt(1:3)` /
      !! `sv(1:3)` are the left column's T / S (top, bottom, mean),
      !! `tt(4:6)` / `sv(4:6)` the right column's.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: e_l, e_r, dz_l, dz_r
      real(wp), intent(in) :: tt(6), sv(6), hw(4)
      logical, intent(in) :: parabolic
      real(wp) :: dpa_face
      integer, parameter :: NPAN = 16
      real(wp), parameter :: X5(5) = [-0.9061798459386640_wp, -0.5384693101056831_wp, &
                                      0.0_wp, 0.5384693101056831_wp, 0.9061798459386640_wp]
      real(wp), parameter :: W5(5) = [0.2369268850561891_wp, 0.4786286704993665_wp, &
                                      0.5688888888888889_wp, 0.4786286704993665_wp, &
                                      0.2369268850561891_wp]
      real(wp) :: hp, w, wl, wr, fl, fr, dpa, intz, acc, q(6)
      integer :: n, m
      hp = 1.0_wp/real(NPAN, wp)
      acc = 0.0_wp
      do n = 1, NPAN
         do m = 1, 5
            w = (real(n, wp) - 0.5_wp)*hp + 0.5_wp*hp*X5(m)
            wr = w
            wl = 1.0_wp - w
            fl = wl*hw(1) + wr*hw(4)
            fr = wl*hw(2) + wr*hw(3)
            q(1:3) = fl*tt(1:3) + fr*tt(4:6)
            q(4:6) = fl*sv(1:3) + fr*sv(4:6)
            call gl_layer_profile(eos, q(1), q(2), q(3), q(4), q(5), q(6), parabolic, &
                                  wl*e_l + wr*e_r, wl*dz_l + wr*dz_r, RHO0, dpa, intz)
            acc = acc + W5(m)*dpa
         end do
      end do
      dpa_face = 0.5_wp*hp*acc
   end function gl_face

   subroutine test_recon_twins(error)
      !! The reconstruct-for-pressure per-EOS twins (Wright
      !! `boole_dpa_intz_layer_wright` / `boole_dpa_face_wright`, Roquet
      !! `roquet_recon_dpa_intz` / `roquet_recon_dpa_face`) on PLM and PPM
      !! profiles: against the generic Boole rule through the `eos_t` handle
      !! (round-off) and against Gauss-Legendre of the generic EOS along the
      !! same profile (the Boole truncation), surface to 6000 m, 1 mm to
      !! 6000 m layers, plus tilted faces.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCASE = 5
      real(wp), parameter :: C(8, NCASE) = reshape([ &
                                                   0.0_wp, 0.5_wp, 2.0_wp, 2.5_wp, 25.0_wp, 18.0_wp, 35.5_wp, 34.0_wp, &
                                                   -200.0_wp, -260.0_wp, 150.0_wp, 90.0_wp, 12.0_wp, 9.0_wp, 35.0_wp, 34.6_wp, &
                                                   -3000.0_wp, -3000.0_wp, 800.0_wp, 30.0_wp, 1.5_wp, 2.5_wp, 34.7_wp, 34.9_wp, &
                                                   -5500.0_wp, -5700.0_wp, 500.0_wp, 300.0_wp, -0.5_wp, 1.0_wp, 34.65_wp, 34.7_wp, &
                                                   -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, 4.0_wp, 1.0_wp, 34.9_wp, 34.7_wp], &
                                                   [8, NCASE])
      real(wp), parameter :: HW_LIN(4) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
      type(eos_t) :: eos
      real(wp) :: t_t, t_b, s_t, s_b, t_t2, t_b2, s_t2, s_b2
      real(wp) :: dpa, intz, dpa_b, intz_b, dpa_q, intz_q, dpa_l, dpa_r, fa, fb, fq
      real(wp) :: scale, err_gen, err_q, err_q_thin, err_face_gen, err_face_q
      integer :: iv, ip, il, is, it, n
      logical :: parabolic, roquet
      character(len=240) :: msg

      err_gen = 0.0_wp
      err_q = 0.0_wp
      err_q_thin = 0.0_wp
      err_face_gen = 0.0_wp
      err_face_q = 0.0_wp
      do iv = 1, 2
         roquet = (iv == 2)
         call make_eos(eos, merge(EOS_VARIANT_ROQUET_SPV, EOS_VARIANT_WRIGHT_97, roquet))
         do ip = 1, 2
            parabolic = (ip == 2)
            do il = 1, N_LAY
               do is = 1, N_S
                  do it = 1, N_T
                     call column_edges(t_t, t_b, T_SET(it), s_t, s_b, S_SET(is))
                     if (roquet) then
                        call roquet_recon_dpa_intz(RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                                   t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                                   parabolic, dpa, intz)
                     else
                        call boole_dpa_intz_layer_wright(RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                                         t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                                         parabolic, dpa, intz)
                     end if
                     call boole_dpa_intz_layer(eos, RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                               t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                               parabolic, dpa_b, intz_b)
                     call gl_layer_profile(eos, t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                           parabolic, E_TOP_SET(il), DZ_SET(il), RHO0, &
                                           dpa_q, intz_q)
                     scale = GRAVITY*RHO0*DZ_SET(il)
                     err_gen = max(err_gen, abs(dpa - dpa_b)/scale, &
                                   abs(intz - intz_b)/(0.5_wp*scale*DZ_SET(il)))
                     err_q = max(err_q, abs(dpa - dpa_q)/scale, &
                                 abs(intz - intz_q)/(0.5_wp*scale*DZ_SET(il)))
                     if (DZ_SET(il) <= 1000.0_wp) then
                        err_q_thin = max(err_q_thin, abs(dpa - dpa_q)/scale, &
                                         abs(intz - intz_q)/(0.5_wp*scale*DZ_SET(il)))
                     end if
                  end do
               end do
            end do
            ! Faces: the five column pairs, each with a PLM/PPM-like profile.
            do n = 1, NCASE
               call column_edges(t_t, t_b, C(5, n), s_t, s_b, C(7, n))
               call column_edges(t_t2, t_b2, C(6, n), s_t2, s_b2, C(8, n))
               call boole_dpa_intz_layer(eos, RHO0, RHO0, C(1, n), C(3, n), t_t, t_b, C(5, n), &
                                         s_t, s_b, C(7, n), parabolic, dpa_l, intz)
               call boole_dpa_intz_layer(eos, RHO0, RHO0, C(2, n), C(4, n), t_t2, t_b2, C(6, n), &
                                         s_t2, s_b2, C(8, n), parabolic, dpa_r, intz)
               if (roquet) then
                  call roquet_recon_dpa_face(RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                             t_t, t_b, C(5, n), t_t2, t_b2, C(6, n), &
                                             s_t, s_b, C(7, n), s_t2, s_b2, C(8, n), &
                                             dpa_l, dpa_r, parabolic, fa)
               else
                  call boole_dpa_face_wright(RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                             t_t, t_b, C(5, n), t_t2, t_b2, C(6, n), &
                                             s_t, s_b, C(7, n), s_t2, s_b2, C(8, n), &
                                             dpa_l, dpa_r, parabolic, fa)
               end if
               call boole_dpa_face(eos, RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                   t_t, t_b, C(5, n), t_t2, t_b2, C(6, n), &
                                   s_t, s_b, C(7, n), s_t2, s_b2, C(8, n), &
                                   dpa_l, dpa_r, parabolic, fb)
               fq = gl_face(eos, C(1, n), C(2, n), C(3, n), C(4, n), &
                            [t_t, t_b, C(5, n), t_t2, t_b2, C(6, n)], &
                            [s_t, s_b, C(7, n), s_t2, s_b2, C(8, n)], HW_LIN, parabolic)
               scale = GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n))
               err_face_gen = max(err_face_gen, abs(fa - fb)/scale)
               err_face_q = max(err_face_q, abs(fa - fq)/scale)
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') &
         "recon twins vs generic Boole: layer ", err_gen, ", face ", err_face_gen, &
         "; vs Gauss-Legendre: layer ", err_q, " (<= 1000 m: ", err_q_thin, "), face ", err_face_q
      call global_logger%info(trim(msg))
      call check(error, err_gen <= TOL_ROUND .and. err_face_gen <= TOL_ROUND, trim(msg))
      if (allocated(error)) return
      call check(error, err_q <= TOL_BOOLE_PROF .and. err_q_thin <= TOL_BOOLE_THIN .and. &
                 err_face_q <= TOL_BOOLE_FACE, trim(msg))
   end subroutine test_recon_twins

   subroutine test_pcm_face_quadrature(error)
      !! The PCM cross-face rules (`roquet_pcm_dpa_face`, and the Wright
      !! closed-form twin `wright_pcm_dpa_face` alongside) against Gauss-
      !! Legendre in both directions of the generic EOS over the same
      !! sub-columns: the error is the Boole truncation across the face.
      !! Plain linear interpolation only: under MOM6 mass weighting the
      !! rule's end points are the columns' OWN integrals while its interior
      !! sub-columns carry the mass-weighted T/S, so the sub-column is not
      !! one smooth curve in `w` and there is no quadrature to compare with
      !! (`roquet_matches_generic_boole` covers mass weighting against the
      !! generic rule).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCASE = 5
      real(wp), parameter :: C(8, NCASE) = reshape([ &
                                                   0.0_wp, 0.5_wp, 2.0_wp, 2.5_wp, 25.0_wp, 18.0_wp, 35.5_wp, 34.0_wp, &
                                                   -200.0_wp, -260.0_wp, 150.0_wp, 90.0_wp, 12.0_wp, 9.0_wp, 35.0_wp, 34.6_wp, &
                                                   -3000.0_wp, -3000.0_wp, 800.0_wp, 30.0_wp, 1.5_wp, 2.5_wp, 34.7_wp, 34.9_wp, &
                                                   -5500.0_wp, -5700.0_wp, 500.0_wp, 300.0_wp, -0.5_wp, 1.0_wp, 34.65_wp, 34.7_wp, &
                                                   -10.0_wp, -4000.0_wp, 5990.0_wp, 2000.0_wp, 4.0_wp, 1.0_wp, 34.9_wp, 34.7_wp], &
                                                   [8, NCASE])
      type(eos_t) :: eos
      real(wp), parameter :: HW(4) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
      real(wp) :: dpa_l, dpa_r, intz, fa, fq, err, max_err
      integer :: iv, n
      character(len=200) :: msg

      max_err = 0.0_wp
      do iv = 1, 2
         call make_eos(eos, merge(EOS_VARIANT_ROQUET_SPV, EOS_VARIANT_WRIGHT_97, iv == 2))
         do n = 1, NCASE
            ! End points: the columns' own integrals (the rule's w = 0, 1).
            call gl_layer(eos, C(5, n), C(7, n), C(1, n), C(3, n), RHO0, dpa_l, intz)
            call gl_layer(eos, C(6, n), C(8, n), C(2, n), C(4, n), RHO0, dpa_r, intz)
            if (iv == 2) then
               call roquet_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                        C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                        HW(1), HW(2), HW(3), HW(4), &
                                        RHO0, RHO0, fa)
            else
               call wright_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                        C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                        HW(1), HW(2), HW(3), HW(4), &
                                        RHO0, RHO0, fa)
            end if
            fq = gl_face(eos, C(1, n), C(2, n), C(3, n), C(4, n), &
                         [C(5, n), C(5, n), C(5, n), C(6, n), C(6, n), C(6, n)], &
                         [C(7, n), C(7, n), C(7, n), C(8, n), C(8, n), C(8, n)], &
                         HW, .false.)
            err = abs(fa - fq)/(GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n)))
            max_err = max(max_err, err)
         end do
      end do
      write (msg, '(a,es10.3)') "PCM face rules vs 2-D Gauss-Legendre: max rel err ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_BOOLE_FACE, trim(msg))
   end subroutine test_pcm_face_quadrature

   subroutine test_edges_layer(error)
      !! The per-layer PLM / PPM edge helpers the reconstruct kernel's 3-D
      !! Pass 0 calls (`plm_edges_layer`, `ppm_edges_layer`) against the
      !! retired per-column routines they replace (`ref_plm_edges_column`,
      !! `ref_ppm_edges_column` below, the verbatim bodies), over columns of
      !! 1 to 12 layers with thicknesses from vanished (1e-4 m) to 800 m and
      !! monotone, non-monotone and extremal T/S profiles.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_SET(8) = [1, 2, 3, 4, 5, 6, 9, 12]
      integer, parameter :: NSEED = 40
      real(wp), allocatable :: h(:), s(:), t(:), st_r(:), sb_r(:), tt_r(:), tb_r(:)
      real(wp) :: st, sb, tt, tb, err, max_err, scale
      integer :: in, nz, iseed, k, km2, km1, kp1, kp2, ip
      integer :: lcg
      character(len=200) :: msg

      max_err = 0.0_wp
      lcg = 12345
      do in = 1, size(NZ_SET)
         nz = NZ_SET(in)
         allocate (h(nz), s(nz), t(nz), st_r(nz), sb_r(nz), tt_r(nz), tb_r(nz))
         do iseed = 1, NSEED
            do k = 1, nz
               h(k) = rand01(lcg)
               if (h(k) < 0.15_wp) then
                  h(k) = 1.0e-4_wp                ! vanished
               else
                  h(k) = 800.0_wp*h(k)**3 + 0.5_wp
               end if
               s(k) = 34.0_wp + 2.0_wp*rand01(lcg)
               t(k) = -1.9_wp + 30.0_wp*rand01(lcg)
            end do
            if (mod(iseed, 3) == 0) then
               ! Monotone (stably stratified) profile.
               do k = 1, nz
                  t(k) = 2.0_wp + 1.5_wp*real(k, wp)
                  s(k) = 35.0_wp - 0.05_wp*real(k, wp)
               end do
            end if
            scale = max(maxval(abs(s)), maxval(abs(t)))
            do ip = 1, 2
               if (ip == 1) then
                  call ref_plm_edges_column(nz, h, s, st_r, sb_r)
                  call ref_plm_edges_column(nz, h, t, tt_r, tb_r)
               else
                  call ref_ppm_edges_column(nz, h, s, st_r, sb_r)
                  call ref_ppm_edges_column(nz, h, t, tt_r, tb_r)
               end if
               do k = 1, nz
                  km2 = max(k - 2, 1)
                  km1 = max(k - 1, 1)
                  kp1 = min(k + 1, nz)
                  kp2 = min(k + 2, nz)
                  if (ip == 1) then
                     call plm_edges_layer(k, nz, h(km1), h(k), h(kp1), s(km1), s(k), s(kp1), st, sb)
                     call plm_edges_layer(k, nz, h(km1), h(k), h(kp1), t(km1), t(k), t(kp1), tt, tb)
                  else
                     call ppm_edges_layer(k, nz, h(km2), h(km1), h(k), h(kp1), h(kp2), &
                                          s(km2), s(km1), s(k), s(kp1), s(kp2), &
                                          t(km2), t(km1), t(k), t(kp1), t(kp2), st, sb, tt, tb)
                  end if
                  err = max(abs(st - st_r(k)), abs(sb - sb_r(k)), abs(tt - tt_r(k)), &
                            abs(tb - tb_r(k)))/scale
                  max_err = max(max_err, err)
               end do
            end do
         end do
         deallocate (h, s, t, st_r, sb_r, tt_r, tb_r)
      end do
      write (msg, '(a,es10.3)') "per-layer PLM/PPM edges vs the per-column routines: max rel diff ", &
         max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ULP_CANCEL, trim(msg))
   end subroutine test_edges_layer

   function rand01(state) result(x)
      !! Park-Miller minimal-standard LCG in [0, 1): deterministic columns.
      integer, intent(inout) :: state
      real(wp) :: x
      state = int(mod(16807_int64*int(state, int64), 2147483647_int64))
      x = real(state, wp)/2147483647.0_wp
   end function rand01

   ! ---- The retired per-column edge routines (verbatim bodies of
   ! `rdb_ocean_pgf_reconstruct :: plm_edges_column / ppm_edges_column /
   ! boundary_edges_linear` before the per-layer split), kept as the oracle.

   pure subroutine ref_boundary_edges_linear(h_self, h_nbr, q_self, dq_up, q_t, q_b)
      real(wp), intent(in)  :: h_self, h_nbr, q_self, dq_up
      real(wp), intent(out) :: q_t, q_b
      real(wp), parameter :: H_TINY = 1.0e-30_wp
      real(wp) :: d
      d = dq_up*h_self/max(h_self + h_nbr, H_TINY)
      d = sign(min(abs(d), abs(dq_up)), d)
      q_t = q_self + d
      q_b = q_self - d
   end subroutine ref_boundary_edges_linear

   pure subroutine ref_plm_edges_column(nz, h, q, q_t, q_b)
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz), q(nz)
      real(wp), intent(out) :: q_t(nz), q_b(nz)
      real(wp) :: slp(nz)
      real(wp) :: h_l, h_c, h_r, sig_c, sig_l, sig_r, slp_max
      real(wp) :: e_t, e_b, q_lo, q_hi
      integer  :: k
      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if
      slp(1) = 0.0_wp
      slp(nz) = 0.0_wp
      do k = 2, nz - 1
         h_l = h(k - 1)
         h_c = h(k)
         h_r = h(k + 1)
         sig_l = q(k) - q(k - 1)
         sig_r = q(k + 1) - q(k)
         if (sig_l*sig_r <= 0.0_wp) then
            slp(k) = 0.0_wp
         else
            sig_c = 2.0_wp*(q(k + 1) - q(k - 1))*h_c/(h_l + 2.0_wp*h_c + h_r)
            slp_max = 2.0_wp*min(abs(sig_l), abs(sig_r))
            slp(k) = sign(min(abs(sig_c), slp_max), sig_c)
         end if
      end do
      call ref_boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call ref_boundary_edges_linear(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                                     q_t(nz), q_b(nz))
      do k = 2, nz - 1
         e_t = q(k) + 0.5_wp*slp(k)
         e_b = q(k) - 0.5_wp*slp(k)
         q_lo = min(q(k), q(k + 1))
         q_hi = max(q(k), q(k + 1))
         q_t(k) = max(q_lo, min(q_hi, e_t))
         q_lo = min(q(k), q(k - 1))
         q_hi = max(q(k), q(k - 1))
         q_b(k) = max(q_lo, min(q_hi, e_b))
      end do
   end subroutine ref_plm_edges_column

   pure subroutine ref_ppm_edges_column(nz, h, q, q_t, q_b)
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz), q(nz)
      real(wp), intent(out) :: q_t(nz), q_b(nz)
      real(wp) :: edge(nz)
      real(wp) :: q_lo, q_hi, ql, qr, qm, dq, dq_l, dq_r, q6
      real(wp) :: h0, h1, h2, h3, hf, h_sum
      real(wp) :: h01, h12, h23, h012, h123, h0123
      real(wp) :: f1, f2, f3, et1, et2, et3
      real(wp), parameter :: H_NEGLECT = 1.0e-30_wp
      real(wp), parameter :: H_MIN_FRAC = 1.0e-5_wp
      integer  :: k
      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if
      if (nz == 2) then
         call ref_boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
         call ref_boundary_edges_linear(h(2), h(1), q(2), q(2) - q(1), q_t(2), q_b(2))
         return
      end if
      do k = 2, nz - 2
         h0 = h(k - 1)
         h1 = h(k)
         h2 = h(k + 1)
         h3 = h(k + 2)
         h_sum = h0 + h1 + h2 + h3
         if (h0 + h1 <= 0.0_wp .or. h1 + h2 <= 0.0_wp .or. h2 + h3 <= 0.0_wp) then
            hf = H_MIN_FRAC*max(H_NEGLECT, h_sum)
            h0 = max(h0, hf)
            h1 = max(h1, hf)
            h2 = max(h2, hf)
            h3 = max(h3, hf)
         end if
         h01 = h0 + h1
         h12 = h1 + h2
         h23 = h2 + h3
         h012 = h0 + h1 + h2
         h123 = h1 + h2 + h3
         h0123 = h0 + h1 + h2 + h3
         f1 = h01*h23/h12
         f2 = h2*q(k) + h1*q(k + 1)
         f3 = 1.0_wp/h012 + 1.0_wp/h123
         et1 = f1*f2*f3
         et2 = (h2*h23/(h012*h01))*((h0 + 2.0_wp*h1)*q(k) - h1*q(k - 1))
         et3 = (h1*h01/(h123*h23))*((2.0_wp*h2 + h3)*q(k + 1) - h2*q(k + 2))
         edge(k) = (et1 + et2 + et3)/h0123
      end do
      edge(1) = (q(1)*h(2) + q(2)*h(1))/(h(1) + h(2))
      edge(nz - 1) = (q(nz - 1)*h(nz) + q(nz)*h(nz - 1))/(h(nz - 1) + h(nz))
      call ref_boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call ref_boundary_edges_linear(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                                     q_t(nz), q_b(nz))
      do k = 2, nz - 1
         qm = q(k)
         ql = edge(k - 1)
         qr = edge(k)
         q_lo = min(q(k - 1), q(k), q(k + 1))
         q_hi = max(q(k - 1), q(k), q(k + 1))
         ql = max(q_lo, min(q_hi, ql))
         qr = max(q_lo, min(q_hi, qr))
         dq = qr - ql
         dq_l = qm - ql
         dq_r = qr - qm
         if (dq_l*dq_r <= 0.0_wp) then
            ql = qm
            qr = qm
         else
            q6 = 6.0_wp*qm - 3.0_wp*(ql + qr)
            if (abs(q6) > abs(dq)) then
               if (q6*dq > 0.0_wp) then
                  ql = 3.0_wp*qm - 2.0_wp*qr
               else
                  qr = 3.0_wp*qm - 2.0_wp*ql
               end if
            end if
         end if
         q_b(k) = ql
         q_t(k) = qr
      end do
   end subroutine ref_ppm_edges_column

end module test_ocean_pgf_eos_fast
