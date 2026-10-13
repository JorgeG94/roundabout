!! Tests for the factored Roquet et al. (2015) in-situ-DENSITY (ΔS = 32,
!! "polyTEOS10_55t") kernels of the FV-MOM6 in-situ PGF: the layer
!! integral `roquet_rho_pcm_dpa_intz` and its cross-face twin
!! `roquet_rho_pcm_dpa_face` (`rdb_roquet_rho.inc`, PGF-side), plus the
!! `rdb_eos` entry points `roquet_rho_ts_coeffs` / `roquet_rho_value` that
!! share the same include. This EOS variant is independent of, and NOT
!! bit-identical to, the already-shipped specific-volume `roquet_spv`
!! variant (`test_ocean_pgf_eos_fast.F90`) -- same paper, ΔS = 24 vs 32,
!! density vs specific volume.
!!
!! The dispatch through `compute_fv_mom6_insitu_pcm_impl` /
!! `compute_fv_mom6_reconstruct_impl` is wired (`use_insitu_pcm` admits
!! `EOS_VARIANT_ROQUET_RHO`, both FV-MOM6 branches carry a third loop copy
!! for it), and the configure-time refusal of `eos='roquet_rho'` with
!! `form='fv_mom6'` is lifted (`rdb_config_validate`,
!! `rdb_ocean_setup_lateral`) -- the combination now validates and runs this
!! branch by default. The full-kernel, device-executing coverage of that
!! namelist-driven path lives in `tests/test_ocean_pgf_insitu.F90`
!! (`thermobaric_face_matches_insitu_oracle_roquet_rho` /
!! `rest_partial_steps_zero_pgf_roquet_rho`); the kernels below, including
!! the reconstruct-for-pressure twins, are tested directly against
!! analytic/quadrature oracles instead:
!!
!! 1. rho_layer_matches_quadrature -- `roquet_rho_pcm_dpa_intz` (5-point
!!    Boole, the (T, S) part of the EOS evaluated once) against a composite
!!    5-point Gauss-Legendre quadrature (64 panels) of the generic
!!    `eos_density_point` (eos%variant = ROQUET_RHO), surface to 6000 m
!!    (p ~ 6.1e7 Pa), 1 mm to 6000 m layers.
!! 2. rho_matches_generic_boole -- the factored rule and its face twin
!!    `roquet_rho_pcm_dpa_face` against the generic Boole path they
!!    replace (`boole_dpa_intz_layer` / `boole_dpa_face_pcm` through the
!!    `eos_t` handle), to round-off: the same five densities, weights and
!!    order.
!! 3. rho_pcm_face_matches_quadrature -- `roquet_rho_pcm_dpa_face` against
!!    a 2-D Gauss-Legendre reference (16 panels across the face x 64 in
!!    the vertical), plain linear interpolation only (see
!!    `test_ocean_pgf_eos_fast.F90`'s twin for why mass weighting has no
!!    quadrature to compare against).
!! 4. rho_inc_matches_eos_point -- `roquet_rho_value(T, S, p)` against
!!    `eos_density_point` at the Roquet et al. (2015) Sec. 5 paper check
!!    point and on a 4x4x5 (T, S, p) grid: bitwise (the inc IS the value
!!    half `eos_density_point` calls for this variant).
!! 5. rho_recon_twins_match_generic_and_quadrature -- the reconstruct-for-
!!    pressure twins `roquet_rho_recon_dpa_intz` / `roquet_rho_recon_dpa_face`
!!    (the in-situ-density value from `rdb_roquet_rho.inc`) on PLM and PPM
!!    profiles: against the generic Boole rule (round-off) and against
!!    Gauss-Legendre of the generic EOS along the same profile (the Boole
!!    truncation), mirroring `test_ocean_pgf_eos_fast.F90`'s
!!    `test_recon_twins`.
module test_ocean_pgf_roquet_rho
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_eos, only: eos_t, eos_density_point, roquet_rho_value, roquet_rho_point, &
                      EOS_VARIANT_ROQUET_RHO
   use rdb_ocean_pressure_force, only: roquet_rho_pcm_dpa_intz, roquet_rho_pcm_dpa_face, &
                                       roquet_rho_recon_dpa_intz, roquet_rho_recon_dpa_face
   use rdb_ocean_pgf_reconstruct, only: boole_dpa_intz_layer, boole_dpa_face_pcm, &
                                        boole_dpa_face
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_ocean_pgf_roquet_rho_tests

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
      !! Round-off bound relative to the full hydrostatic scale `g*rho0*dz`
      !! (copied from `test_ocean_pgf_eos_fast.F90`, never loosened).
   real(wp), parameter :: TOL_BOOLE_ROQ = 2.0e-9_wp
      !! Bound on the 5-point Boole truncation error of the Roquet layer
      !! integral, relative to `g*rho0*dz` (same bound as the SpV fit's
      !! `TOL_BOOLE_ROQ` in `test_ocean_pgf_eos_fast.F90` -- the in-situ
      !! density fit is degree-6 in p too, so the same Boole-exactness
      !! argument applies; the actual max is logged, not asserted tighter).
   real(wp), parameter :: TOL_BOOLE_FACE = 1.0e-10_wp
      !! Bound on the cross-face Boole truncation, relative to
      !! `g*rho0*mean(dz)` (copied from `test_ocean_pgf_eos_fast.F90`).
   real(wp), parameter :: TOL_BOOLE_PROF = 1.0e-8_wp
      !! Bound on the in-layer Boole truncation on a PLM / PPM sub-layer
      !! profile, relative to `g*rho0*dz` (copied from
      !! `test_ocean_pgf_eos_fast.F90`, never loosened).
   real(wp), parameter :: TOL_BOOLE_THIN = 1.0e-9_wp
      !! The same for layers up to 1000 m (copied from
      !! `test_ocean_pgf_eos_fast.F90`).

contains

   subroutine collect_ocean_pgf_roquet_rho_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("rho_layer_matches_quadrature", test_rho_quadrature), &
                  new_unittest("rho_matches_generic_boole", test_rho_generic), &
                  new_unittest("rho_pcm_face_matches_quadrature", test_rho_pcm_face_quadrature), &
                  new_unittest("rho_inc_matches_eos_point", test_rho_inc_matches_point), &
                  new_unittest("rho_recon_twins_match_generic_and_quadrature", test_rho_recon_twins) &
                  ]
   end subroutine collect_ocean_pgf_roquet_rho_tests

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
      !! `rho' = EOS(T, S, g*rho0*(s - e_top)) - rho_ref`. Mirrors
      !! `test_ocean_pgf_eos_fast.F90`'s `gl_layer` exactly.
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

   subroutine test_rho_quadrature(error)
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: rho_ref, dpa, intz, dpa_q, intz_q, err, max_err, max_err_thin
      integer :: it, is, il, ir
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_RHO)
      max_err = 0.0_wp
      max_err_thin = 0.0_wp
      do ir = 1, 2
         rho_ref = merge(RHO0, 0.0_wp, ir == 1)
         do il = 1, N_LAY
            do is = 1, N_S
               do it = 1, N_T
                  call roquet_rho_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                               RHO0, rho_ref, dpa, intz)
                  call gl_layer(eos, T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
                                rho_ref, dpa_q, intz_q)
                  err = max(abs(dpa - dpa_q)/(GRAVITY*RHO0*DZ_SET(il)), &
                            abs(intz - intz_q)/(0.5_wp*GRAVITY*RHO0*DZ_SET(il)**2))
                  max_err = max(max_err, err)
                  if (DZ_SET(il) <= 1000.0_wp) max_err_thin = max(max_err_thin, err)
                  if (err > TOL_BOOLE_ROQ) then
                     write (msg, '(a,2f7.2,2es11.3,a,es10.3)') "roquet_rho layer T,S,e_top,dz=", &
                        T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), " rel err ", err
                     call check(error, .false., trim(msg))
                     return
                  end if
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3)') "roquet_rho layer vs 64-panel GL: max rel err ", max_err, &
         " (layers <= 1000 m: ", max_err_thin
      call global_logger%info(trim(msg)//")")
   end subroutine test_rho_quadrature

   subroutine test_rho_generic(error)
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

      call make_eos(eos, EOS_VARIANT_ROQUET_RHO)
      max_err = 0.0_wp
      do il = 1, N_LAY
         do is = 1, N_S
            do it = 1, N_T
               call roquet_rho_pcm_dpa_intz(T_SET(it), S_SET(is), E_TOP_SET(il), DZ_SET(il), &
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
            call roquet_rho_pcm_dpa_intz(C(5, n), C(7, n), C(1, n), C(3, n), RHO0, RHO0, dpa_l, intz)
            call roquet_rho_pcm_dpa_intz(C(6, n), C(8, n), C(2, n), C(4, n), RHO0, RHO0, dpa_r, intz)
            call roquet_rho_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                         C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                         hw(1, iw), hw(2, iw), hw(3, iw), hw(4, iw), &
                                         RHO0, RHO0, fa)
            call boole_dpa_face_pcm(eos, RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                    C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                    hw(1, iw), hw(2, iw), hw(3, iw), hw(4, iw), fb)
            max_err = max(max_err, abs(fa - fb)/(GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n))))
         end do
      end do
      write (msg, '(a,es10.3)') "roquet_rho factored vs generic Boole: max rel diff ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ROUND, trim(msg))
   end subroutine test_rho_generic

   subroutine test_rho_pcm_face_quadrature(error)
      !! The PCM cross-face rule `roquet_rho_pcm_dpa_face` against
      !! Gauss-Legendre in both directions of the generic EOS over the
      !! same sub-columns: the error is the Boole truncation across the
      !! face. Plain linear interpolation only (see
      !! `test_ocean_pgf_eos_fast.F90`'s `pcm_face_matches_quadrature` for
      !! why mass weighting has no quadrature to compare against).
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
      real(wp) :: dpa_l, dpa_r, intz, fa, fq, err, max_err
      integer :: n
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_RHO)
      max_err = 0.0_wp
      do n = 1, NCASE
         ! End points: the columns' own integrals (the rule's w = 0, 1).
         call gl_layer(eos, C(5, n), C(7, n), C(1, n), C(3, n), RHO0, dpa_l, intz)
         call gl_layer(eos, C(6, n), C(8, n), C(2, n), C(4, n), RHO0, dpa_r, intz)
         call roquet_rho_pcm_dpa_face(C(1, n), C(2, n), C(3, n), C(4, n), &
                                      C(5, n), C(6, n), C(7, n), C(8, n), dpa_l, dpa_r, &
                                      1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
                                      RHO0, RHO0, fa)
         fq = gl_face_pcm(eos, C(1, n), C(2, n), C(3, n), C(4, n), &
                          C(5, n), C(6, n), C(7, n), C(8, n))
         err = abs(fa - fq)/(GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n)))
         max_err = max(max_err, err)
      end do
      write (msg, '(a,es10.3)') "roquet_rho PCM face rule vs 2-D Gauss-Legendre: max rel err ", max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_BOOLE_FACE, trim(msg))
   end subroutine test_rho_pcm_face_quadrature

   function gl_face_pcm(eos, e_l, e_r, dz_l, dz_r, t_l, t_r, s_l, s_r) result(dpa_face)
      !! Along-face mean of the layer `dpa` by composite Gauss-Legendre in
      !! BOTH directions, specialised to a PCM (constant-T/S) column pair:
      !! 16 panels x 5 points across the face (fraction `w` from the left
      !! column), each sub-column integrated in the vertical by `gl_layer`
      !! (64 panels), the sub-column height/thickness/T/S plain-linearly
      !! interpolated between the columns -- the PCM specialisation of
      !! `test_ocean_pgf_eos_fast.F90`'s `gl_face` (no PLM/PPM profile
      !! needed; T/S is already constant through the layer).
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: e_l, e_r, dz_l, dz_r, t_l, t_r, s_l, s_r
      real(wp) :: dpa_face
      integer, parameter :: NPAN = 16
      real(wp), parameter :: X5(5) = [-0.9061798459386640_wp, -0.5384693101056831_wp, &
                                      0.0_wp, 0.5384693101056831_wp, 0.9061798459386640_wp]
      real(wp), parameter :: W5(5) = [0.2369268850561891_wp, 0.4786286704993665_wp, &
                                      0.5688888888888889_wp, 0.4786286704993665_wp, &
                                      0.2369268850561891_wp]
      real(wp) :: hp, w, wl, wr, dpa, intz, acc
      integer :: n, m
      hp = 1.0_wp/real(NPAN, wp)
      acc = 0.0_wp
      do n = 1, NPAN
         do m = 1, 5
            w = (real(n, wp) - 0.5_wp)*hp + 0.5_wp*hp*X5(m)
            wr = w
            wl = 1.0_wp - w
            call gl_layer(eos, wl*t_l + wr*t_r, wl*s_l + wr*s_r, &
                          wl*e_l + wr*e_r, wl*dz_l + wr*dz_r, RHO0, dpa, intz)
            acc = acc + W5(m)*dpa
         end do
      end do
      dpa_face = 0.5_wp*hp*acc
   end function gl_face_pcm

   subroutine test_rho_inc_matches_point(error)
      !! The inc (`roquet_rho_value`, what `eos_density_point` and the PGF
      !! kernels evaluate) against the INDEPENDENT `roquet_rho_point` (the
      !! value-and-derivative routine `eos_specvol_derivs` uses, its own
      !! Horner code): two implementations of the same 55-term table must
      !! agree, so a dropped or mistyped coefficient in either shows up
      !! here. Also locks the `eos_density_point` wiring to the inc bitwise.
      !! The paper check value itself is pinned in `test_ocean_eos_roquet_rho`.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: T_GRID(4) = [-1.9_wp, 2.0_wp, 15.0_wp, 30.0_wp]
      real(wp), parameter :: S_GRID(4) = [0.0_wp, 30.0_wp, 34.7_wp, 37.0_wp]
      real(wp), parameter :: P_GRID(5) = [0.0_wp, 1.0e5_wp, 1.0e6_wp, 2.0e7_wp, 6.0e7_wp]
      type(eos_t) :: eos
      real(wp) :: v_inc, v_point, v_indep, d_dt, d_ds, max_diff, max_rel
      integer :: it, is, ip
      character(len=200) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_RHO)
      max_diff = 0.0_wp
      max_rel = 0.0_wp
      do ip = 1, size(P_GRID)
         do is = 1, size(S_GRID)
            do it = 1, size(T_GRID)
               v_inc = roquet_rho_value(T_GRID(it), S_GRID(is), P_GRID(ip))
               v_point = eos_density_point(eos, T_GRID(it), S_GRID(is), P_GRID(ip))
               max_diff = max(max_diff, abs(v_inc - v_point))
               call roquet_rho_point(T_GRID(it), S_GRID(is), P_GRID(ip), v_indep, d_dt, d_ds)
               max_rel = max(max_rel, abs(v_inc - v_indep)/abs(v_indep))
            end do
         end do
      end do
      write (msg, '(a,es10.3)') "roquet_rho_value vs eos_density_point: max abs diff ", max_diff
      call global_logger%info(trim(msg))
      call check(error, max_diff == 0.0_wp, trim(msg))
      if (allocated(error)) return
      write (msg, '(a,es10.3)') "roquet_rho_value vs roquet_rho_point: max rel diff ", max_rel
      call global_logger%info(trim(msg))
      call check(error, max_rel <= 1.0e-14_wp, trim(msg))
   end subroutine test_rho_inc_matches_point

   subroutine column_edges(t_t, t_b, t_m, s_t, s_b, s_m)
      !! A PLM/PPM-like sub-layer profile: edges off the mean. Copied
      !! verbatim from `test_ocean_pgf_eos_fast.F90`.
      real(wp), intent(in) :: t_m, s_m
      real(wp), intent(out) :: t_t, t_b, s_t, s_b
      t_t = t_m + 0.7_wp
      t_b = t_m - 0.4_wp
      s_t = s_m - 0.05_wp
      s_b = s_m + 0.08_wp
   end subroutine column_edges

   pure subroutine gl_layer_profile(eos, t_t, t_b, t_m, s_t, s_b, s_m, parabolic, &
                                    e_top, dz, rho_ref, dpa, intz_dpa)
      !! `gl_layer` for a layer whose T/S follow the reconstruct-for-pressure
      !! sub-layer profile: `q(wt) = wt*q_t + (1-wt)*q_b + q6*wt*(1-wt)`,
      !! `wt = 1 - s/dz` the fraction of the way UP from the bottom edge,
      !! `q6 = 3*(2*q_m - (q_t + q_b))` when `parabolic` (else 0). Mirrors
      !! `test_ocean_pgf_eos_fast.F90`'s helper of the same name.
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

   pure function gl_face_profile(eos, e_l, e_r, dz_l, dz_r, tt, sv, hw, parabolic) result(dpa_face)
      !! Along-face mean of the layer `dpa` by composite Gauss-Legendre in
      !! BOTH directions, specialised to a PLM/PPM profile pair -- the
      !! `gl_face` of `test_ocean_pgf_eos_fast.F90`. `tt(1:3)` / `sv(1:3)`
      !! are the left column's T / S (top, bottom, mean), `tt(4:6)` /
      !! `sv(4:6)` the right column's; `hw = [1, 0, 1, 0]` is plain linear
      !! interpolation.
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
   end function gl_face_profile

   subroutine test_rho_recon_twins(error)
      !! The reconstruct-for-pressure Roquet in-situ-density twins
      !! (`roquet_rho_recon_dpa_intz` / `roquet_rho_recon_dpa_face`) on PLM
      !! and PPM profiles: against the generic Boole rule through the
      !! `eos_t` handle (round-off) and against Gauss-Legendre of the
      !! generic EOS along the same profile (the Boole truncation), surface
      !! to 6000 m, 1 mm to 6000 m layers, plus tilted faces. Mirrors
      !! `test_ocean_pgf_eos_fast.F90`'s `test_recon_twins`, Roquet-only.
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
      integer :: ip, il, is, it, n
      logical :: parabolic
      character(len=240) :: msg

      call make_eos(eos, EOS_VARIANT_ROQUET_RHO)
      err_gen = 0.0_wp
      err_q = 0.0_wp
      err_q_thin = 0.0_wp
      err_face_gen = 0.0_wp
      err_face_q = 0.0_wp
      do ip = 1, 2
         parabolic = (ip == 2)
         do il = 1, N_LAY
            do is = 1, N_S
               do it = 1, N_T
                  call column_edges(t_t, t_b, T_SET(it), s_t, s_b, S_SET(is))
                  call roquet_rho_recon_dpa_intz(RHO0, RHO0, E_TOP_SET(il), DZ_SET(il), &
                                                 t_t, t_b, T_SET(it), s_t, s_b, S_SET(is), &
                                                 parabolic, dpa, intz)
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
            call roquet_rho_recon_dpa_face(RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                           t_t, t_b, C(5, n), t_t2, t_b2, C(6, n), &
                                           s_t, s_b, C(7, n), s_t2, s_b2, C(8, n), &
                                           dpa_l, dpa_r, parabolic, fa)
            call boole_dpa_face(eos, RHO0, RHO0, C(1, n), C(2, n), C(3, n), C(4, n), &
                                t_t, t_b, C(5, n), t_t2, t_b2, C(6, n), &
                                s_t, s_b, C(7, n), s_t2, s_b2, C(8, n), &
                                dpa_l, dpa_r, parabolic, fb)
            fq = gl_face_profile(eos, C(1, n), C(2, n), C(3, n), C(4, n), &
                                 [t_t, t_b, C(5, n), t_t2, t_b2, C(6, n)], &
                                 [s_t, s_b, C(7, n), s_t2, s_b2, C(8, n)], HW_LIN, parabolic)
            scale = GRAVITY*RHO0*0.5_wp*(C(3, n) + C(4, n))
            err_face_gen = max(err_face_gen, abs(fa - fb)/scale)
            err_face_q = max(err_face_q, abs(fa - fq)/scale)
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') &
         "roquet_rho recon twins vs generic Boole: layer ", err_gen, ", face ", err_face_gen, &
         "; vs Gauss-Legendre: layer ", err_q, " (<= 1000 m: ", err_q_thin, "), face ", err_face_q
      call global_logger%info(trim(msg))
      call check(error, err_gen <= TOL_ROUND .and. err_face_gen <= TOL_ROUND, trim(msg))
      if (allocated(error)) return
      call check(error, err_q <= TOL_BOOLE_PROF .and. err_q_thin <= TOL_BOOLE_THIN .and. &
                 err_face_q <= TOL_BOOLE_FACE, trim(msg))
   end subroutine test_rho_recon_twins

end module test_ocean_pgf_roquet_rho
