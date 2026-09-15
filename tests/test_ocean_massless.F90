!! Unit + analytical tests for the D4 massless-merge column helper
!! (`rdb_massless`) and its kappa-shear end-to-end wiring.
!!
!! Golden values come from the verified clean-room prototype
!! `local_archive/prototypes/d4_massless_prototype.py` (run its
!! `--report` worked example to regenerate).  The prototype is 0-based
!! Python; this test is 1-based Fortran.  TRANSLATION RULE for the maps:
!!   kc_fortran(k) = kc_python[k-1] + 1     (layer/interface -> merged index)
!!   nzc is identical (a count, not an index)
!!   the bed sentinel: python kc[nz] = nzc_count (= nzc, 0-based bed iface)
!!     -> fortran kc(nz+1) = nzc + 1 (1-based bed merged interface in
!!        qc(1:nzc+1))
!!   kf is a fraction (no index translation), kf_fortran(k) = kf_python[k-1]
!!
!! Ordering: all helpers are LOCAL SURFACE-DOWN (local k=1 = surface,
!! k=nz = bed; interface K=1 surface .. K=nz+1 bed).
module test_ocean_massless
   use testdrive, only: new_unittest, unittest_type, error_type, check
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, H_DIV_EPS
#endif
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_massless, only: massless_build_maps, massless_merge_fields, &
                           massless_interp_back
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t, kappa_shear_compute
   use rdb_eos, only: EOS_VARIANT_LINEAR
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 128
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=128).
#endif

   public :: collect_ocean_massless_tests

   integer, parameter :: NZL = NZ_STACK_MAX
   integer, parameter :: NZLI = NZ_STACK_MAX + 1
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0_KS = 1025.0_wp
   real(wp), parameter :: ALPHA_LIN_KS = RHO0_KS*2.0e-4_wp
   real(wp), parameter :: BETA_LIN_KS = RHO0_KS*7.6e-4_wp
   real(wp), parameter :: F_CORIOLIS = 1.0e-4_wp

contains

   subroutine collect_ocean_massless_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("massless_golden_maps", test_golden_maps), &
                  new_unittest("massless_identity_I1", test_identity), &
                  new_unittest("massless_conservation_I2", test_conservation), &
                  new_unittest("massless_kshear_e2e_merge", test_kshear_merge), &
                  new_unittest("massless_kshear_e2e_identity", test_kshear_identity) &
                  ]
   end subroutine collect_ocean_massless_tests

   ! ------------------------------------------------------------------
   ! T1 — golden maps: the prototype's nz=8 worked example.
   ! ------------------------------------------------------------------
   subroutine test_golden_maps(error)
      !! Prototype --report: nz=8, massless at python k=2 and k=5
      !! (0-based) = fortran k=3 and k=6 (1-based surface-down).
      !! Hardcoded expectations are the 1-based translation of the
      !! prototype maps (see module header TRANSLATION RULE).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: h(NZL), u(NZL), v(NZL), t(NZL), s(NZL)
      real(wp) :: hc(NZL), uc(NZL), vc(NZL), tc(NZL), sc(NZL)
      real(wp) :: kf(NZL + 1)
      integer :: kc(NZL + 1)
      integer :: nzc, k
      ! 1-based expected maps (prototype 0-based + 1):
      integer, parameter :: EXP_KC(NZ + 1) = [1, 2, 2, 3, 4, 4, 5, 6, 7]
      integer, parameter :: EXP_NZC = 6
      ! kf nonzero only at fortran K=3 and K=6 (prototype 0-based K=2,5)
      real(wp), parameter :: EXP_KF3 = 0.99999250005624951_wp
      real(wp), parameter :: EXP_KF6 = 0.99999062508788972_wp
      checks: block
         h = 0.0_wp; u = 0.0_wp; v = 0.0_wp; t = 0.0_wp; s = 0.0_wp
         h(1:NZ) = [12.0_wp, 10.0_wp, H_VANISHED*0.5_wp, 15.0_wp, 8.0_wp, &
                    H_VANISHED*0.5_wp, 6.0_wp, 9.0_wp]
         u(1:NZ) = [0.1_wp, 0.2_wp, 0.0_wp, 0.15_wp, 0.05_wp, 0.0_wp, &
                    -0.1_wp, -0.05_wp]
         v(1:NZ) = [0.0_wp, 0.05_wp, 0.0_wp, -0.1_wp, 0.1_wp, 0.0_wp, &
                    0.08_wp, 0.02_wp]
         t(1:NZ) = [22.0_wp, 20.0_wp, 19.5_wp, 18.0_wp, 16.0_wp, 15.5_wp, &
                    14.0_wp, 12.0_wp]
         s(1:NZ) = [35.0_wp, 35.2_wp, 35.3_wp, 35.5_wp, 35.6_wp, 35.65_wp, &
                    35.7_wp, 35.8_wp]

         call massless_build_maps(h, NZ, H_VANISHED, nzc, hc, kc, kf)

         call check(error, nzc == EXP_NZC, "golden: nzc must be 6")
         if (allocated(error)) exit checks
         do k = 1, NZ + 1
            call check(error, kc(k) == EXP_KC(k), "golden: kc mismatch")
            if (allocated(error)) exit checks
         end do
         if (allocated(error)) exit checks

         ! hc: [12, 10+7.5e-5, 15, 8+7.5e-5, 6, 9]
         call check(error, hc(1) == 12.0_wp, "golden: hc(1)")
         if (allocated(error)) exit checks
         call check(error, abs(hc(2) - (10.0_wp + H_VANISHED*0.5_wp)) < 1.0e-15_wp, &
                    "golden: hc(2)")
         if (allocated(error)) exit checks
         call check(error, hc(3) == 15.0_wp, "golden: hc(3)")
         if (allocated(error)) exit checks
         call check(error, abs(hc(4) - (8.0_wp + H_VANISHED*0.5_wp)) < 1.0e-15_wp, &
                    "golden: hc(4)")
         if (allocated(error)) exit checks
         call check(error, hc(5) == 6.0_wp .and. hc(6) == 9.0_wp, "golden: hc(5,6)")
         if (allocated(error)) exit checks

         ! kf: zero everywhere except K=3 and K=6
         call check(error, kf(1) == 0.0_wp .and. kf(2) == 0.0_wp, "golden: kf(1,2)")
         if (allocated(error)) exit checks
         call check(error, abs(kf(3) - EXP_KF3) < 1.0e-14_wp, "golden: kf(3)")
         if (allocated(error)) exit checks
         call check(error, kf(4) == 0.0_wp .and. kf(5) == 0.0_wp, "golden: kf(4,5)")
         if (allocated(error)) exit checks
         call check(error, abs(kf(6) - EXP_KF6) < 1.0e-14_wp, "golden: kf(6)")
         if (allocated(error)) exit checks
         call check(error, kf(7) == 0.0_wp .and. kf(8) == 0.0_wp .and. &
                    kf(9) == 0.0_wp, "golden: kf(7,8,9)")
         if (allocated(error)) exit checks

         ! merged fields: surface/bed pure-massive layers recover exactly
         call massless_merge_fields(h, kc, NZ, nzc, u, v, t, s, uc, vc, tc, sc)
         call check(error, tc(1) == 22.0_wp .and. tc(3) == 18.0_wp .and. &
                    tc(6) == 12.0_wp, "golden: merged T pure layers")
         if (allocated(error)) exit checks
         ! merged layer 2 = [10, 7.5e-5] cluster, thickness-weighted T mean
         call check(error, abs(tc(2) - 19.999996250028122_wp) < 1.0e-12_wp, &
                    "golden: merged T(2)")
      end block checks
   end subroutine test_golden_maps

   ! ------------------------------------------------------------------
   ! T2 — identity invariant I1: healthy column => identity maps + bitwise
   ! interp_back roundtrip.
   ! ------------------------------------------------------------------
   subroutine test_identity(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 12
      real(wp) :: h(NZL), u(NZL), v(NZL), t(NZL), s(NZL)
      real(wp) :: hc(NZL), uc(NZL), vc(NZL), tc(NZL), sc(NZL)
      real(wp) :: kf(NZL + 1)
      integer :: kc(NZL + 1)
      real(wp) :: qc(NZLI), q(NZLI)
      integer :: nzc, k
      logical :: ok
      checks: block
         h = 0.0_wp; u = 0.0_wp; v = 0.0_wp; t = 0.0_wp; s = 0.0_wp
         do k = 1, NZ
            h(k) = 5.0_wp + real(k, wp)        ! all >> H_VANISHED
            u(k) = 0.01_wp*real(k, wp)
            v(k) = -0.02_wp*real(k, wp)
            t(k) = 20.0_wp - 0.5_wp*real(k, wp)
            s(k) = 35.0_wp + 0.01_wp*real(k, wp)
         end do

         call massless_build_maps(h, NZ, H_VANISHED, nzc, hc, kc, kf)

         call check(error, nzc == NZ, "I1: nzc must equal nz")
         if (allocated(error)) exit checks
         ok = .true.
         do k = 1, NZ
            if (kc(k) /= k) ok = .false.
            if (kf(k) /= 0.0_wp) ok = .false.
            if (hc(k) /= h(k)) ok = .false.   ! BITWISE
         end do
         if (kc(NZ + 1) /= NZ + 1) ok = .false.
         if (kf(NZ + 1) /= 0.0_wp) ok = .false.
         call check(error, ok, "I1: identity maps (kc==k, kf==0, hc==h bitwise)")
         if (allocated(error)) exit checks

         ! merged means: single-layer thickness-weighted mean = a*h/h.
         ! IEEE (a*h)/h is not guaranteed bit-identical to a, but a/h*h
         ! variants... we assert to 1 ULP.
         call massless_merge_fields(h, kc, NZ, nzc, u, v, t, s, uc, vc, tc, sc)
         ok = .true.
         do k = 1, NZ
            if (abs(tc(k) - t(k)) > 1.0e-14_wp*max(abs(t(k)), 1.0_wp)) ok = .false.
            if (abs(uc(k) - u(k)) > 1.0e-14_wp*max(abs(u(k)), 1.0_wp)) ok = .false.
         end do
         call check(error, ok, "I1: merged means recover layer values to 1 ULP")
         if (allocated(error)) exit checks

         ! interp_back on a random interface profile: BITWISE identity.
         do k = 1, NZ + 1
            qc(k) = 1.0e-3_wp*real(k*k - 3*k, wp)
         end do
         call massless_interp_back(qc, kc, kf, NZ, q)
         ok = .true.
         do k = 1, NZ + 1
            if (q(k) /= qc(k)) ok = .false.   ! BITWISE
         end do
         call check(error, ok, "I1: interp_back is bitwise identity on a healthy column")
      end block checks
   end subroutine test_identity

   ! ------------------------------------------------------------------
   ! T3 — conservation invariant I2: scattered-massless column conserves
   ! mass + thickness-weighted tracer integrals to round-off.
   ! ------------------------------------------------------------------
   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 11
      real(wp) :: h(NZL), u(NZL), v(NZL), t(NZL), s(NZL)
      real(wp) :: hc(NZL), uc(NZL), vc(NZL), tc(NZL), sc(NZL)
      real(wp) :: kf(NZL + 1)
      integer :: kc(NZL + 1)
      integer :: nzc, k
      real(wp) :: sum_h, sum_hc, sum_th, sum_thc, sum_sh, sum_shc
      checks: block
         h = 0.0_wp; u = 0.0_wp; v = 0.0_wp; t = 0.0_wp; s = 0.0_wp
         ! scattered massless at k=2,5,6,9 (1-based surface-down)
         h(1:NZ) = [10.0_wp, H_VANISHED*0.4_wp, 12.0_wp, 8.0_wp, &
                    H_VANISHED*0.6_wp, H_VANISHED*0.2_wp, 7.0_wp, 9.0_wp, &
                    H_VANISHED*0.3_wp, 5.0_wp, 6.0_wp]
         do k = 1, NZ
            t(k) = 18.0_wp - 0.3_wp*real(k, wp)
            s(k) = 34.5_wp + 0.05_wp*real(k, wp)
         end do

         call massless_build_maps(h, NZ, H_VANISHED, nzc, hc, kc, kf)
         call massless_merge_fields(h, kc, NZ, nzc, u, v, t, s, uc, vc, tc, sc)

         ! mass conservation: sum(hc) == sum(h)  (same accumulation order)
         sum_h = 0.0_wp; sum_th = 0.0_wp; sum_sh = 0.0_wp
         do k = 1, NZ
            sum_h = sum_h + h(k)
            sum_th = sum_th + t(k)*h(k)
            sum_sh = sum_sh + s(k)*h(k)
         end do
         sum_hc = 0.0_wp; sum_thc = 0.0_wp; sum_shc = 0.0_wp
         do k = 1, nzc
            sum_hc = sum_hc + hc(k)
            sum_thc = sum_thc + tc(k)*hc(k)
            sum_shc = sum_shc + sc(k)*hc(k)
         end do

         call check(error, abs(sum_hc - sum_h) <= 1.0e-13_wp*sum_h, &
                    "I2: thickness conserved")
         if (allocated(error)) exit checks
         call check(error, abs(sum_thc - sum_th) <= 1.0e-13_wp*abs(sum_th), &
                    "I2: T integral conserved")
         if (allocated(error)) exit checks
         call check(error, abs(sum_shc - sum_sh) <= 1.0e-13_wp*abs(sum_sh), &
                    "I2: S integral conserved")
      end block checks
   end subroutine test_conservation

   ! ------------------------------------------------------------------
   ! Kappa-shear end-to-end helpers
   ! ------------------------------------------------------------------
   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine setup_ks(ks, grid, nz, merge_on)
      type(ocean_kappa_shear_t), intent(inout) :: ks
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      logical, intent(in) :: merge_on
      call ks%init(grid, nz_ml=nz)
      ks%enable = .true.
      ks%massless_merge = merge_on
      ks%eos%variant = EOS_VARIANT_LINEAR
      ks%eos%alpha_T = ALPHA_LIN_KS
      ks%eos%beta_S = BETA_LIN_KS
      ks%eos%rho0 = RHO0_KS
      ks%rho0 = RHO0_KS
      call ks%set_f_centre(grid, F_CORIOLIS, 0.0_wp, 0.0_wp)
   end subroutine setup_ks

   !! ZSTAR_FULL-like column: bed-side layers vanished, sheared interior.
   !! has_vanish = .false. gives a healthy column (uniform dz) for the
   !! identity end-to-end test.
   subroutine setup_zstar_like_column(ms, nz, has_vanish)
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      logical, intent(in) :: has_vanish
      integer :: k, k_sd, n_vanish
      real(wp) :: dz_surf, z_sd, t_k, u_k, h_k
      real(wp), parameter :: DU = 0.5_wp, SHEAR_DEPTH = 50.0_wp
      real(wp), parameter :: SHEAR_WIDTH = 20.0_wp, STRAT_N2 = 3.0e-5_wp
      real(wp) :: dt_per_m

      dt_per_m = STRAT_N2*RHO0_KS/(GRAVITY*ALPHA_LIN_KS)
      dz_surf = 10.0_wp
      n_vanish = 0
      if (has_vanish) n_vanish = 3   ! 3 bed-side vanished layers

      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      do k = 1, nz
         ! local surface-down index k_sd = nz+1-k (k_sd=1 surface)
         k_sd = nz + 1 - k
         ! bed-side (largest k_sd) layers vanish
         if (k_sd > nz - n_vanish) then
            h_k = H_VANISHED*0.5_wp
         else
            h_k = dz_surf
         end if
         ! approximate depth-below-surface of the layer centre
         z_sd = (real(k_sd - 1, wp) + 0.5_wp)*dz_surf
         t_k = 20.0_wp - dt_per_m*z_sd
         u_k = DU*tanh((SHEAR_DEPTH - z_sd)/SHEAR_WIDTH)

         ms%h_layer(:, :, k) = h_k
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*h_k
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*h_k
         ms%u_face_x_layer(:, :, k) = u_k
      end do
   end subroutine setup_zstar_like_column

   ! ------------------------------------------------------------------
   ! T4a — kappa-shear end-to-end: vanished column, merge ON vs OFF differ.
   ! ------------------------------------------------------------------
   subroutine test_kshear_merge(error)
      !! On a column WITH bed-side vanished layers, massless_merge on vs
      !! off produce DIFFERENT kd_int (the merge engages, eliminating the
      !! blunt 1/H_VANISHED grid spacing at the vanished bed interfaces).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_on, ms_off
      type(ocean_kappa_shear_t) :: ks_on, ks_off
      integer, parameter :: NZ = 12
      real(wp), parameter :: DT = 1800.0_wp
      integer :: ip, jp, k
      real(wp) :: max_diff
      logical :: differ
      checks: block
         call make_grid(grid, 4, 4)
         ms_on%nz_ml = NZ; ms_off%nz_ml = NZ
         call ms_on%init(grid); call ms_off%init(grid)
         call setup_ks(ks_on, grid, NZ, .true.)
         call setup_ks(ks_off, grid, NZ, .false.)
         call setup_zstar_like_column(ms_on, NZ, .true.)
         call setup_zstar_like_column(ms_off, NZ, .true.)

         call kappa_shear_compute(grid, ks_on, ms_on, DT)
         call kappa_shear_compute(grid, ks_off, ms_off, DT)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         ! Both still pin bed + surface to 0.
         call check(error, ks_on%kd_int(ip, jp, 1) == 0.0_wp .and. &
                    ks_on%kd_int(ip, jp, NZ + 1) == 0.0_wp, &
                    "merge-on: bed/surface kd_int still 0")
         if (allocated(error)) exit checks

         ! Merge engages: profiles differ somewhere.
         max_diff = 0.0_wp
         differ = .false.
         do k = 1, NZ + 1
            max_diff = max(max_diff, abs(ks_on%kd_int(ip, jp, k) - &
                                         ks_off%kd_int(ip, jp, k)))
         end do
         differ = max_diff > 1.0e-9_wp
         call check(error, differ, &
                    "merge-on vs off must differ on a vanished column")
      end block checks
      call ks_on%destroy(); call ks_off%destroy()
      call ms_on%destroy(); call ms_off%destroy()
   end subroutine test_kshear_merge

   ! ------------------------------------------------------------------
   ! T4b — kappa-shear end-to-end identity (I1 through the whole kernel):
   ! on a HEALTHY column, merge ON == merge OFF BITWISE.
   ! ------------------------------------------------------------------
   subroutine test_kshear_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_on, ms_off
      type(ocean_kappa_shear_t) :: ks_on, ks_off
      integer, parameter :: NZ = 12
      real(wp), parameter :: DT = 1800.0_wp
      integer :: ip, jp, k
      logical :: bitwise
      checks: block
         call make_grid(grid, 4, 4)
         ms_on%nz_ml = NZ; ms_off%nz_ml = NZ
         call ms_on%init(grid); call ms_off%init(grid)
         call setup_ks(ks_on, grid, NZ, .true.)
         call setup_ks(ks_off, grid, NZ, .false.)
         call setup_zstar_like_column(ms_on, NZ, .false.)   ! healthy
         call setup_zstar_like_column(ms_off, NZ, .false.)

         call kappa_shear_compute(grid, ks_on, ms_on, DT)
         call kappa_shear_compute(grid, ks_off, ms_off, DT)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         bitwise = .true.
         do k = 1, NZ + 1
            if (ks_on%kd_int(ip, jp, k) /= ks_off%kd_int(ip, jp, k)) &
               bitwise = .false.
            if (ks_on%tke_int(ip, jp, k) /= ks_off%tke_int(ip, jp, k)) &
               bitwise = .false.
         end do
         call check(error, bitwise, &
                    "I1 end-to-end: merge on == off bitwise on a healthy column")
      end block checks
      call ks_on%destroy(); call ks_off%destroy()
      call ms_on%destroy(); call ms_off%destroy()
   end subroutine test_kshear_identity

end module test_ocean_massless
