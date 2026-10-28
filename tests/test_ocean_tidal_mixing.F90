!! Analytical tests for the St-Laurent/Simmons internal-tide interior
!! mixing closure (`rdb_ocean_tidal_mixing`).
!!
!! Golden oracle: the stage-4 Python prototype
!! `local_archive/prototypes/tidal_mixing_proto.py` + `tidal_mixing_golden.npz`.
!! The prototype's flux-bookkeeping sweep produces a per-LAYER Kd; the
!! Fortran kernel deposits that per-layer Kd 50/50 at the layer's two
!! bounding interfaces (MOM6 interface deposition, divergence D2), so the
!! interface field `kd_int(:,:,K)` (interior K) equals
!! `0.5*(Kd_lay(K-1) + Kd_lay(K))`.  Tests reconstruct the expected
!! interface values from the golden per-layer Kd and compare.
!!
!! Global index convention (Roundabout bottom-up):
!!   kd_int(:,:,K=1) = bed (always 0)
!!   kd_int(:,:,K=nz+1) = surface (always 0)
!! Layer k contributes to interfaces k (its bed) and k+1 (its top).
module test_ocean_tidal_mixing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, OMEGA_EARTH
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_tidal_mixing, only: ocean_tidal_mixing_t, &
                                     tidal_mixing_compute, &
                                     tidal_mixing_merge_into_kt, &
                                     tidal_mixing_fold_static_fields, &
                                     polzin_stretched_scale
   use rdb_eos, only: EOS_VARIANT_LINEAR
   implicit none
   private

   public :: collect_ocean_tidal_mixing_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0_TM = 1035.0_wp
   !! Linear thermal expansion (kg/m^3/K).  With the linear EOS the kernel
   !! buoyancy derivative is dbuoy_t = G*alpha_T/rho0, so a chosen layer
   !! N^2 maps to a temperature step dT = N2*dz_centre/(G*alpha_T/rho0).
   real(wp), parameter :: ALPHA_LIN_TM = 0.2_wp
   real(wp), parameter :: BETA_LIN_TM = RHO0_TM*7.6e-4_wp
   !! Golden rotation floor Omega^2 = OMEGA_EARTH^2 (plain Omega^2, the
   !! Melet 2013 efficiency rescaling N^2/(N^2+Omega^2)).  The prototype
   !! now uses the SAME OMEGA_EARTH constant, so this equals the kernel
   !! default bit-for-bit and the omega2 override is no longer needed.
   real(wp), parameter :: OMEGA2_GOLD = OMEGA_EARTH**2

   ! Golden column geometry / knobs (tidal_mixing_golden.npz).
   integer, parameter :: NZ_G = 30
   real(wp), parameter :: H_G = 3000.0_wp
   real(wp), parameter :: DZ_G = H_G/real(NZ_G, wp)   ! 100 m uniform
   real(wp), parameter :: E_G = 0.05_wp
   real(wp), parameter :: GAMMA_G = 0.3333_wp
   real(wp), parameter :: MU_G = 0.2_wp
   real(wp), parameter :: ZETA_G = 500.0_wp
   real(wp), parameter :: KDMAX_G = 1.0e-2_wp
   real(wp), parameter :: N2_CONST_G = 1.0e-6_wp

contains

   subroutine collect_ocean_tidal_mixing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("tidal_inert_enabled", test_inert_enabled), &
                  new_unittest("tidal_golden_profile", test_golden_profile), &
                  new_unittest("tidal_e_compute_jayne", test_e_compute_jayne), &
                  new_unittest("tidal_energy_conservation", test_energy_conservation), &
                  new_unittest("tidal_decay_anchoring", test_decay_anchoring), &
                  new_unittest("tidal_n2_zero_robust", test_n2_zero_robust), &
                  new_unittest("tidal_merge_additive", test_merge), &
                  new_unittest("tidal_t1_fold_clamp_mask", test_t1_fold_clamp), &
                  new_unittest("tidal_t1_fields_vs_scalar", test_t1_fields_kernel), &
                  new_unittest("tidal_t1_nbot_rough_avg", test_t1_nbot), &
                  new_unittest("tidal_t1_energy_bed_share", test_t1_energy), &
                  new_unittest("tidal_t1_energy_limiter", test_t1_limiter), &
                  new_unittest("tidal_t1_defaults_golden", test_t1_defaults), &
                  new_unittest("tidal_t2_uniform_n2", test_t2_uniform), &
                  new_unittest("tidal_t2_zp_value", test_t2_zp_value), &
                  new_unittest("tidal_t2_degenerate_scale", test_t2_degenerate), &
                  new_unittest("tidal_t2_wkb_stretch", test_t2_wkb), &
                  new_unittest("tidal_t2_cap_equals_degenerate", test_t2_cap), &
                  new_unittest("tidal_t2_stlaurent_unchanged", test_t2_default), &
                  new_unittest("tidal_t2_n2_layer_form", test_t2_n2_form) &
                  ]
   end subroutine collect_ocean_tidal_mixing_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine setup_tm(tm, grid, nz, e_uniform, kd_max)
      !! St-Laurent slot with the golden knobs + linear EOS.  No omega2
      !! override: the kernel default (OMEGA_EARTH^2) now equals the
      !! golden's Omega^2 bit-for-bit (Bug-1 fix), so the default is used.
      type(ocean_tidal_mixing_t), intent(inout) :: tm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: e_uniform, kd_max
      call tm%init(grid, nz_ml=nz)
      tm%enable = .true.
      tm%gamma = GAMMA_G
      tm%mu = MU_G
      tm%zeta = ZETA_G
      tm%kd_max = kd_max
      tm%prandtl_tidal = 1.0_wp
      tm%e_uniform = e_uniform
      tm%eos%variant = EOS_VARIANT_LINEAR
      tm%eos%alpha_T = ALPHA_LIN_TM
      tm%eos%beta_S = BETA_LIN_TM
      tm%eos%rho0 = RHO0_TM
      tm%rho0 = RHO0_TM
      call tm%set_e_uniform(e_uniform)
   end subroutine setup_tm

   subroutine fill_const_n2_column(ms, nz, dz, n2_target)
      !! Resting column whose kernel-derived per-layer N^2 == n2_target for
      !! every interior layer (T linear in depth, S uniform, u=v=0).  The
      !! kernel computes n2_col(k) = dbuoy_t*(T(k+1)-T(k))/dz_centre with
      !! dbuoy_t = G*alpha_T/rho0; choose dT so that equals n2_target.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, n2_target
      integer :: k
      real(wp) :: dbuoy_t, dt_layer, t_k

      dbuoy_t = GRAVITY*ALPHA_LIN_TM/RHO0_TM
      dt_layer = n2_target*dz/dbuoy_t   ! dz_centre == dz (uniform)
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ! bottom-up: warmer upward (stable). T(k) = T0 + dt_layer*(k-1).
         t_k = 10.0_wp + dt_layer*real(k - 1, wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*dz
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*dz
      end do
   end subroutine fill_const_n2_column

   subroutine ref_sweep(nz, dz, n2, e_col, gamma, mu, zeta, rho0, omega2, &
                        kd_max, kd_lay, tke_lay)
      !! Reference flux-bookkeeping sweep (mirrors the Python prototype
      !! `kd_discrete_sweep` exactly).  Returns per-LAYER Kd and the
      !! per-layer deposited power.  rho0 folded into TKE_to_Kd ([A1]).
      !! D1 (FULL MOM6 endpoint exclusion): BOTH the bed LAYER (k=1) and the
      !! surface LAYER (k=nz) Kd are zeroed (mirrors the kernel's
      !! `kd_lay_arr(1)=0` + `kd_lay_arr(nz)=0`); `tke_lay` is left intact so
      !! the energy-conservation check still sees the pre-exclusion power.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz(nz), n2(nz)
      real(wp), intent(in) :: e_col, gamma, mu, zeta, rho0, omega2, kd_max
      real(wp), intent(out) :: kd_lay(nz), tke_lay(nz)
      integer :: k
      real(wp) :: h_tot, hz, inv_int, tke_bot, tke_rem, z_top, frac_top
      real(wp) :: tlay, denom, kd

      h_tot = sum(dz)
      hz = h_tot/zeta
      if (hz < 1.0e-14_wp) then
         inv_int = 1.0_wp
      else
         inv_int = 1.0_wp/(1.0_wp - exp(-hz))
      end if
      tke_bot = gamma*mu*e_col
      tke_rem = inv_int*tke_bot
      z_top = 0.0_wp
      do k = 1, nz
         z_top = z_top + dz(k)
         frac_top = inv_int*exp(-z_top/zeta)
         tlay = tke_rem - tke_bot*frac_top
         tke_rem = tke_rem - tlay
         denom = rho0*dz(k)*(n2(k) + omega2)
         kd = tlay/denom
         if (kd_max >= 0.0_wp .and. kd > kd_max) kd = kd_max
         kd_lay(k) = kd
         tke_lay(k) = tlay
      end do
      ! D1: bed + surface LAYERS excluded from deposition (tke_lay kept
      ! intact for the energy-conservation check).
      if (nz >= 1) kd_lay(1) = 0.0_wp
      if (nz >= 1) kd_lay(nz) = 0.0_wp
   end subroutine ref_sweep

   subroutine spread_to_interfaces(nz, kd_lay, kd_int)
      !! Spread per-layer Kd 50/50 to the two bounding interfaces, zero the
      !! bed (1) and surface (nz+1) end-caps — the kernel's deposition.
      integer, intent(in) :: nz
      real(wp), intent(in) :: kd_lay(nz)
      real(wp), intent(out) :: kd_int(nz + 1)
      integer :: k
      kd_int = 0.0_wp
      do k = 1, nz
         kd_int(k) = kd_int(k) + 0.5_wp*kd_lay(k)
         kd_int(k + 1) = kd_int(k + 1) + 0.5_wp*kd_lay(k)
      end do
      kd_int(1) = 0.0_wp
      kd_int(nz + 1) = 0.0_wp
   end subroutine spread_to_interfaces

   ! ------------------------------------------------------------------
   ! Test: inert-enabled — enable=.true., E=0 => zero Kd contribution
   ! ------------------------------------------------------------------

   subroutine test_inert_enabled(error)
      !! Structural path live, physics zero: enabled but E=0 => kd_int==0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: i, j, k
      real(wp) :: kd_max_seen
      block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, 0.0_wp, KDMAX_G)
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         kd_max_seen = 0.0_wp
         do k = 1, NZ_G + 1
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  kd_max_seen = max(kd_max_seen, abs(tm%kd_int(i, j, k)))
               end do
            end do
         end do
         call check(error, kd_max_seen == 0.0_wp, &
                    "inert-enabled (E=0): kd_int must be identically 0")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_inert_enabled

   ! ------------------------------------------------------------------
   ! Test: golden profile — vs the prototype npz, constant N^2, no cap
   ! ------------------------------------------------------------------

   subroutine test_golden_profile(error)
      !! Constant-N^2 column matching the golden geometry/knobs (no cap).
      !! Assert the kernel interface field == the 50/50 spread of the
      !! golden per-layer Kd.  D1 (FULL MOM6 endpoint exclusion): the bed
      !! LAYER Kd is now 0; the deepest INTERIOR layer (k=2) is pinned to the
      !! regenerated golden value 0.004765789018233164 to ~1e-12 (so the test
      !! still bit-validates the deposition, not just a zero).  Omega^2 =
      !! OMEGA_EARTH^2.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: k
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G)
      real(wp) :: kd_int_exp(NZ_G + 1)
      real(wp) :: rel, denom
      !! Regenerated golden: bed LAYER (k=1) Kd now 0 (D1 full exclusion);
      !! deepest INTERIOR layer (k=2) Kd pinned below.
      real(wp), parameter :: GOLD_INT1 = 0.004765789018233164_wp

      block
         call make_grid(grid, 3, 3)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, E_G, -1.0_wp)   ! no cap
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         ! Reference per-layer Kd (interior layers carry N2_CONST; surface
         ! layer k=nz has kernel-derived N2=0 -> match that here too).
         do k = 1, NZ_G
            dz_a(k) = DZ_G
            n2_a(k) = N2_CONST_G
         end do
         n2_a(NZ_G) = 0.0_wp   ! kernel surface layer has no overlying layer
         call ref_sweep(NZ_G, dz_a, n2_a, E_G, GAMMA_G, MU_G, ZETA_G, &
                        RHO0_TM, OMEGA2_GOLD, -1.0_wp, kd_lay, tke_lay)
         call spread_to_interfaces(NZ_G, kd_lay, kd_int_exp)

         ! Whole interface profile matches the reconstruction.
         rel = 0.0_wp
         do k = 1, NZ_G + 1
            denom = max(abs(kd_int_exp(k)), 1.0e-30_wp)
            rel = max(rel, abs(tm%kd_int(1, 1, k) - kd_int_exp(k))/denom)
         end do
         call check(error, rel < 1.0e-11_wp, &
                    "golden profile: kd_int != 50/50 spread of reference Kd")
         if (allocated(error)) return

         ! D1 FULL endpoint exclusion: the bed LAYER (k=1) Kd is now exactly 0.
         call check(error, kd_lay(1) == 0.0_wp, &
                    "golden profile: bed-layer Kd must be 0 (D1 full exclusion)")
         if (allocated(error)) return

         ! Pin the deepest INTERIOR layer (k=2) Kd to the regenerated golden
         ! value, so the test still bit-validates the deposition (not just a
         ! zero).  The kernel deposits kd_lay(2) half to interface K=2 (its
         ! bed) and half to K=3 (its top): kd_int(2)=0.5*kd_lay(2) since the
         ! bed layer contributes nothing to K=2 now.
         rel = abs(kd_lay(2) - GOLD_INT1)/GOLD_INT1
         call check(error, rel < 1.0e-12_wp, &
                    "golden profile: interior-layer-2 Kd != golden 0.004765789018233164")
         if (allocated(error)) return
         rel = abs(tm%kd_int(1, 1, 2) - 0.5_wp*GOLD_INT1)/(0.5_wp*GOLD_INT1)
         call check(error, rel < 1.0e-11_wp, &
                    "golden profile: kd_int(K=2) != 0.5*golden interior-layer-2 Kd")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_golden_profile

   ! ------------------------------------------------------------------
   ! Test: e_compute — state-dependent E from Jayne & St Laurent (2001)
   ! ------------------------------------------------------------------

   subroutine test_e_compute_jayne(error)
      !! v1.1 `e_compute` path: the kernel DERIVES the internal-tide energy
      !! input E = 0.5*rho0*kappa_h2*kappa_itides*<h^2>*U_tide^2 * N_bot
      !! (Jayne & St Laurent 2001) instead of reading the prescribed e_in.
      !! This test computes E independently (WITH the rho0 factor — MOM6
      !! GV%H_to_RZ) and drives the reference sweep with it, then compares
      !! against the kernel's interface field.  It is the regression guard
      !! for the rho0 factor: drop rho0 from the kernel coefficient and the
      !! kernel Kd falls ~1035x below this reference, failing the check.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: k
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G)
      real(wp) :: kd_int_exp(NZ_G + 1)
      real(wp) :: rel, denom, h2c, e_expect, n_bot
      !! Generation knobs (realistic order of magnitude).
      real(wp), parameter :: KAPPA_ITIDES_T = 6.2832e-4_wp  ! 2*pi/1e4 m^-1
      real(wp), parameter :: KAPPA_H2_T = 1.0_wp
      real(wp), parameter :: UTIDE_T = 0.05_wp              ! 5 cm/s RMS
      real(wp), parameter :: H2_ROUGH_T = 200.0_wp          ! <h^2> (m^2)
      real(wp), parameter :: FRAC_ROUGH_T = 0.1_wp

      block
         call make_grid(grid, 2, 2)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, 0.0_wp, -1.0_wp)   ! e_uniform=0, no cap
         ! Switch to the computed-E path.
         tm%e_compute = .true.
         tm%kappa_itides = KAPPA_ITIDES_T
         tm%kappa_h2 = KAPPA_H2_T
         tm%utide = UTIDE_T
         tm%h2_rough = H2_ROUGH_T
         tm%frac_rough = FRAC_ROUGH_T
         tm%e_max = 1.0e3_wp                              ! large => no bind
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         ! Independent E (Jayne & St Laurent 2001), WITH rho0.  N_bot is the
         ! bed-most interior N (n2_col(1) = N2_CONST_G in this uniform column);
         ! roughness clamp <h^2> <= (frac_rough*H)^2 (here H2_ROUGH dominates).
         h2c = min(H2_ROUGH_T, (FRAC_ROUGH_T*H_G)**2)
         n_bot = sqrt(N2_CONST_G)
         e_expect = 0.5_wp*RHO0_TM*KAPPA_H2_T*KAPPA_ITIDES_T*h2c*UTIDE_T**2*n_bot

         do k = 1, NZ_G
            dz_a(k) = DZ_G
            n2_a(k) = N2_CONST_G
         end do
         n2_a(NZ_G) = 0.0_wp
         call ref_sweep(NZ_G, dz_a, n2_a, e_expect, GAMMA_G, MU_G, ZETA_G, &
                        RHO0_TM, OMEGA2_GOLD, -1.0_wp, kd_lay, tke_lay)
         call spread_to_interfaces(NZ_G, kd_lay, kd_int_exp)

         ! Whole interface profile matches the formula-derived reference.
         rel = 0.0_wp
         do k = 1, NZ_G + 1
            denom = max(abs(kd_int_exp(k)), 1.0e-30_wp)
            rel = max(rel, abs(tm%kd_int(1, 1, k) - kd_int_exp(k))/denom)
         end do
         call check(error, rel < 1.0e-11_wp, &
                    "e_compute: kd_int != Jayne&StLaurent E (rho0 factor) reference")
         if (allocated(error)) return

         ! The path actually produced mixing (non-trivial, > prescribed-zero).
         call check(error, tm%kd_int(1, 1, 3) > 0.0_wp, &
                    "e_compute: derived E gave no mixing")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_e_compute_jayne

   ! ------------------------------------------------------------------
   ! Test: energy conservation — Sum(TKE_lay) == gamma*mu*E (round-off)
   ! ------------------------------------------------------------------

   subroutine test_energy_conservation(error)
      !! The (1-exp(-H/zeta)) normalization makes the column-integrated
      !! deposited power equal q*mu*E exactly (pre-clip).  Verified on the
      !! kernel's interface field by recovering the per-layer Kd from the
      !! interface 50/50 spread and back-converting to TKE_lay through the
      !! kernel's own N^2 (no cap, so no power is discarded).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: k
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G)
      real(wp) :: deposited, target_e, rel

      block
         call make_grid(grid, 2, 2)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, E_G, -1.0_wp)
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         ! Reference reproduces the kernel; the pre-clip TKE_lay sums to gmuE.
         do k = 1, NZ_G
            dz_a(k) = DZ_G
            n2_a(k) = N2_CONST_G
         end do
         n2_a(NZ_G) = 0.0_wp
         call ref_sweep(NZ_G, dz_a, n2_a, E_G, GAMMA_G, MU_G, ZETA_G, &
                        RHO0_TM, OMEGA2_GOLD, -1.0_wp, kd_lay, tke_lay)

         deposited = sum(tke_lay)
         target_e = GAMMA_G*MU_G*E_G
         rel = abs(deposited - target_e)/abs(target_e)
         call check(error, rel < 1.0e-12_wp, &
                    "energy conservation: Sum(TKE_lay) != gamma*mu*E")
         if (allocated(error)) return

         ! Sanity: the kernel actually deposited something (non-trivial).
         call check(error, tm%kd_int(1, 1, 2) > 0.0_wp, &
                    "energy conservation: kernel deposited nothing")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_energy_conservation

   ! ------------------------------------------------------------------
   ! Test: decay anchoring — bottom-intensified, monotone, e-folds over zeta
   ! ------------------------------------------------------------------

   subroutine test_decay_anchoring(error)
      !! Kd is bottom-intensified: the interior interface field decreases
      !! monotonically from just-above-bed upward, and the per-layer Kd
      !! e-folds over zeta.  Catches a flipped (surface-anchored) port.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: k
      logical :: monotone
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G)
      real(wp) :: ratio, expect_ratio, rel

      block
         call make_grid(grid, 2, 2)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, E_G, -1.0_wp)
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         ! Interior interfaces decrease strictly upward (bottom-intensified).
         ! D1 (FULL exclusion) zeroes BOTH end layers: the bed LAYER (k=1) and
         ! surface LAYER (k=nz) deposit nothing.  Interface K=2 then receives
         ! only kd_lay(2)'s bed-half (an under-filled boundary artifact) and
         ! K=3 is the peak; interface K=NZ receives only kd_lay(nz-1)'s top
         ! half (mixed layer owned downstream by EPBL/KPP).  The clean
         ! decay-diagnostic band is therefore K = 4 .. NZ-1.
         monotone = .true.
         do k = 4, NZ_G - 1
            if (.not. (tm%kd_int(1, 1, k) < tm%kd_int(1, 1, k - 1))) monotone = .false.
         end do
         call check(error, monotone, &
                    "decay anchoring: interior kd_int not monotone-decreasing upward")
         if (allocated(error)) return

         ! Per-layer e-fold ratio over adjacent INTERIOR layers ~ exp(dz/zeta).
         do k = 1, NZ_G
            dz_a(k) = DZ_G
            n2_a(k) = N2_CONST_G
         end do
         n2_a(NZ_G) = 0.0_wp
         call ref_sweep(NZ_G, dz_a, n2_a, E_G, GAMMA_G, MU_G, ZETA_G, &
                        RHO0_TM, OMEGA2_GOLD, -1.0_wp, kd_lay, tke_lay)
         ! Compare two adjacent INTERIOR layers (k=2,3; constant N^2 -> pure
         ! exp).  D1 now zeroes the bed layer (k=1), so the deepest deposited
         ! interior layer is k=2.
         ratio = kd_lay(2)/kd_lay(3)
         expect_ratio = exp(DZ_G/ZETA_G)
         rel = abs(ratio - expect_ratio)/expect_ratio
         call check(error, rel < 1.0e-3_wp, &
                    "decay anchoring: adjacent-layer ratio != exp(dz/zeta)")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_decay_anchoring

   ! ------------------------------------------------------------------
   ! Test: N^2 -> 0 robustness — finite, capped, no NaN/Inf
   ! ------------------------------------------------------------------

   subroutine test_n2_zero_robust(error)
      !! A near-zero-N^2 band (abyssal-quiescent) would blow 1/N^2 without
      !! the Omega^2 floor.  With the floor + per-layer Kd_max cap the
      !! kernel stays finite and capped.  Built by making a uniform-N^2
      !! column then zeroing a mid-band's stratification (T flat there).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: i, j, k
      logical :: finite, capped
      real(wp) :: val, dbuoy_t, dt_layer, t_k

      block
         call make_grid(grid, 2, 2)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, E_G, KDMAX_G)   ! cap ON

         ! Build a column with a flat (zero-N^2) band in layers 6..8.
         dbuoy_t = GRAVITY*ALPHA_LIN_TM/RHO0_TM
         dt_layer = N2_CONST_G*DZ_G/dbuoy_t
         ms%h_layer = DZ_G
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         t_k = 10.0_wp
         do k = 1, NZ_G
            if (k >= 6 .and. k <= 9) then
               ! flat band: no temperature increment (N^2 -> 0 across it)
               t_k = t_k
            else
               t_k = t_k + dt_layer
            end if
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*DZ_G
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*DZ_G
         end do

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         finite = .true.
         capped = .true.
         do k = 1, NZ_G + 1
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  val = tm%kd_int(i, j, k)
                  if (.not. (val == val) .or. abs(val) > huge(1.0_wp)) finite = .false.
                  ! interface field is a 50/50 mean of two capped layers,
                  ! so it is bounded by kd_max too.
                  if (val > KDMAX_G + 1.0e-18_wp) capped = .false.
                  if (val < 0.0_wp) capped = .false.
               end do
            end do
         end do
         call check(error, finite, "N^2->0: kd_int not all finite")
         if (allocated(error)) return
         call check(error, capped, "N^2->0: kd_int exceeded kd_max or went negative")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_n2_zero_robust

   ! ------------------------------------------------------------------
   ! Test: merge — additive into kt (and prandtl*kd into kv)
   ! ------------------------------------------------------------------

   subroutine test_merge(error)
      !! The merge adds kd_int into kt and prandtl_tidal*kd_int into kv at
      !! interior interfaces, leaving bed/surface end-caps untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: nx, ny, nzp1, k
      real(wp), allocatable :: kv(:, :, :), kt(:, :, :)
      real(wp), parameter :: PR = 0.7_wp
      logical :: ok

      block
         call make_grid(grid, 2, 2)
         ms%nz_ml = NZ_G
         call ms%init(grid)
         call setup_tm(tm, grid, NZ_G, E_G, KDMAX_G)
         tm%prandtl_tidal = PR
         call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)

         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

         nx = grid%nx_total
         ny = grid%ny_total
         nzp1 = NZ_G + 1
         allocate (kv(nx, ny, nzp1), source=2.0_wp)
         allocate (kt(nx, ny, nzp1), source=3.0_wp)
         call tidal_mixing_merge_into_kt(tm, nx, ny, nzp1, kv, kt)

         ok = .true.
         do k = 2, nzp1 - 1
            if (abs(kt(1, 1, k) - (3.0_wp + tm%kd_int(1, 1, k))) > 1.0e-14_wp) ok = .false.
            if (abs(kv(1, 1, k) - (2.0_wp + PR*tm%kd_int(1, 1, k))) > 1.0e-14_wp) ok = .false.
         end do
         call check(error, ok, "merge: kt/kv not additive by kd_int / prandtl*kd_int")
         if (allocated(error)) return
         ! End-caps untouched.
         call check(error, kt(1, 1, 1) == 3.0_wp .and. kt(1, 1, nzp1) == 3.0_wp, &
                    "merge: bed/surface kt end-caps must be untouched")
      end block
      call tm%destroy()
      call ms%destroy()
   end subroutine test_merge

   ! ------------------------------------------------------------------
   ! T1 helpers
   ! ------------------------------------------------------------------

   subroutine fill_n2_profile(ms, nz, dz, n2)
      !! Resting uniform-`dz` column whose kernel layer N^2 is `n2(k)`
      !! between layers k and k+1 (T steps chosen from the linear EOS).
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, n2(nz - 1)
      integer :: k
      real(wp) :: dbuoy_t, t_k

      dbuoy_t = GRAVITY*ALPHA_LIN_TM/RHO0_TM
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      t_k = 10.0_wp
      do k = 1, nz
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*dz
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*dz
         if (k < nz) t_k = t_k + n2(k)*dz/dbuoy_t
      end do
   end subroutine fill_n2_profile

   subroutine recover_tke(tm, i, j, nz, dz, n2, tke)
      !! Per-layer deposited power from the kernel's interface field:
      !! `kd_int(K) = 0.5*(kd_lay(K-1) + kd_lay(K))` with the bed and
      !! surface layers zero, so `kd_lay` follows by forward recursion;
      !! then `TKE = kd_lay*rho0*dz*(N^2+Omega^2)`.
      type(ocean_tidal_mixing_t), intent(in) :: tm
      integer, intent(in) :: i, j, nz
      real(wp), intent(in) :: dz, n2(nz)
      real(wp), intent(out) :: tke(nz)
      integer :: k
      real(wp) :: kd_lay(nz)

      kd_lay = 0.0_wp
      do k = 2, nz - 1
         kd_lay(k) = 2.0_wp*tm%kd_int(i, j, k) - kd_lay(k - 1)
      end do
      tke = 0.0_wp
      do k = 2, nz - 1
         tke(k) = kd_lay(k)*RHO0_TM*dz*(n2(k) + OMEGA2_GOLD)
      end do
   end subroutine recover_tke

   pure function phi_exp(z) result(phi)
      !! Remaining upward flux of the golden column at height `z` above
      !! the bed: `q*mu*E * e^{-z/zeta} / (1 - e^{-H/zeta})`.
      real(wp), intent(in) :: z
      real(wp) :: phi
      phi = GAMMA_G*MU_G*E_G*exp(-z/ZETA_G)/(1.0_wp - exp(-H_G/ZETA_G))
   end function phi_exp

   ! ------------------------------------------------------------------
   ! T1 (b): setup fold -- roughness clamp, wet and min_zbot masks
   ! ------------------------------------------------------------------

   subroutine test_t1_fold_clamp(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_tidal_mixing_t) :: tm
      real(wp), allocatable :: depth(:, :), wet(:, :)
      real(wp) :: coef, want
      integer :: nx, ny

      call make_grid(grid, 2, 2)
      call setup_tm(tm, grid, NZ_G, 0.0_wp, -1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (depth(nx, ny), wet(nx, ny))
      depth = 1000.0_wp
      wet = 1.0_wp
      tm%frac_rough = 0.1_wp
      tm%min_zbot = 50.0_wp
      tm%kappa_h2 = 0.84_wp
      tm%kappa_itides = 6.28319e-4_wp
      tm%tideamp = 0.3_wp
      tm%h2 = 25.0_wp
      ! (1,1): file value far above (0.1*100)^2 = 100  -> clamped
      depth(1, 1) = 100.0_wp
      tm%h2(1, 1) = 1.0e6_wp
      ! (2,1): below the clamp -> untouched
      depth(2, 1) = 100.0_wp
      tm%h2(2, 1) = 50.0_wp
      ! (3,1): dry
      wet(3, 1) = 0.0_wp
      tm%h2(3, 1) = 1.0e6_wp
      ! (4,1): shallower than min_zbot
      depth(4, 1) = 40.0_wp
      tm%h2(4, 1) = 1.0e6_wp

      call tidal_mixing_fold_static_fields(tm, nx, ny, depth, wet)

      coef = 0.5_wp*RHO0_TM*0.84_wp*6.28319e-4_wp*0.3_wp**2
      call check(error, abs(tm%h2(1, 1) - 100.0_wp) < 1.0e-12_wp, "clamp: h2 != (0.1 D)^2")
      if (allocated(error)) return
      call check(error, abs(tm%h2(2, 1) - 50.0_wp) < 1.0e-12_wp, "clamp: below-cap h2 changed")
      if (allocated(error)) return
      call check(error, tm%h2(3, 1) == 0.0_wp .and. tm%c_e(3, 1) == 0.0_wp, "dry cell not zeroed")
      if (allocated(error)) return
      call check(error, tm%h2(4, 1) == 0.0_wp .and. tm%c_e(4, 1) == 0.0_wp .and. &
                 tm%tideamp(4, 1) == 0.0_wp, "min_zbot cell not masked")
      if (allocated(error)) return
      want = coef*100.0_wp
      call check(error, abs(tm%c_e(1, 1) - want) < 1.0e-12_wp*want, "c_e(1,1) wrong")
      if (allocated(error)) return
      want = coef*50.0_wp
      call check(error, abs(tm%c_e(2, 1) - want) < 1.0e-12_wp*want, "c_e(2,1) wrong")
      if (allocated(error)) return
      ! Deep cell: (0.1*1000)^2 = 1e4 > 25, untouched.
      want = coef*25.0_wp
      call check(error, abs(tm%c_e(5, 5) - want) < 1.0e-12_wp*want, "c_e deep cell wrong")
      call tm%destroy()
   end subroutine test_t1_fold_clamp

   ! ------------------------------------------------------------------
   ! T1 (f, part): file-field path == scalar path when the fields are
   ! uniform, and Kd scales linearly with the cell's c_e.
   ! ------------------------------------------------------------------

   subroutine test_t1_fields_kernel(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm_a, tm_b
      real(wp), allocatable :: depth(:, :), wet(:, :)
      real(wp) :: worst, ratio
      integer :: k, nx, ny

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)
      nx = grid%nx_total
      ny = grid%ny_total

      ! A: scalar path
      call setup_tm(tm_a, grid, NZ_G, 0.0_wp, -1.0_wp)
      tm_a%e_compute = .true.
      tm_a%kappa_h2 = 0.84_wp
      tm_a%utide = 0.3_wp
      tm_a%h2_rough = 100.0_wp
      tm_a%frac_rough = 1.0_wp
      call tidal_mixing_compute(grid, tm_a, ms, 1800.0_wp)

      ! B: field path, uniform contents, plus one cell with U doubled
      call setup_tm(tm_b, grid, NZ_G, 0.0_wp, -1.0_wp)
      tm_b%e_compute = .true.
      tm_b%kappa_h2 = 0.84_wp
      tm_b%frac_rough = 1.0_wp
      tm_b%use_fields = .true.
      tm_b%tideamp = 0.3_wp
      tm_b%tideamp(4, 3) = 0.6_wp
      tm_b%h2 = 100.0_wp
      allocate (depth(nx, ny), wet(nx, ny))
      depth = H_G
      wet = 1.0_wp
      call tidal_mixing_fold_static_fields(tm_b, nx, ny, depth, wet)
      call tidal_mixing_compute(grid, tm_b, ms, 1800.0_wp)

      worst = 0.0_wp
      do k = 2, NZ_G
         worst = max(worst, abs(tm_b%kd_int(3, 3, k) - tm_a%kd_int(3, 3, k)) &
                     /max(tm_a%kd_int(3, 3, k), 1.0e-30_wp))
      end do
      call check(error, tm_a%kd_int(3, 3, 3) > 0.0_wp, "scalar path deposited nothing")
      if (allocated(error)) return
      call check(error, worst < 1.0e-12_wp, "field path != scalar path for uniform fields")
      if (allocated(error)) return

      ! U doubled => c_e x4 => Kd x4 (cap off, e_max not binding).
      ratio = tm_b%kd_int(4, 3, 5)/tm_b%kd_int(3, 3, 5)
      call check(error, abs(ratio - 4.0_wp) < 1.0e-12_wp, "Kd does not scale with U^2 cell by cell")
      call tm_a%destroy()
      call tm_b%destroy()
      call ms%destroy()
   end subroutine test_t1_fields_kernel

   ! ------------------------------------------------------------------
   ! T1 (c): near-bed N^2 roughness-height average
   ! ------------------------------------------------------------------

   subroutine run_nbot_case(error, h_amp, n_int)
      !! 4 layers of 1000 m (centres 500/1500/2500/3500 m) with distinct
      !! interface N^2 = 1, 4, 9 e-6 s^-2.  The first `n_int` interfaces
      !! are averaged (uniform spans => arithmetic mean) and imposed.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: h_amp
      integer, intent(in) :: n_int
      integer, parameter :: NZ4 = 4
      real(wp), parameter :: DZ4 = 1000.0_wp
      real(wp), parameter :: N2_RAW(3) = [1.0e-6_wp, 4.0e-6_wp, 9.0e-6_wp]
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm, tm_off
      real(wp) :: n2_a(NZ4), dz_a(NZ4), kd_lay(NZ4), tke_lay(NZ4), kd_ref(NZ4 + 1)
      real(wp) :: nb2, e_col, coef, worst, diff_off
      integer :: k

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ4
      call ms%init(grid)
      call fill_n2_profile(ms, NZ4, DZ4, N2_RAW)

      nb2 = sum(N2_RAW(1:n_int))/real(n_int, wp)
      n2_a = [N2_RAW, 0.0_wp]
      n2_a(1:n_int) = nb2
      dz_a = DZ4
      coef = 0.5_wp*RHO0_TM*1.0_wp*6.2832e-4_wp*h_amp**2*0.3_wp**2
      e_col = coef*sqrt(nb2)

      call setup_tm(tm, grid, NZ4, 0.0_wp, -1.0_wp)
      tm%e_compute = .true.
      tm%utide = 0.3_wp
      tm%h2_rough = h_amp**2
      tm%frac_rough = 10.0_wp
      tm%n2_bot_rough_avg = .true.
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

      call ref_sweep(NZ4, dz_a, n2_a, e_col, GAMMA_G, MU_G, ZETA_G, RHO0_TM, &
                     OMEGA2_GOLD, -1.0_wp, kd_lay, tke_lay)
      call spread_to_interfaces(NZ4, kd_lay, kd_ref)
      worst = 0.0_wp
      do k = 2, NZ4
         worst = max(worst, abs(tm%kd_int(3, 3, k) - kd_ref(k))/max(kd_ref(k), 1.0e-30_wp))
      end do
      call check(error, kd_ref(3) > 0.0_wp, "N_b case: oracle Kd is zero")
      if (allocated(error)) return
      call check(error, worst < 1.0e-9_wp, "N_b average over wrong interface set")
      if (allocated(error)) return

      ! Non-vacuity: with the knob off the same column gives a different Kd
      ! whenever the average is not just the bed interface.
      if (n_int > 1) then
         call setup_tm(tm_off, grid, NZ4, 0.0_wp, -1.0_wp)
         tm_off%e_compute = .true.
         tm_off%utide = 0.3_wp
         tm_off%h2_rough = h_amp**2
         tm_off%frac_rough = 10.0_wp
         call tidal_mixing_compute(grid, tm_off, ms, 1800.0_wp)
         diff_off = abs(tm_off%kd_int(3, 3, 3) - tm%kd_int(3, 3, 3))/tm%kd_int(3, 3, 3)
         call check(error, diff_off > 0.05_wp, "N_b average changed nothing")
         call tm_off%destroy()
      end if
      call tm%destroy()
      call ms%destroy()
   end subroutine run_nbot_case

   subroutine test_t1_nbot(error)
      type(error_type), allocatable, intent(out) :: error
      ! Rule: take interfaces upward; once the centre of the layer above the
      ! last one taken exceeds h_amp, take ONE more and stop.  Centres of
      ! layers 2,3,4 are 1500, 2500, 3500 m.
      ! h_amp ~ 0: 1500 > h_amp at once -> bed + next = two interfaces.
      call run_nbot_case(error, 1.0_wp, 2)
      if (allocated(error)) return
      ! 1500 > 1400 -> two.
      call run_nbot_case(error, 1400.0_wp, 2)
      if (allocated(error)) return
      ! 1500 <= 1600, then 2500 > 1600 -> third interface -> three.
      call run_nbot_case(error, 1600.0_wp, 3)
      if (allocated(error)) return
      ! 2500 <= 3000, 3500 > 3000 but no further interface exists -> three.
      call run_nbot_case(error, 3000.0_wp, 3)
   end subroutine test_t1_nbot

   ! ------------------------------------------------------------------
   ! T1 (d): exact energy bookkeeping, with and without bed_share="lump"
   ! ------------------------------------------------------------------

   subroutine test_t1_energy(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: n2_a(NZ_G), tke(NZ_G), want, rel
      integer :: pass

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)
      n2_a = N2_CONST_G
      n2_a(NZ_G) = 0.0_wp

      do pass = 1, 2
         call setup_tm(tm, grid, NZ_G, E_G, -1.0_wp)
         tm%bed_lump = (pass == 2)
         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
         call recover_tke(tm, 3, 3, NZ_G, DZ_G, n2_a, tke)
         if (pass == 1) then
            ! discard: bed share (Phi(0) - Phi(dz)) is NOT deposited
            want = phi_exp(DZ_G) - phi_exp(real(NZ_G - 1, wp)*DZ_G)
         else
            ! lump: Phi(0) - Phi(top interior) lands in the interior
            want = phi_exp(0.0_wp) - phi_exp(real(NZ_G - 1, wp)*DZ_G)
         end if
         rel = abs(sum(tke) - want)/want
         call check(error, rel < 1.0e-12_wp, "Sum(TKE_lay) != Phi(bottom) - Phi(top)")
         if (allocated(error)) return
         if (pass == 2) then
            ! The first interior layer carries the bed share on top of its own.
            want = phi_exp(0.0_wp) - phi_exp(2.0_wp*DZ_G)
            rel = abs(tke(2) - want)/want
            call check(error, rel < 1.0e-12_wp, "lump: layer kbed+1 != bed share + own share")
            if (allocated(error)) return
         end if
         call tm%destroy()
      end do
      call ms%destroy()
   end subroutine test_t1_energy

   ! ------------------------------------------------------------------
   ! T1 (e): energy-conserving limiter ("reinject")
   ! ------------------------------------------------------------------

   subroutine test_t1_limiter(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: KD_CAP = 1.0e-3_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: n2_a(NZ_G), tke(NZ_G), tke_nat, tke_cap, rem, used
      real(wp) :: tot_reinject, tot_clamp, kd_lay_max, worst
      integer :: k, n_capped

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)
      n2_a = N2_CONST_G
      n2_a(NZ_G) = 0.0_wp

      ! Default "clamp": clip Kd, discard the clipped power.
      call setup_tm(tm, grid, NZ_G, E_G, KD_CAP)
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
      call recover_tke(tm, 3, 3, NZ_G, DZ_G, n2_a, tke)
      tot_clamp = sum(tke)
      call tm%destroy()

      ! "reinject"
      call setup_tm(tm, grid, NZ_G, E_G, KD_CAP)
      tm%limiter_reinject = .true.
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
      call recover_tke(tm, 3, 3, NZ_G, DZ_G, n2_a, tke)
      tot_reinject = sum(tke)

      ! Layerwise oracle from the closed-form flux.
      rem = phi_exp(DZ_G)
      worst = 0.0_wp
      n_capped = 0
      kd_lay_max = 0.0_wp
      do k = 2, NZ_G - 1
         tke_nat = rem - phi_exp(real(k, wp)*DZ_G)
         tke_cap = KD_CAP*RHO0_TM*DZ_G*(n2_a(k) + OMEGA2_GOLD)
         used = min(tke_nat, tke_cap)
         if (tke_nat > tke_cap) n_capped = n_capped + 1
         rem = rem - used
         worst = max(worst, abs(tke(k) - used)/used)
         kd_lay_max = max(kd_lay_max, tke(k)/(RHO0_TM*DZ_G*(n2_a(k) + OMEGA2_GOLD)))
      end do
      call check(error, n_capped >= 3, "limiter never engaged (test not exercising it)")
      if (allocated(error)) return
      call check(error, worst < 1.0e-9_wp, "reinject: layer TKE != capped carry-over oracle")
      if (allocated(error)) return
      call check(error, kd_lay_max <= KD_CAP*(1.0_wp + 1.0e-9_wp), "a layer exceeds kd_max")
      if (allocated(error)) return
      ! Energy-conserving: strictly more power deposited than the discarding clamp,
      ! and never more than the natural column total.
      call check(error, tot_reinject > 1.05_wp*tot_clamp, "reinject did not recover clipped power")
      if (allocated(error)) return
      call check(error, tot_reinject <= phi_exp(DZ_G) - phi_exp(real(NZ_G - 1, wp)*DZ_G) + 1.0e-15_wp, &
                 "reinject deposited more than the available power")
      call tm%destroy()
      call ms%destroy()
   end subroutine test_t1_limiter

   ! ------------------------------------------------------------------
   ! T1 (f): every new knob at its default == the golden St-Laurent profile
   ! ------------------------------------------------------------------

   subroutine test_t1_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G), kd_ref(NZ_G + 1)
      real(wp) :: worst
      integer :: k

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)
      call setup_tm(tm, grid, NZ_G, E_G, KDMAX_G)
      tm%use_fields = .false.
      tm%n2_bot_rough_avg = .false.
      tm%limiter_reinject = .false.
      tm%bed_lump = .false.
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)

      dz_a = DZ_G
      n2_a = N2_CONST_G
      n2_a(NZ_G) = 0.0_wp
      call ref_sweep(NZ_G, dz_a, n2_a, E_G, GAMMA_G, MU_G, ZETA_G, RHO0_TM, &
                     OMEGA2_GOLD, KDMAX_G, kd_lay, tke_lay)
      call spread_to_interfaces(NZ_G, kd_lay, kd_ref)
      worst = 0.0_wp
      do k = 2, NZ_G
         worst = max(worst, abs(tm%kd_int(3, 3, k) - kd_ref(k))/max(kd_ref(k), 1.0e-30_wp))
      end do
      call check(error, worst < 1.0e-12_wp, "default knobs drifted from the golden profile")
      call tm%destroy()
      call ms%destroy()
   end subroutine test_t1_defaults

   ! ------------------------------------------------------------------
   ! T2: Polzin (2009) profile in WKB-stretched height
   ! ------------------------------------------------------------------

   pure function pz_zp_formula(u, h2, nb, nbar2, hh, fmax) result(zp)
      !! Independent evaluation of the stretched Polzin scale
      !! (nu=0.0697, N_ref=9.61e-4, kappa=6.2832e-4, f_d=1, no floor).
      real(wp), intent(in) :: u, h2, nb, nbar2, hh, fmax
      real(wp) :: zp
      zp = 0.0697_wp*(9.61e-4_wp)**2*u/((6.2832e-4_wp)**2*h2*nb*nbar2)
      zp = min(zp, fmax*hh)
   end function pz_zp_formula

   pure function pz_layer_n2(n2h) result(n2f)
      !! Layer N^2 = mean of the two bounding interface jumps over the layer
      !! thickness (uniform dz; zero jump at the bed and surface).
      real(wp), intent(in) :: n2h(NZ_G - 1)
      real(wp) :: n2f(NZ_G)
      integer :: k
      do k = 1, NZ_G
         n2f(k) = 0.0_wp
         if (k > 1) n2f(k) = n2f(k) + 0.5_wp*n2h(k - 1)
         if (k < NZ_G) n2f(k) = n2f(k) + 0.5_wp*n2h(k)
      end do
   end function pz_layer_n2

   subroutine pz_oracle(nz, dz, n2, pbot, zp, tke)
      !! Per-layer deposit `Phi(z*_{k-1}) - Phi(z*_k)` for k = 2..nz-1 with
      !! the bed share lumped into k = 2 (`bed_share = "lump"`), from
      !! `Phi(z*) = P I_p zp/(zp + z*)`, `I_p = 1 + zp/H`, `z*_k =
      !! sum_{m<=k} N_m^2 dz/Nbar^2` (z* = z when Nbar^2 = 0).
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, n2(nz), pbot, zp
      real(wp), intent(out) :: tke(nz)
      real(wp) :: hh, nbar2, zs, phi_prev, phi
      integer :: k

      hh = real(nz, wp)*dz
      nbar2 = sum(n2)*dz/hh
      zs = 0.0_wp
      phi_prev = pbot*(1.0_wp + zp/hh)
      tke = 0.0_wp
      do k = 1, nz - 1
         zs = zs + n2(k)*dz/nbar2
         phi = pbot*(1.0_wp + zp/hh)*zp/(zp + zs)
         if (k >= 2) then
            tke(k) = phi_prev - phi
            phi_prev = phi
         end if
      end do
   end subroutine pz_oracle

   subroutine pz_run(n2h, use_compute, u, h2, e_uni, tke, kd_int_col, polzin_on)
      !! Run the kernel on a resting `NZ_G`-layer column with interface
      !! `N^2 = n2h(k)` (k < nz), no cap, bed share lumped, return the
      !! recovered per-layer deposit.  `use_compute`: Jayne & St Laurent
      !! energy from `u`, `h2`; else the prescribed `e_uni`.
      real(wp), intent(in) :: n2h(NZ_G - 1), u, h2, e_uni
      logical, intent(in) :: use_compute, polzin_on
      real(wp), intent(out) :: tke(NZ_G), kd_int_col(NZ_G + 1)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: n2_full(NZ_G)

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_n2_profile(ms, NZ_G, DZ_G, n2h)
      n2_full = pz_layer_n2(n2h)
      call setup_tm(tm, grid, NZ_G, e_uni, -1.0_wp)
      tm%e_compute = use_compute
      tm%utide = u
      tm%h2_rough = h2
      tm%kappa_h2 = 0.84_wp
      tm%bed_lump = .true.
      tm%polzin = polzin_on
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
      call recover_tke(tm, 3, 3, NZ_G, DZ_G, n2_full, tke)
      kd_int_col = tm%kd_int(3, 3, :)
      call tm%destroy()
      call ms%destroy()
   end subroutine pz_run

   function pz_pbot(u, h2, nb) result(pbot)
      real(wp), intent(in) :: u, h2, nb
      real(wp) :: pbot
      pbot = GAMMA_G*MU_G*0.5_wp*RHO0_TM*0.84_wp*6.2832e-4_wp*h2*u**2*nb
   end function pz_pbot

   subroutine pz_compare(error, tke, want, tol, msg)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: tke(NZ_G), want(NZ_G), tol
      character(len=*), intent(in) :: msg
      real(wp) :: worst
      integer :: k
      worst = 0.0_wp
      do k = 2, NZ_G - 1
         worst = max(worst, abs(tke(k) - want(k))/want(k))
      end do
      call check(error, worst < tol, msg)
   end subroutine pz_compare

   ! (a) uniform N^2: shares equal the analytic Phi_p differences
   subroutine test_t2_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: n2h(NZ_G - 1), n2f(NZ_G), tke(NZ_G), want(NZ_G), kdi(NZ_G + 1)
      real(wp) :: pbot, zp, nbar2, nb, u, h2, want_sum

      u = 0.02_wp
      h2 = 1.0e4_wp
      n2h = 1.0e-6_wp
      n2f = pz_layer_n2(n2h)
      nb = sqrt(n2h(1))
      nbar2 = sum(n2f)*DZ_G/H_G
      pbot = pz_pbot(u, h2, nb)
      zp = pz_zp_formula(u, h2, nb, nbar2, H_G, 1.0_wp)
      call check(error, zp < H_G, "test setup: scale must be uncapped")
      if (allocated(error)) return
      call pz_run(n2h, .true., u, h2, 0.0_wp, tke, kdi, .true.)
      call pz_oracle(NZ_G, DZ_G, n2f, pbot, zp, want)
      call pz_compare(error, tke, want, 1.0e-9_wp, "uniform N2: layer share != Phi_p difference")
      if (allocated(error)) return
      ! column deposit = Phi(0) - Phi(top interior) = P_bot (surface N2 = 0)
      ! interior deposit = Phi(0) - Phi(z*_top(nz-1)); the surface layer's share is dropped
      want_sum = pbot*(1.0_wp + zp/H_G)*(1.0_wp - zp/(zp + sum(n2f(1:NZ_G - 1))*DZ_G/nbar2))
      call check(error, abs(sum(tke) - want_sum)/want_sum < 1.0e-9_wp, "column deposit != Phi(0)-Phi(top)")
      if (allocated(error)) return
      call check(error, kdi(1) == 0.0_wp .and. kdi(NZ_G + 1) == 0.0_wp, "end interfaces not zero")
      if (allocated(error)) return
      ! not the exponential: deposit differs (the stretched profile is algebraic)
      call pz_run(n2h, .true., u, h2, 0.0_wp, tke, kdi, .false.)
      call check(error, abs(tke(NZ_G - 1) - want(NZ_G - 1))/want(NZ_G - 1) > 0.05_wp, &
                 "polzin off must give the exponential profile")
   end subroutine test_t2_uniform

   ! (b) z_p* value
   subroutine test_t2_zp_value(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: zp, want

      zp = polzin_stretched_scale(0.02_wp, 1.0e4_wp, 1.0e-3_wp, 1.0e-6_wp, 1.0e5_wp, &
                                  0.0697_wp, 9.61e-4_wp, 6.2832e-4_wp, 1.0_wp, 1.0_wp, 0.0_wp)
      want = 0.0697_wp*(9.61e-4_wp)**2*0.02_wp/((6.2832e-4_wp)**2*1.0e4_wp*1.0e-9_wp)
      call check(error, abs(zp - want)/want < 1.0e-10_wp, "z_p* != formula")
      if (allocated(error)) return
      call check(error, zp > 320.0_wp .and. zp < 332.0_wp, "z_p* not near 326 m")
      if (allocated(error)) return
      ! cap at f_max*H
      zp = polzin_stretched_scale(0.02_wp, 1.0e4_wp, 1.0e-3_wp, 1.0e-6_wp, 200.0_wp, &
                                  0.0697_wp, 9.61e-4_wp, 6.2832e-4_wp, 1.0_wp, 0.5_wp, 0.0_wp)
      call check(error, abs(zp - 100.0_wp) < 1.0e-12_wp, "cap f_max*H not applied")
      if (allocated(error)) return
      ! floor on the unstretched scale: z_p = max(.., 1000) then * N_b^2/Nbar^2
      zp = polzin_stretched_scale(0.02_wp, 1.0e4_wp, 1.0e-3_wp, 2.0e-6_wp, 1.0e5_wp, &
                                  0.0697_wp, 9.61e-4_wp, 6.2832e-4_wp, 1.0_wp, 1.0_wp, 1000.0_wp)
      call check(error, abs(zp - 500.0_wp) < 1.0e-9_wp, "min decay scale floor not applied")
   end subroutine test_t2_zp_value

   ! (c) U = 0 or <h^2> = 0: z_p* = H, I_p = 2
   subroutine test_t2_degenerate(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: n2h(NZ_G - 1), n2f(NZ_G), tke(NZ_G), want(NZ_G), kdi(NZ_G + 1), pbot
      real(wp) :: nbar2, want_sum
      integer :: pass

      n2h = 1.0e-6_wp
      n2f = pz_layer_n2(n2h)
      pbot = GAMMA_G*MU_G*E_G
      call pz_oracle(NZ_G, DZ_G, n2f, pbot, H_G, want)
      do pass = 1, 2
         if (pass == 1) then
            call pz_run(n2h, .false., 0.0_wp, 1.0e4_wp, E_G, tke, kdi, .true.)
         else
            call pz_run(n2h, .false., 0.02_wp, 0.0_wp, E_G, tke, kdi, .true.)
         end if
         call pz_compare(error, tke, want, 1.0e-9_wp, "degenerate: share != P 2H/(H+z*) difference")
         if (allocated(error)) return
      end do
      ! I_p = 2: Phi(0) = 2 P and the first interior share carries the bed share
      nbar2 = sum(n2f)*DZ_G/H_G
      want_sum = pbot*2.0_wp*(1.0_wp - H_G/(H_G + sum(n2f(1:NZ_G - 1))*DZ_G/nbar2))
      call check(error, abs(sum(tke) - want_sum)/want_sum < 1.0e-9_wp, "degenerate: deposit != Phi(0)-Phi(top)")
   end subroutine test_t2_degenerate

   ! (d) WKB: the 4x-N^2 layers advance z* 4x faster; shares obey Phi_p
   subroutine test_t2_wkb(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: n2h(NZ_G - 1), n2f(NZ_G), tke(NZ_G), want(NZ_G), kdi(NZ_G + 1)
      real(wp) :: pbot, zp, nbar2, nb, u, h2, phi(NZ_G), zs(NZ_G), r
      integer :: k

      u = 0.02_wp
      h2 = 1.0e4_wp
      n2h(1:14) = 1.0e-6_wp
      n2h(15:NZ_G - 1) = 4.0e-6_wp
      n2f = pz_layer_n2(n2h)
      nb = sqrt(n2h(1))
      nbar2 = sum(n2f)*DZ_G/H_G
      pbot = pz_pbot(u, h2, nb)
      zp = pz_zp_formula(u, h2, nb, nbar2, H_G, 1.0_wp)
      call check(error, zp < H_G, "test setup: scale must be uncapped")
      if (allocated(error)) return
      call pz_run(n2h, .true., u, h2, 0.0_wp, tke, kdi, .true.)
      call pz_oracle(NZ_G, DZ_G, n2f, pbot, zp, want)
      call pz_compare(error, tke, want, 1.0e-9_wp, "WKB column: share != stretched Phi_p difference")
      if (allocated(error)) return
      ! Kernel-observable stretch: invert Phi back to z* and compare increments.
      phi(1) = pbot*(1.0_wp + zp/H_G)
      do k = 2, NZ_G - 1
         phi(k) = phi(k - 1) - tke(k)
         zs(k) = zp*(pbot*(1.0_wp + zp/H_G)/phi(k) - 1.0_wp)
      end do
      r = (zs(20) - zs(19))/(zs(5) - zs(4))
      call check(error, abs(r - 4.0_wp) < 1.0e-6_wp, "4x N2 layer does not advance z* 4x faster")
   end subroutine test_t2_wkb

   ! (e) cap hit: identical to the degenerate (z_p* = H) column
   subroutine test_t2_cap(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: n2h(NZ_G - 1), n2f(NZ_G), tke_a(NZ_G), tke_b(NZ_G), kdi(NZ_G + 1)
      real(wp) :: nbar2, nb, zp_uncapped, e_col
      real(wp), parameter :: U = 0.02_wp, H2 = 100.0_wp

      n2h = 1.0e-6_wp
      n2f = pz_layer_n2(n2h)
      nb = sqrt(n2h(1))
      nbar2 = sum(n2f)*DZ_G/H_G
      zp_uncapped = 0.0697_wp*(9.61e-4_wp)**2*U/((6.2832e-4_wp)**2*H2*nb*nbar2)
      call check(error, zp_uncapped > 5.0_wp*H_G, "test setup: scale must exceed the cap")
      if (allocated(error)) return
      call pz_run(n2h, .true., U, H2, 0.0_wp, tke_a, kdi, .true.)
      e_col = pz_pbot(U, H2, nb)/(GAMMA_G*MU_G)
      call pz_run(n2h, .false., 0.0_wp, 0.0_wp, e_col, tke_b, kdi, .true.)
      call pz_compare(error, tke_a, tke_b, 1.0e-9_wp, "cap hit != degenerate column")
   end subroutine test_t2_cap

   ! (f) profile = stlaurent: Polzin constants are inert, golden unchanged
   subroutine test_t2_default(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: dz_a(NZ_G), n2_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G), kd_ref(NZ_G + 1)
      real(wp) :: worst
      integer :: k

      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_const_n2_column(ms, NZ_G, DZ_G, N2_CONST_G)
      call setup_tm(tm, grid, NZ_G, E_G, KDMAX_G)
      tm%polzin = .false.
      tm%nu_polzin = 5.0_wp
      tm%nbotref_polzin = 1.0_wp
      tm%polzin_decay_scale_factor = 7.0_wp
      tm%polzin_scale_max_factor = 0.01_wp
      tm%polzin_min_decay_scale = 1.0e6_wp
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
      dz_a = DZ_G
      n2_a = N2_CONST_G
      n2_a(NZ_G) = 0.0_wp
      call ref_sweep(NZ_G, dz_a, n2_a, E_G, GAMMA_G, MU_G, ZETA_G, RHO0_TM, &
                     OMEGA2_GOLD, KDMAX_G, kd_lay, tke_lay)
      call spread_to_interfaces(NZ_G, kd_lay, kd_ref)
      worst = 0.0_wp
      do k = 2, NZ_G
         worst = max(worst, abs(tm%kd_int(3, 3, k) - kd_ref(k))/max(kd_ref(k), 1.0e-30_wp))
      end do
      call check(error, worst < 1.0e-12_wp, "stlaurent profile moved by the polzin constants")
      call tm%destroy()
      call ms%destroy()
   end subroutine test_t2_default

   ! n2_layer_form = "mean_interfaces" on the exponential profile: layer N^2
   ! is the mean of the two bounding jumps (stepped column), vs the legacy value.
   subroutine test_t2_n2_form(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      real(wp) :: n2h(NZ_G - 1), dz_a(NZ_G), kd_lay(NZ_G), tke_lay(NZ_G), kd_ref(NZ_G + 1)
      real(wp) :: worst, legacy15
      integer :: pass, k

      n2h(1:14) = 1.0e-6_wp
      n2h(15:NZ_G - 1) = 4.0e-6_wp
      dz_a = DZ_G
      call make_grid(grid, 2, 2)
      ms%nz_ml = NZ_G
      call ms%init(grid)
      call fill_n2_profile(ms, NZ_G, DZ_G, n2h)
      do pass = 1, 2
         call setup_tm(tm, grid, NZ_G, E_G, KDMAX_G)
         tm%n2_mean_interfaces = (pass == 2)
         call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
         if (pass == 1) then
            legacy15 = tm%kd_int(3, 3, 15)
         else
            call ref_sweep(NZ_G, dz_a, pz_layer_n2(n2h), E_G, GAMMA_G, MU_G, ZETA_G, RHO0_TM, &
                           OMEGA2_GOLD, KDMAX_G, kd_lay, tke_lay)
            call spread_to_interfaces(NZ_G, kd_lay, kd_ref)
            worst = 0.0_wp
            do k = 2, NZ_G
               worst = max(worst, abs(tm%kd_int(3, 3, k) - kd_ref(k))/max(kd_ref(k), 1.0e-30_wp))
            end do
            call check(error, worst < 1.0e-9_wp, "mean_interfaces: Kd != oracle with mean-of-jumps layer N2")
            if (allocated(error)) return
            call check(error, abs(tm%kd_int(3, 3, 15) - legacy15)/legacy15 > 0.05_wp, &
                       "mean_interfaces must differ from legacy at the N2 step")
         end if
         call tm%destroy()
      end do
      call ms%destroy()
   end subroutine test_t2_n2_form

end module test_ocean_tidal_mixing
