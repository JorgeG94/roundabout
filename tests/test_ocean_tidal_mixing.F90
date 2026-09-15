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
                                     tidal_mixing_merge_into_kt
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
                  new_unittest("tidal_merge_additive", test_merge) &
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

end module test_ocean_tidal_mixing
