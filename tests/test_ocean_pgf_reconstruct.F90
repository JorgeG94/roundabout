!! Analytical unit tests for the FV_MOM6 in-layer T/S reconstruction
!! (RECONSTRUCT_FOR_PRESSURE) of the ocean pressure-force kernel.
!!
!! The kernel replaces the layer-mean (PCM) density integral with a
!! 5-point Boole quadrature of a monotone PLM/PPM sub-layer T/S profile
!! (Adcroft, Hallberg & Harrison 2008; White, Adcroft & Hallberg 2009).
!!
!! The oracle uses the LINEAR EOS (rho = rho0 + beta_S(S-S_ref) -
!! alpha_T(T-T_ref), no pressure dependence).  With linear T(z), S(z) the
!! in-situ density rho(z) is linear in z, so:
!!   * the true layer integral int rho' dz and its first moment are
!!     CLOSED-FORM (a linear integrand) — the analytic oracle;
!!   * the 5-point Boole rule is EXACT for a linear integrand, and the
!!     interior-layer PLM/PPM edges of a monotone linear profile recover
!!     the exact endpoint values, so the reconstructed dpa/intz_dpa hit
!!     the analytic answer to round-off;
!!   * the PCM (layer-mean) intz_dpa = 0.5*dpa*h is the moment of a
!!     CONSTANT density, biased on a sloped layer — strictly worse.
!!
!! Tests:
!!   1. flat_strat_consistency  — uniform T/S => rec dpa/intz == PCM ==
!!      analytic to ~1e-12 (incompressible linear EOS; the N8 recipe).
!!   2/3. error_reduction_sloped (PLM, PPM) — thick sloped interior layer,
!!      linear T(z)/S(z): the reconstructed intz_dpa error vs the analytic
!!      oracle is strictly smaller than the PCM error.
!!   4. default_off_bit_identity — reconstruct OFF reproduces the baseline
!!      FV_MOM6 PFu/PFv bit-for-bit.
!!   5. quiescent_uniform_rho   — uniform T/S + sloped bathy + reconstruct
!!      ON => zero PGF (no spurious acceleration on a neutral column).
module test_ocean_pgf_reconstruct
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_reconstruct_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZTEST = 5
      !! Layer count used by every test column (thin bed, thick interior,
      !! three thin surface layers — the thick layer is interior so the
      !! PLM/PPM reconstruction acts on it rather than falling back to PCM).

   ! Linear-EOS coefficients (the test is analytic so values only set scale).
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: RHO_REF = 1030.0_wp
   real(wp), parameter :: T_REF = 10.0_wp, S_REF = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3 per degC
   real(wp), parameter :: BETA_S = 0.78_wp     ! kg/m^3 per PSU

contains

   subroutine collect_ocean_pgf_reconstruct_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("flat_strat_consistency", test_flat_strat), &
                  new_unittest("error_reduction_sloped_plm", test_error_reduction_plm), &
                  new_unittest("error_reduction_sloped_ppm", test_error_reduction_ppm), &
                  new_unittest("ppm_curvature_cubic_vs_plm", test_ppm_curvature), &
                  new_unittest("default_off_bit_identity", test_default_off), &
                  new_unittest("quiescent_uniform_rho_zero_pgf", test_quiescent) &
                  ]
   end subroutine collect_ocean_pgf_reconstruct_tests

   pure function rho_lin(T, S) result(r)
      !! Linear EOS absolute density (mirror of eos_density_point LINEAR).
      real(wp), intent(in) :: T, S
      real(wp) :: r
      r = RHO0 + BETA_S*(S - S_REF) - ALPHA_T*(T - T_REF)
   end function rho_lin

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%T_ref = T_REF
      eos%S_ref = S_REF
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%is_init = .true.
   end subroutine make_eos

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine seed_linear_TS(ms, dTdz, dSdz, b)
      !! Fill h_layer + tracer hTr from a linear continuous T(z)=T_REF+dTdz*z,
      !! S(z)=S_REF+dSdz*z (z surface-relative, <=0), with a per-column
      !! bathymetry b(i,j) so interfaces slope across faces.  Layer thicknesses
      !! are split sigma-style (equal fraction of the local depth), and the
      !! layer-mean T/S equals the mid-depth value (exact mean of a linear
      !! profile).  Two interior thin layers bracket a thick interior layer.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dTdz, dSdz
      real(wp), intent(in) :: b(:, :)
      integer :: i, j, k, nx, ny, nz
      real(wp) :: depth, e_lo, e_hi, z_mid, Tk, Sk
      real(wp) :: frac(NZTEST)
      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      ! Fixed thickness fractions per layer (k=1 bed .. k=nz surface).
      ! A THICK interior layer (k=2) so PLM/PPM act on it (not a boundary).
      ! nz = 5: [thin bed, THICK, thin, thin, thin surface].
      frac(1) = 0.04_wp
      frac(2) = 0.70_wp
      frac(3) = 0.10_wp
      frac(4) = 0.08_wp
      frac(5) = 0.08_wp
      do j = 1, ny
         do i = 1, nx
            depth = b(i, j)
            e_lo = -depth                      ! bed (deepest, k=1 bottom)
            do k = 1, nz
               ms%h_layer(i, j, k) = frac(k)*depth
               e_hi = e_lo + ms%h_layer(i, j, k)
               z_mid = 0.5_wp*(e_lo + e_hi)
               Tk = T_REF + dTdz*z_mid
               Sk = S_REF + dSdz*z_mid
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = Tk*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = Sk*ms%h_layer(i, j, k)
               ms%rho_layer(i, j, k) = rho_lin(Tk, Sk)
               e_lo = e_hi
            end do
         end do
      end do
   end subroutine seed_linear_TS

   subroutine run_pgf(grid, ms, pgf, eos, dx)
      !! One compute pass; pulls results back to host.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: dx
      type(ocean_metrics_t) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data, &
      !$acc&            pgf%intz_dpa%data, pgf%pa%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   pure function analytic_intz_dpa(e_top, dz, T_top, T_bot, S_top, S_bot) result(intz)
      !! Closed-form first moment (from the TOP edge inward) of a LINEAR
      !! density anomaly across a layer.  rho'(s) varies linearly from
      !! rho'_top (s=0, top edge) to rho'_bot (s=1, bottom edge).  The
      !! kernel's intz_dpa = 0.5*g*dz^2 * bracket where bracket =
      !! 2*int_0^1 rho'(s)*(1-s) ds (moment weighted from the bottom edge,
      !! matching the Boole formula).  For linear rho'(s) =
      !! rho'_top + (rho'_bot-rho'_top)*s:
      !!   2*int_0^1 (a + b*s)(1-s) ds = a + b/3   (a=rho'_top, b=rho'_bot-rho'_top)
      !! => intz = 0.5*g*dz^2*(rho'_top + (rho'_bot-rho'_top)/3).
      real(wp), intent(in) :: e_top, dz, T_top, T_bot, S_top, S_bot
      real(wp) :: intz, rp_top, rp_bot, aa, bb
      rp_top = rho_lin(T_top, S_top) - RHO_REF
      rp_bot = rho_lin(T_bot, S_bot) - RHO_REF
      aa = rp_top
      bb = rp_bot - rp_top
      intz = 0.5_wp*GRAVITY*dz*dz*(aa + bb/3.0_wp)
      if (.false.) intz = intz + 0.0_wp*e_top   ! e_top unused (linear EOS)
   end function analytic_intz_dpa

   subroutine test_flat_strat(error)
      !! Uniform T/S in every layer + linear (incompressible) EOS => the
      !! in-layer density is truly constant, so the reconstructed dpa /
      !! intz_dpa must equal the PCM (layer-mean) values AND the analytic
      !! oracle to round-off (the N8 recipe — never assert rec==pcm under
      !! a compressible EOS).
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_rec, pgf_pcm
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 2000.0_wp, DEPTH = 1500.0_wp
      real(wp), parameter :: TUNI = 12.0_wp, SUNI = 34.5_wp
      real(wp), allocatable :: b(:, :), intz_pcm(:, :, :)
      real(wp) :: max_rec_pcm, max_rec_ana, intz_ana
      integer :: i, j, k, nx, ny, nz, ng
      integer :: ksel
      checks: block
         call make_grid(grid, 8, 6, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ng = NGHOST

         allocate (b(nx, ny), source=DEPTH)
         ! Uniform T/S everywhere (override the linear seeding).
         call seed_linear_TS(ms, 0.0_wp, 0.0_wp, b)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = TUNI*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = SUNI*ms%h_layer(i, j, k)
                  ms%rho_layer(i, j, k) = rho_lin(TUNI, SUNI)
               end do
            end do
         end do

         call pgf_rec%init(grid, nz_ml=nz)
         pgf_rec%variant = OPGF_VARIANT_FV_MOM6
         pgf_rec%rho0 = RHO0
         pgf_rec%rho_ref = RHO_REF
         pgf_rec%reconstruct_for_pressure = .true.
         pgf_rec%recon_scheme = 1
         call pgf_rec%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_rec, eos, DX)

         call pgf_pcm%init(grid, nz_ml=nz)
         pgf_pcm%variant = OPGF_VARIANT_FV_MOM6
         pgf_pcm%rho0 = RHO0
         pgf_pcm%rho_ref = RHO_REF
         pgf_pcm%reconstruct_for_pressure = .false.
         call pgf_pcm%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_pcm, eos, DX)

         allocate (intz_pcm, source=pgf_pcm%intz_dpa%data)
         max_rec_pcm = maxval(abs(pgf_rec%intz_dpa%data - intz_pcm))

         ! Analytic: uniform rho' => intz = 0.5*g*dz^2*rho'.  Compare the
         ! thick interior layer k=2 at an interior column.
         ksel = 2
         max_rec_ana = 0.0_wp
         do j = ng + 1, ng + 6
            do i = ng + 1, ng + 8
               intz_ana = 0.5_wp*GRAVITY*ms%h_layer(i, j, ksel)**2 &
                          *(rho_lin(TUNI, SUNI) - RHO_REF)
               max_rec_ana = max(max_rec_ana, &
                                 abs(pgf_rec%intz_dpa%data(i, j, ksel) - intz_ana))
            end do
         end do

         call check(error, max_rec_pcm < 1.0e-9_wp, &
                    "flat strat: reconstructed intz_dpa != PCM")
         if (allocated(error)) exit checks
         call check(error, max_rec_ana < 1.0e-6_wp, &
                    "flat strat: reconstructed intz_dpa != analytic")
      end block checks
      if (allocated(b)) deallocate (b)
      if (allocated(intz_pcm)) deallocate (intz_pcm)
      call pgf_rec%destroy()
      call pgf_pcm%destroy()
      call ms%destroy()
   end subroutine test_flat_strat

   subroutine error_reduction_core(error, scheme)
      !! Thick sloped interior layer with linear T(z)/S(z): the
      !! reconstructed intz_dpa hits the analytic linear-density moment to
      !! round-off, while PCM (constant-density moment) carries an O(slope)
      !! bias => rec error strictly < pcm error.  Headline assertion.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_rec, pgf_pcm
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DTDZ = 0.012_wp, DSDZ = 0.0006_wp
      real(wp), allocatable :: b(:, :)
      real(wp) :: e_lo, e_hi, e_top, dz, z_t, z_b
      real(wp) :: T_top, T_bot, S_top, S_bot, intz_ana
      real(wp) :: err_rec, err_pcm
      integer :: i, j, k, nx, ny, nz, ng, ksel
      checks: block
         call make_grid(grid, 8, 6, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ng = NGHOST

         ! Sloped bathymetry across x: deeper to the east => interfaces tilt.
         allocate (b(nx, ny), source=0.0_wp)
         do j = 1, ny
            do i = 1, nx
               b(i, j) = 1200.0_wp + 40.0_wp*real(i - ng, wp)
            end do
         end do
         call seed_linear_TS(ms, DTDZ, DSDZ, b)

         call pgf_rec%init(grid, nz_ml=nz)
         pgf_rec%variant = OPGF_VARIANT_FV_MOM6
         pgf_rec%rho0 = RHO0
         pgf_rec%rho_ref = RHO_REF
         pgf_rec%reconstruct_for_pressure = .true.
         pgf_rec%recon_scheme = scheme
         call pgf_rec%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_rec, eos, DX)

         call pgf_pcm%init(grid, nz_ml=nz)
         pgf_pcm%variant = OPGF_VARIANT_FV_MOM6
         pgf_pcm%rho0 = RHO0
         pgf_pcm%rho_ref = RHO_REF
         pgf_pcm%reconstruct_for_pressure = .false.
         call pgf_pcm%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_pcm, eos, DX)

         ! Compare intz_dpa on the THICK interior layer (k=2) at interior
         ! columns vs the analytic linear-density moment.
         ksel = 2
         err_rec = 0.0_wp
         err_pcm = 0.0_wp
         do j = ng + 1, ng + 6
            do i = ng + 1, ng + 8
               ! Rebuild the layer geometry for column (i,j).
               e_lo = -b(i, j)
               do k = 1, ksel - 1
                  e_lo = e_lo + ms%h_layer(i, j, k)
               end do
               dz = ms%h_layer(i, j, ksel)
               e_hi = e_lo + dz
               e_top = e_hi          ! shallower interface
               z_t = e_top           ! top edge depth
               z_b = e_lo            ! bottom edge depth
               T_top = T_REF + DTDZ*z_t
               T_bot = T_REF + DTDZ*z_b
               S_top = S_REF + DSDZ*z_t
               S_bot = S_REF + DSDZ*z_b
               intz_ana = analytic_intz_dpa(e_top, dz, T_top, T_bot, S_top, S_bot)
               err_rec = max(err_rec, abs(pgf_rec%intz_dpa%data(i, j, ksel) - intz_ana))
               err_pcm = max(err_pcm, abs(pgf_pcm%intz_dpa%data(i, j, ksel) - intz_ana))
            end do
         end do

         call check(error, err_pcm > 1.0e3_wp, &
                    "error reduction: PCM intz_dpa error should be sizeable")
         if (allocated(error)) exit checks
         call check(error, err_rec < 0.1_wp*err_pcm, &
                    "error reduction: reconstructed intz_dpa not >=10x better than PCM")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf_rec%destroy()
      call pgf_pcm%destroy()
      call ms%destroy()
   end subroutine error_reduction_core

   subroutine test_error_reduction_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call error_reduction_core(error, 1)
   end subroutine test_error_reduction_plm

   subroutine test_error_reduction_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call error_reduction_core(error, 2)
   end subroutine test_error_reduction_ppm

   subroutine seed_cubic_TS(ms, dTdz, dSdz, c2T, c2S, c3T, c3S, depth)
      !! Uniform-thickness column with a CUBIC continuous profile
      !! T(z)=T_REF+dTdz*z+c2T*z^2+c3T*z^3 (z surface-relative, <=0); S alike.
      !! The layer-mean seeded into hTr is the EXACT layer-AVERAGE over
      !! [z_mid-dz/2, z_mid+dz/2]:
      !!   <z^2> = z_mid^2 + dz^2/12,  <z^3> = z_mid^3 + z_mid*dz^2/4.
      !! A CUBIC is the smallest profile that distinguishes PPM from PLM in
      !! the first moment intz_dpa: for a quadratic, PLM's centred slope on
      !! exact cell-means is exactly b+c, which makes the linear first moment
      !! equal the quadratic moment (degenerate).  A cubic breaks that, so
      !! PPM's parabola (the q6 term) is strictly closer to the analytic
      !! cubic moment than PLM's line.  rho is affine in (T,S) under the
      !! linear EOS, so the layer-average density = rho_lin(<T>,<S>).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dTdz, dSdz, c2T, c2S, c3T, c3S, depth
      integer :: i, j, k, nx, ny, nz
      real(wp) :: dz, e_lo, z_mid, Tk, Sk, z2m, z3m
      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      dz = depth/real(nz, wp)            ! uniform layers
      do j = 1, ny
         do i = 1, nx
            e_lo = -depth                  ! bed (k=1 bottom)
            do k = 1, nz
               ms%h_layer(i, j, k) = dz
               z_mid = e_lo + 0.5_wp*dz
               z2m = z_mid*z_mid + dz*dz/12.0_wp            ! exact <z^2>
               z3m = z_mid*z_mid*z_mid + z_mid*dz*dz/4.0_wp ! exact <z^3>
               Tk = T_REF + dTdz*z_mid + c2T*z2m + c3T*z3m
               Sk = S_REF + dSdz*z_mid + c2S*z2m + c3S*z3m
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = Tk*dz
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = Sk*dz
               ms%rho_layer(i, j, k) = rho_lin(Tk, Sk)
               e_lo = e_lo + dz
            end do
         end do
      end do
   end subroutine seed_cubic_TS

   pure function analytic_intz_dpa_cubic(z_t, dz, dTdz, dSdz, c2T, c2S, c3T, c3S) result(intz)
      !! Closed-form kernel intz_dpa = 0.5*g*dz^2 * [2*int_0^1 rho'(s)(1-s)ds]
      !! for a CUBIC-in-z density anomaly.  Under the linear EOS,
      !! rho'(z) = (RHO0-RHO_REF) + p*z + q*z^2 + r*z^3 with
      !!   p = BETA_S*dSdz - ALPHA_T*dTdz,  q = BETA_S*c2S - ALPHA_T*c2T,
      !!   r = BETA_S*c3S - ALPHA_T*c3T.
      !! Map z(s)=z_t-dz*s (s=0 top edge .. s=1 bottom edge), expand to
      !!   rho'(s) = A + B*s + C*s^2 + D*s^3,
      !!   A = rho'(z_t),  B = -dz(p+2q z_t+3r z_t^2),
      !!   C = dz^2(q+3r z_t),  D = -r dz^3.
      !! 2*int_0^1 (A+Bs+Cs^2+Ds^3)(1-s) ds = A + B/3 + C/6 + D/10
      !! (the integrand is quartic in s; the 5-pt Boole rule is exact for it).
      real(wp), intent(in) :: z_t, dz, dTdz, dSdz, c2T, c2S, c3T, c3S
      real(wp) :: intz, p, q, r, AA, BB, CC, DD
      p = BETA_S*dSdz - ALPHA_T*dTdz
      q = BETA_S*c2S - ALPHA_T*c2T
      r = BETA_S*c3S - ALPHA_T*c3T
      AA = (RHO0 - RHO_REF) + p*z_t + q*z_t*z_t + r*z_t*z_t*z_t
      BB = -dz*(p + 2.0_wp*q*z_t + 3.0_wp*r*z_t*z_t)
      CC = dz*dz*(q + 3.0_wp*r*z_t)
      DD = -r*dz*dz*dz
      intz = 0.5_wp*GRAVITY*dz*dz*(AA + BB/3.0_wp + CC/6.0_wp + DD/10.0_wp)
   end function analytic_intz_dpa_cubic

   subroutine test_ppm_curvature(error)
      !! Cubic T(z)/S(z) on a uniform column: PPM (recon_scheme=2) fits a
      !! parabola (using the q6 = 3*(2*mean-(top+bot)) curvature term) and is
      !! strictly closer to the analytic cubic first moment than PLM's line.
      !! This is the regression guard on q6 that the existing LINEAR-profile
      !! tests leave at q6==0 (and which a QUADRATIC profile cannot exercise
      !! either: PLM's centred slope makes its first moment exact for a
      !! quadratic -- see seed_cubic_TS).  Compared on a PURE-interior layer
      !! (k=3 of 5 uniform layers; both neighbours interior).  A broken q6
      !! makes PPM no better than -- or worse than -- PLM, failing the
      !! ratio check.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_ppm, pgf_plm
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 2000.0_wp, DEPTH = 1200.0_wp
      real(wp), parameter :: DTDZ = 0.01_wp, DSDZ = 0.0003_wp
      real(wp), parameter :: C2T = 4.0e-6_wp, C2S = -1.5e-7_wp    ! z^2
      real(wp), parameter :: C3T = 2.0e-9_wp, C3S = -1.0e-10_wp   ! z^3
      real(wp), allocatable :: b(:, :)
      real(wp) :: dz, z_t, intz_ana, err_ppm, err_plm
      integer :: i, j, nx, ny, nz, ng, ksel
      checks: block
         call make_grid(grid, 8, 6, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ng = NGHOST
         allocate (b(nx, ny), source=DEPTH)
         call seed_cubic_TS(ms, DTDZ, DSDZ, C2T, C2S, C3T, C3S, DEPTH)

         call pgf_ppm%init(grid, nz_ml=nz)
         pgf_ppm%variant = OPGF_VARIANT_FV_MOM6
         pgf_ppm%rho0 = RHO0
         pgf_ppm%rho_ref = RHO_REF
         pgf_ppm%reconstruct_for_pressure = .true.
         pgf_ppm%recon_scheme = 2          ! PPM
         call pgf_ppm%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_ppm, eos, DX)

         call pgf_plm%init(grid, nz_ml=nz)
         pgf_plm%variant = OPGF_VARIANT_FV_MOM6
         pgf_plm%rho0 = RHO0
         pgf_plm%rho_ref = RHO_REF
         pgf_plm%reconstruct_for_pressure = .true.
         pgf_plm%recon_scheme = 1          ! PLM
         call pgf_plm%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_plm, eos, DX)

         ! Pure interior layer k=3 (uniform layers => dz = DEPTH/nz).  Top
         ! (shallower) edge z_t = -DEPTH + ksel*dz.
         ksel = 3
         dz = DEPTH/real(nz, wp)
         z_t = -DEPTH + real(ksel, wp)*dz
         intz_ana = analytic_intz_dpa_cubic(z_t, dz, DTDZ, DSDZ, C2T, C2S, C3T, C3S)

         err_ppm = 0.0_wp
         err_plm = 0.0_wp
         do j = ng + 1, ng + 6
            do i = ng + 1, ng + 8
               err_ppm = max(err_ppm, abs(pgf_ppm%intz_dpa%data(i, j, ksel) - intz_ana))
               err_plm = max(err_plm, abs(pgf_plm%intz_dpa%data(i, j, ksel) - intz_ana))
            end do
         end do

         ! PLM carries a genuine O(curvature) bias on the cubic profile.
         call check(error, err_plm > 1.0_wp, &
                    "ppm curvature: PLM error on cubic should be sizeable")
         if (allocated(error)) exit checks
         ! PPM's parabola (the q6 term) is much closer to the cubic moment.
         ! A broken/zeroed q6 collapses PPM toward PLM and fails this.
         call check(error, err_ppm < 0.2_wp*err_plm, &
                    "ppm curvature: PPM (q6) not closer to analytic cubic than PLM")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf_ppm%destroy()
      call pgf_plm%destroy()
      call ms%destroy()
   end subroutine test_ppm_curvature

   subroutine test_default_off(error)
      !! reconstruct_for_pressure = .false. must reproduce the baseline
      !! FV_MOM6 PFu/PFv bit-for-bit (regression gate; the orchestrator
      !! relies on this default-OFF bit-identity).  We compare the
      !! default-constructed pgf (reconstruct off) against an explicitly
      !! PCM pgf — both must give identical fields on a stratified sloped
      !! column.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), allocatable :: b(:, :), dpdx_a(:, :, :), dpdy_a(:, :, :)
      real(wp) :: max_diff
      integer :: i, j, nx, ny, nz, ng
      checks: block
         call make_grid(grid, 8, 6, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ng = NGHOST
         allocate (b(nx, ny), source=0.0_wp)
         do j = 1, ny
            do i = 1, nx
               b(i, j) = 1200.0_wp + 40.0_wp*real(i - ng, wp)
            end do
         end do
         call seed_linear_TS(ms, 0.012_wp, 0.0006_wp, b)

         ! pgf_a: default (reconstruct off).
         call pgf_a%init(grid, nz_ml=nz)
         pgf_a%variant = OPGF_VARIANT_FV_MOM6
         pgf_a%rho0 = RHO0
         pgf_a%rho_ref = RHO_REF
         call pgf_a%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_a, eos, DX)
         allocate (dpdx_a, source=pgf_a%dpdx_face%data)
         allocate (dpdy_a, source=pgf_a%dpdy_face%data)

         ! pgf_b: explicitly off, same setup.
         call pgf_b%init(grid, nz_ml=nz)
         pgf_b%variant = OPGF_VARIANT_FV_MOM6
         pgf_b%rho0 = RHO0
         pgf_b%rho_ref = RHO_REF
         pgf_b%reconstruct_for_pressure = .false.
         call pgf_b%set_bathymetry(b)
         call run_pgf(grid, ms, pgf_b, eos, DX)

         max_diff = max(maxval(abs(pgf_b%dpdx_face%data - dpdx_a)), &
                        maxval(abs(pgf_b%dpdy_face%data - dpdy_a)))
         call check(error, max_diff == 0.0_wp, &
                    "default-off: reconstruct flag not bit-identical to PCM baseline")
      end block checks
      if (allocated(b)) deallocate (b)
      if (allocated(dpdx_a)) deallocate (dpdx_a)
      if (allocated(dpdy_a)) deallocate (dpdy_a)
      call pgf_a%destroy()
      call pgf_b%destroy()
      call ms%destroy()
   end subroutine test_default_off

   subroutine test_quiescent(error)
      !! Uniform T/S (=> uniform rho) + sloped bathy + reconstruct ON.
      !! With a flat free surface (b = sum_k h_layer per column) and no
      !! density structure, the true PGF is identically zero; the
      !! reconstruction must not inject any spurious acceleration on the
      !! neutral column (latent-bug net).
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 2000.0_wp, TUNI = 11.0_wp, SUNI = 35.2_wp
      real(wp), allocatable :: b(:, :)
      real(wp) :: pgf_max
      integer :: i, j, k, nx, ny, nz, ng
      checks: block
         call make_grid(grid, 10, 6, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         nz = ms%nz_ml
         ng = NGHOST

         ! Sloped bathy; uniform T/S everywhere => rho uniform.
         allocate (b(nx, ny), source=0.0_wp)
         do j = 1, ny
            do i = 1, nx
               b(i, j) = 1000.0_wp + 60.0_wp*real(i - ng, wp)
            end do
         end do
         call seed_linear_TS(ms, 0.0_wp, 0.0_wp, b)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = TUNI*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = SUNI*ms%h_layer(i, j, k)
                  ms%rho_layer(i, j, k) = rho_lin(TUNI, SUNI)
               end do
            end do
         end do

         ! b = sum_k h_layer per column => eta = 0 (flat free surface).  The
         ! seeding already used b as the column depth, so sum_k h = b.
         call pgf%init(grid, nz_ml=nz)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = RHO0
         pgf%rho_ref = rho_lin(TUNI, SUNI)   ! anomaly vanishes => no PGF
         pgf%reconstruct_for_pressure = .true.
         pgf%recon_scheme = 1
         call pgf%set_bathymetry(b)
         call run_pgf(grid, ms, pgf, eos, DX)

         pgf_max = max(maxval(abs(pgf%dpdx_face%data)), &
                       maxval(abs(pgf%dpdy_face%data)))
         call check(error, pgf_max < 1.0e-10_wp, &
                    "quiescent: reconstruct ON injected a spurious PGF on a neutral column")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_quiescent

end module test_ocean_pgf_reconstruct
