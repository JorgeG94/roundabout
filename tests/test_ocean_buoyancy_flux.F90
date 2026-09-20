!! Units + cross-scheme consistency tests for the SURFACE BUOYANCY FLUX
!! `B_0` that both boundary-layer schemes build their convective scale
!! from.
!!
!! The convention under test
!! -------------------------
!! Roundabout's linear EOS is the density-ANOMALY form
!!
!!     rho = rho_0 + beta_S*(S - S_ref) - alpha_T*(T - T_ref)
!!
!! so `alpha_T = -d rho/dT` and `beta_S = +d rho/dS` are DIMENSIONAL —
!! kg m^-3 per degC and per psu — NOT the fractional `(1/rho) d rho/dT`
!! coefficients (~2e-4 1/K) that the KPP literature writes as `alpha`.
!! Buoyancy is `b = -g*rho'/rho_0`, so a buoyancy flux built from a
!! DIMENSIONAL alpha costs a `1/rho_0`:
!!
!!     B_0 = (g/rho_0)*(alpha_T*F_T - beta_S*F_S)      [m^2/s^3]
!!
!! with `F_T = Q_heat/(rho_0*cp)` [K m/s] and `F_S = Q_salt/rho_0`
!! [psu m/s].  Dropping that `1/rho_0` makes `B_0` a factor `rho_0`
!! (~1035) too large; `w_* = (-B_0*h)^(1/3)` is then `rho_0^(1/3)`
!! ~ 10.1x too large, and every KPP stability function, boundary-layer
!! depth and non-local gate downstream of it is wrong with it.  The
!! error is INVISIBLE at the historical `alpha_T = 1.7e-4` type default,
!! which is numerically the FRACTIONAL coefficient used in a
!! dimensional slot, and so silently cancels the missing `1/rho_0`.
!!
!! Cases:
!!   * `buoyancy_flux_kpp_matches_epbl` — for one surface heat + salt
!!     flux and one dimensional (alpha, beta), KPP's `B_0` and EPBL's
!!     `b0` (each read back from its OWN scheme's output field) agree to
!!     round-off.  EPBL builds it as `g*rho_0*(dSV/dT*F_T + dSV/dS*F_S)`
!!     with `dSV/dT = +alpha/rho_0^2`, so the two are the same quantity
!!     written two ways; they must never drift apart again.
!!   * `buoyancy_flux_convective_wstar_analytic` — neutrally stratified
!!     resting column, uniform cooling, no wind: the whole column is the
!!     boundary layer and `w_* = (-B_0*H)^(1/3)` exactly.  Asserted
!!     against a hand-computed number (see the derivation in the body).
!!   * `buoyancy_flux_dimensional_fractional_invariance` — the two
!!     spellings of the same physics, `(alpha dimensional, rho_0)` and
!!     `(alpha/rho_0 fractional, rho_0 = 1)`, give the identical `B_0`;
!!     and the pre-fix expression `g*(alpha*F_T - beta*F_S)` is exactly
!!     `rho_0` times too big.
module test_ocean_buoyancy_flux
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: EOS_VARIANT_LINEAR
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, SEAWATER_CP
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_compute, EPBL_MSTAR_CONSTANT
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_apply_kpp_overlay, &
                             kpp_surface_buoyancy_flux
   implicit none
   private

   public :: collect_ocean_buoyancy_flux_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.2_wp   !! kg/m^3/K  (REALISTIC, dimensional)
   real(wp), parameter :: BETA_LIN = 0.8_wp    !! kg/m^3/psu (REALISTIC, dimensional)

contains

   subroutine collect_ocean_buoyancy_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("buoyancy_flux_kpp_matches_epbl", test_kpp_matches_epbl), &
                  new_unittest("buoyancy_flux_convective_wstar_analytic", test_wstar_analytic), &
                  new_unittest("buoyancy_flux_dimensional_fractional_invariance", test_invariance) &
                  ]
   end subroutine collect_ocean_buoyancy_flux_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine setup_column(ms, nz, dz, rho_col, t_col)
      !! Resting column of uniform layers: uniform density (so KPP's
      !! bulk-Ri never crosses and the BL is the whole column) and
      !! uniform T/S (EPBL needs the tracer fields to evaluate its
      !! specific-volume derivatives; the linear EOS makes them
      !! T/S/p-independent anyway).
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, rho_col, t_col
      integer :: k

      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%rho_layer(:, :, k) = rho_col
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_col*dz
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*dz
      end do
   end subroutine setup_column

   subroutine setup_vmix(vmix, grid, nz)
      type(ocean_vmix_t), intent(inout) :: vmix
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call vmix%init(grid, nz_ml=nz)
      vmix%rho0 = RHO0
      vmix%eos%variant = EOS_VARIANT_LINEAR
      vmix%eos%rho0 = RHO0
      vmix%eos%alpha_T = ALPHA_LIN
      vmix%eos%beta_S = BETA_LIN
   end subroutine setup_vmix

   subroutine run_kpp(grid, ms, vmix, ss, sf)
      !! One KPP overlay call on device.  mem:separate discipline: every
      !! state object AND the surface-flux companion is mapped before the
      !! kernel, and only COMPONENT arrays are updated back (never the
      !! aggregate derived type — that would overwrite the host
      !! descriptors with device addresses).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$acc enter data copyin(ms, vmix, ss, sf)
      call ms%enter_data()
      call vmix%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth, vmix%b0)
      call sf%exit_data()
      call ss%exit_data()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
   end subroutine run_kpp

   subroutine run_epbl(grid, ms, epbl, ss, sf, dt)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: dt

      !$acc enter data copyin(ms, epbl, ss, sf)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      call epbl_compute(grid, epbl, ms, ss, dt, sf=sf)
      !$acc update self(epbl%b0, epbl%mld)
      call sf%exit_data()
      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, sf)
   end subroutine run_epbl

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_kpp_matches_epbl(error)
      !! Same grid, same column, same (dimensional) alpha / beta, same
      !! surface heat AND salt flux -> the two boundary-layer schemes
      !! must form the SAME surface buoyancy flux.  Each value is read
      !! back from that scheme's own persisted field (`vmix%b0`,
      !! `epbl%b0`), not re-derived here.
      !!
      !! They are not bit-identical: KPP forms `Q/(rho0*cp)` with a
      !! division and `(g/rho0)*alpha`, EPBL multiplies by the
      !! reciprocal `1/(rho0*cp)` and forms `g*rho0*(alpha/rho0^2)`.
      !! Round-off (1e-12 relative) is the right gate.
      !!
      !! BEFORE the units fix, KPP's b0 was `rho0` = 1035x EPBL's and
      !! this test failed by three orders of magnitude.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp), parameter :: DZ = 20.0_wp
      real(wp), parameter :: Q_HEAT = -150.0_wp   ! W/m^2, upward (cooling)
      real(wp), parameter :: Q_SALT = 2.0e-5_wp   ! kg/m^2/s, salting
      real(wp) :: b0_kpp, b0_epbl, scale
      integer :: i_probe, j_probe

      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_column(ms, NZ, DZ, RHO0, 12.0_wp)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call sf%init(grid)
         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

         call setup_vmix(vmix, grid, NZ)
         call run_kpp(grid, ms, vmix, ss, sf)

         call epbl%init(grid, nz_ml=NZ)
         epbl%enable = .true.
         epbl%mstar_scheme = EPBL_MSTAR_CONSTANT
         epbl%mstar_const = 1.2_wp
         epbl%rho0 = RHO0
         epbl%eos%variant = EOS_VARIANT_LINEAR
         epbl%eos%rho0 = RHO0
         epbl%eos%alpha_T = ALPHA_LIN
         epbl%eos%beta_S = BETA_LIN
         epbl%f_centre = 0.0_wp
         call run_epbl(grid, ms, epbl, ss, sf, 1800.0_wp)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         b0_kpp = vmix%b0(i_probe, j_probe)
         b0_epbl = epbl%b0(i_probe, j_probe)

         call check(error, b0_kpp < 0.0_wp, &
                    "test setup: cooling + salting must give a destabilizing B_0")
         if (allocated(error)) exit checks

         scale = max(abs(b0_kpp), abs(b0_epbl))
         call check(error, abs(b0_kpp - b0_epbl) <= 1.0e-12_wp*scale, &
                    "KPP B_0 and EPBL b0 disagree beyond round-off")
         if (allocated(error)) exit checks

         ! And both must equal the convention the helper documents.
         call check(error, abs(b0_kpp - kpp_surface_buoyancy_flux( &
                               ALPHA_LIN, BETA_LIN, RHO0, &
                               Q_HEAT/(RHO0*SEAWATER_CP), Q_SALT/RHO0)) <= &
                    1.0e-14_wp*scale, &
                    "KPP B_0 is not (g/rho0)*(alpha*F_T - beta*F_S)")
      end block checks

      call sf%destroy(); call ss%destroy()
      call epbl%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_kpp_matches_epbl

   subroutine test_wstar_analytic(error)
      !! ANALYTIC convective case.  Neutrally stratified (uniform rho),
      !! resting, no wind, uniform surface cooling over a flat 200 m
      !! column of 10 x 20 m layers.  With no density contrast the bulk
      !! Richardson number never reaches `ri_crit`, so `h_b` is the full
      !! column; with `u_* = 0` the KPP velocity scale is purely
      !! convective and the overlay reduces to
      !!
      !!     kv(k) = h_b * w_* * G(sigma),  G(s) = s*(1-s)^2
      !!     w_*   = (-B_0 * h_b)^(1/3)
      !!
      !! so `w_*` is recoverable exactly from the scheme's own `kv`.
      !!
      !! Hand derivation (g = 9.80665, cp = 3992, rho_0 = 1035,
      !! alpha_T = 0.2 kg/m^3/K DIMENSIONAL, Q = -100 W/m^2, H = 200 m):
      !!
      !!   F_T = -100/(1035*3992)            = -2.4202995e-5  K m/s
      !!   B_0 = (9.80665/1035)*0.2*F_T      = -4.586480e-8   m^2/s^3
      !!   w_* = (4.586480e-8 * 200)^(1/3)   =  0.0209333     m/s
      !!   kv  = 200 * 0.0209333 * G(0.3)    =  0.615439      m^2/s
      !!         (interface k = 8, depth 60 m, sigma = 0.3, G = 0.147)
      !!
      !! On UNFIXED main the missing `1/rho_0` gives
      !!   B_0 = -4.747e-5, w_* = 0.211748 m/s, kv = 6.2254 m^2/s
      !! i.e. exactly rho_0^(1/3) = 10.1153 too large in w_* and in kv,
      !! which is how this test fails on the pre-fix code.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 10
      integer, parameter :: K_PROBE = 8          ! interface at 60 m depth
      real(wp), parameter :: DZ = 20.0_wp
      real(wp), parameter :: H_COL = real(NZ, wp)*DZ
      real(wp), parameter :: Q_HEAT = -100.0_wp
      real(wp), parameter :: SIGMA_PROBE = 0.3_wp
      real(wp), parameter :: G_PROBE = SIGMA_PROBE*(1.0_wp - SIGMA_PROBE)**2
      !! Hand-computed values quoted in the derivation above; the loose
      !! 1e-5 relative gate here is a check on the COMMENT, while the
      !! tight gate below is the real assertion.
      real(wp), parameter :: WSTAR_HAND = 0.0209333_wp
      real(wp), parameter :: B0_HAND = -4.586480e-8_wp
      real(wp) :: f_t, b0_expect, wstar_expect, kv_expect
      real(wp) :: wstar_run, kv_run
      integer :: i_probe, j_probe

      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_column(ms, NZ, DZ, RHO0, 12.0_wp)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call sf%init(grid)
         call sf%set_surface_flux_const(Q_HEAT, 0.0_wp)
         call setup_vmix(vmix, grid, NZ)
         call run_kpp(grid, ms, vmix, ss, sf)

         f_t = Q_HEAT/(RHO0*SEAWATER_CP)
         b0_expect = (GRAVITY/RHO0)*ALPHA_LIN*f_t
         wstar_expect = (-b0_expect*H_COL)**(1.0_wp/3.0_wp)
         kv_expect = H_COL*wstar_expect*G_PROBE

         ! The hand-written numbers in the docstring must still be true.
         call check(error, abs(b0_expect - B0_HAND) <= 1.0e-5_wp*abs(B0_HAND), &
                    "hand-computed B_0 in the test comment has gone stale")
         if (allocated(error)) exit checks
         call check(error, abs(wstar_expect - WSTAR_HAND) <= 1.0e-5_wp*WSTAR_HAND, &
                    "hand-computed w_* in the test comment has gone stale")
         if (allocated(error)) exit checks

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2

         ! The whole (neutral) column is the boundary layer.
         call check(error, abs(vmix%bl_depth(i_probe, j_probe) - H_COL) <= &
                    1.0e-10_wp*H_COL, &
                    "neutral column: KPP BL depth is not the full column")
         if (allocated(error)) exit checks

         call check(error, abs(vmix%b0(i_probe, j_probe) - b0_expect) <= &
                    1.0e-13_wp*abs(b0_expect), &
                    "KPP B_0 != (g/rho0)*alpha*Q/(rho0*cp)")
         if (allocated(error)) exit checks

         kv_run = vmix%kv(i_probe, j_probe, K_PROBE)
         call check(error, abs(kv_run - kv_expect) <= 1.0e-10_wp*kv_expect, &
                    "convective kv != h_b*w_**G(sigma) for the analytic w_*")
         if (allocated(error)) exit checks

         ! Recover the scheme's own w_* and gate it on the analytic one.
         wstar_run = kv_run/(H_COL*G_PROBE)
         call check(error, abs(wstar_run - wstar_expect) <= 1.0e-10_wp*wstar_expect, &
                    "recovered convective velocity scale != (-B_0*h)^(1/3)")
      end block checks

      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_wstar_analytic

   subroutine test_invariance(error)
      !! Convention documentation, executable.
      !!
      !! (a) The dimensional pair `(alpha_T, beta_S, rho_0)` and the
      !!     fractional pair `(alpha_T/rho_0, beta_S/rho_0, 1)` are two
      !!     spellings of the same physics and give the same `B_0`.
      !! (b) The pre-fix expression `g*(alpha*F_T - beta*F_S)` — which
      !!     is what the kernel used to compute — is exactly `rho_0`
      !!     times the correct value.  Kept as an executable record of
      !!     the magnitude of the bug this test file was added for.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: F_T = -2.4202995e-5_wp   ! K m/s
      real(wp), parameter :: F_S = 1.9323671e-8_wp    ! psu m/s
      real(wp) :: b0_dim, b0_frac, b0_prefix

      checks: block
         b0_dim = kpp_surface_buoyancy_flux(ALPHA_LIN, BETA_LIN, RHO0, F_T, F_S)
         b0_frac = kpp_surface_buoyancy_flux(ALPHA_LIN/RHO0, BETA_LIN/RHO0, &
                                             1.0_wp, F_T, F_S)

         call check(error, b0_dim < 0.0_wp, "test setup: expected a cooling B_0")
         if (allocated(error)) exit checks

         call check(error, abs(b0_dim - b0_frac) <= 1.0e-14_wp*abs(b0_dim), &
                    "dimensional and fractional alpha spellings disagree")
         if (allocated(error)) exit checks

         ! The bug, written out: no 1/rho_0.
         b0_prefix = GRAVITY*(ALPHA_LIN*F_T - BETA_LIN*F_S)
         call check(error, abs(b0_prefix - RHO0*b0_dim) <= 1.0e-13_wp*abs(b0_prefix), &
                    "pre-fix B_0 expression is not exactly rho_0 times the correct one")
      end block checks
   end subroutine test_invariance

end module test_ocean_buoyancy_flux
