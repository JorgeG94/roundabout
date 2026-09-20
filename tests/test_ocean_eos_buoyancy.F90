!! Tests for the EOS-derived thermal-expansion / haline-contraction
!! coefficients (E4) and for the `&ocean_vmix_nml buoyancy_coeffs` knob
!! that routes them into the KPP surface buoyancy flux and the
!! double-diffusion density ratio.
module test_ocean_eos_buoyancy
   !! What is locked here, and why each case exists.
   !!
   !!   (a) **Correctness of the derivatives.** `eos_buoyancy_coeffs` is
   !!       differentiated in closed form per EOS variant, so the one
   !!       failure mode is an algebra slip in a polynomial.  Centred
   !!       finite differences of THE SAME variant's `eos_density_point`
   !!       catch that and nothing else — a shared coefficient typo would
   !!       cancel, which is exactly right: the contract is "the analytic
   !!       derivative of the density this EOS actually returns".  The
   !!       sweep deliberately includes the COLD CORNER (−2…2 degC,
   !!       33–35 PSU, 0–1e7 Pa), the ice-shelf-cavity regime the
   !!       capability exists for and the one the 10 degC-tuned constants
   !!       are worst in.
   !!
   !!       No published α table is used as an oracle: neither Wright
   !!       (1997, JAOT 14:735-740) nor Roquet et al. (2015, Ocean
   !!       Modelling 90:29-43) is available in this tree's papers folder
   !!       (`/home/jorge/nci/cdx/python_prototypes/papers/` carries the
   !!       ice-shelf set only), and this repository does not ship
   !!       unverified constants — so the finite-difference identity plus
   !!       the qualitative trends below are the whole check.  The
   !!       underlying coefficient tables are themselves already locked
   !!       against published values by `test_ocean_eos_handle` (Wright)
   !!       and `test_ocean_eos_roquet` (Roquet, 128/128 vs the verified
   !!       prototype), so this suite only has to lock the DERIVATIVE.
   !!
   !!   (b) **Physical trends.**  α under Wright must FALL toward the
   !!       freezing point and RISE with pressure (thermobaricity).  Those
   !!       two signs are the entire physical argument for the knob, so
   !!       they are asserted directly rather than inferred.
   !!
   !!   (c) **Bit-identity under the linear EOS.**  Structural, not
   !!       coincidental: the linear branch of `eos_buoyancy_coeffs`
   !!       returns the handle members themselves rather than
   !!       `−ρ²·dSV/dX`, so `buoyancy_coeffs = "eos"` and `"constant"`
   !!       must agree to the LAST BIT — asserted with `==` on the
   !!       coefficients and on whole kv/kt fields out of the real
   !!       kernels.
   !!
   !!   (d) **The routed consumers, on the device path.**  The KPP and
   !!       double-diffusion cases drive the PRODUCTION kernels through
   !!       `enter_data` / `!$acc update self` under the `mem:separate`
   !!       contract, which is also the only correct way to test the new
   !!       `!$acc routine seq` point routine from a test executable: the
   !!       `do concurrent` that calls it lives inside the shared library
   !!       (`vmix_kpp_overlay_impl`, `vmix_split_ddiff_eos_impl`), never
   !!       in this translation unit, because nvlink cannot resolve a
   !!       `routine seq` symbol across the library boundary from a test
   !!       executable's own kernel.  That is the pattern vmix and EPBL
   !!       already use for `eos_specvol_derivs`, followed exactly.
   !!
   !!       `kpp_wright_cold_kv_scales_as_cuberoot_alpha` is the
   !!       quantitative one: with zero wind and an unstratified column
   !!       the KPP scale collapses to `w_s = w_* = (−B_0·h_b)^(1/3)` and
   !!       `B_0 ∝ α`, so the overlay's kv must scale as the CUBE ROOT of
   !!       the α ratio — an analytic relation that pins the routing, the
   !!       device data motion and the coefficient value in one number.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_buoyancy_coeffs, eos_density_derivs, &
                      eos_density_point, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, &
                             vmix_apply_kpp_overlay, vmix_split_kd_heat_salt, &
                             parse_buoyancy_coeffs, &
                             BUOY_COEFFS_CONSTANT, BUOY_COEFFS_EOS, &
                             BUOY_COEFFS_INVALID
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   implicit none
   private

   public :: collect_ocean_eos_buoyancy_tests

   integer, parameter :: NGHOST = 2

   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.17_wp
      !! Dimensional linear-EOS α (kg/m³ per degC) — roughly the Wright
      !! value at 10 degC / 35 PSU / surface, so the "constant" setting
      !! in the cold-corner cases below is the realistic open-ocean
      !! stand-in a cavity run would otherwise inherit.
   real(wp), parameter :: BETA_LIN = 0.78_wp
      !! Dimensional linear-EOS β (kg/m³ per PSU).

   ! Cold-corner + mid-ocean (T, S, p) sweep.  The first nine points are
   ! the cavity corner the capability targets; the last three keep an
   ! open-ocean anchor so a change that only fixed the cold end would
   ! still be caught.
   integer, parameter :: N_PTS = 12
   real(wp), parameter :: T_PTS(N_PTS) = [ &
                          -2.0_wp, -1.9_wp, -1.0_wp, 0.0_wp, 1.0_wp, 2.0_wp, &
                          -2.0_wp, -1.9_wp, 2.0_wp, &
                          10.0_wp, 18.0_wp, 25.0_wp]
   real(wp), parameter :: S_PTS(N_PTS) = [ &
                          33.0_wp, 34.5_wp, 34.0_wp, 35.0_wp, 33.5_wp, 34.5_wp, &
                          34.5_wp, 33.0_wp, 35.0_wp, &
                          35.0_wp, 36.0_wp, 34.0_wp]
   real(wp), parameter :: P_PTS(N_PTS) = [ &
                          0.0_wp, 0.0_wp, 2.5e6_wp, 5.0e6_wp, 7.5e6_wp, 1.0e7_wp, &
                          1.0e7_wp, 5.0e6_wp, 0.0_wp, &
                          0.0_wp, 2.0e6_wp, 4.0e7_wp]

   real(wp), parameter :: EPS_T = 1.0e-4_wp
      !! Centred-difference step in degC.  With ρ ~ 1e3 the round-off
      !! floor of a centred quotient is ~eps·ρ/h ≈ 2e-9 kg/m³/degC, i.e.
      !! ~5e-8 relative even against the SMALL cold-corner α ≈ 0.04, while
      !! the h²·ρ'''/6 truncation term is orders below that — so 1e-7
      !! relative is a genuine agreement bound, not a tuned one.
   real(wp), parameter :: EPS_S = 1.0e-4_wp
      !! Centred-difference step in PSU.
   real(wp), parameter :: FD_RTOL = 1.0e-7_wp
      !! Relative agreement demanded of analytic-vs-finite-difference.
   real(wp), parameter :: FD_ATOL = 1.0e-9_wp
      !! Absolute floor, so a point where α passes near zero is judged on
      !! the round-off floor rather than on a meaningless ratio.

contains

   subroutine collect_ocean_eos_buoyancy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("eos_buoy_wright_matches_finite_difference", test_fd_wright), &
                  new_unittest("eos_buoy_roquet_matches_finite_difference", test_fd_roquet), &
                  new_unittest("eos_buoy_linear_matches_finite_difference", test_fd_linear), &
                  new_unittest("eos_buoy_linear_returns_handle_pair_bitwise", test_linear_bitwise), &
                  new_unittest("eos_buoy_density_derivs_is_exact_sign_twin", test_sign_twin), &
                  new_unittest("eos_buoy_wright_alpha_falls_toward_freezing", test_alpha_cold), &
                  new_unittest("eos_buoy_wright_alpha_rises_with_pressure", test_alpha_pressure), &
                  new_unittest("eos_buoy_parse_round_trip_and_fail_loud", test_parse), &
                  new_unittest("kpp_linear_eos_identical_under_both_settings", test_kpp_linear), &
                  new_unittest("kpp_wright_cold_kv_scales_as_cuberoot_alpha", test_kpp_wright), &
                  new_unittest("ddiff_linear_eos_identical_under_both_settings", test_ddiff_linear), &
                  new_unittest("ddiff_wright_cold_changes_the_split", test_ddiff_wright) &
                  ]
   end subroutine collect_ocean_eos_buoyancy_tests

   ! -----------------------------------------------------------------
   ! Handle constructors
   ! -----------------------------------------------------------------

   pure function make_eos(variant) result(eos)
      !! One constructor for all three variants so the linear members are
      !! present (and identical) on every handle — which is what makes
      !! "the linear branch returns THESE numbers" a meaningful claim.
      integer, intent(in) :: variant
      type(eos_t) :: eos
      eos%variant = variant
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_LIN
      eos%beta_S = BETA_LIN
      eos%T_ref = 10.0_wp
      eos%S_ref = 35.0_wp
      eos%p_ref = 0.0_wp
   end function make_eos

   ! -----------------------------------------------------------------
   ! (a) Analytic derivatives vs centred finite differences
   ! -----------------------------------------------------------------

   subroutine check_fd_for(error, eos, rtol, label)
      !! Sweep `T_PTS`/`S_PTS`/`P_PTS` and compare `eos_density_derivs`
      !! against centred differences of the SAME handle's
      !! `eos_density_point`.  Density, not specific volume: dSV is
      !! already covered by `test_ocean_epbl`, and ρ is the quantity the
      !! α/β convention is defined against.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rtol
      character(len=*), intent(in) :: label
      real(wp) :: drho_dt, drho_ds, fd_dt, fd_ds, rho_p, rho_m
      integer :: n
      checks: block
         do n = 1, N_PTS
            call eos_density_derivs(eos, T_PTS(n), S_PTS(n), P_PTS(n), &
                                    drho_dt, drho_ds)
            rho_p = eos_density_point(eos, T_PTS(n) + EPS_T, S_PTS(n), P_PTS(n))
            rho_m = eos_density_point(eos, T_PTS(n) - EPS_T, S_PTS(n), P_PTS(n))
            fd_dt = (rho_p - rho_m)/(2.0_wp*EPS_T)
            rho_p = eos_density_point(eos, T_PTS(n), S_PTS(n) + EPS_S, P_PTS(n))
            rho_m = eos_density_point(eos, T_PTS(n), S_PTS(n) - EPS_S, P_PTS(n))
            fd_ds = (rho_p - rho_m)/(2.0_wp*EPS_S)
            call check(error, abs(drho_dt - fd_dt) <= &
                       max(rtol*abs(fd_dt), FD_ATOL), &
                       label//": d(rho)/dT disagrees with the finite difference")
            if (allocated(error)) exit checks
            call check(error, abs(drho_ds - fd_ds) <= &
                       max(rtol*abs(fd_ds), FD_ATOL), &
                       label//": d(rho)/dS disagrees with the finite difference")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine check_fd_for

   subroutine test_fd_wright(error)
      type(error_type), allocatable, intent(out) :: error
      call check_fd_for(error, make_eos(EOS_VARIANT_WRIGHT_97), FD_RTOL, "wright")
   end subroutine test_fd_wright

   subroutine test_fd_roquet(error)
      !! Roquet carries the extra PT->CT conversion polynomial in the
      !! chain rule, so it gets the same bound rather than a looser one —
      !! if the CT-via-SR term were dropped the error would be ~5e-3, far
      !! outside this.
      type(error_type), allocatable, intent(out) :: error
      call check_fd_for(error, make_eos(EOS_VARIANT_ROQUET_SPV), FD_RTOL, "roquet")
   end subroutine test_fd_roquet

   subroutine test_fd_linear(error)
      !! The linear branch ignores `p` and is constant in (T, S), so the
      !! centred difference is EXACT up to the subtraction's round-off —
      !! this is the case that would catch a sign slip between the
      !! `−alpha_T·(T − T_ref)` in the density and the `−∂ρ/∂T` here.
      type(error_type), allocatable, intent(out) :: error
      call check_fd_for(error, make_eos(EOS_VARIANT_LINEAR), 1.0e-10_wp, "linear")
   end subroutine test_fd_linear

   ! -----------------------------------------------------------------
   ! (c) Exactness of the linear branch + the sign twin
   ! -----------------------------------------------------------------

   subroutine test_linear_bitwise(error)
      !! The structural claim behind the knob's bit-identity: under
      !! `eos = "linear"` the EOS-derived pair IS the handle pair, to the
      !! last bit, at every (T, S, p) — including a non-zero pressure,
      !! which the linear branch must ignore entirely.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: a, b
      integer :: n
      eos = make_eos(EOS_VARIANT_LINEAR)
      checks: block
         do n = 1, N_PTS
            call eos_buoyancy_coeffs(eos, T_PTS(n), S_PTS(n), P_PTS(n), a, b)
            call check(error, a == eos%alpha_T, &
                       "linear alpha_T is not bitwise the handle member")
            if (allocated(error)) exit checks
            call check(error, b == eos%beta_S, &
                       "linear beta_S is not bitwise the handle member")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_linear_bitwise

   subroutine test_sign_twin(error)
      !! `eos_density_derivs` must be exactly `(−alpha, +beta)` — a
      !! negation, never a second evaluation of the polynomial (which
      !! would differ in the last bits and silently cost twice).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: a, b, dt, ds
      integer :: n, v
      integer, parameter :: VARIANTS(3) = [EOS_VARIANT_LINEAR, &
                                           EOS_VARIANT_WRIGHT_97, &
                                           EOS_VARIANT_ROQUET_SPV]
      checks: block
         do v = 1, 3
            eos = make_eos(VARIANTS(v))
            do n = 1, N_PTS
               call eos_buoyancy_coeffs(eos, T_PTS(n), S_PTS(n), P_PTS(n), a, b)
               call eos_density_derivs(eos, T_PTS(n), S_PTS(n), P_PTS(n), dt, ds)
               call check(error, dt == -a, "d(rho)/dT is not bitwise -alpha_T")
               if (allocated(error)) exit checks
               call check(error, ds == b, "d(rho)/dS is not bitwise +beta_S")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
   end subroutine test_sign_twin

   ! -----------------------------------------------------------------
   ! (b) Physical trends under Wright
   ! -----------------------------------------------------------------

   subroutine test_alpha_cold(error)
      !! Thermal expansion collapses toward the freezing point.  Asserted
      !! as STRICT MONOTONICITY of α in T along −2 -> 20 degC at S = 34.5,
      !! surface pressure, plus the headline ratio the capability exists
      !! for: α at the freezing point is several times smaller than at
      !! 10 degC, so a 10 degC-calibrated constant over-sizes the
      !! melt-driven buoyancy flux in a cavity by that factor.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      integer, parameter :: N_T = 23
      real(wp) :: a(N_T), b, t
      integer :: n
      eos = make_eos(EOS_VARIANT_WRIGHT_97)
      checks: block
         do n = 1, N_T
            t = -2.0_wp + real(n - 1, wp)
            call eos_buoyancy_coeffs(eos, t, 34.5_wp, 0.0_wp, a(n), b)
         end do
         do n = 2, N_T
            call check(error, a(n) > a(n - 1), &
                       "Wright alpha_T is not strictly increasing with temperature")
            if (allocated(error)) exit checks
         end do
         ! a(1) is T = -2 degC, a(13) is T = +10 degC.
         call check(error, a(13) > 3.0_wp*a(1), &
                    "Wright alpha_T at 10 degC is not several times the freezing-point value")
         if (allocated(error)) exit checks
         call check(error, a(1) > 0.0_wp, &
                    "Wright alpha_T at the freezing point came out non-positive")
      end block checks
   end subroutine test_alpha_cold

   subroutine test_alpha_pressure(error)
      !! Thermobaricity: at a fixed cold state, α must GROW with
      !! pressure, and appreciably so over a cavity's depth range — this
      !! is the second half of why a single scalar α cannot stand in
      !! under a nonlinear EOS.  1000 dbar = 1e7 Pa.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      integer, parameter :: N_P = 11
      real(wp) :: a(N_P), b, p
      integer :: n
      eos = make_eos(EOS_VARIANT_WRIGHT_97)
      checks: block
         do n = 1, N_P
            p = real(n - 1, wp)*1.0e6_wp
            call eos_buoyancy_coeffs(eos, -1.9_wp, 34.5_wp, p, a(n), b)
         end do
         do n = 2, N_P
            call check(error, a(n) > a(n - 1), &
                       "Wright alpha_T is not strictly increasing with pressure")
            if (allocated(error)) exit checks
         end do
         call check(error, a(N_P) > 1.5_wp*a(1), &
                    "Wright alpha_T barely moved over 1000 dbar (thermobaricity missing)")
      end block checks
   end subroutine test_alpha_pressure

   subroutine test_parse(error)
      !! The enum round-trips, and anything else is INVALID rather than a
      !! silent fallback to the constants (`configure_ocean_vmix` turns
      !! the invalid tag into a fail-loud abort).
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_buoyancy_coeffs("constant") == BUOY_COEFFS_CONSTANT, &
                    "'constant' did not parse to BUOY_COEFFS_CONSTANT")
         if (allocated(error)) exit checks
         call check(error, parse_buoyancy_coeffs("eos") == BUOY_COEFFS_EOS, &
                    "'eos' did not parse to BUOY_COEFFS_EOS")
         if (allocated(error)) exit checks
         call check(error, parse_buoyancy_coeffs("EOS") == BUOY_COEFFS_INVALID, &
                    "a mis-cased spelling silently parsed to a valid tag")
         if (allocated(error)) exit checks
         call check(error, parse_buoyancy_coeffs("active") == BUOY_COEFFS_INVALID, &
                    "an unknown spelling did not come back INVALID")
      end block checks
   end subroutine test_parse

   ! -----------------------------------------------------------------
   ! (d) The routed consumers — production kernels, device path
   ! -----------------------------------------------------------------

   subroutine build_column(ms, h_layer, temp, salt, rho)
      !! An unstratified (in `rho_layer`) column at rest, carrying a
      !! uniform (T, S).  `rho_layer` is set INDEPENDENTLY of the tracers
      !! on purpose: KPP's bulk-Ri sweep reads `rho_layer` while the E4
      !! α/β read the tracers, so decoupling them lets a case fix the
      !! boundary-layer depth analytically while varying the state the
      !! EOS sees.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: h_layer, temp, salt, rho
      ms%h_layer = h_layer
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      ms%rho_layer = rho
      ms%tracers(ms%idx_temperature)%hTr = temp*h_layer
      ms%tracers(ms%idx_salinity)%hTr = salt*h_layer
      ms%p_top = 0.0_wp
   end subroutine build_column

   subroutine run_kpp_case(grid, ms, vmix, ss, sf)
      !! Drive the PRODUCTION PP81 + KPP overlay under the `mem:separate`
      !! contract: every state object AND its scratch companion mapped
      !! before the kernels, the results pulled back COMPONENT-wise
      !! (never the aggregate derived type — an aggregate D->H copy
      !! overwrites the host descriptors with device addresses).
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
      ! `p_top` and the tracer `hTr` are host-set inputs mapped by
      ! `enter_data`; push them so the device reads the values this test
      ! set, not whatever the map left behind.
      !$acc update device(ms%p_top, ms%h_layer, ms%rho_layer)
      !$acc update device(ms%tracers(ms%idx_temperature)%hTr)
      !$acc update device(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update device(sf%Q_heat, sf%Q_salt)
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth, vmix%gamma_t, vmix%gamma_s)
      call sf%exit_data()
      call ss%exit_data()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
   end subroutine run_kpp_case

   subroutine test_kpp_linear(error)
      !! Under `eos = "linear"` the two settings must produce BYTE-EQUAL
      !! kv, kt, bl_depth and the non-local gammas — the structural
      !! bit-identity claim, driven all the way through the real kernel
      !! on the device path rather than argued at the point routine.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp), parameter :: H_LAYER = 20.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_c, vmix_e
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer :: n_diff
      checks: block
         call grid%init(10, 8, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix_c%init(grid, nz_ml=NZ)
         call vmix_e%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call sf%init(grid)
         call ss%set_wind_stress_const(0.08_wp, 0.0_wp)
         call build_column(ms, H_LAYER, -1.9_wp, 34.5_wp, 1027.8_wp)
         sf%Q_heat = -150.0_wp    ! cooling  => destabilizing B_0
         sf%Q_salt = 2.0e-5_wp    ! salting  => also destabilizing
         vmix_c%eos = make_eos(EOS_VARIANT_LINEAR)
         vmix_e%eos = vmix_c%eos
         vmix_c%rho0 = RHO0; vmix_e%rho0 = RHO0
         vmix_c%buoyancy_coeffs = BUOY_COEFFS_CONSTANT
         vmix_e%buoyancy_coeffs = BUOY_COEFFS_EOS

         call run_kpp_case(grid, ms, vmix_c, ss, sf)
         call run_kpp_case(grid, ms, vmix_e, ss, sf)

         n_diff = count(vmix_c%kv /= vmix_e%kv) + count(vmix_c%kt /= vmix_e%kt) &
                  + count(vmix_c%bl_depth /= vmix_e%bl_depth) &
                  + count(vmix_c%gamma_t /= vmix_e%gamma_t) &
                  + count(vmix_c%gamma_s /= vmix_e%gamma_s)
         call check(error, n_diff == 0, &
                    "linear EOS: buoyancy_coeffs='eos' is not bitwise equal to 'constant'")
      end block checks
      call sf%destroy(); call ss%destroy()
      call vmix_e%destroy(); call vmix_c%destroy(); call ms%destroy()
   end subroutine test_kpp_linear

   subroutine test_kpp_wright(error)
      !! The quantitative Wright case, at the cold corner the capability
      !! exists for.
      !!
      !! Configuration chosen so the answer is CLOSED FORM: zero wind
      !! (`u_* = 0`), a uniform `rho_layer` (no bulk-Ri crossing, so
      !! `h_b` is the full column depth in BOTH runs), net cooling and
      !! salting (`B_0 < 0`, destabilizing).  KPP then reduces to
      !!
      !!     w_s = w_* = (−B_0·h_b)^(1/3),   kv = h_b·w_s·G(σ)
      !!
      !! with `B_0 = g·(α·q_T − β·q_S)`, so every interface's kv scales
      !! as the CUBE ROOT of the B_0 ratio between the two settings.
      !! That ratio is computed HOST-SIDE from `eos_buoyancy_coeffs` at
      !! the column's own (−1.9 degC, 34.5 PSU, 0 Pa) and from the handle
      !! constants, and must show up in the kernel's output — which pins
      !! the routing, the α value and the device data motion at once.
      !!
      !! The surface layer is kept thick enough that the shape-function
      !! kv is orders above the PP81 background, so the `max()` in the
      !! overlay never truncates the comparison.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp), parameter :: H_LAYER = 40.0_wp
      real(wp), parameter :: T_COL = -1.9_wp, S_COL = 34.5_wp
      real(wp), parameter :: Q_HEAT = -200.0_wp, Q_SALT = 3.0e-5_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_c, vmix_e
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      real(wp) :: a_eos, b_eos, q_t_kin, q_s_kin, b0_c, b0_e
      real(wp) :: ratio_expect, ratio_obs, kv_c, kv_e
      integer :: ip, jp, k
      checks: block
         call grid%init(10, 8, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix_c%init(grid, nz_ml=NZ)
         call vmix_e%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call sf%init(grid)
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call build_column(ms, H_LAYER, T_COL, S_COL, 1027.8_wp)
         sf%Q_heat = Q_HEAT
         sf%Q_salt = Q_SALT
         vmix_c%eos = make_eos(EOS_VARIANT_WRIGHT_97)
         vmix_e%eos = vmix_c%eos
         vmix_c%rho0 = RHO0; vmix_e%rho0 = RHO0
         vmix_c%buoyancy_coeffs = BUOY_COEFFS_CONSTANT
         vmix_e%buoyancy_coeffs = BUOY_COEFFS_EOS

         call run_kpp_case(grid, ms, vmix_c, ss, sf)
         call run_kpp_case(grid, ms, vmix_e, ss, sf)

         ! Host-side oracle for the B_0 ratio.
         call eos_buoyancy_coeffs(vmix_e%eos, T_COL, S_COL, 0.0_wp, a_eos, b_eos)
         q_t_kin = Q_HEAT/(RHO0*sf%cp)
         q_s_kin = Q_SALT/RHO0
         b0_c = GRAVITY*(ALPHA_LIN*q_t_kin - BETA_LIN*q_s_kin)
         b0_e = GRAVITY*(a_eos*q_t_kin - b_eos*q_s_kin)
         call check(error, b0_c < 0.0_wp .and. b0_e < 0.0_wp, &
                    "test setup: B_0 is not destabilizing in both settings")
         if (allocated(error)) exit checks
         ! The EOS alpha at the freezing point is well below the 10 degC
         ! constant, so the EOS run must be the WEAKER of the two.
         call check(error, a_eos < 0.5_wp*ALPHA_LIN, &
                    "test setup: Wright alpha at -1.9 degC is not far below the constant")
         if (allocated(error)) exit checks
         ratio_expect = (b0_e/b0_c)**(1.0_wp/3.0_wp)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ! Every interior interface inside the BL carries the same ratio;
         ! check them all rather than one probe.
         do k = 2, NZ
            kv_c = vmix_c%kv(ip, jp, k)
            kv_e = vmix_e%kv(ip, jp, k)
            call check(error, kv_c > 1.0e-3_wp, &
                       "test setup: the KPP overlay did not dominate the PP81 background")
            if (allocated(error)) exit checks
            ratio_obs = kv_e/kv_c
            call check(error, abs(ratio_obs - ratio_expect) < 1.0e-10_wp*ratio_expect, &
                       "KPP kv did not scale as the cube root of the alpha-driven B_0 ratio")
            if (allocated(error)) exit checks
         end do
         ! And the two runs really did differ — a routing that silently
         ! fell back to the constants would pass every ratio above only
         ! if ratio_expect were 1, which this guard forbids.
         call check(error, abs(ratio_expect - 1.0_wp) > 0.1_wp, &
                    "test setup: the two settings are indistinguishable")
      end block checks
      call sf%destroy(); call ss%destroy()
      call vmix_e%destroy(); call vmix_c%destroy(); call ms%destroy()
   end subroutine test_kpp_wright

   subroutine run_ddiff_case(grid, ms, vmix)
      !! Drive `vmix_split_kd_heat_salt` (the double-diffusion arm) on the
      !! device path.  `kt` is a host-set input here — the production
      !! caller has every contributor write it first — so it is pushed
      !! explicitly after the map.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      !$acc update device(ms%p_top, ms%h_layer)
      !$acc update device(ms%tracers(ms%idx_temperature)%hTr)
      !$acc update device(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update device(vmix%kt)
      call vmix_split_kd_heat_salt(grid, vmix, ms)
      !$acc update self(vmix%kt, vmix%ks)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_ddiff_case

   subroutine setup_ddiff(ms, vmix, nz, h_layer, variant)
      !! A salt-fingering column: warm+salty over cold+fresh, so
      !! `alpha*dT >= beta*dS > 0` at every interior interface and the
      !! fingering arm of the closure is the one under test.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer
      integer, intent(in) :: variant
      real(wp) :: t_k, s_k
      integer :: k
      ms%h_layer = h_layer
      ms%rho_layer = 1027.0_wp
      ms%p_top = 0.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ! k = 1 is the bed, k = nz the surface: warm/salty at the top.
         t_k = -1.0_wp + 0.45_wp*real(k - 1, wp)
         s_k = 34.0_wp + 0.06_wp*real(k - 1, wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*h_layer
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = s_k*h_layer
      end do
      vmix%eos = make_eos(variant)
      vmix%rho0 = RHO0
      vmix%ddiff_enable = .true.
      vmix%kt = 1.0e-5_wp
   end subroutine setup_ddiff

   subroutine test_ddiff_linear(error)
      !! Double diffusion under `eos = "linear"`: byte-equal kt and ks
      !! between the two settings.  This also exercises the NEW
      !! column-structured `vmix_split_ddiff_eos_impl` against the
      !! untouched collapsed one, so a divergence in the closure's
      !! branch structure (not just in α/β) would show up here.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp), parameter :: H_LAYER = 25.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_c, vmix_e
      integer :: n_diff
      checks: block
         call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix_c%init(grid, nz_ml=NZ)
         call vmix_e%init(grid, nz_ml=NZ)
         call setup_ddiff(ms, vmix_c, NZ, H_LAYER, EOS_VARIANT_LINEAR)
         call setup_ddiff(ms, vmix_e, NZ, H_LAYER, EOS_VARIANT_LINEAR)
         vmix_c%buoyancy_coeffs = BUOY_COEFFS_CONSTANT
         vmix_e%buoyancy_coeffs = BUOY_COEFFS_EOS

         call run_ddiff_case(grid, ms, vmix_c)
         call run_ddiff_case(grid, ms, vmix_e)

         ! The closure must actually have fired, or "identical" is vacuous.
         call check(error, any(vmix_c%ks > 1.0e-5_wp), &
                    "test setup: double diffusion produced no salt-fingering Kd")
         if (allocated(error)) exit checks
         n_diff = count(vmix_c%kt /= vmix_e%kt) + count(vmix_c%ks /= vmix_e%ks)
         call check(error, n_diff == 0, &
                    "linear EOS: ddiff under buoyancy_coeffs='eos' is not bitwise "// &
                    "equal to 'constant'")
      end block checks
      call vmix_e%destroy(); call vmix_c%destroy(); call ms%destroy()
   end subroutine test_ddiff_linear

   subroutine test_ddiff_wright(error)
      !! Under Wright at the cold corner the per-interface α is far below
      !! the 10 degC constant while β moves much less, so the density
      !! ratio `R_rho = α·ΔT / β·ΔS` — the quantity the whole closure
      !! branches on — genuinely moves, and the resulting Kd must differ.
      !! Asserted as "the split changed", not as a number: the closed
      !! form here is the CVMix clamp, which is not worth re-deriving in
      !! a test when the FD suite already pins α itself.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp), parameter :: H_LAYER = 25.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_c, vmix_e
      checks: block
         call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix_c%init(grid, nz_ml=NZ)
         call vmix_e%init(grid, nz_ml=NZ)
         call setup_ddiff(ms, vmix_c, NZ, H_LAYER, EOS_VARIANT_WRIGHT_97)
         call setup_ddiff(ms, vmix_e, NZ, H_LAYER, EOS_VARIANT_WRIGHT_97)
         vmix_c%buoyancy_coeffs = BUOY_COEFFS_CONSTANT
         vmix_e%buoyancy_coeffs = BUOY_COEFFS_EOS

         call run_ddiff_case(grid, ms, vmix_c)
         call run_ddiff_case(grid, ms, vmix_e)

         call check(error, any(vmix_c%ks > 1.0e-5_wp) .or. any(vmix_e%ks > 1.0e-5_wp), &
                    "test setup: double diffusion produced no Kd in either setting")
         if (allocated(error)) exit checks
         call check(error, any(abs(vmix_c%ks - vmix_e%ks) > 1.0e-12_wp), &
                    "Wright EOS: the ddiff split did not respond to the EOS-derived alpha/beta")
      end block checks
      call vmix_e%destroy(); call vmix_c%destroy(); call ms%destroy()
   end subroutine test_ddiff_wright

end module test_ocean_eos_buoyancy
