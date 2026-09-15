!! Analytical unit tests for the shortwave / boundary-layer coupling
!! (PR-21): the KPP `KPP_SHORTWAVE_METHOD` `B_0` correction in
!! `rdb_ocean_vmix` and the EPBL penetrating-SW TKE ledger in
!! `rdb_ocean_epbl`, both computed from the shared two-band
!! `sw_transmission` / `sw_pe_cost_shape` in `rdb_ocean_surface_flux`.
!!
!! Physics (Paulson & Simpson 1977 two-band; Large, McWilliams & Doney
!! 1994 App. B buoyancy decomposition; Reichl & Hallberg 2018 EPBL TKE):
!! only the shortwave ABSORBED INSIDE the boundary layer stabilises it;
!! the fraction that leaks below cannot.  Charging the BL for the full
!! surface deposit spuriously suppresses mixing on sunny days.
!!
!! Bottom-up convention: surface = `k = nz`, bed = `k = 1`.
!!
!! Cases:
!!   1 sw_off_bitident            — the new knobs are inert at sw off (==)
!!   2 phi_recovers_skin_limit    — Phi(inf)=1, Phi(0)=0, Taylor seam
!!   3 epbl_ctke_ledger_partitions_i0 — ledger heat == deposition heat
!!   4 epbl_midday_mld_deepens    — telling EPBL the truth deepens the ML
!!                                  (NEGATIVE CONTROL: neutralise -> fails)
!!   5 kpp_mxl_sw_deepens_under_sunny_cooling — B_0 responds only here
!!   6 sw_bl_gpu_resident         — mem:separate device-residency (CLAUDE.md)
module test_ocean_sw_bl_coupling
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     sw_transmission, sw_pe_cost_shape, SEAWATER_CP
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, &
                             vmix_apply_kpp_overlay, &
                             KPP_SW_ALL, KPP_SW_MXL
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_compute, EPBL_MSTAR_CONSTANT
   implicit none
   private

   public :: collect_ocean_sw_bl_coupling_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.2_wp     ! kg/m^3/degC
   real(wp), parameter :: BETA_LIN = 0.8_wp      ! kg/m^3/PSU
   real(wp), parameter :: R = 0.58_wp
   real(wp), parameter :: ZETA1 = 0.35_wp
   real(wp), parameter :: ZETA2 = 23.0_wp

contains

   subroutine collect_ocean_sw_bl_coupling_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sw_off_bitident", test_sw_off_bitident), &
                  new_unittest("phi_recovers_skin_limit", test_phi_limits), &
                  new_unittest("epbl_ctke_ledger_partitions_i0", test_ledger_partitions), &
                  new_unittest("epbl_midday_mld_deepens", test_mld_deepens), &
                  new_unittest("kpp_mxl_sw_deepens_under_sunny_cooling", test_kpp_mxl_sw), &
                  new_unittest("sw_bl_gpu_resident", test_gpu_resident) &
                  ]
   end subroutine collect_ocean_sw_bl_coupling_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine setup_epbl(epbl, grid, nz)
      type(ocean_epbl_t), intent(inout) :: epbl
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call epbl%init(grid, nz_ml=nz)
      epbl%enable = .true.
      epbl%mstar_scheme = EPBL_MSTAR_CONSTANT
      epbl%mstar_const = 1.2_wp
      epbl%nstar = 0.2_wp
      epbl%eos%variant = EOS_VARIANT_LINEAR
      epbl%eos%alpha_T = ALPHA_LIN
      epbl%eos%beta_S = BETA_LIN
      epbl%eos%rho0 = RHO0
      epbl%rho0 = RHO0
      epbl%f_centre = 0.0_wp
      epbl%tke_diags = .true.
   end subroutine setup_epbl

   !> Fill a uniform-T/S column (no stratification) of NZ layers.
   subroutine set_uniform_column(ms, nz, dz, t_val, s_val)
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, t_val, s_val
      integer :: k
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_val*dz
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = s_val*dz
      end do
   end subroutine set_uniform_column

   ! -----------------------------------------------------------------
   ! Test 1 — the new knobs are inert when SW is off (bit-identity).
   ! -----------------------------------------------------------------

   subroutine test_sw_off_bitident(error)
      !! With `sw_pen_frac = 0` (⇒ `has_sw = .false.`) the SW terms are
      !! gated OUT host-side, so the KPP method and the EPBL ledger knob
      !! must make NO difference: a run with the "physics-parity" defaults
      !! (mxl_sw, epbl_sw_ctke=.true.) is byte-for-byte equal (`==`) to a
      !! run with the legacy knobs (all_sw, epbl_sw_ctke=.false.).  This
      !! is the assertion that catches relying on `x - 0.0*y == x` instead
      !! of gating on `has_sw`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(ocean_epbl_t) :: epbl_a, epbl_b
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      logical :: ok
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_strat_ms(ms, NZ)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.05_wp, 0.0_wp)
         call sf%init(grid)
         ! SW OFF (sw_pen_frac = 0) — set a non-zero heat flux so the
         ! else-branch arithmetic is exercised, but has_sw stays false.
         call sf%set_surface_flux_const(150.0_wp, 0.0_wp)
         call sf%set_sw_penetration(0.0_wp, R, ZETA1, ZETA2)

         ! Run A — parity defaults (ms unchanged by the read-only kernels).
         call vmix_a%init(grid, nz_ml=NZ)
         call setup_eos_vmix(vmix_a)
         vmix_a%kpp_sw_method = KPP_SW_MXL
         call setup_epbl(epbl_a, grid, NZ)
         epbl_a%epbl_sw_ctke = .true.
         call run_kpp(grid, ms, vmix_a, ss, sf)
         call run_epbl(grid, ms, epbl_a, ss, sf)

         ! Run B — legacy knobs.
         call vmix_b%init(grid, nz_ml=NZ)
         call setup_eos_vmix(vmix_b)
         vmix_b%kpp_sw_method = KPP_SW_ALL
         call setup_epbl(epbl_b, grid, NZ)
         epbl_b%epbl_sw_ctke = .false.
         call run_kpp(grid, ms, vmix_b, ss, sf)
         call run_epbl(grid, ms, epbl_b, ss, sf)

         ok = all(vmix_a%kv == vmix_b%kv) .and. all(vmix_a%kt == vmix_b%kt) .and. &
              all(vmix_a%gamma_t == vmix_b%gamma_t) .and. &
              all(vmix_a%bl_depth == vmix_b%bl_depth)
         call check(error, ok, "KPP kv/kt/gamma_t/bl_depth must be == with SW off")
         if (allocated(error)) exit checks
         ok = all(epbl_a%kd_int == epbl_b%kd_int) .and. all(epbl_a%mld == epbl_b%mld)
         call check(error, ok, "EPBL kd_int/mld must be == with SW off")
      end block checks
      call cleanup_vmix(vmix_a); call cleanup_vmix(vmix_b)
      call cleanup_epbl(epbl_a); call cleanup_epbl(epbl_b)
      call ss%destroy(); call sf%destroy(); call ms%destroy()
   end subroutine test_sw_off_bitident

   ! -----------------------------------------------------------------
   ! Test 2 — the in-layer PE-cost shape function Phi(tau).
   ! -----------------------------------------------------------------

   subroutine test_phi_limits(error)
      !! `Phi(inf) = 1` (the shape function reduces to the existing skin
      !! `ctke_sfc` — the correctness anchor for the EPBL half); the
      !! approach is the analytic `Phi(tau) ~ 1 - 2/tau` (a thick layer
      !! relative to the absorption depth deposits all the heat near its
      !! top ⇒ full skin cost), so `Phi(50) = 0.96 < 1` and only very
      !! large `tau` reaches 1.  `Phi(0) = 0` (a uniformly-heated layer
      !! costs nothing to homogenise); the closed form and the Taylor
      !! branch agree across the `tau = 1e-2` seam (catches a wrong Taylor
      !! coefficient); `Phi` monotone increasing on `[1e-4, 50]`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tau, phi_prev, phi_now, closed, em1, phi_seam_taylor, phi_seam_closed
      integer :: n
      logical :: monotone
      checks: block
         call check(error, sw_pe_cost_shape(1.0e7_wp) > 1.0_wp - 1.0e-6_wp, &
                    "Phi(inf) must -> 1 (recovers the skin limit)")
         if (allocated(error)) exit checks
         ! The exact asymptotic Phi ~ 1 - 2/tau is the shape's fingerprint.
         call check(error, abs(sw_pe_cost_shape(1.0e4_wp) - (1.0_wp - 2.0e-4_wp)) < 1.0e-7_wp, &
                    "Phi must approach 1 as 1 - 2/tau")
         if (allocated(error)) exit checks
         call check(error, sw_pe_cost_shape(50.0_wp) < 1.0_wp, &
                    "Phi(tau) must be strictly < 1 for finite tau")
         if (allocated(error)) exit checks
         call check(error, sw_pe_cost_shape(1.0e-6_wp) < 1.0e-6_wp, &
                    "Phi(0) must -> 0")
         if (allocated(error)) exit checks
         ! Seam: at tau = 1e-2 the function uses the Taylor branch (which
         ! is accurate to ~1e-15); the naive closed form suffers ~1e-12 of
         ! catastrophic cancellation there (exactly why the Taylor branch
         ! exists).  They must still agree to <= 1e-10 — a tolerance well
         ! below the ~3e-9 shift a wrong C1_6 / C1_60 coefficient would
         ! introduce, so a mis-typed Taylor coefficient still fails here.
         tau = 1.0e-2_wp
         phi_seam_taylor = sw_pe_cost_shape(tau)
         em1 = 1.0_wp - exp(-tau)
         phi_seam_closed = (tau*(1.0_wp + exp(-tau)) - 2.0_wp*em1)/(tau*em1)
         call check(error, abs(phi_seam_taylor - phi_seam_closed) < 1.0e-10_wp, &
                    "Taylor and closed-form Phi disagree across the tau=1e-2 seam")
         if (allocated(error)) exit checks
         ! Monotone increasing on a log-spaced grid over [1e-4, 50].
         monotone = .true.
         phi_prev = sw_pe_cost_shape(1.0e-4_wp)
         do n = 1, 200
            tau = 1.0e-4_wp*(50.0_wp/1.0e-4_wp)**(real(n, wp)/200.0_wp)
            phi_now = sw_pe_cost_shape(tau)
            if (.not. (phi_now > phi_prev)) monotone = .false.
            phi_prev = phi_now
         end do
         call check(error, monotone, "Phi must be monotone increasing on [1e-4, 50]")
         if (allocated(error)) exit checks
         ! Bounded in (0, 1) for a mid-range tau.
         closed = sw_pe_cost_shape(2.0_wp)
         call check(error, closed > 0.0_wp .and. closed < 1.0_wp, &
                    "Phi(2) must lie strictly in (0, 1)")
      end block checks
   end subroutine test_phi_limits

   ! -----------------------------------------------------------------
   ! Test 3 — the EPBL ledger and the deposition kernel account the
   ! same I0 (single band, uniform linear-EOS column: the prefactor is
   ! constant across k, so Sum ctke_sw recovers Sum heat analytically).
   ! -----------------------------------------------------------------

   subroutine test_ledger_partitions(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 20
      real(wp), parameter :: DZ = 5.0_wp
      real(wp), parameter :: DT = 1800.0_wp
      real(wp), parameter :: QSW = 300.0_wp   ! W/m^2 into the ocean
      real(wp), parameter :: FRAC = 1.0_wp
      real(wp) :: ctke_sum, prefac, sum_heat, i0, expected_heat, phi1
      real(wp) :: dcolht_expect, dcolht_got
      integer :: i, j, k
      checks: block
         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.05_wp, 0.0_wp)
         call set_uniform_column(ms, NZ, DZ, 12.0_wp, 35.0_wp)
         call sf%init(grid)
         call sf%set_surface_flux_const(QSW, 0.0_wp)
         ! Single band (R = 1 ⇒ only zeta1) so the two-band Phi weighting
         ! collapses to one factor; net-heat source (I0 = FRAC*Q_heat > 0).
         call sf%set_sw_penetration(FRAC, 1.0_wp, ZETA1, ZETA2)
         call setup_epbl(epbl, grid, NZ)
         epbl%epbl_sw_ctke = .true.

         call run_epbl(grid, ms, epbl, ss, sf)

         i = grid%nx_total/2
         j = grid%ny_total/2

         ! (a) The dcolht_t identity rho0^2 h dsv_dt == rho0 dcolht_t.
         !     Linear EOS: dsv_dt = alpha_T/rho0^2, so dcolht_t = h*alpha/rho0.
         dcolht_expect = DZ*ALPHA_LIN/RHO0
         dcolht_got = epbl%dcolht_t%data(i, j, NZ)
         call check(error, abs(dcolht_got - dcolht_expect) < 1.0e-12_wp*dcolht_expect, &
                    "dcolht_t identity broken (linear EOS)")
         if (allocated(error)) exit checks

         ! (b) The prefactor -0.5 g rho0 dcolht_t Phi(h/zeta1) is constant
         !     across k (uniform column), so Sum_k ctke_sw = prefac*Sum heat.
         phi1 = sw_pe_cost_shape(DZ/ZETA1)
         prefac = -0.5_wp*GRAVITY*RHO0*dcolht_expect*phi1
         ctke_sum = 0.0_wp
         do k = 1, NZ
            ctke_sum = ctke_sum + epbl%ctke_sw%data(i, j, k)
         end do
         sum_heat = ctke_sum/prefac
         i0 = FRAC*QSW
         expected_heat = i0*DT/(RHO0*SEAWATER_CP)
         call check(error, abs(sum_heat - expected_heat) < 1.0e-9_wp*expected_heat, &
                    "EPBL ledger heat must equal the deposition kernel's I0 dt/(rho0 cp)")
         if (allocated(error)) exit checks

         ! (c) NEGATIVE-CONTROL sanity: the term is actually non-zero.
         call check(error, abs(ctke_sum) > 0.0_wp, &
                    "ctke_sw column must be non-zero when SW is active")
      end block checks
      call cleanup_epbl(epbl)
      call ss%destroy(); call sf%destroy(); call ms%destroy()
   end subroutine test_ledger_partitions

   ! -----------------------------------------------------------------
   ! Test 4 — telling EPBL the truth about where the sunlight went makes
   ! the midday mixed layer DEEPER (the all-skin ledger over-charges and
   ! spuriously suppresses mixing; MOM6's own MXL_SW rationale).
   ! This is the NEGATIVE CONTROL: epbl_sw_ctke=.false. neutralises the
   ! new term; if the ledger contribution were (wrongly) zero in BOTH
   ! runs, mld_true == mld_false and this assertion FAILS.
   ! -----------------------------------------------------------------

   subroutine test_mld_deepens(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl_on, epbl_off
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 30
      real(wp) :: mld_on, mld_off
      integer :: i, j
      checks: block
         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.08_wp, 0.0_wp)   ! moderate wind
         call sf%init(grid)
         call sf%set_surface_flux_const(300.0_wp, 0.0_wp)  ! strong daytime heating
         call sf%set_sw_penetration(1.0_wp, R, ZETA1, ZETA2)

         ! EPBL ledger ON (the truth: SW spread through the column).
         call setup_strat_ms(ms, NZ)
         call setup_epbl(epbl_on, grid, NZ)
         epbl_on%epbl_sw_ctke = .true.
         call run_epbl(grid, ms, epbl_on, ss, sf)

         ! EPBL ledger OFF (the all-skin control).
         call setup_strat_ms(ms, NZ)
         call setup_epbl(epbl_off, grid, NZ)
         epbl_off%epbl_sw_ctke = .false.
         call run_epbl(grid, ms, epbl_off, ss, sf)

         i = grid%nx_total/2
         j = grid%ny_total/2
         mld_on = epbl_on%mld(i, j)
         mld_off = epbl_off%mld(i, j)
         call check(error, mld_on > mld_off + 1.0e-6_wp, &
                    "midday MLD must DEEPEN when the SW ledger is charged truthfully")
      end block checks
      call cleanup_epbl(epbl_on); call cleanup_epbl(epbl_off)
      call ss%destroy(); call sf%destroy(); call ms%destroy()
   end subroutine test_mld_deepens

   ! -----------------------------------------------------------------
   ! Test 5 — KPP B_0 responds only in the sunny-but-net-cooling regime
   ! (Q_heat < 0, q_sw > 0).  With u_star = 0 and c_vt2 = 0 the BL depth
   ! is B_0-independent, so kv ratio at any BL interface is exactly
   ! (w*_mxl / w*_all) = ((-B_0_mxl)/(-B_0_all))^(1/3), pinned against an
   ! independent sw_transmission evaluation.  Under heating (B_0 > 0)
   ! w* = 0, so all_sw and mxl_sw are bit-for-bit identical.
   ! -----------------------------------------------------------------

   subroutine test_kpp_mxl_sw(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_mxl, vmix_all
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 10
      real(wp), parameter :: QCOOL = -250.0_wp   ! net cooling
      real(wp), parameter :: QSW = 400.0_wp      ! sunlight in
      real(wp) :: h_b, i0, t_hb, q_bl_mxl, b0_mxl, b0_all
      real(wp) :: ratio_pred, ratio_obs, kv_mxl, kv_all
      integer :: i, j, k, kmax
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid)
         ! No wind ⇒ u_star = 0 ⇒ w_s = w_* exactly.
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call sf%init(grid)
         call sf%set_components(grid, .true.)         ! allocate q_sw
         call sf%set_surface_flux_const(QCOOL, 0.0_wp)
         ! q_sw source: I0 = frac*q_sw > 0 while Q_heat < 0 (the regime).
         call sf%set_sw_penetration(1.0_wp, R, ZETA1, ZETA2, sw_source="q_sw")
         sf%q_sw = QSW

         i = grid%nx_total/2
         j = grid%ny_total/2

         ! ---- Sunny cooling: MXL vs ALL ----
         ! Unstratified + no wind: h_b = full depth (bumps everywhere) and
         ! u_star = 0 ⇒ w_s = w_* ⇒ kv ratio is the clean cube-root ratio.
         call setup_unstrat_ms(ms, NZ)
         call vmix_mxl%init(grid, nz_ml=NZ)
         vmix_mxl%c_vt2 = 0.0_wp
         vmix_mxl%kpp_sw_method = KPP_SW_MXL
         call setup_eos_vmix(vmix_mxl)
         call run_kpp(grid, ms, vmix_mxl, ss, sf)

         call setup_unstrat_ms(ms, NZ)
         call vmix_all%init(grid, nz_ml=NZ)
         vmix_all%c_vt2 = 0.0_wp
         vmix_all%kpp_sw_method = KPP_SW_ALL
         call setup_eos_vmix(vmix_all)
         call run_kpp(grid, ms, vmix_all, ss, sf)

         ! Both runs share the same BL depth (c_vt2 = 0 ⇒ B_0-independent).
         h_b = vmix_mxl%bl_depth(i, j)
         call check(error, abs(h_b - vmix_all%bl_depth(i, j)) < 1.0e-12_wp*h_b, &
                    "BL depth must be B_0-independent (c_vt2 = 0)")
         if (allocated(error)) exit checks

         ! Predicted B_0 (q_S = 0): B_0 = g*alpha*Q_bl/(rho0*cp).
         i0 = 1.0_wp*QSW
         t_hb = sw_transmission(h_b, R, ZETA1, ZETA2)
         q_bl_mxl = QCOOL - i0*t_hb
         b0_all = GRAVITY*ALPHA_LIN*QCOOL/(RHO0*SEAWATER_CP)
         b0_mxl = GRAVITY*ALPHA_LIN*q_bl_mxl/(RHO0*SEAWATER_CP)
         call check(error, b0_mxl < b0_all .and. b0_all < 0.0_wp, &
                    "expected B_0(mxl) < B_0(all) < 0 under sunny cooling")
         if (allocated(error)) exit checks

         ! Ratio of w_* (= kv, since h_b, G, u_star all match) is exact.
         ratio_pred = ((-b0_mxl)/(-b0_all))**(1.0_wp/3.0_wp)

         ! Probe the interface with the largest kv (well inside the BL).
         kmax = 0
         kv_mxl = 0.0_wp
         do k = 2, NZ
            if (vmix_mxl%kv(i, j, k) > kv_mxl) then
               kv_mxl = vmix_mxl%kv(i, j, k)
               kmax = k
            end if
         end do
         call check(error, kmax > 0 .and. kv_mxl > 0.0_wp, &
                    "expected a non-zero KPP kv bump under sunny cooling")
         if (allocated(error)) exit checks
         kv_all = vmix_all%kv(i, j, kmax)
         ratio_obs = kv_mxl/kv_all
         call check(error, abs(ratio_obs - ratio_pred) < 1.0e-9_wp*ratio_pred, &
                    "kv(mxl)/kv(all) must equal ((-B0mxl)/(-B0all))^(1/3)")
         if (allocated(error)) exit checks
         call check(error, kv_mxl > kv_all, &
                    "mxl_sw must convect harder than all_sw under sunny cooling")
      end block checks
      call cleanup_vmix(vmix_mxl); call cleanup_vmix(vmix_all)
      call ss%destroy(); call sf%destroy(); call ms%destroy()

      if (allocated(error)) return
      ! ---- Heating no-op: B_0 > 0 ⇒ w_* = 0 ⇒ mxl == all bit-for-bit ----
      call kpp_heating_noop(error)
   end subroutine test_kpp_mxl_sw

   subroutine kpp_heating_noop(error)
      !! Under net heating (Q_heat > 0), B_0 > 0 ⇒ w_* = 0, so the SW
      !! method cannot move kv: mxl_sw and all_sw are bit-identical.
      !! Documents the §3.2 no-op honestly.  Needs wind (u_star > 0) so
      !! kv is non-zero and the equality is informative.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_mxl, vmix_all
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 10
      logical :: ok
      checks: block
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.1_wp, 0.0_wp)
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call sf%set_surface_flux_const(200.0_wp, 0.0_wp)   ! heating
         call sf%set_sw_penetration(1.0_wp, R, ZETA1, ZETA2, sw_source="q_sw")
         sf%q_sw = 300.0_wp

         ! Unstratified + wind: h_b = full depth, u_star > 0 ⇒ kv non-zero;
         ! under heating B_0 > 0 ⇒ w_* = 0 ⇒ mxl_sw and all_sw coincide.
         call setup_unstrat_ms(ms, NZ)
         call vmix_mxl%init(grid, nz_ml=NZ)
         vmix_mxl%kpp_sw_method = KPP_SW_MXL
         call setup_eos_vmix(vmix_mxl)
         call run_kpp(grid, ms, vmix_mxl, ss, sf)

         call setup_unstrat_ms(ms, NZ)
         call vmix_all%init(grid, nz_ml=NZ)
         vmix_all%kpp_sw_method = KPP_SW_ALL
         call setup_eos_vmix(vmix_all)
         call run_kpp(grid, ms, vmix_all, ss, sf)

         ok = all(vmix_mxl%kv == vmix_all%kv) .and. all(vmix_mxl%kt == vmix_all%kt)
         call check(error, ok, "under heating (B_0>0) mxl_sw must equal all_sw (w_*=0)")
      end block checks
      call cleanup_vmix(vmix_mxl); call cleanup_vmix(vmix_all)
      call ss%destroy(); call sf%destroy(); call ms%destroy()
   end subroutine kpp_heating_noop

   ! -----------------------------------------------------------------
   ! Test 6 — mem:separate device residency (CLAUDE.md:312).  A missing
   ! ctke_sw enter_data gives a device-side ctke_sw == 0 with no crash;
   ! this test drives EPBL fully on device and asserts the ledger is
   ! non-zero and matches a host reference.  Follows the
   ! test_open_boundary_out_closes template.
   ! -----------------------------------------------------------------

   subroutine test_gpu_resident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 20
      real(wp) :: ctke_sum, mld_probe
      integer :: i, j, k
      checks: block
         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.06_wp, 0.0_wp)
         call setup_strat_ms(ms, NZ)
         call sf%init(grid)
         call sf%set_surface_flux_const(300.0_wp, 0.0_wp)
         call sf%set_sw_penetration(1.0_wp, R, ZETA1, ZETA2)
         call setup_epbl(epbl, grid, NZ)
         epbl%epbl_sw_ctke = .true.

         !$acc enter data copyin(ms, ss, sf, epbl)
         call ms%enter_data()
         call ss%enter_data()
         call sf%enter_data()
         call epbl%enter_data()
         !$acc update device(sf%Q_heat)
         call epbl_compute(grid, epbl, ms, ss, 1800.0_wp, sf=sf)
         ! Pull COMPONENT arrays (never the whole DT) back to host.
         !$acc update self(epbl%mld, epbl%kd_int)
         !$acc update self(epbl%ctke_sw%data)
         call sf%exit_data()
         call ss%exit_data()
         call epbl%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ss, sf, epbl)

         i = grid%nx_total/2
         j = grid%ny_total/2
         ctke_sum = 0.0_wp
         do k = 1, NZ
            ctke_sum = ctke_sum + epbl%ctke_sw%data(i, j, k)
         end do
         mld_probe = epbl%mld(i, j)
         call check(error, abs(ctke_sum) > 0.0_wp, &
                    "device ctke_sw must be non-zero (missing enter_data ⇒ silent 0)")
         if (allocated(error)) exit checks
         call check(error, mld_probe > 0.0_wp, "device EPBL mld must be positive")
      end block checks
      call cleanup_epbl(epbl)
      call ss%destroy(); call sf%destroy(); call ms%destroy()
   end subroutine test_gpu_resident

   ! -----------------------------------------------------------------
   ! Shared run helpers (device-resident, mem:separate discipline).
   ! -----------------------------------------------------------------

   subroutine setup_strat_ms(ms, nz)
      !! Stably stratified column (surface lightest): rho_layer for KPP,
      !! and a matching T/S profile for EPBL (T decreasing downward).
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      integer :: k
      real(wp) :: t_k, z_below
      ms%h_layer = 5.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%rho_layer(:, :, k) = 1030.0_wp - 0.15_wp*real(k - 1, wp)
         z_below = (real(nz - k, wp) + 0.5_wp)*5.0_wp
         t_k = 18.0_wp - 0.05_wp*z_below
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*5.0_wp
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*5.0_wp
      end do
   end subroutine setup_strat_ms

   subroutine setup_unstrat_ms(ms, nz)
      !! Unstratified column (uniform rho, zero velocities): the bulk-Ri
      !! sweep never crosses ri_crit ⇒ h_b = full column depth ⇒ EVERY
      !! interior interface carries a KPP bump.  Density-independent, so
      !! h_b is identical between the mxl_sw and all_sw runs.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      integer :: k
      ms%h_layer = 5.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%rho_layer(:, :, k) = 1030.0_wp
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 15.0_wp*5.0_wp
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*5.0_wp
      end do
   end subroutine setup_unstrat_ms

   subroutine setup_eos_vmix(vmix)
      type(ocean_vmix_t), intent(inout) :: vmix
      vmix%rho0 = RHO0
      vmix%eos%variant = EOS_VARIANT_LINEAR
      vmix%eos%alpha_T = ALPHA_LIN
      vmix%eos%beta_S = BETA_LIN
      vmix%eos%rho0 = RHO0
      ! Kill the PP81 shear-enhanced part so the KPP bump is the clear
      ! signal in kv (background only underneath).
      vmix%pp81_nu0 = 0.0_wp
   end subroutine setup_eos_vmix

   subroutine run_kpp(grid, ms, vmix, ss, sf)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      !$acc enter data copyin(ms, vmix, ss, sf)
      call ms%enter_data()
      call vmix%enter_data()
      call sf%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth, vmix%gamma_t, vmix%gamma_s)
      call sf%exit_data()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
   end subroutine run_kpp

   subroutine run_epbl(grid, ms, epbl, ss, sf)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      !$acc enter data copyin(ms, epbl, ss, sf)
      call ms%enter_data()
      call sf%enter_data()
      call epbl%enter_data()
      call epbl_compute(grid, epbl, ms, ss, 1800.0_wp, sf=sf)
      !$acc update self(epbl%mld, epbl%kd_int)
      !$acc update self(epbl%ctke_sw%data, epbl%dcolht_t%data)
      call epbl%exit_data()
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, sf)
   end subroutine run_epbl

   subroutine cleanup_vmix(vmix)
      type(ocean_vmix_t), intent(inout) :: vmix
      call vmix%destroy()
   end subroutine cleanup_vmix

   subroutine cleanup_epbl(epbl)
      type(ocean_epbl_t), intent(inout) :: epbl
      call epbl%destroy()
   end subroutine cleanup_epbl

end module test_ocean_sw_bl_coupling
