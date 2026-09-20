!! Phase 4b — THE BOUNDARY-LAYER SCHEMES UNDER AN ICE SHELF.
!!
!! Both "surface" schemes in this model (KPP in `rdb_ocean_vmix`, EPBL in
!! `rdb_ocean_epbl`) take their friction velocity from ONE place,
!! `ocean_surface_stress_t`.  Under an ice shelf the wind is masked out of
!! the `tau` pair at configure, so before this PR `stress_mag` — and
!! therefore `u_*` — was EXACTLY ZERO under every covered column while the
!! ice-ocean stress sat unread in the top-drag slot.  A covered column was
!! being mixed as if nothing touched it.
!!
!! `ocean_surface_stress_t%stress_shelf` closes that: it is the
!! cell-centred ice-base stress, always allocated, zero without a cavity,
!! and both schemes now read
!!
!!     u_*^2 = (stress_mag + stress_shelf) / rho_0.
!!
!! What this suite pins:
!!
!!   * `ustar_under_lid_is_sqrt_cd_u` (KPP) and `epbl_ustar_under_lid`
!!     (EPBL) — under a full lid with zero wind, `u_*` is exactly
!!     `sqrt(C_d)*|U|`, the SAME number `cavity_ustar` gives the melt
!!     solve under the one-drag-coefficient rule; and an OPEN column on
!!     the same grid still gets the wind's `sqrt(|tau|/rho_0)`.
!!   * `quiescent_lid_is_the_documented_floor` — zero flow, zero tide,
!!     zero wind under the lid: KPP's `u_*` is 0 and its overlay is
!!     inert (no NaN, no 0/0), EPBL's is its `ustar_min` floor exactly.
!!   * `b0_carries_the_cavity_fluxes_kpp` / `_epbl` — the melt heat and
!!     salt fluxes reach the surface buoyancy flux of BOTH schemes,
!!     through the PRODUCTION assembler, with the analytically expected
!!     magnitude.  Sign convention (both schemes): **B_0 > 0 is
!!     STABILISING**.
!!   * `kpp_bl_deepens_with_ustar` / `epbl_mld_deepens_with_ustar` —
!!     more ice-ocean stress ⇒ deeper boundary layer.
!!   * `epbl_mld_shoals_with_melt_buoyancy` — a stabilising (melt-like)
!!     surface buoyancy flux shoals EPBL's boundary layer, monotonically.
!!   * `kpp_bl_is_insensitive_to_a_stabilising_b0` — the honest
!!     counterpart, and a FINDING rather than an omission: this KPP is
!!     the Large et al. (1994) bulk-Richardson depth, in which `B_0`
!!     enters ONLY through `destabilizing = (B_0 < 0)` and
!!     `w_*^3 = max(0, -B_0)*h_b`.  For any `B_0 >= 0` all three uses are
!!     identical to `B_0 = 0`, so a stabilising flux cannot shoal the KPP
!!     boundary layer WITHIN A CALL — it shoals it over time, through the
!!     stratification the fresh meltwater builds.  Asserted here so the
!!     behaviour is pinned rather than assumed; a future
!!     Monin-Obukhov/Ekman stable-depth limiter is what would change it.
!!   * `tau_writers_do_not_touch_stress_shelf` — `stress_shelf` is NOT
!!     derived from `tau`: the wind setters and the cover mask must leave
!!     it alone, which is the whole reason it is a second field.
!!   * `driver_publishes_shelf_stress_end_to_end` — the wiring gate, run
!!     through the public API: under a full lid with no wind, switching
!!     `&ocean_tdrag_nml` on raises KPP's `kv` by orders of magnitude,
!!     because the boundary layer now feels the ice.
!!   * `melt_only_fallback_publishes_shelf_stress` — the OTHER filler:
!!     with `&ocean_tdrag_nml` OFF and `&ocean_cavity_melt_nml` ON,
!!     `engine_step_finalize` publishes `rho_0*u_*^2` from the melt
!!     slot's own `u_*` instead, one outer step late.
!!   * `epbl_p_top_*` — the EPBL in-situ pressure port
!!     (`&ocean_psurf_nml in_eos`): the knob off is bit-identical even
!!     under a loaded column; on, the pressure that reaches the EOS is
!!     `>= p_top > 0` and increases toward the bed; a UNIFORM `p_top`
!!     under a LINEAR (pressure-independent) EOS is gauge-neutral to a
!!     derived bound; and under the nonlinear Wright EOS the same uniform
!!     load DOES move the answer, through `alpha(p)`/`beta(p)` alone.
!!
!! `mem:separate` discipline throughout: every kernel is given
!! device-present arrays (state objects through their own `enter_data`,
!! bare scratch through its own `!$acc enter data ... exit data`), host-set
!! inputs are pushed with `!$acc update device`, and every read-back is an
!! `!$acc update self` of COMPONENT arrays — never of an aggregate derived
!! type.  All of it is inert on a host build, which is exactly why it is
!! written unconditionally.
!!
!! Tolerances are STATED and derived; nothing here asserts bit-equality
!! across two different code paths.
module test_ocean_bl_under_ice
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, vmix_apply_kpp_overlay
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_compute, EPBL_MSTAR_CONSTANT
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_apply_cover
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_assemble, &
                                     SEAWATER_CP
   use rdb_ocean_top_drag, only: ocean_top_drag_t, TDRAG_QUADRATIC, &
                                 ocean_top_drag_compute_tendencies, &
                                 top_drag_fill_face_cover_impl
   use rdb_ocean_cavity_melt, only: cavity_ustar, CAVITY_MELT_OK
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_kv_ptr, rdb_ocean_set_u, &
                            rdb_ocean_get_grid_info
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_bl_under_ice_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 6
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.2_wp
      !! kg/m^3/K — the linear-EOS thermal contraction the vmix/EPBL
      !! slots are given here (`eos%alpha_T`).
   real(wp), parameter :: BETA_LIN = 0.8_wp
      !! kg/m^3/PSU (`eos%beta_S`).
   real(wp), parameter :: CD = 2.5e-3_wp
      !! ISOMIP+ ice-base drag coefficient (Asay-Davis et al. 2016), the
      !! ONE `C_d` shared by `&ocean_tdrag_nml cd` and
      !! `&ocean_cavity_melt_nml cdrag_top`.
   real(wp), parameter :: PTOP_HK = 40.0_wp
      !! Layer thickness (m) in the EPBL pressure-port cases.
   real(wp), parameter :: PTOP_DT = 1800.0_wp
      !! Thermo step (s) in the EPBL pressure-port cases.
   real(wp), parameter :: U_JET = 0.2_wp
      !! Depth-uniform zonal jet (m/s) seeded into the end-to-end runs —
      !! the stirrer that gives the ice base something to rub against.
   real(wp), parameter :: TOL_REL = 1.0e-11_wp
      !! Relative tolerance for a closed-form comparison.  The chains
      !! here are O(30) flops (a band integral, a square root, a cube
      !! root, a shape function), so the accumulated relative round-off
      !! is O(100*eps) ~ 2e-14 in double; 1e-11 is three decades of
      !! headroom and still several decades below any physical error.
      !! A cube root costs the most precision, which is why this is
      !! looser than the 1e-12 the pure-algebra drag tests use.

contains

   subroutine collect_ocean_bl_under_ice_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ustar_under_lid_is_sqrt_cd_u", test_kpp_ustar_lid), &
                  new_unittest("epbl_ustar_under_lid", test_epbl_ustar_lid), &
                  new_unittest("quiescent_lid_is_the_documented_floor", test_quiescent_floor), &
                  new_unittest("b0_carries_the_cavity_fluxes_kpp", test_b0_kpp), &
                  new_unittest("b0_carries_the_cavity_fluxes_epbl", test_b0_epbl), &
                  new_unittest("kpp_bl_deepens_with_ustar", test_kpp_deepens), &
                  new_unittest("kpp_bl_is_insensitive_to_a_stabilising_b0", test_kpp_stable_b0), &
                  new_unittest("epbl_mld_deepens_with_ustar", test_epbl_deepens), &
                  new_unittest("epbl_mld_shoals_with_melt_buoyancy", test_epbl_shoals), &
                  new_unittest("tau_writers_do_not_touch_stress_shelf", test_tau_writers), &
                  new_unittest("driver_publishes_shelf_stress_end_to_end", test_driver_publish), &
                  new_unittest("epbl_p_top_off_is_bit_identical", test_epbl_ptop_off), &
                  new_unittest("epbl_p_top_sign_and_monotone", test_epbl_ptop_monotone), &
                  new_unittest("epbl_uniform_p_top_linear_eos_is_gauge_neutral", &
                               test_epbl_ptop_gauge), &
                  new_unittest("epbl_p_top_moves_the_nonlinear_eos_only", test_epbl_ptop_wright), &
                  new_unittest("melt_only_fallback_publishes_shelf_stress", test_melt_fallback) &
                  ]
   end subroutine collect_ocean_bl_under_ice_tests

   ! ==================================================================
   ! Harness
   ! ==================================================================

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1000.0_wp, 1000.0_wp)
   end subroutine make_grid

   subroutine build_lid(grid, ms, ss, td, u_flow, h_layer, cover_split)
      !! A domain with a uniform zonal flow `u_flow` in EVERY layer,
      !! uniform layer thickness `h_layer`, and an ice shelf covering
      !! cells `i <= cover_split`.  Wind is left at zero by default; a
      !! caller that wants wind sets `ss%tau_*` and re-applies the cover.
      !!
      !! `htbl = 0` (layer-only top drag) and `drag_bg_vel = 0` make
      !! `|tau_top| = rho_0*C_d*U^2` exactly, with no band integral and
      !! no background-velocity floor in the way — the point of these
      !! tests is the SEAM, not the drag law, which `test_ocean_top_drag`
      !! already pins to its closed form.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_top_drag_t), intent(inout) :: td
      real(wp), intent(in) :: u_flow, h_layer
      integer, intent(in) :: cover_split
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid)
      ss%rho0 = RHO0
      td%enable = .true.
      call td%init(grid, nz_ml=NZ)
      td%variant = TDRAG_QUADRATIC
      td%c_drag = CD
      td%htbl = 0.0_wp
      td%drag_bg_vel = 0.0_wp
      td%rho0 = RHO0

      ms%h_layer = h_layer
      ms%u_face_x_layer = u_flow
      ms%v_face_y_layer = 0.0_wp

      td%cover_t(:, :) = 0.0_wp
      if (cover_split >= 1) td%cover_t(1:min(cover_split, nx), :) = 1.0_wp
      call top_drag_fill_face_cover_impl(td%cover_u, td%cover_v, td%cover_t, nx, ny)
   end subroutine build_lid

   subroutine publish_shelf_stress(ss, td)
      !! The driver's publish, verbatim: `ss%stress_shelf = td%stress_top`
      !! over the whole array.  In production this is an INLINE
      !! `do concurrent` in `run_stage` / `run_stage_split` (a call
      !! handing a state array to an external subroutine pessimises every
      !! `do concurrent` in the caller — CLAUDE.md); here it is a helper
      !! because a test routine has no hot kernels to pessimise.  The
      !! end-to-end gate `driver_publishes_shelf_stress_end_to_end` is
      !! what proves the PRODUCTION copy fires.
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_top_drag_t), intent(in) :: td
      integer :: i, j, nx, ny
      nx = size(ss%stress_shelf, 1)
      ny = size(ss%stress_shelf, 2)
      do concurrent(j=1:ny, i=1:nx)
         ss%stress_shelf(i, j) = td%stress_top(i, j)
      end do
   end subroutine publish_shelf_stress

   subroutine setup_kpp(vmix, grid)
      !! KPP on top of PP81, linear EOS, no V_t^2 (so the bulk-Ri
      !! crossing is the pure shear/buoyancy one and the BL depth is an
      !! analytic function of the column, not of `w_s`), unless a case
      !! turns `c_vt2` back on.
      type(ocean_vmix_t), intent(inout) :: vmix
      type(hgrid_t), intent(in) :: grid
      call vmix%init(grid, nz_ml=NZ)
      vmix%use_closure = .true.
      vmix%use_kpp = .true.
      vmix%rho0 = RHO0
      vmix%eos%variant = EOS_VARIANT_LINEAR
      vmix%eos%alpha_T = ALPHA_LIN
      vmix%eos%beta_S = BETA_LIN
      vmix%eos%rho0 = RHO0
   end subroutine setup_kpp

   subroutine setup_epbl(epbl, grid)
      type(ocean_epbl_t), intent(inout) :: epbl
      type(hgrid_t), intent(in) :: grid
      call epbl%init(grid, nz_ml=NZ)
      epbl%enable = .true.
      epbl%mstar_scheme = EPBL_MSTAR_CONSTANT
      epbl%mstar_const = 1.2_wp
      epbl%nstar = 0.2_wp
      epbl%use_lt = .false.
      epbl%eos%variant = EOS_VARIANT_LINEAR
      epbl%eos%alpha_T = ALPHA_LIN
      epbl%eos%beta_S = BETA_LIN
      epbl%eos%rho0 = RHO0
      epbl%rho0 = RHO0
      epbl%f_centre = 0.0_wp
      epbl%tke_diags = .true.
   end subroutine setup_epbl

   subroutine seed_uniform_TS(ms, h_layer, t_val, s_val)
      !! Uniform T and S — an UNSTRATIFIED column, so KPP finds no
      !! bulk-Ri crossing and `h_b` is the full column depth EXACTLY.
      !! That is what turns the overlay `kv = h_b*w_s*g(sigma)` into a
      !! direct, invertible read-out of `w_s`.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: h_layer, t_val, s_val
      integer :: k
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_val*h_layer
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = s_val*h_layer
         ms%rho_layer(:, :, k) = RHO0
      end do
   end subroutine seed_uniform_TS

   pure function kpp_ws_from_kv(kv_face, h_b, d_face) result(w_s)
      !! Invert the KPP overlay at one interface:
      !!   kv = h_b * w_s * sigma*(1-sigma)^2,   sigma = d_face/h_b
      !! so `w_s = kv / (h_b * g(sigma))`.  With zero surface buoyancy
      !! flux `w_s == u_*`, which is how these tests read `u_*` back out
      !! of the scheme instead of trusting an internal.
      real(wp), intent(in) :: kv_face, h_b, d_face
      real(wp) :: w_s, sigma, g_shape
      sigma = d_face/h_b
      g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
      w_s = kv_face/(h_b*g_shape)
   end function kpp_ws_from_kv

   subroutine run_kpp(grid, ms, vmix, ss, sf)
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
      !$acc update device(ss%stress_mag, ss%stress_shelf)
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth, vmix%gamma_t, vmix%gamma_s)
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
      !$acc update device(ss%stress_mag, ss%stress_shelf)
      call epbl_compute(grid, epbl, ms, ss, dt, sf=sf)
      !$acc update self(epbl%mld, epbl%kd_int, epbl%b0, epbl%tke_wind)
      call sf%exit_data()
      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, sf)
   end subroutine run_epbl

   subroutine run_top_drag(ms, td, ss)
      !! Compute `stress_top` on device and publish it into
      !! `ss%stress_shelf` — the sequence the stage drivers run.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_top_drag_t), intent(inout) :: td
      type(ocean_surface_stress_t), intent(inout) :: ss
      !$acc enter data copyin(ms, td)
      call ms%enter_data()
      call td%enter_data()
      !$acc update device(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, ms%wet_mask)
      !$acc update device(td%cover_u, td%cover_v, td%cover_t)
      call ocean_top_drag_compute_tendencies(td, ms, 1.0_wp)
      !$acc update self(td%stress_top)
      call td%exit_data()
      !$acc exit data delete(td)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call publish_shelf_stress(ss, td)
   end subroutine run_top_drag

   ! ==================================================================
   ! (i) u_* under the lid IS sqrt(C_d)*|U| — KPP
   ! ==================================================================

   subroutine test_kpp_ustar_lid(error)
      !! Full-strength statement of the fix, read back out of KPP itself.
      !!
      !! Unstratified column ⇒ no bulk-Ri crossing ⇒ `h_b` is the exact
      !! column depth.  Zero surface fluxes ⇒ `B_0 = 0` ⇒ `w_* = 0` ⇒
      !! `w_s = u_*`.  So inverting the shape function at any interior
      !! interface returns `u_*` directly.
      !!
      !! West of the calving front the wind is masked (`stress_mag = 0`)
      !! and the lid supplies `sqrt(C_d)*|U|`.  East of it the lid is
      !! absent (`stress_shelf = 0`) and the wind supplies
      !! `sqrt(|tau|/rho_0)`.  BOTH are asserted: a fix that simply
      !! replaced the wind everywhere would pass the first and fail the
      !! second.
      !!
      !! The third assertion is the one-drag-coefficient rule: the melt
      !! solve's own `cavity_ustar` returns the SAME number, so the
      !! boundary layer and the melt rate are being driven by one stress.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: U0 = 0.25_wp
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: TAU = 0.08_wp
      integer :: nx, ny, i_cov, i_open, j_p, k_probe, i_split, ierr_us
      real(wp) :: h_b, d_face, us_cov, us_open, us_melt, expect_cov, expect_open
      checks: block

         call make_grid(grid, 16, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         i_split = nx/2
         call build_lid(grid, ms, ss, td, U0, HK, i_split)
         call setup_kpp(vmix, grid)
         vmix%c_vt2 = 0.0_wp
         call sf%init(grid)

         call seed_uniform_TS(ms, HK, 1.0_wp, 34.5_wp)

         ! Wind everywhere, then masked by the SAME cover the top drag
         ! uses.  `apply_cover` zeroes the `tau` pair on every face
         ! touching a covered cell and refreshes `stress_mag` in the same
         ! call — the production contract.
         call ss%set_wind_stress_const(TAU, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)

         call run_top_drag(ms, td, ss)
         call run_kpp(grid, ms, vmix, ss, sf)

         j_p = ny/2
         i_cov = i_split/2
         i_open = min(i_split + 4, nx - 1)
         k_probe = NZ - 2
         h_b = real(NZ, wp)*HK
         d_face = real(NZ - k_probe + 1, wp)*HK

         call check(error, abs(vmix%bl_depth(i_cov, j_p) - h_b) <= TOL_REL*h_b, &
                    "an unstratified column must give h_b = the full depth "// &
                    "(else the shape-function inversion below is not valid)")
         if (allocated(error)) exit checks

         us_cov = kpp_ws_from_kv(vmix%kv(i_cov, j_p, k_probe), h_b, d_face)
         us_open = kpp_ws_from_kv(vmix%kv(i_open, j_p, k_probe), h_b, d_face)
         expect_cov = sqrt(CD)*U0
         expect_open = sqrt(TAU/RHO0)

         call check(error, abs(us_cov - expect_cov) <= TOL_REL*expect_cov, &
                    "KPP u_* under the lid is not sqrt(C_d)*|U| — the ice-ocean "// &
                    "stress is not reaching the boundary-layer scheme")
         if (allocated(error)) exit checks
         call check(error, abs(us_open - expect_open) <= TOL_REL*expect_open, &
                    "KPP u_* on the OPEN side is not sqrt(|tau|/rho_0) — the fix "// &
                    "must not replace the wind where there is no ice")
         if (allocated(error)) exit checks

         ! ... and the two are genuinely different numbers, else the two
         ! assertions above could both pass on one value.
         call check(error, abs(us_cov - us_open) > 0.1_wp*max(us_cov, us_open), &
                    "the covered and open friction velocities must differ by more "// &
                    "than a tolerance, else this test is vacuous")
         if (allocated(error)) exit checks

         ! The one-C_d rule: the melt solve's u_* is the same number.
         call cavity_ustar(U0, 0.0_wp, CD, 0.0_wp, 0.0_wp, us_melt, ierr_us)
         call check(error, ierr_us == CAVITY_MELT_OK, "cavity_ustar refused a clean input")
         if (allocated(error)) exit checks
         call check(error, abs(us_cov - us_melt) <= TOL_REL*us_melt, &
                    "the boundary-layer u_* and the MELT u_* disagree — the "// &
                    "one-drag-coefficient rule is broken")
         if (allocated(error)) exit checks

         ! Under the lid the wind really is gone: a fix that forgot the
         ! cover mask would show up here.
         call check(error, ss%stress_mag(i_cov, j_p) == 0.0_wp, &
                    "the wind stress magnitude must be EXACTLY zero under cover")
         if (allocated(error)) exit checks
         call check(error, ss%stress_shelf(i_open, j_p) == 0.0_wp, &
                    "the shelf stress must be EXACTLY zero on an open column")

      end block checks
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call vmix%destroy(); call ms%destroy()
   end subroutine test_kpp_ustar_lid

   ! ==================================================================
   ! (i) the same statement for EPBL
   ! ==================================================================

   subroutine test_epbl_ustar_lid(error)
      !! EPBL publishes its wind-TKE input directly:
      !!
      !!     tke_wind = mstar * dt * rho_0 * u_*^3
      !!
      !! and with `mstar_scheme = CONSTANT` and Langmuir off, `mstar` is
      !! the constant we set.  So `u_*` comes straight back out, with no
      !! inversion of the MLD solve at all.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: U0 = 0.25_wp
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: TAU = 0.08_wp
      real(wp), parameter :: DT = 900.0_wp
      integer :: nx, ny, i_cov, i_open, j_p, i_split
      real(wp) :: us_cov, us_open, expect_cov, expect_open, scal
      checks: block

         call make_grid(grid, 16, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         i_split = nx/2
         call build_lid(grid, ms, ss, td, U0, HK, i_split)
         call setup_epbl(epbl, grid)
         call sf%init(grid)
         call seed_uniform_TS(ms, HK, 1.0_wp, 34.5_wp)

         call ss%set_wind_stress_const(TAU, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)
         call run_top_drag(ms, td, ss)
         call run_epbl(grid, ms, epbl, ss, sf, DT)

         j_p = ny/2
         i_cov = i_split/2
         i_open = min(i_split + 4, nx - 1)
         scal = epbl%mstar_const*RHO0

         us_cov = (epbl%tke_wind(i_cov, j_p)/scal)**(1.0_wp/3.0_wp)
         us_open = (epbl%tke_wind(i_open, j_p)/scal)**(1.0_wp/3.0_wp)
         expect_cov = sqrt(CD)*U0
         expect_open = sqrt(TAU/RHO0)

         call check(error, abs(us_cov - expect_cov) <= TOL_REL*expect_cov, &
                    "EPBL u_* under the lid is not sqrt(C_d)*|U|")
         if (allocated(error)) exit checks
         call check(error, abs(us_open - expect_open) <= TOL_REL*expect_open, &
                    "EPBL u_* on the OPEN side is not sqrt(|tau|/rho_0)")
         if (allocated(error)) exit checks
         call check(error, abs(us_cov - us_open) > 0.1_wp*max(us_cov, us_open), &
                    "the covered and open friction velocities must differ, else "// &
                    "this test is vacuous")

      end block checks
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call epbl%destroy(); call ms%destroy()
   end subroutine test_epbl_ustar_lid

   ! ==================================================================
   ! (ii) the quiescent lid — the documented floor, and no NaN
   ! ==================================================================

   subroutine test_quiescent_floor(error)
      !! Zero flow, zero tide, zero wind, full lid.  `|tau_top| = 0`
      !! (`drag_bg_vel = 0`), so:
      !!
      !!   * KPP's `u_*` is 0.  Its floor is not a constant — it is the
      !!     `w_s > 0` guard, which makes the overlay INERT.  `kv` must
      !!     therefore be exactly the PP81 interior value, and finite.
      !!   * EPBL's floor IS a constant, `epbl%ustar_min`, applied by the
      !!     `max()` at the `u_*` site; `tke_wind` must be
      !!     `mstar*dt*rho_0*ustar_min^3` exactly, and `idecay`
      !!     (`tke_decay*|f|/u_*`) must not have divided by zero.
      !!
      !! Everything is checked for finiteness, because the failure this
      !! guards is a silent NaN in `kd_int` that only shows up decades
      !! later as a dead column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix, vmix_ref
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: DT = 900.0_wp
      integer :: nx, ny
      real(wp) :: tke_expect, scal
      checks: block

         call make_grid(grid, 12, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         call build_lid(grid, ms, ss, td, 0.0_wp, HK, nx)
         call setup_kpp(vmix, grid)
         call setup_kpp(vmix_ref, grid)
         vmix_ref%use_kpp = .false.
         call setup_epbl(epbl, grid)
         call epbl%set_f_centre(grid, 1.4e-4_wp, 0.0_wp, 0.0_wp)
         call sf%init(grid)
         call seed_uniform_TS(ms, HK, 1.0_wp, 34.5_wp)

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)
         call run_top_drag(ms, td, ss)

         call check(error, all(ss%stress_shelf == 0.0_wp), &
                    "a motionless lid with no background velocity must exert "// &
                    "EXACTLY zero stress")
         if (allocated(error)) exit checks

         ! --- KPP: overlay must be inert, and finite ---
         call run_kpp(grid, ms, vmix, ss, sf)
         call run_kpp(grid, ms, vmix_ref, ss, sf)
         call check(error, all(ieee_is_finite(vmix%kv)) .and. all(ieee_is_finite(vmix%kt)), &
                    "KPP produced a non-finite diffusivity at u_* = 0")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vmix%kv - vmix_ref%kv)) == 0.0_wp, &
                    "the KPP overlay is not inert at u_* = 0 — its documented "// &
                    "floor is the w_s > 0 guard, not a constant")
         if (allocated(error)) exit checks

         ! --- EPBL: the floor is ustar_min, exactly ---
         call run_epbl(grid, ms, epbl, ss, sf, DT)
         call check(error, all(ieee_is_finite(epbl%kd_int)) .and. &
                    all(ieee_is_finite(epbl%mld)), &
                    "EPBL produced a non-finite kd_int / mld at u_* = 0 — check "// &
                    "the idecay = tke_decay*|f|/u_* division")
         if (allocated(error)) exit checks
         scal = epbl%mstar_const*RHO0
         tke_expect = scal*epbl%ustar_min**3
         call check(error, abs(epbl%tke_wind(nx/2, ny/2) - tke_expect) <= &
                    TOL_REL*max(tke_expect, tiny(1.0_wp)), &
                    "EPBL did not floor u_* at ustar_min under a motionless lid")

      end block checks
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call epbl%destroy(); call vmix_ref%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_quiescent_floor

   ! ==================================================================
   ! (iii) B_0 carries the cavity fluxes — through the PRODUCTION assembler
   ! ==================================================================

   subroutine assemble_cavity_fluxes(grid, ms, sf, q_heat_cav, q_salt_cav)
      !! Put a melt-like heat and salt flux into the two components the
      !! cavity OWNS, then run the production assembler with full cover.
      !! That is exactly `engine_step_finalize`'s sequence, so what comes
      !! out of `sf%Q_heat` / `sf%Q_salt` is what the schemes see in
      !! production.  Full cover ⇒ every ATMOSPHERIC band is multiplied
      !! by `1 - cover_frac = 0`, so `Q_heat == heat_cavity` and
      !! `Q_salt == salt_cavity`.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: q_heat_cav, q_salt_cav
      real(wp), allocatable :: cover(:, :)

      allocate (cover(grid%nx_total, grid%ny_total), source=1.0_wp)
      sf%heat_cavity = q_heat_cav
      sf%salt_cavity = q_salt_cav
      sf%has_heat = .true.
      sf%has_salt = .true.
      !$acc enter data copyin(ms, sf, cover)
      call ms%enter_data()
      call sf%enter_data()
      !$acc update device(sf%heat_cavity, sf%salt_cavity)
      call ocean_surface_flux_assemble(grid, sf, ms, cover_frac=cover)
      !$acc update self(sf%Q_heat, sf%Q_salt)
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf, cover)
      deallocate (cover)
   end subroutine assemble_cavity_fluxes

   subroutine test_b0_kpp(error)
      !! KPP's surface buoyancy flux is
      !!
      !!     B_0 = g*(alpha_T*Q_heat/(rho_0*cp) - beta_S*Q_salt/rho_0),
      !!     B_0 > 0 STABILISING          (rdb_ocean_vmix, KPP pass 1/2)
      !!
      !! and it is built from `sf%Q_heat` / `sf%Q_salt` — the ASSEMBLER's
      !! outputs — so the cavity's `heat_cavity` / `salt_cavity` reach it
      !! with no change to KPP at all.  This test proves that rather than
      !! assuming it, by RECOVERING `B_0` from the scheme's own output.
      !!
      !! To make `B_0` observable the case is DESTABILISING (brine
      !! rejection: a positive salt flux into the ocean), because that is
      !! the branch KPP responds to.  With zero stress, `u_* = 0`, so
      !! `w_s = w_* = (-B_0*h_b)^(1/3)` exactly and the shape-function
      !! inversion returns `w_*`.
      !!
      !! The stabilising (melting) sign is covered by
      !! `kpp_bl_is_insensitive_to_a_stabilising_b0` — a documented
      !! property of this KPP, not an omission.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: Q_HEAT_CAV = -40.0_wp
         !! W/m^2, negative = melting takes heat OUT of the ocean
         !! (`heat_cavity = -q_ocean`).  Cooling the surface is
         !! DEstabilising on its own.
      real(wp), parameter :: Q_SALT_CAV = 2.0e-3_wp
         !! kg/m^2/s of salt INTO the ocean — brine rejection, the
         !! freezing (marine-ice-accretion) sign.  Chosen so the net
         !! `B_0` is clearly destabilising and the recovery below is
         !! well conditioned.
      integer :: nx, ny, i_p, j_p, k_probe
      real(wp) :: h_b, d_face, w_s, b0_got, b0_expect, q_t_kin, q_s_kin
      checks: block

         call make_grid(grid, 12, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         call build_lid(grid, ms, ss, td, 0.0_wp, HK, nx)
         call setup_kpp(vmix, grid)
         vmix%c_vt2 = 0.0_wp
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call seed_uniform_TS(ms, HK, 1.0_wp, 34.5_wp)

         call assemble_cavity_fluxes(grid, ms, sf, Q_HEAT_CAV, Q_SALT_CAV)

         i_p = nx/2
         j_p = ny/2
         call check(error, sf%Q_heat(i_p, j_p) == Q_HEAT_CAV .and. &
                    sf%Q_salt(i_p, j_p) == Q_SALT_CAV, &
                    "under FULL cover the assembler must deliver the cavity "// &
                    "components and nothing else")
         if (allocated(error)) exit checks

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)
         call run_kpp(grid, ms, vmix, ss, sf)

         k_probe = NZ - 2
         h_b = real(NZ, wp)*HK
         d_face = real(NZ - k_probe + 1, wp)*HK
         call check(error, abs(vmix%bl_depth(i_p, j_p) - h_b) <= TOL_REL*h_b, &
                    "an unstratified column must give h_b = the full depth")
         if (allocated(error)) exit checks

         ! u_* = 0 here, so w_s == w_* == (-B_0*h_b)^(1/3).
         w_s = kpp_ws_from_kv(vmix%kv(i_p, j_p, k_probe), h_b, d_face)
         b0_got = -(w_s**3)/h_b

         q_t_kin = sf%Q_heat(i_p, j_p)/(RHO0*sf%cp)
         q_s_kin = sf%Q_salt(i_p, j_p)/RHO0
         b0_expect = GRAVITY*(ALPHA_LIN*q_t_kin - BETA_LIN*q_s_kin)

         call check(error, b0_expect < 0.0_wp, &
                    "the test case must be DESTABILISING, else KPP does not "// &
                    "respond to B_0 at all and the recovery is vacuous")
         if (allocated(error)) exit checks
         call check(error, abs(b0_got - b0_expect) <= TOL_REL*abs(b0_expect), &
                    "KPP's B_0 does not match g*(alpha*Q_heat/(rho0*cp) - "// &
                    "beta*Q_salt/rho0) built from the ASSEMBLED cavity fluxes")
         if (allocated(error)) exit checks

         ! The non-local transport is the other consumer of the same two
         ! kinematic fluxes; it must carry the cavity salt flux too.
         call check(error, abs(vmix%gamma_s(i_p, j_p, k_probe)) > 0.0_wp, &
                    "the KPP non-local term is dead under a destabilising "// &
                    "cavity flux")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%gamma_s(i_p, j_p, k_probe)/q_s_kin - &
                               vmix%gamma_t(i_p, j_p, k_probe)/q_t_kin) <= &
                    TOL_REL*abs(vmix%gamma_s(i_p, j_p, k_probe)/q_s_kin), &
                    "gamma_t and gamma_s must share one shape factor — they are "// &
                    "gamma_factor times the two kinematic surface fluxes")

      end block checks
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call vmix%destroy(); call ms%destroy()
   end subroutine test_b0_kpp

   subroutine test_b0_epbl(error)
      !! EPBL stores its surface buoyancy flux verbatim (`epbl%b0`,
      !! persisted for the Bodner MLE convective velocity scale), so no
      !! inversion is needed:
      !!
      !!     b0 = g*rho_0*(dSV/dT * q_T + dSV/dS * q_S),  > 0 STABILISING
      !!
      !! With a LINEAR EOS the specific-volume derivatives are
      !! `dSV/dT = alpha_T/rho^2` and `dSV/dS = -beta_S/rho^2` at the
      !! surface layer's own density, so the expected value is analytic.
      !!
      !! Here the MELTING sign is used — cold AND fresh — and the
      !! assertion is that the net `b0` comes out POSITIVE: melting under
      !! an ice shelf is stabilising, the freshening beating the cooling.
      !! That is the physical claim of Phase 4b, and it is checked
      !! against a number rather than a sign alone.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: DT = 900.0_wp
      real(wp), parameter :: T_COL = 1.0_wp, S_COL = 34.5_wp
      real(wp), parameter :: Q_HEAT_CAV = -40.0_wp
         !! W/m^2 — melting cools the ocean.
      real(wp), parameter :: Q_SALT_CAV = -3.0e-3_wp
         !! kg/m^2/s — melting FRESHENS (a virtual salt sink).
      integer :: nx, ny, i_p, j_p
      real(wp) :: q_t_kin, q_s_kin, dsv_dt, dsv_ds, rho_s, b0_expect
      checks: block

         call make_grid(grid, 12, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         call build_lid(grid, ms, ss, td, 0.0_wp, HK, nx)
         call setup_epbl(epbl, grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call seed_uniform_TS(ms, HK, T_COL, S_COL)

         call assemble_cavity_fluxes(grid, ms, sf, Q_HEAT_CAV, Q_SALT_CAV)

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)
         call run_epbl(grid, ms, epbl, ss, sf, DT)

         i_p = nx/2
         j_p = ny/2
         q_t_kin = sf%Q_heat(i_p, j_p)/(RHO0*sf%cp)
         q_s_kin = sf%Q_salt(i_p, j_p)/RHO0
         ! Linear EOS: rho = rho0 + beta_S*(S - S_ref) - alpha_T*(T - T_ref),
         ! so `eos_specvol_derivs` returns the CONSTANT sensitivities
         ! dSV/dT = +alpha_T/rho0^2 and dSV/dS = -beta_S/rho0^2, evaluated
         ! at the REFERENCE density (the Boussinesq weights that consume
         ! them are referenced there too) — not at the local density.
         rho_s = epbl%eos%rho0
         dsv_dt = ALPHA_LIN/(rho_s*rho_s)
         dsv_ds = -BETA_LIN/(rho_s*rho_s)
         b0_expect = GRAVITY*RHO0*(dsv_dt*q_t_kin + dsv_ds*q_s_kin)

         call check(error, abs(epbl%b0(i_p, j_p) - b0_expect) <= &
                    TOL_REL*abs(b0_expect), &
                    "EPBL's b0 does not match g*rho0*(dSV/dT*q_T + dSV/dS*q_S) "// &
                    "built from the ASSEMBLED cavity fluxes")
         if (allocated(error)) exit checks
         call check(error, epbl%b0(i_p, j_p) > 0.0_wp, &
                    "basal MELT must be STABILISING (b0 > 0): the freshening "// &
                    "must beat the cooling")
         if (allocated(error)) exit checks

         ! Non-vacuity: both bands actually contribute.
         call check(error, GRAVITY*RHO0*dsv_dt*q_t_kin < 0.0_wp .and. &
                    GRAVITY*RHO0*dsv_ds*q_s_kin > 0.0_wp, &
                    "the melt case must carry a destabilising HEAT band and a "// &
                    "stabilising SALT band, else the sign test is trivial")

      end block checks
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call epbl%destroy(); call ms%destroy()
   end subroutine test_b0_epbl

   ! ==================================================================
   ! (v) monotonicity
   ! ==================================================================

   subroutine kpp_depth_at_flow(u_flow, h_b, error)
      !! Stratified column under a full lid; returns KPP's boundary-layer
      !! depth at the probe column for a given uniform flow speed.  The
      !! flow is uniform in `k`, so it contributes NO resolved shear —
      !! the only way it can change `h_b` is through `u_*` feeding `w_s`
      !! and thence `V_t^2`, which is exactly the coupling under test.
      real(wp), intent(in) :: u_flow
      real(wp), intent(out) :: h_b
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      integer :: nx, ny, k

      call make_grid(grid, 12, 8)
      nx = grid%nx_total
      ny = grid%ny_total
      call build_lid(grid, ms, ss, td, u_flow, HK, nx)
      call setup_kpp(vmix, grid)
      call sf%init(grid)
      ! Stable stratification, bottom-up: k = NZ is the ICE BASE (the
      ! lightest water), k = 1 the bed.
      do k = 1, NZ
         ms%rho_layer(:, :, k) = RHO0 + 0.04_wp*real(NZ - k, wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 1.0_wp*HK
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 34.5_wp*HK
      end do
      ! Two KPP passes so `bl_depth` is no longer at its zero seed when
      ! the V_t^2 term reads the lagged depth.
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      call ocean_surface_stress_apply_cover(ss, td%cover_t)
      call run_top_drag(ms, td, ss)
      call run_kpp(grid, ms, vmix, ss, sf)
      call run_kpp(grid, ms, vmix, ss, sf)
      h_b = vmix%bl_depth(nx/2, ny/2)
      call check(error, ieee_is_finite(h_b) .and. h_b > 0.0_wp, &
                 "KPP returned a non-positive or non-finite boundary-layer depth")
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call vmix%destroy(); call ms%destroy()
   end subroutine kpp_depth_at_flow

   subroutine test_kpp_deepens(error)
      !! More flow under the lid ⇒ more ice-ocean stress ⇒ larger `u_*`
      !! ⇒ larger `w_s` ⇒ larger `V_t^2` ⇒ the bulk Richardson number
      !! crosses `ri_crit` deeper.  A monotone ladder in `U`, not a
      !! two-point comparison, so a non-monotone response cannot hide.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U_LADDER(4) = [0.02_wp, 0.08_wp, 0.20_wp, 0.50_wp]
      real(wp) :: h_b(4)
      integer :: n
      checks: block
         do n = 1, 4
            call kpp_depth_at_flow(U_LADDER(n), h_b(n), error)
            if (allocated(error)) exit checks
         end do
         do n = 2, 4
            call check(error, h_b(n) > h_b(n - 1), &
                       "the KPP boundary layer under an ice shelf must DEEPEN "// &
                       "with the ice-ocean stress")
            if (allocated(error)) exit checks
         end do
         call check(error, h_b(4) > 1.5_wp*h_b(1), &
                    "the deepening across the ladder must be a real effect, not "// &
                    "a round-off-sized one")
      end block checks
   end subroutine test_kpp_deepens

   subroutine test_kpp_stable_b0(error)
      !! THE HONEST COUNTERPART, and a documented finding.
      !!
      !! In the Large et al. (1994) bulk-Richardson boundary-layer depth
      !! this model implements, `B_0` is consumed in exactly three
      !! places, all of them `max(0, -B_0)` or `B_0 < 0`.  A STABILISING
      !! `B_0` (which is what basal melting produces — see
      !! `b0_carries_the_cavity_fluxes_epbl`) is therefore
      !! indistinguishable from no buoyancy flux at all WITHIN one call:
      !! KPP shoals under melt over TIME, through the stratification the
      !! meltwater builds, not instantaneously through `B_0`.
      !!
      !! This is asserted so the property is pinned.  Adding a
      !! Monin-Obukhov / Ekman stable-depth limiter (MOM6 carries one;
      !! this port does not) is what would legitimately change it, and
      !! that change must then fail this test and update it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_zero, vmix_stab
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf_zero, sf_stab
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: U0 = 0.2_wp
      integer :: nx, ny, k
      checks: block

         call make_grid(grid, 12, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         call build_lid(grid, ms, ss, td, U0, HK, nx)
         call setup_kpp(vmix_zero, grid)
         call setup_kpp(vmix_stab, grid)
         call sf_zero%init(grid)
         call sf_stab%init(grid)
         call sf_stab%set_components(grid, .true.)
         do k = 1, NZ
            ms%rho_layer(:, :, k) = RHO0 + 0.04_wp*real(NZ - k, wp)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 1.0_wp*HK
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 34.5_wp*HK
         end do

         ! A strongly STABILISING melt-like flux: cold and very fresh.
         call assemble_cavity_fluxes(grid, ms, sf_stab, -40.0_wp, -2.0e-2_wp)
         call check(error, GRAVITY*(ALPHA_LIN*sf_stab%Q_heat(nx/2, ny/2)/(RHO0*sf_stab%cp) &
                                    - BETA_LIN*sf_stab%Q_salt(nx/2, ny/2)/RHO0) > 0.0_wp, &
                    "the control case must be STABILISING, else this test says "// &
                    "nothing")
         if (allocated(error)) exit checks

         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         call ocean_surface_stress_apply_cover(ss, td%cover_t)
         call run_top_drag(ms, td, ss)
         call run_kpp(grid, ms, vmix_zero, ss, sf_zero)
         call run_kpp(grid, ms, vmix_stab, ss, sf_stab)

         call check(error, vmix_zero%bl_depth(nx/2, ny/2) == &
                    vmix_stab%bl_depth(nx/2, ny/2), &
                    "this KPP's boundary-layer depth responded to a STABILISING "// &
                    "B_0 — it consumes B_0 only through max(0,-B_0) and "// &
                    "(B_0 < 0), so either the kernel changed or a stable-depth "// &
                    "limiter landed; update this test with the change")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vmix_zero%gamma_t)) == 0.0_wp, &
                    "the KPP non-local term must be dead under a stabilising flux")

      end block checks
      call sf_stab%destroy(); call sf_zero%destroy(); call td%destroy()
      call ss%destroy(); call vmix_stab%destroy(); call vmix_zero%destroy()
      call ms%destroy()
   end subroutine test_kpp_stable_b0

   subroutine epbl_mld_case(u_flow, q_heat_cav, q_salt_cav, mld, error)
      !! One EPBL call on a stratified column under a full lid.
      real(wp), intent(in) :: u_flow, q_heat_cav, q_salt_cav
      real(wp), intent(out) :: mld
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: HK = 20.0_wp
      real(wp), parameter :: DT = 1800.0_wp
      integer :: nx, ny, k
      real(wp) :: t_k

      call make_grid(grid, 12, 8)
      nx = grid%nx_total
      ny = grid%ny_total
      call build_lid(grid, ms, ss, td, u_flow, HK, nx)
      call setup_epbl(epbl, grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      ! Stable T stratification: warmest at the ice base (k = NZ).
      do k = 1, NZ
         t_k = 1.0_wp - 0.02_wp*real(NZ - k, wp)*HK
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*HK
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 34.5_wp*HK
         ms%rho_layer(:, :, k) = RHO0
      end do

      call assemble_cavity_fluxes(grid, ms, sf, q_heat_cav, q_salt_cav)
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      call ocean_surface_stress_apply_cover(ss, td%cover_t)
      call run_top_drag(ms, td, ss)
      call run_epbl(grid, ms, epbl, ss, sf, DT)

      mld = epbl%mld(nx/2, ny/2)
      call check(error, ieee_is_finite(mld) .and. mld > 0.0_wp, &
                 "EPBL returned a non-positive or non-finite MLD")
      call sf%destroy(); call td%destroy(); call ss%destroy()
      call epbl%destroy(); call ms%destroy()
   end subroutine epbl_mld_case

   subroutine test_epbl_deepens(error)
      !! More ice-ocean stress ⇒ more mechanical TKE (`mstar*dt*rho_0*
      !! u_*^3`) ⇒ EPBL mixes deeper.  Monotone ladder in `U`.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U_LADDER(4) = [0.02_wp, 0.08_wp, 0.20_wp, 0.50_wp]
      real(wp) :: mld(4)
      integer :: n
      checks: block
         do n = 1, 4
            call epbl_mld_case(U_LADDER(n), 0.0_wp, 0.0_wp, mld(n), error)
            if (allocated(error)) exit checks
         end do
         do n = 2, 4
            call check(error, mld(n) > mld(n - 1), &
                       "EPBL's boundary layer under an ice shelf must DEEPEN "// &
                       "with the ice-ocean stress")
            if (allocated(error)) exit checks
         end do
         call check(error, mld(4) > 1.5_wp*mld(1), &
                    "the deepening across the ladder must be a real effect")
      end block checks
   end subroutine test_epbl_deepens

   subroutine test_epbl_shoals(error)
      !! The melt-buoyancy half of item (v): at FIXED flow, a stronger
      !! stabilising (melt-like) surface buoyancy flux must SHOAL the
      !! EPBL boundary layer, monotonically.  EPBL charges the freshly
      !! applied skin fluxes to its TKE ledger (`ctke_sfc`), and a
      !! stabilising flux makes that a DEBIT clipped against the
      !! mechanical reservoir — the physical statement that fresh, cold
      !! meltwater at the ice base costs energy to mix down.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U0 = 0.20_wp
      real(wp), parameter :: SALT_LADDER(4) = &
                             [0.0_wp, -2.0e-3_wp, -6.0e-3_wp, -2.0e-2_wp]
         !! kg/m^2/s of salt into the ocean; negative = FRESHENING, which
         !! is the melt sign and is stabilising.
      real(wp) :: mld(4)
      integer :: n
      checks: block
         do n = 1, 4
            call epbl_mld_case(U0, 0.0_wp, SALT_LADDER(n), mld(n), error)
            if (allocated(error)) exit checks
         end do
         do n = 2, 4
            call check(error, mld(n) < mld(n - 1), &
                       "EPBL's boundary layer must SHOAL as the stabilising melt "// &
                       "buoyancy flux grows")
            if (allocated(error)) exit checks
         end do
         call check(error, mld(4) < 0.9_wp*mld(1), &
                    "the shoaling across the ladder must be a real effect, not "// &
                    "a round-off-sized one")
      end block checks
   end subroutine test_epbl_shoals

   ! ==================================================================
   ! The separation of the two fields
   ! ==================================================================

   subroutine test_tau_writers(error)
      !! `stress_shelf` is NOT derived from `tau`, and that is the whole
      !! reason it is a second field: a `tau` writer rebuilds
      !! `stress_mag` from scratch, and a folded-in ice stress would be
      !! silently wiped at the next data-forcing bracket or sea-ice
      !! blend.  So: set a shelf stress, then run every `tau` writer this
      !! module exposes, and demand the shelf stress survives BIT for bit
      !! (this is one path leaving an array untouched, not two paths
      !! agreeing, so bit-equality is the right assertion).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), allocatable :: cover(:, :), shelf0(:, :)
      integer :: nx, ny, i, j
      checks: block

         call make_grid(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         call ss%init(grid)

         call check(error, all(ss%stress_shelf == 0.0_wp), &
                    "stress_shelf must be zero-filled at init — that zero is what "// &
                    "makes a cavity-free run bit-identical")
         if (allocated(error)) exit checks

         do j = 1, ny
            do i = 1, nx
               ss%stress_shelf(i, j) = 0.1_wp + 0.01_wp*real(i + j, wp)
            end do
         end do
         allocate (shelf0, source=ss%stress_shelf)

         call ss%set_wind_stress_const(0.05_wp, -0.02_wp)
         call check(error, all(ss%stress_shelf == shelf0), &
                    "a wind-stress setter moved stress_shelf")
         if (allocated(error)) exit checks

         call ss%set_wind_stress_2gyre(grid, 0.1_wp)
         call check(error, all(ss%stress_shelf == shelf0), &
                    "the 2-gyre wind setter moved stress_shelf")
         if (allocated(error)) exit checks

         allocate (cover(nx, ny), source=0.0_wp)
         cover(1:nx/2, :) = 1.0_wp
         call ocean_surface_stress_apply_cover(ss, cover)
         call check(error, all(ss%stress_shelf == shelf0), &
                    "the ice-shelf COVER mask moved stress_shelf — it masks the "// &
                    "wind, it does not own the ice stress")
         if (allocated(error)) exit checks
         ! `2:` because the WEST RIM face (i = 1) carries no cover mask
         ! by design (the face rule is stated once in `rdb_ocean_top_drag`
         ! and fills interior faces only), so cell 1 keeps half its wind.
         call check(error, all(ss%stress_mag(2:nx/2, :) == 0.0_wp), &
                    "the cover mask must still zero the wind magnitude it owns")

      end block checks
      call ss%destroy()
   end subroutine test_tau_writers

   ! ==================================================================
   ! The driver wiring, end to end through the public API
   ! ==================================================================

   subroutine test_driver_publish(error)
      !! THE WIRING GATE.  Everything above calls the kernels directly;
      !! this one runs the real split solver through `rdb_ocean_api` and
      !! asks whether the stage driver actually published the shelf
      !! stress into `ss%stress_shelf` before `vmix_apply_in_stage` read
      !! it.
      !!
      !! Configuration: a flat 300 m ice draft over the WHOLE domain, no
      !! wind (and under full cover there could be none anyway), a seeded
      !! zonal jet, KPP on.  Two runs, identical but for
      !! `&ocean_tdrag_nml enable`.
      !!
      !! With the top drag OFF the covered column has `stress_mag = 0`
      !! and `stress_shelf = 0`, so `u_* = 0` and the KPP overlay is
      !! INERT: `kv` can rise no higher than what the PP81 interior
      !! closure alone produces, which is bounded above by its own
      !! `nu_0 + kv_bg` ceiling — measured `1.01e-2 m^2/s`, i.e. the
      !! closure is already saturated and MORE shear cannot raise it.
      !! That bound is what makes the comparison a test of the seam and
      !! not of the top drag's other effect on the flow.
      !!
      !! With it ON, `u_* = sqrt(C_d)*|U| ~ 1e-2 m/s` over the
      !! boundary layer gives `kv = h_b*u_**g(sigma)` — measured
      !! `2.75e-1 m^2/s`, a factor 27 above the saturated interior
      !! closure.  The assertions below are `kv_off` at or below the
      !! interior ceiling, `kv_on` an order of magnitude past it, and a
      !! ratio of at least 5 (a fifth of the measured margin).
      !!
      !! Measured, gfortran 15.1 Release, 2026-09-20:
      !!   kv_off = 1.010e-02 m^2/s,  kv_on = 2.753e-01 m^2/s  (27.3x)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kv_off, kv_on
      checks: block
         call kv_max_under_lid(.false., kv_off, error)
         if (allocated(error)) exit checks
         call kv_max_under_lid(.true., kv_on, error)
         if (allocated(error)) exit checks

         call check(error, kv_off <= 2.0e-2_wp, &
                    "with no top drag the covered column must see u_* = 0 and "// &
                    "therefore no KPP overlay (kv stays at or below the "// &
                    "saturated PP81 interior ceiling)")
         if (allocated(error)) exit checks
         call check(error, kv_on > 1.0e-1_wp, &
                    "with the top drag on, the stage driver must publish the "// &
                    "ice-base stress into ss%stress_shelf before KPP reads it — "// &
                    "kv shows no boundary-layer overlay at all")
         if (allocated(error)) exit checks
         call check(error, kv_on > 5.0_wp*max(kv_off, tiny(1.0_wp)), &
                    "the ON/OFF separation must be an order of magnitude, not a "// &
                    "shear tweak on a saturated interior closure")
      end block checks
   end subroutine test_driver_publish

   subroutine kv_max_under_lid(tdrag_on, kv_max, error, melt_on, n_steps)
      !! One end-to-end run; returns the largest `kv` over the interior.
      !! Every read goes through `rdb_ocean_refresh_host` first — a getter
      !! alone never triggers a device->host copy, so skipping it would
      !! compare stale host memory on the GPU build.
      logical, intent(in) :: tdrag_on
      real(wp), intent(out) :: kv_max
      type(error_type), allocatable, intent(inout) :: error
      logical, intent(in), optional :: melt_on
      integer, intent(in), optional :: n_steps
         !! Steps to run.  The MELT fallback publishes at the END of a
         !! step, so it needs at least two before the boundary-layer
         !! schemes can have read it — that one-step lag is the documented
         !! cost of the fallback path.
      type(c_ptr) :: handle, ptr
      real(wp), pointer :: kv(:, :, :)
      real(wp), allocatable :: ubuf(:, :, :)
      integer(c_int) :: status, nx, ny, nz, gen, nx_p, ny_p, nz_p, ng
      character(len=:), allocatable :: nml
      integer :: ng_i, n_run

      kv_max = 0.0_wp
      handle = c_null_ptr
      n_run = 3
      if (present(n_steps)) n_run = n_steps
      if (present(melt_on)) then
         nml = lid_namelist(tdrag_on, melt_on=melt_on)
      else
         nml = lid_namelist(tdrag_on)
      end if
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, &
                 "the end-to-end cavity + KPP namelist must configure")
      if (allocated(error)) return

      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      if (status == OCEAN_STATUS_OK) then
         allocate (ubuf(nx_p + 1, ny_p, nz_p), source=U_JET)
         status = rdb_ocean_set_u(handle, ubuf, nx_p, ny_p, nz_p)
         deallocate (ubuf)
      end if
      if (status == OCEAN_STATUS_OK) then
         status = rdb_ocean_step(handle, int(n_run, c_int))
      end if
      if (status == OCEAN_STATUS_OK) status = rdb_ocean_refresh_host(handle)
      if (status == OCEAN_STATUS_OK) then
         status = rdb_ocean_get_kv_ptr(handle, ptr, nx, ny, nz, gen)
      end if
      if (status == OCEAN_STATUS_OK) then
         call c_f_pointer(ptr, kv, [int(nx), int(ny), int(nz)])
         ng_i = int(ng)
         kv_max = maxval(kv(ng_i + 1:ng_i + int(nx_p), ng_i + 1:ng_i + int(ny_p), :))
      end if
      call check(error, status == OCEAN_STATUS_OK, &
                 "the end-to-end run did not step / read back cleanly")
      if (.not. allocated(error)) then
         call check(error, ieee_is_finite(kv_max), "kv came back non-finite")
      end if
      status = rdb_ocean_destroy(handle)
   end subroutine kv_max_under_lid

   function lid_namelist(tdrag_on, melt_on) result(nml)
      !! A flat, fully covered ice shelf over a flat bed, with KPP and
      !! the split solver on its DEFAULT scheme (`pred_corr`) — the scheme
      !! every cavity configuration runs; nothing here is outside its
      !! envelope, so it is not pinned.
      logical, intent(in) :: tdrag_on
      logical, intent(in), optional :: melt_on
         !! Turn on `&ocean_cavity_melt_nml` instead of / as well as the
         !! top drag.  Melt needs the PR-12 component set and the ISOMIP+
         !! liquidus, both already in the base namelist.
      character(len=:), allocatable :: nml
      logical :: melt

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 12, ny = 10, nghost = 2, dx = 5000.0, dy = 5000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 6 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0 /"//new_line("a")// &
            "&ocean_ic_nml rho_0 = 1035.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., "// &
            "n_inner = 12 /"//new_line("a")// &
            "&ocean_vmix_nml use_closure = .true., use_kpp = .true. /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")// &
            "&ocean_eos_nml tfreeze_set = 'isomip' /"//new_line("a")// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', "// &
            "draft_depth = 300.0 /"//new_line("a")
      if (tdrag_on) then
         nml = nml//"&ocean_tdrag_nml enable = .true., cd = 2.5e-3 /"//new_line("a")
      end if
      melt = .false.
      if (present(melt_on)) melt = melt_on
      if (melt) then
         nml = nml//"&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
               "&ocean_cavity_melt_nml enable = .true., cdrag_top = 2.5e-3, "// &
               "u_tide = 0.0 /"//new_line("a")
      end if
   end function lid_namelist

   subroutine test_melt_fallback(error)
      !! THE OTHER FILLER.  `&ocean_tdrag_nml` is OFF here and
      !! `&ocean_cavity_melt_nml` is ON, so the stage drivers publish
      !! nothing and `engine_step_finalize` fills `ss%stress_shelf` from
      !! `rho_0*u_*^2` with the MELT slot's own `u_*` instead — the same
      !! `C_d` under the one-drag-coefficient rule.  Without it a
      !! melt-only cavity would mix its covered columns on `u_* = 0`,
      !! which is the whole defect Phase 4b closes.
      !!
      !! That fill runs at the END of an outer step, so it reaches the
      !! boundary-layer schemes ONE STEP LATE.  The test therefore runs
      !! five steps, and asserts the same decisive separation the
      !! top-drag gate uses: `kv` rises an order of magnitude past the
      !! saturated PP81 interior ceiling.
      !!
      !! `u_tide = 0` so the melt `u_*` is exactly `sqrt(C_d)*|U_far|`,
      !! the same number the top drag would have published.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kv_plain, kv_melt
      checks: block
         call kv_max_under_lid(.false., kv_plain, error, melt_on=.false., n_steps=5)
         if (allocated(error)) exit checks
         call kv_max_under_lid(.false., kv_melt, error, melt_on=.true., n_steps=5)
         if (allocated(error)) exit checks

         call check(error, kv_plain <= 2.0e-2_wp, &
                    "the no-cavity-forcing control must stay at the saturated "// &
                    "PP81 interior ceiling")
         if (allocated(error)) exit checks
         call check(error, kv_melt > 1.0e-1_wp, &
                    "with melt on and the top drag OFF, engine_step_finalize must "// &
                    "still publish an under-ice u_* into ss%stress_shelf — KPP "// &
                    "shows no boundary-layer overlay at all")
         if (allocated(error)) exit checks
         call check(error, kv_melt > 5.0_wp*max(kv_plain, tiny(1.0_wp)), &
                    "the melt-fallback separation must be an order of magnitude")
      end block checks
   end subroutine test_melt_fallback

   ! ==================================================================
   ! The EPBL in-situ pressure port (`&ocean_psurf_nml in_eos`)
   ! ==================================================================

   subroutine epbl_ptop_run(in_eos, wright, p_top_val, epbl, error)
      !! One EPBL call on a stratified, wind-forced column with a uniform
      !! top-of-column load `p_top_val` (Pa).  The slot is returned to
      !! the caller (host arrays already pulled back) so a case can read
      !! `kd_int`, `mld`, and the PE/steric weights `dpe_t`/`dcolht_t`.
      !!
      !! `ms%p_top` is filled unconditionally: it is ALWAYS allocated and
      !! mapped (the `p_top` seam contract), so there is no placeholder to
      !! keep out of a kernel, and whether it is READ is the `in_eos`
      !! gate's business, not the caller's.
      logical, intent(in) :: in_eos, wright
      real(wp), intent(in) :: p_top_val
      type(ocean_epbl_t), intent(inout) :: epbl
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: TAU = 0.1_wp
      integer :: k
      real(wp) :: t_k

      call make_grid(grid, 10, 8)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid)
      ss%rho0 = RHO0
      call sf%init(grid)
      call setup_epbl(epbl, grid)
      epbl%in_eos = in_eos
      if (wright) epbl%eos%variant = EOS_VARIANT_WRIGHT_97

      ms%h_layer = PTOP_HK
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      ms%p_top = p_top_val
      do k = 1, NZ
         t_k = 2.0_wp - 0.01_wp*real(NZ - k, wp)*PTOP_HK
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_k*PTOP_HK
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 34.5_wp*PTOP_HK
         ms%rho_layer(:, :, k) = RHO0
      end do
      call ss%set_wind_stress_const(TAU, 0.0_wp)

      !$acc enter data copyin(ms, epbl, ss, sf)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      !$acc update device(ms%p_top, ss%stress_mag, ss%stress_shelf)
      call epbl_compute(grid, epbl, ms, ss, PTOP_DT, sf=sf)
      !$acc update self(epbl%mld, epbl%kd_int, epbl%b0)
      !$acc update self(epbl%dpe_t%data, epbl%dcolht_t%data)
      call sf%exit_data()
      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, sf)

      call check(error, all(ieee_is_finite(epbl%kd_int)), &
                 "EPBL returned a non-finite kd_int")
      call sf%destroy(); call ss%destroy(); call ms%destroy()
   end subroutine epbl_ptop_run

   subroutine test_epbl_ptop_off(error)
      !! THE BIT-IDENTITY GATE, and it is not the trivial one.  A cavity
      !! fills `ms%p_top` with the ice load whether or not
      !! `&ocean_psurf_nml in_eos` is set, so "the array is zero" does
      !! NOT keep an existing cavity + EPBL run unchanged -- only the
      !! `in_eos` gate does.  So this compares `p_top = 3e6 Pa` with
      !! `in_eos = .false.` against `p_top = 0`, and demands BYTE
      !! equality: one code path either reads the seed or it does not,
      !! which is exactly the case where bit-equality is the right
      !! assertion.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_epbl_t) :: e_zero, e_loaded
      checks: block
         call epbl_ptop_run(.false., .false., 0.0_wp, e_zero, error)
         if (allocated(error)) exit checks
         call epbl_ptop_run(.false., .false., 3.0e6_wp, e_loaded, error)
         if (allocated(error)) exit checks

         call check(error, all(e_zero%kd_int == e_loaded%kd_int), &
                    "in_eos = .false. is not bit-identical under a loaded "// &
                    "column -- the gate is leaking p_top into the EPBL stack")
         if (allocated(error)) exit checks
         call check(error, all(e_zero%mld == e_loaded%mld), &
                    "in_eos = .false. moved the EPBL MLD under a loaded column")
         if (allocated(error)) exit checks
         call check(error, maxval(e_zero%kd_int) > 0.0_wp, &
                    "the reference run must actually mix, else the identity is "// &
                    "vacuous")
      end block checks
      call e_loaded%destroy(); call e_zero%destroy()
   end subroutine test_epbl_ptop_off

   subroutine test_epbl_ptop_monotone(error)
      !! The `p_top` seam contract's standing requirement for a joining
      !! builder: the pressure that reaches the EOS must be `>= p_top > 0`
      !! and must INCREASE toward the bed (`k = 1`).  MOM6 shipped a
      !! NEGATIVE EOS pressure in one path for years, which is why this is
      !! asserted by INVERTING the scheme's own output rather than by
      !! re-deriving the stack.
      !!
      !! The inversion is exact and duplicates no coefficient:
      !!
      !!     dpe_t    = dmass * p_mid * dSV/dT
      !!     dcolht_t = dmass *         dSV/dT
      !!   => p_mid   = dpe_t / dcolht_t.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_epbl_t) :: epbl
      real(wp), parameter :: P_TOP = 3.0e6_wp
      real(wp) :: p_mid(NZ)
      integer :: k, i_p, j_p
      checks: block
         call epbl_ptop_run(.true., .false., P_TOP, epbl, error)
         if (allocated(error)) exit checks

         i_p = size(epbl%mld, 1)/2
         j_p = size(epbl%mld, 2)/2
         do k = 1, NZ
            call check(error, abs(epbl%dcolht_t%data(i_p, j_p, k)) > 0.0_wp, &
                       "the steric weight is zero -- the inversion below would "// &
                       "divide by it")
            if (allocated(error)) exit checks
            p_mid(k) = epbl%dpe_t%data(i_p, j_p, k)/epbl%dcolht_t%data(i_p, j_p, k)
         end do

         do k = 1, NZ
            call check(error, p_mid(k) >= P_TOP, &
                       "an EOS pressure below the top-of-column load -- the stack "// &
                       "is not seeded at p_top")
            if (allocated(error)) exit checks
         end do
         do k = NZ - 1, 1, -1
            call check(error, p_mid(k) > p_mid(k + 1), &
                       "the EOS pressure must INCREASE toward the bed (k = 1)")
            if (allocated(error)) exit checks
         end do
         ! ... and the top layer's pressure is the load plus exactly half
         ! its own hydrostatic weight, which pins the seed itself.
         call check(error, abs(p_mid(NZ) - (P_TOP + 0.5_wp*GRAVITY*RHO0*PTOP_HK)) <= &
                    TOL_REL*p_mid(NZ), &
                    "the surface layer's EOS pressure is not p_top + 0.5*g*rho0*h")
      end block checks
      call epbl%destroy()
   end subroutine test_epbl_ptop_monotone

   subroutine test_epbl_ptop_gauge(error)
      !! THE PE-LEDGER DECISION, tested rather than asserted.
      !!
      !! The seed moves BOTH consumers of the column stack: the in-situ
      !! EOS argument and the PE weight `dpe = dmass*p_mid*dSV`.  That is
      !! deliberate -- the weight is the hydrostatic load a layer's centre
      !! of mass has to lift, and under a floating shelf the ice is part
      !! of that load, so it is the SAME pressure and splitting the two
      !! would put two conventions in one column.
      !!
      !! The risk that decision carries is that a UNIFORM load, which is
      !! pure gauge (it has no gradient and moves no dynamics), would
      !! nevertheless change the mixing energetics.  It must not, and
      !! with a LINEAR EOS -- whose `dSV/dT`, `dSV/dS` are CONSTANTS,
      !! independent of pressure -- the uniform offset `P` enters only
      !! through `dpe -> dpe + P*dcolht`, i.e. only through the column
      !! HEIGHT change of a mixing event.  A linear EOS conserves that
      !! exactly (mixing at fixed mass conserves `sum mass*T` and
      !! `sum mass*S`, and the height is a fixed linear functional of
      !! them), so the offset must cancel to round-off.
      !!
      !! Bound: the cancellation is between terms of size `P*dcolht`
      !! against a `pec_core` of size `p_mid*dcolht`, so the relative
      !! residual is `O(n_ops * eps * P/p_mid)`.  With `P = 3e6 Pa`,
      !! `p_mid ~ O(1e6 Pa)` and `O(100)` operations in the sweep that is
      !! `~ 1e-12` relative -- the bound below is `1e-9`, three decades of
      !! headroom, and still eleven decades below any physical effect.
      !!
      !! MEASURED, gfortran 15.1 Release, 2026-09-20: the difference is
      !! EXACTLY ZERO -- the cancellation is algebraic, not statistical,
      !! because `pec_core -> pec_core + P*colht_core` and `colht_core`
      !! is itself identically zero for constant `dSV/dX`.  The bound is
      !! kept rather than asserting `== 0` because an FMA-contracting
      !! build is free to associate the two expression trees differently.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_epbl_t) :: e_zero, e_loaded
      real(wp) :: scale, worst
      checks: block
         call epbl_ptop_run(.true., .false., 0.0_wp, e_zero, error)
         if (allocated(error)) exit checks
         call epbl_ptop_run(.true., .false., 3.0e6_wp, e_loaded, error)
         if (allocated(error)) exit checks

         scale = maxval(abs(e_zero%kd_int))
         call check(error, scale > 0.0_wp, &
                    "the reference run must actually mix, else this is vacuous")
         if (allocated(error)) exit checks
         worst = maxval(abs(e_loaded%kd_int - e_zero%kd_int))
         call check(error, worst <= 1.0e-9_wp*scale, &
                    "a UNIFORM p_top changed the EPBL diffusivity under a LINEAR "// &
                    "(pressure-independent) EOS by more than round-off -- the PE "// &
                    "ledger is not gauge-neutral")
         if (allocated(error)) exit checks
         worst = maxval(abs(e_loaded%mld - e_zero%mld))
         call check(error, worst <= 1.0e-9_wp*maxval(abs(e_zero%mld)), &
                    "a UNIFORM p_top moved the EPBL MLD under a LINEAR EOS")
      end block checks
      call e_loaded%destroy(); call e_zero%destroy()
   end subroutine test_epbl_ptop_gauge

   subroutine test_epbl_ptop_wright(error)
      !! The other half of the same statement, and what makes the port
      !! worth having: under the NONLINEAR Wright (1997) EOS the same
      !! uniform load DOES move the answer, because `dSV/dT` and `dSV/dS`
      !! are genuine functions of pressure there.  Three megapascals is
      !! 300 m of ice, and the thermal expansion of seawater changes by
      !! several percent over it.
      !!
      !! Without this the gauge test above could pass for the wrong
      !! reason -- a seed that never reached the EOS at all.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_epbl_t) :: e_zero, e_loaded
      real(wp) :: scale, worst
      checks: block
         call epbl_ptop_run(.true., .true., 0.0_wp, e_zero, error)
         if (allocated(error)) exit checks
         call epbl_ptop_run(.true., .true., 3.0e6_wp, e_loaded, error)
         if (allocated(error)) exit checks

         scale = maxval(abs(e_zero%kd_int))
         call check(error, scale > 0.0_wp, "the Wright reference run must mix")
         if (allocated(error)) exit checks
         worst = maxval(abs(e_loaded%kd_int - e_zero%kd_int))
         call check(error, worst > 1.0e-6_wp*scale, &
                    "a 3 MPa load did NOT move the nonlinear-EOS answer -- the "// &
                    "seed is not reaching eos_specvol_derivs at all")
      end block checks
      call e_loaded%destroy(); call e_zero%destroy()
   end subroutine test_epbl_ptop_wright

end module test_ocean_bl_under_ice
