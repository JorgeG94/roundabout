!! THE P5.2 PHYSICS GATE: a stratified cavity under a SLOPING ice shelf,
!! initialised at rest, stays at rest.
module test_ocean_cavity_load
   !! ### What this slice wired, and what has to be true afterwards
   !!
   !! P5.2 completes the load partition:
   !!
   !! ```
   !! ms%p_top   = metrics%p_ice_ref + sf%p_surf     (the pressure)
   !! bt_H_ref   = b - z_draft                       (the barotropic datum)
   !! eta_ib     = -sf%p_surf/(rho0*g_bt)            (the seam: ANOMALY only)
   !! ```
   !!
   !! so `pa(nz+1) = rho_ref*g*eta_geo + p_top` is `0` at rest (the two
   !! terms are the same product with opposite signs), the barotropic mode
   !! sees `bt_eta = 0` and a seam of zero, and the only thing left that
   !! can move the water is the PGF's own TRUNCATION ERROR.
   !!
   !! ### The residual, derived rather than tolerated
   !!
   !! Writing `Pi(z) = P(z) + rho_ref*g*z` for the pressure anomaly
   !! potential, the FV_MOM6 Pass-3 numerator telescopes exactly:
   !!
   !! ```
   !!   numer(k) = G(k) - G(k+1),
   !!   G(K) = integral_{e_L(K)}^{e_R(K)} Pi dz  -  De(K)*intx_pa(K)
   !! ```
   !!
   !! i.e. `G(K)` is the TRAPEZOID error of the face pressure over the
   !! vertical gap `De(K) = e_R(K) - e_L(K)` between the two columns'
   !! K-th interfaces.  It vanishes iff `Pi` is linear there, i.e. iff the
   !! density is constant across the gap; for smooth density
   !! `G(K) = -(De(K)^3/12)*Pi'' = -(De(K)^3/12)*rho_0*N^2`.
   !!
   !! Flat bed + `VCOORD_SIGMA` makes `De(K) = sigma_K*D` with
   !! `D = z_draft,R - z_draft,L = slope*dx`, and what survives the split
   !! solver's depth-mean replacement is
   !!
   !! ```
   !!   PFu(k) - <PFu>_h  =  N^2 * D^3 * (3*sigma_k^2 - 1) / (12*dx*Hbar)
   !! ```
   !!
   !! — SECOND order in `dx`, CUBIC in the draft slope, independent of
   !! `nz`, peaking at the ice base (`sigma = 1`) at
   !! `a_peak = N^2*D^3/(6*dx*Hbar)`.  On an f-plane at `f = 0` the
   !! spurious velocity is then bounded by `a_peak * t`, and THAT is the
   !! bound asserted below.  It is a formula, so a future PGF change that
   !! alters the order of accuracy is caught rather than absorbed.
   !!
   !! ### Measured, gfortran 15.1 Release, 2026-09-20
   !!
   !! ```
   !!   END TO END (40 steps x 300 s = 12000 s, f = 0, no forcing)
   !!     stratified sloping lid   max|u| = 7.083e-08 m/s
   !!                              max|v| = 1.17e-31 m/s  (draft varies in x only)
   !!                              A_PEAK (derived)        8.681e-12 m/s^2
   !!                              A_PEAK*T_TOTAL          1.042e-07 m/s
   !!                              bound (x4)              4.166e-07   -> 5.9x margin
   !!     uniform density          max|u| = 2.578e-13 m/s  (5 decades below)
   !!
   !!   AT THE PGF (unit, D = 2 m, dx = 1000 m)
   !!     residual deviation       -1.1559e-11 m/s^2
   !!     formula prediction        1.1421e-11 m/s^2      -> 1.2 % agreement
   !!     raw |PFu|, load OFF       1.9619e-02 m/s^2  = g*slope (0.03 % off)
   !!     raw |PFu|, load ON        5.8128e-06 m/s^2  = N^2*z_draft*slope
   !! ```
   !!
   !! For scale: the (mis-attributed, see below) literature figure the
   !! Phase-5 design quotes for quiet linear-stratification cavity cases
   !! WITH the sloping-coordinate PGF corrections is `O(1e-9 m/s)`; this
   !! build, with NO such correction, sits at `7.1e-8 m/s` over 12000 s of
   !! a 1e-3 draft slope.  The correction family (top mass-weighting,
   !! reference-interface reset, flattest-interface fallback) is the NEXT
   !! phase's work, and this number is its baseline.  NOTE the design's
   !! own citation correction: `papers/19b_yung2026_isomip2_results.pdf`
   !! is Yung, Asay-Davis, Adcroft et al. (2026), *Results of the second
   !! Ice Shelf-Ocean Model Intercomparison Project*, The Cryosphere 20,
   !! 2053-2088 — an intercomparison with NO PGF equations in it.  The
   !! `1e-9` figure is quoted here because the task asked for it recorded
   !! next to the measurement, not because this file could verify it.
   !!
   !! ### The uniform-density case
   !!
   !! With `N^2 = 0` every `G(K)` is identically zero, so the truncation
   !! residual is not merely small, it is absent: what is left is the
   !! rounding of the geopotential interface stack that closes
   !! `pa(nz+1)`.  That is the case which isolates the LOAD, and it lands
   !! at `2.6e-13 m/s` — five decades under the stratified case, which is
   !! the statement that the O(D^3) term really is what the stratified
   !! number is made of.
   !!
   !! ### The control, and why it lives at the PGF and not end to end
   !!
   !! `load_off_control_at_the_pgf` evaluates the SAME resting state twice,
   !! with `p_top_in_bc` on and off.  With the load off the raw face force
   !! is `g*grad(z_draft) = 1.96e-2 m/s^2` — the design's
   !! `Delta u ~ -g*grad(d)*dt` prediction, exactly, and what the unsplit
   !! driver or any consumer of the raw `pa` stack would feel.  With it on
   !! the raw force falls by `g/(N^2 z_draft) ~ 3.4e3` to `5.81e-6 m/s^2`
   !! — NOT to the truncation floor, and deliberately so: the
   !! Boussinesq-isostatic load is the displaced weight at the REFERENCE
   !! density, and against a stratified column it is short by
   !! `-g*integral(rho_hat - rho_0)`, whose gradient is `N^2 z_draft s`.
   !! The design picks that load anyway (its §3.1, R1 vs R2) because it is
   !! the one that makes the DISCRETE BAROTROPIC state exactly at rest,
   !! and ISOMIP+ prescribes the same constant-density form.
   !!
   !! The same test then asserts the half that makes all of the above
   !! safe: the DEVIATION of `PFu` from its thickness-weighted depth mean
   !! is the same with the load on and off, to 1e-10 relative.  A
   !! depth-uniform top load is baroclinically inert (the theorem in
   !! `compute_fv_mom6_impl`'s docstring) and the split solver replaces
   !! the depth mean with the barotropic solution, so BOTH the missing
   !! load (control (a)) and the load's own `rho_hat` approximation
   !! (control (b)) are annihilated there.  That is why the end-to-end
   !! gate sees `7e-8 m/s` and not `5.8e-6 * 12000 = 7e-2 m/s`, and why
   !! an end-to-end "load off" control would NOT blow up either: the
   !! split path is protected by construction and would only lose 5.4
   !! decades of conditioning (which is what
   !! `test_ocean_cavity_equivalence` measures).  The configure-time
   !! refusal therefore protects the raw stack, the unsplit driver and
   !! `&ocean_bt_nml correction_h_weighted`, not the default path's
   !! stability.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, rdb_ocean_get_b_ptr, &
                            rdb_ocean_get_grid_info, rdb_ocean_set_tracer
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_cavity_load_tests

   ! ---- End-to-end case geometry -------------------------------------
   integer, parameter :: NX_PHYS = 32, NZ_ML = 10
   real(wp), parameter :: DXY = 2000.0_wp
   real(wp), parameter :: BED = 1000.0_wp
      !! Flat bed (m).  Flat on purpose: it puts the WHOLE interface tilt
      !! into the ice base, which is what the residual formula assumes.
   real(wp), parameter :: DRAFT_X0 = -10000.0_wp
      !! Anchor + western edge of the shelf box (m).  Placed well WEST of
      !! the ghost band so the shelf covers the ghosts continuously —
      !! a box edge inside the domain would put a calving-front STEP in
      !! `z_draft` and the residual formula (which assumes a smooth
      !! slope) would no longer describe it.
   real(wp), parameter :: DRAFT_AT_X0 = 190.0_wp
   real(wp), parameter :: DRAFT_SLOPE = 1.0e-3_wp
      !! Ice base deepens eastward: 200 m at x = 0 to 264 m at x = 64 km.
   real(wp), parameter :: N2_TARGET = 1.0e-5_wp
      !! Buoyancy frequency squared (1/s^2) — the stratification the
      !! geopotential T(z) profile below is built to produce.
   real(wp), parameter :: ALPHA_T = 0.2_wp
      !! `&ocean_ic_nml alpha_T`, DIMENSIONAL (kg/m^3 per degC).  Set
      !! explicitly: the shipped default (1.7e-4) is three decades below
      !! a physical thermal expansion and would need a 6 K/m lapse rate
      !! to reach N^2 = 1e-5.
   real(wp), parameter :: RHO_0 = 1035.0_wp
   real(wp), parameter :: T_SURF = 10.0_wp
      !! T at z = 0 (degC) — also `&ocean_ic_nml T_ref`, so the surface
      !! density anomaly is zero.
   real(wp), parameter :: DTDZ = N2_TARGET*RHO_0/(GRAVITY*ALPHA_T)
      !! degC per metre of GEOPOTENTIAL height, from
      !! `N^2 = (g/rho_0)*alpha_T*dT/dz` for the linear EOS.

   integer, parameter :: N_STEPS = 40
   real(wp), parameter :: DT = 300.0_wp
   real(wp), parameter :: T_TOTAL = DT*real(N_STEPS, wp)

   real(wp), parameter :: D_FACE = DRAFT_SLOPE*DXY
      !! Draft step between adjacent columns (m) — the `D` of `G(K)`.
   real(wp), parameter :: H_BAR = BED - (DRAFT_AT_X0 + DRAFT_SLOPE* &
                                         (-DRAFT_X0 + 0.5_wp*real(NX_PHYS, wp)*DXY))
      !! Mean water-column thickness (m) at the domain centre.
   real(wp), parameter :: A_PEAK = N2_TARGET*D_FACE**3/(6.0_wp*DXY*H_BAR)
      !! Peak surviving baroclinic acceleration (m/s^2):
      !! `N^2*D^3*(3*sigma^2-1)/(12*dx*Hbar)` at `sigma = 1`.
   real(wp), parameter :: REST_SAFETY = 4.0_wp
      !! Covers the barotropic adjustment the residual drives on its way
      !! to a balance and the ALE remap's own second-order error.  It is
      !! NOT fitted: the integrated residual `A_PEAK*T_TOTAL` is
      !! `1.04e-7 m/s` and the measurement is `7.08e-8` — 0.68 of the raw
      !! prediction before any safety factor at all.

   ! ---- PGF-level (unit) case ----------------------------------------
   integer, parameter :: NG_U = 2, NZ_U = 4
   real(wp), parameter :: DX_U = 1000.0_wp
   real(wp), parameter :: RHO_REF_U = 1035.0_wp
   real(wp), parameter :: BED_U = 1024.0_wp
   real(wp), parameter :: DRAFT0_U = 256.0_wp
   real(wp), parameter :: DRAFT_STEP_U = 2.0_wp
      !! Per-column draft step (m) for the sloping-lid unit case — `D`.

contains

   subroutine collect_ocean_cavity_load_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_resting_sloping_lid_stratified", test_rest_stratified), &
                  new_unittest("cavity_resting_uniform_density", test_rest_uniform), &
                  new_unittest("flat_lid_stratified_pfu_is_zero", test_flat_lid_pfu), &
                  new_unittest("sloping_lid_residual_matches_formula", test_slope_formula), &
                  new_unittest("load_off_control_at_the_pgf", test_load_off_control) &
                  ]
   end subroutine collect_ocean_cavity_load_tests

   ! ==================================================================
   ! End-to-end: the resting loaded cavity
   ! ==================================================================

   function rest_nml() result(nml)
      !! The resting cavity namelist.  No wind, no heat/salt flux, `f = 0`
      !! (so a spurious acceleration integrates cleanly to `a*t` instead of
      !! turning into a geostrophic balance whose amplitude is `a/f`), a
      !! flat bed, a linear ice draft, sigma, `pred_corr`, a pinned
      !! `n_inner`, and the load wired all the way through.
      character(len=:), allocatable :: nml
      character(len=32) :: t_bot_s

      ! Ghost columns keep whatever the IC seeds (the interior T is
      ! overwritten from GEOPOTENTIAL depth afterwards, and the API's
      ! tracer setter writes the interior only).  With WALL boundaries the
      ! ghosts reach nothing: the physical boundary faces are hard-zeroed
      ! by the wall closure, and the first interior cell's continuity uses
      ! that same zero flux.  Seeding the IC with the SURFACE value at
      ! both ends keeps the ghost columns unstratified and the mismatch
      ! trivially small anyway.
      write (t_bot_s, '(F12.4)') T_SURF
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 32, ny = 6, nghost = 2, dx = 2000.0, dy = 2000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 10 /"//new_line("a")// &
            "&time_nml t_end = 1.0e7, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0, taux_magnitude = 0.0 /"//new_line("a")// &
            "&physics_nml coriolis_f = 0.0 /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, T_init_bottom = "// &
            trim(adjustl(t_bot_s))//", T_init_surface = "// &
            trim(adjustl(t_bot_s))//" /"//new_line("a")// &
            "&ocean_ic_nml alpha_T = 0.2, T_ref = 10.0, S_ref = 35.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .false., "// &
            "n_inner = 24 /"//new_line("a")// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 190.0, draft_slope = 1.0e-3, draft_x0 = -10000.0 /"// &
            new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function rest_nml

   subroutine run_rest_case(stratified, umax, vmax, ok)
      !! Create the resting cavity, optionally impose FLAT ISOPYCNALS by
      !! setting the interior temperature from the layer centre's
      !! GEOPOTENTIAL height (not from its layer index — under a sloping
      !! draft the sigma layers tilt with the ice base, so a profile laid
      !! out per index would tilt the isopycnals with them and the run
      !! would not be at rest for a physical reason), step, and report the
      !! largest interior face speed.
      logical, intent(in) :: stratified
      real(wp), intent(out) :: umax, vmax
      logical, intent(out) :: ok

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen, nxp, nyp, nzp, ngc
      real(wp), pointer :: h3(:, :, :), b2(:, :), u3(:, :, :), v3(:, :, :)
      real(wp), allocatable :: tval(:, :, :)
      character(len=:), allocatable :: nml
      integer :: i, j, k, ng
      real(wp) :: e_low, z_ctr

      ok = .false.
      umax = 0.0_wp
      vmax = 0.0_wp
      handle = c_null_ptr
      nml = rest_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return

      if (stratified) then
         status = rdb_ocean_refresh_host(handle)
         if (status /= OCEAN_STATUS_OK) then
            status = rdb_ocean_destroy(handle); return
         end if
         status = rdb_ocean_get_grid_info(handle, nxp, nyp, nzp, ngc)
         if (status /= OCEAN_STATUS_OK) then
            status = rdb_ocean_destroy(handle); return
         end if
         ng = int(ngc)
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h3, [int(nx), int(ny), int(nz)])
         status = rdb_ocean_get_b_ptr(handle, ptr, nx, ny, gen)
         call c_f_pointer(ptr, b2, [int(nx), int(ny)])

         allocate (tval(int(nxp), int(nyp), int(nzp)))
         do j = 1, int(nyp)
            do i = 1, int(nxp)
               ! Bottom-up stack: e_face(1) = -b, e_face(k+1) = e_face(k) + h(k).
               e_low = -b2(ng + i, ng + j)
               do k = 1, int(nzp)
                  z_ctr = e_low + 0.5_wp*h3(ng + i, ng + j, k)
                  tval(i, j, k) = T_SURF + DTDZ*z_ctr
                  e_low = e_low + h3(ng + i, ng + j, k)
               end do
            end do
         end do
         status = rdb_ocean_set_tracer(handle, "temperature", 11_c_int, tval, &
                                       nxp, nyp, nzp)
         if (status /= OCEAN_STATUS_OK) then
            status = rdb_ocean_destroy(handle); return
         end if
      end if

      status = rdb_ocean_step(handle, int(N_STEPS, c_int))
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_refresh_host(handle)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_get_grid_info(handle, nxp, nyp, nzp, ngc)
      ng = int(ngc)
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, u3, [int(nx), int(ny), int(nz)])
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, v3, [int(nx), int(ny), int(nz)])
      ! STRICTLY interior faces: the wall faces are hard-zeroed and the
      ! ghost band is not a solution anywhere.
      umax = maxval(abs(u3(ng + 2:ng + int(nxp), ng + 1:ng + int(nyp), :)))
      vmax = maxval(abs(v3(ng + 1:ng + int(nxp), ng + 2:ng + int(nyp), :)))
      ok = all(ieee_is_finite(u3)) .and. all(ieee_is_finite(v3))

      status = rdb_ocean_destroy(handle)
      ok = ok .and. (status == OCEAN_STATUS_OK)
   end subroutine run_rest_case

   subroutine test_rest_stratified(error)
      !! THE GATE.  Sloping draft, flat bed, flat isopycnals, at rest, no
      !! forcing: the spurious velocity after `N_STEPS` must stay under
      !! the sigma-PGF truncation bound `REST_SAFETY * A_PEAK * T_TOTAL`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: umax, vmax, bound
      logical :: ok

      call run_rest_case(.true., umax, vmax, ok)
      call check(error, ok, "the resting stratified cavity must run and stay finite")
      if (allocated(error)) return

      bound = REST_SAFETY*A_PEAK*T_TOTAL
      call check(error, umax <= bound, &
                 "the spurious zonal velocity must stay under the derived "// &
                 "sigma-coordinate PGF truncation bound N^2 D^3/(6 dx Hbar) * t")
      if (allocated(error)) return
      call check(error, vmax <= bound, &
                 "the spurious meridional velocity must stay under the same bound")
      if (allocated(error)) return
      ! Non-vacuity in the OTHER direction: the residual is a real,
      ! identified truncation error, not an accidental exact zero.  If it
      ! ever became zero the formula above would no longer be what is
      ! being tested and the bound would stop meaning anything.
      call check(error, umax > 0.0_wp, &
                 "the truncation residual must be present (a bit-zero result here "// &
                 "means the case stopped exercising the sloping-coordinate PGF)")
   end subroutine test_rest_stratified

   subroutine test_rest_uniform(error)
      !! `N^2 = 0` kills every `G(K)`, so the truncation residual is not
      !! small but ABSENT, and what is left is the rounding of the
      !! `pa(nz+1) = rho_ref*g*(-z_draft) + p_ice_ref` cancellation.  This
      !! is the case that isolates the LOAD: it is the same sloping ice
      !! base, and without the load term in the top BC the raw stack would
      !! carry `g*grad(z_draft)`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: umax, vmax, bound
      logical :: ok

      call run_rest_case(.false., umax, vmax, ok)
      call check(error, ok, "the resting uniform-density cavity must run and stay finite")
      if (allocated(error)) return

      ! Round-off floor.  `pa(nz+1) = rho_ref*g*eta_geo + p_ice_ref`
      ! cancels EXACTLY only if `eta_geo` lands on `-z_draft` exactly; it
      ! is accumulated as `-b + sum_k h(k)`, so it carries `~sqrt(nz)`
      ! roundings at the ulp of the BED depth.  The surface BC turns that
      ! into `eps*rho_0*g*BED*sqrt(nz)` Pa, and the Pass-3 face force
      ! turns a `dp` into `dp/(rho_0*dx)` of acceleration, integrated
      ! over `T_TOTAL`.
      bound = REST_SAFETY*(epsilon(1.0_wp)*RHO_0*GRAVITY*BED* &
                           sqrt(real(NZ_ML, wp)))*T_TOTAL/(RHO_0*DXY)
      call check(error, umax <= bound, &
                 "with uniform density the resting cavity must be at the round-off "// &
                 "floor of the pa(nz+1) cancellation")
      if (allocated(error)) return
      call check(error, vmax <= bound, "and the same meridionally")
      if (allocated(error)) return
      ! The point of the case: with the stratification removed the
      ! residual is not merely under a bound, it is DECADES under the
      ! stratified case's own O(D^3) prediction.  That is what says the
      ! stratified number is the truncation term and not something else.
      call check(error, umax <= 1.0e-3_wp*A_PEAK*T_TOTAL, &
                 "removing the stratification must remove the residual, not just "// &
                 "shrink it — N^2 = 0 makes every trapezoid error G(K) vanish")
   end subroutine test_rest_uniform

   ! ==================================================================
   ! PGF-level unit cases
   ! ==================================================================

   subroutine make_grid_u(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NG_U, DX_U, DX_U)
   end subroutine make_grid_u

   subroutine make_pgf_u(grid, pgf, p_top_in_bc)
      type(hgrid_t), intent(in) :: grid
      type(ocean_pressure_force_t), intent(out) :: pgf
      logical, intent(in) :: p_top_in_bc
      call pgf%init(grid, nz_ml=NZ_U)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO_REF_U
      pgf%rho_ref = RHO_REF_U
      pgf%p_top_in_bc = p_top_in_bc
   end subroutine make_pgf_u

   subroutine run_pgf_u(grid, ms, pgf, b)
      !! One FV_MOM6 pass.  `mem:separate` discipline: `ms`, its `p_top`
      !! companion and `pgf` are all mapped before the kernel; `p_top`,
      !! `h_layer` and `rho_layer` are HOST-set inputs so they are pushed
      !! explicitly after `enter_data`; the read-back uses `update self`
      !! on the COMPONENT arrays, never the aggregate derived type.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: b(:, :)
      type(ocean_metrics_t) :: metrics

      call make_cartesian_metrics(metrics, grid)
      call pgf%set_bathymetry(b)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      !$acc update device(ms%p_top, ms%h_layer, ms%rho_layer)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data, pgf%pa%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf_u

   pure subroutine seed_cavity_column(ms, b, draft_step)
      !! Flat bed, a draft stepping by `draft_step` per column, sigma
      !! layers filling `b - z_draft`, a LINEAR-IN-GEOPOTENTIAL-z density
      !! (so `N^2` is a constant) and the matching isostatic load.  This
      !! is the discrete resting state.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(out) :: b(:, :)
      real(wp), intent(in) :: draft_step
      integer :: i, j, k, nx, ny
      real(wp) :: z_draft, water, e_low, z_ctr, drho_dz

      nx = size(ms%p_top, 1)
      ny = size(ms%p_top, 2)
      ! dRho/dz for N^2 = N2_TARGET: N^2 = -(g/rho_0) dRho/dz.
      drho_dz = -N2_TARGET*RHO_REF_U/GRAVITY
      do j = 1, ny
         do i = 1, nx
            b(i, j) = BED_U
            z_draft = DRAFT0_U + draft_step*real(i - 1, wp)
            water = BED_U - z_draft
            e_low = -BED_U
            do k = 1, NZ_U
               ms%h_layer(i, j, k) = water/real(NZ_U, wp)
               z_ctr = e_low + 0.5_wp*ms%h_layer(i, j, k)
               ms%rho_layer(i, j, k) = RHO_REF_U + drho_dz*z_ctr
               e_low = e_low + ms%h_layer(i, j, k)
            end do
            ms%p_top(i, j) = (RHO_REF_U*GRAVITY)*z_draft
         end do
      end do
   end subroutine seed_cavity_column

   subroutine test_flat_lid_pfu(error)
      !! A FLAT lid over a flat bed with a real stratification: every
      !! `De(K)` is zero, so every `G(K)` is zero and the face force is
      !! exactly zero — no truncation term to bound, structurally.  This
      !! is the `D -> 0` end of the formula the next case checks.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)

      checks: block
         call make_grid_u(grid, 10, 8)
         ms%nz_ml = NZ_U
         call ms%init(grid)
         allocate (b(grid%nx_total, grid%ny_total))
         call seed_cavity_column(ms, b, 0.0_wp)
         call make_pgf_u(grid, pgf, .true.)
         call run_pgf_u(grid, ms, pgf, b)

         call check(error, maxval(abs(pgf%dpdx_face%data)) == 0.0_wp, &
                    "a FLAT loaded lid over flat bed with N^2 > 0 must give an "// &
                    "exactly zero zonal face force (every interface gap is zero)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(pgf%dpdy_face%data)) == 0.0_wp, &
                    "and an exactly zero meridional face force")
      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_flat_lid_pfu

   subroutine test_slope_formula(error)
      !! The residual FORMULA, not just a tolerance.  With a linear
      !! density profile, a flat bed and sigma layers, what survives the
      !! depth-mean subtraction must be
      !! `N^2 D^3 (3 sigma_k^2 - 1)/(12 dx Hbar)` at every layer.  Checked
      !! at the peak (`k = NZ`, `sigma ~ 1`) to 20 %: the discrete
      !! `sigma_k` is the layer CENTRE, `(k - 0.5)/nz`, not 1, so the
      !! comparison uses the centre value and the slack covers the
      !! difference between the continuum integral and the discrete stack.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp) :: dev(NZ_U), wsum, hsum, pred, hbar, sig, measured
      integer :: i, j, k, nx, ny

      checks: block
         call make_grid_u(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ_U
         call ms%init(grid)
         allocate (b(nx, ny))
         call seed_cavity_column(ms, b, DRAFT_STEP_U)
         call make_pgf_u(grid, pgf, .true.)
         call run_pgf_u(grid, ms, pgf, b)

         ! One interior face, mid-domain.
         i = nx/2
         j = ny/2
         wsum = 0.0_wp
         hsum = 0.0_wp
         do k = 1, NZ_U
            hbar = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
            wsum = wsum + hbar*pgf%dpdx_face%data(i, j, k)
            hsum = hsum + hbar
         end do
         do k = 1, NZ_U
            dev(k) = pgf%dpdx_face%data(i, j, k) - wsum/hsum
         end do

         sig = (real(NZ_U, wp) - 0.5_wp)/real(NZ_U, wp)
         pred = N2_TARGET*DRAFT_STEP_U**3*(3.0_wp*sig*sig - 1.0_wp)/ &
                (12.0_wp*DX_U*hsum)
         measured = dev(NZ_U)

         call check(error, abs(measured) > 0.0_wp, &
                    "the sloping-lid residual must be non-zero (else the formula "// &
                    "check is vacuous)")
         if (allocated(error)) exit checks
         call check(error, abs(abs(measured) - abs(pred)) <= 0.2_wp*abs(pred), &
                    "the surviving baroclinic residual must match "// &
                    "N^2 D^3 (3 sigma^2 - 1)/(12 dx Hbar) to 20 % — a formula "// &
                    "check, so a change in the PGF's order of accuracy is caught")
      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_slope_formula

   subroutine test_load_off_control(error)
      !! THE CONTROL.  Same resting state, `p_top_in_bc` toggled.  Both
      !! sides are checked against a FORMULA, not against each other.
      !!
      !!   (a) LOAD OFF — the raw face force is the whole unbalanced ice
      !!       base, `g*|grad z_draft| = g*s`.  This is the design's
      !!       `Delta u ~ -g grad(d) dt` per step, and it is what the
      !!       UNSPLIT driver, or any consumer that reads the raw `pa`
      !!       stack, would feel.
      !!   (b) LOAD ON — the raw force collapses by `g/(N^2 z_draft)`, to
      !!       `N^2 * z_draft * s`.  It does NOT collapse to the
      !!       truncation residual, and that is not a defect: the
      !!       Boussinesq-isostatic load `rho_0 g z_draft` is the weight
      !!       of the displaced water AT THE REFERENCE DENSITY, while the
      !!       real column is stratified, so the exact-analytic-rest load
      !!       `g*integral(rho_hat)` differs by `-g*integral(rho_hat -
      !!       rho_0)` and its gradient is precisely `N^2 z_draft s`.
      !!       The design (§3.1, its equations R1 vs R2) chooses
      !!       `rho_0 g z_draft` deliberately, because THAT is the load
      !!       that makes the DISCRETE BAROTROPIC state exactly at rest,
      !!       and ISOMIP+ prescribes the same constant-density form.
      !!   (c) The residual in (b) is DEPTH-UNIFORM — it enters through
      !!       `pa(nz+1)` and nowhere else — so the split solver's
      !!       depth-mean replacement annihilates it, which is why the
      !!       end-to-end gate above measures `7e-8 m/s` (the O(D^3)
      !!       truncation) rather than `5.8e-6 * 12000 = 7e-2 m/s`.  This
      !!       test asserts that inertness directly: the DEVIATION of
      !!       `PFu` from its thickness-weighted depth mean is the same
      !!       with the load on and off to 1e-10 relative.
      !!
      !! (c) is also why an end-to-end "load off" control would NOT blow
      !! up — the split path is protected by construction; it would only
      !! lose 5.4 decades of conditioning (that is what
      !! `test_ocean_cavity_equivalence` measures).  The configure-time
      !! refusal therefore protects the raw stack, the unsplit driver and
      !! `&ocean_bt_nml correction_h_weighted`, not the default path's
      !! stability.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_on, pgf_off
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf_on(:, :, :)
      real(wp) :: raw_on, raw_off, analytic, dev_on, dev_off, hbar
      real(wp) :: wsum_on, wsum_off, hsum, worst_rel, scale
      real(wp) :: mismatch, d_face
      integer :: i, j, k, nx, ny

      checks: block
         call make_grid_u(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ_U
         call ms%init(grid)
         allocate (b(nx, ny))
         call seed_cavity_column(ms, b, DRAFT_STEP_U)

         call make_pgf_u(grid, pgf_on, .true.)
         call run_pgf_u(grid, ms, pgf_on, b)
         allocate (pf_on, source=pgf_on%dpdx_face%data)
         call make_pgf_u(grid, pgf_off, .false.)
         call run_pgf_u(grid, ms, pgf_off, b)

         raw_on = maxval(abs(pf_on(3:nx - 1, 3:ny - 1, :)))
         raw_off = maxval(abs(pgf_off%dpdx_face%data(3:nx - 1, 3:ny - 1, :)))
         ! `-(1/rho_0) d(rho_ref g z_draft)/dx` = `-g * slope`.
         analytic = GRAVITY*DRAFT_STEP_U/DX_U

         call check(error, abs(raw_off - analytic) <= 0.02_wp*analytic, &
                    "with the load OFF the raw face force must be g*grad(z_draft) — "// &
                    "the unbalanced term the top BC exists to cancel")
         if (allocated(error)) exit checks

         ! (b) the load ON leaves exactly the rho_hat-vs-rho_0 mismatch,
         ! `N^2 * z_draft * s`, evaluated at the deepest draft in the
         ! sampled band (the face between the last two sampled columns).
         d_face = DRAFT0_U + DRAFT_STEP_U*(real(nx - 1, wp) - 1.5_wp)
         mismatch = N2_TARGET*d_face*(DRAFT_STEP_U/DX_U)
         call check(error, abs(raw_on - mismatch) <= 0.15_wp*mismatch, &
                    "with the load ON the raw face force must be exactly the "// &
                    "Boussinesq-isostatic load's own approximation error, "// &
                    "N^2*z_draft*slope — not zero, and not the truncation term")
         if (allocated(error)) exit checks
         call check(error, raw_on <= 1.0e-3_wp*raw_off, &
                    "that is a collapse of g/(N^2 z_draft) ~ 3.5e3 against the "// &
                    "unloaded stack")
         if (allocated(error)) exit checks

         ! (b) the depth-mean DEVIATION is the same either way.
         worst_rel = 0.0_wp
         scale = 0.0_wp
         do j = 3, ny - 1
            do i = 3, nx - 1
               wsum_on = 0.0_wp; wsum_off = 0.0_wp; hsum = 0.0_wp
               do k = 1, NZ_U
                  hbar = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
                  wsum_on = wsum_on + hbar*pf_on(i, j, k)
                  wsum_off = wsum_off + hbar*pgf_off%dpdx_face%data(i, j, k)
                  hsum = hsum + hbar
               end do
               do k = 1, NZ_U
                  dev_on = pf_on(i, j, k) - wsum_on/hsum
                  dev_off = pgf_off%dpdx_face%data(i, j, k) - wsum_off/hsum
                  worst_rel = max(worst_rel, abs(dev_on - dev_off))
                  scale = max(scale, abs(dev_on))
               end do
            end do
         end do
         call check(error, scale > 0.0_wp, &
                    "the baroclinic deviation must be non-zero (else (b) is vacuous)")
         if (allocated(error)) exit checks
         call check(error, worst_rel <= 1.0e-10_wp*max(scale, raw_off), &
                    "the load must be baroclinically INERT: the deviation from the "// &
                    "depth mean must be the same with it on and off")
      end block checks
      call pgf_on%destroy(); call pgf_off%destroy(); call ms%destroy()
   end subroutine test_load_off_control

end module test_ocean_cavity_load
