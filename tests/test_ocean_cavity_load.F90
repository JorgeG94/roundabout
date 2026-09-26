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
   !!
   !! ### Under the MOM6 split (`&ocean_bt_nml bc_pgf_forcing`, default)
   !!
   !! Everything in the previous paragraph that says "the split solver
   !! replaces the depth mean" describes the LEGACY split
   !! (`bc_pgf_forcing = .false.`).  The default now forces the barotropic
   !! mode with the depth mean of the full layer PGF, as MOM6 does, so the
   !! load's reference-density shortfall (control (b), `N^2 z_draft s`)
   !! is no longer annihilated: it is a real depth-mean bottom-pressure
   !! gradient the initial state does not balance, and the barotropic mode
   !! adjusts to it (`cavity_sloping_lid_load_shortfall_drives_bt`,
   !! 4.5e-4 m/s).  The truncation gate `cavity_resting_sloping_lid_
   !! stratified` therefore pins the legacy split — the quantity its
   !! formula describes is the deviation from the depth mean.
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
   use rdb_ocean_cavity, only: cavity_trim_eta_linear_impl
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

   ! ---- Trimmed-IC (MOM6 TRIM_IC_FOR_P_SURF) unit case ----------------
   ! The `cavity_sloping_lid_rest` geometry and ISOMIP+ COLD density on
   ! the linear EOS: rho(z) = RHO_S_T + DRHO_DZ_T*z (z positive UP), with
   ! the surface value BELOW rho_ref (the T/S references are not the
   ! surface values), which is what makes the load shortfall non-trivial.
   real(wp), parameter :: RHO_REF_T = 1027.51_wp
   real(wp), parameter :: RHO_S_T = RHO_REF_T + 8.0587609e-1_wp*(33.8_wp - 34.2_wp) &
                          - 3.8356948e-2_wp*(-1.9_wp + 1.0_wp)
      !! `rho_0 + beta_S*(S_s - S_ref) - alpha_T*(T_s - T_ref)`.
   real(wp), parameter :: DRHO_DZ_T = 8.0587609e-1_wp*(-1.0416667e-3_wp)
      !! `beta_S * dS/dz` (kg/m^4); N^2 = -(g/rho_0)*DRHO_DZ_T = 8.0e-6.
   real(wp), parameter :: DX_T = 2000.0_wp
   real(wp), parameter :: BED_T = 720.0_wp
   real(wp), parameter :: DRAFT0_T = 494.0_wp
   real(wp), parameter :: DRAFT_STEP_T = -13.8_wp
      !! The case's own per-face draft step (6.9e-3 * 2 km), shallowing.
   integer, parameter :: NZ_T = 15
      !! The case's own layer count.

contains

   subroutine collect_ocean_cavity_load_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_resting_sloping_lid_stratified", test_rest_stratified), &
                  new_unittest("cavity_sloping_lid_load_shortfall_drives_bt", test_load_shortfall_bt), &
                  new_unittest("cavity_resting_uniform_density", test_rest_uniform), &
                  new_unittest("flat_lid_stratified_pfu_is_zero", test_flat_lid_pfu), &
                  new_unittest("sloping_lid_residual_matches_formula", test_slope_formula), &
                  new_unittest("load_off_control_at_the_pgf", test_load_off_control), &
                  new_unittest("trim_ic_root_is_the_displaced_weight", test_trim_root), &
                  new_unittest("trim_ic_balances_the_depth_mean_pfu", test_trim_balances_pfu) &
                  ]
   end subroutine collect_ocean_cavity_load_tests

   ! ==================================================================
   ! End-to-end: the resting loaded cavity
   ! ==================================================================

   function rest_nml(legacy_split) result(nml)
      !! The resting cavity namelist.  No wind, no heat/salt flux, `f = 0`
      !! (so a spurious acceleration integrates cleanly to `a*t` instead of
      !! turning into a geostrophic balance whose amplitude is `a/f`), a
      !! flat bed, a linear ice draft, sigma, `pred_corr`, a pinned
      !! `n_inner`, and the load wired all the way through.
      logical, intent(in) :: legacy_split
         !! `&ocean_bt_nml bc_pgf_forcing = .not. legacy_split`.
      character(len=:), allocatable :: nml
      character(len=32) :: t_bot_s
      character(len=8) :: bcf

      ! Ghost columns keep whatever the IC seeds (the interior T is
      ! overwritten from GEOPOTENTIAL depth afterwards, and the API's
      ! tracer setter writes the interior only).  With WALL boundaries the
      ! ghosts reach nothing: the physical boundary faces are hard-zeroed
      ! by the wall closure, and the first interior cell's continuity uses
      ! that same zero flux.  Seeding the IC with the SURFACE value at
      ! both ends keeps the ghost columns unstratified and the mismatch
      ! trivially small anyway.
      write (t_bot_s, '(F12.4)') T_SURF
      bcf = ".true."
      if (legacy_split) bcf = ".false."
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
            "n_inner = 24, bc_pgf_forcing = "//trim(bcf)//" /"//new_line("a")// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 190.0, draft_slope = 1.0e-3, draft_x0 = -10000.0 /"// &
            new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function rest_nml

   subroutine run_rest_case(stratified, umax, vmax, ok, ubt_max, ubc_max, legacy_split)
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
      real(wp), intent(out), optional :: ubt_max
         !! Largest interior thickness-weighted DEPTH-MEAN zonal velocity.
      real(wp), intent(out), optional :: ubc_max
         !! Largest interior |u(k) - depth mean| (the baroclinic part).
      logical, intent(in), optional :: legacy_split
         !! Run under the legacy split (`bc_pgf_forcing = .false.`), which
         !! discards the depth-mean layer PGF.  Default: the MOM6 split.

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen, nxp, nyp, nzp, ngc
      real(wp), pointer :: h3(:, :, :), b2(:, :), u3(:, :, :), v3(:, :, :)
      real(wp), allocatable :: tval(:, :, :)
      character(len=:), allocatable :: nml
      integer :: i, j, k, ng
      real(wp) :: e_low, z_ctr, num, den, hf, ubar, ubt_s, ubc_s
      logical :: legacy

      ok = .false.
      umax = 0.0_wp
      vmax = 0.0_wp
      handle = c_null_ptr
      legacy = .false.
      if (present(legacy_split)) legacy = legacy_split
      nml = rest_nml(legacy)
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
      if (present(ubt_max) .or. present(ubc_max)) then
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h3, [int(nx), int(ny), int(nz)])
         ubt_s = 0.0_wp
         ubc_s = 0.0_wp
         do j = ng + 1, ng + int(nyp)
            do i = ng + 2, ng + int(nxp)
               num = 0.0_wp
               den = 0.0_wp
               do k = 1, int(nzp)
                  hf = 0.5_wp*(h3(i - 1, j, k) + h3(i, j, k))
                  num = num + hf*u3(i, j, k)
                  den = den + hf
               end do
               ubar = num/den
               ubt_s = max(ubt_s, abs(ubar))
               ubc_s = max(ubc_s, maxval(abs(u3(i, j, 1:int(nzp)) - ubar)))
            end do
         end do
         if (present(ubt_max)) ubt_max = ubt_s
         if (present(ubc_max)) ubc_max = ubc_s
      end if
      ok = all(ieee_is_finite(u3)) .and. all(ieee_is_finite(v3))

      status = rdb_ocean_destroy(handle)
      ok = ok .and. (status == OCEAN_STATUS_OK)
   end subroutine run_rest_case

   subroutine test_rest_stratified(error)
      !! THE GATE.  Sloping draft, flat bed, flat isopycnals, at rest, no
      !! forcing: the spurious velocity after `N_STEPS` must stay under
      !! the sigma-PGF truncation bound `REST_SAFETY * A_PEAK * T_TOTAL`.
      !!
      !! Run under the LEGACY split (`&ocean_bt_nml bc_pgf_forcing =
      !! .false.`), because what the bound describes is the part of the
      !! truncation that survives a depth-mean REPLACEMENT.  Under the
      !! default MOM6 split the depth-uniform load shortfall reaches the
      !! barotropic mode and dominates by three decades — that is
      !! `cavity_sloping_lid_load_shortfall_drives_bt`, below.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: umax, vmax, bound
      logical :: ok

      call run_rest_case(.true., umax, vmax, ok, legacy_split=.true.)
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

   subroutine test_load_shortfall_bt(error)
      !! The same resting stratified cavity under the DEFAULT split
      !! (`&ocean_bt_nml bc_pgf_forcing`, MOM6 `BT_force`).  The load
      !! `p_ice_ref = rho_0*g*z_draft` is the displaced weight at the
      !! REFERENCE density; the stratified column it floats on is heavier
      !! by `rho_0*N^2*z_draft^2/(2g)`, so the raw face force is short by
      !! the depth-UNIFORM
      !!
      !!     a_0 = N^2 * z_draft * slope          (module header, control (b))
      !!
      !! The legacy split discarded it with the depth mean; the MOM6 split
      !! hands it to the barotropic mode, which is exactly what it is: a
      !! bottom-pressure gradient the initial state does not balance.
      !! The closed basin answers with a gravity-wave adjustment to a
      !! surface tilt `a_0/g`: starting from `eta = 0` the tilt deficit
      !! `a_0*L/(2g)` at the ends rings as the gravest seiche, whose
      !! velocity amplitude is `(c/H)*a_0*L/(2g) = a_0*L/(2c)`.  That is
      !! the bound (with the deepest draft and the shallowest column, so
      !! it is the largest the geometry allows) and, since nothing else
      !! forces the depth mean, a lower bound at a tenth of it keeps the
      !! case from passing vacuously.  Measured (gfortran 15.1): 4.49e-4
      !! m/s against the 9.9e-4 m/s bound, after ~8 seiche periods of
      !! BEBT damping; under the legacy split the same run holds
      !! 7.1e-8 m/s — the truncation floor of the gate above.
      !!
      !! A TRIMMED initial surface (`&ocean_cavity_dyn_nml
      !! trim_ic_for_p_surf`, MOM6 `trim_for_ice`) starts balanced and does
      !! not show this — `trim_ic_balances_the_depth_mean_pfu` below.  This
      !! case seeds its temperature through the API instead of the zinit
      !! overlay the trim needs, so it keeps measuring the untrimmed
      !! adjustment on purpose.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: umax, vmax, ubt_max, ubc_max, draft_max, h_min, a_0, bound
      logical :: ok

      call run_rest_case(.true., umax, vmax, ok, ubt_max, ubc_max)
      call check(error, ok, "the resting stratified cavity must run and stay finite")
      if (allocated(error)) return

      draft_max = DRAFT_AT_X0 + DRAFT_SLOPE*(-DRAFT_X0 + real(NX_PHYS, wp)*DXY)
      h_min = BED - draft_max
      a_0 = N2_TARGET*draft_max*DRAFT_SLOPE
      bound = a_0*(real(NX_PHYS, wp)*DXY)/(2.0_wp*sqrt(GRAVITY*h_min))
      call check(error, ubt_max <= bound, &
                 "the barotropic response to the load shortfall must stay under the "// &
                 "gravest-seiche amplitude a_0*L/(2c)")
      if (allocated(error)) return
      call check(error, ubt_max >= 0.1_wp*bound, &
                 "the depth-mean load shortfall must REACH the barotropic mode under "// &
                 "the MOM6 split (a near-zero depth-mean response means it was discarded)")
      if (allocated(error)) return
      call check(error, vmax <= 1.0e-6_wp*bound, &
                 "the draft varies in x only: no meridional response")
   end subroutine test_load_shortfall_bt

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

   ! ==================================================================
   ! Trimmed initial condition (MOM6 TRIM_IC_FOR_P_SURF)
   ! ==================================================================

   pure function int_rho_linear(s) result(w)
      !! `int_{-s}^{0} rho(z) dz` for the affine trimmed-case profile —
      !! the displaced water's mass per unit area above depth `s`.
      real(wp), intent(in) :: s
      real(wp) :: w
      w = RHO_S_T*s - 0.5_wp*DRHO_DZ_T*s*s
   end function int_rho_linear

   subroutine test_trim_root(error)
      !! The helper's closed-form root IS the trim condition
      !! `g*int_{-s}^{0} rho dz = rho_ref*g*z_draft`, column by column, with
      !! `eta = z_draft - s`; an open-ocean column (`z_draft = 0`) and a
      !! GROUNDED one (`water < h_min_cavity`) are left exactly at 0.  The
      !! magnitude is checked against the leading-order estimate
      !! `eta ~ I(z_draft)/rho_ref`, `I = int (rho - rho_ref)`, so a sign
      !! or factor-of-two slip in the root cannot pass.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NC = 5
      real(wp) :: zd(NC, 1), water(NC, 1), eta(NC, 1)
      real(wp) :: s, resid, lead
      logical :: ok
      integer :: i

      zd(:, 1) = [0.0_wp, 11.0_wp, 342.9_wp, 494.0_wp, 700.0_wp]
      water(:, 1) = BED_T - zd(:, 1)
      call cavity_trim_eta_linear_impl(eta, ok, zd, water, 40.0_wp, RHO_REF_T, &
                                       RHO_S_T, DRHO_DZ_T, NC, 1)
      call check(error, ok, "an ISOMIP+ COLD column must admit a trim depth")
      if (allocated(error)) return
      call check(error, eta(1, 1) == 0.0_wp, "open ocean (z_draft = 0) is not trimmed")
      if (allocated(error)) return
      call check(error, eta(5, 1) == 0.0_wp, &
                 "a GROUNDED column (20 m < h_min_cavity) is not trimmed")
      if (allocated(error)) return
      do i = 2, 4
         s = zd(i, 1) - eta(i, 1)
         resid = int_rho_linear(s) - RHO_REF_T*zd(i, 1)
         call check(error, abs(resid) <= 1.0e-12_wp*RHO_REF_T*zd(i, 1), &
                    "the trimmed column's displaced mass must equal the load's "// &
                    "rho_ref*z_draft (the MOM6 trim_for_ice condition)")
         if (allocated(error)) return
         lead = (int_rho_linear(zd(i, 1)) - RHO_REF_T*zd(i, 1))/RHO_REF_T
         call check(error, abs(eta(i, 1) - lead) <= 0.01_wp*abs(lead), &
                    "eta must match the leading-order I(z_draft)/rho_ref to 1 %")
         if (allocated(error)) return
         call check(error, eta(i, 1) < 0.0_wp, &
                    "water lighter than rho_ref above the ice base needs a DEEPER top")
         if (allocated(error)) return
      end do
      ! The extremum of I sits where rho(-z_draft) = rho_ref, z = 342.9 m:
      ! -4.80e-2 m, the number the case's configure line prints.
      call check(error, abs(eta(3, 1) + 4.80e-2_wp) <= 1.0e-4_wp, &
                 "the deepest trim on the ISOMIP+ COLD profile is -4.80 cm")
   end subroutine test_trim_root

   pure subroutine seed_trim_column(ms, b, water, trimmed)
      !! The `cavity_sloping_lid_rest` resting state across a strip of
      !! columns (draft stepping by `DRAFT_STEP_T`), `NZ_T` sigma layers,
      !! the layer density sampled at each layer centre's GEOPOTENTIAL
      !! height (flat isopycnals), the Boussinesq-isostatic load — and,
      !! when `trimmed`, the column top moved to `-z_draft + eta_trim`.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(out) :: b(:, :)
      real(wp), intent(out) :: water(:, :)
         !! The seeded water column, trim included.
      logical, intent(in) :: trimmed
      integer :: i, j, k, nx, ny
      real(wp), allocatable :: zd(:, :), eta(:, :)
      real(wp) :: e_low
      logical :: ok

      nx = size(ms%p_top, 1)
      ny = size(ms%p_top, 2)
      allocate (zd(nx, ny), eta(nx, ny))
      do j = 1, ny
         do i = 1, nx
            b(i, j) = BED_T
            zd(i, j) = DRAFT0_T + DRAFT_STEP_T*real(i - 1, wp)
            water(i, j) = BED_T - zd(i, j)
         end do
      end do
      eta = 0.0_wp
      if (trimmed) call cavity_trim_eta_linear_impl(eta, ok, zd, water, 40.0_wp, &
                                                    RHO_REF_T, RHO_S_T, DRHO_DZ_T, nx, ny)
      do j = 1, ny
         do i = 1, nx
            water(i, j) = water(i, j) + eta(i, j)
            e_low = -BED_T
            do k = 1, NZ_T
               ms%h_layer(i, j, k) = water(i, j)/real(NZ_T, wp)
               ms%rho_layer(i, j, k) = RHO_S_T + DRHO_DZ_T*(e_low + 0.5_wp*ms%h_layer(i, j, k))
               e_low = e_low + ms%h_layer(i, j, k)
            end do
            ms%p_top(i, j) = (RHO_REF_T*GRAVITY)*zd(i, j)
         end do
      end do
   end subroutine seed_trim_column

   subroutine test_trim_balances_pfu(error)
      !! THE t = 0 MEASUREMENT behind `cavity_sloping_lid_rest`, at the PGF,
      !! on the case's own geometry, stratification and `nz = 15`.
      !!
      !! UNTRIMMED, every level below the ice sits `-g*I(z_draft)` off the
      !! open-ocean hydrostatic pressure (`I = int_{-z_draft}^{0}
      !! (rho - rho_ref) dz`), so the thickness-weighted depth mean of the
      !! FV_MOM6 face force carries the bottom-pressure gradient
      !! `g*(I(z_R) - I(z_L))/(rho_ref*dx)` — up to 1.8e-5 m/s^2 here, what
      !! the MOM6 split (`bc_pgf_forcing`) hands the barotropic mode.
      !!
      !! TRIMMED (`cavity_trim_eta_linear_impl`), every column's INTERFACE
      !! pressures equal the open ocean's at the same z (the layer-midpoint
      !! stack integrates a linear density exactly), and what is left of
      !! the depth mean is the discretisation's own truncation, derived:
      !!
      !!   * the Pass-3 in-layer integral `pa(top)*h + rho'*g*h^2/2` treats
      !!     each layer's density as uniform, short by `g*rho_z*h^3/12` per
      !!     column; sigma layers `h = W/nz` differ across the face, so the
      !!     depth mean carries `-g*rho_z*(W_L^3 - W_R^3)/(12*nz^2*rho_0*dx*Hbar)`
      !!     (~ N^2*W*dW/(4*nz^2*dx), 1-3e-8 m/s^2 here).  Under sigma it is
      !!     DEPTH-UNIFORM, which is why the legacy split never saw it;
      !!   * plus the top interface's trapezoid error `G = -(De^3/12)*rho_0*N^2`
      !!     (flat bed: the depth-summed `G(k) - G(k+1)` telescopes to it),
      !!     `N^2*D^3/(12*dx*Hbar)` ~ 4e-9 m/s^2 — bounded, not predicted.
      !!
      !! Both are checked: the untrimmed mean to 1 % of its peak against
      !! (shortfall + quadrature), the trimmed one against the quadrature
      !! term to within the trapezoid bound.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms_raw, ms_trim
      type(ocean_pressure_force_t) :: pgf_raw, pgf_trim
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), w_raw(:, :), w_trim(:, :)
      real(wp) :: mean_raw, mean_trim, hsum_r, hsum_t, hb, zl, zr
      real(wp) :: pred_load, pred_q_raw, pred_q_trim, g_top, n2
      real(wp) :: err_raw, err_trim, peak_load, peak_trim
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(30, 4, NG_U, DX_T, DX_T)
         nx = grid%nx_total
         ny = grid%ny_total
         ms_raw%nz_ml = NZ_T
         call ms_raw%init(grid)
         ms_trim%nz_ml = NZ_T
         call ms_trim%init(grid)
         allocate (b(nx, ny), w_raw(nx, ny), w_trim(nx, ny))
         call seed_trim_column(ms_raw, b, w_raw, .false.)
         call seed_trim_column(ms_trim, b, w_trim, .true.)
         call make_pgf_t(grid, pgf_raw)
         call run_pgf_u(grid, ms_raw, pgf_raw, b)
         call make_pgf_t(grid, pgf_trim)
         call run_pgf_u(grid, ms_trim, pgf_trim, b)

         n2 = -GRAVITY*DRHO_DZ_T/RHO_REF_T
         err_raw = 0.0_wp
         err_trim = 0.0_wp
         peak_load = 0.0_wp
         peak_trim = 0.0_wp
         do j = 3, ny - 1
            do i = 3, nx - 1
               mean_raw = 0.0_wp; mean_trim = 0.0_wp; hsum_r = 0.0_wp; hsum_t = 0.0_wp
               do k = 1, NZ_T
                  hb = 0.5_wp*(ms_raw%h_layer(i - 1, j, k) + ms_raw%h_layer(i, j, k))
                  mean_raw = mean_raw + hb*pgf_raw%dpdx_face%data(i, j, k)
                  hsum_r = hsum_r + hb
                  hb = 0.5_wp*(ms_trim%h_layer(i - 1, j, k) + ms_trim%h_layer(i, j, k))
                  mean_trim = mean_trim + hb*pgf_trim%dpdx_face%data(i, j, k)
                  hsum_t = hsum_t + hb
               end do
               mean_raw = mean_raw/hsum_r
               mean_trim = mean_trim/hsum_t
               ! u-face i sits between columns i-1 (L) and i (R).
               zl = DRAFT0_T + DRAFT_STEP_T*real(i - 2, wp)
               zr = DRAFT0_T + DRAFT_STEP_T*real(i - 1, wp)
               pred_load = GRAVITY*((int_rho_linear(zr) - RHO_REF_T*zr) - &
                                    (int_rho_linear(zl) - RHO_REF_T*zl))/(RHO_REF_T*DX_T)
               pred_q_raw = -GRAVITY*DRHO_DZ_T*(w_raw(i - 1, j)**3 - w_raw(i, j)**3)/ &
                            (12.0_wp*real(NZ_T*NZ_T, wp)*RHO_REF_T*DX_T*hsum_r)
               pred_q_trim = -GRAVITY*DRHO_DZ_T*(w_trim(i - 1, j)**3 - w_trim(i, j)**3)/ &
                             (12.0_wp*real(NZ_T*NZ_T, wp)*RHO_REF_T*DX_T*hsum_t)
               g_top = n2*abs(DRAFT_STEP_T)**3/(12.0_wp*DX_T*hsum_t)
               err_raw = max(err_raw, abs(mean_raw - (pred_load + pred_q_raw)))
               peak_load = max(peak_load, abs(pred_load))
               err_trim = max(err_trim, abs(mean_trim - pred_q_trim)/(1.5_wp*g_top))
               peak_trim = max(peak_trim, abs(mean_trim))
            end do
         end do
         call check(error, peak_load > 1.0e-5_wp, &
                    "the untrimmed load shortfall must be the 1e-5 m/s^2 force it is "// &
                    "on this geometry (else the comparison below is vacuous)")
         if (allocated(error)) exit checks
         call check(error, err_raw <= 0.01_wp*peak_load, &
                    "untrimmed, the depth-mean face force must be the load "// &
                    "shortfall's bottom-pressure gradient g*dI/dx/rho_ref (plus the "// &
                    "in-layer quadrature term) to 1 % of its peak")
         if (allocated(error)) exit checks
         call check(error, err_trim <= 1.0_wp, &
                    "trimmed, the depth-mean face force must be the in-layer "// &
                    "quadrature truncation to within the top trapezoid error")
         if (allocated(error)) exit checks
         call check(error, peak_trim <= 1.0e-2_wp*peak_load, &
                    "and the trim must remove the load shortfall: two decades or more")
      end block checks
      call pgf_raw%destroy(); call pgf_trim%destroy()
      call ms_raw%destroy(); call ms_trim%destroy()
   end subroutine test_trim_balances_pfu

   subroutine make_pgf_t(grid, pgf)
      !! FV_MOM6 with the load in the top BC, at the trimmed case's rho_ref.
      type(hgrid_t), intent(in) :: grid
      type(ocean_pressure_force_t), intent(out) :: pgf
      call pgf%init(grid, nz_ml=NZ_T)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO_REF_T
      pgf%rho_ref = RHO_REF_T
      pgf%p_top_in_bc = .true.
   end subroutine make_pgf_t

end module test_ocean_cavity_load
