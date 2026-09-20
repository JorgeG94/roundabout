!! Analytical + contract tests for the ice-shelf basal-melt COUPLING
!! (`&ocean_cavity_melt_nml`): the far-field sampler, the masked 2-D
!! driver, the two owned surface-flux components, the budget, the status
!! accounting and the refusal matrix.
module test_ocean_cavity_flux
   !! The kernel itself is tested to a 17-digit oracle in
   !! `test_ocean_cavity_melt`; NOTHING here re-derives three-equation
   !! physics.  What this suite pins is everything BETWEEN the ocean
   !! state and that kernel, which is where a coupling PR goes wrong:
   !!
   !!   * the far field is sampled over METRES, thickness-weighted, with
   !!     a partial last layer and vanished layers skipped;
   !!   * the 2-D driver reproduces scalar `cavity_melt_point` calls on
   !!     the sampled state, and leaves uncovered columns EXACTLY zero;
   !!   * the delivered signs COOL and FRESHEN (a sign error here is a
   !!     plausible, wrong, publishable answer);
   !!   * the virtual salt flux is the dilution identity, so the heat and
   !!     salt budgets close to round-off through the PRODUCTION
   !!     assembler + apply path;
   !!   * deeper ice ⇒ more melt, through `ms%p_top` (the ice pump);
   !!   * a non-finite covered column is FATAL and a refused one is
   !!     counted with zero melt;
   !!   * every v1 gap is a configure REFUSAL, not silence.
   !!
   !! `mem:separate` discipline.  Every kernel here is given
   !! device-present arrays: state objects through their own
   !! `enter_data`, bare test scratch through its own
   !! `!$acc enter data ... exit data`, and every host read-back is an
   !! `!$acc update self` of COMPONENT arrays, never of an aggregate
   !! derived type (the `test_ocean_ice_evp` segfault class).  All of it
   !! is inert on the host build, which is exactly why it is written
   !! unconditionally: a green gfortran run proves nothing about device
   !! data motion.
   !!
   !! Tolerances are STATED and never bit-equality across two different
   !! code paths — the driver and the scalar kernel compute the same
   !! algebra through different expression trees, and an FMA-contracting
   !! build is free to round them differently.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t, eos_apply_tfreeze_set, eos_freezing_point, &
                      TFREEZE_SET_ISOMIP, TFREEZE_SET_SEAICE
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_assemble, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_apply_sw_penetration, &
                                     ocean_surface_restore_apply_tracers
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_apply_cover
   use rdb_ocean_top_drag, only: top_drag_fill_face_cover_impl
   use rdb_ocean_cavity_melt, only: cavity_melt_point, cavity_ustar, &
                                    ocean_cavity_exchange_t, ocean_cavity_ice_t, &
                                    ocean_cavity_const_t, CAVITY_MELT_OK, &
                                    CAVITY_MELT_BAD_INPUT, &
                                    CAVITY_GAMMA_T_ISOMIP, CAVITY_GAMMA_RATIO_ISOMIP
   use rdb_ocean_cavity_flux, only: ocean_cavity_flux_t, ocean_cavity_flux_step, &
                                    cavity_far_field_impl, cavity_flux_fill_impl, &
                                    cavity_status_counts_impl, cavity_melt_status_is_fatal
   use rdb_ocean_setup, only: cavity_resolve_gamma_s
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_cavity_flux_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 4
   integer, parameter :: NY_PHYS = 3
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: DY = 2000.0_wp

   ! A deliberately WARM cavity: +0.5 degC at S = 34.6 under 500 m of ice
   ! melts at a few tens of m/yr, i.e. comfortably above round-off in
   ! every assertion below without being a physically silly state.
   real(wp), parameter :: T_WARM = 0.5_wp
   real(wp), parameter :: S_REF = 34.6_wp
   real(wp), parameter :: P_DRAFT = 4.5e6_wp
      !! ~500 m of ice, `rho_i*g*h` in Pa.
   real(wp), parameter :: H0 = 50.0_wp
   real(wp), parameter :: FAR_DEPTH = 10.0_wp
   real(wp), parameter :: F_COR = -1.4e-4_wp

   real(wp), parameter :: TOL_PATH = 1.0e-12_wp
      !! Relative tolerance for "two code paths, one algebra" (FMA).
   real(wp), parameter :: TOL_EXACT = 1.0e-14_wp
      !! Relative tolerance for a thickness-weighted mean of exact data.

contains

   subroutine collect_ocean_cavity_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_far_field_uniform_is_exact", test_ff_uniform), &
                  new_unittest("cavity_far_field_partial_layer_weighted", test_ff_partial), &
                  new_unittest("cavity_far_field_depth_exceeds_column", test_ff_deep), &
                  new_unittest("cavity_far_field_skips_vanished_layers", test_ff_vanished), &
                  new_unittest("cavity_far_field_velocity_is_centred", test_ff_velocity), &
                  new_unittest("cavity_far_field_inactive_off_cover", test_ff_mask), &
                  new_unittest("cavity_driver_matches_scalar_kernel", test_driver_vs_kernel), &
                  new_unittest("cavity_flux_signs_cool_and_freshen", test_flux_signs), &
                  new_unittest("cavity_virtual_salt_is_dilution_identity", test_virtual_salt), &
                  new_unittest("cavity_budget_closes_to_roundoff", test_budget_closes), &
                  new_unittest("cavity_ice_pump_deeper_melts_more", test_ice_pump), &
                  new_unittest("cavity_status_nonfinite_is_fatal", test_status_nonfinite), &
                  new_unittest("cavity_status_fresh_column_counted", test_status_fresh), &
                  new_unittest("cavity_gamma_s_sentinel_resolves", test_gamma_s_sentinel), &
                  new_unittest("cavity_melt_off_is_placeholder", test_off_placeholder), &
                  new_unittest("cavity_cover_absent_is_unmasked", test_cover_absent), &
                  new_unittest("cavity_cover_stress_zero_wind_under_ice", test_cover_stress), &
                  new_unittest("cavity_cover_face_rule_matches_top_drag", &
                               test_cover_face_rule), &
                  new_unittest("cavity_cover_assembler_masks_atmosphere", test_cover_assemble), &
                  new_unittest("cavity_cover_blocks_sw_and_restore", test_cover_sw_restore), &
                  new_unittest("cavity_cover_budget_closes_both_sources", test_cover_budget), &
                  new_unittest("cavity_melt_validate_refusals", test_validate_refusals) &
                  ]
   end subroutine collect_ocean_cavity_flux_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   function make_eos_isomip() result(eos)
      !! The ISOMIP+ liquidus, through the PRODUCTION handle path — never
      !! a local copy of the coefficients.
      type(eos_t) :: eos
      call eos_apply_tfreeze_set(eos, TFREEZE_SET_ISOMIP)
   end function make_eos_isomip

   pure function default_par() result(par)
      !! The shipped default exchange bundle: `const_gamma` at the
      !! ISOMIP+ coefficients.
      type(ocean_cavity_exchange_t) :: par
      par%gamma_t_coeff = CAVITY_GAMMA_T_ISOMIP
      par%gamma_s_coeff = CAVITY_GAMMA_T_ISOMIP/CAVITY_GAMMA_RATIO_ISOMIP
      par%f_cor = F_COR
   end function default_par

   subroutine sample_one_column(nz, far_depth, h, t, s, u_w, u_e, v_s, v_n, &
                                active, t_far, s_far, u_far, v_far)
      !! Run `cavity_far_field_impl` on a ONE-cell plane with a given
      !! column, with the full `mem:separate` dance (the impl is a device
      !! kernel; its inputs and outputs must be mapped by this caller).
      integer, intent(in) :: nz
      real(wp), intent(in) :: far_depth
      real(wp), intent(in) :: h(nz), t(nz), s(nz)
      real(wp), intent(in) :: u_w(nz), u_e(nz), v_s(nz), v_n(nz)
      real(wp), intent(out) :: active, t_far, s_far, u_far, v_far
      real(wp) :: cover(1, 1), wet(1, 1)
      real(wp) :: h3(1, 1, nz), ht(1, 1, nz), hs(1, 1, nz)
      real(wp) :: uf(2, 1, nz), vf(1, 2, nz)
      real(wp) :: a2(1, 1), t2(1, 1), s2(1, 1), u2(1, 1), v2(1, 1)
      integer :: k

      cover = 1.0_wp
      wet = 1.0_wp
      do k = 1, nz
         h3(1, 1, k) = h(k)
         ht(1, 1, k) = h(k)*t(k)
         hs(1, 1, k) = h(k)*s(k)
         uf(1, 1, k) = u_w(k)
         uf(2, 1, k) = u_e(k)
         vf(1, 1, k) = v_s(k)
         vf(1, 2, k) = v_n(k)
      end do

      !$acc enter data copyin(cover, wet, h3, ht, hs, uf, vf) &
      !$acc&           create(a2, t2, s2, u2, v2)
      call cavity_far_field_impl(1, 1, nz, far_depth, cover, wet, h3, ht, hs, uf, vf, &
                                 a2, t2, s2, u2, v2)
      !$acc update self(a2, t2, s2, u2, v2)
      !$acc exit data delete(cover, wet, h3, ht, hs, uf, vf, a2, t2, s2, u2, v2)

      active = a2(1, 1)
      t_far = t2(1, 1)
      s_far = s2(1, 1)
      u_far = u2(1, 1)
      v_far = v2(1, 1)
   end subroutine sample_one_column

   subroutine build_cavity_plane(grid, metrics, ms, sf, cav, eos, p_top_val, t_val, &
                                 u_val, cover_x_max)
      !! Build a small covered plane: metrics with a cavity, a uniform
      !! stratification-free column, the surface-flux component set, and
      !! a configured melt slot — then map all of it to the device.
      !!
      !! `cover_x_max` is the last PHYSICAL i index that carries ice, so
      !! a single call gives a plane with BOTH covered and uncovered
      !! columns (the calving-front geometry the driver must handle).
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_cavity_flux_t), intent(inout) :: cav
      type(eos_t), intent(out) :: eos
      real(wp), intent(in) :: p_top_val, t_val, u_val
      integer, intent(in) :: cover_x_max
      integer :: i, j, nx, ny

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total

      metrics%use_cavity = .true.
      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = t_val*H0
      ms%tracers(ms%idx_salinity)%hTr = S_REF*H0
      ms%u_face_x_layer = u_val
      ms%v_face_y_layer = 0.0_wp
      ms%wet_mask = 1.0_wp
      ms%p_top = p_top_val

      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      cav%enable = .true.
      call cav%init(grid)
      cav%far_field_depth = FAR_DEPTH
      cav%s_ice = 0.0_wp
      cav%par = default_par()
      cav%f_cor = F_COR
      eos = make_eos_isomip()

      ! Binary cover, exactly as `cavity_fill_cover_frac` builds it, but
      ! stopped at `cover_x_max` so the plane has a calving front.
      metrics%cover_frac = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            if (i - NGHOST >= 1 .and. i - NGHOST <= cover_x_max) then
               metrics%cover_frac(i, j) = 1.0_wp
            end if
         end do
      end do

      ! `metrics%enter_data` already ran inside `make_cartesian_metrics`,
      ! so the cover mask written after it owes an explicit push — the
      ! canonical `mem:separate` trap (2).
      !$acc update device(metrics%cover_frac)
      !$acc enter data copyin(ms, sf, cav)
      call ms%enter_data()
      call sf%enter_data()
      call cav%enter_data()
   end subroutine build_cavity_plane

   subroutine teardown_cavity_plane(metrics, ms, sf, cav)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_cavity_flux_t), intent(inout) :: cav
      call cav%exit_data()
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf, cav)
      call cav%destroy()
      call sf%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown_cavity_plane

   pure function rel_diff(a, b) result(r)
      !! Relative difference with an absolute fallback, so a legitimately
      !! zero expectation does not divide by zero.
      real(wp), intent(in) :: a, b
      real(wp) :: r
      r = abs(a - b)/max(abs(b), 1.0_wp)
   end function rel_diff

   ! ------------------------------------------------------------------
   ! Far-field sampler
   ! ------------------------------------------------------------------

   subroutine test_ff_uniform(error)
      !! A uniform column returns its own values EXACTLY, whatever the
      !! sampling depth — the thickness weights sum to one.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ), t(NZ), s(NZ), uw(NZ), ue(NZ), vs(NZ), vn(NZ)
      real(wp) :: act, tf, sf_, uf, vf

      h = H0
      t = T_WARM
      s = S_REF
      uw = 0.03_wp
      ue = 0.03_wp
      vs = -0.01_wp
      vn = -0.01_wp

      call sample_one_column(NZ, FAR_DEPTH, h, t, s, uw, ue, vs, vn, act, tf, sf_, uf, vf)
      call check(error, act > 0.5_wp, "a covered wet column with mass must be active")
      if (allocated(error)) return
      call check(error, rel_diff(tf, T_WARM) <= TOL_EXACT, "uniform T sampled exactly")
      if (allocated(error)) return
      call check(error, rel_diff(sf_, S_REF) <= TOL_EXACT, "uniform S sampled exactly")
      if (allocated(error)) return
      call check(error, rel_diff(uf, 0.03_wp) <= TOL_EXACT, "uniform u sampled exactly")
      if (allocated(error)) return
      call check(error, rel_diff(vf, -0.01_wp) <= TOL_EXACT, "uniform v sampled exactly")
   end subroutine test_ff_uniform

   subroutine test_ff_partial(error)
      !! A two-layer step profile with `far_field_depth` CUTTING the
      !! second layer must give the exact thickness-weighted mean — this
      !! is the "metres, not layers" contract, and it is the one thing a
      !! "just take layer nz" implementation gets wrong.
      !!
      !! Surface layer (k = nz) is 4 m at T = 1.0; the far field is 10 m,
      !! so 6 m of the next layer (T = -1.0) is included:
      !!   T_far = (4*1.0 + 6*(-1.0))/10 = -0.2 exactly.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ), t(NZ), s(NZ), uw(NZ), ue(NZ), vs(NZ), vn(NZ)
      real(wp) :: act, tf, sf_, uf, vf, expect_t, expect_s

      h = 100.0_wp
      h(NZ) = 4.0_wp
      h(NZ - 1) = 20.0_wp
      t = -1.0_wp
      t(NZ) = 1.0_wp
      s = 34.0_wp
      s(NZ) = 35.0_wp
      uw = 0.0_wp
      ue = 0.0_wp
      vs = 0.0_wp
      vn = 0.0_wp

      expect_t = (4.0_wp*1.0_wp + 6.0_wp*(-1.0_wp))/10.0_wp
      expect_s = (4.0_wp*35.0_wp + 6.0_wp*34.0_wp)/10.0_wp

      call sample_one_column(NZ, 10.0_wp, h, t, s, uw, ue, vs, vn, act, tf, sf_, uf, vf)
      call check(error, act > 0.5_wp, "active")
      if (allocated(error)) return
      call check(error, rel_diff(tf, expect_t) <= TOL_EXACT, &
                 "partial-layer thickness-weighted T must be exact")
      if (allocated(error)) return
      call check(error, rel_diff(sf_, expect_s) <= TOL_EXACT, &
                 "partial-layer thickness-weighted S must be exact")
   end subroutine test_ff_partial

   subroutine test_ff_deep(error)
      !! `far_field_depth` >= the column thickness uses the WHOLE column
      !! and is not an error: the column-mean is the correct limit.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ), t(NZ), s(NZ), uw(NZ), ue(NZ), vs(NZ), vn(NZ)
      real(wp) :: act, tf, sf_, uf, vf, expect_t
      integer :: k

      h = 5.0_wp
      t = [(real(k, wp), k=1, NZ)]
      s = S_REF
      uw = 0.0_wp
      ue = 0.0_wp
      vs = 0.0_wp
      vn = 0.0_wp
      expect_t = sum(t)/real(NZ, wp)

      call sample_one_column(NZ, 1.0e6_wp, h, t, s, uw, ue, vs, vn, act, tf, sf_, uf, vf)
      call check(error, act > 0.5_wp, "active")
      if (allocated(error)) return
      call check(error, rel_diff(tf, expect_t) <= TOL_EXACT, &
                 "far_field_depth beyond the column must average the whole column")
   end subroutine test_ff_deep

   subroutine test_ff_vanished(error)
      !! Vanished layers at the TOP of the stack are SKIPPED, not
      !! clamped: the sample starts at the first massive layer, and the
      !! answer is identical to a stack without them.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ), t(NZ), s(NZ), uw(NZ), ue(NZ), vs(NZ), vn(NZ)
      real(wp) :: act, tf, sf_, uf, vf

      h = 100.0_wp
      h(NZ) = 0.5_wp*H_VANISHED
      h(NZ - 1) = 0.5_wp*H_VANISHED
      t = T_WARM
      t(NZ) = 999.0_wp
      t(NZ - 1) = -999.0_wp
      s = S_REF
      uw = 0.0_wp
      ue = 0.0_wp
      vs = 0.0_wp
      vn = 0.0_wp

      call sample_one_column(NZ, FAR_DEPTH, h, t, s, uw, ue, vs, vn, act, tf, sf_, uf, vf)
      call check(error, act > 0.5_wp, "a column with two vanished top layers is still active")
      if (allocated(error)) return
      call check(error, rel_diff(tf, T_WARM) <= TOL_EXACT, &
                 "vanished layers must contribute NOTHING (not even their tracer value)")
   end subroutine test_ff_vanished

   subroutine test_ff_velocity(error)
      !! The C-grid faces are averaged to the CELL CENTRE per layer, then
      !! averaged vertically.  An asymmetric face pair proves the centring
      !! actually happens rather than one face being picked.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ), t(NZ), s(NZ), uw(NZ), ue(NZ), vs(NZ), vn(NZ)
      real(wp) :: act, tf, sf_, uf, vf

      h = H0
      t = T_WARM
      s = S_REF
      uw = 0.10_wp
      ue = 0.20_wp
      vs = -0.04_wp
      vn = 0.00_wp

      call sample_one_column(NZ, FAR_DEPTH, h, t, s, uw, ue, vs, vn, act, tf, sf_, uf, vf)
      call check(error, rel_diff(uf, 0.15_wp) <= TOL_EXACT, &
                 "u must be the 2-point face average, not a single face")
      if (allocated(error)) return
      call check(error, rel_diff(vf, -0.02_wp) <= TOL_EXACT, &
                 "v must be the 2-point face average")
   end subroutine test_ff_velocity

   subroutine test_ff_mask(error)
      !! Uncovered, dry, and all-vanished columns are INACTIVE with
      !! zeroed outputs — the three ways a column legitimately has no
      !! interface.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: cover(3, 1), wet(3, 1)
      real(wp) :: h3(3, 1, 1), ht(3, 1, 1), hs(3, 1, 1)
      real(wp) :: uf(4, 1, 1), vf(3, 2, 1)
      real(wp) :: a2(3, 1), t2(3, 1), s2(3, 1), u2(3, 1), v2(3, 1)

      cover = reshape([1.0_wp, 0.0_wp, 1.0_wp], [3, 1])
      wet = reshape([0.0_wp, 1.0_wp, 1.0_wp], [3, 1])
      h3 = reshape([0.5_wp*H_VANISHED, H0, H0], [3, 1, 1])
      ht = h3*T_WARM
      hs = h3*S_REF
      uf = 0.0_wp
      vf = 0.0_wp

      !$acc enter data copyin(cover, wet, h3, ht, hs, uf, vf) &
      !$acc&           create(a2, t2, s2, u2, v2)
      call cavity_far_field_impl(3, 1, 1, FAR_DEPTH, cover, wet, h3, ht, hs, uf, vf, &
                                 a2, t2, s2, u2, v2)
      !$acc update self(a2, t2, s2, u2, v2)
      !$acc exit data delete(cover, wet, h3, ht, hs, uf, vf, a2, t2, s2, u2, v2)

      call check(error, a2(1, 1) == 0.0_wp, "a DRY covered column must be inactive")
      if (allocated(error)) return
      call check(error, a2(2, 1) == 0.0_wp, "an UNCOVERED wet column must be inactive")
      if (allocated(error)) return
      call check(error, a2(3, 1) == 1.0_wp, "a covered wet massive column must be active")
      if (allocated(error)) return
      call check(error, t2(1, 1) == 0.0_wp .and. s2(2, 1) == 0.0_wp, &
                 "inactive columns must carry exactly zero, never a stale sample")
   end subroutine test_ff_mask

   ! ------------------------------------------------------------------
   ! Driver vs kernel
   ! ------------------------------------------------------------------

   subroutine test_driver_vs_kernel(error)
      !! On a plane with a calving front, every COVERED column's
      !! `(melt, T_b, S_b)` equals a direct scalar `cavity_melt_point`
      !! call on the sampled far field, and every UNCOVERED column is
      !! exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: const
      real(wp) :: us_ref, tb_ref, sb_ref, m_ref, q_ref
      integer :: ierr_ref, ierr_us, i, j, nx, ny, i_cov, i_open

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, 2)
      nx = grid%nx_total
      ny = grid%ny_total

      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)

      !$acc update self(cav%melt, cav%t_b, cav%s_b, cav%ustar, cav%t_far, cav%s_far, &
      !$acc&            cav%u_far, cav%v_far, cav%active, cav%status)
      !$acc update self(sf%heat_cavity, sf%salt_cavity)

      i_cov = NGHOST + 1
      i_open = NGHOST + 4
      j = NGHOST + 1

      ! The scalar reference: the SAME sampled far field, through the
      ! scalar entry point, with no 2-D machinery in the way.
      call cavity_ustar(cav%u_far(i_cov, j), cav%v_far(i_cov, j), cav%cdrag_top, &
                        cav%u_tide, cav%ustar_min, us_ref, ierr_us)
      call cavity_melt_point(cav%t_far(i_cov, j), cav%s_far(i_cov, j), &
                             ms%p_top(i_cov, j), us_ref, cav%s_ice, cav%par, ice, &
                             eos, const, tb_ref, sb_ref, m_ref, q_ref, ierr_ref)

      call check(error, ierr_us == CAVITY_MELT_OK .and. ierr_ref == CAVITY_MELT_OK, &
                 "the reference scalar solve must succeed")
      if (allocated(error)) go to 900
      call check(error, m_ref > 0.0_wp, "warm water under ice must MELT (m_mass > 0)")
      if (allocated(error)) go to 900
      call check(error, rel_diff(cav%ustar(i_cov, j), us_ref) <= TOL_PATH, "u* matches")
      if (allocated(error)) go to 900
      call check(error, rel_diff(cav%melt(i_cov, j), m_ref) <= TOL_PATH, "melt matches")
      if (allocated(error)) go to 900
      call check(error, rel_diff(cav%t_b(i_cov, j), tb_ref) <= TOL_PATH, "T_b matches")
      if (allocated(error)) go to 900
      call check(error, rel_diff(cav%s_b(i_cov, j), sb_ref) <= TOL_PATH, "S_b matches")
      if (allocated(error)) go to 900

      ! Uncovered columns: EXACTLY zero, and OK status.
      call check(error, cav%active(i_open, j) == 0.0_wp, "beyond the front: inactive")
      if (allocated(error)) go to 900
      call check(error, cav%melt(i_open, j) == 0.0_wp .and. &
                 cav%t_b(i_open, j) == 0.0_wp .and. cav%s_b(i_open, j) == 0.0_wp, &
                 "an uncovered column must be EXACTLY zero, not nearly zero")
      if (allocated(error)) go to 900
      call check(error, cav%status(i_open, j) == CAVITY_MELT_OK, &
                 "an uncovered column is not a solver failure")
      if (allocated(error)) go to 900
      call check(error, sf%heat_cavity(i_open, j) == 0.0_wp .and. &
                 sf%salt_cavity(i_open, j) == 0.0_wp, &
                 "an uncovered column's flux components stay untouched at zero")
      if (allocated(error)) go to 900
      call check(error, cav%n_nonfinite_step == 0 .and. cav%n_not_converged_step == 0 &
                 .and. cav%n_no_root_step == 0 .and. cav%n_other_step == 0, &
                 "a healthy plane reports no refusals")
      if (allocated(error)) go to 900

      ! Every covered column carries the same state, so the whole covered
      ! band must agree with the one reference (a per-column indexing slip
      ! would show up here and nowhere else).
      do j = NGHOST + 1, ny - NGHOST
         do i = NGHOST + 1, NGHOST + 2
            if (rel_diff(cav%melt(i, j), m_ref) > TOL_PATH) then
               call check(error, .false., "every covered column must match the reference")
               go to 900
            end if
         end do
      end do

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_driver_vs_kernel

   ! ------------------------------------------------------------------
   ! Signs
   ! ------------------------------------------------------------------

   subroutine test_flux_signs(error)
      !! Warm water under ice COOLS and FRESHENS the top of the column —
      !! asserted end to end, through the production assembler and apply
      !! path, on the actual tracer loads.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp) :: t_before, s_before, t_after, s_after
      integer :: i, j

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, NX_PHYS)
      i = NGHOST + 1
      j = NGHOST + 1
      t_before = T_WARM
      s_before = S_REF

      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
      call ocean_surface_flux_assemble(grid, sf, ms)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, 3600.0_wp)

      !$acc update self(sf%heat_cavity, sf%salt_cavity, sf%Q_heat, sf%Q_salt)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr, &
      !$acc&            ms%tracers(ms%idx_salinity)%hTr)

      t_after = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H0
      s_after = ms%tracers(ms%idx_salinity)%hTr(i, j, NZ)/H0

      call check(error, sf%heat_cavity(i, j) < 0.0_wp, &
                 "melting must take heat OUT of the ocean (heat_cavity < 0)")
      if (allocated(error)) go to 900
      call check(error, sf%salt_cavity(i, j) < 0.0_wp, &
                 "meltwater must FRESHEN (salt_cavity < 0)")
      if (allocated(error)) go to 900
      call check(error, sf%Q_heat(i, j) < 0.0_wp .and. sf%Q_salt(i, j) < 0.0_wp, &
                 "the assembler must carry both signs through to Q_heat/Q_salt")
      if (allocated(error)) go to 900
      call check(error, t_after < t_before, "the top layer must COOL")
      if (allocated(error)) go to 900
      call check(error, s_after < s_before, "the top layer must FRESHEN")
      if (allocated(error)) go to 900
      call check(error, sf%has_heat .and. sf%has_salt, &
                 "the filler must latch has_heat/has_salt host-side")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_flux_signs

   subroutine test_virtual_salt(error)
      !! The delivered salt component IS the dilution identity
      !! `-m_mass*(S_far - s_ice)`, to round-off.  This is the whole
      !! virtual-salt decision in one assertion: get it wrong and the
      !! meltwater buoyancy — the thing that drives a cavity — is wrong
      !! by a factor, with no other symptom.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp) :: expect
      integer :: i, j

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, NX_PHYS)
      i = NGHOST + 1
      j = NGHOST + 1
      cav%s_ice = 2.0_wp
      !$acc update self(cav%active)
      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
      !$acc update self(cav%melt, cav%s_far, cav%q_ocean)
      !$acc update self(sf%heat_cavity, sf%salt_cavity)

      expect = -cav%melt(i, j)*(cav%s_far(i, j) - cav%s_ice)
      call check(error, rel_diff(sf%salt_cavity(i, j), expect) <= TOL_PATH, &
                 "salt_cavity must be -m_mass*(S_far - s_ice)")
      if (allocated(error)) go to 900
      call check(error, rel_diff(sf%heat_cavity(i, j), -cav%q_ocean(i, j)) <= TOL_PATH, &
                 "heat_cavity must be -q_ocean")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_virtual_salt

   ! ------------------------------------------------------------------
   ! Budget
   ! ------------------------------------------------------------------

   subroutine test_budget_closes(error)
      !! Over several thermo steps the change in the domain heat/salt
      !! integral equals the accumulated surface-budget contributor, to
      !! round-off — i.e. the cavity source is FULLY accounted for by the
      !! path the console reports, with no separate accumulator needed.
      !!
      !! Both sides are reduced with the same stencil the console uses
      !! (physical interior only), so the residual is pure summation
      !! error.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp) :: heat0, salt0, heat1, salt1, bud_h, bud_s, scale
      integer :: nx, ny, step
      integer, parameter :: N_STEP = 5
      real(wp), parameter :: DT = 1800.0_wp

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, NX_PHYS - 1)
      nx = grid%nx_total
      ny = grid%ny_total
      scale = grid%dx*grid%dy*RHO_WATER

      heat0 = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                     NGHOST + 1:ny - NGHOST, :))*scale
      salt0 = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                  NGHOST + 1:ny - NGHOST, :))*scale

      do step = 1, N_STEP
         call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
         call ocean_surface_flux_assemble(grid, sf, ms)
         call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
      end do

      !$acc update self(ms%tracers(ms%idx_temperature)%hTr, &
      !$acc&            ms%tracers(ms%idx_salinity)%hTr, &
      !$acc&            ms%heat_budget_surface, ms%salt_budget_surface)

      heat1 = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                     NGHOST + 1:ny - NGHOST, :))*scale
      salt1 = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                  NGHOST + 1:ny - NGHOST, :))*scale
      bud_h = sum(ms%heat_budget_surface(NGHOST + 1:nx - NGHOST, &
                                         NGHOST + 1:ny - NGHOST, :))*scale
      bud_s = sum(ms%salt_budget_surface(NGHOST + 1:nx - NGHOST, &
                                         NGHOST + 1:ny - NGHOST, :))*scale

      call check(error, abs(bud_h) > 0.0_wp, "the cavity must have moved some heat")
      if (allocated(error)) go to 900
      call check(error, abs(bud_s) > 0.0_wp, "the cavity must have moved some salt")
      if (allocated(error)) go to 900
      ! Normalised the way the console's `relative_drift` normalises:
      ! against the INITIAL TOTAL, not against the source.  That is not a
      ! weakened gate, it is the only meaningful one in double
      ! precision — a five-step cavity source is ~1e-10 of the domain
      ! salt integral, so `salt1 - salt0` cancels fifteen digits before
      ! the budget term is even consulted, and a "1e-12 of the source"
      ! bound would be asserting something finite arithmetic cannot
      ! deliver.  The residual here is pure summation error and comes in
      ! two to three orders BELOW this bound.
      call check(error, abs((heat1 - heat0) - bud_h) <= 1.0e-13_wp*abs(heat0), &
                 "the heat budget must close to round-off with the cavity on")
      if (allocated(error)) go to 900
      call check(error, abs((salt1 - salt0) - bud_s) <= 1.0e-13_wp*abs(salt0), &
                 "the salt budget must close to round-off with the cavity on")
      if (allocated(error)) go to 900
      call check(error, bud_h < 0.0_wp .and. bud_s < 0.0_wp, &
                 "the accumulated cavity source must be a heat AND salt sink")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_budget_closes

   ! ------------------------------------------------------------------
   ! Ice pump
   ! ------------------------------------------------------------------

   subroutine test_ice_pump(error)
      !! THE ice pump, through `ms%p_top`: identical water, a deeper
      !! interface (higher pressure) depresses the freezing point, which
      !! RAISES the thermal driving and therefore the melt rate.  This is
      !! also the end-to-end proof that the melt liquidus reads the
      !! interface pressure at all — with `p_top` ignored the two runs
      !! would be identical.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m_shallow, m_deep

      call melt_at_pressure(1.0e6_wp, m_shallow, error)
      if (allocated(error)) return
      call melt_at_pressure(6.0e6_wp, m_deep, error)
      if (allocated(error)) return

      call check(error, m_shallow > 0.0_wp, "the shallow case must melt at all")
      if (allocated(error)) return
      call check(error, m_deep > m_shallow, &
                 "a deeper interface pressure must melt MORE (the ice pump)")
   end subroutine test_ice_pump

   subroutine melt_at_pressure(p_top_val, m_out, error)
      !! One covered plane at a given interface pressure; returns the
      !! melt rate at a representative covered column.
      real(wp), intent(in) :: p_top_val
      real(wp), intent(out) :: m_out
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, p_top_val, T_WARM, &
                              0.05_wp, NX_PHYS)
      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
      !$acc update self(cav%melt, cav%status)
      m_out = cav%melt(NGHOST + 1, NGHOST + 1)
      call check(error, cav%status(NGHOST + 1, NGHOST + 1) == CAVITY_MELT_OK, &
                 "the ice-pump probe column must solve cleanly")
      call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine melt_at_pressure

   ! ------------------------------------------------------------------
   ! Status handling
   ! ------------------------------------------------------------------

   subroutine test_status_nonfinite(error)
      !! A NaN in one covered column's far field is counted as
      !! NONFINITE and the FAIL-LOUD PREDICATE fires.  The predicate is
      !! asserted, never the `error stop` — the repo's standing pattern
      !! for testing a fail-loud rule without killing the runner.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: act(3, 1)
      integer :: st(3, 1)
      integer :: n_nf, n_nc, n_nr, n_ot

      act = 1.0_wp
      st(1, 1) = CAVITY_MELT_OK
      st(2, 1) = 1   ! CAVITY_MELT_NONFINITE_INPUT
      st(3, 1) = 5   ! CAVITY_MELT_NOT_CONVERGED

      !$acc enter data copyin(act, st)
      call cavity_status_counts_impl(3, 1, act, st, n_nf, n_nc, n_nr, n_ot)
      !$acc exit data delete(act, st)

      call check(error, n_nf == 1, "one non-finite column must be counted")
      if (allocated(error)) return
      call check(error, n_nc == 1, "one non-converged column must be counted")
      if (allocated(error)) return
      call check(error, n_nr == 0 .and. n_ot == 0, "nothing else is counted")
      if (allocated(error)) return
      call check(error, cavity_melt_status_is_fatal(n_nf), &
                 "a non-finite covered column must be FATAL")
      if (allocated(error)) return
      call check(error,.not. cavity_melt_status_is_fatal(0), &
                 "NOT_CONVERGED alone must NOT be fatal (zero melt is the safe state)")
      if (allocated(error)) return

      ! And an INACTIVE column's status is never counted, whatever it says.
      act(2, 1) = 0.0_wp
      !$acc enter data copyin(act, st)
      call cavity_status_counts_impl(3, 1, act, st, n_nf, n_nc, n_nr, n_ot)
      !$acc exit data delete(act, st)
      call check(error, n_nf == 0, "an inactive column is never a solver failure")
   end subroutine test_status_nonfinite

   subroutine test_status_fresh(error)
      !! A covered column whose far-field salinity is at or below the ice
      !! salinity is REFUSED by the kernel (`CAVITY_MELT_BAD_INPUT` — the
      !! root-bracketing argument rests on `S_w > S_i`), counted in the
      !! "other" bucket, and given EXACTLY zero melt and zero flux.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      integer :: i, j

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, NX_PHYS)
      ! Ice salinity above the water's: every covered column is refused.
      cav%s_ice = S_REF + 1.0_wp
      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
      !$acc update self(cav%melt, cav%status)
      !$acc update self(sf%heat_cavity, sf%salt_cavity)

      i = NGHOST + 1
      j = NGHOST + 1
      call check(error, cav%status(i, j) == CAVITY_MELT_BAD_INPUT, &
                 "S_w <= S_i must be refused, not guessed at")
      if (allocated(error)) go to 900
      call check(error, cav%melt(i, j) == 0.0_wp, &
                 "a refused column takes the documented zero-melt safe state")
      if (allocated(error)) go to 900
      call check(error, sf%heat_cavity(i, j) == 0.0_wp .and. &
                 sf%salt_cavity(i, j) == 0.0_wp, &
                 "a refused column delivers exactly zero flux")
      if (allocated(error)) go to 900
      call check(error, cav%n_other_step > 0 .and. cav%n_nonfinite_step == 0, &
                 "a refusal is COUNTED, and it is not the fatal bucket")
      if (allocated(error)) go to 900
      call check(error, cav%n_other_total == int(cav%n_other_step, kind(cav%n_other_total)), &
                 "the running total must absorb the step count")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_status_fresh

   ! ------------------------------------------------------------------
   ! Knobs
   ! ------------------------------------------------------------------

   subroutine test_gamma_s_sentinel(error)
      !! The negative `gamma_s` sentinel resolves to the ISOMIP+
      !! `gamma_t/35`; zero and positive values are taken LITERALLY (zero
      !! is then refused by `validate_config`, not silently defaulted).
      type(error_type), allocatable, intent(out) :: error
      call check(error, rel_diff(cavity_resolve_gamma_s(-1.0_wp, 2.2e-2_wp), &
                                 2.2e-2_wp/35.0_wp) <= TOL_EXACT, &
                 "negative gamma_s must resolve to gamma_t/35")
      if (allocated(error)) return
      call check(error, cavity_resolve_gamma_s(1.0e-3_wp, 2.2e-2_wp) == 1.0e-3_wp, &
                 "a positive gamma_s must be taken literally")
      if (allocated(error)) return
      call check(error, cavity_resolve_gamma_s(0.0_wp, 2.2e-2_wp) == 0.0_wp, &
                 "zero must be taken literally (validate_config then refuses it)")
   end subroutine test_gamma_s_sentinel

   subroutine test_off_placeholder(error)
      !! Knob off ⇒ the slot allocates `(1,1)` placeholders and its
      !! counted footprint is the placeholder size, so a build without a
      !! cavity pays nothing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_cavity_flux_t) :: cav

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      cav%enable = .false.
      call cav%init(grid)
      call check(error, size(cav%melt, 1) == 1 .and. size(cav%melt, 2) == 1, &
                 "knob off must leave a (1,1) placeholder")
      if (allocated(error)) go to 900
      call check(error, cav%is_init, "init must still mark the slot initialised")
      if (allocated(error)) go to 900
      call check(error, cav%bytes() < int(1000, kind(cav%bytes())), &
                 "a gated-off slot must count essentially zero bytes")
900   call cav%destroy()
   end subroutine test_off_placeholder

   ! ------------------------------------------------------------------
   ! P2c: the cover mask on atmospheric forcing
   ! ------------------------------------------------------------------
   !
   ! Under `cover_frac = 1` there is no atmosphere.  What that has to
   ! mean, field by field, is asserted below.  Two properties run through
   ! every one of them:
   !
   !   ABSENT ⇒ UNMASKED.  Every masked kernel is reached through an
   !   OPTIONAL `cover_frac`; with the argument absent the call dispatches
   !   to the untouched `_impl`, so a run without a cavity is byte-
   !   identical by construction.  `test_cover_absent` pins that the
   !   cover-aware twin agrees with it when the cover is all zero, which
   !   is the part construction does NOT give you for free.
   !
   !   THE CAVITY'S OWN FLUXES ARE NOT MASKED.  `heat_cavity` /
   !   `salt_cavity` are the melt interface, not the atmosphere.  A mask
   !   applied after the assembler would scale them away with everything
   !   else, which is why it is applied inside it.

   subroutine test_cover_absent(error)
      !! Cover present-but-zero reproduces cover-absent.  Asserted to a
      !! stated relative tolerance rather than bitwise: the two kernels
      !! are different expression trees (`open_f*(sum)` vs `sum`), and an
      !! FMA-contracting build is free to round them differently — even
      !! though `open_f` is exactly 1 here.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp), allocatable :: zero_cover(:, :), qh_ref(:, :), qs_ref(:, :)
      integer :: nx, ny

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, 0)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (zero_cover(nx, ny), source=0.0_wp)
      allocate (qh_ref(nx, ny), source=0.0_wp)
      allocate (qs_ref(nx, ny), source=0.0_wp)
      ! A live atmosphere so the comparison has something to compare.
      call sf%set_surface_flux_const(40.0_wp, 1.0e-5_wp)
      !$acc update device(sf%Q_heat, sf%Q_salt)
      !$acc enter data copyin(zero_cover) create(qh_ref, qs_ref)

      call ocean_surface_flux_assemble(grid, sf, ms)
      !$acc update self(sf%Q_heat, sf%Q_salt)
      qh_ref = sf%Q_heat
      qs_ref = sf%Q_salt

      call ocean_surface_flux_assemble(grid, sf, ms, cover_frac=zero_cover)
      !$acc update self(sf%Q_heat, sf%Q_salt)

      call check(error, maxval(abs(sf%Q_heat - qh_ref)) <= TOL_PATH*max(maxval(abs(qh_ref)), 1.0_wp), &
                 "a zero cover must reproduce the unmasked Q_heat")
      if (allocated(error)) go to 900
      call check(error, maxval(abs(sf%Q_salt - qs_ref)) <= TOL_PATH*max(maxval(abs(qs_ref)), 1.0_wp), &
                 "a zero cover must reproduce the unmasked Q_salt")
      if (allocated(error)) go to 900
      call check(error, maxval(abs(qh_ref)) > 0.0_wp, &
                 "the comparison must not be 0 == 0")

900   continue
      !$acc exit data delete(zero_cover, qh_ref, qs_ref)
      call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_cover_absent

   subroutine test_cover_stress(error)
      !! No wind acts on a covered cell, and KPP/EPBL see that.
      !!
      !! FACE RULE: a face is closed when EITHER neighbour is covered.
      !! Both of a covered cell's faces therefore close, so its
      !! `stress_mag` — the ONLY source of `u*` for KPP and EPBL — is
      !! EXACTLY zero.  The price, asserted here rather than hidden, is
      !! one face-wide transition: the first OPEN cell at the calving
      !! front keeps one live face and so reads half the open-ocean
      !! `|tau|`, which is the `u*` consistent with the momentum it
      !! actually received.  An open cell with two open neighbours is
      !! untouched, bitwise (its `tau` pair never changed, and
      !! `stress_mag` is recomputed by the same three lines).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), allocatable :: cover(:, :), mag0(:, :)
      real(wp), parameter :: TAUX = 0.2_wp, TAUY = 0.1_wp
      real(wp) :: full_mag
      integer :: nx, ny, i, j, i_ice, i_front, i_open

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      call ss%init(grid, NZ)
      call ss%set_wind_stress_const(TAUX, TAUY)
      full_mag = sqrt(TAUX*TAUX + TAUY*TAUY)

      ! Ice on the two westmost physical columns; open to the east.
      allocate (cover(nx, ny), source=0.0_wp)
      allocate (mag0(nx, ny), source=0.0_wp)
      i_ice = NGHOST + 1
      i_front = NGHOST + 3      ! first OPEN column, shares a face with ice
      i_open = NGHOST + 4       ! open with open neighbours on both sides
      do j = 1, ny
         cover(NGHOST + 1, j) = 1.0_wp
         cover(NGHOST + 2, j) = 1.0_wp
      end do
      mag0 = ss%stress_mag

      call ss%enter_data()
      !$acc enter data copyin(cover)
      call ocean_surface_stress_apply_cover(ss, cover)
      ! Idempotence: the mask multiplies by 0 or 1, so a second pass —
      ! which the per-bracket forcing seam would do — must change nothing.
      call ocean_surface_stress_apply_cover(ss, cover)
      !$acc update self(ss%tau_x, ss%tau_y, ss%stress_mag)

      j = NGHOST + 1
      call check(error, ss%stress_mag(i_ice, j) == 0.0_wp, &
                 "stress_mag must be EXACTLY zero under cover (KPP/EPBL u* = 0)")
      if (allocated(error)) go to 900
      call check(error, ss%tau_x(i_ice, j) == 0.0_wp .and. &
                 ss%tau_x(i_ice + 1, j) == 0.0_wp, &
                 "both u-faces of a covered cell must carry zero stress")
      if (allocated(error)) go to 900
      call check(error, ss%tau_y(i_ice, j) == 0.0_wp .and. &
                 ss%tau_y(i_ice, j + 1) == 0.0_wp, &
                 "both v-faces of a covered cell must carry zero stress")
      if (allocated(error)) go to 900
      call check(error, ss%stress_mag(i_open, j) == mag0(i_open, j), &
                 "an open cell away from the front must be bitwise untouched")
      if (allocated(error)) go to 900
      call check(error, rel_diff(mag0(i_open, j), full_mag) <= TOL_EXACT, &
                 "the untouched open cell must still carry the full |tau|")
      if (allocated(error)) go to 900
      ! The documented one-face transition, asserted as a fact rather
      ! than tolerated silently.
      call check(error, ss%stress_mag(i_front, j) > 0.0_wp .and. &
                 ss%stress_mag(i_front, j) < full_mag, &
                 "the first open cell at the front keeps exactly one live "// &
                 "u-face, so 0 < |tau| < full")
      if (allocated(error)) go to 900
      call check(error, ss%tau_x(i_front, j) == 0.0_wp, &
                 "the face SHARED with the ice must be closed (EITHER-neighbour rule)")

900   continue
      !$acc exit data delete(cover)
      call ss%exit_data()
      call ss%destroy()
   end subroutine test_cover_stress

   subroutine test_cover_face_rule(error)
      !! ONE FACE RULE, NOT TWO.
      !!
      !! The wind mask and the ice-shelf top drag both project the
      !! cell-centred `cover_frac` onto velocity faces, and if they ever
      !! disagreed the calving-front face would take the drag while still
      !! feeling the wind (or the reverse) — a slip/no-slip line that is
      !! exactly one face wide and in exactly the place the cavity
      !! outflow jet leaves.  The rule is STATED once, in
      !! `rdb_ocean_top_drag`'s module docstring; this is what keeps the
      !! two IMPLEMENTATIONS of it honest, the same way
      !! `mirror_of_bottom_drag` keeps the top and bottom drag one
      !! closure rather than two.
      !!
      !! Method: run `top_drag_fill_face_cover_impl` to get the drag's
      !! own face masks, then run the wind mask over a `tau` pair of
      !! ones and read back `1 - tau`, which IS the wind mask's face
      !! cover.  Asserted face-for-face over the whole plane, rim faces
      !! included — those are the only ones where the two could plausibly
      !! have drifted apart, since there is only one neighbour in range.
      !!
      !! A ragged cover (not a clean western block) so the comparison
      !! sees interior fronts in BOTH directions, not just one edge.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), allocatable :: cover(:, :), cu(:, :), cv(:, :)
      integer :: nx, ny, i, j, n_front
      real(wp) :: worst_u, worst_v

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      call ss%init(grid, NZ)
      ! tau = 1 everywhere, so `1 - tau` after masking IS the face cover.
      call ss%set_wind_stress_const(1.0_wp, 1.0_wp)

      allocate (cover(nx, ny), source=0.0_wp)
      allocate (cu(nx + 1, ny), source=-1.0_wp)
      allocate (cv(nx, ny + 1), source=-1.0_wp)
      do j = 1, ny
         do i = 1, nx
            ! Ragged: a western block plus a southern tongue, so fronts
            ! run in x AND y and touch the rim on two sides.
            if (i <= 3 .or. (j <= 2 .and. i <= 5)) cover(i, j) = 1.0_wp
         end do
      end do

      call top_drag_fill_face_cover_impl(cu, cv, cover, nx, ny)
      call ocean_surface_stress_apply_cover(ss, cover)

      worst_u = 0.0_wp
      do j = 1, ny
         do i = 1, nx + 1
            worst_u = max(worst_u, abs((1.0_wp - ss%tau_x(i, j)) - cu(i, j)))
         end do
      end do
      worst_v = 0.0_wp
      do j = 1, ny + 1
         do i = 1, nx
            worst_v = max(worst_v, abs((1.0_wp - ss%tau_y(i, j)) - cv(i, j)))
         end do
      end do

      call check(error, worst_u == 0.0_wp, &
                 "the wind mask and the top drag must agree on EVERY u-face "// &
                 "(one rule, two implementations)")
      if (allocated(error)) go to 900
      call check(error, worst_v == 0.0_wp, &
                 "the same on every v-face")
      if (allocated(error)) go to 900

      ! Non-vacuity: the plane must actually contain frontal faces, i.e.
      ! faces the OR rule closes and an AND rule would have left open.
      n_front = 0
      do j = 1, ny
         do i = 2, nx
            if (cover(i - 1, j) /= cover(i, j)) n_front = n_front + 1
         end do
      end do
      call check(error, n_front > 0, &
                 "the test plane must contain a calving front, else the two "// &
                 "rules agree trivially")
      if (allocated(error)) go to 900
      call check(error, cu(4, NGHOST + 1) == 1.0_wp, &
                 "a frontal u-face must be CLOSED by the shared OR rule")

900   continue
      call ss%destroy()
   end subroutine test_cover_face_rule

   subroutine test_cover_assemble(error)
      !! The assembler masks the ATMOSPHERE and passes the CAVITY
      !! through.  Under cover `Q_heat` must be exactly `heat_cavity` and
      !! `Q_salt` exactly `salt_cavity` — every atmospheric band gone,
      !! the melt flux intact.  In the open cells the answer must be the
      !! no-cavity answer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp), allocatable :: qh_open(:, :), qs_open(:, :)
      integer :: nx, ny, i_ice, i_open, j
      integer, parameter :: COVER_TO = 2

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, COVER_TO)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (qh_open(nx, ny), source=0.0_wp)
      allocate (qs_open(nx, ny), source=0.0_wp)
      i_ice = NGHOST + 1
      i_open = NGHOST + COVER_TO + 1
      j = NGHOST + 1

      ! A full atmosphere: the uniform scalar fill AND two component
      ! bands, so the test would catch a mask that only covered one of
      ! the two routes into Q_heat.
      call sf%set_surface_flux_const(40.0_wp, 2.0e-5_wp)
      sf%heat_added = -15.0_wp
      sf%q_lw = -30.0_wp
      sf%salt_flux = 5.0e-6_wp
      !$acc update device(sf%Q_heat, sf%Q_salt, sf%heat_added, sf%q_lw, sf%salt_flux)

      ! Reference: the SAME state assembled with no cover at all.
      call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
      call ocean_surface_flux_assemble(grid, sf, ms)
      !$acc update self(sf%Q_heat, sf%Q_salt)
      qh_open = sf%Q_heat
      qs_open = sf%Q_salt

      call ocean_surface_flux_assemble(grid, sf, ms, cover_frac=metrics%cover_frac)
      !$acc update self(sf%Q_heat, sf%Q_salt, sf%heat_cavity, sf%salt_cavity)
      !$acc update self(sf%heat_content_massin, sf%heat_content_massout)

      call check(error, sf%heat_cavity(i_ice, j) /= 0.0_wp, &
                 "the covered probe column must actually be melting")
      if (allocated(error)) go to 900
      call check(error, sf%Q_heat(i_ice, j) == sf%heat_cavity(i_ice, j), &
                 "under cover Q_heat must be EXACTLY heat_cavity — every "// &
                 "atmospheric band masked, the melt flux untouched")
      if (allocated(error)) go to 900
      call check(error, sf%Q_salt(i_ice, j) == sf%salt_cavity(i_ice, j), &
                 "under cover Q_salt must be EXACTLY salt_cavity")
      if (allocated(error)) go to 900
      call check(error, sf%heat_content_massin(i_ice, j) == 0.0_wp .and. &
                 sf%heat_content_massout(i_ice, j) == 0.0_wp, &
                 "the mass-enthalpy ledger must read zero under cover, so it "// &
                 "adds up against the masked Q_heat")
      if (allocated(error)) go to 900
      ! Open water: unchanged, to the two-expression-trees tolerance.
      call check(error, sf%heat_cavity(i_open, j) == 0.0_wp, &
                 "an uncovered column has no melt flux to protect")
      if (allocated(error)) go to 900
      call check(error, rel_diff(sf%Q_heat(i_open, j), qh_open(i_open, j)) <= TOL_PATH, &
                 "an open column must get the no-cavity Q_heat")
      if (allocated(error)) go to 900
      call check(error, rel_diff(sf%Q_salt(i_open, j), qs_open(i_open, j)) <= TOL_PATH, &
                 "an open column must get the no-cavity Q_salt")
      if (allocated(error)) go to 900
      call check(error, abs(qh_open(i_open, j)) > 0.0_wp, &
                 "the open-water comparison must not be 0 == 0")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_cover_assemble

   subroutine test_cover_sw_restore(error)
      !! Sunlight and an atmospheric restoring target both stop at the
      !! ice base.  A covered column's `hTr` and its surface budget
      !! contributor must come out EXACTLY unchanged (the increment is a
      !! literal `+ 0.0`, not a small number), while an open column moves
      !! by the same amount it would without a cavity.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp), allocatable :: t_before(:, :, :), t_ref(:, :, :)
      integer :: nx, ny, nzl, i_ice, i_open, j, k
      integer, parameter :: COVER_TO = 2
      real(wp), parameter :: DT = 3600.0_wp

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, COVER_TO)
      nx = grid%nx_total
      ny = grid%ny_total
      nzl = NZ
      i_ice = NGHOST + 1
      i_open = NGHOST + COVER_TO + 1
      j = NGHOST + 1
      allocate (t_before(nx, ny, nzl), source=0.0_wp)
      allocate (t_ref(nx, ny, nzl), source=0.0_wp)

      ! A strong sun and a cold restoring target: either one leaking
      ! under the shelf would be unmissable.
      call sf%set_surface_flux_const(200.0_wp, 0.0_wp)
      call sf%set_sw_penetration(0.5_wp, 0.58_wp, 0.35_wp, 23.0_wp)
      call sf%set_restore(.true., .false., 10.0_wp, 0.0_wp, -20.0_wp, 0.0_wp)
      !$acc update device(sf%Q_heat, sf%Q_salt)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      t_before = ms%tracers(ms%idx_temperature)%hTr

      ! Reference pass with NO cover, on a scratch copy of the state.
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, DT)
      call ocean_surface_restore_apply_tracers(grid, sf, ms, DT)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      t_ref = ms%tracers(ms%idx_temperature)%hTr

      ! Restore the pre-pass state and redo it WITH the cover.
      ms%tracers(ms%idx_temperature)%hTr = t_before
      !$acc update device(ms%tracers(ms%idx_temperature)%hTr)
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, DT, &
                                                   cover_frac=metrics%cover_frac)
      call ocean_surface_restore_apply_tracers(grid, sf, ms, DT, &
                                               cover_frac=metrics%cover_frac)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)

      do k = 1, nzl
         call check(error, ms%tracers(ms%idx_temperature)%hTr(i_ice, j, k) == &
                    t_before(i_ice, j, k), &
                    "a covered column must take EXACTLY zero shortwave and zero "// &
                    "restoring, at every layer")
         if (allocated(error)) go to 900
      end do
      call check(error, t_ref(i_ice, j, nzl) /= t_before(i_ice, j, nzl), &
                 "the unmasked reference must have moved the covered column, "// &
                 "else the assertion above is vacuous")
      if (allocated(error)) go to 900
      call check(error, ms%tracers(ms%idx_temperature)%hTr(i_open, j, nzl) == &
                 t_ref(i_open, j, nzl), &
                 "an open column must get the no-cavity answer, bitwise (the "// &
                 "cover factor there is exactly 1)")
      if (allocated(error)) go to 900
      call check(error, t_ref(i_open, j, nzl) /= t_before(i_open, j, nzl), &
                 "the open column must actually have been forced")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_cover_sw_restore

   subroutine test_cover_budget(error)
      !! Heat and salt still close to round-off with BOTH sources live —
      !! an atmosphere over the open cells and a melting interface under
      !! the covered ones.  This is the composition the lifted refusals
      !! made legal, so it is the one that has to be shown conservative.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_cavity_flux_t) :: cav
      type(eos_t) :: eos
      real(wp) :: heat0, salt0, heat1, salt1, bud_h, bud_s, scale
      real(wp) :: bud_ice, bud_open
      integer :: nx, ny, step, i_ice, i_open, j
      integer, parameter :: N_STEP = 5, COVER_TO = 2
      real(wp), parameter :: DT = 1800.0_wp

      call build_cavity_plane(grid, metrics, ms, sf, cav, eos, P_DRAFT, T_WARM, &
                              0.05_wp, COVER_TO)
      nx = grid%nx_total
      ny = grid%ny_total
      scale = grid%dx*grid%dy*RHO_WATER
      i_ice = NGHOST + 1
      i_open = NGHOST + COVER_TO + 1
      j = NGHOST + 1

      call sf%set_surface_flux_const(60.0_wp, 3.0e-5_wp)
      !$acc update device(sf%Q_heat, sf%Q_salt)

      heat0 = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                     NGHOST + 1:ny - NGHOST, :))*scale
      salt0 = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                  NGHOST + 1:ny - NGHOST, :))*scale

      do step = 1, N_STEP
         call ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf)
         call ocean_surface_flux_assemble(grid, sf, ms, cover_frac=metrics%cover_frac)
         call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
      end do

      !$acc update self(ms%tracers(ms%idx_temperature)%hTr, &
      !$acc&            ms%tracers(ms%idx_salinity)%hTr, &
      !$acc&            ms%heat_budget_surface, ms%salt_budget_surface)

      heat1 = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                     NGHOST + 1:ny - NGHOST, :))*scale
      salt1 = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                  NGHOST + 1:ny - NGHOST, :))*scale
      bud_h = sum(ms%heat_budget_surface(NGHOST + 1:nx - NGHOST, &
                                         NGHOST + 1:ny - NGHOST, :))*scale
      bud_s = sum(ms%salt_budget_surface(NGHOST + 1:nx - NGHOST, &
                                         NGHOST + 1:ny - NGHOST, :))*scale
      bud_ice = sum(ms%heat_budget_surface(i_ice, j, :))
      bud_open = sum(ms%heat_budget_surface(i_open, j, :))

      ! Non-vacuity: the covered column must be COOLING (cavity only) and
      ! the open one WARMING (atmosphere only).  Opposite signs are the
      ! cheapest proof that the mask separated the two.
      call check(error, bud_ice < 0.0_wp, &
                 "the covered column's only heat source is the melt interface, "// &
                 "which cools it")
      if (allocated(error)) go to 900
      call check(error, bud_open > 0.0_wp, &
                 "the open column's only heat source is the +60 W/m^2 atmosphere")
      if (allocated(error)) go to 900
      ! Normalised against the INITIAL TOTAL, exactly as the console's
      ! `relative_drift` does and for the reason documented in
      ! `test_budget_closes`.
      call check(error, abs((heat1 - heat0) - bud_h) <= 1.0e-13_wp*abs(heat0), &
                 "heat must close to round-off with atmosphere AND cavity live")
      if (allocated(error)) go to 900
      call check(error, abs((salt1 - salt0) - bud_s) <= 1.0e-13_wp*abs(salt0), &
                 "salt must close to round-off with atmosphere AND cavity live")

900   call teardown_cavity_plane(metrics, ms, sf, cav)
   end subroutine test_cover_budget

   ! ------------------------------------------------------------------
   ! Refusal matrix
   ! ------------------------------------------------------------------

   subroutine test_validate_refusals(error)
      !! Every v1 restriction fails loud at configure, and the refusal
      !! comes from `validate_config` — NOT from the strict schema
      !! tripping over a malformed test namelist.  Each case is asserted
      !! to PARSE cleanly first, so a silently mis-built namelist cannot
      !! make a refusal test pass for the wrong reason.
      !!
      !! Four entries here used to be refusals and are now ACCEPTANCES:
      !! wind stress, surface restoring, shortwave penetration and the
      !! uniform scalar `q_heat`/`q_salt` were all refused while there
      !! was no per-cell cover mask on the atmospheric forcing.  P2c
      !! ships one, so the verdict flips — and the acceptance is asserted
      !! rather than deleted, so a regression that re-introduced the
      !! refusal has to delete a test to pass.  The MASKING itself is
      !! asserted numerically in the `cavity_cover_*` suite above; this
      !! test only pins the configure-time verdict.
      !!
      !! What remains a refusal, and why, is the file-driven data
      !! override: it rewrites `tau` (and `Q_heat`) per time bracket
      !! without the cover, so it would silently restore the unmasked
      !! atmosphere under the shelf.  The interface pressure is NOT a
      !! gap: the cavity assembles `ms%p_top = p_ice_ref + sf%p_surf`, so
      !! `&ocean_psurf_nml` is asserted ACCEPTED below.
      type(error_type), allocatable, intent(out) :: error

      ! --- the prerequisites ---
      call expect_refused("melt without the cavity geometry group", &
                          melt_nml("enable = .true.", "isomip", .true., "", .false.), error)
      if (allocated(error)) return
      call expect_refused("melt with the sea-ice liquidus", &
                          melt_nml("enable = .true.", "seaice", .true., "", .true.), error)
      if (allocated(error)) return
      call expect_refused("melt without the surface-flux component set", &
                          melt_nml("enable = .true.", "isomip", .false., "", .true.), error)
      if (allocated(error)) return

      ! --- the enum seams (RESERVED must be distinguishable from a typo) ---
      call expect_refused("a RESERVED exchange law", &
                          melt_nml("enable = .true., exchange_law = 'burchard22'", &
                                   "isomip", .true., "", .true.), error)
      if (allocated(error)) return
      call expect_refused("the RESERVED diffusive ice mode", &
                          melt_nml("enable = .true., ice_conduction = 'diffusive'", &
                                   "isomip", .true., "", .true.), error)
      if (allocated(error)) return
      call expect_refused("hj99 on an f = 0 grid", &
                          melt_nml("enable = .true., exchange_law = 'hj99'", &
                                   "isomip", .true., "", .true.), error)
      if (allocated(error)) return

      ! --- ranges ---
      call expect_refused("gamma_s = 0 (the three-equation form divides by it)", &
                          melt_nml("enable = .true., gamma_s = 0.0", &
                                   "isomip", .true., "", .true.), error)
      if (allocated(error)) return

      ! --- the cover mask SHIPS (P2c): four ACCEPTANCES that used to be
      !     refusals.  Kept as explicit acceptances rather than deleted,
      !     so a regression that re-introduced the refusal — or that
      !     quietly dropped the mask and had to put it back — has to
      !     delete a test to pass.  The masking ITSELF is asserted
      !     numerically in `cavity_cover_*` below; this block only pins
      !     the configure-time verdict.
      call expect_accepted("a non-zero wind stress, now cover-masked", &
                           melt_nml("enable = .true.", "isomip", .true., &
                                    "&physics_nml wind_stress_x = 0.05 /", .true.), error)
      if (allocated(error)) return
      call expect_accepted("surface restoring, now cover-masked", &
                           melt_nml("enable = .true.", "isomip", .true., &
                                    "&ocean_restore_nml enable_restore_temp = .true., "// &
                                    "piston_t = 1.0 /", .true.), error)
      if (allocated(error)) return
      call expect_accepted("shortwave penetration, now cover-masked", &
                           melt_nml("enable = .true.", "isomip", .true., &
                                    "&ocean_thermo_nml sw_pen_frac = 0.3 /", .true.), error)
      if (allocated(error)) return
      call expect_accepted("a uniform scalar surface flux, now cover-masked", &
                           melt_nml("enable = .true.", "isomip", .true., &
                                    "&ocean_thermo_nml q_heat = 10.0 /", .true.), error)
      if (allocated(error)) return

      ! --- the ONE atmospheric-forcing path the mask does not reach ---
      ! `ocean_data_forcing_apply` rewrites tau_x/tau_y (and Q_heat) per
      ! time bracket and is handed no metrics slot, so it would restore
      ! the unmasked atmosphere under the shelf.  Refused on the GEOMETRY
      ! group, because that is where the cover comes from.
      call expect_refused("the file-driven data override, which re-writes tau "// &
                          "per bracket with no cover", &
                          melt_nml("enable = .true.", "isomip", .true., &
                                   "&ocean_dataovr_nml enable = .true. /", .true.), error)
      if (allocated(error)) return

      ! --- NOT a gap: the psurf seam composes, it does not clobber ---
      ! `ms%p_top` is assembled as `p_ice_ref + sf%p_surf` by the cavity
      ! (P5.2), so the seam can no longer overwrite the ice load with the
      ! atmospheric one alone and the combination is ACCEPTED.  Asserted
      ! as an acceptance rather than deleted, so a future regression that
      ! re-introduced the clobber would have to delete a test to pass.
      call expect_accepted("cavity melt alongside the psurf seam", &
                           melt_nml("enable = .true.", "isomip", .true., &
                                    "&ocean_psurf_nml enable = .true. /", .true.), error)
      if (allocated(error)) return

      ! ... and the in-envelope configuration is ACCEPTED, so the matrix
      ! above is not simply refusing everything.
      call expect_accepted("the shipped in-envelope cavity-melt run", &
                           melt_nml("enable = .true.", "isomip", .true., "", .true.), error)
      if (allocated(error)) return
      call expect_accepted("the geometry-only cavity run (melt off)", &
                           melt_nml("", "isomip", .true., "", .true.), error)
   end subroutine test_validate_refusals

   function melt_nml(melt_body, tfreeze, components, extra, with_geometry) result(nml)
      !! Build one complete test namelist.  Every group appears EXACTLY
      !! ONCE — the strict schema refuses a duplicate group, and a
      !! refusal test that trips on that would be testing the parser, not
      !! the rule it names.
      character(len=*), intent(in) :: melt_body
         !! `&ocean_cavity_melt_nml` body without the group name or `/`;
         !! empty ⇒ the group is omitted entirely.
      character(len=*), intent(in) :: tfreeze
         !! `&ocean_eos_nml tfreeze_set` value.
      logical, intent(in) :: components
         !! `&ocean_forcing_nml enable_components`.
      character(len=*), intent(in) :: extra
         !! One extra complete group (with name and `/`), or empty.
      logical, intent(in) :: with_geometry
         !! Append `&ocean_cavity_dyn_nml` (the melt group's prerequisite).
      character(len=:), allocatable :: nml, comp

      comp = ".false."
      if (components) comp = ".true."
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 8 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")// &
            "&ocean_eos_nml tfreeze_set = '"//tfreeze//"' /"//new_line("a")// &
            "&ocean_forcing_nml enable_components = "//comp//" /"//new_line("a")
      if (with_geometry) then
         nml = nml//"&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', "// &
               "draft_depth = 300.0 /"//new_line("a")
      end if
      if (len_trim(melt_body) > 0) then
         nml = nml//"&ocean_cavity_melt_nml "//melt_body//" /"//new_line("a")
      end if
      if (len_trim(extra) > 0) nml = nml//extra//new_line("a")
   end function melt_nml

   subroutine expect_refused(what, nml, error)
      !! The namelist must PARSE and then be refused by `validate_config`.
      character(len=*), intent(in) :: what, nml
      type(error_type), allocatable, intent(inout) :: error
      type(config_t) :: cfg
      integer :: ierr

      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 "the refusal case must PARSE cleanly (else it tests the schema, "// &
                 "not the rule): "//what)
      if (allocated(error)) return
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "validate_config must refuse: "//what)
   end subroutine expect_refused

   subroutine expect_accepted(what, nml, error)
      character(len=*), intent(in) :: what, nml
      type(error_type), allocatable, intent(inout) :: error
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "must parse: "//what)
      if (allocated(error)) return
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "must be accepted: "//what)
   end subroutine expect_accepted

end module test_ocean_cavity_flux
