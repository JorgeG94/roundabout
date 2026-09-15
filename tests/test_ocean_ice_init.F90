!! Analytic tests for the sea-ice ANALYTIC initial-condition path (PR 24):
!! `rdb_ice_init` (`ice_ic_parse_conc_config`/`ice_ic_params_from_config`/
!! `ice_ic_target_category`/`ice_init_apply`).
!!
!! Harness mirrors `test_ocean_ice_itd.F90` (grid + `multilayer_state_t`
!! + `ocean_sea_ice_t` built the "driver way": `enable`/`ncat`/`nk_ice` set
!! BEFORE `ice%init`) plus `ocean_test_metrics` for the `ocean_metrics_t`
!! slot (`geolatT`, needed by `"latitudes"`).  `ice_init_apply` itself is
!! HOST-side pure code — never mapped, never called after `enter_data` — so
!! most cases need no `!$acc` discipline at all.  The two cases that DO
!! touch a device kernel (case 4: `ice_adjust_categories`; case 10:
!! `ice_evp_step`) follow the `mem:separate` template
!! (`test_open_boundary_out_closes`, `test_ocean_conservation_salt_heat.F90`):
!! unconditional `!$acc enter data copyin(...)` before the kernel, `!$acc
!! update self` before any host assertion, `exit data` at teardown.
!!
!! `ice%mh_lim` is filled by `ice%init` (`ice_itd_category_bounds`) —
!! PR-58's `&ocean_ice_nml hlim` override is NOT merged on this branch, so
!! case 3 hand-seeds a non-default `ice%mh_lim` after `init` and re-runs the
!! binning to prove `ice_ic_target_category`/`ice_init_apply` read the LIVE
!! `ice%mh_lim` field, never a hardcoded/re-derived table (the PR-58 seam,
!! §13.2 of the plan).
!!
!! Cases (PLAN_PR24_ice_ic_path.md §9):
!!   1  ic_disabled_bitident         — conc_config="zero" touches nothing.
!!   2  ic_uniform_mass_and_area     — mass/area unit conversion, ncat>1.
!!   3  ic_category_binning          — ITD binning + the mh_lim seam.
!!   4  ic_is_itd_fixed_point        — ice_adjust_categories is a no-op.
!!   5  ic_enthalpy_roundtrip        — exact T<->E inversion + liquid guard.
!!   6  ic_land_and_ghosts           — land untouched, ghosts DO seed.
!!   7  ic_ncat1_lumped              — legacy per-cell mode, part_size frozen.
!!   8  ic_latitudes                 — SIS2 polar-cap 0/1 step.
!!   9  ic_restart_roundtrip_uniform — bit-exact restart round-trip.
!!  10  ic_nansen_free_drift_from_ic — IC -> ice_evp_step -> Nansen golden.
module test_ocean_ice_init
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ice_state, only: ocean_sea_ice_t, ice_cell_concentration_impl
   use rdb_ice_itd, only: ice_adjust_categories
   use rdb_ice_evp, only: ice_evp_step, ice_evp_params_t
   use rdb_ice_column, only: ICE_RHO_ICE, ICE_RHO_SNOW, ICE_BULK_SALINITY
   use rdb_ice_enthalpy, only: ice_enth_from_ts, ice_temp_from_en_s, ice_enthalpy_liquid_freeze
   use rdb_ice_init, only: ice_ic_params_t, ice_ic_params_from_config, ice_ic_parse_conc_config, &
                           ice_ic_target_category, ice_init_apply, &
                           ICE_IC_CONC_ZERO, ICE_IC_CONC_UNIFORM, ICE_IC_CONC_LATITUDES, &
                           ICE_IC_CONC_INVALID
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data, &
                              ocean_state_restart_write, ocean_state_restart_read
   use rdb_decomp, only: decomp_t, decomp_init
   implicit none
   private

   public :: collect_ocean_ice_init_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 4
   real(wp), parameter :: DX = 1000.0_wp

contains

   subroutine collect_ocean_ice_init_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ic_disabled_bitident", test_disabled_bitident), &
                  new_unittest("ic_uniform_mass_and_area", test_uniform_mass_and_area), &
                  new_unittest("ic_category_binning", test_category_binning), &
                  new_unittest("ic_is_itd_fixed_point", test_is_itd_fixed_point), &
                  new_unittest("ic_enthalpy_roundtrip", test_enthalpy_roundtrip), &
                  new_unittest("ic_land_and_ghosts", test_land_and_ghosts), &
                  new_unittest("ic_ncat1_lumped", test_ncat1_lumped), &
                  new_unittest("ic_latitudes", test_latitudes), &
                  new_unittest("ic_restart_roundtrip_uniform", test_restart_roundtrip_uniform), &
                  new_unittest("ic_nansen_free_drift_from_ic", test_nansen_free_drift_from_ic) &
                  ]
   end subroutine collect_ocean_ice_init_tests

   ! -----------------------------------------------------------------
   ! Shared setup / teardown
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, ms, ice, ncat, nk_ice)
      !! Tiny all-wet grid + ice slot built the "driver way": `enable`/
      !! `ncat`/`nk_ice` latched BEFORE `ice%init` (mirrors
      !! `test_ocean_ice_itd.F90::setup_state`).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: ncat, nk_ice

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      ms%nz_ml = 1
      call ms%init(grid)
      ice%enable = .true.
      ice%ncat = ncat
      ice%nk_ice = nk_ice
      call ice%init(grid)
   end subroutine setup_state

   subroutine teardown(ms, ice, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_metrics_t), intent(inout), optional :: metrics
      call ice%destroy()
      call ms%destroy()
      if (present(metrics)) call metrics%destroy()
   end subroutine teardown

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_disabled_bitident(error)
      !! Case 1: `conc_config="zero"` (the default) leaves every ice array
      !! exactly at its `ocean_sea_ice_init` value.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      checks: block

         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)

         ! par is default-constructed => conc_config = ICE_IC_CONC_ZERO.
         call check(error, par%conc_config == ICE_IC_CONC_ZERO, &
                    "default ice_ic_params_t must be conc_config=ZERO")
         if (allocated(error)) exit checks

         call ice_init_apply(grid, ms, metrics, ice, par)

         call check(error, all(ice%part_size(:, :, 0) == 1.0_wp), &
                    "zero: part_size(:,:,0) must stay 1 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ice%part_size(:, :, 1:) == 0.0_wp), &
                    "zero: part_size(:,:,1:) must stay 0")
         if (allocated(error)) exit checks
         call check(error, all(ice%m_ice == 0.0_wp), "zero: m_ice must stay 0")
         if (allocated(error)) exit checks
         call check(error, all(ice%m_snow == 0.0_wp), "zero: m_snow must stay 0")
         if (allocated(error)) exit checks
         call check(error, all(ice%enth_ice == 0.0_wp), "zero: enth_ice must stay 0")
         if (allocated(error)) exit checks
         call check(error, all(ice%enth_snow == 0.0_wp), "zero: enth_snow must stay 0")
         if (allocated(error)) exit checks
         call check(error, all(ice%sal_ice == ICE_BULK_SALINITY), &
                    "zero: sal_ice must stay ICE_BULK_SALINITY")

      end block checks
      call teardown(ms, ice, metrics)
   end subroutine test_disabled_bitident

   subroutine test_uniform_mass_and_area(error)
      !! Case 2: ncat=5, h_ice=2.0, h_snow=0.1, conc=0.8 => the per-ICE-area
      !! convention (m_ice = RHO*h_ice, NOT x conc), landing in cat 5 (top,
      !! unbounded bin: m=1810 >= mh_lim(5)=995.5).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      real(wp), allocatable :: ci(:, :), mis(:, :), mice(:, :)
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)
         par = ice_ic_params_from_config("uniform", 0.8_wp, 2.0_wp, 0.1_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)

         call ice_init_apply(grid, ms, metrics, ice, par)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         call check(error, abs(ice%part_size(ip, jp, 0) - 0.2_wp) < 1.0e-14_wp, &
                    "uniform: part_size(0) must be 1-conc")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 5) - 0.8_wp) < 1.0e-14_wp, &
                    "uniform: part_size(5) must be conc (top unbounded bin)")
         if (allocated(error)) exit checks
         call check(error, all(ice%part_size(ip, jp, 1:4) == 0.0_wp), &
                    "uniform: part_size(1:4) must stay 0")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(ip, jp, 5) == ICE_RHO_ICE*2.0_wp, &
                    "uniform: m_ice(5) must equal RHO_ICE*h_ice EXACTLY "// &
                    "(per-ICE-area, NOT x conc)")
         if (allocated(error)) exit checks
         call check(error, all(ice%m_ice(ip, jp, 1:4) == 0.0_wp), &
                    "uniform: m_ice(1:4) must stay 0")
         if (allocated(error)) exit checks
         call check(error, ice%m_snow(ip, jp, 5) == ICE_RHO_SNOW*0.1_wp, &
                    "uniform: m_snow(5) must equal RHO_SNOW*h_snow")

         if (allocated(error)) exit checks

         allocate (ci(grid%nx_total, grid%ny_total))
         allocate (mis(grid%nx_total, grid%ny_total))
         allocate (mice(grid%nx_total, grid%ny_total))
         call ice_cell_concentration_impl(metrics%wet_T, ice%part_size, ice%m_ice, &
                                          ice%m_snow, mis, mice, ci, ice%ncat, &
                                          grid%nx_total, grid%ny_total)
         call check(error, abs(ci(ip, jp) - 0.8_wp) < 1.0e-14_wp, &
                    "uniform: ice_cell_concentration_impl ci must equal conc")
         if (allocated(error)) exit checks
         call check(error, abs(mice(ip, jp) - 0.8_wp*ICE_RHO_ICE*2.0_wp) <= 1.0e-9_wp, &
                    "uniform: ice_cell_concentration_impl mice must recover "// &
                    "the grid-mean mass (0.8*1810)")

      end block checks
      if (allocated(ci)) deallocate (ci)
      if (allocated(mis)) deallocate (mis)
      if (allocated(mice)) deallocate (mice)
      call teardown(ms, ice, metrics)
   end subroutine test_uniform_mass_and_area

   subroutine test_category_binning(error)
      !! Case 3: sweep h_ice at ncat=5 => occupied category {1,2,3,4,5,5};
      !! every occupied category within its `mh_lim` bin. Then (the PR-58
      !! seam): overwrite `ice%mh_lim` with a non-default table and re-bin
      !! the SAME thickness — the target category must move, proving
      !! `ice_ic_target_category`/`ice_init_apply` read the LIVE
      !! `ice%mh_lim`, never a hardcoded/re-derived copy.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NSWEEP = 6
      real(wp), parameter :: H_SWEEP(NSWEEP) = [0.05_wp, 0.2_wp, 0.5_wp, 0.9_wp, 1.3_wp, 2.0_wp]
      integer, parameter :: CAT_EXPECT(NSWEEP) = [1, 2, 3, 4, 5, 5]
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      real(wp) :: m_target, mh_lim_override(6)
      integer :: n, ip, jp, c_star, cat
      checks: block

         do n = 1, NSWEEP
            call setup_state(grid, ms, ice, 5, 2)
            call metrics%init(grid)
            par = ice_ic_params_from_config("uniform", 1.0_wp, H_SWEEP(n), 0.0_wp, &
                                            -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
            call ice_init_apply(grid, ms, metrics, ice, par)

            ip = grid%nx_total/2
            jp = grid%ny_total/2
            m_target = ICE_RHO_ICE*H_SWEEP(n)

            cat = 0
            do c_star = 1, 5
               if (ice%part_size(ip, jp, c_star) > 0.0_wp) cat = c_star
            end do
            call check(error, cat == CAT_EXPECT(n), "category_binning: wrong occupied category")
            if (allocated(error)) exit checks
            call check(error, ice%m_ice(ip, jp, cat) >= ice%mh_lim(cat), &
                       "category_binning: below the lower bin bound")
            if (allocated(error)) exit checks
            if (cat < ice%ncat) then
               call check(error, ice%m_ice(ip, jp, cat) < ice%mh_lim(cat + 1), &
                          "category_binning: above the upper bin bound")
               if (allocated(error)) exit checks
            end if
            ! ice_ic_target_category must agree with what ice_init_apply did.
            call check(error, ice_ic_target_category(m_target, ice%mh_lim, ice%ncat) == cat, &
                       "category_binning: ice_ic_target_category disagrees with ice_init_apply")
            if (allocated(error)) exit checks

            call teardown(ms, ice, metrics)
         end do

         ! ---- The PR-58 seam: hand-seed a NON-DEFAULT mh_lim, prove the
         ! binning reads it live. h=0.9 defaults to cat 4 (mh_lim(4)=0.7*905
         ! =633.5 <= m=814.5 < mh_lim(5)=1.1*905=995.5, see sweep above).
         ! Override mh_lim so the SAME mass now falls in cat 2 instead.
         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)
         m_target = ICE_RHO_ICE*0.9_wp
         ! Default table would give cat 4; override so cat 2 catches it:
         ! mh_lim = [0, 100, 900, 1000, 2000, 3000] (kg/m^2) -- 814.5 in [100,900).
         mh_lim_override = [0.0_wp, 100.0_wp, 900.0_wp, 1000.0_wp, 2000.0_wp, 3000.0_wp]
         ice%mh_lim(:) = mh_lim_override(:)
         par = ice_ic_params_from_config("uniform", 1.0_wp, 0.9_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
         call ice_init_apply(grid, ms, metrics, ice, par)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, ice%part_size(ip, jp, 2) == 1.0_wp .and. &
                    all(ice%part_size(ip, jp, [1, 3, 4, 5]) == 0.0_wp), &
                    "category_binning: with an overridden ice%mh_lim, h=0.9 "// &
                    "must land in cat 2 (NOT cat 4, the default-table answer) "// &
                    "-- the IC must read ice%mh_lim, never a hardcoded table")

      end block checks
      call teardown(ms, ice, metrics)
   end subroutine test_category_binning

   subroutine test_is_itd_fixed_point(error)
      !! Case 4: for the same h_ice sweep, `ice_adjust_categories` (the
      !! PRODUCTION kernel, used as a TEST ORACLE) applied to the IC-seeded
      !! state is a bit-exact no-op. Device kernel => mem:separate
      !! discipline.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NSWEEP = 6
      real(wp), parameter :: H_SWEEP(NSWEEP) = [0.05_wp, 0.2_wp, 0.5_wp, 0.9_wp, 1.3_wp, 2.0_wp]
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      real(wp), allocatable :: ps0(:, :, :), mi0(:, :, :), msn0(:, :, :)
      real(wp), allocatable :: ei0(:, :, :, :), es0(:, :, :, :), si0(:, :, :, :)
      integer :: n
      checks: block

         do n = 1, NSWEEP
            call setup_state(grid, ms, ice, 5, 2)
            call metrics%init(grid)
            par = ice_ic_params_from_config("uniform", 1.0_wp, H_SWEEP(n), 0.05_wp, &
                                            -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
            call ice_init_apply(grid, ms, metrics, ice, par)

            allocate (ps0, source=ice%part_size)
            allocate (mi0, source=ice%m_ice)
            allocate (msn0, source=ice%m_snow)
            allocate (ei0, source=ice%enth_ice)
            allocate (es0, source=ice%enth_snow)
            allocate (si0, source=ice%sal_ice)

            !$acc enter data copyin(ms)
            call ms%enter_data()
            call ice%enter_data()
            call ice_adjust_categories(grid, ms, ice)
            associate (ps => ice%part_size, mi => ice%m_ice, msn => ice%m_snow, &
                       ei => ice%enth_ice, es => ice%enth_snow, si => ice%sal_ice)
               !$acc update self(ps, mi, msn, ei, es, si)
            end associate
            call ice%exit_data()
            call ms%exit_data()
            !$acc exit data delete(ms)

            call check(error, maxval(abs(ice%part_size - ps0)) == 0.0_wp, &
                       "itd_fixed_point: part_size must be bit-exact")
            if (.not. allocated(error)) &
               call check(error, maxval(abs(ice%m_ice - mi0)) == 0.0_wp, &
                          "itd_fixed_point: m_ice must be bit-exact")
            if (.not. allocated(error)) &
               call check(error, maxval(abs(ice%m_snow - msn0)) == 0.0_wp, &
                          "itd_fixed_point: m_snow must be bit-exact")
            if (.not. allocated(error)) &
               call check(error, maxval(abs(ice%enth_ice - ei0)) == 0.0_wp, &
                          "itd_fixed_point: enth_ice must be bit-exact")
            if (.not. allocated(error)) &
               call check(error, maxval(abs(ice%enth_snow - es0)) == 0.0_wp, &
                          "itd_fixed_point: enth_snow must be bit-exact")
            if (.not. allocated(error)) &
               call check(error, maxval(abs(ice%sal_ice - si0)) == 0.0_wp, &
                          "itd_fixed_point: sal_ice must be bit-exact")

            deallocate (ps0, mi0, msn0, ei0, es0, si0)
            call teardown(ms, ice, metrics)
            if (allocated(error)) exit checks
         end do

      end block checks
      if (allocated(ps0)) deallocate (ps0)
      if (allocated(mi0)) deallocate (mi0)
      if (allocated(msn0)) deallocate (msn0)
      if (allocated(ei0)) deallocate (ei0)
      if (allocated(es0)) deallocate (es0)
      if (allocated(si0)) deallocate (si0)
      if (ice%is_init) call teardown(ms, ice, metrics)
   end subroutine test_is_itd_fixed_point

   subroutine test_enthalpy_roundtrip(error)
      !! Case 5: exact T<->E inversion for every seeded (i,j,c,k), plus the
      !! liquid-branch guard: enth_ice must be strictly below the
      !! liquid-freeze enthalpy (i.e. it really is ice, not water).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: T_ICE = -4.0_wp
      real(wp), parameter :: S_ICE = 4.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      integer :: ip, jp, cat, k
      real(wp) :: t_back
      checks: block

         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)
         par = ice_ic_params_from_config("uniform", 1.0_wp, 2.0_wp, 0.5_wp, &
                                         T_ICE, S_ICE, 91.0_wp, -91.0_wp)
         call ice_init_apply(grid, ms, metrics, ice, par)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         cat = 5   ! h=2.0 -> top bin at ncat=5.

         do k = 1, ice%nk_ice
            t_back = ice_temp_from_en_s(ice%enth_ice(ip, jp, cat, k), ice%sal_ice(ip, jp, cat, k))
            call check(error, abs(t_back - T_ICE) < 1.0e-12_wp, &
                       "enthalpy_roundtrip: T<->E inversion must be exact for enth_ice")
            if (allocated(error)) exit checks
         end do
         t_back = ice_temp_from_en_s(ice%enth_snow(ip, jp, cat, 1), 0.0_wp)
         call check(error, abs(t_back - T_ICE) < 1.0e-12_wp, &
                    "enthalpy_roundtrip: T<->E inversion must be exact for enth_snow (fresh)")
         if (allocated(error)) exit checks
         call check(error, all(ice%enth_ice(ip, jp, cat, :) < ice_enthalpy_liquid_freeze(S_ICE)), &
                    "enthalpy_roundtrip: enth_ice must be strictly below the liquid-freeze "// &
                    "enthalpy (this really is ICE, not water)")

      end block checks
      call teardown(ms, ice, metrics)
   end subroutine test_enthalpy_roundtrip

   subroutine test_land_and_ghosts(error)
      !! Case 6: an interior land block stays untouched; ghost rows/columns
      !! ARE seeded identically to the interior (uniform mode) — the
      !! "formula bathymetry setters must fill ghost rows" gotcha wearing a
      !! different hat.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      integer :: il, jl, ig, jg
      checks: block

         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)

         ! Interior land block: one cell, dead centre.
         il = grid%nx_total/2
         jl = grid%ny_total/2
         ms%wet_mask(il, jl) = 0.0_wp

         par = ice_ic_params_from_config("uniform", 1.0_wp, 2.0_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
         call ice_init_apply(grid, ms, metrics, ice, par)

         call check(error, ice%part_size(il, jl, 0) == 1.0_wp .and. &
                    all(ice%part_size(il, jl, 1:) == 0.0_wp), &
                    "land_and_ghosts: land cell must stay open-water")
         if (allocated(error)) exit checks
         call check(error, all(ice%m_ice(il, jl, :) == 0.0_wp) .and. &
                    all(ice%m_snow(il, jl, :) == 0.0_wp), &
                    "land_and_ghosts: land cell must stay massless")
         if (allocated(error)) exit checks

         ! Ghost corner (top-left NGHOSTxNGHOST block) must be seeded
         ! IDENTICALLY to a wet interior cell (uniform mode: spatially
         ! invariant conc/h_ice).
         ig = 1
         jg = 1
         call check(error, ig <= NGHOST .and. jg <= NGHOST, &
                    "land_and_ghosts: (ig,jg) must actually be a ghost cell")
         if (allocated(error)) exit checks
         call check(error, ms%wet_mask(ig, jg) > 0.5_wp, &
                    "land_and_ghosts: the default wet_mask must leave ghosts wet "// &
                    "(the seeder trusts wet_mask, not a ghost-band special case)")
         if (allocated(error)) exit checks
         call check(error, ice%part_size(ig, jg, 5) == 1.0_wp, &
                    "land_and_ghosts: ghost cell part_size(5) must match the interior seed")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(ig, jg, 5) == ICE_RHO_ICE*2.0_wp, &
                    "land_and_ghosts: ghost cell m_ice(5) must match the interior seed")

      end block checks
      call teardown(ms, ice, metrics)
   end subroutine test_land_and_ghosts

   subroutine test_ncat1_lumped(error)
      !! Case 7: ncat=1 (legacy lumped mode) -- m_ice per CELL area,
      !! part_size UNTOUCHED (frozen PR 3b/3c contract).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      real(wp), allocatable :: ci(:, :), mis(:, :), mice(:, :)
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, ice, 1, 1)
         call metrics%init(grid)
         par = ice_ic_params_from_config("uniform", 1.0_wp, 2.0_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
         call ice_init_apply(grid, ms, metrics, ice, par)

         call check(error, all(ice%part_size(:, :, 0) == 1.0_wp), &
                    "ncat1_lumped: part_size(:,:,0) must stay 1 (legacy contract)")
         if (allocated(error)) exit checks
         call check(error, all(ice%m_ice(:, :, 1) == ICE_RHO_ICE*2.0_wp), &
                    "ncat1_lumped: m_ice(:,:,1) must be RHO_ICE*h_ice per CELL area")
         if (allocated(error)) exit checks

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         allocate (ci(grid%nx_total, grid%ny_total))
         allocate (mis(grid%nx_total, grid%ny_total))
         allocate (mice(grid%nx_total, grid%ny_total))
         call ice_cell_concentration_impl(metrics%wet_T, ice%part_size, ice%m_ice, &
                                          ice%m_snow, mis, mice, ci, ice%ncat, &
                                          grid%nx_total, grid%ny_total)
         call check(error, ci(ip, jp) == 1.0_wp, "ncat1_lumped: ci must be 1 (binary)")
         if (allocated(error)) exit checks
         call check(error, mice(ip, jp) == ICE_RHO_ICE*2.0_wp, &
                    "ncat1_lumped: mice must equal RHO_ICE*h_ice")

      end block checks
      if (allocated(ci)) deallocate (ci)
      if (allocated(mis)) deallocate (mis)
      if (allocated(mice)) deallocate (mice)
      call teardown(ms, ice, metrics)
   end subroutine test_ncat1_lumped

   subroutine test_latitudes(error)
      !! Case 8: hand-filled `geolatT` ramp; `arctic_edge=60`,
      !! `antarctic_edge=-91` (default => no Antarctic ice). Cells with
      !! lat > 60 must carry conc==1 (binned); every other cell open water.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: par
      integer :: i, j
      real(wp) :: lat
      checks: block

         call setup_state(grid, ms, ice, 5, 2)
         call metrics%init(grid)

         ! Ramp: lat runs from -90 at j=1 to +90 at j=ny_total.
         do j = 1, grid%ny_total
            lat = -90.0_wp + 180.0_wp*real(j - 1, wp)/real(grid%ny_total - 1, wp)
            do i = 1, grid%nx_total
               metrics%geolatT(i, j) = lat
            end do
         end do

         par = ice_ic_params_from_config("latitudes", 0.0_wp, 2.0_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 60.0_wp, -91.0_wp)
         call ice_init_apply(grid, ms, metrics, ice, par)

         do j = 1, grid%ny_total
            lat = -90.0_wp + 180.0_wp*real(j - 1, wp)/real(grid%ny_total - 1, wp)
            do i = 1, grid%nx_total
               if (lat > 60.0_wp) then
                  call check(error, ice%part_size(i, j, 5) == 1.0_wp, &
                             "latitudes: lat>arctic_edge must be fully ice-covered")
               else
                  call check(error, ice%part_size(i, j, 0) == 1.0_wp .and. &
                             all(ice%m_ice(i, j, :) == 0.0_wp), &
                             "latitudes: lat<=arctic_edge (and >antarctic_edge) "// &
                             "must be open water")
               end if
               if (allocated(error)) exit checks
            end do
            if (allocated(error)) exit checks
         end do

      end block checks
      call teardown(ms, ice, metrics)
   end subroutine test_latitudes

   subroutine test_restart_roundtrip_uniform(error)
      !! Case 9: build a 2 m / conc=1 / ncat=5 pack via `ice_init_apply` on
      !! a live `ocean_state_t`, write a restart, destroy, read it back
      !! into a fresh state -- bit-exact (model on
      !! `test_evp_restart_roundtrip`, `test_ocean_ice_evp.F90:1402`).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NGH = 2
      character(len=*), parameter :: FN = "test_ocean_ice_ic_rt.nc"

      type(hgrid_t) :: grid
      type(ocean_state_t) :: a_state, b_state
      type(decomp_t) :: decomp
      type(ice_ic_params_t) :: par
      real(wp), allocatable :: ps_a(:, :, :), mi_a(:, :, :), msn_a(:, :, :)
      real(wp), allocatable :: ei_a(:, :, :, :), es_a(:, :, :, :), si_a(:, :, :, :)
      real(wp) :: t_read
      integer :: step_read

      call grid%init(NX_PHYS, NY_PHYS, NGH, DX, DX)
      call decomp_init(decomp, NX_PHYS, NY_PHYS, 1, 1, 0)
      call cleanup_test_file(FN)

      a_state%multilayer%nz_ml = 1
      a_state%ice%enable = .true.
      a_state%ice%ncat = 5
      a_state%ice%nk_ice = 2
      call a_state%init(grid)
      call metrics_seed_cartesian(a_state, grid)

      checks: block
         par = ice_ic_params_from_config("uniform", 1.0_wp, 2.0_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
         ! Host-side, BEFORE enter_data -- the driver contract.
         call ice_init_apply(grid, a_state%multilayer, a_state%metrics, a_state%ice, par)

         allocate (ps_a, source=a_state%ice%part_size)
         allocate (mi_a, source=a_state%ice%m_ice)
         allocate (msn_a, source=a_state%ice%m_snow)
         allocate (ei_a, source=a_state%ice%enth_ice)
         allocate (es_a, source=a_state%ice%enth_snow)
         allocate (si_a, source=a_state%ice%sal_ice)

         call check(error, any(ps_a(:, :, 5) > 0.0_wp), &
                    "restart_roundtrip: the IC must actually seed something before the write")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(a_state)
         call ocean_state_restart_write(a_state, grid, decomp, FN, 0.0_wp, 0)
         call ocean_state_exit_data(a_state)
         call a_state%destroy()

         b_state%multilayer%nz_ml = 1
         b_state%ice%enable = .true.
         b_state%ice%ncat = 5
         b_state%ice%nk_ice = 2
         call b_state%init(grid)
         call metrics_seed_cartesian(b_state, grid)
         call ocean_state_restart_read(b_state, grid, decomp, FN, t_read, step_read)
         call ocean_state_enter_data(b_state)
         associate (ps => b_state%ice%part_size, mi => b_state%ice%m_ice, &
                    msn => b_state%ice%m_snow, ei => b_state%ice%enth_ice, &
                    es => b_state%ice%enth_snow, si => b_state%ice%sal_ice)
            !$acc update self(ps, mi, msn, ei, es, si)
         end associate

         call check(error, all(b_state%ice%part_size == ps_a), "restart: part_size mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%m_ice == mi_a), "restart: m_ice mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%m_snow == msn_a), "restart: m_snow mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%enth_ice == ei_a), "restart: enth_ice mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%enth_snow == es_a), "restart: enth_snow mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%sal_ice == si_a), "restart: sal_ice mismatch")

         call ocean_state_exit_data(b_state)
      end block checks

      if (allocated(ps_a)) deallocate (ps_a)
      if (allocated(mi_a)) deallocate (mi_a)
      if (allocated(msn_a)) deallocate (msn_a)
      if (allocated(ei_a)) deallocate (ei_a)
      if (allocated(es_a)) deallocate (es_a)
      if (allocated(si_a)) deallocate (si_a)
      call b_state%destroy()
      call cleanup_test_file(FN)
   end subroutine test_restart_roundtrip_uniform

   subroutine metrics_seed_cartesian(state, grid)
      !! Fill the ocean_state's own metrics slot with a Cartesian fill
      !! (mirrors `test_ocean_ice_evp.F90::metrics_seed_cartesian`).
      use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine metrics_seed_cartesian

   subroutine cleanup_test_file(fn)
      character(len=*), intent(in) :: fn
      integer :: u, ios
      character(len=256) :: msg
      logical :: exists
      inquire (file=fn, exist=exists)
      if (exists) then
         open (newunit=u, file=fn, status="old", action="readwrite", iostat=ios, iomsg=msg)
         if (ios == 0) close (u, status="delete")
      end if
   end subroutine cleanup_test_file

   subroutine test_nansen_free_drift_from_ic(error)
      !! Case 10: build a 3 m / conc=1 pack via `ice_init_apply` on a
      !! double-periodic 4x4, then drive the SAME hard-coded 15-digit
      !! Nansen free-drift goldens as `test_nansen_free_drift`
      !! (`test_ocean_ice_evp.F90:112-114`) through `ice_evp_step` (NOT
      !! `ice_evp_dynamics`) -- exercising `ice_cell_concentration_impl`'s
      !! gather from `part_size x m_ice` for the first time from an IC.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NGH = 3
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DXG = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 48
      integer, parameter :: N_TAU = 4
      real(wp), parameter :: TAU_VALUES(N_TAU) = [0.05_wp, 0.1_wp, 0.2_wp, 0.5_wp]
      real(wp), parameter :: U_ANALYTIC(N_TAU) = [0.122403513677564_wp, 0.173104709124932_wp, &
                                                  0.244807027355129_wp, 0.387073896828676_wp]

      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_metrics_t) :: metrics
      type(ice_ic_params_t) :: ic_par
      type(ice_evp_params_t) :: par
      real(wp), allocatable :: f_corner(:, :)
      integer :: nx, ny, itau, n, i, j
      real(wp) :: u_mean, u_min, u_max

      call grid%init(NXP, NYP, NGH, DXG, DXG)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%nz_ml = 1
      call ms%init(grid)
      ice%enable = .true.
      ice%ncat = 5
      ice%nk_ice = 2
      ice%dynamics = .true.
      call ice%init(grid)

      ! Host-side IC, BEFORE any enter_data -- the driver contract. Uniform
      ! 3 m / conc=1 everywhere (incl. ghosts, sec 7.4) => matches the
      ! golden test's hand-filled mis=mice=3*ICE_RHO_ICE, ci=1 exactly.
      ic_par = ice_ic_params_from_config("uniform", 1.0_wp, 3.0_wp, 0.0_wp, &
                                         -4.0_wp, 4.0_wp, 91.0_wp, -91.0_wp)
      call ice_init_apply(grid, ms, metrics, ice, ic_par)

      allocate (f_corner(nx + 1, ny + 1), source=0.0_wp)
      par%p0 = 2.75e4_wp
      par%c0 = 20.0_wp
      par%ec = 2.0_wp
      par%cdw = 3.24e-3_wp
      par%rho_ocean = 1030.0_wp
      par%del_sh_min_scale = 2.0_wp
      par%tdamp = -0.2_wp
      par%evp_sub_steps = 432

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      !$acc enter data copyin(f_corner)

      checks: block
         do itau = 1, N_TAU
            ice%u_ice = 0.0_wp
            ice%v_ice = 0.0_wp
            ice%str_d = 0.0_wp
            ice%str_t = 0.0_wp
            ice%str_s = 0.0_wp
            ice%fxoc = 0.0_wp
            ice%fyoc = 0.0_wp
            ice%tau_a_x = 0.0_wp
            ice%tau_a_y = 0.0_wp
            do j = NGH + 1, NGH + NYP
               do i = NGH + 1, NGH + NXP + 1
                  ice%tau_a_x(i, j) = TAU_VALUES(itau)
               end do
            end do
            associate (ui => ice%u_ice, vi => ice%v_ice, sd => ice%str_d, &
                       st => ice%str_t, ss => ice%str_s, fx => ice%fxoc, &
                       fy => ice%fyoc, tx => ice%tau_a_x, ty => ice%tau_a_y)
               !$acc update device(ui, vi, sd, st, ss, fx, fy, tx, ty)
            end associate

            do n = 1, N_OUTER
               call ice_evp_step(grid, metrics, f_corner, ice, ms, DT_SLOW, par, &
                                 .true., .true.)
            end do

            associate (ui => ice%u_ice, vi => ice%v_ice)
               !$acc update self(ui, vi)
            end associate

            u_mean = sum(ice%u_ice(NGH + 1:NGH + NXP + 1, NGH + 1:NGH + NYP)) &
                     /real((NXP + 1)*NYP, wp)
            u_min = minval(ice%u_ice(NGH + 1:NGH + NXP + 1, NGH + 1:NGH + NYP))
            u_max = maxval(ice%u_ice(NGH + 1:NGH + NXP + 1, NGH + 1:NGH + NYP))

            call check(error, abs(u_mean - U_ANALYTIC(itau)) <= 1.0e-10_wp*U_ANALYTIC(itau), &
                       "nansen_from_ic: |u| mismatch against the golden")
            if (allocated(error)) exit checks
            call check(error, (u_max - u_min) <= 1.0e-9_wp, "nansen_from_ic: face spread")
            if (allocated(error)) exit checks
            call check(error, all(abs(ice%v_ice(NGH + 1:NGH + NXP, NGH + 1:NGH + NYP + 1)) &
                                  <= 1.0e-14_wp), "nansen_from_ic: v /= 0")
            if (allocated(error)) exit checks
         end do
      end block checks

      !$acc exit data delete(f_corner)
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
      call ice%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      if (allocated(f_corner)) deallocate (f_corner)
   end subroutine test_nansen_free_drift_from_ic

end module test_ocean_ice_init
