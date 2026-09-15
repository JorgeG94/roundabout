!! Analytic tests for the multi-category ice thickness distribution (ITD)
!! restore (sea-ice PR 4a): `ice_itd_category_bounds` (`rdb_ice_state`) and
!! `ice_adjust_categories` (`rdb_ice_itd`).
!!
!! Harness mirrors `test_ocean_ice_driver_column.F90`: grid 6x4, NGHOST=2,
!! NZ=3, H_LAYER=10, S_INIT=35; `dt_therm = 1200 s`. Ice seeded directly onto
!! `part_size`/`m_ice`/`m_snow`/`enth_ice`/`enth_snow`/`sal_ice` (never via the
!! frazil path) so each case isolates the ITD-restore algebra. Seed numbers
!! reproduce `tmp_local_artifacts/proto_itd_adjust.py` exactly (see the
!! literal comments at each assertion — they are the prototype's printed
!! numbers, not independently re-derived).
!!
!! **GPU mem:separate discipline** (per PR-3a/3b/3c precedent): every case
!! that calls `ice_adjust_categories` maps state via `ice%enter_data()` before
!! the kernel and pulls the touched arrays host-ward with
!! `!$acc update self` before any host assertion. Directives are inert
!! no-ops on host builds — written unconditionally; a green multicore run
!! proves nothing about device data motion.
!!
!! Cases (SPEC_ice-pr4a-multicat.md §8):
!!   * `itd_bounds` — `ice_itd_category_bounds` exact literals at ncat=1/5/9,
!!     plus the `ice%init`-computed `h_lim`/`mh_lim` match.
!!   * `adjust_conserves` — the prototype's `mixed_grow_melt` case: five
!!     per-cell conservation totals to <=1e-13 rel, every occupied category
!!     in its bin, Σpart==1, and a second call is a bit-exact no-op
!!     (idempotence).
!!   * `grow_promotes_category` — three sub-cases (promote into an occupied
!!     destination, a 3-boundary cascade, a demote) plus the two cleanup
!!     edges (part>0/m=0 cleanup, phantom m>0/part=0 inert).
!!   * `multicat_thermo_conserves` — the full ncat=5 thermo-window chain
!!     (forcing -> basal -> uptake -> thermo-driver -> couplers -> adjust)
!!     WITH a live frazil bank on the cell: the PR-3c driver-column budget
!!     identities, part-weighted, with the frazil salt term proven live
!!     (not a `0.0_wp` placeholder) alongside the net-melt term.
!!   * `ncat1_bit_identical` — an ncat=1 state is an untouched bit-exact
!!     no-op through `ice_adjust_categories` (the short-circuit contract).
!!   * `disabled_bitident` — a never-`init`'d `ocean_sea_ice_t` (is_init
!!     false) is a guarded no-op.
!!   * PR-58 (`hlim` override, PLAN_PR58_ice_hlim_override.md §9) adds:
!!     `itd_bounds` itself is extended (case 1: `hlim_cfg` UNALLOCATED
!!     reproduces the pre-PR-58 default path exactly — the byte-identity
!!     gate). `itd_bounds_override_full` (case 2), `itd_bounds_override_
!!     partial_extrapolates` (case 3 — the extrapolation-resume-point
!!     pin, the single highest-value case), `itd_bounds_override_ncat1`
!!     (case 4), `itd_adjust_respects_override_bins` (case 5 — the knob
!!     changes which category ice occupies), `itd_override_resolves_
!!     thick_pack` (case 6 — the Antarctic motivation: a thick pack that
!!     collapses into one category under the default bins spreads over
!!     >=3 categories under an override).
module test_ocean_ice_itd
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_apply_tracers
   use rdb_ice_state, only: ocean_sea_ice_t, ice_itd_category_bounds
   use rdb_ice_itd, only: ice_adjust_categories
   use rdb_ice_atm_forcing, only: ice_atm_forcing_restoring
   use rdb_ice_basal_flux, only: ice_compute_basal_flux
   use rdb_ice_frazil_uptake, only: ice_frazil_uptake
   use rdb_ice_thermo_driver, only: ice_thermo_driver_step
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux, ice_ocean_heat_flux
   use rdb_ice_column, only: ICE_BULK_SALINITY, ICE_RHO_ICE
   implicit none
   private

   public :: collect_ocean_ice_itd_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: H_LAYER = 10.0_wp
      !! Uniform layer thickness (m).
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Uniform salinity (PSU).
   real(wp), parameter :: T_DEEP = 1.0_wp
      !! Sub-surface IC (degC) — untouched by every case.
   real(wp), parameter :: DT_THERM = 1200.0_wp
      !! Thermo timestep (s).
   integer, parameter :: N_TOTALS = 5
      !! Number of per-cell conservation totals tracked by `cell_totals`/
      !! `check_conserved` (area, ice mass, snow mass, enthalpy, salt).
   real(wp), parameter :: E_BANK = 3.34e6_wp
      !! Analytic frazil bank (J/m² of cell) for `multicat_thermo_conserves`
      !! — same value as `test_ocean_ice_coupling.F90`'s `E_BANK`.

contains

   subroutine collect_ocean_ice_itd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("itd_bounds", test_itd_bounds), &
                  new_unittest("adjust_conserves", test_adjust_conserves), &
                  new_unittest("grow_promotes_category", test_grow_promotes_category), &
                  new_unittest("multicat_thermo_conserves", test_multicat_thermo_conserves), &
                  new_unittest("ncat1_bit_identical", test_ncat1_bit_identical), &
                  new_unittest("disabled_bitident", test_disabled_bitident), &
                  new_unittest("itd_bounds_override_full", test_itd_bounds_override_full), &
                  new_unittest("itd_bounds_override_partial_extrapolates", &
                               test_itd_bounds_override_partial_extrapolates), &
                  new_unittest("itd_bounds_override_ncat1", test_itd_bounds_override_ncat1), &
                  new_unittest("itd_adjust_respects_override_bins", &
                               test_itd_adjust_respects_override_bins), &
                  new_unittest("itd_override_resolves_thick_pack", &
                               test_itd_override_resolves_thick_pack) &
                  ]
   end subroutine collect_ocean_ice_itd_tests

   ! -----------------------------------------------------------------
   ! Shared setup / run / teardown
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, ms, eos, ice, ncat, t_surface, hlim_cfg)
      !! Tiny ocean state: uniform h + S, sub-surface layers at T_DEEP,
      !! surface layer at `t_surface`. Ice slot enabled the driver way
      !! (enable latched before init) with the requested `ncat`.
      !! PR-58: `hlim_cfg`, when present, is latched onto `ice%hlim_cfg`
      !! BEFORE `ice%init` — mirroring `ocean_state_init_from_config`'s
      !! latch-before-init contract (`rdb_ocean_state.F90`).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: ncat
      real(wp), intent(in) :: t_surface
      real(wp), intent(in), optional :: hlim_cfg(:)

      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ice%enable = .true.
      ice%ncat = ncat
      if (present(hlim_cfg)) ice%hlim_cfg = hlim_cfg
      call ice%init(grid)

      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_surface*H_LAYER
   end subroutine setup_state

   subroutine run_adjust(grid, ms, ice)
      !! GPU mem:separate discipline: map, run `ice_adjust_categories`
      !! only, pull the touched arrays host-ward, unmap.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice

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
   end subroutine run_adjust

   subroutine teardown(ms, eos, ice)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      call ice%destroy()
      call eos%destroy()
      call ms%destroy()
   end subroutine teardown

   pure function cell_totals(ice, ip, jp) result(tot)
      !! Per-cell conservation totals: area, ice mass, snow mass,
      !! enthalpy (ice+snow), salt. Mirrors the prototype's `totals()`.
      type(ocean_sea_ice_t), intent(in) :: ice
      integer, intent(in) :: ip, jp
      real(wp) :: tot(N_TOTALS)
      integer :: cat, k
      real(wp) :: mca

      tot = 0.0_wp
      do cat = 0, ice%ncat
         tot(1) = tot(1) + ice%part_size(ip, jp, cat)
      end do
      do cat = 1, ice%ncat
         tot(2) = tot(2) + ice%part_size(ip, jp, cat)*ice%m_ice(ip, jp, cat)
         tot(3) = tot(3) + ice%part_size(ip, jp, cat)*ice%m_snow(ip, jp, cat)
         mca = ice%part_size(ip, jp, cat)*ice%m_ice(ip, jp, cat)
         do k = 1, ice%nk_ice
            tot(4) = tot(4) + (mca/real(ice%nk_ice, wp))*ice%enth_ice(ip, jp, cat, k)
            tot(5) = tot(5) + (mca/real(ice%nk_ice, wp))*ice%sal_ice(ip, jp, cat, k)
         end do
         tot(4) = tot(4) + ice%part_size(ip, jp, cat)*ice%m_snow(ip, jp, cat) &
                  *ice%enth_snow(ip, jp, cat, 1)
      end do
   end function cell_totals

   subroutine check_conserved(error, before, after, label)
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: before(N_TOTALS), after(N_TOTALS)
      character(*), intent(in) :: label
      integer :: k
      real(wp) :: scale, rel
      character(*), parameter :: names(N_TOTALS) = &
                                 ["area  ", "m_ice ", "m_snow", "enth  ", "salt  "]

      do k = 1, N_TOTALS
         scale = max(abs(before(k)), abs(after(k)), 1.0e-30_wp)
         rel = abs(after(k) - before(k))/scale
         call check(error, rel <= 1.0e-13_wp, &
                    label//": "//names(k)//" not conserved to 1e-13 rel")
         if (allocated(error)) return
      end do
   end subroutine check_conserved

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_itd_bounds(error)
      !! Case 1: `ice_itd_category_bounds` exact literals + the
      !! `ice%init`-computed fields match. PR-58: every call below omits
      !! `hlim_vals`/leaves `ice%hlim_cfg` unallocated, so this is also
      !! the default-off BYTE-IDENTITY gate for the PR-58 override — the
      !! `.not. present(hlim_vals)` branch is, textually, the pre-PR-58
      !! code.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: h_lim(:), mh_lim(:)
      integer :: k
      checks: block

         ! ncat=5: h_lim = [1e-10, 0.1, 0.3, 0.7, 1.1, 1.5]; mh_lim = 905*h_lim.
         allocate (h_lim(6), mh_lim(6))
         call ice_itd_category_bounds(5, h_lim, mh_lim)
         call check(error, h_lim(1) == 1.0e-10_wp, "ncat=5 h_lim(1)")
         if (allocated(error)) exit checks
         call check(error, h_lim(2) == 0.1_wp, "ncat=5 h_lim(2)")
         if (allocated(error)) exit checks
         call check(error, h_lim(3) == 0.3_wp, "ncat=5 h_lim(3)")
         if (allocated(error)) exit checks
         call check(error, h_lim(4) == 0.7_wp, "ncat=5 h_lim(4)")
         if (allocated(error)) exit checks
         call check(error, h_lim(5) == 1.1_wp, "ncat=5 h_lim(5)")
         if (allocated(error)) exit checks
         call check(error, h_lim(6) == 1.5_wp, "ncat=5 h_lim(6)")
         if (allocated(error)) exit checks
         do k = 1, 6
            call check(error, mh_lim(k) == ICE_RHO_ICE*h_lim(k), "ncat=5 mh_lim vs ICE_RHO_ICE*h_lim")
            if (allocated(error)) exit checks
         end do
         deallocate (h_lim, mh_lim)
         if (allocated(error)) exit checks

         ! ncat=1: h_lim = [1e-10, 0.1].
         allocate (h_lim(2), mh_lim(2))
         call ice_itd_category_bounds(1, h_lim, mh_lim)
         call check(error, h_lim(1) == 1.0e-10_wp, "ncat=1 h_lim(1)")
         if (allocated(error)) exit checks
         call check(error, h_lim(2) == 0.1_wp, "ncat=1 h_lim(2)")
         if (allocated(error)) exit checks
         deallocate (h_lim, mh_lim)
         if (allocated(error)) exit checks

         ! ncat=9: last two entries extrapolated [3.0, 3.5] within 1e-14.
         allocate (h_lim(10), mh_lim(10))
         call ice_itd_category_bounds(9, h_lim, mh_lim)
         call check(error, h_lim(8) == 2.5_wp, "ncat=9 h_lim(8) (last default)")
         if (allocated(error)) exit checks
         call check(error, abs(h_lim(9) - 3.0_wp) < 1.0e-14_wp, "ncat=9 h_lim(9) extrapolated")
         if (allocated(error)) exit checks
         call check(error, abs(h_lim(10) - 3.5_wp) < 1.0e-14_wp, "ncat=9 h_lim(10) extrapolated")
         if (allocated(error)) exit checks
         deallocate (h_lim, mh_lim)
         if (allocated(error)) exit checks

         ! After ice%init with ncat=5, ice%h_lim/mh_lim match.
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         allocate (h_lim(6), mh_lim(6))
         call ice_itd_category_bounds(5, h_lim, mh_lim)
         do k = 1, 6
            call check(error, ice%h_lim(k) == h_lim(k), "ice%init h_lim mismatch")
            if (allocated(error)) exit checks
            call check(error, ice%mh_lim(k) == mh_lim(k), "ice%init mh_lim mismatch")
            if (allocated(error)) exit checks
         end do

      end block checks
      if (allocated(h_lim)) deallocate (h_lim)
      if (allocated(mh_lim)) deallocate (mh_lim)
      if (ice%is_init) call teardown(ms, eos, ice)
   end subroutine test_itd_bounds

   subroutine test_adjust_conserves(error)
      !! Case 2: the prototype's `mixed_grow_melt` seed on one wet cell.
      !! cat1 h=0.12/part=0.2/snow=3 (promotes to cat2, occupied dest
      !! cat2 h=0.11/part=0.1); cat4 h=0.65/part=0.15/snow=8 (demotes to
      !! cat3, empty dest). part(0) = 1 - Sum. Assert the five per-cell
      !! totals conserved to <=1e-13 rel, every occupied cat in its bin,
      !! Sum(0..ncat) part == 1 to 1e-14, and idempotence (second call
      !! bit-identical).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: part_sum, h
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      real(wp), allocatable :: ps2(:, :, :), mi2(:, :, :), msn2(:, :, :)
      real(wp), allocatable :: ei2(:, :, :, :), es2(:, :, :, :), si2(:, :, :, :)
      integer :: ip, jp, cat
      checks: block

         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2

         ! cat 1: h=0.12, part=0.2, snow=3, heterogeneous enth/sal.
         ice%part_size(ip, jp, 1) = 0.2_wp
         ice%m_ice(ip, jp, 1) = 0.12_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, 1) = 3.0_wp
         ice%enth_ice(ip, jp, 1, 1) = -3.0e5_wp
         ice%enth_ice(ip, jp, 1, 2) = -2.0e5_wp
         ice%sal_ice(ip, jp, 1, 1) = 4.0_wp
         ice%sal_ice(ip, jp, 1, 2) = 5.0_wp
         ice%enth_snow(ip, jp, 1, 1) = -3.0e5_wp

         ! cat 2: h=0.11, part=0.1 (occupied, stays put).
         ice%part_size(ip, jp, 2) = 0.1_wp
         ice%m_ice(ip, jp, 2) = 0.11_wp*ICE_RHO_ICE

         ! cat 4: h=0.65, part=0.15, snow=8 (demotes to cat3, empty dest).
         ice%part_size(ip, jp, 4) = 0.15_wp
         ice%m_ice(ip, jp, 4) = 0.65_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, 4) = 8.0_wp
         ice%enth_ice(ip, jp, 4, 1) = -2.8e5_wp
         ice%enth_ice(ip, jp, 4, 2) = -1.8e5_wp
         ice%enth_snow(ip, jp, 4, 1) = -3.4e5_wp

         part_sum = 0.0_wp
         do cat = 1, ice%ncat
            part_sum = part_sum + ice%part_size(ip, jp, cat)
         end do
         ice%part_size(ip, jp, 0) = max(1.0_wp - part_sum, 0.0_wp)

         before = cell_totals(ice, ip, jp)

         call run_adjust(grid, ms, ice)

         after = cell_totals(ice, ip, jp)
         call check_conserved(error, before, after, "adjust_conserves")
         if (allocated(error)) exit checks

         ! Every occupied category in its bin.
         do cat = 1, ice%ncat
            if (ice%part_size(ip, jp, cat) > 0.0_wp) then
               call check(error, ice%m_ice(ip, jp, cat) >= ice%mh_lim(cat), &
                          "occupied category below its lower bin bound")
               if (allocated(error)) exit checks
               if (cat < ice%ncat) then
                  call check(error, ice%m_ice(ip, jp, cat) < ice%mh_lim(cat + 1), &
                             "occupied category above its upper bin bound")
                  if (allocated(error)) exit checks
               end if
            end if
         end do

         ! Sum(0..ncat) part == 1 to 1e-14.
         part_sum = 0.0_wp
         do cat = 0, ice%ncat
            part_sum = part_sum + ice%part_size(ip, jp, cat)
         end do
         call check(error, abs(part_sum - 1.0_wp) < 1.0e-14_wp, "Sum part must equal 1")
         if (allocated(error)) exit checks

         ! Exact literals from the prototype (mixed_grow_melt): part(1)=0,
         ! part(2)=0.3 (m_ice=905*105.58333...), part(3)=0.15 (receives the
         ! cat4 demote, h=0.65), part(4)=0 (emptied).
         call check(error, ice%part_size(ip, jp, 1) == 0.0_wp, "cat1 emptied")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 2) - 0.3_wp) < 1.0e-14_wp, "cat2 part after merge")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 3) - 0.15_wp) < 1.0e-14_wp, &
                    "cat3 receives the cat4 demote")
         if (allocated(error)) exit checks
         call check(error, ice%part_size(ip, jp, 4) == 0.0_wp, "cat4 emptied by the demote")
         if (allocated(error)) exit checks
         h = ice%m_ice(ip, jp, 3)/ICE_RHO_ICE
         call check(error, abs(h - 0.65_wp) < 1.0e-13_wp, "cat3 h matches the demoted cat4 thickness")
         if (allocated(error)) exit checks

         ! ---- Idempotence: snapshot, call again, assert bit-identical. ----
         allocate (ps2, source=ice%part_size)
         allocate (mi2, source=ice%m_ice)
         allocate (msn2, source=ice%m_snow)
         allocate (ei2, source=ice%enth_ice)
         allocate (es2, source=ice%enth_snow)
         allocate (si2, source=ice%sal_ice)

         call run_adjust(grid, ms, ice)

         call check(error, maxval(abs(ice%part_size - ps2)) == 0.0_wp, &
                    "idempotence: part_size must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_ice - mi2)) == 0.0_wp, &
                    "idempotence: m_ice must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_snow - msn2)) == 0.0_wp, &
                    "idempotence: m_snow must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_ice - ei2)) == 0.0_wp, &
                    "idempotence: enth_ice must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_snow - es2)) == 0.0_wp, &
                    "idempotence: enth_snow must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%sal_ice - si2)) == 0.0_wp, &
                    "idempotence: sal_ice must be bit-identical")

      end block checks
      if (allocated(ps2)) deallocate (ps2)
      if (allocated(mi2)) deallocate (mi2)
      if (allocated(msn2)) deallocate (msn2)
      if (allocated(ei2)) deallocate (ei2)
      if (allocated(es2)) deallocate (es2)
      if (allocated(si2)) deallocate (si2)
      call teardown(ms, eos, ice)
   end subroutine test_adjust_conserves

   subroutine test_grow_promotes_category(error)
      !! Case 3: three sub-cases + two cleanup edges, each its own
      !! fresh state (own wet cell, cat arrays zeroed by ice%init).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      integer :: ip, jp
      real(wp) :: h
      checks: block

         ! ---- (a) promote 1->2, occupied destination. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 1) = 0.3_wp
         ice%m_ice(ip, jp, 1) = 0.15_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 2) = 0.2_wp
         ice%m_ice(ip, jp, 2) = 0.2_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.3_wp - 0.2_wp

         call run_adjust(grid, ms, ice)

         call check(error, ice%part_size(ip, jp, 1) == 0.0_wp, "3a: cat1 emptied")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 2) - 0.5_wp) < 1.0e-14_wp, &
                    "3a: cat2 part == 0.5")
         if (allocated(error)) exit checks
         ! merged m_ice(2) = (0.3*0.15 + 0.2*0.2)*905/0.5.  Tolerance, not exact
         ! ==: the kernel folds ICE_RHO_ICE into each m_ice BEFORE the area-
         ! weighted sum, so GPU FMA contraction of part*m + part*m rounds a few
         ! ULP away from this factored RHS (value ~154; correct to ~1e-13).
         call check(error, abs(ice%m_ice(ip, jp, 2) - (0.3_wp*0.15_wp + 0.2_wp*0.2_wp) &
                               *ICE_RHO_ICE/0.5_wp) <= 1.0e-9_wp, "3a: merged m_ice(2) area-weighted mean")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- (b) cascade 1->4, three boundaries crossed in one call. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 1) = 0.4_wp
         ice%m_ice(ip, jp, 1) = 0.8_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.4_wp

         call run_adjust(grid, ms, ice)

         call check(error, abs(ice%part_size(ip, jp, 4) - 0.4_wp) < 1.0e-14_wp, &
                    "3b: cascade lands part in cat4")
         if (allocated(error)) exit checks
         h = ice%m_ice(ip, jp, 4)/ICE_RHO_ICE
         call check(error, abs(h - 0.8_wp) <= 1.0e-14_wp, "3b: cat4 h matches 0.8")
         if (allocated(error)) exit checks
         call check(error, ice%part_size(ip, jp, 1) == 0.0_wp .and. &
                    ice%part_size(ip, jp, 2) == 0.0_wp .and. &
                    ice%part_size(ip, jp, 3) == 0.0_wp, &
                    "3b: cats 1-3 emptied by the cascade")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- (c) demote 3->2, empty destination. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 3) = 0.25_wp
         ice%m_ice(ip, jp, 3) = 0.25_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.25_wp

         call run_adjust(grid, ms, ice)

         call check(error, ice%part_size(ip, jp, 3) == 0.0_wp, "3c: cat3 emptied")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 2) - 0.25_wp) < 1.0e-14_wp, &
                    "3c: cat2 receives the demote")
         if (allocated(error)) exit checks
         h = ice%m_ice(ip, jp, 2)/ICE_RHO_ICE
         call check(error, abs(h - 0.25_wp) <= 1.0e-14_wp, "3c: cat2 h matches 0.25")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- Cleanup edge: part>0/m=0 => part->0, open water resumed. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 1) = 0.4_wp
         ice%m_ice(ip, jp, 1) = 0.0_wp
         ice%part_size(ip, jp, 3) = 0.3_wp
         ice%m_ice(ip, jp, 3) = 0.5_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.4_wp - 0.3_wp

         call run_adjust(grid, ms, ice)

         call check(error, ice%part_size(ip, jp, 1) == 0.0_wp, &
                    "cleanup: massless category area returns to open water")
         if (allocated(error)) exit checks
         call check(error, abs(ice%part_size(ip, jp, 0) - 0.7_wp) < 1.0e-14_wp, &
                    "cleanup: open water resumed to 0.7")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- Phantom edge: m>0/part=0 => completely inert. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%m_ice(ip, jp, 2) = 2.0_wp*ICE_RHO_ICE
         ! part_size(*,*,2) stays 0 (never set) -> mca = 0, must be inert.
         ice%part_size(ip, jp, 3) = 0.3_wp
         ice%m_ice(ip, jp, 3) = 0.5_wp*ICE_RHO_ICE
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.3_wp

         block
            real(wp), allocatable :: mi_ic(:, :, :)
            allocate (mi_ic, source=ice%m_ice)
            call run_adjust(grid, ms, ice)
            call check(error, ice%part_size(ip, jp, 2) == 0.0_wp, &
                       "phantom: part stays 0 (never touched)")
            if (allocated(error)) then
               deallocate (mi_ic)
               exit checks
            end if
            call check(error, maxval(abs(ice%m_ice(ip, jp, :) - mi_ic(ip, jp, :))) == 0.0_wp, &
                       "phantom: m_ice bit-identical (mca=0 => no transfer fires)")
            deallocate (mi_ic)
         end block

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_grow_promotes_category

   subroutine test_multicat_thermo_conserves(error)
      !! Case 4: the full ncat=5 thermo window (forcing -> basal ->
      !! uptake -> thermo-driver -> couplers -> adjust) on a small wet
      !! grid, mirroring `test_ocean_ice_driver_column`'s harness. Seed a
      !! multi-cat state (occupied cats 1/2/4, nonzero snow on cat4) plus
      !! a frazil bank on the cell and a warm SST (melt side on the
      !! seeded ice). Assert the PR-3c driver-column identities,
      !! part-weighted, WITH the frazil term proven live
      !! (`m_frozen_diag > 0`) — this is the `ncat>1` compose case, the
      !! only gate on the `+=` at `rdb_ice_thermo_driver.F90`'s multicat
      !! reduce (the `ncat==1` twin is `compose_frazil_and_melt_salt` in
      !! `test_ocean_ice_driver_column.F90`). Note the cell state here is
      !! `SST = 5 degC` (drives `fb > 0`) **and** a live bank — a state
      !! the live driver cannot itself reach (the PR-1 frazil clamp only
      !! banks when `SST < T_f`, warming the surface to `T_f` as it
      !! banks), but a legal unit-test state (the kernels gate on
      !! `frazil_heat > 0`, not on SST) that exercises `fb > 0` and the
      !! bank simultaneously — strictly more coverage than the driver can
      !! reach. `compose_frazil_and_melt_salt` covers the
      !! physically-reachable (`fb == 0`) compose state; this is the
      !! adversarial one.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: m_before, m_after, d_m, expected_d_m, tol
      real(wp) :: sum_heat, sum_ocn, sum_ice_m, m_net, q_heat_expected, q_salt_expected
      real(wp) :: part_sum
      real(wp), allocatable :: part_pre_adjust(:)
      integer :: ip, jp, cat
      checks: block

         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         call sf%init(grid)
         call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2

         ! cat1: occupied, cold.
         ice%part_size(ip, jp, 1) = 0.2_wp
         ice%m_ice(ip, jp, 1) = 0.15_wp*ICE_RHO_ICE
         ice%enth_ice(ip, jp, 1, :) = -3.0e5_wp
         ice%sal_ice(ip, jp, 1, :) = ICE_BULK_SALINITY
         ! cat2: occupied, cold.
         ice%part_size(ip, jp, 2) = 0.15_wp
         ice%m_ice(ip, jp, 2) = 0.35_wp*ICE_RHO_ICE
         ice%enth_ice(ip, jp, 2, :) = -2.8e5_wp
         ice%sal_ice(ip, jp, 2, :) = ICE_BULK_SALINITY
         ! cat4: occupied, cold, snow.
         ice%part_size(ip, jp, 4) = 0.1_wp
         ice%m_ice(ip, jp, 4) = 0.9_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, 4) = 5.0_wp
         ice%enth_ice(ip, jp, 4, :) = -2.5e5_wp
         ice%sal_ice(ip, jp, 4, :) = ICE_BULK_SALINITY
         ice%enth_snow(ip, jp, 4, 1) = -3.3e5_wp

         part_sum = 0.2_wp + 0.15_wp + 0.1_wp
         ice%part_size(ip, jp, 0) = 1.0_wp - part_sum

         ! Warm SST => melt side (drives fb) AND a live frazil bank => the
         ! compose case (module docstring): both the frazil salt term and
         ! the net-melt salt term are live in the same window. The bank
         ! spends via the open-water annexation into k_merge (SIS2 ITD
         ! mode, rdb_ice_frazil_uptake) since part(0) > 0 here.
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = 5.0_wp*H_LAYER
         ice%frazil_heat(ip, jp) = E_BANK

         m_before = 0.0_wp
         do cat = 1, ice%ncat
            m_before = m_before + ice%part_size(ip, jp, cat) &
                       *(ice%m_ice(ip, jp, cat) + ice%m_snow(ip, jp, cat))
         end do

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call ice%enter_data()
         call sf%enter_data()

         call ice_atm_forcing_restoring(ice, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp)
         call ice_compute_basal_flux(grid, eos, ms, ice, DT_THERM)
         call ice_frazil_uptake(grid, eos, ms, ice, DT_THERM)
         call ice_thermo_driver_step(grid, eos, ms, ice, DT_THERM)
         call ice_ocean_brine_flux(sf, ice)
         call ice_ocean_heat_flux(sf, ice)

         ! Snapshot part_size AS THE REDUCE KERNEL SAW IT — the driver's
         ! mandated order runs `ice_adjust_categories` AFTER the couplers
         ! (SPEC §6), so the part-weighted identities below must use the
         ! PRE-adjust partition, not whatever ice_adjust_categories leaves
         ! behind (it may fully collapse a melted-out cell to open water).
         allocate (part_pre_adjust(0:ice%ncat))
         associate (ps_pre => ice%part_size)
            !$acc update self(ps_pre)
         end associate
         part_pre_adjust(:) = ice%part_size(ip, jp, :)

         call ice_adjust_categories(grid, ms, ice)

         associate (ps => ice%part_size, mi => ice%m_ice, msn => ice%m_snow, &
                    ei => ice%enth_ice, es => ice%enth_snow, si => ice%sal_ice, &
                    md => ice%m_melt_diag, mf => ice%m_frozen_diag, &
                    hf => ice%heat_flux_diag, sd => ice%salt_flux_diag, &
                    fb => ice%fb, fps => ice%fb_part_sum, &
                    h2o => ice%h2o_ocn_to_ice, h2i => ice%h2o_ice_to_ocn, &
                    hto => ice%heat_to_ocn, qh => sf%Q_heat, qs => sf%Q_salt, &
                    ssurf => ice%ssurf_seam)
            !$acc update self(ps, mi, msn, ei, es, si, md, mf, hf, sd, fb, fps, h2o, h2i, hto, qh, qs, ssurf)
         end associate

         call sf%exit_data()
         call ice%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         m_after = 0.0_wp
         do cat = 1, ice%ncat
            m_after = m_after + ice%part_size(ip, jp, cat) &
                      *(ice%m_ice(ip, jp, cat) + ice%m_snow(ip, jp, cat))
         end do
         d_m = m_after - m_before

         ! Delta(Sum_c part*(m_ice+m_snow)) == m_frozen_diag - m_melt_diag.
         expected_d_m = ice%m_frozen_diag(ip, jp) - ice%m_melt_diag(ip, jp)
         tol = 1.0e-12_wp*max(abs(expected_d_m), 1.0e-6_wp)
         call check(error, abs(d_m - expected_d_m) <= tol, &
                    "Delta(part-weighted ice+snow mass) must equal m_frozen-m_melt")
         if (allocated(error)) exit checks

         ! Q_heat == Sum_c part*heat_to_ocn/dt - fb_part_sum*fb (recomputed
         ! in the test — the PR-3c identity, part-weighted). Uses the
         ! PRE-adjust partition (the reduce kernel ran before adjust).
         sum_heat = 0.0_wp
         do cat = 1, ice%ncat
            sum_heat = sum_heat + part_pre_adjust(cat)*ice%heat_to_ocn(ip, jp, cat)
         end do
         q_heat_expected = sum_heat/DT_THERM - ice%fb_part_sum(ip, jp)*ice%fb(ip, jp)
         tol = 1.0e-12_wp*max(abs(q_heat_expected), 1.0e-6_wp)
         call check(error, abs(sf%Q_heat(ip, jp) - q_heat_expected) <= tol, &
                    "Q_heat must equal part-weighted heat_to_ocn/dt - fb_part_sum*fb")
         if (allocated(error)) exit checks

         ! Q_salt equals the uptake + net-melt formula with part weights
         ! (PRE-adjust partition, same reasoning as Q_heat above).
         sum_ocn = 0.0_wp
         sum_ice_m = 0.0_wp
         do cat = 1, ice%ncat
            sum_ocn = sum_ocn + part_pre_adjust(cat)*ice%h2o_ocn_to_ice(ip, jp, cat)
            sum_ice_m = sum_ice_m + part_pre_adjust(cat)*ice%h2o_ice_to_ocn(ip, jp, cat)
         end do
         m_net = sum_ocn - sum_ice_m
         ! salt_flux_diag was zeroed+written by the frazil uptake (the LIVE
         ! bank term, part-weighted into m_frozen_diag by the multicat
         ! spend kernel), then the thermo-driver ADDED its net-melt term
         ! on top (the compose contract, rdb_ice_thermo_driver.F90:265).
         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "m_frozen_diag must be strictly positive (frazil term live, not a 0.0 no-op)")
         if (allocated(error)) exit checks
         q_salt_expected = ice%m_frozen_diag(ip, jp) &
                           *(ice%ssurf_seam(ip, jp) - ICE_BULK_SALINITY)/DT_THERM &
                           + m_net*(ice%ssurf_seam(ip, jp) - ICE_BULK_SALINITY)/DT_THERM
         tol = 1.0e-12_wp*max(abs(q_salt_expected), 1.0e-6_wp)
         call check(error, abs(ice%salt_flux_diag(ip, jp) - q_salt_expected) <= tol, &
                    "salt_flux_diag must equal frazil + part-weighted net-melt salt formula")
         if (allocated(error)) exit checks

         ! Sum part == 1.
         part_sum = 0.0_wp
         do cat = 0, ice%ncat
            part_sum = part_sum + ice%part_size(ip, jp, cat)
         end do
         call check(error, abs(part_sum - 1.0_wp) < 1.0e-13_wp, "Sum part must equal 1 after the window")
         if (allocated(error)) exit checks

         ! Occupied categories in their bins after adjust.
         do cat = 1, ice%ncat
            if (ice%part_size(ip, jp, cat) > 0.0_wp) then
               call check(error, ice%m_ice(ip, jp, cat) >= ice%mh_lim(cat), &
                          "post-window occupied category below its lower bin bound")
               if (allocated(error)) exit checks
               if (cat < ice%ncat) then
                  call check(error, ice%m_ice(ip, jp, cat) < ice%mh_lim(cat + 1), &
                             "post-window occupied category above its upper bin bound")
                  if (allocated(error)) exit checks
               end if
            end if
         end do

      end block checks
      call sf%destroy()
      call teardown(ms, eos, ice)
   end subroutine test_multicat_thermo_conserves

   subroutine test_ncat1_bit_identical(error)
      !! Case 5: an ncat=1 state built the PR-3c way (part(0)=1
      !! untouched, m_ice seeded per cell). Snapshot every ice array +
      !! part_size, call `ice_adjust_categories`, assert bit-exact no-op
      !! (the short-circuit contract).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: ps_ic(:, :, :), mi_ic(:, :, :), msn_ic(:, :, :)
      real(wp), allocatable :: ei_ic(:, :, :, :), es_ic(:, :, :, :), si_ic(:, :, :, :)
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, 1, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%m_ice(ip, jp, 1) = 0.4_wp*ICE_RHO_ICE
         ice%enth_ice(ip, jp, 1, :) = -2.0e5_wp
         ice%sal_ice(ip, jp, 1, :) = ICE_BULK_SALINITY

         allocate (ps_ic, source=ice%part_size)
         allocate (mi_ic, source=ice%m_ice)
         allocate (msn_ic, source=ice%m_snow)
         allocate (ei_ic, source=ice%enth_ice)
         allocate (es_ic, source=ice%enth_snow)
         allocate (si_ic, source=ice%sal_ice)

         call run_adjust(grid, ms, ice)

         call check(error, maxval(abs(ice%part_size - ps_ic)) == 0.0_wp, &
                    "ncat=1: part_size must be bit-identical (short-circuit)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_ice - mi_ic)) == 0.0_wp, &
                    "ncat=1: m_ice must be bit-identical (short-circuit)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_snow - msn_ic)) == 0.0_wp, &
                    "ncat=1: m_snow must be bit-identical (short-circuit)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_ice - ei_ic)) == 0.0_wp, &
                    "ncat=1: enth_ice must be bit-identical (short-circuit)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_snow - es_ic)) == 0.0_wp, &
                    "ncat=1: enth_snow must be bit-identical (short-circuit)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%sal_ice - si_ic)) == 0.0_wp, &
                    "ncat=1: sal_ice must be bit-identical (short-circuit)")

      end block checks
      if (allocated(ps_ic)) deallocate (ps_ic)
      if (allocated(mi_ic)) deallocate (mi_ic)
      if (allocated(msn_ic)) deallocate (msn_ic)
      if (allocated(ei_ic)) deallocate (ei_ic)
      if (allocated(es_ic)) deallocate (es_ic)
      if (allocated(si_ic)) deallocate (si_ic)
      call teardown(ms, eos, ice)
   end subroutine test_ncat1_bit_identical

   subroutine test_disabled_bitident(error)
      !! Case 6: `ice%is_init == .false.` (never `init`'d) — assert
      !! `ice_adjust_categories` is a guarded no-op (no touch/crash).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      checks: block

         call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         ! ice never init'd: enable stays .false., is_init stays .false.

         call check(error,.not. ice%is_init, "is_init must be false before the no-op call")
         if (allocated(error)) exit checks

         call ice_adjust_categories(grid, ms, ice)

         call check(error,.not. ice%is_init, "is_init must remain false (guarded no-op)")

      end block checks
      call eos%destroy()
      call ms%destroy()
   end subroutine test_disabled_bitident

   ! -----------------------------------------------------------------
   ! PR-58: ITD category-bound override (`hlim`)
   ! -----------------------------------------------------------------

   subroutine test_itd_bounds_override_full(error)
      !! Case 2 (PLAN_PR58 §9): a fully-specified `hlim_vals` (ncat+1
      !! entries) reproduces itself exactly in `h_lim`, and
      !! `mh_lim = ICE_RHO_ICE*h_lim` exactly. No extrapolation is in
      !! play — a failure here is a pure plumbing failure.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_lim(:), mh_lim(:)
      real(wp), parameter :: hv(6) = [0.05_wp, 0.5_wp, 1.2_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      integer :: k
      checks: block

         allocate (h_lim(6), mh_lim(6))
         call ice_itd_category_bounds(5, h_lim, mh_lim, hlim_vals=hv)
         do k = 1, 6
            call check(error, h_lim(k) == hv(k), "override_full: h_lim mismatch")
            if (allocated(error)) exit checks
            call check(error, mh_lim(k) == ICE_RHO_ICE*h_lim(k), &
                       "override_full: mh_lim vs ICE_RHO_ICE*h_lim")
            if (allocated(error)) exit checks
         end do

      end block checks
      if (allocated(h_lim)) deallocate (h_lim)
      if (allocated(mh_lim)) deallocate (mh_lim)
   end subroutine test_itd_bounds_override_full

   subroutine test_itd_bounds_override_partial_extrapolates(error)
      !! Case 3 (PLAN_PR58 §9) — THE single highest-value case. A
      !! partial `hlim_vals` (3 entries at ncat=5) must extrapolate by
      !! continuing the USER'S constant width from `n+1 = 4`, NOT from a
      !! fixed index 9 (`N_HLIM_DFLT + 1`) and NOT by falling back into
      !! `HLIM_DFLT_TABLE`. Both wrong readings are natural implementation
      !! errors and both produce a plausible, monotone, silently-wrong
      !! ladder that every other case in this file would pass.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_lim(:), mh_lim(:)
      real(wp), parameter :: hv(3) = [0.05_wp, 0.5_wp, 1.2_wp]
         !! Constant width Delta = 1.2 - 0.5 = 0.7, continued upward.
      checks: block

         allocate (h_lim(6), mh_lim(6))
         call ice_itd_category_bounds(5, h_lim, mh_lim, hlim_vals=hv)
         call check(error, h_lim(1) == hv(1), "override_partial: h_lim(1)")
         if (allocated(error)) exit checks
         call check(error, h_lim(2) == hv(2), "override_partial: h_lim(2)")
         if (allocated(error)) exit checks
         call check(error, h_lim(3) == hv(3), "override_partial: h_lim(3)")
         if (allocated(error)) exit checks
         call check(error, abs(h_lim(4) - 1.9_wp) < 1.0e-14_wp, &
                    "override_partial: h_lim(4) must continue the USER width (0.7), not the default table")
         if (allocated(error)) exit checks
         call check(error, abs(h_lim(5) - 2.6_wp) < 1.0e-14_wp, "override_partial: h_lim(5)")
         if (allocated(error)) exit checks
         call check(error, abs(h_lim(6) - 3.3_wp) < 1.0e-14_wp, "override_partial: h_lim(6)")
         if (allocated(error)) exit checks
         call check(error, abs(mh_lim(4) - ICE_RHO_ICE*h_lim(4)) < 1.0e-10_wp, &
                    "override_partial: mh_lim(4) vs ICE_RHO_ICE*h_lim(4)")

      end block checks
      if (allocated(h_lim)) deallocate (h_lim)
      if (allocated(mh_lim)) deallocate (mh_lim)
   end subroutine test_itd_bounds_override_partial_extrapolates

   subroutine test_itd_bounds_override_ncat1(error)
      !! Case 4 (PLAN_PR58 §9): the ncat=1 legacy lumped mode still
      !! accepts and honours the minimum 2-entry list. Guards an
      !! off-by-one in the `min(ncat+1, size(hlim_vals))` clamp that the
      !! ncat=5 cases cannot see.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_lim(:), mh_lim(:)
      real(wp), parameter :: hv(2) = [0.05_wp, 0.6_wp]
      checks: block

         allocate (h_lim(2), mh_lim(2))
         call ice_itd_category_bounds(1, h_lim, mh_lim, hlim_vals=hv)
         call check(error, h_lim(1) == hv(1), "override_ncat1: h_lim(1)")
         if (allocated(error)) exit checks
         call check(error, h_lim(2) == hv(2), "override_ncat1: h_lim(2)")

      end block checks
      if (allocated(h_lim)) deallocate (h_lim)
      if (allocated(mh_lim)) deallocate (mh_lim)
   end subroutine test_itd_bounds_override_ncat1

   subroutine test_itd_adjust_respects_override_bins(error)
      !! Case 5 (PLAN_PR58 §9) — THE PHYSICS GATE. One wet cell, ncat=5,
      !! a single occupied category holding ice at h=1.8m. (a) Under the
      !! default bins (h_lim(5)=1.1, unbounded above) the seed cascades
      !! all the way to cat 5. (b) Under an override
      !! (hlim=[1e-10,0.5,1.2,2.0,3.0,4.0]) the SAME seed lands in cat 3
      !! (1.2 <= 1.8 < 2.0). (c) In BOTH runs the five per-cell
      !! conservation totals hold to <=1e-13 rel — no bounds test on
      !! h_lim alone can show that the knob reaches the re-binning
      !! kernel without breaking `ice_adjust_categories`' conservation.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: hv(6) = [1.0e-10_wp, 0.5_wp, 1.2_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      integer :: ip, jp, cat
      checks: block

         ! ---- (a) default bins: h=1.8 cascades to cat 5. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 1) = 0.3_wp
         ice%m_ice(ip, jp, 1) = 1.8_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, 1) = 4.0_wp
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.3_wp

         before = cell_totals(ice, ip, jp)
         call run_adjust(grid, ms, ice)
         after = cell_totals(ice, ip, jp)
         call check_conserved(error, before, after, "override_bins(default)")
         if (allocated(error)) exit checks

         call check(error, ice%part_size(ip, jp, 5) > 0.0_wp, &
                    "default bins: h=1.8 must land in cat5 (top, unbounded)")
         if (allocated(error)) exit checks
         do cat = 1, 4
            call check(error, ice%part_size(ip, jp, cat) == 0.0_wp, &
                       "default bins: cats 1-4 must be emptied by the cascade to cat5")
            if (allocated(error)) exit checks
         end do
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- (b) override bins: the SAME seed lands in cat 3. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp, hlim_cfg=hv)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%part_size(ip, jp, 1) = 0.3_wp
         ice%m_ice(ip, jp, 1) = 1.8_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, 1) = 4.0_wp
         ice%part_size(ip, jp, 0) = 1.0_wp - 0.3_wp

         before = cell_totals(ice, ip, jp)
         call run_adjust(grid, ms, ice)
         after = cell_totals(ice, ip, jp)
         call check_conserved(error, before, after, "override_bins(override)")
         if (allocated(error)) exit checks

         call check(error, ice%part_size(ip, jp, 3) > 0.0_wp, &
                    "override bins: h=1.8 must land in cat3 (1.2 <= 1.8 < 2.0)")
         if (allocated(error)) exit checks
         call check(error, ice%part_size(ip, jp, 4) == 0.0_wp .and. &
                    ice%part_size(ip, jp, 5) == 0.0_wp, &
                    "override bins: cats 4-5 must stay empty (the override bounds the pack lower)")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_itd_adjust_respects_override_bins

   subroutine test_itd_override_resolves_thick_pack(error)
      !! Case 6 (PLAN_PR58 §9) — THE ANTARCTIC MOTIVATION. One wet cell,
      !! ncat=5, four equal-area category populations at
      !! h=1.2/1.6/2.0/2.4 m (RESUME §10's validated polar freeze-up run
      !! landed a 2.49 m mean pack against a 1.1 m default top edge).
      !! (a) Default bins: the ENTIRE ice area collapses into cat5 —
      !! `part_size(5) == Sum(part_size(1:5))`, `part_size(1:4) == 0`
      !! — reproducing the pathology this PR exists to fix, as a test
      !! that can never be re-forgotten. (b) An override placing bin
      !! edges strictly between each pair of seed thicknesses
      !! (hlim=[1e-10,1.3,1.7,2.1,2.5,3.0]) keeps every seed in its OWN
      !! starting category (no promotion fires — each h_k is already
      !! inside its own bin) — the area spreads over >= 3 distinct
      !! categories (here, all 4).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: hv(6) = [1.0e-10_wp, 1.3_wp, 1.7_wp, 2.1_wp, 2.5_wp, 3.0_wp]
      real(wp), parameter :: seed_h(4) = [1.2_wp, 1.6_wp, 2.0_wp, 2.4_wp]
      real(wp), parameter :: seed_part = 0.1_wp
      real(wp) :: part_sum
      integer :: ip, jp, cat
      checks: block

         ! ---- (a) default bins: total collapse into cat5. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         do cat = 1, 4
            ice%part_size(ip, jp, cat) = seed_part
            ice%m_ice(ip, jp, cat) = seed_h(cat)*ICE_RHO_ICE
         end do
         ice%part_size(ip, jp, 0) = 1.0_wp - 4.0_wp*seed_part

         call run_adjust(grid, ms, ice)

         part_sum = 0.0_wp
         do cat = 1, 5
            part_sum = part_sum + ice%part_size(ip, jp, cat)
         end do
         call check(error, abs(ice%part_size(ip, jp, 5) - part_sum) < 1.0e-14_wp, &
                    "default bins: cat5 must hold the ENTIRE ice area")
         if (allocated(error)) exit checks
         do cat = 1, 4
            call check(error, ice%part_size(ip, jp, cat) == 0.0_wp, &
                       "default bins: cats 1-4 must be emptied (the pathology this PR fixes)")
            if (allocated(error)) exit checks
         end do
         call teardown(ms, eos, ice)
         if (allocated(error)) exit checks

         ! ---- (b) override bins: the pack spreads over >= 3 categories. ----
         call setup_state(grid, ms, eos, ice, 5, 0.0_wp, hlim_cfg=hv)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         do cat = 1, 4
            ice%part_size(ip, jp, cat) = seed_part
            ice%m_ice(ip, jp, cat) = seed_h(cat)*ICE_RHO_ICE
         end do
         ice%part_size(ip, jp, 0) = 1.0_wp - 4.0_wp*seed_part

         call run_adjust(grid, ms, ice)

         call check(error, count(ice%part_size(ip, jp, 1:5) > 0.0_wp) >= 3, &
                    "override bins: the thick pack must resolve into >= 3 categories "// &
                    "(the Antarctic gate — the reason this PR exists)")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_itd_override_resolves_thick_pack

end module test_ocean_ice_itd
