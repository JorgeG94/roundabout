!! Unit tests for the configure-time ocean stability audit
!! (`rdb_ocean_stability_audit.F90`). Motivated by a real failure
!! (`tmp_local_artifacts/global_run/FINDINGS.md`): a global tripolar
!! aquaplanet NaN'd at outer step 7 with NO configure-time diagnostic
!! naming the cause (a 13x viscous-CFL violation). "Things that could be
!! unit tests should be" — these are the tests that should have existed
!! before that debugging session.
!!
!! Two layers:
!!   1. Pure-function tests on the checked-in formulae (no I/O, no
!!      config/grid scaffolding) — cheap, exhaustive edge cases.
!!   2. End-to-end `ocean_stability_audit` tests against a real
!!      `config_t` + `ocean_metrics_t` + `hgrid_t`, including the EXACT
!!      regression numbers from the real failure (nu_h=2e4, dt=900,
!!      dx=3.3km) and a Cartesian control that must NOT trip.
module test_ocean_stability_audit
   use rdb_constants, only: wp, VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_SIGMA, VCOORD_Z_FIXED, VCOORD_RHO
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, make_spherical_metrics, &
                                 destroy_cartesian_metrics
   use rdb_ocean_stability_audit, only: ocean_stability_audit, &
                                        ocean_viscous_cfl_number, &
                                        ocean_viscous_cfl_max_nu_h, &
                                        ocean_viscous_cfl_max_dt, &
                                        ocean_diffusive_number, &
                                        ocean_munk_delta_m, &
                                        ocean_munk_required_nu_h, &
                                        ocean_viscous_cfl_limit, &
                                        ocean_sigma_stiffness, &
                                        ocean_sigma_stiffness_worst, &
                                        ocean_sigma_stiffness_limit, &
                                        ocean_vcoord_is_terrain_following
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_stability_audit_tests

contains

   subroutine collect_ocean_stability_audit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("viscous_cfl_number_formula", test_viscous_cfl_number), &
                  new_unittest("viscous_cfl_fix_formulae", test_viscous_cfl_fix_formulae), &
                  new_unittest("diffusive_number_formula", test_diffusive_number), &
                  new_unittest("munk_delta_m_formula", test_munk_delta_m), &
                  new_unittest("munk_required_nu_h_formula", test_munk_required_nu_h), &
                  new_unittest("audit_viscous_cfl_real_failure_errors", &
                               test_audit_viscous_cfl_real_failure_errors), &
                  new_unittest("audit_viscous_cfl_bound_kh_downgrades", &
                               test_audit_viscous_cfl_bound_kh_downgrades), &
                  new_unittest("audit_viscous_cfl_safe_passes", test_audit_viscous_cfl_safe_passes), &
                  new_unittest("audit_viscous_cfl_ignores_land_cell", &
                               test_audit_viscous_cfl_ignores_land_cell), &
                  new_unittest("audit_kappa_h_spherical_uses_real_metric", &
                               test_audit_kappa_h_spherical_uses_real_metric), &
                  new_unittest("audit_ah_max_below_nu_h_warns_only", &
                               test_audit_ah_max_below_nu_h_warns_only), &
                  new_unittest("sigma_stiffness_formula", test_sigma_stiffness_formula), &
                  new_unittest("sigma_stiffness_vcoord_gate", test_sigma_stiffness_vcoord_gate), &
                  new_unittest("sigma_stiffness_worst_finds_isomip_sidewall", &
                               test_sigma_stiffness_worst_finds_isomip_sidewall), &
                  new_unittest("sigma_stiffness_land_is_not_a_stiff_face", &
                               test_sigma_stiffness_land_is_not_a_stiff_face), &
                  new_unittest("sigma_stiffness_uniform_column_is_zero", &
                               test_sigma_stiffness_uniform_column_is_zero) &
                  ]
   end subroutine collect_ocean_stability_audit_tests

   ! -----------------------------------------------------------------
   ! Layer 1: pure-function formula tests
   ! -----------------------------------------------------------------

   subroutine test_viscous_cfl_number(error)
      !! Pins the exact regression numbers from the real failure
      !! (FINDINGS.md section 6): nu_h=2e4, dt=900s, dx=3300m ->
      !! viscous CFL ~1.65, ~13x the 0.125 limit.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: cfl

      cfl = ocean_viscous_cfl_number(2.0e4_wp, 900.0_wp, 3300.0_wp)
      ! 2e4*900/3300^2 = 1.8e7/1.089e7 = 1.652892...
      call check(error, abs(cfl - 1.652892561983471_wp) < 1.0e-9_wp, &
                 "FINDINGS.md regression number: viscous CFL ~1.65, ~13x the 0.125 limit")
      if (allocated(error)) return
      call check(error, cfl/ocean_viscous_cfl_limit() > 13.0_wp, &
                 "must be >13x the limit, matching the '13x over' finding")
   end subroutine test_viscous_cfl_number

   subroutine test_viscous_cfl_fix_formulae(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: nu_h_max, dt_max

      ! Inverse relationships: max_nu_h(dt, dx, limit) * dt / dx^2 == limit.
      nu_h_max = ocean_viscous_cfl_max_nu_h(900.0_wp, 3300.0_wp, 0.125_wp)
      call check(error, abs(ocean_viscous_cfl_number(nu_h_max, 900.0_wp, 3300.0_wp) - 0.125_wp) < 1.0e-9_wp, &
                 "max_nu_h must saturate the limit exactly")
      if (allocated(error)) return

      dt_max = ocean_viscous_cfl_max_dt(2.0e4_wp, 3300.0_wp, 0.125_wp)
      call check(error, abs(ocean_viscous_cfl_number(2.0e4_wp, dt_max, 3300.0_wp) - 0.125_wp) < 1.0e-9_wp, &
                 "max_dt must saturate the limit exactly")
      if (allocated(error)) return

      ! Degenerate dx_min <= 0 -> no constraint expressible (returns 0, never
      ! spuriously trips a > limit test upstream).
      call check(error, ocean_viscous_cfl_number(2.0e4_wp, 900.0_wp, 0.0_wp) == 0.0_wp, &
                 "dx_min<=0 -> cfl=0 (no constraint)")
   end subroutine test_viscous_cfl_fix_formulae

   subroutine test_diffusive_number(error)
      !! Two-axis form: kappa_h*dt*(2/dx^2) at the same dx on both axes.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dnum

      dnum = ocean_diffusive_number(25.0_wp, 1200.0_wp, 2000.0_wp)
      ! 25 * 1200 * 2 / 2000^2 = 60000/4000000 = 0.015
      call check(error, abs(dnum - 0.015_wp) < 1.0e-12_wp, "diffusive number formula")
      if (allocated(error)) return
      call check(error, ocean_diffusive_number(25.0_wp, 1200.0_wp, 0.0_wp) == 0.0_wp, &
                 "dx_min<=0 -> 0")
   end subroutine test_diffusive_number

   subroutine test_munk_delta_m(error)
      !! delta_M = (nu_h/beta)^(1/3); matches the acc_channel.nml
      !! hand-derived numbers (nu_h=22000, beta~1.62e-11 -> ~119.6 km,
      !! well above the acc_channel dy=55.5km*2=111km 2-cell floor).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: delta_m

      delta_m = ocean_munk_delta_m(22000.0_wp, 1.62e-11_wp)
      ! (22000/1.62e-11)^(1/3) = 110739.499... m.
      call check(error, abs(delta_m - 110739.4992680662_wp) < 1.0e-6_wp, &
                 "acc_channel.nml Munk width ~110.7 km")
      if (allocated(error)) return
      ! f-plane (beta<=0): no Munk boundary layer -> no constraint (huge).
      call check(error, ocean_munk_delta_m(22000.0_wp, 0.0_wp) == huge(1.0_wp), &
                 "beta<=0 -> huge (no constraint)")
      if (allocated(error)) return
      call check(error, ocean_munk_delta_m(0.0_wp, 1.62e-11_wp) == huge(1.0_wp), &
                 "nu_h<=0 -> huge (no constraint)")
   end subroutine test_munk_delta_m

   subroutine test_munk_required_nu_h(error)
      !! acc_channel.nml's own hand-derived number: nu_h >= beta*(2*dy)^3
      !! at beta=1.62e-11, dy=55.5km -> ~2.2e4 (the file's own nu_h).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: nu_h_req

      nu_h_req = ocean_munk_required_nu_h(1.62e-11_wp, 55500.0_wp, 2.0_wp)
      call check(error, nu_h_req > 2.0e4_wp .and. nu_h_req < 2.4e4_wp, &
                 "required nu_h matches acc_channel.nml's own derivation (~2.2e4)")
   end subroutine test_munk_required_nu_h

   ! -----------------------------------------------------------------
   ! Layer 2: end-to-end `ocean_stability_audit` tests
   ! -----------------------------------------------------------------

   subroutine base_config(cfg, dt, nu_h, ah_max, kappa_h, bound_kh, stress_tensor)
      !! Minimal config with just the knobs the audit reads.
      type(config_t), intent(inout) :: cfg
      real(wp), intent(in) :: dt, nu_h, ah_max, kappa_h
      logical, intent(in) :: bound_kh, stress_tensor
      cfg%dt_fixed = dt
      cfg%ocean%hvisc%nu_h = nu_h
      cfg%ocean%hvisc%ah_max = ah_max
      cfg%ocean%hvisc%bound_kh = bound_kh
      cfg%ocean%hvisc%stress_tensor = stress_tensor
      cfg%ocean%hdiff%kappa_h = kappa_h
      cfg%ocean%vmix%dt_therm_ratio = 1
      cfg%ocean%topo%coriolis_beta = 0.0_wp   ! f-plane: Munk check inert
      cfg%ocean%grid%grid_config = "cartesian"
   end subroutine base_config

   subroutine test_audit_viscous_cfl_real_failure_errors(error)
      !! THE regression case: nu_h=2e4, dt=900s, dx=3.3km, bound_kh OFF
      !! (FINDINGS.md section 6, the actual global-aquaplanet failure).
      !! Must be a hard ERROR (OCEAN_STATUS_ERR_SETUP), not a silent
      !! pass-through to a step-7 NaN.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         call grid%init(8, 8, 2, 3300.0_wp, 3300.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call base_config(cfg, 900.0_wp, 2.0e4_wp, 1.0e5_wp, 0.0_wp, &
                          bound_kh=.false., stress_tensor=.false.)

         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)

         call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                    "unprotected nu_h=2e4/dt=900/dx=3.3km must ERROR (the real failure)")
      end block checks
   end subroutine test_audit_viscous_cfl_real_failure_errors

   subroutine test_audit_viscous_cfl_bound_kh_downgrades(error)
      !! The SAME configuration as the real failure, but with bound_kh
      !! enabled (FINDINGS.md's actual fix) — must NOT error: the
      !! runtime per-cell clamp makes the raw nu_h figure informational.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         call grid%init(8, 8, 2, 3300.0_wp, 3300.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call base_config(cfg, 900.0_wp, 2.0e4_wp, 1.0e5_wp, 0.0_wp, &
                          bound_kh=.true., stress_tensor=.false.)

         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)

         call check(error, ierr == OCEAN_STATUS_OK, &
                    "bound_kh=.true. must downgrade the same config to a pass "// &
                    "(the FINDINGS.md fix)")
      end block checks
   end subroutine test_audit_viscous_cfl_bound_kh_downgrades

   subroutine test_audit_viscous_cfl_safe_passes(error)
      !! A config comfortably inside the bound must not trip anything —
      !! the negative control every check above needs.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         call grid%init(8, 8, 2, 3300.0_wp, 3300.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ! nu_h*dt/dx^2 = 1000*900/3300^2 = 0.0826 < 0.125.
         call base_config(cfg, 900.0_wp, 1000.0_wp, 1.0e5_wp, 0.0_wp, &
                          bound_kh=.false., stress_tensor=.false.)

         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)

         call check(error, ierr == OCEAN_STATUS_OK, &
                    "a config comfortably inside the viscous-CFL bound must pass")
      end block checks
   end subroutine test_audit_viscous_cfl_safe_passes

   subroutine test_audit_viscous_cfl_ignores_land_cell(error)
      !! The 1-degree tripolar abort in miniature: a 362 m cell among 33 km
      !! cells.  As a LAND cell it must not trip the viscous-CFL check (no
      !! operator acts there); the SAME cell made wet must.  Both configs
      !! are otherwise identical, so the pair isolates the wet mask.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      real(wp), allocatable :: wm(:, :)
      integer :: ierr, il, jl

      checks: block
         call grid%init(8, 8, 2, 3.3e4_wp, 3.3e4_wp)
         il = grid%nghost + 4
         jl = grid%nghost + 4
         allocate (wm(grid%nx_total, grid%ny_total), source=1.0_wp)
         wm(il, jl) = 0.0_wp
         ! nu_h*dt/dx^2: 33 km -> 1000*900/3.3e4^2 = 8.3e-4 (fine);
         ! 362 m -> 6.9 (55x over 0.125).
         call base_config(cfg, 900.0_wp, 1000.0_wp, 1.0e5_wp, 0.0_wp, &
                          bound_kh=.false., stress_tensor=.false.)

         call make_cartesian_metrics(metrics, grid, wet_mask=wm)
         ! Host-only edit: the audit is a host-side configure scan.
         metrics%dxT(il, jl) = 362.0_wp
         metrics%dyT(il, jl) = 362.0_wp
         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "a tiny LAND cell must not trip the viscous-CFL check")
         if (allocated(error)) exit checks

         call make_cartesian_metrics(metrics, grid)
         metrics%dxT(il, jl) = 362.0_wp
         metrics%dyT(il, jl) = 362.0_wp
         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)
         call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                    "the same tiny cell WET must still trip the viscous-CFL check")
      end block checks
   end subroutine test_audit_viscous_cfl_ignores_land_cell

   subroutine test_audit_kappa_h_spherical_uses_real_metric(error)
      !! The FIXED kappa_h check: on a spherical grid the nominal
      !! &grid_nml dx/dy are DEGREES (meaningless as a length), but the
      !! real per-cell metric (metrics%dxT/dyT, metres) still shrinks
      !! toward the pole exactly like the viscous check needs. Build a
      !! spherical sector reaching high latitude so the ACTUAL minimum
      !! cell is small enough to trip the bound, and confirm the check
      !! actually looks at metrics (not nominal dx/dy, which would be
      !! degrees and never trip a metres-scale bound).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      integer :: ierr
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 0.5_wp, LAT_S = 85.0_wp

      checks: block
         ! 6x8 sector from 85N; physical rows span ~85.25N (south edge)
         ! to ~88.75N (north edge, safely short of the pole -- the plain
         ! lon-lat formula is only valid for |lat|<90). cos(lat) shrinks
         ! fast this close to the pole: dx ~ 9.7 km at the south edge,
         ! ~2.4 km at the north edge -- exactly the FINDINGS.md near-pole
         ! shrinkage regime, and small enough that kappa_h=2000 trips the
         ! 0.5 bound (diffusive number ~2.4) while nowhere in the sector
         ! goes non-positive (which would trigger the audit's own
         ! defensive dx_min<=0 skip).
         call grid%init(6, 8, 2, DLON, DLAT)
         call make_spherical_metrics(metrics, grid, 0.0_wp, LAT_S, DLON, DLAT, 6.371e6_wp)
         call base_config(cfg, 3600.0_wp, 0.0_wp, 1.0e5_wp, 2000.0_wp, &
                          bound_kh=.false., stress_tensor=.false.)
         cfg%ocean%grid%grid_config = "spherical"

         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)

         ! kappa_h=2000 m2/s at dt_therm=3600s and dx~O(10km) gives a
         ! diffusive number well over 0.5 -- must ERROR. The OLD
         ! rdb_config.F90 check would have SKIPPED entirely on this grid
         ! (grid_config /= "cartesian"), silently missing it.
         call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                    "kappa_h diffusive-number check must fire on a spherical grid "// &
                    "using the real per-cell metric")
      end block checks
   end subroutine test_audit_kappa_h_spherical_uses_real_metric

   subroutine test_audit_ah_max_below_nu_h_warns_only(error)
      !! ah_max < nu_h silently clamps nu_h -- a real config bug (wasted
      !! two runs per FINDINGS.md/CLAUDE.md), but a WARNING not an ERROR
      !! (the run is still numerically well-defined, just not what the
      !! user asked for). Must not fail configure.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         call grid%init(8, 8, 2, 10000.0_wp, 10000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ! nu_h=20000 with ah_max=10000: ah_max < nu_h (the CLAUDE.md
         ! "wasted two runs" trap). Viscous CFL itself is tiny at this
         ! dx (10 km), so only check 4 fires.
         call base_config(cfg, 300.0_wp, 20000.0_wp, 10000.0_wp, 0.0_wp, &
                          bound_kh=.false., stress_tensor=.false.)

         call ocean_stability_audit(cfg, metrics, grid, 0, ierr=ierr)
         call destroy_cartesian_metrics(metrics)

         call check(error, ierr == OCEAN_STATUS_OK, &
                    "ah_max < nu_h must warn, not fail configure")
      end block checks
   end subroutine test_audit_ah_max_below_nu_h_warns_only

   ! -----------------------------------------------------------------
   ! Terrain-following stiffness (rx0) — Check 5
   !
   ! Motivating failure: ISOMIP+ Ocean0 on the idealised draft
   ! (`validation_examples/ocean/isomip_plus/ocean0_idealised_draft.nml`)
   ! runs clean for 3 days and then drives the barotropic free surface
   ! below the bed at the trough sidewall, where the ISOMIP+ bathymetry
   ! drops ~122 m across ONE 2 km cell and a 23 m water column sits
   ! beside a 146 m one. Nothing said so at configure time. Every number
   ! below is measured off that configuration's own `water_column`
   ! diagnostic.
   ! -----------------------------------------------------------------

   subroutine test_sigma_stiffness_formula(error)
      !! `rx0 = |H_a-H_b|/(H_a+H_b)`, pinned on the ISOMIP+ Ocean0
      !! sidewall pair, plus the two limits and the degenerate guard.
      type(error_type), allocatable, intent(out) :: error

      ! The worst wet-wet face of ocean0_idealised_draft: 23.09 m beside
      ! 145.55 m across one 2 km face on the trough sidewall.
      call check(error, abs(ocean_sigma_stiffness(23.09_wp, 145.55_wp) &
                            - 0.7261622390891841_wp) < 1.0e-12_wp, &
                 "ISOMIP+ Ocean0 sidewall regression number: rx0 = 0.726")
      if (allocated(error)) return
      call check(error, ocean_sigma_stiffness(23.09_wp, 145.55_wp) &
                 > 3.0_wp*ocean_sigma_stiffness_limit(), &
                 "that face is more than 3x the 0.2 bound")
      if (allocated(error)) return
      ! Symmetric in its arguments.
      call check(error, abs(ocean_sigma_stiffness(145.55_wp, 23.09_wp) &
                            - ocean_sigma_stiffness(23.09_wp, 145.55_wp)) == 0.0_wp, &
                 "rx0 must be symmetric across the face")
      if (allocated(error)) return
      ! Equal columns: no stiffness at all, at any depth.
      call check(error, ocean_sigma_stiffness(720.0_wp, 720.0_wp) == 0.0_wp, &
                 "two equal columns must give exactly zero")
      if (allocated(error)) return
      ! Degenerate pair => "no constraint expressible", never a NaN and
      ! never 1 (which would make every coastline the worst face).
      call check(error, ocean_sigma_stiffness(0.0_wp, 100.0_wp) == 0.0_wp, &
                 "a non-positive column returns 0, not 1 and not NaN")
      if (allocated(error)) return
      call check(error, abs(ocean_sigma_stiffness(20.0_wp, 100.0_wp) &
                            - 0.6666666666666666_wp) < 1.0e-12_wp, &
                 "plain algebra: |20-100|/120 = 2/3")
   end subroutine test_sigma_stiffness_formula

   subroutine test_sigma_stiffness_vcoord_gate(error)
      !! The check must run for every coordinate whose interfaces follow
      !! the topography and for no other. `VCOORD_ZSTAR` is IN because on
      !! the ocean path it shares the `VCOORD_SIGMA` branch of
      !! `ocean_vcoord_compute_target_h`.
      type(error_type), allocatable, intent(out) :: error
      integer :: k
      integer :: following(4), flat(4)

      following = [VCOORD_SIGMA, VCOORD_ZSTAR, VCOORD_ZSIGMA, VCOORD_ZSTAR_SIGMA]
      flat = [VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, VCOORD_Z_FIXED, VCOORD_RHO]

      do k = 1, 4
         call check(error, ocean_vcoord_is_terrain_following(following(k)), &
                    "sigma / zstar / zsigma / zstar_sigma are terrain-following")
         if (allocated(error)) return
      end do
      do k = 1, 4
         call check(error,.not. ocean_vcoord_is_terrain_following(flat(k)), &
                    "lagrangian / eulerian_z / z_fixed / rho do not follow the topography")
         if (allocated(error)) return
      end do
   end subroutine test_sigma_stiffness_vcoord_gate

   subroutine test_sigma_stiffness_worst_finds_isomip_sidewall(error)
      !! An ISOMIP+-shaped sidewall in miniature: one row of thin wet
      !! columns between grounded land and the deep trough. The scan must
      !! find THAT face, report both thicknesses and the axis, and count
      !! the faces over bound.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 9
      real(wp) :: column(NG, NG), wet(NG, NG)
      real(wp) :: rx0_max, h_thin, h_thick
      integer :: i, j, i_at, j_at, n_over, n_face
      logical :: is_x

      ! 5x5 physical interior (nghost = 2 => i,j = 3..7).
      ! j = 3   : land (grounded ice meets the bed)
      ! j = 4   : the thin sidewall row, 23.09 m
      ! j = 5..7: the trough, 145.55 m
      wet = 0.0_wp
      column = 0.0_wp
      do j = 4, 7
         do i = 3, 7
            wet(i, j) = 1.0_wp
            if (j == 4) then
               column(i, j) = 23.09_wp
            else
               column(i, j) = 145.55_wp
            end if
         end do
      end do

      call ocean_sigma_stiffness_worst(NG, NG, 3, 7, 3, 7, column, wet, &
                                       rx0_max, i_at, j_at, is_x, &
                                       h_thin, h_thick, n_over, n_face)

      call check(error, abs(rx0_max - 0.7261622390891841_wp) < 1.0e-12_wp, &
                 "the scan must find the sidewall face, rx0 = 0.726")
      if (allocated(error)) return
      call check(error,.not. is_x, &
                 "the worst face is a y face (the column jumps across j)")
      if (allocated(error)) return
      call check(error, j_at == 4, "reported at the THIN side of the face")
      if (allocated(error)) return
      call check(error, abs(h_thin - 23.09_wp) < 1.0e-12_wp .and. &
                 abs(h_thick - 145.55_wp) < 1.0e-12_wp, &
                 "both column thicknesses must be reported, thin first")
      if (allocated(error)) return
      ! Wet-wet faces: 4 rows x 4 x-faces = 16, plus 3 y-face rows x 5 = 15.
      call check(error, n_face == 31, "every interior wet-wet face is scanned")
      if (allocated(error)) return
      ! Over bound: only the five sidewall y faces (j=4 -> j=5).
      call check(error, n_over == 5, &
                 "exactly the five sidewall faces are over the 0.2 bound")
   end subroutine test_sigma_stiffness_worst_finds_isomip_sidewall

   subroutine test_sigma_stiffness_land_is_not_a_stiff_face(error)
      !! A land T-cell holds `h_layer = H_VANISHED` under the land-state
      !! contract, so its column is ~0 and a thickness-only test would
      !! read `rx0 -> 1` at every coastline. A coastline is a WALL (the
      !! face metrics are zeroed, no pressure gradient is taken), so the
      !! scan must exclude it by the WET MASK and report the interior.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 9
      real(wp) :: column(NG, NG), wet(NG, NG)
      real(wp) :: rx0_max, h_thin, h_thick
      integer :: i, j, i_at, j_at, n_over, n_face
      logical :: is_x

      wet = 0.0_wp
      column = 0.0_wp
      do j = 3, 7
         do i = 3, 7
            wet(i, j) = 1.0_wp
            column(i, j) = 300.0_wp
         end do
      end do
      ! One interior island, held at the land-state marker (36 * 1.5e-4).
      wet(5, 5) = 0.0_wp
      column(5, 5) = 36.0_wp*1.5e-4_wp

      call ocean_sigma_stiffness_worst(NG, NG, 3, 7, 3, 7, column, wet, &
                                       rx0_max, i_at, j_at, is_x, &
                                       h_thin, h_thick, n_over, n_face)

      call check(error, rx0_max == 0.0_wp, &
                 "a land neighbour must not register as a stiff face")
      if (allocated(error)) return
      call check(error, n_over == 0, "and must not be counted over bound")
      if (allocated(error)) return
      ! 5x5 block minus the island: 40 faces total, 4 of them touch the
      ! island and are dropped.
      call check(error, n_face == 36, &
                 "the island's four faces are excluded, the rest are scanned")
   end subroutine test_sigma_stiffness_land_is_not_a_stiff_face

   subroutine test_sigma_stiffness_uniform_column_is_zero(error)
      !! A flat-bed, flat-lid domain is exactly `rx0 = 0` — the control
      !! that says the check cannot fire on the geometry every
      !! bit-identity case in the corpus uses.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 9
      real(wp) :: column(NG, NG), wet(NG, NG)
      real(wp) :: rx0_max, h_thin, h_thick
      integer :: i_at, j_at, n_over, n_face
      logical :: is_x

      wet = 1.0_wp
      column = 4000.0_wp

      call ocean_sigma_stiffness_worst(NG, NG, 3, 7, 3, 7, column, wet, &
                                       rx0_max, i_at, j_at, is_x, &
                                       h_thin, h_thick, n_over, n_face)

      call check(error, rx0_max == 0.0_wp .and. n_over == 0, &
                 "a flat bed is exactly zero stiffness, bit for bit")
      if (allocated(error)) return
      call check(error, n_face == 40, "5x5 interior => 20 x faces + 20 y faces")
   end subroutine test_sigma_stiffness_uniform_column_is_zero

end module test_ocean_stability_audit
