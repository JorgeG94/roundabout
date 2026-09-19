!! The forcing + KPP reference densities must follow the SINGLE configured ρ₀.
!!
!! The defect these pin
!! --------------------
!! `&ocean_ic_nml rho_0` lands on `eos%rho0` — the one ρ₀ of record.  Eleven
!! slots already copy it in their `configure_ocean_*` (EPBL, kappa-shear,
!! tidal mixing, wave speed, the `eta_ib` seam, GM / MEKE / Redi / MLE, the
!! isopycnal slopes, and — since the preceding commit — the PGF).  Four did
!! NOT: `surface_flux%rho0`, `surface_stress%rho0`, `vmix%rho0` and
!! `geothermal%rho0` each declared `= 1035.0_wp` and nothing ever assigned
!! them.  A namelist setting `rho_0 /= 1035` therefore ran its equation of
!! state on the configured density while
!!
!!   * every surface heat source used `dt/(1035·cp)` and every surface salt
!!     source `dt/1035` — including whatever the sea-ice coupler delivered
!!     through `Q_heat`/`Q_salt`;
!!   * every wind-stress acceleration used `τ/(1035·h_top)`;
!!   * KPP's N², u* = √(|τ|/ρ₀) and the kinematic fluxes behind B_0 used 1035;
!!   * the geothermal bed source used `dt·Q_geo/(1035·cp)`
!!
!! — silently, with no warning and no fail-loud.
!!
!! Cases
!! -----
!!   * `configured_rho0_reaches_every_forcing_slot` — the plumbing proper.
!!     A config with `rho_0 = RHO0_CFG` (3.4 % off the default, far outside
!!     any tolerance here) is driven through the production path
!!     (`init_from_config` → `configure_ocean_reference_density`) and all
!!     four slots are asserted equal to `eos%rho0`.  NON-VACUITY: the
!!     configured value is first asserted to differ from the 1035 type
!!     default, so no assertion can pass by both sides being 1035.
!!   * `surface_heat_flux_warming_uses_the_configured_rho0` — the analytical
!!     arm, where ρ₀ is the answer rather than a stored scalar.  The SAME
!!     configured state is given a uniform top-layer temperature and one
!!     constant-`Q_heat` apply step, whose closed form is
!!
!!         Δ(hT)|_{k=nz} = Q · dt / (ρ₀ · cp)
!!
!!     i.e. ΔT = Q·dt/(ρ₀·cp·h_top).  ρ₀ enters ONLY as that divisor, so the
!!     kernel answer is asserted against the closed form at the CONFIGURED
!!     ρ₀ and, for non-vacuity, shown to be resolvably away from the same
!!     form at 1035.  Layers below the surface must not move (bottom-up
!!     convention: surface = `k = nz`).
!!
!! `mem:separate` contract: `run_apply` maps the multilayer slot AND its
!! surface-flux companion before the kernel, pushes the host-set inputs with
!! `!$acc update device`, and pulls the answer back with `!$acc update self`
!! on the COMPONENT arrays (never the aggregate — an aggregate D→H copy
!! overwrites the host descriptors).  All of it is an inert no-op on
!! host / multicore builds.
module test_ocean_forcing_rho_ref
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_reference_density
   use rdb_ocean_geothermal, only: ocean_geothermal_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_apply_tracers
   implicit none
   private

   public :: collect_ocean_forcing_rho_ref_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: RHO0_DEFAULT = 1035.0_wp
      !! The literal every one of these slots defaulted to, and the
      !! `&ocean_ic_nml rho_0` default.  Named here only so the non-vacuity
      !! arms can show the test is not comparing 1035 against 1035.
   real(wp), parameter :: RHO0_CFG = 1000.0_wp
      !! Deliberately NOT the default.
   real(wp), parameter :: TOL_SCALAR = 1.0e-12_wp
   real(wp), parameter :: H_LAYER = 10.0_wp
   real(wp), parameter :: T_IC = 4.0_wp
   real(wp), parameter :: DT = 3600.0_wp
   real(wp), parameter :: Q_HEAT = 100.0_wp
      !! W/m², positive down.  Non-zero is what makes the analytical arm
      !! non-vacuous — `has_heat` only latches on a non-zero fill.

contains

   subroutine collect_ocean_forcing_rho_ref_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("configured_rho0_reaches_every_forcing_slot", &
                               test_forcing_rho0_plumbing), &
                  new_unittest("surface_heat_flux_warming_uses_the_configured_rho0", &
                               test_surface_heat_answer) &
                  ]
   end subroutine collect_ocean_forcing_rho_ref_tests

   ! ---------------------------------------------------------------------
   ! Scaffolding
   ! ---------------------------------------------------------------------

   subroutine build_configured_state(cfg, state, grid, geo)
      !! The production wiring, nothing test-specific: fill a config with
      !! `rho_0 = RHO0_CFG`, allocate the god state through
      !! `init_from_config` (which is what lands `rho_0` on `eos%rho0`),
      !! init the engine-held geothermal slot, then run the one configure
      !! stage that owns these reference densities.
      type(config_t), intent(inout) :: cfg
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(inout) :: grid
      type(ocean_geothermal_t), intent(inout) :: geo

      cfg%sim_type = "ocean"
      cfg%nx = NX_PHYS
      cfg%ny = NY_PHYS
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%nghost = NGHOST
      cfg%ocean%ic%rho_0 = RHO0_CFG

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      call state%init_from_config(cfg, grid)
      call geo%init(grid)
      call configure_ocean_reference_density(state, geo=geo)
   end subroutine build_configured_state

   subroutine run_apply(grid, state)
      !! One surface-flux apply over the configured god state's own slots.
      !! See the module header for the `mem:separate` contract this follows.
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(inout) :: state

      !$acc enter data copyin(state%multilayer, state%surface_flux)
      call state%multilayer%enter_data()
      call state%surface_flux%enter_data()
      !$acc update device(state%multilayer%h_layer, state%surface_flux%Q_heat)
      call ocean_surface_flux_apply_tracers(grid, state%surface_flux, &
                                            state%multilayer, DT)
      associate (hT => state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)
         !$acc update self(hT)
      end associate
      call state%surface_flux%exit_data()
      call state%multilayer%exit_data()
      !$acc exit data delete(state%multilayer, state%surface_flux)
   end subroutine run_apply

   ! ---------------------------------------------------------------------
   ! Cases
   ! ---------------------------------------------------------------------

   subroutine test_forcing_rho0_plumbing(error)
      !! `&ocean_ic_nml rho_0` must reach all four remaining reference
      !! densities, and they must agree with the EOS.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      type(ocean_geothermal_t) :: geo

      checks: block
         call build_configured_state(cfg, state, grid, geo)

         ! Non-vacuity: the configured value must not be the type default,
         ! or every assertion below would hold with the assignments missing.
         call check(error, abs(RHO0_CFG - RHO0_DEFAULT) > 1.0_wp, &
                    "test setup is vacuous: RHO0_CFG equals the 1035 default")
         if (allocated(error)) exit checks

         call check(error, abs(state%eos%rho0 - RHO0_CFG) < TOL_SCALAR, &
                    "eos%rho0 did not follow &ocean_ic_nml rho_0")
         if (allocated(error)) exit checks
         call check(error, abs(state%surface_flux%rho0 - state%eos%rho0) < TOL_SCALAR, &
                    "surface_flux%rho0 did not follow the configured rho_0 — every "// &
                    "surface heat/salt source (sea ice included) is scaled by 1035/rho_0")
         if (allocated(error)) exit checks
         call check(error, abs(state%surface_stress%rho0 - state%eos%rho0) < TOL_SCALAR, &
                    "surface_stress%rho0 did not follow the configured rho_0 — the "// &
                    "wind-stress acceleration tau/(rho0*h) is scaled by 1035/rho_0")
         if (allocated(error)) exit checks
         call check(error, abs(state%vmix%rho0 - state%eos%rho0) < TOL_SCALAR, &
                    "vmix%rho0 did not follow the configured rho_0 — KPP N2, u* and "// &
                    "the kinematic fluxes behind B_0 run on a different rho0 than the EOS")
         if (allocated(error)) exit checks
         call check(error, abs(geo%rho0 - state%eos%rho0) < TOL_SCALAR, &
                    "geothermal%rho0 did not follow the configured rho_0")
      end block checks
      call geo%destroy()
      call state%destroy()
   end subroutine test_forcing_rho0_plumbing

   subroutine test_surface_heat_answer(error)
      !! The ANSWER, not just the scalar.  One constant-`Q_heat` apply step
      !! on a uniform column: Δ(hT) at `k = nz` is exactly `Q·dt/(ρ₀·cp)`,
      !! whose only ρ₀ dependence is that divisor.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      type(ocean_geothermal_t) :: geo
      integer :: i_probe, j_probe, k, idx_T
      real(wp) :: want_cfg, want_default, got, drift, discrimination

      checks: block
         call build_configured_state(cfg, state, grid, geo)

         idx_T = state%multilayer%idx_temperature
         call check(error, idx_T > 0, "no temperature tracer registered")
         if (allocated(error)) exit checks

         state%multilayer%h_layer = H_LAYER
         state%multilayer%tracers(idx_T)%hTr = T_IC*H_LAYER
         call state%surface_flux%set_surface_flux_const(Q_HEAT, 0.0_wp)

         call run_apply(grid, state)

         ! Closed form at the configured rho0 and at the old hard-coded one.
         want_cfg = DT*Q_HEAT/(RHO0_CFG*state%surface_flux%cp)
         want_default = DT*Q_HEAT/(RHO0_DEFAULT*state%surface_flux%cp)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         got = state%multilayer%tracers(idx_T)%hTr(i_probe, j_probe, NZ) - T_IC*H_LAYER

         ! Non-vacuity 1: the step must actually have heated something.
         call check(error, want_cfg > 1.0e-6_wp, &
                    "test setup is vacuous: the expected top-layer warming is zero")
         if (allocated(error)) exit checks

         call check(error, abs(got - want_cfg) < 1.0e-12_wp*want_cfg, &
                    "the top-layer warming does not use the configured rho_0 "// &
                    "(surface_flux%rho0 left at the 1035 type default)")
         if (allocated(error)) exit checks

         ! Non-vacuity 2: the 1035 answer is far outside that tolerance, so
         ! the assertion above genuinely discriminates the two divisors.
         discrimination = abs(want_default - want_cfg)
         call check(error, discrimination > 1.0e-3_wp*want_cfg, &
                    "test is vacuous: the 1035 answer is within tolerance of the "// &
                    "configured-rho_0 answer")
         if (allocated(error)) exit checks

         ! Bottom-up convention: only the surface layer is forced.
         drift = 0.0_wp
         do k = 1, NZ - 1
            drift = max(drift, abs(state%multilayer%tracers(idx_T)%hTr(i_probe, j_probe, k) &
                                   - T_IC*H_LAYER))
         end do
         call check(error, drift < 1.0e-14_wp, &
                    "layers below the surface drifted under a surface-only heat flux")
      end block checks
      call geo%destroy()
      call state%destroy()
   end subroutine test_surface_heat_answer

end module test_ocean_forcing_rho_ref
