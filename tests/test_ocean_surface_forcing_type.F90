!! Unit tests for the PR-12 surface-flux type reshape: the component set
!! on `ocean_surface_flux_t` (`rdb_ocean_surface_flux`), the assembler
!! `ocean_surface_flux_assemble` that derives `Q_heat`/`Q_salt` from it,
!! and the `stress_mag` dedup on `ocean_surface_stress_t`
!! (`rdb_ocean_surface_stress`).
!!
!! Cases (mirroring the PR-12 plan §9):
!!   1. components_off_bitident        — default off ⇒ Q_heat/Q_salt and
!!      `bytes()` untouched, apply_tracers reproduces the pre-PR analytic.
!!   2. assemble_matches_component_sum — plain sum, uniform then
!!      per-column-varying q_sw.
!!   3. isothermal_mass_exchange_adds_no_heat — the enthalpy-of-mass
!!      bookkeeping invariant (§3.1): balanced isothermal exchange adds
!!      exactly zero heat; warm rain adds exactly its own enthalpy.
!!   4. ice_components_bitident        — the ice-coupler retarget (§5.4)
!!      gives bit-identical Q_heat/Q_salt with components on vs off.
!!   5. stress_mag_matches_inline      — the §7.5 dedup is bit-identical.
!!   6. assemble_gpu_resident          — mem:separate discipline: a
!!      post-map host write must be pushed with `update device` or the
!!      assembler silently reads stale (zero) component arrays.
module test_ocean_surface_forcing_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_assemble, SEAWATER_CP
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux, ice_ocean_heat_flux
   implicit none
   private

   public :: collect_ocean_surface_forcing_type_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H_LAYER = 10.0_wp

contains

   subroutine collect_ocean_surface_forcing_type_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("components_off_bitident", test_components_off_bitident), &
                  new_unittest("assemble_matches_component_sum", &
                               test_assemble_matches_component_sum), &
                  new_unittest("isothermal_mass_exchange_adds_no_heat", &
                               test_isothermal_mass_exchange), &
                  new_unittest("ice_components_bitident", test_ice_components_bitident), &
                  new_unittest("stress_mag_matches_inline", test_stress_mag_matches_inline), &
                  new_unittest("assemble_gpu_resident", test_assemble_gpu_resident) &
                  ]
   end subroutine collect_ocean_surface_forcing_type_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine seed_uniform(ms, t0, s0)
      !! Uniform h + T + S across every layer (including ghosts) — the
      !! shared IC for every case below.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: t0, s0
      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = t0*H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = s0*H_LAYER
   end subroutine seed_uniform

   ! -----------------------------------------------------------------
   ! Case 1 — components off ⇒ byte-identical
   ! -----------------------------------------------------------------

   subroutine test_components_off_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf, sf_ref
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_HEAT = 100.0_wp
      real(wp), parameter :: Q_SALT = 1.0e-5_wp
      real(wp), allocatable :: Q_heat_before(:, :), Q_salt_before(:, :)
      real(wp) :: expected_dhT, top_err
      integer :: ip, jp
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call sf%init(grid)
         call sf_ref%init(grid)
         call seed_uniform(ms, 10.0_wp, 35.0_wp)

         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)
         allocate (Q_heat_before, source=sf%Q_heat)
         allocate (Q_salt_before, source=sf%Q_salt)

         ! use_components stays .false. (never touched) — the assembler
         ! must be a plain early return.
         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%Q_heat, sf%Q_salt)
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         call check(error, maxval(abs(sf%Q_heat - Q_heat_before)) == 0.0_wp, &
                    "assemble() with use_components=.false. must not touch Q_heat")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf%Q_salt - Q_salt_before)) == 0.0_wp, &
                    "assemble() with use_components=.false. must not touch Q_salt")
         if (allocated(error)) exit checks

         call check(error, sf%bytes() == sf_ref%bytes(), &
                    "components-off bytes() must equal an init-only slot's bytes()")
         if (allocated(error)) exit checks

         ! The pre-PR-12 apply_tracers analytic must still hold exactly.
         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
         associate (hT => ms%tracers(ms%idx_temperature)%hTr)
            !$acc update self(hT)
         end associate
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         expected_dhT = DT*Q_HEAT/(sf%rho0*sf%cp)
         top_err = abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ) &
                       - (10.0_wp*H_LAYER + expected_dhT))
         call check(error, top_err < 1.0e-10_wp, &
                    "components-off apply_tracers must match the pre-PR-12 analytic")

      end block checks
      if (allocated(Q_heat_before)) deallocate (Q_heat_before)
      if (allocated(Q_salt_before)) deallocate (Q_salt_before)
      call sf%destroy(); call sf_ref%destroy(); call ms%destroy()
   end subroutine test_components_off_bitident

   ! -----------------------------------------------------------------
   ! Case 2 — plain component sum
   ! -----------------------------------------------------------------

   subroutine test_assemble_matches_component_sum(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: Q_HEAT_CONST = 300.0_wp
      real(wp) :: expected_uniform
      integer :: i, j, nx
      real(wp) :: expected_cell
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call sf%set_surface_flux_const(Q_HEAT_CONST, 0.0_wp)
         call seed_uniform(ms, 10.0_wp, 35.0_wp)

         sf%q_sw = 200.0_wp
         sf%q_lw = -60.0_wp
         sf%q_lat = -90.0_wp
         sf%q_sens = -30.0_wp
         sf%heat_added = 5.0_wp
         ! evap / mass fluxes stay zero => heat_content_massin/massout = 0.

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%Q_heat)
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         expected_uniform = Q_HEAT_CONST + 25.0_wp   ! 200-60-90-30+5
         call check(error, maxval(abs(sf%Q_heat - expected_uniform)) < 1.0e-9_wp, &
                    "Q_heat must equal Q_heat_const + sum(components) per cell")
         if (allocated(error)) exit checks

         ! ---- Non-uniform q_sw: per-column read, not a broadcast ----
         nx = grid%nx_total
         do j = 1, grid%ny_total
            do i = 1, nx
               sf%q_sw(i, j) = 10.0_wp*real(i, wp)
            end do
         end do

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%Q_heat)
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         do j = 1, grid%ny_total
            do i = 1, nx
               expected_cell = Q_HEAT_CONST + 10.0_wp*real(i, wp) - 60.0_wp - 90.0_wp - 30.0_wp + 5.0_wp
               call check(error, abs(sf%Q_heat(i, j) - expected_cell) < 1.0e-9_wp, &
                          "Q_heat must track a per-column q_sw field, not a broadcast")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      call sf%destroy(); call ms%destroy()
   end subroutine test_assemble_matches_component_sum

   ! -----------------------------------------------------------------
   ! Case 3 — the enthalpy-of-mass bookkeeping invariant
   ! -----------------------------------------------------------------

   subroutine test_isothermal_mass_exchange(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: T_S = 15.0_wp
      real(wp), parameter :: H_TOP = 1.0_wp   !! exact SST recovery: hTr/h == T_s
      real(wp), parameter :: E = 1.0e-5_wp
      real(wp), parameter :: Q_HEAT_CONST = 42.0_wp
      real(wp) :: massin_expect, massout_expect
      integer :: ip, jp
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call sf%set_surface_flux_const(Q_HEAT_CONST, 0.0_wp)

         ms%h_layer = H_LAYER
         ms%h_layer(:, :, NZ) = H_TOP
         ms%tracers(ms%idx_temperature)%hTr = T_S*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = T_S*H_TOP
         ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*H_LAYER

         sf%evap = -E
         sf%lprec = E
         sf%heat_content_lprec = SEAWATER_CP*T_S*E   !! isothermal: rain at T_s

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%Q_heat, sf%heat_content_massin, sf%heat_content_massout)
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         massin_expect = SEAWATER_CP*T_S*E
         massout_expect = SEAWATER_CP*T_S*(-E)

         call check(error, abs(sf%heat_content_massin(ip, jp) - massin_expect) &
                    < 1.0e-12_wp*abs(massin_expect), &
                    "heat_content_massin must equal heat_content_lprec (only nonzero term)")
         if (allocated(error)) exit checks
         call check(error, abs(sf%heat_content_massout(ip, jp) - massout_expect) &
                    < 1.0e-12_wp*abs(massout_expect), &
                    "heat_content_massout must equal SEAWATER_CP*T_sst*evap")
         if (allocated(error)) exit checks
         call check(error, abs(sf%heat_content_massin(ip, jp) + sf%heat_content_massout(ip, jp)) &
                    < 1.0e-9_wp*abs(massin_expect), &
                    "balanced isothermal exchange must add ~zero net enthalpy")
         if (allocated(error)) exit checks
         call check(error, abs(sf%Q_heat(ip, jp) - Q_HEAT_CONST) < 1.0e-9_wp, &
                    "isothermal balanced exchange must leave Q_heat == Q_heat_const")

      end block checks
      call sf%destroy(); call ms%destroy()
      if (allocated(error)) return

      ! ---- Warm rain: same setup, heat_content_lprec 10 degC warmer ----
      warm: block
         type(hgrid_t) :: grid2
         type(multilayer_state_t) :: ms2
         type(ocean_surface_flux_t) :: sf2
         real(wp), parameter :: DT_WARM = 10.0_wp
         real(wp) :: dq_expect
         integer :: ip2, jp2

         call make_grid(grid2, 6, 4)
         ms2%nz_ml = NZ
         call ms2%init(grid2)
         call sf2%init(grid2)
         call sf2%set_components(grid2, .true.)
         call sf2%set_surface_flux_const(Q_HEAT_CONST, 0.0_wp)

         ms2%h_layer = H_LAYER
         ms2%h_layer(:, :, NZ) = H_TOP
         ms2%tracers(ms2%idx_temperature)%hTr = T_S*H_LAYER
         ms2%tracers(ms2%idx_temperature)%hTr(:, :, NZ) = T_S*H_TOP
         ms2%tracers(ms2%idx_salinity)%hTr = 35.0_wp*H_LAYER

         sf2%evap = -E
         sf2%lprec = E
         sf2%heat_content_lprec = SEAWATER_CP*(T_S + DT_WARM)*E

         !$acc enter data copyin(ms2, sf2)
         call ms2%enter_data()
         call sf2%enter_data()
         call ocean_surface_flux_assemble(grid2, sf2, ms2)
         !$acc update self(sf2%Q_heat)
         call sf2%exit_data()
         call ms2%exit_data()
         !$acc exit data delete(ms2, sf2)

         ip2 = grid2%nx_total/2
         jp2 = grid2%ny_total/2
         dq_expect = SEAWATER_CP*DT_WARM*E
         call check(error, abs((sf2%Q_heat(ip2, jp2) - Q_HEAT_CONST) - dq_expect) &
                    < 1.0e-9_wp*abs(dq_expect), &
                    "warm rain must warm the ocean by exactly the enthalpy it carries")

         call sf2%destroy(); call ms2%destroy()
      end block warm
   end subroutine test_isothermal_mass_exchange

   ! -----------------------------------------------------------------
   ! Case 4 — ice-coupler retarget bit-identity
   ! -----------------------------------------------------------------

   subroutine test_ice_components_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf_off, sf_on
      real(wp), parameter :: Q_HEAT_CONST = 10.0_wp
      real(wp), parameter :: Q_SALT_CONST = 1.0e-6_wp
      real(wp), parameter :: SFLUX = 2.0e-6_wp
      real(wp), parameter :: HFLUX = 50.0_wp
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call seed_uniform(ms, 5.0_wp, 35.0_wp)

         ice%enable = .true.
         call ice%init(grid)
         ice%salt_flux_diag = SFLUX
         ice%heat_flux_diag = HFLUX

         call sf_off%init(grid)
         call sf_off%set_surface_flux_const(Q_HEAT_CONST, Q_SALT_CONST)

         call sf_on%init(grid)
         call sf_on%set_components(grid, .true.)
         call sf_on%set_surface_flux_const(Q_HEAT_CONST, Q_SALT_CONST)

         !$acc enter data copyin(ms, ice, sf_off, sf_on)
         call ms%enter_data()
         call ice%enter_data()
         call sf_off%enter_data()
         call sf_on%enter_data()

         call ice_ocean_brine_flux(sf_off, ice)
         call ice_ocean_heat_flux(sf_off, ice)
         call ice_ocean_brine_flux(sf_on, ice)
         call ice_ocean_heat_flux(sf_on, ice)
         call ocean_surface_flux_assemble(grid, sf_on, ms)

         !$acc update self(sf_off%Q_heat, sf_off%Q_salt, sf_on%Q_heat, sf_on%Q_salt, &
         !$acc&            sf_on%heat_added, sf_on%salt_flux)

         call sf_on%exit_data()
         call sf_off%exit_data()
         call ice%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ice, sf_off, sf_on)

         call check(error, sf_off%has_heat .and. sf_on%has_heat, &
                    "has_heat must be latched true on both paths")
         if (allocated(error)) exit checks
         call check(error, sf_off%has_salt .and. sf_on%has_salt, &
                    "has_salt must be latched true on both paths")
         if (allocated(error)) exit checks

         call check(error, maxval(abs(sf_off%Q_heat - sf_on%Q_heat)) == 0.0_wp, &
                    "components on vs off must give bit-identical Q_heat under ice")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf_off%Q_salt - sf_on%Q_salt)) == 0.0_wp, &
                    "components on vs off must give bit-identical Q_salt under ice")
         if (allocated(error)) exit checks

         call check(error, maxval(abs(sf_on%heat_added - HFLUX)) == 0.0_wp, &
                    "components-on path must write the ice heat into the heat_added COMPONENT")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf_on%salt_flux - SFLUX)) == 0.0_wp, &
                    "components-on path must write the ice salt into the salt_flux COMPONENT")

      end block checks
      call sf_on%destroy(); call sf_off%destroy(); call ice%destroy(); call ms%destroy()
   end subroutine test_ice_components_bitident

   ! -----------------------------------------------------------------
   ! Case 5 — stress_mag dedup bit-identity
   ! -----------------------------------------------------------------

   subroutine test_stress_mag_matches_inline(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAUX_MAG = 0.1_wp
      real(wp) :: tau_x_cell, tau_y_cell, expected
      integer :: i, j
      checks: block

         call make_grid(grid, 8, 6)
         call ss%init(grid)
         call ss%set_wind_stress_2gyre(grid, TAUX_MAG)

         call check(error, ss%has_stress_mag, &
                    "has_stress_mag must be latched true after a wind-stress setter")
         if (allocated(error)) exit checks

         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               tau_x_cell = 0.5_wp*(ss%tau_x(i, j) + ss%tau_x(i + 1, j))
               tau_y_cell = 0.5_wp*(ss%tau_y(i, j) + ss%tau_y(i, j + 1))
               expected = sqrt(tau_x_cell*tau_x_cell + tau_y_cell*tau_y_cell)
               call check(error, ss%stress_mag(i, j) == expected, &
                          "stress_mag must match the inline 3-line computation bit-for-bit")
               if (allocated(error)) exit checks
            end do
            if (allocated(error)) exit checks
         end do

      end block checks
      call ss%destroy()
   end subroutine test_stress_mag_matches_inline

   ! -----------------------------------------------------------------
   ! Case 6 — mem:separate discipline
   ! -----------------------------------------------------------------

   subroutine test_assemble_gpu_resident(error)
      !! Following `test_open_boundary_out_closes`'s template: map,
      !! write a component AFTER the map, push it with `update device`
      !! explicitly (the exact contract every non-reader filler must
      !! follow — module docstring item (b)), run the assembler,
      !! `update self` the result, and assert it is non-zero and
      !! correct.  Without the `update device`, `Q_heat` would come
      !! back bit-identical to `Q_heat_const` (the stale-zero trap) —
      !! this is the gate a missed map/push would fail.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: Q_HEAT_CONST = 0.0_wp
      real(wp), parameter :: Q_SW_VAL = 250.0_wp
      integer :: ip, jp
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)
         call sf%set_surface_flux_const(Q_HEAT_CONST, 0.0_wp)
         call seed_uniform(ms, 10.0_wp, 35.0_wp)

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()

         ! Post-map host write — the exact scenario the contract covers.
         sf%q_sw = Q_SW_VAL
         !$acc update device(sf%q_sw)

         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%Q_heat)

         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, sf%Q_heat(ip, jp) /= 0.0_wp, &
                    "Q_heat must be non-zero — a missed 'update device' would leave it at 0")
         if (allocated(error)) exit checks
         call check(error, abs(sf%Q_heat(ip, jp) - Q_SW_VAL) < 1.0e-9_wp, &
                    "Q_heat must equal the pushed q_sw value (Q_heat_const=0, all else 0)")

      end block checks
      call sf%destroy(); call ms%destroy()
   end subroutine test_assemble_gpu_resident

end module test_ocean_surface_forcing_type
