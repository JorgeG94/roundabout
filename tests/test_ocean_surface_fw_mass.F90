!! Real fresh-water MASS at the ocean surface, `&ocean_forcing_nml freshwater = "mass"`
!! (`ocean_surface_fw_mass_step`, module `rdb_ocean_surface_fw_mass`).  Analytical gates:
!!
!!   1. `rain_grows_column_exactly`: N steps of uniform rain: the top layer grows by N*P*dt/rho0,
!!      T keeps its exact bits, S = S0*h0/h, the content h*S is unchanged, the tracked source is
!!      RHO_WATER * sum(dh * area).
!!   2. `evap_concentrates_salt_t_unchanged`: the same with evaporation: S rises as h0/h.
!!   3. `evap_cfl_passdown`: a 1 mm top layer under a large evaporation ends at 0.2*h0, the rest is
!!      drawn from the next live layer, the column total is exact, every layer stays positive, the
!!      pass-down is counted and flagged per column.
!!   4. `column_exhaustion_is_fatal`: a stack that cannot supply the request is clamped at
!!      `(1-f)*h` per layer, counted over OWNED columns only, and the verdict is fatal.
!!   5. `passive_and_age_diluted` and `pseudo_salt_matches_salinity_bitwise`: dissolved tracers
!!      keep their content; the pseudo-salt mirror equals salinity bit for bit.
!!   6. `b0_identical_virtual_vs_mass`: `Q_salt` (the boundary-layer buoyancy forcing) is bitwise
!!      the same under "virtual" and "mass"; `Q_salt_tracer = Q_salt - salt_fw`.
!!   7. `virtual_vs_mass_second_order`: one step of each: the salinities differ by exactly
!!      `S0 eps^2/(1+eps)` (the cavity identity), not by zero.
!!   8. `mass_budget_closes` / 9. `heat_budget_closes`: the console residual
!!      `(M - M0) + out - src` through the production reducers, ssp_rk2, pred_corr, ratio 2 and the
!!      deferred (`diabatic_first = .false.`) order; the source is non-zero.
!!  10. `uniform_rain_closed_basin_no_flow`: eta rises by P*t/rho0 everywhere, no current.
!!  11. `zstar_interfaces_after_rain`: after the remap the interfaces equal the z* target built from
!!      H and eta = sum(h) - H (a stale `bt_eta` builds it on H + dh).
!!  12. `none_and_virtual_bitidentical`: the assembler under "none" and "virtual" reproduces the
!!      legacy sums bitwise, and leaves `fw_net` / `Q_salt_tracer` untouched.
!!  13. `refusals_by_name`: the configure-time envelope.
module test_ocean_surface_fw_mass
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_efp, only: efp_t, efp_to_real
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_assemble, &
                                     ocean_surface_flux_fw_virtual_step, &
                                     SURF_FW_NONE, SURF_FW_VIRTUAL, SURF_FW_MASS
   use rdb_ocean_surface_fw_mass, only: ocean_surface_fw_mass_step, fw_mass_column_impl, &
                                        fw_mass_thin_is_fatal
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_console_stats, only: compute_total_h, compute_total_tracer, compute_total_content, &
                                      ocean_budget_src, ocean_budget_out, &
                                      ocean_budget_stage_weight, &
                                      ocean_salt_src_sum, ocean_heat_src_sum
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_surface_fw_mass_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   integer, parameter :: NXP = 6
   integer, parameter :: NYP = 4
   real(wp), parameter :: H_LAYER = 10.0_wp
   real(wp), parameter :: DT = 3600.0_wp
   real(wp), parameter :: AREA = 1.0e6_wp
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: T0 = 8.0_wp
   character(len=*), parameter :: NL = new_line("a")
   logical, save :: comm_inited = .false.

contains

   subroutine collect_ocean_surface_fw_mass_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fwm_rain_grows_column_exactly", test_rain_grows_column_exactly), &
                  new_unittest("fwm_evap_concentrates_salt_t_unchanged", test_evap_concentrates), &
                  new_unittest("fwm_evap_cfl_passdown", test_evap_cfl_passdown), &
                  new_unittest("fwm_column_exhaustion_is_fatal", test_column_exhaustion_is_fatal), &
                  new_unittest("fwm_evap_keeps_salt_at_top_of_cfl_range", test_evap_keeps_salt), &
                  new_unittest("fwm_passive_and_age_diluted", test_passive_and_age_diluted), &
                  new_unittest("fwm_pseudo_salt_matches_salinity_bitwise", test_pseudo_salt_bitwise), &
                  new_unittest("fwm_b0_identical_virtual_vs_mass", test_b0_identical), &
                  new_unittest("fwm_virtual_vs_mass_second_order", test_second_order), &
                  new_unittest("fwm_none_and_virtual_bitidentical", test_none_virtual_bitident), &
                  new_unittest("fwm_refusals_by_name", test_refusals), &
                  new_unittest("fwm_mass_budget_closes", test_mass_budget_closes), &
                  new_unittest("fwm_absolute_time_factor", test_absolute_time_factor), &
                  new_unittest("fwm_heat_budget_closes", test_heat_budget_closes), &
                  new_unittest("fwm_uniform_rain_closed_basin_no_flow", test_no_flow), &
                  new_unittest("fwm_zstar_interfaces_after_rain", test_zstar_interfaces), &
                  new_unittest("fwm_diagnostics_registered", test_diagnostics) &
                  ]
   end subroutine collect_ocean_surface_fw_mass_tests

   ! -----------------------------------------------------------------
   ! unit-level scaffolding
   ! -----------------------------------------------------------------

   subroutine setup_column(grid, ms, sf, at, mode, extra_tracers)
      !! A components-on slot over a 6x4 tile, uniform h, stratified T and S.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(out) :: at(:, :)
      integer, intent(in) :: mode
      logical, intent(in), optional :: extra_tracers
      integer :: k, idx

      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid, with_ideal_age=present(extra_tracers))
      if (present(extra_tracers)) then
         call ms%register_passive_tracer(grid, "dye", "1", "passive dye", idx)
         call ms%register_passive_tracer(grid, "pseudo_salt", "psu", "pseudo salt", ms%idx_pseudo_salt)
      end if
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      call sf%set_freshwater(mode)
      allocate (at(grid%nx_total, grid%ny_total), source=AREA)
      ms%h_layer = H_LAYER
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%c(:, :, k) = T0 + 0.5_wp*real(k, wp)
         ms%tracers(ms%idx_salinity)%c(:, :, k) = S0 - 0.1_wp*real(k, wp)
      end do
      if (present(extra_tracers)) then
         do k = 1, size(ms%tracers)
            if (k == ms%idx_temperature .or. k == ms%idx_salinity) cycle
            ms%tracers(k)%c = 1.0_wp + 0.01_wp*real(k, wp)
         end do
         ms%tracers(ms%idx_pseudo_salt)%c = ms%tracers(ms%idx_salinity)%c
      end if
      sf%has_heat = .true.
   end subroutine setup_column

   subroutine map_all(ms, sf, at)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(inout) :: at(:, :)

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      !$acc enter data copyin(at)
   end subroutine map_all

   subroutine unmap_all(ms, sf, at)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(inout) :: at(:, :)

      !$acc exit data delete(at)
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine unmap_all

   subroutine pull_state(ms, sf)
      !! Device -> host for every array the gates read.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer :: it

      !$acc update self(ms%h_layer, ms%heat_budget_surface, ms%salt_budget_surface)
      do it = 1, size(ms%tracers)
         associate (cc => ms%tracers(it)%c)
            !$acc update self(cc)
         end associate
      end do
      !$acc update self(sf%fw_net, sf%Q_salt, sf%Q_salt_tracer, sf%salt_fw, sf%Q_heat)
      !$acc update self(ms%fw_dh, ms%fw_passdown)
   end subroutine pull_state

   subroutine one_thermo_step(grid, ms, sf, at, weight)
      !! The engine's per-thermo-step order for ONE snapshot: the virtual salt writer, the
      !! assembler, the tracer deposit, the mass step.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: at(grid%nx_total, grid%ny_total)
      real(wp), intent(in), optional :: weight
      real(wp) :: w

      w = 1.0_wp
      if (present(weight)) w = weight
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
      call ocean_surface_fw_mass_step(grid, sf, ms, at, DT, w)
   end subroutine one_thermo_step

   function src_mass(ms) result(m)
      !! The tracked mass source (kg), exactly what the console folds in.
      type(multilayer_state_t), intent(in) :: ms
      real(wp) :: m
      type(efp_t) :: e
      e%v = ms%mass_src_efp
      e%poison = ms%mass_src_efp_poison
      m = real(efp_to_real(e), wp)
   end function src_mass

   ! -----------------------------------------------------------------
   ! 1. rain
   ! -----------------------------------------------------------------

   subroutine test_rain_grows_column_exactly(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: P = 1.0e-2_wp
      integer, parameter :: NSTEP = 10
      real(wp) :: dh, expect_h, c_before, c_after, worst_s, worst_c, vol
      real(wp), allocatable :: s_ini(:), t_ini(:)
      integer :: n, i, j, k, nowned

      call setup_column(grid, ms, sf, at, SURF_FW_MASS)
      allocate (s_ini(NZ), t_ini(NZ))
      s_ini = ms%tracers(ms%idx_salinity)%c(NGHOST + 1, NGHOST + 1, :)
      t_ini = ms%tracers(ms%idx_temperature)%c(NGHOST + 1, NGHOST + 1, :)
      sf%lprec = P
      sf%has_mass_flux = .true.
      call map_all(ms, sf, at)
      !$acc update device(sf%lprec)
      do n = 1, NSTEP
         call one_thermo_step(grid, ms, sf, at)
      end do
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at)

      dh = P*DT/sf%rho0
      expect_h = H_LAYER + real(NSTEP, wp)*dh
      worst_s = 0.0_wp
      worst_c = 0.0_wp
      nowned = 0
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            call check(error, abs(ms%h_layer(i, j, NZ) - expect_h) < 1.0e-14_wp*expect_h, &
                       "the top layer grows by N P dt/rho0 (ghosts included)")
            if (allocated(error)) return
            do k = 1, NZ - 1
               call check(error, ms%h_layer(i, j, k) == H_LAYER, "layers below the surface keep h exactly")
               if (allocated(error)) return
            end do
            do k = 1, NZ
               call check(error, ms%tracers(ms%idx_temperature)%c(i, j, k) == t_ini(k), &
                          "T keeps its exact bits (the rain joins at the layer temperature)")
               if (allocated(error)) return
            end do
            worst_s = max(worst_s, abs(ms%tracers(ms%idx_salinity)%c(i, j, NZ) - s_ini(NZ)*H_LAYER/expect_h)/s_ini(NZ))
            c_before = sum(s_ini*H_LAYER)
            c_after = 0.0_wp
            do k = 1, NZ
               c_after = c_after + ms%tracers(ms%idx_salinity)%c(i, j, k)*ms%h_layer(i, j, k)
            end do
            worst_c = max(worst_c, abs(c_after - c_before)/c_before)
         end do
      end do
      call check(error, worst_s < 1.0e-14_wp, "S = S0 h0/h on the surface layer")
      if (allocated(error)) return
      call check(error, worst_c < 1.0e-14_wp, "the salt content h*S of every column is unchanged")
      if (allocated(error)) return
      do k = 1, NZ - 1
         call check(error, ms%tracers(ms%idx_salinity)%c(NGHOST + 2, NGHOST + 2, k) == s_ini(k), &
                    "S below the surface layer is untouched")
         if (allocated(error)) return
      end do
      nowned = NXP*NYP
      vol = real(nowned, wp)*real(NSTEP, wp)*dh*AREA
      call check(error, abs(src_mass(ms) - RHO_WATER*vol) < 1.0e-13_wp*RHO_WATER*vol, &
                 "the tracked source is RHO_WATER * sum(dh * area) over the owned cells")
      if (allocated(error)) return
      call check(error, ms%fw_n_thin_step == 0_int64 .and. ms%fw_n_passdown_step == 0_int64, &
                 "rain never passes down or runs dry")
      call sf%destroy(); call ms%destroy()
   end subroutine test_rain_grows_column_exactly

   ! -----------------------------------------------------------------
   ! 2. evaporation
   ! -----------------------------------------------------------------

   subroutine test_evap_concentrates(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: E = -1.0e-2_wp
      integer, parameter :: NSTEP = 10
      real(wp) :: dh, expect_h, s_top0, t_top0, c_before, c_after
      integer :: n, i, j, k

      call setup_column(grid, ms, sf, at, SURF_FW_MASS)
      s_top0 = ms%tracers(ms%idx_salinity)%c(NGHOST + 1, NGHOST + 1, NZ)
      t_top0 = ms%tracers(ms%idx_temperature)%c(NGHOST + 1, NGHOST + 1, NZ)
      sf%evap = E
      sf%has_mass_flux = .true.
      call map_all(ms, sf, at)
      !$acc update device(sf%evap)
      do n = 1, NSTEP
         call one_thermo_step(grid, ms, sf, at)
      end do
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at)

      dh = E*DT/sf%rho0
      expect_h = H_LAYER + real(NSTEP, wp)*dh
      i = NGHOST + 2
      j = NGHOST + 2
      call check(error, abs(ms%h_layer(i, j, NZ) - expect_h) < 1.0e-14_wp*H_LAYER, "evaporation thins the top layer")
      if (allocated(error)) return
      call check(error, abs(ms%tracers(ms%idx_salinity)%c(i, j, NZ) - s_top0*H_LAYER/expect_h) < 1.0e-14_wp*s_top0, &
                 "evaporation concentrates salt: S = S0 h0/h")
      if (allocated(error)) return
      call check(error, ms%tracers(ms%idx_salinity)%c(i, j, NZ) > s_top0, "S rose")
      if (allocated(error)) return
      call check(error, ms%tracers(ms%idx_temperature)%c(i, j, NZ) == t_top0, "T keeps its exact bits")
      if (allocated(error)) return
      c_before = 0.0_wp
      c_after = 0.0_wp
      do k = 1, NZ
         c_before = c_before + (S0 - 0.1_wp*real(k, wp))*H_LAYER
         c_after = c_after + ms%tracers(ms%idx_salinity)%c(i, j, k)*ms%h_layer(i, j, k)
      end do
      call check(error, abs(c_after - c_before) < 1.0e-14_wp*c_before, "the column salt content is unchanged")
      if (allocated(error)) return
      call check(error, src_mass(ms) < 0.0_wp, "the tracked source is negative (mass left)")
      call sf%destroy(); call ms%destroy()
   end subroutine test_evap_concentrates

   ! -----------------------------------------------------------------
   ! 3. pass-down
   ! -----------------------------------------------------------------

   subroutine test_evap_cfl_passdown(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: H_TOP = 2.0e-3_wp
      real(wp), parameter :: REQ = 0.5_wp
      real(wp) :: e_flux, left, h2_expect, c_ini, c_fin, sum_before, sum_after
      integer :: i, j, k

      call setup_column(grid, ms, sf, at, SURF_FW_MASS)
      ms%h_layer(:, :, NZ) = H_TOP
      e_flux = -REQ*sf%rho0/DT
      sf%evap = e_flux
      sf%has_mass_flux = .true.
      call map_all(ms, sf, at)
      !$acc update device(ms%h_layer)
      !$acc update device(sf%evap)
      call one_thermo_step(grid, ms, sf, at)
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at)

      i = NGHOST + 3
      j = NGHOST + 2
      left = REQ - 0.8_wp*H_TOP
      h2_expect = H_LAYER - left
      call check(error, abs(ms%h_layer(i, j, NZ) - 0.2_wp*H_TOP) < 1.0e-14_wp*H_TOP, &
                 "the thin top layer ends at (1 - 0.8) h0")
      if (allocated(error)) return
      call check(error, abs(ms%h_layer(i, j, NZ - 1) - h2_expect) < 1.0e-13_wp*H_LAYER, &
                 "the remainder is drawn from the next live layer")
      if (allocated(error)) return
      call check(error, ms%h_layer(i, j, 1) == H_LAYER .and. ms%h_layer(i, j, 2) == H_LAYER, &
                 "deeper layers are not touched")
      if (allocated(error)) return
      sum_before = 3.0_wp*H_LAYER + H_TOP
      sum_after = sum(ms%h_layer(i, j, :))
      call check(error, abs((sum_before - sum_after) - REQ) < 1.0e-13_wp, "the column loses exactly the request")
      if (allocated(error)) return
      call check(error, all(ms%h_layer > 0.0_wp) .and. all(ms%h_layer(:, :, NZ) > H_VANISHED), &
                 "every layer stays positive (the top stays live)")
      if (allocated(error)) return
      call check(error, abs(ms%tracers(ms%idx_salinity)%c(i, j, NZ)*ms%h_layer(i, j, NZ) &
                            - (S0 - 0.1_wp*real(NZ, wp))*H_TOP) < 1.0e-14_wp*H_TOP*S0, &
                 "the thin layer keeps its salt content (it concentrates 5x)")
      if (allocated(error)) return
      c_ini = (S0 - 0.1_wp*real(NZ - 1, wp))*H_LAYER
      c_fin = ms%tracers(ms%idx_salinity)%c(i, j, NZ - 1)*ms%h_layer(i, j, NZ - 1)
      call check(error, abs(c_fin - c_ini) < 1.0e-13_wp*c_ini, "the layer drawn from keeps its salt content too")
      if (allocated(error)) return
      call check(error, ms%tracers(ms%idx_temperature)%c(i, j, NZ - 1) == T0 + 0.5_wp*real(NZ - 1, wp), &
                 "T of the layer drawn from keeps its exact bits")
      if (allocated(error)) return
      call check(error, all(ms%fw_passdown(:, :) == 1.0_wp), "every wet column flags exactly one pass-down layer")
      if (allocated(error)) return
      call check(error, ms%fw_n_passdown_step == int(NXP*NYP, int64), "the pass-down count is over OWNED columns only")
      if (allocated(error)) return
      call check(error, ms%fw_n_thin_step == 0_int64, "no column ran dry")
      if (allocated(error)) return
      call check(error, abs(src_mass(ms) + RHO_WATER*REQ*real(NXP*NYP, wp)*AREA) < 1.0e-13_wp*RHO_WATER*REQ*24.0_wp*AREA, &
                 "the tracked source books the APPLIED withdrawal")
      if (allocated(error)) return
      ! heat content booked at the layers' own temperatures
      sum_before = ms%heat_budget_surface(i, j, NZ) + ms%heat_budget_surface(i, j, NZ - 1)
      call check(error, abs(sum_before + 0.8_wp*H_TOP*(T0 + 0.5_wp*real(NZ, wp)) &
                            + left*(T0 + 0.5_wp*real(NZ - 1, wp))) < 1.0e-12_wp, &
                 "the heat content dh*T is booked layer by layer")
      call sf%destroy(); call ms%destroy()
   end subroutine test_evap_cfl_passdown

   ! -----------------------------------------------------------------
   ! 4. exhaustion
   ! -----------------------------------------------------------------

   subroutine test_column_exhaustion_is_fatal(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: H_THIN = 2.0e-3_wp
      integer :: n_thin, n_pd, nx, ny

      call setup_column(grid, ms, sf, at, SURF_FW_MASS)
      ms%h_layer = H_THIN
      sf%fw_net = -1.0_wp*sf%rho0/DT
      nx = grid%nx_total
      ny = grid%ny_total
      call map_all(ms, sf, at)
      !$acc update device(ms%h_layer, sf%fw_net)
      call fw_mass_column_impl(nx, ny, NZ, NGHOST, DT/sf%rho0, sf%fw_evap_cfl_limit, sf%fw_net, ms%wet_mask, &
                               ms%k_top, ms%h_layer, ms%tracers(ms%idx_temperature)%c, ms%heat_budget_surface, &
                               ms%fw_dh, ms%fw_passdown, n_thin, n_pd)
      !$acc update self(ms%h_layer, ms%fw_dh)
      call unmap_all(ms, sf, at)

      call check(error, n_thin == NXP*NYP, "every OWNED column is counted (ghost columns are not)")
      if (allocated(error)) return
      call check(error, fw_mass_thin_is_fatal(int(n_thin, int64)), "a non-zero count is fatal")
      if (allocated(error)) return
      call check(error,.not. fw_mass_thin_is_fatal(0_int64), "a zero count is not")
      if (allocated(error)) return
      call check(error, all(abs(ms%h_layer - 0.2_wp*H_THIN) < 1.0e-16_wp), &
                 "the clamped stack is left at (1 - f) h per layer, still live")
      if (allocated(error)) return
      call check(error, all(ms%h_layer > H_VANISHED), "no layer fell below the vanished marker")
      if (allocated(error)) return
      call check(error, abs(ms%fw_dh(NGHOST + 1, NGHOST + 1) + 0.8_wp*H_THIN*real(NZ, wp)) < 1.0e-15_wp, &
                 "the applied change is what the stack could give")
      call sf%destroy(); call ms%destroy()
   end subroutine test_column_exhaustion_is_fatal

   ! -----------------------------------------------------------------
   ! 4b. no accepted limit may lose dissolved content
   ! -----------------------------------------------------------------

   subroutine test_evap_keeps_salt(error)
      !! A thin top layer, an evaporation larger than it, a limit at the very top of the range (and
      !! at 1.0, which `validate_config` refuses but the kernel must survive): the layer must stay
      !! live, the pass-down must supply the rest, the column's salt is exact and its volume falls
      !! by exactly the request.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: H_TOP = 1.0e-3_wp
      real(wp), parameter :: REQ = 0.5_wp
      real(wp), parameter :: LIMITS(2) = [1.0_wp - 1.0e-9_wp, 1.0_wp]
      real(wp) :: salt0, salt1, vol0
      integer :: i, j, ic

      do ic = 1, size(LIMITS)
         call setup_column(grid, ms, sf, at, SURF_FW_MASS)
         ms%h_layer(:, :, NZ) = H_TOP
         sf%fw_evap_cfl_limit = LIMITS(ic)
         sf%evap = -REQ*sf%rho0/DT
         sf%has_mass_flux = .true.
         i = NGHOST + 3
         j = NGHOST + 2
         salt0 = sum(ms%h_layer(i, j, :)*ms%tracers(ms%idx_salinity)%c(i, j, :))
         vol0 = sum(ms%h_layer(i, j, :))
         call map_all(ms, sf, at)
         !$acc update device(ms%h_layer)
         !$acc update device(sf%evap)
         call one_thermo_step(grid, ms, sf, at)
         call pull_state(ms, sf)
         call unmap_all(ms, sf, at)

         salt1 = sum(ms%h_layer(i, j, :)*ms%tracers(ms%idx_salinity)%c(i, j, :))
         call check(error, abs(salt1 - salt0) < 1.0e-13_wp*salt0, "the column's salt content is exact")
         if (allocated(error)) return
         call check(error, abs((vol0 - sum(ms%h_layer(i, j, :))) - REQ) < 1.0e-13_wp, &
                    "the column loses exactly E*dt (the request)")
         if (allocated(error)) return
         call check(error, ms%h_layer(i, j, NZ) > H_VANISHED, "the thin top layer stays live")
         if (allocated(error)) return
         call check(error, ms%tracers(ms%idx_salinity)%c(i, j, NZ) < 1.0e4_wp, &
                    "its salinity stays physical (no emptied-layer blow-up)")
         if (allocated(error)) return
         call check(error, ms%fw_n_thin_step == 0_int64, "no column ran dry")
         if (allocated(error)) return
         call sf%destroy(); call ms%destroy()
      end do
   end subroutine test_evap_keeps_salt

   ! -----------------------------------------------------------------
   ! 5. passive tracers, age, pseudo-salt
   ! -----------------------------------------------------------------

   subroutine fill_pattern(sf, thin_h, ms)
      !! Rain where `mod(i+j,2)==0`, evaporation elsewhere (a different rate in every column),
      !! a thin top layer in every third column so the pass-down fires, one salt source.
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: thin_h
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j

      do j = 1, size(sf%lprec, 2)
         do i = 1, size(sf%lprec, 1)
            if (mod(i + j, 2) == 0) then
               sf%lprec(i, j) = 4.0e-3_wp*real(1 + mod(i, 3), wp)
            else
               sf%evap(i, j) = -6.0e-3_wp*real(1 + mod(j, 3), wp)
            end if
            if (mod(i, 3) == 0) ms%h_layer(i, j, NZ) = thin_h
         end do
      end do
      sf%salt_flux = 2.0e-6_wp
      sf%has_mass_flux = .true.
   end subroutine fill_pattern

   subroutine test_passive_and_age_diluted(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: P = 5.0e-3_wp
      real(wp) :: dh, worst_dye, worst_age, expect_h, c_dye0, c_age0
      integer :: i, j, k, idx_dye

      call setup_column(grid, ms, sf, at, SURF_FW_MASS, extra_tracers=.true.)
      idx_dye = 4
      c_dye0 = ms%tracers(idx_dye)%c(1, 1, 1)
      c_age0 = ms%tracers(ms%idx_age)%c(1, 1, 1)
      sf%lprec = P
      sf%has_mass_flux = .true.
      call map_all(ms, sf, at)
      !$acc update device(sf%lprec)
      call one_thermo_step(grid, ms, sf, at)
      call one_thermo_step(grid, ms, sf, at)
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at)

      dh = P*DT/sf%rho0
      expect_h = H_LAYER + 2.0_wp*dh
      worst_dye = 0.0_wp
      worst_age = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            worst_dye = max(worst_dye, abs(ms%tracers(idx_dye)%c(i, j, NZ) - c_dye0*H_LAYER/expect_h)/c_dye0)
            worst_age = max(worst_age, abs(ms%tracers(ms%idx_age)%c(i, j, NZ) - c_age0*H_LAYER/expect_h)/c_age0)
            do k = 1, NZ - 1
               call check(error, ms%tracers(idx_dye)%c(i, j, k) == c_dye0 .and. &
                          ms%tracers(ms%idx_age)%c(i, j, k) == c_age0, "lower layers keep their tracer exactly")
               if (allocated(error)) return
            end do
         end do
      end do
      call check(error, worst_dye < 1.0e-14_wp, "a passive tracer dilutes as h0/h (inflow carries none)")
      if (allocated(error)) return
      call check(error, worst_age < 1.0e-12_wp, "ideal age is diluted the same way (meltwater-like: new water)")
      if (allocated(error)) return
      call check(error, ms%tracers(idx_dye)%c(NGHOST + 2, NGHOST + 2, NZ) < c_dye0, "and it fell")
      call sf%destroy(); call ms%destroy()
   end subroutine test_passive_and_age_diluted

   subroutine test_pseudo_salt_bitwise(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      integer :: n, i, j, k, nbad

      call setup_column(grid, ms, sf, at, SURF_FW_MASS, extra_tracers=.true.)
      call fill_pattern(sf, 1.0e-3_wp, ms)
      call map_all(ms, sf, at)
      !$acc update device(ms%h_layer)
      !$acc update device(sf%lprec, sf%evap, sf%salt_flux)
      do n = 1, 6
         call one_thermo_step(grid, ms, sf, at)
      end do
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at)

      nbad = 0
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               if (ms%tracers(ms%idx_pseudo_salt)%c(i, j, k) /= ms%tracers(ms%idx_salinity)%c(i, j, k)) nbad = nbad + 1
            end do
         end do
      end do
      call check(error, nbad == 0, "pseudo-salt equals salinity bit for bit through rain, evaporation and pass-down")
      if (allocated(error)) return
      call check(error, ms%fw_n_passdown_total > 0_int64, "non-vacuous: the pass-down fired")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%c(:, :, NZ) - (S0 - 0.1_wp*real(NZ, wp)))) > 1.0e-3_wp, &
                 "non-vacuous: salinity moved")
      call sf%destroy(); call ms%destroy()
   end subroutine test_pseudo_salt_bitwise

   ! -----------------------------------------------------------------
   ! 6. B_0 identical, 7. second order
   ! -----------------------------------------------------------------

   subroutine assemble_pattern(mode, grid, ms, sf, at)
      integer, intent(in) :: mode
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(out) :: at(:, :)
      integer :: i, j, k

      call setup_column(grid, ms, sf, at, mode)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            sf%lprec(i, j) = 1.0e-5_wp*real(1 + mod(i, 3), wp)
            sf%fprec(i, j) = 2.0e-6_wp*real(mod(j, 2), wp)
            sf%evap(i, j) = -3.0e-6_wp*real(mod(i + j, 4), wp)
            sf%lrunoff(i, j) = 1.0e-6_wp*real(mod(i, 2), wp)
            do k = 1, NZ
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 30.0_wp + 0.5_wp*real(i, wp) + 0.25_wp*real(j, wp) &
                                                        + 0.1_wp*real(k, wp)
            end do
         end do
      end do
      sf%salt_flux = 1.0e-7_wp
      sf%salt_ice = -2.0e-7_wp
      sf%salt_cavity = 3.0e-7_wp
      sf%q_sw = 120.0_wp
      sf%q_lat = -40.0_wp
   end subroutine assemble_pattern

   subroutine test_b0_identical(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid, grid_m
      type(multilayer_state_t) :: ms_v, ms_m
      type(ocean_surface_flux_t) :: sf_v, sf_m
      real(wp), allocatable :: at_v(:, :), at_m(:, :), qs_v(:, :), qh_v(:, :)
      real(wp) :: worst, scale, fnet
      integer :: i, j, ng

      call assemble_pattern(SURF_FW_VIRTUAL, grid, ms_v, sf_v, at_v)
      call map_all(ms_v, sf_v, at_v)
      !$acc update device(sf_v%lprec, sf_v%fprec, sf_v%evap, sf_v%lrunoff, sf_v%salt_flux, sf_v%salt_ice)
      !$acc update device(sf_v%salt_cavity, sf_v%q_sw, sf_v%q_lat, ms_v%tracers(ms_v%idx_salinity)%c)
      call ocean_surface_flux_fw_virtual_step(grid, sf_v, ms_v)
      call ocean_surface_flux_assemble(grid, sf_v, ms_v)
      !$acc update self(sf_v%Q_salt, sf_v%Q_heat, sf_v%salt_fw, sf_v%fw_net, sf_v%Q_salt_tracer)
      allocate (qs_v, source=sf_v%Q_salt)
      allocate (qh_v, source=sf_v%Q_heat)
      call unmap_all(ms_v, sf_v, at_v)

      call assemble_pattern(SURF_FW_MASS, grid_m, ms_m, sf_m, at_m)
      call map_all(ms_m, sf_m, at_m)
      !$acc update device(sf_m%lprec, sf_m%fprec, sf_m%evap, sf_m%lrunoff, sf_m%salt_flux, sf_m%salt_ice)
      !$acc update device(sf_m%salt_cavity, sf_m%q_sw, sf_m%q_lat, ms_m%tracers(ms_m%idx_salinity)%c)
      call ocean_surface_flux_fw_virtual_step(grid_m, sf_m, ms_m)
      call ocean_surface_flux_assemble(grid_m, sf_m, ms_m)
      call pull_state(ms_m, sf_m)
      call unmap_all(ms_m, sf_m, at_m)

      ng = grid%nghost
      call check(error, all(sf_m%Q_salt == qs_v), "B_0 input: Q_salt is BITWISE the same under virtual and mass")
      if (allocated(error)) return
      call check(error, all(sf_m%Q_heat == qh_v), "and so is Q_heat (nothing owes mass enthalpy)")
      if (allocated(error)) return
      worst = 0.0_wp
      scale = maxval(abs(sf_m%Q_salt))
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            worst = max(worst, abs(sf_m%Q_salt_tracer(i, j) - (sf_m%Q_salt(i, j) - sf_m%salt_fw(i, j))))
         end do
      end do
      call check(error, worst < 1.0e-14_wp*scale, "Q_salt_tracer = Q_salt - salt_fw")
      if (allocated(error)) return
      call check(error, maxval(abs(sf_m%salt_fw)) > 0.0_wp, "non-vacuous: the virtual term is non-zero")
      if (allocated(error)) return
      i = ng + 2
      j = ng + 2
      fnet = (1.0e-5_wp*real(1 + mod(i, 3), wp) + 2.0e-6_wp*real(mod(j, 2), wp)) &
             - 3.0e-6_wp*real(mod(i + j, 4), wp) + 1.0e-6_wp*real(mod(i, 2), wp)
      call check(error, abs(sf_m%fw_net(i, j) - fnet) < 1.0e-18_wp, "fw_net is the net flux of the components")
      if (allocated(error)) return
      call check(error, all(sf_v%fw_net == 0.0_wp) .and. all(sf_v%Q_salt_tracer == 0.0_wp), &
                 "virtual leaves fw_net and Q_salt_tracer alone")
      call sf_v%destroy(); call ms_v%destroy(); call sf_m%destroy(); call ms_m%destroy()
   end subroutine test_b0_identical

   subroutine test_second_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :)
      real(wp), parameter :: EPS = 0.05_wp
      real(wp) :: p_rain, s_v, s_m, diff, expect, s_start
      integer :: pass

      do pass = 1, 2
         call setup_column(grid, ms, sf, at, merge(SURF_FW_VIRTUAL, SURF_FW_MASS, pass == 1))
         s_start = ms%tracers(ms%idx_salinity)%c(NGHOST + 2, NGHOST + 2, NZ)
         p_rain = EPS*sf%rho0*H_LAYER/DT
         sf%lprec = p_rain
         sf%has_mass_flux = .true.
         call map_all(ms, sf, at)
         !$acc update device(sf%lprec)
         call one_thermo_step(grid, ms, sf, at)
         call pull_state(ms, sf)
         call unmap_all(ms, sf, at)
         if (pass == 1) s_v = ms%tracers(ms%idx_salinity)%c(NGHOST + 2, NGHOST + 2, NZ)
         if (pass == 2) s_m = ms%tracers(ms%idx_salinity)%c(NGHOST + 2, NGHOST + 2, NZ)
         call sf%destroy(); call ms%destroy()
      end do
      call check(error, abs(s_v - s_start*(1.0_wp - EPS)) < 1.0e-13_wp, "virtual: S0 (1 - eps)")
      if (allocated(error)) return
      call check(error, abs(s_m - s_start/(1.0_wp + EPS)) < 1.0e-13_wp, "mass: S0 /(1 + eps)")
      if (allocated(error)) return
      diff = s_m - s_v
      expect = s_start*EPS*EPS/(1.0_wp + EPS)
      call check(error, abs(diff - expect) < 1.0e-13_wp, "dS_mass - dS_virtual = S0 eps^2/(1+eps), the cavity identity")
      if (allocated(error)) return
      call check(error, diff > 1.0e-2_wp, "non-vacuous: the two forms differ at second order")
   end subroutine test_second_order

   ! -----------------------------------------------------------------
   ! 12. none / virtual bit-identity
   ! -----------------------------------------------------------------

   subroutine test_none_virtual_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: at(:, :), garbage_net(:, :), garbage_tr(:, :)
      real(wp) :: massin, massout, qh, qs, worst_h, worst_s, sst, s_top
      integer :: i, j, mode, mm
      integer, parameter :: MODES(2) = [SURF_FW_NONE, SURF_FW_VIRTUAL]

      do mm = 1, 2
         mode = MODES(mm)
         call setup_column(grid, ms, sf, at, mode)
         sst = T0 + 0.5_wp*real(NZ, wp)
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               sf%lprec(i, j) = 1.0e-5_wp*real(i, wp)
               sf%fprec(i, j) = 2.0e-6_wp*real(j, wp)
               sf%evap(i, j) = -3.0e-6_wp*real(i + j, wp)
               sf%heat_content_lprec(i, j) = 4000.0_wp*sst*sf%lprec(i, j)
               sf%heat_content_fprec(i, j) = 4000.0_wp*sst*sf%fprec(i, j)
               sf%heat_content_vprec(i, j) = 1.5_wp
               sf%heat_content_seaice_melt(i, j) = 0.5_wp
               sf%q_lat(i, j) = -50.0_wp
               sf%q_sw(i, j) = 120.0_wp + real(i, wp)
            end do
         end do
         sf%salt_flux = 1.0e-7_wp
         sf%salt_ice = -2.0e-7_wp
         sf%salt_cavity = 3.0e-7_wp
         sf%fw_net = 777.0_wp
         sf%Q_salt_tracer = 888.0_wp
         allocate (garbage_net, source=sf%fw_net)
         allocate (garbage_tr, source=sf%Q_salt_tracer)
         call map_all(ms, sf, at)
         !$acc update device(sf%lprec, sf%fprec, sf%evap, sf%heat_content_lprec, sf%heat_content_fprec)
         !$acc update device(sf%heat_content_vprec, sf%heat_content_seaice_melt, sf%fw_net, sf%Q_salt_tracer)
         !$acc update device(sf%q_lat, sf%q_sw, sf%salt_flux, sf%salt_ice, sf%salt_cavity)
         call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
         call ocean_surface_flux_assemble(grid, sf, ms)
         !$acc update self(sf%salt_fw, sf%Q_salt, sf%Q_heat, sf%fw_net, sf%Q_salt_tracer)
         call unmap_all(ms, sf, at)

         call check(error, all(sf%fw_net == garbage_net) .and. all(sf%Q_salt_tracer == garbage_tr), &
                    "none/virtual: the mass outputs are never written")
         if (allocated(error)) return
         worst_h = 0.0_wp
         worst_s = 0.0_wp
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               if (mode == SURF_FW_NONE) then
                  massin = sf%heat_content_lprec(i, j) + sf%heat_content_fprec(i, j) + sf%heat_content_vprec(i, j) &
                           + sf%heat_content_lrunoff(i, j) + sf%heat_content_frunoff(i, j) &
                           + sf%heat_content_seaice_melt(i, j)
                  massout = sf%cp*sst*sf%evap(i, j)
                  qs = 0.0_wp + (1.0e-7_wp - 2.0e-7_wp) + 3.0e-7_wp
               else
                  massin = sf%heat_content_vprec(i, j) + sf%heat_content_seaice_melt(i, j)
                  massout = 0.0_wp
                  s_top = S0 - 0.1_wp*real(NZ, wp)
                  qs = 0.0_wp + (1.0e-7_wp - 2.0e-7_wp) + sf%salt_fw(i, j) + 3.0e-7_wp
                  if (i > NGHOST .and. i <= grid%nx_total - NGHOST .and. j > NGHOST .and. &
                      j <= grid%ny_total - NGHOST) then
                     call check(error, abs(sf%salt_fw(i, j) + s_top*(sf%lprec(i, j) + sf%fprec(i, j) + sf%evap(i, j))) &
                                < 1.0e-15_wp, "virtual: salt_fw = -S_surf F")
                     if (allocated(error)) return
                  end if
               end if
               qh = (0.0_wp + (sf%q_sw(i, j) + 0.0_wp) + 0.0_wp + sf%q_lat(i, j) + 0.0_wp &
                     + (0.0_wp + 0.0_wp) + 0.0_wp + massin + massout)
               worst_h = max(worst_h, abs(sf%Q_heat(i, j) - qh)/(abs(qh) + 1.0_wp))
               worst_s = max(worst_s, abs(sf%Q_salt(i, j) - qs)/abs(qs))
            end do
         end do
         call check(error, worst_s < 1.0e-15_wp, "Q_salt is the legacy sum")
         if (allocated(error)) return
         call check(error, worst_h < 1.0e-15_wp, "Q_heat is the legacy sum (mass enthalpy per the mode)")
         if (allocated(error)) return
         deallocate (garbage_net, garbage_tr)
         call sf%destroy(); call ms%destroy()
      end do
   end subroutine test_none_virtual_bitident

   ! -----------------------------------------------------------------
   ! engine-level scaffolding
   ! -----------------------------------------------------------------

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function mass_nml(scheme, therm_ratio, deferred, vcoord, t_top, extra, adv_ratio, forcing_keys, components, &
                     diag_dir, freshwater) result(nml)
      !! A closed 12 x 8 x 4 basin with the real-mass fresh-water mode on.  Every group appears ONCE
      !! (the schema parser refuses a duplicate), so the knobs the refusal cases vary are arguments.
      character(len=*), intent(in) :: scheme, vcoord, extra
      integer, intent(in) :: therm_ratio
      logical, intent(in) :: deferred
      real(wp), intent(in) :: t_top
      integer, intent(in), optional :: adv_ratio
      character(len=*), intent(in), optional :: forcing_keys
      logical, intent(in), optional :: components
      character(len=*), intent(in), optional :: diag_dir, freshwater
      character(len=:), allocatable :: nml, vmix_keys, fw_s, diag_s, out_s
      character(len=32) :: ratio_s, tt_s, adv_s
      logical :: comps

      write (ratio_s, '(I0)') therm_ratio
      write (tt_s, '(F8.3)') t_top
      comps = .true.
      if (present(components)) comps = components
      fw_s = "mass"
      if (present(freshwater)) fw_s = freshwater
      diag_s = "&ocean_diag_nml enabled = .false. /"//NL
      out_s = "&output_nml output_to_file = .false. /"//NL
      if (present(diag_dir)) then
         diag_s = "&ocean_diag_nml enabled = .true., filename = 'd.nc', dt_out = 1200.0 /"//NL
         out_s = "&output_nml output_to_file = .true., output_dir = '"//diag_dir//"' /"//NL
      end if
      vmix_keys = ""
      if (therm_ratio > 1) vmix_keys = vmix_keys//" dt_therm_ratio = "//trim(adjustl(ratio_s))
      if (deferred) vmix_keys = vmix_keys//", diabatic_first = .false."
      if (deferred .or. present(adv_ratio)) then
         if (present(adv_ratio)) then
            write (adv_s, '(I0)') adv_ratio
         else
            write (adv_s, '(I0)') therm_ratio
         end if
         vmix_keys = vmix_keys//", dt_tracer_advect_ratio = "//trim(adjustl(adv_s))
      end if
      if (len(vmix_keys) > 0) then
         if (vmix_keys(1:1) == ",") vmix_keys = vmix_keys(2:)
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 12, ny = 8, nghost = 3, dx = 4000.0, dy = 4000.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 400.0, wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_bottom = 6.0, T_init_surface = "// &
            trim(adjustl(tt_s))//" /"//NL// &
            "&ocean_eos_nml eos = 'linear' /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&ocean_bt_nml split_scheme = '"//scheme//"', auto_n_inner = .true. /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = "//merge(".true. ", ".false.", comps)// &
            ", freshwater = '"//fw_s//"'"
      if (present(forcing_keys)) nml = nml//forcing_keys
      nml = nml//" /"//NL// &
            "&ocean_hvisc_nml nu_h = 100.0 /"//NL//diag_s//out_s
      if (trim(vcoord) == "zstar_list") then
         nml = nml//"&vcoord_nml vcoord_type = 'zstar', z_fixed_profile = 'list', "// &
               "z_fixed_dz = 20.0, 60.0, 120.0, 200.0 /"//NL
      else
         nml = nml//"&vcoord_nml vcoord_type = '"//trim(vcoord)//"' /"//NL
      end if
      if (len(vmix_keys) > 0) nml = nml//"&ocean_vmix_nml "//vmix_keys//" /"//NL
      nml = nml//extra
   end function mass_nml

   subroutine start_engine(nml, engine, cfg, ok)
      character(len=*), intent(in) :: nml
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      integer :: ierr
      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= 0) return
      call validate_config(cfg, ierr)
      if (ierr /= 0) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine start_engine

   subroutine set_forcing(engine, mode)
      !! `mode = 1`: uniform rain.  `mode = 2`: a rain/evaporation pattern, different in every column.
      type(ocean_engine_t), intent(inout) :: engine
      integer, intent(in) :: mode
      integer :: i, j

      associate (sf => engine%state%surface_flux)
         if (mode == 1) then
            sf%lprec = 5.0e-3_wp
         else
            do j = 1, size(sf%lprec, 2)
               do i = 1, size(sf%lprec, 1)
                  if (mod(i + j, 2) == 0) then
                     sf%lprec(i, j) = 4.0e-3_wp*real(1 + mod(i, 3), wp)
                  else
                     sf%evap(i, j) = -3.0e-3_wp*real(1 + mod(j, 3), wp)
                  end if
               end do
            end do
         end if
         sf%has_mass_flux = .true.
         !$acc update device(sf%lprec, sf%evap)
      end associate
   end subroutine set_forcing

   subroutine advance(engine, cfg, t, ierr)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t
      integer, intent(out) :: ierr
      real(wp), parameter :: DTE = 300.0_wp
      call engine_step(engine, DTE, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_ice(engine, cfg, DTE, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_finalize(engine, DTE, t, ierr=ierr)
   end subroutine advance

   subroutine totals(engine, total_mass, total_heat, src_heat, out_heat, mass_out, mass_src, total_salt, src_salt, out_salt)
      !! Exactly what `ocean_console_stats_report` forms, through the production reducers.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: total_mass, total_heat, src_heat, out_heat, mass_out, mass_src
      real(wp), intent(out) :: total_salt, src_salt, out_salt
      real(wp) :: bud_w
      type(efp_t) :: e

      associate (ms => engine%state%multilayer, mt => engine%state%metrics, ng => engine%grid%nghost)
         !$acc update self(ms%h_layer, ms%heat_budget_surface, ms%salt_budget_surface)
         !$acc update self(ms%heat_budget_horiz_adv, ms%heat_budget_hdiff, ms%salt_budget_horiz_adv, ms%salt_budget_hdiff)
         associate (cT => ms%tracers(ms%idx_temperature)%c, cS => ms%tracers(ms%idx_salinity)%c)
            !$acc update self(cT, cS)
         end associate
         bud_w = ocean_budget_stage_weight(engine%state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR)
         total_mass = compute_total_h(ms%h_layer, mt%areaT, ng)*RHO_WATER
         total_heat = compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, mt%areaT, ng)*RHO_WATER
         total_salt = compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, mt%areaT, ng)*RHO_WATER
         src_heat = ocean_budget_src( &
                    ocean_heat_src_sum( &
                    compute_total_content(ms%heat_budget_surface, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_geothermal, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_heat = ocean_budget_out( &
                    compute_total_content(ms%heat_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)
         src_salt = ocean_budget_src( &
                    ocean_salt_src_sum( &
                    compute_total_content(ms%salt_budget_surface, mt%areaT, ng), &
                    compute_total_content(ms%salt_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_salt = ocean_budget_out( &
                    compute_total_content(ms%salt_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_content(ms%salt_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)
         mass_out = ms%mass_out
         e%v = ms%mass_src_efp
         e%poison = ms%mass_src_efp_poison
         mass_src = ms%mass_src + real(efp_to_real(e), wp)
      end associate
   end subroutine totals

   pure function residual(total, ref, out_term, src_term) result(rel)
      real(wp), intent(in) :: total, ref, out_term, src_term
      real(wp) :: rel
      rel = ((total - ref) + out_term - src_term)/ref
   end function residual

   subroutine budget_case(error, scheme, ratio, deferred, label, check_heat)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme, label
      integer, intent(in) :: ratio
      logical, intent(in) :: deferred, check_heat
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr
      real(wp) :: t, m0, h0_, s0_, mt_, ht_, st_, src_h, out_h, m_out, m_src, src_s, out_s
      real(wp) :: worst_m, worst_h, worst_s, src_h_last, src_m_last

      call start_engine(mass_nml(scheme, ratio, deferred, "zstar", 6.0_wp, ""), engine, cfg, ok)
      call check(error, ok, label//": the engine starts")
      if (allocated(error)) return
      gate: block
         call set_forcing(engine, 2)
         call totals(engine, m0, h0_, src_h, out_h, m_out, m_src, s0_, src_s, out_s)
         worst_m = 0.0_wp
         worst_h = 0.0_wp
         worst_s = 0.0_wp
         t = 0.0_wp
         do n = 1, 12
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": the step succeeds")
            if (allocated(error)) exit gate
            t = t + 300.0_wp
            call totals(engine, mt_, ht_, src_h, out_h, m_out, m_src, st_, src_s, out_s)
            worst_m = max(worst_m, abs(residual(mt_, m0, m_out, m_src)))
            worst_h = max(worst_h, abs(residual(ht_, h0_, out_h, src_h)))
            worst_s = max(worst_s, abs(residual(st_, s0_, out_s, src_s)))
         end do
         src_h_last = src_h
         src_m_last = m_src
         write (output_unit, "(a,a,3es11.3,2es13.5)") "  fw mass ", label, worst_m, worst_h, worst_s, m_src, src_h
         call check(error, abs(src_m_last) > 1.0e-6_wp*m0, label//": non-vacuous: the mass source is non-zero")
         if (allocated(error)) exit gate
         call check(error, abs(mt_ - m0) > 1.0e-7_wp*m0, label//": non-vacuous: the domain mass actually changed")
         if (allocated(error)) exit gate
         call check(error, worst_m < 1.0e-13_wp, label//": mass residual (M - M0) + out - src")
         if (allocated(error)) exit gate
         call check(error, worst_s < 1.0e-13_wp, label//": salt residual (the content never moves)")
         if (allocated(error)) exit gate
         if (check_heat) then
            call check(error, abs(src_h_last) > 1.0e-6_wp*abs(h0_), label//": non-vacuous: heat content was booked")
            if (allocated(error)) exit gate
            call check(error, worst_h < 1.0e-12_wp, label//": heat residual (dh*T booked layer by layer)")
            if (allocated(error)) exit gate
         end if
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine budget_case

   subroutine test_mass_budget_closes(error)
      type(error_type), allocatable, intent(out) :: error
      call budget_case(error, "ssp_rk2", 1, .false., "ssp_rk2", .false.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 1, .false., "pred_corr", .false.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 2, .false., "ratio2", .false.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 2, .true., "deferred", .false.)
   end subroutine test_mass_budget_closes

   subroutine test_heat_budget_closes(error)
      type(error_type), allocatable, intent(out) :: error
      call budget_case(error, "ssp_rk2", 1, .false., "ssp_rk2", .true.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 1, .false., "pred_corr", .true.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 2, .false., "ratio2", .true.)
      if (allocated(error)) return
      call budget_case(error, "pred_corr", 2, .true., "deferred", .true.)
   end subroutine test_heat_budget_closes

   subroutine absolute_case(error, scheme, ratio, deferred, label, n_apply)
      !! ABSOLUTE time factor: the budget gates close against the APPLIED dh, so a wrong factor
      !! (outer dt for thermo dt, a doubled call) would still close.  Uniform P on a closed basin:
      !! `M - M0` must equal `n_apply * dt_thermo * F * A * RHO_WATER/rho_0`, `n_apply` an INTEGER
      !! count of thermo intervals, to round-off.  In the IN-STAGE orders (ssp_rk2, pred_corr) the
      !! count is one short of the intervals run: the stage spends the snapshot of the PREVIOUS
      !! assembly, which is zero in the first interval after a cold start, so every later interval
      !! spends exactly one snapshot over the whole interval.  The DEFERRED order assembles and
      !! spends in the same block: no lag, all six intervals.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme, label
      integer, intent(in) :: ratio, n_apply
      logical, intent(in) :: deferred
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr
      real(wp) :: t, m0, mt, dummy(8), expect, area_tot, rho0
      integer, parameter :: NSTEPS = 12
      real(wp), parameter :: P = 5.0e-3_wp, DTE = 300.0_wp

      call start_engine(mass_nml(scheme, ratio, deferred, "zstar", 6.0_wp, ""), engine, cfg, ok)
      call check(error, ok, label//": the engine starts")
      if (allocated(error)) return
      gate: block
         call set_forcing(engine, 1)
         call totals(engine, m0, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8))
         t = 0.0_wp
         do n = 1, NSTEPS
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": the step succeeds")
            if (allocated(error)) exit gate
            t = t + DTE
         end do
         call totals(engine, mt, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8))
         rho0 = engine%state%surface_flux%rho0
         area_tot = 12.0_wp*8.0_wp*4000.0_wp*4000.0_wp
         expect = real(n_apply, wp)*DTE*real(ratio, wp)*P*area_tot*RHO_WATER/rho0
         write (output_unit, "(a,a,2es24.14,f10.5,es11.3)") "  fw abs ", label, mt - m0, expect, &
            (mt - m0)/(P*area_tot*RHO_WATER/rho0*DTE), (mt - m0 - expect)/m0
         ! the tolerance is the round-off of the TOTAL (M ~ 6e14 kg summed over 96 cells), ~1e-14 M,
         ! not of the change: still 1e8 times finer than the smallest wrong factor (one interval in n)
         call check(error, abs((mt - m0) - expect) < 1.0e-13_wp*m0, &
                    label//": M - M0 is exactly the applied intervals times F*A*dt")
         if (allocated(error)) exit gate
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine absolute_case

   subroutine test_absolute_time_factor(error)
      type(error_type), allocatable, intent(out) :: error
      call absolute_case(error, "ssp_rk2", 1, .false., "ssp_rk2", 11)
      if (allocated(error)) return
      call absolute_case(error, "pred_corr", 2, .false., "ratio2", 5)
      if (allocated(error)) return
      call absolute_case(error, "pred_corr", 2, .true., "deferred", 6)
   end subroutine test_absolute_time_factor

   subroutine no_flow_case(error, scheme, label)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme, label
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr, i, j, nsteps
      real(wp) :: t, umax, vmax, dh, eta_expect, eta_min, eta_max, h400
      real(wp), parameter :: P = 5.0e-3_wp

      call start_engine(mass_nml(scheme, 1, .false., "sigma", 6.0_wp, ""), engine, cfg, ok)
      call check(error, ok, label//": the engine starts")
      if (allocated(error)) return
      gate: block
         call set_forcing(engine, 1)
         nsteps = 12
         t = 0.0_wp
         do n = 1, nsteps
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": the step succeeds")
            if (allocated(error)) exit gate
            t = t + 300.0_wp
         end do
         associate (ms => engine%state%multilayer, ng => engine%grid%nghost)
            !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
            umax = maxval(abs(ms%u_face_x_layer))
            vmax = maxval(abs(ms%v_face_y_layer))
            eta_min = huge(1.0_wp)
            eta_max = -huge(1.0_wp)
            do j = ng + 1, engine%grid%ny_total - ng
               do i = ng + 1, engine%grid%nx_total - ng
                  h400 = sum(ms%h_layer(i, j, :))
                  eta_min = min(eta_min, h400)
                  eta_max = max(eta_max, h400)
               end do
            end do
         end associate
         ! the first step runs on the configure-time seed (no snapshot yet): nsteps - 1 applications
         dh = P*300.0_wp/engine%state%surface_flux%rho0
         eta_expect = 400.0_wp + real(nsteps - 1, wp)*dh
         write (output_unit, "(a,a,3es12.4)") "  fw no-flow ", label, umax, vmax, eta_min - 400.0_wp
         call check(error, eta_min - 400.0_wp > 0.5_wp*real(nsteps - 1, wp)*dh, label//": the surface rose (non-vacuous)")
         if (allocated(error)) exit gate
         call check(error, abs(eta_min - eta_expect) < 1.0e-9_wp*400.0_wp .and. abs(eta_max - eta_expect) < 1.0e-9_wp*400.0_wp, &
                    label//": eta rises by P t/rho0 in every column")
         if (allocated(error)) exit gate
         call check(error, umax <= 1.0e-12_wp .and. vmax <= 1.0e-12_wp, label//": no current")
         if (allocated(error)) exit gate
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine no_flow_case

   subroutine test_no_flow(error)
      type(error_type), allocatable, intent(out) :: error
      call no_flow_case(error, "ssp_rk2", "ssp_rk2")
      if (allocated(error)) return
      call no_flow_case(error, "pred_corr", "pred_corr")
   end subroutine test_no_flow

   subroutine test_zstar_interfaces(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr, i, j, k
      real(wp) :: t, worst, sum_h, expect
      real(wp), parameter :: DZ(4) = [20.0_wp, 60.0_wp, 120.0_wp, 200.0_wp]

      ! DZ is surface first; the layer index k = 1 is the bed
      call start_engine(mass_nml("pred_corr", 1, .false., "zstar_list", 6.0_wp, ""), engine, cfg, ok)
      call check(error, ok, "the engine starts")
      if (allocated(error)) return
      gate: block
         call set_forcing(engine, 1)
         engine%state%surface_flux%lprec = 4.0e-2_wp
         !$acc update device(engine%state%surface_flux%lprec)
         t = 0.0_wp
         do n = 1, 6
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "the step succeeds")
            if (allocated(error)) exit gate
            t = t + 300.0_wp
         end do
         associate (ms => engine%state%multilayer, ng => engine%grid%nghost)
            !$acc update self(ms%h_layer)
            worst = 0.0_wp
            do j = ng + 1, engine%grid%ny_total - ng
               do i = ng + 1, engine%grid%nx_total - ng
                  sum_h = sum(ms%h_layer(i, j, :))
                  do k = 1, 4
                     expect = DZ(5 - k)*sum_h/400.0_wp
                     worst = max(worst, abs(ms%h_layer(i, j, k) - expect)/expect)
                  end do
               end do
            end do
          write (output_unit, "(a,2es13.5)") "  fw zstar: worst rel dev, eta:", worst, sum(ms%h_layer(ng + 2, ng + 2, :)) - 400.0_wp
            call check(error, sum(ms%h_layer(ng + 2, ng + 2, :)) - 400.0_wp > 1.0e-3_wp, "non-vacuous: eta rose")
            if (allocated(error)) exit gate
            call check(error, worst < 1.0e-9_wp, &
                       "interfaces = the z* target of H and eta = sum(h) - H (a stale bt_eta builds it on H + dh)")
            if (allocated(error)) exit gate
         end associate
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine test_zstar_interfaces

   ! -----------------------------------------------------------------
   ! diagnostics
   ! -----------------------------------------------------------------

   subroutine test_diagnostics(error)
      !! `net_fw`, `evap`, `lprec`, `fprec`, `fw_passdown` are registered under `freshwater = "mass"` only, and
      !! the registered fill of `net_fw` returns the assembled snapshot.
      type(error_type), allocatable, intent(out) :: error
      call diagnostics_pass(error, .true.)
      if (allocated(error)) return
      call diagnostics_pass(error, .false.)
   end subroutine test_diagnostics

   subroutine diagnostics_pass(error, mass)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: mass
#ifndef RDB_NO_NETCDF
      character(len=*), parameter :: NAMES(5) = [character(len=16) :: "net_fw", "evap", "lprec", "fprec", "fw_passdown"]
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr, k, kk, ip, ngot
      real(wp) :: t
      real(wp), allocatable :: buf(:, :, :)

      call start_engine(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, "", diag_dir="tmp_fwm_diag", &
                                 freshwater=merge("mass", "none", mass)), engine, cfg, ok)
      call check(error, ok, "the engine with a diagnostic stream starts")
      if (allocated(error)) return
      call set_forcing(engine, 1)
      t = 0.0_wp
      do n = 1, 3
         call advance(engine, cfg, t, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "the step succeeds")
         if (allocated(error)) exit
         t = t + 300.0_wp
      end do
      ngot = 0
      kk = 0
      do k = 1, engine%state%diag%nvars
         do ip = 1, size(NAMES)
            if (trim(engine%state%diag%vars(k)%name) == trim(NAMES(ip))) then
               ngot = ngot + 1
               if (ip == 1) kk = k
            end if
         end do
      end do
      if (.not. allocated(error)) then
         if (mass) then
            call check(error, ngot == size(NAMES), "mass: all five fresh-water diagnostics are registered")
            if (.not. allocated(error)) then
               allocate (buf(size(engine%state%metrics%areaT, 1), size(engine%state%metrics%areaT, 2), 1), &
                         source=-777.0_wp)
               !$acc enter data copyin(buf)
               call engine%state%diag%vars(kk)%fill(engine%state, buf)
               !$acc update self(buf)
               !$acc exit data delete(buf)
               call check(error, abs(buf(engine%grid%nghost + 2, engine%grid%nghost + 2, 1) - 5.0e-3_wp) < 1.0e-15_wp, &
                          "net_fw is the assembled snapshot (uniform rain 5e-3 kg m-2 s-1)")
            end if
         else
            call check(error, ngot == 0, "none: no fresh-water diagnostic is registered")
         end if
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
      call execute_command_line("rm -rf tmp_fwm_diag")
#else
      call check(error, mass .or. .true., "diagnostics need NetCDF")
#endif
   end subroutine diagnostics_pass

   ! -----------------------------------------------------------------
   ! 13. configure-time envelope
   ! -----------------------------------------------------------------

   function nml_ierr(nml) result(ierr)
      character(len=*), intent(in) :: nml
      integer :: ierr
      type(config_t) :: cfg
      integer :: ierr_read

      call read_config_from_string(nml, cfg, ierr=ierr_read)
      if (ierr_read /= OCEAN_STATUS_OK) then
         ierr = -ierr_read - 1000
         return
      end if
      call validate_config(cfg, ierr)
   end function nml_ierr

   subroutine test_refusals(error)
      !! Every refusal differs from an accepted baseline by one knob and must be a
      !! `validate_config` verdict (a parse failure would be a different code).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: REFUSED = OCEAN_STATUS_ERR_CONFIG_VALIDATE
      character(len=*), parameter :: ICE = "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2 /"//NL

      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, "")) == OCEAN_STATUS_OK, &
                 "accepted: pred_corr, zstar, ratio 1")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("ssp_rk2", 1, .false., "sigma", 6.0_wp, "")) == OCEAN_STATUS_OK, &
                 "accepted: ssp_rk2, sigma")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 2, .true., "zstar", 6.0_wp, "")) == OCEAN_STATUS_OK, &
                 "accepted: the deferred path (dt_therm_ratio 2, diabatic_first off, windowed advection)")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "z_fixed", 6.0_wp, "")) == OCEAN_STATUS_OK, &
                 "accepted: z_fixed without a cavity")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("ssp_rk2", 1, .false., "eulerian_z", 6.0_wp, "")) == REFUSED, &
                 "refused: eulerian_z")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("ssp_rk2", 1, .false., "lagrangian", 6.0_wp, "")) == REFUSED, &
                 "refused: lagrangian")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("ssp_rk2", 1, .false., "sigma", 6.0_wp, &
                                          "&ocean_wetdry_nml enable = .true. /"//NL)) == REFUSED, &
                 "refused: wet/dry")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, "", adv_ratio=2)) == REFUSED, &
                 "refused: windowed tracer advection with the in-stage order (ratio 1)")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 2, .false., "zstar", 6.0_wp, "", adv_ratio=2)) == REFUSED, &
                 "refused: windowed advection with diabatic_first on")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, ICE)) == REFUSED, &
                 "refused: sea ice (the hand-off has no mass slot yet)")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, "", components=.false.)) == REFUSED, &
                 "refused: without the component set")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, "", &
                                          forcing_keys=", fw_evap_cfl_limit = 1.5")) == REFUSED, &
                 "refused: fw_evap_cfl_limit outside (0, 1)")
      if (allocated(error)) return
      call check(error, nml_ierr(mass_nml("pred_corr", 1, .false., "zstar", 6.0_wp, &
                                          "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
                                          "draft_depth = 100.0, draft_slope = -1.0e-3, draft_x0 = 0.0, "// &
                                          "h_min_cavity = 50.0 /"//NL)) == REFUSED, &
                 "refused: the ice-shelf cavity")
   end subroutine test_refusals

end module test_ocean_surface_fw_mass
