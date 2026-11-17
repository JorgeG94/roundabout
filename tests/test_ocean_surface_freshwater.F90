!! Surface fresh-water forcing, `&ocean_forcing_nml freshwater = "virtual"`: the net open-ocean
!! fresh-water mass flux `F = lprec + fprec + evap + lrunoff + frunoff` becomes the virtual salt
!! flux `salt_fw = -S_surf*F` (`ocean_surface_flux_fw_virtual_step`), summed by the assembler.
!!
!! Cases:
!!   1. `rain_resting_column`: uniform rain on a resting stratified column; the surface salinity
!!      follows the closed form `S0*(1 - a)^n` (`a = dt*P/(rho0*h)`), whose continuous limit is the
!!      exponential `S0*exp(-a n)`; temperature and column mass do not move; with the knob off the
!!      same forcing warms the layer by the analytic `dt*T*P/(rho0*h)` (the control that makes the
!!      heat check able to fail).
!!   2. `components_sign_composition`: every cell carries different components and a different
!!      salinity profile; each component of `F` enters with its own sign, the excluded ones
!!      (`vprec`, `seaice_melt`, the cavity term, land) do not, a vanished top layer reads the first
!!      live one.
!!   3. `budget_closure`: the salt the deposit mirrors into `salt_budget_surface` equals the change
!!      of the content; the same run's budget WITHOUT the source (a control run with the knob off) is
!!      far from closing, so the check would catch a missing source term.
!!   4. `knob_off_bitident`: with the knob off the slot is neither read nor written and the assembler
!!      reproduces the legacy sums (`Q_salt`, `Q_heat` incl. the mass enthalpy).
!!   5. `validate`: `none`/`virtual` accepted, `mass` needs a remapped coordinate, unknown values refused, `virtual` needs the
!!      component set.
!!   6. `cice_no_double_count`: the CICE-port sea ice under the prescribed atmosphere through the
!!      production engine: the ice-covered share of the rain dilutes through `salt_ice`, the open
!!      share through `salt_fw`, and the two add to the whole precipitation.
module test_ocean_surface_freshwater
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_assemble, &
                                     ocean_surface_flux_fw_virtual_step, &
                                     parse_surface_freshwater, SURF_FW_NONE, SURF_FW_VIRTUAL, &
                                     SURF_FW_MASS, SURF_FW_INVALID, SEAWATER_CP
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use test_cice_synth_atm, only: write_synth_atm, remove_synth_atm
   implicit none
   private

   public :: collect_ocean_surface_freshwater_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H_LAYER = 10.0_wp
   real(wp), parameter :: DT = 3600.0_wp
   character(len=*), parameter :: NL = new_line("a")
   character(len=*), parameter :: VFW = "&ocean_forcing_nml enable_components = .true., freshwater = 'virtual' /"//NL
   character(len=*), parameter :: FN_ATM = "test_fw_cice_atm.nc"
   character(len=*), parameter :: FN_LST = "test_fw_cice_atm.txt"
   real(wp), parameter :: DT_ENGINE = 600.0_wp
   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_surface_freshwater_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fw_rain_resting_column", test_rain_resting_column), &
                  new_unittest("fw_components_sign_composition", test_components_sign_composition), &
                  new_unittest("fw_budget_closure", test_budget_closure), &
                  new_unittest("fw_knob_off_bitident", test_knob_off_bitident), &
                  new_unittest("fw_validate", test_validate), &
                  new_unittest("fw_cover_assembler", test_cover_assembler), &
                  new_unittest("fw_cice_salt_closure", test_cice_salt_closure), &
                  new_unittest("fw_cice_no_double_count", test_cice_no_double_count) &
                  ]
   end subroutine collect_ocean_surface_freshwater_tests

   subroutine setup_column(grid, ms, sf, virtual, s_top, t0)
      !! A components-on slot over a 6x4 tile, uniform h, T, S (`s_top` at the surface layer).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      logical, intent(in) :: virtual
      real(wp), intent(in) :: s_top, t0

      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      if (virtual) call sf%set_freshwater(SURF_FW_VIRTUAL)
      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_temperature)%c = t0
      ms%tracers(ms%idx_salinity)%c = 35.0_wp
      ms%tracers(ms%idx_salinity)%c(:, :, NZ) = s_top
      sf%has_heat = .true.
   end subroutine setup_column

   subroutine map_all(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
   end subroutine map_all

   subroutine unmap_all(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine unmap_all

   subroutine one_thermo_step(grid, ms, sf)
      !! The engine's per-thermo-step order: the writer, the assembler, the deposit.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
   end subroutine one_thermo_step

   ! -----------------------------------------------------------------
   ! 1. rain on a resting column
   ! -----------------------------------------------------------------

   subroutine test_rain_resting_column(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms, ms_c
      type(ocean_surface_flux_t) :: sf, sf_c
      real(wp), parameter :: P = 5.0e-4_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 10.0_wp
      integer, parameter :: NSTEP = 100
      real(wp), allocatable :: s_top(:, :), t_all(:, :, :), h_all(:, :, :), s_low(:, :, :)
      real(wp) :: a, closed, expo, dt_ctl, tmax
      integer :: n, ip, jp

      call setup_column(grid, ms, sf, .true., S0, T0)
      sf%lprec = P
      sf%heat_content_lprec = SEAWATER_CP*T0*P
      ! the control: same forcing, knob off
      call setup_column(grid, ms_c, sf_c, .false., S0, T0)
      sf_c%lprec = P
      sf_c%heat_content_lprec = SEAWATER_CP*T0*P

      call map_all(ms, sf)
      !$acc update device(sf%lprec, sf%heat_content_lprec)
      do n = 1, NSTEP
         call one_thermo_step(grid, ms, sf)
      end do
      !$acc update self(ms%tracers(ms%idx_salinity)%c, ms%tracers(ms%idx_temperature)%c, ms%h_layer)
      call unmap_all(ms, sf)

      call map_all(ms_c, sf_c)
      !$acc update device(sf_c%lprec, sf_c%heat_content_lprec)
      call one_thermo_step(grid, ms_c, sf_c)
      !$acc update self(ms_c%tracers(ms_c%idx_temperature)%c)
      call unmap_all(ms_c, sf_c)

      ip = grid%nx_total/2
      jp = grid%ny_total/2
      a = DT*P/(sf%rho0*H_LAYER)
      closed = S0*(1.0_wp - a)**NSTEP
      expo = S0*exp(-a*real(NSTEP, wp))
      allocate (s_top(grid%nx_total, grid%ny_total), source=ms%tracers(ms%idx_salinity)%c(:, :, NZ))
      allocate (s_low(grid%nx_total, grid%ny_total, NZ - 1), &
                source=ms%tracers(ms%idx_salinity)%c(:, :, 1:NZ - 1))
      allocate (t_all, source=ms%tracers(ms%idx_temperature)%c)
      allocate (h_all, source=ms%h_layer)

      ! the fall is of the analytic size and in the right direction
      call check(error, abs(s_top(ip, jp) - closed) < 1.0e-12_wp*S0, &
                 "surface salinity after n steps = S0 (1 - dt P/(rho0 h))^n")
      if (allocated(error)) return
      call check(error, abs(s_top(ip, jp) - expo) < 1.0e-5_wp*S0, &
                 "surface salinity follows S0 exp(-P t/(rho0 h)) to the Euler error a^2 n/2")
      if (allocated(error)) return
      call check(error, S0 - s_top(ip, jp) > 0.05_wp, "the rain freshens the surface layer measurably")
      if (allocated(error)) return
      ! every column alike (uniform forcing, resting)
      call check(error, maxval(abs(s_top(NGHOST + 1:grid%nx_total - NGHOST, NGHOST + 1:grid%ny_total - NGHOST) &
                                   - s_top(ip, jp))) < 1.0e-13_wp, "uniform rain freshens every column alike")
      if (allocated(error)) return
      call check(error, maxval(abs(s_low - 35.0_wp)) == 0.0_wp, "only the surface layer is diluted")
      if (allocated(error)) return
      ! heat unchanged by the salt term, mass unchanged
      tmax = maxval(abs(t_all - T0))
      call check(error, tmax < 1.0e-12_wp, "virtual fresh water leaves the temperature alone (isothermal rain)")
      if (allocated(error)) return
      call check(error, maxval(abs(h_all - H_LAYER)) == 0.0_wp, "no column mass change (thickness untouched)")
      if (allocated(error)) return
      ! the control: with the knob off, the enthalpy of the (unadded) rain warms the layer
      dt_ctl = ms_c%tracers(ms_c%idx_temperature)%c(ip, jp, NZ) - T0
      call check(error, abs(dt_ctl - DT*T0*P/(sf_c%rho0*H_LAYER)) < 1.0e-9_wp*abs(dt_ctl), &
                 "control (knob off): the rain enthalpy warms by dt T P/(rho0 h)")
      if (allocated(error)) return
      call check(error, dt_ctl > 1.0e-3_wp, "the control warming is large enough that the heat check could fail")
      call sf%destroy(); call ms%destroy()
      call sf_c%destroy(); call ms_c%destroy()
   end subroutine test_rain_resting_column

   ! -----------------------------------------------------------------
   ! 2. sign and composition
   ! -----------------------------------------------------------------

   subroutine test_components_sign_composition(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp) :: f, s_read, expect, expect_q, worst, worst_q, worst_heat, qh_expect
      integer :: i, j, k, nxt, nyt, ng, ip, jp, n_nonzero
      integer :: i_zero, j_zero, i_land, j_land, i_thin, j_thin
      real(wp), parameter :: CPT = 11.0_wp

      call setup_column(grid, ms, sf, .true., 35.0_wp, 5.0_wp)
      nxt = grid%nx_total
      nyt = grid%ny_total
      ng = grid%nghost
      do k = 1, NZ
         do j = 1, nyt
            do i = 1, nxt
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 30.0_wp + 0.5_wp*real(i, wp) + 0.25_wp*real(j, wp) &
                                                        + 0.1_wp*real(k, wp)
               sf%lprec(i, j) = 1.0e-5_wp*real(1 + mod(i, 3), wp)
               sf%fprec(i, j) = 2.0e-6_wp*real(mod(j, 2), wp)
               sf%evap(i, j) = -3.0e-6_wp*real(mod(i + j, 4), wp)
               if (i == 4) sf%lrunoff(i, j) = 1.0e-6_wp
               if (j == 3) sf%frunoff(i, j) = 5.0e-7_wp
            end do
         end do
      end do
      ! components that must NOT enter F: set large and everywhere
      sf%vprec = 7.0e-5_wp
      sf%seaice_melt = 9.0e-5_wp
      sf%heat_content_vprec = CPT
      sf%heat_content_seaice_melt = 2.0_wp
      ! other salt sources, summed unchanged
      sf%salt_flux = 1.0e-7_wp
      sf%salt_ice = -2.0e-7_wp
      sf%salt_cavity = 3.0e-7_wp
      ! rain + equal evaporation: exactly zero F
      i_zero = ng + 2
      j_zero = ng + 2
      sf%lprec(i_zero, j_zero) = 4.0e-5_wp
      sf%evap(i_zero, j_zero) = -4.0e-5_wp
      sf%fprec(i_zero, j_zero) = 0.0_wp
      sf%lrunoff(i_zero, j_zero) = 0.0_wp
      sf%frunoff(i_zero, j_zero) = 0.0_wp
      ! a land cell
      i_land = ng + 3
      j_land = ng + 1
      ms%wet_mask(i_land, j_land) = 0.0_wp
      ! a collapsed top layer: the first live layer is read
      i_thin = ng + 5
      j_thin = ng + 3
      ms%h_layer(i_thin, j_thin, NZ) = 0.0_wp

      call map_all(ms, sf)
      !$acc update device(ms%wet_mask, ms%h_layer)
      !$acc update device(sf%lprec, sf%fprec, sf%evap, sf%lrunoff, sf%frunoff, sf%vprec, sf%seaice_melt)
      !$acc update device(sf%heat_content_vprec, sf%heat_content_seaice_melt)
      !$acc update device(sf%salt_flux, sf%salt_ice, sf%salt_cavity)
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      !$acc update self(sf%salt_fw, sf%Q_salt, sf%Q_heat)
      call unmap_all(ms, sf)

      worst = 0.0_wp
      worst_q = 0.0_wp
      worst_heat = 0.0_wp
      n_nonzero = 0
      do j = ng + 1, nyt - ng
         do i = ng + 1, nxt - ng
            f = 1.0e-5_wp*real(1 + mod(i, 3), wp) + 2.0e-6_wp*real(mod(j, 2), wp) &
                - 3.0e-6_wp*real(mod(i + j, 4), wp)
            if (i == 4) f = f + 1.0e-6_wp
            if (j == 3) f = f + 5.0e-7_wp
            if (i == i_zero .and. j == j_zero) f = 0.0_wp
            k = NZ
            if (i == i_thin .and. j == j_thin) k = NZ - 1
            s_read = 30.0_wp + 0.5_wp*real(i, wp) + 0.25_wp*real(j, wp) + 0.1_wp*real(k, wp)
            expect = -s_read*f
            expect_q = expect + 1.0e-7_wp - 2.0e-7_wp + 3.0e-7_wp
            qh_expect = CPT + 2.0_wp
            if (i == i_land .and. j == j_land) then
               expect = 0.0_wp
               expect_q = 0.0_wp
               qh_expect = 0.0_wp
            end if
            if (.not. ieee_is_finite(sf%salt_fw(i, j))) then
               call check(error, .false., "salt_fw is finite")
               return
            end if
            worst = max(worst, abs(sf%salt_fw(i, j) - expect)/(abs(expect) + 1.0e-12_wp))
            worst_q = max(worst_q, abs(sf%Q_salt(i, j) - expect_q)/(abs(expect_q) + 1.0e-12_wp))
            worst_heat = max(worst_heat, abs(sf%Q_heat(i, j) - qh_expect))
            if (abs(expect) > 0.0_wp) n_nonzero = n_nonzero + 1
         end do
      end do
      call check(error, n_nonzero >= 12, "the test states differ from cell to cell and are mostly nonzero")
      if (allocated(error)) return
      call check(error, worst < 1.0e-13_wp, "salt_fw = -S_surf (lprec + fprec + evap + lrunoff + frunoff) per cell")
      if (allocated(error)) return
      call check(error, worst_q < 1.0e-12_wp, &
                 "Q_salt = salt_flux + salt_ice + salt_cavity + salt_fw; vprec and seaice_melt do not enter")
      if (allocated(error)) return
      call check(error, worst_heat < 1.0e-9_wp, &
                 "Q_heat carries only the enthalpy of the mass flux NOT in F (vprec, seaice_melt)")
      if (allocated(error)) return
      ip = i_zero
      jp = j_zero
      call check(error, sf%salt_fw(ip, jp) == 0.0_wp, "rain + equal evaporation gives zero, whatever vprec/seaice_melt hold")
      if (allocated(error)) return
      call check(error, sf%salt_fw(i_land, j_land) == 0.0_wp .and. sf%Q_salt(i_land, j_land) == 0.0_wp, &
                 "a land cell receives nothing")
      if (allocated(error)) return
      call check(error, sf%salt_fw(1, 1) == 0.0_wp, "the writer fills owned cells only (ghosts are exchanged)")
      call sf%destroy(); call ms%destroy()
   end subroutine test_components_sign_composition

   ! -----------------------------------------------------------------
   ! 3. budget closure
   ! -----------------------------------------------------------------

   function owned_content(grid, ms) result(c)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      real(wp) :: c
      integer :: i, j, k, ng

      ng = grid%nghost
      c = 0.0_wp
      do k = 1, NZ
         do j = ng + 1, grid%ny_total - ng
            do i = ng + 1, grid%nx_total - ng
               c = c + ms%h_layer(i, j, k)*ms%tracers(ms%idx_salinity)%c(i, j, k)
            end do
         end do
      end do
   end function owned_content

   function owned_budget(grid, ms) result(b)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      real(wp) :: b
      integer :: i, j, k, ng

      ng = grid%nghost
      b = 0.0_wp
      do k = 1, NZ
         do j = ng + 1, grid%ny_total - ng
            do i = ng + 1, grid%nx_total - ng
               b = b + ms%salt_budget_surface(i, j, k)
            end do
         end do
      end do
   end function owned_budget

   subroutine test_budget_closure(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms, ms_n
      type(ocean_surface_flux_t) :: sf, sf_n
      integer, parameter :: NSTEP = 40
      real(wp) :: c0, c1, bud, bud_missing, dc
      integer :: n, i, j

      call setup_column(grid, ms, sf, .true., 34.0_wp, 8.0_wp)
      call setup_column(grid, ms_n, sf_n, .false., 34.0_wp, 8.0_wp)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            sf%lprec(i, j) = 2.0e-4_wp*real(mod(i, 3), wp)
            sf%evap(i, j) = -1.0e-4_wp*real(mod(j, 2) + 1, wp)
            sf%fprec(i, j) = 5.0e-5_wp
         end do
      end do
      sf_n%lprec = sf%lprec
      sf_n%evap = sf%evap
      sf_n%fprec = sf%fprec
      ! a salt source of another kind alongside, so the closure is of the SUM
      sf%salt_ice = 3.0e-5_wp
      sf_n%salt_ice = 3.0e-5_wp
      sf%has_salt = .true.
      sf_n%has_salt = .true.

      c0 = owned_content(grid, ms)
      call map_all(ms, sf)
      !$acc update device(sf%lprec, sf%evap, sf%fprec, sf%salt_ice)
      do n = 1, NSTEP
         call one_thermo_step(grid, ms, sf)
      end do
      !$acc update self(ms%tracers(ms%idx_salinity)%c, ms%salt_budget_surface)
      call unmap_all(ms, sf)
      c1 = owned_content(grid, ms)
      bud = owned_budget(grid, ms)
      dc = c1 - c0

      call map_all(ms_n, sf_n)
      !$acc update device(sf_n%lprec, sf_n%evap, sf_n%fprec, sf_n%salt_ice)
      do n = 1, NSTEP
         call one_thermo_step(grid, ms_n, sf_n)
      end do
      !$acc update self(ms_n%salt_budget_surface)
      call unmap_all(ms_n, sf_n)
      bud_missing = owned_budget(grid, ms_n)

      write (output_unit, "(a,3es14.5)") "  fw budget: content change, mirrored budget, budget without the source:", &
         dc, bud, bud_missing
      call check(error, abs(dc) > 1.0e-3_wp, "the surface salt content changes measurably")
      if (allocated(error)) return
      call check(error, abs(dc - bud) < 1.0e-11_wp*abs(c0), "content change = mirrored source budget (round-off)")
      if (allocated(error)) return
      call check(error, abs(dc - bud_missing) > 1.0e-3_wp*abs(dc), &
                 "a budget missing the fresh-water source would fail the closure by far more than round-off")
      call sf%destroy(); call ms%destroy()
      call sf_n%destroy(); call ms_n%destroy()
   end subroutine test_budget_closure

   ! -----------------------------------------------------------------
   ! 4. knob off
   ! -----------------------------------------------------------------

   subroutine test_knob_off_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp) :: massin, massout, qh, qs, worst_h, worst_s, sst
      real(wp), allocatable :: fw_before(:, :)
      integer :: i, j, ng

      call setup_column(grid, ms, sf, .false., 35.0_wp, 7.0_wp)
      ng = grid%nghost
      sst = 7.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            sf%lprec(i, j) = 1.0e-5_wp*real(i, wp)
            sf%fprec(i, j) = 2.0e-6_wp*real(j, wp)
            sf%evap(i, j) = -3.0e-6_wp*real(i + j, wp)
            sf%heat_content_lprec(i, j) = SEAWATER_CP*sst*sf%lprec(i, j)
            sf%heat_content_fprec(i, j) = SEAWATER_CP*sst*sf%fprec(i, j)
            sf%q_lat(i, j) = -50.0_wp
            sf%q_sw(i, j) = 120.0_wp + real(i, wp)
         end do
      end do
      sf%salt_flux = 1.0e-7_wp
      sf%salt_ice = -2.0e-7_wp
      sf%salt_cavity = 3.0e-7_wp
      sf%salt_fw = 1.0e3_wp   ! garbage: nobody may read it with the knob off
      allocate (fw_before, source=sf%salt_fw)

      call map_all(ms, sf)
      !$acc update device(sf%lprec, sf%fprec, sf%evap, sf%heat_content_lprec, sf%heat_content_fprec)
      !$acc update device(sf%q_lat, sf%q_sw, sf%salt_flux, sf%salt_ice, sf%salt_cavity, sf%salt_fw)
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      !$acc update self(sf%salt_fw, sf%Q_salt, sf%Q_heat)
      call unmap_all(ms, sf)

      call check(error, all(sf%salt_fw == fw_before), "knob off: the writer does not touch the slot")
      if (allocated(error)) return
      worst_h = 0.0_wp
      worst_s = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            massin = sf%heat_content_lprec(i, j) + sf%heat_content_fprec(i, j) + sf%heat_content_vprec(i, j) &
                     + sf%heat_content_lrunoff(i, j) + sf%heat_content_frunoff(i, j) &
                     + sf%heat_content_seaice_melt(i, j)
            massout = SEAWATER_CP*sst*sf%evap(i, j)
            qh = (0.0_wp + (sf%q_sw(i, j) + 0.0_wp) + 0.0_wp + sf%q_lat(i, j) + 0.0_wp &
                  + (0.0_wp + 0.0_wp) + 0.0_wp + massin + massout)
            qs = 0.0_wp + (1.0e-7_wp - 2.0e-7_wp) + 3.0e-7_wp
            worst_h = max(worst_h, abs(sf%Q_heat(i, j) - qh)/(abs(qh) + 1.0_wp))
            worst_s = max(worst_s, abs(sf%Q_salt(i, j) - qs)/abs(qs))
         end do
      end do
      call check(error, worst_s < 1.0e-15_wp, "knob off: Q_salt is the legacy sum, the garbage slot is not in it")
      if (allocated(error)) return
      call check(error, worst_h < 1.0e-15_wp, "knob off: Q_heat keeps the enthalpy of the mass fluxes")
      call sf%destroy(); call ms%destroy()
   end subroutine test_knob_off_bitident

   ! -----------------------------------------------------------------
   ! 4b. the ice-shelf-cover twin of the assembler
   ! -----------------------------------------------------------------

   subroutine test_cover_assembler(error)
      !! `salt_fw` is an ATMOSPHERIC term: under a shelf it is scaled by the open fraction `1 - cover_frac`,
      !! like `salt_flux`/`salt_ice`, while the cavity's own `salt_cavity`/`heat_cavity` pass unmasked; the
      !! mass enthalpy that stays (`vprec`, `seaice_melt`) is masked too.  Independent values per cell.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: cover(:, :)
      real(wp) :: of, f, s_top, e_s, e_h, worst_s, worst_h
      integer :: i, j, ng

      call setup_column(grid, ms, sf, .true., 35.0_wp, 5.0_wp)
      ng = grid%nghost
      allocate (cover(grid%nx_total, grid%ny_total))
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            cover(i, j) = 0.25_wp*real(mod(i + 2*j, 5), wp)
            ms%tracers(ms%idx_salinity)%c(i, j, NZ) = 32.0_wp + 0.3_wp*real(i, wp) + 0.1_wp*real(j, wp)
            sf%lprec(i, j) = 2.0e-5_wp + 1.0e-6_wp*real(i, wp)
            sf%evap(i, j) = -5.0e-6_wp*real(1 + mod(j, 2), wp)
         end do
      end do
      sf%salt_flux = 1.0e-7_wp
      sf%salt_cavity = 4.0e-7_wp
      sf%heat_cavity = -3.0_wp
      sf%heat_content_vprec = 8.0_wp
      sf%heat_content_seaice_melt = 1.0_wp
      sf%q_sw = 100.0_wp
      sf%heat_content_lprec = 999.0_wp   ! must not enter under the knob

      call map_all(ms, sf)
      !$acc update device(sf%lprec, sf%evap, sf%salt_flux, sf%salt_cavity, sf%heat_cavity, sf%heat_content_vprec)
      !$acc update device(sf%heat_content_seaice_melt, sf%q_sw, sf%heat_content_lprec)
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms, cover_frac=cover)
      !$acc update self(sf%Q_salt, sf%Q_heat, sf%salt_fw)
      call unmap_all(ms, sf)

      worst_s = 0.0_wp
      worst_h = 0.0_wp
      do j = ng + 1, grid%ny_total - ng
         do i = ng + 1, grid%nx_total - ng
            of = 1.0_wp - cover(i, j)
            s_top = 32.0_wp + 0.3_wp*real(i, wp) + 0.1_wp*real(j, wp)
            f = 2.0e-5_wp + 1.0e-6_wp*real(i, wp) - 5.0e-6_wp*real(1 + mod(j, 2), wp)
            e_s = of*(1.0e-7_wp - s_top*f) + 4.0e-7_wp
            e_h = of*(100.0_wp + 8.0_wp + 1.0_wp) - 3.0_wp
            worst_s = max(worst_s, abs(sf%Q_salt(i, j) - e_s)/abs(e_s))
            worst_h = max(worst_h, abs(sf%Q_heat(i, j) - e_h)/abs(e_h))
         end do
      end do
      call check(error, maxval(cover) > 0.5_wp .and. minval(cover) == 0.0_wp, "the cover fraction is nonuniform")
      if (allocated(error)) return
      call check(error, worst_s < 1.0e-12_wp, "cover twin: Q_salt = (1-c)(salt_flux + salt_fw) + salt_cavity")
      if (allocated(error)) return
      call check(error, worst_h < 1.0e-12_wp, "cover twin: Q_heat masks the remaining mass enthalpy, keeps heat_cavity")
      call sf%destroy(); call ms%destroy()
   end subroutine test_cover_assembler

   ! -----------------------------------------------------------------
   ! 5. configure-time refusal
   ! -----------------------------------------------------------------

   function valid_ierr(extra) result(ierr)
      character(len=*), intent(in) :: extra
      integer :: ierr
      type(config_t) :: cfg
      integer :: ierr_read

      call read_config_from_string("&sim_nml sim_type = 'ocean' /"//NL// &
                                   "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
                                   extra//NL// &
                                   "&ocean_diag_nml enabled = .false. /"//NL// &
                                   "&output_nml output_to_file = .false. /"//NL, cfg, ierr=ierr_read)
      if (ierr_read /= OCEAN_STATUS_OK) then
         ierr = ierr_read
         return
      end if
      call validate_config(cfg, ierr)
   end function valid_ierr

   subroutine test_validate(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, parse_surface_freshwater("none") == SURF_FW_NONE .and. &
                 parse_surface_freshwater("virtual") == SURF_FW_VIRTUAL .and. &
                 parse_surface_freshwater("mass") == SURF_FW_MASS .and. &
                 parse_surface_freshwater("bogus") == SURF_FW_INVALID, "the knob values parse")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true. /") == OCEAN_STATUS_OK, &
                 "baseline namelist validates (knob absent)")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true., freshwater = 'none' /") &
                 == OCEAN_STATUS_OK, "freshwater = 'none' is accepted")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true., freshwater = 'virtual' /") &
                 == OCEAN_STATUS_OK, "freshwater = 'virtual' is accepted")
      if (allocated(error)) return
      ! 'mass' ships (test_ocean_surface_fw_mass): the default sigma coordinate remaps, so it is
      ! accepted; `eulerian_z` cannot give the volume back to a remap target and stays refused.
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true., freshwater = 'mass' /"//NL// &
                                   "&ocean_bt_nml auto_n_inner = .true. /") &
                 == OCEAN_STATUS_OK, "freshwater = 'mass' is accepted (default sigma coordinate)")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true., freshwater = 'mass' /") &
                 /= OCEAN_STATUS_OK, "freshwater = 'mass' is refused on the unsplit solver (no remap)")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true., freshwater = 'mass' /"//NL// &
                                   "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
                                   "&vcoord_nml vcoord_type = 'eulerian_z' /") /= OCEAN_STATUS_OK, &
                 "freshwater = 'mass' is refused with eulerian_z")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml freshwater = 'virtual' /") /= OCEAN_STATUS_OK, &
                 "freshwater = 'virtual' without the component set is refused")
      if (allocated(error)) return
      call check(error, parse_surface_freshwater("off") == SURF_FW_INVALID .and. &
                 parse_surface_freshwater("real") == SURF_FW_INVALID, "no aliases beyond the enum")
      if (allocated(error)) return
      ! SIS2-derived ice: the lumped mode and the open-water snowfall path double count
      call check(error, valid_ierr(VFW//"&ocean_ice_nml enable = .true., ncat = 5 /") == OCEAN_STATUS_OK, &
                 "virtual with SIS2-derived ice at ncat = 5 is accepted")
      if (allocated(error)) return
      call check(error, valid_ierr(VFW//"&ocean_ice_nml enable = .true., ncat = 1 /") /= OCEAN_STATUS_OK, &
                 "virtual with SIS2-derived ice at ncat = 1 (lumped, aice0 = 1) is refused")
      if (allocated(error)) return
      call check(error, valid_ierr(VFW//"&ocean_ice_nml enable = .true., ncat = 5, snowfall = 2.0e-5 /") &
                 /= OCEAN_STATUS_OK, "virtual with SIS2-derived ice snowfall /= 0 (counted twice) is refused")
      if (allocated(error)) return
      call check(error, valid_ierr("&ocean_forcing_nml enable_components = .true. /"//NL// &
                                   "&ocean_ice_nml enable = .true., ncat = 1 /") == OCEAN_STATUS_OK, &
                 "the same ice with the knob off is untouched")
   end subroutine test_validate

   ! -----------------------------------------------------------------
   ! 6. CICE-port sea ice: the two shares of the rain
   ! -----------------------------------------------------------------

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function cice_nml(freshwater) result(nml)
      !! A closed spherical high-latitude basin, the ice 0.6 covered, the synthetic prescribed atmosphere.
      character(len=*), intent(in) :: freshwater
      character(len=:), allocatable :: nml

      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0, start_date = '2000-01-01' /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&tracer_nml initial_temperature = -1.0, initial_salinity = 34.0, "// &
            "T_init_surface = -1.0, T_init_bottom = -1.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            "&grid_nml nx = 12, ny = 10, nghost = 3, dx = 0.25, dy = 0.25 /"//NL// &
            "&ocean_grid_nml grid_config = 'spherical', lon_west = 20.0, lat_south = 70.0, "// &
            "rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
            "&physics_nml wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true., freshwater = '"//freshwater//"' /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = 1 /"//NL// &
            "&ocean_ice_nml enable = .true., scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, "// &
            "dynamics = .true., transport = .true. /"//NL// &
            "&ocean_ice_thermo_nml update_ocn_f = .true. /"//NL// &
            "&ocean_ice_dyn_nml advection = 'upwind', ndte = 40 /"//NL// &
            "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.6, h_snow = 0.1, t_ice = -8.0 /"//NL// &
            "&ocean_atm_forcing_nml enable = .true., sw_interp = 'linear', fluxes_regrid = 'bilinear', "// &
            "time_mode = 'interannual', wind_relative = .false., "// &
            "filelist_tas = '"//FN_LST//"', filelist_huss = '"//FN_LST//"', filelist_uas = '"//FN_LST//"', "// &
            "filelist_vas = '"//FN_LST//"', filelist_psl = '"//FN_LST//"', filelist_rlds = '"//FN_LST//"', "// &
            "filelist_rsds = '"//FN_LST//"', filelist_prra = '"//FN_LST//"', filelist_prsn = '"//FN_LST//"' /"//NL// &
            "&ocean_bulk_flux_nml enable = .true. /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function cice_nml

   subroutine run_one_step(fw, ok, salt_ice, salt_fw, q_salt, q_heat, lprec, fprec, evap, s_top, frain, aice_init, &
                           sss_ocn, aice0_open)
      !! One production outer step (`engine_step` -> `engine_step_ice` -> `engine_step_finalize`) of the
      !! `scheme = 'cice'` case, the slots the ocean assembled read back to the host.
      character(len=*), intent(in) :: fw
      logical, intent(out) :: ok
      real(wp), allocatable, intent(out) :: salt_ice(:, :), salt_fw(:, :), q_salt(:, :), q_heat(:, :)
      real(wp), allocatable, intent(out) :: lprec(:, :), fprec(:, :), evap(:, :), s_top(:, :)
      real(wp), allocatable, intent(out) :: frain(:, :), aice_init(:, :), sss_ocn(:, :), aice0_open(:, :)
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t
      integer :: ierr

      ok = .false.
      call write_synth_atm("mixed", FN_ATM, FN_LST)
      call ensure_comm()
      call read_config_from_string(cice_nml(fw), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      call engine_step(engine, DT_ENGINE, t, ierr=ierr)
      if (ierr == OCEAN_STATUS_OK) call engine_step_ice(engine, cfg, DT_ENGINE, t, ierr=ierr)
      if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DT_ENGINE, t, ierr=ierr)
      if (ierr == OCEAN_STATUS_OK) then
         associate (sf => engine%state%surface_flux, ms => engine%state%multilayer, w => engine%state%icec%work, &
                    c => engine%state%ice%cice, o => engine%state%icec%oin)
            !$acc update self(sf%salt_ice, sf%salt_fw, sf%Q_salt, sf%Q_heat, sf%lprec, sf%fprec, sf%evap)
            !$acc update self(w%frain, c%aice_init, c%aice0, o%sss_ocn)
            associate (cS => ms%tracers(ms%idx_salinity)%c)
               !$acc update self(cS)
               allocate (s_top, source=cS(:, :, ms%nz_ml))
            end associate
            allocate (salt_ice, source=sf%salt_ice)
            allocate (salt_fw, source=sf%salt_fw)
            allocate (q_salt, source=sf%Q_salt)
            allocate (q_heat, source=sf%Q_heat)
            allocate (lprec, source=sf%lprec)
            allocate (fprec, source=sf%fprec)
            allocate (evap, source=sf%evap)
            allocate (frain, source=w%frain)
            allocate (aice_init, source=c%aice_init)
            allocate (aice0_open, source=c%aice0)
            allocate (sss_ocn, source=o%sss_ocn)
            ok = .true.
         end associate
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
      call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine run_one_step

   subroutine test_cice_no_double_count(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: si_v(:, :), fw_v(:, :), qs_v(:, :), qh_v(:, :), lp_v(:, :), fp_v(:, :), ev_v(:, :)
      real(wp), allocatable :: st_v(:, :), fr_v(:, :), ai_v(:, :), ss_v(:, :), a0_v(:, :)
      real(wp), allocatable :: si_n(:, :), fw_n(:, :), qs_n(:, :), qh_n(:, :), lp_n(:, :), fp_n(:, :), ev_n(:, :)
      real(wp), allocatable :: st_n(:, :), fr_n(:, :), ai_n(:, :), ss_n(:, :), a0_n(:, :)
      real(wp) :: e_fw, e_q, e_ice, e_tot, term, scl, n_nz
      integer :: i, j, ng, is, ie, js, je
      logical :: ok

      call run_one_step("virtual", ok, si_v, fw_v, qs_v, qh_v, lp_v, fp_v, ev_v, st_v, fr_v, ai_v, ss_v, a0_v)
      call check(error, ok, "virtual engine step ran")
      if (allocated(error)) return
      call run_one_step("none", ok, si_n, fw_n, qs_n, qh_n, lp_n, fp_n, ev_n, st_n, fr_n, ai_n, ss_n, a0_n)
      call check(error, ok, "none engine step ran")
      if (allocated(error)) return

      ng = 3
      is = ng + 1
      ie = ng + 12
      js = ng + 1
      je = ng + 10
      e_fw = 0.0_wp
      e_q = 0.0_wp
      e_ice = 0.0_wp
      e_tot = 0.0_wp
      n_nz = 0.0_wp
      do j = js, je
         do i = is, ie
            ! open share: salt_fw = -S_surf (lprec + fprec + evap), independently from the state
            term = -st_v(i, j)*(lp_v(i, j) + fp_v(i, j) + ev_v(i, j))
            e_fw = max(e_fw, abs(fw_v(i, j) - term)/(abs(term) + 1.0e-12_wp))
            if (abs(term) > 0.0_wp) n_nz = n_nz + 1.0_wp
            ! the whole salt source is the ice's hand-over plus the open-water term, nothing else
            e_q = max(e_q, abs(qs_v(i, j) - (si_v(i, j) + fw_v(i, j)))/(abs(qs_v(i, j)) + 1.0e-12_wp))
            ! knob off: no fresh-water term, the ice's alone
            e_q = max(e_q, abs(qs_n(i, j) - si_n(i, j))/(abs(qs_n(i, j)) + 1.0e-12_wp), abs(fw_n(i, j)))
            ! ice share: the rain through the ice stays in `fresh` under the knob, so the hand-over
            ! differs by exactly -S_ocn * frain * aice_init
            term = -ss_n(i, j)*fr_n(i, j)*ai_n(i, j)
            scl = abs(term) + 1.0e-12_wp
            e_ice = max(e_ice, abs((si_v(i, j) - si_n(i, j)) - term)/scl)
            ! the two shares of the rain add to the whole precipitation (open by aice0, ice by aice_init)
            e_tot = max(e_tot, abs((lp_v(i, j) + fr_v(i, j)*ai_v(i, j)) - fr_v(i, j))/fr_v(i, j))
         end do
      end do
      write (output_unit, "(a,5es11.3)") "  fw cice: worst rel err open-share, Q_salt sum, ice share, shares-sum-to-P, n:", &
         e_fw, e_q, e_ice, e_tot, n_nz
      call check(error, n_nz > 20.0_wp, "the open-water term is nonzero in the cells compared")
      if (allocated(error)) return
      call check(error, e_fw < 1.0e-13_wp, "open share of the rain: salt_fw = -S (lprec + fprec + evap)")
      if (allocated(error)) return
      call check(error, e_q < 1.0e-12_wp, "Q_salt = ice hand-over + salt_fw (knob on), ice hand-over alone (knob off)")
      if (allocated(error)) return
      call check(error, e_ice < 1.0e-9_wp, "the ice-covered share of the rain dilutes through the hand-over")
      if (allocated(error)) return
      call check(error, e_tot < 1.0e-9_wp, "open and ice-covered shares of the rain add to the whole precipitation")
      if (allocated(error)) return
   end subroutine test_cice_no_double_count

   subroutine test_cice_salt_closure(error)
      !! Coupled closure under `freshwater = "virtual"` with the CICE-port ice in a CLOSED basin: the ocean
      !! salt content changes by exactly what the surface deposit mirrored into the salt budget (the ice
      !! hand-over incl. the kept rain-through term, plus `salt_fw`), over several outer steps; and the
      !! knob-off run, whose rain never reaches the salt, ends measurably elsewhere.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: d_v, b_v, d_n, b_n, scale
      logical :: ok

      call cice_content_run("virtual", 6, ok, d_v, b_v, scale)
      call check(error, ok, "virtual closed-basin run")
      if (allocated(error)) return
      call cice_content_run("none", 6, ok, d_n, b_n, scale)
      call check(error, ok, "knob-off closed-basin run")
      if (allocated(error)) return
      write (output_unit, "(a,4es13.5)") "  fw cice closure: d(content) V, budget V, d(content) none, budget none:", &
         d_v, b_v, d_n, b_n
      call check(error, abs(d_v) > 1.0e-9_wp*scale, "the salt content changes")
      if (allocated(error)) return
      call check(error, abs(d_v - b_v) < 1.0e-10_wp*scale + 1.0e-8_wp*abs(d_v), &
                 "virtual: d(salt content) = mirrored surface budget (ice hand-over + salt_fw)")
      if (allocated(error)) return
      call check(error, abs(d_n - b_n) < 1.0e-10_wp*scale + 1.0e-8_wp*abs(d_n), "none: the same closure")
      if (allocated(error)) return
      call check(error, abs(d_v - d_n) > 1.0e-3_wp*abs(d_v), "the fresh-water term moves the salt content")
   end subroutine test_cice_salt_closure

   subroutine cice_content_run(fw, nstep, ok, dcontent, budget, scale)
      character(len=*), intent(in) :: fw
      integer, intent(in) :: nstep
      logical, intent(out) :: ok
      real(wp), intent(out) :: dcontent, budget, scale
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, c0, c1
      integer :: ierr, n

      ok = .false.
      dcontent = 0.0_wp
      budget = 0.0_wp
      scale = 1.0_wp
      call write_synth_atm("mixed", FN_ATM, FN_LST)
      call ensure_comm()
      call read_config_from_string(cice_nml(fw), cfg, ierr=ierr)
      if (ierr == OCEAN_STATUS_OK) call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      call engine_enter_data(engine, cfg)
      c0 = salt_content(engine)
      scale = abs(c0)
      t = 0.0_wp
      do n = 1, nstep
         call engine_step(engine, DT_ENGINE, t, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_ice(engine, cfg, DT_ENGINE, t, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DT_ENGINE, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT_ENGINE
      end do
      if (ierr == OCEAN_STATUS_OK) then
         c1 = salt_content(engine)
         dcontent = c1 - c0
         budget = salt_budget(engine)
         ok = .true.
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
      call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine cice_content_run

   function salt_content(engine) result(c)
      !! `sum(h*S*area)` over the owned cells (device -> host first).
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp) :: c
      integer :: i, j, k, ng

      associate (cS => engine%state%multilayer%tracers(engine%state%multilayer%idx_salinity)%c)
         !$acc update self(cS)
      end associate
      !$acc update self(engine%state%multilayer%h_layer)
      ng = engine%grid%nghost
      c = 0.0_wp
      do k = 1, engine%state%multilayer%nz_ml
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               c = c + engine%state%multilayer%h_layer(i, j, k) &
                   *engine%state%multilayer%tracers(engine%state%multilayer%idx_salinity)%c(i, j, k) &
                   *engine%state%metrics%areaT(i, j)
            end do
         end do
      end do
   end function salt_content

   function salt_budget(engine) result(b)
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp) :: b
      integer :: i, j, k, ng

      !$acc update self(engine%state%multilayer%salt_budget_surface)
      ng = engine%grid%nghost
      b = 0.0_wp
      do k = 1, engine%state%multilayer%nz_ml
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               b = b + engine%state%multilayer%salt_budget_surface(i, j, k)*engine%state%metrics%areaT(i, j)
            end do
         end do
      end do
   end function salt_budget

end module test_ocean_surface_freshwater
