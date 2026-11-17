!! Runoff on the model grid (`&ocean_dataovr_nml lrunoff / frunoff`, slice B of the real fresh-water
!! mass work).  Analytical gates:
!!
!!   1. `runoff_point_source_exact`: a liquid-runoff point source read from a file grows exactly
!!      its own column by `R*dt/rho0`, T keeps its exact bits everywhere, the salt content of the
!!      column is unchanged, and the tracked mass source is `RHO_WATER * dh * area`.
!!   2. `calving_latent_heat_known_answer`: a frozen-runoff point source writes
!!      `q_fus_rof = -L_f*frunoff` at that cell only and cools its top layer by
!!      `dT = -L_f*frunoff*dt/(rho0*cp*h)` (the deposit sees the old thickness), the other layers
!!      and columns keep their exact bits.
!!   3. `runoff_on_land_fatal`: the verdict function and the OWNED-cell counting kernel (a ghost
!!      land cell and a wet cell do not count, a NaN on land does).  The abort itself is not run.
!!   4. `runoff_cyclic_interp`: `time_mode = "cyclic"` blends linearly mid-record, across the year
!!      wrap (last record to first) and reproduces itself one period later.
!!   5. `runoff_config_gates`: accepted under "mass" and "virtual", refused under "none", a blank
!!      variable is refused, and the runoff tags are not among the ones the atmosphere forcing
!!      conflicts with.
!!   6. `runoff_engine_closed_basin`: the engine path (namelist, reader, apply, assembler, mass
!!      step) raises a closed basin's mass by exactly the runoff, and the console residual closes.
module test_ocean_runoff
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER
   use rdb_config, only: config_t, read_config_from_string, validate_config, &
                         ocean_dataovr_config_t, atm_forcing_prereq_code, &
                         ATM_PREREQ_DATAOVR, ATM_PREREQ_BULK_FLUX
   use rdb_grid, only: hgrid_t
   use rdb_efp, only: efp_t, efp_to_real
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_enthalpy, only: ICE_LAT_FUS
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_ocean_data_forcing, only: ocean_data_forcing_t, ocean_data_forcing_configure, &
                                     ocean_data_forcing_apply, runoff_on_land_is_fatal, &
                                     runoff_land_count_impl
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_assemble, &
                                     ocean_surface_flux_fw_virtual_step, &
                                     SURF_FW_MASS
   use rdb_ocean_surface_fw_mass, only: ocean_surface_fw_mass_step
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_console_stats, only: compute_total_h
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            rdb_def_var_1d, nc_enddef, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_runoff_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   integer, parameter :: NXP = 6
   integer, parameter :: NYP = 4
   real(wp), parameter :: H_LAYER = 10.0_wp
   real(wp), parameter :: DT = 3600.0_wp
   real(wp), parameter :: AREA = 1.0e6_wp
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: T0 = 8.0_wp
   real(wp), parameter :: RATE = 2.0e-3_wp
      !! The point-source runoff rate (kg m-2 s-1).
   integer, parameter :: PI_FILE = 3, PJ_FILE = 2
      !! The point source's index in the FILE (physical numbering).
   character(len=*), parameter :: NL = new_line("a")
   logical, save :: comm_inited = .false.

contains

   subroutine collect_ocean_runoff_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("runoff_point_source_exact", test_point_source_exact), &
                  new_unittest("runoff_calving_latent_heat_known_answer", test_calving_latent_heat), &
                  new_unittest("runoff_on_land_fatal", test_runoff_on_land_fatal), &
                  new_unittest("runoff_cyclic_interp", test_cyclic_interp), &
                  new_unittest("runoff_config_gates", test_config_gates), &
                  new_unittest("runoff_engine_closed_basin", test_engine_closed_basin) &
                  ]
   end subroutine collect_ocean_runoff_tests

   ! -----------------------------------------------------------------
   ! scaffolding
   ! -----------------------------------------------------------------

   subroutine write_field_file(filename, nx, ny, nt, t_axis, f)
      !! `field(x, y, time)` and a `time` coordinate, the layout the reader requires.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nt
      real(wp), intent(in) :: t_axis(nt), f(nx, ny, nt)
      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_3d(ncid, "field", [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_field_file

   subroutine write_point_file(filename, nx, ny, pi, pj, val)
      !! A two-record (constant in time) field that is `val` at physical cell `(pi, pj)`, zero elsewhere.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, pi, pj
      real(wp), intent(in) :: val
      real(wp) :: f(nx, ny, 2)

      f = 0.0_wp
      f(pi, pj, :) = val
      call write_field_file(filename, nx, ny, 2, [0.0_wp, 1.0e9_wp], f)
   end subroutine write_point_file

   subroutine setup_unit(grid, ms, sf, ss, at, reader, df, cfg, ierr)
      !! A components-on, `freshwater = "mass"` slot over a 6x4 tile with the runoff tags of `cfg`
      !! registered.  Uniform thickness, stratified T and S.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_surface_stress_t), intent(inout) :: ss
      real(wp), allocatable, intent(out) :: at(:, :)
      type(ocean_data_input_t), intent(inout) :: reader
      type(ocean_data_forcing_t), intent(out) :: df
      type(ocean_dataovr_config_t), intent(in) :: cfg
      integer, intent(out) :: ierr
      type(ocean_bc_state_t) :: bc
      integer :: k

      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid, with_ideal_age=.false.)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      call sf%set_freshwater(SURF_FW_MASS)
      call ss%init(grid)
      allocate (at(grid%nx_total, grid%ny_total), source=AREA)
      ms%h_layer = H_LAYER
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%c(:, :, k) = T0 + 0.5_wp*real(k, wp)
         ms%tracers(ms%idx_salinity)%c(:, :, k) = S0 - 0.1_wp*real(k, wp)
      end do
      sf%has_heat = .true.
      call reader%init()
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df, ierr=ierr)
   end subroutine setup_unit

   subroutine map_all(ms, sf, at, reader)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(inout) :: at(:, :)
      type(ocean_data_input_t), intent(inout) :: reader

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      !$acc enter data copyin(at)
   end subroutine map_all

   subroutine unmap_all(ms, sf, at, reader)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), allocatable, intent(inout) :: at(:, :)
      type(ocean_data_input_t), intent(inout) :: reader

      !$acc exit data delete(at)
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
      call reader%exit_data()
      !$acc exit data delete(reader)
   end subroutine unmap_all

   subroutine pull_state(ms, sf)
      !! Device -> host for every array the gates read.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer :: it

      !$acc update self(ms%h_layer)
      do it = 1, size(ms%tracers)
         associate (cc => ms%tracers(it)%c)
            !$acc update self(cc)
         end associate
      end do
      !$acc update self(sf%fw_net, sf%Q_heat, sf%lrunoff, sf%frunoff, sf%q_fus_rof)
   end subroutine pull_state

   subroutine one_step(grid, ms, sf, ss, at, reader, df, t)
      !! The engine's per-step order for the runoff path: reader, apply, then the thermo step
      !! (virtual-salt writer, assembler, tracer deposit, mass step).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_surface_stress_t), intent(inout) :: ss
      real(wp), intent(in) :: at(grid%nx_total, grid%ny_total)
      type(ocean_data_input_t), intent(inout) :: reader
      type(ocean_data_forcing_t), intent(inout) :: df
      real(wp), intent(in) :: t
      type(ocean_bc_state_t) :: bc

      call ocean_data_input_update_all(reader, t)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, t, wet_mask=ms%wet_mask)
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
      call ocean_surface_fw_mass_step(grid, sf, ms, at, DT, 1.0_wp)
   end subroutine one_step

   function src_mass(ms) result(m)
      type(multilayer_state_t), intent(in) :: ms
      real(wp) :: m
      type(efp_t) :: e
      e%v = ms%mass_src_efp
      e%poison = ms%mass_src_efp_poison
      m = real(efp_to_real(e), wp)
   end function src_mass

   ! -----------------------------------------------------------------
   ! 1. point source, liquid
   ! -----------------------------------------------------------------

   subroutine test_point_source_exact(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FN = "fwb_runoff_point.nc"
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_surface_stress_t) :: ss
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_dataovr_config_t) :: cfg
      real(wp), allocatable :: at(:, :)
      real(wp), allocatable :: t_ini(:), s_ini(:)
      real(wp) :: dh, expect, content0, content1
      integer :: i, j, k, ip, jp, ierr

      call write_point_file(FN, NXP, NYP, PI_FILE, PJ_FILE, RATE)
      cfg%enable = .true.
      cfg%time_mode = "linear"
      cfg%lrunoff%file = FN
      cfg%lrunoff%var = "field"
      call setup_unit(grid, ms, sf, ss, at, reader, df, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. df%id_lrunoff > 0 .and. df%id_frunoff == 0, &
                 "lrunoff registers (and frunoff does not)")
      if (allocated(error)) return
      call check(error, sf%has_mass_flux, "registration latches has_mass_flux")
      if (allocated(error)) return
      ip = NGHOST + PI_FILE
      jp = NGHOST + PJ_FILE
      t_ini = ms%tracers(ms%idx_temperature)%c(ip, jp, :)
      s_ini = ms%tracers(ms%idx_salinity)%c(ip, jp, :)
      call map_all(ms, sf, at, reader)
      call one_step(grid, ms, sf, ss, at, reader, df, 0.0_wp)
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at, reader)

      dh = RATE*DT/sf%rho0
      expect = H_LAYER + dh
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            do k = 1, NZ
               if (i == ip .and. j == jp .and. k == NZ) then
                  call check(error, abs(ms%h_layer(i, j, k) - expect) < 1.0e-14_wp*expect, &
                             "the point-source top layer grows by R dt/rho0")
               else
                  call check(error, ms%h_layer(i, j, k) == H_LAYER, &
                             "every other layer of every other column keeps h exactly")
               end if
               if (allocated(error)) return
               call check(error, ms%tracers(ms%idx_temperature)%c(i, j, k) == T0 + 0.5_wp*real(k, wp), &
                          "T keeps its exact bits (runoff joins at the layer temperature)")
               if (allocated(error)) return
            end do
         end do
      end do
      do k = 1, NZ
         call check(error, ms%tracers(ms%idx_temperature)%c(ip, jp, k) == t_ini(k), "T of the source column is unchanged")
         if (allocated(error)) return
      end do
      content0 = 0.0_wp
      content1 = 0.0_wp
      do k = 1, NZ
         content0 = content0 + s_ini(k)*H_LAYER
         content1 = content1 + ms%tracers(ms%idx_salinity)%c(ip, jp, k)*ms%h_layer(ip, jp, k)
      end do
      call check(error, abs(content1 - content0) < 1.0e-14_wp*content0, "the salt content of the source column is unchanged")
      if (allocated(error)) return
      call check(error, ms%tracers(ms%idx_salinity)%c(ip, jp, NZ) < s_ini(NZ), "the source column is diluted")
      if (allocated(error)) return
      call check(error, abs(src_mass(ms) - RHO_WATER*dh*AREA) < 1.0e-13_wp*RHO_WATER*dh*AREA, &
                 "the tracked source is RHO_WATER * dh * area, once, on the one cell")
      if (allocated(error)) return
      call check(error, sf%lrunoff(ip, jp) == RATE .and. sum(sf%lrunoff) == RATE, &
                 "the reader landed the file's (3,2) at (nghost+3, nghost+2), ghosts untouched")
      call sf%destroy(); call ms%destroy(); call reader%destroy(); call ss%destroy()
   end subroutine test_point_source_exact

   ! -----------------------------------------------------------------
   ! 2. calving latent heat
   ! -----------------------------------------------------------------

   subroutine test_calving_latent_heat(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FN = "fwb_runoff_calving.nc"
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_surface_stress_t) :: ss
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_dataovr_config_t) :: cfg
      real(wp), allocatable :: at(:, :)
      real(wp), allocatable :: t_ini(:)
      real(wp) :: dtemp, dh
      integer :: i, j, k, ip, jp, ierr

      call write_point_file(FN, NXP, NYP, PI_FILE, PJ_FILE, RATE)
      cfg%enable = .true.
      cfg%time_mode = "linear"
      cfg%frunoff%file = FN
      cfg%frunoff%var = "field"
      call setup_unit(grid, ms, sf, ss, at, reader, df, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. df%id_frunoff > 0, "frunoff registers")
      if (allocated(error)) return
      ip = NGHOST + PI_FILE
      jp = NGHOST + PJ_FILE
      t_ini = ms%tracers(ms%idx_temperature)%c(ip, jp, :)
      call map_all(ms, sf, at, reader)
      call one_step(grid, ms, sf, ss, at, reader, df, 0.0_wp)
      call pull_state(ms, sf)
      call unmap_all(ms, sf, at, reader)

      call check(error, sf%q_fus_rof(ip, jp) == -ICE_LAT_FUS*RATE, "q_fus_rof = -L_f*frunoff at the source")
      if (allocated(error)) return
      call check(error, count(sf%q_fus_rof /= 0.0_wp) == 1, "q_fus_rof is non-zero at the source cell only")
      if (allocated(error)) return
      call check(error, sf%Q_heat(ip, jp) == -ICE_LAT_FUS*RATE, "the assembler adds q_fus_rof into Q_heat (nothing else is on)")
      if (allocated(error)) return

      ! The deposit sees the OLD thickness (it runs before the mass step).
      dtemp = -ICE_LAT_FUS*RATE*DT/(sf%rho0*sf%cp*H_LAYER)
      call check(error, abs((ms%tracers(ms%idx_temperature)%c(ip, jp, NZ) - t_ini(NZ)) - dtemp) < 1.0e-12_wp*abs(dtemp), &
                 "dT = -L_f frunoff dt/(rho0 cp h) on the top layer of the source column")
      if (allocated(error)) return
      call check(error, dtemp < 0.0_wp, "non-vacuous: calving cools")
      if (allocated(error)) return
      dh = RATE*DT/sf%rho0
      call check(error, abs(ms%h_layer(ip, jp, NZ) - (H_LAYER + dh)) < 1.0e-14_wp*H_LAYER, &
                 "the calved mass also joins as volume")
      if (allocated(error)) return
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            do k = 1, NZ
               if (i == ip .and. j == jp .and. k == NZ) cycle
               call check(error, ms%tracers(ms%idx_temperature)%c(i, j, k) == T0 + 0.5_wp*real(k, wp), &
                          "every other layer and column keeps its exact T bits")
               if (allocated(error)) return
            end do
         end do
      end do
      call sf%destroy(); call ms%destroy(); call reader%destroy(); call ss%destroy()
   end subroutine test_calving_latent_heat

   ! -----------------------------------------------------------------
   ! 3. runoff on land is fatal
   ! -----------------------------------------------------------------

   subroutine test_runoff_on_land_fatal(error)
      !! The verdict and the counting kernel (an abort is not run).  Owned cells only: a land cell in
      !! the halo is its neighbour's to count, and a wet cell never counts.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = NXP + 2*NGHOST, NY = NYP + 2*NGHOST
      real(wp), allocatable :: wet(:, :), rn(:, :)
      integer :: n

      call check(error,.not. runoff_on_land_is_fatal(0_int64), "no land runoff is fine")
      if (allocated(error)) return
      call check(error, runoff_on_land_is_fatal(1_int64), "one land cell with runoff is fatal")
      if (allocated(error)) return

      allocate (wet(NX, NY), rn(NX, NY))
      wet = 1.0_wp
      rn = 0.0_wp
      wet(NGHOST + 3, NGHOST + 2) = 0.0_wp
      wet(1, 1) = 0.0_wp
      !$acc enter data copyin(wet, rn)

      cases: block
         call runoff_land_count_impl(NX, NY, NGHOST, wet, rn, n)
         call check(error, n == 0, "a land cell with zero runoff does not count")
         if (allocated(error)) exit cases

         rn(NGHOST + 3, NGHOST + 2) = 1.0e-9_wp
         !$acc update device(rn)
         call runoff_land_count_impl(NX, NY, NGHOST, wet, rn, n)
         call check(error, n == 1, "a tiny positive runoff on an owned land cell counts")
         if (allocated(error)) exit cases

         rn(1, 1) = 5.0_wp
         rn(NGHOST + 1, NGHOST + 1) = 7.0_wp
         !$acc update device(rn)
         call runoff_land_count_impl(NX, NY, NGHOST, wet, rn, n)
         call check(error, n == 1, "a ghost land cell and a WET cell with runoff do not count")
         if (allocated(error)) exit cases

         wet(NGHOST + 2, NGHOST + 4) = 0.0_wp
         rn(NGHOST + 2, NGHOST + 4) = ieee_value(1.0_wp, ieee_quiet_nan)
         !$acc update device(wet, rn)
         call runoff_land_count_impl(NX, NY, NGHOST, wet, rn, n)
         call check(error, n == 2, "a NaN on an owned land cell counts")
      end block cases
      !$acc exit data delete(wet, rn)
   end subroutine test_runoff_on_land_fatal

   ! -----------------------------------------------------------------
   ! 4. cyclic time
   ! -----------------------------------------------------------------

   subroutine test_cyclic_interp(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FN = "fwb_runoff_cyclic.nc"
      real(wp), parameter :: PERIOD = 300.0_wp
      real(wp), parameter :: A(3) = [1.0e-3_wp, 2.0e-3_wp, 4.0e-3_wp]
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(ocean_surface_stress_t) :: ss
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: at(:, :)
      real(wp) :: f(NXP, NYP, 3), tq(4), want(4)
      integer :: n, ierr

      do n = 1, 3
         f(:, :, n) = A(n)
      end do
      call write_field_file(FN, NXP, NYP, 3, [0.0_wp, 100.0_wp, 200.0_wp], f)
      cfg%enable = .true.
      cfg%time_mode = "cyclic"
      cfg%cycle_period = PERIOD
      cfg%lrunoff%file = FN
      cfg%lrunoff%var = "field"
      call setup_unit(grid, ms, sf, ss, at, reader, df, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. df%id_lrunoff > 0, "cyclic lrunoff registers")
      if (allocated(error)) return
      call map_all(ms, sf, at, reader)

      ! mid-record (0..100), mid-record (100..200), the WRAP bracket (200 -> record 1 at 300),
      ! and one full period later.
      tq = [50.0_wp, 150.0_wp, 250.0_wp, 350.0_wp]
      want = [0.5_wp*(A(1) + A(2)), 0.5_wp*(A(2) + A(3)), 0.5_wp*(A(3) + A(1)), 0.5_wp*(A(1) + A(2))]
      do n = 1, 4
         call ocean_data_input_update_all(reader, tq(n))
         call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, tq(n), wet_mask=ms%wet_mask)
         !$acc update self(sf%lrunoff)
         call check(error, maxval(abs(sf%lrunoff(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP) - want(n))) < 1.0e-15_wp, &
                    "cyclic linear blend (mid-record, year wrap, next period)")
         if (allocated(error)) exit
      end do
      call unmap_all(ms, sf, at, reader)
      call sf%destroy(); call ms%destroy(); call reader%destroy(); call ss%destroy()
   end subroutine test_cyclic_interp

   ! -----------------------------------------------------------------
   ! 5. configure-time gates
   ! -----------------------------------------------------------------

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function basin_nml(freshwater, extra, tail) result(nml)
      !! A closed 12 x 8 x 4 basin, `ssp_rk2`, sigma, components on, with `freshwater` as given;
      !! `extra` is appended to `&ocean_forcing_nml`'s group as further groups.
      character(len=*), intent(in) :: freshwater, extra, tail
      character(len=:), allocatable :: nml

      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 12, ny = 8, nghost = 3, dx = 4000.0, dy = 4000.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 400.0, wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_bottom = 6.0, T_init_surface = 8.0 /"//NL// &
            "&ocean_eos_nml eos = 'linear' /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&ocean_bt_nml split_scheme = 'ssp_rk2', auto_n_inner = .true. /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true., freshwater = '"//freshwater//"' /"//NL// &
            "&ocean_hvisc_nml nu_h = 100.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL//extra//tail
   end function basin_nml

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

   subroutine test_config_gates(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: REFUSED = OCEAN_STATUS_ERR_CONFIG_VALIDATE
      character(len=*), parameter :: FN = "fwb_runoff_cfg.nc"
      character(len=*), parameter :: LR = "&ocean_dataovr_nml enable = .true., time_mode = 'static', "// &
                                     "lrunoff_file = '"//FN//"', lrunoff_var = 'field' /"//NL
      character(len=*), parameter :: FR = "&ocean_dataovr_nml enable = .true., time_mode = 'static', "// &
                                     "frunoff_file = '"//FN//"', frunoff_var = 'field' /"//NL
      type(config_t) :: cfg

      call ensure_comm()
      call write_point_file(FN, 12, 8, 5, 3, RATE)
      call check(error, nml_ierr(basin_nml("mass", LR, "")) == OCEAN_STATUS_OK, "accepted: lrunoff alone under 'mass'")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("mass", FR, "")) == OCEAN_STATUS_OK, "accepted: frunoff alone under 'mass'")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("virtual", LR, "")) == OCEAN_STATUS_OK, &
                 "accepted: lrunoff under 'virtual' (it enters the virtual salt flux)")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("none", LR, "")) == REFUSED, "refused: lrunoff under freshwater = 'none' (inert)")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("none", FR, "")) == REFUSED, "refused: frunoff under freshwater = 'none' (inert)")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("mass", &
                                           "&ocean_dataovr_nml enable = .true., time_mode = 'static', "// &
                                           "lrunoff_file = '"//FN//"' /"//NL, "")) == REFUSED, &
                 "refused: a runoff file with no variable name")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("mass", &
                                           "&ocean_dataovr_nml enable = .true., time_mode = 'static', "// &
                                           "lrunoff_file = 'fwb_no_such_file.nc', lrunoff_var = 'field' /"//NL, "")) == REFUSED, &
                 "refused: a runoff file that does not exist")
      if (allocated(error)) return
      call check(error, nml_ierr(basin_nml("virtual", LR, &
                                           "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
                                           "draft_depth = 100.0, draft_slope = -1.0e-3, draft_x0 = 0.0, "// &
                                           "h_min_cavity = 50.0 /"//NL)) == REFUSED, &
                 "refused: runoff with an ice-shelf cavity under 'virtual' (the cover twin would drop it)")
      if (allocated(error)) return

      ! Not an atmosphere-forcing prerequisite: the runoff tags are exempt from the
      ! `&ocean_atm_forcing_nml` / `&ocean_dataovr_nml` conflict, tau/heat/evap/lprec are not.
      cfg%start_date = "1958-01-01"
      cfg%ocean%forcing%enable_components = .true.
      cfg%ocean%grid%grid_config = "spherical"
      cfg%ocean%atm_forcing%filelist_tas = "a"
      cfg%ocean%atm_forcing%filelist_huss = "a"
      cfg%ocean%atm_forcing%filelist_uas = "a"
      cfg%ocean%atm_forcing%filelist_vas = "a"
      cfg%ocean%atm_forcing%filelist_psl = "a"
      cfg%ocean%atm_forcing%filelist_rlds = "a"
      cfg%ocean%atm_forcing%filelist_rsds = "a"
      cfg%ocean%atm_forcing%filelist_prra = "a"
      cfg%ocean%atm_forcing%filelist_prsn = "a"
      cfg%ocean%dataovr%lrunoff%file = FN
      cfg%ocean%dataovr%frunoff%file = FN
      call check(error, atm_forcing_prereq_code(cfg) /= ATM_PREREQ_DATAOVR, &
                 "runoff tags do not conflict with the atmosphere forcing")
      if (allocated(error)) return
      call check(error, atm_forcing_prereq_code(cfg) == ATM_PREREQ_BULK_FLUX, &
                 "non-vacuous: the prerequisite chain got past the dataovr check")
      if (allocated(error)) return
      cfg%ocean%dataovr%lprec%file = FN
      call check(error, atm_forcing_prereq_code(cfg) == ATM_PREREQ_DATAOVR, "an lprec file still conflicts")
   end subroutine test_config_gates

   ! -----------------------------------------------------------------
   ! 6. engine path
   ! -----------------------------------------------------------------

   subroutine test_engine_closed_basin(error)
      !! 12 x 8 closed basin, ssp_rk2, `freshwater = "mass"`, a point source of liquid runoff from a
      !! static file.  The basin's mass rises by exactly `(N - 1) * R * dt * area * RHO_WATER/rho0`
      !! (the snapshot is one step old) and the console's mass residual `(M - M0) + out - src` closes.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FN = "fwb_runoff_engine.nc"
      integer, parameter :: NSTEP = 6, IP = 5, JP = 3
      real(wp), parameter :: DTE = 300.0_wp
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      type(efp_t) :: e
      real(wp) :: m0, m1, expect, src, resid, t
      integer :: n, ierr

      call ensure_comm()
      call write_point_file(FN, 12, 8, IP, JP, RATE)
      call read_config_from_string(basin_nml("mass", "&ocean_dataovr_nml enable = .true., time_mode = 'static', "// &
                                             "lrunoff_file = '"//FN//"', lrunoff_var = 'field' /"//NL, ""), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the namelist parses")
      if (allocated(error)) return
      call validate_config(cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the namelist validates")
      if (allocated(error)) return
      call engine_setup(engine, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the engine sets up with the runoff reader")
      if (allocated(error)) return
      call engine_enter_data(engine, cfg)

      gate: block
         m0 = mass_now(engine)
         t = 0.0_wp
         do n = 1, NSTEP
            call engine_step(engine, DTE, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
            if (allocated(error)) exit gate
            call engine_step_ice(engine, cfg, DTE, t, ierr=ierr)
            call engine_step_finalize(engine, DTE, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine_step_finalize")
            if (allocated(error)) exit gate
            t = t + DTE
         end do
         m1 = mass_now(engine)
         ! The thermo step spends the snapshot assembled at the END of the previous step
         ! (`engine_step_finalize`), so the first step spends nothing: NSTEP - 1 applications.
         expect = real(NSTEP - 1, wp)*RATE*DTE*engine%state%metrics%areaT(engine%grid%nghost + IP, engine%grid%nghost + JP) &
                  *RHO_WATER/engine%state%surface_flux%rho0
         ! (the tolerance is the round-off of the two ~1e15 kg totals, not of their difference)
         call check(error, abs((m1 - m0) - expect) < 1.0e-14_wp*m0, "the basin gains exactly the runoff mass")
         if (allocated(error)) exit gate
         e%v = engine%state%multilayer%mass_src_efp
         e%poison = engine%state%multilayer%mass_src_efp_poison
         src = engine%state%multilayer%mass_src + real(efp_to_real(e), wp)
         resid = ((m1 - m0) + engine%state%multilayer%mass_out - src)/m0
         call check(error, abs(resid) < 1.0e-13_wp, "the console mass residual closes")
         if (allocated(error)) exit gate
         call check(error, abs(src) > 0.0_wp, "non-vacuous: the runoff booked a mass source")
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine test_engine_closed_basin

   function mass_now(engine) result(m)
      type(ocean_engine_t), intent(inout) :: engine
      real(wp) :: m

      associate (ms => engine%state%multilayer, mt => engine%state%metrics, ng => engine%grid%nghost)
         !$acc update self(ms%h_layer)
         m = compute_total_h(ms%h_layer, mt%areaT, ng)*RHO_WATER
      end associate
   end function mass_now

end module test_ocean_runoff
