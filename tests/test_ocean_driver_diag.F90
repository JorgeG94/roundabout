!! End-to-end driver hook test for the ocean diag manager.
!!
!! Exercises the same sequence the driver runs:
!!   ocean_state%init  →  register_default_tracers  →
!!   ocean_state_seed_from_cfg  →  register_default_diags  →
!!   open_stream  →  N × (ocean_dyn_step + diag%step)  →  close_stream
!!
!! Then re-opens the NetCDF file and verifies the six default variables
!! (SSH, temperature, salinity, u, v, KE) are present and that the
!! per-variable time dim equals the expected fire count under the
!! cadence.  Zero-forcing, uniform IC → no flow → kernels touch state
!! but the test is structural (variable existence + frame counts), not
!! a physics check.
module test_ocean_driver_diag
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, &
                              ocean_state_exit_data, ocean_state_seed_from_cfg
   use rdb_state, only: register_default_tracers
   use rdb_ocean_dyn, only: ocean_dyn_step
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_diag_fills, only: register_default_diags
   use rdb_ocean_diag_netcdf, only: open_stream, close_stream
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, nc_get_dim_len, &
                            nc_get_varid
   use netcdf, only: nf90_inq_varid, NF90_NOERR
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_driver_diag_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 4
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: H_TOTAL = 100.0_wp
   real(wp), parameter :: DT = 100.0_wp
   integer, parameter :: N_STEPS = 8
   real(wp), parameter :: DT_OUT = 200.0_wp

contains

   subroutine collect_ocean_driver_diag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("driver_emits_default_diag_set", test_emits_default_set), &
                  new_unittest("driver_emits_expected_frame_count", test_frame_count) &
                  ]
   end subroutine collect_ocean_driver_diag_tests

   subroutine cleanup_file(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=filename, exist=exists)
      if (exists) then
         open (newunit=u, file=filename, status='old')
         close (u, status='delete')
      end if
   end subroutine cleanup_file

   subroutine make_cfg(cfg)
      type(config_t), intent(out) :: cfg
      cfg%sim_type = "ocean"
      cfg%nx = NX
      cfg%ny = NY
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%ocean%topo%max_depth = H_TOTAL
      cfg%initial_temperature = 15.0_wp
      cfg%initial_salinity = 35.0_wp
      cfg%ocean%diag%enabled = .true.
      cfg%ocean%diag%dt_out = DT_OUT
   end subroutine make_cfg

   subroutine run_driver_hook_sequence(cfg, filename, state, sf)
      !! Mirror of the diag-relevant init + time-loop hooks in
      !! `driver_run_ocean`.  Leaves `state` and `sf` alive for the
      !! caller's teardown so individual tests can do per-case assertions.
      type(config_t), intent(in) :: cfg
      character(len=*), intent(in) :: filename
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(inout) :: sf

      type(hgrid_t) :: grid
      real(wp) :: t
      integer :: step

      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      ! Fill the metrics slot (cartesian) so the PGF's idxCu/idyCv reads
      ! are non-zero — the driver does this in configure_ocean_metrics;
      ! this test bypasses the driver.
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call ocean_state_seed_from_cfg(state, grid, cfg)
      call register_default_tracers( &
         state%multilayer%tracers(state%multilayer%idx_salinity), &
         state%multilayer%tracers(state%multilayer%idx_temperature), &
         cfg)
      call register_default_diags(state, dt_out=cfg%ocean%diag%dt_out)
      call open_stream(state%diag, filename)

      call sf%init(grid)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      call ocean_state_enter_data(state)
      !$acc enter data copyin(sf)
      call sf%enter_data()

      t = 0.0_wp
      do step = 1, N_STEPS
         call ocean_dyn_step(grid, state%metrics, state%dyn, state%eos, &
                             state%coriolis_adv, state%continuity, &
                             state%pressure_force, state%hvisc, &
                             state%bdrag, state%surface_stress, &
                             state%vert_advect, state%hdiff_tracer, &
                             state%vdiff, state%vmix, &
                             state%multilayer, DT, sf=sf)
         t = t + DT
         call state%diag%step(state, DT, t)
      end do

      call sf%exit_data()
      !$acc exit data delete(sf)
      call ocean_state_exit_data(state)
      call close_stream(state%diag)
   end subroutine run_driver_hook_sequence

   subroutine test_emits_default_set(error)
      !! Run the driver hook sequence, then re-open the NetCDF and
      !! verify that every variable in the default registrar set
      !! (SSH, temperature, salinity, u, v, KE) has a varid.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(ocean_surface_flux_t) :: sf
      character(len=*), parameter :: FN = "test_ocean_driver_diag_set.nc"
      integer :: ncid, varid, ios
      character(len=*), parameter :: names(*) = [character(len=16) :: &
                                                 "SSH", "temperature", "salinity", &
                                                 "u", "v", "KE"]
      integer :: i

      call cleanup_file(FN)
      call make_cfg(cfg)

      checks: block
         call run_driver_hook_sequence(cfg, FN, state, sf)

         call nc_open_read(FN, ncid)
         do i = 1, size(names)
            ios = nf90_inq_varid(ncid, trim(names(i)), varid)
            call check(error, ios == NF90_NOERR, &
                       "default diag var missing from NetCDF: "//trim(names(i)))
            if (allocated(error)) exit checks
         end do
         call nc_close(ncid)
      end block checks

      call sf%destroy()
      call state%destroy()
      call cleanup_file(FN)
   end subroutine test_emits_default_set

   subroutine test_frame_count(error)
      !! Verify the per-variable time-dim length matches what the
      !! cadence dictates: with N_STEPS*DT total time and DT_OUT
      !! cadence, the time op fires `floor(N_STEPS*DT / DT_OUT)`
      !! times (MEAN flushes on cadence boundary, INSTANT on the same
      !! boundary).  Concrete: 8*100 / 200 = 4 frames per variable.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(ocean_surface_flux_t) :: sf
      character(len=*), parameter :: FN = "test_ocean_driver_diag_frames.nc"
      integer :: ncid, nfires
      integer, parameter :: EXPECTED_FRAMES = (N_STEPS*int(DT))/int(DT_OUT)

      call cleanup_file(FN)
      call make_cfg(cfg)

      checks: block
         call run_driver_hook_sequence(cfg, FN, state, sf)

         call nc_open_read(FN, ncid)
         call nc_get_dim_len(ncid, "time_SSH", nfires)
         call check(error, nfires == EXPECTED_FRAMES, &
                    "time_SSH frame count should match cadence")
         if (allocated(error)) exit checks
         call nc_get_dim_len(ncid, "time_temperature", nfires)
         call check(error, nfires == EXPECTED_FRAMES, &
                    "time_temperature frame count should match cadence")
         if (allocated(error)) exit checks
         call nc_get_dim_len(ncid, "time_KE", nfires)
         call check(error, nfires == EXPECTED_FRAMES, &
                    "time_KE frame count should match cadence")
         call nc_close(ncid)
      end block checks

      call sf%destroy()
      call state%destroy()
      call cleanup_file(FN)
   end subroutine test_frame_count

end module test_ocean_driver_diag
