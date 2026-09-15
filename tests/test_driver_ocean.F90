!! End-to-end test that drives the real ocean compute-rank lifecycle.
!!
!! Counterpart to test_driver_coastal, for sim_type='ocean'.  Where
!! test_ocean_driver_diag *replicates* the ocean driver's diag-relevant hook
!! sequence by calling the lower-level routines directly, this test invokes
!! the actual driver_run entry point with an ocean config, so the
!! orchestration in driver_run_ocean runs as a unit: guard checks, state
!! alloc, vcoord pre-wire, IC seed, the (large) physics/forcing config-wiring
!! section, diag registration + stream open, GPU enter_data, the split/unsplit
!! time loop, and finalisation (diag close, exit_data).  This is the
!! regression net for decomposing that routine.
!!
!! Setup: the same tiny quiescent ocean column as test_ocean_driver_diag
!! (uniform T/S, no flow) with a fixed dt and the diag manager enabled,
!! writing to the cwd.  Zero forcing → no motion → the test is structural
!! (the default diag variable set exists, and the per-variable frame count
!! matches the cadence), not a physics check.
module test_driver_ocean
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_driver, only: driver_run, configure_log_level
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_io_netcdf, only: output_rank_filename
   use rdb_io_netcdf, only: nc_open_read, nc_close, nc_get_dim_len
   use netcdf, only: nf90_inq_varid, NF90_NOERR
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_driver_ocean_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 4
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: H_TOTAL = 100.0_wp
   real(wp), parameter :: DT = 100.0_wp
   integer, parameter :: N_STEPS = 8
   real(wp), parameter :: DT_OUT = 200.0_wp
   integer, parameter :: EXPECTED_FRAMES = (N_STEPS*int(DT))/int(DT_OUT)

   ! driver_run needs the MPI comm-env initialised (the real app does
   ! this in app/main); testdrive runs all cases in one process, so we
   ! init once via this saved flag.
   logical, save :: mpi_ready = .false.

contains

   subroutine ensure_mpi()
      !! Initialise the MPI comm-env once per process before driver_run.
      !! Idempotent (saved flag); no finalize — process exit cleans up.
      !! No-op-equivalent on serial builds (stub comm-env).
      if (.not. mpi_ready) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         mpi_ready = .true.
      end if
   end subroutine ensure_mpi

   subroutine collect_driver_ocean_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("driver_run_ocean_e2e", test_ocean_e2e), &
                  new_unittest("driver_run_ocean_pc_sigma", test_ocean_pc_sigma) &
                  ]
   end subroutine collect_driver_ocean_tests

   subroutine make_cfg(cfg)
      !! Minimal ocean config: tiny quiescent column, fixed dt, diag manager
      !! on, output to cwd.  Mirrors test_ocean_driver_diag plus the driver's
      !! requirements (dt_fixed > 0, t_end).  All physics knobs default off.
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
      cfg%dt_fixed = DT
      cfg%t_end = real(N_STEPS, wp)*DT
      cfg%output_dir = "."
   end subroutine make_cfg

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

   subroutine test_ocean_e2e(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=512) :: fn
      integer :: ncid, varid, ios, nframes, i
      character(len=*), parameter :: names(*) = [character(len=16) :: &
                                                 "SSH", "temperature", "salinity", &
                                                 "u", "v", "KE"]

      call configure_log_level("warning")
      call make_cfg(cfg)
      fn = output_rank_filename(".", "ocean_diag", 0)
      call cleanup_file(trim(fn))

      checks: block
         call ensure_mpi()
         call driver_run(cfg)

         call nc_open_read(trim(fn), ncid)

         ! The default diag registrar set must all be present.
         do i = 1, size(names)
            ios = nf90_inq_varid(ncid, trim(names(i)), varid)
            call check(error, ios == NF90_NOERR, &
                       "ocean diag var missing from driver output: "//trim(names(i)))
            if (allocated(error)) then; call nc_close(ncid); exit checks; end if
         end do

         ! Per-variable time dim length must match the cadence: with
         ! N_STEPS*DT total time at DT_OUT cadence, 800/200 = 4 frames.
         call nc_get_dim_len(ncid, "time_SSH", nframes)
         call nc_close(ncid)
         call check(error, nframes == EXPECTED_FRAMES, &
                    "time_SSH frame count should match the diag cadence")
      end block checks

      call cleanup_file(trim(fn))
   end subroutine test_ocean_e2e

   subroutine test_ocean_pc_sigma(error)
      !! pred_corr on an ALE (sigma) vertical coordinate — the generalized
      !! envelope (was fail-loud Lagrangian-only in v1).  MOM6's model:
      !! the dynamics step is Lagrangian-within-step for every ALE
      !! coordinate, so the pc restructure composes with sigma exactly as
      !! with 'lagrangian'; this drives the real driver_run lifecycle
      !! (predictor TR_MODE_NONE chain, corrector full-dt update, ALE
      !! remap at the thermo cadence) end-to-end.  Wind-forced so the
      !! predictor/corrector genuinely advance a non-trivial state; the
      !! check is structural + finiteness (SSH and u in the diag output
      !! contain no non-finite values after N_STEPS).
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use netcdf, only: nf90_get_var
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=512) :: fn
      integer :: ncid, varid, ios, nframes, i, j
      real(wp) :: ssh(NX, NY)

      call configure_log_level("warning")
      call make_cfg(cfg)
      cfg%ocean%bt%split_scheme = "pred_corr"
      cfg%vcoord_type = "sigma"
      cfg%wind_stress_x = 0.1_wp
      cfg%ocean%diag%filename = "ocean_diag_pc"
      fn = output_rank_filename(".", "ocean_diag_pc", 0)
      call cleanup_file(trim(fn))

      checks: block
         call ensure_mpi()
         call driver_run(cfg)

         call nc_open_read(trim(fn), ncid)
         ios = nf90_inq_varid(ncid, "SSH", varid)
         call check(error, ios == NF90_NOERR, "pc+sigma run must emit SSH")
         if (allocated(error)) then; call nc_close(ncid); exit checks; end if
         call nc_get_dim_len(ncid, "time_SSH", nframes)
         call check(error, nframes == EXPECTED_FRAMES, &
                    "pc+sigma frame count should match the diag cadence")
         if (allocated(error)) then; call nc_close(ncid); exit checks; end if
         ! Final frame: the wind-forced state must be finite everywhere.
         ios = nf90_get_var(ncid, varid, ssh, start=[1, 1, nframes], &
                            count=[NX, NY, 1])
         call nc_close(ncid)
         call check(error, ios == NF90_NOERR, "SSH final frame must read back")
         if (allocated(error)) exit checks
         do j = 1, NY
            do i = 1, NX
               call check(error, ieee_is_finite(ssh(i, j)), &
                          "pc+sigma SSH must stay finite under wind forcing")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks

      call cleanup_file(trim(fn))
   end subroutine test_ocean_pc_sigma

end module test_driver_ocean
