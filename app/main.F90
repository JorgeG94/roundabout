!! Roundabout — GPU-native shallow-water / multilayer / regional ocean solver.
program rdb
   !! Thin entry-point wrapper.  Owns:
   !!   * MPI lifecycle (init / role assignment / finalize),
   !!   * command-line argument parsing + config load,
   !!   * output-directory bootstrap,
   !!   * dispatch between I/O-server and compute-rank paths.
   !!
   !! All compute-side logic lives in `rdb_driver%driver_run` so the
   !! same lifecycle is exercised from both `app/main.F90` and from
   !! unit tests (`tests/unit/test_driver_smoke.F90`).  This catches
   !! integration-level bugs (OpenACC partial-presence, init ordering,
   !! missing attaches) that per-kernel unit tests miss.
   use rdb_config, only: config_t, read_config, validate_config
   use rdb_nml_schema, only: nml_schema_t
   use rdb_driver, only: driver_run, configure_log_level
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_abort, &
                           comm_env_compute_size
   use pic_logger, only: logger => global_logger
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=:), allocatable :: doc_dir
   character(len=256) :: input_file
   integer :: arg_len, io_stat
   integer :: mpi_rank
   logical :: input_exists

   ! Phase 1: Basic MPI init (rank/size only, no GPU binding yet)
   call comm_env_init()
   mpi_rank = comm_env_rank()

   ! Read input file from command line (all ranks read independently).
   ! Running without a namelist used to fall back to a hard-coded default
   ! state, which then SIGBUS'd downstream once any output path was hit.
   ! Fail loudly here instead so the user gets a clear "no nml" message
   ! rather than a stale-default crash several seconds in.
   if (command_argument_count() < 1) then
      if (mpi_rank == 0) then
         call logger%error("FATAL: no input file provided.")
         call logger%error("Usage: rdb <input_file.nml>")
      end if
      call comm_env_abort(2)
   end if

   call get_command_argument(1, input_file, arg_len, io_stat)
   inquire (file=trim(input_file), exist=input_exists)
   if (.not. input_exists) then
      if (mpi_rank == 0) then
         call logger%error("FATAL: input file not found: "//trim(input_file))
      end if
      call comm_env_abort(2)
   end if
   ! read_config builds + strict-parses the schema internally (so any
   ! unknown-key / range / enum error error-stops here, with the full
   ! report); `schema` is returned for the rank-0 parameter-doc dumps.
   call read_config(trim(input_file), cfg, schema)
   call validate_config(cfg)
   ! Apply log level from config
   call configure_log_level(cfg%log_level)
   if (mpi_rank == 0) then
      call logger%info("Configuration read from: "//trim(input_file))
   end if

   ! Ensure output directory exists. Fortran's NetCDF/open calls don't create
   ! parent directories — mkdir -p is idempotent and safe to call from every
   ! rank on a shared filesystem.
   if (len_trim(cfg%output_dir) > 0) then
      call execute_command_line('mkdir -p "'//trim(cfg%output_dir)//'"', &
                                wait=.true., exitstat=io_stat)
      if (io_stat /= 0 .and. mpi_rank == 0) then
         call logger%warning("Could not create output directory: "// &
                             trim(cfg%output_dir))
      end if
   end if

   ! Rank 0 emits the MOM6-style parameter-documentation dumps from the
   ! schema returned by read_config (already built + strict-parsed).
   if (mpi_rank == 0) then
      if (len_trim(cfg%output_dir) > 0) then
         doc_dir = trim(cfg%output_dir)//"/"
      else
         doc_dir = ""
      end if
      call schema%write_doc_all(doc_dir//"rdb_parameter_doc.all")
      call schema%write_doc_short(doc_dir//"rdb_parameter_doc.short")
   end if

   ! Phase 2: Assign roles, bind GPU.  The dedicated I/O-server rank was a
   ! coastal-path feature; the ocean core writes its own per-rank diag/restart
   ! streams, so every rank is a compute rank.
   call comm_env_setup_roles(.false.)

   ! Guard: GPU + MPI + multiple compute ranks without CUDA-aware MPI
   ! is a performance trap — host-staged halo exchange destroys scaling.
#if defined(RDB_GPU_OFFLOAD) && !defined(RDB_CUDA_AWARE_MPI)
   if (comm_env_compute_size() > 1) then
      if (mpi_rank == 0) then
         call logger%error("FATAL: multi-GPU MPI build without CUDA-aware MPI.")
         call logger%error("Recompile with -DRDB_CUDA_AWARE_MPI=ON")
      end if
      call comm_env_abort(99)
   end if
#endif

   call driver_run(cfg)

   call comm_env_finalize()

end program rdb
