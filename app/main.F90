!! Roundabout — GPU-native shallow-water / multilayer / regional ocean solver.
program rdb
   !! Thin entry-point wrapper.  Owns:
   !!   * MPI lifecycle (init / role assignment / finalize),
   !!   * command-line argument parsing + config load (and the
   !!     `--validate-only` configure-check-and-exit mode),
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
   use rdb_driver, only: driver_run, driver_validate, configure_log_level
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_abort, &
                           comm_env_compute_size
   use pic_logger, only: logger => global_logger
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=:), allocatable :: doc_dir
   character(len=256) :: input_file, arg
   integer :: arg_len, io_stat, i_arg, n_files
   integer :: mpi_rank, validate_ierr
   logical :: input_exists, validate_only

   ! Phase 1: Basic MPI init (rank/size only, no GPU binding yet)
   call comm_env_init()
   mpi_rank = comm_env_rank()

   ! Read input file from command line (all ranks read independently).
   ! Running without a namelist used to fall back to a hard-coded default
   ! state, which then SIGBUS'd downstream once any output path was hit.
   ! Fail loudly here instead so the user gets a clear "no nml" message
   ! rather than a stale-default crash several seconds in.
   !
   ! `--validate-only` (anywhere on the line) stops after the configure
   ! checks: parse + validate_config + the engine_setup configure stages,
   ! exit 0 if the configuration is accepted, 3 with every refusal reason
   ! logged if not (a stage that still `error stop`s instead of returning
   ! a status exits non-zero with its own message).  No device mapping, no
   ! stepping, no output files.  Run it with `&logging_nml log_level =
   ! "error"` and stdout carries exactly the refusal reasons.
   validate_only = .false.
   n_files = 0
   input_file = ""
   do i_arg = 1, command_argument_count()
      call get_command_argument(i_arg, arg, arg_len, io_stat)
      if (io_stat /= 0) then
         ! -1: the argument is longer than `arg` (truncated); > 0: it could
         ! not be read.  Either way the path below would be wrong, so stop
         ! here instead of reporting a misleading "file not found".
         if (mpi_rank == 0) then
            call logger%error("FATAL: command-line argument longer than 256 characters, "// &
                              "or unreadable: "//trim(arg))
         end if
         call comm_env_abort(2)
      end if
      if (trim(arg) == "--validate-only") then
         validate_only = .true.
      else if (arg(1:1) == "-") then
         if (mpi_rank == 0) then
            call logger%error("FATAL: unknown option: "//trim(arg))
            call logger%error("Usage: rdb [--validate-only] <input_file.nml>")
         end if
         call comm_env_abort(2)
      else
         n_files = n_files + 1
         input_file = arg
      end if
   end do
   if (n_files /= 1) then
      if (mpi_rank == 0) then
         if (n_files == 0) then
            call logger%error("FATAL: no input file provided.")
         else
            call logger%error("FATAL: more than one input file given.")
         end if
         call logger%error("Usage: rdb [--validate-only] <input_file.nml>")
      end if
      call comm_env_abort(2)
   end if

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
   if (validate_only) then
      ! Every validate_config failure is logged (it accumulates them all)
      ! before the rollup, so a refused configuration names each reason.
      call configure_log_level(cfg%log_level)
      call validate_config(cfg, validate_ierr)
      if (validate_ierr /= 0) then
         if (mpi_rank == 0) call logger%error("VALIDATE-ONLY: REFUSED by validate_config: "// &
                                              trim(input_file))
         ! Every rank evaluates the same checks and refuses together, so a
         ! clean finalize (not an abort) keeps the exit code: 3 = refused.
         call comm_env_finalize()
         error stop 3
      end if
      call comm_env_setup_roles(.false.)
      call driver_validate(cfg, validate_ierr)
      if (validate_ierr /= 0) then
         if (mpi_rank == 0) call logger%error("VALIDATE-ONLY: REFUSED by engine_setup: "// &
                                              trim(input_file))
         call comm_env_finalize()
         error stop 3
      end if
      if (mpi_rank == 0) call logger%info("VALIDATE-ONLY: ACCEPTED: "//trim(input_file))
      call comm_env_finalize()
      stop
   end if
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
