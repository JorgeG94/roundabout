!! Communication environment — the ONLY backend.
!!
!! Roundabout no longer carries a hand-written single-process stub beside an
!! MPI implementation.  This module always compiles and always talks to
!! `pic_mpi_lib`; whether that resolves to a real MPI library or to
!! pic-mpi's serial backend is pic-mpi's build-time choice
!! (`PIC_ENABLE_MPI`), invisible here.  See `src/comm/README` on the
!! single-rank contract every collective call site must honour.
module rdb_comm_env
   !! Wraps pic_mpi_lib for MPI init/finalize, rank/size queries,
   !! broadcast, abort, and per-node GPU binding.
   !!
   !! Two-phase initialisation:
   !!   1. comm_env_init()   — basic MPI, cache world rank/size
   !!   2. comm_env_setup_roles(use_io_server) — split comms, assign roles, bind GPU
   !!
   !! When use_io_server is enabled, the last rank on each node becomes
   !! a dedicated I/O server (no GPU, no solver). All solver/decomp/halo
   !! operations use compute_comm (excludes I/O ranks). When disabled,
   !! compute_comm = comm_world (current behavior, no overhead).

   use rdb_constants, only: wp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize, &
                          bcast, abort_comm
#ifdef RDB_HAS_OPENACC_RUNTIME
   ! On NVHPC and Cray, stdpar (do concurrent) dispatches through the
   ! OpenACC runtime regardless of whether -acc= is on, so the openacc
   ! module is always available and is needed to bind the stdpar device
   ! under multi-GPU MPI — even on the OpenMP backend (e.g. -mp=gpu on
   ! NVHPC or -h omp on Cray with -acc= off). Sentinel set by
   ! cmake/compiler_flags.cmake.
   use openacc, only: acc_get_device_num, acc_get_device_type, &
                      acc_set_device_num, acc_device_default, &
                      acc_get_num_devices
#endif
#ifdef _OPENMP
   use omp_lib, only: omp_set_default_device
#endif
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   public :: comm_env_init
   public :: comm_env_setup_roles
   public :: comm_env_finalize
   public :: comm_env_rank
   public :: comm_env_size
   public :: comm_env_bcast_real
   public :: comm_env_abort
   public :: comm_env_is_io_server
   public :: comm_env_io_server_rank
   public :: comm_env_compute_comm
   public :: comm_env_global_comm
   public :: comm_env_compute_rank
   public :: comm_env_compute_size
   public :: comm_env_node_compute_ranks
   public :: comm_env_node_n_compute

   type(comm_t), save :: comm_global
      !! Global MPI communicator (cached after init)
   type(comm_t), save :: comm_compute
      !! Compute-only communicator (excludes I/O ranks)
   logical, save :: env_initialised = .false.
      !! True between comm_env_init and comm_env_finalize.  Makes init
      !! idempotent and finalize a no-op when MPI was never initialised
      !! (e.g. a test binary that never inits the comm-env), so the shared
      !! per-test main can call finalize unconditionally at exit.
   integer, save :: cached_rank = -1
   integer, save :: cached_size = -1
   integer, save :: cached_compute_rank = -1
   integer, save :: cached_compute_size = -1
   logical, save :: cached_is_io_server = .false.
   integer, save :: cached_io_server_rank = -1
      !! World rank of this node's I/O server (-1 if none)
   integer, save :: cached_node_n_compute = 0
      !! Number of compute ranks on this node
   integer, save, allocatable :: cached_node_compute_ranks(:)
      !! World ranks of compute processes on this node

contains

   subroutine comm_env_init()
      !! Phase 1: Initialise MPI, cache world rank/size
      !! Call this before reading config.  Idempotent: a second call is a
      !! no-op, so multiple entry points / testdrive cases can call it.
      if (env_initialised) return
      call pic_mpi_init()
      comm_global = comm_world()
      cached_rank = comm_global%rank()
      cached_size = comm_global%size()

      ! Default: all ranks are compute, compute_comm = comm_world
      comm_compute = comm_global
      cached_compute_rank = cached_rank
      cached_compute_size = cached_size
      env_initialised = .true.
   end subroutine comm_env_init

   subroutine comm_env_setup_roles(use_io_server)
      !! Phase 2: Assign roles (I/O server vs compute), bind GPU
      !! Call this after reading config, before decomp/solver init.
      logical, intent(in) :: use_io_server
         !! Dedicate one rank per node as I/O server

      type(comm_t) :: node_comm
      integer :: node_rank, node_size, gpu_rank, i, n_dev
      integer, allocatable :: world_ranks(:)

      ! Compute node-local rank for GPU binding and I/O server assignment
      node_comm = comm_global%split()
      node_rank = node_comm%rank()
      node_size = node_comm%size()

      if (use_io_server .and. cached_size > 1 .and. node_size > 1) then
         ! Last rank on each node becomes the I/O server
         cached_is_io_server = (node_rank == node_size - 1)

         ! Broadcast I/O server's world rank to all ranks on this node
         cached_io_server_rank = cached_rank
         call bcast(node_comm, cached_io_server_rank, 1, node_size - 1)

         ! Gather world ranks of all node-local processes to build
         ! the compute rank list (exclude the I/O server)
         ! Use bcast from each rank as a simple allgather substitute
         allocate (world_ranks(0:node_size - 1))
         do i = 0, node_size - 1
            world_ranks(i) = cached_rank
            call bcast(node_comm, world_ranks(i), 1, i)
         end do

         cached_node_n_compute = node_size - 1
         allocate (cached_node_compute_ranks(cached_node_n_compute))
         do i = 0, node_size - 2
            cached_node_compute_ranks(i + 1) = world_ranks(i)
         end do
         deallocate (world_ranks)

         ! Create compute communicator (excludes I/O ranks)
         if (cached_is_io_server) then
            comm_compute = comm_global%split_by(1)  ! color 1 = I/O
         else
            comm_compute = comm_global%split_by(0)  ! color 0 = compute
         end if

         cached_compute_rank = comm_compute%rank()
         cached_compute_size = comm_compute%size()

         ! GPU binding: compute ranks bind to their node-local rank
         gpu_rank = node_rank
      else
         ! No I/O server: all ranks are compute
         cached_is_io_server = .false.
         cached_io_server_rank = -1
         comm_compute = comm_global
         cached_compute_rank = cached_rank
         cached_compute_size = cached_size
         cached_node_n_compute = node_size
         allocate (cached_node_compute_ranks(node_size))
         do i = 0, node_size - 1
            cached_node_compute_ranks(i + 1) = cached_rank - node_rank + i
         end do
         gpu_rank = node_rank
      end if

      call node_comm%finalize()

      ! Bind GPU (compute ranks only). Two API calls are needed when both
      ! `do concurrent` (stdpar) and `!$omp target` regions co-exist:
      !
      !   * acc_set_device_num — selects the OpenACC runtime device, which
      !     is what stdpar (NVHPC, Cray) uses to dispatch `do concurrent`
      !     loops. Active whenever the compiler exposes the OpenACC
      !     runtime — see RDB_HAS_OPENACC_RUNTIME in compiler_flags.
      !   * omp_set_default_device — selects the OpenMP target device,
      !     used by `!$omp target` regions on the dc-openmp backend.
      !     Active when -mp= / -fopenmp is on.
      !
      ! On the OpenACC backend (-acc=gpu) only the first fires; on the
      ! OpenMP backend both fire. Either way, every rank lands on the
      ! right GPU under mpirun. acc_device_default lets the runtime
      ! pick the appropriate device kind (NVIDIA, AMD, etc.).
      if (.not. cached_is_io_server) then
#ifdef RDB_HAS_OPENACC_RUNTIME
         ! Clamp the node-local GPU index to the number of devices this
         ! process can actually see. With all GPUs visible this is a no-op
         ! (mod(node_rank, ndev) == node_rank for node_rank < ndev). When the
         ! launcher pins one GPU per rank via CUDA_VISIBLE_DEVICES, only one
         ! device is visible, so this binds device 0 (the rank's own GPU) and
         ! the CUDA-aware-MPI primary context also lands there -- no leftover
         ! context on the global device 0. Also makes ranks > GPUs (over-
         ! subscription) bind round-robin instead of failing.
         !
         ! REVISIT: this clamp may be unnecessary -- NVHPC's acc_set_device_num
         ! might already tolerate an out-of-range device num under one-GPU
         ! CUDA_VISIBLE_DEVICES pinning (clamping internally). It's untested
         ! (the runs that confirmed the fix used this clamp). If a no-clamp
         ! build + pin reports a clean bind, drop this. Kept for now because
         ! it's a no-op when all GPUs are visible and the spec calls
         ! out-of-range device nums undefined, so it's the portable choice.
         n_dev = acc_get_num_devices(acc_device_default)
         if (n_dev > 0) gpu_rank = mod(gpu_rank, n_dev)
         call acc_set_device_num(gpu_rank, acc_device_default)
#endif
#ifdef _OPENMP
         call omp_set_default_device(gpu_rank)
#endif
         call print_gpu_binding(gpu_rank)
      end if

   end subroutine comm_env_setup_roles

   subroutine print_gpu_binding(requested_gpu)
      !! Per-rank one-shot diagnostic: world rank -> requested device num /
      !! actual device num reported by the OpenACC runtime.  Lets us
      !! confirm that mpirun is binding each rank to its own GPU rather
      !! than serialising N ranks on device 0.
      integer, intent(in) :: requested_gpu
      integer :: actual_dev

#ifdef RDB_HAS_OPENACC_RUNTIME
      ! Works on either backend: stdpar+OpenACC and -mp=gpu / -h omp
      ! both keep the OpenACC runtime live, and acc_get_device_num
      ! reports the same device that stdpar/-acc dispatches against.
      actual_dev = acc_get_device_num(acc_get_device_type())
#else
      actual_dev = -1
#endif
      write (output_unit, "(A,I0,A,I0,A,I0,A,I0,A,I0)") &
         "[gpu-bind] world_rank=", cached_rank, &
         "/", cached_size, &
         " compute_rank=", cached_compute_rank, &
         " requested_dev=", requested_gpu, &
         " actual_dev=", actual_dev
      if (cached_rank == 0) then
#ifdef RDB_CUDA_AWARE_MPI
         write (output_unit, "(A)") "[gpu-bind] halo path: GPU-direct (CUDA-aware MPI)"
#else
         write (output_unit, "(A)") "[gpu-bind] halo path: HOST-STAGED (rebuild with -DRDB_CUDA_AWARE_MPI=ON for GPU-direct)"
#endif
      end if
      flush (output_unit)
   end subroutine print_gpu_binding

   subroutine comm_env_finalize()
      !! Finalise MPI.  No-op if the comm-env was never initialised, so it
      !! is safe to call unconditionally at program exit -- the shared
      !! per-test main does this for the few tests that init MPI, and it
      !! costs nothing for the rest (and on serial builds via the stub).
      if (.not. env_initialised) return
      if (allocated(cached_node_compute_ranks)) deallocate (cached_node_compute_ranks)
      call comm_global%finalize()
      call pic_mpi_finalize()
      env_initialised = .false.
   end subroutine comm_env_finalize

   function comm_env_rank() result(rank)
      !! Return this process's world MPI rank
      integer :: rank
      rank = cached_rank
   end function comm_env_rank

   function comm_env_size() result(nprocs)
      !! Return total number of MPI processes (including I/O ranks)
      integer :: nprocs
      nprocs = cached_size
   end function comm_env_size

   function comm_env_compute_rank() result(rank)
      !! Return this process's rank in the compute communicator
      integer :: rank
      rank = cached_compute_rank
   end function comm_env_compute_rank

   function comm_env_compute_size() result(nprocs)
      !! Return number of compute ranks (excludes I/O ranks)
      integer :: nprocs
      nprocs = cached_compute_size
   end function comm_env_compute_size

   function comm_env_compute_comm() result(comm)
      !! Return the compute-only communicator
      !! Falls back to comm_world() if comm_env_init has not been called
      !! (e.g. tests that use raw MPI_Init).
      type(comm_t) :: comm
      if (cached_rank >= 0) then
         comm = comm_compute
      else
         comm = comm_world()
      end if
   end function comm_env_compute_comm

   function comm_env_global_comm() result(comm)
      !! Return the cached global (world) communicator.
      !!
      !! Use this instead of `pic_mpi_lib::comm_world()` from inside
      !! per-step code: `comm_world()` does an `MPI_Comm_dup` on every
      !! call, and that dup is collective on `MPI_COMM_WORLD`.  Calling
      !! it from a routine that runs on only some ranks (e.g. the I/O
      !! client send path, which compute ranks hit but the I/O server
      !! does not) deadlocks.  `comm_global` is duplicated exactly once
      !! by `comm_env_init` on every rank, so handing it out is free.
      type(comm_t) :: comm
      if (cached_rank >= 0) then
         comm = comm_global
      else
         comm = comm_world()
      end if
   end function comm_env_global_comm

   function comm_env_is_io_server() result(is_io)
      !! Return .true. if this rank is a dedicated I/O server
      logical :: is_io
      is_io = cached_is_io_server
   end function comm_env_is_io_server

   function comm_env_io_server_rank() result(rank)
      !! Return world rank of this node's I/O server (-1 if none)
      integer :: rank
      rank = cached_io_server_rank
   end function comm_env_io_server_rank

   function comm_env_node_n_compute() result(n)
      !! Return number of compute ranks on this node
      integer :: n
      n = cached_node_n_compute
   end function comm_env_node_n_compute

   function comm_env_node_compute_ranks() result(ranks)
      !! Return world ranks of compute processes on this node
      integer, allocatable :: ranks(:)
      ranks = cached_node_compute_ranks
   end function comm_env_node_compute_ranks

   subroutine comm_env_bcast_real(val)
      !! Broadcast a single real(wp) scalar from compute rank 0
      !! Uses compute communicator so I/O server ranks don't participate.
      real(wp), intent(inout) :: val
      call bcast(comm_compute, val, 1, 0)
   end subroutine comm_env_bcast_real

   subroutine comm_env_abort(code)
      !! Abort all MPI processes with given error code
      integer, intent(in) :: code
      call abort_comm(comm_global, code)
   end subroutine comm_env_abort

end module rdb_comm_env
