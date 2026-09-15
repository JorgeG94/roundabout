!! MPI integration test for halo exchange
!! Run with: mpirun -np 4 ./test_halo_mpi
!!
!! Tests a 2x2 process grid on a 10x10 global domain (5x5 per rank).
!! Each rank fills its physical cells with its rank number, then after
!! halo exchange, ghost cells should contain the neighbour's rank value.
!! Also tests MPI allreduce for timestep synchronisation.
#ifdef RDB_ENABLE_MPI
program test_halo_mpi
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_init, decomp_rank_from_coords
   use rdb_halo, only: halo_exchange_2d, halo_allreduce_min
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, &
                           comm_env_finalize, comm_env_rank, comm_env_size, &
                           comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   integer :: rank, nprocs
   integer :: nx_global, ny_global, px, py, nghost
   integer :: nx_local, ny_local, nx_total, ny_total
   integer :: i, j, k, n_fail
   type(decomp_t) :: decomp
   type(comm_t) :: comm
   real(wp), allocatable :: fld(:, :)
   real(wp) :: expected, local_dt, global_dt
   integer :: rank_west, rank_east, rank_south, rank_north

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()

   if (nprocs /= 4) then
      if (rank == 0) write (*, *) "ERROR: test requires exactly 4 MPI ranks"
      call comm_env_finalize()
      error stop 1
   end if

   n_fail = 0
   nx_global = 10
   ny_global = 10
   px = 2
   py = 2
   nghost = 2

   call decomp_init(decomp, nx_global, ny_global, px, py, rank)
   nx_local = decomp%nx_local
   ny_local = decomp%ny_local
   nx_total = nx_local + 2*nghost
   ny_total = ny_local + 2*nghost

   allocate (fld(nx_total, ny_total))
   fld = -1.0_wp

   ! Fill physical cells with this rank's ID
   do j = nghost + 1, nghost + ny_local
      do i = nghost + 1, nghost + nx_local
         fld(i, j) = real(rank, wp)
      end do
   end do

   ! Perform halo exchange
   call halo_exchange_2d(fld, decomp, nghost, nx_local, ny_local)

   ! --- Check west ghost cells ---
   if (.not. decomp%has_west) then
      rank_west = decomp_rank_from_coords(px, decomp%rx - 1, decomp%ry)
      expected = real(rank_west, wp)
      do j = nghost + 1, nghost + ny_local
         do k = 1, nghost
            if (abs(fld(k, j) - expected) > 1.0e-10_wp) then
               n_fail = n_fail + 1
               if (n_fail <= 3) write (*, *) "FAIL west ghost: rank=", rank, &
                  " i=", k, " j=", j, " expected=", expected, " got=", fld(k, j)
            end if
         end do
      end do
   end if

   ! --- Check east ghost cells ---
   if (.not. decomp%has_east) then
      rank_east = decomp_rank_from_coords(px, decomp%rx + 1, decomp%ry)
      expected = real(rank_east, wp)
      do j = nghost + 1, nghost + ny_local
         do k = 1, nghost
            if (abs(fld(nghost + nx_local + k, j) - expected) > 1.0e-10_wp) then
               n_fail = n_fail + 1
               if (n_fail <= 3) write (*, *) "FAIL east ghost: rank=", rank, &
                  " i=", nghost + nx_local + k, " j=", j, &
                  " expected=", expected, " got=", fld(nghost + nx_local + k, j)
            end if
         end do
      end do
   end if

   ! --- Check south ghost cells ---
   if (.not. decomp%has_south) then
      rank_south = decomp_rank_from_coords(px, decomp%rx, decomp%ry - 1)
      expected = real(rank_south, wp)
      do k = 1, nghost
         do i = nghost + 1, nghost + nx_local
            if (abs(fld(i, k) - expected) > 1.0e-10_wp) then
               n_fail = n_fail + 1
               if (n_fail <= 3) write (*, *) "FAIL south ghost: rank=", rank, &
                  " i=", i, " j=", k, " expected=", expected, " got=", fld(i, k)
            end if
         end do
      end do
   end if

   ! --- Check north ghost cells ---
   if (.not. decomp%has_north) then
      rank_north = decomp_rank_from_coords(px, decomp%rx, decomp%ry + 1)
      expected = real(rank_north, wp)
      do k = 1, nghost
         do i = nghost + 1, nghost + nx_local
            if (abs(fld(i, nghost + ny_local + k) - expected) > 1.0e-10_wp) then
               n_fail = n_fail + 1
               if (n_fail <= 3) write (*, *) "FAIL north ghost: rank=", rank, &
                  " i=", i, " j=", nghost + ny_local + k, &
                  " expected=", expected, " got=", fld(i, nghost + ny_local + k)
            end if
         end do
      end do
   end if

   ! --- Test allreduce min ---
   local_dt = real(rank + 1, wp)*0.1_wp   ! rank 0->0.1, rank 1->0.2, etc.
   call halo_allreduce_min(local_dt, global_dt)

   if (abs(global_dt - 0.1_wp) > 1.0e-10_wp) then
      n_fail = n_fail + 1
      if (rank == 0) write (*, *) "FAIL allreduce_min: expected 0.1, got", global_dt
   end if

   deallocate (fld)

   ! Gather failures to rank 0
   comm = comm_env_compute_comm()
   call comm%barrier()

   if (n_fail > 0) then
      write (*, *) "Rank", rank, ": ", n_fail, " test(s) FAILED"
   else
      write (*, *) "Rank", rank, ": all tests PASSED"
   end if

   ! Global reduce of failures
   block
      integer :: total_fail
      total_fail = n_fail
      call allreduce(comm, total_fail, MPI_SUM)
      call comm_env_finalize()
      if (total_fail > 0) error stop 1
   end block

end program test_halo_mpi
#endif
