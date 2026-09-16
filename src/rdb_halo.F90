!! Halo exchange for domain decomposition — the ONLY backend.
!! See `rdb_comm_env.F90` for why there is no longer a `_stub` twin.
#ifdef RDB_DOUBLE_PRECISION
#define HALO_ISEND_N comm_isend_real_dp_array_n
#define HALO_IRECV_N comm_irecv_real_dp_array_n
#else
#define HALO_ISEND_N comm_isend_real_sp_array_n
#define HALO_IRECV_N comm_irecv_real_sp_array_n
#endif
module rdb_halo
   !! Exchanges ghost-cell strips between neighbouring MPI ranks
   !!
   !! Uses pic_mpi_lib for non-blocking sends/receives.
   !! Plus-shaped stencil: 4 exchanges (N/S/E/W), no corner exchanges.
   !!
   !! Three exchange modes:
   !!   halo_exchange_2d        — host-staged, blocking (caller copies device<->host)
   !!   halo_exchange_2d_device — GPU-direct, blocking (CUDA-aware MPI)
   !!   halo_exchange_begin/end — async split for comm/compute overlap
   use, intrinsic :: iso_fortran_env, only: real64, int64
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_rank_from_coords
   use pic_mpi_lib, only: comm_t, request_t, MPI_Status, &
                          isend, irecv, waitall, allreduce, MPI_MIN, MPI_MAX, MPI_SUM, &
                          HALO_ISEND_N, HALO_IRECV_N
   use rdb_comm_env, only: comm_env_compute_comm
   use rdb_efp, only: efp_t, EFP_DIGITS, EFP_MAX_RANKS, &
                      efp_to_transport, efp_from_transport, efp_carry, &
                      efp_bin1_within_transport_bound
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: halo_exchange_2d
   public :: halo_exchange_2d_device
   public :: halo_allreduce_min
   public :: halo_allreduce_max
   public :: halo_async_t
   public :: halo_async_init
   public :: halo_async_destroy
   public :: halo_exchange_begin
   public :: halo_exchange_end
   public :: halo_exchange_3d
   public :: halo_exchange_3d_device
   public :: halo_allreduce_sum
   public :: halo_allreduce_sum_i8
   public :: halo_sync_buffers_ensure
   public :: halo_sync_buffers_cleanup
   public :: halo_sync_buffers_ensure_3d
   public :: halo_sync_buffers_cleanup_3d
   public :: halo_allreduce_efp_list

   integer, parameter :: NFIELDS = 4
      !! Number of fields exchanged (h, hu, hv, b)
   integer, parameter :: MAX_REQS = 8
      !! Max MPI requests: 4 directions x 2 (send + recv)

   type :: halo_async_t
      !! Persistent state for split begin/end halo exchange
      !!
      !! Pre-allocates combined send/recv buffers for all 4 fields
      !! in each direction, avoiding per-timestep allocate/deallocate.
      !! Buffers are stored at module level (ha_buf_*) so OpenACC
      !! can resolve present() lookups without derived-type traversal.
      type(request_t) :: reqs(MAX_REQS)
      integer :: nreq = 0
      integer :: strip_ew = 0
      integer :: strip_sn = 0
      integer :: nx_total = 0, ny_total = 0
      integer :: nghost = 0, nx_local = 0, ny_local = 0
      integer :: rank_west = -1, rank_east = -1
      integer :: rank_south = -1, rank_north = -1
      type(decomp_t) :: decomp
      logical :: initialised = .false.
   end type halo_async_t

   ! Module-level send/recv buffers for async halo exchange.
   ! Stored here (not in halo_async_t) because nvhpc OpenACC cannot
   ! resolve present() through derived-type component access.
   real(wp), allocatable :: ha_buf_send_west(:), ha_buf_recv_west(:)
   real(wp), allocatable :: ha_buf_send_east(:), ha_buf_recv_east(:)
   real(wp), allocatable :: ha_buf_send_south(:), ha_buf_recv_south(:)
   real(wp), allocatable :: ha_buf_send_north(:), ha_buf_recv_north(:)

   ! Module-level persistent buffers for the synchronous halo exchange
   ! (`halo_exchange_2d` and `halo_exchange_2d_device`).  Allocated once
   ! by `halo_sync_buffers_ensure` (called from `solver_enter_data`) and
   ! freed by `halo_sync_buffers_cleanup` (from `solver_exit_data`).
   ! Reusing these across the ~50 halo calls per ML step avoids the
   ! cudaMalloc/cudaFree storm that otherwise serialises the runtime
   ! allocator across ranks.  Single-field strips (the async buffers
   ! pack NFIELDS at once and aren't reusable here).
   real(wp), allocatable :: hs_buf_send_west(:), hs_buf_recv_west(:)
   real(wp), allocatable :: hs_buf_send_east(:), hs_buf_recv_east(:)
   real(wp), allocatable :: hs_buf_send_south(:), hs_buf_recv_south(:)
   real(wp), allocatable :: hs_buf_send_north(:), hs_buf_recv_north(:)
   integer :: hs_nghost = 0, hs_nx_total = 0, hs_ny_total = 0

   ! Module-level persistent buffers for the 3D batched halo exchange
   ! (`halo_exchange_3d_device`).  Sized for nghost*ny_total*nz_capacity
   ! (E/W) and nx_total*nghost*nz_capacity (S/N).  Capacity grows on
   ! demand so a single allocation services both ML (nz_ml) and NH
   ! (nz_ml+1 for w) without thrashing.  One MPI message per direction
   ! covering all `nz` layers, in place of the per-layer loop that used
   ! to go through halo_exchange_2d_device.
   real(wp), allocatable :: hs3_buf_send_west(:), hs3_buf_recv_west(:)
   real(wp), allocatable :: hs3_buf_send_east(:), hs3_buf_recv_east(:)
   real(wp), allocatable :: hs3_buf_send_south(:), hs3_buf_recv_south(:)
   real(wp), allocatable :: hs3_buf_send_north(:), hs3_buf_recv_north(:)
   integer :: hs3_nghost = 0, hs3_nx_total = 0, hs3_ny_total = 0
   integer :: hs3_nz_capacity = 0

contains

   subroutine halo_exchange_2d(fld, decomp, nghost, nx_local, ny_local)
      !! Exchange ghost-cell halos for a single 2D field
      !!
      !! The field has dimensions (nx_local + 2*nghost, ny_local + 2*nghost).
      !! Physical cells occupy indices (nghost+1 : nghost+nx_local, nghost+1 : nghost+ny_local).
      !! This routine fills the nghost-wide ghost strips on each side by
      !! sending/receiving from neighbouring ranks.
      real(wp), intent(inout) :: fld(:, :)
         !! 2D field with ghost cells
      type(decomp_t), intent(in) :: decomp
         !! Domain decomposition descriptor
      integer, intent(in) :: nghost
         !! Ghost cell width
      integer, intent(in) :: nx_local
         !! Local physical cells in x
      integer, intent(in) :: ny_local
         !! Local physical cells in y

      type(comm_t) :: comm
      integer :: nx_total, ny_total
      integer :: i, j, k, idx
      integer :: rank_west, rank_east, rank_south, rank_north
      real(wp), allocatable :: send_west(:), send_east(:)
      real(wp), allocatable :: recv_west(:), recv_east(:)
      real(wp), allocatable :: send_south(:), send_north(:)
      real(wp), allocatable :: recv_south(:), recv_north(:)
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      integer :: nreq
      integer :: strip_ew, strip_sn

      comm = comm_env_compute_comm()

      nx_total = nx_local + 2*nghost
      ny_total = ny_local + 2*nghost

      ! East-West strips: nghost columns x ny_total rows
      strip_ew = nghost*ny_total
      ! South-North strips: nx_total columns x nghost rows
      strip_sn = nx_total*nghost

      nreq = 0

      ! --- East/West exchange ---
      if (.not. decomp%has_west) then
         rank_west = decomp_rank_from_coords(decomp%px, decomp%rx - 1, decomp%ry)
         allocate (send_west(strip_ew), recv_west(strip_ew))

         ! Pack west send buffer: columns nghost+1 .. 2*nghost (first nghost physical columns)
         idx = 0
         do j = 1, ny_total
            do k = 1, nghost
               idx = idx + 1
               send_west(idx) = fld(nghost + k, j)
            end do
         end do

         nreq = nreq + 1
         call isend(comm, send_west, rank_west, 1, reqs(nreq))
         nreq = nreq + 1
         call irecv(comm, recv_west, rank_west, 2, reqs(nreq))
      end if

      if (.not. decomp%has_east) then
         rank_east = decomp_rank_from_coords(decomp%px, decomp%rx + 1, decomp%ry)
         allocate (send_east(strip_ew), recv_east(strip_ew))

         ! Pack east send buffer: columns nx_local+1 .. nx_local+nghost (last nghost physical cols)
         idx = 0
         do j = 1, ny_total
            do k = 1, nghost
               idx = idx + 1
               send_east(idx) = fld(nghost + nx_local - nghost + k, j)
            end do
         end do

         nreq = nreq + 1
         call isend(comm, send_east, rank_east, 2, reqs(nreq))
         nreq = nreq + 1
         call irecv(comm, recv_east, rank_east, 1, reqs(nreq))
      end if

      ! --- South/North exchange ---
      if (.not. decomp%has_south) then
         rank_south = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry - 1)
         allocate (send_south(strip_sn), recv_south(strip_sn))

         ! Pack south send buffer: rows nghost+1 .. 2*nghost (first nghost physical rows)
         idx = 0
         do k = 1, nghost
            do i = 1, nx_total
               idx = idx + 1
               send_south(idx) = fld(i, nghost + k)
            end do
         end do

         nreq = nreq + 1
         call isend(comm, send_south, rank_south, 3, reqs(nreq))
         nreq = nreq + 1
         call irecv(comm, recv_south, rank_south, 4, reqs(nreq))
      end if

      if (.not. decomp%has_north) then
         rank_north = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry + 1)
         allocate (send_north(strip_sn), recv_north(strip_sn))

         ! Pack north send buffer: rows ny_local+1 .. ny_local+nghost (last nghost physical rows)
         idx = 0
         do k = 1, nghost
            do i = 1, nx_total
               idx = idx + 1
               send_north(idx) = fld(i, nghost + ny_local - nghost + k)
            end do
         end do

         nreq = nreq + 1
         call isend(comm, send_north, rank_north, 4, reqs(nreq))
         nreq = nreq + 1
         call irecv(comm, recv_north, rank_north, 3, reqs(nreq))
      end if

      ! Wait for all
      if (nreq > 0) then
         call waitall(reqs(1:nreq), stats(1:nreq))
      end if

      ! --- Unpack received data into ghost cells ---
      if (.not. decomp%has_west) then
         idx = 0
         do j = 1, ny_total
            do k = 1, nghost
               idx = idx + 1
               fld(k, j) = recv_west(idx)
            end do
         end do
         deallocate (send_west, recv_west)
      end if

      if (.not. decomp%has_east) then
         idx = 0
         do j = 1, ny_total
            do k = 1, nghost
               idx = idx + 1
               fld(nghost + nx_local + k, j) = recv_east(idx)
            end do
         end do
         deallocate (send_east, recv_east)
      end if

      if (.not. decomp%has_south) then
         idx = 0
         do k = 1, nghost
            do i = 1, nx_total
               idx = idx + 1
               fld(i, k) = recv_south(idx)
            end do
         end do
         deallocate (send_south, recv_south)
      end if

      if (.not. decomp%has_north) then
         idx = 0
         do k = 1, nghost
            do i = 1, nx_total
               idx = idx + 1
               fld(i, nghost + ny_local + k) = recv_north(idx)
            end do
         end do
         deallocate (send_north, recv_north)
      end if

   end subroutine halo_exchange_2d

   subroutine halo_exchange_2d_device(fld, decomp, nghost, nx_local, ny_local)
      !! GPU-direct halo exchange via CUDA-aware MPI
      !!
      !! Pack/unpack buffers live on the device. MPI operates on device
      !! pointers via !$acc host_data use_device, eliminating the
      !! full-array GPU<->host copies required by the host-staged path.
      !!
      !! Uses persistent module-level send/recv buffers
      !! (`hs_buf_*`) so back-to-back halo calls do not hit the
      !! cudaMalloc/cudaFree allocator on every call.  Lazy-allocated by
      !! `halo_sync_buffers_ensure`; the solver normally calls that from
      !! `solver_enter_data`, but a guard here keeps the routine
      !! self-contained for callers that haven't been migrated yet.
      real(wp), intent(inout) :: fld(:, :)
         !! 2D field with ghost cells (present on device)
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nghost
      integer, intent(in) :: nx_local
      integer, intent(in) :: ny_local

      type(comm_t) :: comm
      integer :: nx_total, ny_total
      integer :: i, j, k
      integer :: rank_west, rank_east, rank_south, rank_north
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      integer :: nreq
      integer :: strip_ew, strip_sn

      comm = comm_env_compute_comm()
      nx_total = nx_local + 2*nghost
      ny_total = ny_local + 2*nghost
      strip_ew = nghost*ny_total
      strip_sn = nx_total*nghost

      call halo_sync_buffers_ensure(nghost, nx_total, ny_total)

      ! --- Pack on device ---
      if (.not. decomp%has_west) then
         rank_west = decomp_rank_from_coords(decomp%px, decomp%rx - 1, decomp%ry)
         !$acc parallel loop collapse(2) present(hs_buf_send_west, fld)
         do j = 1, ny_total
            do k = 1, nghost
               hs_buf_send_west((j - 1)*nghost + k) = fld(nghost + k, j)
            end do
         end do
      end if

      if (.not. decomp%has_east) then
         rank_east = decomp_rank_from_coords(decomp%px, decomp%rx + 1, decomp%ry)
         !$acc parallel loop collapse(2) present(hs_buf_send_east, fld)
         do j = 1, ny_total
            do k = 1, nghost
               hs_buf_send_east((j - 1)*nghost + k) = fld(nghost + nx_local - nghost + k, j)
            end do
         end do
      end if

      if (.not. decomp%has_south) then
         rank_south = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry - 1)
         !$acc parallel loop collapse(2) present(hs_buf_send_south, fld)
         do k = 1, nghost
            do i = 1, nx_total
               hs_buf_send_south((k - 1)*nx_total + i) = fld(i, nghost + k)
            end do
         end do
      end if

      if (.not. decomp%has_north) then
         rank_north = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry + 1)
         !$acc parallel loop collapse(2) present(hs_buf_send_north, fld)
         do k = 1, nghost
            do i = 1, nx_total
               hs_buf_send_north((k - 1)*nx_total + i) = fld(i, nghost + ny_local - nghost + k)
            end do
         end do
      end if

      ! --- MPI with device pointers ---
      nreq = 0

      if (.not. decomp%has_west) then
         !$acc host_data use_device(hs_buf_send_west, hs_buf_recv_west)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs_buf_send_west, strip_ew, rank_west, 1, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs_buf_recv_west, strip_ew, rank_west, 2, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_east) then
         !$acc host_data use_device(hs_buf_send_east, hs_buf_recv_east)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs_buf_send_east, strip_ew, rank_east, 2, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs_buf_recv_east, strip_ew, rank_east, 1, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_south) then
         !$acc host_data use_device(hs_buf_send_south, hs_buf_recv_south)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs_buf_send_south, strip_sn, rank_south, 3, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs_buf_recv_south, strip_sn, rank_south, 4, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_north) then
         !$acc host_data use_device(hs_buf_send_north, hs_buf_recv_north)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs_buf_send_north, strip_sn, rank_north, 4, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs_buf_recv_north, strip_sn, rank_north, 3, reqs(nreq))
         !$acc end host_data
      end if

      if (nreq > 0) call waitall(reqs(1:nreq), stats(1:nreq))

      ! --- Unpack on device ---
      if (.not. decomp%has_west) then
         !$acc parallel loop collapse(2) present(hs_buf_recv_west, fld)
         do j = 1, ny_total
            do k = 1, nghost
               fld(k, j) = hs_buf_recv_west((j - 1)*nghost + k)
            end do
         end do
      end if

      if (.not. decomp%has_east) then
         !$acc parallel loop collapse(2) present(hs_buf_recv_east, fld)
         do j = 1, ny_total
            do k = 1, nghost
               fld(nghost + nx_local + k, j) = hs_buf_recv_east((j - 1)*nghost + k)
            end do
         end do
      end if

      if (.not. decomp%has_south) then
         !$acc parallel loop collapse(2) present(hs_buf_recv_south, fld)
         do k = 1, nghost
            do i = 1, nx_total
               fld(i, k) = hs_buf_recv_south((k - 1)*nx_total + i)
            end do
         end do
      end if

      if (.not. decomp%has_north) then
         !$acc parallel loop collapse(2) present(hs_buf_recv_north, fld)
         do k = 1, nghost
            do i = 1, nx_total
               fld(i, nghost + ny_local + k) = hs_buf_recv_north((k - 1)*nx_total + i)
            end do
         end do
      end if

   end subroutine halo_exchange_2d_device

   subroutine halo_sync_buffers_ensure(nghost, nx_total, ny_total)
      !! Lazy-allocate the persistent send/recv buffers used by
      !! `halo_exchange_2d_device` (and `halo_exchange_2d` via the same
      !! pool).  Sized to the per-rank subdomain; resizes on grid change.
      !!
      !! Allocates conservatively in all four directions even if this
      !! rank is at a domain edge -- the pack/unpack guards in
      !! `halo_exchange_2d_device` then skip directions that aren't
      !! used, but the buffers exist so back-to-back calls don't hit
      !! cudaMalloc.  At ~strip_ew + strip_sn doubles per direction the
      !! total residency is negligible compared with the field arrays.
      integer, intent(in) :: nghost, nx_total, ny_total

      integer :: strip_ew, strip_sn
      logical :: needs_realloc

      strip_ew = nghost*ny_total
      strip_sn = nx_total*nghost

      needs_realloc = (.not. allocated(hs_buf_send_west)) .or. &
                      hs_nghost /= nghost .or. &
                      hs_nx_total /= nx_total .or. &
                      hs_ny_total /= ny_total

      if (.not. needs_realloc) return

      call halo_sync_buffers_cleanup()

      allocate (hs_buf_send_west(strip_ew), hs_buf_recv_west(strip_ew))
      allocate (hs_buf_send_east(strip_ew), hs_buf_recv_east(strip_ew))
      allocate (hs_buf_send_south(strip_sn), hs_buf_recv_south(strip_sn))
      allocate (hs_buf_send_north(strip_sn), hs_buf_recv_north(strip_sn))
      !$acc enter data create(hs_buf_send_west, hs_buf_recv_west, &
      !$acc&                  hs_buf_send_east, hs_buf_recv_east, &
      !$acc&                  hs_buf_send_south, hs_buf_recv_south, &
      !$acc&                  hs_buf_send_north, hs_buf_recv_north)

      hs_nghost = nghost
      hs_nx_total = nx_total
      hs_ny_total = ny_total

   end subroutine halo_sync_buffers_ensure

   subroutine halo_sync_buffers_cleanup()
      !! Release the persistent halo buffers.  Idempotent.
      if (.not. allocated(hs_buf_send_west)) return

      !$acc exit data delete(hs_buf_send_west, hs_buf_recv_west, &
      !$acc&                 hs_buf_send_east, hs_buf_recv_east, &
      !$acc&                 hs_buf_send_south, hs_buf_recv_south, &
      !$acc&                 hs_buf_send_north, hs_buf_recv_north)
      deallocate (hs_buf_send_west, hs_buf_recv_west)
      deallocate (hs_buf_send_east, hs_buf_recv_east)
      deallocate (hs_buf_send_south, hs_buf_recv_south)
      deallocate (hs_buf_send_north, hs_buf_recv_north)

      hs_nghost = 0
      hs_nx_total = 0
      hs_ny_total = 0

   end subroutine halo_sync_buffers_cleanup

   ! ================================================================
   ! Async split halo exchange for comm/compute overlap
   ! ================================================================

   subroutine halo_async_init(ha, decomp, nghost, nx_local, ny_local)
      !! Pre-allocate halo buffers for all 4 fields
      type(halo_async_t), intent(out) :: ha
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nghost, nx_local, ny_local

      integer :: nx_total, ny_total

      ha%decomp = decomp
      ha%nghost = nghost
      ha%nx_local = nx_local
      ha%ny_local = ny_local
      nx_total = nx_local + 2*nghost
      ny_total = ny_local + 2*nghost
      ha%nx_total = nx_total
      ha%ny_total = ny_total

      ! Per-field strip sizes, multiplied by NFIELDS for combined buffers
      ha%strip_ew = NFIELDS*nghost*ny_total
      ha%strip_sn = NFIELDS*nx_total*nghost

      if (.not. decomp%has_west) then
         ha%rank_west = decomp_rank_from_coords(decomp%px, decomp%rx - 1, decomp%ry)
         allocate (ha_buf_send_west(ha%strip_ew), ha_buf_recv_west(ha%strip_ew))
         !$acc enter data create( ha_buf_send_west, ha_buf_recv_west)
      end if
      if (.not. decomp%has_east) then
         ha%rank_east = decomp_rank_from_coords(decomp%px, decomp%rx + 1, decomp%ry)
         allocate (ha_buf_send_east(ha%strip_ew), ha_buf_recv_east(ha%strip_ew))
         !$acc enter data create( ha_buf_send_east, ha_buf_recv_east)
      end if
      if (.not. decomp%has_south) then
         ha%rank_south = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry - 1)
         allocate (ha_buf_send_south(ha%strip_sn), ha_buf_recv_south(ha%strip_sn))
         !$acc enter data create( ha_buf_send_south, ha_buf_recv_south)
      end if
      if (.not. decomp%has_north) then
         ha%rank_north = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry + 1)
         allocate (ha_buf_send_north(ha%strip_sn), ha_buf_recv_north(ha%strip_sn))
         !$acc enter data create( ha_buf_send_north, ha_buf_recv_north)
      end if

      ha%initialised = .true.

   end subroutine halo_async_init

   subroutine halo_async_destroy(ha)
      !! Free pre-allocated halo buffers
      type(halo_async_t), intent(inout) :: ha

      if (.not. ha%initialised) return

      if (allocated(ha_buf_send_west)) then
         !$acc exit data delete( ha_buf_send_west, ha_buf_recv_west)
         deallocate (ha_buf_send_west, ha_buf_recv_west)
      end if
      if (allocated(ha_buf_send_east)) then
         !$acc exit data delete( ha_buf_send_east, ha_buf_recv_east)
         deallocate (ha_buf_send_east, ha_buf_recv_east)
      end if
      if (allocated(ha_buf_send_south)) then
         !$acc exit data delete( ha_buf_send_south, ha_buf_recv_south)
         deallocate (ha_buf_send_south, ha_buf_recv_south)
      end if
      if (allocated(ha_buf_send_north)) then
         !$acc exit data delete( ha_buf_send_north, ha_buf_recv_north)
         deallocate (ha_buf_send_north, ha_buf_recv_north)
      end if

      ha%initialised = .false.

   end subroutine halo_async_destroy

   subroutine halo_exchange_begin(ha, h, hu, hv, b_fld)
      !! Pack and post non-blocking MPI sends/recvs for all 4 fields
      !!
      !! After this returns, the interior cells (not touching ghost cells)
      !! can be computed while the exchange is in flight.
      !!
      !! Buffers are module-level arrays (ha_buf_*) so that OpenACC can
      !! resolve present() lookups without derived-type traversal.
      type(halo_async_t), intent(inout) :: ha
      real(wp), intent(in) :: h(:, :), hu(:, :), hv(:, :), b_fld(:, :)

      type(comm_t) :: comm
      integer :: i, j, k, base
      integer :: ng, nxl, nyl, nxt, nyt

      comm = comm_env_compute_comm()
      ng = ha%nghost
      nxl = ha%nx_local
      nyl = ha%ny_local
      nxt = ha%nx_total
      nyt = ha%ny_total
      ha%nreq = 0

      ! --- Pack all 4 fields into combined buffers on device ---
      if (.not. ha%decomp%has_west) then
         !$acc parallel loop collapse(2)
         do j = 1, nyt
            do k = 1, ng
               base = (j - 1)*ng + k
               ha_buf_send_west(base) = h(ng + k, j)
               ha_buf_send_west(base + ng*nyt) = hu(ng + k, j)
               ha_buf_send_west(base + 2*ng*nyt) = hv(ng + k, j)
               ha_buf_send_west(base + 3*ng*nyt) = b_fld(ng + k, j)
            end do
         end do
      end if

      if (.not. ha%decomp%has_east) then
         !$acc parallel loop collapse(2)
         do j = 1, nyt
            do k = 1, ng
               base = (j - 1)*ng + k
               ha_buf_send_east(base) = h(ng + nxl - ng + k, j)
               ha_buf_send_east(base + ng*nyt) = hu(ng + nxl - ng + k, j)
               ha_buf_send_east(base + 2*ng*nyt) = hv(ng + nxl - ng + k, j)
               ha_buf_send_east(base + 3*ng*nyt) = b_fld(ng + nxl - ng + k, j)
            end do
         end do
      end if

      if (.not. ha%decomp%has_south) then
         !$acc parallel loop collapse(2)
         do k = 1, ng
            do i = 1, nxt
               base = (k - 1)*nxt + i
               ha_buf_send_south(base) = h(i, ng + k)
               ha_buf_send_south(base + nxt*ng) = hu(i, ng + k)
               ha_buf_send_south(base + 2*nxt*ng) = hv(i, ng + k)
               ha_buf_send_south(base + 3*nxt*ng) = b_fld(i, ng + k)
            end do
         end do
      end if

      if (.not. ha%decomp%has_north) then
         !$acc parallel loop collapse(2)
         do k = 1, ng
            do i = 1, nxt
               base = (k - 1)*nxt + i
               ha_buf_send_north(base) = h(i, ng + nyl - ng + k)
               ha_buf_send_north(base + nxt*ng) = hu(i, ng + nyl - ng + k)
               ha_buf_send_north(base + 2*nxt*ng) = hv(i, ng + nyl - ng + k)
               ha_buf_send_north(base + 3*nxt*ng) = b_fld(i, ng + nyl - ng + k)
            end do
         end do
      end if

      ! --- Post MPI Isend/Irecv with device pointers ---
      if (.not. ha%decomp%has_west) then
         !$acc host_data use_device(ha_buf_send_west, ha_buf_recv_west)
         ha%nreq = ha%nreq + 1
         call HALO_ISEND_N(comm, ha_buf_send_west, ha%strip_ew, ha%rank_west, 1, ha%reqs(ha%nreq))
         ha%nreq = ha%nreq + 1
         call HALO_IRECV_N(comm, ha_buf_recv_west, ha%strip_ew, ha%rank_west, 2, ha%reqs(ha%nreq))
         !$acc end host_data
      end if

      if (.not. ha%decomp%has_east) then
         !$acc host_data use_device(ha_buf_send_east, ha_buf_recv_east)
         ha%nreq = ha%nreq + 1
         call HALO_ISEND_N(comm, ha_buf_send_east, ha%strip_ew, ha%rank_east, 2, ha%reqs(ha%nreq))
         ha%nreq = ha%nreq + 1
         call HALO_IRECV_N(comm, ha_buf_recv_east, ha%strip_ew, ha%rank_east, 1, ha%reqs(ha%nreq))
         !$acc end host_data
      end if

      if (.not. ha%decomp%has_south) then
         !$acc host_data use_device(ha_buf_send_south, ha_buf_recv_south)
         ha%nreq = ha%nreq + 1
         call HALO_ISEND_N(comm, ha_buf_send_south, ha%strip_sn, ha%rank_south, 3, ha%reqs(ha%nreq))
         ha%nreq = ha%nreq + 1
         call HALO_IRECV_N(comm, ha_buf_recv_south, ha%strip_sn, ha%rank_south, 4, ha%reqs(ha%nreq))
         !$acc end host_data
      end if

      if (.not. ha%decomp%has_north) then
         !$acc host_data use_device(ha_buf_send_north, ha_buf_recv_north)
         ha%nreq = ha%nreq + 1
         call HALO_ISEND_N(comm, ha_buf_send_north, ha%strip_sn, ha%rank_north, 4, ha%reqs(ha%nreq))
         ha%nreq = ha%nreq + 1
         call HALO_IRECV_N(comm, ha_buf_recv_north, ha%strip_sn, ha%rank_north, 3, ha%reqs(ha%nreq))
         !$acc end host_data
      end if

   end subroutine halo_exchange_begin

   subroutine halo_exchange_end(ha, h, hu, hv, b_fld)
      !! Wait for MPI to complete and unpack received ghost cells
      !!
      !! Buffers are module-level arrays (ha_buf_*) — see
      !! halo_exchange_begin for explanation.
      type(halo_async_t), intent(inout) :: ha
      real(wp), intent(inout) :: h(:, :), hu(:, :), hv(:, :), b_fld(:, :)

      type(MPI_Status) :: stats(MAX_REQS)
      integer :: i, j, k, base
      integer :: ng, nxl, nyl, nxt, nyt

      ng = ha%nghost
      nxl = ha%nx_local
      nyl = ha%ny_local
      nxt = ha%nx_total
      nyt = ha%ny_total

      if (ha%nreq > 0) then
         call waitall(ha%reqs(1:ha%nreq), stats(1:ha%nreq))
      end if

      ! --- Unpack all 4 fields from combined buffers on device ---
      if (.not. ha%decomp%has_west) then
         !$acc parallel loop collapse(2)
         do j = 1, nyt
            do k = 1, ng
               base = (j - 1)*ng + k
               h(k, j) = ha_buf_recv_west(base)
               hu(k, j) = ha_buf_recv_west(base + ng*nyt)
               hv(k, j) = ha_buf_recv_west(base + 2*ng*nyt)
               b_fld(k, j) = ha_buf_recv_west(base + 3*ng*nyt)
            end do
         end do
      end if

      if (.not. ha%decomp%has_east) then
         !$acc parallel loop collapse(2)
         do j = 1, nyt
            do k = 1, ng
               base = (j - 1)*ng + k
               h(ng + nxl + k, j) = ha_buf_recv_east(base)
               hu(ng + nxl + k, j) = ha_buf_recv_east(base + ng*nyt)
               hv(ng + nxl + k, j) = ha_buf_recv_east(base + 2*ng*nyt)
               b_fld(ng + nxl + k, j) = ha_buf_recv_east(base + 3*ng*nyt)
            end do
         end do
      end if

      if (.not. ha%decomp%has_south) then
         !$acc parallel loop collapse(2)
         do k = 1, ng
            do i = 1, nxt
               base = (k - 1)*nxt + i
               h(i, k) = ha_buf_recv_south(base)
               hu(i, k) = ha_buf_recv_south(base + nxt*ng)
               hv(i, k) = ha_buf_recv_south(base + 2*nxt*ng)
               b_fld(i, k) = ha_buf_recv_south(base + 3*nxt*ng)
            end do
         end do
      end if

      if (.not. ha%decomp%has_north) then
         !$acc parallel loop collapse(2)
         do k = 1, ng
            do i = 1, nxt
               base = (k - 1)*nxt + i
               h(i, ng + nyl + k) = ha_buf_recv_north(base)
               hu(i, ng + nyl + k) = ha_buf_recv_north(base + nxt*ng)
               hv(i, ng + nyl + k) = ha_buf_recv_north(base + 2*nxt*ng)
               b_fld(i, ng + nyl + k) = ha_buf_recv_north(base + 3*nxt*ng)
            end do
         end do
      end if

   end subroutine halo_exchange_end

   subroutine halo_allreduce_min(local_val, global_val)
      !! MPI_Allreduce with MPI_MIN for global timestep
      real(wp), intent(in) :: local_val
      real(wp), intent(out) :: global_val

      type(comm_t) :: comm

      comm = comm_env_compute_comm()
      ! Single rank: the reduction is the identity, so return this rank's
      ! own contribution without entering a collective.  Not just an
      ! optimisation -- pic-mpi's serial backend (PIC_ENABLE_MPI=OFF)
      ! deliberately `error stop`s in `allreduce`, pushing the size()==1
      ! case onto the caller.  This IS that case.
      if (comm%size() == 1) then
         global_val = local_val
         return
      end if
      call allreduce(comm, local_val, global_val, op=MPI_MIN)

   end subroutine halo_allreduce_min

   subroutine halo_allreduce_max(local_val, global_val)
      !! MPI_Allreduce with MPI_MAX — max-type reductions are exact in FP,
      !! so a global max stays layout-reproducible (ocean-MPI plan D5).
      !! Used for auto_n_inner's global gravity-wave CFL (shared n_inner).
      real(wp), intent(in) :: local_val
      real(wp), intent(out) :: global_val

      type(comm_t) :: comm

      comm = comm_env_compute_comm()
      ! Single rank: the reduction is the identity, so return this rank's
      ! own contribution without entering a collective.  Not just an
      ! optimisation -- pic-mpi's serial backend (PIC_ENABLE_MPI=OFF)
      ! deliberately `error stop`s in `allreduce`, pushing the size()==1
      ! case onto the caller.  This IS that case.
      if (comm%size() == 1) then
         global_val = local_val
         return
      end if
      call allreduce(comm, local_val, global_val, op=MPI_MAX)

   end subroutine halo_allreduce_max

   subroutine halo_exchange_3d(fld, decomp, nghost, nx_local, ny_local, nz)
      !! Exchange ghost-cell halos for a 3D field (all nz layers packed per direction)
      real(wp), intent(inout) :: fld(:, :, :)
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nghost, nx_local, ny_local, nz

      integer :: kk

      do kk = 1, nz
         call halo_exchange_2d(fld(:, :, kk), decomp, nghost, nx_local, ny_local)
      end do

   end subroutine halo_exchange_3d

   subroutine halo_exchange_3d_device(fld, decomp, nghost, nx_local, ny_local, nz)
      !! GPU-direct halo exchange for a 3D field, batched across layers.
      !!
      !! Packs all `nz` layers into one persistent device buffer per
      !! direction, fires one MPI Isend/Irecv per direction, and unpacks
      !! all layers in one kernel.  This replaces the previous
      !! `do kk = 1, nz; call halo_exchange_2d_device(fld(:,:,kk), ...)`
      !! loop, which (a) issued nz×4 MPI messages per call and (b)
      !! tripped NVHPC's per-call array-section descriptor push for each
      !! 2D slice -- both visible as gaps in the timeline.
      !!
      !! `fld` is declared explicit-shape so that NVHPC passes a raw
      !! pointer + bounds rather than a descriptor that has to be
      !! re-attached to the device on each call.
      integer, intent(in) :: nghost, nx_local, ny_local, nz
      real(wp), intent(inout) :: fld(nx_local + 2*nghost, ny_local + 2*nghost, nz)
      type(decomp_t), intent(in) :: decomp

      type(comm_t) :: comm
      integer :: nx_total, ny_total
      integer :: i, j, k, L
      integer :: rank_west, rank_east, rank_south, rank_north
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      integer :: nreq
      integer :: strip_ew, strip_sn

      comm = comm_env_compute_comm()
      nx_total = nx_local + 2*nghost
      ny_total = ny_local + 2*nghost
      strip_ew = nghost*ny_total*nz
      strip_sn = nx_total*nghost*nz

      call halo_sync_buffers_ensure_3d(nghost, nx_total, ny_total, nz)

      ! --- Pack on device (all layers in a single kernel per direction) ---
      if (.not. decomp%has_west) then
         rank_west = decomp_rank_from_coords(decomp%px, decomp%rx - 1, decomp%ry)
         !$acc parallel loop collapse(3) present(hs3_buf_send_west, fld)
         do L = 1, nz
            do j = 1, ny_total
               do k = 1, nghost
                  hs3_buf_send_west(((L - 1)*ny_total + (j - 1))*nghost + k) = fld(nghost + k, j, L)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_east) then
         rank_east = decomp_rank_from_coords(decomp%px, decomp%rx + 1, decomp%ry)
         !$acc parallel loop collapse(3) present(hs3_buf_send_east, fld)
         do L = 1, nz
            do j = 1, ny_total
               do k = 1, nghost
                  hs3_buf_send_east(((L - 1)*ny_total + (j - 1))*nghost + k) = fld(nghost + nx_local - nghost + k, j, L)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_south) then
         rank_south = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry - 1)
         !$acc parallel loop collapse(3) present(hs3_buf_send_south, fld)
         do L = 1, nz
            do k = 1, nghost
               do i = 1, nx_total
                  hs3_buf_send_south(((L - 1)*nghost + (k - 1))*nx_total + i) = fld(i, nghost + k, L)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_north) then
         rank_north = decomp_rank_from_coords(decomp%px, decomp%rx, decomp%ry + 1)
         !$acc parallel loop collapse(3) present(hs3_buf_send_north, fld)
         do L = 1, nz
            do k = 1, nghost
               do i = 1, nx_total
                  hs3_buf_send_north(((L - 1)*nghost + (k - 1))*nx_total + i) = fld(i, nghost + ny_local - nghost + k, L)
               end do
            end do
         end do
      end if

      ! --- One MPI exchange per direction, all layers in one message ---
      nreq = 0

      if (.not. decomp%has_west) then
         !$acc host_data use_device(hs3_buf_send_west, hs3_buf_recv_west)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs3_buf_send_west, strip_ew, rank_west, 1, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs3_buf_recv_west, strip_ew, rank_west, 2, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_east) then
         !$acc host_data use_device(hs3_buf_send_east, hs3_buf_recv_east)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs3_buf_send_east, strip_ew, rank_east, 2, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs3_buf_recv_east, strip_ew, rank_east, 1, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_south) then
         !$acc host_data use_device(hs3_buf_send_south, hs3_buf_recv_south)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs3_buf_send_south, strip_sn, rank_south, 3, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs3_buf_recv_south, strip_sn, rank_south, 4, reqs(nreq))
         !$acc end host_data
      end if

      if (.not. decomp%has_north) then
         !$acc host_data use_device(hs3_buf_send_north, hs3_buf_recv_north)
         nreq = nreq + 1
         call HALO_ISEND_N(comm, hs3_buf_send_north, strip_sn, rank_north, 4, reqs(nreq))
         nreq = nreq + 1
         call HALO_IRECV_N(comm, hs3_buf_recv_north, strip_sn, rank_north, 3, reqs(nreq))
         !$acc end host_data
      end if

      if (nreq > 0) call waitall(reqs(1:nreq), stats(1:nreq))

      ! --- Unpack on device (all layers in a single kernel per direction) ---
      if (.not. decomp%has_west) then
         !$acc parallel loop collapse(3) present(hs3_buf_recv_west, fld)
         do L = 1, nz
            do j = 1, ny_total
               do k = 1, nghost
                  fld(k, j, L) = hs3_buf_recv_west(((L - 1)*ny_total + (j - 1))*nghost + k)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_east) then
         !$acc parallel loop collapse(3) present(hs3_buf_recv_east, fld)
         do L = 1, nz
            do j = 1, ny_total
               do k = 1, nghost
                  fld(nghost + nx_local + k, j, L) = hs3_buf_recv_east(((L - 1)*ny_total + (j - 1))*nghost + k)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_south) then
         !$acc parallel loop collapse(3) present(hs3_buf_recv_south, fld)
         do L = 1, nz
            do k = 1, nghost
               do i = 1, nx_total
                  fld(i, k, L) = hs3_buf_recv_south(((L - 1)*nghost + (k - 1))*nx_total + i)
               end do
            end do
         end do
      end if

      if (.not. decomp%has_north) then
         !$acc parallel loop collapse(3) present(hs3_buf_recv_north, fld)
         do L = 1, nz
            do k = 1, nghost
               do i = 1, nx_total
                  fld(i, nghost + ny_local + k, L) = hs3_buf_recv_north(((L - 1)*nghost + (k - 1))*nx_total + i)
               end do
            end do
         end do
      end if

   end subroutine halo_exchange_3d_device

   subroutine halo_sync_buffers_ensure_3d(nghost, nx_total, ny_total, nz)
      !! Lazy-allocate the persistent 3D send/recv buffers used by
      !! `halo_exchange_3d_device`.  Sized to the per-rank subdomain
      !! and the largest `nz` seen so far -- subsequent calls with a
      !! smaller `nz` reuse the existing buffer (just write fewer
      !! elements), grow-only on `nz` so an ML run that later does an
      !! NH `w` halo (nz_ml+1) doesn't free/reallocate.  Resizes in
      !! full if (nx_total, ny_total, nghost) change.
      integer, intent(in) :: nghost, nx_total, ny_total, nz

      integer :: cap_nz, strip_ew, strip_sn
      logical :: needs_full_realloc, needs_grow

      needs_full_realloc = (.not. allocated(hs3_buf_send_west)) .or. &
                           hs3_nghost /= nghost .or. &
                           hs3_nx_total /= nx_total .or. &
                           hs3_ny_total /= ny_total
      needs_grow = (.not. needs_full_realloc) .and. nz > hs3_nz_capacity

      if (.not. needs_full_realloc .and. .not. needs_grow) return

      call halo_sync_buffers_cleanup_3d()

      cap_nz = nz
      strip_ew = nghost*ny_total*cap_nz
      strip_sn = nx_total*nghost*cap_nz

      allocate (hs3_buf_send_west(strip_ew), hs3_buf_recv_west(strip_ew))
      allocate (hs3_buf_send_east(strip_ew), hs3_buf_recv_east(strip_ew))
      allocate (hs3_buf_send_south(strip_sn), hs3_buf_recv_south(strip_sn))
      allocate (hs3_buf_send_north(strip_sn), hs3_buf_recv_north(strip_sn))
      !$acc enter data create(hs3_buf_send_west, hs3_buf_recv_west, &
      !$acc&                  hs3_buf_send_east, hs3_buf_recv_east, &
      !$acc&                  hs3_buf_send_south, hs3_buf_recv_south, &
      !$acc&                  hs3_buf_send_north, hs3_buf_recv_north)

      hs3_nghost = nghost
      hs3_nx_total = nx_total
      hs3_ny_total = ny_total
      hs3_nz_capacity = cap_nz

   end subroutine halo_sync_buffers_ensure_3d

   subroutine halo_sync_buffers_cleanup_3d()
      !! Release the persistent 3D halo buffers.  Idempotent.
      if (.not. allocated(hs3_buf_send_west)) return

      !$acc exit data delete(hs3_buf_send_west, hs3_buf_recv_west, &
      !$acc&                 hs3_buf_send_east, hs3_buf_recv_east, &
      !$acc&                 hs3_buf_send_south, hs3_buf_recv_south, &
      !$acc&                 hs3_buf_send_north, hs3_buf_recv_north)
      deallocate (hs3_buf_send_west, hs3_buf_recv_west)
      deallocate (hs3_buf_send_east, hs3_buf_recv_east)
      deallocate (hs3_buf_send_south, hs3_buf_recv_south)
      deallocate (hs3_buf_send_north, hs3_buf_recv_north)

      hs3_nghost = 0
      hs3_nx_total = 0
      hs3_ny_total = 0
      hs3_nz_capacity = 0

   end subroutine halo_sync_buffers_cleanup_3d

   subroutine halo_allreduce_sum(local_val, global_val)
      !! MPI_Allreduce with MPI_SUM for CG dot products
      real(wp), intent(in) :: local_val
      real(wp), intent(out) :: global_val

      type(comm_t) :: comm

      comm = comm_env_compute_comm()
      ! Single rank: the reduction is the identity, so return this rank's
      ! own contribution without entering a collective.  Not just an
      ! optimisation -- pic-mpi's serial backend (PIC_ENABLE_MPI=OFF)
      ! deliberately `error stop`s in `allreduce`, pushing the size()==1
      ! case onto the caller.  This IS that case.
      if (comm%size() == 1) then
         global_val = local_val
         return
      end if
      call allreduce(comm, local_val, global_val, op=MPI_SUM)

   end subroutine halo_allreduce_sum

   subroutine halo_allreduce_sum_i8(local_val, global_val)
      !! Cross-rank int64 sum for the decomposition-invariant chksum
      !! bitcount (`rdb_ocean_chksum`).  `pic_mpi_lib` exposes no
      !! `integer(int64)` allreduce overload (MPI_INTEGER8 reaches only
      !! send/recv), so the value rides the EXACT `real64` allreduce: a
      !! per-field POPCNT sum is bounded by (#elements x 64), which stays
      !! FAR below the double-mantissa bound 2**53 for any realistic grid
      !! (2**53/64 ~ 1.4e14 cells), so both the local->double cast and
      !! every partial MPI_SUM are exact — the invariant survives.
      integer(int64), intent(in) :: local_val
      integer(int64), intent(out) :: global_val

      type(comm_t) :: comm
      real(real64) :: acc_local, acc_global

      comm = comm_env_compute_comm()
      ! Single rank: the reduction is the identity, so return this rank's
      ! own contribution without entering a collective.  Not just an
      ! optimisation -- pic-mpi's serial backend (PIC_ENABLE_MPI=OFF)
      ! deliberately `error stop`s in `allreduce`, pushing the size()==1
      ! case onto the caller.  This IS that case.
      if (comm%size() == 1) then
         global_val = local_val
         return
      end if

      acc_local = real(local_val, real64)
      call allreduce(comm, acc_local, acc_global, op=MPI_SUM)
      global_val = int(acc_global, int64)

   end subroutine halo_allreduce_sum_i8

   subroutine halo_allreduce_efp_list(local_list, global_list, nval)
      !! Order-invariant EXACT cross-rank combine of `nval` EFP values in
      !! ONE collective (PR-32).  Replaces N separate scalar
      !! `halo_allreduce_sum` calls with one packed `allreduce`.
      !!
      !! Transport mechanism (see `rdb_efp`'s module docstring for the
      !! full derivation): `pic_mpi_lib` has NO `integer(int64)` allreduce
      !! overload (`MPI_INTEGER8` reaches only send/recv, never
      !! allreduce), so the six int64 bins per value are packed as
      !! EXACTLY-representable `real64` doubles and combined with a plain
      !! `MPI_SUM`.  This is exact, not approximate, PROVIDED every
      !! partial sum MPI could form stays `<= 2**53` (the double mantissa
      !! bound) — enforced below by two fail-loud guards, never a silent
      !! fallback:
      !!
      !!   1. `num_ranks <= EFP_MAX_RANKS` (2**17 = 131072) — the bound
      !!      that makes bins 2..6 (each `< 2**P` after a local carry)
      !!      summable across ranks without exceeding 2**53.
      !!   2. `efp_bin1_within_transport_bound` on every LOCAL value —
      !!      bin 1 is NOT bounded by `efp_carry` (see `rdb_efp`), so it
      !!      needs its own runtime check: `|e(1)| <= 2**53 / num_ranks`.
      !!
      !! COLLECTIVE: every compute rank must call this with the SAME
      !! `nval`.  Both guards are evaluated on rank-uniform data (the
      !! rank count, and each rank's own local value) so a violation
      !! aborts identically on every rank — no rank-dependent branch that
      !! could hang the collective (CLAUDE.md's collective-panic idiom,
      !! `rdb_console_stats.F90:255-280`).
      !!
      !! Non-in-place — `local_list` and `global_list` must be distinct
      !! actual arguments, matching `halo_allreduce_sum`'s aliasing
      !! contract.
      type(efp_t), intent(in) :: local_list(:)
      type(efp_t), intent(out) :: global_list(:)
      integer, intent(in) :: nval

      type(comm_t) :: comm
      integer :: num_ranks, i
      real(real64) :: sendbuf(EFP_DIGITS*nval), recvbuf(EFP_DIGITS*nval)
      logical :: transport_ok

      comm = comm_env_compute_comm()
      num_ranks = comm%size()

      ! Single rank: the reduction is the identity, so return this rank's
      ! own contribution without entering a collective.  Not just an
      ! optimisation -- pic-mpi's serial backend (PIC_ENABLE_MPI=OFF)
      ! deliberately `error stop`s in `allreduce`, pushing the size()==1
      ! case onto the caller.  This IS that case.
      ! The transport-bound guards below police the cross-rank
      ! transport only, so there is nothing for them to check here;
      ! the copy is bit-identical to the serial path by construction.
      if (num_ranks == 1) then
         global_list(1:nval) = local_list(1:nval)
         return
      end if

      if (num_ranks > EFP_MAX_RANKS) then
         call logger%error("============================================")
         call logger%error("[panic] halo_allreduce_efp_list: num_ranks exceeds EFP_MAX_RANKS")
         call logger%error("============================================")
         error stop "halo_allreduce_efp_list: rank count exceeds the exact-transport bound"
      end if

      do i = 1, nval
         if (.not. efp_bin1_within_transport_bound(local_list(i), num_ranks)) then
            call logger%error("============================================")
            call logger%error("[panic] halo_allreduce_efp_list: EFP bin-1 transport bound violated")
            call logger%error("============================================")
            error stop "halo_allreduce_efp_list: bin-1 magnitude exceeds 2**53/num_ranks"
         end if
      end do

      call efp_to_transport(local_list(1:nval), sendbuf)
      call allreduce(comm, sendbuf, recvbuf, EFP_DIGITS*nval, op=MPI_SUM)

      call efp_from_transport(recvbuf, global_list(1:nval), transport_ok)
      if (.not. transport_ok) then
         call logger%error("============================================")
         call logger%error("[panic] halo_allreduce_efp_list: post-combine bin is not an exact integer")
         call logger%error("============================================")
         error stop "halo_allreduce_efp_list: exact-double transport invariant violated"
      end if

   end subroutine halo_allreduce_efp_list

end module rdb_halo
