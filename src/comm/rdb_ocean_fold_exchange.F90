!! Distributed tripolar north-fold exchange (`px > 1`) — MPI backend.
!!
!! ROUTING
!! =======
!! The routing is the pure plan of `rdb_ocean_fold_plan` (read its header
!! first): every north-ghost / fold-line value comes from the tile that OWNS
!! the mirror point, never from a halo copy, so the fold needs no preceding
!! x exchange and is correct even where no halo precedes it.  Only the north
!! rank row takes part (`decomp%has_north` and a folded north edge); every
!! rank of that row reaches the same call sites, so each call is collective
!! over exactly that row.  Other ranks never enter.
!!
!! MESSAGES
!! ========
!! One message per ordered (sender, receiver) pair per call; the self pair
!! is a local copy through the same buffers (no MPI).  A call may carry a
!! GROUP of fields: `ocean_fold_begin` → `ocean_fold_pack` per field →
!! `ocean_fold_exchange` → `ocean_fold_unpack` per field, in the SAME
!! order → `ocean_fold_end`.  Per-peer layout (field-major, then layer,
!! then row, then entry — the halo's layer-outer convention):
!!
!!   idx = (p-1)*cap + offT*n(p,T) + offU*n(p,U)
!!         + ((L-1)*nrow + (r-1))*n(p,fam) + e
!!
!! where `offT`/`offU` count the (row x layer) slabs already packed per
!! column family.  Both ends advance the same cursors over the same field
!! order, and the receiver's `n(p,·)` (its receive counts) equal the
!! sender's send counts for that pair, so no lengths are exchanged.
!!
!! SIGN AND THE FOLD-LINE ROW
!! ==========================
!! The sign is applied at UNPACK (`dst = -buf` for a true vector), so one
!! pack kernel serves every field; IEEE negation is exact, so where it is
!! applied cannot change bits.  For v and corners the fold-line row is
!! always sent (message row 1); the receiver writes it only on its WEST-half
!! columns, zeroes a self-conjugate column for a vector and leaves it for a
!! scalar, and never touches the EAST (authoritative) half — the serial
!! projection rule, as a COPY: one value moves, exactly, so no
!! floating-point expression has to agree across ranks.  If the projection
!! is ever changed to an AVERAGE, the operand order must be canonical (east
!! operand first) on every rank.
!!
!! No source of the fold is ever a destination (sources: rows below the
!! fold line and the east half of the fold row), so pack-all → MPI →
!! unpack-all is hazard-free in place.
!!
!! DEVICE vs HOST
!! ==============
!! Mirrors `rdb_ocean_halo` exactly: the `device_resident` switch selects
!! `!$acc parallel loop` pack/unpack over device-present arrays (default)
!! or plain host loops (setup-time host fills), and the MPI calls pass the
!! buffers through `!$acc host_data use_device` — like the halo, the device
!! path always hands MPI DEVICE addresses (`RDB_CUDA_AWARE_MPI` only selects
!! a log line today; a host-staged branch would be added to both modules
!! together).  The index lists are copied to the device once at init.
!!
!! TAG
!! ===
!!   TAG_OC_FOLD = 15  (halo uses 11-14, `rdb_halo` 1-4).  One message per
!!   ordered pair per call and every call completes with `waitall` before
!!   returning, so (source, tag) is unambiguous, px = 2 included.
#ifdef RDB_DOUBLE_PRECISION
#define FOLD_ISEND_N comm_isend_real_dp_array_n
#define FOLD_IRECV_N comm_irecv_real_dp_array_n
#else
#define FOLD_ISEND_N comm_isend_real_sp_array_n
#define FOLD_IRECV_N comm_irecv_real_sp_array_n
#endif
module rdb_ocean_fold_exchange
   !! Owner-routed north-fold exchange for an east-west split north rank
   !! row, plus the px-dispatching single-field fold entry points
   !! (`ocean_fold_north_*`): `px = 1` keeps the local kernels of
   !! `rdb_ocean_fold` verbatim, `px > 1` always uses the exchange (the
   !! self-mirror tile through a local copy).
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_rank_from_coords
   use pic_mpi_lib, only: comm_t, request_t, MPI_Status, waitall, &
                          FOLD_ISEND_N, FOLD_IRECV_N
   use rdb_comm_env, only: comm_env_compute_comm
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face, fold_north_corner
   use rdb_ocean_fold_plan, only: fold_plan_t, fold_plan_build, fold_stagger_family, &
                                  fold_stagger_nrows, FOLD_PLAN_OK, FOLD_FAM_T, &
                                  FOLD_FAM_U, FOLD_NFAM, FOLD_ROW_WEST, FOLD_ROW_SELF, &
                                  FOLD_STAG_T, FOLD_STAG_U, FOLD_STAG_V, FOLD_STAG_CORNER
   implicit none
   private

   public :: ocean_fold_exchange_init
   public :: ocean_fold_exchange_destroy
   public :: ocean_fold_exchange_reserve
   public :: ocean_fold_is_distributed
   public :: ocean_fold_begin
   public :: ocean_fold_pack
   public :: ocean_fold_exchange
   public :: ocean_fold_unpack
   public :: ocean_fold_end
   public :: ocean_fold_north_centre
   public :: ocean_fold_north_u_face
   public :: ocean_fold_north_v_face
   public :: ocean_fold_north_corner
   public :: FOLD_STAG_T, FOLD_STAG_U, FOLD_STAG_V, FOLD_STAG_CORNER

   ! Rank-generic group primitives: 2D or 3D by the rank of fld.
   interface ocean_fold_pack
      module procedure ocean_fold_pack_2d, ocean_fold_pack_3d
   end interface ocean_fold_pack

   interface ocean_fold_unpack
      module procedure ocean_fold_unpack_2d, ocean_fold_unpack_3d
   end interface ocean_fold_unpack

   ! px-dispatching single-field folds (local kernel at px = 1, exchange
   ! at px > 1).  Same stagger semantics as the `rdb_ocean_fold` kernels.
   interface ocean_fold_north_centre
      module procedure fold_centre_2d, fold_centre_3d
   end interface ocean_fold_north_centre

   interface ocean_fold_north_u_face
      module procedure fold_u_2d, fold_u_3d
   end interface ocean_fold_north_u_face

   interface ocean_fold_north_v_face
      module procedure fold_v_2d, fold_v_3d
   end interface ocean_fold_north_v_face

   interface ocean_fold_north_corner
      module procedure fold_corner_2d
   end interface ocean_fold_north_corner

   integer, parameter :: TAG_OC_FOLD = 15
      !! Fold message tag (halo: 11-14, `rdb_halo`: 1-4).

   ! -----------------------------------------------------------------
   ! Topology
   ! -----------------------------------------------------------------
   logical :: fx_initialised = .false.
      !! True after `ocean_fold_exchange_init`.
   logical :: fx_active = .false.
      !! This rank takes part in a distributed fold (px > 1, folded north
      !! edge, north rank row).
   integer :: fx_ng = 0
      !! Ghost width.
   integer :: fx_nxl = 0
      !! Tile physical width.
   integer :: fx_nyl = 0
      !! Tile physical height (identical on every north-row tile).
   integer :: fx_px = 1
      !! Tiles along x.
   integer :: fx_ry = 0
      !! This rank's y coordinate (the north row).

   ! -----------------------------------------------------------------
   ! Plan, flattened for the kernels (device-resident when active)
   ! -----------------------------------------------------------------
   type(fold_plan_t) :: fx_plan
      !! Host copy of the routing plan.
   integer :: fx_npeer = 0
      !! Peers (self included when it routes to itself).
   integer :: fx_self = 0
      !! Peer index of this tile (0 if not a peer).
   integer :: fx_nmax = 1
      !! Largest (peer, family) entry count.
   integer, allocatable :: fx_peer_rank(:)
      !! Compute-comm rank of each peer.
   integer, allocatable :: fx_sn(:, :)
      !! Send counts (npeer, FOLD_NFAM).
   integer, allocatable :: fx_rn(:, :)
      !! Receive counts (npeer, FOLD_NFAM).
   integer :: fx_s0(FOLD_NFAM) = 0
      !! Start of each family's block in the flat send arrays.
   integer :: fx_r0(FOLD_NFAM) = 0
      !! Start of each family's block in the flat receive arrays.
   integer :: fx_ns(FOLD_NFAM) = 0
      !! Send entries per family.
   integer :: fx_nr(FOLD_NFAM) = 0
      !! Receive entries per family.
   integer, allocatable :: fx_scol(:), fx_speer(:), fx_se(:)
      !! Send entries: source column, peer index, position in its list.
   integer, allocatable :: fx_rcol(:), fx_rpeer(:), fx_re(:), fx_rcls(:)
      !! Receive entries: destination column, peer, position, fold-row class.

   ! -----------------------------------------------------------------
   ! Buffers + group cursors
   ! -----------------------------------------------------------------
   real(wp), allocatable :: fx_sbuf(:), fx_rbuf(:)
      !! Send / receive buffers, `fx_cap` values per peer region.
   integer :: fx_cap = 0
      !! Per-peer capacity (values).
   integer :: fx_slab_cap = 0
      !! Per-peer capacity in (row x layer) slabs of `fx_nmax` entries.
   type(request_t), allocatable :: fx_reqs(:)
      !! MPI requests (2 per peer).
   type(MPI_Status), allocatable :: fx_stats(:)
      !! MPI statuses (2 per peer).
   integer :: fx_poff(FOLD_NFAM) = 0
      !! Pack cursor: slabs packed per family in the open group.
   integer :: fx_uoff(FOLD_NFAM) = 0
      !! Unpack cursor: slabs unpacked per family in the open group.
   logical :: fx_open = .false.
      !! A group is open (between begin and end).
   logical :: fx_sent = .false.
      !! The open group has been exchanged (packing closed).

contains

   ! ==================================================================
   ! Lifecycle
   ! ==================================================================

   subroutine ocean_fold_exchange_init(decomp, nghost, north_fold, ierr)
      !! Record the topology and, on a north-row rank of an east-west split
      !! folded grid, build the routing plan and map it to the device.
      !! Every rank may call it (it is not collective); a rank that does
      !! not fold, or a `px = 1` run, keeps the exchange inactive.
      type(decomp_t), intent(in) :: decomp
         !! Domain decomposition.
      integer, intent(in) :: nghost
         !! Ghost width.
      logical, intent(in) :: north_fold
         !! This rank folds its north edge (`bc%north_fold`, rank-local).
      integer, intent(out), optional :: ierr
         !! `OCEAN_STATUS_ERR_SETUP` when the plan cannot be built (more
         !! tiles than fold-row columns); absent ⇒ `error stop`.

      integer :: status, p, f

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (fx_initialised) call ocean_fold_exchange_destroy()

      fx_ng = nghost
      fx_nxl = decomp%nx_local
      fx_nyl = decomp%ny_local
      fx_px = decomp%px
      fx_ry = decomp%ry
      fx_active = north_fold .and. decomp%px > 1
      fx_initialised = .true.
      if (.not. fx_active) return

      call fold_plan_build(fx_plan, decomp%nx_global, decomp%px, nghost, decomp%rx, status)
      if (status /= FOLD_PLAN_OK) then
         call fail("ocean_fold_exchange_init: cannot build the fold plan for nx = "// &
                   to_string(decomp%nx_global)//", px = "//to_string(decomp%px)// &
                   ", nghost = "//to_string(nghost)//" (need nx >= px, nghost >= 1).", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         fx_active = .false.
         return
      end if

      fx_npeer = fx_plan%npeer
      fx_self = fx_plan%self_peer
      fx_nmax = max(fx_plan%nmax, 1)
      allocate (fx_peer_rank(fx_npeer))
      do p = 1, fx_npeer
         fx_peer_rank(p) = decomp_rank_from_coords(decomp%px, fx_plan%peer_rx(p), decomp%ry)
      end do
      allocate (fx_sn(fx_npeer, FOLD_NFAM), fx_rn(fx_npeer, FOLD_NFAM))
      fx_sn = fx_plan%send_n
      fx_rn = fx_plan%recv_n
      do f = 1, FOLD_NFAM
         fx_ns(f) = fx_plan%nsend(f)
         fx_nr(f) = fx_plan%nrecv(f)
         fx_s0(f) = 0
         fx_r0(f) = 0
         if (fx_npeer > 0) then
            fx_s0(f) = fx_plan%send_start(1, f)
            fx_r0(f) = fx_plan%recv_start(1, f)
         end if
      end do
      fx_scol = fx_plan%send_col
      fx_speer = fx_plan%send_peer
      fx_se = fx_plan%send_e
      fx_rcol = fx_plan%recv_col
      fx_rpeer = fx_plan%recv_peer
      fx_re = fx_plan%recv_e
      fx_rcls = fx_plan%recv_cls
      allocate (fx_reqs(2*max(fx_npeer, 1)), fx_stats(2*max(fx_npeer, 1)))
      !$acc enter data copyin(fx_sn, fx_rn, fx_scol, fx_speer, fx_se, &
      !$acc&                  fx_rcol, fx_rpeer, fx_re, fx_rcls)

      ! Base capacity: one field of the widest stagger, one layer.
      call grow_buffers(fx_ng + 1)
   end subroutine ocean_fold_exchange_init

   subroutine ocean_fold_exchange_destroy()
      !! Release the plan, buffers and topology.  Idempotent.
      if (.not. fx_initialised) return
      if (fx_active) then
         !$acc exit data delete(fx_sn, fx_rn, fx_scol, fx_speer, fx_se, &
         !$acc&                 fx_rcol, fx_rpeer, fx_re, fx_rcls)
         if (allocated(fx_sbuf)) then
            !$acc exit data delete(fx_sbuf, fx_rbuf)
            deallocate (fx_sbuf, fx_rbuf)
         end if
         deallocate (fx_peer_rank, fx_sn, fx_rn, fx_scol, fx_speer, fx_se, &
                     fx_rcol, fx_rpeer, fx_re, fx_rcls, fx_reqs, fx_stats)
      end if
      call fx_plan%destroy()
      fx_cap = 0
      fx_slab_cap = 0
      fx_npeer = 0
      fx_self = 0
      fx_nmax = 1
      fx_active = .false.
      fx_open = .false.
      fx_sent = .false.
      fx_initialised = .false.
   end subroutine ocean_fold_exchange_destroy

   subroutine ocean_fold_exchange_reserve(nslab)
      !! Pre-size the buffers for the largest group of the run, in
      !! (row x layer) slabs summed over its fields — e.g. the
      !! `ml_state` group is `(3 + ntracer)*nz` fields of at most `ng+1`
      !! rows — so no device reallocation fires mid-run.  No-op when
      !! inactive or already large enough.
      integer, intent(in) :: nslab
         !! Total (rows x layers) over the fields of the largest group.
      if (.not. fx_active) return
      if (nslab > fx_slab_cap) call grow_buffers(nslab)
   end subroutine ocean_fold_exchange_reserve

   pure function ocean_fold_is_distributed() result(flag)
      !! True iff this rank folds through the exchange (`px > 1` and a
      !! folded north edge on this rank).  False on every `px = 1` run and
      !! on ranks off the north row.
      logical :: flag
      flag = fx_initialised .and. fx_active
   end function ocean_fold_is_distributed

   subroutine grow_buffers(nslab)
      !! (Re)allocate both buffers for `nslab` slabs per peer
      !! (exit-delete / dealloc / alloc / enter-create, as
      !! `ocean_halo_buffers_ensure_nz`).  Never while a group is open.
      integer, intent(in) :: nslab

      if (fx_open) error stop "rdb_ocean_fold_exchange: buffer growth inside an open group"
      if (allocated(fx_sbuf)) then
         !$acc exit data delete(fx_sbuf, fx_rbuf)
         deallocate (fx_sbuf, fx_rbuf)
      end if
      fx_slab_cap = nslab
      fx_cap = nslab*fx_nmax
      allocate (fx_sbuf(fx_cap*max(fx_npeer, 1)), fx_rbuf(fx_cap*max(fx_npeer, 1)))
      !$acc enter data create(fx_sbuf, fx_rbuf)
   end subroutine grow_buffers

   ! ==================================================================
   ! Group protocol
   ! ==================================================================

   subroutine ocean_fold_begin(nslab)
      !! Open a group of `nslab` (rows x layers) slabs in total (an upper
      !! bound is fine).  Grows the buffers, with a warning, if the run did
      !! not reserve enough.  No-op when inactive.
      integer, intent(in) :: nslab
         !! Sum over the group's fields of rows x layers.

      if (.not. fx_active) return
      if (fx_open) error stop "rdb_ocean_fold_exchange: ocean_fold_begin with a group already open"
      if (nslab > fx_slab_cap) then
         call logger%warning("ocean_fold_begin: growing the fold buffers mid-run ("// &
                             to_string(fx_slab_cap)//" -> "//to_string(nslab)// &
                             " slabs). Call ocean_fold_exchange_reserve at init.")
         call grow_buffers(nslab)
      end if
      fx_poff = 0
      fx_uoff = 0
      fx_open = .true.
      fx_sent = .false.
   end subroutine ocean_fold_begin

   subroutine ocean_fold_end()
      !! Close the group: every packed slab must have been unpacked.
      if (.not. fx_active) return
      if (.not. fx_open) error stop "rdb_ocean_fold_exchange: ocean_fold_end without begin"
      if (any(fx_uoff /= fx_poff)) then
         error stop "rdb_ocean_fold_exchange: group closed with fields left unpacked"
      end if
      fx_open = .false.
      fx_sent = .false.
   end subroutine ocean_fold_end

   subroutine ocean_fold_pack_2d(fld, nxa, nya, stagger, device_resident)
      !! Pack a 2D field of storage shape (nxa, nya) into the open group.
      integer, intent(in) :: nxa, nya
         !! Storage extents of `fld`.
      real(wp), intent(in) :: fld(nxa, nya)
         !! Field (T: nxt x nyt; u: nxt+1 x nyt; v: nxt x nyt+1;
         !! corner: nxt+1 x nyt+1).
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in), optional :: device_resident
         !! `.false.` ⇒ host arrays; default device-resident.

      call ocean_fold_pack_3d(fld, nxa, nya, 1, stagger, device_resident)
   end subroutine ocean_fold_pack_2d

   subroutine ocean_fold_pack_3d(fld, nxa, nya, nz, stagger, device_resident)
      !! Pack every layer of a 3D field into the open group: the sender's
      !! owned mirror-source points, rows below (and, for v / corner, on)
      !! the fold line.
      integer, intent(in) :: nxa, nya, nz
         !! Storage extents of `fld`.
      real(wp), intent(in) :: fld(nxa, nya, nz)
         !! Field; see `ocean_fold_pack_2d` for the per-stagger shape.
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in), optional :: device_resident
         !! `.false.` ⇒ host arrays; default device-resident.

      integer :: fam, nrow, sbase, g, g0, ng_e, L, r, p, cap, oT, oU
      logical :: on_device

      if (.not. fx_active) return
      if (.not. fx_open .or. fx_sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_pack outside an open, unsent group"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      fam = fold_stagger_family(stagger)
      nrow = fold_stagger_nrows(stagger, fx_ng)
      if (sum(fx_poff) + nrow*nz > fx_slab_cap) then
         error stop "rdb_ocean_fold_exchange: group larger than its ocean_fold_begin size"
      end if
      ! Source row of message row r: ng+nyl+1-r (T, u) / ng+nyl+2-r (v,
      ! corner; r = 1 is the fold-line row) — `fold_row_map`.
      sbase = fx_ng + fx_nyl + 1
      if (nrow == fx_ng + 1) sbase = sbase + 1
      g0 = fx_s0(fam)
      ng_e = fx_ns(fam)
      cap = fx_cap
      oT = fx_poff(FOLD_FAM_T)
      oU = fx_poff(FOLD_FAM_U)

      if (ng_e > 0) then
         if (on_device) then
            !$acc parallel loop collapse(3) private(p) &
            !$acc& present(fx_sbuf, fx_scol, fx_speer, fx_se, fx_sn, fld)
            do L = 1, nz
               do r = 1, nrow
                  do g = g0 + 1, g0 + ng_e
                     p = fx_speer(g)
                     fx_sbuf((p - 1)*cap + oT*fx_sn(p, FOLD_FAM_T) + oU*fx_sn(p, FOLD_FAM_U) &
                             + ((L - 1)*nrow + (r - 1))*fx_sn(p, fam) + fx_se(g)) = &
                        fld(fx_scol(g), sbase - r, L)
                  end do
               end do
            end do
         else
            do L = 1, nz
               do r = 1, nrow
                  do g = g0 + 1, g0 + ng_e
                     p = fx_speer(g)
                     fx_sbuf((p - 1)*cap + oT*fx_sn(p, FOLD_FAM_T) + oU*fx_sn(p, FOLD_FAM_U) &
                             + ((L - 1)*nrow + (r - 1))*fx_sn(p, fam) + fx_se(g)) = &
                        fld(fx_scol(g), sbase - r, L)
                  end do
               end do
            end do
         end if
      end if
      fx_poff(fam) = fx_poff(fam) + nrow*nz
   end subroutine ocean_fold_pack_3d

   subroutine ocean_fold_exchange(device_resident)
      !! Move every packed message: one isend + irecv per non-self peer
      !! with a non-empty message, the self pair as a local copy, then
      !! `waitall`.  Collective over the north rank row.
      logical, intent(in), optional :: device_resident
         !! `.false.` ⇒ host buffers; default device-resident.

      integer :: p, n_s, n_r, o, nreq, i
      logical :: on_device
      type(comm_t) :: comm

      if (.not. fx_active) return
      if (.not. fx_open .or. fx_sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_exchange outside an open, unsent group"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident
      comm = comm_env_compute_comm()

      nreq = 0
      do p = 1, fx_npeer
         n_s = fx_poff(FOLD_FAM_T)*fx_sn(p, FOLD_FAM_T) + fx_poff(FOLD_FAM_U)*fx_sn(p, FOLD_FAM_U)
         n_r = fx_poff(FOLD_FAM_T)*fx_rn(p, FOLD_FAM_T) + fx_poff(FOLD_FAM_U)*fx_rn(p, FOLD_FAM_U)
         o = (p - 1)*fx_cap
         if (p == fx_self) then
            ! Self pair: the receive region is the send region (same list).
            if (on_device) then
               !$acc parallel loop present(fx_sbuf, fx_rbuf)
               do i = o + 1, o + n_s
                  fx_rbuf(i) = fx_sbuf(i)
               end do
            else
               fx_rbuf(o + 1:o + n_s) = fx_sbuf(o + 1:o + n_s)
            end if
            cycle
         end if
         if (on_device) then
            !$acc host_data use_device(fx_sbuf, fx_rbuf)
            if (n_r > 0) then
               nreq = nreq + 1
               call FOLD_IRECV_N(comm, fx_rbuf(o + 1:o + n_r), n_r, fx_peer_rank(p), &
                                 TAG_OC_FOLD, fx_reqs(nreq))
            end if
            if (n_s > 0) then
               nreq = nreq + 1
               call FOLD_ISEND_N(comm, fx_sbuf(o + 1:o + n_s), n_s, fx_peer_rank(p), &
                                 TAG_OC_FOLD, fx_reqs(nreq))
            end if
            !$acc end host_data
         else
            if (n_r > 0) then
               nreq = nreq + 1
               call FOLD_IRECV_N(comm, fx_rbuf(o + 1:o + n_r), n_r, fx_peer_rank(p), &
                                 TAG_OC_FOLD, fx_reqs(nreq))
            end if
            if (n_s > 0) then
               nreq = nreq + 1
               call FOLD_ISEND_N(comm, fx_sbuf(o + 1:o + n_s), n_s, fx_peer_rank(p), &
                                 TAG_OC_FOLD, fx_reqs(nreq))
            end if
         end if
      end do
      if (nreq > 0) call waitall(fx_reqs(1:nreq), fx_stats(1:nreq))
      fx_sent = .true.
   end subroutine ocean_fold_exchange

   subroutine ocean_fold_unpack_2d(fld, nxa, nya, stagger, negate, device_resident)
      !! Unpack the next field of the exchanged group into a 2D field.
      integer, intent(in) :: nxa, nya
         !! Storage extents of `fld`.
      real(wp), intent(inout) :: fld(nxa, nya)
         !! Field; same stagger and order as at pack.
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in) :: negate
         !! `.true.` for a true-vector component (sign flips across the fold).
      logical, intent(in), optional :: device_resident
         !! `.false.` ⇒ host arrays; default device-resident.

      call ocean_fold_unpack_3d(fld, nxa, nya, 1, stagger, negate, device_resident)
   end subroutine ocean_fold_unpack_2d

   subroutine ocean_fold_unpack_3d(fld, nxa, nya, nz, stagger, negate, device_resident)
      !! Unpack the next field of the exchanged group: write every north
      !! ghost row (all storage columns) and, for v / corner, the fold-line
      !! row's west half (self-conjugate column → 0 for a vector).
      integer, intent(in) :: nxa, nya, nz
         !! Storage extents of `fld`.
      real(wp), intent(inout) :: fld(nxa, nya, nz)
         !! Field; same stagger and order as at pack.
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in) :: negate
         !! `.true.` for a true-vector component.
      logical, intent(in), optional :: device_resident
         !! `.false.` ⇒ host arrays; default device-resident.

      integer :: fam, nrow, dbase, g, g0, ng_e, L, r, p, cap, oT, oU, ib, cls
      logical :: on_device, has_row
      real(wp) :: sgn, val

      if (.not. fx_active) return
      if (.not. fx_open .or. .not. fx_sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_unpack before ocean_fold_exchange"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      fam = fold_stagger_family(stagger)
      nrow = fold_stagger_nrows(stagger, fx_ng)
      has_row = (nrow == fx_ng + 1)
      ! Destination row of message row r: ng+nyl+r for every stagger
      ! (`fold_row_map`: T/u ng+nyl+d, v/corner ng+nyl+1+d with d = r-1).
      dbase = fx_ng + fx_nyl
      g0 = fx_r0(fam)
      ng_e = fx_nr(fam)
      cap = fx_cap
      oT = fx_uoff(FOLD_FAM_T)
      oU = fx_uoff(FOLD_FAM_U)
      sgn = 1.0_wp
      if (negate) sgn = -1.0_wp

      if (ng_e > 0) then
         if (on_device) then
            !$acc parallel loop collapse(3) private(p, ib, cls, val) &
            !$acc& present(fx_rbuf, fx_rcol, fx_rpeer, fx_re, fx_rcls, fx_rn, fld)
            do L = 1, nz
               do r = 1, nrow
                  do g = g0 + 1, g0 + ng_e
                     p = fx_rpeer(g)
                     ib = (p - 1)*cap + oT*fx_rn(p, FOLD_FAM_T) + oU*fx_rn(p, FOLD_FAM_U) &
                          + ((L - 1)*nrow + (r - 1))*fx_rn(p, fam) + fx_re(g)
                     val = sgn*fx_rbuf(ib)
                     if (has_row .and. r == 1) then
                        cls = fx_rcls(g)
                        if (cls == FOLD_ROW_WEST) then
                           fld(fx_rcol(g), dbase + 1, L) = val
                        else if (cls == FOLD_ROW_SELF .and. negate) then
                           fld(fx_rcol(g), dbase + 1, L) = 0.0_wp
                        end if
                     else
                        fld(fx_rcol(g), dbase + r, L) = val
                     end if
                  end do
               end do
            end do
         else
            do L = 1, nz
               do r = 1, nrow
                  do g = g0 + 1, g0 + ng_e
                     p = fx_rpeer(g)
                     ib = (p - 1)*cap + oT*fx_rn(p, FOLD_FAM_T) + oU*fx_rn(p, FOLD_FAM_U) &
                          + ((L - 1)*nrow + (r - 1))*fx_rn(p, fam) + fx_re(g)
                     val = sgn*fx_rbuf(ib)
                     if (has_row .and. r == 1) then
                        cls = fx_rcls(g)
                        if (cls == FOLD_ROW_WEST) then
                           fld(fx_rcol(g), dbase + 1, L) = val
                        else if (cls == FOLD_ROW_SELF .and. negate) then
                           fld(fx_rcol(g), dbase + 1, L) = 0.0_wp
                        end if
                     else
                        fld(fx_rcol(g), dbase + r, L) = val
                     end if
                  end do
               end do
            end do
         end if
      end if
      fx_uoff(fam) = fx_uoff(fam) + nrow*nz
   end subroutine ocean_fold_unpack_3d

   ! ==================================================================
   ! px-dispatching single-field folds
   ! ==================================================================
   ! Each takes the field's storage extents and the tile geometry the
   ! local kernel needs.  px = 1 (or an inactive rank): the `rdb_ocean_fold`
   ! kernel, verbatim.  px > 1: a one-field group through the exchange.

   subroutine fold_centre_2d(fld, nxt, nyt, nx_phys, ny_phys, ng, device_resident)
      !! Cell-centred 2D north fold (copy).
      integer, intent(in) :: nxt, nyt, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxt, nyt)
         !! Cell-centred field (nx_total, ny_total).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_T, fx_ng))
         call ocean_fold_pack(fld, nxt, nyt, FOLD_STAG_T, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxt, nyt, FOLD_STAG_T, .false., device_resident)
         call ocean_fold_end()
      else
         call fold_north_centre(fld, nxt, nyt, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_centre_2d

   subroutine fold_centre_3d(fld, nxt, nyt, nz, nx_phys, ny_phys, ng, device_resident)
      !! Cell-centred 3D north fold (copy).
      integer, intent(in) :: nxt, nyt, nz, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxt, nyt, nz)
         !! Cell-centred field (nx_total, ny_total, nz).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_T, fx_ng)*nz)
         call ocean_fold_pack(fld, nxt, nyt, nz, FOLD_STAG_T, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxt, nyt, nz, FOLD_STAG_T, .false., device_resident)
         call ocean_fold_end()
      else
         call fold_north_centre(fld, nxt, nyt, nz, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_centre_3d

   subroutine fold_u_2d(fld, nxf, nyt, nx_phys, ny_phys, ng, device_resident)
      !! x-face (u) 2D north fold (negate).
      integer, intent(in) :: nxf, nyt, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxf, nyt)
         !! x-face field (nx_total+1, ny_total).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_U, fx_ng))
         call ocean_fold_pack(fld, nxf, nyt, FOLD_STAG_U, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxf, nyt, FOLD_STAG_U, .true., device_resident)
         call ocean_fold_end()
      else
         call fold_north_u_face(fld, nxf, nyt, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_u_2d

   subroutine fold_u_3d(fld, nxf, nyt, nz, nx_phys, ny_phys, ng, device_resident)
      !! x-face (u) 3D north fold (negate).
      integer, intent(in) :: nxf, nyt, nz, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxf, nyt, nz)
         !! x-face field (nx_total+1, ny_total, nz).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_U, fx_ng)*nz)
         call ocean_fold_pack(fld, nxf, nyt, nz, FOLD_STAG_U, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxf, nyt, nz, FOLD_STAG_U, .true., device_resident)
         call ocean_fold_end()
      else
         call fold_north_u_face(fld, nxf, nyt, nz, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_u_3d

   subroutine fold_v_2d(fld, nxt, nyf, nx_phys, ny_phys, ng, device_resident)
      !! y-face (v) 2D north fold (negate + fold-line projection).
      integer, intent(in) :: nxt, nyf, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxt, nyf)
         !! y-face field (nx_total, ny_total+1).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_V, fx_ng))
         call ocean_fold_pack(fld, nxt, nyf, FOLD_STAG_V, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxt, nyf, FOLD_STAG_V, .true., device_resident)
         call ocean_fold_end()
      else
         call fold_north_v_face(fld, nxt, nyf, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_v_2d

   subroutine fold_v_3d(fld, nxt, nyf, nz, nx_phys, ny_phys, ng, device_resident)
      !! y-face (v) 3D north fold (negate + fold-line projection).
      integer, intent(in) :: nxt, nyf, nz, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxt, nyf, nz)
         !! y-face field (nx_total, ny_total+1, nz).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_V, fx_ng)*nz)
         call ocean_fold_pack(fld, nxt, nyf, nz, FOLD_STAG_V, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxt, nyf, nz, FOLD_STAG_V, .true., device_resident)
         call ocean_fold_end()
      else
         call fold_north_v_face(fld, nxt, nyf, nz, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_v_3d

   subroutine fold_corner_2d(fld, nxf, nyf, nx_phys, ny_phys, ng, negate, device_resident)
      !! SW-corner 2D north fold (+ fold-line projection); `negate` for a
      !! true-vector component, copy for a scalar / vorticity.
      integer, intent(in) :: nxf, nyf, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxf, nyf)
         !! Corner field (nx_total+1, ny_total+1).
      logical, intent(in) :: negate
         !! `.true.` ⇒ negate (true vector); `.false.` ⇒ copy.
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` ⇒ host arrays.

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_CORNER, fx_ng))
         call ocean_fold_pack(fld, nxf, nyf, FOLD_STAG_CORNER, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxf, nyf, FOLD_STAG_CORNER, negate, device_resident)
         call ocean_fold_end()
      else
         call fold_north_corner(fld, nxf, nyf, nx_phys, ny_phys, ng, negate)
      end if
   end subroutine fold_corner_2d

end module rdb_ocean_fold_exchange
