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
!! `do concurrent` pack/unpack over device-present arrays (default)
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
   !!
   !! TWO CONTEXTS.  The routing plan, its device index lists and the
   !! message buffers depend on the ghost width, so the exchange keeps one
   !! context per width: `FOLD_CTX_NORMAL` (`nghost`, every state-level
   !! fold) and `FOLD_CTX_WIDE` (`nghost + bt_halo`, the wide-halo
   !! barotropic march-in, built only when `bt_halo > 0`).  The wide fold is
   !! the SAME mechanism on a plan built for the wider band (plan
   !! `fold_plan_build(ni, px, nghost + bt_halo, rx)`): every value still
   !! comes from the tile that owns the mirror point.  The group protocol
   !! takes an optional `ctx` (default normal); a group is opened, packed,
   !! exchanged, unpacked and closed on ONE context.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_mem_report, only: arr_bytes
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
   public :: ocean_fold_exchange_grow_count
   public :: ocean_fold_exchange_bytes
   public :: ocean_fold_is_distributed
   public :: ocean_fold_ctx_is_distributed
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

   integer, parameter, public :: FOLD_CTX_NORMAL = 1
      !! Context of the `nghost`-wide plan (every state-level fold).
   integer, parameter, public :: FOLD_CTX_WIDE = 2
      !! Context of the `nghost + bt_halo`-wide plan (wide-halo march-in).
   integer, parameter :: NCTX = 2

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
      !! Fold message tag (halo: 11-14, `rdb_halo`: 1-4).  Both contexts
      !! share it: one group is in flight at a time and every call completes
      !! with `waitall` before returning.

   type :: fx_ctx_t
      !! Everything that depends on the ghost width: topology, plan, device
      !! index lists, buffers and the open group's cursors.
      logical :: initialised = .false.
         !! True after the context was initialised.
      logical :: active = .false.
         !! This rank takes part in a distributed fold (px > 1, folded north
         !! edge, north rank row).
      integer :: ng = 0
         !! Ghost width of the plan.
      integer :: nxl = 0
         !! Tile physical width.
      integer :: nyl = 0
         !! Tile physical height (identical on every north-row tile).
      integer :: px = 1
         !! Tiles along x.
      integer :: ry = 0
         !! This rank's y coordinate (the north row).
      type(fold_plan_t) :: plan
         !! Host copy of the routing plan.
      integer :: npeer = 0
         !! Peers (self included when it routes to itself).
      integer :: self = 0
         !! Peer index of this tile (0 if not a peer).
      integer :: nmax = 1
         !! Largest (peer, family) entry count.
      integer, allocatable :: peer_rank(:)
         !! Compute-comm rank of each peer.
      integer, allocatable :: sn(:, :)
         !! Send counts (npeer, FOLD_NFAM).
      integer, allocatable :: rn(:, :)
         !! Receive counts (npeer, FOLD_NFAM).
      integer :: s0(FOLD_NFAM) = 0
         !! Start of each family's block in the flat send arrays.
      integer :: r0(FOLD_NFAM) = 0
         !! Start of each family's block in the flat receive arrays.
      integer :: ns(FOLD_NFAM) = 0
         !! Send entries per family.
      integer :: nr(FOLD_NFAM) = 0
         !! Receive entries per family.
      integer, allocatable :: scol(:), speer(:), se(:)
         !! Send entries: source column, peer index, position in its list.
      integer, allocatable :: rcol(:), rpeer(:), re(:), rcls(:)
         !! Receive entries: destination column, peer, position, fold-row class.
      real(wp), allocatable :: sbuf(:), rbuf(:)
         !! Send / receive buffers, `cap` values per peer region.
      integer :: cap = 0
         !! Per-peer capacity (values).
      integer :: grow_count = 0
         !! Mid-run regrows (`ocean_fold_begin` finding the reserve too small).
      integer :: slab_cap = 0
         !! Per-peer capacity in (row x layer) slabs of `nmax` entries.
      type(request_t), allocatable :: reqs(:)
         !! MPI requests (2 per peer).
      type(MPI_Status), allocatable :: stats(:)
         !! MPI statuses (2 per peer).
      integer :: poff(FOLD_NFAM) = 0
         !! Pack cursor: slabs packed per family in the open group.
      integer :: uoff(FOLD_NFAM) = 0
         !! Unpack cursor: slabs unpacked per family in the open group.
      logical :: open = .false.
         !! A group is open (between begin and end).
      logical :: sent = .false.
         !! The open group has been exchanged (packing closed).
   contains
      procedure :: bytes => fx_ctx_bytes
   end type fx_ctx_t

   type(fx_ctx_t), target, save :: fxc(NCTX)
      !! The contexts, indexed by `FOLD_CTX_*`.

contains

   ! ==================================================================
   ! Lifecycle
   ! ==================================================================

   subroutine ocean_fold_exchange_init(decomp, nghost, north_fold, ierr, nghost_wide)
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
         !! tiles than fold-row columns); absent => `error stop`.
      integer, intent(in), optional :: nghost_wide
         !! `nghost + bt_halo` of the wide-halo march-in: also builds the
         !! `FOLD_CTX_WIDE` plan.  Absent or `<= nghost`: no wide context.

      integer :: ierr_w

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      call ocean_fold_exchange_destroy()
      call init_ctx(fxc(FOLD_CTX_NORMAL), decomp, nghost, north_fold, ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      if (present(nghost_wide)) then
         if (nghost_wide > nghost) then
            ierr_w = OCEAN_STATUS_OK
            call init_ctx(fxc(FOLD_CTX_WIDE), decomp, nghost_wide, north_fold, ierr_w)
            if (present(ierr)) ierr = ierr_w
         end if
      end if
   end subroutine ocean_fold_exchange_init

   subroutine init_ctx(c, decomp, nghost, north_fold, ierr)
      !! Initialise one context (see `ocean_fold_exchange_init`).
      type(fx_ctx_t), intent(inout), target :: c
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nghost
      logical, intent(in) :: north_fold
      integer, intent(inout), optional :: ierr

      integer :: status, p, f

      c%grow_count = 0
      c%ng = nghost
      c%nxl = decomp%nx_local
      c%nyl = decomp%ny_local
      c%px = decomp%px
      c%ry = decomp%ry
      c%active = north_fold .and. decomp%px > 1
      c%initialised = .true.
      if (.not. c%active) return

      call fold_plan_build(c%plan, decomp%nx_global, decomp%px, nghost, decomp%rx, status)
      if (status /= FOLD_PLAN_OK) then
         call fail("ocean_fold_exchange_init: cannot build the fold plan for nx = "// &
                   to_string(decomp%nx_global)//", px = "//to_string(decomp%px)// &
                   ", nghost = "//to_string(nghost)//" (need nx >= px, nghost >= 1).", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         c%active = .false.
         return
      end if

      c%npeer = c%plan%npeer
      c%self = c%plan%self_peer
      c%nmax = max(c%plan%nmax, 1)
      allocate (c%peer_rank(c%npeer))
      do p = 1, c%npeer
         c%peer_rank(p) = decomp_rank_from_coords(decomp%px, c%plan%peer_rx(p), decomp%ry)
      end do
      allocate (c%sn(c%npeer, FOLD_NFAM), c%rn(c%npeer, FOLD_NFAM))
      c%sn = c%plan%send_n
      c%rn = c%plan%recv_n
      do f = 1, FOLD_NFAM
         c%ns(f) = c%plan%nsend(f)
         c%nr(f) = c%plan%nrecv(f)
         c%s0(f) = 0
         c%r0(f) = 0
         if (c%npeer > 0) then
            c%s0(f) = c%plan%send_start(1, f)
            c%r0(f) = c%plan%recv_start(1, f)
         end if
      end do
      c%scol = c%plan%send_col
      c%speer = c%plan%send_peer
      c%se = c%plan%send_e
      c%rcol = c%plan%recv_col
      c%rpeer = c%plan%recv_peer
      c%re = c%plan%recv_e
      c%rcls = c%plan%recv_cls
      allocate (c%reqs(2*max(c%npeer, 1)), c%stats(2*max(c%npeer, 1)))
      !$acc enter data copyin(c%sn, c%rn, c%scol, c%speer, c%se, &
      !$acc&                  c%rcol, c%rpeer, c%re, c%rcls)

      ! Base capacity: one field of the widest stagger, one layer.
      call grow_buffers(c, c%ng + 1)
   end subroutine init_ctx

   subroutine ocean_fold_exchange_destroy()
      !! Release the plans, buffers and topology of every context.  Idempotent.
      integer :: ic
      do ic = 1, NCTX
         call destroy_ctx(fxc(ic))
      end do
   end subroutine ocean_fold_exchange_destroy

   subroutine destroy_ctx(c)
      !! Release one context.
      type(fx_ctx_t), intent(inout), target :: c
      if (.not. c%initialised) return
      if (c%active) then
         !$acc exit data delete(c%sn, c%rn, c%scol, c%speer, c%se, &
         !$acc&                 c%rcol, c%rpeer, c%re, c%rcls)
         if (allocated(c%sbuf)) then
            !$acc exit data delete(c%sbuf, c%rbuf)
            deallocate (c%sbuf, c%rbuf)
         end if
         deallocate (c%peer_rank, c%sn, c%rn, c%scol, c%speer, c%se, &
                     c%rcol, c%rpeer, c%re, c%rcls, c%reqs, c%stats)
      end if
      call c%plan%destroy()
      c%cap = 0
      c%slab_cap = 0
      c%npeer = 0
      c%self = 0
      c%nmax = 1
      c%active = .false.
      c%open = .false.
      c%sent = .false.
      c%initialised = .false.
   end subroutine destroy_ctx

   subroutine ocean_fold_exchange_reserve(nslab, ctx)
      !! Pre-size the buffers for the largest group of the run, in
      !! (row x layer) slabs summed over its fields -- e.g. the
      !! `ml_state` group is `(3 + ntracer)*nz` fields of at most `ng+1`
      !! rows -- so no device reallocation fires mid-run.  No-op when
      !! inactive or already large enough.
      integer, intent(in) :: nslab
         !! Total (rows x layers) over the fields of the largest group.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).
      integer :: ic
      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (nslab > fxc(ic)%slab_cap) call grow_buffers(fxc(ic), nslab)
   end subroutine ocean_fold_exchange_reserve

   pure function ocean_fold_exchange_grow_count() result(n)
      !! Number of mid-run buffer regrows since init (an under-sized
      !! `ocean_fold_exchange_reserve`), over both contexts.  A correctly
      !! reserved run keeps this at 0, so a test can turn the regrow warning
      !! into a failure.
      integer :: n
      n = sum(fxc(:)%grow_count)
   end function ocean_fold_exchange_grow_count

   pure function fx_ctx_bytes(this) result(nbytes)
      !! Counted allocatable footprint of one context (host index lists and
      !! buffers; the device copies are the same sizes).  Zero when inactive.
      class(fx_ctx_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%peer_rank) &
               + arr_bytes(this%sn) &
               + arr_bytes(this%rn) &
               + arr_bytes(this%scol) &
               + arr_bytes(this%speer) &
               + arr_bytes(this%se) &
               + arr_bytes(this%rcol) &
               + arr_bytes(this%rpeer) &
               + arr_bytes(this%re) &
               + arr_bytes(this%rcls) &
               + arr_bytes(this%sbuf) &
               + arr_bytes(this%rbuf)
   end function fx_ctx_bytes

   pure function ocean_fold_exchange_bytes() result(nbytes)
      !! Counted footprint of the fold exchange (both contexts: the plans'
      !! index lists and the message buffers).
      integer(int64) :: nbytes
      integer :: ic
      nbytes = 0_int64
      do ic = 1, NCTX
         nbytes = nbytes + fxc(ic)%bytes()
      end do
   end function ocean_fold_exchange_bytes

   pure function ocean_fold_is_distributed() result(flag)
      !! True iff this rank folds through the exchange (`px > 1` and a
      !! folded north edge on this rank).  False on every `px = 1` run and
      !! on ranks off the north row.
      logical :: flag
      flag = fxc(FOLD_CTX_NORMAL)%initialised .and. fxc(FOLD_CTX_NORMAL)%active
   end function ocean_fold_is_distributed

   pure function ocean_fold_ctx_is_distributed(ctx) result(flag)
      !! `ocean_fold_is_distributed` for one routing context: true iff that
      !! context's plan was built and this rank folds through it.  The
      !! wide-halo march-in asks it for `FOLD_CTX_WIDE`.
      integer, intent(in) :: ctx
         !! `FOLD_CTX_NORMAL` or `FOLD_CTX_WIDE`.
      logical :: flag
      flag = fxc(ctx)%initialised .and. fxc(ctx)%active
   end function ocean_fold_ctx_is_distributed

   subroutine grow_buffers(c, nslab)
      !! (Re)allocate both buffers for `nslab` slabs per peer
      !! (exit-delete / dealloc / alloc / enter-create, as
      !! `ocean_halo_buffers_ensure_nz`).  Never while a group is open.
      type(fx_ctx_t), intent(inout), target :: c
      integer, intent(in) :: nslab

      if (c%open) error stop "rdb_ocean_fold_exchange: buffer growth inside an open group"
      if (allocated(c%sbuf)) then
         !$acc exit data delete(c%sbuf, c%rbuf)
         deallocate (c%sbuf, c%rbuf)
      end if
      c%slab_cap = nslab
      c%cap = nslab*c%nmax
      allocate (c%sbuf(c%cap*max(c%npeer, 1)), c%rbuf(c%cap*max(c%npeer, 1)))
      !$acc enter data create(c%sbuf, c%rbuf)
   end subroutine grow_buffers

   ! ==================================================================
   ! Group protocol
   ! ==================================================================

   subroutine ocean_fold_begin(nslab, ctx)
      !! Open a group of `nslab` (rows x layers) slabs in total (an upper
      !! bound is fine).  Grows the buffers, with a warning, if the run did
      !! not reserve enough.  No-op when inactive.
      integer, intent(in) :: nslab
         !! Sum over the group's fields of rows x layers.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).
      integer :: ic

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (fxc(ic)%open) error stop "rdb_ocean_fold_exchange: ocean_fold_begin with a group already open"
      if (nslab > fxc(ic)%slab_cap) then
         call logger%warning("ocean_fold_begin: growing the fold buffers mid-run ("// &
                             to_string(fxc(ic)%slab_cap)//" -> "//to_string(nslab)// &
                             " slabs). Call ocean_fold_exchange_reserve at init.")
         fxc(ic)%grow_count = fxc(ic)%grow_count + 1
         call grow_buffers(fxc(ic), nslab)
      end if
      fxc(ic)%poff = 0
      fxc(ic)%uoff = 0
      fxc(ic)%open = .true.
      fxc(ic)%sent = .false.
   end subroutine ocean_fold_begin

   subroutine ocean_fold_end(ctx)
      !! Close the group: every packed slab must have been unpacked.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).
      integer :: ic

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (.not. fxc(ic)%open) error stop "rdb_ocean_fold_exchange: ocean_fold_end without begin"
      if (any(fxc(ic)%uoff /= fxc(ic)%poff)) then
         error stop "rdb_ocean_fold_exchange: group closed with fields left unpacked"
      end if
      fxc(ic)%open = .false.
      fxc(ic)%sent = .false.
   end subroutine ocean_fold_end

   subroutine ocean_fold_pack_2d(fld, nxa, nya, stagger, device_resident, ctx)
      !! Pack a 2D field of storage shape (nxa, nya) into the open group.
      integer, intent(in) :: nxa, nya
         !! Storage extents of `fld`.
      real(wp), intent(in) :: fld(nxa, nya)
         !! Field (T: nxt x nyt; u: nxt+1 x nyt; v: nxt x nyt+1;
         !! corner: nxt+1 x nyt+1).
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in), optional :: device_resident
         !! `.false.` => host arrays; default device-resident.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).

      call ocean_fold_pack_3d(fld, nxa, nya, 1, stagger, device_resident, ctx)
   end subroutine ocean_fold_pack_2d

   subroutine ocean_fold_pack_3d(fld, nxa, nya, nz, stagger, device_resident, ctx)
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
         !! `.false.` => host arrays; default device-resident.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).

      integer :: ic, fam, nrow, sbase, g0, ng_e, cap, oT, oU
      logical :: on_device

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (.not. fxc(ic)%open .or. fxc(ic)%sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_pack outside an open, unsent group"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      fam = fold_stagger_family(stagger)
      nrow = fold_stagger_nrows(stagger, fxc(ic)%ng)
      if (sum(fxc(ic)%poff) + nrow*nz > fxc(ic)%slab_cap) then
         error stop "rdb_ocean_fold_exchange: group larger than its ocean_fold_begin size"
      end if
      ! Source row of message row r: ng+nyl+1-r (T, u) / ng+nyl+2-r (v,
      ! corner; r = 1 is the fold-line row) -- `fold_row_map`.
      sbase = fxc(ic)%ng + fxc(ic)%nyl + 1
      if (nrow == fxc(ic)%ng + 1) sbase = sbase + 1
      g0 = fxc(ic)%s0(fam)
      ng_e = fxc(ic)%ns(fam)
      cap = fxc(ic)%cap
      oT = fxc(ic)%poff(FOLD_FAM_T)
      oU = fxc(ic)%poff(FOLD_FAM_U)

      if (ng_e > 0) then
         if (on_device) then
            call pack_dev(fld, nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam, &
                          fxc(ic)%npeer, size(fxc(ic)%scol), size(fxc(ic)%sbuf), &
                          fxc(ic)%scol, fxc(ic)%speer, fxc(ic)%se, fxc(ic)%sn, &
                          fxc(ic)%sbuf)
         else
            call pack_host(fld, nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam, &
                           fxc(ic)%npeer, size(fxc(ic)%scol), size(fxc(ic)%sbuf), &
                           fxc(ic)%scol, fxc(ic)%speer, fxc(ic)%se, fxc(ic)%sn, &
                           fxc(ic)%sbuf)
         end if
      end if
      fxc(ic)%poff(fam) = fxc(ic)%poff(fam) + nrow*nz
   end subroutine ocean_fold_pack_3d

   subroutine pack_dev(fld, nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam, &
                       npeer, nsend, nbuf, scol, speer, se, sn, sbuf)
      !! Device pack kernel (explicit-shape dummies; the context's index
      !! lists are mapped once at init).
      integer, intent(in) :: nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam
      integer, intent(in) :: npeer, nsend, nbuf
      real(wp), intent(in) :: fld(nxa, nya, nz)
      integer, intent(in) :: scol(nsend), speer(nsend), se(nsend), sn(npeer, FOLD_NFAM)
      real(wp), intent(inout) :: sbuf(nbuf)
      integer :: g, L, r, p

      do concurrent(L=1:nz, r=1:nrow, g=g0 + 1:g0 + ng_e) local(p)
         p = speer(g)
         sbuf((p - 1)*cap + oT*sn(p, FOLD_FAM_T) + oU*sn(p, FOLD_FAM_U) &
              + ((L - 1)*nrow + (r - 1))*sn(p, fam) + se(g)) = &
            fld(scol(g), sbase - r, L)
      end do
   end subroutine pack_dev

   subroutine pack_host(fld, nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam, &
                        npeer, nsend, nbuf, scol, speer, se, sn, sbuf)
      !! Host twin of `pack_dev` (setup-time host arrays).
      integer, intent(in) :: nxa, nya, nz, nrow, sbase, g0, ng_e, cap, oT, oU, fam
      integer, intent(in) :: npeer, nsend, nbuf
      real(wp), intent(in) :: fld(nxa, nya, nz)
      integer, intent(in) :: scol(nsend), speer(nsend), se(nsend), sn(npeer, FOLD_NFAM)
      real(wp), intent(inout) :: sbuf(nbuf)
      integer :: g, L, r, p

      do L = 1, nz
         do r = 1, nrow
            do g = g0 + 1, g0 + ng_e
               p = speer(g)
               sbuf((p - 1)*cap + oT*sn(p, FOLD_FAM_T) + oU*sn(p, FOLD_FAM_U) &
                    + ((L - 1)*nrow + (r - 1))*sn(p, fam) + se(g)) = &
                  fld(scol(g), sbase - r, L)
            end do
         end do
      end do
   end subroutine pack_host

   subroutine copy_dev(n, o, sbuf, rbuf, nbuf)
      !! Self pair on the device: the receive region is the send region.
      integer, intent(in) :: n, o, nbuf
      real(wp), intent(in) :: sbuf(nbuf)
      real(wp), intent(inout) :: rbuf(nbuf)
      integer :: i
      do concurrent(i=o + 1:o + n)
         rbuf(i) = sbuf(i)
      end do
   end subroutine copy_dev

   subroutine ocean_fold_exchange(device_resident, ctx)
      !! Move every packed message: one isend + irecv per non-self peer
      !! with a non-empty message, the self pair as a local copy, then
      !! `waitall`.  Collective over the north rank row.
      logical, intent(in), optional :: device_resident
         !! `.false.` => host buffers; default device-resident.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).

      integer :: ic, p, n_s, n_r, o, nreq
      logical :: on_device
      type(comm_t) :: comm

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (.not. fxc(ic)%open .or. fxc(ic)%sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_exchange outside an open, unsent group"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident
      comm = comm_env_compute_comm()

      nreq = 0
      do p = 1, fxc(ic)%npeer
         n_s = fxc(ic)%poff(FOLD_FAM_T)*fxc(ic)%sn(p, FOLD_FAM_T) + &
               fxc(ic)%poff(FOLD_FAM_U)*fxc(ic)%sn(p, FOLD_FAM_U)
         n_r = fxc(ic)%poff(FOLD_FAM_T)*fxc(ic)%rn(p, FOLD_FAM_T) + &
               fxc(ic)%poff(FOLD_FAM_U)*fxc(ic)%rn(p, FOLD_FAM_U)
         o = (p - 1)*fxc(ic)%cap
         if (p == fxc(ic)%self) then
            ! Self pair: the receive region is the send region (same list).
            if (on_device) then
               call copy_dev(n_s, o, fxc(ic)%sbuf, fxc(ic)%rbuf, size(fxc(ic)%sbuf))
            else
               fxc(ic)%rbuf(o + 1:o + n_s) = fxc(ic)%sbuf(o + 1:o + n_s)
            end if
            cycle
         end if
         if (on_device) then
            ! Out of line: `host_data use_device` must name plain (dummy)
            ! arrays, not subscripted derived-type members.
            call post_dev(comm, fxc(ic)%sbuf, fxc(ic)%rbuf, size(fxc(ic)%sbuf), o, n_s, n_r, &
                          fxc(ic)%peer_rank(p), size(fxc(ic)%reqs), fxc(ic)%reqs, nreq)
         else
            if (n_r > 0) then
               nreq = nreq + 1
               call FOLD_IRECV_N(comm, fxc(ic)%rbuf(o + 1:o + n_r), n_r, fxc(ic)%peer_rank(p), &
                                 TAG_OC_FOLD, fxc(ic)%reqs(nreq))
            end if
            if (n_s > 0) then
               nreq = nreq + 1
               call FOLD_ISEND_N(comm, fxc(ic)%sbuf(o + 1:o + n_s), n_s, fxc(ic)%peer_rank(p), &
                                 TAG_OC_FOLD, fxc(ic)%reqs(nreq))
            end if
         end if
      end do
      if (nreq > 0) call waitall(fxc(ic)%reqs(1:nreq), fxc(ic)%stats(1:nreq))
      fxc(ic)%sent = .true.
   end subroutine ocean_fold_exchange

   subroutine post_dev(comm, sbuf, rbuf, nbuf, o, n_s, n_r, peer, nreqmax, reqs, nreq)
      !! Post one peer's receive and send on DEVICE buffers (CUDA-aware MPI
      !! gets the device addresses through `host_data use_device`).
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: nbuf, o, n_s, n_r, peer, nreqmax
      real(wp), intent(inout) :: sbuf(nbuf), rbuf(nbuf)
      type(request_t), intent(inout) :: reqs(nreqmax)
      integer, intent(inout) :: nreq

      !$acc host_data use_device(sbuf, rbuf)
      if (n_r > 0) then
         nreq = nreq + 1
         call FOLD_IRECV_N(comm, rbuf(o + 1:o + n_r), n_r, peer, TAG_OC_FOLD, reqs(nreq))
      end if
      if (n_s > 0) then
         nreq = nreq + 1
         call FOLD_ISEND_N(comm, sbuf(o + 1:o + n_s), n_s, peer, TAG_OC_FOLD, reqs(nreq))
      end if
      !$acc end host_data
   end subroutine post_dev

   subroutine ocean_fold_unpack_2d(fld, nxa, nya, stagger, negate, device_resident, ctx)
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
         !! `.false.` => host arrays; default device-resident.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).

      call ocean_fold_unpack_3d(fld, nxa, nya, 1, stagger, negate, device_resident, ctx)
   end subroutine ocean_fold_unpack_2d

   subroutine ocean_fold_unpack_3d(fld, nxa, nya, nz, stagger, negate, device_resident, ctx)
      !! Unpack the next field of the exchanged group: write every north
      !! ghost row (all storage columns) and, for v / corner, the fold-line
      !! row's west half (self-conjugate column -> 0 for a vector).
      integer, intent(in) :: nxa, nya, nz
         !! Storage extents of `fld`.
      real(wp), intent(inout) :: fld(nxa, nya, nz)
         !! Field; same stagger and order as at pack.
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*`.
      logical, intent(in) :: negate
         !! `.true.` for a true-vector component.
      logical, intent(in), optional :: device_resident
         !! `.false.` => host arrays; default device-resident.
      integer, intent(in), optional :: ctx
         !! Context (default `FOLD_CTX_NORMAL`).

      integer :: ic, fam, nrow, dbase, g0, ng_e, cap, oT, oU
      logical :: on_device, has_row
      real(wp) :: sgn

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (.not. fxc(ic)%active) return
      if (.not. fxc(ic)%open .or. .not. fxc(ic)%sent) then
         error stop "rdb_ocean_fold_exchange: ocean_fold_unpack before ocean_fold_exchange"
      end if
      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      fam = fold_stagger_family(stagger)
      nrow = fold_stagger_nrows(stagger, fxc(ic)%ng)
      has_row = (nrow == fxc(ic)%ng + 1)
      ! Destination row of message row r: ng+nyl+r for every stagger
      ! (`fold_row_map`: T/u ng+nyl+d, v/corner ng+nyl+1+d with d = r-1).
      dbase = fxc(ic)%ng + fxc(ic)%nyl
      g0 = fxc(ic)%r0(fam)
      ng_e = fxc(ic)%nr(fam)
      cap = fxc(ic)%cap
      oT = fxc(ic)%uoff(FOLD_FAM_T)
      oU = fxc(ic)%uoff(FOLD_FAM_U)
      sgn = 1.0_wp
      if (negate) sgn = -1.0_wp

      if (ng_e > 0) then
         if (on_device) then
            call unpack_dev(fld, nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam, &
                            fxc(ic)%npeer, size(fxc(ic)%rcol), size(fxc(ic)%rbuf), &
                            sgn, negate, has_row, &
                            fxc(ic)%rcol, fxc(ic)%rpeer, fxc(ic)%re, fxc(ic)%rcls, &
                            fxc(ic)%rn, fxc(ic)%rbuf)
         else
            call unpack_host(fld, nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam, &
                             fxc(ic)%npeer, size(fxc(ic)%rcol), size(fxc(ic)%rbuf), &
                             sgn, negate, has_row, &
                             fxc(ic)%rcol, fxc(ic)%rpeer, fxc(ic)%re, fxc(ic)%rcls, &
                             fxc(ic)%rn, fxc(ic)%rbuf)
         end if
      end if
      fxc(ic)%uoff(fam) = fxc(ic)%uoff(fam) + nrow*nz
   end subroutine ocean_fold_unpack_3d

   subroutine unpack_dev(fld, nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam, &
                         npeer, nrecv, nbuf, sgn, negate, has_row, &
                         rcol, rpeer, re, rcls, rn, rbuf)
      !! Device unpack kernel (explicit-shape dummies).
      integer, intent(in) :: nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam
      integer, intent(in) :: npeer, nrecv, nbuf
      real(wp), intent(in) :: sgn
      logical, intent(in) :: negate, has_row
      real(wp), intent(inout) :: fld(nxa, nya, nz)
      integer, intent(in) :: rcol(nrecv), rpeer(nrecv), re(nrecv), rcls(nrecv)
      integer, intent(in) :: rn(npeer, FOLD_NFAM)
      real(wp), intent(in) :: rbuf(nbuf)
      integer :: g, L, r, p, ib, cls
      real(wp) :: val

      do concurrent(L=1:nz, r=1:nrow, g=g0 + 1:g0 + ng_e) local(p, ib, cls, val)
         p = rpeer(g)
         ib = (p - 1)*cap + oT*rn(p, FOLD_FAM_T) + oU*rn(p, FOLD_FAM_U) &
              + ((L - 1)*nrow + (r - 1))*rn(p, fam) + re(g)
         val = sgn*rbuf(ib)
         if (has_row .and. r == 1) then
            cls = rcls(g)
            if (cls == FOLD_ROW_WEST) then
               fld(rcol(g), dbase + 1, L) = val
            else if (cls == FOLD_ROW_SELF .and. negate) then
               fld(rcol(g), dbase + 1, L) = 0.0_wp
            end if
         else
            fld(rcol(g), dbase + r, L) = val
         end if
      end do
   end subroutine unpack_dev

   subroutine unpack_host(fld, nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam, &
                          npeer, nrecv, nbuf, sgn, negate, has_row, &
                          rcol, rpeer, re, rcls, rn, rbuf)
      !! Host twin of `unpack_dev` (setup-time host arrays).
      integer, intent(in) :: nxa, nya, nz, nrow, dbase, g0, ng_e, cap, oT, oU, fam
      integer, intent(in) :: npeer, nrecv, nbuf
      real(wp), intent(in) :: sgn
      logical, intent(in) :: negate, has_row
      real(wp), intent(inout) :: fld(nxa, nya, nz)
      integer, intent(in) :: rcol(nrecv), rpeer(nrecv), re(nrecv), rcls(nrecv)
      integer, intent(in) :: rn(npeer, FOLD_NFAM)
      real(wp), intent(in) :: rbuf(nbuf)
      integer :: g, L, r, p, ib, cls
      real(wp) :: val

      do L = 1, nz
         do r = 1, nrow
            do g = g0 + 1, g0 + ng_e
               p = rpeer(g)
               ib = (p - 1)*cap + oT*rn(p, FOLD_FAM_T) + oU*rn(p, FOLD_FAM_U) &
                    + ((L - 1)*nrow + (r - 1))*rn(p, fam) + re(g)
               val = sgn*rbuf(ib)
               if (has_row .and. r == 1) then
                  cls = rcls(g)
                  if (cls == FOLD_ROW_WEST) then
                     fld(rcol(g), dbase + 1, L) = val
                  else if (cls == FOLD_ROW_SELF .and. negate) then
                     fld(rcol(g), dbase + 1, L) = 0.0_wp
                  end if
               else
                  fld(rcol(g), dbase + r, L) = val
               end if
            end do
         end do
      end do
   end subroutine unpack_host

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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_T, fxc(FOLD_CTX_NORMAL)%ng))
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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_T, fxc(FOLD_CTX_NORMAL)%ng)*nz)
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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_U, fxc(FOLD_CTX_NORMAL)%ng))
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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_U, fxc(FOLD_CTX_NORMAL)%ng)*nz)
         call ocean_fold_pack(fld, nxf, nyt, nz, FOLD_STAG_U, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxf, nyt, nz, FOLD_STAG_U, .true., device_resident)
         call ocean_fold_end()
      else
         call fold_north_u_face(fld, nxf, nyt, nz, nx_phys, ny_phys, ng)
      end if
   end subroutine fold_u_3d

   subroutine fold_v_2d(fld, nxt, nyf, nx_phys, ny_phys, ng, device_resident, ctx)
      !! y-face (v) 2D north fold (negate + fold-line projection).
      integer, intent(in) :: nxt, nyf, nx_phys, ny_phys, ng
      real(wp), intent(inout) :: fld(nxt, nyf)
         !! y-face field (nx_total, ny_total+1).
      logical, intent(in), optional :: device_resident
         !! Exchange path only: `.false.` => host arrays.
      integer, intent(in), optional :: ctx
         !! Exchange path only: context of the plan (default normal; the
         !! wide-halo march-in passes `FOLD_CTX_WIDE`, whose plan was built
         !! for `ng`).

      integer :: ic

      ic = FOLD_CTX_NORMAL
      if (present(ctx)) ic = ctx
      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_V, fxc(ic)%ng), ic)
         call ocean_fold_pack(fld, nxt, nyf, FOLD_STAG_V, device_resident, ic)
         call ocean_fold_exchange(device_resident, ic)
         call ocean_fold_unpack(fld, nxt, nyf, FOLD_STAG_V, .true., device_resident, ic)
         call ocean_fold_end(ic)
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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_V, fxc(FOLD_CTX_NORMAL)%ng)*nz)
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
         call ocean_fold_begin(fold_stagger_nrows(FOLD_STAG_CORNER, fxc(FOLD_CTX_NORMAL)%ng))
         call ocean_fold_pack(fld, nxf, nyf, FOLD_STAG_CORNER, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxf, nyf, FOLD_STAG_CORNER, negate, device_resident)
         call ocean_fold_end()
      else
         call fold_north_corner(fld, nxf, nyf, nx_phys, ny_phys, ng, negate)
      end if
   end subroutine fold_corner_2d

end module rdb_ocean_fold_exchange
