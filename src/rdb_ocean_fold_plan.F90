!! Routing plan for the DISTRIBUTED tripolar north fold (`px > 1`).
module rdb_ocean_fold_plan
   !! Pure index math for the owner-routed north-fold exchange — no MPI,
   !! no field data.  Given the global fold-row width `ni`, the x process
   !! count `px`, the ghost width `ng` and a north-row tile `rx`, it lists,
   !! per peer tile and per column FAMILY, which local storage columns this
   !! tile sends and which it receives.  The exchange engine
   !! (`rdb_ocean_fold_exchange`) moves the values; this module only says
   !! where they come from and where they go, so the routing is testable
   !! serially (`tests/test_ocean_fold_plan.F90`).
   !!
   !! ## Mirror (see the `rdb_ocean_fold` header for the derivation)
   !!
   !! Global physical columns, reduced periodically into `1..ni`:
   !!   * T and v (cell columns, family `FOLD_FAM_T`): column `c` mirrors
   !!     to `ni+1-c`.
   !!   * u and corner (WEST-face / SW-vertex columns, family
   !!     `FOLD_FAM_U`): face `f` (`f = 1 ≡ ni+1`) mirrors to `ni+2-f`.
   !! Rows: the receiver's north ghost row `d` (and, for v / corner, the
   !! fold-line row itself) reads the sender's row `d` below the fold line;
   !! both are north-row tiles of the same height, so the row map needs no
   !! global j offset (`fold_row_map`).
   !!
   !! ## Owner routing (plan §2.3)
   !!
   !! Every value comes from the tile that OWNS the mirror point, never
   !! from a halo copy, so the fold needs no preceding x exchange:
   !!   * T / v column `m`: the tile whose cells contain `m`.
   !!   * u / corner face `mf`: the EAST face of cell `mf-1`, which under
   !!     the halo's D1 rule (the west/south rank owns a seam face) belongs
   !!     to the tile holding cell `mf-1`.  The sender reads its owned faces
   !!     `ng+2 .. ng+w+1`, never the west-seam copy at `ng+1`.
   !!
   !! ## Lists
   !!
   !! A receiver's destination columns are EVERY storage column of its
   !! window: `1..w+2ng` (family T) or `1..w+2ng+1` (family U).  The
   !! canonical entry order of a (sender, receiver) pair is the receiver's
   !! destination columns in increasing order, filtered to that sender;
   !! both ends evaluate the same pure function, so their lists agree with
   !! no handshake.  Lists are index lists, not ranges: a peer's columns
   !! can wrap modulo `ni`.
   !!
   !! ## Fold-line row (v and corner)
   !!
   !! Each receive entry carries the fold-row class of its destination
   !! column (the serial projection rule of `fold_north_v_face_*` /
   !! `fold_north_corner_2d`): `FOLD_ROW_WEST` takes the (sign-applied)
   !! mirror value, `FOLD_ROW_SELF` is zeroed for a true vector and left
   !! for a scalar, `FOLD_ROW_EAST` is the authoritative half and is never
   !! written.  The v class uses the T-family columns, the corner class the
   !! U-family columns, so each family carries exactly one class.
   implicit none
   private

   public :: fold_plan_t
   public :: fold_plan_build
   public :: fold_tile_extent
   public :: fold_tile_owner
   public :: fold_receiver_entries
   public :: fold_stagger_family
   public :: fold_stagger_nrows
   public :: fold_row_map

   integer, parameter, public :: FOLD_STAG_T = 1
      !! Cell centre (h, η, tracers): T columns, `ng` halo rows.
   integer, parameter, public :: FOLD_STAG_U = 2
      !! West face (u): U columns, `ng` halo rows.
   integer, parameter, public :: FOLD_STAG_V = 3
      !! South face (v): T columns, `ng` halo rows + the fold-line row.
   integer, parameter, public :: FOLD_STAG_CORNER = 4
      !! SW vertex: U columns, `ng` halo rows + the fold-line row.

   integer, parameter, public :: FOLD_FAM_T = 1
      !! Column family of T and v (mirror `ni+1-c`).
   integer, parameter, public :: FOLD_FAM_U = 2
      !! Column family of u and corner (mirror `ni+2-f`).
   integer, parameter, public :: FOLD_NFAM = 2
      !! Number of column families.

   integer, parameter, public :: FOLD_ROW_EAST = -1
      !! Fold-row slot in the authoritative (east) half: never written.
   integer, parameter, public :: FOLD_ROW_SELF = 0
      !! Self-conjugate fold-row slot: zeroed for vectors, left for scalars.
   integer, parameter, public :: FOLD_ROW_WEST = 1
      !! Fold-row slot in the west half: takes the sign-applied mirror.

   integer, parameter, public :: FOLD_PLAN_OK = 0
      !! `fold_plan_build` status: plan built.
   integer, parameter, public :: FOLD_PLAN_ERR_ARGS = 1
      !! `fold_plan_build` status: `ni < px`, `px < 1`, `ng < 1` or `rx`
      !! outside `0..px-1`.

   type :: fold_plan_t
      !! One north-row tile's routing for the distributed fold.  Entries of
      !! family `f` for peer `p` live at `start(p,f)+1 .. start(p,f)+n(p,f)`
      !! of the flat entry arrays (send and receive separately).
      integer :: ni = 0
         !! Global physical width of the fold row.
      integer :: px = 0
         !! Tiles along the fold row.
      integer :: ng = 0
         !! Ghost width.
      integer :: rx = -1
         !! This tile's x-coordinate (0-based).
      integer :: npeer = 0
         !! Peer tiles this tile sends to or receives from (self included
         !! when its own mirror falls in its window).
      integer :: self_peer = 0
         !! Index into `peer_rx` of this tile itself (0 if not a peer).
      integer :: nmax = 0
         !! Largest entry count of any (peer, family) list, send or receive.
      integer, allocatable :: peer_rx(:)
         !! Peer tile x-coordinates, ascending, shape (npeer).
      integer, allocatable :: send_n(:, :)
         !! Send entry counts, shape (npeer, FOLD_NFAM).
      integer, allocatable :: send_start(:, :)
         !! Offsets into the flat send arrays, shape (npeer, FOLD_NFAM).
      integer, allocatable :: recv_n(:, :)
         !! Receive entry counts, shape (npeer, FOLD_NFAM).
      integer, allocatable :: recv_start(:, :)
         !! Offsets into the flat receive arrays, shape (npeer, FOLD_NFAM).
      integer :: nsend(FOLD_NFAM) = 0
         !! Total send entries per family.
      integer :: nrecv(FOLD_NFAM) = 0
         !! Total receive entries per family.
      integer, allocatable :: send_col(:)
         !! Local storage column this tile reads, flat (all families).
      integer, allocatable :: send_peer(:)
         !! Peer index of each send entry, flat.
      integer, allocatable :: send_e(:)
         !! 1-based position of each send entry within its (peer, family)
         !! list, flat.
      integer, allocatable :: recv_col(:)
         !! Local storage column this tile writes, flat (all families).
      integer, allocatable :: recv_peer(:)
         !! Peer index of each receive entry, flat.
      integer, allocatable :: recv_e(:)
         !! 1-based position within its (peer, family) list, flat.
      integer, allocatable :: recv_cls(:)
         !! Fold-row class of each receive entry (`FOLD_ROW_*`), flat.
   contains
      procedure :: destroy => fold_plan_destroy
   end type fold_plan_t

contains

   pure subroutine fold_tile_extent(ni, px, rx, a, w)
      !! First global cell `a` and width `w` of tile `rx` — the
      !! `decomp_init` split (remainder to the WEST tiles), restated here so
      !! the plan stays pure and dependency-free (cross-checked against
      !! `decomp_init` by the unit test).
      integer, intent(in) :: ni
         !! Global physical width of the fold row (cells).
      integer, intent(in) :: px
         !! Tiles along the fold row.
      integer, intent(in) :: rx
         !! Tile x-coordinate (0-based, `0..px-1`).
      integer, intent(out) :: a
         !! First global cell of the tile (1-based).
      integer, intent(out) :: w
         !! Tile width (cells).

      integer :: base, rem

      base = ni/px
      rem = mod(ni, px)
      if (rx < rem) then
         w = base + 1
         a = rx*(base + 1) + 1
      else
         w = base
         a = rem*(base + 1) + (rx - rem)*base + 1
      end if
   end subroutine fold_tile_extent

   pure function fold_tile_owner(ni, px, c) result(rx)
      !! Tile holding global cell `c` (1..ni), closed form of the
      !! `decomp_init` split.
      integer, intent(in) :: ni
         !! Global physical width of the fold row (cells).
      integer, intent(in) :: px
         !! Tiles along the fold row.
      integer, intent(in) :: c
         !! Global cell index (1-based, `1..ni`).
      integer :: rx
         !! Owning tile x-coordinate (0-based).

      integer :: base, rem

      base = ni/px
      rem = mod(ni, px)
      if (c <= rem*(base + 1)) then
         rx = (c - 1)/(base + 1)
      else
         rx = rem + (c - 1 - rem*(base + 1))/base
      end if
   end function fold_tile_owner

   pure integer function fold_stagger_family(stagger) result(fam)
      !! Column family of a stagger (`FOLD_FAM_T` for T/v, `FOLD_FAM_U` for
      !! u/corner).
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*` stagger.

      if (stagger == FOLD_STAG_U .or. stagger == FOLD_STAG_CORNER) then
         fam = FOLD_FAM_U
      else
         fam = FOLD_FAM_T
      end if
   end function fold_stagger_family

   pure integer function fold_stagger_nrows(stagger, ng) result(nrow)
      !! Rows a stagger moves per column: `ng` (T, u) or `ng+1` (v, corner:
      !! the fold-line row is always sent, the receiver uses it only where
      !! `recv_cls == FOLD_ROW_WEST`).
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*` stagger.
      integer, intent(in) :: ng
         !! Ghost width.

      if (stagger == FOLD_STAG_V .or. stagger == FOLD_STAG_CORNER) then
         nrow = ng + 1
      else
         nrow = ng
      end if
   end function fold_stagger_nrows

   pure subroutine fold_row_map(stagger, ng, nyl, r, src_row, dst_row)
      !! Storage rows of message row `r` (1..`fold_stagger_nrows`) on a
      !! north-row tile of `nyl` physical rows: the sender reads
      !! `src_row`, the receiver writes `dst_row`.
      !!   * T, u:  r = d = 1..ng:   src ng+nyl+1-d, dst ng+nyl+d.
      !!   * v, corner: r = d+1, d = 0..ng: src ng+nyl+1-d, dst ng+nyl+1+d
      !!     (d = 0 is the fold-line row, src = dst = ng+nyl+1).
      integer, intent(in) :: stagger
         !! `FOLD_STAG_*` stagger.
      integer, intent(in) :: ng
         !! Ghost width.
      integer, intent(in) :: nyl
         !! Physical rows of the north-row tile.
      integer, intent(in) :: r
         !! Message row (1-based, `1..fold_stagger_nrows(stagger, ng)`).
      integer, intent(out) :: src_row
         !! Local storage row the sender reads.
      integer, intent(out) :: dst_row
         !! Local storage row the receiver writes.

      integer :: d

      if (stagger == FOLD_STAG_V .or. stagger == FOLD_STAG_CORNER) then
         d = r - 1
         src_row = ng + nyl + 1 - d
         dst_row = ng + nyl + 1 + d
      else
         d = r
         src_row = ng + nyl + 1 - d
         dst_row = ng + nyl + d
      end if
   end subroutine fold_row_map

   pure subroutine fold_receiver_entries(ni, px, ng, fam, rx, ncol, own, src, cls)
      !! Every destination column `i = 1..ncol` of tile `rx`'s window for a
      !! column family: the owning tile `own(i)` of its mirror, the owner's
      !! local storage column `src(i)` to read, and the fold-row class
      !! `cls(i)`.  `ncol` = `w+2ng` (T) or `w+2ng+1` (U); the arrays must
      !! hold at least that many entries.
      integer, intent(in) :: ni
         !! Global physical width of the fold row (cells).
      integer, intent(in) :: px
         !! Tiles along the fold row.
      integer, intent(in) :: ng
         !! Ghost width.
      integer, intent(in) :: fam
         !! Column family (`FOLD_FAM_T` or `FOLD_FAM_U`).
      integer, intent(in) :: rx
         !! Receiving tile x-coordinate (0-based).
      integer, intent(out) :: ncol
         !! Destination columns filled (`w+2ng` for T, `w+2ng+1` for U).
      integer, intent(out) :: own(:)
         !! Owning tile x-coordinate (0-based) of each column's mirror.
      integer, intent(out) :: src(:)
         !! Owner's local storage column to read (1-based, ghosts included).
      integer, intent(out) :: cls(:)
         !! Fold-row class of each column (`FOLD_ROW_*`).

      integer :: a, w, a_o, w_o, i, g, c, m, p, pm

      call fold_tile_extent(ni, px, rx, a, w)
      if (fam == FOLD_FAM_U) then
         ncol = w + 2*ng + 1
      else
         ncol = w + 2*ng
      end if

      do i = 1, ncol
         ! Global column of storage slot i, reduced into 1..ni.
         g = a + i - ng - 1
         c = modulo(g - 1, ni) + 1
         if (fam == FOLD_FAM_U) then
            ! Face c (west face of cell c, c = 1 ≡ ni+1) mirrors to face
            ! m = ni+2-c in 2..ni+1: the EAST face of cell m-1, owned (D1)
            ! by the tile holding that cell.
            m = ni + 2 - c
            own(i) = fold_tile_owner(ni, px, m - 1)
            call fold_tile_extent(ni, px, own(i), a_o, w_o)
            src(i) = ng + m - a_o + 1
            ! Corner projection rule (`fold_north_corner_2d`).
            p = c
            pm = modulo(ni + 1 - p, ni) + 1
         else
            m = ni + 1 - c
            own(i) = fold_tile_owner(ni, px, m)
            call fold_tile_extent(ni, px, own(i), a_o, w_o)
            src(i) = ng + m - a_o + 1
            ! v projection rule (`fold_north_v_face_*`).
            p = c
            pm = ni + 1 - p
         end if
         if (p < pm) then
            cls(i) = FOLD_ROW_WEST
         else if (p == pm) then
            cls(i) = FOLD_ROW_SELF
         else
            cls(i) = FOLD_ROW_EAST
         end if
      end do
   end subroutine fold_receiver_entries

   pure subroutine fold_plan_build(plan, ni, px, ng, rx, status)
      !! Build tile `rx`'s send/receive lists for every peer and both
      !! column families.  Pure: every north-row rank computes every tile's
      !! receive lists itself (`decomp_init` arithmetic), so no handshake.
      type(fold_plan_t), intent(out) :: plan
      integer, intent(in) :: ni
         !! Global physical width of the fold row.
      integer, intent(in) :: px
         !! Tiles along the fold row.
      integer, intent(in) :: ng
         !! Ghost width.
      integer, intent(in) :: rx
         !! This tile's x-coordinate (0-based).
      integer, intent(out) :: status
         !! `FOLD_PLAN_OK`, or `FOLD_PLAN_ERR_ARGS`.

      integer, allocatable :: own(:), src(:), cls(:)
      integer, allocatable :: sn(:, :), rn(:, :), pidx(:)
      integer :: ncap, ncol, fam, r, i, p, k, ks, kr
      integer :: wmax, a, w

      status = FOLD_PLAN_OK
      if (px < 1 .or. ng < 1 .or. ni < px .or. rx < 0 .or. rx >= px) then
         status = FOLD_PLAN_ERR_ARGS
         return
      end if

      plan%ni = ni
      plan%px = px
      plan%ng = ng
      plan%rx = rx

      call fold_tile_extent(ni, px, 0, a, wmax)   ! tile 0 is the widest
      ncap = wmax + 2*ng + 1
      allocate (own(ncap), src(ncap), cls(ncap))
      ! Per-tile counts, indexed by tile x-coordinate 0..px-1.
      allocate (sn(0:px - 1, FOLD_NFAM), rn(0:px - 1, FOLD_NFAM))
      sn = 0
      rn = 0

      ! Pass 1: counts.  Receive: my own window, by owner.  Send: every
      ! receiver's window, the entries whose owner is me.
      do fam = 1, FOLD_NFAM
         do r = 0, px - 1
            call fold_receiver_entries(ni, px, ng, fam, r, ncol, own, src, cls)
            do i = 1, ncol
               if (r == rx) rn(own(i), fam) = rn(own(i), fam) + 1
               if (own(i) == rx) sn(r, fam) = sn(r, fam) + 1
            end do
         end do
      end do

      ! Peers: ascending tile order, union of send and receive partners.
      allocate (pidx(0:px - 1))
      pidx = 0
      plan%npeer = 0
      do r = 0, px - 1
         if (sum(sn(r, :)) + sum(rn(r, :)) > 0) then
            plan%npeer = plan%npeer + 1
            pidx(r) = plan%npeer
         end if
      end do
      allocate (plan%peer_rx(plan%npeer))
      allocate (plan%send_n(plan%npeer, FOLD_NFAM), plan%send_start(plan%npeer, FOLD_NFAM))
      allocate (plan%recv_n(plan%npeer, FOLD_NFAM), plan%recv_start(plan%npeer, FOLD_NFAM))
      do r = 0, px - 1
         if (pidx(r) > 0) then
            plan%peer_rx(pidx(r)) = r
            plan%send_n(pidx(r), :) = sn(r, :)
            plan%recv_n(pidx(r), :) = rn(r, :)
         end if
      end do
      plan%self_peer = pidx(rx)

      ! Offsets: family-major, then peer, into one flat array per side.
      ks = 0
      kr = 0
      do fam = 1, FOLD_NFAM
         do p = 1, plan%npeer
            plan%send_start(p, fam) = ks
            plan%recv_start(p, fam) = kr
            ks = ks + plan%send_n(p, fam)
            kr = kr + plan%recv_n(p, fam)
         end do
         plan%nsend(fam) = sum(plan%send_n(:, fam))
         plan%nrecv(fam) = sum(plan%recv_n(:, fam))
      end do
      plan%nmax = 0
      if (plan%npeer > 0) plan%nmax = max(maxval(plan%send_n), maxval(plan%recv_n))
      allocate (plan%send_col(max(ks, 1)), plan%send_peer(max(ks, 1)), plan%send_e(max(ks, 1)))
      allocate (plan%recv_col(max(kr, 1)), plan%recv_peer(max(kr, 1)), &
                plan%recv_e(max(kr, 1)), plan%recv_cls(max(kr, 1)))

      ! Pass 2: fill, in the canonical order (receiver's destination
      ! columns ascending).  `sn`/`rn` are reused as running cursors.
      sn = 0
      rn = 0
      do fam = 1, FOLD_NFAM
         do r = 0, px - 1
            call fold_receiver_entries(ni, px, ng, fam, r, ncol, own, src, cls)
            do i = 1, ncol
               if (r == rx) then
                  p = pidx(own(i))
                  rn(own(i), fam) = rn(own(i), fam) + 1
                  k = plan%recv_start(p, fam) + rn(own(i), fam)
                  plan%recv_col(k) = i
                  plan%recv_peer(k) = p
                  plan%recv_e(k) = rn(own(i), fam)
                  plan%recv_cls(k) = cls(i)
               end if
               if (own(i) == rx) then
                  p = pidx(r)
                  sn(r, fam) = sn(r, fam) + 1
                  k = plan%send_start(p, fam) + sn(r, fam)
                  plan%send_col(k) = src(i)
                  plan%send_peer(k) = p
                  plan%send_e(k) = sn(r, fam)
               end if
            end do
         end do
      end do
   end subroutine fold_plan_build

   pure subroutine fold_plan_destroy(this)
      !! Free the plan's lists and reset it to the empty default.  Safe on
      !! a plan that was never built.
      class(fold_plan_t), intent(inout) :: this
      if (allocated(this%peer_rx)) deallocate (this%peer_rx)
      if (allocated(this%send_n)) deallocate (this%send_n)
      if (allocated(this%send_start)) deallocate (this%send_start)
      if (allocated(this%recv_n)) deallocate (this%recv_n)
      if (allocated(this%recv_start)) deallocate (this%recv_start)
      if (allocated(this%send_col)) deallocate (this%send_col)
      if (allocated(this%send_peer)) deallocate (this%send_peer)
      if (allocated(this%send_e)) deallocate (this%send_e)
      if (allocated(this%recv_col)) deallocate (this%recv_col)
      if (allocated(this%recv_peer)) deallocate (this%recv_peer)
      if (allocated(this%recv_e)) deallocate (this%recv_e)
      if (allocated(this%recv_cls)) deallocate (this%recv_cls)
      this%ni = 0
      this%px = 0
      this%ng = 0
      this%rx = -1
      this%npeer = 0
      this%self_peer = 0
      this%nmax = 0
      this%nsend = 0
      this%nrecv = 0
   end subroutine fold_plan_destroy

end module rdb_ocean_fold_plan
