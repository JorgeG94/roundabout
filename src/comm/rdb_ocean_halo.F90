!! Ocean-grid staggered halo-exchange — MPI backend.
!!
!! Implements the phase-O1 ocean halo API for multi-rank builds.
!!
!! ARCHITECTURE
!! ============
!! Two-pass X-then-Y ordering: every exchange runs the FULL X pass (E/W
!! messages or local wrap) FIRST, then the FULL Y pass.  The Y-pass pack
!! ranges span the FULL i extent INCLUDING the ghost columns the X pass
!! just filled.  This makes diagonal corner ghosts correct by transitivity
!! with no corner messages.
!!
!! INVARIANT:
!!   Corner ghosts are valid after a full two-pass exchange.
!!   Never call a single-axis partial exchange before a corner-reading kernel.
!!
!! FACE OWNERSHIP (D1)
!! ===================
!! For interior seams: the WEST/SOUTH rank OWNS the shared seam face.
!!
!!   Face-x: rank R owns the seam face at its local column i_own = ng + nxl + 1.
!!     R's east neighbour (rank E) holds a copy at ITS local i_copy = ng + 1.
!!     R overwrites E's copy during every exchange by sending nghost+1 values east.
!!     E does NOT send its i_copy = ng+1 back westward to R.
!!
!!   Face-y: symmetric, south rank owns seam face at j_own = ng + nyl + 1.
!!
!! For periodic wrap links: ownership is UNIFORMLY D1 — the western partner
!!   of every exchanging pair owns the shared seam face.  For the x-wrap pair
!!   (rx = px-1 → rx = 0) the east-most rank IS the western partner, so it
!!   owns the wrap-seam face: it sends nghost+1 values "eastward" (across the
!!   wrap) to rx=0, filling rx=0's west ghosts i=1..ng AND its west-seam copy
!!   i=ng+1; rx=0 sends nghost values back, filling the east-most rank's east
!!   ghost columns only.  NOTE this is the REVERSE of the single-rank wrap
!!   kernel's belt-and-braces copy direction (fld(i_e) := fld(i_w), i.e. the
!!   west-edge copy authoritative) — harmless, because the two copies of the
!!   wrap face are bit-equal whenever the seam invariant holds, and the
!!   single-rank (px==1) path delegates to the wrap kernels verbatim, so
!!   single-rank behaviour is unchanged.
!!
!! FACE-X CENTRE INDEX MAP (equal-split, nghost=ng, nx_local=nxl)
!! ===============================================================
!!   array dims: (nxl + 2*ng + 1) x ny_total  (= nx_total+1 = nxt+1)
!!   physical face cols: i = ng+1 .. ng+nxl+1
!!     i = ng+1:     west-seam copy (owned by WEST neighbour)
!!     i = ng+nxl+1: east seam face (owned by THIS rank)
!!
!!   Send EAST  [ng+1 values per j-row]:
!!     source i = nxl+1 .. ng+nxl+1  (last element is owned seam face)
!!     fills east-nbr's i = 1..ng (west ghosts) + i = ng+1 (seam copy)
!!     (same-global-face identity: my i == east-nbr's i_E + nxl)
!!
!!   Recv FROM EAST  [ng values per j-row]:
!!     east-nbr sends i = ng+2..2*ng+1 → fills OUR i = ng+nxl+2..ng+nxl+ng+1 (east ghosts)
!!
!!   Send WEST  [ng values per j-row]:
!!     source i = ng+2..2*ng+1  (first ng interior faces after west-seam)
!!     fills west-nbr's east ghosts i = ng+nxl+2..ng+nxl+ng+1
!!
!!   Recv FROM WEST  [ng+1 values per j-row]:
!!     west-nbr sends its i = nxl+1..ng+nxl+1 (ng pre-seam faces + owned seam)
!!     fills OUR i = 1..ng+1 (west ghosts + seam copy)
!!
!! FACE-Y CENTRE INDEX MAP: symmetric (transpose i↔j, x↔y).
!!
!! TAG SCHEME (non-overlapping with rdb_halo tags 1-4)
!! ===========================================
!!   TAG_OC_W_TO_E = 11  (rank's eastward send; peer receives as "from west")
!!   TAG_OC_E_TO_W = 12  (rank's westward send; peer receives as "from east")
!!   TAG_OC_S_TO_N = 13  (rank's northward send)
!!   TAG_OC_N_TO_S = 14  (rank's southward send)
!!
!!   2-rank periodic disambiguation (px=2):
!!     rank 0 and rank 1 are each other's east AND west periodic wrap partners.
!!     Each rank posts ONE isend (tag 11 east OR tag 12 west) and ONE irecv
!!     (the complementary tag from the same partner).  With px=2 each rank
!!     sends both eastward (tag 11) AND westward (tag 12) to rank 1-rx.
!!     Four messages per rank, four matching irecvs — no ambiguity because
!!     each (source_rank, tag) pair is posted exactly once.
!!
!! v1: single code path — device-resident under OpenACC with CUDA-aware MPI;
!!     host-staged otherwise.  O2 staging variant deferred to solver wiring.
!!
!! O3 (DONE): centre_3d / face_*_3d are now batched — one message per
!!     direction carrying all nz layers.  Persistent buffers grow on first
!!     call with nz>1 via ocean_halo_buffers_ensure_nz (exit-delete /
!!     dealloc / alloc / enter-create pattern; no per-call allocation).
!!     The 2D primitives use the same buffers unchanged (their 2D capacity
!!     is always <= any nz>=1 grown capacity).
#ifdef RDB_DOUBLE_PRECISION
#define HALO_ISEND_N comm_isend_real_dp_array_n
#define HALO_IRECV_N comm_irecv_real_dp_array_n
#else
#define HALO_ISEND_N comm_isend_real_sp_array_n
#define HALO_IRECV_N comm_irecv_real_sp_array_n
#endif
module rdb_ocean_halo
   !! Ocean staggered halo exchange — the ONLY backend.
   !!
   !! At px==1 / py==1 the X / Y pass never enters a collective: it either
   !! applies `ocean_periodic_wrap_*` locally or does nothing.  That is
   !! what lets this module compile and run against pic-mpi's serial
   !! backend, which `error stop`s on point-to-point calls by design.
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_rank_from_coords
   use pic_mpi_lib, only: comm_t, request_t, MPI_Status, &
                          isend, irecv, waitall, &
                          HALO_ISEND_N, HALO_IRECV_N
   use rdb_comm_env, only: comm_env_compute_comm
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_centre_3d, &
                                 ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_x_3d, &
                                 ocean_periodic_wrap_face_y_2d, &
                                 ocean_periodic_wrap_face_y_3d
   use rdb_ocean_halo_counters, only: oh_count_centre_2d, oh_count_centre_3d, &
                                      oh_count_face_x_2d, oh_count_face_x_3d, &
                                      oh_count_face_y_2d, oh_count_face_y_3d, &
                                      oh_count_bt_group, oh_count_msgs, &
                                      oh_count_suppress_on, oh_count_suppress_off
   implicit none
   private

   public :: ocean_halo_init
   public :: ocean_halo_destroy
   public :: ocean_halo_reserve
   public :: ocean_halo_is_init
   public :: ocean_halo_is_decomposed
   public :: ocean_halo_is_decomposed_x
   public :: ocean_halo_is_decomposed_y
   public :: ocean_halo_centre
   public :: ocean_halo_face_x
   public :: ocean_halo_face_y
   public :: ocean_halo_bt_group_2d

   ! Rank-generic dispatch: ocean_halo_<stagger>(fld[, nz][, device_resident])
   ! selects the 2D or 3D specific by the rank of fld (2 vs 3).  The
   ! specifics stay private -- no caller outside this module names them.
   interface ocean_halo_centre
      module procedure ocean_halo_centre_2d, ocean_halo_centre_3d
   end interface ocean_halo_centre

   interface ocean_halo_face_x
      module procedure ocean_halo_face_x_2d, ocean_halo_face_x_3d
   end interface ocean_halo_face_x

   interface ocean_halo_face_y
      module procedure ocean_halo_face_y_2d, ocean_halo_face_y_3d
   end interface ocean_halo_face_y
   public :: ocean_halo_centre_2d_wide
   public :: ocean_halo_face_x_2d_wide
   public :: ocean_halo_face_y_2d_wide
   public :: ocean_halo_bt_group_2d_wide

   ! -----------------------------------------------------------------
   ! Tags (non-overlapping with coastal halo tags 1-4)
   ! -----------------------------------------------------------------
   integer, parameter :: TAG_OC_W_TO_E = 11
      !! Eastward send (filled into west-recv on peer)
   integer, parameter :: TAG_OC_E_TO_W = 12
      !! Westward send (filled into east-recv on peer)
   integer, parameter :: TAG_OC_S_TO_N = 13
      !! Northward send
   integer, parameter :: TAG_OC_N_TO_S = 14
      !! Southward send

   integer, parameter :: MAX_REQS = 8
      !! Upper bound on requests per exchange (2 per active direction x 2 axes)

   ! -----------------------------------------------------------------
   ! Module-level topology
   ! -----------------------------------------------------------------
   integer :: oh_nghost = 0
      !! Ghost cell width
   integer :: oh_nx_local = 0
      !! Physical cells in x on this subdomain
   integer :: oh_ny_local = 0
      !! Physical cells in y on this subdomain
   integer :: oh_nx_total = 0
      !! nx_local + 2*nghost
   integer :: oh_ny_total = 0
      !! ny_local + 2*nghost
   logical :: oh_periodic_x = .false.
      !! x axis is reentrant periodic
   logical :: oh_periodic_y = .false.
      !! y axis is reentrant periodic
   type(decomp_t) :: oh_decomp
      !! Copy of domain decomposition
   logical :: oh_initialised = .false.
      !! True after ocean_halo_init

   ! -----------------------------------------------------------------
   ! Persistent device-resident send/recv buffers.
   !
   ! Buffer naming convention (from THIS rank's perspective):
   !   send_east : data THIS rank sends to its EAST neighbour
   !   recv_east : data THIS rank receives from its EAST neighbour
   !   send_west : data THIS rank sends to its WEST neighbour
   !   recv_west : data THIS rank receives from its WEST neighbour
   !   (and symmetric for N/S)
   !
   ! Buffer capacity:
   !   E/W east direction (send east, recv from east):
   !     For centre and face-y: nghost*ny_total per j-row
   !     For face-x (send east): (nghost+1)*ny_total — THIS RANK sends seam+nghost
   !     For face-x (recv from east): nghost*ny_total — EAST sends only nghost west
   !     So cap_ew_east = max((ng+1)*ny_total, ng*ny_total) = (ng+1)*ny_total
   !
   !   E/W west direction (send west, recv from west):
   !     For centre and face-y: nghost*ny_total
   !     For face-x (send west): nghost*ny_total
   !     For face-x (recv from west): (nghost+1)*ny_total — WEST sends seam+nghost
   !     So cap_ew_west = (ng+1)*ny_total (for recv_west face-x case)
   !
   !   N/S north direction: (ng+1)*nx_total for face-y send-north case
   !   N/S south direction: (ng+1)*nx_total for face-y recv-from-north case
   !
   !   Use same capacity for all four E/W and all four N/S buffers =
   !   max((ng+1)*ny_total, ng*(ny_total+1)) for E/W
   !   max((ng+1)*nx_total, ng*(nx_total+1)) for N/S
   !   In practice (ng+1)*ny_total >= ng*(ny_total+1) iff ny_total >= ng which
   !   holds for any sensible grid, so use (ng+1)*ny_total and (ng+1)*nx_total.
   !   Also need to handle face-y x-pass: uses ng*(ny_total+1) per direction.
   !   Check: (ng+1)*ny_total >= ng*(ny_total+1) = ng*ny_total + ng
   !          ng*ny_total + ny_total >= ng*ny_total + ng  iff  ny_total >= ng  (true).
   ! -----------------------------------------------------------------

   real(wp), allocatable :: oh_buf_send_east(:), oh_buf_recv_east(:)
      !! Data sent/received to/from the east neighbour (capacity = (ng+1)*ny_total)
   real(wp), allocatable :: oh_buf_send_west(:), oh_buf_recv_west(:)
      !! Data sent/received to/from the west neighbour (capacity = (ng+1)*ny_total)
   real(wp), allocatable :: oh_buf_send_north(:), oh_buf_recv_north(:)
      !! Data sent/received to/from the north neighbour (capacity = (ng+1)*nx_total)
   real(wp), allocatable :: oh_buf_send_south(:), oh_buf_recv_south(:)
      !! Data sent/received to/from the south neighbour (capacity = (ng+1)*nx_total)
   integer :: oh_cap_ew = 0
      !! Allocated capacity for E/W buffers
   integer :: oh_cap_ns = 0
      !! Allocated capacity for N/S buffers

contains

   ! ==================================================================
   ! Lifecycle
   ! ==================================================================

   subroutine ocean_halo_init(decomp, nghost, periodic_x, periodic_y, ierr)
      !! Initialise topology and allocate persistent device-resident buffers.
      !! Idempotent: destroys prior state before re-init on grid change.
      type(decomp_t), intent(in) :: decomp
         !! Domain decomposition descriptor
      integer, intent(in) :: nghost
         !! Ghost cell width
      logical, intent(in) :: periodic_x
         !! Reentrant periodic in x
      logical, intent(in) :: periodic_y
         !! Reentrant periodic in y
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`) on the decomposed-run
         !! nghost<3 seam-stencil guard when present; absent behaves as
         !! today (`error stop`). (F5 residual, P2.4)

      integer :: cap_ew, cap_ns

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      ! O4 seam-stencil guard (FAIL-LOUD): the ocean dyn-core's
      ! continuity/tracer PPM reconstruction is 5-point; a seam face's
      ! ghost-side donor cell needs two neighbours beyond itself, so
      ! full-order reconstruction at a rank-seam face needs THREE ghost
      ! columns.  With nghost < 3 the PPM local-array-edge fallback fires AT
      ! the seam and degrades the face to first order, producing a
      ! seam-inconsistent tracer mass flux: it conserves a UNIFORM tracer
      ! (Σ flux·C cancels) but LEAKS a STRUCTURED one — measured heat
      ! closure grew to +5e-11 by day 2 on seamount_bench_full np2 at
      ! nghost=2 (invisible to the mass/salt gates; the O4 heat-leak
      ! postmortem).  This is an ABORT, not a warning: a silent
      ! decomposition-dependent conservation leak is a correctness bug.
      ! `decomp%px/py` are the TRUE runtime decomposition (namelist px/py
      ! are auto-factored by then), so this gate sees the real seam even
      ! when the namelist left px=py=1.  The default nghost is now 3, so
      ! this only trips an EXPLICIT decomposed override below 3.  The ocean
      ! halo-primitive tests (test_halo_ocean_mpi) already init at nghost=3.
      if ((decomp%px > 1 .or. decomp%py > 1) .and. nghost < 3) then
         call fail("ocean_halo_init: decomposed ocean run (px="// &
                   to_string(decomp%px)//", py="//to_string(decomp%py)// &
                   ") requires nghost >= 3, got nghost = "//to_string(nghost)// &
                   ". The continuity/tracer PPM 5-point stencil degrades to "// &
                   "FIRST ORDER at rank-seam faces, leaking a structured tracer "// &
                   "(heat) while conserving a uniform one (salt). "// &
                   "Set &grid_nml nghost = 3 (the default) or higher.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      if (oh_initialised) call ocean_halo_destroy()

      oh_decomp = decomp
      oh_nghost = nghost
      oh_nx_local = decomp%nx_local
      oh_ny_local = decomp%ny_local
      oh_nx_total = decomp%nx_local + 2*nghost
      oh_ny_total = decomp%ny_local + 2*nghost
      oh_periodic_x = periodic_x
      oh_periodic_y = periodic_y
      oh_initialised = .true.

      ! Capacity: worst case is face-x send-east = (nghost+1)*ny_total for E/W.
      ! For face-y x-pass, ng*(ny_total+1) <= (ng+1)*ny_total when ny_total >= ng.
      cap_ew = (nghost + 1)*oh_ny_total
      cap_ns = (nghost + 1)*oh_nx_total

      allocate (oh_buf_send_east(cap_ew), oh_buf_recv_east(cap_ew))
      allocate (oh_buf_send_west(cap_ew), oh_buf_recv_west(cap_ew))
      allocate (oh_buf_send_north(cap_ns), oh_buf_recv_north(cap_ns))
      allocate (oh_buf_send_south(cap_ns), oh_buf_recv_south(cap_ns))
      !$acc enter data create(oh_buf_send_east, oh_buf_recv_east, &
      !$acc&                  oh_buf_send_west, oh_buf_recv_west, &
      !$acc&                  oh_buf_send_north, oh_buf_recv_north, &
      !$acc&                  oh_buf_send_south, oh_buf_recv_south)

      oh_cap_ew = cap_ew
      oh_cap_ns = cap_ns

   end subroutine ocean_halo_init

   subroutine ocean_halo_destroy()
      !! Release buffers and topology state.  Idempotent.
      if (.not. oh_initialised) return

      !$acc exit data delete(oh_buf_send_east, oh_buf_recv_east, &
      !$acc&                 oh_buf_send_west, oh_buf_recv_west, &
      !$acc&                 oh_buf_send_north, oh_buf_recv_north, &
      !$acc&                 oh_buf_send_south, oh_buf_recv_south)
      deallocate (oh_buf_send_east, oh_buf_recv_east)
      deallocate (oh_buf_send_west, oh_buf_recv_west)
      deallocate (oh_buf_send_north, oh_buf_recv_north)
      deallocate (oh_buf_send_south, oh_buf_recv_south)

      oh_cap_ew = 0
      oh_cap_ns = 0
      oh_nghost = 0
      oh_nx_local = 0
      oh_ny_local = 0
      oh_nx_total = 0
      oh_ny_total = 0
      oh_periodic_x = .false.
      oh_periodic_y = .false.
      oh_initialised = .false.

   end subroutine ocean_halo_destroy

   subroutine ocean_halo_reserve(nz, ng_wide, ierr)
      !! Pre-size the pack/recv buffers to the worst-case capacity for this
      !! run at init time, so that NO device reallocation fires mid-run.
      !!
      !! Worst-case capacity derivation:
      !!
      !!   E/W buffers:
      !!     (a) 2D face-x send-east:  (ng+1)*nyt            (standard init)
      !!     (b) batched 3D face-x:    (ng+1)*nyt*nz          (ensure_nz path)
      !!     (c) wide 2D face-x:       (ng_wide+1)*(ny_local+2*ng_wide) (ensure_wide path)
      !!     worst = max(a, b, c)
      !!
      !!   N/S buffers: symmetric (transpose x↔y).
      !!
      !! Call AFTER ocean_halo_init; calling before init is a programming error
      !! (error stop with a clear message).  Passing nz = 0 or ng_wide = 0
      !! means "not used" — those arms are skipped.  If the computed capacity
      !! is already covered by the current allocation, this is a no-op.
      integer, intent(in) :: nz
         !! Layer count for 3D batched exchanges (0 = not used)
      integer, intent(in) :: ng_wide
         !! Wide-ghost width for the barotropic march-in (0 = not used)
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`) when called before
         !! `ocean_halo_init` when present; absent behaves as today
         !! (`error stop`). (F5 residual, P2.4)

      integer :: need_ew, need_ns
      integer :: w_ew, w_ns
      integer :: ng, nyt, nxt, nyl, nxl

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. oh_initialised) then
         call fail("ocean_halo_reserve: called before ocean_halo_init. "// &
                   "Call ocean_halo_reserve immediately after ocean_halo_init.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      ng = oh_nghost
      nyt = oh_ny_total
      nxt = oh_nx_total
      nyl = oh_ny_local
      nxl = oh_nx_local

      ! (a) 2D base capacity (already set by init, replicated here for clarity)
      need_ew = (ng + 1)*nyt
      need_ns = (ng + 1)*nxt

      ! (b) 3D batched capacity
      if (nz > 0) then
         need_ew = max(need_ew, (ng + 1)*nyt*nz)
         need_ns = max(need_ns, (ng + 1)*nxt*nz)
      end if

      ! (c) Wide-ghost 2D capacity
      if (ng_wide > 0) then
         w_ew = (ng_wide + 1)*(nyl + 2*ng_wide)
         w_ns = (ng_wide + 1)*(nxl + 2*ng_wide)
         need_ew = max(need_ew, w_ew)
         need_ns = max(need_ns, w_ns)
      end if

      if (need_ew <= oh_cap_ew .and. need_ns <= oh_cap_ns) return

      ! Grow E/W if needed (exit-delete / dealloc / alloc / enter-create)
      if (need_ew > oh_cap_ew) then
         !$acc exit data delete(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                 oh_buf_send_west, oh_buf_recv_west)
         deallocate (oh_buf_send_east, oh_buf_recv_east)
         deallocate (oh_buf_send_west, oh_buf_recv_west)
         allocate (oh_buf_send_east(need_ew), oh_buf_recv_east(need_ew))
         allocate (oh_buf_send_west(need_ew), oh_buf_recv_west(need_ew))
         !$acc enter data create(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                  oh_buf_send_west, oh_buf_recv_west)
         oh_cap_ew = need_ew
      end if

      ! Grow N/S if needed
      if (need_ns > oh_cap_ns) then
         !$acc exit data delete(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                 oh_buf_send_south, oh_buf_recv_south)
         deallocate (oh_buf_send_north, oh_buf_recv_north)
         deallocate (oh_buf_send_south, oh_buf_recv_south)
         allocate (oh_buf_send_north(need_ns), oh_buf_recv_north(need_ns))
         allocate (oh_buf_send_south(need_ns), oh_buf_recv_south(need_ns))
         !$acc enter data create(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                  oh_buf_send_south, oh_buf_recv_south)
         oh_cap_ns = need_ns
      end if

   end subroutine ocean_halo_reserve

   subroutine ocean_halo_buffers_ensure_nz(nz)
      !! Grow the persistent send/recv buffers to hold at least nz layers.
      !! Called at the top of each 3D exchange routine before pack/unpack.
      !! The 2D routines use the same buffers and always need <= any 3D
      !! capacity, so they are unaffected.
      !!
      !! Grow policy: exit data delete → deallocate → reallocate at new
      !! capacity → enter data create.  The current 2D capacity (initialised
      !! in ocean_halo_init as (ng+1)*ny_total / (ng+1)*nx_total) covers nz=1;
      !! any call with nz>1 that needs more space triggers a one-time grow.
      integer, intent(in) :: nz
         !! Number of vertical layers required

      integer :: need_ew, need_ns

      ! Worst-case 3D capacity:
      !   E/W: face-x send-east uses (ng+1)*nyt per layer ⇒ (ng+1)*nyt*nz
      !   N/S: face-y send-north uses (ng+1)*nxt per layer ⇒ (ng+1)*nxt*nz
      need_ew = (oh_nghost + 1)*oh_ny_total*nz
      need_ns = (oh_nghost + 1)*oh_nx_total*nz

      if (need_ew <= oh_cap_ew .and. need_ns <= oh_cap_ns) return

      call logger%warning("ocean_halo_buffers_ensure_nz: growing halo buffers "// &
                          "mid-run (E/W: "//to_string(oh_cap_ew)//" -> "// &
                          to_string(need_ew)//", N/S: "//to_string(oh_cap_ns)//" -> "// &
                          to_string(need_ns)//"). "// &
                          "Call ocean_halo_reserve(nz, ng_wide) at init to avoid "// &
                          "a mid-run device reallocation (breaks UCX IPC handle reuse).")

      ! Grow E/W buffers if needed
      if (need_ew > oh_cap_ew) then
         !$acc exit data delete(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                 oh_buf_send_west, oh_buf_recv_west)
         deallocate (oh_buf_send_east, oh_buf_recv_east)
         deallocate (oh_buf_send_west, oh_buf_recv_west)
         allocate (oh_buf_send_east(need_ew), oh_buf_recv_east(need_ew))
         allocate (oh_buf_send_west(need_ew), oh_buf_recv_west(need_ew))
         !$acc enter data create(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                  oh_buf_send_west, oh_buf_recv_west)
         oh_cap_ew = need_ew
      end if

      ! Grow N/S buffers if needed
      if (need_ns > oh_cap_ns) then
         !$acc exit data delete(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                 oh_buf_send_south, oh_buf_recv_south)
         deallocate (oh_buf_send_north, oh_buf_recv_north)
         deallocate (oh_buf_send_south, oh_buf_recv_south)
         allocate (oh_buf_send_north(need_ns), oh_buf_recv_north(need_ns))
         allocate (oh_buf_send_south(need_ns), oh_buf_recv_south(need_ns))
         !$acc enter data create(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                  oh_buf_send_south, oh_buf_recv_south)
         oh_cap_ns = need_ns
      end if

   end subroutine ocean_halo_buffers_ensure_nz

   pure function ocean_halo_is_init() result(flag)
      !! True if ocean_halo_init has been called and not yet destroyed.
      logical :: flag
      flag = oh_initialised
   end function ocean_halo_is_init

   ! These predicates gate on DECOMPOSITION (px/py > 1), never on periodicity:
   ! a single-rank periodic domain wraps its ghosts locally, so there is no
   ! rank seam. Two traps they guard against: (1) a setup-time seam fill (e.g.
   ! the wet_mask land-mask exchange) must skip when px=py=1 — which also keeps
   ! non-MPI unit tests off the MPI runtime; (2) on a split axis the exchange
   ! does the periodic wrap over the rank link, so the solver's local
   ! ocean_periodic_wrap_* must be skipped there or it double-wraps with the
   ! wrong partner columns. (D4)

   pure function ocean_halo_is_decomposed() result(flag)
      !! True iff a genuine multi-rank seam exists (initialised, px>1 or py>1).
      logical :: flag
      flag = oh_initialised .and. (oh_decomp%px > 1 .or. oh_decomp%py > 1)
   end function ocean_halo_is_decomposed

   pure function ocean_halo_is_decomposed_x() result(flag)
      !! True iff the X axis is split across ranks (px>1); skip the local x wrap.
      logical :: flag
      flag = oh_initialised .and. oh_decomp%px > 1
   end function ocean_halo_is_decomposed_x

   pure function ocean_halo_is_decomposed_y() result(flag)
      !! Y analogue of ocean_halo_is_decomposed_x (py>1).
      logical :: flag
      flag = oh_initialised .and. oh_decomp%py > 1
   end function ocean_halo_is_decomposed_y

   ! ==================================================================
   ! Centre-2D
   ! ==================================================================

   subroutine ocean_halo_centre_2d(fld, device_resident)
      !! Two-pass X-then-Y halo exchange for a cell-centred 2D field.
      real(wp), intent(inout) :: fld(oh_nx_total, oh_ny_total)
         !! Cell-centred 2D field (nx_total, ny_total)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      call oh_count_centre_2d()
      call ocean_halo_centre_2d_impl(fld, oh_nghost, device_resident)
   end subroutine ocean_halo_centre_2d

   ! ==================================================================
   ! Centre-3D  (batched: one message per direction carrying all nz layers)
   ! ==================================================================

   subroutine ocean_halo_centre_3d(fld, nz, device_resident)
      !! Two-pass X-then-Y halo exchange for a cell-centred 3D field.
      !! Batched: packs ALL nz layers into one buffer per direction and
      !! posts ONE isend+irecv pair per needed direction (not per layer).
      !! Buffer index formula: ((L-1)*rows + (row-1))*width + k_within_strip
      !! where L=layer, rows=nyt, width=ng (centre X-pass).
      integer, intent(in) :: nz
         !! Number of vertical layers
      real(wp), intent(inout) :: fld(oh_nx_total, oh_ny_total, nz)
         !! Cell-centred 3D field (nx_total, ny_total, nz)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      integer :: ng, nxl, nyl, nxt, nyt, j, k, L
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device
      integer :: strip_ew_3d, strip_ns_3d

      call oh_count_centre_3d()

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      ng = oh_nghost
      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = oh_nx_total
      nyt = oh_ny_total
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! Ensure buffers are large enough for nz layers
      call ocean_halo_buffers_ensure_nz(nz)

      ! Batched strip sizes: ng elements per row per layer
      strip_ew_3d = ng*nyt*nz   ! centre X-pass: ng columns per j-row per layer
      strip_ns_3d = nxt*ng*nz   ! centre Y-pass: ng rows per i-col per layer

      ! ---- X pass ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_centre_3d(fld, nxt, nyt, nz, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         ! Pack: buffer index = ((L-1)*nyt + (j-1))*ng + k
         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_send_east, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_east(((L - 1)*nyt + (j - 1))*ng + k) = fld(nxl + k, j, L)
                     end do
                  end do
               end do
            end if
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_send_west, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_west(((L - 1)*nyt + (j - 1))*ng + k) = fld(ng + k, j, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_e) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_east(((L - 1)*nyt + (j - 1))*ng + k) = fld(nxl + k, j, L)
                     end do
                  end do
               end do
            end if
            if (need_w) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_west(((L - 1)*nyt + (j - 1))*ng + k) = fld(ng + k, j, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, strip_ew_3d, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, strip_ew_3d, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, strip_ew_3d, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, strip_ew_3d, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, strip_ew_3d, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, strip_ew_3d, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, strip_ew_3d, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, strip_ew_3d, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_e + need_w); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         ! Unpack: X pass fills west ghosts (i=1..ng) and east ghosts (i=ng+nxl+1..nxt)
         if (on_device) then
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_west, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(k, j, L) = oh_buf_recv_west(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_east, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(ng + nxl + k, j, L) = oh_buf_recv_east(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         else
            if (need_w) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(k, j, L) = oh_buf_recv_west(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
            if (need_e) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(ng + nxl + k, j, L) = oh_buf_recv_east(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass (spans full i=1..nxt including x-filled ghosts) ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_centre_3d(fld, nxt, nyt, nz, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         ! Pack: buffer index = ((L-1)*ng + (k-1))*nxt + j
         if (on_device) then
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_send_north, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        oh_buf_send_north(((L - 1)*ng + (k - 1))*nxt + j) = fld(j, nyl + k, L)
                     end do
                  end do
               end do
            end if
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_send_south, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        oh_buf_send_south(((L - 1)*ng + (k - 1))*nxt + j) = fld(j, ng + k, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_n) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        oh_buf_send_north(((L - 1)*ng + (k - 1))*nxt + j) = fld(j, nyl + k, L)
                     end do
                  end do
               end do
            end if
            if (need_s) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        oh_buf_send_south(((L - 1)*ng + (k - 1))*nxt + j) = fld(j, ng + k, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, strip_ns_3d, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, strip_ns_3d, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, strip_ns_3d, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, strip_ns_3d, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, strip_ns_3d, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, strip_ns_3d, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, strip_ns_3d, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, strip_ns_3d, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_n + need_s); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         ! Unpack: Y pass fills south ghosts (j=1..ng) and north ghosts (j=ng+nyl+1..nyt)
         if (on_device) then
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_south, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        fld(j, k, L) = oh_buf_recv_south(((L - 1)*ng + (k - 1))*nxt + j)
                     end do
                  end do
               end do
            end if
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_north, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        fld(j, ng + nyl + k, L) = oh_buf_recv_north(((L - 1)*ng + (k - 1))*nxt + j)
                     end do
                  end do
               end do
            end if
         else
            if (need_s) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        fld(j, k, L) = oh_buf_recv_south(((L - 1)*ng + (k - 1))*nxt + j)
                     end do
                  end do
               end do
            end if
            if (need_n) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt
                        fld(j, ng + nyl + k, L) = oh_buf_recv_north(((L - 1)*ng + (k - 1))*nxt + j)
                     end do
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_centre_3d

   ! ==================================================================
   ! Face-x-2D  (west rank owns the seam face)
   ! ==================================================================

   subroutine ocean_halo_face_x_2d(fld, device_resident)
      !! Two-pass X-then-Y halo exchange for a 2D x-face field.
      !!
      !! Index map (ng=nghost, nxl=nx_local):
      !!   Send EAST (ng+1 per j-row): source i = nxl+1..ng+nxl+1
      !!     (last i = ng+nxl+1 is owned seam face; first ng are pre-seam interior)
      !!     → fills east-nbr's i = 1..ng+1 (west ghosts + seam-copy overwrite)
      !!
      !!   Recv FROM EAST (ng per j-row):
      !!     east-nbr sends its i = ng+2..2*ng+1 → fills OUR east ghosts i = ng+nxl+2..ng+nxl+ng+1
      !!
      !!   Send WEST (ng per j-row): source i = ng+2..2*ng+1
      !!     → fills west-nbr's east ghosts
      !!
      !!   Recv FROM WEST (ng+1 per j-row):
      !!     west-nbr sends its owned seam + ng interior → fills OUR i = 1..ng+1
      real(wp), intent(inout) :: fld(oh_nx_total + 1, oh_ny_total)
         !! x-face field (nx_total+1, ny_total)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      call oh_count_face_x_2d()
      call ocean_halo_face_x_2d_impl(fld, oh_nghost, device_resident)
   end subroutine ocean_halo_face_x_2d

   ! ==================================================================
   ! Face-x-3D  (batched: one message per direction carrying all nz layers)
   ! ==================================================================

   subroutine ocean_halo_face_x_3d(fld, nz, device_resident)
      !! Two-pass X-then-Y halo exchange for a 3D x-face field.
      !! Batched: packs ALL nz layers into one buffer per direction.
      !! Face-x ownership: west rank sends (ng+1)*nyt*nz eastward,
      !! east rank sends ng*nyt*nz westward (D1 asymmetry preserved).
      !! Y-pass strip is nxt1*ng*nz (full face-x i extent × ng rows × nz).
      integer, intent(in) :: nz
         !! Number of vertical layers
      real(wp), intent(inout) :: fld(oh_nx_total + 1, oh_ny_total, nz)
         !! x-face field (nx_total+1, ny_total, nz)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      integer :: ng, nxl, nyl, nxt, nxt1, nyt, j, k, L
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device

      call oh_count_face_x_3d()

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      ng = oh_nghost
      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = oh_nx_total
      nxt1 = nxt + 1
      nyt = oh_ny_total
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! Ensure buffers are large enough for nz layers
      call ocean_halo_buffers_ensure_nz(nz)

      ! ---- X pass ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_face_x_3d(fld, nxt1, nyt, nz, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         ! Send east: ng+1 per j-row per layer — my LAST ng+1 physical faces
         ! Pack: buffer index = ((L-1)*nyt + (j-1))*(ng+1) + k
         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_send_east, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng + 1
                        oh_buf_send_east(((L - 1)*nyt + (j - 1))*(ng + 1) + k) = fld(nxl + k, j, L)
                     end do
                  end do
               end do
            end if
            ! Send west: ng per j-row per layer (first ng interior faces after west-seam copy)
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_send_west, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_west(((L - 1)*nyt + (j - 1))*ng + k) = fld(ng + 1 + k, j, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_e) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng + 1
                        oh_buf_send_east(((L - 1)*nyt + (j - 1))*(ng + 1) + k) = fld(nxl + k, j, L)
                     end do
                  end do
               end do
            end if
            if (need_w) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        oh_buf_send_west(((L - 1)*nyt + (j - 1))*ng + k) = fld(ng + 1 + k, j, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, (ng + 1)*nyt*nz, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt*nz, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt*nz, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, (ng + 1)*nyt*nz, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, (ng + 1)*nyt*nz, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt*nz, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt*nz, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, (ng + 1)*nyt*nz, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_e + need_w); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            ! Unpack recv-from-west: fills i=1..ng+1 (west ghosts + seam-copy)
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_west, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng + 1
                        fld(k, j, L) = oh_buf_recv_west(((L - 1)*nyt + (j - 1))*(ng + 1) + k)
                     end do
                  end do
               end do
            end if
            ! Unpack recv-from-east: fills east ghosts i=ng+nxl+2..ng+nxl+ng+1
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_east, fld)
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(ng + nxl + 1 + k, j, L) = oh_buf_recv_east(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         else
            if (need_w) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng + 1
                        fld(k, j, L) = oh_buf_recv_west(((L - 1)*nyt + (j - 1))*(ng + 1) + k)
                     end do
                  end do
               end do
            end if
            if (need_e) then
               do L = 1, nz
                  do j = 1, nyt
                     do k = 1, ng
                        fld(ng + nxl + 1 + k, j, L) = oh_buf_recv_east(((L - 1)*nyt + (j - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass (spans full i=1..nxt1 including x-ghosts) ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_face_x_3d(fld, nxt1, nyt, nz, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         ! Pack Y-pass: buffer index = ((L-1)*ng + (k-1))*nxt1 + j
         if (on_device) then
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_send_north, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        oh_buf_send_north(((L - 1)*ng + (k - 1))*nxt1 + j) = fld(j, nyl + k, L)
                     end do
                  end do
               end do
            end if
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_send_south, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        oh_buf_send_south(((L - 1)*ng + (k - 1))*nxt1 + j) = fld(j, ng + k, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_n) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        oh_buf_send_north(((L - 1)*ng + (k - 1))*nxt1 + j) = fld(j, nyl + k, L)
                     end do
                  end do
               end do
            end if
            if (need_s) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        oh_buf_send_south(((L - 1)*ng + (k - 1))*nxt1 + j) = fld(j, ng + k, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, nxt1*ng*nz, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, nxt1*ng*nz, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, nxt1*ng*nz, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, nxt1*ng*nz, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, nxt1*ng*nz, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, nxt1*ng*nz, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, nxt1*ng*nz, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, nxt1*ng*nz, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_n + need_s); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_south, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        fld(j, k, L) = oh_buf_recv_south(((L - 1)*ng + (k - 1))*nxt1 + j)
                     end do
                  end do
               end do
            end if
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_north, fld)
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        fld(j, ng + nyl + k, L) = oh_buf_recv_north(((L - 1)*ng + (k - 1))*nxt1 + j)
                     end do
                  end do
               end do
            end if
         else
            if (need_s) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        fld(j, k, L) = oh_buf_recv_south(((L - 1)*ng + (k - 1))*nxt1 + j)
                     end do
                  end do
               end do
            end if
            if (need_n) then
               do L = 1, nz
                  do k = 1, ng
                     do j = 1, nxt1
                        fld(j, ng + nyl + k, L) = oh_buf_recv_north(((L - 1)*ng + (k - 1))*nxt1 + j)
                     end do
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_face_x_3d

   ! ==================================================================
   ! Face-y-2D  (south rank owns the seam face)
   ! ==================================================================

   subroutine ocean_halo_face_y_2d(fld, device_resident)
      !! Two-pass X-then-Y halo exchange for a 2D y-face field.
      !! Symmetric to face_x_2d with i↔j, x↔y.
      !!
      !! Index map (ng=nghost, nyl=ny_local):
      !!   Send NORTH (ng+1 per i-col): source j = nyl+1..ng+nyl+1
      !!     (last j = ng+nyl+1 is owned seam face; my j == north-nbr's j_N + nyl)
      !!     → fills north-nbr's j = 1..ng+1 (south ghosts + seam-copy)
      !!
      !!   Recv FROM NORTH (ng per i-col):
      !!     north-nbr sends j = ng+2..2*ng+1 → fills OUR north ghosts j = ng+nyl+2..ng+nyl+ng+1
      !!
      !!   Send SOUTH (ng per i-col): source j = ng+2..2*ng+1
      !!   Recv FROM SOUTH (ng+1 per i-col): fills OUR j = 1..ng+1
      real(wp), intent(inout) :: fld(oh_nx_total, oh_ny_total + 1)
         !! y-face field (nx_total, ny_total+1)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      call oh_count_face_y_2d()
      call ocean_halo_face_y_2d_impl(fld, oh_nghost, device_resident)
   end subroutine ocean_halo_face_y_2d

   ! ==================================================================
   ! Face-y-3D  (batched: one message per direction carrying all nz layers)
   ! ==================================================================

   subroutine ocean_halo_face_y_3d(fld, nz, device_resident)
      !! Two-pass X-then-Y halo exchange for a 3D y-face field.
      !! Batched: packs ALL nz layers into one buffer per direction.
      !! Face-y ownership (D1): south rank sends (ng+1)*nxt*nz northward,
      !! north rank sends ng*nxt*nz southward.
      !! X-pass uses standard centre-style (ng per column per layer over nyt1 rows).
      integer, intent(in) :: nz
         !! Number of vertical layers
      real(wp), intent(inout) :: fld(oh_nx_total, oh_ny_total + 1, nz)
         !! y-face field (nx_total, ny_total+1, nz)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory (no OpenACC directives).
         !! Default (.true.) is the normal device-resident path.

      integer :: ng, nxl, nyl, nxt, nyt, nyt1, i, k, L
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device

      call oh_count_face_y_3d()

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      ng = oh_nghost
      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = oh_nx_total
      nyt = oh_ny_total
      nyt1 = nyt + 1
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! Ensure buffers are large enough for nz layers
      call ocean_halo_buffers_ensure_nz(nz)

      ! ---- X pass: standard centre-style, full y-face extent 1..nyt1 ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_face_y_3d(fld, nxt, nyt1, nz, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         ! Pack X-pass: buffer index = ((L-1)*ng + (k-1))*nyt1 + i
         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_send_east, fld)
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        oh_buf_send_east(((L - 1)*ng + (k - 1))*nyt1 + i) = fld(nxl + k, i, L)
                     end do
                  end do
               end do
            end if
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_send_west, fld)
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        oh_buf_send_west(((L - 1)*ng + (k - 1))*nyt1 + i) = fld(ng + k, i, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_e) then
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        oh_buf_send_east(((L - 1)*ng + (k - 1))*nyt1 + i) = fld(nxl + k, i, L)
                     end do
                  end do
               end do
            end if
            if (need_w) then
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        oh_buf_send_west(((L - 1)*ng + (k - 1))*nyt1 + i) = fld(ng + k, i, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, ng*nyt1*nz, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt1*nz, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt1*nz, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, ng*nyt1*nz, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, ng*nyt1*nz, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt1*nz, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt1*nz, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, ng*nyt1*nz, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_e + need_w); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_w) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_west, fld)
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        fld(k, i, L) = oh_buf_recv_west(((L - 1)*ng + (k - 1))*nyt1 + i)
                     end do
                  end do
               end do
            end if
            if (need_e) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_east, fld)
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        fld(ng + nxl + k, i, L) = oh_buf_recv_east(((L - 1)*ng + (k - 1))*nyt1 + i)
                     end do
                  end do
               end do
            end if
         else
            if (need_w) then
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        fld(k, i, L) = oh_buf_recv_west(((L - 1)*ng + (k - 1))*nyt1 + i)
                     end do
                  end do
               end do
            end if
            if (need_e) then
               do L = 1, nz
                  do k = 1, ng
                     do i = 1, nyt1
                        fld(ng + nxl + k, i, L) = oh_buf_recv_east(((L - 1)*ng + (k - 1))*nyt1 + i)
                     end do
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass: ownership — south rank owns seam face ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_face_y_3d(fld, nxt, nyt1, nz, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         ! Pack Y-pass: buffer index = ((L-1)*nxt + (i-1))*(ng+1) + k (send-north, ng+1 per col)
         ! and ((L-1)*nxt + (i-1))*ng + k (send-south, ng per col)
         if (on_device) then
            ! Send north: ng+1 per i-col per layer
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_send_north, fld)
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng + 1
                        oh_buf_send_north(((L - 1)*nxt + (i - 1))*(ng + 1) + k) = fld(i, nyl + k, L)
                     end do
                  end do
               end do
            end if
            ! Send south: ng per i-col per layer
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_send_south, fld)
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng
                        oh_buf_send_south(((L - 1)*nxt + (i - 1))*ng + k) = fld(i, ng + 1 + k, L)
                     end do
                  end do
               end do
            end if
         else
            if (need_n) then
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng + 1
                        oh_buf_send_north(((L - 1)*nxt + (i - 1))*(ng + 1) + k) = fld(i, nyl + k, L)
                     end do
                  end do
               end do
            end if
            if (need_s) then
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng
                        oh_buf_send_south(((L - 1)*nxt + (i - 1))*ng + k) = fld(i, ng + 1 + k, L)
                     end do
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, (ng + 1)*nxt*nz, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, ng*nxt*nz, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, ng*nxt*nz, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, (ng + 1)*nxt*nz, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, (ng + 1)*nxt*nz, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, ng*nxt*nz, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, ng*nxt*nz, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, (ng + 1)*nxt*nz, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         ! nreq = 2*(need_n + need_s); isends = nreq/2 (paired isend+irecv)
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            ! Recv from south (ng+1 per i-col per layer) → fills south ghosts + seam-copy j=1..ng+1
            if (need_s) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_south, fld)
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng + 1
                        fld(i, k, L) = oh_buf_recv_south(((L - 1)*nxt + (i - 1))*(ng + 1) + k)
                     end do
                  end do
               end do
            end if
            ! Recv from north (ng per i-col per layer) → fills north ghosts j=ng+nyl+2..ng+nyl+ng+1
            if (need_n) then
               !$acc parallel loop collapse(3) present(oh_buf_recv_north, fld)
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng
                        fld(i, ng + nyl + 1 + k, L) = oh_buf_recv_north(((L - 1)*nxt + (i - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         else
            if (need_s) then
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng + 1
                        fld(i, k, L) = oh_buf_recv_south(((L - 1)*nxt + (i - 1))*(ng + 1) + k)
                     end do
                  end do
               end do
            end if
            if (need_n) then
               do L = 1, nz
                  do i = 1, nxt
                     do k = 1, ng
                        fld(i, ng + nyl + 1 + k, L) = oh_buf_recv_north(((L - 1)*nxt + (i - 1))*ng + k)
                     end do
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_face_y_3d

   ! ==================================================================
   ! Grouped barotropic exchange
   ! ==================================================================

   subroutine ocean_halo_bt_group_2d(eta, ubt, vbt, device_resident)
      !! Single-call two-pass exchange of three barotropic fields.
      !! Calls the three individual 2D exchanges in order, which each
      !! run their own full two-pass X-then-Y exchange.
      real(wp), intent(inout) :: eta(oh_nx_total, oh_ny_total)
         !! Barotropic sea-surface height (nx_total x ny_total)
      real(wp), intent(inout) :: ubt(oh_nx_total + 1, oh_ny_total)
         !! Barotropic u-transport (nx_total+1 x ny_total)
      real(wp), intent(inout) :: vbt(oh_nx_total, oh_ny_total + 1)
         !! Barotropic v-transport (nx_total x ny_total+1)
      logical, intent(in), optional :: device_resident
         !! Forwarded to the three individual exchanges.  See centre_2d.

      call oh_count_bt_group()
      call oh_count_suppress_on()
      call ocean_halo_centre_2d(eta, device_resident)
      call ocean_halo_face_x_2d(ubt, device_resident)
      call ocean_halo_face_y_2d(vbt, device_resident)
      call oh_count_suppress_off()
   end subroutine ocean_halo_bt_group_2d

   ! ==================================================================
   ! Wide halo primitives (Phase 3a)
   ! ==================================================================
   !
   ! Each wide routine is a thin wrapper that validates ng_wide, ensures
   ! buffer capacity, bumps the SAME counter category as the matching
   ! normal primitive (no new counter categories), then delegates to the
   ! shared _impl body with ng_wide.  The normal primitives also delegate
   ! to the same _impl body with oh_nghost — ensuring the normal path
   ! executes EXACTLY the same code as before (bit-identity guarantee).
   !
   ! Buffer capacity for a wide exchange (ng = ng_wide):
   !   E/W: (ng+1)*nyt_wide   where nyt_wide = oh_ny_local + 2*ng_wide
   !   N/S: (ng+1)*nxt_wide   where nxt_wide = oh_nx_local + 2*ng_wide
   !
   ! FAIL-LOUD guards (each wide entry point):
   !   ng_wide < oh_nghost        => the ghost band would be narrower than
   !                                 the base halo, leaving stale cells
   !   px > 1 .and. ng_wide > nxl => send strip exceeds neighbour interior
   !   py > 1 .and. ng_wide > nyl => same in y
   ! ==================================================================

   subroutine ocean_halo_buffers_ensure_wide(ng_wide)
      !! Grow the persistent send/recv buffers to hold at least a wide
      !! 2D exchange at ghost width ng_wide.
      !! Uses the same exit-delete / dealloc / alloc / enter-create pattern
      !! as ocean_halo_buffers_ensure_nz.  Safe to call with ng_wide ==
      !! oh_nghost (no-op if already large enough).
      integer, intent(in) :: ng_wide
         !! Wide ghost cell width required

      integer :: nyt_wide, nxt_wide, need_ew, need_ns

      nyt_wide = oh_ny_local + 2*ng_wide
      nxt_wide = oh_nx_local + 2*ng_wide
      ! face-x send-east uses (ng_wide+1)*nyt_wide; face-y send-north uses
      ! (ng_wide+1)*nxt_wide.
      need_ew = (ng_wide + 1)*nyt_wide
      need_ns = (ng_wide + 1)*nxt_wide

      if (need_ew <= oh_cap_ew .and. need_ns <= oh_cap_ns) return

      call logger%warning("ocean_halo_buffers_ensure_wide: growing halo buffers "// &
                          "mid-run (E/W: "//to_string(oh_cap_ew)//" -> "// &
                          to_string(need_ew)//", N/S: "//to_string(oh_cap_ns)//" -> "// &
                          to_string(need_ns)//"). "// &
                          "Call ocean_halo_reserve(nz, ng_wide) at init to avoid "// &
                          "a mid-run device reallocation (breaks UCX IPC handle reuse).")

      if (need_ew > oh_cap_ew) then
         !$acc exit data delete(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                 oh_buf_send_west, oh_buf_recv_west)
         deallocate (oh_buf_send_east, oh_buf_recv_east)
         deallocate (oh_buf_send_west, oh_buf_recv_west)
         allocate (oh_buf_send_east(need_ew), oh_buf_recv_east(need_ew))
         allocate (oh_buf_send_west(need_ew), oh_buf_recv_west(need_ew))
         !$acc enter data create(oh_buf_send_east, oh_buf_recv_east, &
         !$acc&                  oh_buf_send_west, oh_buf_recv_west)
         oh_cap_ew = need_ew
      end if

      if (need_ns > oh_cap_ns) then
         !$acc exit data delete(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                 oh_buf_send_south, oh_buf_recv_south)
         deallocate (oh_buf_send_north, oh_buf_recv_north)
         deallocate (oh_buf_send_south, oh_buf_recv_south)
         allocate (oh_buf_send_north(need_ns), oh_buf_recv_north(need_ns))
         allocate (oh_buf_send_south(need_ns), oh_buf_recv_south(need_ns))
         !$acc enter data create(oh_buf_send_north, oh_buf_recv_north, &
         !$acc&                  oh_buf_send_south, oh_buf_recv_south)
         oh_cap_ns = need_ns
      end if

   end subroutine ocean_halo_buffers_ensure_wide

   ! ------------------------------------------------------------------
   ! Centre-2D wide public entry
   ! ------------------------------------------------------------------

   subroutine ocean_halo_centre_2d_wide(fld, ng_wide, device_resident)
      !! Width-parameterized centre-2D halo exchange.
      !! Identical algorithm to ocean_halo_centre_2d but operates on a
      !! wider ghost band ng_wide (>= oh_nghost).
      !! Fail-loud guards prevent ng_wide from exceeding the neighbour interior.
      integer, intent(in) :: ng_wide
         !! Wide ghost cell width (>= oh_nghost)
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng_wide, oh_ny_local + 2*ng_wide)
         !! Cell-centred 2D field (nx_local + 2*ng_wide, ny_local + 2*ng_wide)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      call oh_count_centre_2d()
      call wide_halo_guards(ng_wide, "ocean_halo_centre_2d_wide")
      call ocean_halo_buffers_ensure_wide(ng_wide)
      call ocean_halo_centre_2d_impl(fld, ng_wide, device_resident)
   end subroutine ocean_halo_centre_2d_wide

   ! ------------------------------------------------------------------
   ! Face-x-2D wide public entry
   ! ------------------------------------------------------------------

   subroutine ocean_halo_face_x_2d_wide(fld, ng_wide, device_resident)
      !! Width-parameterized face-x-2D halo exchange.
      integer, intent(in) :: ng_wide
         !! Wide ghost cell width (>= oh_nghost)
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng_wide + 1, oh_ny_local + 2*ng_wide)
         !! x-face 2D field (nx_local + 2*ng_wide + 1, ny_local + 2*ng_wide)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      call oh_count_face_x_2d()
      call wide_halo_guards(ng_wide, "ocean_halo_face_x_2d_wide")
      call ocean_halo_buffers_ensure_wide(ng_wide)
      call ocean_halo_face_x_2d_impl(fld, ng_wide, device_resident)
   end subroutine ocean_halo_face_x_2d_wide

   ! ------------------------------------------------------------------
   ! Face-y-2D wide public entry
   ! ------------------------------------------------------------------

   subroutine ocean_halo_face_y_2d_wide(fld, ng_wide, device_resident)
      !! Width-parameterized face-y-2D halo exchange.
      integer, intent(in) :: ng_wide
         !! Wide ghost cell width (>= oh_nghost)
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng_wide, oh_ny_local + 2*ng_wide + 1)
         !! y-face 2D field (nx_local + 2*ng_wide, ny_local + 2*ng_wide + 1)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      call oh_count_face_y_2d()
      call wide_halo_guards(ng_wide, "ocean_halo_face_y_2d_wide")
      call ocean_halo_buffers_ensure_wide(ng_wide)
      call ocean_halo_face_y_2d_impl(fld, ng_wide, device_resident)
   end subroutine ocean_halo_face_y_2d_wide

   ! ------------------------------------------------------------------
   ! Grouped barotropic wide exchange
   ! ------------------------------------------------------------------

   subroutine ocean_halo_bt_group_2d_wide(eta, ubt, vbt, ng_wide, device_resident)
      !! Width-parameterized grouped barotropic exchange (eta+ubt+vbt).
      !! Bumps oh_count_bt_group; suppresses the three inner primitive
      !! counters exactly as ocean_halo_bt_group_2d does.
      integer, intent(in) :: ng_wide
         !! Wide ghost cell width (>= oh_nghost)
      real(wp), intent(inout) :: eta(oh_nx_local + 2*ng_wide, oh_ny_local + 2*ng_wide)
         !! Barotropic SSH (nx_local+2*ng_wide, ny_local+2*ng_wide)
      real(wp), intent(inout) :: ubt(oh_nx_local + 2*ng_wide + 1, oh_ny_local + 2*ng_wide)
         !! Barotropic u-transport (nx_local+2*ng_wide+1, ny_local+2*ng_wide)
      real(wp), intent(inout) :: vbt(oh_nx_local + 2*ng_wide, oh_ny_local + 2*ng_wide + 1)
         !! Barotropic v-transport (nx_local+2*ng_wide, ny_local+2*ng_wide+1)
      logical, intent(in), optional :: device_resident
         !! Forwarded to the three individual exchanges.  See centre_2d.

      call oh_count_bt_group()
      call wide_halo_guards(ng_wide, "ocean_halo_bt_group_2d_wide")
      call ocean_halo_buffers_ensure_wide(ng_wide)
      call oh_count_suppress_on()
      call ocean_halo_centre_2d_impl(eta, ng_wide, device_resident)
      call ocean_halo_face_x_2d_impl(ubt, ng_wide, device_resident)
      call ocean_halo_face_y_2d_impl(vbt, ng_wide, device_resident)
      call oh_count_suppress_off()
   end subroutine ocean_halo_bt_group_2d_wide

   ! ==================================================================
   ! Private helpers
   ! ==================================================================

   pure subroutine needs_flags(need_w, need_e, need_s, need_n)
      !! Compute which directions need actual MPI or local-wrap work.
      logical, intent(out) :: need_w, need_e, need_s, need_n

      need_e = (.not. oh_decomp%has_east) .or. (oh_periodic_x .and. oh_decomp%px > 1)
      need_w = (.not. oh_decomp%has_west) .or. (oh_periodic_x .and. oh_decomp%px > 1)
      need_n = (.not. oh_decomp%has_north) .or. (oh_periodic_y .and. oh_decomp%py > 1)
      need_s = (.not. oh_decomp%has_south) .or. (oh_periodic_y .and. oh_decomp%py > 1)
   end subroutine needs_flags

   pure function ew_rank_west() result(r)
      !! West neighbour rank with periodic wrap-around.
      integer :: r
      integer :: rx_w
      rx_w = oh_decomp%rx - 1
      if (rx_w < 0) rx_w = oh_decomp%px - 1
      r = decomp_rank_from_coords(oh_decomp%px, rx_w, oh_decomp%ry)
   end function ew_rank_west

   pure function ew_rank_east() result(r)
      !! East neighbour rank with periodic wrap-around.
      integer :: r
      integer :: rx_e
      rx_e = oh_decomp%rx + 1
      if (rx_e >= oh_decomp%px) rx_e = 0
      r = decomp_rank_from_coords(oh_decomp%px, rx_e, oh_decomp%ry)
   end function ew_rank_east

   pure function ns_rank_south() result(r)
      !! South neighbour rank with periodic wrap-around.
      integer :: r
      integer :: ry_s
      ry_s = oh_decomp%ry - 1
      if (ry_s < 0) ry_s = oh_decomp%py - 1
      r = decomp_rank_from_coords(oh_decomp%px, oh_decomp%rx, ry_s)
   end function ns_rank_south

   pure function ns_rank_north() result(r)
      !! North neighbour rank with periodic wrap-around.
      integer :: r
      integer :: ry_n
      ry_n = oh_decomp%ry + 1
      if (ry_n >= oh_decomp%py) ry_n = 0
      r = decomp_rank_from_coords(oh_decomp%px, oh_decomp%rx, ry_n)
   end function ns_rank_north

   subroutine wide_halo_guards(ng_wide, caller)
      !! Fail-loud validation for wide halo entry points.
      !! Called at the top of every *_wide public routine.
      integer, intent(in) :: ng_wide
         !! Requested wide ghost width
      character(len=*), intent(in) :: caller
         !! Caller name for error messages

      if (ng_wide < oh_nghost) then
         call logger%error(caller//": ng_wide="//to_string(ng_wide)// &
                           " < oh_nghost="//to_string(oh_nghost)// &
                           ". A wide exchange narrower than the base ghost band "// &
                           "would leave stale ghost cells. Use ng_wide >= oh_nghost.")
         error stop "wide halo: ng_wide < oh_nghost"
      end if
      if (oh_decomp%px > 1 .and. ng_wide > oh_nx_local) then
         call logger%error(caller//": ng_wide="//to_string(ng_wide)// &
                           " > nx_local="//to_string(oh_nx_local)// &
                           " with px="//to_string(oh_decomp%px)//". "// &
                           "Wide halo width exceeds local interior; decompose "// &
                           "less or lower bt_halo.")
         error stop "wide halo width exceeds local interior (x)"
      end if
      if (oh_decomp%py > 1 .and. ng_wide > oh_ny_local) then
         call logger%error(caller//": ng_wide="//to_string(ng_wide)// &
                           " > ny_local="//to_string(oh_ny_local)// &
                           " with py="//to_string(oh_decomp%py)//". "// &
                           "Wide halo width exceeds local interior; decompose "// &
                           "less or lower bt_halo.")
         error stop "wide halo width exceeds local interior (y)"
      end if
      ! Single-rank periodic guard: the wrap reads i + nx_phys, so ng_wide
      ! must not exceed nx_phys (= oh_nx_local).
      if (oh_decomp%px == 1 .and. oh_periodic_x .and. ng_wide > oh_nx_local) then
         call logger%error(caller//": ng_wide="//to_string(ng_wide)// &
                           " > nx_phys="//to_string(oh_nx_local)// &
                           " for single-rank periodic-x. "// &
                           "The local wrap reads i+nx_phys which would exceed "// &
                           "the array. Reduce ng_wide.")
         error stop "wide halo: ng_wide > nx_phys for single-rank periodic-x"
      end if
      if (oh_decomp%py == 1 .and. oh_periodic_y .and. ng_wide > oh_ny_local) then
         call logger%error(caller//": ng_wide="//to_string(ng_wide)// &
                           " > ny_phys="//to_string(oh_ny_local)// &
                           " for single-rank periodic-y. "// &
                           "The local wrap reads j+ny_phys which would exceed "// &
                           "the array. Reduce ng_wide.")
         error stop "wide halo: ng_wide > ny_phys for single-rank periodic-y"
      end if
   end subroutine wide_halo_guards

   ! ------------------------------------------------------------------
   ! Centre-2D impl (ng as explicit argument — shared by normal + wide)
   ! ------------------------------------------------------------------

   subroutine ocean_halo_centre_2d_impl(fld, ng, device_resident)
      !! Two-pass X-then-Y halo exchange for a cell-centred 2D field.
      !! Core body shared by ocean_halo_centre_2d (ng=oh_nghost) and
      !! ocean_halo_centre_2d_wide (ng=ng_wide).
      integer, intent(in) :: ng
         !! Ghost cell width for this exchange
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng, oh_ny_local + 2*ng)
         !! Cell-centred 2D field
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      integer :: nxl, nyl, nxt, nyt, j, k
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: strip_ew, strip_ns, nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = nxl + 2*ng
      nyt = nyl + 2*ng
      strip_ew = ng*nyt
      strip_ns = nxt*ng
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! ---- X pass ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_centre_2d(fld, nxt, nyt, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_send_east, fld)
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_east((j - 1)*ng + k) = fld(nxl + k, j)
                  end do
               end do
            end if
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_send_west, fld)
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_west((j - 1)*ng + k) = fld(ng + k, j)
                  end do
               end do
            end if
         else
            if (need_e) then
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_east((j - 1)*ng + k) = fld(nxl + k, j)
                  end do
               end do
            end if
            if (need_w) then
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_west((j - 1)*ng + k) = fld(ng + k, j)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, strip_ew, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, strip_ew, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, strip_ew, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, strip_ew, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, strip_ew, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, strip_ew, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, strip_ew, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, strip_ew, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_west, fld)
               do j = 1, nyt
                  do k = 1, ng
                     fld(k, j) = oh_buf_recv_west((j - 1)*ng + k)
                  end do
               end do
            end if
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_east, fld)
               do j = 1, nyt
                  do k = 1, ng
                     fld(ng + nxl + k, j) = oh_buf_recv_east((j - 1)*ng + k)
                  end do
               end do
            end if
         else
            if (need_w) then
               do j = 1, nyt
                  do k = 1, ng
                     fld(k, j) = oh_buf_recv_west((j - 1)*ng + k)
                  end do
               end do
            end if
            if (need_e) then
               do j = 1, nyt
                  do k = 1, ng
                     fld(ng + nxl + k, j) = oh_buf_recv_east((j - 1)*ng + k)
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass (spans full i=1..nxt including x-filled ghosts) ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_centre_2d(fld, nxt, nyt, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         if (on_device) then
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_send_north, fld)
               do k = 1, ng
                  do j = 1, nxt
                     oh_buf_send_north((k - 1)*nxt + j) = fld(j, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_send_south, fld)
               do k = 1, ng
                  do j = 1, nxt
                     oh_buf_send_south((k - 1)*nxt + j) = fld(j, ng + k)
                  end do
               end do
            end if
         else
            if (need_n) then
               do k = 1, ng
                  do j = 1, nxt
                     oh_buf_send_north((k - 1)*nxt + j) = fld(j, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               do k = 1, ng
                  do j = 1, nxt
                     oh_buf_send_south((k - 1)*nxt + j) = fld(j, ng + k)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, strip_ns, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, strip_ns, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, strip_ns, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, strip_ns, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, strip_ns, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, strip_ns, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, strip_ns, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, strip_ns, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_south, fld)
               do k = 1, ng
                  do j = 1, nxt
                     fld(j, k) = oh_buf_recv_south((k - 1)*nxt + j)
                  end do
               end do
            end if
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_north, fld)
               do k = 1, ng
                  do j = 1, nxt
                     fld(j, ng + nyl + k) = oh_buf_recv_north((k - 1)*nxt + j)
                  end do
               end do
            end if
         else
            if (need_s) then
               do k = 1, ng
                  do j = 1, nxt
                     fld(j, k) = oh_buf_recv_south((k - 1)*nxt + j)
                  end do
               end do
            end if
            if (need_n) then
               do k = 1, ng
                  do j = 1, nxt
                     fld(j, ng + nyl + k) = oh_buf_recv_north((k - 1)*nxt + j)
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_centre_2d_impl

   ! ------------------------------------------------------------------
   ! Face-x-2D impl (ng as explicit argument)
   ! ------------------------------------------------------------------

   subroutine ocean_halo_face_x_2d_impl(fld, ng, device_resident)
      !! Two-pass X-then-Y halo exchange for a 2D x-face field.
      !! Core body shared by ocean_halo_face_x_2d (ng=oh_nghost) and
      !! ocean_halo_face_x_2d_wide (ng=ng_wide).
      integer, intent(in) :: ng
         !! Ghost cell width for this exchange
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng + 1, oh_ny_local + 2*ng)
         !! x-face field (nx_local + 2*ng + 1, ny_local + 2*ng)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      integer :: nxl, nyl, nxt, nxt1, nyt, j, k
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = nxl + 2*ng
      nxt1 = nxt + 1
      nyt = nyl + 2*ng
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! ---- X pass ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_face_x_2d(fld, nxt1, nyt, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_send_east, fld)
               do j = 1, nyt
                  do k = 1, ng + 1
                     oh_buf_send_east((j - 1)*(ng + 1) + k) = fld(nxl + k, j)
                  end do
               end do
            end if
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_send_west, fld)
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_west((j - 1)*ng + k) = fld(ng + 1 + k, j)
                  end do
               end do
            end if
         else
            if (need_e) then
               do j = 1, nyt
                  do k = 1, ng + 1
                     oh_buf_send_east((j - 1)*(ng + 1) + k) = fld(nxl + k, j)
                  end do
               end do
            end if
            if (need_w) then
               do j = 1, nyt
                  do k = 1, ng
                     oh_buf_send_west((j - 1)*ng + k) = fld(ng + 1 + k, j)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, (ng + 1)*nyt, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, (ng + 1)*nyt, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, (ng + 1)*nyt, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, (ng + 1)*nyt, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_west, fld)
               do j = 1, nyt
                  do k = 1, ng + 1
                     fld(k, j) = oh_buf_recv_west((j - 1)*(ng + 1) + k)
                  end do
               end do
            end if
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_east, fld)
               do j = 1, nyt
                  do k = 1, ng
                     fld(ng + nxl + 1 + k, j) = oh_buf_recv_east((j - 1)*ng + k)
                  end do
               end do
            end if
         else
            if (need_w) then
               do j = 1, nyt
                  do k = 1, ng + 1
                     fld(k, j) = oh_buf_recv_west((j - 1)*(ng + 1) + k)
                  end do
               end do
            end if
            if (need_e) then
               do j = 1, nyt
                  do k = 1, ng
                     fld(ng + nxl + 1 + k, j) = oh_buf_recv_east((j - 1)*ng + k)
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass (spans full i=1..nxt1 including x-ghosts) ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_face_x_2d(fld, nxt1, nyt, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         if (on_device) then
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_send_north, fld)
               do k = 1, ng
                  do j = 1, nxt1
                     oh_buf_send_north((k - 1)*nxt1 + j) = fld(j, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_send_south, fld)
               do k = 1, ng
                  do j = 1, nxt1
                     oh_buf_send_south((k - 1)*nxt1 + j) = fld(j, ng + k)
                  end do
               end do
            end if
         else
            if (need_n) then
               do k = 1, ng
                  do j = 1, nxt1
                     oh_buf_send_north((k - 1)*nxt1 + j) = fld(j, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               do k = 1, ng
                  do j = 1, nxt1
                     oh_buf_send_south((k - 1)*nxt1 + j) = fld(j, ng + k)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, nxt1*ng, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, nxt1*ng, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, nxt1*ng, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, nxt1*ng, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, nxt1*ng, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, nxt1*ng, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, nxt1*ng, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, nxt1*ng, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_south, fld)
               do k = 1, ng
                  do j = 1, nxt1
                     fld(j, k) = oh_buf_recv_south((k - 1)*nxt1 + j)
                  end do
               end do
            end if
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_north, fld)
               do k = 1, ng
                  do j = 1, nxt1
                     fld(j, ng + nyl + k) = oh_buf_recv_north((k - 1)*nxt1 + j)
                  end do
               end do
            end if
         else
            if (need_s) then
               do k = 1, ng
                  do j = 1, nxt1
                     fld(j, k) = oh_buf_recv_south((k - 1)*nxt1 + j)
                  end do
               end do
            end if
            if (need_n) then
               do k = 1, ng
                  do j = 1, nxt1
                     fld(j, ng + nyl + k) = oh_buf_recv_north((k - 1)*nxt1 + j)
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_face_x_2d_impl

   ! ------------------------------------------------------------------
   ! Face-y-2D impl (ng as explicit argument)
   ! ------------------------------------------------------------------

   subroutine ocean_halo_face_y_2d_impl(fld, ng, device_resident)
      !! Two-pass X-then-Y halo exchange for a 2D y-face field.
      !! Core body shared by ocean_halo_face_y_2d (ng=oh_nghost) and
      !! ocean_halo_face_y_2d_wide (ng=ng_wide).
      integer, intent(in) :: ng
         !! Ghost cell width for this exchange
      real(wp), intent(inout) :: fld(oh_nx_local + 2*ng, oh_ny_local + 2*ng + 1)
         !! y-face field (nx_local + 2*ng, ny_local + 2*ng + 1)
      logical, intent(in), optional :: device_resident
         !! If .false., operate on host memory.  Default = .true.

      integer :: nxl, nyl, nxt, nyt, nyt1, i, k
      integer :: rk_e, rk_w, rk_n, rk_s
      integer :: nreq
      type(comm_t) :: comm
      type(request_t) :: reqs(MAX_REQS)
      type(MPI_Status) :: stats(MAX_REQS)
      logical :: need_e, need_w, need_n, need_s
      logical :: on_device

      on_device = .true.
      if (present(device_resident)) on_device = device_resident

      nxl = oh_nx_local
      nyl = oh_ny_local
      nxt = nxl + 2*ng
      nyt = nyl + 2*ng
      nyt1 = nyt + 1
      if (oh_decomp%px > 1 .or. oh_decomp%py > 1) comm = comm_env_compute_comm()

      call needs_flags(need_w, need_e, need_s, need_n)

      ! ---- X pass: standard centre-style, full y-face extent 1..nyt1 ----
      if (oh_decomp%px == 1) then
         if (oh_periodic_x) then
            call ocean_periodic_wrap_face_y_2d(fld, nxt, nyt1, nxl, nyl, ng, .true., .false.)
         end if
      else
         if (need_e) rk_e = ew_rank_east()
         if (need_w) rk_w = ew_rank_west()

         if (on_device) then
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_send_east, fld)
               do k = 1, ng
                  do i = 1, nyt1
                     oh_buf_send_east((k - 1)*nyt1 + i) = fld(nxl + k, i)
                  end do
               end do
            end if
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_send_west, fld)
               do k = 1, ng
                  do i = 1, nyt1
                     oh_buf_send_west((k - 1)*nyt1 + i) = fld(ng + k, i)
                  end do
               end do
            end if
         else
            if (need_e) then
               do k = 1, ng
                  do i = 1, nyt1
                     oh_buf_send_east((k - 1)*nyt1 + i) = fld(nxl + k, i)
                  end do
               end do
            end if
            if (need_w) then
               do k = 1, ng
                  do i = 1, nyt1
                     oh_buf_send_west((k - 1)*nyt1 + i) = fld(ng + k, i)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_e) then
               !$acc host_data use_device(oh_buf_send_east, oh_buf_recv_east)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, ng*nyt1, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt1, rk_e, TAG_OC_E_TO_W, reqs(nreq))
               !$acc end host_data
            end if
            if (need_w) then
               !$acc host_data use_device(oh_buf_send_west, oh_buf_recv_west)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt1, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, ng*nyt1, rk_w, TAG_OC_W_TO_E, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_e) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_east, ng*nyt1, rk_e, TAG_OC_W_TO_E, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_east, ng*nyt1, rk_e, TAG_OC_E_TO_W, reqs(nreq))
            end if
            if (need_w) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_west, ng*nyt1, rk_w, TAG_OC_E_TO_W, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_west, ng*nyt1, rk_w, TAG_OC_W_TO_E, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_w) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_west, fld)
               do k = 1, ng
                  do i = 1, nyt1
                     fld(k, i) = oh_buf_recv_west((k - 1)*nyt1 + i)
                  end do
               end do
            end if
            if (need_e) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_east, fld)
               do k = 1, ng
                  do i = 1, nyt1
                     fld(ng + nxl + k, i) = oh_buf_recv_east((k - 1)*nyt1 + i)
                  end do
               end do
            end if
         else
            if (need_w) then
               do k = 1, ng
                  do i = 1, nyt1
                     fld(k, i) = oh_buf_recv_west((k - 1)*nyt1 + i)
                  end do
               end do
            end if
            if (need_e) then
               do k = 1, ng
                  do i = 1, nyt1
                     fld(ng + nxl + k, i) = oh_buf_recv_east((k - 1)*nyt1 + i)
                  end do
               end do
            end if
         end if
      end if

      ! ---- Y pass: ownership — south rank owns seam face ----
      if (oh_decomp%py == 1) then
         if (oh_periodic_y) then
            call ocean_periodic_wrap_face_y_2d(fld, nxt, nyt1, nxl, nyl, ng, .false., .true.)
         end if
      else
         if (need_n) rk_n = ns_rank_north()
         if (need_s) rk_s = ns_rank_south()

         if (on_device) then
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_send_north, fld)
               do i = 1, nxt
                  do k = 1, ng + 1
                     oh_buf_send_north((i - 1)*(ng + 1) + k) = fld(i, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_send_south, fld)
               do i = 1, nxt
                  do k = 1, ng
                     oh_buf_send_south((i - 1)*ng + k) = fld(i, ng + 1 + k)
                  end do
               end do
            end if
         else
            if (need_n) then
               do i = 1, nxt
                  do k = 1, ng + 1
                     oh_buf_send_north((i - 1)*(ng + 1) + k) = fld(i, nyl + k)
                  end do
               end do
            end if
            if (need_s) then
               do i = 1, nxt
                  do k = 1, ng
                     oh_buf_send_south((i - 1)*ng + k) = fld(i, ng + 1 + k)
                  end do
               end do
            end if
         end if

         nreq = 0
         if (on_device) then
            if (need_n) then
               !$acc host_data use_device(oh_buf_send_north, oh_buf_recv_north)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, (ng + 1)*nxt, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, ng*nxt, rk_n, TAG_OC_N_TO_S, reqs(nreq))
               !$acc end host_data
            end if
            if (need_s) then
               !$acc host_data use_device(oh_buf_send_south, oh_buf_recv_south)
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, ng*nxt, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, (ng + 1)*nxt, rk_s, TAG_OC_S_TO_N, reqs(nreq))
               !$acc end host_data
            end if
         else
            if (need_n) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_north, (ng + 1)*nxt, rk_n, TAG_OC_S_TO_N, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_north, ng*nxt, rk_n, TAG_OC_N_TO_S, reqs(nreq))
            end if
            if (need_s) then
               nreq = nreq + 1
               call HALO_ISEND_N(comm, oh_buf_send_south, ng*nxt, rk_s, TAG_OC_N_TO_S, reqs(nreq))
               nreq = nreq + 1
               call HALO_IRECV_N(comm, oh_buf_recv_south, (ng + 1)*nxt, rk_s, TAG_OC_S_TO_N, reqs(nreq))
            end if
         end if
         if (nreq > 0) then
            call oh_count_msgs(nreq/2)
            call waitall(reqs(1:nreq), stats(1:nreq))
         end if

         if (on_device) then
            if (need_s) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_south, fld)
               do i = 1, nxt
                  do k = 1, ng + 1
                     fld(i, k) = oh_buf_recv_south((i - 1)*(ng + 1) + k)
                  end do
               end do
            end if
            if (need_n) then
               !$acc parallel loop collapse(2) present(oh_buf_recv_north, fld)
               do i = 1, nxt
                  do k = 1, ng
                     fld(i, ng + nyl + 1 + k) = oh_buf_recv_north((i - 1)*ng + k)
                  end do
               end do
            end if
         else
            if (need_s) then
               do i = 1, nxt
                  do k = 1, ng + 1
                     fld(i, k) = oh_buf_recv_south((i - 1)*(ng + 1) + k)
                  end do
               end do
            end if
            if (need_n) then
               do i = 1, nxt
                  do k = 1, ng
                     fld(i, ng + nyl + 1 + k) = oh_buf_recv_north((i - 1)*ng + k)
                  end do
               end do
            end if
         end if
      end if

   end subroutine ocean_halo_face_y_2d_impl

end module rdb_ocean_halo
