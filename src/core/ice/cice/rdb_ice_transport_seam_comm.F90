!! Communication side of the CICE-port transport's decomposition seam (roundabout-original).
module rdb_ice_transport_seam_comm
   !! Binds `transport_seam_t` (`rdb_ice_transport_seam`, the standalone ice library) to the ocean's
   !! communication layer and holds the one halo the transport needs that no existing routine
   !! provides as a bundle: the category state ring.  Part of the MAIN library (it uses the ocean's
   !! halo, fold and reductions).
   !!
   !! **`ice_transport_seam_bind`** fills a seam from the rank's boundary state: the collectives
   !! (`halo_allreduce_*` over the compute communicator, an identity on one rank), `fold` from the
   !! rank-local `bc%north_fold`, `fold_dist` from the ocean's fold dispatch
   !! (`ocean_fold_is_distributed`: `px > 1`), and the fold-row pairing of a split fold through the
   !! ocean's owner-routed fold exchange with its `row_shift` option: the ice's north-face flux
   !! `N(i,j)` is stored at the row of its cell `j`, one below the y-face row of the same face, so
   !! the exchange is told `row_shift = -1` and moves the fold-line row `ng+nj` as the v-face fold
   !! line `ng+nj+1`, with the ocean's own rule (east half authoritative, west half the negated
   !! mirror, a self-conjugate column zeroed).  The north ghost rows of the flux buffers are also
   !! written by that exchange; nothing reads them.  The exchange buffers are reserved for the
   !! widest batch at bind time, so the step allocates nothing.
   !!
   !! **`ice_halo_category_state`** refreshes the ghost ring of the category state and of `aice0`
   !! (the cell-centred scalars `aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice0`): the
   !! ocean's two-pass centre exchange on a split axis, the local periodic wrap (the kernel the single-rank
   !! hooks always used) on one rank; a wall edge keeps its ghosts, and then, on the north fold row, the tripolar fold as a plain
   !! copy (scalars: no sign).  **`ice_halo_centre_pair`** is the same refresh for a TRUE-VECTOR pair held at cell
   !! centres (the T-cell wind stress the face hand-off averages): its fold negates both components
   !! (`ocean_fold_wrap_centre_vector_pair`, one exchange group on a split fold row); a scalar copy there would hand
   !! the north-face stress of the fold row the wrong sign.  It is the multi-rank form of the single-rank wrap of
   !! `rdb_ice_cice_hooks` and is what the transport asks for before its call (two rings valid).
   use, intrinsic :: iso_fortran_env, only: int64, int32
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_halo, only: halo_allreduce_sum_vec, halo_allreduce_efp_list, halo_allreduce_min
   use rdb_ocean_halo, only: ocean_halo_centre, ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_3d, ocean_periodic_wrap_centre_2d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_centre_flat, ocean_fold_wrap_centre_vector_pair
   use rdb_ocean_fold_exchange, only: ocean_fold_is_distributed, ocean_fold_begin, ocean_fold_pack, &
                                      ocean_fold_exchange, ocean_fold_unpack, ocean_fold_end, &
                                      ocean_fold_exchange_reserve, FOLD_STAG_V
   use rdb_ice_transport_seam, only: transport_seam_t, seam_min_one
   use rdb_ice_errors, only: ice_err_summary_t, ice_err_set_combine, ice_err_combine_iface, ICE_NERR
   implicit none
   private

   public :: ice_transport_seam_bind
   public :: ice_halo_category_state
   public :: ice_halo_centre_pair
   public :: ice_transport_seam_bind_reductions, ice_err_bind_comm, ice_err_combine_collective

   integer, save :: sc_nghost = 0
      !! Ghost width of the grid the seam was bound for (the fold exchange's row count).

contains

   subroutine ice_transport_seam_bind(seam, grid, bc, nbatch, ntrace)
      !! Bind `seam` for this rank.  Call once at setup, after the ocean's halo and fold exchange are
      !! initialised and the boundary state's edges are set (`bc%north_fold`).  `nbatch`, `ntrace`
      !! size the fold-pairing buffers (the largest category batch and the tracer count of the
      !! transport workspace).  Re-binding is allowed.
      type(transport_seam_t), intent(inout) :: seam
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      integer, intent(in) :: nbatch, ntrace

      sc_nghost = grid%nghost
      seam%multirank = ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()
      seam%fold = bc%north_fold
      seam%fold_dist = bc%north_fold .and. ocean_fold_is_distributed()
      seam%sum_i8 => comm_sum_i8
      seam%sum_efp => comm_sum_efp
      seam%pair_fold_dist => comm_pair_fold
      seam%min_i => comm_min_i
      ! the batch size is the minimum over the ranks: tiles differ, and the batch loop and the fold pairing are
      ! collective per batch
      seam%nbatch = nbatch
      call seam_min_one(seam, seam%nbatch)
      if (seam%fold_dist) call ocean_fold_exchange_reserve((sc_nghost + 1)*seam%nbatch*(1 + ntrace))
   end subroutine ice_transport_seam_bind

   subroutine ice_transport_seam_bind_reductions(seam, grid, bc)
      !! The part of `ice_transport_seam_bind` that is not the fold pairing: the rank-uniform reductions
      !! (`multirank`, `sum_i8`, `sum_efp`, `min_i`) of `seam`, no collective and no buffer reservation.  What the
      !! step sequencer uses to make a hook's failure rank-uniform.
      type(transport_seam_t), intent(inout) :: seam
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc

      sc_nghost = grid%nghost
      seam%multirank = ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()
      seam%fold = bc%north_fold
      seam%fold_dist = bc%north_fold .and. ocean_fold_is_distributed()
      seam%sum_i8 => comm_sum_i8
      seam%sum_efp => comm_sum_efp
      seam%min_i => comm_min_i
   end subroutine ice_transport_seam_bind_reductions

   subroutine ice_err_bind_comm()
      !! Bind the cross-rank combination of `ice_err_reduce` on a decomposed run, unbind it on one rank (the
      !! ice error registry cannot use MPI).  Called on every rank at the top of every ice step.
      procedure(ice_err_combine_iface), pointer :: p

      p => null()
      if (ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()) p => ice_err_combine_collective
      call ice_err_set_combine(p)
   end subroutine ice_err_bind_comm

   subroutine ice_err_combine_collective(summary)
      !! Collective combination of one error summary: ONE sum allreduce carries the 32 bits of the OR'd mask, the
      !! per-bit cell counts and the number of ranks that name a first offender; only when some rank does
      !! (a fatal bit: the abnormal path) a min and a sum more fix the global first offender and its word.
      !! Every rank takes the same branches because they depend on the reduced values only.
      type(ice_err_summary_t), intent(inout) :: summary
      integer, parameter :: NB = 32, NV = NB + ICE_NERR + 1
      integer(int64) :: v(NV)
      real(wp) :: xmin, ymin
      real(wp) :: xw(1), yw(1)
      integer :: b

      do b = 0, NB - 1
         v(b + 1) = merge(1_int64, 0_int64, btest(summary%mask, b))
      end do
      do b = 0, ICE_NERR - 1
         v(NB + b + 1) = int(summary%count(b), int64)
      end do
      v(NV) = merge(1_int64, 0_int64, summary%first_idx > 0_int64)
      call comm_sum_i8(NV, v)
      summary%mask = 0_int32
      do b = 0, NB - 1
         if (v(b + 1) > 0_int64) summary%mask = ibset(summary%mask, b)
      end do
      do b = 0, ICE_NERR - 1
         summary%count(b) = int(v(NB + b + 1), int32)
      end do
      if (v(NV) == 0_int64) then
         summary%first_idx = 0_int64
         summary%first_word = 0_int32
         return
      end if
      ! global first offender: the smallest global index (below 2**53, exact in real64); its owner holds the word
      xmin = huge(1.0_wp)
      if (summary%first_idx > 0_int64) xmin = real(summary%first_idx, wp)
      call halo_allreduce_min(xmin, ymin)
      xw(1) = 0.0_wp
      if (summary%first_idx > 0_int64 .and. real(summary%first_idx, wp) == ymin) xw(1) = real(summary%first_word, wp)
      call halo_allreduce_sum_vec(xw, yw, 1)
      summary%first_idx = nint(ymin, int64)
      summary%first_word = nint(yw(1), int32)
   end subroutine ice_err_combine_collective

   subroutine comm_sum_i8(n, v)
      !! Integer sum over the ranks.  Rides the exact `real64` allreduce (counts far below 2**53).
      integer, intent(in) :: n
      integer(int64), intent(inout) :: v(n)
      real(wp) :: loc(n), glob(n)

      loc = real(v, wp)
      call halo_allreduce_sum_vec(loc, glob, n)
      v = nint(glob, int64)
   end subroutine comm_sum_i8

   subroutine comm_min_i(v)
      integer, intent(inout) :: v
      real(wp) :: x, y

      x = real(v, wp)
      call halo_allreduce_min(x, y)
      v = nint(y)
   end subroutine comm_min_i

   subroutine comm_sum_efp(n, list)
      integer, intent(in) :: n
      type(efp_t), intent(inout) :: list(n)
      type(efp_t) :: glob(n)

      call halo_allreduce_efp_list(list, glob, n)
      list = glob
   end subroutine comm_sum_efp

   subroutine comm_pair_fold(nx, ny, nb, ntrace, with_tracers, mflx_n, mtflx_n)
      !! Pair the fold-row north fluxes over a split north rank row (module docstring).
      integer, intent(in) :: nx, ny, nb, ntrace
      logical, intent(in) :: with_tracers
      real(wp), intent(inout) :: mflx_n(nx, ny, nb), mtflx_n(nx, ny, ntrace, nb)

      call ocean_fold_begin((sc_nghost + 1)*nb*(1 + merge(ntrace, 0, with_tracers)))
      call fold_pack_flat(mflx_n, nx, ny, nb)
      if (with_tracers) call fold_pack_flat(mtflx_n, nx, ny, ntrace*nb)
      call ocean_fold_exchange(.true.)
      call fold_unpack_flat(mflx_n, nx, ny, nb)
      if (with_tracers) call fold_unpack_flat(mtflx_n, nx, ny, ntrace*nb)
      call ocean_fold_end()
   end subroutine comm_pair_fold

   subroutine fold_pack_flat(fld, nx, ny, nz)
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: fld(nx, ny, nz)

      call ocean_fold_pack(fld, nx, ny, nz, FOLD_STAG_V, .true., row_shift=-1)
   end subroutine fold_pack_flat

   subroutine fold_unpack_flat(fld, nx, ny, nz)
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: fld(nx, ny, nz)

      call ocean_fold_unpack(fld, nx, ny, nz, FOLD_STAG_V, .true., .true., row_shift=-1)
   end subroutine fold_unpack_flat

   subroutine ice_halo_centre_pair(grid, bc, nx, ny, fa, fb, device_resident, test_fold)
      !! Refresh the ghost ring of two cell-centred 2-D fields (the T-cell wind stress the face hand-off
      !! averages): the ocean's two-pass centre exchange on a split axis, the local periodic wrap on one
      !! rank, then the tripolar fold on the fold row with the VECTOR sign (the pair is a true vector, both components
      !! negate).  Collective on a split axis.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fa(nx, ny), fb(nx, ny)
      logical, intent(in), optional :: device_resident
      integer, intent(in), optional :: test_fold
         !! Test harness only: 1 = fold the pair as a SCALAR copy (the wrong sign), 2 = no fold.  Absent: the vector fold.
      logical :: split
      integer :: tf

      tf = 0
      if (present(test_fold)) tf = test_fold
      split = ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()
      call refresh_one(fa)
      call refresh_one(fb)
      ! the pair is a TRUE VECTOR: both components change sign across the fold (one exchange group on a split fold row)
      if (tf == 0) then
         call ocean_fold_wrap_centre_vector_pair(grid, bc, fa, fb, device_resident)
      else if (tf == 1) then
         call ocean_fold_wrap_centre_flat(grid, bc, fa, nx, ny, 1, device_resident)
         call ocean_fold_wrap_centre_flat(grid, bc, fb, nx, ny, 1, device_resident)
      end if

   contains

      subroutine refresh_one(fld)
         real(wp), intent(inout) :: fld(nx, ny)

         if (split) then
            call ocean_halo_centre(fld, device_resident)
         else
            call ocean_periodic_wrap_centre_2d(fld, nx, ny, grid%nx_phys, grid%ny_phys, grid%nghost, bc%periodic_x, &
                                               bc%periodic_y)
         end if
      end subroutine refresh_one

   end subroutine ice_halo_centre_pair

   subroutine ice_halo_category_state(grid, bc, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, &
                                      qsnon, aice0, device_resident)
      !! Refresh the ghost ring of the category state and `aice0` (module docstring).  Collective on a
      !! split axis; every rank calls it.  Every array device-present on an offload build.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat*nilyr), sicen(nx, ny, ncat*nilyr), qsnon(nx, ny, ncat*nslyr)
      real(wp), intent(inout) :: aice0(nx, ny, 1)
      logical, intent(in), optional :: device_resident
      logical :: split

      split = ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()
      call refresh_flat(aicen, ncat)
      call refresh_flat(vicen, ncat)
      call refresh_flat(vsnon, ncat)
      call refresh_flat(Tsfcn, ncat)
      call refresh_flat(qicen, ncat*nilyr)
      call refresh_flat(sicen, ncat*nilyr)
      call refresh_flat(qsnon, ncat*nslyr)
      call refresh_flat(aice0, 1)

   contains

      subroutine refresh_flat(fld, nz)
         integer, intent(in) :: nz
         real(wp), intent(inout) :: fld(nx, ny, nz)

         if (split) then
            call ocean_halo_centre(fld, nz, device_resident)
         else
            ! one rank: the local wrap, as the single-rank hooks always did (no halo buffers touched)
            call ocean_periodic_wrap_centre_3d(fld, nx, ny, nz, grid%nx_phys, grid%ny_phys, grid%nghost, bc%periodic_x, &
                                               bc%periodic_y)
         end if
         if (bc%north_fold) call ocean_fold_wrap_centre_flat(grid, bc, fld, nx, ny, nz, device_resident)
      end subroutine refresh_flat

   end subroutine ice_halo_category_state

end module rdb_ice_transport_seam_comm
