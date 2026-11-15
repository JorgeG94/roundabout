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
   !! copy (scalars: no sign).  It is the multi-rank form of the single-rank wrap of
   !! `rdb_ice_cice_hooks` and is what the transport asks for before its call (two rings valid).
   use, intrinsic :: iso_fortran_env, only: int64, int32
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_halo, only: halo_allreduce_sum_vec, halo_allreduce_efp_list
   use rdb_ocean_halo, only: ocean_halo_centre, ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_3d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_centre_flat
   use rdb_ocean_fold_exchange, only: ocean_fold_is_distributed, ocean_fold_begin, ocean_fold_pack, &
                                      ocean_fold_exchange, ocean_fold_unpack, ocean_fold_end, &
                                      ocean_fold_exchange_reserve, FOLD_STAG_V
   use rdb_ice_transport_seam, only: transport_seam_t
   implicit none
   private

   public :: ice_transport_seam_bind
   public :: ice_halo_category_state

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
      if (seam%fold_dist) call ocean_fold_exchange_reserve((sc_nghost + 1)*nbatch*(1 + ntrace))
   end subroutine ice_transport_seam_bind

   subroutine comm_sum_i8(n, v)
      !! Integer sum over the ranks.  Rides the exact `real64` allreduce (counts far below 2**53).
      integer, intent(in) :: n
      integer(int64), intent(inout) :: v(n)
      real(wp) :: loc(n), glob(n)

      loc = real(v, wp)
      call halo_allreduce_sum_vec(loc, glob, n)
      v = nint(glob, int64)
   end subroutine comm_sum_i8

   subroutine comm_sum_efp(n, list)
      integer, intent(in) :: n
      type(efp_t), intent(inout) :: list(n)
      type(efp_t) :: glob(n)

      call halo_allreduce_efp_list(list, glob, n)
      list = glob
   end subroutine comm_sum_efp

   subroutine comm_pair_fold(nx, ny, nb, ntrace, mflx_n, mtflx_n)
      !! Pair the fold-row north fluxes over a split north rank row (module docstring).
      integer, intent(in) :: nx, ny, nb, ntrace
      real(wp), intent(inout) :: mflx_n(nx, ny, nb), mtflx_n(nx, ny, ntrace, nb)

      call ocean_fold_begin((sc_nghost + 1)*nb*(1 + ntrace))
      call fold_pack_flat(mflx_n, nx, ny, nb)
      call fold_pack_flat(mtflx_n, nx, ny, ntrace*nb)
      call ocean_fold_exchange(.true.)
      call fold_unpack_flat(mflx_n, nx, ny, nb)
      call fold_unpack_flat(mtflx_n, nx, ny, ntrace*nb)
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
