!! Seam hooks of the CICE-port sea-ice transport: the rank-uniform reductions and the
!! paired fluxes of the tripolar fold (roundabout-original, no upstream ancestry).
module rdb_ice_transport_seam
   !! Everything the incremental-remapping transport needs from the domain decomposition
   !! and the north fold that is NOT a halo of its input (those are the caller's: the category
   !! state, the face velocities and the metrics arrive with valid ghost rings).  The transport
   !! lives in the standalone ice library, which may not use MPI, so the collectives and the
   !! distributed fold exchange are PROCEDURE POINTERS the main library binds
   !! (`rdb_ice_transport_seam_comm`); this module holds the type, the dispatch and the one
   !! kernel that needs no communication.
   !!
   !! **Rank-uniform decisions.**  A step that stops (a CFL violation, a non-finite value, a
   !! negative area, a geometry error, a refused reconstruction input, a conservation failure)
   !! must stop on EVERY rank, or the ranks that go on enter a collective the others never
   !! reach.  Each decision point of the driver therefore sums its local counts over the ranks
   !! first (`seam_sum`) and decides on the global value, so after the call `res%n_*` and
   !! `res%ierr` are GLOBAL and identical on every rank; the cell error word `ice_err` stays
   !! per cell (the offender's own rank holds the bit).  The remap's clip record
   !! (`res%ledger`) is a plain sum over the OWNED cells of this rank (a diagnostic); the conservation ledger reads the
   !! per-cell planes (`fw%clip`) and reduces them with the EFP sums through `seam_sum_efp`, once per step.
   !! The configuration refusals (`TRANSPORT_CFG_*`: range, workspace, table) depend on the tile shape and the tracer table, so a decomposition whose tiles differ in
   !! them must be refused at configure; the land mask (`REMAP_CFG_MASK`) is checked on the
   !! read ring of each tile and IS reduced.
   !!
   !! **The fold seam.**  On a tripolar grid the north face of row `je` at column `i` and at
   !! its mirror column `i' = ni + 1 - i` are ONE face seen from its two sides: the area and
   !! the tracer content leaving cell `i` northward enter cell `i'`.  Each side computes its
   !! flux from its own (rotated) neighbourhood, so the two agree to round-off, not bitwise,
   !! and the update would create or destroy that difference.  After the flux integrals of a
   !! category batch and before the update the fold-row fluxes are PAIRED as the ocean pairs the
   !! duplicated v-face of the fold line: the east half is authoritative, the west half takes its
   !! negated mirror, and a self-conjugate column (odd `ni`) carries no flux.  Content
   !! then leaves one cell and enters the other in exactly the same bits, whatever the
   !! decomposition (`px = 1`: the local kernel here; `px > 1`: the owner-routed exchange of the
   !! ocean's fold, bound through `pair_fold_dist`).
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t
   implicit none
   private

   public :: transport_seam_t
   public :: seam_sum, seam_sum_efp, seam_pair_fold_flux, seam_sum_int, seam_sum_one, seam_or_bits
   public :: seam_sum_iface, seam_sum_efp_iface, seam_pair_fold_iface, seam_min_iface
   public :: seam_min_one, seam_all_true, seam_max_one

   integer, parameter :: ERR_WORD_BITS = 32
      !! Bits of the 32-bit error word.

   abstract interface
      subroutine seam_sum_iface(n, v)
         !! In-place sum of `n` integers over the ranks (collective).
         import :: int64
         implicit none
         integer, intent(in) :: n
         integer(int64), intent(inout) :: v(n)
      end subroutine seam_sum_iface

      subroutine seam_min_iface(v)
         !! In-place minimum of one integer over the ranks (collective).
         implicit none
         integer, intent(inout) :: v
      end subroutine seam_min_iface

      subroutine seam_sum_efp_iface(n, list)
         !! In-place, order-invariant sum of `n` EFP values over the ranks (collective).
         import :: efp_t
         implicit none
         integer, intent(in) :: n
         type(efp_t), intent(inout) :: list(n)
      end subroutine seam_sum_efp_iface

      subroutine seam_pair_fold_iface(nx, ny, nb, ntrace, with_tracers, mflx_n, mtflx_n)
         !! Pair the fold-row north fluxes across a split north rank row (collective over it).
         import :: wp
         implicit none
         integer, intent(in) :: nx, ny, nb, ntrace
         logical, intent(in) :: with_tracers
         real(wp), intent(inout) :: mflx_n(nx, ny, nb), mtflx_n(nx, ny, ntrace, nb)
      end subroutine seam_pair_fold_iface
   end interface

   type :: transport_seam_t
      !! What the transport's caller tells it about the decomposition (default: one rank, no fold).
      logical :: multirank = .false.
         !! More than one rank: the decisions are reduced through `sum_i8` / `sum_efp`.
      logical :: fold = .false.
         !! This rank holds the tripolar fold row (`bc%north_fold`, rank-local).
      logical :: fold_dist = .false.
         !! The fold row is split east-west (`px > 1`): `pair_fold_dist` does the pairing.
      logical :: pair_fold = .true.
         !! `.false.` leaves the fold-row fluxes UNPAIRED.  A diagnostic switch for the tests'
         !! negative control; the model never clears it.
      logical :: pair_tracers = .true.
         !! `.false.` pairs the area flux only (diagnostic switch of the tests' negative control).
      integer :: nbatch = 0
         !! Category indices per flux batch, the MINIMUM over the ranks (0: not bound, use the local one).  The
         !! batch loop and the fold pairing are collective per batch, so every rank must run the same batches.
      procedure(seam_sum_iface), pointer, nopass :: sum_i8 => null()
      procedure(seam_min_iface), pointer, nopass :: min_i => null()
      procedure(seam_sum_efp_iface), pointer, nopass :: sum_efp => null()
      procedure(seam_pair_fold_iface), pointer, nopass :: pair_fold_dist => null()
   end type transport_seam_t

contains

   subroutine seam_sum(seam, n, v)
      !! Global sum of `n` integers (no-op on one rank).  Collective on `multirank`.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(in) :: n
      integer(int64), intent(inout) :: v(n)

      if (.not. seam%multirank) return
      if (.not. associated(seam%sum_i8)) error stop "transport seam: multirank without a sum_i8 binding"
      call seam%sum_i8(n, v)
   end subroutine seam_sum

   subroutine seam_sum_int(seam, n, v)
      !! Global sum of `n` default integers (counts), in place.  No-op on one rank.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(in) :: n
      integer, intent(inout) :: v(n)
      integer(int64) :: w(n)

      if (.not. seam%multirank) return
      w = int(v, int64)
      call seam_sum(seam, n, w)
      v = int(w)
   end subroutine seam_sum_int

   subroutine seam_sum_one(seam, v)
      !! Global sum of one default integer (a count), in place.  No-op on one rank.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(inout) :: v
      integer :: w(1)

      if (.not. seam%multirank) return
      w(1) = v
      call seam_sum_int(seam, 1, w)
      v = w(1)
   end subroutine seam_sum_one

   subroutine seam_or_bits(seam, word)
      !! Global bitwise OR of an error word (each of its 32 bits summed over the ranks, set where
      !! any rank has it), in place.  No-op on one rank.
      type(transport_seam_t), intent(in) :: seam
      integer(int32), intent(inout) :: word
      integer(int64) :: cnt(0:ERR_WORD_BITS - 1)
      integer :: b

      if (.not. seam%multirank) return
      do b = 0, ERR_WORD_BITS - 1
         cnt(b) = merge(1_int64, 0_int64, btest(word, b))
      end do
      call seam_sum(seam, ERR_WORD_BITS, cnt)
      word = 0_int32
      do b = 0, ERR_WORD_BITS - 1
         if (cnt(b) > 0_int64) word = ibset(word, b)
      end do
   end subroutine seam_or_bits

   subroutine seam_min_one(seam, v)
      !! Global minimum of one integer, in place.  No-op on one rank.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(inout) :: v

      if (.not. seam%multirank) return
      if (.not. associated(seam%min_i)) error stop "transport seam: multirank without a min_i binding"
      call seam%min_i(v)
   end subroutine seam_min_one

   subroutine seam_max_one(seam, v)
      !! Global maximum of one integer, in place.  No-op on one rank.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(inout) :: v

      if (.not. seam%multirank) return
      v = -v
      call seam_min_one(seam, v)
      v = -v
   end subroutine seam_max_one

   function seam_all_true(seam, flag) result(all_ok)
      !! `flag` on EVERY rank (collective on `multirank`): the way a rank-local condition that gates a later
      !! collective is made rank-uniform.
      type(transport_seam_t), intent(in) :: seam
      logical, intent(in) :: flag
      logical :: all_ok
      integer :: v

      v = merge(1, 0, flag)
      call seam_min_one(seam, v)
      all_ok = v == 1
   end function seam_all_true

   subroutine seam_sum_efp(seam, n, list)
      !! Global EFP sum of `n` values (no-op on one rank).  Collective on `multirank`.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(in) :: n
      type(efp_t), intent(inout) :: list(n)

      if (.not. seam%multirank) return
      if (.not. associated(seam%sum_efp)) error stop "transport seam: multirank without a sum_efp binding"
      call seam%sum_efp(n, list)
   end subroutine seam_sum_efp

   subroutine seam_pair_fold_flux(seam, nx, ny, is, ie, je, nb, ntrace, mflx_n, mtflx_n)
      !! Pair the fold-row north fluxes of one category batch (module docstring).  A no-op
      !! on a rank that does not hold the fold row.  `is:ie` is the owned column range, `je`
      !! the last owned row; `mflx_n(:, je, :)` is the fold-line row.
      type(transport_seam_t), intent(in) :: seam
      integer, intent(in) :: nx, ny, is, ie, je, nb, ntrace
      real(wp), intent(inout) :: mflx_n(nx, ny, nb), mtflx_n(nx, ny, ntrace, nb)

      if (.not. seam%fold .or. .not. seam%pair_fold) return
      if (seam%fold_dist) then
         if (.not. associated(seam%pair_fold_dist)) error stop "transport seam: split fold without a binding"
         call seam%pair_fold_dist(nx, ny, nb, ntrace, seam%pair_tracers, mflx_n, mtflx_n)
      else
         call pair_fold_local(nx, ny, is, ie, je, nb, ntrace, seam%pair_tracers, mflx_n, mtflx_n)
      end if
   end subroutine seam_pair_fold_flux

   subroutine pair_fold_local(nx, ny, is, ie, je, nb, ntrace, with_tracers, mflx_n, mtflx_n)
      !! `px = 1`: the whole fold row is local.  Owned column `p` (1-based) is in the west half
      !! when `p < ni + 1 - p`: it takes the negated flux of its mirror, the east half is read
      !! only; the self-conjugate column carries nothing.  Two launches (area, tracers).
      integer, intent(in) :: nx, ny, is, ie, je, nb, ntrace
      logical, intent(in) :: with_tracers
      real(wp), intent(inout) :: mflx_n(nx, ny, nb), mtflx_n(nx, ny, ntrace, nb)
      integer :: i, b, t, p, pm, ni

      ni = ie - is + 1
      !$acc data present(mflx_n, mtflx_n)
      do concurrent(b=1:nb, i=is:ie) local(p, pm)
         p = i - is + 1
         pm = ni + 1 - p
         if (p < pm) then
            mflx_n(i, je, b) = -mflx_n(is - 1 + pm, je, b)
         else if (p == pm) then
            mflx_n(i, je, b) = 0.0_wp
         end if
      end do
      if (with_tracers) then
      do concurrent(b=1:nb, t=1:ntrace, i=is:ie) local(p, pm)
         p = i - is + 1
         pm = ni + 1 - p
         if (p < pm) then
            mtflx_n(i, je, t, b) = -mtflx_n(is - 1 + pm, je, t, b)
         else if (p == pm) then
            mtflx_n(i, je, t, b) = 0.0_wp
         end if
      end do
      end if
      !$acc end data
   end subroutine pair_fold_local

end module rdb_ice_transport_seam
