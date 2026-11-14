module rdb_ice_errors
   !! Device-side error registry of the sea-ice column physics.
   !!
   !! Every ice column routine owns one `int32` word per cell (`ice_err(nx,ny)`,
   !! zeroed at the start of each ice step) and ORs a bit into it where a
   !! reference implementation would abort or warn.  A word is written only by
   !! its own cell's work item, so the kernels stay race-free, and nothing
   !! aborts inside a kernel.  After the step's kernels ONE host call to
   !! `ice_err_reduce` folds the grid into an `ice_err_summary_t`: the OR of all
   !! words, a per-bit cell count, and the first (lowest linear index) cell that
   !! carries a FATAL bit.  The host then either aborts through `ice_err_enforce`
   !! (FATAL bit present: first-offender dump, `error stop`) or prints the
   !! COUNTED bits on the console line (`ice_err_counted_line`).
   !!
   !! Bits are `ishft(1, n)` for `n = 0 .. ICE_NERR-1`; the registry below is the
   !! single authority.  Class is FATAL (F) or COUNTED (C); bit 10
   !! (`ICE_ERR_ITD_HBNEW_ORDER`) is COUNTED unless `&ocean_ice_nml
   !! fatal_on_error`, which `ice_err_fatal_set` folds in.
   !!
   !! Reduction layout.  The grid is reduced in passes no wider than a few
   !! accumulators each (a reduction kernel carrying many accumulators has been seen
   !! to corrupt on the GPU build): pass 1 carries two accumulators
   !! (`ior` mask, `min` first-offender index); the per-bit counters run
   !! `ICE_ERR_PASS_WIDTH` bits per pass and only for chunks whose bits appear
   !! in the mask, so a clean grid costs exactly one pass.
   !!
   !! Multi-rank: `ice_err_reduce` ends in `combine_across_ranks`, today a
   !! no-op.  A later PR fills that one private routine with a collective
   !! (bitwise-OR of `mask`, sum of `count`, min of the GLOBAL first-offender
   !! index) through `pic_mpi_lib`.  The reduce already takes the owned range and
   !! global offsets, and every caller consumes the summary, never the raw grid, so
   !! the only caller-visible change is a communicator argument (or a module-level
   !! one set at setup); see `combine_across_ranks` for the full contract.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use pic_logger, only: global_logger
   implicit none
   private

   public :: ICE_NERR, ICE_ERR_PASS_WIDTH
   public :: ICE_ERR_NONE
   public :: ICE_ERR_TSFC_NOCONV, ICE_ERR_TSN_RANGE, ICE_ERR_SIN_LOW, ICE_ERR_TIN_RANGE
   public :: ICE_ERR_QIN_SQRT, ICE_ERR_NONFINITE, ICE_ERR_THERMO_CONSERV, ICE_ERR_TBOT
   public :: ICE_ERR_NEG_THICK, ICE_ERR_NEW_ICE_AREA, ICE_ERR_ITD_HBNEW_ORDER
   public :: ICE_ERR_ITD_SHIFT_RANGE, ICE_ERR_AICE_RANGE, ICE_ERR_AICEN_NEG
   public :: ICE_ERR_RIDGE_NOCONV, ICE_ERR_RIDGE_AICE0, ICE_ERR_RIDGE_OVERSHOOT
   public :: ICE_ERR_RIDGE_DEGENERATE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_NEG_AREA
   public :: ICE_ERR_TRANSPORT_CONSERV, ICE_ERR_DYN_VEL, ICE_ERR_REMAP_FALLBACK
   public :: ICE_ERR_RIDGE_PASSES, ICE_ERR_FLUX_DROPPED
   public :: ICE_ERR_FATAL_MASK, ICE_ERR_COUNTED_MASK
   public :: ice_err_summary_t
   public :: ice_err_name, ice_err_fatal_set, ice_err_clear, ice_err_reduce
   public :: ice_err_is_fatal, ice_err_dump, ice_err_enforce
   public :: ice_err_counted_line, ice_err_accumulate

   integer, parameter :: ICE_NERR = 25
      !! Number of registered bits (0 .. ICE_NERR-1).
   integer, parameter :: ICE_ERR_PASS_WIDTH = 4
      !! Per-bit counters carried by one reduction kernel.

   integer(int32), parameter :: ICE_ERR_NONE = 0_int32
   integer(int32), parameter :: ICE_ERR_TSFC_NOCONV = ishft(1_int32, 0)
      !! F. Surface-temperature iteration did not converge.
   integer(int32), parameter :: ICE_ERR_TSN_RANGE = ishft(1_int32, 1)
      !! F. Snow-layer temperature outside [Tmin, Tmax].
   integer(int32), parameter :: ICE_ERR_SIN_LOW = ishft(1_int32, 2)
      !! F. Ice salinity below the minimum.
   integer(int32), parameter :: ICE_ERR_TIN_RANGE = ishft(1_int32, 3)
      !! F. Ice-layer temperature above melting or below Tmin.
   integer(int32), parameter :: ICE_ERR_QIN_SQRT = ishft(1_int32, 4)
      !! F. Negative discriminant in the enthalpy-to-temperature inversion.
   integer(int32), parameter :: ICE_ERR_NONFINITE = ishft(1_int32, 5)
      !! F. A kernel's entry/exit finiteness guard tripped.
   integer(int32), parameter :: ICE_ERR_THERMO_CONSERV = ishft(1_int32, 6)
      !! F. Vertical-thermodynamics conservation error above tolerance.
   integer(int32), parameter :: ICE_ERR_TBOT = ishft(1_int32, 7)
      !! F. Bottom temperature above the freezing point.
   integer(int32), parameter :: ICE_ERR_NEG_THICK = ishft(1_int32, 8)
      !! F. Negative layer thickness.
   integer(int32), parameter :: ICE_ERR_NEW_ICE_AREA = ishft(1_int32, 9)
      !! F. New-ice step: a category area below `-puny` at entry, or a negative open-water area at exit.
   integer(int32), parameter :: ICE_ERR_ITD_HBNEW_ORDER = ishft(1_int32, 10)
      !! C (F with `fatal_on_error`). Category-boundary order violated; the cell falls back to a rebin.
   integer(int32), parameter :: ICE_ERR_ITD_SHIFT_RANGE = ishft(1_int32, 11)
      !! F. Category-shift area/volume change outside its tolerance band.
   integer(int32), parameter :: ICE_ERR_AICE_RANGE = ishft(1_int32, 12)
      !! F. Total concentration outside [-puny, 1+puny].
   integer(int32), parameter :: ICE_ERR_AICEN_NEG = ishft(1_int32, 13)
      !! F. Negative category concentration beyond tolerance.
   integer(int32), parameter :: ICE_ERR_RIDGE_NOCONV = ishft(1_int32, 14)
      !! F. Ridging iteration did not converge.
   integer(int32), parameter :: ICE_ERR_RIDGE_AICE0 = ishft(1_int32, 15)
      !! F. Ridging left negative open-water area.
   integer(int32), parameter :: ICE_ERR_RIDGE_OVERSHOOT = ishft(1_int32, 16)
      !! F. Ridged area exceeds the category area.
   integer(int32), parameter :: ICE_ERR_RIDGE_DEGENERATE = ishft(1_int32, 17)
      !! F. Degenerate ridge-participation sum under closing.
   integer(int32), parameter :: ICE_ERR_TRANSPORT_CFL = ishft(1_int32, 18)
      !! F. Departure-point CFL violated.
   integer(int32), parameter :: ICE_ERR_TRANSPORT_NEG_AREA = ishft(1_int32, 19)
      !! F. Negative mass after a transport update.
   integer(int32), parameter :: ICE_ERR_TRANSPORT_CONSERV = ishft(1_int32, 20)
      !! F. Transport conservation check off tolerance.
   integer(int32), parameter :: ICE_ERR_DYN_VEL = ishft(1_int32, 21)
      !! F. Non-finite or over-cap face velocity after the dynamics solve.
   integer(int32), parameter :: ICE_ERR_REMAP_FALLBACK = ishft(1_int32, 22)
      !! C. Cell took the non-remapping fallback.
   integer(int32), parameter :: ICE_ERR_RIDGE_PASSES = ishft(1_int32, 23)
      !! C. Ridging needed more than one pass.
   integer(int32), parameter :: ICE_ERR_FLUX_DROPPED = ishft(1_int32, 24)
      !! C. Ice flux produced in a step that ended with no ice (ledger only).

   integer(int32), parameter :: ICE_ERR_COUNTED_MASK = ior(ior(ICE_ERR_ITD_HBNEW_ORDER, ICE_ERR_REMAP_FALLBACK), &
                                                           ior(ICE_ERR_RIDGE_PASSES, ICE_ERR_FLUX_DROPPED))
      !! Bits that are COUNTED by default (bit 10 included; see `ice_err_fatal_set`).
   integer(int32), parameter :: ICE_ERR_FATAL_MASK = &
                                iand(not(ICE_ERR_COUNTED_MASK), ishft(1_int32, ICE_NERR) - 1_int32)
      !! Bits that are always FATAL (every registered bit not in the counted set).

   type :: ice_err_summary_t
      !! One reduction of an `ice_err` grid.
      integer(int32) :: mask = 0_int32
         !! OR of every cell's word.
      integer(int32) :: count(0:ICE_NERR - 1) = 0_int32
         !! Number of cells carrying each bit (filled only for bits present in `mask`).
      integer(int64) :: first_idx = 0_int64
         !! Lowest GLOBAL linear index `(jg-1)*nx_global + ig` of an owned cell with a
         !! FATAL bit (`ig = (i - is + 1) + i_offset_global`, likewise `jg`; the offsets are the global index of the first owned cell minus 1); 0 when none.
      integer(int32) :: first_word = 0_int32
         !! The error word of that cell.
   end type ice_err_summary_t

contains

   pure function ice_err_name(ibit) result(name)
      !! Registry name of bit `ibit` (0-based), or `UNKNOWN`.
      integer, intent(in) :: ibit
      character(len=26) :: name
      character(len=26), parameter :: table(0:ICE_NERR - 1) = [character(len=26) :: &
                                                               "TSFC_NOCONV", &
                                                               "TSN_RANGE", &
                                                               "SIN_LOW", &
                                                               "TIN_RANGE", &
                                                               "QIN_SQRT", &
                                                               "NONFINITE", &
                                                               "THERMO_CONSERV", &
                                                               "TBOT", &
                                                               "NEG_THICK", &
                                                               "NEW_ICE_AREA", &
                                                               "ITD_HBNEW_ORDER", &
                                                               "ITD_SHIFT_RANGE", &
                                                               "AICE_RANGE", &
                                                               "AICEN_NEG", &
                                                               "RIDGE_NOCONV", &
                                                               "RIDGE_AICE0", &
                                                               "RIDGE_OVERSHOOT", &
                                                               "RIDGE_DEGENERATE", &
                                                               "TRANSPORT_CFL", &
                                                               "TRANSPORT_NEG_AREA", &
                                                               "TRANSPORT_CONSERV", &
                                                               "DYN_VEL", &
                                                               "REMAP_FALLBACK", &
                                                               "RIDGE_PASSES", &
                                                               "FLUX_DROPPED"]
      name = "UNKNOWN"
      if (ibit >= 0 .and. ibit < ICE_NERR) name = table(ibit)
   end function ice_err_name

   pure function ice_err_fatal_set(fatal_on_error) result(m)
      !! The set of bits that abort the run: always `ICE_ERR_FATAL_MASK`, plus
      !! the ITD boundary-order bit when `fatal_on_error`.
      logical, intent(in) :: fatal_on_error
      integer(int32) :: m
      m = ICE_ERR_FATAL_MASK
      if (fatal_on_error) m = ior(m, ICE_ERR_ITD_HBNEW_ORDER)
   end function ice_err_fatal_set

   pure subroutine ice_err_clear(ice_err, nx, ny)
      !! Zero the per-cell error words (start of every ice step).
      integer, intent(in) :: nx, ny
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         ice_err(i, j) = 0_int32
      end do
   end subroutine ice_err_clear

   subroutine ice_err_reduce(ice_err, nx, ny, is, ie, js, je, i_offset_global, j_offset_global, nx_global, &
                             fatal_mask, summary)
      !! Reduce the OWNED cells `is:ie, js:je` of the error grid into `summary`
      !! (host-resident result).  Ghost cells are excluded: under periodic or fold
      !! halos they duplicate owned cells, so counting them would double-count and
      !! could name a ghost as the first offender.  `first_idx` is the GLOBAL
      !! linear index, so a cross-rank MINLOC of it is meaningful.  `ice_err` must
      !! be device-present on an offload build (the kernels assert it).
      integer, intent(in) :: nx, ny, is, ie, js, je
      integer, intent(in) :: i_offset_global, j_offset_global, nx_global
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int32), intent(in) :: fatal_mask
      type(ice_err_summary_t), intent(out) :: summary
      integer :: b0, q
      integer(int32) :: c4(ICE_ERR_PASS_WIDTH)
      integer(int32) :: chunk

      call reduce_mask_first(ice_err, nx, ny, is, ie, js, je, i_offset_global, j_offset_global, &
                             nx_global, fatal_mask, summary%mask, summary%first_idx)
      summary%count = 0_int32
      do b0 = 0, ICE_NERR - 1, ICE_ERR_PASS_WIDTH
         chunk = 0_int32
         do q = b0, min(b0 + ICE_ERR_PASS_WIDTH - 1, ICE_NERR - 1)
            chunk = ior(chunk, ishft(1_int32, q))
         end do
         if (iand(summary%mask, chunk) == 0_int32) cycle
         call reduce_count4(ice_err, nx, ny, is, ie, js, je, b0, c4)
         do q = 1, ICE_ERR_PASS_WIDTH
            if (b0 + q - 1 <= ICE_NERR - 1) summary%count(b0 + q - 1) = c4(q)
         end do
      end do
      summary%first_word = 0_int32
      if (summary%first_idx > 0_int64) then
         call fetch_word(ice_err, nx, ny, is, js, i_offset_global, j_offset_global, nx_global, &
                         summary%first_idx, summary%first_word)
      end if
      call combine_across_ranks(summary)
   end subroutine ice_err_reduce

   subroutine reduce_mask_first(ice_err, nx, ny, is, ie, js, je, i_off, j_off, nxg, fatal_mask, mask, first_idx)
      !! Pass 1: OR of owned words and lowest fatal GLOBAL linear index (2 accumulators).
      integer, intent(in) :: nx, ny, is, ie, js, je, i_off, j_off, nxg
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int32), intent(in) :: fatal_mask
      integer(int32), intent(out) :: mask
      integer(int64), intent(out) :: first_idx
      integer(int32) :: m, e
      integer(int64) :: f
      integer :: i, j
      m = 0_int32
      f = huge(0_int64)
      ! `ice_err` must be device-mapped by the caller.
      do concurrent(j=js:je, i=is:ie) local(e) reduce(ior:m) reduce(min:f)
         e = ice_err(i, j)
         m = ior(m, e)
         if (iand(e, fatal_mask) /= 0_int32) f = min(f, int(j - js + j_off, int64)*nxg + int(i - is + 1 + i_off, int64))
      end do
      mask = m
      first_idx = f
      if (first_idx == huge(0_int64)) first_idx = 0_int64
   end subroutine reduce_mask_first

   subroutine reduce_count4(ice_err, nx, ny, is, ie, js, je, b0, c4)
      !! Counter pass: owned cells carrying each of bits `b0 .. b0+3` (4 accumulators).
      integer, intent(in) :: nx, ny, is, ie, js, je, b0
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int32), intent(out) :: c4(ICE_ERR_PASS_WIDTH)
      integer(int32) :: c1, c2, c3, c4s, e
      integer :: i, j
      c1 = 0_int32
      c2 = 0_int32
      c3 = 0_int32
      c4s = 0_int32
      ! `ice_err` must be device-mapped by the caller.
      do concurrent(j=js:je, i=is:ie) local(e) reduce(+:c1, c2, c3, c4s)
         e = ice_err(i, j)
         if (btest(e, b0)) c1 = c1 + 1_int32
         if (btest(e, b0 + 1)) c2 = c2 + 1_int32
         if (btest(e, b0 + 2)) c3 = c3 + 1_int32
         if (btest(e, b0 + 3)) c4s = c4s + 1_int32
      end do
      c4 = [c1, c2, c3, c4s]
   end subroutine reduce_count4

   subroutine fetch_word(ice_err, nx, ny, is, js, i_off, j_off, nxg, idx, word)
      !! The word of the cell with global linear index `idx`: a one-element reduction
      !! (no whole-grid pass; `ice_err` stays `intent(in)`, no host update of it).
      integer, intent(in) :: nx, ny, is, js, i_off, j_off, nxg
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int64), intent(in) :: idx
      integer(int32), intent(out) :: word
      integer(int32) :: w
      integer :: i, j, k, ig, jg
      jg = int((idx - 1_int64)/int(nxg, int64)) + 1
      ig = int(idx - int(jg - 1, int64)*int(nxg, int64))
      j = jg - j_off - 1 + js
      i = ig - i_off - 1 + is
      w = 0_int32
      do concurrent(k=1:1) reduce(ior:w)
         w = ior(w, ice_err(i, j))
      end do
      word = w
   end subroutine fetch_word

   subroutine combine_across_ranks(summary)
      !! Multi-rank hook, a no-op on one rank.  The later MPI PR (through
      !! `pic_mpi_lib`; no MPI import here) must make this the collective:
      !!   1. BOR-allreduce `summary%mask`;
      !!   2. SUM-allreduce `summary%count`;
      !!   3. MINLOC-allreduce `summary%first_idx` (the global index, 0 = none
      !!      mapped to `huge` for the min), then broadcast the OWNING rank's
      !!      `first_word` so every rank holds the same offender.
      !! Callers must then (a) call `ice_err_enforce` UNCONDITIONALLY on every
      !! rank after this collective so the abort verdict is rank-uniform, and
      !! (b) `ice_err_dump` decodes the global index with the GLOBAL `nx` and
      !! logs on one rank only.  Signatures stay unchanged except that
      !! `ice_err_reduce` gains a communicator argument or reads a module-level
      !! one set at setup; that is the only caller-visible change.
      type(ice_err_summary_t), intent(inout) :: summary
      ! Single rank: nothing to combine; `summary` is already global (the dummy is
      ! deliberately unused until the collective lands).
   end subroutine combine_across_ranks

   pure function ice_err_is_fatal(summary, fatal_mask) result(fatal)
      !! True when any bit of `fatal_mask` is present in the reduced mask.
      type(ice_err_summary_t), intent(in) :: summary
      integer(int32), intent(in) :: fatal_mask
      logical :: fatal
      fatal = iand(summary%mask, fatal_mask) /= 0_int32
   end function ice_err_is_fatal

   subroutine ice_err_dump(summary, nx, fatal_mask)
      !! Log every present bit (name, class, cell count) and the first-offender cell.
      !! `nx` is the GLOBAL x extent (decodes the global first-offender index).
      type(ice_err_summary_t), intent(in) :: summary
      integer, intent(in) :: nx
      integer(int32), intent(in) :: fatal_mask
      character(len=256) :: msg
      character(len=1) :: cls
      integer :: ib, i, j
      do ib = 0, ICE_NERR - 1
         if (.not. btest(summary%mask, ib)) cycle
         cls = "C"
         if (btest(fatal_mask, ib)) cls = "F"
         write (msg, "(a,a,a,a,a,i0,a)") "ice error bit ", trim(ice_err_name(ib)), " [", cls, "] in ", &
            summary%count(ib), " cell(s)"
         call global_logger%error(trim(msg))
      end do
      if (summary%first_idx > 0_int64) then
         j = int((summary%first_idx - 1_int64)/int(nx, int64)) + 1
         i = int(summary%first_idx - int(j - 1, int64)*int(nx, int64))
         write (msg, "(a,i0,a,i0,a,z8.8)") "ice first offender: global cell (i,j) = (", i, ",", j, &
            "), error word 0x", summary%first_word
         call global_logger%error(trim(msg))
      end if
   end subroutine ice_err_dump

   subroutine ice_err_enforce(summary, nx, fatal_mask)
      !! Abort loudly (dump, then `error stop`) when a FATAL bit is present.
      type(ice_err_summary_t), intent(in) :: summary
      integer, intent(in) :: nx
      integer(int32), intent(in) :: fatal_mask
      if (.not. ice_err_is_fatal(summary, fatal_mask)) return
      call ice_err_dump(summary, nx, fatal_mask)
      error stop "rdb_ice_errors: fatal sea-ice error bit set"
   end subroutine ice_err_enforce

   pure function ice_err_counted_line(total, fatal_mask) result(line)
      !! Console fragment for the COUNTED bits with a non-zero run total
      !! (`total` from `ice_err_accumulate`); empty when none.
      integer(int64), intent(in) :: total(0:ICE_NERR - 1)
      integer(int32), intent(in) :: fatal_mask
      character(len=:), allocatable :: line
      character(len=48) :: piece
      integer :: ib
      line = ""
      do ib = 0, ICE_NERR - 1
         if (btest(fatal_mask, ib) .or. total(ib) == 0_int64) cycle
         write (piece, "(a,a,i0)") trim(ice_err_name(ib)), "=", total(ib)
         if (len(line) > 0) line = line//" "
         line = line//trim(piece)
      end do
   end function ice_err_counted_line

   pure subroutine ice_err_accumulate(total, summary)
      !! Add one step's per-bit cell counts into the run totals.
      integer(int64), intent(inout) :: total(0:ICE_NERR - 1)
      type(ice_err_summary_t), intent(in) :: summary
      total = total + int(summary%count, int64)
   end subroutine ice_err_accumulate

end module rdb_ice_errors
