!! Ring buffer of specific error messages for the ocean solver-creation
!! path (Python runtime API plan, P2, `docs/ocean_python_api_plan.md`).
module rdb_error_ring
   !! P0/P0.1 gave every setup procedure an `optional, intent(out) :: ierr`
   !! returning a coarse STAGE code (`rdb_ocean_status`) while the specific
   !! reason went only to `global_logger%error` — reachable from a Fortran
   !! log, not from a C-ABI caller. `RuntimeError: error 3` is exactly the
   !! failure this module exists to prevent.
   !!
   !! A RING, not a single last-error slot: the P0 review (finding F2) found
   !! that at ~23 sites an outer `configure_ocean_*` wrapper calls an inner
   !! procedure with `ierr=` unconditionally, so the inner SPECIFIC sentence
   !! (e.g. a NetCDF dimension-mismatch message) is immediately followed by
   !! an outer GENERIC one. A single slot keeps only the last write and loses
   !! the specific sentence; the ring keeps both, in push order, so a caller
   !! renders index 0 (most recent = deepest/most specific push) as the
   !! headline message.
   !!
   !! `fail(msg, ierr, code)` (below) is the single entry point that replaces
   !! the hand-written `logger%error(msg)` + `if (present(ierr)) then; ierr =
   !! code; return; else; error stop msg; end if` triple that P0 scattered
   !! across ~99 call sites. It ALWAYS pushes `msg` to the ring and logs it
   !! (unchanged behaviour for anyone still reading the buffered log), then:
   !! `ierr` present -> sets `ierr = code` and returns to the caller (which
   !! must itself `return` immediately — this routine cannot return from ITS
   !! caller); `ierr` absent -> `error stop msg`, so the abort text and the
   !! ring both carry the SAME specific message (previously these could
   !! diverge: the `error stop` text was often a shorter, separately
   !! hand-written string).
   !!
   !! Cleared at the top of every `bind(c)` entry point in `rdb_ocean_api`
   !! that begins a fresh operation, so a stale message from a previous call
   !! never masquerades as the cause of a later, unrelated one.
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: fail
   public :: error_ring_push
   public :: error_ring_get
   public :: error_ring_clear
   public :: error_ring_count

   integer, parameter, public :: ERROR_RING_SLOTS = 16
      !! Comfortably deeper than the deepest `ierr=` threading chain on the
      !! create() path (a handful of frames at most).
   integer, parameter, public :: ERROR_RING_MSG_LEN = 512
      !! Matches the C-ABI buffer contract documented at
      !! `rdb_ocean_last_error` (`include/rdb_ocean.h`).

   character(len=ERROR_RING_MSG_LEN), save :: ring(ERROR_RING_SLOTS) = ""
   integer, save :: ring_count = 0
      !! Number of live slots, saturating at ERROR_RING_SLOTS.
   integer, save :: ring_head = 0
      !! 1-based index of the MOST RECENT push; wraps modulo ERROR_RING_SLOTS.

contains

   subroutine error_ring_clear()
      !! Reset the ring to empty.
      ring = ""
      ring_count = 0
      ring_head = 0
   end subroutine error_ring_clear

   subroutine error_ring_push(msg)
      !! Push one message; `error_ring_get(0)` is always the most recent.
      !! Overwrites the oldest slot once full.
      character(len=*), intent(in) :: msg

      ring_head = mod(ring_head, ERROR_RING_SLOTS) + 1
      ring(ring_head) = msg
      if (ring_count < ERROR_RING_SLOTS) ring_count = ring_count + 1
   end subroutine error_ring_push

   pure function error_ring_count() result(n)
      !! Number of live messages in the ring (0..ERROR_RING_SLOTS).
      integer :: n
      n = ring_count
   end function error_ring_count

   pure function error_ring_get(i) result(msg)
      !! Message at ring-relative index `i` (0 = most recent, 1 =
      !! next-most-recent, ...). Empty string if `i` is out of
      !! `[0, error_ring_count())`.
      integer, intent(in) :: i
      character(len=ERROR_RING_MSG_LEN) :: msg

      integer :: slot

      if (i < 0 .or. i >= ring_count) then
         msg = ""
         return
      end if
      slot = ring_head - i
      if (slot < 1) slot = slot + ERROR_RING_SLOTS
      msg = ring(slot)
   end function error_ring_get

   subroutine fail(msg, ierr, code)
      !! Push `msg` to the error ring, log it (unchanged behaviour), then
      !! either set `ierr = code` and return (caller must `return`
      !! immediately after this call) or `error stop msg` when `ierr` is
      !! absent. See the module docstring for the full contract.
      character(len=*), intent(in) :: msg
      integer, intent(out), optional :: ierr
      integer, intent(in) :: code

      call error_ring_push(msg)
      call logger%error(msg)
      if (present(ierr)) then
         ierr = code
      else
         error stop msg
      end if
   end subroutine fail

end module rdb_error_ring
