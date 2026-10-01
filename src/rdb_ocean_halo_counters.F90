!! Ocean halo-exchange counters — backend-agnostic bookkeeping module.
!!
!! Provides deterministic per-category exchange counters that are
!! incremented at SEMANTIC call sites (NOT deep inside primitives) to
!! give a load-independent gate for halo-exchange profiling and testing.
!!
!! DOUBLE-COUNT PREVENTION
!! =======================
!! The 3d wrappers (`centre_3d`, `face_x_3d`, `face_y_3d`) now perform
!! batched single-message exchanges (no per-layer 2d loop).  The grouped
!! `bt_group_2d` and `ocean_halo_exchange_ml_state` wrappers call the
!! primitives internally.
!! Suppression is NESTING-SAFE via a depth counter (`oh_suppress_depth`):
!! each call to `oh_count_suppress_on()` increments the depth, and each
!! `oh_count_suppress_off()` decrements it (floored at 0).  Increment
!! routines check `oh_suppress_depth > 0` before bumping their counter.
!! This handles multi-level nesting: ml_state suppress_on → outer depth
!! still positive, so inner 3d counters stay suppressed inside ml_state.
!!
!! CATEGORIES (9 semantic counters, disjoint and additive):
!!   oh_cnt_bt_group    — ocean_halo_bt_group_2d grouped exchange
!!   oh_cnt_bt_u_mid    — mid-substep bt_ubt face_x exchange (BT Pass-2b→2c)
!!   oh_cnt_ml_state    — ocean_halo_exchange_ml_state whole-ML-state exchange
!!   oh_cnt_centre_2d   — standalone centre_2d primitive
!!   oh_cnt_centre_3d   — standalone centre_3d primitive (batched)
!!   oh_cnt_face_x_2d   — standalone face_x_2d primitive
!!   oh_cnt_face_x_3d   — standalone face_x_3d primitive (batched)
!!   oh_cnt_face_y_2d   — standalone face_y_2d primitive
!!   oh_cnt_face_y_3d   — standalone face_y_3d primitive (batched)
!!
!! MESSAGE COUNTER (oh_cnt_msgs):
!!   Total MPI ISENDs posted by the ocean halo MPI backend.  NOT subject
!!   to the suppress depth — messages are physical MPI operations regardless
!!   of semantic grouping, so suppression does not apply.  Incremented in
!!   `rdb_ocean_halo` on the decomposed axes only; an undecomposed axis posts
!!   no messages, so this stays zero in non-MPI builds.
!!
!! All counters are host-side scalars (never inside a do concurrent /
!! OpenACC region); they are pure integer bookkeeping.
!! No MPI dependency — compiled in BOTH MPI and single-process configs.
module rdb_ocean_halo_counters
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_strings, only: to_string
   implicit none
   private

   public :: oh_count_bt_group
   public :: oh_count_bt_u_mid
   public :: oh_count_ml_state
   public :: oh_count_centre_2d
   public :: oh_count_centre_3d
   public :: oh_count_face_x_2d
   public :: oh_count_face_x_3d
   public :: oh_count_face_y_2d
   public :: oh_count_face_y_3d
   public :: oh_count_msgs
   public :: oh_count_suppress_on
   public :: oh_count_suppress_off
   public :: oh_count_suppressed
   public :: oh_counters_reset
   public :: oh_counters_get
   public :: oh_counters_msgs
   public :: oh_counters_total
   public :: oh_counters_format

   ! ------------------------------------------------------------------
   ! Module-level counters (all initialised 0 at program start).
   ! ------------------------------------------------------------------
   integer(int64), save :: oh_cnt_bt_group = 0_int64
      !! ocean_halo_bt_group_2d grouped exchange calls
   integer(int64), save :: oh_cnt_bt_u_mid = 0_int64
      !! Mid-substep bt_ubt face_x exchange calls (BT Pass-2b→2c u seam)
   integer(int64), save :: oh_cnt_ml_state = 0_int64
      !! ocean_halo_exchange_ml_state whole-ML-state exchange calls
   integer(int64), save :: oh_cnt_centre_2d = 0_int64
      !! Standalone centre_2d primitive calls
   integer(int64), save :: oh_cnt_centre_3d = 0_int64
      !! Standalone centre_3d primitive calls
   integer(int64), save :: oh_cnt_face_x_2d = 0_int64
      !! Standalone face_x_2d primitive calls
   integer(int64), save :: oh_cnt_face_x_3d = 0_int64
      !! Standalone face_x_3d primitive calls
   integer(int64), save :: oh_cnt_face_y_2d = 0_int64
      !! Standalone face_y_2d primitive calls
   integer(int64), save :: oh_cnt_face_y_3d = 0_int64
      !! Standalone face_y_3d primitive calls
   integer(int64), save :: oh_cnt_msgs = 0_int64
      !! Total MPI ISENDs posted by the ocean halo MPI backend.
      !! NOT subject to the suppress depth — messages are physical MPI
      !! operations regardless of semantic grouping, so suppression does
      !! not apply.  Zero on single-process builds (no messages posted).

   ! ------------------------------------------------------------------
   ! Suppress depth counter — incremented by oh_count_suppress_on() and
   ! decremented (floored at 0) by oh_count_suppress_off().  Nesting-safe:
   ! ml_state wraps 3d wrappers; the outer depth stays positive until the
   ! outermost suppress_off(), ensuring no inner 3d counter fires while
   ! inside a grouped wrapper.  oh_cnt_msgs is NOT suppressed.
   ! ------------------------------------------------------------------
   integer, save :: oh_suppress_depth = 0
      !! When > 0 all oh_count_* increment calls (except oh_count_msgs) are no-ops.

contains

   ! ------------------------------------------------------------------
   ! Suppress control
   ! ------------------------------------------------------------------

   subroutine oh_count_suppress_on()
      !! Increment suppress depth so inner primitive increments are no-ops.
      !! Nesting-safe: each on/off pair increments/decrements a counter so
      !! nested suppress pairs (ml_state wrapping 3d wrappers) work correctly.
      oh_suppress_depth = oh_suppress_depth + 1
   end subroutine oh_count_suppress_on

   subroutine oh_count_suppress_off()
      !! Decrement suppress depth (floored at 0) after a wrapper completes.
      oh_suppress_depth = max(oh_suppress_depth - 1, 0)
   end subroutine oh_count_suppress_off

   function oh_count_suppressed() result(flag)
      !! Return .true. if suppress depth > 0 (any suppress scope is active).
      logical :: flag
      flag = (oh_suppress_depth > 0)
   end function oh_count_suppressed

   ! ------------------------------------------------------------------
   ! Increment routines (each checks suppress before bumping).
   ! ------------------------------------------------------------------

   subroutine oh_count_bt_group()
      !! Increment the bt_group counter.
      !! Called once in ocean_halo_bt_group_2d bodies.
      if (oh_suppress_depth > 0) return
      oh_cnt_bt_group = oh_cnt_bt_group + 1_int64
   end subroutine oh_count_bt_group

   subroutine oh_count_bt_u_mid()
      !! Increment the bt_u_mid counter.
      !! Called at the mid-substep bt_ubt face_x exchange site.
      if (oh_suppress_depth > 0) return
      oh_cnt_bt_u_mid = oh_cnt_bt_u_mid + 1_int64
   end subroutine oh_count_bt_u_mid

   subroutine oh_count_ml_state()
      !! Increment the ml_state counter.
      !! Called once in ocean_halo_exchange_ml_state.
      if (oh_suppress_depth > 0) return
      oh_cnt_ml_state = oh_cnt_ml_state + 1_int64
   end subroutine oh_count_ml_state

   subroutine oh_count_centre_2d()
      !! Increment the centre_2d standalone counter.
      !! Called at the top of each ocean_halo_centre_2d body.
      if (oh_suppress_depth > 0) return
      oh_cnt_centre_2d = oh_cnt_centre_2d + 1_int64
   end subroutine oh_count_centre_2d

   subroutine oh_count_centre_3d()
      !! Increment the centre_3d standalone counter.
      !! Called once in ocean_halo_centre_3d (after suppressing inner 2d bumps).
      if (oh_suppress_depth > 0) return
      oh_cnt_centre_3d = oh_cnt_centre_3d + 1_int64
   end subroutine oh_count_centre_3d

   subroutine oh_count_face_x_2d()
      !! Increment the face_x_2d standalone counter.
      !! Called at the top of each ocean_halo_face_x_2d body.
      if (oh_suppress_depth > 0) return
      oh_cnt_face_x_2d = oh_cnt_face_x_2d + 1_int64
   end subroutine oh_count_face_x_2d

   subroutine oh_count_face_x_3d()
      !! Increment the face_x_3d standalone counter.
      !! Called once in ocean_halo_face_x_3d (after suppressing inner 2d bumps).
      if (oh_suppress_depth > 0) return
      oh_cnt_face_x_3d = oh_cnt_face_x_3d + 1_int64
   end subroutine oh_count_face_x_3d

   subroutine oh_count_face_y_2d()
      !! Increment the face_y_2d standalone counter.
      !! Called at the top of each ocean_halo_face_y_2d body.
      if (oh_suppress_depth > 0) return
      oh_cnt_face_y_2d = oh_cnt_face_y_2d + 1_int64
   end subroutine oh_count_face_y_2d

   subroutine oh_count_face_y_3d()
      !! Increment the face_y_3d standalone counter.
      !! Called once in ocean_halo_face_y_3d (batched, no per-layer 2d loop).
      if (oh_suppress_depth > 0) return
      oh_cnt_face_y_3d = oh_cnt_face_y_3d + 1_int64
   end subroutine oh_count_face_y_3d

   subroutine oh_count_msgs(n)
      !! Increment the MPI isend post counter by n.
      !! NOT subject to the suppress depth: messages are physical MPI
      !! operations regardless of semantic grouping.  Called in the MPI
      !! `rdb_ocean_halo` once per exchange after each set of
      !! HALO_ISEND_N calls; n is the number of ISENDs just posted.
      integer, intent(in) :: n
      oh_cnt_msgs = oh_cnt_msgs + int(n, int64)
   end subroutine oh_count_msgs

   ! ------------------------------------------------------------------
   ! Reset
   ! ------------------------------------------------------------------

   subroutine oh_counters_reset()
      !! Reset all 9 semantic counters and the message counter to zero.
      !! Does NOT reset the suppress flag.
      oh_cnt_bt_group = 0_int64
      oh_cnt_bt_u_mid = 0_int64
      oh_cnt_ml_state = 0_int64
      oh_cnt_centre_2d = 0_int64
      oh_cnt_centre_3d = 0_int64
      oh_cnt_face_x_2d = 0_int64
      oh_cnt_face_x_3d = 0_int64
      oh_cnt_face_y_2d = 0_int64
      oh_cnt_face_y_3d = 0_int64
      oh_cnt_msgs = 0_int64
   end subroutine oh_counters_reset

   ! ------------------------------------------------------------------
   ! Getter
   ! ------------------------------------------------------------------

   subroutine oh_counters_get(bt_group, bt_u_mid, ml_state, &
                              centre_2d, centre_3d, &
                              face_x_2d, face_x_3d, &
                              face_y_2d, face_y_3d)
      !! Return all 9 semantic counters in one call.
      integer(int64), intent(out) :: bt_group
         !! bt_group_2d call count
      integer(int64), intent(out) :: bt_u_mid
         !! mid-substep bt_ubt face_x call count
      integer(int64), intent(out) :: ml_state
         !! ml_state exchange call count
      integer(int64), intent(out) :: centre_2d
         !! standalone centre_2d call count
      integer(int64), intent(out) :: centre_3d
         !! standalone centre_3d call count
      integer(int64), intent(out) :: face_x_2d
         !! standalone face_x_2d call count
      integer(int64), intent(out) :: face_x_3d
         !! standalone face_x_3d call count
      integer(int64), intent(out) :: face_y_2d
         !! standalone face_y_2d call count
      integer(int64), intent(out) :: face_y_3d
         !! standalone face_y_3d call count

      bt_group = oh_cnt_bt_group
      bt_u_mid = oh_cnt_bt_u_mid
      ml_state = oh_cnt_ml_state
      centre_2d = oh_cnt_centre_2d
      centre_3d = oh_cnt_centre_3d
      face_x_2d = oh_cnt_face_x_2d
      face_x_3d = oh_cnt_face_x_3d
      face_y_2d = oh_cnt_face_y_2d
      face_y_3d = oh_cnt_face_y_3d
   end subroutine oh_counters_get

   ! ------------------------------------------------------------------
   ! Message counter getter
   ! ------------------------------------------------------------------

   function oh_counters_msgs() result(n)
      !! Return the total MPI ISEND post count.
      !! NOT included in oh_counters_total (separate physical counter).
      integer(int64) :: n
      n = oh_cnt_msgs
   end function oh_counters_msgs

   ! ------------------------------------------------------------------
   ! Total
   ! ------------------------------------------------------------------

   function oh_counters_total() result(total)
      !! Return the sum of all 9 semantic counters.
      integer(int64) :: total
      total = oh_cnt_bt_group + oh_cnt_bt_u_mid + oh_cnt_ml_state + &
              oh_cnt_centre_2d + oh_cnt_centre_3d + &
              oh_cnt_face_x_2d + oh_cnt_face_x_3d + &
              oh_cnt_face_y_2d + oh_cnt_face_y_3d
   end function oh_counters_total

   ! ------------------------------------------------------------------
   ! Format
   ! ------------------------------------------------------------------

   function oh_counters_format() result(str)
      !! Return a one-line counter summary string.
      !!
      !! Output format (single line, <= 240 chars):
      !!   "halo exch: total=NNN bt_group=NN bt_u=NN ml=NN c2d=NN c3d=NN fx2=NN fx3=NN fy2=NN fy3=NN msgs=NN"
      !! msgs is the physical MPI ISEND count (not in total).
      character(len=256) :: str
         !! One-line summary (at most 240 chars; remainder is blank-padded).

      str = "halo exch: total="//to_string(oh_counters_total())// &
            " bt_group="//to_string(oh_cnt_bt_group)// &
            " bt_u="//to_string(oh_cnt_bt_u_mid)// &
            " ml="//to_string(oh_cnt_ml_state)// &
            " c2d="//to_string(oh_cnt_centre_2d)// &
            " c3d="//to_string(oh_cnt_centre_3d)// &
            " fx2="//to_string(oh_cnt_face_x_2d)// &
            " fx3="//to_string(oh_cnt_face_x_3d)// &
            " fy2="//to_string(oh_cnt_face_y_2d)// &
            " fy3="//to_string(oh_cnt_face_y_3d)// &
            " msgs="//to_string(oh_cnt_msgs)
   end function oh_counters_format

end module rdb_ocean_halo_counters
