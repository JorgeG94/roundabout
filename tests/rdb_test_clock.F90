!! Wall-clock seconds for test timing, through pic_timer (no system_clock).
module rdb_test_clock
   use, intrinsic :: iso_fortran_env, only: real64
   use pic_timer, only: timer_type
   implicit none
   private
   public :: wall_seconds

   type(timer_type), save :: epoch
   logical, save :: epoch_started = .false.

contains

   function wall_seconds() result(t)
      !! Wall seconds since the first call; subtract two readings for an interval.
      real(real64) :: t
      if (.not. epoch_started) then
         call epoch%start()
         epoch_started = .true.
      end if
      t = epoch%get_elapsed_time()
   end function wall_seconds

end module rdb_test_clock
