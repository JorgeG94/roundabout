!! Proleptic-Gregorian calendar arithmetic (date<->seconds-since-epoch,
!! CF time-units parsing, repeat-year date-arithmetic wrap).
module rdb_calendar
   !! Pure calendar primitives shared by the astronomical-tide generator
   !! (`rdb_ocean_tide_astro`, which re-exports `gregorian_day_number`/
   !! `days_since_1900` from here rather than carrying its own copy) and
   !! the absolute-date / multi-file time axis of the time-varying
   !! data-input reader (`rdb_ocean_data_input`, OM3 wave-1 PR-C2a).
   !!
   !! Proleptic Gregorian only: every CF `calendar` attribute this tree
   !! reads in practice ("gregorian"/"standard" for any real date, which
   !! this century always is; "proleptic_gregorian" for the same rule
   !! extended before 1582) maps onto ONE day-count rule --
   !! `gregorian_day_number`'s. `noleap`/`360_day` are a genuinely
   !! different day-count and must fail loud
   !! (`calendar_name_is_implemented`), never silently alias here.
   !!
   !! References: Fliegel & Van Flandern (1968), Comm. ACM 11(10) (the
   !! Julian-day-number <-> (y,m,d) algorithm, both directions); CF
   !! Conventions §4.4 (the "<unit> since <reference-date>" time-units
   !! syntax parsed by `parse_time_units`).
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use pic_ascii, only: to_lower
   implicit none
   private

   public :: date_t
   public :: gregorian_day_number, day_number_to_date, days_since_1900
   public :: date_to_seconds_since, seconds_to_date
   public :: parse_date, parse_time_units
   public :: calendar_name_is_implemented
   public :: calendar_ryf_wrap_time

   type :: date_t
      !! Proleptic-Gregorian calendar timestamp. `s` carries fractional
      !! seconds so a CF "seconds since ..." axis round-trips exactly.
      integer :: y = 1970, m = 1, d = 1, h = 0, mi = 0
      real(wp) :: s = 0.0_wp
   end type date_t

contains

   pure function gregorian_day_number(y, m, d) result(jdn)
      !! Proleptic-Gregorian Julian Day Number (integer, at 00:00 UT).
      !! Fliegel & Van Flandern algorithm.
      integer, intent(in) :: y, m, d
      integer(int64) :: jdn
      integer(int64) :: a, yy, mm
      a = int((14 - m)/12, int64)
      yy = int(y, int64) + 4800_int64 - a
      mm = int(m, int64) + 12_int64*a - 3_int64
      jdn = int(d, int64) + (153_int64*mm + 2_int64)/5_int64 + 365_int64*yy &
            + yy/4_int64 - yy/100_int64 + yy/400_int64 - 32045_int64
   end function gregorian_day_number

   pure subroutine day_number_to_date(jdn, y, m, d)
      !! Inverse of `gregorian_day_number` (Fliegel & Van Flandern's
      !! companion integer-arithmetic inverse) -- the piece
      !! `gregorian_day_number` alone never needed, now load-bearing for
      !! `seconds_to_date`.
      integer(int64), intent(in) :: jdn
      integer, intent(out) :: y, m, d
      integer(int64) :: l, n, i, j, k
      l = jdn + 68569_int64
      n = (4_int64*l)/146097_int64
      l = l - (146097_int64*n + 3_int64)/4_int64
      i = (4000_int64*(l + 1_int64))/1461001_int64
      l = l - (1461_int64*i)/4_int64 + 31_int64
      j = (80_int64*l)/2447_int64
      d = int(l - (2447_int64*j)/80_int64)
      k = j/11_int64
      m = int(j + 2_int64 - 12_int64*k)
      y = int(100_int64*(n - 49_int64) + i + k)
   end subroutine day_number_to_date

   pure function days_since_1900(y, m, d) result(dnum)
      !! Days since the astronomical origin 1900-01-01 00:00 UT.
      integer, intent(in) :: y, m, d
      real(wp) :: dnum
      dnum = real(gregorian_day_number(y, m, d) &
                  - gregorian_day_number(1900, 1, 1), wp)
   end function days_since_1900

   pure function date_to_seconds_since(date, epoch) result(t)
      !! Seconds from `epoch` to `date` (negative if `date` precedes
      !! `epoch`). Pure calendar-day arithmetic -- no leap seconds,
      !! matching every CF/UDUNITS clock this tree reads.
      type(date_t), intent(in) :: date, epoch
      real(wp) :: t
      integer(int64) :: ddiff
      ddiff = gregorian_day_number(date%y, date%m, date%d) &
              - gregorian_day_number(epoch%y, epoch%m, epoch%d)
      t = real(ddiff, wp)*86400.0_wp &
          + real(date%h - epoch%h, wp)*3600.0_wp &
          + real(date%mi - epoch%mi, wp)*60.0_wp &
          + (date%s - epoch%s)
   end function date_to_seconds_since

   pure subroutine seconds_to_date(t, epoch, date)
      !! Inverse of `date_to_seconds_since`: `epoch + t` seconds -> date.
      real(wp), intent(in) :: t
      type(date_t), intent(in) :: epoch
      type(date_t), intent(out) :: date
      integer(int64) :: epoch_jdn, day_offset, jdn
      real(wp) :: epoch_sec_of_day, total_sec, sec_of_day

      epoch_jdn = gregorian_day_number(epoch%y, epoch%m, epoch%d)
      epoch_sec_of_day = real(epoch%h, wp)*3600.0_wp + real(epoch%mi, wp)*60.0_wp + epoch%s
      total_sec = epoch_sec_of_day + t
      day_offset = int(floor(total_sec/86400.0_wp), int64)
      sec_of_day = total_sec - real(day_offset, wp)*86400.0_wp
      jdn = epoch_jdn + day_offset

      call day_number_to_date(jdn, date%y, date%m, date%d)
      date%h = int(sec_of_day/3600.0_wp)
      date%mi = int((sec_of_day - real(date%h, wp)*3600.0_wp)/60.0_wp)
      date%s = sec_of_day - real(date%h, wp)*3600.0_wp - real(date%mi, wp)*60.0_wp
   end subroutine seconds_to_date

   pure subroutine parse_date(str, date, ok)
      !! Parse "YYYY-MM-DD" or "YYYY-MM-DD hh:mm:ss" (CF reference-date
      !! syntax). `ok = .false.` on any malformed field; `date` is left
      !! at the `date_t` default (1970-01-01 00:00:00) in that case.
      character(len=*), intent(in) :: str
      type(date_t), intent(out) :: date
      logical, intent(out) :: ok
      character(len=:), allocatable :: t, time_part
      integer :: p1, p2, psp, ios
      character(len=256) :: emsg

      date%y = 1970
      date%m = 1
      date%d = 1
      date%h = 0
      date%mi = 0
      date%s = 0.0_wp
      ok = .false.
      t = trim(adjustl(str))
      if (len(t) == 0) return

      psp = index(t, " ")
      if (psp > 0) then
         time_part = trim(adjustl(t(psp + 1:)))
         t = t(1:psp - 1)
      else
         time_part = ""
      end if

      p1 = index(t, "-")
      if (p1 <= 1) return
      p2 = index(t(p1 + 1:), "-")
      if (p2 <= 1) return
      p2 = p1 + p2
      read (t(1:p1 - 1), *, iostat=ios, iomsg=emsg) date%y
      if (ios /= 0) return
      read (t(p1 + 1:p2 - 1), *, iostat=ios, iomsg=emsg) date%m
      if (ios /= 0) return
      read (t(p2 + 1:), *, iostat=ios, iomsg=emsg) date%d
      if (ios /= 0) return
      if (date%m < 1 .or. date%m > 12 .or. date%d < 1 .or. date%d > 31) return

      if (len(time_part) > 0) then
         call parse_time_of_day(time_part, date%h, date%mi, date%s, ok)
         if (.not. ok) return
      end if
      ok = .true.
   end subroutine parse_date

   pure subroutine parse_time_of_day(str, h, mi, s, ok)
      !! "hh:mm:ss[.sss]" -> (h, mi, s). Helper for `parse_date`.
      character(len=*), intent(in) :: str
      integer, intent(out) :: h, mi
      real(wp), intent(out) :: s
      logical, intent(out) :: ok
      integer :: p1, p2, ios
      character(len=256) :: emsg
      h = 0
      mi = 0
      s = 0.0_wp
      ok = .false.
      p1 = index(str, ":")
      if (p1 <= 1) return
      p2 = index(str(p1 + 1:), ":")
      if (p2 <= 1) return
      p2 = p1 + p2
      read (str(1:p1 - 1), *, iostat=ios, iomsg=emsg) h
      if (ios /= 0) return
      read (str(p1 + 1:p2 - 1), *, iostat=ios, iomsg=emsg) mi
      if (ios /= 0) return
      read (str(p2 + 1:), *, iostat=ios, iomsg=emsg) s
      if (ios /= 0) return
      if (h < 0 .or. h > 23 .or. mi < 0 .or. mi > 59 .or. s < 0.0_wp .or. s >= 60.0_wp) return
      ok = .true.
   end subroutine parse_time_of_day

   pure subroutine parse_time_units(units, scale, epoch, ok)
      !! CF `units` attribute ("<second|minute|hour|day>[s] since
      !! <reference-date>") -> the SI-second multiplier for the raw file
      !! axis + the reference date as a `date_t`. `ok = .false.` on an
      !! unrecognised leading unit word, a missing " since ", or an
      !! unparsable reference date; `scale = 1`, `epoch` left at its
      !! 1970-01-01 default in that case.
      character(len=*), intent(in) :: units
      real(wp), intent(out) :: scale
      type(date_t), intent(out) :: epoch
      logical, intent(out) :: ok
      character(len=:), allocatable :: s, u
      integer :: p

      scale = 1.0_wp
      epoch%y = 1970
      epoch%m = 1
      epoch%d = 1
      epoch%h = 0
      epoch%mi = 0
      epoch%s = 0.0_wp
      ok = .true.

      s = trim(adjustl(units))
      u = to_lower(s)
      if (index(u, "second") == 1 .or. index(u, "sec") == 1) then
         scale = 1.0_wp
      else if (index(u, "minute") == 1 .or. index(u, "min") == 1) then
         scale = 60.0_wp
      else if (index(u, "hour") == 1) then
         scale = 3600.0_wp
      else if (index(u, "day") == 1) then
         scale = 86400.0_wp
      else
         ok = .false.
         return
      end if

      p = index(u, " since ")
      if (p <= 0) then
         ok = .false.
         return
      end if
      call parse_date(trim(adjustl(s(p + 7:))), epoch, ok)
   end subroutine parse_time_units

   pure logical function calendar_name_is_implemented(tag) result(ok)
      !! `.true.` iff `tag` (case-insensitive) names the one calendar
      !! this tree implements -- proleptic Gregorian, under any of its
      !! three CF spellings. `noleap`/`360_day` are a different
      !! day-count rule and must fail loud, never silently alias here.
      character(len=*), intent(in) :: tag
      select case (to_lower(trim(tag)))
      case ("gregorian", "standard", "proleptic_gregorian")
         ok = .true.
      case default
         ok = .false.
      end select
   end function calendar_name_is_implemented

   pure subroutine calendar_ryf_wrap_time(t_query, ryf_start_t, period_seconds, t_wrapped)
      !! Repeat-year-forcing date arithmetic: fold an absolute query
      !! time `t_query` (seconds since the model epoch) into the window
      !! `[ryf_start_t, ryf_start_t + period_seconds)`. `ryf_start_t`
      !! itself comes from calendar date arithmetic
      !! (`date_to_seconds_since` on the RYF start date), which is the
      !! part a plain seconds-only period cannot get right on its own --
      !! once anchored, the fold itself is an ordinary modulo. The
      !! wrapped result still needs a CYCLIC bracket search against the
      !! RYF window's own time axis to interpolate across the seam
      !! between its last and first record -- that composition is the
      !! reader's job (`rdb_ocean_data_input`), not this pure helper's.
      real(wp), intent(in) :: t_query, ryf_start_t, period_seconds
      real(wp), intent(out) :: t_wrapped
      t_wrapped = ryf_start_t + modulo(t_query - ryf_start_t, period_seconds)
   end subroutine calendar_ryf_wrap_time

end module rdb_calendar
