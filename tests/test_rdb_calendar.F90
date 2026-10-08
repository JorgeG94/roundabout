!! Unit tests for `rdb_calendar` (OM3 wave-1 PR-C2a: proleptic-Gregorian
!! calendar arithmetic for the absolute-date / multi-file time axis).
!!
!! Coverage:
!!   * `calendar_day_roundtrip`    — `gregorian_day_number` <->
!!                                   `day_number_to_date`, incl. a leap
!!                                   day (2000-02-29) and the Gregorian
!!                                   century rule (1900-03-01 is day 60
!!                                   of 1900, which is NOT a leap year).
!!   * `calendar_parse_date`       — "YYYY-MM-DD" and "YYYY-MM-DD
!!                                   hh:mm:ss"; malformed inputs report
!!                                   `ok = .false.`.
!!   * `calendar_parse_time_units` — CF `units` attribute parsing
!!                                   ("days since ...", "hours since
!!                                   ..."); an unrecognised leading word
!!                                   or missing "since" is `ok = .false.`.
!!   * `calendar_seconds_date_roundtrip` — `date_to_seconds_since` <->
!!                                   `seconds_to_date`, exact to
!!                                   round-off, across a multi-decade
!!                                   span (the JRA55-do IAF record is
!!                                   67 years, ~2.1e9 s).
!!   * `calendar_name_implemented` — gregorian/standard/
!!                                   proleptic_gregorian accepted;
!!                                   noleap/360_day refused.
!!   * `calendar_ryf_wrap`         — `calendar_ryf_wrap_time`: a query
!!                                   inside the window is unchanged; one
!!                                   full period later wraps back to the
!!                                   same in-window value; several
!!                                   periods and a negative offset both
!!                                   land in `[ryf_start_t, ryf_start_t +
!!                                   period)`.
module test_rdb_calendar
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_calendar, only: date_t, gregorian_day_number, day_number_to_date, &
                           date_to_seconds_since, seconds_to_date, &
                           parse_date, parse_time_units, &
                           calendar_name_is_implemented, calendar_ryf_wrap_time
   implicit none
   private

   public :: collect_rdb_calendar_tests

contains

   subroutine collect_rdb_calendar_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("calendar_day_roundtrip", test_day_roundtrip), &
                  new_unittest("calendar_parse_date", test_parse_date), &
                  new_unittest("calendar_parse_time_units", test_parse_time_units), &
                  new_unittest("calendar_seconds_date_roundtrip", test_seconds_date_roundtrip), &
                  new_unittest("calendar_name_implemented", test_name_implemented), &
                  new_unittest("calendar_ryf_wrap", test_ryf_wrap) &
                  ]
   end subroutine collect_rdb_calendar_tests

   subroutine test_day_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCASE = 4
      integer :: y0(NCASE), m0(NCASE), d0(NCASE)
      integer :: y, m, d, k

      ! 1958-01-01 (JRA55-do IAF epoch), 2000-02-29 (leap day), 1900-03-01
      ! (the Gregorian century rule: 1900 is NOT a leap year, so this is
      ! day 60 of 1900, not day 61), 2024-01-01 (a recent leap year's
      ! start).
      y0 = [1958, 2000, 1900, 2024]
      m0 = [1, 2, 3, 1]
      d0 = [1, 29, 1, 1]

      do k = 1, NCASE
         call day_number_to_date(gregorian_day_number(y0(k), m0(k), d0(k)), y, m, d)
         call check(error, y == y0(k), "day roundtrip: year mismatch")
         if (allocated(error)) return
         call check(error, m == m0(k), "day roundtrip: month mismatch")
         if (allocated(error)) return
         call check(error, d == d0(k), "day roundtrip: day mismatch")
         if (allocated(error)) return
      end do

      ! Century rule, checked directly: 1900-03-01 is day number
      ! 1900-01-01 + 59 (31 Jan + 28 Feb, NOT 29 — 1900 is not a leap
      ! year under the proleptic-Gregorian rule this tree implements).
      call check(error, gregorian_day_number(1900, 3, 1) - gregorian_day_number(1900, 1, 1) == 59, &
                 "century rule: 1900 must not be treated as a leap year")
      if (allocated(error)) return

      ! 2000 IS a leap year (divisible by 400): Feb has 29 days.
      call check(error, gregorian_day_number(2000, 3, 1) - gregorian_day_number(2000, 1, 1) == 60, &
                 "2000 must be treated as a leap year (divisible by 400)")
      if (allocated(error)) return
   end subroutine test_day_roundtrip

   subroutine test_parse_date(error)
      type(error_type), allocatable, intent(out) :: error
      type(date_t) :: date
      logical :: ok

      call parse_date("1958-01-01", date, ok)
      call check(error, ok, "parse_date: plain date should parse")
      if (allocated(error)) return
      call check(error, date%y == 1958 .and. date%m == 1 .and. date%d == 1, &
                 "parse_date: plain date fields wrong")
      if (allocated(error)) return

      call parse_date("1990-05-01 12:30:45", date, ok)
      call check(error, ok, "parse_date: date+time should parse")
      if (allocated(error)) return
      call check(error, date%y == 1990 .and. date%m == 5 .and. date%d == 1 .and. &
                 date%h == 12 .and. date%mi == 30 .and. abs(date%s - 45.0_wp) < 1.0e-10_wp, &
                 "parse_date: date+time fields wrong")
      if (allocated(error)) return

      call parse_date("not-a-date", date, ok)
      call check(error,.not. ok, "parse_date: malformed date must report ok=.false.")
      if (allocated(error)) return

      call parse_date("", date, ok)
      call check(error,.not. ok, "parse_date: empty string must report ok=.false.")
      if (allocated(error)) return
   end subroutine test_parse_date

   subroutine test_parse_time_units(error)
      type(error_type), allocatable, intent(out) :: error
      type(date_t) :: epoch
      real(wp) :: scale
      logical :: ok

      call parse_time_units("days since 1958-01-01 00:00:00", scale, epoch, ok)
      call check(error, ok, "parse_time_units: 'days since ...' should parse")
      if (allocated(error)) return
      call check(error, abs(scale - 86400.0_wp) < 1.0e-10_wp, "parse_time_units: day scale wrong")
      if (allocated(error)) return
      call check(error, epoch%y == 1958 .and. epoch%m == 1 .and. epoch%d == 1, &
                 "parse_time_units: epoch date wrong")
      if (allocated(error)) return

      call parse_time_units("hours since 1900-01-01", scale, epoch, ok)
      call check(error, ok, "parse_time_units: 'hours since ...' (no time-of-day) should parse")
      if (allocated(error)) return
      call check(error, abs(scale - 3600.0_wp) < 1.0e-10_wp, "parse_time_units: hour scale wrong")
      if (allocated(error)) return

      call parse_time_units("furlongs since 1900-01-01", scale, epoch, ok)
      call check(error,.not. ok, "parse_time_units: unrecognised unit must report ok=.false.")
      if (allocated(error)) return

      call parse_time_units("days 1958-01-01", scale, epoch, ok)
      call check(error,.not. ok, "parse_time_units: missing 'since' must report ok=.false.")
      if (allocated(error)) return
   end subroutine test_parse_time_units

   subroutine test_seconds_date_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      type(date_t) :: epoch, mid_date, date_out
      real(wp) :: t, t_out
      integer :: k
      real(wp) :: offsets(5)
      logical :: ok_discard

      epoch%y = 1958; epoch%m = 1; epoch%d = 1; epoch%h = 0; epoch%mi = 0; epoch%s = 0.0_wp

      ! Spans the ~67-year JRA55-do IAF record (1958..2024), incl. a
      ! negative offset (before the epoch).
      offsets = [-86400.0_wp, 0.0_wp, 12345.5_wp, 365.0_wp*86400.0_wp, &
                 67.0_wp*365.25_wp*86400.0_wp]

      do k = 1, size(offsets)
         t = offsets(k)
         call seconds_to_date(t, epoch, mid_date)
         t_out = date_to_seconds_since(mid_date, epoch)
         call check(error, abs(t_out - t) < 1.0e-6_wp, "seconds<->date roundtrip not exact")
         if (allocated(error)) return
      end do

      ! Date -> seconds -> date roundtrip, anchored on a date far from
      ! the epoch (exercises the day-number inverse on a large JDN).
      call parse_date("2024-01-01 06:00:00", mid_date, ok_discard)
      t = date_to_seconds_since(mid_date, epoch)
      call seconds_to_date(t, epoch, date_out)
      call check(error, date_out%y == 2024 .and. date_out%m == 1 .and. date_out%d == 1 .and. &
                 date_out%h == 6, "date->seconds->date roundtrip wrong")
      if (allocated(error)) return
   end subroutine test_seconds_date_roundtrip

   subroutine test_name_implemented(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, calendar_name_is_implemented("gregorian"), "gregorian must be implemented")
      if (allocated(error)) return
      call check(error, calendar_name_is_implemented("STANDARD"), "standard (any case) must be implemented")
      if (allocated(error)) return
      call check(error, calendar_name_is_implemented("proleptic_gregorian"), &
                 "proleptic_gregorian must be implemented")
      if (allocated(error)) return
      call check(error,.not. calendar_name_is_implemented("noleap"), &
                 "noleap must NOT be implemented (different day-count rule)")
      if (allocated(error)) return
      call check(error,.not. calendar_name_is_implemented("360_day"), &
                 "360_day must NOT be implemented (different day-count rule)")
      if (allocated(error)) return
   end subroutine test_name_implemented

   subroutine test_ryf_wrap(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: RYF_START = 1000.0_wp, PERIOD = 400.0_wp*86400.0_wp
      real(wp) :: t_in, t_wrapped

      ! Inside the window: unchanged.
      t_in = RYF_START + 50.0_wp*86400.0_wp
      call calendar_ryf_wrap_time(t_in, RYF_START, PERIOD, t_wrapped)
      call check(error, abs(t_wrapped - t_in) < 1.0e-9_wp, "ryf_wrap: in-window query must be unchanged")
      if (allocated(error)) return

      ! Exactly one period later -> wraps back to the same in-window value.
      call calendar_ryf_wrap_time(t_in + PERIOD, RYF_START, PERIOD, t_wrapped)
      call check(error, abs(t_wrapped - t_in) < 1.0e-6_wp, &
                 "ryf_wrap: one period later must wrap to the same value")
      if (allocated(error)) return

      ! Several periods later, and a negative multiple-of-period offset,
      ! both land in [RYF_START, RYF_START + PERIOD).
      call calendar_ryf_wrap_time(t_in + 7.0_wp*PERIOD, RYF_START, PERIOD, t_wrapped)
      call check(error, abs(t_wrapped - t_in) < 1.0e-6_wp, "ryf_wrap: 7 periods later must wrap identically")
      if (allocated(error)) return
      call check(error, t_wrapped >= RYF_START .and. t_wrapped < RYF_START + PERIOD, &
                 "ryf_wrap: result must land inside the window")
      if (allocated(error)) return

      call calendar_ryf_wrap_time(t_in - 3.0_wp*PERIOD, RYF_START, PERIOD, t_wrapped)
      call check(error, abs(t_wrapped - t_in) < 1.0e-6_wp, &
                 "ryf_wrap: negative multiple-of-period offset must wrap identically")
      if (allocated(error)) return
   end subroutine test_ryf_wrap

end module test_rdb_calendar
