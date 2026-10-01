!! Equilibrium-tide astronomy generator (mean longitudes, equilibrium
!! arguments, nodal corrections, constituent catalog).
module rdb_ocean_tide_astro
   !! Clean-room astronomy for the equilibrium (astronomical) body-force
   !! tide.  Computes the four mean longitudes (moon `s`, sun `h`, lunar
   !! perigee `p`, ascending node `N`) from Schureman's polynomials, the
   !! per-constituent equilibrium argument `V_c` (with the load-bearing
   !! ±pi/2 diurnal signs), and the slowly-varying nodal amplitude/phase
   !! corrections `f_c(N)`, `u_c(N)`.  Host-side scalar generator — called
   !! once at init (to bake the ref-date offset into `phase0`) and once per
   !! outer step (to advance the running argument); NOT a device kernel.
   !!
   !! References (the recipe is built from these, not from any model source):
   !!   * Doodson (1921), Proc. R. Soc. Lond. A 100.
   !!   * Schureman (1958), "Manual of harmonic analysis and prediction of
   !!     tides", US C&GS Spec. Pub. 98 (mean-longitude polynomials; the
   !!     36525-day Julian century).
   !!   * Cartwright & Tayler (1971) / Cartwright & Edden (1973) (amplitudes).
   !!   * Kowalik & Luick (2019), "Modern Theory and Practice of Tide
   !!     Analysis and Prediction" (Tables I.4 argument, I.6 nodal).
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   implicit none
   private

   public :: TIDE_NAME, TIDE_SPECIES, TIDE_AMP, TIDE_LOVE, TIDE_OMEGA
   public :: TIDES_CATALOG_SIZE
   public :: gregorian_day_number, days_since_1900
   public :: mean_longitudes, equilibrium_arguments, nodal_fu
   public :: tide_name_index, parse_date_string

   real(wp), parameter, public :: TIDE_PI = 4.0_wp*atan(1.0_wp)
      !! Pi to working precision.
   real(wp), parameter, public :: TIDE_DEG2RAD = TIDE_PI/180.0_wp
      !! Degrees -> radians.

   integer, parameter :: TIDES_CATALOG_SIZE = 10
      !! Full catalog: 8 default (M2 S2 N2 K2 K1 O1 P1 Q1) + MF MM branch.

   character(len=2), parameter :: TIDE_NAME(TIDES_CATALOG_SIZE) = &
                                  ["M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "MF", "MM"]
      !! Constituent names (catalog order).
   integer, parameter :: TIDE_SPECIES(TIDES_CATALOG_SIZE) = &
                         [2, 2, 2, 2, 1, 1, 1, 1, 3, 3]
      !! Species / structure-slice index: 1 diurnal, 2 semidiurnal,
      !! 3 long-period.  Directly indexes the (nx,ny,3) struct slices.
   real(wp), parameter :: TIDE_AMP(TIDES_CATALOG_SIZE) = &
                          [0.242334_wp, 0.112743_wp, 0.046397_wp, 0.030684_wp, &
                           0.141565_wp, 0.100661_wp, 0.046848_wp, 0.019273_wp, &
                           0.042041_wp, 0.022191_wp]
      !! Equilibrium amplitudes A (m).
   real(wp), parameter :: TIDE_LOVE(TIDES_CATALOG_SIZE) = &
                          [0.693_wp, 0.693_wp, 0.693_wp, 0.693_wp, &
                           0.736_wp, 0.695_wp, 0.706_wp, 0.695_wp, 0.693_wp, 0.693_wp]
      !! Effective Love-number factors (1 + k - h).
   real(wp), parameter :: TIDE_OMEGA(TIDES_CATALOG_SIZE) = &
                          [1.4051890e-4_wp, 1.4544410e-4_wp, 1.3787970e-4_wp, 1.4584234e-4_wp, &
                           0.7292117e-4_wp, 0.6759774e-4_wp, 0.7252295e-4_wp, 0.6495854e-4_wp, &
                           0.053234e-4_wp, 0.026392e-4_wp]
      !! Angular frequencies omega_c (rad/s).

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

   pure function days_since_1900(y, m, d) result(dnum)
      !! Days since the astronomical origin 1900-01-01 00:00 UT.
      integer, intent(in) :: y, m, d
      real(wp) :: dnum
      dnum = real(gregorian_day_number(y, m, d) &
                  - gregorian_day_number(1900, 1, 1), wp)
   end function days_since_1900

   pure subroutine mean_longitudes(dnum, s_deg, h_deg, p_deg, n_deg)
      !! Mean longitudes at day number `dnum` (Schureman polynomials).
      !! Returned in DEGREES, folded to [0,360).  T = dnum/36525 (Julian
      !! centuries, Schureman's 36525-day century).  Multiply by
      !! TIDE_DEG2RAD to get radians for the equilibrium arguments.
      real(wp), intent(in) :: dnum
      real(wp), intent(out) :: s_deg, h_deg, p_deg, n_deg
      real(wp) :: t
      t = dnum/36525.0_wp
      s_deg = wrap360(277.0248_wp + 481267.8906_wp*t + 0.0011_wp*t**2)
      h_deg = wrap360(280.1895_wp + 36000.7689_wp*t + 3.0310e-4_wp*t**2)
      p_deg = wrap360(334.3853_wp + 4069.0340_wp*t - 0.0103_wp*t**2)
      n_deg = wrap360(259.1568_wp - 1934.142_wp*t + 0.0021_wp*t**2)
   end subroutine mean_longitudes

   pure function wrap360(x) result(y)
      !! Fold an angle in degrees to the canonical residue [0,360).
      !! (Fortran `mod` returns the sign of the argument, so a negative
      !! polynomial value needs the +360 canonicalisation.)
      real(wp), intent(in) :: x
      real(wp) :: y
      y = mod(x, 360.0_wp)
      if (y < 0.0_wp) y = y + 360.0_wp
   end function wrap360

   pure subroutine equilibrium_arguments(dnum, v_arg)
      !! Equilibrium argument `V_c` (radians) for every catalog
      !! constituent at day number `dnum`.  s,h,p,N are taken in radians
      !! (deg-mod-360 -> rad); the result is left un-modded (cos/sin are
      !! periodic).  The ±pi/2 diurnal signs are load-bearing.
      real(wp), intent(in) :: dnum
      real(wp), intent(out) :: v_arg(TIDES_CATALOG_SIZE)
      real(wp) :: s_deg, h_deg, p_deg, n_deg
      real(wp) :: s, h, p
      real(wp) :: half_pi
      call mean_longitudes(dnum, s_deg, h_deg, p_deg, n_deg)
      s = s_deg*TIDE_DEG2RAD
      h = h_deg*TIDE_DEG2RAD
      p = p_deg*TIDE_DEG2RAD
      half_pi = 0.5_wp*TIDE_PI
      v_arg(1) = 2.0_wp*(h - s)                    ! M2
      v_arg(2) = 0.0_wp                            ! S2
      v_arg(3) = -3.0_wp*s + 2.0_wp*h + p          ! N2
      v_arg(4) = 2.0_wp*h                          ! K2
      v_arg(5) = h + half_pi                       ! K1  (+pi/2)
      v_arg(6) = -2.0_wp*s + h - half_pi           ! O1  (-pi/2)
      v_arg(7) = -h - half_pi                      ! P1  (-pi/2)
      v_arg(8) = -3.0_wp*s + h + p - half_pi       ! Q1  (-pi/2)
      v_arg(9) = 2.0_wp*s                          ! MF
      v_arg(10) = s - p                            ! MM
   end subroutine equilibrium_arguments

   pure subroutine nodal_fu(dnum, add_nodal, f_nodal, u_nodal)
      !! Nodal amplitude factor `f_c(N)` (nondim) and phase `u_c(N)`
      !! (radians), fixed at the nodal reference date's `N`.
      !! `add_nodal = .false.` returns f=1, u=0 for every constituent.
      real(wp), intent(in) :: dnum
      logical, intent(in) :: add_nodal
      real(wp), intent(out) :: f_nodal(TIDES_CATALOG_SIZE)
      real(wp), intent(out) :: u_nodal(TIDES_CATALOG_SIZE)
      real(wp) :: s_deg, h_deg, p_deg, n_deg
      real(wp) :: cos_n, sin_n
      integer :: c
      if (.not. add_nodal) then
         do c = 1, TIDES_CATALOG_SIZE
            f_nodal(c) = 1.0_wp
            u_nodal(c) = 0.0_wp
         end do
         return
      end if
      call mean_longitudes(dnum, s_deg, h_deg, p_deg, n_deg)
      cos_n = cos(n_deg*TIDE_DEG2RAD)
      sin_n = sin(n_deg*TIDE_DEG2RAD)
      ! M2, N2
      f_nodal(1) = 1.0_wp - 0.037_wp*cos_n
      u_nodal(1) = (-2.1_wp)*TIDE_DEG2RAD*sin_n
      f_nodal(3) = 1.0_wp - 0.037_wp*cos_n
      u_nodal(3) = (-2.1_wp)*TIDE_DEG2RAD*sin_n
      ! K2
      f_nodal(4) = 1.024_wp + 0.286_wp*cos_n
      u_nodal(4) = (-17.7_wp)*TIDE_DEG2RAD*sin_n
      ! K1
      f_nodal(5) = 1.006_wp + 0.115_wp*cos_n
      u_nodal(5) = (-8.9_wp)*TIDE_DEG2RAD*sin_n
      ! O1, Q1  (NOTE +sign on u)
      f_nodal(6) = 1.009_wp + 0.187_wp*cos_n
      u_nodal(6) = (10.8_wp)*TIDE_DEG2RAD*sin_n
      f_nodal(8) = 1.009_wp + 0.187_wp*cos_n
      u_nodal(8) = (10.8_wp)*TIDE_DEG2RAD*sin_n
      ! S2, P1
      f_nodal(2) = 1.0_wp
      u_nodal(2) = 0.0_wp
      f_nodal(7) = 1.0_wp
      u_nodal(7) = 0.0_wp
      ! MF
      f_nodal(9) = 1.043_wp + 0.414_wp*cos_n
      u_nodal(9) = (-23.7_wp)*TIDE_DEG2RAD*sin_n
      ! MM
      f_nodal(10) = 1.0_wp - 0.130_wp*cos_n
      u_nodal(10) = 0.0_wp
   end subroutine nodal_fu

   pure function tide_name_index(name) result(idx)
      !! Catalog index (1..TIDES_CATALOG_SIZE) of a constituent name,
      !! matched case-insensitively; -1 if unknown.
      character(len=*), intent(in) :: name
      integer :: idx
      integer :: c
      character(len=2) :: up
      up = upcase2(adjustl(name))
      idx = -1
      do c = 1, TIDES_CATALOG_SIZE
         if (up == TIDE_NAME(c)) then
            idx = c
            return
         end if
      end do
   end function tide_name_index

   pure function upcase2(s) result(u)
      !! Uppercase the first two characters of a token.
      character(len=*), intent(in) :: s
      character(len=2) :: u
      integer :: k, ic
      u = "  "
      do k = 1, min(2, len_trim(s))
         ic = iachar(s(k:k))
         if (ic >= iachar("a") .and. ic <= iachar("z")) ic = ic - 32
         u(k:k) = achar(ic)
      end do
   end function upcase2

   pure subroutine parse_date_string(str, y, m, d, ok)
      !! Parse a "YYYY-MM-DD" calendar date.  `ok = .false.` on any
      !! malformed field (caller decides fail-loud policy).
      character(len=*), intent(in) :: str
      integer, intent(out) :: y, m, d
      logical, intent(out) :: ok
      character(len=:), allocatable :: t
      character(len=128) :: emsg
      integer :: p1, p2, ios
      y = 0
      m = 0
      d = 0
      ok = .false.
      t = trim(adjustl(str))
      p1 = index(t, "-")
      if (p1 <= 1) return
      p2 = index(t(p1 + 1:), "-")
      if (p2 <= 1) return
      p2 = p1 + p2
      read (t(1:p1 - 1), *, iostat=ios, iomsg=emsg) y
      if (ios /= 0) return
      read (t(p1 + 1:p2 - 1), *, iostat=ios, iomsg=emsg) m
      if (ios /= 0) return
      read (t(p2 + 1:), *, iostat=ios, iomsg=emsg) d
      if (ios /= 0) return
      if (m < 1 .or. m > 12 .or. d < 1 .or. d > 31) return
      ok = .true.
   end subroutine parse_date_string

end module rdb_ocean_tide_astro
