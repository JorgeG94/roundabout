!! Analytic tests for the sea-ice enthalpy library (`rdb_ice_enthalpy`,
!! sea-ice PR 2) — the closed-form T<->E inversions with
!! Cp_brine == Cp_ice.
!!
!! Cases:
!!   * `enthalpy_round_trip_identity` — the headline case: for a grid
!!     of salinities and ice temperatures, ice_temp_from_en_s(
!!     ice_enth_from_ts(T, S), S) == T to 1e-14 relative — a real
!!     algebraic identity, not a loosened solver tolerance.
!!   * `freezing_point_linear` — ice_t_freeze is the exact linear
!!     liquidus.
!!   * `freeze_point_continuity` — the ice branch of ice_enth_from_ts
!!     matches ice_enthalpy_liquid_freeze exactly at T = t_fr.
!!   * `melted_branch_matches_liquid` — above the freezing point,
!!     ice_enth_from_ts collapses to ice_enthalpy_liquid, and the
!!     inverse recovers T.
!!   * `fresh_water_sub_branches` — S = 0 exercises all three
!!     ice_temp_from_en_s fresh-water sub-branches (liquid / mushy
!!     plateau / solid).
!!   * `monotonic_and_latent_jump` — ice_enth_from_ts is monotonically
!!     increasing in T, and the enthalpy deficit at the ice/liquid
!!     midpoint carries the LAT_FUS mushy signature (loose physics
!!     sanity check, not an identity).
!!   * `elemental_dispatch` — whole-array calls match per-element
!!     scalar calls bitwise, guarding the `elemental` attribute.
module test_ocean_ice_enthalpy
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_enthalpy, only: ice_t_freeze, ice_enth_from_ts, &
                               ice_enthalpy_liquid_freeze, ice_enthalpy_liquid, &
                               ice_temp_from_en_s, ICE_LAT_FUS, ICE_CP_ICE, &
                               ICE_CP_WATER, ICE_DTF_DS, ICE_ENTH_LIQ_0
   implicit none
   private

   public :: collect_ocean_ice_enthalpy_tests

contains

   subroutine collect_ocean_ice_enthalpy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("enthalpy_round_trip_identity", test_round_trip), &
                  new_unittest("freezing_point_linear", test_freezing_point_linear), &
                  new_unittest("freeze_point_continuity", test_freeze_point_continuity), &
                  new_unittest("melted_branch_matches_liquid", test_melted_branch), &
                  new_unittest("fresh_water_sub_branches", test_fresh_water_sub_branches), &
                  new_unittest("monotonic_and_latent_jump", test_monotonic_and_latent_jump), &
                  new_unittest("elemental_dispatch", test_elemental_dispatch) &
                  ]
   end subroutine collect_ocean_ice_enthalpy_tests

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_round_trip(error)
      !! Headline case: for S in {0, 5, 10, 20, 30, 35} and, per S, 40 T
      !! values spanning [-40, t_fr) (endpoint strictly below t_fr so
      !! every sample lands in the ice branch), the round trip
      !! ice_temp_from_en_s(ice_enth_from_ts(T, S), S) == T holds to
      !! 1e-14 relative — an algebraic identity (see module docstring
      !! in rdb_ice_enthalpy for why it is exact, not merely tight).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS = 6
      real(wp), parameter :: s_grid(NS) = [0.0_wp, 5.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 35.0_wp]
      integer, parameter :: NT = 40
      real(wp) :: s, t, t_fr, t_hi, en, t_back, tol
      integer :: is, it

      do is = 1, NS
         s = s_grid(is)
         t_fr = ice_t_freeze(s)
         t_hi = t_fr - 1.0e-6_wp
         do it = 1, NT
            t = -40.0_wp + (t_hi - (-40.0_wp))*real(it - 1, wp)/real(NT - 1, wp)
            en = ice_enth_from_ts(t, s)
            t_back = ice_temp_from_en_s(en, s)
            tol = 1.0e-14_wp*(abs(t) + 1.0_wp)
            call check(error, abs(t_back - t) < tol, &
                       "round trip must recover T to 1e-14 relative")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_round_trip

   subroutine test_freezing_point_linear(error)
      !! ice_t_freeze(s) == ICE_DTF_DS*s exactly. Values are routed
      !! through separate intermediate variables (not compared as one
      !! fused inline expression) so `-O3` FMA contraction cannot fuse
      !! the multiply-subtract into a differently-rounded single
      !! instruction across the two independently-evaluated sides.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: s_grid(3) = [0.0_wp, 17.5_wp, 35.0_wp]
      integer :: is
      real(wp) :: lhs, rhs
      real(wp) :: t_zero

      do is = 1, size(s_grid)
         lhs = ice_t_freeze(s_grid(is))
         rhs = ICE_DTF_DS*s_grid(is)
         call check(error, lhs == rhs, &
                    "ice_t_freeze must equal ICE_DTF_DS*s exactly")
         if (allocated(error)) return
      end do

      t_zero = ice_t_freeze(0.0_wp)
      call check(error, t_zero == 0.0_wp, &
                 "ice_t_freeze(0) must be exactly 0")
   end subroutine test_freezing_point_linear

   subroutine test_freeze_point_continuity(error)
      !! At T = t_fr, the ice branch of ice_enth_from_ts matches
      !! ice_enthalpy_liquid_freeze(s).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: s_grid(2) = [5.0_wp, 35.0_wp]
      real(wp) :: s, t_fr, lhs, rhs, tol
      integer :: is

      do is = 1, size(s_grid)
         s = s_grid(is)
         t_fr = ice_t_freeze(s)
         lhs = ice_enth_from_ts(t_fr, s)
         rhs = ice_enthalpy_liquid_freeze(s)
         tol = 1.0e-12_wp*(abs(rhs) + 1.0_wp)
         call check(error, abs(lhs - rhs) <= tol, &
                    "enth_from_ts at t_fr must match enthalpy_liquid_freeze")
         if (allocated(error)) return
      end do
   end subroutine test_freeze_point_continuity

   subroutine test_melted_branch(error)
      !! Above the freezing point, ice_enth_from_ts collapses exactly
      !! to ice_enthalpy_liquid, and the enthalpy inverse recovers T.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: t_grid(3) = [2.0_wp, 0.5_wp, 10.0_wp]
      real(wp), parameter :: s_grid(3) = [35.0_wp, 5.0_wp, 0.0_wp]
      real(wp) :: t, s, en, en_liq, t_back, tol
      integer :: i

      do i = 1, size(t_grid)
         t = t_grid(i)
         s = s_grid(i)
         en = ice_enth_from_ts(t, s)
         en_liq = ice_enthalpy_liquid(t, s)
         tol = 1.0e-14_wp*(abs(en_liq) + 1.0_wp)
         call check(error, abs(en - en_liq) <= tol, &
                    "melted branch must match ice_enthalpy_liquid")
         if (allocated(error)) return

         t_back = ice_temp_from_en_s(en, s)
         tol = 1.0e-14_wp*(abs(t) + 1.0_wp)
         call check(error, abs(t_back - t) < tol, &
                    "melted branch enthalpy must invert back to T")
         if (allocated(error)) return
      end do
   end subroutine test_melted_branch

   subroutine test_fresh_water_sub_branches(error)
      !! S = 0 exercises all three ice_temp_from_en_s fresh-water
      !! sub-branches: liquid (en_j >= 0), mushy plateau
      !! (-LAT_FUS <= en_j < 0), and solid (en_j < -LAT_FUS).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: en, t_back, tol, t_mushy

      ! (a) Liquid leg: T = 5 round-trips.
      en = ice_enth_from_ts(5.0_wp, 0.0_wp)
      t_back = ice_temp_from_en_s(en, 0.0_wp)
      tol = 1.0e-14_wp*(abs(5.0_wp) + 1.0_wp)
      call check(error, abs(t_back - 5.0_wp) < tol, &
                 "fresh-water liquid leg must round-trip")
      if (allocated(error)) return

      ! (b) Mushy plateau leg: en_j in (-LAT_FUS, 0) maps to T = 0.
      t_mushy = ice_temp_from_en_s(-0.5_wp*ICE_LAT_FUS, 0.0_wp)
      call check(error, t_mushy == 0.0_wp, &
                 "fresh-water mushy plateau must map to T = 0")
      if (allocated(error)) return

      ! (c) Solid leg: T = -10 round-trips.
      en = ice_enth_from_ts(-10.0_wp, 0.0_wp)
      t_back = ice_temp_from_en_s(en, 0.0_wp)
      tol = 1.0e-14_wp*(abs(-10.0_wp) + 1.0_wp)
      call check(error, abs(t_back - (-10.0_wp)) < tol, &
                 "fresh-water solid leg must round-trip")
   end subroutine test_fresh_water_sub_branches

   subroutine test_monotonic_and_latent_jump(error)
      !! At S = 35: ice_enth_from_ts is strictly increasing in T over
      !! [-40, +5]. The enthalpy deficit at T = 2*t_fr (brine fraction
      !! 1/2) carries the ~0.5*LAT_FUS mushy latent signature — a loose
      !! physics sanity check, not an identity.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: s = 35.0_wp
      integer, parameter :: NT = 90
      real(wp) :: t_lo, t_hi, t_cur, en_prev, en_cur, t_fr, deficit, tol
      integer :: it

      t_lo = -40.0_wp
      t_hi = 5.0_wp
      en_prev = ice_enth_from_ts(t_lo, s)
      do it = 2, NT
         t_cur = t_lo + (t_hi - t_lo)*real(it - 1, wp)/real(NT - 1, wp)
         en_cur = ice_enth_from_ts(t_cur, s)
         call check(error, en_cur > en_prev, &
                    "ice_enth_from_ts must be strictly increasing in T")
         if (allocated(error)) return
         en_prev = en_cur
      end do

      t_fr = ice_t_freeze(s)
      deficit = ice_enthalpy_liquid_freeze(s) - ice_enth_from_ts(2.0_wp*t_fr, s)
      tol = 0.25_wp*0.5_wp*ICE_LAT_FUS
      call check(error, abs(deficit - 0.5_wp*ICE_LAT_FUS) <= tol, &
                 "mushy latent deficit at T = 2*t_fr must be within 25% of 0.5*LAT_FUS")
   end subroutine test_monotonic_and_latent_jump

   subroutine test_elemental_dispatch(error)
      !! Whole-array elemental calls must match per-element scalar
      !! calls bitwise, guarding the `elemental` attribute.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NV = 5
      real(wp) :: tv(NV), sv(NV), ev(NV), tb(NV)
      real(wp) :: scalar_en, scalar_t
      integer :: i

      tv = [-10.0_wp, -1.0_wp, 0.0_wp, 2.0_wp, 10.0_wp]
      sv = [35.0_wp, 5.0_wp, 0.0_wp, 20.0_wp, 0.0_wp]

      ev = ice_enth_from_ts(tv, sv)
      do i = 1, NV
         scalar_en = ice_enth_from_ts(tv(i), sv(i))
         call check(error, ev(i) == scalar_en, &
                    "elemental ice_enth_from_ts must match scalar dispatch bitwise")
         if (allocated(error)) return
      end do

      tb = ice_temp_from_en_s(ev, sv)
      do i = 1, NV
         scalar_t = ice_temp_from_en_s(ev(i), sv(i))
         call check(error, tb(i) == scalar_t, &
                    "elemental ice_temp_from_en_s must match scalar dispatch bitwise")
         if (allocated(error)) return
      end do
   end subroutine test_elemental_dispatch

end module test_ocean_ice_enthalpy
