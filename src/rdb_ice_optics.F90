!! Sea-ice shortwave optics (SIS2 port, PR 3a).
module rdb_ice_optics
   !! CSIM4 (non-delta-Eddington) branch of `ice_optics_SIS2`
   !! (SIS_optics.F90:371-409, Apache-2.0) — snow/ice albedo + the
   !! Beer's-law vertical partition of absorbed shortwave into the
   !! surface, snow, per-layer ice, and transmitted-to-ocean terms.
   !!
   !! The delta-Eddington branch (SIS_optics.F90 `do_deltaEdd`, spectral
   !! direct/diffuse albedo bands) is deliberately NOT ported — CSIM4 is
   !! SIS2's non-delta-Eddington default and the only branch this PR
   !! needs.
   !!
   !! Everything pure + `!$acc routine seq`: called from inside the
   !! per-column `do concurrent` driver in `rdb_ice_column`.
   use rdb_constants, only: wp
   use rdb_ice_enthalpy, only: ice_t_freeze
   implicit none
   private

   public :: ice_optics_csim4

   ! ---- CSIM4 optics constants (SIS_optics.F90:129-141 defaults) ----
   real(wp), parameter, public :: ICE_ALB_SNOW = 0.85_wp
      !! Cold-snow albedo — SIS2 `SNOW_ALBEDO`.
   real(wp), parameter, public :: ICE_ALB_ICE = 0.5826_wp
      !! Cold bare-ice albedo (non-slab default) — SIS2 `ICE_ALBEDO`.
   real(wp), parameter, public :: ICE_PEN_ICE = 0.3_wp
      !! Fraction of absorbed SW that penetrates below the surface skin
      !! for bare ice — SIS2 `ICE_SW_PEN_FRAC`.
   real(wp), parameter, public :: ICE_OPT_DEP_ICE = 0.67_wp
      !! E-folding optical depth of penetrating SW in ice (m) — SIS2
      !! `ICE_OPTICAL_DEPTH`.
   real(wp), parameter, public :: ICE_T_RANGE_MELT = 1.0_wp
      !! Temperature range (degC) over which the melting-albedo
      !! reduction ramps in — SIS2 `T_RANGE_MELT`.
   real(wp), parameter, public :: ICE_SNOW_PATCH = 0.02_wp
      !! Thin-snow masking depth (m) — inline `0.02` at
      !! SIS_optics.F90:373; kept as a named constant here.

contains

   pure subroutine ice_optics_csim4(nk, hs, hi, ts, sal_ice_top, &
                                    albedo, abs_sfc, abs_snow, abs_ice_lay, &
                                    abs_ocn, abs_int, pen)
      !! CSIM4 albedo + Beer's-law vertical SW partition
      !! (SIS_optics.F90:371-409). Exact port; inline literals
      !! `0.1235`/`0.075` (melt-albedo reductions) and `5.0`/`0.5`
      !! (thin-ice atan ramp) and `0.06` (thin-ice albedo floor) stay
      !! inline per SIS2 (SIS_optics.F90:378-384).
      !!
      !! Partition identity (up to ~1 `exp` round-off): `abs_sfc +
      !! abs_snow + sum(abs_ice_lay) + abs_ocn == 1`, because
      !! `opt_decay_lay**nk == exp(-hi/ICE_OPT_DEP_ICE)`.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(in) :: hs
         !! Snow thickness (m).
      real(wp), intent(in) :: hi
         !! Ice thickness (m).
      real(wp), intent(in) :: ts
         !! Skin/surface temperature stand-in (degC).
      real(wp), intent(in) :: sal_ice_top
         !! Bulk salinity of the top ice layer (PSU) — sets the
         !! melt-onset freezing temperature.
      real(wp), intent(out) :: albedo
         !! Combined snow+ice broadband albedo (nondim).
      real(wp), intent(out) :: abs_sfc
         !! Fraction of absorbed SW deposited at the surface skin.
      real(wp), intent(out) :: abs_snow
         !! Fraction of absorbed SW deposited in the snow (always 0
         !! in the CSIM4 branch — SIS2 keeps the term for symmetry
         !! with the delta-Eddington branch).
      real(wp), intent(out) :: abs_ice_lay(nk)
         !! Fraction of absorbed SW deposited per ice layer, TOP-DOWN
         !! (`abs_ice_lay(1)` = top ice layer).
      real(wp), intent(out) :: abs_ocn
         !! Fraction of absorbed SW transmitted through to the ocean.
      real(wp), intent(out) :: abs_int
         !! Fraction of absorbed SW deposited in the ice interior
         !! (`pen - sw_frac_top` after the Beer's-law drain).
      real(wp), intent(out) :: pen
         !! Fraction of the total (post-albedo) SW that penetrates
         !! below the surface skin.

      real(wp) :: as, ai, snow_cover, temp_ice_freeze, fh, melt_ramp
      real(wp) :: opt_decay_lay, sw_frac_top
      integer :: m

      as = ICE_ALB_SNOW
      ai = ICE_ALB_ICE
      snow_cover = hs/(hs + ICE_SNOW_PATCH)
      temp_ice_freeze = ice_t_freeze(sal_ice_top)

      fh = min(atan(5.0_wp*hi)/atan(5.0_wp*0.5_wp), 1.0_wp)

      if (ts + ICE_T_RANGE_MELT > temp_ice_freeze) then
         ! Reduce albedo for melting, CSIM4 0.53/0.47 vis/ir split
         ! (SIS_optics.F90:378-384).
         melt_ramp = min((ts + ICE_T_RANGE_MELT - temp_ice_freeze)/ICE_T_RANGE_MELT, 1.0_wp)
         as = as - 0.1235_wp*melt_ramp
         ai = ai - 0.075_wp*melt_ramp
      end if
      ai = fh*ai + (1.0_wp - fh)*0.06_wp

      albedo = snow_cover*as + (1.0_wp - snow_cover)*ai

      pen = (1.0_wp - snow_cover)*ICE_PEN_ICE
      opt_decay_lay = exp(-hi/(real(nk, wp)*ICE_OPT_DEP_ICE))
      abs_ocn = pen*exp(-hi/ICE_OPT_DEP_ICE)
      abs_sfc = 1.0_wp - pen
      abs_snow = 0.0_wp

      sw_frac_top = pen
      do m = 1, nk
         abs_ice_lay(m) = sw_frac_top*(1.0_wp - opt_decay_lay)
         sw_frac_top = sw_frac_top*opt_decay_lay
      end do
      abs_int = pen - sw_frac_top
   end subroutine ice_optics_csim4

end module rdb_ice_optics
