!! Sea-ice enthalpy library (SIS2 port, PR 2).
module rdb_ice_enthalpy
   !! Specific-enthalpy <-> temperature inversions for the sea-ice
   !! thermodynamic column model.
   !!
   !! * The thermodynamic prognostic of the ice model is SPECIFIC
   !!   ENTHALPY (J/kg), never temperature; these are the exact T<->E
   !!   inversions every later rung (PR 3 Winton column, PR 4
   !!   transport re-layering) consumes.
   !! * Source: ../SIS2/src/SIS2_ice_thm.F90 (Apache-2.0) — T_Freeze
   !!   (:1622), enth_from_TS (:1647), enthalpy_liquid_freeze (:1679),
   !!   enthalpy_liquid (:1691), Temp_from_En_S (:1830); constant
   !!   defaults from ice_thermo_init (:1534-1617).
   !! * THE SIMPLIFICATION (ratified in PLAN_SEA_ICE.md): Cp_brine ==
   !!   Cp_ice. SIS2's default Cp_brine = Cp_water makes Temp_from_En_S
   !!   transcendental (a T_fr*log(T_fr/T) term) and needs the
   !!   20-iteration Newton/false-position solve at SIS2_ice_thm.F90:
   !!   1876-1933; with Cp_brine = Cp_ice the log term vanishes and the
   !!   T<->E map is a closed-form quadratic, invertible EXACTLY. Only
   !!   the Cp_ice == Cp_brine branches are ported; the else branches
   !!   are deliberately NOT.
   !! * WHY THE ROUND-TRIP IS EXACT: for S > 0, T < t_fr, substituting
   !!   ice_enth_from_ts's En into ice_temp_from_en_s gives
   !!   BB = 0.5*(ICE_LAT_FUS*t_fr/T + ICE_CP_ICE*T), and the quadratic
   !!   ICE_CP_ICE*T**2 - 2*BB*T + ICE_LAT_FUS*t_fr = 0 has T as its
   !!   smaller root by construction — so
   !!   ice_temp_from_en_s(ice_enth_from_ts(T,S), S) == T to machine
   !!   precision, not merely to a solver tolerance. The fresh-water
   !!   and melted branches invert trivially.
   !! * Everything pure elemental + `!$acc routine seq`: host scalars,
   !!   whole-array elemental calls, or per-point inside a
   !!   do concurrent device kernel.
   use rdb_constants, only: wp
   implicit none
   private

   public :: ice_t_freeze, ice_enth_from_ts, ice_enthalpy_liquid_freeze, &
             ice_enthalpy_liquid, ice_temp_from_en_s

   ! SIS2 ice-thermodynamics constants (ice_thermo_init defaults,
   ! SIS2_ice_thm.F90:1534-1617), unscaled SI. Co-located with their
   ! kernel per the rdb_eos TFR_S_COEFF precedent.
   real(wp), parameter, public :: ICE_LAT_FUS = 3.34e5_wp
      !! Latent heat of fusion (J/kg) — SIS2 `LATENT_HEAT_FUSION`.
   real(wp), parameter, public :: ICE_CP_ICE = 2100.0_wp
      !! Heat capacity of fresh ice (J/kg/K) — SIS2 `CP_ICE`.
   real(wp), parameter, public :: ICE_CP_WATER = 4200.0_wp
      !! Heat capacity of seawater as carried by the ICE model
      !! (J/kg/K) — SIS2 `CP_SEAWATER`. Deliberately DISTINCT from the
      !! ocean's `SEAWATER_CP = 3992` (`rdb_ocean_surface_flux`): SIS2's
      !! ice model owns its own Cp_water = 4200 as the enthalpy
      !! zero-of-reference, and reusing 3992 here would break bit-parity
      !! with SIS2's enthalpy and mis-close the PR-3 frazil energy
      !! handshake. Do not "unify" them.
   real(wp), parameter, public :: ICE_CP_BRINE = ICE_CP_ICE
      !! Heat capacity of brine pockets (J/kg/K). THE simplification:
      !! set equal to `ICE_CP_ICE` (SIS2's documented computational-
      !! convenience option for `CP_BRINE`, not its default Cp_water) so
      !! the T<->E map closes in quadratic form with no Newton loop.
   real(wp), parameter, public :: ICE_DTF_DS = -0.054_wp
      !! Liquidus slope dT_f/dS (degC/PSU) — SIS2 `DTFREEZE_DS`. Same
      !! value as the ocean-side `TFR_S_COEFF` (rdb_eos) by
      !! construction — SIS2 parity on both sides of the seam.
      !!
      !! DELIBERATELY NOT switched by `&ocean_eos_nml tfreeze_set`.
      !! That knob selects the OCEAN-side liquidus that
      !! `eos_freezing_point` returns (the sea-surface freezing
      !! temperature the frazil / basal-flux seam works against); THIS
      !! constant is internal to the SIS2 enthalpy relation, where it
      !! fixes the brine-pocket melting temperature inside the ice and
      !! is baked into the closed-form quadratic T<->E map (and into
      !! SIS2 bit-parity). Retuning it is an ice-thermodynamics change,
      !! not an ocean-liquidus one. Consequence, documented rather than
      !! papered over: under `tfreeze_set = "isomip"` the ocean surface
      !! freezing point and the ice-internal brine liquidus disagree by
      !! ~0.03 degC. The ISOMIP+ set exists for ICE-SHELF-CAVITY work,
      !! where the sea-ice column model is normally off; running both at
      !! once is legal but means accepting that offset.
   real(wp), parameter, public :: ICE_ENTH_LIQ_0 = 0.0_wp
      !! Enthalpy of liquid fresh water at 0 degC (J/kg) — SIS2
      !! `ENTHALPY_LIQUID_0`. Zero-point of the enthalpy scale; kept in
      !! the algebra (even though it is 0) so a future non-zero
      !! reference is a one-line change.
   real(wp), parameter, public :: ICE_LIQ_LIM = 0.99_wp
      !! `LIQUID_LIMIT` (PR 3a) — max liquid fraction before the
      !! excess-enthalpy clamp (SIS2_ice_thm.F90:500-524) and the
      !! bottom-freeze `min_dEnth_freeze` floor (:1164-1188). Co-located
      !! here (not `rdb_ice_column`) so both `rdb_ice_column` and
      !! `rdb_ice_mass` can depend on this leaf module without a cycle.
   integer, parameter, public :: ICE_NK_MAX = 8
      !! Compile-time cap on per-column local-array size (kappa-shear
      !! `NZ_STACK_MAX` precedent). `nk_ice <= ICE_NK_MAX` is asserted
      !! at slot init (`ocean_sea_ice_init`, `rdb_ice_state`). Co-located
      !! on this leaf (not `rdb_ice_column`) — like `ICE_LIQ_LIM` — so
      !! `rdb_ice_mass` can size its fixed device-stack locals by it
      !! without a `rdb_ice_column` → `rdb_ice_mass` → `rdb_ice_column`
      !! cycle. `rdb_ice_column` re-exports it for existing importers.

contains

   pure elemental function ice_t_freeze(s) result(t_fr)
      !! Freezing temperature (degC) of the ice brine at bulk salinity
      !! `s` — SIS2 `T_Freeze` (SIS2_ice_thm.F90:1622), linear liquidus.
      !! RAW `s` (no max(0,s) clamp) — faithful to SIS2, which clamps
      !! only in enth_from_TS. Do not unify with ice_enth_from_ts's
      !! internal t_fr.
      !$acc routine seq
      real(wp), intent(in) :: s
         !! Ice bulk salinity (PSU).
      real(wp) :: t_fr

      t_fr = ICE_DTF_DS*s
   end function ice_t_freeze

   pure elemental function ice_enthalpy_liquid_freeze(s) result(enth)
      !! Enthalpy (J/kg) of liquid water at the freezing point for
      !! salinity `s` — SIS2 `enthalpy_liquid_freeze` (SIS2_ice_thm.F90:1679).
      !$acc routine seq
      real(wp), intent(in) :: s
         !! Ice bulk salinity (PSU).
      real(wp) :: enth

      enth = ICE_CP_WATER*(ICE_DTF_DS*s) + ICE_ENTH_LIQ_0
   end function ice_enthalpy_liquid_freeze

   pure elemental function ice_enthalpy_liquid(t, s) result(enth)
      !! Enthalpy (J/kg) of liquid water at temperature `t` — SIS2
      !! `enthalpy_liquid` (SIS2_ice_thm.F90:1691). `s` is unused in
      !! this linear form; the argument is kept for SIS2 call-site
      !! parity (`enthalpy_liquid(T, S, ITV)`).
      !$acc routine seq
      real(wp), intent(in) :: t
         !! Water temperature (degC).
      real(wp), intent(in) :: s
         !! Ice bulk salinity (PSU) — unused, kept for call-site parity.
      real(wp) :: enth

      enth = ICE_ENTH_LIQ_0 + ICE_CP_WATER*t
   end function ice_enthalpy_liquid

   pure elemental function ice_enth_from_ts(t, s) result(enth)
      !! Ice specific enthalpy (J/kg) from temperature + bulk salinity —
      !! SIS2 `enth_from_TS` (SIS2_ice_thm.F90:1647). NOTE: here (and
      !! only here) the freezing point uses max(0, s), per SIS2.
      !! The SIS2 `else` branch at :1667-1671 (Cp_brine /= Cp_ice, the
      !! T_fr*log(T_fr/T) form) is deliberately NOT ported — with
      !! ICE_CP_BRINE == ICE_CP_ICE it is unreachable, and dropping it
      !! keeps the map closed-form-invertible (see module docstring).
      !$acc routine seq
      real(wp), intent(in) :: t
         !! Ice temperature (degC).
      real(wp), intent(in) :: s
         !! Ice bulk salinity (PSU).
      real(wp) :: enth

      real(wp) :: t_fr

      t_fr = ICE_DTF_DS*max(0.0_wp, s)

      if (s == 0.0_wp .and. t <= 0.0_wp) then
         ! Fresh water at/below freezing is assumed all ice, due to the
         ! degeneracy in inverting temperature for enthalpy (SIS2 note).
         enth = (ICE_ENTH_LIQ_0 - ICE_LAT_FUS) + ICE_CP_ICE*t
      else if (t >= t_fr) then
         ! Already melted: just the sensible heat relative to 0 degC.
         enth = ICE_ENTH_LIQ_0 + ICE_CP_WATER*t
      else
         ! Cp_ice == Cp_brine closed form (SIS2 :1664-1666).
         enth = (ICE_ENTH_LIQ_0 - ICE_LAT_FUS*(1.0_wp - t_fr/t)) + &
                (ICE_CP_ICE*t + (ICE_CP_WATER - ICE_CP_ICE)*t_fr)
      end if
   end function ice_enth_from_ts

   pure elemental function ice_temp_from_en_s(en, s) result(t)
      !! Ice temperature (degC) from specific enthalpy + bulk salinity —
      !! SIS2 `Temp_from_En_S` (SIS2_ice_thm.F90:1830), Cp_ice ==
      !! Cp_brine path only: the quadratic
      !!   ICE_CP_ICE*T**2 - 2*BB*T + ICE_LAT_FUS*t_fr = 0
      !! solved for its smaller root. The SIS2 Newton/false-position
      !! refinement (:1876-1933) is only needed when Cp_brine /= Cp_ice
      !! and is deliberately NOT ported. Uses RAW `s` for t_fr (no max),
      !! faithful to SIS2.
      !$acc routine seq
      real(wp), intent(in) :: en
         !! Ice specific enthalpy (J/kg).
      real(wp), intent(in) :: s
         !! Ice bulk salinity (PSU).
      real(wp) :: t

      real(wp) :: t_fr, en_j, bb

      t_fr = ICE_DTF_DS*s
      en_j = en - ICE_ENTH_LIQ_0

      if (s <= 0.0_wp) then
         ! Step function for fresh water: liquid / mushy plateau / solid.
         if (en_j >= 0.0_wp) then
            t = en_j/ICE_CP_WATER
         else if (en_j >= -ICE_LAT_FUS) then
            t = 0.0_wp
         else
            t = (en_j + ICE_LAT_FUS)/ICE_CP_ICE
         end if
      else if (en_j >= t_fr*ICE_CP_WATER) then
         ! Completely melted layer.
         t = en_j/ICE_CP_WATER
      else
         ! Closed-form quadratic root (SIS2 :1868-1874).
         bb = 0.5_wp*((en_j - t_fr*(ICE_CP_WATER - ICE_CP_ICE)) + ICE_LAT_FUS)
         t = (bb - sqrt(bb*bb - t_fr*ICE_CP_ICE*ICE_LAT_FUS))/ICE_CP_ICE
      end if
   end function ice_temp_from_en_s

end module rdb_ice_enthalpy
