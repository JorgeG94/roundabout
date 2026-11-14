!! Uniform-slab seeding of the CICE-port ice state (roundabout-original).
module rdb_ice_seed
   !! `ice_seed_uniform_impl` fills the per-category state arrays with one
   !! uniform slab: every cell (ghosts included) gets concentration `aice`,
   !! mean ice thickness `hi` and mean snow depth `hs`.  It is SETUP code that
   !! runs on the host BEFORE the device map (or is followed by an
   !! `!$acc update device` of the component arrays by the slot wrapper): it
   !! executes once, evaluates a handful of column profiles and broadcasts
   !! them, so a kernel would only add a launch and a mapping requirement.
   !!
   !! Placement.  The slab is one category, the one whose bounds contain `hi`:
   !! `hin_max(n-1) <= hi < hin_max(n)`.  There `aicen = aice`, `vicen =
   !! aice*hi`, `vsnon = aice*hs` (volumes per cell area), so the requested
   !! volume `aice*hi` is conserved exactly by the placement.  All other
   !! categories are empty.
   !!
   !! Profiles (Bitz and Lipscomb 1999 thermodynamics; layer `k = 1` is the TOP,
   !! the documented exception to the ocean's bottom-up order).  With
   !! `zn = (k - 1/2)/nilyr`:
   !!   * salinity `sicen(k) = max(S_prof(zn), min_salin)` (`ice_salinity_profile`;
   !!     all zero when brine pockets are off), in EVERY category, as the
   !!     thermodynamics treats it as a stored tracer;
   !!   * temperature linear in `zn` from the surface value `Ts = min(Tsfc, Timelt)`
   !!     to the ocean freezing point `Tf` at the base, capped at the layer's
   !!     melting point `Tmlt(k) = -depressT*sicen(k)`:
   !!     `T(k) = min(Ts + (Tf - Ts)*zn, Tmlt(k))`;
   !!   * `qicen(k) = ice_q_from_T(T(k), sicen(k))` and
   !!     `qsnon = ice_qsno_init_from_T(Ts)` (the initialisation form, clamped at
   !!     the melting point); `Tsfcn = Ts`.
   !! An EMPTY category carries the isothermal base column: `Ts = min(Tf, Timelt)`,
   !! so every array is finite and brine-consistent, never zero (zero enthalpy is
   !! warmer than any ice).  Only `Tsfcn = Tf` is the upstream value of an empty
   !! cell (the area tracer's empty value); the empty-category LAYER tracers
   !! (`qicen`, `sicen`, `qsnon`) are DON'T-CARE, since every one of them is
   !! volume-weighted and an empty category has zero volume: a golden comparison
   !! against a reference must mask them (compare occupied categories only).
   !! A thickness exactly on a bound, `hi == hin_max(n)`, goes to category n+1
   !! (the bounds are lower-inclusive).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX
   use rdb_ice_params, only: ice_params_t, ice_salinity_profile, ice_q_from_T, ice_qsno_init_from_T
   implicit none
   private

   public :: ice_seed_uniform_impl
   public :: ICE_SEED_OK, ICE_SEED_BAD_ARG, ICE_SEED_BAD_EXTENT

   integer, parameter :: ICE_SEED_OK = 0
      !! Seeded.
   integer, parameter :: ICE_SEED_BAD_ARG = 1
      !! `aice` outside [0,1], negative/non-finite `hi`/`hs`/`Tsfc`/`Tf`, an
      !! occupied slab with `hi = 0`, or `hi` at or above the top category bound.
      !! The arrays are untouched.
   integer, parameter :: ICE_SEED_BAD_EXTENT = 2
      !! `ncat`/`nilyr`/`nslyr` disagree with `p`, or exceed the compile-time bounds.

contains

   subroutine ice_seed_uniform_impl(nx, ny, ncat, nilyr, nslyr, p, hin_max, aice_in, hi, hs, Tsfc, Tf, &
                                    aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr)
      !! Fill the state arrays with the uniform slab described in the module
      !! header.  `ierr` is `ICE_SEED_*`; on a non-zero value nothing is written.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: aice_in, hi, hs, Tsfc, Tf
      real(wp), intent(out) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(out) :: Tsfcn(nx, ny, ncat)
      real(wp), intent(out) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(out) :: qsnon(nx, ny, ncat, nslyr)
      integer, intent(out) :: ierr

      integer :: n, k, ncell
      real(wp) :: sprof(ICE_NILYR_MAX + 1)
      real(wp) :: Ts, T, Tmlt, zn
      logical :: occupied

      ierr = ICE_SEED_OK
      if (p%ncat /= ncat .or. p%nilyr /= nilyr .or. p%nslyr /= nslyr .or. nilyr > ICE_NILYR_MAX &
          .or. nilyr < 1 .or. nslyr < 1) then
         ierr = ICE_SEED_BAD_EXTENT
         return
      end if
      if (.not. (ieee_is_finite(aice_in) .and. ieee_is_finite(hi) .and. ieee_is_finite(hs) &
                 .and. ieee_is_finite(Tsfc) .and. ieee_is_finite(Tf))) then
         ierr = ICE_SEED_BAD_ARG
         return
      end if
      if (aice_in < 0.0_wp .or. aice_in > 1.0_wp .or. hi < 0.0_wp .or. hs < 0.0_wp) then
         ierr = ICE_SEED_BAD_ARG
         return
      end if
      occupied = aice_in > 0.0_wp
      ncell = 0
      if (occupied) then
         if (hi <= 0.0_wp .or. hi >= hin_max(ncat)) then
            ierr = ICE_SEED_BAD_ARG
            return
         end if
         do n = 1, ncat
            if (hi < hin_max(n)) then
               ncell = n
               exit
            end if
         end do
      end if

      call ice_salinity_profile(nilyr, p%c%saltmax, p%c%min_salin, sprof(1:nilyr + 1))

      aicen = 0.0_wp
      vicen = 0.0_wp
      vsnon = 0.0_wp
      do n = 1, ncat
         if (n == ncell) then
            Ts = min(Tsfc, p%c%Timelt)
            aicen(:, :, n) = aice_in
            vicen(:, :, n) = aice_in*hi
            vsnon(:, :, n) = aice_in*hs
         else
            Ts = min(Tf, p%c%Timelt)
         end if
         Tsfcn(:, :, n) = Ts
         do k = 1, nilyr
            zn = (real(k, wp) - 0.5_wp)/real(nilyr, wp)
            if (n == ncell) then
               T = Ts + (Tf - Ts)*zn
            else
               T = Ts
            end if
            if (p%thermo%l_brine) then
               Tmlt = -p%c%depressT*sprof(k)
            else
               Tmlt = 0.0_wp
            end if
            T = min(T, Tmlt)
            sicen(:, :, n, k) = sprof(k)
            qicen(:, :, n, k) = ice_q_from_T(T, sprof(k), p%c, p%thermo%l_brine)
         end do
         qsnon(:, :, n, :) = ice_qsno_init_from_T(Ts, p%c)
      end do
   end subroutine ice_seed_uniform_impl

end module rdb_ice_seed
