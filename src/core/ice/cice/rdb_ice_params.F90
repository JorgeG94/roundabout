!! Sea-ice constants (two selectable sets), parameter types, tracer index model,
!! salinity profile, freezing point and the enthalpy <-> temperature conversions.
module rdb_ice_params
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), modules icepack_parameters,
   ! icepack_therm_shared and icepack_tracers.
   ! Translated to roundabout's do-concurrent/GPU form and modified by the
   ! roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   ! Original authors (from the upstream headers): Elizabeth C. Hunke (LANL),
   ! C. M. Bitz (UW), William H. Lipscomb (LANL).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); Maykut and Untersteiner
   ! (1971); Hunke et al. (2015, CICE documentation, constants table).
   !
   !! **Constants.**  `ice_constants(set)` returns the physical constants of
   !! record as a plain value type, built once on the host (no default
   !! initialisers, so it is safe to pass `intent(in)` into device code).
   !! `ICE_CONSTANTS_ICEPACK` is the upstream code default.
   !! `ICE_CONSTANTS_CESM` is the set a CESM-coupled build installs through
   !! the CESM shared constants (ocean heat capacity 3996, `cp_ice` 2117.27,
   !! `Lfresh` 3.337e5, `zvir` 0.608, `snowpatch` 0.005, `Tocnfrz = -34*0.054`,
   !! ...).  Which of the two a given OM3 configuration actually runs is an
   !! open maintainer item: the default is the upstream set (interim choice).
   !!
   !! **Kernels take constants from `p%c`** (the params bundle), never from
   !! `ice_constants(set)`: configuration writes `ksno` and `snowpatch` into
   !! `p%c`, so a fresh table lookup would silently drop those overrides.
   !!
   !! **Tracer index model.**  The packed column is
   !! `Tsfc | qice(1:nilyr) | qsno(1:nslyr) | sice(1:nilyr)`, top-down in the
   !! ice (`k = 1` top).  Indices are pure functions of the layer counts (no
   !! module-global mutable indices).  `ice_trcr_depend` gives each tracer's
   !! weighting: 0 area, 1 ice volume, 2 snow volume.
   !!
   !! **Thermodynamics.**  Ice enthalpy `q` [J/m3 of ice, negative] and
   !! temperature `T` [degC] are related through the layer's brine-pocket
   !! melting temperature `Tmlt = -depressT*S` (Bitz and Lipscomb 1999 eqs.
   !! 1-3, with the `cp_ocn*Tmlt` offset of the reference state).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_safe_math, only: safe_cos, safe_pow
   implicit none
   private

   public :: ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   public :: ICE_TFRZ_LINEAR_SALT, ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8
   public :: ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71
   public :: ice_constants_t, ice_thermo_params_t, ice_itd_params_t, ice_params_t
   public :: ice_sw_params_t
   public :: ice_ridge_params_t, ice_ridge_params_default, ice_ridge_params_valid
   public :: ice_constants, ice_params_default
   public :: ice_nt_Tsfc, ice_nt_qice, ice_nt_qsno, ice_nt_sice, ice_ntrcr
   public :: ice_trcr_depend, ice_trcr_base, ice_intensive
   public :: ice_salinity_at, ice_salinity_profile, ice_sea_freezing_temperature
   public :: ice_q_from_T, ice_qsno_from_T, ice_qsno_init_from_T, ice_T_from_q, ice_Tsno_from_q

   integer, parameter :: ICE_CONSTANTS_ICEPACK = 0
      !! `constants_set`: upstream Icepack code defaults.
   integer, parameter :: ICE_CONSTANTS_CESM = 1
      !! `constants_set`: the CESM-coupled override set.

   integer, parameter :: ICE_TFRZ_LINEAR_SALT = 1
      !! `tfrz_option = "linear_salt"`: `Tf = -depressT*sss`.
   integer, parameter :: ICE_TFRZ_CONSTANT = 2
      !! `tfrz_option = "constant"`: `Tf = Tocnfrz`.
   integer, parameter :: ICE_TFRZ_MINUS1P8 = 3
      !! `tfrz_option = "minus1p8"`: `Tf = -1.8`.

   integer, parameter :: ICE_CONDUCT_BUBBLY = 1
      !! `conduct = "bubbly"` (Pringle et al. 2007).
   integer, parameter :: ICE_CONDUCT_MU71 = 2
      !! `conduct = "MU71"` (Maykut and Untersteiner 1971).

   real(wp), parameter :: PI_CONST = 3.14159265358979323846_wp
   real(wp), parameter :: SAL_NSAL = 0.407_wp   !! profile exponent numerator (BL99)
   real(wp), parameter :: SAL_MSAL = 0.573_wp   !! profile exponent offset (BL99)

   type :: ice_constants_t
      !! Physical constants of record.  Values only; no default initialisers.
      real(wp) :: rhoi            !! ice density [kg/m3]
      real(wp) :: rhos            !! snow density [kg/m3]
      real(wp) :: rhosi           !! mean sea-ice density [kg/m3]
      real(wp) :: rhow            !! seawater density [kg/m3]
      real(wp) :: rhofresh        !! fresh water density [kg/m3]
      real(wp) :: cp_ice          !! specific heat of fresh ice [J/kg/K]
      real(wp) :: cp_ocn          !! specific heat of seawater [J/kg/K]
      real(wp) :: cp_air          !! specific heat of air [J/kg/K]
      real(wp) :: cp_wv           !! specific heat of water vapour [J/kg/K]
      real(wp) :: Lsub            !! latent heat of sublimation [J/kg]
      real(wp) :: Lvap            !! latent heat of vaporisation [J/kg]
      real(wp) :: Lfresh          !! latent heat of melting, Lsub - Lvap [J/kg]
      real(wp) :: depressT        !! freezing-point depression per ppt [degC/ppt]
      real(wp) :: Tffresh         !! freezing temperature of fresh water [K]
      real(wp) :: Timelt          !! melting temperature, ice top surface [degC]
      real(wp) :: Tsmelt          !! melting temperature, snow top surface [degC]
      real(wp) :: Tocnfrz         !! seawater freezing temperature, "constant" option [degC]
      real(wp) :: kice            !! thermal conductivity of fresh ice [W/m/K]
      real(wp) :: ksno            !! thermal conductivity of snow [W/m/K]
      real(wp) :: hs_min          !! minimum snow thickness for snow temperature [m]
      real(wp) :: snowpatch       !! snow-cover patchiness length [m]
      real(wp) :: saltmax         !! salinity at the ice base, BL99 profile [ppt]
      real(wp) :: min_salin       !! brine-pocket salinity threshold [ppt]
      real(wp) :: puny            !! small number
      real(wp) :: stefan_boltzmann  !! [W/m2/K4]
      real(wp) :: vonkar          !! von Karman constant
      real(wp) :: gravit          !! gravitational acceleration [m/s2]
      real(wp) :: zvir            !! rh2o/rair - 1
      real(wp) :: zref            !! reference height for stability [m]
      real(wp) :: iceruf_ocn      !! under-ice roughness [m]
      real(wp) :: itd_area_min    !! zap residual ice below this area
      real(wp) :: itd_mass_min    !! zap residual ice below this mass [kg/m2]
      real(wp) :: secday          !! seconds per day
      real(wp) :: pi              !! pi
      logical :: bad              !! .true. when `constants_set` was not a known selector
   end type ice_constants_t

   type :: ice_thermo_params_t
      !! Bitz-Lipscomb thermodynamics knobs (upstream code defaults).
      integer :: ktherm               !! 1 = BL99 (only value accepted)
      integer :: conduct              !! ICE_CONDUCT_*
      integer :: tfrz_option          !! ICE_TFRZ_*
      real(wp) :: ice_ref_salinity    !! salt-flux reference salinity [ppt]
      real(wp) :: hfrazilmin          !! minimum new-frazil thickness [m]
      real(wp) :: hi_min              !! minimum ice thickness for thermo [m]
      real(wp) :: ustar_min           !! minimum friction velocity [m/s]
      real(wp) :: floediam            !! floe diameter for lateral melt [m]
      real(wp) :: emissivity          !! surface emissivity
      logical :: update_ocn_f         !! frazil fresh/salt fluxes to the ocean
      logical :: calc_Tsfc            !! solve for the surface temperature
      logical :: l_brine              !! brine pockets active (saltmax > min_salin)
      real(wp) :: floeshape           !! floe shape factor for lateral melt (Steele 1992)
      logical :: congel_one_step      !! congelation freezing in one step (default two-step)
   end type ice_thermo_params_t

   type :: ice_itd_params_t
      !! Thickness-distribution structure knobs.
      integer :: kitd                 !! 1 = linear ITD remap
      integer :: kcatbound            !! category-bound formula (upstream default 1)
   end type ice_itd_params_t

   type :: ice_sw_params_t
      !! Shortwave albedo parameters (CCSM3 albedo; upstream code defaults).
      !! D1.3 owns this type and the `&ocean_ice_sw_nml` group; the shortwave
      !! driver PR (D7.1) EXTENDS both and adds the refusals for unsupported
      !! shortwave options.  Do not create a second type.
      real(wp) :: albicev             !! visible albedo of thick ice
      real(wp) :: albicei             !! near-infrared albedo of thick ice
      real(wp) :: albsnowv            !! visible albedo of cold snow
      real(wp) :: albsnowi            !! near-infrared albedo of cold snow
      real(wp) :: ahmax               !! thickness above which the ice albedo is constant [m]
   end type ice_sw_params_t
   type :: ice_ridge_params_t
      !! Ridging and ice-strength knobs (Thorndike et al. 1975; Rothrock 1975;
      !! Hibler 1979, 1980; Flato and Hibler 1995; Lipscomb et al. 2007).
      !! Scalars only, no default initialisers (safe `intent(in)` into device
      !! code).  Only true knobs live here: the constants the upstream source
      !! fixes are module parameters of `rdb_ice_ridge_itd`, and everything
      !! derived from them or from the physical constants (`Cp`, `1/G*`,
      !! `1/a*`, the exponential normaliser) is computed where it is used, so
      !! no derived value can go stale.
      integer :: kstrength            !! 0 = Hibler (1979), 1 = Rothrock (1975)
      integer :: krdg_partic          !! participation: 0 = Thorndike (1975), 1 = Lipscomb (2007)
      integer :: krdg_redist          !! redistribution: 0 = Hibler (1980), 1 = Lipscomb (2007)
      real(wp) :: mu_rdg              !! e-folding scale of ridged ice, krdg_redist=1 [m^0.5]
      real(wp) :: Cf                  !! frictional-dissipation factor of the strength
      real(wp) :: Pstar               !! Hibler strength constant [N/m]
      real(wp) :: Cstar               !! Hibler concentration constant
   end type ice_ridge_params_t

   type :: ice_params_t
      !! All sea-ice parameters: scalars only, no allocatables, no strings.
      !! Later PRs add their own sub-types in their own files.
      logical :: bad                  !! .true. when the structure or `constants_set` was refused
      integer :: constants_set        !! ICE_CONSTANTS_*
      integer :: ncat                 !! categories
      integer :: nilyr                !! ice layers
      integer :: nslyr                !! snow layers
      type(ice_constants_t) :: c
      type(ice_thermo_params_t) :: thermo
      type(ice_itd_params_t) :: itd
      type(ice_sw_params_t) :: sw
      type(ice_ridge_params_t) :: ridge
   end type ice_params_t

contains

   pure function ice_constants(constants_set) result(c)
      !! Constants of record for `constants_set` (`ICE_CONSTANTS_ICEPACK` or
      !! `ICE_CONSTANTS_CESM`).  Any other value sets `c%bad = .true.` (and
      !! returns the upstream values so nothing downstream sees garbage); the
      !! caller must test the flag and refuse the configuration.
      integer, intent(in) :: constants_set
      type(ice_constants_t) :: c

      ! Upstream default values (Icepack code defaults).
      c%rhoi = 917.0_wp
      c%rhos = 330.0_wp
      c%rhosi = 940.0_wp
      c%rhow = 1026.0_wp
      c%rhofresh = 1000.0_wp
      c%cp_ice = 2106.0_wp
      c%cp_ocn = 4218.0_wp
      c%cp_air = 1005.0_wp
      c%cp_wv = 1.81e3_wp
      c%Lsub = 2.835e6_wp
      c%Lvap = 2.501e6_wp
      c%depressT = 0.054_wp
      c%Tffresh = 273.15_wp
      c%Timelt = 0.0_wp
      c%Tsmelt = 0.0_wp
      c%Tocnfrz = -1.8_wp
      c%kice = 2.03_wp
      c%ksno = 0.30_wp
      c%hs_min = 1.0e-4_wp
      c%snowpatch = 0.02_wp
      c%saltmax = 3.2_wp        ! Bitz and Lipscomb (1999)
      c%min_salin = 0.1_wp
      c%puny = 1.0e-11_wp
      c%stefan_boltzmann = 5.67e-8_wp
      c%vonkar = 0.4_wp
      c%gravit = 9.80616_wp
      c%zvir = 0.606_wp
      c%zref = 10.0_wp
      c%iceruf_ocn = 0.03_wp
      c%itd_area_min = 1.0e-11_wp
      c%itd_mass_min = 1.0e-10_wp
      c%secday = 86400.0_wp
      c%pi = PI_CONST
      c%bad = constants_set /= ICE_CONSTANTS_ICEPACK .and. constants_set /= ICE_CONSTANTS_CESM

      if (constants_set == ICE_CONSTANTS_CESM) then
         ! CESM shared constants (shr_const): the overrides a CESM-coupled
         ! build installs.  Transcribed values; to be confirmed against an
         ! OM3 run log (interim choice).
         c%cp_air = 1004.64_wp
         c%cp_ice = 2117.27_wp
         c%cp_ocn = 3996.0_wp
         c%cp_wv = 1810.0_wp
         c%gravit = 9.80616_wp
         c%zvir = cesm_zvir()
         c%Lsub = 2.834700e6_wp      ! latice 3.337e5 + latvap 2.501e6
         c%Lvap = 2.501e6_wp
         c%Tocnfrz = -34.0_wp*0.054_wp
         c%snowpatch = 0.005_wp
      end if

      c%Lfresh = c%Lsub - c%Lvap
   end function ice_constants

   pure function cesm_zvir() result(z)
      !! CESM `zvir = RWV/RDAIR - 1` from the gas constants, in the same
      !! operation order as the shared constants (so it matches bit for bit).
      real(wp) :: z, rgas, rdair, rwv

      rgas = 6.02214e26_wp*1.38065e-23_wp      ! Avogadro * Boltzmann
      rdair = rgas/28.966_wp                   ! molecular weight of dry air
      rwv = rgas/18.016_wp                     ! molecular weight of water vapour
      z = rwv/rdair - 1.0_wp
   end function cesm_zvir

   pure function ice_params_default(constants_set, ncat, nilyr, nslyr) result(p)
      !! Parameter bundle with upstream code defaults for the given structure,
      !! EXCEPT `tfrz_option`: upstream's default is `mushy` (needs the mushy
      !! thermodynamics, not ported); `linear_salt` is the plan's choice for
      !! OM3 (interim).  OM3 values otherwise live in an OM3 example namelist.
      !! `p%bad` is set when `constants_set` is unknown or a layer/category
      !! count is outside `1..ICE_*_MAX`; the caller must test it.
      integer, intent(in) :: constants_set, ncat, nilyr, nslyr
      type(ice_params_t) :: p

      p%constants_set = constants_set
      p%bad = ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nilyr > ICE_NILYR_MAX &
              .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX
      p%ncat = ncat
      p%nilyr = nilyr
      p%nslyr = nslyr
      p%c = ice_constants(constants_set)
      p%bad = p%bad .or. p%c%bad

      p%thermo%ktherm = 1
      p%thermo%conduct = ICE_CONDUCT_BUBBLY
      p%thermo%tfrz_option = ICE_TFRZ_LINEAR_SALT
      p%thermo%ice_ref_salinity = 4.0_wp
      p%thermo%hfrazilmin = 0.05_wp
      p%thermo%hi_min = 0.01_wp
      p%thermo%ustar_min = 0.005_wp
      p%thermo%floediam = 300.0_wp
      p%thermo%floeshape = 0.66_wp
      p%thermo%emissivity = 0.985_wp
      p%thermo%update_ocn_f = .false.
      p%thermo%calc_Tsfc = .true.
      p%thermo%l_brine = p%c%saltmax > p%c%min_salin
      p%thermo%congel_one_step = .false.

      p%itd%kitd = 1
      p%itd%kcatbound = 1

      p%sw%albicev = 0.78_wp
      p%sw%albicei = 0.36_wp
      p%sw%albsnowv = 0.98_wp
      p%sw%albsnowi = 0.70_wp
      p%sw%ahmax = 0.3_wp
      p%ridge = ice_ridge_params_default()
   end function ice_params_default

   pure function ice_ridge_params_default() result(rp)
      !! Ridging knobs at the upstream code defaults (`kstrength = 1`,
      !! `krdg_partic = 1`, `krdg_redist = 1`, `mu_rdg = 3`, `Cf = 17`,
      !! `Pstar = 2.75e4`, `Cstar = 20`).
      type(ice_ridge_params_t) :: rp

      rp%kstrength = 1
      rp%krdg_partic = 1
      rp%krdg_redist = 1
      rp%mu_rdg = 3.0_wp
      rp%Cf = 17.0_wp
      rp%Pstar = 2.75e4_wp
      rp%Cstar = 20.0_wp
   end function ice_ridge_params_default

   pure function ice_ridge_params_valid(rp) result(ok)
      !! `.true.` when every option of `rp` is one the kernels implement
      !! (`kstrength`, `krdg_partic`, `krdg_redist` each 0 or 1).  Host check
      !! for the drivers and the configure-time refusal.
      type(ice_ridge_params_t), intent(in) :: rp
      logical :: ok

      ok = (rp%kstrength == 0 .or. rp%kstrength == 1) .and. (rp%krdg_partic == 0 .or. rp%krdg_partic == 1) &
           .and. (rp%krdg_redist == 0 .or. rp%krdg_redist == 1)
   end function ice_ridge_params_valid

   ! ------------------------------------------------------------------
   ! Tracer index model (packed column: Tsfc | qice | qsno | sice)
   ! ------------------------------------------------------------------

   pure function ice_nt_Tsfc() result(it)
      !! Index of the surface temperature in the packed column.
      integer :: it
      it = 1
   end function ice_nt_Tsfc

   pure function ice_nt_qice() result(it)
      !! Index of the first ice-enthalpy layer (layers `it .. it+nilyr-1`).
      integer :: it
      it = 2
   end function ice_nt_qice

   pure function ice_nt_qsno(nilyr) result(it)
      !! Index of the first snow-enthalpy layer.
      !$acc routine seq
      integer, intent(in) :: nilyr
      integer :: it
      it = 2 + nilyr
   end function ice_nt_qsno

   pure function ice_nt_sice(nilyr, nslyr) result(it)
      !! Index of the first ice-salinity layer.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      integer :: it
      it = 2 + nilyr + nslyr
   end function ice_nt_sice

   pure function ice_ntrcr(nilyr, nslyr) result(n)
      !! Length of the packed tracer column.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      integer :: n
      n = 1 + 2*nilyr + nslyr
   end function ice_ntrcr

   pure function ice_trcr_depend(it, nilyr, nslyr) result(d)
      !! Weighting of tracer `it`: 0 area (Tsfc), 1 ice volume (qice, sice),
      !! 2 snow volume (qsno); -1 for an index outside the packed column.
      !$acc routine seq
      integer, intent(in) :: it, nilyr, nslyr
      integer :: d

      if (it < 1 .or. it > ice_ntrcr(nilyr, nslyr)) then
         d = -1
      else if (it == 1) then
         d = 0
      else if (it < ice_nt_qsno(nilyr)) then
         d = 1
      else if (it < ice_nt_sice(nilyr, nslyr)) then
         d = 2
      else
         d = 1
      end if
   end function ice_trcr_depend

   pure function ice_trcr_base(depend) result(base)
      !! One-hot weights (area, ice volume, snow volume) for a `depend` code.
      integer, intent(in) :: depend
      real(wp) :: base(3)

      base = 0.0_wp
      if (depend >= 0 .and. depend <= 2) base(depend + 1) = 1.0_wp
   end function ice_trcr_base

   pure function ice_intensive(atr, depend, aicen, vicen, vsnon, puny, empty_value) result(tr)
      !! Intensive tracer from its extensive product `atr` (area-, ice-volume-
      !! or snow-volume-weighted per `depend`).  An area tracer on a category
      !! with `aicen <= puny` returns `empty_value` (the freezing point for
      !! `Tsfc`, 0 for other area tracers).  Volume tracers keep upstream's
      !! test, a strictly positive divisor, with no `puny` threshold: below
      !! it the result is 0.
      !$acc routine seq
      real(wp), intent(in) :: atr, aicen, vicen, vsnon, puny, empty_value
      integer, intent(in) :: depend
      real(wp) :: tr
      real(wp) :: div

      if (depend == 0) then
         if (aicen > puny) then
            tr = atr/aicen
         else
            tr = empty_value
         end if
      else
         if (depend == 1) then
            div = vicen
         else
            div = vsnon
         end if
         if (div > 0.0_wp) then
            tr = atr/div
         else
            tr = 0.0_wp
         end if
      end if
   end function ice_intensive

   ! ------------------------------------------------------------------
   ! Salinity profile and freezing point
   ! ------------------------------------------------------------------

   elemental function ice_salinity_at(zn, saltmax) result(s)
      !! Bitz and Lipscomb (1999) prescribed salinity at normalised depth
      !! `zn` in (0,1): `(saltmax/2)*(1 - cos(pi*zn**(0.407/(0.573+zn))))`.
      !! Host/setup evaluation through the reproducible safe-math library.
      real(wp), intent(in) :: zn, saltmax
      real(wp) :: s

      s = 0.5_wp*saltmax*(1.0_wp - safe_cos(PI_CONST*safe_pow(zn, SAL_NSAL/(SAL_MSAL + zn))))
   end function ice_salinity_at

   pure subroutine ice_salinity_profile(nilyr, saltmax, min_salin, sprofile)
      !! Per-layer salinity `sprofile(1:nilyr)` (floored at `min_salin`,
      !! midpoints `zn = (k-0.5)/nilyr`, k=1 top) and `sprofile(nilyr+1) =
      !! saltmax` (the base value).  With `saltmax <= min_salin` brine pockets
      !! are off and the whole profile is zero.  Constant per configuration:
      !! evaluate once at setup, never per cell.
      integer, intent(in) :: nilyr
      real(wp), intent(in) :: saltmax, min_salin
      real(wp), intent(out) :: sprofile(nilyr + 1)
      integer :: k
      real(wp) :: zn

      if (saltmax > min_salin) then
         do k = 1, nilyr
            zn = (real(k, wp) - 0.5_wp)/real(nilyr, wp)
            sprofile(k) = max(ice_salinity_at(zn, saltmax), min_salin)
         end do
         sprofile(nilyr + 1) = saltmax
      else
         sprofile = 0.0_wp
      end if
   end subroutine ice_salinity_profile

   pure subroutine ice_sea_freezing_temperature(sss, tfrz_option, depressT, Tocnfrz, Tf, bad)
      !! Ocean freezing temperature [degC] for `tfrz_option`
      !! (`ICE_TFRZ_LINEAR_SALT`: `-depressT*sss`; `_CONSTANT`: `Tocnfrz`;
      !! `_MINUS1P8`: -1.8).  An unknown option sets `bad = .true.` and
      !! `Tf = 0`; it never guesses a form.
      !$acc routine seq
      real(wp), intent(in) :: sss, depressT, Tocnfrz
      integer, intent(in) :: tfrz_option
      real(wp), intent(out) :: Tf
      logical, intent(out) :: bad

      bad = .false.
      select case (tfrz_option)
      case (ICE_TFRZ_LINEAR_SALT)
         Tf = -depressT*sss
      case (ICE_TFRZ_CONSTANT)
         Tf = Tocnfrz
      case (ICE_TFRZ_MINUS1P8)
         Tf = -1.8_wp
      case default
         Tf = 0.0_wp
         bad = .true.
      end select
   end subroutine ice_sea_freezing_temperature

   ! ------------------------------------------------------------------
   ! Enthalpy <-> temperature (BL99)
   ! ------------------------------------------------------------------

   pure function ice_q_from_T(T, S, c, l_brine) result(q)
      !! Ice layer enthalpy [J/m3 of ice] (negative) from temperature `T`
      !! [degC] and salinity `S` [ppt] (Bitz and Lipscomb 1999):
      !! `q = -rhoi*(cp_ice*(Tmlt-T) + Lfresh*(1-Tmlt/T) - cp_ocn*Tmlt)`,
      !! `Tmlt = -depressT*S`.  `T` must be negative where `Tmlt /= 0`; with
      !! brine off (`l_brine = .false.`) `Tmlt = 0` and the fresh-ice form
      !! `q = rhoi*(cp_ice*T - Lfresh)` results.
      !$acc routine seq
      real(wp), intent(in) :: T, S
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine
      real(wp) :: q
      real(wp) :: Tmlt, lat

      if (l_brine) then
         Tmlt = -c%depressT*S
      else
         Tmlt = 0.0_wp
      end if
      if (Tmlt /= 0.0_wp) then
         lat = c%Lfresh*(1.0_wp - Tmlt/T)
      else
         lat = c%Lfresh
      end if
      q = -c%rhoi*(c%cp_ice*(Tmlt - T) + lat - c%cp_ocn*Tmlt)
   end function ice_q_from_T

   pure subroutine ice_T_from_q(q, S, c, l_brine, T, bad)
      !! Inverse of `ice_q_from_T`: the root of the BL99 quadratic
      !! `cp_ice*T**2 + b*T + Lfresh*Tmlt = 0` with
      !! `b = (cp_ocn-cp_ice)*Tmlt - q/rhoi - Lfresh`, taking the root below
      !! `Tmlt`: `T = min((-b - sqrt(b**2 - 4*cp_ice*Lfresh*Tmlt))/(2*cp_ice), Tmlt)`.
      !! Fresh ice (`l_brine = .false.`): `T = (Lfresh + q/rhoi)/cp_ice`.
      !! `bad = .true.` flags a negative discriminant or non-finite input
      !! (the caller ORs an error bit); `T` is then 0.  Never aborts.
      !$acc routine seq
      real(wp), intent(in) :: q, S
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine
      real(wp), intent(out) :: T
      logical, intent(out) :: bad
      real(wp) :: Tmlt, aa, bb, cc, disc

      bad = .false.
      T = 0.0_wp
      if (.not. (ieee_is_finite(q) .and. ieee_is_finite(S))) then
         bad = .true.
         return
      end if
      if (l_brine) then
         Tmlt = -c%depressT*S
         aa = c%cp_ice
         bb = (c%cp_ocn - c%cp_ice)*Tmlt - q/c%rhoi - c%Lfresh
         cc = c%Lfresh*Tmlt
         disc = bb*bb - 4.0_wp*aa*cc
         if (disc < 0.0_wp) then
            bad = .true.
         else
            T = min((-bb - sqrt(disc))/(2.0_wp*aa), Tmlt)
         end if
      else
         T = (c%Lfresh + q/c%rhoi)/c%cp_ice
      end if
   end subroutine ice_T_from_q

   pure function ice_qsno_init_from_T(T, c) result(q)
      !! Snow enthalpy [J/m3 of snow] in the INITIALISATION form, which clamps
      !! the temperature at the melting point: `q = -rhos*(Lfresh -
      !! cp_ice*min(0,T))`.  Use only to seed a state from a temperature; the
      !! thermodynamic update uses `ice_qsno_from_T` (no clamp).
      !$acc routine seq
      real(wp), intent(in) :: T
      type(ice_constants_t), intent(in) :: c
      real(wp) :: q

      q = -c%rhos*(c%Lfresh - c%cp_ice*min(0.0_wp, T))
   end function ice_qsno_init_from_T

   pure function ice_qsno_from_T(T, c) result(q)
      !! Snow enthalpy [J/m3 of snow] in the thermodynamic-solver form (Bitz
      !! and Lipscomb 1999), NO clamp: `q = -rhos*(Lfresh - cp_ice*T)`.  The
      !! snow temperature is already capped at the melting point by the
      !! solver, so a clamp here would hide a bug.
      !$acc routine seq
      real(wp), intent(in) :: T
      type(ice_constants_t), intent(in) :: c
      real(wp) :: q

      q = -c%rhos*(c%Lfresh - c%cp_ice*T)
   end function ice_qsno_from_T

   pure function ice_Tsno_from_q(q, c) result(T)
      !! Snow temperature [degC] from enthalpy: `T = (Lfresh + q/rhos)/cp_ice`.
      !$acc routine seq
      real(wp), intent(in) :: q
      type(ice_constants_t), intent(in) :: c
      real(wp) :: T

      T = (c%Lfresh + q/c%rhos)/c%cp_ice
   end function ice_Tsno_from_q

end module rdb_ice_params
