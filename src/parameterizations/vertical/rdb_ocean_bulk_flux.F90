!! Large & Yeager bulk air-sea flux formulae — the CORE/OMIP algorithm
!! (Large & Yeager 2004, NCAR Tech Note TN-460; Large & Yeager 2009,
!! Climate Dynamics 33, "The global climatology of an interannually
!! varying air-sea flux data set") used by ACCESS-OM3's JRA55-do
!! forcing path.  C3 of the OM3 forcing plan: a GPU-safe, `pure`
!! column kernel plus a thin `do concurrent` driver over the 2D
!! model grid.
!!
!! **Scope of this module (C3).**  This ships the PHYSICS KERNEL only.
!! It does NOT read JRA55-do files, does NOT regrid, and does NOT
!! define the `atm_state_t` atmospheric-state provider — those are
!! C2b/C2c of the same plan.  `bulk_flux_driver_2d` is written against
!! plain explicit-shape `(nx,ny)` arrays named after the fields
!! `atm_state_t` will eventually carry (`tas`, `huss`, `uas`, `vas`,
!! `psl`), so that once C2c lands, the driver call site is a pointer
!! substitution, not a rewrite.  See the module-level "Seam for C2c"
!! note at the end of this file.
!!
!! ## Algorithm
!!
!! Inputs: 10 m wind (relative to the surface current), air
!! temperature + specific humidity at a (possibly different) reference
!! height — JRA55-do co-locates both at 2 m — sea-level pressure, and
!! SST.  The column kernel:
!!
!!   1. builds the saturation specific humidity over seawater at SST
!!      (`bulk_flux_qsat_seawater`, the 0.98-reduced Large & Pond 1982
!!      formula used throughout CORE/OMIP);
!!   2. forms a first-guess neutral-at-10m drag coefficient
!!      (`bulk_flux_cdn10_neutral`, Large & Yeager 2009 eq. 11, with the
!!      2009 high-wind correction term and a constant cap above 33 m/s)
!!      and the companion neutral Stanton/Dalton numbers;
!!   3. iterates `n_iter` times (OM3: `flux_max_iteration = 5`) over the
!!      Monin-Obukhov similarity correction: friction/temperature/
!!      humidity turbulent scales, the buoyancy parameter, the stability
!!      parameter `zeta` at each of the wind and temperature/humidity
!!      reference heights, the Businger-Dyer `psi_m`/`psi_h` stability
!!      functions, the 10 m-equivalent wind speed, the reference-height-
!!      adjusted air temperature/humidity, and refreshed transfer
!!      coefficients;
!!   4. forms the final stress, sensible, latent and evaporative fluxes
!!      from the converged coefficients.
!!
!! The iteration count is FIXED (no early-exit convergence test): a
!! data-dependent iteration count is GPU-hostile (divergent warps) and
!! the project convention is a namelist-exposed iteration count, not an
!! adaptive tolerance.
!!
!! ## Sign conventions — READ BEFORE WIRING A CONSUMER
!!
!! This kernel follows the meteorological/CORE convention: `sensible`
!! and `latent` are POSITIVE when they warm the ocean (down into the
!! water, matching `rdb_ocean_surface_flux`'s `q_sens`/`q_lat`
!! component fields — no sign flip needed there).  `evap` is POSITIVE
!! for net evaporation (mass leaving the ocean upward) — this is the
!! OPPOSITE of the `ocean_surface_flux_t%evap` component's MOM6
!! convention (`<= 0` = mass leaving); a caller writing into that
!! field MUST negate it (`sf%evap(i,j) = -evap`).  `taux`/`tauy` are
!! positive in the direction of the wind-minus-current vector, matching
!! `ocean_surface_stress_t%tau_x/tau_y`'s sign convention directly.
!!
!! ## GPU
!!
!! Every procedure is `pure` with an `!$acc routine seq` directive so it
!! can be called from inside a `do concurrent` on the device build;
!! `bulk_flux_driver_2d` takes explicit-shape `(nx,ny)` dummies (no
!! assumed-shape) per the project's `do concurrent` convention.
!!
!! ## Seam for C2c
!!
!! `bulk_flux_driver_2d` produces CELL-CENTRED `taux_cell`/`tauy_cell`
!! (T-point, like `uas`/`vas`), plus `q_sens`/`q_lat`/`evap_massflux`
!! already in `ocean_surface_flux_t`'s sign convention.  Wiring this
!! into the live component set needs, at minimum: (a) the `atm_state_t`
!! provider (C2c) to supply `tas`/`huss`/`uas`/`vas`/`psl` on the model
!! grid; (b) averaging `taux_cell`/`tauy_cell` onto the C-grid
!! `tau_x`/`tau_y` FACES (`ocean_surface_stress_t`, east/north faces,
!! `(nx+1,ny)`/`(nx,ny+1)`) the same way `ocean_surfstress_derived_impl`
!! goes the OTHER direction (face -> centre) — a two-line face-average
!! kernel, not written here because there is no live face-shaped
!! consumer to test it against yet; (c) a real `&ocean_bulk_flux_nml`
!! namelist group threaded through `rdb_config.F90`/`validate_config`
!! (today `bulk_flux_config_t` below is a free-standing config type,
!! default `enable = .false.`, NOT registered on `cfg%ocean` — adding
!! that registration is a `rdb_config.F90`-sized seam in its own right,
!! out of scope per the C3 work-package note "if wiring needs more than
!! a small seam, stop at the kernel + driver"); (d) an engine call site
!! (`rdb_ocean_engine.F90`) that calls the driver once per atmospheric
!! bracket and writes `evap_massflux` into `sf%evap` with the sign
!! FLIPPED (see the sign-convention note above) and `q_sens`/`q_lat`
!! straight through, then re-runs `ocean_surface_flux_assemble`.
module rdb_ocean_bulk_flux
   use rdb_constants, only: wp, GRAVITY
   implicit none
   private

   public :: bulk_flux_config_t
   public :: bulk_flux_config_n_iter_ok
   public :: bulk_flux_qsat_seawater
   public :: bulk_flux_cdn10_neutral
   public :: bulk_flux_longwave_up
   public :: bulk_flux_column
   public :: bulk_flux_driver_2d

   ! ---- Large & Yeager (2004, 2009) / Large & Pond (1982) constants ----
   ! Fixed physical constants of the algorithm (NOT namelist knobs).
   real(wp), parameter, public :: BULK_VON_KARMAN = 0.4_wp
      !! von Karman constant.
   real(wp), parameter, public :: BULK_GAS_CONST_AIR = 287.04_wp
      !! Specific gas constant of dry air, J/(kg*K).
   real(wp), parameter, public :: BULK_CP_AIR = 1000.5_wp
      !! Specific heat capacity of (moist) air at constant pressure,
      !! J/(kg*K) — the Large & Yeager CORE/OMIP constant value (as
      !! opposed to Gill 1982's temperature/humidity-dependent `cp`).
   real(wp), parameter, public :: BULK_LATENT_VAPORIZATION = 2.5e6_wp
      !! Latent heat of vaporization, J/kg (CORE/OMIP constant value).
   real(wp), parameter, public :: BULK_TVQ_AIR = 0.6078_wp
      !! Virtual-temperature humidity factor `1/eps_air - 1`, with
      !! `eps_air = M_water/M_air = 18.016/28.966`.
   real(wp), parameter, public :: BULK_WIND_MIN = 0.3_wp
      !! Minimum relative wind speed (m/s) used inside the transfer
      !! coefficients and the ustar/tstar/qstar scales — the empirical
      !! free-convection floor: at dead calm, turbulent exchange does
      !! not vanish (buoyancy-driven convective cells still ventilate
      !! the interface), so the bulk formulae are not evaluated at a
      !! true zero wind speed.  Does NOT floor the wind-stress VECTOR
      !! components themselves (see `bulk_flux_column`'s zero-wind
      !! test case).
   real(wp), parameter, public :: BULK_T0_KELVIN = 273.15_wp
      !! Celsius -> Kelvin offset.
   real(wp), parameter, public :: BULK_STEFAN_BOLTZMANN = 5.670374419e-8_wp
      !! Stefan-Boltzmann constant, W/(m^2*K^4) (CODATA 2018).
   real(wp), parameter, public :: BULK_EMISSIVITY_SEAWATER = 0.97_wp
      !! Longwave emissivity of seawater.  Commonly cited value for
      !! bulk air-sea flux formulations (e.g. Konda et al. 1994); kept
      !! as a named default, overridable by the caller.
   real(wp), parameter :: BULK_CDN10_HIGH_WIND = 2.34e-3_wp
      !! Constant neutral drag coefficient above the 33 m/s cap
      !! (Large & Yeager 2009 eq. 11b).
   real(wp), parameter :: BULK_CDN10_CAP_WIND = 33.0_wp
      !! Wind speed (m/s) above which `bulk_flux_cdn10_neutral` switches
      !! to the constant high-wind value.
   real(wp), parameter :: BULK_Z_REF10 = 10.0_wp
      !! Fixed 10 m reference height used throughout the Large & Yeager
      !! (2004/2009) neutral-coefficient fits and the height-adjustment
      !! log terms -- distinct from `bulk_flux_config_t%z_wind`, which
      !! is the ACTUAL measurement height of the wind input.
   real(wp), parameter :: BULK_QSAT_SALINITY_FACTOR = 0.98_wp
      !! Large & Pond (1982) salinity-depression factor applied to the
      !! fresh-water saturation vapour pressure to get the value over
      !! seawater.
   real(wp), parameter :: BULK_QSAT_COEF_A = 640380.0_wp
      !! Saturation-humidity prefactor, kg/m^3 (Large & Pond 1982 /
      !! CORE-OMIP `bulk_flux_qsat_seawater` fit).
   real(wp), parameter :: BULK_QSAT_COEF_B = 5107.4_wp
      !! Saturation-humidity Clausius-Clapeyron-like exponent constant,
      !! K (Large & Pond 1982 / CORE-OMIP `bulk_flux_qsat_seawater` fit).
   real(wp), parameter :: BULK_CDN10_FIT_A = 2.7_wp
      !! Large & Yeager (2009) eq. 11a neutral-drag polynomial: 1/u term
      !! coefficient.
   real(wp), parameter :: BULK_CDN10_FIT_B = 0.142_wp
      !! Large & Yeager (2009) eq. 11a neutral-drag polynomial: constant
      !! term.
   real(wp), parameter :: BULK_CDN10_FIT_C = 0.0764_wp
      !! Large & Yeager (2009) eq. 11a neutral-drag polynomial: linear
      !! (`u`) term coefficient.
   real(wp), parameter :: BULK_CDN10_FIT_D = 3.14807e-10_wp
      !! Large & Yeager (2009) eq. 11a neutral-drag polynomial: `u**6`
      !! term coefficient.
   real(wp), parameter :: BULK_CDN10_FIT_SCALE = 1.0e-3_wp
      !! Large & Yeager (2009) eq. 11a neutral-drag polynomial: overall
      !! 1e-3 scale.
   real(wp), parameter :: BULK_CEN10_COEF = 34.6e-3_wp
      !! Large & Yeager (2004) neutral Dalton-number (moisture transfer)
      !! coefficient, `cen10 = BULK_CEN10_COEF * sqrt(cdn10)`.
   real(wp), parameter :: BULK_CTN10_STABLE_COEF = 18.0e-3_wp
      !! Large & Yeager (2004) neutral Stanton-number (heat transfer)
      !! coefficient on the STABLE branch (`dtemp <= 0`).
   real(wp), parameter :: BULK_CTN10_UNSTABLE_COEF = 32.7e-3_wp
      !! Large & Yeager (2004) neutral Stanton-number (heat transfer)
      !! coefficient on the UNSTABLE branch (`dtemp > 0`).

   type :: bulk_flux_config_t
      !! Free-standing configuration bundle for the bulk-flux kernel.
      !! NOT currently registered on `cfg%ocean` (see the module's "Seam
      !! for C2c" note) — a standalone type so a future `&ocean_bulk_
      !! flux_nml` reader has an obvious target, and so a caller that
      !! does wire it up today (e.g. a test or a bench) has one place to
      !! hold the knobs rather than a scatter of bare arguments.
      logical :: enable = .false.
         !! Master switch.  `.false.` (default): nothing in this module
         !! is reachable from the engine — bit-identical to a pre-C3
         !! build by construction (no call site exists yet).
      integer :: n_iter = 5
         !! Monin-Obukhov iteration count.  OM3's CMEPS/CDEPS coupler
         !! configuration value is `flux_max_iteration = 5`.
      real(wp) :: z_wind = 10.0_wp
         !! Wind reference height (m).  JRA55-do: 10 m.
      real(wp) :: z_ta = 2.0_wp
         !! Air temperature / specific humidity reference height (m),
         !! assumed CO-LOCATED (JRA55-do measures both at 2 m).
      real(wp) :: emissivity = BULK_EMISSIVITY_SEAWATER
         !! Longwave emissivity used by `bulk_flux_longwave_up`.
   end type bulk_flux_config_t

contains

   pure function bulk_flux_config_n_iter_ok(cfg) result(ok)
      !! Validates `bulk_flux_config_t%n_iter`: the Monin-Obukhov
      !! iteration count must be `>= 1` (`bulk_flux_column`'s
      !! unchecked precondition -- it runs the fixed-point loop at
      !! least once).  No `&ocean_bulk_flux_nml` reader / `validate_
      !! config` call site exists yet (see the module's "Seam for
      !! C2c" note); this is the check that wiring must call, fail
      !! loud on `.false.`, before ever passing `cfg%n_iter` into
      !! `bulk_flux_column` / `bulk_flux_driver_2d`.  A caller that
      !! builds `cfg` directly (a test or a bench) should call this
      !! too rather than relying on the type's `n_iter = 5` default.
      type(bulk_flux_config_t), intent(in) :: cfg
      logical :: ok
      ok = (cfg%n_iter >= 1)
   end function bulk_flux_config_n_iter_ok

   pure function bulk_flux_qsat_seawater(sst_k, rho_air) result(qsat)
      !! Saturation specific humidity over seawater (kg/kg) — the
      !! Large & Pond (1982) / CORE-OMIP formula, with the 0.98 factor
      !! that reduces the fresh-water saturation value for the salinity
      !! depression of vapour pressure over seawater.
      !$omp declare target
      real(wp), intent(in) :: sst_k
         !! SST in Kelvin.
      real(wp), intent(in) :: rho_air
         !! Ambient air density, kg/m^3.
      real(wp) :: qsat
      qsat = BULK_QSAT_SALINITY_FACTOR*BULK_QSAT_COEF_A/rho_air*exp(-BULK_QSAT_COEF_B/sst_k)
   end function bulk_flux_qsat_seawater

   pure function bulk_flux_cdn10_neutral(u10n) result(cdn10)
      !! Neutral-stability, 10 m drag coefficient (Large & Yeager 2009
      !! eq. 11a-b): a wind-speed polynomial fit below 33 m/s, pinned to
      !! a constant above it.
      !$omp declare target
      real(wp), intent(in) :: u10n
         !! 10 m-equivalent neutral wind speed (m/s), already floored
         !! by the caller at `BULK_WIND_MIN`.
      real(wp) :: cdn10
      if (u10n >= BULK_CDN10_CAP_WIND) then
         cdn10 = BULK_CDN10_HIGH_WIND
      else
         cdn10 = (BULK_CDN10_FIT_A/u10n + BULK_CDN10_FIT_B + BULK_CDN10_FIT_C*u10n &
                  - BULK_CDN10_FIT_D*u10n**6)*BULK_CDN10_FIT_SCALE
      end if
   end function bulk_flux_cdn10_neutral

   pure function bulk_flux_longwave_up(sst_k, emissivity) result(lw_up)
      !! Upward longwave emission from the sea surface: the
      !! Stefan-Boltzmann law, `lw_up = emissivity * sigma * sst_k**4`,
      !! with `sigma = BULK_STEFAN_BOLTZMANN` (W/(m^2*K^4)) and a
      !! seawater `emissivity` (see `BULK_EMISSIVITY_SEAWATER`).
      !!
      !! Units: W/m^2.  Sign convention: POSITIVE, radiating AWAY from
      !! the ocean (upward) — the OPPOSITE of this module's "positive
      !! down into the ocean" convention for `sensible`/`latent`.  A
      !! caller forms a NET longwave flux in the positive-down
      !! convention by subtracting this from the downward longwave
      !! input, `q_lw = rlds - lw_up`, before writing the
      !! `ocean_surface_flux_t%q_lw` component.
      !$omp declare target
      real(wp), intent(in) :: sst_k
         !! SST in Kelvin.
      real(wp), intent(in) :: emissivity
         !! Longwave emissivity of the sea surface (dimensionless,
         !! [0, 1]) -- the result scales linearly with this factor.
      real(wp) :: lw_up
      lw_up = emissivity*BULK_STEFAN_BOLTZMANN*sst_k**4
   end function bulk_flux_longwave_up

   pure subroutine bulk_flux_psi(zeta_in, psi_m, psi_h)
      !! Businger-Dyer flux-profile stability functions (unstable branch
      !! after Paulson 1970 / the KEYPS form used throughout CORE bulk
      !! codes; linear stable branch).  Shared by all three reference
      !! heights (wind, temperature, humidity) in `bulk_flux_column`.
      !$omp declare target
      real(wp), intent(in) :: zeta_in
      real(wp), intent(out) :: psi_m, psi_h
      real(wp) :: zeta, x2, x
      zeta = sign(min(abs(zeta_in), 10.0_wp), zeta_in)
      x2 = sqrt(abs(1.0_wp - 16.0_wp*zeta))
      x2 = max(x2, 1.0_wp)
      x = sqrt(x2)
      if (zeta > 0.0_wp) then
         psi_m = -5.0_wp*zeta
         psi_h = -5.0_wp*zeta
      else
         psi_m = log((1.0_wp + 2.0_wp*x + x2)*(1.0_wp + x2)/8.0_wp) &
                 - 2.0_wp*(atan(x) - atan(1.0_wp))
         psi_h = 2.0_wp*log((1.0_wp + x2)/2.0_wp)
      end if
   end subroutine bulk_flux_psi

   pure subroutine bulk_flux_column(u_wind, v_wind, u_cur, v_cur, t_air, q_air, &
                                    slp, sst_degc, z_wind, z_ta, n_iter, &
                                    taux, tauy, sensible, latent, evap, &
                                    cd, ch, ce, ustar_out)
      !! Large & Yeager (2004, 2009) bulk air-sea flux column kernel.
      !! See the module docstring for the algorithm outline and the
      !! SIGN CONVENTIONS a caller must honour.
      !!
      !! `n_iter` must be `>= 1` (an unchecked precondition — this is a
      !! `pure` hot-loop kernel, not a validating entry point; the
      !! namelist-facing caller validates `bulk_flux_config_t%n_iter`
      !! via `bulk_flux_config_n_iter_ok` before calling here).
      !$omp declare target
      real(wp), intent(in) :: u_wind, v_wind
         !! Wind components at `z_wind` (m/s), earth-relative.
      real(wp), intent(in) :: u_cur, v_cur
         !! Ocean surface current components (m/s).
      real(wp), intent(in) :: t_air
         !! Air temperature at `z_ta` (K).
      real(wp), intent(in) :: q_air
         !! Air specific humidity at `z_ta` (kg/kg).
      real(wp), intent(in) :: slp
         !! Sea-level pressure (Pa).
      real(wp), intent(in) :: sst_degc
         !! Sea surface temperature (degC) — the ocean model's native
         !! temperature unit.
      real(wp), intent(in) :: z_wind, z_ta
         !! Reference heights (m) for the wind and the air temperature
         !! / specific humidity pair (assumed co-located).
      integer, intent(in) :: n_iter
         !! Fixed Monin-Obukhov iteration count (OM3: 5).
      real(wp), intent(out) :: taux, tauy
         !! Wind stress (N/m^2), along the wind-minus-current vector —
         !! same sign convention as `ocean_surface_stress_t%tau_x/y`.
      real(wp), intent(out) :: sensible
         !! Sensible heat flux (W/m^2), POSITIVE DOWN into the ocean.
      real(wp), intent(out) :: latent
         !! Latent heat flux (W/m^2), POSITIVE DOWN into the ocean.
      real(wp), intent(out) :: evap
         !! Evaporative mass flux (kg/m^2/s), POSITIVE = evaporation
         !! (mass leaving the ocean upward) — the OPPOSITE sign of
         !! `ocean_surface_flux_t%evap`'s MOM6 convention; negate at
         !! the call site writing into that component.
      real(wp), intent(out) :: cd, ch, ce
         !! Converged (height- and stability-corrected) drag / Stanton
         !! (heat) / Dalton (moisture) transfer coefficients.
      real(wp), intent(out) :: ustar_out
         !! Friction velocity (m/s), `sqrt(cd)*|relative wind|`.

      real(wp) :: du, dv, wv, sst_k
      real(wp) :: rho_air, qs, dtemp, dqr, tv
      real(wp) :: ta_use, qa_use
      real(wp) :: wv10n, stab
      real(wp) :: cdn10, cdn10_rt, cen10, ctn10
      real(wp) :: cdn, ctn, cen
      real(wp) :: cd_rt, ustar, tstar, qstar, bstar
      real(wp) :: zetau, zetat
      real(wp) :: psi_mu, psi_hu, psi_mt, psi_ht
      real(wp) :: xx
      integer :: n

      sst_k = sst_degc + BULK_T0_KELVIN
      du = u_wind - u_cur
      dv = v_wind - v_cur
      wv = max(sqrt(du*du + dv*dv), BULK_WIND_MIN)

      ! Ambient-level air density + the (fixed) saturation humidity at
      ! SST — both evaluated once, upstream of the height-adjustment
      ! iteration (Large & Yeager hold `qs` fixed through the loop; only
      ! the AMBIENT temperature/humidity get shifted to the wind height).
      rho_air = slp/(BULK_GAS_CONST_AIR*t_air*(1.0_wp + BULK_TVQ_AIR*q_air))
      qs = bulk_flux_qsat_seawater(sst_k, rho_air)

      dtemp = sst_k - t_air
      ! > 0 <=> SST warmer than air <=> convectively UNSTABLE.
      dqr = qs - q_air
      tv = t_air*(1.0_wp + BULK_TVQ_AIR*q_air)
      ta_use = t_air
      qa_use = q_air

      wv10n = wv
      ! `stab` selects the Stanton-number branch from the SIGN of
      ! `dtemp` alone -- never `sign()` of a value that can be a signed
      ! zero (see the git history for the cross-compiler hazard this
      ! replaced).  `merge` puts exact neutrality (`dtemp == 0`, as at
      ! `SST == T_air`) on the STABLE (18e-3) branch, identically on
      ! every compiler.
      stab = merge(0.0_wp, 1.0_wp, dtemp > 0.0_wp)
      cdn10 = bulk_flux_cdn10_neutral(wv10n)
      cdn10_rt = sqrt(cdn10)
      cen10 = BULK_CEN10_COEF*cdn10_rt
      ctn10 = (BULK_CTN10_STABLE_COEF*stab + BULK_CTN10_UNSTABLE_COEF*(1.0_wp - stab))*cdn10_rt
      cdn = cdn10
      ctn = ctn10
      cen = cen10

      do n = 1, n_iter
         cd_rt = sqrt(cdn)
         ustar = cd_rt*wv
         tstar = (ctn/cd_rt)*(-dtemp)
         qstar = (cen/cd_rt)*(-dqr)
         bstar = GRAVITY*(tstar/tv + qstar/(qa_use + 1.0_wp/BULK_TVQ_AIR))

         zetau = BULK_VON_KARMAN*bstar*z_wind/(ustar*ustar)
         call bulk_flux_psi(zetau, psi_mu, psi_hu)
         zetat = BULK_VON_KARMAN*bstar*z_ta/(ustar*ustar)
         call bulk_flux_psi(zetat, psi_mt, psi_ht)
         ! Humidity reference height assumed == z_ta (JRA55-do
         ! co-locates tas/huss at 2 m), so psi_h at the humidity height
         ! is the same psi_ht just computed -- no third stability branch.

         wv10n = wv/(1.0_wp + cdn10_rt*(log(z_wind/BULK_Z_REF10) - psi_mu)/BULK_VON_KARMAN)
         wv10n = max(wv10n, BULK_WIND_MIN)

         ta_use = t_air - tstar*(log(z_ta/z_wind) + psi_hu - psi_ht)/BULK_VON_KARMAN
         qa_use = q_air - qstar*(log(z_ta/z_wind) + psi_hu - psi_ht)/BULK_VON_KARMAN
         tv = ta_use*(1.0_wp + BULK_TVQ_AIR*qa_use)

         cdn10 = bulk_flux_cdn10_neutral(wv10n)
         cdn10_rt = sqrt(cdn10)
         cen10 = BULK_CEN10_COEF*cdn10_rt
         ! Same signed-zero hazard as the pre-loop `stab` above, this
         ! time keyed on the stability parameter `zetau` (MO convention:
         ! `zeta >= 0` <=> stable); `>=` keeps exact neutrality
         ! (`zetau == 0`, e.g. `bstar == 0`) on the STABLE branch too.
         stab = merge(1.0_wp, 0.0_wp, zetau >= 0.0_wp)
         ctn10 = (BULK_CTN10_STABLE_COEF*stab + BULK_CTN10_UNSTABLE_COEF*(1.0_wp - stab))*cdn10_rt

         xx = (log(z_wind/BULK_Z_REF10) - psi_mu)/BULK_VON_KARMAN
         cdn = cdn10/(1.0_wp + cdn10_rt*xx)**2
         xx = (log(z_wind/BULK_Z_REF10) - psi_hu)/BULK_VON_KARMAN
         ctn = ctn10/(1.0_wp + ctn10*xx/cdn10_rt)*sqrt(cdn/cdn10)
         cen = cen10/(1.0_wp + cen10*xx/cdn10_rt)*sqrt(cdn/cdn10)

         dtemp = sst_k - ta_use
         dqr = qs - qa_use
         rho_air = slp/(BULK_GAS_CONST_AIR*tv)
      end do

      cd = cdn
      ch = ctn
      ce = cen
      cd_rt = sqrt(cd)
      ustar_out = cd_rt*wv

      sensible = -rho_air*BULK_CP_AIR*ch*dtemp*wv
      latent = -rho_air*BULK_LATENT_VAPORIZATION*ce*dqr*wv
      evap = rho_air*wv*ce*dqr

      taux = rho_air*cd*wv*du
      tauy = rho_air*cd*wv*dv
   end subroutine bulk_flux_column

   pure subroutine bulk_flux_driver_2d(nx, ny, u_wind, v_wind, u_cur, v_cur, &
                                       t_air, q_air, slp, sst_degc, wet_mask, &
                                       z_wind, z_ta, n_iter, &
                                       taux_cell, tauy_cell, q_sens, q_lat, &
                                       evap_massflux, cd_out, ch_out, ce_out)
      !! Thin `do concurrent` driver — one `bulk_flux_column` call per
      !! CELL-CENTRED grid point.  Explicit-shape dummies throughout
      !! (the project's `do concurrent` convention), no `associate` over
      !! derived-type components, no namelist/engine coupling: a caller
      !! supplies plain model-grid arrays and reads plain outputs.
      !!
      !! Outputs are ALREADY in `ocean_surface_flux_t`'s sign convention
      !! for `q_sens`/`q_lat` (positive down) and for `evap_massflux`
      !! (negated here to the MOM6 `evap <= 0` convention -- see the
      !! module's sign-convention note).  `taux_cell`/`tauy_cell` are
      !! CELL-CENTRED (T-point, like the `uas`/`vas` inputs); averaging
      !! onto the C-grid stress FACES is the caller's job (see the
      !! module's "Seam for C2c" note) -- deliberately not done here, so
      !! this driver has no dependency on `ocean_surface_stress_t`'s
      !! face-shaped arrays.
      !!
      !! `wet_mask` gates every output over land (`wet_mask(i,j) <= 0`):
      !! ALL EIGHT outputs are written as explicit zero on a land cell,
      !! and `bulk_flux_column` is never called there -- land `t_air`/
      !! `q_air`/`slp`/`sst_degc`/`u_wind`/`v_wind` may legitimately be
      !! fill values (including NaN), and MULTIPLYING a column result
      !! by a zero mask does not scrub a NaN (`NaN*0.0 == NaN`), so the
      !! branch must skip the call outright rather than mask its
      !! result.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: u_wind(nx, ny), v_wind(nx, ny)
      real(wp), intent(in) :: u_cur(nx, ny), v_cur(nx, ny)
      real(wp), intent(in) :: t_air(nx, ny), q_air(nx, ny)
      real(wp), intent(in) :: slp(nx, ny), sst_degc(nx, ny)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: z_wind, z_ta
      integer, intent(in) :: n_iter
      real(wp), intent(out) :: taux_cell(nx, ny), tauy_cell(nx, ny)
      real(wp), intent(out) :: q_sens(nx, ny), q_lat(nx, ny)
      real(wp), intent(out) :: evap_massflux(nx, ny)
      real(wp), intent(out) :: cd_out(nx, ny), ch_out(nx, ny), ce_out(nx, ny)

      integer :: i, j
      real(wp) :: taux, tauy, sens, lat, evp, cd, ch, ce, ustar

      do concurrent(j=1:ny, i=1:nx) local(taux, tauy, sens, lat, evp, cd, ch, ce, ustar)
         if (wet_mask(i, j) > 0.0_wp) then
            call bulk_flux_column(u_wind(i, j), v_wind(i, j), u_cur(i, j), v_cur(i, j), &
                                  t_air(i, j), q_air(i, j), slp(i, j), sst_degc(i, j), &
                                  z_wind, z_ta, n_iter, &
                                  taux, tauy, sens, lat, evp, cd, ch, ce, ustar)
            taux_cell(i, j) = taux
            tauy_cell(i, j) = tauy
            q_sens(i, j) = sens
            q_lat(i, j) = lat
            evap_massflux(i, j) = -evp
            cd_out(i, j) = cd
            ch_out(i, j) = ch
            ce_out(i, j) = ce
         else
            taux_cell(i, j) = 0.0_wp
            tauy_cell(i, j) = 0.0_wp
            q_sens(i, j) = 0.0_wp
            q_lat(i, j) = 0.0_wp
            evap_massflux(i, j) = 0.0_wp
            cd_out(i, j) = 0.0_wp
            ch_out(i, j) = 0.0_wp
            ce_out(i, j) = 0.0_wp
         end if
      end do
   end subroutine bulk_flux_driver_2d

end module rdb_ocean_bulk_flux
