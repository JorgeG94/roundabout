!! I1, the Stefan column: one wet column of the CICE-port sea ice (`[ocean.ice] scheme = "cice"`) over a 20 m
!! ocean slab at its freezing point, under a constant -30 C atmosphere, 60 days at dt = 900 s, run through the
!! production engine path (`engine_step` -> `engine_step_ice` -> `engine_step_finalize`).
!!
!! The reference is INDEPENDENT of the model's thermodynamics: `ref_integrate` below is a small one-dimensional
!! enthalpy-method (fixed grid, 2 cm, explicit) solution of the freezing problem of Stefan (1891) for sea ice, written
!! from physical constants and the published constitutive laws (the Pringle et al. 2007 conductivity, the Bitz and
!! Lipscomb 1999 enthalpy of brine-bearing ice at the mean salinity 2.6 ppt of their profile, the freezing point of sea
!! water and its rise with the brine the ice rejects, a surface balance of longwave, sensible and latent heat).  It calls
!! none of the model's thermodynamics.  Its own correctness is pinned by `test_reference_neumann`: with the surface held
!! at a fixed temperature and fresh, constant-property ice it reproduces the exact one-phase solution (Neumann/Stefan:
!! `h = 2 lambda sqrt(kappa t)`, `lambda sqrt(pi) exp(lambda^2) erf(lambda) = c dT / L`) to better than 0.5 percent.
!!
!! SHARED WITH THE MODEL BY DESIGN (the forcing definition, copied from the model's configuration, NOT independent
!! physics, and not checked here): the air density `p/(R (1 + 0.608 q) T)` and the air POTENTIAL temperature referenced
!! to 1000 hPa, `T (1e5/p)^0.286` (0.9 K below `uni_tas`; worth about 1.4 percent of h(60 d)); equal scalar and
!! momentum roughness; the stable-branch profile function (Holtslag and de Bruin; never reached, the surface is warmer
!! than the air throughout, so the layer is unstable and the transfer coefficient about twice the neutral value);
!! the windless sensible-heat term of Jordan et al. (1999), 1 W m-2 K-1; the roughness `z0 = 5e-4 m`, the emissivity
!! 0.985, the 4 ppt salt-flux reference salinity.  A wrong value of any of them in BOTH the model and this file is not
!! seen.
!!
!! WEAKLY SENSITIVE TO THE TURBULENT COEFFICIENTS.  The conductive resistance of the ice dominates the heat transfer
!! (h/k = 0.6 m2 K/W at 1.4 m against about 0.1 for the whole surface), so scaling the model's sensible and latent
!! coefficients by 0.8 moves h(60 d) by only 0.3 percent, below what the gate resolves.  This test does NOT gate the
!! transfer coefficients; the atmosphere kernel's own tests do.
!!
!! THE GATE.  The model's four ice layers (and five-category linear ITD) overshoot the continuum answer by a ONE-SIGNED
!! error: with 2, 4 and 8 layers h(60 d) = 1.3985, 1.3914, 1.3870 m (differences 0.0071, 0.0044, ratio 0.62, Richardson
!! limit 1.3798 m against the reference's 1.3796 m).  So the gate is asymmetric and centred where the 4-layer model
!! sits: `+0.3 percent <= (h - h_ref)/h_ref <= +1.8 percent` on every day 10-60 (measured +0.77 to +1.33; margin 0.47
!! percent to each side).  The ice salinity (0 to 4 ppt) moves the reference by 0.5 percent, inside the margin.  A model
!! error of -1.4 percent in h(60 d) (conductivity x0.95, or Lfresh x1.05) lands at -0.6 percent and fails it.  The heat
!! the surface loses to the atmosphere over the run is gated the same way, `+0.5 .. +2.0 percent` of the reference
!! (measured +1.1).
!!
!! Conservation (independently accumulated): the ocean's heat and salt come from the ocean's own tracer arrays,
!! the ice's from the state totals; the atmosphere's input from the ice's surface-flux terms (`ledger tab`), the
!! hand-over from the ocean assembler's arrays.  The vertical solver's error is a bound
!! (`1.1 ICE_FERRMAX dt sum(A aice_init)`, Bitz and Lipscomb 1999), never a subtracted term.
module test_cice_stefan_column
   use, intrinsic :: iso_fortran_env, only: int32, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_toml_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_errors, only: ICE_ERR_FATAL_MASK
   use rdb_ice_thickness, only: ICE_FERRMAX
   use rdb_ice_ledger_terms, only: ICE_LT_THERM1, ICE_LT_HEAT
   use rdb_ocean_console_stats, only: compute_total_tracer
   implicit none
   private

   public :: collect_cice_stefan_column_tests

   ! ---- the experiment (the example input `cice_stefan_column.toml` carries the same numbers) ----
   real(wp), parameter :: DT = 900.0_wp
   integer, parameter :: NDAYS = 60
   integer, parameter :: NSTEPS = NDAYS*96
      !! 60 days of 900 s steps.
   real(wp), parameter :: TAIR = 243.15_wp
      !! Air temperature, -30 C [K].
   real(wp), parameter :: RLDS = 198.2_wp
      !! Downward longwave [W m-2]: sigma TAIR^4 to 0.01 percent, so that a thin ice surface relaxes to TAIR.
   real(wp), parameter :: HUSS = 1.87e-4_wp
      !! Specific humidity [kg kg-1]: 80 percent of saturation over ice at TAIR.
   real(wp), parameter :: WIND = 5.0_wp
   real(wp), parameter :: PSL = 101325.0_wp
   real(wp), parameter :: H0 = 0.3_wp
      !! Initial ice thickness [m].
   real(wp), parameter :: TICE0 = -10.0_wp
   real(wp), parameter :: SO0 = 34.0_wp
   real(wp), parameter :: HSLAB = 20.0_wp

   ! ---- constants of the reference (physical values; none read from the model) ----
   real(wp), parameter :: SIGSB = 5.670374419e-8_wp
   real(wp), parameter :: RHOI = 917.0_wp
   real(wp), parameter :: LFUS = 3.34e5_wp
   real(wp), parameter :: CPICE = 2106.0_wp
   real(wp), parameter :: MUF = 0.054_wp
      !! Freezing-point slope of sea water [K per ppt] (linear in salinity).
   real(wp), parameter :: RHOW = 1026.0_wp
   real(wp), parameter :: CPW = 3990.0_wp
   real(wp), parameter :: EMISS = 0.985_wp
   real(wp), parameter :: CPAIR = 1005.0_wp
   real(wp), parameter :: LSUB = 2.835e6_wp
   real(wp), parameter :: SALREF = 4.0_wp
   real(wp), parameter :: Z0ICE = 5.0e-4_wp
   real(wp), parameter :: ZREF = 10.0_wp
   real(wp), parameter :: VKARMAN = 0.4_wp

   real(wp), parameter :: SNOWFALL = 1.0e-4_wp
      !! Snowfall of the snow variant [kg m-2 s-1].
   integer, parameter :: SNOW_DAYS = 15
      !! Length of the snow variant [days].
   real(wp), parameter :: SNOW_RHO = 330.0_wp
      !! Density of the fallen snow [kg m-3].
   real(wp), parameter :: GATE_LO = 0.003_wp, GATE_HI = 0.018_wp
      !! Asymmetric gate on `(h - h_ref)/h_ref` (see the header).
   real(wp), parameter :: HEAT_LO = 0.005_wp, HEAT_HI = 0.020_wp

   logical :: comm_inited = .false.

   type :: ref_cfg_t
      !! The reference problem.
      real(wp) :: sice = 2.6_wp
         !! Ice bulk salinity [ppt] (the ice carries it uniformly).
      real(wp) :: z0t_fac = 1.0_wp
         !! Scalar roughness over the momentum roughness.
      real(wp) :: dz = 0.02_wp
      real(wp) :: dt = 100.0_wp
      real(wp) :: so = SO0
         !! Ocean salinity [ppt].
      real(wp) :: hoc = HSLAB
         !! Slab depth [m]; huge freezes the ocean (no salinity feedback).
      real(wp) :: ts_fixed = 1.0e3_wp
         !! A value below 900 holds the surface at that temperature [K] instead of the atmosphere balance.
      real(wp) :: kfix = -1.0_wp
         !! Positive: constant conductivity [W m-1 K-1] (else Pringle et al. 2007).
      real(wp) :: cfix = -1.0_wp
         !! Positive: constant heat capacity [J kg-1 K-1] and no brine (else Bitz and Lipscomb 1999).
      real(wp) :: lam = -1.0_wp
         !! Positive: start from the exact Neumann profile of that root (`kfix`, `cfix`, `ts_fixed` set).
      integer :: ndays = NDAYS
   end type ref_cfg_t

   type :: column_run_t
      !! What one engine run records.
      real(wp) :: h(0:NDAYS) = 0.0_wp
         !! Mean ice thickness at the end of each day [m].
      real(wp) :: hs_end = 0.0_wp, aice_end = 0.0_wp, aice0_end = 0.0_wp, vice_end = 0.0_wp
      real(wp) :: dh_minus = 0.0_wp, dh_plus = 0.0_wp
         !! Largest per-step decrease and increase of the mean thickness.
      real(wp) :: tmargin = huge(1.0_wp)
         !! Lowest `T - Tf(S)` over the ocean cells, layers and steps [K].
      real(wp) :: r_heat = 0.0_wp, bound_heat = 0.0_wp, sc_heat = 0.0_wp, r_salt = 0.0_wp, sc_salt = 0.0_wp
      real(wp) :: r_water = 0.0_wp, sc_water = 0.0_wp, hand_mis = 0.0_wp
      real(wp) :: bl99_ratio = 0.0_wp
         !! Worst per-step `|heat residual in therm1|` over `1.1 ICE_FERRMAX dt sum(A aice_init)` (<= 1).
      real(wp) :: frz_lat = 0.0_wp, atot = 1.0_wp, e_gain = 0.0_wp
         !! New-ice latent heat booked by the ledger [J]; the column area [m2]; heat the ice gained from the atmosphere [J].
      real(wp) :: atm_heat = 0.0_wp
         !! Heat the ice lost to the atmosphere, per unit area [J m-2] (positive = lost).
      real(wp) :: fresh_sum = 0.0_wp
         !! Sum over the run of `fresh` [kg m-2], area-weighted mean.
      real(wp) :: so_end = 0.0_wp, tf_drop = 0.0_wp
      integer :: step_out = 0
         !! First step at which the mean ice thickness is below a nanometre (0: never).
      logical :: finite = .true., fatal = .false., ok = .false.
   end type column_run_t

contains

   subroutine collect_cice_stefan_column_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("reference_neumann", test_reference_neumann), &
                  new_unittest("stefan_growth", test_stefan_growth), &
                  new_unittest("stefan_snow", test_stefan_snow), &
                  new_unittest("melt_out_reversal", test_melt_out) &
                  ]
   end subroutine collect_cice_stefan_column_tests

   ! ------------------------------------------------------------------
   ! the independent reference
   ! ------------------------------------------------------------------

   pure function q_ice(t, s) result(q)
      !! Energy to bring a unit volume of ice of salinity `s` [ppt] at `t` [C] to sea water at its melting point,
      !! [J m-3] (Bitz and Lipscomb 1999, eq. 1): `rho [c (Tm - T) + L (1 - Tm/T)]`, `Tm = -mu s`.
      real(wp), intent(in) :: t, s
      real(wp) :: q, tm

      tm = -MUF*s
      q = RHOI*(CPICE*(tm - t) + LFUS*(1.0_wp - tm/t))
   end function q_ice

   pure function t_of_e(e, s, tf, cfix) result(t)
      !! Temperature [C] of a cell holding enthalpy `e` [J m-3] (zero when melted, `-q_ice` frozen); the cell sits at
      !! the freezing point `tf` while it freezes.  With `cfix > 0` the ice has that constant heat capacity.
      real(wp), intent(in) :: e, s, tf, cfix
      real(wp) :: t, tm, q, b

      if (cfix > 0.0_wp) then
         if (e >= -RHOI*LFUS) then
            t = 0.0_wp
         else
            t = (e + RHOI*LFUS)/(RHOI*cfix)
         end if
         return
      end if
      tm = -MUF*s
      if (e >= -q_ice(tf, s)) then
         t = tf
      else
         q = -e/RHOI
         b = LFUS + CPICE*tm - q
         t = (b - sqrt(b*b - 4.0_wp*CPICE*LFUS*tm))/(2.0_wp*CPICE)
      end if
   end function t_of_e

   pure function qsat_ice(t, p) result(q)
      !! Saturation specific humidity over ice, `t` [K], `p` [Pa] (Murphy and Koop 2005).
      real(wp), intent(in) :: t, p
      real(wp) :: q, es

      es = exp(9.550426_wp - 5723.265_wp/t + 3.53068_wp*log(t) - 0.00728332_wp*t)
      q = 0.622_wp*es/(p - 0.378_wp*es)
   end function qsat_ice

   pure function atm_loss(ts, z0t_fac, qa) result(f)
      !! Net heat the surface at `ts` [K] loses to the atmosphere [W m-2]: longwave (emissivity, downward longwave
      !! `RLDS`) + sensible + latent (sublimation).  The turbulent fluxes are Monin-Obukhov similarity over a rough
      !! surface (`Z0ICE`, scalar roughness `z0t_fac Z0ICE`) at `ZREF`: Paulson (1970) functions for the unstable surface
      !! layer (the ice is warmer than the air at every stage of the growth: the conductive flux needs `Ts > Ta`), iterated
      !! to convergence, the virtual-temperature buoyancy flux, and the windless term of Jordan et al. (1999), `1 W m-2 K-1`.
      !! The air temperature the surface sees is the potential temperature referenced to 1000 hPa, which is how the
      !! model's forcing is defined (`theta = T (1e5/p)^0.286`).
      real(wp), intent(in) :: ts, z0t_fac, qa
      real(wp) :: f, rhoa, theta, qs, ust, tst, qst, zeta, xx, psim, psih, dth, dq, ths, qsurf, vt, lnm, lnh, lnq, tv
      real(wp), parameter :: GRAV = 9.80665_wp
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      integer :: it

      rhoa = PSL/(287.058_wp*(1.0_wp + 0.608_wp*qa)*TAIR)
      theta = TAIR*(1.0e5_wp/PSL)**0.286_wp
      qsurf = qsat_ice(ts, PSL)
      ths = ts
      dth = theta - ths
      dq = qa - qsurf
      lnm = log(ZREF/Z0ICE)
      lnh = log(ZREF/(Z0ICE*z0t_fac))
      lnq = lnh
      ust = VKARMAN*WIND/lnm
      tst = VKARMAN*dth/lnh
      qst = VKARMAN*dq/lnq
      tv = TAIR*(1.0_wp + 0.608_wp*qa)
      do it = 1, 15
         zeta = ZREF*VKARMAN*GRAV*(tst/theta + 0.608_wp*qst/(1.0_wp + 0.608_wp*qa))/ust**2
         zeta = max(min(zeta, 10.0_wp), -10.0_wp)
         if (zeta < 0.0_wp) then
            xx = (1.0_wp - 16.0_wp*zeta)**0.25_wp
            psim = 2.0_wp*log((1.0_wp + xx)/2.0_wp) + log((1.0_wp + xx*xx)/2.0_wp) - 2.0_wp*atan(xx) + 0.5_wp*PI
            psih = 2.0_wp*log((1.0_wp + xx*xx)/2.0_wp)
         else
            ! stable: Holtslag and de Bruin (1988), as used for the stable branch of this surface layer
            psim = -(0.7_wp*zeta + 0.75_wp*(zeta - 14.3_wp)*exp(-0.35_wp*zeta) + 10.7_wp)
            psih = psim
         end if
         ust = VKARMAN*WIND/(lnm - psim)
         tst = VKARMAN*dth/(lnh - psih)
         qst = VKARMAN*dq/(lnq - psih)
      end do
      vt = tv
      f = EMISS*(SIGSB*ts**4 - RLDS) &
          + rhoa*CPAIR*ust*(-tst) + (ts - theta)*1.0_wp &
          + rhoa*LSUB*ust*(-qst)
   end function atm_loss

   subroutine ref_integrate(rc, h_t, e_atm, h_final)
      !! Integrate the reference.  `h_t(d)` is the thickness at the end of day `d`; `e_atm` the heat the surface lost to
      !! the atmosphere (or, for a fixed surface temperature, the conductive heat that reached it) [J m-2].
      type(ref_cfg_t), intent(in) :: rc
      real(wp), intent(out) :: h_t(0:NDAYS), e_atm, h_final
      integer, parameter :: NMAX = 600
      real(wp), parameter :: HINI = H0
      real(wp) :: e(NMAX), temp(NMAX), kk(NMAX), qd(0:NMAX), so, tf, h, hprev, qc, z0f, qa, qthr, dso, qw, lo, hi, mid, glo, ghi, ts_old, tfi
      real(wp) :: z, kappa, ts_s
      integer :: i, n, nstep, istep, it, spd, nd

      n = nint(3.0_wp/rc%dz)
      spd = nint(86400.0_wp/rc%dt)
      nstep = rc%ndays*spd
      so = rc%so
      z0f = rc%z0t_fac
      ts_old = 250.0_wp
      qa = HUSS
      ! the initial ice: H0 thick, linear from TICE0 at the top to the freezing point at the base (as the model seeds it)
      tfi = -MUF*so
      e = 0.0_wp
      nd = nint(HINI/rc%dz)
      do i = 1, nd
         z = (real(i, wp) - 0.5_wp)*rc%dz/HINI
         e(i) = -cell_q(TICE0 + (tfi - TICE0)*z, rc)
      end do
      if (rc%lam > 0.0_wp) then
         ! the exact one-phase solution at the time its front is HINI deep
         kappa = rc%kfix/(RHOI*rc%cfix)
         ts_s = (HINI/(2.0_wp*rc%lam))**2/kappa
         do i = 1, nd
            z = (real(i, wp) - 0.5_wp)*rc%dz
            e(i) = RHOI*rc%cfix*(rc%ts_fixed - 273.15_wp)*(1.0_wp - erf(z/(2.0_wp*sqrt(kappa*ts_s)))/erf(rc%lam)) &
                   - RHOI*LFUS
         end do
      end if
      h_t = 0.0_wp
      h_t(0) = HINI
      hprev = HINI
      e_atm = 0.0_wp
      h = HINI
      do istep = 1, nstep
         tf = -MUF*so
         qthr = cell_q(tf, rc)
         do i = 1, n
            temp(i) = t_of_e(e(i), rc%sice, tf, rc%cfix)
            kk(i) = k_of_t(temp(i), rc)
         end do
         if (rc%ts_fixed < 900.0_wp) then
            qc = kk(1)*(temp(1) - (rc%ts_fixed - 273.15_wp))/(0.5_wp*rc%dz)
         else
            ! secant iteration on the balance `conduction supply(Ts) = atmosphere loss(Ts)`, warm-started from the last step
            lo = ts_old
            hi = ts_old + 0.05_wp
            glo = kk(1)*(temp(1) + 273.15_wp - lo)/(0.5_wp*rc%dz) - atm_loss(lo, z0f, qa)
            do it = 1, 40
               ghi = kk(1)*(temp(1) + 273.15_wp - hi)/(0.5_wp*rc%dz) - atm_loss(hi, z0f, qa)
               if (abs(ghi - glo) < 1.0e-14_wp) exit
               mid = hi - ghi*(hi - lo)/(ghi - glo)
               lo = hi
               glo = ghi
               hi = mid
               if (abs(hi - lo) < 1.0e-7_wp) exit
            end do
            ts_old = hi
            qc = atm_loss(hi, z0f, qa)
         end if
         e_atm = e_atm + qc*rc%dt
         qd(0) = -qc
         do i = 1, n - 1
            qd(i) = -0.5_wp*(kk(i) + kk(i + 1))*(temp(i + 1) - temp(i))/rc%dz
         end do
         qd(n) = 0.0_wp
         do i = 1, n
            e(i) = e(i) + rc%dt*(qd(i - 1) - qd(i))/rc%dz
         end do
         h = 0.0_wp
         do i = 1, n
            if (e(i) <= -qthr) then
               h = h + rc%dz
            else if (e(i) < 0.0_wp) then
               h = h + rc%dz*(-e(i))/qthr
            end if
         end do
         if (mod(istep, spd) == 0) h_t(istep/spd) = h
         ! brine rejection raises the slab salinity, lowering its freezing point: the sensible heat that releases goes to the ice base
         dso = RHOI*(so - SALREF)*(h - hprev)/(RHOW*rc%hoc)
         hprev = h
         so = so + dso
         qw = RHOW*CPW*rc%hoc*MUF*dso
         do i = 1, n
            if (e(i) > -qthr) then
               e(i) = e(i) + qw/rc%dz
               exit
            end if
         end do
      end do
      h_final = h
   end subroutine ref_integrate

   pure function cell_q(t, rc) result(q)
      type(ref_cfg_t), intent(in) :: rc
      real(wp), intent(in) :: t
      real(wp) :: q

      if (rc%cfix > 0.0_wp) then
         q = RHOI*(LFUS - rc%cfix*t)
      else
         q = q_ice(t, rc%sice)
      end if
   end function cell_q

   pure function k_of_t(t, rc) result(k)
      !! Conductivity [W m-1 K-1] of ice at `t` [C] (Pringle et al. 2007, floor 0.1).
      type(ref_cfg_t), intent(in) :: rc
      real(wp), intent(in) :: t
      real(wp) :: k

      if (rc%kfix > 0.0_wp) then
         k = rc%kfix
      else
         k = max(2.11_wp - 0.011_wp*t + 0.09_wp*rc%sice/min(t, -1.0e-3_wp), 0.1_wp)
      end if
   end function k_of_t

   function neumann_lambda(st) result(lam)
      !! Root of `lambda sqrt(pi) exp(lambda^2) erf(lambda) = st` (the one-phase Stefan problem, Stefan 1891).
      real(wp), intent(in) :: st
      real(wp) :: lam, lo, hi, mid, g
      integer :: it
      real(wp), parameter :: PI = 3.14159265358979323846_wp

      lo = 0.0_wp
      hi = 3.0_wp
      do it = 1, 100
         mid = 0.5_wp*(lo + hi)
         g = mid*sqrt(PI)*exp(mid*mid)*erf(mid) - st
         if (g > 0.0_wp) then
            hi = mid
         else
            lo = mid
         end if
      end do
      lam = 0.5_wp*(lo + hi)
   end function neumann_lambda

   subroutine test_reference_neumann(error)
      !! The reference against the exact one-phase Stefan solution: fresh constant-property ice, surface held at
      !! -30 C, water at 0 C (frozen slab).  Started from the exact solution at `t_s`; compared 20 days later.
      type(error_type), allocatable, intent(out) :: error
      type(ref_cfg_t) :: rc
      real(wp) :: h_t(0:NDAYS), ea, hf, k, c, dtemp, st, lam, kappa, ts_s, h_exact, h_stefan0
      real(wp), parameter :: TSURF = -30.0_wp

      k = 2.03_wp
      c = 2106.0_wp
      kappa = k/(RHOI*c)
      dtemp = 0.0_wp - TSURF
      st = c*dtemp/LFUS
      lam = neumann_lambda(st)
      ts_s = (H0/(2.0_wp*lam))**2/kappa
      rc = ref_cfg_t(sice=0.0_wp, dz=0.02_wp, dt=100.0_wp, so=0.0_wp, hoc=1.0e30_wp, &
                     ts_fixed=TSURF + 273.15_wp, kfix=k, cfix=c, lam=lam, ndays=20)
      call ref_integrate(rc, h_t, ea, hf)
      ! exact: h(t) = 2 lam sqrt(kappa (t + t_s)) once the profile has relaxed to the similarity form (a few h0^2/kappa)
      h_exact = 2.0_wp*lam*sqrt(kappa*(20.0_wp*86400.0_wp + ts_s))
      h_stefan0 = sqrt(2.0_wp*k*dtemp*20.0_wp*86400.0_wp/(RHOI*LFUS) + H0**2)
      write (output_unit, "(a,f7.4,a,f7.4,a,f7.4,a,f6.3)") "  Neumann: reference h(20 d) =", hf, ", exact =", h_exact, &
         ", quasi-steady =", h_stefan0, "; Stefan number", st
      call check(error, abs(hf - h_exact)/h_exact < 5.0e-3_wp, "the reference reproduces the exact one-phase Stefan solution")
   end subroutine test_reference_neumann

   ! ------------------------------------------------------------------
   ! the engine run
   ! ------------------------------------------------------------------

   function stefan_nml(tas, rlds_, sw, snow, h_ice, nilyr) result(nml)
      real(wp), intent(in) :: tas, rlds_, sw, snow, h_ice
      integer, intent(in) :: nilyr
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=512) :: a, b, c

      write (a, "(a,f12.3,a,f12.3,a,es12.3,a,f12.3,a)") &
         "[ocean.ice_atm]"//NL//"uni_tas =", tas, NL//"uni_rlds =", rlds_, &
         NL//"uni_huss =", HUSS, NL//"uni_uas =", WIND, NL//"uni_vas = 0.0"//NL//"iceruf = 5.0e-4"//NL
      write (b, "(a,i0,a,f12.2,a,es12.3,a)") &
         "[ocean.ice]"//NL//"enable = true"//NL//'scheme = "cice"'//NL//"ncat = 5"//NL//"nilyr = ", nilyr, &
         NL//"hfreeze = 20.0"//NL//"sw_down =", sw, NL//"snowfall =", snow, &
         NL//"dynamics = true"//NL//"transport = true"//NL//"ledger = true"//NL
      write (c, "(a,f12.3,a,f12.3,a)") &
         "[ocean.ice_ic]"//NL//'conc_config = "uniform"'//NL//"conc = 1.0"//NL//"h_ice =", h_ice, &
         NL//"t_ice =", TICE0, NL//"h_snow = 0.0"//NL
      nml = "[sim]"//NL// &
            'sim_type = "ocean"'//NL// &
            "[time]"//NL// &
            "t_end = 5184000.0"//NL// &
            "dt_fixed = 900.0"//NL// &
            "[nonhydrostatic]"//NL// &
            "nz_layers = 3"//NL// &
            "[tracer]"//NL// &
            "initial_temperature = -1.836"//NL// &
            "initial_salinity = 34.0"//NL// &
            "T_init_surface = -1.836"//NL// &
            "T_init_bottom = -1.836"//NL// &
            "[ocean.bt]"//NL// &
            "auto_n_inner = true"//NL// &
            "[ocean.hvisc]"//NL// &
            "nu_h = 200.0"//NL// &
            "[ocean.diag]"//NL// &
            "enabled = false"//NL// &
            "[grid]"//NL// &
            "nx = 8"//NL// &
            "ny = 8"//NL// &
            "nghost = 3"//NL// &
            "dx = 20000.0"//NL// &
            "dy = 20000.0"//NL// &
            "[physics]"//NL// &
            "coriolis_f = 1.4e-4"//NL// &
            "wind_stress_x = 0.0"//NL// &
            "wind_stress_y = 0.0"//NL// &
            "[vcoord]"//NL// &
            'vcoord_type = "sigma"'//NL// &
            "[ocean.topo]"//NL// &
            'topo_config = "flat"'//NL// &
            "max_depth = 20.0"//NL// &
            "[ocean.thermo]"//NL// &
            "enable_thermodynamics = true"//NL// &
            "[ocean.forcing]"//NL// &
            "enable_components = true"//NL// &
            "[ocean.vmix]"//NL// &
            "dt_therm_ratio = 1"//NL// &
            trim(b)//NL// &
            "[ocean.ice_thermo]"//NL// &
            "update_ocn_f = true"//NL// &
            "emissivity = 0.985"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'advection = "upwind"'//NL// &
            "ndte = 40"//NL// &
            trim(c)//NL// &
            trim(a)//NL// &
            "[ocean.bc]"//NL// &
            'west = "periodic"'//NL// &
            'east = "periodic"'//NL// &
            'south = "periodic"'//NL// &
            'north = "periodic"'//NL// &
            "[output]"//NL// &
            "output_to_file = false"//NL
   end function stefan_nml

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   subroutine pull_state(engine)
      type(ocean_engine_t), intent(inout) :: engine

      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update self(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%vsno)
      !$acc update self(engine%state%ice%cice%aice0, engine%state%ice%cice%aice_init, engine%state%ice%cice%ice_err)
      !$acc update self(engine%state%icec%flux%fresh, engine%state%icec%flux%fsalt)
      !$acc update self(engine%state%icec%flux%fhocn, engine%state%icec%flux%fswthru)
      !$acc update self(engine%state%icec%oin%sss_ocn)
      !$acc update self(engine%state%surface_flux%heat_ice, engine%state%surface_flux%sw_ice)
      !$acc update self(engine%state%surface_flux%q_sw, engine%state%surface_flux%q_lw, engine%state%surface_flux%q_lat)
      !$acc update self(engine%state%surface_flux%q_sens, engine%state%surface_flux%heat_content_lprec)
      !$acc update self(engine%state%surface_flux%q_fus)
      !$acc update self(engine%state%surface_flux%heat_content_fprec, engine%state%surface_flux%evap)
      !$acc update self(engine%state%ice%frazil_heat, engine%state%ice%frzmlt, engine%state%ice%heat_budget_frazil)
      !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c)
      !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_salinity)%c)
   end subroutine pull_state

   subroutine run_column(nml, nsteps, rd)
      !! Run the engine `nsteps` steps; record the thickness curve, the ocean's margin above freezing and the closures.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: nsteps
      type(column_run_t), intent(out) :: rd
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      integer :: ierr, n, i, j, k, ng, is, ie, js, je, nz, spd
      real(wp) :: t, hmean, hprev, atot, rho0, h0o, so0o, h1o, so1o, m0, si0, e0, m1, si1, e1, a_dt
      real(wp) :: icefl_all, ice_app, ice_pend, ice_asm, ice_asm_last, open_app, open_pend, open_forc, icefl
      real(wp) :: sum_ext, w_in, fresh_all, dil_app, dil_pend, dil_n, fsalt_last, fresh_n, sst, tfz, rain_in
      real(wp) :: fsalt_n, atm_sum, expected
      real(wp), allocatable :: area(:, :)

      call ensure_comm()
      call read_config_from_toml_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      nz = size(engine%state%multilayer%h_layer, 3)
      spd = 96
      area = engine%state%metrics%areaT
      atot = sum(area(is:ie, js:je))
      rd%atot = atot
      rho0 = engine%state%eos%rho0
      associate (ms => engine%state%multilayer)
         h0o = engine%state%eos%rho0*engine%state%eos%c_p*compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, &
                                                                               engine%state%metrics%areaT, ng)
         so0o = rho0*compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, engine%state%metrics%areaT, ng)
      end associate
      call engine%state%icec%ledger%totals(engine%state%ice%cice, m0, si0, e0)
      t = 0.0_wp
      hprev = H0
      rd%h(0) = H0
      icefl_all = 0.0_wp
      ice_app = 0.0_wp
      ice_pend = 0.0_wp
      ice_asm_last = 0.0_wp
      open_app = 0.0_wp
      open_pend = 0.0_wp
      dil_app = 0.0_wp
      dil_pend = 0.0_wp
      fsalt_last = 0.0_wp
      fresh_all = 0.0_wp
      sum_ext = 0.0_wp
      w_in = 0.0_wp
      a_dt = 0.0_wp
      atm_sum = 0.0_wp
      rain_in = 0.0_wp
      rd%ok = .true.
      do n = 1, nsteps
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) then
            rd%ok = .false.
            exit
         end if
         t = t + DT
         call pull_state(engine)
         associate (c => engine%state%ice%cice, f => engine%state%icec%flux, sf => engine%state%surface_flux, &
                    ms => engine%state%multilayer, L => engine%state%icec%ledger)
            if (any(iand(c%ice_err, ICE_ERR_FATAL_MASK) /= 0_int32)) rd%fatal = .true.
            if (.not. L%valid .or. .not. L%finite) rd%finite = .false.
            hmean = 0.0_wp
            icefl = 0.0_wp
            ice_asm = 0.0_wp
            open_forc = 0.0_wp
            dil_n = 0.0_wp
            fresh_n = 0.0_wp
            fsalt_n = 0.0_wp
            do j = js, je
               do i = is, ie
                  hmean = hmean + area(i, j)*c%vice(i, j)
                  icefl = icefl + DT*area(i, j)*(f%fhocn(i, j) + f%fswthru(i, j))
                  ice_asm = ice_asm + DT*area(i, j)*(sf%heat_ice(i, j) + sf%sw_ice(i, j))
                  sst = ms%tracers(ms%idx_temperature)%c(i, j, nz)
                  open_forc = open_forc + DT*area(i, j)*(sf%q_sw(i, j) + sf%q_lw(i, j) + sf%q_lat(i, j) + sf%q_sens(i, j) &
                                                         + sf%q_fus(i, j) &
                                                         + sf%heat_content_lprec(i, j) + sf%heat_content_fprec(i, j) &
                                                         + engine%state%eos%c_p*sst*sf%evap(i, j))
                  dil_n = dil_n + DT*area(i, j)*engine%state%icec%oin%sss_ocn(i, j)*f%fresh(i, j)
                  fresh_n = fresh_n + DT*area(i, j)*f%fresh(i, j)
                  fsalt_n = fsalt_n + DT*area(i, j)*f%fsalt(i, j)
                  a_dt = a_dt + DT*area(i, j)*c%aice_init(i, j)
                  do k = 1, nz
                     ! the slab never sits below the freezing point of its own salinity (linear in S)
                     tfz = -MUF*ms%tracers(ms%idx_salinity)%c(i, j, k)
                     rd%tmargin = min(rd%tmargin, ms%tracers(ms%idx_temperature)%c(i, j, k) - tfz)
                  end do
               end do
            end do
            hmean = hmean/atot
            rd%dh_minus = max(rd%dh_minus, hprev - hmean)
            rd%dh_plus = max(rd%dh_plus, hmean - hprev)
            hprev = hmean
            if (rd%step_out == 0 .and. hmean < 1.0e-9_wp) rd%step_out = n
            if (mod(n, spd) == 0 .and. n/spd <= NDAYS) rd%h(n/spd) = hmean
            ! the ocean applies the heat assembled by the PREVIOUS step
            open_app = open_app + open_pend
            ice_app = ice_app + ice_pend
            dil_app = dil_app + dil_pend
            open_pend = open_forc
            ice_pend = ice_asm
            dil_pend = dil_n
            fsalt_last = fsalt_n
            ice_asm_last = ice_asm
            icefl_all = icefl_all + icefl
            fresh_all = fresh_all + fresh_n
            sum_ext = sum_ext + L%tab%lw_in - L%tab%snow_enth + L%tab%atm_rad + L%tab%atm_lat
            atm_sum = atm_sum + L%tab%atm_rad + L%tab%atm_lat + L%tab%lw_in
            w_in = w_in + L%tab%snow + L%tab%evap
            rain_in = rain_in + L%tab%rain
            rd%frz_lat = rd%frz_lat + L%tab%frz_lat
            rd%bl99_ratio = max(rd%bl99_ratio, abs(L%tab%res(ICE_LT_HEAT, ICE_LT_THERM1))/ &
                                (1.1_wp*ICE_FERRMAX*max(sum(DT*area(is:ie, js:je)*c%aice_init(is:ie, js:je)), tiny(1.0_wp))))
         end associate
         if (.not. rd%finite .or. rd%fatal) exit
      end do
      if (rd%ok) then
         associate (ms => engine%state%multilayer, c => engine%state%ice%cice)
            h1o = engine%state%eos%rho0*engine%state%eos%c_p*compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, &
                                                                                  engine%state%metrics%areaT, ng)
            so1o = rho0*compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, engine%state%metrics%areaT, ng)
            call engine%state%icec%ledger%totals(c, m1, si1, e1)
            rd%hs_end = sum(area(is:ie, js:je)*c%vsno(is:ie, js:je))/atot
            rd%aice_end = sum(area(is:ie, js:je)*c%aice(is:ie, js:je))/atot
            rd%aice0_end = minval(c%aice0(is:ie, js:je))
            rd%vice_end = hprev
         end associate
         ! heat: dH_ocean + dE_ice = atmosphere -> ice + the open water's share + the hand-over applied - the hand-over booked
         expected = sum_ext + open_app + ice_app - icefl_all
         rd%r_heat = (h1o - h0o) + (e1 - e0) - expected
         rd%bound_heat = 1.1_wp*ICE_FERRMAX*a_dt
         rd%sc_heat = max(abs(h1o - h0o), abs(e1 - e0), abs(sum_ext), abs(icefl_all), 1.0_wp)
         rd%hand_mis = abs(ice_app + ice_asm_last - icefl_all)/max(abs(icefl_all), 1.0_wp)
         ! salt [kg]: ice + (ocean salt content)/1000 + the salt flux in transit + the applied virtual-salt dilution
         rd%r_salt = (si1 - si0) + (so1o - so0o)/1000.0_wp + fsalt_last + dil_app/1000.0_wp
         rd%sc_salt = max(abs(si1 - si0), abs(so1o - so0o)/1000.0_wp, abs(dil_app)/1000.0_wp, 1.0e-30_wp)
         ! water [kg]: the ice's mass change + what it booked to the ocean (`fresh`) = what the atmosphere delivered to it
         rd%r_water = (m1 - m0) + fresh_all - w_in - rain_in
         rd%sc_water = max(abs(m1 - m0), abs(fresh_all), abs(w_in), 1.0e-30_wp)
         rd%atm_heat = -atm_sum/atot
         rd%fresh_sum = fresh_all/atot
         rd%so_end = so1o/(rho0*HSLAB*atot)
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_column

   subroutine print_curve(label, rd, rh, dmax)
      character(len=*), intent(in) :: label
      type(column_run_t), intent(in) :: rd
      real(wp), intent(in), optional :: rh(0:NDAYS)
      integer, intent(in), optional :: dmax
      integer :: d, d1

      write (output_unit, "(a,a)") "  ", label
      d1 = NDAYS
      if (present(dmax)) d1 = dmax
      do d = 5, d1, 5
         if (present(rh)) then
            write (output_unit, "(4x,a,i3,a,f8.4,a,f8.4,a,f8.3,a)") "day", d, ": model", rd%h(d), ", reference", rh(d), &
               ", difference", 100.0_wp*(rd%h(d) - rh(d))/rh(d), " %"
         else
            write (output_unit, "(4x,a,i3,a,f8.4)") "day", d, ": model", rd%h(d)
         end if
      end do
   end subroutine print_curve

   subroutine check_closures(error, rd)
      type(error_type), allocatable, intent(out) :: error
      type(column_run_t), intent(in) :: rd

      call check(error, rd%ok .and. rd%finite .and. .not. rd%fatal, "ran to the end: finite ledger, no fatal ice error bit")
      if (allocated(error)) return
      write (output_unit, "(a,es10.3,a,es10.3,a,es10.3)") "    heat residual", rd%r_heat, " J, solver bound", rd%bound_heat, &
         " J, scale", rd%sc_heat
      write (output_unit, "(a,es10.3,a,es10.3,a,es10.3)") "    salt residual (relative)", abs(rd%r_salt)/rd%sc_salt, &
         "; water residual (relative)", abs(rd%r_water)/rd%sc_water, "; hand-over mismatch", rd%hand_mis
     write (output_unit, "(a,es10.3)") "    worst per-step heat residual in therm1 over the solver tolerance (<= 1):", rd%bl99_ratio
      call check(error, rd%bl99_ratio <= 1.0_wp, "every step's heat residual in therm1 is within the solver tolerance")
      if (allocated(error)) return
      call check(error, abs(rd%r_heat) <= rd%bound_heat + 1.0e-11_wp*rd%sc_heat, &
                 "ocean + ice heat closes within the vertical solver's tolerance bound")
      if (allocated(error)) return
      call check(error, abs(rd%r_salt)/rd%sc_salt < 1.0e-10_wp, "ocean + ice salt (with the virtual-salt dilution) closes")
      if (allocated(error)) return
      call check(error, abs(rd%r_water)/rd%sc_water < 1.0e-12_wp, "ice mass change = atmosphere exchange - booked fresh water")
      if (allocated(error)) return
      call check(error, rd%hand_mis < 1.0e-12_wp, "the ice hand-over: applied + in transit == what the ice booked")
   end subroutine check_closures

   subroutine test_stefan_growth(error)
      !! I1 gates 1-3: thickness against the independent reference (asymmetric gate, see the module header), monotone
      !! growth, closures, slab T >= Tf, no new-ice (frazil) heat, and the heat the surface lost to the atmosphere.
      type(error_type), allocatable, intent(out) :: error
      type(column_run_t) :: rd
      type(ref_cfg_t) :: rc
      real(wp) :: hc(0:NDAYS), ea, hf, dev, devmin, devmax, dstefan, edev
      integer :: d

      call run_column(stefan_nml(TAIR, RLDS, 0.0_wp, 0.0_wp, H0, 4), NSTEPS, rd)
      call check(error, rd%ok, "the engine ran 60 days")
      if (allocated(error)) return
      rc = ref_cfg_t()
      call ref_integrate(rc, hc, ea, hf)
      call print_curve("ice thickness [m]: model vs reference:", rd, hc)
      devmin = huge(1.0_wp)
      devmax = -huge(1.0_wp)
      do d = 10, NDAYS
         dev = (rd%h(d) - hc(d))/hc(d)
         devmin = min(devmin, dev)
         devmax = max(devmax, dev)
      end do
      write (output_unit, "(a,f8.4,a,f8.4)") "    day 60: model", rd%h(NDAYS), ", reference", hc(NDAYS)
      write (output_unit, "(a,f8.4,a,f8.4,a,f6.3,a,f6.3)") "    days 10-60, relative deviation from", devmin, " to", devmax, &
         "; gate", GATE_LO, " to", GATE_HI
      write (output_unit, "(a,f8.4,a,f8.4,a,f8.4)") "    margins: lower", devmin - GATE_LO, ", upper", GATE_HI - devmax
      ! the zero-surface-resistance Stefan law of the plan (h^2 = h0^2 + 2 k dT t/(rho L), k = 2.35 above any ice value)
      dstefan = sqrt(H0**2 + 2.0_wp*2.35_wp*(30.0_wp - 1.836_wp)*real(NDAYS, wp)*86400.0_wp/(RHOI*LFUS))
      call check(error, rd%h(NDAYS) < dstefan, "h(60 d) is below the zero-surface-resistance Stefan law (k = 2.35)")
      if (allocated(error)) return
      call check(error, devmin >= GATE_LO, "h(t) not below the reference + 0.3 percent (days 10-60)")
      if (allocated(error)) return
      call check(error, devmax <= GATE_HI, "h(t) not above the reference + 1.8 percent (days 10-60)")
      if (allocated(error)) return
      call check(error, rd%dh_minus < 1.0e-6_wp, "growth is monotone (no step lost more than 1 micron)")
      if (allocated(error)) return
      call check(error, rd%tmargin >= -1.0e-12_wp, "slab temperature stays at or above the freezing point")
      if (allocated(error)) return
      write (output_unit, "(a,es10.3,a,es10.3)") "    new-ice latent heat booked [J]:", rd%frz_lat, ", atmosphere heat loss x area", rd%atm_heat*rd%atot
      ! A trace of frazil ice forms (brine rejection and the one-step-old freezing point); its latent heat is paid by the
      ! ocean's own tracers, which the closure already holds, so it is not an extra source.  It must stay far below the
      ! throughput, and the closure must hold within the solver bound with it added in either sign.
   call check(error, abs(rd%frz_lat) <= 2.0e-6_wp*abs(rd%atm_heat)*rd%atot, "new-ice latent heat below 2e-6 of the heat throughput")
      if (allocated(error)) return
      call check(error, abs(abs(rd%r_heat) + abs(rd%frz_lat)) <= rd%bound_heat, "closure holds within the solver bound even with the new-ice heat added")
      if (allocated(error)) return
      call check_closures(error, rd)
      if (allocated(error)) return
      edev = (rd%atm_heat - ea)/ea
 write (output_unit, "(a,es10.3,a,es10.3,a,f8.4,a,f6.3,a,f6.3)") "    atmosphere heat loss over 60 d [J m-2]: model", rd%atm_heat, &
         ", reference", ea, ", deviation", edev, "; gate", HEAT_LO, " to", HEAT_HI
call check(error, edev >= HEAT_LO .and. edev <= HEAT_HI, "heat lost to the atmosphere within +0.5 .. +2.0 percent of the reference")
   end subroutine test_stefan_growth

   subroutine test_stefan_snow(error)
      !! Variant with snowfall (1e-4 kg m-2 s-1, 130 kg m-2 over 15 days): the snow insulates (the ice is more than 10
      !! percent thinner than the snow-free reference's at day 15); the snow that floods through the negative freeboard
      !! becomes ice (less than 80 percent of the fallen snow mass is still snow; sublimation is two orders smaller and
      !! nothing melts), and the mass, salt and heat closures hold through that conversion.  No independent reference for
      !! the snow itself.
      type(error_type), allocatable, intent(out) :: error
      type(column_run_t) :: rs
      type(ref_cfg_t) :: rc
      real(wp) :: hc(0:NDAYS), ea, hf

      call run_column(stefan_nml(TAIR, RLDS, 0.0_wp, SNOWFALL, H0, 4), SNOW_DAYS*96, rs)
      call check(error, rs%ok, "the run completed")
      if (allocated(error)) return
      rc = ref_cfg_t(ndays=SNOW_DAYS)
      call ref_integrate(rc, hc, ea, hf)
      call print_curve("snow variant, ice thickness [m]:", rs, hc, SNOW_DAYS)
      write (output_unit, "(a,f8.4,a,f8.4,a,f8.4,a,f8.4)") "    h(15 d) with snow", rs%h(SNOW_DAYS), ", snow-free reference", &
 hc(SNOW_DAYS), ", snow depth", rs%hs_end, ", snow mass fraction left", SNOW_RHO*rs%hs_end/(SNOWFALL*real(SNOW_DAYS, wp)*86400.0_wp)
      call check(error, rs%hs_end > 0.05_wp, "snow accumulated (more than 5 cm)")
      if (allocated(error)) return
      call check(error, SNOW_RHO*rs%hs_end < 0.8_wp*SNOWFALL*real(SNOW_DAYS, wp)*86400.0_wp, &
                 "snow-ice: a fifth or more of the fallen snow was converted to ice through the freeboard")
      if (allocated(error)) return
      call check(error, rs%h(SNOW_DAYS) < 0.9_wp*hc(SNOW_DAYS), "snow insulates: the ice is more than 10 percent thinner than the snow-free reference")
      if (allocated(error)) return
      call check(error, rs%dh_minus < 1.0e-6_wp, "growth is monotone")
      if (allocated(error)) return
      call check_closures(error, rs)
   end subroutine test_stefan_snow

   subroutine test_melt_out(error)
      !! Reversal: warm air and shortwave melt a 0.1 m ice out (day 8.1 measured); meltwater goes to the ocean (the slab
      !! freshens), the cell returns to open water (`aice0 = 1`, to a nanometre of residue: a 1e-12 threshold was fragile, a sub-nanometre sliver
      !! outlives melt-out for a while when Lfresh is 5 percent larger) and stays there, and the closures hold.  The melt-out day
      !! (6.9-9.4) is a regression pin on the measured 8.1, not independent physics; the independent part is the energy bound.
      type(error_type), allocatable, intent(out) :: error
      type(column_run_t) :: rd
      real(wp) :: eneed
      integer :: i

      call run_column(stefan_nml(278.15_wp, 330.0_wp, 250.0_wp, 0.0_wp, 0.1_wp, 4), 10*96, rd)
      call check(error, rd%ok, "the engine ran")
      if (allocated(error)) return
      write (output_unit, "(a,f7.3,a,es10.3,a,es10.3,a,f8.3)") "    melt-out at day", real(rd%step_out, wp)/96.0_wp, ", aice0", rd%aice0_end, &
         ", fresh [kg m-2]", rd%fresh_sum, ", slab salinity", rd%so_end
      call check(error, rd%vice_end < 1.0e-9_wp .and. rd%aice0_end > 1.0_wp - 1.0e-9_wp, "melted out: open water, aice0 = 1 (no residue above a nanometre)")
      if (allocated(error)) return
      call check(error, rd%step_out > 0 .and. rd%step_out < 10*96, "the ice melted out before the end, and stayed out")
      if (allocated(error)) return
      eneed = 0.0_wp
      do i = 1, 40
         eneed = eneed + 0.1_wp*q_ice(TICE0 + (-MUF*SO0 - TICE0)*(real(i, wp) - 0.5_wp)/40.0_wp, 2.6_wp)/40.0_wp
      end do
      write (output_unit, "(a,es10.3,a,es10.3,a,f7.3)") "    energy to melt the seed [J m-2]:", eneed, "; heat the ice gained from the atmosphere:", &
         -rd%atm_heat, "; ratio", -rd%atm_heat/eneed
      ! the ice cannot melt without receiving its enthalpy (ratio >= 1); the excess is the shortwave that passes through the thin
      ! ice into the ocean, which the same ledger term counts (at most a third of 250 W m-2 over the 8 days)
      call check(error, -rd%atm_heat/eneed >= 1.0_wp .and. -rd%atm_heat/eneed < 1.7_wp, &
                 "the ice gained at least the energy its enthalpy needs to melt, and not 70 percent more")
      if (allocated(error)) return
      call check(error, real(rd%step_out, wp)/96.0_wp > 6.9_wp .and. real(rd%step_out, wp)/96.0_wp < 9.4_wp, &
                 "melt-out day within 15 percent of the 8.1 measured")
      if (allocated(error)) return
      call check(error, rd%fresh_sum > 0.0_wp, "meltwater: fresh has the sign opposite to growth")
      if (allocated(error)) return
      call check(error, rd%so_end < SO0, "the ocean received the meltwater: the slab is fresher")
      if (allocated(error)) return
      call check_closures(error, rd)
   end subroutine test_melt_out

end module test_cice_stefan_column
