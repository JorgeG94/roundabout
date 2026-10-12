module rdb_teos10_convert
   !! Pure TEOS-10 conversion primitives (no EOS state, no namelist knobs).
   !!
   !! Salinity-scale and temperature-scale converters between the
   !! "model" variables this tree has used historically (practical
   !! salinity SP, potential temperature PT) and the TEOS-10 variables
   !! (Reference/Absolute Salinity SA/SR, Conservative Temperature CT),
   !! plus the TEOS-10 polynomial freezing point.
   !!
   !! Refs: IOC/SCOR/IAPSO (2010), "The international thermodynamic
   !! equation of seawater - 2010 (TEOS-10)", IOC Manuals and Guides
   !! No. 56 (the "TEOS-10 manual"); McDougall & Barker (2011), "Getting
   !! started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic
   !! Toolbox", SCOR/IAPSO WG127 (the "GSW toolbox" paper, which the
   !! `pt_from_ct` Newton-McDougall iteration and the `ct_from_pt`
   !! surface polynomial are from).
   !!
   !! OM3 parity override (no delta-SA atlas): this tree sets
   !! SA == SR (Reference Salinity from practical salinity alone, no
   !! spatial Absolute-Salinity anomaly atlas) -- the same simplification
   !! the production `EOS_VARIANT_ROQUET_*` branches already make.  Every
   !! routine below that takes an "sa" argument therefore expects SR
   !! (or, under the cons_abs convention a later PR adds, the model's SA
   !! tracer taken as-is); there is no separate `sa_from_sp`.
   !!
   !! All point routines are `pure elemental` and `!$acc routine seq` --
   !! parameters only, no allocation, fixed iteration counts -- so they
   !! are device-callable from inside a `do concurrent` kernel in this
   !! library. The explicit-shape 1-D `*_impl` kernels below wrap them for
   !! callers outside the library (a test executable's own kernel cannot
   !! resolve a library `routine seq` symbol at device link).
   !!
   !! NOT YET WIRED into `rdb_eos`: `roquet_spv_point` / `roquet_rho_point`
   !! carry their OWN copy of a similarly-shaped PT->CT surface
   !! polynomial, independently transcribed, and are NOT bit-identical to
   !! `ct_from_pt` here (see that module's docstrings) -- a later PR either
   !! reconciles the two or documents the gap; this module does not
   !! change `rdb_eos`'s existing behaviour.
   use rdb_constants, only: wp
   implicit none
   private

   public :: sr_from_sp
   public :: ct_from_pt_impl, pt_from_ct_impl, ct_freezing_poly_impl
   public :: sp_from_sr
   public :: ct_from_pt
   public :: pt_from_ct
   public :: ct_freezing_poly
   public :: pt_freezing_poly

   real(wp), parameter, public :: CP0 = 3991.86795711963_wp
      !! TEOS-10 "Conservative Temperature" heat-capacity constant
      !! (TEOS-10 manual eq. 3.3.3; McDougall 2003) -- the reference
      !! value that makes CT a conserved-like quantity: potential
      !! enthalpy = CP0 * CT.  Also the GSW toolbox's `gsw_cp0`.

   real(wp), parameter, public :: SR_FACTOR = 35.16504_wp/35.0_wp
      !! SR = SP * SR_FACTOR, in (g/kg) per unit Practical Salinity
      !! (Reference Salinity from Practical Salinity; TEOS-10 manual
      !! eq. 2.4.1 at the global-average composition, no spatial anomaly
      !! -- see the OM3-parity note above). The single copy:
      !! `rdb_eos`'s `ROQ_SR_FACTOR` is defined from it.

   real(wp), parameter :: CT_SFAC = 0.0248826675584615_wp
      !! Salinity normalisation, in kg/g: `x2 = CT_SFAC*SA` (SA in g/kg)
      !! makes `x2` dimensionless in the `ct_from_pt` and
      !! `gibbs_pt0_pt0` polynomials. Equals 1/(40·35.16504/35)
      !! (McDougall & Barker 2011; GSW `gsw_sfac`).

   ! --- TEOS-10 polynomial freezing-point coefficients (manual §3.33;
   ! McDougall & Barker 2011, the `ct_freezing_poly` fit).  `a`, `b` are
   ! the dissolved-air correction constants; `a = 0.502500117621/SSO`.
   real(wp), parameter :: FZ_SSO = 35.16504_wp
   real(wp), parameter :: FZ_C0 = 0.017947064327968736_wp
   real(wp), parameter :: FZ_C1 = -6.076099099929818_wp
   real(wp), parameter :: FZ_C2 = 4.883198653547851_wp
   real(wp), parameter :: FZ_C3 = -11.88081601230542_wp
   real(wp), parameter :: FZ_C4 = 13.34658511480257_wp
   real(wp), parameter :: FZ_C5 = -8.722761043208607_wp
   real(wp), parameter :: FZ_C6 = 2.082038908808201_wp
   real(wp), parameter :: FZ_C7 = -7.389420998107497_wp
   real(wp), parameter :: FZ_C8 = -2.110913185058476_wp
   real(wp), parameter :: FZ_C9 = 0.2295491578006229_wp
   real(wp), parameter :: FZ_C10 = -0.9891538123307282_wp
   real(wp), parameter :: FZ_C11 = -0.08987150128406496_wp
   real(wp), parameter :: FZ_C12 = 0.3831132432071728_wp
   real(wp), parameter :: FZ_C13 = 1.054318231187074_wp
   real(wp), parameter :: FZ_C14 = 1.065556599652796_wp
   real(wp), parameter :: FZ_C15 = -0.7997496801694032_wp
   real(wp), parameter :: FZ_C16 = 0.3850133554097069_wp
   real(wp), parameter :: FZ_C17 = -2.078616693017569_wp
   real(wp), parameter :: FZ_C18 = 0.8756340772729538_wp
   real(wp), parameter :: FZ_C19 = -2.079022768390933_wp
   real(wp), parameter :: FZ_C20 = 1.596435439942262_wp
   real(wp), parameter :: FZ_C21 = 0.1338002171109174_wp
   real(wp), parameter :: FZ_C22 = 1.242891021876471_wp
   real(wp), parameter :: FZ_A = 0.014289763856964_wp
   real(wp), parameter :: FZ_B = 0.057000649899720_wp

   ! --- pt_from_ct rational initial guess (McDougall & Barker 2011):
   ! pt0 = (PA0 + PA1*s1 + PA2*s1**2 + CT*(PA3 + PA4*s1 + PA5*CT))
   !       / (PB0 + PB1*s1 + CT*(PB2 + PB3*CT)), s1 = SA/SR_FACTOR (SP-like),
   ! CT and pt0 in degC. Coefficients as shipped in the GSW toolbox. ---
   real(wp), parameter :: PA0 = -1.446013646344788e-2_wp
      !! Constant term of the numerator (degC).
   real(wp), parameter :: PA1 = -3.305308995852924e-3_wp
      !! Numerator coefficient of s1, s1 = SA/SR_FACTOR (degC per unit s1).
   real(wp), parameter :: PA2 = 1.062415929128982e-4_wp
      !! Numerator coefficient of s1² (degC per unit s1²).
   real(wp), parameter :: PA3 = 9.477566673794488e-1_wp
      !! Numerator coefficient of CT (dimensionless).
   real(wp), parameter :: PA4 = 2.166591947736613e-3_wp
      !! Numerator coefficient of s1·CT (per unit s1).
   real(wp), parameter :: PA5 = 3.828842955039902e-3_wp
      !! Numerator coefficient of CT² (1/degC).
   real(wp), parameter :: PB0 = 1.0_wp
      !! Constant term of the denominator (dimensionless).
   real(wp), parameter :: PB1 = 6.506097115635800e-4_wp
      !! Denominator coefficient of s1 (per unit s1).
   real(wp), parameter :: PB2 = 3.830289486850898e-3_wp
      !! Denominator coefficient of CT (1/degC).
   real(wp), parameter :: PB3 = 1.247811760368034e-6_wp
      !! Denominator coefficient of CT² (1/degC²).
   real(wp), parameter :: T0_KELVIN = 273.15_wp
      !! Celsius -> Kelvin offset (TEOS-10 manual eq. 2.2).

contains

   pure elemental function sr_from_sp(sp) result(sr)
      !! Reference Salinity from Practical Salinity, OM3-parity
      !! simplification SA == SR (no delta-SA atlas -- see the module
      !! header note). TEOS-10 manual eq. 2.4.1 at global-mean
      !! composition.
      !$acc routine seq
      real(wp), intent(in) :: sp  !! Practical Salinity (PSU).
      real(wp) :: sr  !! Reference Salinity (g/kg).
      sr = sp*SR_FACTOR
   end function sr_from_sp

   pure elemental function sp_from_sr(sr) result(sp)
      !! Inverse of `sr_from_sp`.
      !$acc routine seq
      real(wp), intent(in) :: sr  !! Reference Salinity (g/kg).
      real(wp) :: sp  !! Practical Salinity (PSU).
      sp = sr/SR_FACTOR
   end function sp_from_sr

   pure elemental function ct_from_pt(sa, pt) result(ct)
      !! Conservative Temperature from potential temperature (McDougall
      !! & Barker 2011; GSW `gsw_ct_from_pt`) -- the degree-7-in-PT,
      !! degree-7-in-sqrt(SA) surface polynomial evaluated at p = 0.
      !! `sa` is Reference (or Absolute) Salinity in g/kg, `pt`
      !! potential temperature in degC referenced to p = 0.
      !$acc routine seq
      real(wp), intent(in) :: sa  !! Reference/Absolute Salinity (g/kg).
      real(wp), intent(in) :: pt  !! Potential temperature, p_ref=0 (degC).
      real(wp) :: ct  !! Conservative Temperature (degC).
      real(wp) :: x2, x, y, pot_enthalpy

      x2 = CT_SFAC*max(sa, 0.0_wp)
      x = sqrt(x2)
      y = pt*0.025_wp

      pot_enthalpy = 61.01362420681071_wp &
                     + y*(168776.46138048015_wp &
                          + y*(-2735.2785605119625_wp &
                               + y*(2574.2164453821433_wp &
                                    + y*(-1536.6644434977543_wp &
                                         + y*(545.7340497931629_wp &
                                              + (-50.91091728474331_wp - 18.30489878927802_wp*y)*y))))) &
                     + x2*(268.5520265845071_wp &
                           + y*(-12019.028203559312_wp &
                                + y*(3734.858026725145_wp &
                                     + y*(-2046.7671145057618_wp &
                                          + y*(465.28655623826234_wp &
                                               + (-0.6370820302376359_wp - 10.650848542359153_wp*y)*y)))) &
                           + x*(937.2099110620707_wp &
                                + y*(588.1802812170108_wp &
                                     + y*(248.39476522971285_wp &
                                          + (-3.871557904936333_wp - 2.6268019854268356_wp*y)*y)) &
                                + x*(-1687.914374187449_wp &
                                     + x*(246.9598888781377_wp &
                                          + x*(123.59576582457964_wp - 48.5891069025409_wp*x)) &
                                     + y*(936.3206544460336_wp &
                                          + y*(-942.7827304544439_wp &
                                               + y*(369.4389437509002_wp &
                                                    + (-33.83664947895248_wp - 9.987880382780322_wp*y)*y))))))
      ct = pot_enthalpy/CP0
   end function ct_from_pt

   pure elemental function gibbs_pt0_pt0(sa, pt0) result(g_tt)
      !! Second temperature derivative of the Gibbs function at p = 0,
      !! g_TT(SA, pt0, 0) in J/(kg K²) (McDougall & Barker 2011; GSW
      !! `gsw_gibbs_pt0_pt0`) -- the one extra primitive the modified
      !! Newton-McDougall step in `pt_from_ct` needs beyond `ct_from_pt`.
      !$acc routine seq
      real(wp), intent(in) :: sa   !! Reference/Absolute Salinity (g/kg).
      real(wp), intent(in) :: pt0  !! Potential temperature, p_ref = 0 (degC).
      real(wp) :: g_tt  !! d²g/dT² at p = 0 (J/(kg K²)).
      real(wp) :: x2, x, y, g03, g08

      x2 = CT_SFAC*max(sa, 0.0_wp)
      x = sqrt(x2)
      y = pt0*0.025_wp

      g03 = -24715.571866078_wp &
            + y*(4420.4472249096725_wp &
                 + y*(-1778.231237203896_wp &
                      + y*(1160.5182516851419_wp &
                           + y*(-569.531539542516_wp + y*128.13429152494615_wp))))

      g08 = x2*(1760.062705994408_wp &
                + x*(-86.1329351956084_wp &
                     + x*(-137.1145018408982_wp &
                          + y*(296.20061691375236_wp + y*(-205.67709290374563_wp + 49.9394019139016_wp*y))) &
                     + y*(-60.136422517125_wp + y*10.50720794170734_wp)) &
                + y*(-1351.605895580406_wp &
                     + y*(1097.1125373015109_wp + y*(-433.20648175062206_wp + 63.905091254154904_wp*y))))

      g_tt = (g03 + g08)*0.000625_wp
   end function gibbs_pt0_pt0

   pure elemental function pt_from_ct(sa, ct) result(pt)
      !! Potential temperature (p_ref=0) from Conservative Temperature --
      !! the modified Newton-Raphson ("Newton-McDougall") method, 1.5
      !! iterations fixed (McDougall & Barker 2011; GSW `gsw_pt_from_ct`):
      !! a rational initial guess, one full Newton step with the cheap
      !! analytic `dpt/dct` from the guess itself, then one more step
      !! using the expensive-but-exact `dpt/dct = -cp0/((pt+T0)*gibbs_tt)`
      !! evaluated at the midpoint.  Fixed iteration count, no
      !! data-dependent loop -- device-callable.
      !$acc routine seq
      real(wp), intent(in) :: sa  !! Reference/Absolute Salinity (g/kg).
      real(wp), intent(in) :: ct  !! Conservative Temperature (degC).
      real(wp) :: pt  !! Potential temperature, p_ref=0 (degC).
      real(wp) :: s1, a5ct, b3ct, ct_factor, pt_num, pt_recden
      real(wp) :: pt_old, ptm, dpt_dct, ct_diff

      s1 = sa/SR_FACTOR

      a5ct = PA5*ct
      b3ct = PB3*ct
      ct_factor = PA3 + PA4*s1 + a5ct
      pt_num = PA0 + s1*(PA1 + PA2*s1) + ct*ct_factor
      pt_recden = 1.0_wp/(PB0 + PB1*s1 + ct*(PB2 + b3ct))
      pt = pt_num*pt_recden
      dpt_dct = (ct_factor + a5ct - (PB2 + b3ct + b3ct)*pt)*pt_recden

      ! --- 1st (cheap) Newton-McDougall step ---
      ct_diff = ct_from_pt(sa, pt) - ct
      pt_old = pt
      pt = pt_old - ct_diff*dpt_dct
      ptm = 0.5_wp*(pt + pt_old)

      ! --- 0.5 (exact) Newton-McDougall step at the midpoint ---
      dpt_dct = -CP0/((ptm + T0_KELVIN)*gibbs_pt0_pt0(sa, ptm))
      pt = pt_old - ct_diff*dpt_dct
      ct_diff = ct_from_pt(sa, pt) - ct
      pt_old = pt
      pt = pt_old - ct_diff*dpt_dct
   end function pt_from_ct

   pure elemental function ct_freezing_poly(sa, p, sat_frac) result(ct_f)
      !! Conservative Temperature at the freezing point of seawater
      !! (TEOS-10 manual §3.33; McDougall & Barker 2011; GSW
      !! `gsw_ct_freezing_poly`) -- polynomial fit, not the exact
      !! Newton-Raphson chemical-potential solve.  `p` is SEA PRESSURE
      !! IN DBAR (the TEOS-10/GSW convention this routine and its
      !! coefficients are fit in -- NOT the Pa convention `rdb_eos`'s
      !! `p_ref` uses; a caller crossing that boundary must convert).
      !! `sat_frac` is the saturation fraction of dissolved air in
      !! seawater (0 = no dissolved air, 1 = fully air-saturated).
      !$acc routine seq
      real(wp), intent(in) :: sa  !! Reference/Absolute Salinity (g/kg).
      real(wp), intent(in) :: p  !! Sea pressure (dbar).
      real(wp), intent(in) :: sat_frac  !! Dissolved-air saturation fraction.
      real(wp) :: ct_f  !! CT at the freezing point (degC).
      real(wp) :: sa_r, x, p_r

      sa_r = sa*1.0e-2_wp
      x = sqrt(max(sa_r, 0.0_wp))
      p_r = p*1.0e-4_wp

      ct_f = FZ_C0 + sa_r*(FZ_C1 + x*(FZ_C2 + x*(FZ_C3 + x*(FZ_C4 + x*(FZ_C5 + FZ_C6*x))))) &
             + p_r*(FZ_C7 + p_r*(FZ_C8 + FZ_C9*p_r)) &
             + sa_r*p_r*(FZ_C10 + p_r*(FZ_C12 + p_r*(FZ_C15 + FZ_C21*sa_r)) &
                         + sa_r*(FZ_C13 + FZ_C17*p_r + FZ_C19*sa_r) &
                         + x*(FZ_C11 + p_r*(FZ_C14 + FZ_C18*p_r) + sa_r*(FZ_C16 + FZ_C20*p_r + FZ_C22*sa_r)))

      ! Dissolved-air depression of the freezing point.
      ct_f = ct_f - sat_frac*1.0e-3_wp*(2.4_wp - FZ_A*sa)*(1.0_wp + FZ_B*(1.0_wp - sa/FZ_SSO))
   end function ct_freezing_poly

   pure elemental function pt_freezing_poly(sa, p, sat_frac) result(pt_f)
      !! Potential-temperature freezing point, by composition:
      !! `pt_from_ct(sa, ct_freezing_poly(sa, p, sat_frac))` -- there is
      !! no separate GSW potential-temperature freezing polynomial, the
      !! toolbox gets it from the CT one the same way.  Same pressure
      !! convention as `ct_freezing_poly` (sea pressure in DBAR).
      !$acc routine seq
      real(wp), intent(in) :: sa, p, sat_frac
      real(wp) :: pt_f
      pt_f = pt_from_ct(sa, ct_freezing_poly(sa, p, sat_frac))
   end function pt_freezing_poly

   pure subroutine ct_from_pt_impl(n, sa, pt, ct)
      !! `ct = ct_from_pt(sa, pt)` over `n` points, on the device when the
      !! arrays are mapped.
      integer, intent(in) :: n
      real(wp), intent(in) :: sa(n), pt(n)
      real(wp), intent(out) :: ct(n)
      integer :: i
      do concurrent(i=1:n)
         ct(i) = ct_from_pt(sa(i), pt(i))
      end do
   end subroutine ct_from_pt_impl

   pure subroutine pt_from_ct_impl(n, sa, ct, pt)
      !! `pt = pt_from_ct(sa, ct)` over `n` points.
      integer, intent(in) :: n
      real(wp), intent(in) :: sa(n), ct(n)
      real(wp), intent(out) :: pt(n)
      integer :: i
      do concurrent(i=1:n)
         pt(i) = pt_from_ct(sa(i), ct(i))
      end do
   end subroutine pt_from_ct_impl

   pure subroutine ct_freezing_poly_impl(n, sa, p, sat_frac, ct_f)
      !! `ct_f = ct_freezing_poly(sa, p, sat_frac)` over `n` points.
      integer, intent(in) :: n
      real(wp), intent(in) :: sa(n), p(n), sat_frac(n)
      real(wp), intent(out) :: ct_f(n)
      integer :: i
      do concurrent(i=1:n)
         ct_f(i) = ct_freezing_poly(sa(i), p(i), sat_frac(i))
      end do
   end subroutine ct_freezing_poly_impl

end module rdb_teos10_convert
