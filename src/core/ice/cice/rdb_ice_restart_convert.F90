!! One-shot SIS2 -> Icepack ice-state conversion (roundabout-original; host code).
module rdb_ice_restart_convert
   !! `ice_restart_convert_sis2` turns the SIS2-derived category state
   !! (`rdb_ice_state`: `part_size`, `m_ice`, `m_snow`, `enth_ice`, `sal_ice`,
   !! `enth_snow`) into the Icepack state of `ice_cice_state_t` (`aicen`,
   !! `vicen`, `vsnon`, `Tsfcn`, `qicen`, `sicen`, `qsnon`).  It is plain-array,
   !! HOST-ONLY set-up code that runs ONCE at start-up before the device map
   !! (decision M-D1-c: it lives until one release after the default flips), so
   !! nothing here is a kernel: no `do concurrent`, no directives.
   !!
   !! **What the SIS2 slot stores** (read from `rdb_ice_state`, not assumed).
   !! `part_size(0:ncat)` (category 0 = open water; NOT maintained when
   !! `ncat = 1`, the legacy lumped mode, where a cell is ice-covered with
   !! concentration 1 iff `m_ice(1) > 0`); `m_ice`, `m_snow` [kg/m2] per unit
   !! ICE-covered area (per CELL area in the lumped mode, which is the same
   !! thing at concentration 1); `enth_ice`, `enth_snow` [J/kg] specific
   !! enthalpy, layers BOTTOM-UP (`k = 1` at the ice base), equal-mass layers;
   !! `sal_ice` [ppt] bulk salinity, same order.  No surface temperature is
   !! stored (`tsurf_out` is scratch and not in the restart file).
   !!
   !! The SIS2 enthalpy side is Winton (2000, J. Atmos. Oceanic Technol. 17, 525-531).
   !!
   !! **Lumped mode.**  `ncat = 1` masses are per CELL area.  A film with
   !! `m/rhoi < hi_min` (`p%thermo%hi_min`) becomes a partial cover
   !! `aicen = m/(rhoi*hi_min)` of `hi_min`-thick ice, mass exact; thicker slabs
   !! are concentration 1.
   !!
   !! **Dropped mass.**  Snow on a category with no ice, and categories below
   !! `puny` after merging, are not carried; their mass is in
   !! `mass_*_dropped`, so `in = out + dropped` holds to rounding.
   !!
   !! **Rules (MASTER_PLAN B.5, decision M-D1-g).**
   !!  1. MASS is conserved: `vicen = a*m_ice/rhoi`, `vsnon = a*m_snow/rhos`
   !!     with the Icepack densities of `p%c` (917 and 330), so the ice THICKNESS
   !!     changes by `905/917 - 1 = -1.3 %` while `rhoi*vicen = a*m_ice` holds
   !!     to rounding (within 2 ulp: one division and one multiplication).
   !!  2. Layers are reversed to the Icepack TOP-DOWN order.  When the layer count
   !!     differs the specific enthalpy and the salinity are remapped
   !!     conservatively over normalised depth (overlap-weighted means of
   !!     equal-mass layers, `ice_restart_remap_layers`; the layer mean is
   !!     preserved exactly) BEFORE the temperature is recovered.
   !!  3. Layer TEMPERATURE is preserved, not enthalpy: `T = ice_temp_from_en_s`
   !!     of the (remapped) SIS2 enthalpy and bulk salinity, capped at the new
   !!     layer's melting point, then `qicen = ice_q_from_T(T, sicen)` with the
   !!     Icepack (Bitz-Lipscomb) relation and the Bitz-Lipscomb salinity
   !!     profile (the SIS2 bulk salinity is discarded).  Snow likewise, with
   !!     `ice_qsno_init_from_T`.  The heat content therefore CHANGES (different
   !!     heat capacity, density, liquidus) and the relative change is REPORTED
   !!     (`ice_convert_report_t`), never hidden.  Interim choice (the plan only
   !!     says "T from SIS2 enthalpy, then the Icepack formula").
   !!  4. `Tsfcn` = the snow temperature capped at `Timelt` where the category
   !!     carries snow, else the top ice-layer temperature capped at `Timelt` (the
   !!     surface temperature is not in the SIS2 restart; area-weighted on merge).
   !!     SIS2 snow enthalpy is the Winton (2000) form (`enth_from_TS` at zero
   !!     salinity); partially melted snow (the 0 degC plateau) becomes fully
   !!     FROZEN snow at 0 degC under `ice_qsno_init_from_T`, a heat change that
   !!     is part of the reported figure.
   !!  5. Categories are re-binned on the NEW bounds `hin_max` by their mean
   !!     thickness `h = vicen/aicen` (lower-inclusive bounds; above the top
   !!     bound goes to the top category; a thickness exactly ON a bound goes to
   !!     the UPPER bin, the bounds being lower-inclusive).  Old categories landing in one bin are
   !!     merged: areas and volumes add, layer enthalpies and salinities are
   !!     volume-weighted, `Tsfcn` area-weighted, so volume and energy are
   !!     conserved by the merge.
   !!  6. Total area above 1 (round-off or legacy) is rescaled to 1 keeping the
   !!     volumes (thickness rises), counted in the report.  Categories with
   !!     `aicen <= puny` or `vicen <= puny` after merging are cleaned to the
   !!     empty value, counted as dropped; every empty category (and every land
   !!     or open-water cell) carries EXACT zero `aicen/vicen/vsnon` and the same
   !!     don't-care tracers as `ice_seed_uniform_impl` (it is called with
   !!     `aice = 0` to fill them), so the two paths agree.
   !!
   !! Non-finite or negative input is REFUSED (`ierr`, outputs untouched).  All
   !! cells, ghosts included, are converted cell by cell, so the result does not
   !! depend on a decomposition.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX
   use rdb_ice_params, only: ice_params_t, ice_salinity_profile, ice_q_from_T, ice_qsno_init_from_T
   use rdb_ice_seed, only: ice_seed_uniform_impl, ICE_SEED_OK
   use rdb_ice_enthalpy, only: ice_temp_from_en_s
   implicit none
   private

   public :: ice_restart_convert_sis2, ice_restart_remap_layers, ice_convert_report_t
   public :: ICE_CONV_OK, ICE_CONV_BAD_EXTENT, ICE_CONV_NONFINITE, ICE_CONV_NEGATIVE

   integer, parameter :: ICE_CONV_OK = 0
      !! Converted.
   integer, parameter :: ICE_CONV_BAD_EXTENT = 1
      !! Layer/category counts disagree with `p`, exceed the compile-time bounds,
      !! or `hin_max` is not increasing.  Outputs untouched.
   integer, parameter :: ICE_CONV_NONFINITE = 2
      !! A non-finite input value (or `Tf`).  Outputs untouched.
   integer, parameter :: ICE_CONV_NEGATIVE = 3
      !! A negative `part_size`, `m_ice` or `m_snow`.  Outputs untouched.

   type :: ice_convert_report_t
      !! Totals over ALL cells of the arrays (per unit cell area, summed; cell
      !! areas cancel in the relative figures).
      real(wp) :: mass_ice_in = 0.0_wp      !! `sum a*m_ice` [kg/m2]
      real(wp) :: mass_ice_out = 0.0_wp     !! `sum rhoi*vicen`
      real(wp) :: mass_snow_in = 0.0_wp     !! `sum a*m_snow`
      real(wp) :: mass_snow_out = 0.0_wp    !! `sum rhos*vsnon`
      real(wp) :: mass_ice_dropped = 0.0_wp  !! ice mass cleaned away (categories below `puny`); `in = out + dropped`
      real(wp) :: mass_snow_dropped = 0.0_wp  !! snow mass lost: snow on a category with no ice, or cleaned categories
      real(wp) :: heat_in = 0.0_wp          !! SIS2 enthalpy content `sum a*(m_ice*mean_k(enth) + m_snow*enth_snow)` [J/m2]
      real(wp) :: heat_out = 0.0_wp         !! Icepack content `sum (vicen*mean_k(qicen) + vsnon*qsnon)` [J/m2]
      real(wp) :: heat_rel_change = 0.0_wp  !! `(heat_out - heat_in)/|heat_in|` (0 when `heat_in = 0`)
      integer :: n_rescaled = 0             !! cells whose total area exceeded 1 and was rescaled
      integer :: n_dropped = 0              !! occupied categories cleaned to empty (below `puny`)
      integer :: n_t_capped = 0             !! ice layers whose temperature was capped at the new melting point
   end type ice_convert_report_t

contains

   pure subroutine ice_restart_remap_layers(nki, nko, xin, xout)
      !! Conservative resampling of a layer-mean quantity of `nki` equal layers
      !! onto `nko` equal layers over the same normalised depth: each output is
      !! the overlap-weighted mean of the input layers it covers (integer
      !! overlaps in units of `1/(nki*nko)`), so `mean(xout) = mean(xin)` to
      !! rounding and a uniform input stays uniform.  `nki = nko` is an exact copy.
      integer, intent(in) :: nki, nko
      real(wp), intent(in) :: xin(nki)
      real(wp), intent(out) :: xout(nko)
      integer :: i, j, ov
      real(wp) :: s

      if (nki == nko) then
         xout = xin
         return
      end if
      do j = 1, nko
         s = 0.0_wp
         do i = 1, nki
            ov = min(i*nko, j*nki) - max((i - 1)*nko, (j - 1)*nki)
            if (ov > 0) s = s + real(ov, wp)*xin(i)
         end do
         xout(j) = s/real(nki, wp)
      end do
   end subroutine ice_restart_remap_layers

   pure subroutine ice_restart_convert_sis2(nx, ny, ncat_o, nk_o, ncat, nilyr, nslyr, p, hin_max, Tf, &
                                            part_size, m_ice, m_snow, enth_ice, sal_ice, enth_snow, &
                                            aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, rep, ierr)
      !! Convert the SIS2 state (`ncat_o` categories, `nk_o` layers, bottom-up) to
      !! the Icepack state (`ncat`, `nilyr`, `nslyr = 1`, top-down) with the rules
      !! of the module docstring.  `Tf` [degC] is the ocean freezing point used
      !! for the empty categories (as in the uniform seed).  `ierr` is
      !! `ICE_CONV_*`; non-zero leaves the outputs untouched.
      integer, intent(in) :: nx, ny, ncat_o, nk_o, ncat, nilyr, nslyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: Tf
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat_o)
      real(wp), intent(in) :: m_ice(nx, ny, ncat_o), m_snow(nx, ny, ncat_o)
      real(wp), intent(in) :: enth_ice(nx, ny, ncat_o, nk_o), sal_ice(nx, ny, ncat_o, nk_o)
      real(wp), intent(in) :: enth_snow(nx, ny, ncat_o, 1)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: qsnon(nx, ny, ncat, nslyr)
      type(ice_convert_report_t), intent(out) :: rep
      integer, intent(out) :: ierr

      integer, parameter :: NKO_MAX = 64
      integer :: i, j, c, n, k, ierr_seed
      real(wp) :: a_o, mi, ms, vi, vs, h, asum, scale, T, Tmlt, tsno
      real(wp) :: eb(NKO_MAX), sb(NKO_MAX), et(NKO_MAX), st(NKO_MAX), er(ICE_NILYR_MAX), sr(ICE_NILYR_MAX)
      real(wp) :: tl(ICE_NILYR_MAX), sprof(ICE_NILYR_MAX + 1)
      real(wp) :: acc_a(ncat), acc_v(ncat), acc_s(ncat), acc_ts(ncat)
      real(wp) :: acc_q(ncat, ICE_NILYR_MAX), acc_sal(ncat, ICE_NILYR_MAX), acc_qs(ncat)
      real(wp) :: a_c(ncat_o), v_c(ncat_o), s_c(ncat_o), mic(ncat_o), msc(ncat_o)
      real(wp) :: tl_c(ncat_o, ICE_NILYR_MAX), sl_c(ncat_o, ICE_NILYR_MAX), qs_c(ncat_o), ts_c(ncat_o)
      integer :: bin_c(ncat_o)
      logical :: lumped

      ierr = ICE_CONV_OK
      rep = ice_convert_report_t()
      if (p%ncat /= ncat .or. p%nilyr /= nilyr .or. p%nslyr /= nslyr .or. nslyr /= 1 .or. &
          nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nk_o < 1 .or. nk_o > NKO_MAX .or. ncat_o < 1 .or. ncat < 1) then
         ierr = ICE_CONV_BAD_EXTENT
         return
      end if
      do n = 1, ncat
         if (.not. (hin_max(n) > hin_max(n - 1))) ierr = ICE_CONV_BAD_EXTENT
      end do
      if (ierr /= ICE_CONV_OK) return
      if (.not. (ieee_is_finite(Tf) .and. all(ieee_is_finite(part_size)) .and. all(ieee_is_finite(m_ice)) &
                 .and. all(ieee_is_finite(m_snow)) .and. all(ieee_is_finite(enth_ice)) &
                 .and. all(ieee_is_finite(sal_ice)) .and. all(ieee_is_finite(enth_snow)))) then
         ierr = ICE_CONV_NONFINITE
         return
      end if
      if (any(part_size < 0.0_wp) .or. any(m_ice < 0.0_wp) .or. any(m_snow < 0.0_wp)) then
         ierr = ICE_CONV_NEGATIVE
         return
      end if

      ! Empty value everywhere (also fills the don't-care tracers like the seed).
      call ice_seed_uniform_impl(nx, ny, ncat, nilyr, nslyr, p, hin_max, 0.0_wp, 0.0_wp, 0.0_wp, &
                                 min(Tf, p%c%Timelt), Tf, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                 ierr_seed)
      if (ierr_seed /= ICE_SEED_OK) then
         ierr = ICE_CONV_BAD_EXTENT
         return
      end if
      call ice_salinity_profile(nilyr, p%c%saltmax, p%c%min_salin, sprof(1:nilyr + 1))
      lumped = ncat_o == 1

      do j = 1, ny
         do i = 1, nx
            ! ---- per old category: area, volumes, layer temperatures (top-down) ----
            asum = 0.0_wp
            do c = 1, ncat_o
               a_c(c) = 0.0_wp
               v_c(c) = 0.0_wp
               s_c(c) = 0.0_wp
               mic(c) = 0.0_wp
               msc(c) = 0.0_wp
               bin_c(c) = 0
               mi = m_ice(i, j, c)
               ms = m_snow(i, j, c)
               if (lumped) then
                  ! per-CELL masses: a film thinner than hi_min is a partial cover of hi_min ice
                  a_o = merge(min(1.0_wp, mi/(p%c%rhoi*p%thermo%hi_min)), 0.0_wp, mi > 0.0_wp)
                  mic(c) = mi
                  msc(c) = ms
               else
                  a_o = part_size(i, j, c)
                  mic(c) = a_o*mi
                  msc(c) = a_o*ms
               end if
               rep%mass_ice_in = rep%mass_ice_in + mic(c)
               rep%mass_snow_in = rep%mass_snow_in + msc(c)
               if (.not. (a_o > 0.0_wp .and. mi > 0.0_wp)) then
                  rep%mass_snow_dropped = rep%mass_snow_dropped + msc(c)   ! snow without ice
                  cycle
               end if
               a_c(c) = a_o
               asum = asum + a_o
               rep%heat_in = rep%heat_in + mic(c)*sum(enth_ice(i, j, c, :))/real(nk_o, wp) &
                             + msc(c)*enth_snow(i, j, c, 1)
            end do
            if (asum <= 0.0_wp) cycle
            scale = 1.0_wp
            if (asum > 1.0_wp) then
               scale = 1.0_wp/asum
               rep%n_rescaled = rep%n_rescaled + 1
            end if
            do c = 1, ncat_o
               if (a_c(c) <= 0.0_wp) cycle
               a_o = a_c(c)*scale
               vi = mic(c)/p%c%rhoi      ! volumes keep the ORIGINAL mass (exact)
               vs = msc(c)/p%c%rhos
               ! enthalpy and salinity, bottom-up -> top-down, then resample.
               do k = 1, nk_o
                  eb(k) = enth_ice(i, j, c, nk_o + 1 - k)
                  sb(k) = sal_ice(i, j, c, nk_o + 1 - k)
               end do
               call ice_restart_remap_layers(nk_o, nilyr, eb(1:nk_o), er(1:nilyr))
               call ice_restart_remap_layers(nk_o, nilyr, sb(1:nk_o), sr(1:nilyr))
               do k = 1, nilyr
                  tl_c(c, k) = ice_temp_from_en_s(er(k), sr(k))
               end do
               tsno = ice_temp_from_en_s(enth_snow(i, j, c, 1), 0.0_wp)
               qs_c(c) = ice_qsno_init_from_T(tsno, p%c)
               if (vs > 0.0_wp) then
                  ts_c(c) = min(tsno, p%c%Timelt)      ! under snow the surface is the snow
               else
                  ts_c(c) = min(tl_c(c, 1), p%c%Timelt)
               end if
               a_c(c) = a_o
               v_c(c) = vi
               s_c(c) = vs
               h = vi/a_o
               bin_c(c) = ncat
               do n = 1, ncat
                  if (h < hin_max(n)) then
                     bin_c(c) = n
                     exit
                  end if
               end do
            end do
            ! ---- merge into the new categories ----
            acc_a = 0.0_wp
            acc_v = 0.0_wp
            acc_s = 0.0_wp
            acc_ts = 0.0_wp
            acc_qs = 0.0_wp
            acc_q = 0.0_wp
            acc_sal = 0.0_wp
            do c = 1, ncat_o
               n = bin_c(c)
               if (n == 0) cycle
               acc_a(n) = acc_a(n) + a_c(c)
               acc_v(n) = acc_v(n) + v_c(c)
               acc_s(n) = acc_s(n) + s_c(c)
               acc_ts(n) = acc_ts(n) + a_c(c)*ts_c(c)
               acc_qs(n) = acc_qs(n) + s_c(c)*qs_c(c)
               do k = 1, nilyr
                  ! new-layer temperature capped at the new melting point, then enthalpy
                  if (p%thermo%l_brine) then
                     Tmlt = -p%c%depressT*sprof(k)
                  else
                     Tmlt = 0.0_wp
                  end if
                  T = tl_c(c, k)
                  if (T > Tmlt) then
                     T = Tmlt
                     rep%n_t_capped = rep%n_t_capped + 1
                  end if
                  acc_q(n, k) = acc_q(n, k) + v_c(c)*ice_q_from_T(T, sprof(k), p%c, p%thermo%l_brine)
                  acc_sal(n, k) = acc_sal(n, k) + v_c(c)*sprof(k)
               end do
            end do
            do n = 1, ncat
               if (acc_a(n) <= 0.0_wp) cycle
               if (acc_a(n) <= p%c%puny .or. acc_v(n) <= p%c%puny) then
                  rep%n_dropped = rep%n_dropped + 1
                  rep%mass_ice_dropped = rep%mass_ice_dropped + p%c%rhoi*acc_v(n)
                  rep%mass_snow_dropped = rep%mass_snow_dropped + p%c%rhos*acc_s(n)
                  cycle
               end if
               aicen(i, j, n) = acc_a(n)
               vicen(i, j, n) = acc_v(n)
               vsnon(i, j, n) = acc_s(n)
               Tsfcn(i, j, n) = acc_ts(n)/acc_a(n)
               do k = 1, nilyr
                  qicen(i, j, n, k) = acc_q(n, k)/acc_v(n)
                  sicen(i, j, n, k) = acc_sal(n, k)/acc_v(n)
               end do
               if (acc_s(n) > 0.0_wp) then
                  qsnon(i, j, n, 1) = acc_qs(n)/acc_s(n)
               else
                  qsnon(i, j, n, 1) = ice_qsno_init_from_T(Tsfcn(i, j, n), p%c)
               end if
            end do
         end do
      end do

      ! ---- report ----
      do n = 1, ncat
         do j = 1, ny
            do i = 1, nx
               rep%mass_ice_out = rep%mass_ice_out + p%c%rhoi*vicen(i, j, n)
               rep%mass_snow_out = rep%mass_snow_out + p%c%rhos*vsnon(i, j, n)
               if (aicen(i, j, n) > 0.0_wp) then
                  rep%heat_out = rep%heat_out + vicen(i, j, n)*sum(qicen(i, j, n, :))/real(nilyr, wp) &
                                 + vsnon(i, j, n)*qsnon(i, j, n, 1)
               end if
            end do
         end do
      end do
      if (rep%heat_in /= 0.0_wp) rep%heat_rel_change = (rep%heat_out - rep%heat_in)/abs(rep%heat_in)
   end subroutine ice_restart_convert_sis2

end module rdb_ice_restart_convert
