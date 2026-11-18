!! Shortwave drivers of the CICE-port sea ice: band split, per-category step, cell aggregate (roundabout-original).
module rdb_ice_radiation
   !! The grid-sweeping drivers around the CCSM3 column kernels of `rdb_ice_shortwave`
   !! (Briegleb et al. 2004, NCAR/TN-463+STR).  Three pieces:
   !!
   !! 1. **Band split** (`ice_sw_band_split`).  A forced run has one broadband downward
   !!    shortwave `sw_down`; the visible share is the knob `sw_vis_frac` (default 0.52):
   !!    `swv = sw*sw_vis_frac`, `swn = sw - swv`.  The direct and diffuse albedos of this
   !!    scheme are identical, so only the visible and near-IR TOTALS matter and no
   !!    direct/diffuse split exists here (a consumer that needs one, e.g. a delta-Eddington
   !!    scheme, would add it).  `swn` is formed by subtraction so the two bands add back to
   !!    `sw` to round-off by construction.
   !! 2. **Category step** (`ice_step_radiation_impl`, ONE `do concurrent (n, j, i)` launch):
   !!    albedos and absorbed shortwave per (i, j, n), written ONLY to per-category outputs
   !!    `(nx,ny,ncat[,k])`; nothing is accumulated across `n`.  A category with
   !!    `aicen <= puny` gets the defaults: both albedos `ICE_ALBOCN`, every absorbed
   !!    flux and layer absorption zero.  A category with a non-finite `aicen, vicen,
   !!    vsnon, Tsfcn` or `sw_down` gets the same defaults and `ICE_ERR_NONFINITE` in its own
   !!    `errn(i,j,n)` word (an ice-free category with a NaN is also flagged: the state is
   !!    corrupt).  The tests come FIRST and each output is then written exactly once.
   !!    No other cell or category is touched.
   !! 3. **Cell aggregate** (`ice_radiation_aggregate_impl`, ONE `do concurrent (j, i)` launch):
   !!    `alv_ai = sum_n aicen*alvn` over categories with `aicen > puny`, in the fixed order
   !!    `n = 1..ncat`, and likewise `ali_ai`.  This is the PARTIAL sum over the ice
   !!    categories: NOT divided by the ice area and WITHOUT open water.  A grid-cell albedo
   !!    is `alv_ai + aice0*(ocean albedo)`; an ice-only mean is `alv_ai/aice`.  The
   !!    per-category error words are ORed into `ice_err(i,j)`.
   !!
   !! **Call order (`sw_lag = .false.`, the only mode).**  Once per ice step, AFTER the
   !! forcing import and `ice_save_init`, BEFORE therm1 and before anything that changes
   !! the category state (new ice, ITD, ridging, transport).  `aicen` is then equal to
   !! `aicen_init`, so the albedo weights are the same start-of-step areas that therm1's
   !! flux merge uses; the call is also correct before `ice_save_init`.  The state is
   !! recomputed from scratch every step: nothing is read from the previous step.
   !!
   !! **Fail loud.**  The host shim returns `ierr /= ICE_RAD_OK` (and runs nothing) for an
   !! uninitialised or unconfigured slot, or extents that differ from the slot's; the caller
   !! MUST treat that as fatal, never as "no shortwave".
   !!
   !! **Not here.**  The `fswthru` / `fswthru_vis` flux producer (`sum_n aicen*fswthrun`
   !! into the flux slot) and feeding therm1 are the wiring PR's.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_sw_params_t
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_shortwave, only: ice_albedos_ccsm3, ice_absorbed_ccsm3, ICE_ALBOCN, ICE_SW_PUNY
   implicit none
   private

   public :: ice_sw_band_split, ice_radiation_cell
   public :: ice_step_radiation_impl, ice_radiation_aggregate_impl, ice_step_radiation
   public :: ICE_RAD_OK, ICE_RAD_ERR_NOT_INIT, ICE_RAD_ERR_NOT_CONFIGURED, ICE_RAD_ERR_SHAPE

   integer, parameter :: ICE_RAD_OK = 0
      !! `ice_step_radiation` ran.
   integer, parameter :: ICE_RAD_ERR_NOT_INIT = 1
      !! The slot is not initialised; nothing ran.
   integer, parameter :: ICE_RAD_ERR_NOT_CONFIGURED = 2
      !! `rad%configure` was not called; nothing ran.
   integer, parameter :: ICE_RAD_ERR_SHAPE = 3
      !! Extents (or the cell range) disagree with the slot; nothing ran.

contains

   pure subroutine ice_sw_band_split(sw, vis_frac, swv, swn)
      !! Visible and near-IR totals of the broadband downward shortwave `sw` [W m-2]:
      !! `swv = sw*vis_frac`, `swn = sw - swv`.
      !$acc routine seq
      real(wp), intent(in) :: sw
         !! Broadband downward shortwave (W/m2).
      real(wp), intent(in) :: vis_frac
         !! Visible share of `sw` (`sw_vis_frac`, in [0, 1]).
      real(wp), intent(out) :: swv, swn
         !! Visible / near-IR shortwave (W/m2).
      swv = sw*vis_frac
      swn = sw - swv
   end subroutine ice_sw_band_split

   pure subroutine ice_radiation_cell(i, j, n, nx, ny, ncat, nilyr, &
                                      albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, vis_frac, &
                                      aicen, vicen, vsnon, tsfcn, sw_down, &
                                      alvn, alin, fswsfcn, fswintn, fswthrun, iswabsn, errn)
      !! One (i, j, n) work item: albedos and absorbed shortwave of category `n` in cell
      !! `(i, j)`.  Reads the state in place and writes ONLY `(i, j, n[, k])` of the outputs,
      !! each exactly once.  The indices are passed BY VALUE (nvfortran `do concurrent`
      !! index-by-reference bug).
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nilyr
      real(wp), intent(in) :: albicev, albicei, albsnowv, albsnowi
         !! Albedo knobs (`ice_sw_params_t`).
      real(wp), intent(in) :: fhtan
         !! `ice_albedo_fhtan(ahmax)`, once per configuration.
      real(wp), intent(in) :: snowpatch
         !! Snow-cover length scale (m).
      real(wp), intent(in) :: vis_frac
         !! `sw_vis_frac`.
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: tsfcn(nx, ny, ncat), sw_down(nx, ny)
      real(wp), intent(inout) :: alvn(nx, ny, ncat), alin(nx, ny, ncat)
      real(wp), intent(inout) :: fswsfcn(nx, ny, ncat), fswintn(nx, ny, ncat), fswthrun(nx, ny, ncat)
      real(wp), intent(inout) :: iswabsn(nx, ny, ncat, nilyr)
      integer(int32), intent(inout) :: errn(nx, ny, ncat)

      real(wp) :: a, hi, hs, sw, swv, swn
      real(wp) :: alv_i, ali_i, alv_s, ali_s, asnow, alv, ali
      real(wp) :: fsfc, fint, fthru
      real(wp) :: isw(ICE_NILYR_MAX)
      integer :: k
      logical :: fin

      a = aicen(i, j, n)
      sw = sw_down(i, j)
      fin = ieee_is_finite(a + vicen(i, j, n) + vsnon(i, j, n) + tsfcn(i, j, n) + sw)
      if (fin .and. a > ICE_SW_PUNY) then
         hi = vicen(i, j, n)/a
         hs = vsnon(i, j, n)/a
         call ice_albedos_ccsm3(hi, hs, tsfcn(i, j, n), albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, &
                                alv_i, ali_i, alv_s, ali_s, asnow, alv, ali)
         call ice_sw_band_split(sw, vis_frac, swv, swn)
         call ice_absorbed_ccsm3(nilyr, hi, asnow, swv, swn, alv_i, ali_i, alv_s, ali_s, fsfc, fint, fthru, isw)
         alvn(i, j, n) = alv
         alin(i, j, n) = ali
         fswsfcn(i, j, n) = fsfc
         fswintn(i, j, n) = fint
         fswthrun(i, j, n) = fthru
         do k = 1, nilyr
            iswabsn(i, j, n, k) = isw(k)
         end do
      else
         ! ice-free or corrupt input: the upstream defaults (ocean albedo, zero absorption)
         alvn(i, j, n) = ICE_ALBOCN
         alin(i, j, n) = ICE_ALBOCN
         fswsfcn(i, j, n) = 0.0_wp
         fswintn(i, j, n) = 0.0_wp
         fswthrun(i, j, n) = 0.0_wp
         do k = 1, nilyr
            iswabsn(i, j, n, k) = 0.0_wp
         end do
      end if
      errn(i, j, n) = merge(0_int32, ICE_ERR_NONFINITE, fin)
   end subroutine ice_radiation_cell

   subroutine ice_step_radiation_impl(nx, ny, ncat, nilyr, is, ie, js, je, &
                                      albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, vis_frac, &
                                      aicen, vicen, vsnon, tsfcn, sw_down, &
                                      alvn, alin, fswsfcn, fswintn, fswthrun, iswabsn, errn)
      !! ONE category-parallel launch over the cells `is:ie, js:je` (module docstring, piece 2).
      !! Every array must be device-present on an offload build.
      integer, intent(in) :: nx, ny, ncat, nilyr, is, ie, js, je
      real(wp), intent(in) :: albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, vis_frac
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: tsfcn(nx, ny, ncat), sw_down(nx, ny)
      real(wp), intent(inout) :: alvn(nx, ny, ncat), alin(nx, ny, ncat)
      real(wp), intent(inout) :: fswsfcn(nx, ny, ncat), fswintn(nx, ny, ncat), fswthrun(nx, ny, ncat)
      real(wp), intent(inout) :: iswabsn(nx, ny, ncat, nilyr)
      integer(int32), intent(inout) :: errn(nx, ny, ncat)
      integer :: i, j, n

      !$acc data present(aicen, vicen, vsnon, tsfcn, sw_down, alvn, alin) &
      !$acc&     present(fswsfcn, fswintn, fswthrun, iswabsn, errn)
      do concurrent(n=1:ncat, j=js:je, i=is:ie)
         call ice_radiation_cell(i, j, n, nx, ny, ncat, nilyr, &
                                 albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, vis_frac, &
                                 aicen, vicen, vsnon, tsfcn, sw_down, &
                                 alvn, alin, fswsfcn, fswintn, fswthrun, iswabsn, errn)
      end do
      !$acc end data
   end subroutine ice_step_radiation_impl

   subroutine ice_radiation_aggregate_impl(nx, ny, ncat, is, ie, js, je, aicen, alvn, alin, errn, &
                                           alv_ai, ali_ai, ice_err)
      !! ONE per-cell launch: the partial sums `alv_ai = sum_n aicen*alvn`, `ali_ai =
      !! sum_n aicen*alin` over categories with `aicen > puny` (fixed order; NOT divided by
      !! the ice area, WITHOUT open water: a grid-cell albedo is `alv_ai + aice0*(ocean
      !! albedo)`, an ice-only mean `alv_ai/aice`) and the merge of the per-category error
      !! words into `ice_err`.  Each work item writes only its own cell.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: aicen(nx, ny, ncat), alvn(nx, ny, ncat), alin(nx, ny, ncat)
      integer(int32), intent(in) :: errn(nx, ny, ncat)
      real(wp), intent(inout) :: alv_ai(nx, ny), ali_ai(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j, n
      real(wp) :: a, sv, si
      integer(int32) :: e

      !$acc data present(aicen, alvn, alin, errn, alv_ai, ali_ai, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, a, sv, si, e)
         sv = 0.0_wp
         si = 0.0_wp
         e = 0_int32
         do n = 1, ncat
            a = aicen(i, j, n)
            if (a > ICE_SW_PUNY) then
               sv = sv + alvn(i, j, n)*a
               si = si + alin(i, j, n)*a
            end if
            e = ior(e, errn(i, j, n))
         end do
         alv_ai(i, j) = sv
         ali_ai(i, j) = si
         ice_err(i, j) = ior(ice_err(i, j), e)
      end do
      !$acc end data
   end subroutine ice_radiation_aggregate_impl

   subroutine ice_step_radiation(rad, swp, is, ie, js, je, nx, ny, ncat, nilyr, aicen, vicen, vsnon, &
                                 tsfcn, sw_down, ice_err, ierr)
      !! Host shim: recompute the radiation slot `rad` from the current category state and the
      !! broadband downward shortwave `sw_down` over the cells `is:ie, js:je` (two launches).
      !! `swp` holds the albedo knobs and `sw_vis_frac`; `fhtan` and `snowpatch` come from
      !! `rad%configure`.  FAILS LOUD: `ierr /= ICE_RAD_OK` (and nothing runs) if the slot is not
      !! initialised or configured, or `nx, ny, ncat, nilyr` or the cell range disagree with the
      !! slot.  The slot and every input must be device-mapped on an offload build.
      type(ice_rad_state_t), intent(inout) :: rad
      type(ice_sw_params_t), intent(in) :: swp
      integer, intent(in) :: is, ie, js, je, nx, ny, ncat, nilyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: tsfcn(nx, ny, ncat), sw_down(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: ierr

      ierr = ICE_RAD_OK
      if (.not. rad%is_init) then
         ierr = ICE_RAD_ERR_NOT_INIT
         return
      end if
      if (.not. rad%configured) then
         ierr = ICE_RAD_ERR_NOT_CONFIGURED
         return
      end if
      if (nx /= rad%nx .or. ny /= rad%ny .or. ncat /= rad%ncat .or. nilyr /= rad%nilyr &
          .or. is < 1 .or. ie > nx .or. js < 1 .or. je > ny) then
         ierr = ICE_RAD_ERR_SHAPE
         return
      end if
      call ice_step_radiation_impl(nx, ny, ncat, nilyr, is, ie, js, je, &
                                   swp%albicev, swp%albicei, swp%albsnowv, swp%albsnowi, rad%fhtan, rad%snowpatch, &
                                   swp%sw_vis_frac, aicen, vicen, vsnon, tsfcn, sw_down, &
                                   rad%alvn, rad%alin, rad%fswsfcn, rad%fswintn, rad%fswthrun, rad%Iswabsn, rad%errn)
      call ice_radiation_aggregate_impl(nx, ny, ncat, is, ie, js, je, aicen, rad%alvn, rad%alin, rad%errn, &
                                        rad%alv_ai, rad%ali_ai, ice_err)
   end subroutine ice_step_radiation

end module rdb_ice_radiation
