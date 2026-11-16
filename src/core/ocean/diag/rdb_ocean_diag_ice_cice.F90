!! CICE-named sea-ice history diagnostics of the CICE/Icepack-port scheme (roundabout-original).
module rdb_ocean_diag_ice_cice
   !! The fills behind the CICE history names under `&ocean_ice_nml scheme = "cice"`.  They are STEP
   !! quantities: each fill reads what the LAST ice step left in the flux slot (`ice_flux_t`: the
   !! merged fluxes, the per-step `[m per step]` thickness bookkeeping, and the post-thermodynamics
   !! aggregate snapshot `aice_th, vice_th`), the radiation slot (albedo partial sums) and the
   !! end-of-step category state, and converts it to the CICE history unit.  The WINDOW MEAN over the
   !! output interval is the diagnostic manager's own `DIAG_OP_MEAN` fold (every fill runs every outer
   !! step, weighted by the step), so there is no second accumulator: nothing here is state, nothing is
   !! restart-carried, and nothing is written but the diagnostic buffer (the slots are read only).
   !! The ice steps every outer step (`ndtd = 1`, ice step = outer step), so the mean over the output
   !! window is the exact time mean of the per-step values.
   !!
   !! **Restart.**  The first output window after a resume is built from the steps taken after the
   !! resume only (the diagnostic manager's accumulators are not checkpointed, as for every `mean`
   !! diagnostic).  Each per-step value is a pure function of restart state, so the values themselves
   !! continue bitwise (`test_cice_diags`, `restart_fills_equal`).
   !!
   !! **Normalisation (the convention of every field here and of the SIS2 path's rates).**  Per unit
   !! GRID-CELL area, as CICE reports `congel`, `meltt`, `dvidtt` ...: a category quantity is weighted by
   !! the category's start-of-step area when the cell is merged (`ice_flux_merge_categories`).  The rates
   !! are in `cm/day` of ice (or snow/water) thickness per cell area, `[m per step] * 8.64e6 / dt`.
   !! The heat and mass fluxes are the cell sums.  CICE names the cell-weighted flux fields with the
   !! `_ai` suffix (its un-suffixed `fsens`, `flat`, ... are per unit ICE area, a divide by `aice` that
   !! cannot be taken on a window mean), so the flux diagnostics carry the `_ai` names.
   !!
   !! **Definitions** (Hunke et al., CICE documentation, "History"; Icepack documentation, "Thermodynamics"):
   !! | name | definition | unit, sign |
   !! |---|---|---|
   !! | `congel` | basal congelation ice growth | cm/day, >= 0 |
   !! | `frazil` | new ice made from the frazil bank in open water (`frazil_vol`) | cm/day, >= 0 |
   !! | `snoice` | snow-ice (flooding) formation | cm/day, >= 0 |
   !! | `meltt` | top ice melt | cm/day, >= 0 (melt is positive) |
   !! | `meltb` | basal ice melt | cm/day, >= 0 |
   !! | `meltl` | lateral ice melt | cm/day, >= 0 |
   !! | `dvidtt` | `(vice_after_thermo - vice_init)/dt`: thermodynamic ice-volume tendency | cm/day |
   !! | `dvidtd` | `(vice_end - vice_after_thermo)/dt`: dynamic (advection, ridging, rebin) tendency | cm/day |
   !! | `daidtt`, `daidtd` | the same for the cell ice area `aice` | %/day |
   !! | `Tsfc` | ice-area weighted category surface temperature plus the open water at the freezing point | degC |
   !! | `fswabs_ai` | shortwave absorbed by snow, ice AND the part passed to the ocean (`fswsfc + fswint + fswthru`) | W m-2, + down |
   !! | `fswthru_ai` | shortwave through the ice into the ocean | W m-2, + down |
   !! | `albsni` | broadband albedo, ice-area weighted `sum_n aicen (f_vis alvn + (1-f_vis) alin)` with `sw_vis_frac`; divide by the `aice` mean for the ice-only albedo | % |
   !! | `flwup_ai` | outgoing longwave net of the reflected incoming part (`-eps sigma T^4` form, NEGATIVE) | W m-2, + down |
   !! | `fsens_ai`, `flat_ai` | sensible and latent heat, atmosphere to surface | W m-2, + down |
   !! | `fhocn_ai` | heat from ice to ocean (never frazil: that is the ocean's bank) | W m-2, + into the ocean |
   !! | `fresh_ai` | fresh-water flux ice to ocean, water equivalent (`rhofresh = 1000`); contains `evap` | cm/day, + into the ocean |
   !! | `fsalt_ai` | salt flux ice to ocean | kg m-2 s-1, + into the ocean |
   !! | `evap_ai` | evaporative water flux, water equivalent | cm/day, + CONDENSATION (`< 0` is sublimation) |
   !!
   !! The thermodynamic volume identity is `dvidtt = congel + frazil + snoice - meltt - meltb - meltl`
   !! up to the ICE SUBLIMATION term (`evap` over ice, not a named melt in the history set) and the
   !! thin-category zap residual of the cleanup; the budget test in `test_cice_diags` measures the
   !! difference against the sublimation it should equal.  Not provided, because the column
   !! thermodynamics does not export them: `melts` (snow melt; computed per column, discarded before
   !! the merge) and `albice` (needs the bare-ice albedo split).
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_ocean_state, only: ocean_state_t
   implicit none
   private

   public :: cice_fill, cice_scheme_active
   public :: CD_CONGEL, CD_FRAZIL, CD_SNOICE, CD_MELTT, CD_MELTB, CD_MELTL, CD_DVIDTT, CD_DVIDTD
   public :: CD_DAIDTT, CD_DAIDTD, CD_TSFC, CD_FSWABS, CD_FLWUP, CD_FSENS, CD_FLAT, CD_FHOCN
   public :: CD_FSWTHRU, CD_FRESH, CD_FSALT, CD_EVAP, CD_ALBSNI
   public :: cice_scaled_impl, cice_dvidtt_impl, cice_diff_impl, cice_tsfc_impl, cice_albsni_impl

   integer, parameter :: CD_CONGEL = 1, CD_FRAZIL = 2, CD_SNOICE = 3, CD_MELTT = 4, CD_MELTB = 5
   integer, parameter :: CD_MELTL = 6, CD_DVIDTT = 7, CD_DVIDTD = 8, CD_DAIDTT = 9, CD_DAIDTD = 10
   integer, parameter :: CD_TSFC = 11, CD_FSWABS = 12, CD_FLWUP = 13, CD_FSENS = 14, CD_FLAT = 15
   integer, parameter :: CD_FHOCN = 16, CD_FSWTHRU = 17, CD_FRESH = 18, CD_FSALT = 19, CD_EVAP = 20
   integer, parameter :: CD_ALBSNI = 21
      !! Which diagnostic `cice_fill` computes.

   real(wp), parameter :: MPS_TO_CMPDAY = 8.64e6_wp
      !! m/s -> cm/day (86400 s/day x 100 cm/m); CICE `mps_to_cmpdy`.  Also `secday*c100`, the %/day
      !! factor of the area tendencies.
   real(wp), parameter :: RHOFRESH = 1000.0_wp
      !! kg m-3, the fresh-water density CICE uses to turn a mass flux into cm/day.

contains

   pure function cice_scheme_active(state) result(on)
      !! True when the CICE-port scheme is the running ice (its flux slot exists and is initialised).
      type(ocean_state_t), intent(in) :: state
      logical :: on
      on = state%icec%is_init .and. state%icec%flux%is_init .and. state%ice%cice%is_init
   end function cice_scheme_active

   subroutine cice_fill(state, id, buf)
      !! Fill `buf(:, :, 1)` with diagnostic `id` (`CD_*`) from the slots of the last ice step; NaN on
      !! land and ghosts, 0 before the first ice step (`dt_step = 0`).  Outer shim of the flat kernels.
      type(ocean_state_t), intent(in) :: state
      integer, intent(in) :: id
      real(wp), intent(inout) :: buf(:, :, :)
      real(wp) :: rate, wrate

      call cice_zero_impl(buf)
      if (.not. cice_scheme_active(state)) return
      rate = 0.0_wp
      if (state%icec%flux%dt_step > 0.0_wp) rate = MPS_TO_CMPDAY/state%icec%flux%dt_step
      wrate = 0.0_wp
      if (state%icec%flux%dt_step > 0.0_wp) wrate = MPS_TO_CMPDAY/RHOFRESH
      associate (f => state%icec%flux, c => state%ice%cice, w => state%metrics%wet_T)
         select case (id)
         case (CD_CONGEL)
            call cice_scaled_impl(f%nx, f%ny, w, f%congel, rate, buf)
         case (CD_FRAZIL)
            call cice_scaled_impl(f%nx, f%ny, w, f%frazil_vol, rate, buf)
         case (CD_SNOICE)
            call cice_scaled_impl(f%nx, f%ny, w, f%snoice, rate, buf)
         case (CD_MELTT)
            call cice_scaled_impl(f%nx, f%ny, w, f%meltt, rate, buf)
         case (CD_MELTB)
            call cice_scaled_impl(f%nx, f%ny, w, f%meltb, rate, buf)
         case (CD_MELTL)
            call cice_scaled_impl(f%nx, f%ny, w, f%meltl, rate, buf)
         case (CD_DVIDTT)
            call cice_dvidtt_impl(f%nx, f%ny, f%ncat, w, c%vicen_init, f%vice_th, rate, buf)
         case (CD_DVIDTD)
            call cice_diff_impl(f%nx, f%ny, w, c%vice, f%vice_th, rate, buf)
         case (CD_DAIDTT)
            call cice_diff_impl(f%nx, f%ny, w, f%aice_th, c%aice_init, rate, buf)
         case (CD_DAIDTD)
            call cice_diff_impl(f%nx, f%ny, w, c%aice, f%aice_th, rate, buf)
         case (CD_TSFC)
            call cice_tsfc_impl(f%nx, f%ny, f%ncat, w, c%p%c%puny, c%aicen, c%Tsfcn, state%icec%work%Tf, buf)
         case (CD_FSWABS)
            call cice_scaled_impl(f%nx, f%ny, w, f%fswabs, 1.0_wp, buf)
         case (CD_FLWUP)
            call cice_scaled_impl(f%nx, f%ny, w, f%flwout, 1.0_wp, buf)
         case (CD_FSENS)
            call cice_scaled_impl(f%nx, f%ny, w, f%fsens, 1.0_wp, buf)
         case (CD_FLAT)
            call cice_scaled_impl(f%nx, f%ny, w, f%flat, 1.0_wp, buf)
         case (CD_FHOCN)
            call cice_scaled_impl(f%nx, f%ny, w, f%fhocn, 1.0_wp, buf)
         case (CD_FSWTHRU)
            call cice_scaled_impl(f%nx, f%ny, w, f%fswthru, 1.0_wp, buf)
         case (CD_FRESH)
            ! a mass flux is already per second: no 1/dt
            call cice_scaled_impl(f%nx, f%ny, w, f%fresh, wrate, buf)
         case (CD_FSALT)
            call cice_scaled_impl(f%nx, f%ny, w, f%fsalt, 1.0_wp, buf)
         case (CD_EVAP)
            call cice_scaled_impl(f%nx, f%ny, w, f%evap, wrate, buf)
         case (CD_ALBSNI)
            call cice_albsni_impl(f%nx, f%ny, w, c%p%sw%sw_vis_frac, state%icec%rad%alv_ai, state%icec%rad%ali_ai, buf)
         case default
            continue   ! unknown id: the buffer stays zero
         end select
      end associate
   end subroutine cice_fill

   pure subroutine cice_zero_impl(buf)
      !! Device-side zero of the whole buffer (a host array assignment would miss the device copy).
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k

      do concurrent(k=1:size(buf, 3), j=1:size(buf, 2), i=1:size(buf, 1))
         buf(i, j, k) = 0.0_wp
      end do
   end subroutine cice_zero_impl

   pure subroutine cice_scaled_impl(nx, ny, wet_T, src, scale, buf)
      !! `buf = scale*src` on wet cells, NaN on land/ghost.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: wet_T(nx, ny), src(nx, ny)
      real(wp), intent(in) :: scale
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nxl, nyl
      real(wp) :: qnan

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:nyl, i=1:nxl)
         if (wet_T(i, j) > 0.5_wp) then
            buf(i, j, 1) = scale*src(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cice_scaled_impl

   pure subroutine cice_diff_impl(nx, ny, wet_T, a, b, scale, buf)
      !! `buf = scale*(a - b)` on wet cells, NaN on land/ghost.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: wet_T(nx, ny), a(nx, ny), b(nx, ny)
      real(wp), intent(in) :: scale
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nxl, nyl
      real(wp) :: qnan

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:nyl, i=1:nxl)
         if (wet_T(i, j) > 0.5_wp) then
            buf(i, j, 1) = scale*(a(i, j) - b(i, j))
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cice_diff_impl

   pure subroutine cice_dvidtt_impl(nx, ny, ncat, wet_T, vicen_init, vice_th, scale, buf)
      !! `buf = scale*(vice_th - sum_n vicen_init)`: the thermodynamic volume change of the step against
      !! the start-of-step category volumes (fixed category order).
      integer, intent(in) :: nx, ny, ncat
      real(wp), intent(in) :: wet_T(nx, ny), vicen_init(nx, ny, ncat), vice_th(nx, ny)
      real(wp), intent(in) :: scale
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, n, nxl, nyl
      real(wp) :: qnan, v0

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:nyl, i=1:nxl) local(n, v0)
         if (wet_T(i, j) > 0.5_wp) then
            v0 = 0.0_wp
            do n = 1, ncat
               v0 = v0 + vicen_init(i, j, n)
            end do
            buf(i, j, 1) = scale*(vice_th(i, j) - v0)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cice_dvidtt_impl

   pure subroutine cice_tsfc_impl(nx, ny, ncat, wet_T, puny, aicen, Tsfcn, Tf, buf)
      !! `Tsfc = sum_n aicen Tsfcn + max(0, 1 - sum_n aicen) Tf` over the categories with `aicen > puny`
      !! (CICE: "averaged with Tf if no ice is present"), degC; NaN on land/ghost.
      integer, intent(in) :: nx, ny, ncat
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: aicen(nx, ny, ncat), Tsfcn(nx, ny, ncat), Tf(nx, ny)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, n, nxl, nyl
      real(wp) :: qnan, tsum, asum

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:nyl, i=1:nxl) local(n, tsum, asum)
         if (wet_T(i, j) > 0.5_wp) then
            tsum = 0.0_wp
            asum = 0.0_wp
            do n = 1, ncat
               if (aicen(i, j, n) > puny) then
                  tsum = tsum + aicen(i, j, n)*Tsfcn(i, j, n)
                  asum = asum + aicen(i, j, n)
               end if
            end do
            buf(i, j, 1) = tsum + max(0.0_wp, 1.0_wp - asum)*Tf(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cice_tsfc_impl

   pure subroutine cice_albsni_impl(nx, ny, wet_T, vis_frac, alv_ai, ali_ai, buf)
      !! `buf = 100 (vis_frac alv_ai + (1 - vis_frac) ali_ai)`, percent, ice-area weighted (the radiation
      !! slot's partial sums `sum_n aicen albedo_n`); NaN on land/ghost.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: vis_frac
      real(wp), intent(in) :: wet_T(nx, ny), alv_ai(nx, ny), ali_ai(nx, ny)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nxl, nyl
      real(wp) :: qnan

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:nyl, i=1:nxl)
         if (wet_T(i, j) > 0.5_wp) then
            buf(i, j, 1) = 100.0_wp*(vis_frac*alv_ai(i, j) + (1.0_wp - vis_frac)*ali_ai(i, j))
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cice_albsni_impl

end module rdb_ocean_diag_ice_cice
