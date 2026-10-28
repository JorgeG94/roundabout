!! Ocean-side frazil accumulator (SIS2 port, PR 1).
module rdb_ice_frazil
   !! Clamp the ocean SURFACE layer at the seawater freezing point and
   !! BANK the removed supercooling deficit for the sea-ice model
   !! (`PLAN_SEA_ICE.md` "Prerequisites" #2 — the MOM6/SIS2 frazil
   !! limiter on the ocean side of the coupling seam).
   !!
   !! Physics: after a full outer dynamics step the surface layer
   !! (k = nz, Roundabout's bottom-up convention) can be advected/mixed
   !! below the salinity-dependent freezing point `T_f =
   !! eos_freezing_point(S, p)`.  Liquid seawater cannot supercool at
   !! leading order — the deficit freezes out as frazil crystals whose
   !! latent-heat release warms the water back to T_f.  This kernel
   !! applies exactly that: it raises the surface tracer temperature to
   !! T_f and deposits the heat spent,
   !!
   !!   deficit = ρ·Cp·h_surf·(T_f − T)⁺   [J/m²],
   !!
   !! into the persistent `frazil_heat` bank on the ice slot.  ENERGY
   !! CONSERVING by construction: the sensible heat ADDED to the ocean
   !! equals the latent heat the (future, PR 3) ice model will spend
   !! forming new ice from the same bank — banked, never discarded.
   !!
   !! The matching hTr increment `h·(T_f − T)` (K·m) is mirrored into
   !! the slot's `heat_budget_frazil` contributor so the global heat
   !! budget (`rdb_ocean_budgets`) still closes with ice on.
   !!
   !! Runs ONCE per outer step on the FINAL (post-RK2-average,
   !! post-remap) tracer state — clamping inside the RK2 stages would
   !! not bound the averaged result.  Gated at the driver on
   !! `state%ice%enable` (default off ⇒ byte-identical).
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED, H_DIV_EPS, GRAVITY, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   implicit none
   private

   public :: ice_frazil_accumulate
   public :: ice_frazil_accumulate_column

contains

   pure subroutine ice_frazil_accumulate(grid, eos, ms, frazil_heat, &
                                         heat_budget_frazil)
      !! Outer shim (outer-shim + flat-impl pattern): dereference the
      !! tracer registry (`ms%tracers(idx)%hTr`) on the HOST and forward
      !! bare arrays to the device kernel — NVHPC stdpar cannot follow
      !! the array-of-derived-types indirection inside a do-concurrent
      !! body.  No-op when either S or T is unregistered.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(inout) :: frazil_heat(:, :)
         !! Persistent supercooling bank (J/m²) on the ice slot —
         !! `ocean_sea_ice_t%frazil_heat`, device-mapped by the slot's
         !! `enter_data`.
      real(wp), intent(inout) :: heat_budget_frazil(:, :, :)
         !! Heat-budget contributor accumulator (K·m per cell) —
         !! `ocean_sea_ice_t%heat_budget_frazil`.

      integer :: idx_T, idx_S

      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (idx_T <= 0 .or. idx_S <= 0) return

      call ice_frazil_accumulate_impl(ms%tracers(idx_T)%hTr, &
                                      ms%tracers(idx_S)%hTr, &
                                      ms%h_layer, ms%wet_mask, &
                                      frazil_heat, heat_budget_frazil, &
                                      eos, grid%nghost, &
                                      ms%nz_ml, grid%nx_total, grid%ny_total)
   end subroutine ice_frazil_accumulate

   pure subroutine ice_frazil_accumulate_impl(hTr_T, hTr_S, h_layer, wet_mask, &
                                              frazil_heat, heat_budget_frazil, &
                                              eos, nghost, nz, nx, ny)
      !! Device kernel over PHYSICAL surface cells (ghosts excluded —
      !! wall ghosts are inert and seam ghosts are rebuilt by the wrap /
      !! exchange, so clamping them would double-count the bank in any
      !! area integral).  Per wet cell at k = nz:
      !!
      !!   T = hTr_T/h,  S = hTr_S/h,  T_f = eos_freezing_point(S, 0)
      !!   (OCEAN side: under `ts_convention = "cons_abs"` T, S, T_f are
      !!   CT, SA and the CT liquidus -- no conversion here; `Cp` is
      !!   `eos%c_p`, `&ocean_eos_nml c_p`, default 3992.)
      !!   if T < T_f:  frazil_heat += ρ·Cp·h·(T_f − T)  and
      !!                hTr_T = h·T_f   (exact clamp)
      !!
      !! Surface pressure p = 0 (surface-relative hydrostatic zero —
      !! same convention as the surface-flux kernels; the in-situ
      !! layer-centre depression over h_surf/2 is O(1 mK) and ignored,
      !! matching MOM6's p=0 surface-frazil evaluation).  Vanished
      !! columns (`h <= H_VANISHED`) and land (`wet_mask = 0`) are
      !! skipped — the division for a PHYSICAL T/S gates on
      !! `H_VANISHED` per the thin-layer taxonomy.  The conditional
      !! write uses an inner `if`, NOT a `do concurrent` mask header
      !! (masked DC headers fail tools/dc_audit.py --strict).
      integer, intent(in) :: nghost, nz, nx, ny
      real(wp), intent(inout) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(inout) :: frazil_heat(nx, ny)
      real(wp), intent(inout) :: heat_budget_frazil(nx, ny, 1)
      type(eos_t), intent(in) :: eos

      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h, t_sfc, s_sfc, t_f, d_htr

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(h, t_sfc, s_sfc, t_f, d_htr)
         h = h_layer(i, j, nz)
         if (wet_mask(i, j) > 0.5_wp .and. h > H_VANISHED) then
            t_sfc = hTr_T(i, j, nz)/h
            s_sfc = hTr_S(i, j, nz)/h
            t_f = eos_freezing_point(eos, s_sfc, 0.0_wp)
            if (t_sfc < t_f) then
               d_htr = h*(t_f - t_sfc)
               frazil_heat(i, j) = frazil_heat(i, j) &
                                   + RHO_WATER*eos%c_p*d_htr
               heat_budget_frazil(i, j, 1) = heat_budget_frazil(i, j, 1) + d_htr
               ! Exact clamp: T_after = (h·T_f)/h == T_f.
               hTr_T(i, j, nz) = h*t_f
            end if
         end if
      end do
   end subroutine ice_frazil_accumulate_impl

   pure subroutine ice_frazil_accumulate_column(grid, eos, ms, frazil_heat, &
                                                heat_budget_frazil, &
                                                pressure_dependent, reclaim)
      !! Full-column, optionally pressure-dependent frazil with reclaim
      !! (`&ocean_ice_nml frazil_column`; OM3 `FRAZIL` +
      !! `PRESSURE_DEPENDENT_FRAZIL` + `RECLAIM_FRAZIL`).  Outer shim:
      !! dereferences the tracer registry on the host and forwards bare
      !! arrays to `ice_frazil_accumulate_column_impl`.  No-op when S or
      !! T is unregistered.  The surface-only `ice_frazil_accumulate` is
      !! untouched and remains the default path.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(inout) :: frazil_heat(:, :)
         !! Persistent frazil bank F (J/m2), `ocean_sea_ice_t%frazil_heat`.
      real(wp), intent(inout) :: heat_budget_frazil(:, :, :)
         !! Heat-budget contributor (K.m per cell).
      logical, intent(in) :: pressure_dependent
         !! Evaluate T_f at the layer-centre hydrostatic pressure.
      logical, intent(in) :: reclaim
         !! Spend the standing bank F on a warm surface layer first.

      integer :: idx_T, idx_S

      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (idx_T <= 0 .or. idx_S <= 0) return

      call ice_frazil_accumulate_column_impl(ms%tracers(idx_T)%hTr, &
                                             ms%tracers(idx_S)%hTr, &
                                             ms%h_layer, ms%wet_mask, &
                                             frazil_heat, heat_budget_frazil, &
                                             eos, pressure_dependent, reclaim, &
                                             grid%nghost, ms%nz_ml, &
                                             grid%nx_total, grid%ny_total)
   end subroutine ice_frazil_accumulate_column

   pure subroutine ice_frazil_accumulate_column_impl(hTr_T, hTr_S, h_layer, wet_mask, &
                                                     frazil_heat, heat_budget_frazil, &
                                                     eos, pressure_dependent, reclaim, &
                                                     nghost, nz, nx, ny)
      !! Device kernel, one thread per PHYSICAL wet column, serial in k,
      !! scalars only.  T, S, T_f are CT, SA and the CT liquidus under
      !! `cons_abs` (potential T/practical S otherwise); `hc = rho c_p h`
      !! with `rho = RHO_WATER`, `c_p = eos%c_p` (the same heat weight as
      !! the surface-only kernel and the console `src` term).
      !!
      !! Pressure.  With `pressure_dependent`, the layer-centre pressure
      !! is `p_k = g rho0 (sum_{m>k} h_m + h_k/2)` (Pa, surface pressure
      !! zero, `rho0 = eos%rho0`), summed top-down exactly as
      !! `p = g rho0 h_nz/2`, then `+= g rho0 (h_{m+1} + h_m)/2`.
      !! Otherwise p = 0.
      !!
      !! 1. Reclaim (`reclaim`).  If the bank `F = frazil_heat > 0` and the
      !!    surface layer is above its T_f: with `hc` of the surface
      !!    layer, `F <= hc (T - T_f)` cools it by `F/hc` and zeroes F;
      !!    otherwise `F -= hc (T - T_f)` and `T = T_f`.
      !! 2. Bottom-up sweep (k = 1 bed ... nz surface) with a carry
      !!    `c = 0` (J/m2).  A layer is visited when `T < 1 degC` (the
      !!    cheap pre-check; every liquidus used is below that) or `c > 0`.
      !!    A VANISHED layer (`h <= H_VANISHED`, the project dynamic-vanish
      !!    threshold; layers at `h <= H_DIV_EPS` are skipped outright)
      !!    can only freeze: `T < T_f` adds `hc (T_f - T)` to c and sets
      !!    `T = T_f`.  Any other layer, if `c + hc (T_f - T) < 0` (warm
      !!    enough to absorb the carry), cools by `c/hc` and zeroes c;
      !!    else `c += hc (T_f - T)` and `T = T_f`.  Frazil formed at depth
      !!    therefore rises and is re-melted by warmer layers above; only
      !!    what survives the top leaves the column.
      !! 3. Hand-off.  `F += c` (J/m2, positive = heat the ice model must
      !!    spend, same sign and units as the surface-only kernel).
      !!
      !! Conservation.  The OCEAN heat change of the column is exactly
      !! `c_final - F_reclaimed` (J/m2); the kernel accumulates the true
      !! `sum(Delta hTr_T)` (K.m) into `heat_budget_frazil`, so the
      !! console `src` term (`rho * budget`) matches the bank change.
      integer, intent(in) :: nghost, nz, nx, ny
      real(wp), intent(inout) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(inout) :: frazil_heat(nx, ny)
      real(wp), intent(inout) :: heat_budget_frazil(nx, ny, 1)
      type(eos_t), intent(in) :: eos
      logical, intent(in) :: pressure_dependent, reclaim

      integer :: i, j, k, m, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h, t_k, s_k, t_f, p_k, hc, c, f_bank, dsum, gr, new_t

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      gr = GRAVITY*eos%rho0

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(k, m, h, t_k, s_k, t_f, p_k, hc, c, f_bank, dsum, new_t)
         if (wet_mask(i, j) > 0.5_wp .and. rdb_vl_is_live(h_layer(i, j, nz))) then
            c = 0.0_wp
            dsum = 0.0_wp
            f_bank = frazil_heat(i, j)

            ! ---- 1. reclaim on the surface layer ----
            if (reclaim .and. f_bank > 0.0_wp) then
               h = h_layer(i, j, nz)
               t_k = rdb_vl_conc(hTr_T(i, j, nz), h)
               s_k = rdb_vl_conc(hTr_S(i, j, nz), h)
               p_k = 0.0_wp
               if (pressure_dependent) p_k = 0.5_wp*gr*h
               t_f = eos_freezing_point(eos, s_k, p_k)
               if (t_k > t_f) then
                  hc = RHO_WATER*eos%c_p*h
                  if (f_bank <= hc*(t_k - t_f)) then
                     new_t = t_k - f_bank/hc
                     f_bank = 0.0_wp
                  else
                     f_bank = f_bank - hc*(t_k - t_f)
                     new_t = t_f
                  end if
                  dsum = dsum + h*(new_t - t_k)
                  hTr_T(i, j, nz) = h*new_t
               end if
            end if

            ! ---- 2. bottom-up sweep with the carry ----
            do k = 1, nz
               h = h_layer(i, j, k)
               if (h <= H_DIV_EPS) cycle
               t_k = rdb_vl_conc(hTr_T(i, j, k), h)
               if (t_k < 1.0_wp .or. c > 0.0_wp) then
                  s_k = rdb_vl_conc(hTr_S(i, j, k), h)
                  p_k = 0.0_wp
                  if (pressure_dependent) then
                     p_k = 0.5_wp*gr*h_layer(i, j, nz)
                     do m = nz - 1, k, -1
                        p_k = p_k + 0.5_wp*gr*(h_layer(i, j, m + 1) + h_layer(i, j, m))
                     end do
                  end if
                  t_f = eos_freezing_point(eos, s_k, p_k)
                  hc = RHO_WATER*eos%c_p*h
                  if (.not. rdb_vl_is_live(h)) then
                     if (t_k < t_f) then
                        c = c + hc*(t_f - t_k)
                        dsum = dsum + h*(t_f - t_k)
                        hTr_T(i, j, k) = h*t_f
                     end if
                  else if (c > 0.0_wp .or. t_k < t_f) then
                     if (c + hc*(t_f - t_k) < 0.0_wp) then
                        new_t = t_k - c/hc
                        c = 0.0_wp
                     else
                        c = c + hc*(t_f - t_k)
                        new_t = t_f
                     end if
                     dsum = dsum + h*(new_t - t_k)
                     hTr_T(i, j, k) = h*new_t
                  end if
               end if
            end do

            ! ---- 3. hand-off ----
            frazil_heat(i, j) = f_bank + c
            heat_budget_frazil(i, j, 1) = heat_budget_frazil(i, j, 1) + dsum
         end if
      end do
   end subroutine ice_frazil_accumulate_column_impl

#include "rdb_vanished_layer.inc"

end module rdb_ice_frazil
