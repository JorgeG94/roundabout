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
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   implicit none
   private

   public :: ice_frazil_accumulate

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
      !! (masked DC headers break tools/dc_to_omp.py).
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
                                   + RHO_WATER*SEAWATER_CP*d_htr
               heat_budget_frazil(i, j, 1) = heat_budget_frazil(i, j, 1) + d_htr
               ! Exact clamp: T_after = (h·T_f)/h == T_f.
               hTr_T(i, j, nz) = h*t_f
            end if
         end if
      end do
   end subroutine ice_frazil_accumulate_impl

end module rdb_ice_frazil
