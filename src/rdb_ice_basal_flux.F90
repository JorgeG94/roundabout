!! Ocean -> ice basal heat flux fb from above-freezing surface heat (PR 3c).
module rdb_ice_basal_flux
   !! fb = RHO_WATER*SEAWATER_CP*max(0, SST - T_f)*h_top/dt_therm [W/m^2],
   !! the complement of the frazil bank: a warm (above-freezing) ocean surface
   !! under ice melts the base (fb > 0 -> bmelt in ice_temp_sis2 TRAP #3);
   !! a supercooled surface grows it (that path is the frazil bank, PR 1/3b).
   !! Same RHO_WATER*SEAWATER_CP convention as ice_frazil_accumulate so growth
   !! and melt share one energy scale.
   !!
   !! NO CAP (v1): fb is the full above-freezing flux. The column's bottom-melt
   !! peel is self-limiting (clamps to available ice mass, spills the remainder
   !! to heat_to_ocn), so nothing is discarded and the heat budget closes.
   !!
   !! ICE-PRESENCE GATE (crucial): fb is ZERO on ice-free cells — no ice base,
   !! no basal flux. The gate mirrors ice_thermo_columns' own ice threshold
   !! `sum_cat m_ice > ICE_RHO_ICE*H_VANISHED` so the column and the coupler
   !! agree EXACTLY on which cells exchange. Without it, a warm ice-free ocean
   !! cell (the normal open-ocean state, SST > T_f) would compute a large
   !! fb > 0 that the column then never consumes (m_ice = 0 => heat_to_ocn = 0),
   !! and the melt-side reduce kernel would inject a spurious Q_heat = -fb,
   !! cooling the open ocean by several degC per thermo step. The sample seam
   !! (sst_seam/ssurf_seam/tfw_seam) is STILL filled on all wet cells (harmless
   !! — the column only reads it where it has ice).
   !!
   !! One-step lag: reads the outer step's FINAL surface state (called post-dyn,
   !! pre-column, like ice_frazil_accumulate). Physical cells only; wet +
   !! non-vanished gate. Outer-shim + flat-impl (registry deref on host).
   !!
   !! ALSO fills the sample seam (sst_seam/ssurf_seam/tfw_seam) the column
   !! driver (rdb_ice_thermo_driver) reuses — one SST/SSS/T_f sample serves
   !! both fb and the column's ocean-side inputs, keeping them at the same
   !! one-step-lagged snapshot (physically consistent, PLAN_ICE_PR3c
   !! §"Ocean -> ice basal heat flux").
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private
   public :: ice_compute_basal_flux
contains
   pure subroutine ice_compute_basal_flux(grid, eos, ms, ice, dt_therm)
      !! Outer shim (outer-shim + flat-impl pattern): dereference the
      !! tracer registry (`ms%tracers(idx)%hTr`) on the HOST and forward
      !! bare arrays to the device kernel — same rule as `ice_frazil_accumulate`
      !! / `ice_frazil_uptake`. No-op when either S or T is unregistered.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(in) :: ms
         !! READ-ONLY: SST/SSS are sampled here, never written.
      type(ocean_sea_ice_t), intent(inout) :: ice
         !! Writes ice%fb and the sst_seam/ssurf_seam/tfw_seam sample seam.
      real(wp), intent(in) :: dt_therm
         !! Effective thermo timestep (s) — `ocean_dyn_t%therm_dt(dt)`.
      integer :: idx_T, idx_S
      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (idx_T <= 0 .or. idx_S <= 0) return
      call ice_compute_basal_flux_impl(ms%tracers(idx_T)%hTr, ms%tracers(idx_S)%hTr, &
                                       ms%h_layer, ms%wet_mask, ice%m_ice, eos, &
                                       ice%fb, ice%sst_seam, ice%ssurf_seam, ice%tfw_seam, &
                                       dt_therm, grid%nghost, ice%ncat, ms%nz_ml, &
                                       grid%nx_total, grid%ny_total)
   end subroutine ice_compute_basal_flux

   pure subroutine ice_compute_basal_flux_impl(hTr_T, hTr_S, h_layer, wet_mask, m_ice, eos, &
                                               fb, sst_seam, ssurf_seam, tfw_seam, &
                                               dt_therm, nghost, ncat, nz, nx, ny)
      !! Device kernel over PHYSICAL cells (ghosts excluded — same
      !! physical-cells-only contract as `ice_frazil_accumulate_impl`).
      !! `fb`/`sst_seam`/`ssurf_seam`/`tfw_seam` are zeroed unconditionally
      !! first. The SAMPLE seam (sst_seam/ssurf_seam/tfw_seam) is filled on
      !! every wet, non-vanished cell (harmless — the column reads it only
      !! where it has ice). But `fb` is filled ONLY where BOTH the cell is
      !! wet+non-vanished AND it carries ice (`sum_cat m_ice >
      !! ICE_RHO_ICE*H_VANISHED`, exactly ice_thermo_columns' own per-cat
      !! ice threshold, summed): no ice base ⇒ no basal flux. This keeps the
      !! coupler and the column in lockstep on which cells exchange, so an
      !! ice-free warm ocean cell (SST > T_f, the normal open-ocean state)
      !! reports fb = 0 and the melt-side reduce kernel's `-fb` term is a
      !! harmless subtraction of zero (rather than a spurious ocean-cooling
      !! Q_heat = -fb).
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if` gate only (never a masked
      !! `do concurrent` header).
      integer, intent(in) :: nghost, ncat, nz, nx, ny
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: fb(nx, ny)
      real(wp), intent(inout) :: sst_seam(nx, ny)
      real(wp), intent(inout) :: ssurf_seam(nx, ny)
      real(wp), intent(inout) :: tfw_seam(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h, sst, s_surf, tfw, m_ice_tot

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(h, sst, s_surf, tfw, cat, m_ice_tot)
         fb(i, j) = 0.0_wp
         sst_seam(i, j) = 0.0_wp
         ssurf_seam(i, j) = 0.0_wp
         tfw_seam(i, j) = 0.0_wp
         h = h_layer(i, j, nz)
         if (wet_mask(i, j) > 0.5_wp .and. h > H_VANISHED) then
            sst = hTr_T(i, j, nz)/h
            s_surf = hTr_S(i, j, nz)/h
            tfw = eos_freezing_point(eos, s_surf, 0.0_wp)
            sst_seam(i, j) = sst
            ssurf_seam(i, j) = s_surf
            tfw_seam(i, j) = tfw
            m_ice_tot = 0.0_wp
            do cat = 1, ncat
               m_ice_tot = m_ice_tot + m_ice(i, j, cat)
            end do
            if (m_ice_tot > ICE_RHO_ICE*H_VANISHED) then
               fb(i, j) = RHO_WATER*SEAWATER_CP*max(0.0_wp, sst - tfw)*h/dt_therm
            end if
         end if
      end do
   end subroutine ice_compute_basal_flux_impl

end module rdb_ice_basal_flux
