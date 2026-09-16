!! Driver-facing ice thermodynamics: run the column step, reduce the
!! per-(cell,cat) exchange into the ocean-coupling seam diags (PR 3c).
module rdb_ice_thermo_driver
   !! Bridges ice_thermo_columns (PR 3a) into the ocean-coupling seam:
   !! runs the column, then reduces h2o_ocn_to_ice/h2o_ice_to_ocn/heat_to_ocn
   !! over categories into the per-cell heat_flux_diag/m_melt_diag and the
   !! net-melt salt contribution added into salt_flux_diag.
   !!   Q_heat  = sum_cat(heat_to_ocn)/dt_therm - fb        [W/m^2, +down]
   !!   m_net   = sum_cat(h2o_ocn_to_ice - h2o_ice_to_ocn)  [kg/m^2, +freeze]
   !!   salt   += m_net*(s_surf - ICE_BULK_SALINITY)/dt_therm  (ADD, not set)
   !!   sw_thru_diag = sum_cat(sw_thru)                     [W/m^2, +down]
   !! PR 31: `sw_thru` (the shortwave penetrating the ice to the water) is
   !! now reduced to the per-cell `sw_thru_diag` and coupled to the ocean
   !! by `ice_ocean_sw_flux` (`rdb_ice_ocean_coupler`).  It is DELIBERATELY
   !! kept OUT of the `heat_flux_diag` sum above — `heat_flux_diag` carries
   !! only the NON-shortwave heat; the coupler delivers the shortwave via a
   !! separate path (the `q_sw` component, or a direct `Q_heat` add) so the
   !! energy reaches `Q_heat` exactly once in either component mode.
   !!
   !! **Ordering contract** (driver-enforced, see rdb_driver.F90): this
   !! module must run AFTER `ice_frazil_uptake` in the thermo-cadence
   !! block. `ice_frazil_uptake_impl` OVERWRITES `salt_flux_diag` at loop
   !! top (unconditional zero, then a gated write) — this module's
   !! reduction kernel ADDS its net-melt term on top of that write, never
   !! zeroing `salt_flux_diag` itself. `heat_flux_diag`/`m_melt_diag` ARE
   !! zeroed unconditionally here (this module owns them outright).
   !!
   !! **fb double-count guard**: the column already folds `fb` into
   !! `bmelt` (rdb_ice_column TRAP #3) and reports `heat_to_ocn`/
   !! `h2o_ice_to_ocn` net of that. The only place `fb` re-enters the
   !! energy accounting is the `Q_heat` formula above — never add it
   !! anywhere else.
   !!
   !! Reads the sst/s_surf/tfw sample seam filled by
   !! `rdb_ice_basal_flux%ice_compute_basal_flux` earlier in the same
   !! thermo-cadence block (one sample, reused — see that module's
   !! docstring).
   !!
   !! **PR 4a multicat dispatch**: `ncat == 1` runs the EXISTING
   !! `ice_thermo_driver_reduce_impl` byte-for-byte UNCHANGED (bit-identity
   !! contract). `ncat > 1` additionally fills `ice%fb_part_sum` — the
   !! fb-charged ice-cover-fraction snapshot — in a small DC kernel BEFORE
   !! `ice_thermo_columns` runs (using the SAME entry gate the column
   !! itself uses, `m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`, evaluated
   !! PRE-column so a category that melts out entirely this step is still
   !! counted in the weight it was actually charged `fb` under — counting
   !! POST-column would leak `part*fb*dt_therm` of energy from the books
   !! for exactly that category), then dispatches to
   !! `ice_thermo_driver_reduce_multicat_impl`, which part-weights the
   !! per-category sums and subtracts `fb_part_sum*fb` (rather than the
   !! bare `fb` the ncat==1 path subtracts) from the heat diag.
   !!
   !! **PR 26 snowfall**: when `ice%has_snowfall`, `ice_snow_part_ocn_
   !! fill_impl` runs (same PRE-column-snapshot contract and reasoning as
   !! `ice_fb_part_sum_fill_impl` — a category that melts out this window
   !! must still be counted in the cover it actually caught snow under)
   !! for BOTH `ncat==1` and `ncat>1` (unlike `fb_part_sum`, which is
   !! `ncat>1`-only: at `ncat==1` the ice-free share is a binary 0/1, not
   !! a part-weighted sum, so there is no "skip at ncat==1" case here).
   !! `ice%atm_fprec` is threaded into `ice_thermo_columns` unconditionally
   !! (harmless zeros when `has_snowfall` is false). The two reduce impls
   !! above are NOT touched by PR 26 — the ocean-bound share is delivered
   !! separately by `rdb_ice_snow%ice_snowfall_ocean_share`, called by the
   !! driver AFTER this step (mandated order, `rdb_driver.F90`).
   !!
   !! **PR 27 snow-ice flooding**: `ice%snow_ice` is threaded into
   !! `ice_thermo_columns` as `do_snow_ice` unconditionally (harmless
   !! `.false.` when the knob is off); the flood is entirely
   !! column-internal (SIS2's `ice_resize_SIS2`, no ocean mass/heat/salt
   !! exchange — see `rdb_ice_mass%ice_snow_ice_flood`'s docstring) and
   !! therefore does NOT enter the `Q_heat`/`m_net`/`salt` reduction
   !! above; `ice%snow_to_ice` is filled but not yet reduced/coupled
   !! (same "filled, not consumed" contract as `sw_thru`).
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ice_column, only: ice_thermo_columns, ICE_BULK_SALINITY, ICE_RHO_ICE
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private
   public :: ice_thermo_driver_step
contains

   pure subroutine ice_thermo_driver_step(grid, eos, ms, ice, dt_therm)
      !! Outer shim (outer-shim + flat-impl pattern): `ms`/`eos` are
      !! accepted for call-site parity with the other thermo-cadence
      !! kernels (`ice_compute_basal_flux`, `ice_frazil_uptake`) even
      !! though this step reads its ocean-surface sample from the
      !! `ice%sst_seam`/`ssurf_seam`/`tfw_seam` scratch (filled by
      !! `ice_compute_basal_flux` earlier in the same window) rather
      !! than re-deriving it from `ms`/`eos` directly — `eos` is unused
      !! here but kept in the signature so a future revision that DOES
      !! need a fresh EOS evaluation does not have to change every call
      !! site.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt_therm
         !! Effective thermo timestep (s) — `ocean_dyn_t%therm_dt(dt)`.

      associate (unused => eos)
      end associate

      if (ice%ncat > 1) then
         call ice_fb_part_sum_fill_impl(ms%wet_mask, ice%part_size, ice%m_ice, ice%fb_part_sum, &
                                        grid%nghost, ice%ncat, grid%nx_total, grid%ny_total)
      end if

      if (ice%has_snowfall) then
         call ice_snow_part_ocn_fill_impl(ms%wet_mask, ice%part_size, ice%m_ice, &
                                          ice%snow_part_ocn, grid%nghost, ice%ncat, &
                                          grid%nx_total, grid%ny_total)
      end if

      call ice_thermo_columns(grid%nghost, grid%nx_total, grid%ny_total, &
                              ice%ncat, ice%nk_ice, dt_therm, ice%snow_ice, ms%wet_mask, &
                              ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                              ice%sal_ice, ice%atm_sf0, ice%atm_dsfdt, ice%atm_sw_dn, &
                              ice%atm_fprec, ice%tfw_seam, ice%fb, ice%sst_seam, &
                              ice%ssurf_seam, ice%tsurf_out, ice%h2o_ocn_to_ice, &
                              ice%h2o_ice_to_ocn, ice%heat_to_ocn, ice%sw_thru, &
                              ice%snow_to_ice)

      if (ice%ncat == 1) then
         call ice_thermo_driver_reduce_impl(ms%wet_mask, ice%h2o_ocn_to_ice, &
                                            ice%h2o_ice_to_ocn, ice%heat_to_ocn, &
                                            ice%sw_thru, ice%ssurf_seam, ice%fb, &
                                            ice%heat_flux_diag, ice%sw_thru_diag, &
                                            ice%m_melt_diag, &
                                            ice%salt_flux_diag, dt_therm, grid%nghost, &
                                            ice%ncat, grid%nx_total, grid%ny_total)
      else
         call ice_thermo_driver_reduce_multicat_impl(ms%wet_mask, ice%part_size, &
                                                     ice%h2o_ocn_to_ice, ice%h2o_ice_to_ocn, &
                                                     ice%heat_to_ocn, ice%sw_thru, &
                                                     ice%ssurf_seam, ice%fb, &
                                                     ice%fb_part_sum, &
                                                     ice%heat_flux_diag, ice%sw_thru_diag, &
                                                     ice%m_melt_diag, &
                                                     ice%salt_flux_diag, dt_therm, grid%nghost, &
                                                     ice%ncat, grid%nx_total, grid%ny_total)
      end if
   end subroutine ice_thermo_driver_step

   pure subroutine ice_fb_part_sum_fill_impl(wet_mask, part_size, m_ice, fb_part_sum, &
                                             nghost, ncat, nx, ny)
      !! Fill `fb_part_sum(i,j) = Σ_c part_size(i,j,c)` restricted to
      !! categories passing the column's OWN entry gate
      !! (`m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`) — MUST run BEFORE
      !! `ice_thermo_columns` mutates `m_ice` (module docstring). Device
      !! kernel over PHYSICAL cells, inner `if`/serial `do cat` — never a
      !! masked `do concurrent` header.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them.
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: fb_part_sum(nx, ny)

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(cat)
         fb_part_sum(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            do cat = 1, ncat
               if (m_ice(i, j, cat) > ICE_RHO_ICE*H_VANISHED) then
                  fb_part_sum(i, j) = fb_part_sum(i, j) + part_size(i, j, cat)
               end if
            end do
         end if
      end do
   end subroutine ice_fb_part_sum_fill_impl

   pure subroutine ice_snow_part_ocn_fill_impl(wet_mask, part_size, m_ice, snow_part_ocn, &
                                               nghost, ncat, nx, ny)
      !! PR 26: fill `snow_part_ocn(i,j)` — the ice-FREE share of the cell
      !! that a uniform snowfall lands on with no ice underneath — under
      !! the SAME PRE-column-snapshot contract as `ice_fb_part_sum_fill_
      !! impl` (module docstring, this module docstring's PR-26
      !! paragraph): MUST run BEFORE `ice_thermo_columns` mutates `m_ice`,
      !! using the column's OWN entry gate
      !! (`m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`).
      !!
      !! Two-mode dispatch (`rdb_ice_state`'s per-CELL vs per-ICE-AREA
      !! convention, load-bearing): at `ncat==1` (legacy lumped, `part_size`
      !! never maintained) the ice-free share is BINARY — 0 if the cell's
      !! one category has ice, 1 if it does not (`ice_cell_concentration_
      !! impl`'s own `ci = 1 iff m_ice(1) > 0` convention, mirrored here).
      !! At `ncat>1` (SIS2 ITD, `part_size` live) it is `1 -
      !! Σ_{c : gate passes} part_size(c)`, clamped at 0 (categories that
      !! fail the gate but still carry `part_size>0` — SIS2's orphan-snow
      !! trap, PLAN_PR26_snowfall.md §11.1 — count toward the ocean share,
      !! never toward the column). The `ncat` branch is loop-invariant
      !! (same for every cell this call) so it stays INSIDE the one `do
      !! concurrent` (`dc-uniform-case` — splitting it would double the
      !! launch count for no benefit).
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if`/serial `do cat` — never a
      !! masked `do concurrent` header.
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: snow_part_ocn(nx, ny)

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi
      real(wp) :: cover

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(cat, cover)
         snow_part_ocn(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            if (ncat == 1) then
               if (m_ice(i, j, 1) > ICE_RHO_ICE*H_VANISHED) then
                  snow_part_ocn(i, j) = 0.0_wp
               else
                  snow_part_ocn(i, j) = 1.0_wp
               end if
            else
               cover = 0.0_wp
               do cat = 1, ncat
                  if (m_ice(i, j, cat) > ICE_RHO_ICE*H_VANISHED) then
                     cover = cover + part_size(i, j, cat)
                  end if
               end do
               snow_part_ocn(i, j) = max(0.0_wp, 1.0_wp - cover)
            end if
         end if
      end do
   end subroutine ice_snow_part_ocn_fill_impl

   pure subroutine ice_thermo_driver_reduce_impl(wet_mask, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                                 heat_to_ocn, sw_thru, ssurf_seam, fb, &
                                                 heat_flux_diag, sw_thru_diag, m_melt_diag, &
                                                 salt_flux_diag, &
                                                 dt_therm, nghost, ncat, nx, ny)
      !! Device kernel over PHYSICAL cells (ghosts excluded). Reduces the
      !! per-category column outputs into the per-cell coupling diags.
      !! `heat_flux_diag`/`sw_thru_diag`/`m_melt_diag` are zeroed
      !! unconditionally at loop top (this module owns them outright);
      !! `salt_flux_diag` is NOT zeroed — `ice_frazil_uptake` (which runs
      !! first this window, see the module docstring) already zeroed +
      !! wrote it, and this kernel ADDS its net-melt contribution on top,
      !! gated on wet to avoid land noise (an unwetted cell adds exactly 0).
      !!
      !! PR 31: `sw_thru_diag(i,j) = Σ_cat sw_thru(i,j,cat)` (W/m^2, +down)
      !! — the lumped ncat==1 reduction (no part-weight, matching
      !! `heat_flux_diag`). It is written to its OWN field, NEVER folded
      !! into `heat_flux_diag`: `heat_flux_diag` stays the non-shortwave
      !! heat share, and `ice_ocean_sw_flux` delivers `sw_thru_diag` to
      !! `Q_heat` exactly once (see that routine for the no-double-count
      !! argument across both component modes).
      !!
      !! ICE-FREE CELLS ARE A CLEAN NO-OP: on an ice-free wet cell the column
      !! zeroes all its per-cat outputs (h2o_*/heat_to_ocn = 0), and
      !! `ice_compute_basal_flux` gates `fb = 0` there (no ice base ⇒ no basal
      !! flux — the same ice-presence threshold the column uses). So
      !! `heat_flux_diag = 0/dt - 0 = 0`, `m_melt_diag = 0`, and the net-melt
      !! salt add is `0*(...) = 0`. The whole open-ocean surface therefore sees
      !! Q_heat = Q_salt = 0 from the ice path, as physics requires.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if`/serial `do cat` — never a
      !! masked `do concurrent` header.
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: h2o_ocn_to_ice(nx, ny, ncat)
      real(wp), intent(in) :: h2o_ice_to_ocn(nx, ny, ncat)
      real(wp), intent(in) :: heat_to_ocn(nx, ny, ncat)
      real(wp), intent(in) :: sw_thru(nx, ny, ncat)
      real(wp), intent(in) :: ssurf_seam(nx, ny)
      real(wp), intent(in) :: fb(nx, ny)
      real(wp), intent(inout) :: heat_flux_diag(nx, ny)
      real(wp), intent(inout) :: sw_thru_diag(nx, ny)
      real(wp), intent(inout) :: m_melt_diag(nx, ny)
      real(wp), intent(inout) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi
      real(wp) :: sum_ocn, sum_ice, sum_heat, sum_sw, m_net

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(cat, sum_ocn, sum_ice, sum_heat, sum_sw, m_net)
         heat_flux_diag(i, j) = 0.0_wp
         sw_thru_diag(i, j) = 0.0_wp
         m_melt_diag(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            sum_ocn = 0.0_wp
            sum_ice = 0.0_wp
            sum_heat = 0.0_wp
            sum_sw = 0.0_wp
            do cat = 1, ncat
               sum_ocn = sum_ocn + h2o_ocn_to_ice(i, j, cat)
               sum_ice = sum_ice + h2o_ice_to_ocn(i, j, cat)
               sum_heat = sum_heat + heat_to_ocn(i, j, cat)
               sum_sw = sum_sw + sw_thru(i, j, cat)
            end do
            m_net = sum_ocn - sum_ice
            m_melt_diag(i, j) = sum_ice - sum_ocn
            heat_flux_diag(i, j) = sum_heat/dt_therm - fb(i, j)
            sw_thru_diag(i, j) = sum_sw
            salt_flux_diag(i, j) = salt_flux_diag(i, j) &
                                   + m_net*(ssurf_seam(i, j) - ICE_BULK_SALINITY)/dt_therm
         end if
      end do
   end subroutine ice_thermo_driver_reduce_impl

   pure subroutine ice_thermo_driver_reduce_multicat_impl(wet_mask, part_size, h2o_ocn_to_ice, &
                                                          h2o_ice_to_ocn, heat_to_ocn, sw_thru, &
                                                          ssurf_seam, fb, fb_part_sum, &
                                                          heat_flux_diag, sw_thru_diag, &
                                                          m_melt_diag, &
                                                          salt_flux_diag, dt_therm, nghost, &
                                                          ncat, nx, ny)
      !! ncat>1 SIS2 ITD-mode reduce (module docstring). Same zero/gate/
      !! ordering contract as `ice_thermo_driver_reduce_impl` — the only
      !! change is PART-WEIGHTING the per-category sums (the column's
      !! `h2o_*`/`heat_to_ocn`/`sw_thru` outputs are per unit ICE-COVERED
      !! area in this mode, so a per-cell total needs the `part_size`
      !! weight) and subtracting `fb_part_sum*fb` instead of the bare `fb`
      !! the ncat==1 path subtracts (`fb` is a per-cell flux; `fb_part_sum`
      !! is the fraction of the cell it was actually charged against —
      !! see the module docstring and `rdb_ice_state%fb_part_sum`).
      !!
      !! PR 31: `sw_thru_diag(i,j) = Σ_cat part_size(i,j,cat)*sw_thru(i,j,cat)`
      !! (W/m^2, +down) — the area-weighted twin of the ncat==1 reduce's
      !! lumped sum, matching `heat_to_ocn`'s own part-weighting here.
      !! Same OWN-field / never-folded-into-`heat_flux_diag` contract as
      !! the ncat==1 path.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if`/serial `do cat` — never a
      !! masked `do concurrent` header.
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: h2o_ocn_to_ice(nx, ny, ncat)
      real(wp), intent(in) :: h2o_ice_to_ocn(nx, ny, ncat)
      real(wp), intent(in) :: heat_to_ocn(nx, ny, ncat)
      real(wp), intent(in) :: sw_thru(nx, ny, ncat)
      real(wp), intent(in) :: ssurf_seam(nx, ny)
      real(wp), intent(in) :: fb(nx, ny)
      real(wp), intent(in) :: fb_part_sum(nx, ny)
      real(wp), intent(inout) :: heat_flux_diag(nx, ny)
      real(wp), intent(inout) :: sw_thru_diag(nx, ny)
      real(wp), intent(inout) :: m_melt_diag(nx, ny)
      real(wp), intent(inout) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi
      real(wp) :: sum_ocn, sum_ice, sum_heat, sum_sw, m_net

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(cat, sum_ocn, sum_ice, sum_heat, sum_sw, m_net)
         heat_flux_diag(i, j) = 0.0_wp
         sw_thru_diag(i, j) = 0.0_wp
         m_melt_diag(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            sum_ocn = 0.0_wp
            sum_ice = 0.0_wp
            sum_heat = 0.0_wp
            sum_sw = 0.0_wp
            do cat = 1, ncat
               sum_ocn = sum_ocn + part_size(i, j, cat)*h2o_ocn_to_ice(i, j, cat)
               sum_ice = sum_ice + part_size(i, j, cat)*h2o_ice_to_ocn(i, j, cat)
               sum_heat = sum_heat + part_size(i, j, cat)*heat_to_ocn(i, j, cat)
               sum_sw = sum_sw + part_size(i, j, cat)*sw_thru(i, j, cat)
            end do
            m_net = sum_ocn - sum_ice
            m_melt_diag(i, j) = sum_ice - sum_ocn
            heat_flux_diag(i, j) = sum_heat/dt_therm - fb_part_sum(i, j)*fb(i, j)
            sw_thru_diag(i, j) = sum_sw
            salt_flux_diag(i, j) = salt_flux_diag(i, j) &
                                   + m_net*(ssurf_seam(i, j) - ICE_BULK_SALINITY)/dt_therm
         end if
      end do
   end subroutine ice_thermo_driver_reduce_multicat_impl

end module rdb_ice_thermo_driver
