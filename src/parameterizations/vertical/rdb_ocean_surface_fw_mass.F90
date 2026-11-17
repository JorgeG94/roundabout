!! Real fresh-water MASS at the ocean surface (`&ocean_forcing_nml freshwater = "mass"`).
!!
!! The Boussinesq column carries VOLUME.  This module spends the one assembled snapshot of the
!! net open-ocean fresh-water mass flux `sf%fw_net` (kg m-2 s-1, into the ocean; see
!! `ocean_surface_flux_assemble`) as a thickness change of the first live layer, with the
!! evaporation pass-down, and dilutes or concentrates the dissolved tracers accordingly.  It is
!! the surface twin of `ocean_cavity_mass_step` (`rdb_ocean_cavity_flux`), and the same
!! derivation applies; this docstring states only what is different.
!!
!! ## Arithmetic (per wet column, over the whole array, ghosts included)
!!
!!     dh = fw_net * dt / rho_0         (rho_0 converts kg to m, as everywhere in the model)
!!
!! * `dh > 0` (rain, snow, runoff, ice melt): `h(kt) += dh`, `kt` the first LIVE layer at or
!!   below `k_top` (`rdb_vl_first_live`).
!! * `dh < 0` (evaporation, freezing): take `min(left, f*h_k, h_k - 2*H_VANISHED)` from each live layer going DOWN
!!   from `kt`, `f = fw_evap_cfl_limit` (0.8), until the request is met.  The unmet remainder of a
!!   thin top layer is thus passed down the column ("the unused mass loss is passed down through
!!   the column", the usual `EVAP_CFL_LIMIT` rule).  Each layer also keeps at least `2*H_VANISHED` (`fw_mass_take`), so no
!!   withdrawal empties a layer and strands its dissolved content.  A column whose whole live
!!   stack cannot supply the request is CLAMPED (every live layer left at `(1-f)*h` or `2*H_VANISHED`), COUNTED, and the count makes the
!!   step fail loud, rank-uniformly (one int64 reduction), because a clamped withdrawal no longer
!!   matches the tracked mass source.  A non-finite `dh` is counted the same way.
!! * The mass step reads `fw_net`, never the raw components: bulk forcing rewrites `lprec`/`evap`
!!   at the start of every outer step while the assembly happens at the end of the thermo step, so
!!   a kernel reading the components would spend step n+1's P - E against step n's salt flux and
!!   enthalpy.  Reading the snapshot (and its seam-exchanged ghosts) is what keeps `h` in the
!!   ghost columns bit-identical to the neighbour's owned ones.
!!
!! ## Tracers
!!
!! * **Temperature is not touched.**  Every mass flux enters or leaves at the temperature of the
!!   layer it touches (rain, snow and liquid runoff at the sea-surface temperature, evaporation
!!   and the pass-down withdrawal at the layer's own `T`), so `T` keeps its exact bits and the
!!   heat CONTENT `dh_k * T_k` is booked into `ms%heat_budget_surface`.  The assembler drops
!!   `heat_content_massin/massout` from `Q_heat` under mass so the enthalpy is not counted twice;
!!   the latent heat of fusion (`q_fus`) stays in `Q_heat`.  The sea-surface-temperature reference
!!   for the enthalpy of a mass flux (`SST*cp`) is the usual default of the large-scale ocean
!!   models, and the layer-`T` rule avoids the "concentrate T" absurdity of a layer that loses
!!   80 % of its volume keeping all its heat.
!! * **Salinity, pseudo-salt, every passive tracer, ideal age** (dissolved tracers stay behind):
!!   inflow carries zero concentration and withdrawal takes no tracer, so the content `h*c` is
!!   unchanged and `c <- c*h0/hn` (`rdb_conc_add_source_mass(c, h0, hn, 0)`): rain dilutes,
!!   evaporation concentrates.  Salt content is unchanged, so nothing is booked for salt.  The
!!   pseudo-salt is given exactly salinity's treatment through the same expression, so the
!!   mirror holds bitwise.  The surface salt flux of the deposit (`Q_salt_tracer`) is applied
!!   BEFORE this step, at the old thickness, by `ocean_surface_flux_apply_tracers`.
!! * **Boundary-layer buoyancy.**  `Q_salt` keeps the virtual term `salt_fw = -S*F` so KPP/EPBL
!!   see the fresh-water buoyancy forcing; only the tracer deposit is switched to
!!   `Q_salt_tracer = Q_salt - salt_fw`.  The virtual term is thereby removed from the tracers
!!   exactly once, structurally, in the assembler (Griffies et al. 2001; Huang 1993).
!!
!! ## Placement and weights
!!
!! Called immediately after every `ocean_surface_flux_apply_tracers`: `run_stage` (weight 0.5 per
!! SSP-RK2 stage), the split stage (`chain_weight`: 0.5 for ssp_rk2, 0 for the pred_corr predictor
!! (which never runs it: thermo is off there), 1 for the corrector) and `apply_deferred_thermo`
!! (weight 1; `diabatic_first = .false.`, `dt_therm_ratio > 1`).  The barotropic mode sees the
!! change through `derive_bt_from_layers` at the next stage top; the ALE remap takes
!! `h_ref = sum(h) - bt_eta`, so the engine re-derives `bt_eta = sum(h) - bt_H_ref` right before
!! the remap (`ocean_dyn_step_split`), after the stage average, or the remap would build its
!! target on `H + dh`.
!!
!! ## Budgets
!!
!! The applied `dh` (never the request) times `areaT * RHO_WATER * weight` is summed over this
!! tile's OWNED cells as EFP bins into `ms%mass_src_efp` (no collective; the console sums ranks),
!! so `(M - M0) + mass_out - mass_src` closes to round-off and is the same on every rank count.
!!
!! ## What this does not do
!! * No sea ice: under `scheme = "cice"` or the old ice scheme the mode is refused at configure
!!   (the ice hand-off has no mass slot yet).
!! * No dynamic wet/dry, no `eulerian_z` / `lagrangian` coordinate (no remap to give the volume
!!   back to the target), no cavity: refused by name in `validate_config`.
module rdb_ocean_surface_fw_mass
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS, RHO_WATER, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_efp, only: efp_t, EFP_DIGITS, EFP_MAX_SUMMANDS, efp_carry
   use rdb_ocean_console_stats, only: efp_decompose_impl
   use rdb_grid_reduce, only: grid_allreduce_sum_i8
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   implicit none
   private

   public :: ocean_surface_fw_mass_step
   public :: fw_mass_column_impl
   public :: fw_mass_dilute_impl
   public :: fw_mass_volume_efp
   public :: fw_mass_thin_is_fatal
   public :: fw_mass_rederive_bt_eta

contains

   pure function fw_mass_thin_is_fatal(n_thin) result(is_fatal)
      !! The verdict on the global count of columns that could not give up the requested
      !! thickness (or fed the step a non-finite flux): any is fatal.
      integer(int64), intent(in) :: n_thin
      logical :: is_fatal
      is_fatal = (n_thin > 0_int64)
   end function fw_mass_thin_is_fatal

   subroutine ocean_surface_fw_mass_step(grid, sf, ms, areaT, dt, weight, active)
      !! The real fresh-water mass update, once per thermo application of the surface fluxes,
      !! immediately after `ocean_surface_flux_apply_tracers` (see the module docstring).
      !!
      !! NON-`pure`: it logs, fails loud, writes `ms%mass_src_efp` and the counters of `ms`, and
      !! issues ONE integer reduction.  Its kernels are `pure`.  `mem:separate`: every array is
      !! mapped by its owning slot; the tracer-registry derefs stay on the host before each
      !! kernel (outer-shim rule).
      use pic_logger, only: global_logger
      use pic_strings, only: to_string
      use rdb_error_ring, only: fail
      use rdb_ocean_status, only: OCEAN_STATUS_ERR_SETUP
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in) :: sf
         !! Reads `fw_net`, `fw_evap_cfl_limit`, `rho0` (read-only: the dynamics hand it down so).
      type(multilayer_state_t), intent(inout) :: ms
         !! Writes `h_layer`, the tracer loads, `heat_budget_surface`, `mass_src_efp`,
         !! `fw_dh`, `fw_passdown` and the counters.
      real(wp), intent(in) :: areaT(grid%nx_total, grid%ny_total)
         !! Cell area `metrics%areaT` (explicit-shape from `grid`, so no repack copy can
         !! sit between the device-resident array and the reduction's `present` clause).
      real(wp), intent(in) :: dt
         !! Thermo timestep (s): the one `ocean_surface_flux_apply_tracers` was given.
      real(wp), intent(in) :: weight
         !! Per-call weight of the `mass_src` accumulator (`ocean_accumulate_mass_out`'s).
      logical, intent(in), optional :: active
         !! Thermo-cadence gate; present-and-false => early return.

      integer :: nx, ny, nz, idx_t, idx_s, it, n_thin, n_pd
      real(wp) :: dt_over_rho0
      type(efp_t) :: tot
      integer(int64) :: loc, glob
      integer(int64), parameter :: PACK = 4294967296_int64

      if (.not. sf%fw_mass) return
      if (.not. sf%use_components) return
      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      if (idx_t <= 0 .or. idx_s <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      dt_over_rho0 = dt/sf%rho0

      ! Every dissolved tracer first (they need the pre-source thickness): S, the pseudo-salt
      ! mirror, passive tracers, ideal age.  Temperature is not diluted (see the docstring).
      do it = 1, size(ms%tracers)
         if (it == idx_t) cycle
         if (.not. allocated(ms%tracers(it)%c)) cycle
         call fw_mass_dilute_impl(nx, ny, nz, dt_over_rho0, sf%fw_evap_cfl_limit, sf%fw_net, &
                                  ms%wet_mask, ms%k_top, ms%h_layer, ms%tracers(it)%c)
      end do
      call fw_mass_column_impl(nx, ny, nz, grid%nghost, dt_over_rho0, sf%fw_evap_cfl_limit, &
                               sf%fw_net, ms%wet_mask, ms%k_top, ms%h_layer, &
                               ms%tracers(idx_t)%c, ms%heat_budget_surface, &
                               ms%fw_dh, ms%fw_passdown, n_thin, n_pd)

      ! The tracked source: APPLIED volume x weight x RHO_WATER, EFP over owned cells.
      tot = fw_mass_volume_efp(nx, ny, grid%nghost, weight*RHO_WATER, ms%fw_dh, areaT)
      ms%mass_src_efp = ms%mass_src_efp + tot%v
      ms%mass_src_efp_poison = ms%mass_src_efp_poison + tot%poison
      call efp_carry(ms%mass_src_efp)

      ! Rank-uniform verdict, one packed reduction (each count is far below 2**32).
      loc = int(n_thin, int64) + int(n_pd, int64)*PACK
      call grid_allreduce_sum_i8(grid, loc, glob)
      ms%fw_n_thin_step = modulo(glob, PACK)
      ms%fw_n_passdown_step = glob/PACK
      ms%fw_n_passdown_total = ms%fw_n_passdown_total + ms%fw_n_passdown_step
      if (fw_mass_thin_is_fatal(ms%fw_n_thin_step)) then
         call global_logger%error("surface fresh-water mass: "//to_string(int(ms%fw_n_thin_step))// &
                                  " column(s) could not give up the requested thickness")
         call fail("&ocean_forcing_nml freshwater='mass': "//to_string(int(ms%fw_n_thin_step))// &
                   " column(s) cannot supply the evaporative withdrawal from their whole live "// &
                   "stack at fw_evap_cfl_limit (or carry a non-finite fresh-water flux).  The "// &
                   "withdrawal was clamped so the state stays finite, but a clamped withdrawal "// &
                   "no longer matches the tracked mass source.  Use thicker layers, a shorter "// &
                   "dt, or a smaller evaporation.", code=OCEAN_STATUS_ERR_SETUP)
      end if
   end subroutine ocean_surface_fw_mass_step

   pure subroutine fw_mass_column_impl(nx, ny, nz, nghost, dt_over_rho0, cfl, fw_net, wet_mask, &
                                       k_top, h_layer, c_T, heat_budget, fw_dh, fw_passdown, &
                                       n_thin, n_pd)
      !! The thickness update, the heat-content booking and the counts (see the module
      !! docstring).  The withdrawal profile is recomputed identically by `fw_mass_dilute_impl`
      !! from the pre-update thickness: keep the two loops textually in step.
      integer, intent(in) :: nx, ny, nz, nghost
      real(wp), intent(in) :: dt_over_rho0, cfl
      real(wp), intent(in) :: fw_net(nx, ny), wet_mask(nx, ny)
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(inout) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: c_T(nx, ny, nz)
         !! Temperature (concentration): unchanged, read for the heat content.
      real(wp), intent(inout) :: heat_budget(nx, ny, nz)
         !! `ms%heat_budget_surface` (content `h*T` units).
      real(wp), intent(inout) :: fw_dh(nx, ny), fw_passdown(nx, ny)
      integer, intent(out) :: n_thin
         !! Owned columns that could not be satisfied (or non-finite flux).
      integer, intent(out) :: n_pd
         !! Owned columns that drew from a layer below the first live one.
      integer :: i, j, k, kt, c_thin, c_pd, npd_col
      real(wp) :: dh, left, take, h0, dh_app
      logical :: owned

      c_thin = 0
      c_pd = 0
      do concurrent(j=1:ny, i=1:nx) local(k, kt, dh, left, take, h0, dh_app, npd_col, owned) &
         reduce(+:c_thin, c_pd)
         fw_dh(i, j) = 0.0_wp
         fw_passdown(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            owned = (i > nghost .and. i <= nx - nghost .and. j > nghost .and. j <= ny - nghost)
            dh = fw_net(i, j)*dt_over_rho0
            kt = rdb_vl_first_live(nx, ny, nz, h_layer, i, j, k_top(i, j))
            if (.not. ieee_is_finite(dh)) then
               if (owned) c_thin = c_thin + 1
            else if (dh > 0.0_wp .and. h_layer(i, j, kt) > H_VANISHED) then
               h_layer(i, j, kt) = h_layer(i, j, kt) + dh
               heat_budget(i, j, kt) = heat_budget(i, j, kt) + dh*c_T(i, j, kt)
               fw_dh(i, j) = dh
            else if (dh < 0.0_wp .and. h_layer(i, j, kt) > H_VANISHED) then
               left = -dh
               dh_app = 0.0_wp
               npd_col = 0
               do k = kt, 1, -1
                  h0 = h_layer(i, j, k)
                  if (h0 > H_VANISHED) then
                     take = fw_mass_take(left, cfl, h0)
                     if (take > 0.0_wp) then
                        h_layer(i, j, k) = h0 - take
                        heat_budget(i, j, k) = heat_budget(i, j, k) - take*c_T(i, j, k)
                        dh_app = dh_app - take
                        left = left - take
                        if (k < kt) npd_col = npd_col + 1
                     end if
                     if (left <= 0.0_wp) exit
                  end if
               end do
               if (left > 0.0_wp .and. owned) c_thin = c_thin + 1
               if (npd_col > 0 .and. owned) c_pd = c_pd + 1
               fw_dh(i, j) = dh_app
               fw_passdown(i, j) = real(npd_col, wp)
            end if
         end if
      end do
      n_thin = c_thin
      n_pd = c_pd
   end subroutine fw_mass_column_impl

   pure function fw_mass_take(left, cfl, h0) result(take)
      !! The thickness one withdrawal takes from a live layer of thickness `h0`:
      !! `min(left, cfl*h0, h0 - 2*H_VANISHED)`, never negative.  The last bound keeps the layer
      !! LIVE after the withdrawal, so `c <- c*h0/hn` always runs with `hn > H_VANISHED` and the
      !! dissolved content `h*c` stays in the column: evaporation removes fresh water only, and a
      !! layer emptied to `hn <= H_DIV_EPS` would keep its `c` and lose `h0*c` (the reason
      !! `fw_evap_cfl_limit = 1` is refused).  A layer thinner than `2*H_VANISHED` gives nothing and
      !! the request passes down.  `fw_mass_column_impl` and `fw_mass_dilute_impl` both call it.
      !$acc routine seq
      real(wp), intent(in) :: left, cfl, h0
      real(wp) :: take
      take = max(0.0_wp, min(left, cfl*h0, h0 - 2.0_wp*H_VANISHED))
   end function fw_mass_take

   pure subroutine fw_mass_dilute_impl(nx, ny, nz, dt_over_rho0, cfl, fw_net, wet_mask, k_top, &
                                       h_layer, c)
      !! Dilution / concentration of ONE dissolved tracer by the volume change of
      !! `fw_mass_column_impl`: `c <- c*h0/hn` on every layer that gains or loses water, the
      !! content `h*c` unchanged.  Reads the thickness BEFORE the update (call first).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt_over_rho0, cfl
      real(wp), intent(in) :: fw_net(nx, ny), wet_mask(nx, ny)
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
         !! Thickness before the update.
      real(wp), intent(inout) :: c(nx, ny, nz)
         !! Tracer concentration.
      integer :: i, j, k, kt
      real(wp) :: dh, left, take, h0

      do concurrent(j=1:ny, i=1:nx) local(k, kt, dh, left, take, h0)
         if (wet_mask(i, j) > 0.5_wp) then
            dh = fw_net(i, j)*dt_over_rho0
            kt = rdb_vl_first_live(nx, ny, nz, h_layer, i, j, k_top(i, j))
            if (.not. ieee_is_finite(dh)) then
               continue
            else if (dh > 0.0_wp .and. h_layer(i, j, kt) > H_VANISHED) then
               h0 = h_layer(i, j, kt)
               c(i, j, kt) = rdb_conc_add_source_mass(c(i, j, kt), h0, h0 + dh, 0.0_wp)
            else if (dh < 0.0_wp .and. h_layer(i, j, kt) > H_VANISHED) then
               left = -dh
               do k = kt, 1, -1
                  h0 = h_layer(i, j, k)
                  if (h0 > H_VANISHED) then
                     take = fw_mass_take(left, cfl, h0)
                     if (take > 0.0_wp) then
                        c(i, j, k) = rdb_conc_add_source_mass(c(i, j, k), h0, h0 - take, 0.0_wp)
                        left = left - take
                     end if
                     if (left <= 0.0_wp) exit
                  end if
               end do
            end if
         end if
      end do
   end subroutine fw_mass_dilute_impl

   function fw_mass_volume_efp(nx, ny, nghost, scale, dh, areaT) result(tot)
      !! RANK-LOCAL tracked mass source: `sum(scale * dh * areaT)` over the OWNED cells
      !! (`nghost+1 .. n-nghost`) as an EFP reproducing sum (layout-independent once the ranks'
      !! values are added).  Device: one small `!$acc parallel loop` reduction with the console's
      !! proven seven-accumulator footprint.  A NaN/Inf summand sets `poison`.
      integer, intent(in) :: nx, ny, nghost
      real(wp), intent(in) :: scale
         !! `weight*RHO_WATER`.
      real(wp), intent(in) :: dh(nx, ny)
         !! Applied thickness change per column (m).
      real(wp), intent(in) :: areaT(nx, ny)
      type(efp_t) :: tot
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      integer(int64) :: e1, e2, e3, e4, e5, e6, epoison
      integer(int64) :: d1, d2, d3, d4, d5, d6, dpoison
      integer(int64) :: bins(EFP_DIGITS)
      real(real64) :: val

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      if (int(i_hi - i_lo + 1, int64)*int(j_hi - j_lo + 1, int64) > EFP_MAX_SUMMANDS) then
         error stop "fw_mass_volume_efp: cell count exceeds EFP_MAX_SUMMANDS"
      end if
      e1 = 0_int64
      e2 = 0_int64
      e3 = 0_int64
      e4 = 0_int64
      e5 = 0_int64
      e6 = 0_int64
      epoison = 0_int64
      !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6,epoison) &
      !$acc&         private(val, d1, d2, d3, d4, d5, d6, dpoison) present(dh, areaT)
      do j = j_lo, j_hi
         do i = i_lo, i_hi
            val = real(scale, real64)*real(dh(i, j), real64)*real(areaT(i, j), real64)
            call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6, dpoison)
            e1 = e1 + d1
            e2 = e2 + d2
            e3 = e3 + d3
            e4 = e4 + d4
            e5 = e5 + d5
            e6 = e6 + d6
            epoison = epoison + dpoison
         end do
      end do
      bins = [e1, e2, e3, e4, e5, e6]
      call efp_carry(bins)
      tot%v = bins
      tot%poison = epoison
   end function fw_mass_volume_efp

   pure subroutine fw_mass_rederive_bt_eta(nx, ny, nz, h_layer, bt_H_ref, bt_eta)
      !! `bt_eta = sum_k h - bt_H_ref` over the full array: the barotropic free surface the ALE
      !! remap reads (`h_ref = sum_k h - bt_eta`) once a real volume source has moved the column
      !! sum.  Without it the z*/z_fixed target is built on `H + dh` instead of `H` (the column
      !! sum is still conserved, so no conservation test would notice).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz), bt_H_ref(nx, ny)
      real(wp), intent(inout) :: bt_eta(nx, ny)
      integer :: i, j, k
      real(wp) :: total_h

      do concurrent(j=1:ny, i=1:nx) local(k, total_h)
         total_h = 0.0_wp
         do k = 1, nz
            total_h = total_h + h_layer(i, j, k)
         end do
         bt_eta(i, j) = total_h - bt_H_ref(i, j)
      end do
   end subroutine fw_mass_rederive_bt_eta

#include "rdb_vanished_layer.inc"
#include "rdb_tracer_conc.inc"

end module rdb_ocean_surface_fw_mass
