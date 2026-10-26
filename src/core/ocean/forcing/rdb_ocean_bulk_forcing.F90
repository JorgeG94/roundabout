!! OCEAN consumer of the JRA55-do-style prescribed atmospheric state
!! (`rdb_atm_state`) -- OM3 forcing wave-3, C2c-3a.
module rdb_ocean_bulk_forcing
   !! Turns the nine-field `atm_state_t` into the ocean's own surface
   !! forcing slots (`ocean_surface_stress_t%tau_x/tau_y`,
   !! `ocean_surface_flux_t`'s component set): rotate the prescribed
   !! wind onto the model grid axes, form the relative wind against the
   !! ocean's own surface current, run the Large & Yeager bulk-flux
   !! column kernel (`rdb_ocean_bulk_flux`) cell by cell, write the
   !! shortwave / longwave / precipitation / pressure-loading
   !! components, apply the sea-ice open-water weighting when a
   !! sea-ice component is configured, and finish by turning the
   !! cell-centred stress into the C-grid face pair.
   !!
   !! **No ice seam here.**  The atmosphere -> sea-ice forcing
   !! dispatcher (`ice_atm_forcing_from_atm`) is a separate chunk
   !! (C2c-3b); this module only READS the sea-ice state's already-live
   !! `part_size` category array to scale the ocean's own open-water
   !! fraction -- it never writes into `ocean_sea_ice_t`.
   !!
   !! **Why the bespoke per-cell kernel instead of calling
   !! `bulk_flux_driver_2d` over the whole `(nx_total, ny_total)`
   !! tile.**  `atm_state_t`'s own docstring is explicit: every field
   !! array is ghosted but only the PHYSICAL interior is ever written
   !! -- "no consumer in this chunk reads a ghost cell" -- so a ghost
   !! cell's `tas`/`huss`/`psl` sits at its `source=0.0_wp` allocation
   !! value forever.  Calling the bulk-flux column kernel there would
   !! divide by `t_air = 0` (`rho_air = slp/(R*t_air*(...))`) and most
   !! namelists this provider is gated on (`grid_config` must be
   !! spherical/supergrid/tripolar, and a tripolar cap REQUIRES
   !! periodic east/west) have genuinely WET ghost cells, so this is
   !! not a hypothetical edge case.  `bulk_forcing_core_impl` below
   !! loops the physical interior only (the `ii = nghost + i` offset
   !! idiom `rdb_atm_state`'s own `atm_state_scatter_impl` uses) and
   !! calls `bulk_flux_column` -- the same per-column kernel
   !! `bulk_flux_driver_2d` itself calls -- directly, so every physical
   !! cell's `q_sens`/`q_lat`/`evap`/stress is bit-for-bit what a direct
   !! `bulk_flux_column` call produces (what the C2c tests assert).
   !! Ghost cells of every output are intentionally left UNTOUCHED by
   !! this kernel, matching `rdb_ocean_data_forcing`'s "ghosts are
   !! EXCHANGED, never extrapolated" contract; the stress pair's
   !! ghosts are filled afterwards, see below.
   !!
   !! **Stress: cell-centred, then a scalar halo copy, then the face
   !! average, then the REAL vector seam refresh.**  `bulk_flux_column`
   !! naturally produces a CELL-CENTRED stress (same point as `uas`/
   !! `vas`).  Averaging that onto the C-grid `tau_x`/`tau_y` FACES
   !! needs the neighbour cell on both sides of every face, including
   !! the faces at the tile edge -- which needs this rank's cell-centred
   !! ghost ring to hold the NEIGHBOUR's real physical value, not a
   !! locally-recomputed (zero-atm-input) one.  `ocean_halo_centre`
   !! (the 2D specific) supplies exactly that: an MPI exchange (or local periodic wrap)
   !! of the plain SCALAR field, no vector sign handling -- which is
   !! why the tripolar fold row comes out of it with the wrong sign
   !! (a scalar copy, not a vector one) and is deliberately not relied
   !! on: `ocean_seam_refresh_surface_stress`, called last, re-does the
   !! face pair's own proper vector exchange + periodic wrap + fold and
   !! is what leaves `tau_x`/`tau_y` (and `stress_mag`) correct
   !! everywhere a kernel reads them.
   !!
   !! **Open-water weighting.**  `ocean_sea_ice_t%part_size(:,:,0)` IS
   !! the open-water fraction already -- the category partition sums
   !! to 1 by construction (`part_size(:,:,0) = 1 - conc` at the ITD
   !! seed, `rdb_ice_init`), so `1 - Sum_{c=1}^{ncat} part_size(c)`
   !! and `part_size(0)` are the same number; this module reads the
   !! field directly rather than re-summing it. `this%open_water` is
   !! seeded to 1 at configure and is the array `bulk_forcing_core_impl`
   !! always multiplies by; it is refreshed from `ice%part_size(:,:,0)`
   !! only when `ice%enable`, so a run with no sea-ice component never
   !! touches an `ocean_sea_ice_t` array (CLAUDE.md's "a host-gated call
   !! that hands a state array to an external subroutine still costs" --
   !! here the guard is at the CALL SITE, so the array is never even
   !! named when ice is off).
   !!
   !! **Component ownership.**  Every field this module writes is a
   !! PR-12 component (`q_sens`, `q_lat`, `q_sw`, `q_lw`, `evap`,
   !! `lprec`, `fprec`, `p_surf_atm`) -- never `Q_heat`/`Q_salt`
   !! directly (`ocean_surface_flux_assemble` rebuilds those every
   !! thermo step).  `lprec`/`fprec` are mass-flux components carrying
   !! enthalpy, so this module also fills their `heat_content_<flux>`
   !! companions -- isothermal at the ocean's own SST, the convention
   !! `rdb_ocean_surface_flux`'s own component contract and
   !! `tests/test_ocean_surface_forcing_type.F90` already use ("rain at
   !! T_s"). No salt component is written (v1 bulk fluxes carry no
   !! virtual-salt term), so `has_salt` is left untouched.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t, ocean_atm_forcing_config_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_atm_state, only: atm_state_t
   use rdb_teos10_convert, only: pt_from_ct
   use rdb_ocean_bulk_flux, only: bulk_flux_config_t, bulk_flux_column, &
                                  bulk_flux_longwave_up, BULK_T0_KELVIN
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress
   use rdb_ocean_halo, only: ocean_halo_centre
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: bulk_forcing_t
   public :: bulk_forcing_configure
   public :: bulk_forcing_apply

   type :: bulk_forcing_t
      !! Persistent per-step scratch for the ocean bulk-flux consumer.
      !! `taux_c`/`tauy_c` MUST be module-/state-level (never a local
      !! automatic array inside `bulk_forcing_apply`): they are written
      !! by `bulk_forcing_core_impl`, then mutated IN PLACE by the halo
      !! exchange, then read by the face-average kernel -- three
      !! separate calls sharing one buffer every step, and a per-step
      !! local-allocate of a `(nx_total,ny_total)` array stalls the
      !! device on entry/exit (CLAUDE.md "Persistent kernel
      !! workspaces").
      logical :: enable = .false.
         !! Mirrors `&ocean_bulk_flux_nml enable` (`cfg%ocean%bulk_flux
         !! %enable`). `.false.` (default): every procedure below is an
         !! immediate no-op -- bit-identical to a pre-C2c-3a build.
      logical :: on_device = .false.
      logical :: cons_abs = .false.
         !! Latched `&ocean_eos_nml ts_convention == "cons_abs"`: the model
         !! T/S are CT/SA, so the SST fed to the bulk formulae is
         !! `pt_from_ct(SA, CT)` of the top layer (IOC/SCOR/IAPSO 2010;
         !! McDougall & Barker 2011).  `.false.` (default) is the original
         !! expression, bit-identical.
      integer :: nx_total = 0, ny_total = 0
      type(ocean_atm_forcing_config_t) :: cfg_atm
         !! Latched copy of `cfg%ocean%atm_forcing` -- `wind_relative`/
         !! `ocean_albedo` -- so `bulk_forcing_apply` needs no `cfg`
         !! argument (mirrors `ocean_engine_t%diag_enabled`'s "latch at
         !! setup" idiom, scoped to this slot instead of the engine).
      type(bulk_flux_config_t) :: cfg_bulk
         !! Latched copy of `cfg%ocean%bulk_flux` -- `n_iter`/`z_wind`/
         !! `z_ta`/`emissivity`.
      real(wp), allocatable :: taux_c(:, :), tauy_c(:, :)
         !! Cell-centred stress (T point, same shape as `uas`/`vas`),
         !! `(nx_total, ny_total)`.
      real(wp), allocatable :: open_water(:, :)
         !! Open-water fraction, `(nx_total, ny_total)`. Seeded to 1 at
         !! configure; refreshed from `ice%part_size(:,:,0)` once per
         !! step only when `ice%enable` -- otherwise stays 1 forever,
         !! so a run without sea ice never reads an ice array.
   contains
      procedure, non_overridable :: enter_data => bulk_forcing_enter_data
      procedure, non_overridable :: exit_data => bulk_forcing_exit_data
      procedure, non_overridable :: bytes => bulk_forcing_bytes
      procedure, non_overridable :: destroy => bulk_forcing_destroy
   end type bulk_forcing_t

contains

   ! ======================================================================
   ! Lifecycle
   ! ======================================================================

   subroutine bulk_forcing_enter_data(this)
      class(bulk_forcing_t), intent(inout) :: this
      if (.not. this%enable) return
      !$acc enter data copyin(this%taux_c, this%tauy_c, this%open_water)
      this%on_device = .true.
   end subroutine bulk_forcing_enter_data

   subroutine bulk_forcing_exit_data(this)
      class(bulk_forcing_t), intent(inout) :: this
      if (.not. this%enable) return
      !$acc exit data delete(this%taux_c, this%tauy_c, this%open_water)
      this%on_device = .false.
   end subroutine bulk_forcing_exit_data

   pure function bulk_forcing_bytes(this) result(nbytes)
      class(bulk_forcing_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (.not. this%enable) return
      nbytes = arr_bytes(this%taux_c) + arr_bytes(this%tauy_c) + arr_bytes(this%open_water)
   end function bulk_forcing_bytes

   subroutine bulk_forcing_destroy(this)
      class(bulk_forcing_t), intent(inout) :: this
      if (this%on_device) call this%exit_data()
      if (allocated(this%taux_c)) deallocate (this%taux_c)
      if (allocated(this%tauy_c)) deallocate (this%tauy_c)
      if (allocated(this%open_water)) deallocate (this%open_water)
      this%enable = .false.
      this%nx_total = 0
      this%ny_total = 0
   end subroutine bulk_forcing_destroy

   ! ======================================================================
   ! Configure (once, before ocean_state_enter_data)
   ! ======================================================================

   subroutine bulk_forcing_configure(cfg, grid, sf, this)
      !! A no-op when `&ocean_bulk_flux_nml enable = .false.` (the
      !! default) -- nothing allocates and `bulk_forcing_apply` stays a
      !! no-op: bit-identical. `validate_config` ties `enable` to
      !! `&ocean_atm_forcing_nml enable` (each requires the other), so
      !! by the time this runs with `enable = .true.` the atm provider
      !! is configured too.
      !!
      !! Sets the host-side component latches
      !! (`has_heat`/`has_mass_flux`/`has_q_sw`) ONCE here rather than
      !! every step -- same contract `rdb_ocean_data_forcing`'s own
      !! configure uses for its tags, and required for
      !! `ocean_surface_flux_apply_tracers` to ever read the `Q_heat`
      !! the assembler derives from the components this module fills
      !! (its own guard is `if (.not. has_heat .and. .not. has_salt)
      !! return`).
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(bulk_forcing_t), intent(out) :: this

      this%enable = cfg%ocean%bulk_flux%enable
      if (.not. this%enable) return

      this%cfg_atm = cfg%ocean%atm_forcing
      this%cfg_bulk = cfg%ocean%bulk_flux
      this%cons_abs = trim(adjustl(cfg%ocean%eos%ts_convention)) == "cons_abs"
      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      allocate (this%taux_c(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%tauy_c(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%open_water(this%nx_total, this%ny_total), source=1.0_wp)

      sf%has_heat = .true.
      sf%has_mass_flux = .true.
      sf%has_q_sw = .true.
   end subroutine bulk_forcing_configure

   ! ======================================================================
   ! Per-step apply (immediately after atm_state_update)
   ! ======================================================================

   subroutine bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, this, ss, sf)
      !! No-op unless both `this%enable` (the bulk-flux consumer) and
      !! `atm%enable` (the atmospheric-state provider) hold -- the
      !! second check is paranoia, since `validate_config` already
      !! requires one iff the other. Reads `this%cfg_atm`/
      !! `this%cfg_bulk` (latched at `bulk_forcing_configure`) rather
      !! than taking a `cfg` argument -- the engine's per-step call
      !! site has no `cfg` in scope (mirrors `ocean_engine_t
      !! %diag_enabled`'s "latch at setup" idiom).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_bc_state_t), intent(in) :: bc
      type(atm_state_t), intent(in) :: atm
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_sea_ice_t), intent(in) :: ice
      type(bulk_forcing_t), intent(inout) :: this
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf

      integer :: nx, ny, ng, nxp, nyp, nz, idx_T, idx_Sarg
      logical :: cons_abs

      if (.not. this%enable) return
      if (.not. atm%enable) return
      if (.not. allocated(ms%tracers)) return
      idx_T = ms%idx_temperature
      if (idx_T <= 0) return

      ! TEOS-10: export theta = pt_from_ct(SA, CT) to the bulk formulae.
      ! The salinity tracer is the SA source; without one there is nothing
      ! to convert with, so the model T is passed as is.
      cons_abs = this%cons_abs .and. ms%idx_salinity > 0
      idx_Sarg = idx_T
      if (cons_abs) idx_Sarg = ms%idx_salinity

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost
      nxp = grid%nx_phys
      nyp = grid%ny_phys
      nz = ms%nz_ml

      ! Open-water fraction: only touch the sea-ice state when it is
      ! actually enabled (module docstring) -- `this%open_water` stays
      ! at its configure-time 1.0 seed otherwise.
      if (ice%enable) then
         call bulk_forcing_open_water_impl(nx, ny, ng, nxp, nyp, ice%ncat, &
                                           ice%part_size, this%open_water)
      end if

      call bulk_forcing_core_impl(nx, ny, ng, nxp, nyp, nz, this%cfg_bulk%n_iter, &
                                  this%cfg_atm%wind_relative, this%cfg_atm%ocean_albedo, &
                                  this%cfg_bulk%z_wind, this%cfg_bulk%z_ta, this%cfg_bulk%emissivity, &
                                  metrics%wet_T, metrics%angle_dx, &
                                  atm%uas, atm%vas, atm%tas, atm%huss, atm%psl, &
                                  atm%rlds, atm%rsds, atm%prra, atm%prsn, &
                                  ms%u_face_x_layer, ms%v_face_y_layer, &
                                  ms%tracers(idx_T)%hTr, ms%tracers(idx_Sarg)%hTr, ms%h_layer, &
                                  cons_abs, this%open_water, &
                                  this%taux_c, this%tauy_c, &
                                  sf%q_sens, sf%q_lat, sf%evap, sf%q_sw, sf%q_lw, &
                                  sf%lprec, sf%fprec, sf%p_surf_atm, &
                                  sf%heat_content_lprec, sf%heat_content_fprec, sf%cp)

      ! Scalar halo of the cell-centred stress -- see the module
      ! docstring for why this (not the vector seam refresh) runs
      ! first: the face average below needs a real neighbour-cell
      ! value in every ghost cell, including at tile/periodic seams.
      call ocean_halo_centre(this%taux_c)
      call ocean_halo_centre(this%tauy_c)

      call bulk_forcing_stress_face_impl(nx, ny, this%taux_c, this%tauy_c, ss%tau_x, ss%tau_y)

      ! The REAL vector exchange + periodic wrap + tripolar fold, and
      ! the `stress_mag` refresh KPP/EPBL read for u_* -- overwrites
      ! whatever sign the scalar halo above left at a fold row.
      call ocean_seam_refresh_surface_stress(ss, grid, bc)
   end subroutine bulk_forcing_apply

   pure subroutine bulk_forcing_open_water_impl(nx, ny, ng, nxp, nyp, ncat, part_size, open_water)
      !! `open_water(i,j) = part_size(i,j,0)` over the physical
      !! interior -- the open-water fraction IS category 0 by
      !! construction (module docstring); no re-summation needed.
      !! Explicit-shape, integer dims first (decl-order).
      integer, intent(in) :: nx, ny, ng, nxp, nyp, ncat
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: open_water(nx, ny)
      integer :: i, j, ii, jj
      do concurrent(j=1:nyp, i=1:nxp) local(ii, jj)
         ii = ng + i
         jj = ng + j
         open_water(ii, jj) = part_size(ii, jj, 0)
      end do
   end subroutine bulk_forcing_open_water_impl

   pure subroutine bulk_forcing_core_impl(nx, ny, ng, nxp, nyp, nz, n_iter, &
                                          wind_relative, ocean_albedo, z_wind, z_ta, emissivity, &
                                          wet_T, angle_dx, uas, vas, tas, huss, psl, &
                                          rlds, rsds, prra, prsn, &
                                          u_face_x_layer, v_face_y_layer, hTr_T, hTr_S, h_layer, &
                                          cons_abs, open_water, &
                                          taux_c, tauy_c, q_sens, q_lat, evap, q_sw, q_lw, &
                                          lprec, fprec, p_surf_atm, &
                                          heat_content_lprec, heat_content_fprec, cp)
      !! One `bulk_flux_column` call per PHYSICAL T-cell (see the
      !! module docstring for why ghost cells are never touched here).
      !! Land / dry cells (`wet_T <= 0.5`) write explicit zero to every
      !! output, matching `bulk_flux_driver_2d`'s own NaN-safe `if`
      !! convention (a multiply by a zero mask would not scrub a NaN
      !! land input). Explicit-shape throughout, integer dims first.
      integer, intent(in) :: nx, ny, ng, nxp, nyp, nz, n_iter
      logical, intent(in) :: wind_relative
      real(wp), intent(in) :: ocean_albedo, z_wind, z_ta, emissivity
      real(wp), intent(in) :: wet_T(nx, ny), angle_dx(nx, ny)
      real(wp), intent(in) :: uas(nx, ny), vas(nx, ny), tas(nx, ny), huss(nx, ny), psl(nx, ny)
      real(wp), intent(in) :: rlds(nx, ny), rsds(nx, ny), prra(nx, ny), prsn(nx, ny)
      real(wp), intent(in) :: u_face_x_layer(nx + 1, ny, nz), v_face_y_layer(nx, ny + 1, nz)
      real(wp), intent(in) :: hTr_T(nx, ny, nz), hTr_S(nx, ny, nz), h_layer(nx, ny, nz)
      real(wp), intent(in) :: open_water(nx, ny)
      real(wp), intent(inout) :: taux_c(nx, ny), tauy_c(nx, ny)
      real(wp), intent(inout) :: q_sens(nx, ny), q_lat(nx, ny), evap(nx, ny)
      real(wp), intent(inout) :: q_sw(nx, ny), q_lw(nx, ny)
      real(wp), intent(inout) :: lprec(nx, ny), fprec(nx, ny), p_surf_atm(nx, ny)
      real(wp), intent(inout) :: heat_content_lprec(nx, ny), heat_content_fprec(nx, ny)
      logical, intent(in) :: cons_abs
         !! `.true.` under `&ocean_eos_nml ts_convention = "cons_abs"`
         !! (with a salinity tracer): the SST handed to the bulk formulae
         !! and to the longwave emission is `pt_from_ct(SA, CT)` of the top
         !! layer, SA = `hTr_S/h`.  `.false.` takes the model T as is
         !! (bit-identical default).
      real(wp), intent(in) :: cp
         !! Specific heat (J/kg/K) of the heat content of precipitation:
         !! `sf%cp` (`&ocean_eos_nml c_p`, default 3992).

      integer :: i, j, ii, jj
      real(wp) :: ow, cosA, sinA, u_g, v_g, u_c, v_c, sst, sst_f, sst_k, lw_up
      real(wp) :: taux, tauy, sens, lat, evp, cd, ch, ce, ustar

      do concurrent(j=1:nyp, i=1:nxp) &
         local(ii, jj, ow, cosA, sinA, u_g, v_g, u_c, v_c, sst, sst_f, sst_k, lw_up, &
               taux, tauy, sens, lat, evp, cd, ch, ce, ustar)
         ii = ng + i
         jj = ng + j
         if (wet_T(ii, jj) > 0.5_wp) then
            ow = open_water(ii, jj)
            cosA = cos(angle_dx(ii, jj))
            sinA = sin(angle_dx(ii, jj))
            u_g = cosA*uas(ii, jj) + sinA*vas(ii, jj)
            v_g = -sinA*uas(ii, jj) + cosA*vas(ii, jj)
            if (wind_relative) then
               u_c = 0.5_wp*(u_face_x_layer(ii, jj, nz) + u_face_x_layer(ii + 1, jj, nz))
               v_c = 0.5_wp*(v_face_y_layer(ii, jj, nz) + v_face_y_layer(ii, jj + 1, nz))
            else
               u_c = 0.0_wp
               v_c = 0.0_wp
            end if
            sst = hTr_T(ii, jj, nz)/(h_layer(ii, jj, nz) + H_DIV_EPS)
            ! Bulk formulae and longwave emission take a POTENTIAL
            ! temperature; the model T is CT under cons_abs.  The heat
            ! content of the precipitation below keeps the model T
            ! (ocean-frame enthalpy with c_p).
            sst_f = sst
            if (cons_abs) then
               sst_f = pt_from_ct(hTr_S(ii, jj, nz)/(h_layer(ii, jj, nz) + H_DIV_EPS), sst)
            end if
            call bulk_flux_column(u_g, v_g, u_c, v_c, tas(ii, jj), huss(ii, jj), psl(ii, jj), sst_f, &
                                  z_wind, z_ta, n_iter, &
                                  taux, tauy, sens, lat, evp, cd, ch, ce, ustar)
            sst_k = sst_f + BULK_T0_KELVIN
            lw_up = bulk_flux_longwave_up(sst_k, emissivity)

            taux_c(ii, jj) = ow*taux
            tauy_c(ii, jj) = ow*tauy
            q_sens(ii, jj) = ow*sens
            q_lat(ii, jj) = ow*lat
            evap(ii, jj) = -ow*evp
            q_sw(ii, jj) = ow*(1.0_wp - ocean_albedo)*rsds(ii, jj)
            q_lw(ii, jj) = ow*(rlds(ii, jj) - lw_up)
            lprec(ii, jj) = ow*prra(ii, jj)
            fprec(ii, jj) = ow*prsn(ii, jj)
            p_surf_atm(ii, jj) = psl(ii, jj)
            heat_content_lprec(ii, jj) = cp*sst*lprec(ii, jj)
            heat_content_fprec(ii, jj) = cp*sst*fprec(ii, jj)
         else
            taux_c(ii, jj) = 0.0_wp
            tauy_c(ii, jj) = 0.0_wp
            q_sens(ii, jj) = 0.0_wp
            q_lat(ii, jj) = 0.0_wp
            evap(ii, jj) = 0.0_wp
            q_sw(ii, jj) = 0.0_wp
            q_lw(ii, jj) = 0.0_wp
            lprec(ii, jj) = 0.0_wp
            fprec(ii, jj) = 0.0_wp
            p_surf_atm(ii, jj) = 0.0_wp
            heat_content_lprec(ii, jj) = 0.0_wp
            heat_content_fprec(ii, jj) = 0.0_wp
         end if
      end do
   end subroutine bulk_forcing_core_impl

   pure subroutine bulk_forcing_stress_face_impl(nx, ny, taux_c, tauy_c, tau_x, tau_y)
      !! Cell-centred stress -> C-grid face average, ALL faces
      !! (interior and the two tile-edge faces per axis), inverting
      !! `ocean_surfstress_derived_impl`'s own face -> centre average
      !! (`rdb_ocean_surface_stress`) so the two are each other's
      !! round-trip on a uniform field. Edge faces (`i=1`/`nx+1`,
      !! `j=1`/`ny+1`) are one-sided (only one neighbouring cell exists
      !! in the array) -- excluded from the interior loop by RANGE, per
      !! CLAUDE.md's vectoriser-edge-bug rule, never by an in-loop
      !! `i==1 .or. i==nx+1` guard.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: taux_c(nx, ny), tauy_c(nx, ny)
      real(wp), intent(inout) :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=2:nx)
         tau_x(i, j) = 0.5_wp*(taux_c(i - 1, j) + taux_c(i, j))
      end do
      do concurrent(j=1:ny)
         tau_x(1, j) = taux_c(1, j)
         tau_x(nx + 1, j) = taux_c(nx, j)
      end do

      do concurrent(j=2:ny, i=1:nx)
         tau_y(i, j) = 0.5_wp*(tauy_c(i, j - 1) + tauy_c(i, j))
      end do
      do concurrent(i=1:nx)
         tau_y(i, 1) = tauy_c(i, 1)
         tau_y(i, ny + 1) = tauy_c(i, ny)
      end do
   end subroutine bulk_forcing_stress_face_impl

end module rdb_ocean_bulk_forcing
