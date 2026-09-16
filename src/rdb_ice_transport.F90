!! Category ice/snow transport + compress_ice (SIS2 port, PR 4b).
module rdb_ice_transport
   !! Port of SIS2's DEFAULT (velocity, non-merged) `ice_cat_transport`
   !! (`SIS_transport.F90:127`) + `finish_ice_transport` (`:255`) +
   !! `compress_ice` (`:898`), grounded per `SPEC_ice-pr4b-transport.md`.
   !! Source citations below are `SIS_transport.F90` / `SIS_continuity.F90`
   !! / `SIS_tracer_advect.F90` unless noted; when this module and the
   !! SIS2 source disagree, the source (re-verified against the running
   !! SIS2 tree for this port) wins.
   !!
   !! **What SIS2 actually runs** (`MERGED_CONTINUITY` defaults `.false.`,
   !! `SIS_dyn_trans.F90:2377`): per-medium continuity driven by FACE
   !! VELOCITIES (`uc`/`vc`), not the merged/proportionate-split variant.
   !! `ice_continuity` (`SIS_continuity.F90:69`) reconstructs a PPM
   !! parabola on the CATEGORY-SUMMED mass, computes ONE total-mass face
   !! transport from that parabola, then splits it into per-category
   !! flux proportionally to each category's share of the summed mass
   !! (`SIS_continuity.F90:1184-1215`) — there is no separate "split"
   !! pass. Snow rides its OWN independent PPM solve on summed snow mass,
   !! then is masked to zero wherever the CO-LOCATED ice flux is exactly
   !! zero (`SIS_continuity.F90:1210-1215`; the `frac_neglect` second
   !! masking clause is dead at SIS2's own default `frac_neglect=0`, so it
   !! is not ported). `x_first` is fixed by `FIRST_DIRECTION` (default 0);
   !! the second (y) pass reads the ALREADY-UPDATED post-x-pass masses
   !! (`SIS_continuity.F90:170-216`).
   !!
   !! Algorithm per outer call (`ice_transport_step`):
   !!   Phase 0 — gates (`is_init`, `ncat==1`) + velocity sampling +
   !!     zero-velocity exact no-op.
   !!   Phase 1 — IST -> CAS (`ice_state_to_cell_ave_state`, `:465`):
   !!     `mca_ice(c) = part_size(c)*m_ice(c)`, ditto snow.
   !!   Phase 2 — `adv_substeps` iterations, each an x-pass then y-pass;
   !!     each pass does (a) total-mass PPM + proportionate ice-flux
   !!     split, (b) same for snow + the co-located mask, (c) PCM tracer
   !!     riding (mass-weighted, BEFORE the mass update, with SIS2's
   !!     `H_NEGLECT` conditioning guard for thin remainders), (d) the
   !!     mass update, then a fail-loud positivity/orphan-snow reduction.
   !!   Phase 3 — CAS -> IST (`cell_ave_state_to_ice_state`, `:540`):
   !!     re-derive `part_size` by division, with the pre-floor + optional
   !!     thin-ice rolling + general floor.  `part_size(0)` may go
   !!     negative here — compress (Phase 4) is what fixes it.
   !!   Phase 4 — `compress_ice` (`:898`, no-ridge default, ponds
   !!     dropped): thinnest-first cascade that returns `part_size(0)` to
   !!     >= 0 by compacting / promoting categories, tracer-merging on
   !!     transfer.
   !!   Phase 5 — `ice_adjust_categories` (PR 4a, unchanged): restore the
   !!     ITD partition.
   !!
   !! **GPU race-hazard note (Phase 2c, tracer riding).** A per-cell
   !! update that reads a NEIGHBOUR cell's intensive value (the upwind
   !! donor of the west/east face) while another loop iteration may be
   !! WRITING that same neighbour is a `do concurrent` race — iteration
   !! order is unspecified.  This module never reads a neighbour's `val`
   !! directly: it first runs a FACE-indexed gather kernel
   !! (`ice_gather_flux_{x,y}[_layer]_impl`) that computes `tr_flux(I) =
   !! uh(I)*val(donor of I)` — reading `val` ONLY at that face's own donor
   !! cell, writing to the SEPARATE `tr_flux_x_work`/`tr_flux_y_work` face
   !! buffer (never `val` itself) — then a CELL-indexed update kernel
   !! (`ice_ride_update_{x,y}[_layer]_impl`) that reads only the face
   !! buffer plus its OWN cell's prior mass/value and writes
   !! `val(i,j,c[,l])`.  This mirrors `rdb_continuity`'s `Tr_face_left_x`
   !! gather/scatter split (`tracer_advect_zonal_one_impl`) verbatim in
   !! spirit.
   !!
   !! **Reuse contract.** The five SIS2-equivalent PPM helpers
   !! (`ppm_mirror_h`, `volcfl_face`, `ppm_limited_slope`,
   !! `ppm_cell_limiter`, `ppm_limit_pos`) are promoted to production
   !! public in `rdb_continuity` and called verbatim here — their bodies
   !! are NOT duplicated (SPEC §2).
   !!
   !! **Documented divergences from SIS2** (kept in sync with
   !! `SPEC_ice-pr4b-transport.md` §10):
   !!   D1: continuity scheme fixed to Roundabout's PPM (H3-style limited
   !!     edges + CW84 + limit_pos) rather than SIS2's in-code legacy
   !!     default `UPWIND_2D`; matches SIS2's modern PPM configs in
   !!     structure.
   !!   D2: no massless-category `mH` fill (SIS2 `:509-520`) — moot under
   !!     PCM tracer riding (a massless category's scalar value is never
   !!     read: its face fluxes are always zero because `mca=0`).
   !!   D3: tracer riding is PCM (piecewise-constant upwind), one of
   !!     SIS2's four schemes (`SIS_TRACER_ADVECTION_SCHEME=PCM`) — not
   !!     the PPM:H3 modern configs use.  Upgrade path open.
   !!   D4: no melt ponds (Roundabout carries none).
   !!   D5: fixed x-first directional split (SIS2 `FIRST_DIRECTION`
   !!     default 0; no alternation).
   !!   D6: zero-velocity early-exit and `ncat==1` early-exit are Roundabout
   !!     bit-identity contract additions; SIS2 has neither.
   !!   D7: fail-loud negativity/orphan-snow detection via a post-pass
   !!     device REDUCTION + driver abort, instead of SIS2's in-loop
   !!     FATALs (GPU portability — a `do concurrent` kernel cannot abort
   !!     mid-loop).  Also STRICTLY SAFER than SIS2 at the top-category
   !!     compaction: `part(ncat) <= excess` (negative- or, at the exact
   !!     tie, zero-denominator `f`) is routed into the fail-loud branch
   !!     (SIS2 has that hole).
   !!   D8: `compress_ice` tracer-merges INSIDE the thinnest-first cascade
   !!     at each transfer site, whereas SIS2 defers to
   !!     `advect_tracers_thicker` AFTER the category k-loop — equivalent
   !!     to round-off (each boundary visited once in a fixed order with
   !!     running masses; same argument as the PR-4a `ice_adjust_categories`
   !!     merge).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_itd, only: ice_adjust_categories
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_continuity, only: ppm_mirror_h, volcfl_face, ppm_limited_slope, &
                             ppm_cell_limiter, ppm_limit_pos
   implicit none
   private

   public :: ice_transport_step
   public :: ice_transport_compress_cell
   public :: H_NEGLECT_ICE_TRANSPORT
   public :: MASS_NEGLECT_ICE_TRANSPORT

   real(wp), parameter :: H_NEGLECT_ICE_TRANSPORT = 1.0e-30_wp
      !! SIS2 `IG%H_subroundoff` role (`SIS_tracer_advect.F90`
      !! `advect_scalar_x` conditioning guard): a per-CELL mass floor
      !! (kg/m² of cell area) below which the tracer-riding division is
      !! reconditioned rather than divided directly.  Deliberately
      !! distinct from `rdb_constants`' `H_VANISHED`/`H_DIV_EPS` — this
      !! guard's algebra (proportional `h_add` redistribution across the
      !! old mass + both face transports) is SIS2-specific and belongs
      !! with the kernel that uses it.
   real(wp), parameter :: MASS_NEGLECT_ICE_TRANSPORT = 1.0e-60_wp
      !! SIS2 `compress_ice`'s `mass_neglect` (`SIS_transport.F90:959`) —
      !! the cell-average-mass gate below which a category is treated as
      !! contributing nothing to a transfer (guards the thickness-merge
      !! division, not a physical floor).

contains

   subroutine ice_transport_step(grid, metrics, ms, ice, dt, adv_substeps, roll_factor, ok)
      !! Entry point (host orchestration; kernels inside).  Outer-shim +
      !! flat-impl: every phase below dispatches to a `pure` `_impl`
      !! kernel; this routine only sequences them and owns the
      !! host-visible `ok` fail-loud signal.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
         !! READ-ONLY: `wet_mask` (physical-cell gate) + the surface-layer
         !! face velocities the v1 sampler reads.
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt
         !! The thermo-step dt (cadence = the ice slow step).
      integer, intent(in) :: adv_substeps
         !! SIS2 `NSTEPS_ADV` (`&ocean_ice_nml adv_substeps`), >= 1.
      real(wp), intent(in) :: roll_factor
         !! SIS2 `SEA_ICE_ROLL_FACTOR` (`&ocean_ice_nml roll_factor`).
         !! 0 disables rolling.
      logical, intent(out) :: ok
         !! `.false.` => the caller (driver) must abort fail-loud
         !! (conservation/positivity violation, or a compress-time
         !! consistency failure SIS2 would FATAL on).

      integer :: nx, ny, nz, nghost, n
      real(wp) :: dt_adv, vmax

      ok = .true.
      if (.not. ice%is_init) return
      if (ice%ncat == 1) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nghost = grid%nghost

      ! ---- Phase 0: sample the ocean surface-layer velocity onto the
      ! ice C-grid faces (v1 interim filler).  PR 5: skipped when EVP
      ! dynamics is on — ice_evp_step (driver, runs BEFORE this call)
      ! already wrote u_ice/v_ice directly; the sampler would clobber
      ! them with the ocean surface velocity. ----
      if (.not. ice%dynamics) then
         call ice_sample_velocity_impl(ms%u_face_x_layer(:, :, nz), ms%v_face_y_layer(:, :, nz), &
                                       metrics%wet_u, metrics%wet_v, ice%u_ice, ice%v_ice, nx, ny)
      end if

      ! ---- Zero-velocity exact no-op (D6): the CAS<->IST round trip
      ! re-derives part_size by division and would otherwise perturb
      ! last bits even under zero flow. ----
      call ice_max_speed_impl(ice%u_ice, ice%v_ice, nghost, grid%nx_phys, grid%ny_phys, &
                              nx, ny, vmax)
      if (vmax == 0.0_wp) return

      ! ---- Phase 1: IST -> CAS ----
      call ice_ist_to_cas_impl(ms%wet_mask, ice%part_size, ice%m_ice, ice%m_snow, &
                               ice%mca_ice, ice%mca_snow, nghost, ice%ncat, nx, ny)

      ! ---- Phase 2: adv_substeps advective iterations, x-pass then y-pass ----
      dt_adv = dt/real(max(adv_substeps, 1), wp)
      do n = 1, adv_substeps
         call ice_pass_x(grid, metrics, ice, dt_adv, ok)
         if (.not. ok) return
         call ice_pass_y(grid, metrics, ice, dt_adv, ok)
         if (.not. ok) return
      end do

      ! ---- Phase 3: CAS -> IST ----
      call ice_cas_to_ist_impl(ms%wet_mask, metrics%areaT, ice%mca_ice, ice%mca_snow, &
                               ice%part_size, ice%m_ice, ice%m_snow, ice%mh_lim, roll_factor, &
                               nghost, ice%ncat, nx, ny)

      ! ---- Phase 4: compress_ice ----
      call ice_compress_impl(ms%wet_mask, ice%part_size, ice%m_ice, ice%m_snow, &
                             ice%enth_ice, ice%enth_snow, ice%sal_ice, ice%mh_lim, &
                             nghost, ice%ncat, ice%nk_ice, nx, ny, ok)
      if (.not. ok) return

      ! ---- Phase 5: recategorize (PR 4a, unchanged) ----
      call ice_adjust_categories(grid, ms, ice)
   end subroutine ice_transport_step

   ! ======================================================================
   ! Phase 0: velocity sampling + zero-velocity gate
   ! ======================================================================

   pure subroutine ice_sample_velocity_impl(u_surf, v_surf, wet_u, wet_v, u_ice, v_ice, nx, ny)
      !! `u_ice(i,j) = u_surf(i,j)*wet_u(i,j)`, ditto v.  v1 interim
      !! filler (SPEC §5 Phase 0) — PR 5 EVP replaces this call with a
      !! dynamics solve writing the SAME `u_ice`/`v_ice` faces.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: u_surf(nx + 1, ny)
      real(wp), intent(in) :: v_surf(nx, ny + 1)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(in) :: wet_v(nx, ny + 1)
      real(wp), intent(out) :: u_ice(nx + 1, ny)
      real(wp), intent(out) :: v_ice(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         u_ice(i, j) = u_surf(i, j)*wet_u(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         v_ice(i, j) = v_surf(i, j)*wet_v(i, j)
      end do
   end subroutine ice_sample_velocity_impl

   pure subroutine ice_max_speed_impl(u_ice, v_ice, nghost, nx_phys, ny_phys, nx, ny, vmax)
      !! Max |u_ice|/|v_ice| over PHYSICAL faces only (the array-edge
      !! ghost faces carry no meaning here).  `!$acc parallel loop
      !! reduction` (inert comment on non-OpenACC compilers) — GPU-safe
      !! max reduction for the zero-velocity exact no-op gate.
      integer, intent(in) :: nghost, nx_phys, ny_phys, nx, ny
      real(wp), intent(in) :: u_ice(nx + 1, ny)
      real(wp), intent(in) :: v_ice(nx, ny + 1)
      real(wp), intent(out) :: vmax
      integer :: i, j, i_lo, i_hi, j_lo, j_hi

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      vmax = 0.0_wp
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi + 1) reduce(max:vmax)
         vmax = max(vmax, abs(u_ice(i, j)))
      end do
      do concurrent(j=j_lo:j_hi + 1, i=i_lo:i_hi) reduce(max:vmax)
         vmax = max(vmax, abs(v_ice(i, j)))
      end do
   end subroutine ice_max_speed_impl

   ! ======================================================================
   ! Phase 1: IST -> CAS
   ! ======================================================================

   pure subroutine ice_ist_to_cas_impl(wet_mask, part_size, m_ice, m_snow, mca_ice, mca_snow, &
                                       nghost, ncat, nx, ny)
      !! SIS2 `ice_state_to_cell_ave_state` (`:465`): `mca(c) =
      !! part_size(c)*m(c)` per category, physical cells only.  Ghost
      !! cells are zeroed (never read as donors in the flux kernels below
      !! — wet mirroring + wall zeroing exclude them).
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      real(wp), intent(in) :: m_snow(nx, ny, ncat)
      real(wp), intent(out) :: mca_ice(nx, ny, ncat)
      real(wp), intent(out) :: mca_snow(nx, ny, ncat)
      integer :: i, j, c, i_lo, i_hi, j_lo, j_hi

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=1:ny, i=1:nx, c=1:ncat)
         mca_ice(i, j, c) = 0.0_wp
         mca_snow(i, j, c) = 0.0_wp
      end do
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi, c=1:ncat)
         if (wet_mask(i, j) > 0.5_wp) then
            mca_ice(i, j, c) = part_size(i, j, c)*m_ice(i, j, c)
            mca_snow(i, j, c) = part_size(i, j, c)*m_snow(i, j, c)
         end if
      end do
   end subroutine ice_ist_to_cas_impl

   ! ======================================================================
   ! Phase 2: one directional pass (x then y)
   ! ======================================================================

   pure subroutine ice_pass_x(grid, metrics, ice, dt_adv, ok)
      !! One zonal pass: total-mass PPM + proportionate ice-flux split,
      !! the snow twin + co-located mask, PCM tracer riding (gather then
      !! cell update), the mass update (AFTER every ride reads the
      !! pre-update mass), and the post-pass validity reduction.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt_adv
      logical, intent(out) :: ok
      integer :: nx, ny, nghost, l

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      call ice_cat_flux_x_impl(metrics%wet_T, metrics%dy_cu, metrics%idxT, ice%u_ice, &
                               ice%mca_ice, ice%htot_work, ice%hl_x_work, ice%hr_x_work, &
                               ice%uhtot_work, ice%uh_ice, dt_adv, nghost, grid%nx_phys, &
                               ice%ncat, nx, ny)
      call ice_cat_flux_x_impl(metrics%wet_T, metrics%dy_cu, metrics%idxT, ice%u_ice, &
                               ice%mca_snow, ice%htot_work, ice%hl_x_work, ice%hr_x_work, &
                               ice%uhtot_work, ice%uh_snow, dt_adv, nghost, grid%nx_phys, &
                               ice%ncat, nx, ny)
      call ice_mask_snow_by_ice_impl(ice%uh_ice, ice%uh_snow, nx + 1, ny, ice%ncat)

      ! ---- Tracer riding (Phase 2c), BEFORE the mass update ----
      call ice_gather_flux_x_impl(ice%uh_ice, ice%m_ice, ice%tr_flux_x_work, ice%ncat, nx, ny)
      call ice_ride_update_x_impl(metrics%iareaT, ice%uh_ice, ice%tr_flux_x_work, ice%mca_ice, &
                                  ice%m_ice, dt_adv, ice%ncat, nx, ny)
      do l = 1, ice%nk_ice
         call ice_gather_flux_x_layer_impl(ice%uh_ice, ice%enth_ice, l, ice%tr_flux_x_work, &
                                           ice%ncat, ice%nk_ice, nx, ny)
         call ice_ride_update_x_layer_impl(metrics%iareaT, ice%uh_ice, ice%tr_flux_x_work, &
                                           ice%mca_ice, ice%enth_ice, l, dt_adv, ice%ncat, &
                                           ice%nk_ice, nx, ny)
         call ice_gather_flux_x_layer_impl(ice%uh_ice, ice%sal_ice, l, ice%tr_flux_x_work, &
                                           ice%ncat, ice%nk_ice, nx, ny)
         call ice_ride_update_x_layer_impl(metrics%iareaT, ice%uh_ice, ice%tr_flux_x_work, &
                                           ice%mca_ice, ice%sal_ice, l, dt_adv, ice%ncat, &
                                           ice%nk_ice, nx, ny)
      end do
      call ice_gather_flux_x_layer_impl(ice%uh_snow, ice%enth_snow, 1, ice%tr_flux_x_work, &
                                        ice%ncat, 1, nx, ny)
      call ice_ride_update_x_layer_impl(metrics%iareaT, ice%uh_snow, ice%tr_flux_x_work, &
                                        ice%mca_snow, ice%enth_snow, 1, dt_adv, ice%ncat, &
                                        1, nx, ny)

      ! ---- Mass update (Phase 2d) ----
      call ice_mass_update_x_impl(metrics%iareaT, ice%uh_ice, ice%mca_ice, dt_adv, &
                                  nghost, grid%nx_phys, ice%ncat, nx, ny)
      call ice_mass_update_x_impl(metrics%iareaT, ice%uh_snow, ice%mca_snow, dt_adv, &
                                  nghost, grid%nx_phys, ice%ncat, nx, ny)

      call ice_validity_reduce_impl(ice%mca_ice, ice%mca_snow, nghost, grid%nx_phys, &
                                    grid%ny_phys, ice%ncat, nx, ny, ok)
   end subroutine ice_pass_x

   pure subroutine ice_pass_y(grid, metrics, ice, dt_adv, ok)
      !! Meridional twin of `ice_pass_x`.  Reads the POST-x-pass masses
      !! (SIS2 `SIS_continuity.F90:170-216`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt_adv
      logical, intent(out) :: ok
      integer :: nx, ny, nghost, l

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      call ice_cat_flux_y_impl(metrics%wet_T, metrics%dx_cv, metrics%idyT, ice%v_ice, &
                               ice%mca_ice, ice%htot_work, ice%hl_y_work, ice%hr_y_work, &
                               ice%vhtot_work, ice%vh_ice, dt_adv, nghost, grid%ny_phys, &
                               ice%ncat, nx, ny)
      call ice_cat_flux_y_impl(metrics%wet_T, metrics%dx_cv, metrics%idyT, ice%v_ice, &
                               ice%mca_snow, ice%htot_work, ice%hl_y_work, ice%hr_y_work, &
                               ice%vhtot_work, ice%vh_snow, dt_adv, nghost, grid%ny_phys, &
                               ice%ncat, nx, ny)
      call ice_mask_snow_by_ice_impl(ice%vh_ice, ice%vh_snow, nx, ny + 1, ice%ncat)

      call ice_gather_flux_y_impl(ice%vh_ice, ice%m_ice, ice%tr_flux_y_work, ice%ncat, nx, ny)
      call ice_ride_update_y_impl(metrics%iareaT, ice%vh_ice, ice%tr_flux_y_work, ice%mca_ice, &
                                  ice%m_ice, dt_adv, ice%ncat, nx, ny)
      do l = 1, ice%nk_ice
         call ice_gather_flux_y_layer_impl(ice%vh_ice, ice%enth_ice, l, ice%tr_flux_y_work, &
                                           ice%ncat, ice%nk_ice, nx, ny)
         call ice_ride_update_y_layer_impl(metrics%iareaT, ice%vh_ice, ice%tr_flux_y_work, &
                                           ice%mca_ice, ice%enth_ice, l, dt_adv, ice%ncat, &
                                           ice%nk_ice, nx, ny)
         call ice_gather_flux_y_layer_impl(ice%vh_ice, ice%sal_ice, l, ice%tr_flux_y_work, &
                                           ice%ncat, ice%nk_ice, nx, ny)
         call ice_ride_update_y_layer_impl(metrics%iareaT, ice%vh_ice, ice%tr_flux_y_work, &
                                           ice%mca_ice, ice%sal_ice, l, dt_adv, ice%ncat, &
                                           ice%nk_ice, nx, ny)
      end do
      call ice_gather_flux_y_layer_impl(ice%vh_snow, ice%enth_snow, 1, ice%tr_flux_y_work, &
                                        ice%ncat, 1, nx, ny)
      call ice_ride_update_y_layer_impl(metrics%iareaT, ice%vh_snow, ice%tr_flux_y_work, &
                                        ice%mca_snow, ice%enth_snow, 1, dt_adv, ice%ncat, &
                                        1, nx, ny)

      call ice_mass_update_y_impl(metrics%iareaT, ice%vh_ice, ice%mca_ice, dt_adv, &
                                  nghost, grid%ny_phys, ice%ncat, nx, ny)
      call ice_mass_update_y_impl(metrics%iareaT, ice%vh_snow, ice%mca_snow, dt_adv, &
                                  nghost, grid%ny_phys, ice%ncat, nx, ny)

      call ice_validity_reduce_impl(ice%mca_ice, ice%mca_snow, nghost, grid%nx_phys, &
                                    grid%ny_phys, ice%ncat, nx, ny, ok)
   end subroutine ice_pass_y

   ! ----------------------------------------------------------------------
   ! (a) Total-mass PPM face transport + proportionate category split
   ! ----------------------------------------------------------------------

   pure subroutine ice_cat_flux_x_impl(wet_T, dy_cu, idxT, u_ice, mca, htot_work, hl_x_work, &
                                       hr_x_work, uhtot_work, uh_out, dt_adv, nghost, nx_phys, &
                                       ncat, nx, ny)
      !! SIS2 `zonal_mass_flux` (`SIS_continuity.F90:1064`): PPM
      !! reconstruction of the category-SUMMED mass `htot`, ONE total
      !! face transport `uhtot` from the swept-volume parabola integral
      !! (`volcfl_face`, bit-for-bit the SIS2 face expression), then the
      !! PROPORTIONATE split `uh(c) = uhtot*mca(donor,c)*I_htot(donor)`
      !! (`SIS_continuity.F90:1199-1205`, Adcroft reciprocal — `I_htot=0`
      !! when `htot(donor)<=0`).  Stencil + edge STORAGE CONVENTION +
      !! swept-face orientation copied verbatim from `continuity_zonal_flux`
      !! (`rdb_continuity.F90:871-948`): `hl_x_work(i)` == `h_face_left_x(i)`
      !! is the value AT east face `i` from the LEFT cell `i-1` (that
      !! cell's OWN downwind edge); `hr_x_work(i)` == `h_face_right_x(i)`
      !! is from the RIGHT cell `i` (its OWN left edge).  H3-style limited
      !! edges + CW84 + `ppm_limit_pos` (SIS2 runs `PPM_limit_pos`
      !! UNCONDITIONALLY on the PD scheme, so it is not optional here).
      !! CFL metric: SIS2's shipped default `vol_CFL=.false.` uses
      !! `CFL = |u|*dt*IdxT(donor)` (`SIS_continuity.F90:1165`) — the
      !! T-cell inverse spacing, NOT the `dy_cu*iareaT` swept-area ratio
      !! (== SIS2's `vol_CFL=.true.` variant).  Identical on uniform
      !! Cartesian; correct on spherical / anisotropic.
      integer, intent(in) :: nghost, nx_phys, ncat, nx, ny
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: idxT(nx, ny)
      real(wp), intent(in) :: u_ice(nx + 1, ny)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: htot_work(nx, ny)
      real(wp), intent(inout) :: hl_x_work(nx + 1, ny)
      real(wp), intent(inout) :: hr_x_work(nx + 1, ny)
      real(wp), intent(inout) :: uhtot_work(nx + 1, ny)
      real(wp), intent(out) :: uh_out(nx + 1, ny, ncat)
      real(wp), intent(in) :: dt_adv

      integer :: i, j, c
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right
      real(wp) :: hm2, hm1, h0, hp1, hp2
      real(wp) :: u, cfl, dh_d, curv3_d, h_face, i_htot

      ! ---- Category-summed mass ----
      do concurrent(j=1:ny, i=1:nx)
         htot_work(i, j) = sum(mca(i, j, :))
      end do

      ! ---- PPM reconstruction (interior 5-point stencil) ----
      ! Cell i's LEFT edge -> hr_x_work(i) (right-cell state at face i);
      ! cell i's RIGHT edge -> hl_x_work(i+1) (left-cell state at face
      ! i+1).  Exactly `continuity_compute_fluxes`'s write pattern.
      do concurrent(j=1:ny, i=3:nx - 2) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         h0 = htot_work(i, j)
         hm1 = ppm_mirror_h(htot_work(i - 1, j), h0, wet_T(i - 1, j))
         hp1 = ppm_mirror_h(htot_work(i + 1, j), h0, wet_T(i + 1, j))
         hm2 = ppm_mirror_h(htot_work(i - 2, j), hm1, wet_T(i - 2, j))
         hp2 = ppm_mirror_h(htot_work(i + 2, j), hp1, wet_T(i + 2, j))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*wet_T(i - 1, j)*wet_T(i, j)*wet_T(i + 1, j)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         call ppm_limit_pos(h0, h_left, h_right, 0.0_wp)
         hr_x_work(i, j) = h_left
         hl_x_work(i + 1, j) = h_right
      end do
      ! First-order fallback within 2 cells of either array edge
      ! (face-shaped writes, in bounds: `nx+1` is the last valid face).
      do concurrent(j=1:ny)
         hl_x_work(1, j) = htot_work(1, j)
         hr_x_work(1, j) = htot_work(1, j)
         hl_x_work(2, j) = htot_work(1, j)
         hr_x_work(2, j) = htot_work(2, j)
         hl_x_work(3, j) = htot_work(2, j)
         hr_x_work(3, j) = htot_work(2, j)
         hr_x_work(nx - 1, j) = htot_work(nx - 1, j)
         hl_x_work(nx, j) = htot_work(nx - 1, j)
         hr_x_work(nx, j) = htot_work(nx, j)
         hl_x_work(nx + 1, j) = htot_work(nx, j)
         hr_x_work(nx + 1, j) = htot_work(nx, j)
      end do

      ! ---- Total-mass face transport (swept-volume PPM integral) ----
      ! Donor-edge orientation copied verbatim from
      ! `continuity_zonal_flux` (`rdb_continuity.F90:923-948`):
      !   u>0 (donor i-1): edge = hl_x_work(i) [cell i-1's downwind
      !     edge]; dh = hr_x_work(i-1) - hl_x_work(i); curv3 =
      !     hr_x_work(i-1) + hl_x_work(i) - 2*htot(i-1).
      !   u<0 (donor i):   edge = hr_x_work(i) [cell i's downwind edge];
      !     dh = hl_x_work(i+1) - hr_x_work(i); curv3 = hr_x_work(i) +
      !     hl_x_work(i+1) - 2*htot(i).
      ! No `max(...,0)` clamp: a `ppm_limit_pos`'d parabola has a
      ! non-negative swept mean for CFL <= 1, so `volcfl_face` stays >= 0
      ! here — the earlier negative values were the SIGNATURE of the
      ! wrong (now-fixed) donor-edge stencil feeding mixed-cell inputs,
      ! not a property of the scheme (neither SIS2 nor
      ! continuity_zonal_flux clamps).
      do concurrent(j=1:ny, i=2:nx) local(u, cfl, dh_d, curv3_d, h_face)
         u = u_ice(i, j)
         if (u > 0.0_wp) then
            cfl = u*dt_adv*idxT(i - 1, j)
            h_face = hl_x_work(i, j)
            dh_d = hr_x_work(i - 1, j) - h_face
            curv3_d = hr_x_work(i - 1, j) + h_face - 2.0_wp*htot_work(i - 1, j)
            uhtot_work(i, j) = dy_cu(i, j)*u*volcfl_face(h_face, dh_d, curv3_d, cfl)
         else if (u < 0.0_wp) then
            cfl = (-u)*dt_adv*idxT(i, j)
            h_face = hr_x_work(i, j)
            dh_d = hl_x_work(i + 1, j) - h_face
            curv3_d = h_face + hl_x_work(i + 1, j) - 2.0_wp*htot_work(i, j)
            uhtot_work(i, j) = dy_cu(i, j)*u*volcfl_face(h_face, dh_d, curv3_d, cfl)
         else
            uhtot_work(i, j) = 0.0_wp
         end if
      end do
      do concurrent(j=1:ny)
         uhtot_work(1, j) = 0.0_wp
         uhtot_work(nx + 1, j) = 0.0_wp
      end do
      ! Physical wall faces (not just array edges) — mirrors
      ! `continuity_zonal_flux`'s wall-zeroing rationale: a moving ocean
      ! surface layer can leave a nonzero sampled velocity at the
      ! interior physical wall even on a closed-boundary configuration.
      do concurrent(j=1:ny)
         uhtot_work(nghost + 1, j) = 0.0_wp
         uhtot_work(nghost + nx_phys + 1, j) = 0.0_wp
      end do

      ! ---- Proportionate category split (Adcroft reciprocal) ----
      do concurrent(j=1:ny, i=1:nx + 1, c=1:ncat) local(i_htot)
         if (uhtot_work(i, j) == 0.0_wp) then
            uh_out(i, j, c) = 0.0_wp
         else if (u_ice(i, j) >= 0.0_wp) then
            if (htot_work(i - 1, j) > 0.0_wp) then
               i_htot = 1.0_wp/htot_work(i - 1, j)
            else
               i_htot = 0.0_wp
            end if
            uh_out(i, j, c) = uhtot_work(i, j)*mca(i - 1, j, c)*i_htot
         else
            if (htot_work(i, j) > 0.0_wp) then
               i_htot = 1.0_wp/htot_work(i, j)
            else
               i_htot = 0.0_wp
            end if
            uh_out(i, j, c) = uhtot_work(i, j)*mca(i, j, c)*i_htot
         end if
      end do
   end subroutine ice_cat_flux_x_impl

   pure subroutine ice_cat_flux_y_impl(wet_T, dx_cv, idyT, v_ice, mca, htot_work, hl_y_work, &
                                       hr_y_work, vhtot_work, vh_out, dt_adv, nghost, ny_phys, &
                                       ncat, nx, ny)
      !! Meridional twin of `ice_cat_flux_x_impl`.  Same face-indexed edge
      !! STORAGE + swept orientation as `continuity_meridional_flux`
      !! (`rdb_continuity.F90:1159-1202`): `hl_y_work(i,j)` == north-face
      !! `h_face_left_y` (from the SOUTH cell `j-1`), `hr_y_work(i,j)` ==
      !! `h_face_right_y` (from the NORTH cell `j`).  CFL uses
      !! `IdyT(donor)` (SIS2 `vol_CFL=.false.` default).
      integer, intent(in) :: nghost, ny_phys, ncat, nx, ny
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idyT(nx, ny)
      real(wp), intent(in) :: v_ice(nx, ny + 1)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: htot_work(nx, ny)
      real(wp), intent(inout) :: hl_y_work(nx, ny + 1)
      real(wp), intent(inout) :: hr_y_work(nx, ny + 1)
      real(wp), intent(inout) :: vhtot_work(nx, ny + 1)
      real(wp), intent(out) :: vh_out(nx, ny + 1, ncat)
      real(wp), intent(in) :: dt_adv

      integer :: i, j, c
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right
      real(wp) :: hm2, hm1, h0, hp1, hp2
      real(wp) :: v, cfl, dh_d, curv3_d, h_face, i_htot

      do concurrent(j=1:ny, i=1:nx)
         htot_work(i, j) = sum(mca(i, j, :))
      end do

      do concurrent(j=3:ny - 2, i=1:nx) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         h0 = htot_work(i, j)
         hm1 = ppm_mirror_h(htot_work(i, j - 1), h0, wet_T(i, j - 1))
         hp1 = ppm_mirror_h(htot_work(i, j + 1), h0, wet_T(i, j + 1))
         hm2 = ppm_mirror_h(htot_work(i, j - 2), hm1, wet_T(i, j - 2))
         hp2 = ppm_mirror_h(htot_work(i, j + 2), hp1, wet_T(i, j + 2))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*wet_T(i, j - 1)*wet_T(i, j)*wet_T(i, j + 1)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         call ppm_limit_pos(h0, h_left, h_right, 0.0_wp)
         hr_y_work(i, j) = h_left
         hl_y_work(i, j + 1) = h_right
      end do
      do concurrent(i=1:nx)
         hl_y_work(i, 1) = htot_work(i, 1)
         hr_y_work(i, 1) = htot_work(i, 1)
         hl_y_work(i, 2) = htot_work(i, 1)
         hr_y_work(i, 2) = htot_work(i, 2)
         hl_y_work(i, 3) = htot_work(i, 2)
         hr_y_work(i, 3) = htot_work(i, 2)
         hr_y_work(i, ny - 1) = htot_work(i, ny - 1)
         hl_y_work(i, ny) = htot_work(i, ny - 1)
         hr_y_work(i, ny) = htot_work(i, ny)
         hl_y_work(i, ny + 1) = htot_work(i, ny)
         hr_y_work(i, ny + 1) = htot_work(i, ny)
      end do

      ! Donor-edge orientation from `continuity_meridional_flux`:
      !   v>0 (donor j-1): edge = hl_y_work(i,j); dh = hr_y_work(i,j-1)
      !     - hl_y_work(i,j); curv3 = hr_y_work(i,j-1) + hl_y_work(i,j)
      !     - 2*htot(i,j-1).
      !   v<0 (donor j):   edge = hr_y_work(i,j); dh = hl_y_work(i,j+1)
      !     - hr_y_work(i,j); curv3 = hr_y_work(i,j) + hl_y_work(i,j+1)
      !     - 2*htot(i,j).
      ! No `max(...,0)` clamp (see the x-pass twin).
      do concurrent(j=2:ny, i=1:nx) local(v, cfl, dh_d, curv3_d, h_face)
         v = v_ice(i, j)
         if (v > 0.0_wp) then
            cfl = v*dt_adv*idyT(i, j - 1)
            h_face = hl_y_work(i, j)
            dh_d = hr_y_work(i, j - 1) - h_face
            curv3_d = hr_y_work(i, j - 1) + h_face - 2.0_wp*htot_work(i, j - 1)
            vhtot_work(i, j) = dx_cv(i, j)*v*volcfl_face(h_face, dh_d, curv3_d, cfl)
         else if (v < 0.0_wp) then
            cfl = (-v)*dt_adv*idyT(i, j)
            h_face = hr_y_work(i, j)
            dh_d = hl_y_work(i, j + 1) - h_face
            curv3_d = h_face + hl_y_work(i, j + 1) - 2.0_wp*htot_work(i, j)
            vhtot_work(i, j) = dx_cv(i, j)*v*volcfl_face(h_face, dh_d, curv3_d, cfl)
         else
            vhtot_work(i, j) = 0.0_wp
         end if
      end do
      do concurrent(i=1:nx)
         vhtot_work(i, 1) = 0.0_wp
         vhtot_work(i, ny + 1) = 0.0_wp
      end do
      do concurrent(i=1:nx)
         vhtot_work(i, nghost + 1) = 0.0_wp
         vhtot_work(i, nghost + ny_phys + 1) = 0.0_wp
      end do

      do concurrent(j=1:ny + 1, i=1:nx, c=1:ncat) local(i_htot)
         if (vhtot_work(i, j) == 0.0_wp) then
            vh_out(i, j, c) = 0.0_wp
         else if (v_ice(i, j) >= 0.0_wp) then
            if (htot_work(i, j - 1) > 0.0_wp) then
               i_htot = 1.0_wp/htot_work(i, j - 1)
            else
               i_htot = 0.0_wp
            end if
            vh_out(i, j, c) = vhtot_work(i, j)*mca(i, j - 1, c)*i_htot
         else
            if (htot_work(i, j) > 0.0_wp) then
               i_htot = 1.0_wp/htot_work(i, j)
            else
               i_htot = 0.0_wp
            end if
            vh_out(i, j, c) = vhtot_work(i, j)*mca(i, j, c)*i_htot
         end if
      end do
   end subroutine ice_cat_flux_y_impl

   ! ----------------------------------------------------------------------
   ! (b) Snow masked by the co-located ice flux
   ! ----------------------------------------------------------------------

   pure subroutine ice_mask_snow_by_ice_impl(uh_ice, uh_snow, nfi, nfj, ncat)
      !! SIS2 `masking_uh=uh_ice` (`SIS_continuity.F90:1210-1215`):
      !! `uh_snow(I,c) = 0` wherever `uh_ice(I,c) == 0`.  Face-shaped
      !! array, works identically for the x-face `(nx+1,ny,ncat)` and
      !! y-face `(nx,ny+1,ncat)` layouts (caller passes the right
      !! extents).  `frac_neglect` (SIS2's second masking clause) is dead
      !! at SIS2's own default `frac_neglect=0` — not ported (D-noted in
      !! the module docstring).
      integer, intent(in) :: nfi, nfj, ncat
      real(wp), intent(in) :: uh_ice(nfi, nfj, ncat)
      real(wp), intent(inout) :: uh_snow(nfi, nfj, ncat)
      integer :: i, j, c

      do concurrent(j=1:nfj, i=1:nfi, c=1:ncat)
         if (uh_ice(i, j, c) == 0.0_wp) uh_snow(i, j, c) = 0.0_wp
      end do
   end subroutine ice_mask_snow_by_ice_impl

   ! ----------------------------------------------------------------------
   ! (c) PCM tracer riding: gather (face) then cell update
   ! ----------------------------------------------------------------------

   pure subroutine ice_gather_flux_x_impl(uh, val, tr_flux_x_work, ncat, nx, ny)
      !! Gather pass (race-free): `tr_flux_x_work(I,c) = val(donor of
      !! I,c)`, donor by the sign of the FLUX.  `uh` is exactly 0.0 at
      !! the array-edge faces `I=1` and `I=nx+1` (the flux kernel zeros
      !! them unconditionally), so the `I=1` donor-by-sign branch
      !! (`uh>=0`) would read the out-of-bounds `val(0,j,c)` — guarded
      !! explicitly below (`I==1`/`I==nx+1` fall back to the IN-BOUNDS
      !! neighbour; the value is never actually consumed by the update
      !! kernel there since `uh==0` at both those faces makes them
      !! inert, but the gather must still avoid the invalid index).
      !! Reads `val` ONLY at a donor cell — never writes `val` — so this
      !! kernel has no race with any other iteration.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: uh(nx + 1, ny, ncat)
      real(wp), intent(in) :: val(nx, ny, ncat)
      real(wp), intent(out) :: tr_flux_x_work(nx + 1, ny, ncat)
      integer :: i, j, c

      do concurrent(j=1:ny, i=1:nx + 1, c=1:ncat)
         if (i == 1) then
            tr_flux_x_work(i, j, c) = val(1, j, c)
         else if (i == nx + 1) then
            tr_flux_x_work(i, j, c) = val(nx, j, c)
         else if (uh(i, j, c) >= 0.0_wp) then
            tr_flux_x_work(i, j, c) = val(i - 1, j, c)
         else
            tr_flux_x_work(i, j, c) = val(i, j, c)
         end if
      end do
   end subroutine ice_gather_flux_x_impl

   pure subroutine ice_gather_flux_x_layer_impl(uh, val4, layer, tr_flux_x_work, ncat, nk, nx, ny)
      !! Layer-indexed twin of `ice_gather_flux_x_impl` for
      !! `enth_ice`/`sal_ice`/`enth_snow` (shape `(nx,ny,ncat,nk)`).
      integer, intent(in) :: layer, ncat, nk, nx, ny
      real(wp), intent(in) :: uh(nx + 1, ny, ncat)
      real(wp), intent(in) :: val4(nx, ny, ncat, nk)
      real(wp), intent(out) :: tr_flux_x_work(nx + 1, ny, ncat)
      integer :: i, j, c

      do concurrent(j=1:ny, i=1:nx + 1, c=1:ncat)
         if (i == 1) then
            tr_flux_x_work(i, j, c) = val4(1, j, c, layer)
         else if (i == nx + 1) then
            tr_flux_x_work(i, j, c) = val4(nx, j, c, layer)
         else if (uh(i, j, c) >= 0.0_wp) then
            tr_flux_x_work(i, j, c) = val4(i - 1, j, c, layer)
         else
            tr_flux_x_work(i, j, c) = val4(i, j, c, layer)
         end if
      end do
   end subroutine ice_gather_flux_x_layer_impl

   pure subroutine ice_ride_update_x_impl(iareaT, uh, tr_flux_x_work, mca, val, dt_adv, &
                                          ncat, nx, ny)
      !! Cell-update pass (race-free): reads `tr_flux_x_work` (the
      !! GATHERED donor VALUE at each face, from `ice_gather_flux_x_impl`)
      !! + its OWN cell's `mca`/`val` (never a neighbour's `val`), applies
      !! SIS2's `advect_scalar_x` flux-form update + the `H_NEGLECT`
      !! conditioning guard (`SIS_tracer_advect.F90:735-760`), writes
      !! `val(i,j,c)`.  `hnew` is the SAME expression the mass-update
      !! kernel (d) uses, so the implied masses agree bitwise.
      !! Bitwise no-op when both faces carry zero flux (`F_W==F_E==0`).
      !!
      !! Per-area transcription of SIS2's cell-integrated algebra
      !! (SIS2 works in `hprev = mca_old*areaT`, `uhh = flux*dt`; dividing
      !! every SIS2 quantity by `areaT` gives the per-area form below —
      !! `dtI = dt_adv*iareaT` plays the role of SIS2's `dt/areaT`):
      !!   `hlst = max(mca_old, 0)`; `hnew = mca_old - dtI*(fe-fw)`.
      !!   `hnew <= 0`            : hold (mass gone).
      !!   `0 < hnew < H_NEGLECT`  : `h_add = H_NEGLECT - hnew`;
      !!     `I_htot = 1/(hlst + dtI*(|fe|+|fw|))` (0 if the denominator
      !!     is 0); `hlst_adj = hlst + h_add*hlst*I_htot`;
      !!     `haddE = h_add*dtI*|fe|*I_htot`, `haddW = h_add*dtI*|fw|*I_htot`
      !!     (SIS2's sign convention: `haddE` is SUBTRACTED from the east
      !!     outflow term, `haddW` is ADDED to the west inflow term — both
      !!     push the effective in/outflow toward "more mass stays");
      !!     `val_new = (val*hlst_adj - ((fe*dtI-haddE)*val_e -
      !!     (fw*dtI+haddW)*val_w)) / H_NEGLECT`.
      !!   `hnew >= H_NEGLECT`     : plain flux form,
      !!     `val_new = (val*mca_old - dtI*(fe*val_e-fw*val_w))/hnew`.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: uh(nx + 1, ny, ncat)
      real(wp), intent(in) :: tr_flux_x_work(nx + 1, ny, ncat)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: val(nx, ny, ncat)
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c
      real(wp) :: fw, fe, val_w, val_e, mca_old, hlst, dti, hnew, h_add, denom, i_htot
      real(wp) :: hlst_adj, haddw, hadde, fw_term, fe_term

      do concurrent(j=1:ny, i=1:nx, c=1:ncat) &
         local(fw, fe, val_w, val_e, mca_old, hlst, dti, hnew, h_add, denom, i_htot, &
               hlst_adj, haddw, hadde, fw_term, fe_term)
         fw = uh(i, j, c)
         fe = uh(i + 1, j, c)
         if (fw /= 0.0_wp .or. fe /= 0.0_wp) then
            val_w = tr_flux_x_work(i, j, c)
            val_e = tr_flux_x_work(i + 1, j, c)
            mca_old = mca(i, j, c)
            hlst = max(mca_old, 0.0_wp)
            dti = dt_adv*iareaT(i, j)
            hnew = mca_old - dti*(fe - fw)
            if (hnew <= 0.0_wp) then
               ! mass gone; hold val (inert — PR 4a massless convention)
               continue
            else if (hnew < H_NEGLECT_ICE_TRANSPORT) then
               h_add = H_NEGLECT_ICE_TRANSPORT - hnew
               denom = hlst + dti*(abs(fe) + abs(fw))
               if (denom > 0.0_wp) then
                  i_htot = 1.0_wp/denom
               else
                  i_htot = 0.0_wp
               end if
               hlst_adj = hlst + h_add*hlst*i_htot
               haddw = h_add*dti*abs(fw)*i_htot
               hadde = h_add*dti*abs(fe)*i_htot
               fe_term = fe*dti - hadde
               fw_term = fw*dti + haddw
               val(i, j, c) = (val(i, j, c)*hlst_adj - (fe_term*val_e - fw_term*val_w)) &
                              /H_NEGLECT_ICE_TRANSPORT
            else
               val(i, j, c) = (val(i, j, c)*mca_old - dti*(fe*val_e - fw*val_w))/hnew
            end if
         end if
      end do
   end subroutine ice_ride_update_x_impl

   pure subroutine ice_ride_update_x_layer_impl(iareaT, uh, tr_flux_x_work, mca, val4, layer, &
                                                dt_adv, ncat, nk, nx, ny)
      !! Layer-indexed twin of `ice_ride_update_x_impl` for
      !! `enth_ice`/`sal_ice`/`enth_snow` (shape `(nx,ny,ncat,nk)`).  Same
      !! algebra, applied to `val4(:,:,:,layer)`.
      integer, intent(in) :: ncat, nk, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: uh(nx + 1, ny, ncat)
      real(wp), intent(in) :: tr_flux_x_work(nx + 1, ny, ncat)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: val4(nx, ny, ncat, nk)
      integer, intent(in) :: layer
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c
      real(wp) :: fw, fe, val_w, val_e, mca_old, hlst, dti, hnew, h_add, denom, i_htot
      real(wp) :: hlst_adj, haddw, hadde, fw_term, fe_term

      do concurrent(j=1:ny, i=1:nx, c=1:ncat) &
         local(fw, fe, val_w, val_e, mca_old, hlst, dti, hnew, h_add, denom, i_htot, &
               hlst_adj, haddw, hadde, fw_term, fe_term)
         fw = uh(i, j, c)
         fe = uh(i + 1, j, c)
         if (fw /= 0.0_wp .or. fe /= 0.0_wp) then
            val_w = tr_flux_x_work(i, j, c)
            val_e = tr_flux_x_work(i + 1, j, c)
            mca_old = mca(i, j, c)
            hlst = max(mca_old, 0.0_wp)
            dti = dt_adv*iareaT(i, j)
            hnew = mca_old - dti*(fe - fw)
            if (hnew <= 0.0_wp) then
               continue
            else if (hnew < H_NEGLECT_ICE_TRANSPORT) then
               h_add = H_NEGLECT_ICE_TRANSPORT - hnew
               denom = hlst + dti*(abs(fe) + abs(fw))
               if (denom > 0.0_wp) then
                  i_htot = 1.0_wp/denom
               else
                  i_htot = 0.0_wp
               end if
               hlst_adj = hlst + h_add*hlst*i_htot
               haddw = h_add*dti*abs(fw)*i_htot
               hadde = h_add*dti*abs(fe)*i_htot
               fe_term = fe*dti - hadde
               fw_term = fw*dti + haddw
               val4(i, j, c, layer) = (val4(i, j, c, layer)*hlst_adj - (fe_term*val_e - fw_term*val_w)) &
                                      /H_NEGLECT_ICE_TRANSPORT
            else
               val4(i, j, c, layer) = (val4(i, j, c, layer)*mca_old - dti*(fe*val_e - fw*val_w))/hnew
            end if
         end if
      end do
   end subroutine ice_ride_update_x_layer_impl

   pure subroutine ice_gather_flux_y_impl(vh, val, tr_flux_y_work, ncat, nx, ny)
      !! Meridional twin of `ice_gather_flux_x_impl`.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: vh(nx, ny + 1, ncat)
      real(wp), intent(in) :: val(nx, ny, ncat)
      real(wp), intent(out) :: tr_flux_y_work(nx, ny + 1, ncat)
      integer :: i, j, c

      do concurrent(j=1:ny + 1, i=1:nx, c=1:ncat)
         if (j == 1) then
            tr_flux_y_work(i, j, c) = val(i, 1, c)
         else if (j == ny + 1) then
            tr_flux_y_work(i, j, c) = val(i, ny, c)
         else if (vh(i, j, c) >= 0.0_wp) then
            tr_flux_y_work(i, j, c) = val(i, j - 1, c)
         else
            tr_flux_y_work(i, j, c) = val(i, j, c)
         end if
      end do
   end subroutine ice_gather_flux_y_impl

   pure subroutine ice_gather_flux_y_layer_impl(vh, val4, layer, tr_flux_y_work, ncat, nk, nx, ny)
      !! Meridional twin of `ice_gather_flux_x_layer_impl`.
      integer, intent(in) :: layer, ncat, nk, nx, ny
      real(wp), intent(in) :: vh(nx, ny + 1, ncat)
      real(wp), intent(in) :: val4(nx, ny, ncat, nk)
      real(wp), intent(out) :: tr_flux_y_work(nx, ny + 1, ncat)
      integer :: i, j, c

      do concurrent(j=1:ny + 1, i=1:nx, c=1:ncat)
         if (j == 1) then
            tr_flux_y_work(i, j, c) = val4(i, 1, c, layer)
         else if (j == ny + 1) then
            tr_flux_y_work(i, j, c) = val4(i, ny, c, layer)
         else if (vh(i, j, c) >= 0.0_wp) then
            tr_flux_y_work(i, j, c) = val4(i, j - 1, c, layer)
         else
            tr_flux_y_work(i, j, c) = val4(i, j, c, layer)
         end if
      end do
   end subroutine ice_gather_flux_y_layer_impl

   pure subroutine ice_ride_update_y_impl(iareaT, vh, tr_flux_y_work, mca, val, dt_adv, &
                                          ncat, nx, ny)
      !! Meridional twin of `ice_ride_update_x_impl`.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: vh(nx, ny + 1, ncat)
      real(wp), intent(in) :: tr_flux_y_work(nx, ny + 1, ncat)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: val(nx, ny, ncat)
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c
      real(wp) :: fs, fn, val_s, val_n, mca_old, hlst, dti, hnew, h_add, denom, i_htot
      real(wp) :: hlst_adj, hadds, haddn, fs_term, fn_term

      do concurrent(j=1:ny, i=1:nx, c=1:ncat) &
         local(fs, fn, val_s, val_n, mca_old, hlst, dti, hnew, h_add, denom, i_htot, &
               hlst_adj, hadds, haddn, fs_term, fn_term)
         fs = vh(i, j, c)
         fn = vh(i, j + 1, c)
         if (fs /= 0.0_wp .or. fn /= 0.0_wp) then
            val_s = tr_flux_y_work(i, j, c)
            val_n = tr_flux_y_work(i, j + 1, c)
            mca_old = mca(i, j, c)
            hlst = max(mca_old, 0.0_wp)
            dti = dt_adv*iareaT(i, j)
            hnew = mca_old - dti*(fn - fs)
            if (hnew <= 0.0_wp) then
               continue
            else if (hnew < H_NEGLECT_ICE_TRANSPORT) then
               h_add = H_NEGLECT_ICE_TRANSPORT - hnew
               denom = hlst + dti*(abs(fn) + abs(fs))
               if (denom > 0.0_wp) then
                  i_htot = 1.0_wp/denom
               else
                  i_htot = 0.0_wp
               end if
               hlst_adj = hlst + h_add*hlst*i_htot
               hadds = h_add*dti*abs(fs)*i_htot
               haddn = h_add*dti*abs(fn)*i_htot
               fn_term = fn*dti - haddn
               fs_term = fs*dti + hadds
               val(i, j, c) = (val(i, j, c)*hlst_adj - (fn_term*val_n - fs_term*val_s)) &
                              /H_NEGLECT_ICE_TRANSPORT
            else
               val(i, j, c) = (val(i, j, c)*mca_old - dti*(fn*val_n - fs*val_s))/hnew
            end if
         end if
      end do
   end subroutine ice_ride_update_y_impl

   pure subroutine ice_ride_update_y_layer_impl(iareaT, vh, tr_flux_y_work, mca, val4, layer, &
                                                dt_adv, ncat, nk, nx, ny)
      !! Layer-indexed twin of `ice_ride_update_y_impl` for
      !! `enth_ice`/`sal_ice`/`enth_snow` (shape `(nx,ny,ncat,nk)`).
      integer, intent(in) :: ncat, nk, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: vh(nx, ny + 1, ncat)
      real(wp), intent(in) :: tr_flux_y_work(nx, ny + 1, ncat)
      real(wp), intent(in) :: mca(nx, ny, ncat)
      real(wp), intent(inout) :: val4(nx, ny, ncat, nk)
      integer, intent(in) :: layer
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c
      real(wp) :: fs, fn, val_s, val_n, mca_old, hlst, dti, hnew, h_add, denom, i_htot
      real(wp) :: hlst_adj, hadds, haddn, fs_term, fn_term

      do concurrent(j=1:ny, i=1:nx, c=1:ncat) &
         local(fs, fn, val_s, val_n, mca_old, hlst, dti, hnew, h_add, denom, i_htot, &
               hlst_adj, hadds, haddn, fs_term, fn_term)
         fs = vh(i, j, c)
         fn = vh(i, j + 1, c)
         if (fs /= 0.0_wp .or. fn /= 0.0_wp) then
            val_s = tr_flux_y_work(i, j, c)
            val_n = tr_flux_y_work(i, j + 1, c)
            mca_old = mca(i, j, c)
            hlst = max(mca_old, 0.0_wp)
            dti = dt_adv*iareaT(i, j)
            hnew = mca_old - dti*(fn - fs)
            if (hnew <= 0.0_wp) then
               continue
            else if (hnew < H_NEGLECT_ICE_TRANSPORT) then
               h_add = H_NEGLECT_ICE_TRANSPORT - hnew
               denom = hlst + dti*(abs(fn) + abs(fs))
               if (denom > 0.0_wp) then
                  i_htot = 1.0_wp/denom
               else
                  i_htot = 0.0_wp
               end if
               hlst_adj = hlst + h_add*hlst*i_htot
               hadds = h_add*dti*abs(fs)*i_htot
               haddn = h_add*dti*abs(fn)*i_htot
               fn_term = fn*dti - haddn
               fs_term = fs*dti + hadds
               val4(i, j, c, layer) = (val4(i, j, c, layer)*hlst_adj - (fn_term*val_n - fs_term*val_s)) &
                                      /H_NEGLECT_ICE_TRANSPORT
            else
               val4(i, j, c, layer) = (val4(i, j, c, layer)*mca_old - dti*(fn*val_n - fs*val_s))/hnew
            end if
         end if
      end do
   end subroutine ice_ride_update_y_layer_impl

   ! ----------------------------------------------------------------------
   ! (d) Mass update
   ! ----------------------------------------------------------------------

   pure subroutine ice_mass_update_x_impl(iareaT, uh, mca, dt_adv, nghost, nx_phys, ncat, nx, ny)
      !! SIS2 `:183-191`: `mca(c) -= dt_adv*iareaT*(uh(I)-uh(I-1))`,
      !! physical cells only.  No race (reads faces, writes cells).
      integer, intent(in) :: nghost, nx_phys, ncat, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: uh(nx + 1, ny, ncat)
      real(wp), intent(inout) :: mca(nx, ny, ncat)
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c, i_lo, i_hi

      i_lo = nghost + 1
      i_hi = nghost + nx_phys

      do concurrent(j=1:ny, i=i_lo:i_hi, c=1:ncat)
         mca(i, j, c) = mca(i, j, c) - dt_adv*iareaT(i, j)*(uh(i + 1, j, c) - uh(i, j, c))
      end do
   end subroutine ice_mass_update_x_impl

   pure subroutine ice_mass_update_y_impl(iareaT, vh, mca, dt_adv, nghost, ny_phys, ncat, nx, ny)
      !! Meridional twin of `ice_mass_update_x_impl`.
      integer, intent(in) :: nghost, ny_phys, ncat, nx, ny
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: vh(nx, ny + 1, ncat)
      real(wp), intent(inout) :: mca(nx, ny, ncat)
      real(wp), intent(in) :: dt_adv
      integer :: i, j, c, j_lo, j_hi

      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      do concurrent(j=j_lo:j_hi, i=1:nx, c=1:ncat)
         mca(i, j, c) = mca(i, j, c) - dt_adv*iareaT(i, j)*(vh(i, j + 1, c) - vh(i, j, c))
      end do
   end subroutine ice_mass_update_y_impl

   ! ----------------------------------------------------------------------
   ! Post-pass validity reduction (D7): fail-loud negativity/orphan-snow
   ! detection, GPU-safe (device reduction instead of an in-loop FATAL).
   ! ----------------------------------------------------------------------

   pure subroutine ice_validity_reduce_impl(mca_ice, mca_snow, nghost, nx_phys, ny_phys, &
                                            ncat, nx, ny, ok)
      !! `ok = .false.` when `min(mca_ice) < 0`, `min(mca_snow) < 0`, or
      !! an "orphan snow" cell exists (`mca_snow > 0` where `mca_ice <= 0`
      !! by more than `H_NEGLECT_ICE_TRANSPORT`) — SIS2 FATALs on any of
      !! these; ported as a post-pass reduction (GPU-safe fail-loud, D7)
      !! rather than an in-loop abort.  Physical cells only.
      integer, intent(in) :: nghost, nx_phys, ny_phys, ncat, nx, ny
      real(wp), intent(in) :: mca_ice(nx, ny, ncat)
      real(wp), intent(in) :: mca_snow(nx, ny, ncat)
      logical, intent(out) :: ok
      integer :: i, j, c, i_lo, i_hi, j_lo, j_hi
      real(wp) :: min_ice, min_snow, max_orphan

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      min_ice = 0.0_wp
      min_snow = 0.0_wp
      max_orphan = 0.0_wp
      do concurrent(c=1:ncat, j=j_lo:j_hi, i=i_lo:i_hi) reduce(min:min_ice, min_snow)
         min_ice = min(min_ice, mca_ice(i, j, c))
         min_snow = min(min_snow, mca_snow(i, j, c))
      end do
      do concurrent(c=1:ncat, j=j_lo:j_hi, i=i_lo:i_hi) reduce(max:max_orphan)
         if (mca_ice(i, j, c) <= 0.0_wp) then
            max_orphan = max(max_orphan, mca_snow(i, j, c))
         end if
      end do

      ok = (min_ice >= 0.0_wp) .and. (min_snow >= 0.0_wp) &
           .and. (max_orphan <= H_NEGLECT_ICE_TRANSPORT)
   end subroutine ice_validity_reduce_impl

   ! ======================================================================
   ! Phase 3: CAS -> IST
   ! ======================================================================

   pure subroutine ice_cas_to_ist_impl(wet_mask, areaT, mca_ice, mca_snow, part_size, m_ice, &
                                       m_snow, mh_lim, roll_factor, nghost, ncat, nx, ny)
      !! SIS2 `cell_ave_state_to_ice_state` (`:540`).  Per cell, per
      !! category: pre-floor category 1, optional thin-ice rolling
      !! (`roll_factor > 0`), general floor, then re-derive
      !! `part_size(c) = mca_ice(c)/m_ice(c)` and
      !! `m_snow(c) = m_ice(c)*(mca_snow(c)/mca_ice(c))` (per-ICE-area).
      !! `part_size(0) = 1 - Sum_c part_size(c)` — MAY be negative here;
      !! `ice_compress_impl` (Phase 4) is what restores >= 0.
      !!
      !! `roll_factor` rolling test — VERIFIED against the literal SIS2
      !! source (`SIS_transport.F90:572-581`, `L_to_H = US%L_to_Z*Rho_ice`;
      !! Roundabout carries no L_to_Z length rescaling, so `L_to_H =
      !! ICE_RHO_ICE`): a VOLUME-vs-VOLUME comparison, NOT a per-area one
      !! — `areaT` does NOT cancel (only the RHS carries it explicitly in
      !! SIS2's actual, non-simplified code):
      !!   `h_eff = m_ice(c)/ICE_RHO_ICE` (pre-roll per-ICE-area thickness,
      !!   a length);
      !!   `if (roll_factor*h_eff**3 > (mca_ice(c)/ICE_RHO_ICE)*areaT(i,j))`
      !!   — LHS a volume (h_eff cubed), RHS `(mca_ice/ICE_RHO_ICE)` is a
      !!   per-CELL-area ice-volume-per-area (a length) times `areaT` =
      !!   the cell's total ice volume;
      !!   then `m_ice(c) = max(mh_lim(1), ICE_RHO_ICE*sqrt((mca_ice(c)*
      !!   areaT(i,j)/ICE_RHO_ICE) / (roll_factor*h_eff)))` — SIS2's sqrt
      !!   argument is `(m_ice*areaT)/(roll_factor*mH_ice)` in ITS
      !!   pre-division units; converting consistently (dividing the
      !!   `mca_ice*areaT` term by `ICE_RHO_ICE` once, matching the
      !!   `(mca_ice/ICE_RHO_ICE)*areaT` volume on the LHS test) gives the
      !!   expression coded below.
      integer, intent(in) :: nghost, ncat, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: mca_ice(nx, ny, ncat)
      real(wp), intent(in) :: mca_snow(nx, ny, ncat)
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(in) :: mh_lim(ncat + 1)
      real(wp), intent(in) :: roll_factor
      integer :: i, j, c, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h_eff, ice_vol, part_sum

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi, c=1:ncat) local(h_eff, ice_vol)
         if (wet_mask(i, j) > 0.5_wp) then
            if (c == 1 .and. mca_ice(i, j, 1) > 0.0_wp .and. m_ice(i, j, 1) < mh_lim(1)) then
               m_ice(i, j, 1) = mh_lim(1)
            end if
            if (mca_ice(i, j, c) > 0.0_wp) then
               if (roll_factor > 0.0_wp) then
                  h_eff = m_ice(i, j, c)/ICE_RHO_ICE
                  ice_vol = (mca_ice(i, j, c)/ICE_RHO_ICE)*areaT(i, j)
                  if (roll_factor*h_eff**3 > ice_vol) then
                     m_ice(i, j, c) = max(mh_lim(1), &
                                          ICE_RHO_ICE*sqrt(ice_vol/(roll_factor*h_eff)))
                  end if
               end if
               if (m_ice(i, j, c) < mh_lim(1)) m_ice(i, j, c) = mh_lim(1)
               part_size(i, j, c) = mca_ice(i, j, c)/m_ice(i, j, c)
               m_snow(i, j, c) = m_ice(i, j, c)*(mca_snow(i, j, c)/mca_ice(i, j, c))
            else
               part_size(i, j, c) = 0.0_wp
               m_ice(i, j, c) = 0.0_wp
               m_snow(i, j, c) = 0.0_wp
               ! enth/sal held (inert — massless convention)
            end if
         end if
      end do
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(part_sum, c)
         if (wet_mask(i, j) > 0.5_wp) then
            part_sum = 0.0_wp
            do c = 1, ncat
               part_sum = part_sum + part_size(i, j, c)
            end do
            part_size(i, j, 0) = 1.0_wp - part_sum
         end if
      end do
   end subroutine ice_cas_to_ist_impl

   ! ======================================================================
   ! Phase 4: compress_ice
   ! ======================================================================

   pure subroutine ice_compress_impl(wet_mask, part_size, m_ice, m_snow, enth_ice, enth_snow, &
                                     sal_ice, mh_lim, nghost, ncat, nk, nx, ny, ok)
      !! Outer per-cell dispatch: ONE `do concurrent(j,i)` over physical
      !! cells (ghosts excluded), `wet_mask > 0.5` inner gate, serial in
      !! category within the cell (same shape as
      !! `ice_adjust_categories_impl`).  The per-cell algebra is the SAME
      !! algorithm `ice_transport_compress_cell` implements (public,
      !! directly unit-testable on small standalone arrays — SPEC §7 test
      !! 3's single-cell hand-check) — but this production kernel does
      !! NOT call it with a derived slice: `part_size(i,j,:)` etc. are
      !! NON-CONTIGUOUS sections (the category axis is not the fastest
      !! dimension), which a device `!$acc routine seq` call cannot take
      !! safely (would force a compiler temporary — the array-of-
      !! derived-type/strided-slice indirection trap, CLAUDE.md memory).
      !! Instead the algorithm is INLINED here operating on `(i,j,c)`
      !! triples directly, mirroring `ice_adjust_categories_impl`'s
      !! established pattern.  `ok` is a per-cell-then-reduced flag
      !! (D7: SIS2 FATALs on the top-category overflow inconsistency;
      !! ported as a `!$acc parallel loop reduction` instead, since
      !! `do concurrent` cannot itself carry a boolean/min reduction in
      !! the house style — CLAUDE.md "Reductions use `!$acc parallel loop
      !! reduction(...)`").
      integer, intent(in) :: nghost, ncat, nk, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: enth_snow(nx, ny, ncat, 1)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
      real(wp), intent(in) :: mh_lim(ncat + 1)
      logical, intent(out) :: ok
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: ok_acc
      logical :: cell_ok

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      ok_acc = 1.0_wp
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(cell_ok) reduce(min:ok_acc)
         if (wet_mask(i, j) > 0.5_wp) then
            call ice_compress_cell_inline(part_size, m_ice, m_snow, enth_ice, &
                                          enth_snow, sal_ice, mh_lim, i, j, ncat, &
                                          nk, nx, ny, cell_ok)
            if (.not. cell_ok) ok_acc = min(ok_acc, 0.0_wp)
         end if
      end do
      ok = (ok_acc > 0.5_wp)
   end subroutine ice_compress_impl

   pure subroutine ice_compress_cell_inline(part_size, m_ice, m_snow, enth_ice, enth_snow, &
                                            sal_ice, mh_lim, i, j, ncat, nk, nx, ny, ok)
      !$acc routine seq
      !! KEEP IN SYNC with `ice_transport_compress_cell` (the HOST tested
      !! seam twin, directly below). Same excess/ratio/compaction algorithm;
      !! this twin exists only for a different argument shape — full
      !! device-present arrays + scalar `(i,j)` indices for the fused device
      !! path (fixed-size device locals), vs the twin's small standalone
      !! per-cell arrays for the unit test. Any change to the compaction
      !! logic here MUST be mirrored there (the test only drives the twin).
      !!
      !! Device-safe per-cell compress body: takes the FULL state arrays
      !! plus scalar `(i,j)` indices (no array-section slicing — every
      !! access below is a direct `(i,j,c[,l])` element read/write, so
      !! this is safe to call from `!$acc parallel loop`/`do concurrent`
      !! with the full arrays already device-present).  Same algorithm as
      !! `ice_transport_compress_cell` (see that routine's docstring for
      !! the algorithm narrative); duplicated rather than shared because
      !! the two need different argument shapes (device-safe
      !! full-array-plus-index here; small standalone per-cell arrays
      !! there for the unit test) — CLAUDE.md "duplicate explicitly" (no
      !! include-style body sharing).  Conservation is exact and the two
      !! are SIS2-roundoff-equivalent, but NOT bitwise identical across
      !! CHAINED cascades: this device body recomputes `mca_this =
      !! part_size(c)*m_ice(c)` fresh per iteration, so after a c-1 -> c
      !! transfer (which just re-derived `m_ice(c)` by division) its
      !! `mca_this` differs at round-off from the standalone's running
      !! `mca_c(c)` array (and from SIS2's own running array) — the
      !! transferred/removed masses still telescope to the same total, so
      !! `Σ_c part*m` is conserved regardless.
      integer, intent(in) :: i, j, ncat, nk, nx, ny
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: enth_snow(nx, ny, ncat, 1)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
      real(wp), intent(in) :: mh_lim(ncat + 1)
      logical, intent(out) :: ok

      integer :: c, l
      real(wp) :: excess, ratio, f, part_sum
      real(wp) :: mca_this, msnow_this, mca_next, msnow_next
      real(wp) :: mca_old, trans, mnew, msnew

      ok = .true.
      if (part_size(i, j, 0) >= 0.0_wp) return

      excess = -part_size(i, j, 0)
      part_size(i, j, 0) = 0.0_wp

      ! No pre-materialized per-category mca/msnow arrays (ncat is not a
      ! compile-time bound — CLAUDE.md fixed-size-locals rule for
      ! `!$acc routine seq`).  Each iteration recomputes `mca_this`
      ! (category c) and `mca_next` (category c+1) FRESH from
      ! `part_size*m_ice`, BEFORE this iteration's own mutations — since
      ! the cascade only ever reads/writes categories c and c+1 at step
      ! c, and category c+1 has not yet been touched by any earlier
      ! iteration, this reproduces the array-based algorithm exactly
      ! (verified against the Python prototype, both branches, to full
      ! precision).
      do c = 1, ncat - 1
         mca_this = part_size(i, j, c)*m_ice(i, j, c)
         if (excess > 0.0_wp .and. mca_this > 0.0_wp) then
            msnow_this = part_size(i, j, c)*m_snow(i, j, c)
            ratio = m_ice(i, j, c)/mh_lim(c + 1)
            if (part_size(i, j, c)*(1.0_wp - ratio) >= excess) then
               f = part_size(i, j, c)/(part_size(i, j, c) - excess)
               m_ice(i, j, c) = m_ice(i, j, c)*f
               m_snow(i, j, c) = m_snow(i, j, c)*f
               part_size(i, j, c) = part_size(i, j, c) - excess
               excess = 0.0_wp
            else
               excess = excess - part_size(i, j, c)*(1.0_wp - ratio)
               if (mca_this > MASS_NEGLECT_ICE_TRANSPORT) then
                  mca_next = part_size(i, j, c + 1)*m_ice(i, j, c + 1)
                  msnow_next = part_size(i, j, c + 1)*m_snow(i, j, c + 1)
                  part_size(i, j, c + 1) = part_size(i, j, c + 1) + part_size(i, j, c)*ratio
                  mca_old = mca_next
                  trans = mca_this
                  mnew = mca_next + mca_this
                  if (part_size(i, j, c + 1) > MASS_NEGLECT_ICE_TRANSPORT) then
                     m_ice(i, j, c + 1) = mnew/part_size(i, j, c + 1)
                  else if (trans > mca_old) then
                     part_size(i, j, c + 1) = mnew/mh_lim(c + 1)
                     m_ice(i, j, c + 1) = mh_lim(c + 1)
                  else
                     part_size(i, j, c + 1) = mnew/m_ice(i, j, c + 1)
                  end if
                  msnew = msnow_next + msnow_this
                  if (part_size(i, j, c + 1) > 0.0_wp) then
                     m_snow(i, j, c + 1) = msnew/part_size(i, j, c + 1)
                  else
                     m_snow(i, j, c + 1) = 0.0_wp
                  end if
                  if (trans > 0.0_wp .and. mnew > 0.0_wp) then
                     do l = 1, nk
                        enth_ice(i, j, c + 1, l) = (trans*enth_ice(i, j, c, l) &
                                                    + mca_old*enth_ice(i, j, c + 1, l))/mnew
                        sal_ice(i, j, c + 1, l) = (trans*sal_ice(i, j, c, l) &
                                                   + mca_old*sal_ice(i, j, c + 1, l))/mnew
                     end do
                  end if
                  if (msnow_this > 0.0_wp .and. msnew > 0.0_wp) then
                     enth_snow(i, j, c + 1, 1) = ((msnew - msnow_this)*enth_snow(i, j, c + 1, 1) &
                                                  + msnow_this*enth_snow(i, j, c, 1))/msnew
                  end if
               end if
               m_ice(i, j, c) = 0.0_wp
               m_snow(i, j, c) = 0.0_wp
               part_size(i, j, c) = 0.0_wp
            end if
         end if
      end do

      if (excess > 0.0_wp) then
         c = ncat
         ! Fail-loud (D7) when the top category cannot absorb the leftover
         ! excess: (a) SIS2's own consistency check (`part(ncat) <= 1 .and.
         ! excess > 2*ncat*eps`), PLUS (b) `part(ncat) <= excess` — the
         ! in-place compaction `f = part/(part-excess)` has a NEGATIVE (or,
         ! at the exact `part(ncat) == excess` tie, ZERO) denominator there
         ! and would flip `m_ice(ncat)` large-negative or divide by zero
         ! (SIS2 has the same hole; guarding it is strictly safer and the
         ! Phase-2 mca-only reduction cannot catch a compress-produced
         ! negative).
         if ((part_size(i, j, c) <= 1.0_wp .and. excess > 2.0_wp*real(ncat, wp)*epsilon(1.0_wp)) &
             .or. part_size(i, j, c) <= excess) then
            ok = .false.
            return
         end if
         f = part_size(i, j, c)/(part_size(i, j, c) - excess)
         m_ice(i, j, c) = m_ice(i, j, c)*f
         m_snow(i, j, c) = m_snow(i, j, c)*f
         part_size(i, j, c) = part_size(i, j, c) - excess
      end if

      part_sum = 0.0_wp
      do c = 1, ncat
         part_sum = part_sum + part_size(i, j, c)
      end do
      part_size(i, j, 0) = max(1.0_wp - part_sum, 0.0_wp)
   end subroutine ice_compress_cell_inline

   pure subroutine ice_transport_compress_cell(part_size, m_ice, m_snow, enth_ice, enth_snow, &
                                               sal_ice, mh_lim, ncat, nk, ok)
      !! KEEP IN SYNC with `ice_compress_cell_inline` (the DEVICE production
      !! twin, directly above). Same excess/ratio/compaction algorithm; this
      !! HOST twin exists only as the directly unit-testable seam (small
      !! standalone per-cell arrays, no `!$acc routine seq`), while the
      !! device twin takes full device-present arrays + a scalar `(i,j)`
      !! index. Any change to the compaction logic here MUST be mirrored
      !! there (the fused device path is what production runs).
      !!
      !! SIS2 `compress_ice` (`SIS_transport.F90:898-1100`), no-ridge
      !! default, ponds dropped — single-CELL kernel (public: directly
      !! unit-testable, SPEC §7 test 3).  `part_size(0:ncat)`,
      !! `m_ice`/`m_snow`/`enth_snow` `(ncat)`, `enth_ice`/`sal_ice`
      !! `(ncat, nk)`.  `part_size(0)` may be negative on entry (the
      !! open-water deficit from Phase 3); no-op (returns `ok=.true.`)
      !! when it is already `>= 0`.
      !!
      !! HOST-ONLY test entry point (no `!$acc routine seq` — this is
      !! deliberately NOT called from device kernels; a dummy-`ncat`-sized
      !! automatic array (`mca_c`/`msnow_c` below) is fine on the host but
      !! would violate the fixed-size-locals rule for a device routine).
      !! The PRODUCTION per-cell body is `ice_compress_cell_inline`
      !! (device-safe: full arrays + scalar `(i,j)` index, no per-category
      !! `mca_c`/`msnow_c` snapshot).  The two are CONSERVING +
      !! SIS2-roundoff-equivalent, but NOT bitwise identical across
      !! CHAINED cascades: this host body carries a pre-materialized
      !! `mca_c(ncat)` array (SIS2's own running array), whereas the
      !! device body recomputes `mca_this = part*m_ice` fresh per
      !! iteration — after a c-1 -> c transfer they differ at round-off,
      !! though both conserve `Σ_c part*m` exactly.  On a SINGLE transfer
      !! (the unit-test hand-check inputs) they ARE bit-identical and
      !! reproduce the Python prototype exactly.  Duplicated rather than
      !! shared, per CLAUDE.md "duplicate explicitly".
      !!
      !! Documented divergence (D8): the tracer merge is done INSIDE the
      !! thinnest-first cascade at each transfer site, whereas SIS2 defers
      !! it to `advect_tracers_thicker` AFTER the category k-loop
      !! (`SIS_tracer_advect.F90`) — equivalent to round-off since each
      !! boundary is visited once in a fixed order with running masses
      !! (same argument as the PR-4a `ice_adjust_categories` merge).
      !!
      !! Algorithm (thinnest-first, `c = 1..ncat-1`):
      !!   `excess = -part_size(0)`; `part_size(0) = 0`.
      !!   `mca_c(c) = part_size(c)*m_ice(c)`, `msnow_c(c) =
      !!   part_size(c)*m_snow(c)` (the CAS-absent branch, SIS2 `:986-992`
      !!   — recomputed here, not reused from Phase 2's `mca_ice`; SIS2
      !!   documents the same roundoff freedom, `:986`).
      !!   Per category `c` (while `excess > 0 .and. mca_c(c) > 0`):
      !!     `ratio = m_ice(c)/mh_lim(c+1)`.
      !!     If `part_size(c)*(1-ratio) >= excess`: IN-PLACE compaction —
      !!       `f = part_size(c)/(part_size(c)-excess)`;
      !!       `m_ice(c) *= f`; `m_snow(c) *= f`; `part_size(c) -= excess`;
      !!       `excess = 0`.
      !!     Else: `excess -= part_size(c)*(1-ratio)`; if
      !!       `mca_c(c) > MASS_NEGLECT_ICE_TRANSPORT`: transfer c -> c+1
      !!       (`part_size(c+1) += part_size(c)*ratio`; mass-weighted
      !!       merge of `mca_c`, `msnow_c`, mass-weighted tracer merge of
      !!       `enth_ice`/`sal_ice`/`enth_snow`; three-branch thickness
      !!       underflow guard on `m_ice(c+1)`/`part_size(c+1)`); zero
      !!       category `c`.
      !!   Top category (`c = ncat`): if `excess > 0`, consistency check
      !!   (`part_size(ncat) <= 1 .and. excess > 2*ncat*epsilon` => FATAL
      !!   in SIS2; here `ok = .false.`), else in-place compaction (same
      !!   formula, no transfer branch).
      !!   Cell epilogue: `part_size(0) = max(1 - Sum_c part_size(c), 0)`.
      integer, intent(in) :: ncat, nk
      real(wp), intent(inout) :: part_size(0:ncat)
      real(wp), intent(inout) :: m_ice(ncat)
      real(wp), intent(inout) :: m_snow(ncat)
      real(wp), intent(inout) :: enth_ice(ncat, nk)
      real(wp), intent(inout) :: enth_snow(ncat, 1)
      real(wp), intent(inout) :: sal_ice(ncat, nk)
      real(wp), intent(in) :: mh_lim(ncat + 1)
      logical, intent(out) :: ok

      integer :: c, l
      real(wp) :: excess, ratio, f, part_sum
      real(wp) :: mca_c(ncat), msnow_c(ncat)
      real(wp) :: mca_old, trans, mnew, msnew

      ok = .true.
      if (part_size(0) >= 0.0_wp) return

      excess = -part_size(0)
      part_size(0) = 0.0_wp
      do c = 1, ncat
         mca_c(c) = part_size(c)*m_ice(c)
         msnow_c(c) = part_size(c)*m_snow(c)
      end do

      do c = 1, ncat - 1
         if (excess > 0.0_wp .and. mca_c(c) > 0.0_wp) then
            ratio = m_ice(c)/mh_lim(c + 1)
            if (part_size(c)*(1.0_wp - ratio) >= excess) then
               f = part_size(c)/(part_size(c) - excess)
               m_ice(c) = m_ice(c)*f
               m_snow(c) = m_snow(c)*f
               part_size(c) = part_size(c) - excess
               excess = 0.0_wp
            else
               excess = excess - part_size(c)*(1.0_wp - ratio)
               if (mca_c(c) > MASS_NEGLECT_ICE_TRANSPORT) then
                  part_size(c + 1) = part_size(c + 1) + part_size(c)*ratio
                  mca_old = mca_c(c + 1)
                  trans = mca_c(c)
                  mca_c(c + 1) = mca_c(c + 1) + mca_c(c)
                  if (part_size(c + 1) > MASS_NEGLECT_ICE_TRANSPORT) then
                     m_ice(c + 1) = mca_c(c + 1)/part_size(c + 1)
                  else if (trans > mca_old) then
                     part_size(c + 1) = mca_c(c + 1)/mh_lim(c + 1)
                     m_ice(c + 1) = mh_lim(c + 1)
                  else
                     part_size(c + 1) = mca_c(c + 1)/m_ice(c + 1)
                  end if
                  msnow_c(c + 1) = msnow_c(c + 1) + msnow_c(c)
                  if (part_size(c + 1) > 0.0_wp) then
                     m_snow(c + 1) = msnow_c(c + 1)/part_size(c + 1)
                  else
                     m_snow(c + 1) = 0.0_wp
                  end if
                  mnew = trans + mca_old
                  if (trans > 0.0_wp .and. mnew > 0.0_wp) then
                     do l = 1, nk
                        enth_ice(c + 1, l) = (trans*enth_ice(c, l) + mca_old*enth_ice(c + 1, l))/mnew
                        sal_ice(c + 1, l) = (trans*sal_ice(c, l) + mca_old*sal_ice(c + 1, l))/mnew
                     end do
                  end if
                  msnew = msnow_c(c + 1)
                  if (msnow_c(c) > 0.0_wp .and. msnew > 0.0_wp) then
                     enth_snow(c + 1, 1) = ((msnew - msnow_c(c))*enth_snow(c + 1, 1) &
                                            + msnow_c(c)*enth_snow(c, 1))/msnew
                  end if
               end if
               mca_c(c) = 0.0_wp
               msnow_c(c) = 0.0_wp
               m_ice(c) = 0.0_wp
               m_snow(c) = 0.0_wp
               part_size(c) = 0.0_wp
            end if
         end if
      end do

      if (excess > 0.0_wp) then
         c = ncat
         ! Same top-category fail-loud as `ice_compress_cell_inline`,
         ! incl. the `part(ncat) <= excess` negative-or-zero-denominator guard
         ! (the `==` tie is the exact zero-denominator divide).
         if ((part_size(c) <= 1.0_wp .and. excess > 2.0_wp*real(ncat, wp)*epsilon(1.0_wp)) &
             .or. part_size(c) <= excess) then
            ok = .false.
            return
         end if
         f = part_size(c)/(part_size(c) - excess)
         m_ice(c) = m_ice(c)*f
         m_snow(c) = m_snow(c)*f
         part_size(c) = part_size(c) - excess
      end if

      part_sum = 0.0_wp
      do c = 1, ncat
         part_sum = part_sum + part_size(c)
      end do
      part_size(0) = max(1.0_wp - part_sum, 0.0_wp)
   end subroutine ice_transport_compress_cell

end module rdb_ice_transport
