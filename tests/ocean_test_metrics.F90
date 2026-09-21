!! Shared test helper: build a uniform-Cartesian `ocean_metrics_t` from a
!! grid and map it to the device.  Used by the ocean unit tests that call
!! kernels which now take the metrics slot (PGF, dyn steps).  On uniform
!! Cartesian every metric equals the legacy scalar bitwise, so this keeps
!! the analytical tests unchanged.
!!
!! `make_anisotropic_metrics` fills a PURELY-j-varying dx so that the
!! ratio bundle is non-trivial in the x-direction: dxT(i,j) grows with j,
!! dy is uniform.  Designed as a regression sentinel for the corner
!! shear-strain ratio-bundle fix (Finding 1) — the old
!! `(v_diff)*idxCu + (u_diff)*idyCv` form produces a DIFFERENT D_S when
!! dx varies with j because it ignores the per-face metric asymmetry.
module ocean_test_metrics
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_fill_spherical, metrics_fill_tripolar, &
                                metrics_finalize, metrics_apply_land_mask, &
                                metrics_closed_faces_alloc, &
                                adcroft_recip
   implicit none
   private

   public :: make_cartesian_metrics
   public :: make_spherical_metrics
   public :: make_tripolar_metrics
   public :: make_anisotropic_metrics
   public :: destroy_cartesian_metrics

contains

   subroutine make_cartesian_metrics(metrics, grid, wet_mask, nz_closed)
      !! Init + cartesian-fill + finalize [+ land-mask] + device-map a
      !! metrics slot.
      !!
      !! Both optionals exist for the same reason and are orthogonal: each
      !! edits the slot BEFORE the device map, which is the only order that
      !! is correct under `-gpu=...,mem:separate`.
      !!
      !! `wet_mask` (optional, T-cell wet=1 / land=0, `(nx_total,
      !! ny_total)`) is applied by `metrics_apply_land_mask` BEFORE the
      !! device map — the ONLY order that is correct under
      !! `-gpu=...,mem:separate`.  The mask routine is a host-loop
      !! setup-time editor with no `!$acc update`, so a caller that maps
      !! first and masks afterwards leaves the device holding the
      !! ALL-WET metrics: every kernel then transports straight across
      !! the coast on the GPU while the host-side assertions, reading
      !! the masked host copy, look perfectly fine.  Taking the mask
      !! here makes that ordering unskippable.
      !!
      !! `nz_closed` (optional) grows the z-level closed-face masks
      !! `open_u`/`open_v` from their `(1,1,1)` placeholder to full face
      !! size BEFORE the device map.  A test that wants them must pass it
      !! here and must NOT call `metrics_closed_faces_alloc` itself after
      !! this routine: that deallocates the 8-byte placeholders while
      !! they are still mapped, so the device table keeps two stale
      !! entries pointing into freed host memory.  Whatever the host heap
      !! hands out next then reads as "partially present" — a FATAL
      !! runtime error on `-gpu=mem:separate`, at an unrelated array,
      !! whose identity depends only on the allocation order.  The
      !! `metrics_closed_faces_alloc` docstring states the rule; this is
      !! how a test obeys it.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in), optional :: wet_mask(:, :)
         !! T-cell wet/land mask.  Absent ⇒ the masker is not called at
         !! all and the slot keeps its unmasked fill (`init` sources
         !! `wet_*` to 1), which is what every pre-existing caller gets.
      integer, intent(in), optional :: nz_closed
      call metrics%init(grid)
      if (present(nz_closed)) call metrics_closed_faces_alloc(metrics, grid, nz_closed)
      call metrics_fill_cartesian(metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(metrics)
      if (present(wet_mask)) then
         call metrics_apply_land_mask(metrics, wet_mask, grid, &
                                      periodic_x=.false., periodic_y=.false., &
                                      north_fold=.false.)
      end if
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine make_cartesian_metrics

   subroutine make_spherical_metrics(metrics, grid, lon_west, lat_south, &
                                     dlon_deg, dlat_deg, rad_earth)
      !! Init + spherical-fill + finalize + device-map a metrics slot.
      !! Builds a genuinely non-uniform lon-lat sector (dx varies with
      !! latitude), so kernels that secretly still use `1/dx` or
      !! `dx*dy` fail conservation here while passing on Cartesian.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg, rad_earth
      call metrics%init(grid)
      call metrics_fill_spherical(metrics, grid, lon_west, lat_south, &
                                  dlon_deg, dlat_deg, rad_earth)
      call metrics_finalize(metrics)
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine make_spherical_metrics

   subroutine make_tripolar_metrics(metrics, grid, lon_west, lat_south, &
                                    dlon_deg, dlat_deg, rad_earth, phi_join, lon_pole)
      !! Init + tripolar-fill (Murray bipolar cap + seam ghost fold) +
      !! finalize + device-map.  The fill routine folds the north ghosts
      !! and periodic-wraps the east/west ghosts internally.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg
      real(wp), intent(in) :: rad_earth, phi_join, lon_pole
      call metrics%init(grid)
      call metrics_fill_tripolar(metrics, grid, lon_west, lat_south, &
                                 dlon_deg, dlat_deg, rad_earth, phi_join, lon_pole)
      call metrics_finalize(metrics)
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine make_tripolar_metrics

   subroutine destroy_cartesian_metrics(metrics)
      !! Device-unmap + deallocate a metrics slot.
      type(ocean_metrics_t), intent(inout) :: metrics
      call metrics%exit_data()
      !$acc exit data delete(metrics)
      call metrics%destroy()
   end subroutine destroy_cartesian_metrics

   subroutine make_anisotropic_metrics(metrics, grid, dx0, dy, amp)
      !! Build non-uniform metrics where dxT varies with j only:
      !!   dxT(i,j) = dx0 * (1 + amp * (j-1) / ny)
      !! dy is uniform.  All other staggers derived consistently; the ratio
      !! bundle dy_dxBu / dx_dyBu is therefore non-trivial in x — on uniform
      !! Cartesian these are both 1, which makes the D_S old and new forms
      !! algebraically identical.  This metric fill is the regression sentinel
      !! for the corner shear-strain ratio-bundle fix (Finding 1).
      !!
      !! Stagger convention (mirrors metrics_fill_spherical):
      !!   T  at (i,j): lat = (j - 0.5)/ny, dx proportional to (1+amp·lat)
      !!   Cu at (i,j): same latitude row as T, west face of T(i,j)
      !!   Cv at (i,j): corner latitude row, south face of T(i,j)
      !!   Bu at (i,j): SW corner of T(i,j), same latitude as Cv
      !! dyCv / dyBu = dy everywhere (meridional spacing is uniform).
      !! dxCv = dxBu = dx at the corner-latitude row
      !!   = dx0 * (1 + amp * (j-1)/ny)  [corner j-row is between T-rows j-1 and j]
      !! dxCu = dxT (faces share the T-row latitude).
      !! Areas = dx * dy for each stagger.
      !!
      !! Caller must NOT call `metrics_finalize` — this routine does it
      !! internally to keep the helper self-contained.  Device-map is done
      !! via the standard enter_data path.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: dx0
         !! Base zonal spacing (m) at j = 0 latitude.
      real(wp), intent(in) :: dy
         !! Uniform meridional spacing (m).
      real(wp), intent(in) :: amp
         !! Fractional amplitude of the j-variation: dx(j) = dx0*(1+amp*(j-1)/ny).
         !! Typical test value 0.1 gives ~10 % variation across the domain.

      integer :: i, j, nx, ny
      real(wp) :: dx_t, dx_b

      call metrics%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! T-cell and Cu-face metrics (share T-row latitude)
      do j = 1, ny
         dx_t = dx0*(1.0_wp + amp*real(j - 1, wp)/real(ny, wp))
         do i = 1, nx
            metrics%dxT(i, j) = dx_t
            metrics%dyT(i, j) = dy
            metrics%areaT(i, j) = dx_t*dy
         end do
         do i = 1, nx + 1
            metrics%dxCu(i, j) = dx_t
            metrics%dyCu(i, j) = dy
            metrics%dy_cu(i, j) = dy
            metrics%areaCu(i, j) = dx_t*dy
         end do
      end do

      ! Cv-face and Bu-corner metrics (corner latitude rows)
      ! Corner row j sits between T-rows j-1 and j.  Use j-1 as the
      ! latitude index so Bu(i,1) uses dx_b = dx0 (the j=0 corner row).
      do j = 1, ny + 1
         dx_b = dx0*(1.0_wp + amp*real(j - 1, wp)/real(ny, wp))
         do i = 1, nx
            metrics%dxCv(i, j) = dx_b
            metrics%dyCv(i, j) = dy
            metrics%dx_cv(i, j) = dx_b
            metrics%areaCv(i, j) = dx_b*dy
         end do
         do i = 1, nx + 1
            metrics%dxBu(i, j) = dx_b
            metrics%dyBu(i, j) = dy
            metrics%areaBu(i, j) = dx_b*dy
         end do
      end do

      ! Geography left at zero (Cartesian, no lat/lon needed for this helper).
      metrics%geolatT = 0.0_wp
      metrics%geolonT = 0.0_wp
      metrics%geolatBu = 0.0_wp
      metrics%geolonBu = 0.0_wp

      ! Single-source inverses + ratio bundle via finalize.
      call metrics_finalize(metrics)

      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine make_anisotropic_metrics

end module ocean_test_metrics
