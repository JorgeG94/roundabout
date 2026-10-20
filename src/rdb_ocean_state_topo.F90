submodule(rdb_ocean_state) rdb_ocean_state_topo
   !! R2b split: host-only formula-bathymetry setter bodies, moved
   !! verbatim from rdb_ocean_state.F90 (see
   !! src/core/ocean/state/rdb_ocean_state.F90 for the public interface).
   implicit none

contains

   module subroutine set_bathymetry_island(b, grid, max_depth, half_frac)
      !! Flat-bottom basin at `max_depth` everywhere, except a central
      !! square of LAND (`b = 0`, well below `LAND_DEPTH_THRESHOLD`).  The
      !! land square is centred on the physical domain and spans the
      !! central `2*half_frac` fraction of each axis (e.g. `half_frac=0.2`
      !! ⇒ the middle 40 % is land).  `half_frac` is taken from
      !! `&ocean_topo_nml slope_scale` at the call site (no new knob).
      !!
      !! Fills the FULL array including ghost rows (the formula-bathy
      !! ghost-fill gotcha): a ghost `b=0` at a wall-adjacent face sends
      !! the EOS into its ρ=ρ₀ fallback → spurious density jump → blowup.
      !! Ghost rows here inherit the flat `max_depth` (the land square is
      !! strictly interior), so every ghost is ocean.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so the land square is centred on the
      !! WHOLE domain and each rank carves only the part of it that falls
      !! inside its own tile.  On a single rank the offsets are 0 and
      !! global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, half_frac
      integer :: i, j, i_phys, j_phys, ng, ci0, ci1, cj0, cj1, half_i, half_j
      integer :: ioff, joff, nxg, nyg
      real(wp) :: hf

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      ng = grid%nghost
      hf = half_frac
      if (hf <= 0.0_wp) hf = 0.2_wp     ! sensible default if knob left 0
      if (hf > 0.49_wp) hf = 0.49_wp    ! keep at least one wet ring inside

      ! Land-square physical-index bounds (centred on the GLOBAL domain), inclusive.
      half_i = nint(hf*real(nxg, wp))
      half_j = nint(hf*real(nyg, wp))
      ci0 = nxg/2 - half_i + 1
      ci1 = nxg/2 + half_i
      cj0 = nyg/2 - half_j + 1
      cj1 = nyg/2 + half_j

      ! Flat basin everywhere (incl. ghosts), then carve the interior land.
      ! The global physical index of local cell (i,j) is (i_phys + ioff).
      b = max_depth
      do j = 1, size(b, 2)
         j_phys = j - ng
         do i = 1, size(b, 1)
            i_phys = i - ng
            if ((i_phys + ioff) >= ci0 .and. (i_phys + ioff) <= ci1 .and. &
                (j_phys + joff) >= cj0 .and. (j_phys + joff) <= cj1) then
               b(i, j) = 0.0_wp   ! LAND (< LAND_DEPTH_THRESHOLD)
            end if
         end do
      end do
   end subroutine set_bathymetry_island

   module subroutine set_bathymetry_double_drake(b, grid, max_depth, half_frac)
      !! "Double Drake" idealised supercontinent (Ferreira, Marshall &
      !! Campin 2010, *J. Climate*): a flat-bottom global ocean at
      !! `max_depth` with TWO thin meridional wall-continents 90° of
      !! longitude apart, each running from the north pole down to a
      !! southern-channel latitude, leaving a reentrant circumpolar
      !! channel (Drake-Passage analogue) to the south.  Used as the
      !! static-land-mask integration showcase on a spherical periodic-x
      !! sector grid.
      !!
      !! Index-space geometry (mapped to the Ferreira config by the
      !! caller's grid lon/lat extent): physical-x is the 360° periodic
      !! longitude; the two walls sit at the western seam (`i_phys=1`)
      !! and a quarter of the way across (`i_phys≈nx_phys/4`, ≡ 90°
      !! apart on a 360° domain).  Each wall is one cell wide and spans
      !! the NORTHERN `(1-2*half_frac)` fraction of the physical-y axis
      !! (north = pole), leaving the southern `2*half_frac` band fully
      !! open as the reentrant channel.  `half_frac` is taken from
      !! `&ocean_topo_nml slope_scale` (no new knob), default 0.2
      !! ⇒ southern 40 % of the basin is the open channel.
      !!
      !! Fills the FULL array incl. ghosts (formula-bathy ghost-fill
      !! gotcha): the walls extend through the y-ghost rows so the
      !! periodic-x seam wall is land in the halo too; the southern
      !! channel + interior basin stay ocean.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so the two wall meridians and the
      !! channel latitude are located on the WHOLE domain and each rank
      !! carves only the part inside its own tile.  On a single rank the
      !! offsets are 0 and global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, half_frac
      integer :: i, j, i_phys, j_phys, ng, iw1, iw2, j_chan_top
      integer :: ioff, joff, nxg, nyg
      real(wp) :: hf

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      ng = grid%nghost
      hf = half_frac
      if (hf <= 0.0_wp) hf = 0.2_wp
      if (hf > 0.49_wp) hf = 0.49_wp

      ! Two one-cell-wide meridional walls, 90° apart on the GLOBAL 360° axis.
      iw1 = 1                   ! western-seam global wall at global i_phys = 1
      iw2 = max(2, nxg/4)       ! a quarter across ≡ 90° of longitude globally
      ! Channel occupies the southern `2*hf` fraction of the GLOBAL domain;
      ! walls cover the rest (the northern band up to the pole).
      j_chan_top = nint(2.0_wp*hf*real(nyg, wp))

      b = max_depth
      do j = 1, size(b, 2)
         j_phys = j - ng
         do i = 1, size(b, 1)
            i_phys = i - ng
            ! Wall cells: on either meridional line AND north of the
            ! reentrant channel.  Compare GLOBAL physical indices.
            ! j-ghosts north of the channel are wall (walls reach the pole),
            ! south ghosts stay ocean (open channel).
            if (((i_phys + ioff) == iw1 .or. (i_phys + ioff) == iw2) .and. &
                (j_phys + joff) > j_chan_top) then
               b(i, j) = 0.0_wp   ! LAND (< LAND_DEPTH_THRESHOLD)
            end if
         end do
      end do
   end subroutine set_bathymetry_double_drake

   module subroutine set_bathymetry_spoon(b, grid, max_depth, edge_depth, slope_scale)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Fill `b(:,:)` with the MOM6 spoon bathymetry.  In Cartesian
      !! terms (we collapse MOM6's lat/lon factors of `R_earth · π / 180`
      !! into a direct meters scale), the local depth is
      !!
      !!     D(i,j) = D_edge
      !!            + D_0 · sin(π · x_phys / x_len)
      !!                  · (1 − exp((y_phys − y_len) / slope_scale))
      !!
      !! with `D_0 = (max_depth − D_edge) / (1 − exp(−0.5 · y_len /
      !! slope_scale))²`.  Zero at east/west walls, full depth at the
      !! south wall, exponentially decaying to `D_edge` at the north
      !! wall.  Physical interior only; ghost cells retain whatever
      !! value `b` had on entry.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so `x_len` / `y_len` and the
      !! index→position map describe the WHOLE domain on every rank.  On a
      !! single rank the offsets are 0 and global == local, so the fill is
      !! byte-identical to the undecomposed formula.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, edge_depth, slope_scale
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp) :: x_len, y_len, x_phys, y_phys, D_0, denom, D_local
      integer :: i, j, i_phys, j_phys, ng
      integer :: ioff, joff

      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ng = grid%nghost
      x_len = real(grid%nx_global, wp)*grid%dx
      y_len = real(grid%ny_global, wp)*grid%dy
      denom = 1.0_wp - exp(-0.5_wp*y_len/slope_scale)
      D_0 = (max_depth - edge_depth)/(denom*denom)

      ! Fill the full array including ghost rows.  Ghost cells left at
      ! the alloc-time zero cause `h_layer = 0` there, which sends the
      ! EOS into its vanishing-layer fallback (`rho_layer = rho_0`).
      ! That puts a non-physical density jump at every wall-adjacent
      ! face whenever the IC has `T_init /= T_ref` or `S_init /= S_ref`,
      ! and the BPG operator picks it up as a spurious horizontal
      ! pressure gradient.  Confirmed via the seamount test.
      do j = 1, size(b, 2)
         j_phys = j - ng
         ! Global physical y position: local j_phys shifted by joff.
         y_phys = (real(j_phys + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, size(b, 1)
            i_phys = i - ng
            ! Global physical x position: local i_phys shifted by ioff.
            x_phys = (real(i_phys + ioff, wp) - 0.5_wp)*grid%dx
            D_local = edge_depth &
                      + D_0*sin(PI*x_phys/x_len) &
                      *(1.0_wp - exp((y_phys - y_len)/slope_scale))
            ! Safety: the sin term goes negative outside [0, x_len],
            ! and the exp term saturates the depth at the north wall.
            ! Clamp to [edge_depth, max_depth].
            if (D_local < edge_depth) D_local = edge_depth
            if (D_local > max_depth) D_local = max_depth
            b(i, j) = D_local
         end do
      end do
   end subroutine set_bathymetry_spoon

   module subroutine set_bathymetry_seamount(b, grid, max_depth, peak_depth, half_width)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Fill `b(:,:)` with a centred Gaussian seamount bathymetry:
      !!
      !!     D(i,j) = max_depth − (max_depth − peak_depth)
      !!                          · exp(−((x − xc)² + (y − yc)²) / L²)
      !!
      !! with `(xc, yc)` at the centre of the physical domain and
      !! `L = half_width` (gaussian e-folding distance) — expressed in the
      !! GRID coordinate units (metres on Cartesian, degrees on spherical/
      !! curvilinear).  The production dispatch converts the metres
      !! `slope_scale` knob via `topo_length_to_grid_units`; callers that
      !! pass `half_width` directly must match the grid's units.
      !! Reaches `peak_depth` at the bump centre, asymptotes to
      !! `max_depth` far from the bump.
      !!
      !! Designed as a Tier-1.5 bridge test between flat-bottom
      !! analytical setups and full real-bathymetry regional runs.
      !! Tests σ-coord pressure-gradient and Coriolis-advection
      !! kernels under a varying h_layer without the confounds of
      !! a 50× depth ratio + IC stratification + closed-basin
      !! resonance that `tasman_2km.nml` introduces.
      !!
      !! Standard test protocol: uniform T,S (no APE), zero
      !! Coriolis or f-plane only, no wind, walls all sides.  The
      !! expected steady state is **identically zero motion** — any
      !! non-zero u, v, η is a σ-coord-related numerical artefact.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so `(xc, yc)` sits at the centre of the
      !! WHOLE domain on every rank rather than the centre of each tile.
      !! On a single rank the offsets are 0 and global == local, so the fill
      !! is byte-identical to the undecomposed formula.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, peak_depth, half_width
      real(wp) :: x_len, y_len, xc, yc, x_phys, y_phys, r2, depression
      integer :: i, j, i_phys, j_phys, ng, nx_total, ny_total
      integer :: ioff, joff

      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ng = grid%nghost
      x_len = real(grid%nx_global, wp)*grid%dx
      y_len = real(grid%ny_global, wp)*grid%dy
      xc = 0.5_wp*x_len
      yc = 0.5_wp*y_len
      depression = max_depth - peak_depth
      nx_total = size(b, 1)
      ny_total = size(b, 2)

      ! Fill the full array including ghost rows so the EOS, BPG, and
      ! mass-flux operators see a smooth bathymetry on every face they
      ! touch.  Leaving ghosts at zero injects a spurious density jump
      ! at wall-adjacent faces (EOS falls back to `rho_0` where
      ! `h_layer = 0`) — the smoking-gun bug the seamount test was built
      ! to expose.
      do j = 1, ny_total
         j_phys = j - ng
         ! Global physical y position: local j_phys shifted by joff.
         y_phys = (real(j_phys + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, nx_total
            i_phys = i - ng
            ! Global physical x position: local i_phys shifted by ioff.
            x_phys = (real(i_phys + ioff, wp) - 0.5_wp)*grid%dx
            r2 = (x_phys - xc)**2 + (y_phys - yc)**2
            b(i, j) = max_depth - depression*exp(-r2/(half_width*half_width))
         end do
      end do
   end subroutine set_bathymetry_seamount

   module subroutine set_bathymetry_isomip_plus(b, grid, max_depth, x_origin, m_per_grid)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !!
      !! Fill `b(:,:)` with the MISMIP+ / ISOMIP+ analytic bedrock,
      !! Asay-Davis et al. (2016) Eqs. (1)-(4) + Table 1:
      !!
      !!     z_b(x,y) = max( Bx(x) + By(y), z_b,deep )      [Eq. (1)]
      !!
      !! `z_b` is an ELEVATION (positive up, sea level at 0) and is
      !! negative throughout the ISOMIP+ box; Roundabout's `b` is a
      !! DEPTH (positive down), so the last line is `b = -z_b`, floored
      !! at 0 so that a bed which the formula puts ABOVE sea level (it
      !! does for `x < ~140 km`, outside the ISOMIP+ box but reachable if
      !! a caller sets a smaller `x_origin`) is reported as dry land
      !! rather than as a negative depth.  `b = 0` is below
      !! `LAND_DEPTH_THRESHOLD`, so `seed_wet_mask_impl` masks the column
      !! out through the ordinary land path — there is no ISOMIP+ branch
      !! anywhere downstream.
      !!
      !! `max_depth` is the deep clip, i.e. `-z_b,deep`; the protocol
      !! value is `ISOMIP_ZB_DEEP` ⇒ `&ocean_topo_nml max_depth = 720.0`.
      !!
      !! MINIMUM WATER COLUMN.  The protocol (their Sect. 3.1.5) asks for
      !! "the minimum ocean column as thin as can reasonably be achieved"
      !! and leaves the value to the modeller, with the choice being
      !! either to modify the topography or to mark the column land.
      !! Roundabout takes the second option and it is NOT this routine's
      !! job: `&ocean_cavity_dyn_nml h_min_cavity` is the threshold and
      !! `seed_wet_mask_impl(water, land_cutoff=h_min_cavity)` is where it
      !! bites, on `water = b - z_draft`.
      !!
      !! UNITS.  Every Table-1 constant is METRES, as printed.  Grid
      !! positions are in GRID coordinate units (metres on Cartesian,
      !! DEGREES on spherical/curvilinear), so `m_per_grid` converts them
      !! to metres before the formula sees them — the inverse of the
      !! `topo_length_to_grid_units` conversion the spoon/seamount
      !! dispatch applies to `slope_scale`, and the same trap.  On a
      !! Cartesian grid `m_per_grid = 1` exactly and this is the identity.
      !! (The protocol prescribes a Cartesian box; the conversion exists
      !! so a curvilinear caller degrades predictably rather than
      !! silently collapsing the basin flat.)
      !!
      !! Fills the FULL array INCLUDING ghost rows by evaluating the
      !! formula at the ghost index — the CLAUDE.md rule every formula
      !! bathymetry setter follows; a ghost row left at the alloc-time
      !! zero sends the EOS into its `rho_0` vanishing-layer fallback and
      !! puts a spurious density jump at every wall-adjacent face.
      !!
      !! MPI: positions come off the GLOBAL index offsets + extents, so
      !! each rank fills its window of ONE global bed.  Single rank ⇒
      !! offsets 0 ⇒ byte-identical to the undecomposed formula.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth
         !! Deep clip (m, positive down) = `-z_b,deep`.  Protocol: 720.
      real(wp), intent(in) :: x_origin
         !! Absolute MISMIP+ x (m) of the domain's west edge.  Protocol
         !! (ISOMIP+): 320e3.
      real(wp), intent(in) :: m_per_grid
         !! Metres per grid coordinate unit (1 on Cartesian).

      real(wp) :: y_len, x_m, y_m, by, zb
      integer :: i, j, ng, nxt, nyt, ioff, joff

      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxt = size(b, 1)
      nyt = size(b, 2)
      y_len = real(grid%ny_global, wp)*grid%dy*m_per_grid

      do j = 1, nyt
         y_m = (real(j - ng + joff, wp) - 0.5_wp)*grid%dy*m_per_grid
         ! `By` depends on y alone — hoisted out of the i loop.
         by = isomip_plus_by(y_m, y_len)
         do i = 1, nxt
            x_m = x_origin + (real(i - ng + ioff, wp) - 0.5_wp)*grid%dx*m_per_grid
            zb = isomip_plus_bx(x_m) + by
            if (zb < -max_depth) zb = -max_depth     ! Eq. (1) deep clip
            b(i, j) = max(-zb, 0.0_wp)               ! elevation -> depth
         end do
      end do
   end subroutine set_bathymetry_isomip_plus

   pure function nw2_cosbell(x, L) result(c)
      !! Cosine-bell kernel for the Neverworld2 basin: `0.5·(1 + cos(π·min(|x/L|,1)))`.
      !! Peaks at 1 for x=0, decays smoothly to 0 at |x|=L.  Re-derived from
      !! Marques et al. (2022, GMD) "Neverworld2"; MOM6-inspired.
      real(wp), intent(in) :: x, L
      real(wp) :: c
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      c = 0.5_wp*(1.0_wp + cos(PI*min(abs(x/L), 1.0_wp)))
   end function nw2_cosbell

   pure function nw2_spike(x, L) result(s)
      !! Sin-spike kernel for the Neverworld2 basin: `1 − sin(π·min(|x/L|,0.5))`.
      !! Equals 1 at x=0 and drops to 0 at |x|=L/2 (the 0.5 cap stops the sine
      !! re-ascending past its first zero).  Re-derived from Marques et al.
      !! (2022, GMD); MOM6-inspired.
      real(wp), intent(in) :: x, L
      real(wp) :: s
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      s = 1.0_wp - sin(PI*min(abs(x/L), 0.5_wp))
   end function nw2_spike

   module subroutine set_bathymetry_neverworld2(b, grid, max_depth, nl_continent_amp, nl_roughness_amp, min_depth)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !!
      !! Fill `b(:,:)` with the **Neverworld2** idealized-basin bathymetry
      !! (Marques et al. 2022, GMD; MOM6-inspired).  A single Pangaea-style
      !! basin spanning a 60°×140° spherical sector with a re-entrant
      !! (periodic-x) southern channel — the Drake-Passage analog.  Depth is
      !! built in **normalized coordinates**
      !!
      !!     x = (i_phys − 0.5)/nx_phys ∈ [0,1],   y = (j_phys − 0.5)/ny_phys ∈ [0,1]
      !!
      !! which equal MOM6's `(lon − west)/len_lon` and `(lat − south)/len_lat`
      !! on a uniform grid, so no metrics access is needed.  The fractional
      !! depth is
      !!
      !!     D_frac(x,y) = 1
      !!        − 1.1·spike(y−1, 0.12)              ! great northern wall
      !!        − 1.1·spike(y,   0.12)              ! Antarctica (south wall)
      !!        − A_c·[ continents + ridges ]        ! A_c = nl_continent_amp
      !!        − A_r·cos(14πx)·sin(14πy)            ! A_r = nl_roughness_amp
      !!        − A_r·cos(20πx)·cos(20πy)
      !!
      !! clamped `D_frac = max(D_frac, 0)`, then `D = D_frac·max_depth`.  The
      !! continent/ridge block carves two meridional barriers (S-America at the
      !! x-edges, Africa mid-basin) that wall the northern gyres, plus the
      !! Drake/Scotia ridge system that sets the channel sill; `A_c=0` gives an
      !! aquaplanet with the southern channel only (useful for bring-up).
      !!
      !! Sign convention: `b` is bottom depth positive-down (matching the spoon
      !! + flat branches).  Land emerges where D_frac→0 (the wet/dry mask is
      !! seeded downstream from `b ≥ LAND_DEPTH_THRESHOLD`).
      !!
      !! GHOST-ROW FOOTGUN (load-bearing): the loop runs over the FULL array
      !! incl. ghost rows.  Leaving ghosts at the alloc-time zero gives
      !! `h_layer = 0` there, sending the EOS into its vanishing-layer fallback
      !! (`rho_layer = rho_0`) → a spurious density jump at wall-adjacent faces
      !! → ~12-h e-fold blowup (the same trap the spoon/seamount setters guard).
      !! Host only — call `enter_data` afterwards.
      !!
      !! MPI: the normalised coordinates divide by the GLOBAL extents
      !! `grid%nx_global` / `grid%ny_global`, and local cell indices are
      !! mapped to global physical indices with `grid%i_offset_global` /
      !! `grid%j_offset_global`, so `x`/`y` sweep [0,1] across the WHOLE
      !! basin.  Normalised by the LOCAL extents instead, every rank would
      !! rebuild the entire Pangaea basin — walls, continents and all —
      !! inside its own tile.  On a single rank the offsets are 0 and
      !! global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, nl_continent_amp, nl_roughness_amp
      real(wp), intent(in) :: min_depth
         !! Floor on the resulting depth (m).  The fractional-depth formula
         !! produces b in [0, ~1.1*max_depth]; the `1.1*spike` walls and the
         !! continent terms drive b to 0 (land) at the basin edges.  The ocean
         !! C-grid dyn-core does not yet carry a robust wet/dry path, and true
         !! zero-depth land + the resulting sub-metre surface layers blow up
         !! within ~12 h.  Flooring every cell to `min_depth` (MOM6's
         !! MINIMUM_DEPTH approach) turns the continents into shallow shelves
         !! that steer — rather than hard-block — the flow, keeping the basin
         !! all-wet and stable.  Lower values give stronger topographic
         !! steering at the cost of thinner layers (min_depth/nz); true land
         !! barriers wait on the wet/dry-plumbing follow-up.
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp) :: x, y, d_frac, nxp, nyp
      integer :: i, j, i_phys, j_phys, ng, ioff, joff

      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxp = real(grid%nx_global, wp)
      nyp = real(grid%ny_global, wp)

      do j = 1, size(b, 2)
         ! GLOBAL physical indices of this local cell.
         j_phys = j - ng + joff
         y = (real(j_phys, wp) - 0.5_wp)/nyp
         do i = 1, size(b, 1)
            i_phys = i - ng + ioff
            x = (real(i_phys, wp) - 0.5_wp)/nxp
            d_frac = 1.0_wp &
                     - 1.1_wp*nw2_spike(y - 1.0_wp, 0.12_wp) &   ! great northern wall
                     - 1.1_wp*nw2_spike(y, 0.12_wp) &            ! Antarctica (south wall)
                     - nl_continent_amp*( &
                     (1.2_wp*nw2_spike(x, 0.2_wp) + 1.2_wp*nw2_spike(x - 1.0_wp, 0.2_wp)) &
                     *nw2_spike(min(0.0_wp, y - 0.3_wp), 0.2_wp) &                          ! South America
                     + 1.2_wp*nw2_spike(x - 0.5_wp, 0.2_wp) &
                     *nw2_spike(min(0.0_wp, y - 0.55_wp), 0.2_wp) &                         ! Africa
                     + 1.2_wp*(nw2_spike(x, 0.12_wp) + nw2_spike(x - 1.0_wp, 0.12_wp)) &
                     *nw2_spike(max(0.0_wp, y - 0.06_wp), 0.12_wp) &                        ! Antarctic Peninsula
                     + 0.1_wp*(nw2_cosbell(x, 0.1_wp) + nw2_cosbell(x - 1.0_wp, 0.1_wp)) &  ! Drake Passage ridge
                     + 0.5_wp*nw2_cosbell(x - 0.16_wp, 0.05_wp)*(nw2_cosbell(y - 0.18_wp, 0.13_wp)**0.4_wp) &  ! Scotia Arc E
                     + 0.4_wp*(nw2_cosbell(x - 0.09_wp, 0.08_wp)**0.4_wp)*nw2_cosbell(y - 0.26_wp, 0.05_wp) &  ! Scotia Arc N
                     + 0.4_wp*(nw2_cosbell(x - 0.08_wp, 0.08_wp)**0.4_wp)*nw2_cosbell(y - 0.1_wp, 0.05_wp)) &  ! Scotia Arc S
                     - nl_roughness_amp*cos(14.0_wp*PI*x)*sin(14.0_wp*PI*y) &  ! roughness
                     - nl_roughness_amp*cos(20.0_wp*PI*x)*cos(20.0_wp*PI*y)
            if (d_frac < 0.0_wp) d_frac = 0.0_wp
            b(i, j) = max(d_frac*max_depth, min_depth)
         end do
      end do
   end subroutine set_bathymetry_neverworld2

end submodule rdb_ocean_state_topo
