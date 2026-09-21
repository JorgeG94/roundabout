!! Porous barriers — subgrid sill/strait blocking of C-grid face widths.
module rdb_ocean_porous
   !! Adcroft (2013) three-parameter porous-barrier fit for the ocean
   !! C-grid dyn-core.
   !!
   !! A model face is a straight segment of length `dy_cu` (u-face) or
   !! `dx_cv` (v-face), but the real seafloor along that segment is not
   !! flat: a strait or a sill leaves only PART of the segment open at a
   !! given depth.  Porous barriers replace the binary open/closed face
   !! with a depth-dependent OPEN FRACTION, so a deep sill blocks the
   !! bottom layers while the surface layers stay fully open.  The face
   !! stays a single C-grid degree of freedom — only the transport
   !! width narrows.
   !!
   !! The along-face seafloor is summarised by three numbers per face,
   !! all TOPOGRAPHIC HEIGHTS (positive up, negative below the sea
   !! surface) so they compare directly against interface heights.  NOTE
   !! the ocean path's `barotropic%b` is the opposite sign — a reference
   !! column DEPTH, positive down (`SSH = sum(h_layer) - b`) — so the
   !! caller negates it once when filling `metrics%por_bed`:
   !!   * `d_min` — deepest along-face point (most negative),
   !!   * `d_max` — shallowest along-face point (least negative),
   !!   * `d_avg` — mean along-face height, `d_min <= d_avg <= d_max`.
   !!
   !! Adcroft's fit picks the one-parameter family of monotone profiles
   !! that reproduces exactly those three numbers.  With
   !!   `m    = (d_avg - d_min) / (d_max - d_min)`   (nondim, in [0,1])
   !!   `zeta = (eta   - d_min) / (d_max - d_min)`   (nondim)
   !!   `a    = (1 - m) / m`
   !! the OPEN WIDTH FRACTION at an interface of height `eta` is
   !!   `w(eta) = 0`                        for `eta <= d_min`,
   !!   `w(eta) = zeta**(1/a)`              for `m < 1/2`,
   !!   `w(eta) = zeta`                     for `m = 1/2`,
   !!   `w(eta) = 1 - (1 - zeta)**a`        for `m > 1/2`,
   !!   `w(eta) = 1`                        for `eta > d_max`,
   !! and its vertical integral from the bottom (units of length) is
   !!   `A(eta) = 0`                                      `eta <= d_min`,
   !!   `A(eta) = (d_max-d_min) * (1-m) * zeta**(1/(1-m))`   `m < 1/2`,
   !!   `A(eta) = (d_max-d_min) * 0.5*zeta*zeta`             `m = 1/2`,
   !!   `A(eta) = (d_max-d_min) * (zeta - m + m*(1-zeta)**(1/m))`  `m > 1/2`,
   !!   `A(eta) = eta - d_avg`                             `eta > d_max`.
   !! `dA/d(eta) = w(eta)` identically, and both branches join
   !! continuously at `eta = d_max` where `A = (d_max-d_min)*(1-m)`.
   !!
   !! The LAYER-AVERAGED open AREA fraction of layer `k`, whose lower and
   !! upper interfaces sit at face heights `eta_lo` / `eta_hi`, is the
   !! mean of `w` over the layer, obtained exactly from the cumulative
   !! integral:
   !!   `por_face_area(k) = min(1, (A(eta_hi) - A(eta_lo))/(eta_hi - eta_lo))`
   !! with a zero fallback for a vanishing layer.  It multiplies the face
   !! width in the transport: `uh = u * h_face * dy_cu * por_face_area_u`.
   !!
   !! Reference: Adcroft, A. (2013), "Representation of topography by
   !! porous barriers and objective interpolation of topographic data",
   !! Ocean Modelling 67, 13-27.  MOM6's `MOM_porous_barriers` was the
   !! inspiration for the discrete layer-averaging and the eta-at-velocity
   !! interpolation options; this is an independent implementation
   !! re-derived from the paper's fit.
   !!
   !! SUBGRID DATA CAVEAT.  The fit needs min/max/mean of the
   !! HIGH-RESOLUTION seafloor along each face — a statistic that can only
   !! come from a bathymetry dataset finer than the model grid (MOM6 reads
   !! it from an offline-generated `topog_edge.nc`).  Roundabout has no such
   !! file plumbing yet, so the shipped source is `POROUS_SOURCE_RESOLVED`:
   !! min/max/mean of the RESOLVED bathymetry sampled at three along-face
   !! points (south corner, midpoint, north corner for a u-face).  That is
   !! a genuine along-face statistic of the data we have, NOT true subgrid
   !! information — it captures along-face slope but cannot see structure
   !! below the grid scale.  `POROUS_SOURCE_FILE` is the honest
   !! MOM6-parity path and fails loud pending the file-forcing backend
   !! (the same gate the barotropic wave-drag `form="file"` waits on).
   !!
   !! WHAT THE PROXY ACTUALLY DOES — read this before trusting it.  All
   !! three samples are CELL-CENTRE AVERAGES, so the statistic only sees
   !! variation ALONG the face.
   !!   * Bathymetry uniform ALONG the face — a ridge or shelf break that
   !!     runs parallel to it, the common case — collapses all three
   !!     samples onto one value.  The fit would then be a STEP at the
   !!     two-cell mean height: a hard wall on every layer below it, on a
   !!     face the grid already resolves as open.  That wall is
   !!     manufactured by the fill routine, not information, so
   !!     `porous_update_face_areas` treats a DEGENERATE face
   !!     (`d_max <= d_min`) as FULLY OPEN and blocks nothing.  A flat
   !!     seafloor is the same degenerate case, which is also what makes
   !!     the resolved source a literal no-op on a flat basin.
   !!   * Where the bathymetry DOES vary along the face, the narrowing is
   !!     set by the spread of the three corner-mean samples.  A face
   !!     whose along-face spread is small compared with the water depth
   !!     still yields a near-step over that narrow spread — the fit is
   !!     only ever as smooth as the sampled spread.  Genuinely graded
   !!     narrowing needs real subgrid data (`source="file"`).
   !!
   !! WET-CELL GATING.  A corner sample averages four cells; if any of
   !! them is LAND its elevation would pull `d_max` up and manufacture
   !! blockage on a face the grid resolves as fully open.
   !! `porous_fill_stats_resolved` therefore drops a corner sample whose
   !! four-cell stencil is not entirely wet, falling back to the two-cell
   !! face midpoint for that sample.
   !!
   !! VANISHING-LAYER THRESHOLD.  A layer whose face thickness is at or
   !! below `H_VANISHED` (1.5e-4 m) is given `por = 0`.  MOM6 instead uses
   !! `Angstrom_Z` (1e-10 m), 1.5 million times smaller, so a layer
   !! between the two thresholds gets a real open fraction there and zero
   !! here.  The consequence is that `por_col` (the column-integrated
   !! fraction the barotropic width carries) is the thickness-weighted
   !! mean of the per-layer fractions only over the NON-vanished layers.
   !! `H_VANISHED` is the repo-wide dynamic-vanish constant
   !! (`rdb_constants`), not a porous-barrier choice.
   !!
   !! NOT NARROWED.  `areaCu` / `areaCv` keep their full geometric value:
   !! only the TRANSPORT widths are narrowed (`dy_cu` / `dx_cv` per layer
   !! and `dy_cu_bt` / `dx_cv_bt` for the barotropic mode).  The
   !! barotropic Coriolis and KE terms therefore combine a narrowed
   !! transport with an un-narrowed cell area — the same split MOM6 has,
   !! but worth knowing before reading a BT energy budget under a strong
   !! barrier.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS
   implicit none
   private

   public :: porous_open_width, porous_cum_area, porous_eta_face
   public :: porous_fill_stats_resolved
   public :: porous_stats_are_ordered
   public :: porous_update_face_areas
   public :: porous_narrow_3d
   public :: closed_faces_update_bt_widths
   public :: parse_porous_source, parse_porous_eta_interp
   public :: POROUS_SOURCE_RESOLVED, POROUS_SOURCE_FILE
   public :: POROUS_ETA_MAX, POROUS_ETA_MIN, POROUS_ETA_ARITH, POROUS_ETA_HARM

   ! ---- Subgrid-statistics source enum (`&ocean_porous_nml source`) ----
   integer, parameter :: POROUS_SOURCE_RESOLVED = 0
      !! Along-face min/max/mean of the RESOLVED bathymetry (three
      !! samples: the two face corners + the face midpoint).  A proxy
      !! for true subgrid statistics — see the module caveat.
   integer, parameter :: POROUS_SOURCE_FILE = 1
      !! Offline subgrid-bathymetry file (MOM6 `topog_edge.nc`).  Not
      !! implemented — fails loud at configure.

   ! ---- Interface-height-at-velocity-point interpolation enum ----
   integer, parameter :: POROUS_ETA_MAX = 0
      !! Higher (shallower) of the two adjacent interface heights — the
      !! default, and the LEAST blocking of the four rules: `w` is
      !! monotone increasing in the interface height, so raising both
      !! interfaces raises the layer-averaged open fraction.
   integer, parameter :: POROUS_ETA_MIN = 1
      !! Lower (deeper) of the two adjacent interface heights — the MOST
      !! blocking rule, by the same monotonicity.
   integer, parameter :: POROUS_ETA_ARITH = 2
      !! Arithmetic mean of the two adjacent interface heights.
   integer, parameter :: POROUS_ETA_HARM = 3
      !! Harmonic mean of the two adjacent interface heights.

contains

   ! =================================================================
   ! Enum parsing
   ! =================================================================

   pure function parse_porous_source(name) result(src)
      !! Map `&ocean_porous_nml source` onto `POROUS_SOURCE_*`.  An
      !! unrecognised string returns -1 so the caller can fail loud.
      character(len=*), intent(in) :: name
      integer :: src
      select case (trim(adjustl(name)))
      case ("resolved")
         src = POROUS_SOURCE_RESOLVED
      case ("file")
         src = POROUS_SOURCE_FILE
      case default
         src = -1
      end select
   end function parse_porous_source

   pure function parse_porous_eta_interp(name) result(interp)
      !! Map `&ocean_porous_nml eta_interp` onto `POROUS_ETA_*`.  An
      !! unrecognised string returns -1 so the caller can fail loud.
      character(len=*), intent(in) :: name
      integer :: interp
      select case (trim(adjustl(name)))
      case ("max")
         interp = POROUS_ETA_MAX
      case ("min")
         interp = POROUS_ETA_MIN
      case ("arithmetic")
         interp = POROUS_ETA_ARITH
      case ("harmonic")
         interp = POROUS_ETA_HARM
      case default
         interp = -1
      end select
   end function parse_porous_eta_interp

   ! =================================================================
   ! Adcroft (2013) three-parameter fit
   ! =================================================================

   pure function porous_open_width(d_min, d_max, d_avg, eta) result(w)
      !! Open WIDTH fraction of a face at interface height `eta`
      !! (dimensionless, in `[0, 1]`).  Zero when the interface is at or
      !! below the deepest along-face point, one when it is above the
      !! shallowest.  A degenerate face (`d_max <= d_min`, i.e. a flat
      !! along-face seafloor) reduces to the binary open/closed step, so
      !! the fit never divides by zero.
      real(wp), intent(in) :: d_min
         !! Deepest along-face topographic height (m, positive up).
      real(wp), intent(in) :: d_max
         !! Shallowest along-face topographic height (m, positive up).
      real(wp), intent(in) :: d_avg
         !! Mean along-face topographic height (m, positive up).
      real(wp), intent(in) :: eta
         !! Interface height at the face (m, positive up).
      real(wp) :: w
      real(wp) :: m, a, zeta, drange

      if (eta <= d_min) then
         w = 0.0_wp
      else if (eta > d_max) then
         w = 1.0_wp
      else if (d_avg <= d_min) then
         ! Degenerate m = 0: the profile collapses to a step at d_min, so
         ! everything above it is open.  Guarded explicitly because
         ! `a = (1-m)/m` would divide by zero here.
         w = 1.0_wp
      else if (d_avg >= d_max) then
         ! Degenerate m = 1: a step at d_max — nothing below it is open.
         w = 0.0_wp
      else
         drange = d_max - d_min
         m = (d_avg - d_min)/drange
         zeta = (eta - d_min)/drange
         if (m < 0.5_wp) then
            ! a = (1-m)/m, so 1/a = m/(1-m) = (d_avg-d_min)/(d_max-d_avg).
            a = (1.0_wp - m)/m
            w = zeta**(1.0_wp/a)
         else if (m > 0.5_wp) then
            a = (1.0_wp - m)/m
            w = 1.0_wp - (1.0_wp - zeta)**a
         else
            w = zeta
         end if
      end if
   end function porous_open_width

   pure function porous_cum_area(d_min, d_max, d_avg, eta) result(area)
      !! Cumulative open area of a face from the deepest along-face point
      !! up to interface height `eta`, per unit face length (m — it is
      !! the vertical integral of `porous_open_width`).  Layer-averaged
      !! open fractions are differences of this function divided by the
      !! layer thickness, which is exact (no quadrature error) because
      !! `d(area)/d(eta) = porous_open_width(eta)` identically.
      real(wp), intent(in) :: d_min
         !! Deepest along-face topographic height (m, positive up).
      real(wp), intent(in) :: d_max
         !! Shallowest along-face topographic height (m, positive up).
      real(wp), intent(in) :: d_avg
         !! Mean along-face topographic height (m, positive up).
      real(wp), intent(in) :: eta
         !! Interface height at the face (m, positive up).
      real(wp) :: area
      real(wp) :: m, zeta, drange

      if (eta <= d_min) then
         area = 0.0_wp
      else if (eta > d_max) then
         ! Above the sill the face is fully open, so the integral grows
         ! linearly; the offset `d_avg` is what makes the mean height of
         ! the fit equal the prescribed `d_avg`.
         area = eta - d_avg
      else if (d_avg <= d_min) then
         ! Degenerate m = 0 (step at d_min): fully open above it, so the
         ! integral is the same linear form the above-d_max branch uses.
         area = eta - d_avg
      else if (d_avg >= d_max) then
         ! Degenerate m = 1 (step at d_max): nothing open below it.
         area = 0.0_wp
      else
         drange = d_max - d_min
         m = (d_avg - d_min)/drange
         zeta = (eta - d_min)/drange
         if (m < 0.5_wp) then
            area = drange*((1.0_wp - m)*zeta**(1.0_wp/(1.0_wp - m)))
         else if (m > 0.5_wp) then
            area = drange*(zeta - m + m*((1.0_wp - zeta)**(1.0_wp/m)))
         else
            area = drange*(0.5_wp*zeta*zeta)
         end if
      end if
   end function porous_cum_area

   pure function porous_eta_face(z_a, z_b, interp) result(eta)
      !! Interface height at a velocity point from the two adjacent
      !! cell-centre interface heights.  MOM6's `PORBAR_ETA_INTERP`
      !! options.  `POROUS_ETA_MAX` (the higher, i.e. shallower,
      !! interface) is the default and the LEAST blocking: the open width
      !! `w` is monotone increasing in the interface height, so the rule
      !! that returns the larger height leaves the most of the face open.
      !! `POROUS_ETA_MIN` is the most blocking.
      real(wp), intent(in) :: z_a
         !! Interface height in the first adjacent cell (m, positive up).
      real(wp), intent(in) :: z_b
         !! Interface height in the second adjacent cell (m, positive up).
      integer, intent(in) :: interp
         !! One of the `POROUS_ETA_*` enum values.
      real(wp) :: eta

      select case (interp)
      case (POROUS_ETA_MAX)
         eta = max(z_a, z_b)
      case (POROUS_ETA_MIN)
         eta = min(z_a, z_b)
      case (POROUS_ETA_ARITH)
         eta = 0.5_wp*(z_a + z_b)
      case (POROUS_ETA_HARM)
         ! Only meaningful when both interfaces are on the SAME side of
         ! the datum (the ocean case, both negative).  A straddling pair
         ! can drive the denominator through zero, where the H_DIV_EPS
         ! armour prevents the division but not a large result.
         eta = 2.0_wp*(z_a*z_b)/(z_a + z_b + H_DIV_EPS)
      case default
         ! UNREACHABLE.  MOM6 issues a FATAL here; this is a `pure`
         ! device-side function and cannot, so the fail-loud lives at the
         ! host boundary instead: `parse_porous_eta_interp` returns -1 for
         ! any unrecognised string and `configure_ocean_porous`
         ! `error stop`s on it, so no other value can reach this kernel.
         ! Falling back to MAX keeps the function total.
         eta = max(z_a, z_b)
      end select
   end function porous_eta_face

   ! =================================================================
   ! Setup: along-face statistics from the resolved bathymetry
   ! =================================================================

   pure subroutine porous_fill_stats_resolved(nx, ny, b, wet_t, &
                                              dmin_u, dmax_u, davg_u, &
                                              dmin_v, dmax_v, davg_v)
      !! Fill the per-face `d_min` / `d_max` / `d_avg` from the RESOLVED
      !! bottom elevation, sampled at three points ALONG each face: the
      !! two end corners and the midpoint.  A corner sample is the mean
      !! of the four cells around it, the midpoint the mean of the two
      !! cells the face separates — all three are genuine points on the
      !! resolved seafloor along the face segment.
      !!
      !! WET GATING.  A corner sample is used ONLY when all four cells in
      !! its stencil are wet; otherwise it falls back to the two-cell face
      !! midpoint.  Without this a wet-wet face with one diagonal LAND
      !! neighbour picks the land elevation up into `d_max` and blocks a
      !! face the grid fully resolves as open (measured: a single 4000 m
      !! land diagonal on a 4000 m column gives ~37% spurious blockage on
      !! the deepest layers).  When both corners are gated out all three
      !! samples coincide, the face statistic is degenerate, and
      !! `porous_update_face_areas` leaves it fully open.
      !!
      !! This is a PROXY for true subgrid statistics (see the module
      !! caveat): it resolves along-face slope but is blind to structure
      !! below the grid scale.  A flat seafloor — or ANY seafloor uniform
      !! along the face — gives `d_min = d_max = d_avg`, the degenerate
      !! case the kernel leaves fully open.  That is the bit-identity
      !! property flat-bottom configurations rely on, and the reason a
      !! ridge running parallel to the face does not become a wall.
      !!
      !! The outermost face columns (`i = 1`, `i = nx+1` for u;
      !! `j = 1`, `j = ny+1` for v) are left at zero: they have no adjacent
      !! cell pair, `porous_update_face_areas` skips them, and their
      !! metric width is already zero.
      !!
      !! Setup-phase host code: plain sequential loops, NOT
      !! `do concurrent` (this runs before `enter_data`, where a DC loop
      !! would round-trip the unmapped arrays through the device).
      integer, intent(in) :: nx, ny
         !! `grid%nx_total`, `grid%ny_total`.
      real(wp), intent(in) :: b(nx, ny)
         !! Bottom elevation at cell centres (m, positive up).
      real(wp), intent(in) :: wet_t(nx, ny)
         !! T-cell wet (1) / land (0) mask (`ocean_metrics_t%wet_T`).
         !! All-ones on a run with no land, where every corner sample is
         !! kept and the statistic is the un-gated one.
      real(wp), intent(out) :: dmin_u(nx + 1, ny), dmax_u(nx + 1, ny), davg_u(nx + 1, ny)
         !! u-face along-face deepest / shallowest / mean height (m).
      real(wp), intent(out) :: dmin_v(nx, ny + 1), dmax_v(nx, ny + 1), davg_v(nx, ny + 1)
         !! v-face along-face deepest / shallowest / mean height (m).

      integer :: i, j
      real(wp) :: s_lo, s_mid, s_hi

      ! ---- u-faces: face i separates cells (i-1,j) and (i,j) ----
      ! Samples run south -> north along the face segment.
      dmin_u = 0.0_wp
      dmax_u = 0.0_wp
      davg_u = 0.0_wp
      do j = 2, ny - 1
         do i = 2, nx
            s_mid = 0.5_wp*(b(i - 1, j) + b(i, j))
            if (corner_is_wet(wet_t(i - 1, j - 1), wet_t(i, j - 1), &
                              wet_t(i - 1, j), wet_t(i, j))) then
               s_lo = 0.25_wp*(b(i - 1, j - 1) + b(i, j - 1) + b(i - 1, j) + b(i, j))
            else
               s_lo = s_mid
            end if
            if (corner_is_wet(wet_t(i - 1, j + 1), wet_t(i, j + 1), &
                              wet_t(i - 1, j), wet_t(i, j))) then
               s_hi = 0.25_wp*(b(i - 1, j + 1) + b(i, j + 1) + b(i - 1, j) + b(i, j))
            else
               s_hi = s_mid
            end if
            dmin_u(i, j) = min(s_lo, s_mid, s_hi)
            dmax_u(i, j) = max(s_lo, s_mid, s_hi)
            ! Simpson (1,4,1)/6 — the consistent along-face MEAN of three
            ! evenly spaced samples.  A plain (1,1,1)/3 average would
            ! under-weight the midpoint and bias `d_avg` toward the ends.
            !
            ! Bracketed because the invariant `d_min <= d_avg <= d_max` is
            ! exact only in exact arithmetic: the weights are a convex
            ! combination, but three near-equal samples can round the sum
            ! an ulp outside the range, and `m = (d_avg-d_min)/(d_max-d_min)`
            ! must stay in [0,1].  The bracket moves nothing else.
            davg_u(i, j) = min(dmax_u(i, j), max(dmin_u(i, j), &
                                                 (s_lo + 4.0_wp*s_mid + s_hi)/6.0_wp))
         end do
      end do
      ! Outermost rows have no cross-face neighbour for the corner
      ! samples: fall back to the two-cell midpoint (a flat face, so the
      ! fit degenerates to the open/closed step and blocks nothing extra).
      do i = 2, nx
         s_mid = 0.5_wp*(b(i - 1, 1) + b(i, 1))
         dmin_u(i, 1) = s_mid
         dmax_u(i, 1) = s_mid
         davg_u(i, 1) = s_mid
         s_mid = 0.5_wp*(b(i - 1, ny) + b(i, ny))
         dmin_u(i, ny) = s_mid
         dmax_u(i, ny) = s_mid
         davg_u(i, ny) = s_mid
      end do

      ! ---- v-faces: face j separates cells (i,j-1) and (i,j) ----
      ! Samples run west -> east along the face segment.
      dmin_v = 0.0_wp
      dmax_v = 0.0_wp
      davg_v = 0.0_wp
      do j = 2, ny
         do i = 2, nx - 1
            s_mid = 0.5_wp*(b(i, j - 1) + b(i, j))
            if (corner_is_wet(wet_t(i - 1, j - 1), wet_t(i - 1, j), &
                              wet_t(i, j - 1), wet_t(i, j))) then
               s_lo = 0.25_wp*(b(i - 1, j - 1) + b(i - 1, j) + b(i, j - 1) + b(i, j))
            else
               s_lo = s_mid
            end if
            if (corner_is_wet(wet_t(i + 1, j - 1), wet_t(i + 1, j), &
                              wet_t(i, j - 1), wet_t(i, j))) then
               s_hi = 0.25_wp*(b(i + 1, j - 1) + b(i + 1, j) + b(i, j - 1) + b(i, j))
            else
               s_hi = s_mid
            end if
            dmin_v(i, j) = min(s_lo, s_mid, s_hi)
            dmax_v(i, j) = max(s_lo, s_mid, s_hi)
            ! Bracketed — see the zonal twin.
            davg_v(i, j) = min(dmax_v(i, j), max(dmin_v(i, j), &
                                                 (s_lo + 4.0_wp*s_mid + s_hi)/6.0_wp))
         end do
      end do
      do j = 2, ny
         s_mid = 0.5_wp*(b(1, j - 1) + b(1, j))
         dmin_v(1, j) = s_mid
         dmax_v(1, j) = s_mid
         davg_v(1, j) = s_mid
         s_mid = 0.5_wp*(b(nx, j - 1) + b(nx, j))
         dmin_v(nx, j) = s_mid
         dmax_v(nx, j) = s_mid
         davg_v(nx, j) = s_mid
      end do
   end subroutine porous_fill_stats_resolved

   pure function corner_is_wet(w1, w2, w3, w4) result(ok)
      !! `.true.` iff all four cells contributing to a corner sample are
      !! wet.  The masks are real 0/1 (`ocean_metrics_t%wet_T`), so the
      !! test is a mid-point comparison rather than an equality.
      real(wp), intent(in) :: w1, w2, w3, w4
      logical :: ok
      ok = (w1 > 0.5_wp) .and. (w2 > 0.5_wp) .and. &
           (w3 > 0.5_wp) .and. (w4 > 0.5_wp)
   end function corner_is_wet

   pure function porous_stats_are_ordered(n1, n2, dmin, dmax, davg) result(ok)
      !! `.true.` iff every face satisfies `d_min <= d_avg <= d_max`, the
      !! invariant the whole fit rests on (`m = (d_avg-d_min)/(d_max-d_min)`
      !! must lie in `[0,1]`).
      !!
      !! `porous_fill_stats_resolved` cannot violate it — `d_avg` is a
      !! convex combination of the same three samples `d_min`/`d_max`
      !! bracket.  A FILE-backed source can, so this is the assertion the
      !! reader boundary owes: `configure_ocean_porous` calls it after
      !! filling, whichever source produced the numbers, and fails loud.
      integer, intent(in) :: n1, n2
         !! Face-array extents.
      real(wp), intent(in) :: dmin(n1, n2), dmax(n1, n2), davg(n1, n2)
         !! Along-face deepest / shallowest / mean height (m, positive up).
      logical :: ok
      integer :: i, j
      ok = .true.
      do j = 1, n2
         do i = 1, n1
            if (davg(i, j) < dmin(i, j) .or. davg(i, j) > dmax(i, j)) then
               ok = .false.
               return
            end if
         end do
      end do
   end function porous_stats_are_ordered

   ! =================================================================
   ! Applying the open fractions to a face-staggered 3D field
   ! =================================================================

   pure subroutine porous_narrow_3d(n1, n2, nz, por, arr)
      !! Multiply a face-staggered per-layer field by the open-area
      !! fraction: `arr <- arr * por`.
      !!
      !! Deliberately a SEPARATE pass rather than a factor folded into the
      !! transport loops.  Folding it in would put
      !! `metrics%por_face_area_u(i,j,k)` inside a `do concurrent` whose
      !! implicit data clause is generated from the LOOP bounds, so the
      !! `(1,1,1)` placeholder a knob-off run carries would be reported
      !! partially present and abort under `mem:separate`.  With the
      !! multiply hoisted behind a host-side `if (metrics%use_porous)`
      !! there is no kernel launch at all when the knob is off — which
      !! also leaves the transport loops textually untouched, so
      !! bit-identity is by construction rather than by argument.  The
      !! extra pass is bandwidth-bound and only runs when the (opt-in)
      !! scheme is active.
      integer, intent(in) :: n1, n2, nz
         !! Face-array extents (`nx+1, ny, nz` for u; `nx, ny+1, nz` for v).
      real(wp), intent(in) :: por(n1, n2, nz)
         !! Layer-averaged open-area fraction (nondim, `[0,1]`).
      real(wp), intent(inout) :: arr(n1, n2, nz)
         !! Face-staggered field to narrow (a mass transport).
      integer :: i, j, k

      do concurrent(k=1:nz, j=1:n2, i=1:n1)
         arr(i, j, k) = arr(i, j, k)*por(i, j, k)
      end do
   end subroutine porous_narrow_3d

   ! =================================================================
   ! Per-stage device recompute of the layer-averaged open fractions
   ! =================================================================

   pure subroutine porous_update_face_areas(nx, ny, nz, interp, mask_depth, &
                                            b, h_layer, &
                                            dmin_u, dmax_u, davg_u, &
                                            dmin_v, dmax_v, davg_v, &
                                            dy_cu, dx_cv, &
                                            por_u, por_v, dy_cu_bt, dx_cv_bt)
      !! Recompute the layer-averaged open-area fractions from the
      !! CURRENT layer thicknesses.  Interface-height dependent, so this
      !! runs once per RK2 stage on the device.
      !!
      !! Column recurrence (bottom-up, `k=1` is the bed): the face
      !! interface height is interpolated from the two adjacent columns'
      !! running interface heights, the cumulative open area is evaluated
      !! there, and the layer fraction is the increment over the layer
      !! divided by the layer's face thickness.  Everything is a running
      !! scalar — no per-column workspace, so there is nothing to
      !! allocate and nothing to map.
      !!
      !! Three faces are left FULLY OPEN (`por = 1`, BT width untouched):
      !!   * array-edge faces (no adjacent cell pair; metric width is
      !!     already zero),
      !!   * faces whose mean along-face height is at or above
      !!     `mask_depth` — MOM6's `PORBAR_MASKING_DEPTH` gate; porous
      !!     barriers are a deep-sill parameterization and should not
      !!     narrow shelf faces the grid already resolves,
      !!   * DEGENERATE faces, `d_max <= d_min`.  The fit's limit there is
      !!     a step at `d_min`, which on a seafloor uniform ALONG the face
      !!     puts a hard wall at the two-cell mean height — an artifact of
      !!     the resolved-bathymetry proxy, not subgrid information (see
      !!     the module docstring).  Leaving them open also makes a flat
      !!     basin a literal no-op, which is what the bit-identity
      !!     argument for `source="resolved"` rests on.
      !!
      !! A wholly VANISHED column (top interface within `H_VANISHED` of
      !! the bed) gets `por = 0` on every layer AND a zero BT width, so
      !! the barotropic mode and the layers agree that nothing can move
      !! through it.  Handing the BT the un-narrowed width there would let
      !! it transport at full width across a face every layer has blocked.
      integer, intent(in) :: nx, ny, nz
         !! `grid%nx_total`, `grid%ny_total`, number of layers.
      integer, intent(in) :: interp
         !! Interface-at-velocity-point rule, a `POROUS_ETA_*` value.
      real(wp), intent(in) :: mask_depth
         !! Gate height (m, positive up, `<= 0`): faces with
         !! `d_avg >= mask_depth` stay fully open.
      real(wp), intent(in) :: b(nx, ny)
         !! Bottom elevation at cell centres (m, positive up).
      real(wp), intent(in) :: h_layer(nx, ny, nz)
         !! Layer thicknesses (m), `k=1` the bed layer.
      real(wp), intent(in) :: dmin_u(nx + 1, ny), dmax_u(nx + 1, ny), davg_u(nx + 1, ny)
         !! u-face along-face deepest / shallowest / mean height (m).
      real(wp), intent(in) :: dmin_v(nx, ny + 1), dmax_v(nx, ny + 1), davg_v(nx, ny + 1)
         !! v-face along-face deepest / shallowest / mean height (m).
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
         !! Un-narrowed open u-face width for transport (m).
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
         !! Un-narrowed open v-face width for transport (m).
      real(wp), intent(out) :: por_u(nx + 1, ny, nz)
         !! u-face layer-averaged open-area fraction (nondim, `[0,1]`).
      real(wp), intent(out) :: por_v(nx, ny + 1, nz)
         !! v-face layer-averaged open-area fraction (nondim, `[0,1]`).
      real(wp), intent(out) :: dy_cu_bt(nx + 1, ny)
         !! u-face width the BAROTROPIC substep transports on (m):
         !! `dy_cu` scaled by the COLUMN-INTEGRATED open fraction.
      real(wp), intent(out) :: dx_cv_bt(nx, ny + 1)
         !! v-face twin (m).

      integer :: i, j, k
      real(wp) :: z_a, z_b, e_lo, e_hi, a_lo, a_hi, dz
      real(wp) :: e_bed, a_bed, por_col

      ! ---- u-faces ----
      do concurrent(j=1:ny, i=1:nx + 1) &
         local(k, z_a, z_b, e_lo, e_hi, a_lo, a_hi, dz, e_bed, a_bed, por_col)
         if (i < 2 .or. i > nx) then
            ! Array-edge faces have no adjacent cell pair; their metric
            ! width is already zero, so leave them fully open.
            do k = 1, nz
               por_u(i, j, k) = 1.0_wp
            end do
            dy_cu_bt(i, j) = dy_cu(i, j)
         else if (davg_u(i, j) >= mask_depth) then
            do k = 1, nz
               por_u(i, j, k) = 1.0_wp
            end do
            dy_cu_bt(i, j) = dy_cu(i, j)
         else if (dmax_u(i, j) <= dmin_u(i, j)) then
            ! Degenerate face: the along-face statistic carries no
            ! information, so block nothing (see the routine docstring).
            do k = 1, nz
               por_u(i, j, k) = 1.0_wp
            end do
            dy_cu_bt(i, j) = dy_cu(i, j)
         else
            z_a = b(i - 1, j)
            z_b = b(i, j)
            e_lo = porous_eta_face(z_a, z_b, interp)
            a_lo = porous_cum_area(dmin_u(i, j), dmax_u(i, j), davg_u(i, j), e_lo)
            e_bed = e_lo
            a_bed = a_lo
            do k = 1, nz
               z_a = z_a + h_layer(i - 1, j, k)
               z_b = z_b + h_layer(i, j, k)
               e_hi = porous_eta_face(z_a, z_b, interp)
               a_hi = porous_cum_area(dmin_u(i, j), dmax_u(i, j), davg_u(i, j), e_hi)
               dz = e_hi - e_lo
               if (dz > H_VANISHED) then
                  por_u(i, j, k) = clamp_fraction((a_hi - a_lo)/dz)
               else
                  ! Vanished layer at this face: no open area to speak of.
                  por_u(i, j, k) = 0.0_wp
               end if
               e_lo = e_hi
               a_lo = a_hi
            end do
            ! Column-integrated open fraction: the SAME layer-average
            ! formula applied over the whole water column, which is
            ! identically the thickness-weighted mean of the per-layer
            ! fractions.  This is what makes the barotropic transport see
            ! the barrier (MOM6 folds the equivalent quantity into its
            ! BT_cont face area).  A wholly VANISHED column has every
            ! layer blocked above, so the BT width goes to zero too — the
            ! two modes must agree about a face nothing can pass through.
            if (e_lo - e_bed > H_VANISHED) then
               por_col = clamp_fraction((a_lo - a_bed)/(e_lo - e_bed))
               dy_cu_bt(i, j) = dy_cu(i, j)*por_col
            else
               dy_cu_bt(i, j) = 0.0_wp
            end if
         end if
      end do

      ! ---- v-faces ----
      do concurrent(j=1:ny + 1, i=1:nx) &
         local(k, z_a, z_b, e_lo, e_hi, a_lo, a_hi, dz, e_bed, a_bed, por_col)
         if (j < 2 .or. j > ny) then
            do k = 1, nz
               por_v(i, j, k) = 1.0_wp
            end do
            dx_cv_bt(i, j) = dx_cv(i, j)
         else if (davg_v(i, j) >= mask_depth) then
            do k = 1, nz
               por_v(i, j, k) = 1.0_wp
            end do
            dx_cv_bt(i, j) = dx_cv(i, j)
         else if (dmax_v(i, j) <= dmin_v(i, j)) then
            ! Degenerate face — see the zonal twin.
            do k = 1, nz
               por_v(i, j, k) = 1.0_wp
            end do
            dx_cv_bt(i, j) = dx_cv(i, j)
         else
            z_a = b(i, j - 1)
            z_b = b(i, j)
            e_lo = porous_eta_face(z_a, z_b, interp)
            a_lo = porous_cum_area(dmin_v(i, j), dmax_v(i, j), davg_v(i, j), e_lo)
            e_bed = e_lo
            a_bed = a_lo
            do k = 1, nz
               z_a = z_a + h_layer(i, j - 1, k)
               z_b = z_b + h_layer(i, j, k)
               e_hi = porous_eta_face(z_a, z_b, interp)
               a_hi = porous_cum_area(dmin_v(i, j), dmax_v(i, j), davg_v(i, j), e_hi)
               dz = e_hi - e_lo
               if (dz > H_VANISHED) then
                  por_v(i, j, k) = clamp_fraction((a_hi - a_lo)/dz)
               else
                  por_v(i, j, k) = 0.0_wp
               end if
               e_lo = e_hi
               a_lo = a_hi
            end do
            ! Column-integrated open fraction — see the zonal twin.
            if (e_lo - e_bed > H_VANISHED) then
               por_col = clamp_fraction((a_lo - a_bed)/(e_lo - e_bed))
               dx_cv_bt(i, j) = dx_cv(i, j)*por_col
            else
               dx_cv_bt(i, j) = 0.0_wp
            end if
         end if
      end do
   end subroutine porous_update_face_areas

   pure function clamp_fraction(x) result(f)
      !! Clamp an open-area fraction into `[0,1]`, NaN-safely.
      !!
      !! The bracket is unreachable in exact arithmetic — `A` is monotone
      !! with slope in `[0,1]`, so `(A(hi)-A(lo))/(hi-lo)` already lies in
      !! `[0,1]` — but it is NOT unreachable in floating point: the
      !! above-`d_max` branch evaluates `A = eta - d_avg`, and the
      !! difference of two such values can exceed `hi - lo` by an ulp, so
      !! the `min` is what makes a fully submerged sill give back EXACTLY
      !! 1.  The explicit non-finite test is the repo's clamp rule
      !! (CLAUDE.md): under nvfortran's relaxed FP an `if/else` clamp
      !! lowers to a NaN-blind min/max select, which would launder a NaN
      !! into a plausible 0 or 1 instead of a visibly blocked face.
      real(wp), intent(in) :: x
      real(wp) :: f
      if (.not. ieee_is_finite(x)) then
         f = 0.0_wp
      else
         f = max(0.0_wp, min(1.0_wp, x))
      end if
   end function clamp_fraction

   pure subroutine closed_faces_update_bt_widths(nx, ny, nz, use_por, &
                                                 dy_cu, dx_cv, h_layer, &
                                                 por_u, por_v, open_u, open_v, &
                                                 dy_cu_bt, dx_cv_bt)
      !! Refresh the BAROTROPIC face widths from the LIVE layer
      !! thicknesses when `&vcoord_nml zfixed_closed_faces` is on.
      !!
      !! The barotropic substep transports on `ubt * FA * dy_cu_bt` with
      !! `FA = sum_k h_face` — the FULL column.  A closed layer carries no
      !! transport, so the BT solve has to see the OPEN depth of the face
      !! or it will hand the blocked transport straight back through the
      !! per-layer renormalisation.  Narrowing the WIDTH by the open
      !! fraction and reducing the DEPTH to the open depth are the same
      !! number here, so:
      !!
      !! ```
      !! dy_cu_bt(I,j) = dy_cu(I,j) * (sum_k h_face_k * por_k * open_k)
      !!                            / (sum_k h_face_k)
      !! ```
      !!
      !! WRITE, not multiply: when porous barriers are also on,
      !! `porous_update_face_areas` has just written
      !! `dy_cu_bt = dy_cu * por_col` with its own thickness-weighted mean
      !! `por_col`, and the expression above already contains that mean.
      !! Multiplying would count the porous fraction twice.  This routine
      !! therefore SUPERSEDES the porous write and must run after it —
      !! which is exactly the order `ocean_porous_refresh` calls them in.
      !!
      !! Per OUTER step (MOM6's porous cadence), from the live `h`, not
      !! frozen at `eta = 0`: under `z_fixed` the `eta = 0` value is only
      !! `O(eta/H) ~ 1E-4` off because eta lands in the first LIVE layer,
      !! but "only a small error" is not a reason to carry one.
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: use_por
         !! Porous barriers also active.  `.false.` => `por_u`/`por_v`
         !! are never indexed, and the caller must hand over a full-size,
         !! device-present stand-in rather than the `(1,1,1)` porous
         !! placeholder: nvfortran builds the `do concurrent` data clause
         !! from the LOOP BOUNDS, not from the descriptor, so a
         !! placeholder aborts under `mem:separate` ("variable in data
         !! clause is partially present") even though the branch that
         !! indexes it is never taken.  `ocean_porous_refresh` passes
         !! `open_u`/`open_v` themselves — right shape, already mapped,
         !! `intent(in)` on both dummies so the double association is not
         !! aliasing.  (Found by the GPU build; both CPU builds were
         !! silently happy.)
      real(wp), intent(in) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: por_u(nx + 1, ny, nz), por_v(nx, ny + 1, nz)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz), open_v(nx, ny + 1, nz)
      real(wp), intent(inout) :: dy_cu_bt(nx + 1, ny), dx_cv_bt(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: h_face, sum_all, sum_open, wk

      do concurrent(j=1:ny, i=2:nx) local(k, h_face, sum_all, sum_open, wk)
         sum_all = 0.0_wp
         sum_open = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            wk = open_u(i, j, k)
            if (use_por) wk = wk*por_u(i, j, k)
            sum_all = sum_all + h_face
            sum_open = sum_open + h_face*wk
         end do
         if (sum_all > 0.0_wp) then
            dy_cu_bt(i, j) = dy_cu(i, j)*(sum_open/sum_all)
         else
            dy_cu_bt(i, j) = 0.0_wp
         end if
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, h_face, sum_all, sum_open, wk)
         sum_all = 0.0_wp
         sum_open = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            wk = open_v(i, j, k)
            if (use_por) wk = wk*por_v(i, j, k)
            sum_all = sum_all + h_face
            sum_open = sum_open + h_face*wk
         end do
         if (sum_all > 0.0_wp) then
            dx_cv_bt(i, j) = dx_cv(i, j)*(sum_open/sum_all)
         else
            dx_cv_bt(i, j) = 0.0_wp
         end if
      end do
   end subroutine closed_faces_update_bt_widths

end module rdb_ocean_porous
