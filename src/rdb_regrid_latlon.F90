!! Regular lat-lon source -> model T-grid regridding (in-model weights).
module rdb_regrid_latlon
   !! Grid-agnostic regridding library: builds bilinear or first-order
   !! conservative weights from a regular lat-lon SOURCE grid onto the
   !! model's T-grid, stored in an ELL (padded-stencil) layout, and applies
   !! them with one device kernel. No namelist knob, no engine wiring here
   !! -- bit-identity by construction; the engine-side selector knobs
   !! (`states_regrid`/`fluxes_regrid`/...) are added by a later PR.
   !!
   !! References: Jones, P. W. (1999), "First- and Second-Order Conservative
   !! Remapping Schemes for Grids in Spherical Coordinates", Monthly Weather
   !! Review, 127, 2204-2210; Sutherland, I. E. and Hodgman, G. W. (1974),
   !! "Reentrant Polygon Clipping", Communications of the ACM, 17(1), 32-42;
   !! Murray, R. J. (1996), "Explicit Generation of Orthogonal Grids for
   !! Ocean Models", Journal of Computational Physics, 126, 251-273.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, PI, TWO_PI, DEG2RAD
   implicit none
   private

   public :: regrid_weights_t
   public :: REGRID_NMAX_CAP
   public :: regrid_bnds_to_edges
   public :: regrid_weights_bilinear
   public :: regrid_weights_conservative
   public :: regrid_apply_2d_impl
   public :: regrid_apply_2d
   public :: regrid_window

   integer, parameter :: REGRID_NMAX_CAP = 64
      !! Fail-loud ceiling on the per-cell stencil size `nmax`: a 4-degree
      !! target on a TL319-like (~0.5625-degree) source needs ~64 source
      !! cells inside one conservative polygon's bounding box; a 1-degree
      !! target needs 9-12. Exceeding it is always `ierr`, never a silent
      !! truncation of the accumulated area.

   integer, parameter :: REGRID_POLY_MAXV = 24
      !! Max vertex count carried through the Sutherland-Hodgman clip
      !! (quad or the 7-vertex pole-cap polygon, clipped by one rectangle:
      !! at most 7+4 = 11 vertices in practice).
   real(wp), parameter :: POLE_TOL_DEG = 1.0e-9_wp
      !! A target corner within this of +-90 deg is treated as AT the pole.

   type :: regrid_weights_t
      !! ELL/padded regridding weights, physical target cells only.
      integer :: nx = 0, ny = 0
         !! Target tile size, PHYSICAL cells only (no ghosts).
      integer :: nmax = 0
         !! Stencil width (padded): 4 for bilinear; up to
         !! `REGRID_NMAX_CAP` for conservative.
      integer :: nsx = 0
         !! Source window width in longitude -- ALWAYS the full lon
         !! circle (`i_lo = 1` always; `idx_i` is a GLOBAL column index).
      integer :: j_lo = 0, nsy = 0
         !! Source window in latitude: GLOBAL source rows
         !! `[j_lo, j_lo + nsy - 1]`, clamped to `[1, nlat]`.
      integer, allocatable :: idx_i(:, :, :), idx_j(:, :, :)
         !! `(nx, ny, nmax)`. `idx_i` is a GLOBAL source column
         !! (`1..nsx`); `idx_j` is WINDOW-LOCAL
         !! (`idx_j = j_global - j_lo + 1`), so it indexes straight into
         !! a `(nsx, nsy)` windowed source array handed in by the caller.
      real(wp), allocatable :: w(:, :, :)
         !! `(nx, ny, nmax)` interpolation/area weights. Padding slots
         !! (beyond a cell's actual stencil count): `w = 0`,
         !! `idx_i = idx_j = 1` (always-valid indices, zero contribution).
      real(wp), allocatable :: area_t(:, :)
         !! Host-only diagnostic (conservative only): target-cell area in
         !! (lon [rad], sin(lat)) units -- the sum of the clipped source
         !! fragment areas, which equals the full target-cell area when
         !! coverage is complete.
      integer :: n_fallback = 0
         !! Host-only diagnostic: count of degenerate conservative target
         !! cells that fell back to the bilinear-at-centre 4-point stencil.
      logical :: is_init = .false.
      logical :: on_device = .false.
   contains
      procedure, non_overridable :: enter_data => regrid_weights_enter_data
      procedure, non_overridable :: exit_data => regrid_weights_exit_data
      procedure, non_overridable :: bytes => regrid_weights_bytes
      procedure, non_overridable :: destroy => regrid_weights_destroy
   end type regrid_weights_t

contains

   ! =====================================================================
   ! Type-bound: device mapping, footprint, teardown.
   ! =====================================================================

   subroutine regrid_weights_enter_data(this)
      class(regrid_weights_t), intent(inout) :: this
      !$acc enter data copyin(this%idx_i, this%idx_j, this%w)
      this%on_device = .true.
   end subroutine regrid_weights_enter_data

   subroutine regrid_weights_exit_data(this)
      class(regrid_weights_t), intent(inout) :: this
      !$acc exit data delete(this%idx_i, this%idx_j, this%w)
      this%on_device = .false.
   end subroutine regrid_weights_exit_data

   pure function regrid_weights_bytes(this) result(nbytes)
      class(regrid_weights_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (allocated(this%idx_i)) nbytes = nbytes + int(size(this%idx_i), int64)*4_int64
      if (allocated(this%idx_j)) nbytes = nbytes + int(size(this%idx_j), int64)*4_int64
      if (allocated(this%w)) nbytes = nbytes + int(size(this%w), int64)*8_int64
      if (allocated(this%area_t)) nbytes = nbytes + int(size(this%area_t), int64)*8_int64
   end function regrid_weights_bytes

   subroutine regrid_weights_destroy(this)
      class(regrid_weights_t), intent(inout) :: this
      if (this%on_device) call this%exit_data()
      if (allocated(this%idx_i)) deallocate (this%idx_i)
      if (allocated(this%idx_j)) deallocate (this%idx_j)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%area_t)) deallocate (this%area_t)
      this%nx = 0
      this%ny = 0
      this%nmax = 0
      this%nsx = 0
      this%j_lo = 0
      this%nsy = 0
      this%n_fallback = 0
      this%is_init = .false.
   end subroutine regrid_weights_destroy

   ! =====================================================================
   ! CF lon_bnds/lat_bnds -> monotone edge array.
   ! =====================================================================

   pure subroutine regrid_bnds_to_edges(n, bnds, edges)
      !! CF convention: `bnds(1,k)`/`bnds(2,k)` are cell k's two edges,
      !! adjacent cells sharing an edge (`bnds(2,k) == bnds(1,k+1)`).
      !! `edges(1) = bnds(1,1)`, `edges(n+1) = bnds(2,n)`.
      integer, intent(in) :: n
      real(wp), intent(in) :: bnds(2, n)
      real(wp), intent(out) :: edges(n + 1)
      integer :: k
      edges(1) = bnds(1, 1)
      do k = 1, n
         edges(k + 1) = bnds(2, k)
      end do
   end subroutine regrid_bnds_to_edges

   ! =====================================================================
   ! Shared host-side locate helpers (bisection on a monotone-increasing
   ! array; the lon/lat axis orientation is normalised at each call site).
   ! =====================================================================

   pure subroutine bisect_bracket(n, x_inc, x, i0)
      !! `i0` in `[1, n-1]` with `x_inc(i0) <= x <= x_inc(i0+1)`; clamps
      !! to the end interval when `x` is outside the array's range.
      integer, intent(in) :: n
      real(wp), intent(in) :: x_inc(n)
      real(wp), intent(in) :: x
      integer, intent(out) :: i0
      integer :: lo, hi, mid
      if (x <= x_inc(1)) then
         i0 = 1
      else if (x >= x_inc(n)) then
         i0 = n - 1
      else
         lo = 1
         hi = n
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (x_inc(mid) <= x) then
               lo = mid
            else
               hi = mid
            end if
         end do
         i0 = lo
      end if
   end subroutine bisect_bracket

   pure subroutine lon_bisect(nlon, lon_c, x, i0, fx)
      !! Locates the longitude cell bracketing `x` (already unwrapped
      !! into `[0, 360)` relative to `lon_c(1)`), by bisection on `lon_c`
      !! -- works for a non-uniform axis. `i0` in `[1, nlon]`, `fx` the
      !! fractional position within the bracket `[i0, i0+1]` (periodic).
      integer, intent(in) :: nlon
      real(wp), intent(in) :: lon_c(nlon)
      real(wp), intent(in) :: x
      integer, intent(out) :: i0
      real(wp), intent(out) :: fx
      real(wp) :: lon_e_rel(nlon + 1)
      integer :: k
      do k = 1, nlon
         lon_e_rel(k) = modulo(lon_c(k) - lon_c(1), 360.0_wp)
      end do
      lon_e_rel(nlon + 1) = 360.0_wp
      call bisect_bracket(nlon + 1, lon_e_rel, x, i0)
      if (i0 > nlon) i0 = nlon
      fx = (x - lon_e_rel(i0))/(lon_e_rel(i0 + 1) - lon_e_rel(i0))
   end subroutine lon_bisect

   pure function lon_is_uniform(nlon, lon_c, dlon) result(uniform)
      integer, intent(in) :: nlon
      real(wp), intent(in) :: lon_c(nlon)
      real(wp), intent(in) :: dlon
      logical :: uniform
      real(wp) :: maxdev
      integer :: k
      maxdev = 0.0_wp
      do k = 1, nlon - 1
         maxdev = max(maxdev, abs((lon_c(k + 1) - lon_c(k)) - dlon))
      end do
      uniform = maxdev < 1.0e-9_wp*dlon
   end function lon_is_uniform

   ! =====================================================================
   ! Bilinear weight builder.
   ! =====================================================================

   pure subroutine regrid_weights_bilinear(wts, nlon, nlat, lon_c, lat_c, &
                                           nx_phys, ny_phys, nghost, n1, n2, &
                                           geolonT, geolatT, ierr)
      !! `n1,n2 = nx_total,ny_total`; reads physical T cells
      !! `(nghost+1 : nghost+nx_phys, nghost+1 : nghost+ny_phys)` only.
      integer, intent(in) :: nlon, nlat, nx_phys, ny_phys, nghost, n1, n2
      real(wp), intent(in) :: lon_c(nlon), lat_c(nlat)
      real(wp), intent(in) :: geolonT(n1, n2), geolatT(n1, n2)
      type(regrid_weights_t), intent(out) :: wts
      integer, intent(out) :: ierr

      integer, parameter :: NMAX_BILINEAR = 4
      logical :: lat_flip, lon_uniform
      real(wp) :: dlon, lon_t, lat_t, x, fx, fy
      real(wp) :: lat_inc(nlat)
      integer :: i, j, ip, jp, i0, i1, j0f, j1f, j0, j1
      integer :: jmin, jmax

      ierr = 0
      dlon = 360.0_wp/real(nlon, wp)
      lon_uniform = lon_is_uniform(nlon, lon_c, dlon)

      lat_flip = lat_c(1) > lat_c(nlat)
      if (lat_flip) then
         lat_inc = lat_c(nlat:1:-1)
      else
         lat_inc = lat_c
      end if

      wts%nx = nx_phys
      wts%ny = ny_phys
      wts%nmax = NMAX_BILINEAR
      wts%nsx = nlon
      allocate (wts%idx_i(nx_phys, ny_phys, NMAX_BILINEAR), source=1)
      allocate (wts%idx_j(nx_phys, ny_phys, NMAX_BILINEAR), source=1)
      allocate (wts%w(nx_phys, ny_phys, NMAX_BILINEAR), source=0.0_wp)

      jmin = nlat
      jmax = 1

      do j = 1, ny_phys
         jp = j + nghost
         do i = 1, nx_phys
            ip = i + nghost
            lon_t = geolonT(ip, jp)
            lat_t = geolatT(ip, jp)

            x = modulo(lon_t - lon_c(1), 360.0_wp)
            if (lon_uniform) then
               i0 = min(int(x/dlon) + 1, nlon)
               fx = x/dlon - real(i0 - 1, wp)
            else
               call lon_bisect(nlon, lon_c, x, i0, fx)
            end if
            i1 = mod(i0, nlon) + 1

            if (lat_t <= lat_inc(1)) then
               j0f = 1
               j1f = 1
               fy = 0.0_wp
            else if (lat_t >= lat_inc(nlat)) then
               j0f = nlat
               j1f = nlat
               fy = 0.0_wp
            else
               call bisect_bracket(nlat, lat_inc, lat_t, j0f)
               j1f = j0f + 1
               fy = (lat_t - lat_inc(j0f))/(lat_inc(j1f) - lat_inc(j0f))
            end if

            if (lat_flip) then
               j0 = nlat + 1 - j0f
               j1 = nlat + 1 - j1f
            else
               j0 = j0f
               j1 = j1f
            end if

            wts%idx_i(i, j, 1) = i0
            wts%idx_j(i, j, 1) = j0
            wts%idx_i(i, j, 2) = i1
            wts%idx_j(i, j, 2) = j0
            wts%idx_i(i, j, 3) = i0
            wts%idx_j(i, j, 3) = j1
            wts%idx_i(i, j, 4) = i1
            wts%idx_j(i, j, 4) = j1
            wts%w(i, j, 1) = (1.0_wp - fx)*(1.0_wp - fy)
            wts%w(i, j, 2) = fx*(1.0_wp - fy)
            wts%w(i, j, 3) = (1.0_wp - fx)*fy
            wts%w(i, j, 4) = fx*fy

            jmin = min(jmin, min(j0, j1))
            jmax = max(jmax, max(j0, j1))
         end do
      end do

      wts%j_lo = max(1, jmin - 1)
      wts%nsy = min(nlat, jmax + 1) - wts%j_lo + 1
      wts%idx_j = wts%idx_j - wts%j_lo + 1

      wts%is_init = .true.
   end subroutine regrid_weights_bilinear

   pure subroutine bilinear_point_from_edges(nlon, nlat, lon_e, lat_e, lon_t, lat_t, &
                                             ii, jj, ww)
      !! Fallback bilinear-at-centre 4-point stencil for a degenerate
      !! conservative cell: cell centres taken as edge midpoints
      !! (`lon_c = 0.5*(lon_e(i)+lon_e(i+1))`, likewise for lat).
      integer, intent(in) :: nlon, nlat
      real(wp), intent(in) :: lon_e(nlon + 1), lat_e(nlat + 1)
      real(wp), intent(in) :: lon_t, lat_t
      integer, intent(out) :: ii(4), jj(4)
      real(wp), intent(out) :: ww(4)

      real(wp) :: lon_c(nlon), lat_c(nlat), lat_inc(nlat)
      integer :: k, i0, i1, j0, j1
      real(wp) :: x, fx, fy
      logical :: lat_flip

      do k = 1, nlon
         lon_c(k) = 0.5_wp*(lon_e(k) + lon_e(k + 1))
      end do
      do k = 1, nlat
         lat_c(k) = 0.5_wp*(lat_e(k) + lat_e(k + 1))
      end do

      lat_flip = lat_c(1) > lat_c(nlat)
      if (lat_flip) then
         lat_inc = lat_c(nlat:1:-1)
      else
         lat_inc = lat_c
      end if

      x = modulo(lon_t - lon_c(1), 360.0_wp)
      call lon_bisect(nlon, lon_c, x, i0, fx)
      i1 = mod(i0, nlon) + 1

      if (lat_t <= lat_inc(1)) then
         j0 = 1
         j1 = 1
         fy = 0.0_wp
      else if (lat_t >= lat_inc(nlat)) then
         j0 = nlat
         j1 = nlat
         fy = 0.0_wp
      else
         call bisect_bracket(nlat, lat_inc, lat_t, j0)
         j1 = j0 + 1
         fy = (lat_t - lat_inc(j0))/(lat_inc(j1) - lat_inc(j0))
      end if
      if (lat_flip) then
         j0 = nlat + 1 - j0
         j1 = nlat + 1 - j1
      end if

      ii = [i0, i1, i0, i1]
      jj = [j0, j0, j1, j1]
      ww(1) = (1.0_wp - fx)*(1.0_wp - fy)
      ww(2) = fx*(1.0_wp - fy)
      ww(3) = (1.0_wp - fx)*fy
      ww(4) = fx*fy
   end subroutine bilinear_point_from_edges

   ! =====================================================================
   ! Sutherland-Hodgman polygon clip (subject clipped by a convex region)
   ! and shoelace area, used by the conservative weight builder.
   ! =====================================================================

   pure subroutine clip_halfplane(n_in, px, py, a, b, c, n_out, qx, qy)
      !! Clips polygon `(px,py)` (length `n_in`) to the half-plane
      !! `a*x + b*y <= c`, writing the result into `(qx,qy)` / `n_out`.
      integer, intent(in) :: n_in
      real(wp), intent(in) :: px(REGRID_POLY_MAXV), py(REGRID_POLY_MAXV)
      real(wp), intent(in) :: a, b, c
      integer, intent(out) :: n_out
      real(wp), intent(out) :: qx(REGRID_POLY_MAXV), qy(REGRID_POLY_MAXV)
      integer :: k
      real(wp) :: d_cur, d_prev, t, xcur, ycur, xprev, yprev

      n_out = 0
      if (n_in < 1) return
      xprev = px(n_in)
      yprev = py(n_in)
      d_prev = a*xprev + b*yprev - c
      do k = 1, n_in
         xcur = px(k)
         ycur = py(k)
         d_cur = a*xcur + b*ycur - c
         if (d_cur <= 0.0_wp) then
            if (d_prev > 0.0_wp .and. n_out < REGRID_POLY_MAXV) then
               t = d_prev/(d_prev - d_cur)
               n_out = n_out + 1
               qx(n_out) = xprev + t*(xcur - xprev)
               qy(n_out) = yprev + t*(ycur - yprev)
            end if
            if (n_out < REGRID_POLY_MAXV) then
               n_out = n_out + 1
               qx(n_out) = xcur
               qy(n_out) = ycur
            end if
         else if (d_prev <= 0.0_wp .and. n_out < REGRID_POLY_MAXV) then
            t = d_prev/(d_prev - d_cur)
            n_out = n_out + 1
            qx(n_out) = xprev + t*(xcur - xprev)
            qy(n_out) = yprev + t*(ycur - yprev)
         end if
         xprev = xcur
         yprev = ycur
         d_prev = d_cur
      end do
   end subroutine clip_halfplane

   pure subroutine clip_box(n_in, px, py, xmin, xmax, ymin, ymax, n_out, qx, qy)
      !! Clips polygon `(px,py)` to the axis-aligned box
      !! `[xmin,xmax] x [ymin,ymax]` by composing four half-plane clips.
      integer, intent(in) :: n_in
      real(wp), intent(in) :: px(REGRID_POLY_MAXV), py(REGRID_POLY_MAXV)
      real(wp), intent(in) :: xmin, xmax, ymin, ymax
      integer, intent(out) :: n_out
      real(wp), intent(out) :: qx(REGRID_POLY_MAXV), qy(REGRID_POLY_MAXV)
      real(wp) :: ax(REGRID_POLY_MAXV), ay(REGRID_POLY_MAXV)
      real(wp) :: bx(REGRID_POLY_MAXV), by(REGRID_POLY_MAXV)
      integer :: na, nb

      call clip_halfplane(n_in, px, py, 1.0_wp, 0.0_wp, xmax, na, ax, ay)
      call clip_halfplane(na, ax, ay, -1.0_wp, 0.0_wp, -xmin, nb, bx, by)
      call clip_halfplane(nb, bx, by, 0.0_wp, 1.0_wp, ymax, na, ax, ay)
      call clip_halfplane(na, ax, ay, 0.0_wp, -1.0_wp, -ymin, n_out, qx, qy)
   end subroutine clip_box

   pure function segments_cross(ax, ay, bx, by, cx, cy, dx, dy) result(cross)
      !! True when the open segments AB and CD properly cross (each pair of
      !! endpoints strictly on opposite sides of the other segment's line).
      !! The relative tolerance keeps touching / collinear / coincident
      !! endpoints (a collapsed cell edge) from counting as a crossing.
      real(wp), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy
      logical :: cross
      real(wp) :: o1, o2, o3, o4, tol

      tol = 1.0e-12_wp*max(abs(bx - ax) + abs(by - ay), abs(dx - cx) + abs(dy - cy))**2
      o1 = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
      o2 = (bx - ax)*(dy - ay) - (by - ay)*(dx - ax)
      o3 = (dx - cx)*(ay - cy) - (dy - cy)*(ax - cx)
      o4 = (dx - cx)*(by - cy) - (dy - cy)*(bx - cx)
      cross = ((o1 > tol .and. o2 < -tol) .or. (o1 < -tol .and. o2 > tol)) .and. &
              ((o3 > tol .and. o4 < -tol) .or. (o3 < -tol .and. o4 > tol))
   end function segments_cross

   pure function poly_area2(n, px, py) result(a2)
      !! Twice the SIGNED shoelace area (CCW positive).
      integer, intent(in) :: n
      real(wp), intent(in) :: px(REGRID_POLY_MAXV), py(REGRID_POLY_MAXV)
      real(wp) :: a2
      integer :: k, kn
      a2 = 0.0_wp
      if (n < 3) return
      do k = 1, n
         kn = k + 1
         if (kn > n) kn = 1
         a2 = a2 + px(k)*py(kn) - px(kn)*py(k)
      end do
   end function poly_area2

   ! =====================================================================
   ! Conservative (first-order, Jones 1999) weight builder.
   ! =====================================================================

   pure subroutine regrid_weights_conservative(wts, nlon, nlat, lon_e, lat_e, &
                                               nx_phys, ny_phys, nghost, n1, n2, &
                                               geolonBu, geolatBu, geolonT, geolatT, &
                                               ierr, nmax_cap_override)
      !! Per-target-cell area-weighted overlap in the (lon [rad],
      !! mu = sin(lat)) plane, where a lat-lon source box's spherical
      !! area is exactly `R^2 . d(lon) . d(mu)` -- so clipping is planar
      !! there. Bilinear-at-centre fallback for degenerate quads.
      integer, intent(in) :: nlon, nlat, nx_phys, ny_phys, nghost, n1, n2
      real(wp), intent(in) :: lon_e(nlon + 1), lat_e(nlat + 1)
      real(wp), intent(in) :: geolonBu(n1 + 1, n2 + 1), geolatBu(n1 + 1, n2 + 1)
      real(wp), intent(in) :: geolonT(n1, n2), geolatT(n1, n2)
      type(regrid_weights_t), intent(out) :: wts
      integer, intent(out) :: ierr
      integer, intent(in), optional :: nmax_cap_override

      integer, parameter :: LOCAL_MAXCAND = 512

      integer :: cap, i, j, ip, jp, k, kc
      integer :: ncur, nmax_used, jmin_g, jmax_g
      logical :: lat_flip, is_cap, degenerate
      logical :: at_pole(4)
      integer :: kref
      real(wp) :: lat_e_rad(nlat + 1), mu_e_inc(nlat + 1)
      real(wp) :: lon_e_rad(nlon + 1)
      real(wp) :: c_lon0(4), c_lat0(4), lam(4), mu(4)
      real(wp) :: sorted_lon(4), sorted_mu(4)
      real(wp) :: px(REGRID_POLY_MAXV), py(REGRID_POLY_MAXV)
      real(wp) :: qx(REGRID_POLY_MAXV), qy(REGRID_POLY_MAXV)
      integer :: npoly
      real(wp) :: span, lam1, mu1
      real(wp) :: lam_min, lam_max, mu_min, mu_max
      real(wp) :: cand_area(LOCAL_MAXCAND)
      integer :: cand_i(LOCAL_MAXCAND), cand_j(LOCAL_MAXCAND)
      integer :: n_cand, i_cur, jf, jf_lo, jf_hi, j_orig
      real(wp) :: box_lo, box_hi, mu_lo, mu_hi, total_area
      real(wp) :: lon_t, lat_t
      integer :: ii4(4), jj4(4)
      real(wp) :: ww4(4)
      integer, allocatable :: idx_i_tmp(:, :, :), idx_j_tmp(:, :, :)
      real(wp), allocatable :: w_tmp(:, :, :)
      integer, allocatable :: ncell(:, :)

      ierr = 0
      if (present(nmax_cap_override)) then
         cap = nmax_cap_override
      else
         cap = REGRID_NMAX_CAP
      end if

      do k = 1, nlon + 1
         lon_e_rad(k) = lon_e(k)*DEG2RAD
      end do
      do k = 1, nlat + 1
         lat_e_rad(k) = lat_e(k)*DEG2RAD
      end do

      lat_flip = lat_e_rad(1) > lat_e_rad(nlat + 1)
      if (lat_flip) then
         mu_e_inc = sin(lat_e_rad(nlat + 1:1:-1))
      else
         mu_e_inc = sin(lat_e_rad)
      end if

      wts%nx = nx_phys
      wts%ny = ny_phys
      wts%nsx = nlon
      allocate (idx_i_tmp(nx_phys, ny_phys, LOCAL_MAXCAND), source=1)
      allocate (idx_j_tmp(nx_phys, ny_phys, LOCAL_MAXCAND), source=1)
      allocate (w_tmp(nx_phys, ny_phys, LOCAL_MAXCAND), source=0.0_wp)
      allocate (ncell(nx_phys, ny_phys), source=0)
      allocate (wts%area_t(nx_phys, ny_phys), source=0.0_wp)
      wts%n_fallback = 0
      nmax_used = 0
      jmin_g = nlat
      jmax_g = 1

      do j = 1, ny_phys
         jp = j + nghost
         do i = 1, nx_phys
            ip = i + nghost

            ! ---- build the (lam, mu) target polygon ----
            c_lon0(1) = geolonBu(ip, jp)
            c_lat0(1) = geolatBu(ip, jp)
            c_lon0(2) = geolonBu(ip + 1, jp)
            c_lat0(2) = geolatBu(ip + 1, jp)
            c_lon0(3) = geolonBu(ip + 1, jp + 1)
            c_lat0(3) = geolatBu(ip + 1, jp + 1)
            c_lon0(4) = geolonBu(ip, jp + 1)
            c_lat0(4) = geolatBu(ip, jp + 1)
            ! A corner AT a geographic pole has no meaningful longitude:
            ! unwrap around the first non-pole corner and leave pole
            ! corners out of the span (they are expanded below).
            do k = 1, 4
               at_pole(k) = abs(c_lat0(k)) >= 90.0_wp - POLE_TOL_DEG
            end do
            kref = 1
            do k = 4, 1, -1
               if (.not. at_pole(k)) kref = k
            end do
            do k = 1, 4
               c_lon0(k) = c_lon0(kref) + (modulo(c_lon0(k) - c_lon0(kref) + 180.0_wp, &
                                                  360.0_wp) - 180.0_wp)
            end do
            do k = 1, 4
               lam(k) = c_lon0(k)*DEG2RAD
               mu(k) = sin(c_lat0(k)*DEG2RAD)
            end do
            span = maxval(c_lon0, mask=.not. at_pole) - minval(c_lon0, mask=.not. at_pole)
            is_cap = span > 180.0_wp

            degenerate = .false.
            if (.not. is_cap) then
               ! In (lam, mu) a pole is the whole line mu = +-1, so a pole
               ! corner becomes two points on it at its neighbours'
               ! longitudes (meridian edges stay vertical, so the cells
               ! sharing the pole still tile the plane exactly).
               npoly = 0
               do k = 1, 4
                  if (at_pole(k)) then
                     npoly = npoly + 1
                     px(npoly) = lam(mod(k + 2, 4) + 1)
                     py(npoly) = sign(1.0_wp, c_lat0(k))
                     npoly = npoly + 1
                     px(npoly) = lam(mod(k, 4) + 1)
                     py(npoly) = sign(1.0_wp, c_lat0(k))
                  else
                     npoly = npoly + 1
                     px(npoly) = lam(k)
                     py(npoly) = mu(k)
                  end if
               end do
               ! Only a SELF-INTERSECTING (bow-tie) quad is degenerate:
               ! Sutherland-Hodgman against a convex box is exact for any
               ! simple polygon, concave ones included (cells touching the
               ! pole), and a collapsed edge (bipole cells) is harmless.
               if (npoly == 4) degenerate = segments_cross(px(1), py(1), px(2), py(2), &
                                                           px(3), py(3), px(4), py(4)) .or. &
                                            segments_cross(px(2), py(2), px(3), py(3), &
                                                           px(4), py(4), px(1), py(1))
            else
               sorted_lon = c_lon0
               sorted_mu = mu
               call sort4(sorted_lon, sorted_mu)
               lam1 = sorted_lon(1)*DEG2RAD
               mu1 = sorted_mu(1)
               npoly = 7
               do k = 1, 4
                  px(k) = sorted_lon(k)*DEG2RAD
                  py(k) = sorted_mu(k)
               end do
               px(5) = lam1 + TWO_PI
               py(5) = mu1
               px(6) = lam1 + TWO_PI
               py(6) = 1.0_wp
               px(7) = lam1
               py(7) = 1.0_wp
            end if

            lon_t = geolonT(ip, jp)
            lat_t = geolatT(ip, jp)

            if (degenerate) then
               call bilinear_point_from_edges(nlon, nlat, lon_e, lat_e, lon_t, lat_t, &
                                              ii4, jj4, ww4)
               wts%n_fallback = wts%n_fallback + 1
               do k = 1, 4
                  idx_i_tmp(i, j, k) = ii4(k)
                  idx_j_tmp(i, j, k) = jj4(k)
                  w_tmp(i, j, k) = ww4(k)
               end do
               ncell(i, j) = 4
               wts%area_t(i, j) = 0.0_wp
               jmin_g = min(jmin_g, minval(jj4))
               jmax_g = max(jmax_g, maxval(jj4))
               nmax_used = max(nmax_used, 4)
               cycle
            end if

            lam_min = minval(px(1:npoly))
            lam_max = maxval(px(1:npoly))
            mu_min = minval(py(1:npoly))
            mu_max = maxval(py(1:npoly))

            ! ---- j (latitude) candidate range: contiguous in the
            ! increasing-mu space, mapped back to the original row index.
            call bisect_bracket(nlat + 1, mu_e_inc, mu_min, jf_lo)
            call bisect_bracket(nlat + 1, mu_e_inc, mu_max, jf_hi)
            jf_lo = max(1, jf_lo)
            jf_hi = min(nlat, jf_hi)

            ! ---- i (longitude) candidate walk: start at the cell
            ! bracketing lam_min in the SAME absolute (unwrapped) frame
            ! as the target polygon, then walk east until lam_max is
            ! covered (periodic wrap in i, incrementing a 2*pi shift).
            n_cand = 0
            call lon_cell_shifted(nlon, lon_e_rad, lam_min, i_cur, box_lo, box_hi)
            do kc = 1, nlon + 2
               do jf = jf_lo, jf_hi
                  j_orig = merge(nlat + 1 - jf, jf, lat_flip)
                  mu_lo = mu_e_inc(jf)
                  mu_hi = mu_e_inc(jf + 1)
                  call clip_box(npoly, px, py, box_lo, box_hi, mu_lo, mu_hi, &
                                ncur, qx, qy)
                  if (ncur >= 3) then
                     if (n_cand < LOCAL_MAXCAND) then
                        n_cand = n_cand + 1
                        cand_i(n_cand) = i_cur
                        cand_j(n_cand) = j_orig
                        cand_area(n_cand) = 0.5_wp*abs(poly_area2(ncur, qx, qy))
                     end if
                  end if
               end do
               if (box_hi >= lam_max) exit
               i_cur = i_cur + 1
               if (i_cur > nlon) then
                  i_cur = 1
                  box_lo = box_hi
                  box_hi = box_lo + (lon_e_rad(2) - lon_e_rad(1))
               else
                  box_lo = box_hi
                  box_hi = box_lo + (lon_e_rad(i_cur + 1) - lon_e_rad(i_cur))
               end if
            end do

            total_area = 0.0_wp
            do k = 1, n_cand
               total_area = total_area + cand_area(k)
            end do

            if (n_cand == 0 .or. total_area < 1.0e-14_wp*max(1.0e-30_wp, &
                                                             (lam_max - lam_min)* &
                                                             max(1.0e-30_wp, mu_max - mu_min))) &
               then
               call bilinear_point_from_edges(nlon, nlat, lon_e, lat_e, lon_t, lat_t, &
                                              ii4, jj4, ww4)
               wts%n_fallback = wts%n_fallback + 1
               do k = 1, 4
                  idx_i_tmp(i, j, k) = ii4(k)
                  idx_j_tmp(i, j, k) = jj4(k)
                  w_tmp(i, j, k) = ww4(k)
               end do
               ncell(i, j) = 4
               wts%area_t(i, j) = 0.0_wp
               jmin_g = min(jmin_g, minval(jj4))
               jmax_g = max(jmax_g, maxval(jj4))
               nmax_used = max(nmax_used, 4)
               cycle
            end if

            if (n_cand > cap) ierr = 1
            ncell(i, j) = min(n_cand, LOCAL_MAXCAND)
            do k = 1, ncell(i, j)
               idx_i_tmp(i, j, k) = cand_i(k)
               idx_j_tmp(i, j, k) = cand_j(k)
               w_tmp(i, j, k) = cand_area(k)/total_area
               jmin_g = min(jmin_g, cand_j(k))
               jmax_g = max(jmax_g, cand_j(k))
            end do
            wts%area_t(i, j) = total_area
            nmax_used = max(nmax_used, ncell(i, j))
         end do
      end do

      wts%nmax = nmax_used
      wts%j_lo = max(1, jmin_g - 1)
      wts%nsy = min(nlat, jmax_g + 1) - wts%j_lo + 1

      allocate (wts%idx_i(nx_phys, ny_phys, wts%nmax), source=1)
      allocate (wts%idx_j(nx_phys, ny_phys, wts%nmax), source=1)
      allocate (wts%w(nx_phys, ny_phys, wts%nmax), source=0.0_wp)
      do j = 1, ny_phys
         do i = 1, nx_phys
            do k = 1, ncell(i, j)
               wts%idx_i(i, j, k) = idx_i_tmp(i, j, k)
               wts%idx_j(i, j, k) = idx_j_tmp(i, j, k) - wts%j_lo + 1
               wts%w(i, j, k) = w_tmp(i, j, k)
            end do
         end do
      end do

      wts%is_init = .true.
   end subroutine regrid_weights_conservative

   pure subroutine lon_cell_shifted(nlon, lon_e_rad, lam, i0, box_lo, box_hi)
      !! Locates the source longitude cell bracketing `lam` (an arbitrary
      !! real, possibly outside `lon_e_rad`'s own `[.,.+2*pi)` range) in
      !! `lam`'s OWN absolute frame: returns the GLOBAL column `i0` and
      !! that cell's edges shifted into the same frame as `lam`.
      integer, intent(in) :: nlon
      real(wp), intent(in) :: lon_e_rad(nlon + 1)
      real(wp), intent(in) :: lam
      integer, intent(out) :: i0
      real(wp), intent(out) :: box_lo, box_hi
      real(wp) :: x, raw, shift
      real(wp) :: lon_e_rel(nlon + 1)
      integer :: k

      do k = 1, nlon
         lon_e_rel(k) = modulo(lon_e_rad(k) - lon_e_rad(1), TWO_PI)
      end do
      lon_e_rel(nlon + 1) = TWO_PI
      raw = lam - lon_e_rad(1)
      x = modulo(raw, TWO_PI)
      shift = (raw - x)/TWO_PI

      call bisect_bracket(nlon + 1, lon_e_rel, x, i0)
      if (i0 > nlon) i0 = nlon
      box_lo = lon_e_rad(1) + shift*TWO_PI + lon_e_rel(i0)
      box_hi = lon_e_rad(1) + shift*TWO_PI + lon_e_rel(i0 + 1)
   end subroutine lon_cell_shifted

   pure subroutine sort4(a, b)
      !! Insertion sort of 4 (key, paired-value) pairs by key ascending.
      real(wp), intent(inout) :: a(4), b(4)
      integer :: i, j
      real(wp) :: ka, kb
      do i = 2, 4
         ka = a(i)
         kb = b(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= ka) exit
            a(j + 1) = a(j)
            b(j + 1) = b(j)
            j = j - 1
         end do
         a(j + 1) = ka
         b(j + 1) = kb
      end do
   end subroutine sort4

   ! =====================================================================
   ! Apply: one device kernel + host conveniences.
   ! =====================================================================

   pure subroutine regrid_apply_2d_impl(nx, ny, nmax, nsx, nsy, idx_i, idx_j, w, src, dst)
      !! The ONE device kernel: physical-only target, explicit shapes.
      integer, intent(in) :: nx, ny, nmax, nsx, nsy
      integer, intent(in) :: idx_i(nx, ny, nmax), idx_j(nx, ny, nmax)
      real(wp), intent(in) :: w(nx, ny, nmax)
      real(wp), intent(in) :: src(nsx, nsy)
      real(wp), intent(inout) :: dst(nx, ny)
      integer :: i, j, k
      real(wp) :: s
      do concurrent(j=1:ny, i=1:nx)
         s = 0.0_wp
         do k = 1, nmax
            s = s + w(i, j, k)*src(idx_i(i, j, k), idx_j(i, j, k))
         end do
         dst(i, j) = s
      end do
   end subroutine regrid_apply_2d_impl

   subroutine regrid_apply_2d(wts, nsx, nsy, src, n1, n2, nghost, dst)
      !! Host shim: applies `wts` and writes the physical interior of a
      !! `(n1,n2)` ghosted destination. Ghost fill is the consumer's.
      type(regrid_weights_t), intent(in) :: wts
      integer, intent(in) :: nsx, nsy, n1, n2, nghost
      real(wp), intent(in) :: src(nsx, nsy)
      real(wp), intent(inout) :: dst(n1, n2)
      real(wp) :: dst_phys(wts%nx, wts%ny)
      integer :: i, j

      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, nsx, nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst_phys)
      do j = 1, wts%ny
         do i = 1, wts%nx
            dst(nghost + i, nghost + j) = dst_phys(i, j)
         end do
      end do
   end subroutine regrid_apply_2d

   pure subroutine regrid_window(wts, j_lo, nsy)
      !! Getter for the C2c reader hyperslab `(1, j_lo, nlon, nsy)`.
      type(regrid_weights_t), intent(in) :: wts
      integer, intent(out) :: j_lo, nsy
      j_lo = wts%j_lo
      nsy = wts%nsy
   end subroutine regrid_window

end module rdb_regrid_latlon
