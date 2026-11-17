!! Lon-lat runoff product -> model-grid runoff file (offline helper).
program om_runoff_to_model
   !! Offline preprocessing helper for `&ocean_dataovr_nml lrunoff / frunoff`.  Not part of the
   !! model, which reads runoff only from a file ALREADY on its grid with every value on a wet
   !! cell (a value on a land cell aborts, because that mass would vanish silently).
   !!
   !! Given a regular lon-lat runoff field in kg m-2 s-1 (JRA55-do `friver` or `licalvf`; Tsujino
   !! et al., 2018), it
   !!
   !! 1. converts every source cell to a mass rate (kg/s) with its exact spherical area
   !!    `R^2 dlon (sin(lat_n) - sin(lat_s))`;
   !! 2. remaps it conservatively onto the model's T cells: each source cell is cut into
   !!    `msub x msub` equal-AREA sub-cells whose centres are located in the model grid (the
   !!    nearest cell whose quadrilateral, projected gnomonically on its own tangent plane,
   !!    contains the point, else the nearest centre when the point lies within two cell radii of it),
   !!    and every sub-cell hands its equal share of the source mass to the cell it falls in.  A
   !!    sub-cell farther than that from every model cell is OUTSIDE the model domain (a regional
   !!    grid, the polar hole of a lon-lat mosaic): its share is DISCARDED, never dumped on an edge
   !!    cell, and its total (kg/s) and sub-cell count are printed per record.  The total check is
   !!    `in-domain after = source - discarded`; `--global` makes any discard fatal (use it for a
   !!    global grid, where nothing may be out of domain).  A point poleward of the mosaic is
   !!    outside; a point within two cell radii of the fold, the poles of a tripolar cap or the
   !!    dateline seam of a periodic grid is never discarded (the periodic i-wrap is searched);
   !! 3. moves every share that landed on a model LAND cell to the nearest wet COASTAL cell (a
   !!    wet cell with a land 4-neighbour; any wet cell if the model has none) by great-circle
   !!    distance between cell centres, the role of the coupled models' runoff-to-coast map;
   !! 3b. spreads every cell's mass over the wet cells reachable from it through wet cells only (4-neighbour
   !!    flood fill, periodic i-wrap, north fold) within `--spread-radius-km` (default 500) with weights
   !!    exp(-d/`--spread-efold-km`) (default 150), normalised per source, so a river is not released into one
   !!    thin coastal column and cannot cross a land bridge; radius 0 turns it off (bitwise the unspread
   !!    answer).  `--spread-max-rate` (m/day) widens a source that stays above it up to `--spread-hard-km`.
   !!    The operator is static and built once per source cell; the global total is asserted unchanged.
   !! 4. divides by the target cell area and writes `runoff(x, y, time)` (kg m-2 s-1, float64) with
   !!    the source's time axis and `units` attribute, the layout the reader requires.
   !!
   !! It asserts that the total after (sum of out * target area) equals the total before (sum of
   !! source * source area) minus the discarded out-of-domain mass to 1e-12 relative, that no land
   !! cell is non-zero, and prints both totals and the discard per record.  A missing value (fill, NaN, |x| > 1e10) in the source counts as zero
   !! and is reported.  No third-party package; the self-test builds a synthetic mosaic, depth map
   !! and runoff product in the working directory and checks the same assertions.
   !!
   !! Arguments (positional): hgrid.nc depth.nc src.nc src_var out.nc
   !! Options: --depth-var NAME (default depth; wet where depth > 0)  --out-var NAME (default runoff)
   !!          --lon-name NAME (default lon) --lat-name NAME (default lat) --time-name NAME (default time)
   !!          --msub N (default 4)  --min-depth M (wet where depth > M, default 0)  --global
   !!          --selftest   (no other arguments)
   use, intrinsic :: iso_fortran_env, only: real64, error_unit, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use netcdf, only: nf90_open, nf90_close, nf90_create, nf90_enddef, nf90_noerr, &
                     nf90_nowrite, nf90_clobber, nf90_64bit_offset, nf90_double, &
                     nf90_global, nf90_inq_varid, nf90_inq_dimid, &
                     nf90_inquire_dimension, nf90_inquire_variable, nf90_get_var, &
                     nf90_put_var, nf90_def_dim, nf90_def_var, nf90_put_att, nf90_get_att, &
                     nf90_strerror, nf90_inquire_attribute
   implicit none

   integer, parameter :: dp = real64
   real(dp), parameter :: PI = acos(-1.0_dp)
   real(dp), parameter :: DEG = PI/180.0_dp
   real(dp), parameter :: R_EARTH = 6.378e6_dp
      !! Earth radius (m) of the source-cell areas.
   real(dp), parameter :: FILL_LIMIT = 1.0e10_dp
      !! Any |value| above this (or non-finite) is a missing value, taken as zero.
   real(dp), parameter :: TOL_TOTAL = 1.0e-12_dp
   real(dp), parameter :: RHO_FW = 1000.0_dp
      !! Fresh-water density (kg/m3) of the m/day rates the spreading reports and limits.

   type :: point_index_t
      !! Nearest-neighbour index over unit vectors: a cubic bin grid (no pole or dateline
      !! singularity), CSR bin lists.
      integer :: nb = 0
      integer :: n = 0
      integer, allocatable :: bstart(:), ids(:)
      real(dp), allocatable :: v(:, :)
   end type point_index_t

   integer, parameter :: NXS = 72, NYS = 32, NSX = 180, NSY = 90
      !! The five-degree mosaic and two-degree source of the spreading self-test.
   integer, parameter :: SPREAD_MAXLEV = 8
      !! Radius levels of the spreading operator (level 1 = `--spread-radius-km`, each next one x1.5).

   type :: spread_t
      !! The static spreading operator, built lazily one source cell and radius level at a time and
      !! kept: for source cell k at level l, `tgt/wgt(start(k,l) : start(k,l)+cnt(k,l)-1)` are the
      !! wet cells it feeds and the shares (they sum to 1).
      integer :: nl = 1
      real(dp) :: rad(SPREAD_MAXLEV) = 0.0_dp
      real(dp) :: efold = 100.0_dp
      integer :: used = 0
      integer, allocatable :: start(:, :), cnt(:, :), tgt(:), stamp(:), queue(:)
      real(dp), allocatable :: wgt(:)
      integer :: stampn = 0
   end type spread_t

   character(len=1024) :: arg
   character(len=1024) :: p_hgrid, p_depth, p_src, p_out
   character(len=64) :: src_var, depth_var, out_var, lon_name, lat_name, time_name
   integer :: msub, iarg, nargs
   real(dp) :: min_depth
   logical :: strict_global
   real(dp) :: spread_r, spread_e, spread_max, spread_hard

   nargs = command_argument_count()
   depth_var = "depth"
   out_var = "runoff"
   lon_name = "lon"
   lat_name = "lat"
   time_name = "time"
   msub = 4
   min_depth = 0.0_dp
   strict_global = .false.
   spread_r = 500.0_dp
   spread_e = 150.0_dp
   spread_max = 0.0_dp
   spread_hard = 1500.0_dp

   if (nargs == 1) then
      call get_command_argument(1, arg)
      if (trim(arg) == "--selftest") then
         call selftest()
         stop 0
      end if
   end if
   if (nargs < 5) then
      write (error_unit, "(a)") "usage: om_runoff_to_model hgrid.nc depth.nc src.nc src_var out.nc "// &
         "[--depth-var NAME] [--out-var NAME] [--lon-name NAME] [--lat-name NAME] [--time-name NAME] "// &
         "[--msub N] [--min-depth M] [--global] [--spread-radius-km R] [--spread-efold-km L] "// &
         "[--spread-max-rate M_PER_DAY] [--spread-hard-km R] | --selftest"
      stop 1
   end if
   call get_command_argument(1, p_hgrid)
   call get_command_argument(2, p_depth)
   call get_command_argument(3, p_src)
   call get_command_argument(4, src_var)
   call get_command_argument(5, p_out)
   iarg = 6
   do while (iarg <= nargs)
      call get_command_argument(iarg, arg)
      select case (trim(arg))
      case ("--depth-var")
         call get_command_argument(iarg + 1, depth_var)
         iarg = iarg + 2
      case ("--out-var")
         call get_command_argument(iarg + 1, out_var)
         iarg = iarg + 2
      case ("--lon-name")
         call get_command_argument(iarg + 1, lon_name)
         iarg = iarg + 2
      case ("--lat-name")
         call get_command_argument(iarg + 1, lat_name)
         iarg = iarg + 2
      case ("--time-name")
         call get_command_argument(iarg + 1, time_name)
         iarg = iarg + 2
      case ("--msub")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) msub
         iarg = iarg + 2
      case ("--min-depth")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) min_depth
         iarg = iarg + 2
      case ("--spread-radius-km")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) spread_r
         iarg = iarg + 2
      case ("--spread-efold-km")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) spread_e
         iarg = iarg + 2
      case ("--spread-max-rate")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) spread_max
         iarg = iarg + 2
      case ("--spread-hard-km")
         call get_command_argument(iarg + 1, arg)
         read (arg, *) spread_hard
         iarg = iarg + 2
      case ("--global")
         strict_global = .true.
         iarg = iarg + 1
      case default
         write (error_unit, "(a)") "unknown option: "//trim(arg)
         stop 1
      end select
   end do
   if (msub < 1) then
      write (error_unit, "(a)") "--msub must be >= 1"
      stop 1
   end if
   call convert(trim(p_hgrid), trim(p_depth), trim(depth_var), trim(p_src), trim(src_var), &
                trim(p_out), trim(out_var), trim(lon_name), trim(lat_name), trim(time_name), msub, &
                min_depth, strict_global, spread_radius_km=spread_r, spread_efold_km=spread_e, &
                spread_max_rate=spread_max, spread_hard_km=spread_hard)

contains

   subroutine check(status, what)
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a)") "om_runoff_to_model: "//what//": "//trim(nf90_strerror(status))
         stop 1
      end if
   end subroutine check

   function itoa(n) result(str)
      integer, intent(in) :: n
      character(len=:), allocatable :: str
      character(len=16) :: buf
      write (buf, "(i0)") n
      str = trim(buf)
   end function itoa

   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      write (error_unit, "(a)") "om_runoff_to_model: "//msg
      stop 1
   end subroutine fail

   ! ====================================================================
   ! geometry
   ! ====================================================================

   pure function unit_vec(lon, lat) result(v)
      !! Unit vector of a (lon, lat) in degrees.
      real(dp), intent(in) :: lon, lat
      real(dp) :: v(3)
      v = [cos(lat*DEG)*cos(lon*DEG), cos(lat*DEG)*sin(lon*DEG), sin(lat*DEG)]
   end function unit_vec

   pure function cross3(a, b) result(c)
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross3

   pure function bin_of(x, nb) result(b)
      real(dp), intent(in) :: x
      integer, intent(in) :: nb
      integer :: b
      b = min(nb, max(1, int(floor((x + 1.0_dp)*0.5_dp*real(nb, dp))) + 1))
   end function bin_of

   subroutine build_index(idx, n, v)
      !! Index `n` unit vectors `v(3, n)`.
      type(point_index_t), intent(out) :: idx
      integer, intent(in) :: n
      real(dp), intent(in) :: v(3, n)
      integer :: nb, k, b, ix, iy, iz
      integer, allocatable :: cnt(:), pos(:)

      idx%n = n
      nb = max(2, min(64, nint((real(max(n, 1), dp)/4.0_dp)**(1.0_dp/3.0_dp)) + 1))
      idx%nb = nb
      allocate (idx%v(3, n), idx%bstart(nb*nb*nb + 1), idx%ids(n), cnt(nb*nb*nb), pos(nb*nb*nb))
      idx%v = v
      cnt = 0
      do k = 1, n
         ix = bin_of(v(1, k), nb)
         iy = bin_of(v(2, k), nb)
         iz = bin_of(v(3, k), nb)
         b = ix + (iy - 1)*nb + (iz - 1)*nb*nb
         cnt(b) = cnt(b) + 1
      end do
      idx%bstart(1) = 1
      do b = 1, nb*nb*nb
         idx%bstart(b + 1) = idx%bstart(b) + cnt(b)
      end do
      pos = idx%bstart(1:nb*nb*nb)
      do k = 1, n
         ix = bin_of(v(1, k), nb)
         iy = bin_of(v(2, k), nb)
         iz = bin_of(v(3, k), nb)
         b = ix + (iy - 1)*nb + (iz - 1)*nb*nb
         idx%ids(pos(b)) = k
         pos(b) = pos(b) + 1
      end do
   end subroutine build_index

   function nearest_point(idx, q) result(ibest)
      !! Index (into the points handed to `build_index`) of the point nearest to unit vector `q`.
      type(point_index_t), intent(in) :: idx
      real(dp), intent(in) :: q(3)
      integer :: ibest
      integer :: nb, bx, by, bz, r, dx, dy, dz, ix, iy, iz, b, m, k
      real(dp) :: best, d2, binw

      nb = idx%nb
      binw = 2.0_dp/real(nb, dp)
      bx = bin_of(q(1), nb)
      by = bin_of(q(2), nb)
      bz = bin_of(q(3), nb)
      ibest = 0
      best = huge(1.0_dp)
      do r = 0, nb
         if (ibest > 0) then
            if (real(r - 1, dp)*binw > sqrt(best)) exit
         end if
         do dz = -r, r
            iz = bz + dz
            if (iz < 1 .or. iz > nb) cycle
            do dy = -r, r
               iy = by + dy
               if (iy < 1 .or. iy > nb) cycle
               do dx = -r, r
                  if (max(abs(dx), abs(dy), abs(dz)) /= r) cycle
                  ix = bx + dx
                  if (ix < 1 .or. ix > nb) cycle
                  b = ix + (iy - 1)*nb + (iz - 1)*nb*nb
                  do m = idx%bstart(b), idx%bstart(b + 1) - 1
                     k = idx%ids(m)
                     d2 = sum((idx%v(:, k) - q)**2)
                     if (d2 < best) then
                        best = d2
                        ibest = k
                     end if
                  end do
               end do
            end do
         end do
      end do
   end function nearest_point

   function quad_contains(c4, q) result(inside)
      !! `.true.` when unit vector `q` lies inside the quadrilateral of the four unit-vector corners
      !! `c4(:, 1:4)` (in order around the cell), judged in the gnomonic projection onto the tangent
      !! plane at the corners' mean direction.
      real(dp), intent(in) :: c4(3, 4), q(3)
      logical :: inside
      real(dp) :: c(3), e1(3), e2(3), p(2), a(2, 4), ref(3)
      real(dp) :: s, cr
      integer :: k, k2, npos, nneg

      c = sum(c4, dim=2)
      c = c/sqrt(sum(c**2))
      inside = .false.
      s = dot_product(q, c)
      if (s <= 0.0_dp) return
      if (abs(c(3)) < 0.9_dp) then
         ref = [0.0_dp, 0.0_dp, 1.0_dp]
      else
         ref = [1.0_dp, 0.0_dp, 0.0_dp]
      end if
      e1 = cross3(ref, c)
      e1 = e1/sqrt(sum(e1**2))
      e2 = cross3(c, e1)
      do k = 1, 4
         s = dot_product(c4(:, k), c)
         a(1, k) = dot_product(c4(:, k), e1)/s
         a(2, k) = dot_product(c4(:, k), e2)/s
      end do
      s = dot_product(q, c)
      p(1) = dot_product(q, e1)/s
      p(2) = dot_product(q, e2)/s
      npos = 0
      nneg = 0
      do k = 1, 4
         k2 = mod(k, 4) + 1
         cr = (a(1, k2) - a(1, k))*(p(2) - a(2, k)) - (a(2, k2) - a(2, k))*(p(1) - a(1, k))
         if (cr > 0.0_dp) npos = npos + 1
         if (cr < 0.0_dp) nneg = nneg + 1
      end do
      inside = (npos == 0 .or. nneg == 0)
   end function quad_contains

   ! ====================================================================
   ! the model grid and the point location
   ! ====================================================================

   subroutine cell_corners(sg_x, sg_y, nx, ny, i, j, c4)
      !! The four corners (SW, SE, NE, NW) of T cell `(i, j)` as unit vectors, from the supergrid
      !! nodes `(2i - 1 .. 2i + 1, 2j - 1 .. 2j + 1)`.
      integer, intent(in) :: nx, ny, i, j
      real(dp), intent(in) :: sg_x(2*nx + 1, 2*ny + 1), sg_y(2*nx + 1, 2*ny + 1)
      real(dp), intent(out) :: c4(3, 4)
      c4(:, 1) = unit_vec(sg_x(2*i - 1, 2*j - 1), sg_y(2*i - 1, 2*j - 1))
      c4(:, 2) = unit_vec(sg_x(2*i + 1, 2*j - 1), sg_y(2*i + 1, 2*j - 1))
      c4(:, 3) = unit_vec(sg_x(2*i + 1, 2*j + 1), sg_y(2*i + 1, 2*j + 1))
      c4(:, 4) = unit_vec(sg_x(2*i - 1, 2*j + 1), sg_y(2*i - 1, 2*j + 1))
   end subroutine cell_corners

   function locate(idx, sg_x, sg_y, nx, ny, q, periodic) result(cell)
      !! Flat index (`i + (j - 1) nx`) of the model T cell holding unit vector `q`: the nearest
      !! centre's quadrilateral, then its eight grid neighbours' (wrapping in i when `periodic`),
      !! else that nearest centre if `q` is within two cell radii of it; `0` (OUTSIDE the domain)
      !! beyond that.  The radius test, not the quadrilateral test, decides "outside", so a point
      !! the planar quadrilaterals miss near a fold, a pole or a seam is still placed.
      type(point_index_t), intent(in) :: idx
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: sg_x(2*nx + 1, 2*ny + 1), sg_y(2*nx + 1, 2*ny + 1)
      real(dp), intent(in) :: q(3)
      logical, intent(in) :: periodic
      integer :: cell
      integer :: k0, i0, j0, di, dj, i, j, k
      real(dp) :: rfar2, d2
      real(dp) :: c4(3, 4), ctr(3)

      k0 = nearest_point(idx, q)
      cell = k0
      i0 = mod(k0 - 1, nx) + 1
      j0 = (k0 - 1)/nx + 1
      call cell_corners(sg_x, sg_y, nx, ny, i0, j0, c4)
      if (quad_contains(c4, q)) return
      do dj = -1, 1
         do di = -1, 1
            if (di == 0 .and. dj == 0) cycle
            i = i0 + di
            j = j0 + dj
            if (periodic) i = modulo(i - 1, nx) + 1
            if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
            call cell_corners(sg_x, sg_y, nx, ny, i, j, c4)
            if (quad_contains(c4, q)) then
               cell = i + (j - 1)*nx
               return
            end if
         end do
      end do
      call cell_corners(sg_x, sg_y, nx, ny, i0, j0, c4)
      ctr = unit_vec(sg_x(2*i0, 2*j0), sg_y(2*i0, 2*j0))
      rfar2 = 0.0_dp
      do k = 1, 4
         rfar2 = max(rfar2, sum((c4(:, k) - ctr)**2))
      end do
      d2 = sum((q - ctr)**2)
      if (d2 > 4.0_dp*rfar2) cell = 0
   end function locate

   function grid_is_periodic(sg_x, sg_y, nxp, nyp) result(per)
      !! The mosaic closes on itself in i (first and last node columns coincide on the sphere).
      integer, intent(in) :: nxp, nyp
      real(dp), intent(in) :: sg_x(nxp, nyp), sg_y(nxp, nyp)
      logical :: per
      integer :: j
      per = .true.
      do j = 1, nyp
         if (sum((unit_vec(sg_x(1, j), sg_y(1, j)) - unit_vec(sg_x(nxp, j), sg_y(nxp, j)))**2) > 1.0e-18_dp) per = .false.
      end do
   end function grid_is_periodic

   function grid_has_fold(sg_x, sg_y, nxp, nyp) result(fold)
      !! A tripolar north fold: the top node row is its own mirror (`i <-> nxp + 1 - i`) and is
      !! not a degenerate pole row (all nodes equal).
      integer, intent(in) :: nxp, nyp
      real(dp), intent(in) :: sg_x(nxp, nyp), sg_y(nxp, nyp)
      logical :: fold
      integer :: i
      logical :: distinct
      fold = .true.
      distinct = .false.
      do i = 1, nxp
         if (sum((unit_vec(sg_x(i, nyp), sg_y(i, nyp)) - unit_vec(sg_x(nxp + 1 - i, nyp), sg_y(nxp + 1 - i, nyp)))**2) &
             > 1.0e-18_dp) fold = .false.
         if (i < nxp) then
            if (sum((unit_vec(sg_x(i, nyp), sg_y(i, nyp)) - unit_vec(sg_x(i + 1, nyp), sg_y(i + 1, nyp)))**2) &
                > 1.0e-18_dp) distinct = .true.
         end if
      end do
      fold = fold .and. distinct
   end function grid_has_fold

   ! ====================================================================
   ! the conversion
   ! ====================================================================

   subroutine convert(p_hgrid, p_depth, depth_var, p_src, src_var, p_out, out_var, lon_name, lat_name, time_name, msub, &
                      min_depth, strict_global, discarded, n_redirected, max_redirect_km, &
                      spread_radius_km, spread_efold_km, spread_max_rate, spread_hard_km)
      character(len=*), intent(in) :: p_hgrid, p_depth, depth_var, p_src, src_var, p_out, out_var
      character(len=*), intent(in) :: lon_name, lat_name, time_name
      integer, intent(in) :: msub
      real(dp), intent(in) :: min_depth
         !! A cell is wet where `depth > min_depth`.
      logical, intent(in) :: strict_global
         !! Any out-of-domain source mass is fatal.
      real(dp), allocatable, intent(out), optional :: discarded(:)
         !! Per record, the discarded out-of-domain mass (kg/s).
      integer, intent(out), optional :: n_redirected
         !! Land cells that carried mass in some record and were moved to the coast.
      real(dp), intent(out), optional :: max_redirect_km
         !! The farthest such move (great-circle, centre to centre).
      real(dp), intent(in), optional :: spread_radius_km, spread_efold_km, spread_max_rate, spread_hard_km
         !! Spreading (see the program docstring); radius 0 or absent = none.

      integer :: nx, ny, nxp, nyp, nlon, nlat, nt
      real(dp), allocatable :: sg_x(:, :), sg_y(:, :), area_sg(:, :), area_t(:, :), depth(:, :)
      real(dp), allocatable :: slon(:), slat(:), stime(:), asrc(:, :), f(:, :), mass(:), out(:, :, :)
      real(dp), allocatable :: cv(:, :), wetv(:, :)
      integer, allocatable :: tidx(:, :, :), tidx0(:, :, :), redir(:), coast_id(:)
      logical, allocatable :: wet(:, :), hit(:), coastal(:, :), moved(:)
      real(dp) :: dkm, dkm_max, moved_mass, d_dot
      integer :: n_moved_used, k_far
      type(spread_t) :: sp
      real(dp), allocatable :: mass2(:), rate(:)
      integer, allocatable :: lev(:)
      real(dp) :: sp_r, sp_hard, sp_max, rate_b, rate_a, tot_after_spread
      integer :: iter, m, n_b, n_a, k_b, k_a, n_recv_b, n_recv_a, n_over
      logical :: do_spread, changed
      character(len=256) :: tunits, tcal
      logical :: has_cal
      type(point_index_t) :: cidx, ridx
      integer :: i, j, k, t, a, b, nwet, ncoast, nmoved_cells, nmiss, n_out_sub
      logical :: periodic, fold
      real(dp) :: before, after, lon0, lat0, dlon, dlat, disc
      real(dp) :: v(3)
      real(dp) :: lonw, sins, sinn, w_sub
      integer :: ncid, vid

      call read_supergrid(p_hgrid, nx, ny, sg_x, sg_y, area_sg)
      nxp = 2*nx + 1
      nyp = 2*ny + 1
      allocate (area_t(nx, ny))
      do j = 1, ny
         do i = 1, nx
            area_t(i, j) = area_sg(2*i - 1, 2*j - 1) + area_sg(2*i, 2*j - 1) &
                           + area_sg(2*i - 1, 2*j) + area_sg(2*i, 2*j)
         end do
      end do
      allocate (depth(nx, ny), wet(nx, ny))
      call check(nf90_open(p_depth, nf90_nowrite, ncid), "open "//p_depth)
      call check(nf90_inq_varid(ncid, depth_var, vid), "variable "//depth_var)
      call check(nf90_get_var(ncid, vid, depth), "read "//depth_var)
      call check(nf90_close(ncid), "close depth")
      wet = depth > min_depth
      nwet = count(wet)
      if (nwet == 0) call fail("the depth map has no wet cell")

      ! source axes and areas
      call read_source_axes(p_src, src_var, lon_name, lat_name, time_name, nlon, nlat, nt, slon, slat, stime, tunits, tcal, has_cal)
      allocate (asrc(nlon, nlat))
      dlon = slon(2) - slon(1)
      if (nlon < 2 .or. nlat < 2) call fail("the source grid must have at least 2 points per axis")
      if (maxval(abs(slon(2:) - slon(:nlon - 1) - dlon)) > 1.0e-6_dp*abs(dlon)) call fail("source longitudes must be regular")
      dlat = slat(2) - slat(1)
      if (dlat <= 0.0_dp .or. maxval(abs(slat(2:) - slat(:nlat - 1) - dlat)) > 1.0e-6_dp*dlat) then
         call fail("source latitudes must be regular and increasing")
      end if
      do j = 1, nlat
         sins = sin(max(-90.0_dp, slat(j) - 0.5_dp*dlat)*DEG)
         sinn = sin(min(90.0_dp, slat(j) + 0.5_dp*dlat)*DEG)
         do i = 1, nlon
            asrc(i, j) = R_EARTH**2*abs(dlon)*DEG*(sinn - sins)
         end do
      end do

      ! target centres -> index
      allocate (cv(3, nx*ny))
      do j = 1, ny
         do i = 1, nx
            cv(:, i + (j - 1)*nx) = unit_vec(sg_x(2*i, 2*j), sg_y(2*i, 2*j))
         end do
      end do
      call build_index(cidx, nx*ny, cv)
      periodic = grid_is_periodic(sg_x, sg_y, nxp, nyp)
      fold = grid_has_fold(sg_x, sg_y, nxp, nyp)

      ! the point map: msub x msub equal-area sub-cells of every source cell
      allocate (tidx(msub*msub, nlon, nlat), hit(nx*ny))
      hit = .false.
      do j = 1, nlat
         sins = sin(max(-90.0_dp, slat(j) - 0.5_dp*dlat)*DEG)
         sinn = sin(min(90.0_dp, slat(j) + 0.5_dp*dlat)*DEG)
         do i = 1, nlon
            lonw = slon(i) - 0.5_dp*dlon
            do b = 1, msub
               lat0 = asin(sins + (real(b, dp) - 0.5_dp)/real(msub, dp)*(sinn - sins))/DEG
               do a = 1, msub
                  lon0 = lonw + (real(a, dp) - 0.5_dp)/real(msub, dp)*dlon
                  v = unit_vec(lon0, lat0)
                  k = locate(cidx, sg_x, sg_y, nx, ny, v, periodic)
                  tidx(a + (b - 1)*msub, i, j) = k
                  if (k > 0) hit(k) = .true.
               end do
            end do
         end do
      end do

      ! land -> nearest wet coastal cell
      allocate (coastal(nx, ny))
      coastal = .false.
      do j = 1, ny
         do i = 1, nx
            if (.not. wet(i, j)) cycle
            ! the 4 neighbours: wrapping in i on a periodic mosaic, mirrored across a north fold
            if (i > 1) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(i - 1, j)
            else if (periodic) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(nx, j)
            end if
            if (i < nx) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(i + 1, j)
            else if (periodic) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(1, j)
            end if
            if (j > 1) coastal(i, j) = coastal(i, j) .or. .not. wet(i, j - 1)
            if (j < ny) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(i, j + 1)
            else if (fold) then
               coastal(i, j) = coastal(i, j) .or. .not. wet(nx + 1 - i, j)
            end if
         end do
      end do
      ncoast = count(coastal)
      if (ncoast == 0) then
         coastal = wet
         ncoast = nwet
      end if
      allocate (coast_id(ncoast), wetv(3, ncoast))
      k = 0
      do j = 1, ny
         do i = 1, nx
            if (coastal(i, j)) then
               k = k + 1
               coast_id(k) = i + (j - 1)*nx
               wetv(:, k) = cv(:, i + (j - 1)*nx)
            end if
         end do
      end do
      call build_index(ridx, ncoast, wetv)
      allocate (redir(nx*ny))
      nmoved_cells = 0
      do k = 1, nx*ny
         redir(k) = k
         i = mod(k - 1, nx) + 1
         j = (k - 1)/nx + 1
         if (hit(k) .and. .not. wet(i, j)) then
            redir(k) = coast_id(nearest_point(ridx, cv(:, k)))
            nmoved_cells = nmoved_cells + 1
         end if
      end do
      allocate (tidx0(msub*msub, nlon, nlat), moved(nx*ny))
      tidx0 = tidx
      moved = .false.
      do j = 1, nlat
         do i = 1, nlon
            do a = 1, msub*msub
               if (tidx(a, i, j) > 0) tidx(a, i, j) = redir(tidx(a, i, j))
            end do
         end do
      end do

      ! the spreading operator (static, built lazily per source cell)
      sp_r = 0.0_dp
      sp_hard = 1500.0_dp
      sp_max = 0.0_dp
      sp%efold = 100.0_dp
      if (present(spread_radius_km)) sp_r = spread_radius_km
      if (present(spread_efold_km)) sp%efold = spread_efold_km
      if (present(spread_max_rate)) sp_max = spread_max_rate
      if (present(spread_hard_km)) sp_hard = spread_hard_km
      do_spread = sp_r > 0.0_dp
      if (sp_r < 0.0_dp .or. sp%efold <= 0.0_dp .or. sp_max < 0.0_dp) call fail("bad spreading option")
      if (do_spread) then
         sp%nl = 1
         sp%rad(1) = sp_r
         if (sp_max > 0.0_dp) then
            do while (sp%nl < SPREAD_MAXLEV .and. sp%rad(sp%nl)*1.5_dp <= max(sp_hard, sp_r))
               sp%rad(sp%nl + 1) = sp%rad(sp%nl)*1.5_dp
               sp%nl = sp%nl + 1
            end do
         end if
         allocate (sp%start(nx*ny, sp%nl), sp%cnt(nx*ny, sp%nl), sp%stamp(nx*ny), sp%queue(nx*ny))
         sp%start = 0
         sp%cnt = 0
         sp%stamp = 0
         allocate (sp%tgt(1024), sp%wgt(1024), mass2(nx*ny), rate(nx*ny), lev(nx*ny))
         write (output_unit, "(a,f7.1,a,f7.1,a)") "spreading: radius ", sp_r, " km, e-folding ", sp%efold, &
            " km, through wet cells only"
      end if

      ! the records
      allocate (f(nlon, nlat), mass(nx*ny), out(nx, ny, max(nt, 1)))
      if (present(discarded)) allocate (discarded(max(nt, 1)))
      n_out_sub = count(tidx == 0)
      w_sub = 1.0_dp/real(msub*msub, dp)
      write (output_unit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "target ", nx, " x ", ny, " (", nwet, " wet, ", ncoast, &
         " coastal); source ", nlon
      write (output_unit, "(a,i0,a,i0,a,i0)") "  x ", nlat, ", records ", nt, "; land cells hit by any source cell moved to the coast: ", &
         nmoved_cells
      do t = 1, max(nt, 1)
         call read_source_record(p_src, src_var, nlon, nlat, nt, t, f, nmiss)
         mass = 0.0_dp
         before = 0.0_dp
         disc = 0.0_dp
         moved_mass = 0.0_dp
         do j = 1, nlat
            do i = 1, nlon
               if (f(i, j) == 0.0_dp) cycle
               before = before + f(i, j)*asrc(i, j)
               do a = 1, msub*msub
                  k = tidx(a, i, j)
                  if (tidx0(a, i, j) /= k) then
                     moved(tidx0(a, i, j)) = .true.
                     moved_mass = moved_mass + f(i, j)*asrc(i, j)*w_sub
                  end if
                  if (k == 0) then
                     disc = disc + f(i, j)*asrc(i, j)*w_sub
                  else
                     mass(k) = mass(k) + f(i, j)*asrc(i, j)*w_sub
                  end if
               end do
            end do
         end do
         if (do_spread) then
            rate_b = 0.0_dp
            k_b = 0
            n_recv_b = count(mass /= 0.0_dp)
            do k = 1, nx*ny
               if (mass(k) /= 0.0_dp) then
                  rate_a = mass(k)/area_t(mod(k - 1, nx) + 1, (k - 1)/nx + 1)/RHO_FW*86400.0_dp
                  if (rate_a > rate_b) then
                     rate_b = rate_a
                     k_b = k
                  end if
               end if
            end do
            n_b = count_over(mass, area_t, nx, ny, 1.0_dp)
            lev = 1
            do iter = 1, 30
               call spread_apply(sp, mass, mass2, lev, wet, cv, nx, ny, periodic, fold)
               if (sp_max <= 0.0_dp) exit
               rate = 0.0_dp
               do k = 1, nx*ny
                  rate(k) = mass2(k)/area_t(mod(k - 1, nx) + 1, (k - 1)/nx + 1)/RHO_FW*86400.0_dp
               end do
               changed = .false.
               do k = 1, nx*ny
                  if (mass(k) == 0.0_dp .or. lev(k) >= sp%nl) cycle
                  do m = 0, sp%cnt(k, lev(k)) - 1
                     if (rate(sp%tgt(sp%start(k, lev(k)) + m)) > sp_max) then
                        lev(k) = lev(k) + 1
                        changed = .true.
                        exit
                     end if
                  end do
               end do
               if (.not. changed) exit
            end do
            tot_after_spread = sum(mass2)
            if (abs(tot_after_spread - sum(mass)) > TOL_TOTAL*max(abs(sum(mass)), tiny(1.0_dp))) then
               call fail("spreading did not conserve the global runoff")
            end if
            rate_a = 0.0_dp
            k_a = 0
            do k = 1, nx*ny
               if (mass2(k) /= 0.0_dp) then
                  rate = 0.0_dp
                  rate(1) = mass2(k)/area_t(mod(k - 1, nx) + 1, (k - 1)/nx + 1)/RHO_FW*86400.0_dp
                  if (rate(1) > rate_a) then
                     rate_a = rate(1)
                     k_a = k
                  end if
               end if
            end do
            n_recv_a = count(mass2 /= 0.0_dp)
            n_a = count_over(mass2, area_t, nx, ny, 1.0_dp)
            n_over = 0
            if (sp_max > 0.0_dp) n_over = count_over(mass2, area_t, nx, ny, sp_max)
            write (output_unit, "(a,i0,a,i0,a,f9.3,a,i0,a,i0,a,f9.3,a,i0,a,i0,a,i0)") &
               "    spread: cells with runoff ", n_recv_b, " -> ", n_recv_a, "; max fresh-water rate ", rate_b, &
               " m/day (cell ", mod(k_b - 1, nx) + 1, ",", (k_b - 1)/nx + 1, ") -> ", rate_a, " m/day (cell ", &
               mod(k_a - 1, nx) + 1, ",", (k_a - 1)/nx + 1
            write (output_unit, "(a,i0,a,i0,a,i0)") "    spread: cells above 1 m/day ", n_b, " -> ", n_a, &
               "; sources widened beyond level 1: ", count(lev > 1 .and. mass /= 0.0_dp)
            if (sp_max > 0.0_dp) write (output_unit, "(a,es10.3,a,i0)") "    spread: cells still above ", sp_max, " m/day: ", n_over
            mass = mass2
         end if
         after = 0.0_dp
         do j = 1, ny
            do i = 1, nx
               k = i + (j - 1)*nx
               if (.not. wet(i, j) .and. mass(k) /= 0.0_dp) call fail("internal: mass left on a land cell")
               out(i, j, t) = mass(k)/area_t(i, j)
               after = after + out(i, j, t)*area_t(i, j)
            end do
         end do
         write (output_unit, "(a,i0,a,es22.14,a,es22.14,a,i0,a)") "  record ", t, ": total before ", before, " kg/s, after ", &
            after, " kg/s  (", nmiss, " missing values taken as zero)"
         write (output_unit, "(a,es22.14,a,i0,a)") "    discarded (outside the model domain): ", disc, " kg/s in ", &
            n_out_sub, " source sub-cells"
         write (output_unit, "(a,es22.14,a)") "    redirected from land to the coast: ", moved_mass, " kg/s"
         if (present(discarded)) discarded(t) = disc
         if (strict_global .and. disc /= 0.0_dp) then
            call fail("--global: record "//itoa(t)//" has source mass outside the model domain")
         end if
         if (abs(after + disc - before) > TOL_TOTAL*max(abs(before), tiny(1.0_dp))) then
            call fail("the total was not preserved: in-domain + discarded /= source")
         end if
      end do
      ! how far did the mass travel?  Only land cells that actually carried mass in some record.
      n_moved_used = count(moved)
      dkm_max = 0.0_dp
      k_far = 0
      do k = 1, nx*ny
         if (.not. moved(k)) cycle
         d_dot = min(1.0_dp, max(-1.0_dp, dot_product(cv(:, k), cv(:, redir(k)))))
         dkm = 6371.0_dp*acos(d_dot)
         if (dkm > dkm_max) then
            dkm_max = dkm
            k_far = k
         end if
      end do
      if (present(n_redirected)) n_redirected = n_moved_used
      if (present(max_redirect_km)) max_redirect_km = dkm_max
      write (output_unit, "(a,i0,a,f8.1,a)") "land cells that carried mass and were redirected: ", n_moved_used, &
         "; farthest redirect ", dkm_max, " km"
      if (k_far > 0) then
        write (output_unit, "(a,i0,a,i0,a,i0,a,i0,a)") "  farthest: land cell (", mod(k_far - 1, nx) + 1, ",", (k_far - 1)/nx + 1, &
            ") -> coast cell (", mod(redir(k_far) - 1, nx) + 1, ",", (redir(k_far) - 1)/nx + 1, ")"
      end if
      write (output_unit, "(a)") "ok: in-domain total + discarded = source total, no land cell non-zero"
      call write_output(p_out, out_var, nx, ny, max(nt, 1), stime, tunits, tcal, has_cal, out, p_src)
   end subroutine convert

   ! ====================================================================
   ! I/O
   ! ====================================================================

   pure integer function count_over(mass, area_t, nx, ny, lim) result(n)
      !! Cells whose fresh-water rate (m/day) exceeds `lim`.
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: lim
      real(dp), intent(in) :: mass(nx*ny), area_t(nx, ny)
      integer :: k
      n = 0
      do k = 1, nx*ny
         if (mass(k) /= 0.0_dp) then
            if (mass(k)/area_t(mod(k - 1, nx) + 1, (k - 1)/nx + 1)/RHO_FW*86400.0_dp > lim) n = n + 1
         end if
      end do
   end function count_over

   pure subroutine spread_nbrs(i, j, nx, ny, periodic, fold, n, ni, nj)
      !! The 4-neighbours of (i,j) a flood fill may step to: wrapping in i on a periodic mosaic and
      !! mirrored (i -> nx+1-i, same row) across a north fold; no wrap otherwise.
      integer, intent(in) :: i, j, nx, ny
      logical, intent(in) :: periodic, fold
      integer, intent(out) :: n
      integer, intent(out) :: ni(4), nj(4)
      n = 0
      if (i > 1) then
         n = n + 1
         ni(n) = i - 1
         nj(n) = j
      else if (periodic) then
         n = n + 1
         ni(n) = nx
         nj(n) = j
      end if
      if (i < nx) then
         n = n + 1
         ni(n) = i + 1
         nj(n) = j
      else if (periodic) then
         n = n + 1
         ni(n) = 1
         nj(n) = j
      end if
      if (j > 1) then
         n = n + 1
         ni(n) = i
         nj(n) = j - 1
      end if
      if (j < ny) then
         n = n + 1
         ni(n) = i
         nj(n) = j + 1
      else if (fold) then
         n = n + 1
         ni(n) = nx + 1 - i
         nj(n) = j
      end if
   end subroutine spread_nbrs

   subroutine spread_build(sp, k, l, wet, cv, nx, ny, periodic, fold)
      !! Flood fill from wet cell k through WET cells only, within the great-circle radius of level l,
      !! weights exp(-d/efold) normalised to sum to 1.  A land bridge therefore cannot be crossed.
      type(spread_t), intent(inout) :: sp
      integer, intent(in) :: k, l, nx, ny
      logical, intent(in) :: periodic, fold
      logical, intent(in) :: wet(nx, ny)
      real(dp), intent(in) :: cv(3, nx*ny)
      integer :: head, tail, c, ci, cj, n, q, m
      integer :: ni(4), nj(4)
      integer, allocatable :: itmp(:)
      real(dp), allocatable :: rtmp(:)
      real(dp) :: d, wsum, dd
      sp%stampn = sp%stampn + 1
      head = 1
      tail = 1
      sp%queue(1) = k
      sp%stamp(k) = sp%stampn
      do while (head <= tail)
         c = sp%queue(head)
         head = head + 1
         ci = mod(c - 1, nx) + 1
         cj = (c - 1)/nx + 1
         call spread_nbrs(ci, cj, nx, ny, periodic, fold, n, ni, nj)
         do m = 1, n
            q = ni(m) + (nj(m) - 1)*nx
            if (.not. wet(ni(m), nj(m)) .or. sp%stamp(q) == sp%stampn) cycle
            dd = min(1.0_dp, max(-1.0_dp, dot_product(cv(:, k), cv(:, q))))
            if (6371.0_dp*acos(dd) > sp%rad(l)) cycle
            sp%stamp(q) = sp%stampn
            tail = tail + 1
            sp%queue(tail) = q
         end do
      end do
      if (sp%used + tail > size(sp%tgt)) then
         allocate (itmp(2*(sp%used + tail)), rtmp(2*(sp%used + tail)))
         itmp(1:sp%used) = sp%tgt(1:sp%used)
         rtmp(1:sp%used) = sp%wgt(1:sp%used)
         call move_alloc(itmp, sp%tgt)
         call move_alloc(rtmp, sp%wgt)
      end if
      wsum = 0.0_dp
      do m = 1, tail
         dd = min(1.0_dp, max(-1.0_dp, dot_product(cv(:, k), cv(:, sp%queue(m)))))
         d = 6371.0_dp*acos(dd)
         sp%tgt(sp%used + m) = sp%queue(m)
         sp%wgt(sp%used + m) = exp(-d/sp%efold)
         wsum = wsum + sp%wgt(sp%used + m)
      end do
      sp%wgt(sp%used + 1:sp%used + tail) = sp%wgt(sp%used + 1:sp%used + tail)/wsum
      sp%start(k, l) = sp%used + 1
      sp%cnt(k, l) = tail
      sp%used = sp%used + tail
   end subroutine spread_build

   subroutine spread_apply(sp, mass, mass2, lev, wet, cv, nx, ny, periodic, fold)
      !! mass2 = the spreading operator applied to mass (kg/s per cell), source cell k at radius level lev(k).
      type(spread_t), intent(inout) :: sp
      integer, intent(in) :: nx, ny
      integer, intent(in) :: lev(nx*ny)
      real(dp), intent(in) :: mass(nx*ny), cv(3, nx*ny)
      real(dp), intent(out) :: mass2(nx*ny)
      logical, intent(in) :: periodic, fold
      logical, intent(in) :: wet(nx, ny)
      integer :: k, m, l
      mass2 = 0.0_dp
      do k = 1, nx*ny
         if (mass(k) == 0.0_dp) cycle
         l = lev(k)
         if (sp%start(k, l) == 0) call spread_build(sp, k, l, wet, cv, nx, ny, periodic, fold)
         do m = sp%start(k, l), sp%start(k, l) + sp%cnt(k, l) - 1
            mass2(sp%tgt(m)) = mass2(sp%tgt(m)) + mass(k)*sp%wgt(m)
         end do
      end do
   end subroutine spread_apply

   subroutine read_supergrid(path, nx, ny, sg_x, sg_y, area_sg)
      character(len=*), intent(in) :: path
      integer, intent(out) :: nx, ny
      real(dp), allocatable, intent(out) :: sg_x(:, :), sg_y(:, :), area_sg(:, :)
      integer :: nc, vid, did, nxp, nyp

      call check(nf90_open(path, nf90_nowrite, nc), "open "//path)
      call check(nf90_inq_dimid(nc, "nxp", did), "nxp")
      call check(nf90_inquire_dimension(nc, did, len=nxp), "nxp len")
      call check(nf90_inq_dimid(nc, "nyp", did), "nyp")
      call check(nf90_inquire_dimension(nc, did, len=nyp), "nyp len")
      nx = (nxp - 1)/2
      ny = (nyp - 1)/2
      allocate (sg_x(nxp, nyp), sg_y(nxp, nyp), area_sg(2*nx, 2*ny))
      call check(nf90_inq_varid(nc, "x", vid), "x")
      call check(nf90_get_var(nc, vid, sg_x), "read x")
      call check(nf90_inq_varid(nc, "y", vid), "y")
      call check(nf90_get_var(nc, vid, sg_y), "read y")
      call check(nf90_inq_varid(nc, "area", vid), "area (the supergrid cell areas, m2)")
      call check(nf90_get_var(nc, vid, area_sg), "read area")
      call check(nf90_close(nc), "close hgrid")
   end subroutine read_supergrid

   subroutine read_source_axes(path, var, lon_name, lat_name, time_name, nlon, nlat, nt, slon, slat, stime, tunits, tcal, has_cal)
      character(len=*), intent(in) :: path, var, lon_name, lat_name, time_name
      integer, intent(out) :: nlon, nlat, nt
      real(dp), allocatable, intent(out) :: slon(:), slat(:), stime(:)
      character(len=256), intent(out) :: tunits, tcal
      logical, intent(out) :: has_cal
      integer :: nc, vid, did, vt, st, ndims
      integer :: dimids(4)

      call check(nf90_open(path, nf90_nowrite, nc), "open "//path)
      call check(nf90_inq_varid(nc, var, vid), "variable "//var)
      call check(nf90_inquire_variable(nc, vid, ndims=ndims, dimids=dimids), "inquire "//var)
      if (ndims /= 2 .and. ndims /= 3) call fail("the source variable must be (lon, lat) or (lon, lat, time)")
      call check(nf90_inq_dimid(nc, lon_name, did), "lon dim")
      call check(nf90_inquire_dimension(nc, did, len=nlon), "lon len")
      call check(nf90_inq_dimid(nc, lat_name, did), "lat dim")
      call check(nf90_inquire_dimension(nc, did, len=nlat), "lat len")
      allocate (slon(nlon), slat(nlat))
      call check(nf90_inq_varid(nc, lon_name, vt), "lon var")
      call check(nf90_get_var(nc, vt, slon), "read lon")
      call check(nf90_inq_varid(nc, lat_name, vt), "lat var")
      call check(nf90_get_var(nc, vt, slat), "read lat")
      tunits = ""
      tcal = ""
      has_cal = .false.
      if (ndims == 3) then
         call check(nf90_inq_dimid(nc, time_name, did), "time dim")
         call check(nf90_inquire_dimension(nc, did, len=nt), "time len")
         allocate (stime(nt))
         call check(nf90_inq_varid(nc, time_name, vt), "time var")
         call check(nf90_get_var(nc, vt, stime), "read time")
         call check(nf90_get_att(nc, vt, "units", tunits), "time units (required: the reader needs them)")
         st = nf90_get_att(nc, vt, "calendar", tcal)
         has_cal = (st == nf90_noerr)
      else
         nt = 0
         allocate (stime(1))
         stime = 0.0_dp
      end if
      call check(nf90_close(nc), "close source")
   end subroutine read_source_axes

   subroutine read_source_record(path, var, nlon, nlat, nt, t, f, nmiss)
      !! Record `t` (the whole field when the source has no time axis) with missing values zeroed.
      character(len=*), intent(in) :: path, var
      integer, intent(in) :: nlon, nlat, nt, t
      real(dp), intent(out) :: f(nlon, nlat)
      integer, intent(out) :: nmiss
      integer :: nc, vid, i, j

      call check(nf90_open(path, nf90_nowrite, nc), "open "//path)
      call check(nf90_inq_varid(nc, var, vid), "variable "//var)
      if (nt > 0) then
         call check(nf90_get_var(nc, vid, f, start=[1, 1, t], count=[nlon, nlat, 1]), "read record")
      else
         call check(nf90_get_var(nc, vid, f), "read field")
      end if
      call check(nf90_close(nc), "close source")
      nmiss = 0
      do j = 1, nlat
         do i = 1, nlon
            if (.not. ieee_is_finite(f(i, j))) then
               f(i, j) = 0.0_dp
               nmiss = nmiss + 1
            else if (abs(f(i, j)) > FILL_LIMIT) then
               f(i, j) = 0.0_dp
               nmiss = nmiss + 1
            end if
         end do
      end do
   end subroutine read_source_record

   subroutine write_output(path, var, nx, ny, nt, stime, tunits, tcal, has_cal, out, src_path)
      character(len=*), intent(in) :: path, var, tunits, tcal, src_path
      integer, intent(in) :: nx, ny, nt
      real(dp), intent(in) :: stime(:), out(nx, ny, nt)
      logical, intent(in) :: has_cal
      integer :: nc, d_x, d_y, d_t, v_t, v_f, k0

      call check(nf90_create(path, ior(nf90_clobber, nf90_64bit_offset), nc), "create "//path)
      call check(nf90_def_dim(nc, "x", nx, d_x), "def x")
      call check(nf90_def_dim(nc, "y", ny, d_y), "def y")
      call check(nf90_def_dim(nc, "time", nt, d_t), "def time")
      call check(nf90_def_var(nc, "time", nf90_double, [d_t], v_t), "def time var")
      if (len_trim(tunits) > 0) call check(nf90_put_att(nc, v_t, "units", trim(tunits)), "att")
      if (has_cal) call check(nf90_put_att(nc, v_t, "calendar", trim(tcal)), "att")
      call check(nf90_def_var(nc, var, nf90_double, [d_x, d_y, d_t], v_f), "def "//var)
      call check(nf90_put_att(nc, v_f, "units", "kg m-2 s-1"), "att")
      call check(nf90_put_att(nc, v_f, "long_name", "runoff mass flux into the ocean on model wet cells"), "att")
      call check(nf90_put_att(nc, nf90_global, "source", trim(src_path)), "att")
      call check(nf90_put_att(nc, nf90_global, "method", &
                              "conservative sub-cell remap, land values moved to the nearest wet coastal cell"), "att")
      call check(nf90_put_att(nc, nf90_global, "generated_by", "roundabout tools/om_runoff_to_model.f90"), "att")
      call check(nf90_enddef(nc), "enddef")
      if (size(stime) >= nt) then
         call check(nf90_put_var(nc, v_t, stime(1:nt)), "put time")
      else
         call check(nf90_put_var(nc, v_t, [(real(k0, dp), k0=0, nt - 1)]), "put time")
      end if
      call check(nf90_put_var(nc, v_f, out), "put field")
      call check(nf90_close(nc), "close out")
   end subroutine write_output

   ! ====================================================================
   ! self-test
   ! ====================================================================

   subroutine write_mosaic(path, nx, ny, lon0, lat0, dl)
      !! A lon-lat supergrid of `nx x ny` T cells of `dl` degrees from `(lon0, lat0)`.
      character(len=*), intent(in) :: path
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: lon0, lat0, dl
      real(dp), allocatable :: sg_x(:, :), sg_y(:, :), area(:, :)
      real(dp) :: sn, ss, hd
      integer :: i, j, nc, d1, d2, d3, v1, v2, v3

      allocate (sg_x(2*nx + 1, 2*ny + 1), sg_y(2*nx + 1, 2*ny + 1), area(2*nx, 2*ny))
      hd = 0.5_dp*dl
      do j = 1, 2*ny + 1
         do i = 1, 2*nx + 1
            sg_x(i, j) = lon0 + real(i - 1, dp)*hd
            sg_y(i, j) = lat0 + real(j - 1, dp)*hd
         end do
      end do
      do j = 1, 2*ny
         ss = sin(sg_y(1, j)*DEG)
         sn = sin(sg_y(1, j + 1)*DEG)
         do i = 1, 2*nx
            area(i, j) = R_EARTH**2*hd*DEG*(sn - ss)
         end do
      end do
      call check(nf90_create(path, ior(nf90_clobber, nf90_64bit_offset), nc), "create hgrid")
      call check(nf90_def_dim(nc, "nxp", 2*nx + 1, d1), "d")
      call check(nf90_def_dim(nc, "nyp", 2*ny + 1, d2), "d")
      call check(nf90_def_dim(nc, "nx", 2*nx, d3), "d")
      call check(nf90_def_var(nc, "x", nf90_double, [d1, d2], v1), "v")
      call check(nf90_def_var(nc, "y", nf90_double, [d1, d2], v2), "v")
      call check(nf90_def_dim(nc, "ny", 2*ny, d1), "d")
      call check(nf90_def_var(nc, "area", nf90_double, [d3, d1], v3), "v")
      call check(nf90_enddef(nc), "e")
      call check(nf90_put_var(nc, v1, sg_x), "p")
      call check(nf90_put_var(nc, v2, sg_y), "p")
      call check(nf90_put_var(nc, v3, area), "p")
      call check(nf90_close(nc), "c")
   end subroutine write_mosaic

   subroutine write_depth(path, depth)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: depth(:, :)
      integer :: nc, d1, d2, v1
      call check(nf90_create(path, ior(nf90_clobber, nf90_64bit_offset), nc), "create depth")
      call check(nf90_def_dim(nc, "x", size(depth, 1), d1), "d")
      call check(nf90_def_dim(nc, "y", size(depth, 2), d2), "d")
      call check(nf90_def_var(nc, "depth", nf90_double, [d1, d2], v1), "v")
      call check(nf90_enddef(nc), "e")
      call check(nf90_put_var(nc, v1, depth), "p")
      call check(nf90_close(nc), "c")
   end subroutine write_depth

   subroutine write_source(path, slon, slat, stime, src)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: slon(:), slat(:), stime(:), src(:, :, :)
      integer :: nc, d1, d2, d3, v1, v2, v3, vt
      call check(nf90_create(path, ior(nf90_clobber, nf90_64bit_offset), nc), "create src")
      call check(nf90_def_dim(nc, "lon", size(slon), d1), "d")
      call check(nf90_def_dim(nc, "lat", size(slat), d2), "d")
      call check(nf90_def_dim(nc, "time", size(stime), d3), "d")
      call check(nf90_def_var(nc, "lon", nf90_double, [d1], v1), "v")
      call check(nf90_def_var(nc, "lat", nf90_double, [d2], v2), "v")
      call check(nf90_def_var(nc, "time", nf90_double, [d3], vt), "v")
      call check(nf90_put_att(nc, vt, "units", "days since 1900-01-01 00:00:00"), "a")
      call check(nf90_def_var(nc, "friver", nf90_double, [d1, d2, d3], v3), "v")
      call check(nf90_enddef(nc), "e")
      call check(nf90_put_var(nc, v1, slon), "p")
      call check(nf90_put_var(nc, v2, slat), "p")
      call check(nf90_put_var(nc, vt, stime), "p")
      call check(nf90_put_var(nc, v3, src), "p")
      call check(nf90_close(nc), "c")
   end subroutine write_source

   subroutine read_out(path, nx, ny, nt, out)
      character(len=*), intent(in) :: path
      integer, intent(in) :: nx, ny, nt
      real(dp), allocatable, intent(out) :: out(:, :, :)
      integer :: nc, v1
      allocate (out(nx, ny, nt))
      call check(nf90_open(path, nf90_nowrite, nc), "open out")
      call check(nf90_inq_varid(nc, "runoff", v1), "runoff")
      call check(nf90_get_var(nc, v1, out), "read runoff")
      call check(nf90_close(nc), "close")
   end subroutine read_out

   subroutine placement_selftest()
      !! Known-answer placement on a 5-degree mosaic (continent lon [100,160) x lat [-20,45)):
      !! (a) an open-ocean source cell lands ENTIRELY in its own model cell; (b) a source cell over
      !! land next to the west coast lands entirely on the adjacent wet coastal cell; (c) on a
      !! REGIONAL 0..90E, 40S..40N mosaic a source cell in the domain is kept and one at 200E is
      !! discarded, the discarded amount being exactly its mass; `--global` refuses the latter.
      integer, parameter :: NX = 72, NY = 32, NSLON = 180, NSLAT = 90
      character(len=*), parameter :: F_HG = "om_runoff_selftest_hgrid.nc"
      character(len=*), parameter :: F_DP = "om_runoff_selftest_depth.nc"
      character(len=*), parameter :: F_SR = "om_runoff_selftest_src.nc"
      character(len=*), parameter :: F_OUT = "om_runoff_selftest_out.nc"
      real(dp), allocatable :: depth(:, :), slon(:), slat(:), src(:, :, :), out(:, :, :), disc(:)
      real(dp) :: m_src, area_c, expect_disc
      integer :: i, j, ii, jj, nreg_x, nreg_y
      logical :: land

      allocate (slon(NSLON), slat(NSLAT), src(NSLON, NSLAT, 1), depth(NX, NY))
      do i = 1, NSLON
         slon(i) = (real(i, dp) - 0.5_dp)*2.0_dp
      end do
      do j = 1, NSLAT
         slat(j) = -90.0_dp + (real(j, dp) - 0.5_dp)*2.0_dp
      end do
      call write_mosaic(F_HG, NX, NY, 0.0_dp, -80.0_dp, 5.0_dp)
      do j = 1, NY
         do i = 1, NX
            land = cell_is_land((real(i, dp) - 0.5_dp)*5.0_dp, -80.0_dp + (real(j, dp) - 0.5_dp)*5.0_dp)
            depth(i, j) = merge(0.0_dp, 4000.0_dp, land)
         end do
      end do
      call write_depth(F_DP, depth)

      ! (a) source cell lon [10,12], lat [-10,-8]  ->  model cell (3, 15) = lon [10,15], lat [-10,-5]
      src = 0.0_dp
      ii = 6
      jj = 41
      src(ii, jj, 1) = 3.0e-5_dp
      call write_source(F_SR, slon, slat, [1.0_dp], src)
      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4, 0.0_dp, .true., disc)
      call read_out(F_OUT, NX, NY, 1, out)
      m_src = src(ii, jj, 1)*R_EARTH**2*(2.0_dp*DEG)*(sin(-8.0_dp*DEG) - sin(-10.0_dp*DEG))
      area_c = R_EARTH**2*(5.0_dp*DEG)*(sin(-5.0_dp*DEG) - sin(-10.0_dp*DEG))
      if (count(out(:, :, 1) /= 0.0_dp) /= 1) call fail("selftest(a): the open-ocean source did not land in ONE cell")
      if (abs(out(3, 15, 1)*area_c - m_src) > 1.0e-10_dp*m_src) call fail("selftest(a): wrong cell or wrong mass")

      ! (b) source cell lon [100,102], lat [10,12]: over land cell (21, 19); nearest wet coastal cell (20, 19)
      src = 0.0_dp
      ii = 51
      jj = 51
      src(ii, jj, 1) = 3.0e-5_dp
      call write_source(F_SR, slon, slat, [1.0_dp], src)
      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4, 0.0_dp, .true., disc)
      call read_out(F_OUT, NX, NY, 1, out)
      m_src = src(ii, jj, 1)*R_EARTH**2*(2.0_dp*DEG)*(sin(12.0_dp*DEG) - sin(10.0_dp*DEG))
      area_c = R_EARTH**2*(5.0_dp*DEG)*(sin(15.0_dp*DEG) - sin(10.0_dp*DEG))
      if (depth(21, 19) > 0.0_dp .or. depth(20, 19) <= 0.0_dp) call fail("selftest(b): the geography is not as assumed")
      if (count(out(:, :, 1) /= 0.0_dp) /= 1) call fail("selftest(b): the land source did not land in ONE cell")
      if (abs(out(20, 19, 1)*area_c - m_src) > 1.0e-10_dp*m_src) call fail("selftest(b): not on the adjacent coastal cell")

      ! (c) regional mosaic lon [0,90), lat [-40,40): one source cell inside, one at lon 200
      nreg_x = 18
      nreg_y = 16
      call write_mosaic(F_HG, nreg_x, nreg_y, 0.0_dp, -40.0_dp, 5.0_dp)
      deallocate (depth)
      allocate (depth(nreg_x, nreg_y), source=4000.0_dp)
      call write_depth(F_DP, depth)
      src = 0.0_dp
      src(6, 46, 1) = 2.0e-5_dp
      src(101, 46, 1) = 5.0e-5_dp
      call write_source(F_SR, slon, slat, [1.0_dp], src)
      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4, 0.0_dp, .false., disc)
      expect_disc = src(101, 46, 1)*R_EARTH**2*(2.0_dp*DEG)*(sin((slat(46) + 1.0_dp)*DEG) - sin((slat(46) - 1.0_dp)*DEG))
      if (abs(disc(1) - expect_disc) > 1.0e-10_dp*expect_disc) call fail("selftest(c): wrong discarded mass")
      call read_out(F_OUT, nreg_x, nreg_y, 1, out)
      m_src = src(6, 46, 1)*R_EARTH**2*(2.0_dp*DEG)*(sin((slat(46) + 1.0_dp)*DEG) - sin((slat(46) - 1.0_dp)*DEG))
      area_c = R_EARTH**2*(5.0_dp*DEG)*(sin(5.0_dp*DEG) - sin(0.0_dp))
      if (count(out(:, :, 1) /= 0.0_dp) /= 1) call fail("selftest(c): the in-domain source did not land in ONE cell")
      if (abs(maxval(out(:, :, 1))*area_c - m_src) > 1.0e-10_dp*m_src) call fail("selftest(c): the in-domain mass was not kept")
      write (output_unit, "(a)") "selftest placement passed (open ocean, coast, regional discard)"
   end subroutine placement_selftest

   subroutine spread_case(depth, isrc, jsrc, rad, efold, max_rate, out, area_t)
      !! One spreading case on the 72 x 32 five-degree mosaic: a single 2-degree source cell inside model
      !! cell (isrc, jsrc), converted with the given spreading options; returns the model-grid runoff
      !! and the cell areas.
      real(dp), intent(in) :: depth(NXS, NYS)
      real(dp), intent(in) :: rad, efold, max_rate
      integer, intent(in) :: isrc, jsrc
      real(dp), intent(out) :: out(NXS, NYS), area_t(NXS, NYS)
      character(len=*), parameter :: F_HG = "om_runoff_selftest_hgrid.nc"
      character(len=*), parameter :: F_DP = "om_runoff_selftest_depth.nc"
      character(len=*), parameter :: F_SR = "om_runoff_selftest_src.nc"
      character(len=*), parameter :: F_OUT = "om_runoff_selftest_out.nc"
      real(dp) :: slon(NSX), slat(NSY), src(NSX, NSY, 1)
      real(dp), allocatable :: o3(:, :, :)
      integer :: i, j
      do i = 1, NSX
         slon(i) = (real(i, dp) - 0.5_dp)*2.0_dp
      end do
      do j = 1, NSY
         slat(j) = -90.0_dp + (real(j, dp) - 0.5_dp)*2.0_dp
      end do
      src = 0.0_dp
      src((isrc - 1)*5/2 + 2, ((jsrc - 1)*5 + 10)/2 + 2, 1) = 1.0e-5_dp
      ! (isrc-1)*5/2 and the latitude index are integer divisions: the 2-degree cell lies inside the 5-degree one
      call write_mosaic(F_HG, NXS, NYS, 0.0_dp, -80.0_dp, 5.0_dp)
      call write_depth(F_DP, depth)
      call write_source(F_SR, slon, slat, [15.0_dp], src)
      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4, 0.0_dp, .true., &
                   spread_radius_km=rad, spread_efold_km=efold, spread_max_rate=max_rate, spread_hard_km=3000.0_dp)
      call read_out(F_OUT, NXS, NYS, 1, o3)
      out = o3(:, :, 1)
      do j = 1, NYS
         do i = 1, NXS
            area_t(i, j) = R_EARTH**2*(5.0_dp*DEG)*(sin((-80.0_dp + 5.0_dp*real(j, dp))*DEG) &
                                                    - sin((-80.0_dp + 5.0_dp*real(j - 1, dp))*DEG))
         end do
      end do
   end subroutine spread_case

   subroutine spread_selftest()
      !! Spreading: an open coast spreads symmetrically and conserves; a full land wall is not crossed;
      !! a closed one-cell bay keeps everything; the periodic seam is crossed; radius 0 is the unspread
      !! answer bit for bit; the max-rate option widens; the flood-fill neighbours wrap and fold.
      real(dp) :: m0, m1, m2
      real(dp) :: depth(NXS, NYS), o0(NXS, NYS), o1(NXS, NYS), a(NXS, NYS)
      integer :: n
      integer :: ni(4), nj(4)
      logical :: ok

      ! A: coast at i = 11 (land i <= 10); source in (11,16)
      depth = 4000.0_dp
      depth(1:10, :) = 0.0_dp
      call spread_case(depth, 11, 16, 0.0_dp, 300.0_dp, 0.0_dp, o0, a)
      call spread_case(depth, 11, 16, 1500.0_dp, 300.0_dp, 0.0_dp, o1, a)
      m0 = sum(o0*a)
      m1 = sum(o1*a)
      if (abs(m1 - m0) > TOL_TOTAL*m0) call fail("selftest spread: the total changed")
      if (count(o1 /= 0.0_dp) < 5) call fail("selftest spread: the open-coast source did not spread")
      if (abs(o1(11, 15)*a(11, 15) - o1(11, 17)*a(11, 17)) > 1.0e-9_dp*m0) then
         call fail("selftest spread: the coastal spread is not symmetric")
      end if
      if (o1(11, 16) >= o0(11, 16)) call fail("selftest spread: the source cell was not relieved")
      if (any(o1(1:10, :) /= 0.0_dp)) call fail("selftest spread: mass reached land")

      ! radius 0 is bitwise the unspread answer (the option absent)
      call spread_case(depth, 11, 16, 0.0_dp, 300.0_dp, 0.0_dp, o0, a)
      call spread_case(depth, 11, 16, -0.0_dp, 123.0_dp, 0.0_dp, o1, a)
      if (any(o0 /= o1)) call fail("selftest spread: radius 0 is not the unspread answer")

      ! B: full land walls at i = 20 and i = 60 (the world is periodic): the 21..59 basin is a separate sea
      depth = 4000.0_dp
      depth(20, :) = 0.0_dp
      depth(60, :) = 0.0_dp
      call spread_case(depth, 19, 16, 3000.0_dp, 400.0_dp, 0.0_dp, o1, a)
      if (any(o1(21:59, :) /= 0.0_dp)) call fail("selftest spread: runoff crossed a land wall")
      if (o1(18, 16) <= 0.0_dp) call fail("selftest spread: no spread on the source's own side of the wall")

      ! C: a closed one-cell bay keeps its runoff
      depth = 4000.0_dp
      depth(39, 16) = 0.0_dp
      depth(41, 16) = 0.0_dp
      depth(40, 15) = 0.0_dp
      depth(40, 17) = 0.0_dp
      call spread_case(depth, 40, 16, 0.0_dp, 300.0_dp, 0.0_dp, o0, a)
      call spread_case(depth, 40, 16, 2000.0_dp, 300.0_dp, 0.0_dp, o1, a)
      if (count(o1 /= 0.0_dp) /= 1) call fail("selftest spread: a closed bay leaked")
      m0 = o0(40, 16)*a(40, 16)
      if (abs(o1(40, 16)*a(40, 16) - m0) > TOL_TOTAL*m0) call fail("selftest spread: the bay lost mass")

      ! D: the periodic seam: a source in column 1 feeds column 72 as it feeds column 2
      depth = 4000.0_dp
      call spread_case(depth, 1, 16, 1000.0_dp, 300.0_dp, 0.0_dp, o1, a)
      if (o1(72, 16) <= 0.0_dp) call fail("selftest spread: no spread across the periodic seam")
      if (abs(o1(72, 16) - o1(2, 16)) > 1.0e-9_dp*o1(2, 16)) call fail("selftest spread: the seam is not symmetric")

      ! E: the max-rate option widens a source until the peak falls
      depth = 4000.0_dp
      depth(1:10, :) = 0.0_dp
      call spread_case(depth, 11, 16, 600.0_dp, 800.0_dp, 0.0_dp, o0, a)
      call spread_case(depth, 11, 16, 600.0_dp, 800.0_dp, maxval(o0)*86400.0_dp/RHO_FW*0.25_dp, o1, a)
      if (maxval(o1) >= 0.9_dp*maxval(o0)) call fail("selftest spread: --spread-max-rate did not widen")
      m1 = sum(o1*a)
      m2 = sum(o0*a)
      if (abs(m1 - m2) > TOL_TOTAL*m2) call fail("selftest spread: widening changed the total")

      ! F: flood-fill neighbours: wrap in i only when periodic, mirror at the north row only with a fold
      call spread_nbrs(1, 5, 10, 8, .true., .false., n, ni, nj)
      ok = any(ni(1:n) == 10 .and. nj(1:n) == 5)
      call spread_nbrs(1, 5, 10, 8, .false., .false., n, ni, nj)
      ok = ok .and. .not. any(ni(1:n) == 10)
      call spread_nbrs(3, 8, 10, 8, .true., .true., n, ni, nj)
      ok = ok .and. any(ni(1:n) == 8 .and. nj(1:n) == 8)
      call spread_nbrs(3, 8, 10, 8, .true., .false., n, ni, nj)
      ok = ok .and. n == 3
      if (.not. ok) call fail("selftest spread: wrong flood-fill neighbours")
      write (output_unit, "(a)") "selftest spreading passed"
   end subroutine spread_selftest

   subroutine selftest()
      !! A synthetic 5-degree mosaic (lon 0..360, lat -80..80) with a continent and an island, a
      !! 2-degree source that puts runoff on the continent, on its coast, in the open ocean and in
      !! the polar cap beyond the mosaic (discarded, and checked against its known mass), the
      !! conversion's own assertions plus an independent re-read of the output, then the
      !! known-answer placement cases.
      integer, parameter :: NX = 72, NY = 32, NSLON = 180, NSLAT = 90, NT = 3
      character(len=*), parameter :: F_HG = "om_runoff_selftest_hgrid.nc"
      character(len=*), parameter :: F_DP = "om_runoff_selftest_depth.nc"
      character(len=*), parameter :: F_SR = "om_runoff_selftest_src.nc"
      character(len=*), parameter :: F_OUT = "om_runoff_selftest_out.nc"
      real(dp), allocatable :: depth(:, :), area_t(:, :)
      real(dp), allocatable :: slon(:), slat(:), stime(:), src(:, :, :), out(:, :, :), disc(:)
      real(dp) :: before, after, a_src, cap
      integer :: i, j, t, n_redir
      real(dp) :: redir_km
      logical :: land

      call write_mosaic(F_HG, NX, NY, 0.0_dp, -80.0_dp, 5.0_dp)
      allocate (depth(NX, NY), area_t(NX, NY))
      do j = 1, NY
         do i = 1, NX
            land = cell_is_land((real(i, dp) - 0.5_dp)*5.0_dp, -80.0_dp + (real(j, dp) - 0.5_dp)*5.0_dp)
            depth(i, j) = merge(0.0_dp, 4000.0_dp, land)
            area_t(i, j) = R_EARTH**2*(5.0_dp*DEG)*(sin((-80.0_dp + 5.0_dp*real(j, dp))*DEG) &
                                                    - sin((-80.0_dp + 5.0_dp*real(j - 1, dp))*DEG))
         end do
      end do
      call write_depth(F_DP, depth)

      ! source: 2-degree, three records, positive on the continent, its coast, an island, open ocean, the cap
      allocate (slon(NSLON), slat(NSLAT), stime(NT), src(NSLON, NSLAT, NT))
      do i = 1, NSLON
         slon(i) = (real(i, dp) - 0.5_dp)*2.0_dp
      end do
      do j = 1, NSLAT
         slat(j) = -90.0_dp + (real(j, dp) - 0.5_dp)*2.0_dp
      end do
      stime = [15.0_dp, 45.0_dp, 75.0_dp]
      src = 0.0_dp
      do t = 1, NT
         do j = 1, NSLAT
            do i = 1, NSLON
            if (cell_is_land(slon(i), slat(j)) .or. (abs(slon(i) - 130.0_dp) < 40.0_dp .and. abs(slat(j) - 12.0_dp) < 40.0_dp)) then
               src(i, j, t) = 1.0e-5_dp*real(t, dp)*real(1 + mod(3*i + 7*j, 5), dp)
            end if
            if (abs(slat(j)) > 85.0_dp .and. mod(i, 9) == 0) src(i, j, t) = 2.0e-5_dp*real(t, dp)
            end do
         end do
      end do
      src(10, 20, 2) = 1.0e30_dp
      call write_source(F_SR, slon, slat, stime, src)

      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4, 0.0_dp, .false., disc, &
                   n_redir, redir_km)
      ! the continent carries mass, so land cells were redirected; the nearest coast of a 5-degree mosaic's
      ! continent is never farther than the continent is wide (about 40 degrees = 4500 km), and is not zero
      if (n_redir < 10) call fail("selftest: too few land cells were redirected")
      if (redir_km <= 0.0_dp .or. redir_km > 4500.0_dp) call fail("selftest: implausible farthest redirect distance")

      ! independent verification from the files
      call read_out(F_OUT, NX, NY, NT, out)
      do t = 1, NT
         after = 0.0_dp
         before = 0.0_dp
         cap = 0.0_dp
         do j = 1, NSLAT
            do i = 1, NSLON
               if (abs(src(i, j, t)) < FILL_LIMIT) then
                  a_src = src(i, j, t)*R_EARTH**2*(2.0_dp*DEG)* &
                          (sin(min(90.0_dp, slat(j) + 1.0_dp)*DEG) - sin(max(-90.0_dp, slat(j) - 1.0_dp)*DEG))
                  before = before + a_src
                  ! the polar caps lie poleward of the 80-degree mosaic edge: the known discard
                  if (abs(slat(j)) > 85.0_dp) cap = cap + a_src
               end if
            end do
         end do
         do j = 1, NY
            do i = 1, NX
               if (depth(i, j) <= 0.0_dp .and. out(i, j, t) /= 0.0_dp) call fail("selftest: a land cell is non-zero")
               if (out(i, j, t) < 0.0_dp) call fail("selftest: negative runoff")
               after = after + out(i, j, t)*area_t(i, j)
            end do
         end do
         if (abs(after - (before - cap)) > TOL_TOTAL*before) call fail("selftest: re-read total /= source - polar cap")
         if (abs(disc(t) - cap) > TOL_TOTAL*before) call fail("selftest: the discarded mass is not the polar-cap mass")
         write (output_unit, "(a,i0,a,es22.14,a,es22.14,a,es22.14)") "selftest record ", t, ": source ", before, &
            ", re-read output ", after, ", discarded (polar cap) ", cap
      end do
      if (count(out(:, :, 1) > 0.0_dp) < 20) call fail("selftest: the runoff reached suspiciously few cells")
      call placement_selftest()
      call spread_selftest()
      write (output_unit, "(a)") "selftest passed"
   end subroutine selftest

   pure function cell_is_land(lon, lat) result(land)
      !! The synthetic geography of the self-test (degrees).
      real(dp), intent(in) :: lon, lat
      logical :: land
      land = (lon >= 100.0_dp .and. lon < 160.0_dp .and. lat >= -20.0_dp .and. lat < 45.0_dp) &
             .or. (abs(lon - 252.5_dp) < 2.6_dp .and. abs(lat - 12.5_dp) < 2.6_dp)
   end function cell_is_land

end program om_runoff_to_model
