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
   !!    contains the point; the nearest centre when none does), and every sub-cell hands its
   !!    equal share of the source mass to the cell it falls in.  Every share goes somewhere, so
   !!    the global total is preserved exactly (to round-off) whatever the geometry, including a
   !!    tripolar cap or a regional domain (mass outside the domain is handed to the nearest cell
   !!    and reported);
   !! 3. moves every share that landed on a model LAND cell to the nearest wet COASTAL cell (a
   !!    wet cell with a land 4-neighbour; any wet cell if the model has none) by great-circle
   !!    distance between cell centres, the role of the coupled models' runoff-to-coast map;
   !! 4. divides by the target cell area and writes `runoff(x, y, time)` (kg m-2 s-1, float64) with
   !!    the source's time axis and `units` attribute, the layout the reader requires.
   !!
   !! It asserts that the global total after (sum of out * target area) equals the total before
   !! (sum of source * source area) to 1e-12 relative, that no land cell is non-zero, and prints
   !! both totals per record.  A missing value (fill, NaN, |x| > 1e10) in the source counts as zero
   !! and is reported.  No third-party package; the self-test builds a synthetic mosaic, depth map
   !! and runoff product in the working directory and checks the same assertions.
   !!
   !! Arguments (positional): hgrid.nc depth.nc src.nc src_var out.nc
   !! Options: --depth-var NAME (default depth; wet where depth > 0)  --out-var NAME (default runoff)
   !!          --lon-name NAME (default lon) --lat-name NAME (default lat) --time-name NAME (default time)
   !!          --msub N (default 4)
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

   type :: point_index_t
      !! Nearest-neighbour index over unit vectors: a cubic bin grid (no pole or dateline
      !! singularity), CSR bin lists.
      integer :: nb = 0
      integer :: n = 0
      integer, allocatable :: bstart(:), ids(:)
      real(dp), allocatable :: v(:, :)
   end type point_index_t

   character(len=1024) :: arg
   character(len=1024) :: p_hgrid, p_depth, p_src, p_out
   character(len=64) :: src_var, depth_var, out_var, lon_name, lat_name, time_name
   integer :: msub, iarg, nargs

   nargs = command_argument_count()
   depth_var = "depth"
   out_var = "runoff"
   lon_name = "lon"
   lat_name = "lat"
   time_name = "time"
   msub = 4

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
         "[--msub N] | --selftest"
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
                trim(p_out), trim(out_var), trim(lon_name), trim(lat_name), trim(time_name), msub)

contains

   subroutine check(status, what)
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a)") "om_runoff_to_model: "//what//": "//trim(nf90_strerror(status))
         stop 1
      end if
   end subroutine check

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

   function locate(idx, sg_x, sg_y, nx, ny, q) result(cell)
      !! Flat index (`i + (j - 1) nx`) of the model T cell holding unit vector `q`: the nearest
      !! centre's quadrilateral, then its eight grid neighbours', else that nearest centre.
      type(point_index_t), intent(in) :: idx
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: sg_x(2*nx + 1, 2*ny + 1), sg_y(2*nx + 1, 2*ny + 1)
      real(dp), intent(in) :: q(3)
      integer :: cell
      integer :: k0, i0, j0, di, dj, i, j
      real(dp) :: c4(3, 4)

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
            if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
            call cell_corners(sg_x, sg_y, nx, ny, i, j, c4)
            if (quad_contains(c4, q)) then
               cell = i + (j - 1)*nx
               return
            end if
         end do
      end do
   end function locate

   ! ====================================================================
   ! the conversion
   ! ====================================================================

   subroutine convert(p_hgrid, p_depth, depth_var, p_src, src_var, p_out, out_var, lon_name, lat_name, time_name, msub)
      character(len=*), intent(in) :: p_hgrid, p_depth, depth_var, p_src, src_var, p_out, out_var
      character(len=*), intent(in) :: lon_name, lat_name, time_name
      integer, intent(in) :: msub

      integer :: nx, ny, nxp, nyp, nlon, nlat, nt
      real(dp), allocatable :: sg_x(:, :), sg_y(:, :), area_sg(:, :), area_t(:, :), depth(:, :)
      real(dp), allocatable :: slon(:), slat(:), stime(:), asrc(:, :), f(:, :), mass(:), out(:, :, :)
      real(dp), allocatable :: cv(:, :), wetv(:, :)
      integer, allocatable :: tidx(:, :, :), redir(:), coast_id(:)
      logical, allocatable :: wet(:, :), hit(:), coastal(:, :)
      character(len=256) :: tunits, tcal
      logical :: has_cal
      type(point_index_t) :: cidx, ridx
      integer :: i, j, k, t, a, b, nwet, ncoast, nmoved_cells, nmiss
      real(dp) :: before, after, lon0, lat0, dlon, dlat
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
      wet = depth > 0.0_dp
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
                  k = locate(cidx, sg_x, sg_y, nx, ny, v)
                  tidx(a + (b - 1)*msub, i, j) = k
                  hit(k) = .true.
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
            if (i > 1) coastal(i, j) = coastal(i, j) .or. .not. wet(i - 1, j)
            if (i < nx) coastal(i, j) = coastal(i, j) .or. .not. wet(i + 1, j)
            if (j > 1) coastal(i, j) = coastal(i, j) .or. .not. wet(i, j - 1)
            if (j < ny) coastal(i, j) = coastal(i, j) .or. .not. wet(i, j + 1)
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
      do j = 1, nlat
         do i = 1, nlon
            do a = 1, msub*msub
               tidx(a, i, j) = redir(tidx(a, i, j))
            end do
         end do
      end do

      ! the records
      allocate (f(nlon, nlat), mass(nx*ny), out(nx, ny, max(nt, 1)))
      w_sub = 1.0_dp/real(msub*msub, dp)
      write (output_unit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "target ", nx, " x ", ny, " (", nwet, " wet, ", ncoast, &
         " coastal); source ", nlon
      write (output_unit, "(a,i0,a,i0,a,i0)") "  x ", nlat, ", records ", nt, "; land cells receiving mass moved to the coast: ", &
         nmoved_cells
      do t = 1, max(nt, 1)
         call read_source_record(p_src, src_var, nlon, nlat, nt, t, f, nmiss)
         mass = 0.0_dp
         before = 0.0_dp
         do j = 1, nlat
            do i = 1, nlon
               if (f(i, j) == 0.0_dp) cycle
               before = before + f(i, j)*asrc(i, j)
               do a = 1, msub*msub
                  k = tidx(a, i, j)
                  mass(k) = mass(k) + f(i, j)*asrc(i, j)*w_sub
               end do
            end do
         end do
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
         if (abs(after - before) > TOL_TOTAL*max(abs(before), tiny(1.0_dp))) then
            call fail("the global total was not preserved")
         end if
      end do
      write (output_unit, "(a)") "ok: global total preserved, no land cell non-zero"
      call write_output(p_out, out_var, nx, ny, max(nt, 1), stime, tunits, tcal, has_cal, out, p_src)
   end subroutine convert

   ! ====================================================================
   ! I/O
   ! ====================================================================

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

   subroutine selftest()
      !! A synthetic 5-degree mosaic (lon 0..360, lat -80..80) with a continent and an island, a
      !! 2-degree source that puts runoff on the continent, on its coast, in the open ocean and in
      !! the polar cap beyond the mosaic, and the conversion's own assertions plus an independent
      !! re-read of the output.
      integer, parameter :: NX = 72, NY = 32, NSLON = 180, NSLAT = 90, NT = 3
      character(len=*), parameter :: F_HG = "om_runoff_selftest_hgrid.nc"
      character(len=*), parameter :: F_DP = "om_runoff_selftest_depth.nc"
      character(len=*), parameter :: F_SR = "om_runoff_selftest_src.nc"
      character(len=*), parameter :: F_OUT = "om_runoff_selftest_out.nc"
      real(dp), allocatable :: sg_x(:, :), sg_y(:, :), area(:, :), depth(:, :)
      real(dp), allocatable :: slon(:), slat(:), stime(:), src(:, :, :), out(:, :, :), area_t(:, :)
      real(dp) :: dl, sn, ss, before, after
      integer :: i, j, t, nc, d1, d2, d3, v1, v2, v3, vt
      logical :: land

      ! mosaic
      allocate (sg_x(2*NX + 1, 2*NY + 1), sg_y(2*NX + 1, 2*NY + 1), area(2*NX, 2*NY))
      dl = 360.0_dp/real(2*NX, dp)
      do j = 1, 2*NY + 1
         do i = 1, 2*NX + 1
            sg_x(i, j) = real(i - 1, dp)*dl
            sg_y(i, j) = -80.0_dp + real(j - 1, dp)*160.0_dp/real(2*NY, dp)
         end do
      end do
      do j = 1, 2*NY
         ss = sin(sg_y(1, j)*DEG)
         sn = sin(sg_y(1, j + 1)*DEG)
         do i = 1, 2*NX
            area(i, j) = R_EARTH**2*dl*DEG*(sn - ss)
         end do
      end do
      call check(nf90_create(F_HG, ior(nf90_clobber, nf90_64bit_offset), nc), "create hgrid")
      call check(nf90_def_dim(nc, "nxp", 2*NX + 1, d1), "d")
      call check(nf90_def_dim(nc, "nyp", 2*NY + 1, d2), "d")
      call check(nf90_def_dim(nc, "nx", 2*NX, d3), "d")
      call check(nf90_def_var(nc, "x", nf90_double, [d1, d2], v1), "v")
      call check(nf90_def_var(nc, "y", nf90_double, [d1, d2], v2), "v")
      call check(nf90_def_dim(nc, "ny", 2*NY, d1), "d")
      call check(nf90_def_var(nc, "area", nf90_double, [d3, d1], v3), "v")
      call check(nf90_enddef(nc), "e")
      call check(nf90_put_var(nc, v1, sg_x), "p")
      call check(nf90_put_var(nc, v2, sg_y), "p")
      call check(nf90_put_var(nc, v3, area), "p")
      call check(nf90_close(nc), "c")

      ! depth: a continent [100,160) x [-20,45), an island at lon 250 / lat 10, ocean 4000 m elsewhere
      allocate (depth(NX, NY), area_t(NX, NY))
      do j = 1, NY
         do i = 1, NX
            land = cell_is_land((real(i, dp) - 0.5_dp)*5.0_dp, -80.0_dp + (real(j, dp) - 0.5_dp)*5.0_dp)
            depth(i, j) = merge(0.0_dp, 4000.0_dp, land)
            area_t(i, j) = area(2*i - 1, 2*j - 1) + area(2*i, 2*j - 1) + area(2*i - 1, 2*j) + area(2*i, 2*j)
         end do
      end do
      call check(nf90_create(F_DP, ior(nf90_clobber, nf90_64bit_offset), nc), "create depth")
      call check(nf90_def_dim(nc, "x", NX, d1), "d")
      call check(nf90_def_dim(nc, "y", NY, d2), "d")
      call check(nf90_def_var(nc, "depth", nf90_double, [d1, d2], v1), "v")
      call check(nf90_enddef(nc), "e")
      call check(nf90_put_var(nc, v1, depth), "p")
      call check(nf90_close(nc), "c")

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
      call check(nf90_create(F_SR, ior(nf90_clobber, nf90_64bit_offset), nc), "create src")
      call check(nf90_def_dim(nc, "lon", NSLON, d1), "d")
      call check(nf90_def_dim(nc, "lat", NSLAT, d2), "d")
      call check(nf90_def_dim(nc, "time", NT, d3), "d")
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

      call convert(F_HG, F_DP, "depth", F_SR, "friver", F_OUT, "runoff", "lon", "lat", "time", 4)

      ! independent verification from the files
      allocate (out(NX, NY, NT))
      call check(nf90_open(F_OUT, nf90_nowrite, nc), "open out")
      call check(nf90_inq_varid(nc, "runoff", v1), "runoff")
      call check(nf90_get_var(nc, v1, out), "read runoff")
      call check(nf90_close(nc), "close")
      do t = 1, NT
         after = 0.0_dp
         before = 0.0_dp
         do j = 1, NSLAT
            do i = 1, NSLON
               if (abs(src(i, j, t)) < FILL_LIMIT) then
                  before = before + src(i, j, t)*R_EARTH**2*(2.0_dp*DEG)* &
                           (sin(min(90.0_dp, slat(j) + 1.0_dp)*DEG) - sin(max(-90.0_dp, slat(j) - 1.0_dp)*DEG))
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
         if (abs(after - before) > TOL_TOTAL*before) call fail("selftest: re-read total differs from the source total")
         write (output_unit, "(a,i0,a,es22.14,a,es22.14)") "selftest record ", t, ": source ", before, ", re-read output ", after
      end do
      if (count(out(:, :, 1) > 0.0_dp) < 20) call fail("selftest: the runoff reached suspiciously few cells")
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
