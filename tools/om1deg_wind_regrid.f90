!! JRA55-do near-surface wind -> OM_1deg C-grid wind stress (offline helper).
program om1deg_wind_regrid
   !! Offline preprocessing helper for
   !! `validation_examples/ocean/global_1deg/global_1deg_wind.nml`, driven by
   !! `tools/om1deg_prepare_wind.py` (which builds it with the system Fortran
   !! compiler + `nf-config` and passes the arguments).  Not part of the model.
   !!
   !! For every 3-hourly JRA55-do record of the requested year:
   !!
   !! 1. bilinear interpolation of `uas`/`vas` (10 m wind, m/s, east/north,
   !!    on the TL319 640x320 lon/Gaussian-lat grid) to the model's C-grid
   !!    u-face points (supergrid node `(2i-1, 2j)`, the WEST face of T-cell
   !!    `(i,j)`, i = 1..nx+1) and v-face points (node `(2i, 2j-1)`, the SOUTH
   !!    face, j = 1..ny+1);
   !! 2. quadratic bulk stress in geographic axes,
   !!    `tau = rho_air * C_d(|U|) * |U| * (u_E, v_N)`, with `rho_air = 1.22`
   !!    kg/m^3 and C_d either the Large & Yeager (2004) neutral 10 m drag
   !!    coefficient `C_d = (2.7/U + 0.142 + 0.0764 U) * 1e-3`, U floored at
   !!    0.5 m/s (the NCAR/FMS `ncar_ocean_fluxes` floor), or a constant;
   !!    the ocean surface velocity is ignored (absolute, not relative, wind);
   !! 3. rotation onto the grid axes by the grid angle `a` at that face
   !!    (`u_grid = cos a u_E + sin a v_N`, `v_grid = -sin a u_E + cos a v_N`,
   !!    the `ocean_metrics_t%angle_dx` convention), `a` from the node
   !!    geography (`face_angle` says why not from the mosaic's `angle_dx`);
   !! 4. accumulation into `avg_hours`-long bins with trapezoidal weights: a
   !!    record exactly on a bin edge counts half to each side, so a daily bin
   !!    is centred on 12:00 when the next day's 00:00 record exists (the
   !!    `padded` JRA files carry one record past the year end).
   !!
   !! Output: NetCDF classic (64-bit offset) in the layout
   !! `&ocean_dataovr_nml` reads — `taux(xf=nx+1, y=ny, time)`,
   !! `tauy(x=nx, yf=ny+1, time)`, Fortran order, float32, `time` in
   !! "days since <year>-01-01 00:00:00" at the bin centres.
   !!
   !! Arguments (positional):
   !!   hgrid uas_file vas_file out_file year cd_form avg_hours
   !! with `cd_form` either `ly04` or a constant drag coefficient (e.g. 1.2e-3).
   use, intrinsic :: iso_fortran_env, only: real32, real64, error_unit, output_unit
   use netcdf, only: nf90_open, nf90_close, nf90_create, nf90_enddef, nf90_noerr, &
                     nf90_nowrite, nf90_clobber, nf90_64bit_offset, nf90_float, &
                     nf90_double, nf90_global, nf90_inq_varid, nf90_inq_dimid, &
                     nf90_inquire_dimension, nf90_get_var, nf90_put_var, nf90_def_dim, &
                     nf90_def_var, nf90_put_att, nf90_get_att, nf90_strerror
   implicit none

   integer, parameter :: dp = real64
   real(dp), parameter :: RHO_AIR = 1.22_dp
      !! Air density (kg/m^3) for the bulk stress.
   real(dp), parameter :: U_FLOOR = 0.5_dp
      !! Low-wind floor (m/s) of the LY04 drag law, as in NCAR/FMS.
   real(dp), parameter :: DEG = acos(-1.0_dp)/180.0_dp
   real(dp), parameter :: FILL_LIMIT = 1.0e10_dp
      !! Any |value| above this is a missing value: abort, never interpolate it.

   character(len=1024) :: hgrid_path, uas_path, vas_path, out_path, arg
   character(len=64) :: cd_form
   integer :: year, avg_hours
   real(dp) :: cd_const
   logical :: use_ly04

   ! Supergrid
   integer :: nxp, nyp, nx, ny
   real(dp), allocatable :: sg_x(:, :), sg_y(:, :), sg_ang(:, :)

   ! JRA source grid
   integer :: nlon, nlat, ntime
   real(dp), allocatable :: jlon(:), jlat(:), jtime(:)
   real(real32), allocatable :: ua(:, :), va(:, :)

   ! Interpolation stencils for the u-face and v-face target sets
   integer, allocatable :: ui0(:, :), ui1(:, :), uj0(:, :), uj1(:, :)
   real(dp), allocatable :: uwx(:, :), uwy(:, :), ucos(:, :), usin(:, :)
   integer, allocatable :: vi0(:, :), vi1(:, :), vj0(:, :), vj1(:, :)
   real(dp), allocatable :: vwx(:, :), vwy(:, :), vcos(:, :), vsin(:, :)

   ! Accumulators, one slab per output bin
   integer :: nbins
   real(dp), allocatable :: taux(:, :, :), tauy(:, :, :), wsum(:)

   integer :: ncu, ncv, varu, varv, rec, b, nused
   real(dp) :: t_year0, t_rel, bin_len, pos, frac, taumax, t_first, t_last

   ! ------------------------------------------------------------------
   ! Arguments
   ! ------------------------------------------------------------------
   if (command_argument_count() /= 7) then
      write (error_unit, "(a)") "usage: om1deg_wind_regrid hgrid uas vas out year cd_form avg_hours"
      stop 2
   end if
   call get_command_argument(1, hgrid_path)
   call get_command_argument(2, uas_path)
   call get_command_argument(3, vas_path)
   call get_command_argument(4, out_path)
   call get_command_argument(5, arg)
   read (arg, *) year
   call get_command_argument(6, cd_form)
   call get_command_argument(7, arg)
   read (arg, *) avg_hours
   use_ly04 = (trim(cd_form) == "ly04")
   cd_const = 0.0_dp
   if (.not. use_ly04) read (cd_form, *) cd_const
   if (avg_hours <= 0 .or. mod(24, avg_hours) /= 0) then
      write (error_unit, "(a)") "avg_hours must divide 24"
      stop 2
   end if

   ! ------------------------------------------------------------------
   ! Supergrid: node lon/lat + angle_dx
   ! ------------------------------------------------------------------
   call read_supergrid(trim(hgrid_path))
   nx = (nxp - 1)/2
   ny = (nyp - 1)/2
   write (output_unit, "(a,i0,a,i0,a,i0,a,i0)") "supergrid ", nxp, " x ", nyp, &
      " nodes -> model ", nx, " x ", ny

   ! ------------------------------------------------------------------
   ! JRA grid + time axis (uas and vas must share them)
   ! ------------------------------------------------------------------
   call open_jra(trim(uas_path), "uas", ncu, varu)
   call open_jra(trim(vas_path), "vas", ncv, varv)
   write (output_unit, "(a,i0,a,i0,a,i0,a)") "JRA grid ", nlon, " x ", nlat, ", ", ntime, " records"

   ! ------------------------------------------------------------------
   ! Stencils: u faces (i = 1..nx+1, node 2i-1) x (j = 1..ny, node 2j);
   !           v faces (i = 1..nx, node 2i) x (j = 1..ny+1, node 2j-1)
   ! ------------------------------------------------------------------
   call build_stencil(nx + 1, ny, 1, 2, ui0, ui1, uj0, uj1, uwx, uwy, ucos, usin)
   call build_stencil(nx, ny + 1, 2, 1, vi0, vi1, vj0, vj1, vwx, vwy, vcos, vsin)

   ! ------------------------------------------------------------------
   ! Bins over the year
   ! ------------------------------------------------------------------
   t_year0 = real(days_from_civil(year, 1, 1) - days_from_civil(1900, 1, 1), dp)
   nbins = (days_from_civil(year + 1, 1, 1) - days_from_civil(year, 1, 1))*(24/avg_hours)
   bin_len = real(avg_hours, dp)/24.0_dp
   allocate (taux(nx + 1, ny, nbins), tauy(nx, ny + 1, nbins), wsum(nbins))
   taux = 0.0_dp
   tauy = 0.0_dp
   wsum = 0.0_dp
   allocate (ua(nlon, nlat), va(nlon, nlat))

   nused = 0
   t_first = huge(1.0_dp)
   t_last = -huge(1.0_dp)
   taumax = 0.0_dp
   do rec = 1, ntime
      t_rel = jtime(rec) - t_year0
      if (t_rel < -1.0e-6_dp .or. t_rel > nbins*bin_len + 1.0e-6_dp) cycle
      call check(nf90_get_var(ncu, varu, ua, start=[1, 1, rec], count=[nlon, nlat, 1]), "read uas")
      call check(nf90_get_var(ncv, varv, va, start=[1, 1, rec], count=[nlon, nlat, 1]), "read vas")
      if (maxval(abs(ua)) > FILL_LIMIT .or. maxval(abs(va)) > FILL_LIMIT) then
         write (error_unit, "(a,i0)") "missing values in JRA record ", rec
         stop 1
      end if
      nused = nused + 1
      t_first = min(t_first, jtime(rec))
      t_last = max(t_last, jtime(rec))
      pos = t_rel/bin_len
      b = nint(pos)
      if (abs(pos - b) < 1.0e-6_dp) then
         ! On a bin edge: half to the bin that ends here, half to the one that starts.
         if (b >= 1) call accumulate(b, 0.5_dp)
         if (b + 1 <= nbins) call accumulate(b + 1, 0.5_dp)
      else
         b = int(floor(pos)) + 1
         frac = 1.0_dp
         if (b >= 1 .and. b <= nbins) call accumulate(b, frac)
      end if
   end do
   call check(nf90_close(ncu), "close uas")
   call check(nf90_close(ncv), "close vas")

   write (output_unit, "(a,i0,a,f12.4,a,f12.4,a)") "used ", nused, " records (days since 1900: ", &
      t_first, " .. ", t_last, ")"
   if (any(wsum <= 0.0_dp)) then
      write (error_unit, "(a,i0,a)") "bins without data: ", count(wsum <= 0.0_dp), &
         " — the JRA files do not cover the requested year"
      stop 1
   end if
   write (output_unit, "(a,f8.3,a,f8.3)") "weight per bin: min ", minval(wsum), " max ", maxval(wsum)
   do b = 1, nbins
      taux(:, :, b) = taux(:, :, b)/wsum(b)
      tauy(:, :, b) = tauy(:, :, b)/wsum(b)
   end do
   write (output_unit, "(a,es12.4,a)") "max |tau| of any 3-hourly sample: ", taumax, " Pa"
   write (output_unit, "(a,es12.4,a,es12.4,a)") "bin-mean taux range ", minval(taux), " .. ", &
      maxval(taux), " Pa"
   write (output_unit, "(a,es12.4,a,es12.4,a)") "bin-mean tauy range ", minval(tauy), " .. ", &
      maxval(tauy), " Pa"

   call write_output(trim(out_path))
   write (output_unit, "(a)") "wrote "//trim(out_path)

contains

   subroutine check(status, what)
      !! Abort with the NetCDF message on any error.
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a)") "om1deg_wind_regrid: "//what//": "//trim(nf90_strerror(status))
         stop 1
      end if
   end subroutine check

   pure integer function days_from_civil(y, m, d) result(n)
      !! Days since 1970-01-01 of a proleptic-Gregorian date (H. Hinnant's algorithm).
      integer, intent(in) :: y, m, d
      integer :: yy, era, yoe, doy, doe
      yy = y
      if (m <= 2) yy = y - 1
      if (yy >= 0) then
         era = yy/400
      else
         era = (yy - 399)/400
      end if
      yoe = yy - era*400
      if (m > 2) then
         doy = (153*(m - 3) + 2)/5 + d - 1
      else
         doy = (153*(m + 9) + 2)/5 + d - 1
      end if
      doe = yoe*365 + yoe/4 - yoe/100 + doy
      n = era*146097 + doe - 719468
   end function days_from_civil

   subroutine read_supergrid(path)
      !! Node longitude/latitude (degrees) and `angle_dx` (degrees) of the mosaic.
      character(len=*), intent(in) :: path
      integer :: nc, vid, did
      call check(nf90_open(path, nf90_nowrite, nc), "open "//path)
      call check(nf90_inq_dimid(nc, "nxp", did), "nxp")
      call check(nf90_inquire_dimension(nc, did, len=nxp), "nxp len")
      call check(nf90_inq_dimid(nc, "nyp", did), "nyp")
      call check(nf90_inquire_dimension(nc, did, len=nyp), "nyp len")
      allocate (sg_x(nxp, nyp), sg_y(nxp, nyp), sg_ang(nxp, nyp))
      call check(nf90_inq_varid(nc, "x", vid), "x")
      call check(nf90_get_var(nc, vid, sg_x), "read x")
      call check(nf90_inq_varid(nc, "y", vid), "y")
      call check(nf90_get_var(nc, vid, sg_y), "read y")
      call check(nf90_inq_varid(nc, "angle_dx", vid), "angle_dx (read for the cross-check)")
      call check(nf90_get_var(nc, vid, sg_ang), "read angle_dx")
      call check(nf90_close(nc), "close hgrid")
   end subroutine read_supergrid

   subroutine open_jra(path, var, nc, vid)
      !! Open one JRA file; read (or cross-check) the lon/lat/time axes.
      character(len=*), intent(in) :: path, var
      integer, intent(out) :: nc, vid
      integer :: did, vt, nlon_f, nlat_f, nt_f
      real(dp), allocatable :: lon_f(:), lat_f(:), t_f(:)
      character(len=128) :: units
      call check(nf90_open(path, nf90_nowrite, nc), "open "//path)
      call check(nf90_inq_varid(nc, var, vid), "variable "//var)
      call check(nf90_inq_dimid(nc, "lon", did), "lon dim")
      call check(nf90_inquire_dimension(nc, did, len=nlon_f), "lon len")
      call check(nf90_inq_dimid(nc, "lat", did), "lat dim")
      call check(nf90_inquire_dimension(nc, did, len=nlat_f), "lat len")
      call check(nf90_inq_dimid(nc, "time", did), "time dim")
      call check(nf90_inquire_dimension(nc, did, len=nt_f), "time len")
      allocate (lon_f(nlon_f), lat_f(nlat_f), t_f(nt_f))
      call check(nf90_inq_varid(nc, "lon", vt), "lon var")
      call check(nf90_get_var(nc, vt, lon_f), "read lon")
      call check(nf90_inq_varid(nc, "lat", vt), "lat var")
      call check(nf90_get_var(nc, vt, lat_f), "read lat")
      call check(nf90_inq_varid(nc, "time", vt), "time var")
      call check(nf90_get_var(nc, vt, t_f), "read time")
      units = ""
      call check(nf90_get_att(nc, vt, "units", units), "time units")
      if (index(units, "days since 1900-01-01") /= 1) then
         write (error_unit, "(a)") "unexpected JRA time units: "//trim(units)
         stop 1
      end if
      if (.not. allocated(jlon)) then
         nlon = nlon_f
         nlat = nlat_f
         ntime = nt_f
         jlon = lon_f
         jlat = lat_f
         jtime = t_f
         if (any(jlat(2:) <= jlat(:nlat - 1))) then
            write (error_unit, "(a)") "JRA latitude must be strictly increasing"
            stop 1
         end if
      else
         if (nlon_f /= nlon .or. nlat_f /= nlat .or. nt_f /= ntime) then
            write (error_unit, "(a)") "uas and vas grids/time axes differ"
            stop 1
         end if
         if (maxval(abs(lon_f - jlon)) > 1.0e-9_dp .or. maxval(abs(lat_f - jlat)) > 1.0e-9_dp &
             .or. maxval(abs(t_f - jtime)) > 1.0e-9_dp) then
            write (error_unit, "(a)") "uas and vas grids/time axes differ"
            stop 1
         end if
      end if
   end subroutine open_jra

   subroutine build_stencil(ni, nj, oi, oj, i0, i1, j0, j1, wx, wy, ca, sa)
      !! Bilinear stencil on the JRA grid for target node `(2i - 2 + oi, 2j - 2 + oj)`,
      !! i = 1..ni, j = 1..nj, plus cos/sin of the grid angle there.  Longitude is
      !! periodic (JRA lon is uniform); latitude beyond the outermost Gaussian
      !! row takes that row.
      integer, intent(in) :: ni, nj, oi, oj
      integer, allocatable, intent(out) :: i0(:, :), i1(:, :), j0(:, :), j1(:, :)
      real(dp), allocatable, intent(out) :: wx(:, :), wy(:, :), ca(:, :), sa(:, :)
      integer :: i, j, si, sj, k, lo, hi, mid
      real(dp) :: dlon, x, lat, ang, dev, ang_dev
      allocate (i0(ni, nj), i1(ni, nj), j0(ni, nj), j1(ni, nj))
      allocate (wx(ni, nj), wy(ni, nj), ca(ni, nj), sa(ni, nj))
      ang_dev = 0.0_dp
      dlon = 360.0_dp/nlon
      if (maxval(abs(jlon(2:) - jlon(:nlon - 1) - dlon)) > 1.0e-6_dp) then
         write (error_unit, "(a)") "JRA longitude is not uniform over 360 degrees"
         stop 1
      end if
      do j = 1, nj
         do i = 1, ni
            si = 2*i - 2 + oi
            sj = 2*j - 2 + oj
            x = modulo(sg_x(si, sj) - jlon(1), 360.0_dp)/dlon
            k = int(floor(x))
            wx(i, j) = x - k
            i0(i, j) = modulo(k, nlon) + 1
            i1(i, j) = modulo(k + 1, nlon) + 1
            lat = sg_y(si, sj)
            if (lat <= jlat(1)) then
               j0(i, j) = 1
               j1(i, j) = 1
               wy(i, j) = 0.0_dp
            else if (lat >= jlat(nlat)) then
               j0(i, j) = nlat
               j1(i, j) = nlat
               wy(i, j) = 0.0_dp
            else
               lo = 1
               hi = nlat
               do while (hi - lo > 1)
                  mid = (lo + hi)/2
                  if (jlat(mid) <= lat) then
                     lo = mid
                  else
                     hi = mid
                  end if
               end do
               j0(i, j) = lo
               j1(i, j) = hi
               wy(i, j) = (lat - jlat(lo))/(jlat(hi) - jlat(lo))
            end if
            ang = face_angle(si, sj)
            ca(i, j) = cos(ang)
            sa(i, j) = sin(ang)
            if (sg_y(si, sj) < 59.0_dp) then
               dev = abs(modulo(ang/DEG - sg_ang(si, sj) + 180.0_dp, 360.0_dp) - 180.0_dp)
               ang_dev = max(ang_dev, dev)
            end if
         end do
      end do
      write (output_unit, "(a,f8.3,a)") "face rotation vs the mosaic angle_dx south of the cap (lat < 59): "// &
         "max |difference| ", ang_dev, " degrees"
   end subroutine build_stencil

   real(dp) function face_angle(si, sj) result(ang)
      !! Grid rotation (radians, counter-clockwise from true east) at supergrid node
      !! `(si, sj)`: the chord from node `si-1` to node `si+1` along the row
      !! (periodic in i), projected on the local east/north tangent plane in 3-D,
      !! so it stays exact next to the pole.  Computed here rather than read from
      !! the mosaic's `angle_dx`, which is not usable for this on OM_1deg's
      !! `ocean_hgrid.nc`: (a) in the bipolar Arctic cap it is the heading in
      !! DEGREE space, `atan2(dlat, dlon)` without the `cos(lat)` metric (e.g.
      !! 10.4 against the true 32.7 degrees at 73 N) — the error MOM6's
      !! `GRID_ROTATION_ANGLE_BUGS = False` exists to avoid; (b) its edge nodes
      !! are one-sided: column 1 differs from its periodic twin, column nxp,
      !! and the fold row is not antisymmetric about the fold, so the same
      !! physical face would get two different stresses.  South of the cap
      !! (the plain lon-lat band) both are zero.
      integer, intent(in) :: si, sj
      integer :: sw, se
      real(dp) :: pw(3), pe(3), d(3), east(3), north(3)
      real(dp) :: lam, phi
      sw = si - 1
      se = si + 1
      if (sw < 1) sw = nxp - 1
      if (se > nxp) se = 2
      pw = unit_vec(sg_x(sw, sj), sg_y(sw, sj))
      pe = unit_vec(sg_x(se, sj), sg_y(se, sj))
      d = pe - pw
      lam = sg_x(si, sj)*DEG
      phi = sg_y(si, sj)*DEG
      east = [-sin(lam), cos(lam), 0.0_dp]
      north = [-sin(phi)*cos(lam), -sin(phi)*sin(lam), cos(phi)]
      ang = atan2(dot_product(d, north), dot_product(d, east))
   end function face_angle

   pure function unit_vec(lon, lat) result(p)
      !! Unit-sphere Cartesian position of (lon, lat) in degrees.
      real(dp), intent(in) :: lon, lat
      real(dp) :: p(3)
      p = [cos(lat*DEG)*cos(lon*DEG), cos(lat*DEG)*sin(lon*DEG), sin(lat*DEG)]
   end function unit_vec

   pure real(dp) function drag_coeff(u) result(cd)
      !! Neutral 10 m drag coefficient: Large & Yeager (2004) eq. 6a
      !! (`(2.7/U + 0.142 + 0.0764 U) * 1e-3`, U floored at 0.5 m/s), or the constant.
      real(dp), intent(in) :: u
      real(dp) :: uf
      if (use_ly04) then
         uf = max(u, U_FLOOR)
         cd = (2.7_dp/uf + 0.142_dp + 0.0764_dp*uf)*1.0e-3_dp
      else
         cd = cd_const
      end if
   end function drag_coeff

   subroutine accumulate(bin, w)
      !! Add weight `w` of the current record's stress (both face sets) to `bin`.
      integer, intent(in) :: bin
      real(dp), intent(in) :: w
      integer :: i, j
      real(dp) :: ue, vn, spd, k, te, tn
      do j = 1, ny
         do i = 1, nx + 1
            call interp(ui0(i, j), ui1(i, j), uj0(i, j), uj1(i, j), uwx(i, j), uwy(i, j), ue, vn)
            spd = sqrt(ue*ue + vn*vn)
            k = RHO_AIR*drag_coeff(spd)*spd
            te = k*ue
            tn = k*vn
            taumax = max(taumax, k*spd)
            taux(i, j, bin) = taux(i, j, bin) + w*(ucos(i, j)*te + usin(i, j)*tn)
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            call interp(vi0(i, j), vi1(i, j), vj0(i, j), vj1(i, j), vwx(i, j), vwy(i, j), ue, vn)
            spd = sqrt(ue*ue + vn*vn)
            k = RHO_AIR*drag_coeff(spd)*spd
            te = k*ue
            tn = k*vn
            tauy(i, j, bin) = tauy(i, j, bin) + w*(-vsin(i, j)*te + vcos(i, j)*tn)
         end do
      end do
      wsum(bin) = wsum(bin) + w
   end subroutine accumulate

   subroutine interp(a0, a1, b0, b1, fx, fy, ue, vn)
      !! Bilinear value of the current `ua`/`va` record at one target.
      integer, intent(in) :: a0, a1, b0, b1
      real(dp), intent(in) :: fx, fy
      real(dp), intent(out) :: ue, vn
      ue = (1.0_dp - fy)*((1.0_dp - fx)*ua(a0, b0) + fx*ua(a1, b0)) &
           + fy*((1.0_dp - fx)*ua(a0, b1) + fx*ua(a1, b1))
      vn = (1.0_dp - fy)*((1.0_dp - fx)*va(a0, b0) + fx*va(a1, b0)) &
           + fy*((1.0_dp - fx)*va(a0, b1) + fx*va(a1, b1))
   end subroutine interp

   subroutine write_output(path)
      !! The `&ocean_dataovr_nml` layout (see the program docstring).
      character(len=*), intent(in) :: path
      integer :: nc, d_x, d_y, d_xf, d_yf, d_t, v_t, v_tx, v_ty, b2
      character(len=64) :: tunits, cdtxt
      call check(nf90_create(path, ior(nf90_clobber, nf90_64bit_offset), nc), "create "//path)
      call check(nf90_def_dim(nc, "time", nbins, d_t), "def time")
      call check(nf90_def_dim(nc, "y", ny, d_y), "def y")
      call check(nf90_def_dim(nc, "x", nx, d_x), "def x")
      call check(nf90_def_dim(nc, "xf", nx + 1, d_xf), "def xf")
      call check(nf90_def_dim(nc, "yf", ny + 1, d_yf), "def yf")
      call check(nf90_def_var(nc, "time", nf90_double, [d_t], v_t), "def time var")
      write (tunits, "(a,i4.4,a)") "days since ", year, "-01-01 00:00:00"
      call check(nf90_put_att(nc, v_t, "units", trim(tunits)), "att")
      call check(nf90_put_att(nc, v_t, "calendar", "gregorian"), "att")
      call check(nf90_put_att(nc, v_t, "long_name", "centre of the averaging bin"), "att")
      call check(nf90_def_var(nc, "taux", nf90_float, [d_xf, d_y, d_t], v_tx), "def taux")
      call check(nf90_put_att(nc, v_tx, "units", "Pa"), "att")
      call check(nf90_put_att(nc, v_tx, "long_name", &
                              "wind stress along the grid +i axis on u faces (west face of T-cell i)"), "att")
      call check(nf90_def_var(nc, "tauy", nf90_float, [d_x, d_yf, d_t], v_ty), "def tauy")
      call check(nf90_put_att(nc, v_ty, "units", "Pa"), "att")
      call check(nf90_put_att(nc, v_ty, "long_name", &
                              "wind stress along the grid +j axis on v faces (south face of T-cell j)"), "att")
      if (use_ly04) then
         cdtxt = "Large & Yeager (2004) neutral 10 m, U >= 0.5 m/s"
      else
         write (cdtxt, "(a,es10.3)") "constant ", cd_const
      end if
      call check(nf90_put_att(nc, nf90_global, "title", &
                              "JRA55-do 10 m wind -> OM_1deg C-grid bulk wind stress"), "att")
      call check(nf90_put_att(nc, nf90_global, "source", trim(uas_path)//" , "//trim(vas_path)), "att")
      call check(nf90_put_att(nc, nf90_global, "drag_coefficient", trim(cdtxt)), "att")
      call check(nf90_put_att(nc, nf90_global, "rho_air", RHO_AIR), "att")
      call check(nf90_put_att(nc, nf90_global, "method", &
                              "bilinear uas/vas to each face, tau = rho_air*Cd(|U|)*|U|*U per 3-hourly "// &
                              "record (ocean velocity ignored), rotated to the grid axes (angle from node geography), "// &
                              "trapezoidal bin mean"), "att")
      call check(nf90_put_att(nc, nf90_global, "generated_by", "roundabout tools/om1deg_prepare_wind.py"), "att")
      call check(nf90_enddef(nc), "enddef")
      call check(nf90_put_var(nc, v_t, [((b2 - 0.5_dp)*bin_len, b2=1, nbins)]), "put time")
      call check(nf90_put_var(nc, v_tx, real(taux, real32)), "put taux")
      call check(nf90_put_var(nc, v_ty, real(tauy, real32)), "put tauy")
      call check(nf90_close(nc), "close out")
   end subroutine write_output

end program om1deg_wind_regrid
