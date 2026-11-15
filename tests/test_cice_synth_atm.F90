!! The synthetic prescribed atmosphere of the CICE-port forcing and closure tests (test-only).
!! Nine fields, each affine in (lon, lat, time) so the bilinear regrid and the time blend reproduce the
!! closed form exactly; written as one small NetCDF file plus the file list the forcing reader takes.
module test_cice_synth_atm
   use rdb_constants, only: wp
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            nc_def_var_2d, rdb_def_var_1d, nc_enddef, rdb_put_var_1d, nc_put_att
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: write_synth_atm, remove_synth_atm

   integer, parameter :: NLON_SRC = 36, NLAT_SRC = 18

contains

   subroutine uniform_axis(n, lo, hi, c, bnds)
      integer, intent(in) :: n
      real(wp), intent(in) :: lo, hi
      real(wp), intent(out) :: c(n), bnds(2, n)
      real(wp) :: dx
      integer :: k

      dx = (hi - lo)/real(n, wp)
      do k = 1, n
         bnds(1, k) = lo + real(k - 1, wp)*dx
         bnds(2, k) = lo + real(k, wp)*dx
         c(k) = 0.5_wp*(bnds(1, k) + bnds(2, k))
      end do
   end subroutine uniform_axis

   subroutine atm_coefs(scn, coefs)
      !! `f = c1 + c2*lon + c3*lat + c4*t[days]` of the nine fields (tas, huss, uas, vas, psl, rlds, rsds,
      !! prra, prsn).  `warm`: air above freezing, no ice forms.  `cold`: ice persists and grows.
      character(len=*), intent(in) :: scn
      real(wp), intent(out) :: coefs(4, 9)

      coefs(:, 1) = [255.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]               ! tas
      coefs(:, 2) = [4.0e-4_wp, 0.0_wp, 0.0_wp, 0.0_wp]              ! huss
      coefs(:, 3) = [3.0_wp, 0.8_wp, 0.5_wp, 30.0_wp]                ! uas
      coefs(:, 4) = [-1.0_wp, 0.3_wp, 0.4_wp, -20.0_wp]              ! vas
      coefs(:, 5) = [101000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]            ! psl
      coefs(:, 6) = [170.0_wp, 0.5_wp, 0.0_wp, 20.0_wp]              ! rlds
      coefs(:, 7) = [80.0_wp, 0.0_wp, 0.0_wp, 120.0_wp]              ! rsds
      coefs(:, 8) = [6.0e-5_wp, 0.0_wp, 0.0_wp, 4.0e-5_wp]           ! prra
      coefs(:, 9) = [3.0e-5_wp, 0.0_wp, 0.0_wp, 0.0_wp]              ! prsn
      select case (scn)
      case ("warm")
         coefs(1, 1) = 285.0_wp
         coefs(1, 7) = 250.0_wp
      case ("cold")
         continue
      case ("mixed_dry")
         coefs(1, 1) = 250.0_wp - 14.0_wp*70.0_wp
         coefs(3, 1) = 14.0_wp
         coefs(1, 7) = 150.0_wp
         coefs(:, 8) = 0.0_wp
         coefs(:, 9) = 0.0_wp
      case ("mixed")
         ! a latitude gradient of the air temperature: freezing in the south, melting in the north of the basin
         coefs(1, 1) = 250.0_wp - 14.0_wp*70.0_wp
         coefs(3, 1) = 14.0_wp
         coefs(1, 7) = 150.0_wp
      case default
         error stop "test_cice_forcing: unknown atmosphere scenario"
      end select
   end subroutine atm_coefs

   subroutine write_synth_atm(scn, fn_atm, fn_lst)
      !! Write the atmosphere `fn_atm` and its file list `fn_lst` for scenario `scn` (`warm`, `cold`, `mixed`).
      character(len=*), intent(in) :: scn, fn_atm, fn_lst
      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", "rlds", "rsds", "prra", "prsn"]
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC), lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(3), coefs(4, 9)
      real(wp), allocatable :: f(:, :, :)
      integer :: ncid, dim_lon, dim_lat, dim_time, dim_bnds
      integer :: vid_lon, vid_lat, vid_lonb, vid_latb, vid_time, vid_f(9)
      integer :: v, i, j, k, ierr, u

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 1.0_wp, 2.0_wp]
      call atm_coefs(scn, coefs)
      allocate (f(NLON_SRC, NLAT_SRC, 3))
      call nc_create_file(fn_atm, ncid)
      call nc_def_dim(ncid, "lon", NLON_SRC, dim_lon)
      call nc_def_dim(ncid, "lat", NLAT_SRC, dim_lat)
      call nc_def_dim(ncid, "bnds", 2, dim_bnds)
      call nc_def_dim(ncid, "time", 3, dim_time)
      call rdb_def_var_1d(ncid, "lon", dim_lon, vid_lon)
      call rdb_def_var_1d(ncid, "lat", dim_lat, vid_lat)
      call nc_def_var_2d(ncid, "lon_bnds", [dim_bnds, dim_lon], vid_lonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dim_bnds, dim_lat], vid_latb)
      call rdb_def_var_1d(ncid, "time", dim_time, vid_time)
      call nc_put_att(ncid, vid_time, "units", "days since 2000-01-01 00:00:00")
      do v = 1, 9
         call nc_def_var_3d(ncid, trim(VARNAMES(v)), [dim_lon, dim_lat, dim_time], vid_f(v))
      end do
      call nc_enddef(ncid)
      call rdb_put_var_1d(ncid, vid_lon, lon_c)
      call rdb_put_var_1d(ncid, vid_lat, lat_c)
      ierr = nf90_put_var(ncid, vid_lonb, lon_bnds)
      ierr = nf90_put_var(ncid, vid_latb, lat_bnds)
      call rdb_put_var_1d(ncid, vid_time, t_axis)
      do v = 1, 9
         do k = 1, 3
            do j = 1, NLAT_SRC
               do i = 1, NLON_SRC
                  f(i, j, k) = coefs(1, v) + coefs(2, v)*lon_c(i) + coefs(3, v)*lat_c(j) + coefs(4, v)*t_axis(k)
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vid_f(v), f)
      end do
      call nc_close(ncid)
      open (newunit=u, file=fn_lst, status="replace", action="write")
      write (u, '(a)') fn_atm
      close (u)
   end subroutine write_synth_atm

   subroutine remove_synth_atm(fn_atm, fn_lst)
      character(len=*), intent(in) :: fn_atm, fn_lst
      integer :: u, ios

      open (newunit=u, file=fn_atm, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
      open (newunit=u, file=fn_lst, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine remove_synth_atm

end module test_cice_synth_atm
