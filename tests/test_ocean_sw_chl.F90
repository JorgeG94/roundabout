!! Analytical tests for the chlorophyll-dependent N-band penetrating
!! shortwave (`&ocean_thermo_nml sw_scheme = "manizza05"`, Manizza et al.
!! 2005, GRL 32, L05603): `rdb_ocean_sw_chl` (opacity + reader consumer),
!! `sw_band_fractions` / `sw_transmission_nb` / the N-band deposition
!! kernels in `rdb_ocean_surface_flux`, and the KPP `B_0` + EPBL ledger
!! consumers.  Every reference value below is computed from the closed
!! forms in the test itself, never by calling the routine under test.
!!
!! Cases (NetCDF I/O -- run with OMP_NUM_THREADS=1):
!!   1 kappa_values            kappa_n = a_n + b_n chl^p_n at chl = 0.1, 1, 10
!!                             (+ chl floor, + land value) through the fill kernel
!!   2 column_absorption       per-band layer absorption = e^{-k z_top} - e^{-k z_bot},
!!                             total + bed residual = SW_total (conservation)
!!   3 band_split              blue = red = 0.21 SW, NIR = 0.58 SW; NIR split over n>3
!!   4 time_interpolation      2-record cyclic file -> mid-interval linear blend
!!   5 kpp_consistency         KPP B_0 follows the N-band transmission; off == default
!!   6 epbl_consistency        EPBL SW ledger follows the N-band profile; off == default
!!   7 configure_refusals      missing file / variable, conflicting sw_pen_frac
module test_ocean_sw_chl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t, ocean_thermo_config_t, read_config_from_string, &
                         validate_config
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: EOS_VARIANT_LINEAR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_sw_penetration, &
                                     sw_band_fractions, sw_opacity_manizza, &
                                     sw_pe_cost_shape, SEAWATER_CP
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, vmix_apply_kpp_overlay, &
                             KPP_SW_ALL, KPP_SW_MXL
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_compute, EPBL_MSTAR_CONSTANT
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_ocean_sw_chl, only: ocean_sw_chl_t, sw_chl_configure, sw_chl_apply, &
                               sw_chl_fill_kappa, SW_CHL_FLOOR, SW_CHL_YEAR_SECONDS
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            rdb_def_var_1d, nc_enddef, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_sw_chl_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.2_wp
   real(wp), parameter :: BETA_LIN = 0.8_wp
   ! Manizza et al. (2005) coefficients, restated independently of the code.
   real(wp), parameter :: A1 = 0.0232_wp, B1 = 0.074_wp, P1 = 0.674_wp
   real(wp), parameter :: A2 = 0.225_wp, B2 = 0.037_wp, P2 = 0.629_wp
   real(wp), parameter :: K3 = 2.86_wp
   ! Two-band comparison parameters (Jerlov I).
   real(wp), parameter :: R2 = 0.58_wp, ZETA1 = 0.35_wp, ZETA2 = 23.0_wp

contains

   subroutine collect_ocean_sw_chl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kappa_values", test_kappa_values), &
                  new_unittest("column_absorption", test_column_absorption), &
                  new_unittest("band_split", test_band_split), &
                  new_unittest("time_interpolation", test_time_interp), &
                  new_unittest("kpp_consistency", test_kpp_consistency), &
                  new_unittest("epbl_consistency", test_epbl_consistency), &
                  new_unittest("configure_refusals", test_refusals) &
                  ]
   end subroutine collect_ocean_sw_chl_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   pure function kappa_ref(n, chl) result(k)
      !! Independent restatement of the Manizza (2005) opacity.
      integer, intent(in) :: n
      real(wp), intent(in) :: chl
      real(wp) :: k
      select case (n)
      case (1)
         k = A1 + B1*chl**P1
      case (2)
         k = A2 + B2*chl**P2
      case default
         k = K3
      end select
   end function kappa_ref

   pure function trans_ref(d, nb, kap, frac) result(t)
      !! sum_n frac_n exp(-kappa_n d), the N-band transmission.
      real(wp), intent(in) :: d
      integer, intent(in) :: nb
      real(wp), intent(in) :: kap(nb), frac(nb)
      real(wp) :: t
      integer :: n
      t = 0.0_wp
      do n = 1, nb
         t = t + frac(n)*exp(-kap(n)*d)
      end do
   end function trans_ref

   subroutine thermo_cfg(cfg, file, var)
      type(ocean_thermo_config_t), intent(out) :: cfg
      character(len=*), intent(in) :: file, var
      cfg%sw_scheme = "manizza05"
      cfg%sw_chl_file = file
      cfg%sw_chl_var = var
   end subroutine thermo_cfg

   subroutine write_chl_file(filename, nx, ny, nt, t_axis, amp, varname)
      !! `varname(x, y, time)` uniform in space, `amp(k)` at record `k`.
      character(len=*), intent(in) :: filename, varname
      integer, intent(in) :: nx, ny, nt
      real(wp), intent(in) :: t_axis(nt), amp(nt)
      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr, i, j, k
      real(wp) :: f(nx, ny, nt)
      do k = 1, nt
         do j = 1, ny
            do i = 1, nx
               f(i, j, k) = amp(k)
            end do
         end do
      end do
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_3d(ncid, varname, [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_chl_file

   ! -----------------------------------------------------------------
   ! Test 1 -- opacity formulas
   ! -----------------------------------------------------------------

   subroutine test_kappa_values(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sw_chl_t) :: sc
      real(wp), parameter :: CHLS(3) = [0.1_wp, 1.0_wp, 10.0_wp]
      real(wp), allocatable :: wet(:, :), chl(:, :), kap(:, :, :)
      real(wp) :: want, got
      integer :: n, m, j

      call make_grid(grid, 5, 3)
      sc%nb = 3
      sc%a(1:3) = [A1, A2, K3]
      sc%b(1:3) = [B1, B2, 0.0_wp]
      sc%p(1:3) = [P1, P2, 0.0_wp]
      sc%land = 10.0_wp
      allocate (wet(grid%nx_total, grid%ny_total), source=1.0_wp)
      allocate (chl(grid%nx_total, grid%ny_total), source=1.0_wp)
      allocate (kap(grid%nx_total, grid%ny_total, 3), source=-1.0_wp)
      ! Columns 1..3 carry the three test chlorophylls, 4 a negative fill
      ! value (floored), 5 is land.
      do j = 1, grid%ny_total
         do m = 1, 3
            chl(m, j) = CHLS(m)
         end do
         chl(4, j) = -1.0e34_wp
         wet(5, j) = 0.0_wp
      end do

      !$acc enter data copyin(sc, wet, chl, kap)
      call sw_chl_fill_kappa(sc, grid%nx_total, grid%ny_total, wet, chl, kap)
      !$acc update self(kap)
      !$acc exit data delete(sc, wet, chl, kap)

      j = 2
      do m = 1, 3
         do n = 1, 3
            want = kappa_ref(n, CHLS(m))
            got = kap(m, j, n)
            call check(error, abs(got - want) < 1.0e-14_wp, &
                       "kappa must match a + b*chl**p at chl = 0.1, 1, 10")
            if (allocated(error)) return
         end do
      end do
      ! chl = 1: kappa_1 = a1 + b1 exactly, kappa_2 = a2 + b2.
      call check(error, abs(kap(2, j, 1) - (A1 + B1)) < 1.0e-14_wp .and. &
                 abs(kap(2, j, 2) - (A2 + B2)) < 1.0e-14_wp, "chl = 1 gives a + b")
      if (allocated(error)) return
      ! Fill value (negative) -> floored chlorophyll.
      do n = 1, 3
         want = kappa_ref(n, SW_CHL_FLOOR)
         call check(error, abs(kap(4, j, n) - want) < 1.0e-14_wp, &
                    "negative / fill chlorophyll must be floored")
         if (allocated(error)) return
      end do
      ! Land columns carry the land opacity in every band.
      do n = 1, 3
         call check(error, kap(5, j, n) == 10.0_wp, "land opacity must be 10")
         if (allocated(error)) return
      end do
      ! Negative control: opacity grows with chlorophyll (b1, b2 > 0).
      call check(error, kap(3, j, 1) > kap(2, j, 1) .and. kap(2, j, 1) > kap(1, j, 1) .and. &
                 kap(3, j, 2) > kap(1, j, 2), "opacity must increase with chlorophyll")
      if (allocated(error)) return
      call check(error, abs(sw_opacity_manizza(0.1_wp, A1, B1, P1, SW_CHL_FLOOR) - &
                            kappa_ref(1, 0.1_wp)) < 1.0e-14_wp, "pure function agrees")
   end subroutine test_kappa_values

   ! -----------------------------------------------------------------
   ! Test 2 -- column absorption + conservation
   ! -----------------------------------------------------------------

   subroutine deposit_column(grid, sf, nz, dz, absorbed_frac)
      !! Run the SW deposition on a hTr = 0 uniform column and return the
      !! fraction of I0 absorbed in each layer (surface lump added back).
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz
      real(wp), intent(out) :: absorbed_frac(nz)
      type(multilayer_state_t) :: ms
      real(wp), parameter :: DT = 1800.0_wp
      integer :: i, j, k
      real(wp) :: inv_scale, i0

      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = dz
      ms%tracers(ms%idx_temperature)%hTr = 0.0_wp
      i = grid%nx_total/2
      j = grid%ny_total/2
      i0 = sf%Q_heat(i, j)
      inv_scale = DT/(sf%rho0*sf%cp)

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, DT)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr)
         !$acc update self(hT)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)

      do k = 1, nz
         absorbed_frac(k) = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/(inv_scale*i0)
      end do
      absorbed_frac(nz) = absorbed_frac(nz) + 1.0_wp   ! put the surface lump back
      call ms%destroy()
   end subroutine deposit_column

   subroutine test_column_absorption(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 10
      real(wp), parameter :: DZ = 10.0_wp, QSW = 300.0_wp
      real(wp), parameter :: KB(3) = [0.05_wp, 0.02_wp, 0.011_wp]
      real(wp) :: got(NZ), want, ztop, zbot, tot, resid, kap(3), frac(3)
      integer :: n, k

      call make_grid(grid, 6, 4)

      ! ---- (a) each band alone: unit-vector band share ----
      do n = 1, 3
         call sf%init(grid)
         call sf%set_surface_flux_const(QSW, 0.0_wp)
         call sf%set_sw_penetration(0.0_wp, R2, ZETA1, ZETA2)
         call sf%set_sw_chl(3, 0.5_wp)
         sf%sw_frac = 0.0_wp
         sf%sw_frac(n) = 1.0_wp
         do k = 1, 3
            sf%sw_kappa(:, :, k) = KB(k)
         end do
         call deposit_column(grid, sf, NZ, DZ, got)
         do k = 1, NZ
            ztop = real(NZ - k, wp)*DZ
            zbot = ztop + DZ
            if (k > 1) then
               want = exp(-KB(n)*ztop) - exp(-KB(n)*zbot)
            else
               want = exp(-KB(n)*ztop)      ! bed residual lands in the bottom layer
            end if
            call check(error, abs(got(k) - want) < 1.0e-13_wp, &
                       "per-band layer absorption must be e^{-k z_top} - e^{-k z_bot}")
            if (allocated(error)) return
         end do
         call check(error, abs(sum(got) - 1.0_wp) < 1.0e-13_wp, &
                    "absorbed + bed residual must equal SW_total (one band)")
         if (allocated(error)) return
         call sf%destroy()
      end do

      ! ---- (b) chlorophyll bands, default split, conservation ----
      kap = [kappa_ref(1, 0.5_wp), kappa_ref(2, 0.5_wp), kappa_ref(3, 0.5_wp)]
      call sf%init(grid)
      call sf%set_surface_flux_const(QSW, 0.0_wp)
      call sf%set_sw_penetration(0.0_wp, R2, ZETA1, ZETA2)
      call sf%set_sw_chl(3, 0.5_wp)
      frac = [0.21_wp, 0.21_wp, 0.58_wp]
      do k = 1, 3
         sf%sw_kappa(:, :, k) = kap(k)
      end do
      call deposit_column(grid, sf, NZ, DZ, got)
      tot = 0.0_wp
      do k = 1, NZ
         ztop = real(NZ - k, wp)*DZ
         zbot = ztop + DZ
         if (k > 1) then
            want = trans_ref(ztop, 3, kap, frac) - trans_ref(zbot, 3, kap, frac)
         else
            want = trans_ref(ztop, 3, kap, frac)
         end if
         call check(error, abs(got(k) - want) < 1.0e-13_wp, &
                    "3-band layer absorption must equal sum_n f_n (e^{-k z_top} - e^{-k z_bot})")
         if (allocated(error)) return
         tot = tot + got(k)
      end do
      resid = trans_ref(real(NZ - 1, wp)*DZ, 3, kap, frac)   ! what reaches the bed layer top
      call check(error, abs(tot - 1.0_wp) < 1.0e-13_wp, &
                 "total absorbed (incl. bed residual) must equal SW_total")
      if (allocated(error)) return
      call check(error, got(1) > 0.0_wp .and. abs(got(1) - resid) < 1.0e-13_wp, &
                 "bed layer must absorb the whole remaining flux")
      if (allocated(error)) return
      ! Negative control: the two-band profile absorbs differently.
      call check(error, abs(got(NZ) - (R2*(1.0_wp - exp(-DZ/ZETA1)) + &
                                       (1.0_wp - R2)*(1.0_wp - exp(-DZ/ZETA2)))) > 1.0e-3_wp, &
                 "N-band top layer must differ from the two-band top layer")
      call sf%destroy()
   end subroutine test_column_absorption

   ! -----------------------------------------------------------------
   ! Test 3 -- band split
   ! -----------------------------------------------------------------

   subroutine test_band_split(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: f3(3), f5(5), f3b(3)
      real(wp), parameter :: SW = 400.0_wp

      f3 = sw_band_fractions(3, 0.5_wp)
      call check(error, abs(f3(1)*SW - 0.21_wp*SW) < 1.0e-13_wp .and. &
                 abs(f3(2)*SW - 0.21_wp*SW) < 1.0e-13_wp .and. &
                 abs(f3(3)*SW - 0.58_wp*SW) < 1.0e-13_wp, &
                 "blue = red = 0.21 SW, NIR = 0.58 SW")
      if (allocated(error)) return
      call check(error, abs(sum(f3)*SW - SW) < 1.0e-13_wp, "bands must sum to SW_total")
      if (allocated(error)) return

      f3b = sw_band_fractions(3, 0.3_wp)
      call check(error, abs(f3b(1) - 0.126_wp) < 1.0e-15_wp .and. &
                 abs(f3b(2) - 0.294_wp) < 1.0e-15_wp .and. abs(f3b(3) - 0.58_wp) < 1.0e-15_wp, &
                 "blue_frac = 0.3 gives 0.126 / 0.294 / 0.58")
      if (allocated(error)) return

      f5 = sw_band_fractions(5, 0.5_wp)
      call check(error, abs(f5(3) - 0.58_wp/3.0_wp) < 1.0e-15_wp .and. &
                 abs(f5(4) - f5(3)) < 1.0e-15_wp .and. abs(f5(5) - f5(3)) < 1.0e-15_wp, &
                 "NIR must be split evenly over bands 3..n")
      if (allocated(error)) return
      call check(error, abs(sum(f5) - 1.0_wp) < 1.0e-14_wp, "5 bands sum to 1")
   end subroutine test_band_split

   ! -----------------------------------------------------------------
   ! Test 4 -- cyclic time interpolation through the shared reader
   ! -----------------------------------------------------------------

   subroutine test_time_interp(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 3, NT = 2
      character(len=*), parameter :: FNAME = "test_sw_chl_interp.nc"
      real(wp), parameter :: C0 = 0.2_wp, C1 = 1.4_wp
      type(ocean_thermo_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_sw_chl_t) :: sc
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      real(wp), allocatable :: wet(:, :)
      real(wp) :: t_axis(NT), tq(3), mid, want
      integer :: ierr, iq, i, j, n

      t_axis = [0.0_wp, 0.5_wp*SW_CHL_YEAR_SECONDS]
      call write_chl_file(FNAME, NXP, NYP, NT, t_axis, [C0, C1], "CHL")

      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      call sf%init(grid)
      call reader%init()
      call thermo_cfg(cfg, FNAME, "CHL")
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. sc%active, "configure must succeed")
      if (allocated(error)) return
      call check(error, sf%sw_chl .and. sf%sw_nb == 3 .and. sf%has_sw .and. sf%sw_pen_frac == 1.0_wp, &
                 "configure must switch the flux slot to the N-band path with I0 = 100% of SW")
      if (allocated(error)) return

      allocate (wet(grid%nx_total, grid%ny_total), source=1.0_wp)
      wet(NGHOST + 1, NGHOST + 1) = 0.0_wp    ! one land column

      mid = 0.5_wp*(C0 + C1)
      ! Mid-interval of records 1 and 2, the mid-interval of the wrap
      ! segment (record 2 -> record 1 of the next year), and the same
      ! point one full year later.
      tq = [0.25_wp*SW_CHL_YEAR_SECONDS, 0.75_wp*SW_CHL_YEAR_SECONDS, &
            1.25_wp*SW_CHL_YEAR_SECONDS]

      !$acc enter data copyin(reader)
      call reader%enter_data()
      call sf%enter_data()
      !$acc enter data copyin(wet)
      do iq = 1, 3
         call ocean_data_input_update_all(reader, tq(iq))
         call sw_chl_apply(sc, reader, grid, wet, sf, tq(iq))
         !$acc update self(sf%sw_chl_field, sf%sw_kappa)
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               if (wet(i, j) < 0.5_wp) then
                  call check(error, sf%sw_kappa(i, j, 1) == 10.0_wp, "land column keeps land opacity")
                  if (allocated(error)) exit
                  cycle
               end if
               call check(error, abs(sf%sw_chl_field(i, j) - mid) < 1.0e-13_wp, &
                          "chlorophyll must be the mid-interval linear blend")
               if (allocated(error)) exit
               do n = 1, 3
                  want = kappa_ref(n, mid)
                  call check(error, abs(sf%sw_kappa(i, j, n) - want) < 1.0e-13_wp, &
                             "kappa must follow the blended chlorophyll")
                  if (allocated(error)) exit
               end do
               if (allocated(error)) exit
            end do
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      !$acc exit data delete(wet)
      call sf%exit_data()
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()
      call sf%destroy()
   end subroutine test_time_interp

   ! -----------------------------------------------------------------
   ! Test 5 -- KPP B_0 consumes the same N-band transmission
   ! -----------------------------------------------------------------

   subroutine kpp_case(grid, ss, sf, method, vmix)
      !! Unstratified 8 x 5 m column (no wind): h_b = full depth,
      !! B_0-independent (c_vt2 = 0), so only the SW term moves kv.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: method
      type(ocean_vmix_t), intent(inout) :: vmix
      type(multilayer_state_t) :: ms
      integer, parameter :: NZ = 8
      integer :: k
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 5.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%rho_layer(:, :, k) = 1030.0_wp
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 15.0_wp*5.0_wp
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*5.0_wp
      end do
      call vmix%init(grid, nz_ml=NZ)
      vmix%c_vt2 = 0.0_wp
      vmix%kpp_sw_method = method
      vmix%rho0 = RHO0
      vmix%eos%variant = EOS_VARIANT_LINEAR
      vmix%eos%alpha_T = ALPHA_LIN
      vmix%eos%beta_S = BETA_LIN
      vmix%eos%rho0 = RHO0
      vmix%pp81_nu0 = 0.0_wp

      !$acc enter data copyin(ms, vmix, ss, sf)
      call ms%enter_data()
      call vmix%enter_data()
      call sf%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth)
      call sf%exit_data()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
      call ms%destroy()
   end subroutine kpp_case

   subroutine seed_sf(grid, sf, qsw, qcool)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: qsw, qcool
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_surface_flux_const(qcool, 0.0_wp)
      call sf%set_sw_penetration(1.0_wp, R2, ZETA1, ZETA2, sw_source="q_sw")
      sf%q_sw = qsw
   end subroutine seed_sf

   subroutine test_kpp_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf_def, sf_chl, sf_off
      type(ocean_vmix_t) :: v_def, v_chl_mxl, v_chl_all, v_off
      real(wp), parameter :: QSW = 400.0_wp, QCOOL = -50.0_wp
      real(wp), parameter :: KAP(3) = [0.1_wp, 0.3_wp, 2.86_wp]
      real(wp), parameter :: FR(3) = [0.21_wp, 0.21_wp, 0.58_wp]
      real(wp) :: h_b, q_bl, b0_all, b0_mxl, ratio_pred, ratio_obs, kv_m, kv_a, kv_def_m
      integer :: i, j, k, kmax, n

      call make_grid(grid, 8, 6)
      call ss%init(grid)
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)    ! u_star = 0 => w_s = w_*
      i = grid%nx_total/2
      j = grid%ny_total/2

      call seed_sf(grid, sf_def, QSW, QCOOL)
      call seed_sf(grid, sf_chl, QSW, QCOOL)
      call sf_chl%set_sw_chl(3, 0.5_wp)
      do n = 1, 3
         sf_chl%sw_kappa(:, :, n) = KAP(n)
      end do
      call seed_sf(grid, sf_off, QSW, QCOOL)
      call sf_off%set_sw_chl(3, 0.5_wp)
      do n = 1, 3
         sf_off%sw_kappa(:, :, n) = KAP(n)
      end do
      sf_off%sw_chl = .false.        ! scheme off, arrays populated

      call kpp_case(grid, ss, sf_def, KPP_SW_MXL, v_def)
      call kpp_case(grid, ss, sf_chl, KPP_SW_MXL, v_chl_mxl)
      call kpp_case(grid, ss, sf_chl, KPP_SW_ALL, v_chl_all)
      call kpp_case(grid, ss, sf_off, KPP_SW_MXL, v_off)

      checks: block
         ! (a) scheme off == default, bit for bit.
         call check(error, all(v_off%kv == v_def%kv) .and. all(v_off%kt == v_def%kt), &
                    "sw_chl = .false. must leave KPP identical to the two-band default")
         if (allocated(error)) exit checks

         ! (b) with the scheme on, B_0 uses the N-band T(h_b): the kv ratio
         !     mxl/all is the cube root of the B_0 ratio.
         h_b = v_chl_mxl%bl_depth(i, j)
         call check(error, abs(h_b - v_chl_all%bl_depth(i, j)) < 1.0e-12_wp*h_b, &
                    "BL depth must be B_0-independent (c_vt2 = 0)")
         if (allocated(error)) exit checks
         q_bl = QCOOL - QSW*trans_ref(h_b, 3, KAP, FR)
         b0_all = GRAVITY*ALPHA_LIN*QCOOL/(RHO0*SEAWATER_CP)
         b0_mxl = GRAVITY*ALPHA_LIN*q_bl/(RHO0*SEAWATER_CP)
         call check(error, b0_mxl < b0_all .and. b0_all < 0.0_wp, "expected B_0(mxl) < B_0(all) < 0")
         if (allocated(error)) exit checks
         ratio_pred = ((-b0_mxl)/(-b0_all))**(1.0_wp/3.0_wp)
         kmax = 0
         kv_m = 0.0_wp
         do k = 2, 8
            if (v_chl_mxl%kv(i, j, k) > kv_m) then
               kv_m = v_chl_mxl%kv(i, j, k)
               kmax = k
            end if
         end do
         call check(error, kmax > 0, "expected a KPP kv bump")
         if (allocated(error)) exit checks
         kv_a = v_chl_all%kv(i, j, kmax)
         ratio_obs = kv_m/kv_a
         call check(error, abs(ratio_obs - ratio_pred) < 1.0e-9_wp*ratio_pred, &
                    "kv(mxl)/kv(all) must follow the N-band transmission at h_b")
         if (allocated(error)) exit checks

         ! (c) it differs from the two-band default (a column where it matters).
         kv_def_m = v_def%kv(i, j, kmax)
         call check(error, abs(kv_m - kv_def_m) > 1.0e-3_wp*kv_def_m, &
                    "N-band KPP kv must differ from the two-band default")
      end block checks
      call v_def%destroy(); call v_chl_mxl%destroy(); call v_chl_all%destroy(); call v_off%destroy()
      call sf_def%destroy(); call sf_chl%destroy(); call sf_off%destroy(); call ss%destroy()
   end subroutine test_kpp_consistency

   ! -----------------------------------------------------------------
   ! Test 6 -- EPBL SW TKE ledger consumes the N-band profile
   ! -----------------------------------------------------------------

   subroutine epbl_case(grid, ss, sf, nz, dz, epbl)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz
      type(ocean_epbl_t), intent(inout) :: epbl
      type(multilayer_state_t) :: ms
      integer :: k
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 12.0_wp*dz
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*dz
      end do
      call epbl%init(grid, nz_ml=nz)
      epbl%enable = .true.
      epbl%mstar_scheme = EPBL_MSTAR_CONSTANT
      epbl%mstar_const = 1.2_wp
      epbl%nstar = 0.2_wp
      epbl%eos%variant = EOS_VARIANT_LINEAR
      epbl%eos%alpha_T = ALPHA_LIN
      epbl%eos%beta_S = BETA_LIN
      epbl%eos%rho0 = RHO0
      epbl%rho0 = RHO0
      epbl%f_centre = 0.0_wp
      epbl%tke_diags = .true.
      epbl%epbl_sw_ctke = .true.

      !$acc enter data copyin(ms, epbl, ss, sf)
      call ms%enter_data()
      call sf%enter_data()
      call epbl%enter_data()
      call epbl_compute(grid, epbl, ms, ss, 1800.0_wp, sf=sf)
      !$acc update self(epbl%ctke_sw%data, epbl%dcolht_t%data)
      call epbl%exit_data()
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, sf)
      call ms%destroy()
   end subroutine epbl_case

   subroutine test_epbl_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf_def, sf_chl, sf_off
      type(ocean_epbl_t) :: e_def, e_chl, e_off
      integer, parameter :: NZ = 20
      real(wp), parameter :: DZ = 5.0_wp, DT = 1800.0_wp, QSW = 300.0_wp
      real(wp), parameter :: KAP(3) = [0.08_wp, 0.25_wp, 2.86_wp]
      real(wp), parameter :: FR(3) = [0.21_wp, 0.21_wp, 0.58_wp]
      real(wp) :: want(NZ), prefac, ztop, zbot, heat, scale, dcol, diff_def
      integer :: i, j, k, n

      call make_grid(grid, 6, 4)
      call ss%init(grid)
      call ss%set_wind_stress_const(0.05_wp, 0.0_wp)
      i = grid%nx_total/2
      j = grid%ny_total/2

      call sf_def%init(grid)
      call sf_def%set_surface_flux_const(QSW, 0.0_wp)
      call sf_def%set_sw_penetration(1.0_wp, R2, ZETA1, ZETA2)
      call sf_chl%init(grid)
      call sf_chl%set_surface_flux_const(QSW, 0.0_wp)
      call sf_chl%set_sw_penetration(0.0_wp, R2, ZETA1, ZETA2)
      call sf_chl%set_sw_chl(3, 0.5_wp)
      call sf_off%init(grid)
      call sf_off%set_surface_flux_const(QSW, 0.0_wp)
      call sf_off%set_sw_penetration(1.0_wp, R2, ZETA1, ZETA2)
      call sf_off%set_sw_chl(3, 0.5_wp)
      do n = 1, 3
         sf_chl%sw_kappa(:, :, n) = KAP(n)
         sf_off%sw_kappa(:, :, n) = KAP(n)
      end do
      sf_off%sw_chl = .false.

      call epbl_case(grid, ss, sf_def, NZ, DZ, e_def)
      call epbl_case(grid, ss, sf_chl, NZ, DZ, e_chl)
      call epbl_case(grid, ss, sf_off, NZ, DZ, e_off)

      checks: block
         call check(error, all(e_off%ctke_sw%data == e_def%ctke_sw%data), &
                    "sw_chl = .false. must leave the EPBL SW ledger identical to the default")
         if (allocated(error)) exit checks

         ! Analytical ledger: uniform linear-EOS column => constant prefactor.
         dcol = DZ*ALPHA_LIN/RHO0
         prefac = -0.5_wp*GRAVITY*RHO0*dcol
         scale = 0.0_wp
         do k = 1, NZ
            ztop = real(NZ - k, wp)*DZ
            zbot = ztop + DZ
            want(k) = 0.0_wp
            do n = 1, 3
               if (k > 1) then
                  heat = QSW*FR(n)*DT/(RHO0*SEAWATER_CP)*(exp(-KAP(n)*ztop) - exp(-KAP(n)*zbot))
               else
                  heat = QSW*FR(n)*DT/(RHO0*SEAWATER_CP)*exp(-KAP(n)*ztop)
               end if
               want(k) = want(k) + prefac*heat*sw_pe_cost_shape(DZ*KAP(n))
            end do
            scale = max(scale, abs(want(k)))
         end do
         do k = 1, NZ
            call check(error, abs(e_chl%ctke_sw%data(i, j, k) - want(k)) < 1.0e-9_wp*scale, &
                       "EPBL ledger must equal sum_n heat_n Phi(h kappa_n) layer by layer")
            if (allocated(error)) exit checks
         end do

         ! The default two-band ledger is measurably different.
         diff_def = maxval(abs(e_chl%ctke_sw%data(i, j, :) - e_def%ctke_sw%data(i, j, :)))
         call check(error, diff_def > 1.0e-3_wp*scale, &
                    "N-band ledger must differ from the two-band ledger")
      end block checks
      call e_def%destroy(); call e_chl%destroy(); call e_off%destroy()
      call sf_def%destroy(); call sf_chl%destroy(); call sf_off%destroy(); call ss%destroy()
   end subroutine test_epbl_consistency

   ! -----------------------------------------------------------------
   ! Test 7 -- configure-time refusals
   ! -----------------------------------------------------------------

   function validate_thermo(body) result(ierr)
      character(len=*), intent(in) :: body
      integer :: ierr
      type(config_t) :: cfg
      ierr = -999
      call read_config_from_string( &
         "&sim_nml sim_type = 'ocean' /"//new_line('a')// &
         "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"// &
         new_line('a')//"&ocean_thermo_nml "//body//" /"//new_line('a'), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr=ierr)
   end function validate_thermo

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 3
      character(len=*), parameter :: FNAME = "test_sw_chl_refuse.nc"
      type(ocean_thermo_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_sw_chl_t) :: sc
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      integer :: ierr

      ! ---- validate_config: namelist level ----
      call check(error, validate_thermo("sw_scheme = 'manizza05', sw_chl_file = 'x.nc', "// &
                                        "sw_chl_var = 'chl_a'") == OCEAN_STATUS_OK, &
                 "a complete manizza05 block must validate (and its keys must parse)")
      if (allocated(error)) return
      call check(error, validate_thermo("sw_scheme = 'manizza05', sw_chl_file = 'x.nc', "// &
                                        "sw_pen_frac = 0.5") /= OCEAN_STATUS_OK, &
                 "a conflicting sw_pen_frac must be refused")
      if (allocated(error)) return
      call check(error, validate_thermo("sw_scheme = 'manizza05'") /= OCEAN_STATUS_OK, &
                 "manizza05 without sw_chl_file must be refused")
      if (allocated(error)) return
      call check(error, validate_thermo("sw_scheme = 'bogus'") /= OCEAN_STATUS_OK, &
                 "an unknown sw_scheme must be refused")
      if (allocated(error)) return
      call check(error, validate_thermo("sw_scheme = 'default', sw_pen_frac = 0.5") == OCEAN_STATUS_OK, &
                 "the default scheme keeps its sw_pen_frac semantics")
      if (allocated(error)) return

      ! ---- sw_chl_configure ----
      call write_chl_file(FNAME, NXP, NYP, 2, [0.0_wp, 1.0e6_wp], [0.2_wp, 0.4_wp], "CHL")
      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      call sf%init(grid)
      call reader%init()

      call thermo_cfg(cfg, "no_such_chl_file.nc", "CHL")
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK .and. .not. sc%active .and. .not. sf%sw_chl, &
                 "a missing file must be refused and leave the flux slot untouched")
      if (allocated(error)) return

      call thermo_cfg(cfg, FNAME, "not_a_variable")
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK .and. .not. sf%sw_chl, &
                 "a missing variable must be refused")
      if (allocated(error)) return

      call thermo_cfg(cfg, "", "CHL")
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, "an empty sw_chl_file must be refused")
      if (allocated(error)) return

      call thermo_cfg(cfg, FNAME, "CHL")
      cfg%sw_pen_frac = 0.5_wp
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK .and. .not. sf%sw_chl, &
                 "a conflicting sw_pen_frac must be refused at configure")
      if (allocated(error)) return

      call thermo_cfg(cfg, FNAME, "CHL")
      cfg%sw_nbands = 4                       ! band 4 has no coefficients
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, "a band without an opacity must be refused")
      if (allocated(error)) return

      ! Positive control: the same file + variable configures cleanly.
      call thermo_cfg(cfg, FNAME, "CHL")
      call sw_chl_configure(cfg, reader, grid, sf, sc, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. sc%active .and. sf%sw_chl, &
                 "a valid configuration must configure")
      call reader%destroy()
      call sf%destroy()
   end subroutine test_refusals

end module test_ocean_sw_chl
