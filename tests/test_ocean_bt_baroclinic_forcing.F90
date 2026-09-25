!! THE BAROCLINIC FORCING OF THE BAROTROPIC MODE: a horizontal density
!! gradient at rest, with a flat free surface, accelerates the depth mean
!! at `-(g H / 2 rho_0) d(rho)/dx` — and the split solver must deliver it.
module test_ocean_bt_baroclinic_forcing
   !! ### The defect this gates
   !!
   !! The split solver built the barotropic substep's frozen forcing as
   !! `F_bt_fast = F_bt - <PGF>_h`: the WHOLE thickness-weighted depth mean
   !! of the layer pressure force was taken out, on the reading that the
   !! substep's own `-g*grad(eta)` replaces it.  That is only true for a
   !! barotropic pressure field.  The layer PGF of a sloping density field
   !! has a depth mean of its own — the bottom-pressure gradient, whose
   !! curl over topography is JEBAR — and `apply_bt_correction` subtracts
   !! `dt*F_bt` from every layer, so it was removed from the layers too:
   !! the depth mean ended every stage at `u_bt^end`, which never felt it.
   !! Every baroclinic test still passed (the layer SHEAR was right) and
   !! every flat-bottom rest test still passed (no density gradient, no
   !! depth-mean PGF).  On the global 1-degree WOA13 spin-up it held Drake
   !! Passage at 0.2 Sv on day 1 where MOM6, forced by the same density
   !! field alone, surges to -188 Sv and settles near 160 Sv.
   !!
   !! MOM6 keeps the full `PFu` in `BT_force` (`MOM_barotropic.F90`,
   !! `BT_force_u += wt_u*bc_accel_u`, `bc_accel = CAu + PFu + diffu` from
   !! `MOM_dynamics_split_RK2.F90`) and lets the barotropic loop see only
   !! the ANOMALY `-gtot*grad(eta - eta_PF)` about the free surface the
   !! slow PGF was built on (`btloop_find_PF`).  `&ocean_bt_nml
   !! bc_pgf_forcing` (default on) is that split: the forcing sheds only
   !! the free-surface term the slow PGF itself carries (none for the
   !! surface-relative Montgomery / FV-lite forms, `g*grad(eta_PF)` for
   !! FV-MOM6), see `set_fast_forcing_eta_pf`.
   !!
   !! ### The case, and the analytic answer
   !!
   !! A flat-bottomed (H = 1000 m) re-entrant channel, periodic in x,
   !! walled in y, f = 0, no wind, no drag, no viscosity, at rest with
   !! `eta = 0`, sigma layers, and a temperature
   !! `T = T_ref + dT*cos(2*pi*x/L)` that is UNIFORM IN DEPTH — so the
   !! density anomaly `rho'(x)` is the same in every layer and there is no
   !! vertical stratification for anything to adjust against.  The
   !! hydrostatic pressure anomaly is `p'(z) = -g*rho'*z` and its depth
   !! mean gradient
   !!
   !!     dU/dt = -(1/rho_0) (1/H) int_{-H}^{0} d(p')/dx dz
   !!           = -(g H / (2 rho_0)) d(rho')/dx
   !!
   !! is the barotropic acceleration at t = 0, evaluated here on the
   !! model's own discrete `rho_layer` (face difference over `dx`).  After
   !! ONE outer step of `dt` the depth-mean face velocity is `a*dt` up to
   !! the free surface the acceleration itself raises: `eta ~ -H*da/dx*dt^2/2`
   !! pushes back with `g*grad(eta)*dt/2`, a relative `(c*k*dt)^2/4` with
   !! `c = sqrt(g H)` and `k = 2*pi/L` (0.14 % here).  The asserted
   !! tolerance is `(c*k*dt)^2 + 1e-3`; the legacy split gives a ratio of
   !! ZERO to round-off, which is the fails-before half, kept as a control.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_rho_layer_ptr, rdb_ocean_get_grid_info, &
                            rdb_ocean_set_tracer
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_bt_baroclinic_forcing_tests

   integer, parameter :: NX_PHYS = 32
   real(wp), parameter :: DXY = 31250.0_wp
      !! Channel length `L = NX_PHYS*DXY = 1000 km`: one cosine wavelength.
   real(wp), parameter :: H_BED = 1000.0_wp
   real(wp), parameter :: RHO_0 = 1035.0_wp
      !! `&ocean_ic_nml rho_0` (its default, set explicitly below).
   real(wp), parameter :: T_REF = 10.0_wp
   real(wp), parameter :: DT_T = 2.0_wp
      !! Temperature amplitude (degC).  With `alpha_T = 0.2` kg/m^3/degC the
      !! density anomaly is +-0.4 kg/m^3 and the peak acceleration
      !! `(g H/2 rho_0)*0.4*k` is 1.2e-5 m/s^2.
   real(wp), parameter :: DT = 120.0_wp
   real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)

contains

   subroutine collect_ocean_bt_baroclinic_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bc_pgf_accelerates_depth_mean_pred_corr_mont", test_pc_mont), &
                  new_unittest("bc_pgf_accelerates_depth_mean_ssp_rk2_mont", test_ssp_mont), &
                  new_unittest("bc_pgf_accelerates_depth_mean_pred_corr_fv_mom6", test_pc_fv_mom6), &
                  new_unittest("bc_pgf_accelerates_depth_mean_ssp_rk2_fv_mom6", test_ssp_fv_mom6), &
                  new_unittest("legacy_split_discards_depth_mean_bc_pgf", test_legacy_control) &
                  ]
   end subroutine collect_ocean_bt_baroclinic_forcing_tests

   function case_nml(scheme, pgf_form, bc_pgf_forcing) result(nml)
      character(len=*), intent(in) :: scheme, pgf_form
      logical, intent(in) :: bc_pgf_forcing
      character(len=:), allocatable :: nml
      character(len=8) :: flag

      flag = ".true."
      if (.not. bc_pgf_forcing) flag = ".false."
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 32, ny = 4, nghost = 3, dx = 31250.0, dy = 31250.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 10 /"//new_line("a")// &
            "&time_nml t_end = 1.0e7, dt_fixed = 120.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0, taux_magnitude = 0.0 /"//new_line("a")// &
            "&physics_nml coriolis_f = 0.0 /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, T_init_bottom = 10.0, "// &
            "T_init_surface = 10.0 /"//new_line("a")// &
            "&ocean_ic_nml alpha_T = 0.2, T_ref = 10.0, S_ref = 35.0, rho_0 = 1035.0 /"// &
            new_line("a")// &
            "&ocean_pgf_nml form = '"//pgf_form//"' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic' /"//new_line("a")// &
            "&ocean_bdrag_nml form = 'linear', r = 0.0 /"//new_line("a")// &
            "&ocean_bt_nml split_scheme = '"//scheme//"', auto_n_inner = .false., "// &
            "n_inner = 10, bc_pgf_forcing = "//trim(flag)//" /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function case_nml

   subroutine run_case(scheme, pgf_form, bc_pgf_forcing, ratio_min, ratio_max, amax, ok)
      !! Create the channel, impose the depth-uniform `T(x)`, read the
      !! model's own `rho_layer` for the analytic acceleration, take ONE
      !! outer step, and return the range of `U/(a*dt)` over the faces
      !! where `|a|` is at least half its peak (away from the sine's zeros,
      !! where the ratio is 0/0).
      character(len=*), intent(in) :: scheme, pgf_form
      logical, intent(in) :: bc_pgf_forcing
      real(wp), intent(out) :: ratio_min, ratio_max, amax
      logical, intent(out) :: ok

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen, nxp, nyp, nzp, ngc
      real(wp), pointer :: r3(:, :, :), h3(:, :, :), u3(:, :, :)
      real(wp), allocatable :: tval(:, :, :), accel(:, :)
      character(len=:), allocatable :: nml
      integer :: i, j, k, ng
      real(wp) :: num, den, hf, ubar, r

      ok = .false.
      ratio_min = huge(1.0_wp)
      ratio_max = -huge(1.0_wp)
      amax = 0.0_wp
      handle = c_null_ptr
      nml = case_nml(scheme, pgf_form, bc_pgf_forcing)
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return
      status = rdb_ocean_get_grid_info(handle, nxp, nyp, nzp, ngc)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      ng = int(ngc)

      allocate (tval(int(nxp), int(nyp), int(nzp)))
      do k = 1, int(nzp)
         do j = 1, int(nyp)
            do i = 1, int(nxp)
               tval(i, j, k) = T_REF + DT_T*cos(TWO_PI*(real(i, wp) - 0.5_wp)/real(nxp, wp))
            end do
         end do
      end do
      status = rdb_ocean_set_tracer(handle, "temperature", 11_c_int, tval, nxp, nyp, nzp)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if

      ! The analytic acceleration on the model's own discrete density
      ! (`rdb_ocean_set_tracer` recomputes `rho_layer` and pulls it back
      ! to the host).  Layer 1 stands for every layer: T is depth-uniform.
      status = rdb_ocean_get_rho_layer_ptr(handle, ptr, nx, ny, nz, gen)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      call c_f_pointer(ptr, r3, [int(nx), int(ny), int(nz)])
      allocate (accel(int(nxp) + 1, int(nyp)))
      accel = 0.0_wp
      do j = 1, int(nyp)
         do i = 2, int(nxp)
            accel(i, j) = -GRAVITY*H_BED/(2.0_wp*RHO_0)* &
                          (r3(ng + i, ng + j, 1) - r3(ng + i - 1, ng + j, 1))/DXY
         end do
      end do
      amax = maxval(abs(accel))

      status = rdb_ocean_step(handle, 1_c_int)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_refresh_host(handle)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, h3, [int(nx), int(ny), int(nz)])
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, u3, [int(nx), int(ny), int(nz)])
      ok = all(ieee_is_finite(u3)) .and. all(ieee_is_finite(h3))

      ! Face `i` (physical) sits between cells `i-1` and `i`; the array
      ! face index is `ng + i`.  Interior faces only (the seam faces are
      ! periodic copies of these).
      do j = 1, int(nyp)
         do i = 2, int(nxp)
            if (abs(accel(i, j)) < 0.5_wp*amax) cycle
            num = 0.0_wp
            den = 0.0_wp
            do k = 1, int(nzp)
               hf = 0.5_wp*(h3(ng + i - 1, ng + j, k) + h3(ng + i, ng + j, k))
               num = num + hf*u3(ng + i, ng + j, k)
               den = den + hf
            end do
            ubar = num/den
            r = ubar/(accel(i, j)*DT)
            ratio_min = min(ratio_min, r)
            ratio_max = max(ratio_max, r)
         end do
      end do

      status = rdb_ocean_destroy(handle)
      ok = ok .and. (status == OCEAN_STATUS_OK)
   end subroutine run_case

   pure function ratio_tol() result(tol)
      !! `(c*k*dt)^2 + 1e-3`: the free surface the acceleration raises
      !! within the step pushes back by a relative `(c*k*dt)^2/4`; the
      !! factor 4 and the 1e-3 cover the scheme's own O(dt^2) stage
      !! weighting.  5.6e-3 + 1e-3 = 6.6e-3 here; measured |1 - ratio| is
      !! 9.4e-4 (pred_corr) and 3.7e-3 (ssp_rk2), both PGF forms.
      real(wp) :: tol
      tol = (sqrt(GRAVITY*H_BED)*(TWO_PI/(real(NX_PHYS, wp)*DXY))*DT)**2 + 1.0e-3_wp
   end function ratio_tol

   subroutine check_accelerates(error, scheme, pgf_form)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme, pgf_form
      real(wp) :: rmin, rmax, amax
      logical :: ok
      character(len=256) :: msg

      call run_case(scheme, pgf_form, .true., rmin, rmax, amax, ok)
      call check(error, ok, "the channel must run one step and stay finite")
      if (allocated(error)) return
      call check(error, amax > 1.0e-6_wp, "the imposed density gradient must be non-vacuous")
      if (allocated(error)) return
      write (msg, '("depth-mean U/(a*dt) in [",f10.6,",",f10.6,"], analytic 1 +- ",es9.2, &
            &" (",a,", ",a,"): the barotropic mode must feel the depth-mean baroclinic PGF")') &
         rmin, rmax, ratio_tol(), scheme, pgf_form
      call check(error, abs(rmin - 1.0_wp) <= ratio_tol() .and. abs(rmax - 1.0_wp) <= ratio_tol(), &
                 trim(msg))
   end subroutine check_accelerates

   subroutine test_pc_mont(error)
      type(error_type), allocatable, intent(out) :: error
      call check_accelerates(error, "pred_corr", "mont")
   end subroutine test_pc_mont

   subroutine test_ssp_mont(error)
      type(error_type), allocatable, intent(out) :: error
      call check_accelerates(error, "ssp_rk2", "mont")
   end subroutine test_ssp_mont

   subroutine test_pc_fv_mom6(error)
      !! FV-MOM6 closes its pressure stack with `rho_ref*g*eta`, so its
      !! layer PGF CARRIES a free-surface term (unlike Montgomery): the
      !! forcing must shed exactly that and nothing else.
      type(error_type), allocatable, intent(out) :: error
      call check_accelerates(error, "pred_corr", "fv_mom6")
   end subroutine test_pc_fv_mom6

   subroutine test_ssp_fv_mom6(error)
      type(error_type), allocatable, intent(out) :: error
      call check_accelerates(error, "ssp_rk2", "fv_mom6")
   end subroutine test_ssp_fv_mom6

   subroutine test_legacy_control(error)
      !! The defect, kept as a control: `bc_pgf_forcing = .false.` removes
      !! the whole depth-mean PGF from the barotropic forcing and the depth
      !! mean does not move (measured ratio 4e-7, second-order advection,
      !! where the analytic answer is 1).  This arm runs the pre-fix code
      !! path verbatim, so it is the fails-before half of the gate.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: rmin, rmax, amax
      logical :: ok

      call run_case("pred_corr", "mont", .false., rmin, rmax, amax, ok)
      call check(error, ok, "the legacy-split channel must run one step and stay finite")
      if (allocated(error)) return
      call check(error, max(abs(rmin), abs(rmax)) <= 1.0e-3_wp, &
                 "the legacy split must leave the depth mean at rest (it discards the "// &
                 "depth-mean baroclinic PGF) — if this moved, the control no longer "// &
                 "documents the defect")
   end subroutine test_legacy_control

end module test_ocean_bt_baroclinic_forcing
