!! VCOORD_HYCOM (hybrid z*/isopycnal) ALE-regrid tests.
!!
!! HYCOM runs the exact VCOORD_RHO density-space inversion
!! (`ocean_vcoord_compute_target_h_rho`, `hybrid=.true.`) plus two
!! deltas: a bottom-up density monotonize before the inversion and a
!! z* nominal-floor sweep after it.  The net effect is a fixed-
!! resolution near-surface z* band (no surface collapse) over an
!! isopycnal interior.
!!
!! These tests target the stage-5 blind spots that the RHO suite (and
!! the η=0 prototype) cannot see:
!!   - hycom_surface_band_zstar — weak surface ML + stratified interior:
!!     the top interface sits at the z* floor; floor invariant holds.
!!   - hycom_nonuniform_dsig_flip — the UNCONFIGURED-slot fallback
!!     (`z_fixed_h_ref = 0`, the historical column-fraction floor) on a
!!     NON-uniform dsig: the surface layer uses the FINE dsig (catches a
!!     surface<->bed flip that uniform dsig hides).
!!   - hycom_nonuniform_profile_flip — the same flip check on the
!!     production floor: a stretched z* profile in METRES
!!     (`ocean_vcoord_set_z_fixed_profile`, as `z_fixed_profile` installs).
!!   - hycom_zstar_floor_metres_shallow_column — a 300 m column under a
!!     600 m profile, weakly stratified: the interfaces sit at the profile
!!     DEPTHS (10, 30, 70, 150 m), not at k/nz of the column (the sigma
!!     floor would make every layer 50 m) — audit finding H1.
!!   - hycom_strong_strat_density_sets — a strongly stratified column
!!     whose isopycnal interfaces all lie below the floor: HYCOM's grid is
!!     the pure RHO grid bit-for-bit (density sets every interface).
!!   - hycom_free_surface_stretching — η /= 0: the z* band scales as
!!     dsig*(H+η), NOT dsig*(H+η)^2/H (catches the stretching-factor
!!     trap the prototype was blind to).
!!   - hycom_conserves — T·h, S·h, total to round-off (single + multi).
!!   - hycom_unstratified_surface_protected — pure-unstratified column:
!!     the surface band is protected (top layers at z*) though the deep
!!     interior collapses.
!!   - hycom_on_device — the full do-concurrent kernel through a device
!!     enter_data round-trip stays finite + conservative.
!!   - hycom_tanh_profile_engine / hycom_refuses_* — the namelist path:
!!     `z_fixed_profile = "tanh"` and `rho_target_profile = "list"` reach
!!     the kernel (a target list lighter than the whole column leaves the
!!     floor binding everywhere, so the surface layer is `z_fixed_dz_top`),
!!     and `validate_config` refuses a bad `rho_target_list` or one set on a
!!     coordinate that never reads it.
!!   - hycom_floor_stress_sweep — 25 344 columns (depth vs z* profile
!!     depth, eta, stratification, target placement, collapsed sources,
!!     three floor sources, both families): no target layer below the
!!     inflation floor, column sum conserved to 1e-12 relative.
!!
!! Linear EOS so layer densities are an exact closed form of (T, S):
!! with uniform S, rho = rho0 - alpha_T*(T - T_ref).
module test_ocean_vcoord_hycom
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp, REMAP_PPM, VCOORD_HYCOM, VCOORD_RHO, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_LINEAR
   use rdb_ocean_vcoord, only: ocean_vcoord_t, ocean_vcoord_set_z_fixed_profile
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_grid_info
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_remap, only: ocean_apply_ale_remap_step
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_hycom_tests

   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3 per degC
   real(wp), parameter :: BETA_S = 0.8_wp      ! kg/m^3 per PSU
   real(wp), parameter :: T_REF = 10.0_wp
   real(wp), parameter :: S_REF = 35.0_wp

contains

   subroutine collect_ocean_vcoord_hycom_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("hycom_surface_band_zstar", test_surface_band_zstar), &
                  new_unittest("hycom_nonuniform_dsig_flip", test_nonuniform_dsig_flip), &
                  new_unittest("hycom_free_surface_stretching", test_free_surface_stretching), &
                  new_unittest("hycom_conserves", test_conserves), &
                  new_unittest("hycom_unstratified_surface_protected", &
                               test_unstratified_surface_protected), &
                  new_unittest("hycom_on_device", test_on_device), &
                  new_unittest("hycom_nonuniform_profile_flip", test_nonuniform_profile_flip), &
                  new_unittest("hycom_zstar_floor_metres_shallow_column", &
                               test_zstar_floor_metres_shallow), &
                  new_unittest("hycom_strong_strat_density_sets", test_strong_strat_density), &
                  new_unittest("hycom_tanh_profile_engine", test_tanh_profile_engine), &
                  new_unittest("hycom_refuses_rho_list_length", test_refuses_rho_list_length), &
                  new_unittest("hycom_refuses_rho_list_not_increasing", &
                               test_refuses_rho_list_order), &
                  new_unittest("hycom_refuses_rho_list_on_sigma", test_refuses_rho_list_sigma), &
                  new_unittest("hycom_refuses_rho_list_without_profile", &
                               test_refuses_rho_list_no_profile), &
                  new_unittest("hycom_floor_stress_sweep", test_floor_stress_sweep) &
                  ]
   end subroutine collect_ocean_vcoord_hycom_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T_REF
      eos%S_ref = S_REF
      eos%is_init = .true.
   end subroutine make_eos

   pure function lin_rho(T, S) result(r)
      real(wp), intent(in) :: T, S
      real(wp) :: r
      r = RHO0 + BETA_S*(S - S_REF) - ALPHA_T*(T - T_REF)
   end function lin_rho

   subroutine setup_column(grid, ms, nz, h_lay, T_lay, S_lay)
      !! Build a uniform-over-(i,j) multilayer state from per-layer
      !! bottom-up (k=1 bed .. k=nz surface) thickness / T / S profiles.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_lay(nz), T_lay(nz), S_lay(nz)
      integer :: i, j, k
      call grid%init(3, 3, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               ms%h_layer(i, j, k) = h_lay(k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = T_lay(k)*h_lay(k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S_lay(k)*h_lay(k)
            end do
         end do
      end do
   end subroutine setup_column

   subroutine run_remap_host(grid, vc, ms, eos, eta_val)
      !! Host-side remap (gfortran test build runs `do concurrent` on
      !! the CPU).  `eta_val` sets a uniform free-surface anomaly; the
      !! column total (sum of h_layer) is the live H+eta the new grid
      !! must span, and bt_H_ref = total - eta is the H reference the
      !! z*-floor sweep stretches.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vcoord_t), intent(inout) :: vc
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: eta_val
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: Htot
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      Htot = sum(ms%h_layer(1, 1, :))
      allocate (bt_eta(nx, ny), source=eta_val)
      allocate (bt_H_ref(nx, ny), source=(Htot - eta_val))
      call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                      method=REMAP_PPM, eos=eos)
      deallocate (bt_eta, bt_H_ref)
   end subroutine run_remap_host

   pure function floor_invariant(h_new, dsig, col_extent, nz) result(ok)
      !! HYCOM floor invariant: cumulative interface depth (bottom-up
      !! state, surface at k=nz) >= cumulative z* = sum dsig*(H+eta).
      !! Both walked surface->bed.  z(state) interfaces accumulate from
      !! the surface (k=nz) down to the bed (k=1).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_new(nz), dsig(nz)
      real(wp), intent(in) :: col_extent
      logical :: ok
      real(wp) :: z, nom
      integer :: k
      ok = .true.
      z = 0.0_wp
      nom = 0.0_wp
      ! Walk surface->bed: state index k = nz, nz-1, ..., 1.
      do k = nz, 1, -1
         z = z + h_new(k)
         nom = nom + dsig(k)*col_extent
         if (z < nom - 1.0e-6_wp) ok = .false.
      end do
   end function floor_invariant

   ! -----------------------------------------------------------------
   ! 1. surface band z* — weak surface ML + stratified interior
   ! -----------------------------------------------------------------
   subroutine test_surface_band_zstar(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_lay(NZ)
      real(wp) :: rmin, rmax, H, hsurf
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         ! bottom-up: k=NZ, NZ-1 surface ML (uniform T), interior cools
         ! toward the bed (k=1 coldest/densest).
         T_lay = [4.0_wp, 6.0_wp, 9.0_wp, 11.0_wp, 12.0_wp, 12.0_wp]
         S_lay = S_REF
         do k = 1, NZ
            h_lay(k) = H/real(NZ, wp)
            rho_lay(k) = lin_rho(T_lay(k), S_REF)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Production floor: uniform z* in METRES (setup writes max_depth
         ! here); max_depth = H, so the nominal layer is H/NZ.
         vc%z_fixed_h_ref = H
         ! uniform dsig (1/NZ) -> z* nominal layer = H/NZ = 100 m.
         rmin = rho_lay(NZ) - 0.2_wp
         rmax = rho_lay(1) + 0.2_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos, 0.0_wp)

         ! Surface layer (k=NZ) must be >= its z* floor (100 m); the weak
         ! ML cannot collapse it.
         hsurf = ms%h_layer(1, 1, NZ)
         call check(error, hsurf >= H/real(NZ, wp) - 1.0e-3_wp, &
                    "surface-band: top layer >= z* floor (H/NZ)")
         if (allocated(error)) exit checks
         call check(error, floor_invariant(ms%h_layer(1, 1, :), vc%dsig, H, NZ), &
                    "surface-band: floor invariant z(k) >= cumulative-z* holds")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "surface-band: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_surface_band_zstar

   ! -----------------------------------------------------------------
   ! 2. NON-uniform dsig — surface uses the fine dsig (flip check)
   ! -----------------------------------------------------------------
   subroutine test_nonuniform_dsig_flip(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: dsig_td(NZ)
      real(wp) :: H, rmin, rmax, hsurf
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         ! Unstratified column so the floor binds on every layer -> the
         ! grid is pure z*, exposing the dsig orientation directly.
         T_lay = 12.0_wp
         S_lay = S_REF
         do k = 1, NZ
            h_lay(k) = H/real(NZ, wp)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! NON-uniform dsig, TOP-DOWN (surface->bed): fine surface, coarse
         ! bed.  The slot stores dsig BOTTOM-UP (dsig(NZ)=surface), so flip.
         dsig_td = [0.05_wp, 0.05_wp, 0.10_wp, 0.15_wp, 0.25_wp, 0.40_wp]
         do k = 1, NZ
            vc%dsig(k) = dsig_td(NZ - k + 1)   ! bottom-up storage
         end do
         rmin = lin_rho(12.0_wp, S_REF) - 1.0_wp
         rmax = lin_rho(12.0_wp, S_REF) + 1.0_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos, 0.0_wp)

         ! Surface layer (k=NZ) must use the FINE surface dsig (0.05*600 =
         ! 30 m), NOT the coarse bed value (0.40*600 = 240 m).
         hsurf = ms%h_layer(1, 1, NZ)
         call check(error, abs(hsurf - dsig_td(1)*H) < 1.0_wp, &
                    "nonuniform-dsig: surface uses FINE dsig (~30 m), not coarse bed")
         if (allocated(error)) exit checks
         call check(error, floor_invariant(ms%h_layer(1, 1, :), vc%dsig, H, NZ), &
                    "nonuniform-dsig: floor invariant holds")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "nonuniform-dsig: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_nonuniform_dsig_flip

   ! -----------------------------------------------------------------
   ! 3. free surface (eta /= 0) — z* band scales as dsig*(H+eta)
   ! -----------------------------------------------------------------
   subroutine test_free_surface_stretching(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: H, eta, col_extent, rmin, rmax, hsurf, expect_correct, expect_trap
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         eta = 60.0_wp     ! 10% free-surface rise
         col_extent = H + eta
         ! Unstratified -> floor binds everywhere -> pure z* grid; the
         ! surface layer exposes the stretching factor exactly.
         T_lay = 12.0_wp
         S_lay = S_REF
         ! Column total = H + eta (the live thickness the IC carries).
         do k = 1, NZ
            h_lay(k) = col_extent/real(NZ, wp)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Production floor: uniform z* in METRES (setup writes max_depth
         ! here); max_depth = H, so the nominal layer is H/NZ.
         vc%z_fixed_h_ref = H
         ! uniform dsig = 1/NZ.
         rmin = lin_rho(12.0_wp, S_REF) - 1.0_wp
         rmax = lin_rho(12.0_wp, S_REF) + 1.0_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos, eta)

         ! CORRECT: surface layer ~ dsig*(H+eta) = (1/NZ)*660 = 110 m.
         ! TRAP   : dsig*(H+eta)^2/H = (1/NZ)*660*1.1 = 121 m.
         hsurf = ms%h_layer(1, 1, NZ)
         expect_correct = col_extent/real(NZ, wp)
         expect_trap = col_extent/real(NZ, wp)*(col_extent/H)
         call check(error, abs(hsurf - expect_correct) < 0.5_wp, &
                    "free-surface: surface layer ~ dsig*(H+eta) (correct factor)")
         if (allocated(error)) exit checks
         call check(error, abs(hsurf - expect_trap) > 5.0_wp, &
                    "free-surface: surface layer is NOT dsig*(H+eta)^2/H (over-stretch trap)")
         if (allocated(error)) exit checks
         ! Floor invariant uses col_extent = H + eta.
         call check(error, floor_invariant(ms%h_layer(1, 1, :), vc%dsig, col_extent, NZ), &
                    "free-surface: floor invariant z(k) >= cumulative dsig*(H+eta)")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - col_extent) < 1.0e-8_wp, &
                    "free-surface: column total conserved (= H+eta)")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_free_surface_stretching

   ! -----------------------------------------------------------------
   ! 4. conservation — T*h, S*h, total to round-off (single + multi)
   ! -----------------------------------------------------------------
   subroutine test_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_lay(NZ)
      real(wp) :: rmin, rmax, H
      real(wp) :: Th0, Sh0, Th1, Sh1, Th2, Sh2
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         T_lay = [4.0_wp, 6.0_wp, 8.0_wp, 10.0_wp, 12.0_wp, 12.0_wp]
         S_lay = [34.0_wp, 34.5_wp, 35.0_wp, 35.0_wp, 35.5_wp, 36.0_wp]
         do k = 1, NZ
            h_lay(k) = H/real(NZ, wp)
            rho_lay(k) = lin_rho(T_lay(k), S_lay(k))
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Production floor: uniform z* in METRES (setup writes max_depth
         ! here); max_depth = H, so the nominal layer is H/NZ.
         vc%z_fixed_h_ref = H
         rmin = rho_lay(NZ) - 0.5_wp
         rmax = rho_lay(1) + 0.5_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         ! single regrid
         call run_remap_host(grid, vc, ms, eos, 0.0_wp)
         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh1 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-9_wp, "conserves: T*h (single regrid)")
         if (allocated(error)) exit checks
         call check(error, abs(Sh1 - Sh0) < 1.0e-9_wp, "conserves: S*h (single regrid)")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "conserves: total (single regrid)")
         if (allocated(error)) exit checks

         ! second regrid (the collapsed/floored column is now h_old)
         call run_remap_host(grid, vc, ms, eos, 0.0_wp)
         Th2 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh2 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th2 - Th0) < 1.0e-9_wp, "conserves: T*h (multi regrid)")
         if (allocated(error)) exit checks
         call check(error, abs(Sh2 - Sh0) < 1.0e-9_wp, "conserves: S*h (multi regrid)")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "conserves: total (multi regrid)")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_conserves

   ! -----------------------------------------------------------------
   ! 5. unstratified column — surface band protected though deep collapses
   ! -----------------------------------------------------------------
   subroutine test_unstratified_surface_protected(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_mean, H, zstar
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         T_lay = 12.0_wp     ! pure unstratified
         S_lay = S_REF
         do k = 1, NZ
            h_lay(k) = H/real(NZ, wp)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Production floor: uniform z* in METRES (setup writes max_depth
         ! here); max_depth = H, so the nominal layer is H/NZ.
         vc%z_fixed_h_ref = H
         ! Targets bracketing the (uniform) column density: RHO alone would
         ! collapse every interior interface to the surface or bed.
         rho_mean = lin_rho(12.0_wp, S_REF)
         do k = 0, NZ
            vc%rho_target(k) = rho_mean - 1.0_wp + 2.0_wp*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos, 0.0_wp)

         ! Surface band protected: top two layers at the z* floor (~100 m).
         zstar = H/real(NZ, wp)
         call check(error, abs(ms%h_layer(1, 1, NZ) - zstar) < 1.0_wp, &
                    "unstratified: top layer at z* floor (surface protected)")
         if (allocated(error)) exit checks
         call check(error, abs(ms%h_layer(1, 1, NZ - 1) - zstar) < 1.0_wp, &
                    "unstratified: 2nd layer at z* floor (surface band protected)")
         if (allocated(error)) exit checks
         call check(error, floor_invariant(ms%h_layer(1, 1, :), vc%dsig, H, NZ), &
                    "unstratified: floor invariant holds (no surface collapse)")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "unstratified: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_unstratified_surface_protected

   ! -----------------------------------------------------------------
   ! 6. on-device round-trip — kernel runs through an enter_data round-trip
   ! -----------------------------------------------------------------
   subroutine test_on_device(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_lay(NZ)
      real(wp) :: rmin, rmax, H, Hsum, Th0, Th1
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      integer :: k, nx, ny
      logical :: finite
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         T_lay = [4.0_wp, 6.0_wp, 9.0_wp, 11.0_wp, 12.0_wp, 12.0_wp]
         S_lay = S_REF
         do k = 1, NZ
            h_lay(k) = H/real(NZ, wp)
            rho_lay(k) = lin_rho(T_lay(k), S_REF)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         nx = grid%nx_total
         ny = grid%ny_total

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Production floor: uniform z* in METRES (setup writes max_depth
         ! here); max_depth = H, so the nominal layer is H/NZ.
         vc%z_fixed_h_ref = H
         rmin = rho_lay(NZ) - 0.2_wp
         rmax = rho_lay(1) + 0.2_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         allocate (bt_eta(nx, ny), source=0.0_wp)
         allocate (bt_H_ref(nx, ny), source=H)

         !$omp target enter data map(to: ms, vc, bt_eta, bt_H_ref)
         call ms%enter_data()
         call vc%enter_data()
         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                         method=REMAP_PPM, eos=eos)
         !$omp target update from(ms%h_layer, ms%tracers(ms%idx_temperature)%hTr)
         !$omp target update from(bt_eta)
         call vc%exit_data()
         call ms%exit_data()
         !$omp target exit data map(delete: ms, vc, bt_eta, bt_H_ref)

         finite = .true.
         do k = 1, NZ
            if (ms%h_layer(1, 1, k) /= ms%h_layer(1, 1, k)) finite = .false.
            if (ms%h_layer(1, 1, k) < 0.0_wp) finite = .false.
         end do
         call check(error, finite, "device: thicknesses finite + non-negative")
         if (allocated(error)) exit checks
         Hsum = sum(ms%h_layer(1, 1, :))
         call check(error, abs(Hsum - H) < 1.0e-7_wp, "device: column total conserved")
         if (allocated(error)) exit checks
         call check(error, floor_invariant(ms%h_layer(1, 1, :), vc%dsig, H, NZ), &
                    "device: floor invariant holds")
         if (allocated(error)) exit checks
         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-7_wp, "device: T*h conserved")
      end block checks
      if (allocated(bt_eta)) deallocate (bt_eta)
      if (allocated(bt_H_ref)) deallocate (bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_on_device

   ! -----------------------------------------------------------------
   ! 7. NON-uniform z* PROFILE in metres — surface uses the fine entry
   ! -----------------------------------------------------------------
   subroutine test_nonuniform_profile_flip(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: dz_td(NZ), H, rmin, rmax, hsurf
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         ! Unstratified: the floor binds on every layer -> pure z* grid.
         T_lay = 12.0_wp
         S_lay = S_REF
         h_lay = H/real(NZ, wp)
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         ! Surface-first nominal thicknesses (m): fine surface, coarse bed.
         dz_td = [30.0_wp, 30.0_wp, 60.0_wp, 90.0_wp, 150.0_wp, 240.0_wp]
         call ocean_vcoord_set_z_fixed_profile(vc, dz_td)
         ! Every target lighter than the column: each interior interface
         ! inverts to the surface and the floor alone places it.
         rmin = lin_rho(12.0_wp, S_REF) - 5.0_wp
         rmax = lin_rho(12.0_wp, S_REF) - 4.0_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos, 0.0_wp)

         hsurf = ms%h_layer(1, 1, NZ)
         call check(error, abs(hsurf - dz_td(1)) < 1.0e-6_wp, &
                    "profile: surface layer is the FINE 30 m entry, not the 240 m bed one")
         if (allocated(error)) exit checks
         do k = 1, NZ
            call check(error, abs(ms%h_layer(1, 1, k) - dz_td(NZ - k + 1)) < 1.0e-6_wp, &
                       "profile: every layer sits on its nominal z* thickness")
            if (allocated(error)) exit checks
         end do
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-8_wp, &
                    "profile: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_nonuniform_profile_flip

   ! -----------------------------------------------------------------
   ! 8. z* floor in METRES on a column shallower than the profile (H1)
   ! -----------------------------------------------------------------
   subroutine test_zstar_floor_metres_shallow(error)
      !! A 300 m column under a 600 m surface-first profile
      !! 10/20/40/80/150/300 m, weakly stratified (0.05 degC over the
      !! column, every target far lighter than it).  The density inversion
      !! leaves every interior interface at the surface, so the z* floor
      !! places them: at 10, 30, 70, 150, 300 (= the bed) m — layers
      !! 10/20/40/80/150 m and a collapsed bed layer.  The pre-fix sigma
      !! floor (k/nz of the column) made every layer 50 m.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ), dz_td(NZ)
      real(wp) :: H, rho_c, h_floor_eff
      integer :: k
      checks: block
         call make_eos(eos)
         H = 300.0_wp
         do k = 1, NZ
            T_lay(k) = 10.0_wp + 0.01_wp*real(k - 1, wp)   ! bed coldest
         end do
         S_lay = S_REF
         h_lay = H/real(NZ, wp)
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_HYCOM
         dz_td = [10.0_wp, 20.0_wp, 40.0_wp, 80.0_wp, 150.0_wp, 300.0_wp]
         call ocean_vcoord_set_z_fixed_profile(vc, dz_td)
         ! Every target LIGHTER than the column: the inversion puts every
         ! interior interface at the surface and the floor sets them all.
         rho_c = lin_rho(T_lay(NZ), S_REF)
         do k = 0, NZ
            vc%rho_target(k) = rho_c - 5.0_wp + 0.1_wp*real(k, wp)
         end do
         h_floor_eff = max(vc%zstar_h_min, 3.0e-4_wp)

         call run_remap_host(grid, vc, ms, eos, 0.0_wp)

         ! state k = NZ is the surface; layers NZ..3 carry dz_td(1..4).
         do k = 1, 4
            call check(error, abs(ms%h_layer(1, 1, NZ - k + 1) - dz_td(k)) < 1.0e-9_wp, &
                       "shallow: interface at the z* profile DEPTH (metres), not k/nz of H")
            if (allocated(error)) exit checks
         end do
         call check(error, abs(ms%h_layer(1, 1, NZ) - H/real(NZ, wp)) > 1.0_wp, &
                    "shallow: the surface layer is NOT the sigma floor's H/nz = 50 m")
         if (allocated(error)) exit checks
         ! The bed layer collapsed and was inflated; the 150 m layer paid.
         call check(error, abs(ms%h_layer(1, 1, 1) - h_floor_eff) < 1.0e-12_wp, &
                    "shallow: the layer below the bed clamp sits at the inflation floor")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H) < 1.0e-9_wp, &
                    "shallow: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_zstar_floor_metres_shallow

   ! -----------------------------------------------------------------
   ! 9. strongly stratified column — density sets every interface
   ! -----------------------------------------------------------------
   subroutine test_strong_strat_density(error)
      !! Linear stable stratification (2 degC per 100 m layer), targets
      !! spanning the column, z* floor 1 m per layer: every isopycnal
      !! interface lies far below the floor and the monotonize is a no-op
      !! on the stable column, so HYCOM must reproduce the pure RHO grid
      !! bit-for-bit — the floor only ever acts where density does not.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_h, ms_r
      type(ocean_vcoord_t) :: vc_h, vc_r
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: H, rtop, rbed, rmov
      integer :: k
      checks: block
         call make_eos(eos)
         H = 600.0_wp
         do k = 1, NZ
            T_lay(k) = 2.0_wp + 2.0_wp*real(k - 1, wp)   ! bed coldest
         end do
         S_lay = S_REF
         h_lay = H/real(NZ, wp)
         rtop = lin_rho(T_lay(NZ), S_REF)
         rbed = lin_rho(T_lay(1), S_REF)
         call setup_column(grid, ms_h, NZ, h_lay, T_lay, S_lay)
         call setup_column(grid, ms_r, NZ, h_lay, T_lay, S_lay)
         call vc_h%init(grid, nz_ml=NZ)
         call vc_r%init(grid, nz_ml=NZ)
         vc_h%coord_type = VCOORD_HYCOM
         vc_r%coord_type = VCOORD_RHO
         vc_h%z_fixed_h_ref = real(NZ, wp)      ! 1 m nominal layers
         do k = 0, NZ
            ! Targets shifted off the layer densities so the interfaces
            ! move off the layer boundaries (non-vacuous).
            rmov = rtop + (rbed - rtop)*(real(k, wp) - 0.33_wp)/real(NZ - 1, wp)
            vc_h%rho_target(k) = rmov
            vc_r%rho_target(k) = rmov
         end do

         call run_remap_host(grid, vc_h, ms_h, eos, 0.0_wp)
         call run_remap_host(grid, vc_r, ms_r, eos, 0.0_wp)

         call check(error, maxval(abs(ms_r%h_layer(1, 1, :) - h_lay)) > 1.0_wp, &
                    "strong-strat: the RHO grid moved (test is not vacuous)")
         if (allocated(error)) exit checks
         call check(error, all(ms_h%h_layer(1, 1, :) == ms_r%h_layer(1, 1, :)), &
                    "strong-strat: HYCOM == RHO bit-for-bit (density sets every interface)")
      end block checks
      call vc_h%destroy()
      call vc_r%destroy()
      call ms_h%destroy()
      call ms_r%destroy()
   end subroutine test_strong_strat_density

   ! -----------------------------------------------------------------
   ! 10. namelist path — tanh z* floor + target list reach the kernel
   ! -----------------------------------------------------------------
   function engine_nml(vcoord, extra) result(txt)
      !! Flat 600 m, 6 layers, Wright EOS (sigma-0 coordinate), salinity
      !! stratified by the `&tracer_nml` linear-in-layer IC (no NetCDF).
      character(len=*), intent(in) :: vcoord, extra
      character(len=:), allocatable :: txt
      character(len=*), parameter :: nl = new_line("a")
      txt = '&sim_nml sim_type = "ocean" /'//nl// &
            "&grid_nml nx = 6, ny = 4, dx = 2000.0, dy = 2000.0, nghost = 2 /"//nl// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//nl// &
            "&physics_nml coriolis_f = -1.409e-4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 6 /"//nl// &
            '&vcoord_nml vcoord_type = "'//vcoord//'", rho_ref_pressure = 0.0, '// &
            extra//" /"//nl// &
            '&ocean_topo_nml topo_config = "flat", max_depth = 600.0 /'//nl// &
            '&ocean_pgf_nml form = "fv_mom6", reconstruct_for_pressure = .true. /'//nl// &
            '&ocean_eos_nml eos = "wright" /'//nl// &
            "&tracer_nml initial_temperature = -1.9, initial_salinity = 33.8, "// &
            "S_init_surface = 33.8, S_init_bottom = 34.55 /"//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function engine_nml

   subroutine test_tanh_profile_engine(error)
      !! `z_fixed_profile = "tanh"` from a 10 m surface layer and a target
      !! LIST lighter than the whole column (sigma-0 ~1027.2-1027.9 here):
      !! every interior interface inverts to the surface, the floor sets all
      !! of them, so after a regrid the surface layer is 10 m (x (H+eta)/H,
      !! eta ~ 0 at rest) — the profile reached the kernel through setup.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :)
      real(wp) :: err_top
      character(len=120) :: msg
      nml = engine_nml("hycom", "z_fixed_profile = 'tanh', z_fixed_dz_top = 10.0, "// &
                       "rho_target_profile = 'list', rho_target_list = 1020.0, 1020.5, "// &
                       "1021.0, 1021.5, 1022.0, 1022.5, 1023.0")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "tanh+list hycom case must build")
      if (allocated(error)) return
      body: block
         status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
         status = rdb_ocean_step(handle, 3_c_int)
         call check(error, status == OCEAN_STATUS_OK, "tanh+list hycom case must step")
         if (allocated(error)) exit body
         status = rdb_ocean_refresh_host(handle)
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h, [nx, ny, nz])
         err_top = maxval(abs(h(ng + 1:ng + nx_p, ng + 1:ng + ny_p, nz) - 10.0_wp))
         write (msg, "(a,es10.3,a)") "surface layer is z_fixed_dz_top = 10 m (max err ", &
            err_top, " m)"
         call check(error, err_top < 1.0e-3_wp, trim(msg))
      end block body
      status = rdb_ocean_destroy(handle)
   end subroutine test_tanh_profile_engine

   subroutine expect_refused(error, vcoord, extra, what)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: vcoord, extra, what
      character(len=:), allocatable :: nml
      type(c_ptr) :: handle
      integer(c_int) :: status
      nml = engine_nml(vcoord, extra)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status /= OCEAN_STATUS_OK, what//" built")
      if (status == OCEAN_STATUS_OK) status = rdb_ocean_destroy(handle)
   end subroutine expect_refused

   subroutine test_refuses_rho_list_length(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "hycom", "rho_target_profile = 'list', "// &
                          "rho_target_list = 1020.0, 1021.0, 1022.0", &
                          "a 3-entry rho_target_list for 6 layers")
   end subroutine test_refuses_rho_list_length

   subroutine test_refuses_rho_list_order(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "rho", "rho_target_profile = 'list', "// &
                          "rho_target_list = 1020.0, 1021.0, 1022.0, 1022.0, 1023.0, "// &
                          "1024.0, 1025.0", "a non-increasing rho_target_list")
   end subroutine test_refuses_rho_list_order

   subroutine test_refuses_rho_list_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "sigma", "rho_target_profile = 'list', "// &
                          "rho_target_list = 1020.0, 1021.0, 1022.0, 1023.0, 1024.0, "// &
                          "1025.0, 1026.0", "rho_target_profile='list' on sigma")
   end subroutine test_refuses_rho_list_sigma

   subroutine test_refuses_rho_list_no_profile(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, "hycom", "rho_target_list = 1020.0, 1021.0, 1022.0, "// &
                          "1023.0, 1024.0, 1025.0, 1026.0", &
                          "rho_target_list under rho_target_profile='uniform'")
   end subroutine test_refuses_rho_list_no_profile

   ! -----------------------------------------------------------------
   ! 11. floor stress sweep — no configuration goes below the floor
   ! -----------------------------------------------------------------
   subroutine test_floor_stress_sweep(error)
      !! Sweep the RHO/HYCOM target builder over every shape that could
      !! defeat the z* floor sweep + inflation: column depth against the
      !! profile depth (just below / at / just above `nz·floor` and the
      !! profile's interface depths, up to deeper than the profile), eta
      !! of 0 / +-0.5 / +3 (scaled to the column), uniform / stable /
      !! unstable / two-layer columns, targets lighter than / bracketing /
      !! denser than / inside the column, uniform / bed-collapsed /
      !! surface-collapsed source layers; uniform-metres, stretched and
      !! unconfigured floors; both families.  Every target layer must be at
      !! or above the inflation floor (or the column left at h_old — the
      !! too-thin guard) and the column sum conserved to 1e-12 relative.
      !! 25 344 columns; the target builder alone, no remap.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 15
      type(hgrid_t) :: grid
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      real(wp), allocatable :: eta(:, :)
      real(wp) :: hs(22), etas(4), dz(NZ), hsrc(NZ), t_col(NZ)
      real(wp) :: f, tot, hmin, rlo, rhi, hc, e
      integer :: ih, ie, it, ip, ig, k, nbad, ncase, prof, coord
      character(len=240) :: first_bad
      hs = [0.003_wp, 0.0045_wp, 0.006_wp, 1.0_wp, 4.0_wp, 9.99_wp, 10.0_wp, 10.01_wp, &
            30.0_wp, 70.0_wp, 100.0_wp, 149.9_wp, 150.0_wp, 300.0_wp, 449.0_wp, 450.0_wp, &
            451.0_wp, 600.0_wp, 900.0_wp, 1049.0_wp, 1050.0_wp, 2000.0_wp]
      etas = [0.0_wp, 0.5_wp, -0.5_wp, 3.0_wp]
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = 1027.51_wp
      eos%alpha_T = 3.8356948e-2_wp
      eos%beta_S = 8.0587609e-1_wp
      eos%T_ref = -1.0_wp
      eos%S_ref = 34.2_wp
      eos%is_init = .true.
      call grid%init(1, 1, 0, 1.0_wp, 1.0_wp)
      allocate (eta(grid%nx_total, grid%ny_total))
      nbad = 0
      ncase = 0
      first_bad = ""
      do coord = 1, 2
         do prof = 1, 3
            call vc%init(grid, nz_ml=NZ)
            vc%coord_type = merge(VCOORD_HYCOM, VCOORD_RHO, coord == 1)
            vc%zstar_h_min = 1.0e-4_wp
            select case (prof)
            case (1)      ! uniform metres: 70 m layers over 1050 m
               vc%z_fixed_h_ref = 1050.0_wp
            case (2)      ! stretched: 2 m x 1.35**k, surface first
               do k = 1, NZ
                  dz(k) = 2.0_wp*1.35_wp**(k - 1)
               end do
               call ocean_vcoord_set_z_fixed_profile(vc, dz)
            case default  ! unconfigured column-fraction fallback
               vc%z_fixed_h_ref = 0.0_wp
            end select
            f = max(vc%zstar_h_min, 2.0_wp*H_VANISHED)
            do ih = 1, size(hs)
               do ie = 1, size(etas)
                  do it = 1, 4
                     do ip = 1, 4
                        do ig = 1, 3
                           hc = hs(ih)
                           e = etas(ie)*min(1.0_wp, hc)
                           if (hc + e <= 0.0_wp) cycle
                           select case (ig)
                           case (1)
                              hsrc = (hc + e)/real(NZ, wp)
                           case (2)      ! bed-collapsed (state k=1 = bed)
                              hsrc = 3.0e-4_wp
                              hsrc(NZ) = max(hc + e - real(NZ - 1, wp)*3.0e-4_wp, 3.0e-4_wp)
                           case default  ! surface-collapsed
                              hsrc = 3.0e-4_wp
                              hsrc(1) = max(hc + e - real(NZ - 1, wp)*3.0e-4_wp, 3.0e-4_wp)
                           end select
                           do k = 1, NZ
                              select case (it)
                              case (1)
                                 t_col(k) = -1.9_wp
                              case (2)
                                 t_col(k) = -1.9_wp + 0.2_wp*real(k - 1, wp)
                              case (3)
                                 t_col(k) = 1.0_wp - 0.2_wp*real(k - 1, wp)
                              case default
                                 t_col(k) = merge(-1.9_wp, 1.0_wp, k <= NZ/2)
                              end select
                           end do
                           select case (ip)
                           case (1)
                              rlo = 1020.0_wp
                              rhi = 1021.0_wp
                           case (2)
                              rlo = 1027.3_wp
                              rhi = 1027.8_wp
                           case (3)
                              rlo = 1035.0_wp
                              rhi = 1036.0_wp
                           case default
                              rlo = 1027.55_wp
                              rhi = 1027.62_wp
                           end select
                           do k = 0, NZ
                              vc%rho_target(k) = rlo + (rhi - rlo)*real(k, wp)/real(NZ, wp)
                           end do
                           do k = 1, NZ
                              vc%remap_h_old(:, :, k) = hsrc(k)
                              vc%remap_conc_t(:, :, k) = t_col(k)
                           end do
                           vc%remap_conc_s = 34.2_wp
                           tot = sum(hsrc)
                           eta = e
                           vc%remap_h_ref = tot - e
                           call vc%compute_target_h_rho(vc%remap_h_ref, eta, vc%remap_conc_t, &
                                                        vc%remap_conc_s, eos, &
                                                        hybrid=(coord == 1))
                           ncase = ncase + 1
                           hmin = minval(vc%target_h(1, 1, :))
                           if ((hmin < f*(1.0_wp - 1.0e-12_wp) .and. &
                                any(vc%target_h(1, 1, :) /= hsrc)) .or. &
                               abs(sum(vc%target_h(1, 1, :)) - tot) > &
                               1.0e-12_wp*max(tot, 1.0_wp)) then
                              nbad = nbad + 1
                              if (nbad == 1) then
                                 write (first_bad, '(a,i0,a,i0,a,es10.3,a,es10.3,3(a,i0),a,es11.3,a,es11.3)') &
                                    "family ", coord, " profile ", prof, " H ", hc, " eta ", e, &
                                    " T ", it, " targets ", ip, " source ", ig, " hmin ", hmin, &
                                    " dsum ", sum(vc%target_h(1, 1, :)) - tot
                              end if
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            call vc%destroy()
         end do
      end do
      call check(error, ncase > 20000, "sweep: the case grid ran (non-vacuous)")
      if (allocated(error)) return
      call check(error, nbad == 0, "sweep: a target layer went below the floor or the "// &
                 "column sum moved; first: "//trim(first_bad))
   end subroutine test_floor_stress_sweep

end module test_ocean_vcoord_hycom
