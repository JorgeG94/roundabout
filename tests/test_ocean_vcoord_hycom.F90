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
!!   - hycom_nonuniform_dsig_flip — NON-uniform dsig (fine surface,
!!     coarse bed): the surface layer uses the FINE dsig (catches a
!!     surface<->bed flip that uniform dsig hides).
!!   - hycom_free_surface_stretching — η /= 0: the z* band scales as
!!     dsig*(H+η), NOT dsig*(H+η)^2/H (catches the stretching-factor
!!     trap the prototype was blind to).
!!   - hycom_conserves — T·h, S·h, total to round-off (single + multi).
!!   - hycom_unstratified_surface_protected — pure-unstratified column:
!!     the surface band is protected (top layers at z*) though the deep
!!     interior collapses.
!!   - hycom_on_device — the full do-concurrent kernel through a device
!!     enter_data round-trip stays finite + conservative.
!!
!! Linear EOS so layer densities are an exact closed form of (T, S):
!! with uniform S, rho = rho0 - alpha_T*(T - T_ref).
module test_ocean_vcoord_hycom
   use rdb_constants, only: wp, REMAP_PPM, VCOORD_HYCOM
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_LINEAR
   use rdb_ocean_vcoord, only: ocean_vcoord_t
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
                  new_unittest("hycom_on_device", test_on_device) &
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
         rmin = rho_lay(NZ) - 0.2_wp
         rmax = rho_lay(1) + 0.2_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         allocate (bt_eta(nx, ny), source=0.0_wp)
         allocate (bt_H_ref(nx, ny), source=H)

         !$acc enter data copyin(ms, vc, bt_eta, bt_H_ref)
         call ms%enter_data()
         call vc%enter_data()
         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                         method=REMAP_PPM, eos=eos)
         !$acc update self(ms%h_layer, ms%tracers(ms%idx_temperature)%hTr)
         !$acc update self(bt_eta)
         call vc%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vc, bt_eta, bt_H_ref)

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

end module test_ocean_vcoord_hycom
