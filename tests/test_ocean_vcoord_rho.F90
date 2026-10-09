!! VCOORD_RHO (isopycnal) ALE-regrid tests.
!!
!! Exercises the density-space interface placement in
!! `ocean_vcoord_compute_target_h_rho` + the conservative remap it
!! drives through `ocean_apply_ale_remap_step`.  The coordinate places
!! layer interfaces on prescribed potential-density surfaces
!! `rho_target(0:nz)` (lightest = surface, densest = bed) by inverting a
!! PPM reconstruction of the column density profile.
!!
!! Cases (analytical — the high-leverage kind):
!!   1. two_layer_interface_lands_on_jump — a 2-layer warm-over-cold
!!      column with the interior target at the density step: the new
!!      interface lands on the layer boundary and T·h / S·h conserve.
!!   2. linear_strat_targets_recovered — linearly stratified column,
!!      targets spanning the range: the post-remap layer densities track
!!      `rho_target` and the column integrals conserve.
!!   3. unstable_column_graceful — near-neutral + a small inversion:
!!      no NaN, finite thicknesses summing to the column total, conserve.
!!   4. regrid_runs_on_device — the full `do concurrent` inversion runs
!!      through a device enter_data round-trip (GPU exercise) and stays
!!      finite + conservative.
!!   5. fully_vanished_column — all but one source layer below the
!!      min-thickness floor: the `<=1`-survivor fast path keeps
!!      `h_new = h_old` exactly.
!!   7. thin_column_{land,minted,shared_debit}_{rho,hycom} — columns too
!!      thin for every layer to sit at the inflation floor (`nz` x
!!      `2*H_VANISHED`), and a column just thick enough whose thickest
!!      layer alone cannot pay the inflation debit.  The step-5 inflation
!!      used to write a NEGATIVE thickness (a 50 x 1.5e-4 m land column
!!      collapsed into one layer gave -7.2e-3 m, the 1-degree Southern
!!      Ocean `rho` crash) or mint mass (`ns == 0`).  Asserts h >= 0,
!!      every layer at or above the floor (or untouched), and h, T*h, S*h
!!      conserved to round-off, on both density families.
!!
!! Uses the linear EOS so layer densities are an exact closed form of
!! (T, S): with uniform S, rho = rho0 - alpha_T*(T - T_ref).
module test_ocean_vcoord_rho
   use rdb_constants, only: wp, REMAP_PPM, VCOORD_RHO, VCOORD_HYCOM, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_LINEAR
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_remap, only: ocean_apply_ale_remap_step
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_rho_tests

   ! Linear-EOS coefficients used by every test: realistic thermal
   ! expansion so density steps are large relative to PPM round-off;
   ! salinity contraction off (S held uniform) so rho is a clean
   ! function of T alone.
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3 per degC
   real(wp), parameter :: BETA_S = 0.8_wp      ! kg/m^3 per PSU
   real(wp), parameter :: T_REF = 10.0_wp
   real(wp), parameter :: S_REF = 35.0_wp

contains

   subroutine collect_ocean_vcoord_rho_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("two_layer_interface_lands_on_jump", test_two_layer_jump), &
                  new_unittest("linear_strat_targets_recovered", test_linear_strat), &
                  new_unittest("unstable_column_graceful", test_unstable), &
                  new_unittest("regrid_runs_on_device", test_on_device), &
                  new_unittest("fully_vanished_column", test_vanished), &
                  new_unittest("multi_regrid_conserves", test_multi_regrid), &
                  new_unittest("thin_column_land_rho", test_thin_land_rho), &
                  new_unittest("thin_column_land_hycom", test_thin_land_hycom), &
                  new_unittest("thin_column_minted_rho", test_thin_minted_rho), &
                  new_unittest("thin_column_minted_hycom", test_thin_minted_hycom), &
                  new_unittest("thin_column_shared_debit_rho", test_shared_debit_rho) &
                  ]
   end subroutine collect_ocean_vcoord_rho_tests

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
      !! Mirror of the linear EOS for the test's analytic predictions.
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

   subroutine run_remap_host(grid, vc, ms, eos)
      !! Host-side remap (gfortran test build runs `do concurrent` on
      !! the CPU; no enter_data needed).
      type(hgrid_t), intent(in) :: grid
      type(ocean_vcoord_t), intent(inout) :: vc
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: Htot
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      Htot = sum(ms%h_layer(1, 1, :))
      allocate (bt_eta(nx, ny), source=0.0_wp)
      allocate (bt_H_ref(nx, ny), source=Htot)
      call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                      method=REMAP_PPM, eos=eos)
      deallocate (bt_eta, bt_H_ref)
   end subroutine run_remap_host

   ! -----------------------------------------------------------------
   ! 1. two-layer column — interface lands on the density jump
   ! -----------------------------------------------------------------
   subroutine test_two_layer_jump(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 2
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_cold, rho_warm, Th0, Sh0, Th1, Sh1
      real(wp) :: hsurf
      checks: block
         call make_eos(eos)
         ! bottom-up: k=1 bed = cold/dense, k=2 surface = warm/light.
         h_lay = [40.0_wp, 60.0_wp]
         T_lay = [4.0_wp, 18.0_wp]
         S_lay = [S_REF, S_REF]
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         rho_warm = lin_rho(18.0_wp, S_REF)
         rho_cold = lin_rho(4.0_wp, S_REF)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         ! rho_target(0)=light(surface) .. rho_target(NZ)=dense(bed).
         ! interior target = midpoint of the jump -> interface on z=h_surf.
         vc%rho_target(0) = rho_warm - 1.0_wp
         vc%rho_target(1) = 0.5_wp*(rho_warm + rho_cold)
         vc%rho_target(2) = rho_cold + 1.0_wp

         call run_remap_host(grid, vc, ms, eos)

         ! Surface (k=NZ) thickness should be ~ the warm layer thickness
         ! (60 m): the midpoint target lands on the warm/cold boundary.
         hsurf = ms%h_layer(1, 1, NZ)
         call check(error, abs(hsurf - 60.0_wp) < 1.0e-6_wp, &
                    "two-layer: surface layer should equal the warm-layer thickness")
         if (allocated(error)) exit checks

         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh1 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-9_wp, "two-layer: T*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(Sh1 - Sh0) < 1.0e-9_wp, "two-layer: S*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - 100.0_wp) < 1.0e-9_wp, &
                    "two-layer: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_two_layer_jump

   ! -----------------------------------------------------------------
   ! 2. linearly stratified column — target densities recovered
   ! -----------------------------------------------------------------
   subroutine test_linear_strat(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 6
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rho_lay(NZ)
      real(wp) :: Th0, Sh0, Th1, Sh1
      real(wp) :: rho_new, Tk, Sk, rmin, rmax
      integer :: k
      checks: block
         call make_eos(eos)
         ! Cools with depth (k=1 bed coldest) -> denser at the bed.
         do k = 1, NZ
            h_lay(k) = 100.0_wp/real(NZ, wp)
            ! T increases toward the surface (k=NZ).
            T_lay(k) = 4.0_wp + 14.0_wp*real(k - 1, wp)/real(NZ - 1, wp)
            S_lay(k) = S_REF
            rho_lay(k) = lin_rho(T_lay(k), S_REF)
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         ! Targets span the column's density range (lightest=surface).
         rmin = rho_lay(NZ)   ! surface (warmest, lightest)
         rmax = rho_lay(1)    ! bed (coldest, densest)
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos)

         ! Recovered layer densities must (a) stay in the target range,
         ! (b) be MONOTONE — denser toward the bed (k=1) — and (c) track
         ! the target-interface midpoints, so a constant/flat column would
         ! FAIL (the weak "within range" check alone is tautological).
         do k = NZ, 1, -1           ! walk bed-up: density must not increase
            Tk = ms%tracers(ms%idx_temperature)%hTr(1, 1, k)/ms%h_layer(1, 1, k)
            Sk = ms%tracers(ms%idx_salinity)%hTr(1, 1, k)/ms%h_layer(1, 1, k)
            rho_new = lin_rho(Tk, Sk)
            call check(error, rho_new >= rmin - 1.0e-6_wp .and. rho_new <= rmax + 1.0e-6_wp, &
                       "linear strat: recovered density within target range")
            if (allocated(error)) exit checks
            ! Surface↔bed flip: bottom-up layer k (k=1 bed) sits between
            ! target interfaces NZ-k (lighter, above) and NZ-k+1 (denser,
            ! below).  Its mean density should land near that midpoint.
            call check(error, abs(rho_new &
                                  - 0.5_wp*(vc%rho_target(NZ - k) + vc%rho_target(NZ - k + 1))) &
                       < 0.15_wp*(rmax - rmin), &
                       "linear strat: recovered density tracks the target-interface midpoint")
            if (allocated(error)) exit checks
         end do
         ! Real stratification preserved (a flat column collapses this spread).
         Tk = ms%tracers(ms%idx_temperature)%hTr(1, 1, 1)/ms%h_layer(1, 1, 1)
         Sk = ms%tracers(ms%idx_salinity)%hTr(1, 1, 1)/ms%h_layer(1, 1, 1)
         rho_new = lin_rho(Tk, Sk)
         Tk = ms%tracers(ms%idx_temperature)%hTr(1, 1, NZ)/ms%h_layer(1, 1, NZ)
         Sk = ms%tracers(ms%idx_salinity)%hTr(1, 1, NZ)/ms%h_layer(1, 1, NZ)
         call check(error, (rho_new - lin_rho(Tk, Sk)) > 0.5_wp*(rmax - rmin), &
                    "linear strat: bed-to-surface density spread preserved")
         if (allocated(error)) exit checks

         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh1 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-9_wp, "linear strat: T*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(Sh1 - Sh0) < 1.0e-9_wp, "linear strat: S*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - 100.0_wp) < 1.0e-8_wp, &
                    "linear strat: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_linear_strat

   ! -----------------------------------------------------------------
   ! 3. near-neutral + small inversion — graceful, finite, conserves
   ! -----------------------------------------------------------------
   subroutine test_unstable(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 5
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rmin, rmax, Th0, Sh0, Th1, Sh1
      integer :: k
      logical :: finite
      checks: block
         call make_eos(eos)
         h_lay = 20.0_wp
         ! near-neutral with a tiny inversion (k=3 slightly warm)
         T_lay = [10.0_wp, 10.0_wp, 10.01_wp, 9.99_wp, 10.0_wp]
         S_lay = S_REF
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         rmin = lin_rho(10.01_wp, S_REF) - 0.2_wp
         rmax = lin_rho(9.99_wp, S_REF) + 0.2_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos)

         finite = .true.
         do k = 1, NZ
            if (ms%h_layer(1, 1, k) /= ms%h_layer(1, 1, k)) finite = .false.   ! NaN
            if (ms%h_layer(1, 1, k) < 0.0_wp) finite = .false.
         end do
         call check(error, finite, "unstable: thicknesses finite + non-negative")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - 100.0_wp) < 1.0e-7_wp, &
                    "unstable: column total conserved")
         if (allocated(error)) exit checks

         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh1 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-8_wp, "unstable: T*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(Sh1 - Sh0) < 1.0e-8_wp, "unstable: S*h conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_unstable

   ! -----------------------------------------------------------------
   ! 4. on-device round-trip — exercise the do concurrent kernel on GPU
   ! -----------------------------------------------------------------
   subroutine test_on_device(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 4
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rmin, rmax, Hsum, Th0, Th1
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      integer :: k, nx, ny
      logical :: finite
      checks: block
         call make_eos(eos)
         do k = 1, NZ
            h_lay(k) = 25.0_wp
            T_lay(k) = 4.0_wp + 12.0_wp*real(k - 1, wp)/real(NZ - 1, wp)
            S_lay(k) = S_REF
         end do
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         nx = grid%nx_total
         ny = grid%ny_total

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         rmin = lin_rho(16.0_wp, S_REF)
         rmax = lin_rho(4.0_wp, S_REF)
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         allocate (bt_eta(nx, ny), source=0.0_wp)
         allocate (bt_H_ref(nx, ny), source=real(NZ, wp)*25.0_wp)

         ! Device round-trip: map the state + vcoord, remap on-device,
         ! pull the results back.  (eos is a flat POD passed by value —
         ! copied into the kernel automatically, no enter_data needed.)
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
         call check(error, abs(Hsum - 100.0_wp) < 1.0e-7_wp, &
                    "device: column total conserved")
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
   ! 5. fully-vanished column — <=1 survivor fast path keeps h unchanged
   ! -----------------------------------------------------------------
   subroutine test_vanished(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 4
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: h_before(NZ)
      real(wp) :: max_dh
      integer :: k
      checks: block
         call make_eos(eos)
         ! Only k=NZ (surface) holds water; the rest are below h_min.
         h_lay = [1.0e-6_wp, 1.0e-6_wp, 1.0e-6_wp, 100.0_wp]
         T_lay = [4.0_wp, 6.0_wp, 8.0_wp, 12.0_wp]
         S_lay = S_REF
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         do k = 1, NZ
            h_before(k) = ms%h_layer(1, 1, k)
         end do

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         vc%zstar_h_min = 1.0e-4_wp   ! floor above the 1e-6 vanished layers
         do k = 0, NZ
            vc%rho_target(k) = RHO0 - 2.0_wp + 4.0_wp*real(k, wp)/real(NZ, wp)
         end do

         call run_remap_host(grid, vc, ms, eos)

         ! Fast path: <=1 survivor -> target_h = remap_h_old (the live h).
         max_dh = 0.0_wp
         do k = 1, NZ
            max_dh = max(max_dh, abs(ms%h_layer(1, 1, k) - h_before(k)))
         end do
         call check(error, max_dh < 1.0e-12_wp, &
                    "vanished: <=1-survivor fast path keeps h_new = h_old")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_vanished

   ! -----------------------------------------------------------------
   ! 6. multi-regrid conservation — the inflation floor must sit STRICTLY
   !    above the remap drain threshold (H_FLOOR), or a collapsed column's
   !    min-thickness layers get their tracer mass zeroed when fed back as
   !    `h_old` on the next regrid.  A single regrid never exposes this
   !    (the inflated layers are h_new, not h_old); the SECOND regrid does.
   !    Uses a near-neutral column that collapses to a thick layer + floor
   !    layers, then regrids again with the same targets.
   ! -----------------------------------------------------------------
   subroutine test_multi_regrid(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      integer, parameter :: NZ = 5
      real(wp) :: h_lay(NZ), T_lay(NZ), S_lay(NZ)
      real(wp) :: rmin, rmax, Th0, Sh0, Th2, Sh2
      integer :: k
      checks: block
         call make_eos(eos)
         h_lay = 20.0_wp
         ! Near-neutral: the targets won't bracket the (flat) column, so it
         ! collapses to one thick layer + min-thickness floor layers.
         T_lay = [10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp]
         S_lay = S_REF
         call setup_column(grid, ms, NZ, h_lay, T_lay, S_lay)
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_RHO
         rmin = lin_rho(10.0_wp, S_REF) - 1.0_wp
         rmax = lin_rho(10.0_wp, S_REF) + 1.0_wp
         do k = 0, NZ
            vc%rho_target(k) = rmin + (rmax - rmin)*real(k, wp)/real(NZ, wp)
         end do

         ! First regrid: collapses to a thick layer + floor layers.
         call run_remap_host(grid, vc, ms, eos)
         ! Second regrid: the collapsed column (with floor layers) is now
         ! h_old.  If the floor <= H_FLOOR, the drain zeroes those layers'
         ! tracer mass here and conservation breaks.
         call run_remap_host(grid, vc, ms, eos)

         Th2 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh2 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th2 - Th0) < 1.0e-9_wp, &
                    "multi-regrid: T*h conserved across two regrids")
         if (allocated(error)) exit checks
         call check(error, abs(Sh2 - Sh0) < 1.0e-9_wp, &
                    "multi-regrid: S*h conserved across two regrids")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - 100.0_wp) < 1.0e-8_wp, &
                    "multi-regrid: column total conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine test_multi_regrid

   ! -----------------------------------------------------------------
   ! 7. thin columns — the step-5 inflation must neither write a negative
   !    thickness nor mint mass when the column cannot hold nz floors.
   ! -----------------------------------------------------------------
   subroutine check_thin_column(error, coord, label, h_lay, T_lay, rho_lo, rho_hi, all_vanished)
      !! Regrid one column with targets linspace(rho_lo, rho_hi) and
      !! assert: every h >= 0; every h either >= the floor or unchanged
      !! from h_old (the too-thin column keeps h_old); sum(h), sum(T*h),
      !! sum(S*h) conserved to round-off.
      !!
      !! `all_vanished`: every source layer is at `H_VANISHED` (a land
      !! column).  Such a column has no live layer, so the remap's I1'
      !! content rule zeroes its tracer content whatever the coordinate
      !! does (`rdb_vl_merge_content`: `hTr = 0` in a column with no live
      !! layer) — the T/S checks are replaced by "h is left bit-for-bit".
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: coord
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: h_lay(:), T_lay(:)
      real(wp), intent(in) :: rho_lo, rho_hi
      logical, intent(in) :: all_vanished
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      real(wp), allocatable :: S_lay(:), h_before(:)
      real(wp) :: H0, Th0, Sh0, floor_eff, hmin
      integer :: k, nz
      logical :: floored
      nz = size(h_lay)
      allocate (S_lay(nz), source=S_REF)
      checks: block
         call make_eos(eos)
         call setup_column(grid, ms, nz, h_lay, T_lay, S_lay)
         h_before = ms%h_layer(1, 1, :)
         H0 = sum(ms%h_layer(1, 1, :))
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call vc%init(grid, nz_ml=nz)
         vc%coord_type = coord
         do k = 0, nz
            vc%rho_target(k) = rho_lo + (rho_hi - rho_lo)*real(k, wp)/real(nz, wp)
         end do
         floor_eff = max(vc%zstar_h_min, 2.0_wp*H_VANISHED)

         call run_remap_host(grid, vc, ms, eos)

         hmin = minval(ms%h_layer(1, 1, :))
         call check(error, hmin >= 0.0_wp, label//": no negative thickness (min h = "// &
                    real_str(hmin)//")")
         if (allocated(error)) exit checks
         floored = .true.
         do k = 1, nz
            if (ms%h_layer(1, 1, k) < floor_eff*(1.0_wp - 1.0e-9_wp) .and. &
                ms%h_layer(1, 1, k) /= h_before(k)) floored = .false.
         end do
         call check(error, floored, label//": every layer >= the floor or untouched")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%h_layer(1, 1, :)) - H0) <= 1.0e-14_wp*max(H0, 1.0_wp), &
                    label//": column total conserved (no minted mass)")
         if (allocated(error)) exit checks
         if (all_vanished) then
            call check(error, all(ms%h_layer(1, 1, :) == h_before), &
                       label//": too-thin column left bit-for-bit")
            exit checks
         end if
         call check(error, abs(sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :)) - Th0) &
                    <= 1.0e-13_wp*max(abs(Th0), 1.0_wp), label//": T*h conserved")
         if (allocated(error)) exit checks
         call check(error, abs(sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :)) - Sh0) &
                    <= 1.0e-13_wp*max(abs(Sh0), 1.0_wp), label//": S*h conserved")
      end block checks
      call vc%destroy()
      call ms%destroy()
   end subroutine check_thin_column

   function real_str(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(es12.4)') x
      s = trim(adjustl(buf))
   end function real_str

   subroutine test_thin_land_rho(error)
      !! A land column (50 x H_VANISHED = 7.5e-3 m), uniform and lighter
      !! than every target: the inversion puts the whole column in the
      !! surface layer, the inflation needs 49 x 3e-4 = 1.47e-2 m — main
      !! debited it from the one survivor: 7.5e-3 - 1.47e-2 = -7.2e-3 m.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(50), T(50)
      h = H_VANISHED
      T = 10.0_wp
      call check_thin_column(error, VCOORD_RHO, "land/rho", h, T, &
                             RHO0 + 1.0_wp, RHO0 + 5.0_wp, .true.)
   end subroutine test_thin_land_rho

   subroutine test_thin_land_hycom(error)
      !! Same land column on HYCOM: the (sigma) floor spreads it into 50
      !! layers of 1.5e-4 <= floor, so `ns == 0` set every layer to 3e-4 —
      !! minting 7.5e-3 m of water (and T, S) per land column per regrid.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(50), T(50)
      h = H_VANISHED
      T = 10.0_wp
      call check_thin_column(error, VCOORD_HYCOM, "land/hycom", h, T, &
                             RHO0 + 1.0_wp, RHO0 + 5.0_wp, .true.)
   end subroutine test_thin_land_hycom

   subroutine test_thin_minted_rho(error)
      !! A thin wet sliver (50 x 2e-4 = 1e-2 m < 50 x 3e-4), stratified so
      !! the inversion spreads it: every layer <= floor => `ns == 0` minted
      !! 5e-3 m, or a partial spread debited a survivor negative.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(50), T(50)
      integer :: k
      h = 2.0e-4_wp
      do k = 1, 50
         T(k) = 4.0_wp + 0.2_wp*real(k - 1, wp)
      end do
      call check_thin_column(error, VCOORD_RHO, "minted/rho", h, T, &
                             lin_rho(T(50), S_REF), lin_rho(T(1), S_REF), .false.)
   end subroutine test_thin_minted_rho

   subroutine test_thin_minted_hycom(error)
      !! As `test_thin_minted_rho` on HYCOM.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(50), T(50)
      integer :: k
      h = 2.0e-4_wp
      do k = 1, 50
         T(k) = 4.0_wp + 0.2_wp*real(k - 1, wp)
      end do
      call check_thin_column(error, VCOORD_HYCOM, "minted/hycom", h, T, &
                             lin_rho(T(50), S_REF), lin_rho(T(1), S_REF), .false.)
   end subroutine test_thin_minted_hycom

   subroutine test_shared_debit_rho(error)
      !! A column just thick enough (50 x 3.2e-4 = 1.6e-2 m >= 50 x 3e-4)
      !! that is a two-layer step (warm 25 over cold 25) with the first
      !! interior target on the step and every other target denser than
      !! the column: the inversion leaves two ~8e-3 m layers and 48 empty
      !! ones.  The 48 x 3e-4 = 1.44e-2 m debit exceeds either survivor —
      !! main drove the thickest to -6.4e-3 m; the fix shares the debit.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(50), T(50), rho_w, rho_c
      h = 3.2e-4_wp
      T(1:25) = 4.0_wp
      T(26:50) = 12.0_wp
      rho_w = lin_rho(12.0_wp, S_REF)
      rho_c = lin_rho(4.0_wp, S_REF)
      ! target(1) = midway between the two layer densities; target(k>=2)
      ! denser than the cold layer => interfaces 2..49 land on the bed.
      call check_thin_column(error, VCOORD_RHO, "shared-debit/rho", h, T, &
                             rho_w - 0.5_wp*(rho_c - rho_w), &
                             rho_w - 0.5_wp*(rho_c - rho_w) + 50.0_wp*(rho_c - rho_w), .false.)
   end subroutine test_shared_debit_rho

end module test_ocean_vcoord_rho
