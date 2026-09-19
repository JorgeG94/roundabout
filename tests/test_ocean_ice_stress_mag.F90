!! Gates for the cell-centred surface-stress magnitude `stress_mag` under
!! sea ice.
!!
!! `stress_mag` is the ONLY source KPP (`rdb_ocean_vmix`) and EPBL
!! (`rdb_ocean_epbl`) take the surface friction velocity from —
!! `u_* = sqrt(|tau| / rho0)`.  The sea-ice coupler
!! (`ice_ocean_stress_flux`, `rdb_ice_ocean_coupler`) OVERWRITES
!! `tau_x`/`tau_y` on the device every outer step with the
!! concentration-weighted blend of the atmospheric-stress snapshot and the
!! EVP ice->ocean drag, so a `stress_mag` that is only ever filled at
!! configure time leaves both boundary-layer schemes mixing on the
!! configure-time WIND — zero, under full ice cover with no wind at all.
!!
!! These gates pin the refresh at the blend:
!!   1. `full_cover_tracks_drag`  — zero wind, full ice cover, non-zero
!!      ice->ocean drag: `stress_mag == |(fxoc, fyoc)|`.  Without the
!!      refresh this reads 0 everywhere (the configure-time wind).
!!   2. `ice_free_wind_bitident`  — ice present but zero concentration:
!!      the blend reproduces the wind exactly, so `stress_mag` must come
!!      back BIT-IDENTICAL to its configure-time value.
!!   3. `partial_cover_consistent` — mixed cover: `stress_mag` matches the
!!      cell-centred magnitude of whatever `tau` pair the blend actually
!!      produced, and is distinct from both the pure-wind and the pure-drag
!!      magnitude (so gate 1 cannot pass by accident on a stale field).
!!
!! GPU `mem:separate` discipline throughout (canonical template:
!! `test_open_boundary_out_closes` in `test_ocean_conservation_salt_heat.F90`):
!! every object AND its scratch companion is mapped before the kernel call,
!! host-set inputs ride in on `copyin`, and every host assertion is preceded
!! by an `!$acc update self` of the COMPONENT arrays (never the aggregate
!! derived type — that would overwrite the host descriptors with device
!! addresses).  The coupler's own per-step blend scratch is module-level and
!! released by `ice_ocean_stress_cleanup`.
module test_ocean_ice_stress_mag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux, ice_ocean_stress_cleanup
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   implicit none
   private

   public :: collect_ocean_ice_stress_mag_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 4, NYP = 4
   real(wp), parameter :: DX = 2000.0_wp

   ! Ice->ocean drag (Pa) the EVP core would have left in `fxoc`/`fyoc`.
   real(wp), parameter :: FXOC = 0.7_wp
   real(wp), parameter :: FYOC = -0.2_wp
   ! Atmospheric stress snapshot (Pa) — gate 2/3 only; gate 1 runs windless.
   real(wp), parameter :: WIND_X = 0.3_wp
   real(wp), parameter :: WIND_Y = 0.1_wp

   real(wp), parameter :: TOL = 8.0_wp*epsilon(1.0_wp)

contains

   subroutine collect_ocean_ice_stress_mag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("full_cover_tracks_drag", test_full_cover_tracks_drag), &
                  new_unittest("ice_free_wind_bitident", test_ice_free_wind_bitident), &
                  new_unittest("partial_cover_consistent", test_partial_cover_consistent) &
                  ]
   end subroutine collect_ocean_ice_stress_mag_tests

   ! =====================================================================
   ! Gate 1: full ice cover, zero wind -> stress_mag == |ice-ocean drag|
   ! =====================================================================

   subroutine test_full_cover_tracks_drag(error)
      !! The discriminating case.  Wind is IDENTICALLY ZERO, so the
      !! configure-time `stress_mag` is 0 everywhere; ice covers every cell
      !! (ghosts included, so `a_u == a_v == 1` at every face and the blend
      !! is the exact identity `tau = fxoc`).  A correct `stress_mag` is
      !! `sqrt(FXOC^2 + FYOC^2)`; the pre-fix code leaves 0, which is
      !! `u_* = 0` — KPP/EPBL wind mixing silently switched off under ice.
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_stress_t) :: stress
      real(wp) :: expected

      expected = sqrt(FXOC*FXOC + FYOC*FYOC)

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      call make_ice(ice, grid)
      call stress%init(grid, nz_ml=1)

      ! Windless: this also seeds `stress_mag = 0` host-side.
      call stress%set_wind_stress_const(0.0_wp, 0.0_wp)
      ice%tau_a_x = 0.0_wp
      ice%tau_a_y = 0.0_wp
      ice%fxoc = FXOC
      ice%fyoc = FYOC
      ! Full cover EVERYWHERE, ghost rows included: `ncat == 1` is the
      ! lumped mode, where any positive `m_ice` on a wet cell gives ci = 1.
      ice%m_ice = 3.0_wp*ICE_RHO_ICE
      ice%m_snow = 0.0_wp

      checks: block
         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(stress)
         call stress%enter_data()
         ! `metrics%wet_T` is read on-device by `ice_cell_concentration_impl`
         ! inside `ice_ocean_stress_flux`, and `make_cartesian_metrics`
         ! already mapped it (parent + arrays) — do NOT re-map here, or the
         ! paired unmap below would leave it at presentcount 1 on-device.

         call ice_ocean_stress_flux(metrics, stress, ice)

         ! Component arrays only — never `update self(stress)`.
         !$acc update self(stress%tau_x, stress%tau_y, stress%stress_mag)

         ! Physical faces / cells only.  The blend gives the two ARRAY-edge
         ! faces (i = 1, i = nx+1) a deliberate half weight — they have only
         ! one abutting cell — so the full-cover identity `tau == fxoc` is a
         ! statement about the physical band, not about the padding.
         call check(error, all(abs(stress%tau_x(NGHOST + 1:NGHOST + NXP + 1, &
                                                NGHOST + 1:NGHOST + NYP) - FXOC) <= TOL), &
                    "full_cover_tracks_drag: blend did not hand the ocean fxoc")
         if (allocated(error)) exit checks
         call check(error, all(abs(stress%stress_mag(NGHOST + 1:NGHOST + NXP, &
                                                     NGHOST + 1:NGHOST + NYP) - expected) &
                               <= TOL), &
                    "full_cover_tracks_drag: stress_mag is stale (KPP/EPBL u_* "// &
                    "still on the configure-time wind)")
      end block checks

      call teardown(ice, stress, metrics)
   end subroutine test_full_cover_tracks_drag

   ! =====================================================================
   ! Gate 2: ice present but zero concentration -> wind, bit-identically
   ! =====================================================================

   subroutine test_ice_free_wind_bitident(error)
      !! With `ci == 0` the blend is `tau = 1*tau_a + 0*fxoc`, which
      !! reproduces the wind snapshot exactly — so re-deriving `stress_mag`
      !! from it must reproduce the configure-time field BIT-FOR-BIT.  This
      !! is the "the refresh cannot move an ice-free answer" gate: the same
      !! three lines, the same operands, the same FP op order.
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_stress_t) :: stress
      real(wp), allocatable :: mag_configure(:, :)

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      call make_ice(ice, grid)
      call stress%init(grid, nz_ml=1)

      call stress%set_wind_stress_const(WIND_X, WIND_Y)
      ice%tau_a_x = WIND_X
      ice%tau_a_y = WIND_Y
      ice%fxoc = FXOC
      ice%fyoc = FYOC
      ice%m_ice = 0.0_wp
      ice%m_snow = 0.0_wp
      allocate (mag_configure, source=stress%stress_mag)

      checks: block
         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(stress)
         call stress%enter_data()

         call ice_ocean_stress_flux(metrics, stress, ice)

         !$acc update self(stress%tau_x, stress%stress_mag)

         call check(error, all(stress%tau_x == WIND_X), &
                    "ice_free_wind_bitident: zero-cover blend perturbed tau_x")
         if (allocated(error)) exit checks
         call check(error, all(stress%stress_mag == mag_configure), &
                    "ice_free_wind_bitident: refresh moved the ice-free stress_mag")
      end block checks

      call teardown(ice, stress, metrics)
   end subroutine test_ice_free_wind_bitident

   ! =====================================================================
   ! Gate 3: partial cover -> stress_mag consistent with the blended tau
   ! =====================================================================

   subroutine test_partial_cover_consistent(error)
      !! Half the physical domain iced.  `stress_mag` must equal the
      !! cell-centred magnitude of the `tau` pair the blend ACTUALLY left
      !! behind, recomputed here on the host, everywhere — including the
      !! partially covered faces where `a_u` is 0.5.  The second pair of
      !! checks pins that an iced interior cell's value is neither the
      !! pure-wind nor the pure-drag magnitude, so gate 1 cannot be
      !! satisfied by a field that merely happens to be constant.
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_stress_t) :: stress
      real(wp) :: mag_host, tau_x_cell, tau_y_cell, worst
      real(wp) :: mag_wind, mag_drag
      integer :: i, j, nx, ny, i_probe, j_probe

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx = grid%nx_total
      ny = grid%ny_total
      call make_cartesian_metrics(metrics, grid)
      call make_ice(ice, grid)
      call stress%init(grid, nz_ml=1)

      call stress%set_wind_stress_const(WIND_X, WIND_Y)
      ice%tau_a_x = WIND_X
      ice%tau_a_y = WIND_Y
      ice%fxoc = FXOC
      ice%fyoc = FYOC
      ice%m_snow = 0.0_wp
      ice%m_ice = 0.0_wp
      ! Ice on the western half of the physical domain only.
      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP/2
            ice%m_ice(i, j, 1) = 3.0_wp*ICE_RHO_ICE
         end do
      end do
      ! Probe the western-most iced column: its two v-faces are fully
      ! ice-mediated (both abutting cells iced) while its western u-face
      ! straddles the ice edge (a_u = 0.5), so the cell-centred |tau| lands
      ! strictly BETWEEN the pure-wind and the pure-drag magnitude — a value
      ! no stale field can produce.
      i_probe = NGHOST + 1
      j_probe = NGHOST + 2

      mag_wind = sqrt(WIND_X*WIND_X + WIND_Y*WIND_Y)
      mag_drag = sqrt(FXOC*FXOC + FYOC*FYOC)

      checks: block
         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(stress)
         call stress%enter_data()

         call ice_ocean_stress_flux(metrics, stress, ice)

         !$acc update self(stress%tau_x, stress%tau_y, stress%stress_mag)

         worst = 0.0_wp
         do j = 1, ny
            do i = 1, nx
               tau_x_cell = 0.5_wp*(stress%tau_x(i, j) + stress%tau_x(i + 1, j))
               tau_y_cell = 0.5_wp*(stress%tau_y(i, j) + stress%tau_y(i, j + 1))
               mag_host = sqrt(tau_x_cell*tau_x_cell + tau_y_cell*tau_y_cell)
               worst = max(worst, abs(stress%stress_mag(i, j) - mag_host))
            end do
         end do
         call check(error, worst <= TOL, &
                    "partial_cover_consistent: stress_mag does not match the |tau| "// &
                    "of the blend it was derived from")
         if (allocated(error)) exit checks

         call check(error, abs(stress%stress_mag(i_probe, j_probe) - mag_wind) > 1.0e-3_wp, &
                    "partial_cover_consistent: iced cell still reads the wind magnitude")
         if (allocated(error)) exit checks
         call check(error, stress%stress_mag(i_probe, j_probe) > mag_wind, &
                    "partial_cover_consistent: ice mediation did not raise |tau| "// &
                    "towards the drag")
         if (allocated(error)) exit checks
         call check(error, stress%stress_mag(i_probe, j_probe) <= mag_drag + TOL, &
                    "partial_cover_consistent: blended |tau| exceeds the pure-drag bound")
      end block checks

      call teardown(ice, stress, metrics)
   end subroutine test_partial_cover_consistent

   ! =====================================================================
   ! Shared harness
   ! =====================================================================

   subroutine make_ice(ice, grid)
      !! Minimal lumped (`ncat == 1`) sea-ice slot: the coupler only ever
      !! reads `part_size`/`m_ice`/`m_snow` (through
      !! `ice_cell_concentration_impl`) plus `tau_a_*`/`fxoc`/`fyoc`, and
      !! `dynamics` stays `.false.` so no EVP workspace is allocated or
      !! mapped.
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(hgrid_t), intent(in) :: grid
      ice%enable = .true.
      ice%ncat = 1
      ice%nk_ice = 1
      call ice%init(grid)
   end subroutine make_ice

   subroutine teardown(ice, stress, metrics)
      !! Unmap in the reverse of the map order, then release the coupler's
      !! module-level blend scratch (it is sized lazily per grid and lives
      !! across calls — the production teardown does the same next to
      !! `ice_evp_cleanup`).  `metrics` was mapped ONCE by
      !! `make_cartesian_metrics`; `destroy_cartesian_metrics` is its one
      !! matching unmap.
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_stress_t), intent(inout) :: stress
      type(ocean_metrics_t), intent(inout) :: metrics
      !$acc exit data delete(stress)
      call stress%exit_data()
      !$acc exit data delete(ice)
      call ice%exit_data()
      call ice_ocean_stress_cleanup()
      call ice%destroy()
      call stress%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown

end module test_ocean_ice_stress_mag
