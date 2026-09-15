!! Unit tests for the KPP non-local (counter-gradient) tracer
!! transport (`vmix_apply_nonlocal_tendencies` after
!! `vmix_apply_kpp_overlay` populates `gamma_t` /
!! `gamma_s`).
!!
!! Non-local γ is non-zero only inside the surface boundary layer
!! under destabilizing flux (B_0 < 0).  Each interface k gets
!!   γ_T(k) = C_s · σ(1-σ)² · (Q_heat / (ρ_0 · cp))
!! with σ = d_face(k) / h_BL.  The applied tendency is the flux
!! divergence at layer k:
!!   hT(k) += dt · (γ_T(k+1) - γ_T(k))
!! (k=1 bed, k=nz surface — ROMS-style).
!!
!! Cases:
!!   * stabilizing flux (Q_heat > 0) → γ everywhere zero.  hT
!!     is bit-identical between a no-γ baseline and a destabilizing-
!!     flux KPP overlay that didn't trigger γ.
!!   * cooling redistributes — Q_heat < 0, a column with no surface
!!     T anomaly.  After one apply the *interior* (k < nz) layers
!!     gain heat at the expense of the surface layer.  Subsurface
!!     hT delta is non-zero and has the same sign as the surface
!!     loss.
!!   * column conservation — total column hT is conserved to round
!!     off (γ_T at bed and surface is zero by construction, so the
!!     summed divergence is identically zero).
!!   * γ vanishes below the BL — for cells where `bl_depth < d_face`,
!!     gamma_t(:, :, k) must stay at zero even under cooling.
module test_ocean_kpp_nonlocal
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_vmix, only: ocean_vmix_t, &
                             vmix_apply_kpp_overlay, &
                             vmix_apply_nonlocal_tendencies
   implicit none
   private

   public :: collect_ocean_kpp_nonlocal_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_kpp_nonlocal_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nonlocal_zero_under_stabilizing_flux", &
                               test_stabilizing_no_op), &
                  new_unittest("nonlocal_cooling_redistributes_T", &
                               test_cooling_redistributes), &
                  new_unittest("nonlocal_column_heat_conserved", &
                               test_column_conserved), &
                  new_unittest("nonlocal_zero_outside_BL", &
                               test_zero_outside_bl) &
                  ]
   end subroutine collect_ocean_kpp_nonlocal_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine setup_state(grid, ms, vmix, ss, sf, nz)
      !! Mixed-layer-over-thermocline IC matching `test_ocean_kpp_-
      !! convective` — strong shear in the top 3 layers, density
      !! step at the thermocline, scalar wind stress.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: nz
      integer, parameter :: NZ_ML = 3
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: U_ML = 0.5_wp
      real(wp), parameter :: S_REF = 35.0_wp
      real(wp), parameter :: T_REF = 18.0_wp
      integer :: k

      ms%nz_ml = nz
      call ms%init(grid)
      call vmix%init(grid, nz_ml=nz)
      call ss%init(grid)
      call sf%init(grid)

      ms%h_layer = H_LAYER
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         if (k > nz - NZ_ML) then
            ms%u_face_x_layer(:, :, k) = U_ML
            ms%rho_layer(:, :, k) = 1027.0_wp
         else
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
         end if
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_REF*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_REF*H_LAYER
      end do
   end subroutine setup_state

   subroutine run_overlay(grid, ms, vmix, ss, sf, dt, do_apply)
      !! Run the KPP overlay (populates γ_T / γ_S) and optionally
      !! apply the non-local tendency on the tracers.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: dt
      logical, intent(in) :: do_apply
      !$acc enter data copyin(ms, vmix, ss, sf)
      call ms%enter_data()
      call vmix%enter_data()
      call sf%enter_data()
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf=sf)
      if (do_apply) then
         call vmix_apply_nonlocal_tendencies(grid, vmix, ms, dt)
      end if
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth, &
      !$acc&            vmix%gamma_t, vmix%gamma_s)
      call sf%exit_data()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
   end subroutine run_overlay

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_stabilizing_no_op(error)
      !! Q_heat > 0 (warming) → B_0 > 0 → γ = 0 everywhere.  hT and
      !! hS after the apply must be bit-identical to the IC.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 100.0_wp
      real(wp), allocatable :: hT_ic(:, :, :)
      real(wp) :: max_gamma, max_dT
      checks: block

         call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
         call setup_state(grid, ms, vmix, ss, sf, NZ)
         call ss%set_wind_stress_const(0.1_wp, 0.0_wp)
         call sf%set_surface_flux_const(500.0_wp, 0.0_wp)    ! W/m² downward → warming

         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call run_overlay(grid, ms, vmix, ss, sf, DT, do_apply=.true.)

         max_gamma = maxval(abs(vmix%gamma_t)) + maxval(abs(vmix%gamma_s))
         max_dT = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_gamma < 1.0e-15_wp, &
                    "stabilizing flux: γ_T or γ_S non-zero")
         if (allocated(error)) exit checks
         call check(error, max_dT < 1.0e-12_wp, &
                    "stabilizing flux: hT changed despite γ = 0")

      end block checks
      deallocate (hT_ic)
      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_stabilizing_no_op

   subroutine test_cooling_redistributes(error)
      !! Q_heat < 0 (cooling) → B_0 < 0 → γ_T < 0 inside the BL
      !! (downward γ negative = upward heat flux toward the surface
      !! to compensate for the lost surface heat that the
      !! `surface_flux_apply` removes).  Run the non-local apply
      !! in *isolation* (no surface flux apply) so the only term
      !! moving heat is γ.  The divergence then warms the surface
      !! layer (`γ(nz+1)=0 - γ(nz)<0 = positive`) and cools the
      !! subsurface (`γ(nz)<γ(nz-1)<0 ⇒ negative`).  Probes the
      !! sign convention of the LMD94 form.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 100.0_wp
      real(wp), allocatable :: hT_ic(:, :, :)
      real(wp) :: min_gamma, dT_surface, dT_subsurface
      integer :: i_probe, j_probe
      checks: block

         call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
         call setup_state(grid, ms, vmix, ss, sf, NZ)
         call ss%set_wind_stress_const(0.1_wp, 0.0_wp)
         call sf%set_surface_flux_const(-500.0_wp, 0.0_wp)     ! cooling

         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call run_overlay(grid, ms, vmix, ss, sf, DT, do_apply=.true.)

         min_gamma = minval(vmix%gamma_t)
         call check(error, min_gamma < 0.0_wp, &
                    "cooling: γ_T did not become negative inside BL")
         if (allocated(error)) exit checks

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         dT_surface = ms%tracers(ms%idx_temperature)%hTr(i_probe, j_probe, NZ) - &
                      hT_ic(i_probe, j_probe, NZ)
         dT_subsurface = ms%tracers(ms%idx_temperature)%hTr(i_probe, j_probe, NZ - 1) - &
                         hT_ic(i_probe, j_probe, NZ - 1)

         ! Non-local in isolation (no surface_flux_apply): the
         ! divergence at the surface layer is `γ(nz+1) - γ(nz) =
         ! 0 - (negative) = positive` → surface hT INCREASES.  The
         ! subsurface sees `γ(nz) - γ(nz-1) = larger_negative -
         ! smaller_negative = negative` → subsurface hT DECREASES.
         call check(error, dT_surface > 0.0_wp, &
                    "cooling-only γ apply: surface hT did not increase")
         if (allocated(error)) exit checks
         call check(error, dT_subsurface < 0.0_wp, &
                    "cooling-only γ apply: subsurface hT did not decrease")

      end block checks
      deallocate (hT_ic)
      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_cooling_redistributes

   subroutine test_column_conserved(error)
      !! γ_T at the bed (k=1) and free surface (k=nz+1) is identically
      !! zero, so `Σ_k (γ(k+1) - γ(k)) = 0` per column.  Total column
      !! hT is therefore conserved by the non-local apply to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 100.0_wp
      real(wp), allocatable :: col_before(:, :), col_after(:, :)
      real(wp) :: max_drift

      call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
      call setup_state(grid, ms, vmix, ss, sf, NZ)
      call ss%set_wind_stress_const(0.1_wp, 0.0_wp)
      call sf%set_surface_flux_const(-500.0_wp, 0.0_wp)

      allocate (col_before(grid%nx_total, grid%ny_total))
      allocate (col_after(grid%nx_total, grid%ny_total))
      col_before = sum(ms%tracers(ms%idx_temperature)%hTr, dim=3)

      call run_overlay(grid, ms, vmix, ss, sf, DT, do_apply=.true.)

      col_after = sum(ms%tracers(ms%idx_temperature)%hTr, dim=3)
      max_drift = maxval(abs(col_after - col_before))

      call check(error, max_drift < 1.0e-9_wp, &
                 "non-local apply did not conserve column total hT")

      deallocate (col_before, col_after)
      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_column_conserved

   subroutine test_zero_outside_bl(error)
      !! For every cell, γ_T at interface k must be zero when the
      !! interface lies BELOW the BL (i.e. `d_face(k) >= h_BL`).
      !! Walks each column, accumulates `d_face`, and confirms the
      !! γ values at deeper interfaces are exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 100.0_wp
      real(wp) :: d_face, h_b, gamma_below
      integer :: i, j, k, nx, ny

      call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
      call setup_state(grid, ms, vmix, ss, sf, NZ)
      call ss%set_wind_stress_const(0.1_wp, 0.0_wp)
      call sf%set_surface_flux_const(-500.0_wp, 0.0_wp)

      call run_overlay(grid, ms, vmix, ss, sf, DT, do_apply=.false.)

      nx = grid%nx_total
      ny = grid%ny_total
      gamma_below = 0.0_wp
      do j = grid%nghost + 1, ny - grid%nghost
         do i = grid%nghost + 1, nx - grid%nghost
            h_b = vmix%bl_depth(i, j)
            d_face = 0.0_wp
            ! Walk from surface (k=nz) down to bed; interface k
            ! sits at the bottom of layer k.  γ is only populated
            ! at k = 2..nz (the kernel skips k=1).  So we check
            ! interfaces from nz down to 2.
            do k = NZ, 2, -1
               d_face = d_face + ms%h_layer(i, j, k)
               if (d_face >= h_b) then
                  gamma_below = max(gamma_below, abs(vmix%gamma_t(i, j, k)))
                  gamma_below = max(gamma_below, abs(vmix%gamma_s(i, j, k)))
               end if
            end do
         end do
      end do

      call check(error, gamma_below < 1.0e-15_wp, &
                 "γ non-zero at interface below the BL")

      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_zero_outside_bl

end module test_ocean_kpp_nonlocal
