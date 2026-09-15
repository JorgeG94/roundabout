!! Unit tests for the shear-driven KPP boundary-layer overlay
!! (`vmix_apply_kpp_overlay` in `rdb_ocean_vmix`).
!!
!! KPP runs on top of an interior closure (PP81 in our case) and
!! overlays a shape-function-shaped kv inside the surface boundary
!! layer.  Phase 1 (this branch) handles the shear-driven case
!! only — surface heat-flux coupling, V_t², and the non-local
!! transport term land in a follow-up branch.
!!
!! Cases:
!!   * No wind stress → u_star = 0 → overlay is a no-op.  kv stays
!!     at whatever the interior closure produced.
!!   * Strongly stratified column + wind → finite h_b that's clearly
!!     less than the total column depth and clearly greater than the
!!     top layer thickness.  The shape function dumps a recognisable
!!     bump into the upper interfaces.
!!   * Unstratified column + wind → no Ri-bulk crossing → h_b = full
!!     column depth.  The entire column gets the cubic shape overlay.
!!   * Shape-function check at a known h_b: probe several interfaces
!!     and confirm `kv_kpp = h_b · u_star · σ(1-σ)²`.
module test_ocean_kpp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, &
                             vmix_compute_pp81, &
                             vmix_apply_kpp_overlay
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   implicit none
   private

   public :: collect_ocean_kpp_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_kpp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("kpp_no_wind_is_no_op", test_no_wind_noop), &
                  new_unittest("kpp_stratified_bl_depth_in_column", test_strat_bl_depth), &
                  new_unittest("kpp_unstratified_bl_fills_column", test_unstrat_fills_column), &
                  new_unittest("kpp_shape_function_recovers_analytic", test_shape_function), &
                  new_unittest("kpp_vt2_deepens_bl", test_vt2_deepens_bl), &
                  new_unittest("kpp_vt2_no_forcing_invariant", test_vt2_no_forcing_invariant), &
                  new_unittest("kpp_vt2_monotonic_in_coefficient", &
                               test_vt2_monotonic_in_coefficient), &
                  new_unittest("kpp_half_domain_cooling_bl_depth_differs", &
                               test_half_domain_cooling) &
                  ]
   end subroutine collect_ocean_kpp_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine run_kpp(grid, ms, vmix, ss)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      ! Zero-flux surface-flux slot: sf is now REQUIRED by the overlay
      ! (A7 per-column B_0); zero fields = the old shear-only behaviour.
      type(ocean_surface_flux_t) :: sf_zero
      call sf_zero%init(grid)
      !$acc enter data copyin(ms, vmix, ss)
      call ms%enter_data()
      call vmix%enter_data()
      !$acc enter data copyin(sf_zero)
      call sf_zero%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf_zero)
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth)
      call sf_zero%exit_data()
      !$acc exit data delete(sf_zero)
      call sf_zero%destroy()
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss)
   end subroutine run_kpp

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_no_wind_noop(error)
      !! Zero wind stress → u_star = 0 → KPP overlay must not write
      !! anything (the `max()` against the interior closure stays at
      !! the interior value).  Verify by comparing a pure-PP81 run
      !! vs a PP81+KPP run with zero stress.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_pp81, vmix_kpp
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NZ = 5
      real(wp), parameter :: H_LAYER = 5.0_wp
      real(wp) :: max_kv_diff, max_kt_diff
      integer :: k
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix_pp81%init(grid, nz_ml=NZ)
         call vmix_kpp%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
         end do

         ! PP81-only reference run
         !$acc enter data copyin(ms, vmix_pp81)
         call ms%enter_data()
         call vmix_pp81%enter_data()
         call vmix_compute_pp81(grid, vmix_pp81, ms)
         !$acc update self(vmix_pp81%kv, vmix_pp81%kt)
         call vmix_pp81%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vmix_pp81)

         ! PP81 + KPP overlay
         call run_kpp(grid, ms, vmix_kpp, ss)

         max_kv_diff = maxval(abs(vmix_kpp%kv - vmix_pp81%kv))
         max_kt_diff = maxval(abs(vmix_kpp%kt - vmix_pp81%kt))

         call check(error, max_kv_diff < 1.0e-14_wp, &
                    "no-wind KPP modified kv vs PP81")
         if (allocated(error)) exit checks
         call check(error, max_kt_diff < 1.0e-14_wp, &
                    "no-wind KPP modified kt vs PP81")

      end block checks
      call ss%destroy(); call vmix_kpp%destroy(); call vmix_pp81%destroy(); call ms%destroy()
   end subroutine test_no_wind_noop

   subroutine test_strat_bl_depth(error)
      !! Mixed-layer-over-thermocline setup: top 3 layers have
      !! uniform velocity + uniform ρ (well-mixed), bottom 5 layers
      !! are quiescent + strongly stratified.  Bulk Ri stays near
      !! zero across the mixed layer (no Δρ, no ΔV²) and jumps
      !! sharply at the thermocline base.  The BL should land
      !! somewhere in or just below the mixed layer — not collapse
      !! to a sub-layer thickness and not run all the way to the bed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NZ = 8
      integer, parameter :: NZ_ML = 3              ! top 3 layers form the mixed layer
      real(wp), parameter :: H_LAYER = 25.0_wp     ! 8 * 25 m = 200 m
      real(wp), parameter :: TAU = 0.1_wp
      real(wp), parameter :: U_ML = 0.5_wp         ! mixed-layer velocity
      real(wp) :: total_H, h_b_probe, h_ml_top
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call ss%set_wind_stress_const(TAU, 0.0_wp)

         ms%h_layer = H_LAYER
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            if (k > NZ - NZ_ML) then
               ! Mixed layer: uniform velocity, uniform ρ.
               ms%u_face_x_layer(:, :, k) = U_ML
               ms%rho_layer(:, :, k) = 1027.0_wp
            else
               ! Thermocline + abyss: at rest, stratified.
               ms%u_face_x_layer(:, :, k) = 0.0_wp
               ms%rho_layer(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
            end if
         end do

         call run_kpp(grid, ms, vmix, ss)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         h_b_probe = vmix%bl_depth(i_probe, j_probe)
         total_H = NZ*H_LAYER
         h_ml_top = real(NZ_ML, wp)*H_LAYER         ! mixed-layer thickness from surface

         ! BL must extend below the topmost layer (the bulk-Ri sweep is
         ! seeded from the top layer centre, so anything shallower
         ! than H_LAYER means it terminated at the very first
         ! comparison — a stencil bug).  Linear bulk-Ri interpolation
         ! lands h_b somewhere inside the mixed layer rather than at
         ! the thermocline base — same known behaviour as MOM6's
         ! CVMix_KPP — so we don't pin it to `h_ml_top` exactly.
         call check(error, h_b_probe > H_LAYER, &
                    "stratified BL collapsed to a sub-layer thickness")
         if (allocated(error)) exit checks
         call check(error, h_b_probe < total_H, &
                    "stratified BL extended to full column (Ri never crossed?)")
         if (allocated(error)) exit checks
         ! KPP overlay at the mid-BL interfaces should lift kv above
         ! the PP81 quiescent baseline (ν_bg).  We don't compare to the
         ! shear-driven PP81 cap because the mixed layer has no shear.
         call check(error, &
                    maxval(vmix%kv(i_probe, j_probe, 2:NZ)) > 5.0_wp*vmix%pp81_nu_bg, &
                    "KPP overlay didn't lift kv above background")

      end block checks
      call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_strat_bl_depth

   subroutine test_unstrat_fills_column(error)
      !! Unstratified column (ρ uniform) + wind → bulk Ri stays at
      !! 0 throughout (no buoyancy gradient).  Crossing never
      !! happens, so h_b falls through to the total column depth.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 20.0_wp
      real(wp), parameter :: TAU = 0.05_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      real(wp) :: total_H, h_b_probe
      integer :: i_probe, j_probe

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZ)
      call ss%init(grid)
      call ss%set_wind_stress_const(TAU, 0.0_wp)

      ms%h_layer = H_LAYER
      ms%rho_layer = 1030.0_wp    ! uniform
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      call run_kpp(grid, ms, vmix, ss)

      i_probe = grid%nx_total/2
      j_probe = grid%ny_total/2
      h_b_probe = vmix%bl_depth(i_probe, j_probe)
      total_H = NZ*H_LAYER

      call check(error, abs(h_b_probe - total_H) < TOL, &
                 "unstratified column: BL depth didn't fall through to total H")

      call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_unstrat_fills_column

   subroutine test_shape_function(error)
      !! Confirm the cubic shape function `G(σ) = σ(1-σ)²` shows up
      !! in `kv` at the right interface depths.  We use the
      !! unstratified case so h_b = total H is exactly known.
      !! Probe a few interfaces, compute expected kv, compare.
      !!
      !! Note: the closure uses `max(kv_kpp, kv_interior)`, so we
      !! must place the probe where kv_kpp > kv_interior.  PP81 on
      !! a uniform-ρ column with no shear ⇒ kv_pp81 = ν_bg + ν_0
      !! (Ri = 0 ⇒ full mixing).  KPP at σ near the peak gives
      !! 4/27·h_b·u* ≈ 0.148·h_b·u*.  For 80 m column and τ = 0.05:
      !! u* ≈ 0.00695 m/s, kv_kpp_peak ≈ 0.082 m²/s — comfortably
      !! larger than the PP81 0.0101 m²/s.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 20.0_wp
      real(wp), parameter :: TAU = 0.05_wp, RHO0 = 1035.0_wp
      real(wp) :: total_H, u_star, kv_expected, sigma, g_shape, d_face_k
      real(wp) :: kv_obs, kv_pp81_max
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call ss%set_wind_stress_const(TAU, 0.0_wp)

         ms%h_layer = H_LAYER
         ms%rho_layer = 1030.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call run_kpp(grid, ms, vmix, ss)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         total_H = NZ*H_LAYER
         u_star = sqrt(TAU/RHO0)

         ! Reference interior value PP81 produces under uniform ρ + no shear:
         kv_pp81_max = vmix%pp81_nu_bg + vmix%pp81_nu0

         ! Probe interior interfaces (k = 2..NZ).  Interface k sits at
         ! depth = (NZ - k + 1) * H_LAYER from the surface.
         do k = 2, NZ
            d_face_k = real(NZ - k + 1, wp)*H_LAYER
            sigma = d_face_k/total_H
            g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
            kv_expected = max(kv_pp81_max, total_H*u_star*g_shape)
            kv_obs = vmix%kv(i_probe, j_probe, k)
            call check(error, abs(kv_obs - kv_expected) < 1.0e-9_wp, &
                       "KPP shape function: kv didn't match analytic G(σ)")
            if (allocated(error)) exit checks
         end do

      end block checks
      call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_shape_function

   ! -----------------------------------------------------------------
   ! V_t² (LMD94 eq 23) — unresolved-turbulence term in bulk-Ri.
   ! -----------------------------------------------------------------

   subroutine run_kpp_with_flux(grid, ms, vmix, ss, sf)
      !! Like `run_kpp` but threads a surface-flux slot through so the
      !! V_t² code path sees a non-zero B_0.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
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
   end subroutine run_kpp_with_flux

   subroutine setup_stratified_column(ms, NZ)
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: NZ
      integer :: k
      ms%h_layer = 5.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      ! Linear stratification: surface (k = NZ) lightest, bed (k = 1)
      ! densest.  Δρ = 1 kg/m^3 across 5 layers (~0.2 kg/m^3 per
      ! layer) — strong enough to clearly bracket the bulk-Ri
      ! crossing.
      do k = 1, NZ
         ms%rho_layer(:, :, k) = 1030.0_wp - 0.2_wp*real(k - 1, wp)
      end do
   end subroutine setup_stratified_column

   subroutine test_vt2_deepens_bl(error)
      !! V_t² (LMD94 eq 23) adds unresolved turbulence to the bulk-Ri
      !! denominator → smaller Ri_bulk at the same depth → the
      !! crossing of `ri_crit` shifts deeper.  Quantitative
      !! discriminator: same wind + stratification, two runs differing
      !! only in `c_vt2`, the V_t²-on run must produce a strictly
      !! greater `bl_depth` at the probe.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_off, vmix_on
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 6
      real(wp), parameter :: TAU = 0.05_wp
      real(wp), parameter :: Q_HEAT = -200.0_wp  ! W/m² — cooling
      real(wp) :: h_b_off, h_b_on
      integer :: i_probe, j_probe

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid)
      call ss%set_wind_stress_const(TAU, 0.0_wp)
      call sf%init(grid)
      call sf%set_surface_flux_const(Q_HEAT, 0.0_wp)

      checks: block
         ! ---- V_t² off (c_vt2 = 0) ----
         call vmix_off%init(grid, nz_ml=NZ)
         vmix_off%c_vt2 = 0.0_wp
         call setup_stratified_column(ms, NZ)
         call run_kpp_with_flux(grid, ms, vmix_off, ss, sf)

         ! ---- V_t² on (c_vt2 = 1.8) ----
         call vmix_on%init(grid, nz_ml=NZ)
         vmix_on%c_vt2 = 1.8_wp
         call setup_stratified_column(ms, NZ)
         call run_kpp_with_flux(grid, ms, vmix_on, ss, sf)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         h_b_off = vmix_off%bl_depth(i_probe, j_probe)
         h_b_on = vmix_on%bl_depth(i_probe, j_probe)

         ! V_t²-on must produce a deeper or equal BL.  Equality is
         ! permitted because the discrete sweep walks layer centres, so
         ! V_t² of a few cm² s⁻² may not be enough to shift the crossing
         ! into a new layer.  The strict inequality is the usual case
         ! under the cooling + wind setup here.
         call check(error, h_b_on >= h_b_off, &
                    "V_t²-on BL must be at least as deep as V_t²-off")
         if (allocated(error)) exit checks
         call check(error, h_b_on > h_b_off + 1.0e-6_wp, &
                    "V_t²-on did not visibly deepen BL "// &
                    "(check that wind + stratification are non-trivial)")
      end block checks

      call sf%destroy(); call ss%destroy()
      call vmix_on%destroy(); call vmix_off%destroy(); call ms%destroy()
   end subroutine test_vt2_deepens_bl

   subroutine test_vt2_no_forcing_invariant(error)
      !! When there's no wind and no surface buoyancy flux,
      !! u_* = 0, w_* = 0, w_s_col = 0 → V_t² = 0 regardless of
      !! `c_vt2`.  Two runs differing only in `c_vt2` must produce
      !! identical `bl_depth` to round-off.  Catches a regression
      !! where V_t² is computed without a w_s_col guard.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 6
      real(wp), parameter :: TOL = 1.0e-13_wp
      real(wp) :: max_diff
      integer :: i_probe, j_probe

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid)
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      call sf%init(grid)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      call vmix_a%init(grid, nz_ml=NZ)
      vmix_a%c_vt2 = 0.0_wp
      call setup_stratified_column(ms, NZ)
      call run_kpp_with_flux(grid, ms, vmix_a, ss, sf)

      call vmix_b%init(grid, nz_ml=NZ)
      vmix_b%c_vt2 = 10.0_wp  ! exaggerated; should still produce same h_b
      call setup_stratified_column(ms, NZ)
      call run_kpp_with_flux(grid, ms, vmix_b, ss, sf)

      max_diff = maxval(abs(vmix_a%bl_depth - vmix_b%bl_depth))
      i_probe = grid%nx_total/2
      j_probe = grid%ny_total/2

      call check(error, max_diff < TOL, &
                 "V_t² changed BL depth with no wind / no flux "// &
                 "(w_s_col guard missing)")

      call sf%destroy(); call ss%destroy()
      call vmix_b%destroy(); call vmix_a%destroy(); call ms%destroy()
   end subroutine test_vt2_no_forcing_invariant

   subroutine test_vt2_monotonic_in_coefficient(error)
      !! Increasing `c_vt2` strictly deepens (or holds) the BL depth.
      !! V_t² appears positive-definite in the denominator, so the
      !! mapping `c_vt2 → bl_depth` is monotonic.  Three c_vt2
      !! sample points; each step must produce h_b ≥ the previous
      !! sample, with at least one strict increase.  Confirms the
      !! sign + dependence of V_t² rather than just that "something
      !! happens".
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_lo, vmix_md, vmix_hi
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 6
      real(wp), parameter :: TAU = 0.05_wp
      real(wp), parameter :: Q_HEAT = -150.0_wp
      real(wp) :: h_lo, h_md, h_hi
      integer :: i_probe, j_probe

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid)
      call ss%set_wind_stress_const(TAU, 0.0_wp)
      call sf%init(grid)
      call sf%set_surface_flux_const(Q_HEAT, 0.0_wp)

      checks: block
         call vmix_lo%init(grid, nz_ml=NZ)
         vmix_lo%c_vt2 = 0.0_wp
         call setup_stratified_column(ms, NZ)
         call run_kpp_with_flux(grid, ms, vmix_lo, ss, sf)

         call vmix_md%init(grid, nz_ml=NZ)
         vmix_md%c_vt2 = 1.8_wp
         call setup_stratified_column(ms, NZ)
         call run_kpp_with_flux(grid, ms, vmix_md, ss, sf)

         call vmix_hi%init(grid, nz_ml=NZ)
         vmix_hi%c_vt2 = 5.0_wp
         call setup_stratified_column(ms, NZ)
         call run_kpp_with_flux(grid, ms, vmix_hi, ss, sf)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         h_lo = vmix_lo%bl_depth(i_probe, j_probe)
         h_md = vmix_md%bl_depth(i_probe, j_probe)
         h_hi = vmix_hi%bl_depth(i_probe, j_probe)

         call check(error, h_md >= h_lo, &
                    "h_b non-monotonic between c_vt2 = 0 and 1.8")
         if (allocated(error)) exit checks
         call check(error, h_hi >= h_md, &
                    "h_b non-monotonic between c_vt2 = 1.8 and 5.0")
         if (allocated(error)) exit checks
         call check(error, h_hi > h_lo + 1.0e-6_wp, &
                    "h_b did not visibly grow from c_vt2 = 0 to 5.0")
      end block checks

      call sf%destroy(); call ss%destroy()
      call vmix_hi%destroy(); call vmix_md%destroy(); call vmix_lo%destroy()
      call ms%destroy()
   end subroutine test_vt2_monotonic_in_coefficient

   subroutine test_half_domain_cooling(error)
      !! Half-domain cooling: Q_heat < 0 in the west half (x ≤ nx/2),
      !! Q_heat = 0 in the east half.  After TWO consecutive KPP calls
      !! (so the lagged w* from the first call seeds the second), the
      !! BL depth must be strictly deeper in the cooled west half than
      !! in the east half.  This catches scalar-broadcast regressions in
      !! the closure consumers — if BL depth is spatially uniform
      !! despite a spatial Q pattern, the 2D field is being ignored.
      !!
      !! Note: KPP uses a lagged w* (`h_b_lagged`); a single call with
      !! bl_depth=0 IC produces w*=0 everywhere regardless of B_0 (the
      !! product `(-B_0)·h_b_lagged` is zero).  Two calls are needed:
      !! the first seeding bl_depth; the second using it to build w*.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 6
      integer, parameter :: NX_PHYS = 12, NY_PHYS = 8
      real(wp), parameter :: TAU = 0.05_wp      ! uniform wind
      real(wp), parameter :: Q_COOL = -300.0_wp ! W/m² cooling in west half
      real(wp) :: h_b_west, h_b_east
      integer :: i, j_mid, i_west, i_east, nx_tot

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZ)
      call ss%init(grid)
      call ss%set_wind_stress_const(TAU, 0.0_wp)
      call sf%init(grid)
      ! Seed with zero first, then overwrite west half with cooling.
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      ! Cooling only in the west half (i ≤ nx_total/2).
      nx_tot = grid%nx_total
      do i = 1, nx_tot/2
         sf%Q_heat(i, :) = Q_COOL
      end do
      sf%has_heat = .true.

      ! First call: seeds bl_depth so the second call has a non-zero
      ! h_b_lagged to build w* from.
      call setup_stratified_column(ms, NZ)
      call run_kpp_with_flux(grid, ms, vmix, ss, sf)
      ! Second call: uses the seeded bl_depth to produce a non-zero w*
      ! in the cooled west half — should push BL deeper than the east.
      call run_kpp_with_flux(grid, ms, vmix, ss, sf)

      j_mid = grid%ny_total/2
      i_west = max(1, nx_tot/2 - 1)      ! probe inside the cooled region
      i_east = min(nx_tot, nx_tot/2 + 2) ! probe inside the unforced region
      h_b_west = vmix%bl_depth(i_west, j_mid)
      h_b_east = vmix%bl_depth(i_east, j_mid)

      call check(error, h_b_west > h_b_east + 1.0e-6_wp, &
                 "half-domain cooling: west BL not deeper than east BL "// &
                 "(scalar-broadcast regression in KPP?)")

      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_half_domain_cooling
end module test_ocean_kpp
