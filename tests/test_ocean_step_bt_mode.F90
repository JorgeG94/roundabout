!! STEPPED-TERRAIN BAROTROPIC NON-GROWTH — the gate for finding B of the
!! vertical-coordinate stability matrix (`renorm_consistent_flux`).
!!
!! ## What is guarded
!!
!! Under `pred_corr` the slow continuity renormalises its per-layer mass
!! fluxes onto the barotropic substep's time-mean transport `uhbt`
!! (`renormalise_zonal_flux_to_uhbt`), so that the layer free surface
!! `Σh − H` ends the step equal to the substep's `η_end`.  The Newton solve
!! re-picks each layer's upwind donor at the corrected velocity but kept
!! the OLD donor's `u0·h_old` in its flux model, `flux0 + du·h_new`, which
!! jumps by `u0·(h_new − h_old)·w` where the donor flips.  Over a
!! bathymetric step (sigma layers of different thickness on the two
!! sides, both PPM edges flattened by the limiter) a target of the
!! opposite sign to the layer velocity lands in that jump whenever
!! `|uhbt| < |u0|·Δh·w` — Newton has no root, cycles, and returns a
!! layer transport of the WRONG SIGN.  The layer `η` then differs from
!! `η_end` by O(η) in the two step columns, every such step.
!!
!! Under `pred_corr` that is frequent: the layer velocity entering the
!! corrector continuity is the END-of-step barotropic velocity and the
!! target is the TIME-MEAN transport; for a grid-scale barotropic gravity
!! mode (`ω·dt ≈ 50`) the two are uncorrelated in sign and the mean is
!! ~`sinc(ω·dt/2) ≈ 2 %` of the end value — inside the gap about half
!! the time.  The η kick pumps the mode, which nothing else damps under
!! `pred_corr` at `bebt = 0` (the lateral viscosity sees only the time
!! mean `u_av`), so it grows exponentially: measured 8 % per outer step
!! in amplitude at rx0 = 0.1, 21 % at rx0 = 0.3, on ONE layer, at f = 0,
!! with no viscosity and no remap.  In the 48 × 6 × 15 matrix cell
!! `vcmv_rx0_010_sigma` it is the domain-wide 2Δx barotropic mode that
!! drives a layer to −255 m at step 2087.
!!
!! ## Why this shape of test
!!
!! Flat-bed control is exactly neutral, so the property is energetic: the
!! basin is closed and unforced, so `KE + PE` may not grow.  The seed is
!! a y-uniform barotropic 2Δx velocity (the growing mode's own shape), at
!! 1e-6 m/s so the whole run is linear.  Two assertions:
!!
!!   * `η` CONSISTENCY after every step, `max|Σ_k h − H − η_end|`, the
!!     invariant the defect breaks directly (measured on this case: knob
!!     off 7.2e-4 m, knob on 2.0e-13 m, round-off);
!!   * `KE + PE` ratio after 150 outer steps below `ENERGY_GROWTH_BAR`.
module test_ocean_step_bt_mode
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_LITE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_step_bt_mode_tests

   ! The failing matrix cell's horizontal geometry, time step and rotation
   ! (tests/regression/vcoord_matrix.py): 2 km, dt = 600 s, f-plane at 75 S,
   ! rx0 = 0.1 step (825 m | 675 m) at mid-channel, and its viscous leg's
   ! Laplacian (nu_h = 160 m2/s, stress-divergence form).  Two cells across
   ! the channel is enough: the growing mode is y-uniform.
   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 48
   integer, parameter :: NYP = 2
   integer, parameter :: NZ = 2
   integer, parameter :: N_INNER = 59            !! the cell's auto_n_inner
   integer, parameter :: N_STEPS = 150
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: H_BAR = 750.0_wp
   real(wp), parameter :: RX0 = 0.1_wp
   real(wp), parameter :: F0 = -1.409e-4_wp
   real(wp), parameter :: DT = 600.0_wp
   real(wp), parameter :: NU_H = 160.0_wp
   real(wp), parameter :: U_SEED = 1.0e-6_wp
   real(wp), parameter :: GRAV = 9.80665_wp

   real(wp), parameter :: ENERGY_GROWTH_BAR = 4.0_wp
      !! Bar on `(KE+PE)_end / (KE+PE)_0` after `N_STEPS`.  Measured on this
      !! case (gfortran host build): knob OFF 6.7e+02 and still exponential;
      !! knob ON 1.08.  A neutral scheme sloshes KE <-> PE: the same seed on
      !! a 6-row, 1-layer version of this channel wanders 0.86-1.8 over 600
      !! steps with the knob on.  The bar sits 2.2x above that envelope and
      !! 170x below the defect at this length.  Never widen it to make this
      !! pass.
   real(wp), parameter :: ETA_CONSISTENCY_BAR = 1.0e-10_wp
      !! Bar on `max|Σh − H − η_end|` (m).  Round-off at H = 825 m is
      !! ~1e-13; the defect is O(η) ~ 1e-5.

contains

   subroutine collect_ocean_step_bt_mode_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("pred_corr_step_channel_no_bt_growth", test_step_channel) &
                  ]
   end subroutine collect_ocean_step_bt_mode_tests

   subroutine test_step_channel(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc
      integer :: i, j, k, i0, i1, j0, j1, step
      real(wp) :: hcol, e0, e1, eta, eta_err, ratio
      logical :: finite
      character(len=320) :: msg

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      ct%renorm_consistent_flux = .true.
      cor%f_0 = F0
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_FV_LITE
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      ! bebt = 0 (pure forward-backward, the pre-2026-09-22 default): the
      ! guard must see the renormaliser alone.  MOM6's default BEBT = 0.1
      ! damps the pumped 2Δx mode by itself and would let the historical
      ! discontinuous flux model pass this gate.
      dyn%bt_work%bebt = 0.0_wp
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_LAGRANGIAN
      dyn%split_scheme = SPLIT_SCHEME_PRED_CORR
      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = NU_H
      hv%stress_tensor = .true.
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

      i0 = NGHOST + 1; i1 = NGHOST + NXP
      j0 = NGHOST + 1; j1 = NGHOST + NYP
      ! Uniform density, sigma-like equal layers, a single step face at
      ! mid-channel (ghosts included, so the edge stencils see the step too).
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            hcol = H_BAR*(1.0_wp + RX0)
            if (i > NGHOST + NXP/2) hcol = H_BAR*(1.0_wp - RX0)
            do k = 1, NZ
               ms%h_layer(i, j, k) = hcol/real(NZ, wp)
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            ms%rho_layer(i, j, :) = eos%rho0
            dyn%bt_work%bt_H_ref(i, j) = hcol
         end do
      end do
      ! Barotropic, y-uniform 2Δx seed on the interior faces; walls closed.
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do j = j0, j1
         do i = i0 + 1, i1
            ms%u_face_x_layer(i, j, :) = U_SEED*real(1 - 2*mod(i, 2), wp)
         end do
      end do
      e0 = channel_energy(ms, dyn%bt_work%bt_H_ref, i0, i1, j0, j1)

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      eta_err = 0.0_wp
      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         !$acc update self(ms%h_layer, dyn%bt_work%bt_eta_end)
         do j = j0, j1
            do i = i0, i1
               eta = sum(ms%h_layer(i, j, 1:NZ)) - dyn%bt_work%bt_H_ref(i, j)
               eta_err = max(eta_err, abs(eta - dyn%bt_work%bt_eta_end(i, j)))
            end do
         end do
      end do

      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      finite = .true.
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) finite = .false.
            end do
         end do
      end do
      e1 = channel_energy(ms, dyn%bt_work%bt_H_ref, i0, i1, j0, j1)
      ratio = huge(1.0_wp)
      if (finite .and. e0 > 0.0_wp .and. ieee_is_finite(e1)) ratio = e1/e0

      call vc%exit_data()
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data()
      call ss%exit_data(); call bd%exit_data(); call hv%exit_data()
      call pgf%exit_data(); call cor%exit_data(); call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()

      call check(error, finite, "stepped channel: went non-finite (the barotropic mode blew up)")
      if (allocated(error)) return
      write (msg, '("stepped channel: max|sum(h) - H - eta_end| = ", es11.3, &
            &" m; the layer free surface left the barotropic one — suspect the", &
            &" uhbt renormalisation at the step face (renorm_consistent_flux)")') eta_err
      call check(error, eta_err <= ETA_CONSISTENCY_BAR, trim(msg))
      if (allocated(error)) return
      write (msg, '("stepped channel: KE+PE ratio = ", es11.3, " after ", i0, &
            &" steps; closed + unforced, so growth is manufactured energy")') ratio, N_STEPS
      call check(error, ratio < ENERGY_GROWTH_BAR, trim(msg))
   end subroutine test_step_channel

   pure function channel_energy(ms, h_ref, i0, i1, j0, j1) result(energy)
      !! `Σ_k h·(u² + v²)/2 + g·η²/2` over the interior, per unit density and
      !! area, `η = Σ_k h − H` — a monotone energy diagnostic (face
      !! velocities squared in place), which is all a non-growth assertion
      !! needs.
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: h_ref(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: energy
      integer :: i, j, k
      real(wp) :: eta
      energy = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            eta = -h_ref(i, j)
            do k = 1, NZ
               eta = eta + ms%h_layer(i, j, k)
               energy = energy + 0.5_wp*ms%h_layer(i, j, k)* &
                        (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
            end do
            energy = energy + 0.5_wp*GRAV*eta*eta
         end do
      end do
   end function channel_energy
end module test_ocean_step_bt_mode
