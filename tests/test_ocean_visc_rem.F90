!! Unit tests for the `visc_rem` (viscous remnant) PRODUCER — PR-19.
!!
!! `bt_work%visc_rem_u/v` — the per-face, per-layer γ_k ≡
!! (1/Δt)·∂u_k^{n+1}/∂Ā, the sensitivity of the post-friction layer
!! velocity to a uniform barotropic acceleration — is now FILLED by
!! `vdiff_apply_momentum` (optional `visc_rem_u`/`visc_rem_v` dummies),
!! which solves the SAME already-factorized momentum tridiagonal a
!! second time with RHS ≡ 1 (Roundabout's rows are pre-normalized by h_k,
!! unlike MOM6's un-normalized `h_u(k)` RHS convention — see
!! `rdb_ocean_vdiff.F90`'s `diffuse_velocity_columns_impl` docstring).
!!
!! Convention reminder (Roundabout is bottom-up): k=1 is the BED (γ →
!! smallest), k=nz is the SURFACE (γ → 1) — the MIRROR of MOM6's
!! top-down k=1/k=nz roles.
!!
!! Cases (§9 of PLAN_PR19_visc_rem.md):
!!   1. no_drag_remnant_is_exactly_one — the row-sum invariant: a
!!      no-flux-top/no-flux-bottom viscous operator cannot remove a
!!      uniform acceleration, so γ ≡ 1 to round-off for ANY interior
!!      viscosity/thickness/dt as long as the operator carries no
!!      drag.
!!   2. strong_drag_remnant_matches_analytic_profile — with strong
!!      bottom drag (near no-slip), γ(z) ≈ 1 − exp(−(z+D)/√(νΔt))
!!      in the interior (excludes the near-bed band where the
!!      finite-λ discretization departs from the continuum no-slip
!!      limit, and the top layer where the finite-D free-surface BC
!!      matters).
!!   3. remnant_bounded_and_monotone — the discrete maximum principle:
!!      0 < γ ≤ 1, non-decreasing bed→surface, no NaN.
!!   4. remnant_absent_is_bit_identical — the default-off contract:
!!      the producer is a pure add-on and cannot perturb the momentum
!!      solve it shares a kernel with.
!!   5. remnant_consistent_with_momentum_response — γ IS (not merely
!!      approximates) (1/Δt)·∂u^{n+1}/∂Ā; verified by a finite
!!      difference against the EXACT same linear operator (no
!!      discretization error since the operator is linear).
!!   6. producer_then_corrector_biases_against_bbl — the end-to-end
!!      wire: `vdiff_apply_momentum` (the real producer) fills
!!      `bt_work%visc_rem_u/v`, and `apply_bt_correction` (the real
!!      consumer) reads that SAME field and visibly biases the
!!      BT-corrector's Δu against a frictional bottom boundary layer
!!      while preserving the column's depth-mean Δu — the trap-#1
!!      defence (a producer and consumer that both pass while the
!!      wire between them is dead).
module test_ocean_visc_rem
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_momentum
   use rdb_barotropic_coupling, only: apply_bt_correction
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_visc_rem_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_visc_rem_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_drag_remnant_is_exactly_one", &
                               test_no_drag_remnant_is_one), &
                  new_unittest("strong_drag_remnant_matches_analytic_profile", &
                               test_analytic_profile), &
                  new_unittest("remnant_bounded_and_monotone", &
                               test_bounded_and_monotone), &
                  new_unittest("remnant_absent_is_bit_identical", &
                               test_absent_is_bit_identical), &
                  new_unittest("remnant_consistent_with_momentum_response", &
                               test_consistent_with_definition), &
                  new_unittest("producer_then_corrector_biases_against_bbl", &
                               test_end_to_end_wire) &
                  ]
   end subroutine collect_ocean_visc_rem_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vd)
      call vd%enter_data()
   end subroutine map_in

   subroutine map_out(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      call vd%exit_data()
      !$acc exit data delete(vd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------

   subroutine test_no_drag_remnant_is_one(error)
      !! Stratified h, large constant K_v, NO drag: the operator's row
      !! sums are exactly 1 (a no-flux/no-flux viscous operator cannot
      !! remove a uniform acceleration), so visc_rem must equal 1 to
      !! round-off — for any viscosity, any layer thicknesses, any dt.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 5
      real(wp), parameter :: DT = 900.0_wp
      integer :: k
      real(wp) :: max_dev
      checks: block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 3.5_wp        ! large interior viscosity
         vd%implicit_drag = .false.      ! NO drag: row sums stay exactly 1

         do k = 1, NZ
            ms%h_layer(:, :, k) = real(k, wp)*17.0_wp   ! stratified, non-uniform
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         max_dev = max(maxval(abs(bt%visc_rem_u - 1.0_wp)), &
                       maxval(abs(bt%visc_rem_v - 1.0_wp)))
         call check(error, max_dev < 1.0e-12_wp, &
                    "no-drag row-sum invariant: visc_rem must equal 1 to round-off")
      end block checks
      call bt%destroy(); call vd%destroy(); call ms%destroy()
   end subroutine test_no_drag_remnant_is_one

   subroutine test_analytic_profile(error)
      !! Single-column analytic target (§3.4): constant ν, uniform dz,
      !! deep water (D >> δ=√(νΔt)), near-no-slip bed (huge λ_bot).
      !! γ(z) ≈ 1 − exp(−(z+D)/δ) with z measured from the bed at
      !! z=-D (so z+D = height above bed).  Excludes the near-bed band
      !! (finite-λ departs from the continuum Dirichlet limit — the
      !! producer/consumer wire test 6 exercises that regime directly)
      !! and the top layer (finite-D free-surface BC).  Tolerance
      !! chosen empirically from the discrete solve (§9.2 spec quotes
      !! ~2-3% — the achieved max abs error here is ~2.4%).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 40
      real(wp), parameter :: DZ = 2.0_wp, NU = 0.1_wp, DT = 1000.0_wp
      real(wp), parameter :: LAM_HUGE = 1.0e6_wp
      real(wp), parameter :: DELTA = sqrt(NU*DT)   ! = 10 m
      integer, parameter :: K_LO = 8   ! exclude k < K_LO (near-bed band)
      real(wp), parameter :: TOL = 0.03_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      real(wp) :: z_above_bed, analytic, max_abs_err
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      checks: block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = NU
         vd%implicit_drag = .true.
         nu_face = size(ms%u_face_x_layer, 1)
         nv_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         nv_face = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         allocate (lam_u(nu_face, nv_uface), source=LAM_HUGE)
         allocate (lam_v(nx_vface, nv_face), source=LAM_HUGE)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         max_abs_err = 0.0_wp
         do k = K_LO, NZ - 1
            z_above_bed = (real(k, wp) - 0.5_wp)*DZ
            analytic = 1.0_wp - exp(-z_above_bed/DELTA)
            max_abs_err = max(max_abs_err, abs(bt%visc_rem_u(ig, jg, k) - analytic))
         end do
         call check(error, max_abs_err < TOL, &
                    "strong-drag remnant profile deviates from 1-exp(-(z+D)/sqrt(nu*dt))")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call bt%destroy(); call vd%destroy(); call ms%destroy()
   end subroutine test_analytic_profile

   subroutine test_bounded_and_monotone(error)
      !! Discrete maximum principle: for several ν (including ν=0) with
      !! strong bottom drag and stratified h, 0 < γ ≤ 1, γ is
      !! non-decreasing bed(k=1)→surface(k=nz), γ_1 < γ_nz strictly
      !! (λ_bot > 0), and no NaN.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 6
      real(wp), parameter :: DT = 450.0_wp, LAM = 0.08_wp
      real(wp) :: nu_list(3)
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      integer :: k, trial, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      logical :: bounded, monotone, no_nan, strict

      nu_list = [0.0_wp, 0.02_wp, 2.0_wp]
      do trial = 1, size(nu_list)
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = nu_list(trial)
         vd%implicit_drag = .true.
         nu_face = size(ms%u_face_x_layer, 1)
         nv_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         nv_face = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         do k = 1, NZ
            ms%h_layer(:, :, k) = real(k, wp)*12.0_wp
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         allocate (lam_u(nu_face, nv_uface), source=LAM)
         allocate (lam_v(nx_vface, nv_face), source=LAM)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         bounded = .true.
         monotone = .true.
         no_nan = .true.
         do k = 1, NZ
            if (bt%visc_rem_u(ig, jg, k) <= 0.0_wp .or. &
                bt%visc_rem_u(ig, jg, k) > 1.0_wp + 1.0e-12_wp) bounded = .false.
            if (bt%visc_rem_u(ig, jg, k) /= bt%visc_rem_u(ig, jg, k)) no_nan = .false.
            if (k > 1) then
               if (bt%visc_rem_u(ig, jg, k) < bt%visc_rem_u(ig, jg, k - 1) - 1.0e-13_wp) &
                  monotone = .false.
            end if
         end do
         strict = bt%visc_rem_u(ig, jg, 1) < bt%visc_rem_u(ig, jg, NZ) - 1.0e-10_wp

         deallocate (lam_u, lam_v)
         call bt%destroy(); call vd%destroy(); call ms%destroy()

         call check(error, bounded, "visc_rem: violated 0 < gamma <= 1")
         if (allocated(error)) return
         call check(error, no_nan, "visc_rem: NaN detected")
         if (allocated(error)) return
         call check(error, monotone, "visc_rem: not non-decreasing bed->surface")
         if (allocated(error)) return
         call check(error, strict, &
                    "visc_rem: bed value should be strictly less than surface (lambda_bot>0)")
         if (allocated(error)) return
      end do
   end subroutine test_bounded_and_monotone

   subroutine test_absent_is_bit_identical(error)
      !! Producer is a pure add-on: two identical states, A solved
      !! WITHOUT visc_rem_u/v, B solved WITH.  u_face/v_face_layer
      !! must be bit-identical — the remnant reads the surviving
      !! factorization and writes only visc_rem_out, so it cannot
      !! perturb the momentum answer it shares a kernel with (also
      !! catches an NVHPC FMA/codegen perturbation from the new
      !! in-kernel branch).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(barotropic_workstate_t) :: bt_b
      integer, parameter :: NZ = 6
      real(wp), parameter :: DT = 300.0_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      real(wp) :: max_du
      integer :: k, nu_face, nv_uface, nx_vface, nv_face
      checks: block
         call make_grid(grid, 5, 4)
         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid); call ms_b%init(grid)
         call vd_a%init(grid, nz_ml=NZ); call vd_b%init(grid, nz_ml=NZ)
         call bt_b%init(grid, nz_ml=NZ)
         vd_a%K_v_momentum = 2.0e-2_wp
         vd_b%K_v_momentum = 2.0e-2_wp
         vd_a%implicit_drag = .true.
         vd_b%implicit_drag = .true.
         nu_face = size(ms_a%u_face_x_layer, 1)
         nv_uface = size(ms_a%u_face_x_layer, 2)
         nx_vface = size(ms_a%v_face_y_layer, 1)
         nv_face = size(ms_a%v_face_y_layer, 2)
         allocate (lam_u(nu_face, nv_uface), source=0.02_wp)
         allocate (lam_v(nx_vface, nv_face), source=0.02_wp)

         do k = 1, NZ
            ms_a%h_layer(:, :, k) = real(k, wp)*8.0_wp
            ms_b%h_layer(:, :, k) = real(k, wp)*8.0_wp
            ms_a%u_face_x_layer(:, :, k) = 0.2_wp*real(k, wp)
            ms_b%u_face_x_layer(:, :, k) = 0.2_wp*real(k, wp)
            ms_a%v_face_y_layer(:, :, k) = -0.1_wp*real(k, wp)
            ms_b%v_face_y_layer(:, :, k) = -0.1_wp*real(k, wp)
         end do

         call map_in(ms_a, vd_a)
         call vdiff_apply_momentum(grid, vd_a, ms_a, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         call map_out(ms_a, vd_a)

         call map_in(ms_b, vd_b)
         call vdiff_apply_momentum(grid, vd_b, ms_b, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt_b%visc_rem_u, visc_rem_v=bt_b%visc_rem_v)
         call map_out(ms_b, vd_b)

         max_du = max(maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer)), &
                      maxval(abs(ms_a%v_face_y_layer - ms_b%v_face_y_layer)))
         call check(error, max_du == 0.0_wp, &
                    "visc_rem absent vs present: momentum solve must be bit-identical")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call bt_b%destroy()
      call vd_a%destroy(); call vd_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_absent_is_bit_identical

   subroutine test_consistent_with_definition(error)
      !! gamma_k IS (1/dt)*d(u_k^{n+1})/d(Abar) by construction: solve
      !! A with u^n = u0, B with u^n = u0 + dt*Abar (a spatially-uniform-
      !! in-k barotropic acceleration applied over the SAME dt as a
      !! uniform RHS-seed increment).  Because the backward-Euler
      !! operator is exactly linear, (u_B,k - u_A,k)/(dt*Abar) equals
      !! visc_rem_k to machine precision — no discretization error, a
      !! sharper check than the analytic-profile test.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(barotropic_workstate_t) :: bt_a
      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 240.0_wp, ABAR = 0.37_wp
      real(wp), parameter :: DELTA_U = DT*ABAR
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      real(wp) :: max_diff, fd_gamma
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      checks: block
         call make_grid(grid, 4, 4)
         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid); call ms_b%init(grid)
         call vd_a%init(grid, nz_ml=NZ); call vd_b%init(grid, nz_ml=NZ)
         call bt_a%init(grid, nz_ml=NZ)
         vd_a%K_v_momentum = 0.15_wp
         vd_b%K_v_momentum = 0.15_wp
         vd_a%implicit_drag = .true.
         vd_b%implicit_drag = .true.
         nu_face = size(ms_a%u_face_x_layer, 1)
         nv_uface = size(ms_a%u_face_x_layer, 2)
         nx_vface = size(ms_a%v_face_y_layer, 1)
         nv_face = size(ms_a%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost
         allocate (lam_u(nu_face, nv_uface), source=0.03_wp)
         allocate (lam_v(nx_vface, nv_face), source=0.03_wp)

         do k = 1, NZ
            ms_a%h_layer(:, :, k) = real(k, wp)*9.0_wp
            ms_b%h_layer(:, :, k) = real(k, wp)*9.0_wp
            ! Non-trivial, non-monotone u^n so the check isn't hiding a
            ! symmetry accident.
            ms_a%u_face_x_layer(:, :, k) = 0.3_wp*sin(real(k, wp))
            ms_b%u_face_x_layer(:, :, k) = 0.3_wp*sin(real(k, wp)) + DELTA_U
            ms_a%v_face_y_layer(:, :, k) = 0.0_wp
            ms_b%v_face_y_layer(:, :, k) = 0.0_wp
         end do

         call map_in(ms_a, vd_a)
         call vdiff_apply_momentum(grid, vd_a, ms_a, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt_a%visc_rem_u, visc_rem_v=bt_a%visc_rem_v)
         call map_out(ms_a, vd_a)

         call map_in(ms_b, vd_b)
         call vdiff_apply_momentum(grid, vd_b, ms_b, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         call map_out(ms_b, vd_b)

         max_diff = 0.0_wp
         do k = 1, NZ
            fd_gamma = (ms_b%u_face_x_layer(ig, jg, k) - ms_a%u_face_x_layer(ig, jg, k))/DELTA_U
            max_diff = max(max_diff, abs(fd_gamma - bt_a%visc_rem_u(ig, jg, k)))
         end do
         call check(error, max_diff < 1.0e-9_wp, &
                    "visc_rem must equal the finite-difference momentum response to Abar")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call bt_a%destroy()
      call vd_a%destroy(); call vd_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_consistent_with_definition

   subroutine test_end_to_end_wire(error)
      !! The trap-#1 defence: run the REAL producer
      !! (vdiff_apply_momentum) into bt_work%visc_rem_u/v, then feed
      !! that SAME field into the REAL consumer (apply_bt_correction,
      !! use_h_weighted + use_visc_rem).  Assert (a) visc_rem is no
      !! longer identically 1 after the producer call; (b) the bed
      !! layer's BT-corrector Δu is strictly smaller than the h-only
      !! path's (a barotropic acceleration no longer pushes water
      !! inside the frictional BBL as hard); (c) the column depth-mean
      !! Δu is unchanged to round-off (the mass-flux invariant survives
      !! the visc_rem re-weighting) in BOTH runs.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NZ = 4
      real(wp), parameter :: DT = 600.0_wp, LAM_STRONG = 0.05_wp
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      real(wp) :: max_dev, du_bed_a, du_bed_b, mean_a, mean_b
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      call make_grid(grid, 4, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      call bt%init(grid, nz_ml=NZ)
      ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
      call ms_a%init(grid); call ms_b%init(grid)
      call make_cartesian_metrics(metrics, grid)
      checks: block
         vd%K_v_momentum = 5.0e-2_wp
         vd%implicit_drag = .true.
         nu_face = size(ms%u_face_x_layer, 1)
         nv_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         nv_face = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost
         allocate (lam_u(nu_face, nv_uface), source=LAM_STRONG)
         allocate (lam_v(nx_vface, nv_face), source=LAM_STRONG)

         do k = 1, NZ
            ms%h_layer(:, :, k) = 25.0_wp   ! uniform h: isolate visc_rem's effect
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! ---- Step 9 (producer): the real vdiff kernel fills bt%visc_rem ----
         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         max_dev = maxval(abs(bt%visc_rem_u - 1.0_wp))
         call check(error, max_dev > 1.0e-3_wp, &
                    "producer: visc_rem_u should no longer be identically 1 after the "// &
                    "implicit-drag vdiff call")
         if (allocated(error)) exit checks

         ! ---- Step 7 (consumer, on TWO independent states): h-only vs
         ! h*visc_rem, using the SAME producer-filled bt%visc_rem ----
         ms_a%h_layer = 25.0_wp; ms_b%h_layer = 25.0_wp
         ms_a%u_face_x_layer = U0; ms_b%u_face_x_layer = U0
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp

         bt%bt_ubt_end = 1.0_wp; bt%bt_vbt_end = 0.0_wp
         bt%ubt_at_n = U0; bt%vbt_at_n = 0.0_wp
         bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 100.0_wp; bt%bt_eta_end = 0.0_wp

         call apply_bt_correction(bt, ms_a, DT, metrics, skip_h_rescale=.true., &
                                  use_h_weighted=.true., use_visc_rem=.false.)
         call apply_bt_correction(bt, ms_b, DT, metrics, skip_h_rescale=.true., &
                                  use_h_weighted=.true., use_visc_rem=.true.)

         du_bed_a = ms_a%u_face_x_layer(ig, jg, 1) - U0
         du_bed_b = ms_b%u_face_x_layer(ig, jg, 1) - U0
         call check(error, du_bed_b < du_bed_a - 1.0e-8_wp, &
                    "consumer: visc_rem-weighted bed Delta-u must be smaller than the "// &
                    "h-only path's (biased against the frictional BBL)")
         if (allocated(error)) exit checks

         mean_a = sum(ms_a%u_face_x_layer(ig, jg, :))/real(NZ, wp)
         mean_b = sum(ms_b%u_face_x_layer(ig, jg, :))/real(NZ, wp)
         call check(error, abs(mean_a - 1.0_wp) < 1.0e-12_wp, &
                    "h-only run: depth-mean must equal bt_ubt_end")
         if (allocated(error)) exit checks
         call check(error, abs(mean_b - 1.0_wp) < 1.0e-12_wp, &
                    "visc_rem run: depth-mean must STILL equal bt_ubt_end")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call destroy_cartesian_metrics(metrics)
      call bt%destroy()
      call vd%destroy(); call ms%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_end_to_end_wire

end module test_ocean_visc_rem
