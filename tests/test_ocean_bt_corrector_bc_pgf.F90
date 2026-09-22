!! Unit tests for the bc-PGF per-layer BT corrector
!! (`apply_bt_correction` with `use_bc_pgf = .true.`, plus the
!! supporting `compute_pbce / compute_gtot_faces / compute_e_anom`
!! kernels in `rdb_barotropic_coupling`).
!!
!! The bc-PGF correction is MOM6's `btstep_layer_accel` analogue
!! (`MOM_barotropic.F90:3720-3789`).  It adds a per-layer
!! baroclinic-PGF retro-correction on top of the uniform / h-weighted
!! BT Δu, accounting for the η evolution during the BT substep that
!! the slow PGF didn't see.
!!
!! Tests:
!!   1. `pbce_uniform_density_equals_g_scaled`
!!      — uniform rho_layer ⇒ pbce(k) = g · rho_ref/rho_0 for all k.
!!      ⇒ pbce − gtot ≡ 0 ⇒ bc-PGF Δ ≡ 0 everywhere.
!!   2. `pbce_two_layer_matches_closed_form`
!!      — 2-layer column (rho_bed=1036, rho_surf=1035, rho_ref=1035)
!!      with prescribed thicknesses; verify pbce follows the
!!      surface-to-bed recursion to roundoff.
!!   3. `gtot_depth_mean_zero_invariant`
!!      — for arbitrary pbce + h_layer, Σ_k h_face(k)·(pbce(k) −
!!      gtot_face) = 0 to roundoff (the mass-conservation guarantee).
!!   4. `bc_pgf_no_op_when_e_anom_zero`
!!      — with `use_bc_pgf=.true.` but `e_anom ≡ 0`, the resulting
!!      u_face_x_layer must be bit-identical to the `use_bc_pgf=.false.`
!!      run.
!!   5. `bc_pgf_depth_mean_preserved`
!!      — non-trivial e_anom + stratified pbce + stratified h; the
!!      column-mean of u_face_x_layer must still equal `bt_ubt_end`.
!!   6. `bc_pgf_opposite_sign_at_bed_vs_surf`
!!      — heavier bed (pbce(bed) > pbce(surf)) + positive e_anom on
!!      east cell, zero on west cell: the bc-PGF Δ at the bed has
!!      OPPOSITE sign from the Δ at the surface.  Qualitative
!!      direction-of-damping check.
module test_ocean_bt_corrector_bc_pgf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: apply_bt_correction, &
                                      compute_pbce, compute_gtot_faces, &
                                      compute_e_anom, snapshot_eta_PF
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       OPGF_VARIANT_FV_MOM6
   use rdb_scratch_3d, only: scratch_3d_buffer_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_bt_corrector_bc_pgf_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_bt_corrector_bc_pgf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("pbce_uniform_density_equals_g_scaled", test_pbce_uniform), &
                  new_unittest("pbce_two_layer_matches_closed_form", test_pbce_two_layer), &
                  new_unittest("gtot_depth_mean_zero_invariant", test_gtot_depth_mean_zero), &
                  new_unittest("bc_pgf_no_op_when_e_anom_zero", test_no_op_e_anom_zero), &
                  new_unittest("bc_pgf_depth_mean_preserved", test_depth_mean_preserved), &
                  new_unittest("bc_pgf_opposite_sign_at_bed_vs_surf", test_opposite_sign) &
                  ]
   end subroutine collect_bt_corrector_bc_pgf_tests

   subroutine build_state(grid, ms, bt_work, pgf, nz_ml)
      !! Build a minimal grid + state + workstate + pgf, with the FV_MOM6
      !! variant pre-selected on the pgf (so compute_pbce's assert passes).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      type(ocean_pressure_force_t), intent(out) :: pgf
      integer, intent(in) :: nz_ml
      call grid%init(4, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_ml
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=nz_ml)
      call pgf%init(grid, nz_ml=nz_ml)
      pgf%variant = OPGF_VARIANT_FV_MOM6
   end subroutine build_state

   subroutine cleanup(ms, bt_work, pgf)
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(ocean_pressure_force_t), intent(inout) :: pgf
      call pgf%destroy()
      call bt_work%destroy()
      call ms%destroy()
   end subroutine cleanup

   subroutine seed_e_face_uniform(pgf, ms, b)
      !! Populate pgf%e_face by the same recursion the FV_MOM6 PGF
      !! does: e(1) = -b; e(k+1) = e(k) + h_layer(k).
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: b
      integer :: i, j, k, nx, ny
      nx = size(pgf%e_face%data, 1)
      ny = size(pgf%e_face%data, 2)
      do j = 1, ny
         do i = 1, nx
            pgf%e_face%data(i, j, 1) = -b
            do k = 1, ms%nz_ml
               pgf%e_face%data(i, j, k + 1) = pgf%e_face%data(i, j, k) + ms%h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine seed_e_face_uniform

   ! -----------------------------------------------------------------
   ! Tests
   ! -----------------------------------------------------------------

   subroutine test_pbce_uniform(error)
      !! Uniform rho_layer everywhere ⇒ pbce(k) = g · rho_ref/rho_0
      !! for all k.  With rho_ref = rho_0 = 1035, pbce ≡ g everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: RHO = 1035.0_wp
      real(wp) :: pbce_min, pbce_max, expected
      integer :: i, j, k
      checks: block
         call build_state(grid, ms, bt, pgf, 3)
         call make_cartesian_metrics(metrics, grid)
         pgf%rho0 = RHO; pgf%rho_ref = RHO
         ms%rho_layer = RHO
         ms%h_layer = 100.0_wp
         call seed_e_face_uniform(pgf, ms, 300.0_wp)

         call compute_pbce(grid, bt, pgf, ms)

         pbce_min = minval(bt%pbce); pbce_max = maxval(bt%pbce)
         expected = GRAVITY*RHO/RHO
         call check(error, abs(pbce_min - expected) < 1.0e-10_wp, &
                    "uniform rho: pbce_min should equal g·rho_ref/rho_0")
         if (allocated(error)) exit checks
         call check(error, abs(pbce_max - expected) < 1.0e-10_wp, &
                    "uniform rho: pbce_max should equal g·rho_ref/rho_0")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt, pgf)
   end subroutine test_pbce_uniform

   subroutine test_pbce_two_layer(error)
      !! Two-layer 1036/1035 (rho_ref = 1035) with h_bed = 500, h_surf
      !! = 1000.  Closed form (our k=1 bed, k=nz=2 surf):
      !!   pbce(surf=2) = g · rho_ref / rho_0 = g
      !!   g_prime_at_K = g · (rho_layer(2) − rho_layer(1)) / rho_0
      !!                = g · (1035 − 1036) / 1035 = -g/1035
      !!   e_above_bed = e_face(2) = -b + h_bed   = -1500 + 500 = -1000
      !!   e_bed       = e_face(1) = -b                       = -1500
      !!   pbce(bed=1) = pbce(2) + g_prime_at_K · (e_above_bed − e_bed) / H_col
      !!              = g + (-g/1035) · 500/1500
      !!              = g · (1 − 1/(1035·3))
      !!              ≈ g · 0.999678
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: RHO_BED = 1036.0_wp, RHO_SURF = 1035.0_wp
      real(wp), parameter :: H_BED = 500.0_wp, H_SURF = 1000.0_wp
      real(wp) :: pbce_bed_expected, pbce_surf_expected, g_prime_K
      integer :: i, j
      checks: block
         call build_state(grid, ms, bt, pgf, 2)
         call make_cartesian_metrics(metrics, grid)
         pgf%rho0 = RHO_SURF; pgf%rho_ref = RHO_SURF

         ms%h_layer(:, :, 1) = H_BED
         ms%h_layer(:, :, 2) = H_SURF
         ms%rho_layer(:, :, 1) = RHO_BED
         ms%rho_layer(:, :, 2) = RHO_SURF
         call seed_e_face_uniform(pgf, ms, H_BED + H_SURF)

         call compute_pbce(grid, bt, pgf, ms)

         pbce_surf_expected = GRAVITY*RHO_SURF/RHO_SURF
         g_prime_K = GRAVITY*(RHO_SURF - RHO_BED)/RHO_SURF
         pbce_bed_expected = pbce_surf_expected &
                             + g_prime_K*H_BED/(H_BED + H_SURF)

         i = NGHOST + 2; j = NGHOST + 2
         call check(error, abs(bt%pbce(i, j, 2) - pbce_surf_expected) < 1.0e-12_wp, &
                    "two-layer: pbce(surf) deviates from closed form")
         if (allocated(error)) exit checks
         call check(error, abs(bt%pbce(i, j, 1) - pbce_bed_expected) < 1.0e-12_wp, &
                    "two-layer: pbce(bed) deviates from closed form")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt, pgf)
   end subroutine test_pbce_two_layer

   subroutine test_gtot_depth_mean_zero(error)
      !! For any (pbce, h), the column-mean of `pbce(k) − gtot_face`
      !! must be zero by construction of `gtot_face`.  Mass-flux
      !! invariant for the bc-PGF correction.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf
      real(wp) :: dev_sum, h_face, h_col, dev_max
      integer :: i, j, k, nz
      checks: block
         call build_state(grid, ms, bt, pgf, 3)
         call make_cartesian_metrics(metrics, grid)
         pgf%rho0 = 1035.0_wp; pgf%rho_ref = 1035.0_wp
         ms%rho_layer(:, :, 1) = 1036.0_wp     ! bed
         ms%rho_layer(:, :, 2) = 1035.5_wp     ! middle
         ms%rho_layer(:, :, 3) = 1035.0_wp     ! surf
         ms%h_layer(:, :, 1) = 50.0_wp
         ms%h_layer(:, :, 2) = 200.0_wp
         ms%h_layer(:, :, 3) = 1000.0_wp
         call seed_e_face_uniform(pgf, ms, 1250.0_wp)

         call compute_pbce(grid, bt, pgf, ms)
         call compute_gtot_faces(grid, bt, ms)

         nz = ms%nz_ml
         dev_max = 0.0_wp
         ! Check at one interior u-face (east face of cell (i,j))
         do j = NGHOST + 1, NGHOST + 4
            do i = NGHOST + 1, NGHOST + 3
               h_col = 0.0_wp
               dev_sum = 0.0_wp
               do k = 1, nz
                  h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i + 1, j, k))
                  h_col = h_col + h_face
                  dev_sum = dev_sum + h_face*(bt%pbce(i, j, k) - bt%gtot_E(i, j))
               end do
               dev_max = max(dev_max, abs(dev_sum/max(h_col, 1.0e-12_wp)))
            end do
         end do
         call check(error, dev_max < 1.0e-12_wp, &
                    "gtot_E: column-mean of (pbce − gtot) should be zero")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt, pgf)
   end subroutine test_gtot_depth_mean_zero

   subroutine test_no_op_e_anom_zero(error)
      !! With `e_anom = 0` everywhere, the bc-PGF Δ is zero regardless
      !! of pbce/gtot — so `use_bc_pgf=.true.` must produce bit-identical
      !! u_face_x_layer to `use_bc_pgf=.false.`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt_a, bt_b
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      real(wp) :: max_diff
      checks: block
         call build_state(grid, ms_a, bt_a, pgf_a, 2)
         call make_cartesian_metrics(metrics, grid)
         call build_state(grid, ms_b, bt_b, pgf_b, 2)
         pgf_a%rho0 = 1035.0_wp; pgf_a%rho_ref = 1035.0_wp
         pgf_b%rho0 = 1035.0_wp; pgf_b%rho_ref = 1035.0_wp

         ms_a%h_layer = 500.0_wp; ms_b%h_layer = 500.0_wp
         ms_a%rho_layer(:, :, 1) = 1036.0_wp; ms_a%rho_layer(:, :, 2) = 1035.0_wp
         ms_b%rho_layer(:, :, 1) = 1036.0_wp; ms_b%rho_layer(:, :, 2) = 1035.0_wp
         ms_a%u_face_x_layer = 0.1_wp; ms_b%u_face_x_layer = 0.1_wp
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp
         call seed_e_face_uniform(pgf_a, ms_a, 1000.0_wp)
         call seed_e_face_uniform(pgf_b, ms_b, 1000.0_wp)

         ! Both branches: a non-trivial BT correction, but e_anom = 0.
         bt_a%bt_ubt_end = 0.5_wp; bt_a%ubt_at_n = 0.1_wp; bt_a%F_bt_u = 0.0_wp
         bt_a%bt_vbt_end = 0.0_wp; bt_a%vbt_at_n = 0.0_wp; bt_a%F_bt_v = 0.0_wp
         bt_a%bt_H_ref = 1000.0_wp; bt_a%bt_eta_end = 0.0_wp
         bt_b%bt_ubt_end = 0.5_wp; bt_b%ubt_at_n = 0.1_wp; bt_b%F_bt_u = 0.0_wp
         bt_b%bt_vbt_end = 0.0_wp; bt_b%vbt_at_n = 0.0_wp; bt_b%F_bt_v = 0.0_wp
         bt_b%bt_H_ref = 1000.0_wp; bt_b%bt_eta_end = 0.0_wp

         call compute_pbce(grid, bt_b, pgf_b, ms_b)
         call compute_gtot_faces(grid, bt_b, ms_b)
         bt_b%e_anom = 0.0_wp     ! the critical setting

         call apply_bt_correction(bt_a, ms_a, 1.0_wp, metrics, skip_h_rescale=.true.)
         call apply_bt_correction(bt_b, ms_b, 1.0_wp, skip_h_rescale=.true., &
                                  grid=grid, use_bc_pgf=.true., metrics=metrics)

         max_diff = maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer))
         call check(error, max_diff < 1.0e-12_wp, &
                    "e_anom = 0: use_bc_pgf must be bit-identical to no-correction")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms_a, bt_a, pgf_a)
      call cleanup(ms_b, bt_b, pgf_b)
   end subroutine test_no_op_e_anom_zero

   subroutine test_depth_mean_preserved(error)
      !! After the bc-PGF correction, the depth-mean velocity must
      !! still equal `bt_ubt_end` to roundoff.  Stratified ρ, stratified
      !! h, non-trivial e_anom.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf
      real(wp) :: u_mean, h_total, err_max
      integer :: i, j, k, nz
      checks: block
         call build_state(grid, ms, bt, pgf, 3)
         call make_cartesian_metrics(metrics, grid)
         pgf%rho0 = 1035.0_wp; pgf%rho_ref = 1035.0_wp
         nz = ms%nz_ml

         ms%h_layer(:, :, 1) = 50.0_wp        ! thin bed
         ms%h_layer(:, :, 2) = 200.0_wp
         ms%h_layer(:, :, 3) = 1000.0_wp      ! thick surf
         ms%rho_layer(:, :, 1) = 1036.0_wp
         ms%rho_layer(:, :, 2) = 1035.5_wp
         ms%rho_layer(:, :, 3) = 1035.0_wp
         ms%u_face_x_layer = 0.0_wp; ms%v_face_y_layer = 0.0_wp
         call seed_e_face_uniform(pgf, ms, 1250.0_wp)

         bt%bt_ubt_end = 0.3_wp; bt%ubt_at_n = 0.0_wp; bt%F_bt_u = 0.0_wp
         bt%bt_vbt_end = 0.0_wp; bt%vbt_at_n = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 1250.0_wp; bt%bt_eta_end = 0.0_wp

         call compute_pbce(grid, bt, pgf, ms)
         call compute_gtot_faces(grid, bt, ms)
         ! Mimic a non-trivial SSH anomaly that varies in x.
         do j = 1, size(bt%e_anom, 2)
            do i = 1, size(bt%e_anom, 1)
               bt%e_anom(i, j) = 0.01_wp*real(i, wp)
            end do
         end do

         call apply_bt_correction(bt, ms, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true., grid=grid, &
                                  use_bc_pgf=.true., metrics=metrics)

         err_max = 0.0_wp
         do j = NGHOST + 1, NGHOST + 4
            do i = NGHOST + 2, NGHOST + 3
               u_mean = 0.0_wp; h_total = 0.0_wp
               do k = 1, nz
                  u_mean = u_mean + ms%u_face_x_layer(i, j, k)* &
                           0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
                  h_total = h_total + 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
               end do
               u_mean = u_mean/h_total
               err_max = max(err_max, abs(u_mean - bt%bt_ubt_end(i, j)))
            end do
         end do
         call check(error, err_max < 1.0e-12_wp, &
                    "depth-mean of u_layer must still equal bt_ubt_end")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt, pgf)
   end subroutine test_depth_mean_preserved

   subroutine test_opposite_sign(error)
      !! Heavier bed (pbce(bed) > pbce(surf)) + positive e_anom on
      !! east cell, zero on west cell.  The bc-PGF Δ should have
      !! OPPOSITE signs at the bed vs the surface layer (because
      !! pbce(k) − gtot is sign-opposite at bed vs surface).  This is
      !! the qualitative direction-of-damping check.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      type(ocean_pressure_force_t) :: pgf
      real(wp) :: du_bed, du_surf
      integer :: i, j
      checks: block
         call build_state(grid, ms, bt, pgf, 2)
         call make_cartesian_metrics(metrics, grid)
         pgf%rho0 = 1035.0_wp; pgf%rho_ref = 1035.0_wp
         ms%h_layer(:, :, 1) = 500.0_wp     ! bed
         ms%h_layer(:, :, 2) = 1000.0_wp    ! surf
         ms%rho_layer(:, :, 1) = 1036.0_wp  ! heavier
         ms%rho_layer(:, :, 2) = 1035.0_wp  ! lighter
         ms%u_face_x_layer = 0.0_wp; ms%v_face_y_layer = 0.0_wp
         call seed_e_face_uniform(pgf, ms, 1500.0_wp)

         ! Zero BT correction so the only u change is the bc-PGF term
         bt%bt_ubt_end = 0.0_wp; bt%ubt_at_n = 0.0_wp; bt%F_bt_u = 0.0_wp
         bt%bt_vbt_end = 0.0_wp; bt%vbt_at_n = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 1500.0_wp; bt%bt_eta_end = 0.0_wp

         call compute_pbce(grid, bt, pgf, ms)
         call compute_gtot_faces(grid, bt, ms)

         ! e_anom: positive only on the east cell (i=NGHOST+3), zero
         ! elsewhere.  At the u-face between (i=NGHOST+2) and
         ! (i=NGHOST+3), this gives a non-trivial Δ.
         bt%e_anom = 0.0_wp
         bt%e_anom(NGHOST + 3, NGHOST + 2) = 0.05_wp

         call apply_bt_correction(bt, ms, 1.0_wp, skip_h_rescale=.true., &
                                  grid=grid, use_bc_pgf=.true., metrics=metrics)

         i = NGHOST + 3; j = NGHOST + 2
         du_bed = ms%u_face_x_layer(i, j, 1)
         du_surf = ms%u_face_x_layer(i, j, 2)
         call check(error, du_bed*du_surf < 0.0_wp, &
                    "bc-PGF: bed and surface Δu must have opposite signs")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt, pgf)
   end subroutine test_opposite_sign

end module test_ocean_bt_corrector_bc_pgf
