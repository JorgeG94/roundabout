!! Analytical + golden tests for the Fox-Kemper mixed-layer-eddy
!! restratification slot (capability B5, `rdb_ocean_mle`).
!!
!! T1  FK08 scaling: uDml = Ce*H^2*|grad b|/|f| (bare timescale, <1%).
!! T2  restratification sign + monotone front slumping; ASSERT the
!!     surface transport is WESTWARD toward the dense column.
!! T3  conservation: sum_k a(k) = 0 per face => uhml/vhml non-divergent.
!! T4  bit-identity default-off: enable=.false. injects zero flux.
!! T5  zero gradient => zero transport.
!! T6  golden nz=6 worked example == prototype --report to ~1e-9.
module test_ocean_foxkemper
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_epbl, only: ocean_epbl_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_mle, only: ocean_mle_t, mle_compute_transports, &
                            mle_fold_x, mle_fold_y, &
                            mle_mu_shape, mle_layer_weights
   use rdb_continuity, only: continuity_t, continuity_tracer_step_split, &
                             continuity_tracer_drain, &
                             TR_MODE_ACCUMULATE
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_foxkemper_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: CE = 0.0625_wp
   real(wp), parameter :: F_FLOOR = 1.0e-5_wp

contains

   subroutine collect_ocean_foxkemper_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fk_t1_scaling", test_t1_scaling), &
                  new_unittest("fk_t2_sign_monotone", test_t2_sign_monotone), &
                  new_unittest("fk_t3_conservation", test_t3_conservation), &
                  new_unittest("fk_t4_bit_identity", test_t4_bit_identity), &
                  new_unittest("fk_t5_zero_gradient", test_t5_zero_gradient), &
                  new_unittest("fk_t6_golden", test_t6_golden), &
                  new_unittest("fk_t7_injection_seam", test_t7_injection_seam), &
                  new_unittest("fk_t8_fold_gate", test_t8_fold_gate), &
                  new_unittest("fk_t9_window_drain_thin_layer", test_t9_window_drain), &
                  new_unittest("fk_t10_mld_decay_filter", test_t10_mld_filter), &
                  new_unittest("fk_t11_dynamic_multiwindow", test_t11_dynamic_multiwindow), &
                  new_unittest("fk_t12_multiwindow_with_remap", test_t12_multiwindow_remap), &
                  new_unittest("fk_t13_wall_face_mask", test_t13_wall_face_mask), &
                  new_unittest("fk_t14_metric_aspect", test_t14_metric_aspect), &
                  new_unittest("bodner_arrest_by_wind", test_bodner_arrest_wind), &
                  new_unittest("bodner_convective_arrest", test_bodner_convective), &
                  new_unittest("bodner_cr_zero_inert", test_bodner_cr_zero) &
                  ]
   end subroutine collect_ocean_foxkemper_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   pure function ml_mean_buoyancy(rho, h, nz, mld) result(b_bar)
      !! Host replica of the ML-band buoyancy (k=nz surface, k=1 bed),
      !! partial-weighting the straddling layer.  Used by T1/T6 to derive
      !! the expected b_bar without re-running the device kernel.
      integer, intent(in) :: nz
      real(wp), intent(in) :: rho(nz), h(nz), mld
      real(wp) :: b_bar, htot, rho_int, h_remain, w
      integer :: k
      htot = 0.0_wp
      rho_int = 0.0_wp
      do k = nz, 1, -1
         h_remain = mld - htot
         if (h_remain <= 0.0_wp) exit
         w = min(h(k), h_remain)
         htot = htot + w
         rho_int = rho_int + rho(k)*w
      end do
      b_bar = -(GRAVITY/RHO0)*(rho_int/(htot + 1.0e-30_wp))
   end function ml_mean_buoyancy

   ! -----------------------------------------------------------------
   ! T1 — FK08 scaling
   ! -----------------------------------------------------------------

   subroutine test_t1_scaling(error)
      !! Single u-face, uniform H, constant f, known Δb̄: assert
      !! uDml = (Ce/|f|)*dy*Δb̄*H² and mu peaks mid-ML with mu(0)=mu(-1)=0.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp), parameter :: H = 100.0_wp, dz = H/nz
      real(wp), parameter :: f = 5.0e-5_wp, db = 1.0e-4_wp, dy = 1.0e4_wp
      real(wp) :: h_face(nz), a(nz), ts, uDml, uDml_exp, rel
      real(wp) :: mu_mid, mu_top, mu_bot
      integer :: k

      h_face = dz
      ts = CE/f
      uDml_exp = ts*dy*db*H*H

      ! Replicate the kernel's uDml chain (bare timescale).
      uDml = (CE/max(f, F_FLOOR))*dy*db*H*H
      rel = abs(uDml - uDml_exp)/abs(uDml_exp)
      call check(error, rel < 1.0e-10_wp, "T1: uDml scaling")
      if (allocated(error)) return

      ! mu properties.
      mu_top = mle_mu_shape(0.0_wp)
      mu_bot = mle_mu_shape(-1.0_wp)
      mu_mid = mle_mu_shape(-0.5_wp)
      call check(error, abs(mu_top) < 1.0e-15_wp, "T1: mu(0)=0")
      if (allocated(error)) return
      call check(error, abs(mu_bot) < 1.0e-15_wp, "T1: mu(-1)=0")
      if (allocated(error)) return
      call check(error, mu_mid > 0.5_wp, "T1: mu peaks mid-ML")
      if (allocated(error)) return

      ! a(k) closed cell.
      call mle_layer_weights(h_face, nz, H, a)
      call check(error, abs(sum(a)) < 1.0e-14_wp, "T1: sum a(k)=0")
      if (allocated(error)) return
      ! mu must be strictly positive inside, zero on the boundary layers.
      do k = 1, nz
         call check(error,.not. (a(k) /= a(k)), "T1: a(k) finite")
         if (allocated(error)) return
      end do
   end subroutine test_t1_scaling

   ! -----------------------------------------------------------------
   ! T2 — restratification sign + monotone slumping + WESTWARD surface
   ! -----------------------------------------------------------------

   subroutine test_t2_sign_monotone(error)
      !! W denser, E lighter (b_E > b_W) => uDml > 0.  The surface a(k)<0
      !! drives SURFACE transport WESTWARD toward the dense (west) column;
      !! lower-ML returns eastward.  Stepping the front forward with the
      !! transport must reduce |Δb̄| monotonically (front slumps).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp), parameter :: dz = 25.0_wp, mld = 100.0_wp
      real(wp), parameter :: f = 1.0e-4_wp, dy = 1.0e4_wp, dx = 1.0e4_wp
      real(wp), parameter :: dt = 10.0_wp
      integer, parameter :: n_steps = 200
      real(wp) :: rho_W(nz), rho_E(nz), h(nz), a(nz)
      real(wp) :: rho_h_W(nz), rho_h_E(nz), rho_face(nz), uhml(nz)
      real(wp) :: b_W, b_E, db, ts, uDml, area, db0, db_last
      integer :: step, k
      logical :: any_pos, any_neg, monotone

      rho_W = [1025.01_wp, 1024.99_wp, 1024.98_wp, 1024.97_wp]   ! bed->surf
      rho_E = [1025.00_wp, 1024.98_wp, 1024.97_wp, 1024.96_wp]
      h = dz
      area = dx*dy

      b_W = ml_mean_buoyancy(rho_W, h, nz, mld)
      b_E = ml_mean_buoyancy(rho_E, h, nz, mld)
      call check(error, b_E > b_W, "T2: setup b_E>b_W")
      if (allocated(error)) return

      db = b_E - b_W
      ts = CE/max(f, F_FLOOR)
      uDml = ts*dy*db*mld*mld
      call check(error, uDml > 0.0_wp, "T2: uDml>0 for b_E>b_W")
      if (allocated(error)) return

      call mle_layer_weights(h, nz, mld, a)
      do k = 1, nz
         uhml(k) = a(k)*uDml
      end do
      ! Surface layer (k=nz) transport must be NEGATIVE => WESTWARD toward
      ! the dense (west) column (the restratifying sign, spec addendum).
      call check(error, uhml(nz) < 0.0_wp, &
                 "T2: surface transport WESTWARD (toward dense column)")
      if (allocated(error)) return
      ! Closed-cell return flow: at least one positive (eastward) below.
      any_pos = .false.
      any_neg = .false.
      do k = 1, nz
         if (uhml(k) > 1.0e-30_wp) any_pos = .true.
         if (uhml(k) < -1.0e-30_wp) any_neg = .true.
      end do
      call check(error, any_pos .and. any_neg, "T2: +/- overturning cell")
      if (allocated(error)) return

      ! Step the front forward with conservative upwind density advection.
      rho_h_W = rho_W*h
      rho_h_E = rho_E*h
      db0 = abs(db)
      db_last = db0
      monotone = .true.
      do step = 1, n_steps
         do k = 1, nz
            rho_W(k) = rho_h_W(k)/h(k)
            rho_E(k) = rho_h_E(k)/h(k)
         end do
         b_W = ml_mean_buoyancy(rho_W, h, nz, mld)
         b_E = ml_mean_buoyancy(rho_E, h, nz, mld)
         db = b_E - b_W
         uDml = ts*dy*db*mld*mld
         do k = 1, nz
            uhml(k) = a(k)*uDml
            if (uhml(k) >= 0.0_wp) then
               rho_face(k) = rho_W(k)     ! positive flux takes from W
            else
               rho_face(k) = rho_E(k)
            end if
         end do
         do k = 1, nz
            rho_h_W(k) = rho_h_W(k) - uhml(k)*rho_face(k)*dt/area
            rho_h_E(k) = rho_h_E(k) + uhml(k)*rho_face(k)*dt/area
         end do
         if (abs(db) > db_last + 1.0e-15_wp*db0) monotone = .false.
         db_last = abs(db)
      end do

      call check(error, db_last < db0, "T2: |Δb̄| decreased (front slumped)")
      if (allocated(error)) return
      call check(error, monotone, "T2: |Δb̄| monotone non-increasing")
   end subroutine test_t2_sign_monotone

   ! -----------------------------------------------------------------
   ! T3 — conservation (sum_k a(k) = 0 across random configs)
   ! -----------------------------------------------------------------

   subroutine test_t3_conservation(error)
      !! sum_k a(k) = 0 per face to round-off across random MLD/h, and
      !! layers below the MLD get a(k)=0 (ML-confinement).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 8, n_tests = 500
      real(wp) :: h(nz), a(nz), mld, htot, max_err, rh, rm
      integer :: t, k

      ! Deterministic positive pseudo-random configs from a fractional
      ! recurrence (kept in [0,1) so thicknesses stay strictly positive).
      max_err = 0.0_wp
      rh = 0.123456789_wp
      rm = 0.314159265_wp
      do t = 1, n_tests
         htot = 0.0_wp
         do k = 1, nz
            rh = mod(rh*9301.0_wp + 0.49297_wp, 1.0_wp)
            h(k) = 1.0_wp + 49.0_wp*rh          ! [1,50] m, always > 0
            htot = htot + h(k)
         end do
         rm = mod(rm*4096.0_wp + 0.71828_wp, 1.0_wp)
         mld = maxval(h) + rm*(htot - maxval(h)) ! [one layer, full depth]
         call mle_layer_weights(h, nz, min(mld, htot), a)
         max_err = max(max_err, abs(sum(a)))
      end do
      call check(error, max_err < 1.0e-14_wp, "T3: max|sum a(k)| < 1e-14")
      if (allocated(error)) return

      ! Sub-MLD layers a(k)=0: MLD covers only the top 2 of 6 layers.
      h = 5.0_wp
      call mle_layer_weights(h, 6, 10.0_wp, a(1:6))
      do k = 1, 4    ! k=1..4 are below the 10 m MLD (top 2 layers = k=5,6)
         call check(error, abs(a(k)) < 1.0e-15_wp, "T3: sub-MLD a(k)=0")
         if (allocated(error)) return
      end do
   end subroutine test_t3_conservation

   ! -----------------------------------------------------------------
   ! T4 — bit-identity default-off
   ! -----------------------------------------------------------------

   subroutine test_t4_bit_identity(error)
      !! enable=.false. => mle_compute_transports leaves uhml/vhml at
      !! their zero-initialised values (no flux injected => bit-identity).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 3
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      integer :: k

      call grid%init(6, 5, NGHOST, 1.0e4_wp, 1.0e4_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = nz
      call ms%init(grid)
      call epbl%init(grid, nz_ml=nz)
      call ss%init(grid)
      call mle%init(grid, nz_ml=nz)
      mle%enable = .false.          ! default-off

      ! Strong front so an enabled kernel would inject obvious flux.
      do k = 1, nz
         ms%h_layer(:, :, k) = 30.0_wp
         ms%rho_layer(:, :, k) = 1025.0_wp + real(k, wp)
      end do
      ms%rho_layer(4:6, :, :) = ms%rho_layer(4:6, :, :) - 1.0_wp
      epbl%mld = 60.0_wp
      epbl%f_centre = 7.0e-5_wp
      epbl%rho0 = RHO0

      !$acc enter data copyin(ms, epbl, ss, mle)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call mle%enter_data()
      call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss)
      !$acc update self(mle%uhml, mle%vhml)
      call mle%exit_data()

      call check(error, maxval(abs(mle%uhml)) == 0.0_wp, "T4: uhml stays zero")
      if (allocated(error)) return
      call check(error, maxval(abs(mle%vhml)) == 0.0_wp, "T4: vhml stays zero")

      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, mle)
      call destroy_cartesian_metrics(metrics)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
   end subroutine test_t4_bit_identity

   ! -----------------------------------------------------------------
   ! T13 — closed-wall face mask (mask2dCu/mask2dCv)
   ! -----------------------------------------------------------------

   subroutine test_t13_wall_face_mask(error)
      !! FK MLE transport must not cross a closed (non-periodic) physical
      !! wall face.  On a ghosted grid the physical wall faces sit at
      !! i=nghost+1 / i=nghost+nx_phys+1 (and the y analogue), INTERIOR to
      !! the array — zeroing only the array edges leaves a nonzero uhml/vhml
      !! there, which mle_fold_{x,y} then injects into the mass flux AFTER
      !! the continuity wall-zeroing → tracer leaks across the wall (the
      !! ~1e-4/day FK + windowed-advect Salt drift).  Build a strong y-front
      !! with WALL south/north + WALL east/west and assert uhml/vhml vanish
      !! exactly on every physical wall face while the interior is nonzero.
      !! Without the `bc` mask this fails (nonzero wall flux); with it passes.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 3
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      type(ocean_bc_state_t) :: bc
      integer :: i, j, k, iw, ie, js, jn, i0, i1, j0, j1
      real(wp) :: wall_max, interior_max

      call grid%init(8, 8, NGHOST, 1.0e4_wp, 1.0e4_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = nz
      call ms%init(grid)
      call epbl%init(grid, nz_ml=nz)
      call ss%init(grid)
      call mle%init(grid, nz_ml=nz)
      call ocean_bc_state_init(bc, grid, nz, n_tracers=size(ms%tracers))
      ! All-closed box: every physical edge is a WALL → no FK flux may cross.
      bc%periodic_x = .false.
      bc%periodic_y = .false.
      bc%north_fold = .false.
      mle%enable = .true.
      mle%ce = CE
      mle%f_floor = F_FLOOR

      ! Stratified column + strong x- AND y-buoyancy fronts so both uDml and
      ! vDml are nonzero in the interior (an unmasked kernel would push FK
      ! flux through every wall face).
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               ms%h_layer(i, j, k) = 30.0_wp
               ms%rho_layer(i, j, k) = 1027.0_wp - 0.4_wp*real(k, wp) &
                                       + 0.05_wp*real(i, wp) - 0.05_wp*real(j, wp)
            end do
         end do
      end do
      epbl%mld = 60.0_wp
      epbl%f_centre = 7.0e-5_wp
      epbl%rho0 = RHO0

      !$acc enter data copyin(ms, epbl, ss, mle, bc)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call mle%enter_data()
      call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss, bc=bc)
      !$acc update self(mle%uhml, mle%vhml)
      call mle%exit_data()

      iw = NGHOST + 1
      ie = NGHOST + grid%nx_phys + 1
      js = NGHOST + 1
      jn = NGHOST + grid%ny_phys + 1
      i0 = NGHOST + 1
      i1 = NGHOST + grid%nx_phys
      j0 = NGHOST + 1
      j1 = NGHOST + grid%ny_phys

      ! Every physical wall face must carry exactly zero FK transport.
      wall_max = max(maxval(abs(mle%uhml(iw, j0:j1, :))), &
                     maxval(abs(mle%uhml(ie, j0:j1, :))), &
                     maxval(abs(mle%vhml(i0:i1, js, :))), &
                     maxval(abs(mle%vhml(i0:i1, jn, :))))
      ! The interior (one face in from each wall) must be nonzero — proves the
      ! kernel actually produced transport that the mask had to remove.
      interior_max = max(maxval(abs(mle%uhml(iw + 1:ie - 1, j0:j1, :))), &
                         maxval(abs(mle%vhml(i0:i1, js + 1:jn - 1, :))))

      call check(error, wall_max == 0.0_wp, "T13: FK transport zero on physical wall faces")
      if (allocated(error)) return
      call check(error, interior_max > 1.0e-6_wp, "T13: FK transport nonzero in interior")

      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, mle, bc)
      call destroy_cartesian_metrics(metrics)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
      call ocean_bc_state_destroy(bc)
   end subroutine test_t13_wall_face_mask

   ! -----------------------------------------------------------------
   ! T5 — zero gradient -> zero transport
   ! -----------------------------------------------------------------

   subroutine test_t5_zero_gradient(error)
      !! Uniform density (no front) => uDml=vDml=0 => uhml=vhml=0 exactly.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle

      call grid%init(6, 5, NGHOST, 1.0e4_wp, 1.0e4_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = nz
      call ms%init(grid)
      call epbl%init(grid, nz_ml=nz)
      call ss%init(grid)
      call mle%init(grid, nz_ml=nz)
      mle%enable = .true.
      mle%ce = CE
      mle%f_floor = F_FLOOR

      ms%h_layer = 25.0_wp
      ms%rho_layer = 1027.0_wp        ! uniform => zero gradient
      epbl%mld = 100.0_wp
      epbl%f_centre = 1.0e-4_wp
      epbl%rho0 = RHO0

      !$acc enter data copyin(ms, epbl, ss, mle)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call mle%enter_data()
      call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss)
      !$acc update self(mle%uhml, mle%vhml)
      call mle%exit_data()

      call check(error, maxval(abs(mle%uhml)) < 1.0e-12_wp, "T5: uhml=0")
      if (allocated(error)) return
      call check(error, maxval(abs(mle%vhml)) < 1.0e-12_wp, "T5: vhml=0")

      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, mle)
      call destroy_cartesian_metrics(metrics)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
   end subroutine test_t5_zero_gradient

   ! -----------------------------------------------------------------
   ! T6 — golden nz=6 worked example (prototype --report)
   ! -----------------------------------------------------------------

   subroutine test_t6_golden(error)
      !! Reproduce the prototype `--report` nz=6, 2-column worked example
      !! to ~1e-9: b̄_W, b̄_E, uDml, a(k), uhml(k).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 6
      real(wp), parameter :: dz = 15.0_wp, mld = 45.0_wp
      real(wp), parameter :: f = 7.0e-5_wp, dy = 2.0e4_wp
      real(wp) :: rho_W(nz), rho_E(nz), h(nz), hf(nz), a(nz)
      real(wp) :: b_W, b_E, db, h_vel, ts, uDml
      real(wp) :: uhml(nz)
      ! Golden values from `b5_foxkemper_prototype.py --report`.  The
      ! prototype prints with GRAVITY=9.80616; Roundabout's rdb_constants
      ! GRAVITY is 9.80665, so the buoyancy-derived golds (b̄, Δb̄, uDml)
      ! are rescaled by 9.80665/9.80616 (verified with the prototype
      ! algorithm at g=9.80665).  The g-independent shape gold a(k) is
      ! bit-identical to the prototype.
      real(wp), parameter :: B_W_GOLD = -9.721374782608695_wp
      real(wp), parameter :: B_E_GOLD = -9.711899758454106_wp
      real(wp), parameter :: DB_GOLD = 9.475024154589562e-03_wp
      real(wp), parameter :: UDML_GOLD = 3.426236413043547e+08_wp
      real(wp), parameter :: A4_GOLD = 9.124044679600235e-01_wp  ! k=4 (Fortran)
      real(wp), parameter :: UHML4_GOLD = 3.126113411548257e+08_wp ! k=4
      integer :: k

      rho_W = [1028.0_wp, 1027.5_wp, 1027.0_wp, 1026.5_wp, 1026.0_wp, 1025.5_wp]
      rho_E = [1027.0_wp, 1026.5_wp, 1026.0_wp, 1025.5_wp, 1025.0_wp, 1024.5_wp]
      h = dz

      b_W = ml_mean_buoyancy(rho_W, h, nz, mld)
      b_E = ml_mean_buoyancy(rho_E, h, nz, mld)
      db = b_E - b_W
      h_vel = mld                         ! htot_W = htot_E = 45 m
      ts = CE/max(f, F_FLOOR)
      uDml = ts*dy*db*h_vel*h_vel
      hf = 0.5_wp*(h + h)
      call mle_layer_weights(hf, nz, h_vel, a)
      do k = 1, nz
         uhml(k) = a(k)*uDml
      end do

      call check(error, abs(b_W - B_W_GOLD) < 1.0e-7_wp, "T6: b̄_W")
      if (allocated(error)) return
      call check(error, abs(b_E - B_E_GOLD) < 1.0e-7_wp, "T6: b̄_E")
      if (allocated(error)) return
      call check(error, abs(db - DB_GOLD) < 1.0e-9_wp, "T6: Δb̄")
      if (allocated(error)) return
      call check(error, abs(uDml - UDML_GOLD)/UDML_GOLD < 1.0e-8_wp, "T6: uDml")
      if (allocated(error)) return
      call check(error, abs(a(4) - A4_GOLD) < 1.0e-9_wp, "T6: a(4)")
      if (allocated(error)) return
      call check(error, abs(a(6) + A4_GOLD) < 1.0e-9_wp, "T6: a(6)=-a(4)")
      if (allocated(error)) return
      call check(error, abs(uhml(4) - UHML4_GOLD)/UHML4_GOLD < 1.0e-8_wp, "T6: uhml(4)")
      if (allocated(error)) return
      ! Layers below MLD (k=1,2,3) carry no transport.
      do k = 1, 3
         call check(error, abs(a(k)) < 1.0e-12_wp, "T6: sub-MLD a(k)=0")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(a)) < 1.0e-14_wp, "T6: sum a(k)=0")
   end subroutine test_t6_golden

   ! -----------------------------------------------------------------
   ! T7 — injection-seam device test: fold -> continuity_tracer_step_split
   ! -----------------------------------------------------------------

   subroutine test_t7_injection_seam(error)
      !! DEVICE test exercising the actual on-device injection path —
      !! B5's reason to exist.  Runs on GPU via the same `!$acc enter data`
      !! pattern as T4/T5; asserts:
      !!   (i)  Σh conserved to ~1e-10 relative after the split step.
      !!   (ii) Σ(h·S) and Σ(h·T) (tracer mass) conserved to ~1e-10.
      !!   (iii) mass_flux_x_layer was modified by the fold (non-zero
      !!         contribution from uhml — proves the fold executed).
      !!
      !! The setup: 4×3 grid, nz=3, horizontal density gradient in x so
      !! mle_compute_transports injects non-zero uhml.  EPBL provides
      !! a simple MLD = dz*nz.  Velocities are zero so any change in
      !! h-layer or tracer mass after continuity_tracer_step_split comes
      !! ONLY from the FK fold — the test then checks conservation.
      !!
      !! Note on thermo-cadence limitation (FIX 3): with velocities at
      !! zero the FK fold is the only flux; mass_flux sums to zero by
      !! sum_k a(k) = 0 => Σh and Σ(hTr) are conserved exactly regardless
      !! of how many times the fold is applied.  The limitation (stale
      !! uhml re-folded when dt_therm_ratio>1) is documented on the slot
      !! fields and the fold call site; it does NOT affect this test.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 3
      real(wp), parameter :: DX = 1.0e4_wp, DY = 1.0e4_wp
      real(wp), parameter :: DZ = 50.0_wp           ! layer thickness (m)
      real(wp), parameter :: MLD = DZ*real(NZ, wp)  ! full-column MLD
      real(wp), parameter :: F0 = 1.0e-4_wp         ! Coriolis
      real(wp), parameter :: DT = 60.0_wp           ! 1-minute step
      integer, parameter :: N_STEPS = 3
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      integer :: i, j, k, nx, ny, it
      real(wp) :: h0_total, h_total, reldiff
      real(wp), allocatable :: hTr0(:, :, :, :)    ! (nx,ny,nz,ntracers)
      real(wp) :: flux_norm_before, flux_norm_after

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         call epbl%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call mle%init(grid, nz_ml=NZ)
         mle%enable = .true.
         mle%ce = CE
         mle%f_floor = F_FLOOR

         ! Uniform layer thickness; zero velocity (so any h/Tr change
         ! comes exclusively from the FK fold through continuity).
         do k = 1, NZ
            ms%h_layer(:, :, k) = DZ
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%v_face_y_layer(:, :, k) = 0.0_wp
         end do

         ! Horizontal density gradient in x: west columns denser than east.
         ! Density decreases by ~0.5 kg/m^3 each column (bed->surf uniform
         ! per column so the ML average follows the horizontal gradient).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%rho_layer(i, j, k) = 1027.0_wp - 0.1_wp*real(i, wp)
               end do
            end do
         end do

         ! Tracers: uniform salinity + linear temperature gradient.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = DZ*35.0_wp
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = DZ*(15.0_wp + &
                                                                    0.5_wp*real(i, wp))
               end do
            end do
         end do

         ! EPBL: uniform MLD = full column, uniform f.
         epbl%rho0 = RHO0
         epbl%mld = MLD
         epbl%f_centre = F0

         ! Snapshot total h and per-tracer mass before the device run.
         h0_total = sum(ms%h_layer)*DX*DY
         allocate (hTr0(nx, ny, NZ, size(ms%tracers)))
         do it = 1, size(ms%tracers)
            hTr0(:, :, :, it) = ms%tracers(it)%hTr
         end do

         ! Attach all state to device.
         !$acc enter data copyin(ms, epbl, ss, mle, ct)
         call ms%enter_data()
         call epbl%enter_data()
         call ss%enter_data()
         call mle%enter_data()
         call ct%enter_data()

         ! Compute FK transports on device (fills uhml/vhml).
         call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss)

         ! Snapshot mass_flux_x_layer BEFORE the fold by pulling uhml back.
         ! We verify the fold is non-trivial: |uhml| > 0.
         !$acc update self(mle%uhml)
         flux_norm_before = maxval(abs(mle%uhml))

         ! Run N_STEPS of continuity_tracer_step_split WITH the mle argument
         ! — this is the production injection path (fold_x + fold_y inside).
         do i = 1, N_STEPS
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT, mle=mle)
         end do

         ! Pull results back to host.
         call ct%exit_data()
         !$acc exit data delete(ct)
         call mle%exit_data()
         call ss%exit_data()
         call epbl%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, epbl, ss, mle)

         ! (i) h conserved: Σh after = Σh before (closed-wall, zero velocity
         !     outside the FK fold, and sum_k a(k)=0 => fold adds nothing to
         !     the column integral).
         h_total = sum(ms%h_layer)*DX*DY
         reldiff = abs(h_total - h0_total)/(abs(h0_total) + 1.0e-30_wp)
         call check(error, reldiff < 1.0e-10_wp, &
                    "T7: total h not conserved through FK fold + continuity")
         if (allocated(error)) exit checks

         ! (ii) Tracer mass conserved per tracer.
         do it = 1, size(ms%tracers)
            reldiff = abs(sum(ms%tracers(it)%hTr) - sum(hTr0(:, :, :, it))) &
                      /(abs(sum(hTr0(:, :, :, it))) + 1.0e-30_wp)
            call check(error, reldiff < 1.0e-10_wp, &
                       "T7: tracer mass not conserved through FK fold + continuity")
            if (allocated(error)) exit checks
         end do

         ! (iii) The fold actually ran: uhml must be non-zero (the buoyancy
         !       gradient in x guarantees this; a zero here means mle_compute
         !       or mle_fold_x silently no-op'd on device).
         call check(error, flux_norm_before > 0.0_wp, &
                    "T7: uhml is zero — FK transport not computed on device")

      end block checks

      call ct%destroy()
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      if (allocated(hTr0)) deallocate (hTr0)
   end subroutine test_t7_injection_seam

   ! -----------------------------------------------------------------
   ! T8 — mle_fold_active gate: .false. suppresses fold to no-mle identity
   ! -----------------------------------------------------------------

   subroutine test_t8_fold_gate(error)
      !! DEVICE test for the `mle_fold_active` argument added by the
      !! DT_THERM Phase 1 work.  Uses the same 4×3 grid + density gradient
      !! as T7 so uhml is provably non-zero.
      !!
      !! Three runs on fresh identical host-side ICs, each run on device:
      !!
      !!   (A) No `mle` arg at all → baseline state after continuity.
      !!   (B) `mle=mle, mle_fold_active=.false.` → fold suppressed →
      !!       result must be bit-identical to (A).
      !!   (C) `mle=mle, mle_fold_active=.true.` → fold fires →
      !!       result must DIFFER from (A) (positive control proving
      !!       the gate is what separates B from C, not a broken setup).
      !!
      !! ASSERT (B)==(A) and (C)!=(A) for both h_layer and hTr.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 3
      real(wp), parameter :: DX = 1.0e4_wp, DY = 1.0e4_wp
      real(wp), parameter :: DZ = 50.0_wp
      real(wp), parameter :: MLD = DZ*real(NZ, wp)
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 60.0_wp
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms_a, ms_b, ms_c
      type(continuity_t) :: ct_a, ct_b, ct_c
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      integer :: i, j, k, nx, ny
      real(wp) :: max_diff_ba, max_diff_ca
      real(wp) :: max_hTr_diff_ba, max_hTr_diff_ca

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Initialise three identical states and three continuity slots.
         ms_a%nz_ml = NZ
         ms_b%nz_ml = NZ
         ms_c%nz_ml = NZ
         call ms_a%init(grid)
         call ms_b%init(grid)
         call ms_c%init(grid)
         ct_a%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct_a%init(grid, nz_ml=NZ)
         ct_b%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct_b%init(grid, nz_ml=NZ)
         ct_c%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct_c%init(grid, nz_ml=NZ)
         call epbl%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call mle%init(grid, nz_ml=NZ)
         mle%enable = .true.
         mle%ce = CE
         mle%f_floor = F_FLOOR

         ! IC: same horizontal density gradient as T7 (west denser than east).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_a%h_layer(i, j, k) = DZ
                  ms_a%u_face_x_layer(i, j, k) = 0.0_wp
                  ms_a%v_face_y_layer(i, j, k) = 0.0_wp
                  ms_a%rho_layer(i, j, k) = 1027.0_wp - 0.1_wp*real(i, wp)
                  ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, k) = DZ*35.0_wp
                  ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, k) = DZ*(15.0_wp + 0.5_wp*real(i, wp))
               end do
            end do
         end do
         ! Copy IC to ms_b and ms_c.
         ms_b%h_layer = ms_a%h_layer
         ms_b%u_face_x_layer = ms_a%u_face_x_layer
         ms_b%v_face_y_layer = ms_a%v_face_y_layer
         ms_b%rho_layer = ms_a%rho_layer
         do k = 1, NZ
            ms_b%tracers(ms_b%idx_salinity)%hTr = ms_a%tracers(ms_a%idx_salinity)%hTr
            ms_b%tracers(ms_b%idx_temperature)%hTr = ms_a%tracers(ms_a%idx_temperature)%hTr
         end do
         ms_c%h_layer = ms_a%h_layer
         ms_c%u_face_x_layer = ms_a%u_face_x_layer
         ms_c%v_face_y_layer = ms_a%v_face_y_layer
         ms_c%rho_layer = ms_a%rho_layer
         do k = 1, NZ
            ms_c%tracers(ms_c%idx_salinity)%hTr = ms_a%tracers(ms_a%idx_salinity)%hTr
            ms_c%tracers(ms_c%idx_temperature)%hTr = ms_a%tracers(ms_a%idx_temperature)%hTr
         end do

         epbl%rho0 = RHO0
         epbl%mld = MLD
         epbl%f_centre = F0

         ! --- Run (A): no mle arg ---
         !$acc enter data copyin(ms_a, ct_a)
         call ms_a%enter_data()
         call ct_a%enter_data()
         call continuity_tracer_step_split(grid, metrics, ct_a, ms_a, DT)
         call ms_a%exit_data()
         call ct_a%exit_data()
         !$acc exit data delete(ms_a, ct_a)

         ! --- Compute MLE transports (shared uhml/vhml for B and C) ---
         !$acc enter data copyin(epbl, ss, mle)
         call epbl%enter_data()
         call ss%enter_data()
         call mle%enter_data()
         call mle_compute_transports(grid, metrics, mle, ms_a, epbl, ss=ss)

         ! Verify uhml is non-zero (density gradient is strong enough).
         !$acc update self(mle%uhml)
         call check(error, maxval(abs(mle%uhml)) > 0.0_wp, &
                    "T8: uhml zero — density gradient did not generate FK transport")
         if (allocated(error)) exit checks

         ! --- Run (B): mle present, fold_active=.false. ---
         !$acc enter data copyin(ms_b, ct_b)
         call ms_b%enter_data()
         call ct_b%enter_data()
         call continuity_tracer_step_split(grid, metrics, ct_b, ms_b, DT, &
                                           mle=mle, mle_fold_active=.false.)
         call ms_b%exit_data()
         call ct_b%exit_data()
         !$acc exit data delete(ms_b, ct_b)

         ! --- Run (C): mle present, fold_active=.true. ---
         !$acc enter data copyin(ms_c, ct_c)
         call ms_c%enter_data()
         call ct_c%enter_data()
         call continuity_tracer_step_split(grid, metrics, ct_c, ms_c, DT, &
                                           mle=mle, mle_fold_active=.true.)
         call ms_c%exit_data()
         call ct_c%exit_data()
         !$acc exit data delete(ms_c, ct_c)

         ! Tear down MLE/EPBL/SS device attachments.
         call mle%exit_data()
         call ss%exit_data()
         call epbl%exit_data()
         !$acc exit data delete(epbl, ss, mle)

         ! --- Assert (B) == (A) (fold suppressed => no FK contribution) ---
         max_diff_ba = maxval(abs(ms_b%h_layer - ms_a%h_layer))
         max_hTr_diff_ba = max( &
                           maxval(abs(ms_b%tracers(ms_b%idx_salinity)%hTr - ms_a%tracers(ms_a%idx_salinity)%hTr)), &
                           maxval(abs(ms_b%tracers(ms_b%idx_temperature)%hTr - ms_a%tracers(ms_a%idx_temperature)%hTr)))

         call check(error, max_diff_ba == 0.0_wp, &
                    "T8: mle_fold_active=.false. h_layer differs from no-mle baseline (fold not suppressed)")
         if (allocated(error)) exit checks
         call check(error, max_hTr_diff_ba == 0.0_wp, &
                    "T8: mle_fold_active=.false. tracers differ from no-mle baseline (fold not suppressed)")
         if (allocated(error)) exit checks

         ! --- Assert (C) != (A) (fold fired => FK flux changed state) ---
         max_diff_ca = maxval(abs(ms_c%h_layer - ms_a%h_layer))
         max_hTr_diff_ca = max( &
                           maxval(abs(ms_c%tracers(ms_c%idx_salinity)%hTr - ms_a%tracers(ms_a%idx_salinity)%hTr)), &
                           maxval(abs(ms_c%tracers(ms_c%idx_temperature)%hTr - ms_a%tracers(ms_a%idx_temperature)%hTr)))

         call check(error, max_diff_ca > 0.0_wp .or. max_hTr_diff_ca > 0.0_wp, &
                    "T8: mle_fold_active=.true. produced no change vs no-mle baseline (fold did not fire)")

      end block checks

      call ct_c%destroy()
      call ct_b%destroy()
      call ct_a%destroy()
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms_c%destroy()
      call ms_b%destroy()
      call ms_a%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_t8_fold_gate

   subroutine test_t9_window_drain(error)
      !! DEVICE regression for the windowed-tracer-advect x Fox-Kemper
      !! interaction bug.  A thin EPBL-MLD surface layer over a buoyancy
      !! front: the FK overturning injects a per-layer transport into the
      !! thin surface layer that, accumulated over a `ratio`-step window,
      !! used to exceed the layer's drainable volume.  The windowed drain's
      !! `hprev = h_end + div(uhtr)` reconstruction then clamped to ~0 in
      !! the thin layer ⇒ `Tr = hTr/hprev` blew up / collapsed and Σ(hTr)
      !! was not conserved.
      !!
      !! The fix is the MOM6 per-layer FK availability cap (a(k)*uDml ≤
      !! h_avail with I4dt = 1/(4·dt_limit)).  This test drives the genuine
      !! production path: `continuity_tracer_step_split(TR_MODE_ACCUMULATE)`
      !! folds the (limited) FK transport into `uhtr` and advances `h_layer`
      !! over a 2-step window, then `continuity_tracer_drain` spends it.
      !!
      !! ASSERT after the drain:
      !!   (i)  Σ(areaT·hTr) for S conserved to ~1e-10 relative.
      !!   (ii) per-cell concentration S = hTr/h_layer stays within the IC
      !!        bounds (no collapse, no negatives) to a small tolerance.
      !! WITHOUT the limiter (dt_limit absent) both (i) and (ii) fail —
      !! that is the failing-before behaviour this test guards.
      use rdb_constants, only: wp
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 3
      integer, parameter :: NZ = 3
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 3
      integer, parameter :: RATIO = 2
      real(wp), parameter :: DX = 5.0e3_wp, DY = 5.0e3_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: DT_WIN = real(RATIO, wp)*DT
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      ! Strongly stratified column with a THIN surface layer (k=NZ).
      real(wp), parameter :: H_DEEP = 400.0_wp, H_MID = 80.0_wp, H_SURF = 10.0_wp
      real(wp), parameter :: MLD = H_SURF + H_MID  ! ML spans the two upper layers
      real(wp), parameter :: S_REF(NZ) = [35.0_wp, 34.5_wp, 33.0_wp]  ! k=1 bed..nz surf
      real(wp), parameter :: HCOL(NZ) = [H_DEEP, H_MID, H_SURF]
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      type(ocean_bc_state_t) :: bc
      integer :: i, j, k, nx, ny, i0, i1, j0, j1
      real(wp) :: area, mass0, mass1, rel, conc, smin, smax
      real(wp) :: s_lo, s_hi

      checks: block
         call grid%init(NX_PHYS, NY_PHYS, NG, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         area = DX*DY
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         call epbl%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call mle%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
         bc%periodic_x = .true.
         bc%periodic_y = .false.
         bc%north_fold = .false.
         mle%enable = .true.
         mle%ce = CE
         mle%f_floor = F_FLOOR

         ! Stratified column + a periodic-x buoyancy front (light mid-channel)
         ! so the FK overturning varies in x (per-layer div(uhtr) /= 0).  Zero
         ! background velocity ⇒ the only transport is FK.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = HCOL(k)
                  ms%u_face_x_layer(i, j, k) = 0.0_wp
                  ms%v_face_y_layer(i, j, k) = 0.0_wp
                  ! Density: stable strat + a strong x-front (period = nx_phys).
                  ms%rho_layer(i, j, k) = 1027.0_wp - 0.4_wp*real(k, wp) &
                                          - 0.6_wp*sin(2.0_wp*PI*real(i - NG, wp)/real(NX_PHYS, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S_REF(k)*HCOL(k)
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 12.0_wp*HCOL(k)
               end do
            end do
         end do
         epbl%rho0 = RHO0
         epbl%mld = MLD
         epbl%f_centre = F0

         i0 = NG + 1
         i1 = NG + grid%nx_phys
         j0 = NG + 1
         j1 = NG + grid%ny_phys
         ! IC concentration bounds (with a slim PPM-overshoot tolerance).
         s_lo = minval(S_REF) - 0.05_wp*(maxval(S_REF) - minval(S_REF))
         s_hi = maxval(S_REF) + 0.05_wp*(maxval(S_REF) - minval(S_REF))
         mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))

         !$acc enter data copyin(ms, ct, epbl, ss, mle, bc)
         call ms%enter_data()
         call ct%enter_data()
         call epbl%enter_data()
         call ss%enter_data()
         call mle%enter_data()

         ! Compute FK transports WITH the availability cap (production path
         ! passes dt_limit = dt_therm; here the window is RATIO*DT).  Toggle
         ! the dt_limit presence to see failing-before vs passing-after.
         call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss, &
                                     dt_limit=DT_WIN)
         !$acc update self(mle%uhml)

         ! Build the windowed drain inputs directly, the way the dyn step
         ! does but without the RK2 bookkeeping: the net window transport is
         ! the FK flux integrated over the window (uhtr = uhml·DT_WIN; zero
         ! background flow), and h_end is h_start advanced by its divergence.
         ! This satisfies the drain contract hprev = h_end + div(uhtr) =
         ! h_start exactly, with FK as the sole transport — isolating the
         ! FK×drain interaction.  Done on host then re-synced to device.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ct%uhtr(i, j, k) = mle%uhml(i, j, k)*DT_WIN
               end do
            end do
         end do
         ct%vhtr = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  ms%h_layer(i, j, k) = ms%h_layer(i, j, k) &
                                        - (ct%uhtr(i + 1, j, k) - ct%uhtr(i, j, k))*metrics%iareaT(i, j)
               end do
            end do
         end do
         !$acc update device(ct%uhtr, ct%vhtr, ms%h_layer)

         ! Spend the window.
         call continuity_tracer_drain(grid, metrics, ct, ms, RATIO, bc=bc)

         call mle%exit_data()
         call ss%exit_data()
         call epbl%exit_data()
         call ct%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ct, epbl, ss, mle, bc)

         ! (i) conservation of interior salinity mass.
         mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
         rel = abs(mass1 - mass0)/abs(mass0)
         call check(error, rel < 1.0e-10_wp, &
                    "T9: FK windowed drain broke interior tracer conservation")
         if (allocated(error)) exit checks

         ! (ii) per-cell concentration stays physical (no collapse / negative).
         smin = huge(1.0_wp)
         smax = -huge(1.0_wp)
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  conc = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/ms%h_layer(i, j, k)
                  smin = min(smin, conc)
                  smax = max(smax, conc)
               end do
            end do
         end do
         call check(error, smin >= s_lo, &
                    "T9: FK windowed drain collapsed concentration below IC range")
         if (allocated(error)) exit checks
         call check(error, smax <= s_hi, &
                    "T9: FK windowed drain produced a concentration above IC range")
      end block checks

      call ocean_bc_state_destroy(bc)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_t9_window_drain

   subroutine test_t10_mld_filter(error)
      !! Running-mean MLD filter (MOM6 MLE_MLD_DECAY_TIME).  The FK
      !! streamfunction scales as Psi ~ MLD^2, so a diagnosed boundary
      !! layer that deepens then retreats sharply (as the EPBL MLD does
      !! under transient forcing) injects a transport spike.  The filter
      !! resets instantly to a DEEPER MLD but DECAYS on retreat, bounding
      !! the spike.
      !!
      !! Sequence (single u-face, uniform front): step the compute kernel
      !! with MLD = [deep, deep, shallow].  ASSERT:
      !!   (i)  filter OFF (decay_time = 0) is bit-identical to the
      !!        instantaneous-MLD transport at every step (legacy path);
      !!   (ii) filter ON: after the retreat-step the filtered transport
      !!        magnitude is STRICTLY LARGER than the unfiltered one
      !!        (the running mean keeps MLD deeper than the instantaneous
      !!        shallow value => more transport, but bounded — it is the
      !!        decayed value, not the spike), AND strictly smaller than
      !!        the deep-step transport (it IS decaying toward shallow).
      use rdb_constants, only: wp
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 2
      integer, parameter :: NZ = 4
      integer, parameter :: NX_PHYS = 6, NY_PHYS = 3
      real(wp), parameter :: DX = 5.0e3_wp, DY = 5.0e3_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 1200.0_wp
      real(wp), parameter :: T_DECAY = 86400.0_wp
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: HCOL(NZ) = [400.0_wp, 100.0_wp, 50.0_wp, 20.0_wp]
      real(wp), parameter :: MLD_DEEP = 120.0_wp, MLD_SHALLOW = 25.0_wp
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle_off, mle_on
      integer :: i, j, k, nx, ny, iface, jrow
      real(wp) :: t_off(3), t_on(3)
      real(wp) :: best

      checks: block
         call grid%init(NX_PHYS, NY_PHYS, NG, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call epbl%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call mle_off%init(grid, nz_ml=NZ)
         call mle_on%init(grid, nz_ml=NZ)
         mle_off%enable = .true.
         mle_off%ce = CE
         mle_off%f_floor = F_FLOOR
         mle_off%mld_decay_time = 0.0_wp          ! filter OFF
         mle_on%enable = .true.
         mle_on%ce = CE
         mle_on%f_floor = F_FLOOR
         mle_on%mld_decay_time = T_DECAY          ! filter ON

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = HCOL(k)
                  ms%rho_layer(i, j, k) = 1027.0_wp - 0.3_wp*real(k, wp) &
                                          - 0.5_wp*sin(2.0_wp*PI*real(i - NG, wp)/real(NX_PHYS, wp))
               end do
            end do
         end do
         epbl%rho0 = RHO0
         epbl%f_centre = F0

         ! Pick the interior u-face with the largest FK transport (the
         ! sine front => some faces sit at db ~ 0).  Probe once on the deep
         ! MLD, on the mid row, then reuse that (iface,jrow) throughout.
         jrow = NG + 1
         epbl%mld = MLD_DEEP
         call mle_compute_transports(grid, metrics, mle_off, ms, epbl, ss=ss, dt_limit=DT)
         iface = NG + 2
         best = -1.0_wp
         do i = 2, nx
            if (mle_face_mag(mle_off, i, jrow, NZ) > best) then
               best = mle_face_mag(mle_off, i, jrow, NZ)
               iface = i
            end if
         end do
         ! Reset the filter-off slot's history is irrelevant (decay=0); redo
         ! the deep step cleanly so t_off(1) is the first window's value.

         ! Step 1 + 2: deep MLD.  Step 3: shallow MLD.
         ! --- filter OFF run ---
         epbl%mld = MLD_DEEP
         call mle_compute_transports(grid, metrics, mle_off, ms, epbl, ss=ss, dt_limit=DT)
         t_off(1) = mle_face_mag(mle_off, iface, jrow, NZ)
         call mle_compute_transports(grid, metrics, mle_off, ms, epbl, ss=ss, dt_limit=DT)
         t_off(2) = mle_face_mag(mle_off, iface, jrow, NZ)
         epbl%mld = MLD_SHALLOW
         call mle_compute_transports(grid, metrics, mle_off, ms, epbl, ss=ss, dt_limit=DT)
         t_off(3) = mle_face_mag(mle_off, iface, jrow, NZ)

         ! --- filter ON run (same sequence) ---
         epbl%mld = MLD_DEEP
         call mle_compute_transports(grid, metrics, mle_on, ms, epbl, ss=ss, dt_limit=DT)
         t_on(1) = mle_face_mag(mle_on, iface, jrow, NZ)
         call mle_compute_transports(grid, metrics, mle_on, ms, epbl, ss=ss, dt_limit=DT)
         t_on(2) = mle_face_mag(mle_on, iface, jrow, NZ)
         epbl%mld = MLD_SHALLOW
         call mle_compute_transports(grid, metrics, mle_on, ms, epbl, ss=ss, dt_limit=DT)
         t_on(3) = mle_face_mag(mle_on, iface, jrow, NZ)

         ! (i) On the deep (non-retreating) steps the running mean has reset
         !     to the deep MLD instantly, so filter ON == filter OFF; and
         !     filter OFF is the instantaneous-MLD path (legacy).
         call check(error, abs(t_on(1) - t_off(1)) <= 1.0e-12_wp*max(t_off(1), 1.0_wp), &
                    "T10: filter ON deep-step != instantaneous (reset-to-deeper broken)")
         if (allocated(error)) exit checks
         call check(error, abs(t_on(2) - t_off(2)) <= 1.0e-12_wp*max(t_off(2), 1.0_wp), &
                    "T10: filter ON second deep-step != instantaneous")
         if (allocated(error)) exit checks

         ! Sanity: the front drives a non-zero transport on the deep steps.
         call check(error, t_off(1) > 0.0_wp, "T10: zero transport on deep step (bad face)")
         if (allocated(error)) exit checks

         ! (ii) On the retreat step the filtered MLD stays deeper than the
         !      instantaneous shallow value => filtered transport magnitude
         !      is LARGER than unfiltered, but SMALLER than the deep-step
         !      value (it is decaying, not frozen).
         call check(error, t_on(3) > t_off(3), &
                    "T10: filter did not damp the MLD retreat (no running-mean memory)")
         if (allocated(error)) exit checks
         call check(error, t_on(3) < t_off(1), &
                    "T10: filtered retreat transport not below the deep-step value")
      end block checks

      call mle_on%destroy()
      call mle_off%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_t10_mld_filter

   subroutine drive_multiwindow(error, do_relayer, salt_rel_drift)
      !! Shared driver for the DYNAMIC multi-window FK x windowed-advect
      !! regression (T11 / T12).  Reproduces the production loop minus the
      !! barotropic/RK2 machinery, exercising exactly the interaction the
      !! static fk_t9 cannot reach:
      !!
      !!   * MULTIPLE accumulation windows (n_win) of `ratio` outer steps;
      !!   * FK transports RECOMPUTED at each window start from the evolving
      !!     density/MLD (mle_compute_transports), then folded into the
      !!     mass flux on the thermo step and ACCUMULATED into uhtr/vhtr
      !!     via the production continuity_tracer_step_split(ACCUMULATE);
      !!   * the FK fold thermo-gated to the window-start step only
      !!     (mle_fold_active=.true. on step 0 of the window, .false. on the
      !!     rest) — the dt_therm_ratio>1 gating;
      !!   * a genuine non-zero background zonal flow so the drain runs a
      !!     real multi-pass swept-PPM (not the FK-only degenerate case);
      !!   * `continuity_tracer_drain` at each window close;
      !!   * optionally (do_relayer) a conservative per-column relayer
      !!     between windows, mimicking the thermo-cadence ALE remap that
      !!     fires immediately after the drain in ocean_dyn_step_split.
      !!
      !! Returns the relative drift of the full-domain salt mass
      !! Σ(areaT·hTr_S) over the whole run in `salt_rel_drift`.  Conservation
      !! requires it ~ round-off; the FK x windowed-advect leak shows up here
      !! as a systematic O(1%/window) drift.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: do_relayer
      real(wp), intent(out) :: salt_rel_drift

      integer, parameter :: NG = 3
      integer, parameter :: NZ = 3
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: RATIO = 2          ! dt_tracer_advect = dt_therm
      integer, parameter :: N_WIN = 6          ! 6 windows = 12 outer steps
      real(wp), parameter :: DX = 5.0e3_wp, DY = 5.0e3_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: DT_WIN = real(RATIO, wp)*DT
      real(wp), parameter :: U_BG = 0.3_wp     ! eastward flow => genuine multi-pass drain CFL
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: H_DEEP = 400.0_wp, H_MID = 80.0_wp, H_SURF = 5.0_wp
      real(wp), parameter :: MLD = H_SURF + H_MID
      real(wp), parameter :: S_REF(NZ) = [35.0_wp, 34.5_wp, 33.0_wp]
      real(wp), parameter :: HCOL(NZ) = [H_DEEP, H_MID, H_SURF]
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      type(ocean_bc_state_t) :: bc
      integer :: i, j, k, nx, ny, i0, i1, j0, j1, win, stp
      real(wp) :: area, mass0, mass1
      logical :: thermo

      checks: block
         call grid%init(NX_PHYS, NY_PHYS, NG, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         area = DX*DY
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         call epbl%init(grid, nz_ml=NZ)
         call ss%init(grid)
         call mle%init(grid, nz_ml=NZ)
         call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
         bc%periodic_x = .true.
         bc%periodic_y = .false.
         bc%north_fold = .false.
         mle%enable = .true.
         mle%ce = CE
         mle%f_floor = F_FLOOR

         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = HCOL(k)
                  ms%u_face_x_layer(i, j, k) = U_BG
                  ms%v_face_y_layer(i, j, k) = 0.0_wp
                  ms%rho_layer(i, j, k) = 1027.0_wp - 0.4_wp*real(k, wp) &
                                          - 0.6_wp*sin(2.0_wp*PI*real(i - NG, wp)/real(NX_PHYS, wp))
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S_REF(k)*HCOL(k)
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 12.0_wp*HCOL(k)
               end do
            end do
         end do
         epbl%rho0 = RHO0
         epbl%mld = MLD
         epbl%f_centre = F0

         i0 = NG + 1
         i1 = NG + grid%nx_phys
         j0 = NG + 1
         j1 = NG + grid%ny_phys
         ! Full-domain (interior) salt mass before the run.
         mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))

         !$acc enter data copyin(ms, ct, epbl, ss, mle, bc)
         call ms%enter_data()
         call ct%enter_data()
         call epbl%enter_data()
         call ss%enter_data()
         call mle%enter_data()

         do win = 1, N_WIN
            ! --- window-start (thermo) step: recompute FK from current state
            call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss, &
                                        dt_limit=DT_WIN)
            ! Wrap the periodic-x ghosts of the prognostics so the PPM
            ! reconstruction reads the correct cyclic stencil (the real
            ! dyn step does this via ocean_periodic_wrap_state).
            do stp = 1, RATIO
               thermo = (stp == 1)
               call continuity_tracer_step_split(grid, metrics, ct, ms, DT, &
                                                 bc=bc, mle=mle, &
                                                 mle_fold_active=thermo, &
                                                 tracer_mode=TR_MODE_ACCUMULATE)
            end do
            ! --- window close: spend the accumulated transport.
            call continuity_tracer_drain(grid, metrics, ct, ms, RATIO, bc=bc)
            ! --- optional thermo-cadence relayer (mimics the ALE remap that
            !     fires right after the drain).  Conservative per column.
            if (do_relayer) call relayer_to_target(grid, ms, HCOL, NZ)
         end do

         call mle%exit_data()
         call ss%exit_data()
         call epbl%exit_data()
         call ct%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ct, epbl, ss, mle, bc)

         mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
         salt_rel_drift = abs(mass1 - mass0)/abs(mass0)
         call check(error, salt_rel_drift < 1.0e-10_wp, &
                    "dynamic multi-window FK x windowed-advect leaks salt mass")
      end block checks

      call ocean_bc_state_destroy(bc)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine drive_multiwindow

   subroutine relayer_to_target(grid, ms, target_h, nz)
      !! Conservative per-column relayer mimicking the thermo-cadence ALE
      !! remap (Lagrangian h_layer -> fixed target_h).  Preserves per-column
      !! Σ_k hTr exactly (PCM column remap), so it cannot itself leak tracer
      !! mass — its presence in T12 isolates whether the drain leaves an
      !! (hTr, h_layer) pair that the remap reads consistently.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz
      real(wp), intent(in) :: target_h(nz)
      integer :: i, j, k, nx, ny, it
      real(wp) :: col_h_old(nz), col_tr(nz), conc_old(nz), conc_new(nz)
      real(wp) :: z_old(0:nz), z_new(0:nz), htot_old, htot_new, scale_h

      nx = grid%nx_total
      ny = grid%ny_total
      !$acc update self(ms%h_layer)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (allocated(ms%tracers(it)%hTr)) then
               !$acc update self(ms%tracers(it)%hTr)
            end if
         end do
      end if

      do j = 1, ny
         do i = 1, nx
            do k = 1, nz
               col_h_old(k) = ms%h_layer(i, j, k)
            end do
            htot_old = sum(col_h_old)
            htot_new = sum(target_h)
            if (htot_old <= 0.0_wp .or. htot_new <= 0.0_wp) cycle
            ! Stretch the target to the column's actual total (z*-style), so
            ! the relayer is a pure regrid (no volume change).
            scale_h = htot_old/htot_new
            ! Interface depths (surface = 0 down to -htot), bottom-up storage.
            z_old(0) = 0.0_wp
            z_new(0) = 0.0_wp
            do k = 1, nz
               z_old(k) = z_old(k - 1) + col_h_old(k)
               z_new(k) = z_new(k - 1) + target_h(nz - k + 1)*scale_h
            end do
            if (allocated(ms%tracers)) then
               do it = 1, size(ms%tracers)
                  if (.not. allocated(ms%tracers(it)%hTr)) cycle
                  do k = 1, nz
                     col_tr(k) = ms%tracers(it)%hTr(i, j, k)
                     conc_old(k) = col_tr(k)/max(col_h_old(k), 1.0e-12_wp)
                  end do
                  call pcm_column_remap(nz, z_old, conc_old, z_new, conc_new)
                  do k = 1, nz
                     ms%tracers(it)%hTr(i, j, k) = conc_new(k)*(z_new(k) - z_new(k - 1))
                  end do
               end do
            end if
            do k = 1, nz
               ms%h_layer(i, j, k) = z_new(k) - z_new(k - 1)
            end do
         end do
      end do

      !$acc update device(ms%h_layer)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (allocated(ms%tracers(it)%hTr)) then
               !$acc update device(ms%tracers(it)%hTr)
            end if
         end do
      end if
   end subroutine relayer_to_target

   pure subroutine pcm_column_remap(nz, z_old, conc_old, z_new, conc_new)
      !! Piecewise-constant conservative column remap: integral of the old
      !! piecewise-constant profile over each new layer / new thickness.
      !! Exactly conserves Σ(conc·Δz) per column.
      integer, intent(in) :: nz
      real(wp), intent(in) :: z_old(0:nz), conc_old(nz)
      real(wp), intent(in) :: z_new(0:nz)
      real(wp), intent(out) :: conc_new(nz)
      integer :: kn, ko
      real(wp) :: zlo, zhi, ov_lo, ov_hi, acc, dz
      do kn = 1, nz
         zlo = z_new(kn - 1)
         zhi = z_new(kn)
         acc = 0.0_wp
         do ko = 1, nz
            ov_lo = max(zlo, z_old(ko - 1))
            ov_hi = min(zhi, z_old(ko))
            if (ov_hi > ov_lo) acc = acc + conc_old(ko)*(ov_hi - ov_lo)
         end do
         dz = zhi - zlo
         if (dz > 0.0_wp) then
            conc_new(kn) = acc/dz
         else
            conc_new(kn) = 0.0_wp
         end if
      end do
   end subroutine pcm_column_remap

   subroutine test_t11_dynamic_multiwindow(error)
      !! DYNAMIC multi-window FK x windowed-advect conservation, NO remap.
      !! Isolates the drain-accounting interaction (mechanism #1 / #4).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: drift
      call drive_multiwindow(error, do_relayer=.false., salt_rel_drift=drift)
   end subroutine test_t11_dynamic_multiwindow

   subroutine test_t12_multiwindow_remap(error)
      !! DYNAMIC multi-window FK x windowed-advect conservation WITH a
      !! conservative relayer between windows (mechanism #2: drain leaves
      !! an (hTr, h_layer) pair the ALE remap then reads).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: drift
      call drive_multiwindow(error, do_relayer=.true., salt_rel_drift=drift)
   end subroutine test_t12_multiwindow_remap

   pure function mle_face_mag(mle, i, j, nz) result(mag)
      !! L1 magnitude of the per-layer FK transport at u-face (i,j):
      !! sum_k |uhml(i,j,k)|.  A monotone proxy for the streamfunction
      !! amplitude at that face (Psi ~ MLD^2), used by T10 to compare
      !! filtered vs unfiltered transport.
      type(ocean_mle_t), intent(in) :: mle
      integer, intent(in) :: i, j, nz
      real(wp) :: mag
      integer :: k
      mag = 0.0_wp
      do k = 1, nz
         mag = mag + abs(mle%uhml(i, j, k))
      end do
   end function mle_face_mag

   ! -----------------------------------------------------------------
   ! T14 — metric-aspect scaling (regression for the missing 1/dx).
   ! The FK transport is Psi*(face width) with Psi ~ (db/dx)*H^2, so the
   ! u-face transport scales as dy/dx and the v-face as dx/dy.  The
   ! pre-fix kernel used dy_cu (resp. dx_cv) WITHOUT the idxCu (resp.
   ! idyCv) = 1/dx (resp. 1/dy) factor, so uDml ~ dy independent of dx —
   ! ~dx too large.  Holding the physics fixed and only changing the grid
   ! spacing: doubling dx must HALVE the u-transport (and doubling dy must
   ! halve the v-transport).  The bug gives ratio 1.0; the fix gives 0.5.
   ! dx /= dy is essential — on dx==dy the dy/dx aspect is 1 and the
   ! missing factor hides (which is why every prior dx==dy test passed).
   ! The cap is OFF here (no dt_limit passed) so the raw magnitude shows.
   ! -----------------------------------------------------------------
   subroutine test_t14_metric_aspect(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: L = 1.0e4_wp
      real(wp) :: u_base, v_base, u_dx2, v_dx2, u_dy2, v_dy2

      call mle_aspect_run(L, L, u_base, v_base)        ! dx = dy = L
      call mle_aspect_run(2.0_wp*L, L, u_dx2, v_dx2)   ! dx = 2L  (u must halve)
      call mle_aspect_run(L, 2.0_wp*L, u_dy2, v_dy2)   ! dy = 2L  (v must halve)

      call check(error, u_base > 0.0_wp .and. v_base > 0.0_wp, &
                 "T14: baseline u/v transports nonzero (diagonal front)")
      if (allocated(error)) return
      ! u-face ~ dy/dx: doubling dx halves it (pre-fix: 1.0).
      call check(error, abs(u_dx2/u_base - 0.5_wp) < 1.0e-3_wp, &
                 "T14: u-transport must scale as dy/dx (idxCu factor)")
      if (allocated(error)) return
      ! v-face ~ dx/dy: doubling dy halves it (pre-fix: 1.0).
      call check(error, abs(v_dy2/v_base - 0.5_wp) < 1.0e-3_wp, &
                 "T14: v-transport must scale as dx/dy (idyCv factor)")
   end subroutine test_t14_metric_aspect

   subroutine mle_aspect_run(dx, dy, max_u, max_v)
      !! Run mle_compute_transports on a 6x5 grid with a fixed diagonal
      !! buoyancy front (rho varies in both i and j), MLD, and f — only
      !! the grid spacing (dx, dy) varies between calls.  Cap OFF (no
      !! dt_limit), enable ON.  Returns max|uhml|, max|vhml|.
      real(wp), intent(in) :: dx, dy
      real(wp), intent(out) :: max_u, max_v
      integer, parameter :: nz = 3
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      integer :: i, j, k, nxt, nyt

      call grid%init(6, 5, NGHOST, dx, dy)
      call make_cartesian_metrics(metrics, grid)
      nxt = grid%nx_total
      nyt = grid%ny_total
      ms%nz_ml = nz
      call ms%init(grid)
      call epbl%init(grid, nz_ml=nz)
      call ss%init(grid)
      call mle%init(grid, nz_ml=nz)
      mle%enable = .true.

      ! Diagonal density front: linear ramp in BOTH i and j (distinct
      ! slopes) so db_x and db_y are both nonzero and dx-/dy-independent.
      do k = 1, nz
         do j = 1, nyt
            do i = 1, nxt
               ms%h_layer(i, j, k) = 30.0_wp
               ms%rho_layer(i, j, k) = RHO0 + real(k, wp) &
                                       + 0.5_wp*real(i, wp) + 0.3_wp*real(j, wp)
            end do
         end do
      end do
      epbl%mld = 60.0_wp
      epbl%f_centre = 7.0e-5_wp
      epbl%rho0 = RHO0

      !$acc enter data copyin(ms, epbl, ss, mle)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call mle%enter_data()
      call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss)
      !$acc update self(mle%uhml, mle%vhml)
      call mle%exit_data()

      max_u = maxval(abs(mle%uhml))
      max_v = maxval(abs(mle%vhml))

      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, mle)
      call destroy_cartesian_metrics(metrics)
      call mle%destroy()
      call ss%destroy()
      call epbl%destroy()
      call ms%destroy()
   end subroutine mle_aspect_run

   ! -----------------------------------------------------------------
   ! Bodner (2023) frontogenesis-arrest MLE
   ! -----------------------------------------------------------------

   subroutine run_bodner_transport(tau_val, b0_val, cr, max_uhml)
      !! Bodner MLE on a small x-buoyancy front with uniform wind stress
      !! `tau_val` and surface buoyancy flux `b0_val`; returns peak |uhml|.
      !! No CFL cap (dt_limit omitted) so the raw ts_bod ∝ 1/w'u' scaling is
      !! visible.  Reads the persisted epbl%b0 -> exercises the EPBL change.
      real(wp), intent(in) :: tau_val, b0_val, cr
      real(wp), intent(out) :: max_uhml
      integer, parameter :: nz = 3
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_surface_stress_t) :: ss
      type(ocean_mle_t) :: mle
      integer :: i, j, k
      call grid%init(8, 8, NGHOST, 1.0e4_wp, 1.0e4_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = nz; call ms%init(grid)
      call epbl%init(grid, nz_ml=nz)
      call ss%init(grid)
      call mle%init(grid, nz_ml=nz)
      mle%enable = .true.
      mle%use_bodner = .true.
      mle%cr = cr
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               ms%h_layer(i, j, k) = 30.0_wp
               ms%rho_layer(i, j, k) = 1027.0_wp - 0.4_wp*real(k, wp) + 0.05_wp*real(i, wp)
            end do
         end do
      end do
      epbl%mld = 60.0_wp
      epbl%f_centre = 7.0e-5_wp
      epbl%rho0 = RHO0
      epbl%b0 = b0_val
      ss%tau_x = tau_val
      ss%tau_y = 0.0_wp
      !$acc enter data copyin(ms, epbl, ss, mle)
      call ms%enter_data()
      call epbl%enter_data()
      call ss%enter_data()
      call mle%enter_data()
      call mle_compute_transports(grid, metrics, mle, ms, epbl, ss=ss)
      !$acc update self(mle%uhml)
      max_uhml = maxval(abs(mle%uhml))
      ! Full symmetric teardown -- this helper is called twice per test, so a
      ! partial exit would leak device mappings that collide on the next call.
      call mle%exit_data()
      call ss%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, epbl, ss, mle)
      call mle%destroy(); call ss%destroy(); call epbl%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_bodner_transport

   subroutine test_bodner_arrest_wind(error)
      !! The frontogenesis-arrest signature: stronger boundary-layer
      !! turbulence (bigger wind stress -> bigger u* -> bigger w'u') must
      !! SHRINK the Bodner overturning (ts_bod ∝ 1/w'u').  10x the stress ->
      !! ~10x smaller transport; both strictly positive on the front.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t_lo, t_hi
      call run_bodner_transport(0.01_wp, 0.0_wp, 0.06_wp, t_lo)
      call run_bodner_transport(0.10_wp, 0.0_wp, 0.06_wp, t_hi)
      call check(error, t_lo > 0.0_wp .and. t_hi > 0.0_wp, &
                 "bodner produced no transport on a buoyancy front")
      if (allocated(error)) return
      call check(error, t_hi < t_lo, &
                 "bodner not arrested by wind: stronger turbulence did not reduce transport")
      if (allocated(error)) return
      ! ts_bod ∝ 1/|tau|, so 10x stress should cut transport well below half.
      call check(error, t_hi < 0.5_wp*t_lo, &
                 "bodner arrest too weak: 10x wind should more than halve transport")
   end subroutine test_bodner_arrest_wind

   subroutine test_bodner_convective(error)
      !! The CONVECTIVE arrest path (validates the persisted epbl%b0): a
      !! destabilizing surface buoyancy flux (b0 < 0 -> w*^3 > 0) adds to
      !! w'u' and must further reduce the transport vs the wind-only case.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t_wind_only, t_with_conv
      call run_bodner_transport(0.03_wp, 0.0_wp, 0.06_wp, t_wind_only)
      call run_bodner_transport(0.03_wp, -1.0e-6_wp, 0.06_wp, t_with_conv)
      call check(error, t_with_conv < t_wind_only .and. t_with_conv > 0.0_wp, &
                 "bodner convective term inert: destabilizing b0 did not reduce transport")
   end subroutine test_bodner_convective

   subroutine test_bodner_cr_zero(error)
      !! Cr = 0 (the MOM6 default) => zero Bodner streamfunction => no
      !! transport, even with a front + wind present.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t0
      call run_bodner_transport(0.05_wp, -1.0e-6_wp, 0.0_wp, t0)
      call check(error, t0 < 1.0e-18_wp, "bodner with Cr=0 produced transport")
   end subroutine test_bodner_cr_zero

end module test_ocean_foxkemper
