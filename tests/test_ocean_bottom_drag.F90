!! Unit tests for the ocean bottom-drag kernel
!! (rdb_ocean_bottom_drag).  Two variants: linear Rayleigh and
!! quadratic log-layer.  Drag acts only on the bed-most layer
!! (k = 1 under the ROMS-style convention).
!!
!! Cases:
!!   * Linear decay rate — uniform u at k=1, k>=2 zero.  After N
!!     steps the bed-layer velocity must equal (1 - dt*r)^N * U0,
!!     and the layers above must remain at zero.  Pinpoint check
!!     on the analytic rate constant.
!!   * Bed-only support — non-trivial multilayer u, v.  After
!!     one drag step every layer EXCEPT k=1 must remain at its
!!     IC.  Catches accidental whole-column writes.
!!   * Zero-coefficient short-circuit — `r_linear = 0` AND
!!     `c_drag = 0` is a no-op.  Lets the driver wire the kernel
!!     unconditionally and gate via the namelist.
!!   * Quadratic stability — Non-trivial multilayer u, v with
!!     bottom-layer flow.  After several steps the bed-layer KE
!!     must decrease and stay finite (the `h_min` floor keeps
!!     the kernel bounded under vanishing layers).
module test_ocean_bottom_drag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, &
                                    BDRAG_LINEAR, BDRAG_QUADRATIC, &
                                    ocean_bottom_drag_compute_tendencies, &
                                    ocean_bottom_drag_apply_tendencies, &
                                    ocean_channel_drag_compute_tendencies, &
                                    ocean_channel_drag_apply_tendencies
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_finalize
   implicit none
   private

   public :: collect_ocean_bottom_drag_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_bottom_drag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("linear_decay_rate", test_linear_decay), &
                  new_unittest("bed_only_support", test_bed_only), &
                  new_unittest("zero_coeff_no_op", test_zero_coeff), &
                  new_unittest("quadratic_KE_monotone", test_quadratic_ke), &
                  new_unittest("bed_factor_unity_bit_identical", test_bed_factor_unity), &
                  new_unittest("bed_factor_scales_bed_only", test_bed_factor_scales_bed), &
                  new_unittest("channel_drag_off_no_op", test_channel_drag_off), &
                  new_unittest("channel_drag_allwet_no_op", test_channel_drag_allwet), &
                  new_unittest("channel_drag_rate_closed_form", test_channel_drag_rate), &
                  new_unittest("channel_drag_implicit_update", test_channel_drag_implicit), &
                  new_unittest("implicit_thin_layer_stable", test_implicit_thin_layer) &
                  ]
   end subroutine collect_ocean_bottom_drag_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, bd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_bottom_drag_t), intent(inout) :: bd
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(bd)
      call bd%enter_data()
   end subroutine map_in

   subroutine map_out(ms, bd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_bottom_drag_t), intent(inout) :: bd
      call bd%exit_data()
      !$acc exit data delete(bd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_linear_decay(error)
      !! Uniform u = U0 at k=1, zero elsewhere.  Linear drag with
      !! coefficient r predicts u_N = (1 - dt*r)^N * U0 in the
      !! bottom layer (away from walls where the kernel forces
      !! zero tendency).  After N steps the realised value must
      !! match the analytic to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: U0 = 0.4_wp
      real(wp), parameter :: R = 0.1_wp
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_STEPS = 12
      integer :: i, j, k, step, nx, ny
      real(wp) :: u_expected, u_obs, max_upper
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         bd%variant = BDRAG_LINEAR
         bd%r_linear = R
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               ms%u_face_x_layer(i, j, 1) = U0
            end do
         end do

         call map_in(ms, bd)
         do step = 1, N_STEPS
            call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
            call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
         end do
         call map_out(ms, bd)

         u_expected = U0*(1.0_wp - DT*R)**N_STEPS
         u_obs = ms%u_face_x_layer(nx/2, ny/2, 1)
         max_upper = max(maxval(abs(ms%u_face_x_layer(:, :, 2:))), &
                         maxval(abs(ms%v_face_y_layer(:, :, 2:))))

         call check(error, abs(u_obs - u_expected) < 1.0e-12_wp, &
                    "linear decay: bed-layer u off analytic")
         if (allocated(error)) exit checks
         call check(error, max_upper < 1.0e-12_wp, &
                    "linear decay: drag leaked into layers k>=2")

      end block checks
      call bd%destroy(); call ms%destroy()
   end subroutine test_linear_decay

   subroutine test_bed_only(error)
      !! Non-trivial multilayer u, v.  After one drag step, every
      !! layer EXCEPT k=1 must be byte-identical to the IC.  Guards
      !! against accidentally writing to k>=2.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      real(wp) :: max_du_upper, max_dv_upper
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         bd%variant = BDRAG_LINEAR
         bd%r_linear = 0.5_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i + k, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = cos(2.0_wp*PI*real(j + k, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, bd)
         call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
         call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
         call map_out(ms, bd)

         max_du_upper = maxval(abs(ms%u_face_x_layer(:, :, 2:) - u_ic(:, :, 2:)))
         max_dv_upper = maxval(abs(ms%v_face_y_layer(:, :, 2:) - v_ic(:, :, 2:)))

         call check(error, max_du_upper < 1.0e-12_wp, &
                    "bed-only: u in layers k>=2 was touched")
         if (allocated(error)) exit checks
         call check(error, max_dv_upper < 1.0e-12_wp, &
                    "bed-only: v in layers k>=2 was touched")

      end block checks
      deallocate (u_ic, v_ic)
      call bd%destroy(); call ms%destroy()
   end subroutine test_bed_only

   subroutine test_zero_coeff(error)
      !! `r_linear = c_drag = 0` → tendency stays zero regardless
      !! of variant; apply is a no-op.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         bd%variant = BDRAG_QUADRATIC
         bd%r_linear = 0.0_wp
         bd%c_drag = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = 0.2_wp*sin(PI*real(i, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = 0.1_wp*cos(PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, bd)
         call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
         call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
         call map_out(ms, bd)

         call check(error, maxval(abs(ms%u_face_x_layer - u_ic)) < 1.0e-12_wp, &
                    "zero coeff: u changed")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%v_face_y_layer - v_ic)) < 1.0e-12_wp, &
                    "zero coeff: v changed")

      end block checks
      deallocate (u_ic, v_ic)
      call bd%destroy(); call ms%destroy()
   end subroutine test_zero_coeff

   subroutine test_quadratic_ke(error)
      !! Quadratic-drag closure on a non-trivial bottom-layer flow.
      !! Bed-layer KE must decrease step-over-step (drag is
      !! dissipative; the apply step removes energy without
      !! re-injecting any).  Tail check: the field stays finite —
      !! `h_min` keeps the `u/h_bot` division well-behaved.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.1_wp
      integer, parameter :: N_STEPS = 6
      integer :: i, j, k, step, nx, ny
      real(wp) :: ke_prev, ke_now, ke_initial, max_finite
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         bd%variant = BDRAG_QUADRATIC
         bd%c_drag = 2.5e-3_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, 1) = 0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, 1) = 0.3_wp*cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
            end do
         end do

         ke_initial = sum(ms%u_face_x_layer(:, :, 1)**2) + sum(ms%v_face_y_layer(:, :, 1)**2)
         ke_prev = ke_initial
         call map_in(ms, bd)
         do step = 1, N_STEPS
            call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
            call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
            !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
            ke_now = sum(ms%u_face_x_layer(:, :, 1)**2) + sum(ms%v_face_y_layer(:, :, 1)**2)
            if (ke_now >= ke_prev + 1.0e-14_wp) exit
            ke_prev = ke_now
         end do
         call map_out(ms, bd)

         max_finite = maxval(abs(ms%u_face_x_layer)) + maxval(abs(ms%v_face_y_layer))

         call check(error, ke_now < ke_prev + 1.0e-12_wp, &
                    "quadratic drag: bed-KE not monotone decreasing")
         if (allocated(error)) exit checks
         call check(error, ke_now < ke_initial, &
                    "quadratic drag: no net KE decrease")
         if (allocated(error)) exit checks
         call check(error, max_finite < 1.0e6_wp, &
                    "quadratic drag: field went non-finite")

      end block checks
      call bd%destroy(); call ms%destroy()
   end subroutine test_quadratic_ke

   ! -----------------------------------------------------------------
   ! bed_factor (layer-specific drag) cases
   ! -----------------------------------------------------------------
   !
   ! `bdrag%bed_factor` multiplies the bed-layer (k=1) drag tendency
   ! only — layers k>=2 stay at the nominal rate.  Default 1.0 keeps
   ! the historical drag bit-identical; values > 1 strengthen bed
   ! damping without altering the rest of the BBL.

   subroutine test_bed_factor_unity(error)
      !! `bed_factor = 1.0` (default) must give the same drag tendency
      !! as not setting the knob at all.  Tests the no-op guarantee:
      !! existing namelists that don't mention `bed_factor` get the
      !! pre-knob behaviour to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_bottom_drag_t) :: bd_a, bd_b
      real(wp), parameter :: U0 = 0.4_wp
      real(wp), parameter :: R = 0.1_wp
      real(wp), parameter :: HBBL = 50.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp) :: max_diff
      integer :: i, j, nx, ny
      checks: block
         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid); call ms_b%init(grid)
         call bd_a%init(grid, nz_ml=NZ); call bd_b%init(grid, nz_ml=NZ)
         bd_a%variant = BDRAG_LINEAR; bd_a%r_linear = R; bd_a%hbbl = HBBL
         bd_b%variant = BDRAG_LINEAR; bd_b%r_linear = R; bd_b%hbbl = HBBL
         ! Engage the HBBL-distributed path (this is where bed_factor lives).
         bd_a%bed_factor = 1.0_wp
         ! bd_b%bed_factor untouched ⇒ uses its type-default 1.0
         nx = grid%nx_total; ny = grid%ny_total

         ms_a%h_layer = 10.0_wp; ms_b%h_layer = 10.0_wp
         ms_a%u_face_x_layer = 0.0_wp; ms_b%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               ms_a%u_face_x_layer(i, j, 1) = U0
               ms_b%u_face_x_layer(i, j, 1) = U0
            end do
         end do

         call map_in(ms_a, bd_a); call map_in(ms_b, bd_b)
         call ocean_bottom_drag_compute_tendencies(grid, bd_a, ms_a, DT)
         call ocean_bottom_drag_apply_tendencies(bd_a, ms_a, DT)
         call ocean_bottom_drag_compute_tendencies(grid, bd_b, ms_b, DT)
         call ocean_bottom_drag_apply_tendencies(bd_b, ms_b, DT)
         call map_out(ms_a, bd_a); call map_out(ms_b, bd_b)

         max_diff = max(maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer)), &
                        maxval(abs(ms_a%v_face_y_layer - ms_b%v_face_y_layer)))
         call check(error, max_diff < 1.0e-14_wp, &
                    "bed_factor=1.0 must reproduce default-knob path bit-identically")
      end block checks
      call bd_a%destroy(); call bd_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_bed_factor_unity

   subroutine test_bed_factor_scales_bed(error)
      !! With h_face_k > hbbl on the bed layer alone, the linear
      !! distributed drag reduces to `du/dt = -r · (hbbl/h_face) · u`
      !! on k=1 (bed gets full BBL coverage, k>=2 is outside the band
      !! and stays at zero tendency).  Scaling by `bed_factor = 3`
      !! multiplies the bed tendency by 3 exactly; nothing else.
      !!
      !! Two runs from identical IC: `bf=1` and `bf=3`.  After one
      !! step the relative |Δu| change on the bed must be
      !! `3·R·DT·(hbbl/h_face)` — and layers k>=2 must remain at U0
      !! (drag never reaches them through the BBL band) so the diff
      !! is exactly zero there.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_one, ms_three
      type(ocean_bottom_drag_t) :: bd_one, bd_three
      real(wp), parameter :: U0 = 0.4_wp
      real(wp), parameter :: R = 0.1_wp
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: HBBL = 8.0_wp      ! < H_LAYER ⇒ bed-only BBL
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: BF = 3.0_wp
      real(wp) :: du_bed_one, du_bed_three, max_upper_diff
      real(wp) :: expected_ratio, observed_ratio
      integer :: i, j, k, nx, ny
      checks: block
         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms_one%nz_ml = NZ; ms_three%nz_ml = NZ
         call ms_one%init(grid); call ms_three%init(grid)
         call bd_one%init(grid, nz_ml=NZ); call bd_three%init(grid, nz_ml=NZ)
         bd_one%variant = BDRAG_LINEAR; bd_one%r_linear = R; bd_one%hbbl = HBBL
         bd_three%variant = BDRAG_LINEAR; bd_three%r_linear = R; bd_three%hbbl = HBBL
         bd_one%bed_factor = 1.0_wp
         bd_three%bed_factor = BF
         nx = grid%nx_total; ny = grid%ny_total

         ms_one%h_layer = H_LAYER; ms_three%h_layer = H_LAYER
         ms_one%u_face_x_layer = 0.0_wp; ms_three%u_face_x_layer = 0.0_wp
         ms_one%v_face_y_layer = 0.0_wp; ms_three%v_face_y_layer = 0.0_wp
         ! Fill every layer with U0 so we can also check that k>=2
         ! is unchanged in BOTH runs.
         do k = 1, NZ
            do j = 1, ny
               do i = 2, nx
                  ms_one%u_face_x_layer(i, j, k) = U0
                  ms_three%u_face_x_layer(i, j, k) = U0
               end do
            end do
         end do

         call map_in(ms_one, bd_one); call map_in(ms_three, bd_three)
         call ocean_bottom_drag_compute_tendencies(grid, bd_one, ms_one, DT)
         call ocean_bottom_drag_apply_tendencies(bd_one, ms_one, DT)
         call ocean_bottom_drag_compute_tendencies(grid, bd_three, ms_three, DT)
         call ocean_bottom_drag_apply_tendencies(bd_three, ms_three, DT)
         call map_out(ms_one, bd_one); call map_out(ms_three, bd_three)

         ! Bed-layer Δu sample at an interior face
         du_bed_one = U0 - ms_one%u_face_x_layer(nx/2, ny/2, 1)
         du_bed_three = U0 - ms_three%u_face_x_layer(nx/2, ny/2, 1)
         observed_ratio = du_bed_three/du_bed_one
         expected_ratio = BF

         call check(error, abs(observed_ratio - expected_ratio) < 1.0e-10_wp, &
                    "bed_factor: Δu(bed) ratio must equal bed_factor")
         if (allocated(error)) exit checks

         ! Layers k>=2 should be untouched by drag in BOTH runs — the
         ! BBL band stops at hbbl < h_face_k of layer 1, so layers
         ! k>=2 never enter the drag loop.  Compare run-vs-run rather
         ! than against U0 (the IC fill leaves wall faces at 0 which
         ! would confuse an absolute-equality check).  If bed_factor
         ! incorrectly leaked into upper layers, the bf=1 vs bf=3 runs
         ! would diverge at k>=2; they don't.
         max_upper_diff = maxval(abs(ms_one%u_face_x_layer(:, :, 2:) - &
                                     ms_three%u_face_x_layer(:, :, 2:)))
         call check(error, max_upper_diff < 1.0e-12_wp, &
                    "bed_factor: must not affect layers k>=2")
      end block checks
      call bd_one%destroy(); call bd_three%destroy()
      call ms_one%destroy(); call ms_three%destroy()
   end subroutine test_bed_factor_scales_bed

   ! -----------------------------------------------------------------
   ! Channel (side-wall) drag cases
   ! -----------------------------------------------------------------
   !
   ! `bdrag%channel_drag` adds a per-layer lateral Rayleigh drag with
   ! rate `lambda = cdrag_side * |U_face| * f_blocked / cell_width`,
   ! applied implicitly `u <- u/(1+dt*lambda)`.  f_blocked is the
   ! fraction of the cross-stream perimeter blocked by land (`wet_q==0`)
   ! or by a vanished neighbour layer.  All-wet / flat-bottom ⇒
   ! f_blocked==0 ⇒ exact no-op; default-off short-circuits to zero.

   subroutine make_metrics(metrics, grid)
      !! Cartesian metrics WITHOUT the device map yet — caller may edit
      !! wet_q / face widths before mapping.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      call metrics%init(grid)
      call metrics_fill_cartesian(metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(metrics)
   end subroutine make_metrics

   subroutine map_in_cd(ms, bd, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_metrics_t), intent(inout) :: metrics
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(bd)
      call bd%enter_data()
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine map_in_cd

   subroutine map_out_cd(ms, bd, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_metrics_t), intent(inout) :: metrics
      call metrics%exit_data()
      !$acc exit data delete(metrics)
      call bd%exit_data()
      !$acc exit data delete(bd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_cd

   subroutine test_channel_drag_off(error)
      !! `channel_drag = .false.` (default) ⇒ even with land present and
      !! a non-trivial flow, the side-drag apply leaves u/v unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      integer :: nx, ny
      checks: block
         call make_grid(grid, 12, 10, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         call make_metrics(metrics, grid)
         ! channel_drag left .false. (default), but give a real coeff +
         ! a land block to prove the GATE not the coefficient is what
         ! short-circuits.
         bd%cdrag_side = 2.5e-3_wp
         nx = grid%nx_total; ny = grid%ny_total
         metrics%wet_q(nx/2, ny/2) = 0.0_wp   ! a blocked corner
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.5_wp
         ms%v_face_y_layer = 0.3_wp
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in_cd(ms, bd, metrics)
         call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
         call ocean_channel_drag_apply_tendencies(bd, ms, DT)
         call map_out_cd(ms, bd, metrics)

         call check(error, maxval(abs(ms%u_face_x_layer - u_ic)) < 1.0e-14_wp, &
                    "channel_drag off: u changed")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%v_face_y_layer - v_ic)) < 1.0e-14_wp, &
                    "channel_drag off: v changed")
      end block checks
      if (allocated(u_ic)) deallocate (u_ic, v_ic)
      call metrics%destroy(); call bd%destroy(); call ms%destroy()
   end subroutine test_channel_drag_off

   subroutine test_channel_drag_allwet(error)
      !! `channel_drag = .true.` but ALL-WET + flat-bottom ⇒ every
      !! `wet_q == 1` and no layer vanishes ⇒ f_blocked == 0 ⇒ exact
      !! no-op.  The bit-identity guarantee for the realistic-bathy knob
      !! on a clean (no-land) domain.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 100.0_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      checks: block
         call make_grid(grid, 12, 10, 1000.0_wp, 1000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         call make_metrics(metrics, grid)
         bd%channel_drag = .true.
         bd%cdrag_side = 2.5e-3_wp
         ! metrics%wet_q stays 1.0 everywhere (init default); h all thick.
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.5_wp
         ms%v_face_y_layer = 0.3_wp
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in_cd(ms, bd, metrics)
         call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
         call ocean_channel_drag_apply_tendencies(bd, ms, DT)
         call map_out_cd(ms, bd, metrics)

         call check(error, maxval(abs(ms%u_face_x_layer - u_ic)) < 1.0e-14_wp, &
                    "channel_drag all-wet: u changed (should be zero blockage)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%v_face_y_layer - v_ic)) < 1.0e-14_wp, &
                    "channel_drag all-wet: v changed")
      end block checks
      if (allocated(u_ic)) deallocate (u_ic, v_ic)
      call metrics%destroy(); call bd%destroy(); call ms%destroy()
   end subroutine test_channel_drag_allwet

   subroutine test_channel_drag_rate(error)
      !! Analytic blocked-perimeter rate.  Block BOTH corners flanking a
      !! single interior u-face (south corner (i,j) and north corner
      !! (i,j+1)) by setting `wet_q == 0` there ⇒ f_blocked == 1.  With
      !! v == 0 everywhere the face speed is |u|, so the closed-form rate
      !! is `lambda = cdrag_side * |u| * 1 / dyCu`.  The implicit update
      !! gives `u_new = u/(1 + dt*lambda)`; back out lambda and compare.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: CD = 2.5e-3_wp
      real(wp), parameter :: DX = 1000.0_wp, DY = 800.0_wp
      real(wp), parameter :: DT = 60.0_wp
      integer :: ii, jj, nx, ny, k
      real(wp) :: lambda_expected, u_new, lambda_obs
      checks: block
         call make_grid(grid, 12, 10, DX, DY)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         call make_metrics(metrics, grid)
         bd%channel_drag = .true.
         bd%cdrag_side = CD
         nx = grid%nx_total; ny = grid%ny_total
         ii = nx/2; jj = ny/2
         ! Block both flanking corners of u-face (ii,jj): wet_q(ii,jj) and
         ! wet_q(ii,jj+1).  All other corners stay wet.
         metrics%wet_q(ii, jj) = 0.0_wp
         metrics%wet_q(ii, jj + 1) = 0.0_wp
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%u_face_x_layer(ii, jj, k) = U0
         end do

         call map_in_cd(ms, bd, metrics)
         call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
         call ocean_channel_drag_apply_tendencies(bd, ms, DT)
         call map_out_cd(ms, bd, metrics)

         ! f_blocked = 1, v_at_u = 0 ⇒ speed = U0, width = dyCu = DY.
         lambda_expected = CD*U0*1.0_wp/DY
         u_new = ms%u_face_x_layer(ii, jj, 1)
         ! Back out lambda from the implicit update: u_new = U0/(1+dt*lam).
         lambda_obs = (U0/u_new - 1.0_wp)/DT

         call check(error, abs(lambda_obs - lambda_expected) < 1.0e-12_wp, &
                    "channel_drag: per-layer rate off closed form")
         if (allocated(error)) exit checks
         ! Same on every layer k (per-layer, not bed-only).
         call check(error, abs(ms%u_face_x_layer(ii, jj, NZ) - u_new) < 1.0e-14_wp, &
                    "channel_drag: surface layer rate differs from bed")
         if (allocated(error)) exit checks
         ! A wet-corner face nearby is untouched.
         call check(error, abs(ms%u_face_x_layer(ii + 2, jj, 1)) < 1.0e-14_wp, &
                    "channel_drag: leaked to an unblocked face")
      end block checks
      call metrics%destroy(); call bd%destroy(); call ms%destroy()
   end subroutine test_channel_drag_rate

   subroutine test_channel_drag_implicit(error)
      !! One corner blocked ⇒ f_blocked = 0.5.  Verify the implicit
      !! update form exactly: with v == 0, u_new == u/(1+dt*lambda) where
      !! lambda = cd*|u|*0.5/dyCu.  Distinct from the rate test: this
      !! pins the 0.5 half-perimeter weight + the implicit denominator.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: U0 = 0.4_wp
      real(wp), parameter :: CD = 3.0e-3_wp
      real(wp), parameter :: DX = 1000.0_wp, DY = 500.0_wp
      real(wp), parameter :: DT = 120.0_wp
      integer :: ii, jj, nx, ny, k
      real(wp) :: lambda, u_expected, u_obs
      checks: block
         call make_grid(grid, 12, 10, DX, DY)
         ms%nz_ml = NZ
         call ms%init(grid)
         call bd%init(grid, nz_ml=NZ)
         call make_metrics(metrics, grid)
         bd%channel_drag = .true.
         bd%cdrag_side = CD
         nx = grid%nx_total; ny = grid%ny_total
         ii = nx/2; jj = ny/2
         metrics%wet_q(ii, jj) = 0.0_wp   ! south corner only
         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%u_face_x_layer(ii, jj, k) = U0
         end do

         call map_in_cd(ms, bd, metrics)
         call ocean_channel_drag_compute_tendencies(grid, metrics, bd, ms)
         call ocean_channel_drag_apply_tendencies(bd, ms, DT)
         call map_out_cd(ms, bd, metrics)

         lambda = CD*U0*0.5_wp/DY
         u_expected = U0/(1.0_wp + DT*lambda)
         u_obs = ms%u_face_x_layer(ii, jj, 1)
         call check(error, abs(u_obs - u_expected) < 1.0e-12_wp, &
                    "channel_drag: implicit update off u/(1+dt*lambda)")
      end block checks
      call metrics%destroy(); call bd%destroy(); call ms%destroy()
   end subroutine test_channel_drag_implicit

   subroutine test_implicit_thin_layer(error)
      !! Implicit (backward-Euler) bottom drag is unconditionally stable on
      !! a THIN bottom layer where the explicit form overshoots.  Setup: bed
      !! layer h=5 m, |U|=8 m/s, C_d=2.5e-3, dt=300 ⇒ explicit drag CFL
      !! `dt·C_d·|U|/h = 1.2 > 1`, so one explicit step flips the sign
      !! (u → U0·(1−1.2) = −0.2·U0, the runaway signature).  The implicit
      !! step must instead give the closed form `u = U0/(1+dt·C_d·|U|/h)`
      !! (monotone decay, same sign, |u|<|U0|) to round-off.  Bed-only mode
      !! (hbbl=0) for a clean single-face analytic.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_i, ms_e
      type(ocean_bottom_drag_t) :: bd_i, bd_e
      real(wp), parameter :: U0 = 8.0_wp, CD = 2.5e-3_wp, HBOT = 5.0_wp
      real(wp), parameter :: DT = 300.0_wp
      integer :: i, j, nx, ny
      real(wp) :: lam, u_imp_exact, u_obs_i, u_obs_e
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ! ---- implicit run ----
         ms_i%nz_ml = NZ
         call ms_i%init(grid)
         call bd_i%init(grid, nz_ml=NZ)
         bd_i%variant = BDRAG_QUADRATIC
         bd_i%c_drag = CD
         bd_i%implicit = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         ms_i%h_layer = HBOT
         ms_i%u_face_x_layer = 0.0_wp
         ms_i%v_face_y_layer = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               ms_i%u_face_x_layer(i, j, 1) = U0
            end do
         end do
         call map_in(ms_i, bd_i)
         call ocean_bottom_drag_compute_tendencies(grid, bd_i, ms_i, DT)
         call ocean_bottom_drag_apply_tendencies(bd_i, ms_i, DT)
         call map_out(ms_i, bd_i)
         u_obs_i = ms_i%u_face_x_layer(nx/2, ny/2, 1)

         ! ---- explicit run (same state) ----
         ms_e%nz_ml = NZ
         call ms_e%init(grid)
         call bd_e%init(grid, nz_ml=NZ)
         bd_e%variant = BDRAG_QUADRATIC
         bd_e%c_drag = CD
         bd_e%implicit = .false.
         ms_e%h_layer = HBOT
         ms_e%u_face_x_layer = 0.0_wp
         ms_e%v_face_y_layer = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               ms_e%u_face_x_layer(i, j, 1) = U0
            end do
         end do
         call map_in(ms_e, bd_e)
         call ocean_bottom_drag_compute_tendencies(grid, bd_e, ms_e, DT)
         call ocean_bottom_drag_apply_tendencies(bd_e, ms_e, DT)
         call map_out(ms_e, bd_e)
         u_obs_e = ms_e%u_face_x_layer(nx/2, ny/2, 1)

         lam = CD*U0/HBOT                       ! drag rate (speed=|U|=U0, v=0)
         u_imp_exact = U0/(1.0_wp + DT*lam)     ! backward-Euler closed form

         ! Implicit matches the closed form to round-off.
         call check(error, abs(u_obs_i - u_imp_exact) < 1.0e-12_wp, &
                    "implicit drag /= u/(1+dt*Cd*|U|/h)")
         if (allocated(error)) exit checks
         ! Implicit is a monotone decay: same sign, magnitude reduced.
         call check(error, u_obs_i > 0.0_wp .and. u_obs_i < U0, &
                    "implicit drag overshot / changed sign")
         if (allocated(error)) exit checks
         ! Explicit, on the same thin-layer state, overshoots to NEGATIVE u
         ! (the conditional-instability signature) — proving the knob acts.
         call check(error, u_obs_e < 0.0_wp, &
                    "explicit drag did not overshoot (test setup below CFL)")
         if (allocated(error)) exit checks
         call check(error, abs(u_obs_i - u_obs_e) > 0.1_wp*U0, &
                    "implicit and explicit drag indistinguishable")

      end block checks
      call bd_i%destroy(); call ms_i%destroy()
      call bd_e%destroy(); call ms_e%destroy()
   end subroutine test_implicit_thin_layer

end module test_ocean_bottom_drag
