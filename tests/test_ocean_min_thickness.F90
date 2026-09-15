module test_ocean_min_thickness
   !! Analytical conservation tests for the conservative minimum-thickness
   !! borrow (`rdb_ocean_min_thickness`).  Proves, on a hand-built column with a
   !! layer driven below the floor:
   !!   * total thickness  Sigma h_k          conserved to round-off,
   !!   * every layer >= floor afterwards,
   !!   * every tracer mass Sigma h_k.Tr_k     conserved to round-off,
   !!   * per-face momentum  Sigma h_face.u    conserved to round-off,
   !! and that a column with all layers already >= floor is byte-unchanged
   !! (strict no-op), including the default (knob-off) path.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_min_thickness, only: ocean_apply_conservative_min_thickness, &
                                      min_thickness_target_column

   implicit none
   private

   public :: collect_ocean_min_thickness_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: FLOOR = 0.1_wp

contains

   subroutine collect_ocean_min_thickness_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("column_target_conserves_and_floors", &
                               test_column_target), &
                  new_unittest("field_conserves_h_mom_tracers", &
                               test_field_conservation), &
                  new_unittest("thick_column_is_byte_noop", &
                               test_thick_noop) &
                  ]
   end subroutine collect_ocean_min_thickness_tests

   ! ---- Test 1: per-column target builder (pure, no state) --------------------
   subroutine test_column_target(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h_old(NZ), h_new(NZ), h_flat(NZ), h_flat_new(NZ)
      real(wp) :: sum_old, sum_new
      logical :: grounded
      integer :: k

      ! A column with layer 2 well below the floor.
      h_old = [50.0_wp, 0.02_wp, 30.0_wp, 20.0_wp]
      call min_thickness_target_column(NZ, h_old, FLOOR, h_new, grounded)

      call check(error, grounded, "grounded column must report grounded=.true.")
      if (allocated(error)) return

      sum_old = sum(h_old)
      sum_new = sum(h_new)
      call check(error, abs(sum_new - sum_old) <= 1.0e-12_wp*sum_old, &
                 "column total thickness must be conserved")
      if (allocated(error)) return

      do k = 1, NZ
         call check(error, h_new(k) >= FLOOR - 1.0e-14_wp, &
                    "every layer must meet the floor after borrow")
         if (allocated(error)) return
      end do

      ! An already-thick column: strict no-op (grounded=.false., exact copy).
      h_flat = [25.0_wp, 25.0_wp, 25.0_wp, 25.0_wp]
      call min_thickness_target_column(NZ, h_flat, FLOOR, h_flat_new, grounded)
      call check(error,.not. grounded, "thick column must report grounded=.false.")
      if (allocated(error)) return
      do k = 1, NZ
         call check(error, h_flat_new(k) == h_flat(k), &
                    "thick column target must be a byte-exact copy")
         if (allocated(error)) return
      end do
   end subroutine test_column_target

   ! ---- Test 2: full-field conservation through the state routine -------------
   subroutine test_field_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp), allocatable :: h_new(:, :, :)
      real(wp), allocatable :: gmask(:, :, :)
      real(wp), allocatable :: h0(:, :, :), u0(:, :, :), v0(:, :, :)
      real(wp), allocatable :: hTrS0(:, :, :), hTrT0(:, :, :)
      integer :: nx, ny, i, j, k, ic, jc, ifar, jfar
      real(wp) :: s_old, s_new, mom0, mom1, err_x, err_y

      call grid%init(3, 3, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)

      allocate (h_new(nx, ny, NZ))
      allocate (gmask(nx, ny, 1))

      ! HOST setup (plain do, NOT do concurrent): ms is not device-mapped yet,
      ! and under -stdpar=gpu a do concurrent here launches on the device
      ! against unmapped host memory -> CUDA_ERROR_ILLEGAL_ADDRESS. Host-fill
      ! first; enter_data (copyin) below makes ms device-resident.
      ! Thick everywhere (incl. ghosts) so only our test column grounds.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 25.0_wp
            end do
         end do
      end do
      ! Distinct per-layer tracer concentrations and face velocities.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = (35.0_wp + real(k, wp))*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = (10.0_wp + 2.0_wp*real(k, wp))*ms%h_layer(i, j, k)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = 0.05_wp*real(k, wp) + 0.01_wp*real(i, wp)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = -0.03_wp*real(k, wp) + 0.02_wp*real(j, wp)
            end do
         end do
      end do

      ! Interior physical cell with layer 2 driven below floor; keep the same
      ! column total by pushing the removed volume into the other layers.
      ic = NGHOST + 2
      jc = NGHOST + 2
      ms%h_layer(ic, jc, 1) = 50.0_wp
      ms%h_layer(ic, jc, 2) = 0.02_wp
      ms%h_layer(ic, jc, 3) = 29.99_wp
      ms%h_layer(ic, jc, 4) = 19.99_wp
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(ic, jc, k) = (35.0_wp + real(k, wp))*ms%h_layer(ic, jc, k)
         ms%tracers(ms%idx_temperature)%hTr(ic, jc, k) = (10.0_wp + 2.0_wp*real(k, wp))*ms%h_layer(ic, jc, k)
      end do

      ! Snapshots (host, before the device map).
      h0 = ms%h_layer
      u0 = ms%u_face_x_layer
      v0 = ms%v_face_y_layer
      hTrS0 = ms%tracers(ms%idx_salinity)%hTr
      hTrT0 = ms%tracers(ms%idx_temperature)%hTr

      ! Make ms device-resident (copyin of h_layer/faces/tracers%hTr) + map the
      ! caller-owned scratch h_new; the kernel runs on-device. exit_data copies
      ! h_layer/faces/tracers%hTr back to the host for the checks below.
      call ms%enter_data()
      !$acc enter data create(h_new, gmask)

      call ocean_apply_conservative_min_thickness(grid, ms, h_new, gmask, FLOOR)

      !$acc exit data delete(h_new, gmask)
      call ms%exit_data()

      ! (a) Total thickness of the grounded column conserved.
      s_old = sum(h0(ic, jc, :))
      s_new = sum(ms%h_layer(ic, jc, :))
      call check(error, abs(s_new - s_old) <= 1.0e-12_wp*s_old, &
                 "grounded-column total thickness must be conserved")
      if (allocated(error)) return

      ! (b) Every layer meets the floor.
      do k = 1, NZ
         call check(error, ms%h_layer(ic, jc, k) >= FLOOR - 1.0e-13_wp, &
                    "every layer of the grounded column must meet the floor")
         if (allocated(error)) return
      end do

      ! (c) Tracer mass conserved (both prognostic tracers).
      s_old = sum(hTrS0(ic, jc, :))
      s_new = sum(ms%tracers(ms%idx_salinity)%hTr(ic, jc, :))
      call check(error, abs(s_new - s_old) <= 1.0e-12_wp*abs(s_old), &
                 "salinity mass must be conserved in the grounded column")
      if (allocated(error)) return
      s_old = sum(hTrT0(ic, jc, :))
      s_new = sum(ms%tracers(ms%idx_temperature)%hTr(ic, jc, :))
      call check(error, abs(s_new - s_old) <= 1.0e-12_wp*abs(s_old), &
                 "temperature mass must be conserved in the grounded column")
      if (allocated(error)) return

      ! (d) Per-face momentum Sigma_k h_face.u conserved on every face.
      !     before uses h0 face-avg & u0; after uses new h_layer & u.
      !     Loop bodies must stay branch-free (worst-defect reduction, ONE
      !     check per direction): nvfortran 25.9/26.5 ICEs ("flowgraph:
      !     node is zero") at -O2+ on the CPU path with check/early-return
      !     inside these nested face loops.
      err_x = 0.0_wp
      do j = 1, ny
         do i = 1, nx + 1
            mom0 = 0.0_wp
            mom1 = 0.0_wp
            do k = 1, NZ
               mom0 = mom0 + face_h(h0, nx, i, j, k)*u0(i, j, k)
               mom1 = mom1 + face_h(ms%h_layer, nx, i, j, k)*ms%u_face_x_layer(i, j, k)
            end do
            err_x = max(err_x, abs(mom1 - mom0)/max(abs(mom0), 1.0_wp))
         end do
      end do
      call check(error, err_x <= 1.0e-11_wp, &
                 "x-face momentum must be conserved")
      if (allocated(error)) return
      err_y = 0.0_wp
      do j = 1, ny + 1
         do i = 1, nx
            mom0 = 0.0_wp
            mom1 = 0.0_wp
            do k = 1, NZ
               mom0 = mom0 + face_hy(h0, ny, i, j, k)*v0(i, j, k)
               mom1 = mom1 + face_hy(ms%h_layer, ny, i, j, k)*ms%v_face_y_layer(i, j, k)
            end do
            err_y = max(err_y, abs(mom1 - mom0)/max(abs(mom0), 1.0_wp))
         end do
      end do
      call check(error, err_y <= 1.0e-11_wp, &
                 "y-face momentum must be conserved")
      if (allocated(error)) return

      ! (e) A far ungrounded cell is byte-unchanged (no isopycnal pinning).
      ifar = NGHOST + 1
      jfar = NGHOST + 1
      do k = 1, NZ
         call check(error, ms%h_layer(ifar, jfar, k) == h0(ifar, jfar, k), &
                    "ungrounded cell h_layer must be byte-unchanged")
         if (allocated(error)) return
         call check(error, ms%tracers(ms%idx_salinity)%hTr(ifar, jfar, k) == hTrS0(ifar, jfar, k), &
                    "ungrounded cell salinity must be byte-unchanged")
         if (allocated(error)) return
      end do

      call ms%destroy()
   end subroutine test_field_conservation

   ! ---- Test 3: strict no-op when every layer already meets the floor ---------
   subroutine test_thick_noop(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp), allocatable :: h_new(:, :, :)
      real(wp), allocatable :: gmask(:, :, :)
      real(wp), allocatable :: h0(:, :, :), u0(:, :, :), v0(:, :, :)
      real(wp), allocatable :: hTrS0(:, :, :), hTrT0(:, :, :)
      integer :: nx, ny, i, j, k

      call grid%init(3, 3, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      allocate (h_new(nx, ny, NZ))
      allocate (gmask(nx, ny, 1))

      ! HOST setup (plain do; enter_data maps ms to the device below).
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 10.0_wp + real(k, wp)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 34.0_wp*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 12.0_wp*ms%h_layer(i, j, k)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = 0.2_wp*real(k, wp)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = 0.3_wp*real(k, wp)
            end do
         end do
      end do

      h0 = ms%h_layer
      u0 = ms%u_face_x_layer
      v0 = ms%v_face_y_layer
      hTrS0 = ms%tracers(ms%idx_salinity)%hTr
      hTrT0 = ms%tracers(ms%idx_temperature)%hTr

      call ms%enter_data()
      !$acc enter data create(h_new, gmask)

      call ocean_apply_conservative_min_thickness(grid, ms, h_new, gmask, FLOOR)

      !$acc exit data delete(h_new, gmask)
      call ms%exit_data()

      call check(error, all(ms%h_layer == h0), "h_layer must be byte-identical (no-op)")
      if (allocated(error)) return
      call check(error, all(ms%u_face_x_layer == u0), "u faces must be byte-identical (no-op)")
      if (allocated(error)) return
      call check(error, all(ms%v_face_y_layer == v0), "v faces must be byte-identical (no-op)")
      if (allocated(error)) return
      call check(error, all(ms%tracers(ms%idx_salinity)%hTr == hTrS0), &
                 "salinity must be byte-identical (no-op)")
      if (allocated(error)) return
      call check(error, all(ms%tracers(ms%idx_temperature)%hTr == hTrT0), &
                 "temperature must be byte-identical (no-op)")
      if (allocated(error)) return

      call ms%destroy()
   end subroutine test_thick_noop

   ! ---- Face-thickness helpers (arithmetic mean, wall = single cell) ----------
   pure function face_h(h, nx, i, j, k) result(hf)
      !! East-face thickness at x-face i (1..nx+1) using the same
      !! arithmetic-mean / wall convention as the remap.
      real(wp), intent(in) :: h(:, :, :)
      integer, intent(in) :: nx, i, j, k
      real(wp) :: hf
      if (i == 1) then
         hf = h(1, j, k)
      else if (i == nx + 1) then
         hf = h(nx, j, k)
      else
         hf = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
      end if
   end function face_h

   pure function face_hy(h, ny, i, j, k) result(hf)
      !! North-face thickness at y-face j (1..ny+1).
      real(wp), intent(in) :: h(:, :, :)
      integer, intent(in) :: ny, i, j, k
      real(wp) :: hf
      if (j == 1) then
         hf = h(i, 1, k)
      else if (j == ny + 1) then
         hf = h(i, ny, k)
      else
         hf = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
      end if
   end function face_hy

end module test_ocean_min_thickness
