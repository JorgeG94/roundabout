!! Unit tests for the C-grid barotropic + multilayer state types
!! introduced in Phase 1 of the ocean dynamical core (see
!! `docs/ROADMAP_OCEAN.md`).  Checks:
!!
!!   * Allocation shapes (centre + east-face + north-face stagger)
!!   * Zero-fill on init
!!   * Symmetric destroy
!!   * GPU device-mapping round trip: host → device → mutate on
!!     device → device → host, asserting the host sees the
!!     device-mutated values.
!!
!! The round-trip cases construct the state types standalone (not
!! through `state_t`) so they exercise just the bound `enter_data` /
!! `exit_data` plumbing without dragging in the rest of the solver.
!! Mirrors the device-cycle pattern in `test_state` and the
!! `!$acc parallel loop present(state)` pattern in
!! `test_ml_dynamics_unstr`.
module test_state_roundtrip
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   implicit none
   private

   public :: collect_state_roundtrip_tests

   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_state_roundtrip_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("barotropic_init_shapes", test_barotropic_shapes), &
                  new_unittest("barotropic_init_zeroes", test_barotropic_zeroes), &
                  new_unittest("barotropic_destroy", test_barotropic_destroy), &
                  new_unittest("barotropic_gpu_roundtrip", test_barotropic_roundtrip), &
                  new_unittest("multilayer_init_shapes", test_multilayer_shapes), &
                  new_unittest("multilayer_init_zeroes", test_multilayer_zeroes), &
                  new_unittest("multilayer_tracers_registered", test_multilayer_tracers), &
                  new_unittest("multilayer_destroy", test_multilayer_destroy), &
                  new_unittest("multilayer_gpu_roundtrip", test_multilayer_roundtrip) &
                  ]
   end subroutine collect_state_roundtrip_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   ! -----------------------------------------------------------------
   ! Barotropic C-grid state
   ! -----------------------------------------------------------------

   subroutine test_barotropic_shapes(error)
      !! Centre arrays are (nx_total, ny_total); east-face arrays carry
      !! one extra i-column; north-face arrays carry one extra j-row.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      integer :: nx, ny

      call make_grid(grid)
      call bs%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call check(error, bs%is_init, "is_init should be true after init")
      if (allocated(error)) return

      call check(error, size(bs%h, 1) == nx .and. size(bs%h, 2) == ny, &
                 "h should be (nx_total, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%b, 1) == nx .and. size(bs%b, 2) == ny, &
                 "b should be (nx_total, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%u_face_x, 1) == nx + 1 .and. size(bs%u_face_x, 2) == ny, &
                 "u_face_x should be (nx_total+1, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%hu_face_x, 1) == nx + 1 .and. size(bs%hu_face_x, 2) == ny, &
                 "hu_face_x should be (nx_total+1, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%mass_flux_x, 1) == nx + 1 .and. size(bs%mass_flux_x, 2) == ny, &
                 "mass_flux_x should be (nx_total+1, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%v_face_y, 1) == nx .and. size(bs%v_face_y, 2) == ny + 1, &
                 "v_face_y should be (nx_total, ny_total+1)")
      if (allocated(error)) return
      call check(error, size(bs%hv_face_y, 1) == nx .and. size(bs%hv_face_y, 2) == ny + 1, &
                 "hv_face_y should be (nx_total, ny_total+1)")
      if (allocated(error)) return
      call check(error, size(bs%mass_flux_y, 1) == nx .and. size(bs%mass_flux_y, 2) == ny + 1, &
                 "mass_flux_y should be (nx_total, ny_total+1)")
      if (allocated(error)) return
      call check(error, size(bs%u_face_x0, 1) == nx + 1 .and. size(bs%u_face_x0, 2) == ny, &
                 "u_face_x0 should be (nx_total+1, ny_total)")
      if (allocated(error)) return
      call check(error, size(bs%v_face_y0, 1) == nx .and. size(bs%v_face_y0, 2) == ny + 1, &
                 "v_face_y0 should be (nx_total, ny_total+1)")

      call bs%destroy()
   end subroutine test_barotropic_shapes

   subroutine test_barotropic_zeroes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs

      call make_grid(grid)
      call bs%init(grid)

      call check(error, maxval(abs(bs%h)) < 1.0e-15_wp, "h zero")
      if (allocated(error)) return
      call check(error, maxval(abs(bs%u_face_x)) < 1.0e-15_wp, "u_face_x zero")
      if (allocated(error)) return
      call check(error, maxval(abs(bs%v_face_y)) < 1.0e-15_wp, "v_face_y zero")
      if (allocated(error)) return
      call check(error, maxval(abs(bs%mass_flux_x)) < 1.0e-15_wp, "mass_flux_x zero")
      if (allocated(error)) return
      call check(error, maxval(abs(bs%mass_flux_y)) < 1.0e-15_wp, "mass_flux_y zero")

      call bs%destroy()
   end subroutine test_barotropic_zeroes

   subroutine test_barotropic_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs

      call make_grid(grid)
      call bs%init(grid)
      call bs%destroy()

      call check(error,.not. bs%is_init, "is_init false after destroy")
      if (allocated(error)) return
      call check(error,.not. allocated(bs%h), "h deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(bs%u_face_x), "u_face_x deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(bs%v_face_y), "v_face_y deallocated")
   end subroutine test_barotropic_destroy

   subroutine test_barotropic_roundtrip(error)
      !! init → enter_data → mutate on device → exit_data → host check.
      !! Tagging each cell-centre, east-face, and north-face slot with
      !! a distinct value verifies the stagger ships across the
      !! transfer (a swapped or off-by-one face stride would land
      !! values in the wrong cell).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      integer :: i, j, nx, ny
      real(wp) :: expected

      call make_grid(grid)
      call bs%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! Map parent + components.  Parent first per the
      ! derived-type-GPU-mapping rule (see CLAUDE.md memory).
      !$acc enter data copyin(bs)
      call bs%enter_data()

      ! Mutate every cell-centred slot, every east face, every north
      ! face — three loops with distinct shapes so a stride bug in
      ! any of them stands out.
      !$acc parallel loop collapse(2) present(bs)
      do j = 1, ny
         do i = 1, nx
            bs%h(i, j) = real(100*i + j, wp)
         end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(2) present(bs)
      do j = 1, ny
         do i = 1, nx + 1
            bs%u_face_x(i, j) = real(1000*i + j, wp)
         end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(2) present(bs)
      do j = 1, ny + 1
         do i = 1, nx
            bs%v_face_y(i, j) = real(10000*i + j, wp)
         end do
      end do
      !$acc end parallel loop

      call bs%exit_data()
      !$acc exit data delete(bs)

      ! Host should now see the device-set tags.
      do j = 1, ny
         do i = 1, nx
            expected = real(100*i + j, wp)
            call check(error, abs(bs%h(i, j) - expected) < 1.0e-12_wp, &
                       "h round-trip mismatch")
            if (allocated(error)) then
               call bs%destroy()
               return
            end if
         end do
      end do

      do j = 1, ny
         do i = 1, nx + 1
            expected = real(1000*i + j, wp)
            call check(error, abs(bs%u_face_x(i, j) - expected) < 1.0e-12_wp, &
                       "u_face_x round-trip mismatch (east-face stagger)")
            if (allocated(error)) then
               call bs%destroy()
               return
            end if
         end do
      end do

      do j = 1, ny + 1
         do i = 1, nx
            expected = real(10000*i + j, wp)
            call check(error, abs(bs%v_face_y(i, j) - expected) < 1.0e-12_wp, &
                       "v_face_y round-trip mismatch (north-face stagger)")
            if (allocated(error)) then
               call bs%destroy()
               return
            end if
         end do
      end do

      call bs%destroy()
   end subroutine test_barotropic_roundtrip

   ! -----------------------------------------------------------------
   ! Multilayer C-grid state
   ! -----------------------------------------------------------------

   subroutine test_multilayer_shapes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: nx, ny

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call check(error, ms%is_init, "is_init should be true after init")
      if (allocated(error)) return

      call check(error, all(shape(ms%h_layer) == [nx, ny, NZ]), &
                 "h_layer should be (nx_total, ny_total, nz_ml)")
      if (allocated(error)) return
      call check(error, all(shape(ms%u_face_x_layer) == [nx + 1, ny, NZ]), &
                 "u_face_x_layer should be (nx_total+1, ny_total, nz_ml)")
      if (allocated(error)) return
      call check(error, all(shape(ms%v_face_y_layer) == [nx, ny + 1, NZ]), &
                 "v_face_y_layer should be (nx_total, ny_total+1, nz_ml)")
      if (allocated(error)) return
      call check(error, all(shape(ms%mass_flux_x_layer) == [nx + 1, ny, NZ]), &
                 "mass_flux_x_layer should be (nx_total+1, ny_total, nz_ml)")
      if (allocated(error)) return
      call check(error, all(shape(ms%mass_flux_y_layer) == [nx, ny + 1, NZ]), &
                 "mass_flux_y_layer should be (nx_total, ny_total+1, nz_ml)")
      if (allocated(error)) return
      call check(error, all(shape(ms%w_interface) == [nx, ny, NZ + 1]), &
                 "w_interface should be (nx_total, ny_total, nz_ml+1)")

      call ms%destroy()
   end subroutine test_multilayer_shapes

   subroutine test_multilayer_zeroes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)

      call check(error, maxval(abs(ms%h_layer)) < 1.0e-15_wp, "h_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%u_face_x_layer)) < 1.0e-15_wp, "u_face_x_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%v_face_y_layer)) < 1.0e-15_wp, "v_face_y_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%w_interface)) < 1.0e-15_wp, "w_interface zero")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr)) < 1.0e-15_wp, &
                 "salinity hTr zero")
      if (allocated(error)) return
      call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr)) < 1.0e-15_wp, &
                 "temperature hTr zero")

      call ms%destroy()
   end subroutine test_multilayer_zeroes

   subroutine test_multilayer_tracers(error)
      !! Salinity + temperature registered with the expected metadata.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)

      call check(error, allocated(ms%tracers), "tracers allocated")
      if (allocated(error)) return
      call check(error, size(ms%tracers) == 2, "two tracers registered (S + T)")
      if (allocated(error)) return
      call check(error, ms%idx_salinity == 1, "salinity at index 1")
      if (allocated(error)) return
      call check(error, ms%idx_temperature == 2, "temperature at index 2")
      if (allocated(error)) return
      call check(error, trim(ms%tracers(ms%idx_salinity)%name) == "salinity", &
                 "salinity tracer named")
      if (allocated(error)) return
      call check(error, trim(ms%tracers(ms%idx_temperature)%name) == "temperature", &
                 "temperature tracer named")

      call ms%destroy()
   end subroutine test_multilayer_tracers

   subroutine test_multilayer_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ms%destroy()

      call check(error,.not. ms%is_init, "is_init false after destroy")
      if (allocated(error)) return
      call check(error,.not. allocated(ms%h_layer), "h_layer deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(ms%u_face_x_layer), "u_face_x_layer deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(ms%tracers), "tracers deallocated")
      if (allocated(error)) return
      call check(error, ms%idx_salinity == 0, "idx_salinity reset")
      if (allocated(error)) return
      call check(error, ms%idx_temperature == 0, "idx_temperature reset")
   end subroutine test_multilayer_destroy

   subroutine test_multilayer_roundtrip(error)
      !! 3D round trip: prognostic centre + east-face + north-face +
      !! interface arrays + a tracer hTr.  Each tagged with a
      !! distinct (i, j, k)-encoded value so a stride/stagger bug at
      !! any axis stands out.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: i, j, k, nx, ny, idx_s
      real(wp) :: expected

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      idx_s = ms%idx_salinity

      !$acc enter data copyin(ms)
      call ms%enter_data()

      ! h_layer at cell centres, shape (nx, ny, nz)
      !$acc parallel loop collapse(3) present(ms)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = real(10000*k + 100*i + j, wp)
            end do
         end do
      end do
      !$acc end parallel loop

      ! u_face_x_layer at east faces, shape (nx+1, ny, nz)
      !$acc parallel loop collapse(3) present(ms)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = real(20000*k + 100*i + j, wp)
            end do
         end do
      end do
      !$acc end parallel loop

      ! v_face_y_layer at north faces, shape (nx, ny+1, nz)
      !$acc parallel loop collapse(3) present(ms)
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = real(30000*k + 100*i + j, wp)
            end do
         end do
      end do
      !$acc end parallel loop

      ! w_interface at layer interfaces, shape (nx, ny, nz+1)
      !$acc parallel loop collapse(3) present(ms)
      do k = 1, NZ + 1
         do j = 1, ny
            do i = 1, nx
               ms%w_interface(i, j, k) = real(40000*k + 100*i + j, wp)
            end do
         end do
      end do
      !$acc end parallel loop

      ! Salinity hTr at cell centres, shape (nx, ny, nz).  Exercises
      ! the two-step tracer-registry attach (parent array descriptor
      ! + per-element hTr) at the deepest level the ocean ml state
      ! reaches.
      !$acc parallel loop collapse(3) present(ms)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(idx_s)%hTr(i, j, k) = real(50000*k + 100*i + j, wp)
            end do
         end do
      end do
      !$acc end parallel loop

      call ms%exit_data()
      !$acc exit data delete(ms)

      ! Verify
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               expected = real(10000*k + 100*i + j, wp)
               call check(error, abs(ms%h_layer(i, j, k) - expected) < 1.0e-12_wp, &
                          "h_layer round-trip mismatch")
               if (allocated(error)) then
                  call ms%destroy(); return
               end if
            end do
         end do
      end do

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               expected = real(20000*k + 100*i + j, wp)
               call check(error, abs(ms%u_face_x_layer(i, j, k) - expected) < 1.0e-12_wp, &
                          "u_face_x_layer round-trip mismatch")
               if (allocated(error)) then
                  call ms%destroy(); return
               end if
            end do
         end do
      end do

      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               expected = real(30000*k + 100*i + j, wp)
               call check(error, abs(ms%v_face_y_layer(i, j, k) - expected) < 1.0e-12_wp, &
                          "v_face_y_layer round-trip mismatch")
               if (allocated(error)) then
                  call ms%destroy(); return
               end if
            end do
         end do
      end do

      do k = 1, NZ + 1
         do j = 1, ny
            do i = 1, nx
               expected = real(40000*k + 100*i + j, wp)
               call check(error, abs(ms%w_interface(i, j, k) - expected) < 1.0e-12_wp, &
                          "w_interface round-trip mismatch")
               if (allocated(error)) then
                  call ms%destroy(); return
               end if
            end do
         end do
      end do

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               expected = real(50000*k + 100*i + j, wp)
               call check(error, &
                          abs(ms%tracers(idx_s)%hTr(i, j, k) - expected) < 1.0e-12_wp, &
                          "salinity hTr round-trip mismatch")
               if (allocated(error)) then
                  call ms%destroy(); return
               end if
            end do
         end do
      end do

      call ms%destroy()
   end subroutine test_multilayer_roundtrip

end module test_state_roundtrip
