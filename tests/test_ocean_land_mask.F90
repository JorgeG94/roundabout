!! Unit tests for the ocean static land-mask foundation (CHUNK A) in
!! `rdb_ocean_metrics` (`metrics_apply_land_mask`).
!!
!! Coverage (spec §14 C3 + §9.1):
!!   * face/corner mask derivation: `wet_u = wet_T(i-1)*wet_T(i)`,
!!     `wet_v = wet_T(j-1)*wet_T(j)`, `wet_q = product of the 4 corner
!!     T-cells`, with an interior island.
!!   * `metrics_apply_land_mask` zeros EXACTLY the 6 face metrics
!!     (`dy_cu, idxCu, dxCu` at land u-faces; `dx_cv, idyCv, dyCv` at
!!     land v-faces) and leaves the DO-NOT-ZERO set + every wet face
!!     untouched.
!!   * all-wet bit-identity: `wet_mask≡1` ⇒ masks≡1, metrics byte-unchanged.
module test_ocean_land_mask
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_finalize, &
                                metrics_fill_cartesian, metrics_apply_land_mask
   use rdb_continuity, only: ppm_mirror_h
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_dyn, only: mask_layer_velocities
   implicit none
   private

   public :: collect_ocean_land_mask_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 250.0_wp, DY = 400.0_wp

contains

   subroutine collect_ocean_land_mask_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("land_mask_face_corner_derivation", test_face_corner), &
                  new_unittest("land_mask_zeros_exactly_six_metrics", test_six_metrics), &
                  new_unittest("land_mask_donotzero_untouched", test_donotzero), &
                  new_unittest("land_mask_all_wet_bit_identical", test_all_wet), &
                  new_unittest("land_mask_wet_T_stored_for_ppm_mirror", test_wet_T_stored), &
                  new_unittest("ppm_mirror_h_reflects_land_neighbour", test_ppm_mirror), &
                  new_unittest("wall_velocity_mask_zeros_wall_faces", test_wall_velocity_mask), &
                  new_unittest("wall_velocity_mask_off_bit_identical", test_wall_velocity_off) &
                  ]
   end subroutine collect_ocean_land_mask_tests

   subroutine test_ppm_mirror(error)
      !! The mirror-h helper (spec §14 C2): a wet neighbour (w=1) is a
      !! literal no-op (returns the neighbour, byte-exact ⇒ all-wet
      !! bit-identity), a land neighbour (w=0) is mirrored to the local
      !! cell's value so the PPM parabola sees a flat coast.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HNBR = 1234.5_wp, HLOC = 67.0_wp
      ! Wet neighbour: byte-identical passthrough (the bit-identity hinge).
      call check(error, ppm_mirror_h(HNBR, HLOC, 1.0_wp) == HNBR, &
                 "mirror w=1 must return the neighbour byte-exactly")
      if (allocated(error)) return
      ! Land neighbour: reflected to the local cell ⇒ zero stencil gradient.
      call check(error, ppm_mirror_h(HNBR, HLOC, 0.0_wp) == HLOC, &
                 "mirror w=0 must return the local cell (reflected coast)")
   end subroutine test_ppm_mirror

   function make_grid(nx_phys, ny_phys) result(g)
      integer, intent(in) :: nx_phys, ny_phys
      type(hgrid_t) :: g
      call g%init(nx_phys, ny_phys, NGHOST, DX, DY)
   end function make_grid

   !! Build a wet_mask (nx_total, ny_total) all-wet except a 2x2 interior
   !! island.  The island occupies T-cells (ci, cj)..(ci+1, cj+1).
   subroutine build_island_mask(g, ci, cj, wm)
      type(hgrid_t), intent(in) :: g
      integer, intent(in) :: ci, cj
      real(wp), allocatable, intent(out) :: wm(:, :)
      allocate (wm(g%nx_total, g%ny_total), source=1.0_wp)
      wm(ci, cj) = 0.0_wp
      wm(ci + 1, cj) = 0.0_wp
      wm(ci, cj + 1) = 0.0_wp
      wm(ci + 1, cj + 1) = 0.0_wp
   end subroutine build_island_mask

   ! ================================================================
   ! T1: face/corner mask derivation vs the analytic product rule
   ! ================================================================
   subroutine test_face_corner(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)
      integer :: i, j, ci, cj
      real(wp) :: expect_u, expect_v, expect_q

      g = make_grid(6, 6)
      ci = 4
      cj = 4   ! interior 2x2 island, clear of the ghost ring
      call build_island_mask(g, ci, cj, wm)

      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)

      ! wet_u(i,j) = wet_T(i-1,j)*wet_T(i,j), i in [2,nx]; outer walls 0.
      do j = 1, g%ny_total
         call check(error, m%wet_u(1, j) == 0.0_wp, "wet_u west wall not 0")
         if (allocated(error)) return
         call check(error, m%wet_u(g%nx_total + 1, j) == 0.0_wp, "wet_u east wall not 0")
         if (allocated(error)) return
         do i = 2, g%nx_total
            expect_u = wm(i - 1, j)*wm(i, j)
            call check(error, m%wet_u(i, j) == expect_u, "wet_u /= wet_T(i-1)*wet_T(i)")
            if (allocated(error)) return
         end do
      end do

      ! wet_v(i,j) = wet_T(i,j-1)*wet_T(i,j), j in [2,ny]; outer walls 0.
      do i = 1, g%nx_total
         call check(error, m%wet_v(i, 1) == 0.0_wp, "wet_v south wall not 0")
         if (allocated(error)) return
         call check(error, m%wet_v(i, g%ny_total + 1) == 0.0_wp, "wet_v north wall not 0")
         if (allocated(error)) return
         do j = 2, g%ny_total
            expect_v = wm(i, j - 1)*wm(i, j)
            call check(error, m%wet_v(i, j) == expect_v, "wet_v /= wet_T(j-1)*wet_T(j)")
            if (allocated(error)) return
         end do
      end do

      ! wet_q(i,j) = product of the 4 surrounding T-cells (interior only).
      do j = 2, g%ny_total
         do i = 2, g%nx_total
            expect_q = wm(i - 1, j - 1)*wm(i, j - 1)*wm(i - 1, j)*wm(i, j)
            call check(error, m%wet_q(i, j) == expect_q, "wet_q /= 4-cell product")
            if (allocated(error)) return
         end do
      end do

      ! Spot check: the 4 corners touching the island are land (0).
      ! Corner (ci, cj) is SW of island cell (ci,cj): T-cells (ci-1,cj-1),
      ! (ci,cj-1),(ci-1,cj),(ci,cj) — (ci,cj) is land ⇒ wet_q=0.
      call check(error, m%wet_q(ci, cj) == 0.0_wp, "SW island corner not land")
      if (allocated(error)) return
      call check(error, m%wet_q(ci + 2, cj + 2) == 0.0_wp, "NE island corner not land")
      if (allocated(error)) return
      ! The corner one cell SW of the island touches only wet T-cells ⇒ wet.
      call check(error, m%wet_q(ci - 1, cj - 1) == 1.0_wp, "far SW corner not wet")
   end subroutine test_face_corner

   ! ================================================================
   ! T2: exactly-6 metrics zeroed at land faces, scaled by mask elsewhere
   ! ================================================================
   subroutine test_six_metrics(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)
      integer :: i, j, ci, cj

      g = make_grid(6, 6)
      ci = 4
      cj = 4
      call build_island_mask(g, ci, cj, wm)

      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)

      ! u-face metrics: dy_cu, idxCu, dxCu == base * wet_u everywhere.
      do j = 1, g%ny_total
         do i = 1, g%nx_total + 1
            call check(error, m%dy_cu(i, j) == DY*m%wet_u(i, j), "dy_cu not masked")
            if (allocated(error)) return
            call check(error, m%idxCu(i, j) == (1.0_wp/DX)*m%wet_u(i, j), "idxCu not masked")
            if (allocated(error)) return
            call check(error, m%dxCu(i, j) == DX*m%wet_u(i, j), "dxCu not masked")
            if (allocated(error)) return
         end do
      end do

      ! v-face metrics: dx_cv, idyCv, dyCv == base * wet_v everywhere.
      do j = 1, g%ny_total + 1
         do i = 1, g%nx_total
            call check(error, m%dx_cv(i, j) == DX*m%wet_v(i, j), "dx_cv not masked")
            if (allocated(error)) return
            call check(error, m%idyCv(i, j) == (1.0_wp/DY)*m%wet_v(i, j), "idyCv not masked")
            if (allocated(error)) return
            call check(error, m%dyCv(i, j) == DY*m%wet_v(i, j), "dyCv not masked")
            if (allocated(error)) return
         end do
      end do

      ! At least one land u-face really hit zero (not a vacuous test):
      ! the west face of island cell (ci,cj) = u-face ci is land.
      call check(error, m%dy_cu(ci, cj) == 0.0_wp, "expected land u-face dy_cu=0")
      if (allocated(error)) return
      call check(error, m%idxCu(ci, cj) == 0.0_wp, "expected land u-face idxCu=0")
      if (allocated(error)) return
      ! And a known-wet u-face is unchanged.
      call check(error, m%dy_cu(2, 2) == DY, "wet u-face dy_cu altered")
   end subroutine test_six_metrics

   ! ================================================================
   ! T3: the DO-NOT-ZERO metrics are left byte-identical to the
   !     pre-mask (uniform-Cartesian) values, even at land faces/cells.
   ! ================================================================
   subroutine test_donotzero(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)
      real(wp) :: area, iarea_t, iarea_bu

      g = make_grid(6, 6)
      call build_island_mask(g, 4, 4, wm)

      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)

      area = DX*DY
      iarea_t = 1.0_wp/area
      iarea_bu = 1.0_wp/area
      ! iareaT / areaT / areaCu / areaCv / iareaBu untouched (every cell).
      call check(error, all(m%areaT == area), "areaT altered")
      if (allocated(error)) return
      call check(error, all(m%areaCu == area), "areaCu altered")
      if (allocated(error)) return
      call check(error, all(m%areaCv == area), "areaCv altered")
      if (allocated(error)) return
      call check(error, all(m%iareaT == iarea_t), "iareaT altered")
      if (allocated(error)) return
      call check(error, all(m%iareaBu == iarea_bu), "iareaBu altered")
   end subroutine test_donotzero

   ! ================================================================
   ! T4: all-wet domain ⇒ masks≡1, the 6 metrics byte-unchanged.
   ! ================================================================
   subroutine test_all_wet(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)

      g = make_grid(6, 6)
      allocate (wm(g%nx_total, g%ny_total), source=1.0_wp)

      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)

      ! Interior face masks all 1 (only the outer walls are forced to 0).
      call check(error, all(m%wet_u(2:g%nx_total, :) == 1.0_wp), "interior wet_u /= 1")
      if (allocated(error)) return
      call check(error, all(m%wet_v(:, 2:g%ny_total) == 1.0_wp), "interior wet_v /= 1")
      if (allocated(error)) return
      ! Interior corners (need all 4 T-cells, i,j in [2,n]) all 1.
      call check(error, all(m%wet_q(2:g%nx_total, 2:g%ny_total) == 1.0_wp), &
                 "interior wet_q /= 1")
      if (allocated(error)) return

      ! The 6 face metrics: identical to the unmasked cartesian fill at
      ! every interior face (where wet=1), byte-for-byte.
      call check(error, all(m%dy_cu(2:g%nx_total, :) == DY), "dy_cu changed all-wet")
      if (allocated(error)) return
      call check(error, all(m%dxCu(2:g%nx_total, :) == DX), "dxCu changed all-wet")
      if (allocated(error)) return
      call check(error, all(m%idxCu(2:g%nx_total, :) == 1.0_wp/DX), "idxCu changed all-wet")
      if (allocated(error)) return
      call check(error, all(m%dx_cv(:, 2:g%ny_total) == DX), "dx_cv changed all-wet")
      if (allocated(error)) return
      call check(error, all(m%dyCv(:, 2:g%ny_total) == DY), "dyCv changed all-wet")
      if (allocated(error)) return
      call check(error, all(m%idyCv(:, 2:g%ny_total) == 1.0_wp/DY), "idyCv changed all-wet")
   end subroutine test_all_wet

   !! CHUNK B (C2): the halo-valid T-cell mask `wet_T` must be stored on
   !! `ocean_metrics_t` so the continuity + tracer PPM reconstruction can
   !! mirror a land neighbour's thickness.  Checks: (a) all-wet ⇒ wet_T≡1
   !! (mirror never triggers, bit-identity), (b) island cells are 0 and
   !! their wet neighbours are 1.
   subroutine test_wet_T_stored(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)
      integer :: ci, cj

      ! ---- (a) all-wet ⇒ wet_T ≡ 1 ----
      g = make_grid(6, 6)
      allocate (wm(g%nx_total, g%ny_total), source=1.0_wp)
      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)
      call check(error, all(m%wet_T == 1.0_wp), "wet_T not all 1 for all-wet mask")
      if (allocated(error)) return
      deallocate (wm)

      ! ---- (b) island block: land cells 0, wet cells 1 ----
      ci = 4
      cj = 4
      call build_island_mask(g, ci, cj, wm)
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.false., periodic_y=.false., &
                                   north_fold=.false.)
      call check(error, m%wet_T(ci, cj) == 0.0_wp, "wet_T island SW cell not land")
      if (allocated(error)) return
      call check(error, m%wet_T(ci + 1, cj + 1) == 0.0_wp, "wet_T island NE cell not land")
      if (allocated(error)) return
      call check(error, m%wet_T(ci - 1, cj) == 1.0_wp, "wet_T west-of-island cell not wet")
      if (allocated(error)) return
      ! wet_T must equal the (interior) input mask exactly.
      call check(error, all(m%wet_T == wm), "wet_T /= halo-valid input mask")
   end subroutine test_wet_T_stored

   ! ================================================================
   ! Solid-wall velocity masking (mask_wall_velocity opt-in).
   ! Flat all-wet channel: periodic-x + WALL north/south.  The knob
   ! zeros the ghost wet_T beyond the N/S walls so the wall v-face masks
   ! vanish and mask_layer_velocities clears the wall-normal velocity —
   ! MOM6's mask-in-the-update, no bespoke velocity BC.
   ! ================================================================
   subroutine test_wall_velocity_mask(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      type(multilayer_state_t) :: ms
      real(wp), allocatable :: wm(:, :)
      integer :: ng, ni, nj, i, j, k, js_wall, jn_wall

      g = make_grid(6, 6)
      ng = g%nghost
      ni = g%nx_phys
      nj = g%ny_phys
      js_wall = ng + 1        ! v-face on the south wall (south of 1st phys cell)
      jn_wall = ng + nj + 1   ! v-face on the north wall (north of last phys cell)

      ! Flat all-wet channel (wet_mask ≡ 1 everywhere, incl. wall ghosts).
      allocate (wm(g%nx_total, g%ny_total), source=1.0_wp)
      call m%init(g)
      call metrics_fill_cartesian(m, g, DX, DY)
      call metrics_finalize(m)

      ! Knob ON: periodic-x, WALL n/s.
      call metrics_apply_land_mask(m, wm, g, &
                                   periodic_x=.true., periodic_y=.false., &
                                   north_fold=.false., &
                                   mask_wall_velocity=.true., &
                                   wall_west=.false., wall_east=.false., &
                                   wall_south=.true., wall_north=.true.)

      ! (a) wet_v == 0 at the N/S wall v-faces (every physical column).
      do i = ng + 1, ng + ni
         call check(error, m%wet_v(i, js_wall) == 0.0_wp, &
                    "wet_v not zeroed at south wall v-face")
         if (allocated(error)) return
         call check(error, m%wet_v(i, jn_wall) == 0.0_wp, &
                    "wet_v not zeroed at north wall v-face")
         if (allocated(error)) return
      end do
      ! Interior v-faces (between two wet physical cells) stay open.
      do i = ng + 1, ng + ni
         do j = js_wall + 1, jn_wall - 1
            call check(error, m%wet_v(i, j) == 1.0_wp, &
                       "interior v-face wrongly masked")
            if (allocated(error)) return
         end do
      end do
      ! wet_u stays 1 at the periodic-x seam u-faces (physical rows) — the
      ! periodic wrap keeps the west/east ghosts wet.
      do j = ng + 1, ng + nj
         do i = ng + 1, ng + ni + 1
            call check(error, m%wet_u(i, j) == 1.0_wp, &
                       "periodic-x u-face wrongly masked")
            if (allocated(error)) return
         end do
      end do

      ! (b) seed a non-zero v everywhere, run mask_layer_velocities, and
      ! confirm the wall-normal v is exactly 0 while the interior survives.
      ms%nz_ml = 3
      call ms%init(g)
      ms%v_face_y_layer = 1.7_wp
      ms%u_face_x_layer = 0.9_wp
      call mask_layer_velocities(g, m, ms)
      do k = 1, ms%nz_ml
         do i = ng + 1, ng + ni
            call check(error, ms%v_face_y_layer(i, js_wall, k) == 0.0_wp, &
                       "south wall-normal v not zeroed by mask_layer_velocities")
            if (allocated(error)) return
            call check(error, ms%v_face_y_layer(i, jn_wall, k) == 0.0_wp, &
                       "north wall-normal v not zeroed by mask_layer_velocities")
            if (allocated(error)) return
            call check(error, ms%v_face_y_layer(i, js_wall + 1, k) == 1.7_wp, &
                       "interior v wrongly zeroed by mask_layer_velocities")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_wall_velocity_mask

   subroutine test_wall_velocity_off(error)
      !! (c) With the knob OFF (and the legacy absent-arg call) the masks
      !! are byte-identical: the flat all-wet channel keeps every wet_* ≡ 1
      !! regardless of the WALL edges (nothing zeroed) ⇒ bit-identity.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m_off, m_legacy
      type(hgrid_t) :: g
      real(wp), allocatable :: wm(:, :)

      g = make_grid(6, 6)
      allocate (wm(g%nx_total, g%ny_total), source=1.0_wp)

      ! Legacy path (no wall-mask args at all).
      call m_legacy%init(g)
      call metrics_fill_cartesian(m_legacy, g, DX, DY)
      call metrics_finalize(m_legacy)
      call metrics_apply_land_mask(m_legacy, wm, g, &
                                   periodic_x=.true., periodic_y=.false., &
                                   north_fold=.false.)

      ! Knob explicitly OFF but WALL n/s tagged — must match the legacy path.
      call m_off%init(g)
      call metrics_fill_cartesian(m_off, g, DX, DY)
      call metrics_finalize(m_off)
      call metrics_apply_land_mask(m_off, wm, g, &
                                   periodic_x=.true., periodic_y=.false., &
                                   north_fold=.false., &
                                   mask_wall_velocity=.false., &
                                   wall_west=.false., wall_east=.false., &
                                   wall_south=.true., wall_north=.true.)

      call check(error, all(m_off%wet_v == m_legacy%wet_v), &
                 "knob-off wet_v differs from legacy path")
      if (allocated(error)) return
      call check(error, all(m_off%wet_u == m_legacy%wet_u), &
                 "knob-off wet_u differs from legacy path")
      if (allocated(error)) return
      call check(error, all(m_off%wet_q == m_legacy%wet_q), &
                 "knob-off wet_q differs from legacy path")
      if (allocated(error)) return
      ! Interior (physical) wet_v on a flat all-wet channel is a literal 1
      ! (the wall ghosts are untouched ⇒ the wall v-face product stays 1).
      call check(error, all(m_off%wet_v(g%nghost + 1:g%nghost + g%nx_phys, &
                                        g%nghost + 1:g%nghost + g%ny_phys + 1) == 1.0_wp), &
                 "knob-off wall v-face masked (should be untouched ⇒ 1.0)")
   end subroutine test_wall_velocity_off

end module test_ocean_land_mask
