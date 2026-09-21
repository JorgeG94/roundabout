!! `k_top` — the shared index of the first LIVE layer, counting down
!! from the top (P6.3).
!!
!! Under a quasi-geopotential coordinate beneath an ice shelf
!! (`vcoord_type = "z_fixed"` with `vcoord%z_top > 0`) the layers whose
!! nominal range lies inside the draft are inert fillers at
!! `zstar_h_min`, so on an ice-covered column `k = nz` is NOT the
!! ice-adjacent layer.  `ocean_vcoord_k_top_from_target` is the ONE
!! producer of the index every top-side consumer reads instead.
!!
!! Three things are asserted here, and the third is the one that keeps
!! the tree honest:
!!
!!   1. **Fallback.** No top-side filler ⇒ `k_top ≡ nz` on centres AND
!!      on both face staggers.  That is what makes routing every
!!      consumer through `k_top` bit-identical on sigma / z*-lite / every
!!      family shipped today.
!!   2. **Placement.** Under a rigid top the index lands on the PARTIAL
!!      TOP CELL — the layer the ice base cuts — with every layer above
!!      it at or below `H_VANISHED`.
!!   3. **Agreement with the closed-face mask.** Both are built from the
!!      same `z_fixed` target at `eta = 0` by the same definition of
!!      "live", so `open_u(i,j,k)` must be exactly zero for every
!!      `k > k_top_u(i,j)`, and the face index must be the `min` of its
!!      two columns — never the `max`, which is the bug that would put
!!      the ice-ocean drag on a row that is a filler on one side.
module test_ocean_ktop
   use, intrinsic :: iso_fortran_env, only: real64
   use rdb_constants, only: wp, H_VANISHED
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_k_top_from_target, &
                               ocean_vcoord_closed_face_masks
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_ktop_tests

   ! A 1000 m bed, nz = 20 ⇒ 50 m nominal spacing, filler 1.0e-4 m
   ! (strictly below H_VANISHED = 1.5e-4, the inert-filler contract).
   integer, parameter :: NZ = 20
   integer, parameter :: NX = 6
   integer, parameter :: NY = 4
   real(wp), parameter :: H_REF = 1000.0_wp
   real(wp), parameter :: H_NOM = H_REF/real(NZ, wp)
   real(wp), parameter :: H_MIN = 1.0e-4_wp

contains

   subroutine collect_ocean_ktop_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("no_filler_falls_back_to_nz", test_fallback_nz), &
                  new_unittest("flat_lid_lands_on_the_partial_top_cell", test_flat_lid), &
                  new_unittest("faces_take_the_min_of_their_two_columns", test_face_min), &
                  new_unittest("agrees_with_the_closed_face_mask", test_mask_agreement), &
                  new_unittest("land_column_reads_nz", test_land_column) &
                  ]
   end subroutine collect_ocean_ktop_tests

   subroutine build(z_top_val, total_h_val, tgt, k_top, k_top_u, k_top_v)
      !! Uniform geometry over the whole array: one `z_top`, one column
      !! thickness, `eta = 0` (the configure-time datum `k_top` is built
      !! from in production).
      real(wp), intent(in) :: z_top_val, total_h_val
      real(wp), intent(out) :: tgt(NX, NY, NZ)
      integer, intent(out) :: k_top(NX, NY)
      integer, intent(out) :: k_top_u(NX + 1, NY)
      integer, intent(out) :: k_top_v(NX, NY + 1)
      real(wp) :: total_h(NX, NY), eta(NX, NY), z_top(NX, NY)

      total_h = total_h_val
      eta = 0.0_wp
      z_top = z_top_val
      call ocean_vcoord_z_fixed_target(tgt, total_h, eta, z_top, &
                                       NX, NY, NZ, H_NOM, H_MIN)
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, &
                                          NX, NY, NZ, H_VANISHED)
   end subroutine build

   subroutine test_fallback_nz(error)
      !! No rigid top: the whole 1000 m column is live, nothing vanishes
      !! against the surface, and every index reads `nz`.  This is the
      !! bit-identity gate for the entire slice — every consumer rewritten
      !! from `nz` to `k_top(i,j)` reads the same memory here.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(NX, NY, NZ)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)

      call build(0.0_wp, H_REF, tgt, k_top, k_top_u, k_top_v)

      call check(error, all(k_top == NZ), "k_top is nz with no rigid top")
      if (allocated(error)) return
      call check(error, all(k_top_u == NZ), "k_top_u is nz with no rigid top")
      if (allocated(error)) return
      call check(error, all(k_top_v == NZ), "k_top_v is nz with no rigid top")
      if (allocated(error)) return
      call check(error, tgt(3, 2, NZ) > H_VANISHED, &
                 "and the top layer really is live")
   end subroutine test_fallback_nz

   subroutine test_flat_lid(error)
      !! A 500 m flat lid over a 1000 m bed: the nominal 50 m stack is
      !! cut at the ice base, layers 11..20 collapse to the inert filler
      !! and layer 10 is the partial top cell.  `k_top` must name layer
      !! 10 — the first layer that carries mass — and every layer above
      !! it must be at or below the vanish marker.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(NX, NY, NZ)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)
      integer :: k
      logical :: ok

      call build(500.0_wp, 500.0_wp, tgt, k_top, k_top_u, k_top_v)

      call check(error, all(k_top == 10), "k_top names the partial top cell")
      if (allocated(error)) return
      call check(error, tgt(3, 2, 10) > H_VANISHED, "the partial top cell is live")
      if (allocated(error)) return

      ok = .true.
      do k = 11, NZ
         if (tgt(3, 2, k) > H_VANISHED) ok = .false.
      end do
      call check(error, ok, "every layer above k_top is an inert filler")
      if (allocated(error)) return

      ! The whole point: the column sum is still exact, so the fillers
      ! are accounted for and nothing has been invented or lost.
      call check(error, abs(sum(tgt(3, 2, :)) - 500.0_wp) < 1.0e-10_wp, &
                 "Sum target_h = H + eta with the fillers included")
   end subroutine test_flat_lid

   subroutine test_face_min(error)
      !! Two columns whose drafts straddle a nominal level: one carries
      !! more fillers than the other, so their `k_top` differ by one.
      !! A velocity face carries water in layer `k` only where BOTH sides
      !! do, so the face index is the DEEPER of the two tops — the
      !! smaller index, `min`.  `max` would hand the ice-ocean drag and
      !! the implicit stress fold a row that is a filler on one side.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(NX, NY, NZ), total_h(NX, NY), eta(NX, NY), z_top(NX, NY)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)
      integer :: i, j
      logical :: ok

      eta = 0.0_wp
      ! A staircase in i: deeper draft to the east ⇒ fewer live layers.
      do j = 1, NY
         do i = 1, NX
            z_top(i, j) = 100.0_wp*real(i - 1, wp)
            total_h(i, j) = H_REF - z_top(i, j)
         end do
      end do
      call ocean_vcoord_z_fixed_target(tgt, total_h, eta, z_top, &
                                       NX, NY, NZ, H_NOM, H_MIN)
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, &
                                          NX, NY, NZ, H_VANISHED)

      ! The staircase is real (otherwise the test proves nothing).
      call check(error, k_top(1, 2) > k_top(NX, 2), &
                 "the draft staircase does move k_top")
      if (allocated(error)) return

      ok = .true.
      do j = 1, NY
         do i = 2, NX
            if (k_top_u(i, j) /= min(k_top(i - 1, j), k_top(i, j))) ok = .false.
            if (k_top_u(i, j) > k_top(i - 1, j)) ok = .false.
            if (k_top_u(i, j) > k_top(i, j)) ok = .false.
         end do
      end do
      call check(error, ok, "k_top_u is the min of its two columns")
      if (allocated(error)) return

      ! The v faces run along the draft contour here, so both columns
      ! agree and the min is that common value — still the min rule.
      ok = .true.
      do j = 2, NY
         do i = 1, NX
            if (k_top_v(i, j) /= min(k_top(i, j - 1), k_top(i, j))) ok = .false.
         end do
      end do
      call check(error, ok, "k_top_v is the min of its two columns")
   end subroutine test_face_min

   subroutine test_mask_agreement(error)
      !! `k_top` and `metrics%open_u/open_v` are two readings of ONE
      !! live/filler pattern, built from the same target by the same
      !! strict `> H_VANISHED` test.  If they ever disagree, the
      !! closed-face wall and the forcing row are talking about different
      !! layers.  Asserted on the staircase geometry, where the two
      !! columns of a face genuinely differ.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(NX, NY, NZ), total_h(NX, NY), eta(NX, NY), z_top(NX, NY)
      real(wp) :: open_u(NX + 1, NY, NZ), open_v(NX, NY + 1, NZ)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)
      integer :: i, j, k
      logical :: ok

      eta = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            z_top(i, j) = 100.0_wp*real(i - 1, wp)
            total_h(i, j) = H_REF - z_top(i, j)
         end do
      end do
      call ocean_vcoord_z_fixed_target(tgt, total_h, eta, z_top, &
                                       NX, NY, NZ, H_NOM, H_MIN)
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, &
                                          NX, NY, NZ, H_VANISHED)
      call ocean_vcoord_closed_face_masks(open_u, open_v, tgt, &
                                          NX, NY, NZ, H_VANISHED)

      ! Above the face's first live layer the mask MUST be closed.
      ok = .true.
      do k = 1, NZ
         do j = 1, NY
            do i = 2, NX
               if (k > k_top_u(i, j) .and. open_u(i, j, k) /= 0.0_wp) ok = .false.
            end do
         end do
      end do
      call check(error, ok, "open_u is closed above k_top_u")
      if (allocated(error)) return

      ! And AT it the mask must be open — that is what "first live layer
      ! of the face" means.  (Both columns are far thicker than the bed
      ! filler band on this geometry, so the overlap is non-empty.)
      ok = .true.
      do j = 1, NY
         do i = 2, NX
            if (open_u(i, j, k_top_u(i, j)) /= 1.0_wp) ok = .false.
         end do
      end do
      call check(error, ok, "open_u is open AT k_top_u")
      if (allocated(error)) return

      ok = .true.
      do k = 1, NZ
         do j = 2, NY
            do i = 1, NX
               if (k > k_top_v(i, j) .and. open_v(i, j, k) /= 0.0_wp) ok = .false.
            end do
         end do
      end do
      call check(error, ok, "open_v is closed above k_top_v")
   end subroutine test_mask_agreement

   subroutine test_land_column(error)
      !! A column with no live layer at all — the land-state contract
      !! holds every layer AT the marker — falls back to `nz`, which is
      !! the row those consumers index today and which `wet_mask` then
      !! zeroes.  A fallback of 0 would be an out-of-bounds write on the
      !! GPU with no diagnostic.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(NX, NY, NZ)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)

      tgt = H_VANISHED
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, &
                                          NX, NY, NZ, H_VANISHED)
      call check(error, all(k_top == NZ), "a dead column falls back to nz")
      if (allocated(error)) return
      call check(error, all(k_top_u == NZ) .and. all(k_top_v == NZ), &
                 "and so do its faces")
   end subroutine test_land_column

end module test_ocean_ktop
