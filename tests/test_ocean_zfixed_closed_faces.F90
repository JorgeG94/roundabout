!! The partial-step z-level FACE-CLOSURE mask builder
!! (`&vcoord_nml zfixed_closed_faces`; Adcroft, Hill & Marshall 1997;
!! Losch 2008 §2.1 for the ice-shelf cavity).
!!
!! `ocean_vcoord_closed_face_masks` is a `pure` procedure over explicit
!! arrays, so it can be tested against a HAND-BUILT staircase with no
!! grid, no state and no engine — which is the point of owning it in the
!! vcoord module rather than inlining it into a `configure_*`.
!!
!! What each case pins:
!!
!!   * `staircase_three_columns` — the rule itself, on a 3-column bed
!!     staircase carrying a LEDGE (an isolated live cell whose four
!!     own-layer faces are all closed).  Every u- and v-face of the
!!     interior is enumerated by hand.
!!   * `periodic_wrap_seam` — the seam of a periodic axis is an INTERIOR
!!     index (`nghost+1`), not the array edge, so it is covered by the
!!     builder's `2:nx` sweep PROVIDED the caller filled the ghosts.
!!     This asserts exactly that: a ghost column carrying the wrapped
!!     partner's target closes the seam face iff the partner is a filler.
!!   * `array_edge_faces_stay_open` — `I = 1` / `I = nx+1` are left at 1,
!!     the same convention the porous kernel uses; their `dy_cu` is
!!     already zero and the continuity wall zeroing owns them.
!!   * `all_live_is_all_open` — a column with no fillers anywhere leaves
!!     every interior face open, i.e. the mask is inert where the
!!     coordinate has no staircase (this is what makes the knob
!!     bit-identical on a deep-water configuration).
!!   * `mask_matches_the_z_fixed_target` — the mask is built from the
!!     SAME `ocean_vcoord_z_fixed_target` the ALE regrid and the IC seed
!!     use, so "live" has one definition.  A stepped bed under a flat lid
!!     is laid with that kernel and the mask cross-checked against the
!!     target it came from.
!!   * `refuses_correction_bc_pgf` / `refuses_substep_drag` /
!!     `refuses_wave_drag` — the three barotropic paths that still weight
!!     by the FULL column (`compute_pbce` + `compute_gtot_faces` + the
!!     bc-PGF `du_bc` block; `compute_bt_rem`; `compute_bt_rem_wave_drag`)
!!     are FAIL-LOUD at configure under the knob, and the SAME request
!!     with the knob off is accepted (so the refusal is the knob's, not
!!     the path's).
module test_ocean_zfixed_closed_faces
   use rdb_constants, only: wp, H_VANISHED
   use rdb_ocean_vcoord, only: ocean_vcoord_closed_face_masks, &
                               ocean_vcoord_count_ledges, &
                               ocean_vcoord_z_fixed_target
   use rdb_ocean_vcoord, only: VCOORD_Z_FIXED
   use rdb_config, only: config_t, read_config_from_string
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_closed_faces
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: error_ring_get, error_ring_clear
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_zfixed_closed_faces_tests

   real(wp), parameter :: FILLER = 1.0e-4_wp
      !! `zstar_h_min` at its default — an inert filler, below H_VANISHED.
   real(wp), parameter :: LIVE = 50.0_wp
      !! An ordinary live layer.

contains

   subroutine collect_ocean_zfixed_closed_faces_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("staircase_three_columns", test_staircase), &
                  new_unittest("periodic_wrap_seam", test_periodic_seam), &
                  new_unittest("array_edge_faces_stay_open", test_array_edges), &
                  new_unittest("all_live_is_all_open", test_all_live), &
                  new_unittest("mask_matches_the_z_fixed_target", test_from_target), &
                  new_unittest("refuses_correction_bc_pgf", test_refuses_bc_pgf), &
                  new_unittest("refuses_substep_drag", test_refuses_substep_drag), &
                  new_unittest("refuses_wave_drag", test_refuses_wave_drag) &
                  ]
   end subroutine collect_ocean_zfixed_closed_faces_tests

   subroutine test_staircase(error)
      !! A hand-built 3-column × 3-row × 3-layer staircase.
      !!
      !! Liveness by layer (`k = 1` is the BED, `k = 3` the top), on the
      !! middle row `j = 2`:
      !!
      !! ```
      !!            i = 1      i = 2      i = 3
      !!   k = 3    live       live       live
      !!   k = 2    FILLER     live       live      <- the step
      !!   k = 1    FILLER     FILLER     live      <- the deeper step
      !! ```
      !!
      !! plus a LEDGE at `(2,2,1)`: that cell is made live while all four
      !! of its own-layer neighbours stay fillers, so the mask isolates
      !! it.  Rows `j = 1` and `j = 3` are all-filler at `k <= 2`, which
      !! is what closes the ledge's v-faces.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3, NZ = 3
      real(wp) :: tgt(NX, NY, NZ)
      real(wp) :: ou(NX + 1, NY, NZ), ov(NX, NY + 1, NZ)
      integer :: n_ledge

      ! Everything a filler, then open up the staircase.
      tgt = FILLER
      ! k = 3 (top): live everywhere.
      tgt(:, :, 3) = LIVE
      ! k = 2: live from i = 2 eastwards, on every row.
      tgt(2:3, :, 2) = LIVE
      ! k = 1 (bed): live only in the easternmost column.
      tgt(3, :, 1) = LIVE
      ! The LEDGE: (2,2,1) is live, but (1,2,1), (3,2,1) is live... so to
      ! isolate it we make its EAST neighbour a filler too.  Redo the bed
      ! row: only (2,2,1) is live on j = 2, and nothing on j = 1 / j = 3.
      tgt(:, :, 1) = FILLER
      tgt(2, 2, 1) = LIVE

      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NX, NY, NZ, H_VANISHED)

      ! ---- k = 3: every interior face open ----
      call check(error, ou(2, 2, 3) == 1.0_wp, "k=3 u-face I=2 must be open")
      if (allocated(error)) return
      call check(error, ou(3, 2, 3) == 1.0_wp, "k=3 u-face I=3 must be open")
      if (allocated(error)) return
      call check(error, ov(2, 2, 3) == 1.0_wp, "k=3 v-face J=2 must be open")
      if (allocated(error)) return
      call check(error, ov(2, 3, 3) == 1.0_wp, "k=3 v-face J=3 must be open")
      if (allocated(error)) return

      ! ---- k = 2: the step.  Face I=2 separates the FILLER column i=1
      ! from the live column i=2 => CLOSED.  Face I=3 separates two live
      ! columns => open.  Every v-face is live/live => open.
      call check(error, ou(2, 2, 2) == 0.0_wp, &
                 "k=2 u-face I=2 straddles a filler (i=1) and must be CLOSED")
      if (allocated(error)) return
      call check(error, ou(3, 2, 2) == 1.0_wp, &
                 "k=2 u-face I=3 has water on both sides and must be open")
      if (allocated(error)) return
      call check(error, ov(2, 2, 2) == 1.0_wp, "k=2 v-face J=2 must be open")
      if (allocated(error)) return
      call check(error, ov(3, 3, 2) == 1.0_wp, "k=2 v-face J=3 at i=3 must be open")
      if (allocated(error)) return
      ! i = 1 is a filler on every row at k = 2, so its own v-faces close.
      call check(error, ov(1, 2, 2) == 0.0_wp, &
                 "k=2 v-face J=2 at i=1 is filler/filler and must be CLOSED")
      if (allocated(error)) return

      ! ---- k = 1: the ledge.  (2,2,1) is the ONLY live cell, so all
      ! four of its own-layer faces are closed.
      call check(error, ou(2, 2, 1) == 0.0_wp, "ledge west u-face must be CLOSED")
      if (allocated(error)) return
      call check(error, ou(3, 2, 1) == 0.0_wp, "ledge east u-face must be CLOSED")
      if (allocated(error)) return
      call check(error, ov(2, 2, 1) == 0.0_wp, "ledge south v-face must be CLOSED")
      if (allocated(error)) return
      call check(error, ov(2, 3, 1) == 0.0_wp, "ledge north v-face must be CLOSED")
      if (allocated(error)) return

      ! ...and the ledge census must SEE it (exactly one, at (2,2,1)).
      n_ledge = ocean_vcoord_count_ledges(ou, ov, tgt, NX, NY, NZ, H_VANISHED)
      call check(error, n_ledge == 1, "the ledge census must report exactly 1 isolated cell")
      if (allocated(error)) return

      ! A closed face is EXACTLY zero and an open one EXACTLY one — the
      ! mask is multiplied into transports, so anything in between would
      ! be a silent partial leak rather than a wall.
      call check(error, all(ou == 0.0_wp .or. ou == 1.0_wp), &
                 "open_u must be exactly 0 or exactly 1")
      if (allocated(error)) return
      call check(error, all(ov == 0.0_wp .or. ov == 1.0_wp), &
                 "open_v must be exactly 0 or exactly 1")
   end subroutine test_staircase

   subroutine test_periodic_seam(error)
      !! PERIODIC wrap.  `nghost = 1`, `nx_phys = 3`, so the array is
      !! `nx_total = 5` and the physical seam u-face is `I = nghost+1 = 2`
      !! — an INTERIOR index, covered by the builder's `2:nx` sweep.  The
      !! west ghost column `i = 1` carries the EAST physical column's
      !! target (that is what `ocean_periodic_wrap_centre_2d` +
      !! `ocean_halo_centre` leave behind before this runs).
      !!
      !! Case: the east physical column is a FILLER at `k = 1`, the west
      !! physical column is live.  Across the seam that face must CLOSE —
      !! which it only does if the ghost was filled.  The same face with
      !! a live ghost must stay open.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 1, NZ = 2
      real(wp) :: tgt(NX, NY, NZ)
      real(wp) :: ou(NX + 1, NY, NZ), ov(NX, NY + 1, NZ)

      ! Physical columns are i = 2,3,4; ghosts are i = 1 and i = 5.
      tgt = LIVE
      ! East physical column (i = 4) is a filler at the bed.
      tgt(4, 1, 1) = FILLER
      ! Periodic wrap: west ghost i = 1 mirrors the east physical i = 4.
      tgt(1, 1, 1) = tgt(4, 1, 1)
      ! East ghost i = 5 mirrors the west physical i = 2.
      tgt(5, 1, 1) = tgt(2, 1, 1)

      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NX, NY, NZ, H_VANISHED)

      call check(error, ou(2, 1, 1) == 0.0_wp, &
                 "the periodic SEAM face (I = nghost+1) must close when the "// &
                 "wrapped partner is a filler — if this passes as 1 the ghost "// &
                 "was not filled before the mask was built")
      if (allocated(error)) return
      call check(error, ou(2, 1, 2) == 1.0_wp, &
                 "the seam face must stay open at a layer that is live on both sides")
      if (allocated(error)) return
      ! The east ghost is live, so the face at I = 5 (between i=4 filler
      ! and i=5 live) also closes: the mask is symmetric about the seam.
      call check(error, ou(5, 1, 1) == 0.0_wp, &
                 "the east side of the seam must agree with the west side")
   end subroutine test_periodic_seam

   subroutine test_array_edges(error)
      !! Array-edge faces are left fully OPEN (1), exactly as the porous
      !! kernel leaves them: `dy_cu` is already zero there and the
      !! continuity wall zeroing owns them.  Pinning this stops a future
      !! "tidy up the edges" change from silently making the mask the
      !! only thing holding a wall shut.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 2, NZ = 2
      real(wp) :: tgt(NX, NY, NZ)
      real(wp) :: ou(NX + 1, NY, NZ), ov(NX, NY + 1, NZ)

      tgt = FILLER
      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NX, NY, NZ, H_VANISHED)

      call check(error, all(ou(1, :, :) == 1.0_wp), "u-face I=1 must stay open")
      if (allocated(error)) return
      call check(error, all(ou(NX + 1, :, :) == 1.0_wp), "u-face I=nx+1 must stay open")
      if (allocated(error)) return
      call check(error, all(ov(:, 1, :) == 1.0_wp), "v-face J=1 must stay open")
      if (allocated(error)) return
      call check(error, all(ov(:, NY + 1, :) == 1.0_wp), "v-face J=ny+1 must stay open")
      if (allocated(error)) return
      ! ...while every INTERIOR face of an all-filler stack is closed.
      call check(error, all(ou(2:NX, :, :) == 0.0_wp), &
                 "every interior u-face of an all-filler stack must be CLOSED")
      if (allocated(error)) return
      call check(error, all(ov(:, 2:NY, :) == 0.0_wp), &
                 "every interior v-face of an all-filler stack must be CLOSED")
   end subroutine test_array_edges

   subroutine test_all_live(error)
      !! No fillers anywhere ⇒ no closed faces ⇒ the mask is the identity
      !! and every consumer's multiply is a multiply by 1.  This is the
      !! structural reason a deep-water `z_fixed` configuration is
      !! unaffected by the knob.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 4, NZ = 5
      real(wp) :: tgt(NX, NY, NZ)
      real(wp) :: ou(NX + 1, NY, NZ), ov(NX, NY + 1, NZ)

      tgt = LIVE
      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NX, NY, NZ, H_VANISHED)
      call check(error, all(ou == 1.0_wp), "an all-live stack must leave open_u ≡ 1")
      if (allocated(error)) return
      call check(error, all(ov == 1.0_wp), "an all-live stack must leave open_v ≡ 1")
      if (allocated(error)) return
      call check(error, ocean_vcoord_count_ledges(ou, ov, tgt, NX, NY, NZ, H_VANISHED) == 0, &
                 "an all-live stack must isolate nothing")
   end subroutine test_all_live

   subroutine test_from_target(error)
      !! End-to-end against the REAL target kernel: a stepped bed under
      !! no lid, laid by `ocean_vcoord_z_fixed_target` exactly as
      !! `configure_ocean_closed_faces` lays it (η = 0, `z_top = 0`).
      !!
      !! Geometry: `h_nominal = 100 m`, `nz = 4` ⇒ nominal interfaces at
      !! 0/100/200/300/400 m.  Column depths 400, 250, 120 m.  The 400 m
      !! column is fully live; the 250 m column loses its bed layer; the
      !! 120 m column loses two.  The mask must then close the bed face
      !! between columns 1 and 2, and the two bed-side faces between 2
      !! and 3 — and leave the surface layer open everywhere.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 1, NZ = 4
      real(wp), parameter :: H_NOM = 100.0_wp
      real(wp) :: tgt(NX, NY, NZ), total_h(NX, NY), eta(NX, NY), z_top(NX, NY)
      real(wp) :: ou(NX + 1, NY, NZ), ov(NX, NY + 1, NZ)
      integer :: k
      logical :: live1, live2, live3

      total_h(1, 1) = 400.0_wp
      total_h(2, 1) = 250.0_wp
      total_h(3, 1) = 120.0_wp
      eta = 0.0_wp
      z_top = 0.0_wp
      call ocean_vcoord_z_fixed_target(tgt, total_h, eta, z_top, &
                                       NX, NY, NZ, H_NOM, FILLER)
      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NX, NY, NZ, H_VANISHED)

      ! The target must have produced the staircase this case is about.
      call check(error, tgt(1, 1, 1) > H_VANISHED, &
                 "the 400 m column's bed layer must be live (check the target, not the mask)")
      if (allocated(error)) return
      call check(error, tgt(3, 1, 1) <= H_VANISHED, &
                 "the 120 m column's bed layer must be a filler")
      if (allocated(error)) return

      ! Cross-check every interior face against the liveness of its two
      ! columns, read from the SAME target.  This is the invariant, not a
      ! transcription of the expected answer.
      do k = 1, NZ
         live1 = tgt(1, 1, k) > H_VANISHED
         live2 = tgt(2, 1, k) > H_VANISHED
         live3 = tgt(3, 1, k) > H_VANISHED
         call check(error, (ou(2, 1, k) == 1.0_wp) .eqv. (live1 .and. live2), &
                    "u-face I=2 must be open exactly when both columns are live")
         if (allocated(error)) return
         call check(error, (ou(3, 1, k) == 1.0_wp) .eqv. (live2 .and. live3), &
                    "u-face I=3 must be open exactly when both columns are live")
         if (allocated(error)) return
      end do

      ! The surface layer is live in every column, so its faces are open.
      call check(error, ou(2, 1, NZ) == 1.0_wp .and. ou(3, 1, NZ) == 1.0_wp, &
                 "the surface layer must stay open across the whole staircase")
      if (allocated(error)) return
      ! ...and the bed layer's face into the shallowest column is closed.
      call check(error, ou(3, 1, 1) == 0.0_wp, &
                 "the bed face into the 120 m column must be CLOSED")
   end subroutine test_from_target

   subroutine test_refuses_bc_pgf(error)
      !! `correction_bc_pgf`: `compute_pbce`, `compute_gtot_faces` and the
      !! bc-PGF `du_bc` block weight by the FULL column, so the
      !! correction's depth-mean-zero identity is not the OPEN column's.
      type(error_type), allocatable, intent(out) :: error
      call check_refusal(error, "&ocean_bt_nml correction_bc_pgf = .true. /", &
                         "correction_bc_pgf")
   end subroutine test_refuses_bc_pgf

   subroutine test_refuses_substep_drag(error)
      !! `substep_drag`: `compute_bt_rem` damps on the FULL-column depth.
      type(error_type), allocatable, intent(out) :: error
      call check_refusal(error, "&ocean_bt_nml substep_drag = .true. /"// &
                         new_line("a")//'&ocean_bdrag_nml form = "linear", '// &
                         "r = 1.0e-4, hbbl = 10.0 /", "substep_drag")
   end subroutine test_refuses_substep_drag

   subroutine test_refuses_wave_drag(error)
      !! `wave_drag`: `compute_bt_rem_wave_drag` damps on the FULL-column
      !! depth.
      type(error_type), allocatable, intent(out) :: error
      call check_refusal(error, "&ocean_bt_nml wave_drag = .true., "// &
                         "wave_drag_r_uniform = 1.0e-3 /", "wave_drag")
   end subroutine test_refuses_wave_drag

   subroutine check_refusal(error, extra, knob)
      !! Configure the closed faces on a bare state that passes every
      !! EARLIER guard of `configure_ocean_closed_faces` (multilayer
      !! initialised, `z_fixed`, a resolved `z_fixed_h_ref`, no BT
      !! workstate yet) with `extra` on top, and assert:
      !!
      !!   1. with `zfixed_closed_faces = .true.` the call is REFUSED with
      !!      `OCEAN_STATUS_ERR_SETUP` and the message names `knob`;
      !!   2. with the knob off the same request is accepted (the routine
      !!      returns before touching the grid, which is why a bare
      !!      `hgrid_t` suffices for both legs).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: extra, knob
      type(config_t) :: cfg
      type(ocean_state_t) :: st
      type(hgrid_t) :: grid
      integer :: ierr
      character(len=:), allocatable :: msg

      st%multilayer%is_init = .true.
      st%vcoord%coord_type = VCOORD_Z_FIXED
      st%vcoord%z_fixed_h_ref = 400.0_wp

      call parse_case(cfg, extra, .true.)
      call error_ring_clear()
      call configure_ocean_closed_faces(cfg, st, grid, 1, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                 "zfixed_closed_faces with "//knob//" must be REFUSED at configure")
      if (allocated(error)) return
      msg = trim(error_ring_get(0))
      call check(error, index(msg, "zfixed_closed_faces") > 0 .and. &
                 index(msg, knob) > 0, &
                 "the refusal must name both zfixed_closed_faces and "//knob// &
                 "; got: "//msg)
      if (allocated(error)) return
      call check(error,.not. st%metrics%use_closed_faces, &
                 "a refused configure must not latch use_closed_faces")
      if (allocated(error)) return

      call parse_case(cfg, extra, .false.)
      call configure_ocean_closed_faces(cfg, st, grid, 1, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 knob//" without zfixed_closed_faces must be accepted here")
   end subroutine check_refusal

   subroutine parse_case(cfg, extra, closed)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: extra
      logical, intent(in) :: closed
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            '&vcoord_nml vcoord_type = "z_fixed", zfixed_closed_faces = '// &
            merge(".true. ", ".false.", closed)//" /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
   end subroutine parse_case

end module test_ocean_zfixed_closed_faces
