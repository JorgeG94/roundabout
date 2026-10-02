!! The partial-step face closure (`&vcoord_nml zfixed_closed_faces`) on
!! `vcoord_type = "zstar_full"`.
!!
!! `build_zref_full` lays a z-level FINE zone (`zstar_h_surf_target`-thick
!! layers from the surface, `nz/3` of them on the ocean path) over a
!! terrain-following COARSE zone.  A column shallower than the fine zone
!! ends in a partial cell and every layer below it is a `zstar_h_min`
!! FILLER — the `z_fixed` bed staircase, with the same
!! thickness-independent FV-PGF defect across it (the 1-degree Southern
!! Ocean dies at step 10 of it, `vcoord_audit` H2).  The knob builds the
!! same static mask from the `ZSTAR_FULL` target at `eta = 0` and reuses
!! every consumer.
!!
!! What each case pins:
!!
!!   * `mask_matches_the_zstar_full_target` — the mask is built from the
!!     SAME target the ALE regrid computes (`ocean_vcoord_eta0_target` ->
!!     the `ZSTAR_FULL` branch of the geometric kernel, over the `z_ref`
!!     table `build_zref_full` lays): every interior face is open exactly
!!     when both columns are live, the fine zone stays open across a
!!     shelf break and the coarse layers facing a filler are closed.
!!   * `pattern_static_for_positive_eta` — the property that makes a
!!     configure-time mask sound: under `eta >= 0` the `ZSTAR_FULL`
!!     target puts the whole anomaly in the surface layer, so no layer
!!     changes between live and filler.
!!   * `refuses_without_a_fine_zone` / `still_refuses_sigma` — the
!!     generalised coordinate guard: `zstar_full` with
!!     `zstar_h_surf_target <= 0` (uniform sigma, nothing to close) and
!!     every non-geometric family stay fail-loud, naming the knob.
!!   * `resting_staircase_stays_at_rest` — END TO END through the C API:
!!     a stratified resting ocean over a 200 m / 30 m shelf break, seeded
!!     ON the coordinate.  With the mask the velocity stays at round-off
!!     (every open face pairs the same density stack on both sides; the
!!     coarse layers that face the shelf's fillers are walls).  The SAME
!!     state with the knob off (thickness and tracers copied across) is
!!     driven from rest by the staircase pressure gradient — the test
!!     fails without the mask.
module test_ocean_zstar_full_closed_faces
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer, c_double
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, ocean_vcoord_closed_face_masks, &
                               ocean_vcoord_eta0_target
   use rdb_constants, only: VCOORD_ZSTAR_FULL, VCOORD_SIGMA
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_closed_faces
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: error_ring_get, error_ring_clear
   use rdb_ocean_api, only: rdb_ocean_create_pending, rdb_ocean_create_finalize, &
                            rdb_ocean_stage_bathymetry, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_h_layer_ptr, &
                            rdb_ocean_get_tracer_ptr, &
                            rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, &
                            rdb_ocean_set_h, rdb_ocean_set_tracer
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_zstar_full_closed_faces_tests

   integer, parameter :: NZ = 6
      !! `nz/3 = 2` fine layers on the ocean path.
   real(wp), parameter :: H_SURF = 20.0_wp
      !! Fine-zone layer thickness ⇒ a 40 m fine zone.
   real(wp), parameter :: FILLER = 1.0e-4_wp
      !! `zstar_h_min`.
   real(wp), parameter :: DEEP = 200.0_wp
      !! Off-shelf depth: 2 x 20 m fine + 4 x 40 m coarse, no filler.
   real(wp), parameter :: SHELF = 30.0_wp
      !! Shelf depth: 20 m + a 10 m partial cell, 4 coarse fillers.
   integer(c_int), parameter :: NXP = 8, NYP = 4
      !! Physical grid of the end-to-end case (shelf break at i = 4|5).
   integer, parameter :: N_STEPS = 24
      !! Outer steps (2 hours at dt = 300 s).
   real(wp), parameter :: REST_TOL = 1.0e-12_wp
      !! Round-off ceiling on |u| with the mask (m/s).
   real(wp), parameter :: DRIVEN_MIN = 1.0e-4_wp
      !! Floor on |u| the knob-OFF twin must exceed (m/s).
   integer(c_int), parameter :: BATHY_DEPTH_POSITIVE_DOWN = 1

contains

   subroutine collect_ocean_zstar_full_closed_faces_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("mask_matches_the_zstar_full_target", test_from_target), &
                  new_unittest("pattern_static_for_positive_eta", test_static_eta), &
                  new_unittest("refuses_without_a_fine_zone", test_refuses_no_fine_zone), &
                  new_unittest("still_refuses_sigma", test_refuses_sigma), &
                  new_unittest("resting_staircase_stays_at_rest", test_resting_staircase) &
                  ]
   end subroutine collect_ocean_zstar_full_closed_faces_tests

   subroutine lay_three_columns(vc, grid, b)
      !! A 3-column shelf break (deep | shelf | very shallow shelf),
      !! ghosts copied from the edge columns, with the `zstar_full` slot
      !! configured exactly as `engine_setup` does and its `z_ref` table
      !! laid by `build_zref_full`.
      type(ocean_vcoord_t), intent(out) :: vc
      type(hgrid_t), intent(out) :: grid
      real(wp), allocatable, intent(out) :: b(:, :)
      integer :: i
      real(wp), parameter :: depth(3) = [DEEP, SHELF, 12.0_wp]

      call grid%init(3, 1, 1, 1000.0_wp, 1000.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = H_SURF
      vc%zstar_h_min = FILLER
      allocate (b(grid%nx_total, grid%ny_total))
      do i = 1, grid%nx_total
         b(i, :) = depth(min(3, max(1, i - 1)))
      end do
      call vc%build_zref_full(b)
   end subroutine lay_three_columns

   subroutine test_from_target(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), tgt(:, :, :), ou(:, :, :), ov(:, :, :)
      integer :: nx, ny, i, k
      logical :: lw, le
      character(len=120) :: msg

      call lay_three_columns(vc, grid, b)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (tgt(nx, ny, NZ), source=0.0_wp)
      allocate (ou(nx + 1, ny, NZ), ov(nx, ny + 1, NZ))
      call ocean_vcoord_eta0_target(vc, tgt, b, nx, ny, NZ)
      call ocean_vcoord_closed_face_masks(ou, ov, tgt, nx, ny, NZ, H_VANISHED)

      ! The target must have produced the staircase this case is about
      ! (physical columns are i = 2, 3, 4).
      call check(error, all(tgt(2, 1, :) > H_VANISHED), &
                 "the 200 m column must be live in every layer (check the target)")
      if (allocated(error)) return
      call check(error, all(tgt(3, 1, 1:4) <= H_VANISHED) .and. &
                 all(tgt(3, 1, 5:6) > H_VANISHED), &
                 "the 30 m column must carry 4 coarse fillers under 2 live fine layers")
      if (allocated(error)) return
      call check(error, all(tgt(4, 1, 1:5) <= H_VANISHED) .and. tgt(4, 1, 6) > H_VANISHED, &
                 "the 12 m column must be live in its surface layer only")
      if (allocated(error)) return
      do i = 2, 4
         write (msg, "(a,i0,a,es12.5)") "column ", i, " must sum to its depth; sum = ", &
            sum(tgt(i, 1, :))
         call check(error, abs(sum(tgt(i, 1, :)) - b(i, 1)) <= 1.0e-10_wp*b(i, 1), trim(msg))
         if (allocated(error)) return
      end do

      ! Every interior face against the liveness of its two columns —
      ! the invariant, not a transcription of the expected answer.
      do k = 1, NZ
         do i = 2, nx
            lw = tgt(i - 1, 1, k) > H_VANISHED
            le = tgt(i, 1, k) > H_VANISHED
            call check(error, (ou(i, 1, k) == 1.0_wp) .eqv. (lw .and. le), &
                       "a u-face must be open exactly when both columns are live")
            if (allocated(error)) return
         end do
      end do
      ! And the physics of it: the fine zone is open across the shelf
      ! break, the deep column's coarse layers are walls facing rock.
      call check(error, all(ou(3, 1, 5:6) == 1.0_wp), &
                 "the fine zone must stay open across the 200|30 m shelf break")
      if (allocated(error)) return
      call check(error, all(ou(3, 1, 1:4) == 0.0_wp), &
                 "the coarse layers facing the shelf's fillers must be CLOSED")
      if (allocated(error)) return
      call check(error, ou(4, 1, 6) == 1.0_wp .and. ou(4, 1, 5) == 0.0_wp, &
                 "30|12 m: the surface layer open, the 30 m column's partial cell closed")
   end subroutine test_from_target

   subroutine test_static_eta(error)
      !! `eta >= 0` lands in the surface layer: the live/filler pattern of
      !! the running target equals the configure-time one.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), tgt0(:, :, :), eta(:, :)
      integer :: n
      real(wp), parameter :: etas(3) = [1.0e-6_wp, 0.5_wp, 3.0_wp]

      call lay_three_columns(vc, grid, b)
      allocate (tgt0(grid%nx_total, grid%ny_total, NZ), source=0.0_wp)
      call ocean_vcoord_eta0_target(vc, tgt0, b, grid%nx_total, grid%ny_total, NZ)
      allocate (eta(grid%nx_total, grid%ny_total))
      do n = 1, size(etas)
         eta = etas(n)
         call vc%compute_target_h(b, eta)
         call check(error, all((vc%target_h > H_VANISHED) .eqv. (tgt0 > H_VANISHED)), &
                    "a positive eta must not flip any layer between live and filler")
         if (allocated(error)) return
      end do
   end subroutine test_static_eta

   subroutine test_refuses_no_fine_zone(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refusal(error, VCOORD_ZSTAR_FULL, 0.0_wp, "zstar_h_surf_target")
   end subroutine test_refuses_no_fine_zone

   subroutine test_refuses_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refusal(error, VCOORD_SIGMA, H_SURF, "zstar_full")
   end subroutine test_refuses_sigma

   subroutine check_refusal(error, coord, h_surf, needle)
      !! Configure the closed faces on a bare state with `coord` and
      !! `h_surf`, and assert the call is REFUSED naming both the knob and
      !! `needle`, without latching the mask.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: coord
      real(wp), intent(in) :: h_surf
      character(len=*), intent(in) :: needle
      type(config_t) :: cfg
      type(ocean_state_t) :: st
      type(hgrid_t) :: grid
      integer :: ierr
      character(len=:), allocatable :: msg, nml

      st%multilayer%is_init = .true.
      st%vcoord%coord_type = coord
      st%vcoord%zstar_h_surf_target = h_surf
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            '&vcoord_nml vcoord_type = "zstar_full", zfixed_closed_faces = .true. /'// &
            new_line("a")
      call read_config_from_string(nml, cfg)
      call error_ring_clear()
      call configure_ocean_closed_faces(cfg, st, grid, 1, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                 "zfixed_closed_faces must be REFUSED here (needle: "//needle//")")
      if (allocated(error)) return
      msg = trim(error_ring_get(0))
      call check(error, index(msg, "zfixed_closed_faces") > 0 .and. index(msg, needle) > 0, &
                 "the refusal must name zfixed_closed_faces and "//needle//"; got: "//msg)
      if (allocated(error)) return
      call check(error,.not. st%metrics%use_closed_faces, &
                 "a refused configure must not latch use_closed_faces")
   end subroutine check_refusal

   function staircase_nml(closed) result(nml)
      !! A resting, stably stratified (per coordinate layer, linear EOS)
      !! ocean over the staged 200 m / 30 m shelf break.  Every process
      !! that could move a resting column is off or inert: no wind, no
      !! buoyancy forcing, no rotation needed, no vertical mixing (it would
      !! diffuse the 20 m deep fine layer and the 10 m shelf partial cell
      !! differently and build a REAL horizontal gradient).
      logical, intent(in) :: closed
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 8, ny = 4, nghost = 2, dx = 2000.0, dy = 2000.0 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0, cfl_interval = 1 /"//nl// &
            "&physics_nml coriolis_f = 1.0e-4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 6 /"//nl// &
            "&vcoord_nml vcoord_type = 'zstar_full', zstar_h_surf_target = 20.0, "// &
            "zstar_h_min = 1.0e-4, check_vanished_content = .true., "// &
            "zfixed_closed_faces = "//merge(".true. ", ".false.", closed)//" /"//nl// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 200.0 /"//nl// &
            "&tracer_nml initial_salinity = 35.0, T_init_surface = 15.0, "// &
            "T_init_bottom = 5.0 /"//nl// &
            "&ocean_eos_nml eos = 'linear' /"//nl// &
            "&ocean_ic_nml alpha_T = 0.2, beta_S = 0.8, T_ref = 10.0, S_ref = 35.0, "// &
            "rho_0 = 1027.0 /"//nl// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny_energy' /"//nl// &
            "&ocean_hvisc_nml nu_h = 10.0 /"//nl// &
            "&ocean_bdrag_nml form = 'linear', r = 1.0e-4, hbbl = 10.0 /"//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false., "// &
            "pp81_nu_bg = 0.0, pp81_kappa_bg = 0.0 /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function staircase_nml

   subroutine create_staircase(handle, closed, status)
      type(c_ptr), intent(out) :: handle
      logical, intent(in) :: closed
      integer(c_int), intent(out) :: status
      character(len=:), allocatable :: nml
      real(c_double) :: bathy(NXP, NYP)

      bathy(1:4, :) = real(DEEP, c_double)
      bathy(5:8, :) = real(SHELF, c_double)
      nml = staircase_nml(closed)
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return
      status = rdb_ocean_stage_bathymetry(handle, bathy, NXP, NYP, BATHY_DEPTH_POSITIVE_DOWN)
      if (status /= OCEAN_STATUS_OK) return
      status = rdb_ocean_create_finalize(handle)
   end subroutine create_staircase

   function max_speed(handle) result(umax)
      !! max |u|, |v| over the physical interior faces, every layer.
      type(c_ptr), intent(in) :: handle
      real(wp) :: umax
      type(c_ptr) :: ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: u(:, :, :), v(:, :, :)
      status = rdb_ocean_refresh_host(handle)
      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, u, [nx, ny, nz])
      umax = maxval(abs(u(ng + 1:ng + nx_p + 1, ng + 1:ng + ny_p, :)))
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, v, [nx, ny, nz])
      umax = max(umax, maxval(abs(v(ng + 1:ng + nx_p, ng + 1:ng + ny_p + 1, :))))
   end function max_speed

   subroutine test_resting_staircase(error)
      !! The C API holds ONE ocean at a time, so the two legs run in
      !! sequence: the masked leg first (its on-coordinate seed is read
      !! back before stepping), then the knob-off twin, overwritten with
      !! that seed so the two runs differ in the mask only.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: on, off, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :)
      real(c_double), allocatable :: h_seed(:, :, :), t_seed(:, :, :), s_seed(:, :, :)
      real(wp) :: u_on, u_off
      integer :: i, j, k, kb
      character(len=200) :: msg

      ! ---- Leg 1: zstar_full + the mask ----
      call create_staircase(on, .true., status)
      call check(error, status == OCEAN_STATUS_OK, &
                 "zstar_full + zfixed_closed_faces must build an ocean")
      if (.not. allocated(error)) then
         leg_on: block
            status = rdb_ocean_refresh_host(on)
            status = rdb_ocean_get_grid_info(on, nx_p, ny_p, nz_p, ng)
            status = rdb_ocean_get_h_layer_ptr(on, ptr, nx, ny, nz, gen)
            call c_f_pointer(ptr, h, [nx, ny, nz])
            ! Explicit allocate + element copy, not an allocate-on-assign
            ! from the pointer section: nvfortran (25.5, -fast) built
            ! `h_seed` with a scrambled layout from that form — the same
            ! 64 fillers, at the wrong (i, j, k).
            allocate (h_seed(nx_p, ny_p, nz_p))
            do k = 1, nz_p
               do j = 1, ny_p
                  do i = 1, nx_p
                     h_seed(i, j, k) = real(h(ng + i, ng + j, k), c_double)
                  end do
               end do
            end do
            call check(error, count(h_seed <= H_VANISHED) == (NXP/2)*NYP*4, &
                       "the on-target seed must give every shelf column 4 coarse fillers")
            if (allocated(error)) exit leg_on
            ! The stratification, PER COORDINATE LAYER (5, 7, ..., 15 degC
            ! bed to surface), with every filler holding its DONOR's value
            ! (I1′).  Written here rather than taken from `&tracer_nml`:
            ! a per-index IC gives the shelf's coarse fillers 5-11 degC of
            ! FOREIGN content, which the seed's I1′ pooling then folds into
            ! the live partial cell above (12.9998 instead of 13 degC, a
            ! real density step at an OPEN face).  This state has no
            ! horizontal density difference at any open face.
            allocate (t_seed(nx_p, ny_p, nz_p), s_seed(nx_p, ny_p, nz_p))
            s_seed = 35.0_c_double
            ! The fillers sit at the bed, so a filler's donor is the
            ! column's lowest LIVE layer `kb`: T(k) = T_index(max(k, kb)).
            do j = 1, ny_p
               do i = 1, nx_p
                  kb = nz_p
                  do k = 1, nz_p
                     if (h_seed(i, j, k) > H_VANISHED) then
                        kb = k
                        exit
                     end if
                  end do
                  do k = 1, nz_p
                     t_seed(i, j, k) = 5.0_c_double + 2.0_c_double*real(max(k, kb) - 1, c_double)
                  end do
               end do
            end do

            status = rdb_ocean_set_tracer(on, "temperature", 11_c_int, t_seed, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, "set_tracer(T) on the masked leg")
            if (allocated(error)) exit leg_on
            status = rdb_ocean_set_tracer(on, "salinity", 8_c_int, s_seed, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, "set_tracer(S) on the masked leg")
            if (allocated(error)) exit leg_on

            status = rdb_ocean_step(on, int(N_STEPS, c_int))
            call check(error, status == OCEAN_STATUS_OK, "the masked run must step cleanly")
            if (allocated(error)) exit leg_on
            u_on = max_speed(on)
         end block leg_on
      end if
      status = rdb_ocean_destroy(on)
      if (allocated(error)) return

      ! ---- Leg 2: the same state, knob off ----
      call create_staircase(off, .false., status)
      call check(error, status == OCEAN_STATUS_OK, "the knob-off twin must build an ocean")
      if (.not. allocated(error)) then
         leg_off: block
            status = rdb_ocean_set_h(off, h_seed, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, "set_h on the knob-off twin")
            if (allocated(error)) exit leg_off
            status = rdb_ocean_set_tracer(off, "temperature", 11_c_int, t_seed, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, "set_tracer(T) on the knob-off twin")
            if (allocated(error)) exit leg_off
            status = rdb_ocean_set_tracer(off, "salinity", 8_c_int, s_seed, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, "set_tracer(S) on the knob-off twin")
            if (allocated(error)) exit leg_off
            status = rdb_ocean_step(off, int(N_STEPS, c_int))
            call check(error, status == OCEAN_STATUS_OK, "the knob-off twin must step")
            if (allocated(error)) exit leg_off
            u_off = max_speed(off)

            write (msg, "(a,es10.3,a,es10.3,a)") "max|u| with the mask = ", u_on, &
               " m/s, without = ", u_off, " m/s"
            call check(error, u_on <= REST_TOL, trim(msg)// &
                       " — the masked resting staircase must stay at rest")
            if (allocated(error)) exit leg_off
            call check(error, u_off >= DRIVEN_MIN, trim(msg)// &
                       " — without the mask the staircase PGF must drive the flow "// &
                       "(else this test proves nothing)")
         end block leg_off
      end if
      status = rdb_ocean_destroy(off)
   end subroutine test_resting_staircase

end module test_ocean_zstar_full_closed_faces
