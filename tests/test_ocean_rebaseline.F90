!! `rdb_ocean_rebaseline`: re-deriving the setup-time snapshots after a
!! post-create overwrite of the initial state.
!!
!! Cases:
!!   (A5) the SAME initial condition built (1) staged in the pending-create
!!        window and (2) by create + post-create `set_h/u/v/tracer` +
!!        `rebaseline` agrees on the wet cells: sponge reference (T, S,
!!        pseudo-salt, u, v), u, v, h, T, S, pseudo-salt, rho_layer, bt_eta;
!!        then 10 steps stay bit-identical.  The namelist enables the sponge
!!        (`target_source = "ic"`, strong enough to matter within 10 steps),
!!        pseudo-salt, an island and a periodic x axis, so each step of the
!!        rebaseline is exercised.  Without `rebaseline` the sponge target and
!!        the pseudo-salt are still the built-in uniform ones (checked, so the
!!        test cannot pass vacuously).
!!   (A8) `rebaseline` after a step fails `OCEAN_STATUS_ERR_SETUP`.
!!   (P)  a pseudo-salt the caller set itself survives the rebaseline; the
!!        flag is consumed by it.
!!
!! The compared fields are read from the live handle after `refresh_host`, so
!! on the GPU build the host values are the device ones; the 10-step
!! comparison is what proves the sponge reference and the masked state reached
!! the device.  Own ctest binary: single-live-handle module state.
module test_ocean_rebaseline
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_handle, only: ocean_handle_t, handle_check, HANDLE_OK
   use rdb_ocean_api, only: rdb_ocean_create_pending, rdb_ocean_create_finalize, &
                            rdb_ocean_stage_bathymetry, rdb_ocean_stage_h, rdb_ocean_stage_u, &
                            rdb_ocean_stage_v, rdb_ocean_stage_tracer, rdb_ocean_set_h, &
                            rdb_ocean_set_u, rdb_ocean_set_v, rdb_ocean_set_tracer, &
                            rdb_ocean_rebaseline, rdb_ocean_step, rdb_ocean_destroy, &
                            rdb_ocean_refresh_host
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   implicit none
   private

   public :: collect_ocean_rebaseline_tests

   integer, parameter :: NX = 10, NY = 8, NZ = 3, NG = 3
   real(wp), parameter :: DEPTH = 300.0_wp
   integer, parameter :: LI = 4, LJ = 3     ! the island cell
   real(wp), parameter :: ETA_SLOPE = 1.0e-3_wp

   type :: snap_t
      real(wp), allocatable :: h(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: hT(:, :, :), hS(:, :, :), hP(:, :, :)
      real(wp), allocatable :: rho(:, :, :), eta(:, :)
      real(wp), allocatable :: refT(:, :, :), refS(:, :, :), refP(:, :, :)
      real(wp), allocatable :: uref(:, :, :), vref(:, :, :)
      real(wp), allocatable :: wet(:, :), wet_u(:, :), wet_v(:, :)
   end type snap_t

contains

   subroutine collect_ocean_rebaseline_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("rebaseline_matches_staged_and_steps_bit_identically", &
                               test_matches_staged), &
                  new_unittest("rebaseline_refused_after_step", test_refused_after_step), &
                  new_unittest("rebaseline_keeps_user_set_pseudo_salt", test_user_pseudo_salt) &
                  ]
   end subroutine collect_ocean_rebaseline_tests

   function make_nml() result(nml)
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 10, ny = 8, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 300.0 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'sponge', sponge_width = 3, sponge_strength = 2.0e-3, "// &
            "sponge_relax_tracers = .true. /"//NL// &
            "&ocean_sponge_nml enable = .true., damp_source = 'band', "// &
            "target_source = 'ic', relax_uv = .true., relax_tracers = .true. /"//NL// &
            "&ocean_tracers_nml enable_pseudo_salt = .true. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function make_nml

   pure function bathymetry() result(b)
      real(wp) :: b(NX, NY)
      b = DEPTH
      b(LI, LJ) = 0.0_wp
   end function bathymetry

   pure function field3(base, di, dj, dk, ni) result(f)
      real(wp), intent(in) :: base, di, dj, dk
      integer, intent(in) :: ni
      real(wp) :: f(ni, NY, NZ)
      integer :: i, j, k
      do k = 1, NZ
         do j = 1, NY
            do i = 1, ni
               f(i, j, k) = base + di*real(i, wp) + dj*real(j, wp) + dk*real(k, wp)
            end do
         end do
      end do
   end function field3

   pure function thickness() result(h)
      real(wp) :: h(NX, NY, NZ)
      integer :: i, k
      do k = 1, NZ
         h(:, :, k) = DEPTH*real(k, wp)/6.0_wp
      end do
      h(3, 4, :) = h(3, 4, :)*[1.2_wp, 0.9_wp, 0.95_wp]
      h(3, 4, :) = h(3, 4, :)*DEPTH/sum(h(3, 4, :))
      ! A free-surface bulge growing with i: bt_eta = DEPTH*ETA_SLOPE*i.
      do i = 1, NX
         h(i, :, :) = h(i, :, :)*(1.0_wp + ETA_SLOPE*real(i, wp))
      end do
   end function thickness

   function v_field() result(v)
      real(wp) :: v(NX, NY + 1, NZ)
      integer :: i, j, k
      do k = 1, NZ
         do j = 1, NY + 1
            do i = 1, NX
               v(i, j, k) = 0.004_wp*real(j, wp) - 0.001_wp*real(i, wp) + 0.02_wp*real(k, wp)
            end do
         end do
      end do
   end function v_field

   function open_pending(c_handle) result(status)
      !! Pending handle with the island staged.
      type(c_ptr), intent(out) :: c_handle
      integer :: status
      character(len=:), allocatable :: nml
      real(wp) :: b(NX, NY)
      nml = make_nml()
      status = int(rdb_ocean_create_pending(nml, int(len(nml), c_int), c_handle))
      if (status /= OCEAN_STATUS_OK) return
      b = bathymetry()
      status = int(rdb_ocean_stage_bathymetry(c_handle, b, int(NX, c_int), int(NY, c_int), 1_c_int))
   end function open_pending

   function create_staged(c_handle, with_pseudo) result(status)
      type(c_ptr), intent(out) :: c_handle
      logical, intent(in) :: with_pseudo
      integer :: status
      real(wp) :: h(NX, NY, NZ), u(NX + 1, NY, NZ), v(NX, NY + 1, NZ)
      real(wp) :: tc(NX, NY, NZ), sc(NX, NY, NZ), pc(NX, NY, NZ)
      status = open_pending(c_handle)
      if (status /= OCEAN_STATUS_OK) return
      h = thickness()
      u = field3(0.01_wp, 0.01_wp, -0.002_wp, 0.01_wp, NX + 1)
      v = v_field()
      tc = field3(8.0_wp, 0.1_wp, 0.03_wp, 1.7_wp, NX)
      sc = field3(33.0_wp, 0.05_wp, 0.01_wp, 0.4_wp, NX)
      pc = field3(20.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, NX)
      status = int(rdb_ocean_stage_h(c_handle, h, int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_stage_u(c_handle, u, int(NX + 1, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_stage_v(c_handle, v, int(NX, c_int), int(NY + 1, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_stage_tracer(c_handle, "temperature", 11_c_int, tc, &
                                                   int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_stage_tracer(c_handle, "salinity", 8_c_int, sc, &
                                                   int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      if (with_pseudo) status = status + int(rdb_ocean_stage_tracer(c_handle, "pseudo_salt", 11_c_int, pc, &
                                                                    int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      if (status /= OCEAN_STATUS_OK) return
      status = int(rdb_ocean_create_finalize(c_handle))
   end function create_staged

   function create_then_set(c_handle, with_pseudo) result(status)
      !! The same initial condition through the post-create setters.  h
      !! first, tracers after (set_tracer forms hTr against the live h).
      type(c_ptr), intent(out) :: c_handle
      logical, intent(in) :: with_pseudo
      integer :: status
      real(wp) :: h(NX, NY, NZ), u(NX + 1, NY, NZ), v(NX, NY + 1, NZ)
      real(wp) :: tc(NX, NY, NZ), sc(NX, NY, NZ), pc(NX, NY, NZ)
      status = open_pending(c_handle)
      if (status /= OCEAN_STATUS_OK) return
      status = int(rdb_ocean_create_finalize(c_handle))
      if (status /= OCEAN_STATUS_OK) return
      h = thickness()
      u = field3(0.01_wp, 0.01_wp, -0.002_wp, 0.01_wp, NX + 1)
      v = v_field()
      tc = field3(8.0_wp, 0.1_wp, 0.03_wp, 1.7_wp, NX)
      sc = field3(33.0_wp, 0.05_wp, 0.01_wp, 0.4_wp, NX)
      pc = field3(20.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, NX)
      status = int(rdb_ocean_set_h(c_handle, h, int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_set_u(c_handle, u, int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_set_v(c_handle, v, int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_set_tracer(c_handle, "temperature", 11_c_int, tc, &
                                                 int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      status = status + int(rdb_ocean_set_tracer(c_handle, "salinity", 8_c_int, sc, &
                                                 int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
      if (with_pseudo) status = status + int(rdb_ocean_set_tracer(c_handle, "pseudo_salt", 11_c_int, pc, &
                                                                  int(NX, c_int), int(NY, c_int), int(NZ, c_int)))
   end function create_then_set

   subroutine take_snapshot(c_handle, s)
      type(c_ptr), intent(in) :: c_handle
      type(snap_t), intent(out) :: s
      type(ocean_handle_t), pointer :: h
      integer :: st, is, it, ip

      st = int(rdb_ocean_refresh_host(c_handle))
      st = handle_check(c_handle, h)
      is = h%state%multilayer%idx_salinity
      it = h%state%multilayer%idx_temperature
      ip = h%state%multilayer%idx_pseudo_salt
      s%h = h%state%multilayer%h_layer
      s%u = h%state%multilayer%u_face_x_layer
      s%v = h%state%multilayer%v_face_y_layer
      s%hT = h%state%multilayer%tracers(it)%hTr
      s%hS = h%state%multilayer%tracers(is)%hTr
      s%hP = h%state%multilayer%tracers(ip)%hTr
      s%rho = h%state%multilayer%rho_layer
      s%eta = h%state%dyn%bt_work%bt_eta
      s%refT = h%state%sponge%ref_tracer(:, :, :, it)
      s%refS = h%state%sponge%ref_tracer(:, :, :, is)
      s%refP = h%state%sponge%ref_tracer(:, :, :, ip)
      s%uref = h%state%sponge%u_ref
      s%vref = h%state%sponge%v_ref
      s%wet = h%state%multilayer%wet_mask
      s%wet_u = h%state%metrics%wet_u
      s%wet_v = h%state%metrics%wet_v
   end subroutine take_snapshot

   function max_diff_wet(a, b, wet) result(d)
      !! Largest |a - b| over the interior wet T-cells (3-D, cell-centred).
      real(wp), intent(in) :: a(:, :, :), b(:, :, :), wet(:, :)
      real(wp) :: d
      integer :: i, j
      d = 0.0_wp
      do j = NG + 1, NG + NY
         do i = NG + 1, NG + NX
            if (wet(i, j) <= 0.5_wp) cycle
            d = max(d, maxval(abs(a(i, j, :) - b(i, j, :))))
         end do
      end do
   end function max_diff_wet

   function max_diff_interior(a, b, ni, nj) result(d)
      !! Largest |a - b| over the interior physical block (faces: ni = NX+1 ...).
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      integer, intent(in) :: ni, nj
      real(wp) :: d
      d = maxval(abs(a(NG + 1:NG + ni, NG + 1:NG + nj, :) - b(NG + 1:NG + ni, NG + 1:NG + nj, :)))
   end function max_diff_interior

   function max_diff_faces(a, b, wet_f, ni, nj) result(d)
      !! Largest |a - b| over the open (wet_f > 0.5) interior faces.
      real(wp), intent(in) :: a(:, :, :), b(:, :, :), wet_f(:, :)
      integer, intent(in) :: ni, nj
      real(wp) :: d
      integer :: i, j
      d = 0.0_wp
      do j = NG + 1, NG + nj
         do i = NG + 1, NG + ni
            if (wet_f(i, j) <= 0.5_wp) cycle
            d = max(d, maxval(abs(a(i, j, :) - b(i, j, :))))
         end do
      end do
   end function max_diff_faces

   subroutine compare_all(error, a, b, tag, ref_too, eta_too)
      type(error_type), allocatable, intent(inout) :: error
      type(snap_t), intent(in) :: a, b
      character(len=*), intent(in) :: tag
      logical, intent(in) :: ref_too
      logical, intent(in) :: eta_too
      real(wp) :: d

      d = max_diff_wet(a%h, b%h, a%wet)
      call check(error, d == 0.0_wp, tag//": h differs on wet cells")
      if (allocated(error)) return
      d = max_diff_wet(a%hT, b%hT, a%wet)
      call check(error, d == 0.0_wp, tag//": hTr(temperature) differs on wet cells")
      if (allocated(error)) return
      d = max_diff_wet(a%hS, b%hS, a%wet)
      call check(error, d == 0.0_wp, tag//": hTr(salinity) differs on wet cells")
      if (allocated(error)) return
      d = max_diff_wet(a%hP, b%hP, a%wet)
      call check(error, d == 0.0_wp, tag//": hTr(pseudo_salt) differs on wet cells")
      if (allocated(error)) return
      d = max_diff_interior(a%u, b%u, NX + 1, NY)
      call check(error, d == 0.0_wp, tag//": u differs (land faces must be zero in both)")
      if (allocated(error)) return
      d = max_diff_interior(a%v, b%v, NX, NY + 1)
      call check(error, d == 0.0_wp, tag//": v differs (land faces must be zero in both)")
      if (allocated(error)) return
      d = max_diff_wet(a%rho, b%rho, a%wet)
      call check(error, d == 0.0_wp, tag//": rho_layer differs on wet cells")
      if (allocated(error)) return
      if (eta_too) then
         d = max_diff_wet(spread(a%eta, 3, 1), spread(b%eta, 3, 1), a%wet)
         call check(error, d == 0.0_wp, tag//": bt_eta differs on wet cells")
         if (allocated(error)) return
      end if
      if (ref_too) then
         d = max_diff_wet(a%refT, b%refT, a%wet)
         call check(error, d == 0.0_wp, tag//": sponge ref_tracer(T) differs on wet cells")
         if (allocated(error)) return
         d = max_diff_wet(a%refS, b%refS, a%wet)
         call check(error, d == 0.0_wp, tag//": sponge ref_tracer(S) differs on wet cells")
         if (allocated(error)) return
         d = max_diff_wet(a%refP, b%refP, a%wet)
         call check(error, d == 0.0_wp, tag//": sponge ref_tracer(pseudo_salt) differs on wet cells")
         if (allocated(error)) return
         ! Open faces only, and not the periodic duplicate seam face i = NX+1:
         ! the staged path snapshots BEFORE the land-face mask and the seam
         ! wrap (a staged u on a closed face, and the east copy of the seam
         ! face, stay as staged in u_ref); the rebaseline masks and wraps
         ! first, so u_ref is 0 on closed faces and the wrapped value at the
         ! seam face.
         d = max_diff_faces(a%uref, b%uref, a%wet_u, NX, NY)
         call check(error, d == 0.0_wp, tag//": sponge u_ref differs on open faces")
         if (allocated(error)) return
         d = max_diff_faces(a%vref, b%vref, a%wet_v, NX, NY + 1)
         call check(error, d == 0.0_wp, tag//": sponge v_ref differs on open faces")
      end if
   end subroutine compare_all

   subroutine test_matches_staged(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hnd
      type(snap_t) :: s0_staged, s10_staged, s0_re, s10_re, s0_stale
      integer :: st, i, j
      real(wp) :: worst

      ! (1) staged
      st = create_staged(hnd, .false.)
      call check(error, st, OCEAN_STATUS_OK, "staged create failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s0_staged)
      st = int(rdb_ocean_step(hnd, 10_c_int))
      call check(error, st, OCEAN_STATUS_OK, "staged step failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s10_staged)
      st = int(rdb_ocean_destroy(hnd))

      ! (2a) the stale state WITHOUT rebaseline: the sponge target and the
      ! pseudo-salt are still the built-in uniform ones (anti-vacuity).
      st = create_then_set(hnd, .false.)
      call check(error, st, OCEAN_STATUS_OK, "post-create sets failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s0_stale)
      call check(error, max_diff_wet(s0_stale%refT, s0_staged%refT, s0_staged%wet) > 1.0e-3_wp, &
                 "without rebaseline the sponge reference should still be the built-in IC")
      if (allocated(error)) return
      call check(error, max_diff_wet(s0_stale%hP, s0_staged%hP, s0_staged%wet) > 1.0e-3_wp, &
                 "without rebaseline pseudo-salt should still be the built-in S")
      if (allocated(error)) return

      ! (2b) rebaseline
      st = int(rdb_ocean_rebaseline(hnd))
      call check(error, st, OCEAN_STATUS_OK, "rebaseline failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s0_re)
      call compare_all(error, s0_staged, s0_re, "t=0", .true., .false.)
      if (allocated(error)) return
      ! The eta getter is current straight after the rebaseline: sum(h) - H.
      worst = 0.0_wp
      do j = NG + 1, NG + NY
         do i = NG + 1, NG + NX
            if (s0_re%wet(i, j) <= 0.5_wp) cycle
            worst = max(worst, abs(s0_re%eta(i, j) - DEPTH*ETA_SLOPE*real(i - NG, wp)))
         end do
      end do
      call check(error, worst < 1.0e-9_wp, "bt_eta is not sum(h) - H after the rebaseline")
      if (allocated(error)) return
      st = int(rdb_ocean_step(hnd, 10_c_int))
      call check(error, st, OCEAN_STATUS_OK, "rebaselined step failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s10_re)
      call compare_all(error, s10_staged, s10_re, "after 10 steps", .false., .true.)
      st = int(rdb_ocean_destroy(hnd))
   end subroutine test_matches_staged

   subroutine test_refused_after_step(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hnd
      integer :: st

      st = create_then_set(hnd, .false.)
      call check(error, st, OCEAN_STATUS_OK, "post-create sets failed")
      if (allocated(error)) return
      st = int(rdb_ocean_step(hnd, 1_c_int))
      st = int(rdb_ocean_rebaseline(hnd))
      call check(error, st, OCEAN_STATUS_ERR_SETUP, "rebaseline after step(1) must fail ERR_SETUP")
      st = int(rdb_ocean_destroy(hnd))
   end subroutine test_refused_after_step

   subroutine test_user_pseudo_salt(error)
      !! A pseudo-salt the caller wrote survives; it is NOT the S re-seed.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hnd
      type(snap_t) :: s_staged, s_re
      integer :: st

      st = create_staged(hnd, .true.)
      call check(error, st, OCEAN_STATUS_OK, "staged create (pseudo staged) failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s_staged)
      st = int(rdb_ocean_destroy(hnd))

      st = create_then_set(hnd, .true.)
      call check(error, st, OCEAN_STATUS_OK, "post-create sets failed")
      if (allocated(error)) return
      st = int(rdb_ocean_rebaseline(hnd))
      call check(error, st, OCEAN_STATUS_OK, "rebaseline failed")
      if (allocated(error)) return
      call take_snapshot(hnd, s_re)
      call compare_all(error, s_staged, s_re, "user pseudo-salt", .true., .false.)
      if (allocated(error)) return
      call check(error, max_diff_wet(s_re%hP, s_re%hS, s_re%wet) > 1.0e-3_wp, &
                 "pseudo-salt was re-seeded from S although the caller set it")
      st = int(rdb_ocean_destroy(hnd))
   end subroutine test_user_pseudo_salt

end module test_ocean_rebaseline
