!! `VCOORD_Z_FIXED` on a STRETCHED nominal profile
!! (`&vcoord_nml z_fixed_profile = "list" | "tanh"`).
!!
!! Until this knob `z_fixed` only knew UNIFORM nominal layers
!! (`h_nominal = max_depth/nz`, i.e. 130 m each at 50 levels over a 6500 m
!! ocean), which a realistic global grid cannot use.  A profile replaces the
!! one spacing with a per-layer table, and every consumer of "where are the
!! nominal interfaces" must read it: the target builder (bed partial cell,
!! partial top cell under a draft, the fillers at both ends), and through
!! the target the closed-face mask and `k_top`.  This suite pins each of
!! them on hand-computable columns, then runs the model:
!!
!!   * `tanh_profile_shape` / `list_profile_*` — the profile builder
!!     (`rdb_vcoord :: z_fixed_nominal_dz`): exact surface thickness,
!!     monotone, sums to `max_depth`; a list of the wrong length refused.
!!   * `stretched_flat_bed` / `stretched_partial_bed` /
!!     `stretched_partial_top_under_draft` / `stretched_top_sliver_merges`
!!     — the target on a 6-layer 2/4/8/16/30/40 m column, every thickness
!!     derived by hand in the test.
!!   * `stretched_closed_faces_staircase` — three columns stepping down a
!!     stretched staircase; the mask closes exactly the filler layers.
!!   * `list_of_equal_dz_matches_uniform` — a list of `nz` equal entries is
!!     the uniform coordinate (to round-off: the two tables are built by
!!     different arithmetic, so this is NOT a bitwise claim).
!!   * `flat_lid_rest_is_exact_zero` — the full model: a resting stratified
!!     cavity under a FLAT ice shelf on a tanh profile keeps every face
!!     velocity at exactly 0.0 (the uniform-profile gate
!!     `cavity_flat_lid_rest_zfixed.nml`, on a stretched stack).
!!   * `refuses_*` — `validate_config` fails loud on a list of the wrong
!!     length and on a profile set on another coordinate.
!!
!! Knob OFF (the default `"uniform"`) takes the unchanged uniform branch of
!! the same kernel; that byte-identity is carried by every pre-existing
!! `z_fixed` test and golden, not re-asserted here.
module test_ocean_zfixed_stretched
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp, H_VANISHED
   use rdb_vcoord, only: z_fixed_nominal_dz, ZFIXED_PROFILE_LIST, ZFIXED_PROFILE_TANH, &
                         ZFIXED_DZ_OK, ZFIXED_DZ_ERR_COUNT, ZFIXED_DZ_ERR_TOO_DEEP
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks, &
                               ocean_vcoord_k_top_from_target
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_h_layer_ptr, &
                            rdb_ocean_get_wet_t_ptr, rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_zfixed_stretched_tests

   integer, parameter :: NZ = 6
      !! Layers of the hand-computable column.
   real(wp), parameter :: DZ_SURF_FIRST(NZ) = [2.0_wp, 4.0_wp, 8.0_wp, 16.0_wp, 30.0_wp, 40.0_wp]
      !! Nominal thicknesses, surface first (total 100 m).
   real(wp), parameter :: H_MIN = 1.0e-4_wp
      !! Inert-filler thickness (strictly below `H_VANISHED`).
   real(wp), parameter :: TOL = 1.0e-12_wp
      !! Absolute tolerance on a hand-derived thickness (m).

contains

   subroutine collect_ocean_zfixed_stretched_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("tanh_profile_shape", test_tanh_shape), &
                  new_unittest("list_profile_count_refused", test_list_count), &
                  new_unittest("tanh_profile_too_deep_refused", test_tanh_too_deep), &
                  new_unittest("stretched_flat_bed", test_flat_bed), &
                  new_unittest("stretched_partial_bed", test_partial_bed), &
                  new_unittest("stretched_partial_top_under_draft", test_partial_top), &
                  new_unittest("stretched_top_sliver_merges", test_top_sliver), &
                  new_unittest("stretched_closed_faces_staircase", test_closed_faces), &
                  new_unittest("list_of_equal_dz_matches_uniform", test_list_matches_uniform), &
                  new_unittest("flat_lid_rest_is_exact_zero", test_flat_lid_rest), &
                  new_unittest("refuses_list_of_wrong_length", test_refuses_list_length), &
                  new_unittest("refuses_profile_off_z_fixed", test_refuses_off_z_fixed) &
                  ]
   end subroutine collect_ocean_zfixed_stretched_tests

   pure subroutine tables(dz_sf, n, zi, dz)
      !! Bottom-up tables from a surface-first list — the same flip and
      !! accumulation `ocean_vcoord_set_z_fixed_profile` performs.
      integer, intent(in) :: n
      real(wp), intent(in) :: dz_sf(n)
      real(wp), intent(out) :: zi(0:n), dz(n)
      integer :: k
      do k = 1, n
         dz(k) = dz_sf(n - k + 1)
      end do
      zi(n) = 0.0_wp
      do k = n, 1, -1
         zi(k - 1) = zi(k) + dz(k)
      end do
   end subroutine tables

   subroutine column(tgt, total_h, z_top)
      !! One-column (plus ring) stretched target at `eta = 0`.
      real(wp), intent(out) :: tgt(3, 3, NZ)
      real(wp), intent(in) :: total_h, z_top
      real(wp) :: th(3, 3), eta(3, 3), zt(3, 3), zi(0:NZ), dz(NZ)
      th = total_h
      eta = 0.0_wp
      zt = z_top
      call tables(DZ_SURF_FIRST, NZ, zi, dz)
      call ocean_vcoord_z_fixed_target(tgt, th, eta, zt, 3, 3, NZ, 0.0_wp, .true., &
                                       zi, dz, H_MIN)
   end subroutine column

   subroutine expect(error, got, want, what)
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: got(:), want(:)
      character(len=*), intent(in) :: what
      integer :: k
      character(len=160) :: msg
      do k = 1, size(want)
         write (msg, "(a,a,i0,a,es22.15,a,es22.15)") what, ": k=", k, " got ", got(k), &
            " want ", want(k)
         call check(error, abs(got(k) - want(k)) <= TOL, trim(msg))
         if (allocated(error)) return
      end do
   end subroutine expect

   subroutine test_tanh_shape(error)
      !! 50 levels over 6500 m from a 2 m surface layer: exact surface
      !! thickness, monotone, closes on max_depth, ~200+ m at depth.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 50
      real(wp) :: dz(N), unused(1)
      integer :: ierr, k
      unused = -1.0_wp
      call z_fixed_nominal_dz(ZFIXED_PROFILE_TANH, N, 6500.0_wp, unused, 2.0_wp, &
                              0.5_wp, 0.25_wp, dz, ierr)
      call check(error, ierr == ZFIXED_DZ_OK, "tanh profile did not build")
      if (allocated(error)) return
      call check(error, dz(1) == 2.0_wp, "surface layer is not exactly dz_top")
      if (allocated(error)) return
      do k = 2, N
         call check(error, dz(k) >= dz(k - 1), "tanh profile is not monotone")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(dz) - 6500.0_wp) <= 1.0e-9_wp*6500.0_wp, &
                 "tanh profile does not sum to max_depth")
      if (allocated(error)) return
      call check(error, dz(N) > 200.0_wp .and. dz(N) < 400.0_wp, &
                 "tanh bed layer is not O(200+ m)")
   end subroutine test_tanh_shape

   subroutine test_list_count(error)
      !! A list shorter than nz, or with a gap, is refused; the exact list
      !! is taken verbatim.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz(NZ), lst(10)
      integer :: ierr
      lst = -1.0_wp
      lst(1:NZ) = DZ_SURF_FIRST
      call z_fixed_nominal_dz(ZFIXED_PROFILE_LIST, NZ, 0.0_wp, lst, 0.0_wp, 0.0_wp, &
                              0.0_wp, dz, ierr)
      call check(error, ierr == ZFIXED_DZ_OK .and. all(dz == DZ_SURF_FIRST), &
                 "exact list not taken verbatim")
      if (allocated(error)) return
      lst(NZ) = -1.0_wp
      call z_fixed_nominal_dz(ZFIXED_PROFILE_LIST, NZ, 0.0_wp, lst, 0.0_wp, 0.0_wp, &
                              0.0_wp, dz, ierr)
      call check(error, ierr == ZFIXED_DZ_ERR_COUNT, "short list accepted")
      if (allocated(error)) return
      lst(1:NZ) = DZ_SURF_FIRST
      lst(3) = -1.0_wp
      call z_fixed_nominal_dz(ZFIXED_PROFILE_LIST, NZ, 0.0_wp, lst, 0.0_wp, 0.0_wp, &
                              0.0_wp, dz, ierr)
      call check(error, ierr == ZFIXED_DZ_ERR_COUNT, "list with a gap accepted")
   end subroutine test_list_count

   subroutine test_tanh_too_deep(error)
      !! nz*dz_top >= max_depth leaves nothing to stretch into.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz(10), unused(1)
      integer :: ierr
      unused = -1.0_wp
      call z_fixed_nominal_dz(ZFIXED_PROFILE_TANH, 10, 100.0_wp, unused, 10.0_wp, &
                              0.5_wp, 0.25_wp, dz, ierr)
      call check(error, ierr == ZFIXED_DZ_ERR_TOO_DEEP, "tanh with no room accepted")
   end subroutine test_tanh_too_deep

   subroutine test_flat_bed(error)
      !! A column exactly as deep as the profile: every layer at its
      !! nominal thickness, bottom-up.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(3, 3, NZ)
      call column(tgt, 100.0_wp, 0.0_wp)
      call expect(error, tgt(2, 2, :), [40.0_wp, 30.0_wp, 16.0_wp, 8.0_wp, 4.0_wp, 2.0_wp], &
                  "flat bed")
   end subroutine test_flat_bed

   subroutine test_partial_bed(error)
      !! Bed at 40 m: the top four layers (2+4+8+16 = 30 m) are nominal,
      !! the 30 m layer is the partial BOTTOM cell holding the remaining
      !! 10 m less the one bed filler's debt, and the 40 m layer is a filler.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(3, 3, NZ)
      call column(tgt, 40.0_wp, 0.0_wp)
      call expect(error, tgt(2, 2, :), &
                  [H_MIN, 10.0_wp - H_MIN, 16.0_wp, 8.0_wp, 4.0_wp, 2.0_wp], "partial bed")
      if (allocated(error)) return
      call check(error, abs(sum(tgt(2, 2, :)) - 40.0_wp) <= TOL, "partial bed: column not closed")
   end subroutine test_partial_bed

   subroutine test_partial_top(error)
      !! Flat bed at 100 m under a 5 m draft (water column 95 m).  The 2 m
      !! surface layer lies inside the ice (filler); the 4 m layer [2, 6]
      !! STRADDLES the ice base and keeps 6 - 5 = 1 m less the filler's
      !! debt.  1 m clears ITS OWN threshold 0.1*4 = 0.4 m — the uniform
      !! spacing (100/6 m) would have put the threshold at 1.67 m and
      !! merged it, so this also pins the per-layer partial-top rule.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(3, 3, NZ)
      integer :: k_top(3, 3), k_top_u(4, 3), k_top_v(3, 4)
      call column(tgt, 95.0_wp, 5.0_wp)
      call expect(error, tgt(2, 2, :), &
                  [40.0_wp, 30.0_wp, 16.0_wp, 8.0_wp, 1.0_wp - H_MIN, H_MIN], &
                  "partial top under a draft")
      if (allocated(error)) return
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, 3, 3, NZ, H_VANISHED)
      call check(error, k_top(2, 2) == NZ - 1, "k_top is not the partial top cell")
   end subroutine test_partial_top

   subroutine test_top_sliver(error)
      !! Draft 5.8 m: the 4 m layer would keep only 0.2 m < 0.4 m, so the
      !! sliver merges into the 8 m layer below (which then holds
      !! 14 - 5.8 = 8.2 m less two fillers' debt) and both upper layers are
      !! fillers.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tgt(3, 3, NZ)
      call column(tgt, 100.0_wp - 5.8_wp, 5.8_wp)
      call expect(error, tgt(2, 2, :), &
                  [40.0_wp, 30.0_wp, 16.0_wp, (14.0_wp - 5.8_wp) - 2.0_wp*H_MIN, H_MIN, H_MIN], &
                  "top sliver")
   end subroutine test_top_sliver

   subroutine test_closed_faces(error)
      !! Three columns stepping down the stretched staircase — beds at 100,
      !! 40 and 20 m.  Column B vanishes the 40 m bed layer (k = 1);
      !! column C also vanishes the 30 m layer (k = 2, nominal range
      !! [30, 60] entirely below its 20 m bed) and keeps 6 m of the 16 m
      !! one.  The A|B face must close k = 1 only; the B|C face k = 1, 2.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXC = 5, NYC = 3
      real(wp) :: tgt(NXC, NYC, NZ), th(NXC, NYC), eta(NXC, NYC), zt(NXC, NYC)
      real(wp) :: ou(NXC + 1, NYC, NZ), ov(NXC, NYC + 1, NZ), zi(0:NZ), dz(NZ)
      real(wp) :: want_ab(NZ), want_bc(NZ)
      integer :: k
      character(len=64) :: msg
      th(1:2, :) = 100.0_wp
      th(3, :) = 40.0_wp
      th(4:5, :) = 20.0_wp
      eta = 0.0_wp
      zt = 0.0_wp
      call tables(DZ_SURF_FIRST, NZ, zi, dz)
      call ocean_vcoord_z_fixed_target(tgt, th, eta, zt, NXC, NYC, NZ, 0.0_wp, .true., &
                                       zi, dz, H_MIN)
      call check(error, abs(tgt(4, 2, 3) - (6.0_wp - 2.0_wp*H_MIN)) <= TOL, &
                 "column C partial bed cell is not 6 m less two fillers")
      if (allocated(error)) return
      call ocean_vcoord_closed_face_masks(ou, ov, tgt, NXC, NYC, NZ, H_VANISHED)
      want_ab = [0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      want_bc = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      do k = 1, NZ
         write (msg, "(a,i0)") "closed-face mask wrong at k=", k
         ! u-face I is the west face of cell I: A|B is I = 3, B|C is I = 4.
         call check(error, ou(3, 2, k) == want_ab(k) .and. ou(4, 2, k) == want_bc(k) .and. &
                    ou(2, 2, k) == 1.0_wp .and. ou(5, 2, k) == want_bc(k), trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_closed_faces

   subroutine test_list_matches_uniform(error)
      !! `nz` equal entries of `H/nz` reproduce the uniform coordinate on a
      !! column with BOTH a partial bed and a partial top cell, to
      !! round-off.  (The two paths build the interface depths by
      !! different arithmetic — a product vs a running sum — so bitwise
      !! agreement is neither expected nor claimed.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 12
      real(wp), parameter :: HREF = 720.0_wp
      real(wp) :: t_uni(3, 3, N), t_lst(3, 3, N), th(3, 3), eta(3, 3), zt(3, 3)
      real(wp) :: zi(0:N), dz(N), sf(N)
      th = 437.0_wp - 123.0_wp
      eta = 0.25_wp
      zt = 123.0_wp
      sf = HREF/real(N, wp)
      call tables(sf, N, zi, dz)
      call ocean_vcoord_z_fixed_target_uniform(t_uni, th, eta, zt, 3, 3, N, HREF/real(N, wp), H_MIN)
      call ocean_vcoord_z_fixed_target(t_lst, th, eta, zt, 3, 3, N, 0.0_wp, .true., zi, dz, H_MIN)
      call check(error, maxval(abs(t_uni - t_lst)) <= 1.0e-10_wp, &
                 "a list of equal thicknesses is not the uniform coordinate")
   end subroutine test_list_matches_uniform

   function nml_flat_lid(extra) result(nml)
      !! `cavity_flat_lid_rest_zfixed.nml`, reduced (24 x 6), on a tanh
      !! profile from a 10 m surface layer over 720 m / 15 levels.
      character(len=*), intent(in) :: extra
         !! Appended `&vcoord_nml` keys.
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 24, ny = 6, nghost = 2, dx = 2000.0, dy = 2000.0 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 600.0, cfl_interval = 1 /"//nl// &
            "&physics_nml coriolis_f = -1.409e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 15 /"//nl// &
            "&vcoord_nml "//extra//" /"//nl// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 720.0, wind_config = 'constant', "// &
            "taux_magnitude = 0.0, coriolis_beta = 0.0, coriolis_y_ref = 0.0 /"//nl// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', draft_depth = 300.0, "// &
            "h_min_cavity = 96.0 /"//nl// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true., maxvel = 2.0 /"//nl// &
            "&ocean_eos_nml eos = 'linear' /"//nl// &
            "&ocean_ic_nml alpha_T = 3.8356948e-2, beta_S = 8.0587609e-1, T_ref = -1.0, "// &
            "S_ref = 34.2, rho_0 = 1027.51 /"//nl// &
            "&ocean_zinit_nml enable = .true., source = 'linear', lin_t_ref = -1.9, "// &
            "lin_dt_dz = 0.0, lin_s_ref = 33.8, lin_ds_dz = -1.0416667e-3 /"//nl// &
            "&tracer_nml initial_temperature = -1.9, initial_salinity = 34.175 /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny' /"//nl// &
            "&ocean_hvisc_nml nu_h = 0.0, nu_4 = 0.0 /"//nl// &
            "&ocean_bdrag_nml form = 'quadratic', cd = 0.0 /"//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function nml_flat_lid

   subroutine test_flat_lid_rest(error)
      !! Resting, stratified, flat ice shelf, stretched `z_fixed` stack with
      !! closed faces: every column vanishes the SAME layers and cuts its
      !! partial top cell at the SAME depth, so every interface offset is
      !! identically zero and the answer must be BIT-ZERO — the bar the
      !! uniform-profile twin meets.  Non-vacuity: the stack really is
      !! stretched (two live layers of different thickness) and really has
      !! fillers under the ice.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :), u(:, :, :), v(:, :, :), wet(:, :)
      real(wp) :: umax, vmax, hmin_live, hmax_live
      integer :: n_fill, i, j, k
      character(len=200) :: msg

      nml = nml_flat_lid("vcoord_type = 'z_fixed', zfixed_closed_faces = .true., "// &
                         "check_vanished_content = .true., z_fixed_profile = 'tanh', "// &
                         "z_fixed_dz_top = 10.0, z_fixed_tanh_center = 0.5, "// &
                         "z_fixed_tanh_width = 0.25")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "stretched flat-lid case must build")
      if (allocated(error)) return

      body: block
         status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
         status = rdb_ocean_step(handle, 60_c_int)
         call check(error, status == OCEAN_STATUS_OK, "stretched flat-lid case must step")
         if (allocated(error)) exit body
         status = rdb_ocean_refresh_host(handle)
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h, [nx, ny, nz])
         status = rdb_ocean_get_wet_t_ptr(handle, ptr, nx, ny, gen)
         call c_f_pointer(ptr, wet, [nx, ny])
         status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, u, [nx, ny, nz])
         umax = maxval(abs(u))
         status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, v, [nx, ny, nz])
         vmax = maxval(abs(v))

         n_fill = 0
         hmin_live = huge(1.0_wp)
         hmax_live = 0.0_wp
         do k = 1, nz_p
            do j = ng + 1, ng + ny_p
               do i = ng + 1, ng + nx_p
                  if (wet(i, j) <= 0.5_wp) cycle
                  if (h(i, j, k) <= H_VANISHED) then
                     n_fill = n_fill + 1
                  else if (k > 1) then
                     ! k = 1 is the bed cell, which absorbs the residual;
                     ! the interior live layers carry nominal thicknesses.
                     hmin_live = min(hmin_live, h(i, j, k))
                     hmax_live = max(hmax_live, h(i, j, k))
                  end if
               end do
            end do
         end do
         call check(error, n_fill > 0, "no filler under the ice — the test is vacuous")
         if (allocated(error)) exit body
         call check(error, hmax_live > 1.5_wp*hmin_live, &
                    "live layers are not stretched — the test is vacuous")
         if (allocated(error)) exit body
         write (msg, "(a,es10.3,a,es10.3)") "flat lid at rest moved: max|u| = ", umax, &
            ", max|v| = ", vmax
         call check(error, umax == 0.0_wp .and. vmax == 0.0_wp, trim(msg))
      end block body
      status = rdb_ocean_destroy(handle)
   end subroutine test_flat_lid_rest

   subroutine test_refuses_list_length(error)
      !! `list` with 3 entries for 15 layers must not build.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml
      type(c_ptr) :: handle
      integer(c_int) :: status
      nml = nml_flat_lid("vcoord_type = 'z_fixed', zfixed_closed_faces = .true., "// &
                         "z_fixed_profile = 'list', z_fixed_dz = 10.0, 20.0, 30.0")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status /= OCEAN_STATUS_OK, "a z_fixed_dz list of the wrong length built")
      if (status == OCEAN_STATUS_OK) status = rdb_ocean_destroy(handle)
   end subroutine test_refuses_list_length

   subroutine test_refuses_off_z_fixed(error)
      !! A profile on a coordinate that never reads it must not build.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml
      type(c_ptr) :: handle
      integer(c_int) :: status
      nml = nml_flat_lid("vcoord_type = 'sigma', z_fixed_profile = 'tanh'")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status /= OCEAN_STATUS_OK, "z_fixed_profile on sigma built")
      if (status == OCEAN_STATUS_OK) status = rdb_ocean_destroy(handle)
   end subroutine test_refuses_off_z_fixed

end module test_ocean_zfixed_stretched
