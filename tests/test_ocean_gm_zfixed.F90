!! Gent-McWilliams thickness diffusion on `z_fixed` with partial-step
!! CLOSED faces (`&vcoord_nml zfixed_closed_faces`).  Every case runs the
!! production device kernels — `ocean_slopes_compute` +
!! `gm_compute_transports` (+ the GM operator `continuity_gm_apply` and
!! the I1′ enforcement point for the multi-step case) — over a
!! hand-specified staircase whose layer thicknesses AND closed-face masks
!! come from the SAME production builders a configured run uses
!! (`ocean_vcoord_z_fixed_target_uniform` → `ocean_vcoord_closed_face_masks`),
!! so "filler" and "closed" have their one production definition.
!!
!! Bottom-up (k = 1 bed, k = nz top).  Linear EOS, T-only stratification.
!!
!!  1. `rest_flat_staircase_zero_transport` — z-level (layer-uniform)
!!     stratification over a staircase bed with partial bottom cells and
!!     bed fillers: the slope is 0 at every interface and GM moves NOTHING
!!     (to round-off: |S| <= 1e-14, |uhD| <= 1e-12·KhTh·dx, gm_src <=
!!     1e-20 — the T = hTr/h recovery is exact only to an ulp on a partial
!!     cell carrying the fillers' h_min debt).  Fails on the pre-port
!!     slopes, whose interface-tilt term differences the heights ABOVE THE
!!     LOCAL BED and so reads the bathymetry step as an isopycnal slope.
!!  2. `tilted_staircase_open_column` — tilted isopycnals over the same
!!     staircase, NSTEP GM-operator steps: a closed face-layer carries
!!     EXACTLY zero GM flux, every face's column-integrated transport is 0
!!     to round-off, no field is non-finite, the fillers keep their
!!     thickness exactly, I1′ holds after every step, and mass + T content
!!     are conserved.  Fails on the pre-port GM (flux through closed
!!     faces / into fillers via the full-column closure).
!!  3. `ice_top_closure_in_open_column` — fillers ABOVE the open column
!!     (an ice draft, `z_top > 0`): the non-divergence closure must land in
!!     the topmost OPEN layer, never in the `k = nz` filler.  The kernel
!!     property only — GM under a cavity is still refused at configure
!!     (the slopes slot is), this pins the closure rule it will need.
module test_ocean_gm_zfixed
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports
   use rdb_continuity, only: continuity_t, continuity_gm_apply, TR_MODE_ADVECT
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_gm_zfixed_tests

   integer, parameter :: NG = 2
   integer, parameter :: NXP = 8, NYP = 4, NZ = 8
   real(wp), parameter :: DX = 20000.0_wp
   real(wp), parameter :: H_NOM = 50.0_wp
      !! Nominal z_fixed spacing (m): NZ*H_NOM = 400 m full depth.
   real(wp), parameter :: H_MIN = 1.0e-4_wp
      !! `zstar_h_min` at its default: the inert filler thickness.
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: KHTH = 1000.0_wp
   real(wp), parameter :: GZ = 0.02_wp
      !! Vertical T gradient (K/m), warm above: stable.

contains

   subroutine collect_ocean_gm_zfixed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("rest_flat_staircase_zero_transport", test_rest_flat), &
                  new_unittest("tilted_staircase_open_column", test_tilted_open_column), &
                  new_unittest("ice_top_closure_in_open_column", test_ice_top_closure) &
                  ]
   end subroutine collect_ocean_gm_zfixed_tests

   ! ------------------------------------------------------------------
   ! Setup helpers
   ! ------------------------------------------------------------------

   subroutine build_case(grid, metrics, ms, sl, gm, eos, tgt, draft, gx, gy)
      !! Grid + closed-face metrics + a z_fixed staircase state.
      !!
      !! Bed depth steps down eastward AND northward (30 m per column, 20 m
      !! per row — off the 50 m nominal spacing, so the bottom live layer
      !! is a PARTIAL cell and the columns carry 0..3 bed fillers); ghost
      !! columns replicate the nearest physical column.  `draft` (m, >= 0)
      !! puts an ice base over the WHOLE domain for case 3.  The thickness
      !! field is the production `z_fixed` target at eta = 0 and the masks
      !! are built from that same target by the production builder.
      !!
      !! T is a function of NOMINAL layer depth plus a horizontal tilt
      !! `gx*x + gy*y`: with `gx = gy = 0` every column carries the same T
      !! in the same layer — a z-level (flat) stratification.  Fillers are
      !! then given their donor's concentration (I1′) by the host twin of
      !! the enforcement point.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(eos_t), intent(out) :: eos
      real(wp), allocatable, intent(out) :: tgt(:, :, :)
      real(wp), intent(in) :: draft, gx, gy
      real(wp), allocatable :: tot_h(:, :), eta0(:, :), z_top(:, :)
      integer :: i, j, k, ni, nj, ip, jp
      real(wp) :: zc, t_val

      call grid%init(NXP, NYP, NG, DX, DX)
      ni = grid%nx_total
      nj = grid%ny_total

      ! `nz_closed` grows open_u/open_v BEFORE the device map (see the
      ! helper's docstring).
      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)

      allocate (tot_h(ni, nj), eta0(ni, nj), z_top(ni, nj))
      allocate (tgt(ni, nj, NZ), source=0.0_wp)
      eta0 = 0.0_wp
      z_top = draft
      do j = 1, nj
         jp = min(max(j - NG, 1), NYP)
         do i = 1, ni
            ip = min(max(i - NG, 1), NXP)
            ! The target is laid in "depth below the column top"; the
            ! column (bed to ice base) holds `bed depth - draft`.
            tot_h(i, j) = (real(NZ, wp)*H_NOM - 30.0_wp*real(ip - 1, wp) &
                           - 20.0_wp*real(jp - 1, wp)) - draft
         end do
      end do
      call ocean_vcoord_z_fixed_target_uniform(tgt, tot_h, eta0, z_top, &
                                               ni, nj, NZ, H_NOM, H_MIN)
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                          tgt, ni, nj, NZ, H_VANISHED)
      metrics%use_closed_faces = .true.
      ! The masks are device-present (grown before the map), so the builder
      ! wrote the device copy; pull it back for the host-side assertions.
      !$omp target update from(metrics%open_u, metrics%open_v)

      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = tgt
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         zc = (real(k, wp) - 0.5_wp)*H_NOM          ! nominal height above 400 m
         do j = 1, nj
            jp = min(max(j - NG, 1), NYP)
            do i = 1, ni
               ip = min(max(i - NG, 1), NXP)
               t_val = T0 + GZ*zc + gx*real(ip, wp)*DX + gy*real(jp, wp)*DX
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*tgt(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*tgt(i, j, k)
            end do
         end do
      end do
      ! Hand each filler its donor's T BEFORE the I1' sweep, so the sweep's
      ! pooling mixes equal concentrations and does not shift the donor
      ! (a filler seeded at its own nominal-depth T would leak ~1e-4/h of
      ! that difference into the bottom live layer — a real along-layer
      ! gradient between columns with and without fillers).
      do j = 1, nj
         do i = 1, ni
            t_val = -1.0_wp
            do k = NZ, 1, -1
               if (tgt(i, j, k) > H_VANISHED) then
                  t_val = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/tgt(i, j, k)
               else if (t_val >= 0.0_wp) then
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*tgt(i, j, k)
               end if
            end do
         end do
      end do
      call ms%enforce_vanished_content_host(ni, nj)

      call sl%init(grid, nz_ml=NZ)
      ! Bed datum: the staircase bed depth (the column holds `bed - draft`).
      call sl%set_bathymetry(tot_h + draft)
      sl%enable = .true.
      sl%rho0 = RHO0
      ! No vert-fill smoothing: the smoothing couples each live layer to
      ! the one below it, and columns of different depth then disagree in
      ! the fourth digit of a z-level T — a real (tiny) slope, which would
      ! blur the exact-zero assertion of case 1.  Production keeps 1e-6.
      sl%kd_smooth = 0.0_wp
      sl%min_dz_for_n2 = 1.0_wp

      call gm%init(grid, nz_ml=NZ)
      gm%enable = .true.
      gm%khth = KHTH
      gm%khth_slope_max = 0.01_wp
      gm%rho0 = RHO0

      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine build_case

   subroutine map_in(ms, sl, gm, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(continuity_t), intent(inout), optional :: ct
      !$omp target enter data map(to: ms)
      call ms%enter_data()
      !$omp target enter data map(to: sl)
      call sl%enter_data()
      !$omp target enter data map(to: gm)
      call gm%enter_data()
      if (present(ct)) then
         !$omp target enter data map(to: ct)
         call ct%enter_data()
      end if
   end subroutine map_in

   subroutine pull_back(ms, sl, gm)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      !$omp target update from(gm%uhD, gm%vhD, gm%gm_src)
      !$omp target update from(sl%slope_x, sl%slope_y, sl%n2_u, sl%n2_v)
      !$omp target update from(ms%h_layer)
      !$omp target update from(ms%tracers(ms%idx_temperature)%hTr)
      !$omp target update from(ms%tracers(ms%idx_salinity)%hTr)
   end subroutine pull_back

   subroutine map_out(ms, sl, gm, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(continuity_t), intent(inout), optional :: ct
      if (present(ct)) then
         call ct%exit_data()
         !$omp target exit data map(delete: ct)
      end if
      call gm%exit_data()
      !$omp target exit data map(delete: gm)
      call sl%exit_data()
      !$omp target exit data map(delete: sl)
      call ms%exit_data()
      !$omp target exit data map(delete: ms)
   end subroutine map_out

   pure logical function all_finite_3d(a)
      real(wp), intent(in) :: a(:, :, :)
      all_finite_3d = all(ieee_is_finite(a))
   end function all_finite_3d

   pure logical function all_finite_2d(a)
      real(wp), intent(in) :: a(:, :)
      all_finite_2d = all(ieee_is_finite(a))
   end function all_finite_2d

   subroutine closed_face_census(metrics, ms, gm, ni, nj, n_closed, worst_closed, &
                                 worst_colsum, max_flux)
      !! Host census over EVERY face of the array: the largest |GM flux| on
      !! a face-layer that is closed (`open == 0`) or touches a filler on
      !! either side, the largest |Sum_k flux| of any face column, and the
      !! largest |flux| overall (the scale for the round-off assertion).
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_gm_t), intent(in) :: gm
      integer, intent(in) :: ni, nj
      integer, intent(out) :: n_closed
      real(wp), intent(out) :: worst_closed, worst_colsum, max_flux
      integer :: i, j, k
      real(wp) :: csum
      logical :: shut

      n_closed = 0
      worst_closed = 0.0_wp
      worst_colsum = 0.0_wp
      max_flux = max(maxval(abs(gm%uhD)), maxval(abs(gm%vhD)))
      do j = 1, nj
         do i = 2, ni
            csum = 0.0_wp
            do k = 1, NZ
               csum = csum + gm%uhD(i, j, k)
               shut = metrics%open_u(i, j, k) < 0.5_wp .or. &
                      ms%h_layer(i - 1, j, k) <= H_VANISHED .or. &
                      ms%h_layer(i, j, k) <= H_VANISHED
               if (shut) then
                  n_closed = n_closed + 1
                  worst_closed = max(worst_closed, abs(gm%uhD(i, j, k)))
               end if
            end do
            worst_colsum = max(worst_colsum, abs(csum))
         end do
      end do
      do j = 2, nj
         do i = 1, ni
            csum = 0.0_wp
            do k = 1, NZ
               csum = csum + gm%vhD(i, j, k)
               shut = metrics%open_v(i, j, k) < 0.5_wp .or. &
                      ms%h_layer(i, j - 1, k) <= H_VANISHED .or. &
                      ms%h_layer(i, j, k) <= H_VANISHED
               if (shut) then
                  n_closed = n_closed + 1
                  worst_closed = max(worst_closed, abs(gm%vhD(i, j, k)))
               end if
            end do
            worst_colsum = max(worst_colsum, abs(csum))
         end do
      end do
   end subroutine closed_face_census

   real(wp) function content_sum(ms, idx) result(s)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx
      s = sum(ms%tracers(idx)%hTr)
   end function content_sum

   ! ------------------------------------------------------------------
   ! Case 1: resting flat stratification — GM identically zero
   ! ------------------------------------------------------------------
   subroutine test_rest_flat(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :)
      integer :: n_fill
      real(wp) :: smax, fmax
      character(len=96) :: msg

      call build_case(grid, metrics, ms, sl, gm, eos, tgt, 0.0_wp, 0.0_wp, 0.0_wp)
      checks: block
         ! Non-vacuous: the staircase must actually carry fillers and
         ! partial cells, and the masks closed faces.
         n_fill = count(tgt <= H_VANISHED)
         call check(error, n_fill > 0, "the staircase must carry bed fillers")
         if (allocated(error)) exit checks
         call check(error, count(metrics%open_u < 0.5_wp) > 0 .and. &
                    count(metrics%open_v < 0.5_wp) > 0, &
                    "the staircase must close u- AND v-face layers")
         if (allocated(error)) exit checks

         call map_in(ms, sl, gm)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call pull_back(ms, sl, gm)
         call map_out(ms, sl, gm)

         smax = max(maxval(abs(sl%slope_x)), maxval(abs(sl%slope_y)))
         fmax = max(maxval(abs(gm%uhD)), maxval(abs(gm%vhD)))
         write (msg, '(a,es10.3,a,es10.3)') "max|S| = ", smax, ", max|uhD,vhD| = ", fmax
         ! Zero to ROUND-OFF, not bitwise: the slope pass recovers T as
         ! hTr/h, which returns the seeded T only to an ulp when h is a
         ! partial cell carrying the fillers' `h_min` debt, so two columns'
         ! T at one level can differ by ~1e-16 relative.  The bounds are
         ! ~12 decades below the pre-port answer (max|S| ~ 5e-2, the
         ! bathymetry step read as a slope; max|uhD| ~ 1e4 m^3/s).
         call check(error, smax <= 1.0e-14_wp, &
                    "a z-level stratification over a staircase has ZERO isopycnal "// &
                    "slope (the bathymetry step is not an isopycnal tilt): "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, fmax <= 1.0e-12_wp*KHTH*DX, &
                    "GM must move nothing in a resting flat stratification: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, maxval(abs(gm%gm_src)) <= 1.0e-20_wp, &
                    "no slope => no GM PE release (to round-off)")
         if (allocated(error)) exit checks
         call check(error, all_finite_3d(sl%n2_u) .and. all_finite_3d(sl%n2_v), &
                    "N^2 must be finite everywhere, fillers included")
      end block checks
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_rest_flat

   ! ------------------------------------------------------------------
   ! Case 2: tilted isopycnals — open-column GM over NSTEP steps
   ! ------------------------------------------------------------------
   subroutine test_tilted_open_column(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :), h_fill0(:, :, :)
      integer, parameter :: NSTEP = 12
      integer :: it, ni, nj, n_closed, n_bad, n_bad_max
      real(wp) :: worst_closed, worst_colsum, max_flux, worst_i1
      real(wp) :: wc_max, cs_rel_max, mf_max, m0, m1, t0s, t1s
      logical :: finite_all, fill_h_kept
      character(len=96) :: msg

      ! Tilt ~1e-3: warm east / north at 2e-5 K/m against GZ = 0.02 K/m.
      call build_case(grid, metrics, ms, sl, gm, eos, tgt, 0.0_wp, 2.0e-5_wp, 1.0e-5_wp)
      ni = grid%nx_total
      nj = grid%ny_total
      call ct%init(grid, nz_ml=NZ)
      allocate (h_fill0, source=ms%h_layer)
      m0 = sum(ms%h_layer)
      t0s = content_sum(ms, ms%idx_temperature)

      wc_max = 0.0_wp
      cs_rel_max = 0.0_wp
      mf_max = 0.0_wp
      n_bad_max = 0
      finite_all = .true.
      fill_h_kept = .true.
      checks: block
         call map_in(ms, sl, gm, ct)
         do it = 1, NSTEP
            call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
            call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
            call continuity_gm_apply(grid, metrics, ct, ms, gm, DT, 1.0_wp, TR_MODE_ADVECT)
            ! The production step tail: restore I1′, then the tripwire scan.
            call ms%enforce_vanished_content(ni, nj)
            call ms%scan_vanished_content(ni, nj, n_bad, worst_i1)
            n_bad_max = max(n_bad_max, n_bad)
            call pull_back(ms, sl, gm)
            call closed_face_census(metrics, ms, gm, ni, nj, n_closed, worst_closed, &
                                    worst_colsum, max_flux)
            wc_max = max(wc_max, worst_closed)
            mf_max = max(mf_max, max_flux)
            if (max_flux > 0.0_wp) cs_rel_max = max(cs_rel_max, worst_colsum/max_flux)
            finite_all = finite_all .and. all_finite_3d(gm%uhD) .and. &
                         all_finite_3d(gm%vhD) .and. all_finite_2d(gm%gm_src) .and. &
                         all_finite_3d(sl%slope_x) .and. all_finite_3d(sl%slope_y) .and. &
                         all_finite_3d(sl%n2_u) .and. all_finite_3d(sl%n2_v) .and. &
                         all_finite_3d(ms%h_layer) .and. &
                         all_finite_3d(ms%tracers(ms%idx_temperature)%hTr)
            fill_h_kept = fill_h_kept .and. &
                          all(merge(ms%h_layer == h_fill0, .true., h_fill0 <= H_VANISHED))
         end do
         call map_out(ms, sl, gm, ct)

         call check(error, n_closed > 0, "the census must see closed face-layers")
         if (allocated(error)) exit checks
         call check(error, mf_max > 0.0_wp, &
                    "tilted isopycnals must drive a non-zero GM transport")
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3,a,es10.3,a)') "max closed |flux| = ", wc_max, &
            " m^3/s (open max ", mf_max, ")"
         call check(error, wc_max == 0.0_wp, &
                    "a CLOSED or filler face-layer must carry EXACTLY zero GM flux: "// &
                    trim(msg))
         if (allocated(error)) exit checks
         call check(error, cs_rel_max < 1.0e-13_wp, &
                    "every face's column-integrated GM transport must be 0 to round-off")
         if (allocated(error)) exit checks
         call check(error, finite_all, &
                    "slopes, N^2, uhD/vhD, gm_src, h, hTr must stay finite everywhere")
         if (allocated(error)) exit checks
         call check(error, fill_h_kept, &
                    "a filler's thickness must not change (no GM mass into or out of it)")
         if (allocated(error)) exit checks
         call check(error, n_bad_max == 0, "I1' must hold after every step")
         if (allocated(error)) exit checks
         m1 = sum(ms%h_layer)
         t1s = content_sum(ms, ms%idx_temperature)
         call check(error, abs(m1 - m0) <= 1.0e-13_wp*m0, "GM must conserve mass")
         if (allocated(error)) exit checks
         call check(error, abs(t1s - t0s) <= 1.0e-13_wp*abs(t0s), &
                    "GM must conserve temperature content")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_tilted_open_column

   ! ------------------------------------------------------------------
   ! Case 3: fillers above the open column (ice draft) — top closure
   ! ------------------------------------------------------------------
   subroutine test_ice_top_closure(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :)
      integer :: ni, nj, n_closed
      real(wp) :: worst_closed, worst_colsum, max_flux
      character(len=96) :: msg

      ! A 120 m draft everywhere: the top two nominal layers are inside the
      ! ice (fillers), the third is a partial top cell.
      call build_case(grid, metrics, ms, sl, gm, eos, tgt, 120.0_wp, 2.0e-5_wp, 1.0e-5_wp)
      ni = grid%nx_total
      nj = grid%ny_total
      checks: block
         call check(error, all(tgt(:, :, NZ) <= H_VANISHED), &
                    "the k = nz layer must be a filler under the draft")
         if (allocated(error)) exit checks

         call map_in(ms, sl, gm)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call pull_back(ms, sl, gm)
         call map_out(ms, sl, gm)

         call closed_face_census(metrics, ms, gm, ni, nj, n_closed, worst_closed, &
                                 worst_colsum, max_flux)
         call check(error, max_flux > 0.0_wp, "the tilted column must carry GM")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(gm%uhD(:, :, NZ))) == 0.0_wp .and. &
                    maxval(abs(gm%vhD(:, :, NZ))) == 0.0_wp, &
                    "the closure must not land in the k = nz filler under the ice")
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3,a,es10.3,a)') "max closed |flux| = ", worst_closed, &
            " m^3/s (open max ", max_flux, ")"
         call check(error, worst_closed == 0.0_wp, &
                    "no GM flux on any closed / filler face-layer under the ice: "// &
                    trim(msg))
         if (allocated(error)) exit checks
         call check(error, worst_colsum <= 1.0e-13_wp*max_flux, &
                    "column-integrated GM transport must still vanish")
      end block checks
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_ice_top_closure

end module test_ocean_gm_zfixed
