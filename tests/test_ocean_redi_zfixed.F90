!! Redi continuous neutral diffusion on `z_fixed` with partial-step CLOSED
!! faces (`&vcoord_nml zfixed_closed_faces`).  Every case runs the
!! production device kernels — `redi_calc_coeffs` + `redi_apply_flux` (+ the
!! I1′ enforcement point and its scan for the multi-step case) — over a
!! hand-specified staircase whose layer thicknesses AND closed-face masks come
!! from the SAME production builders a configured run uses
!! (`ocean_vcoord_z_fixed_target_uniform` → `ocean_vcoord_closed_face_masks`),
!! so "filler" and "closed" have their one production definition.
!!
!! Bottom-up (k = 1 bed, k = nz top).  Linear EOS.  The bed deepens WESTWARD
!! and SOUTHWARD in 30 m / 20 m steps (off the 50 m nominal spacing: partial
!! bottom cells, 0..3 bed fillers per column).
!!
!! The face-level property is asserted on the Phase-A COEFFICIENTS, which is
!! where Redi decides what it pairs: every neutral sublayer with a non-zero
!! effective thickness scatters (Phase B) into native layers `nz+1-KoL` /
!! `nz+1-KoR`, and both must be open at the face and live on both sides
!! (`pairing_census`).  Phase B is then checked on the tracer field itself
!! (fillers untouched, content, extrema, finiteness, I1′).
!!
!!  1. `uniform_ts_zero_flux` — uniform T/S over the staircase (every
!!     neutral-surface comparison is a tie): the Redi update is identically
!!     zero, everything finite.
!!  2. `tilted_isopycnals_open_window` — T warm over the DEEP side, with a
!!     half-compensating S gradient (so T and S vary along the neutral
!!     surfaces) and an isopycnal offset of ~50 m per column — more than a
!!     bed step, so a deep-side layer that is CLOSED at the face (below the
!!     shallow side's bed) has the density of live water across it: NSTEP
!!     steps, no pairing touches a closed face-layer or a filler, fillers'
!!     content untouched by Redi, T and S content conserved to round-off, no
!!     new extrema, all finite, I1′ after every step.  Fails on the pre-port
!!     Redi (full-column pairing: closed face-layers paired with hEff up to
!!     the layer thickness).
!!  3. `ice_draft_open_window` — fillers ABOVE the open column (an ice
!!     draft, `z_top > 0`): the window's top is the topmost open layer, no
!!     pairing reaches a draft filler.  The kernel property only — Redi
!!     under a cavity is still refused at configure.
!!  4. `all_open_matches_full_column` — flat bed, nothing closed: the knob
!!     ON (window = the whole column) gives the knob-OFF answer BITWISE.
!!  5. `open_steps_drained_cell_bounded` — the same tilted staircase with
!!     closed faces OFF (`zstar`'s open steps, MOM6-style: the window is the
!!     whole column, so fillers are paired), and one live partial bed cell
!!     OVER-DRAINED by the continuity step to `h = -8.2e-4 m`, its content
!!     advected with it (`hTr = c·h`) — the cell kind, and the value, the
!!     1-degree Southern Ocean reached at step 12 (a thin partial cell
!!     beside two open filler faces, drained by the GM bolus flux folded into
!!     the same continuity sweeps as the resolved flux).  Redi must read that
!!     layer's concentration by the I1′ rule (`rdb_vl_column_conc`: a
!!     non-live layer reads its donor's), so every live cell stays inside the
!!     initial live range, no cell's content exceeds 2·max|c|·H_NOM, and
!!     everything stays finite.  Phase A (`redi_calc_coeffs`) runs on the
!!     healthy state and Phase B (`redi_apply_flux`) on the drained one, as
!!     in the model (coefficients at the start of the thermo step, the flux
!!     after continuity).  Fails on the floored read `hTr/max(h, 1e-20)`:
!!     measured hT = 2.4e18 in the drained cell and a live T of -3.0e16.
module test_ocean_redi_zfixed
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_redi, only: ocean_redi_t, redi_calc_coeffs, redi_apply_flux
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_redi_zfixed_tests

   integer, parameter :: NG = 2
   integer, parameter :: NXP = 8, NYP = 4, NZ = 8
   real(wp), parameter :: DX = 5000.0_wp
   real(wp), parameter :: H_NOM = 50.0_wp
      !! Nominal z_fixed spacing (m): NZ*H_NOM = 400 m full depth.
   real(wp), parameter :: H_MIN = 1.0e-4_wp
      !! `zstar_h_min` at its default: the inert filler thickness.
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: KHTR = 600.0_wp
      !! `KHTR*DT/DX^2 = 0.043`: well inside the explicit limit.
   real(wp), parameter :: GZ = 0.01_wp
      !! Vertical T gradient (K/m), warm above: stable.
   real(wp), parameter :: GX = -2.0e-4_wp, GY = -1.0e-4_wp
      !! Horizontal T gradient (K/m): warm WEST and SOUTH, over the deep side.
   real(wp), parameter :: SX = 0.5_wp*ALPHA_T*GX/BETA_S
   real(wp), parameter :: SY = 0.5_wp*ALPHA_T*GY/BETA_S
      !! S gradient compensating HALF the T density gradient: the isopycnal
      !! slope is `0.5*|GX|/GZ = 1e-2` (50 m per 5 km column, more than the
      !! 30 m bed step), and T, S both vary along the neutral surfaces.

contains

   subroutine collect_ocean_redi_zfixed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("uniform_ts_zero_flux", test_uniform_ts), &
                  new_unittest("tilted_isopycnals_open_window", test_tilted_open_window), &
                  new_unittest("ice_draft_open_window", test_ice_draft), &
                  new_unittest("all_open_matches_full_column", test_all_open_bitwise), &
                  new_unittest("open_steps_drained_cell_bounded", test_open_drained_cell) &
                  ]
   end subroutine collect_ocean_redi_zfixed_tests

   ! ------------------------------------------------------------------
   ! Setup helpers
   ! ------------------------------------------------------------------

   subroutine build_case(grid, metrics, ms, rd, eos, tgt, draft, gx_, gy_, gz_, &
                         sx_, sy_, staircase, closed)
      !! Grid + (closed-face) metrics + a z_fixed state + a Redi slot.
      !!
      !! `staircase`: the bed steps 30 m per column eastward and 20 m per row
      !! northward (else a flat 400 m bed); ghost columns replicate the
      !! nearest physical column.  `draft` (m, >= 0) puts an ice base over
      !! the whole domain.  The thickness field is the production `z_fixed`
      !! target at eta = 0; with `closed` the masks are built from that same
      !! target by the production builder, else the metrics keep the
      !! `(1,1,1)` placeholders and the knob is off.
      !!
      !! T = T0 + gz*z + gx*x + gy*y and S = S0 + sx*x + sy*y at the NOMINAL
      !! layer height; fillers are then given their donor's concentration
      !! and the host twin of the enforcement point establishes I1′.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_redi_t), intent(inout) :: rd
      type(eos_t), intent(out) :: eos
      real(wp), allocatable, intent(out) :: tgt(:, :, :)
      real(wp), intent(in) :: draft, gx_, gy_, gz_, sx_, sy_
      logical, intent(in) :: staircase, closed
      real(wp), allocatable :: tot_h(:, :), eta0(:, :), z_top(:, :)
      integer :: i, j, k, ni, nj, ip, jp, it
      real(wp) :: zc, x, y, c_val

      call grid%init(NXP, NYP, NG, DX, DX)
      ni = grid%nx_total
      nj = grid%ny_total

      if (closed) then
         ! `nz_closed` grows open_u/open_v BEFORE the device map.
         call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
      else
         call make_cartesian_metrics(metrics, grid)
      end if

      allocate (tot_h(ni, nj), eta0(ni, nj), z_top(ni, nj))
      allocate (tgt(ni, nj, NZ), source=0.0_wp)
      eta0 = 0.0_wp
      z_top = draft
      do j = 1, nj
         jp = min(max(j - NG, 1), NYP)
         do i = 1, ni
            ip = min(max(i - NG, 1), NXP)
            tot_h(i, j) = real(NZ, wp)*H_NOM - draft
            if (staircase) tot_h(i, j) = tot_h(i, j) - 30.0_wp*real(ip - 1, wp) &
                                         - 20.0_wp*real(jp - 1, wp)
         end do
      end do
      call ocean_vcoord_z_fixed_target_uniform(tgt, tot_h, eta0, z_top, &
                                               ni, nj, NZ, H_NOM, H_MIN)
      if (closed) then
         call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                             tgt, ni, nj, NZ, H_VANISHED)
         metrics%use_closed_faces = .true.
         ! The masks are device-present (grown before the map), so the
         ! builder wrote the device copy; pull it back for the host census.
         !$omp target update from(metrics%open_u, metrics%open_v)
      end if

      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = tgt
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         zc = (real(k, wp) - 0.5_wp)*H_NOM          ! nominal height above 400 m
         do j = 1, nj
            jp = min(max(j - NG, 1), NYP)
            y = real(jp, wp)*DX
            do i = 1, ni
               ip = min(max(i - NG, 1), NXP)
               x = real(ip, wp)*DX
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  (T0 + gz_*zc + gx_*x + gy_*y)*tgt(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = (S0 + sx_*x + sy_*y)*tgt(i, j, k)
            end do
         end do
      end do
      ! Hand each filler its donor's concentration BEFORE the I1' sweep, so
      ! the sweep's pooling mixes equal concentrations and does not shift
      ! the donor (cf. test_ocean_gm_zfixed).
      do it = 1, 2
         associate (q => ms%tracers(merge(ms%idx_temperature, ms%idx_salinity, it == 1))%hTr)
            do j = 1, nj
               do i = 1, ni
                  c_val = -huge(1.0_wp)
                  do k = NZ, 1, -1
                     if (tgt(i, j, k) > H_VANISHED) then
                        c_val = q(i, j, k)/tgt(i, j, k)
                     else if (c_val > -huge(1.0_wp)) then
                        q(i, j, k) = c_val*tgt(i, j, k)
                     end if
                  end do
               end do
            end do
         end associate
      end do
      call ms%enforce_vanished_content_host(ni, nj)

      call rd%init(grid, nz_ml=NZ)
      rd%enable = .true.
      rd%continuous = .true.
      rd%khtr = KHTR

      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine build_case

   subroutine map_in(ms, rd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_redi_t), intent(inout) :: rd
      !$omp target enter data map(to: ms)
      call ms%enter_data()
      !$omp target enter data map(to: rd)
      call rd%enter_data()
   end subroutine map_in

   subroutine pull_back(ms, rd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_redi_t), intent(inout) :: rd
      !$omp target update from(rd%uKoL, rd%uKoR, rd%uhEff, rd%uPoL, rd%uPoR)
      !$omp target update from(rd%vKoL, rd%vKoR, rd%vhEff, rd%vPoL, rd%vPoR)
      !$omp target update from(ms%tracers(ms%idx_temperature)%hTr)
      !$omp target update from(ms%tracers(ms%idx_salinity)%hTr)
   end subroutine pull_back

   subroutine map_out(ms, rd, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_redi_t), intent(inout) :: rd
      type(ocean_metrics_t), intent(inout) :: metrics
      call rd%exit_data()
      !$omp target exit data map(delete: rd)
      call ms%exit_data()
      !$omp target exit data map(delete: ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   pure logical function all_finite_3d(a)
      real(wp), intent(in) :: a(:, :, :)
      all_finite_3d = all(ieee_is_finite(a))
   end function all_finite_3d

   subroutine pairing_census(metrics, ms, rd, ni, nj, n_pair, n_bad, worst_heff)
      !! Host census of the Phase-A pairing over EVERY interior face of the
      !! array: each neutral sublayer with `hEff /= 0` scatters into native
      !! layers `NZ+1-KoL` (left column) and `NZ+1-KoR` (right column); each
      !! must be open at the face and live on BOTH sides.  `n_pair` counts
      !! the active sublayers (non-vacuity), `n_bad` the violations and
      !! `worst_heff` the largest effective thickness (m) of a violating one.
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_redi_t), intent(in) :: rd
      integer, intent(in) :: ni, nj
      integer, intent(out) :: n_pair, n_bad
      real(wp), intent(out) :: worst_heff
      integer :: i, j, ks, kl, kr

      n_pair = 0
      n_bad = 0
      worst_heff = 0.0_wp
      do j = 1, nj
         do i = 2, ni
            do ks = 1, rd%nsurf - 1
               if (rd%uhEff(i, j, ks) == 0.0_wp) cycle
               n_pair = n_pair + 1
               kl = NZ + 1 - rd%uKoL(i, j, ks)
               kr = NZ + 1 - rd%uKoR(i, j, ks)
               if (.not. (face_layer_open(metrics%open_u(i, j, kl), ms%h_layer(i - 1, j, kl), &
                                          ms%h_layer(i, j, kl)) .and. &
                          face_layer_open(metrics%open_u(i, j, kr), ms%h_layer(i - 1, j, kr), &
                                          ms%h_layer(i, j, kr)))) then
                  n_bad = n_bad + 1
                  worst_heff = max(worst_heff, abs(rd%uhEff(i, j, ks)))
               end if
            end do
         end do
      end do
      do j = 2, nj
         do i = 1, ni
            do ks = 1, rd%nsurf - 1
               if (rd%vhEff(i, j, ks) == 0.0_wp) cycle
               n_pair = n_pair + 1
               kl = NZ + 1 - rd%vKoL(i, j, ks)
               kr = NZ + 1 - rd%vKoR(i, j, ks)
               if (.not. (face_layer_open(metrics%open_v(i, j, kl), ms%h_layer(i, j - 1, kl), &
                                          ms%h_layer(i, j, kl)) .and. &
                          face_layer_open(metrics%open_v(i, j, kr), ms%h_layer(i, j - 1, kr), &
                                          ms%h_layer(i, j, kr)))) then
                  n_bad = n_bad + 1
                  worst_heff = max(worst_heff, abs(rd%vhEff(i, j, ks)))
               end if
            end do
         end do
      end do
   end subroutine pairing_census

   pure logical function face_layer_open(open_f, h_a, h_b)
      !! A face-layer a Redi flux may cross: open in the static mask and
      !! live (`> H_VANISHED`) on both sides.
      real(wp), intent(in) :: open_f, h_a, h_b
      face_layer_open = open_f > 0.5_wp .and. h_a > H_VANISHED .and. h_b > H_VANISHED
   end function face_layer_open

   real(wp) function phys_content(ms, idx) result(s)
      !! Content summed over the PHYSICAL cells (the ghost columns trade
      !! among themselves; the physical-domain walls are closed to Redi).
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx
      s = sum(ms%tracers(idx)%hTr(NG + 1:NG + NXP, NG + 1:NG + NYP, :))
   end function phys_content

   subroutine live_range(ms, idx, cmin, cmax)
      !! Min / max concentration over the LIVE physical cells.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx
      real(wp), intent(out) :: cmin, cmax
      integer :: i, j, k
      real(wp) :: c
      cmin = huge(1.0_wp)
      cmax = -huge(1.0_wp)
      do k = 1, NZ
         do j = NG + 1, NG + NYP
            do i = NG + 1, NG + NXP
               if (ms%h_layer(i, j, k) <= H_VANISHED) cycle
               c = ms%tracers(idx)%hTr(i, j, k)/ms%h_layer(i, j, k)
               cmin = min(cmin, c)
               cmax = max(cmax, c)
            end do
         end do
      end do
   end subroutine live_range

   ! ------------------------------------------------------------------
   ! Case 1: uniform T/S — Redi moves nothing
   ! ------------------------------------------------------------------
   subroutine test_uniform_ts(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :), t_in(:, :, :), s_in(:, :, :)
      real(wp) :: dmax
      character(len=96) :: msg

      call build_case(grid, metrics, ms, rd, eos, tgt, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                      0.0_wp, 0.0_wp, .true., .true.)
      checks: block
         call check(error, count(tgt <= H_VANISHED) > 0, "the staircase must carry bed fillers")
         if (allocated(error)) exit checks
         call check(error, count(metrics%open_u < 0.5_wp) > 0 .and. &
                    count(metrics%open_v < 0.5_wp) > 0, &
                    "the staircase must close u- AND v-face layers")
         if (allocated(error)) exit checks
         allocate (t_in, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (s_in, source=ms%tracers(ms%idx_salinity)%hTr)

         call map_in(ms, rd)
         call redi_calc_coeffs(grid, metrics, eos, rd, ms)
         call redi_apply_flux(grid, metrics, rd, ms, DT)
         call pull_back(ms, rd)
         call map_out(ms, rd, metrics)

         call check(error, all_finite_3d(ms%tracers(ms%idx_temperature)%hTr) .and. &
                    all_finite_3d(ms%tracers(ms%idx_salinity)%hTr) .and. &
                    all_finite_3d(rd%uhEff) .and. all_finite_3d(rd%vhEff) .and. &
                    all_finite_3d(rd%uPoL) .and. all_finite_3d(rd%vPoR), &
                    "every field must stay finite (all-tie neutral sweep)")
         if (allocated(error)) exit checks
         dmax = max(maxval(abs(ms%tracers(ms%idx_temperature)%hTr - t_in)), &
                    maxval(abs(ms%tracers(ms%idx_salinity)%hTr - s_in)))
         write (msg, '(a,es10.3)') "max|d hTr| = ", dmax
         call check(error, dmax == 0.0_wp, &
                    "uniform T/S: the Redi update must be identically zero: "//trim(msg))
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_uniform_ts

   ! ------------------------------------------------------------------
   ! Case 2: tilted isopycnals over the staircase — NSTEP steps
   ! ------------------------------------------------------------------
   subroutine test_tilted_open_window(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :), t_prev(:, :, :), s_prev(:, :, :)
      real(wp), allocatable :: t_init(:, :, :)
      integer, parameter :: NSTEP = 10
      integer :: it, ni, nj, n_pair, n_bad, n_pair_tot, n_bad_tot, n_i1, n_i1_max
      real(wp) :: worst_heff, wh_max, worst_i1, t0s, s0s, t1s, s1s
      real(wp) :: tmin0, tmax0, smin0, smax0, tmin, tmax, smin, smax, ext, dtmax
      logical :: finite_all, fill_kept
      character(len=128) :: msg

      call build_case(grid, metrics, ms, rd, eos, tgt, 0.0_wp, GX, GY, GZ, SX, SY, &
                      .true., .true.)
      ni = grid%nx_total
      nj = grid%ny_total
      t0s = phys_content(ms, ms%idx_temperature)
      s0s = phys_content(ms, ms%idx_salinity)
      call live_range(ms, ms%idx_temperature, tmin0, tmax0)
      call live_range(ms, ms%idx_salinity, smin0, smax0)
      allocate (t_init, source=ms%tracers(ms%idx_temperature)%hTr)

      n_pair_tot = 0
      n_bad_tot = 0
      n_i1_max = 0
      wh_max = 0.0_wp
      ext = 0.0_wp
      finite_all = .true.
      fill_kept = .true.
      checks: block
         call map_in(ms, rd)
         do it = 1, NSTEP
            allocate (t_prev, source=ms%tracers(ms%idx_temperature)%hTr)
            allocate (s_prev, source=ms%tracers(ms%idx_salinity)%hTr)
            call redi_calc_coeffs(grid, metrics, eos, rd, ms)
            call redi_apply_flux(grid, metrics, rd, ms, DT)
            call pull_back(ms, rd)
            ! Redi alone, before the step tail: a filler's content is
            ! untouched (no flux into or out of it).
            fill_kept = fill_kept .and. &
                        all(merge(ms%tracers(ms%idx_temperature)%hTr == t_prev, .true., &
                                  tgt <= H_VANISHED)) .and. &
                        all(merge(ms%tracers(ms%idx_salinity)%hTr == s_prev, .true., &
                                  tgt <= H_VANISHED))
            call pairing_census(metrics, ms, rd, ni, nj, n_pair, n_bad, worst_heff)
            n_pair_tot = n_pair_tot + n_pair
            n_bad_tot = n_bad_tot + n_bad
            wh_max = max(wh_max, worst_heff)
            ! The production step tail: restore I1′ (Redi moved the donors'
            ! concentration), then the tripwire scan.
            call ms%enforce_vanished_content(ni, nj)
            call ms%scan_vanished_content(ni, nj, n_i1, worst_i1)
            n_i1_max = max(n_i1_max, n_i1)
            !$omp target update from(ms%tracers(ms%idx_temperature)%hTr)
            !$omp target update from(ms%tracers(ms%idx_salinity)%hTr)
            finite_all = finite_all .and. &
                         all_finite_3d(ms%tracers(ms%idx_temperature)%hTr) .and. &
                         all_finite_3d(ms%tracers(ms%idx_salinity)%hTr) .and. &
                         all_finite_3d(rd%uhEff) .and. all_finite_3d(rd%vhEff) .and. &
                         all_finite_3d(rd%uPoL) .and. all_finite_3d(rd%uPoR) .and. &
                         all_finite_3d(rd%vPoL) .and. all_finite_3d(rd%vPoR)
            call live_range(ms, ms%idx_temperature, tmin, tmax)
            call live_range(ms, ms%idx_salinity, smin, smax)
            ext = max(ext, tmax - tmax0, tmin0 - tmin, smax - smax0, smin0 - smin)
            deallocate (t_prev, s_prev)
         end do
         call map_out(ms, rd, metrics)

         call check(error, n_pair_tot > 0, "the census must see active neutral sublayers")
         if (allocated(error)) exit checks
         write (msg, '(a,i0,a,i0,a,es10.3,a)') "violations ", n_bad_tot, " of ", n_pair_tot, &
            " active sublayers, worst hEff ", wh_max, " m"
         call check(error, n_bad_tot == 0, &
                    "no neutral sublayer may pair a CLOSED face-layer or a filler: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, fill_kept, &
                    "Redi must not change a filler's content (no flux into or out of it)")
         if (allocated(error)) exit checks
         call check(error, finite_all, "hTr and the Redi coefficients must stay finite")
         if (allocated(error)) exit checks
         call check(error, n_i1_max == 0, "I1' must hold after every step")
         if (allocated(error)) exit checks
         dtmax = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - t_init))
         call check(error, dtmax > 1.0e-3_wp, &
                    "tilted isopycnals with spice must drive a non-trivial Redi flux")
         if (allocated(error)) exit checks
         t1s = phys_content(ms, ms%idx_temperature)
         s1s = phys_content(ms, ms%idx_salinity)
         write (msg, '(a,es10.3,a,es10.3)') "rel dT-content ", abs(t1s - t0s)/abs(t0s), &
            ", rel dS-content ", abs(s1s - s0s)/abs(s0s)
         call check(error, abs(t1s - t0s) <= 1.0e-13_wp*abs(t0s) .and. &
                    abs(s1s - s0s) <= 1.0e-13_wp*abs(s0s), &
                    "Redi must conserve T and S content: "//trim(msg))
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3)') "largest excursion beyond the initial range ", ext
         call check(error, ext <= 1.0e-12_wp*max(tmax0, smax0), &
                    "Redi must not create new extrema (live T and S stay in range): "// &
                    trim(msg))
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_tilted_open_window

   ! ------------------------------------------------------------------
   ! Case 3: fillers above the open column (ice draft)
   ! ------------------------------------------------------------------
   subroutine test_ice_draft(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :), t_prev(:, :, :)
      integer :: ni, nj, n_pair, n_bad
      real(wp) :: worst_heff, t0s, t1s
      character(len=128) :: msg

      ! A 120 m draft everywhere: the top two nominal layers are inside the
      ! ice (fillers), the third is a partial top cell.
      call build_case(grid, metrics, ms, rd, eos, tgt, 120.0_wp, GX, GY, GZ, SX, SY, &
                      .true., .true.)
      ni = grid%nx_total
      nj = grid%ny_total
      checks: block
         call check(error, all(tgt(:, :, NZ) <= H_VANISHED), &
                    "the k = nz layer must be a filler under the draft")
         if (allocated(error)) exit checks
         allocate (t_prev, source=ms%tracers(ms%idx_temperature)%hTr)
         t0s = phys_content(ms, ms%idx_temperature)

         call map_in(ms, rd)
         call redi_calc_coeffs(grid, metrics, eos, rd, ms)
         call redi_apply_flux(grid, metrics, rd, ms, DT)
         call pull_back(ms, rd)
         ! Census BEFORE the unmap: `map_out` destroys the metrics (masks).
         call pairing_census(metrics, ms, rd, ni, nj, n_pair, n_bad, worst_heff)
         call map_out(ms, rd, metrics)

         call check(error, n_pair > 0, "the open window under the ice must carry Redi")
         if (allocated(error)) exit checks
         write (msg, '(a,i0,a,i0,a,es10.3,a)') "violations ", n_bad, " of ", n_pair, &
            " active sublayers, worst hEff ", worst_heff, " m"
         call check(error, n_bad == 0, &
                    "no neutral sublayer may reach a draft or bed filler: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, all(merge(ms%tracers(ms%idx_temperature)%hTr == t_prev, .true., &
                                     tgt <= H_VANISHED)), &
                    "Redi must not change a draft filler's content")
         if (allocated(error)) exit checks
         call check(error, all_finite_3d(ms%tracers(ms%idx_temperature)%hTr), &
                    "T must stay finite under the draft")
         if (allocated(error)) exit checks
         t1s = phys_content(ms, ms%idx_temperature)
         call check(error, abs(t1s - t0s) <= 1.0e-13_wp*abs(t0s), &
                    "Redi must conserve T content under the draft")
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_ice_draft

   ! ------------------------------------------------------------------
   ! Case 4: nothing closed — knob on is the knob-off answer, bitwise
   ! ------------------------------------------------------------------
   subroutine test_all_open_bitwise(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :), t_off(:, :, :), s_off(:, :, :), t_in(:, :, :)
      integer :: pass
      logical :: closed

      do pass = 1, 2
         closed = pass == 2
         call build_case(grid, metrics, ms, rd, eos, tgt, 0.0_wp, GX, GY, GZ, SX, SY, &
                         .false., closed)
         if (pass == 1) allocate (t_in, source=ms%tracers(ms%idx_temperature)%hTr)
         call map_in(ms, rd)
         call redi_calc_coeffs(grid, metrics, eos, rd, ms)
         call redi_apply_flux(grid, metrics, rd, ms, DT)
         call pull_back(ms, rd)
         call map_out(ms, rd, metrics)
         if (pass == 1) then
            allocate (t_off, source=ms%tracers(ms%idx_temperature)%hTr)
            allocate (s_off, source=ms%tracers(ms%idx_salinity)%hTr)
         end if
         if (pass == 2) exit
         call rd%destroy()
         call ms%destroy()
      end do
      checks: block
         call check(error, maxval(abs(t_off - t_in)) > 0.0_wp, &
                    "the flat-bed case must carry a Redi flux (non-vacuous)")
         if (allocated(error)) exit checks
         call check(error, all(ms%tracers(ms%idx_temperature)%hTr == t_off) .and. &
                    all(ms%tracers(ms%idx_salinity)%hTr == s_off), &
                    "all-open window must reproduce the full-column Redi bitwise")
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_all_open_bitwise

   ! ------------------------------------------------------------------
   ! Case 5: open steps, an over-drained live partial cell
   ! ------------------------------------------------------------------
   subroutine test_open_drained_cell(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      real(wp), allocatable :: tgt(:, :, :)
      real(wp) :: tmin0, tmax0, smin0, smax0, tmin, tmax, smin, smax, tol, c_t, c_s
      integer :: i, j, k, id, jd, kd
      character(len=200) :: msg

      call build_case(grid, metrics, ms, rd, eos, tgt, 0.0_wp, GX, GY, GZ, SX, SY, &
                      .true., .false.)
      checks: block
         call live_range(ms, ms%idx_temperature, tmin0, tmax0)
         call live_range(ms, ms%idx_salinity, smin0, smax0)
         ! The live partial bed cell sitting directly on the topmost bed
         ! filler of an interior column (the thin cell GM drains through its
         ! open filler faces).
         id = 0
         do j = NG + 2, NG + NYP - 1
            do i = NG + 2, NG + NXP - 1
               do k = NZ - 1, 1, -1
                  if (id == 0 .and. tgt(i, j, k) <= H_VANISHED .and. &
                      tgt(i, j, k + 1) > H_VANISHED) then
                     id = i; jd = j; kd = k
                  end if
               end do
            end do
         end do
         call check(error, id > 0, "the open staircase must carry an interior bed filler")
         if (allocated(error)) exit checks
         kd = kd + 1
         c_t = ms%tracers(ms%idx_temperature)%hTr(id, jd, kd)/ms%h_layer(id, jd, kd)
         c_s = ms%tracers(ms%idx_salinity)%hTr(id, jd, kd)/ms%h_layer(id, jd, kd)

         ! Phase A on the healthy state (the model runs it at the start of
         ! the thermo step) ...
         call map_in(ms, rd)
         call redi_calc_coeffs(grid, metrics, eos, rd, ms)
         ! ... then continuity over-drains the cell to the observed -8.2e-4 m,
         ! its content advected with it, and Phase B runs on that state.
         !$omp target update from(ms%h_layer, ms%tracers(ms%idx_temperature)%hTr)
         !$omp target update from(ms%tracers(ms%idx_salinity)%hTr)
         ms%h_layer(id, jd, kd) = -8.2e-4_wp
         ms%tracers(ms%idx_temperature)%hTr(id, jd, kd) = c_t*ms%h_layer(id, jd, kd)
         ms%tracers(ms%idx_salinity)%hTr(id, jd, kd) = c_s*ms%h_layer(id, jd, kd)
         !$omp target update to(ms%h_layer, ms%tracers(ms%idx_temperature)%hTr)
         !$omp target update to(ms%tracers(ms%idx_salinity)%hTr)
         call redi_apply_flux(grid, metrics, rd, ms, DT)
         call pull_back(ms, rd)
         call map_out(ms, rd, metrics)

         call check(error, all_finite_3d(ms%tracers(ms%idx_temperature)%hTr) .and. &
                    all_finite_3d(ms%tracers(ms%idx_salinity)%hTr), &
                    "open steps, drained cell: the Redi update must stay finite")
         if (allocated(error)) exit checks
         ! Content bound on EVERY cell, the drained one included: twice the
         ! largest concentration times the thickest layer.  The initial state
         ! sits ON `max|c|*H_NOM` (a full layer at the warmest T), so the bound
         ! needs room for round-off (GPU FMA); the defect is 1e18.
         write (msg, '(a,2es12.4)') "max|hT|, max|hS| = ", &
            maxval(abs(ms%tracers(ms%idx_temperature)%hTr)), &
            maxval(abs(ms%tracers(ms%idx_salinity)%hTr))
         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr)) <= &
                    2.0_wp*max(abs(tmin0), abs(tmax0))*H_NOM .and. &
                    maxval(abs(ms%tracers(ms%idx_salinity)%hTr)) <= &
                    2.0_wp*max(abs(smin0), abs(smax0))*H_NOM, &
                    "open steps, drained cell: unbounded content: "//trim(msg))
         if (allocated(error)) exit checks
         call live_range(ms, ms%idx_temperature, tmin, tmax)
         call live_range(ms, ms%idx_salinity, smin, smax)
         tol = 1.0e-9_wp
         write (msg, '(a,4es12.4,a,4es12.4)') "T before/after ", tmin0, tmax0, tmin, tmax, &
            "  S before/after ", smin0, smax0, smin, smax
         call check(error, tmin >= tmin0 - tol .and. tmax <= tmax0 + tol .and. &
                    smin >= smin0 - tol .and. smax <= smax0 + tol, &
                    "open steps, drained cell: Redi made new live extrema: "//trim(msg))
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_open_drained_cell

end module test_ocean_redi_zfixed
