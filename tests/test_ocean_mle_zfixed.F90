!! Fox-Kemper mixed-layer-eddy restratification on `z_fixed` with
!! partial-step CLOSED faces (`&vcoord_nml zfixed_closed_faces`).  Every
!! case runs the production device kernels — `mle_compute_transports`
!! (+ `continuity_tracer_step_split(mle=mle)` and the I1′ enforcement point
!! for the multi-step case) — over a hand-specified staircase whose layer
!! thicknesses AND closed-face masks come from the SAME production builders
!! a configured run uses (`ocean_vcoord_z_fixed_target_uniform` →
!! `ocean_vcoord_closed_face_masks`), so "filler" and "closed" have their
!! one production definition.
!!
!! Bottom-up (k = 1 bed, k = nz top).  `rho_layer` is a linear function of
!! T on live layers; fillers carry a deliberately WRONG density (the EOS
!! reference stands in there in production), so a walk that read them
!! would show.  The MLD (200 m) is deeper than the shallow columns of the
!! staircase, so the open-column `H_vel` clamp is exercised.
!!
!!  1. `closed_faces_zero_flux` — an ML buoyancy front in x AND y over the
!!     staircase: a closed or filler face-layer carries EXACTLY zero MLE
!!     transport, every face column sums to zero to round-off, the
!!     transport is non-trivial.  Fails on the pre-port kernel (full-column
!!     weights put O(1e4) m^3/s on closed face-layers).
!!  2. `ice_top_open_column` — fillers ABOVE the open column (a 120 m ice
!!     draft): sigma = 0 at the ice base, nothing on the `k = nz` filler,
!!     closed layers zero, column sums zero.  The kernel property only
!!     (MLE needs EPBL, which is refused under a cavity).
!!  3. `restrat_front_staircase` — NSTEP MLE + continuity steps from a
!!     vertically uniform ML front: the ML stratification increases from
!!     exactly zero, no new T extrema, mass + T content conserved, filler
!!     thickness unchanged, closed face-layers exactly zero every step, I1′
!!     after every step.
!!  4. `flat_bed_knob_on_equals_off` — no closed face anywhere: the open
!!     path is the full-column form BIT FOR BIT (`uhml`/`vhml` identical
!!     with `use_closed_faces` on and off).
module test_ocean_mle_zfixed
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_epbl, only: ocean_epbl_t
   use rdb_ocean_mle, only: ocean_mle_t, mle_compute_transports
   use rdb_continuity, only: continuity_t, continuity_tracer_step_split
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_mle_zfixed_tests

   integer, parameter :: NG = 2
   integer, parameter :: NXP = 8, NYP = 4, NZ = 8
   real(wp), parameter :: DX = 20000.0_wp
   real(wp), parameter :: H_NOM = 50.0_wp
      !! Nominal z_fixed spacing (m): NZ*H_NOM = 400 m full depth.
   real(wp), parameter :: H_MIN = 1.0e-4_wp
      !! `zstar_h_min` at its default: the inert filler thickness.
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_RHO = 0.2_wp
      !! Linear EOS slope (kg/m^3/K): rho = RHO0 - ALPHA_RHO*(T - T0).
   real(wp), parameter :: RHO_FILLER = 1000.0_wp
      !! A wrong density on fillers: reading it would show in b_bar.
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: MLD = 200.0_wp
   real(wp), parameter :: F0 = 1.0e-4_wp
   real(wp), parameter :: CE = 0.0625_wp

contains

   subroutine collect_ocean_mle_zfixed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("closed_faces_zero_flux", test_closed_zero), &
                  new_unittest("ice_top_open_column", test_ice_top), &
                  new_unittest("restrat_front_staircase", test_restrat), &
                  new_unittest("flat_bed_knob_on_equals_off", test_flat_bitid) &
                  ]
   end subroutine collect_ocean_mle_zfixed_tests

   ! ------------------------------------------------------------------
   ! Setup helpers
   ! ------------------------------------------------------------------

   subroutine build_case(grid, metrics, ms, epbl, mle, tgt, draft, flat)
      !! Grid + closed-face metrics + a z_fixed staircase state with an ML
      !! front.  Bed depth steps down eastward AND northward (30 m per
      !! column, 20 m per row, off the 50 m spacing: partial bottom cells,
      !! 0..5 bed fillers); `flat` keeps every column 400 m deep (no
      !! filler, no closed face).  `draft` (m) puts an ice base over the
      !! whole domain.  Ghost columns replicate the nearest physical one.
      !!
      !! T: in the top MLD of NOMINAL depth a vertically uniform front,
      !! 12 + 3·[ip > 4] + 1.5·[jp > 2]; below it a stable profile colder
      !! than any ML water.  Fillers take their donor's T (I1′).
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_mle_t), intent(inout) :: mle
      real(wp), allocatable, intent(out) :: tgt(:, :, :)
      real(wp), intent(in) :: draft
      logical, intent(in) :: flat
      real(wp), allocatable :: tot_h(:, :), eta0(:, :), z_top(:, :)
      integer :: i, j, k, ni, nj, ip, jp
      real(wp) :: zdep, t_val

      call grid%init(NXP, NYP, NG, DX, DX)
      ni = grid%nx_total
      nj = grid%ny_total
      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)

      allocate (tot_h(ni, nj), eta0(ni, nj), z_top(ni, nj))
      allocate (tgt(ni, nj, NZ), source=0.0_wp)
      eta0 = 0.0_wp
      z_top = draft
      do j = 1, nj
         jp = min(max(j - NG, 1), NYP)
         do i = 1, ni
            ip = min(max(i - NG, 1), NXP)
            if (flat) then
               tot_h(i, j) = real(NZ, wp)*H_NOM - draft
            else
               tot_h(i, j) = (real(NZ, wp)*H_NOM - 30.0_wp*real(ip - 1, wp) &
                              - 20.0_wp*real(jp - 1, wp)) - draft
            end if
         end do
      end do
      call ocean_vcoord_z_fixed_target_uniform(tgt, tot_h, eta0, z_top, &
                                               ni, nj, NZ, H_NOM, H_MIN)
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                          tgt, ni, nj, NZ, H_VANISHED)
      metrics%use_closed_faces = .true.
      !$omp target update from(metrics%open_u, metrics%open_v)

      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = tgt
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         zdep = (real(NZ - k, wp) + 0.5_wp)*H_NOM    ! nominal centre depth
         do j = 1, nj
            jp = min(max(j - NG, 1), NYP)
            do i = 1, ni
               ip = min(max(i - NG, 1), NXP)
               if (zdep < MLD) then
                  t_val = 12.0_wp
                  if (ip > NXP/2) t_val = t_val + 3.0_wp
                  if (jp > NYP/2) t_val = t_val + 1.5_wp
               else
                  t_val = 11.0_wp - 0.005_wp*(zdep - MLD)
               end if
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*tgt(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*tgt(i, j, k)
            end do
         end do
      end do
      ! Fillers take their donor's T before the I1' sweep (see GM twin).
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
      call set_rho_host(ms)

      call epbl%init(grid, nz_ml=NZ)
      epbl%rho0 = RHO0
      epbl%mld = MLD
      epbl%f_centre = F0

      call mle%init(grid, nz_ml=NZ)
      mle%enable = .true.
      mle%ce = CE
      mle%f_floor = 1.0e-5_wp
   end subroutine build_case

   subroutine set_rho_host(ms)
      !! Linear EOS on live layers; RHO_FILLER on fillers.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j, k
      real(wp) :: t
      do k = 1, size(ms%h_layer, 3)
         do j = 1, size(ms%h_layer, 2)
            do i = 1, size(ms%h_layer, 1)
               if (ms%h_layer(i, j, k) > H_VANISHED) then
                  t = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ms%h_layer(i, j, k)
                  ms%rho_layer(i, j, k) = RHO0 - ALPHA_RHO*(t - T0)
               else
                  ms%rho_layer(i, j, k) = RHO_FILLER
               end if
            end do
         end do
      end do
   end subroutine set_rho_host

   subroutine map_in(ms, epbl, mle, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_mle_t), intent(inout) :: mle
      type(continuity_t), intent(inout), optional :: ct
      !$omp target enter data map(to: ms, epbl, mle)
      call ms%enter_data()
      call epbl%enter_data()
      call mle%enter_data()
      if (present(ct)) then
         !$omp target enter data map(to: ct)
         call ct%enter_data()
      end if
   end subroutine map_in

   subroutine map_out(ms, epbl, mle, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_mle_t), intent(inout) :: mle
      type(continuity_t), intent(inout), optional :: ct
      if (present(ct)) then
         call ct%exit_data()
         !$omp target exit data map(delete: ct)
      end if
      call mle%exit_data()
      call epbl%exit_data()
      call ms%exit_data()
      !$omp target exit data map(delete: ms, epbl, mle)
   end subroutine map_out

   subroutine teardown(metrics, ms, epbl, mle)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_epbl_t), intent(inout) :: epbl
      type(ocean_mle_t), intent(inout) :: mle
      call mle%destroy()
      call epbl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown

   subroutine census(metrics, ms, mle, ni, nj, n_closed, worst_closed, &
                     worst_colsum, max_flux)
      !! Host census over EVERY interior face: the largest |MLE transport|
      !! on a face-layer that is closed (`open == 0`) or touches a filler on
      !! either side, the largest |Sum_k transport| of any face column, and
      !! the largest |transport| overall.
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_mle_t), intent(in) :: mle
      integer, intent(in) :: ni, nj
      integer, intent(out) :: n_closed
      real(wp), intent(out) :: worst_closed, worst_colsum, max_flux
      integer :: i, j, k
      real(wp) :: csum
      logical :: shut

      n_closed = 0
      worst_closed = 0.0_wp
      worst_colsum = 0.0_wp
      max_flux = max(maxval(abs(mle%uhml)), maxval(abs(mle%vhml)))
      do j = 1, nj
         do i = 2, ni
            csum = 0.0_wp
            do k = 1, NZ
               csum = csum + mle%uhml(i, j, k)
               shut = metrics%open_u(i, j, k) < 0.5_wp .or. &
                      ms%h_layer(i - 1, j, k) <= H_VANISHED .or. &
                      ms%h_layer(i, j, k) <= H_VANISHED
               if (shut) then
                  n_closed = n_closed + 1
                  worst_closed = max(worst_closed, abs(mle%uhml(i, j, k)))
               end if
            end do
            worst_colsum = max(worst_colsum, abs(csum))
         end do
      end do
      do j = 2, nj
         do i = 1, ni
            csum = 0.0_wp
            do k = 1, NZ
               csum = csum + mle%vhml(i, j, k)
               shut = metrics%open_v(i, j, k) < 0.5_wp .or. &
                      ms%h_layer(i, j - 1, k) <= H_VANISHED .or. &
                      ms%h_layer(i, j, k) <= H_VANISHED
               if (shut) then
                  n_closed = n_closed + 1
                  worst_closed = max(worst_closed, abs(mle%vhml(i, j, k)))
               end if
            end do
            worst_colsum = max(worst_colsum, abs(csum))
         end do
      end do
   end subroutine census

   real(wp) function ml_strat(ms, ni, nj) result(s)
      !! Mean over physical columns of T(top live layer) - T(the live layer
      !! 3 below it, i.e. ~150 m down inside the 200 m ML).  Exactly 0 on
      !! the vertically uniform initial ML.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: ni, nj
      integer :: i, j, k, kt, n
      real(wp) :: tt, tb
      s = 0.0_wp
      n = 0
      do j = NG + 1, nj - NG
         do i = NG + 1, ni - NG
            kt = 0
            do k = NZ, 1, -1
               if (ms%h_layer(i, j, k) > H_VANISHED) then
                  kt = k
                  exit
               end if
            end do
            if (kt < 4) cycle
            if (ms%h_layer(i, j, kt - 3) <= H_VANISHED) cycle
            tt = ms%tracers(ms%idx_temperature)%hTr(i, j, kt)/ms%h_layer(i, j, kt)
            tb = ms%tracers(ms%idx_temperature)%hTr(i, j, kt - 3)/ms%h_layer(i, j, kt - 3)
            s = s + (tt - tb)
            n = n + 1
         end do
      end do
      if (n > 0) s = s/real(n, wp)
   end function ml_strat

   subroutine ml_walk_check(ms, mle, ni, nj, worst_b, worst_h)
      !! Host replica of the ML walk over LIVE layers only (surface down to
      !! MLD, the straddling layer partial-weighted); the worst relative
      !! b_bar and absolute htot mismatch against the device result.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_mle_t), intent(in) :: mle
      integer, intent(in) :: ni, nj
      real(wp), intent(out) :: worst_b, worst_h
      integer :: i, j, k
      real(wp) :: htot, rint, w, b
      worst_b = 0.0_wp
      worst_h = 0.0_wp
      do j = 1, nj
         do i = 1, ni
            htot = 0.0_wp
            rint = 0.0_wp
            do k = NZ, 1, -1
               if (MLD - htot <= 0.0_wp) exit
               if (ms%h_layer(i, j, k) <= H_VANISHED) cycle
               w = min(ms%h_layer(i, j, k), MLD - htot)
               htot = htot + w
               rint = rint + ms%rho_layer(i, j, k)*w
            end do
            b = -(GRAVITY/RHO0)*(rint/(htot + 1.0e-30_wp))
            worst_b = max(worst_b, abs(mle%b_ml(i, j) - b)/abs(b))
            worst_h = max(worst_h, abs(mle%htot_ml(i, j) - htot))
         end do
      end do
   end subroutine ml_walk_check

   subroutine live_t_range(ms, ni, nj, tmin, tmax)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: ni, nj
      real(wp), intent(out) :: tmin, tmax
      integer :: i, j, k
      real(wp) :: t
      tmin = huge(1.0_wp)
      tmax = -huge(1.0_wp)
      do k = 1, NZ
         do j = NG + 1, nj - NG
            do i = NG + 1, ni - NG
               if (ms%h_layer(i, j, k) <= H_VANISHED) cycle
               t = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ms%h_layer(i, j, k)
               tmin = min(tmin, t)
               tmax = max(tmax, t)
            end do
         end do
      end do
   end subroutine live_t_range

   ! ------------------------------------------------------------------
   ! Case 1: front over the staircase — closed face-layers exactly zero
   ! ------------------------------------------------------------------
   subroutine test_closed_zero(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_mle_t) :: mle
      real(wp), allocatable :: tgt(:, :, :)
      integer :: ni, nj, n_closed
      real(wp) :: worst_closed, worst_colsum, max_flux, worst_b, worst_h
      character(len=128) :: msg

      call build_case(grid, metrics, ms, epbl, mle, tgt, 0.0_wp, .false.)
      ni = grid%nx_total
      nj = grid%ny_total
      checks: block
         call check(error, count(tgt <= H_VANISHED) > 0 .and. &
                    count(metrics%open_u < 0.5_wp) > 0 .and. &
                    count(metrics%open_v < 0.5_wp) > 0, &
                    "the staircase must carry fillers and close u- AND v-face layers")
         if (allocated(error)) exit checks
         call check(error, minval(sum(tgt, dim=3)) < MLD, &
                    "some column must be shallower than the MLD (clamp exercised)")
         if (allocated(error)) exit checks

         call map_in(ms, epbl, mle)
         call mle_compute_transports(grid, metrics, mle, ms, epbl, dt_limit=DT)
         !$omp target update from(mle%uhml, mle%vhml, mle%b_ml, mle%htot_ml)
         call map_out(ms, epbl, mle)

         call census(metrics, ms, mle, ni, nj, n_closed, worst_closed, &
                     worst_colsum, max_flux)
         call check(error, n_closed > 0, "the census must see closed face-layers")
         if (allocated(error)) exit checks
         call check(error, max_flux > 1.0e3_wp, &
                    "the ML front must drive a non-trivial MLE transport")
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3,a,es10.3,a)') "max closed |uhml| = ", worst_closed, &
            " m^3/s (open max ", max_flux, ")"
         call check(error, worst_closed == 0.0_wp, &
                    "a CLOSED or filler face-layer must carry EXACTLY zero MLE "// &
                    "transport: "//trim(msg))
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3)') "max |Sum_k| / max = ", worst_colsum/max_flux
         call check(error, worst_colsum <= 1.0e-13_wp*max_flux, &
                    "every face's column-integrated MLE transport must vanish: "//trim(msg))
         if (allocated(error)) exit checks
         ! No filler entered the ML walk: b_bar / htot match a host walk
         ! over the LIVE layers only (a filler at RHO_FILLER with weight
         ! h_min would shift b_bar by ~1e-8 relative).
         call ml_walk_check(ms, mle, ni, nj, worst_b, worst_h)
         write (msg, '(a,es10.3,a,es10.3)') "rel |db_ml| = ", worst_b, &
            ", |dhtot| = ", worst_h
         call check(error, worst_b <= 1.0e-13_wp .and. worst_h <= 1.0e-9_wp, &
                    "the ML walk must skip fillers: "//trim(msg))
      end block checks
      call teardown(metrics, ms, epbl, mle)
   end subroutine test_closed_zero

   ! ------------------------------------------------------------------
   ! Case 2: fillers above the open column (ice draft)
   ! ------------------------------------------------------------------
   subroutine test_ice_top(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_mle_t) :: mle
      real(wp), allocatable :: tgt(:, :, :)
      integer :: ni, nj, n_closed
      real(wp) :: worst_closed, worst_colsum, max_flux
      character(len=128) :: msg

      ! 120 m draft: the top two nominal layers are fillers inside the ice,
      ! the third a partial top cell.
      call build_case(grid, metrics, ms, epbl, mle, tgt, 120.0_wp, .false.)
      ni = grid%nx_total
      nj = grid%ny_total
      checks: block
         call check(error, all(tgt(:, :, NZ) <= H_VANISHED), &
                    "the k = nz layer must be a filler under the draft")
         if (allocated(error)) exit checks
         call map_in(ms, epbl, mle)
         call mle_compute_transports(grid, metrics, mle, ms, epbl, dt_limit=DT)
         !$omp target update from(mle%uhml, mle%vhml)
         call map_out(ms, epbl, mle)

         call census(metrics, ms, mle, ni, nj, n_closed, worst_closed, &
                     worst_colsum, max_flux)
         call check(error, max_flux > 1.0e3_wp, "the front under the ice must carry MLE")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(mle%uhml(:, :, NZ - 1:NZ))) == 0.0_wp .and. &
                    maxval(abs(mle%vhml(:, :, NZ - 1:NZ))) == 0.0_wp, &
                    "nothing may land in the fillers inside the ice")
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3,a,es10.3,a)') "max closed |uhml| = ", worst_closed, &
            " m^3/s (open max ", max_flux, ")"
         call check(error, worst_closed == 0.0_wp, &
                    "no MLE transport on any closed / filler face-layer: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, worst_colsum <= 1.0e-13_wp*max_flux, &
                    "column-integrated MLE transport must still vanish")
      end block checks
      call teardown(metrics, ms, epbl, mle)
   end subroutine test_ice_top

   ! ------------------------------------------------------------------
   ! Case 3: restratification of the front over NSTEP steps
   ! ------------------------------------------------------------------
   subroutine test_restrat(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_mle_t) :: mle
      type(continuity_t) :: ct
      real(wp), allocatable :: tgt(:, :, :), h_fill0(:, :, :)
      integer, parameter :: NSTEP = 24
      integer :: it, ni, nj, n_closed, n_bad, n_bad_max
      real(wp) :: worst_closed, worst_colsum, max_flux, worst_i1
      real(wp) :: wc_max, cs_rel_max, mf_max, m0, m1, t0s, t1s
      real(wp) :: s0, s1, tmin0, tmax0, tmin, tmax, tmin_run, tmax_run
      logical :: finite_all, fill_h_kept
      character(len=128) :: msg

      call build_case(grid, metrics, ms, epbl, mle, tgt, 0.0_wp, .false.)
      ni = grid%nx_total
      nj = grid%ny_total
      call ct%init(grid, nz_ml=NZ)
      allocate (h_fill0, source=ms%h_layer)
      m0 = sum(ms%h_layer(NG + 1:ni - NG, NG + 1:nj - NG, :))
      t0s = sum(ms%tracers(ms%idx_temperature)%hTr(NG + 1:ni - NG, NG + 1:nj - NG, :))
      s0 = ml_strat(ms, ni, nj)
      call live_t_range(ms, ni, nj, tmin0, tmax0)
      tmin_run = tmin0
      tmax_run = tmax0

      wc_max = 0.0_wp
      cs_rel_max = 0.0_wp
      mf_max = 0.0_wp
      n_bad_max = 0
      finite_all = .true.
      fill_h_kept = .true.
      checks: block
         call map_in(ms, epbl, mle, ct)
         do it = 1, NSTEP
            call mle_compute_transports(grid, metrics, mle, ms, epbl, dt_limit=DT)
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT, mle=mle)
            call ms%enforce_vanished_content(ni, nj)
            call ms%scan_vanished_content(ni, nj, n_bad, worst_i1)
            n_bad_max = max(n_bad_max, n_bad)
            !$omp target update from(mle%uhml, mle%vhml, ms%h_layer)
            !$omp target update from(ms%tracers(ms%idx_temperature)%hTr)
            call census(metrics, ms, mle, ni, nj, n_closed, worst_closed, &
                        worst_colsum, max_flux)
            wc_max = max(wc_max, worst_closed)
            mf_max = max(mf_max, max_flux)
            if (max_flux > 0.0_wp) cs_rel_max = max(cs_rel_max, worst_colsum/max_flux)
            finite_all = finite_all .and. all(ieee_is_finite(mle%uhml)) .and. &
                         all(ieee_is_finite(mle%vhml)) .and. &
                         all(ieee_is_finite(ms%h_layer)) .and. &
                         all(ieee_is_finite(ms%tracers(ms%idx_temperature)%hTr))
            fill_h_kept = fill_h_kept .and. &
                          all(merge(ms%h_layer == h_fill0, .true., h_fill0 <= H_VANISHED))
            call live_t_range(ms, ni, nj, tmin, tmax)
            tmin_run = min(tmin_run, tmin)
            tmax_run = max(tmax_run, tmax)
            ! The density follows T for the next step's b_bar.
            call set_rho_host(ms)
            !$omp target update to(ms%rho_layer)
         end do
         call map_out(ms, epbl, mle, ct)

         call check(error, n_closed > 0, "the census must see closed face-layers")
         if (allocated(error)) exit checks
         call check(error, mf_max > 1.0e3_wp, "the front must drive MLE transport")
         if (allocated(error)) exit checks
         write (msg, '(a,es10.3,a,es10.3,a)') "max closed |flux| = ", wc_max, &
            " m^3/s (open max ", mf_max, ")"
         call check(error, wc_max == 0.0_wp, &
                    "a CLOSED or filler face-layer must carry EXACTLY zero MLE flux: "// &
                    trim(msg))
         if (allocated(error)) exit checks
         call check(error, cs_rel_max < 1.0e-13_wp, &
                    "every face's column-integrated MLE transport must be 0 to round-off")
         if (allocated(error)) exit checks
         call check(error, finite_all, "uhml/vhml, h, hTr must stay finite")
         if (allocated(error)) exit checks
         call check(error, fill_h_kept, &
                    "a filler's thickness must not change (no MLE mass into it)")
         if (allocated(error)) exit checks
         call check(error, n_bad_max == 0, "I1' must hold after every step")
         if (allocated(error)) exit checks
         m1 = sum(ms%h_layer(NG + 1:ni - NG, NG + 1:nj - NG, :))
         t1s = sum(ms%tracers(ms%idx_temperature)%hTr(NG + 1:ni - NG, NG + 1:nj - NG, :))
         call check(error, abs(m1 - m0) <= 1.0e-13_wp*m0, "MLE must conserve mass")
         if (allocated(error)) exit checks
         call check(error, abs(t1s - t0s) <= 1.0e-13_wp*abs(t0s), &
                    "MLE must conserve temperature content")
         if (allocated(error)) exit checks
         write (msg, '(a,2f12.8,a,2f12.8)') "T range ", tmin0, tmax0, " -> ", &
            tmin_run, tmax_run
         call check(error, tmin_run >= tmin0 - 1.0e-10_wp .and. &
                    tmax_run <= tmax0 + 1.0e-10_wp, &
                    "MLE + continuity must create no new T extrema: "//trim(msg))
         if (allocated(error)) exit checks
         s1 = ml_strat(ms, ni, nj)
         write (msg, '(a,es11.4,a,es11.4,a)') "ML strat ", s0, " -> ", s1, " K"
         call check(error, abs(s0) <= 1.0e-12_wp .and. s1 > 1.0e-3_wp, &
                    "the front must RESTRATIFY the ML (warm over cold): "//trim(msg))
      end block checks
      call ct%destroy()
      call teardown(metrics, ms, epbl, mle)
   end subroutine test_restrat

   ! ------------------------------------------------------------------
   ! Case 4: nothing closed => the open path IS the full-column form
   ! ------------------------------------------------------------------
   subroutine test_flat_bitid(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_epbl_t) :: epbl
      type(ocean_mle_t) :: mle
      real(wp), allocatable :: tgt(:, :, :), u_on(:, :, :), v_on(:, :, :)

      call build_case(grid, metrics, ms, epbl, mle, tgt, 0.0_wp, .true.)
      checks: block
         call check(error, count(tgt <= H_VANISHED) == 0 .and. &
                    count(metrics%open_u < 0.5_wp) == 0, &
                    "the flat case must have no filler and no closed face")
         if (allocated(error)) exit checks
         call map_in(ms, epbl, mle)
         call mle_compute_transports(grid, metrics, mle, ms, epbl, dt_limit=DT)
         !$omp target update from(mle%uhml, mle%vhml)
         allocate (u_on, source=mle%uhml)
         allocate (v_on, source=mle%vhml)
         metrics%use_closed_faces = .false.
         call mle_compute_transports(grid, metrics, mle, ms, epbl, dt_limit=DT)
         !$omp target update from(mle%uhml, mle%vhml)
         metrics%use_closed_faces = .true.
         call map_out(ms, epbl, mle)

         call check(error, maxval(abs(u_on)) > 1.0e3_wp, "the flat front must carry MLE")
         if (allocated(error)) exit checks
         call check(error, all(u_on == mle%uhml) .and. all(v_on == mle%vhml), &
                    "with nothing closed the open path must be bit-identical "// &
                    "to the full-column (knob-off) path")
      end block checks
      call teardown(metrics, ms, epbl, mle)
   end subroutine test_flat_bitid

end module test_ocean_mle_zfixed
