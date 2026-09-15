!! Unit tests for the Q6 ocean windowed-drain WENO tracer reconstruction.
!!
!! The ocean windowed horizontal tracer-advection drain
!! (`continuity_tracer_drain`, the `dt_tracer_advect_ratio > 1` path) gains
!! an alternative face scheme: instead of the CW-PPM swept parabola it can
!! reconstruct the swept-average donor concentration with the WENO5/7/9-Z
!! ladder reused verbatim from the coastal multilayer path.  These tests
!! exercise that path ON DEVICE:
!!
!!   * uniform invariance (weno5/7/9): a horizontally uniform tracer stays
!!     uniform under a genuinely divergent windowed drain — the reconstruction
!!     reproduces constants to round-off (catches a broken stencil / mirror /
!!     rung ladder).
!!   * front / peak retention: advecting a narrow blob many windows, WENO5
!!     retains a HIGHER peak than the monotone CW-PPM (less numerical
!!     diffusion) while conserving mass and staying bounded (no runaway
!!     extrema).
!!   * u<0 coverage: the uniform-invariance test is repeated with NEGATIVE
!!     transport, and a mirror-symmetry test advects the same symmetric blob
!!     rightward and leftward and demands the leftward field equal the
!!     rightward field reflected to round-off — the strongest guard against a
!!     transposition in the u<0 branch (its swapped avail ladder + reversed
!!     stencil), which every positive-only test would pass silently.
!!   * wall-axis + interior-land: nonzero MERIDIONAL transport on the
!!     non-periodic y-axis exercises weno_face_conc_y's position-ladder
!!     degradation, and a single interior land T-cell inside a live blob
!!     exercises the ppm_mirror_h reflection.  Both assert conservation +
!!     finite + positive.
!!
!! Default off (tracer_recon = ppm) is covered by the existing
!! test_ocean_tracer_advect_window suite (unchanged ⇒ byte-identical).
module test_ocean_tracer_weno
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy
   use rdb_continuity, only: continuity_t, continuity_tracer_drain
   use rdb_recon_weno, only: TRACER_RECON_PPM, TRACER_RECON_WENO5, &
                             TRACER_RECON_WENO7, TRACER_RECON_WENO9
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_tracer_weno_tests

   integer, parameter :: NGHOST = 5   !! covers weno5 (>=3), weno7 (>=4), weno9 (>=5)
   integer, parameter :: NZ = 2
   integer, parameter :: NXP = 32, NYP = 6
   real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp

contains

   subroutine collect_ocean_tracer_weno_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("weno5_uniform_preserved", test_uniform_weno5), &
                  new_unittest("weno7_uniform_preserved", test_uniform_weno7), &
                  new_unittest("weno9_uniform_preserved", test_uniform_weno9), &
                  new_unittest("weno5_sharper_than_ppm_conservative", test_front_sharpening), &
                  new_unittest("weno5_uniform_preserved_leftward", test_uniform_weno5_leftward), &
                  new_unittest("weno5_mirror_symmetry_u_lt_0", test_leftward_mirror), &
                  new_unittest("weno_wall_axis_meridional", test_wall_axis_meridional), &
                  new_unittest("weno_interior_land_mirror", test_interior_land_mirror) &
                  ]
   end subroutine collect_ocean_tracer_weno_tests

   subroutine build_channel(grid, metrics, ms, ct, bc, nx, ny, recon)
      !! Periodic-x channel, uniform thickness H0, mapped state.  Sets the
      !! drain face-reconstruction scheme to `recon`.  hTr is left for the
      !! caller to seed on the host before map-in.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(continuity_t), intent(out) :: ct
      type(ocean_bc_state_t), intent(out) :: bc
      integer, intent(out) :: nx, ny
      integer, intent(in) :: recon
      call grid%init(NXP, NYP, NGHOST, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
      call ct%init(grid, nz_ml=NZ)
      ct%tracer_recon = recon
      call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
      bc%periodic_x = .true.
      bc%periodic_y = .false.
      bc%north_fold = .false.
      ms%h_layer = H0
   end subroutine build_channel

   subroutine teardown(ms, ct, metrics, bc)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_bc_state_t), intent(inout) :: bc
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
   end subroutine teardown

   subroutine one_window(grid, metrics, ct, ms, ratio, bc)
      !! Map ms+ct+bc onto the device, run one drain (metrics is already
      !! device-resident from make_cartesian_metrics), unmap.  hTr round-trips
      !! host<->device so the caller can re-seed uhtr between windows on host.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: ratio
      type(ocean_bc_state_t), intent(inout) :: bc
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
      !$acc enter data copyin(bc)
      call continuity_tracer_drain(grid, metrics, ct, ms, ratio, bc=bc)
      !$acc exit data delete(bc)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine one_window

   ! ---------------------------------------------------------------------
   ! Uniform-invariance: WENO must reproduce a constant to round-off.
   ! ---------------------------------------------------------------------

   subroutine uniform_invariance(error, recon, usign)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon
      real(wp), intent(in) :: usign
         !! Sign of the accumulated zonal transport: +1 exercises the drain
         !! u>0 branch, -1 the u<0 branch (its own swapped avail_up/avail_down
         !! ladder + reversed stencil gather).  |uhtr| is identical either way.
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: C0 = 34.7_wp, U0 = 0.3_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, nxp
      real(wp) :: cmax_err, conc
      call build_channel(grid, metrics, ms, ct, bc, nx, ny, recon)
      i0 = NGHOST + 1
      i1 = NGHOST + grid%nx_phys
      j0 = NGHOST + 1
      j1 = NGHOST + grid%ny_phys
      nxp = grid%nx_phys
      ! NON-uniform, exactly periodic accumulated transport (div /= 0): a
      ! constant scheme must still hold C0 (the drain reconstructs hprev =
      ! h_end + div(uhtr)/area, so Tr = hTr/hprev = C0 ⇒ flat reconstruction).
      ! `usign` flips the whole field: (1.5 + sin) > 0 always, so usign = -1
      ! makes uhtr < 0 at every face ⇒ the drain takes the u<0 branch.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ct%uhtr(i, j, k) = usign*U0*H0*DY*DT_STEP* &
                                  (1.5_wp + sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)))
            end do
         end do
      end do
      ct%vhtr = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = H0 - (ct%uhtr(i + 1, j, k) - ct%uhtr(i, j, k))/(DX*DY)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = C0*H0
            end do
         end do
      end do

      call one_window(grid, metrics, ct, ms, 1, bc)

      cmax_err = 0.0_wp
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1
               conc = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/ms%h_layer(i, j, k)
               cmax_err = max(cmax_err, abs(conc - C0))
            end do
         end do
      end do
      call check(error, cmax_err < 1.0e-12_wp*C0, &
                 "WENO drain created a spurious gradient from a uniform tracer")
      call teardown(ms, ct, metrics, bc)
   end subroutine uniform_invariance

   subroutine test_uniform_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call uniform_invariance(error, TRACER_RECON_WENO5, 1.0_wp)
   end subroutine test_uniform_weno5

   subroutine test_uniform_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call uniform_invariance(error, TRACER_RECON_WENO7, 1.0_wp)
   end subroutine test_uniform_weno7

   subroutine test_uniform_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call uniform_invariance(error, TRACER_RECON_WENO9, 1.0_wp)
   end subroutine test_uniform_weno9

   subroutine test_uniform_weno5_leftward(error)
      !! u<0 uniform invariance: NEGATIVE accumulated transport everywhere
      !! drives the drain u<0 branch (own swapped avail ladder + reversed
      !! stencil gather).  A constant tracer must stay constant to round-off.
      type(error_type), allocatable, intent(out) :: error
      call uniform_invariance(error, TRACER_RECON_WENO5, -1.0_wp)
   end subroutine test_uniform_weno5_leftward

   ! ---------------------------------------------------------------------
   ! Front / peak retention: WENO5 less diffusive than monotone CW-PPM.
   ! ---------------------------------------------------------------------

   subroutine advect_blob(grid, metrics, ms, ct, bc, nx, ny, recon, &
                          peak, tv_ratio, mass_rel, cmin)
      !! Advect a narrow Gaussian salinity blob N_WIN windows of uniform
      !! rightward transport (div = 0 ⇒ hprev = h_end = H0, pure translation +
      !! numerical diffusion).  Returns the final interior peak concentration,
      !! the retained total-variation fraction, the relative mass change, and
      !! the min concentration.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(continuity_t), intent(out) :: ct
      type(ocean_bc_state_t), intent(out) :: bc
      integer, intent(out) :: nx, ny
      integer, intent(in) :: recon
      real(wp), intent(out) :: peak, tv_ratio, mass_rel, cmin
      integer, parameter :: N_WIN = 12, RATIO = 3
      real(wp), parameter :: U0 = 0.2_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: AMP = 4.0_wp, BG = 1.0_wp, RSUP = 8.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: i0, i1, j0, j1, i, j, k, w
      real(wp) :: xc, uvol, mass0, mass1, conc, tv, tv0, area, r
      call build_channel(grid, metrics, ms, ct, bc, nx, ny, recon)
      i0 = NGHOST + 1
      i1 = NGHOST + grid%nx_phys
      j0 = NGHOST + 1
      j1 = NGHOST + grid%ny_phys
      area = DX*DY
      xc = real(NGHOST + grid%nx_phys/2, wp)
      ! Compact raised-cosine bump: EXACTLY BG (constant) outside [xc±RSUP],
      ! so the periodic-seam faces always reconstruct the constant background
      ! ⇒ the interior flux telescopes exactly and conservation is a clean
      ! machine-precision gate even though the bump translates.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               r = real(i, wp) - xc
               if (abs(r) < RSUP) then
                  conc = BG + AMP*0.5_wp*(1.0_wp + cos(PI*r/RSUP))
               else
                  conc = BG
               end if
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = conc*H0
            end do
         end do
      end do
      mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      ! Initial interior total variation (concentration) along i.
      tv0 = 0.0_wp
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1 - 1
               tv0 = tv0 + abs(ms%tracers(ms%idx_salinity)%hTr(i + 1, j, k) &
                               - ms%tracers(ms%idx_salinity)%hTr(i, j, k))/H0
            end do
         end do
      end do

      uvol = real(RATIO, wp)*U0*H0*DY*DT_STEP
      do w = 1, N_WIN
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp
         call one_window(grid, metrics, ct, ms, RATIO, bc)
      end do

      mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      mass_rel = abs(mass1 - mass0)/abs(mass0)
      peak = 0.0_wp
      cmin = huge(1.0_wp)
      tv = 0.0_wp
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1
               conc = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/H0
               peak = max(peak, conc)
               cmin = min(cmin, conc)
            end do
            do i = i0, i1 - 1
               tv = tv + abs(ms%tracers(ms%idx_salinity)%hTr(i + 1, j, k) &
                             - ms%tracers(ms%idx_salinity)%hTr(i, j, k))/H0
            end do
         end do
      end do
      tv_ratio = tv/tv0
   end subroutine advect_blob

   subroutine test_front_sharpening(error)
      !! WENO5 must (a) conserve interior mass, (b) stay bounded (no runaway
      !! extrema, positive), and (c) retain a HIGHER peak + MORE total
      !! variation than the monotone CW-PPM — i.e. sharpen the structure the
      !! PPM limiter diffuses away.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_p, grid_w
      type(multilayer_state_t) :: ms_p, ms_w
      type(continuity_t) :: ct_p, ct_w
      type(ocean_metrics_t) :: met_p, met_w
      type(ocean_bc_state_t) :: bc_p, bc_w
      real(wp), parameter :: AMP = 4.0_wp, BG = 1.0_wp
      integer :: nx, ny
      real(wp) :: peak_p, peak_w, tvr_p, tvr_w, mass_p, mass_w, cmin_p, cmin_w
      checks: block
         call advect_blob(grid_p, met_p, ms_p, ct_p, bc_p, nx, ny, &
                          TRACER_RECON_PPM, peak_p, tvr_p, mass_p, cmin_p)
         call teardown(ms_p, ct_p, met_p, bc_p)
         call advect_blob(grid_w, met_w, ms_w, ct_w, bc_w, nx, ny, &
                          TRACER_RECON_WENO5, peak_w, tvr_w, mass_w, cmin_w)
         call teardown(ms_w, ct_w, met_w, bc_w)

         ! (a) conservation — uniform transport ⇒ div = 0 ⇒ exact translation.
         call check(error, mass_w < 1.0e-10_wp, &
                    "WENO drain broke interior tracer conservation")
         if (allocated(error)) exit checks
         call check(error, mass_p < 1.0e-10_wp, &
                    "PPM drain broke interior tracer conservation (harness sanity)")
         if (allocated(error)) exit checks
         ! (b) bounded: positive, no large new maximum (allow mild WENO Gibbs).
         call check(error, cmin_w > 0.0_wp, "WENO drain lost positivity")
         if (allocated(error)) exit checks
         call check(error, peak_w <= (BG + AMP) + 0.05_wp*AMP, &
                    "WENO drain produced a runaway maximum")
         if (allocated(error)) exit checks
         ! (c) less diffusive than the monotone CW-PPM.
         call check(error, peak_w > peak_p + 1.0e-3_wp*AMP, &
                    "WENO drain not sharper than CW-PPM (peak retention)")
         if (allocated(error)) exit checks
         call check(error, tvr_w > tvr_p, &
                    "WENO drain retained less total variation than CW-PPM")
      end block checks
   end subroutine test_front_sharpening

   ! ---------------------------------------------------------------------
   ! u<0 branch: leftward transport + mirror symmetry (transposition guard).
   ! ---------------------------------------------------------------------

   subroutine run_blob_case(recon, usign, conc, mass_rel, cmin, has_nan)
      !! Advect a compact raised-cosine salinity blob N_WIN windows of
      !! SPATIALLY-UNIFORM transport of sign `usign` (>0 rightward via the
      !! drain u>0 branch, <0 leftward via the u<0 branch) through the WENO
      !! drain on a periodic-x channel (div = 0 ⇒ hprev = H0, pure
      !! translation).  The blob is symmetric about the face between the two
      !! central physical cells, so the exact solution is equivariant under
      !! (i -> mirror, u -> -u): a CORRECT u<0 branch reproduces the rightward
      !! field reflected.  Returns the final concentration field (Tr = hTr/H0),
      !! the relative interior-mass change, the interior min concentration, and
      !! a finiteness flag.
      integer, intent(in) :: recon
      real(wp), intent(in) :: usign
      real(wp), allocatable, intent(out) :: conc(:, :, :)
      real(wp), intent(out) :: mass_rel, cmin
      logical, intent(out) :: has_nan
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_bc_state_t) :: bc
      integer, parameter :: N_WIN = 8, RATIO = 3
      real(wp), parameter :: U0 = 0.15_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: AMP = 4.0_wp, BG = 1.0_wp, RSUP = 6.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, w
      real(wp) :: xc, uvol, mass0, mass1, r, c, val
      call build_channel(grid, metrics, ms, ct, bc, nx, ny, recon)
      i0 = NGHOST + 1
      i1 = NGHOST + grid%nx_phys
      j0 = NGHOST + 1
      j1 = NGHOST + grid%ny_phys
      ! Centre on the face between the two central physical cells:
      ! xc = i0 + (nx_phys-1)/2 (half-integer for even nx_phys) so the IC is
      ! EXACTLY symmetric under the reflection i -> 2*i0 + nx_phys - 1 - i.
      xc = real(i0, wp) + real(grid%nx_phys - 1, wp)*0.5_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               r = real(i, wp) - xc
               if (abs(r) < RSUP) then
                  c = BG + AMP*0.5_wp*(1.0_wp + cos(PI*r/RSUP))
               else
                  c = BG
               end if
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = c*H0
            end do
         end do
      end do
      mass0 = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      ! Spatially uniform transport ⇒ div = 0 ⇒ hprev = H0; usign < 0 ⇒
      ! uhtr < 0 at every face ⇒ the drain takes the u<0 branch.
      uvol = usign*real(RATIO, wp)*U0*H0*DY*DT_STEP
      do w = 1, N_WIN
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp
         call one_window(grid, metrics, ct, ms, RATIO, bc)
      end do
      mass1 = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      mass_rel = abs(mass1 - mass0)/abs(mass0)
      allocate (conc(nx, ny, NZ))
      cmin = huge(1.0_wp)
      has_nan = .false.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               val = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/H0
               conc(i, j, k) = val
               if (.not. (abs(val) <= huge(1.0_wp))) has_nan = .true.
               if (i >= i0 .and. i <= i1 .and. j >= j0 .and. j <= j1) then
                  cmin = min(cmin, val)
               end if
            end do
         end do
      end do
      call teardown(ms, ct, metrics, bc)
   end subroutine run_blob_case

   subroutine test_leftward_mirror(error)
      !! STRONG u<0 guard.  Advect the SAME symmetric blob rightward (u>0) and
      !! leftward (u<0); the leftward field must equal the rightward field
      !! reflected about the domain centre to round-off.  A transposition in
      !! the u<0 branch (swapped avail ladder / reversed stencil) breaks this
      !! equality while every positive-transport test stays green.  Also gates
      !! leftward conservation + positivity + finiteness.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: cR(:, :, :), cL(:, :, :)
      real(wp) :: massR, massL, cminR, cminL, err_mir
      logical :: nanR, nanL
      integer :: i0, i1, j0, j1, i, j, k, imir
      checks: block
         call run_blob_case(TRACER_RECON_WENO5, 1.0_wp, cR, massR, cminR, nanR)
         call run_blob_case(TRACER_RECON_WENO5, -1.0_wp, cL, massL, cminL, nanL)
         i0 = NGHOST + 1
         i1 = NGHOST + NXP
         j0 = NGHOST + 1
         j1 = NGHOST + NYP
         call check(error,.not. nanR .and. .not. nanL, "WENO blob produced NaN/Inf")
         if (allocated(error)) exit checks
         call check(error, massL < 1.0e-10_wp, "leftward WENO drain broke conservation")
         if (allocated(error)) exit checks
         call check(error, cminL > 0.0_wp, "leftward WENO drain lost positivity")
         if (allocated(error)) exit checks
         ! Reflection i -> 2*i0 + nx_phys - 1 - i (a bijection on [i0,i1] for
         ! even nx_phys): the leftward field must equal the rightward field
         ! there.  Mirrored arithmetic ⇒ bit-identical when the branch is right.
         err_mir = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  imir = 2*i0 + NXP - 1 - i
                  err_mir = max(err_mir, abs(cL(i, j, k) - cR(imir, j, k)))
               end do
            end do
         end do
         call check(error, err_mir < 1.0e-10_wp, &
                    "u<0 WENO branch is not the mirror of u>0 (transposition)")
      end block checks
      if (allocated(cR)) deallocate (cR)
      if (allocated(cL)) deallocate (cL)
   end subroutine test_leftward_mirror

   ! ---------------------------------------------------------------------
   ! Wall-axis meridional reconstruction + interior-land mirror.
   ! ---------------------------------------------------------------------

   subroutine test_wall_axis_meridional(error)
      !! Exercise weno_face_conc_y on the NON-PERIODIC (wall) y-axis with
      !! NONZERO meridional transport — the position ladder degrades toward
      !! the wall rows (the periodic-x tests never enter this branch because
      !! vhtr ≡ 0).  A smooth meridional bump is advected a few windows of
      !! uniform northward transport; the boundary faces carry no flux so the
      !! FULL-array tracer mass must be conserved, the field finite + positive.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_bc_state_t) :: bc
      integer, parameter :: N_WIN = 6, RATIO = 2
      real(wp), parameter :: V0 = 0.1_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: AMP = 3.0_wp, BG = 1.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i, j, k, w
      real(wp) :: mass0, mass1, vvol, val, cmin, mass_rel
      logical :: has_nan
      call build_channel(grid, metrics, ms, ct, bc, nx, ny, TRACER_RECON_WENO5)
      ! y is a wall axis (build_channel sets periodic_y = .false.).  Smooth
      ! meridional cosine across the WHOLE array so the ladder degradation
      ! near the wall rows is exercised without needing a y seam wrap.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  (BG + AMP*0.5_wp*(1.0_wp + cos(2.0_wp*PI*real(j - 1, wp)/real(ny, wp))))*H0
            end do
         end do
      end do
      mass0 = sum(ms%tracers(ms%idx_salinity)%hTr)
      vvol = real(RATIO, wp)*V0*H0*DX*DT_STEP
      do w = 1, N_WIN
         ct%uhtr = 0.0_wp
         ct%vhtr = vvol   ! uniform northward; drain zeros the outer faces
         call one_window(grid, metrics, ct, ms, RATIO, bc)
      end do
      mass1 = sum(ms%tracers(ms%idx_salinity)%hTr)
      mass_rel = abs(mass1 - mass0)/abs(mass0)
      cmin = huge(1.0_wp)
      has_nan = .false.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               val = ms%tracers(ms%idx_salinity)%hTr(i, j, k)
               if (.not. (abs(val) <= huge(1.0_wp))) has_nan = .true.
               cmin = min(cmin, val)
            end do
         end do
      end do
      call teardown(ms, ct, metrics, bc)
      checks: block
         call check(error,.not. has_nan, "wall-axis WENO drain produced NaN/Inf")
         if (allocated(error)) exit checks
         call check(error, mass_rel < 1.0e-10_wp, &
                    "wall-axis WENO drain broke total tracer conservation")
         if (allocated(error)) exit checks
         call check(error, cmin >= 0.0_wp, "wall-axis WENO drain lost positivity")
      end block checks
   end subroutine test_wall_axis_meridional

   subroutine test_interior_land_mirror(error)
      !! Exercise the interior-land reflection (ppm_mirror_h with w_nbr = 0)
      !! INSIDE a live nonzero-transport WENO stencil.  A single interior
      !! T-cell is flagged land (wet_T = 0) at the centre of an advected blob;
      !! the reconstruction must mirror the local value across it.  Only wet_T
      !! (not the area metrics) is masked, so the swept-flux divergence still
      !! telescopes and total mass stays conserved; assert conservation +
      !! finite + positive.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_bc_state_t) :: bc
      integer, parameter :: N_WIN = 5, RATIO = 2
      real(wp), parameter :: U0 = 0.15_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: AMP = 3.0_wp, BG = 1.0_wp, RSUP = 6.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, w, iland, jland
      real(wp) :: xc, uvol, mass0, mass1, mass_rel, r, c, val, cmin
      logical :: has_nan
      call build_channel(grid, metrics, ms, ct, bc, nx, ny, TRACER_RECON_WENO5)
      i0 = NGHOST + 1
      i1 = NGHOST + grid%nx_phys
      j0 = NGHOST + 1
      j1 = NGHOST + grid%ny_phys
      ! Flag one interior physical T-cell as land and push it to the device
      ! (the WENO reconstruction reads metrics%wet_T for the mirror).
      iland = i0 + grid%nx_phys/2
      jland = j0 + grid%ny_phys/2
      metrics%wet_T(iland, jland) = 0.0_wp
      !$acc update device(metrics%wet_T)
      xc = real(iland, wp)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               r = real(i, wp) - xc
               if (abs(r) < RSUP) then
                  c = BG + AMP*0.5_wp*(1.0_wp + cos(PI*r/RSUP))
               else
                  c = BG
               end if
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = c*H0
            end do
         end do
      end do
      mass0 = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      uvol = real(RATIO, wp)*U0*H0*DY*DT_STEP
      do w = 1, N_WIN
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp
         call one_window(grid, metrics, ct, ms, RATIO, bc)
      end do
      mass1 = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
      mass_rel = abs(mass1 - mass0)/abs(mass0)
      cmin = huge(1.0_wp)
      has_nan = .false.
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1
               val = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/H0
               if (.not. (abs(val) <= huge(1.0_wp))) has_nan = .true.
               cmin = min(cmin, val)
            end do
         end do
      end do
      call teardown(ms, ct, metrics, bc)
      checks: block
         call check(error,.not. has_nan, "interior-land WENO drain produced NaN/Inf")
         if (allocated(error)) exit checks
         call check(error, mass_rel < 1.0e-10_wp, &
                    "interior-land WENO drain broke tracer conservation")
         if (allocated(error)) exit checks
         call check(error, cmin > 0.0_wp, "interior-land WENO drain lost positivity")
      end block checks
   end subroutine test_interior_land_mirror

end module test_ocean_tracer_weno
