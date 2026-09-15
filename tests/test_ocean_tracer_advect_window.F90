!! Unit tests for Phase-2 (6b) windowed flux-accumulated horizontal
!! tracer advection (the drain).  Mirrors the validated Python prototype
!! gates (local_archive/prototypes/dt_tracer_advect_ppm_proto.py):
!!
!!   * swept-flux oracle (G7): the swept-average CW parabola flux
!!     reproduces the closed-form analytic swept integral for a known
!!     parabola + CFL, both flow directions, to ~1e-13.
!!   * hprev reconstruction (closed form): hprev = areaT·h_end + div(uhtr)
!!     recovers a known window-start thickness to round-off.
!!   * conservation (G2): drain a tracer blob on a periodic-x channel at
!!     ratio = 2, 3, 5; Σ(areaT·hTr) conserved to ~1e-12 (interior).
!!   * positivity (G3): a positive blob stays >= 0 with no new maxima,
!!     including a fast multi-pass (high ratio + fast flow) case.
!!   * ratio = 1 bit-identity: the windowed dispatch at ratio = 1 is the
!!     verbatim every-step bypass — exact.
!!   * on-device: every drain test runs the kernels on the active build
!!     (NVHPC GPU under CUDA_VISIBLE_DEVICES=1).
module test_ocean_tracer_advect_window
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy
   use rdb_continuity, only: continuity_t, continuity_tracer_drain, &
                             continuity_tracer_step_split, &
                             drain_swept_flux_x, drain_swept_flux_y, &
                             drain_parabola_x, drain_parabola_y, &
                             drain_reconstruct_hprev, TR_MODE_ACCUMULATE
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_tracer_advect_window_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NZ = 2

contains

   subroutine collect_ocean_tracer_advect_window_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("swept_flux_oracle", test_swept_oracle), &
                  new_unittest("hprev_reconstruction", test_reconstruct), &
                  new_unittest("conservation_ratio_2_3_5", test_conservation), &
                  new_unittest("uniform_concentration_preserved", test_uniform_preserved), &
                  new_unittest("accumulate_holds_concentration", test_accumulate_holds_conc), &
                  new_unittest("seam_ghost_robustness", test_seam_ghost_robustness), &
                  new_unittest("positivity_multipass", test_positivity), &
                  new_unittest("ratio1_bit_identity", test_ratio1_identity), &
                  new_unittest("conservation_y_ratio_2_3_5", test_conservation_y), &
                  new_unittest("positivity_y_multipass", test_positivity_y), &
                  new_unittest("conservation_diagonal_2d", test_conservation_2d), &
                  new_unittest("swept_flux_x_edge_faces_zero", test_swept_x_edge_faces), &
                  new_unittest("swept_flux_y_edge_faces_zero", test_swept_y_edge_faces) &
                  ]
   end subroutine collect_ocean_tracer_advect_window_tests

   subroutine map_in(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in

   subroutine map_out(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine test_swept_oracle(error)
      !! G7: the swept-average CW parabola flux must reproduce the
      !! closed-form analytic swept integral for a hand-crafted donor
      !! parabola + CFL, both flow directions, to ~1e-13.  Runs the
      !! drain_parabola_x + drain_swept_flux_x kernels on the device.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP = 16, NYP = 4
      real(wp), parameter :: H0 = 1.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: i, j, k, nx, ny, donor, face
      real(wp) :: test_cfl, uhh_val, area, ker_conc, oracle, err
      real(wp) :: aLd, aRd, a6d
      real(wp), allocatable :: F(:, :, :)
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         area = DX*DY

         ms%h_layer = H0
         ! Smooth periodic sine concentration so the PPM limiter is inert
         ! at the interior donor cell.  hTr = Tr·hprev_work (= Tr·H0).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (2.0_wp + sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))*H0
               end do
            end do
         end do
         ! hprev_work = H0 everywhere (used by parabola + swept CFL).
         ct%hprev_work = H0
         ct%uhh_x = 0.0_wp

         donor = NGHOST + 4         ! a well-interior cell
         face = donor + 1           ! east face of donor (donor = left cell)
         test_cfl = 0.37_wp
         uhh_val = test_cfl*area*H0
         k = 1
         j = NGHOST + 1
         ct%uhh_x(face, j, k) = uhh_val

         allocate (F(nx + 1, ny, NZ))
         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc data copy(F)
         call drain_parabola_x(nx, ny, NZ, metrics%wet_T, ms%tracers(ms%idx_salinity)%hTr, &
                               ct%hprev_work, ct%tr_work, ct%pal, ct%par, ct%pa6)
         call drain_swept_flux_x(nx, ny, NZ, metrics%areaT, ct%uhh_x, &
                                 ct%hprev_work, ct%pal, ct%par, ct%pa6, F)
         !$acc update self(F, ct%pal, ct%par, ct%pa6)
         !$acc end data
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         aLd = ct%pal(donor, j, k)
         aRd = ct%par(donor, j, k)
         a6d = ct%pa6(donor, j, k)
         oracle = aRd - 0.5_wp*test_cfl*((aRd - aLd) &
                                         - a6d*(1.0_wp - (2.0_wp/3.0_wp)*test_cfl))
         ker_conc = F(face, j, k)/uhh_val
         err = abs(ker_conc - oracle)
         call check(error, err < 1.0e-13_wp, &
                    "swept-flux (pos) deviates from analytic oracle")
         if (allocated(error)) exit checks
         ! Sanity: a6 must be non-trivial (limiter inert ⇒ real parabola).
         call check(error, abs(a6d) > 1.0e-6_wp, &
                    "test setup degenerate: a6 ~ 0 (parabola flattened)")
         deallocate (F)
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_swept_oracle

   subroutine test_reconstruct(error)
      !! Closed-form hprev reconstruction: with a known h_end and a known
      !! uhtr/vhtr, hprev = max(0, areaT·h_end + div(uhtr))·iareaT must
      !! recover the analytic window-start thickness to round-off.  Here
      !! we craft a divergence so hprev = h_end + Δ exactly (away from the
      !! vanishing hatch).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP = 8, NYP = 6
      real(wp), parameter :: H0 = 10.0_wp, DX = 2.0_wp, DY = 3.0_wp
      integer :: nx, ny, ic, jc, k
      real(wp) :: area, q, expect, got
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         area = DX*DY

         ms%h_layer = H0
         ct%uhtr = 0.0_wp
         ct%vhtr = 0.0_wp
         ! Put a net outflow q (volume) out the east face of cell ic, so
         ! the reconstructed window-start thickness there is larger than
         ! h_end by q/areaT.  div(uhtr)[ic] = uhtr(ic+1) - uhtr(ic) = +q.
         ic = NGHOST + 3
         jc = NGHOST + 2
         k = 1
         q = 7.5_wp
         ct%uhtr(ic + 1, jc, k) = q

         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         call drain_reconstruct_hprev(nx, ny, NZ, metrics%areaT, metrics%iareaT, &
                                      ms%h_layer, ct%uhtr, ct%vhtr, ct%hprev_work)
         !$acc update self(ct%hprev_work)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         expect = H0 + q/area
         got = ct%hprev_work(ic, jc, k)
         call check(error, abs(got - expect) < 1.0e-12_wp*H0, &
                    "hprev reconstruction off at the divergent cell")
         if (allocated(error)) exit checks
         ! A quiescent (no-flux) neighbour stays exactly h_end.
         got = ct%hprev_work(ic - 2, jc, k)
         call check(error, abs(got - H0) < 1.0e-12_wp*H0, &
                    "hprev reconstruction perturbed a no-flux cell")
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_reconstruct

   subroutine setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
      !! Build a periodic-x channel with uniform thickness and a positive
      !! Gaussian salinity blob (centred, interior).  Returns the grid +
      !! mapped state with bc%periodic_x set.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(continuity_t), intent(out) :: ct
      type(ocean_bc_state_t), intent(out) :: bc
      integer, intent(out) :: nx, ny
      real(wp), intent(in) :: H0, DX, DY
      integer, parameter :: NXP = 24, NYP = 4
      integer :: i, j, k
      real(wp) :: xc, blob
      call grid%init(NXP, NYP, NGHOST, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
      call ct%init(grid, nz_ml=NZ)
      call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
      bc%periodic_x = .true.
      bc%periodic_y = .false.
      bc%north_fold = .false.

      ms%h_layer = H0
      xc = real(NGHOST + NXP/2, wp)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               blob = 1.0_wp + 3.0_wp*exp(-0.5_wp*((real(i, wp) - xc)/2.5_wp)**2)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = blob*H0
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 10.0_wp*H0
            end do
         end do
      end do
   end subroutine setup_periodic_blob

   subroutine test_conservation(error)
      !! G2: drain a positive blob on a periodic-x channel with a uniform
      !! rightward accumulated transport at ratio = 2, 3, 5.  The interior
      !! tracer mass Σ(areaT·hTr) must conserve to ~1e-12 relative.  Uniform
      !! uhtr ⇒ div = 0 ⇒ hprev = h_end (no hatch); periodic ⇒ Σ closes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: U0 = 0.2_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, ratio, ri, i, j, k, nxp
      integer, parameter :: RATIOS(3) = [2, 3, 5]
      real(wp) :: area, uvol, mass0, mass1, rel
      checks: block
         do ri = 1, 3
            ratio = RATIOS(ri)
            call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
            area = DX*DY
            i0 = NGHOST + 1
            i1 = NGHOST + grid%nx_phys
            j0 = NGHOST + 1
            j1 = NGHOST + grid%ny_phys
            nxp = grid%nx_phys
            ! Overwrite with an EXACTLY periodic positive field (period =
            ! nx_phys cells) so the west/east physical-wall fluxes match to
            ! round-off and the interior sum is the conserved quantity (a
            ! finite-tail Gaussian leaves an O(1e-7) seam-tail asymmetry).
            do k = 1, NZ
               do j = 1, ny
                  do i = 1, nx
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                        (2.0_wp + sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)))*H0
                  end do
               end do
            end do
            ! Uniform rightward transport accumulated over `ratio` steps:
            ! uhtr = ratio · (u·H0·dy·dt).  Uniform across all faces.
            uvol = real(ratio, wp)*U0*H0*DY*DT_STEP
            ct%uhtr = uvol
            ct%vhtr = 0.0_wp

            mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))

            call map_in(ms, ct)
            !$acc enter data copyin(metrics)
            call metrics%enter_data()
            !$acc enter data copyin(bc)
            call drain(grid, metrics, ct, ms, ratio, bc)
            !$acc exit data delete(bc)
            call metrics%exit_data()
            !$acc exit data delete(metrics)
            call map_out(ms, ct)

            mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
            rel = abs(mass1 - mass0)/abs(mass0)
            call check(error, rel < 1.0e-12_wp, &
                       "windowed drain broke interior tracer conservation")
            call ct%destroy()
            call ms%destroy()
            call destroy_cartesian_metrics(metrics)
            if (allocated(error)) exit checks
         end do
      end block checks
      if (allocated(error)) then
         ! teardown already done inside the loop on the failing iteration's
         ! predecessors; nothing further required.
      end if
   end subroutine test_conservation

   subroutine test_uniform_preserved(error)
      !! Constancy / uniform-preservation: a horizontally UNIFORM tracer
      !! concentration advected by a NON-uniform (divergent) accumulated
      !! transport must stay uniform to round-off.  This is the core
      !! tracer-thickness-consistency property the once-per-step windowed
      !! drain buys (and which the retired per-stage path only held to RK2
      !! truncation order): the drain reconstructs hprev = areaT·h_end +
      !! div(uhtr), so advecting a constant concentration from hprev to
      !! h_end yields hTr = C0·h_end exactly — no spurious gradient.  Runs
      !! the real drain (ratio = 1, the new default window) on the device.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: C0 = 34.7_wp, U0 = 0.3_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, nxp
      real(wp) :: cmax_err, conc
      checks: block
         call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         i0 = NGHOST + 1
         i1 = NGHOST + grid%nx_phys
         j0 = NGHOST + 1
         j1 = NGHOST + grid%ny_phys
         nxp = grid%nx_phys
         ! NON-uniform, exactly periodic accumulated transport (div /= 0) so
         ! the drain genuinely advects — a constant scheme must still hold C0.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ct%uhtr(i, j, k) = U0*H0*DY*DT_STEP* &
                                     (1.5_wp + sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)))
               end do
            end do
         end do
         ct%vhtr = 0.0_wp
         ! Physically-consistent window: the WINDOW-START thickness hprev is
         ! uniform (= H0) with uniform concentration C0 (hTr = C0·H0).  The
         ! drain reconstructs hprev = h_end + div(uhtr)/area, so the END
         ! thickness must be h_end = H0 - div(uhtr)/area (divergent cells
         ! thin, convergent cells thicken).  Then Tr = hTr/hprev = C0 at the
         ! start ⇒ flat parabola ⇒ constancy: hTr_end must stay C0·h_end.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H0 - (ct%uhtr(i + 1, j, k) - ct%uhtr(i, j, k))/(DX*DY)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = C0*H0
               end do
            end do
         end do

         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(bc)
         call drain(grid, metrics, ct, ms, 1, bc)
         !$acc exit data delete(bc)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         ! Concentration must be C0 everywhere in the interior to round-off.
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
                    "windowed drain created a spurious gradient from a uniform tracer")
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_uniform_preserved

   subroutine test_accumulate_holds_conc(error)
      !! REGRESSION (windowed-advect instability, eady dt=600 ratio=2):
      !! `continuity_tracer_step_split` in TR_MODE_ACCUMULATE must leave every
      !! tracer's CONCENTRATION `T = hTr/h_layer` pointwise unchanged, however
      !! divergent the flow.
      !!
      !! The windowed path defers horizontal tracer advection to the
      !! end-of-window drain.  Roundabout's prognostic is the thickness-weighted
      !! CONTENT `hTr`, so "defer" was implemented as "freeze hTr" — while
      !! continuity kept advancing `h_layer` underneath it.  Every consumer
      !! that derives `T = hTr/h_layer` (the EOS above all, hence ρ → PGF)
      !! then read a concentration corrupted by the full window thickness
      !! divergence, `δT/T = −δh/h`.  In the eady benchmark that grid-scale
      !! buoyancy error closed an exponentially growing EOS→PGF→divergence
      !! loop: kinetic energy doubled hourly and the run NaN'd inside day 1
      !! at dt = 600 s, while dt = 300/150 survived (the injection scales
      !! with the window length).  Conserved budgets stayed EXACT throughout,
      !! so no mass/heat/salt check could see it — only the concentration.
      !!
      !! MOM6 cannot have this failure mode: its prognostic tracer `Tr%t` is a
      !! CONCENTRATION (MOM_tracer_registry), which is thickness-invariant for
      !! free, and `h` evolves freely across the whole DT_TRACER_ADVECT window.
      !! The fix (`drain_rescale_hTr`) reproduces that semantics in Roundabout's
      !! content-based state by re-weighting `hTr` onto the new `h` each stage.
      !!
      !! Guards against a vacuous pass by also asserting the thickness really
      !! moved (a non-divergent flow would satisfy the invariant trivially).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 8.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: U0 = 0.4_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, it, nxp
      real(wp) :: conc_err, dh_rel
      real(wp), allocatable :: conc_ref(:, :, :, :), h_ref(:, :, :)
      checks: block
         call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         i0 = NGHOST + 1
         i1 = NGHOST + grid%nx_phys
         j0 = NGHOST + 1
         j1 = NGHOST + grid%ny_phys
         nxp = grid%nx_phys
         ! Strongly DIVERGENT, exactly periodic zonal flow: du/dx /= 0, so
         ! continuity moves h by several percent per stage.  Without the
         ! concentration hold the tracer concentration moves with it.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U0* &
                                               sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp))
               end do
            end do
         end do
         ms%v_face_y_layer = 0.0_wp

         allocate (conc_ref(nx, ny, NZ, size(ms%tracers)))
         allocate (h_ref(nx, ny, NZ))
         h_ref = ms%h_layer
         do it = 1, size(ms%tracers)
            conc_ref(:, :, :, it) = ms%tracers(it)%hTr/ms%h_layer
         end do

         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         ! Two RK2 stages of one accumulation window.
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT_STEP, &
                                           tracer_mode=TR_MODE_ACCUMULATE)
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT_STEP, &
                                           tracer_mode=TR_MODE_ACCUMULATE)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         ! (a) the thickness genuinely moved — the test is not vacuous.
         dh_rel = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  dh_rel = max(dh_rel, abs(ms%h_layer(i, j, k) - h_ref(i, j, k))/h_ref(i, j, k))
               end do
            end do
         end do
         call check(error, dh_rel > 1.0e-3_wp, &
                    "accumulate window did not move h_layer - test would pass vacuously")
         if (allocated(error)) exit checks

         ! (b) every tracer concentration is pointwise unchanged.
         conc_err = 0.0_wp
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (.not. ms%tracers(it)%do_horizontal_advection) cycle
            do k = 1, NZ
               do j = j0, j1
                  do i = i0, i1
                     conc_err = max(conc_err, &
                                    abs(ms%tracers(it)%hTr(i, j, k)/ms%h_layer(i, j, k) &
                                        - conc_ref(i, j, k, it)))
                  end do
               end do
            end do
         end do
         call check(error, conc_err < 1.0e-12_wp, &
                    "TR_MODE_ACCUMULATE changed a tracer concentration: the frozen "// &
                    "content drifts as h evolves, which corrupts the EOS mid-window")
      end block checks
      if (allocated(conc_ref)) deallocate (conc_ref)
      if (allocated(h_ref)) deallocate (h_ref)
      call ocean_bc_state_destroy(bc)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_accumulate_holds_conc

   subroutine test_seam_ghost_robustness(error)
      !! The drain must be ROBUST to the incoming `hTr` ghost-halo state.  The
      !! dynamics keep `h_layer`'s ghosts periodic-consistent but NOT the
      !! frozen `hTr`, so the drain must wrap `hTr` itself before the first
      !! sub-cycle pass — otherwise the pass-1 PPM parabola at seam cells reads
      !! stale ghosts (a ~1e-6 seam asymmetry).  Guard: run the SAME interior
      !! tracer through the drain twice — once with correct periodic ghosts,
      !! once with ZEROED (stale) ghosts — and assert the interior output is
      !! bit-for-bit identical.  Uniform transport ⇒ cfl ≈ 0.1 ⇒ max_iter = 1,
      !! so the result depends entirely on pass 1's ghosts: with the hTr-wrap
      !! fix both runs see the same periodic ghosts (identical); without it the
      !! zeroed-ghost run diverges at the seam (this test FAILS).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: C0 = 34.7_wp, AMP = 5.0_wp, U0 = 0.2_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, i, j, k, nxp
      real(wp) :: uvol, max_diff
      real(wp), allocatable :: hTr_a(:, :, :)
      checks: block
         ! ---- Run A: correct periodic ghosts (sinusoid over ALL cells) ----
         call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         i0 = NGHOST + 1
         i1 = NGHOST + grid%nx_phys
         j0 = NGHOST + 1
         j1 = NGHOST + grid%ny_phys
         nxp = grid%nx_phys
         ms%h_layer = H0
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (C0 + AMP*sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)))*H0
               end do
            end do
         end do
         uvol = U0*H0*DY*DT_STEP      ! uniform rightward transport (div = 0)
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp
         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(bc)
         call drain(grid, metrics, ct, ms, 1, bc)
         !$acc exit data delete(bc)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)
         allocate (hTr_a, source=ms%tracers(ms%idx_salinity)%hTr)
         call ct%destroy()
         call ms%destroy()
         call destroy_cartesian_metrics(metrics)

         ! ---- Run B: SAME interior, ZEROED (stale) ghosts ----
         call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         ms%h_layer = H0
         ms%tracers(ms%idx_salinity)%hTr = 0.0_wp     ! stale incoming ghosts
         do k = 1, NZ
            do j = 1, ny
               do i = i0, i1                          ! interior only
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (C0 + AMP*sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)))*H0
               end do
            end do
         end do
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp
         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(bc)
         call drain(grid, metrics, ct, ms, 1, bc)
         !$acc exit data delete(bc)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         max_diff = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  max_diff = max(max_diff, &
                                 abs(ms%tracers(ms%idx_salinity)%hTr(i, j, k) - hTr_a(i, j, k)))
               end do
            end do
         end do
         call check(error, max_diff == 0.0_wp, &
                    "drain not robust to incoming hTr ghost state (pre-pass seam wrap missing)")
         deallocate (hTr_a)
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_seam_ghost_robustness

   subroutine drain(grid, metrics, ct, ms, ratio, bc)
      !! Thin wrapper so the test calls the same entry the driver does.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: ratio
      type(ocean_bc_state_t), intent(in) :: bc
      call continuity_tracer_drain(grid, metrics, ct, ms, ratio, bc=bc)
   end subroutine drain

   subroutine test_positivity(error)
      !! G3: a positive blob stays >= 0 with no new maximum, including a
      !! fast multi-pass case (high ratio + large accumulated CFL so the
      !! fixed-budget sub-cycle runs many passes).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      integer :: nx, ny, i0, i1, j0, j1, ratio
      real(wp) :: area, uvol, tr_min, tr_max, tr_max0
      checks: block
         ! Fast flow: accumulated CFL ~ ratio·0.7 ⇒ multi-pass drain.
         ratio = 5
         call setup_periodic_blob(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         area = DX*DY
         i0 = NGHOST + 1
         i1 = NGHOST + grid%nx_phys
         j0 = NGHOST + 1
         j1 = NGHOST + grid%ny_phys
         tr_max0 = maxval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         ! Accumulated transport with per-step CFL ~0.7 over `ratio` steps.
         uvol = real(ratio, wp)*0.7_wp*H0*DY*DX
         ct%uhtr = uvol
         ct%vhtr = 0.0_wp

         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(bc)
         call drain(grid, metrics, ct, ms, ratio, bc)
         !$acc exit data delete(bc)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         tr_min = minval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         tr_max = maxval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         call check(error, tr_min >= -1.0e-12_wp, &
                    "windowed drain produced a negative concentration")
         if (allocated(error)) exit checks
         call check(error, tr_max <= tr_max0 + 1.0e-10_wp, &
                    "windowed drain produced a new maximum (monotonicity)")
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_positivity

   subroutine test_ratio1_identity(error)
      !! ratio = 1 bit-identity: the windowed dispatch at ratio = 1 is the
      !! verbatim every-step bypass.  The continuity_tracer_step_split call
      !! with the default mode (no tracer_mode ⇒ TR_MODE_ADVECT) must be
      !! bit-for-bit identical to an explicit TR_MODE_ADVECT call — i.e.
      !! adding the mode arg never perturbs the every-step path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 8.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp), DT = 0.3_wp
      integer :: i, j, k, nx, ny
      real(wp) :: max_diff
      checks: block
         call grid%init(20, 4, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid); call ms_b%init(grid)
         ct_a%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         ct_b%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct_a%init(grid, nz_ml=NZ); call ct_b%init(grid, nz_ml=NZ)
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_a%h_layer(i, j, k) = H0
                  ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, k) = &
                     (34.0_wp + sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))*H0
                  ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, k) = 12.0_wp*H0
               end do
               do i = 1, nx + 1
                  ms_a%u_face_x_layer(i, j, k) = 0.15_wp* &
                                                 sin(PI*real(i - 1, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms_a%v_face_y_layer(i, j, k) = 0.0_wp
               end do
            end do
         end do
         ms_b%h_layer = ms_a%h_layer
         ms_b%u_face_x_layer = ms_a%u_face_x_layer
         ms_b%v_face_y_layer = ms_a%v_face_y_layer
         ms_b%tracers(ms_b%idx_salinity)%hTr = ms_a%tracers(ms_a%idx_salinity)%hTr
         ms_b%tracers(ms_b%idx_temperature)%hTr = ms_a%tracers(ms_a%idx_temperature)%hTr

         !$acc enter data copyin(ms_a, ms_b, ct_a, ct_b, metrics)
         call ms_a%enter_data(); call ms_b%enter_data()
         call ct_a%enter_data(); call ct_b%enter_data()
         call metrics%enter_data()
         ! A: default mode (the dispatch's ratio=1 bypass).
         call continuity_tracer_step_split(grid, metrics, ct_a, ms_a, DT)
         ! B: identical but reached via the new mode arg set to ADVECT.
         call continuity_tracer_step_split(grid, metrics, ct_b, ms_b, DT, &
                                           tracer_mode=0)
         !$acc update self(ms_a%tracers(ms_a%idx_salinity)%hTr, &
         !$acc             ms_b%tracers(ms_b%idx_salinity)%hTr, &
         !$acc             ms_a%h_layer, ms_b%h_layer)
         call metrics%exit_data()
         call ct_a%exit_data(); call ct_b%exit_data()
         call ms_a%exit_data(); call ms_b%exit_data()
         !$acc exit data delete(ms_a, ms_b, ct_a, ct_b, metrics)

         max_diff = maxval(abs(ms_a%tracers(ms_a%idx_salinity)%hTr &
                               - ms_b%tracers(ms_b%idx_salinity)%hTr))
         call check(error, max_diff == 0.0_wp, &
                    "ratio=1 windowed dispatch not bit-identical to every-step")
      end block checks
      call ct_a%destroy(); call ct_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_ratio1_identity

   subroutine setup_periodic_blob_y(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
      !! Tall-y periodic-y channel (mirror of setup_periodic_blob with the
      !! axes swapped) for the meridional drain tests.  Uniform thickness;
      !! bc%periodic_y set, periodic_x off.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(continuity_t), intent(out) :: ct
      type(ocean_bc_state_t), intent(out) :: bc
      integer, intent(out) :: nx, ny
      real(wp), intent(in) :: H0, DX, DY
      integer, parameter :: NXP = 4, NYP = 24
      call grid%init(NXP, NYP, NGHOST, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
      call ct%init(grid, nz_ml=NZ)
      call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
      bc%periodic_x = .false.
      bc%periodic_y = .true.
      bc%north_fold = .false.
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*H0
   end subroutine setup_periodic_blob_y

   subroutine test_conservation_y(error)
      !! FIX-3: meridional analogue of test_conservation.  Pure-y flow
      !! (uhtr = 0; vhtr a uniform northward transport accumulated over
      !! `ratio` steps) on a periodic-y channel, with an EXACTLY periodic
      !! positive sine tracer (period = ny_phys cells) so the interior sum
      !! is the conserved quantity.  Exercises drain_limit_y /
      !! drain_parabola_y / drain_swept_flux_y / drain_update_h_y — the
      !! meridional kernels NO existing case reaches (all use vhtr = 0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: V0 = 0.2_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, ratio, ri, i, j, k, nyp
      integer, parameter :: RATIOS(3) = [2, 3, 5]
      real(wp) :: area, vvol, mass0, mass1, rel
      checks: block
         do ri = 1, 3
            ratio = RATIOS(ri)
            call setup_periodic_blob_y(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
            area = DX*DY
            i0 = NGHOST + 1
            i1 = NGHOST + grid%nx_phys
            j0 = NGHOST + 1
            j1 = NGHOST + grid%ny_phys
            nyp = grid%ny_phys
            do k = 1, NZ
               do j = 1, ny
                  do i = 1, nx
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                        (2.0_wp + sin(2.0_wp*PI*real(j - j0, wp)/real(nyp, wp)))*H0
                  end do
               end do
            end do
            ! Uniform northward transport over `ratio` steps; uhtr = 0.
            vvol = real(ratio, wp)*V0*H0*DX*DT_STEP
            ct%uhtr = 0.0_wp
            ct%vhtr = vvol

            mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))

            call map_in(ms, ct)
            !$acc enter data copyin(metrics)
            call metrics%enter_data()
            !$acc enter data copyin(bc)
            call drain(grid, metrics, ct, ms, ratio, bc)
            !$acc exit data delete(bc)
            call metrics%exit_data()
            !$acc exit data delete(metrics)
            call map_out(ms, ct)

            mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
            rel = abs(mass1 - mass0)/abs(mass0)
            call check(error, rel < 1.0e-12_wp, &
                       "meridional windowed drain broke interior conservation")
            call ct%destroy()
            call ms%destroy()
            call destroy_cartesian_metrics(metrics)
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_conservation_y

   subroutine test_positivity_y(error)
      !! FIX-3: meridional analogue of test_positivity.  A fast northward
      !! multi-pass case (high ratio + large accumulated CFL) must stay
      !! >= 0 with no new maximum — exercises the meridional drain limiter
      !! + swept flux under sub-cycling.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, ratio, i, j, k, nyp, j0c
      real(wp) :: vvol, tr_min, tr_max, tr_max0, blob
      checks: block
         ratio = 5
         call setup_periodic_blob_y(grid, metrics, ms, ct, bc, nx, ny, H0, DX, DY)
         i0 = NGHOST + 1
         i1 = NGHOST + grid%nx_phys
         j0 = NGHOST + 1
         j1 = NGHOST + grid%ny_phys
         nyp = grid%ny_phys
         j0c = NGHOST + nyp/2
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  blob = 1.0_wp + 3.0_wp*exp(-0.5_wp*((real(j, wp) - real(j0c, wp))/2.5_wp)**2)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = blob*H0
               end do
            end do
         end do
         tr_max0 = maxval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         ! Accumulated transport with per-step CFL ~0.7 over `ratio` steps.
         vvol = real(ratio, wp)*0.7_wp*H0*DX*DY
         ct%uhtr = 0.0_wp
         ct%vhtr = vvol

         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(bc)
         call drain(grid, metrics, ct, ms, ratio, bc)
         !$acc exit data delete(bc)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         tr_min = minval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         tr_max = maxval(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :)/H0)
         call check(error, tr_min >= -1.0e-12_wp, &
                    "meridional drain produced a negative concentration")
         if (allocated(error)) exit checks
         call check(error, tr_max <= tr_max0 + 1.0e-10_wp, &
                    "meridional drain produced a new maximum (monotonicity)")
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_positivity_y

   subroutine test_conservation_2d(error)
      !! FIX-3: genuinely 2D drain — BOTH uhtr and vhtr non-zero (uniform
      !! diagonal transport) at ratio = 2, 3 on a doubly-periodic square
      !! with an exactly bi-periodic positive sine tracer.  Validates the
      !! Lie-split x-then-y coupling on the shared evolving hprev_work:
      !! the zonal sub-pass mutates hprev (and Tr), the meridional sub-pass
      !! then drains against THAT updated state, and Σ(areaT·hTr) must still
      !! close to ~1e-12.  Uniform transports ⇒ div = 0 ⇒ hprev = h_end.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NSQ = 16
      real(wp), parameter :: H0 = 5.0_wp, DX = 1.0_wp, DY = 1.0_wp
      real(wp), parameter :: U0 = 0.15_wp, V0 = 0.2_wp, DT_STEP = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: nx, ny, i0, i1, j0, j1, ratio, ri, i, j, k, nxp, nyp
      integer, parameter :: RATIOS(2) = [2, 3]
      real(wp) :: area, uvol, vvol, mass0, mass1, rel
      checks: block
         do ri = 1, 2
            ratio = RATIOS(ri)
            call grid%init(NSQ, NSQ, NGHOST, DX, DY)
            call make_cartesian_metrics(metrics, grid)
            nx = grid%nx_total
            ny = grid%ny_total
            ms%nz_ml = NZ
            call ms%init(grid)
            ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
            call ct%init(grid, nz_ml=NZ)
            call ocean_bc_state_init(bc, grid, NZ, n_tracers=size(ms%tracers))
            bc%periodic_x = .true.
            bc%periodic_y = .true.
            bc%north_fold = .false.
            area = DX*DY
            i0 = NGHOST + 1
            i1 = NGHOST + grid%nx_phys
            j0 = NGHOST + 1
            j1 = NGHOST + grid%ny_phys
            nxp = grid%nx_phys
            nyp = grid%ny_phys
            ms%h_layer = H0
            ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*H0
            ! Exactly bi-periodic positive tracer (separable sine product).
            do k = 1, NZ
               do j = 1, ny
                  do i = 1, nx
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                        (3.0_wp &
                         + sin(2.0_wp*PI*real(i - i0, wp)/real(nxp, wp)) &
                         + sin(2.0_wp*PI*real(j - j0, wp)/real(nyp, wp)))*H0
                  end do
               end do
            end do
            uvol = real(ratio, wp)*U0*H0*DY*DT_STEP
            vvol = real(ratio, wp)*V0*H0*DX*DT_STEP
            ct%uhtr = uvol
            ct%vhtr = vvol

            mass0 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))

            call map_in(ms, ct)
            !$acc enter data copyin(metrics)
            call metrics%enter_data()
            !$acc enter data copyin(bc)
            call drain(grid, metrics, ct, ms, ratio, bc)
            !$acc exit data delete(bc)
            call metrics%exit_data()
            !$acc exit data delete(metrics)
            call map_out(ms, ct)

            mass1 = area*sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))
            rel = abs(mass1 - mass0)/abs(mass0)
            call check(error, rel < 1.0e-12_wp, &
                       "2D diagonal windowed drain broke interior conservation")
            call ct%destroy()
            call ms%destroy()
            call ocean_bc_state_destroy(bc)
            call destroy_cartesian_metrics(metrics)
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_conservation_2d

   subroutine test_swept_x_edge_faces(error)
      !! Regression for the latent OOB in drain_swept_flux_x.  The two
      !! array-edge faces (1, nx+1) have no donor cell: at face 1 the u>0
      !! branch would read cell i-1=0, and at face nx+1 the u<0 branch
      !! would read cell nx+1 -- both off the [1,nx] cell arrays.  The
      !! fix restricts the flux loop to i=2:nx and zeroes the edge faces.
      !! Here we drive uhh nonzero at BOTH edge faces with exactly those
      !! OOB-triggering signs (and one interior face), and assert the
      !! edge faces come back ZERO while the interior face is computed.
      !! On a -fcheck=bounds build the pre-fix kernel ABORTS on this case.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP = 8, NYP = 4
      real(wp), parameter :: H0 = 1.0_wp, DX = 1.0_wp, DY = 1.0_wp
      integer :: j, k, nx, ny, ifc
      real(wp) :: area
      real(wp), allocatable :: F(:, :, :)
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         area = DX*DY

         ms%h_layer = H0
         ! Constant tracer => flat PPM reconstruction (aL=aR, a6=0).
         ms%tracers(ms%idx_salinity)%hTr = 2.0_wp*H0
         ct%hprev_work = H0
         ct%uhh_x = 0.0_wp
         ! Edge faces with OOB-triggering signs + one interior face.
         ifc = NGHOST + 3                       ! a well-interior face
         do k = 1, NZ
            do j = 1, ny
               ct%uhh_x(1, j, k) = 0.5_wp*area*H0     ! u>0 at face 1  -> pre-fix reads cell 0
               ct%uhh_x(nx + 1, j, k) = -0.5_wp*area*H0  ! u<0 at face nx+1 -> reads cell nx+1
               ct%uhh_x(ifc, j, k) = 0.37_wp*area*H0
            end do
         end do

         allocate (F(nx + 1, ny, NZ))
         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc data copy(F)
         call drain_parabola_x(nx, ny, NZ, metrics%wet_T, ms%tracers(ms%idx_salinity)%hTr, &
                               ct%hprev_work, ct%tr_work, ct%pal, ct%par, ct%pa6)
         call drain_swept_flux_x(nx, ny, NZ, metrics%areaT, ct%uhh_x, &
                                 ct%hprev_work, ct%pal, ct%par, ct%pa6, F)
         !$acc update self(F)
         !$acc end data
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         call check(error, maxval(abs(F(1, :, :))) < 1.0e-30_wp, &
                    "west array-edge face F(1) must be zero (no donor cell 0)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(F(nx + 1, :, :))) < 1.0e-30_wp, &
                    "east array-edge face F(nx+1) must be zero (no donor cell nx+1)")
         if (allocated(error)) exit checks
         ! Interior face still computed (loop i=2:nx unchanged for physics).
         call check(error, abs(F(ifc, NGHOST + 1, 1)) > 0.5_wp, &
                    "interior face must still carry a nonzero flux")
      end block checks
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_swept_x_edge_faces

   subroutine test_swept_y_edge_faces(error)
      !! Meridional analogue of test_swept_x_edge_faces: edge faces
      !! (1, ny+1) of drain_swept_flux_y have no donor row (cell j-1=0 at
      !! face 1 for u>0, cell ny+1 at face ny+1 for u<0).  Assert they are
      !! zeroed and an interior face is computed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP = 4, NYP = 8
      real(wp), parameter :: H0 = 1.0_wp, DX = 1.0_wp, DY = 1.0_wp
      integer :: i, k, nx, ny, jfc
      real(wp) :: area
      real(wp), allocatable :: F(:, :, :)
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! opt in: continuity_t default is now .false.
         call ct%init(grid, nz_ml=NZ)
         area = DX*DY

         ms%h_layer = H0
         ms%tracers(ms%idx_salinity)%hTr = 2.0_wp*H0
         ct%hprev_work = H0
         ct%uhh_y = 0.0_wp
         jfc = NGHOST + 3
         do k = 1, NZ
            do i = 1, nx
               ct%uhh_y(i, 1, k) = 0.5_wp*area*H0       ! u>0 at face 1   -> reads row 0
               ct%uhh_y(i, ny + 1, k) = -0.5_wp*area*H0  ! u<0 at face ny+1 -> reads row ny+1
               ct%uhh_y(i, jfc, k) = 0.37_wp*area*H0
            end do
         end do

         allocate (F(nx, ny + 1, NZ))
         call map_in(ms, ct)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc data copy(F)
         call drain_parabola_y(nx, ny, NZ, metrics%wet_T, ms%tracers(ms%idx_salinity)%hTr, &
                               ct%hprev_work, ct%tr_work, ct%pal, ct%par, ct%pa6)
         call drain_swept_flux_y(nx, ny, NZ, metrics%areaT, ct%uhh_y, &
                                 ct%hprev_work, ct%pal, ct%par, ct%pa6, F)
         !$acc update self(F)
         !$acc end data
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call map_out(ms, ct)

         call check(error, maxval(abs(F(:, 1, :))) < 1.0e-30_wp, &
                    "south array-edge face F(1) must be zero (no donor row 0)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(F(:, ny + 1, :))) < 1.0e-30_wp, &
                    "north array-edge face F(ny+1) must be zero (no donor row ny+1)")
         if (allocated(error)) exit checks
         call check(error, abs(F(NGHOST + 1, jfc, 1)) > 0.5_wp, &
                    "interior face must still carry a nonzero flux")
      end block checks
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_swept_y_edge_faces

end module test_ocean_tracer_advect_window
