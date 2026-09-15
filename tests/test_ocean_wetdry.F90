!! Unit tests for ocean dynamic wetting/drying (docs/ocean_wetdry_plan.md).
!!
!! Drives `barotropic_substep_nonlinear` directly (the level the wet/dry
!! machinery lives at) with the wd_* workspaces allocated the way
!! `configure_ocean_wetdry` does.  Gates mirror the validated Python
!! prototype (`local_archive/prototypes/ocean_wetdry_proto.py`):
!!
!!   1. `wetdry_allwet_inert` — deep basin, knob ON: the limiter
!!      (theta) and the bed-blocking gate (wd_open) must be EXACTLY
!!      inert (== 1 everywhere), depth stays positive, and interior
!!      mass is conserved to round-off.  (The knob-OFF path is
!!      byte-identical by construction — it runs the untouched
!!      pre-change code, anchored by the whole existing ocean suite.)
!!
!!   2. `wetdry_thacker` — the Thacker (1981) parabolic-basin planar
!!      solution: analytic moving shoreline.  Gates (prototype numbers:
!!      period −0.89%, shoreline 3.6 cells, mass 0.0, min D = 0):
!!      period error < 2%, first-period mean shoreline error < 4 cells,
!!      min total depth >= 0, interior mass to round-off.
!!
!!   3. `wetdry_drying_slosh` — sloping-beach seiche that dries and
!!      rewets the shallow bank: positivity, round-off mass
!!      conservation, and the hysteresis mask actually toggles
!!      (dries then rewets — the machinery demonstrably engaged).
module test_ocean_wetdry
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, H_VANISHED, VCOORD_SIGMA
   use rdb_grid, only: hgrid_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_config, only: config_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg, &
                              seed_h_layer_uniform_impl
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_wetdry_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_wetdry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wetdry_allwet_inert", test_allwet_inert), &
                  new_unittest("wetdry_thacker", test_thacker), &
                  new_unittest("wetdry_drying_slosh", test_drying_slosh), &
                  new_unittest("wetdry_positivity_substep", test_positivity_substep), &
                  new_unittest("wetdry_starts_dry", test_starts_dry), &
                  new_unittest("wetdry_seed_criterion", test_seed_criterion), &
                  new_unittest("wetdry_emerged_seed_production", &
                               test_emerged_seed_production), &
                  new_unittest("wetdry_emerged_beach_multilayer_step", &
                               test_emerged_beach_multilayer_step), &
                  new_unittest("wetdry_emerged_tracer_gradient", &
                               test_emerged_tracer_gradient), &
                  new_unittest("wetdry_vanished_column_flood", &
                               test_vanished_column_flood) &
                  ]
   end subroutine collect_ocean_wetdry_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine enable_wetdry(bt_work, grid, dry_depth, rewet_depth)
      !! Mirror `configure_ocean_wetdry`: knobs + lazy wd_* allocations +
      !! hysteresis-mask seed from the current (bt_H_ref + bt_eta) depth.
      use rdb_barotropic_workstate, only: barotropic_workstate_t
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: dry_depth, rewet_depth
      integer :: nx, ny, i, j

      nx = grid%nx_total
      ny = grid%ny_total
      bt_work%wetdry_enable = .true.
      bt_work%wd_dry_depth = dry_depth
      bt_work%wd_rewet_depth = rewet_depth
      allocate (bt_work%wd_wet_dyn(nx, ny), source=1.0_wp)
      allocate (bt_work%wd_theta(nx, ny), source=1.0_wp)
      allocate (bt_work%wd_flux_x(nx + 1, ny), source=0.0_wp)
      allocate (bt_work%wd_flux_y(nx, ny + 1), source=0.0_wp)
      allocate (bt_work%wd_open_u(nx + 1, ny), source=1.0_wp)
      allocate (bt_work%wd_open_v(nx, ny + 1), source=1.0_wp)
      do j = 1, ny
         do i = 1, nx
            if (bt_work%bt_H_ref(i, j) + bt_work%bt_eta(i, j) < dry_depth) then
               bt_work%wd_wet_dyn(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine enable_wetdry

   subroutine run_chunk(grid, metrics, dyn, cor, fu, fv, n_steps, dt_inner)
      !! One substep-chunk with full device round-trip: run, pull the
      !! end-of-chunk state + the wet/dry diagnostics host-ward, restore
      !! the prognostics from the end-state (the substep exits with the
      !! TIME-MEAN in bt_eta/bt_ubt/bt_vbt), release the device copies.
      !! The persistent hysteresis mask must be pulled BEFORE exit_data
      !! (exit deletes, not copyout) or the GPU run loses the front state
      !! between chunks.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      real(wp), intent(in) :: fu(:, :), fv(:, :)
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      !$acc enter data copyin(dyn, cor, fu, fv)
      call dyn%enter_data()
      call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                        fu, fv, &
                                        n_steps, dt_inner, &
                                        bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                        eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                        bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                        uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                        vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn%bt_work%bt_eta_end, dyn%bt_work%bt_ubt_end, dyn%bt_work%bt_vbt_end)
      if (allocated(dyn%bt_work%wd_wet_dyn)) then
         !$acc update self(dyn%bt_work%wd_wet_dyn, dyn%bt_work%wd_theta)
         !$acc update self(dyn%bt_work%wd_open_u, dyn%bt_work%wd_open_v)
      end if
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv)
      ! Continue the trajectory from the end-of-chunk snapshot.
      dyn%bt_work%bt_eta = dyn%bt_work%bt_eta_end
      dyn%bt_work%bt_ubt = dyn%bt_work%bt_ubt_end
      dyn%bt_work%bt_vbt = dyn%bt_work%bt_vbt_end
   end subroutine run_chunk

   pure function interior_mass(grid, h_ref, eta) result(m)
      !! Total interior water volume (m³, uniform-area cells).  Physical
      !! cells only — ghost dynamics are wall-isolated from the interior.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h_ref(:, :), eta(:, :)
      real(wp) :: m
      integer :: i, j
      m = 0.0_wp
      do j = grid%nghost + 1, grid%nghost + grid%ny_phys
         do i = grid%nghost + 1, grid%nghost + grid%nx_phys
            m = m + (h_ref(i, j) + eta(i, j))*grid%dx*grid%dy
         end do
      end do
   end function interior_mass

   subroutine test_allwet_inert(error)
      !! Deep basin (H = 200 m), gravity-wave IC, wet/dry knob ON: the
      !! machinery must be exactly inert — theta == 1 and wd_open == 1
      !! everywhere, min depth positive, interior mass to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 200.0_wp
      real(wp), parameter :: DX = 100.0_wp
      real(wp) :: dt_inner, m0, m1, pi_l
      integer :: i, j, nx, ny

      checks: block
         call make_grid(grid, 24, 20, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total
         dyn%bt_work%bt_H_ref = H_REF
         pi_l = 4.0_wp*atan(1.0_wp)
         do j = 1, ny
            do i = 1, nx
               dyn%bt_work%bt_eta(i, j) = &
                  0.5_wp*cos(pi_l*real(i - NGHOST, wp)/real(grid%nx_phys, wp))
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, 0.05_wp, 0.10_wp)
         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)
         dt_inner = 0.3_wp*DX/sqrt(GRAVITY*H_REF)
         m0 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)

         call run_chunk(grid, metrics, dyn, cor, fu, fv, 200, dt_inner)

         m1 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         call check(error, maxval(abs(dyn%bt_work%wd_theta - 1.0_wp)) == 0.0_wp, &
                    "all-wet: limiter theta engaged on deep water (must be exactly 1)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%wd_open_u - 1.0_wp)) == 0.0_wp .and. &
                    maxval(abs(dyn%bt_work%wd_open_v - 1.0_wp)) == 0.0_wp, &
                    "all-wet: bed-blocking gate closed a deep face (must be exactly 1)")
         if (allocated(error)) exit checks
         call check(error, minval(dyn%bt_work%bt_H_ref + dyn%bt_work%bt_eta) > 0.0_wp, &
                    "all-wet: negative total depth")
         if (allocated(error)) exit checks
         call check(error, abs(m1 - m0) <= 1.0e-12_wp*abs(m0), &
                    "all-wet: interior mass not conserved to round-off")
      end block checks
      if (allocated(fu)) deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_allwet_inert

   subroutine test_thacker(error)
      !! Thacker (1981) 1-D paraboloid planar solution (prototype §(a)):
      !!   H_ref(x) = D0·(1 − x²/L²),  ω = √(2·g·D0)/L
      !!   u(t) = −B·ω·sin(ωt)  (spatially uniform ⇒ the KE-gradient and
      !!   ζ advection terms vanish on the exact solution)
      !!   shoreline x_s(t) = B·cos(ωt) ± L
      !! Gates (prototype: period −0.89%, shoreline 3.6 cells, mass 0):
      !! period < 2%, first-period mean shoreline < 4 cells, min D >= 0
      !! (to round-off), interior mass to round-off over 2 periods.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: D0 = 10.0_wp
      real(wp), parameter :: L_BASIN = 4000.0_wp
      real(wp), parameter :: B_AMP = 1500.0_wp
      real(wp), parameter :: DRY_D = 0.02_wp, REWET_D = 0.05_wp
      integer, parameter :: NXP = 160, NYP = 4
      integer, parameter :: CHUNK = 10
      real(wp) :: lx, dx, dt_inner, omega, period, x_c
      real(wp) :: m0, m1, min_d, t_now, t_prev, u_prev, u_now
      real(wp) :: xs_ana, xs_num, sl_err_sum, d_col
      real(wp), allocatable :: crossings(:)
      integer :: i, j, nx, ny, n_chunks, ic, n_cross, sl_n, i_umid
      logical :: in_p1

      checks: block
         lx = 2.9_wp*L_BASIN
         dx = lx/real(NXP, wp)
         call make_grid(grid, NXP, NYP, dx, dx)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total
         omega = sqrt(2.0_wp*GRAVITY*D0)/L_BASIN
         period = 8.0_wp*atan(1.0_wp)/omega

         ! Paraboloid bathymetry + exact planar IC (u = 0 at t = 0);
         ! banks start truly dry (D = 0 exactly — no film).
         do j = 1, ny
            do i = 1, nx
               x_c = (real(i - NGHOST, wp) - 0.5_wp)*dx - 0.5_wp*lx
               dyn%bt_work%bt_H_ref(i, j) = D0*(1.0_wp - (x_c/L_BASIN)**2)
               dyn%bt_work%bt_eta(i, j) = max( &
                                          (D0*B_AMP/L_BASIN**2)*(2.0_wp*x_c - B_AMP), &
                                          -dyn%bt_work%bt_H_ref(i, j))
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)
         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)
         dt_inner = 0.4_wp*dx/sqrt(GRAVITY*D0)
         n_chunks = int(2.0_wp*period/(dt_inner*real(CHUNK, wp))) + 1
         allocate (crossings(n_chunks), source=0.0_wp)

         m0 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         min_d = huge(1.0_wp)
         sl_err_sum = 0.0_wp
         sl_n = 0
         n_cross = 0
         t_prev = 0.0_wp
         u_prev = 0.0_wp
         i_umid = NGHOST + NXP/2  ! u-face near the basin centre
         in_p1 = .true.

         do ic = 1, n_chunks
            call run_chunk(grid, metrics, dyn, cor, fu, fv, CHUNK, dt_inner)
            t_now = real(ic*CHUNK, wp)*dt_inner
            min_d = min(min_d, minval(dyn%bt_work%bt_H_ref(NGHOST + 1:NGHOST + NXP, :) + &
                                      dyn%bt_work%bt_eta(NGHOST + 1:NGHOST + NXP, :)))
            ! upward zero crossings of u(t) at the basin centre -> period
            u_now = dyn%bt_work%bt_ubt(i_umid, NGHOST + 2)
            if (ic > 1 .and. u_prev < 0.0_wp .and. u_now >= 0.0_wp) then
               n_cross = n_cross + 1
               ! linear interpolation of the crossing instant
               crossings(n_cross) = t_prev + (t_now - t_prev)*(-u_prev)/(u_now - u_prev)
            end if
            u_prev = u_now
            t_prev = t_now
            ! first-period shoreline error (right shore, cells)
            in_p1 = t_now <= period
            if (in_p1) then
               xs_ana = B_AMP*cos(omega*t_now) + L_BASIN
               xs_num = -huge(1.0_wp)
               do i = NGHOST + NXP, NGHOST + 1, -1
                  d_col = dyn%bt_work%bt_H_ref(i, NGHOST + 2) + &
                          dyn%bt_work%bt_eta(i, NGHOST + 2)
                  if (d_col > DRY_D) then
                     xs_num = (real(i - NGHOST, wp))*dx - 0.5_wp*lx
                     exit
                  end if
               end do
               if (xs_num > -huge(1.0_wp)) then
                  sl_err_sum = sl_err_sum + abs(xs_ana - xs_num)
                  sl_n = sl_n + 1
               end if
            end if
         end do

         m1 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         call check(error, min_d >= -1.0e-9_wp, &
                    "thacker: negative total depth (positivity broken)")
         if (allocated(error)) exit checks
         call check(error, abs(m1 - m0) <= 1.0e-12_wp*abs(m0), &
                    "thacker: interior mass not conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, sl_n > 0 .and. sl_err_sum/real(max(sl_n, 1), wp) < 4.0_wp*dx, &
                    "thacker: first-period mean shoreline error >= 4 cells")
         if (allocated(error)) exit checks
         call check(error, n_cross >= 2, "thacker: too few u zero-crossings for a period fit")
         if (allocated(error)) exit checks
         call check(error, abs((crossings(n_cross) - crossings(1))/ &
                               real(n_cross - 1, wp) - period)/period < 0.02_wp, &
                    "thacker: period error >= 2%")
      end block checks
      if (allocated(fu)) deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_thacker

   subroutine test_drying_slosh(error)
      !! Sloping-beach seiche: a tilted surface sloshes over a linear
      !! beach whose top stands above still-water level, so the shallow
      !! bank dries and rewets each cycle.  Gates: positivity, interior
      !! mass to round-off, and the hysteresis mask actually toggles
      !! (some cell dries at some point AND is wet again by another) —
      !! i.e. the machinery demonstrably engaged, unlike the deep test.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: DRY_D = 0.02_wp, REWET_D = 0.05_wp
      real(wp), parameter :: H_DEEP = 8.0_wp, H_TOP = -1.0_wp
      integer, parameter :: NXP = 80, NYP = 4
      integer, parameter :: CHUNK = 20
      real(wp) :: lx, dx, dt_inner, m0, m1, min_d, frac
      integer :: i, j, nx, ny, ic, n_chunks
      logical :: saw_dry, saw_rewet
      logical, allocatable :: was_dry(:)

      checks: block
         lx = 8000.0_wp
         dx = lx/real(NXP, wp)
         call make_grid(grid, NXP, NYP, dx, dx)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ! linear beach: deep at west, bank above still water at east
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               dyn%bt_work%bt_H_ref(i, j) = H_DEEP + (H_TOP - H_DEEP)*frac
               ! tilted IC: high water at the west, dry bank at the east
               dyn%bt_work%bt_eta(i, j) = max(0.8_wp*(0.5_wp - frac), &
                                              -dyn%bt_work%bt_H_ref(i, j))
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)
         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)
         allocate (was_dry(nx), source=.false.)
         dt_inner = 0.35_wp*dx/sqrt(GRAVITY*H_DEEP)
         m0 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         min_d = huge(1.0_wp)
         saw_dry = .false.
         saw_rewet = .false.
         n_chunks = 150

         do ic = 1, n_chunks
            call run_chunk(grid, metrics, dyn, cor, fu, fv, CHUNK, dt_inner)
            min_d = min(min_d, minval(dyn%bt_work%bt_H_ref(NGHOST + 1:NGHOST + NXP, :) + &
                                      dyn%bt_work%bt_eta(NGHOST + 1:NGHOST + NXP, :)))
            do i = NGHOST + 1, NGHOST + NXP
               if (dyn%bt_work%wd_wet_dyn(i, NGHOST + 2) < 0.5_wp) then
                  ! only count cells that were WET in the IC (i.e. a
                  ! genuine drying event, not the never-wet bank top)
                  if (dyn%bt_work%bt_H_ref(i, NGHOST + 2) > REWET_D) then
                     saw_dry = .true.
                     was_dry(i) = .true.
                  end if
               else if (was_dry(i)) then
                  saw_rewet = .true.
               end if
            end do
         end do

         m1 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         call check(error, min_d >= -1.0e-9_wp, &
                    "drying slosh: negative total depth (positivity broken)")
         if (allocated(error)) exit checks
         call check(error, abs(m1 - m0) <= 1.0e-12_wp*abs(m0), &
                    "drying slosh: interior mass not conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, saw_dry, &
                    "drying slosh: no initially-wet cell ever dried (machinery inert?)")
         if (allocated(error)) exit checks
         call check(error, saw_rewet, &
                    "drying slosh: no dried cell ever rewetted (hysteresis stuck)")
      end block checks
      if (allocated(fu)) deallocate (fu, fv)
      if (allocated(was_dry)) deallocate (was_dry)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_drying_slosh

   subroutine test_positivity_substep(error)
      !! Substep-granularity positivity.  The thacker / drying_slosh
      !! tests sample `min D` once per 10-20-substep chunk, so an
      !! intra-chunk negative depth could hide.  This drives the same
      !! drying-beach geometry with CHUNK = 1 — the total depth is
      !! pulled host-ward after EVERY barotropic substep and asserted
      !! non-negative (to round-off).  Also confirms the front actually
      !! dries within the window, so the check is not vacuous.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: DRY_D = 0.02_wp, REWET_D = 0.05_wp
      real(wp), parameter :: H_DEEP = 8.0_wp, H_TOP = -1.0_wp
      integer, parameter :: NXP = 40, NYP = 4
      integer, parameter :: N_SUBSTEPS = 800
      real(wp) :: lx, dx, dt_inner, m0, m1, min_d, frac
      integer :: i, j, nx, ny, is
      logical :: saw_dry

      checks: block
         lx = 8000.0_wp
         dx = lx/real(NXP, wp)
         call make_grid(grid, NXP, NYP, dx, dx)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               dyn%bt_work%bt_H_ref(i, j) = H_DEEP + (H_TOP - H_DEEP)*frac
               dyn%bt_work%bt_eta(i, j) = max(0.8_wp*(0.5_wp - frac), &
                                              -dyn%bt_work%bt_H_ref(i, j))
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)
         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)
         dt_inner = 0.35_wp*dx/sqrt(GRAVITY*H_DEEP)
         m0 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         min_d = huge(1.0_wp)
         saw_dry = .false.

         ! CHUNK = 1: one barotropic substep per device round-trip, so
         ! `bt_eta` / `wd_wet_dyn` are inspected at substep granularity.
         do is = 1, N_SUBSTEPS
            call run_chunk(grid, metrics, dyn, cor, fu, fv, 1, dt_inner)
            min_d = min(min_d, minval(dyn%bt_work%bt_H_ref(NGHOST + 1:NGHOST + NXP, :) + &
                                      dyn%bt_work%bt_eta(NGHOST + 1:NGHOST + NXP, :)))
            do i = NGHOST + 1, NGHOST + NXP
               if (dyn%bt_work%wd_wet_dyn(i, NGHOST + 2) < 0.5_wp .and. &
                   dyn%bt_work%bt_H_ref(i, NGHOST + 2) > REWET_D) saw_dry = .true.
            end do
         end do

         m1 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         call check(error, min_d >= -1.0e-12_wp, &
                    "substep positivity: negative total depth between chunks")
         if (allocated(error)) exit checks
         call check(error, abs(m1 - m0) <= 1.0e-12_wp*abs(m0), &
                    "substep positivity: interior mass not conserved to round-off")
         if (allocated(error)) exit checks
         call check(error, saw_dry, &
                    "substep positivity: front never dried (check is vacuous)")
      end block checks
      if (allocated(fu)) deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_positivity_substep

   subroutine test_starts_dry(error)
      !! v2 intertidal test: cells that START DRY (wd_wet_dyn=0 at IC over
      !! the emerged bank) breathe from BELOW — they flood then re-dry.
      !!
      !! Geometry: 1-D channel (nx=40, ny=4).  Linear bed from H_DEEP=5 m
      !! at the western (deep) end to H_TOP=-1 m at the eastern (emerged)
      !! end.  Initial eta = 0 at the bank, ramping up to ETA_AMP=1.0 m at
      !! the deep end.  Total depth D0 = H_ref + eta: emerged cells (H_ref<0)
      !! have D0 = max(H_ref + eta, 0) ≈ 0 < dry_depth=0.05 ⇒ wd_wet_dyn=0.
      !!
      !! Sanity numbers (used to set chunk / run-time parameters):
      !!   emerged cells  = 7 of 40   (H_ref < 0)
      !!   starts-dry     = 7 of 40   (D0 < dry_depth = 0.05 m)
      !!   c_deep         ≈ 7.0 m/s   (sqrt(g * 5))
      !!   dt_inner       ≈ 5.0 s     (0.35 * dx / c_deep, dx=100 m)
      !!   travel time    ≈ 571 s     (lx / c_deep)
      !!   total run time ≈ 20000 s   (200 chunks × 20 steps × dt_inner)
      !!   est. seiche quarter-period ≈ 2285 s  (4*lx/c_deep)
      !! The run covers ~8.75 quarter-periods, so the wave reaches the bank
      !! multiple times and the emerged cells flood and re-dry.
      !!
      !! Gates:
      !!   1. Precondition: >=1 interior shelf column has wd_wet_dyn=0 at IC.
      !!   2. Flood: >=1 IC-dry column becomes wd_wet_dyn=1 during the run.
      !!   3. Re-dry: >=1 flooded column becomes wd_wet_dyn=0 again.
      !!   4. Positivity: min total depth over interior >= -1e-9 throughout.
      !!   5. Mass: interior mass conserved to round-off (1e-9 * |m0|).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: DRY_D = 0.05_wp, REWET_D = 0.10_wp
      real(wp), parameter :: H_DEEP = 5.0_wp, H_TOP = -1.0_wp
      real(wp), parameter :: ETA_AMP = 1.0_wp  ! IC wave amplitude at deep end
      integer, parameter :: NXP = 40, NYP = 4
      integer, parameter :: CHUNK = 20
      real(wp) :: lx, dx, dt_inner, m0, m1, min_d, frac, H_r, eta_ic
      integer :: i, j, nx, ny, ic, n_chunks
      integer :: n_starts_dry, n_shelf
      logical :: saw_flood, saw_redry
      logical, allocatable :: started_dry(:), ever_flooded(:)

      checks: block
         lx = 4000.0_wp
         dx = lx/real(NXP, wp)
         call make_grid(grid, NXP, NYP, dx, dx)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Linear bed: H_DEEP at west (frac=0), H_TOP at east (frac=1).
         ! IC eta: linear ramp from ETA_AMP at west to 0 at east, clamped
         ! so D = H_ref + eta >= 0.  Emerged bank cells (H_ref < 0) get
         ! eta clamped to -H_ref so D = 0 → wd_wet_dyn=0 at IC.
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               H_r = H_DEEP + (H_TOP - H_DEEP)*frac
               eta_ic = ETA_AMP*max(1.0_wp - 2.0_wp*frac, 0.0_wp)
               ! Clamp so D = H_ref + eta >= 0.
               eta_ic = max(eta_ic, -H_r)
               dyn%bt_work%bt_H_ref(i, j) = H_r
               dyn%bt_work%bt_eta(i, j) = eta_ic
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)

         ! Count interior shelf cells that START DRY (wd_wet_dyn == 0).
         n_shelf = 0
         n_starts_dry = 0
         do i = NGHOST + 1, NGHOST + NXP
            ! All j share same x-profile; check one representative row.
            if (dyn%bt_work%bt_H_ref(i, NGHOST + 2) < 0.0_wp) n_shelf = n_shelf + 1
            if (dyn%bt_work%wd_wet_dyn(i, NGHOST + 2) < 0.5_wp) n_starts_dry = n_starts_dry + 1
         end do
         ! Gate 1: precondition — at least one cell starts dry.
         call check(error, n_starts_dry >= 1, &
                    "starts_dry: no interior column has wd_wet_dyn==0 at IC (precondition failed)")
         if (allocated(error)) exit checks

         allocate (started_dry(nx), source=.false.)
         allocate (ever_flooded(nx), source=.false.)
         do i = 1, nx
            if (dyn%bt_work%wd_wet_dyn(i, NGHOST + 2) < 0.5_wp) started_dry(i) = .true.
         end do

         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)
         dt_inner = 0.35_wp*dx/sqrt(GRAVITY*H_DEEP)
         n_chunks = 200
         m0 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         min_d = huge(1.0_wp)
         saw_flood = .false.
         saw_redry = .false.

         do ic = 1, n_chunks
            call run_chunk(grid, metrics, dyn, cor, fu, fv, CHUNK, dt_inner)
            ! Gate 4: positivity throughout.
            min_d = min(min_d, minval(dyn%bt_work%bt_H_ref(NGHOST + 1:NGHOST + NXP, :) + &
                                      dyn%bt_work%bt_eta(NGHOST + 1:NGHOST + NXP, :)))
            ! Track flood and re-dry of originally-dry cells.
            do i = NGHOST + 1, NGHOST + NXP
               if (started_dry(i)) then
                  if (dyn%bt_work%wd_wet_dyn(i, NGHOST + 2) > 0.5_wp) then
                     saw_flood = .true.
                     ever_flooded(i) = .true.
                  else if (ever_flooded(i)) then
                     saw_redry = .true.
                  end if
               end if
            end do
         end do

         m1 = interior_mass(grid, dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta)
         ! Gate 4: positivity.
         call check(error, min_d >= -1.0e-9_wp, &
                    "starts_dry: negative total depth (positivity broken)")
         if (allocated(error)) exit checks
         ! Gate 5: mass conservation.
         call check(error, abs(m1 - m0) <= 1.0e-9_wp*abs(m0), &
                    "starts_dry: interior mass not conserved to round-off")
         if (allocated(error)) exit checks
         ! Gate 2: flood — some IC-dry cell became wet.
         call check(error, saw_flood, &
                    "starts_dry: no IC-dry cell ever flooded (v2 intertidal not working?)")
         if (allocated(error)) exit checks
         ! Gate 3: re-dry — some flooded cell dried again.
         call check(error, saw_redry, &
                    "starts_dry: no flooded cell ever re-dried (hysteresis stuck one way?)")
      end block checks
      if (allocated(fu)) deallocate (fu, fv)
      if (allocated(started_dry)) deallocate (started_dry)
      if (allocated(ever_flooded)) deallocate (ever_flooded)
      call cor%destroy()
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_starts_dry

   subroutine test_seed_criterion(error)
      !! v2 static-mask seed criterion (plan §10.1), exercised through the
      !! REAL production path (`ocean_state_seed_from_cfg` — the same call
      !! the driver makes).  Flat bathymetry b = max_depth everywhere, so
      !! one classification per case:
      !!
      !!   (i)  knob OFF, b = 1 m (< LAND_DEPTH_THRESHOLD = 2 m)
      !!        ⇒ wet_mask ≡ 0 — the v1 rest-depth land test is untouched.
      !!   (ii) knob ON, land_margin = 5 m, b = 1 m
      !!        ⇒ wet_mask ≡ 1 — intertidal terrain participates (THE v2
      !!        point: the same column that was static land in (i)).
      !!  (iii) knob ON, land_margin = 5 m, b = −6 m (bed 6 m above rest
      !!        MSL, above the flood headroom)
      !!        ⇒ wet_mask ≡ 0 — true land is still static.
      !!
      !! Host-only classification test: no metrics, no enter_data, no
      !! dynamics.  Case (iii) seeds negative h_layer/hTr — irrelevant
      !! here, only wet_mask is read before the state is destroyed.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DXS = 100.0_wp

      call seed_case(error, enable=.false., depth=1.0_wp, expect_wet=.false., &
                     what="(i) knob-off b=1 m must stay static land (v1 seed)")
      if (allocated(error)) return
      call seed_case(error, enable=.true., depth=1.0_wp, expect_wet=.true., &
                     what="(ii) knob-on b=1 m intertidal must stay dynamic (wet_mask=1)")
      if (allocated(error)) return
      call seed_case(error, enable=.true., depth=-6.0_wp, expect_wet=.false., &
                     what="(iii) knob-on bed above land_margin must stay static land")

   contains

      subroutine seed_case(err, enable, depth, expect_wet, what)
         type(error_type), allocatable, intent(out) :: err
         logical, intent(in) :: enable, expect_wet
         real(wp), intent(in) :: depth
         character(len=*), intent(in) :: what
         type(hgrid_t) :: grid
         type(ocean_state_t) :: state
         type(config_t) :: cfg
         logical :: all_wet, all_land

         cfg%sim_type = "ocean"
         cfg%nx = 4
         cfg%ny = 4
         cfg%dx = DXS
         cfg%dy = DXS
         cfg%nz_layers = 2
         cfg%nghost = NGHOST
         cfg%ocean%topo%max_depth = depth
         cfg%ocean%wetdry%enable = enable
         cfg%ocean%wetdry%land_margin = 5.0_wp

         call grid%init(4, 4, NGHOST, DXS, DXS)
         state%multilayer%nz_ml = 2
         call state%init(grid)
         call ocean_state_seed_from_cfg(state, grid, cfg)
         all_wet = all(state%multilayer%wet_mask == 1.0_wp)
         all_land = all(state%multilayer%wet_mask == 0.0_wp)
         if (expect_wet) then
            call check(err, all_wet, what)
         else
            call check(err, all_land, what)
         end if
         call state%destroy()
      end subroutine seed_case

   end subroutine test_seed_criterion

   subroutine test_emerged_seed_production(error)
      !! v2 emerged-bed SEED consistency through the REAL production path
      !! (`ocean_state_seed_from_cfg` — the same call the driver makes).
      !! This is the regime EVERY prior test avoids: `seed_criterion`
      !! ignores h_layer, `starts_dry` is barotropic-only (never seeds
      !! h_layer), and the restart test hand-mirrors the config with a
      !! POSITIVE bed.  A seamount with `edge_depth = peak < 0` pierces
      !! rest MSL, so the crest is EMERGED (b < 0, wet_mask=1 intertidal
      !! under the wetdry-aware `-land_margin` cutoff) while the flanks
      !! stay deep — a genuine mixed intertidal/deep column set built by
      !! the production topo dispatch.
      !!
      !! Gates:
      !!  1. No negative layer thickness anywhere (THE fix — unfixed code
      !!     seeds `h_layer = b/nz_ml < 0` on the emerged crest).
      !!  2. Every column `D = Σ h_layer >= 0` (`derive_bt_from_layers`
      !!     makes the barotropic depth `D = Σ h_layer`, so this is the
      !!     `H_ref + eta >= 0` guarantee the limiter needs at IC).
      !!  3. Emerged crest cells are LIVE (`wet_mask=1`) with each layer
      !!     floored to exactly `2·H_VANISHED > 0` (so the EOS `h>H_VANISHED`
      !!     branch fires and the layer survives the strict remap-drain gate).
      !!  4. The floor is SCOPED: every cell equals `max(b/nz_ml,
      !!     2·H_VANISHED)`, so deep flanks keep the un-floored `b/nz_ml` —
      !!     the fix only touches the near-zero band.
      !!  5. Tracers stay non-negative and recover the seeded S/T exactly
      !!     (`hTr/h` finite on both a deep and an emerged column).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DXS = 100.0_wp
      real(wp), parameter :: S0 = 35.0_wp, T0 = 10.0_wp
      real(wp), parameter :: MAXD = 5.0_wp, PEAK = -1.0_wp
      integer, parameter :: NXP = 16, NYP = 16, NZ_ML = 2
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(config_t) :: cfg
      integer :: i, j, k, nx, ny, idxS, idxT, n_emerged, n_deep
      real(wp) :: b_ij, expect_h, col_sum
      real(wp) :: minh, min_col, min_htr, max_scope_err
      integer :: i_deep, j_deep, i_emrg, j_emrg

      checks: block
         cfg%sim_type = "ocean"
         cfg%nx = NXP; cfg%ny = NYP; cfg%dx = DXS; cfg%dy = DXS
         cfg%nz_layers = NZ_ML; cfg%nghost = NGHOST
         cfg%initial_salinity = S0
         cfg%initial_temperature = T0
         cfg%ocean%topo%topo_config = "seamount"
         cfg%ocean%topo%max_depth = MAXD
         cfg%ocean%topo%edge_depth = PEAK   ! seamount crest, emerged (< 0)
         cfg%ocean%topo%slope_scale = 300.0_wp
         cfg%ocean%wetdry%enable = .true.
         cfg%ocean%wetdry%land_margin = 5.0_wp
         cfg%ocean%wetdry%dry_depth = 0.05_wp
         cfg%ocean%wetdry%rewet_depth = 0.10_wp

         call grid%init(NXP, NYP, NGHOST, DXS, DXS)
         state%multilayer%nz_ml = NZ_ML
         call state%init(grid)
         call ocean_state_seed_from_cfg(state, grid, cfg)

         nx = grid%nx_total; ny = grid%ny_total
         idxS = state%multilayer%idx_salinity
         idxT = state%multilayer%idx_temperature

         ! Sweep every column: positivity, D>=0, scope, and locate a
         ! representative deep + emerged column for the tracer check.
         minh = huge(1.0_wp)
         min_col = huge(1.0_wp)
         min_htr = huge(1.0_wp)
         max_scope_err = 0.0_wp
         n_emerged = 0
         n_deep = 0
         i_deep = 0; j_deep = 0; i_emrg = 0; j_emrg = 0
         do j = 1, ny
            do i = 1, nx
               b_ij = state%barotropic%b(i, j)
               col_sum = 0.0_wp
               do k = 1, NZ_ML
                  minh = min(minh, state%multilayer%h_layer(i, j, k))
                  min_htr = min(min_htr, &
                                state%multilayer%tracers(idxS)%hTr(i, j, k))
                  ! Gate 4: floor is exactly max(b/nz, 2·H_VANISHED).
                  expect_h = max(b_ij/real(NZ_ML, wp), 2.0_wp*H_VANISHED)
                  max_scope_err = max(max_scope_err, &
                                      abs(state%multilayer%h_layer(i, j, k) - expect_h))
                  col_sum = col_sum + state%multilayer%h_layer(i, j, k)
               end do
               min_col = min(min_col, col_sum)
               if (b_ij < 0.0_wp) then
                  n_emerged = n_emerged + 1
                  i_emrg = i; j_emrg = j
               else if (b_ij > 1.0_wp) then
                  n_deep = n_deep + 1
                  i_deep = i; j_deep = j
               end if
            end do
         end do

         ! Precondition: the seamount really did produce a MIX.
         call check(error, n_emerged > 0 .and. n_deep > 0, &
                    "emerged_seed: seamount did not produce a mixed emerged/deep domain")
         if (allocated(error)) exit checks
         ! Gate 1: no negative layer thickness.
         call check(error, minh >= 0.0_wp, &
                    "emerged_seed: negative h_layer seeded on the emerged crest (the fix)")
         if (allocated(error)) exit checks
         ! Gate 2: D = Σ h_layer >= 0 everywhere.
         call check(error, min_col >= 0.0_wp, &
                    "emerged_seed: negative barotropic depth D = Σ h_layer at IC")
         if (allocated(error)) exit checks
         ! Gate 3: emerged crest is LIVE and floored to 2·H_VANISHED.
         call check(error, state%multilayer%wet_mask(i_emrg, j_emrg) == 1.0_wp, &
                    "emerged_seed: emerged crest not classified intertidal (wet_mask/=1)")
         if (allocated(error)) exit checks
         call check(error, state%multilayer%h_layer(i_emrg, j_emrg, 1) == 2.0_wp*H_VANISHED, &
                    "emerged_seed: emerged crest layer not floored to 2*H_VANISHED")
         if (allocated(error)) exit checks
         ! Gate 4: floor scoped — every cell == max(b/nz, 2·H_VANISHED).
         call check(error, max_scope_err == 0.0_wp, &
                    "emerged_seed: floor not scoped to max(b/nz, 2*H_VANISHED)")
         if (allocated(error)) exit checks
         ! Gate 5: tracers non-negative + S/T recovered on deep AND emerged.
         call check(error, min_htr >= 0.0_wp, &
                    "emerged_seed: negative tracer mass hTr seeded")
         if (allocated(error)) exit checks
         call check(error, &
                    abs(state%multilayer%tracers(idxS)%hTr(i_emrg, j_emrg, 1)/ &
                        state%multilayer%h_layer(i_emrg, j_emrg, 1) - S0) < 1.0e-10_wp .and. &
                    abs(state%multilayer%tracers(idxT)%hTr(i_emrg, j_emrg, 1)/ &
                        state%multilayer%h_layer(i_emrg, j_emrg, 1) - T0) < 1.0e-10_wp, &
                    "emerged_seed: emerged-crest T/S not recovered from hTr/h (EOS would read garbage)")
         if (allocated(error)) exit checks
         call check(error, &
                    abs(state%multilayer%tracers(idxS)%hTr(i_deep, j_deep, 1)/ &
                        state%multilayer%h_layer(i_deep, j_deep, 1) - S0) < 1.0e-10_wp, &
                    "emerged_seed: deep-column S not recovered from hTr/h")
      end block checks
      call state%destroy()
   end subroutine test_emerged_seed_production

   subroutine test_emerged_beach_multilayer_step(error)
      !! v2 emerged-bed SEED-AND-STEP through the full multilayer split
      !! driver (`ocean_dyn_step_split`) — the coverage gap the blocker
      !! flagged: the existing wetdry tests never seed h_layer on an
      !! emerged bed and never run EOS/PGF/continuity over it.
      !!
      !! Geometry mirrors `test_starts_dry` (a linear beach, H_DEEP=5 m
      !! west → H_TOP=-1 m emerged east, tilted IC) but lifts it to the
      !! multilayer: `h_layer` is seeded through the REAL fixed
      !! `seed_h_layer_uniform_impl` (apply_wetdry_floor=.true.) from the
      !! clamped total depth D = max(H_ref + eta, 0), and uniform S/T so
      !! the flood is barotropic (density is flat ⇒ PGF≈0, but the EOS +
      !! PGF + continuity kernels still EXECUTE over the vanished
      !! intertidal columns and must not NaN).  bt_H_ref carries the true
      !! (partly negative) bed; `derive_bt_from_layers` re-derives bt_eta
      !! from Σ h_layer each stage, so D = Σ h_layer throughout.
      !!
      !! Gates:
      !!  1. IC: no negative h_layer (the fix — unfixed seeds b/nz < 0).
      !!  2. All h_layer finite (no NaN/Inf) every step — EOS/PGF/
      !!     continuity survive the vanished columns.
      !!  3. Positivity: min column depth Σ h_layer >= -1e-6 throughout.
      !!  4. Mass: interior total layer volume conserved (drift < 1e-7).
      !!  5. Flood + re-dry: an IC-dry emerged column wets, then dries.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DX = 100.0_wp
      real(wp), parameter :: H_DEEP = 5.0_wp, H_TOP = -1.0_wp
      real(wp), parameter :: ETA_AMP = 1.0_wp
      real(wp), parameter :: DRY_D = 0.05_wp, REWET_D = 0.10_wp
      integer, parameter :: NXP = 40, NYP = 4, NZ2 = 2
      ! N_INNER = 1 ⇒ the slow (outer) step runs at the BT-safe dt_inner,
      ! matching the stable test_starts_dry cadence: the sloshing layers +
      ! ALE remap tolerate the fine slow step where a coarse split (dt =
      ! N_INNER·dt_inner) over-CFLs the run-up front.
      integer, parameter :: N_INNER = 1, N_STEPS = 800
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(eos_t) :: eos
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc
      real(wp), allocatable :: b_field(:, :)
      logical, allocatable :: started_dry(:), ever_flooded(:)
      real(wp) :: dt_inner, dt, frac, H_r, eta_ic, csum
      real(wp) :: m0, m1, min_col, ic_minh
      integer :: i, j, k, nx, ny, step, jr
      logical :: all_finite, saw_flood, saw_redry

      checks: block
         call make_grid(grid, NXP, NYP, DX, DX)
         ms%nz_ml = NZ2
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ2)
         ! Wet/dry REQUIRES the per-layer PPM positivity guard (the BT
         ! limiter bounds the column; ppm_limit_pos bounds each layer
         ! under drying transports — validate_config enforces this).
         ct%use_ppm_limit_pos = .true.
         call cor%init(grid, nz_ml=NZ2)
         call eos%init(grid)
         call pgf%init(grid, nz_ml=NZ2)
         call hv%init(grid, nz_ml=NZ2)
         call bd%init(grid, nz_ml=NZ2)
         call ss%init(grid, nz_ml=NZ2)
         call va%init(grid, nz_ml=NZ2)
         call hd%init(grid, nz_ml=NZ2)
         call vd%init(grid, nz_ml=NZ2)
         call vmix%init(grid, nz_ml=NZ2)
         call dyn%init(grid, nz_ml=NZ2)
         ! Wet/dry is validated with the sigma ALE path: on sigma the
         ! layers scale to zero TOGETHER as D -> 0 (a dry column is
         ! "all layers vanished"), so the ALE remap keeps the two-layer
         ! partition and no single layer pinches out mid-column.  Without
         ! it the layers drift apart and a thin layer triggers a 1/h
         ! blowup.  Keep the interior closures off so the flood is a
         ! clean barotropic gravity-wave run-up over the emerged bank.
         call vc%init(grid, nz_ml=NZ2)
         vc%coord_type = VCOORD_SIGMA
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total
         jr = NGHOST + 2   ! representative interior row

         ! Linear emerged beach + tilted IC ⇒ total depth field D>=0.
         allocate (b_field(nx, ny))
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               H_r = H_DEEP + (H_TOP - H_DEEP)*frac
               eta_ic = ETA_AMP*max(1.0_wp - 2.0_wp*frac, 0.0_wp)
               ! bt_H_ref carries the true (partly negative) bed.
               dyn%bt_work%bt_H_ref(i, j) = H_r
               ! Seed the multilayer from the CLAMPED total depth D>=0.
               b_field(i, j) = max(H_r + eta_ic, 0.0_wp)
            end do
         end do

         ! Seed h_layer via the REAL fixed routine (floors emerged to
         ! H_VANISHED); tracers uniform so density is flat (barotropic
         ! flood) but EOS/PGF still run over the vanished columns.
         call seed_h_layer_uniform_impl(ms%h_layer, b_field, NZ2, &
                                        apply_wetdry_floor=.true.)
         do k = 1, NZ2
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     ms%h_layer(i, j, k)*eos%S_ref
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     ms%h_layer(i, j, k)*eos%T_ref
               end do
            end do
         end do

         ! Consistent bt_eta = Σ h_layer − bt_H_ref so the wd_wet_dyn seed
         ! (from bt_H_ref + bt_eta = D) classifies the emerged bank dry.
         do j = 1, ny
            do i = 1, nx
               csum = 0.0_wp
               do k = 1, NZ2
                  csum = csum + ms%h_layer(i, j, k)
               end do
               dyn%bt_work%bt_eta(i, j) = csum - dyn%bt_work%bt_H_ref(i, j)
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)

         ! Gate 1: IC no negative layer thickness.
         ic_minh = minval(ms%h_layer)
         call check(error, ic_minh >= 0.0_wp, &
                    "emerged_step: negative h_layer at IC (the fix)")
         if (allocated(error)) exit checks

         ! Record IC-dry emerged columns + IC interior mass.
         allocate (started_dry(nx), source=.false.)
         allocate (ever_flooded(nx), source=.false.)
         do i = 1, nx
            if (dyn%bt_work%wd_wet_dyn(i, jr) < 0.5_wp) started_dry(i) = .true.
         end do
         m0 = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  m0 = m0 + ms%h_layer(i, j, k)
               end do
            end do
         end do

         dt_inner = 0.35_wp*DX/sqrt(GRAVITY*H_DEEP)
         dt = real(N_INNER, wp)*dt_inner
         min_col = huge(1.0_wp)
         all_finite = .true.
         saw_flood = .false.
         saw_redry = .false.

         !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ms%enter_data(); call ct%enter_data(); call cor%enter_data()
         call pgf%enter_data(); call hv%enter_data(); call bd%enter_data()
         call ss%enter_data(); call va%enter_data(); call hd%enter_data()
         call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
         call vc%enter_data()
         call make_cartesian_metrics(metrics, grid)

         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, &
                                      hv, bd, ss, va, hd, vd, vmix, ms, dt, N_INNER, &
                                      vcoord=vc)
            !$acc update self(ms%h_layer, dyn%bt_work%wd_wet_dyn)
            ! Gate 2/3: finite + positivity over interior columns.
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP
                  csum = 0.0_wp
                  do k = 1, NZ2
                     if (.not. ieee_is_finite(ms%h_layer(i, j, k))) all_finite = .false.
                     csum = csum + ms%h_layer(i, j, k)
                  end do
                  min_col = min(min_col, csum)
               end do
            end do
            ! Gate 5: track flood + re-dry of IC-dry emerged columns.
            do i = NGHOST + 1, NGHOST + NXP
               if (started_dry(i)) then
                  if (dyn%bt_work%wd_wet_dyn(i, jr) > 0.5_wp) then
                     saw_flood = .true.
                     ever_flooded(i) = .true.
                  else if (ever_flooded(i)) then
                     saw_redry = .true.
                  end if
               end if
            end do
         end do

         !$acc update self(ms%h_layer)
         call destroy_cartesian_metrics(metrics)
         call vc%exit_data()
         call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
         call hd%exit_data(); call va%exit_data(); call ss%exit_data()
         call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
         call cor%exit_data(); call ct%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         m1 = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  m1 = m1 + ms%h_layer(i, j, k)
               end do
            end do
         end do

         ! Gate 2: no NaN/Inf ever.
         call check(error, all_finite, &
                    "emerged_step: h_layer went non-finite (EOS/PGF NaN on vanished column?)")
         if (allocated(error)) exit checks
         ! Gate 3: positivity throughout.
         call check(error, min_col >= -1.0e-6_wp, &
                    "emerged_step: column depth Σ h_layer went negative")
         if (allocated(error)) exit checks
         ! Gate 4: interior mass conserved.
         call check(error, abs(m1 - m0) <= 1.0e-7_wp*abs(m0), &
                    "emerged_step: interior layer volume not conserved")
         if (allocated(error)) exit checks
         ! Gate 5: an emerged column flooded and later re-dried.
         call check(error, saw_flood, &
                    "emerged_step: no IC-dry emerged column ever flooded")
         if (allocated(error)) exit checks
         call check(error, saw_redry, &
                    "emerged_step: no flooded emerged column ever re-dried")
      end block checks
      if (allocated(b_field)) deallocate (b_field)
      if (allocated(started_dry)) deallocate (started_dry, ever_flooded)
      call vc%destroy()
      call dyn%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy()
      call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy()
      call pgf%destroy(); call eos%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine test_emerged_beach_multilayer_step

   subroutine test_emerged_tracer_gradient(error)
      !! Fix-2 (blocker) DISCRIMINATING test: the same emerged-beach flood as
      !! `test_emerged_beach_multilayer_step`, but with a NON-UNIFORM salinity
      !! gradient so the ALE remap genuinely TRANSPORTS a tracer through the
      !! drain band.  The prior test seeds UNIFORM S/T and gates only VOLUME
      !! (Σ h_layer), so it is blind to the seed-floor × remap-drain-gate
      !! collision (a layer seeded at exactly H_VANISHED fails the strict
      !! `h_old > H_VANISHED` drain gate ⇒ its concentration is zeroed ⇒ the
      !! seeded hS is DESTROYED on the first regrid).  With the 2·H_VANISHED
      !! seed floor the seeded layer clears the gate and Σ hS / Σ hT conserve.
      !!
      !! On the UN-FIXED (exact-H_VANISHED) seed this test FAILS: the emerged
      !! columns lose their seeded salt content on step 1 → Σ hS drops by the
      !! seeded emerged mass → the conservation gate trips.  (Verified by
      !! stashing Fix 1 and re-running — the discriminating check.)
      !!
      !! Gates: no NaN/Inf; Σ hS and Σ hT (all layers+interior cells)
      !! conserved to ≤ 1e-12 relative through flood → drain → reflood;
      !! positivity (hS, hT >= 0) throughout.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DX = 100.0_wp
      real(wp), parameter :: H_DEEP = 5.0_wp, H_TOP = -1.0_wp
      real(wp), parameter :: ETA_AMP = 1.0_wp
      real(wp), parameter :: DRY_D = 0.05_wp, REWET_D = 0.10_wp
      real(wp), parameter :: S_LO = 30.0_wp, S_HI = 38.0_wp, T0 = 10.0_wp
      integer, parameter :: NXP = 40, NYP = 4, NZ2 = 2
      integer, parameter :: N_INNER = 1, N_STEPS = 800
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(eos_t) :: eos
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc
      real(wp), allocatable :: b_field(:, :)
      real(wp) :: dt_inner, dt, frac, H_r, eta_ic, csum, s_col_val
      real(wp) :: hs0, hs1, ht0, ht1, rel_s, rel_t, min_tr
      integer :: i, j, k, nx, ny, step, idxS, idxT

      checks: block
         call make_grid(grid, NXP, NYP, DX, DX)
         ms%nz_ml = NZ2
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ2)
         ct%use_ppm_limit_pos = .true.
         call cor%init(grid, nz_ml=NZ2)
         call eos%init(grid)
         call pgf%init(grid, nz_ml=NZ2)
         call hv%init(grid, nz_ml=NZ2)
         call bd%init(grid, nz_ml=NZ2)
         call ss%init(grid, nz_ml=NZ2)
         call va%init(grid, nz_ml=NZ2)
         call hd%init(grid, nz_ml=NZ2)
         call vd%init(grid, nz_ml=NZ2)
         call vmix%init(grid, nz_ml=NZ2)
         call dyn%init(grid, nz_ml=NZ2)
         call vc%init(grid, nz_ml=NZ2)
         vc%coord_type = VCOORD_SIGMA
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total
         idxS = ms%idx_salinity
         idxT = ms%idx_temperature

         allocate (b_field(nx, ny))
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               H_r = H_DEEP + (H_TOP - H_DEEP)*frac
               eta_ic = ETA_AMP*max(1.0_wp - 2.0_wp*frac, 0.0_wp)
               dyn%bt_work%bt_H_ref(i, j) = H_r
               b_field(i, j) = max(H_r + eta_ic, 0.0_wp)
            end do
         end do

         call seed_h_layer_uniform_impl(ms%h_layer, b_field, NZ2, &
                                        apply_wetdry_floor=.true.)
         ! Non-uniform S: a west→east gradient so the remap genuinely
         ! transports the tracer.  T uniform (keeps the flood barotropic-ish
         ! while still routing hT through the remap).
         do k = 1, NZ2
            do j = 1, ny
               do i = 1, nx
                  frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
                  s_col_val = S_LO + (S_HI - S_LO)*max(0.0_wp, min(1.0_wp, frac))
                  ms%tracers(idxS)%hTr(i, j, k) = ms%h_layer(i, j, k)*s_col_val
                  ms%tracers(idxT)%hTr(i, j, k) = ms%h_layer(i, j, k)*T0
               end do
            end do
         end do

         do j = 1, ny
            do i = 1, nx
               csum = 0.0_wp
               do k = 1, NZ2
                  csum = csum + ms%h_layer(i, j, k)
               end do
               dyn%bt_work%bt_eta(i, j) = csum - dyn%bt_work%bt_H_ref(i, j)
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)

         ! IC interior tracer content.
         hs0 = 0.0_wp; ht0 = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  hs0 = hs0 + ms%tracers(idxS)%hTr(i, j, k)
                  ht0 = ht0 + ms%tracers(idxT)%hTr(i, j, k)
               end do
            end do
         end do

         dt_inner = 0.35_wp*DX/sqrt(GRAVITY*H_DEEP)
         dt = real(N_INNER, wp)*dt_inner
         min_tr = huge(1.0_wp)

         !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ms%enter_data(); call ct%enter_data(); call cor%enter_data()
         call pgf%enter_data(); call hv%enter_data(); call bd%enter_data()
         call ss%enter_data(); call va%enter_data(); call hd%enter_data()
         call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
         call vc%enter_data()
         call make_cartesian_metrics(metrics, grid)

         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, &
                                      hv, bd, ss, va, hd, vd, vmix, ms, dt, N_INNER, &
                                      vcoord=vc)
         end do

         !$acc update self(ms%h_layer)
         !$acc update self(ms%tracers(idxS)%hTr, ms%tracers(idxT)%hTr)
         call destroy_cartesian_metrics(metrics)
         call vc%exit_data()
         call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
         call hd%exit_data(); call va%exit_data(); call ss%exit_data()
         call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
         call cor%exit_data(); call ct%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         hs1 = 0.0_wp; ht1 = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  hs1 = hs1 + ms%tracers(idxS)%hTr(i, j, k)
                  ht1 = ht1 + ms%tracers(idxT)%hTr(i, j, k)
                  min_tr = min(min_tr, ms%tracers(idxS)%hTr(i, j, k))
                  min_tr = min(min_tr, ms%tracers(idxT)%hTr(i, j, k))
                  if (.not. ieee_is_finite(ms%tracers(idxS)%hTr(i, j, k))) &
                     min_tr = -huge(1.0_wp)
                  if (.not. ieee_is_finite(ms%tracers(idxT)%hTr(i, j, k))) &
                     min_tr = -huge(1.0_wp)
               end do
            end do
         end do

         rel_s = abs(hs1 - hs0)/abs(hs0)
         rel_t = abs(ht1 - ht0)/abs(ht0)

         ! Positivity + finiteness.
         call check(error, min_tr >= 0.0_wp, &
                    "emerged_tracer_grad: tracer went negative/non-finite")
         if (allocated(error)) exit checks
         ! THE discriminating gate: Σ hS conserved to round-off (fails on the
         ! un-fixed exact-H_VANISHED seed, which destroys emerged salt at step 1).
         call check(error, rel_s <= 1.0e-12_wp, &
                    "emerged_tracer_grad: Σ hS not conserved (seed-floor/drain-gate collision?)")
         if (allocated(error)) exit checks
         call check(error, rel_t <= 1.0e-12_wp, &
                    "emerged_tracer_grad: Σ hT not conserved")
      end block checks
      if (allocated(b_field)) deallocate (b_field)
      call vc%destroy()
      call dyn%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy()
      call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy()
      call pgf%destroy(); call eos%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine test_emerged_tracer_gradient

   subroutine test_vanished_column_flood(error)
      !! Fix-3 (Finding A / zero-depth plan Phase 1) controlled test: one
      !! interior column seeded FULLY VANISHED (all layers at the 2·H_VANISHED
      !! floor ⇒ D ≈ nz·2·H_VANISHED, well below dry_depth) sits adjacent to
      !! deep wet columns carrying a NON-UNIFORM salinity gradient, under
      !! gentle sub-CFL forcing (an SSH tilt that overtops the vanished column
      !! within the run).  Steps the FULL split driver.  This is the load-
      !! bearing D≈0 experiment the plan calls for — it separates a real
      !! zero-depth kernel bug from the forcing-CFL confound (the tilt is
      !! sub-CFL by construction).
      !!
      !! Gates: no NaN/Inf ever (EOS/PGF survive the vanished column with a
      !! density gradient next door); Σ hS / Σ hT conserved to round-off;
      !! positivity throughout; the vanished column FLOODS (wd_wet_dyn → 1).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DX = 100.0_wp
      real(wp), parameter :: H_DEEP = 6.0_wp
      real(wp), parameter :: DRY_D = 0.05_wp, REWET_D = 0.10_wp
      real(wp), parameter :: S_LO = 30.0_wp, S_HI = 38.0_wp, T0 = 10.0_wp
      real(wp), parameter :: ETA_TILT = 1.5_wp
      integer, parameter :: NXP = 20, NYP = 4, NZ2 = 2
      integer, parameter :: N_INNER = 1, N_STEPS = 600
      integer, parameter :: I_VANISH = NGHOST + 10   ! interior vanished column
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(eos_t) :: eos
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(ocean_dyn_t) :: dyn
      type(ocean_vcoord_t) :: vc
      real(wp), allocatable :: b_field(:, :)
      real(wp) :: dt_inner, dt, csum, s_col_val, frac
      real(wp) :: hs0, hs1, ht0, ht1, rel_s, rel_t, min_tr
      integer :: i, j, k, nx, ny, step, idxS, idxT
      logical :: all_finite, vanish_flooded

      checks: block
         call make_grid(grid, NXP, NYP, DX, DX)
         ms%nz_ml = NZ2
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ2)
         ct%use_ppm_limit_pos = .true.
         call cor%init(grid, nz_ml=NZ2)
         call eos%init(grid)
         call pgf%init(grid, nz_ml=NZ2)
         call hv%init(grid, nz_ml=NZ2)
         call bd%init(grid, nz_ml=NZ2)
         call ss%init(grid, nz_ml=NZ2)
         call va%init(grid, nz_ml=NZ2)
         call hd%init(grid, nz_ml=NZ2)
         call vd%init(grid, nz_ml=NZ2)
         call vmix%init(grid, nz_ml=NZ2)
         call dyn%init(grid, nz_ml=NZ2)
         call vc%init(grid, nz_ml=NZ2)
         vc%coord_type = VCOORD_SIGMA
         vmix%use_closure = .false.
         vmix%use_kpp = .false.
         vd%K_v_tracer = 0.0_wp
         vd%K_v_momentum = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total
         idxS = ms%idx_salinity
         idxT = ms%idx_temperature

         ! Flat deep bed everywhere EXCEPT the one interior vanished column.
         ! A gentle SSH tilt (sub-CFL) drives water toward the vanished
         ! column so it overtops and floods.
         allocate (b_field(nx, ny))
         do j = 1, ny
            do i = 1, nx
               frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
               dyn%bt_work%bt_H_ref(i, j) = H_DEEP
               ! IC total depth: deep + a linear SSH tilt (max ETA_TILT at west).
               b_field(i, j) = max(H_DEEP + ETA_TILT*(1.0_wp - frac), 0.0_wp)
            end do
         end do
         ! Force the target column vanished: seed its D at ~0 (bed at MSL, no
         ! water).  bt_H_ref stays H_DEEP so overtopping is reachable.
         do j = 1, ny
            b_field(I_VANISH, j) = 0.0_wp
         end do

         call seed_h_layer_uniform_impl(ms%h_layer, b_field, NZ2, &
                                        apply_wetdry_floor=.true.)
         ! Non-uniform S gradient across the domain (the vanished column too).
         do k = 1, NZ2
            do j = 1, ny
               do i = 1, nx
                  frac = (real(i - NGHOST, wp) - 0.5_wp)/real(NXP, wp)
                  s_col_val = S_LO + (S_HI - S_LO)*max(0.0_wp, min(1.0_wp, frac))
                  ms%tracers(idxS)%hTr(i, j, k) = ms%h_layer(i, j, k)*s_col_val
                  ms%tracers(idxT)%hTr(i, j, k) = ms%h_layer(i, j, k)*T0
               end do
            end do
         end do

         do j = 1, ny
            do i = 1, nx
               csum = 0.0_wp
               do k = 1, NZ2
                  csum = csum + ms%h_layer(i, j, k)
               end do
               dyn%bt_work%bt_eta(i, j) = csum - dyn%bt_work%bt_H_ref(i, j)
            end do
         end do
         call enable_wetdry(dyn%bt_work, grid, DRY_D, REWET_D)

         hs0 = 0.0_wp; ht0 = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  hs0 = hs0 + ms%tracers(idxS)%hTr(i, j, k)
                  ht0 = ht0 + ms%tracers(idxT)%hTr(i, j, k)
               end do
            end do
         end do

         dt_inner = 0.30_wp*DX/sqrt(GRAVITY*H_DEEP)
         dt = real(N_INNER, wp)*dt_inner
         all_finite = .true.
         vanish_flooded = .false.

         !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ms%enter_data(); call ct%enter_data(); call cor%enter_data()
         call pgf%enter_data(); call hv%enter_data(); call bd%enter_data()
         call ss%enter_data(); call va%enter_data(); call hd%enter_data()
         call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
         call vc%enter_data()
         call make_cartesian_metrics(metrics, grid)

         do step = 1, N_STEPS
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, &
                                      hv, bd, ss, va, hd, vd, vmix, ms, dt, N_INNER, &
                                      vcoord=vc)
            !$acc update self(ms%h_layer, dyn%bt_work%wd_wet_dyn)
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP
                  do k = 1, NZ2
                     if (.not. ieee_is_finite(ms%h_layer(i, j, k))) all_finite = .false.
                  end do
               end do
            end do
            if (dyn%bt_work%wd_wet_dyn(I_VANISH, NGHOST + 2) > 0.5_wp) &
               vanish_flooded = .true.
         end do

         !$acc update self(ms%h_layer)
         !$acc update self(ms%tracers(idxS)%hTr, ms%tracers(idxT)%hTr)
         call destroy_cartesian_metrics(metrics)
         call vc%exit_data()
         call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
         call hd%exit_data(); call va%exit_data(); call ss%exit_data()
         call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
         call cor%exit_data(); call ct%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         hs1 = 0.0_wp; ht1 = 0.0_wp; min_tr = huge(1.0_wp)
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               do k = 1, NZ2
                  hs1 = hs1 + ms%tracers(idxS)%hTr(i, j, k)
                  ht1 = ht1 + ms%tracers(idxT)%hTr(i, j, k)
                  min_tr = min(min_tr, ms%tracers(idxS)%hTr(i, j, k))
                  min_tr = min(min_tr, ms%tracers(idxT)%hTr(i, j, k))
               end do
            end do
         end do
         rel_s = abs(hs1 - hs0)/abs(hs0)
         rel_t = abs(ht1 - ht0)/abs(ht0)

         call check(error, all_finite, &
                    "vanished_flood: h_layer went non-finite (zero-depth slow-path NaN?)")
         if (allocated(error)) exit checks
         call check(error, min_tr >= 0.0_wp, &
                    "vanished_flood: tracer went negative")
         if (allocated(error)) exit checks
         call check(error, rel_s <= 1.0e-12_wp, &
                    "vanished_flood: Σ hS not conserved across the vanished column")
         if (allocated(error)) exit checks
         call check(error, rel_t <= 1.0e-12_wp, &
                    "vanished_flood: Σ hT not conserved")
         if (allocated(error)) exit checks
         call check(error, vanish_flooded, &
                    "vanished_flood: the seeded-vanished interior column never flooded")
      end block checks
      if (allocated(b_field)) deallocate (b_field)
      call vc%destroy()
      call dyn%destroy(); call vmix%destroy(); call vd%destroy(); call hd%destroy()
      call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy()
      call pgf%destroy(); call eos%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine test_vanished_column_flood

end module test_ocean_wetdry
