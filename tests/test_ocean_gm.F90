!! Analytical + device tests for Gent-McWilliams thickness diffusion
!! (capability [2], `rdb_ocean_gm`).  All cases RUN THE DEVICE KERNELS via
!! `gm_compute_transports` (the production entry) with a linear EOS and the
!! slopes slot supplying the stored isopycnal slope.
!!
!! Bottom-up convention: k=1 bed, k=nz surface; interface K=1 bed,
!! K=nz+1 surface (both zero slope).  GM transports `uhD`/`vhD` satisfy
!! `Sum_k uhD = 0` per face by the surface/bed BCs (the column recurrence
!! closes at the surface via `uhD(nz) = -uhtot`).
!!
!!  1. gm_conserves_mass     — Sum_{i,j,k} h unchanged after a GM step.
!!  2. gm_flattens_isopycnal — a tilted 2-layer interface relaxes toward
!!                             flat under GM-only (tilt monotone-decreasing).
!!  3. gm_tracer_conserves   — tracer content conserved through the GM
!!                             operator (`continuity_gm_apply`).
!!  4. gm_slope_limiter      — S >> slope_max ⇒ bounded Psi, h >= H_VANISHED.
!!  5. gm_gm_src_sign        — gm_src >= 0 (PE release) for a stable column.
!!  6. gm_sequential_partial_cell — the 1-degree Southern Ocean failure: an
!!                             8.6 cm partial bed cell next to fillers under a
!!                             dome, drained at the cap on all four faces.
!!                             GM computed from the thickness the dynamics
!!                             LEFT keeps it >= H_VANISHED and conserves
!!                             volume + T/S content; the same transport
!!                             computed from the stage-ENTRY thickness (the
!!                             old fold) takes it negative.
module test_ocean_gm
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports
   use rdb_continuity, only: continuity_t, continuity_gm_apply, TR_MODE_ADVECT
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_gm_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3/degC
   real(wp), parameter :: BETA_S = 0.78_wp     ! kg/m^3/PSU
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: KHTH = 1000.0_wp

contains

   subroutine collect_ocean_gm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("gm_conserves_mass", test_conserves_mass), &
                  new_unittest("gm_flattens_isopycnal", test_flattens), &
                  new_unittest("gm_tracer_conserves", test_tracer_conserves), &
                  new_unittest("gm_wall_no_leak", test_wall_no_leak), &
                  new_unittest("gm_slope_limiter", test_slope_limiter), &
                  new_unittest("gm_gm_src_sign", test_gm_src_sign), &
                  new_unittest("gm_sequential_partial_cell", test_sequential_partial_cell) &
                  ]
   end subroutine collect_ocean_gm_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine make_eos

   subroutine setup_slopes(sl, grid, nz)
      type(ocean_slopes_t), intent(inout) :: sl
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call sl%init(grid, nz_ml=nz)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%kd_smooth = 0.0_wp
      sl%min_dz_for_n2 = 1.0_wp
   end subroutine setup_slopes

   subroutine setup_gm(gm, grid, nz, slope_max)
      type(ocean_gm_t), intent(inout) :: gm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in), optional :: slope_max
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%khth = KHTH
      gm%khth_max_cfl = 1.0e6_wp   ! effectively no CFL clamp for the analytic tests
      gm%khth_slope_max = 0.01_wp
      if (present(slope_max)) gm%khth_slope_max = slope_max
      gm%rho0 = RHO0
   end subroutine setup_gm

   subroutine map_in(ms, sl, gm)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      !$acc enter data copyin(ms)
      call ms%enter_data()
      ! Bed datum of the slopes' geopotential interface heights: every
      ! case here has a flat bed under a constant-depth column at rest
      ! (eta = 0), so D = Sum_k h.  Host-side, BEFORE the slot's map.
      call sl%set_bathymetry(sum(ms%h_layer, dim=3))
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
   end subroutine map_in

   subroutine map_out(ms, sl, gm)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      !$acc update self(gm%uhD, gm%vhD, gm%gm_src, gm%khth_u, gm%khth_v)
      !$acc update self(sl%slope_x, sl%slope_y, sl%n2_u, sl%n2_v)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine fill_tilted_TS(ms, grid, nz, dz, gx)
      !! T tilts in x (warmer east), S uniform, uniform thickness: a stable
      !! column (warm/light over cold/dense via GZ) with a non-zero
      !! horizontal density gradient ⇒ a tilted isopycnal.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, gx
      real(wp), parameter :: GZ = 0.02_wp
      integer :: i, j, k, ni, nj
      real(wp) :: z_c, x_c, t_val
      ni = grid%nx_total
      nj = grid%ny_total
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         z_c = (real(k, wp) - 0.5_wp)*dz
         do j = 1, nj
            do i = 1, ni
               x_c = real(i, wp)*grid%dx
               t_val = T0 + GZ*z_c + gx*x_c
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*dz
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*dz
            end do
         end do
      end do
   end subroutine fill_tilted_TS

   ! ------------------------------------------------------------------
   ! Test 1: column-sum conservation Sum_k uhD = 0 / Sum_k vhD = 0
   ! ------------------------------------------------------------------
   subroutine test_conserves_mass(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 2000.0_wp, DZ = 50.0_wp
      real(wp) :: csum, mx
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         call fill_tilted_TS(ms, grid, NZ, DZ, 1.0e-4_wp)

         call map_in(ms, sl, gm)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call map_out(ms, sl, gm)

         ! Sum_k uhD must be 0 on every u-face; same for vhD.
         mx = 0.0_wp
         do j = 1, nj
            do i = 1, ni + 1
               csum = 0.0_wp
               do k = 1, NZ
                  csum = csum + gm%uhD(i, j, k)
               end do
               mx = max(mx, abs(csum))
            end do
         end do
         do j = 1, nj + 1
            do i = 1, ni
               csum = 0.0_wp
               do k = 1, NZ
                  csum = csum + gm%vhD(i, j, k)
               end do
               mx = max(mx, abs(csum))
            end do
         end do
         call check(error, mx < 1.0e-6_wp, &
                    "GM column transport must close (Sum_k uhD = 0)")
         if (allocated(error)) exit checks
         ! Non-trivial: at least one face must carry a finite transport.
         call check(error, maxval(abs(gm%uhD)) > 0.0_wp, &
                    "GM must produce a non-zero bolus transport for a tilted column")
      end block checks
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_conserves_mass

   ! ------------------------------------------------------------------
   ! Test 4: slope limiter — S >> slope_max ⇒ bounded Psi, h stays >= floor
   ! ------------------------------------------------------------------
   subroutine test_slope_limiter(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 6
      real(wp), parameter :: DX = 2000.0_wp, DZ = 50.0_wp
      real(wp) :: mx_uhD, area_budget
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ, slope_max=0.01_wp)
         call make_eos(eos)
         ! Huge horizontal T gradient ⇒ steep (near-vertical) isopycnal,
         ! S = drdx/sqrt(...) ~ 1 >> slope_max.
         call fill_tilted_TS(ms, grid, NZ, DZ, 5.0e-2_wp)

         call map_in(ms, sl, gm)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call map_out(ms, sl, gm)

         ! The mass-availability limiter caps each layer transport at the
         ! donor budget areaT*(h-H_VANISHED)/(4 dt).  No uhD may exceed it.
         area_budget = (DX*DX)*(DZ - H_VANISHED)/(4.0_wp*DT)
         mx_uhD = maxval(abs(gm%uhD))
         call check(error, mx_uhD <= area_budget*(1.0_wp + 1.0e-9_wp), &
                    "GM transport must stay within the donor mass budget")
         if (allocated(error)) exit checks
         ! And it must be finite (no NaN/Inf runaway).
         call check(error, mx_uhD == mx_uhD .and. mx_uhD < huge(1.0_wp), &
                    "GM transport must be finite under a near-vertical slope")
      end block checks
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_slope_limiter

   ! ------------------------------------------------------------------
   ! Test 5: gm_src >= 0 (PE release) for a stably stratified tilted column
   ! ------------------------------------------------------------------
   subroutine test_gm_src_sign(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 2000.0_wp, DZ = 50.0_wp
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         call make_eos(eos)
         call fill_tilted_TS(ms, grid, NZ, DZ, 1.0e-4_wp)

         call map_in(ms, sl, gm)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call map_out(ms, sl, gm)

         ! gm_src >= 0 alone is vacuous (it is a sum of squares).  The
         ! discriminating checks: positive where tilted (needs KH·S²·N²·h>0),
         ! finite everywhere (no NaN/Inf), and positive across the BULK of the
         ! tilted interior (not a single stray cell — guards a mis-indexed sum).
         call check(error, maxval(gm%gm_src) > 0.0_wp, &
                    "gm_src must be positive where the isopycnal is tilted")
         if (allocated(error)) exit checks
         call check(error, (.not. any(gm%gm_src /= gm%gm_src)) .and. &
                    maxval(gm%gm_src) < 1.0e30_wp, "gm_src must be finite (no NaN/Inf)")
         if (allocated(error)) exit checks
         call check(error, count(gm%gm_src > 0.0_wp) >= (grid%nx_phys*grid%ny_phys)/2, &
                    "gm_src must be positive across the tilted interior, not a stray cell")
      end block checks
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_gm_src_sign

   ! ------------------------------------------------------------------
   ! Test 2: a tilted 2-layer interface flattens under GM-only
   ! ------------------------------------------------------------------
   subroutine test_flattens(error)
      !! Two layers (cold/dense at the bed k=1, warm/light at the surface
      !! k=nz=2), uniform T/S horizontally but a TILTED interface (the
      !! bottom layer is thicker to the west).  The slope-formula's
      !! interface-tilt term gives a non-zero neutral slope; GM should move
      !! mass to flatten the interface — the west-vs-east bottom-layer
      !! thickness difference must shrink monotonically.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      type(eos_t) :: eos
      integer, parameter :: NX = 8, NY = 3, NZ = 2
      real(wp), parameter :: DX = 2000.0_wp, HTOT = 100.0_wp
      real(wp), parameter :: DTL = 600.0_wp
      integer, parameter :: NSTEP = 40
      real(wp) :: tilt0, tilt, h1w, h1e, prev
      integer :: i, j, k, ni, nj, it, iw, ie
      logical :: monotone
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         iw = NGHOST + 1
         ie = NGHOST + grid%nx_phys

         ! Tilted interface: bottom layer thicker to the west, surface layer
         ! takes up the rest so the column total is HTOT everywhere.
         ! Stable 2-layer: cold (dense) bottom, warm (light) top, uniform
         ! horizontally so the slope is driven purely by the interface tilt.
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do j = 1, nj
            do i = 1, ni
               h1w = 0.5_wp*HTOT + 20.0_wp*(real(grid%nx_phys + 1, wp)*0.5_wp - real(i - NGHOST, wp)) &
                     /real(grid%nx_phys, wp)
               if (h1w < 10.0_wp) h1w = 10.0_wp
               if (h1w > HTOT - 10.0_wp) h1w = HTOT - 10.0_wp
               ms%h_layer(i, j, 1) = h1w
               ms%h_layer(i, j, 2) = HTOT - h1w
               ! cold bottom (k=1), warm surface (k=2); uniform in x.
               ms%tracers(ms%idx_temperature)%hTr(i, j, 1) = 8.0_wp*ms%h_layer(i, j, 1)
               ms%tracers(ms%idx_temperature)%hTr(i, j, 2) = 16.0_wp*ms%h_layer(i, j, 2)
               ms%tracers(ms%idx_salinity)%hTr(i, j, 1) = S0*ms%h_layer(i, j, 1)
               ms%tracers(ms%idx_salinity)%hTr(i, j, 2) = S0*ms%h_layer(i, j, 2)
            end do
         end do

         call map_in_ct(ms, sl, gm, ct)

         j = NGHOST + 1
         !$acc update self(ms%h_layer)
         tilt0 = abs(ms%h_layer(iw, j, 1) - ms%h_layer(ie, j, 1))
         tilt = tilt0
         prev = tilt0 + 1.0_wp
         monotone = .true.
         do it = 1, NSTEP
            call ocean_slopes_compute(grid, metrics, eos, sl, ms, DTL)
            call gm_compute_transports(grid, metrics, gm, sl, ms, DTL)
            call continuity_gm_apply(grid, metrics, ct, ms, gm, DTL, 1.0_wp, TR_MODE_ADVECT)
            !$acc update self(ms%h_layer)
            tilt = abs(ms%h_layer(iw, j, 1) - ms%h_layer(ie, j, 1))
            if (tilt > prev + 1.0e-9_wp) monotone = .false.
            prev = tilt
         end do

         call map_out_ct(ms, sl, gm, ct)

         call check(error, tilt < tilt0, &
                    "GM must reduce the interface tilt (flatten isopycnals)")
         if (allocated(error)) exit checks
         call check(error, monotone, "tilt must decrease monotonically under GM")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_flattens

   ! ------------------------------------------------------------------
   ! Test 3: tracer content conserved through the GM operator
   ! ------------------------------------------------------------------
   subroutine test_tracer_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 6
      real(wp), parameter :: DX = 2000.0_wp, DZ = 50.0_wp
      real(wp) :: hsum0, hsum, tsum0, tsum, ssum0, ssum, relh, relt, rels
      integer :: i, j, k, i0, i1, j0, j1
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         call make_eos(eos)
         call fill_tilted_TS(ms, grid, NZ, DZ, 1.0e-4_wp)

         ! GM is conservative over the CLOSED discrete domain: the operator
         ! zeroes the bolus transport on every non-periodic physical edge
         ! face (here all four are walls), so the full-array sum and the
         ! physical-domain sum are both invariant; the full array is the
         ! stricter bookkeeping (nothing may appear in the halo either).
         i0 = 1
         i1 = grid%nx_total
         j0 = 1
         j1 = grid%ny_total

         hsum0 = sum_phys_h(ms, i0, i1, j0, j1, NZ)
         tsum0 = sum_phys_tr(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)
         ssum0 = sum_phys_tr(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)

         call map_in_ct(ms, sl, gm, ct)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call continuity_gm_apply(grid, metrics, ct, ms, gm, DT, 1.0_wp, TR_MODE_ADVECT)
         call map_out_ct(ms, sl, gm, ct)

         hsum = sum_phys_h(ms, i0, i1, j0, j1, NZ)
         tsum = sum_phys_tr(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)
         ssum = sum_phys_tr(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)

         relh = abs(hsum - hsum0)/abs(hsum0)
         relt = abs(tsum - tsum0)/abs(tsum0)
         rels = abs(ssum - ssum0)/abs(ssum0)
         call check(error, relh < 1.0e-10_wp, "GM: column mass conserved through continuity")
         if (allocated(error)) exit checks
         call check(error, relt < 1.0e-10_wp, "GM: temperature content conserved")
         if (allocated(error)) exit checks
         call check(error, rels < 1.0e-10_wp, "GM: salinity content conserved")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_tracer_conserves

   ! ------------------------------------------------------------------
   ! Test 3b: no bolus leak across a no-normal-flow WALL.
   ! Regression for the GM/MLE wall-closure bug found by benchmark_ALE:
   ! the old fold added uhD/vhD at the physical WALL faces (which the
   ! resolved flux had zeroed), so the bolus carried tracer mass into the
   ! ghost halo and the PHYSICAL-domain sum(hTr) drifted (~1e-8/step).
   ! With no `bc` argument the GM operator treats every edge as a wall and
   ! zeroes the bolus there, so the physical-domain content must conserve
   ! to round-off.
   ! ------------------------------------------------------------------
   subroutine test_wall_no_leak(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 6
      real(wp), parameter :: DX = 2000.0_wp, DZ = 50.0_wp
      real(wp) :: tsum0, tsum, ssum0, ssum, relt, rels
      integer :: i0, i1, j0, j1
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         call make_eos(eos)
         call fill_tilted_TS(ms, grid, NZ, DZ, 1.0e-4_wp)

         ! PHYSICAL interior only (excludes the ghost halo).
         i0 = grid%nghost + 1
         i1 = grid%nghost + grid%nx_phys
         j0 = grid%nghost + 1
         j1 = grid%nghost + grid%ny_phys

         tsum0 = sum_phys_tr(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)
         ssum0 = sum_phys_tr(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)

         call map_in_ct(ms, sl, gm, ct)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call continuity_gm_apply(grid, metrics, ct, ms, gm, DT, 1.0_wp, TR_MODE_ADVECT)
         call map_out_ct(ms, sl, gm, ct)

         tsum = sum_phys_tr(ms, ms%idx_temperature, i0, i1, j0, j1, NZ)
         ssum = sum_phys_tr(ms, ms%idx_salinity, i0, i1, j0, j1, NZ)

         relt = abs(tsum - tsum0)/abs(tsum0)
         rels = abs(ssum - ssum0)/abs(ssum0)
         call check(error, relt < 1.0e-12_wp, "GM: no temperature leak across WALL")
         if (allocated(error)) exit checks
         call check(error, rels < 1.0e-12_wp, "GM: no salinity leak across WALL")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_wall_no_leak

   ! ------------------------------------------------------------------
   ! Test 6: the sequential operator cannot take a partial cell negative
   ! ------------------------------------------------------------------
   subroutine test_sequential_partial_cell(error)
      !! The 1-degree Southern Ocean failure, reduced to one column (z* open
      !! steps, nothing closed).  The centre column's bed layer is an 8.6 cm
      !! PARTIAL cell; in its four neighbours the same layer is an inert
      !! filler (1e-4 m).  The stored isopycnals are DOMED over the centre,
      !! so GM drains its deep water outward on all four faces — at the
      !! availability cap `A·(h − H_VANISHED)/(4·dt)` with this slope.  The
      !! dynamics of the step has already taken 60 % of the partial cell.
      !!
      !!   * FOLDED (the pre-2026-10 path): the transport computed from the
      !!     stage-ENTRY thickness (8.6 cm) and applied to what the dynamics
      !!     left (3.44 cm) — the cap bounds the wrong `h`, and the cell goes
      !!     NEGATIVE.  Asserted, as the witness that the case bites.
      !!   * SEQUENTIAL (now): the transport computed from the thickness the
      !!     dynamics left — the cap bounds what is there, the cell ends at
      !!     >= H_VANISHED, every layer stays >= min(its h, H_VANISHED), and
      !!     total volume and T/S content are conserved to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      integer, parameter :: NX = 5, NY = 5, NZ = 4
      real(wp), parameter :: DX = 2000.0_wp, HL = 1000.0_wp
      real(wp), parameter :: PARTIAL = 0.086_wp, FILL = 1.0e-4_wp, KEEP = 0.4_wp
      real(wp), parameter :: SLOPE0 = 5.0e-3_wp
      real(wp), allocatable :: h_left(:, :, :), h_entry(:, :, :)
      real(wp) :: vol0, vol, t0s, ts, s0s, ss, h_fold, h_seq, cap, gain
      integer :: i, j, k, ic, jc, ni, nj, it_t, it_s
      logical :: floor_ok
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         ni = grid%nx_total
         nj = grid%ny_total
         ic = NGHOST + (NX + 1)/2
         jc = NGHOST + (NY + 1)/2
         it_t = ms%idx_temperature
         it_s = ms%idx_salinity

         ! Stage-ENTRY thickness: uniform 1000 m layers, the centre bed layer
         ! an 8.6 cm partial cell, its four neighbours' bed layer a filler.
         allocate (h_entry(ni, nj, NZ))
         h_entry = HL
         h_entry(ic - 1, jc, 1) = FILL
         h_entry(ic + 1, jc, 1) = FILL
         h_entry(ic, jc - 1, 1) = FILL
         h_entry(ic, jc + 1, 1) = FILL
         h_entry(ic, jc, 1) = PARTIAL
         ! What the dynamics LEFT: 60 % of the partial cell moved into the
         ! four neighbouring fillers.
         allocate (h_left(ni, nj, NZ))
         h_left = h_entry
         gain = (1.0_wp - KEEP)*PARTIAL/4.0_wp
         h_left(ic, jc, 1) = KEEP*PARTIAL
         h_left(ic - 1, jc, 1) = FILL + gain
         h_left(ic + 1, jc, 1) = FILL + gain
         h_left(ic, jc - 1, 1) = FILL + gain
         h_left(ic, jc + 1, 1) = FILL + gain

         ! Domed stored slopes: + west/south of the centre, - east/north;
         ! zero at the bed (K=1) and the surface (K=NZ+1).
         sl%slope_x = 0.0_wp
         sl%slope_y = 0.0_wp
         sl%n2_u = 1.0e-6_wp
         sl%n2_v = 1.0e-6_wp
         do k = 2, NZ
            do j = 1, nj
               do i = 1, ni + 1
                  sl%slope_x(i, j, k) = merge(SLOPE0, -SLOPE0, i <= ic)
               end do
            end do
            do j = 1, nj + 1
               do i = 1, ni
                  sl%slope_y(i, j, k) = merge(SLOPE0, -SLOPE0, j <= jc)
               end do
            end do
         end do

         ! ---- FOLDED: transport from the stage-entry h, applied to h_left.
         call fill_layers(ms, h_entry, ni, nj, NZ)
         call map_in_ct(ms, sl, gm, ct)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         !$acc update self(gm%uhD)
         cap = DX*DX*(PARTIAL - H_VANISHED)/(4.0_wp*DT)
         call check(error, abs(gm%uhD(ic + 1, jc, 1) - cap) <= 1.0e-12_wp*cap .and. &
                    abs(-gm%uhD(ic, jc, 1) - cap) <= 1.0e-12_wp*cap, &
                    "the dome must drain the partial cell AT the cap (case strength)")
         if (allocated(error)) exit checks
         ms%h_layer = h_left
         !$acc update device(ms%h_layer)
         call continuity_gm_apply(grid, metrics, ct, ms, gm, DT, 1.0_wp, TR_MODE_ADVECT)
         !$acc update self(ms%h_layer)
         h_fold = ms%h_layer(ic, jc, 1)
         call map_out_ct(ms, sl, gm, ct)
         call check(error, h_fold < 0.0_wp, &
                    "witness: the stage-entry (folded) transport drives the partial cell negative")
         if (allocated(error)) exit checks

         ! ---- SEQUENTIAL: transport from the thickness the dynamics left.
         call fill_layers(ms, h_left, ni, nj, NZ)
         vol0 = sum_phys_h(ms, 1, ni, 1, nj, NZ)
         t0s = sum_phys_tr(ms, it_t, 1, ni, 1, nj, NZ)
         s0s = sum_phys_tr(ms, it_s, 1, ni, 1, nj, NZ)
         call map_in_ct(ms, sl, gm, ct)
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
         call continuity_gm_apply(grid, metrics, ct, ms, gm, DT, 1.0_wp, TR_MODE_ADVECT)
         call map_out_ct(ms, sl, gm, ct)
         h_seq = ms%h_layer(ic, jc, 1)
         floor_ok = .true.
         do k = 1, NZ
            do j = 1, nj
               do i = 1, ni
                  if (ms%h_layer(i, j, k) < min(h_left(i, j, k), H_VANISHED)*(1.0_wp - 1.0e-12_wp)) &
                     floor_ok = .false.
               end do
            end do
         end do
         vol = sum_phys_h(ms, 1, ni, 1, nj, NZ)
         ts = sum_phys_tr(ms, it_t, 1, ni, 1, nj, NZ)
         ss = sum_phys_tr(ms, it_s, 1, ni, 1, nj, NZ)

         call check(error, h_seq >= H_VANISHED*(1.0_wp - 1.0e-12_wp), &
                    "sequential GM must leave the partial cell >= H_VANISHED")
         if (allocated(error)) exit checks
         call check(error, h_seq < KEEP*PARTIAL, &
                    "sequential GM must still drain the partial cell (non-trivial)")
         if (allocated(error)) exit checks
         call check(error, floor_ok, &
                    "sequential GM must keep every layer >= min(h, H_VANISHED)")
         if (allocated(error)) exit checks
         call check(error, abs(vol - vol0) <= 1.0e-14_wp*vol0, &
                    "sequential GM must conserve volume")
         if (allocated(error)) exit checks
         call check(error, abs(ts - t0s) <= 1.0e-13_wp*abs(t0s), &
                    "sequential GM must conserve temperature content")
         if (allocated(error)) exit checks
         call check(error, abs(ss - s0s) <= 1.0e-13_wp*abs(s0s), &
                    "sequential GM must conserve salinity content")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_sequential_partial_cell

   subroutine fill_layers(ms, h, ni, nj, nz)
      !! Set `h_layer` and T/S content (T warm over cold, S uniform; a
      !! filler carries its donor's concentration, I1').
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: ni, nj, nz
      real(wp), intent(in) :: h(ni, nj, nz)
      integer :: k
      ms%h_layer = h
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (2.0_wp + 4.0_wp*real(k, wp))*h(:, :, k)
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S0*h(:, :, k)
      end do
      ! I1': the bed fillers carry the layer above's temperature.
      where (h(:, :, 1) <= H_VANISHED) &
         ms%tracers(ms%idx_temperature)%hTr(:, :, 1) = (2.0_wp + 8.0_wp)*h(:, :, 1)
   end subroutine fill_layers

   ! ------------------------------------------------------------------
   ! continuity-aware map helpers + physical-domain sums
   ! ------------------------------------------------------------------
   subroutine map_in_ct(ms, sl, gm, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      ! Bed datum of the slopes' geopotential interface heights: every
      ! case here has a flat bed under a constant-depth column at rest
      ! (eta = 0), so D = Sum_k h.  Host-side, BEFORE the slot's map.
      call sl%set_bathymetry(sum(ms%h_layer, dim=3))
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in_ct

   subroutine map_out_ct(ms, sl, gm, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(continuity_t), intent(inout) :: ct
      !$acc update self(ms%h_layer)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_ct

   function sum_phys_h(ms, i0, i1, j0, j1, nz) result(s)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: i0, i1, j0, j1, nz
      real(wp) :: s
      integer :: i, j, k
      s = 0.0_wp
      do k = 1, nz
         do j = j0, j1
            do i = i0, i1
               s = s + ms%h_layer(i, j, k)
            end do
         end do
      end do
   end function sum_phys_h

   function sum_phys_tr(ms, idx, i0, i1, j0, j1, nz) result(s)
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx, i0, i1, j0, j1, nz
      real(wp) :: s
      integer :: i, j, k
      s = 0.0_wp
      do k = 1, nz
         do j = j0, j1
            do i = i0, i1
               s = s + ms%tracers(idx)%hTr(i, j, k)
            end do
         end do
      end do
   end function sum_phys_tr

end module test_ocean_gm
