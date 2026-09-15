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
!!  3. gm_tracer_conserves   — column tracer content conserved through the
!!                             fold + continuity drain.
!!  4. gm_slope_limiter      — S >> slope_max ⇒ bounded Psi, h >= H_VANISHED.
!!  5. gm_gm_src_sign        — gm_src >= 0 (PE release) for a stable column.
module test_ocean_gm
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports, gm_fold_x, gm_fold_y
   use rdb_continuity, only: continuity_t, continuity_tracer_step_split
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
                  new_unittest("gm_gm_src_sign", test_gm_src_sign) &
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
            call continuity_tracer_step_split(grid, metrics, ct, ms, DTL, gm=gm)
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
   ! Test 3: tracer content conserved through fold + continuity
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

         ! GM is conservative over the CLOSED discrete domain.  Bolus flux
         ! vanishes only at the ARRAY edges (i=1/nx+1, j=1/ny+1), where the
         ! slopes are zeroed; the physical sub-region exchanges tracer with
         ! the ghost halo across its inner faces (the tilted-T edge fluxes
         ! don't cancel even though uniform-h mass flux nets to zero).  So
         ! the conservation invariant is the FULL-array sum.
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
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT, gm=gm)
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
   ! the fold added uhD/vhD at the physical WALL faces (which the resolved
   ! flux had zeroed), so the bolus carried tracer mass into the ghost halo
   ! and the PHYSICAL-domain sum(hTr) drifted (~1e-8/step) even though the
   ! FULL-array sum stayed closed.  With no `bc` argument the continuity
   ! step defaults every edge to OBC_WALL, so the physical-domain content
   ! must now conserve to round-off.  Resolved velocity is zero (tilted-T
   ! only) ⇒ the bolus fold is the SOLE transport, isolating the wall path.
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
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT, gm=gm)
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
   ! continuity-aware map helpers + physical-domain sums
   ! ------------------------------------------------------------------
   subroutine map_in_ct(ms, sl, gm, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
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
