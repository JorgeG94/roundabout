!! Unit tests for the ocean vertical-advection kernel
!! (rdb_ocean_vertical_advection).  Two pieces:
!!
!!   1. `compute_w_from_continuity` — diagnoses
!!      `ms%w_interface` from `ms%flux_h_layer` by integrating
!!      upward under the bed BC `w(k=1) = 0`.
!!   2. `tracer_advect_vertical` — first-order upwind
!!      vertical tracer advection consuming `ms%w_interface`.
!!
!! Cases:
!!   * Constancy — uniform tracer stays uniform under any w field.
!!     The CWC-like preservation property in the vertical.
!!   * Zero w → no-op — `w_interface = 0` everywhere leaves tracers
!!     bit-for-bit at IC.
!!   * Conservation — sum(hTr) is preserved by vertical advection
!!     across the column (vertical fluxes telescope: F(k) - F(k+1)
!!     summed over k gives 0 because bed and surface fluxes are
!!     forced to zero by the kernel).
!!   * Downward transport — prescribe uniform w < 0 in the interior
!!     (with bed/surface zero).  A surface step in T descends into
!!     the layers below; net mass moves downward.
!!   * w-from-continuity diagnostic — set a known `flux_h_layer`
!!     pattern, check the integrated `w_interface` matches the
!!     analytic cumulative sum.
module test_ocean_vertical_advection
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vertical_advection, only: &
      ocean_vertical_advection_t, &
      compute_w_from_continuity, &
      tracer_advect_vertical
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_HEAT_TOTAL, BUDGET_SALT_TOTAL
   implicit none
   private

   public :: collect_ocean_vertical_advection_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_vertical_advection_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("vert_adv_constancy", test_constancy), &
                  new_unittest("vert_adv_zero_w_noop", test_zero_w_noop), &
                  new_unittest("vert_adv_conservation", test_conservation), &
                  new_unittest("vert_adv_downward_transport", test_downward), &
                  new_unittest("compute_w_from_continuity", test_w_from_continuity), &
                  new_unittest("vert_adv_budget_contributor_telescopes", &
                               test_vert_adv_budget_contributor) &
                  ]
   end subroutine collect_ocean_vertical_advection_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms, vadv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vertical_advection_t), intent(inout) :: vadv
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vadv)
      call vadv%enter_data()
   end subroutine map_in

   subroutine map_out(ms, vadv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vertical_advection_t), intent(inout) :: vadv
      call vadv%exit_data()
      !$acc exit data delete(vadv)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_constancy(error)
      !! Uniform T across all (i, j, k).  Any w field should leave
      !! T uniform.  Strongest discriminator for vertical-CWC bugs.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: T0 = 25.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: max_T_dev, T_obs
      integer :: i, j, k, nx, ny

      call make_grid(grid, 8, 6)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vadv%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = T0*H0
      ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*H0

      ! Prescribe a non-trivial w (interior + bed/surface zero)
      ms%w_interface = 0.0_wp
      do k = 2, NZ
         do j = 1, ny
            do i = 1, nx
               ms%w_interface(i, j, k) = 0.1_wp* &
                                         sin(PI*real(k - 1, wp)/real(NZ, wp))
            end do
         end do
      end do

      call map_in(ms, vadv)
      call tracer_advect_vertical(grid, vadv, ms, DT)
      call map_out(ms, vadv)

      ! h_layer has now been updated by the vertical mass divergence,
      ! so recover T from the post-update hTr/h pair.  CWC says T
      ! must remain uniform = T0 regardless of how non-uniform w made
      ! the layer thicknesses.
      max_T_dev = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (ms%h_layer(i, j, k) > 0.0_wp) then
                  T_obs = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ &
                          ms%h_layer(i, j, k)
                  max_T_dev = max(max_T_dev, abs(T_obs - T0))
               end if
            end do
         end do
      end do
      call check(error, max_T_dev < 1.0e-12_wp, &
                 "vert advect: uniform T not preserved under w divergence")

      call vadv%destroy(); call ms%destroy()
   end subroutine test_constancy

   subroutine test_zero_w_noop(error)
      !! w_interface = 0 everywhere ⇒ no flux ⇒ tracers unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), allocatable :: hT_ic(:, :, :)
      real(wp) :: max_diff
      integer :: i, j, k, nx, ny

      call make_grid(grid, 8, 6)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vadv%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
            real(k, wp)*H0
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do
      allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

      ms%w_interface = 0.0_wp

      call map_in(ms, vadv)
      call tracer_advect_vertical(grid, vadv, ms, DT)
      call map_out(ms, vadv)

      max_diff = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))
      call check(error, max_diff < 1.0e-12_wp, &
                 "vert advect: w=0 not a no-op")

      deallocate (hT_ic)
      call vadv%destroy(); call ms%destroy()
   end subroutine test_zero_w_noop

   subroutine test_conservation(error)
      !! Column-integrated tracer mass `sum_k hTr(:,:,k)` must be
      !! preserved.  Per (i, j) the vertical-flux contributions
      !! telescope: ∑_k (F(k) - F(k+1)) = F(1) - F(nz+1) = 0 by
      !! the bed/surface BC.  Test on a non-trivial stratification
      !! with a non-trivial w field.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: total_before, total_after, drift
      integer :: i, j, k, nx, ny

      call make_grid(grid, 10, 8)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vadv%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  (10.0_wp + 2.0_wp*real(k, wp) + &
                   0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))*H0
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 35.0_wp*H0
            end do
         end do
      end do

      ms%w_interface = 0.0_wp
      do k = 2, NZ
         do j = 1, ny
            do i = 1, nx
               ms%w_interface(i, j, k) = 0.05_wp* &
                                         sin(PI*real(k - 1, wp)/real(NZ, wp))* &
                                         cos(2.0_wp*PI*real(i, wp)/real(nx, wp))
            end do
         end do
      end do

      total_before = sum(ms%tracers(ms%idx_temperature)%hTr)
      call map_in(ms, vadv)
      call tracer_advect_vertical(grid, vadv, ms, DT)
      call map_out(ms, vadv)
      total_after = sum(ms%tracers(ms%idx_temperature)%hTr)

      drift = abs(total_after - total_before)/abs(total_before)
      call check(error, drift < 1.0e-12_wp, &
                 "vert advect: total tracer mass drifted > 1e-12")

      call vadv%destroy(); call ms%destroy()
   end subroutine test_conservation

   subroutine test_downward(error)
      !! Step IC in T: surface layer hot, lower layers cool.  Apply
      !! a uniform downward w in the interior interfaces.  After
      !! one step the surface-layer tracer must decrease and the
      !! layer-below tracer must increase (heat moves down).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: W_DOWN = -0.5_wp
      real(wp) :: hT_surf_before, hT_surf_after
      real(wp) :: hT_below_before, hT_below_after
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vadv%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ms%tracers(ms%idx_temperature)%hTr = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*H0
         ! Step: surface layer hot, others cold
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = 20.0_wp*H0

         ms%w_interface = 0.0_wp
         do k = 2, NZ
            ms%w_interface(:, :, k) = W_DOWN
         end do

         hT_surf_before = ms%tracers(ms%idx_temperature)%hTr(nx/2, ny/2, NZ)
         hT_below_before = ms%tracers(ms%idx_temperature)%hTr(nx/2, ny/2, NZ - 1)

         call map_in(ms, vadv)
         call tracer_advect_vertical(grid, vadv, ms, DT)
         call map_out(ms, vadv)

         hT_surf_after = ms%tracers(ms%idx_temperature)%hTr(nx/2, ny/2, NZ)
         hT_below_after = ms%tracers(ms%idx_temperature)%hTr(nx/2, ny/2, NZ - 1)

         call check(error, hT_surf_after < hT_surf_before, &
                    "downward w: surface tracer did not decrease")
         if (allocated(error)) exit checks
         call check(error, hT_below_after > hT_below_before, &
                    "downward w: layer-below tracer did not increase")

      end block checks
      call vadv%destroy(); call ms%destroy()
   end subroutine test_downward

   subroutine test_w_from_continuity(error)
      !! Prescribe a known `flux_h_layer` pattern and check the
      !! diagnosed `w_interface` matches the analytic cumulative
      !! integral:
      !!   w_interface(k+1) = w_interface(k) - flux_h_layer(k)
      !!   w_interface(1)   = 0
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp) :: w_expected, w_obs, max_diff
      real(wp) :: flux_k(NZ), cum_sum
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vadv%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ! Layer-dependent constant flux divergence pattern.
         flux_k = [0.5_wp, -0.2_wp, 0.3_wp, -0.6_wp]
         do k = 1, NZ
            ms%flux_h_layer(:, :, k) = flux_k(k)
         end do
         ms%w_interface = 0.0_wp

         call map_in(ms, vadv)
         ! flux_h_layer is `create`-mapped (it's scratch in the
         ! normal pipeline — continuity writes it on-device); push
         ! the prescribed host values to device so the kernel reads
         ! them, not the initial uninitialised buffer.
         !$acc update device(ms%flux_h_layer)
         !$acc enter data copyin(vadv)
         call compute_w_from_continuity(grid, vadv, ms)
         !$acc exit data delete(vadv)
         call map_out(ms, vadv)

         max_diff = 0.0_wp
         cum_sum = 0.0_wp
         do k = 1, NZ
            ! w(k+1) = w(k) - flux_k(k), w(1) = 0
            cum_sum = cum_sum - flux_k(k)
            w_expected = cum_sum
            w_obs = ms%w_interface(nx/2, ny/2, k + 1)
            max_diff = max(max_diff, abs(w_obs - w_expected))
         end do
         call check(error, max_diff < 1.0e-12_wp, &
                    "w-from-continuity: diagnosed w off analytic")
         if (allocated(error)) exit checks
         call check(error, abs(ms%w_interface(nx/2, ny/2, 1)) < 1.0e-12_wp, &
                    "w-from-continuity: bed BC w(k=1) != 0")

      end block checks
      call vadv%destroy(); call ms%destroy()
   end subroutine test_w_from_continuity

   subroutine test_vert_adv_budget_contributor(error)
      !! Phase D v2 tracer-vert-adv kernel patch: the heat + salt
      !! contributor slots must telescope to ~0 over the spatial
      !! integral for closed BC vertical advection (bed and surface
      !! fluxes zero ⇒ Σ_k dt · (F(k) - F(k+1)) = 0 per column).
      !! Mirrors `test_conservation` but registers contributors on a
      !! local `ocean_budgets_t`; drains AFTER advection (with device
      !! data still mapped) and asserts contributor totals are bounded
      !! by the total tracer mass round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vertical_advection_t) :: vadv
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: total_T_before, total_S_before, total_T_after, total_S_after
      real(wp) :: rel_T, rel_S
      integer :: i, j, k, nx, ny, idx_heat, idx_salt
      checks: block

         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vadv%init(grid, nz_ml=NZ)
         call budgets%init(grid)
         call budgets%register_contributor("vert_adv_heat", BUDGET_HEAT_TOTAL, &
                                           ms%heat_budget_vert_adv, &
                                           device_resident=.true.)
         idx_heat = budgets%n_contributors
         call budgets%register_contributor("vert_adv_salt", BUDGET_SALT_TOTAL, &
                                           ms%salt_budget_vert_adv, &
                                           device_resident=.true.)
         idx_salt = budgets%n_contributors

         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = H0
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     (10.0_wp + 2.0_wp*real(k, wp) + &
                      0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp)))*H0
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (35.0_wp + 0.3_wp*cos(2.0_wp*PI*real(j, wp)/real(ny, wp)))*H0
               end do
            end do
         end do
         ms%w_interface = 0.0_wp
         do k = 2, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%w_interface(i, j, k) = 0.05_wp* &
                                            sin(PI*real(k - 1, wp)/real(NZ, wp))* &
                                            cos(2.0_wp*PI*real(i, wp)/real(nx, wp))
               end do
            end do
         end do

         total_T_before = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_before = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         call map_in(ms, vadv)
         call tracer_advect_vertical(grid, vadv, ms, DT)
         ! Drain BEFORE map_out so the contributor's device buffer is
         ! still mapped (mirror of the surface-flux test).
         call budgets%drain_contributors()
         call map_out(ms, vadv)

         total_T_after = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_after = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         ! Per-column telescope: bed + surface flux = 0 ⇒ contributor
         ! sums to zero.  Compare to the absolute mass for tolerance.
         rel_T = abs(budgets%contributors(idx_heat)%total_integrated)/abs(total_T_before)
         rel_S = abs(budgets%contributors(idx_salt)%total_integrated)/abs(total_S_before)
         call check(error, rel_T < 1.0e-12_wp, &
                    "vert-adv heat contributor should telescope to FP")
         if (allocated(error)) exit checks
         call check(error, rel_S < 1.0e-12_wp, &
                    "vert-adv salt contributor should telescope to FP")
         if (allocated(error)) exit checks

         ! Closure: LHS (total drift) = RHS (contributor).  Vert-adv is
         ! the only kernel touching hTr here, so LHS = RHS to FP.
         call check(error, abs((total_T_after - total_T_before) - &
                               budgets%contributors(idx_heat)%total_integrated) &
                    < 1.0e-12_wp*abs(total_T_before), &
                    "heat LHS = RHS for closed vert-adv")
         if (allocated(error)) exit checks
         call check(error, abs((total_S_after - total_S_before) - &
                               budgets%contributors(idx_salt)%total_integrated) &
                    < 1.0e-12_wp*abs(total_S_before), &
                    "salt LHS = RHS for closed vert-adv")
      end block checks
      call budgets%destroy()
      call vadv%destroy(); call ms%destroy()
   end subroutine test_vert_adv_budget_contributor

end module test_ocean_vertical_advection
