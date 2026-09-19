!! Unit tests for the FV_WRIGHT variant of the ocean pressure-force
!! kernel (rdb_ocean_pressure_force with
!! `variant = OPGF_VARIANT_FV_WRIGHT`).
!!
!! FV_WRIGHT layers full Wright-along-the-path density on top of the
!! FV_LITE z-correction.  The column sweep helper does a single
!! Picard step using `ms%rho_layer` (= Wright at p_ref=0) as the
!! half-layer pressure seed, then re-runs Wright at that seed to
!! produce `pgf%rho_insitu`.  Pressure stack + `rho_face` averaging
!! both read from the in-situ field rather than from `ms%rho_layer`.
!!
!! Cases:
!!   * No lateral gradient → no PGF, regardless of column depth.
!!     Sanity check that compressibility along the path doesn't
!!     manufacture a horizontal acceleration when none should exist.
!!   * Column sweep recovers Wright(T, S, p_centre) — direct unit
!!     test on the flat-impl: set up a single-column state, run the
!!     sweep, compare rho_insitu against the analytic Wright eval
!!     at the corresponding half-layer pressure.
!!   * Thermobaric amplification at depth — lateral T gradient over
!!     a deep aligned-column setup.  Cold water compresses more than
!!     warm at depth, so the bed-layer dpdx under FV_WRIGHT must
!!     exceed the FV_LITE response by a detectable margin.
module test_ocean_pgf_wright
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_wright_pgf_column_sweep_impl
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_LITE, OPGF_VARIANT_FV_WRIGHT
   implicit none
   private

   public :: collect_ocean_pgf_wright_tests

   integer, parameter :: NGHOST = 2

   ! Wright (1997) coefficients, mirrored from rdb_eos for the
   ! direct sweep unit test.  Same values; copied locally so the test
   ! doesn't depend on the EOS module exposing them.
   real(wp), parameter :: WA0 = 7.057924e-4_wp, WA1 = 3.480336e-7_wp, WA2 = -1.112733e-7_wp
   real(wp), parameter :: WB0 = 5.790749e8_wp, WB1 = 3.516535e6_wp
   real(wp), parameter :: WB2 = -4.002714e4_wp, WB3 = 2.084372e2_wp
   real(wp), parameter :: WB4 = 5.944068e5_wp, WB5 = -9.643486e3_wp
   real(wp), parameter :: WC0 = 1.704853e5_wp, WC1 = 7.904722e2_wp
   real(wp), parameter :: WC2 = -7.984422e0_wp, WC3 = 5.140652e-2_wp
   real(wp), parameter :: WC4 = -2.302158e2_wp, WC5 = -3.079464e0_wp

contains

   subroutine collect_ocean_pgf_wright_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("wright_no_lateral_gradient_no_pgf", &
                               test_no_lateral_no_pgf), &
                  new_unittest("wright_column_sweep_recovers_wright_at_pressure", &
                               test_column_sweep_unit), &
                  new_unittest("wright_thermobaric_amplification_at_depth", &
                               test_thermobaric_amplification) &
                  ]
   end subroutine collect_ocean_pgf_wright_tests

   pure function wright_rho_at(T, S, P) result(rho)
      !! Reference Wright eval for the unit-test analytic.
      real(wp), intent(in) :: T, S, P
      real(wp) :: rho
      real(wp) :: T_sq, T_cu, alpha_0, p_0, lambda, p_plus_p0
      T_sq = T*T
      T_cu = T_sq*T
      alpha_0 = WA0 + WA1*T + WA2*S
      p_0 = WB0 + WB1*T + WB2*T_sq + WB3*T_cu + WB4*S + WB5*S*T
      lambda = WC0 + WC1*T + WC2*T_sq + WC3*T_cu + WC4*S + WC5*S*T
      p_plus_p0 = P + p_0
      rho = p_plus_p0/(lambda + alpha_0*p_plus_p0)
   end function wright_rho_at

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine run_pgf(grid, ms, pgf)
      !! Map ms + pgf to the device, run compute, pull dpdx/dpdy and
      !! rho_insitu back to the host for inspection.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_metrics_t) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data, pgf%rho_insitu%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_no_lateral_no_pgf(error)
      !! Uniform T, S, uniform h_layer, deep column.  No lateral
      !! density variation → both FV_LITE and FV_WRIGHT must produce
      !! a zero PGF.  Compressibility along the column does change
      !! `rho_insitu` away from `rho_layer`, but with no x-dependence
      !! the horizontal gradient is identically zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      integer, parameter :: NZ = 5
      real(wp), parameter :: H_LAYER = 200.0_wp        ! 5 * 200 m = 1 km deep
      real(wp), parameter :: T_UNI = 10.0_wp, S_UNI = 35.0_wp
      real(wp) :: max_dpdx, max_dpdy, max_rho_insitu
      checks: block

         call make_grid(grid, 12, 10, 1.0e3_wp, 1.0e3_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_FV_WRIGHT

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = S_UNI*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = T_UNI*H_LAYER
         ! Seed rho_layer with the surface Wright value so the Picard
         ! seed is sane.
         ms%rho_layer = wright_rho_at(T_UNI, S_UNI, 0.0_wp)

         call run_pgf(grid, ms, pgf)

         max_dpdx = maxval(abs(pgf%dpdx_face%data))
         max_dpdy = maxval(abs(pgf%dpdy_face%data))
         max_rho_insitu = maxval(pgf%rho_insitu%data)

         call check(error, max_dpdx < 1.0e-12_wp, &
                    "uniform T,S: FV_WRIGHT produced spurious dpdx")
         if (allocated(error)) exit checks
         call check(error, max_dpdy < 1.0e-12_wp, &
                    "uniform T,S: FV_WRIGHT produced spurious dpdy")
         if (allocated(error)) exit checks
         ! Sanity: rho_insitu should be larger than surface Wright due to
         ! compressibility (we have ~1 km of water).
         call check(error, max_rho_insitu > ms%rho_layer(1, 1, 1), &
                    "FV_WRIGHT rho_insitu didn't pick up any compressibility")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_no_lateral_no_pgf

   subroutine test_column_sweep_unit(error)
      !! Direct unit test on the column-sweep flat-impl.  Build a
      !! tiny multilayer state, call the sweep, compare each layer's
      !! `rho_insitu(k)` against the analytic Wright eval at the
      !! corresponding half-layer pressure (`p_above + 0.5*g*rho_seed*h`).
      !!
      !! Both `p_edge` and `rho_insitu` are checked column-by-column.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 3, NZ = 4
      real(wp), parameter :: H = 250.0_wp           ! 4*250 m = 1 km
      real(wp), parameter :: T_UNI = 8.0_wp, S_UNI = 34.5_wp
      real(wp), parameter :: RHO_SEED = 1027.5_wp   ! a stand-in for surface ρ
      real(wp), parameter :: RHO_0 = 1035.0_wp
      real(wp), parameter :: TOL_RHO = 5.0e-3_wp   ! 0.005 kg/m^3 tolerance vs analytic
      real(wp), parameter :: TOL_PEDGE = 50.0_wp   ! 50 Pa over ~1e7 Pa is < 1e-5 relative
      real(wp), allocatable :: h_layer(:, :, :), hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: rho_layer_seed(:, :, :)
      real(wp), allocatable :: p_top(:, :)
      real(wp), allocatable :: p_edge_out(:, :, :), rho_insitu_out(:, :, :)
      real(wp) :: p_above, p_centre_seed, rho_analytic
      real(wp) :: max_rho_err, max_pedge_err, p_edge_expected
      integer :: i, j, k
      checks: block

         allocate (h_layer(NX, NY, NZ), source=H)
         allocate (hS(NX, NY, NZ), source=S_UNI*H)
         allocate (hT(NX, NY, NZ), source=T_UNI*H)
         allocate (rho_layer_seed(NX, NY, NZ), source=RHO_SEED)
         allocate (p_top(NX, NY), source=0.0_wp)
         allocate (p_edge_out(NX, NY, NZ + 1))
         allocate (rho_insitu_out(NX, NY, NZ))

         !$acc enter data copyin(h_layer, hS, hT, rho_layer_seed, p_top)
         !$acc enter data create(p_edge_out, rho_insitu_out)
         call eos_wright_pgf_column_sweep_impl( &
            h_layer, hS, hT, rho_layer_seed, p_top, &
            p_edge_out, rho_insitu_out, &
            GRAVITY, RHO_0, NX, NY, NZ)
         !$acc update self(p_edge_out, rho_insitu_out)
         !$acc exit data delete(p_edge_out, rho_insitu_out)
         !$acc exit data delete(h_layer, hS, hT, rho_layer_seed, p_top)

         max_rho_err = 0.0_wp
         max_pedge_err = 0.0_wp
         do j = 1, NY
            do i = 1, NX
               p_above = 0.0_wp
               do k = NZ, 1, -1
                  p_centre_seed = p_above + 0.5_wp*GRAVITY*RHO_SEED*H
                  rho_analytic = wright_rho_at(T_UNI, S_UNI, p_centre_seed)
                  max_rho_err = max(max_rho_err, &
                                    abs(rho_insitu_out(i, j, k) - rho_analytic))
                  p_above = p_above + GRAVITY*rho_analytic*H
                  p_edge_expected = p_above
                  max_pedge_err = max(max_pedge_err, &
                                      abs(p_edge_out(i, j, k) - p_edge_expected))
               end do
            end do
         end do

         call check(error, max_rho_err < TOL_RHO, &
                    "column sweep rho_insitu doesn't match Wright at p_centre")
         if (allocated(error)) exit checks
         call check(error, max_pedge_err < TOL_PEDGE, &
                    "column sweep p_edge accumulation off")
         if (allocated(error)) exit checks
         call check(error, abs(p_edge_out(1, 1, NZ + 1)) < 1.0e-12_wp, &
                    "surface p_edge should be exactly zero")

      end block checks
      deallocate (h_layer, hS, hT, rho_layer_seed, p_edge_out, rho_insitu_out)
   end subroutine test_column_sweep_unit

   subroutine test_thermobaric_amplification(error)
      !! Deep aligned-column setup with a lateral T gradient and
      !! S = 35 throughout.  Cold water compresses more than warm at
      !! pressure (Wright's `a1` and `c1` terms shift `α_0` and `λ`
      !! together), so the lateral density gradient is amplified
      !! near the bed under FV_WRIGHT relative to FV_LITE.
      !!
      !! Both variants run on the SAME `ms%rho_layer` (surface Wright
      !! at p_ref=0).  FV_LITE keeps using that throughout; FV_WRIGHT
      !! re-evaluates at the in-situ pressure.  We compare the
      !! bed-layer (k=1) east-face dpdx between the two runs:
      !! `max |dpdx_wright[k=1]| > max |dpdx_lite[k=1]|` by a clear
      !! margin (we ask for ≥ 5% amplification).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_lite, pgf_wright
      integer, parameter :: NZ = 5
      integer, parameter :: NX_PHYS = 24, NY_PHYS = 8
      real(wp), parameter :: H_LAYER = 800.0_wp     ! 5 * 800 m = 4 km deep
      real(wp), parameter :: DX = 5.0e3_wp, DY = 5.0e3_wp
      real(wp), parameter :: S_UNI = 35.0_wp
      real(wp), parameter :: T_REF = 6.0_wp, T_GRAD = 0.5_wp ! 0.5 °C / cell
      real(wp) :: T_at_i, max_lite, max_wright
      real(wp), allocatable :: dpdx_lite_bed(:, :)
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_lite%init(grid, nz_ml=NZ)
         call pgf_wright%init(grid, nz_ml=NZ)
         pgf_lite%variant = OPGF_VARIANT_FV_LITE
         pgf_wright%variant = OPGF_VARIANT_FV_WRIGHT
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H_LAYER

         ! Lateral T gradient: T(i, j, k) = T_REF + T_GRAD * (i - nx/2).
         ! Uniform in y and k.  Aligned columns so z_R == z_L; any PGF
         ! comes from the lateral ρ gradient (not from the z-correction).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  T_at_i = T_REF + T_GRAD*real(i - nx/2, wp)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S_UNI*H_LAYER
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = T_at_i*H_LAYER
                  ms%rho_layer(i, j, k) = wright_rho_at(T_at_i, S_UNI, 0.0_wp)
               end do
            end do
         end do

         call run_pgf(grid, ms, pgf_lite)
         allocate (dpdx_lite_bed, source=pgf_lite%dpdx_face%data(:, :, 1))
         max_lite = maxval(abs(dpdx_lite_bed))

         call run_pgf(grid, ms, pgf_wright)
         max_wright = maxval(abs(pgf_wright%dpdx_face%data(:, :, 1)))

         ! FV_LITE responds to the surface-ρ gradient; FV_WRIGHT
         ! responds to the in-situ ρ gradient which is amplified by
         ! compressibility at the bed.  Demand at least 5% larger.
         call check(error, max_wright > 1.05_wp*max_lite, &
                    "FV_WRIGHT failed to amplify bed-layer PGF relative to FV_LITE")
         if (allocated(error)) exit checks
         ! Sanity that the FV_LITE run already produced a non-trivial signal
         call check(error, max_lite > 1.0e-6_wp, &
                    "FV_LITE bed-layer PGF too small — test setup didn't engage")

      end block checks
      deallocate (dpdx_lite_bed)
      call pgf_wright%destroy(); call pgf_lite%destroy(); call ms%destroy()
   end subroutine test_thermobaric_amplification

end module test_ocean_pgf_wright
