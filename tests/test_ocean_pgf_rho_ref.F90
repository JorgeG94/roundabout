!! The PGF reference densities must follow the SINGLE configured ρ₀.
!!
!! The defect these pin
!! --------------------
!! `ocean_pressure_force_t` carries two reference densities — `rho0`, the
!! Boussinesq divisor in `du/dt = −(1/ρ₀)·∂p/∂x`, and `rho_ref`, the baseline
!! subtracted from layer densities when building the FV_MOM6 `pa` anomaly
!! stack (also the `g·ρ_ref/ρ₀` surface value in `compute_pbce`).  Both were
!! declared `= 1035.0_wp` and `configure_ocean_pgf` assigned NEITHER, while
!! the EOS took its ρ₀ from `&ocean_ic_nml rho_0`.  A namelist setting
!! `rho_0 /= 1035` therefore ran an equation of state and a pressure gradient
!! on two different reference densities, silently — no warning, no fail-loud,
!! just a pressure gradient scaled by `1035/ρ₀_configured`.
!!
!! Cases
!! -----
!!   * `configured_rho0_reaches_both_pgf_reference_densities` — the plumbing
!!     proper.  A config with `rho_0 = 1000` (3.5 % off the default, so a
!!     dropped assignment cannot hide in round-off) is driven through the
!!     production path (`init_from_config` → `configure_ocean_pgf`), and both
!!     `rho0` and `rho_ref` are asserted equal to `eos%rho0` — the single ρ₀
!!     of record.  NON-VACUITY: the configured value is asserted to differ
!!     from the type default first, so the check cannot pass by both sides
!!     being 1035.
!!   * `two_column_pgf_answer_uses_the_configured_rho0` — the analytical arm.
!!     The SAME configured state is given a two-column hydrostatic field
!!     (flat bed, every column the same layer thickness so the `z_centre`
!!     correction vanishes identically, a horizontal density contrast across
!!     each face) and run through `ocean_pressure_force_compute`.  With
!!     uniform thickness the FV_LITE face expression collapses to the exact
!!     two-column hydrostatic form
!!
!!         dpdx(i,k) = −(1/ρ₀) · (p_c(i,k) − p_c(i−1,k)) / dx
!!         p_c(k) = ½(p_edge(k) + p_edge(k+1)),  p_edge(nz+1) = 0,
!!         p_edge(k) = p_edge(k+1) + g·ρ(k)·h(k)
!!
!!     whose ONLY ρ₀ dependence is the `1/ρ₀` prefactor — `p_edge` is built
!!     from `ρ(k)·h(k)` alone.  The kernel answer is asserted against that
!!     closed form evaluated at the CONFIGURED ρ₀, to round-off.
!!     NON-VACUITY: the same closed form evaluated at the old hard-coded
!!     1035 is asserted to be resolvably different, so the tolerance is
!!     proven to discriminate the two.
module test_ocean_pgf_rho_ref
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_pgf
   use rdb_ocean_pressure_force, only: ocean_pressure_force_compute
   implicit none
   private

   public :: collect_ocean_pgf_rho_ref_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 10000.0_wp
   real(wp), parameter :: H_LAYER = 100.0_wp
      !! Uniform in every column AND every layer — that is what makes the
      !! FV_LITE `z_correction` vanish exactly (both columns share a
      !! `z_centre`), leaving the pure two-column hydrostatic expression.
   real(wp), parameter :: RHO0_DEFAULT = 1035.0_wp
      !! The literal the PGF slot defaults to, and the `&ocean_ic_nml rho_0`
      !! default.  Named here only so the non-vacuity arms can show the test
      !! is not comparing 1035 against 1035.
   real(wp), parameter :: RHO0_CFG = 1000.0_wp
      !! Deliberately NOT the default: 3.5 % apart, far above any tolerance
      !! used below.
   real(wp), parameter :: RHO_BASE = 1028.0_wp
   real(wp), parameter :: RHO_STRAT = 0.4_wp
      !! Per-layer density increment towards the bed (k=1 is the bed).
   real(wp), parameter :: RHO_STEP = 0.6_wp
      !! Total west-to-east density contrast — what actually drives the
      !! gradient; without it every face answer would be zero and both arms
      !! would pass vacuously.

contains

   subroutine collect_ocean_pgf_rho_ref_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("configured_rho0_reaches_both_pgf_reference_densities", test_rho_ref_plumbing), &
                  new_unittest("two_column_pgf_answer_uses_the_configured_rho0", test_rho_ref_answer) &
                  ]
   end subroutine collect_ocean_pgf_rho_ref_tests

   ! ---------------------------------------------------------------------
   ! Scaffolding
   ! ---------------------------------------------------------------------

   subroutine build_configured_state(cfg, state, grid)
      !! The production wiring, nothing test-specific: fill a config with
      !! `rho_0 = RHO0_CFG` and `form = "fv_lite"`, allocate the god state
      !! through `init_from_config` (which is what lands `rho_0` on
      !! `eos%rho0`), then run the one configure stage that owns the PGF
      !! slot's scalars.
      type(config_t), intent(inout) :: cfg
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(inout) :: grid

      cfg%sim_type = "ocean"
      cfg%nx = NX_PHYS
      cfg%ny = NY_PHYS
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%nghost = NGHOST
      cfg%ocean%ic%rho_0 = RHO0_CFG
      cfg%ocean%pgf%form = "fv_lite"

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      call state%init_from_config(cfg, grid)
      call configure_ocean_pgf(cfg, state, compute_rank=0)
   end subroutine build_configured_state

   subroutine build_two_column_field(state, grid)
      !! Flat bed, uniform layer thickness, a west-to-east density ramp.
      !! Every adjacent column pair is an independent two-column hydrostatic
      !! problem; tiling the ramp across the row exercises all of them
      !! instead of a single face.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      integer :: i, j, k, nx, ny
      real(wp) :: ramp

      nx = grid%nx_total
      ny = grid%ny_total
      state%multilayer%h_layer(:, :, :) = H_LAYER
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ramp = real(i - 1, wp)/real(nx - 1, wp)
               state%multilayer%rho_layer(i, j, k) = RHO_BASE &
                                                     + RHO_STRAT*real(NZ - k, wp) &
                                                     + RHO_STEP*ramp
            end do
         end do
      end do
   end subroutine build_two_column_field

   pure subroutine host_two_column_dpdx(rho_layer, h_layer, rho0, dpdx)
      !! The closed form the FV_LITE face pass reduces to at uniform layer
      !! thickness.  `rho0` enters ONLY as the `1/ρ₀` prefactor — the whole
      !! point of the test, so it is an explicit argument and the two arms
      !! call this with the two candidate divisors.
      real(wp), intent(in) :: rho_layer(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: rho0
      real(wp), intent(out) :: dpdx(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp), allocatable :: p_edge(:, :, :)
      real(wp) :: p_c_left, p_c_right

      nx = size(rho_layer, 1)
      ny = size(rho_layer, 2)
      allocate (p_edge(nx, ny, NZ + 1))
      do j = 1, ny
         do i = 1, nx
            p_edge(i, j, NZ + 1) = 0.0_wp
            do k = NZ, 1, -1
               p_edge(i, j, k) = p_edge(i, j, k + 1) &
                                 + GRAVITY*rho_layer(i, j, k)*h_layer(i, j, k)
            end do
         end do
      end do
      dpdx(:, :, :) = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 2, nx
               p_c_right = 0.5_wp*(p_edge(i, j, k) + p_edge(i, j, k + 1))
               p_c_left = 0.5_wp*(p_edge(i - 1, j, k) + p_edge(i - 1, j, k + 1))
               dpdx(i, j, k) = -(1.0_wp/rho0)*(p_c_right - p_c_left)/DX
            end do
         end do
      end do
   end subroutine host_two_column_dpdx

   subroutine run_pgf(state, grid)
      !! One compute pass over the configured god state's own PGF + multilayer
      !! slots, with the face buffers pulled back for the host comparison.
      !!
      !! `mem:separate` contract: EVERY object the kernel touches is mapped
      !! before the call — the multilayer slot AND its PGF scratch companion
      !! AND the metrics — and `dpdx_face` is `create`-mapped scratch, so the
      !! `!$acc update self` is what makes the host read meaningful rather
      !! than stale.  Components are named individually, never the aggregate
      !! (an aggregate D→H copy would overwrite the host descriptors).  All of
      !! it is an inert no-op on host / multicore builds.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t) :: metrics

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(state%multilayer, state%pressure_force)
      call state%multilayer%enter_data()
      call state%pressure_force%enter_data()
      !$acc update device(state%multilayer%h_layer, state%multilayer%rho_layer)
      call ocean_pressure_force_compute(grid, metrics, state%pressure_force, &
                                        state%multilayer)
      !$acc update self(state%pressure_force%dpdx_face%data)
      call state%pressure_force%exit_data()
      call state%multilayer%exit_data()
      !$acc exit data delete(state%multilayer, state%pressure_force)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   ! ---------------------------------------------------------------------
   ! Cases
   ! ---------------------------------------------------------------------

   subroutine test_rho_ref_plumbing(error)
      !! `&ocean_ic_nml rho_0` must reach BOTH PGF reference densities.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid

      checks: block
         call build_configured_state(cfg, state, grid)

         ! Non-vacuity: the configured value must not be the type default,
         ! or every assertion below would hold with the assignments missing.
         call check(error, abs(RHO0_CFG - RHO0_DEFAULT) > 1.0_wp, &
                    "test setup is vacuous: RHO0_CFG equals the 1035 default")
         if (allocated(error)) exit checks

         call check(error, abs(state%eos%rho0 - RHO0_CFG) < 1.0e-12_wp, &
                    "eos%rho0 did not follow &ocean_ic_nml rho_0")
         if (allocated(error)) exit checks
         call check(error, abs(state%pressure_force%rho0 - RHO0_CFG) < 1.0e-12_wp, &
                    "pressure_force%rho0 did not follow the configured rho_0 "// &
                    "(configure_ocean_pgf left the 1035 type default in place)")
         if (allocated(error)) exit checks
         call check(error, abs(state%pressure_force%rho_ref - RHO0_CFG) < 1.0e-12_wp, &
                    "pressure_force%rho_ref did not follow the configured rho_0")
         if (allocated(error)) exit checks
         ! The single-ρ₀-of-record invariant, stated directly.
         call check(error, abs(state%pressure_force%rho0 - state%eos%rho0) < 1.0e-12_wp, &
                    "the EOS and the PGF are running on different reference densities")
      end block checks
      call state%destroy()
   end subroutine test_rho_ref_plumbing

   subroutine test_rho_ref_answer(error)
      !! The ANSWER, not just the scalar: a two-column hydrostatic field whose
      !! face acceleration is exactly `−(1/ρ₀)·Δp_c/dx`, checked against the
      !! configured ρ₀ and shown to be resolvably different from the 1035 one.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      real(wp), allocatable :: want_cfg(:, :, :), want_default(:, :, :)
      integer :: nx, ny
      real(wp) :: err_cfg, err_default, scale

      checks: block
         call build_configured_state(cfg, state, grid)
         call build_two_column_field(state, grid)

         nx = grid%nx_total
         ny = grid%ny_total
         allocate (want_cfg(nx, ny, NZ))
         allocate (want_default(nx, ny, NZ))
         call host_two_column_dpdx(state%multilayer%rho_layer, &
                                   state%multilayer%h_layer, RHO0_CFG, want_cfg)
         call host_two_column_dpdx(state%multilayer%rho_layer, &
                                   state%multilayer%h_layer, RHO0_DEFAULT, want_default)

         call run_pgf(state, grid)

         scale = maxval(abs(want_cfg(2:nx, :, :)))
         ! Non-vacuity 1: the field must actually be driven.
         call check(error, scale > 1.0e-8_wp, &
                    "two-column density contrast produced no pressure gradient")
         if (allocated(error)) exit checks

         err_cfg = maxval(abs(state%pressure_force%dpdx_face%data(2:nx, 1:ny, 1:NZ) &
                              - want_cfg(2:nx, :, :)))
         call check(error, err_cfg < 1.0e-12_wp*scale, &
                    "the two-column PGF answer does not use the configured rho_0")
         if (allocated(error)) exit checks

         ! Non-vacuity 2: the 1035 answer is far outside that tolerance, so
         ! the assertion above genuinely discriminates the two divisors.
         err_default = maxval(abs(want_default(2:nx, :, :) - want_cfg(2:nx, :, :)))
         call check(error, err_default > 1.0e-2_wp*scale, &
                    "test is vacuous: the 1035 answer is within tolerance of "// &
                    "the configured-rho_0 answer")
      end block checks
      call state%destroy()
   end subroutine test_rho_ref_answer

end module test_ocean_pgf_rho_ref
