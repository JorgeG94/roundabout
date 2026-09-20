!! Analytical unit tests for `&ocean_pgf_nml p_top_in_bc` — the
!! top-of-column load `multilayer_state_t%p_top` (Pa) in the FV_MOM6
!! pressure-stack surface boundary condition,
!! `pa(nz+1) = rho_ref*g*eta_geo + p_top`.
!!
!! The knob exists because under an ice-shelf draft the load is `O(5e6 Pa)`
!! while `pa` is an ANOMALY stack about `rho_ref*g*z`: cancel the two
!! inside `pa(nz+1)` and the whole stack stays `O(1e4 Pa)`, which shrinks
!! every downstream `h_neglect` divisor leak by the same 500x.  What makes
!! that legitimate — and not a double count of the barotropic
!! `eta_forcing` seam — is the theorem the tests below pin down.
!!
!! **Theorem.** Perturb the top BC only, `pa(.,nz+1) -> pa(.,nz+1) + dp`
!! with `dp(i,j)` independent of `k`.  Every `dpa`, `intz_dpa` and
!! `intx_dpa` is unchanged and `intx_pa(K) -> intx_pa(K) + 0.5*(dp_L+dp_R)`
!! for EVERY `K`, so the Pass-3 numerator moves by
!! `0.5*(h_L+h_R)*(dp_L-dp_R)` and
!!
!!    dPFu(k) = -(1/rho_0)*d(dp)/dx * (h_L+h_R)/(h_L+h_R+h_neglect)
!!
!! — the SAME acceleration in every layer, up to the `h_neglect` divisor.
!! The split solver replaces the depth mean of the layer PGF with the
!! barotropic solution, so a depth-uniform `dPFu` is annihilated in the
!! baroclinic part and the seam keeps sole ownership of the barotropic
!! response.
!!
!! Cases:
!!   1. `uniform_p_top_bit_identical` — a constant has no gradient: a
!!      uniform 5e6 Pa load over a stratified column set leaves the face
!!      PGF BIT-identical to `p_top = 0`.
!!   2. `theorem_depth_uniform_pfu` — the theorem itself, on a stratified,
!!      non-uniform-thickness stack with flat isopycnals and flat eta:
!!      `PFu(k)` equals `-(1/rho_0)*grad p_top` at EVERY layer, and the
!!      layer PGF minus its thickness-weighted depth mean is zero to
!!      round-off.  Also asserts the `h_neglect` leak bound.
!!   3. `theorem_depth_uniform_pfu_reconstruct` — case 2 through the
!!      `reconstruct_for_pressure` branch, which seeds the SAME `pa(nz+1)`
!!      and reuses the same face assembly.
!!   4. `isostatic_rest_sloping_load` — the physics gate.  Flat bed,
!!      uniform density, a SLOPING draft with `sum(h) = b - z_draft` and
!!      `p_top = (rho_ref*GRAVITY)*z_draft` built with the same product as
!!      the `pa` seed: the column is in discrete hydrostatic balance and
!!      every face force is zero to the rounding of `p_top`.  With the `+ p_top` line removed the
!!      residual is `g*grad(z_draft)` ~ 0.1 m/s^2 — the step-1 blow-up the
!!      knob exists to prevent.
!!   5. `knob_off_ignores_p_top` — a ramped, non-zero `p_top` with the knob
!!      OFF is bit-identical to `p_top = 0` (the default-off contract).
!!   6. `pa_sign_and_monotone` — `pa(nz+1) = p_top >= 0` at `eta = 0`, and
!!      `pa` increases DOWNWARD (toward `k = 1`, the bed) for water denser
!!      than `rho_ref`.
!!   7. `validate_config_p_top_in_bc` — the envelope: fv_mom6 accepted
!!      (with and without the reconstruction), every other form refused.
module test_ocean_pgf_p_top_bc
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   use rdb_config, only: config_t, read_config_from_string, validate_config, &
                         p_top_has_producer
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

   public :: collect_ocean_pgf_p_top_bc_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 1000.0_wp

   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: RHO_REF = 1035.0_wp

   real(wp), parameter :: H_NEGLECT = 1.0e-10_wp
      !! Mirrors `ocean_pressure_force_t%h_neglect` — the face-divisor
      !! regulariser whose leak case 2 bounds.

   ! Exactly-representable layer thicknesses (k = 1 bed .. k = NZ surface),
   ! so `e_face(nz+1) = -b + sum(h)` is exactly 0 when `b = sum(h)` and the
   ! `p_top`/`pa` cancellation in case 4 is bit-exact rather than merely
   ! small.
   real(wp), parameter :: H_STACK(NZ) = [256.0_wp, 128.0_wp, 512.0_wp, 64.0_wp]
   real(wp), parameter :: H_TOTAL = 960.0_wp

   real(wp), parameter :: P_TOP_UNIFORM = 5.0e6_wp
      !! ~500 m of Boussinesq-isostatic ice draft.
   real(wp), parameter :: DPDX_TOP = 1.0_wp
      !! Pa per metre — the ramped load's gradient in cases 2/3/5.

contains

   subroutine collect_ocean_pgf_p_top_bc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("uniform_p_top_bit_identical", test_uniform_bitident), &
                  new_unittest("theorem_depth_uniform_pfu", test_theorem_pcm), &
                  new_unittest("theorem_depth_uniform_pfu_reconstruct", test_theorem_recon), &
                  new_unittest("isostatic_rest_sloping_load", test_isostatic_rest), &
                  new_unittest("knob_off_ignores_p_top", test_knob_off), &
                  new_unittest("pa_sign_and_monotone", test_pa_sign), &
                  new_unittest("validate_config_p_top_in_bc", test_validate_config), &
                  new_unittest("p_top_producer_is_psurf_or_cavity", test_p_top_producer) &
                  ]
   end subroutine collect_ocean_pgf_p_top_bc_tests

   ! ------------------------------------------------------------------
   ! Harness
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, DX, DX)
   end subroutine make_grid

   subroutine run_pgf(grid, ms, pgf, b, eos)
      !! One FV_MOM6 compute pass against a caller-supplied bathymetry.
      !!
      !! `mem:separate` discipline (CLAUDE.md): `ms`, its `p_top`
      !! companion and `pgf` are ALL mapped before the kernel; `p_top` is
      !! a HOST-set input so it is pushed with an explicit
      !! `!$acc update device` after `enter_data` (it is `copyin`-mapped,
      !! so this is belt-and-braces, but it is the rule and the directive
      !! is an inert no-op on the host build); the read-back uses
      !! `update self` on the COMPONENT arrays, never the aggregate
      !! derived type.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: b(:, :)
      type(eos_t), intent(in), optional :: eos
      type(ocean_metrics_t) :: metrics

      call make_cartesian_metrics(metrics, grid)
      call pgf%set_bathymetry(b)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      !$acc update device(ms%p_top, ms%h_layer, ms%rho_layer)
      if (present(eos)) then
         call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      else
         call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      end if
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data, pgf%pa%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   pure subroutine seed_stratified(ms)
      !! Horizontally uniform stack: flat isopycnals, `H_STACK` thicknesses,
      !! density decreasing upward.  With `b = H_TOTAL` this is a resting
      !! stratified column set, so the unloaded PGF is zero and any signal
      !! the tests see comes from `p_top` alone.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: k
      do k = 1, NZ
         ms%h_layer(:, :, k) = H_STACK(k)
         ms%rho_layer(:, :, k) = RHO_REF + 1.2_wp - 0.3_wp*real(k - 1, wp)
      end do
   end subroutine seed_stratified

   pure function x_of(i) result(x)
      !! Cell-centre x (m) from the TOTAL index, ghosts included.
      integer, intent(in) :: i
      real(wp) :: x
      x = (real(i, wp) - 0.5_wp)*DX
   end function x_of

   pure subroutine seed_ramped_p_top(ms)
      !! `p_top = P_TOP_UNIFORM + DPDX_TOP*x`, filled over the WHOLE array
      !! (ghosts included) so every u-face the kernel touches sees the
      !! same linear gradient.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: i, j
      do j = 1, size(ms%p_top, 2)
         do i = 1, size(ms%p_top, 1)
            ms%p_top(i, j) = P_TOP_UNIFORM + DPDX_TOP*x_of(i)
         end do
      end do
   end subroutine seed_ramped_p_top

   subroutine make_pgf(grid, pgf, p_top_in_bc, reconstruct)
      type(hgrid_t), intent(in) :: grid
      type(ocean_pressure_force_t), intent(out) :: pgf
      logical, intent(in) :: p_top_in_bc
      logical, intent(in), optional :: reconstruct
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO0
      pgf%rho_ref = RHO_REF
      pgf%p_top_in_bc = p_top_in_bc
      if (present(reconstruct)) pgf%reconstruct_for_pressure = reconstruct
   end subroutine make_pgf

   pure subroutine make_linear_eos(eos)
      !! Pressure-independent linear EOS, so the reconstruction branch's
      !! own in-layer pressure cannot contaminate the BC test.
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%T_ref = 10.0_wp
      eos%S_ref = 35.0_wp
      eos%alpha_T = 0.2_wp
      eos%beta_S = 0.78_wp
      eos%is_init = .true.
   end subroutine make_linear_eos

   ! ------------------------------------------------------------------
   ! Case 1 — a constant has no gradient
   ! ------------------------------------------------------------------

   subroutine test_uniform_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_on, pgf_off
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), dpdx_off(:, :, :), dpdy_off(:, :, :)

      checks: block
         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call seed_stratified(ms)
         allocate (b(grid%nx_total, grid%ny_total), source=H_TOTAL)

         call make_pgf(grid, pgf_off, .false.)
         call run_pgf(grid, ms, pgf_off, b)
         allocate (dpdx_off, source=pgf_off%dpdx_face%data)
         allocate (dpdy_off, source=pgf_off%dpdy_face%data)

         ms%p_top = P_TOP_UNIFORM
         call make_pgf(grid, pgf_on, .true.)
         call run_pgf(grid, ms, pgf_on, b)

         call check(error, all(pgf_on%dpdx_face%data == dpdx_off), &
                    "uniform p_top must leave PFu bit-identical (no gradient)")
         if (allocated(error)) exit checks
         call check(error, all(pgf_on%dpdy_face%data == dpdy_off), &
                    "uniform p_top must leave PFv bit-identical (no gradient)")
      end block checks
      call pgf_on%destroy(); call pgf_off%destroy(); call ms%destroy()
   end subroutine test_uniform_bitident

   ! ------------------------------------------------------------------
   ! Cases 2/3 — the theorem
   ! ------------------------------------------------------------------

   subroutine test_theorem_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      call theorem_body(error, use_reconstruct=.false.)
   end subroutine test_theorem_pcm

   subroutine test_theorem_recon(error)
      type(error_type), allocatable, intent(out) :: error
      call theorem_body(error, use_reconstruct=.true.)
   end subroutine test_theorem_recon

   subroutine theorem_body(error, use_reconstruct)
      !! Flat isopycnals, flat eta, a linear-in-x `p_top`.  Asserts, at
      !! every interior u-face and every layer:
      !!
      !!   (a) `PFu(k)` equals the analytic `-(1/rho_0)*dp_top/dx`;
      !!   (b) the LOAD-INDUCED part `PFu_on - PFu_off` equals it to the
      !!       `h_neglect` bound `h_neglect/(h_L+h_R)` — a *formula* check
      !!       on the leak, not just a tolerance;
      !!   (c) that part is DEPTH-UNIFORM: its deviation from its own
      !!       thickness-weighted depth mean is zero to round-off, which
      !!       is exactly what the split solver's depth-mean replacement
      !!       annihilates.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_reconstruct
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_on, pgf_off
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), dpdx_off(:, :, :)
      real(wp) :: analytic, dpfu(NZ), h_face(NZ), wmean, sum_h
      real(wp) :: max_rel_direct, max_rel_leak, max_rel_uniform, leak_bound
      integer :: i, j, k, nx, ny

      checks: block
         call make_grid(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call seed_stratified(ms)
         call make_linear_eos(eos)
         call seed_TS_uniform(ms)
         allocate (b(nx, ny), source=H_TOTAL)

         ! Reference run: same state, knob OFF, p_top still zero.
         call make_pgf(grid, pgf_off, .false., reconstruct=use_reconstruct)
         if (use_reconstruct) then
            call run_pgf(grid, ms, pgf_off, b, eos=eos)
         else
            call run_pgf(grid, ms, pgf_off, b)
         end if
         allocate (dpdx_off, source=pgf_off%dpdx_face%data)

         call seed_ramped_p_top(ms)
         call make_pgf(grid, pgf_on, .true., reconstruct=use_reconstruct)
         if (use_reconstruct) then
            call run_pgf(grid, ms, pgf_on, b, eos=eos)
         else
            call run_pgf(grid, ms, pgf_on, b)
         end if

         analytic = -DPDX_TOP/RHO0
         ! Flat horizontal stack => h_L = h_R = H_STACK(k) at every face.
         do k = 1, NZ
            h_face(k) = 2.0_wp*H_STACK(k)
         end do
         leak_bound = H_NEGLECT/(2.0_wp*minval(H_STACK))

         max_rel_direct = 0.0_wp
         max_rel_leak = 0.0_wp
         max_rel_uniform = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               do k = 1, NZ
                  dpfu(k) = pgf_on%dpdx_face%data(i, j, k) - dpdx_off(i, j, k)
                  max_rel_direct = max(max_rel_direct, &
                                       abs(pgf_on%dpdx_face%data(i, j, k) - analytic))
                  max_rel_leak = max(max_rel_leak, abs(dpfu(k) - analytic))
               end do
               ! Thickness-weighted depth mean of the load-induced part.
               wmean = 0.0_wp
               sum_h = 0.0_wp
               do k = 1, NZ
                  wmean = wmean + h_face(k)*dpfu(k)
                  sum_h = sum_h + h_face(k)
               end do
               wmean = wmean/sum_h
               do k = 1, NZ
                  max_rel_uniform = max(max_rel_uniform, abs(dpfu(k) - wmean))
               end do
            end do
         end do
         max_rel_direct = max_rel_direct/abs(analytic)
         max_rel_leak = max_rel_leak/abs(analytic)
         max_rel_uniform = max_rel_uniform/abs(analytic)

         call check(error, max_rel_direct < 1.0e-9_wp, &
                    "PFu under a ramped p_top must equal -(1/rho_0)*grad p_top "// &
                    "at EVERY layer")
         if (allocated(error)) exit checks
         ! A formula check, not a tolerance: the only departure from the
         ! analytic is the `h_neglect` divisor, worst at the THINNEST
         ! layer.  The 5 % head-room covers nothing but rounding of the
         ! bound itself (the leak is 7.8e-13, roundoff 1e-16).
         call check(error, max_rel_leak <= 1.05_wp*leak_bound, &
                    "the load-induced PFu must match the analytic to within the "// &
                    "h_neglect face-divisor bound h_neglect/(h_L+h_R)")
         if (allocated(error)) exit checks
         call check(error, max_rel_uniform < 1.0e-11_wp, &
                    "the load-induced PFu must be DEPTH-UNIFORM (zero deviation "// &
                    "from its thickness-weighted depth mean) — the piece the "// &
                    "split solver annihilates")
      end block checks
      call pgf_on%destroy(); call pgf_off%destroy(); call ms%destroy()
   end subroutine theorem_body

   pure subroutine seed_TS_uniform(ms)
      !! Uniform T/S consistent with `seed_stratified`'s `rho_layer` only
      !! in the sense that neither varies horizontally — the reconstruct
      !! branch builds its own density from T/S, and the theorem is about
      !! the `pa(nz+1)` seed, which is common to both branches.
      type(multilayer_state_t), intent(inout) :: ms
      integer :: k
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
            (10.0_wp + 0.5_wp*real(k - 1, wp))*ms%h_layer(:, :, k)
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*ms%h_layer(:, :, k)
      end do
   end subroutine seed_TS_uniform

   ! ------------------------------------------------------------------
   ! Case 4 — isostatic rest under a sloping load
   ! ------------------------------------------------------------------

   subroutine test_isostatic_rest(error)
      !! Flat bed `b`, uniform density `rho = rho_ref`, a draft sloping in
      !! x, the water column thinned to `b - z_draft` so
      !! `eta_geo = e_face(nz+1) = -z_draft` exactly, and the load built
      !! with the SAME product as the `pa` seed,
      !! `p_top = (rho_ref*GRAVITY)*z_draft`.  Then
      !!
      !!   pa(nz+1) = (rho_ref*GRAVITY)*(-z_draft) + (rho_ref*GRAVITY)*z_draft
      !!            = 0   when the product is rounded before the add; under
      !!                FMA contraction, the rounding error of `p_top`
      !!
      !! and `dpa = (rho - rho_ref)*g*h = 0`, so the WHOLE `pa` stack is
      !! zero to that rounding and every face force is at round-off.  Without the
      !! `+ p_top` term `pa(nz+1) = -rho_ref*g*z_draft` varies across the
      !! face and the PGF is `g*grad(z_draft)` ~ 1e-1 m/s^2.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: BED = 1024.0_wp
      real(wp), parameter :: DRAFT0 = 256.0_wp, DRAFT_STEP = 0.5_wp
      real(wp), parameter :: PGF_REST_TOL = 1.0e-13_wp
         !! m/s^2.  ~1e12 below the `g*grad(z_draft)` ~ 1e-1 m/s^2 the
         !! stack carries without the load term.
      real(wp) :: z_draft, water, pgf_max, pa_tol
      integer :: i, j, k, nx, ny

      checks: block
         call make_grid(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         allocate (b(nx, ny), source=BED)

         do j = 1, ny
            do i = 1, nx
               ! Draft, bed and quarter-column thickness all land on exact
               ! binary fractions, so `e_face(nz+1) = -z_draft` EXACTLY and
               ! the `pa(nz+1)` cancellation is bit-exact, not merely small.
               z_draft = DRAFT0 + DRAFT_STEP*real(i - 1, wp)
               water = BED - z_draft
               do k = 1, NZ
                  ms%h_layer(i, j, k) = water*0.25_wp
                  ms%rho_layer(i, j, k) = RHO_REF
               end do
               ms%p_top(i, j) = (RHO_REF*GRAVITY)*z_draft
            end do
         end do

         call make_pgf(grid, pgf, .true.)
         call run_pgf(grid, ms, pgf, b)

         pgf_max = max(maxval(abs(pgf%dpdx_face%data)), &
                       maxval(abs(pgf%dpdy_face%data)))
         ! NOT `== 0`: the cancellation `(rho_ref*g)*(-z_draft) + p_top` is
         ! bit-exact only when the product is ROUNDED before the add.  A
         ! toolchain that contracts the kernel's `rho_ref*g*eta + p_top` into
         ! a fused multiply-add (nvfortran on the device; any `-mfma` host
         ! build) keeps the product exact, so what survives is the rounding
         ! error of the HOST-built `p_top` itself: |pa| <= ulp(p_top)/2.
         ! That is the honest bound -- ~2e-10 Pa on a 2.6e6 Pa load, a face
         ! force ~1e-16 m/s^2 against the 1e-1 m/s^2 of the un-loaded stack
         ! -- so assert THAT, on every toolchain, rather than a bit-zero only
         ! some code generators can deliver.  (Verified: gfortran gives
         ! exactly 0; the V100 build gives a non-zero value inside the bound.)
         pa_tol = spacing(maxval(ms%p_top))
         call check(error, maxval(abs(pgf%pa%data)) <= pa_tol, &
                    "the whole pa stack must vanish to the rounding of p_top "// &
                    "at isostatic rest")
         if (allocated(error)) exit checks
         call check(error, pgf_max <= PGF_REST_TOL, &
                    "isostatic rest under a sloping load must give a face "// &
                    "force at round-off (the pa stack cancels)")
      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_isostatic_rest

   ! ------------------------------------------------------------------
   ! Case 5 — the default-off contract
   ! ------------------------------------------------------------------

   subroutine test_knob_off(error)
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_zero, pgf_loaded
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), dpdx_zero(:, :, :), dpdy_zero(:, :, :)

      checks: block
         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call seed_stratified(ms)
         allocate (b(grid%nx_total, grid%ny_total), source=H_TOTAL)

         call make_pgf(grid, pgf_zero, .false.)
         call run_pgf(grid, ms, pgf_zero, b)
         allocate (dpdx_zero, source=pgf_zero%dpdx_face%data)
         allocate (dpdy_zero, source=pgf_zero%dpdy_face%data)

         ! A large, strongly sloping load — with the knob OFF none of it
         ! may reach the pressure stack.
         call seed_ramped_p_top(ms)
         call make_pgf(grid, pgf_loaded, .false.)
         call run_pgf(grid, ms, pgf_loaded, b)

         call check(error, all(pgf_loaded%dpdx_face%data == dpdx_zero), &
                    "p_top_in_bc=.false. with a non-zero p_top must leave PFu "// &
                    "bit-identical")
         if (allocated(error)) exit checks
         call check(error, all(pgf_loaded%dpdy_face%data == dpdy_zero), &
                    "p_top_in_bc=.false. with a non-zero p_top must leave PFv "// &
                    "bit-identical")
      end block checks
      call pgf_zero%destroy(); call pgf_loaded%destroy(); call ms%destroy()
   end subroutine test_knob_off

   ! ------------------------------------------------------------------
   ! Case 6 — sign + monotonicity of the loaded pressure stack
   ! ------------------------------------------------------------------

   subroutine test_pa_sign(error)
      !! MOM6 shipped a NEGATIVE EOS pressure in one path for years, so
      !! the sign convention gets its own assertion here too: with
      !! `eta = 0` and `p_top >= 0` the stack must start at `p_top` and
      !! grow DOWNWARD (toward `k = 1`, the bed), never upward.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp) :: max_seed_err, min_step
      integer :: i, j, k, nx, ny

      checks: block
         call make_grid(grid, 10, 8)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call seed_stratified(ms)     ! rho_layer > RHO_REF in every layer
         call seed_ramped_p_top(ms)
         allocate (b(nx, ny), source=H_TOTAL)   ! => eta_geo = 0 exactly

         call make_pgf(grid, pgf, .true.)
         call run_pgf(grid, ms, pgf, b)

         max_seed_err = 0.0_wp
         min_step = huge(1.0_wp)
         do j = 1, ny
            do i = 1, nx
               max_seed_err = max(max_seed_err, &
                                  abs(pgf%pa%data(i, j, NZ + 1) - ms%p_top(i, j)))
               do k = NZ, 1, -1
                  min_step = min(min_step, pgf%pa%data(i, j, k) - pgf%pa%data(i, j, k + 1))
               end do
            end do
         end do

         call check(error, max_seed_err == 0.0_wp, &
                    "at eta = 0 the loaded surface BC must be exactly pa(nz+1) = p_top")
         if (allocated(error)) exit checks
         call check(error, minval(pgf%pa%data(:, :, NZ + 1)) > 0.0_wp, &
                    "pa(nz+1) must be >= 0 for p_top >= 0 at eta = 0")
         if (allocated(error)) exit checks
         call check(error, min_step > 0.0_wp, &
                    "pa must INCREASE downward (toward k = 1) for water denser "// &
                    "than rho_ref")
      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_pa_sign

   ! ------------------------------------------------------------------
   ! Case 7 — the envelope
   ! ------------------------------------------------------------------

   subroutine test_validate_config(error)
      !! Only the FV_MOM6 family builds a `pa` stack: `mont` hard-zeroes
      !! `M(nz)` and `fv_lite`/`fv_wright` seed `p_edge(nz+1) = 0`, so
      !! there is no boundary condition anywhere else to inject the load
      !! into.  Refuse rather than run a silently inert knob on a config
      !! the user believes is carrying an ice load.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true. /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "p_top_in_bc with form='fv_mom6' must be accepted")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true., '// &
                         'reconstruct_for_pressure = .true. /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "p_top_in_bc with the FV_MOM6 reconstruction branch must be "// &
                    "accepted (same pa(nz+1) seed)")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "mont", p_top_in_bc = .true. /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "p_top_in_bc with form='mont' must be refused (no pa stack)")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_wright", p_top_in_bc = .true. /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "p_top_in_bc with form='fv_wright' must be refused "// &
                    "(p_edge(nz+1) = 0, not a pa stack)")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_lite", p_top_in_bc = .true. /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "p_top_in_bc with form='fv_lite' must be refused")
         if (allocated(error)) exit checks

         ! Default off must still validate for every form.
         call parse_case(cfg, '&ocean_pgf_nml form = "mont" /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "the default p_top_in_bc=.false. must validate for any form")
      end block checks
   end subroutine test_validate_config

   subroutine test_p_top_producer(error)
      !! WHO WRITES `ms%p_top` — the predicate behind the `p_top_in_bc`
      !! inert-knob warning.
      !!
      !! `p_top = metrics%p_ice_ref + sf%p_surf` has TWO producers.  The
      !! warning used to name only the `&ocean_psurf_nml` seam, so a
      !! configuration with an ice-shelf cavity — where `p_top` carries
      !! the full `rho_ref*g*z_draft` ice load, and where `p_top_in_bc`
      !! is not merely live but REQUIRED for a varying draft — was told
      !! its load was being ignored.  That is worse than a missing
      !! warning: it invites the user to turn OFF the knob that is
      !! carrying the physics.
      !!
      !! The predicate is tested rather than the log line because the
      !! predicate is the thing with a contract; the sentence is its
      !! rendering.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      checks: block
         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true. /')
         call check(error,.not. p_top_has_producer(cfg), &
                    "with neither psurf nor a cavity, p_top really is the zero "// &
                    "array and the knob really is inert")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true. /'// &
                         new_line("a")//"&ocean_psurf_nml enable = .true. /")
         call check(error, p_top_has_producer(cfg), &
                    "the surface-pressure seam is a producer (it fills sf%p_surf)")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true. /'// &
                         new_line("a")//"&ocean_cavity_dyn_nml enable = .true., "// &
                         "draft_config = 'flat', draft_depth = 300.0 /")
         call check(error, p_top_has_producer(cfg), &
                    "an ice-shelf cavity is a producer too — p_ice_ref = "// &
                    "rho_ref*g*z_draft is assembled into p_top by "// &
                    "configure_ocean_cavity, and p_top_in_bc is REQUIRED there")
         if (allocated(error)) exit checks

         call parse_case(cfg, '&ocean_pgf_nml form = "fv_mom6", p_top_in_bc = .true. /'// &
                         new_line("a")//"&ocean_psurf_nml enable = .true. /"// &
                         new_line("a")//"&ocean_cavity_dyn_nml enable = .true., "// &
                         "draft_config = 'flat', draft_depth = 300.0 /")
         call check(error, p_top_has_producer(cfg), &
                    "both at once is still a producer (p_top is their SUM)")
      end block checks
   end subroutine test_p_top_producer

   subroutine parse_case(cfg, extra)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: extra
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_bt_nml n_inner = 8 /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
   end subroutine parse_case

end module test_ocean_pgf_p_top_bc
