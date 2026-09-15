!! Unit tests for the ocean pressure-gradient kernel
!! (`rdb_ocean_pressure_force`): the rest-state / analytic physics cases on
!! ALIGNED columns, where every variant must agree, plus the per-variant
!! scratch ALLOCATION GATE.
!!
!! The physics cases run on the slot default (`mont`).  On aligned columns
!! that form is algebraically identical to `fv_lite`, so these assertions
!! bind every variant, not just one — the variant-specific behaviour lives in
!! `test_ocean_pgf_mont`, `test_ocean_pgf_fv`, `test_ocean_pgf_fv_mom6` and
!! `test_ocean_pgf_grounded`.
!!
!! Cases:
!!   * Uniform density, uniform thickness, zero velocity — the
!!     pressure gradient must vanish identically.  Constancy
!!     preservation for the unstratified rest state.
!!   * Stratified rest state — different density per layer but
!!     uniform across columns; horizontal pressure gradient still
!!     zero in the domain interior.  Catches per-layer integration
!!     bugs that uniform-density tests miss.
!!   * Horizontal density gradient in one layer — analytic check
!!     against the closed-form hydrostatic pressure formula.  Asserts
!!     (a) acceleration in the gradient layer matches analytic
!!     to round-off, (b) layers above stay quiescent (no spurious
!!     baroclinic propagation upward through the column).
!!   * The `scratch_gated` allocation map, per variant.
module test_ocean_pgf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       ocean_pressure_force_apply, &
                                       OPGF_VARIANT_MONT, &
                                       OPGF_VARIANT_FV_LITE, &
                                       OPGF_VARIANT_FV_WRIGHT, &
                                       OPGF_VARIANT_GPRIME, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_pgf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_rest_state", test_uniform_rest_state), &
                  new_unittest("stratified_rest_state", test_stratified_rest_state), &
                  new_unittest("horizontal_density_gradient", test_horizontal_gradient), &
                  new_unittest("uniform_rho_varying_bathy_cancels", test_uniform_rho_varying_bathy), &
                  new_unittest("scratch_gate_off_allocates_everything", test_gate_off_allocates_all), &
                  new_unittest("scratch_gate_per_variant", test_gate_per_variant), &
                  new_unittest("scratch_gate_recon", test_gate_recon), &
                  new_unittest("scratch_gate_bytes_track_allocation", test_gate_bytes) &
                  ]
   end subroutine collect_ocean_pgf_tests

   ! =====================================================================
   ! Allocation-gate tests (`ocean_pressure_force_t%scratch_gated`).
   !
   ! These assert the SHAPE of the gate, not the physics: which buffers
   ! exist for which `variant`.  A gate that is too tight leaves a buffer
   ! unallocated on a path that reaches it — and under
   ! `-gpu=...,mem:separate` that is a silent wrong answer, not a crash —
   ! so the per-variant map is pinned here rather than left to review.
   ! =====================================================================

   subroutine gate_flags(pgf, has_p_edge, has_z_centre, has_rho_insitu, &
                         has_fv_mom6, has_recon, has_out, has_mont_M)
      !! Which buffer groups the slot actually allocated.
      type(ocean_pressure_force_t), intent(in) :: pgf
      logical, intent(out) :: has_p_edge, has_z_centre, has_rho_insitu
      logical, intent(out) :: has_fv_mom6, has_recon, has_out
      logical, intent(out), optional :: has_mont_M
      has_p_edge = allocated(pgf%p_edge%data)
      has_z_centre = allocated(pgf%z_centre%data)
      has_rho_insitu = allocated(pgf%rho_insitu%data)
      if (present(has_mont_M)) has_mont_M = allocated(pgf%mont_M%data)
      ! The FV_MOM6 stack is all-or-nothing; require every member so a
      ! partial gate cannot pass.
      has_fv_mom6 = allocated(pgf%e_face%data) .and. allocated(pgf%pa%data) &
                    .and. allocated(pgf%intz_dpa%data) &
                    .and. allocated(pgf%intx_pa%data) .and. allocated(pgf%inty_pa%data) &
                    .and. allocated(pgf%intx_dpa%data) .and. allocated(pgf%inty_dpa%data)
      has_recon = allocated(pgf%recon_T_t%data) .and. allocated(pgf%recon_T_b%data) &
                  .and. allocated(pgf%recon_S_t%data) .and. allocated(pgf%recon_S_b%data)
      has_out = allocated(pgf%dpdx_face%data) .and. allocated(pgf%dpdy_face%data)
   end subroutine gate_flags

   subroutine test_gate_off_allocates_all(error)
      !! `scratch_gated = .false.` (the type default, and what every direct
      !! `pgf%init(...)` call site relies on) must allocate all 16 buffers
      !! whatever `variant` says — those call sites set `variant` AFTER
      !! `init`, so the gate would decide on a stale value.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_pressure_force_t) :: pgf
      logical :: hp, hz, hr, hm, hc, ho, hM_

      call make_grid(grid, 6, 5, 1.0e3_wp, 1.0e3_wp)
      call pgf%init(grid, nz_ml=NZ)          ! scratch_gated left .false.
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hp .and. hz .and. hr .and. hm .and. hc .and. ho .and. hM_, &
                 "ungated init must allocate every PGF scratch buffer")
      call pgf%destroy()
   end subroutine test_gate_off_allocates_all

   subroutine test_gate_per_variant(error)
      !! With the gate on, each variant gets exactly the buffers its
      !! branch of `ocean_pressure_force_compute` can reach.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_pressure_force_t) :: pgf
      logical :: hp, hz, hr, hm, hc, ho, hM_

      call make_grid(grid, 6, 5, 1.0e3_wp, 1.0e3_wp)

      ! MONT: interface heights + the Montgomery potential.  It builds `M`
      ! straight from `rho_layer` and never forms a pressure stack, so
      ! `p_edge` must NOT appear — the buffer the old no-geopotential
      ! arithmetic used is precisely the one the correct form does not need.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_MONT
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hz .and. hM_ .and. ho .and. .not. (hp .or. hr .or. hm .or. hc), &
                 "MONT must gate to z_centre + mont_M + the two output face "// &
                 "buffers, and must NOT allocate p_edge")
      call pgf%destroy()
      if (allocated(error)) return

      ! FV_LITE: adds z_centre; rho comes from ms%rho_layer, not rho_insitu.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_FV_LITE
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hp .and. hz .and. ho .and. .not. (hr .or. hm .or. hc .or. hM_), &
                 "FV_LITE must gate to p_edge + z_centre (no rho_insitu, no "// &
                 "FV_MOM6, no mont_M)")
      call pgf%destroy()
      if (allocated(error)) return

      ! FV_WRIGHT: adds the Picard in-situ density.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_FV_WRIGHT
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hp .and. hz .and. hr .and. ho .and. .not. (hm .or. hc .or. hM_), &
                 "FV_WRIGHT must gate to p_edge + z_centre + rho_insitu")
      call pgf%destroy()
      if (allocated(error)) return

      ! GPRIME returns before Pass 1 — only the outputs are ever written.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_GPRIME
      call pgf%init(grid, nz_ml=2)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, ho .and. .not. (hp .or. hz .or. hr .or. hm .or. hc .or. hM_), &
                 "GPRIME must gate to the output face buffers only")
      call pgf%destroy()
      if (allocated(error)) return

      ! FV_MOM6 takes the layer-integrated stack and NOT the FV p_edge.
      ! Its face assembly never reads `z_centre`, so with the grounded-layer
      ! gate off that buffer must stay unallocated — the whole point of the
      ! gate is not to carry a 3D array nobody reads.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_FV_MOM6
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hm .and. ho .and. .not. (hp .or. hz .or. hr .or. hc .or. hM_), &
                 "FV_MOM6 must gate to the e_face/pa/int* stack, not p_edge")
      call pgf%destroy()
      if (allocated(error)) return

      ! FV_MOM6 + grounded-layer gate: `z_centre` is the gate's ONLY input on
      ! this path, so it must appear exactly when `skip_nonoverlap` is on.
      ! That flag is config-time (vcoord + namelist), which is why the ocean
      ! driver latches it BEFORE `init` alongside `variant` — if it arrived
      ! only at configure time, Pass 4 would read an unallocated buffer, and
      ! under `-gpu=...,mem:separate` that is a silent wrong answer.
      ! Kept last: `destroy` does not restore component defaults, so the flag
      ! set here would leak into any block that followed.
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%skip_nonoverlap = .true.
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho, hM_)
      call check(error, hm .and. hz .and. ho .and. .not. (hp .or. hr .or. hc .or. hM_), &
                 "FV_MOM6 with skip_nonoverlap must additionally allocate z_centre")
      call pgf%destroy()
   end subroutine test_gate_per_variant

   subroutine test_gate_recon(error)
      !! The in-layer reconstruction scratch keys off
      !! `reconstruct_for_pressure`, not off `variant` (configure fails
      !! loud if the knob is set with any variant other than FV_MOM6).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_pressure_force_t) :: pgf
      logical :: hp, hz, hr, hm, hc, ho

      call make_grid(grid, 6, 5, 1.0e3_wp, 1.0e3_wp)
      pgf%scratch_gated = .true.
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%reconstruct_for_pressure = .true.
      call pgf%init(grid, nz_ml=NZ)
      call gate_flags(pgf, hp, hz, hr, hm, hc, ho)
      call check(error, hm .and. hc .and. ho, &
                 "reconstruct_for_pressure must add the recon_{T,S}_{t,b} scratch")
      call pgf%destroy()
   end subroutine test_gate_recon

   subroutine test_gate_bytes(error)
      !! `bytes()` must track the gate: a gated-off buffer contributes 0,
      !! so the counted footprint of a gated slot is strictly smaller than
      !! the ungated one and equals the sum of what it actually allocated.
      !! (If it did not, the startup budget would over-report by exactly
      !! the memory the gate just saved.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_pressure_force_t) :: pgf_all, pgf_gated
      integer(kind=8) :: b_all, b_gated, b_expect

      call make_grid(grid, 6, 5, 1.0e3_wp, 1.0e3_wp)

      call pgf_all%init(grid, nz_ml=NZ)
      b_all = pgf_all%bytes()

      pgf_gated%scratch_gated = .true.
      pgf_gated%variant = OPGF_VARIANT_GPRIME
      call pgf_gated%init(grid, nz_ml=NZ)
      b_gated = pgf_gated%bytes()

      ! GPRIME keeps only b + the two face outputs.
      b_expect = int(size(pgf_gated%b), 8)*8_8 &
                 + int(size(pgf_gated%dpdx_face%data), 8)*8_8 &
                 + int(size(pgf_gated%dpdy_face%data), 8)*8_8
      call check(error, b_gated == b_expect, &
                 "gated bytes() must equal the sum of the buffers it allocated")
      if (.not. allocated(error)) &
         call check(error, b_gated < b_all, &
                    "gated bytes() must be strictly smaller than ungated")
      call pgf_all%destroy()
      call pgf_gated%destroy()
   end subroutine test_gate_bytes

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(grid, ms, pgf, metrics)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_metrics_t), intent(inout) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(pgf)
      call pgf%enter_data()
   end subroutine map_in

   subroutine map_out(ms, pgf, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_metrics_t), intent(inout) :: metrics
      call pgf%exit_data()
      !$acc exit data delete(pgf)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_uniform_rest_state(error)
      !! Uniform rho_layer (= rho0), uniform h_layer, zero velocity
      !! initial state.  After one step of PGF compute+apply, u and
      !! v must remain at round-off.  The simplest hydrostatic
      !! balance check.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: max_u, max_v
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         ms%h_layer = H0
         ms%rho_layer = pgf%rho0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(grid, ms, pgf, metrics)
         call ocean_pressure_force_compute(grid, metrics, pgf, ms)
         call ocean_pressure_force_apply(pgf, ms, DT)
         call map_out(ms, pgf, metrics)

         max_u = maxval(abs(ms%u_face_x_layer))
         max_v = maxval(abs(ms%v_face_y_layer))
         call check(error, max_u < 1.0e-12_wp, &
                    "uniform rest state: PGF injected u")
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0e-12_wp, &
                    "uniform rest state: PGF injected v")

      end block checks
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_uniform_rest_state

   subroutine test_stratified_rest_state(error)
      !! Same h_layer per layer (uniform across the column), but
      !! different rho_layer in each layer (stratified rest state).
      !! Because the column structure is identical at every (i, j),
      !! the horizontal pressure gradient must vanish.  Catches a
      !! per-layer integration bug that uniform-density misses.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: rho_k(NZ), max_u, max_v
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         ms%h_layer = H0
         ! Heaviest at bed (k=1), lightest at surface (k=NZ)
         rho_k = [pgf%rho0 + 2.0_wp, pgf%rho0, pgf%rho0 - 2.0_wp]
         do k = 1, NZ
            ms%rho_layer(:, :, k) = rho_k(k)
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(grid, ms, pgf, metrics)
         call ocean_pressure_force_compute(grid, metrics, pgf, ms)
         call ocean_pressure_force_apply(pgf, ms, DT)
         call map_out(ms, pgf, metrics)

         max_u = maxval(abs(ms%u_face_x_layer))
         max_v = maxval(abs(ms%v_face_y_layer))
         call check(error, max_u < 1.0e-12_wp, &
                    "stratified rest state: PGF injected u")
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0e-12_wp, &
                    "stratified rest state: PGF injected v")

      end block checks
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_stratified_rest_state

   subroutine test_horizontal_gradient(error)
      !! Density varies linearly with i in the bed layer; uniform in
      !! the layers above.  Analytic Mont-form prediction:
      !!
      !! Per-column:
      !!   p_edge(NZ+1) = 0
      !!   p_edge(k)    = p_edge(k+1) + g * rho_layer(k) * h_layer(k)
      !!   p_centre(k)  = 0.5 * (p_edge(k) + p_edge(k+1))
      !! Face acceleration:
      !!   dpdx_face(i, j, k) = -(1/rho0) * (p_centre(i, j, k)
      !!                                    - p_centre(i-1, j, k)) / dx
      !!
      !! For the gradient ONLY in layer 1: layer 1 sees a full PGF;
      !! layers 2 and NZ see zero (rho there is uniform across i, so
      !! their p_centre is the same in adjacent cells).
      !!
      !! Direction sanity: rho increases with i (east is denser),
      !! so pressure on the east side of the face is higher,
      !! du/dt < 0 (flow accelerates westward, away from high p).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DRHO_X = 0.05_wp     ! kg/m^3 per cell
      real(wp), parameter :: DX = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      integer :: i, j, k, nx, ny, i_c
      real(wp) :: expected_du_dt, max_diff, max_u_above
      checks: block

         call make_grid(grid, 8, 6, DX, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total
         i_c = nx/2

         ms%h_layer = H0
         ms%rho_layer = pgf%rho0
         ! Bed layer (k=1) carries the gradient
         do j = 1, ny
            do i = 1, nx
               ms%rho_layer(i, j, 1) = pgf%rho0 + DRHO_X*(real(i, wp) - real(i_c, wp))
            end do
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(grid, ms, pgf, metrics)
         call ocean_pressure_force_compute(grid, metrics, pgf, ms)
         call ocean_pressure_force_apply(pgf, ms, DT)
         call map_out(ms, pgf, metrics)

         ! Analytic Mont du/dt in the gradient layer at an interior face.
         ! Within layer 1: p_centre(1, i) = 0.5 * [sum from surface of
         ! g*rho_layer*h, bracketing layer 1's centre].  Only the
         ! gradient layer's rho varies with i, so:
         !   p_centre(1, i) - p_centre(1, i-1) = 0.5 * g * DRHO_X * H0
         !   (since DRHO_X is per unit i and h_layer = H0)
         ! du/dt = -(1/rho0) * (p_centre(1, i) - p_centre(1, i-1)) / dx
         expected_du_dt = -(1.0_wp/pgf%rho0)*0.5_wp*GRAVITY*DRHO_X*H0/DX

         max_diff = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               max_diff = max(max_diff, &
                              abs(ms%u_face_x_layer(i, j, 1) - DT*expected_du_dt))
            end do
         end do
         call check(error, max_diff < 1.0e-12_wp, &
                    "bed-layer PGF u increment deviates from analytic")
         if (allocated(error)) exit checks

         ! Sign sanity
         call check(error, expected_du_dt < 0.0_wp, &
                    "expected_du_dt sign wrong (rho east > rho west should push u west)")
         if (allocated(error)) exit checks

         ! Layers above the gradient layer must stay quiescent
         max_u_above = max(maxval(abs(ms%u_face_x_layer(:, :, 2))), &
                           maxval(abs(ms%u_face_x_layer(:, :, NZ))))
         call check(error, max_u_above < 1.0e-12_wp, &
                    "PGF leaked into upper layers from a bed-only gradient")

      end block checks
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_horizontal_gradient

   subroutine test_uniform_rho_varying_bathy(error)
      !! Regression for the FV-PGF bed→surface z_centre bug
      !! (commit 45421d4): uniform rho_layer over a horizontally-
      !! varying h_layer must produce zero dpdx_face on every face.
      !! Pre-fix the Jacobian was missing the `−g·ρ·dH/dx` term and
      !! injected `g·dH/dx` per face — `O(1 m/s²)` at a typical
      !! shelf-break slope.  Post-fix the Jacobian uses surface-
      !! relative `z_centre` and the residual is FP roundoff only.
      !!
      !! Test setup: h_layer is built from a per-column `b(i)` that
      !! varies linearly with `i`.  All other fields are uniform
      !! (no Coriolis, no shear, no stratification — just bathy and
      !! density).  Use OPGF_VARIANT_FV_LITE explicitly because the
      !! Montgomery variant doesn't apply the Jacobian and would
      !! always fail this test by design.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: B_FLAT = 200.0_wp
      real(wp), parameter :: B_SLOPE = 50.0_wp  ! per cell
      real(wp), parameter :: TOL = 1.0e-10_wp
      real(wp) :: max_dpdx, max_dpdy, b_i
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_FV_LITE
         nx = grid%nx_total
         ny = grid%ny_total

         ! Linear bathymetry: b(i) = B_FLAT + B_SLOPE * i (varies
         ! horizontally).  Distribute uniformly across NZ layers.
         do i = 1, nx
            b_i = B_FLAT + B_SLOPE*real(i, wp)
            do j = 1, ny
               do k = 1, NZ
                  ms%h_layer(i, j, k) = b_i/real(NZ, wp)
               end do
            end do
         end do

         ms%rho_layer = pgf%rho0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(grid, ms, pgf, metrics)
         call ocean_pressure_force_compute(grid, metrics, pgf, ms)
         !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
         max_dpdx = maxval(abs(pgf%dpdx_face%data))
         max_dpdy = maxval(abs(pgf%dpdy_face%data))
         call map_out(ms, pgf, metrics)

         call check(error, max_dpdx < TOL, &
                    "uniform_rho_varying_bathy: spurious dpdx at slope")
         if (allocated(error)) exit checks
         call check(error, max_dpdy < TOL, &
                    "uniform_rho_varying_bathy: spurious dpdy at slope")

      end block checks
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_uniform_rho_varying_bathy

end module test_ocean_pgf
