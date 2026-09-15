!! Unit tests for the ocean salt/heat conservation budget feed.
!!
!! Verifies that:
!!   (1) `horiz_adv_accumulation` — `tracer_advect_zonal`/`_meridional`
!!       correctly fill `*_budget_horiz_adv` with the right sign/magnitude,
!!       and that the S/T dispatch gate is respected (passive tracers do
!!       NOT write into either array).
!!   (2) `salt_heat_error_closes` — algebraic fill-consistency identity on a
!!       closed basin with surface flux (total_change = Σsurface + Σadv).
!!   (3) `console_residual_closes` — a FAITHFUL RK2 step + manual average that
!!       exercises the real console arithmetic via the shared pure helpers
!!       (the 0.5 weight, the leading minus on `out`, the ρ scaling); the
!!       residual closes to round-off while the raw drift is O(src).
!!   (4) `passive_tracer_no_budget_leak` — an ideal-age passive tracer routes
!!       through the wrapper's `else` (no-`budget_adv`) branch; its advection
!!       must not leak into the salt/heat accumulators (proves the
!!       present(budget_adv)-ABSENT path is taken and isolated).
!!   (5) `windowed_advect_falls_back` — the `dt_tracer_advect_ratio>1` gate
!!       (`ocean_budget_is_active`) falls back to raw drift, plus a direct
!!       pin on the `ocean_budget_src`/`_out` weight+sign+ρ arithmetic.
!!   (6) `geothermal_folds_into_heat_src` — the geothermal bottom-heat source
!!       is folded into `heat_src` via `ocean_heat_src_sum` (heat-only); the
!!       heat residual closes only because of the fold.
!!   (7) `open_boundary_out_closes` — a uniform through-flow + linear tracer
!!       drives a genuinely NONZERO boundary transport `out`, pinning
!!       `ocean_budget_out` end-to-end on an OPEN-BC residual.
!!   (8) `prognostic_state_byte_identical` — drives the advect IMPLs WITH vs
!!       WITHOUT `budget_adv` on identical inputs; the `hTr` update must be
!!       bit-identical (a `hTr`-corrupting mutation in the budget-fill loop
!!       makes the WITH run diverge from the clean WITHOUT run).
!!   (9) `windowed_drain_with_surface_flux` — the windowed tracer-advect
!!       drain (`dt_tracer_advect_ratio > 1`) driven TOGETHER WITH an active
!!       surface flux, the combination no test covered before 2026-09-12.
!!       The budget must close mid-window and after the drain.
module test_ocean_conservation_salt_heat
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, &
                             continuity_zonal_flux, &
                             continuity_meridional_flux, &
                             tracer_advect_zonal, &
                             tracer_advect_meridional, &
                             tracer_advect_zonal_one_impl, &
                             tracer_advect_meridional_one_impl, &
                             continuity_tracer_step_split, &
                             continuity_tracer_drain, &
                             TR_MODE_ACCUMULATE
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers
   use rdb_ocean_console_stats, only: ocean_budget_src, ocean_budget_out, &
                                      ocean_budget_is_active, ocean_heat_src_sum
   use rdb_ocean_geothermal, only: ocean_geothermal_t, &
                                   ocean_geothermal_apply_tracers
   implicit none
   private

   public :: collect_ocean_conservation_salt_heat_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_conservation_salt_heat_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("horiz_adv_accumulation", &
                               test_horiz_adv_accumulation), &
                  new_unittest("salt_heat_error_closes", &
                               test_salt_heat_error_closes), &
                  new_unittest("console_residual_closes", &
                               test_console_residual_closes), &
                  new_unittest("passive_tracer_no_budget_leak", &
                               test_passive_tracer_no_budget_leak), &
                  new_unittest("windowed_advect_falls_back", &
                               test_windowed_advect_falls_back), &
                  new_unittest("geothermal_folds_into_heat_src", &
                               test_geothermal_folds_into_heat_src), &
                  new_unittest("open_boundary_out_closes", &
                               test_open_boundary_out_closes), &
                  new_unittest("prognostic_state_byte_identical", &
                               test_prognostic_byte_identical), &
                  new_unittest("windowed_drain_with_surface_flux", &
                               test_windowed_drain_with_surface_flux) &
                  ]
   end subroutine collect_ocean_conservation_salt_heat_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

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
      call pull_budgets(ms)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine pull_budgets(ms)
      !! Sync the device-resident conservation-budget accumulators back to
      !! the host before the host reads them.  The budget scratch is filled
      !! on-device by the advection/flux `do concurrent` kernels; the state's
      !! `exit_data` DELETEs it (production drains it via the budget manager's
      !! own `update self`), so a test that reads the raw accumulators must
      !! pull them itself.  `if_present` ⇒ inert when nothing is mapped (host
      !! build / un-mapped subtest) and on non-OpenACC compilers.
      type(multilayer_state_t), intent(inout) :: ms
      !$acc update self(ms%mass_budget_continuity, &
      !$acc&            ms%heat_budget_surface, ms%salt_budget_surface, &
      !$acc&            ms%heat_budget_geothermal, &
      !$acc&            ms%heat_budget_vert_adv, ms%salt_budget_vert_adv, &
      !$acc&            ms%heat_budget_vdiff, ms%salt_budget_vdiff, &
      !$acc&            ms%heat_budget_hdiff, ms%salt_budget_hdiff, &
      !$acc&            ms%heat_budget_horiz_adv, ms%salt_budget_horiz_adv, &
      !$acc&            ms%mass_budget_remap, ms%heat_budget_remap, &
      !$acc&            ms%salt_budget_remap) if_present
   end subroutine pull_budgets

   pure function area_sum(field, nx, ny, area_w) result(s)
      !! Host-side interior (ghost-excluded) area-weighted reduction of a
      !! per-cell budget array — mirrors `compute_total_tracer`'s stencil for
      !! the Cartesian test grid (uniform per-cell area `area_w = dx·dy`), so
      !! the result is the PRE-weight, PRE-ρ sum the `ocean_budget_*` helpers
      !! expect.  Kept in the test so `console_residual_closes` feeds the
      !! production helpers exactly the argument they get in the reporter.
      real(wp), intent(in) :: field(:, :, :)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: area_w
      real(wp) :: s
      s = sum(field(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*area_w
   end function area_sum

   ! ------------------------------------------------------------------
   ! §5.1  Analytical accumulation test
   ! ------------------------------------------------------------------

   subroutine test_horiz_adv_accumulation(error)
      !! Directly exercises `tracer_advect_zonal` (and `_meridional`)
      !! with a known mass-flux pattern and asserts the accumulator
      !! integrates the correct divergence with the right sign.
      !!
      !! Part A — uniform tracer under uniform flux:
      !!   mass_flux_x = M (constant), hTr = T0·h (uniform field)
      !!   ⇒ Tr_face = T0 everywhere ⇒ Tr_face_left = M·T0 uniform
      !!   ⇒ divergence per interior cell = 0 ⇒ budget_horiz_adv ≈ 0.
      !!   (CWC zero check + sign check.)
      !!
      !! Part B — linearly ramped mass flux:
      !!   mass_flux_x(i,j,k) = a·i (ramp), hTr = T0·h (uniform)
      !!   ⇒ Tr_face_left(i) = a·i·T0
      !!   ⇒ per-cell divergence (F(i+1)-F(i)) = a·T0 (constant)
      !!   ⇒ budget_horiz_adv(i,j,k) += -dt·a·T0·iareaT(i,j)
      !!   Assert to ~1e-13 relative.
      !!
      !! Part C — S/T gate: salt budget fills, heat budget fills, a
      !!   passive tracer (registered as 3rd) does NOT land in either.
      !!   We also assert heat_budget_horiz_adv is distinct from
      !!   salt_budget_horiz_adv (tracers carry different values).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: S0 = 35.0_wp    !! salinity (PSU)
      real(wp), parameter :: T0 = 12.0_wp    !! temperature (°C)
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: A_RAMP = 0.001_wp  !! flux ramp coefficient
      integer :: i, j, k, nx, ny, nx_phys, ny_phys
      real(wp) :: expected_adv, iareaT_val
      real(wp) :: max_salt_zero, max_heat_zero
      real(wp) :: max_salt_err, max_heat_err

      nx_phys = 12
      ny_phys = 8
      checks: block

         call make_grid(grid, nx_phys, ny_phys, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ms%h_layer = H0
         ms%tracers(ms%idx_salinity)%hTr = S0*H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0

         ! --- Part A: uniform mass flux → zero divergence → zero budget ---
         ! mass_flux uniform ⇒ Tr_face_left uniform ⇒ divergence 0 ⇒ budget ≈ 0
         ms%mass_flux_x_layer = A_RAMP
         ms%mass_flux_y_layer = A_RAMP

         call map_in(ms, ct)
         ! `enter_data` maps mass_flux with `create` (production recomputes it
         ! on-device every step), so the host-set flux above is NOT on the
         ! device yet — push it (else the -gpu=mem:separate kernel reads
         ! uninitialised device flux and the budget never fills).
         !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
         call tracer_advect_zonal(grid, metrics, ct, ms, DT)
         call tracer_advect_meridional(grid, metrics, ct, ms, DT)
         call map_out(ms, ct)

         max_salt_zero = maxval(abs(ms%salt_budget_horiz_adv))
         max_heat_zero = maxval(abs(ms%heat_budget_horiz_adv))
         call check(error, max_salt_zero < 1.0e-12_wp, &
                    "uniform-flux salt budget_horiz_adv should be ~0 (CWC)")
         if (allocated(error)) exit checks
         call check(error, max_heat_zero < 1.0e-12_wp, &
                    "uniform-flux heat budget_horiz_adv should be ~0 (CWC)")
         if (allocated(error)) exit checks

         ! --- Part B: ramped mass flux → known constant divergence ---
         ! Reset accumulators to isolate Part B
         ms%salt_budget_horiz_adv = 0.0_wp
         ms%heat_budget_horiz_adv = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = S0*H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0

         ! mass_flux_x(i,j,k) = A_RAMP·i ⇒ Tr_face_left(i) = A_RAMP·i·Tr
         ! per-cell divergence F(i+1)−F(i) = A_RAMP·Tr ⇒ budget += −dt·A_RAMP·Tr·iareaT
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%mass_flux_x_layer(i, j, k) = A_RAMP*real(i, wp)
               end do
            end do
         end do
         ms%mass_flux_y_layer = 0.0_wp

         iareaT_val = metrics%iareaT(nx_phys/2 + NGHOST, ny_phys/2 + NGHOST)

         call map_in(ms, ct)
         !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
         call tracer_advect_zonal(grid, metrics, ct, ms, DT)
         call map_out(ms, ct)

         ! Expected: budget_salt(i,j,k) = -dt · A_RAMP · S0 · iareaT
         ! (boundary cells fall back to first-order ⇒ skip i=1..3 and i=nx-1..nx)
         expected_adv = -DT*A_RAMP*S0*iareaT_val
         max_salt_err = 0.0_wp
         do k = 1, NZ
            do j = NGHOST + 1, ny - NGHOST
               do i = NGHOST + 3, nx - NGHOST - 2
                  max_salt_err = max(max_salt_err, &
                                     abs(ms%salt_budget_horiz_adv(i, j, k) - expected_adv))
               end do
            end do
         end do
         call check(error, max_salt_err <= 1.0e-13_wp*abs(expected_adv), &
                    "salt budget_horiz_adv ramp test: wrong magnitude or sign")
         if (allocated(error)) exit checks

         expected_adv = -DT*A_RAMP*T0*iareaT_val
         max_heat_err = 0.0_wp
         do k = 1, NZ
            do j = NGHOST + 1, ny - NGHOST
               do i = NGHOST + 3, nx - NGHOST - 2
                  max_heat_err = max(max_heat_err, &
                                     abs(ms%heat_budget_horiz_adv(i, j, k) - expected_adv))
               end do
            end do
         end do
         call check(error, max_heat_err <= 1.0e-13_wp*abs(expected_adv), &
                    "heat budget_horiz_adv ramp test: wrong magnitude or sign")
         if (allocated(error)) exit checks

         ! --- Part C: S and T accumulators are distinct ---
         ! (S0 ≠ T0 ⇒ the per-cell budgets differ)
         call check(error, abs(S0 - T0) > 1.0_wp, &
                    "test setup: S0 and T0 must differ for Part C to be meaningful")
         if (allocated(error)) exit checks
         call check(error, &
                    maxval(abs(ms%salt_budget_horiz_adv - ms%heat_budget_horiz_adv)) > 0.0_wp, &
                    "salt and heat budget_horiz_adv should differ (distinct tracer values)")

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_horiz_adv_accumulation

   ! ------------------------------------------------------------------
   ! §5.2  Closure test
   ! ------------------------------------------------------------------

   subroutine test_salt_heat_error_closes(error)
      !! Closed-basin (all-wall) kernel-level conservation closure test.
      !!
      !! Verifies the algebraic identity:
      !!   total_change = sum_budget_surface + sum_budget_adv
      !! over all interior cells (ghosts excluded), where both sides are
      !! in hTr·area·ρ₀ units.
      !!
      !! On a closed (all-wall) basin the advection budget telescopes to
      !! zero (net open-boundary flux = 0), so the identity reduces to
      !! `total_change = sum_budget_surface`, directly testing the surface
      !! fill convention.  Including the adv term makes the test general.
      !!
      !! Two tracer-advection + surface-flux stages are run per outer step
      !! (matching the production RK2 pattern without the rk2_average that
      !! the driver applies).  The budget arrays accumulate the raw sum over
      !! all stages; the algebraic identity holds exactly whether or not the
      !! rk2_average has been applied.
      !!
      !! Also asserts raw drift > 0 to confirm the surface flux ran.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 15.0_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: Q_HEAT = 50.0_wp   !! W/m² surface heat flux
      real(wp), parameter :: Q_SALT = 1.0e-6_wp  !! kg/(m²·s) surface salt flux
      real(wp), parameter :: TOL = 1.0e-10_wp
      integer, parameter :: N_OUTER = 4
      integer, parameter :: NX_PHYS = 8
      integer, parameter :: NY_PHYS = 6
      integer :: nx, ny, outer
      real(wp) :: ref_salt, ref_heat, total_salt, total_heat
      real(wp) :: sum_salt_bud_surface, sum_salt_bud_adv
      real(wp) :: sum_heat_bud_surface, sum_heat_bud_adv
      real(wp) :: err_salt, err_heat

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call sf%init(grid)

         ms%h_layer = H0
         ms%tracers(ms%idx_salinity)%hTr = S0*H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0
         ms%u_face_x_layer = 0.02_wp
         ms%v_face_y_layer = 0.01_wp
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_y_layer = 0.0_wp

         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

         ! Initial totals (host-side, before enter_data)
         ref_salt = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                        NGHOST + 1:ny - NGHOST, :)) &
                    *grid%dx*grid%dy*RHO_WATER
         ref_heat = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                           NGHOST + 1:ny - NGHOST, :)) &
                    *grid%dx*grid%dy*RHO_WATER

         !$acc enter data copyin(ms, sf, ct)
         call ms%enter_data()
         call sf%enter_data()
         call ct%enter_data()

         do outer = 1, N_OUTER
            ! Two stages per outer step (production RK2 pattern, no rk2_average)
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
            call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
            call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
         end do

         call pull_budgets(ms)
         call ct%exit_data()
         call ms%exit_data()
         call sf%exit_data()
         !$acc exit data delete(ms, sf, ct)

         ! Post-run totals (host-side, ghosts excluded)
         total_salt = sum(ms%tracers(ms%idx_salinity)%hTr(NGHOST + 1:nx - NGHOST, &
                                                          NGHOST + 1:ny - NGHOST, :)) &
                      *grid%dx*grid%dy*RHO_WATER
         total_heat = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                             NGHOST + 1:ny - NGHOST, :)) &
                      *grid%dx*grid%dy*RHO_WATER

         ! Budget term reductions (mirroring compute_total_tracer stencil)
         sum_salt_bud_surface = sum(ms%salt_budget_surface(NGHOST + 1:nx - NGHOST, &
                                                           NGHOST + 1:ny - NGHOST, :)) &
                                *grid%dx*grid%dy*RHO_WATER
         sum_heat_bud_surface = sum(ms%heat_budget_surface(NGHOST + 1:nx - NGHOST, &
                                                           NGHOST + 1:ny - NGHOST, :)) &
                                *grid%dx*grid%dy*RHO_WATER
         sum_salt_bud_adv = sum(ms%salt_budget_horiz_adv(NGHOST + 1:nx - NGHOST, &
                                                         NGHOST + 1:ny - NGHOST, :)) &
                            *grid%dx*grid%dy*RHO_WATER
         sum_heat_bud_adv = sum(ms%heat_budget_horiz_adv(NGHOST + 1:nx - NGHOST, &
                                                         NGHOST + 1:ny - NGHOST, :)) &
                            *grid%dx*grid%dy*RHO_WATER

         ! Algebraic identity: total_change = sum_budget_surface + sum_budget_adv
         err_salt = abs((total_salt - ref_salt) - (sum_salt_bud_surface + sum_salt_bud_adv))
         err_heat = abs((total_heat - ref_heat) - (sum_heat_bud_surface + sum_heat_bud_adv))

         ! (1) Budget identity closes to round-off
         call check(error, err_salt <= TOL*abs(ref_salt), &
                    "salt: total_change must equal sum_budget_surface + sum_budget_adv")
         if (allocated(error)) exit checks
         call check(error, err_heat <= TOL*abs(ref_heat), &
                    "heat: total_change must equal sum_budget_surface + sum_budget_adv")
         if (allocated(error)) exit checks

         ! (2) Surface budget is positive — proves the flux kernel actually ran
         ! and wrote into budget_surface.  We check the budget directly rather
         ! than a relative-drift of the total tracer mass (which is too small to
         ! exceed a fixed threshold when H0 is large relative to Q·dt).
         call check(error, sum_salt_bud_surface > 0.0_wp, &
                    "salt surface budget should be positive (surface flux ran)")
         if (allocated(error)) exit checks
         call check(error, sum_heat_bud_surface > 0.0_wp, &
                    "heat surface budget should be positive (surface flux ran)")

      end block checks
      call ct%destroy()
      call sf%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_salt_heat_error_closes

   ! ------------------------------------------------------------------
   ! §5.2b  Faithful-RK2 console-residual closure test
   ! ------------------------------------------------------------------

   subroutine test_console_residual_closes(error)
      !! Exercises the ACTUAL console conservation arithmetic — the 0.5 RK2
      !! weight, the leading minus on `out`, and the RHO_WATER scaling — as
      !! opposed to `salt_heat_error_closes` which only checks the fill
      !! consistency (an identity that stays true even if the console weight
      !! or sign is wrong).
      !!
      !! A FAITHFUL RK2 outer step is run: two stages from the same starting
      !! state, then a manual RK2 average of the tracers (mirroring
      !! `rk2_average_field_3d`, which is not public).  Because we average,
      !!   total_change = 0.5·Σ_stages  and  src = 0.5·Σ_surface
      !! so they balance ONLY if the console weight is 0.5.  Flip the console
      !! to 1.0 (or drop the leading minus on `out`) and check (a) fails.
      !!
      !! Console residual (see `emit_conservation_line` in rdb_console_stats):
      !!   residual = (total − ref) + out − src
      !!   src = +RK2_W·Σ(budget_surface)·dx·dy·RHO_WATER
      !!   out = −RK2_W·(Σ(budget_horiz_adv) + Σ(budget_hdiff))·dx·dy·RHO_WATER
      !! On a closed (all-wall, no bc arg) basin the advection budget
      !! telescopes to ~0, so `out ≈ 0` and the residual reduces to
      !! `(total − ref) − src` — directly testing that the console's src term
      !! (weight + sign + scaling) matches the state change.  hdiff is 0
      !! (kappa_h default 0); the term is kept in the formula for generality.
      !!
      !! Asserts:
      !!   (a) |residual| ≤ 1e-10·|ref|  for S and T          (the closure).
      !!   (b) |total_change| > 1e3·|residual|  for at least T (proves the
      !!       correction MATTERS — without it the Error would read the
      !!       O(src) drift, not the round-off residual).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 15.0_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: Q_HEAT = 200.0_wp   !! W/m² surface heat flux
      real(wp), parameter :: Q_SALT = 1.0e-4_wp  !! kg/(m²·s) surface salt flux
      real(wp), parameter :: TOL = 1.0e-10_wp
      integer, parameter :: NX_PHYS = 8
      integer, parameter :: NY_PHYS = 6
      integer :: nx, ny
      integer :: is, it
      real(wp), allocatable :: hS0(:, :, :), hT0(:, :, :)
      real(wp) :: ref_salt, ref_heat, total_salt, total_heat
      real(wp) :: change_salt, change_heat
      real(wp) :: src_salt, src_heat, out_salt, out_heat
      real(wp) :: res_salt, res_heat
      real(wp) :: scale_v, area_w

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call sf%init(grid)
         is = ms%idx_salinity
         it = ms%idx_temperature
         area_w = grid%dx*grid%dy      ! per-cell area (Cartesian, uniform)
         scale_v = area_w*RHO_WATER

         ms%h_layer = H0
         ms%tracers(is)%hTr = S0*H0
         ms%tracers(it)%hTr = T0*H0
         ! Small nonzero velocity so advection actually runs (closed basin ⇒
         ! interior advection telescopes to ~0, but the budget IS filled).
         ms%u_face_x_layer = 0.02_wp
         ms%v_face_y_layer = 0.01_wp
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_y_layer = 0.0_wp

         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

         ! Save the pre-step tracer fields (host) to form the manual RK2
         ! average AND to reference the interior total.
         allocate (hS0, source=ms%tracers(is)%hTr)
         allocate (hT0, source=ms%tracers(it)%hTr)

         ref_salt = sum(hS0(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*scale_v
         ref_heat = sum(hT0(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*scale_v

         !$acc enter data copyin(ms, sf, ct)
         call ms%enter_data()
         call sf%enter_data()
         call ct%enter_data()

         ! ---- FAITHFUL RK2 outer step ----
         ! Stage 1 (from state^n):
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
         call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
         ! Stage 2 (from the stage-1 state, as the production dyn step does):
         call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
         call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)

         call pull_budgets(ms)
         call ct%exit_data()
         call ms%exit_data()
         call sf%exit_data()
         !$acc exit data delete(ms, sf, ct)

         ! ---- Manual RK2 average of the tracers (mirrors rk2_average_field_3d) ----
         ! Literal 0.5 = the physics RK2 scheme.  The console weight lives in
         ! the production `ocean_budget_src`/`ocean_budget_out` helpers (called
         ! below), so the average here is INDEPENDENT of the console
         ! convention — that separation is what makes this a genuine guard:
         ! mutate RK2_STAGE_WEIGHT in the module and closure (a) breaks,
         ! because total_change stays 0.5·Σ_stages while src/out shift.
         ms%tracers(is)%hTr = 0.5_wp*(hS0 + ms%tracers(is)%hTr)
         ms%tracers(it)%hTr = 0.5_wp*(hT0 + ms%tracers(it)%hTr)

         ! ---- Post-average interior totals ----
         total_salt = sum(ms%tracers(is)%hTr(NGHOST + 1:nx - NGHOST, &
                                             NGHOST + 1:ny - NGHOST, :))*scale_v
         total_heat = sum(ms%tracers(it)%hTr(NGHOST + 1:nx - NGHOST, &
                                             NGHOST + 1:ny - NGHOST, :))*scale_v
         change_salt = total_salt - ref_salt
         change_heat = total_heat - ref_heat

         ! ---- Console src/out via the SAME production pure helpers ----
         ! The helpers own the weight + sign + ρ scaling; we feed the
         ! area-weighted interior budget sums PRE-weight, PRE-ρ (i.e. Σ·dx·dy).
         ! heat src folds surface + geothermal via ocean_heat_src_sum exactly
         ! like the reporter (geothermal is 0 here — no geo forcing — so it is
         ! a no-op fold, but it pins the reporter's argument structure).
         src_salt = ocean_budget_src(area_sum(ms%salt_budget_surface, nx, ny, area_w))
         src_heat = ocean_budget_src(ocean_heat_src_sum( &
                                     area_sum(ms%heat_budget_surface, nx, ny, area_w), &
                                     area_sum(ms%heat_budget_geothermal, nx, ny, area_w), &
                                     0.0_wp))
         out_salt = ocean_budget_out(area_sum(ms%salt_budget_horiz_adv, nx, ny, area_w), &
                                     area_sum(ms%salt_budget_hdiff, nx, ny, area_w))
         out_heat = ocean_budget_out(area_sum(ms%heat_budget_horiz_adv, nx, ny, area_w), &
                                     area_sum(ms%heat_budget_hdiff, nx, ny, area_w))

         ! residual = (total − ref) + out − src   (matches emit_conservation_line)
         res_salt = change_salt + out_salt - src_salt
         res_heat = change_heat + out_heat - src_heat

         ! (a) Console residual closes to round-off (S and T).
         call check(error, abs(res_salt) <= TOL*abs(ref_salt), &
                    "salt console residual must close to round-off")
         if (allocated(error)) exit checks
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "heat console residual must close to round-off")
         if (allocated(error)) exit checks

         ! (b) The correction matters — the raw drive (change ≈ src) is much
         !     larger than the residual.  Assert for heat (Q_HEAT is large
         !     enough to clear the ratio comfortably).
         call check(error, abs(change_heat) > 1.0e3_wp*abs(res_heat), &
                    "heat total_change must dominate the residual (correction matters)")
         if (allocated(error)) exit checks
         call check(error, abs(change_heat - src_heat) <= TOL*abs(ref_heat), &
                    "heat change must equal src on a closed basin (out ~ 0)")

      end block checks
      if (allocated(hS0)) deallocate (hS0)
      if (allocated(hT0)) deallocate (hT0)
      call ct%destroy()
      call sf%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_console_residual_closes

   ! ------------------------------------------------------------------
   ! §5.2c  Passive-tracer (else / present(budget_adv)-ABSENT) path
   ! ------------------------------------------------------------------

   subroutine test_passive_tracer_no_budget_leak(error)
      !! Exercises the wrapper's `else` branch — the `tracer_advect_*` call
      !! that passes NO `budget_adv` arg (so the impl runs the
      !! `.not. present(budget_adv)` path).  The S/T-only test states never
      !! reach it; a passive tracer does.
      !!
      !! Method: two identical advection runs, differing ONLY by the presence
      !! of a third (ideal-age) passive tracer.
      !!   Run A — S + T only.  Snapshot salt/heat budget_horiz_adv.
      !!   Run B — S + T + age (age carries a spatially-varying field).
      !! Asserts:
      !!   (1) the AGE tracer's hTr CHANGED under advection (proves the else
      !!       branch actually ran on it — the advect call took the no-budget
      !!       path), AND
      !!   (2) salt/heat budget_horiz_adv are BIT-IDENTICAL between A and B
      !!       (the passive tracer's advection did NOT leak into either S/T
      !!       accumulator — the dispatch routes age → else, never S/T).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      real(wp), parameter :: H0 = 20.0_wp
      real(wp), parameter :: S0 = 34.5_wp
      real(wp), parameter :: T0 = 10.0_wp
      real(wp), parameter :: U0 = 0.1_wp
      real(wp), parameter :: DT = 300.0_wp
      integer :: nx, ny, i, j, k, ia
      real(wp), allocatable :: age_before(:, :, :)
      real(wp), allocatable :: salt_bud_a(:, :, :), heat_bud_a(:, :, :)
      real(wp) :: max_age_change, max_salt_leak, max_heat_leak

      checks: block

         call make_grid(grid, 12, 10, 2000.0_wp, 2000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ms_a%nz_ml = NZ
         ms_b%nz_ml = NZ
         call ms_a%init(grid)                          ! S + T only
         call ms_b%init(grid, with_ideal_age=.true.)   ! S + T + age (idx 3)
         call ct_a%init(grid, nz_ml=NZ)
         call ct_b%init(grid, nz_ml=NZ)
         ia = ms_b%idx_age

         call check(error, ia > 0, "ideal-age tracer must register at idx_age > 0")
         if (allocated(error)) exit checks

         ! Identical S/T/h/velocity ICs on both states.
         ms_a%h_layer = H0
         ms_b%h_layer = H0
         ms_a%tracers(ms_a%idx_salinity)%hTr = S0*H0
         ms_b%tracers(ms_b%idx_salinity)%hTr = S0*H0
         ms_a%tracers(ms_a%idx_temperature)%hTr = T0*H0
         ms_b%tracers(ms_b%idx_temperature)%hTr = T0*H0
         ! Age: spatially-varying so advection actually moves it (else the
         ! "age changed" assertion could pass trivially for the wrong reason).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_b%tracers(ia)%hTr(i, j, k) = H0*(1.0_wp + 0.1_wp*real(i + j, wp))
               end do
            end do
         end do
         allocate (age_before, source=ms_b%tracers(ia)%hTr)

         ! Same smooth velocity field on both.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms_a%u_face_x_layer(i, j, k) = U0*sin(real(i, wp)*0.3_wp)
                  ms_b%u_face_x_layer(i, j, k) = ms_a%u_face_x_layer(i, j, k)
               end do
               do i = 1, nx
                  ms_a%v_face_y_layer(i, j, k) = U0*cos(real(j, wp)*0.4_wp)
                  ms_b%v_face_y_layer(i, j, k) = ms_a%v_face_y_layer(i, j, k)
               end do
            end do
         end do

         ! Compute fluxes then advect (both directions) on each state.
         call continuity_zonal_flux(grid, metrics, ct_a, ms_a, DT)
         call continuity_zonal_flux(grid, metrics, ct_b, ms_b, DT)

         call map_in(ms_a, ct_a)
         call map_in(ms_b, ct_b)

         call tracer_advect_zonal(grid, metrics, ct_a, ms_a, DT)
         call tracer_advect_zonal(grid, metrics, ct_b, ms_b, DT)
         call continuity_meridional_flux(grid, metrics, ct_a, ms_a, DT)
         call continuity_meridional_flux(grid, metrics, ct_b, ms_b, DT)
         call tracer_advect_meridional(grid, metrics, ct_a, ms_a, DT)
         call tracer_advect_meridional(grid, metrics, ct_b, ms_b, DT)

         call map_out(ms_a, ct_a)
         call map_out(ms_b, ct_b)

         ! (1) The age tracer moved ⇒ the else (no-budget) advect path ran.
         max_age_change = maxval(abs(ms_b%tracers(ia)%hTr - age_before))
         call check(error, max_age_change > 0.0_wp, &
                    "age tracer hTr must change under advection (else-branch path ran)")
         if (allocated(error)) exit checks

         ! (2) The passive tracer's presence must not perturb the S/T budgets.
         allocate (salt_bud_a, source=ms_a%salt_budget_horiz_adv)
         allocate (heat_bud_a, source=ms_a%heat_budget_horiz_adv)
         max_salt_leak = maxval(abs(ms_b%salt_budget_horiz_adv - salt_bud_a))
         max_heat_leak = maxval(abs(ms_b%heat_budget_horiz_adv - heat_bud_a))
         call check(error, max_salt_leak == 0.0_wp, &
                    "age advection must NOT leak into salt_budget_horiz_adv")
         if (allocated(error)) exit checks
         call check(error, max_heat_leak == 0.0_wp, &
                    "age advection must NOT leak into heat_budget_horiz_adv")
         if (allocated(error)) exit checks

         ! Sanity: the S/T budgets themselves are nonzero (advection ran).
         call check(error, maxval(abs(salt_bud_a)) > 0.0_wp, &
                    "salt_budget_horiz_adv should be nonzero (S/T advection ran)")

      end block checks
      if (allocated(age_before)) deallocate (age_before)
      if (allocated(salt_bud_a)) deallocate (salt_bud_a)
      if (allocated(heat_bud_a)) deallocate (heat_bud_a)
      call ct_a%destroy()
      call ct_b%destroy()
      call ms_a%destroy()
      call ms_b%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_passive_tracer_no_budget_leak

   ! ------------------------------------------------------------------
   ! §5.2d  Windowed-advection fallback gate
   ! ------------------------------------------------------------------

   subroutine test_windowed_advect_falls_back(error)
      !! Pins the console fallback: when `horiz_adv_budget_valid = .false.`
      !! -- i.e. some horizontal transport moved `hTr` without filling
      !! `*_budget_horiz_adv` -- the closed budget must be suppressed so the
      !! console reverts to raw drift rather than print a wrong "closed"
      !! residual.  Tests the production gate helper `ocean_budget_is_active`
      !! directly (the reporter uses the same call).
      !!
      !! NOTE (2026-09-12): the windowed tracer-advect path
      !! (`dt_tracer_advect_ratio > 1`) used to be the gate's only live
      !! producer and no longer is -- it is fully instrumented now, see
      !! `windowed_drain_with_surface_flux`.  The gate itself stays, for the
      !! next un-instrumented transport path.
      type(error_type), allocatable, intent(out) :: error

      checks: block
         ! Direct pin on the src/out helper arithmetic (weight + sign + ρ).
         ! The closed-basin `console_residual_closes` cannot distinguish the
         ! SIGN of `out` (it telescopes to ~0 there), so pin it point-blank
         ! here with a known input: src is +0.5·s·ρ, out is −0.5·(a+d)·ρ.
         call check(error, &
                    abs(ocean_budget_src(4.0_wp) - 0.5_wp*4.0_wp*RHO_WATER) &
                    <= 1.0e-9_wp*RHO_WATER, &
                    "ocean_budget_src must be +0.5·source·RHO_WATER")
         if (allocated(error)) exit checks
         call check(error, &
                    abs(ocean_budget_out(3.0_wp, 1.0_wp) + 0.5_wp*4.0_wp*RHO_WATER) &
                    <= 1.0e-9_wp*RHO_WATER, &
                    "ocean_budget_out must be −0.5·(adv+hdiff)·RHO_WATER (note sign)")
         if (allocated(error)) exit checks

         ! Fused path (budget complete) ⇒ a registered tracer is active.
         call check(error, ocean_budget_is_active(1, .true.), &
                    "registered tracer with complete budget must be active")
         if (allocated(error)) exit checks
         call check(error, ocean_budget_is_active(2, .true.), &
                    "registered tracer (idx 2) with complete budget must be active")
         if (allocated(error)) exit checks

         ! Windowed path (budget incomplete) ⇒ fall back to raw drift.
         call check(error,.not. ocean_budget_is_active(1, .false.), &
                    "salt must fall back to raw drift when horiz_adv budget invalid")
         if (allocated(error)) exit checks
         call check(error,.not. ocean_budget_is_active(2, .false.), &
                    "heat must fall back to raw drift when horiz_adv budget invalid")
         if (allocated(error)) exit checks

         ! Unregistered tracer ⇒ never active, regardless of budget validity.
         call check(error,.not. ocean_budget_is_active(0, .true.), &
                    "unregistered tracer (idx 0) must never be active")
         if (allocated(error)) exit checks
         call check(error,.not. ocean_budget_is_active(0, .false.), &
                    "unregistered tracer (idx 0) must never be active (invalid budget)")
      end block checks
   end subroutine test_windowed_advect_falls_back

   ! ------------------------------------------------------------------
   ! §5.2e  Geothermal bottom-heat source folded into heat src
   ! ------------------------------------------------------------------

   subroutine test_geothermal_folds_into_heat_src(error)
      !! Verifies the geothermal bottom-heat flux is accounted in the Heat
      !! Error (MF1): the reporter folds `heat_budget_geothermal` into
      !! `heat_src`, so an otherwise-unaccounted bottom source does NOT show
      !! up as a spurious leak.
      !!
      !! Runs the real geothermal kernel (no advection, no surface flux —
      !! isolates the geothermal limb) over a faithful RK2 outer step +
      !! manual RK2 average, then closes the heat budget using the SAME
      !! production helper the reporter uses, with the src fed
      !! surface(=0) + geothermal exactly as the reporter does:
      !!   residual = (total − ref) + out − src,  src = 0.5·(0 + Σgeo)·ρ
      !! Asserts (a) residual closes to round-off, and (b) dropping the
      !! geothermal term from src would leave an O(src) "leak" — i.e. the
      !! fold is load-bearing (change_heat is dominated by the geothermal
      !! source).  Salt is untouched by geothermal (heat-only) — asserted.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_geothermal_t) :: geo
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 15.0_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: Q_GEO = 0.5_wp   !! W/m² (exaggerated for signal)
      real(wp), parameter :: TOL = 1.0e-10_wp
      integer, parameter :: NX_PHYS = 8
      integer, parameter :: NY_PHYS = 6
      integer :: nx, ny, is, it
      real(wp), allocatable :: hT0(:, :, :), hS0(:, :, :)
      real(wp) :: ref_heat, ref_salt, total_heat, total_salt
      real(wp) :: change_heat, change_salt
      real(wp) :: src_heat, out_heat, res_heat
      real(wp) :: geo_sum, surf_sum, scale_v, area_w
      real(wp) :: src_surf_only, src_with_geo

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call geo%init(grid)
         is = ms%idx_salinity
         it = ms%idx_temperature
         area_w = grid%dx*grid%dy
         scale_v = area_w*RHO_WATER

         ms%h_layer = H0
         ms%tracers(is)%hTr = S0*H0
         ms%tracers(it)%hTr = T0*H0

         geo%enable = .true.
         geo%q_geo_const = Q_GEO

         allocate (hS0, source=ms%tracers(is)%hTr)
         allocate (hT0, source=ms%tracers(it)%hTr)
         ref_heat = sum(hT0(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*scale_v
         ref_salt = sum(hS0(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*scale_v

         !$acc enter data copyin(ms, geo)
         call ms%enter_data()

         ! Faithful RK2 outer step: geothermal each stage (fills hTr +
         ! heat_budget_geothermal, summed over both stages).
         call ocean_geothermal_apply_tracers(grid, geo, ms, DT)
         call ocean_geothermal_apply_tracers(grid, geo, ms, DT)

         call pull_budgets(ms)
         call ms%exit_data()
         !$acc exit data delete(ms, geo)

         ! Manual RK2 average (physics 0.5).
         ms%tracers(is)%hTr = 0.5_wp*(hS0 + ms%tracers(is)%hTr)
         ms%tracers(it)%hTr = 0.5_wp*(hT0 + ms%tracers(it)%hTr)

         total_heat = sum(ms%tracers(it)%hTr(NGHOST + 1:nx - NGHOST, &
                                             NGHOST + 1:ny - NGHOST, :))*scale_v
         total_salt = sum(ms%tracers(is)%hTr(NGHOST + 1:nx - NGHOST, &
                                             NGHOST + 1:ny - NGHOST, :))*scale_v
         change_heat = total_heat - ref_heat
         change_salt = total_salt - ref_salt

         ! Heat src via the SAME production helpers as the reporter:
         ! ocean_budget_src(ocean_heat_src_sum(surface, geothermal)).  out = 0
         ! (no advection / hdiff).  Calling ocean_heat_src_sum (rather than a
         ! local a+b) is what pins the reporter's geothermal FOLD: deleting the
         ! `+ geothermal_sum` in the helper fails assertion (a) below.
         surf_sum = area_sum(ms%heat_budget_surface, nx, ny, area_w)
         geo_sum = area_sum(ms%heat_budget_geothermal, nx, ny, area_w)
         src_heat = ocean_budget_src(ocean_heat_src_sum(surf_sum, geo_sum, 0.0_wp))
         out_heat = ocean_budget_out( &
                    area_sum(ms%heat_budget_horiz_adv, nx, ny, area_w), &
                    area_sum(ms%heat_budget_hdiff, nx, ny, area_w))
         res_heat = change_heat + out_heat - src_heat

         ! (a) Heat closure holds ONLY because geothermal is folded into src.
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "heat residual must close with geothermal folded into src")
         if (allocated(error)) exit checks

         ! (b) The fold helper is exact: ocean_heat_src_sum == surface + geo,
         !     and adding geothermal raises heat_src by EXACTLY the budget-src
         !     scaling of the geothermal part (pins the fold, not just a+b).
         call check(error, ocean_heat_src_sum(surf_sum, geo_sum, 0.0_wp) == surf_sum + geo_sum, &
                    "ocean_heat_src_sum must equal surface_sum + geothermal_sum")
         if (allocated(error)) exit checks
         src_surf_only = ocean_budget_src(ocean_heat_src_sum(surf_sum, 0.0_wp, 0.0_wp))
         src_with_geo = ocean_budget_src(ocean_heat_src_sum(surf_sum, geo_sum, 0.0_wp))
         call check(error, &
                    abs((src_with_geo - src_surf_only) - ocean_budget_src(geo_sum)) &
                    <= TOL*abs(ref_heat), &
                    "geothermal must raise heat_src by exactly budget_src(geo_sum)")
         if (allocated(error)) exit checks

         ! (c) The geothermal source is load-bearing: change_heat ≈ src_heat,
         !     and the geothermal reduction itself is nonzero.
         call check(error, ocean_budget_src(geo_sum) > 0.0_wp, &
                    "geothermal src contribution must be positive (source into ocean)")
         if (allocated(error)) exit checks
         call check(error, abs(change_heat) > 1.0e3_wp*max(abs(res_heat), tiny(1.0_wp)), &
                    "heat change must dominate residual (geothermal fold matters)")
         if (allocated(error)) exit checks

         ! (d) Salt is untouched — geothermal is heat-only.
         call check(error, abs(change_salt) == 0.0_wp, &
                    "salt must be unchanged by geothermal (heat-only source)")

      end block checks
      if (allocated(hS0)) deallocate (hS0)
      if (allocated(hT0)) deallocate (hT0)
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_geothermal_folds_into_heat_src

   ! ------------------------------------------------------------------
   ! §5.2f  Open-boundary through-flow: out ≠ 0 closure
   ! ------------------------------------------------------------------

   subroutine test_open_boundary_out_closes(error)
      !! The headline scenario: a NON-closed budget where the boundary
      !! transport `out` is genuinely nonzero — proving `ocean_budget_out`
      !! actually closes an open-BC residual (every other closure case is a
      !! closed basin where `out` telescopes to ~0 and can't distinguish the
      !! term).
      !!
      !! Setup drives `tracer_advect_zonal_one_impl` DIRECTLY with a UNIFORM
      !! through-flow mass flux `M > 0` at EVERY x-face (including the two
      !! domain-boundary faces — bypassing the wrapper's wall closure) and a
      !! LINEAR tracer `T(i) = a + b·i` (so the per-face tracer values differ
      !! ⇒ the interior divergence does NOT telescope to zero ⇒ net advective
      !! transport crosses the boundary).  Faithful 2-stage RK2 + manual
      !! average, no surface flux (src = 0), so:
      !!   residual = (total − ref) + out − 0,   out = ocean_budget_out(Σadv, 0)
      !! Because `budget_adv` IS exactly the per-stage hTr change, the identity
      !! closes by construction — but with a genuinely NONZERO `out`, which is
      !! the point: it pins `ocean_budget_out` end-to-end (weight + sign + ρ)
      !! against a real through-flow.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 5.0_wp
      real(wp), parameter :: DT = 200.0_wp
      real(wp), parameter :: MFLUX = 50.0_wp   !! uniform through-flow mass flux
      real(wp), parameter :: A_TR = 10.0_wp    !! tracer intercept
      real(wp), parameter :: B_TR = 2.0_wp     !! tracer slope (per i)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer, parameter :: NX_PHYS = 10
      integer, parameter :: NY_PHYS = 6
      integer :: nx, ny, i, j, k
      real(wp), allocatable :: hTr0(:, :, :), mfx(:, :, :), budget(:, :, :)
      real(wp), allocatable :: fxl(:, :, :), fxr(:, :, :)
      real(wp) :: ref_heat, total_heat, change_heat, out_heat, res_heat
      real(wp) :: area_w, scale_v

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 3000.0_wp, 3000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         area_w = grid%dx*grid%dy
         scale_v = area_w*RHO_WATER

         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = H0
         ! Linear tracer field T(i) = A_TR + B_TR·i ⇒ hTr = T·h.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     (A_TR + B_TR*real(i, wp))*H0
               end do
            end do
         end do

         ! Uniform through-flow: mass flux M at EVERY x-face, including the two
         ! boundary faces (i=1 and i=nx+1) so tracer enters/leaves the domain.
         allocate (mfx(nx + 1, ny, NZ), source=MFLUX)
         allocate (hTr0, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (fxl(nx + 1, ny, NZ), source=0.0_wp)
         allocate (fxr(nx + 1, ny, NZ), source=0.0_wp)
         allocate (budget(nx, ny, NZ), source=0.0_wp)

         ref_heat = sum(hTr0(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*scale_v

         ! Faithful 2-stage RK2 on the impl (uniform flux ⇒ same mfx each
         ! stage); budget accumulates over both stages via +=.
         ! GPU (-gpu=mem:separate): the impl's `do concurrent` gets NO implicit
         ! copies, so every argument must be device-present.  ms (h_layer +
         ! tracers%hTr) goes through the production deep-map; the local scratch
         ! (mfx / fxl / fxr / budget) is mapped explicitly.  metrics is already
         ! device-resident (make_cartesian_metrics).  All inert on host builds.
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(mfx) create(fxl, fxr) copyin(budget)
         call tracer_advect_zonal_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                           metrics%wet_T, ms%h_layer, &
                                           ms%tracers(ms%idx_temperature)%hTr, mfx, &
                                           fxl, fxr, budget_adv=budget)
         fxl = 0.0_wp
         fxr = 0.0_wp
         call tracer_advect_zonal_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                           metrics%wet_T, ms%h_layer, &
                                           ms%tracers(ms%idx_temperature)%hTr, mfx, &
                                           fxl, fxr, budget_adv=budget)
         !$acc update self(budget)
         !$acc exit data delete(mfx, fxl, fxr, budget)
         call ms%exit_data()
         !$acc exit data delete(ms)
         ! Manual RK2 average (physics 0.5).
         ms%tracers(ms%idx_temperature)%hTr = &
            0.5_wp*(hTr0 + ms%tracers(ms%idx_temperature)%hTr)

         total_heat = sum(ms%tracers(ms%idx_temperature)%hTr(NGHOST + 1:nx - NGHOST, &
                                                             NGHOST + 1:ny - NGHOST, :))*scale_v
         change_heat = total_heat - ref_heat

         ! out via the production helper (weight + sign + ρ), no hdiff, src = 0.
         out_heat = ocean_budget_out(area_sum(budget, nx, ny, area_w), 0.0_wp)
         res_heat = change_heat + out_heat

         ! (1) out is strictly nonzero — the through-flow populates the term.
         ! (A small FRACTION of the total heat, since the boundary transport
         ! over one step is tiny vs the resident column heat — so the honest
         ! guarantee is "nonzero AND dominates the residual", checked jointly
         ! with (3) below, not a fixed fraction-of-ref threshold.)
         call check(error, abs(out_heat) > 0.0_wp, &
                    "boundary transport out must be nonzero (through-flow populates it)")
         if (allocated(error)) exit checks

         ! (2) The residual closes to round-off WITH the nonzero out term.
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "open-BC residual (change + out) must close to round-off")
         if (allocated(error)) exit checks

         ! (3) out dominates the residual — the term is load-bearing, not noise.
         call check(error, abs(out_heat) > 1.0e3_wp*max(abs(res_heat), tiny(1.0_wp)), &
                    "out must dominate the residual (boundary term is load-bearing)")

      end block checks
      if (allocated(mfx)) deallocate (mfx)
      if (allocated(hTr0)) deallocate (hTr0)
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_open_boundary_out_closes

   ! ------------------------------------------------------------------
   ! §5.3  Byte-identity test
   ! ------------------------------------------------------------------

   subroutine test_prognostic_byte_identical(error)
      !! Prove the budget-fill loop does NOT perturb the prognostic hTr —
      !! by CONTRASTING the `present(budget_adv)`-PRESENT vs -ABSENT paths of
      !! the impls directly on identical inputs.
      !!
      !! The wrapper-level S/T tests can't do this: they always pass
      !! `budget_adv` for both salinity and temperature, so a mutation that
      !! wrote into `hTr` inside the budget loop would corrupt BOTH runs
      !! equally and slip past a state-vs-state comparison.  Here we drive
      !! `tracer_advect_zonal_one_impl` (and `_meridional`) TWICE on two
      !! copies of the same `hTr`: run A WITH `budget_adv`, run B WITHOUT.
      !! The prognostic update loop is shared and runs before the guarded
      !! fill, so the two must be bit-identical iff the fill never touches
      !! `hTr` — a `hTr`-corrupting mutation in the budget loop makes run A
      !! diverge from the clean run B and fails the ==0.0 assertion.
      !!
      !! Each impl mutates its face buffers (they double as the tracer-flux
      !! scratch), so runs A and B get their OWN face buffers.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 20.0_wp
      real(wp), parameter :: T0 = 10.0_wp
      real(wp), parameter :: U0 = 0.1_wp
      real(wp), parameter :: DT = 300.0_wp
      integer :: nx, ny, i, j, k
      real(wp), allocatable :: hTr0(:, :, :), hTr_a(:, :, :), hTr_b(:, :, :)
      real(wp), allocatable :: mfx(:, :, :), mfy(:, :, :)
      real(wp), allocatable :: fxl_a(:, :, :), fxr_a(:, :, :)
      real(wp), allocatable :: fxl_b(:, :, :), fxr_b(:, :, :)
      real(wp), allocatable :: fyl_a(:, :, :), fyr_a(:, :, :)
      real(wp), allocatable :: fyl_b(:, :, :), fyr_b(:, :, :)
      real(wp), allocatable :: budget(:, :, :)
      real(wp) :: max_diff

      checks: block

         call make_grid(grid, 12, 10, 2000.0_wp, 2000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ms%h_layer = H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0
         ! Spatially-varying velocity so the mass fluxes (and hence the tracer
         ! advection) are nontrivial across the domain.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U0*sin(real(i, wp)*0.3_wp)
               end do
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U0*cos(real(j, wp)*0.4_wp)
               end do
            end do
         end do

         ! Compute mass fluxes on the device, then snapshot the exact impl
         ! inputs to the host.  Under -gpu=mem:separate the `do concurrent`
         ! kernels get NO implicit copies, so ms/ct must be device-present for
         ! the flux compute, and the local scratch (mfx / hTr_* / fx*_* /
         ! budget) must be mapped explicitly around the impl calls.  All inert
         ! on host builds.
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(ct)
         call ct%enter_data()

         call continuity_zonal_flux(grid, metrics, ct, ms, DT)
         !$acc update self(ms%mass_flux_x_layer)
         allocate (mfx, source=ms%mass_flux_x_layer)

         allocate (hTr0, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hTr_a, source=hTr0)
         allocate (hTr_b, source=hTr0)
         allocate (fxl_a(nx + 1, ny, NZ), source=0.0_wp)
         allocate (fxr_a(nx + 1, ny, NZ), source=0.0_wp)
         allocate (fxl_b(nx + 1, ny, NZ), source=0.0_wp)
         allocate (fxr_b(nx + 1, ny, NZ), source=0.0_wp)
         allocate (budget(nx, ny, NZ), source=0.0_wp)

         ! Zonal: run A WITH budget_adv, run B WITHOUT — identical everything
         ! else, independent face buffers.
         !$acc enter data copyin(mfx, hTr_a, hTr_b, budget) &
         !$acc            create(fxl_a, fxr_a, fxl_b, fxr_b)
         call tracer_advect_zonal_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                           metrics%wet_T, ms%h_layer, hTr_a, mfx, &
                                           fxl_a, fxr_a, budget_adv=budget)
         call tracer_advect_zonal_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                           metrics%wet_T, ms%h_layer, hTr_b, mfx, &
                                           fxl_b, fxr_b)
         !$acc update self(hTr_a, hTr_b, budget)
         !$acc exit data delete(mfx, hTr_a, hTr_b, budget, fxl_a, fxr_a, fxl_b, fxr_b)

         max_diff = maxval(abs(hTr_a - hTr_b))
         call check(error, max_diff == 0.0_wp, &
                    "zonal hTr must be bit-identical WITH vs WITHOUT budget_adv fill")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(budget)) > 0.0_wp, &
                    "budget scratch must be nonzero after the WITH call (fill ran)")
         if (allocated(error)) exit checks

         ! Meridional: apply the zonal-updated hTr, recompute the meridional
         ! mass flux, then repeat the WITH-vs-WITHOUT contrast.  Reuse hTr_a as
         ! the common input for both runs (it equals hTr_b bit-for-bit here).
         ms%tracers(ms%idx_temperature)%hTr = hTr_a
         !$acc update device(ms%tracers(ms%idx_temperature)%hTr)
         call continuity_meridional_flux(grid, metrics, ct, ms, DT)
         !$acc update self(ms%mass_flux_y_layer)
         allocate (mfy, source=ms%mass_flux_y_layer)
         deallocate (hTr_b)
         allocate (hTr_b, source=hTr_a)
         allocate (fyl_a(nx, ny + 1, NZ), source=0.0_wp)
         allocate (fyr_a(nx, ny + 1, NZ), source=0.0_wp)
         allocate (fyl_b(nx, ny + 1, NZ), source=0.0_wp)
         allocate (fyr_b(nx, ny + 1, NZ), source=0.0_wp)
         budget = 0.0_wp

         !$acc enter data copyin(mfy, hTr_a, hTr_b, budget) &
         !$acc            create(fyl_a, fyr_a, fyl_b, fyr_b)
         call tracer_advect_meridional_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                                metrics%wet_T, ms%h_layer, hTr_a, mfy, &
                                                fyl_a, fyr_a, budget_adv=budget)
         call tracer_advect_meridional_one_impl(nx, ny, NZ, DT, metrics%iareaT, &
                                                metrics%wet_T, ms%h_layer, hTr_b, mfy, &
                                                fyl_b, fyr_b)
         !$acc update self(hTr_a, hTr_b, budget)
         !$acc exit data delete(mfy, hTr_a, hTr_b, budget, fyl_a, fyr_a, fyl_b, fyr_b)

         max_diff = maxval(abs(hTr_a - hTr_b))
         call check(error, max_diff == 0.0_wp, &
                    "meridional hTr must be bit-identical WITH vs WITHOUT budget_adv fill")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(budget)) > 0.0_wp, &
                    "meridional budget scratch must be nonzero after the WITH call")

      end block checks
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_prognostic_byte_identical

   pure subroutine copy_3d(nx, ny, nz_lev, src, dst)
      !! Device-side `dst = src` (explicit-shape, so NVHPC does not walk a
      !! descriptor per launch).  Used to take the top-of-outer-step save the
      !! manual RK2 average below needs, without a host round trip.
      integer, intent(in) :: nx, ny, nz_lev
      real(wp), intent(in) :: src(nx, ny, nz_lev)
      real(wp), intent(inout) :: dst(nx, ny, nz_lev)
      integer :: i, j, k
      do concurrent(k=1:nz_lev, j=1:ny, i=1:nx)
         dst(i, j, k) = src(i, j, k)
      end do
   end subroutine copy_3d

   pure subroutine rk2_avg_3d(nx, ny, nz_lev, f0, f)
      !! Device-side `f = 0.5·(f0 + f)` — the test's stand-in for the
      !! production `rk2_average` / `rk2_average_field_3d`, which are private
      !! to `rdb_ocean_dyn`.  The literal 0.5 is the RK2 SCHEME, deliberately
      !! independent of the console's `RK2_STAGE_WEIGHT` (see the note in
      !! `test_console_residual_closes`).
      integer, intent(in) :: nx, ny, nz_lev
      real(wp), intent(in) :: f0(nx, ny, nz_lev)
      real(wp), intent(inout) :: f(nx, ny, nz_lev)
      integer :: i, j, k
      do concurrent(k=1:nz_lev, j=1:ny, i=1:nx)
         f(i, j, k) = 0.5_wp*(f0(i, j, k) + f(i, j, k))
      end do
   end subroutine rk2_avg_3d

   subroutine test_windowed_drain_with_surface_flux(error)
      !! REGRESSION (2026-09-12): the windowed tracer-advection drain
      !! (`&ocean_vmix_nml dt_tracer_advect_ratio > 1`) TOGETHER WITH an
      !! active surface tracer flux.  That exact combination was uncovered:
      !! every windowed test ran adiabatically and every surface-flux test ran
      !! at ratio = 1, so nothing drove both.
      !!
      !! What went wrong.  In `TR_MODE_ACCUMULATE` the fused
      !! `tracer_advect_*` kernels are skipped, so nothing filled
      !! `ms%*_budget_horiz_adv`; the driver therefore declared the closed
      !! budget invalid and the console fell back to RAW DRIFT
      !! `(Q − Q0)/Q0`.  Raw drift cannot subtract a source, so a perfectly
      !! conservative run with `q_heat = −40 W/m²` reported the heat the flux
      !! legitimately added as a 5.6e-5 "leak" — 1e5x the corpus, on four
      !! shipped `acc_channel` namelists.
      !!
      !! The fix instruments the windowed path end to end: the drain
      !! sub-cycle AND both halves of the concentration hold accumulate into
      !! the same `*_budget_horiz_adv` arrays (see
      !! `DRAIN_BUDGET_POST_AVERAGE_WEIGHT` in `rdb_continuity`), so the
      !! closed budget is valid at ratio > 1 and the gate is gone.
      !!
      !! Two outer steps in accumulate mode with a surface heat + salt flux
      !! on every stage, a faithful RK2 average of `h` and `hTr` between
      !! them, then the drain.  Asserts:
      !!
      !!   (a) MID-WINDOW (after outer step 1, before the drain) the console
      !!       residual `(total − ref) + out − src` closes to round-off.  The
      !!       per-stage concentration hold `hTr := Tr·h` moves
      !!       `Σ areaT·hTr` by `Σ areaT·Tr·δh ≠ 0` for a non-uniform tracer,
      !!       so this only closes if the HOLD is instrumented too.
      !!   (b) mid-window `|out|` is materially non-zero — the accumulator is
      !!       carrying that correction rather than being all zeros (drop the
      !!       hold instrumentation and this is the assertion that dies
      !!       loudly instead of silently).
      !!   (c) at the END of the window, after the drain, the residual still
      !!       closes to round-off for BOTH salt and heat.
      !!   (d) the correction matters: the raw drift the old fall-back would
      !!       have printed is >= 1e3x the closed residual — i.e. this test
      !!       reproduces the 1e5x-over-budget signature if the fix is
      !!       reverted.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: S0 = 35.0_wp
      real(wp), parameter :: T0 = 15.0_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), parameter :: Q_HEAT = 200.0_wp   !! W/m² surface heat flux
      real(wp), parameter :: Q_SALT = 1.0e-4_wp  !! kg/(m²·s) surface salt flux
      real(wp), parameter :: TOL = 1.0e-10_wp
      integer, parameter :: NX_PHYS = 8
      integer, parameter :: NY_PHYS = 6
      integer, parameter :: RATIO = 2
      integer :: nx, ny, i, j, k
      integer :: is, it
      real(wp), allocatable :: h_save(:, :, :), hs_save(:, :, :), ht_save(:, :, :)
      real(wp) :: ref_salt, ref_heat, total_salt, total_heat
      real(wp) :: src_salt, src_heat, out_salt, out_heat
      real(wp) :: res_salt, res_heat, drift_heat
      real(wp) :: scale_v, area_w, ang

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         ct%windowed_advection = .true.   ! allocates uhtr/vhtr + drain scratch
         call ct%init(grid, nz_ml=NZ)
         call sf%init(grid)
         is = ms%idx_salinity
         it = ms%idx_temperature
         area_w = grid%dx*grid%dy
         scale_v = area_w*RHO_WATER

         ! NON-UNIFORM tracers and a DIVERGENT flow: both are load-bearing.
         ! A uniform tracer makes the concentration hold exactly
         ! content-conserving (Σ Tr·δh = Tr·Σ δh = 0) and assertion (b)
         ! degenerates to 0 > 0.
         ms%h_layer = H0
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ang = 0.7_wp*real(i, wp) + 0.4_wp*real(j, wp) + 0.3_wp*real(k, wp)
                  ms%tracers(is)%hTr(i, j, k) = (S0 + 0.8_wp*sin(ang))*H0
                  ms%tracers(it)%hTr(i, j, k) = (T0 + 2.0_wp*cos(ang))*H0
                  ms%u_face_x_layer(i, j, k) = 0.05_wp*sin(0.5_wp*real(i, wp))
                  ms%v_face_y_layer(i, j, k) = 0.04_wp*cos(0.3_wp*real(j, wp))
               end do
            end do
         end do
         ms%mass_flux_x_layer = 0.0_wp
         ms%mass_flux_y_layer = 0.0_wp

         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

         allocate (h_save(nx, ny, NZ), source=0.0_wp)
         allocate (hs_save(nx, ny, NZ), source=0.0_wp)
         allocate (ht_save(nx, ny, NZ), source=0.0_wp)

         ref_salt = area_sum(ms%tracers(is)%hTr, nx, ny, area_w)*RHO_WATER
         ref_heat = area_sum(ms%tracers(it)%hTr, nx, ny, area_w)*RHO_WATER

         !$acc enter data copyin(ms, sf, ct)
         call ms%enter_data()
         call sf%enter_data()
         call ct%enter_data()
         !$acc enter data create(h_save, hs_save, ht_save)

         ! ---- Outer step 1 of the window (accumulate; no drain) ----
         call window_outer_step(grid, metrics, ct, ms, sf, DT, nx, ny, NZ, is, it, &
                                h_save, hs_save, ht_save)

         ! (a)/(b) MID-WINDOW closure: hTr still carries the concentration
         ! hold, and `uhtr` still carries an unspent window.
         call window_totals(ms, nx, ny, is, it, area_w, total_salt, total_heat, &
                            src_salt, src_heat, out_salt, out_heat)
         res_salt = (total_salt - ref_salt) + out_salt - src_salt
         res_heat = (total_heat - ref_heat) + out_heat - src_heat
         call check(error, abs(res_salt) <= TOL*abs(ref_salt), &
                    "mid-window salt residual must close (concentration hold instrumented)")
         if (allocated(error)) exit checks
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "mid-window heat residual must close (concentration hold instrumented)")
         if (allocated(error)) exit checks
         call check(error, abs(out_heat) > 1.0e3_wp*abs(res_heat), &
                    "mid-window heat `out` must carry the concentration-hold correction")
         if (allocated(error)) exit checks

         ! ---- Outer step 2 of the window, then the drain ----
         call window_outer_step(grid, metrics, ct, ms, sf, DT, nx, ny, NZ, is, it, &
                                h_save, hs_save, ht_save)
         call continuity_tracer_drain(grid, metrics, ct, ms, RATIO)

         call window_totals(ms, nx, ny, is, it, area_w, total_salt, total_heat, &
                            src_salt, src_heat, out_salt, out_heat)

         !$acc exit data delete(h_save, hs_save, ht_save)
         call ct%exit_data()
         call ms%exit_data()
         call sf%exit_data()
         !$acc exit data delete(ms, sf, ct)

         res_salt = (total_salt - ref_salt) + out_salt - src_salt
         res_heat = (total_heat - ref_heat) + out_heat - src_heat
         drift_heat = total_heat - ref_heat

         ! (c) The closed budget still closes after the drain has spent the
         !     window -- with a surface flux running through it.
         call check(error, abs(res_salt) <= TOL*abs(ref_salt), &
                    "windowed-drain salt residual must close with a surface flux active")
         if (allocated(error)) exit checks
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "windowed-drain heat residual must close with a surface flux active")
         if (allocated(error)) exit checks

         ! (d) ... and the closure is doing real work: the raw drift the old
         !     fall-back printed is orders of magnitude larger.  This is the
         !     shipped-namelist signature (5.6e-5 drift vs 3e-14 residual).
         call check(error, abs(drift_heat) > 1.0e3_wp*abs(res_heat), &
                    "raw heat drift must dominate the closed residual (the fall-back was wrong)")

      end block checks
      if (allocated(h_save)) deallocate (h_save)
      if (allocated(hs_save)) deallocate (hs_save)
      if (allocated(ht_save)) deallocate (ht_save)
      call ct%destroy()
      call sf%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_windowed_drain_with_surface_flux

   subroutine window_outer_step(grid, metrics, ct, ms, sf, dt, nx, ny, nz_lev, is, it, &
                                h_save, hs_save, ht_save)
      !! One faithful RK2 outer step in `TR_MODE_ACCUMULATE` with the surface
      !! flux applied on BOTH stages, mirroring `run_stage` + `rk2_average` in
      !! `rdb_ocean_dyn`.  Averaging `h_layer` is not optional: the drain
      !! reconstructs `hprev = areaT·h_end + div(uhtr)` against the
      !! RK2-AVERAGED thickness (`uhtr` carries the 0.5 stage weight), so
      !! skipping the average breaks the telescoping the drain is built on.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: dt
      integer, intent(in) :: nx, ny, nz_lev, is, it
      real(wp), intent(inout) :: h_save(nx, ny, nz_lev)
      real(wp), intent(inout) :: hs_save(nx, ny, nz_lev)
      real(wp), intent(inout) :: ht_save(nx, ny, nz_lev)

      call copy_3d(nx, ny, nz_lev, ms%h_layer, h_save)
      call copy_3d(nx, ny, nz_lev, ms%tracers(is)%hTr, hs_save)
      call copy_3d(nx, ny, nz_lev, ms%tracers(it)%hTr, ht_save)

      call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                        tracer_mode=TR_MODE_ACCUMULATE)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, dt)
      call continuity_tracer_step_split(grid, metrics, ct, ms, dt, &
                                        tracer_mode=TR_MODE_ACCUMULATE)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, dt)

      call rk2_avg_3d(nx, ny, nz_lev, h_save, ms%h_layer)
      call rk2_avg_3d(nx, ny, nz_lev, hs_save, ms%tracers(is)%hTr)
      call rk2_avg_3d(nx, ny, nz_lev, ht_save, ms%tracers(it)%hTr)
   end subroutine window_outer_step

   subroutine window_totals(ms, nx, ny, is, it, area_w, total_salt, total_heat, &
                            src_salt, src_heat, out_salt, out_heat)
      !! Pull the device state back far enough to form the console's own
      !! `total` / `src` / `out` terms through the PRODUCTION helpers.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, is, it
      real(wp), intent(in) :: area_w
      real(wp), intent(out) :: total_salt, total_heat
      real(wp), intent(out) :: src_salt, src_heat, out_salt, out_heat

      call pull_budgets(ms)
      associate (hs => ms%tracers(is)%hTr, ht => ms%tracers(it)%hTr)
         !$acc update self(hs, ht) if_present
      end associate

      total_salt = area_sum(ms%tracers(is)%hTr, nx, ny, area_w)*RHO_WATER
      total_heat = area_sum(ms%tracers(it)%hTr, nx, ny, area_w)*RHO_WATER
      src_salt = ocean_budget_src(area_sum(ms%salt_budget_surface, nx, ny, area_w))
      src_heat = ocean_budget_src(ocean_heat_src_sum( &
                                  area_sum(ms%heat_budget_surface, nx, ny, area_w), &
                                  area_sum(ms%heat_budget_geothermal, nx, ny, area_w), &
                                  0.0_wp))
      out_salt = ocean_budget_out(area_sum(ms%salt_budget_horiz_adv, nx, ny, area_w), &
                                  area_sum(ms%salt_budget_hdiff, nx, ny, area_w))
      out_heat = ocean_budget_out(area_sum(ms%heat_budget_horiz_adv, nx, ny, area_w), &
                                  area_sum(ms%heat_budget_hdiff, nx, ny, area_w))
   end subroutine window_totals

end module test_ocean_conservation_salt_heat
