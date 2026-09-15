!! Unit tests for PR-28: `multilayer_state_t%register_passive_tracer`
!! and its first consumer, the pseudo-salt verification tracer
!! (`rdb_ocean_pseudo_salt`).
!!
!! Cases:
!!   * `register_passive_tracer_refuses_when_not_init` — the registry
!!     guard: calling on an uninitialised state returns idx=0.
!!   * `register_passive_tracer_refuses_when_locked` — the enter_data
!!     lock guard: registering after `enter_data` returns idx=0 and
!!     leaves the registry untouched (the mem:separate device-map
!!     foot-gun this PR closes).
!!   * `pseudo_salt_registers_and_seeds` — registration metadata
!!     (name/units/long_name, eos_coeff=0, budget_id=NONE) and the
!!     ghost-inclusive seed (hTr_ps == hTr_S everywhere, incl. ghosts).
!!   * `pseudo_salt_default_off_leaves_idx_zero` — knob-off bit-identity
!!     precondition: idx_pseudo_salt stays 0, registry unchanged.
!!   * `pseudo_salt_bc_ordering_gives_real_clamped_slot` — the §2.5
!!     registration-order contract: registering BEFORE
!!     `ocean_bc_state_init` gives pseudo-salt a real per-edge
!!     `clamped_tracer` slot.
!!   * `pseudo_salt_tracks_salinity_under_advection` — GPU canary: drive
!!     `continuity_tracer_step_split` for 50 steps on a divergent,
!!     sheared flow with a non-uniform S field; pseudo-salt (seeded to
!!     S) must track S to round-off — the "does a newly-registered
!!     passive tracer just work" proof, exercising the device-mapped
!!     array-of-DT registry indirection end to end.
!!   * `pseudo_salt_receives_identical_surface_salt_flux` — the surface
!!     salt-flux mirror is bit-identical to salinity's increment, and
!!     `salt_budget_surface` is unaffected by pseudo-salt's presence.
!!   * `pseudo_salt_no_budget_contribution` — the budget_id dispatch
!!     refactor is inert: heat/salt horizontal-advection budgets are
!!     bit-identical with vs without an extra NONE-budget tracer
!!     registered.
!!   * `pseudo_salt_conflict_predicates` — the three §5.6 fail-loud
!!     guard predicates, unit-tested without booting a config.
module test_ocean_pseudo_salt
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_NONE
   use rdb_ocean_pseudo_salt, only: ocean_pseudo_salt_register, ocean_pseudo_salt_seed, &
                                    pseudo_salt_conflicts_restore, &
                                    pseudo_salt_conflicts_ice, &
                                    pseudo_salt_needs_thermo_warning
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_continuity, only: continuity_t, continuity_tracer_step_split
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_pseudo_salt_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_pseudo_salt_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("register_passive_tracer_refuses_when_not_init", &
                               test_refuses_not_init), &
                  new_unittest("register_passive_tracer_refuses_when_locked", &
                               test_refuses_locked), &
                  new_unittest("pseudo_salt_registers_and_seeds", &
                               test_registers_and_seeds), &
                  new_unittest("pseudo_salt_default_off_leaves_idx_zero", &
                               test_default_off), &
                  new_unittest("pseudo_salt_bc_ordering_gives_real_clamped_slot", &
                               test_bc_ordering), &
                  new_unittest("pseudo_salt_tracks_salinity_under_advection", &
                               test_tracks_salinity_under_advection), &
                  new_unittest("pseudo_salt_receives_identical_surface_salt_flux", &
                               test_identical_surface_salt_flux), &
                  new_unittest("pseudo_salt_no_budget_contribution", &
                               test_no_budget_contribution), &
                  new_unittest("pseudo_salt_conflict_predicates", &
                               test_conflict_predicates) &
                  ]
   end subroutine collect_ocean_pseudo_salt_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   ! ------------------------------------------------------------------
   ! Registry API guards
   ! ------------------------------------------------------------------

   subroutine test_refuses_not_init(error)
      !! An uninitialised state (is_init = .false.) refuses registration.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: idx

      call make_grid(grid, 4, 4, 1.0_wp, 1.0_wp)
      idx = -1
      call ms%register_passive_tracer(grid, "dye", "1", "Dye tracer", idx)
      call check(error, idx == 0, "register_passive_tracer must refuse (idx=0) "// &
                 "when the registry is not initialised")
   end subroutine test_refuses_not_init

   subroutine test_refuses_locked(error)
      !! After enter_data (registry_locked = .true.), registration
      !! refuses AND leaves the registry untouched — the mem:separate
      !! device-map foot-gun this PR closes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: idx, n_before

      checks: block
         call make_grid(grid, 4, 4, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         n_before = size(ms%tracers)

         !$acc enter data copyin(ms)
         call ms%enter_data()
         call check(error, ms%registry_locked, "enter_data must set registry_locked")
         if (allocated(error)) exit checks

         idx = -1
         call ms%register_passive_tracer(grid, "dye", "1", "Dye tracer", idx)
         call check(error, idx == 0, "register_passive_tracer must refuse (idx=0) "// &
                    "once the registry is locked")
         if (allocated(error)) exit checks
         call check(error, size(ms%tracers) == n_before, &
                    "a refused registration must not grow the registry")
      end block checks
      call ms%exit_data()
      !$acc exit data delete(ms)
      if (.not. allocated(error)) then
         call check(error,.not. ms%registry_locked, "exit_data must clear registry_locked")
      end if
      call ms%destroy()
   end subroutine test_refuses_locked

   ! ------------------------------------------------------------------
   ! Registration + seed
   ! ------------------------------------------------------------------

   subroutine test_registers_and_seeds(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: nx, ny, i, j, k

      checks: block
         call make_grid(grid, 5, 4, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Non-uniform S, INCLUDING the ghost band (MOM6 isd:ied parity —
         ! the seed must cover it, not just the physical interior).
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (30.0_wp + 0.1_wp*real(i, wp) + 0.01_wp*real(j, wp) + real(k, wp))*10.0_wp
               end do
            end do
         end do

         call ocean_pseudo_salt_register(ms, grid)
         call check(error, ms%idx_pseudo_salt == size(ms%tracers), &
                    "idx_pseudo_salt must be the last (newest) slot")
         if (allocated(error)) exit checks
         call check(error, ms%idx_pseudo_salt == 3, "pseudo-salt should land at index 3 "// &
                    "(after default S=1, T=2)")
         if (allocated(error)) exit checks
         call check(error, ms%tracers(ms%idx_pseudo_salt)%eos_coeff == 0.0_wp, &
                    "pseudo-salt must be forced passive (eos_coeff = 0)")
         if (allocated(error)) exit checks
         call check(error, ms%tracers(ms%idx_pseudo_salt)%budget_id == TRACER_BUDGET_NONE, &
                    "pseudo-salt must register with budget_id = NONE")
         if (allocated(error)) exit checks
         call check(error, trim(ms%tracers(ms%idx_pseudo_salt)%name) == "pseudo_salt", &
                    "pseudo-salt name mismatch")
         if (allocated(error)) exit checks
         call check(error, trim(ms%tracers(ms%idx_pseudo_salt)%units) == "psu", &
                    "pseudo-salt units mismatch")
         if (allocated(error)) exit checks
         call check(error, &
                    trim(ms%tracers(ms%idx_pseudo_salt)%long_name) == "Pseudo salt passive tracer", &
                    "pseudo-salt long_name mismatch")
         if (allocated(error)) exit checks
         ! S/T keep their indices (§5.2 contract).
         call check(error, ms%idx_salinity == 1 .and. ms%idx_temperature == 2, &
                    "S/T indices must be unchanged by registration")
         if (allocated(error)) exit checks

         call ocean_pseudo_salt_seed(ms)
         call check(error, maxval(abs(ms%tracers(ms%idx_pseudo_salt)%hTr - &
                                      ms%tracers(ms%idx_salinity)%hTr)) < tiny(1.0_wp), &
                    "seed must set hTr_ps == hTr_S EXACTLY, including ghosts")
      end block checks
      call ms%destroy()
   end subroutine test_registers_and_seeds

   subroutine test_default_off(error)
      !! Without registering pseudo-salt, idx_pseudo_salt stays 0 and
      !! the registry is unchanged — the bit-identity precondition.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      checks: block
         call make_grid(grid, 4, 4, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call check(error, ms%idx_pseudo_salt == 0, &
                    "idx_pseudo_salt must be 0 when never registered")
         if (allocated(error)) exit checks
         call check(error, size(ms%tracers) == 2, "registry must stay S+T only")
      end block checks
      call ms%destroy()
   end subroutine test_default_off

   subroutine test_bc_ordering(error)
      !! Registering BEFORE ocean_bc_state_init (the §2.5 contract) gives
      !! pseudo-salt a real per-edge clamped_tracer slot: the size of
      !! bc%west%clamped_tracer must equal size(ms%tracers).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc

      checks: block
         call make_grid(grid, 4, 4, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ocean_pseudo_salt_register(ms, grid)

         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=size(ms%tracers))

         call check(error, allocated(bc%west%clamped_tracer), &
                    "clamped_tracer must be allocated when n_tracers > 0")
         if (allocated(error)) exit checks
         call check(error, size(bc%west%clamped_tracer) == size(ms%tracers), &
                    "clamped_tracer must size to the FULL registry, incl. pseudo-salt")
         if (allocated(error)) exit checks
         call check(error, ms%idx_pseudo_salt <= size(bc%west%clamped_tracer), &
                    "pseudo-salt's slot must be addressable in the OBC reservoir")
      end block checks
      call ms%destroy()
   end subroutine test_bc_ordering

   ! ------------------------------------------------------------------
   ! GPU canary: passive tracer rides the real advection kernel
   ! ------------------------------------------------------------------

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

   subroutine test_tracks_salinity_under_advection(error)
      !! Register + seed pseudo-salt to a NON-uniform S field, drive 50
      !! steps of the split continuity+tracer kernel on a divergent,
      !! sheared flow (test_split_mass_conservation's IC), and assert
      !! pseudo-salt tracks salinity to round-off.  Because both tracers
      !! start bit-identical and every subsequent kernel call applies the
      !! SAME deterministic per-tracer operator (same h_layer, same mass
      !! flux) to each registry slot independently, this is a strong
      !! canary for "does a newly-registered passive tracer ride the
      !! device-mapped array-of-DT registry for free" — exactly the
      !! seam contract this PR owns.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      integer, parameter :: N_STEPS = 50
      real(wp) :: layer_phase, max_dev, max_s
      integer :: i, j, k, nx, ny, step, idx_ps

      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            layer_phase = real(k - 1, wp)*PI/real(NZ, wp)
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = H_BASE*( &
                                                             35.0_wp + 0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp) + layer_phase))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp) + layer_phase)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                               sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp) + layer_phase)
               end do
            end do
         end do

         ! Register + seed AFTER S is filled (so the seed is non-trivial),
         ! BEFORE map_in/enter_data (the lock guard would refuse otherwise).
         call ocean_pseudo_salt_register(ms, grid)
         idx_ps = ms%idx_pseudo_salt
         call ocean_pseudo_salt_seed(ms)

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call continuity_tracer_step_split(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         max_dev = maxval(abs(ms%tracers(idx_ps)%hTr - ms%tracers(ms%idx_salinity)%hTr))
         max_s = maxval(abs(ms%tracers(ms%idx_salinity)%hTr))
         call check(error, max_dev < 1.0e-12_wp*max_s, &
                    "pseudo-salt must track salinity to round-off under advection")

      end block checks
      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_tracks_salinity_under_advection

   ! ------------------------------------------------------------------
   ! Surface salt-flux mirror
   ! ------------------------------------------------------------------

   subroutine test_identical_surface_salt_flux(error)
      !! Pseudo-salt's surface salt-flux increment must be EXACTLY
      !! (== 0.0_wp difference) equal to salinity's, and
      !! salt_budget_surface must be bit-identical to a run without
      !! pseudo-salt registered (the NOBUDGET twin must not leak into
      !! the budget accumulator).
      use rdb_ocean_surface_flux, only: ocean_surface_flux_apply_tracers
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_surface_flux_t) :: sf_a, sf_b
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_SALT_CONST = 5.0e-5_wp
      integer :: nx, ny, i, j, idx_ps
      real(wp) :: delta_ps, delta_s

      checks: block

         call make_grid(grid, 6, 4, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         ! ---- Run A: WITH pseudo-salt registered ----
         ms_a%nz_ml = NZ
         call ms_a%init(grid)
         ms_a%h_layer = H_LAYER
         ms_a%tracers(ms_a%idx_salinity)%hTr = 35.0_wp*H_LAYER
         call ocean_pseudo_salt_register(ms_a, grid)
         idx_ps = ms_a%idx_pseudo_salt
         call ocean_pseudo_salt_seed(ms_a)
         call sf_a%init(grid)
         call sf_a%set_surface_flux_const(0.0_wp, Q_SALT_CONST)
         ! Non-uniform Q_salt so the mirror is tested against a real
         ! per-column field, not a broadcast scalar.
         do j = 1, ny
            do i = 1, nx
               sf_a%Q_salt(i, j) = Q_SALT_CONST*(1.0_wp + 0.1_wp*real(i, wp))
            end do
         end do

         !$acc enter data copyin(ms_a, sf_a)
         call ms_a%enter_data()
         call sf_a%enter_data()
         call ocean_surface_flux_apply_tracers(grid, sf_a, ms_a, DT)
         associate (hS => ms_a%tracers(ms_a%idx_salinity)%hTr, &
                    hPs => ms_a%tracers(idx_ps)%hTr, &
                    budget => ms_a%salt_budget_surface)
            !$acc update self(hS, hPs, budget)
         end associate
         call sf_a%exit_data()
         call ms_a%exit_data()
         !$acc exit data delete(ms_a, sf_a)

         ! ---- Run B: WITHOUT pseudo-salt (reference for the budget check) ----
         ms_b%nz_ml = NZ
         call ms_b%init(grid)
         ms_b%h_layer = H_LAYER
         ms_b%tracers(ms_b%idx_salinity)%hTr = 35.0_wp*H_LAYER
         call sf_b%init(grid)
         call sf_b%set_surface_flux_const(0.0_wp, Q_SALT_CONST)
         do j = 1, ny
            do i = 1, nx
               sf_b%Q_salt(i, j) = Q_SALT_CONST*(1.0_wp + 0.1_wp*real(i, wp))
            end do
         end do
         !$acc enter data copyin(ms_b, sf_b)
         call ms_b%enter_data()
         call sf_b%enter_data()
         call ocean_surface_flux_apply_tracers(grid, sf_b, ms_b, DT)
         associate (budget => ms_b%salt_budget_surface)
            !$acc update self(budget)
         end associate
         call sf_b%exit_data()
         call ms_b%exit_data()
         !$acc exit data delete(ms_b, sf_b)

         ! (1) Exact match: pseudo-salt's increment == salinity's increment.
         delta_s = ms_a%tracers(ms_a%idx_salinity)%hTr(3, 2, NZ) - 35.0_wp*H_LAYER
         delta_ps = ms_a%tracers(idx_ps)%hTr(3, 2, NZ) - 35.0_wp*H_LAYER
         call check(error, delta_s == delta_ps, &
                    "pseudo-salt surface flux increment must EXACTLY equal salinity's")
         if (allocated(error)) exit checks

         ! (2) salt_budget_surface unaffected by pseudo-salt's presence.
         call check(error, maxval(abs(ms_a%salt_budget_surface - ms_b%salt_budget_surface)) &
                    < tiny(1.0_wp), &
                    "salt_budget_surface must be bit-identical with/without pseudo-salt")

      end block checks
      call sf_a%destroy(); call ms_a%destroy()
      call sf_b%destroy(); call ms_b%destroy()
   end subroutine test_identical_surface_salt_flux

   ! ------------------------------------------------------------------
   ! Budget-dispatch refactor bit-identity
   ! ------------------------------------------------------------------

   subroutine test_no_budget_contribution(error)
      !! The budget_id `select case` refactor (§5.3) is inert: running
      !! the split tracer-advection step WITH an extra NONE-budget
      !! (pseudo-salt) tracer registered must leave heat/salt
      !! horizontal-advection budgets bit-identical to a run without it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(ocean_metrics_t) :: metrics_a, metrics_b
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      integer :: nx, ny, i, j, k

      checks: block
         call make_grid(grid, 12, 8, 1.0_wp, 1.0_wp)

         ! ---- Run A: WITH pseudo-salt ----
         call make_cartesian_metrics(metrics_a, grid)
         ms_a%nz_ml = NZ
         call ms_a%init(grid)
         call ct_a%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_a%h_layer(i, j, k) = H_BASE
                  ms_a%tracers(ms_a%idx_salinity)%hTr(i, j, k) = H_BASE*(35.0_wp + 0.1_wp*real(i, wp))
                  ms_a%tracers(ms_a%idx_temperature)%hTr(i, j, k) = H_BASE*(15.0_wp + 0.1_wp*real(j, wp))
               end do
               do i = 1, nx + 1
                  ms_a%u_face_x_layer(i, j, k) = 0.1_wp*sin(real(i, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms_a%v_face_y_layer(i, j, k) = 0.1_wp*cos(real(j, wp))
               end do
            end do
         end do
         call ocean_pseudo_salt_register(ms_a, grid)
         call ocean_pseudo_salt_seed(ms_a)
         ms_a%heat_budget_horiz_adv = 0.0_wp
         ms_a%salt_budget_horiz_adv = 0.0_wp

         call map_in(ms_a, ct_a)
         call continuity_tracer_step_split(grid, metrics_a, ct_a, ms_a, DT)
         associate (hb => ms_a%heat_budget_horiz_adv, sb => ms_a%salt_budget_horiz_adv)
            !$acc update self(hb, sb)
         end associate
         call map_out(ms_a, ct_a)

         ! ---- Run B: WITHOUT pseudo-salt (identical IC otherwise) ----
         call make_cartesian_metrics(metrics_b, grid)
         ms_b%nz_ml = NZ
         call ms_b%init(grid)
         call ct_b%init(grid, nz_ml=NZ)
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_b%h_layer(i, j, k) = H_BASE
                  ms_b%tracers(ms_b%idx_salinity)%hTr(i, j, k) = H_BASE*(35.0_wp + 0.1_wp*real(i, wp))
                  ms_b%tracers(ms_b%idx_temperature)%hTr(i, j, k) = H_BASE*(15.0_wp + 0.1_wp*real(j, wp))
               end do
               do i = 1, nx + 1
                  ms_b%u_face_x_layer(i, j, k) = 0.1_wp*sin(real(i, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms_b%v_face_y_layer(i, j, k) = 0.1_wp*cos(real(j, wp))
               end do
            end do
         end do
         ms_b%heat_budget_horiz_adv = 0.0_wp
         ms_b%salt_budget_horiz_adv = 0.0_wp

         call map_in(ms_b, ct_b)
         call continuity_tracer_step_split(grid, metrics_b, ct_b, ms_b, DT)
         associate (hb => ms_b%heat_budget_horiz_adv, sb => ms_b%salt_budget_horiz_adv)
            !$acc update self(hb, sb)
         end associate
         call map_out(ms_b, ct_b)

         call check(error, maxval(abs(ms_a%heat_budget_horiz_adv - ms_b%heat_budget_horiz_adv)) &
                    < tiny(1.0_wp), &
                    "heat_budget_horiz_adv must be bit-identical with/without pseudo-salt")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms_a%salt_budget_horiz_adv - ms_b%salt_budget_horiz_adv)) &
                    < tiny(1.0_wp), &
                    "salt_budget_horiz_adv must be bit-identical with/without pseudo-salt")

      end block checks
      call ct_a%destroy(); call ms_a%destroy(); call destroy_cartesian_metrics(metrics_a)
      call ct_b%destroy(); call ms_b%destroy(); call destroy_cartesian_metrics(metrics_b)
   end subroutine test_no_budget_contribution

   ! ------------------------------------------------------------------
   ! Fail-loud guard predicates
   ! ------------------------------------------------------------------

   subroutine test_conflict_predicates(error)
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, pseudo_salt_conflicts_restore(.true., .true.), &
                    "pseudo-salt + SSS restore must conflict")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_conflicts_restore(.true., .false.), &
                    "pseudo-salt alone must not conflict with restore")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_conflicts_restore(.false., .true.), &
                    "restore alone (pseudo-salt off) must not conflict")
         if (allocated(error)) exit checks

         call check(error, pseudo_salt_conflicts_ice(.true., .true.), &
                    "pseudo-salt + sea-ice must conflict")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_conflicts_ice(.true., .false.), &
                    "pseudo-salt alone must not conflict with ice")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_conflicts_ice(.false., .true.), &
                    "ice alone (pseudo-salt off) must not conflict")
         if (allocated(error)) exit checks

         call check(error, pseudo_salt_needs_thermo_warning(.true., .false.), &
                    "pseudo-salt with thermo off must warn")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_needs_thermo_warning(.true., .true.), &
                    "pseudo-salt with thermo on must not warn")
         if (allocated(error)) exit checks
         call check(error,.not. pseudo_salt_needs_thermo_warning(.false., .false.), &
                    "pseudo-salt off must never warn")
      end block checks
   end subroutine test_conflict_predicates

end module test_ocean_pseudo_salt
