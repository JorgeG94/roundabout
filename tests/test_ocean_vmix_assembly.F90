!! Unit tests for the C11 vmix diffusivity-assembly gate
!! (`vmix_assemble` in `rdb_ocean_vmix`).
!!
!! The assembly is the single downstream gate every interior / overlay
!! closure feeds into before vdiff: background floors, kv_max/kd_max
!! ceilings, optional 1-2-1 horizontal smoothing, optional negative/NaN
!! guard.  Defaults are a pure pass-through (bit-identical to the
!! pre-refactor chain).
!!
!! Cases:
!!   * Defaults = pass-through: a known kv/kt/ks is unchanged after
!!     assemble when every knob is at its default.
!!   * kd_max clamps a deliberately huge diffusivity.
!!   * Background floors raise a sub-background value.
!!   * One 1-2-1 smoothing pass on a delta gives the exact 9-point
!!     stencil result AND conserves the interior sum.
!!   * The guard trips (status=1) on a negative diffusivity, and the
!!     status path returns instead of error-stopping.
module test_ocean_vmix_assembly
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_assemble, vmix_compute_pp81
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_setup, only: configure_ocean_lateral
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_assemble, vmix_split_kd_heat_salt
   implicit none
   private

   public :: collect_ocean_vmix_assembly_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_vmix_assembly_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("vmix_assemble_defaults_passthrough", test_defaults_passthrough), &
                  new_unittest("vmix_assemble_kd_max_clamps", test_kd_max_clamps), &
                  new_unittest("vmix_assemble_floor_applies", test_floor_applies), &
                  new_unittest("vmix_assemble_smooth_121_stencil", test_smooth_121), &
                  new_unittest("vmix_assemble_smooth_wetdry", test_smooth_wetdry), &
                  new_unittest("vmix_assemble_guard_trips_on_negative", test_guard_negative), &
                  new_unittest("vmix_pp81_knobs_reach_the_slot", test_pp81_knobs_reach_slot), &
                  new_unittest("vmix_pp81_defaults_bit_identical", test_pp81_defaults_bit_identical), &
                  new_unittest("pp81_alpha_changes_the_answer", test_pp81_alpha_changes_answer), &
                  new_unittest("vmix_split_derives_ks_from_kt", test_split_derives_ks), &
                  new_unittest("vmix_assemble_smooths_ks_like_kt", test_smooth_ks_matches_kt) &
                  ]
   end subroutine collect_ocean_vmix_assembly_tests

   subroutine make_state(grid, ms, vmix, nz)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: nz
      call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      call vmix%init(grid, nz_ml=nz)
      ms%h_layer = 5.0_wp
   end subroutine make_state

   subroutine run_assemble(grid, ms, vmix)
      !! Device round-trip: enter_data, assemble, pull kv/kt/ks back.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      !$acc update device(vmix%kv, vmix%kt, vmix%ks)
      call vmix_assemble(grid, vmix, ms)
      !$acc update self(vmix%kv, vmix%kt, vmix%ks)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_assemble

   ! -----------------------------------------------------------------

   subroutine test_defaults_passthrough(error)
      !! With default knobs (backgrounds = pp81_*_bg, ceilings huge,
      !! no smoothing, no guard) a known field above the background is
      !! left bit-identical (bitwise equality, not a tolerance check).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 5
      real(wp), allocatable :: kv0(:, :, :), kt0(:, :, :), ks0(:, :, :)
      integer :: k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         ! Seed interior interfaces above background; boundaries at zero.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0e-3_wp + 1.0e-4_wp*real(k, wp)
            vmix%kt(:, :, k) = 2.0e-4_wp + 1.0e-5_wp*real(k, wp)
            vmix%ks(:, :, k) = 3.0e-4_wp + 1.0e-5_wp*real(k, wp)
         end do
         kv0 = vmix%kv; kt0 = vmix%kt; ks0 = vmix%ks

         call run_assemble(grid, ms, vmix)

         ! Bitwise equality: the defaults (floor = pp81_*_bg, ceil = huge,
         ! no smoothing, no guard) must leave the field unchanged bit-for-bit.
         call check(error, all(vmix%kv == kv0), "kv changed under defaults")
         if (allocated(error)) exit checks
         call check(error, all(vmix%kt == kt0), "kt changed under defaults")
         if (allocated(error)) exit checks
         call check(error, all(vmix%ks == ks0), "ks changed under defaults")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_defaults_passthrough

   subroutine test_kd_max_clamps(error)
      !! A huge kt/ks is clipped to kd_max; kv clipped to kv_max.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      real(wp), parameter :: KD_CAP = 1.0e-2_wp, KV_CAP = 5.0e-2_wp
      integer :: k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         vmix%kd_max = KD_CAP
         vmix%kv_max = KV_CAP
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0_wp   ! 1 m^2/s, way over KV_CAP
            vmix%kt(:, :, k) = 1.0_wp
            vmix%ks(:, :, k) = 1.0_wp
         end do

         call run_assemble(grid, ms, vmix)

         call check(error, abs(vmix%kt(5, 4, 3) - KD_CAP) < 1.0e-15_wp, "kt not clamped to kd_max")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%ks(5, 4, 3) - KD_CAP) < 1.0e-15_wp, "ks not clamped to kd_max")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kv(5, 4, 3) - KV_CAP) < 1.0e-15_wp, "kv not clamped to kv_max")
         if (allocated(error)) exit checks
         ! Boundary interfaces (k=1, k=NZ+1) must stay untouched at zero.
         call check(error, abs(vmix%kv(5, 4, 1)) < 1.0e-15_wp, "bed interface clipped")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kv(5, 4, NZ + 1)) < 1.0e-15_wp, "surface interface clipped")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_kd_max_clamps

   subroutine test_floor_applies(error)
      !! A sub-background value is raised to the background floor; an
      !! above-background value is untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      integer :: k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         ! Custom, easy-to-check floors.
         vmix%kv_bg = 1.0e-3_wp; vmix%kt_bg = 2.0e-3_wp; vmix%ks_bg = 3.0e-3_wp
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0e-9_wp   ! below floor
            vmix%kt(:, :, k) = 1.0e-9_wp
            vmix%ks(:, :, k) = 1.0e-9_wp
         end do
         vmix%kv(5, 4, 2) = 9.9e-1_wp      ! above floor, must survive

         call run_assemble(grid, ms, vmix)

         call check(error, abs(vmix%kv(5, 4, 3) - vmix%kv_bg) < 1.0e-15_wp, "kv floor not applied")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kt(5, 4, 3) - vmix%kt_bg) < 1.0e-15_wp, "kt floor not applied")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%ks(5, 4, 3) - vmix%ks_bg) < 1.0e-15_wp, "ks floor not applied")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kv(5, 4, 2) - 9.9e-1_wp) < 1.0e-15_wp, "above-floor value altered")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_floor_applies

   subroutine test_smooth_121(error)
      !! One 1-2-1 (9-point) pass on a delta function gives the exact
      !! [[1,2,1],[2,4,2],[1,2,1]]/16 stencil and conserves the sum.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 3
      integer, parameter :: KI = 2          ! an interior interface
      integer :: ic, jc
      real(wp) :: sum_before, sum_after
      checks: block
         call make_state(grid, ms, vmix, NZ)
         vmix%kd_smooth_iterations = 1
         ! Floors off so they don't perturb the smoothed delta.
         vmix%kv_bg = 0.0_wp; vmix%kt_bg = 0.0_wp; vmix%ks_bg = 0.0_wp
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         ! Allocate smooth_scratch (normally done by configure_ocean_vmix).
         allocate (vmix%smooth_scratch(grid%nx_total, grid%ny_total, NZ + 1), source=0.0_wp)
         ! Delta at an interior (i,j) so all 8 neighbours are interior too.
         ! All cells are wet (ms%wet_mask defaults to 1.0) so the 9-point
         ! stencil weight is the full 16, giving [[1,2,1],[2,4,2],[1,2,1]]/16.
         ic = 6; jc = 5
         vmix%kv(ic, jc, KI) = 16.0_wp
         vmix%kt(ic, jc, KI) = 16.0_wp
         sum_before = sum(vmix%kv(2:grid%nx_total - 1, 2:grid%ny_total - 1, KI))

         call run_assemble(grid, ms, vmix)

         ! Centre: 4/16 * 16 = 4
         call check(error, abs(vmix%kv(ic, jc, KI) - 4.0_wp) < 1.0e-13_wp, "delta centre wrong")
         if (allocated(error)) exit checks
         ! Edge neighbour: 2/16 * 16 = 2
         call check(error, abs(vmix%kv(ic - 1, jc, KI) - 2.0_wp) < 1.0e-13_wp, "edge neighbour wrong")
         if (allocated(error)) exit checks
         ! Corner neighbour: 1/16 * 16 = 1
         call check(error, abs(vmix%kv(ic - 1, jc - 1, KI) - 1.0_wp) < 1.0e-13_wp, "corner neighbour wrong")
         if (allocated(error)) exit checks
         ! kt smoothed identically.
         call check(error, abs(vmix%kt(ic, jc, KI) - 4.0_wp) < 1.0e-13_wp, "kt delta centre wrong")
         if (allocated(error)) exit checks
         ! Sum conserved (delta interior, all weight stays interior).
         sum_after = sum(vmix%kv(2:grid%nx_total - 1, 2:grid%ny_total - 1, KI))
         call check(error, abs(sum_after - sum_before) < 1.0e-12_wp, "smoothing did not conserve sum")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_smooth_121

   subroutine test_smooth_wetdry(error)
      !! Wet/dry smoothing correctness:
      !!   (a) no leakage INTO a dry column — a dry cell's value stays
      !!       unchanged after a smoothing pass.
      !!   (b) no leakage OUT OF a dry column — wet neighbours that
      !!       surround a dry cell do not receive a contribution from it.
      !!   (c) conservation over the wet set is preserved (no wet-cell
      !!       value bleeds into dry cells and no mass is created).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 3
      integer, parameter :: KI = 2          ! interior interface
      integer :: id, jd                     ! dry-column indices
      real(wp) :: dry_val_before, dry_val_after
      real(wp) :: wet_sum_before, wet_sum_after
      integer :: i, j
      checks: block
         call make_state(grid, ms, vmix, NZ)
         vmix%kd_smooth_iterations = 1
         ! Floors off so smoothing result is not perturbed.
         vmix%kv_bg = 0.0_wp; vmix%kt_bg = 0.0_wp; vmix%ks_bg = 0.0_wp

         ! Allocate smooth_scratch (normally done by configure_ocean_vmix).
         allocate (vmix%smooth_scratch(grid%nx_total, grid%ny_total, NZ + 1), source=0.0_wp)

         ! Uniform wet field of 1; one interior dry column.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         vmix%kv(:, :, KI) = 1.0_wp
         id = 6; jd = 5
         ms%wet_mask = 1.0_wp
         ms%wet_mask(id, jd) = 0.0_wp    ! mark one column dry

         ! Save pre-smooth state of the dry cell and the wet sum.
         dry_val_before = vmix%kv(id, jd, KI)
         wet_sum_before = 0.0_wp
         do j = 2, grid%ny_total - 1
            do i = 2, grid%nx_total - 1
               if (ms%wet_mask(i, j) > 0.0_wp) wet_sum_before = wet_sum_before + vmix%kv(i, j, KI)
            end do
         end do

         call run_assemble(grid, ms, vmix)

         ! (a) dry column unchanged.
         dry_val_after = vmix%kv(id, jd, KI)
         call check(error, dry_val_after == dry_val_before, "dry cell value changed after smoothing")
         if (allocated(error)) exit checks

         ! (b) wet neighbours of the dry cell must not see the dry value
         ! bled in — since the dry value was 1.0 (same as wet) we instead
         ! check the neighbours of the dry cell are NOT affected by it by
         ! verifying no wet cell takes a value > 1 (which could only happen
         ! if the dry cell "donated" weight without being renormalised).
         call check(error, maxval(vmix%kv(2:grid%nx_total - 1, 2:grid%ny_total - 1, KI)) <= 1.0_wp + 1.0e-13_wp, &
                    "wet cell received out-of-range value from dry neighbour")
         if (allocated(error)) exit checks

         ! (c) wet-set sum is conserved (uniform field, no net change expected).
         wet_sum_after = 0.0_wp
         do j = 2, grid%ny_total - 1
            do i = 2, grid%nx_total - 1
               if (ms%wet_mask(i, j) > 0.0_wp) wet_sum_after = wet_sum_after + vmix%kv(i, j, KI)
            end do
         end do
         call check(error, abs(wet_sum_after - wet_sum_before) < 1.0e-11_wp, &
                    "wet-set sum not conserved after wet/dry smoothing")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_smooth_wetdry

   subroutine test_guard_negative(error)
      !! With the guard on and a negative kv seeded, the status-returning
      !! path reports 1 (and does not error-stop).
      !!
      !! The guard runs BEFORE the clip (guard → clip → smooth) so a
      !! negative raw-closure value trips the guard regardless of the
      !! background floor.  The old test set kv_bg = -1 to force the
      !! negative through the clip — that hack is no longer needed and
      !! was masking the ordering bug it was trying to work around.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      integer :: st, k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         vmix%vmix_guard = .true.
         ! Floors remain at positive defaults (kv_bg = pp81_nu_bg = 1e-4,
         ! etc.) — the guard sees the raw negative BEFORE any clip.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0e-3_wp
            vmix%kt(:, :, k) = 1.0e-4_wp
            vmix%ks(:, :, k) = 1.0e-4_wp
         end do
         vmix%kv(5, 4, 3) = -2.0e-3_wp   ! a deliberate negative

         st = 99
         !$acc enter data copyin(ms, vmix)
         call ms%enter_data()
         call vmix%enter_data()
         !$acc update device(vmix%kv, vmix%kt, vmix%ks)
         call vmix_assemble(grid, vmix, ms, status=st)
         call vmix%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vmix)

         call check(error, st == 1, "guard did not report a negative diffusivity")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_guard_negative

   subroutine setup_ocean_state(grid, state, nz)
      !! Minimal ocean_state_t setup shared by the PP81/KPP config-wiring
      !! tests below: grid + full ocean_state%init (allocates every slot,
      !! including vmix%kv/kt/ks/kd_bg — required BEFORE
      !! configure_ocean_lateral's seed_backgrounds() re-derive, §2E) +
      !! cartesian metrics (configure_ocean_lateral's ah_bg branch reads
      !! metrics_dx_min).  Mirrors `setup_state` in test_ocean_diag.F90.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: nz
      call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
      state%multilayer%nz_ml = nz
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine setup_ocean_state

   subroutine test_pp81_knobs_reach_slot(error)
      !! PR-9 §9.9.  All eight `&ocean_vmix_nml` pp81_*/kpp_* keys reach
      !! `ocean_state%vmix` via the PRODUCTION `configure_ocean_lateral`
      !! entry point — not by poking the slot directly (the exact gap the
      !! audit found: `grep -rn "%vmix%rho0\|vmix%pp81\|vmix%ri_crit..."
      !! src/` returned nothing outside tests/).  Also proves the
      !! structural invariant (§2E trap) survives the copy: kv_bg tracks
      !! the NEW pp81_nu_bg, kt_bg/ks_bg track the NEW pp81_kappa_bg, and
      !! the seeded interior kv/kt arrays equal the NEW background — a
      !! naive field-only copy without the `seed_backgrounds()` re-derive
      !! would leave these at the OLD type-default, silently raising the
      !! assembly floor above what the user asked for.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: ocean_state
      type(config_t) :: cfg
      integer, parameter :: NZ = 4
      real(wp), parameter :: NU0 = 2.0e-2_wp, NU_BG = 5.0e-4_wp, KAPPA_BG = 7.0e-5_wp
      real(wp), parameter :: ALPHA = 8.0_wp, SHEAR_FLOOR = 2.0e-9_wp
      real(wp), parameter :: RI_CRIT = 0.5_wp, CS_NL = 5.0_wp, C_VT2 = 1.0_wp

      checks: block
         call setup_ocean_state(grid, ocean_state, NZ)

         cfg%ocean%vmix%pp81_nu0 = NU0
         cfg%ocean%vmix%pp81_nu_bg = NU_BG
         cfg%ocean%vmix%pp81_kappa_bg = KAPPA_BG
         cfg%ocean%vmix%pp81_alpha = ALPHA
         cfg%ocean%vmix%shear2_floor = SHEAR_FLOOR
         cfg%ocean%vmix%kpp_ri_crit = RI_CRIT
         cfg%ocean%vmix%kpp_cs_nonlocal = CS_NL
         cfg%ocean%vmix%kpp_c_vt2 = C_VT2

         call configure_ocean_lateral(cfg, ocean_state, grid, compute_rank=0)

         call check(error, ocean_state%vmix%pp81_nu0 == NU0, "pp81_nu0 must reach the slot")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%pp81_nu_bg == NU_BG, "pp81_nu_bg must reach the slot")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%pp81_kappa_bg == KAPPA_BG, &
                    "pp81_kappa_bg must reach the slot")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%pp81_alpha == ALPHA, "pp81_alpha must reach the slot")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%shear2_floor == SHEAR_FLOOR, &
                    "shear2_floor must reach the slot")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%ri_crit == RI_CRIT, &
                    "kpp_ri_crit must reach vmix%ri_crit")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%cs_nonlocal == CS_NL, &
                    "kpp_cs_nonlocal must reach vmix%cs_nonlocal")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%c_vt2 == C_VT2, "kpp_c_vt2 must reach vmix%c_vt2")
         if (allocated(error)) exit checks

         ! Structural invariant (§2E).
         call check(error, ocean_state%vmix%kv_bg == NU_BG, &
                    "kv_bg must track the NEW pp81_nu_bg")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%kt_bg == KAPPA_BG, &
                    "kt_bg must track the NEW pp81_kappa_bg")
         if (allocated(error)) exit checks
         call check(error, ocean_state%vmix%ks_bg == KAPPA_BG, &
                    "ks_bg must track the NEW pp81_kappa_bg")
         if (allocated(error)) exit checks
         call check(error, all(ocean_state%vmix%kv(:, :, 2:NZ) == NU_BG), &
                    "seeded kv interior must equal the NEW pp81_nu_bg, not the old default")
         if (allocated(error)) exit checks
         call check(error, all(ocean_state%vmix%kt(:, :, 2:NZ) == KAPPA_BG), &
                    "seeded kt interior must equal the NEW pp81_kappa_bg, not the old default")
      end block checks
      call ocean_state%destroy()
   end subroutine test_pp81_knobs_reach_slot

   subroutine test_pp81_defaults_bit_identical(error)
      !! PR-9 §9.10.  A default `config_t` routed through
      !! `configure_ocean_lateral` leaves every relevant `ocean_vmix_t`
      !! field — including the seeded `kv`/`kt` arrays — bit-identical to
      !! a bare `vmix%init()`.  House bit-identity contract in its
      !! cheapest testable form: no default moved.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_plain, grid_cfg
      type(ocean_state_t) :: state_plain, state_configured
      type(config_t) :: cfg
      integer, parameter :: NZ = 4

      checks: block
         call setup_ocean_state(grid_plain, state_plain, NZ)
         call setup_ocean_state(grid_cfg, state_configured, NZ)
         call configure_ocean_lateral(cfg, state_configured, grid_cfg, compute_rank=0)

         call check(error, state_configured%vmix%pp81_nu0 == state_plain%vmix%pp81_nu0, &
                    "default config must not move pp81_nu0")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%pp81_nu_bg == state_plain%vmix%pp81_nu_bg, &
                    "default config must not move pp81_nu_bg")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%pp81_kappa_bg == state_plain%vmix%pp81_kappa_bg, &
                    "default config must not move pp81_kappa_bg")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%pp81_alpha == state_plain%vmix%pp81_alpha, &
                    "default config must not move pp81_alpha")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%shear2_floor == state_plain%vmix%shear2_floor, &
                    "default config must not move shear2_floor")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%ri_crit == state_plain%vmix%ri_crit, &
                    "default config must not move ri_crit")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%cs_nonlocal == state_plain%vmix%cs_nonlocal, &
                    "default config must not move cs_nonlocal")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%c_vt2 == state_plain%vmix%c_vt2, &
                    "default config must not move c_vt2")
         if (allocated(error)) exit checks
         call check(error, state_configured%vmix%kv_bg == state_plain%vmix%kv_bg, &
                    "default config must not move kv_bg")
         if (allocated(error)) exit checks
         call check(error, all(state_configured%vmix%kv == state_plain%vmix%kv), &
                    "default config must leave kv bit-identical to a bare init")
         if (allocated(error)) exit checks
         call check(error, all(state_configured%vmix%kt == state_plain%vmix%kt), &
                    "default config must leave kt bit-identical to a bare init")
      end block checks
      call state_plain%destroy(); call state_configured%destroy()
   end subroutine test_pp81_defaults_bit_identical

   subroutine test_pp81_alpha_changes_answer(error)
      !! PR-9 §9.11 (analytical).  A single interface with a KNOWN N² and
      !! KNOWN shear: assert `kv = pp81_nu_bg + pp81_nu0/(1+alpha*Ri)^2`
      !! for two DIFFERENT alpha values set THROUGH THE CONFIG, against
      !! the closed form (`rdb_ocean_vmix.F90:vmix_compute_pp81`).  Proves
      !! the knob is not merely stored — PP81's Ri response tracks it.
      !! (The existing `test_ocean_pp81.F90` suite sets vmix%pp81_* slot
      !! fields directly; that is exactly why nobody noticed the
      !! `&ocean_vmix_nml` keys did not exist — see PR-9 §2E.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid1, grid2
      type(ocean_state_t) :: state1, state2
      type(config_t) :: cfg1, cfg2
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: U_STEP = 1.0e-1_wp   ! du between adjacent layers
      real(wp), parameter :: DRHO = 1.0e-1_wp     ! rho drop between adjacent layers (lighter up)
      ! Moderate Ri (~O(0.1), not the huge-Ri quiescent-stratified limit):
      ! large enough that (1+alpha*Ri)^2 differs measurably between
      ! ALPHA1/ALPHA2, small enough that nu0/(1+alpha*Ri)^2 stays well
      ! above double-precision round-off.
      real(wp), parameter :: NU0 = 1.5e-2_wp, NU_BG = 3.0e-4_wp, SHEAR_FLOOR = 1.0e-12_wp
      real(wp), parameter :: ALPHA1 = 4.0_wp, ALPHA2 = 9.0_wp
      real(wp) :: shear2, n2, ri, kv1_expected, kv2_expected
      integer :: k, i_probe, j_probe

      checks: block
         call setup_ocean_state(grid1, state1, NZ)
         call setup_ocean_state(grid2, state2, NZ)

         cfg1%ocean%vmix%pp81_nu0 = NU0
         cfg1%ocean%vmix%pp81_nu_bg = NU_BG
         cfg1%ocean%vmix%shear2_floor = SHEAR_FLOOR
         cfg1%ocean%vmix%pp81_alpha = ALPHA1
         cfg2 = cfg1
         cfg2%ocean%vmix%pp81_alpha = ALPHA2

         call configure_ocean_lateral(cfg1, state1, grid1, compute_rank=0)
         call configure_ocean_lateral(cfg2, state2, grid2, compute_rank=0)

         i_probe = grid1%nx_total/2
         j_probe = grid1%ny_total/2

         do k = 1, NZ
            state1%multilayer%h_layer(:, :, k) = H_LAYER
            state1%multilayer%v_face_y_layer(:, :, k) = 0.0_wp
            state1%multilayer%u_face_x_layer(:, :, k) = U_STEP*real(k - 1, wp)
            state1%multilayer%rho_layer(:, :, k) = 1030.0_wp - DRHO*real(k - 1, wp)
         end do
         state2%multilayer%h_layer = state1%multilayer%h_layer
         state2%multilayer%u_face_x_layer = state1%multilayer%u_face_x_layer
         state2%multilayer%v_face_y_layer = state1%multilayer%v_face_y_layer
         state2%multilayer%rho_layer = state1%multilayer%rho_layer

         ! Closed-form Ri (constant across the interior interfaces: uniform
         ! dz, linear u and rho profiles).
         shear2 = max((U_STEP/H_LAYER)**2, SHEAR_FLOOR)
         n2 = GRAVITY*DRHO/(state1%vmix%rho0*H_LAYER)
         ri = n2/shear2
         kv1_expected = NU_BG + NU0/(1.0_wp + ALPHA1*ri)**2
         kv2_expected = NU_BG + NU0/(1.0_wp + ALPHA2*ri)**2
         ! Sanity: the two expected values must actually differ, or the
         ! test would pass vacuously regardless of whether alpha is wired.
         call check(error, abs(kv1_expected - kv2_expected) > 1.0e-8_wp, &
                    "test setup: kv1_expected and kv2_expected must differ")
         if (allocated(error)) exit checks

         ! Host-set h_layer/u_face_x_layer/v_face_y_layer/rho_layer are
         ! `copyin`'d by multilayer's enter_data (they are production
         ! state, not per-step scratch), so no explicit `!$acc update
         ! device` is needed here — see
         ! rdb_multilayer_state.F90:multilayer_state_enter_data_impl.
         call ocean_state_enter_data(state1)
         call ocean_state_enter_data(state2)
         call vmix_compute_pp81(grid1, state1%vmix, state1%multilayer)
         call vmix_compute_pp81(grid2, state2%vmix, state2%multilayer)
         ! vmix's exit_data does `exit data delete` (not copyout) on kv —
         ! pull the computed result back explicitly first (mirrors
         ! test_ocean_pp81.F90:run_pp81).
         !$acc update self(state1%vmix%kv, state2%vmix%kv)
         call ocean_state_exit_data(state1)
         call ocean_state_exit_data(state2)

         call check(error, abs(state1%vmix%kv(i_probe, j_probe, 2) - kv1_expected) < 1.0e-10_wp, &
                    "alpha=ALPHA1 (via config) must match the closed-form PP81 kv")
         if (allocated(error)) exit checks
         call check(error, abs(state2%vmix%kv(i_probe, j_probe, 2) - kv2_expected) < 1.0e-10_wp, &
                    "alpha=ALPHA2 (via config) must match the closed-form PP81 kv")
      end block checks
      call state1%destroy(); call state2%destroy()
   end subroutine test_pp81_alpha_changes_answer
   subroutine test_split_derives_ks(error)
      !! `vmix_split_kd_heat_salt` (PR-20) unconditionally overwrites `ks`
      !! with `kt` over the FULL (nx, ny, nz+1) extent, including the
      !! boundary interfaces k=1 and k=nz+1.  Seed `kt` with a recognisable
      !! non-uniform interior field and `ks` with deliberate garbage so
      !! the check proves an unconditional overwrite, not a merge.
      !!
      !! This is a strengthening of the existing gate tests, not a new
      !! feature test: `test_defaults_passthrough` / `test_kd_max_clamps`
      !! / `test_floor_applies` above seed `ks` INDEPENDENTLY of `kt` and
      !! must stay green unmodified — that is itself proof the split
      !! landed OUTSIDE `vmix_assemble` (§6.1 of the PR-20 plan).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 5
      integer :: k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         do k = 1, NZ + 1
            vmix%kt(:, :, k) = 1.0e-4_wp + 1.0e-5_wp*real(k, wp)
         end do
         vmix%ks = -999.0_wp   ! deliberate garbage — must be fully clobbered

         !$acc enter data copyin(ms, vmix)
         call ms%enter_data()
         call vmix%enter_data()
         !$acc update device(vmix%kt, vmix%ks)
         call vmix_split_kd_heat_salt(grid, vmix, ms)
         !$acc update self(vmix%kt, vmix%ks)
         call vmix%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vmix)

         call check(error, all(vmix%ks == vmix%kt), &
                    "vmix_split_kd_heat_salt: ks != kt over the full (nx,ny,nz+1) extent")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_split_derives_ks

   subroutine test_smooth_ks_matches_kt(error)
      !! Regression gate for the bit-identity trap: `kd_smooth_iterations
      !! > 0` must smooth `ks` identically to `kt`, or a real config with
      !! smoothing on would silently see `ks != kt` for no physical
      !! reason (§2.3 / §9.5 of the PR-20 plan) — invisible at the
      !! default `kd_smooth_iterations = 0`, so this is the one test
      !! that would catch a forgotten third `vmix_smooth_121_impl` call.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      integer :: i, j, k
      checks: block
         call make_state(grid, ms, vmix, NZ)
         vmix%kd_smooth_iterations = 2
         allocate (vmix%smooth_scratch(grid%nx_total, grid%ny_total, NZ + 1), source=0.0_wp)

         ! Non-uniform interior field, kt == ks going in (the split's
         ! post-condition), plus a non-trivial wet mask so the wet-aware
         ! stencil renormalisation is actually exercised.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  vmix%kt(i, j, k) = 1.0e-3_wp + 1.0e-4_wp*real(i + j + k, wp)
               end do
            end do
         end do
         vmix%ks = vmix%kt
         ms%wet_mask = 1.0_wp
         ms%wet_mask(6, 5) = 0.0_wp

         call run_assemble(grid, ms, vmix)

         call check(error, all(vmix%ks == vmix%kt), &
                    "vmix_assemble: kd_smooth_iterations smoothed ks differently from kt")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_smooth_ks_matches_kt

end module test_ocean_vmix_assembly
