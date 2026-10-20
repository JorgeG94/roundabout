!! Tripolar north-fold regression tests for C-grid EVP sea-ice dynamics
!! (A6, #164).
!!
!! **What this guards.**  Before A6, `ice_evp_step` pinned the north ghost
!! band to land under EVERY edge policy including a tripolar fold
!! (`land_n = bc%has_north .and. .not. bc%periodic_y`, true under a fold)
!! -- `evp_build_masks_impl` walled the fold row off, so EVP's ice
!! momentum never crossed the seam even though the category state (mass,
!! concentration) WAS folded (`ocean_fold_wrap_centre_flat`).  A6 (a) lets
!! the fold row flow (`land_n` now excludes `bc%north_fold`) and (b) folds
!! `ui`/`vi` (true-vector, `ocean_fold_wrap_face_pair`) every subcycle
!! plus the restart-carried stresses once before the loop
!! (`fold_north_centre`/`fold_north_corner`).
!!
!! Harness: `ocean_test_metrics::make_tripolar_metrics` (Murray bipolar
!! cap, periodic west/east, `phi_join=60`) + raw test-allocated EVP
!! arrays, mirroring `test_ocean_ice_evp`'s harness.  No land (`wet_T`
!! defaults all-wet) -- "uniform ice over all wet cells" below means
!! every cell.  GPU mem:separate discipline throughout: unconditional
!! `!$acc enter data copyin(...)` before every kernel call, `!$acc update
!! self` before every host assertion, `exit data` at teardown (canonical
!! template: `test_open_boundary_out_closes`,
!! `test_ocean_conservation_salt_heat.F90`).
!!
!! Cases:
!!   * `evp_fold_halo_mirrors_exactly` -- one outer EVP step; assert the
!!     north-halo `ui`/`vi` rows are the EXACT (bitwise) fold image of
!!     their physical source; the stresses (`str_d`/`str_t` in the halo,
!!     `str_s` on the fold row) are fold-symmetric to round-off.
!!   * `evp_fold_seam_row_antisymmetric` -- the on-line `vi` projection at
!!     the fold row is bitwise antisymmetric, and the cross-seam ice mass
!!     flux (uniform `mi_v`, true-vector `vi`, fold-symmetric `dx_cv`)
!!     telescopes to round-off -- the EVP analogue of the ocean's B1 leak
!!     check.
!!   * `evp_fold_run_no_seam_stripe` -- 48 outer steps; every entry stays
!!     finite, bounded, no more CFL truncation than the wall run, and
!!     no seam-row stripe artifact
!!     (the fold-row inter-row |u| jump is not an outlier against the
!!     rest of the domain).  Then the LIVE-change proof: an identical
!!     setup with `north="wall"` instead of the fold must give
!!     `max|vi|` == 0 on the boundary face (walled), while the fold run
!!     gives > 0 (the row actually flows) -- confirms `land_n`'s flip is
!!     what is under test, not a vacuous no-op.
!!   Deviation from the plan: no `ice_transport_step` leg here (transport
!!   across the fold is already gated by
!!   `test_ocean_ice_fold::ice_fold_seam_conserves_and_bounded`; wiring a
!!   full `ocean_sea_ice_t` category-state transport run into this raw-
!!   array EVP harness would roughly double the file for no new seam
!!   coverage). The coarse test grid's near-pole cells trip the CFL clip
!!   under either north boundary, so truncation is compared against the
!!   wall run rather than asserted zero.
module test_ocean_ice_evp_fold
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, OMEGA_EARTH, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_tripolar_metrics, destroy_cartesian_metrics
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_evp, only: ice_evp_dynamics, ice_evp_params_t
   use rdb_ice_state, only: evp_workspace_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_PERIODIC, &
                                       OBC_TRIPOLAR_FOLD, ocean_bc_validate_periodic, &
                                       ocean_bc_validate_fold
   implicit none
   private

   public :: collect_ocean_ice_evp_fold_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 16, NYP = 10
      !! NXP even (no self-conjugate v column).
   real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 50.0_wp
   real(wp), parameter :: DLON = 360.0_wp/real(NXP, wp)
   real(wp), parameter :: DLAT = 3.0_wp
   real(wp), parameter :: REARTH = 6.378e6_wp
   real(wp), parameter :: PHI_JOIN = 60.0_wp, LON_POLE = 0.0_wp

   real(wp), parameter :: MIS_CONST = 900.0_wp, MICE_CONST = 900.0_wp, CI_CONST = 1.0_wp
   real(wp), parameter :: TAU_AX_CONST = 0.1_wp
   real(wp), parameter :: DT_SLOW = 1800.0_wp

contains

   subroutine collect_ocean_ice_evp_fold_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("evp_fold_halo_mirrors_exactly", test_halo_mirrors), &
                  new_unittest("evp_fold_seam_row_antisymmetric", test_seam_antisym), &
                  new_unittest("evp_fold_run_no_seam_stripe", test_no_seam_stripe) &
                  ]
   end subroutine collect_ocean_ice_evp_fold_tests

   ! -----------------------------------------------------------------
   ! Shared setup
   ! -----------------------------------------------------------------

   subroutine make_grid_metrics(grid, metrics)
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
      call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                 REARTH, PHI_JOIN, LON_POLE)
   end subroutine make_grid_metrics

   subroutine make_bc(bc, grid, fold)
      type(ocean_bc_state_t), intent(out) :: bc
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: fold
      call ocean_bc_state_init(bc, grid, nz_ml=0, n_tracers=0)
      bc%west%bc_type = OBC_PERIODIC
      bc%east%bc_type = OBC_PERIODIC
      call ocean_bc_validate_periodic(bc)
      if (fold) then
         bc%north%bc_type = OBC_TRIPOLAR_FOLD
         call ocean_bc_validate_fold(bc)
      end if
      ! south stays the `ocean_bc_state_init` default (OBC_WALL).
   end subroutine make_bc

   pure function default_params() result(par)
      type(ice_evp_params_t) :: par
      par%p0 = 2.75e4_wp
      par%c0 = 20.0_wp
      par%ec = 2.0_wp
      par%cdw = 3.24e-3_wp
      par%rho_ocean = 1030.0_wp
      par%del_sh_min_scale = 2.0_wp
      par%tdamp = -0.2_wp
      par%evp_sub_steps = 60
      par%a_face_stress = .false.
      par%cfl_trunc = 0.5_wp
         !! SIS2 default.
      par%cfl_trunc_dyn_its = .false.
      par%project_ci = .false.
   end function default_params

   subroutine alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                               str_d, str_t, str_s, fxoc, fyoc, f_corner)
      integer, intent(in) :: nx, ny
      real(wp), allocatable, intent(out) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), allocatable, intent(out) :: uo(:, :), vo(:, :)
      real(wp), allocatable, intent(out) :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable, intent(out) :: ui(:, :), vi(:, :)
      real(wp), allocatable, intent(out) :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable, intent(out) :: fxoc(:, :), fyoc(:, :)
      real(wp), allocatable, intent(out) :: f_corner(:, :)

      allocate (mis(nx, ny), mice(nx, ny), ci(nx, ny))
      allocate (uo(nx + 1, ny), vo(nx, ny + 1))
      allocate (tau_ax(nx + 1, ny), tau_ay(nx, ny + 1))
      allocate (ui(nx + 1, ny), vi(nx, ny + 1))
      allocate (str_d(nx, ny), str_t(nx, ny), str_s(nx + 1, ny + 1))
      allocate (fxoc(nx + 1, ny), fyoc(nx, ny + 1))
      allocate (f_corner(nx + 1, ny + 1))
   end subroutine alloc_evp_arrays

   subroutine seed_uniform_state(metrics, mis, mice, ci, uo, vo, tau_ax, tau_ay, &
                                 ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      !! Uniform ice over the whole (all-wet) domain, ocean at rest, a
      !! uniform zonal wind, EVP state seeded at rest -- see module header.
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(inout) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), intent(inout) :: uo(:, :), vo(:, :)
      real(wp), intent(inout) :: tau_ax(:, :), tau_ay(:, :)
      real(wp), intent(inout) :: ui(:, :), vi(:, :)
      real(wp), intent(inout) :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), intent(inout) :: fxoc(:, :), fyoc(:, :)
      real(wp), intent(inout) :: f_corner(:, :)
      integer :: i, j

      mis = MIS_CONST; mice = MICE_CONST; ci = CI_CONST
      uo = 0.0_wp; vo = 0.0_wp
      tau_ax = TAU_AX_CONST; tau_ay = 0.0_wp
      ui = 0.0_wp; vi = 0.0_wp
      str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
      fxoc = 0.0_wp; fyoc = 0.0_wp

      ! Geometric Coriolis parameter -- `geolatBu` is already fold-correct
      ! (the bipolar cap's own reflection symmetry preserves latitude),
      ! same construction `test_ocean_tripolar` uses; no explicit fold
      ! call needed.
      do j = 1, size(f_corner, 2)
         do i = 1, size(f_corner, 1)
            f_corner(i, j) = 2.0_wp*OMEGA_EARTH*sin(metrics%geolatBu(i, j)*PI/180.0_wp)
         end do
      end do
   end subroutine seed_uniform_state

   ! -----------------------------------------------------------------
   ! Case 1: halo == exact fold image after one outer step.
   ! -----------------------------------------------------------------

   subroutine test_halo_mirrors(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j, im, jm, ng, ni, nj
      real(wp) :: max_abs_sd, worst_corner

      call make_grid_metrics(grid, metrics)
      call make_bc(bc, grid, fold=.true.)
      nx = grid%nx_total; ny = grid%ny_total
      ng = NGHOST; ni = NXP; nj = NYP
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call seed_uniform_state(metrics, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                              str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%init(nx, ny)
      call ws%enter_data()

      !$omp target enter data map(to: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target enter data map(to: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

      call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                            tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                            fxoc, fyoc, DT_SLOW, par, .true., .false., ws, &
                            land_w=.false., land_e=.false., land_s=.true., &
                            land_n=.false., bc=bc)

      !$omp target update from(ui, vi, str_d, str_t, str_s)
      !$omp target exit data map(delete: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target exit data map(delete: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%exit_data()
      call ws%destroy()

      ! ui: halo rows j >= ng+nj+1, fsum=2ng+ni+2, jsum=2ng+2nj+1.
      do j = ng + nj + 1, ny
         do i = 1, nx + 1
            im = 2*ng + ni + 2 - i
            jm = 2*ng + 2*nj + 1 - j
            call check(error, ui(i, j), -ui(im, jm), "ui fold-halo mismatch")
            if (allocated(error)) return
         end do
      end do

      ! vi: halo rows j > ng+nj+1, isum=2ng+ni+1, jsum=2ng+2nj+2.
      do j = ng + nj + 2, ny + 1
         do i = 1, nx
            im = 2*ng + ni + 1 - i
            jm = 2*ng + 2*nj + 2 - j
            call check(error, vi(i, j), -vi(im, jm), "vi fold-halo mismatch")
            if (allocated(error)) return
         end do
      end do

      ! str_d / str_t halo rows equal the interior centre-fold image to
      ! ROUND-OFF, not bitwise: the subcycle loop recomputes them in the
      ! halo from folded `ui`/`vi`, summing the mirrored operands in the
      ! opposite order (they are folded bitwise only once, before the
      ! loop, and again at the next call). The outermost halo row is
      ! excluded: its stencil runs off the array and no physical cell
      ! reads it.
      max_abs_sd = max(maxval(abs(str_d)), maxval(abs(str_t)), 1.0_wp)
      do j = ng + nj + 1, ny - 1
         do i = 1, nx
            im = 2*ng + ni + 1 - i
            jm = 2*ng + 2*nj + 1 - j
            call check(error, abs(str_d(i, j) - str_d(im, jm)) <= 1.0e-13_wp*max_abs_sd, &
                       "str_d fold-halo mismatch")
            if (allocated(error)) return
            call check(error, abs(str_t(i, j) - str_t(im, jm)) <= 1.0e-13_wp*max_abs_sd, &
                       "str_t fold-halo mismatch")
            if (allocated(error)) return
         end do
      end do

      ! str_s / sh_ds on the fold row stay parity-symmetric to round-off
      ! (NOT bitwise -- the subcycle loop evolves str_s from folded
      ! `ui`/`vi` afresh every step, it is not re-folded itself).  str_s's
      ! fold uses `negate=.false.` (parity +1, a pseudoscalar like
      ! vorticity) so the projection COPIES, not negates: the correct
      ! relation is EQUALITY, a(i) == a(mirror), not antisymmetry.
      max_abs_sd = maxval(abs(str_s))
      worst_corner = 0.0_wp
      do i = ng + 1, ng + ni
         im = 2*ng + ni + 2 - i  ! corners mirror like u faces (Bu = SW corner of T)
         worst_corner = max(worst_corner, abs(str_s(i, ng + nj + 1) - str_s(im, ng + nj + 1)))
      end do
      call check(error, worst_corner <= 1.0e-13_wp*max(max_abs_sd, 1.0_wp), &
                 "str_s fold-row parity off by more than round-off")
      if (allocated(error)) return

      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
   end subroutine test_halo_mirrors

   ! -----------------------------------------------------------------
   ! Case 2: on-line vi projection antisymmetric + mass flux telescopes.
   ! -----------------------------------------------------------------

   subroutine test_seam_antisym(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, im, ng, ni, nj, j_fold
      real(wp) :: flux_sum, flux_scale

      call make_grid_metrics(grid, metrics)
      call make_bc(bc, grid, fold=.true.)
      nx = grid%nx_total; ny = grid%ny_total
      ng = NGHOST; ni = NXP; nj = NYP
      j_fold = ng + nj + 1
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call seed_uniform_state(metrics, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                              str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%init(nx, ny)
      call ws%enter_data()

      !$omp target enter data map(to: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target enter data map(to: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

      call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                            tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                            fxoc, fyoc, DT_SLOW, par, .true., .false., ws, &
                            land_w=.false., land_e=.false., land_s=.true., &
                            land_n=.false., bc=bc)

      !$omp target update from(vi)
      !$omp target exit data map(delete: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target exit data map(delete: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%exit_data()
      call ws%destroy()

      ! On-line projection: isum = 2*ng+ni+1.
      flux_sum = 0.0_wp
      flux_scale = 0.0_wp
      do i = ng + 1, ng + ni
         im = 2*ng + ni + 1 - i
         call check(error, vi(i, j_fold) + vi(im, j_fold), 0.0_wp, &
                    "vi fold-row projection not antisymmetric")
         if (allocated(error)) return
         ! Uniform mi_v (== MIS_CONST, all-wet) * true-vector vi *
         ! fold-symmetric dx_cv: the ice analogue of the ocean B1 leak.
         flux_sum = flux_sum + MIS_CONST*vi(i, j_fold)*metrics%dx_cv(i, j_fold)
         flux_scale = flux_scale + MIS_CONST*abs(vi(i, j_fold))*metrics%dx_cv(i, j_fold)
      end do
      call check(error, abs(flux_sum) <= 1.0e-12_wp*max(flux_scale, 1.0_wp), &
                 "cross-seam ice mass flux does not telescope to round-off")
      if (allocated(error)) return

      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
   end subroutine test_seam_antisym

   ! -----------------------------------------------------------------
   ! Case 3: many-step stability + no seam stripe + wall<->fold proof.
   ! -----------------------------------------------------------------

   subroutine run_n_steps(grid, metrics, bc, par, n_outer, land_n, use_bc, &
                          ui, vi, n_trunc_total)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_bc_state_t), intent(in) :: bc
      type(ice_evp_params_t), intent(in) :: par
      integer, intent(in) :: n_outer
      logical, intent(in) :: land_n, use_bc
      real(wp), allocatable, intent(out) :: ui(:, :), vi(:, :)
      integer, intent(out) :: n_trunc_total

      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, n_out

      nx = grid%nx_total; ny = grid%ny_total
      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call seed_uniform_state(metrics, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                              str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%init(nx, ny)
      call ws%enter_data()

      !$omp target enter data map(to: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target enter data map(to: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

      n_trunc_total = 0
      do n = 1, n_outer
         if (use_bc) then
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .false., ws, &
                                  n_trunc=n_out, &
                                  land_w=.false., land_e=.false., land_s=.true., &
                                  land_n=land_n, bc=bc)
         else
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .false., ws, &
                                  n_trunc=n_out, &
                                  land_w=.false., land_e=.false., land_s=.true., &
                                  land_n=land_n)
         end if
         n_trunc_total = n_trunc_total + n_out
      end do

      !$omp target update from(ui, vi)
      !$omp target exit data map(delete: mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$omp target exit data map(delete: ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%exit_data()
      call ws%destroy()
   end subroutine run_n_steps

   subroutine test_no_seam_stripe(error)
      !! Deviation from the plan: no `ice_transport_step` leg -- see the
      !! module header.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc_fold, bc_wall
      type(ice_evp_params_t) :: par
      real(wp), allocatable :: ui(:, :), vi(:, :), ui_wall(:, :), vi_wall(:, :)
      integer :: n_trunc, n_trunc_wall
      integer, parameter :: N_OUTER = 48
      integer :: ng, ni, nj, i, j_last, j_prev
      real(wp) :: max_u, jump_boundary, max_v_fold, max_v_wall
      real(wp), allocatable :: jumps(:)

      call make_grid_metrics(grid, metrics)
      call make_bc(bc_fold, grid, fold=.true.)
      par = default_params()
      ng = NGHOST; ni = NXP; nj = NYP

      call run_n_steps(grid, metrics, bc_fold, par, N_OUTER, land_n=.false., &
                       use_bc=.true., ui=ui, vi=vi, n_trunc_total=n_trunc)

      call check(error, all(ieee_is_finite(ui)), "ui not finite after 48 steps")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(vi)), "vi not finite after 48 steps")
      if (allocated(error)) return

      max_u = maxval(abs(ui(ng + 1:ng + ni + 1, ng + 1:ng + nj)))
      call check(error, max_u < 1.0_wp, "max|u_ice| >= 1 m/s after 48 steps")
      if (allocated(error)) return

      ! The coarse test grid's tiny near-pole cells trip the CFL clip
      ! whatever the north boundary is; the fold must add none of its own.
      call make_bc(bc_wall, grid, fold=.false.)
      call run_n_steps(grid, metrics, bc_wall, par, N_OUTER, land_n=.true., &
                       use_bc=.false., ui=ui_wall, vi=vi_wall, n_trunc_total=n_trunc_wall)
      call check(error, n_trunc <= n_trunc_wall, "fold adds CFL truncation over the wall run")
      if (allocated(error)) return

      ! Stripe detector on ui: the inter-row jump at the fold boundary row
      ! must not be an outlier against the typical interior inter-row
      ! jump.
      j_last = ng + nj
      j_prev = j_last - 1
      allocate (jumps(j_prev - (ng + 1) + 1))
      do j_prev = ng + 1, j_last - 1
         jumps(j_prev - ng) = maxval(abs(ui(ng + 1:ng + ni + 1, j_prev + 1) - &
                                         ui(ng + 1:ng + ni + 1, j_prev)))
      end do
      jump_boundary = maxval(abs(ui(ng + 1:ng + ni + 1, j_last) - &
                                 ui(ng + 1:ng + ni + 1, j_last - 1)))
      call check(error, jump_boundary <= 3.0_wp*median(jumps), &
                 "seam-row stripe: boundary inter-row jump is an outlier")
      if (allocated(error)) return

      ! Live-change proof: an identical setup with a WALL instead of the
      ! fold must give exactly zero v on the boundary face (walled); the
      ! fold run must give nonzero v there (it actually flows).
      max_v_fold = maxval(abs(vi(ng + 1:ng + ni, ng + nj + 1)))
      max_v_wall = maxval(abs(vi_wall(ng + 1:ng + ni, ng + nj + 1)))
      call check(error, max_v_wall, 0.0_wp, "wall run: boundary v is not exactly zero")
      if (allocated(error)) return
      call check(error, max_v_fold > 0.0_wp, "fold run: boundary v is zero (fold not live)")
      if (allocated(error)) return

      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc_fold)
      call ocean_bc_state_destroy(bc_wall)
   end subroutine test_no_seam_stripe

   pure function median(x) result(m)
      real(wp), intent(in) :: x(:)
      real(wp) :: m
      real(wp) :: sorted(size(x))
      integer :: n
      sorted = x
      call sort_ascending(sorted)
      n = size(sorted)
      if (mod(n, 2) == 1) then
         m = sorted((n + 1)/2)
      else
         m = 0.5_wp*(sorted(n/2) + sorted(n/2 + 1))
      end if
   end function median

   pure subroutine sort_ascending(x)
      !! Insertion sort -- `x` is tiny (O(ny_phys)).
      real(wp), intent(inout) :: x(:)
      integer :: i, j
      real(wp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_ascending

end module test_ocean_ice_evp_fold
