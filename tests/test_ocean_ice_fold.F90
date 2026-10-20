!! Tripolar north-fold regression tests for the sea-ice category state
!! (fold-seam fix, `fix/ice-fold-seam`).
!!
!! **Root cause this guards.**  `ocean_halo_exchange_ice_state` /
!! `ocean_halo_exchange_ice_transport` refreshed the ice category state's
!! seam ghosts through the plain MPI/periodic halo primitives only -- on a
!! `north = 'tripolar_fold'` grid those primitives fill every ghost band
!! EXCEPT the north cap, which needs the 180-degree-rotated mirror
!! (`rdb_ocean_fold`).  The ocean's own prognostics (`h_layer`,
!! `u_face_x_layer`, ...) got that fold through `ocean_fold_wrap_state`;
!! the ice category state (`part_size`/`m_ice`/`m_snow`/`enth_ice`/
!! `sal_ice`/`enth_snow`/`mca_ice`/`mca_snow`) did not.  Every reader one
!! cell into that band (the transport's 5-point PPM stencil, first and
!! foremost) therefore saw a stale/zero mirror instead of the correct
!! fold image, which is how ice piled up without bound on the fold-seam
!! row (max 1.4 m -> 195 m over a year on the global 1-degree ice case).
!!
!! Cases:
!!   * `ice_state_fold_mirrors_exactly` / `ice_transport_fold_mirrors_exactly`
!!     -- seed every category field with a j/i-dependent (NOT already
!!     fold-symmetric) pattern on physical cells only, call the exchange
!!     host-side (`device_resident=.false.`), and assert every north-halo
!!     ghost cell equals the EXACT `fold_north_centre` mirror of its
!!     physical source cell, for every field the two exchanges carry.
!!     FAILS before the fix (the ghost band stays at its `init` value --
!!     0 for everything but `part_size(:,:,0)`, which stays 1 -- while the
!!     true mirror is the seeded nonzero pattern).
!!   * `ice_fold_seam_conserves_and_bounded` -- a full `ice_transport_step`
!!     integration: an i-mirror-SYMMETRIC category IC (uniform across
!!     physical i) advected by a meridional velocity at the exact fold
!!     face that is i-ANTIsymmetric (step function, `v(p) = -v(mirror
!!     p)`, the only velocity pattern consistent with the fold's own
!!     v-face contract) over N steps.  The combination is a symmetry of
!!     the transport equations, so a correctly fold-exchanged run must
!!     keep the category state exactly i-mirror-symmetric and conserve
!!     total ice mass to round-off; the fold-row thickness, bounded by
!!     that conserved mass, cannot blow up.  Before the fix the PPM
!!     stencil at the seam face reads the wrong ghost, mass conservation
!!     and the mirror symmetry both break, and the fold-row thickness
!!     grows well past the bound this test enforces.
module test_ocean_ice_fold
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_tripolar_metrics, destroy_cartesian_metrics
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_transport, only: ice_transport_step
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ice_state, &
                                   ocean_halo_exchange_ice_transport
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_PERIODIC, &
                                       OBC_TRIPOLAR_FOLD, ocean_bc_validate_periodic, &
                                       ocean_bc_validate_fold
   implicit none
   private

   public :: collect_ocean_ice_fold_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NZ = 2
   integer, parameter :: NXP = 16, NYP = 10
      !! NXP even (no self-conjugate v column); NYP >= NGHOST+1.
   integer, parameter :: NCAT = 2
   integer, parameter :: NK_ICE = 2
   real(wp), parameter :: DLON = 360.0_wp/real(NXP, wp)
   real(wp), parameter :: DLAT = 3.0_wp
   real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 50.0_wp
   real(wp), parameter :: REARTH = 6.378e6_wp
   real(wp), parameter :: PHI_JOIN = 60.0_wp, LON_POLE = 0.0_wp
   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_ice_fold_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ice_state_fold_mirrors_exactly", test_state_fold_mirrors), &
                  new_unittest("ice_transport_fold_mirrors_exactly", &
                               test_transport_fold_mirrors), &
                  new_unittest("ice_fold_seam_conserves_and_bounded", &
                               test_seam_conserves_and_bounded) &
                  ]
   end subroutine collect_ocean_ice_fold_tests

   ! -----------------------------------------------------------------
   ! Shared setup
   ! -----------------------------------------------------------------

   subroutine make_grid_bc(grid, bc)
      type(hgrid_t), intent(out) :: grid
      type(ocean_bc_state_t), intent(out) :: bc
      call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
      call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=0)
      bc%west%bc_type = OBC_PERIODIC
      bc%east%bc_type = OBC_PERIODIC
      bc%north%bc_type = OBC_TRIPOLAR_FOLD
      call ocean_bc_validate_periodic(bc)
      call ocean_bc_validate_fold(bc)
         !! `has_north` defaults `.true.` (single-rank), so this alone
         !! sets `bc%north_fold = .true.` -- matches `test_ocean_tripolar`
         !! `make_bc_fold`, no `ocean_bc_state_set_edges` call needed.
   end subroutine make_grid_bc

   subroutine make_ice(grid, ice)
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(inout) :: ice
      ice%enable = .true.
      ice%ncat = NCAT
      ice%nk_ice = NK_ICE
      ice%transport = .true.
      call ice%init(grid)
   end subroutine make_ice

   !! Fold-image indices for a T-stagger (centre) field -- the EXACT
   !! formula `fold_north_centre` applies (see `rdb_ocean_fold`'s header
   !! table): i' = 2*nghost+nx_phys+1-i, j' = 2*nghost+2*ny_phys+1-j.
   pure subroutine t_mirror(i, j, grid, im, jm)
      integer, intent(in) :: i, j
      type(hgrid_t), intent(in) :: grid
      integer, intent(out) :: im, jm
      im = 2*grid%nghost + grid%nx_phys + 1 - i
      jm = 2*grid%nghost + 2*grid%ny_phys + 1 - j
   end subroutine t_mirror

   ! -----------------------------------------------------------------
   ! Case 1 + 2: ghost == exact mirror after the exchange.
   ! -----------------------------------------------------------------

   subroutine test_state_fold_mirrors(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      type(ocean_sea_ice_t) :: ice
      integer :: i, j, im, jm, c, l

      call make_grid_bc(grid, bc)
      call make_ice(grid, ice)

      ! Seed a j/i-varying, non-fold-symmetric pattern on PHYSICAL cells
      ! only (ghosts stay at `init`'s zero / part_size(:,:,0)=1).
      do j = grid%nghost + 1, grid%nghost + grid%ny_phys
         do i = grid%nghost + 1, grid%nghost + grid%nx_phys
            do c = 1, ice%ncat
               ice%part_size(i, j, c) = 0.1_wp*real(c, wp) + 0.01_wp*real(i + 7*j, wp)
               ice%m_ice(i, j, c) = 100.0_wp*real(c, wp) + real(i + 3*j, wp)
               ice%m_snow(i, j, c) = 10.0_wp*real(c, wp) + real(2*i + j, wp)
               do l = 1, ice%nk_ice
                  ice%enth_ice(i, j, c, l) = -1.0e5_wp + real(i + j + l, wp)
                  ice%sal_ice(i, j, c, l) = 4.0_wp + 0.01_wp*real(i + l, wp)
               end do
               ice%enth_snow(i, j, c, 1) = -5.0e4_wp + real(i - j, wp)
            end do
            ice%part_size(i, j, 0) = 1.0_wp - sum(ice%part_size(i, j, 1:ice%ncat))
         end do
      end do

      call ocean_halo_exchange_ice_state(ice, grid, bc, device_resident=.false.)

      ! Every north-halo T-row must now equal the exact centre-fold
      ! mirror of its (physical) source cell, for every field the
      ! exchange carries.
      do j = grid%nghost + grid%ny_phys + 1, grid%ny_total
         do i = 1, grid%nx_total
            call t_mirror(i, j, grid, im, jm)
            do c = 1, ice%ncat
               call check(error, ice%part_size(i, j, c), ice%part_size(im, jm, c), &
                          "part_size fold mismatch")
               if (allocated(error)) return
               call check(error, ice%m_ice(i, j, c), ice%m_ice(im, jm, c), &
                          "m_ice fold mismatch")
               if (allocated(error)) return
               call check(error, ice%m_snow(i, j, c), ice%m_snow(im, jm, c), &
                          "m_snow fold mismatch")
               if (allocated(error)) return
               do l = 1, ice%nk_ice
                  call check(error, ice%enth_ice(i, j, c, l), ice%enth_ice(im, jm, c, l), &
                             "enth_ice fold mismatch")
                  if (allocated(error)) return
                  call check(error, ice%sal_ice(i, j, c, l), ice%sal_ice(im, jm, c, l), &
                             "sal_ice fold mismatch")
                  if (allocated(error)) return
               end do
               call check(error, ice%enth_snow(i, j, c, 1), ice%enth_snow(im, jm, c, 1), &
                          "enth_snow fold mismatch")
               if (allocated(error)) return
            end do
            call check(error, ice%part_size(i, j, 0), ice%part_size(im, jm, 0), &
                       "part_size(0) fold mismatch")
            if (allocated(error)) return
         end do
      end do

      call ice%destroy()
      call ocean_bc_state_destroy(bc)
   end subroutine test_state_fold_mirrors

   subroutine test_transport_fold_mirrors(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      type(ocean_sea_ice_t) :: ice
      integer :: i, j, im, jm, c, l

      call make_grid_bc(grid, bc)
      call make_ice(grid, ice)

      do j = grid%nghost + 1, grid%nghost + grid%ny_phys
         do i = grid%nghost + 1, grid%nghost + grid%nx_phys
            do c = 1, ice%ncat
               ice%mca_ice(i, j, c) = 90.0_wp*real(c, wp) + real(i + 5*j, wp)
               ice%mca_snow(i, j, c) = 9.0_wp*real(c, wp) + real(i + j, wp)
               ice%m_ice(i, j, c) = 100.0_wp*real(c, wp) + real(i + 3*j, wp)
               do l = 1, ice%nk_ice
                  ice%enth_ice(i, j, c, l) = -1.0e5_wp + real(i + j + l, wp)
                  ice%sal_ice(i, j, c, l) = 4.0_wp + 0.01_wp*real(i + l, wp)
               end do
               ice%enth_snow(i, j, c, 1) = -5.0e4_wp + real(i - j, wp)
            end do
         end do
      end do

      call ocean_halo_exchange_ice_transport(ice, grid, bc, device_resident=.false.)

      do j = grid%nghost + grid%ny_phys + 1, grid%ny_total
         do i = 1, grid%nx_total
            call t_mirror(i, j, grid, im, jm)
            do c = 1, ice%ncat
               call check(error, ice%mca_ice(i, j, c), ice%mca_ice(im, jm, c), &
                          "mca_ice fold mismatch")
               if (allocated(error)) return
               call check(error, ice%mca_snow(i, j, c), ice%mca_snow(im, jm, c), &
                          "mca_snow fold mismatch")
               if (allocated(error)) return
               call check(error, ice%m_ice(i, j, c), ice%m_ice(im, jm, c), &
                          "m_ice fold mismatch (transport exchange)")
               if (allocated(error)) return
            end do
         end do
      end do

      call ice%destroy()
      call ocean_bc_state_destroy(bc)
   end subroutine test_transport_fold_mirrors

   ! -----------------------------------------------------------------
   ! Case 3: full transport integration across the fold.
   ! -----------------------------------------------------------------

   subroutine test_seam_conserves_and_bounded(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(decomp_t) :: decomp
      integer, parameter :: NSTEPS = 20
      real(wp), parameter :: DT = 1800.0_wp
      real(wp), parameter :: V0 = 0.05_wp
         !! m/s at the exact fold face -- generous but CFL-safe against a
         !! ~3-degree (~330 km) row width over DT=1800s (~90 m/step).
      real(wp), parameter :: MASS_REL_TOL = 1.0e-10_wp
      real(wp), parameter :: SYMM_TOL = 1.0e-10_wp
      real(wp) :: mass0, mass1, h_max0, h_max1
      integer :: i, j, c, n, i_lo, i_hi, j_lo, j_hi, ng, p, pm, fold_j
      logical :: ok

      call make_grid_bc(grid, bc)
      call make_ice(grid, ice)
      call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, REARTH, &
                                 PHI_JOIN, LON_POLE)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 500.0_wp

      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
      call decomp_init(decomp, grid%nx_phys, grid%ny_phys, 1, 1, 0)
      call ocean_halo_init(decomp, grid%nghost, bc%periodic_x, bc%periodic_y)

      ng = grid%nghost
      i_lo = ng + 1
      i_hi = ng + grid%nx_phys
      j_lo = ng + 1
      j_hi = ng + grid%ny_phys
      fold_j = ng + grid%ny_phys + 1
         !! The exact fold-line v-face row (south face of the first ghost
         !! row == north face of the last physical T-row).

      ! i-mirror-SYMMETRIC category IC: uniform across physical i (any
      ! uniform-in-i field is trivially its own mirror).
      do j = j_lo, j_hi
         do i = i_lo, i_hi
            do c = 1, ice%ncat
               ice%part_size(i, j, c) = 0.3_wp
               ice%m_ice(i, j, c) = 300.0_wp*real(c, wp)
               ice%enth_ice(i, j, c, :) = -2.0e5_wp
               ice%sal_ice(i, j, c, :) = 4.0_wp
            end do
            ice%part_size(i, j, 0) = 1.0_wp - sum(ice%part_size(i, j, 1:ice%ncat))
         end do
      end do

      ! i-ANTIsymmetric velocity at the exact fold face -- the only
      ! pattern consistent with the fold's own v-face contract
      ! (`v(p, fold) = -v(mirror p, fold)`): a step function, +V0 on the
      ! west half of the cap and -V0 on its mirror (east) half.  Zero
      ! everywhere else (u, and v off the fold face) keeps the problem
      ! purely 1-D across the seam.
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do i = i_lo, i_hi
         p = i - ng
         pm = grid%nx_phys + 1 - p
         if (p < pm) then
            ms%v_face_y_layer(i, fold_j, NZ) = V0
            ms%v_face_y_layer(ng + pm, fold_j, NZ) = -V0
         end if
      end do

      mass0 = global_ice_mass(grid, ice, i_lo, i_hi, j_lo, j_hi)
      h_max0 = max_fold_row_thickness(grid, ice, i_lo, i_hi, j_hi)

      ! GPU mem:separate discipline (CLAUDE.md "Writing GPU tests"): every
      ! array `ice_transport_step` touches must be device-present BEFORE
      ! the first call, or its `do concurrent`/`!$acc` kernels silently
      ! read/write stale host memory -- no crash, just a wrong answer
      ! (this is exactly what happened here before this map was added:
      ! host-only, the case passed; on the GPU build it failed on mass
      ! conservation, not because the fix is wrong, but because the test
      ! never mapped `ms`/`ice`). Map once, after the IC/velocity seed
      ! above (setup-time host edits must precede the map), pull the
      ! touched arrays back every step for the per-step symmetry check,
      ! unmap once at the end.
      !$omp target enter data map(to: ms)
      call ms%enter_data()
      call ice%enter_data()

      do n = 1, NSTEPS
         call ice_transport_step(grid, metrics, ms, ice, DT, 1, 0.0_wp, ok, bc)
         call check(error, ok, "ice_transport_step reported not-ok")
         if (.not. allocated(error)) then
            associate (ps => ice%part_size, mi => ice%m_ice)
               !$omp target update from(ps, mi)
            end associate
         end if
         if (allocated(error)) exit

         ! Mirror symmetry: c(p, j) == c(mirror p, j) for EVERY physical
         ! cell, every step -- a genuine symmetry of this IC + velocity
         ! pair, and exactly what a correctly fold-exchanged run preserves.
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               p = i - ng
               pm = grid%nx_phys + 1 - p
               do c = 1, ice%ncat
                  call check(error, ice%m_ice(i, j, c), ice%m_ice(ng + pm, j, c), &
                             "m_ice lost i-mirror symmetry crossing the fold", &
                             rel=.true., thr=SYMM_TOL)
                  if (allocated(error)) exit
               end do
               if (allocated(error)) exit
            end do
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      if (allocated(error)) then
         call ice%exit_data(); call ms%exit_data()
         !$omp target exit data map(delete: ms)
         call ice%destroy(); call ms%destroy(); call ocean_bc_state_destroy(bc)
         call destroy_cartesian_metrics(metrics); call ocean_halo_destroy()
         return
      end if

      call ice%exit_data()
      call ms%exit_data()
      !$omp target exit data map(delete: ms)

      mass1 = global_ice_mass(grid, ice, i_lo, i_hi, j_lo, j_hi)
      h_max1 = max_fold_row_thickness(grid, ice, i_lo, i_hi, j_hi)

      call check(error, abs(mass1 - mass0) <= MASS_REL_TOL*max(abs(mass0), 1.0_wp), &
                 "total ice mass not conserved across the fold seam")
      if (.not. allocated(error)) then
         ! The bug this guards grew the fold-row max 1.4 m -> 195 m over a
         ! year (140x) on the production case; 20 CFL-modest steps with a
         ! conserved, finite mass supply must stay WELL inside a
         ! generous 5x bound -- a real blow-up trips this hard.
         call check(error, h_max1 <= 5.0_wp*max(h_max0, 1.0e-6_wp), &
                    "fold-row ice thickness grew unbounded")
      end if

      call ice%destroy()
      call ms%destroy()
      call ocean_bc_state_destroy(bc)
      call destroy_cartesian_metrics(metrics)
      call ocean_halo_destroy()
   end subroutine test_seam_conserves_and_bounded

   function global_ice_mass(grid, ice, i_lo, i_hi, j_lo, j_hi) result(m)
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(in) :: ice
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi
      real(wp) :: m
      integer :: i, j, c
      m = 0.0_wp
      do j = j_lo, j_hi
         do i = i_lo, i_hi
            do c = 1, ice%ncat
               m = m + ice%part_size(i, j, c)*ice%m_ice(i, j, c)
            end do
         end do
      end do
   end function global_ice_mass

   function max_fold_row_thickness(grid, ice, i_lo, i_hi, j_hi) result(h)
      !! Category-summed ice "thickness" proxy (mass / rho) on the last
      !! PHYSICAL row (`j_hi`), the fold-seam row the production bug
      !! piled up on.  `ICE_RHO_ICE` is not imported here -- the test
      !! only needs a monotone proxy for "did this blow up", so a plain
      !! mass/area number (kg/m^2) serves identically; the 5x bound in
      !! the caller is scale-invariant.
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(in) :: ice
      integer, intent(in) :: i_lo, i_hi, j_hi
      real(wp) :: h
      integer :: i, c
      real(wp) :: cell
      h = 0.0_wp
      do i = i_lo, i_hi
         cell = 0.0_wp
         do c = 1, ice%ncat
            cell = cell + ice%part_size(i, j_hi, c)*ice%m_ice(i, j_hi, c)
         end do
         h = max(h, cell)
      end do
   end function max_fold_row_thickness

end module test_ocean_ice_fold
