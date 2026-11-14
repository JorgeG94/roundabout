!! C-grid ice metrics map, corner-wrap helper and face-pair halo wrapper
!! (CICE6 port D9.0).  No physics.
!!
!! Cases:
!!   * `cgrid_index_map_*` -- every ice name equals its ocean source at the
!!     documented offset (E = Cu(i+1,j), N = Cv(i,j+1), U = Bu(i+1,j+1)),
!!     bitwise, over the WHOLE array, on anisotropic-Cartesian, lon-lat and
!!     tripolar metrics; masks follow `wet_u/v/q`.
!!   * `cgrid_ratios` -- `ratiodxN(i,j)*dxN(i,j) == -dxN(i+1,j)`,
!!     `ratiodyE*dyE == -dyE(i,j+1)`, reciprocals multiply to 1 where
!!     defined, last column/row and zero denominators are 0.
!!   * `cgrid_masks_land` -- an interior island closes the faces/corners
!!     touching it through `epm/npm/uvm` (hand-derived from the wet cells).
!!   * `cgrid_fold_ghost_metrics` -- on the tripolar cap the north-ghost
!!     metrics of the ice names are finite, non-negative, and are the exact
!!     fold image of the interior (`dxN(:,nj+2)`-type rows, `earea` ghosts).
!!   * `corner_wrap_periodic` -- the promoted `evp_wrap_corner_impl` makes
!!     every ghost corner the periodic image and the east seam corner equal
!!     the west one.
!!   * `halo_faces_equals_ocean_*` -- `ice_halo_faces` is BITWISE the
!!     ocean's own periodic wrap + fold kernels on (periodic, wall,
!!     tripolar) edges; the fold ghosts are the negated mirror image and the
!!     on-line v row is antisymmetric.
!!   * `cgrid_device_kernel` -- the mapped ratio/mask arrays are read by a
!!     device kernel (GPU `mem:separate` residency of `enter_data`).
module test_cice_cgrid_metrics
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, make_anisotropic_metrics, &
                                 make_spherical_metrics, make_tripolar_metrics, &
                                 destroy_cartesian_metrics
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, OBC_OPEN, &
                                       OBC_PERIODIC, OBC_TRIPOLAR_FOLD, &
                                       ocean_bc_validate_periodic, ocean_bc_validate_fold
   use rdb_ocean_periodic, only: ocean_periodic_wrap_face_x_2d, ocean_periodic_wrap_face_y_2d
   use rdb_ocean_fold, only: fold_north_u_face, fold_north_v_face
   use rdb_ice_evp, only: evp_wrap_corner_impl
   use rdb_ice_cgrid_metrics, only: ice_cgrid_metrics_t, ice_halo_faces
   implicit none
   private

   public :: collect_cice_cgrid_metrics_tests

   integer, parameter :: NG = 3
   integer, parameter :: NXP = 12, NYP = 10
      !! NXP even (no self-conjugate v column under the fold).
   real(wp), parameter :: DLON = 360.0_wp/real(NXP, wp)
   real(wp), parameter :: DLAT = 3.0_wp
   real(wp), parameter :: REARTH = 6.378e6_wp

contains

   subroutine collect_cice_cgrid_metrics_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cgrid_index_map_cartesian", test_map_cartesian), &
                  new_unittest("cgrid_index_map_spherical", test_map_spherical), &
                  new_unittest("cgrid_index_map_tripolar", test_map_tripolar), &
                  new_unittest("cgrid_ratios", test_ratios), &
                  new_unittest("cgrid_masks_land", test_masks_land), &
                  new_unittest("cgrid_masks_edges", test_masks_edges), &
                  new_unittest("cgrid_stretched_lengths", test_stretched_lengths), &
                  new_unittest("cgrid_fold_ghost_metrics", test_fold_ghost_metrics), &
                  new_unittest("corner_wrap_periodic", test_corner_wrap), &
                  new_unittest("halo_faces_equals_ocean_periodic", test_halo_periodic), &
                  new_unittest("halo_faces_equals_ocean_wall", test_halo_wall), &
                  new_unittest("halo_faces_equals_ocean_fold", test_halo_fold), &
                  new_unittest("cgrid_device_kernel", test_device_kernel) &
                  ]
   end subroutine collect_cice_cgrid_metrics_tests

   ! -----------------------------------------------------------------
   ! helpers
   ! -----------------------------------------------------------------

   subroutine make_bc(bc, grid, periodic_x, fold)
      type(ocean_bc_state_t), intent(out) :: bc
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: periodic_x, fold
      call ocean_bc_state_init(bc, grid, nz_ml=0, n_tracers=0)
      if (periodic_x) then
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
      end if
      if (fold) then
         bc%north%bc_type = OBC_TRIPOLAR_FOLD
         call ocean_bc_validate_fold(bc)
      end if
   end subroutine make_bc

   subroutine check_map(error, cm, m, nx, ny)
      !! Whole-array bitwise index-map check.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cgrid_metrics_t), intent(in) :: cm
      type(ocean_metrics_t), intent(in) :: m
      integer, intent(in) :: nx, ny
      integer :: i, j
      logical :: bad
      bad = .false.
      do j = 1, ny
         do i = 1, nx
            if (cm%dxE(i, j) /= m%dxCu(i + 1, j)) bad = .true.
            if (cm%dyE(i, j) /= m%dyCu(i + 1, j)) bad = .true.
            if (cm%dxN(i, j) /= m%dxCv(i, j + 1)) bad = .true.
            if (cm%dyN(i, j) /= m%dyCv(i, j + 1)) bad = .true.
            if (cm%dxU(i, j) /= m%dxBu(i + 1, j + 1)) bad = .true.
            if (cm%dyU(i, j) /= m%dyBu(i + 1, j + 1)) bad = .true.
            if (cm%uarea(i, j) /= m%areaBu(i + 1, j + 1)) bad = .true.
            if (cm%earea(i, j) /= m%areaCu(i + 1, j)) bad = .true.
            if (cm%narea(i, j) /= m%areaCv(i, j + 1)) bad = .true.
            if (cm%HTE(i, j) /= m%dyCu(i + 1, j)) bad = .true.
            if (cm%HTN(i, j) /= m%dxCv(i, j + 1)) bad = .true.
            ! Masks follow the ocean masks strictly inside the physical
            ! domain (the pinned edge ghost band is `cgrid_masks_edges`).
            if (i >= NG + 1 .and. i <= NG + NXP - 1 .and. j >= NG + 1 .and. j <= NG + NYP - 1) then
               if (cm%epm(i, j) /= m%wet_u(i + 1, j)) bad = .true.
               if (cm%npm(i, j) /= m%wet_v(i, j + 1)) bad = .true.
               if (cm%uvm(i, j) /= m%wet_q(i + 1, j + 1)) bad = .true.
            end if
            if (cm%earea(i, j) > 0.0_wp) then
               if (abs(cm%earear(i, j)*cm%earea(i, j) - 1.0_wp) > 1.0e-14_wp) bad = .true.
            end if
            if (cm%narea(i, j) > 0.0_wp) then
               if (abs(cm%narear(i, j)*cm%narea(i, j) - 1.0_wp) > 1.0e-14_wp) bad = .true.
            end if
         end do
      end do
      call check(error,.not. bad, "ice C-grid name map differs from ocean source")
   end subroutine check_map

   subroutine test_map_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      call grid%init(NXP, NYP, NG, 1.0e4_wp, 1.0e4_wp)
      call make_anisotropic_metrics(m, grid, 1.0e4_wp, 1.2e4_wp, 0.5_wp)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      call check_map(error, cm, m, grid%nx_total, grid%ny_total)
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_map_cartesian

   subroutine test_map_spherical(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_spherical_metrics(m, grid, 0.0_wp, 20.0_wp, DLON, DLAT, REARTH)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      call check_map(error, cm, m, grid%nx_total, grid%ny_total)
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_map_spherical

   subroutine test_map_tripolar(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_tripolar_metrics(m, grid, 0.0_wp, 50.0_wp, DLON, DLAT, REARTH, 60.0_wp, 0.0_wp)
      call make_bc(bc, grid, .true., .true.)
      call cm%init(grid, m, bc)
      call check_map(error, cm, m, grid%nx_total, grid%ny_total)
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_map_tripolar

   subroutine test_ratios(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      integer :: i, j, nx, ny
      real(wp) :: err, tol
      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_spherical_metrics(m, grid, 0.0_wp, 20.0_wp, DLON, DLAT, REARTH)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      nx = grid%nx_total
      ny = grid%ny_total
      err = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            err = max(err, abs(cm%ratiodxN(i, j)*cm%dxN(i, j) + cm%dxN(i + 1, j))/cm%dxN(i + 1, j))
            err = max(err, abs(cm%ratiodxN(i, j)*cm%ratiodxNr(i, j) - 1.0_wp))
         end do
         call check(error, cm%ratiodxN(nx, j), 0.0_wp, "ratiodxN last column must be 0")
         if (allocated(error)) return
      end do
      do j = 1, ny - 1
         do i = 1, nx
            err = max(err, abs(cm%ratiodyE(i, j)*cm%dyE(i, j) + cm%dyE(i, j + 1))/cm%dyE(i, j + 1))
            err = max(err, abs(cm%ratiodyE(i, j)*cm%ratiodyEr(i, j) - 1.0_wp))
         end do
      end do
      do i = 1, nx
         call check(error, cm%ratiodyE(i, ny), 0.0_wp, "ratiodyE last row must be 0")
         if (allocated(error)) return
      end do
      tol = 1.0e-13_wp
      call check(error, err < tol, "ratio identity violated")
      if (allocated(error)) return
      ! lon-lat: dxN varies with latitude only => in-row ratio is exactly -1
      call check(error, abs(cm%ratiodxN(NG + 3, NG + 3) + 1.0_wp) < 1.0e-12_wp, &
                 "lon-lat ratiodxN in a row should be -1")
      call cm%destroy()
      call destroy_cartesian_metrics(m)
      ! Uniform Cartesian grid: the ratios are EXACTLY -1 (single division).
      call make_cartesian_metrics(m, grid)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      do j = 1, ny - 1
         do i = 1, nx - 1
            if (cm%ratiodxN(i, j) /= -1.0_wp .or. cm%ratiodyE(i, j) /= -1.0_wp .or. &
                cm%ratiodxNr(i, j) /= -1.0_wp .or. cm%ratiodyEr(i, j) /= -1.0_wp) then
               call check(error, .false., "uniform-grid ratios must be exactly -1")
               return
            end if
         end do
      end do
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_ratios

   subroutine test_masks_land(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: wet(:, :)
      integer :: i, j, ic, jc
      logical :: bad
      call grid%init(NXP, NYP, NG, 1.0e4_wp, 1.0e4_wp)
      allocate (wet(grid%nx_total, grid%ny_total))
      wet = 1.0_wp
      ic = NG + 5
      jc = NG + 4
      wet(ic, jc) = 0.0_wp
      call make_cartesian_metrics(m, grid, wet_mask=wet)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      bad = .false.
      do j = NG + 1, NG + NYP - 1
         do i = NG + 1, NG + NXP - 1
            ! E(i,j) is open iff T(i,j) and T(i+1,j) are wet; same for N, U.
            if (cm%epm(i, j) /= wet(i, j)*wet(i + 1, j)) bad = .true.
            if (cm%npm(i, j) /= wet(i, j)*wet(i, j + 1)) bad = .true.
            if (cm%uvm(i, j) /= wet(i, j)*wet(i + 1, j)*wet(i, j + 1)*wet(i + 1, j + 1)) bad = .true.
         end do
      end do
      call check(error,.not. bad, "epm/npm/uvm do not follow the T wet mask")
      if (allocated(error)) return
      call check(error, cm%epm(ic, jc) == 0.0_wp .and. cm%epm(ic - 1, jc) == 0.0_wp .and. &
                 cm%npm(ic, jc) == 0.0_wp .and. cm%npm(ic, jc - 1) == 0.0_wp .and. &
                 cm%uvm(ic, jc) == 0.0_wp .and. cm%uvm(ic - 1, jc - 1) == 0.0_wp, &
                 "island faces/corners must be closed")
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_masks_land

   subroutine edge_expectations(error, cm, seam_east)
      !! Hand-derived: Cartesian NXP x NYP domain, all-wet T, edges not
      !! periodic.  Face E(i) is open iff T(i) and T(i+1) are open T cells;
      !! a pinned ghost band is land.  `seam_east`: the east edge is an MPI
      !! seam (has_east false) and follows wet_T (open).
      type(error_type), allocatable, intent(out) :: error
      type(ice_cgrid_metrics_t), intent(in) :: cm
      logical, intent(in) :: seam_east
      integer :: i, j
      logical :: bad
      real(wp) :: want
      bad = .false.
      do j = NG + 1, NG + NYP - 1
         if (cm%epm(NG, j) /= 0.0_wp .or. cm%epm(NG - 1, j) /= 0.0_wp) bad = .true.
         if (cm%uvm(NG, j) /= 0.0_wp) bad = .true.
         do i = NG + 1, NG + NXP - 1
            if (cm%epm(i, j) /= 1.0_wp .or. cm%uvm(i, j) /= 1.0_wp) bad = .true.
         end do
         want = merge(1.0_wp, 0.0_wp, seam_east)
         if (cm%epm(NG + NXP, j) /= want) bad = .true.
      end do
      do i = NG + 1, NG + NXP - 1
         if (cm%npm(i, NG) /= 0.0_wp .or. cm%npm(i, NG + NYP) /= 0.0_wp) bad = .true.
         if (cm%npm(i, NG - 1) /= 0.0_wp) bad = .true.
         if (cm%uvm(i, NG) /= 0.0_wp) bad = .true.
         do j = NG + 1, NG + NYP - 1
            if (cm%npm(i, j) /= 1.0_wp) bad = .true.
         end do
      end do
      call check(error,.not. bad, "edge-pinned masks differ from the hand-derived pattern")
   end subroutine edge_expectations

   subroutine test_masks_edges(error)
      !! Wall edges with `mask_wall_velocity` both ways (ocean wet_u at the
      !! wall face is 1 without it, 0 with it), an OPEN east edge, an MPI
      !! seam (has_east=.false.: not pinned), and a periodic x edge (not
      !! pinned).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      integer :: icase, j
      logical :: solid

      call grid%init(NXP, NYP, NG, 1.0e4_wp, 1.0e4_wp)
      do icase = 1, 5
         solid = (icase == 2)
         call make_cartesian_metrics(m, grid, solid_walls=solid)
         call make_bc(bc, grid, icase == 5, .false.)
         if (icase == 3) bc%east%bc_type = OBC_OPEN
         if (icase == 4) bc%has_east = .false.
         call cm%init(grid, m, bc)
         select case (icase)
         case (1, 2)
            call check(error, m%wet_u(NG + NXP + 1, NG + 2) == merge(0.0_wp, 1.0_wp, solid), &
                       "precondition: ocean wall-face mask")
            if (.not. allocated(error)) call edge_expectations(error, cm, .false.)
         case (3)
            call check(error, m%wet_u(NG + NXP + 1, NG + 2) == 1.0_wp, &
                       "precondition: open edge keeps ocean mask 1")
            if (.not. allocated(error)) call edge_expectations(error, cm, .false.)
         case (4)
            call check(error, m%wet_u(NG + NXP + 1, NG + 2) == 1.0_wp, &
                       "precondition: seam face mask 1")
            if (.not. allocated(error)) call edge_expectations(error, cm, .true.)
         case (5)
            ! periodic x: both x edges follow the wrapped wet_T (open); y walls pinned.
            do j = NG + 1, NG + NYP - 1
               if (allocated(error)) exit
               call check(error, cm%epm(NG, j) == 1.0_wp .and. cm%epm(NG + NXP, j) == 1.0_wp, &
                          "periodic x edge must not be pinned")
            end do
            if (.not. allocated(error)) then
               call check(error, cm%npm(NG + 2, NG) == 0.0_wp .and. cm%npm(NG + 2, NG + NYP) == 0.0_wp, &
                          "y walls must still be pinned under periodic x")
            end if
         end select
         call cm%destroy()
         call destroy_cartesian_metrics(m)
         if (allocated(error)) return
      end do
   end subroutine test_masks_edges

   subroutine test_stretched_lengths(error)
      !! Independent length check: stretched grid; dxE(i) is the
      !! centre-to-centre distance T(i)->T(i+1) from COORDINATES, dyN(j)
      !! likewise in y.  The ocean's face lengths are rebuilt from those
      !! coordinates here, so the test pins the DIRECTION of the shift.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: xc(:), yc(:)
      integer :: i, j, nx, ny
      real(wp) :: err

      call grid%init(NXP, NYP, NG, 1.0e4_wp, 1.0e4_wp)
      call make_cartesian_metrics(m, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (xc(nx), yc(ny))
      do j = 1, ny
         do i = 1, nx
            m%dxT(i, j) = 1.0e4_wp*(1.0_wp + 0.1_wp*real(i, wp))
            m%dyT(i, j) = 1.0e4_wp*(1.0_wp + 0.07_wp*real(j, wp))
         end do
      end do
      xc(1) = 0.5_wp*m%dxT(1, 1)
      do i = 2, nx
         xc(i) = xc(i - 1) + 0.5_wp*(m%dxT(i - 1, 1) + m%dxT(i, 1))
      end do
      yc(1) = 0.5_wp*m%dyT(1, 1)
      do j = 2, ny
         yc(j) = yc(j - 1) + 0.5_wp*(m%dyT(1, j - 1) + m%dyT(1, j))
      end do
      do j = 1, ny
         do i = 2, nx
            m%dxCu(i, j) = xc(i) - xc(i - 1)
         end do
         m%dxCu(1, j) = m%dxCu(2, j)
         m%dxCu(nx + 1, j) = m%dxCu(nx, j)
      end do
      do j = 2, ny
         do i = 1, nx
            m%dyCv(i, j) = yc(j) - yc(j - 1)
         end do
      end do
      do i = 1, nx
         m%dyCv(i, 1) = m%dyCv(i, 2)
         m%dyCv(i, ny + 1) = m%dyCv(i, ny)
      end do
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      err = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            err = max(err, abs(cm%dxE(i, j) - (xc(i + 1) - xc(i)))/cm%dxE(i, j))
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            err = max(err, abs(cm%dyN(i, j) - (yc(j + 1) - yc(j)))/cm%dyN(i, j))
         end do
      end do
      call check(error, err < 1.0e-12_wp, "dxE/dyN is not the centre-to-centre distance to i+1/j+1")
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_stretched_lengths

   subroutine test_fold_ghost_metrics(error)
      !! North-ghost metrics of the ice names: finite, >= 0, and the exact
      !! fold image of the interior.  Cv fold (scalar): i' = 2ng+ni+1-i,
      !! row r' = 2ng+2nj+2-r with r = j+1 the Cv row => N-row
      !! j' = 2ng+2nj-j.  Cu fold: i' (Cu) = 2ng+ni+2-c, row j' = 2ng+2nj+1-j
      !! => E-col i' = 2ng+ni-i.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      integer :: i, j, im, jm, nx, ny, nghost_rows
      logical :: fin, mir
      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_tripolar_metrics(m, grid, 0.0_wp, 50.0_wp, DLON, DLAT, REARTH, 60.0_wp, 0.0_wp)
      call make_bc(bc, grid, .true., .true.)
      call cm%init(grid, m, bc)
      nx = grid%nx_total
      ny = grid%ny_total
      fin = .true.
      mir = .true.
      nghost_rows = 0
      ! N rows strictly above the fold line (r = ng+nj+1 <=> j = ng+nj)
      do j = NG + NYP + 1, ny
         nghost_rows = nghost_rows + 1
         do i = 1, nx
            if (.not. (ieee_is_finite(cm%dxN(i, j)) .and. ieee_is_finite(cm%dyN(i, j)) .and. &
                       ieee_is_finite(cm%narea(i, j)) .and. ieee_is_finite(cm%ratiodxN(i, j)))) fin = .false.
            if (cm%dxN(i, j) < 0.0_wp .or. cm%narea(i, j) < 0.0_wp) fin = .false.
            im = 2*NG + NXP + 1 - i
            jm = 2*NG + 2*NYP - j
            if (im < 1 .or. im > nx) cycle
            if (cm%dxN(i, j) /= cm%dxN(im, jm)) mir = .false.
            if (cm%dyN(i, j) /= cm%dyN(im, jm)) mir = .false.
            if (cm%narea(i, j) /= cm%narea(im, jm)) mir = .false.
            if (cm%HTN(i, j) /= cm%HTN(im, jm)) mir = .false.
         end do
      end do
      call check(error, nghost_rows >= 2, "need >= 2 north N ghost rows")
      if (allocated(error)) return
      call check(error, fin, "non-finite or negative north-ghost N metric")
      if (allocated(error)) return
      call check(error, mir, "north-ghost N metric is not the fold image")
      if (allocated(error)) return
      ! E and U/earea ghosts (rows above the last physical row).
      fin = .true.
      mir = .true.
      do j = NG + NYP + 1, ny
         do i = 1, nx
            if (.not. (ieee_is_finite(cm%dxE(i, j)) .and. ieee_is_finite(cm%dyE(i, j)) .and. &
                       ieee_is_finite(cm%earea(i, j)) .and. ieee_is_finite(cm%uarea(i, j)) .and. &
                       ieee_is_finite(cm%dxU(i, j)))) fin = .false.
            if (cm%earea(i, j) < 0.0_wp .or. cm%uarea(i, j) < 0.0_wp) fin = .false.
            im = 2*NG + NXP - i
            jm = 2*NG + 2*NYP + 1 - j
            if (im < 1 .or. im > nx .or. jm < 1) cycle
            if (cm%dyE(i, j) /= cm%dyE(im, jm)) mir = .false.
            if (cm%dxE(i, j) /= cm%dxE(im, jm)) mir = .false.
            if (cm%earea(i, j) /= cm%earea(im, jm)) mir = .false.
            if (cm%HTE(i, j) /= cm%HTE(im, jm)) mir = .false.
         end do
      end do
      ! U corners: Bu image c' = 2ng+ni+2-c, J' = 2(ng+nj+1)-J
      ! => U-col i' = 2ng+ni-i, U-row j' = 2ng+2nj-j (rows j >= ng+nj+1).
      do j = NG + NYP + 1, ny
         do i = 1, nx
            im = 2*NG + NXP - i
            jm = 2*NG + 2*NYP - j
            if (im < 1 .or. im > nx .or. jm < 1) cycle
            if (cm%dxU(i, j) /= cm%dxU(im, jm)) mir = .false.
            if (cm%dyU(i, j) /= cm%dyU(im, jm)) mir = .false.
            if (cm%uarea(i, j) /= cm%uarea(im, jm)) mir = .false.
         end do
      end do
      call check(error, fin, "non-finite or negative north-ghost E/U metric")
      if (allocated(error)) return
      call check(error, mir, "north-ghost E/U metric is not the fold image")
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_fold_ghost_metrics

   subroutine test_corner_wrap(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: f(:, :)
      integer :: nxf, nyf, i, j, ii, jj
      logical :: bad
      nxf = NXP + 2*NG + 1
      nyf = NYP + 2*NG + 1
      allocate (f(nxf, nyf))
      f = -1.0e30_wp
      ! interior corners i = ng+1..ng+nx (the seam corner ng+nx+1 is set by the wrap)
      do j = NG + 1, NG + NYP
         do i = NG + 1, NG + NXP
            f(i, j) = real(100*(i - NG) + (j - NG), wp)
         end do
      end do
      !$acc enter data copyin(f)
      call evp_wrap_corner_impl(f, nxf, nyf, NXP, NYP, NG, .true., .false.)
      !$acc update self(f)
      !$acc exit data delete(f)
      bad = .false.
      do j = NG + 1, NG + NYP
         do i = 1, nxf
            ii = NG + 1 + modulo(i - NG - 1, NXP)
            if (f(i, j) /= f(ii, j)) bad = .true.
         end do
      end do
      call check(error,.not. bad, "corner x-wrap is not the periodic image")
      if (allocated(error)) return
      ! both axes
      f = -1.0e30_wp
      do j = NG + 1, NG + NYP
         do i = NG + 1, NG + NXP
            f(i, j) = real(100*(i - NG) + (j - NG), wp)
         end do
      end do
      !$acc enter data copyin(f)
      call evp_wrap_corner_impl(f, nxf, nyf, NXP, NYP, NG, .true., .true.)
      !$acc update self(f)
      !$acc exit data delete(f)
      bad = .false.
      do j = 1, nyf
         do i = 1, nxf
            ii = NG + 1 + modulo(i - NG - 1, NXP)
            jj = NG + 1 + modulo(j - NG - 1, NYP)
            if (f(i, j) /= real(100*(ii - NG) + (jj - NG), wp)) bad = .true.
         end do
      end do
      call check(error,.not. bad, "corner xy-wrap is not the periodic image")
   end subroutine test_corner_wrap

   ! -----------------------------------------------------------------
   ! halo wrapper == the ocean's own kernels, bitwise
   ! -----------------------------------------------------------------

   subroutine seed_faces(u, v, nx, ny)
      !! Index-encoded values everywhere (ghosts included, non-periodic on
      !! purpose so a missing wrap/fold shows).
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: u(nx + 1, ny), v(nx, ny + 1)
      integer :: i, j
      do j = 1, ny
         do i = 1, nx + 1
            u(i, j) = 1.0_wp + real(i, wp)*0.37_wp + real(j, wp)*0.011_wp
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            v(i, j) = -2.0_wp + real(i, wp)*0.013_wp + real(j, wp)*0.29_wp
         end do
      end do
   end subroutine seed_faces

   subroutine run_halo_case(error, periodic_x, fold)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: periodic_x, fold
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: u(:, :), v(:, :), u2(:, :), v2(:, :)
      integer :: nx, ny, i, j, im, jm
      logical :: bad

      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_bc(bc, grid, periodic_x, fold)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (u(nx + 1, ny), v(nx, ny + 1), u2(nx + 1, ny), v2(nx, ny + 1))
      call seed_faces(u, v, nx, ny)
      u2 = u
      v2 = v

      !$acc enter data copyin(u, v, u2, v2)
      call ice_halo_faces(grid, bc, u, v)
      ! The ocean's own sequence, spelled out kernel by kernel.
      call ocean_periodic_wrap_face_x_2d(u2, nx + 1, ny, NXP, NYP, NG, bc%periodic_x, bc%periodic_y)
      call ocean_periodic_wrap_face_y_2d(v2, nx, ny + 1, NXP, NYP, NG, bc%periodic_x, bc%periodic_y)
      if (bc%north_fold) then
         call fold_north_u_face(u2, nx + 1, ny, NXP, NYP, NG)
         call fold_north_v_face(v2, nx, ny + 1, NXP, NYP, NG)
      end if
      !$acc update self(u, v, u2, v2)
      !$acc exit data delete(u, v, u2, v2)

      call check(error, all(u == u2) .and. all(v == v2), &
                 "ice_halo_faces differs from the ocean's wrap + fold")
      if (allocated(error)) return

      if (periodic_x) then
         bad = .false.
         do j = NG + 1, NG + NYP
            do i = 1, NG
               if (u(i, j) /= u(i + NXP, j)) bad = .true.
               if (v(i, j) /= v(i + NXP, j)) bad = .true.
            end do
         end do
         call check(error,.not. bad, "periodic ghost faces are not the interior image")
         if (allocated(error)) return
      else
         call seed_faces(u2, v2, nx, ny)
         call check(error, all(u == u2) .and. all(v == v2), "wall case modified the faces")
         if (allocated(error)) return
      end if

      if (fold) then
         bad = .false.
         do j = NG + NYP + 1, ny
            do i = 1, nx + 1
               im = 2*NG + NXP + 2 - i
               jm = 2*NG + 2*NYP + 1 - j
               if (im >= 1 .and. im <= nx + 1) then
                  if (u(i, j) /= -u(im, jm)) bad = .true.
               end if
            end do
         end do
         call check(error,.not. bad, "u fold ghosts are not the negated mirror")
         if (allocated(error)) return
         bad = .false.
         j = NG + NYP + 1
         do i = 1, nx
            im = 2*NG + NXP + 1 - i
            if (im >= 1 .and. im <= nx) then
               if (v(i, j) /= -v(im, j)) bad = .true.
            end if
         end do
         call check(error,.not. bad, "v fold-row is not antisymmetric")
         if (allocated(error)) return
         ! v ghost rows ABOVE the fold line: v(i,r) = -v(2ng+ni+1-i, 2(ng+nj+1)-r)
         bad = .false.
         do j = NG + NYP + 2, ny + 1
            do i = 1, nx
               im = 2*NG + NXP + 1 - i
               jm = 2*(NG + NYP + 1) - j
               if (im >= 1 .and. im <= nx .and. jm >= 1) then
                  if (v(i, j) /= -v(im, jm)) bad = .true.
               end if
            end do
         end do
         call check(error,.not. bad, "v fold ghost rows are not the negated mirror")
      end if
   end subroutine run_halo_case

   subroutine test_halo_periodic(error)
      type(error_type), allocatable, intent(out) :: error
      call run_halo_case(error, periodic_x=.true., fold=.false.)
   end subroutine test_halo_periodic

   subroutine test_halo_wall(error)
      type(error_type), allocatable, intent(out) :: error
      call run_halo_case(error, periodic_x=.false., fold=.false.)
   end subroutine test_halo_wall

   subroutine test_halo_fold(error)
      type(error_type), allocatable, intent(out) :: error
      call run_halo_case(error, periodic_x=.true., fold=.true.)
   end subroutine test_halo_fold

   ! -----------------------------------------------------------------
   ! device residency of the mapped metric arrays
   ! -----------------------------------------------------------------

   subroutine ratio_kernel(nx, ny, rx, rxr, ry, ryr, epm, npm, uvm, out)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: rx(nx, ny), rxr(nx, ny), ry(nx, ny), ryr(nx, ny)
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny), uvm(nx, ny)
      real(wp), intent(inout) :: out(nx, ny)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         out(i, j) = rx(i, j)*rxr(i, j) + ry(i, j)*ryr(i, j) + epm(i, j) + npm(i, j) + uvm(i, j)
      end do
   end subroutine ratio_kernel

   subroutine test_device_kernel(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ice_cgrid_metrics_t) :: cm
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: out(:, :), ref(:, :)
      integer :: nx, ny
      call grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_spherical_metrics(m, grid, 0.0_wp, 20.0_wp, DLON, DLAT, REARTH)
      call make_bc(bc, grid, .false., .false.)
      call cm%init(grid, m, bc)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (out(nx, ny), ref(nx, ny))
      out = -7.0_wp
      ref = cm%ratiodxN*cm%ratiodxNr + cm%ratiodyE*cm%ratiodyEr + cm%epm + cm%npm + cm%uvm
      call cm%enter_data()
      !$acc enter data copyin(out)
      call ratio_kernel(nx, ny, cm%ratiodxN, cm%ratiodxNr, cm%ratiodyE, cm%ratiodyEr, &
                        cm%epm, cm%npm, cm%uvm, out)
      !$acc update self(out)
      !$acc exit data delete(out)
      call cm%exit_data()
      ! Sum of products: GPU builds contract FMA, so compare to a tolerance.
      call check(error, maxval(abs(out - ref)) < 1.0e-13_wp, "device kernel disagrees with host")
      if (allocated(error)) return
      call check(error, cm%bytes() == int(20*nx*ny*8, kind(cm%bytes())), "bytes() count")
      call cm%destroy()
      call destroy_cartesian_metrics(m)
   end subroutine test_device_kernel

end module test_cice_cgrid_metrics
