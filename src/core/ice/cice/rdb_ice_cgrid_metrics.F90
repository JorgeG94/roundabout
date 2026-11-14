module rdb_ice_cgrid_metrics
   !! C-grid sea-ice dynamics metrics: the name map from the ocean's
   !! `ocean_metrics_t` to the ice C-grid names (Lemieux et al. 2024, the
   !! C-grid EVP formulation; Bouillon et al. 2013 for the EVP rheology
   !! being discretised), the `ratiodxN`/`ratiodyE` lateral-boundary ratios,
   !! the `epm`/`npm`/`uvm` masks, and the ONE-exchange-per-subcycle halo
   !! wrapper `ice_halo_faces` for the face velocity pair.
   !!
   !! Roundabout-original; no physics lives here (D9.1 added the T-cell names
   !! `dxT, dyT, tarea, tarear, hm` the EVP needs).  It needs the ocean's
   !! grid metrics, halo and fold, so it is part of the MAIN library (not
   !! the standalone `rdb_cice_objs`).
   !!
   !! **Index map (load-bearing, pinned by `test_cice_cgrid_metrics`).**
   !! Ice-side names are `(nx_total, ny_total)` arrays addressed by the cell
   !! `(i,j)` they belong to:
   !!   E(i,j) = EAST face of T(i,j)   = ocean x-face (Cu)  (i+1, j)
   !!   N(i,j) = NORTH face of T(i,j)  = ocean y-face (Cv)  (i, j+1)
   !!   U(i,j) = NE corner of T(i,j)   = ocean corner (Bu)  (i+1, j+1)
   !! (the ocean's `Cu(i,j)` is the WEST face of cell i, `Cv(i,j)` the SOUTH
   !! face, `Bu(i,j)` the SW corner).  So the ice velocity `u_ice`/`v_ice`
   !! IS the ocean face arrays: ice `uvelE(i,j) = u_ice(i+1,j)`,
   !! `vvelN(i,j) = v_ice(i,j+1)`; the arrays are not copied or shifted.
   !!
   !! | ice name           | ocean source                          |
   !! |--------------------|---------------------------------------|
   !! | `dxE`, `dyE`       | `dxCu`, `dyCu` at (i+1,j)             |
   !! | `dxN`, `dyN`       | `dxCv`, `dyCv` at (i,j+1)             |
   !! | `dxU`, `dyU`       | `dxBu`, `dyBu` at (i+1,j+1)           |
   !! | `dxT`, `dyT`        | `dxT`, `dyT` at (i,j)                 |
   !! | `tarea`, `tarear`  | `areaT` at (i,j) and its reciprocal   |
   !! | `hm`               | pinned T mask (real 0/1)              |
   !! | `uarea`            | `areaBu` at (i+1,j+1) (weight only)   |
   !! | `earea`, `narea`   | `areaCu` (i+1,j), `areaCv` (i,j+1)    |
   !! | `HTE`, `HTN`       | `dyCu` (i+1,j), `dxCv` (i,j+1)        |
   !! | `epm`              | `wet_u` (i+1,j) x pinned T mask       |
   !! | `npm`              | `wet_v` (i,j+1) x pinned T mask       |
   !! | `uvm`              | `wet_q` (i+1,j+1) x pinned T mask     |
   !!
   !! `ratiodxN(i,j) = -dxN(i+1,j)/dxN(i,j)` and
   !! `ratiodyE(i,j) = -dyE(i,j+1)/dyE(i,j)` (and their reciprocals
   !! `ratiodxNr`, `ratiodyEr`; one guarded division each, so a uniform grid
   !! gives exactly -1) carry the lateral boundary condition of the
   !! C-grid strain rates (no-slip: the mirrored tangential velocity makes the
   !! wall-parallel flow vanish at a straight wall).  They are defined for
   !! `i <= nx_total-1` (resp. `j <= ny_total-1`); the last column/row is
   !! 0, as is any entry with a zero denominator (Adcroft convention) --
   !! with `nghost >= 2` every entry a ghost-ring-1 stencil reads is valid.
   !!
   !! **Halo and fold.**  `ice_halo_faces` refreshes `(u_ice, v_ice)`: on a
   !! locally wrapped periodic axis the DUPLICATE seam face (`u(is,:)`,
   !! `v(:,js)`) is first copied from the face the solver computes (`u(ie+1,:)`,
   !! `v(:,je+1)`; the ocean's wrap goes the other way), then the ocean's local
   !! periodic wrap fills the ghost faces, then the tripolar fold (true
   !! vector: negate, v fold row projected antisymmetric) runs, a no-op off the
   !! north rank row (`bc%north_fold` already carries `has_north`).  A
   !! decomposed axis stops (`error stop`) until the multi-rank PR adds the
   !! shifted one-face exchange.  The corner-staggered
   !! persistent stress is wrapped with the public `evp_wrap_corner_impl`
   !! of `rdb_ice_evp`.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_y_2d
   use rdb_ocean_halo, only: ocean_halo_face_x, ocean_halo_face_y, &
                             ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_face_pair
   implicit none
   private

   public :: ice_cgrid_metrics_t
   public :: ice_halo_faces

   type :: ice_cgrid_metrics_t
      !! Ice C-grid metric set, all `(nx_total, ny_total)`, built once.
      logical :: is_init = .false.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data`.
      real(wp), allocatable :: dxE(:, :), dyE(:, :)
         !! E-face lengths (m).
      real(wp), allocatable :: dxN(:, :), dyN(:, :)
         !! N-face lengths (m).
      real(wp), allocatable :: dxU(:, :), dyU(:, :)
         !! U-corner lengths (m).
      real(wp), allocatable :: dxT(:, :), dyT(:, :)
         !! T-cell widths through the middle (m).
      real(wp), allocatable :: tarea(:, :), tarear(:, :)
         !! T-cell area (m^2) and its Adcroft reciprocal.
      real(wp), allocatable :: hm(:, :)
         !! T open mask (real 0/1; the ghost band of a wall edge pinned to land).
      real(wp), allocatable :: uarea(:, :)
         !! U-corner area (m^2).
      real(wp), allocatable :: earea(:, :), narea(:, :)
         !! E / N face-cell areas (m^2).
      real(wp), allocatable :: earear(:, :), narear(:, :)
         !! 1/earea, 1/narea (Adcroft reciprocal).
      real(wp), allocatable :: HTE(:, :), HTN(:, :)
         !! Face widths (m): `HTE = dyCu` (E face length), `HTN = dxCv`.
      real(wp), allocatable :: ratiodxN(:, :), ratiodxNr(:, :)
         !! `-dxN(i+1,j)/dxN(i,j)` and its reciprocal.
      real(wp), allocatable :: ratiodyE(:, :), ratiodyEr(:, :)
         !! `-dyE(i,j+1)/dyE(i,j)` and its reciprocal.
      real(wp), allocatable :: epm(:, :), npm(:, :), uvm(:, :)
         !! E / N / U open masks (real 0/1).
   contains
      procedure :: init => ice_cgrid_metrics_init
      procedure :: destroy => ice_cgrid_metrics_destroy
      procedure :: bytes => ice_cgrid_metrics_bytes
      procedure :: enter_data => ice_cgrid_metrics_enter_data
      procedure :: exit_data => ice_cgrid_metrics_exit_data
   end type ice_cgrid_metrics_t

contains

   pure function safe_recip(x) result(r)
      !! Adcroft reciprocal: 0 for a zero denominator.
      real(wp), intent(in) :: x
      real(wp) :: r
      if (x /= 0.0_wp) then
         r = 1.0_wp/x
      else
         r = 0.0_wp
      end if
   end function safe_recip

   subroutine ice_cgrid_metrics_init(this, grid, metrics, bc, ierr)
      !! Build the ice names from the ocean metrics (host, configure time).
      !! Call after the ocean metrics are final (land mask applied, ghosts
      !! wrapped/folded) and BEFORE `enter_data`.
      !!
      !! Re-initialising is refused: `init` on an object that is already
      !! initialised stops with a message (call `exit_data` then `destroy`
      !! first), so a device-mapped array is never deallocated under the
      !! device.  `nghost < 2` is refused (the ring-1 strain stencil needs
      !! it): `ierr = 1` if present, else `error stop`.
      !!
      !! Masks: built from a T mask whose ghost band is pinned to LAND on a
      !! physical, non-periodic, non-fold edge (`bc%has_*`, so an MPI seam
      !! is NOT pinned and follows `wet_T`), independent of
      !! `mask_wall_velocity` and of what an OPEN edge's ghost `wet_T`
      !! holds: `epm = wet_u * Tm(i,j)*Tm(i+1,j)`, `npm` likewise in y,
      !! `uvm = wet_q * (4 surrounding Tm)`.  The ice never writes ghosts at
      !! a wall/open edge, so those faces must read closed.
      !!
      !! `uarea` is the local-lengths product from the ocean's corner area
      !! (`areaBu`).  The upstream doubling of the corner sub-cell area on
      !! the fold row is deliberately NOT reproduced (plan decision M-D9-1,
      !! MASTER_PLAN B.7).
      class(ice_cgrid_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_bc_state_t), intent(in) :: bc
      integer, intent(out), optional :: ierr

      integer :: nx, ny, i, j, ng, ni, nj
      logical :: land_w, land_e, land_s, land_n, wrap_x, wrap_y
      real(wp), allocatable :: tm(:, :)

      if (present(ierr)) ierr = 0
      if (this%is_init) then
         error stop "ice_cgrid_metrics_t%init: already initialised; exit_data + destroy first"
      end if
      if (grid%nghost < 2) then
         if (present(ierr)) then
            ierr = 1
            return
         end if
         error stop "ice_cgrid_metrics_t%init: C-grid ice dynamics requires nghost >= 2"
      end if

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost
      ni = grid%nx_phys
      nj = grid%ny_phys
      allocate (this%dxE(nx, ny), this%dyE(nx, ny), this%dxN(nx, ny), this%dyN(nx, ny))
      allocate (this%dxU(nx, ny), this%dyU(nx, ny), this%uarea(nx, ny))
      allocate (this%dxT(nx, ny), this%dyT(nx, ny), this%tarea(nx, ny), this%tarear(nx, ny), this%hm(nx, ny))
      allocate (this%earea(nx, ny), this%narea(nx, ny))
      allocate (this%earear(nx, ny), this%narear(nx, ny))
      allocate (this%HTE(nx, ny), this%HTN(nx, ny))
      allocate (this%ratiodxN(nx, ny), this%ratiodxNr(nx, ny))
      allocate (this%ratiodyE(nx, ny), this%ratiodyEr(nx, ny))
      allocate (this%epm(nx, ny), this%npm(nx, ny), this%uvm(nx, ny))

      ! T mask with the physical-edge ghost band pinned to land.
      land_w = bc%has_west .and. .not. bc%periodic_x
      land_e = bc%has_east .and. .not. bc%periodic_x
      land_s = bc%has_south .and. .not. bc%periodic_y
      land_n = bc%has_north .and. .not. bc%periodic_y .and. .not. bc%north_fold
      allocate (tm(nx, ny))
      do j = 1, ny
         do i = 1, nx
            if ((i <= ng .and. land_w) .or. (i > ng + ni .and. land_e) .or. &
                (j <= ng .and. land_s) .or. (j > ng + nj .and. land_n)) then
               tm(i, j) = 0.0_wp
            else
               tm(i, j) = merge(1.0_wp, 0.0_wp, metrics%wet_T(i, j) > 0.5_wp)
            end if
         end do
      end do
      wrap_x = bc%periodic_x .and. .not. ocean_halo_is_decomposed_x()
      wrap_y = bc%periodic_y .and. .not. ocean_halo_is_decomposed_y()
      call ocean_periodic_wrap_centre_2d(tm, nx, ny, ni, nj, ng, wrap_x, wrap_y)

      this%epm = 0.0_wp
      this%npm = 0.0_wp
      this%uvm = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            this%dxE(i, j) = metrics%dxCu(i + 1, j)
            this%dyE(i, j) = metrics%dyCu(i + 1, j)
            this%dxN(i, j) = metrics%dxCv(i, j + 1)
            this%dyN(i, j) = metrics%dyCv(i, j + 1)
            this%dxU(i, j) = metrics%dxBu(i + 1, j + 1)
            this%dyU(i, j) = metrics%dyBu(i + 1, j + 1)
            this%uarea(i, j) = metrics%areaBu(i + 1, j + 1)
            this%dxT(i, j) = metrics%dxT(i, j)
            this%dyT(i, j) = metrics%dyT(i, j)
            this%tarea(i, j) = metrics%areaT(i, j)
            this%tarear(i, j) = safe_recip(this%tarea(i, j))
            this%hm(i, j) = tm(i, j)
            this%earea(i, j) = metrics%areaCu(i + 1, j)
            this%narea(i, j) = metrics%areaCv(i, j + 1)
            this%earear(i, j) = safe_recip(this%earea(i, j))
            this%narear(i, j) = safe_recip(this%narea(i, j))
            this%HTE(i, j) = metrics%dyCu(i + 1, j)
            this%HTN(i, j) = metrics%dxCv(i, j + 1)
         end do
      end do
      ! The last column / row has no east / north T neighbour in the array:
      ! those faces and corners stay closed (0).
      do j = 1, ny
         do i = 1, nx - 1
            this%epm(i, j) = metrics%wet_u(i + 1, j)*tm(i, j)*tm(i + 1, j)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            this%npm(i, j) = metrics%wet_v(i, j + 1)*tm(i, j)*tm(i, j + 1)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx - 1
            this%uvm(i, j) = metrics%wet_q(i + 1, j + 1)* &
                             tm(i, j)*tm(i + 1, j)*tm(i, j + 1)*tm(i + 1, j + 1)
         end do
      end do

      ! Single guarded division: a uniform grid gives exactly -1.
      this%ratiodxN = 0.0_wp
      this%ratiodxNr = 0.0_wp
      this%ratiodyE = 0.0_wp
      this%ratiodyEr = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            if (this%dxN(i, j) /= 0.0_wp) then
               this%ratiodxN(i, j) = -this%dxN(i + 1, j)/this%dxN(i, j)
            end if
            if (this%ratiodxN(i, j) /= 0.0_wp) then
               this%ratiodxNr(i, j) = 1.0_wp/this%ratiodxN(i, j)
            end if
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            if (this%dyE(i, j) /= 0.0_wp) then
               this%ratiodyE(i, j) = -this%dyE(i, j + 1)/this%dyE(i, j)
            end if
            if (this%ratiodyE(i, j) /= 0.0_wp) then
               this%ratiodyEr(i, j) = 1.0_wp/this%ratiodyE(i, j)
            end if
         end do
      end do
      this%is_init = .true.
   end subroutine ice_cgrid_metrics_init

   subroutine ice_cgrid_metrics_destroy(this)
      class(ice_cgrid_metrics_t), intent(inout) :: this
      if (this%mapped) error stop "ice_cgrid_metrics_t%destroy: exit_data first"
      if (allocated(this%dxE)) deallocate (this%dxE, this%dyE, this%dxN, this%dyN)
      if (allocated(this%dxU)) deallocate (this%dxU, this%dyU, this%uarea)
      if (allocated(this%dxT)) deallocate (this%dxT, this%dyT, this%tarea, this%tarear, this%hm)
      if (allocated(this%earea)) deallocate (this%earea, this%narea, this%earear, this%narear)
      if (allocated(this%HTE)) deallocate (this%HTE, this%HTN)
      if (allocated(this%ratiodxN)) deallocate (this%ratiodxN, this%ratiodxNr, &
                                                this%ratiodyE, this%ratiodyEr)
      if (allocated(this%epm)) deallocate (this%epm, this%npm, this%uvm)
      this%is_init = .false.
   end subroutine ice_cgrid_metrics_destroy

   pure function arr_bytes(a) result(nb)
      real(wp), allocatable, intent(in) :: a(:, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes

   pure function ice_cgrid_metrics_bytes(this) result(nbytes)
      !! Counted allocatable footprint (0 when unallocated).  Add a term when
      !! a new allocatable joins the type.
      class(ice_cgrid_metrics_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%dxE) &
               + arr_bytes(this%dyE) &
               + arr_bytes(this%dxN) &
               + arr_bytes(this%dyN) &
               + arr_bytes(this%dxU) &
               + arr_bytes(this%dyU) &
               + arr_bytes(this%uarea) &
               + arr_bytes(this%dxT) &
               + arr_bytes(this%dyT) &
               + arr_bytes(this%tarea) &
               + arr_bytes(this%tarear) &
               + arr_bytes(this%hm) &
               + arr_bytes(this%earea) &
               + arr_bytes(this%narea) &
               + arr_bytes(this%earear) &
               + arr_bytes(this%narear) &
               + arr_bytes(this%HTE) &
               + arr_bytes(this%HTN) &
               + arr_bytes(this%ratiodxN) &
               + arr_bytes(this%ratiodxNr) &
               + arr_bytes(this%ratiodyE) &
               + arr_bytes(this%ratiodyEr) &
               + arr_bytes(this%epm) &
               + arr_bytes(this%npm) &
               + arr_bytes(this%uvm)
   end function ice_cgrid_metrics_bytes

   subroutine ice_cgrid_metrics_enter_data(this)
      !! Map every component array to the device (`copyin`).
      class(ice_cgrid_metrics_t), intent(inout) :: this
      !$acc enter data copyin(this%dxE, this%dyE, this%dxN, this%dyN, this%dxU, this%dyU)
      !$acc enter data copyin(this%uarea, this%earea, this%narea, this%earear, this%narear)
      !$acc enter data copyin(this%HTE, this%HTN)
      !$acc enter data copyin(this%dxT, this%dyT, this%tarea, this%tarear, this%hm)
      !$acc enter data copyin(this%ratiodxN, this%ratiodxNr, this%ratiodyE, this%ratiodyEr)
      !$acc enter data copyin(this%epm, this%npm, this%uvm)
      this%mapped = .true.
   end subroutine ice_cgrid_metrics_enter_data

   subroutine ice_cgrid_metrics_exit_data(this)
      class(ice_cgrid_metrics_t), intent(inout) :: this
      !$acc exit data delete(this%dxE, this%dyE, this%dxN, this%dyN, this%dxU, this%dyU)
      !$acc exit data delete(this%uarea, this%earea, this%narea, this%earear, this%narear)
      !$acc exit data delete(this%HTE, this%HTN)
      !$acc exit data delete(this%dxT, this%dyT, this%tarea, this%tarear, this%hm)
      !$acc exit data delete(this%ratiodxN, this%ratiodxNr, this%ratiodyE, this%ratiodyEr)
      !$acc exit data delete(this%epm, this%npm, this%uvm)
      this%mapped = .false.
   end subroutine ice_cgrid_metrics_exit_data

   subroutine ice_halo_faces(grid, bc, u, v, device_resident)
      !! One seam refresh of the ice face velocity pair `(u, v)` =
      !! `(u_ice, v_ice)`: duplicate -> periodic wrap -> tripolar fold (true
      !! vector).  The ocean's own sequence (exchange -> wrap -> fold) is
      !! reused, with ONE difference forced by the ice solver.
      !!
      !! **Which face is authoritative.**  The ice solver computes the faces
      !! `is+1 .. ie+1` (the east face of each owned cell, `E(is:ie)`) and the
      !! north ones `js+1 .. je+1`; the face `u(is, :)` (`v(:, js)`) is the
      !! DUPLICATE of the east (north) seam face.  The ocean's periodic wrap
      !! takes the opposite view (`u(ie+1) := u(is)`), so applied as is it
      !! would overwrite the freshly computed seam face with a value nothing
      !! updates.  On a locally wrapped axis this helper therefore first copies
      !! the computed face onto its duplicate (`u(is,:) = u(ie+1,:)`,
      !! `v(:,js) = v(:,je+1)`), after which the ocean's wrap is a no-op on the
      !! seam and fills the ghost faces.
      !!
      !! A DECOMPOSED axis needs the shifted one-face exchange of the
      !! multi-rank PR (the ocean exchange duplicates the other way too) and
      !! stops here until then.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(inout) :: u(grid%nx_total + 1, grid%ny_total)
         !! x-face field.
      real(wp), intent(inout) :: v(grid%nx_total, grid%ny_total + 1)
         !! y-face field.
      logical, intent(in), optional :: device_resident
         !! `.false.` for host-side calls (default: device-resident).

      logical :: halo_x, halo_y, wrap_x, wrap_y
      integer :: nx, ny, i, j, is, ie, js, je

      nx = grid%nx_total
      ny = grid%ny_total
      is = grid%nghost + 1
      ie = grid%nghost + grid%nx_phys
      js = grid%nghost + 1
      je = grid%nghost + grid%ny_phys
      halo_x = ocean_halo_is_decomposed_x()
      halo_y = ocean_halo_is_decomposed_y()
      if (halo_x .or. halo_y) then
         error stop "ice_halo_faces: a decomposed axis needs the shifted one-face exchange "// &
            "of the multi-rank PR (not available yet)"
      end if
      wrap_x = bc%periodic_x
      wrap_y = bc%periodic_y
      if (wrap_x) then
         do concurrent(j=1:ny)
            u(is, j) = u(ie + 1, j)
         end do
      end if
      if (wrap_y) then
         do concurrent(i=1:nx)
            v(i, js) = v(i, je + 1)
         end do
      end if
      call ocean_periodic_wrap_face_x_2d(u, nx + 1, ny, grid%nx_phys, grid%ny_phys, &
                                         grid%nghost, wrap_x, wrap_y)
      call ocean_periodic_wrap_face_y_2d(v, nx, ny + 1, grid%nx_phys, grid%ny_phys, &
                                         grid%nghost, wrap_x, wrap_y)
      call ocean_fold_wrap_face_pair(grid, bc, u, v, device_resident)
   end subroutine ice_halo_faces

end module rdb_ice_cgrid_metrics
