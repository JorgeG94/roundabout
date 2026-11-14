!! The CICE6 C-grid EVP (`rdb_ice_evp_c`) under the tripolar north fold and on
!! the sphere (D9.3), ice-only: raw ice-name arrays and analytic metrics, no
!! ocean.  Hunke and Dukowicz (1997, 2002), Kimmritz et al. (2016) and Lemieux
!! et al. (2024) define the scheme; Murray (1996) the fold.
!!
!! **The fold as a symmetry.**  The folded domain (`nj` physical rows, the north
!! edge glued to itself by `(x, y) -> (ni - x, 2 nj - y)`) is the quotient of a
!! DOUBLED domain (`2 nj` rows, walls south and north, periodic in x) by the
!! 180-degree rotation about `(ni/2, nj)`.  So a state of the doubled domain
!! that is invariant under the rotation (scalars copy, vectors negate) must stay
!! invariant, and the folded run must reproduce the lower half of the doubled
!! run.  Every array of the fixtures is built by ONE rule (`rot_fill`): the rows
!! up to the fold line come from an analytic formula, the rest is the exact
!! rotated copy of them (the fold-line row of a stagger that has one: the east
!! half is the formula, the west half its image), so the doubled state and the
!! folded ghost rows are bit-for-bit images.  The map indices are an
!! INDEPENDENT host statement of the fold (a re-derivation from the continuous
!! coordinates, not a call of the code under test).
!!
!! Cases:
!!   * `fold_kernel_matches_host_map`    the fold kernel (1 launch) equals the
!!                                       host rotation on random faces, bitwise;
!!                                       idempotent; the fold-line row is
!!                                       antisymmetric, the self-conjugate column of an odd ni is 0
!!                                       (ni = 12 and 13).
!!   * `doubled_domain_symmetry`         rotation-invariant state, doubled domain,
!!                                       N steps: asymmetry measured (bitwise or the
!!                                       round-off bound, and which operation).
!!   * `fold_equals_doubled_domain`      folded run == lower half of the doubled run
!!                                       (u, v, the five stresses, strength,
!!                                       deformations), non-uniform metrics, an
!!                                       island next to the fold.
!!   * `corner_storage_mirror_formulas`  the mirror index of every corner storage (stress12U SW storage,
!!                                       shearU/deltaU NE storage) and its fold row, from the continuous map.
!!   * `zero_width_face_pole`            a zero-length closed face between two wet T cells (a bipole of
!!                                       the tripolar grid, where the ocean metrics are degenerate):
!!                                       finite stresses and velocities (was 0/0 in the T strain rates).
!!   * `ice_crossing_fold`               a patch driven across the fold by a vortex
!!                                       wind: flux through the fold telescopes,
!!                                       ghost rows are exact images, the ghost-ring
!!                                       stresses agree with their mirrored owned
!!                                       values, subcycle by subcycle.
!!   * `fold_line_dofs_consistent`       after EVERY subcycle: the duplicated v row
!!                                       is antisymmetric (exactly), the fold is
!!                                       idempotent, the duplicated corner stress
!!                                       slots agree.
!!   * `solid_body_rotation_sphere`      lon-lat metrics (converging meridians): the
!!                                       rigid rotation has exactly zero strain and
!!                                       zero stress divergence with metric-consistent
!!                                       corner lengths; with the natural spherical
!!                                       lengths the corner shear converges at O(h^2).
!!   * `free_drift_sphere`               latitude-varying Coriolis parameter: the
!!                                       scheme's free-drift velocity against the local
!!                                       analytic balance, second order in the cell size.
module test_cice_evp_c_fold
   use, intrinsic :: iso_fortran_env, only: int32, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_VISC_AVG_STRENGTH
   use rdb_ice_evp_c_prep, only: ice_evp_c_const_t, ice_evp_c_work_t, ice_evp_c_constants, ice_evp_c_prep
   use rdb_ice_evp_c, only: ice_evp_c_step, ice_evp_c_substep, ice_evp_c_average, ice_evp_c_finish, &
                            ice_evp_c_wrap_faces, ice_evp_c_fold_faces, ice_evp_c_copy_faces, &
                            ice_evp_c_div_stress, evp_strain_u_impl, evp_stress_t_impl, evp_stress_u_impl
   implicit none
   private

   public :: collect_cice_evp_c_fold_tests

   integer, parameter :: NG = 3
   real(wp), parameter :: PI = 3.14159265358979323846_wp
   real(wp), parameter :: TWO_PI = 2.0_wp*PI
   real(wp), parameter :: OM_EARTH = 7.292e-5_wp
   real(wp), parameter :: DELTAMIN = 1.0e-11_wp
      !! `deltaminEVP`, the default of `ice_dyn_params_default`.

   type :: ft_t
      !! Ice-name fixture (the arrays of `ice_evp_c_step`), `ni` columns, `nyp` rows.
      integer :: ni = 0, nj = 0, nyp = 0, nx = 0, ny = 0, is = 0, ie = 0, js = 0, je = 0
      real(wp) :: dt = 900.0_wp
      real(wp), allocatable :: dxE(:, :), dyE(:, :), dxN(:, :), dyN(:, :), dxU(:, :), dyU(:, :), dxT(:, :), dyT(:, :)
      real(wp), allocatable :: uarea(:, :), tarea(:, :), tarear(:, :), earea(:, :), narea(:, :)
      real(wp), allocatable :: earear(:, :), narear(:, :), rdxN(:, :), rdxNr(:, :), rdyE(:, :), rdyEr(:, :)
      real(wp), allocatable :: epm(:, :), npm(:, :), uvm(:, :), hm(:, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), aice_init(:, :), aice0(:, :)
      real(wp), allocatable :: aicen(:, :, :), vicen(:, :, :)
      real(wp), allocatable :: fcorE(:, :), fcorN(:, :), strairxE(:, :), strairyN(:, :), ss_tltxE(:, :), ss_tltyN(:, :)
      real(wp), allocatable :: uocnE(:, :), vocnE(:, :), uocnN(:, :), vocnN(:, :), u(:, :), v(:, :)
      integer(int32), allocatable :: ice_err(:, :)
      type(ice_evp_c_work_t) :: w
      type(ice_params_t) :: p
      type(ice_evp_c_const_t) :: k
      logical :: mapped = .false.
   end type ft_t

contains

   subroutine collect_cice_evp_c_fold_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fold_kernel_matches_host_map", test_fold_kernel), &
                  new_unittest("doubled_domain_symmetry", test_doubled_symmetry), &
                  new_unittest("fold_equals_doubled_domain", test_fold_equals_doubled), &
                  new_unittest("corner_storage_mirror_formulas", test_corner_formulas), &
                  new_unittest("zero_width_face_pole", test_zero_width_face), &
                  new_unittest("ice_crossing_fold", test_crossing), &
                  new_unittest("fold_line_dofs_consistent", test_fold_line), &
                  new_unittest("solid_body_rotation_sphere", test_solid_body_sphere), &
                  new_unittest("free_drift_sphere", test_free_drift_sphere) &
                  ]
   end subroutine collect_cice_evp_c_fold_tests

   ! ------------------------------------------------------------------
   ! fixture plumbing
   ! ------------------------------------------------------------------

   subroutine ft_alloc(ft, ni, nj, nyp, ndte)
      !! Allocate and zero a fixture of `ni x nyp` owned cells with `NG` ghost
      !! rings; `nj` is the number of rows below the fold line (== `nyp` for a
      !! folded domain, `nyp / 2` for its doubled twin).
      type(ft_t), intent(inout) :: ft
      integer, intent(in) :: ni, nj, nyp, ndte
      integer :: nx, ny

      ft%ni = ni
      ft%nj = nj
      ft%nyp = nyp
      nx = ni + 2*NG
      ny = nyp + 2*NG
      ft%nx = nx
      ft%ny = ny
      ft%is = NG + 1
      ft%ie = NG + ni
      ft%js = NG + 1
      ft%je = NG + nyp
      allocate (ft%dxE(nx, ny), ft%dyE(nx, ny), ft%dxN(nx, ny), ft%dyN(nx, ny), ft%dxU(nx, ny), ft%dyU(nx, ny))
      allocate (ft%dxT(nx, ny), ft%dyT(nx, ny), ft%uarea(nx, ny), ft%tarea(nx, ny), ft%tarear(nx, ny))
      allocate (ft%earea(nx, ny), ft%narea(nx, ny), ft%earear(nx, ny), ft%narear(nx, ny))
      allocate (ft%rdxN(nx, ny), ft%rdxNr(nx, ny), ft%rdyE(nx, ny), ft%rdyEr(nx, ny))
      allocate (ft%epm(nx, ny), ft%npm(nx, ny), ft%uvm(nx, ny), ft%hm(nx, ny))
      allocate (ft%aice(nx, ny), ft%vice(nx, ny), ft%vsno(nx, ny), ft%aice_init(nx, ny), ft%aice0(nx, ny))
      allocate (ft%aicen(nx, ny, 1), ft%vicen(nx, ny, 1))
      allocate (ft%fcorE(nx, ny), ft%fcorN(nx, ny), ft%strairxE(nx, ny), ft%strairyN(nx, ny))
      allocate (ft%ss_tltxE(nx, ny), ft%ss_tltyN(nx, ny))
      allocate (ft%uocnE(nx, ny), ft%vocnE(nx, ny), ft%uocnN(nx, ny), ft%vocnN(nx, ny))
      allocate (ft%u(nx + 1, ny), ft%v(nx, ny + 1), ft%ice_err(nx, ny))
      ft%dxE = 1.0_wp; ft%dxN = 1.0_wp; ft%dxU = 1.0_wp; ft%dxT = 1.0_wp
      ft%dyE = 1.0_wp; ft%dyN = 1.0_wp; ft%dyU = 1.0_wp; ft%dyT = 1.0_wp
      ft%uarea = 1.0_wp; ft%tarea = 1.0_wp; ft%earea = 1.0_wp; ft%narea = 1.0_wp
      ft%tarear = 1.0_wp; ft%earear = 1.0_wp; ft%narear = 1.0_wp
      ft%rdxN = -1.0_wp; ft%rdxNr = -1.0_wp; ft%rdyE = -1.0_wp; ft%rdyEr = -1.0_wp
      ft%epm = 0.0_wp; ft%npm = 0.0_wp; ft%uvm = 0.0_wp; ft%hm = 1.0_wp
      ft%aice = 0.0_wp; ft%vice = 0.0_wp; ft%vsno = 0.0_wp; ft%aice_init = 0.0_wp; ft%aice0 = 1.0_wp
      ft%aicen = 0.0_wp; ft%vicen = 0.0_wp
      ft%fcorE = 0.0_wp; ft%fcorN = 0.0_wp; ft%strairxE = 0.0_wp; ft%strairyN = 0.0_wp
      ft%ss_tltxE = 0.0_wp; ft%ss_tltyN = 0.0_wp
      ft%uocnE = 0.0_wp; ft%vocnE = 0.0_wp; ft%uocnN = 0.0_wp; ft%vocnN = 0.0_wp
      ft%u = 0.0_wp; ft%v = 0.0_wp; ft%ice_err = 0_int32
      ft%p = ice_params_default(ICE_CONSTANTS_ICEPACK, 1, 4, 1)
      ft%p%ridge%kstrength = 0
      ft%p%dyn%ndte = ndte
      ft%dt = 900.0_wp
      ft%k = ice_evp_c_constants(ft%p, ft%dt)
      call ft%w%init(nx, ny)
      ft%w%iceEmask = .true.
      ft%w%iceNmask = .true.
   end subroutine ft_alloc

   subroutine ft_masks(ft)
      !! `epm/npm/uvm` from `hm` (a pure product of 0/1 values: exactly rotation-symmetric).
      type(ft_t), intent(inout) :: ft
      integer :: i, j, nx, ny
      nx = ft%nx
      ny = ft%ny
      ft%epm = 0.0_wp; ft%npm = 0.0_wp; ft%uvm = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            ft%epm(i, j) = ft%hm(i, j)*ft%hm(i + 1, j)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            ft%npm(i, j) = ft%hm(i, j)*ft%hm(i, j + 1)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx - 1
            ft%uvm(i, j) = ft%hm(i, j)*ft%hm(i + 1, j)*ft%hm(i, j + 1)*ft%hm(i + 1, j + 1)
         end do
      end do
   end subroutine ft_masks

   subroutine ft_ratios(ft)
      !! `ratiodxN/ratiodyE` and their reciprocals from the face lengths (as the
      !! production metrics build them; only land faces read them).
      type(ft_t), intent(inout) :: ft
      integer :: i, j
      ft%rdxN = 0.0_wp; ft%rdxNr = 0.0_wp; ft%rdyE = 0.0_wp; ft%rdyEr = 0.0_wp
      do j = 1, ft%ny
         do i = 1, ft%nx - 1
            if (ft%dxN(i, j) /= 0.0_wp) ft%rdxN(i, j) = -ft%dxN(i + 1, j)/ft%dxN(i, j)
            if (ft%rdxN(i, j) /= 0.0_wp) ft%rdxNr(i, j) = 1.0_wp/ft%rdxN(i, j)
         end do
      end do
      do j = 1, ft%ny - 1
         do i = 1, ft%nx
            if (ft%dyE(i, j) /= 0.0_wp) ft%rdyE(i, j) = -ft%dyE(i, j + 1)/ft%dyE(i, j)
            if (ft%rdyE(i, j) /= 0.0_wp) ft%rdyEr(i, j) = 1.0_wp/ft%rdyE(i, j)
         end do
      end do
   end subroutine ft_ratios

   subroutine ft_map(ft)
      type(ft_t), intent(inout) :: ft
      call ft%w%enter_data()
      !$acc enter data copyin(ft%dxE, ft%dyE, ft%dxN, ft%dyN, ft%dxU, ft%dyU, ft%dxT, ft%dyT)
      !$acc enter data copyin(ft%uarea, ft%tarea, ft%tarear, ft%earea, ft%narea, ft%earear, ft%narear)
      !$acc enter data copyin(ft%rdxN, ft%rdxNr, ft%rdyE, ft%rdyEr, ft%epm, ft%npm, ft%uvm, ft%hm)
      !$acc enter data copyin(ft%aice, ft%vice, ft%vsno, ft%aice_init, ft%aice0, ft%aicen, ft%vicen)
      !$acc enter data copyin(ft%fcorE, ft%fcorN, ft%strairxE, ft%strairyN, ft%ss_tltxE, ft%ss_tltyN)
      !$acc enter data copyin(ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN, ft%u, ft%v, ft%ice_err)
      !$acc enter data copyin(ft%p)
      ft%mapped = .true.
   end subroutine ft_map

   subroutine ft_unmap(ft)
      type(ft_t), intent(inout) :: ft
      !$acc exit data delete(ft%dxE, ft%dyE, ft%dxN, ft%dyN, ft%dxU, ft%dyU, ft%dxT, ft%dyT)
      !$acc exit data delete(ft%uarea, ft%tarea, ft%tarear, ft%earea, ft%narea, ft%earear, ft%narear)
      !$acc exit data delete(ft%rdxN, ft%rdxNr, ft%rdyE, ft%rdyEr, ft%epm, ft%npm, ft%uvm, ft%hm)
      !$acc exit data delete(ft%aice, ft%vice, ft%vsno, ft%aice_init, ft%aice0, ft%aicen, ft%vicen)
      !$acc exit data delete(ft%fcorE, ft%fcorN, ft%strairxE, ft%strairyN, ft%ss_tltxE, ft%ss_tltyN)
      !$acc exit data delete(ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN, ft%u, ft%v, ft%ice_err)
      !$acc exit data delete(ft%p)
      call ft%w%exit_data()
      ft%mapped = .false.
   end subroutine ft_unmap

   subroutine ft_pull(ft)
      !! Device -> host: everything a test inspects.
      type(ft_t), intent(inout) :: ft
      !$acc update self(ft%u, ft%v, ft%ice_err)
      !$acc update self(ft%w%stresspT, ft%w%stressmT, ft%w%zetax2T, ft%w%etax2T, ft%w%stress12U)
      !$acc update self(ft%w%strength, ft%w%aiE, ft%w%aiN, ft%w%emassdti, ft%w%nmassdti)
      !$acc update self(ft%w%strintxE, ft%w%strintyN, ft%w%divu, ft%w%shear, ft%w%vort)
      !$acc update self(ft%w%shearU, ft%w%uvelN, ft%w%vvelE, ft%w%iceEmask, ft%w%iceNmask)
   end subroutine ft_pull

   subroutine ft_pull_pair2(ft)
      !! Device -> host: the ping-pong partner pair.
      type(ft_t), intent(inout) :: ft
      !$acc update self(ft%w%u2, ft%w%v2)
   end subroutine ft_pull_pair2

   subroutine ft_destroy(ft)
      type(ft_t), intent(inout) :: ft
      if (ft%mapped) call ft_unmap(ft)
      call ft%w%destroy()
      if (allocated(ft%dxE)) deallocate (ft%dxE, ft%dyE, ft%dxN, ft%dyN, ft%dxU, ft%dyU, ft%dxT, ft%dyT)
      if (allocated(ft%uarea)) deallocate (ft%uarea, ft%tarea, ft%tarear, ft%earea, ft%narea, ft%earear, ft%narear)
      if (allocated(ft%rdxN)) deallocate (ft%rdxN, ft%rdxNr, ft%rdyE, ft%rdyEr, ft%epm, ft%npm, ft%uvm, ft%hm)
      if (allocated(ft%aice)) deallocate (ft%aice, ft%vice, ft%vsno, ft%aice_init, ft%aice0, ft%aicen, ft%vicen)
      if (allocated(ft%fcorE)) deallocate (ft%fcorE, ft%fcorN, ft%strairxE, ft%strairyN, ft%ss_tltxE, ft%ss_tltyN)
      if (allocated(ft%uocnE)) deallocate (ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN, ft%u, ft%v, ft%ice_err)
   end subroutine ft_destroy

   subroutine ft_step(ft, fold)
      !! One full step through the public driver (periodic in x, walls or fold in y).
      type(ft_t), intent(inout) :: ft
      logical, intent(in) :: fold
      call ice_evp_c_step(ft%w, ft%k, ft%nx, ft%ny, 1, ft%is, ft%ie, ft%js, ft%je, ft%p%ridge, ft%p%c, .true., .false., &
                          ft%dxE, ft%dyE, ft%dxN, ft%dyN, ft%dxU, ft%dyU, ft%dxT, ft%dyT, ft%uarea, ft%tarea, &
                          ft%tarear, ft%earea, ft%narea, ft%earear, ft%narear, ft%rdxN, ft%rdxNr, ft%rdyE, &
                          ft%rdyEr, ft%epm, ft%npm, ft%uvm, ft%hm, ft%aice, ft%vice, ft%vsno, ft%aice_init, &
                          ft%aice0, ft%aicen, ft%vicen, ft%fcorE, ft%fcorN, ft%strairxE, ft%strairyN, &
                          ft%ss_tltxE, ft%ss_tltyN, ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN, ft%u, ft%v, ft%ice_err, fold)
   end subroutine ft_step

   subroutine ft_prep(ft)
      type(ft_t), intent(inout) :: ft
      call ice_evp_c_prep(ft%w, ft%k, ft%nx, ft%ny, 1, ft%is, ft%ie, ft%js, ft%je, ft%p%ridge, ft%p%c, ft%hm, &
                          ft%aice, ft%vice, ft%vsno, ft%aice_init, ft%aice0, ft%aicen, ft%vicen, ft%tarea, &
                          ft%epm, ft%npm, ft%fcorE, ft%fcorN, ft%strairxE, ft%strairyN, ft%ss_tltxE, &
                          ft%ss_tltyN, ft%uocnE, ft%vocnN, ft%u, ft%v, ft%ice_err)
   end subroutine ft_prep

   subroutine ft_sub(ft, uo, vo, un, vn)
      type(ft_t), intent(inout) :: ft
      real(wp), intent(in) :: uo(ft%nx + 1, ft%ny), vo(ft%nx, ft%ny + 1)
      real(wp), intent(inout) :: un(ft%nx + 1, ft%ny), vn(ft%nx, ft%ny + 1)
      call ice_evp_c_substep(ft%w, ft%k, ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, uo, vo, un, vn, &
                             ft%dxE, ft%dyE, ft%dxN, ft%dyN, ft%dxU, ft%dyU, ft%dxT, ft%dyT, ft%uarea, ft%tarea, &
                             ft%earea, ft%narea, ft%earear, ft%narear, ft%rdxN, ft%rdxNr, ft%rdyE, ft%rdyEr, &
                             ft%epm, ft%npm, ft%uvm, ft%hm, ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN)
   end subroutine ft_sub

   subroutine ft_k123(ft)
      !! The three stress kernels of one subcycle on the fixture's own `(u, v)`.
      type(ft_t), intent(inout) :: ft
      call evp_strain_u_impl(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%k%e_factor, ft%k%avg_strength, ft%u, &
                             ft%v, ft%earea, ft%narea, ft%dxE, ft%dyN, ft%dxU, ft%dyU, &
                             ft%rdxN, ft%rdxNr, ft%rdyE, ft%rdyEr, ft%epm, ft%npm, ft%uvm, ft%w%imask, &
                             ft%w%shearU, ft%w%deltaU)
      call evp_stress_t_impl(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%k%arlx1i, ft%k%revp, ft%k%denom1, &
                             ft%k%epp2i, ft%k%e_factor, ft%k%Ktens, ft%k%deltamin, ft%k%cap_max, ft%u, ft%v, ft%dxN, &
                             ft%dyE, ft%dxT, ft%dyT, ft%uarea, ft%tarea, ft%w%strength, ft%w%shearU, ft%w%imask, &
                             ft%w%zetax2T, ft%w%etax2T, ft%w%stresspT, ft%w%stressmT)
      call evp_stress_u_impl(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%k%arlx1i, ft%k%revp, ft%k%denom1, &
                             ft%k%epp2i, ft%k%Ktens, ft%k%deltamin, ft%k%cap_max, ft%k%avg_strength, ft%hm, &
                             ft%tarea, ft%uarea, ft%w%etax2T, ft%w%strengthU, ft%w%deltaU, ft%w%shearU, &
                             ft%w%imask, ft%w%stress12U)
   end subroutine ft_k123

   ! ------------------------------------------------------------------
   ! the rotation (host) and the analytic fields
   ! ------------------------------------------------------------------

   subroutine rot_fill(a, ni, isum, jsum, sgn, dup)
      !! Make `a` rotation-covariant: rows whose rotated row `jp = jsum - j` is
      !! smaller (and, on the fold-line row, the west half) become `sgn * a(rotated)`
      !! (a self-conjugate point of a vector becomes 0), the others keep their
      !! (formula) values; then the ghost columns are re-wrapped.  `isum`/`jsum`
      !! are the stagger's index sums (T: `2NG+ni+1`, `2NG+2nj+1`; the table is in
      !! `rot_T` ... below).  The map is the continuous rotation
      !! `(x, y) -> (ni - x, 2 nj - y)` written in storage indices, derived here
      !! independently of the fold kernel.  `dup`: the array carries the duplicate
      !! column `ie+1` of a wrapped face/corner stagger.
      real(wp), intent(inout) :: a(:, :)
      integer, intent(in) :: ni, isum, jsum
      real(wp), intent(in) :: sgn
      logical, intent(in) :: dup
      integer :: i, j, ip, jp, n2

      n2 = size(a, 2)
      do j = 1, n2
         jp = jsum - j
         if (jp < 1 .or. jp > n2) cycle
         do i = NG + 1, NG + ni
            ip = NG + 1 + modulo(isum - i - NG - 1, ni)
            if (jp < j .or. (jp == j .and. ip > i)) then
               a(i, j) = sgn*a(ip, jp)
            else if (jp == j .and. ip == i .and. sgn < 0.0_wp) then
               a(i, j) = 0.0_wp
            end if
         end do
      end do
      call wrap_cols(a, ni, dup)
   end subroutine rot_fill

   subroutine wrap_cols(a, ni, dup)
      !! Periodic ghost columns from the interior ones (`dup`: the stagger stores the
      !! duplicate column `ie+1`, a copy of `is`).
      real(wp), intent(inout) :: a(:, :)
      integer, intent(in) :: ni
      logical, intent(in) :: dup
      integer :: j, kk
      do j = 1, size(a, 2)
         do kk = 1, NG
            a(NG + 1 - kk, j) = a(NG + ni + 1 - kk, j)
            if (dup) then
               a(NG + ni + 1 + kk, j) = a(NG + 1 + kk, j)
            else
               a(NG + ni + kk, j) = a(NG + kk, j)
            end if
         end do
         if (dup) a(NG + ni + 1, j) = a(NG + 1, j)
      end do
   end subroutine wrap_cols

   pure function wrapc(i, ni) result(iw)
      !! Interior column equivalent of storage column `i` (first interior column is `NG+1`).
      integer, intent(in) :: i, ni
      integer :: iw
      iw = NG + 1 + modulo(i - NG - 1, ni)
   end function wrapc

   pure function fld(fid, x, y, ni) result(f)
      !! Smooth non-symmetric analytic fields of position, periodic in x (period `ni`).
      integer, intent(in) :: fid, ni
      real(wp), intent(in) :: x, y
      real(wp) :: f, th
      th = TWO_PI*x/real(ni, wp)
      select case (fid)
      case (1)          ! concentration
         f = 0.62_wp + 0.28_wp*sin(th + 0.4_wp)*cos(0.55_wp*y)
      case (2)          ! thickness
         f = 1.1_wp + 0.4_wp*sin(2.0_wp*th + 1.0_wp)*sin(0.31_wp*y + 0.2_wp)
      case (3)          ! u
         f = 0.04_wp*sin(th + 0.3_wp)*cos(0.4_wp*y)
      case (4)          ! v
         f = 0.03_wp*cos(th)*sin(0.5_wp*y + 1.0_wp)
      case (5)          ! wind stress x
         f = 0.07_wp + 0.04_wp*sin(th + 2.0_wp)*cos(0.2_wp*y)
      case (6)          ! wind stress y
         f = 0.05_wp*cos(th + 0.7_wp)*sin(0.35_wp*y)
      case (7)          ! Coriolis parameter
         f = 1.3e-4_wp*(1.0_wp + 0.3_wp*sin(th)*cos(0.3_wp*y))
      case (8)          ! tilt x
         f = 2.0e-6_wp*sin(th + 1.1_wp)*cos(0.3_wp*y)
      case (9)          ! tilt y
         f = 2.0e-6_wp*cos(th + 0.2_wp)*sin(0.45_wp*y)
      case (10)         ! ocean u
         f = 0.05_wp*sin(th + 0.9_wp)*sin(0.25_wp*y + 0.3_wp)
      case (11)         ! ocean v
         f = 0.04_wp*cos(th + 1.9_wp)*cos(0.33_wp*y)
      case default
         f = 0.0_wp
      end select
   end function fld

   subroutine fill_formula(a, stag, fid, ni, scale)
      !! `a(i,j) = scale * fld(fid, x, y)` at the stagger `stag` ('T', 'E', 'N', 'U' ice-name;
      !! 'u', 'v' face storage), every column through its interior equivalent (so the
      !! ghost columns are the wrapped copy); rows above the fold line are overwritten by `rot_fill`.
      real(wp), intent(inout) :: a(:, :)
      character(len=1), intent(in) :: stag
      integer, intent(in) :: fid, ni
      real(wp), intent(in) :: scale
      integer :: i, j, iw, ip, jp
      real(wp) :: x, y

      do j = 1, size(a, 2)
         jp = j - NG
         do i = 1, size(a, 1)
            iw = wrapc(i, ni)
            ip = iw - NG
            select case (stag)
            case ('T')
               x = real(ip, wp) - 0.5_wp; y = real(jp, wp) - 0.5_wp
            case ('E', 'u')
               x = real(ip, wp); y = real(jp, wp) - 0.5_wp
               if (stag == 'u') x = real(ip - 1, wp)
            case ('N')
               x = real(ip, wp) - 0.5_wp; y = real(jp, wp)
            case ('v')
               x = real(ip, wp) - 0.5_wp; y = real(jp - 1, wp)
            case default                  ! 'U'
               x = real(ip, wp); y = real(jp, wp)
            end select
            a(i, j) = scale*fld(fid, x, y, ni)
         end do
      end do
   end subroutine fill_formula

   function hxf(x, y, ni, nj) result(h)
      !! Zonal scale factor [m]: even about the fold line `y = nj`, periodic-even in x.
      real(wp), intent(in) :: x, y
      integer, intent(in) :: ni, nj
      real(wp) :: h
      h = 2.5e4_wp*(0.72_wp + 0.28_wp*cos(PI*(y - real(nj, wp))/(2.0_wp*real(nj, wp))))* &
          (1.0_wp + 0.08_wp*cos(TWO_PI*x/real(ni, wp)))
   end function hxf

   function hyf(x, y, ni, nj) result(h)
      !! Meridional scale factor [m], same symmetry.
      real(wp), intent(in) :: x, y
      integer, intent(in) :: ni, nj
      real(wp) :: h
      h = 2.0e4_wp*(0.85_wp + 0.15_wp*cos(PI*(y - real(nj, wp))/(2.0_wp*real(nj, wp))))* &
          (1.0_wp + 0.1_wp*cos(TWO_PI*x/real(ni, wp)))
   end function hyf

   subroutine metric_formula(a, stag, which, ni, nj)
      !! `which = 'x'|'y'`: scale factor at the stagger `stag` (as `fill_formula`).
      real(wp), intent(inout) :: a(:, :)
      character(len=1), intent(in) :: stag, which
      integer, intent(in) :: ni, nj
      integer :: i, j, iw, ip, jp
      real(wp) :: x, y

      do j = 1, size(a, 2)
         jp = j - NG
         do i = 1, size(a, 1)
            iw = wrapc(i, ni)
            ip = iw - NG
            select case (stag)
            case ('T')
               x = real(ip, wp) - 0.5_wp; y = real(jp, wp) - 0.5_wp
            case ('E')
               x = real(ip, wp); y = real(jp, wp) - 0.5_wp
            case ('N')
               x = real(ip, wp) - 0.5_wp; y = real(jp, wp)
            case default
               x = real(ip, wp); y = real(jp, wp)
            end select
            if (which == 'x') then
               a(i, j) = hxf(x, y, ni, nj)
            else
               a(i, j) = hyf(x, y, ni, nj)
            end if
         end do
      end do
   end subroutine metric_formula

   subroutine sym_geometry(ft, with_island)
      !! Rotation-covariant geometry and metrics about the fold line `y = nj`:
      !! analytic non-uniform scale factors (formula up to the fold line, exact
      !! rotated copy above), the walled south edge (and, in the doubled domain, the
      !! image wall in the north), an optional island next to the fold line.
      type(ft_t), intent(inout) :: ft
      logical, intent(in) :: with_island
      integer :: ni, nj, i, j, ip, jp
      ni = ft%ni
      nj = ft%nj
      ! T mask: wet except the south ghost band and the island; north by rotation.
      do j = 1, ft%ny
         jp = j - NG
         do i = 1, ft%nx
            ip = wrapc(i, ni) - NG
            ft%hm(i, j) = 1.0_wp
            if (jp < 1) ft%hm(i, j) = 0.0_wp
            if (with_island .and. (ip == 4 .or. ip == 5) .and. jp == nj - 1) ft%hm(i, j) = 0.0_wp
         end do
      end do
      call rot_fill(ft%hm, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call ft_masks(ft)
      ! Scale factors and the derived areas, one rotation sum per stagger.
      call metric_formula(ft%dxT, 'T', 'x', ni, nj); call rot_fill(ft%dxT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call metric_formula(ft%dyT, 'T', 'y', ni, nj); call rot_fill(ft%dyT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call metric_formula(ft%dxE, 'E', 'x', ni, nj); call rot_fill(ft%dxE, ni, 2*NG + ni, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call metric_formula(ft%dyE, 'E', 'y', ni, nj); call rot_fill(ft%dyE, ni, 2*NG + ni, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call metric_formula(ft%dxN, 'N', 'x', ni, nj); call rot_fill(ft%dxN, ni, 2*NG + ni + 1, 2*NG + 2*nj, 1.0_wp, .false.)
      call metric_formula(ft%dyN, 'N', 'y', ni, nj); call rot_fill(ft%dyN, ni, 2*NG + ni + 1, 2*NG + 2*nj, 1.0_wp, .false.)
      call metric_formula(ft%dxU, 'U', 'x', ni, nj); call rot_fill(ft%dxU, ni, 2*NG + ni, 2*NG + 2*nj, 1.0_wp, .false.)
      call metric_formula(ft%dyU, 'U', 'y', ni, nj); call rot_fill(ft%dyU, ni, 2*NG + ni, 2*NG + 2*nj, 1.0_wp, .false.)
      ft%tarea = ft%dxT*ft%dyT
      ft%earea = ft%dxE*ft%dyE
      ft%narea = ft%dxN*ft%dyN
      ft%uarea = ft%dxU*ft%dyU
      ft%tarear = 1.0_wp/ft%tarea
      ft%earear = 1.0_wp/ft%earea
      ft%narear = 1.0_wp/ft%narea
      call ft_ratios(ft)
   end subroutine sym_geometry

   subroutine sym_state(ft, patch)
      !! Rotation-covariant ice state, velocity and forcing (non-symmetric formulas
      !! below the fold line, their exact images above).  `patch`: ice only in a patch
      !! straddling the fold (columns 4..13, the last three rows) with a vortex wind
      !! about `(ni/2, nj)`; else ice everywhere wet and the generic forcing.
      type(ft_t), intent(inout) :: ft
      logical, intent(in) :: patch
      integer :: ni, nj, i, j, ip, jp
      real(wp) :: x, y, tau0, rscale

      ni = ft%ni
      nj = ft%nj
      ! ---- T state ----
      call fill_formula(ft%aice, 'T', 1, ni, 1.0_wp)
      call fill_formula(ft%vice, 'T', 2, ni, 1.0_wp)
      do j = 1, ft%ny
         jp = j - NG
         do i = 1, ft%nx
            ip = wrapc(i, ni) - NG
            ft%vice(i, j) = ft%aice(i, j)*ft%vice(i, j)
            if (patch) then
               if (.not. (ip >= 4 .and. ip <= 13 .and. jp >= nj - 2)) then
                  ft%aice(i, j) = 0.0_wp
                  ft%vice(i, j) = 0.0_wp
               else
                  ft%aice(i, j) = 0.9_wp
                  ft%vice(i, j) = 0.9_wp*1.5_wp
               end if
            end if
            ft%aice(i, j) = ft%aice(i, j)*ft%hm(i, j)
            ft%vice(i, j) = ft%vice(i, j)*ft%hm(i, j)
         end do
      end do
      call rot_fill(ft%aice, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call rot_fill(ft%vice, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      ft%vsno = 0.0_wp
      ft%aice_init = ft%aice
      ft%aice0 = 1.0_wp - ft%aice
      ft%aicen(:, :, 1) = ft%aice
      ft%vicen(:, :, 1) = ft%vice
      ! ---- the scalar Coriolis parameter (copy under rotation) ----
      call fill_formula(ft%fcorE, 'E', 7, ni, 1.0_wp); call rot_fill(ft%fcorE, ni, 2*NG + ni, 2*NG + 2*nj + 1, 1.0_wp, .false.)
      call fill_formula(ft%fcorN, 'N', 7, ni, 1.0_wp); call rot_fill(ft%fcorN, ni, 2*NG + ni + 1, 2*NG + 2*nj, 1.0_wp, .false.)
      ! ---- vectors at E (x components) and N (y components): negate under rotation ----
      if (patch) then
         tau0 = 0.12_wp
         rscale = 6.0_wp
         do j = 1, ft%ny
            do i = 1, ft%nx
               ip = wrapc(i, ni) - NG
               jp = j - NG
               ! vortex stress about (ni/2, nj): tau = tau0 (-(y - nj), (x - ni/2)) / rscale
               x = real(ip, wp); y = real(jp, wp) - 0.5_wp
               ft%strairxE(i, j) = -tau0*(y - real(nj, wp))/rscale
               x = real(ip, wp) - 0.5_wp; y = real(jp, wp)
               ft%strairyN(i, j) = tau0*(x - 0.5_wp*real(ni, wp))/rscale
            end do
         end do
         ft%ss_tltxE = 0.0_wp
         ft%ss_tltyN = 0.0_wp
         ft%uocnE = 0.0_wp; ft%vocnE = 0.0_wp; ft%uocnN = 0.0_wp; ft%vocnN = 0.0_wp
      else
         call fill_formula(ft%strairxE, 'E', 5, ni, 1.0_wp)
         call fill_formula(ft%strairyN, 'N', 6, ni, 1.0_wp)
         call fill_formula(ft%ss_tltxE, 'E', 8, ni, 1.0_wp)
         call fill_formula(ft%ss_tltyN, 'N', 9, ni, 1.0_wp)
         call fill_formula(ft%uocnE, 'E', 10, ni, 1.0_wp)
         call fill_formula(ft%vocnE, 'E', 11, ni, 1.0_wp)
         call fill_formula(ft%uocnN, 'N', 10, ni, 1.0_wp)
         call fill_formula(ft%vocnN, 'N', 11, ni, 1.0_wp)
      end if
      call mask_E(ft%strairxE, ft); call mask_E(ft%ss_tltxE, ft); call mask_E(ft%uocnE, ft); call mask_E(ft%vocnE, ft)
      call mask_N(ft%strairyN, ft); call mask_N(ft%ss_tltyN, ft); call mask_N(ft%uocnN, ft); call mask_N(ft%vocnN, ft)
      call rot_fill(ft%strairxE, ni, 2*NG + ni, 2*NG + 2*nj + 1, -1.0_wp, .false.)
      call rot_fill(ft%ss_tltxE, ni, 2*NG + ni, 2*NG + 2*nj + 1, -1.0_wp, .false.)
      call rot_fill(ft%uocnE, ni, 2*NG + ni, 2*NG + 2*nj + 1, -1.0_wp, .false.)
      call rot_fill(ft%vocnE, ni, 2*NG + ni, 2*NG + 2*nj + 1, -1.0_wp, .false.)
      call rot_fill(ft%strairyN, ni, 2*NG + ni + 1, 2*NG + 2*nj, -1.0_wp, .false.)
      call rot_fill(ft%ss_tltyN, ni, 2*NG + ni + 1, 2*NG + 2*nj, -1.0_wp, .false.)
      call rot_fill(ft%uocnN, ni, 2*NG + ni + 1, 2*NG + 2*nj, -1.0_wp, .false.)
      call rot_fill(ft%vocnN, ni, 2*NG + ni + 1, 2*NG + 2*nj, -1.0_wp, .false.)
      ! ---- the face velocity pair ----
      call fill_formula(ft%u, 'u', 3, ni, 1.0_wp)
      call fill_formula(ft%v, 'v', 4, ni, 1.0_wp)
      if (patch) then
         ft%u = 0.0_wp
         ft%v = 0.0_wp
      end if
      do j = 1, ft%ny
         do i = 1, ft%nx
            ft%u(i + 1, j) = ft%u(i + 1, j)*ft%epm(i, j)
            ft%v(i, j + 1) = ft%v(i, j + 1)*ft%npm(i, j)
         end do
      end do
      call rot_fill(ft%u, ni, 2*NG + ni + 2, 2*NG + 2*nj + 1, -1.0_wp, .true.)
      call rot_fill(ft%v, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, .false.)
   end subroutine sym_state

   subroutine mask_E(a, ft)
      type(ft_t), intent(in) :: ft
      real(wp), intent(inout) :: a(:, :)
      a = a*ft%epm
   end subroutine mask_E

   subroutine mask_N(a, ft)
      type(ft_t), intent(in) :: ft
      real(wp), intent(inout) :: a(:, :)
      a = a*ft%npm
   end subroutine mask_N

   function maxabs_diff(a, b, i0, i1, j0, j1) result(d)
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: d
      integer :: i, j
      d = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            if (.not. (abs(a(i, j) - b(i, j)) <= d)) d = abs(a(i, j) - b(i, j))
         end do
      end do
   end function maxabs_diff

   function maxabs(a, i0, i1, j0, j1) result(d)
      real(wp), intent(in) :: a(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: d
      integer :: i, j
      d = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            if (.not. (abs(a(i, j)) <= d)) d = abs(a(i, j))
         end do
      end do
   end function maxabs

   function asym(a, ni, isum, jsum, sgn, i0, i1, j0, j1) result(d)
      !! Largest `|a(i,j) - sgn a(rot(i,j))|` over the window (rotation defect).
      real(wp), intent(in) :: a(:, :)
      integer, intent(in) :: ni, isum, jsum, i0, i1, j0, j1
      real(wp), intent(in) :: sgn
      real(wp) :: d
      integer :: i, j, ip, jp
      d = 0.0_wp
      do j = j0, j1
         jp = jsum - j
         if (jp < 1 .or. jp > size(a, 2)) cycle
         do i = i0, i1
            ip = NG + 1 + modulo(isum - i - NG - 1, ni)
            if (.not. (abs(a(i, j) - sgn*a(ip, jp)) <= d)) d = abs(a(i, j) - sgn*a(ip, jp))
         end do
      end do
   end function asym

   ! ------------------------------------------------------------------
   ! 1. the fold kernel against the host rotation
   ! ------------------------------------------------------------------

   subroutine test_fold_kernel(error)
      type(error_type), allocatable, intent(out) :: error
      call fold_kernel_case(error, 12)
      if (allocated(error)) return
      ! an odd column count has a self-conjugate v column on the fold row (zeroed): the branch only runs here
      call fold_kernel_case(error, 13)
   end subroutine test_fold_kernel

   subroutine fold_kernel_case(error, ni)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ni
      integer, parameter :: nj = 5
      integer :: nx, ny, is, ie, js, je, i, j, jf, seed
      real(wp), allocatable :: u(:, :), v(:, :), ur(:, :), vr(:, :), u0(:, :), v0(:, :)
      real(wp) :: dmax, rnd

      nx = ni + 2*NG
      ny = nj + 2*NG
      is = NG + 1; ie = NG + ni; js = NG + 1; je = NG + nj
      jf = je + 1
      allocate (u(nx + 1, ny), v(nx, ny + 1), ur(nx + 1, ny), vr(nx, ny + 1), u0(nx + 1, ny), v0(nx, ny + 1))
      seed = 12345
      do j = 1, ny
         do i = 1, nx + 1
            seed = mod(seed*1103515245 + 12345, 2147483647)
            rnd = real(seed, wp)/2147483647.0_wp
            u(i, j) = rnd - 0.5_wp
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            seed = mod(seed*1103515245 + 12345, 2147483647)
            rnd = real(seed, wp)/2147483647.0_wp
            v(i, j) = rnd - 0.5_wp
         end do
      end do
      ! periodic columns first (the kernel's contract)
      call wrap_cols(u, ni, .true.)
      call wrap_cols(v, ni, .false.)
      u0 = u
      v0 = v
      ! host reference: the rotation on the rows above the fold line
      ur = u
      vr = v
      call rot_fill(ur, ni, 2*NG + ni + 2, 2*NG + 2*nj + 1, -1.0_wp, .true.)
      call rot_fill(vr, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, .false.)
      !$acc enter data copyin(u, v)
      call ice_evp_c_fold_faces(nx, ny, is, ie, js, je, u, v)
      !$acc update self(u, v)
      dmax = max(maxabs_diff(u, ur, 1, nx, jf, ny), maxabs_diff(v, vr, 1, nx, jf, ny + 1))
      write (error_unit, '(a,es10.2)') "  [evp_c_fold] fold kernel vs host rotation (ni odd/even), max |diff| = ", dmax
      call check(error, dmax == 0.0_wp, "the fold kernel is bit-for-bit the host rotation")
      if (.not. allocated(error)) then
         ! rows below the fold line are untouched; the east half of the fold row too
         call check(error, maxabs_diff(u, u0, 1, nx + 1, 1, jf - 1) == 0.0_wp, "u below the fold line untouched")
      end if
      if (.not. allocated(error)) then
         call check(error, maxabs_diff(v, v0, 1, nx, 1, jf - 1) == 0.0_wp, "v below the fold row untouched")
      end if
      if (.not. allocated(error)) then
         ! fold-line row: v(i) = -v(i') exactly, every storage column
         call check(error, asym(v, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, 1, nx, jf, jf) == 0.0_wp, &
                    "the fold-line row of v is exactly antisymmetric")
      end if
      if (.not. allocated(error) .and. mod(ni, 2) == 1) then
         ! self-conjugate column (v = -v): zero, ghost copies included
         do i = 1, nx
            if (modulo(i - NG - 1, ni) + 1 == (ni + 1)/2) then
               if (v(i, jf) /= 0.0_wp) call check(error, .false., "the self-conjugate fold-row column is not zero")
            end if
         end do
      end if
      if (.not. allocated(error)) then
         ! idempotence
         u0 = u
         v0 = v
         !$acc update device(u, v)
         call ice_evp_c_fold_faces(nx, ny, is, ie, js, je, u, v)
         !$acc update self(u, v)
         call check(error, maxabs_diff(u, u0, 1, nx, 1, ny) == 0.0_wp .and. &
                    maxabs_diff(v, v0, 1, nx, 1, ny + 1) == 0.0_wp, "the fold is idempotent")
      end if
      !$acc exit data delete(u, v)
   end subroutine fold_kernel_case

   ! ------------------------------------------------------------------
   ! 2. rotation symmetry of the doubled domain
   ! ------------------------------------------------------------------

   subroutine set_variant(ft, variant)
      !! Rheology variants: 1 classic EVP, `avg_zeta`; 2 revised EVP with a tensile
      !! strength and an odd subcycle count (the ping-pong ends on the partner pair);
      !! 3 `avg_strength` viscosity with the summed (Kreyscher) cap and a larger `Ktens`.
      type(ft_t), intent(inout) :: ft
      integer, intent(in) :: variant
      select case (variant)
      case (2)
         ft%p%dyn%ndte = 7
         ft%p%dyn%revised_evp = .true.
         ft%p%dyn%Ktens = 0.2_wp
      case (3)
         ft%p%dyn%visc_method = ICE_VISC_AVG_STRENGTH
         ft%p%dyn%capping_max = .false.
         ft%p%dyn%Ktens = 0.3_wp
      case default
         continue
      end select
      ft%k = ice_evp_c_constants(ft%p, ft%dt)
   end subroutine set_variant

   subroutine make_pair(fdom, ddom, ni, nj, ndte, patch, variant)
      !! The folded fixture and its doubled twin on the same formulas.
      type(ft_t), intent(inout) :: fdom, ddom
      integer, intent(in) :: ni, nj, ndte, variant
      logical, intent(in) :: patch
      call ft_alloc(fdom, ni, nj, nj, ndte)
      call set_variant(fdom, variant)
      call sym_geometry(fdom, .true.)
      call sym_state(fdom, patch)
      call ft_alloc(ddom, ni, nj, 2*nj, ndte)
      call set_variant(ddom, variant)
      call sym_geometry(ddom, .true.)
      call sym_state(ddom, patch)
   end subroutine make_pair

   subroutine doubled_defects(ddom, du, dv, dsp, dsm, ds12, scale_u, scale_v, scale_s)
      !! Rotation defects of the doubled state over the owned rows.
      type(ft_t), intent(in) :: ddom
      real(wp), intent(out) :: du, dv, dsp, dsm, ds12, scale_u, scale_v, scale_s
      integer :: ni, nj
      ni = ddom%ni
      nj = ddom%nj
      du = asym(ddom%u, ni, 2*NG + ni + 2, 2*NG + 2*nj + 1, -1.0_wp, ddom%is, ddom%ie, ddom%js, ddom%je)
      dv = asym(ddom%v, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, ddom%is, ddom%ie, ddom%js + 1, ddom%je + 1)
      dsp = asym(ddom%w%stresspT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, ddom%is, ddom%ie, ddom%js, ddom%je)
      dsm = asym(ddom%w%stressmT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, ddom%is, ddom%ie, ddom%js, ddom%je)
      ds12 = asym(ddom%w%stress12U, ni, 2*NG + ni + 2, 2*NG + 2*nj + 2, 1.0_wp, ddom%is + 1, ddom%ie + 1, ddom%js, ddom%je + 1)
      scale_u = maxabs(ddom%u, ddom%is, ddom%ie + 1, ddom%js, ddom%je)
      scale_v = maxabs(ddom%v, ddom%is, ddom%ie, ddom%js + 1, ddom%je + 1)
      scale_s = max(maxabs(ddom%w%stresspT, ddom%is, ddom%ie, ddom%js, ddom%je), &
                    maxabs(ddom%w%stressmT, ddom%is, ddom%ie, ddom%js, ddom%je))
   end subroutine doubled_defects

   subroutine test_doubled_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      type(ft_t) :: fdom, ddom
      integer :: n
      real(wp) :: du, dv, dsp, dsm, ds12, su, sv, ss, rel

      call make_pair(fdom, ddom, 16, 6, 8, .false., 1)
      call ft_destroy(fdom)
      call ft_map(ddom)
      do n = 1, 4
         call ft_step(ddom, .false.)
         call ft_pull(ddom)
         call doubled_defects(ddom, du, dv, dsp, dsm, ds12, su, sv, ss)
         rel = max(du/su, dv/sv, dsp/ss, dsm/ss, ds12/ss)
         write (error_unit, '(a,i2,a,5es10.2,a,3es10.2)') "  [evp_c_fold] doubled step ", n, ": rotation defect u v s+ s- s12 = ", &
            du, dv, dsp, dsm, ds12, "   scales u v s = ", su, sv, ss
         call check(error, rel < 1.0e-12_wp, "the doubled-domain state stays rotation-invariant to round-off")
         if (allocated(error)) exit
      end do
      call ft_destroy(ddom)
   end subroutine test_doubled_symmetry

   ! ------------------------------------------------------------------
   ! 3. folded == lower half of the doubled domain
   ! ------------------------------------------------------------------

   subroutine test_fold_equals_doubled(error)
      type(error_type), allocatable, intent(out) :: error
      type(ft_t) :: fdom, ddom
      integer :: n, is, ie, js, je, variant
      real(wp) :: du, dv, dp, dm, d12, dstr, ddiv, dshear, su, sv, ss, rel
      character(len=24), parameter :: tag(3) = [character(len=24) :: "classic", "revised, odd ndte, Ktens", "avg_strength, sum cap"]

      do variant = 1, 3
         call make_pair(fdom, ddom, 16, 6, 8, .false., variant)
         call ft_map(fdom)
         call ft_map(ddom)
         is = fdom%is; ie = fdom%ie; js = fdom%js; je = fdom%je
         do n = 1, 4
            call ft_step(fdom, .true.)
            call ft_step(ddom, .false.)
            call ft_pull(fdom)
            call ft_pull(ddom)
            du = maxabs_diff(fdom%u, ddom%u, is, ie + 1, js, je)
            dv = maxabs_diff(fdom%v, ddom%v, is, ie, js + 1, je + 1)
            dp = maxabs_diff(fdom%w%stresspT, ddom%w%stresspT, is, ie, js, je)
            dm = maxabs_diff(fdom%w%stressmT, ddom%w%stressmT, is, ie, js, je)
            d12 = maxabs_diff(fdom%w%stress12U, ddom%w%stress12U, is, ie + 1, js, je + 1)
            dstr = maxabs_diff(fdom%w%strength, ddom%w%strength, is, ie, js, je)
            ddiv = maxabs_diff(fdom%w%divu, ddom%w%divu, is, ie, js, je)
            dshear = maxabs_diff(fdom%w%shear, ddom%w%shear, is, ie, js, je)
            su = maxabs(fdom%u, is, ie + 1, js, je)
            sv = maxabs(fdom%v, is, ie, js + 1, je + 1)
            ss = max(maxabs(fdom%w%stresspT, is, ie, js, je), maxabs(fdom%w%stressmT, is, ie, js, je))
        write (error_unit, '(a,a,a,i2,a,3es10.2,a,3es10.2)') "  [evp_c_fold] fold vs doubled (", trim(tag(variant)), ") step ", n, &
               ": |du| |dv| |ds12| = ", du, dv, d12, "   |dsp| |dsm| |dstrength| = ", dp, dm, dstr
      write (error_unit, '(a,2es10.2,a,3es10.2)') "      |d divu| |d shear| = ", ddiv, dshear, "   scales u v stress = ", su, sv, ss
            rel = max(du/su, dv/sv, dp/ss, dm/ss, d12/ss)
            call check(error, dstr == 0.0_wp, "the strength (a function of the T state only) is identical")
            if (.not. allocated(error)) call check(error, rel < 1.0e-12_wp, &
                                                   "the folded run equals the lower half of the doubled run to round-off")
            if (allocated(error)) exit
         end do
         call ft_destroy(fdom)
         call ft_destroy(ddom)
         if (allocated(error)) return
      end do
   end subroutine test_fold_equals_doubled

   ! ------------------------------------------------------------------
   ! 4. ice crossing the fold, subcycle by subcycle
   ! ------------------------------------------------------------------

   subroutine fold_line_report_arrays(ft, u, v, dv_row, du_ghost, dv_ghost, d12_slots, dsp_ghost, dsm_ghost, ss)
      type(ft_t), intent(in) :: ft
      real(wp), intent(in) :: u(:, :), v(:, :)
      real(wp), intent(out) :: dv_row, du_ghost, dv_ghost, d12_slots, dsp_ghost, dsm_ghost, ss
      integer :: ni, nj, nx, ny, jf
      ni = ft%ni
      nj = ft%nj
      nx = ft%nx
      ny = ft%ny
      jf = NG + nj + 1
      dv_row = asym(v, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, 1, nx, jf, jf)
      du_ghost = asym(u, ni, 2*NG + ni + 2, 2*NG + 2*nj + 1, -1.0_wp, 1, nx, jf, ny)
      dv_ghost = asym(v, ni, 2*NG + ni + 1, 2*NG + 2*nj + 2, -1.0_wp, 1, nx, jf + 1, ny + 1)
      ! corner stress, fold-line row of corners = storage row jf (NE corners of the last T row)
      d12_slots = asym(ft%w%stress12U, ni, 2*NG + ni + 2, 2*NG + 2*nj + 2, 1.0_wp, ft%is, ft%ie, jf, jf)
      ! T stress ring 1 above the fold line against the mirrored owned row
      dsp_ghost = asym(ft%w%stresspT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, ft%is, ft%ie, jf, jf)
      dsm_ghost = asym(ft%w%stressmT, ni, 2*NG + ni + 1, 2*NG + 2*nj + 1, 1.0_wp, ft%is, ft%ie, jf, jf)
      ss = max(maxabs(ft%w%stresspT, ft%is, ft%ie, ft%js, ft%je), maxabs(ft%w%stressmT, ft%is, ft%ie, ft%js, ft%je))
   end subroutine fold_line_report_arrays

   subroutine test_corner_formulas(error)
      !! The header table, from the continuous rotation `(x, y) -> (ni - x, 2 nj - y)`: a marked
      !! value stored at the corner of position `(x, y)` lands, after mirroring, at the storage
      !! index the table gives, in both corner storages, and the fold-line row is the row whose
      !! mirror is itself.  Positions: `stress12U(c, r)` at `x = c - NG - 1`, `y = r - NG - 1` (SW
      !! storage); `shearU(i, j)` at `x = i - NG`, `y = j - NG` (NE storage).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ni = 12, nj = 5
      integer :: c, r, cm, rm, nbad, nfold_s, nfold_u
      real(wp) :: x, y, xm, ym

      nbad = 0
      nfold_s = 0
      nfold_u = 0
      do r = NG + 1, NG + nj + 1 + NG
         do c = NG + 1, NG + ni
            x = real(c - NG - 1, wp); y = real(r - NG - 1, wp)
            xm = real(ni, wp) - x; ym = real(2*nj, wp) - y
            cm = nint(xm) + NG + 1; rm = nint(ym) + NG + 1
            if (cm /= 2*NG + ni + 2 - c .and. modulo(cm - NG - 1, ni) /= modulo(2*NG + ni + 2 - c - NG - 1, ni)) nbad = nbad + 1
            if (rm /= 2*NG + 2*nj + 2 - r) nbad = nbad + 1
            if (rm == r) then
               nfold_s = nfold_s + 1
               if (r /= NG + nj + 1) nbad = nbad + 1
            end if
            ! NE storage, same physical corner
            x = real(c - NG, wp); y = real(r - NG, wp)
            xm = real(ni, wp) - x; ym = real(2*nj, wp) - y
            cm = nint(xm) + NG; rm = nint(ym) + NG
            if (modulo(cm - NG, ni) /= modulo(2*NG + ni - c - NG, ni)) nbad = nbad + 1
            if (rm /= 2*NG + 2*nj - r) nbad = nbad + 1
            if (rm == r) then
               nfold_u = nfold_u + 1
               if (r /= NG + nj) nbad = nbad + 1
            end if
         end do
      end do
      call check(error, nbad == 0, "corner storage mirror formulas of the header table")
      if (.not. allocated(error)) call check(error, nfold_s == ni .and. nfold_u == ni, &
                                           "exactly one fold-line row per corner storage (stress12U row ng+nj+1, shearU row ng+nj)")
   end subroutine test_corner_formulas

   subroutine test_zero_width_face(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: variant
      do variant = 1, 3, 2           ! classic avg_zeta, and avg_strength (the corner viscosity divides by deltamin * uarea)
         call zero_width_case(error, variant)
         if (allocated(error)) return
      end do
   end subroutine test_zero_width_face

   subroutine zero_width_case(error, variant)
      !! The tripolar grids of this tree carry bipoles in open water: the faces there have zero
      !! length, the ocean closes them (`epm = 0`) and they carry no velocity.  The T strain rates
      !! divide the face velocity by the face length; the quotient of a closed zero-length face is
      !! 0/0 unless guarded (`evp_quot`).  One such face (and the N face above a second one) inside
      !! the folded domain: every velocity and stress must stay finite, the error word clear.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: variant
      type(ft_t) :: ft
      integer :: n, i0, j0, i1, j1, i2, j2
      logical :: fin

      call ft_alloc(ft, 16, 6, 6, 8)
      call set_variant(ft, variant)
      call sym_geometry(ft, .false.)
      call sym_state(ft, .false.)
      i0 = NG + 8
      j0 = NG + 3
      i1 = NG + 11
      j1 = NG + 4
      ft%dxE(i0, j0) = 0.0_wp
      ft%dyE(i0, j0) = 0.0_wp
      ft%epm(i0, j0) = 0.0_wp
      ft%u(i0 + 1, j0) = 0.0_wp
      ft%dxN(i1, j1) = 0.0_wp
      ft%dyN(i1, j1) = 0.0_wp
      ft%npm(i1, j1) = 0.0_wp
      ft%v(i1, j1 + 1) = 0.0_wp
      ! a zero-area corner (a pole vertex) with its ice mask closed, as `ice_cgrid_metrics_t%init` leaves it
      i2 = NG + 5
      j2 = NG + 2
      ft%uarea(i2, j2) = 0.0_wp
      ft%dxU(i2, j2) = 0.0_wp
      ft%dyU(i2, j2) = 0.0_wp
      ft%uvm(i2, j2) = 0.0_wp
      call ft_ratios(ft)
      call ft_map(ft)
      do n = 1, 3
         call ft_step(ft, .true.)
      end do
      call ft_pull(ft)
      fin = all(ieee_is_finite(ft%u)) .and. all(ieee_is_finite(ft%v)) .and. all(ieee_is_finite(ft%w%stresspT)) .and. &
            all(ieee_is_finite(ft%w%stressmT)) .and. all(ieee_is_finite(ft%w%stress12U))
      call check(error, fin, "velocities and stresses stay finite next to a zero-length closed face")
      if (.not. allocated(error)) call check(error, all(ft%ice_err == 0_int32), "no error bit next to a zero-length closed face")
      if (.not. allocated(error)) call check(error, ft%u(i0 + 1, j0) == 0.0_wp .and. ft%v(i1, j1 + 1) == 0.0_wp, &
                                             "the closed zero-length faces carry no velocity")
      call ft_destroy(ft)
   end subroutine zero_width_case

   subroutine test_crossing(error)
      type(error_type), allocatable, intent(out) :: error
      type(ft_t) :: ft
      integer :: n, i, ni, nj, je
      real(wp) :: dv_row, du_g, dv_g, d12, dsp, dsm, ss, flux, flux_abs, vmax
      real(wp) :: stored(3)

      ni = 16
      nj = 6
      je = NG + nj
      call ft_alloc(ft, ni, nj, nj, 10)
      call sym_geometry(ft, .false.)
      call sym_state(ft, .true.)
      call ft_map(ft)
      do n = 1, 5
         call ft_step(ft, .true.)
         call ft_pull(ft)
         call fold_line_report_arrays(ft, ft%u, ft%v, dv_row, du_g, dv_g, d12, dsp, dsm, ss)
         ! mass flux through the fold line: sum over the fold row of (mass/dt)*dt*v*dx telescopes pairwise
         flux = 0.0_wp
         flux_abs = 0.0_wp
         do i = ft%is, ft%ie
            flux = flux + ft%dxN(i, je)*ft%w%nmassdti(i, je)*ft%dt*ft%v(i, je + 1)
            flux_abs = flux_abs + abs(ft%dxN(i, je)*ft%w%nmassdti(i, je)*ft%dt*ft%v(i, je + 1))
         end do
         vmax = maxabs(ft%v, ft%is, ft%ie, je + 1, je + 1)
         write (error_unit, '(a,i2,a,es10.2,a,es10.2)') "  [evp_c_fold] crossing step ", n, ": max|v| on the fold line ", vmax, &
            "  net / gross mass flux through it ", abs(flux)/max(flux_abs, tiny(1.0_wp))
         write (error_unit, '(a,3es10.2)') "      defects: v row, u ghost, v ghost = ", dv_row, du_g, dv_g
         write (error_unit, '(a,3es10.2,a,es10.2)') "      ghost-ring |dsp| |dsm| and corner-slot |ds12| = ", dsp, dsm, d12, "  stress scale ", ss
         stored = [dv_row, du_g, dv_g]
         call check(error, all(stored == 0.0_wp), "ghost rows and the fold-line row are exact rotated images")
         if (.not. allocated(error)) call check(error, vmax > 1.0e-4_wp, "the ice does cross the fold line (v /= 0 on it)")
         if (.not. allocated(error)) call check(error, abs(flux) <= 1.0e-12_wp*flux_abs, &
                                                "the mass flux through the fold line telescopes")
         if (.not. allocated(error)) call check(error, max(dsp, dsm, d12) <= 1.0e-12_wp*ss, &
                                                "ghost-ring stresses agree with the mirrored owned values")
         if (allocated(error)) exit
      end do
      call ft_destroy(ft)
   end subroutine test_crossing

   ! ------------------------------------------------------------------
   ! 5. fold-line DOFs after every subcycle (manual subcycle loop)
   ! ------------------------------------------------------------------

   subroutine test_fold_line(error)
      type(error_type), allocatable, intent(out) :: error
      type(ft_t) :: ft, fs
      integer :: ksub, ni, nj, nx, ny, is, ie, js, je, jf
      real(wp) :: dv_row, du_g, dv_g, d12, dsp, dsm, ss, w12, wsp, dstep
      real(wp), allocatable :: u1(:, :), v1(:, :)

      ni = 16
      nj = 6
      call ft_alloc(ft, ni, nj, nj, 6)
      call sym_geometry(ft, .true.)
      call sym_state(ft, .true.)
      call ft_alloc(fs, ni, nj, nj, 6)
      call sym_geometry(fs, .true.)
      call sym_state(fs, .true.)
      nx = ft%nx; ny = ft%ny; is = ft%is; ie = ft%ie; js = ft%js; je = ft%je
      jf = je + 1
      allocate (u1(nx + 1, ny), v1(nx, ny + 1))
      call ft_map(ft)
      call ft_map(fs)
      ! the reference: one driver step with the fold
      call ft_step(fs, .true.)
      call ft_pull(fs)
      ! the same step by hand, inspecting the new pair after every subcycle
      call ft_prep(ft)
      call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, .true., .false., ft%u, ft%v)
      call ice_evp_c_fold_faces(nx, ny, is, ie, js, je, ft%u, ft%v)
      call ice_evp_c_copy_faces(nx, ny, ft%u, ft%v, ft%w%u2, ft%w%v2)
      w12 = 0.0_wp
      wsp = 0.0_wp
      do ksub = 1, ft%k%ndte
         if (mod(ksub, 2) == 1) then
            call ft_sub(ft, ft%u, ft%v, ft%w%u2, ft%w%v2)
            call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, .true., .false., ft%w%u2, ft%w%v2)
            call ice_evp_c_fold_faces(nx, ny, is, ie, js, je, ft%w%u2, ft%w%v2)
         else
            call ft_sub(ft, ft%w%u2, ft%w%v2, ft%u, ft%v)
            call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, .true., .false., ft%u, ft%v)
            call ice_evp_c_fold_faces(nx, ny, is, ie, js, je, ft%u, ft%v)
         end if
         call ft_pull(ft)
         call ft_pull_pair2(ft)
         if (mod(ksub, 2) == 1) then
            u1 = ft%w%u2
            v1 = ft%w%v2
         else
            u1 = ft%u
            v1 = ft%v
         end if
         call fold_line_report_arrays(ft, u1, v1, dv_row, du_g, dv_g, d12, dsp, dsm, ss)
         w12 = max(w12, d12/max(ss, tiny(1.0_wp)))
         wsp = max(wsp, max(dsp, dsm)/max(ss, tiny(1.0_wp)))
         write (error_unit, '(a,i2,a,3es10.2,a,3es10.2)') "  [evp_c_fold] subcycle ", ksub, ": v row, u ghost, v ghost defects ", &
            dv_row, du_g, dv_g, "   relative corner-slot / ghost-ring stress defects ", d12/max(ss, tiny(1.0_wp)), &
            dsp/max(ss, tiny(1.0_wp)), dsm/max(ss, tiny(1.0_wp))
         call check(error, dv_row == 0.0_wp .and. du_g == 0.0_wp .and. dv_g == 0.0_wp, &
                    "the new face pair is exactly fold-consistent after every subcycle")
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         call check(error, w12 < 1.0e-12_wp .and. wsp < 1.0e-12_wp, "duplicated corner stress slots and ghost-ring stresses agree to round-off")
      end if
      if (.not. allocated(error)) then
         ! the by-hand loop is the driver step
         if (mod(ft%k%ndte, 2) == 1) call ice_evp_c_copy_faces(nx, ny, ft%w%u2, ft%w%v2, ft%u, ft%v)
         call ft_pull(ft)
         dstep = max(maxabs_diff(ft%u, fs%u, 1, nx + 1, 1, ny), maxabs_diff(ft%v, fs%v, 1, nx, 1, ny + 1), &
                     maxabs_diff(ft%w%stress12U, fs%w%stress12U, 1, nx + 1, 1, ny + 1))
         call check(error, dstep == 0.0_wp, "the by-hand subcycle loop reproduces ice_evp_c_step(fold) bitwise")
      end if
      deallocate (u1, v1)
      call ft_destroy(ft)
      call ft_destroy(fs)
   end subroutine test_fold_line

   ! ------------------------------------------------------------------
   ! 6. the sphere
   ! ------------------------------------------------------------------

   subroutine sphere_geometry(ft, dlam, dphi, phi0, consistent_dxu)
      !! Lon-lat metrics on a sphere of radius `R`: `dx = R cos(phi) dlam` at the latitude
      !! of each stagger (T and E at the cell-centre latitude, N and U at the face latitude),
      !! `dy = R dphi`, the exact band areas.  `consistent_dxu`: the corner length is the
      !! earea-weighted mean of the two E lengths it sits between (the identity the solver's
      !! metric term of the corner shear needs for an exactly zero rigid-rotation strain).
      type(ft_t), intent(inout) :: ft
      real(wp), intent(in) :: dlam, dphi, phi0
      logical, intent(in) :: consistent_dxu
      real(wp), parameter :: rr = 6.371e6_wp
      integer :: i, j
      real(wp) :: pc, pn, ps, pnn

      do j = 1, ft%ny
         pc = phi0 + (real(j - ft%js, wp) + 0.5_wp)*dphi
         pn = pc + 0.5_wp*dphi
         ps = pc - 0.5_wp*dphi
         pnn = pc + dphi
         do i = 1, ft%nx
            ft%dxT(i, j) = rr*cos(pc)*dlam
            ft%dxE(i, j) = rr*cos(pc)*dlam
            ft%dxN(i, j) = rr*cos(pn)*dlam
            ft%dxU(i, j) = rr*cos(pn)*dlam
            ft%dyT(i, j) = rr*dphi
            ft%dyE(i, j) = rr*dphi
            ft%dyN(i, j) = rr*dphi
            ft%dyU(i, j) = rr*dphi
            ft%tarea(i, j) = rr*rr*dlam*(sin(pn) - sin(ps))
            ft%earea(i, j) = ft%tarea(i, j)
            ft%narea(i, j) = rr*rr*dlam*(sin(pnn) - sin(pc))
            ft%uarea(i, j) = ft%narea(i, j)
            ft%tarear(i, j) = 1.0_wp/ft%tarea(i, j)
            ft%earear(i, j) = 1.0_wp/ft%earea(i, j)
            ft%narear(i, j) = 1.0_wp/ft%narea(i, j)
         end do
      end do
      if (consistent_dxu) then
         do j = 1, ft%ny - 1
            do i = 1, ft%nx
               ft%dxU(i, j) = (ft%dxE(i, j)*ft%earea(i, j) + ft%dxE(i, j + 1)*ft%earea(i, j + 1))/ &
                              (ft%earea(i, j) + ft%earea(i, j + 1))
            end do
         end do
      end if
      ft%hm = 1.0_wp
      call ft_masks(ft)
      call ft_ratios(ft)
   end subroutine sphere_geometry

   subroutine sphere_rotation_run(nn, consistent, shr_rel, divmax, stressmax, strintmax, p_scale)
      !! Rigid rotation about the pole, `u = Omega R cos(phi)` at E, `v = 0`, on `nn x nn`
      !! cells of a lon-lat sector: the corner shear rate of the stress path
      !! (`shearU / uarea`, 1/s) relative to `Omega`, the T divergence, the stresses after
      !! one evaluation and the stress divergence at E/N relative to `P / dx`.
      integer, intent(in) :: nn
      logical, intent(in) :: consistent
      real(wp), intent(out) :: shr_rel, divmax, stressmax, strintmax, p_scale
      type(ft_t) :: ft
      real(wp), parameter :: rr = 6.371e6_wp, om = OM_EARTH, phi0 = 0.3_wp, span = 0.5_wp
      real(wp) :: dphi, pc, c
      integer :: i, j

      dphi = span/real(nn, wp)
      call ft_alloc(ft, nn, nn, nn, 4)
      call sphere_geometry(ft, dphi, dphi, phi0, consistent)
      ft%aice = 0.9_wp
      ft%aice_init = ft%aice
      ft%vice = 0.9_wp
      ft%aice0 = 0.1_wp
      ft%aicen(:, :, 1) = ft%aice
      ft%vicen(:, :, 1) = ft%vice
      do j = 1, ft%ny
         pc = phi0 + (real(j - ft%js, wp) + 0.5_wp)*dphi
         c = om/dphi
         do i = 1, ft%nx + 1
            ft%u(i, j) = c*ft%dxE(1, j)
         end do
      end do
      ft%v = 0.0_wp
      call ft_map(ft)
      call ft_prep(ft)
      !$acc update self(ft%w%strength)
      ft%w%strength = 1.0e4_wp
      !$acc update device(ft%w%strength)
      call ice_evp_c_average(ft%w, ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%u, ft%v, ft%earea, ft%narea, ft%epm, ft%npm, ft%uvm)
      call ft_k123(ft)
      call ice_evp_c_div_stress(ft%w, ft%k, ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%dxE, ft%dyE, ft%dxN, ft%dyN, &
                                ft%dxU, ft%dyU, ft%dxT, ft%dyT, ft%earear, ft%narear)
      call ice_evp_c_finish(ft%w, ft%k, ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, ft%u, ft%v, ft%dxN, ft%dyE, ft%dxT, ft%dyT, &
                            ft%uarea, ft%tarear, ft%uocnE, ft%vocnE, ft%uocnN, ft%vocnN, ft%ice_err)
      call ft_pull(ft)
      ! corner shear of the stress path, in 1/s: shearU is the rate times the corner area
      shr_rel = 0.0_wp
      do j = ft%js, ft%je
         do i = ft%is, ft%ie
            shr_rel = max(shr_rel, abs(ft%w%shearU(i, j))/ft%uarea(i, j))
         end do
      end do
      shr_rel = shr_rel/om
      divmax = maxabs(ft%w%divu, ft%is, ft%ie, ft%js, ft%je)/om
      stressmax = max(maxabs(ft%w%stresspT, ft%is, ft%ie, ft%js, ft%je), maxabs(ft%w%stressmT, ft%is, ft%ie, ft%js, ft%je), &
                      maxabs(ft%w%stress12U, ft%is, ft%ie, ft%js, ft%je))
      p_scale = 1.0e4_wp
      strintmax = max(maxabs(ft%w%strintxE, ft%is, ft%ie, ft%js, ft%je), maxabs(ft%w%strintyN, ft%is, ft%ie, ft%js, ft%je))
      strintmax = strintmax*ft%dxT(ft%is, ft%js + nn/2)/p_scale
      call ft_destroy(ft)
   end subroutine sphere_rotation_run

   subroutine test_solid_body_sphere(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: sh(3), dv(3), st(3), si(3), pp, s2(3), d2(3), t2(3), i2(3), bound
      integer :: k

      do k = 1, 3
         call sphere_rotation_run(8*2**(k - 1), .true., sh(k), dv(k), st(k), si(k), pp)
         call sphere_rotation_run(8*2**(k - 1), .false., s2(k), d2(k), t2(k), i2(k), pp)
      end do
    write (error_unit, '(a,3es10.2)') "  [evp_c_fold] sphere, metric-consistent corner lengths: corner shear/Omega (N=8,16,32) ", sh
      write (error_unit, '(a,3es10.2)') "      divergence/Omega ", dv
      write (error_unit, '(a,3es10.2)') "      max stress component (P = 1e4 N/m) ", st
      write (error_unit, '(a,3es10.2)') "      stress divergence * dx / P ", si
      write (error_unit, '(a,3es10.2,a,2f7.3)') "  [evp_c_fold] sphere, natural lengths: corner shear/Omega ", s2, " ratios ", s2(1)/s2(2), s2(2)/s2(3)
      write (error_unit, '(a,3es10.2)') "      divergence/Omega ", d2
      write (error_unit, '(a,3es10.2)') "      stress divergence * dx / P (plastic response to the O(h^2) strain) ", i2
      call check(error, maxval(dv) < 1.0e-12_wp .and. maxval(d2) < 1.0e-12_wp, "the rigid rotation has zero divergence")
      if (.not. allocated(error)) call check(error, maxval(sh) < 1.0e-12_wp, &
                                             "metric-consistent corner lengths: zero strain to round-off")
      ! The residual strain of the metric-consistent case is machine round-off (CPU: relative 1e-15 of Omega; the
      ! device contracts multiply-adds so the exact zero of the divergence becomes 1e-15 Omega too).  The rheology
      ! turns a strain `eps` below the floor into the stress `P eps / Delta_min`, so the stress relative to `P` is
      ! bounded by `eps_round / Delta_min` with `eps_round = 50 epsilon Omega`: no tolerance fitted to the output.
      bound = 50.0_wp*epsilon(1.0_wp)*OM_EARTH/DELTAMIN
      if (.not. allocated(error)) call check(error, maxval(si) < bound .and. maxval(st) < bound*1.0e4_wp, &
                                         "metric-consistent: stress and stress divergence at the round-off / viscosity-floor level")
      if (.not. allocated(error)) call check(error, s2(1)/s2(2) > 3.5_wp .and. s2(1)/s2(2) < 4.5_wp .and. &
                                             s2(2)/s2(3) > 3.5_wp .and. s2(2)/s2(3) < 4.5_wp, &
                                             "natural spherical lengths: the corner shear residual is second order")
   end subroutine test_solid_body_sphere

   subroutine drift_solution(f, cw, rhow, a, cosw, sinw, mass, taux, tauy, u, v)
      !! Fixed point of the local free-drift balance (no ocean motion, no tilt, no inertia):
      !! `(vr cosw) U + (f m + vr sinw) k x U = F` rotated, `vr = a rho_w cw |U|`
      !! (the complex form of `test_cice_evp_c`).
      real(wp), intent(in) :: f, cw, rhow, a, cosw, sinw, mass, taux, tauy
      real(wp), intent(out) :: u, v
      real(wp) :: ur, vr, vrel, cr, ci, nr, ni, den
      integer :: m

      ur = 0.0_wp
      vr = 0.0_wp
      do m = 1, 6000
         vrel = a*rhow*cw*sqrt(ur*ur + vr*vr)
         cr = vrel*cosw
         ci = f*mass + vrel*sinw
         nr = taux
         ni = tauy
         den = cr*cr + ci*ci
         ur = 0.5_wp*ur + 0.5_wp*((nr*cr + ni*ci)/den)
         vr = 0.5_wp*vr + 0.5_wp*((ni*cr - nr*ci)/den)
      end do
      u = ur
      v = vr
   end subroutine drift_solution

   subroutine drift_error(nn, err_u, err_v)
      !! Free drift on a lon-lat sector (zero strength, ocean at rest): iterate the
      !! subcycle to its fixed point with Dirichlet ghost rows at the analytic solution
      !! and compare the owned faces to the local analytic balance (the Coriolis
      !! parameter of each face's latitude, the face's own wind).
      integer, intent(in) :: nn
      real(wp), intent(out) :: err_u, err_v
      type(ft_t) :: ft
      real(wp), parameter :: om = OM_EARTH, phi0 = 0.9_wp, span = 0.3_wp
      real(wp), parameter :: taux = 0.12_wp, tauy = -0.05_wp
      real(wp) :: dphi, pc, pn, mass, ua, va, cosw, sinw, a0
      integer :: i, j, it, nit
      real(wp), allocatable :: ue(:), vn(:)

      dphi = span/real(nn, wp)
      call ft_alloc(ft, 8, nn, nn, 2)
      ft%p%dyn%cosw = cos(25.0_wp*PI/180.0_wp)
      ft%p%dyn%sinw = sin(25.0_wp*PI/180.0_wp)
      ft%p%dyn%arlx = 20.0_wp
      ft%p%dyn%brlx = 2.0_wp
      ft%k = ice_evp_c_constants(ft%p, ft%dt)
      cosw = ft%p%dyn%cosw
      sinw = ft%p%dyn%sinw
      call sphere_geometry(ft, dphi, dphi, phi0, .false.)
      a0 = 0.9_wp
      ft%aice = a0
      ft%aice_init = ft%aice
      ft%vice = a0*1.0_wp
      ft%aice0 = 1.0_wp - a0
      ft%aicen(:, :, 1) = ft%aice
      ft%vicen(:, :, 1) = ft%vice
      ft%strairxE = taux
      ft%strairyN = tauy
      mass = ft%p%c%rhoi*a0*1.0_wp
      ! Coriolis parameter at the faces' own latitudes, analytic velocity everywhere (ghost rows too)
      allocate (ue(ft%ny), vn(ft%ny + 1))
      do j = 1, ft%ny
         pc = phi0 + (real(j - ft%js, wp) + 0.5_wp)*dphi
         pn = pc + 0.5_wp*dphi
         ft%fcorE(:, j) = 2.0_wp*om*sin(pc)
         ft%fcorN(:, j) = 2.0_wp*om*sin(pn)
         call drift_solution(2.0_wp*om*sin(pc), ft%p%dyn%cdn_ocn, ft%p%c%rhow, a0, cosw, sinw, mass, taux, tauy, ua, va)
         ue(j) = ua
         call drift_solution(2.0_wp*om*sin(pn), ft%p%dyn%cdn_ocn, ft%p%c%rhow, a0, cosw, sinw, mass, taux, tauy, ua, va)
         ft%u(:, j) = 0.0_wp
         ft%v(:, j + 1) = va
         vn(j + 1) = va
      end do
      do j = 1, ft%ny
         ft%u(:, j) = ue(j)
      end do
      ft%v(:, 1) = ft%v(:, 2)
      call ft_map(ft)
      call ft_prep(ft)
      !$acc update self(ft%w%strength)
      ft%w%strength = 0.0_wp
      !$acc update device(ft%w%strength)
      nit = 4000
      call ice_evp_c_wrap_faces(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, .true., .false., ft%u, ft%v)
      call ice_evp_c_copy_faces(ft%nx, ft%ny, ft%u, ft%v, ft%w%u2, ft%w%v2)
      do it = 1, nit
         if (mod(it, 2) == 1) then
            call ft_sub(ft, ft%u, ft%v, ft%w%u2, ft%w%v2)
            call ice_evp_c_wrap_faces(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, .true., .false., ft%w%u2, ft%w%v2)
         else
            call ft_sub(ft, ft%w%u2, ft%w%v2, ft%u, ft%v)
            call ice_evp_c_wrap_faces(ft%nx, ft%ny, ft%is, ft%ie, ft%js, ft%je, .true., .false., ft%u, ft%v)
         end if
      end do
      if (mod(nit, 2) == 1) call ice_evp_c_copy_faces(ft%nx, ft%ny, ft%w%u2, ft%w%v2, ft%u, ft%v)
      call ft_pull(ft)
      err_u = 0.0_wp
      err_v = 0.0_wp
      do j = ft%js + 1, ft%je - 1
         pc = phi0 + (real(j - ft%js, wp) + 0.5_wp)*dphi
         pn = pc + 0.5_wp*dphi
         call drift_solution(2.0_wp*om*sin(pc), ft%p%dyn%cdn_ocn, ft%p%c%rhow, a0, cosw, sinw, mass, taux, tauy, ua, va)
         do i = ft%is, ft%ie
            err_u = max(err_u, abs(ft%u(i + 1, j) - ua))
         end do
         call drift_solution(2.0_wp*om*sin(pn), ft%p%dyn%cdn_ocn, ft%p%c%rhow, a0, cosw, sinw, mass, taux, tauy, ua, va)
         do i = ft%is, ft%ie
            err_v = max(err_v, abs(ft%v(i, j + 1) - va))
         end do
      end do
      deallocate (ue, vn)
      call ft_destroy(ft)
   end subroutine drift_error

   subroutine test_free_drift_sphere(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: eu(3), ev(3)
      integer :: k

      do k = 1, 3
         call drift_error(6*2**(k - 1), eu(k), ev(k))
      end do
      write (error_unit, '(a,3es10.2,a,2f7.3)') "  [evp_c_fold] sphere free drift, max |u - u_analytic| (N=6,12,24) ", eu, &
         " ratios ", eu(1)/eu(2), eu(2)/eu(3)
      write (error_unit, '(a,3es10.2,a,2f7.3)') "      max |v - v_analytic| ", ev, " ratios ", ev(1)/ev(2), ev(2)/ev(3)
      call check(error, all(eu < 5.0e-3_wp) .and. all(ev < 5.0e-3_wp), "free drift on the sphere is near the local balance")
      if (.not. allocated(error)) call check(error, eu(1)/eu(2) > 3.0_wp .and. eu(2)/eu(3) > 3.0_wp .and. &
                                             ev(1)/ev(2) > 3.0_wp .and. ev(2)/ev(3) > 3.0_wp, &
                                             "the departure from the local balance is second order in the cell size")
   end subroutine test_free_drift_sphere

end module test_cice_evp_c_fold
