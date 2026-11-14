!! C-grid EVP dynamics (CICE6 port D9.1): analytic, symmetry, mask, convergence
!! and device-parity tests of `rdb_ice_evp_c` / `rdb_ice_evp_c_prep`.
!!
!! Links only the ice core library.  Every array a kernel touches is mapped
!! (`mem:separate`), results are pulled back with `update self`.  Expected
!! numbers are closed forms written out here, a scalar re-transcription of the
!! per-point formulas on ice-name arrays (`ref_run`, host loops), and the
!! numbers of the stdlib-only Python transcript
!! (`python_prototypes/cice_evp_c/evp_c_transcript.py`).
!!
!! Cases:
!!   * `dyn_params`        defaults, every refusal code, derived constants.
!!   * `uniform_velocity`  uniform u,v: zero strain, zero stress, zero divergence.
!!   * `solid_body_rotation` zero divergence/shear, vorticity 2*Omega.
!!   * `linear_field_*`    exact divu/shear/vort/rdg_*, viscosities, replacement
!!                         pressure and the first stress update in closed form
!!                         (classic + revised EVP, capping max/sum, Ktens,
!!                         avg_zeta/avg_strength).
!!   * `div_stress_linear` stress divergence of prescribed linear stresses;
!!                         `free_drift_approx` gating.
!!   * `free_drift_*`      momentum balance of wind, tilt, ocean drag (with a
!!                         turning angle) and Coriolis, closed form.
!!   * `mirror_x/y`        mirrored inputs give the mirrored solution.
!!   * `periodic_shift`    translation invariance of the periodic channel.
!!   * `closed_basin_*`    onshore wind: no flow through walls, bounded stress.
!!   * `ndte_convergence`  residual against a long-subcycle reference decreases.
!!   * `zero_ice_cells`    ice-free cells untouched.
!!   * `nan_*`             non-finite input is flagged, never laundered.
!!   * `python_transcript` T stress and one momentum update vs the Python numbers.
!!   * `reference_*`       every output on every cell vs the scalar re-transcription
!!                         on a random masked state (classic, revised, avg_strength).
!!   * `throughput`        opt-in (`RDB_CICE_EVP_BENCH=1`).
module test_cice_evp_c
   use, intrinsic :: iso_fortran_env, only: int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_errors, only: ICE_ERR_DYN_VEL
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_dyn_params_t, ice_dyn_params_default, &
                             ice_dyn_params_check, ICE_CONSTANTS_ICEPACK, ICE_VISC_AVG_ZETA, &
                             ICE_VISC_AVG_STRENGTH, ICE_DYN_NMSG
   use rdb_ice_evp_c_prep, only: ice_evp_c_const_t, ice_evp_c_work_t, ice_evp_c_constants, ice_evp_c_prep, &
                                 ice_evp_c_ghost_ok
   use rdb_ice_evp_c, only: ice_evp_c_step, ice_evp_c_substep, ice_evp_c_average, ice_evp_c_finish, &
                            ice_evp_c_wrap_faces, ice_evp_c_div_stress, evp_strain_u_impl, evp_stress_t_impl, &
                            evp_stress_u_impl
   implicit none
   private

   public :: collect_cice_evp_c_tests

   integer, parameter :: NG = 3
   real(wp), parameter :: PI = 3.14159265358979323846_wp
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: TOL_REF = 1.0e-12_wp   ! FMA on host code too; measured worst 1e-15 relative (V100)
   real(wp), parameter :: TOL_CF = 1.0e-11_wp
#else
   real(wp), parameter :: TOL_REF = 1.0e-13_wp   ! measured worst 2e-16 relative
   real(wp), parameter :: TOL_CF = 1.0e-12_wp
#endif

   type :: fx_t
      !! Cartesian fixture: ice-name arrays, an `ice_evp_c_work_t` and the constants.
      integer :: nxp = 0, nyp = 0, nx = 0, ny = 0, is = 0, ie = 0, js = 0, je = 0
      real(wp) :: dx = 0.0_wp, dy = 0.0_wp, dt = 900.0_wp
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
   end type fx_t

contains

   subroutine collect_cice_evp_c_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("dyn_params", test_dyn_params), &
                  new_unittest("ghost_width_check", test_ghost_ok), &
                  new_unittest("uniform_velocity", test_uniform_velocity), &
                  new_unittest("solid_body_rotation", test_solid_body_rotation), &
                  new_unittest("linear_field_classic", test_linear_classic), &
                  new_unittest("linear_field_sum_ktens", test_linear_sum_ktens), &
                  new_unittest("linear_field_revised", test_linear_revised), &
                  new_unittest("linear_field_avg_strength", test_linear_avg_strength), &
                  new_unittest("div_stress_linear", test_div_stress_linear), &
                  new_unittest("free_drift_balance", test_free_drift_balance), &
                  new_unittest("free_drift_tilt_sign", test_tilt_sign), &
                  new_unittest("mirror_x", test_mirror_x), &
                  new_unittest("mirror_y", test_mirror_y), &
                  new_unittest("periodic_shift", test_periodic_shift), &
                  new_unittest("closed_basin_wind", test_closed_basin), &
                  new_unittest("ndte_convergence", test_ndte_convergence), &
                  new_unittest("zero_ice_cells", test_zero_ice), &
                  new_unittest("nan_strength", test_nan_strength), &
                  new_unittest("nan_velocity", test_nan_velocity), &
                  new_unittest("python_transcript", test_python_transcript), &
                  new_unittest("coastal_corner_no_slip", test_coastal_corner), &
                  new_unittest("spherical_rigid_rotation", test_spherical_rotation), &
                  new_unittest("reference_classic", test_reference_classic), &
                  new_unittest("reference_revised", test_reference_revised), &
                  new_unittest("reference_avg_strength_sum", test_reference_avg_strength), &
                  new_unittest("throughput", test_throughput) &
                  ]
   end subroutine collect_cice_evp_c_tests

   ! ------------------------------------------------------------------
   ! fixture
   ! ------------------------------------------------------------------

   subroutine fx_build(fx, nxp, nyp, wall, dx, dy, ndte)
      !! Allocate, zero and fill the metrics/masks of an `nxp x nyp` Cartesian
      !! domain with `NG` ghost rings.  `wall`: the ghost band is land (a closed
      !! box); otherwise every cell and face is open.  Ice is NOT set (call
      !! `fx_ice`); the work arrays are mapped by `fx_map`.
      type(fx_t), intent(inout) :: fx
      integer, intent(in) :: nxp, nyp, ndte
      logical, intent(in) :: wall
      real(wp), intent(in) :: dx, dy
      integer :: nx, ny, i, j
      real(wp), allocatable :: tm(:, :)

      fx%nxp = nxp
      fx%nyp = nyp
      nx = nxp + 2*NG
      ny = nyp + 2*NG
      fx%nx = nx
      fx%ny = ny
      fx%is = NG + 1
      fx%ie = NG + nxp
      fx%js = NG + 1
      fx%je = NG + nyp
      fx%dx = dx
      fx%dy = dy
      allocate (fx%dxE(nx, ny), fx%dyE(nx, ny), fx%dxN(nx, ny), fx%dyN(nx, ny), fx%dxU(nx, ny), fx%dyU(nx, ny))
      allocate (fx%dxT(nx, ny), fx%dyT(nx, ny), fx%uarea(nx, ny), fx%tarea(nx, ny), fx%tarear(nx, ny))
      allocate (fx%earea(nx, ny), fx%narea(nx, ny), fx%earear(nx, ny), fx%narear(nx, ny))
      allocate (fx%rdxN(nx, ny), fx%rdxNr(nx, ny), fx%rdyE(nx, ny), fx%rdyEr(nx, ny))
      allocate (fx%epm(nx, ny), fx%npm(nx, ny), fx%uvm(nx, ny), fx%hm(nx, ny))
      allocate (fx%aice(nx, ny), fx%vice(nx, ny), fx%vsno(nx, ny), fx%aice_init(nx, ny), fx%aice0(nx, ny))
      allocate (fx%aicen(nx, ny, 1), fx%vicen(nx, ny, 1))
      allocate (fx%fcorE(nx, ny), fx%fcorN(nx, ny), fx%strairxE(nx, ny), fx%strairyN(nx, ny))
      allocate (fx%ss_tltxE(nx, ny), fx%ss_tltyN(nx, ny))
      allocate (fx%uocnE(nx, ny), fx%vocnE(nx, ny), fx%uocnN(nx, ny), fx%vocnN(nx, ny))
      allocate (fx%u(nx + 1, ny), fx%v(nx, ny + 1), fx%ice_err(nx, ny))
      fx%dxE = dx; fx%dxN = dx; fx%dxU = dx; fx%dxT = dx
      fx%dyE = dy; fx%dyN = dy; fx%dyU = dy; fx%dyT = dy
      fx%uarea = dx*dy; fx%tarea = dx*dy; fx%earea = dx*dy; fx%narea = dx*dy
      fx%tarear = 1.0_wp/(dx*dy); fx%earear = fx%tarear; fx%narear = fx%tarear
      fx%rdxN = -1.0_wp; fx%rdxNr = -1.0_wp; fx%rdyE = -1.0_wp; fx%rdyEr = -1.0_wp
      fx%rdxN(nx, :) = 0.0_wp; fx%rdxNr(nx, :) = 0.0_wp
      fx%rdyE(:, ny) = 0.0_wp; fx%rdyEr(:, ny) = 0.0_wp
      allocate (tm(nx, ny))
      tm = 1.0_wp
      if (wall) then
         do j = 1, ny
            do i = 1, nx
               if (i < fx%is .or. i > fx%ie .or. j < fx%js .or. j > fx%je) tm(i, j) = 0.0_wp
            end do
         end do
      end if
      fx%hm = tm
      fx%epm = 0.0_wp; fx%npm = 0.0_wp; fx%uvm = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            fx%epm(i, j) = tm(i, j)*tm(i + 1, j)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            fx%npm(i, j) = tm(i, j)*tm(i, j + 1)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx - 1
            fx%uvm(i, j) = tm(i, j)*tm(i + 1, j)*tm(i, j + 1)*tm(i + 1, j + 1)
         end do
      end do
      deallocate (tm)
      fx%aice = 0.0_wp; fx%vice = 0.0_wp; fx%vsno = 0.0_wp; fx%aice_init = 0.0_wp; fx%aice0 = 1.0_wp
      fx%aicen = 0.0_wp; fx%vicen = 0.0_wp
      fx%fcorE = 0.0_wp; fx%fcorN = 0.0_wp; fx%strairxE = 0.0_wp; fx%strairyN = 0.0_wp
      fx%ss_tltxE = 0.0_wp; fx%ss_tltyN = 0.0_wp
      fx%uocnE = 0.0_wp; fx%vocnE = 0.0_wp; fx%uocnN = 0.0_wp; fx%vocnN = 0.0_wp
      fx%u = 0.0_wp; fx%v = 0.0_wp; fx%ice_err = 0_int32
      fx%p = ice_params_default(ICE_CONSTANTS_ICEPACK, 1, 4, 1)
      fx%p%ridge%kstrength = 0
      fx%p%dyn%ndte = ndte
      fx%dt = 900.0_wp
      fx%k = ice_evp_c_constants(fx%p, fx%dt)
      call fx%w%init(nx, ny)
   end subroutine fx_build

   subroutine fx_consts(fx)
      !! Re-derive the constants after editing `fx%p%dyn`.
      type(fx_t), intent(inout) :: fx
      fx%k = ice_evp_c_constants(fx%p, fx%dt)
   end subroutine fx_consts

   subroutine fx_ice(fx, a, h)
      !! Uniform ice of concentration `a` and thickness `h` on every wet cell.
      type(fx_t), intent(inout) :: fx
      real(wp), intent(in) :: a, h
      fx%aice = a*fx%hm
      fx%aice_init = fx%aice
      fx%vice = a*h*fx%hm
      fx%vsno = 0.0_wp
      fx%aicen(:, :, 1) = fx%aice
      fx%vicen(:, :, 1) = fx%vice
      fx%aice0 = 1.0_wp - fx%aice
   end subroutine fx_ice

   subroutine fx_map(fx)
      !! Map every array to the device (current host values).
      type(fx_t), intent(inout) :: fx
      call fx%w%enter_data()
      !$acc enter data copyin(fx%dxE, fx%dyE, fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT)
      !$acc enter data copyin(fx%uarea, fx%tarea, fx%tarear, fx%earea, fx%narea, fx%earear, fx%narear)
      !$acc enter data copyin(fx%rdxN, fx%rdxNr, fx%rdyE, fx%rdyEr, fx%epm, fx%npm, fx%uvm, fx%hm)
      !$acc enter data copyin(fx%aice, fx%vice, fx%vsno, fx%aice_init, fx%aice0, fx%aicen, fx%vicen)
      !$acc enter data copyin(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc enter data copyin(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%u, fx%v, fx%ice_err)
      !$acc enter data copyin(fx%p)
      fx%mapped = .true.
   end subroutine fx_map

   subroutine fx_unmap(fx)
      type(fx_t), intent(inout) :: fx
      !$acc exit data delete(fx%dxE, fx%dyE, fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT)
      !$acc exit data delete(fx%uarea, fx%tarea, fx%tarear, fx%earea, fx%narea, fx%earear, fx%narear)
      !$acc exit data delete(fx%rdxN, fx%rdxNr, fx%rdyE, fx%rdyEr, fx%epm, fx%npm, fx%uvm, fx%hm)
      !$acc exit data delete(fx%aice, fx%vice, fx%vsno, fx%aice_init, fx%aice0, fx%aicen, fx%vicen)
      !$acc exit data delete(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc exit data delete(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%u, fx%v, fx%ice_err)
      !$acc exit data delete(fx%p)
      call fx%w%exit_data()
      fx%mapped = .false.
   end subroutine fx_unmap

   subroutine fx_push_in(fx)
      !! Host -> device: the inputs and the velocity.
      type(fx_t), intent(inout) :: fx
      !$acc update device(fx%aice, fx%vice, fx%vsno, fx%aice_init, fx%aice0, fx%aicen, fx%vicen)
      !$acc update device(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc update device(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%u, fx%v, fx%ice_err)
   end subroutine fx_push_in

   subroutine fx_push_w(fx)
      !! Host -> device: the persistent state and masks of the work type.
      type(fx_t), intent(inout) :: fx
      !$acc update device(fx%w%stresspT, fx%w%stressmT, fx%w%zetax2T, fx%w%etax2T, fx%w%stress12U)
      !$acc update device(fx%w%iceEmask, fx%w%iceNmask, fx%w%strength, fx%w%strengthU)
   end subroutine fx_push_w

   subroutine fx_pull(fx)
      !! Device -> host: everything a test inspects.
      type(fx_t), intent(inout) :: fx
      !$acc update self(fx%u, fx%v, fx%ice_err)
      !$acc update self(fx%w%stresspT, fx%w%stressmT, fx%w%zetax2T, fx%w%etax2T, fx%w%stress12U)
      !$acc update self(fx%w%iceEmask, fx%w%iceNmask, fx%w%tmass, fx%w%tmphm, fx%w%iceTmask, fx%w%iceUmask)
      !$acc update self(fx%w%strength, fx%w%strengthU, fx%w%aiE, fx%w%aiN, fx%w%emassdti, fx%w%nmassdti)
      !$acc update self(fx%w%fmE, fx%w%fmN, fx%w%forcexE, fx%w%forceyN, fx%w%uinit, fx%w%vinit)
      !$acc update self(fx%w%uvelN, fx%w%vvelE, fx%w%uvelU, fx%w%vvelU, fx%w%shearU, fx%w%deltaU)
      !$acc update self(fx%w%strintxE, fx%w%strintyN, fx%w%divu, fx%w%shear, fx%w%vort)
      !$acc update self(fx%w%rdg_conv, fx%w%rdg_shear, fx%w%strocnxE, fx%w%strocnyN)
   end subroutine fx_pull

   subroutine fx_destroy(fx)
      type(fx_t), intent(inout) :: fx
      if (fx%mapped) call fx_unmap(fx)
      call fx%w%destroy()
      if (allocated(fx%dxE)) deallocate (fx%dxE, fx%dyE, fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT)
      if (allocated(fx%uarea)) deallocate (fx%uarea, fx%tarea, fx%tarear, fx%earea, fx%narea, fx%earear, fx%narear)
      if (allocated(fx%rdxN)) deallocate (fx%rdxN, fx%rdxNr, fx%rdyE, fx%rdyEr, fx%epm, fx%npm, fx%uvm, fx%hm)
      if (allocated(fx%aice)) deallocate (fx%aice, fx%vice, fx%vsno, fx%aice_init, fx%aice0, fx%aicen, fx%vicen)
      if (allocated(fx%fcorE)) deallocate (fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      if (allocated(fx%uocnE)) deallocate (fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%u, fx%v, fx%ice_err)
   end subroutine fx_destroy

   subroutine fx_prep(fx)
      type(fx_t), intent(inout) :: fx
      call ice_evp_c_prep(fx%w, fx%k, fx%nx, fx%ny, 1, fx%is, fx%ie, fx%js, fx%je, fx%p%ridge, fx%p%c, fx%hm, &
                          fx%aice, fx%vice, fx%vsno, fx%aice_init, fx%aice0, fx%aicen, fx%vicen, fx%tarea, &
                          fx%epm, fx%npm, fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, &
                          fx%ss_tltyN, fx%uocnE, fx%vocnN, fx%u, fx%v, fx%ice_err)
   end subroutine fx_prep

   subroutine fx_avg(fx)
      type(fx_t), intent(inout) :: fx
      call ice_evp_c_average(fx%w, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%u, fx%v, fx%earea, fx%narea, &
                             fx%epm, fx%npm, fx%uvm)
   end subroutine fx_avg

   subroutine fx_sub(fx)
      type(fx_t), intent(inout) :: fx
      call ice_evp_c_substep(fx%w, fx%k, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%u, fx%v, fx%dxE, fx%dyE, &
                             fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT, fx%uarea, fx%tarea, fx%earear, &
                             fx%narear, fx%rdxN, fx%rdxNr, fx%rdyE, fx%rdyEr, fx%epm, fx%npm, fx%hm, &
                             fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN)
   end subroutine fx_sub

   subroutine fx_k123(fx)
      !! The three stress kernels of one subcycle without the momentum update
      !! (the velocity stays the prescribed field).
      type(fx_t), intent(inout) :: fx
      call evp_strain_u_impl(fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%k%e_factor, fx%k%avg_strength, fx%u, &
                             fx%v, fx%w%uvelN, fx%w%vvelE, fx%w%uvelU, fx%w%vvelU, fx%dxE, fx%dyN, fx%dxU, fx%dyU, &
                             fx%rdxN, fx%rdxNr, fx%rdyE, fx%rdyEr, fx%epm, fx%npm, fx%w%iceUmask, fx%w%shearU, &
                             fx%w%deltaU)
      call evp_stress_t_impl(fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%k%arlx1i, fx%k%revp, fx%k%denom1, &
                             fx%k%epp2i, fx%k%e_factor, fx%k%Ktens, fx%k%deltamin, fx%k%cap_max, fx%u, fx%v, fx%dxN, &
                             fx%dyE, fx%dxT, fx%dyT, fx%uarea, fx%tarea, fx%w%strength, fx%w%shearU, fx%w%iceTmask, &
                             fx%w%zetax2T, fx%w%etax2T, fx%w%stresspT, fx%w%stressmT)
      call evp_stress_u_impl(fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%k%arlx1i, fx%k%revp, fx%k%denom1, &
                             fx%k%epp2i, fx%k%Ktens, fx%k%deltamin, fx%k%cap_max, fx%k%avg_strength, fx%hm, &
                             fx%tarea, fx%uarea, fx%w%etax2T, fx%w%strengthU, fx%w%deltaU, fx%w%shearU, &
                             fx%w%iceUmask, fx%w%stress12U)
   end subroutine fx_k123

   subroutine fx_fin(fx)
      type(fx_t), intent(inout) :: fx
      call ice_evp_c_finish(fx%w, fx%k, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%u, fx%v, fx%dxN, fx%dyE, &
                            fx%dxT, fx%dyT, fx%uarea, fx%tarear, fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%ice_err)
   end subroutine fx_fin

   subroutine fx_wrap(fx, wx, wy)
      type(fx_t), intent(inout) :: fx
      logical, intent(in) :: wx, wy
      call ice_evp_c_wrap_faces(fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, wx, wy, fx%u, fx%v)
   end subroutine fx_wrap

   subroutine fx_step(fx, wx, wy)
      !! One full step through the public driver.
      type(fx_t), intent(inout) :: fx
      logical, intent(in) :: wx, wy
      call ice_evp_c_step(fx%w, fx%k, fx%nx, fx%ny, 1, fx%is, fx%ie, fx%js, fx%je, fx%p%ridge, fx%p%c, wx, wy, &
                          fx%dxE, fx%dyE, fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT, fx%uarea, fx%tarea, &
                          fx%tarear, fx%earea, fx%narea, fx%earear, fx%narear, fx%rdxN, fx%rdxNr, fx%rdyE, &
                          fx%rdyEr, fx%epm, fx%npm, fx%uvm, fx%hm, fx%aice, fx%vice, fx%vsno, fx%aice_init, &
                          fx%aice0, fx%aicen, fx%vicen, fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, &
                          fx%ss_tltxE, fx%ss_tltyN, fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN, fx%u, fx%v, fx%ice_err)
   end subroutine fx_step

   subroutine fx_set_strength(fx, pval)
      !! Overwrite the strength after `fx_prep` (uniform `pval` on the ice mask).
      type(fx_t), intent(inout) :: fx
      real(wp), intent(in) :: pval
      !$acc update self(fx%w%iceTmask)
      fx%w%strength = 0.0_wp
      where (fx%w%iceTmask) fx%w%strength = pval
      !$acc update device(fx%w%strength)
   end subroutine fx_set_strength

   subroutine fx_strength_u(fx)
      !! Host-side T2US of the strength (when a test replaces the strength
      !! after prep and runs `avg_strength`).
      type(fx_t), intent(inout) :: fx
      integer :: i, j
      real(wp) :: w1, w2, w3, w4, wt
      fx%w%strengthU = 0.0_wp
      do j = 1, fx%ny - 1
         do i = 1, fx%nx - 1
            w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i + 1, j)*fx%tarea(i + 1, j)
            w3 = fx%hm(i, j + 1)*fx%tarea(i, j + 1); w4 = fx%hm(i + 1, j + 1)*fx%tarea(i + 1, j + 1)
            wt = w1 + w2 + w3 + w4
            if (wt /= 0.0_wp) fx%w%strengthU(i, j) = (w1*fx%w%strength(i, j) + w2*fx%w%strength(i + 1, j) &
                                                      + w3*fx%w%strength(i, j + 1) + w4*fx%w%strength(i + 1, j + 1))/wt
         end do
      end do
      !$acc update device(fx%w%strengthU)
   end subroutine fx_strength_u

   subroutine near(error, got, want, rtol, atol, msg)
      !! |got - want| <= atol + rtol |want|, NaN-proof.
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: got, want, rtol, atol
      character(len=*), intent(in) :: msg
      character(len=160) :: buf
      if (allocated(error)) return
      if (.not. (ieee_is_finite(got) .and. ieee_is_finite(want))) then
         call test_failed(error, msg//": non-finite value")
         return
      end if
      if (abs(got - want) > atol + rtol*abs(want)) then
         write (buf, '(a,a,es23.15,a,es23.15)') msg, ": got ", got, " want ", want
         call test_failed(error, trim(buf))
      end if
   end subroutine near

   subroutine field_diff(a, b, i0, i1, j0, j1, worst, scale, nbad)
      !! Largest |a-b| over the window, the largest |b|, and the count of
      !! non-finite entries of either (NaN-proof).
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      real(wp), intent(out) :: worst, scale
      integer, intent(out) :: nbad
      integer :: i, j
      worst = 0.0_wp
      scale = 0.0_wp
      nbad = 0
      do j = j0, j1
         do i = i0, i1
            if (.not. (ieee_is_finite(a(i, j)) .and. ieee_is_finite(b(i, j)))) then
               nbad = nbad + 1
            else
               worst = max(worst, abs(a(i, j) - b(i, j)))
               scale = max(scale, abs(b(i, j)))
            end if
         end do
      end do
   end subroutine field_diff

   subroutine expect_field(error, name, a, b, i0, i1, j0, j1, rtol)
      !! `a == b` to `rtol * max|b|` (+ tiny) over a window; prints the measured worst.
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      real(wp), intent(in) :: rtol
      real(wp) :: worst, scale
      integer :: nbad
      character(len=160) :: buf
      if (allocated(error)) return
      call field_diff(a, b, i0, i1, j0, j1, worst, scale, nbad)
      write (error_unit, '(a,a,a,es9.2,a,es9.2)') "  [evp_c] ", name, " worst ", worst, " scale ", scale
      if (nbad > 0) then
         write (buf, '(a,a,i0,a)') name, ": ", nbad, " non-finite entries"
         call test_failed(error, trim(buf))
      else if (worst > rtol*scale + 1.0e-300_wp) then
         write (buf, '(a,a,es10.3,a,es10.3)') name, ": worst ", worst, " exceeds tol*scale ", rtol*scale
         call test_failed(error, trim(buf))
      end if
   end subroutine expect_field

   subroutine fill_linear(fx, a, b, c, d, u0, v0)
      !! `u = a x + b y + u0` at the x faces, `v = c x + d y + v0` at the y faces
      !! (x faces at `x = (i - is) dx`, cell centres at `(j - js + 1/2) dy`).
      type(fx_t), intent(inout) :: fx
      real(wp), intent(in) :: a, b, c, d, u0, v0
      integer :: i, j
      do j = 1, fx%ny
         do i = 1, fx%nx + 1
            fx%u(i, j) = a*real(i - fx%is, wp)*fx%dx + b*(real(j - fx%js, wp) + 0.5_wp)*fx%dy + u0
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx
            fx%v(i, j) = c*(real(i - fx%is, wp) + 0.5_wp)*fx%dx + d*real(j - fx%js, wp)*fx%dy + v0
         end do
      end do
   end subroutine fill_linear

   ! ------------------------------------------------------------------
   ! tests
   ! ------------------------------------------------------------------

   subroutine test_dyn_params(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_dyn_params_t) :: dp
      type(ice_params_t) :: p
      type(ice_evp_c_const_t) :: k
      integer :: codes(ICE_DYN_NMSG), n, c

      dp = ice_dyn_params_default()
      call check(error, dp%ndte, 120)
      if (allocated(error)) return
      call near(error, dp%elasticDamp, 0.36_wp, 0.0_wp, 0.0_wp, "elasticDamp")
      call near(error, dp%e_yieldcurve, 2.0_wp, 0.0_wp, 0.0_wp, "e_yieldcurve")
      call near(error, dp%e_plasticpot, 2.0_wp, 0.0_wp, 0.0_wp, "e_plasticpot")
      call near(error, dp%deltaminEVP, 1.0e-11_wp, 0.0_wp, 0.0_wp, "deltaminEVP")
      call near(error, dp%dyn_area_min, 1.0e-11_wp, 0.0_wp, 0.0_wp, "dyn_area_min")
      call near(error, dp%dyn_mass_min, 1.0e-10_wp, 0.0_wp, 0.0_wp, "dyn_mass_min")
      call near(error, dp%rheo_area_min, 1.0e-3_wp, 0.0_wp, 0.0_wp, "rheo_area_min")
      call near(error, dp%cdn_ocn, 0.00536_wp, 0.0_wp, 0.0_wp, "cdn_ocn")
      call near(error, dp%arlx, 300.0_wp, 0.0_wp, 0.0_wp, "arlx")
      call near(error, dp%brlx, 300.0_wp, 0.0_wp, 0.0_wp, "brlx")
      if (allocated(error)) return
      call check(error, dp%revised_evp .or. .not. dp%capping_max .or. .not. dp%free_drift_approx .or. &
                 dp%visc_method /= ICE_VISC_AVG_ZETA, .false., "default flags")
      if (allocated(error)) return
      call ice_dyn_params_check(dp, codes, n)
      call check(error, n, 0, "defaults are accepted")
      if (allocated(error)) return
      ! every refusal code, one at a time (+ NaN)
      do c = 1, 14
         dp = ice_dyn_params_default()
         select case (c)
         case (1); dp%ndte = 0
         case (2); dp%elasticDamp = 0.0_wp
         case (3); dp%revised_evp = .true.; dp%arlx = -1.0_wp
         case (4); dp%revised_evp = .true.; dp%brlx = 0.0_wp
         case (5); dp%e_yieldcurve = -2.0_wp
         case (6); dp%e_plasticpot = ieee_value(1.0_wp, ieee_quiet_nan)
         case (7); dp%Ktens = 1.0_wp
         case (8); dp%deltaminEVP = 0.0_wp
         case (9); dp%dyn_area_min = -1.0e-3_wp
         case (10); dp%dyn_mass_min = -1.0_wp
         case (11); dp%rheo_area_min = -1.0_wp
         case (12); dp%visc_method = 7
         case (13); dp%cdn_ocn = -0.1_wp
         case (14); dp%vel_cap = 0.0_wp
         end select
         call ice_dyn_params_check(dp, codes, n)
         call check(error, n == 1 .and. codes(1) == c, .true., "refusal code isolated")
         if (allocated(error)) return
      end do
      ! two at once are both reported
      dp = ice_dyn_params_default()
      dp%ndte = -3
      dp%cdn_ocn = -1.0_wp
      call ice_dyn_params_check(dp, codes, n)
      call check(error, n == 2 .and. codes(1) == 1 .and. codes(2) == 13, .true., "all refusals reported")
      if (allocated(error)) return
      ! derived constants: classic (Hunke 2001) and revised (Bouillon et al. 2013)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 7, 1)
      k = ice_evp_c_constants(p, 900.0_wp)
      call near(error, k%arlx1i, 1.0_wp/(2.0_wp*0.36_wp*120.0_wp), 1.0e-15_wp, 0.0_wp, "classic 1/arlx")
      call near(error, k%brlx, 120.0_wp, 0.0_wp, 0.0_wp, "classic brlx = ndte")
      call near(error, k%denom1, 1.0_wp/(1.0_wp + k%arlx1i), 1.0e-15_wp, 0.0_wp, "classic denom1")
      call near(error, k%revp, 0.0_wp, 0.0_wp, 0.0_wp, "revp 0")
      call near(error, k%epp2i, 0.25_wp, 0.0_wp, 0.0_wp, "epp2i")
      call near(error, k%e_factor, 0.25_wp, 0.0_wp, 0.0_wp, "e_factor")
      p%dyn%revised_evp = .true.
      k = ice_evp_c_constants(p, 900.0_wp)
      call near(error, k%arlx1i, 1.0_wp/300.0_wp, 1.0e-15_wp, 0.0_wp, "revised 1/arlx")
      call near(error, k%brlx, 300.0_wp, 0.0_wp, 0.0_wp, "revised brlx")
      call near(error, k%denom1, 1.0_wp, 0.0_wp, 0.0_wp, "revised denom1")
      call near(error, k%revp, 1.0_wp, 0.0_wp, 0.0_wp, "revp 1")
   end subroutine test_dyn_params

   subroutine test_ghost_ok(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, ice_evp_c_ghost_ok(14, 12, 4, 11, 4, 9), .true., "3 ghost cells accepted")
      if (allocated(error)) return
      call check(error,.not. ice_evp_c_ghost_ok(12, 12, 3, 10, 4, 9), .true., "2 ghost cells refused")
   end subroutine test_ghost_ok

   subroutine test_uniform_velocity(error)
      !! Uniform velocity in an open domain: no strain, no stress, no divergence.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j
      real(wp) :: mx

      call fx_build(fx, 8, 7, .false., 2.0e4_wp, 1.5e4_wp, 4)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      fx%u = 0.13_wp
      fx%v = -0.07_wp
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_prep(fx)
      call fx_avg(fx)
      call fx_k123(fx)
      call fx_fin(fx)
      call ice_evp_c_div_stress(fx%w, fx%k, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%dxE, fx%dyE, fx%dxN, &
                                fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT, fx%earear, fx%narear)
      call fx_pull(fx)
      mx = 0.0_wp
      do j = fx%js - 2, fx%je + 1
         do i = fx%is - 2, fx%ie + 1
            mx = max(mx, abs(fx%w%shearU(i, j)))
         end do
      end do
      call near(error, mx, 0.0_wp, 0.0_wp, 1.0e-9_wp, "uniform: shearU")
      call near(error, maxval(abs(fx%w%stresspT(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-9_wp, "uniform: sp")
      call near(error, maxval(abs(fx%w%stressmT(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-9_wp, "uniform: sm")
      call near(error, maxval(abs(fx%w%divu(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-18_wp, "uniform: divu")
      call near(error, maxval(abs(fx%w%shear(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-18_wp, "uniform: shear")
      call near(error, maxval(abs(fx%w%vort(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-18_wp, "uniform: vort")
      call near(error, maxval(abs(fx%w%strintxE(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-15_wp, "uniform: strintx")
      call near(error, maxval(abs(fx%w%strintyN(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 1.0e-15_wp, "uniform: strinty")
      call fx_destroy(fx)
   end subroutine test_uniform_velocity

   subroutine test_solid_body_rotation(error)
      !! u = -Omega y, v = Omega x: div = shear = 0, vorticity = 2 Omega.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j
      real(wp), parameter :: om = 1.0e-4_wp
      real(wp) :: wd, ws, wv

      call fx_build(fx, 9, 8, .false., 2.0e4_wp, 1.5e4_wp, 4)
      call fx_ice(fx, 0.8_wp, 1.5_wp)
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      call fill_linear(fx, 0.0_wp, -om, om, 0.0_wp, 0.0_wp, 0.0_wp)
      call fx_map(fx)
      call fx_prep(fx)
      call fx_avg(fx)
      call fx_fin(fx)
      call fx_pull(fx)
      wd = 0.0_wp
      ws = 0.0_wp
      wv = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            wd = max(wd, abs(fx%w%divu(i, j)))
            ws = max(ws, abs(fx%w%shear(i, j)))
            wv = max(wv, abs(fx%w%vort(i, j) - 2.0_wp*om))
         end do
      end do
      write (error_unit, '(a,3es10.2)') "  [evp_c] rotation: |divu|, |shear|, |vort-2W| = ", wd, ws, wv
      call near(error, wd, 0.0_wp, 0.0_wp, 1.0e-12_wp*om, "rotation divu")
      call near(error, ws, 0.0_wp, 0.0_wp, 1.0e-12_wp*om, "rotation shear")
      call near(error, wv, 0.0_wp, 0.0_wp, 1.0e-12_wp*om, "rotation vort")
      call fx_destroy(fx)
   end subroutine test_solid_body_rotation

   subroutine linear_case(error, tag, revised, cap_max, ktens, avg_strength, e_yc, e_pp)
      !! A linear velocity field in an open domain with uniform strength `P`:
      !! divu, shear, vort, rdg_conv, rdg_shear in closed form, and the first
      !! stress update of one subcycle (viscosities, replacement pressure).
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: tag
      logical, intent(in) :: revised, cap_max, avg_strength
      real(wp), intent(in) :: ktens, e_yc, e_pp
      type(fx_t) :: fx
      integer :: i0, j0, sgn
      real(wp), parameter :: a = 2.0e-6_wp, b = -1.0e-6_wp, c = 1.5e-6_wp, d = 0.5e-6_wp
      real(wp), parameter :: pval = 2.5e4_wp, sp_init = 123.0_wp, sm_init = -45.0_wp, s12_init = 7.0_wp
      real(wp) :: dd, tt, ss, area, ef, ep, dn, delta, dmin, tmp, zeta2, rep, eta2, keep, a1i, den
      real(wp) :: ex_sp, ex_sm, ex_s12, dpp2, rdg_s

      do sgn = 1, -1, -2
         call fx_build(fx, 10, 9, .false., 2.0e4_wp, 1.5e4_wp, 8)
         fx%p%dyn%revised_evp = revised
         fx%p%dyn%capping_max = cap_max
         fx%p%dyn%Ktens = ktens
         fx%p%dyn%e_yieldcurve = e_yc
         fx%p%dyn%e_plasticpot = e_pp
         if (avg_strength) fx%p%dyn%visc_method = ICE_VISC_AVG_STRENGTH
         call fx_consts(fx)
         call fx_ice(fx, 0.9_wp, 1.0_wp)
         fx%w%iceEmask = .true.
         fx%w%iceNmask = .true.
         call fill_linear(fx, real(sgn, wp)*a, real(sgn, wp)*b, real(sgn, wp)*c, real(sgn, wp)*d, 0.01_wp, -0.02_wp)
         fx%w%stresspT = sp_init
         fx%w%stressmT = sm_init
         fx%w%stress12U = s12_init
         call fx_map(fx)
         call fx_prep(fx)
         call fx_set_strength(fx, pval)
         if (avg_strength) call fx_strength_u(fx)
         ! prep zeroed the stresses off the ice mask only; re-impose the initial stresses
         fx%w%stresspT = sp_init
         fx%w%stressmT = sm_init
         fx%w%stress12U = s12_init
         call fx_push_w(fx)
         call fx_avg(fx)
         call fx_k123(fx)
         call fx_fin(fx)
         call fx_pull(fx)
         i0 = fx%is + 4
         j0 = fx%js + 3
         area = fx%dx*fx%dy
         dd = real(sgn, wp)*(a + d)
         tt = real(sgn, wp)*(a - d)
         ss = real(sgn, wp)*(b + c)
         ef = e_yc*e_yc/(e_pp**4)
         ep = 1.0_wp/(e_pp*e_pp)
         dn = sqrt(dd*dd + ef*(tt*tt + ss*ss))
         delta = area*dn
         dmin = fx%p%dyn%deltaminEVP*area
         if (cap_max) then
            tmp = pval/max(delta, dmin)
         else
            tmp = pval/(delta + dmin)
         end if
         zeta2 = (1.0_wp + ktens)*tmp
         rep = (1.0_wp - ktens)*tmp*delta
         eta2 = ep*zeta2
         keep = 1.0_wp - fx%k%arlx1i*fx%k%revp
         a1i = fx%k%arlx1i
         den = fx%k%denom1
         ex_sp = (sp_init*keep + a1i*(zeta2*area*dd - rep))*den
         ex_sm = (sm_init*keep + a1i*eta2*area*tt)*den
         ex_s12 = (s12_init*keep + a1i*0.5_wp*eta2*area*ss)*den
         call near(error, fx%w%zetax2T(i0, j0), zeta2, TOL_CF, 0.0_wp, tag//" zeta2")
         call near(error, fx%w%etax2T(i0, j0), eta2, TOL_CF, 0.0_wp, tag//" eta2")
         call near(error, fx%w%stresspT(i0, j0), ex_sp, TOL_CF, 0.0_wp, tag//" stresspT")
         call near(error, fx%w%stressmT(i0, j0), ex_sm, TOL_CF, 0.0_wp, tag//" stressmT")
         call near(error, fx%w%stress12U(i0 + 1, j0 + 1), ex_s12, TOL_CF, 0.0_wp, tag//" stress12U")
         ! diagnostics (closed forms in 1/s)
         dpp2 = dd*dd + ep*(tt*tt + ss*ss)
         rdg_s = 0.5_wp*(dpp2/dn - abs(dd))
         call near(error, fx%w%divu(i0, j0), dd, TOL_CF, 0.0_wp, tag//" divu")
         call near(error, fx%w%shear(i0, j0), sqrt(tt*tt + ss*ss), TOL_CF, 0.0_wp, tag//" shear")
         call near(error, fx%w%vort(i0, j0), real(sgn, wp)*(c - b), TOL_CF, 0.0_wp, tag//" vort")
         call near(error, fx%w%rdg_conv(i0, j0), max(-dd, 0.0_wp), TOL_CF, 0.0_wp, tag//" rdg_conv")
         call near(error, fx%w%rdg_shear(i0, j0), rdg_s, TOL_CF, 0.0_wp, tag//" rdg_shear")
         call fx_destroy(fx)
         if (allocated(error)) return
      end do
   end subroutine linear_case

   subroutine test_linear_classic(error)
      type(error_type), allocatable, intent(out) :: error
      call linear_case(error, "linear classic", .false., .true., 0.0_wp, .false., 2.0_wp, 2.5_wp)
   end subroutine test_linear_classic

   subroutine test_linear_sum_ktens(error)
      type(error_type), allocatable, intent(out) :: error
      call linear_case(error, "linear sum Ktens", .false., .false., 0.15_wp, .false., 2.0_wp, 2.0_wp)
   end subroutine test_linear_sum_ktens

   subroutine test_linear_revised(error)
      type(error_type), allocatable, intent(out) :: error
      call linear_case(error, "linear revised", .true., .true., 0.0_wp, .false., 1.8_wp, 2.2_wp)
   end subroutine test_linear_revised

   subroutine test_linear_avg_strength(error)
      type(error_type), allocatable, intent(out) :: error
      call linear_case(error, "linear avg_strength", .false., .true., 0.1_wp, .true., 2.0_wp, 2.0_wp)
   end subroutine test_linear_avg_strength

   subroutine test_div_stress_linear(error)
      !! Divergence of prescribed linear stresses:
      !!   fx = (1/2) d(sp)/dx + (1/2) d(sm)/dx + d(s12)/dy,
      !!   fy = (1/2) d(sp)/dy - (1/2) d(sm)/dy + d(s12)/dx.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j, ii, jj
      real(wp), parameter :: ap = 3.0e-3_wp, bp = -2.0e-3_wp, am = 1.0e-3_wp, bm = 4.0e-3_wp
      real(wp), parameter :: cx = 1.5e-3_wp, cy = -2.5e-3_wp
      real(wp) :: ex, ey, xc, yc, xu, yu

      call fx_build(fx, 10, 9, .false., 2.0e4_wp, 1.5e4_wp, 4)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      ! sp = ap x + bp y, sm = am x + bm y at the T centres; s12 = cx y + cy x at the corners
      do j = 1, fx%ny
         do i = 1, fx%nx
            xc = (real(i - fx%is, wp) + 0.5_wp)*fx%dx
            yc = (real(j - fx%js, wp) + 0.5_wp)*fx%dy
            fx%w%stresspT(i, j) = ap*xc + bp*yc
            fx%w%stressmT(i, j) = am*xc + bm*yc
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx + 1
            ! corner U(i-1,j-1) sits at x = (i-1-is+1) dx, y = (j-1-js+1) dy
            xu = real(i - 1 - fx%is + 1, wp)*fx%dx
            yu = real(j - 1 - fx%js + 1, wp)*fx%dy
            fx%w%stress12U(i, j) = cx*yu + cy*xu
         end do
      end do
      call fx_map(fx)
      call fx_prep(fx)
      fx%w%stresspT = 0.0_wp
      do j = 1, fx%ny
         do i = 1, fx%nx
            xc = (real(i - fx%is, wp) + 0.5_wp)*fx%dx
            yc = (real(j - fx%js, wp) + 0.5_wp)*fx%dy
            fx%w%stresspT(i, j) = ap*xc + bp*yc
            fx%w%stressmT(i, j) = am*xc + bm*yc
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx + 1
            xu = real(i - fx%is, wp)*fx%dx
            yu = real(j - fx%js, wp)*fx%dy
            fx%w%stress12U(i, j) = cx*yu + cy*xu
         end do
      end do
      call fx_push_w(fx)
      call ice_evp_c_div_stress(fx%w, fx%k, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%dxE, fx%dyE, fx%dxN, &
                                fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT, fx%earear, fx%narear)
      call fx_pull(fx)
      ex = 0.5_wp*ap + 0.5_wp*am + cx
      ey = 0.5_wp*bp - 0.5_wp*bm + cy
      do jj = fx%js, fx%je
         do ii = fx%is, fx%ie
            call near(error, fx%w%strintxE(ii, jj), ex, 1.0e-11_wp, 0.0_wp, "strintx")
            call near(error, fx%w%strintyN(ii, jj), ey, 1.0e-11_wp, 0.0_wp, "strinty")
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      ! free_drift_approx: below rheo_area_min the rheology is dropped
      if (.not. allocated(error)) then
         fx%aice_init = 5.0e-4_wp*fx%hm
         fx%p%dyn%free_drift_approx = .true.
         call fx_consts(fx)
         call fx_unmap(fx)
         call fx_map(fx)
         call fx_prep(fx)
         fx%w%stresspT = 0.0_wp
         fx%w%stresspT(:, :) = ap*real(1, wp)
         call fx_push_w(fx)
         call ice_evp_c_div_stress(fx%w, fx%k, fx%nx, fx%ny, fx%is, fx%ie, fx%js, fx%je, fx%dxE, fx%dyE, &
                                   fx%dxN, fx%dyN, fx%dxU, fx%dyU, fx%dxT, fx%dyT, fx%earear, fx%narear)
         call fx_pull(fx)
         call near(error, maxval(abs(fx%w%strintxE(fx%is:fx%ie, fx%js:fx%je))), 0.0_wp, 0.0_wp, 0.0_wp, &
                   "free drift: no internal stress below rheo_area_min")
      end if
      call fx_destroy(fx)
   end subroutine test_div_stress_linear

   subroutine free_drift_run(fx, f, uo, vo, taux, tauy, tx, ty, cosw, sinw, revised, niter)
      !! Iterate the momentum equation to its fixed point in a periodic open
      !! domain with no internal stress (strength 0); uniform fields.
      type(fx_t), intent(inout) :: fx
      real(wp), intent(in) :: f, uo, vo, taux, tauy, tx, ty, cosw, sinw
      logical, intent(in) :: revised
      integer, intent(in) :: niter
      integer :: it

      call fx_build(fx, 8, 8, .false., 2.0e4_wp, 1.5e4_wp, 2)
      fx%p%dyn%cosw = cosw
      fx%p%dyn%sinw = sinw
      fx%p%dyn%revised_evp = revised
      fx%p%dyn%arlx = 20.0_wp
      fx%p%dyn%brlx = 2.0_wp
      call fx_consts(fx)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      fx%fcorE = f; fx%fcorN = f
      fx%strairxE = taux; fx%strairyN = tauy
      fx%ss_tltxE = tx; fx%ss_tltyN = ty
      fx%uocnE = uo; fx%vocnE = vo; fx%uocnN = uo; fx%vocnN = vo
      call fx_map(fx)
      call fx_prep(fx)
      call fx_set_strength(fx, 0.0_wp)
      call fx_wrap(fx, .true., .true.)
      call fx_avg(fx)
      do it = 1, niter
         call fx_sub(fx)
         call fx_wrap(fx, .true., .true.)
         call fx_avg(fx)
      end do
      call fx_pull(fx)
   end subroutine free_drift_run

   subroutine test_free_drift_balance(error)
      !! Wind + sea-surface tilt + ocean drag (turning angle) + Coriolis:
      !! the converged velocity satisfies the 2x2 balance of the scheme's fixed
      !! point (classic EVP: no inertia; revised EVP: the implicit-Euler
      !! inertia `m (u - u_start)/dt`, Bouillon et al. 2013), and equals an
      !! independent complex fixed-point solution.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp), parameter :: f = 1.3e-4_wp, uo = 0.05_wp, vo = -0.03_wp
      real(wp), parameter :: taux = 0.12_wp, tauy = -0.05_wp, tx = 2.0e-5_wp, ty = -1.0e-5_wp
      real(wp) :: cosw, sinw, mass, g, vrel, fxx, fyy, rx, ry, ui, vi, ur, vr, cr, ci, nr, ni, den, scale, inert
      integer :: it, i0, j0, m
      logical :: revised
      character(len=20) :: tag

      cosw = cos(25.0_wp*PI/180.0_wp)
      sinw = sin(25.0_wp*PI/180.0_wp)
      do it = 1, 2
         revised = it == 2
         tag = merge("revised EVP", "classic EVP", revised)
         call free_drift_run(fx, f, uo, vo, taux, tauy, tx, ty, cosw, sinw, revised, 900)
         i0 = fx%is + 3
         j0 = fx%js + 3
         mass = fx%p%c%rhoi*0.9_wp*1.0_wp
         g = fx%p%c%gravit
         inert = merge(mass/fx%dt, 0.0_wp, revised)
         ui = fx%u(i0 + 1, j0)
         vi = fx%v(i0, j0 + 1)
         vrel = 0.9_wp*fx%p%c%rhow*fx%p%dyn%cdn_ocn*sqrt((uo - ui)**2 + (vo - vi)**2)
         fxx = taux - g*mass*tx
         fyy = tauy - g*mass*ty
         rx = inert*(ui - uo) + vrel*cosw*ui - (f*mass + vrel*sinw)*vi - (fxx + vrel*(uo*cosw - vo*sinw))
         ry = inert*(vi - vo) + vrel*cosw*vi + (f*mass + vrel*sinw)*ui - (fyy + vrel*(vo*cosw + uo*sinw))
         scale = abs(fxx) + abs(fyy)
         write (error_unit, '(a,a,2es10.2,a,2f10.6)') "  [evp_c] free drift ", trim(tag), rx/scale, ry/scale, &
            "  (u,v) = ", ui, vi
         call near(error, rx/scale, 0.0_wp, 0.0_wp, 1.0e-9_wp, trim(tag)//" x balance")
         call near(error, ry/scale, 0.0_wp, 0.0_wp, 1.0e-9_wp, trim(tag)//" y balance")
         ! independent solution of ((inert + vr cosw) + i (f m + vr sinw)) U
         !   = F + vr (cosw + i sinw)(uo + i vo) + inert (uo + i vo), vr = a rho_w cd |Uo - U|
         ur = 0.0_wp
         vr = 0.0_wp
         do m = 1, 4000
            vrel = 0.9_wp*fx%p%c%rhow*fx%p%dyn%cdn_ocn*sqrt((uo - ur)**2 + (vo - vr)**2)
            cr = inert + vrel*cosw
            ci = f*mass + vrel*sinw
            nr = fxx + vrel*(uo*cosw - vo*sinw) + inert*uo
            ni = fyy + vrel*(vo*cosw + uo*sinw) + inert*vo
            den = cr*cr + ci*ci
            ur = 0.5_wp*ur + 0.5_wp*((nr*cr + ni*ci)/den)
            vr = 0.5_wp*vr + 0.5_wp*((ni*cr - nr*ci)/den)
         end do
         call near(error, ui, ur, 1.0e-8_wp, 1.0e-10_wp, trim(tag)//" u vs complex fixed point")
         call near(error, vi, vr, 1.0e-8_wp, 1.0e-10_wp, trim(tag)//" v vs complex fixed point")
         call fx_destroy(fx)
         if (allocated(error)) return
      end do
   end subroutine test_free_drift_balance

   subroutine test_tilt_sign(error)
      !! Sea surface rising toward +x pushes the ice toward -x (force = -g m slope).
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp) :: mass, vrel, want

      call free_drift_run(fx, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0e-5_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
                          .false., 700)
      mass = fx%p%c%rhoi*0.9_wp
      call check(error, fx%u(fx%is + 3, fx%js + 3) < 0.0_wp, .true., "tilt pushes the ice downhill (u < 0)")
      if (allocated(error)) then
         call fx_destroy(fx)
         return
      end if
      ! with f = 0, no wind, no ocean current: vrel |u| = g m s  =>  0.9 rhow cd u^2 = g m s
      vrel = 0.9_wp*fx%p%c%rhow*fx%p%dyn%cdn_ocn
      want = -sqrt(fx%p%c%gravit*mass*3.0e-5_wp/vrel)
      call near(error, fx%u(fx%is + 3, fx%js + 3), want, 1.0e-8_wp, 0.0_wp, "tilt: quadratic-drag equilibrium")
      call near(error, fx%v(fx%is + 3, fx%js + 3), 0.0_wp, 0.0_wp, 1.0e-10_wp, "tilt: no v")
      call fx_destroy(fx)
   end subroutine test_tilt_sign

   subroutine random_state(fx, seed, with_island)
      !! Random ice, forcing and velocity on a (wall) fixture; deterministic.
      type(fx_t), intent(inout) :: fx
      integer(int64), intent(in) :: seed
      logical, intent(in) :: with_island
      integer(int64) :: s
      integer :: i, j
      real(wp) :: a, h, r

      s = lcg_next(seed)
      do j = 1, fx%ny
         do i = 1, fx%nx
            if (fx%hm(i, j) < 0.5_wp) cycle
            s = lcg_next(s); r = lcg_unit(s)
            if (r < 0.2_wp) then
               a = 0.0_wp
            else if (r < 0.3_wp) then
               a = 4.0e-4_wp           ! below rheo_area_min, above dyn_area_min: free drift
            else
               s = lcg_next(s); a = 0.3_wp + 0.7_wp*lcg_unit(s)
            end if
            s = lcg_next(s); h = 0.2_wp + 2.8_wp*lcg_unit(s)
            fx%aice(i, j) = a
            fx%aice_init(i, j) = a
            fx%vice(i, j) = a*h
            fx%aicen(i, j, 1) = a
            fx%vicen(i, j, 1) = a*h
            fx%aice0(i, j) = 1.0_wp - a
            s = lcg_next(s); fx%vsno(i, j) = a*0.1_wp*lcg_unit(s)
         end do
      end do
      if (with_island) then
         i = fx%is + fx%nxp/2
         j = fx%js + fx%nyp/2
         fx%hm(i:i + 1, j) = 0.0_wp
         call rebuild_masks(fx)
         fx%aice(i:i + 1, j) = 0.0_wp; fx%aice_init(i:i + 1, j) = 0.0_wp; fx%vice(i:i + 1, j) = 0.0_wp
         fx%aicen(i:i + 1, j, 1) = 0.0_wp; fx%vicen(i:i + 1, j, 1) = 0.0_wp
         fx%aice0(i:i + 1, j) = 1.0_wp; fx%vsno(i:i + 1, j) = 0.0_wp
      end if
      do j = 1, fx%ny
         do i = 1, fx%nx
            s = lcg_next(s); fx%fcorE(i, j) = 1.4e-4_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%fcorN(i, j) = 1.4e-4_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%strairxE(i, j) = 0.2_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%strairyN(i, j) = 0.2_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%ss_tltxE(i, j) = 2.0e-5_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%ss_tltyN(i, j) = 2.0e-5_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%uocnE(i, j) = 0.1_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%vocnE(i, j) = 0.1_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%uocnN(i, j) = 0.1_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%vocnN(i, j) = 0.1_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
         end do
      end do
      do j = 1, fx%ny
         do i = 1, fx%nx + 1
            s = lcg_next(s); fx%u(i, j) = 0.2_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx
            s = lcg_next(s); fx%v(i, j) = 0.2_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
         end do
      end do
      do j = 1, fx%ny
         do i = 1, fx%nx
            s = lcg_next(s); fx%w%stresspT(i, j) = 50.0_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%w%stressmT(i, j) = 50.0_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
            s = lcg_next(s); fx%w%zetax2T(i, j) = lcg_unit(s)
            s = lcg_next(s); fx%w%etax2T(i, j) = lcg_unit(s)
            s = lcg_next(s); fx%w%iceEmask(i, j) = lcg_unit(s) > 0.5_wp
            s = lcg_next(s); fx%w%iceNmask(i, j) = lcg_unit(s) > 0.5_wp
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx + 1
            s = lcg_next(s); fx%w%stress12U(i, j) = 50.0_wp*(2.0_wp*lcg_unit(s) - 1.0_wp)
         end do
      end do
   end subroutine random_state

   subroutine rebuild_masks(fx)
      !! `epm/npm/uvm` from `hm` (a closed box + any islands).
      type(fx_t), intent(inout) :: fx
      integer :: i, j, nx, ny
      nx = fx%nx
      ny = fx%ny
      fx%epm = 0.0_wp; fx%npm = 0.0_wp; fx%uvm = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            fx%epm(i, j) = fx%hm(i, j)*fx%hm(i + 1, j)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            fx%npm(i, j) = fx%hm(i, j)*fx%hm(i, j + 1)
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx - 1
            fx%uvm(i, j) = fx%hm(i, j)*fx%hm(i + 1, j)*fx%hm(i, j + 1)*fx%hm(i + 1, j + 1)
         end do
      end do
   end subroutine rebuild_masks

   ! ---- mirror symmetry --------------------------------------------------

   subroutine mirror_into(src, dst, axis)
      !! `dst` = the mirror image of `src`'s inputs (axis 1: x, 2: y).
      type(fx_t), intent(in) :: src
      type(fx_t), intent(inout) :: dst
      integer, intent(in) :: axis
      integer :: i, j, nx, ny, im, jm

      nx = src%nx
      ny = src%ny
      do j = 1, ny
         do i = 1, nx
            if (axis == 1) then
               im = nx + 1 - i; jm = j
            else
               im = i; jm = ny + 1 - j
            end if
            dst%aice(i, j) = src%aice(im, jm); dst%aice_init(i, j) = src%aice_init(im, jm)
            dst%vice(i, j) = src%vice(im, jm); dst%vsno(i, j) = src%vsno(im, jm)
            dst%aice0(i, j) = src%aice0(im, jm); dst%aicen(i, j, 1) = src%aicen(im, jm, 1)
            dst%vicen(i, j, 1) = src%vicen(im, jm, 1)
            ! N-type (index nx+1-i for x, ny-j for y) and E-type arrays
            if (axis == 1) then
               dst%strairyN(i, j) = src%strairyN(im, jm)
               dst%ss_tltyN(i, j) = src%ss_tltyN(im, jm)
               dst%uocnN(i, j) = -src%uocnN(im, jm)
               dst%vocnN(i, j) = src%vocnN(im, jm)
               if (i < nx) then
                  dst%strairxE(i, j) = -src%strairxE(nx - i, j)
                  dst%ss_tltxE(i, j) = -src%ss_tltxE(nx - i, j)
                  dst%uocnE(i, j) = -src%uocnE(nx - i, j)
                  dst%vocnE(i, j) = src%vocnE(nx - i, j)
               end if
            else
               dst%strairxE(i, j) = src%strairxE(im, jm)
               dst%ss_tltxE(i, j) = src%ss_tltxE(im, jm)
               dst%uocnE(i, j) = src%uocnE(im, jm)
               dst%vocnE(i, j) = -src%vocnE(im, jm)
               if (j < ny) then
                  dst%strairyN(i, j) = -src%strairyN(i, ny - j)
                  dst%ss_tltyN(i, j) = -src%ss_tltyN(i, ny - j)
                  dst%uocnN(i, j) = src%uocnN(i, ny - j)
                  dst%vocnN(i, j) = -src%vocnN(i, ny - j)
               end if
            end if
         end do
      end do
   end subroutine mirror_into

   subroutine mirror_case(error, axis)
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: axis
      type(fx_t) :: a, b
      integer :: i, j, nx, ny, nxp, nyp, k
      real(wp) :: worst, scale, e
      character(len=24) :: nm

      nxp = 11
      nyp = 9
      do k = 1, 2
         call fx_build(a, nxp, nyp, .true., 2.0e4_wp, 1.5e4_wp, 12)
         call fx_build(b, nxp, nyp, .true., 2.0e4_wp, 1.5e4_wp, 12)
         if (k == 2) then
            a%p%dyn%revised_evp = .true.
            b%p%dyn%revised_evp = .true.
            call fx_consts(a)
            call fx_consts(b)
         end if
         nx = a%nx
         ny = a%ny
         call random_state(a, 20240307_int64, .false.)
         ! f = 0, no turning, zero start (symmetric); keep the random thickness / forcing
         a%fcorE = 0.0_wp; a%fcorN = 0.0_wp
         a%u = 0.0_wp; a%v = 0.0_wp
         a%w%stresspT = 0.0_wp; a%w%stressmT = 0.0_wp; a%w%zetax2T = 0.0_wp; a%w%etax2T = 0.0_wp
         a%w%stress12U = 0.0_wp; a%w%iceEmask = .false.; a%w%iceNmask = .false.
         call mirror_into(a, b, axis)
         b%u = 0.0_wp; b%v = 0.0_wp
         call fx_map(a); call fx_map(b)
         call fx_step(a, .false., .false.)
         call fx_step(b, .false., .false.)
         call fx_pull(a); call fx_pull(b)
         worst = 0.0_wp
         scale = 0.0_wp
         do j = a%js, a%je + 1
            do i = a%is, a%ie + 1
               if (axis == 1) then
                  e = abs(b%u(i, j) + a%u(nx + 2 - i, j))
               else
                  e = abs(b%u(i, j) - a%u(i, ny + 1 - j))
               end if
               if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
               worst = max(worst, e)
               scale = max(scale, abs(a%u(i, j)))
            end do
         end do
         do j = a%js, a%je + 1
            do i = a%is, a%ie + 1
               if (axis == 1) then
                  e = abs(b%v(i, j) - a%v(nx + 1 - i, j))
               else
                  e = abs(b%v(i, j) + a%v(i, ny + 2 - j))
               end if
               if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
               worst = max(worst, e)
               scale = max(scale, abs(a%v(i, j)))
            end do
         end do
         do j = a%js, a%je
            do i = a%is, a%ie
               if (axis == 1) then
                  e = abs(b%w%stresspT(i, j) - a%w%stresspT(nx + 1 - i, j))
                  e = max(e, abs(b%w%stressmT(i, j) - a%w%stressmT(nx + 1 - i, j)))
                  e = max(e, abs(b%w%divu(i, j) - a%w%divu(nx + 1 - i, j)))
               else
                  e = abs(b%w%stresspT(i, j) - a%w%stresspT(i, ny + 1 - j))
                  e = max(e, abs(b%w%stressmT(i, j) - a%w%stressmT(i, ny + 1 - j)))
                  e = max(e, abs(b%w%divu(i, j) - a%w%divu(i, ny + 1 - j)))
               end if
               if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
               worst = max(worst, e)
            end do
         end do
         write (nm, '(a,i0,a,i0)') "mirror axis ", axis, " rev ", k - 1
         write (error_unit, '(a,a,a,es10.2,a,es10.2)') "  [evp_c] ", trim(nm), " worst ", worst, " |u|max ", scale
         call check(error, scale > 1.0e-3_wp, .true., "mirror case moved (non-trivial)")
         if (.not. allocated(error)) call near(error, worst, 0.0_wp, 0.0_wp, 1.0e-9_wp*scale, trim(nm))
         call fx_destroy(a); call fx_destroy(b)
         if (allocated(error)) return
      end do
   end subroutine mirror_case

   subroutine test_mirror_x(error)
      type(error_type), allocatable, intent(out) :: error
      call mirror_case(error, 1)
   end subroutine test_mirror_x

   subroutine test_mirror_y(error)
      type(error_type), allocatable, intent(out) :: error
      call mirror_case(error, 2)
   end subroutine test_mirror_y

   subroutine test_periodic_shift(error)
      !! x-periodic channel (walls north and south): evaluating every periodic
      !! field one cell further along gives the solution one cell further
      !! along (checks the wrap and the translation invariance of every
      !! stencil).
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: a, b
      integer :: i, j, nxp, nyp
      real(wp) :: worst, scale, e
      real(wp), allocatable :: ra(:), rb(:), rc(:), rd(:), re(:)
      integer(int64) :: s

      nxp = 12
      nyp = 8
      allocate (ra(nxp), rb(nxp), rc(nxp), rd(nxp), re(nxp))
      s = lcg_next(777_int64)
      do i = 1, nxp
         s = lcg_next(s); ra(i) = 0.4_wp + 0.6_wp*lcg_unit(s)
         s = lcg_next(s); rb(i) = 0.3_wp + 1.5_wp*lcg_unit(s)
         s = lcg_next(s); rc(i) = 0.2_wp*(lcg_unit(s) - 0.5_wp)
         s = lcg_next(s); rd(i) = 0.2_wp*(lcg_unit(s) - 0.5_wp)
         s = lcg_next(s); re(i) = 0.1_wp*(lcg_unit(s) - 0.5_wp)
      end do
      call fx_build(a, nxp, nyp, .false., 2.0e4_wp, 1.5e4_wp, 10)
      call fx_build(b, nxp, nyp, .false., 2.0e4_wp, 1.5e4_wp, 10)
      call periodic_fill(a, 0)
      call periodic_fill(b, 1)
      call fx_map(a); call fx_map(b)
      call fx_step(a, .true., .false.)
      call fx_step(b, .true., .false.)
      call fx_pull(a); call fx_pull(b)
      worst = 0.0_wp
      scale = 0.0_wp
      do j = a%js, a%je
         do i = a%is + 1, a%ie
            e = abs(b%u(i, j) - a%u(i + 1, j))
            if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
            worst = max(worst, e)
            scale = max(scale, abs(a%u(i + 1, j)))
            e = abs(b%v(i, j + 1) - a%v(i + 1, j + 1))
            if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
            worst = max(worst, e)
            e = abs(b%w%stresspT(i, j) - a%w%stresspT(i + 1, j))*1.0e-4_wp
            if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
            worst = max(worst, e)
         end do
      end do
      write (error_unit, '(a,es10.2,a,es10.2)') "  [evp_c] periodic shift worst ", worst, " |u|max ", scale
      call check(error, scale > 1.0e-4_wp, .true., "channel moved")
      if (.not. allocated(error)) call near(error, worst, 0.0_wp, 0.0_wp, 1.0e-9_wp*scale, "periodic shift")
      call fx_destroy(a); call fx_destroy(b)
   contains

      integer function pidx(i, sh)
         integer, intent(in) :: i, sh
         pidx = modulo(i - a%is + sh, nxp) + 1
      end function pidx

      subroutine periodic_fill(fxx, sh)
         type(fx_t), intent(inout) :: fxx
         integer, intent(in) :: sh
         integer :: ii, jj
         fxx%hm = 1.0_wp
         do jj = 1, fxx%ny
            if (jj < fxx%js .or. jj > fxx%je) fxx%hm(:, jj) = 0.0_wp
         end do
         call rebuild_masks(fxx)
         do jj = fxx%js, fxx%je
            do ii = 1, fxx%nx
               fxx%aice(ii, jj) = ra(pidx(ii, sh))
               fxx%aice_init(ii, jj) = fxx%aice(ii, jj)
               fxx%vice(ii, jj) = ra(pidx(ii, sh))*rb(pidx(ii, sh))
               fxx%aicen(ii, jj, 1) = fxx%aice(ii, jj)
               fxx%vicen(ii, jj, 1) = fxx%vice(ii, jj)
               fxx%aice0(ii, jj) = 1.0_wp - fxx%aice(ii, jj)
               fxx%strairxE(ii, jj) = rd(pidx(ii + 1, sh))
               fxx%strairyN(ii, jj) = re(pidx(ii, sh))
               fxx%fcorE(ii, jj) = 1.2e-4_wp
               fxx%fcorN(ii, jj) = 1.2e-4_wp
            end do
            do ii = 1, fxx%nx + 1
               fxx%u(ii, jj) = rc(pidx(ii, sh))*(1.0_wp + 0.1_wp*real(jj, wp))
            end do
         end do
         do jj = fxx%js + 1, fxx%je
            do ii = 1, fxx%nx
               fxx%v(ii, jj) = re(pidx(ii, sh))*(1.0_wp - 0.05_wp*real(jj, wp))
            end do
         end do
         fxx%w%iceEmask = .true.
         fxx%w%iceNmask = .true.
      end subroutine periodic_fill

   end subroutine test_periodic_shift

   subroutine test_closed_basin(error)
      !! Onshore wind piles the ice against the east wall: no velocity through
      !! a wall face, bounded stress, convergence at the wall.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: istep, i, j
      real(wp) :: smax, pval, ucore, cmax

      call fx_build(fx, 16, 10, .true., 2.0e4_wp, 2.0e4_wp, 60)
      call fx_ice(fx, 0.95_wp, 1.5_wp)
      fx%strairxE = 0.15_wp
      fx%fcorE = 0.0_wp; fx%fcorN = 0.0_wp
      call fx_map(fx)
      smax = 0.0_wp
      do istep = 1, 12
         call fx_step(fx, .false., .false.)
      end do
      call fx_pull(fx)
      pval = fx%p%ridge%Pstar*1.5_wp*0.95_wp*exp(-fx%p%ridge%Cstar*(1.0_wp - 0.95_wp))
      do j = 1, fx%ny
         do i = 1, fx%nx
            smax = max(smax, abs(fx%w%stresspT(i, j)))
         end do
      end do
      ucore = fx%u(fx%is + 8, fx%js + 5)
      cmax = maxval(fx%w%rdg_conv(fx%is:fx%ie, fx%js:fx%je))
      write (error_unit, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "  [evp_c] closed basin: max|sp| ", smax, &
         " P ", pval, " u(mid) ", ucore, " max rdg_conv ", cmax
      call check(error, all(fx%u(fx%is, :) == 0.0_wp) .and. all(fx%u(fx%ie + 1, :) == 0.0_wp), .true., &
                 "no flow through the west/east wall")
      if (allocated(error)) return
      call check(error, all(fx%v(:, fx%js) == 0.0_wp) .and. all(fx%v(:, fx%je + 1) == 0.0_wp), .true., &
                 "no flow through the south/north wall")
      if (allocated(error)) return
      call check(error, ucore > 0.0_wp, .true., "ice is pushed toward the east")
      if (allocated(error)) return
      call check(error, cmax > 0.0_wp, .true., "convergence against the wall")
      if (allocated(error)) return
      call check(error, all(fx%ice_err == 0_int32), .true., "no error bits")
      if (allocated(error)) return
      ! sigma_1 = sigma_11 + sigma_22 = -2 P at the pure-convergence vertex of the yield curve
      call near(error, smax, 0.0_wp, 0.0_wp, 2.0_wp*pval, "|stresspT| bounded by 2 P (sigma_11 = sigma_22 = -P)")
      call check(error, smax > 0.5_wp*pval, .true., "the ice is actually compressed")
      call fx_destroy(fx)
   end subroutine test_closed_basin

   subroutine test_ndte_convergence(error)
      !! The solution after one step approaches a long-subcycle reference as
      !! ndte grows (classic and revised EVP); the residual decreases.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns(5) = [10, 20, 40, 80, 160]
      type(fx_t) :: fx
      real(wp) :: res(5, 2), ref(2), uref(8, 6), diff
      integer :: im, jn, i, j
      logical :: revised

      do im = 1, 2
         revised = im == 2
         do jn = 0, 5
            call fx_build(fx, 8, 6, .true., 2.0e4_wp, 2.0e4_wp, merge(1280, ns(max(jn, 1)), jn == 0))
            fx%p%dyn%revised_evp = revised
            fx%p%dyn%arlx = 300.0_wp
            fx%p%dyn%brlx = 300.0_wp
            call fx_consts(fx)
            call fx_ice(fx, 0.9_wp, 1.0_wp)
            fx%strairxE = 0.2_wp
            fx%fcorE = 1.0e-4_wp; fx%fcorN = 1.0e-4_wp
            call fx_map(fx)
            call fx_step(fx, .false., .false.)
            call fx_pull(fx)
            if (jn == 0) then
               uref = fx%u(fx%is + 1:fx%is + 8, fx%js:fx%js + 5)
            else
               diff = 0.0_wp
               do j = 1, 6
                  do i = 1, 8
                     diff = max(diff, abs(fx%u(fx%is + i, fx%js + j - 1) - uref(i, j)))
                  end do
               end do
               res(jn, im) = diff
            end if
            call fx_destroy(fx)
         end do
         ref(im) = maxval(abs(uref))
         write (error_unit, '(a,l1,a,5es10.2,a,es10.2)') "  [evp_c] ndte convergence revised=", revised, &
            " residuals ", res(:, im), "  |u|max ", ref(im)
         if (revised) then
            ! revised EVP: the iteration count, not dt/ndte, controls convergence; 10 -> 160 must improve
            call check(error, res(5, im) < res(1, im), .true., "revised: residual at 160 below that at 10")
         else
            call check(error, all(res(2:5, im) < res(1:4, im)), .true., "classic: residual decreases monotonically")
         end if
         if (allocated(error)) return
         call check(error, ref(im) > 1.0e-3_wp, .true., "the reference moved")
         if (allocated(error)) return
      end do
   end subroutine test_ndte_convergence

   subroutine test_zero_ice(error)
      !! Ice in the west half only: the ice-free east half keeps zero velocity,
      !! zero stress and zero deformation; garbage stresses there are reset.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j, imid, nbad

      call fx_build(fx, 16, 8, .true., 2.0e4_wp, 2.0e4_wp, 10)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      imid = fx%is + 7
      fx%aice(imid + 1:, :) = 0.0_wp; fx%aice_init(imid + 1:, :) = 0.0_wp; fx%vice(imid + 1:, :) = 0.0_wp
      fx%aicen(imid + 1:, :, 1) = 0.0_wp; fx%vicen(imid + 1:, :, 1) = 0.0_wp; fx%aice0(imid + 1:, :) = 1.0_wp
      fx%strairxE = 0.1_wp
      fx%w%stresspT = 9.0_wp; fx%w%stressmT = -9.0_wp; fx%w%stress12U = 5.0_wp
      fx%u = 0.3_wp; fx%v = 0.3_wp
      fx%w%iceEmask = .true.; fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_step(fx, .false., .false.)
      call fx_pull(fx)
      nbad = 0
      do j = fx%js, fx%je
         do i = imid + 3, fx%ie                 ! beyond the ice edge + the dilation ring
            if (fx%u(i + 1, j) /= 0.0_wp .or. fx%v(i, j + 1) /= 0.0_wp) nbad = nbad + 1
            if (fx%w%stresspT(i, j) /= 0.0_wp .or. fx%w%stressmT(i, j) /= 0.0_wp) nbad = nbad + 1
            if (fx%w%stress12U(i + 1, j + 1) /= 0.0_wp) nbad = nbad + 1
            if (fx%w%divu(i, j) /= 0.0_wp .or. fx%w%shear(i, j) /= 0.0_wp .or. fx%w%vort(i, j) /= 0.0_wp .or. &
                fx%w%rdg_conv(i, j) /= 0.0_wp .or. fx%w%rdg_shear(i, j) /= 0.0_wp) nbad = nbad + 1
            if (fx%w%strocnxE(i, j) /= 0.0_wp .or. fx%w%strocnyN(i, j) /= 0.0_wp) nbad = nbad + 1
         end do
      end do
      call check(error, nbad, 0, "ice-free cells carry no velocity, stress, deformation or ocean stress")
      if (allocated(error)) return
      call check(error, any(abs(fx%u(fx%is + 1:imid, fx%js:fx%je)) > 0.0_wp), .true., "the ice half moves")
      call fx_destroy(fx)
   end subroutine test_zero_ice

   subroutine test_nan_strength(error)
      !! A NaN strength in one cell is flagged (ICE_ERR_DYN_VEL) and the
      !! velocity near it is non-finite; nothing finite and plausible comes out.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i0, j0
      logical :: any_nan

      call fx_build(fx, 12, 10, .true., 2.0e4_wp, 2.0e4_wp, 6)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      fx%strairxE = 0.1_wp
      i0 = fx%is + 5
      j0 = fx%js + 4
      fx%aice(i0, j0) = ieee_value(1.0_wp, ieee_quiet_nan)
      call fx_map(fx)
      call fx_step(fx, .false., .false.)
      call fx_pull(fx)
      any_nan = any(.not. ieee_is_finite(fx%u(fx%is:fx%ie + 1, fx%js:fx%je))) .or. &
                any(.not. ieee_is_finite(fx%v(fx%is:fx%ie, fx%js:fx%je + 1)))
      call check(error, any(iand(fx%ice_err, ICE_ERR_DYN_VEL) /= 0_int32), .true., "NaN input raises ICE_ERR_DYN_VEL")
      if (allocated(error)) return
      call check(error, any_nan, .true., "NaN is not turned into a finite velocity")
      call fx_destroy(fx)
   end subroutine test_nan_strength

   subroutine test_nan_velocity(error)
      !! A NaN face velocity entering the step is flagged; a velocity above
      !! `vel_cap` is flagged without being clamped.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i0, j0, ii, jj, nlaund, nnan

      call fx_build(fx, 12, 10, .true., 2.0e4_wp, 2.0e4_wp, 6)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      i0 = fx%is + 5
      j0 = fx%js + 4
      fx%u(i0, j0) = ieee_value(1.0_wp, ieee_quiet_nan)
      fx%w%iceEmask = .true.; fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_step(fx, .false., .false.)
      call fx_pull(fx)
      call check(error, any(iand(fx%ice_err, ICE_ERR_DYN_VEL) /= 0_int32), .true., "NaN velocity raises ICE_ERR_DYN_VEL")
      if (.not. allocated(error)) then
         ! ridging reads rdg_conv: a non-finite divergence must never give a finite rdg_conv
         nlaund = 0
         nnan = 0
         do jj = fx%js, fx%je
            do ii = fx%is, fx%ie
               if (.not. ieee_is_finite(fx%w%divu(ii, jj))) then
                  nnan = nnan + 1
                  if (ieee_is_finite(fx%w%rdg_conv(ii, jj))) nlaund = nlaund + 1
               end if
            end do
         end do
         call check(error, nnan > 0 .and. nlaund == 0, .true., "NaN divu is not laundered into a finite rdg_conv")
      end if
      call fx_destroy(fx)
      if (allocated(error)) return
      ! over the cap: the wind alone drives a thin ice sheet past 0.5 m/s
      call fx_build(fx, 12, 10, .false., 2.0e4_wp, 2.0e4_wp, 6)
      fx%p%dyn%vel_cap = 0.05_wp
      call fx_consts(fx)
      call fx_ice(fx, 0.9_wp, 0.5_wp)
      fx%strairxE = 0.3_wp
      call fx_map(fx)
      call fx_step(fx, .true., .true.)
      call fx_pull(fx)
      call check(error, any(iand(fx%ice_err, ICE_ERR_DYN_VEL) /= 0_int32) .and. &
                 maxval(abs(fx%u(fx%is:fx%ie + 1, fx%js:fx%je))) > 0.05_wp, .true., "over-cap speed flagged, not clamped")
      call fx_destroy(fx)
   end subroutine test_nan_velocity

   subroutine test_python_transcript(error)
      !! One subcycle on an open 8x7 domain against the Python transcript
      !! (evp_c_transcript.py): the T stress, the corner stress, the stress
      !! divergence and the new u at one point.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j, i0, j0
      real(wp), parameter :: uo = 0.05_wp, vo = -0.02_wp, uon = 0.03_wp, von = 0.04_wp
      real(wp), parameter :: py_sp = -10264.906866703866_wp, py_sm = -300.2983875877536_wp
      real(wp), parameter :: py_eta = 6.796367424627854_wp, py_s12 = -39.704514032527676_wp
      real(wp), parameter :: py_unew = 0.044395397682571315_wp

      call fx_build(fx, 8, 7, .false., 2.0e4_wp, 1.5e4_wp, 4)
      call fx_ice(fx, 0.9_wp, 1.2_wp)
      fx%strairxE = 0.07_wp
      fx%strairyN = -0.03_wp
      fx%fcorE = 1.4e-4_wp
      fx%fcorN = 1.4e-4_wp
      fx%uocnE = uo; fx%vocnE = vo; fx%uocnN = uon; fx%vocnN = von
      do j = 1, fx%ny
         do i = 1, fx%nx + 1
            fx%u(i, j) = 0.05_wp*cos(0.7_wp*real(i, wp) + 0.3_wp*real(j, wp)) + 0.02_wp
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx
            fx%v(i, j) = 0.04_wp*sin(0.5_wp*real(i, wp) - 0.4_wp*real(j, wp)) - 0.01_wp
         end do
      end do
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_prep(fx)
      call fx_set_strength(fx, 2.0e4_wp)
      call fx_avg(fx)
      call fx_sub(fx)
      call fx_pull(fx)
      i0 = 7
      j0 = 7
      call near(error, fx%w%stresspT(i0, j0), py_sp, 1.0e-11_wp, 0.0_wp, "python stresspT")
      call near(error, fx%w%stressmT(i0, j0), py_sm, 1.0e-11_wp, 0.0_wp, "python stressmT")
      call near(error, fx%w%etax2T(i0, j0), py_eta, 1.0e-11_wp, 0.0_wp, "python etax2T")
      call near(error, fx%w%stress12U(i0 + 1, j0 + 1), py_s12, 1.0e-11_wp, 0.0_wp, "python stress12U")
      ! the new u uses the OLD tangential average (the substep reads w%vvelE from the previous K5)
      call near(error, fx%u(i0 + 1, j0), py_unew, 1.0e-10_wp, 0.0_wp, "python u at E(7,7)")
      call fx_destroy(fx)
   end subroutine test_python_transcript

   subroutine test_coastal_corner(error)
      !! Wall channel (walls south and north), uniform along-wall velocity U0:
      !! the coastal corner (two wet T cells) carries the no-slip shear
      !! `+/- 2 dxU U0` (the mirrored tangential value is -U0) and a non-zero
      !! corner stress, equal to the Python transcript's number.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j, jj
      real(wp), parameter :: u0 = 0.2_wp, dxx = 2.0e4_wp
      real(wp), parameter :: py_s12 = 1822.440157697287_wp, py_eta = 1.7677669529663687_wp

      call fx_build(fx, 10, 8, .false., dxx, 1.5e4_wp, 4)
      do jj = 1, fx%ny
         if (jj < fx%js .or. jj > fx%je) fx%hm(:, jj) = 0.0_wp
      end do
      call rebuild_masks(fx)
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      fx%u = 0.0_wp
      fx%v = 0.0_wp
      do j = fx%js, fx%je
         fx%u(:, j) = u0
      end do
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_prep(fx)
      call fx_set_strength(fx, 2.0e4_wp)
      call fx_avg(fx)
      call fx_k123(fx)
      call fx_pull(fx)
      i = fx%is + 4
      call check(error, fx%w%iceUmask(i, fx%js - 1) .and. fx%w%iceUmask(i, fx%je), .true., "coastal corners are stress points")
      call near(error, fx%w%shearU(i, fx%js - 1), 2.0_wp*dxx*u0, 1.0e-12_wp, 0.0_wp, "south wall corner shear = +2 dxU U0")
      call near(error, fx%w%shearU(i, fx%je), -2.0_wp*dxx*u0, 1.0e-12_wp, 0.0_wp, "north wall corner shear = -2 dxU U0")
      call near(error, fx%w%shearU(i, fx%js + 2), 0.0_wp, 0.0_wp, 1.0e-9_wp, "interior corner shear is zero")
      call near(error, fx%w%etax2T(i, fx%js), py_eta, 1.0e-11_wp, 0.0_wp, "wall-row T viscosity (python)")
      call near(error, fx%w%stress12U(i + 1, fx%js), py_s12, 1.0e-11_wp, 0.0_wp, "south wall corner stress (python)")
      call near(error, fx%w%stress12U(i + 1, fx%je + 1), -py_s12, 1.0e-11_wp, 0.0_wp, "north wall corner stress")
      call fx_destroy(fx)
   end subroutine test_coastal_corner

   subroutine spherical_error(nn, e_shear, e_vort, e_div)
      !! Rigid rotation about the pole on a lon-lat sector of `nn x nn` cells
      !! (u = Omega R cos(phi), v = 0): normalised residuals of the shear
      !! diagnostic and of vorticity - 2 Omega sin(phi), and the divergence.
      integer, intent(in) :: nn
      real(wp), intent(out) :: e_shear, e_vort, e_div
      type(fx_t) :: fx
      real(wp), parameter :: rr = 6.371e6_wp, om = 7.292e-5_wp, phi0 = 0.3_wp, span = 0.5_wp
      real(wp) :: h, pc, pn, ps, pnn
      integer :: i, j, ic, jc

      h = span/real(nn, wp)
      call fx_build(fx, nn, nn, .false., 1.0_wp, 1.0_wp, 4)
      do j = 1, fx%ny
         pc = phi0 + (real(j - fx%js, wp) + 0.5_wp)*h
         pn = pc + 0.5_wp*h
         ps = pc - 0.5_wp*h
         pnn = pc + h
         do i = 1, fx%nx
            fx%dxT(i, j) = rr*cos(pc)*h
            fx%dxE(i, j) = rr*cos(pc)*h
            fx%dxN(i, j) = rr*cos(pn)*h
            fx%dxU(i, j) = rr*cos(pn)*h
            fx%dyT(i, j) = rr*h; fx%dyE(i, j) = rr*h; fx%dyN(i, j) = rr*h; fx%dyU(i, j) = rr*h
            fx%tarea(i, j) = rr*rr*h*(sin(pn) - sin(ps))
            fx%earea(i, j) = fx%tarea(i, j)
            fx%narea(i, j) = rr*rr*h*(sin(pnn) - sin(pc))
            fx%uarea(i, j) = fx%narea(i, j)
            fx%tarear(i, j) = 1.0_wp/fx%tarea(i, j)
            fx%earear(i, j) = 1.0_wp/fx%earea(i, j)
            fx%narear(i, j) = 1.0_wp/fx%narea(i, j)
         end do
      end do
      call fx_ice(fx, 0.9_wp, 1.0_wp)
      do j = 1, fx%ny
         pc = phi0 + (real(j - fx%js, wp) + 0.5_wp)*h
         fx%u(:, j) = om*rr*cos(pc)
      end do
      fx%v = 0.0_wp
      fx%w%iceEmask = .true.
      fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_prep(fx)
      call fx_set_strength(fx, 1.0e4_wp)
      call fx_avg(fx)
      call fx_k123(fx)
      call fx_fin(fx)
      call fx_pull(fx)
      ic = fx%is + nn/2
      jc = fx%js + nn/2
      pc = phi0 + (real(jc - fx%js, wp) + 0.5_wp)*h
      e_shear = abs(fx%w%shear(ic, jc))/om
      e_vort = abs(fx%w%vort(ic, jc) - 2.0_wp*om*sin(pc))/om
      e_div = abs(fx%w%divu(ic, jc))/om
      call fx_destroy(fx)
   end subroutine spherical_error

   subroutine test_spherical_rotation(error)
      !! Metrics that VARY in space (a lon-lat sector): the rigid rotation has
      !! zero divergence exactly and zero shear / vorticity 2 Omega sin(phi)
      !! to second order, so the error drops by ~4 per halving of the cell
      !! size.  A swapped dxT / dxN / dxE or a wrong area weight would not.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: es(3), ev(3), ed(3), r1, r2, r3, r4
      integer :: k

      do k = 1, 3
         call spherical_error(8*2**(k - 1), es(k), ev(k), ed(k))
      end do
      r1 = es(1)/es(2); r2 = es(2)/es(3); r3 = ev(1)/ev(2); r4 = ev(2)/ev(3)
      write (error_unit, '(a,3es10.2,a,2f7.3)') "  [evp_c] spherical shear err/Omega (N=8,16,32) ", es, " ratios ", r1, r2
      write (error_unit, '(a,3es10.2,a,2f7.3)') "  [evp_c] spherical vort err/Omega ", ev, " ratios ", r3, r4
      call check(error, maxval(ed) < 1.0e-12_wp, .true., "divergence of the rigid rotation vanishes")
      if (allocated(error)) return
      call check(error, r1 > 3.5_wp .and. r1 < 4.5_wp .and. r2 > 3.5_wp .and. r2 < 4.5_wp, .true., "shear residual is second order")
      if (allocated(error)) return
  call check(error, r3 > 3.5_wp .and. r3 < 4.5_wp .and. r4 > 3.5_wp .and. r4 < 4.5_wp, .true., "vorticity residual is second order")
      if (allocated(error)) return
      call check(error, ev(3) < 1.0e-3_wp .and. es(3) < 1.0e-3_wp, .true., "residuals are small against Omega")
   end subroutine test_spherical_rotation

   ! ------------------------------------------------------------------
   ! scalar re-transcription on ice-name arrays
   ! ------------------------------------------------------------------

   subroutine ref_run(fx, u0, v0, iceE0, iceN0, sp0, sm0, s120, uout, vout, spo, smo, zo, eo, s12o, divu, shear, vort, &
                      rc, rs, sox, soy, iceT, iceE, iceN, iceU)
      !! Straight-loop transcription of one dynamics step (prep, average,
      !! `ndte` subcycles, finish) on ICE-NAME arrays: `uE(i,j) = u(i+1,j)`,
      !! `vN(i,j) = v(i,j+1)`, corner stress at `(i,j)`.  Independent of the
      !! library's kernels (own loops, own index shifts).  The strength is the
      !! library's (`fx%w%strength`, verified by the ridge tests).
      type(fx_t), intent(in) :: fx
      real(wp), intent(in) :: u0(:, :), v0(:, :), sp0(:, :), sm0(:, :), s120(:, :)
      logical, intent(in) :: iceE0(:, :), iceN0(:, :)
      real(wp), intent(out) :: uout(:, :), vout(:, :), spo(:, :), smo(:, :), zo(:, :), eo(:, :), s12o(:, :)
      real(wp), intent(out) :: divu(:, :), shear(:, :), vort(:, :), rc(:, :), rs(:, :), sox(:, :), soy(:, :)
      logical, intent(out) :: iceT(:, :), iceE(:, :), iceN(:, :), iceU(:, :)
      integer :: nx, ny, is, ie, js, je, i, j, ks
      real(wp), allocatable :: uE(:, :), vN(:, :), tmass(:, :), uvN(:, :), vvE(:, :), uvU(:, :), vvU(:, :), shU(:, :), &
                               dlU(:, :), aE(:, :), aN(:, :), mE(:, :), mN(:, :), fmE(:, :), fmN(:, :), fxE(:, :), &
                               fyN(:, :), ui(:, :), vi(:, :), sU(:, :)
      logical, allocatable :: tph(:, :)
      real(wp) :: w1, w2, w3, w4, wt, m, a, k_arlx1i, revp, denom1, brlx, ep, ef, kt, dmin, rho_i, rho_s, rho_w, g, cw
      real(wp) :: cosw, sinw, dtl, dx, dy, ar, divT, tenT, shrT, ssq, dlt, zeta, eta, rep, tmp, keep, vrel, sg
      real(wp) :: strint, uold, vold, cca, ccb, cc, rheo, vcap, dpp2, dv
      logical :: capmax, avgs, fd
      type(ice_params_t) :: p

      p = fx%p
      nx = fx%nx; ny = fx%ny; is = fx%is; ie = fx%ie; js = fx%js; je = fx%je
      dx = fx%dx; dy = fx%dy; ar = dx*dy
      allocate (uE(nx, ny), vN(nx, ny), tmass(nx, ny), uvN(nx, ny), vvE(nx, ny), uvU(nx, ny), vvU(nx, ny), shU(nx, ny), &
                dlU(nx, ny), aE(nx, ny), aN(nx, ny), mE(nx, ny), mN(nx, ny), fmE(nx, ny), fmN(nx, ny), fxE(nx, ny), &
                fyN(nx, ny), ui(nx, ny), vi(nx, ny), sU(nx, ny), tph(nx, ny))
      k_arlx1i = fx%k%arlx1i; revp = fx%k%revp; denom1 = fx%k%denom1; brlx = fx%k%brlx
      ep = fx%k%epp2i; ef = fx%k%e_factor; kt = fx%k%Ktens; dmin = fx%k%deltamin*ar
      rho_i = p%c%rhoi; rho_s = p%c%rhos; rho_w = p%c%rhow; g = p%c%gravit
      cw = p%dyn%cdn_ocn; cosw = p%dyn%cosw; sinw = p%dyn%sinw; dtl = fx%dt; vcap = p%dyn%vel_cap
      capmax = p%dyn%capping_max; avgs = p%dyn%visc_method == ICE_VISC_AVG_STRENGTH; fd = p%dyn%free_drift_approx
      keep = 1.0_wp - k_arlx1i*revp
      do j = 1, ny
         do i = 1, nx
            uE(i, j) = u0(i + 1, j)
            vN(i, j) = v0(i, j + 1)
         end do
      end do
      spo = sp0; smo = sm0; zo = 0.0_wp; eo = 0.0_wp
      ! ---- T masks ----
      tmass = 0.0_wp
      tph = .false.
      do j = 1, ny
         do i = 1, nx
            if (fx%hm(i, j) > 0.5_wp) then
               tmass(i, j) = rho_i*fx%vice(i, j) + rho_s*fx%vsno(i, j)
               tph(i, j) = fx%aice(i, j) > p%dyn%dyn_area_min .and. tmass(i, j) > p%dyn%dyn_mass_min
            end if
         end do
      end do
      iceT = .false.
      do j = js - 1, je + 1
         do i = is - 1, ie + 1
            iceT(i, j) = fx%hm(i, j) > 0.5_wp .and. any(tph(i - 1:i + 1, j - 1:j + 1))
            if (.not. iceT(i, j)) then
               spo(i, j) = 0.0_wp
               smo(i, j) = 0.0_wp
            end if
         end do
      end do
      ! ---- E, N, U points ----
      iceE = .false.; iceN = .false.; iceU = .false.
      aE = 0.0_wp; aN = 0.0_wp; mE = 0.0_wp; mN = 0.0_wp; fmE = 0.0_wp; fmN = 0.0_wp; fxE = 0.0_wp; fyN = 0.0_wp
      ui = 0.0_wp; vi = 0.0_wp
      do j = js, je
         do i = is, ie
            w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i + 1, j)*fx%tarea(i + 1, j); wt = w1 + w2
            m = 0.0_wp; a = 0.0_wp
            if (wt /= 0.0_wp) then
               m = (w1*tmass(i, j) + w2*tmass(i + 1, j))/wt
               a = (w1*fx%aice_init(i, j) + w2*fx%aice_init(i + 1, j))/wt
            end if
            iceE(i, j) = fx%epm(i, j) > 0.5_wp .and. a > p%dyn%dyn_area_min .and. m > p%dyn%dyn_mass_min
            if (iceE(i, j)) then
               if (.not. iceE0(i, j)) uE(i, j) = fx%uocnE(i, j)
               aE(i, j) = a; mE(i, j) = m/dtl; fmE(i, j) = fx%fcorE(i, j)*m
               fxE(i, j) = fx%strairxE(i, j) - g*m*fx%ss_tltxE(i, j)
               ui(i, j) = uE(i, j)
            else
               uE(i, j) = 0.0_wp
            end if
            w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i, j + 1)*fx%tarea(i, j + 1); wt = w1 + w2
            m = 0.0_wp; a = 0.0_wp
            if (wt /= 0.0_wp) then
               m = (w1*tmass(i, j) + w2*tmass(i, j + 1))/wt
               a = (w1*fx%aice_init(i, j) + w2*fx%aice_init(i, j + 1))/wt
            end if
            iceN(i, j) = fx%npm(i, j) > 0.5_wp .and. a > p%dyn%dyn_area_min .and. m > p%dyn%dyn_mass_min
            if (iceN(i, j)) then
               if (.not. iceN0(i, j)) vN(i, j) = fx%vocnN(i, j)
               aN(i, j) = a; mN(i, j) = m/dtl; fmN(i, j) = fx%fcorN(i, j)*m
               fyN(i, j) = fx%strairyN(i, j) - g*m*fx%ss_tltyN(i, j)
               vi(i, j) = vN(i, j)
            else
               vN(i, j) = 0.0_wp
            end if
         end do
      end do
      s12o = s120
      do j = js - 2, je + 1
         do i = is - 2, ie + 1
            w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i + 1, j)*fx%tarea(i + 1, j)
            w3 = fx%hm(i, j + 1)*fx%tarea(i, j + 1); w4 = fx%hm(i + 1, j + 1)*fx%tarea(i + 1, j + 1)
            wt = w1 + w2 + w3 + w4
            m = 0.0_wp; a = 0.0_wp
            if (wt /= 0.0_wp) then
               m = (w1*tmass(i, j) + w2*tmass(i + 1, j) + w3*tmass(i, j + 1) + w4*tmass(i + 1, j + 1))/wt
               a = (w1*fx%aice_init(i, j) + w2*fx%aice_init(i + 1, j) + w3*fx%aice_init(i, j + 1) &
                    + w4*fx%aice_init(i + 1, j + 1))/wt
            end if
            iceU(i, j) = (fx%hm(i, j) + fx%hm(i + 1, j) + fx%hm(i, j + 1) + fx%hm(i + 1, j + 1)) > 1.5_wp &
                         .and. a > p%dyn%dyn_area_min .and. m > p%dyn%dyn_mass_min
            if (.not. iceU(i, j)) s12o(i + 1, j + 1) = 0.0_wp
         end do
      end do
      ! strength at U (avg_strength)
      sU = 0.0_wp
      do j = js - 1, je
         do i = is - 1, ie
            w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i + 1, j)*fx%tarea(i + 1, j)
            w3 = fx%hm(i, j + 1)*fx%tarea(i, j + 1); w4 = fx%hm(i + 1, j + 1)*fx%tarea(i + 1, j + 1)
            wt = w1 + w2 + w3 + w4
            if (wt /= 0.0_wp) sU(i, j) = (w1*fx%w%strength(i, j) + w2*fx%w%strength(i + 1, j) &
                                          + w3*fx%w%strength(i, j + 1) + w4*fx%w%strength(i + 1, j + 1))/wt
         end do
      end do
      call ref_avg()
      do ks = 1, fx%k%ndte
         ! K1
         do j = js - 2, je + 1
            do i = is - 2, ie + 1
               shU(i, j) = 0.0_wp
               dlU(i, j) = 0.0_wp
               if (.not. iceU(i, j)) cycle
               call ref_corner(i, j, shU(i, j), dlU(i, j))
            end do
         end do
         ! K2
         do j = js - 1, je + 1
            do i = is - 1, ie + 1
               if (.not. iceT(i, j)) cycle
               divT = dy*(uE(i, j) - uE(i - 1, j)) + dx*(vN(i, j) - vN(i, j - 1))
               tenT = dy*dy*((uE(i, j) - uE(i - 1, j))/dy) - dx*dx*((vN(i, j) - vN(i, j - 1))/dx)
               ssq = 0.25_wp*(shU(i, j)**2 + shU(i, j - 1)**2 + shU(i - 1, j - 1)**2 + shU(i - 1, j)**2)
               dlt = sqrt(divT**2 + ef*(tenT**2 + ssq))
               call ref_visc(fx%w%strength(i, j), dmin, dlt, zeta, eta, rep)
               zo(i, j) = zeta; eo(i, j) = eta
               spo(i, j) = (spo(i, j)*keep + k_arlx1i*(zeta*divT - rep))*denom1
               smo(i, j) = (smo(i, j)*keep + k_arlx1i*eta*tenT)*denom1
            end do
         end do
         ! K3
         do j = js - 1, je
            do i = is - 1, ie
               if (.not. iceU(i, j)) cycle
               if (avgs) then
                  call ref_visc(sU(i, j), fx%k%deltamin*ar, dlU(i, j), zeta, eta, rep)
               else
                  w1 = fx%hm(i, j)*fx%tarea(i, j); w2 = fx%hm(i + 1, j)*fx%tarea(i + 1, j)
                  w3 = fx%hm(i, j + 1)*fx%tarea(i, j + 1); w4 = fx%hm(i + 1, j + 1)*fx%tarea(i + 1, j + 1)
                  eta = (w1*eo(i, j) + w2*eo(i + 1, j) + w3*eo(i, j + 1) + w4*eo(i + 1, j + 1))/(w1 + w2 + w3 + w4)
               end if
               s12o(i + 1, j + 1) = (s12o(i + 1, j + 1)*keep + k_arlx1i*0.5_wp*eta*shU(i, j))*denom1
            end do
         end do
         ! K4
         do j = js, je
            do i = is, ie
               if (iceE(i, j)) then
                  rheo = merge(0.0_wp, 1.0_wp, fd .and. .not. (aE(i, j) > p%dyn%rheo_area_min))
                  strint = rheo/ar*(0.5_wp*dy*(spo(i + 1, j) - spo(i, j)) + (0.5_wp/dy)*dy*dy*(smo(i + 1, j) - smo(i, j)) &
                                    + (1.0_wp/dx)*dx*dx*(s12o(i + 1, j + 1) - s12o(i + 1, j)))
                  uold = uE(i, j); vold = vvE(i, j)
                  vrel = aE(i, j)*rho_w*cw*sqrt((fx%uocnE(i, j) - uold)**2 + (fx%vocnE(i, j) - vold)**2)
                  sg = sign(1.0_wp, fmE(i, j))
                  cca = (brlx + revp)*mE(i, j) + vrel*cosw
                  ccb = fmE(i, j) + sg*vrel*sinw
                  cc = strint + fxE(i, j) + vrel*(fx%uocnE(i, j)*cosw - fx%vocnE(i, j)*sinw*sg) &
                       + mE(i, j)*(brlx*uold + revp*ui(i, j))
                  uE(i, j) = (ccb*vold + cc)/cca
               end if
               if (iceN(i, j)) then
                  rheo = merge(0.0_wp, 1.0_wp, fd .and. .not. (aN(i, j) > p%dyn%rheo_area_min))
                  strint = rheo/ar*(0.5_wp*dx*(spo(i, j + 1) - spo(i, j)) - (0.5_wp/dx)*dx*dx*(smo(i, j + 1) - smo(i, j)) &
                                    + (1.0_wp/dy)*dy*dy*(s12o(i + 1, j + 1) - s12o(i, j + 1)))
                  uold = uvN(i, j); vold = vN(i, j)
                  vrel = aN(i, j)*rho_w*cw*sqrt((fx%uocnN(i, j) - uold)**2 + (fx%vocnN(i, j) - vold)**2)
                  sg = sign(1.0_wp, fmN(i, j))
                  cca = (brlx + revp)*mN(i, j) + vrel*cosw
                  ccb = fmN(i, j) + sg*vrel*sinw
                  cc = strint + fyN(i, j) + vrel*(fx%vocnN(i, j)*cosw + fx%uocnN(i, j)*sinw*sg) &
                       + mN(i, j)*(brlx*vold + revp*vi(i, j))
                  vN(i, j) = (-ccb*uold + cc)/cca
               end if
            end do
         end do
         call ref_avg()
      end do
      ! finish
      divu = 0.0_wp; shear = 0.0_wp; vort = 0.0_wp; rc = 0.0_wp; rs = 0.0_wp; sox = 0.0_wp; soy = 0.0_wp
      do j = js, je
         do i = is, ie
            if (iceT(i, j)) then
               divT = dy*(uE(i, j) - uE(i - 1, j)) + dx*(vN(i, j) - vN(i, j - 1))
               tenT = dy*dy*((uE(i, j) - uE(i - 1, j))/dy) - dx*dx*((vN(i, j) - vN(i, j - 1))/dx)
               shrT = dx*dx*((uvN(i, j) - uvN(i, j - 1))/dx) + dy*dy*((vvE(i, j) - vvE(i - 1, j))/dy)
               ssq = 0.25_wp*(shU(i, j)**2 + shU(i, j - 1)**2 + shU(i - 1, j - 1)**2 + shU(i - 1, j)**2)
               dlt = sqrt(divT**2 + ef*(tenT**2 + ssq))
               dpp2 = divT**2 + ep*(tenT**2 + ssq)
               dv = divT/ar
               divu(i, j) = dv
               tmp = dpp2/max(dlt, 1.0e-13_wp)/ar
               rc(i, j) = -min(dv, 0.0_wp)
               rs(i, j) = 0.5_wp*(tmp - abs(dv))
               shear(i, j) = sqrt(tenT**2 + shrT**2)/ar
               vort(i, j) = ((dy*vvE(i, j) - dy*vvE(i - 1, j)) - (dx*uvN(i, j) - dx*uvN(i, j - 1)))/ar
            end if
            if (iceE(i, j)) then
               vrel = aE(i, j)*rho_w*cw*sqrt((fx%uocnE(i, j) - uE(i, j))**2 + (fx%vocnE(i, j) - vvE(i, j))**2)
               sg = sign(1.0_wp, fmE(i, j))
               sox(i, j) = vrel*((fx%uocnE(i, j) - uE(i, j))*cosw - (fx%vocnE(i, j) - vvE(i, j))*sinw*sg)
            end if
            if (iceN(i, j)) then
               vrel = aN(i, j)*rho_w*cw*sqrt((fx%uocnN(i, j) - uvN(i, j))**2 + (fx%vocnN(i, j) - vN(i, j))**2)
               sg = sign(1.0_wp, fmN(i, j))
               soy(i, j) = vrel*((fx%vocnN(i, j) - vN(i, j))*cosw + (fx%uocnN(i, j) - uvN(i, j))*sinw*sg)
            end if
         end do
      end do
      uout = u0
      vout = v0
      do j = 1, ny
         do i = 1, nx
            uout(i + 1, j) = merge(uE(i, j), u0(i + 1, j), i >= is .and. i <= ie .and. j >= js .and. j <= je)
            vout(i, j + 1) = merge(vN(i, j), v0(i, j + 1), i >= is .and. i <= ie .and. j >= js .and. j <= je)
         end do
      end do
   contains

      subroutine ref_avg()
         integer :: ia, ja
         real(wp) :: q1, q2, q3, q4, qt, s
         do ja = js - 2, je + 2
            do ia = is - 2, ie + 2
               q1 = fx%earea(ia - 1, ja); q2 = fx%earea(ia, ja); q3 = fx%earea(ia - 1, ja + 1); q4 = fx%earea(ia, ja + 1)
               qt = q1 + q2 + q3 + q4
               s = (uE(ia - 1, ja)*q1 + uE(ia, ja)*q2 + uE(ia - 1, ja + 1)*q3 + uE(ia, ja + 1)*q4)/qt
               uvN(ia, ja) = s*fx%npm(ia, ja)
               q1 = fx%narea(ia, ja - 1); q2 = fx%narea(ia + 1, ja - 1); q3 = fx%narea(ia, ja); q4 = fx%narea(ia + 1, ja)
               qt = q1 + q2 + q3 + q4
               s = (vN(ia, ja - 1)*q1 + vN(ia + 1, ja - 1)*q2 + vN(ia, ja)*q3 + vN(ia + 1, ja)*q4)/qt
               vvE(ia, ja) = s*fx%epm(ia, ja)
               q1 = fx%earea(ia, ja); q2 = fx%earea(ia, ja + 1); qt = q1 + q2
               uvU(ia, ja) = (uE(ia, ja)*q1 + uE(ia, ja + 1)*q2)/qt*fx%uvm(ia, ja)
               q1 = fx%narea(ia, ja); q2 = fx%narea(ia + 1, ja); qt = q1 + q2
               vvU(ia, ja) = (vN(ia, ja)*q1 + vN(ia + 1, ja)*q2)/qt*fx%uvm(ia, ja)
            end do
         end do
      end subroutine ref_avg

      subroutine ref_corner(ic, jc, sh, dl)
         !! Corner strain rate with the free-slip ratio extrapolation, written
         !! with explicit east/north neighbour values.
         integer, intent(in) :: ic, jc
         real(wp), intent(out) :: sh, dl
         real(wp) :: e1, e0, n1, n0, ue1, ue0, vn1, vn0, dv2, tn
         ! tangential-for-shear: u at E(ic,jc+1), E(ic,jc); v at N(ic+1,jc), N(ic,jc)
         ue1 = uE(ic, jc + 1)*fx%epm(ic, jc + 1) + (fx%epm(ic, jc) - fx%epm(ic, jc + 1))*fx%epm(ic, jc)*fx%rdyE(ic, jc)*uE(ic, jc)
         ue0 = uE(ic, jc)*fx%epm(ic, jc) + (fx%epm(ic, jc + 1) - fx%epm(ic, jc))*fx%epm(ic, jc + 1)*fx%rdyEr(ic, jc)*uE(ic, jc + 1)
         vn1 = vN(ic + 1, jc)*fx%npm(ic + 1, jc) + (fx%npm(ic, jc) - fx%npm(ic + 1, jc))*fx%npm(ic, jc)*fx%rdxN(ic, jc)*vN(ic, jc)
         vn0 = vN(ic, jc)*fx%npm(ic, jc) + (fx%npm(ic + 1, jc) - fx%npm(ic, jc))*fx%npm(ic + 1, jc)*fx%rdxNr(ic, jc)*vN(ic + 1, jc)
         sh = dx*(ue1 - ue0) + dy*(vn1 - vn0)
         e1 = uvN(ic + 1, jc)*fx%npm(ic + 1, jc) + (fx%npm(ic, jc) - fx%npm(ic + 1, jc))*fx%npm(ic, jc)*fx%rdxN(ic, jc)*uvN(ic, jc)
         e0 = uvN(ic, jc)*fx%npm(ic, jc) + (fx%npm(ic + 1, jc) - fx%npm(ic, jc))*fx%npm(ic + 1, jc)*fx%rdxNr(ic, jc)*uvN(ic + 1, jc)
         n1 = vvE(ic, jc + 1)*fx%epm(ic, jc + 1) + (fx%epm(ic, jc) - fx%epm(ic, jc + 1))*fx%epm(ic, jc)*fx%rdyE(ic, jc)*vvE(ic, jc)
         n0 = vvE(ic, jc)*fx%epm(ic, jc) + (fx%epm(ic, jc + 1) - fx%epm(ic, jc))*fx%epm(ic, jc + 1)*fx%rdyEr(ic, jc)*vvE(ic, jc + 1)
         dv2 = dy*(e1 - e0) + dx*(n1 - n0)
         tn = dy*(e1 - e0) - dx*(n1 - n0)
         dl = sqrt(dv2*dv2 + ef*(tn*tn + sh*sh))
      end subroutine ref_corner

      subroutine ref_visc(pp, dm, dl, zeta, eta, rep)
         real(wp), intent(in) :: pp, dm, dl
         real(wp), intent(out) :: zeta, eta, rep
         real(wp) :: t
         if (capmax) then
            t = pp/max(dl, dm)
         else
            t = pp/(dl + dm)
         end if
         zeta = (1.0_wp + kt)*t
         rep = (1.0_wp - kt)*t*dl
         eta = ep*zeta
      end subroutine ref_visc

   end subroutine ref_run

   subroutine reference_case(error, tag, revised, avg_strength, cap_sum)
      !! The public step against the scalar transcription on a random masked
      !! state (wall box + island, free-drift patches, Coriolis of both signs,
      !! a turning angle): every output on every cell, and the ring outside the
      !! computed windows untouched (bitwise).
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: tag
      logical, intent(in) :: revised, avg_strength, cap_sum
      type(fx_t) :: fx
      integer :: nx, ny, is, ie, js, je, i, j, nbad, nstep
      real(wp), allocatable :: u0(:, :), v0(:, :), sp0(:, :), sm0(:, :), s120(:, :), u1(:, :), v1(:, :), &
                               spo(:, :), smo(:, :), zo(:, :), eo(:, :), s12o(:, :), divu(:, :), shear(:, :), &
                               vort(:, :), rc(:, :), rs(:, :), sox(:, :), soy(:, :)
      logical, allocatable :: iceE0(:, :), iceN0(:, :), mT(:, :), mE(:, :), mN(:, :), mU(:, :)
      real(wp) :: tol

      tol = TOL_REF
      call fx_build(fx, 14, 12, .true., 2.0e4_wp, 1.5e4_wp, 8)
      fx%p%dyn%revised_evp = revised
      fx%p%dyn%arlx = 40.0_wp
      fx%p%dyn%brlx = 20.0_wp
      fx%p%dyn%cosw = cos(20.0_wp*PI/180.0_wp)
      fx%p%dyn%sinw = sin(20.0_wp*PI/180.0_wp)
      if (avg_strength) fx%p%dyn%visc_method = ICE_VISC_AVG_STRENGTH
      if (cap_sum) then
         fx%p%dyn%capping_max = .false.
         fx%p%dyn%Ktens = 0.1_wp
      end if
      fx%p%dyn%e_plasticpot = 2.3_wp
      call fx_consts(fx)
      call random_state(fx, 99173_int64, .true.)
      nx = fx%nx; ny = fx%ny; is = fx%is; ie = fx%ie; js = fx%js; je = fx%je
      allocate (u0(nx + 1, ny), v0(nx, ny + 1), sp0(nx, ny), sm0(nx, ny), s120(nx + 1, ny + 1), u1(nx + 1, ny), v1(nx, ny + 1))
      allocate (spo(nx, ny), smo(nx, ny), zo(nx, ny), eo(nx, ny), s12o(nx + 1, ny + 1), divu(nx, ny), shear(nx, ny))
      allocate (vort(nx, ny), rc(nx, ny), rs(nx, ny), sox(nx, ny), soy(nx, ny))
      allocate (iceE0(nx, ny), iceN0(nx, ny), mT(nx, ny), mE(nx, ny), mN(nx, ny), mU(nx, ny))
      u0 = fx%u; v0 = fx%v; sp0 = fx%w%stresspT; sm0 = fx%w%stressmT; s120 = fx%w%stress12U
      iceE0 = fx%w%iceEmask; iceN0 = fx%w%iceNmask
      call fx_map(fx)
      do nstep = 1, 2                    ! two steps: the persistent state is carried
         if (nstep == 2) then
            u0 = fx%u; v0 = fx%v; sp0 = fx%w%stresspT; sm0 = fx%w%stressmT; s120 = fx%w%stress12U
            iceE0 = fx%w%iceEmask; iceN0 = fx%w%iceNmask
         end if
         call fx_step(fx, .false., .false.)
         call fx_pull(fx)
         call ref_run(fx, u0, v0, iceE0, iceN0, sp0, sm0, s120, u1, v1, spo, smo, zo, eo, s12o, divu, shear, vort, rc, rs, &
                      sox, soy, mT, mE, mN, mU)
         ! outputs on the computed windows
         call expect_field(error, tag//" u", fx%u, u1, 1, nx + 1, 1, ny, tol)
         call expect_field(error, tag//" v", fx%v, v1, 1, nx, 1, ny + 1, tol)
         call expect_field(error, tag//" stresspT", fx%w%stresspT, spo, is - 1, ie + 1, js - 1, je + 1, tol)
         call expect_field(error, tag//" stressmT", fx%w%stressmT, smo, is - 1, ie + 1, js - 1, je + 1, tol)
         call expect_field(error, tag//" zetax2T", fx%w%zetax2T, zo, is - 1, ie + 1, js - 1, je + 1, tol)
         call expect_field(error, tag//" etax2T", fx%w%etax2T, eo, is - 1, ie + 1, js - 1, je + 1, tol)
         call expect_field(error, tag//" stress12U", fx%w%stress12U, s12o, is - 1, ie + 1, js - 1, je + 1, tol)
         call expect_field(error, tag//" divu", fx%w%divu, divu, is, ie, js, je, tol)
         call expect_field(error, tag//" shear", fx%w%shear, shear, is, ie, js, je, tol)
         call expect_field(error, tag//" vort", fx%w%vort, vort, is, ie, js, je, tol)
         call expect_field(error, tag//" rdg_conv", fx%w%rdg_conv, rc, is, ie, js, je, tol)
         call expect_field(error, tag//" rdg_shear", fx%w%rdg_shear, rs, is, ie, js, je, tol)
         call expect_field(error, tag//" strocnxE", fx%w%strocnxE, sox, is, ie, js, je, tol)
         call expect_field(error, tag//" strocnyN", fx%w%strocnyN, soy, is, ie, js, je, tol)
         if (allocated(error)) exit
         ! masks (exact)
         nbad = 0
         do j = js - 1, je + 1
            do i = is - 1, ie + 1
               if (fx%w%iceTmask(i, j) .neqv. mT(i, j)) nbad = nbad + 1
               if (i >= is .and. i <= ie .and. j >= js .and. j <= je) then
                  if (fx%w%iceEmask(i, j) .neqv. mE(i, j)) nbad = nbad + 1
                  if (fx%w%iceNmask(i, j) .neqv. mN(i, j)) nbad = nbad + 1
               end if
            end do
         end do
         do j = js - 2, je + 1
            do i = is - 2, ie + 1
               if (fx%w%iceUmask(i, j) .neqv. mU(i, j)) nbad = nbad + 1
            end do
         end do
         call check(error, nbad, 0, tag//": ice masks")
         if (allocated(error)) exit
         ! outside the written windows the input is untouched, bitwise
         nbad = 0
         do j = 1, ny
            do i = 1, nx + 1
               if (i <= is .or. i >= ie + 2 .or. j < js .or. j > je) then
                  if (fx%u(i, j) /= u0(i, j)) nbad = nbad + 1
               end if
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               if (j <= js .or. j >= je + 2 .or. i < is .or. i > ie) then
                  if (fx%v(i, j) /= v0(i, j)) nbad = nbad + 1
               end if
            end do
         end do
         do j = 1, ny
            do i = 1, nx
               if (i < is - 1 .or. i > ie + 1 .or. j < js - 1 .or. j > je + 1) then
                  if (fx%w%stresspT(i, j) /= sp0(i, j) .or. fx%w%stressmT(i, j) /= sm0(i, j)) nbad = nbad + 1
               end if
            end do
         end do
         call check(error, nbad, 0, tag//": cells outside the computed windows are untouched")
         if (allocated(error)) exit
      end do
      call fx_destroy(fx)
   end subroutine reference_case

   subroutine test_reference_classic(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, "ref classic", .false., .false., .false.)
   end subroutine test_reference_classic

   subroutine test_reference_revised(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, "ref revised", .true., .false., .false.)
   end subroutine test_reference_revised

   subroutine test_reference_avg_strength(error)
      type(error_type), allocatable, intent(out) :: error
      call reference_case(error, "ref avg_strength+sum", .false., .true., .true.)
   end subroutine test_reference_avg_strength

   subroutine test_throughput(error)
      !! Opt-in (`RDB_CICE_EVP_BENCH=1`): wall time of one step at the OM3 grid
      !! (1440 x 1152) with the CICE default ndte = 120.
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      character(len=16) :: env
      integer :: stat, nxp, nyp, c0, c1, cr, irep, nrep
      real(wp) :: tstep

      call get_environment_variable("RDB_CICE_EVP_BENCH", env, status=stat)
      if (stat /= 0) return
      if (trim(env) /= "1") return
      nxp = 1440
      nyp = 1152
      call fx_build(fx, nxp, nyp, .true., 2.5e4_wp, 2.5e4_wp, 120)
      call random_state(fx, 4242_int64, .true.)
      fx%w%iceEmask = .true.; fx%w%iceNmask = .true.
      call fx_map(fx)
      call fx_step(fx, .false., .false.)                 ! warm-up
      call fx_pull(fx)
      nrep = 3
      call system_clock(c0, cr)
      do irep = 1, nrep
         call fx_step(fx, .false., .false.)
      end do
      call fx_pull(fx)                                   ! synchronise
      call system_clock(c1)
      tstep = real(c1 - c0, wp)/real(cr, wp)/real(nrep, wp)
      write (error_unit, '(a,i0,a,i0,a,f9.4,a,f9.3,a)') "  [evp_c] ", nxp, " x ", nyp, " ndte=120: ", tstep, &
         " s per step, ", tstep/120.0_wp*1.0e3_wp, " ms per subcycle"
      call check(error, tstep > 0.0_wp, .true., "timed")
      call fx_destroy(fx)
   end subroutine test_throughput

end module test_cice_evp_c
