!! `&ocean_coriolis_nml form = "sadourny_hk"` next to a `z_fixed` STAIRCASE
!! (`&vcoord_nml zfixed_closed_faces`): thin live partial bottom cells,
!! inert fillers and closed face-layers.
!!
!! ### The defect this file guards
!!
!! Every Arakawa-Hsu term is `coef·transport`, `coef` the sum of the PVs
!! `q = (f+ζ)/h_corner` at three corners of the cell the two faces share.
!! For two of those three corners the far cell of one of the two
!! transports lies OUTSIDE the corner's 4-cell thickness mean (the
!! "cross" pairs), so the factor `h_face/h_corner` the PV·transport
!! product carries over `(f+ζ)·v` is UNBOUNDED.  A `z_fixed` staircase
!! realises it: the live partial bottom cell may be as thin as
!! `H_VANISHED` against a full-depth neighbour across an OPEN face.  On
!! the 1-degree Southern Ocean (50 tanh z-levels, closed faces, pred_corr)
!! the Coriolis tendency next to such cells grew 0.5 → 5.7 → 977 m/s per
!! stage in ten steps and the run was NaN at step 11, while `sadourny`
!! (velocity form, never divides by h) and `sadourny_energy` (each corner
!! only meets transports whose two cells are inside it ⇒
!! `h_corner >= h_face/2`) ran clean.  The fix gives the HK pairs the
!! energy form's bound: a PV is evaluated at a corner thickness of at
!! least half the larger face thickness of its pair.
!!
!! ### The cases
!!
!!   * `rest_is_exactly_zero` — a resting staircase (u = v = 0) gives an
!!     exactly-zero HK tendency on every face, with no non-finite value
!!     (the fillers' `(f+ζ)/h_corner` is large, but it must never leak).
!!   * `geostrophic_step_is_bounded` — uniform `u = U0`, `v = V0` on every
!!     OPEN face of the staircase (zero on the closed ones, as
!!     `mask_layer_velocities` leaves them).  The PV part of the tendency
!!     (`pv_flux + ∇KE`) must satisfy the energy form's bound
!!     `|CA| <= 2·max|f+ζ|·max(|U0|,|V0|)` on every face, and reproduce
!!     `f·V0` / `−f·U0` far from the step.  FAILS on the unfixed kernel:
!!     the bound is exceeded by the thickness ratio across the step.
!!   * `floor_inert_without_contrast` — on a smoothly varying, all-open
!!     column field the closed-face branch returns the original HK
!!     tendencies BIT FOR BIT (the floor is a no-op where no cell
!!     outweighs the other three of its corner).
!!   * `open_step_latch_is_bounded` — the same staircase with closed faces
!!     OFF (`zstar` / `z_fixed` / `zstar_full` over a stepped bed: every face
!!     open, fillers next to live cells) and flow on every face.  With the
!!     `coriolis_adv_t%hk_pair_floor` latch the energy-form bound holds;
!!     without it the same state exceeds the bound (non-vacuous: the
!!     unguarded kernel is what took the compat matrix's
!!     `zstar x sadourny_hk` cell to a negative layer at step 2).
!!   * `pair_antisymmetry_kept` — the Coriolis work
!!     `Σ_U uh·CAu/IdxCu + Σ_V vh·CAv/IdyCv` vanishes to round-off on the
!!     staircase with arbitrary velocities: the floor is a property of the
!!     PAIR, so the coefficient the u-tendency applies to `vh_V` is the one
!!     the v-tendency applies to `uh_U` (HK energy conservation).
module test_ocean_coriolis_hk_vanished
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_coriolis_adv, only: coriolis_adv_t, coriolis_adv_compute_tendencies_hk, &
                               PV_VARIANT_SADOURNY_HK
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_vcoord, only: ocean_vcoord_closed_face_masks
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_coriolis_hk_vanished_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 12, NYP = 12, NZ = 2
   real(wp), parameter :: DX = 1.0e3_wp
   real(wp), parameter :: F0 = 1.0e-4_wp
   real(wp), parameter :: FILLER = 1.0e-4_wp
      !! `zstar_h_min` at its default — an inert filler, below H_VANISHED.
   real(wp), parameter :: THIN = 1.0e-2_wp
      !! A LIVE partial bottom cell (above H_VANISHED), the staircase step.
   real(wp), parameter :: THICK = 50.0_wp
      !! A full z-level.

contains

   subroutine collect_ocean_coriolis_hk_vanished_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("rest_is_exactly_zero", test_rest), &
                  new_unittest("geostrophic_step_is_bounded", test_geostrophic_step), &
                  new_unittest("floor_inert_without_contrast", test_floor_inert), &
                  new_unittest("pair_antisymmetry_kept", test_antisymmetry), &
                  new_unittest("open_step_latch_is_bounded", test_open_step_latch) &
                  ]
   end subroutine collect_ocean_coriolis_hk_vanished_tests

   ! -----------------------------------------------------------------
   ! Fixtures
   ! -----------------------------------------------------------------

   subroutine build_staircase(grid, ms, metrics)
      !! Two layers (`k = 1` the bed).  Layer 2 is a full level everywhere.
      !! Layer 1 is a FILLER in the south-west patch, a THIN live partial
      !! cell along a west strip and a south strip, and a full level in the
      !! north-east: an x-step and a y-step of thickness ratio 5000, plus a
      !! filler patch whose faces the z-level mask closes.  The mask comes
      !! from the production builder, so "closed" has its one definition.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      integer :: i, j, ig, nx_t, ny_t

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ig = NGHOST
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      ! Masks sized BEFORE the device map (see make_cartesian_metrics).
      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
      do j = 1, ny_t
         do i = 1, nx_t
            ms%h_layer(i, j, 2) = THICK
            if (i <= ig + 3 .and. j <= ig + 3) then
               ms%h_layer(i, j, 1) = FILLER
            else if (i <= ig + 6 .or. j <= ig + 6) then
               ms%h_layer(i, j, 1) = THIN
            else
               ms%h_layer(i, j, 1) = THICK
            end if
         end do
      end do
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, ms%h_layer, &
                                          nx_t, ny_t, NZ, H_VANISHED)
      !$omp target update to(metrics%open_u, metrics%open_v)
      metrics%use_closed_faces = .true.
   end subroutine build_staircase

   subroutine set_open_velocities(grid, ms, metrics, u0, v0, wobble)
      !! `u0`/`v0` on every OPEN interior face, zero on closed faces and on
      !! the array-edge faces (no transport leaves the patch).  `wobble`
      !! adds a deterministic face-to-face variation (arbitrary flow).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: u0, v0
      logical, intent(in) :: wobble
      integer :: i, j, k, nx_t, ny_t
      real(wp) :: w

      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_t
            do i = 2, nx_t
               w = 1.0_wp
               if (wobble) w = sin(1.3_wp*real(i, wp) + 0.7_wp*real(j, wp) + real(k, wp))
               ms%u_face_x_layer(i, j, k) = u0*w*metrics%open_u(i, j, k)
            end do
         end do
         do j = 2, ny_t
            do i = 1, nx_t
               w = 1.0_wp
               if (wobble) w = cos(0.9_wp*real(i, wp) - 1.1_wp*real(j, wp) + real(k, wp))
               ms%v_face_y_layer(i, j, k) = v0*w*metrics%open_v(i, j, k)
            end do
         end do
      end do
   end subroutine set_open_velocities

   subroutine run_hk(grid, metrics, ms, pvx, pvy, kec, mfu, mfv, pair_floor)
      !! One HK evaluation on the host state `ms`, every buffer mapped and
      !! read back (`mem:separate`-safe).  Returns the tendencies, the
      !! centre KE and the (masked) transports the kernel used.
      !! `pair_floor` sets the `hk_pair_floor` latch (default off).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), allocatable, intent(out) :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp), allocatable, intent(out) :: mfu(:, :, :), mfv(:, :, :)
      logical, intent(in), optional :: pair_floor
      type(coriolis_adv_t) :: cor

      cor%f_0 = F0
      cor%pv_variant = PV_VARIANT_SADOURNY_HK
      if (present(pair_floor)) cor%hk_pair_floor = pair_floor
      call cor%init(grid, nz_ml=NZ)
      !$omp target enter data map(to: ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_hk(grid, metrics, cor, ms, &
                                              ms%u_face_x_layer, ms%v_face_y_layer, &
                                              ms%h_layer)
      !$omp target update from(cor%pv_flux_x%data, cor%pv_flux_y%data, cor%ke_centre%data)
      !$omp target update from(cor%mass_flux_u%data, cor%mass_flux_v%data)
      pvx = cor%pv_flux_x%data
      pvy = cor%pv_flux_y%data
      kec = cor%ke_centre%data
      mfu = cor%mass_flux_u%data
      mfv = cor%mass_flux_v%data
      call cor%exit_data(); call ms%exit_data()
      !$omp target exit data map(delete: ms, cor)
      call cor%destroy()
   end subroutine run_hk

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_rest(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: mfu(:, :, :), mfv(:, :, :)

      call build_staircase(grid, ms, metrics)
      call set_open_velocities(grid, ms, metrics, 0.0_wp, 0.0_wp, .false.)
      call run_hk(grid, metrics, ms, pvx, pvy, kec, mfu, mfv)

      call check(error, all(ieee_is_finite(pvx)) .and. all(ieee_is_finite(pvy)), &
                 "resting staircase: non-finite HK tendency")
      if (allocated(error)) goto 99
      call check(error, maxval(abs(pvx)) == 0.0_wp .and. maxval(abs(pvy)) == 0.0_wp, &
                 "resting staircase: HK tendency must be EXACTLY zero")

99    call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_rest

   subroutine test_geostrophic_step(error)
      !! Measured on the unfixed kernel: max|CA_pv| / bound = 1.30e2 (the
      !! cross pair `q_S·vh_N` across the y-step, `h_face/h_corner ≈ 25/0.01`
      !! damped by the 1/12 weight); with the pair floor the ratio is <= 1.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: mfu(:, :, :), mfv(:, :, :)
      real(wp), parameter :: U0 = 0.02_wp, V0 = 0.01_wp
      real(wp) :: bound, worst, pv, abs_vort_max
      integer :: i, j, k, nx_t, ny_t, ifar, jfar
      character(len=200) :: msg

      call build_staircase(grid, ms, metrics)
      call set_open_velocities(grid, ms, metrics, U0, V0, .false.)
      call run_hk(grid, metrics, ms, pvx, pvy, kec, mfu, mfv)
      nx_t = grid%nx_total
      ny_t = grid%ny_total

      ! |ζ| at a corner is at most (2·|V0|·dy + 2·|U0|·dx)/(dx·dy): only the
      ! closed faces (zero velocity) next to open ones (U0, V0) shear.
      abs_vort_max = abs(F0) + 2.0_wp*(abs(U0) + abs(V0))/DX
      bound = 2.0_wp*abs_vort_max*max(abs(U0), abs(V0))*(1.0_wp + 1.0e-12_wp)

      worst = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_t
            do i = 2, nx_t
               pv = pvx(i, j, k) + (kec(i, j, k) - kec(i - 1, j, k))/DX
               worst = max(worst, abs(pv))
            end do
         end do
         do j = 2, ny_t
            do i = 1, nx_t
               pv = pvy(i, j, k) + (kec(i, j, k) - kec(i, j - 1, k))/DX
               worst = max(worst, abs(pv))
            end do
         end do
      end do
      write (msg, "(a,es12.4,a,es12.4,a,es12.4)") "max|CA_pv| = ", worst, &
         "  bound 2|f+zeta||v| = ", bound, "  ratio = ", worst/bound
      call check(error, all(ieee_is_finite(pvx)) .and. all(ieee_is_finite(pvy)), &
                 "staircase flow: non-finite HK tendency")
      if (allocated(error)) goto 99
      call check(error, worst <= bound, &
                 "HK PV tendency exceeds the energy-form bound next to the step: "//trim(msg))
      if (allocated(error)) goto 99

      ! Non-vacuous: deep inside the thick, uniform region the HK form is
      ! plain geostrophy, CAu = f·V0 and CAv = −f·U0.
      ifar = NGHOST + 10
      jfar = NGHOST + 10
      call check(error, abs(pvx(ifar, jfar, 1) - F0*V0) <= 1.0e-12_wp*F0*V0, &
                 "far from the step CAu must be f*V0")
      if (allocated(error)) goto 99
      call check(error, abs(pvy(ifar, jfar, 1) + F0*U0) <= 1.0e-12_wp*F0*U0, &
                 "far from the step CAv must be -f*U0")

99    call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_geostrophic_step

   subroutine test_floor_inert(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: pvx_f(:, :, :), pvy_f(:, :, :), pvx_o(:, :, :), pvy_o(:, :, :)
      real(wp), allocatable :: kec(:, :, :), mfu(:, :, :), mfv(:, :, :)
      integer :: i, j, nx_t, ny_t

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
      ! All open (the masks stay at 1): only the floor differs between runs.
      do j = 1, ny_t
         do i = 1, nx_t
            ms%h_layer(i, j, 1) = 40.0_wp + 1.5_wp*real(i, wp) + 0.5_wp*real(j, wp)
            ms%h_layer(i, j, 2) = 60.0_wp - 0.7_wp*real(i, wp) + 1.1_wp*real(j, wp)
         end do
      end do
      call set_open_velocities(grid, ms, metrics, 0.3_wp, -0.2_wp, .true.)

      metrics%use_closed_faces = .true.
      call run_hk(grid, metrics, ms, pvx_f, pvy_f, kec, mfu, mfv)
      metrics%use_closed_faces = .false.
      call run_hk(grid, metrics, ms, pvx_o, pvy_o, kec, mfu, mfv)

      call check(error, all(pvx_f == pvx_o) .and. all(pvy_f == pvy_o), &
                 "closed-face HK branch must be bit-identical to the original "// &
                 "where no thickness contrast engages the pair floor")

      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_floor_inert

   subroutine test_antisymmetry(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: mfu(:, :, :), mfv(:, :, :)
      real(wp) :: work, scale, t
      integer :: i, j, k, nx_t, ny_t
      character(len=160) :: msg

      call build_staircase(grid, ms, metrics)
      call set_open_velocities(grid, ms, metrics, 0.2_wp, 0.15_wp, .true.)
      call run_hk(grid, metrics, ms, pvx, pvy, kec, mfu, mfv)
      nx_t = grid%nx_total
      ny_t = grid%ny_total

      ! Coriolis work of the PV part (the KE gradient is not part of the
      ! antisymmetric operator).  Cartesian: IdxCu = IdyCv = 1/DX.
      work = 0.0_wp
      scale = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_t
            do i = 2, nx_t
               t = mfu(i, j, k)*(pvx(i, j, k) + (kec(i, j, k) - kec(i - 1, j, k))/DX)*DX
               work = work + t
               scale = scale + abs(t)
            end do
         end do
         do j = 2, ny_t
            do i = 1, nx_t
               t = mfv(i, j, k)*(pvy(i, j, k) + (kec(i, j, k) - kec(i, j - 1, k))/DX)*DX
               work = work + t
               scale = scale + abs(t)
            end do
         end do
      end do
      write (msg, "(a,es12.4,a,es12.4)") "sum = ", work, "  sum|terms| = ", scale
      call check(error, scale > 0.0_wp, "antisymmetry check is vacuous (no work terms)")
      if (allocated(error)) goto 99
      call check(error, abs(work) <= 1.0e-12_wp*scale, &
                 "HK Coriolis work must vanish (pair antisymmetry): "//trim(msg))

99    call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_antisymmetry

   subroutine test_open_step_latch(error)
      !! Closed faces OFF over the staircase: every face open, so a live
      !! cell faces a FILLER (1e-4 m) across the patch edge.  Uniform flow
      !! on every interior face.  The pair-floor latch must hold the
      !! energy-form bound; the unlatched kernel must not (the defect).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: mfu(:, :, :), mfv(:, :, :)
      real(wp), parameter :: U0 = 0.02_wp, V0 = 0.01_wp
      real(wp) :: bound, worst(2), abs_vort_max
      integer :: ilatch
      character(len=200) :: msg

      do ilatch = 1, 2
         call build_staircase(grid, ms, metrics)
         ! Open steps: no z-level mask anywhere.
         metrics%use_closed_faces = .false.
         metrics%open_u = 1.0_wp
         metrics%open_v = 1.0_wp
         !$omp target update to(metrics%open_u, metrics%open_v)
         call set_open_velocities(grid, ms, metrics, U0, V0, .false.)
         call run_hk(grid, metrics, ms, pvx, pvy, kec, mfu, mfv, pair_floor=(ilatch == 1))
         worst(ilatch) = max_pv_tendency(grid, pvx, pvy, kec)
         call check(error, all(ieee_is_finite(pvx)) .and. all(ieee_is_finite(pvy)), &
                    "open staircase: non-finite HK tendency")
         call ms%destroy()
         call destroy_cartesian_metrics(metrics)
         if (allocated(error)) return
      end do
      abs_vort_max = abs(F0) + 2.0_wp*(abs(U0) + abs(V0))/DX
      bound = 2.0_wp*abs_vort_max*max(abs(U0), abs(V0))*(1.0_wp + 1.0e-12_wp)
      write (msg, "(a,es12.4,a,es12.4,a,es12.4)") "latched ", worst(1), &
         "  unlatched ", worst(2), "  bound ", bound
      call check(error, worst(1) <= bound, &
                 "open staircase, hk_pair_floor: HK exceeds the energy-form bound: "//trim(msg))
      if (allocated(error)) return
      call check(error, worst(2) > bound, &
                 "open staircase WITHOUT the latch should exceed the bound (vacuous test?): " &
                 //trim(msg))
   end subroutine test_open_step_latch

   function max_pv_tendency(grid, pvx, pvy, kec) result(worst)
      !! max |pv_flux + grad KE| over every face: the PV part of the HK
      !! tendency, the quantity the energy-form bound limits.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: pvx(:, :, :), pvy(:, :, :), kec(:, :, :)
      real(wp) :: worst
      integer :: i, j, k
      worst = 0.0_wp
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 2, grid%nx_total
               worst = max(worst, abs(pvx(i, j, k) + (kec(i, j, k) - kec(i - 1, j, k))/DX))
            end do
         end do
         do j = 2, grid%ny_total
            do i = 1, grid%nx_total
               worst = max(worst, abs(pvy(i, j, k) + (kec(i, j, k) - kec(i, j - 1, k))/DX))
            end do
         end do
      end do
   end function max_pv_tendency

end module test_ocean_coriolis_hk_vanished
