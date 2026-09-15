!! Coverage for the `&ocean_vdiff_nml hvel_mom6` + `bbl_glue`
!! momentum-vdiff path (the MOM6 HARMONIC_VISC face-thickness build and
!! the rest-state BBL-glue viscous absorber, `diffuse_velocity_columns_impl`).
!!
!! `test_ocean_visc_rem` already exercises the plain `implicit_drag` +
!! `lambda_bot` remnant path thoroughly, but every one of its cases runs
!! with `hvel_mom6 = .false.` and `bbl_glue = .false.`, so the harmonic
!! `hvel` / `h_shear` build, the near-bed upwind blend, the `zint`/`botfn`
!! bookkeeping, and the piston bed-coupling all had ZERO coverage.
!!
!! Convention reminder (Roundabout is bottom-up): `k = 1` is the BED
!! (γ → smallest), `k = nz` is the SURFACE (γ → 1).
!!
!! Cases:
!!   1. hvel_mom6_visc_rem_bounded — with the harmonic hvel/h_shear build
!!      on (drag on, stratified physical column), the viscous remnant γ
!!      obeys the discrete maximum principle (finite, 0 < γ ≤ 1,
!!      non-decreasing bed→surface).  Exercises the whole hvel_mom6
!!      branch that the existing suite never enters.
!!   2. bbl_glue_visc_rem_finite_thin_bed — the full glue path
!!      (hvel_mom6 + bbl_glue + implicit_drag) on a grounded-stack column
!!      (a thin bed sliver thickening upward): γ must stay FINITE (no NaN
!!      from an unguarded `hf_k` / `dz` division — the failure mode the
!!      pre-merge review flagged) and bounded 0 < γ ≤ 1.
!!   3. bbl_glue_damps_near_bed — the absorber does its job: with a thin
!!      bed layer and a uniform velocity, turning bbl_glue ON (piston bed
!!      coupling `kv_bbl/(hf₁·min(hf₁/2, bbl_thick))`, which diverges as
!!      the bed layer thins) damps the near-bed velocity strictly harder
!!      than the plain-`lambda` bed drag with everything else identical.
!!   4. face_thick_blend_pointwise — a direct pointwise check of the
!!      public face-thickness blend for known `h_a`/`h_b`: arithmetic
!!      mean when `use_harmonic = .false.`, harmonic mean when `.true.`.
!!   5. hvel_upwind_gate — the near-bed upwind (arithmetic-donor) blend
!!      is gated EXACTLY on `u_face·(h_r − h_l) < 0`: with no thick→thin
!!      flow (uniform layers) `hvel_upwind` on vs off is bit-identical
!!      (blend never triggers); with a thick→thin bed column and the
!!      matching velocity sign it triggers and the two differ.
!!   6. collapsed_interior_layer_gamma_finite — REGRESSION for the
!!      momentum-vdiff `hf_k` floor: an EXACTLY-collapsed interior layer
!!      (both neighbour cells `h_layer = 0` ⇒ `hvel = 0` at that face)
!!      with `kv > 0` at that interface.  Without flooring the face
!!      thicknesses at `H_VANISHED`, the α/β denominators `hf_k·dz` go to
!!      0 ⇒ Inf diagonal ⇒ NaN viscous remnant γ.  With the floor γ stays
!!      finite and bounded 0 < γ ≤ 1 (the tracer path was already
!!      guarded; the momentum path was the gap the review flagged).
module test_ocean_vdiff_bbl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_momentum, face_thick
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_vdiff_bbl_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_vdiff_bbl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("hvel_mom6_visc_rem_bounded", test_hvel_mom6_bounded), &
                  new_unittest("bbl_glue_visc_rem_finite_thin_bed", test_bbl_glue_finite), &
                  new_unittest("bbl_glue_damps_near_bed", test_bbl_glue_damps), &
                  new_unittest("face_thick_blend_pointwise", test_face_thick_pointwise), &
                  new_unittest("hvel_upwind_gate", test_hvel_upwind_gate), &
                  new_unittest("collapsed_interior_layer_gamma_finite", &
                               test_collapsed_layer_finite) &
                  ]
   end subroutine collect_ocean_vdiff_bbl_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vd)
      call vd%enter_data()
   end subroutine map_in

   subroutine map_out(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      call vd%exit_data()
      !$acc exit data delete(vd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------

   subroutine test_hvel_mom6_bounded(error)
      !! hvel_mom6 harmonic face-thickness build + strong bottom drag,
      !! stratified physical column: γ finite, 0 < γ ≤ 1, non-decreasing
      !! bed→surface, strictly γ_bed < γ_surface (λ_bot > 0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 6
      real(wp), parameter :: DT = 450.0_wp, LAM = 0.08_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      logical :: bounded, monotone, no_nan, strict
      checks: block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.2_wp
         vd%implicit_drag = .true.
         vd%hvel_mom6 = .true.        ! <-- the previously-uncovered branch
         nu_face = size(ms%u_face_x_layer, 1)
         nv_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         nv_face = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         do k = 1, NZ
            ms%h_layer(:, :, k) = real(k, wp)*12.0_wp   ! stratified, thickening up
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         allocate (lam_u(nu_face, nv_uface), source=LAM)
         allocate (lam_v(nx_vface, nv_face), source=LAM)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         bounded = .true.; monotone = .true.; no_nan = .true.
         do k = 1, NZ
            if (.not. ieee_is_finite(bt%visc_rem_u(ig, jg, k))) no_nan = .false.
            if (bt%visc_rem_u(ig, jg, k) <= 0.0_wp .or. &
                bt%visc_rem_u(ig, jg, k) > 1.0_wp + 1.0e-12_wp) bounded = .false.
            if (k > 1) then
               if (bt%visc_rem_u(ig, jg, k) < bt%visc_rem_u(ig, jg, k - 1) - 1.0e-13_wp) &
                  monotone = .false.
            end if
         end do
         strict = bt%visc_rem_u(ig, jg, 1) < bt%visc_rem_u(ig, jg, NZ) - 1.0e-10_wp

         call check(error, no_nan, "hvel_mom6: visc_rem NaN/Inf")
         if (allocated(error)) exit checks
         call check(error, bounded, "hvel_mom6: violated 0 < gamma <= 1")
         if (allocated(error)) exit checks
         call check(error, monotone, "hvel_mom6: gamma not non-decreasing bed->surface")
         if (allocated(error)) exit checks
         call check(error, strict, "hvel_mom6: bed gamma should be strictly < surface (lambda>0)")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call bt%destroy(); call vd%destroy(); call ms%destroy()
   end subroutine test_hvel_mom6_bounded

   subroutine test_bbl_glue_finite(error)
      !! The full BBL-glue path on a grounded-stack column (thin bed
      !! sliver thickening upward): the piston bed coupling divides by
      !! `hf₁·min(hf₁/2, bbl_thick)`, so a thin bed layer is exactly the
      !! regime where an unguarded division would produce a NaN.  γ must
      !! stay finite and bounded 0 < γ ≤ 1.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 4
      real(wp), parameter :: DT = 600.0_wp, LAM = 0.02_wp
      real(wp) :: h_prof(NZ)
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      logical :: bounded, no_nan
      checks: block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.1_wp
         vd%implicit_drag = .true.
         vd%hvel_mom6 = .true.
         vd%bbl_glue = .true.         ! <-- the full absorber path
         nu_face = size(ms%u_face_x_layer, 1)
         nv_uface = size(ms%u_face_x_layer, 2)
         nx_vface = size(ms%v_face_y_layer, 1)
         nv_face = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ! Grounded stack: thin bed sliver, thickening upward.
         h_prof = [0.05_wp, 1.0_wp, 10.0_wp, 25.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_prof(k)
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         allocate (lam_u(nu_face, nv_uface), source=LAM)
         allocate (lam_v(nx_vface, nv_face), source=LAM)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         bounded = .true.; no_nan = .true.
         do k = 1, NZ
            if (.not. ieee_is_finite(bt%visc_rem_u(ig, jg, k)) .or. &
                .not. ieee_is_finite(bt%visc_rem_v(ig, jg, k))) no_nan = .false.
            if (bt%visc_rem_u(ig, jg, k) <= 0.0_wp .or. &
                bt%visc_rem_u(ig, jg, k) > 1.0_wp + 1.0e-12_wp) bounded = .false.
         end do

         call check(error, no_nan, &
                    "bbl_glue: visc_rem NaN/Inf on a thin-bed column (hf_k division unguarded?)")
         if (allocated(error)) exit checks
         call check(error, bounded, "bbl_glue: violated 0 < gamma <= 1 on the thin-bed column")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call bt%destroy(); call vd%destroy(); call ms%destroy()
   end subroutine test_bbl_glue_finite

   subroutine test_bbl_glue_damps(error)
      !! Two momentum solves from an IDENTICAL uniform-velocity IC with a
      !! thin bed layer, differing ONLY in bbl_glue.  With glue OFF the
      !! bed drag is `dt·lambda`; with glue ON it is the piston
      !! `dt·kv_bbl/(hf₁·min(hf₁/2, bbl_thick))`, which for a thin bed
      !! layer is far larger — so the near-bed velocity must be damped
      !! strictly harder with glue on.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      integer, parameter :: NZ = 4
      real(wp), parameter :: DT = 300.0_wp, LAM = 0.01_wp, U0 = 0.5_wp
      ! Near-rest interior viscosity (< kv_bbl = bbl_piston*hbbl_visc =
      ! 3e-3): exactly the regime the glue targets, so its bed piston is
      ! the dominant near-bed sink rather than the interior coupling.
      real(wp), parameter :: KV_SMALL = 1.0e-4_wp
      real(wp) :: h_prof(NZ)
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :)
      integer :: k, nu_face, nv_uface, nx_vface, nv_face, ig, jg
      real(wp) :: ubed_a, ubed_b
      checks: block
         call make_grid(grid, 4, 4)
         ms_a%nz_ml = NZ; ms_b%nz_ml = NZ
         call ms_a%init(grid); call ms_b%init(grid)
         call vd_a%init(grid, nz_ml=NZ); call vd_b%init(grid, nz_ml=NZ)
         ! A: hvel_mom6 + plain lambda bed drag (glue OFF).
         vd_a%K_v_momentum = KV_SMALL
         vd_a%implicit_drag = .true.
         vd_a%hvel_mom6 = .true.
         vd_a%bbl_glue = .false.
         ! B: identical, but glue ON.
         vd_b%K_v_momentum = KV_SMALL
         vd_b%implicit_drag = .true.
         vd_b%hvel_mom6 = .true.
         vd_b%bbl_glue = .true.
         nu_face = size(ms_a%u_face_x_layer, 1)
         nv_uface = size(ms_a%u_face_x_layer, 2)
         nx_vface = size(ms_a%v_face_y_layer, 1)
         nv_face = size(ms_a%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         h_prof = [0.2_wp, 5.0_wp, 15.0_wp, 25.0_wp]
         do k = 1, NZ
            ms_a%h_layer(:, :, k) = h_prof(k)
            ms_b%h_layer(:, :, k) = h_prof(k)
         end do
         ms_a%u_face_x_layer = U0; ms_b%u_face_x_layer = U0
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp
         allocate (lam_u(nu_face, nv_uface), source=LAM)
         allocate (lam_v(nx_vface, nv_face), source=LAM)

         call map_in(ms_a, vd_a)
         call vdiff_apply_momentum(grid, vd_a, ms_a, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         call map_out(ms_a, vd_a)

         call map_in(ms_b, vd_b)
         call vdiff_apply_momentum(grid, vd_b, ms_b, DT, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         call map_out(ms_b, vd_b)

         ! Bed layer is k = 1.  Both start at U0 > 0 and can only be
         ! damped toward 0, so smaller residual == stronger damping.
         ubed_a = ms_a%u_face_x_layer(ig, jg, 1)
         ubed_b = ms_b%u_face_x_layer(ig, jg, 1)

         call check(error, ubed_a > 0.0_wp .and. ubed_a < U0, &
                    "damps: glue-off bed velocity should be partially damped, still positive")
         if (allocated(error)) exit checks
         call check(error, ubed_b < ubed_a - 1.0e-6_wp, &
                    "damps: bbl_glue ON must damp the near-bed velocity harder than plain drag")
      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v)
      call vd_a%destroy(); call vd_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_bbl_glue_damps

   subroutine test_face_thick_pointwise(error)
      !! Direct pointwise check of the public face-thickness blend for
      !! known h_a, h_b: arithmetic mean vs harmonic mean.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HA = 2.0_wp, HB = 8.0_wp
      real(wp), parameter :: TOL = 1.0e-14_wp
      real(wp) :: dz_arith, dz_harm
      checks: block
         dz_arith = face_thick(HA, HB, .false.)
         dz_harm = face_thick(HA, HB, .true.)
         ! arithmetic: 0.5*(2+8) = 5
         call check(error, abs(dz_arith - 5.0_wp) < TOL, &
                    "face_thick arithmetic != 0.5*(h_a+h_b)")
         if (allocated(error)) exit checks
         ! harmonic: 2*2*8/(2+8) = 32/10 = 3.2
         call check(error, abs(dz_harm - 3.2_wp) < TOL, &
                    "face_thick harmonic != 2*h_a*h_b/(h_a+h_b)")
         if (allocated(error)) exit checks
         ! equal thicknesses: both means coincide.
         call check(error, abs(face_thick(4.0_wp, 4.0_wp, .true.) - 4.0_wp) < TOL, &
                    "face_thick harmonic of equal thicknesses != the common value")
      end block checks
   end subroutine test_face_thick_pointwise

   subroutine test_hvel_upwind_gate(error)
      !! The near-bed upwind (arithmetic-donor) blend inside the
      !! hvel_mom6 build is gated on `u_face·(h_r − h_l) < 0`.  (a) With
      !! uniform layers (h_r = h_l everywhere, h_delta = 0) the blend can
      !! never trigger, so hvel_upwind ON vs OFF is bit-identical.
      !! (b) With a thick→thin bed column across the face AND a velocity
      !! whose sign satisfies the trigger, ON vs OFF must differ (blend
      !! wired and gated exactly as documented).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: du_uniform, du_thickthin

      ! (a) Uniform layers -> h_delta = 0 -> no trigger -> bit-identical
      !     (the sheared u makes the solve itself non-trivial, so this is
      !     a real bit-identity, not a 0-vs-0 accident).
      call run_upwind_pair(uniform=.true., max_du=du_uniform)
      call check(error, du_uniform == 0.0_wp, &
                 "hvel_upwind: ON vs OFF must be bit-identical when no face is thick->thin")
      if (allocated(error)) return

      ! (b) Thick->thin bed face + positive bed u -> u*h_delta < 0 -> trigger -> differ.
      call run_upwind_pair(uniform=.false., max_du=du_thickthin)
      call check(error, du_thickthin > 1.0e-12_wp, &
                 "hvel_upwind: ON vs OFF must differ once a bed face runs thick->thin (blend wired)")
   end subroutine test_hvel_upwind_gate

   subroutine run_upwind_pair(uniform, max_du)
      !! Solve the momentum vdiff twice (hvel_upwind ON vs OFF, all else
      !! equal, hvel_mom6 on) and return the max |Δu| between them.  A
      !! SHEARED u profile (positive at the bed) drives a non-trivial
      !! viscous solve either way.  When `uniform`, every cell has the
      !! same thickness (h_delta = 0 on every face) so the near-bed
      !! upwind blend can never trigger.  Otherwise a west-thick /
      !! east-thin bed layer makes `h_delta = h_r − h_l < 0` at the
      !! interior x-faces, and the positive bed u there gives
      !! `u·h_delta < 0` — firing the blend at the bed layer only.
      logical, intent(in) :: uniform
      real(wp), intent(out) :: max_du

      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_on, ms_off
      type(ocean_vdiff_t) :: vd_on, vd_off
      integer, parameter :: NZ = 3
      real(wp), parameter :: DT = 300.0_wp
      real(wp), parameter :: U_PROF(NZ) = [0.4_wp, 0.1_wp, -0.2_wp]  !! sheared, bed +ve
      integer :: i, j, k, nx, ny
      real(wp) :: hbed

      call make_grid(grid, 6, 4)
      nx = grid%nx_total
      ny = grid%ny_total
      ms_on%nz_ml = NZ; ms_off%nz_ml = NZ
      call ms_on%init(grid); call ms_off%init(grid)
      call vd_on%init(grid, nz_ml=NZ); call vd_off%init(grid, nz_ml=NZ)

      vd_on%K_v_momentum = 0.3_wp
      vd_on%hvel_mom6 = .true.
      vd_on%hvel_upwind = .true.
      vd_off%K_v_momentum = 0.3_wp
      vd_off%hvel_mom6 = .true.
      vd_off%hvel_upwind = .false.

      do j = 1, ny
         do i = 1, nx
            ! Upper layers: fixed, uniform.
            ms_on%h_layer(i, j, 2) = 12.0_wp
            ms_on%h_layer(i, j, 3) = 20.0_wp
            ! Bed layer: uniform, or west-thick/east-thin ramp in i.
            if (uniform) then
               hbed = 8.0_wp
            else
               hbed = max(0.5_wp, 16.0_wp - 2.0_wp*real(i, wp))
            end if
            ms_on%h_layer(i, j, 1) = hbed
            do k = 1, NZ
               ms_off%h_layer(i, j, k) = ms_on%h_layer(i, j, k)
            end do
         end do
      end do
      do k = 1, NZ
         ms_on%u_face_x_layer(:, :, k) = U_PROF(k)
         ms_off%u_face_x_layer(:, :, k) = U_PROF(k)
      end do
      ms_on%v_face_y_layer = 0.0_wp
      ms_off%v_face_y_layer = 0.0_wp

      call map_in(ms_on, vd_on)
      call vdiff_apply_momentum(grid, vd_on, ms_on, DT)
      call map_out(ms_on, vd_on)

      call map_in(ms_off, vd_off)
      call vdiff_apply_momentum(grid, vd_off, ms_off, DT)
      call map_out(ms_off, vd_off)

      max_du = maxval(abs(ms_on%u_face_x_layer - ms_off%u_face_x_layer))

      call vd_on%destroy(); call vd_off%destroy()
      call ms_on%destroy(); call ms_off%destroy()
   end subroutine run_upwind_pair

   subroutine test_collapsed_layer_finite(error)
      !! REGRESSION for the momentum-vdiff hf_k floor.  Build a column with
      !! an EXACTLY-collapsed interior layer (h_layer = 0 at k = 2 in BOTH
      !! neighbour cells ⇒ hvel(2) = 0 at the face) and a positive uniform
      !! kv, then run the momentum vdiff solve requesting the viscous
      !! remnant γ.  Without the `hf_k = max(hvel, H_VANISHED)` floor the
      !! α/β denominators `hf_k·dz` collapse to 0 → Inf diagonal → NaN γ
      !! (the tracer path was already floored; the momentum path was not —
      !! the gap the pre-merge review flagged).  With the floor γ must stay
      !! finite and bounded 0 < γ ≤ 1 for EVERY layer, including the
      !! collapsed one.  Runs the default (non-mom6) hvel build so hvel(2)
      !! is an EXACT 0.5·(0 + 0); no drag/stress needed — the interior
      !! coupling alone hits the zero denominator.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(barotropic_workstate_t) :: bt
      integer, parameter :: NZ = 4
      real(wp), parameter :: DT = 450.0_wp
      real(wp) :: h_prof(NZ)
      integer :: k, ig, jg
      logical :: bounded, no_nan
      checks: block
         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         call bt%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.2_wp    ! kv > 0 at every interface (incl. the collapse)
         vd%hvel_mom6 = .false.      ! default build ⇒ hvel(2) = 0.5*(0+0) = EXACT 0
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ! Interior layer k = 2 exactly collapsed; thick neighbours around it.
         h_prof = [10.0_wp, 0.0_wp, 12.0_wp, 20.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_prof(k)
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, rho0=RHO0, &
                                   visc_rem_u=bt%visc_rem_u, visc_rem_v=bt%visc_rem_v)
         call map_out(ms, vd)

         bounded = .true.; no_nan = .true.
         do k = 1, NZ
            if (.not. ieee_is_finite(bt%visc_rem_u(ig, jg, k)) .or. &
                .not. ieee_is_finite(bt%visc_rem_v(ig, jg, k))) no_nan = .false.
            if (bt%visc_rem_u(ig, jg, k) <= 0.0_wp .or. &
                bt%visc_rem_u(ig, jg, k) > 1.0_wp + 1.0e-12_wp) bounded = .false.
            if (bt%visc_rem_v(ig, jg, k) <= 0.0_wp .or. &
                bt%visc_rem_v(ig, jg, k) > 1.0_wp + 1.0e-12_wp) bounded = .false.
         end do

         call check(error, no_nan, &
                    "collapsed layer: visc_rem NaN/Inf (hf_k denominator unfloored?)")
         if (allocated(error)) exit checks
         call check(error, bounded, &
                    "collapsed layer: violated 0 < gamma <= 1")
      end block checks
      call bt%destroy(); call vd%destroy(); call ms%destroy()
   end subroutine test_collapsed_layer_finite

end module test_ocean_vdiff_bbl
