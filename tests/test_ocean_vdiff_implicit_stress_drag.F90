!! Unit tests for folding the surface wind stress + bottom drag into the
!! ocean backward-Euler vertical-friction tridiagonal
!! (`rdb_ocean_vdiff` implicit_stress / implicit_drag).
!!
!! Convention reminder (Roundabout is bottom-up): k=1 is the BED (drag → k=1
!! diagonal), k=nz is the SURFACE (wind stress → k=nz RHS) — the MIRROR of
!! MOM6's top-down k=1/k=nz roles.
!!
!! Cases:
!!   (a) THIN-LAYER DRAG STABILITY — thin bed layer + strong quadratic drag
!!       at a dt where the explicit forward-Euler factor (1 - dt·c_d|U|/h)
!!       reverses sign and amplifies (dt·c_d|U|/h = 6).  The implicit fold
!!       must keep the bed velocity in (0, uⁿ) and decay monotonically over
!!       many steps; the explicit reference blows up / flips sign.
!!   (b) WELL-RESOLVED EQUIVALENCE — thick layer + small drag/stress: the
!!       implicit boundary terms reduce to the explicit increments in the
!!       zero-coupling / small-dt·λ limit.  Match to tolerance.
!!   (c) DEFAULT-OFF BIT-IDENTITY — knobs off ⇒ the momentum solve is
!!       bit-for-bit identical whether or not the fold args are supplied.
!!   (d) STRESS PLACEMENT + SIGN — eastward τ_x with implicit_stress lands
!!       the momentum at the SURFACE (k=nz), not the bed, with the right
!!       sign, and (with interior K_v) diffuses downward.
module test_ocean_vdiff_implicit_stress_drag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_momentum
   implicit none
   private

   public :: collect_ocean_vdiff_implicit_stress_drag_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 6
   real(wp), parameter :: RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_vdiff_implicit_stress_drag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("vdiff_implicit_drag_thin_layer_bounded", &
                               test_thin_layer_drag_bounded), &
                  new_unittest("vdiff_implicit_well_resolved_equivalence", &
                               test_well_resolved_equivalence), &
                  new_unittest("vdiff_implicit_default_off_bit_identity", &
                               test_default_off_bit_identity), &
                  new_unittest("vdiff_implicit_stress_placement_and_sign", &
                               test_stress_placement_and_sign), &
                  new_unittest("vdiff_implicit_stress_land_masked", &
                               test_stress_land_masked) &
                  ]
   end subroutine collect_ocean_vdiff_implicit_stress_drag_tests

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

   subroutine test_thin_layer_drag_bounded(error)
      !! Thin bed layer + strong quadratic drag: dt·c_d·|U|/h_1 = 6.  The
      !! explicit forward-Euler bed update u·(1 - 6) = -5u flips sign and
      !! grows; the implicit diagonal fold u/(1 + 6) decays monotonically
      !! and never leaves (0, uⁿ).  Drive vdiff_apply_momentum with
      !! K_v = 0 (no interior viscosity) so only the bed-drag diagonal acts,
      !! refreshing the Rayleigh rate λ = c_d·|U|/h each step (Picard).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H1 = 2.0_wp, CD = 2.5e-3_wp, U0 = 8.0_wp
      real(wp), parameter :: DT = 600.0_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :), tau_u(:, :), tau_v(:, :)
      real(wp) :: u_bed, u_prev, speed, lam_explicit, u_explicit
      integer :: nx, ny, nu, nv_u, nx_v, nv, step, ig, jg
      logical :: monotone, bounded, explicit_blew_up
      checks: block

         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.0_wp        ! isolate the bed-drag diagonal
         vd%implicit_drag = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         nu = size(ms%u_face_x_layer, 1)
         nv_u = size(ms%u_face_x_layer, 2)
         nx_v = size(ms%v_face_y_layer, 1)
         nv = size(ms%v_face_y_layer, 2)

         ms%h_layer = H1
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp

         allocate (lam_u(nu, nv_u), source=0.0_wp)
         allocate (lam_v(nx_v, nv), source=0.0_wp)
         allocate (tau_u(nu, nv_u), source=0.0_wp)
         allocate (tau_v(nx_v, nv), source=0.0_wp)

         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ! Reference: explicit forward-Euler bed factor at step 1.
         lam_explicit = CD*U0/H1            ! = 0.01 (1/s)
         u_explicit = U0*(1.0_wp - DT*lam_explicit)   ! = 8·(1-6) = -40
         explicit_blew_up = (u_explicit < 0.0_wp) .and. (abs(u_explicit) > U0)
         call check(error, explicit_blew_up, &
                    "test setup: explicit forward-Euler should reverse + amplify here")
         if (allocated(error)) exit checks

         call map_in(ms, vd)
         monotone = .true.
         bounded = .true.
         u_prev = U0
         do step = 1, 50
            ! Picard: freeze |U| at the current bed velocity (v=0 here).
            !$acc update self(ms%u_face_x_layer)
            speed = abs(ms%u_face_x_layer(ig, jg, 1))
            lam_u = CD*speed/H1
            ! lam_u/lam_v/tau_u/tau_v are unmapped local host arrays — the
            ! bare `do concurrent` in diffuse_velocity_columns_impl auto-copies
            ! them in per launch (stdpar), picking up this iteration's fresh
            ! Picard lam_u.  No `!$acc update device` (they are not in the
            ! present table; an update on unmapped data is a fatal error).
            call vdiff_apply_momentum(grid, vd, ms, DT, &
                                      tau_u=tau_u, tau_v=tau_v, &
                                      lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
            !$acc update self(ms%u_face_x_layer)
            u_bed = ms%u_face_x_layer(ig, jg, 1)
            if (u_bed < 0.0_wp .or. u_bed > U0 + 1.0e-12_wp) bounded = .false.
            if (u_bed > u_prev + 1.0e-12_wp) monotone = .false.
            if (u_bed /= u_bed) bounded = .false.   ! NaN
            u_prev = u_bed
         end do
         call map_out(ms, vd)

         call check(error, bounded, &
                    "implicit drag: bed velocity left (0, uⁿ) (overshoot / sign flip / NaN)")
         if (allocated(error)) exit checks
         call check(error, monotone, &
                    "implicit drag: bed velocity not monotonically decaying")
         if (allocated(error)) exit checks
         call check(error, u_prev > 0.0_wp .and. u_prev < U0, &
                    "implicit drag: final bed velocity should be decayed but positive")

      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v, tau_u, tau_v)
      call vd%destroy()
      call ms%destroy()
   end subroutine test_thin_layer_drag_bounded

   subroutine test_well_resolved_equivalence(error)
      !! Thick layer (h=100 m), small drag (dt·c_d|U|/h = 1e-3) and small
      !! stress: the implicit boundary terms agree with the explicit
      !! increments to O((dt·λ)²) for drag and EXACTLY for stress (K_v = 0 ⇒
      !! the surface row is an identity operator on the RHS add).  One step.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 100.0_wp, U0 = 0.5_wp, DT = 1.0_wp
      real(wp), parameter :: LAM = 1.0e-3_wp        ! dt·λ = 1e-3
      real(wp), parameter :: TAUX = 0.1_wp          ! N/m^2
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :), tau_u(:, :), tau_v(:, :)
      real(wp) :: u_bed_imp, u_surf_imp, u_bed_exp, u_surf_exp
      integer :: nu, nv_u, nx_v, nv, ig, jg
      checks: block

         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.0_wp
         vd%implicit_drag = .true.
         vd%implicit_stress = .true.
         nu = size(ms%u_face_x_layer, 1)
         nv_u = size(ms%u_face_x_layer, 2)
         nx_v = size(ms%v_face_y_layer, 1)
         nv = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ms%h_layer = H0
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp

         allocate (lam_u(nu, nv_u), source=LAM)       ! λ = 1e-3 1/s
         allocate (lam_v(nx_v, nv), source=LAM)
         allocate (tau_u(nu, nv_u), source=TAUX)
         allocate (tau_v(nx_v, nv), source=0.0_wp)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         !$acc update self(ms%u_face_x_layer)
         u_bed_imp = ms%u_face_x_layer(ig, jg, 1)
         u_surf_imp = ms%u_face_x_layer(ig, jg, NZ)
         call map_out(ms, vd)

         ! Explicit references.
         u_bed_exp = U0*(1.0_wp - DT*LAM)                       ! drag
         u_surf_exp = U0 + DT*TAUX/(RHO0*H0)                    ! stress

         call check(error, abs(u_bed_imp - u_bed_exp) < 1.0e-6_wp*abs(U0), &
                    "implicit drag: not equivalent to explicit in the small-dt·λ limit")
         if (allocated(error)) exit checks
         call check(error, abs(u_surf_imp - u_surf_exp) < 1.0e-12_wp, &
                    "implicit stress: should match explicit increment exactly at K_v=0")

      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v, tau_u, tau_v)
      call vd%destroy()
      call ms%destroy()
   end subroutine test_well_resolved_equivalence

   subroutine test_default_off_bit_identity(error)
      !! Knobs off ⇒ the matrix build + solve is unchanged whether or not
      !! the fold arrays are supplied.  Two states: one solved with NO fold
      !! args, one with the fold args present but vd%implicit_* = .false.
      !! (and non-trivial λ/τ that MUST be ignored).  Results bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :), tau_u(:, :), tau_v(:, :)
      real(wp) :: max_du
      integer :: k, nu, nv_u, nx_v, nv
      checks: block

         call make_grid(grid, 6, 4)
         ms_a%nz_ml = NZ
         ms_b%nz_ml = NZ
         call ms_a%init(grid)
         call ms_b%init(grid)
         call vd_a%init(grid, nz_ml=NZ)
         call vd_b%init(grid, nz_ml=NZ)
         vd_a%K_v_momentum = 1.0e-2_wp
         vd_b%K_v_momentum = 1.0e-2_wp
         vd_a%implicit_drag = .false.
         vd_a%implicit_stress = .false.
         vd_b%implicit_drag = .false.       ! OFF — fold args must be ignored
         vd_b%implicit_stress = .false.
         nu = size(ms_a%u_face_x_layer, 1)
         nv_u = size(ms_a%u_face_x_layer, 2)
         nx_v = size(ms_a%v_face_y_layer, 1)
         nv = size(ms_a%v_face_y_layer, 2)

         ms_a%h_layer = H0
         ms_b%h_layer = H0
         do k = 1, NZ
            ms_a%u_face_x_layer(:, :, k) = 0.1_wp*real(k, wp)
            ms_b%u_face_x_layer(:, :, k) = 0.1_wp*real(k, wp)
            ms_a%v_face_y_layer(:, :, k) = -0.05_wp*real(k, wp)
            ms_b%v_face_y_layer(:, :, k) = -0.05_wp*real(k, wp)
         end do

         allocate (lam_u(nu, nv_u), source=5.0_wp)   ! huge λ/τ — MUST be ignored
         allocate (lam_v(nx_v, nv), source=5.0_wp)
         allocate (tau_u(nu, nv_u), source=10.0_wp)
         allocate (tau_v(nx_v, nv), source=10.0_wp)

         call map_in(ms_a, vd_a)
         call map_in(ms_b, vd_b)
         call vdiff_apply_momentum(grid, vd_a, ms_a, DT)
         call vdiff_apply_momentum(grid, vd_b, ms_b, DT, &
                                   tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         call map_out(ms_a, vd_a)
         call map_out(ms_b, vd_b)

         max_du = maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer)) + &
                  maxval(abs(ms_a%v_face_y_layer - ms_b%v_face_y_layer))
         call check(error, max_du == 0.0_wp, &
                    "default-off: fold args changed the solve (must be bit-identical)")

      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v, tau_u, tau_v)
      call vd_a%destroy()
      call vd_b%destroy()
      call ms_a%destroy()
      call ms_b%destroy()
   end subroutine test_default_off_bit_identity

   subroutine test_stress_placement_and_sign(error)
      !! Eastward τ_x with implicit_stress: the momentum source lands at the
      !! SURFACE (k=nz) with a POSITIVE increment.  With K_v = 0 the bed
      !! (k=1) is untouched (pins placement at k=nz, not k=1).  With interior
      !! K_v > 0 the surface momentum diffuses downward (surface still gains
      !! the most), confirming the RHS term feeds the implicit operator.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 600.0_wp, TAUX = 0.2_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :), tau_u(:, :), tau_v(:, :)
      real(wp) :: u_bed, u_surf, u_mid
      integer :: nu, nv_u, nx_v, nv, ig, jg
      checks: block

         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.0_wp
         vd%implicit_stress = .true.
         nu = size(ms%u_face_x_layer, 1)
         nv_u = size(ms%u_face_x_layer, 2)
         nx_v = size(ms%v_face_y_layer, 1)
         nv = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         allocate (lam_u(nu, nv_u), source=0.0_wp)
         allocate (lam_v(nx_v, nv), source=0.0_wp)
         allocate (tau_u(nu, nv_u), source=TAUX)
         allocate (tau_v(nx_v, nv), source=0.0_wp)

         ! Pass 1: K_v = 0 — momentum only at the surface row.
         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         !$acc update self(ms%u_face_x_layer)
         u_bed = ms%u_face_x_layer(ig, jg, 1)
         u_surf = ms%u_face_x_layer(ig, jg, NZ)
         call map_out(ms, vd)

         call check(error, u_surf > 0.0_wp, &
                    "implicit stress: eastward τ_x must give positive surface u")
         if (allocated(error)) exit checks
         call check(error, abs(u_bed) < 1.0e-14_wp, &
                    "implicit stress: bed (k=1) must be untouched at K_v=0 (placement at k=nz)")
         if (allocated(error)) exit checks
         call check(error, abs(u_surf - DT*TAUX/(RHO0*H0)) < 1.0e-12_wp, &
                    "implicit stress: surface increment ≠ dt·τ/(ρ₀·h)")
         if (allocated(error)) exit checks

         ! Pass 2: K_v > 0 — momentum diffuses downward; surface > mid > 0.
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         vd%K_v_momentum = 1.0e-1_wp
         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         !$acc update self(ms%u_face_x_layer)
         u_surf = ms%u_face_x_layer(ig, jg, NZ)
         u_mid = ms%u_face_x_layer(ig, jg, NZ/2)
         u_bed = ms%u_face_x_layer(ig, jg, 1)
         call map_out(ms, vd)

         call check(error, u_surf > u_mid .and. u_mid > 0.0_wp .and. u_bed > 0.0_wp, &
                    "implicit stress: momentum should diffuse downward (surface > mid > 0)")

      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v, tau_u, tau_v)
      call vd%destroy()
      call ms%destroy()
   end subroutine test_stress_placement_and_sign

   subroutine test_stress_land_masked(error)
      !! No-normal-flow land fidelity: the implicit wind-stress fold must
      !! mask `tau_face` by the face `min` of the two bounding cells' wet
      !! mask (as the explicit surface-stress kernel does), so a land face
      !! receives NO stress.  Mark cell (ig,jg) land; with K_v = 0 and a
      !! uniform eastward τ_x, the U-face at i=ig (cells ig-1 | ig) has
      !! face-mask min(1,0)=0 ⇒ zero surface velocity, while the all-wet
      !! U-face at i=ig-1 (cells ig-2 | ig-1) gets the full dt·τ/(ρ₀·h).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 50.0_wp, TAUX = 0.1_wp, DT = 1.0_wp
      real(wp), allocatable :: lam_u(:, :), lam_v(:, :), tau_u(:, :), tau_v(:, :)
      real(wp) :: u_land, u_wet
      integer :: nu, nv_u, nx_v, nv, ig, jg
      checks: block

         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.0_wp
         vd%implicit_stress = .true.
         nu = size(ms%u_face_x_layer, 1)
         nv_u = size(ms%u_face_x_layer, 2)
         nx_v = size(ms%v_face_y_layer, 1)
         nv = size(ms%v_face_y_layer, 2)
         ig = 1 + grid%nghost
         jg = 1 + grid%nghost

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Cell (ig,jg) is land; the U-face at i=ig (between cells ig-1 and
         ! ig) is therefore a land face.
         ms%wet_mask = 1.0_wp
         ms%wet_mask(ig, jg) = 0.0_wp

         allocate (lam_u(nu, nv_u), source=0.0_wp)
         allocate (lam_v(nx_v, nv), source=0.0_wp)
         allocate (tau_u(nu, nv_u), source=TAUX)
         allocate (tau_v(nx_v, nv), source=0.0_wp)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT, &
                                   tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=lam_u, lambda_bot_v=lam_v, rho0=RHO0)
         !$acc update self(ms%u_face_x_layer)
         u_land = ms%u_face_x_layer(ig, jg, NZ)
         u_wet = ms%u_face_x_layer(ig - 1, jg, NZ)
         call map_out(ms, vd)

         call check(error, abs(u_land) < 1.0e-14_wp, &
                    "implicit stress: land face (masked) must receive no stress")
         if (allocated(error)) exit checks
         call check(error, abs(u_wet - DT*TAUX/(RHO0*H0)) < 1.0e-12_wp, &
                    "implicit stress: all-wet face must receive full dt·τ/(ρ₀·h)")

      end block checks
      if (allocated(lam_u)) deallocate (lam_u, lam_v, tau_u, tau_v)
      call vd%destroy()
      call ms%destroy()
   end subroutine test_stress_land_masked

end module test_ocean_vdiff_implicit_stress_drag
