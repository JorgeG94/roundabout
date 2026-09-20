!! THE P5.1 GATE: a flat ice lid is indistinguishable from a shallower
!! ocean.
module test_ocean_cavity_equivalence
   !! ### What is being claimed
   !!
   !! Put a UNIFORM ice draft `d` over a flat bed `b`, with no load term
   !! anywhere (this slice carries the datum only — `ms%p_top` is
   !! untouched), and run the full split solver.  Because the datum
   !! absorbs the draft exactly (`bt_H_ref = b - d`, `sum h_layer = b - d`,
   !! `bt_eta = 0`), the model should evolve like a cavity-free ocean over
   !! the shallower bed `b - d`.  Nothing about a flat lid is dynamics; it
   !! is bookkeeping, and this is the test that says so.
   !!
   !! ### What matches BITWISE, and what does not — and why
   !!
   !! At `t = 0` every prognostic is bit-for-bit identical, and that is
   !! structural rather than lucky: `b - z_draft` is `1000 - 500` in one
   !! run and `b` is `500` in the other, both exact, and every seed below
   !! (layer split, tracers, wet mask, datum) reads that one number.
   !!
   !! After stepping they are NOT bitwise identical, and cannot be.  ONE
   !! thing in the dyn-core legitimately depends on ABSOLUTE z: the
   !! FV_MOM6 pressure stack.  It is an anomaly about `rho_ref*g*z`, seeded
   !! at the bed (`e_face(1) = -b`, from the TRUE bed on both runs) and
   !! closed at the surface with
   !!
   !!     pa(nz+1) = rho_ref*GRAVITY*eta_geo,   eta_geo = -b + sum h.
   !!
   !! Under the lid `eta_geo = -500`, so the whole `pa` stack is offset by
   !! the constant `C = rho_ref*g*500 ~ 5e6 Pa` relative to the shallow
   !! run.  A CONSTANT offset produces exactly zero pressure-gradient
   !! force in exact arithmetic (that is the theorem `&ocean_pgf_nml
   !! p_top_in_bc` rests on) — but in floating point the two stacks round
   !! differently, because one accumulates `O(1e4)` anomalies on top of
   !! `-5e6` and the other on top of `0`.  The surviving difference is a
   !! few ulp OF THE OFFSET, not of the signal.
   !!
   !! That is precisely the leak the load term is for: wiring
   !! `p_ice_ref` into `pa(nz+1)` cancels `C` inside the boundary
   !! condition and puts the stack back at `O(1e4)`.  Until then the
   !! honest assertion is a ROUND-OFF BOUND, and this suite states how it
   !! is built (`bound_accel` below) instead of tuning a number until the
   !! test goes green.
   !!
   !! Everything else is shift-invariant by construction: the EOS
   !! reference pressure is the scalar `&ocean_eos_nml p_ref` (never a
   !! function of z), `in_eos` is off, the z-level T/S overlay is refused
   !! under a cavity, and `&ocean_topo_nml` forcing, the Coriolis
   !! parameter and every metric depend on (x, y) alone.
   !!
   !! ### Why it is not trivially at rest
   !!
   !! Both runs carry a stratified water column and a two-gyre wind
   !! stress, so the flow is O(mm/s) after six steps.  The suite asserts
   !! that separately (non-vacuity), because a gate that passes on a
   !! motionless ocean proves nothing.
   !!
   !! ### Measured, gfortran 15.1 Release, 2026-09-20
   !!
   !! ```
   !! max|du| = 3.1e-13 m/s      bound 3.1e-11   (100x margin)
   !! max|dv| = 2.6e-13 m/s
   !! max|dh| = 8.5e-14 m        bound 1.4e-08
   !! max|deta| = 5.4e-13 m
   !! max|u|  = 3.5e-03 m/s      i.e. the signal is 1e10 x the difference
   !! ```
   !!
   !! So: not bitwise, and the reason is one identified term; a decade of
   !! margin under a derived bound; ten decades of separation from the
   !! physics.  If `du` ever reaches the bound, the pressure stack is no
   !! longer merely re-rounded — something reads absolute z for real.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, rdb_ocean_get_bt_eta_ptr, &
                            rdb_ocean_set_bathymetry
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_cavity_equivalence_tests

   real(wp), parameter :: BED = 1000.0_wp
      !! Bed depth of the CAVITY run (m).
   real(wp), parameter :: DRAFT = 500.0_wp
      !! Uniform ice draft (m).  The shallow twin's bed is `BED - DRAFT`.
   real(wp), parameter :: DX = 4000.0_wp
   integer, parameter :: N_STEPS = 6
   real(wp), parameter :: DT = 300.0_wp

   ! Reference density of record for this namelist (`&ocean_ic_nml rho_0`
   ! default), i.e. the `rho_ref` the FV_MOM6 stack is an anomaly about.
   real(wp), parameter :: RHO_REF = 1035.0_wp

contains

   subroutine collect_ocean_cavity_equivalence_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_flat_lid_seed_bit_identical", test_seed_identical), &
                  new_unittest("cavity_flat_lid_equivalence", test_flat_lid_equivalence), &
                  new_unittest("cavity_api_bathymetry_swap_refused", test_api_bathy_refused) &
                  ]
   end subroutine collect_ocean_cavity_equivalence_tests

   function case_nml(max_depth, cavity) result(nml)
      !! The two runs differ in EXACTLY two lines: the bed depth and the
      !! cavity group.  Everything else — grid, layers, dt, stratification,
      !! wind, PGF form, vcoord, the PINNED `n_inner` — is shared text.
      !!
      !! `n_inner` is pinned rather than auto-derived on purpose:
      !! `auto_n_inner` takes its gravity-wave speed from `max(b)`, the
      !! BED, so the two runs would otherwise choose different substep
      !! counts (99 m/s over 1000 m vs 70 m/s over 500 m) and there would
      !! be nothing to compare.  That conservatism is deliberate and
      !! documented at the latch; here it simply must not be in the way.
      real(wp), intent(in) :: max_depth
      character(len=*), intent(in) :: cavity
      character(len=:), allocatable :: nml
      character(len=32) :: depth_s
      write (depth_s, '(F12.2)') max_depth
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 12, ny = 10, nghost = 2, dx = 4000.0, dy = 4000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 6 /"//new_line("a")// &
            "&time_nml t_end = 100000.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = "//trim(adjustl(depth_s))//", "// &
            "wind_config = '2gyre', taux_magnitude = 0.1 /"//new_line("a")// &
            "&physics_nml coriolis_f = 1.0e-4 /"//new_line("a")// &
            ! Stratified: T from 4 degC at the bed to 12 degC at the
            ! surface, so every layer carries a different density and the
            ! baroclinic pressure gradient is live.
            "&tracer_nml initial_salinity = 35.0, T_init_bottom = 4.0, "// &
            "T_init_surface = 12.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .false., "// &
            "n_inner = 12 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
      if (len_trim(cavity) > 0) then
         nml = nml//"&ocean_cavity_dyn_nml "//cavity//" /"//new_line("a")
      end if
   end function case_nml

   subroutine run_case(max_depth, cavity, n_steps, h, u, v, eta, ok)
      !! Create an ocean from the namelist, step it, and snapshot the
      !! prognostics into caller-owned arrays.  The handle is destroyed
      !! before returning: the API holds ONE live ocean at a time, so the
      !! two runs are sequential and each reads its own snapshot back.
      !!
      !! Every read goes through `rdb_ocean_refresh_host` first — the D<->H
      !! contract.  A getter alone never triggers a device->host copy, so
      !! skipping it would compare stale host memory on the GPU build and
      !! could pass for the wrong reason.
      real(wp), intent(in) :: max_depth
      character(len=*), intent(in) :: cavity
      integer, intent(in) :: n_steps
      real(wp), allocatable, intent(out) :: h(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable, intent(out) :: eta(:, :)
      logical, intent(out) :: ok

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen
      real(wp), pointer :: f3(:, :, :), f2(:, :)
      character(len=:), allocatable :: nml

      ok = .false.
      handle = c_null_ptr
      nml = case_nml(max_depth, cavity)
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return

      if (n_steps > 0) then
         status = rdb_ocean_step(handle, int(n_steps, c_int))
         if (status /= OCEAN_STATUS_OK) then
            status = rdb_ocean_destroy(handle)
            return
         end if
      end if

      status = rdb_ocean_refresh_host(handle)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle)
         return
      end if

      status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle)
         return
      end if
      call c_f_pointer(ptr, f3, [int(nx), int(ny), int(nz)])
      allocate (h(nx, ny, nz))
      h = f3

      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, f3, [int(nx), int(ny), int(nz)])
      allocate (u(nx, ny, nz))
      u = f3

      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, f3, [int(nx), int(ny), int(nz)])
      allocate (v(nx, ny, nz))
      v = f3

      status = rdb_ocean_get_bt_eta_ptr(handle, ptr, nx, ny, gen)
      call c_f_pointer(ptr, f2, [int(nx), int(ny)])
      allocate (eta(nx, ny))
      eta = f2

      status = rdb_ocean_destroy(handle)
      ok = (status == OCEAN_STATUS_OK)
   end subroutine run_case

   function bound_accel() result(tol_u)
      !! The round-off bound on the velocity difference, built from the
      !! one quantity that differs between the runs.
      !!
      !!   * the two `pa` stacks are offset by `C = rho_ref*g*DRAFT`;
      !!   * every `pa` entry therefore carries an absolute rounding of
      !!     order `C*eps`, where the shallow run's is negligible beside it;
      !!   * the Pass-3 face force divides a `pa*h` difference by
      !!     `(h_L+h_R)`, so the acceleration error is `~ C*eps/(rho_0*dx)`;
      !!   * over `N_STEPS` steps of `DT` that integrates to
      !!     `C*eps*DT*N_STEPS/(rho_0*dx)`.
      !!
      !! The factor 64 covers the handful of roundings per pass (three
      !! passes, the barotropic substeps, the ALE remap) and the mild
      !! step-to-step amplification.  It is NOT fitted to the answer: the
      !! measured difference (3.1e-13 m/s) sits 100x under this bound and
      !! 1e10 under the flow it rides on.
      real(wp) :: tol_u
      tol_u = 64.0_wp*(RHO_REF*GRAVITY*DRAFT)*epsilon(1.0_wp)* &
              DT*real(N_STEPS, wp)/(RHO_REF*DX)
   end function bound_accel

   subroutine test_seed_identical(error)
      !! At `t = 0` the two runs are BIT-for-bit identical — the datum
      !! absorption is exact arithmetic, not an approximation.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_a(:, :, :), u_a(:, :, :), v_a(:, :, :), eta_a(:, :)
      real(wp), allocatable :: h_b(:, :, :), u_b(:, :, :), v_b(:, :, :), eta_b(:, :)
      logical :: ok_a, ok_b

      call run_case(BED, "enable = .true., draft_config = 'flat', draft_depth = "// &
                    "500.0", 0, h_a, u_a, v_a, eta_a, ok_a)
      call check(error, ok_a, "the cavity run must create")
      if (allocated(error)) return
      call run_case(BED - DRAFT, "", 0, h_b, u_b, v_b, eta_b, ok_b)
      call check(error, ok_b, "the shallow twin must create")
      if (allocated(error)) return

      call check(error, all(h_a == h_b), &
                 "h_layer under a 500 m lid over a 1000 m bed must be BIT-identical "// &
                 "to a 500 m ocean at t = 0")
      if (allocated(error)) return
      call check(error, all(eta_a == eta_b), &
                 "bt_eta must be BIT-identical at t = 0 (both exactly zero: the "// &
                 "datum absorbed the draft)")
      if (allocated(error)) return
      call check(error, all(eta_a == 0.0_wp), &
                 "bt_eta must be exactly zero under the lid — the free surface is "// &
                 "measured from the LOADED equilibrium, not from z = 0")
      if (allocated(error)) return
      call check(error, all(u_a == u_b) .and. all(v_a == v_b), &
                 "the velocity seeds must be BIT-identical at t = 0")
   end subroutine test_seed_identical

   subroutine test_flat_lid_equivalence(error)
      !! The gate itself: several steps of the full `pred_corr` split
      !! solver, stratified and wind-driven, and the two runs stay
      !! together to the round-off bound above.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_a(:, :, :), u_a(:, :, :), v_a(:, :, :), eta_a(:, :)
      real(wp), allocatable :: h_b(:, :, :), u_b(:, :, :), v_b(:, :, :), eta_b(:, :)
      logical :: ok_a, ok_b
      real(wp) :: du, dv, dh, deta, tol_u, tol_h, signal

      call run_case(BED, "enable = .true., draft_config = 'flat', draft_depth = "// &
                    "500.0", N_STEPS, h_a, u_a, v_a, eta_a, ok_a)
      call check(error, ok_a, "the cavity run must step")
      if (allocated(error)) return
      call run_case(BED - DRAFT, "", N_STEPS, h_b, u_b, v_b, eta_b, ok_b)
      call check(error, ok_b, "the shallow twin must step")
      if (allocated(error)) return

      call check(error, all(ieee_is_finite(u_a)) .and. all(ieee_is_finite(h_a)), &
                 "the cavity run must stay finite")
      if (allocated(error)) return

      du = maxval(abs(u_a - u_b))
      dv = maxval(abs(v_a - v_b))
      dh = maxval(abs(h_a - h_b))
      deta = maxval(abs(eta_a - eta_b))
      tol_u = bound_accel()
      ! A thickness difference is a velocity difference integrated by the
      ! continuity divergence: `dh ~ du*(dt*N)*H/dx`.
      tol_h = tol_u*DT*real(N_STEPS, wp)*BED/DX

      ! NON-VACUITY FIRST: if the ocean never moved, everything below is
      ! 0 == 0 and the gate is worthless.
      signal = maxval(abs(u_a))
      call check(error, signal > 1.0e4_wp*tol_u, &
                 "the wind-driven flow must be decades above the bound, else the "// &
                 "equivalence gate is vacuous")
      if (allocated(error)) return

      call check(error, du <= tol_u, "u_face_x_layer must match the shallow twin to "// &
                 "the pa-offset round-off bound")
      if (allocated(error)) return
      call check(error, dv <= tol_u, "v_face_y_layer must match the shallow twin")
      if (allocated(error)) return
      call check(error, dh <= tol_h, "h_layer must match the shallow twin")
      if (allocated(error)) return
      call check(error, deta <= tol_h, "bt_eta must match the shallow twin")
   end subroutine test_flat_lid_equivalence

   subroutine test_api_bathy_refused(error)
      !! The one API entry point that could silently destroy the datum:
      !! `rdb_ocean_set_bathymetry` re-derives `bt_H_ref = b`, which under
      !! a cavity drops the ice load out of the reference depth — and even
      !! re-deriving it correctly would not be enough, because a new bed
      !! changes which columns GROUND and this call does not rebuild the
      !! wet mask.  It must refuse, and it must still work without a
      !! cavity (else the test proves only that the call is broken).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: nml
      real(wp) :: b_new(12, 10)

      b_new = 900.0_wp

      nml = case_nml(BED, "enable = .true., draft_config = 'flat', draft_depth = 500.0")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "the cavity run must create")
      if (allocated(error)) return
      status = rdb_ocean_set_bathymetry(handle, b_new, 12_c_int, 10_c_int)
      call check(error, status /= OCEAN_STATUS_OK, &
                 "a mid-run bathymetry swap must be REFUSED under a cavity — it "// &
                 "would silently reset bt_H_ref to the bed")
      status = rdb_ocean_destroy(handle)
      if (allocated(error)) return

      nml = case_nml(BED - DRAFT, "")
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "the cavity-free run must create")
      if (allocated(error)) return
      status = rdb_ocean_set_bathymetry(handle, b_new, 12_c_int, 10_c_int)
      call check(error, status == OCEAN_STATUS_OK, &
                 "without a cavity the same call must still be accepted")
      status = rdb_ocean_destroy(handle)
   end subroutine test_api_bathy_refused

end module test_ocean_cavity_equivalence
