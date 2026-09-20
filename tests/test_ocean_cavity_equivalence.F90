!! THE CAVITY EQUIVALENCE GATE: a flat ice lid is indistinguishable from
!! a shallower ocean — unloaded (P5.1) and loaded (P5.2).
module test_ocean_cavity_equivalence
   !! ### What is being claimed
   !!
   !! Put a UNIFORM ice draft `d` over a flat bed `b` and run the full
   !! split solver.  Because the datum
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
   !! condition and puts the stack back on its anomaly scale.  The
   !! honest assertion either way is a ROUND-OFF BOUND, and this suite
   !! states how each one is built (`bound_accel`,
   !! `bound_accel_loaded`) instead of tuning a number until the test
   !! goes green.  BOTH variants ship: the unloaded one documents the
   !! offset (a UNIFORM draft is the one cavity configuration that does
   !! not require `&ocean_pgf_nml p_top_in_bc`, because a load with no
   !! gradient is provably inert in the top BC), the loaded one is the
   !! P5.2 gate.
   !!
   !! Everything else is shift-invariant by construction: the EOS
   !! reference pressure is the scalar `&ocean_eos_nml p_ref` (never a
   !! function of z), `in_eos` is off, the z-level T/S overlay is refused
   !! under a cavity, and `&ocean_topo_nml` forcing, the Coriolis
   !! parameter and every metric depend on (x, y) alone.
   !!
   !! ### Why it is not trivially at rest
   !!
   !! Both runs carry a stratified water column and an identical INITIAL
   !! zonal jet (`U_SEED`, written through `rdb_ocean_set_u`), so the
   !! flow is O(mm/s) through the whole window — nine decades above the
   !! tighter of the two bounds.  The suite asserts that separately
   !! (non-vacuity), because a gate that passes on a motionless ocean
   !! proves nothing.
   !!
   !! The stirrer used to be a two-gyre WIND stress and can no longer
   !! be: P2c masks every atmospheric forcing term with
   !! `1 - cover_frac`, a uniform flat lid covers the whole domain, and
   !! so the cavity run correctly feels no wind while its uncovered twin
   !! feels all of it.  See `U_SEED`.
   !!
   !! ### Measured, gfortran 15.1 Release, 2026-09-20
   !!
   !! (Re-measured when the stirrer changed from the two-gyre wind to the
   !! `U_SEED` initial jet — see below.  Same shape, same decade counts;
   !! the LOADED margin tightened from 24x to 6x because a jet dropped
   !! into an unbalanced stratified column excites the pressure stack
   !! somewhat harder than a six-step wind spin-up did.)
   !!
   !! ```
   !!                      UNLOADED (P5.1)        LOADED (P5.2)
   !!   max|du|            2.854e-13 m/s          4.987e-18 m/s
   !!   max|dv|            2.645e-13 m/s          3.578e-18 m/s
   !!   max|dh|            5.684e-14 m            0 (exactly)
   !!   max|deta|          4.547e-13 m            0 (exactly)
   !!   bound              3.136e-11              3.09e-17
   !!   margin                 110x                    6x
   !!   max|u| (signal)    1.559e-03 m/s          (same run)
   !! ```
   !!
   !! **The load buys 5.7e4x — 4.8 decades — on the velocity difference,
   !! and exact agreement on thickness and SSH.**  Inverting the bound's
   !! own algebra, the effective pressure-stack discrepancy falls from
   !! `7.2e-10 Pa` (= `eps * rho_ref*g*DRAFT`, i.e. exactly the offset's
   !! own rounding, which is the term the load cancels) to `3.0e-15 Pa`
   !! (= a few `eps * PA_ANOM`, the anomaly stack's own rounding — the
   !! floor, with no offset left to remove).
   !!
   !! So: still not bitwise, and the reason is still one identified term;
   !! a derived bound with margin in both variants; and after the load,
   !! fifteen decades of separation from the physics.  If `du` in the
   !! LOADED variant ever reaches its bound, the cancellation inside
   !! `pa(nz+1)` has stopped happening.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, rdb_ocean_get_bt_eta_ptr, &
                            rdb_ocean_set_bathymetry, rdb_ocean_set_u, &
                            rdb_ocean_get_grid_info
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
   real(wp), parameter :: ALPHA_T = 1.7e-4_wp
      !! `&ocean_ic_nml alpha_T` at its default (kg/m^3 per degC) — the
      !! namelist below does not set it.
   real(wp), parameter :: T_HALF_RANGE = 6.0_wp
      !! Largest `|T - T_ref|` in the IC: T runs 4..12 degC about the
      !! default `T_ref = 10`.  S is uniform at `S_ref`, so the density
      !! anomaly is thermal only.
   real(wp), parameter :: U_SEED = 1.0e-3_wp
      !! Amplitude of the INITIAL zonal jet (m/s) both runs are seeded
      !! with; the `1 - cos` profile peaks at `2*U_SEED`.
      !!
      !! This suite used to stir the basin with a `'2gyre'` wind stress.
      !! It cannot any more, and the reason is the feature under test one
      !! phase later: P2c masks every atmospheric forcing term with
      !! `1 - cover_frac`, and a uniform flat lid covers the WHOLE
      !! domain, so the cavity run correctly feels no wind at all while
      !! its shallow twin feels the full stress.  The gate would then be
      !! comparing two different experiments — and with the wind removed
      !! from both it would compare two oceans at rest, which is the
      !! vacuity the suite explicitly refuses to accept.
      !!
      !! An INITIAL velocity is the right stirrer for a datum-equivalence
      !! gate: it is a prognostic, not a forcing, so no cover mask can
      !! touch it; it is written through `rdb_ocean_set_u` with the SAME
      !! array in both runs, so the two experiments still differ in
      !! exactly the bed depth and the cavity group; and it excites the
      !! same baroclinic pressure stack the wind used to, because the jet
      !! is in no balance with the stratified column it is dropped into.
   real(wp), parameter :: P_SURF_SEAM = 2000.0_wp
      !! Atmospheric load (Pa) for the seam cases — ~20 hPa, i.e. a
      !! 0.197 m inverse-barometer elevation.  Small enough that the
      !! 5.08e6 Pa ice load would stand out by 2500x in `eta_ib` if it
      !! ever leaked onto the seam.

   ! Reference density of record for this namelist (`&ocean_ic_nml rho_0`
   ! default), i.e. the `rho_ref` the FV_MOM6 stack is an anomaly about.
   real(wp), parameter :: RHO_REF = 1035.0_wp

   real(wp), parameter :: PA_ANOM = (ALPHA_T*T_HALF_RANGE)*GRAVITY*(BED - DRAFT)
      !! Bound on `|pa|` once the load is cancelled at `pa(nz+1)`: the
      !! density anomaly integrated over the water column, ~5 Pa — the
      !! anomaly scale the FV_MOM6 stack is SUPPOSED to live on, against
      !! the 5.08e6 Pa the unloaded run carries.

contains

   subroutine collect_ocean_cavity_equivalence_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_flat_lid_seed_bit_identical", test_seed_identical), &
                  new_unittest("cavity_flat_lid_equivalence", test_flat_lid_equivalence), &
                  new_unittest("cavity_flat_lid_equivalence_loaded", &
                               test_flat_lid_equivalence_loaded), &
                  new_unittest("cavity_seam_zero_p_surf_inert", test_seam_zero_p_surf_inert), &
                  new_unittest("cavity_seam_matches_cavity_free", &
                               test_seam_matches_cavity_free), &
                  new_unittest("cavity_api_bathymetry_swap_refused", test_api_bathy_refused) &
                  ]
   end subroutine collect_ocean_cavity_equivalence_tests

   function case_nml(max_depth, cavity, loaded, psurf) result(nml)
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
      logical, intent(in), optional :: loaded
         !! `&ocean_pgf_nml p_top_in_bc` — the P5.2 load in the FV_MOM6
         !! `pa(nz+1)` surface BC.  Default `.false.` = the P5.1
         !! datum-only run this suite was written for.  A UNIFORM draft
         !! is the one cavity configuration that does not REQUIRE it (a
         !! load with no gradient is provably inert in the top BC), which
         !! is exactly why both variants are expressible here.
      real(wp), intent(in), optional :: psurf
         !! When present: `&ocean_psurf_nml enable` + the PR-12 component
         !! set, with `p_surf_const` set to this value (Pa), so the
         !! `eta_ib = -p_surf/(rho0 g_bt)` seam is LIVE.  Absent ⇒ no
         !! psurf group at all and no seam.
      character(len=:), allocatable :: nml
      character(len=32) :: depth_s, psurf_s
      character(len=:), allocatable :: pgf_l
      logical :: want_psurf

      pgf_l = "&ocean_pgf_nml form = 'fv_mom6' /"
      if (present(loaded)) then
         if (loaded) pgf_l = "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"
      end if
      want_psurf = present(psurf)
      psurf_s = "0.0"
      if (want_psurf) write (psurf_s, '(ES16.8)') psurf

      write (depth_s, '(F12.2)') max_depth
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 12, ny = 10, nghost = 2, dx = 4000.0, dy = 4000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 6 /"//new_line("a")// &
            "&time_nml t_end = 100000.0, dt_fixed = 300.0 /"//new_line("a")// &
            ! No wind: under a full-domain lid the P2c cover mask zeroes
            ! it in the cavity run and not in the twin, which would make
            ! the two runs different experiments.  The stirrer is the
            ! `U_SEED` initial jet `run_case` writes into both — see
            ! that parameter's docstring.
            "&ocean_topo_nml max_depth = "//trim(adjustl(depth_s))//" /"//new_line("a")// &
            "&physics_nml coriolis_f = 1.0e-4 /"//new_line("a")// &
            ! Stratified: T from 4 degC at the bed to 12 degC at the
            ! surface, so every layer carries a different density and the
            ! baroclinic pressure gradient is live.
            "&tracer_nml initial_salinity = 35.0, T_init_bottom = 4.0, "// &
            "T_init_surface = 12.0 /"//new_line("a")// &
            pgf_l//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .false., "// &
            "n_inner = 12 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
      if (want_psurf) then
         nml = nml//"&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
               "&ocean_psurf_nml enable = .true., p_surf_const = "// &
               trim(adjustl(psurf_s))//" /"//new_line("a")
      end if
      if (len_trim(cavity) > 0) then
         nml = nml//"&ocean_cavity_dyn_nml "//cavity//" /"//new_line("a")
      end if
   end function case_nml

   subroutine run_case(max_depth, cavity, n_steps, h, u, v, eta, ok, loaded, psurf)
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
      logical, intent(in), optional :: loaded
      real(wp), intent(in), optional :: psurf

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen
      real(wp), pointer :: f3(:, :, :), f2(:, :)
      character(len=:), allocatable :: nml

      ok = .false.
      handle = c_null_ptr
      if (present(psurf)) then
         nml = case_nml(max_depth, cavity, loaded=loaded, psurf=psurf)
      else
         nml = case_nml(max_depth, cavity, loaded=loaded)
      end if
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return

      ! The stirrer: an identical initial zonal jet in BOTH runs, written
      ! through the API (which narrow-pushes it to the device, so this is
      ! correct on the GPU build too).  See `U_SEED`.
      status = seed_initial_jet(handle)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle)
         return
      end if

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

   function seed_initial_jet(handle) result(status)
      !! Write the `U_SEED` initial zonal jet into `u_face_x_layer` over
      !! the physical interior, depth-uniform and meridionally shaped
      !! like the two-gyre stress this suite used to blow:
      !!
      !!     u(i, j, k) = U_SEED * (1 - cos(2*pi*(j_phys - 0.5)/ny_phys))
      !!
      !! The profile is a pure function of the PHYSICAL index, so the
      !! cavity run and its shallow twin (same `nx`/`ny`/`nz`, different
      !! bed) receive bit-identical arrays — which is what
      !! `cavity_flat_lid_seed_bit_identical` then asserts.
      !!
      !! `rdb_ocean_set_u` narrow-pushes with `!$acc update device`, so
      !! the seed reaches the device under `mem:separate`; it does NOT
      !! re-derive `hu_face_x_layer`, which the first step recomputes.
      type(c_ptr), intent(in) :: handle
      integer(c_int) :: status
      integer(c_int) :: nx_p, ny_p, nz_p, ng
      real(wp), allocatable :: ubuf(:, :, :)
      real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)
      integer :: i, j, k

      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      if (status /= OCEAN_STATUS_OK) return

      allocate (ubuf(nx_p + 1, ny_p, nz_p))
      do k = 1, nz_p
         do j = 1, ny_p
            do i = 1, nx_p + 1
               ubuf(i, j, k) = U_SEED* &
                               (1.0_wp - cos(TWO_PI*(real(j, wp) - 0.5_wp)/real(ny_p, wp)))
            end do
         end do
      end do
      status = rdb_ocean_set_u(handle, ubuf, nx_p, ny_p, nz_p)
   end function seed_initial_jet

   pure function bound_accel_from(dp_scale) result(tol_u)
      !! Velocity-difference bound from a PRESSURE-difference scale.
      !!
      !!   * the Pass-3 face force divides a `pa*h` difference by
      !!     `(h_L+h_R)`, so a `dp_scale` (Pa) discrepancy between the two
      !!     runs' pressure stacks becomes an acceleration discrepancy
      !!     `~ dp_scale/(rho_0*dx)`;
      !!   * over `N_STEPS` steps of `DT` that integrates to
      !!     `dp_scale*DT*N_STEPS/(rho_0*dx)`.
      !!
      !! The factor 64 covers the handful of roundings per pass (three
      !! passes, the barotropic substeps, the ALE remap) and the mild
      !! step-to-step amplification.  It is the SAME factor in both
      !! variants, so the two bounds differ only by their physics.
      real(wp), intent(in) :: dp_scale
      real(wp) :: tol_u
      tol_u = 64.0_wp*dp_scale*DT*real(N_STEPS, wp)/(RHO_REF*DX)
   end function bound_accel_from

   pure function bound_accel() result(tol_u)
      !! UNLOADED variant.  The two `pa` stacks are offset by the constant
      !! `C = rho_ref*g*DRAFT ~ 5e6 Pa`, so every `pa` entry in the cavity
      !! run carries an absolute rounding of order `C*eps` that the
      !! shallow run's `O(1e4 Pa)` stack does not.
      real(wp) :: tol_u
      tol_u = bound_accel_from((RHO_REF*GRAVITY*DRAFT)*epsilon(1.0_wp))
   end function bound_accel

   pure function bound_accel_loaded() result(tol_u)
      !! LOADED variant.  `p_top = p_ice_ref` cancels `C` inside
      !! `pa(nz+1) = rho_ref*g*eta_geo + p_top`, so the `C*eps` term the
      !! unloaded bound is made of is simply GONE: `pa(nz+1)` is zero to
      !! the rounding of one product, and the stack that marches down from
      !! it is nothing but the density anomaly it was always meant to be.
      !!
      !! So the scale that replaces `C` is the size of that anomaly stack:
      !!
      !!     |pa| <= |rho_layer - rho_ref| * g * (water column)
      !!           = (ALPHA_T * T_HALF_RANGE) * GRAVITY * (BED - DRAFT)
      !!
      !! with `ALPHA_T` the `&ocean_ic_nml` linear-EOS coefficient this
      !! namelist leaves at its default and `T_HALF_RANGE` the largest
      !! `|T - T_ref|` its 4..12 degC profile reaches about `T_ref = 10`.
      !! Salinity is uniform at `S_ref`, so it contributes nothing.
      !!
      !! `PA_ANOM` here is ~5 Pa against the unloaded `C` of 5.08e6 Pa, so
      !! this bound is SIX DECADES tighter, and it is what the measured
      !! numbers in the module header sit under.
      real(wp) :: tol_u
      tol_u = bound_accel_from(epsilon(1.0_wp)*PA_ANOM)
   end function bound_accel_loaded

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
                 "the seeded jet must keep the flow decades above the bound, else "// &
                 "the equivalence gate is vacuous")
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

   subroutine test_flat_lid_equivalence_loaded(error)
      !! P5.2: the SAME pair with the load wired through
      !! `&ocean_pgf_nml p_top_in_bc`, so `pa(nz+1) = rho_ref*g*eta_geo +
      !! p_ice_ref` cancels the `C = rho_ref*g*DRAFT ~ 5.08e6 Pa` offset
      !! at the top of the stack instead of carrying it through every
      !! entry.  The `dp_scale` the bound is built from loses that term
      !! and keeps only the irreducible one (the two runs stack their
      !! interfaces from different beds) — see `bound_accel_loaded`.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_a(:, :, :), u_a(:, :, :), v_a(:, :, :), eta_a(:, :)
      real(wp), allocatable :: h_b(:, :, :), u_b(:, :, :), v_b(:, :, :), eta_b(:, :)
      logical :: ok_a, ok_b
      real(wp) :: du, dv, dh, deta, tol_u, tol_h, signal

      call run_case(BED, "enable = .true., draft_config = 'flat', draft_depth = "// &
                    "500.0", N_STEPS, h_a, u_a, v_a, eta_a, ok_a, loaded=.true.)
      call check(error, ok_a, "the LOADED cavity run must step")
      if (allocated(error)) return
      call run_case(BED - DRAFT, "", N_STEPS, h_b, u_b, v_b, eta_b, ok_b, loaded=.true.)
      call check(error, ok_b, "the shallow twin must step")
      if (allocated(error)) return

      call check(error, all(ieee_is_finite(u_a)) .and. all(ieee_is_finite(h_a)), &
                 "the loaded cavity run must stay finite")
      if (allocated(error)) return

      du = maxval(abs(u_a - u_b))
      dv = maxval(abs(v_a - v_b))
      dh = maxval(abs(h_a - h_b))
      deta = maxval(abs(eta_a - eta_b))
      tol_u = bound_accel_loaded()
      tol_h = tol_u*DT*real(N_STEPS, wp)*BED/DX

      signal = maxval(abs(u_a))
      call check(error, signal > 1.0e4_wp*tol_u, &
                 "the seeded jet must keep the flow decades above the bound, else "// &
                 "the loaded equivalence gate is vacuous")
      if (allocated(error)) return

      call check(error, du <= tol_u, "u_face_x_layer must match the shallow twin to "// &
                 "the TIGHTER load-cancelled round-off bound")
      if (allocated(error)) return
      call check(error, dv <= tol_u, "v_face_y_layer must match the shallow twin")
      if (allocated(error)) return
      call check(error, dh <= tol_h, "h_layer must match the shallow twin")
      if (allocated(error)) return
      call check(error, deta <= tol_h, "bt_eta must match the shallow twin")
      if (allocated(error)) return

      ! The point of the slice: cancelling the load at the top of the
      ! stack must make the pair AGREE BETTER, not merely still agree.
      call check(error, tol_u < bound_accel(), &
                 "the loaded bound must be tighter than the unloaded one — else "// &
                 "the load bought nothing")
   end subroutine test_flat_lid_equivalence_loaded

   subroutine test_seam_zero_p_surf_inert(error)
      !! SEAM PROOF, part 1 — a cavity with `p_surf = 0` sends the
      !! `eta_forcing` seam NOTHING.
      !!
      !! Run the same loaded cavity twice: once with no `&ocean_psurf_nml`
      !! group at all (no seam array, `eta_forcing` absent from the
      !! barotropic substep's argument list) and once with the seam LIVE
      !! at `p_surf_const = 0`.  If the static ice load had been assembled
      !! into `sf%p_surf` — the array `eta_ib = -p_surf/(rho0 g_bt)` is
      !! built from — the second run would carry a 493 m seam elevation
      !! and could not possibly reproduce the first bit-for-bit.
      !!
      !! It does, because `p_ice_ref` goes to `ms%p_top` and to the datum
      !! and nowhere else: `eta_ib = -0/(rho0 g_bt)` is `-0.0`, and
      !! `eta - (-0.0)` is `eta` exactly under IEEE-754.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_a(:, :, :), u_a(:, :, :), v_a(:, :, :), eta_a(:, :)
      real(wp), allocatable :: h_b(:, :, :), u_b(:, :, :), v_b(:, :, :), eta_b(:, :)
      logical :: ok_a, ok_b
      character(len=*), parameter :: LID = &
                                     "enable = .true., draft_config = 'flat', draft_depth = 500.0"

      call run_case(BED, LID, N_STEPS, h_a, u_a, v_a, eta_a, ok_a, loaded=.true.)
      call check(error, ok_a, "the seam-free cavity run must step")
      if (allocated(error)) return
      call run_case(BED, LID, N_STEPS, h_b, u_b, v_b, eta_b, ok_b, &
                    loaded=.true., psurf=0.0_wp)
      call check(error, ok_b, "the cavity run with a live zero seam must step")
      if (allocated(error)) return

      call check(error, maxval(abs(u_a)) > 1.0e-6_wp, &
                 "the wind-driven flow must be non-trivial, else the seam gate is "// &
                 "vacuous")
      if (allocated(error)) return
      call check(error, all(u_a == u_b) .and. all(v_a == v_b), &
                 "turning the eta_forcing seam ON at p_surf = 0 must be BIT-identical "// &
                 "under a cavity — the static ice load never reaches sf%p_surf")
      if (allocated(error)) return
      call check(error, all(h_a == h_b) .and. all(eta_a == eta_b), &
                 "h_layer and bt_eta must be BIT-identical too")
   end subroutine test_seam_zero_p_surf_inert

   subroutine test_seam_matches_cavity_free(error)
      !! SEAM PROOF, part 2 — with `p_surf /= 0` the seam the cavity run
      !! feeds the barotropic substep is the one a cavity-free run feeds
      !! it, so the loaded pair still tracks.
      !!
      !! The bound here is the OFFSET-scale `bound_accel`, not the tighter
      !! `bound_accel_loaded`, and that is a derived statement about this
      !! configuration rather than a concession: with an atmospheric load
      !! present the cavity run assembles
      !! `p_top = p_ice_ref + sf%p_surf = 5.0777e6 + 2.0e3`, and the small
      !! addend inherits the large one's ulp — one rounding at
      !! `eps*rho_ref*g*DRAFT`, exactly the `dp_scale` the unloaded bound
      !! is built from.  The cavity-free twin adds `0 + 2000` exactly and
      !! pays nothing.  Measured: `du = 1.836e-14 m/s` against the
      !! `3.136e-11` bound — 1700x of margin, and 5 decades under the
      !! run's own `3.5e-3 m/s`.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: h_a(:, :, :), u_a(:, :, :), v_a(:, :, :), eta_a(:, :)
      real(wp), allocatable :: h_b(:, :, :), u_b(:, :, :), v_b(:, :, :), eta_b(:, :)
      logical :: ok_a, ok_b
      real(wp) :: du, tol_u

      call run_case(BED, "enable = .true., draft_config = 'flat', draft_depth = "// &
                    "500.0", N_STEPS, h_a, u_a, v_a, eta_a, ok_a, &
                    loaded=.true., psurf=P_SURF_SEAM)
      call check(error, ok_a, "the loaded cavity run with a live seam must step")
      if (allocated(error)) return
      call run_case(BED - DRAFT, "", N_STEPS, h_b, u_b, v_b, eta_b, ok_b, &
                    loaded=.true., psurf=P_SURF_SEAM)
      call check(error, ok_b, "the shallow twin with the same seam must step")
      if (allocated(error)) return

      du = maxval(abs(u_a - u_b))
      tol_u = bound_accel()
      call check(error, maxval(abs(u_a)) > 1.0e4_wp*tol_u, &
                 "the flow must be decades above the bound, else the seam gate is "// &
                 "vacuous")
      if (allocated(error)) return
      call check(error, du <= tol_u, &
                 "with the SAME p_surf the cavity run and the cavity-free twin must "// &
                 "still agree to the load-cancelled bound — i.e. the cavity put "// &
                 "nothing of its own on the eta_ib seam")
   end subroutine test_seam_matches_cavity_free

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
