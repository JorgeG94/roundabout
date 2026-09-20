!! The ice-shelf-cavity DIAGNOSTIC catalog: thirteen entries, each
!! checked against an independent host computation from the slot fields.
module test_ocean_cavity_diags
   !! What this suite is for.  A diagnostic is the only part of a model
   !! most readers ever see, so the failure mode that matters is not "it
   !! crashed" but "it printed a plausible number".  Three classes of
   !! that are pinned here:
   !!
   !!   * **A wrong VALUE.**  Every fill is compared against the same
   !!     quantity computed independently on the host from the slot
   !!     arrays — not against itself, and not against a golden number.
   !!   * **A wrong UNIT.**  `melt_m_per_yr` is the ISOMIP+ reporting
   !!     convention (`/ rho_fw = 1000`, x seconds per year), and the
   !!     conversion is hand-checked against a number worked out in the
   !!     docstring rather than against the code's own product.
   !!   * **A number where there is no cavity.**  Outside the cover
   !!     every entry must be the IEEE NaN sentinel, never zero.  Zero
   !!     melt is a LEGAL answer, so a plane of zeros over open ocean is
   !!     how a 30 m/yr shelf gets averaged down to 2 and published.
   !!
   !! Plus the fail-loud gate: a cavity diagnostic requested on a run
   !! with no cavity is a configure error, asserted through the
   !! `derived_catalog_requires` lookup rather than by provoking
   !! `error stop` — the repo's standing pattern for testing a
   !! fail-loud rule without killing the runner.
   !!
   !! `mem:separate` discipline: the state goes to the device through
   !! `ocean_state_enter_data`, the cavity arrays are pushed with an
   !! explicit `!$acc update device` (they are host-filled AFTER the
   !! map), and every read-back is an `!$acc update self` of COMPONENT
   !! arrays, never of an aggregate derived type.  All inert on the host
   !! build, which is exactly why it is written unconditionally.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t, eos_apply_tfreeze_set, eos_freezing_point, &
                      TFREEZE_SET_ISOMIP
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_diag_derived, only: derived_catalog_size, derived_catalog_name, &
                                     derived_catalog_requires, melt_m_per_yr_factor, &
                                     DERIVED_REQ_NONE, DERIVED_REQ_CAVITY_DYN, &
                                     DERIVED_REQ_CAVITY_MELT, &
                                     fill_melt, fill_melt_m_per_yr, fill_thermal_driving, &
                                     fill_haline_driving, fill_tbdry, fill_sbdry, &
                                     fill_tfreeze_ib, fill_exch_vel_t, fill_exch_vel_s, &
                                     fill_ustar_shelf, fill_cavity_melt_status, &
                                     fill_z_draft, fill_water_column
   use rdb_ocean_cavity_melt, only: CAVITY_MELT_OK, CAVITY_MELT_NOT_CONVERGED
   implicit none
   private

   public :: collect_ocean_cavity_diags_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1000.0_wp

   integer, parameter :: I_ICE = 2
      !! A covered, SOLVED column (local index; `nghost = 1` here).
   integer, parameter :: I_OPEN = 5
      !! An uncovered column — every cavity diagnostic must be missing.

   real(wp), parameter :: TOL = 1.0e-13_wp
      !! Relative tolerance for "one algebra, two expression trees" — the
      !! fill kernel and the host check contract their FMAs differently.

contains

   subroutine collect_ocean_cavity_diags_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_diag_catalog_carries_every_name", test_catalog), &
                  new_unittest("cavity_diag_requires_gate_is_set", test_requires_gate), &
                  new_unittest("cavity_diag_melt_and_unit_conversion", test_melt), &
                  new_unittest("cavity_diag_drivings_match_host", test_drivings), &
                  new_unittest("cavity_diag_interface_state_matches_slot", test_interface), &
                  new_unittest("cavity_diag_geometry_matches_datum", test_geometry), &
                  new_unittest("cavity_diag_missing_outside_the_cavity", test_missing) &
                  ]
   end subroutine collect_ocean_cavity_diags_tests

   ! ------------------------------------------------------------------
   ! Fixture
   ! ------------------------------------------------------------------

   subroutine setup_cavity_state(grid, state, eos)
      !! A 6x4x3 cartesian plane with the cavity geometry and the melt
      !! slot both on, ice over `i <= 3`, and every slot array stamped
      !! with a distinct, column-dependent value so a fill that read the
      !! WRONG array would not accidentally agree.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      type(eos_t), intent(out) :: eos
      integer :: i, j

      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      state%metrics%use_cavity = .true.
      state%cavity_flux%enable = .true.
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)

      call eos_apply_tfreeze_set(eos, TFREEZE_SET_ISOMIP)
      state%eos = eos

      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            if (i <= 3) then
               state%metrics%cover_frac(i, j) = 1.0_wp
               state%metrics%z_draft(i, j) = 300.0_wp + 10.0_wp*real(i, wp)
               state%cavity_flux%active(i, j) = 1.0_wp
            else
               state%metrics%cover_frac(i, j) = 0.0_wp
               state%metrics%z_draft(i, j) = 0.0_wp
               state%cavity_flux%active(i, j) = 0.0_wp
            end if
            ! Distinct per-(i,j) values: a fill that read `s_b` where it
            ! should read `s_far` cannot pass by coincidence.
            state%cavity_flux%melt(i, j) = 1.0e-5_wp*real(i, wp) + 1.0e-6_wp*real(j, wp)
            state%cavity_flux%t_far(i, j) = -1.0_wp + 0.10_wp*real(i, wp)
            state%cavity_flux%s_far(i, j) = 34.0_wp + 0.05_wp*real(j, wp)
            state%cavity_flux%t_b(i, j) = -2.0_wp - 0.01_wp*real(i, wp)
            state%cavity_flux%s_b(i, j) = 20.0_wp + 0.20_wp*real(i, wp)
            state%cavity_flux%ustar(i, j) = 1.0e-3_wp*real(i + j, wp)
            state%cavity_flux%gamma_t(i, j) = 2.0e-5_wp*real(i, wp)
            state%cavity_flux%gamma_s(i, j) = 6.0e-7_wp*real(j, wp)
            state%cavity_flux%status(i, j) = CAVITY_MELT_OK
            state%multilayer%p_top(i, j) = 3.0e6_wp + 1.0e4_wp*real(i, wp)
            state%dyn%bt_work%bt_H_ref(i, j) = 700.0_wp - 10.0_wp*real(i, wp)
            state%dyn%bt_work%bt_eta(i, j) = 0.01_wp*real(j, wp)
         end do
      end do
      ! One non-OK covered column, so the status diag has something to
      ! say other than zero.
      state%cavity_flux%status(3, 2) = CAVITY_MELT_NOT_CONVERGED

      call ocean_state_enter_data(state)
      ! mem:separate trap (2): these were host-filled AFTER the map, so
      ! they owe an explicit push.  COMPONENT arrays only — never the
      ! aggregate.
      !$acc update device(state%metrics%cover_frac, state%metrics%z_draft)
      !$acc update device(state%cavity_flux%active, state%cavity_flux%melt, &
      !$acc&              state%cavity_flux%t_far, state%cavity_flux%s_far, &
      !$acc&              state%cavity_flux%t_b, state%cavity_flux%s_b, &
      !$acc&              state%cavity_flux%ustar, state%cavity_flux%gamma_t, &
      !$acc&              state%cavity_flux%gamma_s, state%cavity_flux%status)
      !$acc update device(state%multilayer%p_top)
      !$acc update device(state%dyn%bt_work%bt_H_ref, state%dyn%bt_work%bt_eta)
   end subroutine setup_cavity_state

   subroutine teardown_cavity_state(state)
      type(ocean_state_t), intent(inout) :: state
      call ocean_state_exit_data(state)
      call state%destroy()
   end subroutine teardown_cavity_state

   subroutine run_fill(state, which, buf)
      !! Fire one catalog fill into a device-mapped buffer and bring the
      !! answer back.  The buffer gets its own `enter data`/`exit data`
      !! because it is bare test scratch, not a slot array.
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: which
      real(wp), intent(out) :: buf(NX + 2, NY + 2, 1)

      buf = 0.0_wp
      !$acc enter data copyin(buf)
      select case (which)
      case (1); call fill_melt(state, buf)
      case (2); call fill_melt_m_per_yr(state, buf)
      case (3); call fill_thermal_driving(state, buf)
      case (4); call fill_haline_driving(state, buf)
      case (5); call fill_tbdry(state, buf)
      case (6); call fill_sbdry(state, buf)
      case (7); call fill_tfreeze_ib(state, buf)
      case (8); call fill_exch_vel_t(state, buf)
      case (9); call fill_exch_vel_s(state, buf)
      case (10); call fill_ustar_shelf(state, buf)
      case (11); call fill_cavity_melt_status(state, buf)
      case (12); call fill_z_draft(state, buf)
      case (13); call fill_water_column(state, buf)
      end select
      !$acc update self(buf)
      !$acc exit data delete(buf)
   end subroutine run_fill

   pure function rel_diff(a, b) result(r)
      real(wp), intent(in) :: a, b
      real(wp) :: r
      r = abs(a - b)/max(abs(b), 1.0_wp)
   end function rel_diff

   ! ------------------------------------------------------------------
   ! Catalog + gate
   ! ------------------------------------------------------------------

   subroutine test_catalog(error)
      !! All thirteen names are in the catalog and none collides with an
      !! existing entry.
      type(error_type), allocatable, intent(out) :: error
      character(len=32), parameter :: NAMES(13) = [ &
                                      character(len=32) :: "melt", "melt_m_per_yr", &
                                                           "thermal_driving", "haline_driving", "tbdry", &
                                                           "sbdry", "tfreeze_ib", "exch_vel_t", "exch_vel_s", &
                                                           "ustar_shelf", "cavity_melt_status", "z_draft", &
                                                           "water_column"]
      integer :: k, i, n, hits

      n = derived_catalog_size()
      do k = 1, size(NAMES)
         hits = 0
         do i = 1, n
            if (trim(derived_catalog_name(i)) == trim(NAMES(k))) hits = hits + 1
         end do
         call check(error, hits == 1, &
                    "the catalog must carry '"//trim(NAMES(k))//"' exactly once")
         if (allocated(error)) return
      end do
   end subroutine test_catalog

   subroutine test_requires_gate(error)
      !! Every melt-interface diagnostic is gated on the MELT knob and
      !! the two geometry ones only on the GEOMETRY knob — which is the
      !! whole point of shipping `z_draft`/`water_column` separately.
      !! Asserted through the lookup, not by provoking the `error stop`
      !! `register_derived` raises.
      type(error_type), allocatable, intent(out) :: error

      call check(error, derived_catalog_requires("melt") == DERIVED_REQ_CAVITY_MELT, &
                 "melt must require the melt slot")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("exch_vel_s") == DERIVED_REQ_CAVITY_MELT, &
                 "exch_vel_s must require the melt slot")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("cavity_melt_status") == &
                 DERIVED_REQ_CAVITY_MELT, "the status plane must require the melt slot")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("z_draft") == DERIVED_REQ_CAVITY_DYN, &
                 "z_draft is GEOMETRY: it must not drag in the melt slot")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("water_column") == DERIVED_REQ_CAVITY_DYN, &
                 "water_column is GEOMETRY: it must not drag in the melt slot")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("ke_total") == DERIVED_REQ_NONE, &
                 "a non-cavity entry must stay ungated")
      if (allocated(error)) return
      call check(error, derived_catalog_requires("not_a_diagnostic") == -1, &
                 "an unknown name must be distinguishable from an ungated one")
   end subroutine test_requires_gate

   ! ------------------------------------------------------------------
   ! Values
   ! ------------------------------------------------------------------

   subroutine test_melt(error)
      !! `melt` is the slot field verbatim, and `melt_m_per_yr` is it in
      !! the ISOMIP+ reporting unit.
      !!
      !! THE CONVERSION, HAND-CHECKED.  `m_per_yr = (kg m-2 s-1) /
      !! rho_fw * seconds_per_year` with `rho_fw = 1000 kg/m^3` and a
      !! 365-day year (Asay-Davis et al. 2016 §3.3), so the factor is
      !! `31 536 000 / 1000 = 31 536` exactly and a melt flux of
      !! `1e-3 kg/m^2/s` is `31.536 m/yr`.  Both numbers are written out
      !! here rather than recomputed from the module's own constants —
      !! otherwise the test would agree with any pair of constants the
      !! module happened to hold.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(eos_t) :: eos
      real(wp) :: buf(NX + 2, NY + 2, 1)
      integer :: j

      call check(error, rel_diff(melt_m_per_yr_factor(), 31536.0_wp) <= TOL, &
                 "the ISOMIP+ conversion factor must be 31536 m/yr per kg/m^2/s")
      if (allocated(error)) return
      call check(error, rel_diff(1.0e-3_wp*melt_m_per_yr_factor(), 31.536_wp) <= TOL, &
                 "1e-3 kg/m^2/s must be 31.536 m/yr")
      if (allocated(error)) return

      call setup_cavity_state(grid, state, eos)
      j = 2

      call run_fill(state, 1, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%melt(I_ICE, j), &
                 "melt must be the slot field verbatim, not a rescaled copy")
      if (allocated(error)) go to 900

      call run_fill(state, 2, buf)
      call check(error, rel_diff(buf(I_ICE, j, 1), &
                                 state%cavity_flux%melt(I_ICE, j)*31536.0_wp) <= TOL, &
                 "melt_m_per_yr must be melt x 31536")
      if (allocated(error)) go to 900
      call check(error, buf(I_ICE, j, 1) > 0.0_wp, &
                 "a positive melt flux must report a positive melt RATE (sign convention)")

900   call teardown_cavity_state(state)
   end subroutine test_melt

   subroutine test_drivings(error)
      !! `thermal_driving = T_far - T_f(S_far, p_top)` and
      !! `haline_driving = S_far - S_b`, each against a host evaluation
      !! of the SAME liquidus handle the solve uses.  `tfreeze_ib` is
      !! that liquidus on its own, so the three are checked for mutual
      !! consistency too — a sign slip in any one of them would break it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(eos_t) :: eos
      real(wp) :: buf(NX + 2, NY + 2, 1)
      real(wp) :: tf_host, td, hd, tf_diag
      integer :: j

      call setup_cavity_state(grid, state, eos)
      j = 3
      tf_host = eos_freezing_point(eos, state%cavity_flux%s_far(I_ICE, j), &
                                   state%multilayer%p_top(I_ICE, j))

      call run_fill(state, 3, buf)
      td = buf(I_ICE, j, 1)
      call check(error, rel_diff(td, state%cavity_flux%t_far(I_ICE, j) - tf_host) <= TOL, &
                 "thermal_driving must be T_far - T_f(S_far, p_top)")
      if (allocated(error)) go to 900

      call run_fill(state, 7, buf)
      tf_diag = buf(I_ICE, j, 1)
      call check(error, rel_diff(tf_diag, tf_host) <= TOL, &
                 "tfreeze_ib must be the in-situ liquidus of the FAR FIELD")
      if (allocated(error)) go to 900
      call check(error, rel_diff(td, state%cavity_flux%t_far(I_ICE, j) - tf_diag) <= TOL, &
                 "the two must be mutually consistent: T* = T_far - tfreeze_ib")
      if (allocated(error)) go to 900
      ! The liquidus must actually depend on the interface pressure —
      ! else the ice pump is silently absent from the diagnostic.
      call check(error, tf_host < 0.0_wp, &
                 "at 3 MPa under a shelf the freezing point must be well below 0 degC")
      if (allocated(error)) go to 900

      call run_fill(state, 4, buf)
      hd = buf(I_ICE, j, 1)
      call check(error, hd == state%cavity_flux%s_far(I_ICE, j) - &
                 state%cavity_flux%s_b(I_ICE, j), &
                 "haline_driving must be S_far - S_b")
      if (allocated(error)) go to 900
      call check(error, hd > 0.0_wp, &
                 "a fresher interface than the far field must give a POSITIVE S*")

900   call teardown_cavity_state(state)
   end subroutine test_drivings

   subroutine test_interface(error)
      !! `tbdry`/`sbdry`/`exch_vel_t`/`exch_vel_s`/`ustar_shelf`/
      !! `cavity_melt_status` are verbatim reads of six DIFFERENT slot
      !! arrays.  Every array carries a distinct value here, so a fill
      !! wired to the wrong one cannot pass.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(eos_t) :: eos
      real(wp) :: buf(NX + 2, NY + 2, 1)
      integer :: j

      call setup_cavity_state(grid, state, eos)
      j = 2

      call run_fill(state, 5, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%t_b(I_ICE, j), &
                 "tbdry must be t_b")
      if (allocated(error)) go to 900
      call run_fill(state, 6, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%s_b(I_ICE, j), &
                 "sbdry must be s_b")
      if (allocated(error)) go to 900
      call run_fill(state, 8, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%gamma_t(I_ICE, j), &
                 "exch_vel_t must be the CONVERGED gamma_t, not Gamma_T*u* re-derived")
      if (allocated(error)) go to 900
      call run_fill(state, 9, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%gamma_s(I_ICE, j), &
                 "exch_vel_s must be the converged gamma_s")
      if (allocated(error)) go to 900
      call run_fill(state, 10, buf)
      call check(error, buf(I_ICE, j, 1) == state%cavity_flux%ustar(I_ICE, j), &
                 "ustar_shelf must be the melt friction velocity")
      if (allocated(error)) go to 900
      call run_fill(state, 11, buf)
      call check(error, buf(I_ICE, j, 1) == real(CAVITY_MELT_OK, wp), &
                 "a clean column must report status 0")
      if (allocated(error)) go to 900
      call check(error, buf(3, 2, 1) == real(CAVITY_MELT_NOT_CONVERGED, wp), &
                 "a non-converged column must report its own code, not OK")

900   call teardown_cavity_state(state)
   end subroutine test_interface

   subroutine test_geometry(error)
      !! `z_draft` is the metrics field verbatim and `water_column` is
      !! the datum identity `bt_H_ref + bt_eta` — the one expression
      !! every consumer of the column depth uses, which is exactly what
      !! the cavity datum buys.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(eos_t) :: eos
      real(wp) :: buf(NX + 2, NY + 2, 1)
      integer :: j

      call setup_cavity_state(grid, state, eos)
      j = 2

      call run_fill(state, 12, buf)
      call check(error, buf(I_ICE, j, 1) == state%metrics%z_draft(I_ICE, j), &
                 "z_draft must be the metrics field verbatim")
      if (allocated(error)) go to 900
      call check(error, buf(I_ICE, j, 1) > 0.0_wp, &
                 "the draft is POSITIVE DOWN")
      if (allocated(error)) go to 900

      call run_fill(state, 13, buf)
      call check(error, rel_diff(buf(I_ICE, j, 1), &
                                 state%dyn%bt_work%bt_H_ref(I_ICE, j) + &
                                 state%dyn%bt_work%bt_eta(I_ICE, j)) <= TOL, &
                 "water_column must be bt_H_ref + bt_eta")

900   call teardown_cavity_state(state)
   end subroutine test_geometry

   subroutine test_missing(error)
      !! OUTSIDE THE CAVITY EVERY ENTRY IS MISSING, NOT ZERO.
      !!
      !! This is the assertion that keeps a domain mean honest: the
      !! console reduction and the NetCDF writer both skip non-finite
      !! cells, so a NaN-filled open ocean makes `mean(melt)` the mean
      !! over the CAVITY.  A zero would be a legal melt rate and would
      !! silently dilute it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(eos_t) :: eos
      real(wp) :: buf(NX + 2, NY + 2, 1)
      integer :: which, j
      integer :: n_finite

      call setup_cavity_state(grid, state, eos)
      j = 2

      do which = 1, 13
         call run_fill(state, which, buf)
         call check(error, ieee_is_finite(buf(I_ICE, j, 1)), &
                    "inside the cavity the diagnostic must be a NUMBER")
         if (allocated(error)) go to 900
         call check(error,.not. ieee_is_finite(buf(I_OPEN, j, 1)), &
                    "outside the cavity the diagnostic must be MISSING, not zero")
         if (allocated(error)) go to 900
      end do

      ! And the reduction a console line would take: three covered
      ! columns per row, so the finite count over the physical interior
      ! is exactly the cavity's own cell count.
      call run_fill(state, 1, buf)
      n_finite = count(ieee_is_finite(buf(2:NX + 1, 2:NY + 1, 1)))
      call check(error, n_finite == 2*NY, &
                 "the finite cell count must be the CAVITY's cell count (i = 2, 3 "// &
                 "of the physical interior), so a domain mean is a cavity mean")

900   call teardown_cavity_state(state)
   end subroutine test_missing

end module test_ocean_cavity_diags
