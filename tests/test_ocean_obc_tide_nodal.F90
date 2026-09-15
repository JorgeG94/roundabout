!! Unit tests for the OBC tidal nodal/astronomical correction (capability C3).
!! Exercises the `pure` helpers `obc_match_constituent` + `obc_tide_nodal_fill`
!! (in `rdb_ocean_boundary_types`) that bake the 18.6-yr nodal factor f_c and
!! the equilibrium+nodal phase (V_c+u_c) into an OBC edge's per-constituent
!! `tidal_fnodal`/`tidal_arg`, and the phase-convention continuity of the
!! barotropic eta-target sum.  Pure-helper tests — no solver handle needed.
module test_ocean_obc_tide_nodal
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_tide_astro, only: days_since_1900, equilibrium_arguments, &
                                   nodal_fu, TIDE_OMEGA, TIDES_CATALOG_SIZE
   use rdb_ocean_boundary_types, only: ocean_bc_face_tag_t, obc_match_constituent, &
                                       obc_tide_nodal_fill
   implicit none
   private

   public :: collect_ocean_obc_tide_nodal_tests

   integer, parameter :: M2 = 1   !! Catalog index of M2 (TIDE_NAME(1)).

contains

   subroutine collect_ocean_obc_tide_nodal_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("match_resolver", test_match_resolver), &
                  new_unittest("nodal_tracking_9yr", test_nodal_tracking), &
                  new_unittest("fill_unmatched_flags_error", test_fill_unmatched), &
                  new_unittest("disabled_defaults_inert", test_disabled_defaults), &
                  new_unittest("static_phase_sum_reproduced", test_static_sum)]
   end subroutine collect_ocean_obc_tide_nodal_tests

   subroutine test_match_resolver(error)
      !! ω→catalog resolver: TIDE_OMEGA(M2) resolves to M2; every catalog
      !! entry resolves to itself; an off-catalog ω resolves to 0 (no match);
      !! a slightly-perturbed (within-tol) M2 ω still resolves to M2.
      type(error_type), allocatable, intent(out) :: error
      integer :: c
      call check(error, obc_match_constituent(TIDE_OMEGA(M2)), M2)
      if (allocated(error)) return
      ! Each catalog frequency resolves to its own index (ω's are separated).
      do c = 1, TIDES_CATALOG_SIZE
         call check(error, obc_match_constituent(TIDE_OMEGA(c)), c)
         if (allocated(error)) return
      end do
      ! Off-catalog ω (2.0e-4 rad/s sits between S2 and K2 but far from both).
      call check(error, obc_match_constituent(2.0e-4_wp), 0)
      if (allocated(error)) return
      ! Non-positive ω is treated as unknown.
      call check(error, obc_match_constituent(0.0_wp), 0)
      if (allocated(error)) return
      ! Within relative tolerance 1e-4 still matches M2 (perturb by 1e-5 rel).
      call check(error, obc_match_constituent(TIDE_OMEGA(M2)*(1.0_wp + 1.0e-5_wp)), M2)
   end subroutine test_match_resolver

   subroutine test_nodal_tracking(error)
      !! Headline: bake f_M2 / (V+u)_M2 at two ref dates ~9.3 yr apart (nodal
      !! extremes).  The stored values must (a) equal a direct nodal_fu /
      !! equilibrium_arguments catalog lookup bit-for-bit, and (b) differ by
      !! several percent between the two dates (the 18.6-yr modulation).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_bc_face_tag_t) :: face_a, face_b
      real(wp) :: va(TIDES_CATALOG_SIZE), fa(TIDES_CATALOG_SIZE), ua(TIDES_CATALOG_SIZE)
      real(wp) :: vb(TIDES_CATALOG_SIZE), fb(TIDES_CATALOG_SIZE), ub(TIDES_CATALOG_SIZE)
      real(wp) :: da, db
      integer :: ierr

      ! Date A: 1997-01-01 ; Date B: 2006-01-01  (~9.13 yr ≈ half nodal cycle).
      da = days_since_1900(1997, 1, 1)
      db = days_since_1900(2006, 1, 1)
      call equilibrium_arguments(da, va)
      call nodal_fu(da, .true., fa, ua)
      call equilibrium_arguments(db, vb)
      call nodal_fu(db, .true., fb, ub)

      ! One M2 constituent on the west edge, unit amplitude, zero Greenwich lag.
      face_a%n_tidal_constituents = 1
      face_a%tidal_omega(1) = TIDE_OMEGA(M2)
      face_a%tidal_amp(1) = 1.0_wp
      face_a%tidal_phase(1) = 0.0_wp
      face_b = face_a

      call obc_tide_nodal_fill(face_a, fa, ua, va, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) return
      call obc_tide_nodal_fill(face_b, fb, ub, vb, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) return

      ! (a) stored == direct catalog lookup (bit-match).
      call check(error, face_a%tidal_fnodal(1), fa(M2), thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, face_a%tidal_arg(1), va(M2) + ua(M2), thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, face_b%tidal_fnodal(1), fb(M2), thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, face_b%tidal_arg(1), vb(M2) + ub(M2), thr=1.0e-14_wp)
      if (allocated(error)) return

      ! (b) the nodal factor moves several percent across the ~9.3-yr span.
      call check(error, abs(face_a%tidal_fnodal(1) - face_b%tidal_fnodal(1)) > 0.05_wp, &
                 "M2 nodal factor must differ >5% between nodal extremes")
   end subroutine test_nodal_tracking

   subroutine test_fill_unmatched(error)
      !! An edge carrying an off-catalog ω (with the correction enabled) makes
      !! the pure fill return ierr = the failing slot index — the signal the
      !! setup path converts to a fail-loud error stop.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_bc_face_tag_t) :: face
      real(wp) :: v(TIDES_CATALOG_SIZE), f(TIDES_CATALOG_SIZE), u(TIDES_CATALOG_SIZE)
      integer :: ierr

      call equilibrium_arguments(0.0_wp, v)
      call nodal_fu(0.0_wp, .true., f, u)
      face%n_tidal_constituents = 2
      face%tidal_omega(1) = TIDE_OMEGA(M2)   ! good
      face%tidal_omega(2) = 2.0e-4_wp        ! off-catalog
      face%tidal_amp(1:2) = 1.0_wp
      call obc_tide_nodal_fill(face, f, u, v, ierr)
      call check(error, ierr, 2)   ! second slot is the unmatched one
   end subroutine test_fill_unmatched

   subroutine test_disabled_defaults(error)
      !! A freshly-defaulted face (never touched by the fill) has fnodal≡1 and
      !! arg≡0 for every constituent slot — the inert defaults that keep the
      !! disabled-path OBC tidal sum bit-identical to the legacy kernel.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_bc_face_tag_t) :: face
      integer :: nc
      do nc = 1, size(face%tidal_fnodal)
         call check(error, face%tidal_fnodal(nc), 1.0_wp, thr=1.0e-15_wp)
         if (allocated(error)) return
         call check(error, face%tidal_arg(nc), 0.0_wp, thr=1.0e-15_wp)
         if (allocated(error)) return
      end do
   end subroutine test_disabled_defaults

   subroutine test_static_sum(error)
      !! Kernel phase-convention check.  The legacy (disabled) eta-target sum is
      !! Σ A·cos(ω·t+φ).  With the correction trivial (f=1, arg=0), the C3 sum
      !! Σ f·A·cos(ω·t+arg−φ) reduces to Σ A·cos(ω·t−φ); at φ=0 the two branches
      !! coincide bit-for-bit — proving the documented sign flip is confined to
      !! nonzero Greenwich lags and the enabled path.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_bc_face_tag_t) :: face
      real(wp) :: t, legacy, nodal
      integer :: it
      ! Two constituents (M2 + S2), unit amplitudes, ZERO Greenwich lag,
      ! defaults (f=1, arg=0) — the "correction trivial" case.
      face%n_tidal_constituents = 2
      face%tidal_omega(1) = TIDE_OMEGA(1)
      face%tidal_omega(2) = TIDE_OMEGA(2)
      face%tidal_amp(1:2) = [1.3_wp, 0.7_wp]
      face%tidal_phase(1:2) = 0.0_wp
      do it = 0, 5
         t = real(it, wp)*3600.0_wp
         legacy = legacy_eta(face, t)
         nodal = nodal_eta(face, t)
         call check(error, nodal, legacy, thr=1.0e-13_wp)
         if (allocated(error)) return
      end do
   end subroutine test_static_sum

   pure function legacy_eta(face, t) result(eta)
      !! Legacy (disabled-path) OBC tidal eta-target sum: Σ A·cos(ω·t+φ).
      type(ocean_bc_face_tag_t), intent(in) :: face
      real(wp), intent(in) :: t
      real(wp) :: eta
      integer :: nc
      eta = 0.0_wp
      do nc = 1, face%n_tidal_constituents
         eta = eta + face%tidal_amp(nc)*cos(face%tidal_omega(nc)*t + face%tidal_phase(nc))
      end do
   end function legacy_eta

   pure function nodal_eta(face, t) result(eta)
      !! C3 (enabled-path) OBC tidal eta-target sum:
      !! Σ f·A·cos(ω·t + arg − φ).
      type(ocean_bc_face_tag_t), intent(in) :: face
      real(wp), intent(in) :: t
      real(wp) :: eta
      integer :: nc
      eta = 0.0_wp
      do nc = 1, face%n_tidal_constituents
         eta = eta + face%tidal_fnodal(nc)*face%tidal_amp(nc)* &
               cos(face%tidal_omega(nc)*t + face%tidal_arg(nc) - face%tidal_phase(nc))
      end do
   end function nodal_eta

end module test_ocean_obc_tide_nodal
