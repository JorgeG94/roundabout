!! Unit tests for the unified `parse_vcoord_type` parser.
!!
!! Covers:
!!   1. Every recognised string maps to the right `VCOORD_*` code.
!!   2. Alias strings (`zsigma` / `z-sigma` / `z_sigma`, etc.) all map
!!      to the same code as the canonical form.
!!   3. Case-insensitivity for the canonical lowercase variants by
!!      also accepting UPPER (where the parser ships explicit
!!      alternatives).
!!   4. Unrecognised strings fall back to the user-supplied
!!      `default_code` (and to `VCOORD_SIGMA` when no default is
!!      passed — coastal-path historical default).
!!   5. `parse_ocean_vcoord_type` pins the fallback to
!!      `VCOORD_EULERIAN_Z` regardless of caller default.
!!
!! Note: `scoord` / `VCOORD_SCOORD` removed (D7); "scoord" now falls
!! back to the default (SIGMA for parse_vcoord_type).  If the S&H
!! stretch is later exposed as a `stretch_type="song_haidvogel"` option
!! on the existing sigma path, a new test can be added there.
module test_vcoord_parse
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, VCOORD_SIGMA, &
                            VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, VCOORD_Z_FIXED
   use rdb_vcoord, only: parse_vcoord_type
   use rdb_ocean_vcoord, only: parse_ocean_vcoord_type
   implicit none
   private

   public :: collect_vcoord_parse_tests

contains

   subroutine collect_vcoord_parse_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("canonical_strings", test_canonical), &
                  new_unittest("alias_strings", test_aliases), &
                  new_unittest("ocean_alias_gprime", test_ocean_gprime_alias), &
                  new_unittest("unknown_default_sigma", test_unknown_default_sigma), &
                  new_unittest("unknown_default_supplied", test_unknown_default_supplied), &
                  new_unittest("ocean_wrapper_pins_eulerian_z", test_ocean_wrapper_default), &
                  new_unittest("eulerian_z_recognised", test_eulerian_z), &
                  new_unittest("lagrangian_isopycnal", test_lagrangian) &
                  ]
   end subroutine collect_vcoord_parse_tests

   subroutine test_lagrangian(error)
      type(error_type), allocatable, intent(out) :: error
      ! "lagrangian" / "isopycnal" (+ UPPER) all map to VCOORD_LAGRANGIAN.
      call check(error, parse_vcoord_type("lagrangian") == VCOORD_LAGRANGIAN, &
                 "lagrangian not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("isopycnal") == VCOORD_LAGRANGIAN, &
                 "isopycnal not parsed to VCOORD_LAGRANGIAN")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("LAGRANGIAN") == VCOORD_LAGRANGIAN, &
                 "LAGRANGIAN (upper) not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("ISOPYCNAL") == VCOORD_LAGRANGIAN, &
                 "ISOPYCNAL (upper) not parsed")
      if (allocated(error)) return
      ! And via the ocean wrapper.
      call check(error, parse_ocean_vcoord_type("isopycnal") == VCOORD_LAGRANGIAN, &
                 "ocean wrapper failed on isopycnal")
   end subroutine test_lagrangian

   subroutine test_canonical(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_vcoord_type("sigma") == VCOORD_SIGMA, &
                 "sigma not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zsigma") == VCOORD_ZSIGMA, &
                 "zsigma not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar") == VCOORD_ZSTAR, &
                 "zstar not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar_full") == VCOORD_ZSTAR_FULL, &
                 "zstar_full not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar_sigma") == VCOORD_ZSTAR_SIGMA, &
                 "zstar_sigma not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z_fixed") == VCOORD_Z_FIXED, &
                 "z_fixed not parsed")
   end subroutine test_canonical

   subroutine test_aliases(error)
      type(error_type), allocatable, intent(out) :: error
      ! Hyphen + underscore variants should map identically.
      call check(error, parse_vcoord_type("z-sigma") == VCOORD_ZSIGMA, &
                 "z-sigma alias not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z_sigma") == VCOORD_ZSIGMA, &
                 "z_sigma alias not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z-star") == VCOORD_ZSTAR, &
                 "z-star alias not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar_lite") == VCOORD_ZSTAR, &
                 "zstar_lite alias not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z-star-full") == VCOORD_ZSTAR_FULL, &
                 "z-star-full alias not parsed")
      if (allocated(error)) return
      ! UPPERCASE recognised forms.
      call check(error, parse_vcoord_type("SIGMA") == VCOORD_SIGMA, &
                 "SIGMA (upper) not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("ZSTAR") == VCOORD_ZSTAR, &
                 "ZSTAR (upper) not parsed")
   end subroutine test_aliases

   subroutine test_ocean_gprime_alias(error)
      type(error_type), allocatable, intent(out) :: error
      ! MOM6 `COORD_CONFIG="gprime"` analogue should land on Z_FIXED.
      call check(error, parse_vcoord_type("gprime") == VCOORD_Z_FIXED, &
                 "gprime alias not parsed to z_fixed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z_levels") == VCOORD_Z_FIXED, &
                 "z_levels alias not parsed to z_fixed")
   end subroutine test_ocean_gprime_alias

   subroutine test_unknown_default_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      ! No default supplied ⇒ historical coastal fallback SIGMA.
      call check(error, parse_vcoord_type("not_a_real_coord") == VCOORD_SIGMA, &
                 "unknown string should fall back to VCOORD_SIGMA")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("") == VCOORD_SIGMA, &
                 "empty string should fall back to VCOORD_SIGMA")
   end subroutine test_unknown_default_sigma

   subroutine test_unknown_default_supplied(error)
      type(error_type), allocatable, intent(out) :: error
      ! Caller can override the fallback.
      call check(error, &
                 parse_vcoord_type("nonsense", default_code=VCOORD_ZSTAR) == VCOORD_ZSTAR, &
                 "default_code override ignored")
      if (allocated(error)) return
      call check(error, &
                 parse_vcoord_type("garbage", default_code=VCOORD_EULERIAN_Z) == VCOORD_EULERIAN_Z, &
                 "default_code = VCOORD_EULERIAN_Z not honored")
   end subroutine test_unknown_default_supplied

   subroutine test_ocean_wrapper_default(error)
      type(error_type), allocatable, intent(out) :: error
      ! The thin ocean wrapper must pin the fallback to EULERIAN_Z.
      call check(error, parse_ocean_vcoord_type("bogus") == VCOORD_EULERIAN_Z, &
                 "parse_ocean_vcoord_type fallback should be VCOORD_EULERIAN_Z")
      if (allocated(error)) return
      ! And still parse known cases correctly.
      call check(error, parse_ocean_vcoord_type("zstar_full") == VCOORD_ZSTAR_FULL, &
                 "ocean wrapper failed on zstar_full")
   end subroutine test_ocean_wrapper_default

   subroutine test_eulerian_z(error)
      type(error_type), allocatable, intent(out) :: error
      ! "eulerian_z" (and short aliases "z", "Z") map to the
      ! VCOORD_EULERIAN_Z constant (no longer ocean-only).
      call check(error, parse_vcoord_type("eulerian_z") == VCOORD_EULERIAN_Z, &
                 "eulerian_z not parsed")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z") == VCOORD_EULERIAN_Z, &
                 "z alias not parsed to VCOORD_EULERIAN_Z")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("Z") == VCOORD_EULERIAN_Z, &
                 "Z alias not parsed to VCOORD_EULERIAN_Z")
   end subroutine test_eulerian_z

end module test_vcoord_parse
