!! Schema-engine behaviour, driven through the TOML front end.
module test_config_schema_toml
   !! The strict schema (`rdb_nml_schema`) is the one validation engine behind the
   !! TOML reader.  These tests build a small purpose-made schema (not the Roundabout
   !! one) and drive it with TOML text, pinning the engine rules a TOML user meets:
   !! the range minimum AND maximum, required keys, default tracking after a parse,
   !! the per-group `cross_check` callback, the doc writers, a two-element real array,
   !! retired enum spellings and the did-you-mean hints.  Table-syntax refusals and
   !! the group <-> table mapping are pinned in `test_toml_config`.
   use rdb_constants, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_nml_schema, only: nml_schema_t, nml_group_t, &
                             nml_real, nml_int, nml_logical, nml_string, &
                             nml_enum, nml_real_array
   use rdb_toml_read, only: toml_schema_parse_string
   implicit none
   private

   integer, parameter :: rk = wp
      !! The schema's working real kind, so the pointer targets below associate
      !! in either a single- or double-precision build.

   public :: collect_config_schema_toml_tests

   character(len=*), parameter :: NL = achar(10)

   ! Config slots the schema points at; each test resets them through build_schema.
   real(rk), target :: cfg_ri_crit
   integer, target :: cfg_n_iter
   logical, target :: cfg_use_lt
   character(len=32), target :: cfg_form
   character(len=16), target :: cfg_variant
   real(rk), target :: cfg_layers(4)
   real(rk), target :: cfg_lo, cfg_hi
      !! Cross-check slots.

contains

   subroutine collect_config_schema_toml_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("schema_toml_happy_path_all_types", test_happy_path), &
                  new_unittest("schema_toml_unknown_key_and_table_hints", test_unknown), &
                  new_unittest("schema_toml_range_min_max_enum_strlen", test_range_enum_strlen), &
                  new_unittest("schema_toml_required_key", test_required), &
                  new_unittest("schema_toml_defaults_tracking", test_defaults), &
                  new_unittest("schema_toml_doc_writers", test_doc_writers), &
                  new_unittest("schema_toml_cross_check_propagation", test_cross_check), &
                  new_unittest("schema_toml_retired_enum_spellings", test_retired_enum), &
                  new_unittest("schema_toml_duplicate_table_refused", test_duplicates) &
                  ]
   end subroutine collect_config_schema_toml_tests

   ! ---- helpers --------------------------------------------------------

   subroutine reset_targets()
      !! Restore the config slots to the defaults the schema captures at registration.
      cfg_ri_crit = 0.25_rk
      cfg_n_iter = 10
      cfg_use_lt = .false.
      cfg_form = "sadourny"
      cfg_variant = "om4"
      cfg_layers = [1000.0_rk, 1000.0_rk, 0.0_rk, 0.0_rk]
   end subroutine reset_targets

   function str_ptr_form() result(p)
      !! Deferred-length pointer associated with the fixed-length form slot.
      character(len=:), pointer :: p
      p => cfg_form
   end function str_ptr_form

   function str_ptr_variant() result(p)
      !! Deferred-length pointer associated with the fixed-length variant slot.
      character(len=:), pointer :: p
      p => cfg_variant
   end function str_ptr_variant

   subroutine build_schema(schema)
      !! A representative schema over the module slots: `[ocean.kappa_shear]` (every
      !! key type) and `[ocean.epbl]` (an enum with two retired spellings).
      type(nml_schema_t), intent(out) :: schema
      type(nml_group_t) :: g, g2
      real(rk), pointer :: pr, pa(:)
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps, ps2

      call reset_targets()

      g%name = "ocean_kappa_shear"
      g%doc = "Shear-driven mixing knobs."
      pr => cfg_ri_crit
      call g%add(nml_real("ri_crit", pr, "Critical Richardson number", &
                          units="nondim", min=0.0_rk, max=2.0_rk))
      pi => cfg_n_iter
      call g%add(nml_int("n_iter", pi, "Solver iterations", min=1, max=100))
      pl => cfg_use_lt
      call g%add(nml_logical("use_lt", pl, "Enable Langmuir"))
      ps => str_ptr_form()
      call g%add(nml_string("form", ps, "Coriolis form"))
      pa => cfg_layers
      call g%add(nml_real_array("layers", pa, "Layer thicknesses", units="m"))
      call schema%add_group(g)

      g2%name = "ocean_epbl"
      g2%doc = "Energetics PBL knobs."
      ps2 => str_ptr_variant()
      call g2%add(nml_enum("mstar_variant", ps2, "mstar scheme", &
                           allowed=[character(len=8) :: "constant", "om4", "rh18"], &
                           retired=[character(len=8) :: "om3", "om3b"], &
                           retired_hint="renamed to 'om4'"))
      call schema%add_group(g2)
   end subroutine build_schema

   logical function any_contains(arr, n, needle)
      !! True if any of arr(1:n) contains the substring `needle`.
      character(len=*), intent(in) :: arr(:)
      integer, intent(in) :: n
      character(len=*), intent(in) :: needle
      integer :: i
      any_contains = .false.
      do i = 1, min(n, size(arr))
         if (index(arr(i), needle) > 0) then
            any_contains = .true.
            return
         end if
      end do
   end function any_contains

   logical function file_contains(path, needle)
      !! True if any line of `path` contains `needle`.
      character(len=*), intent(in) :: path, needle
      integer :: u, ios
      character(len=512) :: buf
      file_contains = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(A)', iostat=ios) buf
         if (ios /= 0) exit
         if (index(buf, needle) > 0) then
            file_contains = .true.
            exit
         end if
      end do
      close (u)
   end function file_contains

   subroutine rm_file(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open (newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close (u, status='delete')
   end subroutine rm_file

   ! ---- tests ----------------------------------------------------------

   subroutine test_happy_path(error)
      !! Every key type parses; the enum canonicalises; a two-element real array fills
      !! the first two slots of the target and leaves the rest alone.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, &
                                    "[ocean.kappa_shear]"//NL// &
                                    "ri_crit = 0.5  # inline comment"//NL// &
                                    "n_iter = 20"//NL// &
                                    "use_lt = true"//NL// &
                                    "form = 'leith'"//NL// &
                                    "layers = [500.0, 250.0]"//NL// &
                                    "[ocean.epbl]"//NL//"mstar_variant = 'RH18'"//NL, &
                                    status=stat, errors=errs)

      call check(error, stat == 0, "expected no errors")
      if (allocated(error)) return
      call check(error, abs(cfg_ri_crit - 0.5_rk) < 1.0e-12_rk, "ri_crit")
      if (allocated(error)) return
      call check(error, cfg_n_iter == 20, "n_iter")
      if (allocated(error)) return
      call check(error, cfg_use_lt, "use_lt")
      if (allocated(error)) return
      call check(error, trim(cfg_form) == "leith", "form")
      if (allocated(error)) return
      call check(error, abs(cfg_layers(1) - 500.0_rk) < 1.0e-9_rk, "layer 1")
      if (allocated(error)) return
      call check(error, abs(cfg_layers(2) - 250.0_rk) < 1.0e-9_rk, "layer 2")
      if (allocated(error)) return
      call check(error, abs(cfg_layers(3)) < 1.0e-12_rk .and. abs(cfg_layers(4)) < 1.0e-12_rk, &
                 "slots past the supplied elements are untouched")
      if (allocated(error)) return
      call check(error, trim(cfg_variant) == "rh18", "enum canonical spelling from the allowed list")
   end subroutine test_happy_path

   subroutine test_unknown(error)
      !! An unknown key and an unknown table each error with a did-you-mean hint.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shear]"//NL//"ri_crt = 0.3"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown key"), "unknown-key message")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "ri_crit"), "suggestion present")
      if (allocated(error)) return

      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shr]"//NL//"ri_crit = 0.3"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat >= 1, "expected an error for the table")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown table"), "unknown-table message")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "ocean.kappa_shear"), "table suggestion")
   end subroutine test_unknown

   subroutine test_range_enum_strlen(error)
      !! Below the minimum, above the maximum, an over-long string and an enum value
      !! outside the allowed set all error (one run, every error collected).
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, &
                                    "[ocean.kappa_shear]"//NL// &
                                    "ri_crit = -1.0"//NL// &
                                    "n_iter = 999"//NL// &
                                    "form = '"//repeat('x', 40)//"'"//NL// &  ! longer than cfg_form (32)
                                    "[ocean.epbl]"//NL//"mstar_variant = 'bogus'"//NL, &
                                    status=stat, errors=errs)
      call check(error, any_contains(errs, stat, "below min"), "min violation")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "above max"), "max violation (integer key)")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "exceeds target length"), "string length")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "not in allowed set"), "enum")
      if (allocated(error)) return

      ! The max side of a REAL key, on its own.
      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shear]"//NL//"ri_crit = 2.5"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat == 1 .and. any_contains(errs, stat, "above max"), "max violation (real key)")
      if (allocated(error)) return
      ! The bounds themselves are accepted.
      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shear]"//NL//"ri_crit = 2.0"//NL//"n_iter = 1"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat == 0, "values on the bounds are accepted")
   end subroutine test_range_enum_strlen

   subroutine test_retired_enum(error)
      !! EVERY spelling listed in `retired=` is rejected with the migration hint, not
      !! with the generic allowed-set message and not silently.  The list carries two
      !! entries on purpose: a key can outlive more than one rename, so the match must
      !! walk the whole list.
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      character(len=4), parameter :: dead(2) = [character(len=4) :: "om3", "om3b"]

      do i = 1, size(dead)
         call check_retired_spelling(error, trim(dead(i)))
         if (allocated(error)) return
      end do
   end subroutine test_retired_enum

   subroutine check_retired_spelling(error, spelling)
      !! Shared body of `test_retired_enum`, one retired spelling.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: spelling
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.epbl]"//NL//"mstar_variant = '"//spelling//"'"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat >= 1, "retired value must error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "RETIRED"), "says retired")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "renamed to"), "quotes the hint")
      if (allocated(error)) return
      call check(error,.not. any_contains(errs, stat, "not in allowed set"), &
                 "must not ALSO emit the generic allowed-set message")
      if (allocated(error)) return
      call check(error, trim(cfg_variant) == "om4", "target left at its default")
   end subroutine check_retired_spelling

   subroutine test_duplicates(error)
      !! A table given twice is refused, as is a key given twice.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, &
                                    "[ocean.kappa_shear]"//NL//"ri_crit = 0.3"//NL// &
                                    "[ocean.kappa_shear]"//NL//"n_iter = 5"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat >= 1, "a repeated table must be refused")
      if (allocated(error)) return

      call build_schema(schema)
      call toml_schema_parse_string(schema, &
                                    "[ocean.kappa_shear]"//NL//"ri_crit = 0.3"//NL//"ri_crit = 0.4"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat >= 1, "a repeated key must be refused")
   end subroutine test_duplicates

   subroutine test_required(error)
      !! A required key errors when its table is present but the key absent; an absent
      !! table is fine (defaults rule).
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      type(nml_group_t) :: g
      integer :: stat
      character(len=:), allocatable :: errs(:)
      real(rk), pointer :: pr

      call reset_targets()
      g%name = "reqgrp"
      pr => cfg_ri_crit
      call g%add(nml_real("must_set", pr, "Mandatory knob", required=.true.))
      call schema%add_group(g)

      ! Table present, key absent -> error naming the key.
      call toml_schema_parse_string(schema, "[reqgrp]"//NL, status=stat, errors=errs)
      call check(error, any_contains(errs, stat, "required key"), "required key fires")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "must_set"), "message names the key")
      if (allocated(error)) return

      ! Table present with the key -> accepted.
      call reset_targets()
      call toml_schema_parse_string(schema, "[reqgrp]"//NL//"must_set = 1.5"//NL, status=stat, errors=errs)
      call check(error, stat == 0, "required key present is accepted")
      if (allocated(error)) return

      ! Table absent entirely -> OK.
      call toml_schema_parse_string(schema, "# nothing here"//NL, status=stat, errors=errs)
      call check(error, stat == 0, "absent table is accepted")
   end subroutine test_required

   subroutine test_defaults(error)
      !! Untouched keys keep their defaults and report `is_default`; after setting one
      !! knob only that key reports non-default.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shear]"//NL//"n_iter = 42"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat == 0, "no errors")
      if (allocated(error)) return
      call check(error, cfg_n_iter == 42, "n_iter set")
      if (allocated(error)) return
      call check(error, abs(cfg_ri_crit - 0.25_rk) < 1.0e-12_rk, "ri_crit default kept")
      if (allocated(error)) return
      call check(error, schema%groups(1)%keys(1)%key%is_default(), "ri_crit is_default")
      if (allocated(error)) return
      call check(error,.not. schema%groups(1)%keys(2)%key%is_default(), "n_iter non-default")
      if (allocated(error)) return
      call check(error, schema%groups(2)%keys(1)%key%is_default(), "other group untouched")
   end subroutine test_defaults

   subroutine test_doc_writers(error)
      !! `write_doc_all`, `write_doc_short` and `render_markdown` after a parse.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: af = "scratch_schema_all.txt"
      character(len=*), parameter :: sf = "scratch_schema_short.txt"
      character(len=*), parameter :: mf = "scratch_schema_md.md"

      call build_schema(schema)
      call toml_schema_parse_string(schema, "[ocean.kappa_shear]"//NL//"n_iter = 77"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat == 0, "parse ok")
      if (allocated(error)) return

      call schema%write_doc_all(af)
      call schema%write_doc_short(sf)
      call schema%render_markdown(mf)

      ! doc_all lists every key; an untouched one is marked (default).
      call check(error, file_contains(af, "ri_crit"), "all has ri_crit")
      if (allocated(error)) return
      call check(error, file_contains(af, "(default)"), "all marks defaults")
      if (allocated(error)) return
      ! doc_short lists only the non-default key.
      call check(error, file_contains(sf, "n_iter"), "short has n_iter")
      if (allocated(error)) return
      call check(error,.not. file_contains(sf, "ri_crit"), "short omits defaults")
      if (allocated(error)) return
      ! The dumps are TOML: table headers, bare TOML values, `#` comments.
      call check(error, file_contains(af, "[ocean.kappa_shear]") .and. file_contains(af, "[ocean.epbl]"), &
                 "all carries every [ocean.x] table header")
      if (allocated(error)) return
      call check(error, file_contains(af, "ri_crit = 0.25  # [nondim]"), "all: real as a short TOML float")
      if (allocated(error)) return
      call check(error, file_contains(af, "n_iter = 77  # "), "all: integer")
      if (allocated(error)) return
      call check(error, file_contains(af, "use_lt = false  # "), "all: logical true/false")
      if (allocated(error)) return
      call check(error, file_contains(af, 'form = "sadourny"  # '), "all: quoted string")
      if (allocated(error)) return
      call check(error, file_contains(af, "layers = [1000.0, 1000.0, 0.0, 0.0]  # "), "all: array")
      if (allocated(error)) return
      call check(error, file_contains(sf, "[ocean.kappa_shear]") .and. file_contains(sf, "n_iter = 77"), &
                 "short carries the table header of a group with a non-default key")
      if (allocated(error)) return
      call check(error,.not. file_contains(sf, "[ocean.epbl]"), "short omits an all-default group")
      if (allocated(error)) return
      ! The markdown reference carries its table header and every knob.
      call check(error, file_contains(mf, "| Knob |"), "markdown table header")
      if (allocated(error)) return
      call check(error, file_contains(mf, "mstar_variant"), "markdown lists every group")
      if (allocated(error)) return
      call check(error, file_contains(mf, "| `ri_crit` | `0.25` |"), "markdown shows the same value text")

      call rm_file(af)
      call rm_file(sf)
      call rm_file(mf)
   end subroutine test_doc_writers

   function lo_hi_cross_check() result(msg)
      !! Cross-check: `lo` must not exceed `hi`.
      character(len=:), allocatable :: msg
      msg = ''
      if (cfg_lo > cfg_hi) msg = "lo exceeds hi"
   end function lo_hi_cross_check

   subroutine test_cross_check(error)
      !! A registered `cross_check` error is collected after the keys apply, and a
      !! consistent pair passes.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      type(nml_group_t) :: g
      integer :: stat
      character(len=:), allocatable :: errs(:)
      real(rk), pointer :: plo, phi

      cfg_lo = 0.0_rk
      cfg_hi = 10.0_rk
      g%name = "xcgrp"
      plo => cfg_lo
      phi => cfg_hi
      call g%add(nml_real("lo", plo, "low"))
      call g%add(nml_real("hi", phi, "high"))
      g%cross_check => lo_hi_cross_check
      call schema%add_group(g)

      call toml_schema_parse_string(schema, "[xcgrp]"//NL//"lo = 5.0"//NL//"hi = 1.0"//NL, &
                                    status=stat, errors=errs)
      call check(error, any_contains(errs, stat, "lo exceeds hi"), "cross_check message")
      if (allocated(error)) return

      call toml_schema_parse_string(schema, "[xcgrp]"//NL//"lo = 1.0"//NL//"hi = 5.0"//NL, &
                                    status=stat, errors=errs)
      call check(error, stat == 0, "a consistent pair passes the cross-check")
   end subroutine test_cross_check

end module test_config_schema_toml
