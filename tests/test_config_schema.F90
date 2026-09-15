!! Unit tests for rdb_config_schema: the central Roundabout namelist schema.
module test_config_schema
   use rdb_constants, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_config, only: config_t
   use rdb_config_schema, only: build_rdb_schema
   use rdb_nml_schema, only: nml_schema_t
   implicit none
   private

   integer, parameter :: rk = wp
      !! Compare parsed config values (which are `real(wp)`) against
      !! `_rk` literals of the same kind, so the tolerance checks hold in
      !! a single-precision build (`0.3_wp` matches the parsed `0.3`).

   public :: collect_config_schema_tests

contains

   subroutine collect_config_schema_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("happy_path_kshear_epbl", test_happy), &
                  new_unittest("typo_key_suggestion", test_typo), &
                  new_unittest("out_of_range", test_range), &
                  new_unittest("doc_short_non_default", test_doc_short), &
                  new_unittest("legacy_group_external", test_external), &
                  new_unittest("ocean_bc_migrated", test_ocean_bc), &
                  new_unittest("retired_split_scheme_names_successor", &
                               test_retired_mom6_pc), &
                  new_unittest("retired_split_rk2_names_successor", &
                               test_retired_split_rk2), &
                  new_unittest("split_scheme_default_is_pred_corr", &
                               test_split_scheme_default) &
                  ]
   end subroutine collect_config_schema_tests

   subroutine write_file(path, body)
      character(len=*), intent(in) :: path, body
      integer :: u, i, start
      open (newunit=u, file=path, status='replace', action='write')
      start = 1
      do i = 1, len(body)
         if (body(i:i) == achar(10)) then
            write (u, '(A)') body(start:i - 1)
            start = i + 1
         end if
      end do
      if (start <= len(body)) write (u, '(A)') body(start:)
      close (u)
   end subroutine write_file

   subroutine rm_file(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open (newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close (u, status='delete')
   end subroutine rm_file

   logical function any_contains(arr, n, needle)
      character(len=*), intent(in) :: arr(:)
      integer, intent(in) :: n
      character(len=*), intent(in) :: needle
      integer :: i
      any_contains = .false.
      do i = 1, n
         if (index(arr(i), needle) > 0) then
            any_contains = .true.
            return
         end if
      end do
   end function any_contains

   logical function file_contains(path, needle)
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

   subroutine test_happy(error)
      !! Default config + a file setting kshear enable/ri_crit + epbl
      !! mstar_scheme="rh18": fields set, zero errors, enum canonicalized.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_happy.nml"

      call build_rdb_schema(cfg, schema)
      call write_file(fn, &
                      "&ocean_kappa_shear_nml"//achar(10)// &
                      "  enable = .true."//achar(10)// &
                      "  ri_crit = 0.3"//achar(10)// &
                      "/"//achar(10)// &
                      "&ocean_epbl_nml  mstar_scheme = 'RH18' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat == 0, "expected no errors")
      if (allocated(error)) return
      call check(error, cfg%ocean%kshear%enable, "kshear enable set")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%kshear%ri_crit - 0.3_rk) < 1.0e-12_rk, "ri_crit set")
      if (allocated(error)) return
      ! enum canonical spelling is the lowercase "rh18" from the allowed set.
      call check(error, trim(cfg%ocean%epbl%mstar_scheme) == "rh18", "enum canonical")
   end subroutine test_happy

   subroutine test_typo(error)
      !! Typo'd kshear key -> error with a did-you-mean suggestion.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_typo.nml"

      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_kappa_shear_nml  ri_crt = 0.3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown key"), "unknown-key msg")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "ri_crit"), "suggestion present")
   end subroutine test_typo

   subroutine test_range(error)
      !! Out-of-range kshear integer (max_inner_it < 1) -> error.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_range.nml"

      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_kappa_shear_nml  max_inner_it = 0 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "below min"), "range violation")
   end subroutine test_range

   subroutine test_retired_mom6_pc(error)
      !! `&ocean_bt_nml split_scheme = "mom6_pc"` -- the spelling this
      !! scheme shipped under until 2026-09-14 -- must fail LOUD and name
      !! its successor.  The two ways this could go wrong are both silent
      !! and both worse than an error: falling through to the default
      !! (the run would be right by accident today and wrong the next
      !! time a default moves), or the generic "not in allowed set {...}"
      !! message, which tells the reader the value is wrong without
      !! telling them what it is now called.
      type(error_type), allocatable, intent(out) :: error

      call check_retired_split_scheme(error, "mom6_pc", "scratch_cs_retired_mp.nml")
   end subroutine test_retired_mom6_pc

   subroutine test_retired_split_rk2(error)
      !! `split_rk2` was the SECOND spelling of the same scheme, held
      !! briefly between the `mom6_pc` rename and the `pred_corr` one.
      !! It collided with `dynamics/split_rk2/`, the directory holding
      !! the outer-loop machinery of BOTH schemes, so the name described
      !! the family as well as one member.  Both dead spellings are
      !! retired, and both must name `pred_corr` rather than fall through
      !! to the default.
      type(error_type), allocatable, intent(out) :: error

      call check_retired_split_scheme(error, "split_rk2", "scratch_cs_retired_sr.nml")
   end subroutine test_retired_split_rk2

   subroutine check_retired_split_scheme(error, spelling, fn)
      !! Shared body: `spelling` must be rejected with the migration
      !! message naming `pred_corr`, never the generic allowed-set one,
      !! and must not reach the config.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: spelling, fn
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_bt_nml  split_scheme = '"//spelling//"' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "retired split_scheme value must be an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "RETIRED"), &
                 "error must say the value is retired")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "pred_corr"), &
                 "error must name the new spelling")
      if (allocated(error)) return
      call check(error,.not. any_contains(errs, stat, "not in allowed set"), &
                 "must not ALSO emit the generic allowed-set message")
      if (allocated(error)) return
      ! And it must NOT have been accepted into the config.
      call check(error, trim(cfg%ocean%bt%split_scheme) /= spelling, &
                 "retired value must not reach the config")
   end subroutine check_retired_split_scheme

   subroutine test_split_scheme_default(error)
      !! The registered default is the MOM6 predictor-corrector.  Guards
      !! the half of the pair that lives in `rdb_config.F90`; the other
      !! half (`ocean_dyn_t%split_scheme`) is guarded in
      !! `test_ocean_pred_corr`.  They have diverged before.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema

      call build_rdb_schema(cfg, schema)
      call check(error, trim(cfg%ocean%bt%split_scheme) == "pred_corr", &
                 "&ocean_bt_nml split_scheme default must be 'pred_corr'")
   end subroutine test_split_scheme_default

   subroutine test_doc_short(error)
      !! write_doc_short with ONE non-default knob contains exactly it.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_ds.nml"
      character(len=*), parameter :: sf = "scratch_cs_short.nml"

      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_kappa_shear_nml  shearmix_rate = 0.5 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "parse ok")
      if (allocated(error)) return

      call schema%write_doc_short(sf)
      call check(error, file_contains(sf, "shearmix_rate"), "short has the non-default knob")
      if (allocated(error)) return
      ! ri_crit is still default -> must be omitted from the short doc.
      call check(error,.not. file_contains(sf, "ri_crit"), "short omits default knob")
      call rm_file(sf)
   end subroutine test_doc_short

   subroutine test_external(error)
      !! Every group is now migrated onto the strict schema — there are
      !! no external (skipped) groups left (P4.5 closed the last one,
      !! ocean_bc).  An ocean group that used to be external
      !! (ocean_topo_nml) is now strictly validated: a known key applies,
      !! an unknown key errors, and a bad enum value errors.  A migrated
      !! coastal group (physics_nml) likewise rejects garbage.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat, n_ext
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_ext.nml"

      call build_rdb_schema(cfg, schema)

      ! No external (skipped) groups remain.
      n_ext = 0
      if (allocated(schema%external_names)) n_ext = size(schema%external_names)
      call check(error, n_ext == 0, "no transitional external groups remain")
      if (allocated(error)) return

      ! A clean ocean_topo_nml validates and applies.
      call write_file(fn, "&ocean_topo_nml  max_depth = 4000.0  topo_config = 'spoon' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "clean ocean_topo_nml validates")
      if (allocated(error)) return

      ! ocean_topo_nml is now strict — an unknown key must error.
      call write_file(fn, "&ocean_topo_nml  max_depth = 4000.0  some_garbage = 3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "migrated ocean_topo_nml rejects unknown key")
      if (allocated(error)) return

      ! A bad enum value (topo_config) must error too.
      call write_file(fn, "&ocean_topo_nml  topo_config = 'nonsense' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "ocean_topo_nml rejects bad topo_config enum")
      if (allocated(error)) return

      ! physics_nml is a strict group — an unknown key must error.
      call write_file(fn, "&physics_nml  manning_n = 0.025  some_garbage = 3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "migrated physics_nml rejects unknown key")
   end subroutine test_external

   subroutine test_ocean_bc(error)
      !! P4.5: `&ocean_bc_nml` is now a strictly-validated schema group
      !! (no more hand-rolled `read_ocean_bc_nml`).  Covers: defaults
      !! reproduce a closed-wall run, a clean block applies + canonicalises
      !! its enums, an unknown key errors, a bad edge-type enum errors
      !! naming the allowed set, and a tidal-constituent array knob applies.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_bc.nml"

      ! Pristine defaults reproduce a closed-wall run.
      call build_rdb_schema(cfg, schema)
      call check(error, trim(cfg%ocean%bc%west) == "wall", "default west = wall")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%bc%north) == "wall", "default north = wall")
      if (allocated(error)) return

      ! A clean block applies and canonicalises.
      call write_file(fn, "&ocean_bc_nml  west = 'PERIODIC'  east = 'periodic' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "clean ocean_bc_nml validates")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%bc%west) == "periodic", "west canonicalised to lowercase")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%bc%east) == "periodic", "east applied")
      if (allocated(error)) return

      ! An unknown key must error (was silently ignored by the native
      ! namelist reader before P4.5).
      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_bc_nml  west = 'wall'  some_garbage = 3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat >= 1, "migrated ocean_bc_nml rejects unknown key")
      if (allocated(error)) return

      ! A bad edge-type enum must error naming the allowed set (was a
      ! silent fall-back to 'wall' before P4.5 — the NVHPC in-memory bug
      ! this migration fixes).
      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_bc_nml  west = 'not-a-real-bc' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat >= 1, "bad west enum value errors")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "allowed set"), "error names the allowed set")
      if (allocated(error)) return

      ! Tidal-constituent array knob (fixed-size real_array).
      call build_rdb_schema(cfg, schema)
      call write_file(fn, "&ocean_bc_nml  west_n_tidal = 1  west_tidal_amp = 1.5 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "tidal array block validates")
      if (allocated(error)) return
      call check(error, cfg%ocean%bc%west_n_tidal == 1, "west_n_tidal applied")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%bc%west_tidal_amp(1) - 1.5_rk) < 1.0e-12_rk, &
                 "west_tidal_amp(1) applied")
   end subroutine test_ocean_bc

end module test_config_schema
