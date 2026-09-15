!! PR-8 dead-knob amnesty: retirement tests.
module test_ocean_knob_retirement
   !! Roundabout's config schema is strict — `rdb_config.F90` `error stop`s on
   !! any unknown key.  PR-8 deleted five registered `&ocean_*` keys that
   !! validated and then did nothing (`&ocean_bt_nml correction_visc_rem`
   !! is the one exception — kept, per its own coordination note in the
   !! plan, and NOT covered here).  This suite pins the contract that
   !! makes the retirement real: setting a retired key in a `.nml` must be
   !! REJECTED with an "unknown key" error, not silently ignored — the
   !! exact defect class this PR closes.  Cloned from
   !! `tests/test_config_schema.F90:test_typo`
   !! (itself modelled on `tests/test_nml_schema.F90:test_unknown_key`).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_config, only: config_t
   use rdb_config_schema, only: build_rdb_schema
   use rdb_nml_schema, only: nml_schema_t
   implicit none
   private

   public :: collect_ocean_knob_retirement_tests

contains

   subroutine collect_ocean_knob_retirement_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("redi_ref_pres_retired", test_redi_ref_pres), &
                  new_unittest("tides_use_equilibrium_retired", test_tides_use_equilibrium), &
                  new_unittest("tides_use_eq_phase_retired", test_tides_use_eq_phase), &
                  new_unittest("foxkemper_apply_cfl_limit_retired", test_foxkemper_apply_cfl_limit) &
                  ]
   end subroutine collect_ocean_knob_retirement_tests

   subroutine write_file(path, body)
      character(len=*), intent(in) :: path, body
      integer :: u
      open (newunit=u, file=path, status='replace', action='write')
      write (u, '(A)') body
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

   subroutine assert_retired(error, fn, group_line, key_line)
      !! Shared body: write a scratch `.nml` with ONE retired key set
      !! inside its real group, parse against the real production
      !! schema, and assert the strict engine rejects it as unknown.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: fn, group_line, key_line
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)

      call build_rdb_schema(cfg, schema)
      call write_file(fn, trim(group_line)//achar(10)//trim(key_line)//achar(10)//"/")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "retired key should be rejected, not ignored")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown key"), "unknown-key msg")
   end subroutine assert_retired

   subroutine test_redi_ref_pres(error)
      !! `&ocean_redi_nml ref_pres` — accepted-then-unread (§2.5 item 2):
      !! copied into `ocean_redi_t%ref_pres` but never referenced in
      !! either `redi_calc_coeffs_x`/`_y` body.  Deleted whole (field,
      !! knob, copy site, kernel dummy).
      type(error_type), allocatable, intent(out) :: error
      call assert_retired(error, "scratch_knobret_ref_pres.nml", &
                          "&ocean_redi_nml", "  ref_pres = 2.0e7")
   end subroutine test_redi_ref_pres

   subroutine test_tides_use_equilibrium(error)
      !! `&ocean_tides_nml use_equilibrium` — accepted, copied to the
      !! slot, never read; `enable=.true.` always applies the full body
      !! tide regardless.
      type(error_type), allocatable, intent(out) :: error
      call assert_retired(error, "scratch_knobret_use_equilibrium.nml", &
                          "&ocean_tides_nml", "  use_equilibrium = .false.")
   end subroutine test_tides_use_equilibrium

   subroutine test_tides_use_eq_phase(error)
      !! `&ocean_tides_nml use_eq_phase` — accepted, never even copied to
      !! a slot; `rdb_ocean_tides.F90` unconditionally bakes V_c into
      !! phase0 regardless of the setting.
      type(error_type), allocatable, intent(out) :: error
      call assert_retired(error, "scratch_knobret_use_eq_phase.nml", &
                          "&ocean_tides_nml", "  use_eq_phase = .false.")
   end subroutine test_tides_use_eq_phase

   subroutine test_foxkemper_apply_cfl_limit(error)
      !! `&ocean_foxkemper_nml apply_cfl_limit` — accepted then hard
      !! `error stop`s at configure for asking for a thing the code
      !! already always does (the availability cap in
      !! `mle_compute_transports` is unconditional in production).
      !! Deleted as incoherent, not merely deferred.
      type(error_type), allocatable, intent(out) :: error
      call assert_retired(error, "scratch_knobret_apply_cfl_limit.nml", &
                          "&ocean_foxkemper_nml", "  apply_cfl_limit = .true.")
   end subroutine test_foxkemper_apply_cfl_limit

end module test_ocean_knob_retirement
