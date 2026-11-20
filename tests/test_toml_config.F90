!! Unit tests for the TOML configuration reader (rdb_toml_read).
module test_toml_config
   !! These tests pin (1) that a configuration spelled in TOML resolves to the
   !! values it states -- exactly, bit for bit for reals -- (2) every refusal
   !! class of the TOML syntax, with the file/line and key in the message, and
   !! (3) the group <-> table mapping.  The resolved dumps (`write_doc_all`,
   !! `write_resolved`) serve as the oracle where a whole configuration is
   !! compared.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: DEG2RAD
   use rdb_config, only: config_t, read_config, read_config_from_toml_string, build_rdb_schema
   use rdb_nml_schema, only: nml_schema_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_PARSE
   use rdb_toml_read, only: group_to_table, table_to_group, toml_schema_parse_string, &
                            toml_schema_parse_file
   implicit none
   private

   public :: collect_toml_config_tests

   character(len=*), parameter :: NL = achar(10)

contains

   subroutine collect_toml_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("all_key_types_resolve", test_all_key_types), &
                  new_unittest("double_gyre_resolves", test_double_gyre), &
                  new_unittest("file_entry_and_non_toml_refused", test_file_entry), &
                  new_unittest("oracle_sees_one_ulp", test_oracle_sensitivity), &
                  new_unittest("real_stringify_is_exact", test_real_exact), &
                  new_unittest("empty_and_missing_input", test_empty_missing), &
                  new_unittest("refuse_unknown_table", test_unknown_table), &
                  new_unittest("refuse_unknown_key_and_case", test_unknown_key_case), &
                  new_unittest("refuse_outside_table_and_non_table", test_outside_non_table), &
                  new_unittest("refuse_wrong_types", test_wrong_types), &
                  new_unittest("refuse_range_enum_array", test_range_enum_array), &
                  new_unittest("syntax_error_contract", test_syntax_error), &
                  new_unittest("errors_accumulate", test_accumulate), &
                  new_unittest("scalar_for_array_key", test_scalar_for_array), &
                  new_unittest("group_table_mapping", test_mapping), &
                  new_unittest("integer_overflow_not_wrapped", test_int_overflow), &
                  new_unittest("number_spellings", test_number_spellings), &
                  new_unittest("nonfinite_and_ranges", test_nonfinite), &
                  new_unittest("float_overflow_is_deterministic", test_float_overflow), &
                  new_unittest("duplicate_group_routes", test_duplicate_group), &
                  new_unittest("dotted_keys_and_inline_tables", test_dotted_inline), &
                  new_unittest("refuse_exotic_values", test_exotic_values), &
                  new_unittest("string_contents_and_length", test_strings), &
                  new_unittest("unreadable_file_label", test_unreadable) &
                  ]
   end subroutine collect_toml_config_tests

   ! ---- helpers --------------------------------------------------------

   subroutine slurp(path, text)
      !! Whole file into one string (newline-joined), then delete it.
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: text
      character(len=2048) :: buf
      integer :: u, ios

      text = ""
      open (newunit=u, file=path, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      do
         read (u, "(A)", iostat=ios) buf
         if (ios /= 0) exit
         text = text//trim(buf)//NL
      end do
      close (u, status="delete")
   end subroutine slurp

   subroutine dumps_of(text, dump_all, dump_exact, ierr)
      !! Resolved dumps (doc + full-precision) of a config held in `text`.
      character(len=*), intent(in) :: text
      character(len=:), allocatable, intent(out) :: dump_all, dump_exact
      integer, intent(out) :: ierr
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema

      call read_config_from_toml_string(text, cfg, schema, ierr)
      dump_all = ""
      dump_exact = ""
      if (ierr /= OCEAN_STATUS_OK) return
      call schema%write_doc_all("test_toml_dump.all")
      call schema%write_resolved("test_toml_dump.exact")
      call slurp("test_toml_dump.all", dump_all)
      call slurp("test_toml_dump.exact", dump_exact)
   end subroutine dumps_of

   subroutine check_real_value(error, token, expected)
      !! `[physics] coriolis_f = <token>` reads as exactly `expected`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: token
      real(real64), intent(in) :: expected
      type(config_t), target :: cfg
      integer :: ierr

      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = "//token//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "coriolis_f = "//token//" accepted")
      if (allocated(error)) return
      call check(error, real(cfg%coriolis_f, real64) == expected, "coriolis_f = "//token//" reads exactly")
   end subroutine check_real_value

   subroutine check_int_value(error, token, expected)
      !! `[grid] nx = <token>` reads as exactly `expected`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: token
      integer, intent(in) :: expected
      type(config_t), target :: cfg
      integer :: ierr

      call read_config_from_toml_string("[grid]"//NL//"nx = "//token//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "nx = "//token//" accepted")
      if (allocated(error)) return
      call check(error, cfg%nx, expected, "nx = "//token//" reads exactly")
   end subroutine check_int_value

   subroutine errors_of(text, status, errors)
      !! The schema's own error list for a TOML text, against a fresh schema.
      character(len=*), intent(in) :: text
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: errors(:)
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema

      call build_rdb_schema(cfg, schema)
      call toml_schema_parse_string(schema, text, status=status, errors=errors)
   end subroutine errors_of

   function joined(errors) result(s)
      character(len=:), allocatable, intent(in) :: errors(:)
      character(len=:), allocatable :: s
      integer :: i
      s = ""
      if (.not. allocated(errors)) return
      do i = 1, size(errors)
         s = s//trim(errors(i))//NL
      end do
   end function joined

   subroutine check_refused(error, text, needles)
      !! `text` is refused and the message contains every needle.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: text
      character(len=*), intent(in) :: needles(:)
      integer :: status, i
      character(len=:), allocatable :: errors(:), msg

      call errors_of(text, status, errors)
      msg = joined(errors)
      call check(error, status > 0, "refused: "//text)
      if (allocated(error)) return
      do i = 1, size(needles)
         call check(error, index(msg, trim(needles(i))) > 0, &
                    "message names '"//trim(needles(i))//"'; got: "//msg)
         if (allocated(error)) return
      end do
   end subroutine check_refused

   ! ---- equivalence ----------------------------------------------------

   subroutine test_all_key_types(error)
      !! Every key type, a real array, an integer for a real, an enum in
      !! non-canonical case, nested `ocean.*` tables and non-ocean tables.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      integer :: ierr

      call read_config_from_toml_string( &
         "[sim]"//NL//'sim_type = "ocean"'//NL// &
         "[grid]"//NL//"nx = 24"//NL//"ny = 16"//NL//"nghost = 3"//NL// &
         "[time]"//NL//"t_end = 2.5"//NL//'time_unit = "day"'//NL//"dt_fixed = 600"//NL// &
         "[physics]"//NL//"coriolis_f = 9.4e-5"//NL//"wind_stress_x = 0.1"//NL// &
         "[ocean.hvisc]"//NL//"nu_h = 10000"//NL//"smag_ah = true"//NL// &
         'lateral_closure = "Smagorinsky"'//NL//"c_smag = 0.15"//NL// &
         "[ocean.ic]"//NL//"layer_rho_init = [1036.0, 1035.0]"//NL// &
         "[ocean.vmix]"//NL//"dt_therm_ratio = 2"//NL// &
         "[ocean.diag]"//NL//'filename = "run_a"'//NL//"enabled = false"//NL// &
         "[ocean.ice_sw]"//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "accepted")
      if (allocated(error)) return
      call check(error, trim(cfg%sim_type) == "ocean", "string")
      if (allocated(error)) return
      call check(error, cfg%nx == 24 .and. cfg%ny == 16 .and. cfg%nghost == 3, "integers")
      if (allocated(error)) return
      call check(error, real(cfg%t_end, real64) == 2.5_real64*86400.0_real64, "t_end scaled by time_unit = day")
      if (allocated(error)) return
      call check(error, real(cfg%dt_fixed, real64) == 600.0_real64, "integer for a real key")
      if (allocated(error)) return
      call check(error, real(cfg%coriolis_f, real64) == 9.4e-5_real64, "real")
      if (allocated(error)) return
      call check(error, real(cfg%ocean%hvisc%nu_h, real64) == 10000.0_real64, "integer for a real key (ocean table)")
      if (allocated(error)) return
      call check(error, cfg%ocean%hvisc%smag_ah, "boolean true")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%hvisc%lateral_closure) == "smagorinsky", &
                 "enum canonicalised to its registered spelling")
      if (allocated(error)) return
      call check(error, real(cfg%ocean%ic%layer_rho_init(1), real64) == 1036.0_real64 .and. &
                 real(cfg%ocean%ic%layer_rho_init(2), real64) == 1035.0_real64, "real array")
      if (allocated(error)) return
      call check(error, cfg%ocean%vmix%dt_therm_ratio == 2, "integer in an ocean table")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%diag%filename) == "run_a" .and. .not. cfg%ocean%diag%enabled, &
                 "string and boolean false")
   end subroutine test_all_key_types

   subroutine test_double_gyre(error)
      !! The canonical double-gyre reference input (`validation_examples/ocean/
      !! double_gyre/double_gyre_mom6.toml`), comments dropped: every table class
      !! the reference uses resolves to the stated values.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      integer :: ierr

      call read_config_from_toml_string( &
         "[sim]"//NL//'sim_type = "ocean"'//NL// &
         "[grid]"//NL//"nx = 44"//NL//"ny = 40"//NL//"nghost = 3"//NL// &
         "[ocean.grid]"//NL//'axis_units = "degrees"'//NL//"len_lon = 22.0"//NL//"len_lat = 20.0"//NL// &
         "rad_earth = 6.378e6"//NL// &
         "[time]"//NL//"t_end = 10.0"//NL//'time_unit = "day"'//NL//"dt_fixed = 1200.0"//NL// &
         "cfl_interval = 1"//NL// &
         "[physics]"//NL//"coriolis_f = 9.4e-5"//NL//"wind_stress_x = 0.0"//NL//"wind_stress_y = 0.0"//NL// &
         "[nonhydrostatic]"//NL//"nz_layers = 2"//NL// &
         "[tracer]"//NL//"initial_temperature = 15.0"//NL//"initial_salinity = 35.0"//NL// &
         "[ocean.ic]"//NL//"layer_rho_init = [1036.0, 1035.0]"//NL// &
         "[ocean.coriolis]"//NL//'form = "sadourny_energy"'//NL// &
         "[ocean.thermo]"//NL//"enable_thermodynamics = false"//NL// &
         "[ocean.topo]"//NL//'topo_config = "spoon"'//NL//"max_depth = 2000.0"//NL// &
         "edge_depth = 100.0"//NL//"slope_scale = 400000.0"//NL//'wind_config = "2gyre"'//NL// &
         "taux_magnitude = 0.1"//NL//"coriolis_beta = 1.76e-11"//NL//"coriolis_y_ref = 1113200.0"//NL// &
         "[ocean.pgf]"//NL//'form = "gprime"'//NL//"gprime_gfs = 0.98"//NL//"gprime_gint = 0.0098"//NL// &
         "maxvel = 6.0"//NL// &
         "[ocean.bdrag]"//NL//'form = "linear"'//NL//"r = 2.5e-5"//NL//"hbbl = 10.0"//NL// &
         "bg_vel = 0.1"//NL//"bbl_thick_min = 0.1"//NL// &
         "[ocean.hvisc]"//NL//"nu_h = 10000.0"//NL//"nu_4 = 0.0"//NL//"smag_ah = true"//NL// &
         "smag_bi_const = 0.06"//NL//"nu_4_bg = 0.0"//NL//"nu_4_max = 1.0e12"//NL// &
         'lateral_closure = "smagorinsky"'//NL//"c_smag = 0.15"//NL//"kh_vel_scale = 3.0e-3"//NL// &
         "ah_max = 1.0e4"//NL// &
         "[ocean.vmix]"//NL//"kv_ml_invz2 = 0.01"//NL//"hmix_fixed = 20.0"//NL// &
         "direct_stress = false"//NL//"hmix_stress = 20.0"//NL//"harmonic_visc = true"//NL// &
         "dt_therm_ratio = 2"//NL// &
         "[ocean.bt]"//NL//"auto_n_inner = true"//NL//"cfl_bt_safety = 0.65"//NL//"bebt = 0.2"//NL// &
         "[vcoord]"//NL//'vcoord_type = "zstar_full"'//NL//"zstar_h_surf_target = 1000.0"//NL// &
         "zstar_h_min = 1.5e-4"//NL// &
         "[ocean.diag]"//NL//"enabled = true"//NL//'filename = "double_gyre"'//NL//"dt_out = 1.0"//NL// &
         "[output]"//NL//'output_dir = "./out_double_gyre_mom6"'//NL//"output_to_file = false"//NL// &
         "[logging]"//NL//"status_interval = 1.0"//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "accepted")
      if (allocated(error)) return
      call check(error, cfg%nx == 44 .and. cfg%ny == 40 .and. cfg%nghost == 3, "grid size")
      if (allocated(error)) return
      ! len_lon/len_lat in degrees derive dx/dy by the arc-length formula.
      call check(error, abs(real(cfg%dx, real64) - 6.378e6_real64*22.0_real64*DEG2RAD/44.0_real64) < 1.0e-6_real64, &
                 "dx derived from len_lon")
      if (allocated(error)) return
      call check(error, abs(real(cfg%dy, real64) - 6.378e6_real64*20.0_real64*DEG2RAD/40.0_real64) < 1.0e-6_real64, &
                 "dy derived from len_lat")
      if (allocated(error)) return
      call check(error, real(cfg%dt_fixed, real64) == 1200.0_real64 .and. cfg%ocean%vmix%dt_therm_ratio == 2, &
                 "time stepping")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%pgf%form) == "gprime" .and. trim(cfg%ocean%coriolis%form) == "sadourny_energy" &
                 .and. trim(cfg%ocean%bdrag%form) == "linear", "the three `form` keys land in their own tables")
      if (allocated(error)) return
      call check(error, real(cfg%ocean%bt%bebt, real64) == 0.2_real64 .and. cfg%ocean%bt%auto_n_inner, "bt table")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%hvisc%lateral_closure) == "smagorinsky" .and. cfg%ocean%hvisc%smag_ah &
                 .and. real(cfg%ocean%hvisc%nu_4_max, real64) == 1.0e12_real64, "hvisc table")
      if (allocated(error)) return
      call check(error, trim(cfg%vcoord_type) == "zstar_full" .and. trim(cfg%output_dir) == "./out_double_gyre_mom6" &
                 .and. .not. cfg%output_to_file, "vcoord and output tables")
      if (allocated(error)) return
      call check(error, real(cfg%ocean%ic%layer_rho_init(2), real64) == 1035.0_real64, "layer_rho_init")
   end subroutine test_double_gyre

   subroutine test_file_entry(error)
      !! `read_config` reads a `.toml` file; any other name is refused with a
      !! message naming TOML (never silently defaults), and a namelist file is
      !! recognised as such.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: c1
      type(nml_schema_t) :: s1
      character(len=:), allocatable :: a1
      integer :: u, i1

      open (newunit=u, file="test_toml_in.toml", status="replace", action="write")
      write (u, "(A)") "[grid]"
      write (u, "(A)") "nx = 30"
      write (u, "(A)") "[ocean.hvisc]"
      write (u, "(A)") "nu_h = 250.5"
      close (u)
      call read_config("test_toml_in.toml", c1, s1, i1)
      open (newunit=u, file="test_toml_in.toml", status="old")
      close (u, status="delete")
      call check(error, i1, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, c1%nx == 30 .and. real(c1%ocean%hvisc%nu_h, real64) == 250.5_real64, "file values resolved")
      if (allocated(error)) return
      call s1%write_resolved("test_toml_f1.exact")
      call slurp("test_toml_f1.exact", a1)
      call check(error, len(a1) > 0, "resolved dump written")
      if (allocated(error)) return

      ! A namelist file (and any other non-.toml name) is refused, not defaulted.
      open (newunit=u, file="test_toml_in.nml", status="replace", action="write")
      write (u, "(A)") "! a namelist"
      write (u, "(A)") "&grid_nml"
      write (u, "(A)") " nx = 30"
      write (u, "(A)") "/"
      close (u)
      call read_config("test_toml_in.nml", c1, s1, i1)
      open (newunit=u, file="test_toml_in.nml", status="old")
      close (u, status="delete")
      call check(error, i1, OCEAN_STATUS_ERR_CONFIG_PARSE, "a .nml file is refused")
      if (allocated(error)) return
      call check(error, c1%nx == 200, "refusal leaves the defaults untouched")
      if (allocated(error)) return
      call read_config("no_such_file_anywhere.nml", c1, s1, i1)
      call check(error, i1, OCEAN_STATUS_ERR_CONFIG_PARSE, "a missing non-.toml name is refused too")
      if (allocated(error)) return
      call read_config("input.txt", c1, s1, i1)
      call check(error, i1, OCEAN_STATUS_ERR_CONFIG_PARSE, "any non-.toml name is refused")
   end subroutine test_file_entry

   subroutine test_oracle_sensitivity(error)
      !! The equality oracle is not vacuous: the doc dump rounds reals to 11
      !! digits and so cannot tell two values one ulp apart, the full-precision
      !! dump can.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      integer :: i1, i2

      call dumps_of("[physics]"//NL//"coriolis_f = 1.0"//NL, a1, e1, i1)
      call dumps_of("[physics]"//NL//"coriolis_f = 1.0000000000000002"//NL, a2, e2, i2)
      call check(error, i1 == OCEAN_STATUS_OK .and. i2 == OCEAN_STATUS_OK, "both accepted")
      if (allocated(error)) return
      call check(error, e1 /= e2, "full-precision dump separates one ulp")
      if (allocated(error)) return
      call check(error, a1 == a2, "doc dump (11 digits) cannot")
   end subroutine test_oracle_sensitivity

   subroutine test_real_exact(error)
      !! A float written by the reader is read back bit for bit: sweep a few
      !! thousand bit patterns (every exponent, every mantissa shape) through
      !! TOML text -> toml-f -> stringify -> `read` and compare with the value
      !! the text was written from; then a table of known values against the
      !! compiler's own literals.
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: state, bits
      real(real64) :: x
      character(len=40) :: lit
      character(len=:), allocatable :: txt
      type(config_t), target :: cfg
      integer :: i, ierr, nbad, ntest

      state = 88172645463325252_int64
      nbad = 0
      ntest = 0
      do i = 1, 3000
         ! xorshift64
         state = ieor(state, ishft(state, 13))
         state = ieor(state, ishft(state, -7))
         state = ieor(state, ishft(state, 17))
         bits = state
         x = transfer(bits, x)
         ! Every 7th sample a subnormal, every 11th one near huge().
         if (mod(i, 7) == 0) x = transfer(iand(bits, 4503599627370495_int64), x)
         if (mod(i, 11) == 0) x = transfer(ior(iand(bits, 4503599627370495_int64), &
                                               9218868437227405312_int64 - 4503599627370496_int64), x)
         if (.not. ieee_is_finite(x)) cycle
         write (lit, "(ES30.17E3)") x
         txt = "[physics]"//NL//"coriolis_f = "//trim(adjustl(lit))//NL
         call read_config_from_toml_string(txt, cfg, ierr=ierr)
         ntest = ntest + 1
         if (ierr /= OCEAN_STATUS_OK) then
            nbad = nbad + 1
         else if (real(cfg%coriolis_f, real64) /= x) then
            nbad = nbad + 1
         end if
      end do
      call check(error, ntest > 2000, "enough finite samples")
      if (allocated(error)) return
      call check(error, nbad, 0, "every sampled real reads back bit for bit")
      if (allocated(error)) return

      ! Known values, against the compiler's own literal parser.
      call check_real_value(error, "0.1", 0.1_real64)
      if (allocated(error)) return
      call check_real_value(error, "0.3333333333333333", 0.3333333333333333_real64)
      if (allocated(error)) return
      call check_real_value(error, "1.7976931348623157e308", huge(1.0_real64))
      if (allocated(error)) return
      call check_real_value(error, "2.2250738585072014e-308", tiny(1.0_real64))
      if (allocated(error)) return
      call check_real_value(error, "4.9406564584124654e-324", transfer(1_int64, 1.0_real64))
      if (allocated(error)) return
      call check_real_value(error, "-6.378e6", -6.378e6_real64)
      if (allocated(error)) return
      call check_real_value(error, "9007199254740993", 9007199254740992.0_real64)
      if (allocated(error)) return
      call check_real_value(error, "1.0000000000000002", 1.0_real64 + epsilon(1.0_real64))
   end subroutine test_real_exact

   subroutine test_empty_missing(error)
      !! Empty text, a comment-only document and a missing `.toml` file keep the
      !! defaults.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: i1, i2, i3

      call dumps_of("", a1, e1, i1)
      call check(error, i1, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, len(a1) > 0 .and. len(e1) > 0, "the defaults dump is non-empty")
      if (allocated(error)) return
      call dumps_of("# only a comment"//NL, a2, e2, i2)
      call check(error, i2 == OCEAN_STATUS_OK .and. a2 == a1 .and. e2 == e1, "comment-only toml keeps defaults")
      if (allocated(error)) return
      call read_config("no_such_file_anywhere.toml", cfg, schema, i3)
      call check(error, i3, OCEAN_STATUS_OK, "missing .toml keeps defaults")
      if (allocated(error)) return
      call check(error, cfg%nx == 200, "missing .toml leaves the default nx")
   end subroutine test_empty_missing

   ! ---- refusals -------------------------------------------------------

   subroutine test_unknown_table(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[grid]"//NL//"nx = 3"//NL//NL//"[ocean.hvizc]"//NL//"nu_h = 1.0"//NL, &
                         [character(len=24) :: "<config-string>:4:", "unknown table [ocean.hvizc]", &
                          "ocean.hvisc"])
      if (allocated(error)) return
      call check_refused(error, "[nonesuch]"//NL//"a = 1"//NL, &
                         [character(len=24) :: ":1:", "unknown table [nonesuch]"])
      if (allocated(error)) return
      ! The group's one spelling is [ocean.hvisc]; the flat form is refused.
      call check_refused(error, "[ocean_hvisc]"//NL//"nu_h = 1.0"//NL, &
                         [character(len=24) :: ":1:", "unknown table [ocean_hvisc]", "ocean.hvisc"])
      if (allocated(error)) return
      ! A `_nml` suffix is not part of a table name.
      call check_refused(error, "[grid_nml]"//NL//"nx = 3"//NL, &
                         [character(len=24) :: "unknown table [grid_nml]"])
      if (allocated(error)) return
      ! An unknown ocean group and an unknown sub-level.
      call check_refused(error, "[ocean.hvisc.extra]"//NL//"nu_h = 1.0"//NL, &
                         [character(len=24) :: "unexpected sub-table 'extra'"])
   end subroutine test_unknown_table

   subroutine test_unknown_key_case(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[ocean.hvisc]"//NL//"nu_hh = 1.0"//NL, &
                         [character(len=24) :: ":2:", "unknown key 'nu_hh'", "[ocean.hvisc]", "nu_h"])
      if (allocated(error)) return
      ! Keys are exact-case: `NU_H` is NOT a key here.
      call check_refused(error, "[ocean.hvisc]"//NL//"NU_H = 1.0"//NL, &
                         [character(len=24) :: ":2:", "unknown key 'NU_H'", "did you mean 'nu_h'"])
      if (allocated(error)) return
      ! ... and neither is a wrong-case table.
      call check_refused(error, "[Grid]"//NL//"nx = 3"//NL, &
                         [character(len=24) :: "unknown table [Grid]", "did you mean 'grid'"])
   end subroutine test_unknown_key_case

   subroutine test_outside_non_table(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "nx = 3"//NL, &
                         [character(len=24) :: ":1:", "key 'nx' is outside any table"])
      if (allocated(error)) return
      call check_refused(error, "grid = 5"//NL, &
                         [character(len=24) :: "key 'grid' is outside any table"])
      if (allocated(error)) return
      call check_refused(error, "ocean = 3"//NL, &
                         [character(len=24) :: "key 'ocean' is outside any table"])
      if (allocated(error)) return
      call check_refused(error, "[ocean]"//NL//"nu_h = 1.0"//NL, &
                         [character(len=24) :: ":2:", "key 'nu_h' in [ocean] is not a table"])
      if (allocated(error)) return
      call check_refused(error, "[[grid]]"//NL//"nx = 3"//NL, &
                         [character(len=24) :: "key 'grid' is outside any table"])
   end subroutine test_outside_non_table

   subroutine test_wrong_types(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[grid]"//NL//'nx = "44"'//NL, &
                         [character(len=24) :: ":2:", "key 'nx'", "expects an integer, got a string"])
      if (allocated(error)) return
      call check_refused(error, "[grid]"//NL//"nx = 44.0"//NL, &
                         [character(len=24) :: "expects an integer, got a float"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//'smag_ah = "true"'//NL, &
                         [character(len=24) :: "expects a boolean (true/false), got a string"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//"smag_ah = 1"//NL, &
                         [character(len=24) :: "expects a boolean (true/false), got an integer"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//'nu_h = "big"'//NL, &
                         [character(len=24) :: "expects a number, got a string"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//"lateral_closure = 3"//NL, &
                         [character(len=24) :: "expects a string, got an integer"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//"nu_h = [1.0, 2.0]"//NL, &
                         [character(len=24) :: "expects a number, got an array"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//'layer_rho_init = [1.0, "x"]'//NL, &
                         [character(len=24) :: "element 2", "expects a number, got a string"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//"layer_rho_init = []"//NL, &
                         [character(len=24) :: "has no value"])
   end subroutine test_wrong_types

   subroutine test_range_enum_array(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[grid]"//NL//"nx = 0"//NL, &
                         [character(len=24) :: ":2:", "key 'nx'", "below min"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//'lateral_closure = "nonsense"'//NL, &
                         [character(len=24) :: ":2:", "not in allowed set"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.bt]"//NL//'split_scheme = "mom6_pc"'//NL, &
                         [character(len=24) :: ":2:", "RETIRED", "pred_corr"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//"layer_rho_init = ["// &
                         repeat("1.0, ", 70)//"1.0]"//NL, &
                         [character(len=24) :: "accepts at most"])
   end subroutine test_range_enum_array

   subroutine test_syntax_error(error)
      !! A TOML syntax error comes from toml-f, through the same status/ierr
      !! contract; nothing else is reported with it.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: ierr, status
      character(len=:), allocatable :: errors(:)

      call errors_of("[grid"//NL//"nx = 3"//NL, status, errors)
      call check(error, status, 1, "exactly the syntax error")
      if (allocated(error)) return
      call check(error, index(joined(errors), "TOML syntax error") > 0, joined(errors))
      if (allocated(error)) return
      call errors_of("[grid]"//NL//"nx = 3"//NL//"nx = 4"//NL, status, errors)
      call check(error, status > 0 .and. index(joined(errors), "TOML syntax error") > 0, &
                 "duplicate key is a toml-f error: "//joined(errors))
      if (allocated(error)) return
      call read_config_from_toml_string("[grid"//NL, cfg, schema, ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_CONFIG_PARSE, "ierr contract")
   end subroutine test_syntax_error

   subroutine test_accumulate(error)
      !! Several problems in one document are all reported.
      type(error_type), allocatable, intent(out) :: error
      integer :: status
      character(len=:), allocatable :: errors(:), msg

      call errors_of("[grid]"//NL//"nx = 0"//NL//"nz = 3"//NL//"[nonesuch]"//NL//"a = 1"//NL// &
                     "[ocean.hvisc]"//NL//"smag_ah = 2"//NL, status, errors)
      msg = joined(errors)
      call check(error, status, 4, "four independent refusals: "//msg)
      if (allocated(error)) return
      call check(error, index(msg, ":2:") > 0 .and. index(msg, ":3:") > 0 .and. &
                 index(msg, ":4:") > 0 .and. index(msg, ":7:") > 0, "each at its own line: "//msg)
   end subroutine test_accumulate

   subroutine test_scalar_for_array(error)
      !! A lone number is a legal value for an array key (it fills element 1 and
      !! leaves the rest alone); a one-element and a two-element array likewise.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg, def
      integer :: ierr

      call read_config_from_toml_string("[ocean.ic]"//NL//"layer_rho_init = 1030.5"//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, real(cfg%ocean%ic%layer_rho_init(1), real64) == 1030.5_real64, "scalar fills element 1")
      if (allocated(error)) return
      call check(error, cfg%ocean%ic%layer_rho_init(2) == def%ocean%ic%layer_rho_init(2), "element 2 untouched")
      if (allocated(error)) return
      call read_config_from_toml_string("[ocean.ic]"//NL//"layer_rho_init = [1030.5, 1031]"//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, real(cfg%ocean%ic%layer_rho_init(1), real64) == 1030.5_real64 .and. &
                 real(cfg%ocean%ic%layer_rho_init(2), real64) == 1031.0_real64, "two-element array")
   end subroutine test_scalar_for_array

   subroutine test_mapping(error)
      !! The group <-> table mapping round-trips on every registered group and
      !! collides on none.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      character(len=:), allocatable :: t
      integer :: i, j, n_ocean

      call build_rdb_schema(cfg, schema)
      call check(error, size(schema%groups) > 60, "registered groups")
      if (allocated(error)) return
      n_ocean = 0
      do i = 1, size(schema%groups)
         t = group_to_table(schema%groups(i)%name)
         call check(error, table_to_group(t) == schema%groups(i)%name, &
                    "round trip: "//schema%groups(i)%name//" -> "//t)
         if (allocated(error)) return
         if (index(t, "ocean.") == 1) n_ocean = n_ocean + 1
         ! No table path may be the prefix-parent of another (a table that
         ! holds both keys and sub-tables): the only parent is `ocean`.
         do j = 1, size(schema%groups)
            if (i == j) cycle
            call check(error, t /= group_to_table(schema%groups(j)%name), "unique table paths")
            if (allocated(error)) return
            call check(error, index(group_to_table(schema%groups(j)%name), t//".") /= 1, &
                       "no table nests under another group's table: "//t)
            if (allocated(error)) return
         end do
      end do
      call check(error, n_ocean > 40, "ocean_* groups map under [ocean]")
      if (allocated(error)) return
      call check(error, group_to_table("ocean_ice_sw") == "ocean.ice_sw")
      if (allocated(error)) return
      call check(error, group_to_table("ocean_cavity_dyn") == "ocean.cavity_dyn")
      if (allocated(error)) return
      call check(error, group_to_table("grid") == "grid" .and. group_to_table("vcoord") == "vcoord")
      if (allocated(error)) return
      call check(error, table_to_group("ocean") == "" .and. table_to_group("ocean_grid") == "" &
                 .and. table_to_group("a.b") == "" .and. table_to_group("ocean.a.b") == "", &
                 "non-canonical paths name no group")
   end subroutine test_mapping

   ! ---- review round ----------------------------------------------------

   subroutine test_int_overflow(error)
      !! toml-f wraps an overflowing integer silently; the reader hands the
      !! schema the source digits instead, so the schema sees (and refuses) the
      !! real magnitude for an integer key and reads the real value for a real
      !! key, identically on every compiler.
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[grid]"//NL//"nx = 18446744073709551660"//NL, &
                         [character(len=24) :: ":2:", "key 'nx'", "cannot parse integer"])
      if (allocated(error)) return
      call check_refused(error, "[grid]"//NL//"nx = 99999999999"//NL, &
                         [character(len=24) :: ":2:", "key 'nx'"])
      if (allocated(error)) return
      ! An integer literal for a REAL key is read as a real on every compiler
      ! (nvfortran's `read` refuses a bare integer beyond int64 for a real, gfortran
      ! reads it as 1e20; the reader hands over `<digits>.0`) -- never toml-f's
      ! wrapped value.
      call check_real_value(error, "100000000000000000000", 1.0e20_real64)
      if (allocated(error)) return
      call check_real_value(error, "-9223372036854775809", -9.223372036854775808e18_real64)
      if (allocated(error)) return
      call check_real_value(error, "9223372036854775807", 9.223372036854775808e18_real64)
      if (allocated(error)) return
      call check_real_value(error, "1_000_000", 1.0e6_real64)
      if (allocated(error)) return
      call check_real_value(error, "-7", -7.0_real64)
   end subroutine test_int_overflow

   subroutine test_number_spellings(error)
      !! `-0` keeps its sign for a real key; underscores, a leading `+` and an
      !! exponent without a dot read as the digits say; hex integers are refused.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      integer :: ierr

      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = -0"//NL, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, sign(1.0_real64, real(cfg%coriolis_f, real64)) < 0.0_real64, "integer -0 for a real keeps its sign")
      if (allocated(error)) return
      call check_int_value(error, "1_000", 1000)
      if (allocated(error)) return
      call check_real_value(error, "1_500.2_5", 1500.25_real64)
      if (allocated(error)) return
      call check_real_value(error, "+1.5E-4", 1.5e-4_real64)
      if (allocated(error)) return
      call check_real_value(error, "3e3", 3000.0_real64)
      if (allocated(error)) return
      call check_refused(error, "[grid]"//NL//"nx = 0x10"//NL, &
                         [character(len=24) :: ":2:", "hexadecimal"])
   end subroutine test_number_spellings

   subroutine test_nonfinite(error)
      !! inf/nan read as those values for an unranged key and are refused for a
      !! key with a declared range.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      integer :: ierr
      real(real64) :: v

      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = inf"//NL, cfg, ierr=ierr)
      v = cfg%coriolis_f
      call check(error, ierr == OCEAN_STATUS_OK .and. .not. ieee_is_finite(v) .and. v > 0.0_real64, "inf")
      if (allocated(error)) return
      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = -inf"//NL, cfg, ierr=ierr)
      v = cfg%coriolis_f
      call check(error, ierr == OCEAN_STATUS_OK .and. .not. ieee_is_finite(v) .and. v < 0.0_real64, "-inf")
      if (allocated(error)) return
      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = +inf"//NL, cfg, ierr=ierr)
      v = cfg%coriolis_f
      call check(error, ierr == OCEAN_STATUS_OK .and. .not. ieee_is_finite(v) .and. v > 0.0_real64, "+inf")
      if (allocated(error)) return
      call read_config_from_toml_string("[physics]"//NL//"coriolis_f = nan"//NL, cfg, ierr=ierr)
      v = cfg%coriolis_f
      call check(error, ierr == OCEAN_STATUS_OK .and. ieee_is_nan(v), "nan")
      if (allocated(error)) return
      call check_refused(error, "[ocean.bt]"//NL//"pc_be = nan"//NL, &
                         [character(len=24) :: ":2:", "key 'pc_be'", "not finite"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.bt]"//NL//"pc_be = inf"//NL, &
                         [character(len=24) :: "not finite"])
   end subroutine test_nonfinite

   subroutine test_float_overflow(error)
      !! A literal beyond the real range gives the same answer on every
      !! compiler: both runs of the reader agree (the same value, or both
      !! refuse), and a refusal carries the key and line.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      integer :: i1, i2

      call dumps_of("[physics]"//NL//"coriolis_f = 1e400"//NL, a1, e1, i1)
      call dumps_of("[physics]"//NL//"coriolis_f = 1e400"//NL, a2, e2, i2)
      call check(error, i1 == i2, "deterministic accept/refuse verdict for 1e400")
      if (allocated(error)) return
      if (i1 == OCEAN_STATUS_OK) then
         call check(error, e1 == e2, "deterministic value for 1e400")
      else
         call check_refused(error, "[physics]"//NL//"coriolis_f = 1e400"//NL, &
                            [character(len=24) :: ":2:", "key 'coriolis_f'"])
      end if
   end subroutine test_float_overflow

   subroutine test_duplicate_group(error)
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, '["ocean.hvisc"]'//NL//"nu_h = 1.0"//NL, &
                         [character(len=24) :: ":1:", "contains a dot", "[ocean.<group>]"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//"nu_h = 1.0"//NL//'["ocean.hvisc"]'//NL// &
                         "nu_h = 2.0"//NL, [character(len=24) :: ":3:", "contains a dot"])
      if (allocated(error)) return
      call check_refused(error, '["grid.x"]'//NL//"a = 1"//NL, &
                         [character(len=24) :: "contains a dot"])
   end subroutine test_duplicate_group

   subroutine test_dotted_inline(error)
      !! Dotted keys and inline tables build the same tables as headers do.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      integer :: i1, i2

      call dumps_of("[ocean.hvisc]"//NL//"nu_h = 5.0"//NL//"[grid]"//NL//"nx = 20"//NL, a1, e1, i1)
      call check(error, i1, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call dumps_of("[ocean]"//NL//"hvisc.nu_h = 5.0"//NL//"[grid]"//NL//"nx = 20"//NL, a2, e2, i2)
      call check(error, i2 == OCEAN_STATUS_OK .and. a2 == a1 .and. e2 == e1, "dotted key under [ocean]")
      if (allocated(error)) return
      call dumps_of("grid = { nx = 20 }"//NL//"[ocean]"//NL//"hvisc = { nu_h = 5.0 }"//NL, a2, e2, i2)
      call check(error, i2 == OCEAN_STATUS_OK .and. a2 == a1 .and. e2 == e1, "inline tables")
      if (allocated(error)) return
      call dumps_of("ocean.hvisc.nu_h = 5.0"//NL//"grid.nx = 20"//NL, a2, e2, i2)
      call check(error, i2 == OCEAN_STATUS_OK .and. a2 == a1 .and. e2 == e1, "fully dotted keys")
      if (allocated(error)) return
      call check_int_value(error, "20", 20)
   end subroutine test_dotted_inline

   subroutine test_exotic_values(error)
      !! Datetimes, nested / mixed / empty arrays are refused with a line.
      type(error_type), allocatable, intent(out) :: error
      call check_refused(error, "[grid]"//NL//"nx = 1979-05-27"//NL, &
                         [character(len=24) :: ":2:", "datetime"])
      if (allocated(error)) return
      call check_refused(error, "[physics]"//NL//"coriolis_f = 1979-05-27T07:32:00Z"//NL, &
                         [character(len=24) :: ":2:", "datetime"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//"layer_rho_init = [[1.0], [2.0]]"//NL, &
                         [character(len=24) :: ":2:", "element 1"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//'layer_rho_init = [1.0, "a", 2.0]'//NL, &
                         [character(len=24) :: ":2:", "element 2"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//"layer_rho_init = []"//NL, &
                         [character(len=24) :: ":2:", "has no value"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.ic]"//NL//"layer_rho_init = [1979-05-27]"//NL, &
                         [character(len=24) :: "element 1", "datetime"])
   end subroutine test_exotic_values

   subroutine test_strings(error)
      !! Special characters survive; control characters and over-long strings
      !! are refused.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      character(len=:), allocatable :: want, toml, longs
      character(len=1), parameter :: BS = achar(92)
      integer :: ierr

      want = " lead, a/b & c! d# 'q' e"
      toml = "[output]"//NL//'output_dir = " lead, a/b & c! d# ''q'' e"'//NL
      call read_config_from_toml_string(toml, cfg, schema, ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "special characters accepted")
      if (allocated(error)) return
      call check(error, trim(cfg%output_dir) == want, "special characters preserved: '"// &
                 trim(cfg%output_dir)//"'")
      if (allocated(error)) return
      call check_refused(error, "[output]"//NL//'output_dir = "a'//BS//'nb"'//NL, &
                         [character(len=24) :: ":2:", "control character"])
      if (allocated(error)) return
      call check_refused(error, "[output]"//NL//'output_dir = "a'//BS//'tb"'//NL, &
                         [character(len=24) :: ":2:", "control character"])
      if (allocated(error)) return
      call check_refused(error, "[ocean.hvisc]"//NL//'lateral_closure = "none'//BS//'n"'//NL, &
                         [character(len=24) :: "control character"])
      if (allocated(error)) return
      longs = repeat("x", 300)
      call check_refused(error, "[output]"//NL//'output_dir = "'//longs//'"'//NL, &
                         [character(len=24) :: ":2:", "exceeds target length"])
      if (allocated(error)) return
      longs = repeat("y", 200)
      call read_config_from_toml_string("[output]"//NL//'output_dir = "'//longs//'"'//NL, cfg, schema, ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. trim(cfg%output_dir) == longs, "200 chars fit")
      if (allocated(error)) return
   end subroutine test_strings

   subroutine test_unreadable(error)
      !! A path that cannot be opened is not reported as a TOML syntax error.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      character(len=:), allocatable :: errors(:)
      integer :: status

      call build_rdb_schema(cfg, schema)
      call toml_schema_parse_file(schema, ".", status=status, errors=errors)
      call check(error, status, 1)
      if (allocated(error)) return
      call check(error, index(joined(errors), "syntax") == 0, "not labelled a syntax error: "//joined(errors))
   end subroutine test_unreadable

end module test_toml_config
