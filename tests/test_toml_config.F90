!! Unit tests for the TOML configuration front end (rdb_toml_read).
module test_toml_config
   !! The TOML reader is a second front end on the strict namelist schema.  These
   !! tests pin (1) that a configuration spelled in TOML resolves to EXACTLY the
   !! configuration the namelist spells -- compared through the schema's
   !! resolved dumps, including a full-precision one -- (2) every refusal class
   !! of the TOML syntax, with the file/line and key in the message, and (3) the
   !! group <-> table mapping.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_config, only: config_t, read_config, read_config_from_string, &
                         read_config_from_toml_string, build_rdb_schema
   use rdb_nml_schema, only: nml_schema_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_PARSE
   use rdb_toml_read, only: group_to_table, table_to_group, toml_schema_parse_string
   implicit none
   private

   public :: collect_toml_config_tests

   character(len=*), parameter :: NL = achar(10)

contains

   subroutine collect_toml_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("equivalent_all_key_types", test_equivalent_all_types), &
                  new_unittest("equivalent_double_gyre", test_equivalent_double_gyre), &
                  new_unittest("equivalent_file_entry", test_equivalent_file), &
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
                  new_unittest("group_table_mapping", test_mapping) &
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

   subroutine dumps_of(text, as_toml, dump_all, dump_exact, ierr)
      !! Resolved dumps (doc + full-precision) of a config held in `text`,
      !! read as TOML or as namelist.
      character(len=*), intent(in) :: text
      logical, intent(in) :: as_toml
      character(len=:), allocatable, intent(out) :: dump_all, dump_exact
      integer, intent(out) :: ierr
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema

      if (as_toml) then
         call read_config_from_toml_string(text, cfg, schema, ierr)
      else
         call read_config_from_string(text, cfg, schema, ierr)
      end if
      dump_all = ""
      dump_exact = ""
      if (ierr /= OCEAN_STATUS_OK) return
      call schema%write_doc_all("test_toml_dump.all")
      call schema%write_resolved("test_toml_dump.exact")
      call slurp("test_toml_dump.all", dump_all)
      call slurp("test_toml_dump.exact", dump_exact)
   end subroutine dumps_of

   subroutine check_equivalent(error, nml_text, toml_text)
      !! Both texts must resolve to byte-identical dumps.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml_text, toml_text
      character(len=:), allocatable :: a1, e1, a2, e2
      integer :: i1, i2

      call dumps_of(nml_text, .false., a1, e1, i1)
      call dumps_of(toml_text, .true., a2, e2, i2)
      call check(error, i1, OCEAN_STATUS_OK, "namelist text accepted")
      if (allocated(error)) return
      call check(error, i2, OCEAN_STATUS_OK, "toml text accepted")
      if (allocated(error)) return
      call check(error, len(a1) > 0 .and. len(e1) > 0, "dumps are non-empty")
      if (allocated(error)) return
      call check(error, a2 == a1, "doc dumps identical")
      if (allocated(error)) return
      call check(error, e2 == e1, "full-precision dumps identical")
   end subroutine check_equivalent

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

   subroutine test_equivalent_all_types(error)
      !! Every key type, a real array, an integer for a real, an enum in
      !! non-canonical case, nested `ocean.*` tables and non-ocean tables.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: n, t

      n = "&sim_nml"//NL//' sim_type = "ocean"'//NL//"/"//NL// &
          "&grid_nml"//NL//" nx = 24"//NL//" ny = 16"//NL//" nghost = 3"//NL//"/"//NL// &
          "&time_nml"//NL//" t_end = 2.5"//NL//' time_unit = "day"'//NL// &
          " dt_fixed = 600"//NL//"/"//NL// &
          "&physics_nml"//NL//" coriolis_f = 9.4e-5"//NL//" wind_stress_x = 0.1"//NL//"/"//NL// &
          "&ocean_hvisc_nml"//NL//" nu_h = 10000"//NL//" smag_ah = .true."//NL// &
          ' lateral_closure = "Smagorinsky"'//NL//" c_smag = 0.15"//NL//"/"//NL// &
          "&ocean_ic_nml"//NL//" layer_rho_init = 1036.0, 1035.0"//NL//"/"//NL// &
          "&ocean_vmix_nml"//NL//" dt_therm_ratio = 2"//NL//"/"//NL// &
          "&ocean_diag_nml"//NL//' filename = "run_a"'//NL//" enabled = .false."//NL//"/"//NL// &
          "&ocean_ice_sw_nml"//NL//"/"//NL
      t = "[sim]"//NL//'sim_type = "ocean"'//NL// &
          "[grid]"//NL//"nx = 24"//NL//"ny = 16"//NL//"nghost = 3"//NL// &
          "[time]"//NL//"t_end = 2.5"//NL//'time_unit = "day"'//NL//"dt_fixed = 600"//NL// &
          "[physics]"//NL//"coriolis_f = 9.4e-5"//NL//"wind_stress_x = 0.1"//NL// &
          "[ocean.hvisc]"//NL//"nu_h = 10000"//NL//"smag_ah = true"//NL// &
          'lateral_closure = "Smagorinsky"'//NL//"c_smag = 0.15"//NL// &
          "[ocean.ic]"//NL//"layer_rho_init = [1036.0, 1035.0]"//NL// &
          "[ocean.vmix]"//NL//"dt_therm_ratio = 2"//NL// &
          "[ocean.diag]"//NL//'filename = "run_a"'//NL//"enabled = false"//NL// &
          "[ocean.ice_sw]"//NL
      call check_equivalent(error, n, t)
   end subroutine test_equivalent_all_types

   subroutine test_equivalent_double_gyre(error)
      !! `validation_examples/ocean/double_gyre/double_gyre_mom6.nml`,
      !! hand-converted (comments dropped).
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: n, t

      n = &
         '&sim_nml sim_type = "ocean" /'//NL// &
         "&grid_nml nx = 44 ny = 40 nghost = 3 /"//NL// &
         '&ocean_grid_nml axis_units = "degrees" len_lon = 22.0 len_lat = 20.0 rad_earth = 6.378e6 /'//NL// &
         '&time_nml t_end = 10.0 time_unit = "day" dt_fixed = 1200.0 cfl_interval = 1 /'//NL// &
         "&physics_nml coriolis_f = 9.4e-5 wind_stress_x = 0.0 wind_stress_y = 0.0 /"//NL// &
         "&nonhydrostatic_nml nz_layers = 2 /"//NL// &
         "&tracer_nml initial_temperature = 15.0 initial_salinity = 35.0 /"//NL// &
         "&ocean_ic_nml layer_rho_init = 1036.0, 1035.0 /"//NL// &
         '&ocean_coriolis_nml form = "sadourny_energy" /'//NL// &
         "&ocean_thermo_nml enable_thermodynamics = .false. /"//NL// &
         '&ocean_topo_nml topo_config = "spoon" max_depth = 2000.0 edge_depth = 100.0 '// &
         'slope_scale = 400000.0 wind_config = "2gyre" taux_magnitude = 0.1 '// &
         "coriolis_beta = 1.76e-11 coriolis_y_ref = 1113200.0 /"//NL// &
         '&ocean_pgf_nml form = "gprime" gprime_gfs = 0.98 gprime_gint = 0.0098 maxvel = 6.0 /'//NL// &
         '&ocean_bdrag_nml form = "linear" r = 2.5e-5 hbbl = 10.0 bg_vel = 0.1 bbl_thick_min = 0.1 /'//NL// &
         "&ocean_hvisc_nml nu_h = 10000.0 nu_4 = 0.0 smag_ah = .true. smag_bi_const = 0.06 "// &
         'nu_4_bg = 0.0 nu_4_max = 1.0e12 lateral_closure = "smagorinsky" c_smag = 0.15 '// &
         "kh_vel_scale = 3.0e-3 ah_max = 1.0e4 /"//NL// &
         "&ocean_vmix_nml kv_ml_invz2 = 0.01 hmix_fixed = 20.0 direct_stress = .false. "// &
         "hmix_stress = 20.0 harmonic_visc = .true. dt_therm_ratio = 2 /"//NL// &
         "&ocean_bt_nml auto_n_inner = .true. cfl_bt_safety = 0.65 bebt = 0.2 /"//NL// &
         '&vcoord_nml vcoord_type = "zstar_full" zstar_h_surf_target = 1000.0 zstar_h_min = 1.5e-4 /'//NL// &
         '&ocean_diag_nml enabled = .true. filename = "double_gyre" dt_out = 1.0 /'//NL// &
         '&output_nml output_dir = "./out_double_gyre_mom6" output_to_file = .false. /'//NL// &
         "&logging_nml status_interval = 1.0 /"//NL
      t = &
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
         "[logging]"//NL//"status_interval = 1.0"//NL
      call check_equivalent(error, n, t)
   end subroutine test_equivalent_double_gyre

   subroutine test_equivalent_file(error)
      !! `read_config` dispatches on the `.toml` extension; the file form
      !! resolves like the namelist file form.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: c1, c2
      type(nml_schema_t) :: s1, s2
      character(len=:), allocatable :: a1, a2
      integer :: u, i1, i2

      open (newunit=u, file="test_toml_in.nml", status="replace", action="write")
      write (u, "(A)") "&grid_nml"
      write (u, "(A)") " nx = 30"
      write (u, "(A)") "/"
      write (u, "(A)") "&ocean_hvisc_nml"
      write (u, "(A)") " nu_h = 250.5"
      write (u, "(A)") "/"
      close (u)
      open (newunit=u, file="test_toml_in.toml", status="replace", action="write")
      write (u, "(A)") "[grid]"
      write (u, "(A)") "nx = 30"
      write (u, "(A)") "[ocean.hvisc]"
      write (u, "(A)") "nu_h = 250.5"
      close (u)

      call read_config("test_toml_in.nml", c1, s1, i1)
      call read_config("test_toml_in.toml", c2, s2, i2)
      call check(error, i1, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, i2, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call s1%write_resolved("test_toml_f1.exact")
      call s2%write_resolved("test_toml_f2.exact")
      call slurp("test_toml_f1.exact", a1)
      call slurp("test_toml_f2.exact", a2)
      call check(error, len(a1) > 0 .and. a1 == a2, "file readers resolve identically")
      open (newunit=u, file="test_toml_in.nml", status="old")
      close (u, status="delete")
      open (newunit=u, file="test_toml_in.toml", status="old")
      close (u, status="delete")
   end subroutine test_equivalent_file

   subroutine test_oracle_sensitivity(error)
      !! The equality oracle is not vacuous: the doc dump rounds reals to 11
      !! digits and so cannot tell two values one ulp apart, the full-precision
      !! dump can.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      integer :: i1, i2

      call dumps_of("[physics]"//NL//"coriolis_f = 1.0"//NL, .true., a1, e1, i1)
      call dumps_of("[physics]"//NL//"coriolis_f = 1.0000000000000002"//NL, .true., a2, e2, i2)
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
      !! the namelist path reads from the SAME decimal text.
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: state, bits
      real(real64) :: x, got_toml, got_nml
      character(len=40) :: lit
      character(len=:), allocatable :: txt
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
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
         if (abs(x) > 1.0e300_real64 .or. abs(x) < 1.0e-300_real64) cycle
         if (x /= x) cycle
         write (lit, "(ES30.17E3)") x
         txt = "[physics]"//NL//"coriolis_f = "//trim(adjustl(lit))//NL
         call read_config_from_toml_string(txt, cfg, schema, ierr)
         if (ierr /= OCEAN_STATUS_OK) then
            nbad = nbad + 1
            cycle
         end if
         got_toml = cfg%coriolis_f
         txt = "&physics_nml coriolis_f = "//trim(adjustl(lit))//" /"//NL
         call read_config_from_string(txt, cfg, schema, ierr)
         got_nml = cfg%coriolis_f
         ntest = ntest + 1
         if (got_toml /= got_nml) nbad = nbad + 1
      end do
      call check(error, ntest > 2000, "enough finite samples")
      if (allocated(error)) return
      call check(error, nbad, 0, "every sampled real is identical across the readers")
   end subroutine test_real_exact

   subroutine test_empty_missing(error)
      !! Empty text and a missing file keep the defaults, as the namelist
      !! path does.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: a1, e1, a2, e2
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: i1, i2, i3

      call dumps_of("", .false., a1, e1, i1)
      call dumps_of("", .true., a2, e2, i2)
      call check(error, i1, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, i2, OCEAN_STATUS_OK)
      if (allocated(error)) return
      call check(error, a1 == a2 .and. e1 == e2, "empty text keeps defaults in both")
      if (allocated(error)) return
      call dumps_of("# only a comment"//NL, .true., a2, e2, i2)
      call check(error, i2 == OCEAN_STATUS_OK .and. a2 == a1, "comment-only toml keeps defaults")
      if (allocated(error)) return
      call read_config("no_such_file_anywhere.toml", cfg, schema, i3)
      call check(error, i3, OCEAN_STATUS_OK, "missing .toml keeps defaults")
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
      ! A namelist-style suffix is not a table name either.
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
      ! Keys are exact-case: the namelist spelling `NU_H` is NOT a key here.
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
      !! A lone number is a legal value for an array key, as a single value
      !! is in a namelist (it fills element 1 and leaves the rest alone).
      type(error_type), allocatable, intent(out) :: error
      call check_equivalent(error, "&ocean_ic_nml layer_rho_init = 1030.5 /"//NL, &
                            "[ocean.ic]"//NL//"layer_rho_init = 1030.5"//NL)
      if (allocated(error)) return
      call check_equivalent(error, "&ocean_ic_nml layer_rho_init = 1030.5, 1031 /"//NL, &
                            "[ocean.ic]"//NL//"layer_rho_init = [1030.5, 1031]"//NL)
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

end module test_toml_config
