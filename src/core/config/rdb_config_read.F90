submodule(rdb_config) rdb_config_read
   !! TOML file / string reading entry points + post-parse cascades.
   use rdb_toml_read, only: toml_schema_parse_file, toml_schema_parse_string
   implicit none

   integer, parameter :: SRC_TOML_FILE = 1
   integer, parameter :: SRC_TOML_STRING = 2

contains

   module subroutine read_config(filename, cfg, schema, ierr)
      !! Read the simulation configuration from the TOML file `filename`.
      !! The name must end in `.toml` (case-insensitive): any other name is
      !! refused with a message naming TOML (and, when the file is a namelist,
      !! saying so) -- it never falls back to defaults.  The document is
      !! validated and applied through the strict schema (`rdb_nml_schema`)
      !! built from `cfg`'s pristine field initialisers, so `error stop` on
      !! any unknown key / range / enum violation fires for every caller
      !! (driver, FFI, unit tests), not just the production entry points.
      !!
      !! When `schema` is present the built schema is returned for the
      !! caller's parameter-doc dumps.  Its key pointers alias into
      !! `cfg`; both entry points keep their `cfg` alive (and `target`),
      !! so returning it is safe.
      character(len=*), intent(in) :: filename
         !! Path to the TOML input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out), optional :: schema
         !! Built + parsed schema (aliases into `cfg`); for doc dumps.
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure (unknown group/key,
         !! range/enum violation) when present; absent behaves as today
         !! (`error stop`).

      ! Dispatch so the schema is built DIRECTLY into the caller's
      ! variable when requested — never copied.  Intrinsic assignment of
      ! nml_schema_t deep-copies polymorphic allocatable key boxes, which
      ! NVHPC miscompiles (heap corruption in the FFI tests); building in
      ! place sidesteps it and is cheaper anyway.
      if (present(schema)) then
         call read_config_impl(filename, cfg, schema, ierr=ierr)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_impl(filename, cfg, sch, ierr=ierr)
         end block
      end if

   end subroutine read_config

   module subroutine read_config_impl(filename, cfg, sch, ierr)
      !! Body of [[read_config]]: build the strict schema (capturing
      !! defaults from the pristine `cfg` field initialisers), parse the
      !! file (validate + apply straight into `cfg`), run the post-parse
      !! time-unit cascade.
      character(len=*), intent(in) :: filename
         !! Path to the TOML input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out) :: sch
         !! Built + parsed schema (aliases into `cfg`)
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure or a missing input file when present; absent
         !! behaves as today (`error stop`).

      ! `cfg` arrives at its pristine field-initialiser defaults.  Build the
      ! strict schema (its key constructors capture those defaults), then parse
      ! -- which validates AND applies every override straight into `cfg` via
      ! the schema's key pointers.
      call build_rdb_schema(cfg, sch)

      ! Only TOML is read.  Refuse any other name up front, whether or not the
      ! file exists, rather than guessing a format or silently keeping defaults.
      if (.not. is_toml_name(filename)) then
         call refuse_non_toml(filename, ierr)
         return
      end if

      ! A missing / unopenable input file is an ERROR (never "keep defaults"):
      ! refused through the parse-failure contract -- `ierr` =
      ! `OCEAN_STATUS_ERR_CONFIG_PARSE` plus a message naming the path when
      ! `ierr` is present, `error stop` otherwise.  A caller that wants the
      ! defaults uses `read_config_from_toml_string("")` (an empty TOML document
      ! is valid and means "all defaults") or default-constructs `config_t`.
      block
         logical :: file_exists
         character(len=:), allocatable :: msg
         inquire (file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            msg = trim(filename)//": cannot read the input file (it does not exist or is not readable)"
            call error_ring_push(msg)
            call logger%error(msg)
            if (present(ierr)) then
               ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
               return
            end if
            error stop "read_config: the input file does not exist"
         end if
      end block

      ! Strict schema parse: validate + APPLY every group straight into
      ! `cfg`; routing of `ierr` / `error stop` is in `parse_into_schema`.
      block
         logical :: failed
         call parse_into_schema(sch, SRC_TOML_FILE, trim(filename), ierr, failed)
         if (failed) return
      end block

      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK

   end subroutine read_config_impl

   module subroutine read_config_from_toml_string(text, cfg, schema, ierr)
      !! In-memory sibling of [[read_config]]: the TOML document is held in
      !! `text` and nothing touches the filesystem.  Same `ierr` contract and
      !! post-parse cascades as the file reader.
      character(len=*), intent(in) :: text
         !! Whole TOML document as one string.
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out), optional :: schema
      integer, intent(out), optional :: ierr
         !! Non-zero on a refused document when present; absent behaves as
         !! the file reader does (`error stop`).

      if (present(schema)) then
         call read_config_from_toml_string_impl(text, cfg, schema, ierr=ierr)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_from_toml_string_impl(text, cfg, sch, ierr=ierr)
         end block
      end if
   end subroutine read_config_from_toml_string

   module subroutine read_config_from_toml_string_impl(text, cfg, sch, ierr)
      !! Body of [[read_config_from_toml_string]].
      character(len=*), intent(in) :: text
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out) :: sch
      integer, intent(out), optional :: ierr

      logical :: failed

      call build_rdb_schema(cfg, sch)
      if (len_trim(text) == 0) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if
      call parse_into_schema(sch, SRC_TOML_STRING, text, ierr, failed)
      if (failed) return
      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine read_config_from_toml_string_impl

   subroutine refuse_non_toml(filename, ierr)
      !! A non-`.toml` input name: log why (naming TOML; when the file reads as a
      !! namelist, saying namelist input is gone), then return
      !! `OCEAN_STATUS_ERR_CONFIG_PARSE` when `ierr` is present, else `error stop`.
      character(len=*), intent(in) :: filename
      integer, intent(out), optional :: ierr

      character(len=:), allocatable :: msg

      msg = trim(filename)//": the input file must be TOML (a '.toml' file name)"
      if (looks_like_namelist(filename)) then
         msg = msg//"; this file is a Fortran namelist (&..._nml groups), which is no longer read "// &
               "and has no converter -- rewrite it as TOML: &ocean_<x>_nml becomes [ocean.<x>], "// &
               "any other &<x>_nml becomes [<x>] (see 'Input file format' in docs/REFERENCE.md)"
      end if
      call error_ring_push(msg)
      call logger%error(msg)
      if (present(ierr)) then
         ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
      else
         error stop "read_config: the input file must be TOML"
      end if
   end subroutine refuse_non_toml

   function looks_like_namelist(filename) result(yes)
      !! True when the first non-blank, non-comment line of `filename` opens a
      !! namelist group (`&<name>_nml`).  False for an unreadable file.
      character(len=*), intent(in) :: filename
      logical :: yes

      character(len=512) :: line
      character(len=128) :: imsg
      integer :: u, ios, i, k

      yes = .false.
      open (newunit=u, file=trim(filename), status="old", action="read", iostat=ios, iomsg=imsg)
      if (ios /= 0) return
      do k = 1, 200
         read (u, "(A)", iostat=ios, iomsg=imsg) line
         if (ios /= 0) exit
         i = verify(line, " "//achar(9))
         if (i == 0) cycle
         if (line(i:i) == "!" .or. line(i:i) == "#") cycle
         yes = line(i:i) == "&" .and. index(line(i:), "_nml") > 0
         exit
      end do
      close (u)
   end function looks_like_namelist

   pure function is_toml_name(filename) result(yes)
      !! True when `filename` ends in `.toml` (any case).
      character(len=*), intent(in) :: filename
      logical :: yes
      integer :: n

      n = len_trim(filename)
      yes = .false.
      if (n > 5) yes = to_lower(filename(n - 4:n)) == ".toml"
   end function is_toml_name

   subroutine parse_into_schema(sch, src, text, ierr, failed)
      !! Parse + validate + apply `text` (a file path for the `*_FILE` kinds,
      !! the document for the `*_STRING` kinds) into `sch` through the front
      !! end `src` selects.  When `ierr` is present a refusal is routed
      !! through the schema's own non-aborting `status`/`errors` pair: every
      !! message is logged and pushed to the error ring, `ierr` is set to
      !! `OCEAN_STATUS_ERR_CONFIG_PARSE` and `failed` is true.  Absent
      !! `ierr` keeps the bare strict form (log + `error stop`), for every
      !! existing caller (driver, benchmarks, unit tests).
      type(nml_schema_t), intent(inout) :: sch
      integer, intent(in) :: src
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: ierr
      logical, intent(out) :: failed

      integer :: parse_status, e
      character(len=:), allocatable :: parse_errors(:)

      failed = .false.
      if (present(ierr)) then
         if (src == SRC_TOML_FILE) then
            call toml_schema_parse_file(sch, text, status=parse_status, errors=parse_errors)
         else
            call toml_schema_parse_string(sch, text, status=parse_status, errors=parse_errors)
         end if
         if (parse_status /= 0) then
            if (allocated(parse_errors)) then
               do e = 1, size(parse_errors)
                  call error_ring_push(trim(parse_errors(e)))
                  call logger%error(trim(parse_errors(e)))
               end do
            end if
            ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
            failed = .true.
         end if
      else
         if (src == SRC_TOML_FILE) then
            call toml_schema_parse_file(sch, text)
         else
            call toml_schema_parse_string(sch, text)
         end if
      end if
   end subroutine parse_into_schema

   module subroutine apply_time_unit_cascade(cfg)
      !! Post-parse time-unit fixup shared by the file + string config
      !! paths.  `t_end`, `status_interval`, `ocean_diag%dt_out` are given
      !! in `time_unit` from `[time]`; multiply through to SI seconds.
      !! `dt_fixed`/`dt_max` stay in seconds; default `time_unit="s"`
      !! (and it is reset to "s" once the conversion is done)
      !! gives a 1× factor (bit-identical to before).
      type(config_t), intent(inout) :: cfg

      real(wp) :: time_factor
      logical :: unit_ok

      unit_ok = .true.
      select case (trim(cfg%time_unit))
      case ("s", "sec", "second", "seconds")
         time_factor = 1.0_wp
      case ("min", "minute", "minutes")
         time_factor = 60.0_wp
      case ("hr", "hour", "hours", "h")
         time_factor = 3600.0_wp
      case ("day", "days", "d")
         time_factor = 86400.0_wp
      case ("year", "years", "yr")
         ! Julian year (365.25 d) — matches MOM6 / UDUNITS.
         time_factor = 365.25_wp*86400.0_wp
      case default
         call logger%error("Invalid time_unit = '"//trim(cfg%time_unit)// &
                           "': must be one of s/min/hr/day/year")
         time_factor = 1.0_wp  ! fallthrough so validation reports it
         unit_ok = .false.
      end select
      cfg%t_end = cfg%t_end*time_factor
      cfg%status_interval = cfg%status_interval*time_factor
      cfg%ocean%diag%dt_out = cfg%ocean%diag%dt_out*time_factor
      ! The three quantities are now in seconds, so the resolved configuration
      ! says so: `time_unit` reads "s".  Otherwise a doc dump (`write_doc_all`,
      ! `write_doc_short`) that states the converted seconds next to the original
      ! unit would be scaled a second time when it is read back.  Nothing else
      ! reads `time_unit`.
      if (unit_ok) cfg%time_unit = "s"
   end subroutine apply_time_unit_cascade

   module subroutine apply_cartesian_degrees(cfg)
      !! Post-parse Cartesian domain sizing (MOM6 GRID_CONFIG="cartesian" +
      !! AXIS_UNITS / LENLON / LENLAT).  When `[ocean.grid] len_lon`/`len_lat`
      !! are set (> 0) on an ocean Cartesian grid, DERIVE the uniform
      !! `[grid] dx`/`dy` (metres) from the domain extent, the physical cell
      !! count (`nx`/`ny` = MOM6 NIGLOBAL/NJGLOBAL) and `axis_units`, exactly as
      !! MOM6's `set_grid_metrics_cartesian`:
      !!   degrees : dx = rad_earth · len_lon · π/180 / nx  (arc length, no cos(lat))
      !!   km      : dx = 1000 · len_lon / nx
      !!   meters  : dx =        len_lon / nx
      !! and likewise dy from len_lat / ny.  This runs BEFORE grid init, so the
      !! metrics, barotropic CFL and Coriolis all see the derived metres.
      !! Default (len_lon <= 0) leaves `dx`/`dy` untouched ⇒ bit-identical.
      type(config_t), intent(inout) :: cfg

      real(wp) :: dx_scale, dy_scale

      ! Opt-in: only engage when the extent knobs are set.
      if (cfg%ocean%grid%len_lon <= 0.0_wp .and. cfg%ocean%grid%len_lat <= 0.0_wp) return

      if (trim(cfg%sim_type) /= "ocean") then
         call logger%error("len_lon/len_lat are [ocean.grid] knobs; "// &
                           "set sim_type='ocean' to use them")
         return
      end if
      if (trim(cfg%ocean%grid%grid_config) /= "cartesian") then
         call logger%error("len_lon/len_lat (MOM6 LENLON/LENLAT) apply to "// &
                           "grid_config='cartesian' only (got '"// &
                           trim(cfg%ocean%grid%grid_config)//"'); the spherical/"// &
                           "tripolar paths already read [grid] dx/dy as degrees")
         return
      end if
      if (cfg%ocean%grid%len_lon <= 0.0_wp .or. cfg%ocean%grid%len_lat <= 0.0_wp) then
         call logger%error("len_lon AND len_lat must both be > 0 to size the "// &
                           "Cartesian domain (got len_lon="// &
                           to_string(cfg%ocean%grid%len_lon)//", len_lat="// &
                           to_string(cfg%ocean%grid%len_lat)//")")
         return
      end if

      ! Per-unit metres-per-degree/km/metre scale (applied to len/N below).
      select case (trim(cfg%ocean%grid%axis_units))
      case ("meters", "meter", "m")
         dx_scale = 1.0_wp
         dy_scale = 1.0_wp
      case ("km", "kilometer", "kilometers")
         dx_scale = 1000.0_wp
         dy_scale = 1000.0_wp
      case ("degrees", "degree", "deg")
         ! Arc length on a sphere of radius rad_earth; NO cos(lat) — MOM6's
         ! Cartesian grid is a flat plane, unlike its spherical grid.
         dx_scale = cfg%ocean%grid%rad_earth*DEG2RAD
         dy_scale = cfg%ocean%grid%rad_earth*DEG2RAD
      case default
         call logger%error("Invalid axis_units = '"// &
                           trim(cfg%ocean%grid%axis_units)// &
                           "': must be meters/degrees/km")
         return
      end select

      cfg%dx = dx_scale*cfg%ocean%grid%len_lon/real(cfg%nx, wp)
      cfg%dy = dy_scale*cfg%ocean%grid%len_lat/real(cfg%ny, wp)

      call logger%info("Cartesian domain sized from extent (MOM6-style): len_lon="// &
                       to_string(cfg%ocean%grid%len_lon)//", len_lat="// &
                       to_string(cfg%ocean%grid%len_lat)//" "// &
                       trim(cfg%ocean%grid%axis_units)//" over "// &
                       to_string(cfg%nx)//"x"//to_string(cfg%ny)//" cells => dx="// &
                       to_string(cfg%dx)//" m, dy="//to_string(cfg%dy)//" m")
   end subroutine apply_cartesian_degrees

end submodule rdb_config_read
