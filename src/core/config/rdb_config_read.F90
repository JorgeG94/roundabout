submodule(rdb_config) rdb_config_read
   !! Namelist / string reading entry points + post-parse cascades.
   implicit none

contains

   module subroutine read_config(filename, cfg, schema, ierr)
      !! Read simulation configuration from a namelist file.
      !!
      !! Groups still read natively (`read(nml=)`) are seeded from the
      !! `cfg` defaults, read, and copied back as before.  Groups that
      !! have migrated to the strict schema (`rdb_nml_schema`) are
      !! validated and applied here, AFTER the native reads, by building
      !! a local schema (defaults captured from `cfg`'s pristine field
      !! initialisers) and strict-parsing the file — so `error stop` on
      !! any unknown key / range / enum violation fires for every caller
      !! (driver, FFI, unit tests), not just the production entry points.
      !!
      !! When `schema` is present the built schema is returned for the
      !! caller's parameter-doc dumps.  Its key pointers alias into
      !! `cfg`; both entry points keep their `cfg` alive (and `target`),
      !! so returning it is safe.
      character(len=*), intent(in) :: filename
         !! Path to the namelist input file
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
         !! Path to the namelist input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out) :: sch
         !! Built + parsed schema (aliases into `cfg`)
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`).

      ! `cfg` arrives at its pristine field-initialiser defaults.  Every
      ! namelist group is now on the strict schema: build it (capturing
      ! those defaults), then parse — which validates AND applies every
      ! override straight into `cfg` via the schema's key pointers.  No
      ! native readers, no flat shadow locals, no copy-back remain.

      ! Build the strict schema so its key constructors capture defaults
      ! from the pristine `cfg` field initialisers.
      call build_rdb_schema(cfg, sch)

      ! A missing namelist file is not an error — the pristine `cfg`
      ! defaults stand (production main.F90 fails loudly on a missing
      ! file earlier; tests + the FFI rely on this defaults-on-absent
      ! behaviour).  Guard the strict parse so it never error-stops on
      ! "cannot open file".
      block
         logical :: file_exists
         inquire (file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            call logger%debug("Namelist file not found, keeping defaults: "//trim(filename))
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if
      end block

      ! Strict schema parse: validate + APPLY every group straight into
      ! `cfg`.  When `ierr` is present, route through the schema's own
      ! non-aborting `status`/`errors` pair (already implemented in
      ! `rdb_nml_schema`) instead of the bare (strict, error-stopping)
      ! form — every collected message is still logged, just without the
      ! `error stop`.  Absent `ierr` ⇒ unchanged bare call, for EVERY
      ! existing caller (driver, benchmarks, unit tests).
      if (present(ierr)) then
         block
            integer :: parse_status, e
            character(len=:), allocatable :: parse_errors(:)
            call sch%parse(trim(filename), status=parse_status, errors=parse_errors)
            if (parse_status /= 0) then
               if (allocated(parse_errors)) then
                  do e = 1, size(parse_errors)
                     call error_ring_push(trim(parse_errors(e)))
                     call logger%error(trim(parse_errors(e)))
                  end do
               end if
               ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
               return
            end if
         end block
      else
         call sch%parse(trim(filename))
      end if

      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK

   end subroutine read_config_impl

   module subroutine read_config_from_string(text, cfg, schema, ierr)
      !! In-memory sibling of [[read_config]]: build + strict-parse +
      !! apply a namelist held entirely in the `text` buffer (newline-
      !! separated), with NO filesystem touch — no temp file, no `chdir`,
      !! no cleanup.  Same validation behaviour as the file path (the
      !! schema parser and the native &ocean_bc_nml read both run from
      !! the in-memory line array); same `error stop` default and same
      !! optional non-aborting `ierr` as [[read_config]].
      character(len=*), intent(in) :: text
         !! Whole namelist as one string (records separated by '\n').
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out), optional :: schema
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`). The natural FFI/create()
         !! door: a bad namelist string becomes a returned status
         !! instead of aborting the host process.

      if (present(schema)) then
         call read_config_from_string_impl(text, cfg, schema, ierr=ierr)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_from_string_impl(text, cfg, sch, ierr=ierr)
         end block
      end if
   end subroutine read_config_from_string

   module subroutine read_config_from_string_impl(text, cfg, sch, ierr)
      !! Body of [[read_config_from_string]] — mirrors read_config_impl
      !! with the line array replacing the file unit.
      character(len=*), intent(in) :: text
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out) :: sch
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`).

      character(len=:), allocatable :: lines(:)
      integer :: n_lines

      call build_rdb_schema(cfg, sch)

      ! Empty buffer → keep pristine defaults (mirrors the missing-file
      ! branch of read_config_impl).
      if (len_trim(text) == 0) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      call split_to_lines(text, lines, n_lines)

      ! Strict schema parse from the same lines: validate + APPLY every
      ! migrated group straight into `cfg`.  Same status/errors routing
      ! as read_config_impl — see its comment for the rationale.
      if (present(ierr)) then
         block
            integer :: parse_status, e
            character(len=:), allocatable :: parse_errors(:)
            call sch%parse_lines(lines, n_lines, status=parse_status, errors=parse_errors)
            if (parse_status /= 0) then
               if (allocated(parse_errors)) then
                  do e = 1, size(parse_errors)
                     call error_ring_push(trim(parse_errors(e)))
                     call logger%error(trim(parse_errors(e)))
                  end do
               end if
               ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
               return
            end if
         end block
      else
         call sch%parse_lines(lines, n_lines)
      end if

      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine read_config_from_string_impl

   module subroutine split_to_lines(text, lines, n_lines)
      !! Split a newline-separated buffer into a fixed-len character array
      !! (one record per line, trailing CR stripped) for internal-file
      !! reads + the schema line walker.
      character(len=*), intent(in) :: text
      character(len=:), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines

      integer :: i, start, maxlen, cur, nl
      character(len=1) :: c

      ! First pass: count records + longest line.
      nl = 1
      maxlen = 0
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            nl = nl + 1
            maxlen = max(maxlen, cur)
            cur = 0
         else if (c /= achar(13)) then
            cur = cur + 1
         end if
      end do
      maxlen = max(maxlen, cur, 1)

      allocate (character(len=maxlen) :: lines(nl))
      ! Blank-fill via a SECTION, not `lines = ""`.  `lines` has a DEFERRED
      ! length, so whole-variable intrinsic assignment from a zero-length
      ! scalar re-allocates it with len = 0 (F2018 10.2.1.3p3) — ifx does
      ! exactly that, and every record then parses as empty.  A section
      ! designator is not an allocatable variable, so it blank-pads in place.
      lines(:) = ""

      ! Second pass: copy each record.
      n_lines = 0
      start = 1
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            n_lines = n_lines + 1
            if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
            start = i + 1
            cur = 0
         else if (c /= achar(13)) then
            ! CR is skipped (a trailing CR before LF is excluded because
            ! `cur` counts only content chars and `start` marks the first).
            if (cur == 0) start = i
            cur = cur + 1
         end if
      end do
      ! Trailing record (no final newline).
      n_lines = n_lines + 1
      if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
   end subroutine split_to_lines

   module subroutine apply_time_unit_cascade(cfg)
      !! Post-parse time-unit fixup shared by the file + string config
      !! paths.  `t_end`, `status_interval`, `ocean_diag%dt_out` are given
      !! in `time_unit` from `&time_nml`; multiply through to SI seconds.
      !! `dt_fixed`/`dt_max` stay in seconds; default `time_unit="s"`
      !! gives a 1× factor (bit-identical to before).
      type(config_t), intent(inout) :: cfg

      real(wp) :: time_factor

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
      end select
      cfg%t_end = cfg%t_end*time_factor
      cfg%status_interval = cfg%status_interval*time_factor
      cfg%ocean%diag%dt_out = cfg%ocean%diag%dt_out*time_factor
   end subroutine apply_time_unit_cascade

   module subroutine apply_cartesian_degrees(cfg)
      !! Post-parse Cartesian domain sizing (MOM6 GRID_CONFIG="cartesian" +
      !! AXIS_UNITS / LENLON / LENLAT).  When `&ocean_grid_nml len_lon`/`len_lat`
      !! are set (> 0) on an ocean Cartesian grid, DERIVE the uniform
      !! `&grid_nml dx`/`dy` (metres) from the domain extent, the physical cell
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
         call logger%error("len_lon/len_lat are &ocean_grid_nml knobs; "// &
                           "set sim_type='ocean' to use them")
         return
      end if
      if (trim(cfg%ocean%grid%grid_config) /= "cartesian") then
         call logger%error("len_lon/len_lat (MOM6 LENLON/LENLAT) apply to "// &
                           "grid_config='cartesian' only (got '"// &
                           trim(cfg%ocean%grid%grid_config)//"'); the spherical/"// &
                           "tripolar paths already read &grid_nml dx/dy as degrees")
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
