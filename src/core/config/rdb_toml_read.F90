!! TOML front end for the strict configuration schema.
module rdb_toml_read
   !! Reads a TOML document (toml-f) and applies it to an `nml_schema_t` through
   !! the schema's per-key validators: one set of range / enum / length /
   !! required-key / cross-check rules and one error list (`path:line:`
   !! prefixed, accumulating).
   !!
   !! Layout (the group <-> table mapping lives in [[group_to_table]] /
   !! [[table_to_group]] and nowhere else):
   !!
   !!   group `ocean_<x>`  <->  `[ocean.<x>]`   (only the LEADING `ocean_` splits)
   !!   group `<x>`        <->  `[<x>]`
   !!
   !! Key names are EXACT-CASE: a key or table spelled in a different case is
   !! unknown (and the did-you-mean points at the registered spelling).  Value
   !! types are TOML's own: reals are floats or integers, integers are
   !! integers, logicals are booleans, strings and enums are strings, real
   !! arrays are arrays of numbers (a lone number is accepted for an array key).
   !!
   !! Values reach the key validators as TEXT (`parse_tokens`): a number as its
   !! SOURCE spelling (underscores and a leading `+` dropped; an integer for a
   !! real key gets `.0` appended), which list-directed `read` returns bit for
   !! bit -- `tests/test_toml_config.F90` proves it over a sweep of bit patterns
   !! -- so stringifying loses nothing.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use rdb_nml_schema, only: nml_schema_t, nml_group_t, nml_real_key_t, &
                             nml_int_key_t, nml_logical_key_t, nml_string_key_t, &
                             nml_enum_key_t, nml_real_array_key_t, nml_key_t, add_err, &
                             suggest, group_key_list, make_prefix, reset_state, &
                             post_parse_validate
   use tomlf, only: toml_table, toml_array, toml_keyval, toml_value, toml_key, &
                    toml_context, toml_error, toml_load, toml_loads, get_value, len
   use tomlf_constants, only: toml_type
   implicit none
   private

   public :: group_to_table, table_to_group
   public :: toml_schema_parse_file, toml_schema_parse_string

   integer, parameter :: WANT_NUMBER = 1
   integer, parameter :: WANT_INTEGER = 2
   integer, parameter :: WANT_BOOLEAN = 3
   integer, parameter :: WANT_STRING = 4

   character(len=*), parameter :: OCEAN_PREFIX = "ocean_"
   character(len=*), parameter :: OCEAN_TABLE = "ocean"

contains

   ! ------------------------------------------------------------------
   ! The group <-> table mapping (the only place it is spelled out).
   ! ------------------------------------------------------------------

   pure function group_to_table(group) result(table)
      !! Table path of a registered group name (the bare name, as stored in
      !! the schema: `ocean_hvisc`, `grid`).  `ocean_<x>` -> `ocean.<x>`;
      !! anything else is its own top-level table.
      character(len=*), intent(in) :: group
      character(len=:), allocatable :: table
      integer :: n

      n = len(OCEAN_PREFIX)
      if (len_trim(group) > n) then
         if (group(1:n) == OCEAN_PREFIX) then
            table = OCEAN_TABLE//"."//trim(group(n + 1:))
            return
         end if
      end if
      table = trim(group)
   end function group_to_table

   pure function table_to_group(table) result(group)
      !! Group name of a table path, or an EMPTY string when `table` is not the
      !! canonical path of any group: the bare container `ocean`, a path with a
      !! dot other than the `ocean.` one, and a top-level `ocean_<x>` (the
      !! group's one spelling is `ocean.<x>`).  Inverse of [[group_to_table]]
      !! on every name it accepts.
      character(len=*), intent(in) :: table
      character(len=:), allocatable :: group
      character(len=:), allocatable :: t
      integer :: n

      t = trim(table)
      group = ""
      n = len(OCEAN_TABLE) + 1
      if (len(t) == 0) return
      if (len(t) > n) then
         if (t(1:n) == OCEAN_TABLE//".") then
            if (index(t(n + 1:), ".") == 0) group = OCEAN_PREFIX//t(n + 1:)
            return
         end if
      end if
      if (index(t, ".") > 0) return
      if (t == OCEAN_TABLE) return
      if (len(t) >= len(OCEAN_PREFIX)) then
         if (t(1:len(OCEAN_PREFIX)) == OCEAN_PREFIX) return
      end if
      group = t
   end function table_to_group

   ! ------------------------------------------------------------------
   ! Entry points (same error-return convention as `nml_schema_t%parse`).
   ! ------------------------------------------------------------------

   subroutine toml_schema_parse_file(schema, path, status, errors)
      !! Parse + validate + apply the TOML file at `path`.  If BOTH `status`
      !! and `errors` are present the collected messages return in them and
      !! nothing stops; otherwise each is logged and the routine `error stop`s.
      type(nml_schema_t), intent(inout) :: schema
      character(len=*), intent(in) :: path
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)

      type(toml_table), allocatable :: root
      type(toml_error), allocatable :: terr
      type(toml_context) :: ctx
      character(len=:), allocatable :: err_buf(:)
      integer :: n_err
      logical :: syntax_ok

      n_err = 0
      call reset_state(schema)
      call toml_load(root, path, context=ctx, error=terr)
      syntax_ok = .not. allocated(terr)
      if (syntax_ok) then
         call apply_root(schema, root, ctx, path, err_buf, n_err)
      else
         call add_err(err_buf, n_err, syntax_message(terr, path))
      end if
      call post_parse_validate(schema, path, err_buf, n_err, status, errors, &
                               skip_checks=.not. syntax_ok)
   end subroutine toml_schema_parse_file

   subroutine toml_schema_parse_string(schema, text, status, errors)
      !! In-memory sibling of [[toml_schema_parse_file]]; errors are labelled
      !! `<config-string>`.
      type(nml_schema_t), intent(inout) :: schema
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)

      character(len=*), parameter :: LABEL = "<config-string>"
      type(toml_table), allocatable :: root
      type(toml_error), allocatable :: terr
      type(toml_context) :: ctx
      character(len=:), allocatable :: err_buf(:)
      integer :: n_err
      logical :: syntax_ok

      n_err = 0
      call reset_state(schema)
      call toml_loads(root, text, context=ctx, error=terr)
      syntax_ok = .not. allocated(terr)
      if (syntax_ok) then
         call apply_root(schema, root, ctx, LABEL, err_buf, n_err)
      else
         call add_err(err_buf, n_err, syntax_message(terr, LABEL))
      end if
      call post_parse_validate(schema, LABEL, err_buf, n_err, status, errors, &
                               skip_checks=.not. syntax_ok)
   end subroutine toml_schema_parse_string

   function syntax_message(terr, label) result(msg)
      !! toml-f's first syntax error, tagged with the source label.
      type(toml_error), intent(in) :: terr
      character(len=*), intent(in) :: label
      character(len=:), allocatable :: msg
      if (index(terr%message, "Could not open") > 0) then
         ! toml-f's lexer could not read the source (a directory, a file the
         ! process may not read): not a syntax problem.
         msg = label//": cannot open file ("//trim(terr%message)//")"
      else
         msg = label//": TOML syntax error: "//terr%message
      end if
   end function syntax_message

   ! ------------------------------------------------------------------
   ! The walk.
   ! ------------------------------------------------------------------

   pure function line_of(ctx, origin) result(line)
      !! 1-based source line of a toml-f origin token (0 when unknown).
      type(toml_context), intent(in) :: ctx
      integer, intent(in) :: origin
      integer :: line, i, first

      line = 0
      if (origin <= 0 .or. origin > ctx%top) return
      if (.not. allocated(ctx%source)) return
      first = min(ctx%token(origin)%first, len(ctx%source) + 1)
      line = 1
      do i = 1, first - 1
         if (ctx%source(i:i) == achar(10)) line = line + 1
      end do
   end function line_of

   pure function loc(ctx, label, origin) result(prefix)
      !! `label:line: ` for an origin token.
      type(toml_context), intent(in) :: ctx
      character(len=*), intent(in) :: label
      integer, intent(in) :: origin
      character(len=:), allocatable :: prefix
      prefix = make_prefix(label, max(line_of(ctx, origin), 1))
   end function loc

   subroutine apply_root(schema, root, ctx, label, errors, n_err)
      !! Walk the top level: `[ocean]` is the container of the `ocean_*`
      !! groups, every other table is a group, anything else is refused.
      type(nml_schema_t), intent(inout) :: schema
      type(toml_table), intent(inout) :: root
      type(toml_context), intent(in) :: ctx
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      type(toml_key), allocatable :: keys(:)
      class(toml_value), pointer :: v
      integer :: i

      call root%get_keys(keys)
      do i = 1, size(keys)
         call root%get(keys(i)%key, v)
         select type (v)
         type is (toml_table)
            if (index(keys(i)%key, ".") > 0) then
               call add_err(errors, n_err, loc(ctx, label, keys(i)%origin)//"table name '"// &
                            keys(i)%key//"' contains a dot; write nested tables as "// &
                            "[ocean.<group>], not as a quoted name")
            else if (keys(i)%key == OCEAN_TABLE) then
               call apply_ocean_container(schema, v, ctx, label, errors, n_err)
            else
               call apply_group(schema, v, keys(i)%key, keys(i)%origin, ctx, label, errors, n_err)
            end if
         class default
            call add_err(errors, n_err, loc(ctx, label, keys(i)%origin)//"key '"// &
                         keys(i)%key//"' is outside any table; keys belong inside a "// &
                         "[group] table (see docs: [ocean.<group>] or [<group>])")
         end select
      end do
   end subroutine apply_root

   subroutine apply_ocean_container(schema, tbl, ctx, label, errors, n_err)
      !! `[ocean]` holds only sub-tables: `[ocean.<x>]` is group `ocean_<x>`.
      type(nml_schema_t), intent(inout) :: schema
      type(toml_table), intent(inout) :: tbl
      type(toml_context), intent(in) :: ctx
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      type(toml_key), allocatable :: keys(:)
      class(toml_value), pointer :: v
      integer :: i

      call tbl%get_keys(keys)
      do i = 1, size(keys)
         call tbl%get(keys(i)%key, v)
         select type (v)
         type is (toml_table)
            call apply_group(schema, v, OCEAN_TABLE//"."//keys(i)%key, keys(i)%origin, &
                             ctx, label, errors, n_err)
         class default
            call add_err(errors, n_err, loc(ctx, label, keys(i)%origin)//"key '"// &
                         keys(i)%key//"' in [ocean] is not a table; ocean settings "// &
                         "belong in [ocean.<group>] tables")
         end select
      end do
   end subroutine apply_ocean_container

   subroutine apply_group(schema, tbl, tpath, torigin, ctx, label, errors, n_err)
      !! One group table: resolve it EXACT-case, then every key through the
      !! schema's own validators.
      type(nml_schema_t), intent(inout) :: schema
      type(toml_table), intent(inout) :: tbl
      character(len=*), intent(in) :: tpath
      integer, intent(in) :: torigin
      type(toml_context), intent(in) :: ctx
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      type(toml_key), allocatable :: keys(:)
      class(toml_value), pointer :: v
      character(len=:), allocatable :: gname
      character(len=:), allocatable :: names(:)
      integer :: gidx, i

      gname = table_to_group(tpath)
      gidx = 0
      if (len(gname) > 0) then
         gidx = schema%find_group(gname)
      end if
      if (gidx == 0) then
         names = table_name_list(schema)
         call add_err(errors, n_err, loc(ctx, label, torigin)//"unknown table ["// &
                      tpath//"]"//suggest(tpath, names))
         return
      end if

      if (schema%groups(gidx)%found) then
         call add_err(errors, n_err, loc(ctx, label, torigin)//"duplicate occurrence of group '"// &
                      schema%groups(gidx)%name//"' (table ["//tpath//"])")
         return
      end if
      schema%groups(gidx)%found = .true.
      call tbl%get_keys(keys)
      do i = 1, size(keys)
         call tbl%get(keys(i)%key, v)
         call apply_key(schema%groups(gidx), tpath, keys(i), v, ctx, label, errors, n_err)
      end do
   end subroutine apply_group

   function table_name_list(schema) result(names)
      !! The canonical table path of every registered group.
      type(nml_schema_t), intent(in) :: schema
      character(len=:), allocatable :: names(:)
      integer :: i, mlen

      mlen = 1
      do i = 1, size(schema%groups)
         mlen = max(mlen, len(group_to_table(schema%groups(i)%name)))
      end do
      allocate (character(len=mlen) :: names(size(schema%groups)))
      do i = 1, size(schema%groups)
         names(i) = group_to_table(schema%groups(i)%name)
      end do
   end function table_name_list

   subroutine apply_key(g, tpath, key, v, ctx, label, errors, n_err)
      !! One `key = value` of group `g`: resolve (exact-case), type-check,
      !! stringify and hand to the key's `parse_tokens`.
      type(nml_group_t), intent(inout) :: g
      character(len=*), intent(in) :: tpath
      type(toml_key), intent(in) :: key
      class(toml_value), pointer, intent(in) :: v
      type(toml_context), intent(in) :: ctx
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      character(len=:), allocatable :: prefix, emsg, why
      character(len=:), allocatable :: vals(:), names(:)
      integer :: kidx, want
      logical :: array_ok, ok

      prefix = loc(ctx, label, key%origin)

      select type (v)
      type is (toml_table)
         call add_err(errors, n_err, prefix//"unexpected sub-table '"//key%key//"' in ["// &
                      tpath//"]; a group table holds keys only")
         return
      end select

      kidx = g%find_key(key%key)
      if (kidx > 0) then
         if (g%keys(kidx)%key%name /= key%key) kidx = 0
      end if
      if (kidx == 0) then
         names = group_key_list(g)
         call add_err(errors, n_err, prefix//"unknown key '"//key%key//"' in table ["// &
                      tpath//"]"//suggest(key%key, names))
         return
      end if
      g%keys(kidx)%key%found = .true.

      call key_expectation(g%keys(kidx)%key, want, array_ok)
      call to_values(v, ctx, want, array_ok, vals, ok, why)
      if (.not. ok) then
         call add_err(errors, n_err, prefix//"key '"//key%key//"' in ["//tpath//"] "//why)
         return
      end if
      call g%keys(kidx)%key%parse_tokens(vals, emsg)
      if (allocated(emsg)) call add_err(errors, n_err, prefix//trim(emsg))
   end subroutine apply_key

   pure subroutine key_expectation(k, want, array_ok)
      !! What TOML type a registered key takes.
      class(nml_key_t), intent(in) :: k
      integer, intent(out) :: want
      logical, intent(out) :: array_ok

      array_ok = .false.
      select type (k)
      type is (nml_real_array_key_t)
         want = WANT_NUMBER
         array_ok = .true.
      type is (nml_real_key_t)
         want = WANT_NUMBER
      type is (nml_int_key_t)
         want = WANT_INTEGER
      type is (nml_logical_key_t)
         want = WANT_BOOLEAN
      type is (nml_enum_key_t)
         want = WANT_STRING
      type is (nml_string_key_t)
         want = WANT_STRING
      class default
         want = WANT_STRING
      end select
   end subroutine key_expectation

   subroutine to_values(v, ctx, want, array_ok, vals, ok, why)
      !! TOML value -> the text tokens `parse_tokens` takes; `why` is the
      !! refusal (ending the sentence "key 'k' in [t] ...") when not `ok`.
      class(toml_value), pointer, intent(in) :: v
      type(toml_context), intent(in) :: ctx
      integer, intent(in) :: want
      logical, intent(in) :: array_ok
      character(len=:), allocatable, intent(out) :: vals(:)
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: why

      character(len=:), allocatable :: one
      character(len=:), allocatable :: tmp(:)
      class(toml_value), pointer :: ev
      integer :: n, j, mlen

      ok = .false.
      select type (v)
      type is (toml_keyval)
         call scalar_text(v, ctx, want, one, why)
         if (.not. allocated(one)) return
         allocate (character(len=max(len(one), 1)) :: vals(1))
         vals(1) = one
         ok = .true.
      type is (toml_array)
         if (.not. array_ok) then
            why = "expects "//want_name(want)//", got an array"
            return
         end if
         n = len(v)
         if (n < 1) then
            why = "has no value (an empty array)"
            return
         end if
         mlen = 1
         allocate (character(len=512) :: tmp(n))
         do j = 1, n
            call v%get(j, ev)
            select type (ev)
            type is (toml_keyval)
               call scalar_text(ev, ctx, want, one, why)
               if (.not. allocated(one)) then
                  why = "element "//itoa(j)//": "//why
                  return
               end if
               tmp(j) = one
               mlen = max(mlen, len(one))
            class default
               why = "element "//itoa(j)//" must be a number, not a nested array or table"
               return
            end select
         end do
         allocate (character(len=mlen) :: vals(n))
         do j = 1, n
            vals(j) = tmp(j)
         end do
         ok = .true.
      class default
         why = "has an unsupported value"
      end select
   end subroutine to_values

   subroutine scalar_text(kv, ctx, want, text, why)
      !! One TOML scalar as text, or `why` when its type is not `want`.
      !! NUMBERS are handed over as their SOURCE spelling (the value token
      !! in `ctx`), never re-formatted from toml-f's parsed value: toml-f's
      !! integer extraction wraps silently on overflow, so only the raw digits
      !! let the schema's own `read` see the real magnitude (and refuse it).
      type(toml_keyval), intent(inout) :: kv
      type(toml_context), intent(in) :: ctx
      integer, intent(in) :: want
      character(len=:), allocatable, intent(out) :: text
      character(len=:), allocatable, intent(out) :: why

      logical :: lval
      character(len=:), allocatable :: sval
      integer :: t, stat, i

      t = kv%get_type()
      select case (want)
      case (WANT_NUMBER, WANT_INTEGER)
         if (t == toml_type%int .or. (t == toml_type%float .and. want == WANT_NUMBER)) then
            call raw_number(ctx, kv%origin_value, t == toml_type%int, text, why)
            ! An integer literal given for a REAL key reaches the real parser as a
            ! real literal (`.0` appended): list-directed `read` of a bare integer
            ! longer than int64 into a real is refused by some compilers (nvfortran)
            ! and read as 1e20 by others (gfortran); the real spelling reads the same
            ! value everywhere, and is bit-identical to the bare spelling wherever
            ! both read.
            if (want == WANT_NUMBER .and. t == toml_type%int .and. allocated(text)) text = text//".0"
            return
         end if
      case (WANT_BOOLEAN)
         if (t == toml_type%boolean) then
            call get_value(kv, lval, stat=stat)
            if (lval) then
               text = ".true."
            else
               text = ".false."
            end if
            return
         end if
      case (WANT_STRING)
         if (t == toml_type%string) then
            call get_value(kv, sval, stat=stat)
            if (.not. allocated(sval)) sval = ""
            do i = 1, len(sval)
               if (iachar(sval(i:i)) < 32 .or. iachar(sval(i:i)) == 127) then
                  why = "is a string with a control character (newline, tab, ...) at "// &
                        "position "//itoa(i)//"; not allowed in a configuration value"
                  return
               end if
            end do
            text = sval
            return
         end if
      case default
         continue
      end select
      why = "expects "//want_name(want)//", got "//type_name(t)
   end subroutine scalar_text

   subroutine raw_number(ctx, origin, is_int, text, why)
      !! The source spelling of a numeric value token, normalised to what the
      !! schema's list-directed `read` takes: underscores dropped, a leading
      !! `+` dropped, `inf`/`nan` lower-cased.  Hexadecimal / octal / binary
      !! integers are refused (write them in decimal).
      type(toml_context), intent(in) :: ctx
      integer, intent(in) :: origin
      logical, intent(in) :: is_int
      character(len=:), allocatable, intent(out) :: text
      character(len=:), allocatable, intent(out) :: why

      character(len=:), allocatable :: raw
      integer :: i, f, l

      if (origin <= 0 .or. origin > ctx%top .or. .not. allocated(ctx%source)) then
         why = "has a numeric value whose source text cannot be recovered"
         return
      end if
      f = ctx%token(origin)%first
      l = ctx%token(origin)%last
      if (f < 1 .or. l > len(ctx%source) .or. l < f) then
         why = "has a numeric value whose source text cannot be recovered"
         return
      end if
      raw = ""
      do i = f, l
         if (ctx%source(i:i) /= "_") raw = raw//ctx%source(i:i)
      end do
      if (len(raw) > 0) then
         if (raw(1:1) == "+") raw = raw(2:)
      end if
      if (len(raw) == 0) then
         why = "has an empty numeric value"
         return
      end if
      if (is_int .and. len(raw) > 1) then
         if (raw(1:1) == "0" .and. scan(raw(2:2), "xXoObB") > 0) then
            why = "is a hexadecimal/octal/binary integer; write the value in decimal"
            return
         end if
      end if
      if (.not. is_int) then
         if (scan(raw, "iInN") > 0) raw = lowercase(raw)
      end if
      text = raw
   end subroutine raw_number

   pure function lowercase(s) result(o)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: o
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            o(i:i) = achar(c + 32)
         else
            o(i:i) = s(i:i)
         end if
      end do
   end function lowercase

   pure function want_name(want) result(s)
      integer, intent(in) :: want
      character(len=:), allocatable :: s
      select case (want)
      case (WANT_NUMBER)
         s = "a number"
      case (WANT_INTEGER)
         s = "an integer"
      case (WANT_BOOLEAN)
         s = "a boolean (true/false)"
      case default
         s = "a string"
      end select
   end function want_name

   pure function type_name(t) result(s)
      integer, intent(in) :: t
      character(len=:), allocatable :: s
      if (t == toml_type%int) then
         s = "an integer"
      else if (t == toml_type%float) then
         s = "a float"
      else if (t == toml_type%boolean) then
         s = "a boolean"
      else if (t == toml_type%string) then
         s = "a string"
      else if (t == toml_type%datetime) then
         s = "a datetime"
      else
         s = "an unsupported value"
      end if
   end function type_name

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write (buf, "(I0)") i
      s = trim(buf)
   end function itoa

end module rdb_toml_read
