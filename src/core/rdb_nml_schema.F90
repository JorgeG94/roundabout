!! Strict, schema-driven configuration validation engine.
module rdb_nml_schema
   !! A host-only schema of configuration groups and keys with the per-key
   !! validators (range / enum / length / required / cross-check) the TOML
   !! front end (`rdb_toml_read`) applies a parsed document through,
   !! collecting ALL errors (with `path:line:` prefixes and did-you-mean
   !! suggestions).  The module keeps its historical `nml_*` names; a group
   !! `ocean_<x>` is the TOML table `[ocean.<x>]`, any other group `<x>` the
   !! table `[<x>]`.
   !!
   !! Defaults contract (load-bearing): every key constructor captures its
   !! default from the CURRENT target value at registration; the schema
   !! NEVER takes a default argument. Register AFTER config defaults exist,
   !! BEFORE parsing.
   !!
   !! Required-key semantics: a `required=.true.` key must appear whenever
   !! its TABLE appears; an absent group is fine (defaults rule).
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use pic_logger, only: global_logger
   use rdb_error_ring, only: error_ring_push
   implicit none
   private

   integer, parameter :: rk = wp
      !! Working real kind — tracks the solver `wp` so `real(wp)` config
      !! targets associate in either single- or double-precision builds.

   integer, parameter :: MAX_SUGGEST_LEN = 64
      !! Cap on names fed to the Levenshtein suggester.

   public :: nml_key_t
   public :: nml_real_key_t, nml_int_key_t, nml_logical_key_t
   public :: nml_string_key_t, nml_enum_key_t, nml_real_array_key_t
   public :: nml_group_t, nml_schema_t, nml_key_box_t
   public :: nml_real, nml_int, nml_logical, nml_string, nml_enum, nml_real_array
   public :: group_check_iface
   ! Building blocks for the front end (rdb_toml_read) that applies a parsed
   ! document through the key validators and the error list.
   public :: add_err, suggest, group_key_list, make_prefix
   public :: reset_state, post_parse_validate

   type, abstract :: nml_key_t
      !! Abstract base for a single namelist key descriptor.
      character(len=:), allocatable :: name
         !! Key name (matched case-insensitively).
      character(len=:), allocatable :: doc
         !! Human-readable description (FORD-style synthesis doc).
      character(len=:), allocatable :: units
         !! Optional units string, default ''.
      logical :: required = .false.
         !! If true, the key must appear whenever its group appears.
      logical :: found = .false.
         !! Parse-state: set true when the key is seen in the file.
      character(len=:), allocatable :: dead_reason
         !! Unallocated/empty = live. Non-empty = this knob is registered
         !! (accepted, type/range-validated) but does NOT reach any
         !! ocean-path behaviour -- e.g. `&vcoord_nml zstar_stretching`
         !! is captured on `cfg` but never copied onto `ocean_state%vcoord`
         !! (rdb_ocean_vcoord.F90:179,182). Carried into the JSON emitter
         !! (`render_json`) so a generated caller-side layer can warn on
         !! assignment instead of presenting a dead knob as live.
   contains
      procedure(parse_tokens_iface), deferred :: parse_tokens
      procedure(value_string_iface), deferred :: value_string
      procedure(default_string_iface), deferred :: default_string
      procedure(is_default_iface), deferred :: is_default
   end type nml_key_t

   abstract interface
      subroutine parse_tokens_iface(this, values, err_msg)
         !! Validate `values` and assign through the target pointer.
         !! On success `err_msg` is left unallocated.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(inout) :: this
         character(len=*), intent(in) :: values(:)
         character(len=:), allocatable, intent(out) :: err_msg
      end subroutine parse_tokens_iface

      function value_string_iface(this) result(str)
         !! Current target value, formatted for namelist output.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         character(len=:), allocatable :: str
      end function value_string_iface

      function default_string_iface(this) result(str)
         !! Captured default value, formatted for namelist output.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         character(len=:), allocatable :: str
      end function default_string_iface

      function is_default_iface(this) result(yes)
         !! True if the current target value equals the captured default.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         logical :: yes
      end function is_default_iface
   end interface

   abstract interface
      function group_check_iface() result(msg)
         !! Cross-key validation callback.  Empty string = OK; a
         !! non-empty string is appended to the collected errors.
         implicit none
         character(len=:), allocatable :: msg
      end function group_check_iface
   end interface

   type :: nml_key_box_t
      !! Boxes a polymorphic key so groups can hold a growable array.
      class(nml_key_t), allocatable :: key
   end type nml_key_box_t

   type, extends(nml_key_t) :: nml_real_key_t
      !! Real-valued key with optional [min,max] range validation.
      real(rk), pointer :: tgt => null()
      real(rk) :: default = 0.0_rk
      real(rk) :: vmin = 0.0_rk, vmax = 0.0_rk
      logical :: has_min = .false., has_max = .false.
   contains
      procedure :: parse_tokens => real_parse
      procedure :: value_string => real_value_string
      procedure :: default_string => real_default_string
      procedure :: is_default => real_is_default
   end type nml_real_key_t

   type, extends(nml_key_t) :: nml_int_key_t
      !! Integer-valued key with optional [min,max] range validation.
      integer, pointer :: tgt => null()
      integer :: default = 0
      integer :: vmin = 0, vmax = 0
      logical :: has_min = .false., has_max = .false.
   contains
      procedure :: parse_tokens => int_parse
      procedure :: value_string => int_value_string
      procedure :: default_string => int_default_string
      procedure :: is_default => int_is_default
   end type nml_int_key_t

   type, extends(nml_key_t) :: nml_logical_key_t
      !! Logical-valued key.
      logical, pointer :: tgt => null()
      logical :: default = .false.
   contains
      procedure :: parse_tokens => logical_parse
      procedure :: value_string => logical_value_string
      procedure :: default_string => logical_default_string
      procedure :: is_default => logical_is_default
   end type nml_logical_key_t

   type, extends(nml_key_t) :: nml_string_key_t
      !! String key.  Deferred-len char pointer associates with a
      !! fixed-len target and assumes its length; parsing errors if the
      !! parsed value is LONGER than len(tgt) (no silent truncation).
      character(len=:), pointer :: tgt => null()
      character(len=:), allocatable :: default
   contains
      procedure :: parse_tokens => string_parse
      procedure :: value_string => string_value_string
      procedure :: default_string => string_default_string
      procedure :: is_default => string_is_default
   end type nml_string_key_t

   type, extends(nml_string_key_t) :: nml_enum_key_t
      !! Enum key: a string restricted to an allowed list, validated
      !! case-insensitively and assigned its canonical spelling.
      character(len=:), allocatable :: allowed(:)
      character(len=:), allocatable :: retired(:)
         !! Spellings that USED to be legal and no longer are (renamed or
         !! withdrawn).  A namelist still carrying one must not fall
         !! through to the bare "not in allowed set {…}" message — that
         !! tells the reader the value is wrong but not what replaced it.
         !! Matching one of these fails loud with `retired_hint` instead.
      character(len=:), allocatable :: retired_hint
         !! What to say when a retired spelling is used.  Name the
         !! replacement; this string IS the migration instruction.
   contains
      procedure :: parse_tokens => enum_parse
   end type nml_enum_key_t

   type, extends(nml_key_t) :: nml_real_array_key_t
      !! Real-array key: accepts 1..size(tgt) values, fills from
      !! element 1, errors if more than size(tgt) values are given.
      real(rk), pointer :: tgt(:) => null()
      real(rk), allocatable :: default(:)
   contains
      procedure :: parse_tokens => real_array_parse
      procedure :: value_string => real_array_value_string
      procedure :: default_string => real_array_default_string
      procedure :: is_default => real_array_is_default
   end type nml_real_array_key_t

   type :: nml_group_t
      !! A namelist group: a name, doc, growable key list, parse-state,
      !! and an optional cross-key validation callback.
      character(len=:), allocatable :: name
         !! Group name (matched case-insensitively).
      character(len=:), allocatable :: doc
         !! Group description.
      type(nml_key_box_t), allocatable :: keys(:)
         !! Registered keys.
      logical :: found = .false.
         !! Parse-state: true when the group appears in the file.
      procedure(group_check_iface), pointer, nopass :: cross_check => null()
         !! Optional cross-key validator (empty string = OK).
   contains
      procedure :: add => group_add
      procedure :: find_key => group_find_key
   end type nml_group_t

   type :: nml_schema_t
      !! The top-level schema: the validated groups.
      type(nml_group_t), allocatable :: groups(:)
         !! Registered, validated groups.
   contains
      procedure :: add_group => schema_add_group
      procedure :: write_doc_all => schema_write_doc_all
      procedure :: write_doc_short => schema_write_doc_short
      procedure :: write_resolved => schema_write_resolved
      procedure :: render_markdown => schema_render_markdown
      procedure :: render_json => schema_render_json
      procedure :: find_group => schema_find_group
   end type nml_schema_t
contains

   ! ------------------------------------------------------------------
   ! Constructors.  Each captures `default = current target value`.
   ! ------------------------------------------------------------------

   function nml_real(name, tgt, doc, units, required, min, max, dead_on_ocean_path) result(key)
      !! Construct a real key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      real(rk), pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      real(rk), intent(in), optional :: min, max
      character(len=*), intent(in), optional :: dead_on_ocean_path
      type(nml_real_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
      if (present(min)) then
         key%vmin = min
         key%has_min = .true.
      end if
      if (present(max)) then
         key%vmax = max
         key%has_max = .true.
      end if
   end function nml_real

   function nml_int(name, tgt, doc, units, required, min, max, dead_on_ocean_path) result(key)
      !! Construct an integer key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      integer, pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      integer, intent(in), optional :: min, max
      character(len=*), intent(in), optional :: dead_on_ocean_path
      type(nml_int_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
      if (present(min)) then
         key%vmin = min
         key%has_min = .true.
      end if
      if (present(max)) then
         key%vmax = max
         key%has_max = .true.
      end if
   end function nml_int

   function nml_logical(name, tgt, doc, units, required, dead_on_ocean_path) result(key)
      !! Construct a logical key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      logical, pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      character(len=*), intent(in), optional :: dead_on_ocean_path
      type(nml_logical_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
   end function nml_logical

   function nml_string(name, tgt, doc, units, required, dead_on_ocean_path) result(key)
      !! Construct a string key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      character(len=:), pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      character(len=*), intent(in), optional :: dead_on_ocean_path
      type(nml_string_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
   end function nml_string

   function nml_enum(name, tgt, doc, allowed, units, required, dead_on_ocean_path, &
                     retired, retired_hint) result(key)
      !! Construct an enum key restricted to `allowed`, capturing the
      !! current target as default.
      !!
      !! `retired` lists spellings that were legal in an earlier release.
      !! They are NOT accepted — they fail loud, quoting `retired_hint`,
      !! so a stale namelist is told what the value is called now instead
      !! of being handed the generic allowed-set message.
      character(len=*), intent(in) :: name, doc
      character(len=:), pointer, intent(in) :: tgt
      character(len=*), intent(in) :: allowed(:)
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      character(len=*), intent(in), optional :: dead_on_ocean_path
      character(len=*), intent(in), optional :: retired(:)
      character(len=*), intent(in), optional :: retired_hint
      type(nml_enum_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
      ! Explicit allocate + element copy on purpose: the auto-allocating
      ! `key%allowed = allowed` makes NVHPC's str_copy write one byte past
      ! the allocation (heap corruption once enough enums register).
      block
         integer :: i
         allocate (character(len=len(allowed)) :: key%allowed(size(allowed)))
         do i = 1, size(allowed)
            key%allowed(i) = allowed(i)
         end do
      end block
      if (present(retired)) then
         block
            integer :: i
            allocate (character(len=len(retired)) :: key%retired(size(retired)))
            do i = 1, size(retired)
               key%retired(i) = retired(i)
            end do
         end block
      end if
      if (present(retired_hint)) key%retired_hint = retired_hint
   end function nml_enum

   function nml_real_array(name, tgt, doc, units, required, dead_on_ocean_path) result(key)
      !! Construct a real-array key, capturing the current target as
      !! default (whole-array compare for is_default).
      character(len=*), intent(in) :: name, doc
      real(rk), pointer, intent(in) :: tgt(:)
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      character(len=*), intent(in), optional :: dead_on_ocean_path
      type(nml_real_array_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%dead_reason = ""
      if (present(dead_on_ocean_path)) key%dead_reason = dead_on_ocean_path
      key%tgt => tgt
      key%default = tgt
   end function nml_real_array
   ! ------------------------------------------------------------------
   ! Formatting helpers.
   ! ------------------------------------------------------------------

   pure function fmt_real(v) result(str)
      !! Format a real in a stable, namelist-valid exponential form.
      real(rk), intent(in) :: v
      character(len=:), allocatable :: str
      character(len=32) :: buf
      write (buf, "(E18.10)") v
      str = trim(adjustl(buf))
   end function fmt_real

   pure function fmt_int(v) result(str)
      !! Format an integer.
      integer, intent(in) :: v
      character(len=:), allocatable :: str
      character(len=32) :: buf
      write (buf, "(I0)") v
      str = trim(adjustl(buf))
   end function fmt_int

   pure function fmt_logical(v) result(str)
      !! Format a logical as namelist `.true.`/`.false.`.
      logical, intent(in) :: v
      character(len=:), allocatable :: str
      if (v) then
         str = ".true."
      else
         str = ".false."
      end if
   end function fmt_logical

   pure function lower(s) result(out)
      !! Lowercase an ASCII string.
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            out(i:i) = achar(c + 32)
         else
            out(i:i) = s(i:i)
         end if
      end do
   end function lower

   ! ------------------------------------------------------------------
   ! Real key methods.
   ! ------------------------------------------------------------------

   subroutine real_parse(this, values, err_msg)
      class(nml_real_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      real(rk) :: v
      integer :: ios
      character(len=128) :: imsg

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      read (values(1), *, iostat=ios, iomsg=imsg) v
      if (ios /= 0) then
         err_msg = "key '"//this%name//"': cannot parse real from '"//trim(values(1))//"'"
         return
      end if
      ! A declared range is a promise that the value is a usable number:
      ! `v < vmin` / `v > vmax` are both FALSE for NaN, so a NaN (and, with
      ! only one bound declared, an infinity on the other side) must be
      ! refused by name, never by comparison.
      if ((this%has_min .or. this%has_max) .and. .not. ieee_is_finite(v)) then
         err_msg = "key '"//this%name//"' = "//fmt_real(v)//" is not finite (the key has a declared range)"
         return
      end if
      if (this%has_min .and. v < this%vmin) then
         err_msg = "key '"//this%name//"' = "//fmt_real(v)//" below min "//fmt_real(this%vmin)
         return
      end if
      if (this%has_max .and. v > this%vmax) then
         err_msg = "key '"//this%name//"' = "//fmt_real(v)//" above max "//fmt_real(this%vmax)
         return
      end if
      this%tgt = v
   end subroutine real_parse

   function real_value_string(this) result(str)
      class(nml_real_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_real(this%tgt)
   end function real_value_string

   function real_default_string(this) result(str)
      class(nml_real_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_real(this%default)
   end function real_default_string

   function real_is_default(this) result(yes)
      class(nml_real_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt == this%default)
   end function real_is_default

   ! ------------------------------------------------------------------
   ! Integer key methods.
   ! ------------------------------------------------------------------

   subroutine int_parse(this, values, err_msg)
      class(nml_int_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      integer :: v
      integer(int64) :: wide
      integer :: ios
      character(len=128) :: imsg

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      read (values(1), *, iostat=ios, iomsg=imsg) wide
      if (ios /= 0) then
         err_msg = "key '"//this%name//"': cannot parse integer from '"//trim(values(1))//"'"
         return
      end if
      ! Read WIDE and range-check: a narrow `read` of an out-of-range value
      ! is an error on gfortran but wraps silently on nvfortran.
      if (wide > int(huge(v), int64) .or. wide < -int(huge(v), int64) - 1_int64) then
         err_msg = "key '"//this%name//"': integer '"//trim(values(1))//"' out of range"
         return
      end if
      v = int(wide)
      if (this%has_min .and. v < this%vmin) then
         err_msg = "key '"//this%name//"' = "//fmt_int(v)//" below min "//fmt_int(this%vmin)
         return
      end if
      if (this%has_max .and. v > this%vmax) then
         err_msg = "key '"//this%name//"' = "//fmt_int(v)//" above max "//fmt_int(this%vmax)
         return
      end if
      this%tgt = v
   end subroutine int_parse

   function int_value_string(this) result(str)
      class(nml_int_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_int(this%tgt)
   end function int_value_string

   function int_default_string(this) result(str)
      class(nml_int_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_int(this%default)
   end function int_default_string

   function int_is_default(this) result(yes)
      class(nml_int_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt == this%default)
   end function int_is_default
   ! ------------------------------------------------------------------
   ! Logical key methods.
   ! ------------------------------------------------------------------

   subroutine logical_parse(this, values, err_msg)
      class(nml_logical_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: tok

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      tok = lower(trim(adjustl(values(1))))
      select case (tok)
      case (".true.", "t", "true")
         this%tgt = .true.
      case (".false.", "f", "false")
         this%tgt = .false.
      case default
         err_msg = "key '"//this%name//"': cannot parse logical from '"//trim(values(1))//"'"
      end select
   end subroutine logical_parse

   function logical_value_string(this) result(str)
      class(nml_logical_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_logical(this%tgt)
   end function logical_value_string

   function logical_default_string(this) result(str)
      class(nml_logical_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_logical(this%default)
   end function logical_default_string

   function logical_is_default(this) result(yes)
      class(nml_logical_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt .eqv. this%default)
   end function logical_is_default

   ! ------------------------------------------------------------------
   ! String key methods.
   ! ------------------------------------------------------------------

   subroutine string_parse(this, values, err_msg)
      class(nml_string_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: v

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      v = trim(values(1))
      if (len(v) > len(this%tgt)) then
         err_msg = "key '"//this%name//"': value '"//v//"' ("//fmt_int(len(v))// &
                   " chars) exceeds target length "//fmt_int(len(this%tgt))
         return
      end if
      this%tgt = v
   end subroutine string_parse

   function string_value_string(this) result(str)
      class(nml_string_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = '"'//trim(this%tgt)//'"'
   end function string_value_string

   function string_default_string(this) result(str)
      class(nml_string_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = '"'//trim(this%default)//'"'
   end function string_default_string

   function string_is_default(this) result(yes)
      class(nml_string_key_t), intent(in) :: this
      logical :: yes
      yes = (trim(this%tgt) == trim(this%default))
   end function string_is_default

   ! ------------------------------------------------------------------
   ! Enum key methods (extends string).
   ! ------------------------------------------------------------------

   subroutine enum_parse(this, values, err_msg)
      class(nml_enum_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: v, list
      integer :: i

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      v = trim(adjustl(values(1)))
      do i = 1, size(this%allowed)
         if (lower(v) == lower(trim(this%allowed(i)))) then
            if (len(trim(this%allowed(i))) > len(this%tgt)) then
               err_msg = "key '"//this%name//"': canonical value '"//trim(this%allowed(i))// &
                         "' exceeds target length "//fmt_int(len(this%tgt))
               return
            end if
            this%tgt = trim(this%allowed(i))
            return
         end if
      end do
      ! A spelling that WAS legal gets the migration message, never the
      ! generic allowed-set one: "not in allowed set" tells the reader the
      ! value is wrong but not what it is called now, which for a rename is
      ! the only thing they need.
      if (allocated(this%retired)) then
         do i = 1, size(this%retired)
            if (lower(v) == lower(trim(this%retired(i)))) then
               err_msg = "key '"//this%name//"': '"//v//"' was RETIRED"
               if (allocated(this%retired_hint)) then
                  err_msg = err_msg//" — "//this%retired_hint
               end if
               return
            end if
         end do
      end if
      list = ""
      do i = 1, size(this%allowed)
         if (i > 1) list = list//", "
         list = list//"'"//trim(this%allowed(i))//"'"
      end do
      err_msg = "key '"//this%name//"': '"//v//"' not in allowed set {"//list//"}"
   end subroutine enum_parse

   ! ------------------------------------------------------------------
   ! Real-array key methods.
   ! ------------------------------------------------------------------

   subroutine real_array_parse(this, values, err_msg)
      class(nml_real_array_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      real(rk) :: v
      integer :: i, ios
      character(len=128) :: imsg

      if (size(values) > size(this%tgt)) then
         err_msg = "key '"//this%name//"' accepts at most "//fmt_int(size(this%tgt))// &
                   " values, got "//fmt_int(size(values))
         return
      end if
      if (size(values) < 1) then
         err_msg = "key '"//this%name//"' expects at least 1 value"
         return
      end if
      do i = 1, size(values)
         read (values(i), *, iostat=ios, iomsg=imsg) v
         if (ios /= 0) then
            err_msg = "key '"//this%name//"': cannot parse real from '"//trim(values(i))// &
                      "' (element "//fmt_int(i)//")"
            return
         end if
         this%tgt(i) = v
      end do
   end subroutine real_array_parse

   function real_array_value_string(this) result(str)
      class(nml_real_array_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(this%tgt)
         if (i > 1) str = str//", "
         str = str//fmt_real(this%tgt(i))
      end do
   end function real_array_value_string

   function real_array_default_string(this) result(str)
      class(nml_real_array_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(this%default)
         if (i > 1) str = str//", "
         str = str//fmt_real(this%default(i))
      end do
   end function real_array_default_string

   function real_array_is_default(this) result(yes)
      class(nml_real_array_key_t), intent(in) :: this
      logical :: yes
      yes = .false.
      if (size(this%tgt) /= size(this%default)) return
      yes = all(this%tgt == this%default)
   end function real_array_is_default
   ! ------------------------------------------------------------------
   ! Group management.
   ! ------------------------------------------------------------------

   subroutine group_add(this, key)
      !! Add a key to the group (sourced-allocate into a box).  Error
      !! stop on duplicate key name (programming error).
      class(nml_group_t), intent(inout) :: this
      class(nml_key_t), intent(in) :: key
      type(nml_key_box_t), allocatable :: tmp(:)
      integer :: n, i

      if (.not. allocated(this%keys)) allocate (this%keys(0))
      n = size(this%keys)
      do i = 1, n
         if (lower(this%keys(i)%key%name) == lower(key%name)) then
            call global_logger%error("duplicate key '"//key%name//"' in group '"//this%name//"'")
            error stop "rdb_nml_schema: duplicate key registration"
         end if
      end do
      allocate (tmp(n + 1))
      do i = 1, n
         call move_alloc(this%keys(i)%key, tmp(i)%key)
      end do
      allocate (tmp(n + 1)%key, source=key)
      call move_alloc(tmp, this%keys)
   end subroutine group_add

   function group_find_key(this, name) result(idx)
      !! Index of `name` in the group's key list, 0 if absent.
      class(nml_group_t), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: idx, i
      idx = 0
      if (.not. allocated(this%keys)) return
      do i = 1, size(this%keys)
         if (lower(this%keys(i)%key%name) == lower(name)) then
            idx = i
            return
         end if
      end do
   end function group_find_key

   ! ------------------------------------------------------------------
   ! Schema management.
   ! ------------------------------------------------------------------

   subroutine schema_add_group(this, g)
      !! MOVE a locally-built group into the schema — `g` is consumed
      !! (allocatables transferred). Move semantics on purpose: intrinsic
      !! assignment deep-copies the polymorphic key boxes, which NVHPC
      !! miscompiles (heap corruption); move_alloc transfers descriptors
      !! only. Error stop on duplicate name.
      class(nml_schema_t), intent(inout) :: this
      type(nml_group_t), intent(inout) :: g
      type(nml_group_t), allocatable :: tmp(:)
      integer :: n, i

      if (.not. allocated(this%groups)) allocate (this%groups(0))
      n = size(this%groups)
      do i = 1, n
         if (lower(this%groups(i)%name) == lower(g%name)) then
            call global_logger%error("duplicate group '"//g%name//"' in schema")
            error stop "rdb_nml_schema: duplicate group registration"
         end if
      end do
      allocate (tmp(n + 1))
      do i = 1, n
         call move_group(this%groups(i), tmp(i))
      end do
      call move_group(g, tmp(n + 1))
      call move_alloc(tmp, this%groups)
   end subroutine schema_add_group

   subroutine move_group(src, dst)
      !! Transfer a group's contents without copying any polymorphic
      !! key box (descriptor moves only — see [[schema_add_group]]).
      type(nml_group_t), intent(inout) :: src
      type(nml_group_t), intent(inout) :: dst
      if (allocated(src%name)) call move_alloc(src%name, dst%name)
      if (allocated(src%doc)) call move_alloc(src%doc, dst%doc)
      if (allocated(src%keys)) call move_alloc(src%keys, dst%keys)
      dst%found = src%found
      dst%cross_check => src%cross_check
      src%cross_check => null()
   end subroutine move_group

   function schema_find_group(this, name) result(idx)
      !! Index of `name` in the schema's groups (exact, case-sensitive), 0 if absent.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: idx, i
      idx = 0
      if (.not. allocated(this%groups)) return
      do i = 1, size(this%groups)
         if (this%groups(i)%name == name) then
            idx = i
            return
         end if
      end do
   end function schema_find_group

   ! ------------------------------------------------------------------
   ! Levenshtein + did-you-mean.
   ! ------------------------------------------------------------------

   pure function levenshtein(a, b) result(dist)
      !! Levenshtein edit distance between two (capped) strings.
      character(len=*), intent(in) :: a, b
      integer :: dist
      integer :: la, lb, i, j
      integer, allocatable :: prev(:), cur(:)
      character(len=MAX_SUGGEST_LEN) :: aa, bb

      aa = a
      bb = b
      la = min(len_trim(a), MAX_SUGGEST_LEN)
      lb = min(len_trim(b), MAX_SUGGEST_LEN)
      allocate (prev(0:lb), cur(0:lb))
      do j = 0, lb
         prev(j) = j
      end do
      do i = 1, la
         cur(0) = i
         do j = 1, lb
            if (aa(i:i) == bb(j:j)) then
               cur(j) = prev(j - 1)
            else
               cur(j) = 1 + min(prev(j), cur(j - 1), prev(j - 1))
            end if
         end do
         prev = cur
      end do
      dist = prev(lb)
   end function levenshtein

   function suggest(name, candidates) result(msg)
      !! Build a " (did you mean 'x'?)" fragment if a close candidate
      !! exists (distance <= max(2, len/3)); else empty.
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: candidates(:)
      character(len=:), allocatable :: msg
      integer :: i, d, best_d, best_i, thresh

      msg = ""
      best_d = huge(0)
      best_i = 0
      thresh = max(2, len_trim(name)/3)
      do i = 1, size(candidates)
         d = levenshtein(lower(trim(name)), lower(trim(candidates(i))))
         if (d < best_d) then
            best_d = d
            best_i = i
         end if
      end do
      if (best_i > 0 .and. best_d <= thresh) then
         msg = " (did you mean '"//trim(candidates(best_i))//"'?)"
      end if
   end function suggest
   ! ------------------------------------------------------------------
   ! Documentation writers.
   ! ------------------------------------------------------------------

   subroutine write_key_line(unit, box, short_only, wrote_any)
      !! Write one aligned `key = value  ! [units] doc` line.  When
      !! `short_only`, skip default-valued keys (and leave wrote_any).
      integer, intent(in) :: unit
      type(nml_key_box_t), intent(in) :: box
      logical, intent(in) :: short_only
      logical, intent(inout) :: wrote_any
      character(len=:), allocatable :: lhs, comment, ustr
      logical :: is_def

      is_def = box%key%is_default()
      if (short_only .and. is_def) return
      wrote_any = .true.

      lhs = "  "//box%key%name//" = "//box%key%value_string()
      ustr = ""
      if (len_trim(box%key%units) > 0) ustr = "["//trim(box%key%units)//"] "
      if (is_def) then
         comment = "  ! "//ustr//trim(box%key%doc)//" (default)"
      else
         comment = "  ! "//ustr//trim(box%key%doc)// &
                   " *** non-default (default: "//box%key%default_string()//")"
      end if
      write (unit, "(A)") lhs//comment
   end subroutine write_key_line

   subroutine schema_write_doc_all(this, path)
      !! Write every group and key as valid namelist with aligned docs.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik
      logical :: dummy

      open (newunit=u, file=path, status="replace", action="write")
      do ig = 1, size(this%groups)
         write (u, "(A)") "&"//this%groups(ig)%name//"_nml"
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               dummy = .false.
               call write_key_line(u, this%groups(ig)%keys(ik), .false., dummy)
            end do
         end if
         write (u, "(A)") "/"
      end do
      close (u)
   end subroutine schema_write_doc_all

   subroutine schema_write_doc_short(this, path)
      !! Write only non-default keys (MOM6 parameter_doc.short).  Groups
      !! with no non-default keys are omitted entirely.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik
      logical :: wrote_any

      open (newunit=u, file=path, status="replace", action="write")
      do ig = 1, size(this%groups)
         wrote_any = .false.
         if (allocated(this%groups(ig)%keys)) then
            ! First pass: does any non-default key exist in this group?
            do ik = 1, size(this%groups(ig)%keys)
               if (.not. this%groups(ig)%keys(ik)%key%is_default()) wrote_any = .true.
            end do
         end if
         if (.not. wrote_any) cycle
         write (u, "(A)") "&"//this%groups(ig)%name//"_nml"
         do ik = 1, size(this%groups(ig)%keys)
            call write_key_line(u, this%groups(ig)%keys(ik), .true., wrote_any)
         end do
         write (u, "(A)") "/"
      end do
      close (u)
   end subroutine schema_write_doc_short

   subroutine schema_write_resolved(this, path)
      !! Write every key as `group.key = value` with reals at full precision
      !! (`ES30.17E3`: 18 significant digits, an exact binary64 round trip),
      !! where the doc dumps round them to 11.  A function of the resolved
      !! configuration only -- it is the equality oracle between input
      !! formats (`rdb_cfg_dump`).
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik

      open (newunit=u, file=path, status="replace", action="write")
      do ig = 1, size(this%groups)
         if (.not. allocated(this%groups(ig)%keys)) cycle
         do ik = 1, size(this%groups(ig)%keys)
            write (u, "(A)") this%groups(ig)%name//"."// &
               this%groups(ig)%keys(ik)%key%name//" = "// &
               exact_value_string(this%groups(ig)%keys(ik)%key)
         end do
      end do
      close (u)
   end subroutine schema_write_resolved

   function exact_value_string(key) result(str)
      !! `value_string` with reals at full precision.
      class(nml_key_t), intent(in) :: key
      character(len=:), allocatable :: str
      character(len=32) :: buf
      integer :: i

      select type (key)
      type is (nml_real_key_t)
         write (buf, "(ES30.17E3)") key%tgt
         str = trim(adjustl(buf))
      type is (nml_real_array_key_t)
         str = ""
         do i = 1, size(key%tgt)
            write (buf, "(ES30.17E3)") key%tgt(i)
            if (i > 1) str = str//", "
            str = str//trim(adjustl(buf))
         end do
      class default
         str = key%value_string()
      end select
   end function exact_value_string

   subroutine schema_render_markdown(this, path)
      !! Render the schema as Markdown: per group a heading, doc line,
      !! and a knob table.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik
      character(len=:), allocatable :: ustr

      open (newunit=u, file=path, status="replace", action="write")
      do ig = 1, size(this%groups)
         ! Blank separator BETWEEN groups only: a trailing blank line after
         ! the last group gets re-trimmed by the end-of-file pre-commit hook
         ! on every regeneration (a recurring commit-abort nuisance).
         if (ig > 1) write (u, "(A)") ""
         write (u, "(A)") "### &"//this%groups(ig)%name//"_nml"
         write (u, "(A)") ""
         if (allocated(this%groups(ig)%doc)) then
            if (len_trim(this%groups(ig)%doc) > 0) then
               write (u, "(A)") trim(this%groups(ig)%doc)
               write (u, "(A)") ""
            end if
         end if
         write (u, "(A)") "| Knob | Default | Units | Description |"
         write (u, "(A)") "|------|---------|-------|-------------|"
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               associate (k => this%groups(ig)%keys(ik)%key)
                  ustr = ""
                  if (len_trim(k%units) > 0) ustr = trim(k%units)
                  write (u, "(A)") "| `"//k%name//"` | `"//k%default_string()// &
                     "` | "//ustr//" | "//trim(k%doc)//" |"
               end associate
            end do
         end if
      end do
      close (u)
   end subroutine schema_render_markdown

   ! ------------------------------------------------------------------
   ! JSON emitter (P4: python/rdb/_config_generated.py generator input).
   ! ------------------------------------------------------------------

   pure function json_escape(s) result(out)
      !! Minimal JSON string escaping: backslash, double-quote, and the
      !! common control characters.  Doc/name/units strings in this
      !! codebase are plain single-line ASCII, so this is deliberately
      !! not a general-purpose escaper.
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: out
      integer :: i
      character :: c
      out = ""
      do i = 1, len(s)
         c = s(i:i)
         select case (c)
         case ('"')
            out = out//'\"'
         case ("\")
            out = out//"\\"
         case (achar(10))
            out = out//"\n"
         case (achar(9))
            out = out//"\t"
         case default
            if (iachar(c) >= 32) out = out//c
         end select
      end do
   end function json_escape

   pure function fmt_json_bool(v) result(str)
      !! Format a logical as a JSON `true`/`false` literal (as opposed to
      !! `fmt_logical`'s Fortran `.true.`/`.false.`).
      logical, intent(in) :: v
      character(len=:), allocatable :: str
      if (v) then
         str = "true"
      else
         str = "false"
      end if
   end function fmt_json_bool

   pure function json_fmt_real(v) result(str)
      !! Real -> JSON number, robust to extreme magnitude (e.g. `huge(wp)`
      !! floor/ceiling sentinels like `kv_max`'s default): `fmt_real`'s
      !! `E18.10` is too narrow for a 3-digit exponent and gfortran drops
      !! the `E` literal to fit ("0.1797693135+309"), which is not valid
      !! JSON. `ES.E3` explicitly requests a 3-digit exponent field, which
      !! forces the `E` to stay, at a width sized for the worst case. 16
      !! fraction digits (17 significant figures) is real64's round-trip
      !! guarantee -- one digit short (as an earlier version of this
      !! function used) rounds `huge(wp)`'s text UP past the true maximum
      !! finite double, so re-parsing it (e.g. Python's `json.load`)
      !! overflows to infinity, which is not valid JSON and breaks the
      !! generator (`tools/gen_python_config.py`) that reads this file.
      real(rk), intent(in) :: v
      character(len=:), allocatable :: str
      character(len=40) :: buf
      write (buf, "(ES26.16E3)") v
      str = trim(adjustl(buf))
   end function json_fmt_real

   pure function json_string_array(arr) result(str)
      !! Render `arr(:)` as a comma-separated list of JSON-quoted strings
      !! (no enclosing brackets -- the caller supplies those).
      character(len=*), intent(in) :: arr(:)
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(arr)
         if (i > 1) str = str//", "
         str = str//'"'//json_escape(trim(arr(i)))//'"'
      end do
   end function json_string_array

   pure function json_real_array(arr) result(str)
      !! Render `arr(:)` as a comma-separated list of JSON numbers (no
      !! enclosing brackets).
      real(rk), intent(in) :: arr(:)
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(arr)
         if (i > 1) str = str//", "
         str = str//json_fmt_real(arr(i))
      end do
   end function json_real_array

   subroutine write_key_json(u, box, is_last)
      !! Write one key as a JSON object.  A `select type` over the
      !! concrete key recovers the kind-specific fields (vmin/vmax,
      !! allowed(:), array size) that `nml_key_t`'s abstract interface
      !! does not carry -- the one place this walk needs the concrete
      !! type (`schema_render_markdown` does not, since `default_string()`
      !! is polymorphic).
      integer, intent(in) :: u
      type(nml_key_box_t), intent(in) :: box
      logical, intent(in) :: is_last
      character(len=:), allocatable :: dead_str

      dead_str = ""
      if (allocated(box%key%dead_reason)) dead_str = trim(box%key%dead_reason)

      write (u, "(A)") "        {"
      write (u, "(A)") '          "name": "'//json_escape(box%key%name)//'",'
      write (u, "(A)") '          "doc": "'//json_escape(trim(box%key%doc))//'",'
      write (u, "(A)") '          "units": "'//json_escape(trim(box%key%units))//'",'
      write (u, "(A)") '          "required": '//fmt_json_bool(box%key%required)//","
      write (u, "(A)") '          "dead_on_ocean_path": "'//json_escape(dead_str)//'",'

      select type (k => box%key)
      type is (nml_enum_key_t)
         write (u, "(A)") '          "kind": "enum",'
         write (u, "(A)") '          "default": "'//json_escape(trim(k%default))//'",'
         write (u, "(A)") '          "allowed": ['//json_string_array(k%allowed)//"]"
      type is (nml_string_key_t)
         write (u, "(A)") '          "kind": "string",'
         write (u, "(A)") '          "max_len": '//fmt_int(len(k%tgt))//","
         write (u, "(A)") '          "default": "'//json_escape(trim(k%default))//'"'
      type is (nml_real_array_key_t)
         write (u, "(A)") '          "kind": "real_array",'
         write (u, "(A)") '          "size": '//fmt_int(size(k%tgt))//","
         write (u, "(A)") '          "default": ['//json_real_array(k%default)//"]"
      type is (nml_real_key_t)
         write (u, "(A)") '          "kind": "real",'
         write (u, "(A)") '          "default": '//json_fmt_real(k%default)//","
         write (u, "(A)") '          "has_min": '//fmt_json_bool(k%has_min)//","
         write (u, "(A)") '          "vmin": '//json_fmt_real(k%vmin)//","
         write (u, "(A)") '          "has_max": '//fmt_json_bool(k%has_max)//","
         write (u, "(A)") '          "vmax": '//json_fmt_real(k%vmax)
      type is (nml_int_key_t)
         write (u, "(A)") '          "kind": "int",'
         write (u, "(A)") '          "default": '//fmt_int(k%default)//","
         write (u, "(A)") '          "has_min": '//fmt_json_bool(k%has_min)//","
         write (u, "(A)") '          "vmin": '//fmt_int(k%vmin)//","
         write (u, "(A)") '          "has_max": '//fmt_json_bool(k%has_max)//","
         write (u, "(A)") '          "vmax": '//fmt_int(k%vmax)
      type is (nml_logical_key_t)
         write (u, "(A)") '          "kind": "logical",'
         write (u, "(A)") '          "default": '//fmt_json_bool(k%default)
      class default
         write (u, "(A)") '          "kind": "unknown"'
      end select

      if (is_last) then
         write (u, "(A)") "        }"
      else
         write (u, "(A)") "        },"
      end if
   end subroutine write_key_json

   subroutine schema_render_json(this, path)
      !! Emit the schema as JSON: one object per group, nesting one object
      !! per key.  This is the ONLY thing `tools/gen_python_config.py`
      !! reads -- it parses no Fortran source, mirroring the precedent
      !! `schema_render_markdown` already set for `docs/generated_nml_knobs.md`.
      !! `&ocean_bc_nml` is NOT here -- it is registered via
      !! `add_external_group` (see `rdb_config.F90`), so it carries no keys
      !! on this schema; python/rdb/_config_bc.py is a hand-written
      !! stub with its own drift test (see docs D5.7).
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik, nk

      open (newunit=u, file=path, status="replace", action="write")
      write (u, "(A)") "{"
      write (u, "(A)") '  "groups": ['
      do ig = 1, size(this%groups)
         write (u, "(A)") "    {"
         write (u, "(A)") '      "name": "'//json_escape(this%groups(ig)%name)//'",'
         write (u, "(A)") '      "doc": "'//json_escape(trim(this%groups(ig)%doc))//'",'
         write (u, "(A)") '      "keys": ['
         nk = 0
         if (allocated(this%groups(ig)%keys)) nk = size(this%groups(ig)%keys)
         do ik = 1, nk
            call write_key_json(u, this%groups(ig)%keys(ik), ik == nk)
         end do
         write (u, "(A)") "      ]"
         if (ig < size(this%groups)) then
            write (u, "(A)") "    },"
         else
            write (u, "(A)") "    }"
         end if
      end do
      write (u, "(A)") "  ]"
      write (u, "(A)") "}"
      close (u)
   end subroutine schema_render_json

   ! ------------------------------------------------------------------
   ! Error collection helper.
   ! ------------------------------------------------------------------

   subroutine add_err(errors, n_err, msg)
      !! Append a message to a growable error list.
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err
      character(len=*), intent(in) :: msg
      character(len=:), allocatable :: tmp(:)
      integer :: cap, newlen, i

      if (.not. allocated(errors)) allocate (character(len=256) :: errors(8))
      cap = size(errors)
      if (n_err >= cap .or. len(msg) > len(errors)) then
         newlen = max(len(errors), len(msg))
         allocate (character(len=newlen) :: tmp(max(2*cap, n_err + 1)))
         do i = 1, n_err
            tmp(i) = errors(i)
         end do
         call move_alloc(tmp, errors)
      end if
      n_err = n_err + 1
      errors(n_err) = msg
   end subroutine add_err

   ! ------------------------------------------------------------------
   ! Post-parse passes shared by the front ends.
   ! ------------------------------------------------------------------

   subroutine post_parse_validate(this, label, err_buf, n_err, status, errors, skip_checks)
      !! Required-key + cross-check pass shared by the file and buffer
      !! parsers, then the error-return / error-stop finalization.
      !! `skip_checks` (default false) omits the two passes: a source that
      !! never parsed (a TOML syntax error) has applied nothing, so they
      !! would only add noise to the one real error.
      type(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: err_buf(:)
      integer, intent(inout) :: n_err
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)
      logical, intent(in), optional :: skip_checks

      integer :: i, ig
      logical :: testable, do_checks

      testable = present(status) .and. present(errors)
      do_checks = .true.
      if (present(skip_checks)) do_checks = .not. skip_checks

      ! Required-key check (only for groups present in the source).
      do ig = 1, size(this%groups)
         if (.not. do_checks) exit
         if (.not. this%groups(ig)%found) cycle
         if (.not. allocated(this%groups(ig)%keys)) cycle
         do i = 1, size(this%groups(ig)%keys)
            associate (k => this%groups(ig)%keys(i)%key)
               if (k%required .and. .not. k%found) then
                  call add_err(err_buf, n_err, label//": required key '"//k%name// &
                               "' missing in group '"//this%groups(ig)%name//"'")
               end if
            end associate
         end do
      end do

      ! Cross-checks: run for every registered group regardless of presence.
      do ig = 1, size(this%groups)
         if (.not. do_checks) exit
         if (associated(this%groups(ig)%cross_check)) then
            block
               character(len=:), allocatable :: cmsg
               cmsg = this%groups(ig)%cross_check()
               if (len_trim(cmsg) > 0) then
                  call add_err(err_buf, n_err, label//": group '"// &
                               this%groups(ig)%name//"': "//trim(cmsg))
               end if
            end block
         end if
      end do

      if (testable) then
         status = n_err
         if (n_err > 0) then
            allocate (character(len=len(err_buf)) :: errors(n_err))
            do i = 1, n_err
               errors(i) = err_buf(i)
            end do
         end if
         return
      end if

      if (n_err > 0) then
         do i = 1, n_err
            call error_ring_push(trim(err_buf(i)))
            call global_logger%error(trim(err_buf(i)))
         end do
         error stop "rdb_nml_schema: namelist validation failed"
      end if
   end subroutine post_parse_validate

   subroutine reset_state(this)
      !! Reset parse-state (found flags) before a parse.
      type(nml_schema_t), intent(inout) :: this
      integer :: ig, ik
      do ig = 1, size(this%groups)
         this%groups(ig)%found = .false.
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               this%groups(ig)%keys(ik)%key%found = .false.
            end do
         end if
      end do
   end subroutine reset_state

   pure function make_prefix(path, line) result(p)
      !! Build a `path:line: ` error prefix.
      character(len=*), intent(in) :: path
      integer, intent(in) :: line
      character(len=:), allocatable :: p
      p = path//":"//fmt_int(line)//": "
   end function make_prefix

   function group_key_list(g) result(names)
      !! All key names of a group as a char array.
      type(nml_group_t), intent(in) :: g
      character(len=:), allocatable :: names(:)
      integer :: nk, i, mlen

      nk = 0
      if (allocated(g%keys)) nk = size(g%keys)
      mlen = 1
      do i = 1, nk
         mlen = max(mlen, len(g%keys(i)%key%name))
      end do
      allocate (character(len=mlen) :: names(nk))
      do i = 1, nk
         names(i) = g%keys(i)%key%name
      end do
   end function group_key_list

end module rdb_nml_schema
