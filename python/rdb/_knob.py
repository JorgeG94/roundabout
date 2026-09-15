"""P4 -- the descriptor machinery behind the generated config layer.

Hand-written and stable; ``python/rdb/_config_generated.py`` (generated
by ``tools/gen_python_config.py`` from the live ``nml_schema_t``, which
now covers every namelist group including ``&ocean_bc_nml`` since its
P4.5 schema migration) imports from here.

Two rules make this non-drifting, straight from
``tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`` D2.0
and D5.3:

1. **Absence is the default, not the Fortran value.** Every knob starts
   unset; only an explicit assignment marks it "set". ``Group._namelist_lines()``
   emits ONLY set knobs -- a freshly constructed :class:`Config` therefore
   serialises to nothing, which is what gives bit-identity to existing
   namelists (``read_config_from_string`` on an empty buffer keeps the
   pristine ``config_t`` defaults).
2. **Validation on assignment mirrors ``parse_tokens`` arm-for-arm**
   (``rdb_nml_schema.F90``): range checks for real/int, case-insensitive
   membership for enum (stored canonical), a length cap for string/array.
   It does NOT re-implement ``validate_config``'s cross-knob rules or the
   per-group ``cross_check`` callbacks -- those stay in Fortran and surface
   as :class:`~rdb._errors.ConfigValidationError` through the error ring
   when the composed namelist reaches ``read_config_from_string``.

Stdlib only -- no numpy (P3/P4 rule; see ``_ffi.py``'s module docstring).
"""

from __future__ import annotations

import warnings

from ._errors import ConfigParseError


class RdbDeadKnobWarning(UserWarning):
    """A knob was assigned that is registered on the schema (accepted,
    type/range-validated, round-trips through ``to_namelist()``) but does
    NOT reach any ocean-path behaviour -- see each knob's
    ``dead_on_ocean_path`` docstring note for the specific reason. Found by
    the P4 generator sweep (2026-09-10); see
    ``tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`` D5.6.
    """


class _KnobDescriptor:
    """Base data descriptor for one namelist key.

    Per-instance state lives in the OWNING INSTANCE's ``__dict__`` under a
    private key (``self._attr``), never under the public knob name --
    that name is permanently occupied by this descriptor on the class, so
    ``Group.ocean_epbl.mstar`` and (via the class) ``Config.ocean_epbl.mstar``
    both resolve here.
    """

    kind = "unknown"

    def __init__(self, name, *, doc="", units="", default=None,
                 required=False, dead_on_ocean_path=None):
        self.name = name
        self.doc = doc
        self.units = units
        self.default = default
        self.required = required
        self.dead_on_ocean_path = dead_on_ocean_path
        self._attr = "_v_" + name

    def __set_name__(self, owner, attr_name):
        self._pyattr = attr_name

    def __get__(self, obj, objtype=None):
        if obj is None:
            # Accessed on the CLASS (or via Config's class-level group
            # attributes, see _config.py) -- return the descriptor itself
            # so `Config.ocean_epbl.mstar.default` etc. work.
            return self
        return obj.__dict__.get(self._attr, MISSING)

    def __set__(self, obj, value):
        value = self._validate(value)
        if self.dead_on_ocean_path:
            warnings.warn(
                f"{self.name!r} ({self.doc}) has no effect on the ocean "
                f"path: {self.dead_on_ocean_path}",
                RdbDeadKnobWarning, stacklevel=2)
        obj.__dict__[self._attr] = value
        explicit = obj.__dict__.setdefault("_explicit", [])
        if self.name not in explicit:
            explicit.append(self.name)

    def _validate(self, value):  # pragma: no cover - overridden
        raise NotImplementedError

    def _format(self, value) -> str:  # pragma: no cover - overridden
        raise NotImplementedError


class _Missing:
    """Sentinel type for an unset knob. A single instance (`MISSING`) is
    used everywhere; identity comparison (`is MISSING`), never equality."""

    def __repr__(self):
        return "MISSING"

    def __bool__(self):
        return False


MISSING = _Missing()


def _require_number(name, value, kind):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigParseError(
            f"key '{name}': cannot parse {kind} from {value!r} "
            f"(got {type(value).__name__})")


class Real(_KnobDescriptor):
    """Mirrors `nml_real_key_t` / `real_parse`."""

    kind = "real"

    def __init__(self, name, *, doc="", units="", default=None,
                 required=False, has_min=False, vmin=None, has_max=False,
                 vmax=None, dead_on_ocean_path=None):
        super().__init__(name, doc=doc, units=units, default=default,
                          required=required,
                          dead_on_ocean_path=dead_on_ocean_path)
        self.has_min, self.vmin = has_min, vmin
        self.has_max, self.vmax = has_max, vmax

    def _validate(self, value):
        _require_number(self.name, value, "real")
        v = float(value)
        if self.has_min and v < self.vmin:
            raise ConfigParseError(
                f"key '{self.name}' = {v!r} below min {self.vmin!r}")
        if self.has_max and v > self.vmax:
            raise ConfigParseError(
                f"key '{self.name}' = {v!r} above max {self.vmax!r}")
        return v

    def _format(self, value) -> str:
        return repr(float(value))


class Int(_KnobDescriptor):
    """Mirrors `nml_int_key_t` / `int_parse`."""

    kind = "int"

    def __init__(self, name, *, doc="", units="", default=None,
                 required=False, has_min=False, vmin=None, has_max=False,
                 vmax=None, dead_on_ocean_path=None):
        super().__init__(name, doc=doc, units=units, default=default,
                          required=required,
                          dead_on_ocean_path=dead_on_ocean_path)
        self.has_min, self.vmin = has_min, vmin
        self.has_max, self.vmax = has_max, vmax

    def _validate(self, value):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ConfigParseError(
                f"key '{self.name}': cannot parse int from {value!r} "
                f"(got {type(value).__name__})")
        v = int(value)
        if self.has_min and v < self.vmin:
            raise ConfigParseError(
                f"key '{self.name}' = {v} below min {self.vmin}")
        if self.has_max and v > self.vmax:
            raise ConfigParseError(
                f"key '{self.name}' = {v} above max {self.vmax}")
        return v

    def _format(self, value) -> str:
        return str(int(value))


class Bool(_KnobDescriptor):
    """Mirrors `nml_logical_key_t` / `logical_parse`."""

    kind = "logical"

    def _validate(self, value):
        if not isinstance(value, bool):
            raise ConfigParseError(
                f"key '{self.name}': cannot parse logical from {value!r} "
                f"(got {type(value).__name__}, expected bool)")
        return value

    def _format(self, value) -> str:
        return ".true." if value else ".false."


class Str(_KnobDescriptor):
    """Mirrors `nml_string_key_t` / `string_parse` -- errors (no silent
    truncation) if the value is longer than the Fortran target's fixed
    length, when known."""

    kind = "string"

    def __init__(self, name, *, doc="", units="", default=None,
                 required=False, max_len=None, dead_on_ocean_path=None):
        super().__init__(name, doc=doc, units=units, default=default,
                          required=required,
                          dead_on_ocean_path=dead_on_ocean_path)
        self.max_len = max_len

    def _validate(self, value):
        if not isinstance(value, str):
            raise ConfigParseError(
                f"key '{self.name}': cannot parse string from {value!r} "
                f"(got {type(value).__name__})")
        if self.max_len is not None and len(value) > self.max_len:
            raise ConfigParseError(
                f"key '{self.name}': value {value!r} ({len(value)} chars) "
                f"exceeds target length {self.max_len}")
        if '"' in value:
            raise ConfigParseError(
                f"key '{self.name}': value {value!r} contains a double "
                f"quote, which the namelist tokenizer cannot escape")
        return value

    def _format(self, value) -> str:
        return f'"{value}"'


class Enum(_KnobDescriptor):
    """Mirrors `nml_enum_key_t` / `enum_parse`: case-insensitive member
    check, stored in the canonical (schema) spelling."""

    kind = "enum"

    def __init__(self, name, *, doc="", units="", default=None,
                 required=False, allowed=(), dead_on_ocean_path=None):
        super().__init__(name, doc=doc, units=units, default=default,
                          required=required,
                          dead_on_ocean_path=dead_on_ocean_path)
        self.allowed = tuple(allowed)

    def _validate(self, value):
        if not isinstance(value, str):
            raise ConfigParseError(
                f"key '{self.name}': cannot parse enum from {value!r} "
                f"(got {type(value).__name__})")
        lv = value.strip().lower()
        for canonical in self.allowed:
            if canonical.lower() == lv:
                return canonical
        allowed_list = ", ".join(repr(a) for a in self.allowed)
        raise ConfigParseError(
            f"key '{self.name}': {value!r} not in allowed set "
            f"{{{allowed_list}}}")

    def _format(self, value) -> str:
        return f'"{value}"'


class RealArray(_KnobDescriptor):
    """Mirrors `nml_real_array_key_t` / `real_array_parse`: 1..size values,
    filled from element 1 (elements beyond what was set stay at whatever
    the Fortran target already held -- default, on a fresh config_t)."""

    kind = "real_array"

    def __init__(self, name, *, doc="", units="", default=(),
                 required=False, size=None, dead_on_ocean_path=None):
        super().__init__(name, doc=doc, units=units, default=tuple(default),
                          required=required,
                          dead_on_ocean_path=dead_on_ocean_path)
        self.size = size

    def _validate(self, value):
        try:
            values = list(value)
        except TypeError:
            raise ConfigParseError(
                f"key '{self.name}': expected a sequence of floats, got "
                f"{value!r}") from None
        if len(values) < 1:
            raise ConfigParseError(
                f"key '{self.name}' expects at least 1 value")
        if self.size is not None and len(values) > self.size:
            raise ConfigParseError(
                f"key '{self.name}' accepts at most {self.size} values, "
                f"got {len(values)}")
        out = []
        for i, v in enumerate(values):
            _require_number(self.name, v, "real")
            out.append(float(v))
        return tuple(out)

    def _format(self, value) -> str:
        return ", ".join(repr(float(v)) for v in value)


class Group:
    """Base class for one namelist group (generated or hand-written).

    Subclasses set ``_nml_name`` (the bare group name, no leading ``&``,
    no ``_nml`` suffix -- the same spelling `nml_schema_t` uses) and
    declare one descriptor attribute per knob.
    """

    _nml_name: str = ""

    def _explicit_names(self):
        """Knob names assigned on THIS instance, in assignment order."""
        return list(self.__dict__.get("_explicit", ()))

    def _namelist_lines(self):
        """One `  key = value` line per explicitly-set knob."""
        lines = []
        for name in self._explicit_names():
            desc = getattr(type(self), name)
            value = self.__dict__[desc._attr]
            lines.append(f"  {name} = {desc._format(value)}")
        return lines

    def __repr__(self):
        names = self._explicit_names()
        if not names:
            return f"{type(self).__name__}()"
        body = ", ".join(f"{n}={getattr(self, n)!r}" for n in names)
        return f"{type(self).__name__}({body})"
