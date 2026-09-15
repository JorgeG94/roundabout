"""P5/P7 D2.14 -- output & diagnostics. `&ocean_diag_nml` + `&output_nml`.

Follows the P5 pattern exactly (D2.0 rule 1): a writer over the generated
config, setting namelist knobs and nothing else. Exists because P5 deferred
this object, which left NO way to turn diagnostics off from the curated
layer except the raw escape hatch (``model.config.ocean_diag.enabled =
False``) -- and a user who does not know that escape hatch exists gets the
Fortran defaults: file output ON, writing a relative ``./output/...nc``
path that does not exist on a fresh checkout (P7 F1, now fixed so it
raises instead of aborting -- but "raises on every first run" is still a
worse default than "off unless asked for").
"""

from __future__ import annotations

_TIME_OPS = {"instant", "mean", "max", "min", "integral"}
_COORDS = {"layer", "z", "zstar", "z*", "sigma", "density", "rho"}
_CADENCE_UNITS = frozenset("smhd")

_VGRIDS = ("layer", "z_fixed", "sigma", "zstar", "density")
_PRECISIONS = ("double", "single")
_LOG_LEVELS = ("debug", "verbose", "info", "performance", "warning", "error")
#: Mirrors `apply_time_unit_cascade`'s own case list (rdb_config.F90:2886)
#: -- the SAME conversion `Duration.time_unit` feeds.
_TIME_UNITS = ("s", "sec", "second", "seconds", "min", "minute", "minutes",
               "hr", "hour", "hours", "h", "day", "days", "d",
               "year", "years", "yr")

#: vgrid -> (array knob, count knob) on `&ocean_diag_nml`, for `levels=`.
_LEVEL_FIELDS = {
    "z_fixed": ("z_levels", "n_z_levels"),
    "sigma": ("sigma_levels", "n_sigma_levels"),
    "zstar": ("zstar_levels", "n_zstar_levels"),
    "density": ("rho_levels", "n_rho_levels"),
}


def _is_cadence(attr: str) -> bool:
    """`<int|float><unit>`, unit in s/m/h/d -- mirrors
    `parse_cadence_attr` (`rdb_ocean_diag.F90`)."""
    if len(attr) < 2 or attr[-1] not in _CADENCE_UNITS:
        return False
    try:
        val = float(attr[:-1])
    except ValueError:
        return False
    return val > 0


def _validate_token(token: str) -> str:
    """Mirror `parse_one_spec_token`'s grammar (`rdb_ocean_diag.F90`) so a
    malformed `&ocean_diag_nml diags` token is a Python ``ValueError`` at
    ``Diagnostics(...)`` construction, not a Fortran ``error stop`` inside
    ``create()`` -- `parse_one_spec_token`'s unknown-attribute branch is
    one of the `configure_ocean_diag` paths P0.1 explicitly did NOT
    convert to a returnable status (see its "Residual scope" note).  This
    duplicates the GRAMMAR only (never the dispatch/registration logic,
    which stays exclusively in Fortran, per D5.5) -- returns the token
    unchanged, it does not rewrite it.
    """
    parts = token.split(":")
    if not parts[0]:
        raise ValueError(
            f"diagnostic spec token {token!r} has no name before the "
            f"first ':'")
    for attr in parts[1:]:
        if attr == "":
            continue
        low = attr.lower()
        if low == "off" or low in _TIME_OPS or low in _COORDS:
            continue
        if _is_cadence(low):
            continue
        raise ValueError(
            f"unknown diagnostic attribute {attr!r} in spec token "
            f"{token!r}; valid: off | instant|mean|max|min|integral | "
            f"layer|z|zstar|sigma|density|rho | <int>s/m/h/d")
    return token


def _tokens(names) -> str:
    toks = names.split() if isinstance(names, str) else [str(n) for n in names]
    return " ".join(_validate_token(t) for t in toks)


class Diagnostics:
    """``Diagnostics(names=, every=, vgrid=, levels=, filename=,
    precision=, mask_vanished=, reproducing_sums=, deflate=,
    output_dir=, status_interval=, log_level=, enabled=True)`` ->
    ``&ocean_diag_nml`` (+ ``&output_nml`` / ``&logging_nml`` for
    ``deflate``/``output_dir``/``status_interval``/``log_level``).

    ``output_dir=``, ``status_interval=`` and ``log_level=`` are
    applied UNCONDITIONALLY (before the ``enabled=False`` early
    return below) -- they are console/output plumbing, not NetCDF
    diagnostic selection, so turning diagnostics off should not also
    silently drop them (e.g. a restart file still wants
    ``output_dir=``; the console status line still wants
    ``status_interval=`` regardless of whether the diag registry is
    writing anything). Every OTHER keyword below is diagnostic-only
    and is skipped when ``enabled=False``.

    - ``enabled=False`` is the explicit OFF switch this object exists to
      provide: sets ``&ocean_diag_nml enabled = .false.`` (the actual
      gate `engine_configure_diag` checks -- ``&output_nml
      output_to_file`` is ALSO set for documentation parity with the
      escape hatch, but is a no-op on the ocean path per the P4
      dead-knob sweep: it is validated but read nowhere in ``src/``).
      When ``enabled=False``, every other keyword is ignored (nothing
      to configure) rather than silently building unreachable knobs.
    - ``names=`` -- the unified ``&ocean_diag_nml diags`` token-list
      selection: a single pre-built token string, OR an iterable of
      tokens (``"temperature"``, ``"KE:off"``,
      ``"vorticity_z:z:1d"``, ``"temperature:6h:mean"``, ...). Each
      token is validated against the SAME grammar
      ``parse_one_spec_token`` (`rdb_ocean_diag.F90`) enforces, so a
      typo raises ``ValueError`` here rather than aborting `create()`.
    - ``every=`` -- global output cadence in seconds, ``dt_out``.
    - ``vgrid=`` -- default output vertical grid for layered
      diagnostics: one of ``"layer"``/``"z_fixed"``/``"sigma"``/
      ``"zstar"``/``"density"``.
    - ``levels=`` -- output levels for whichever ``vgrid=`` needs them
      (z-levels/sigma fractions/zstar reference depths/density bin
      edges); routed to the matching ``&ocean_diag_nml *_levels`` array
      (+ its companion count knob) by ``vgrid``. Requires ``vgrid=`` to
      be one of the four leveled families (not ``"layer"``, which has
      no vertical remap).
    - ``filename=`` -- output basename (per-rank suffix appended).
    - ``precision=`` -- ``"double"`` (default) or ``"single"``
      (halves the bytes written).
    - ``mask_vanished=`` -- mask below-bottom/pinched remap target
      cells to the missing-value sentinel instead of 0 (non-LAYER
      vgrids only -- see ``rdb_ocean_diag_fills``'s
      ``fill_tracer_impl`` docstring for the SEPARATE, unconditional
      NaN-at-LAYER-vgrid fix this does not affect).
    - ``reproducing_sums=`` -- order-invariant EFP console totals.
    - ``deflate=`` -- an int 0-9 (0 = off); ``&output_nml
      compress_output``/``compress_level``.
    - ``output_dir=`` -- ``&output_nml output_dir``, the directory
      every file this run writes (diag NetCDF, restart) lands in.
      Applied even when ``enabled=False``.
    - ``status_interval=`` -- ``&logging_nml status_interval``, the
      console ``[stats]`` print cadence, in ``time_unit`` (see
      :class:`Duration`) -- NOT the same cadence as ``every=``
      (``&ocean_diag_nml dt_out``, always seconds). Applied even when
      ``enabled=False``.
    - ``log_level=`` -- ``&logging_nml log_level``, one of ``"debug"``/
      ``"verbose"``/``"info"``/``"performance"``/``"warning"``/
      ``"error"``. Applied even when ``enabled=False``.
    """

    def __init__(self, names=None, *, every=None, vgrid=None, levels=None,
                 filename=None, precision=None, mask_vanished=None,
                 reproducing_sums=None, deflate=None, output_dir=None,
                 status_interval=None, log_level=None, enabled=True):
        if vgrid is not None and vgrid not in _VGRIDS:
            raise ValueError(f"vgrid must be one of {_VGRIDS}, got {vgrid!r}")
        if precision is not None and precision not in _PRECISIONS:
            raise ValueError(
                f"precision must be one of {_PRECISIONS}, got {precision!r}")
        if levels is not None and vgrid not in _LEVEL_FIELDS:
            raise ValueError(
                f"levels= needs vgrid= one of {sorted(_LEVEL_FIELDS)} "
                f"(a vertical remap target) -- got vgrid={vgrid!r}")
        if log_level is not None and log_level not in _LOG_LEVELS:
            raise ValueError(
                f"log_level must be one of {_LOG_LEVELS}, got {log_level!r}")
        if names is not None:
            names = _tokens(names)
        self.names = names
        self.every = every
        self.vgrid = vgrid
        self.levels = None if levels is None else list(levels)
        self.filename = filename
        self.precision = precision
        self.mask_vanished = mask_vanished
        self.reproducing_sums = reproducing_sums
        self.deflate = deflate
        self.output_dir = output_dir
        self.status_interval = status_interval
        self.log_level = log_level
        self.enabled = bool(enabled)

    def apply(self, config):
        if self.output_dir is not None:
            config.output.output_dir = self.output_dir
        if self.status_interval is not None:
            config.logging.status_interval = float(self.status_interval)
        if self.log_level is not None:
            config.logging.log_level = self.log_level
        config.ocean_diag.enabled = self.enabled
        config.output.output_to_file = self.enabled
        if not self.enabled:
            return
        if self.names is not None:
            config.ocean_diag.diags = self.names
        if self.every is not None:
            config.ocean_diag.dt_out = float(self.every)
        if self.vgrid is not None:
            config.ocean_diag.vgrid = self.vgrid
        if self.filename is not None:
            config.ocean_diag.filename = self.filename
        if self.precision is not None:
            config.ocean_diag.output_precision = self.precision
        if self.mask_vanished is not None:
            config.ocean_diag.mask_vanished_layers = bool(self.mask_vanished)
        if self.reproducing_sums is not None:
            config.ocean_diag.reproducing_sums = bool(self.reproducing_sums)
        if self.levels is not None:
            arr_name, n_name = _LEVEL_FIELDS[self.vgrid]
            setattr(config.ocean_diag, arr_name, self.levels)
            setattr(config.ocean_diag, n_name, len(self.levels))
        if self.deflate is not None:
            config.output.compress_output = bool(self.deflate)
            config.output.compress_level = max(1, int(self.deflate))


class Duration:
    """``Duration(t_end, time_unit="s")`` -> ``&time_nml t_end`` +
    ``time_unit``. Pass as ``Model(..., duration=Duration(...))``.

    `t_end`/`time_unit` were previously reachable ONLY through the
    escape hatch (`model.config.time.t_end = ...`). :meth:`Model.run`
    drives the outer loop directly in raw simulated seconds via
    ``until=``/``steps=`` and never reads either knob -- but `t_end`
    still gates ``create()``: the Fortran's own cross-check
    (`rdb_config.F90:3252`) refuses a ``Diagnostics(every=...)``
    cadence (``&ocean_diag_nml dt_out``) that exceeds `t_end` once BOTH
    are scaled through `time_unit`, and `t_end` DEFAULTS to 1.0
    **second** -- almost any real diagnostic cadence exceeds that, so
    "turn diagnostics on, get no curated route to raise t_end" was a
    first-run stumbling block.

    `time_unit` ALSO scales `status_interval`
    (``Diagnostics(status_interval=...)``) the same way
    (`apply_time_unit_cascade`, `rdb_config.F90:2886`, multiplies
    `t_end`/`status_interval`/`ocean_diag.dt_out` through by the SAME
    factor) -- set `time_unit=` here to whatever unit every other
    time-ish curated value in this `Model(...)` is meant to be read in.
    """

    def __init__(self, t_end, time_unit="s"):
        if time_unit not in _TIME_UNITS:
            raise ValueError(
                f"time_unit must be one of {_TIME_UNITS!r} (mirrors "
                f"apply_time_unit_cascade's own case list, "
                f"rdb_config.F90:2886), got {time_unit!r}")
        self.t_end, self.time_unit = t_end, time_unit

    def apply(self, config):
        config.time.t_end = float(self.t_end)
        config.time.time_unit = self.time_unit


class Restart:
    """``Restart(every=, file=)`` -> ``&output_nml restart_interval`` /
    ``restart_file``.

    Writing a restart from a live handle has no C-ABI entry point yet
    (:meth:`~rdb._model.Model.write_restart` raises, naming the
    same gap) -- this object still lets a curated ``Model`` schedule
    periodic NAMELIST-DRIVEN restart writes without dropping to the
    escape hatch (``model.config.output.restart_interval = ...``).
    """

    def __init__(self, every=None, file=None):
        self.every = every
        self.file = file

    def apply(self, config):
        if self.every is not None:
            config.output.restart_interval = float(self.every)
        if self.file is not None:
            config.output.restart_file = self.file
