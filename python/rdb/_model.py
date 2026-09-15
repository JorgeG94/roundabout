"""The Model: lifecycle, single-handle guard, field table.

Config is a raw namelist string OR a :class:`rdb._config.Config` (P4,
the generated typed layer -- ``model.config.ocean_epbl.mstar = 1.2``).
Both go through the SAME path: a ``Config`` is serialised by
``to_namelist()`` into exactly the namelist text a hand-written string
would have been, then handed to ``rdb_ocean_create_from_string`` --
there is no second, mirrored entry point. Two-phase create (P2.5) is
exposed directly: ``Model(nml_text_or_config)`` is the one-shot path
(``rdb_ocean_create_from_string``); ``Model.pending(nml_text_or_config)``
returns a handle in the PENDING window so ``stage_bathymetry`` /
``stage_metrics`` / ``stage_topology`` can inject geometry before
``finalize()`` runs ``rdb_ocean_create_finalize``.
"""

from __future__ import annotations

import ctypes
import weakref
from pathlib import Path

from . import _compose, _ffi
from ._config import Config
from ._errors import (AlreadyExistsError, NotInitialisedError,
                       RdbReadOnlyError, RdbUnsupportedError, check)
from ._field import DerivedField, DiagnosticField, Field, Loc
from ._nml_parse import (config_from_namelist_text, merge_explicit,
                          resolve_relative_file_paths)

C, F = Loc.CENTER, Loc.FACE


def _resolve_namelist_text(namelist_text_or_config):
    """Accept a raw namelist string OR a `Config` (P4); return the exact
    text that crosses the ABI, and the `Config` object to remember as
    `model.config` (a fresh, empty one when a raw string was given -- see
    `Model.config`'s docstring for why)."""
    if isinstance(namelist_text_or_config, Config):
        return namelist_text_or_config.to_namelist(), namelist_text_or_config
    if isinstance(namelist_text_or_config, str):
        return namelist_text_or_config, None
    raise TypeError(
        "Model(...) / Model.pending(...) accept a namelist string or a "
        f"rdb.Config, got {type(namelist_text_or_config).__name__}")

_live_model_ref = None  # module-global single-live-handle guard (mirrors
# the Fortran g_handle_live invariant, but checked client-side so the
# error can name the model still open -- info Fortran cannot have).


def _check_no_live_model():
    global _live_model_ref
    if _live_model_ref is not None:
        existing = _live_model_ref()
        if existing is not None:
            raise AlreadyExistsError(
                "a Roundabout ocean Model is already open in this process "
                "(multi-instance is not supported -- see CLAUDE.md / "
                "docs/ocean_python_api_plan.md S2.1); close() or destroy "
                "it before creating another"
            )
        _live_model_ref = None  # stale: prior model was GC'd without close()


def _claim_live_model(model):
    global _live_model_ref
    _live_model_ref = weakref.ref(model)


def _release_live_model(model):
    global _live_model_ref
    if _live_model_ref is not None and _live_model_ref() is model:
        _live_model_ref = None


# Field descriptor table: name -> (getter attr, ndim, location, units,
# long_name, setter attr or None). Shapes and staggering follow D3's
# table in 06_python_surface_design.md / sketch/03_field_semantics.py.
_FIELDS_3D = {
    "h": ("rdb_ocean_get_h_layer_ptr", (C, C, C), "m",
          "layer thickness", "rdb_ocean_set_h"),
    "u": ("rdb_ocean_get_u_face_x_layer_ptr", (F, C, C), "m/s",
          "west-face x-velocity", "rdb_ocean_set_u"),
    "v": ("rdb_ocean_get_v_face_y_layer_ptr", (C, F, C), "m/s",
          "south-face y-velocity", "rdb_ocean_set_v"),
    "hu": ("rdb_ocean_get_hu_ptr", (F, C, C), "m^2/s",
           "west-face x transport (h*u)", None),
    "hv": ("rdb_ocean_get_hv_ptr", (C, F, C), "m^2/s",
           "south-face y transport (h*v)", None),
    "w": ("rdb_ocean_get_w_interface_ptr", (C, C, F), "m/s",
          "vertical velocity at layer interfaces", None),
    "rho": ("rdb_ocean_get_rho_layer_ptr", (C, C, C), "kg/m^3",
            "in-situ density", None),
    "kv": ("rdb_ocean_get_kv_ptr", (C, C, F), "m^2/s",
           "vertical viscosity", None),
    "kt": ("rdb_ocean_get_kt_ptr", (C, C, F), "m^2/s",
           "vertical heat diffusivity", None),
    "ks": ("rdb_ocean_get_ks_ptr", (C, C, F), "m^2/s",
           "vertical salt (+ passive tracer) diffusivity", None),
}

_FIELDS_2D = {
    "bathymetry": ("rdb_ocean_get_b_ptr", (C, C), "m",
                   "bathymetry (positive down)", "rdb_ocean_set_bathymetry"),
    "eta": ("rdb_ocean_get_bt_eta_ptr", (C, C), "m",
            "sea-surface height (diagnostic, sum_k h - b)", None),
    "tau_x": ("rdb_ocean_get_tau_x_ptr", (F, C), "N/m^2",
              "east-face wind stress", None),
    "tau_y": ("rdb_ocean_get_tau_y_ptr", (C, F), "N/m^2",
              "north-face wind stress", None),
    "Q_heat": ("rdb_ocean_get_q_heat_ptr", (C, C), "W/m^2",
               "net surface heat flux", "rdb_ocean_set_heat_flux"),
    "Q_salt": ("rdb_ocean_get_q_salt_ptr", (C, C), "",
               "net surface salt flux", "rdb_ocean_set_salt_flux"),
    "wet_mask": ("rdb_ocean_get_wet_t_ptr", (C, C), "",
                 "wet mask at T points (1=wet, 0=land)", None),
}


#: Curated top-level keywords that mean "compose from P5 objects" --
#: naming them explicitly (rather than a generic **kwargs) keeps
#: Model.__init__'s signature self-documenting and lets a typo raise
#: TypeError immediately instead of silently doing nothing.
_CURATED_KWARGS = (
    "grid", "bathymetry", "timestep", "coriolis", "momentum_advection",
    "buoyancy", "pressure", "vcoord", "closures", "boundaries", "forcing",
    "tracers", "initial_condition", "barotropic", "continuity",
    "diagnostics", "restart", "duration",
)


class _DiagnosticsView:
    """``model.diagnostics`` (P7) — introspection over what CAN be
    requested vs what IS registered.

    ``.available`` is the same static catalog as the module-level
    :func:`rdb.available_diagnostics` (canonical + derived names;
    reachable without a live model at all). ``.selected`` is the live,
    actually-registered set on THIS instance — canonical names that
    were gated off (e.g. ``temperature`` with thermodynamics disabled)
    or never mentioned in ``&ocean_diag_nml diags`` do not appear here
    even though they are in ``.available``.
    """

    __slots__ = ("_model",)

    def __init__(self, model):
        self._model = model

    @property
    def available(self):
        return _ffi.available_diagnostics()

    @property
    def selected(self):
        m = self._model
        m._check_alive()
        count = ctypes.c_int()
        check(m._lib.rdb_ocean_get_diag_count(
            m._handle, ctypes.byref(count)), "rdb_ocean_get_diag_count")
        names = []
        buf = ctypes.create_string_buffer(64)
        for i in range(count.value):
            n = m._lib.rdb_ocean_list_diags(m._handle, i, buf, 64)
            names.append(buf.raw[:n].decode("utf-8"))
        return names

    def __repr__(self):
        return f"Diagnostics(available={len(self.available)} names, " \
               f"selected={self.selected!r})"


class Model:
    """A live ocean solver handle.

    Two construction paths:

    - ``Model(namelist_text)`` / ``Model(config=a_Config)`` -- the P3/P4
      one-shot path: creates immediately, matching the original contract.
      For pre-create geometry injection on THIS path use :meth:`pending`
      -> :meth:`stage_bathymetry` / :meth:`stage_metrics` /
      :meth:`stage_topology` -> :meth:`finalize` directly.
    - ``Model(grid=..., bathymetry=..., closures=[...], ...)`` -- P5/P6:
      curated physics objects compose onto a fresh `Config` (D2.0 rule 1).
      Creation is DEFERRED to :meth:`create` / the first ``with model:``
      so the escape hatch (``model.config.<group>.<knob> = ...``) can
      still edit the composed config -- see the worked example
      (`sketch/01_regional_run.py`): it mutates ``model.config`` and
      prints ``to_namelist()`` between ``Model(...)`` and ``with
      model:``. Pass ``defer=True`` to require an EXPLICIT
      :meth:`create` call even inside ``with`` (mirrors
      `Model.from_namelist`'s same knob).

    **P7 update:** a curated ``Diagnostics``/``Restart`` object now
    exists (``rdb.diagnostics``) -- pass ``diagnostics=Diagnostics(...)``
    (or ``diagnostics=Diagnostics(enabled=False)`` to turn output off
    entirely) as a curated keyword. With NO ``diagnostics=`` given at
    all, ``&ocean_diag_nml``/``&output_nml`` stay at the Fortran's OWN
    defaults, which write a relative-path NetCDF file
    (``./output/ocean_diag_rank_000000.nc``) -- a missing ``./output/``
    directory on a first run no longer ABORTS the interpreter (P7 F1:
    ``nc_create_file`` now threads an ``ierr`` and the directory is
    created if missing), but a genuinely unwritable/locked path still
    raises ``RdbIOError``. The escape hatch
    (``model.config.ocean_diag.enabled = False``) still works and is
    what ``Diagnostics(enabled=False).apply()`` does under the hood.

    **P8 update:** pass ``duration=Duration(t_end=..., time_unit=...)``
    (``rdb.diagnostics``) to reach ``&time_nml t_end``/``time_unit``
    -- previously escape-hatch only. ``t_end`` does not bound
    :meth:`run`'s own loop (``until=``/``steps=`` do that directly),
    but it DOES gate ``create()``: the Fortran cross-check refuses a
    ``Diagnostics(every=...)`` cadence that exceeds ``t_end`` once both
    are scaled through ``time_unit``, and ``t_end`` defaults to 1.0
    SECOND -- almost any real diagnostic cadence exceeds that, so
    diagnostics + no ``duration=`` was a first-run trap.
    """

    def __init__(self, namelist_text: str | Config | None = None, *,
                 grid=None, bathymetry=None, timestep=None, coriolis=None,
                 momentum_advection=None, buoyancy=None, pressure=None,
                 vcoord=None, closures=(), boundaries=None, forcing=(),
                 tracers=(), initial_condition=None, barotropic=None,
                 continuity=None, diagnostics=None, restart=None,
                 duration=None, config=None, defer=False, z=None):
        if z is not None:
            raise RdbUnsupportedError(
                "z= is not supported: Roundabout has no explicit vertical-"
                "grid input in any form (no interface list) -- only "
                "nz_layers + vcoord_type + scalar family parameters. Use "
                "nz= (on the grid's size=) and vcoord= instead. See D6.4.")
        self._bootstrap()

        _frame_locals = locals()
        curated_given = any(
            _frame_locals[k] not in (None, (), []) for k in _CURATED_KWARGS)
        if curated_given:
            if namelist_text is not None or config is not None:
                raise TypeError(
                    "Model(...) accepts EITHER a namelist string / "
                    "Config (the positional argument / config=) OR "
                    "curated keyword arguments (grid=, bathymetry=, "
                    "closures=, ...), not both -- build a Config with "
                    "the curated objects, or use the escape hatch "
                    "(model.config.<group>.<knob> = ...) after "
                    "construction, but not a mix of both entry points")
            composed = _compose.compose(
                grid=grid, bathymetry=bathymetry, timestep=timestep,
                coriolis=coriolis, momentum_advection=momentum_advection,
                buoyancy=buoyancy, pressure=pressure, vcoord=vcoord,
                closures=closures, boundaries=boundaries, forcing=forcing,
                tracers=tracers, initial_condition=initial_condition,
                barotropic=barotropic, continuity=continuity,
                diagnostics=diagnostics, restart=restart,
                duration=duration)
            self._adopt_composed(composed, defer=defer)
            return

        src = config if config is not None else namelist_text
        if src is not None:
            self._create_from_string(src)

    def _bootstrap(self):
        """Shared base state for every construction path (`__init__`,
        `pending`, `from_namelist`)."""
        self._lib = _ffi.get_lib()
        self._wp_bytes = self._lib.wp_bytes
        self._wp_ctype = self._lib.wp_ctype
        self._wp_typestr = self._lib.wp_typestr
        self._handle = ctypes.c_void_p()
        self._epoch = 0
        self._closed = True
        self._pending = False
        self._uncreated = False
        self._auto_create_on_enter = True
        self._nx = self._ny = self._nz = self._nghost = None
        self._config = None
        self._composed = None
        self._flux_components_enabled = False

    def _adopt_composed(self, composed, *, defer):
        """Record a `_compose.Composed` result WITHOUT creating the
        handle yet (P5/P6: creation defers to :meth:`create` / the first
        ``with model:`` so the escape hatch can still edit
        `model.config`)."""
        self._config = composed.config
        self._composed = composed
        self._flux_components_enabled = composed.flux_components_enabled
        self._uncreated = True
        self._auto_create_on_enter = not defer

    def create(self):
        """Perform the deferred `create()` for a curated
        (``Model(grid=...)``) or ``defer=True`` construction. Idempotent
        error if this handle is not in the uncreated/deferred state
        (already created, or built from a raw namelist/Config, which
        creates eagerly in `__init__`)."""
        if not self._uncreated:
            raise NotInitialisedError(
                "create() called on a Model that is not in the deferred "
                "(uncreated) state -- either it is already created, or "
                "it was built from a raw namelist string / Config, which "
                "creates eagerly in __init__ and has no separate "
                "create() step")
        self._create_from_composed(self._composed)
        self._uncreated = False

    def _create_from_composed(self, composed):
        if composed.needs_pending:
            text = composed.config.to_namelist()
            _check_no_live_model()
            buf = text.encode("utf-8")
            status = self._lib.rdb_ocean_create_pending(
                buf, len(buf), ctypes.byref(self._handle))
            check(status, "rdb_ocean_create_pending")
            self._closed = False
            self._pending = True
            _claim_live_model(self)
            try:
                if composed.bathymetry_stage is not None:
                    depth, convention = composed.bathymetry_stage
                    self.stage_bathymetry(depth, convention)
                if composed.periodic_x or composed.periodic_y:
                    self.stage_topology(composed.periodic_x,
                                         composed.periodic_y)
                self.finalize()
            except Exception:
                self.close()
                raise
        else:
            self._create_from_string(composed.config)

    @classmethod
    def from_namelist(cls, path, *, defer=False, **overrides):
        """Load a namelist FILE -- read in Python, with relative paths
        inside it resolved against ITS OWN directory before anything is
        sent (retires the recovered wrapper's ``os.chdir()``). Same C
        entry point, same validation as ``Model(text)``: the resulting
        ``.config`` is the SAME TYPED `Config` a hand-written script
        would build (not a second representation), so a user migrates
        one knob at a time: ``model.config.ocean_hvisc.nu_h = 5000.0``.

        ``**overrides`` are the SAME curated keyword arguments
        ``Model(...)`` accepts (``grid=``, ``closures=``, ...); composing
        them onto the file's config raises ``ConfigConflictError`` if the
        file and a curated object disagree on the same knob, rather than
        silently letting one win.
        """
        p = Path(path)
        text = p.read_text()
        file_cfg = config_from_namelist_text(text)
        resolve_relative_file_paths(file_cfg, p.parent)

        if overrides:
            composed = _compose.compose(config=Config(), **overrides)
            merge_explicit(file_cfg, composed.config)
            needs_pending = composed.needs_pending
            bathymetry_stage = composed.bathymetry_stage
            periodic_x, periodic_y = composed.periodic_x, composed.periodic_y
            flux_components_enabled = composed.flux_components_enabled
        else:
            needs_pending = False
            bathymetry_stage = None
            periodic_x = periodic_y = False
            flux_components_enabled = False

        final = _compose.Composed(file_cfg, needs_pending, bathymetry_stage,
                                   periodic_x, periodic_y,
                                   flux_components_enabled)
        self = cls.__new__(cls)
        self._bootstrap()
        self._adopt_composed(final, defer=defer)
        return self

    # ---- construction -------------------------------------------------
    def _create_from_string(self, namelist_text_or_config):
        text, config = _resolve_namelist_text(namelist_text_or_config)
        _check_no_live_model()
        buf = text.encode("utf-8")
        status = self._lib.rdb_ocean_create_from_string(
            buf, len(buf), ctypes.byref(self._handle))
        check(status, "rdb_ocean_create_from_string")
        self._closed = False
        self._pending = False
        self._config = config
        _claim_live_model(self)
        self._refresh_grid_info()

    @classmethod
    def pending(cls, namelist_text: str | Config) -> "Model":
        """Phase 1 of 2 (P2.5): build + validate a config and allocate a
        handle without running setup or mapping the device. Call
        :meth:`stage_bathymetry` / :meth:`stage_metrics` /
        :meth:`stage_topology` (any subset, any order), then
        :meth:`finalize`."""
        text, config = _resolve_namelist_text(namelist_text)
        self = cls.__new__(cls)
        self._bootstrap()
        self._config = config

        _check_no_live_model()
        buf = text.encode("utf-8")
        status = self._lib.rdb_ocean_create_pending(
            buf, len(buf), ctypes.byref(self._handle))
        check(status, "rdb_ocean_create_pending")
        self._closed = False
        self._pending = True
        _claim_live_model(self)
        return self

    def stage_bathymetry(self, depth, convention: int):
        """Stage an interior-sized bathymetry array on a pending handle.

        ``depth`` is a nested ``[nx][ny]`` list/tuple (or anything with
        ``__array_interface__``, e.g. a numpy array). ``convention`` is
        REQUIRED -- ``_ffi.BATHY_DEPTH_POSITIVE_DOWN`` (Roundabout's own, a
        4000 m-deep cell is +4000) or ``_ffi.BATHY_HEIGHT_POSITIVE_UP``
        (GEBCO/ETOPO/Oceananigans convention, the same cell is -4000).
        There is no default: D6.2 -- getting this wrong silently makes
        every cell read as land and produces a clean, crash-free,
        entirely wrong quiescent run.
        """
        self._require_pending()
        nx_p, ny_p = _shape_of(depth, ndim=2)
        flat = _ffi.flatten_nested(depth, (nx_p, ny_p))
        arr = _ffi.make_double_array(flat)
        status = self._lib.rdb_ocean_stage_bathymetry(
            self._handle, arr, nx_p, ny_p, int(convention))
        check(status, "rdb_ocean_stage_bathymetry")

    def stage_metrics(self, x, y, dx, dy, area, nxp, nyp, nx, ny):
        """Stage MOM6-style supergrid metric arrays on a pending handle.
        See ``rdb_ocean_stage_metrics`` in ``include/rdb_ocean.h``
        for the exact shapes."""
        self._require_pending()
        x_arr = _ffi.make_double_array(_ffi.flatten_nested(x, (nxp, nyp)))
        y_arr = _ffi.make_double_array(_ffi.flatten_nested(y, (nxp, nyp)))
        dx_arr = _ffi.make_double_array(_ffi.flatten_nested(dx, (nx, nyp)))
        dy_arr = _ffi.make_double_array(_ffi.flatten_nested(dy, (nxp, ny)))
        area_arr = _ffi.make_double_array(_ffi.flatten_nested(area, (nx, ny)))
        status = self._lib.rdb_ocean_stage_metrics(
            self._handle, x_arr, y_arr, dx_arr, dy_arr, area_arr,
            nxp, nyp, nx, ny)
        check(status, "rdb_ocean_stage_metrics")

    def stage_topology(self, periodic_x: bool, periodic_y: bool):
        """Stage per-dimension periodicity on a pending handle."""
        self._require_pending()
        status = self._lib.rdb_ocean_stage_topology(
            self._handle, int(bool(periodic_x)), int(bool(periodic_y)))
        check(status, "rdb_ocean_stage_topology")

    def finalize(self):
        """Phase 2 of 2: complete a handle started by :meth:`pending`."""
        self._require_pending()
        status = self._lib.rdb_ocean_create_finalize(
            ctypes.byref(self._handle))
        if status != _ffi.OK:
            # F9 contract: the whole handle was destroyed on failure.
            self._closed = True
            self._pending = False
            _release_live_model(self)
            self._epoch += 1
        check(status, "rdb_ocean_create_finalize")
        self._pending = False
        self._refresh_grid_info()

    def _require_pending(self):
        if self._closed:
            raise NotInitialisedError(
                "model has been closed; nothing to stage/finalize")
        if not self._pending:
            raise NotInitialisedError(
                "model is not in the pending (staging) window -- call "
                "Model.pending(...) first, or this handle was already "
                "finalized")

    # ---- lifecycle ------------------------------------------------------
    def close(self):
        """Idempotent. A second close() (or closing a never-created
        model) is a clean no-op -- safe for __del__ to call blind."""
        if self._closed:
            return
        status = self._lib.rdb_ocean_destroy(ctypes.byref(self._handle))
        self._closed = True
        self._pending = False
        _release_live_model(self)
        self._epoch += 1  # invalidate every outstanding Field
        check(status, "rdb_ocean_destroy")

    def __enter__(self):
        if self._uncreated and self._auto_create_on_enter:
            self.create()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        # Backstop, not policy (D6.1): __del__ timing is not a resource
        # policy and a GPU context should not be left to the collector,
        # but a leaked handle is worse than a best-effort cleanup here.
        try:
            self.close()
        except Exception:
            pass

    def destroy(self):
        """Alias for :meth:`close` (matches the C ABI's naming)."""
        self.close()

    def _check_alive(self):
        if self._uncreated:
            raise NotInitialisedError(
                "model has not been created yet -- this Model was built "
                "from curated keyword arguments (or defer=True), which "
                "defers create() to the first `with model:` or an "
                "explicit model.create() call so the escape hatch "
                "(model.config.<group>.<knob> = ...) can still edit the "
                "composed config beforehand")
        if self._closed:
            raise NotInitialisedError("model has been closed")

    # ---- grid / introspection --------------------------------------------
    def _refresh_grid_info(self):
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        nz = ctypes.c_int()
        ng = ctypes.c_int()
        status = self._lib.rdb_ocean_get_grid_info(
            self._handle, ctypes.byref(nx), ctypes.byref(ny),
            ctypes.byref(nz), ctypes.byref(ng))
        check(status, "rdb_ocean_get_grid_info")
        self._nx, self._ny, self._nz, self._nghost = (
            nx.value, ny.value, nz.value, ng.value)

    @property
    def config(self) -> Config:
        """The P4 typed config (`rdb.Config`) this model was built
        from -- `model.config.ocean_epbl.mstar` etc. Post-create this is a
        RECORD, not a live handle: mutating it further does not reach a
        running solver (there is no reconfigure-in-place entry point).

        When the model was constructed from a raw namelist string (or
        `Model.from_namelist` in a later phase) rather than a `Config`,
        this lazily returns a fresh, ALL-UNSET `Config` -- it does not
        parse the string back into typed knobs, so it will not reflect
        what was actually set. Read `model.grid_info` / the field
        accessors for the ground truth on a string-built model.
        """
        if self._config is None:
            self._config = Config()
        return self._config

    @property
    def nghost(self):
        return self._nghost

    @property
    def grid_info(self):
        """``{'nx', 'ny', 'nz', 'nghost'}`` -- physical (interior) shape
        plus ghost width."""
        self._check_alive()
        return {"nx": self._nx, "ny": self._ny, "nz": self._nz,
                "nghost": self._nghost}

    @property
    def time(self) -> float:
        self._check_alive()
        t = ctypes.c_double()
        check(self._lib.rdb_ocean_get_time(self._handle, ctypes.byref(t)),
              "rdb_ocean_get_time")
        return t.value

    def _step_count(self) -> int:
        self._check_alive()
        n = ctypes.c_int()
        check(self._lib.rdb_ocean_get_step_count(
            self._handle, ctypes.byref(n)), "rdb_ocean_get_step_count")
        return n.value

    @property
    def step_count(self) -> int:
        return self._step_count()

    @property
    def total_mass(self) -> float:
        self._check_alive()
        m = ctypes.c_double()
        check(self._lib.rdb_ocean_get_total_mass(
            self._handle, ctypes.byref(m)), "rdb_ocean_get_total_mass")
        return m.value

    @property
    def kinetic_energy(self) -> float:
        self._check_alive()
        ke = ctypes.c_double()
        check(self._lib.rdb_ocean_get_kinetic_energy(
            self._handle, ctypes.byref(ke)), "rdb_ocean_get_kinetic_energy")
        return ke.value

    # ---- stepping ---------------------------------------------------------
    def step(self, n: int = 1):
        """Advance ``n`` fixed-dt outer steps."""
        self._check_alive()
        check(self._lib.rdb_ocean_step(self._handle, int(n)),
              "rdb_ocean_step")

    def run(self, *, until=None, steps=None, callback=None, every=None):
        """Advance the model. Exactly one of ``until=`` (simulated
        seconds) or ``steps=`` (outer steps) must be given.
        ``callback(model)`` fires every ``every`` -- seconds of simulated
        time under ``until=``, outer steps under ``steps=``; omit
        ``every`` to fire on every outer step.

        NO ``sync()`` CALL ANYWHERE IN THIS LOOP: ``step()`` only
        invalidates the generation counter; every `Field` the callback
        (or the caller, after `run()` returns) touches re-syncs lazily
        on that access (D3.4). This method steps one outer step at a
        time -- simple and obviously correct rather than batched, which
        is the right tradeoff for a callback-driven loop; a
        ``callback=None`` run over a very large ``steps=``/``until=`` may
        prefer a direct :meth:`step` call for fewer ctypes round-trips.
        """
        self._check_alive()
        if (until is None) == (steps is None):
            raise ValueError("run(...) needs exactly one of until= or "
                              "steps=")
        if steps is not None:
            every_steps = 1 if every is None else int(every)
            for i in range(1, int(steps) + 1):
                self.step(1)
                if callback is not None and i % every_steps == 0:
                    callback(self)
            return
        next_fire = self.time + every if every is not None else None
        while self.time < until - 1e-9:
            self.step(1)
            if callback is None:
                continue
            if every is None:
                callback(self)
            elif self.time >= next_fire - 1e-9:
                callback(self)
                next_fire += every

    def write_diagnostics(self):
        """Fire the registered diagnostic set now, out of cadence.

        NOT REACHABLE on the current C ABI: ``rdb_ocean_step`` already
        drives the diag manager's cadence dispatch internally (device-
        side fill, H<-D on the cadence fire, per ``&ocean_diag_nml
        dt_out``), but ``include/rdb_ocean.h`` exposes no call to
        force an OUT-of-cadence write. This is a P2/P2.5 Fortran-ABI
        gap, not something P6 (pure Python, no Fortran changes) can add
        -- raises rather than silently doing nothing.
        """
        self._check_alive()
        raise RdbUnsupportedError(
            "write_diagnostics() has no reachable entry point: the diag "
            "manager's cadence dispatch already runs inside "
            "rdb_ocean_step (driven by &ocean_diag_nml dt_out), but "
            "include/rdb_ocean.h exposes no call to force an out-of-"
            "cadence write today. Set dt_out small enough that the "
            "cadence itself covers what you need, or extend the C ABI "
            "with a force-write entry point.")

    def write_restart(self, path=None):
        """NOT REACHABLE on the current C ABI. ``rdb_ocean_api.F90``
        says so explicitly in its own comment: 'no restart on the C ABI
        yet (a future phase's concern, not silently dropped --
        engine_setup's restart_file is simply [unused by the handle
        path])'. Restart is namelist-driven only today (`&output_nml
        restart_file` / `restart_interval`) -- set those on
        `model.config.output` before `create()`.
        """
        self._check_alive()
        raise RdbUnsupportedError(
            "write_restart() has no reachable entry point: the C ABI "
            "has no restart call yet (rdb_ocean_api.F90's own comment: "
            "'no restart on the C ABI yet -- a future phase's concern, "
            "not silently dropped'). Restart is namelist-driven only "
            "today -- set model.config.output.restart_file / "
            "restart_interval before create().")

    # ---- fields -------------------------------------------------------
    def _get_field(self, name, getter_name, ndim, location, units,
                    long_name, setter_name):
        self._check_alive()
        getter = getattr(self._lib, getter_name)
        ptr = ctypes.c_void_p()
        dims = [ctypes.c_int() for _ in range(ndim)]
        gen = ctypes.c_int()
        args = [self._handle, ctypes.byref(ptr)] + \
            [ctypes.byref(d) for d in dims] + [ctypes.byref(gen)]
        status = getter(*args)
        check(status, getter_name)
        full_shape = tuple(d.value for d in dims)
        setter = None
        if name in ("Q_heat", "Q_salt") and self._flux_components_enabled:
            def setter(handle, arr, *shape, _name=name):
                raise RdbReadOnlyError(
                    f"{_name!r} is DERIVED while "
                    f"SurfaceFluxComponents() is enabled -- "
                    f"ocean_surface_flux_assemble computes it every "
                    f"thermo step from the component fluxes, so a "
                    f"direct write would be silently overwritten at the "
                    f"next thermo step (src/core/ocean/README.md:138). "
                    f"Set the component fluxes instead, or drop "
                    f"SurfaceFluxComponents() from forcing=[...].")
        elif setter_name is not None:
            fn = getattr(self._lib, setter_name)

            def setter(handle, arr, *shape, _fn=fn):
                return _fn(handle, arr, *shape)

        return Field(self, name, ptr.value, full_shape, gen.value, location,
                     units=units, long_name=long_name, setter=setter)

    def __getattr__(self, name):
        # Declarative field table -- __getattr__ only fires when the
        # attribute wasn't found normally, so this never shadows real
        # methods/properties above.
        if name in _FIELDS_3D:
            getter_name, loc, units, long_name, setter_name = _FIELDS_3D[name]
            return self._get_field(name, getter_name, 3, loc, units,
                                    long_name, setter_name)
        if name in _FIELDS_2D:
            getter_name, loc, units, long_name, setter_name = _FIELDS_2D[name]
            return self._get_field(name, getter_name, 2, loc, units,
                                    long_name, setter_name)
        raise AttributeError(
            f"{type(self).__name__!r} object has no attribute {name!r}")

    @property
    def ssh(self):
        """Sea-surface height, DERIVED as sum_k h - b (same definition the
        diag manager uses). eta is the raw bt_work%bt_eta diagnostic;
        prefer this for "is the free surface where I expect".
        """
        h = self.h.copy()
        b = self.bathymetry.copy()
        nx, ny, nz = self.h.shape
        out = [[sum(h[i][j][k] for k in range(nz)) - b[i][j]
                for j in range(ny)] for i in range(nx)]
        return out

    writable_fields = frozenset({
        "h", "u", "v", "bathymetry", "Q_heat", "Q_salt",
        "temperature", "salinity",
    })

    def set_wind(self, tau_x, tau_y):
        """Overwrite wind stress. Two arrays together -- the ABI has no
        single-component wind setter, so `tau_x`/`tau_y` are read-only
        Fields and this is the writer."""
        self._check_alive()
        nx_p, ny_p = self._nx, self._ny
        taux_flat = _ffi.flatten_nested(tau_x, (nx_p + 1, ny_p))
        tauy_flat = _ffi.flatten_nested(tau_y, (nx_p, ny_p + 1))
        taux_arr = _ffi.make_double_array(taux_flat)
        tauy_arr = _ffi.make_double_array(tauy_flat)
        status = self._lib.rdb_ocean_set_wind(
            self._handle, taux_arr, tauy_arr, nx_p, ny_p)
        check(status, "rdb_ocean_set_wind")

    # ---- tracers --------------------------------------------------------
    @property
    def tracer_names(self):
        self._check_alive()
        count = ctypes.c_int()
        check(self._lib.rdb_ocean_get_tracer_count(
            self._handle, ctypes.byref(count)), "rdb_ocean_get_tracer_count")
        names = []
        buf = ctypes.create_string_buffer(64)
        for i in range(count.value):
            n = self._lib.rdb_ocean_list_tracers(self._handle, i, buf, 64)
            names.append(buf.raw[:n].decode("utf-8"))
        return names

    def tracer(self, name: str) -> DerivedField:
        """Concentration (its own units: degC/PSU/...) for the tracer
        `name`, by NAME -- never index. Raises TracerNotFoundError if no
        such tracer is registered."""
        htr_field = self._raw_tracer_field(name)
        h_field = self.h
        units = "degC" if name == "temperature" else (
            "PSU" if name == "salinity" else "")

        def setter(handle, arr, *shape, _name=name):
            enc = _name.encode("utf-8")
            return self._lib.rdb_ocean_set_tracer(
                handle, enc, len(enc), arr, *shape)

        return DerivedField(self, name, htr_field, h_field, units=units,
                             long_name=f"tracer {name!r} concentration",
                             setter=setter)

    def _raw_tracer_field(self, name: str) -> Field:
        self._check_alive()
        enc = name.encode("utf-8")
        ptr = ctypes.c_void_p()
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        nz = ctypes.c_int()
        gen = ctypes.c_int()
        status = self._lib.rdb_ocean_get_tracer_ptr(
            self._handle, enc, len(enc), ctypes.byref(ptr),
            ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(nz),
            ctypes.byref(gen))
        check(status, "rdb_ocean_get_tracer_ptr")
        return Field(self, name + ".hTr", ptr.value,
                     (nx.value, ny.value, nz.value), gen.value, (C, C, C))

    # ---- diagnostics (P7) -----------------------------------------------
    @property
    def diagnostics(self) -> _DiagnosticsView:
        """``.available`` (static catalog) vs ``.selected`` (what THIS
        instance registered). See :func:`rdb.available_diagnostics`
        for the module-level, no-model-required equivalent of
        ``.available``."""
        return _DiagnosticsView(self)

    def diagnostic(self, name: str) -> DiagnosticField:
        """Read-only in-memory view onto the diag manager's own
        ``output_buffer`` for the REGISTERED diagnostic ``name`` — no
        NetCDF round-trip. This is most of the point of driving the
        model in-process: otherwise you write a file and read it back
        to get what memory already holds.

        Shape is ``(nx, ny, nz)`` (interior, per :class:`Field`'s usual
        convention) — ``nz`` is the layer count for a layered diagnostic
        at the default LAYER output grid, ``1`` for a 2D diagnostic
        (e.g. ``"SSH"``), or the remapped level count when the
        diagnostic's own ``output_vgrid`` is not LAYER. The generation
        is the diagnostic's own per-fire counter (see
        :class:`~rdb._field.DiagnosticField`) — reading before the
        diagnostic has ever fired returns its allocation-time value
        (0.0 for every canonical/derived diagnostic here, matching the
        device-mapped initial value).

        Raises :class:`~rdb._errors.TracerNotFoundError`
        (``OCEAN_STATUS_ERR_NOT_FOUND``) if ``name`` is not CURRENTLY
        registered on this instance — a name may be legal on
        ``diagnostics.available`` yet gated off or unselected; check
        ``diagnostics.selected`` to see what is actually live.
        """
        self._check_alive()
        enc = name.encode("utf-8")
        ptr = ctypes.c_void_p()
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        nz = ctypes.c_int()
        gen = ctypes.c_int()
        status = self._lib.rdb_ocean_get_diagnostic_ptr(
            self._handle, enc, len(enc), ctypes.byref(ptr),
            ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(nz),
            ctypes.byref(gen))
        check(status, "rdb_ocean_get_diagnostic_ptr")
        return DiagnosticField(self, name, ptr.value,
                               (nx.value, ny.value, nz.value), gen.value,
                               (C, C, C), long_name=f"diagnostic {name!r}")

    @property
    def temperature(self) -> DerivedField:
        """Layer temperature (degC), DERIVED from hT/h (D3.2)."""
        return self.tracer("temperature")

    @property
    def salinity(self) -> DerivedField:
        """Layer salinity (PSU), DERIVED from hS/h (D3.2)."""
        return self.tracer("salinity")

    def __repr__(self):
        if self._uncreated:
            state = "uncreated"
        else:
            state = "closed" if self._closed else (
                "pending" if self._pending else "live")
        grid = (f"{self._nx}x{self._ny}x{self._nz}"
                if self._nx is not None else "?")
        return f"Model({state}, grid={grid})"


def _shape_of(nested, ndim):
    shape = []
    v = nested
    for _ in range(ndim):
        shape.append(len(v))
        v = v[0]
    return tuple(shape)
