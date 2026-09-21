#!/usr/bin/env python3
"""Drift gate for docs/CLOSURE_MATRIX.md.

Fails (exit 1) when the closure documentation has drifted from the code:

  1. Every closure-selector knob in src/core/rdb_config.F90 must appear in
     docs/CLOSURE_MATRIX.md.  The `vmix_use_*` family is *discovered* from the
     config, so a new vertical-mixing scheme cannot be added without a matrix
     row.  The other closure/forcing toggles are checked by name.
  2. Every `test_*` named in the matrix must be registered in
     tests/CMakeLists.txt.
  3. The knob-driven ALE-remap default in the config (`remap_method = "ppm"`)
     must match the default the matrix states.  This guards the Tracer
     Advection & Reconstruction section the same way the vmix checks guard the
     vertical-mixing section.
     The HARDCODED ocean rows (PPM horizontal advect, upwind-in-z) have no
     config symbol and are intentionally NOT auto-checked here — see the
     comment in check_tracer_defaults().

Wire into pre-commit + CI.  The matrix exists because prose drifts silently:
CAPABILITIES_AND_LIMITATIONS.md still claimed PP81 was the only vertical mixing
scheme three closures after KPP landed.  This makes that fail the build.

Run from anywhere; paths are resolved relative to this file.
"""
import re
import sys
from pathlib import Path
from typing import List

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "src" / "core" / "rdb_config.F90"
MATRIX = ROOT / "docs" / "CLOSURE_MATRIX.md"
TESTS_CMAKE = ROOT / "tests" / "CMakeLists.txt"

# Named closure / forcing selector knobs (besides the auto-discovered
# vmix_use_* family).  Each must exist in the config AND have a matrix cell.
CURATED_KNOBS = [
    "ocean_hvisc_nml",
    "ocean_bdrag_nml",
    "ocean_kappa_shear_nml",
    # Static ice-shelf cavity geometry + the barotropic datum that absorbs
    # it.  Listed here because the datum silently changes what `bt_eta`
    # MEANS (deviation from the loaded equilibrium, not from z = 0), which
    # is exactly the kind of thing prose forgets.
    "ocean_cavity_dyn_nml",
    # Ice-shelf basal-melt thermodynamics.  Listed because the capability
    # is delivered as a VIRTUAL salt flux with no mass, and because the
    # v1 gaps (no cover mask on atmospheric forcing, no top drag) are
    # configure REFUSALS — exactly the kind of scoping prose forgets.
    "ocean_cavity_melt_nml",
    # Ice-shelf TOP drag.  Listed because the FACE cover rule (OR, so the
    # calving-front face is dragged) and the ONE-C_d agreement rule with
    # &ocean_cavity_melt_nml cdrag_top are both deliberate choices that a
    # reader cannot recover from the code without being told.
    "ocean_tdrag_nml",
]

CURATED_KEYS = [
    # Individual keys (not whole groups) that select a CLOSURE BEHAVIOUR and
    # so must be findable in the matrix.  A plain numeric knob does not
    # belong here; a switch that changes which physics runs does.
    #
    # `buoyancy_coeffs` decides whether the KPP surface buoyancy flux and
    # the double-diffusion density ratio take alpha/beta from the scalar
    # linear-EOS pair or from the ACTIVE equation of state.  Listed because
    # the answer is invisible in the output (both settings produce a
    # plausible Kd) and because it is a no-op under eos="linear" but a
    # factor-of-several change in an ice-shelf cavity under Wright.
    "buoyancy_coeffs",
    #
    # `zfixed_closed_faces` decides whether a `z_fixed` face whose layer is
    # an inert filler on one side is a thin passage or a z-LEVEL WALL.  It
    # is listed because it is invisible in the output (both settings
    # produce a plausible field), because it changes SEVEN kernels at once
    # (continuity, the BT renormaliser, transport Coriolis, the layer
    # velocity mask, kappa_h, the harmonic viscosity's slip condition and
    # the momentum vdiff coupling) and because its fail-loud exclusion
    # list -- GM / Redi / MLE / biharmonic / stress_tensor / wet-dry /
    # bt_halo -- is exactly the kind of scoping prose forgets.
    "zfixed_closed_faces",
]

# Named sea-ice selector knobs (&ocean_ice_nml).  The shipped subsystem —
# Winton column thermodynamics, multi-category ITD, category transport +
# compress_ice, and C-grid EVP dynamics — is gated by these.  Each distinctive
# token must exist in the config AND appear in the sea-ice section of the matrix
# so the whole subsystem can't drift out of the enabled-closure ground truth
# (the group name + one knob per shipped component/family).
ICE_KNOBS = [
    "ocean_ice_nml",
    "nk_ice",
    "adv_substeps",
    "evp_sub_steps",
    "del_sh_min_scale",
    "tdamp",
]


def check_tracer_defaults(config: str, matrix: str) -> List[str]:
    """Guard the knob-driven ALE-remap default.

    Same mechanism as the vmix default checks: grep the default straight out
    of the config_t declaration, then assert the doc states the same default.
    If the code and the doc disagree, the doc is the bug (per CLAUDE.md: code
    is authority) — fail loudly so it's fixed in the same PR.

    NOT covered here (no config symbol to grep, so a future maintainer must not
    assume they're guarded):
      * ocean HORIZONTAL tracer advection is HARDCODED PPM — no namelist knob
        (tracer_advect_{zonal,meridional}_one_impl in rdb_continuity.F90).
      * ocean VERTICAL (in-z) tracer advection is HARDCODED 1st-order upwind
        (tracer_advect_vertical_one_impl in rdb_ocean_vertical_advection.F90).
    These are documented in the matrix's honesty notes but cannot be
    config-grepped; verifying them needs source inspection, not this gate.
    """
    errors: List[str] = []

    # remap_method — character default in config_t (the `= "value"` form).
    m = re.search(r'character\(len=\d+\)\s*::\s*remap_method\s*=\s*"(\w+)"', config)
    if not m:
        errors.append(
            "remap_method declaration not found in rdb_config.F90 — "
            "renamed/removed in code, update this script."
        )
    else:
        default = m.group(1)
        if default != "ppm":
            errors.append(
                f'remap_method default in rdb_config.F90 is "{default}", '
                'expected "ppm" — update the matrix and this script together.'
            )
        # The doc must state the same default for this knob.
        if f'`"{default}"`' not in matrix:
            errors.append(
                'docs/CLOSURE_MATRIX.md does not state remap_method\'s default '
                f'("{default}") as `"{default}"`.'
            )

    return errors


def main() -> int:
    for path in (CONFIG, MATRIX, TESTS_CMAKE):
        if not path.exists():
            print(f"check_closure_matrix: FAIL — missing {path}", file=sys.stderr)
            return 1

    config = CONFIG.read_text()
    matrix = MATRIX.read_text()
    cmake = TESTS_CMAKE.read_text()

    # 1a. discover vmix_use_* toggles *declared* in config_t (the `= default`
    #     form, so comments and local `logical :: a, b` decls don't match).
    vmix_knobs = sorted(set(re.findall(r"logical\s*::\s*(vmix_use_\w+)\s*=", config)))
    knobs = vmix_knobs + CURATED_KNOBS + CURATED_KEYS + ICE_KNOBS

    errors: List[str] = []

    # the curated names must actually exist in the config (catch a rename here)
    phantom = [k for k in CURATED_KNOBS + CURATED_KEYS + ICE_KNOBS if k not in config]
    if phantom:
        errors.append(
            "curated knob(s) not found in rdb_config.F90 — renamed/removed in "
            f"code, update this script: {phantom}"
        )

    # 1b. every selector knob must have a row in the matrix
    missing_in_matrix = [k for k in knobs if k not in matrix]
    if missing_in_matrix:
        errors.append(
            "closure knob(s) in the config with no row in docs/CLOSURE_MATRIX.md: "
            f"{missing_in_matrix}"
        )

    # 2. every test named in the matrix must be registered
    matrix_tests = sorted(set(re.findall(r"\btest_[a-z0-9_]+", matrix)))
    missing_tests = [t for t in matrix_tests if t not in cmake]
    if missing_tests:
        errors.append(
            "test(s) named in docs/CLOSURE_MATRIX.md but not in "
            f"tests/CMakeLists.txt: {missing_tests}"
        )

    # 3. the knob-driven ALE-remap default must match
    errors.extend(check_tracer_defaults(config, matrix))

    if errors:
        for e in errors:
            print(f"check_closure_matrix: FAIL — {e}", file=sys.stderr)
        return 1

    print(
        f"check_closure_matrix: OK — {len(knobs)} closure knobs "
        f"({len(vmix_knobs)} vmix discovered) + {len(matrix_tests)} tests "
        "reconciled + ALE-remap default (remap_method) verified"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
