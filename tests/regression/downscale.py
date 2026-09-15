#!/usr/bin/env python3
"""Dimensionless-number rules for building a TIER-2 downscaled twin.

Why this file exists
====================
The obvious way to make a validation case cheap is to shrink the grid. That
is also the fastest way to build a test that no longer tests anything: the
resolved physics of an ocean case is set by a handful of DIMENSIONLESS
numbers, and a naive `nx/2, ny/2` changes every one of them at once.

Three real defects found on 2026-09-11 all hid behind exactly this:

  * `acc_channel` shipped a domain 10x too large for its `nu_h`, so the Munk
    layer was under one cell and the western boundary was grid-scale noise
    (`tmp_local_artifacts/global_run/FINDINGS.md`).
  * `eady.nml` shipped a front 4x weaker than its header described under
    `nu_h = 100`, and decayed where theory says it grows (re-baselined
    2026-09-13; the nml header carries the sweep).
  * A thermo-cadence instability was survivable only by viscosity margin --
    the same margin a downscale silently destroys
    (`tmp_local_artifacts/eady_hunt/FINDINGS_windowed_advect.md`).

So a tier-2 twin is not accepted because it is small and fast. It is accepted
because it still satisfies the constraints below, CHECKED, with the numbers
printed. `check_twin()` is run for every tier-2 case before the case runs; a
violating twin is a manifest bug and fails loudly rather than quietly testing
different physics.

The rules
=========

R1  MUNK LAYER >= 2 CELLS.      nu_h >= beta * (MUNK_CELLS * dy)**3
    The Munk boundary-layer width is delta_M = (nu_h/beta)**(1/3). Under two
    cells the western boundary current is not resolved and you get grid-scale
    wall noise (a non-monotonic velocity profile into the wall) that looks
    like physics. At -45 deg S with dy = 0.5 deg this floor is nu_h = 2.2e4.
    Only applies to cases with a beta plane and a meridional wall.

R2  ah_max CLAMPS nu_h.         ah_max >= nu_h (when ah_max is set at all)
    `&ocean_hvisc_nml ah_max` is a magnitude CEILING applied to the assembled
    viscosity. Raising `nu_h` to satisfy R1 without raising `ah_max` is inert
    -- the clamp silently throws the increase away and the twin runs at the
    old, too-small viscosity.

R3  VISCOUS CFL.                nu_h * dt / min(dx,dy)**2 <= VISC_CFL (0.125)
    Explicit Laplacian friction is forward-Euler here. Satisfying R1 by
    raising nu_h on a coarse twin can trip this; satisfying R3 by shrinking
    dx can trip R1. When both cannot hold, either `bound_kh = .true.` (the
    per-cell CFL clamp) must be on, or the case is tier-1 only. Note ah_max
    does NOT help: it bounds magnitude, not the CFL.

R4  DEFORMATION RADIUS.         Rd / dx >= RD_CELLS  (4 by default)
    An eddy-resolving case needs ~4-7 cells per first-baroclinic Rd. Below
    ~2 there are no eddies at any run length, so an "eddy" twin that violates
    this is testing advection of a laminar flow.

R5  AVAILABLE POTENTIAL ENERGY. a case with uniform density cannot make
    eddies at ANY resolution. This is a property of the IC, not the grid, so
    the rule here is a flag: a twin declared `eddying` must declare where its
    APE comes from. (`acc_channel` as shipped has none -- it never would have
    produced eddies however long it ran.)

R6  NGHOST vs SCHEME.           weno5 -> nghost >= 3, weno7 -> nghost >= 4,
    periodic-x -> nghost >= 3. A downscale must not quietly drop nghost to
    save cells.

R7  ROSSBY / DOMAIN ASPECT.     L_domain / Rd >= DOMAIN_RD (8 by default) for
    an eddying case -- the box must hold several eddies, or the twin measures
    the box, not the turbulence.

Nothing here is tuned to make a particular twin pass. Where a case cannot be
downscaled without breaking a rule, the manifest marks it `tier1_only` and
says which rule blocks it -- a valid and expected outcome.

Stdlib only (the repo forbids `pip install`).
"""

import math

# --- physical / numerical constants -----------------------------------------
OMEGA_EARTH = 7.2921e-5     # rad/s
R_EARTH = 6.371e6           # m
G_ACCEL = 9.80665           # m/s^2

MUNK_CELLS = 2.0            # R1: minimum Munk-layer width in cells
VISC_CFL = 0.125            # R3: forward-Euler Laplacian stability margin
RD_CELLS = 4.0              # R4: minimum cells per deformation radius
DOMAIN_RD = 8.0             # R7: minimum deformation radii across the domain

# nghost floors keyed by the scheme knob that demands them (R6).
NGHOST_FLOOR = {
    "weno5": 3, "weno7": 4, "weno3": 2,
    "periodic": 3,
}


def beta_at(lat_deg):
    """Planetary vorticity gradient beta = df/dy at a latitude (1/(m s))."""
    return 2.0 * OMEGA_EARTH * math.cos(math.radians(lat_deg)) / R_EARTH


def f_at(lat_deg):
    """Coriolis parameter at a latitude (1/s)."""
    return 2.0 * OMEGA_EARTH * math.sin(math.radians(lat_deg))


def munk_width(nu_h, beta):
    """Munk western-boundary-layer width delta_M = (nu_h/beta)**(1/3), m."""
    if beta <= 0.0:
        return float("inf")
    return (max(nu_h, 0.0) / beta) ** (1.0 / 3.0)


def munk_nu_floor(beta, dy, cells=MUNK_CELLS):
    """R1: the smallest nu_h that resolves the Munk layer over `cells` cells."""
    return beta * (cells * dy) ** 3


def deformation_radius(g_prime, depth, f):
    """First-baroclinic Rossby radius Rd = sqrt(g' H)/|f|, m.

    `g_prime` is the reduced gravity g*drho/rho0 for the dominant interface
    (use G_ACCEL for the barotropic/external radius).
    """
    if f == 0.0:
        return float("inf")
    return math.sqrt(max(g_prime, 0.0) * max(depth, 0.0)) / abs(f)


def viscous_cfl(nu_h, dt, dx, dy=None):
    """R3: the forward-Euler Laplacian-friction CFL number (dimensionless)."""
    h = min(dx, dy) if dy else dx
    if h <= 0.0:
        return float("inf")
    return nu_h * dt / (h * h)


def gravity_wave_cfl(depth, dt, dx):
    """External gravity-wave CFL sqrt(gH)*dt/dx -- the barotropic substep's
    own limit; informational, since `auto_n_inner` derives n_inner from it."""
    if dx <= 0.0:
        return float("inf")
    return math.sqrt(G_ACCEL * max(depth, 0.0)) * dt / dx


# ---------------------------------------------------------------------------
# Twin validation
# ---------------------------------------------------------------------------
class RuleResult(object):
    """One rule's verdict on one twin: pass/fail plus the numbers behind it."""

    __slots__ = ("rule", "ok", "detail", "waived", "waiver")

    def __init__(self, rule, ok, detail, waived=False, waiver=""):
        self.rule = rule
        self.ok = ok
        self.detail = detail
        self.waived = waived
        self.waiver = waiver

    def __str__(self):
        if self.waived:
            return "  {:<22} WAIVED  {}  [{}]".format(
                self.rule, self.detail, self.waiver)
        return "  {:<22} {}  {}".format(
            self.rule, "ok  " if self.ok else "FAIL", self.detail)


def check_twin(spec):
    """Validate one tier-2 downscale spec. Returns a list of RuleResult.

    `spec` is the manifest's `downscale` block, already merged with the
    parent case's resolved numbers. Recognised keys (all optional -- a rule
    whose inputs are absent is skipped, never silently passed):

        dx, dy        twin grid spacing, m
        dt            twin outer timestep, s
        nu_h          twin constant lateral viscosity, m^2/s
        ah_max        twin viscosity ceiling, m^2/s (R2)
        bound_kh      True if the per-cell CFL clamp is on (waives R3)
        lat_deg       representative latitude for beta/f (R1, R4)
        beta          explicit beta, overrides lat_deg (beta-plane cases)
        f0            explicit f, overrides lat_deg (f-plane cases)
        has_western_boundary   True => R1 applies
        eddying       True => R4 + R7 apply
        g_prime, depth         inputs to Rd (R4, R7)
        nx, ny        twin cell counts (R7)
        nghost        twin halo width (R6)
        schemes       iterable of scheme tokens demanding a halo (R6)
        waivers       {rule: "reason"} -- an EXPLICIT, documented exemption
    """
    out = []
    waivers = spec.get("waivers", {}) or {}

    def emit(rule, ok, detail):
        if rule in waivers:
            out.append(RuleResult(rule, True, detail, True, waivers[rule]))
        else:
            out.append(RuleResult(rule, ok, detail))

    dx = spec.get("dx")
    dy = spec.get("dy", dx)
    dt = spec.get("dt")
    nu_h = spec.get("nu_h")
    beta = spec.get("beta")
    if beta is None and spec.get("lat_deg") is not None:
        beta = beta_at(spec["lat_deg"])
    f = spec.get("f0")
    if f is None and spec.get("lat_deg") is not None:
        f = f_at(spec["lat_deg"])

    # -- R1 Munk layer --------------------------------------------------
    if spec.get("has_western_boundary") and nu_h is not None and dy and beta:
        floor = munk_nu_floor(beta, dy)
        width_cells = munk_width(nu_h, beta) / dy
        emit("R1-munk-layer", nu_h >= floor,
             "nu_h={:.3g} vs floor {:.3g} m2/s (delta_M = {:.2f} cells, "
             "need >= {:.1f})".format(nu_h, floor, width_cells, MUNK_CELLS))

    # -- R2 ah_max clamps nu_h ------------------------------------------
    ah_max = spec.get("ah_max")
    if ah_max is not None and nu_h is not None:
        emit("R2-ah_max-clamp", ah_max >= nu_h,
             "ah_max={:.3g} vs nu_h={:.3g} m2/s (ah_max CLAMPS nu_h; a lower "
             "ceiling makes the nu_h increase inert)".format(ah_max, nu_h))

    # -- R3 viscous CFL --------------------------------------------------
    if nu_h is not None and dt and dx:
        cfl = viscous_cfl(nu_h, dt, dx, dy)
        ok = cfl <= VISC_CFL or bool(spec.get("bound_kh"))
        note = "" if cfl <= VISC_CFL else (
            " (bound_kh ON -> per-cell clamp covers it)"
            if spec.get("bound_kh") else
            " -- raise dx, lower dt/nu_h, or set bound_kh=.true.")
        emit("R3-viscous-cfl", ok,
             "nu_h*dt/dx^2 = {:.4f} vs limit {:.3f}{}".format(
                 cfl, VISC_CFL, note))

    # -- R4/R7 deformation radius ---------------------------------------
    if spec.get("eddying"):
        gp, depth = spec.get("g_prime"), spec.get("depth")
        if gp is not None and depth is not None and f:
            rd = deformation_radius(gp, depth, f)
            if dx:
                emit("R4-Rd-resolution", rd / dx >= RD_CELLS,
                     "Rd = {:.1f} km = {:.2f} cells (need >= {:.1f}; below "
                     "~2 no eddies form at ANY run length)".format(
                         rd / 1e3, rd / dx, RD_CELLS))
            nx, ny = spec.get("nx"), spec.get("ny")
            if nx and ny and dx and dy:
                lmin = min(nx * dx, ny * dy)
                emit("R7-domain-Rd", lmin / rd >= DOMAIN_RD,
                     "smallest domain side = {:.2f} Rd (need >= {:.1f}; a "
                     "smaller box measures the box, not the turbulence)".format(
                         lmin / rd, DOMAIN_RD))
        # -- R5 APE ------------------------------------------------------
        ape = spec.get("ape_source")
        emit("R5-ape-source", bool(ape),
             "APE source = {} (a uniform-density IC yields NO eddies at any "
             "resolution)".format(ape if ape else "NONE DECLARED"))

    # -- R6 nghost vs scheme ---------------------------------------------
    nghost = spec.get("nghost")
    schemes = [s for s in (spec.get("schemes") or []) if s in NGHOST_FLOOR]
    if nghost is not None and schemes:
        need = max(NGHOST_FLOOR[s] for s in schemes)
        who = ",".join(sorted(set(schemes)))
        emit("R6-nghost", nghost >= need,
             "nghost={} vs {} required by [{}]".format(nghost, need, who))

    return out


def twin_ok(results):
    """True when no rule failed (waived rules count as passing)."""
    return all(r.ok for r in results)


def format_report(name, results):
    lines = ["{}:".format(name)]
    if not results:
        lines.append("  (no rule inputs declared -- nothing checked)")
    lines.extend(str(r) for r in results)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Standalone self-check: validate every tier-2 twin in the manifest.
# ---------------------------------------------------------------------------
def main(argv=None):
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import manifest

    bad = 0
    checked = 0
    for case in manifest.STABILITY_CASES:
        t2 = case.get("tier2")
        if not t2 or t2.get("skip"):
            continue
        spec = dict(t2.get("dimensionless") or {})
        if not spec:
            continue
        checked += 1
        res = check_twin(spec)
        if not twin_ok(res):
            bad += 1
            print(format_report(case["name"], res))
            print("")
    print("downscale rule check: {} twins checked, {} violating".format(
        checked, bad))
    return 1 if bad else 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
