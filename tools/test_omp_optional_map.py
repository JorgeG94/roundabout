"""Unit tests for omp_optional_map.py — the absent-OPTIONAL map() post-pass.

Run with: pytest tools/test_omp_optional_map.py

Background: LLVM Flang emits an unconditional implicit map for an `optional`,
EXPLICIT-SHAPE array dummy referenced inside an `!$omp target` region. The map
size comes from the declared bounds, but an absent argument's base address is
NULL, so libomptarget calls `hsa_amd_memory_lock(0x0, size)` and the run dies.
Naming the argument in an explicit `map()` clause takes flang's
presence-guarded path. Confirmed on Frontier (MI250X, amdflang 23.0.0git).

The cases below are the ones where a plausible-looking implementation goes
wrong, rather than a restatement of the happy path:

  * scoping — an optional declared in a SIBLING procedure must not be
    attributed to this one. A file-global scan claimed 81 sites where the
    real number was 28, because sibling procedures declare same-named dummies
    that are NOT optional.
  * `!$omp declare target` is a device-routine declaration, not a region, and
    a naive `!$omp target` prefix match swallows it.
  * declarations split across `&` continuations must be seen whole, or the
    `optional` attribute and the entity list end up on different "lines".
  * idempotence — regen re-runs the pass over an already-processed tree.
"""
import importlib.util
import pathlib
import sys

_HERE = pathlib.Path(__file__).parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    sys.path.insert(0, str(_HERE))
    spec.loader.exec_module(mod)
    return mod


m = _load("omp_optional_map")


def _xform(src):
    lines, n = m.transform_lines(src.splitlines())
    return "\n".join(lines), n


def _proc(body, name="kern", decls=""):
    return (f"  subroutine {name}(a, opt)\n"
            f"{decls}"
            f"{body}"
            f"  end subroutine {name}\n")


# --------------------------------------------------------------------------
# The core behaviour
# --------------------------------------------------------------------------

def test_explicit_shape_optional_mapped_by_intent():
    src = _proc(
        decls=("    real, intent(inout) :: a(nx, ny)\n"
               "    real, intent(in), optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1
    assert "map(to: opt)" in out, out


def test_intent_inout_gets_tofrom():
    src = _proc(
        decls=("    integer, intent(in) :: nx, ny\n"
               "    real, intent(inout), optional :: budget(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      budget(i,j) = 0.0\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1
    assert "map(tofrom: budget)" in out, out


def test_unreferenced_optional_not_mapped():
    """Mapping something the region never touches is pure noise."""
    src = _proc(
        decls=("    real, intent(inout) :: a(nx, ny)\n"
               "    real, intent(in), optional :: unused(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = 0.0\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0
    assert "map(" not in out


def test_non_optional_dummy_not_mapped():
    src = _proc(
        decls=("    real, intent(in) :: plain(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = plain(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0


# --------------------------------------------------------------------------
# Shapes flang already handles: must NOT be mapped
# --------------------------------------------------------------------------

def test_assumed_shape_optional_skipped():
    """Descriptor-passed: flang null-checks the box, so the bug can't occur."""
    src = _proc(
        decls=("    real, intent(in), optional :: opt(:,:)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0, out


def test_assumed_size_optional_skipped():
    src = _proc(
        decls=("    real, intent(in), optional :: opt(*)\n"),
        body=("    !$omp target teams distribute parallel do\n"
              "    do i = 1, nx\n"
              "      a(i,1) = opt(i)\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0, out


def test_scalar_optional_skipped():
    """A scalar optional is implicitly firstprivate, not mapped."""
    src = _proc(
        decls=("    real, intent(in), optional :: s\n"),
        body=("    !$omp target teams distribute parallel do\n"
              "    do i = 1, nx\n"
              "      a(i,1) = s\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0, out


# --------------------------------------------------------------------------
# The scoping trap: 81 claimed sites vs 28 real ones
# --------------------------------------------------------------------------

def test_optional_in_sibling_procedure_not_attributed():
    """`res` is optional in helper() but a PLAIN dummy in kern().

    A file-global scan maps it in kern() too. That is the exact mistake that
    inflated a hand-scan from 28 sites to 81.
    """
    src = (
        "module mod_x\n"
        "contains\n"
        + _proc(name="helper",
                decls=("    real, intent(in), optional :: res(nx, ny)\n"),
                body="    x = 1\n")
        + _proc(name="kern",
                decls=("    real, intent(in) :: res(nx, ny)\n"),
                body=("    !$omp target teams distribute parallel do collapse(2)\n"
                      "    do j = 1, ny\n"
                      "    do i = 1, nx\n"
                      "      a(i,j) = res(i,j)\n"
                      "    end do\n"
                      "    end do\n"))
        + "end module mod_x\n")
    out, n = _xform(src)
    assert n == 0, "res is not optional in the procedure that owns the region"
    assert "map(" not in out


def test_innermost_enclosing_procedure_wins():
    """With nested/contained procedures, the region belongs to the inner one."""
    src = (
        "  subroutine outer(a, opt)\n"
        "    real, intent(in), optional :: opt(nx, ny)\n"
        "  contains\n"
        "    subroutine inner(a, opt2)\n"
        "      real, intent(in), optional :: opt2(nx, ny)\n"
        "      !$omp target teams distribute parallel do collapse(2)\n"
        "      do j = 1, ny\n"
        "      do i = 1, nx\n"
        "        a(i,j) = opt2(i,j)\n"
        "      end do\n"
        "      end do\n"
        "    end subroutine inner\n"
        "  end subroutine outer\n")
    out, n = _xform(src)
    assert n == 1
    assert "map(to: opt2)" in out
    assert "opt)" not in out.split("!$omp target")[1].split("\n")[0]


# --------------------------------------------------------------------------
# Things that look like a region but are not
# --------------------------------------------------------------------------

def test_declare_target_is_not_a_region():
    """`!$omp declare target` declares a device routine; it maps nothing.

    A loop referencing the optional FOLLOWS the directive on purpose: without
    it `region_end` finds no loop and the test passes no matter what the
    region regex does. With it, a regex that swallows `declare target` appends
    a map() clause to a declarative directive, which is what we're guarding.
    """
    src = _proc(
        decls=("    real, intent(in), optional :: opt(nx, ny)\n"
               "    !$omp declare target\n"),
        body=("    do i = 1, nx\n"
              "      a(i,1) = opt(i,1)\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 0, out
    declare = [l for l in out.splitlines() if "declare target" in l][0]
    assert "map(" not in declare, declare


def test_target_data_directives_are_not_compute_regions():
    """Data-mapping directives take no map() of ours.

    As above, each is followed by a loop touching the optional so the negative
    lookaheads in the region regex are actually exercised.
    """
    for d in ("!$omp target enter data map(to: a)",
              "!$omp target exit data map(delete: a)",
              "!$omp target update from(a)"):
        src = _proc(
            decls=("    real, intent(in), optional :: opt(nx, ny)\n"),
            body=(f"    {d}\n"
                  "    do i = 1, nx\n"
                  "      a(i,1) = opt(i,1)\n"
                  "    end do\n"))
        out, n = _xform(src)
        assert n == 0, f"{d} -> {out}"


# --------------------------------------------------------------------------
# Continuations
# --------------------------------------------------------------------------

def test_continued_declaration_seen_whole():
    """`optional` and the entity list on different physical lines."""
    src = _proc(
        decls=("    real, intent(in), &\n"
               "      optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1, out
    assert "map(to: opt)" in out


def test_clause_appended_to_last_continuation_line():
    """A directive already split across `&` gets the clause on its LAST line."""
    src = _proc(
        decls=("    real, intent(in), optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2) &\n"
              "    !$omp private(tmp)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1, out
    lines = out.splitlines()
    # the directive must stay syntactically one logical line
    dir_lines = [l for l in lines if "!$omp" in l and "map(to: opt)" in l]
    assert dir_lines, out
    assert "private(tmp)" in out
    # every directive line but the last must end in the continuation marker
    di = [i for i, l in enumerate(lines) if l.strip().startswith("!$omp")]
    for i in di[:-1]:
        if lines[i + 1].strip().startswith("!$omp"):
            assert lines[i].rstrip().endswith("&"), lines[i]


def test_long_directive_wraps_rather_than_exceeding_line_limit():
    names = ", ".join(f"v{k}" for k in range(12))
    src = _proc(
        decls=(f"    real, intent(in), optional :: opt(nx, ny)\n"),
        body=(f"    !$omp target teams distribute parallel do collapse(2) private({names})\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1
    for l in out.splitlines():
        assert len(l.rstrip()) <= m.MAX_LINE, f"{len(l)} chars: {l}"


# --------------------------------------------------------------------------
# Grouping and idempotence
# --------------------------------------------------------------------------

def test_same_map_type_grouped_into_one_clause():
    src = _proc(
        decls=("    real, intent(in), optional :: p(nx, ny), q(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = p(i,j) + q(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 1, "one clause, not one per name"
    assert "map(to: p, q)" in out, out


def test_both_map_types_emitted_separately():
    src = _proc(
        decls=("    real, intent(in), optional :: p(nx, ny)\n"
               "    real, intent(inout), optional :: w(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      w(i,j) = p(i,j)\n"
              "    end do\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 2
    assert "map(to: p)" in out and "map(tofrom: w)" in out, out


def test_idempotent():
    """regen re-runs this pass over trees it has already processed."""
    src = _proc(
        decls=("    real, intent(in), optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do collapse(2)\n"
              "    do j = 1, ny\n"
              "    do i = 1, nx\n"
              "      a(i,j) = opt(i,j)\n"
              "    end do\n"
              "    end do\n"))
    once, n1 = _xform(src)
    twice, n2 = m.transform_lines(once.splitlines())
    assert n1 == 1 and n2 == 0, f"second pass added {n2}"
    assert "\n".join(twice) == once
    assert once.count("map(to: opt)") == 1


def test_two_regions_in_one_procedure_both_mapped():
    src = _proc(
        decls=("    real, intent(in), optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do\n"
              "    do i = 1, nx\n"
              "      a(i,1) = opt(i,1)\n"
              "    end do\n"
              "    !$omp target teams distribute parallel do\n"
              "    do i = 1, nx\n"
              "      a(i,2) = opt(i,2)\n"
              "    end do\n"))
    out, n = _xform(src)
    assert n == 2, out
    assert out.count("map(to: opt)") == 2


def test_region_scope_does_not_leak_past_its_loop():
    """An optional used only AFTER the region must not be mapped into it."""
    src = _proc(
        decls=("    real, intent(in), optional :: opt(nx, ny)\n"),
        body=("    !$omp target teams distribute parallel do\n"
              "    do i = 1, nx\n"
              "      a(i,1) = 0.0\n"
              "    end do\n"
              "    x = opt(1,1)\n"))
    out, n = _xform(src)
    assert n == 0, out
