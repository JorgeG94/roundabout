"""Unit tests for CPP locality-macro unwrapping (`DO_LOCALITY(local(...))`).

Run with: pytest tools/test_locality_macro.py

Regression guard for a WRONG-ANSWER bug, not a cosmetic one: the translators
read unpreprocessed source, so a `local(...)` hidden behind MOM6's
`DO_LOCALITY` macro is invisible to a naive clause parser and the emitted
OpenMP loop silently loses its `private(...)`, making the variables shared
across threads.
"""
import importlib.util
import pathlib
import sys

_HERE = pathlib.Path(__file__).parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    # Register before exec: @dataclass resolves annotations via
    # sys.modules[cls.__module__], which is None for an unregistered module.
    sys.modules[name] = mod
    sys.path.insert(0, str(_HERE))
    spec.loader.exec_module(mod)
    return mod


lm = _load("_locality_macro")
u = lm.unwrap_locality_macros


def test_simple_unwrap():
    assert u("do concurrent (i=1:n) DO_LOCALITY(local(a, b))") == \
        "do concurrent (i=1:n) local(a, b)"


def test_nested_parens_in_payload():
    # A mask argument inside the header must survive untouched.
    assert u("do concurrent (I=is:ie, do_I(I)) DO_LOCALITY(local(u_new, duhdu))") == \
        "do concurrent (I=is:ie, do_I(I)) local(u_new, duhdu)"


def test_multiple_specifiers_in_one_macro():
    assert u("do concurrent (i=1:n) DO_LOCALITY(local(x) shared(y))") == \
        "do concurrent (i=1:n) local(x) shared(y)"


def test_false_branch_semicolon_dropped():
    # `#define DO_LOCALITY(X) ;` — the header as the unsupported path sees it.
    assert u("do concurrent (i=1:n) ;").rstrip() == "do concurrent (i=1:n)"


def test_bare_specifier_untouched():
    assert u("do concurrent (i=1:n) local(z)") == "do concurrent (i=1:n) local(z)"


def test_unbalanced_macro_left_alone():
    # A continuation the caller failed to join must not be corrupted.
    src = "do concurrent (i=1:n) DO_LOCALITY(local(a,"
    assert u(src) == src


def test_custom_macro_name():
    assert u("do concurrent (i=1:n) MY_LOC(local(a))", names=("MY_LOC",)) == \
        "do concurrent (i=1:n) local(a)"


def test_statement_separator_semicolon_preserved():
    # A real `;` between two statements is not the macro's false branch.
    src = "a = 1; b = 2"
    assert u(src) == src


def test_dc_to_omp_emits_private_through_macro(tmp_path):
    dc_to_omp = _load("dc_to_omp")
    f = tmp_path / "m.F90"
    f.write_text(
        "subroutine k(n, a)\n"
        "  integer :: i, n\n"
        "  real :: a(n), tmp\n"
        "  do concurrent (i=1:n) DO_LOCALITY(local(tmp))\n"
        "    tmp = a(i)\n"
        "    a(i) = tmp\n"
        "  end do\n"
        "end subroutine k\n")
    lines, n = dc_to_omp.transform_lines(f.read_text().splitlines(), "gpu")
    out = "\n".join(lines)
    assert n == 1
    assert "private(tmp)" in out
