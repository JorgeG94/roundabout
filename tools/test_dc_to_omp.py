"""Structural regression guards for dc_to_omp.py's block finder.

Every test here corresponds to a bug found by running the translator over
MOM6, whose `;`-compound statement style ("endif ; enddo",
"if (m==1) then ; do concurrent (...) ; press(i,j) = 0.0 ; enddo") breaks
line-anchored parsing in ways that are silent rather than loud.
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


d = _load("dc_to_omp")


def _xform(src):
    lines, n = d.transform_lines(src.splitlines(), "gpu")
    return "\n".join(lines), n


def test_unconvertible_loop_does_not_abandon_rest_of_file():
    """A skipped loop must not end the scan — None means 'no more loops'."""
    src = (
        "do concurrent (i)\n"          # no range: unconvertible
        "  a(i) = 1\n"
        "end do\n"
        "do concurrent (j=1:m)\n"
        "  b(j) = 2\n"
        "end do\n")
    out, n = _xform(src)
    assert n == 1, "the loop after the unconvertible one must still convert"
    assert "!$omp target teams distribute parallel do" in out


def test_mask_lowered_to_guard():
    """`do concurrent (i=lo:hi, MASK)` -> loop + `if (MASK) then` inside."""
    src = (
        "do concurrent (i=1:n, mask(i))\n"
        "  a(i) = 1\n"
        "end do\n")
    out, n = _xform(src)
    assert n == 1
    body = out.splitlines()
    assert body[0].strip() == "!$omp target teams distribute parallel do"
    assert body[1].strip() == "do i = 1, n"
    assert body[2].strip() == "if ((mask(i))) then"
    assert body[4].strip() == "end if"
    assert body[5].strip() == "end do"


def test_multiple_masks_conjoined():
    src = ("do concurrent (I=is:ie, do_I(I), flag)\n"
           "  a(I) = 1\n"
           "end do\n")
    out, n = _xform(src)
    assert n == 1
    assert "if ((do_I(I)) .and. (flag)) then" in out


def test_mask_keeps_locality_private():
    src = ("do concurrent (I=is:ie, do_I(I)) DO_LOCALITY(local(u_new))\n"
           "  u_new = 1\n"
           "  x(I) = u_new\n"
           "end do\n")
    out, n = _xform(src)
    assert n == 1
    assert "private(u_new)" in out and "if ((do_I(I))) then" in out


def test_string_literal_is_not_a_loop():
    """A character literal naming the construct must not derail the scan."""
    src = (
        'call check(error, ok, "do concurrent reduce(+) must be exact")\n'
        "do concurrent (i=1:n)\n"
        "  a(i) = 1\n"
        "end do\n")
    out, n = _xform(src)
    assert n == 1
    assert 'call check(error, ok, "do concurrent reduce(+) must be exact")' in out


def test_compound_enddo_closes_depth():
    """`endif ; enddo` closes an inner do; the outer end do is the match."""
    lines = [
        "do concurrent (i=1:n)",          # 0  header
        "  do k=1,nk ; if (c(k)) then",   # 1  opens a do
        "    x = 1",                      # 2
        "  endif ; enddo",                # 3  closes it — NOT at line start
        "end do",                         # 4  the true match
    ]
    assert d.find_matching_end_do(lines, 0) == 4


def test_compound_do_is_counted():
    """A `do` after `;` must open depth, or an outer end do is stolen."""
    lines = [
        "do concurrent (i=1:n)",              # 0  header
        "  if (m==1) then ; do k=1,nk",       # 1  opens a do, not at line start
        "    x = 1",                          # 2
        "  enddo ; endif",                    # 3  closes it
        "end do",                             # 4  the true match
    ]
    # The old line-anchored matcher missed the `do` on line 1 and returned 3.
    assert d.find_matching_end_do(lines, 0) == 4


def test_suffix_compound_header_converted_keeping_trailing_statement():
    """MOM6's `do concurrent (...) ; if (...) then` — the `if` is body, not clause."""
    src = (
        "do concurrent (j=js:je, I=is-1:ie) ; if (mask(I,j) > 0.0) then\n"
        "  a(I,j) = 1\n"
        "endif ; enddo\n")
    out, n = _xform(src)
    assert n == 1
    assert "if (mask(I,j) > 0.0) then" in out, "trailing statement must survive"
    assert "endif" in out
    assert out.count("end do") == 2 and "collapse(2)" in out


def test_prefix_compound_header_converted():
    """`do k=1,nz ; do concurrent (i=is:ie)` — the outer do must survive."""
    src = ("do k=1,nz ; do concurrent (i=is:ie)\n"
           "  b(i) = 2\n"
           "enddo ; enddo\n")
    out, n = _xform(src)
    lines = [l.strip() for l in out.splitlines()]
    assert n == 1
    assert lines[0] == "do k=1,nz"
    assert lines[1] == "!$omp target teams distribute parallel do"
    assert lines[-1] == "enddo", "the outer enddo must survive"


def test_prefix_if_compound_header_converted():
    src = ("if (use_temperature) then ; do concurrent (j=js:je, i=is:ie)\n"
           "  c(i,j) = 3\n"
           "enddo ; endif\n")
    out, n = _xform(src)
    lines = [l.strip() for l in out.splitlines()]
    assert n == 1
    assert lines[0] == "if (use_temperature) then"
    assert lines[-1] == "endif"


def test_one_line_loop_converted():
    """`do concurrent (i=is:ie) ; press(i,j) = 0.0 ; enddo` — body must survive."""
    src = "do concurrent (i=is:ie) ; press(i,j) = 0.0 ; enddo\n"
    out, n = _xform(src)
    lines = [l.strip() for l in out.splitlines()]
    assert n == 1
    assert lines == ["!$omp target teams distribute parallel do",
                     "do i = is, ie",
                     "press(i,j) = 0.0",
                     "end do",
                     "!$omp end target teams distribute parallel do"]


def test_normalisation_leaves_unrelated_compound_lines_alone():
    """Only header / end-do lines are split; the rest keeps MOM6's style."""
    src = ("a = 1 ; b = 2\n"
           "do concurrent (i=1:n)\n"
           "  x = 1 ; y = 2\n"
           "end do\n")
    out, _ = _xform(src)
    assert "a = 1 ; b = 2" in out
    assert "  x = 1 ; y = 2" in out


def test_normalisation_skips_continuations_and_comments():
    src = ("! a comment ; with a semicolon\n"
           "call f(a, &\n"
           "       b) ; call g()\n"
           "do concurrent (i=1:n)\n"
           "  z = 1\n"
           "end do\n")
    out, _ = _xform(src)
    assert "! a comment ; with a semicolon" in out
    assert "       b) ; call g()" in out


def test_shared_closing_line_refused():
    """A closing `end do` sharing its line must not be swallowed by the footer."""
    lines = [
        "do concurrent (i=1:n)",
        "  x = 1",
        "end do ; y = 2",
    ]
    assert d.find_matching_end_do(lines, 0) is None


def test_skip_message_reports_original_line():
    """Line numbers must survive the drift from preceding rewrites."""
    import io
    import contextlib
    src = (
        "do concurrent (i=1:n)\n"      # 1  converts: 1 line -> 3
        "  a(i) = 1\n"                 # 2
        "end do\n"                     # 3
        "do concurrent (j)\n"          # 4  unconvertible — must report line 4
        "  b(j) = 2\n"
        "end do\n")
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        _xform(src)
    assert "line 4" in err.getvalue(), err.getvalue()


def test_locality_macro_survives_to_private():
    src = (
        "do concurrent (j=1:ny, i=1:nx) DO_LOCALITY(local(tmp))\n"
        "  tmp = a(i,j)\n"
        "  b(i,j) = tmp\n"
        "end do\n")
    out, n = _xform(src)
    assert n == 1
    assert "private(tmp)" in out


def test_nested_loops_get_one_directive():
    """An inner `do concurrent` must NOT get its own directive.

    OpenMP forbids a `target` region inside a `target` region, and the inner
    loop is already covered by the outer construct's parallelism.
    """
    src = ("do concurrent (k=1:nz)\n"
           "  do concurrent (j=js:je, i=is:ie)\n"
           "    x(i,j,k) = 1.0\n"
           "  enddo\n"
           "enddo\n")
    out, n = _xform(src)
    assert n == 2, "both loops are rewritten"
    assert out.count("!$omp target teams distribute parallel do") == 1
    assert out.count("!$omp end target teams distribute parallel do") == 1
    assert "do j = js, je" in out and "do i = is, ie" in out


def test_nested_locality_hoisted_to_outer_private():
    """An inner loop's `local(...)` becomes the OUTER directive's `private(...)`."""
    src = ("do concurrent (k=1:nz)\n"
           "  do concurrent (i=is:ie) DO_LOCALITY(local(tmp))\n"
           "    tmp = 0.0\n"
           "    a(i,k) = tmp\n"
           "  enddo\n"
           "enddo\n")
    out, _ = _xform(src)
    directive = [l for l in out.splitlines() if l.strip().startswith("!$omp target")][0]
    assert "private(tmp)" in directive


def test_privates_deduplicated():
    src = ("do concurrent (k=1:nz) DO_LOCALITY(local(tmp))\n"
           "  do concurrent (i=is:ie) DO_LOCALITY(local(tmp))\n"
           "    tmp = 0.0\n"
           "    a(i,k) = tmp\n"
           "  enddo\n"
           "enddo\n")
    out, _ = _xform(src)
    directive = [l for l in out.splitlines() if l.strip().startswith("!$omp target")][0]
    assert directive.count("tmp") == 1, directive


def test_nested_mask_lowered_without_directive():
    src = ("do concurrent (k=1:nz)\n"
           "  do concurrent (i=is:ie, do_I(i))\n"
           "    a(i,k) = 1.0\n"
           "  enddo\n"
           "enddo\n")
    out, _ = _xform(src)
    assert out.count("!$omp target teams distribute parallel do") == 1
    assert "if ((do_I(i))) then" in out


def test_prefix_statement_on_continued_header_line():
    """`do K=2,nz ; do concurrent(I=is-1:ie) &` — the outer `do` must survive.

    Regression: the normaliser used to skip any line ending in `&`, so this
    prefix was never split off and the emitter replaced the whole line,
    deleting `do K=2,nz` while its `enddo` remained.
    """
    src = ("do K=2,nz ; do concurrent(I=is-1:ie) &\n"
           "    DO_LOCALITY(local(Hdn))\n"
           "  Hdn = 1.0\n"
           "  a(I,K) = Hdn\n"
           "enddo ; enddo\n")
    out, n = _xform(src)
    lines = [l.strip() for l in out.splitlines()]
    assert n == 1
    assert lines[0] == "do K=2,nz", lines
    assert "private(Hdn)" in out
    assert lines.count("enddo") + lines.count("end do") == 2


def test_leading_code_never_silently_deleted():
    """Whatever the normaliser declines to split must be refused, not dropped."""
    import io
    import contextlib
    # A continuation line the normaliser will not cut.
    src = ("call setup(a, &\n"
           "     b) ; do concurrent (i=1:n)\n"
           "  x(i) = 1\n"
           "end do\n")
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        out, n = _xform(src)
    assert "call setup(a, &" in out and "b) ;" in out, "leading code preserved"


def test_loop_inside_handwritten_target_region_gets_no_directive():
    """OpenMP forbids a target region inside a target region.

    MOM6's GPU branch wraps `do j=...` in a hand-written `!$omp target teams
    loop` and puts `do concurrent (i=...)` inside it. That inner loop must
    become a plain nest.
    """
    src = ("!$omp target teams loop private(depth)\n"
           "do j=js,je\n"
           "  do concurrent (i=is:ie)\n"
           "    depth(i) = 0.0\n"
           "  enddo\n"
           "enddo\n")
    out, n = _xform(src)
    assert n == 1
    assert out.count("!$omp target teams loop") == 1
    assert "distribute parallel do" not in out
    assert "do i = is, ie" in out


def test_loop_after_handwritten_region_still_gets_a_directive():
    """The suppression must end with the region, not swallow the rest of the file."""
    src = ("!$omp target teams loop\n"
           "do j=js,je\n"
           "  do concurrent (i=is:ie)\n"
           "    a(i,j) = 0.0\n"
           "  enddo\n"
           "enddo\n"
           "do concurrent (i=is:ie)\n"
           "  b(i) = 1.0\n"
           "enddo\n")
    out, n = _xform(src)
    assert n == 2
    assert out.count("!$omp target teams distribute parallel do") == 1, out
