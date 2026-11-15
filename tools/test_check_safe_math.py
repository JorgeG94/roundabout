"""Unit tests for tools/check_safe_math.py (standard-library unittest; pytest runs them too).

Run:  python3 -m unittest tools/test_check_safe_math.py   (from the repo root)
"""
import importlib.util
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

_HERE = pathlib.Path(__file__).parent
_spec = importlib.util.spec_from_file_location("check_safe_math", _HERE / "check_safe_math.py")
sm = importlib.util.module_from_spec(_spec)
sys.modules["check_safe_math"] = sm
_spec.loader.exec_module(sm)


def flagged(src: str):
    """[(line, kind, what)] for a Fortran snippet."""
    return [(h.line, h.kind, h.what) for h in sm.scan_text(src)]


class TruePositives(unittest.TestCase):
    def test_every_intrinsic_call(self):
        for name in sm.INTRINSICS:
            with self.subTest(name=name):
                self.assertEqual(flagged(f"y = {name}(x)\n"), [(1, "call", name)])

    def test_case_and_spacing(self):
        self.assertEqual(flagged("y = EXP (x)\n"), [(1, "call", "exp")])

    def test_two_on_a_line(self):
        self.assertEqual(len(flagged("y = exp(a) + cos(b)\n")), 2)

    def test_after_relational_dot_operator(self):
        self.assertEqual(len(flagged("if (a .lt. log(b)) y = 1\n")), 1)

    def test_real_exponents(self):
        for expr in ("x**0.5_wp", "x**(1.0_wp/3.0_wp)", "x**(-0.5)", "x**1.1", "x**2e0",
                     "x**(2/3)", "x**(1.0_wp/real(n,wp))", "x**(kappa*c)", "x**real(n,wp)",
                     "x**-0.25_wp", "x **0.5_wp"):
            with self.subTest(expr=expr):
                self.assertEqual([k for _, k, _ in flagged(f"y = {expr}\n")], ["pow"])

    def test_continuation_flags_the_call_line(self):
        src = "y = a + &\n    exp(b)\n"
        self.assertEqual(flagged(src), [(2, "call", "exp")])


class FalsePositiveClasses(unittest.TestCase):
    def clean(self, src):
        self.assertEqual(flagged(src), [], src)

    def test_integer_powers(self):
        for expr in ("x**2", "x**3", "x**n", "x**(n-1)", "x**(2)", "x**(-2)", "x**2_i4", "x**ip",
                     "10**nk", "x**(i+1)", "x**-2"):
            with self.subTest(expr=expr):
                self.clean(f"y = {expr}\n")

    def test_comments(self):
        self.clean("! y = exp(x) and x**0.5\n")
        self.clean("   !! uses cos(theta)\n")
        self.clean("y = z  ! was exp(x)\n")

    def test_strings(self):
        self.clean("write_msg = 'exp(x) and 2**0.5'\n")
        self.clean('msg = "cos(theta) ! not a comment start"\n')

    def test_string_with_bang_then_real_call(self):
        # the `!` is inside a string, so the call after it is code
        self.assertEqual(len(flagged("msg = 'a!b'; y = exp(x)\n")), 1)

    def test_names_containing_words(self):
        for line in ("y = cosine_x(a)", "y = explicit(a)", "y = safe_exp(a)", "y = rdb_sm_exp(a)",
                     "y = my_log(a)", "y = exp_decay(a)", "y = tanh2(a)", "call log_message(a)",
                     "y = logistic(a)", "y = alog(a)", "y = dexp (a)"):
            with self.subTest(line=line):
                self.clean(line + "\n")

    def test_type_bound_and_component(self):
        self.clean("y = obj%exp(a)\n")
        self.clean("call logger%log(msg)\n")

    def test_shadowing_array_or_dummy(self):
        self.clean("real(wp) :: exp(n)\ny = exp(i) * 2\n")
        self.clean("real(wp), intent(in) :: a, log(3)\ny = log(1)\n")
        self.clean("real(wp) function sin(x)\nend function sin\ny = sin(1.0)\n")
        self.clean("interface cos\nend interface\ny = cos(1.0)\n")

    def test_declaration_initialiser_is_not_a_shadow(self):
        # `x = atan(1.0)` on a declaration line must not hide later atan calls
        src = "real(wp) :: x = atan(1.0_wp)\ny = atan(z)\n"
        self.assertIn((2, "call", "atan"), flagged(src))

    def test_parameter_initialisers(self):
        self.clean("real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)\n")
        self.clean("real(wp), parameter, public :: TWO_PI = 8*atan(1.0_wp), E = exp(1.0_wp)\n")
        self.clean("real(wp), parameter :: R = 2.0_wp**0.5_wp\n")

    def test_preprocessor_lines(self):
        self.clean("#define LOGC log(1.0)\n")

    def test_exempt_library_paths(self):
        self.assertTrue(sm.exempt("src/framework/rdb_safe_math.F90"))
        self.assertTrue(sm.exempt("src/shared_module_utilities/rdb_safe_math_trig.inc"))
        self.assertFalse(sm.exempt("src/shared_module_utilities/rdb_frhat_face.inc"))
        self.assertFalse(sm.exempt("src/core/ocean/x/rdb_safe_math.F90"))


class Waiver(unittest.TestCase):
    def test_reason_waives(self):
        self.assertEqual(flagged("y = exp(x)  ! safe-math-ok: setup only, one rank value\n"), [])

    def test_reason_required(self):
        self.assertEqual(flagged("y = exp(x)  ! safe-math-ok:\n"), [(1, "waiver", "safe-math-ok needs a reason")])

    def test_waiver_is_per_line(self):
        src = "! safe-math-ok: nope, wrong line\ny = exp(x)\n"
        self.assertEqual(flagged(src), [(2, "call", "exp")])

    def test_waiver_covers_pow(self):
        self.assertEqual(flagged("y = x**0.5_wp ! safe-math-ok: sqrt-exact\n"), [])


class Context(unittest.TestCase):
    def test_procedure_tracking(self):
        src = (
            "module m\ncontains\n"
            "  subroutine alpha(x)\n    y = exp(x)\n  end subroutine alpha\n"
            "  pure real(wp) function beta(x)\n    beta = cos(x)\n  end function beta\n"
            "end module m\n"
        )
        self.assertEqual([h.proc for h in sm.scan_text(src)], ["alpha", "beta"])

    def test_inventory_counts(self):
        hits = sm.scan_text("subroutine s()\n y = exp(a) + exp(b) + x**0.5_wp\nend subroutine s\n", "src/core/ocean/a.F90")
        out = sm.inventory(hits)
        self.assertIn("3 bare sites in 1 files", out)
        self.assertIn("2  exp", out)
        self.assertIn("src/core/ocean/a.F90 :: s", out)


DIFF = """\
diff --git a/src/a.F90 b/src/a.F90
--- a/src/a.F90
+++ b/src/a.F90
@@ -3 +3 @@ x
-old
+new
@@ -10,0 +11,2 @@
+a
+b
diff --git a/src/gone.F90 b/src/gone.F90
--- a/src/gone.F90
+++ /dev/null
@@ -1,2 +0,0 @@
-x
-y
"""


class DiffMode(unittest.TestCase):
    def test_parse_unified0(self):
        got = sm.parse_unified0(DIFF)
        self.assertEqual(got["src/a.F90"], {3, 11, 12})
        self.assertNotIn("src/gone.F90", got)

    def test_git_diff_only_new_lines_are_errors(self):
        # throw-away repo: a file with one old bare exp; add a second; hook-mode flags only the new one
        with tempfile.TemporaryDirectory() as d:
            repo = pathlib.Path(d)
            # pre-commit exports GIT_DIR / GIT_INDEX_FILE: they would redirect git away from the throw-away repo
            env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
            g = lambda *a: subprocess.run(["git", "-C", d, *a], check=True, capture_output=True, text=True, env=env)
            g("init", "-q", "-b", "main")
            g("config", "user.email", "t@t")
            g("config", "user.name", "t")
            (repo / "src").mkdir()
            f = repo / "src" / "k.F90"
            f.write_text("subroutine s()\n  y = exp(a)\nend subroutine s\n")
            g("add", "-A")
            g("commit", "-q", "-m", "base")
            f.write_text("subroutine s()\n  y = exp(a)\n  z = cos(b)\n  w = sin(c) ! safe-math-ok: test\nend subroutine s\n")
            g("add", "-A")
            old_root = sm.ROOT
            sm.ROOT = repo
            try:
                with mock.patch.dict(os.environ, env, clear=True):
                    changed = sm.changed_lines("main")
                hits = sm.scan_text(f.read_text(), "src/k.F90")
                new = [(h.line, h.what) for h in hits if h.line in changed.get("src/k.F90", ())]
                self.assertEqual(new, [(3, "cos")])  # exp on line 2 is pre-existing; sin is waived
                self.assertEqual(len([h for h in hits]), 2)  # the full scan still sees both unwaived sites
            finally:
                sm.ROOT = old_root


if __name__ == "__main__":
    unittest.main()
