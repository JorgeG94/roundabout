"""Unit tests for tools/check_no_vector_math.py (standard-library unittest).

Run:  python3 -m unittest tools/test_check_no_vector_math.py   (from the repo root)
"""
import importlib.util
import pathlib
import sys
import unittest

_HERE = pathlib.Path(__file__).parent
_spec = importlib.util.spec_from_file_location("check_no_vector_math", _HERE / "check_no_vector_math.py")
nv = importlib.util.module_from_spec(_spec)
sys.modules["check_no_vector_math"] = nv
_spec.loader.exec_module(nv)


def matches(sym):
    return any(p.match(sym) for p in nv.PATTERNS)


class Patterns(unittest.TestCase):
    def test_glibc_and_nvhpc_vector_entries(self):
        for s in ("_ZGVdN4v_exp", "_ZGVbN2v_cos", "_ZGVbN2vv_pow", "_ZGVeN8v_sin",
                  "__fd_exp_4_avx2", "__fd_sin_2_avx2", "__fd_exp_4_mn", "__fs_sin_8_avx2",
                  "__pd_exp_2", "__ps_cos_4", "__gvd_exp_4", "__gvs_sin_8", "__mth_i_dpowd"):
            with self.subTest(sym=s):
                self.assertTrue(matches(s))

    def test_ordinary_symbols_do_not_match(self):
        for s in ("exp", "cos", "sqrt", "__gfortran_st_write", "_gfortran_pow_r8_i4", "memcpy",
                  "__rdb_safe_math_MOD_safe_exp", "__fdsomething", "_ZN4pic"):
            with self.subTest(sym=s):
                self.assertFalse(matches(s))


class Compare(unittest.TestCase):
    def test_new_stale_missing(self):
        found = {"a.o": ["_ZGVdN4v_exp"], "c.o": ["_ZGVdN4v_cos"]}
        allowed = {"a.o", "b.o", "d.o"}
        present = {"a.o", "b.o", "c.o"}
        new, stale, missing = nv.compare(found, allowed, present)
        self.assertEqual(new, ["c.o"])      # offender not listed
        self.assertEqual(stale, ["b.o"])    # listed, built, clean
        self.assertEqual(missing, ["d.o"])  # listed, not built: warning only

    def test_clean_when_equal(self):
        self.assertEqual(nv.compare({"a.o": ["x"]}, {"a.o"}, {"a.o"}), ([], [], []))

    def test_allowlist_parse_and_roundtrip(self):
        found = {"src/x.F90.o": ["_ZGVbN2v_exp", "_ZGVdN4v_exp"]}
        text = nv.format_allowlist(found)
        self.assertEqual(nv.parse_allowlist(text), {"src/x.F90.o"})

    def test_object_key_strips_target_dir(self):
        b = pathlib.Path("/b")
        self.assertEqual(nv.object_key(b / "CMakeFiles/core_rdb_objs.dir/src/a/x.F90.o", b), "src/a/x.F90.o")
        self.assertEqual(nv.object_key(b / "app/CMakeFiles/rdb.dir/app/main.F90.o", b), "app/main.F90.o")


if __name__ == "__main__":
    unittest.main()
