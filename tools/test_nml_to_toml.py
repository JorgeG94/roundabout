"""Unit tests for tools/nml_to_toml.py (standard-library unittest).

Run:  python3 -m unittest tools/test_nml_to_toml.py   (from the repo root)
"""
import importlib.util
import pathlib
import sys
import unittest

_HERE = pathlib.Path(__file__).parent
_spec = importlib.util.spec_from_file_location("nml_to_toml", _HERE / "nml_to_toml.py")
nt = importlib.util.module_from_spec(_spec)
sys.modules["nml_to_toml"] = nt
_spec.loader.exec_module(nt)

SCHEMA = nt.Schema({"groups": [
    {"name": "time", "keys": [
        {"name": "t_end", "kind": "real"}, {"name": "dt_fixed", "kind": "real"},
        {"name": "time_unit", "kind": "string"}, {"name": "n", "kind": "int"}]},
    {"name": "ocean_ic", "keys": [
        {"name": "alpha_T", "kind": "real"}, {"name": "flag", "kind": "logical"},
        {"name": "layer_rho_init", "kind": "real_array"},
        {"name": "mode", "kind": "enum"}]},
    {"name": "ocean_bt", "keys": [{"name": "auto", "kind": "logical"}]},
]})


def conv(text):
    return nt.convert(text, SCHEMA)


class Numbers(unittest.TestCase):
    def test_real_forms(self):
        for src, want in [("1.", "1.0"), (".5", "0.5"), ("-.5", "-0.5"),
                          ("1d3", "1e3"), ("1.d3", "1.0e3"), ("2.5D-3", "2.5e-3"),
                          ("1.e5", "1.0e5"), ("1.0e-2", "1.0e-2"),
                          ("3.8356948e-2", "3.8356948e-2"), ("+1.5", "1.5"),
                          ("007.5", "7.5"), ("5", "5"), ("0.0", "0.0")]:
            with self.subTest(src=src):
                self.assertEqual(nt.fmt_number(src, False), want)

    def test_int(self):
        self.assertEqual(nt.fmt_number("+7", True), "7")
        self.assertEqual(nt.fmt_number("007", True), "7")
        with self.assertRaises(nt.ConvertError):
            nt.fmt_number("1.0", True)

    def test_garbage(self):
        for bad in ("abc", "1.0_wp", "3*1.0", "."):
            with self.subTest(bad=bad), self.assertRaises(nt.ConvertError):
                nt.fmt_number(bad, False)

    def test_bool(self):
        for t in (".true.", ".TRUE.", ".t.", "T", "true"):
            self.assertEqual(nt.fmt_bool(t), "true")
        for t in (".false.", ".F.", "f"):
            self.assertEqual(nt.fmt_bool(t), "false")


class Strings(unittest.TestCase):
    def test_quotes(self):
        self.assertEqual(nt.fmt_string("'ocean'"), '"ocean"')
        self.assertEqual(nt.fmt_string("'it''s'"), '"it\'s"')
        self.assertEqual(nt.fmt_string('"a""b"'), '"a\\"b"')
        self.assertEqual(nt.fmt_string("'a\\b'"), '"a\\\\b"')


class Convert(unittest.TestCase):
    def test_basic_layout_and_comments(self):
        got = conv("! top\n\n&time_nml\n   t_end = 10.   ! ten\n"
                   "   ! note\n   time_unit = 'day'\n/\n")
        self.assertEqual(got, "# top\n\n[time]\nt_end = 10.0  # ten\n"
                              "# note\ntime_unit = \"day\"\n")

    def test_bang_inside_string_is_kept(self):
        got = conv("&time_nml\n time_unit = 'a!b' ! c\n/\n")
        self.assertIn('time_unit = "a!b"  # c', got)

    def test_one_line_groups_and_multi_key(self):
        got = conv("&time_nml t_end = 1.0, dt_fixed = 2d0 /\n"
                   "&ocean_bt_nml auto = .true. /\n")
        self.assertEqual(got, "[time]\nt_end = 1.0\ndt_fixed = 2e0\n"
                              "[ocean.bt]\nauto = true\n")

    def test_mapping_and_case(self):
        got = conv("&OCEAN_IC_NML\n ALPHA_t = 0.2\n LAYER_RHO_INIT = 1036., 1035.\n/\n")
        self.assertEqual(got, "[ocean.ic]\nalpha_T = 0.2\n"
                              "layer_rho_init = [1036.0, 1035.0]\n")

    def test_single_value_array_stays_array(self):
        self.assertIn("layer_rho_init = [5.0]",
                      conv("&ocean_ic_nml layer_rho_init = 5. /\n"))

    def test_placeholder_passthrough(self):
        got = conv("&time_nml\n   t_end = {{T}}\n{{EXTRA}}\n/\n")
        self.assertEqual(got, "[time]\nt_end = {{T}}\n{{EXTRA}}\n")

    def test_errors(self):
        for bad in ("&nope_nml\n/\n", "&time_nml\n nope = 1\n/\n",
                    "&time_nml\n t_end = 1\n t_end = 2\n/\n",
                    "&time_nml\n t_end = 3*1.0\n/\n",
                    "&time_nml\n/\n&time_nml\n/\n", "&time_nml\n t_end = 1\n"):
            with self.subTest(bad=bad), self.assertRaises(nt.ConvertError):
                conv(bad)

    def test_table_mapping(self):
        self.assertEqual(nt.group_to_table("ocean_hvisc"), "ocean.hvisc")
        self.assertEqual(nt.group_to_table("grid"), "grid")
        self.assertEqual(nt.group_to_table("ocean_"), "ocean_")


if __name__ == "__main__":
    unittest.main()
