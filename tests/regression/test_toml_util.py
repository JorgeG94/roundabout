"""Unit tests for tests/regression/toml_util.py (standard-library unittest).

Run:  python3 -m unittest tests/regression/test_toml_util.py   (from the repo root)

Covers the three patch operations the drivers rely on (replace a key in a
table, add a key to an existing table, add a new table) plus rendering and
the group -> table mapping.  No `tomllib`: the module must run under 3.10.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toml_util as tu  # noqa: E402

SRC = """\
# header comment
[time]
t_end = 10.0   # trailing
time_unit = "day"
dt_fixed = 1200.0

[ocean.bt]
auto_n_inner = true

[output]
output_dir = "./out#1"
"""


class Mapping(unittest.TestCase):
    def test_group_to_table(self):
        self.assertEqual(tu.group_to_table("time_nml"), "time")
        self.assertEqual(tu.group_to_table("ocean_bt_nml"), "ocean.bt")
        self.assertEqual(tu.group_to_table("ocean_hvisc"), "ocean.hvisc")
        # only the LEADING ocean_ nests
        self.assertEqual(tu.group_to_table("ocean_grid_nml"), "ocean.grid")
        self.assertEqual(tu.group_to_table("grid_nml"), "grid")


class Render(unittest.TestCase):
    def test_values(self):
        self.assertEqual(tu.fmt_value(True), "true")
        self.assertEqual(tu.fmt_value(False), "false")
        self.assertEqual(tu.fmt_value(3), "3")
        self.assertEqual(tu.fmt_value(1200.0), "1200.0")
        self.assertEqual(tu.fmt_value(1e-5), "1e-05")
        self.assertEqual(tu.fmt_value('a"b\\c'), '"a\\"b\\\\c"')
        self.assertEqual(tu.fmt_value([1.0, 2.5]), "[1.0, 2.5]")
        with self.assertRaises(TypeError):
            tu.fmt_value(None)

    def test_render(self):
        txt = tu.render({"ocean_bt_nml": {"b": 1, "a": True},
                         "time_nml": {"t": 2.0}}, header="hello\n\nworld")
        self.assertEqual(
            txt, "# hello\n#\n# world\n[ocean.bt]\na = true\nb = 1\n\n"
                 "[time]\nt = 2.0\n")


class Patch(unittest.TestCase):
    def test_get(self):
        self.assertEqual(tu.get(SRC, "time_nml", "t_end"), "10.0")
        self.assertEqual(tu.get_float(SRC, "time_nml", "dt_fixed"), 1200.0)
        self.assertEqual(tu.get(SRC, "output_nml", "output_dir"), '"./out#1"')
        self.assertIsNone(tu.get(SRC, "time_nml", "nope"))
        self.assertIsNone(tu.get(SRC, "ocean_diag_nml", "enabled"))
        self.assertTrue(tu.has_group(SRC, "ocean_bt_nml"))
        self.assertFalse(tu.has_group(SRC, "ocean_diag_nml"))

    def test_replace_key(self):
        out = tu.patch(SRC, {"time_nml": {"t_end": 5.0, "time_unit": "second"}})
        self.assertIn("t_end = 5.0\n", out)
        self.assertIn('time_unit = "second"\n', out)
        self.assertEqual(out.count("t_end"), 1)
        self.assertIn("# header comment", out)
        self.assertEqual(tu.get_float(out, "time_nml", "t_end"), 5.0)

    def test_add_key_to_existing_table(self):
        out = tu.patch(SRC, {"ocean_bt_nml": {"bt_halo": 0}})
        lines = out.splitlines()
        i = lines.index("[ocean.bt]")
        self.assertEqual(lines[i + 1:i + 3], ["auto_n_inner = true", "bt_halo = 0"])
        # the next table is untouched and still the next table
        self.assertEqual(tu.get(out, "output_nml", "output_dir"), '"./out#1"')

    def test_add_new_table(self):
        out = tu.patch(SRC, {"mpi_nml": {"px": 2, "py": 1}})
        self.assertTrue(out.endswith('output_dir = "./out#1"\n\n[mpi]\npx = 2\npy = 1\n'))
        out2 = tu.patch(SRC, {"ocean_diag_nml": {"enabled": False}})
        self.assertIn("\n[ocean.diag]\nenabled = false\n", out2)

    def test_table_from_empty(self):
        self.assertEqual(tu.patch("", {"time_nml": {"a": 1}}), "[time]\na = 1\n")

    def test_idempotent(self):
        ov = {"time_nml": {"t_end": 7.0}, "ocean_bt_nml": {"bt_halo": 0},
              "mpi_nml": {"px": 1}}
        once = tu.patch(SRC, ov)
        self.assertEqual(tu.patch(once, ov), once)

    def test_multiline_array_refused(self):
        with self.assertRaises(ValueError):
            tu.patch("[a]\nx = [1,\n 2]\n", {"a": {"y": 1}})


if __name__ == "__main__":
    unittest.main()
