import contextlib
import io
import unittest
from unittest.mock import mock_open, patch

from run_simulation import main


class FragmentationConfigTests(unittest.TestCase):
    def run_wrapper(self, yaml_text):
        with patch("builtins.open", mock_open(read_data=yaml_text)), \
             patch("sys.argv", ["run_simulation.py", "-c", "test.yaml"]), \
             patch("run_simulation.run_c_program", return_value=(True, 0)) as run, \
             contextlib.redirect_stdout(io.StringIO()):
            main()
            return run.call_args

    def test_explicit_factor_mapping(self):
        for section in ("simulation_parameters", "dust_parameters"):
            with self.subTest(section=section):
                call = self.run_wrapper(
                    f"{section}:\n  fragmentation_factor: 0.52\n"
                )
                self.assertEqual(call.args[1]["fragmentation_factor"], 0.52)
                self.assertEqual(call.args[2]["fragmentation_factor"], "-ffrag")

    def test_omitted_factor_preserves_binary_default(self):
        call = self.run_wrapper(
            "simulation_parameters:\n  enable_dust_growth: true\n"
        )
        self.assertNotIn("fragmentation_factor", call.args[1])


if __name__ == "__main__":
    unittest.main()
