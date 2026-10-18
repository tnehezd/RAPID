import contextlib
import io
import unittest
from unittest.mock import mock_open, patch

from run_simulation import main, validate_thermal_config


class ThermalConfigTests(unittest.TestCase):
    def test_default_and_explicit_lociso(self):
        for config in ({}, {"thermodynamics_parameters": {}},
                       {"thermodynamics_parameters": {"model": "lociso"}},
                       {"thermodynamics_parameters": {"model": "locally_isotherm"}}):
            with self.subTest(config=config):
                self.assertEqual(validate_thermal_config(config), "lociso")

    def test_invalid_settings(self):
        configs = (
            [],
            {"thermal": {"model": "radiative"}},
            {"thermodynamics_parameters": None},
            {"thermodynamics_parameters": "radiative"},
            {"thermodynamics_parameters": []},
            {"thermodynamics_parameters": {"model": None}},
            {"thermodynamics_parameters": {"model": True}},
            {"thermodynamics_parameters": {"model": "radative"}},
            {"thermodynamics_parameters": {"modle": "radiative"}},
            {"thermodynamics_parameters": {"model": "lociso", "opacity": "lin"}},
        )
        for config in configs:
            with self.subTest(config=config), self.assertRaises(ValueError):
                validate_thermal_config(config)

    def test_polytropic_is_explicit_error(self):
        with self.assertRaisesRegex(ValueError, "not implemented"):
            validate_thermal_config({"thermodynamics_parameters": {"model": "polytropic"}})

    def test_radiative_parameters(self):
        self.assertEqual(validate_thermal_config({
            "thermodynamics_parameters": {"model": "radiative"}
        }), "radiative")
        for value in (0, -1, True, None, "10", float("nan"), float("inf")):
            for key in ("mean_molecular_weight", "background_temperature_k",
                        "minimum_temperature_k", "maximum_temperature_k"):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    validate_thermal_config({
                        "thermodynamics_parameters": {"model": "radiative", key: value}
                    })
        for settings in (
            {"minimum_temperature_k": 1500},
            {"maximum_temperature_k": 5},
            {"maximum_temperature_k": 1e9},
        ):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                validate_thermal_config({
                    "thermodynamics_parameters": {"model": "radiative", **settings}
                })

    def test_energy_mode_and_gamma(self):
        call = self.run_wrapper(
            "disk_parameters:\n  star_mass_msun: 1\n"
            "thermodynamics_parameters:\n  model: radiative_energy\n  gamma: 1.6\n"
        )
        self.assertEqual(call.args[1]["thermal_model"], "radiative_energy")
        self.assertEqual(call.args[1]["gamma"], 1.6)
        self.assertEqual(call.args[2]["gamma"], "-thermal_gamma")
        for value in (0, 1, -1, True, "1.4", float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_thermal_config({
                    "thermodynamics_parameters": {"model": "radiative_energy", "gamma": value}
                })

    def run_wrapper(self, yaml_text):
        with patch("builtins.open", mock_open(read_data=yaml_text)), \
             patch("sys.argv", ["run_simulation.py", "-c", "test.yaml"]), \
             patch("run_simulation.run_c_program", return_value=(True, 0)) as run, \
             contextlib.redirect_stdout(io.StringIO()):
            main()
            return run.call_args

    def test_irradiation_settings(self):
        for settings in (
            {"irradiation_luminosity_lsun": 0, "irradiation_factor": 0},
            {"irradiation_luminosity_lsun": 1, "irradiation_factor": 0.05},
            {"irradiation_factor": 1},
        ):
            self.assertEqual(validate_thermal_config({
                "thermodynamics_parameters": {"model": "radiative", **settings}
            }), "radiative")
        for key in ("irradiation_luminosity_lsun", "irradiation_factor"):
            for value in (-1, True, None, "1", float("nan"), float("inf")):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    validate_thermal_config({
                        "thermodynamics_parameters": {"model": "radiative", key: value}
                    })
        with self.assertRaises(ValueError):
            validate_thermal_config({"thermodynamics_parameters": {"irradiation_factor": 1.1}})

    def test_lociso_preserves_binary_arguments(self):
        config = "disk_parameters:\n  star_mass_msun: 1.0\n"
        original = self.run_wrapper(config)
        self.assertIsNotNone(original)
        for model in ("lociso", "locally_isotherm"):
            with self.subTest(model=model):
                explicit = self.run_wrapper(
                    config + f"thermodynamics_parameters:\n  model: {model}\n"
                )
                self.assertEqual(original, explicit)
                self.assertEqual(explicit.args[1], {"star_val": 1.0})
                with_radiative_settings = self.run_wrapper(
                    config + f"thermodynamics_parameters:\n  model: {model}\n"
                    "  mean_molecular_weight: 2.34\n  background_temperature_k: 10\n"
                    "  irradiation_luminosity_lsun: 1\n  irradiation_factor: 0.05\n"
                )
                self.assertEqual(original, with_radiative_settings)

    def test_invalid_mode_never_launches_binary(self):
        for model in ("polytropic", "typo"):
            with self.subTest(model=model), \
                 patch("builtins.open", mock_open(
                     read_data=f"disk_parameters:\n  star_mass_msun: 1\nthermodynamics_parameters:\n  model: {model}\n"
                 )), \
                 patch("sys.argv", ["run_simulation.py"]), \
                 patch("run_simulation.run_c_program") as run, \
                 contextlib.redirect_stderr(io.StringIO()) as errors:
                with self.assertRaises(SystemExit) as exit_context:
                    main()
                self.assertEqual(exit_context.exception.code, 2)
                self.assertIn("thermodynamics_parameters.model", errors.getvalue())
                run.assert_not_called()

    def test_radiative_arguments_reach_binary(self):
        call = self.run_wrapper(
            "disk_parameters:\n  star_mass_msun: 1\n"
            "thermodynamics_parameters:\n  model: radiative\n"
            "  mean_molecular_weight: 2.34\n  background_temperature_k: 20\n"
            "  minimum_temperature_k: 10\n  maximum_temperature_k: 1200\n"
            "  irradiation_luminosity_lsun: 1\n  irradiation_factor: 0.05\n"
        )
        self.assertEqual(call.args[1]["thermal_model"], "radiative")
        for key, flag in (
            ("thermal_model", "-thermal_model"),
            ("mean_molecular_weight", "-thermal_mu"),
            ("background_temperature_k", "-thermal_background"),
            ("minimum_temperature_k", "-thermal_tmin"),
            ("maximum_temperature_k", "-thermal_tmax"),
            ("irradiation_luminosity_lsun", "-irradiation_luminosity"),
            ("irradiation_factor", "-irradiation_factor"),
        ):
            self.assertIn(key, call.args[1])
            self.assertEqual(call.args[2][key], flag)


if __name__ == "__main__":
    unittest.main()
