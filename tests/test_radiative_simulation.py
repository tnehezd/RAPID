"""End-to-end tests; build bin/simulation with make all first."""

import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "bin" / "simulation"


class RadiativeSimulationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.counter = 0

    def run_simulation(self, *extra, binary=BINARY):
        self.counter += 1
        output = self.directory / f"run_{self.counter}"
        mass_arguments = [] if "-disk_mass" in extra else ["-sigma0_init", "1e-5"]
        command = [
            str(binary), "-drift", "0", "-growth", "0", "-twopop", "0",
            "-n", "48", "-ri", "1", "-ro", "10", *mass_arguments,
            "-alpha_init", "0.003", "-rdzei", "2.7", "-rdzeo", "6",
            "-drdzei", "1", "-drdzeo", "1", "-amod", "0.01",
            "-density_floor", "1e-20", "-tmax", "10", "-outfreq", "5",
            "-o", str(output), "-no_panels", *extra,
        ]
        process = subprocess.run(
            command, cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": "2"},
            capture_output=True, text=True, timeout=120,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        return output

    def check_thermal_snapshots(self, output):
        thermal_files = sorted((output / "LOGS").glob("temperature_profile_*.dat"))
        self.assertGreaterEqual(len(thermal_files), 2)
        for path in thermal_files:
            data = np.loadtxt(path)
            self.assertEqual(data.shape, (48, 9))
            self.assertTrue(np.all(np.isfinite(data)))
            self.assertTrue(np.all(data[:, 1] >= 10.0))
            self.assertTrue(np.all(data[:, 2:7] > 0.0))
            np.testing.assert_allclose(data[:, 7], data[:, 8], rtol=1e-8, atol=1e-20)
        return np.loadtxt(output / "LOGS" / "temperature_profile_0.dat")

    def test_radiative_evolves_on_both_grids_with_and_without_sg(self):
        for grid in ("linear", "logarithmic"):
            for sg in ("0", "1"):
                with self.subTest(grid=grid, sg=sg):
                    output = self.run_simulation(
                        "-thermal_model", "radiative", "-grid_type", grid, "-selfgravity", sg
                    )
                    initial = self.check_thermal_snapshots(output)
                    logs = output / "LOGS"
                    initial_gas = np.loadtxt(logs / "density_profile_00000000.dat")
                    final_gas = np.loadtxt(logs / "density_profile_00000010.dat")
                    self.assertTrue(np.all(final_gas[:, 2] > 0))
                    self.assertGreater(np.max(np.abs(final_gas[:, 2] / initial_gas[:, 2] - 1)),
                                       1e-8)
                    self.assertGreater(np.max(np.abs(
                        np.loadtxt(logs / "temperature_profile_10.dat")[:, 1] /
                        initial[:, 1] - 1)), 1e-9)
                    # Independent CGS -> code conversion and pressure reconstruction.
                    time_unit = 365.242199 * 86400 / (2 * math.pi)
                    radius = initial[:, 0]
                    height = initial[:, 2] * radius
                    cs_squared_cgs = 1.380649e-16 * initial[:, 1] / (2.34 * 1.6735575e-24)
                    cs_squared_code = cs_squared_cgs * (time_unit / 1.496e13) ** 2
                    expected_pressure = initial_gas[:, 2] * cs_squared_code / (
                        math.sqrt(2 * math.pi) * height
                    )
                    np.testing.assert_allclose(initial_gas[:, 3], expected_pressure, rtol=2e-9)
                    sigma_cgs = initial_gas[:, 2] * 1.989e33 / 1.496e13 ** 2
                    np.testing.assert_allclose(initial[:, 3] * math.sqrt(2 * math.pi) *
                                               height * 1.496e13, sigma_cgs, rtol=2e-9)
                    alpha = .003 * (1 - .5 * .99 * (
                        np.tanh(radius - 2.7) + np.tanh(6 - radius)
                    ))
                    np.testing.assert_allclose(
                        initial[:, 4], alpha * np.sqrt(cs_squared_cgs) * height * 1.496e13,
                        rtol=1e-12,
                    )
                    omega = np.sqrt(cs_squared_code) / height
                    extra = omega ** 2 - 1.0 / radius ** 3
                    derivative = np.gradient(extra, radius, edge_order=2)
                    shear = (-3.0 / radius ** 3 + radius * derivative) / (2 * omega ** 2)
                    expected_heating = initial[:, 4] * sigma_cgs * (
                        shear * omega / time_unit
                    ) ** 2
                    np.testing.assert_allclose(
                        initial[1:-1, 7], expected_heating[1:-1], rtol=1e-10
                    )

    def test_sg_changes_radiative_temperature(self):
        without = self.run_simulation("-thermal_model", "radiative", "-selfgravity", "0")
        with_sg = self.run_simulation("-thermal_model", "radiative", "-selfgravity", "1")
        left = self.check_thermal_snapshots(without)
        right = self.check_thermal_snapshots(with_sg)
        self.assertGreater(np.max(np.abs(right[:, 1] / left[:, 1] - 1)), 1e-6)

    def test_hdf5_run_also_writes_thermal_diagnostics(self):
        output = self.run_simulation("-thermal_model", "radiative", "-output_format", "hdf5")
        self.check_thermal_snapshots(output)
        files = sorted((output / "LOGS").glob("snapshot_*.h5"))
        self.assertGreaterEqual(len(files), 2)
        self.assertEqual(files[-1].read_bytes()[:8], b"\x89HDF\r\n\x1a\n")

    def test_zero_alpha_returns_background_and_no_evolution(self):
        output = self.run_simulation("-thermal_model", "radiative", "-alpha_init", "0")
        for path in (output / "LOGS").glob("temperature_profile_*.dat"):
            data = np.loadtxt(path)
            np.testing.assert_array_equal(data[:, 1], 10.0)
            np.testing.assert_array_equal(data[:, 4], 0.0)
            np.testing.assert_array_equal(data[:, 7:9], 0.0)

    def test_lociso_aliases_are_identical(self):
        outputs = [self.run_simulation(*args) for args in (
            (), ("-thermal_model", "lociso"), ("-thermal_model", "locally_isotherm"),
        )]
        reference = np.loadtxt(outputs[0] / "LOGS" / "density_profile_00000010.dat")
        for output in outputs:
            np.testing.assert_array_equal(
                reference, np.loadtxt(output / "LOGS" / "density_profile_00000010.dat")
            )
            self.assertEqual(list((output / "LOGS").glob("temperature_profile_*.dat")), [])

    def test_radiative_with_dust_and_photoevaporation(self):
        for extra in (
            ("-drift", "1", "-growth", "1", "-ndust", "32", "-largesize", "0.01"),
            ("-photoevap", "1", "-photoevap_mode", "Owen", "-xray_luminosity", "1e30"),
        ):
            with self.subTest(extra=extra):
                output = self.run_simulation("-thermal_model", "radiative", *extra)
                self.check_thermal_snapshots(output)

    def test_low_disk_mass_approaches_no_sg(self):
        profiles = []
        for sg in ("0", "1"):
            output = self.run_simulation(
                "-thermal_model", "radiative", "-sigma0_init", "1e-10", "-selfgravity", sg
            )
            profiles.append(self.check_thermal_snapshots(output))
        np.testing.assert_allclose(profiles[0][:, 1], profiles[1][:, 1], rtol=1e-6)

    def test_yaml_wrapper_runs_radiative_model(self):
        config = yaml.safe_load((ROOT / "configs" / "config_thermal.yaml").read_text())
        config["disk_parameters"]["number_of_grid_points"] = 64
        config["disk_parameters"]["radial_grid_type"] = "logarithmic"
        config["simulation_parameters"]["omp_num_threads"] = 2
        config["time_parameters"] = {"total_simulation_time": 2, "output_write_frequency": 1}
        output = self.directory / "yaml_run"
        config["output_parameters"]["output_directory_name"] = str(output)
        config["log_parameters"] = {"disable_terminal_panels": True}
        filename = self.directory / "radiative.yaml"
        filename.write_text(yaml.safe_dump(config))
        process = subprocess.run(
            [sys.executable, str(ROOT / "run_simulation.py"), "-c", str(filename)],
            cwd=ROOT, capture_output=True, text=True, timeout=120,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        data = np.loadtxt(output / "LOGS" / "temperature_profile_2.dat")
        self.assertEqual(data.shape, (64, 9))

    def test_bad_cli_and_unbracketed_equilibrium_fail(self):
        for arguments in (
            ("-thermal_model", "polytropic"),
            ("-thermal_model", "typo"),
            ("-thermal_mu", "nan"),
            ("-thermal_mu", "2garbage"),
            ("-thermal_background", "2000"),
            ("-thermal_tmin", "1500"),
            ("-thermal_model", "radiative", "-alpha_init", "-1"),
        ):
            with self.subTest(arguments=arguments):
                process = subprocess.run(
                    [str(BINARY), *arguments], cwd=ROOT, capture_output=True, text=True,
                    timeout=30,
                )
                self.assertNotEqual(process.returncode, 0)
        output = self.directory / "bad_bracket"
        process = subprocess.run(
            [str(BINARY), "-thermal_model", "radiative", "-thermal_tmax", "10.00001",
             "-drift", "0", "-growth", "0", "-twopop", "0", "-n", "48", "-ri", "1",
             "-ro", "10", "-disk_mass", "0.01", "-o", str(output), "-no_panels"],
            cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": "2"},
            capture_output=True, text=True, timeout=30,
        )
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("Radiative equilibrium failed", process.stdout + process.stderr)


if __name__ == "__main__":
    unittest.main()
