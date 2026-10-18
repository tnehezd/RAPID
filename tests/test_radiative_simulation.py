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
            self.assertTrue(np.all(data[:, [2, 3, 5, 6]] > 0.0))
            self.assertTrue(np.all(data[:, 4] >= 0.0))
            np.testing.assert_allclose(data[:, 7], data[:, 8], rtol=1e-8, atol=1e-20)
        for path in thermal_files:
            self.assertRegex(path.name, r"^temperature_profile_\d{8}\.dat$")
        return np.loadtxt(output / "LOGS" / "temperature_profile_00000000.dat")

    def test_thermodynamics_startup_panel(self):
        for model, sg in (("lociso", "0"), ("lociso", "1"),
                          ("radiative", "0"), ("radiative", "1")):
            with self.subTest(model=model, sg=sg):
                output = self.directory / f"panel_{model}_{sg}"
                process = subprocess.run(
                    [str(BINARY), "-drift", "0", "-growth", "0", "-twopop", "0",
                     "-n", "16", "-ri", "1", "-ro", "10", "-sigma0_init", "1e-5",
                     "-tmax", "0", "-o", str(output),
                     "-thermal_model", model, "-selfgravity", sg,
                     "-thermal_mu", "2.5", "-thermal_background", "12",
                     "-thermal_tmin", "2", "-thermal_tmax", "1400"],
                    cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": "2"},
                    capture_output=True, text=True, timeout=120,
                )
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                panel = process.stderr.split("THERMODYNAMICS", 1)[1].split(
                    "SELF-GRAVITY", 1
                )[0]
                self.assertNotIn("SG COUPLING", panel)
                sg_panel = process.stderr.split("SELF-GRAVITY", 1)[1].split(
                    "DUST SMOOTHING METHOD", 1
                )[0]
                self.assertIn(f"SELF-GRAVITY:   [{'ON' if sg == '1' else 'OFF'}]",
                              sg_panel)
                if sg == "1":
                    self.assertIn("[Axisymmetric disk rotation correction]", sg_panel)
                    self.assertIn("VERTICAL SG:    [Not included]", sg_panel)
                if sg == "1" and model == "radiative":
                    self.assertIn("[Rotation and shear; H-dependent softening]", sg_panel)
                else:
                    self.assertNotIn("THERMAL COUPL.", sg_panel)
                if model == "radiative":
                    for text in (
                        "[Radiative equilibrium]", "[Instantaneous; no energy evolution]",
                        "[Viscous]", "[Opacity-dependent radiation]",
                        "[Smoothed Lin-Papaloizou]", "STELLAR IRRAD.: [OFF]",
                        "[2.5 m_H]", "[12 K]", "[2 -- 1400 K]",
                    ):
                        self.assertIn(text, panel)
                    self.assertNotIn("H/R AT 1 AU", panel)
                else:
                    self.assertIn("[Locally isothermal]", panel)
                    self.assertIn("H/R AT 1 AU", panel)
                    self.assertIn("FLARING INDEX", panel)
                    self.assertIn("THERMAL SOLVER: [OFF]", panel)
                    self.assertNotIn("BACKGROUND T", panel)

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
                        np.loadtxt(logs / "temperature_profile_00000010.dat")[:, 1] /
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

    def test_irradiation_only_matches_prescribed_temperature(self):
        for grid in ("linear", "logarithmic"):
            for sg in ("0", "1"):
                with self.subTest(grid=grid, sg=sg):
                    output = self.run_simulation(
                        "-thermal_model", "radiative", "-alpha_init", "0",
                        "-irradiation_luminosity", "1", "-irradiation_factor", "0.05",
                        "-grid_type", grid, "-selfgravity", sg,
                    )
                    data = self.check_thermal_snapshots(output)
                    irradiated_t4 = .05 * 3.828e33 / (
                        4 * math.pi * 5.670374419e-5 * (data[:, 0] * 1.496e13)**2
                    )
                    expected = (10**4 + irradiated_t4)**.25
                    np.testing.assert_allclose(data[:, 1], expected, rtol=1e-14)
                    np.testing.assert_array_equal(data[:, 7:9], 0.0)
                    self.assertIn("Irradiation:", (
                        output / "LOGS" / "temperature_profile_00000000.dat"
                    ).read_text())

    def test_irradiated_viscous_equilibrium_and_off_compatibility(self):
        baseline = self.run_simulation("-thermal_model", "radiative")
        for disabled in (("-irradiation_luminosity", "0"),
                         ("-irradiation_luminosity", "1", "-irradiation_factor", "0")):
            output = self.run_simulation("-thermal_model", "radiative", *disabled)
            for name in ("density_profile_00000010.dat", "temperature_profile_00000010.dat"):
                np.testing.assert_array_equal(
                    np.loadtxt(baseline / "LOGS" / name),
                    np.loadtxt(output / "LOGS" / name),
                )
        for sg in ("0", "1"):
            output = self.run_simulation(
                "-thermal_model", "radiative", "-selfgravity", sg,
                "-irradiation_luminosity", "1", "-irradiation_factor", "0.05",
            )
            data = self.check_thermal_snapshots(output)
            t_rad4 = 10**4 + .05 * 3.828e33 / (
                4 * math.pi * 5.670374419e-5 * (data[:, 0] * 1.496e13)**2
            )
            self.assertTrue(np.all(data[:, 1] > t_rad4**.25))
            tau = data[:, 6]
            tau_eff = .375 * tau + math.sqrt(3) / 4 + .25 / tau
            emitted_net = 2 * 5.670374419e-5 * (data[:, 1]**4 - t_rad4) / tau_eff
            np.testing.assert_allclose(data[:, 7], emitted_net, rtol=1e-5)
            self.assertTrue(np.all(data[:, 1] > np.loadtxt(
                baseline / "LOGS" / "temperature_profile_00000000.dat"
            )[:, 1]))

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
        config["thermodynamics_parameters"]["irradiation_luminosity_lsun"] = 1
        config["thermodynamics_parameters"]["irradiation_factor"] = 0.05
        filename = self.directory / "radiative.yaml"
        filename.write_text(yaml.safe_dump(config))
        process = subprocess.run(
            [sys.executable, str(ROOT / "run_simulation.py"), "-c", str(filename)],
            cwd=ROOT, capture_output=True, text=True, timeout=120,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        data = np.loadtxt(output / "LOGS" / "temperature_profile_00000002.dat")
        self.assertEqual(data.shape, (64, 9))
        self.assertIn("STELLAR IRRAD.: [ON / prescribed effective temperature]",
                      process.stdout + process.stderr)
        self.assertIn("Irradiation:       ENABLED", (output / "current_run_info.dat").read_text())

    def test_bad_cli_and_unbracketed_equilibrium_fail(self):
        for arguments in (
            ("-thermal_model", "polytropic"),
            ("-thermal_model", "typo"),
            ("-thermal_mu", "nan"),
            ("-thermal_mu", "2garbage"),
            ("-thermal_background", "2000"),
            ("-thermal_tmin", "1500"),
            ("-thermal_model", "radiative", "-alpha_init", "-1"),
            ("-irradiation_luminosity", "-1"),
            ("-irradiation_luminosity", "nan"),
            ("-irradiation_factor", "1.1"),
            ("-irradiation_factor", "0.1junk"),
            ("-irradiation_factor", "-0.05"),
            ("-irradiation_luminosity",),
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
