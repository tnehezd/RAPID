"""End-to-end stationary drag tests; build bin/simulation with make all first."""

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


class BackreactionSimulationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.counter = 0

    def run_case(self, *extra, enabled=True, threads=2, success=True):
        self.counter += 1
        output = self.directory / f"run_{self.counter}"
        command = [
            str(BINARY), "-no_panels", "-drift", "1", "-growth", "0", "-twopop", "0",
            "-evol", "1", "-n", "32", "-ndust", "64", "-ri", "1", "-ro", "10",
            "-sigma0_init", "1e-5", "-alpha_init", ".003", "-rdzei", "2.7",
            "-rdzeo", "6", "-drdzei", "1", "-drdzeo", "1", "-amod", ".01",
            "-density_floor", "1e-20", "-dust_smoothing", "cic", "-eps", ".01",
            "-largesize", "10", "-tmax", "10", "-outfreq", "5",
            "-inner_bc", "reflecting", "-outer_bc", "reflecting",
            "-o", str(output),
        ]
        if enabled:
            command.extend(["-backreaction", "1"])
        process = subprocess.run(
            [*command, *extra], cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": str(threads)},
            capture_output=True, text=True, timeout=60,
        )
        if success:
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            return output
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("ERROR", process.stdout + process.stderr)
        return process

    def final_gas(self, output):
        files = sorted((output / "LOGS").glob("density_profile_*.dat"))
        self.assertGreaterEqual(len(files), 2)
        return np.loadtxt(files[-1])

    def test_default_off_preserves_outputs(self):
        default = self.run_case(enabled=False)
        explicit = self.run_case("-backreaction", "0", enabled=False)
        for prefix in ("density_profile", "dust_size_evolution"):
            left = sorted((default / "LOGS").glob(f"{prefix}_*.dat"))
            right = sorted((explicit / "LOGS").glob(f"{prefix}_*.dat"))
            self.assertEqual(len(left), len(right))
            for a, b in zip(left, right):
                np.testing.assert_array_equal(np.loadtxt(a), np.loadtxt(b))

    def test_drag_changes_gas_transport_on_both_grids(self):
        for grid in ("linear", "logarithmic"):
            with self.subTest(grid=grid):
                light = self.run_case("-grid_type", grid, "-eps", "1e-6")
                loaded = self.run_case("-grid_type", grid, "-eps", "1")
                light_gas = self.final_gas(light)
                loaded_gas = self.final_gas(loaded)
                self.assertTrue(np.all(np.isfinite(loaded_gas)))
                self.assertTrue(np.all(loaded_gas[:, 2] >= 1e-20))
                self.assertGreater(np.max(abs(light_gas[:, 2] - loaded_gas[:, 2])), 1e-12)
                low_drift = np.loadtxt(light / "drift_timecale.dat")
                high_drift = np.loadtxt(loaded / "drift_timecale.dat")
                self.assertLess(np.median(abs(high_drift[4:-4, 5])),
                                np.median(abs(low_drift[4:-4, 5])))
                info = (loaded / "current_run_info.dat").read_text()
                self.assertIn("Dust Backreaction", info)
                self.assertIn("Stationary mutual drag", info)

    def test_dust_free_velocity_limit(self):
        passive = self.run_case("-eps", "1e-10", enabled=False)
        feedback = self.run_case("-eps", "1e-10")
        before = np.loadtxt(passive / "drift_timecale.dat")
        after = np.loadtxt(feedback / "drift_timecale.dat")
        # Legacy startup reloads rounded pressure; feedback reconstructs it from
        # Sigma. Differentiation amplifies the seven-significant-digit input error.
        np.testing.assert_allclose(after[4:-4, 5], before[4:-4, 5], rtol=3e-6)

    def test_growth_and_thread_determinism(self):
        for count in ("16", "32", "96"):
            with self.subTest(count=count):
                one = self.run_case("-growth", "1", "-ndust", count, threads=1)
                four = self.run_case("-growth", "1", "-ndust", count, threads=4)
                np.testing.assert_array_equal(self.final_gas(one), self.final_gas(four))
                path = sorted((four / "LOGS").glob("dust_size_evolution_*.dat"))[-1]
                grains = np.loadtxt(path)
                active = grains[:, 2] > 0
                self.assertTrue(np.any(active))
                self.assertTrue(np.all(grains[active, 3] > 0))
                self.assertTrue(np.all(np.isfinite(grains)))
                dust_path = sorted((four / "LOGS").glob("dust_density_profile_*.dat"))[-1]
                dust = np.loadtxt(dust_path)
                self.assertEqual(dust.shape, (32, 3))
                np.testing.assert_allclose(dust[:, 1], self.final_gas(four)[:, 1],
                                           rtol=1e-5, atol=1e-5)

    def test_inviscid_pressure_feedback(self):
        for grid in ("linear", "logarithmic"):
            with self.subTest(grid=grid):
                output = self.run_case("-alpha_init", "0", "-eps", "1",
                                       "-grid_type", grid, "-tmax", "100")
                files = sorted((output / "LOGS").glob("density_profile_*.dat"))
                initial = np.loadtxt(files[0])
                final = self.final_gas(output)
                self.assertTrue(np.all(np.isfinite(final)))
                self.assertTrue(np.all(final[:, 2] >= 1e-20))
                self.assertGreater(np.max(abs(final[:, 2] - initial[:, 2])), 1e-12)

    def test_rejects_unsupported_configurations(self):
        for arguments in (
            ("-thermal_model", "radiative"), ("-thermal_model", "radiative_energy"),
            ("-selfgravity", "1"), ("-photoevap", "1"), ("-twopop", "1"),
            ("-drift", "0"), ("-evol", "0"), ("-dust_smoothing", "gaussian"),
            ("-inner_bc", "absorbing"), ("-outer_bc", "fixed_flux"),
            ("-alpha_init", "-.01"), ("-amod", "-1"), ("-pdensity", "0"),
            ("-eps", "-.1"), ("-h_init", "0"), ("-density_floor", "nan"),
            ("-backreaction", "nan"), ("-backreaction", "2"), ("-backreaction",),
        ):
            with self.subTest(arguments=arguments):
                self.run_case(*arguments, success=False)

    def test_yaml_wrapper_and_hdf5_smoke(self):
        output = self.directory / "yaml_output"
        config = yaml.safe_load((ROOT / "configs" / "config_backreaction.yaml").read_text())
        config["simulation_parameters"].update(number_of_dust_particles=64, omp_num_threads=2)
        config["disk_parameters"].update(number_of_grid_points=32, outer_radius_au=10)
        config["boundary_conditions"] = {"inner_boundary_condition": "reflecting",
                                         "outer_boundary_condition": "reflecting"}
        config["time_parameters"] = {"total_simulation_time": 1, "output_write_frequency": 1}
        config["output_parameters"] = {"output_directory_name": str(output), "output_format": "hdf5"}
        path = self.directory / "feedback.yaml"
        path.write_text(yaml.safe_dump(config))
        process = subprocess.run(
            [sys.executable, str(ROOT / "run_simulation.py"), "-c", str(path)],
            cwd=ROOT, capture_output=True, text=True, timeout=60,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        self.assertIn("-backreaction 1", process.stdout)
        self.assertTrue(any(output.rglob("*.hdf5")) or any(output.rglob("*.h5")))
        self.assertEqual(np.loadtxt(output / "drift_timecale.dat").shape, (64, 8))


if __name__ == "__main__":
    unittest.main()
