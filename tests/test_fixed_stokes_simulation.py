import bisect
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "bin" / "release" / "simulation"
AU_CM = 1.496e13
SOLAR_MASS = 1.989e33


def rows(path):
    with path.open() as stream:
        return [[float(value) for value in line.split()] for line in stream
                if line.strip() and not line.startswith("#")]


def interpolate(x, xs, ys):
    right = min(max(bisect.bisect_right(xs, x), 1), len(xs) - 1)
    left = right - 1
    return ys[left] + (ys[right] - ys[left]) * (x - xs[left]) / (xs[right] - xs[left])


class FixedStokesSimulationTests(unittest.TestCase):
    def run_binary(self, *args):
        return subprocess.run(
            [str(BINARY), *map(str, args)], cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1"}, timeout=30,
        )

    def check_stokes(self, gas, particles, tolerance=2e-6):
        radial = [row[0] for row in gas]
        density = [row[1] for row in gas]
        checked = 0
        for position, size_cm in particles:
            if position <= radial[0] or position >= radial[-1]:
                continue
            sigma = interpolate(position, radial, density)
            stokes = (1.6 * AU_CM ** 3 / SOLAR_MASS * (size_cm / AU_CM)
                      * math.pi / (2 * sigma))
            self.assertAlmostEqual(stokes, 0.1, delta=tolerance)
            checked += 1
        self.assertGreater(checked, 5)

    def test_ascii_evolving_and_static(self):
        for evolution in (0, 1):
            with self.subTest(evolution=evolution), tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / "run"
                result = self.run_binary(
                    "-n", 30, "-ndust", 20, "-ri", 1, "-ro", 5,
                    "-disk_mass", 0.01, "-alpha_init", 0.001,
                    "-stokes", 0.1, "-drift", 1, "-growth", 0, "-twopop", 0,
                    "-evol", evolution, "-tmax", 2, "-outfreq", 2,
                    "-no_panels", "-o", output,
                )
                self.assertEqual(result.returncode, 0, result.stdout)
                initial = rows(output / "config" / "initial_dust_profile.dat")
                self.assertTrue(all(row[2] > 0 and row[3] == 0 for row in initial))
                diagnostic = rows(output / "drift_timecale.dat")
                self.assertTrue(all(abs(row[4] - 0.1) < 1e-12 for row in diagnostic))
                self.assertIn("Stokes number:     0.10000000000000001",
                              (output / "current_run_info.dat").read_text())
                initial_gas = rows(output / "config" / "initial_gas_profile.dat")
                for time in (0, 2):
                    if evolution:
                        gas_rows = rows(output / "LOGS" / f"density_profile_{time:08d}.dat")
                    else:
                        gas_rows = initial_gas
                    gas = [(row[1], row[2]) for row in gas_rows]
                    dust = rows(output / "LOGS" / f"dust_size_evolution_{time:08d}.dat")
                    self.check_stokes(gas, [(row[2], row[3]) for row in dust])
                before = rows(output / "LOGS" / "dust_size_evolution_00000000.dat")
                after = rows(output / "LOGS" / "dust_size_evolution_00000002.dat")
                before_by_id = {row[1]: row for row in before}
                survivors = {int(row[1]) for row in after}
                self.assertLess(len(survivors), len(initial))
                expected_mass = sum(row[2] for row in initial if int(row[0]) in survivors)
                density = rows(output / "LOGS" / "dust_density_profile_00000002.dat")
                dr = 4 / 29
                deposited_mass = sum(row[2] * 2 * math.pi * row[1] * dr for row in density)
                self.assertAlmostEqual(deposited_mass / expected_mass, 1, delta=3e-6)
                self.assertTrue(any(abs(row[2] - before_by_id[row[1]][2]) > 1e-5
                                    and abs(row[3] - before_by_id[row[1]][3]) > 1e-5
                                    for row in after if row[1] in before_by_id))

    def test_rejected_options(self):
        common = ("-stokes", "0.1", "-drift", "1", "-growth", "0", "-twopop", "0")
        for args in (
            ("-selfgravity", "1"), ("-backreaction", "1"),
            ("-thermal_model", "radiative"), ("-thermal_model", "radiative_energy"),
            ("-selfgravity", "true"), ("-stokes", "nan"), ("-stokes", "0"),
            (*common, "-growth", "1"), (*common, "-twopop", "1"),
            (*common, "-drift", "0"), (*common, "-largesize", "1"),
        ):
            with self.subTest(args=args):
                result = self.run_binary(*args)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn("ERROR", result.stdout)

    def test_yaml_forwards_stokes(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "run"
            config = Path(temp) / "input.yaml"
            config.write_text(json.dumps({
                "simulation_parameters": {
                    "enable_dust_drift": True, "enable_dust_growth": False,
                    "enable_two_dust_populations": False, "enable_gas_evolution": True,
                    "number_of_dust_particles": 20,
                },
                "disk_parameters": {
                    "number_of_grid_points": 30, "inner_radius_au": 1,
                    "outer_radius_au": 5, "disk_mass": .01,
                },
                "dust_parameters": {"constant_stokes_number": .1},
                "time_parameters": {"total_simulation_time": 2, "output_write_frequency": 2},
                "output_parameters": {"output_directory_name": str(output)},
                "log_parameters": {"disable_terminal_panels": True},
            }))
            result = subprocess.run(
                [sys.executable, str(ROOT / "run_simulation.py"), "-c", str(config),
                 "-exec", str(BINARY)],
                cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": "1"},
                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("-stokes 0.1", result.stdout)
            self.assertTrue(all(abs(row[4] - .1) < 1e-12
                                for row in rows(output / "drift_timecale.dat")))

    def test_yaml_rejects_disabled_physics(self):
        configs = (
            {"simulation_parameters": {"enable_selfgravity": True}},
            {"simulation_parameters": {"enable_backreaction": True}},
            {"thermodynamics_parameters": {"model": "radiative"}},
            {"thermodynamics_parameters": {"model": "radiative_energy"}},
        )
        for config in configs:
            with self.subTest(config=config), tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / "input.yaml"
                config["disk_parameters"] = {"star_mass_msun": 1}
                path.write_text(json.dumps(config))
                result = subprocess.run(
                    [sys.executable, str(ROOT / "run_simulation.py"), "-c", str(path),
                     "-exec", str(BINARY)],
                    cwd=ROOT, env={**os.environ, "OMP_NUM_THREADS": "1"},
                    text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30,
                )
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn("pending physical validation", result.stdout.lower())

    @unittest.skipUnless(shutil.which("h5dump"), "h5dump is required to inspect HDF5")
    def test_hdf5_physical_grid_positions_and_sizes(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "run"
            result = self.run_binary(
                "-n", 16, "-ndust", 12, "-ri", 1, "-ro", 5,
                "-disk_mass", .01, "-alpha_init", .001, "-stokes", .1,
                "-drift", 1, "-growth", 0, "-twopop", 0, "-evol", 1,
                "-tmax", 2, "-outfreq", 2, "-output_format", "hdf5",
                "-no_panels", "-o", output,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            snapshot = output / "LOGS" / "snapshot_00000002.h5"

            def dataset(name):
                text = subprocess.check_output(
                    ["h5dump", "-y", "-w", "0", "-m", "%.17g", "-d", name, str(snapshot)],
                    text=True,
                )
                data = re.search(r"DATA\s*\{([^}]*)\}", text)
                self.assertIsNotNone(data, text)
                return [float(value) for value in data.group(1).replace("\n", "").split(",")]

            grid = dataset("/gas_grid/radial_grid")
            sigma = dataset("/gas_grid/surface_density")
            positions = dataset("/particles/position")
            sizes = dataset("/particles/size")
            self.assertAlmostEqual(grid[0], 1)
            self.assertAlmostEqual(grid[-1], 5)
            self.assertEqual(len(positions), 12)
            self.check_stokes(list(zip(grid, sigma)), list(zip(positions, sizes)))


if __name__ == "__main__":
    unittest.main()
