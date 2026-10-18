import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "benchmark_runner", ROOT / "test_examples" / "generate_configs.py"
)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class BenchmarkRunnerTests(unittest.TestCase):
    def test_full_matrix_without_templates(self):
        with tempfile.TemporaryDirectory() as temporary:
            paths = runner.generate_configs(Path(temporary))
            self.assertEqual(len(paths), 80)
            self.assertEqual(len(set(paths)), 80)
            combinations = set()
            outputs = set()
            for path in paths:
                config = yaml.safe_load(path.read_text())
                simulation = config["simulation_parameters"]
                thermal = config["thermodynamics_parameters"]
                model = thermal["model"]
                combinations.add((
                    simulation["enable_selfgravity"],
                    simulation["enable_photoevaporation"],
                    simulation["enable_dust_drift"],
                    config["disk_parameters"]["radial_grid_type"],
                    model, thermal.get("irradiation_model"),
                ))
                outputs.add(config["output_parameters"]["output_directory_name"])
                self.assertEqual(simulation["omp_num_threads"], 4)
                self.assertFalse(simulation["enable_two_dust_populations"])
                self.assertEqual(simulation["enable_dust_growth"],
                                 simulation["enable_dust_drift"])
                if simulation["enable_dust_drift"]:
                    self.assertEqual(simulation["number_of_dust_particles"], 1000)
                    self.assertEqual(config["dust_parameters"]["one_size_particle_value_cm"], .01)
                boundary = "outflow" if model == "radiative_energy" else "zero_gradient"
                self.assertEqual(set(config["boundary_conditions"].values()), {boundary})
            self.assertEqual(len(combinations), 80)
            self.assertEqual(len(outputs), 80)

    def test_launch_environment_and_logs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "case.yaml"

            def launch(command, **kwargs):
                self.assertEqual(command, [sys.executable, str(root / "run_simulation.py"),
                                           "-c", str(path)])
                for variable in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT"):
                    self.assertEqual(kwargs["env"][variable], "4")
                for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
                    self.assertEqual(kwargs["env"][variable], "1")
                kwargs["stdout"].write("simulator diagnostic\n")
                self.assertEqual(kwargs["stderr"], subprocess.STDOUT)
                return subprocess.CompletedProcess(command, 7)

            with patch.object(runner.subprocess, "run", side_effect=launch):
                name, code, seconds, log = runner.run_case(path, root, root)
            self.assertEqual((name, code), ("case.yaml", 7))
            self.assertGreaterEqual(seconds, 0)
            self.assertIn("simulator diagnostic", log.read_text())

    def test_batch_thread_bound_and_failure_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "test_examples"
            directory.mkdir()
            (directory.parent / "bin").mkdir()
            (directory.parent / "bin" / "simulation").touch()
            paths = [directory / f"case_{i}.yaml" for i in range(25)]
            lock = threading.Lock()
            active = peak = 0

            def run_case(path, root, logs):
                nonlocal active, peak
                with lock:
                    active += 1
                    peak = max(peak, active)
                time.sleep(.02)
                with lock:
                    active -= 1
                return path.name, int(path == paths[0]), .02, logs / f"{path.stem}.log"

            with patch.object(runner, "__file__", str(directory / "generate_configs.py")), \
                    patch.object(runner, "generate_configs", return_value=paths), \
                    patch.object(runner, "run_case", side_effect=run_case), \
                    patch.object(sys, "argv", ["generate_configs.py"]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(runner.main(), 1)
            self.assertGreater(peak, 1)
            self.assertLessEqual(peak, 10)
            rows = (directory.parent / "benchmark_logs" / "timings.tsv").read_text().splitlines()
            self.assertEqual(len(rows), 26)
            self.assertEqual(sum(row.split("\t")[1] == "1" for row in rows[1:]), 1)

    def test_generate_only_does_not_launch(self):
        with patch.object(runner, "generate_configs", return_value=[]), \
                patch.object(runner, "run_case") as launch, \
                patch.object(sys, "argv", ["generate_configs.py", "--generate-only"]), \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(runner.main(), 0)
        launch.assert_not_called()

    def test_rejects_more_than_ten_jobs(self):
        with patch.object(sys, "argv", ["generate_configs.py", "--jobs", "11"]), \
                contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            runner.main()
        self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
