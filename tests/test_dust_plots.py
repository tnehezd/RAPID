import tempfile
from pathlib import Path
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_dust import plot_dust


class DustPlotTests(unittest.TestCase):
    def write_snapshots(self, run):
        logs = run / "LOGS"
        logs.mkdir()
        for label, time in ((0, 0), (10, 10.25)):
            gas = [[1, 1, 1e-5, 1e-8, -1e-9], [2, 2, 2e-5, 2e-8, -1e-9]]
            np.savetxt(logs / f"density_profile_{label:08d}.dat", gas,
                       header=f"Gas at {time} years")
            np.savetxt(logs / f"dust_density_profile_{label:08d}.dat",
                       [[label, 1, 0], [label, 2, 1e-7]])
            particles = [[label, 7, 1.5, .1], [label, 2, 1.8, .2]]
            np.savetxt(logs / f"dust_size_evolution_{label:08d}.dat", particles,
                       header=f"Particles at {time / (2 * np.pi)} years")

    def test_profiles_and_consistent_particle_id_colours(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for linear in (False, True):
                with self.subTest(linear=linear):
                    figure = plot_dust(run, linear=linear)
                    try:
                        self.assertEqual(figure.axes[0].get_yscale(), "linear" if linear else "log")
                        points = figure.axes[2].collections
                        np.testing.assert_array_equal(points[0].get_array(), [7, 2])
                        np.testing.assert_array_equal(points[1].get_offsets()[:, 1], [10.25, 10.25])
                        self.assertIs(points[0].norm, points[1].norm)
                        figure.savefig(run / "overview.png")
                        self.assertGreater((run / "overview.png").stat().st_size, 0)
                    finally:
                        plt.close(figure)

    def test_empty_particles_and_shifted_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for path in (run / "LOGS").glob("dust_size_evolution_*.dat"):
                path.write_text("# No active particles\n")
            with self.assertWarnsRegex(UserWarning, "No active particles"):
                figure = plot_dust(run)
            plt.close(figure)
            path = run / "LOGS" / "dust_density_profile_00000010.dat"
            np.savetxt(path, [[10, .9, 1e-7], [10, 1, 2e-7]])
            with self.assertRaisesRegex(ValueError, "grids differ"):
                plot_dust(run)


if __name__ == "__main__":
    unittest.main()
