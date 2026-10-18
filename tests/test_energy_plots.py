import tempfile
from pathlib import Path
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np

from plot_temperature import plot_temperature
from plot_thermal_heatmap import plot_heatmaps


class EnergyPlotTests(unittest.TestCase):
    def test_signed_cooling_and_equilibrium_compatibility(self):
        for energy_mode in (False, True):
            with self.subTest(energy_mode=energy_mode), tempfile.TemporaryDirectory() as directory:
                run = Path(directory)
                logs = run / "LOGS"
                logs.mkdir()
                for time in (0, 1):
                    thermal = np.array([
                        [1, 40, .03, 1e-10, 0, 1, 10, 0, -2],
                        [2, 30, .04, 2e-10, 0, 2, 20, 0, 1],
                    ], dtype=float)
                    if not energy_mode:
                        thermal[:, 7:9] = 2
                    mode = "energy evolution" if energy_mode else "equilibrium"
                    np.savetxt(
                        logs / f"temperature_profile_{time:08d}.dat", thermal,
                        header=f"Radiative {mode} at {time} years; gamma=1.4",
                    )
                    gas = np.array([[1, 1, 1e-5, 1e-8, -1e-9],
                                    [2, 2, 2e-5, 2e-8, 1e-9]])
                    np.savetxt(logs / f"density_profile_{time:08d}.dat", gas,
                               header=f"Gas at {time} years")
                profile = plot_temperature(logs / "temperature_profile_00000001.dat")
                heatmap = plot_heatmaps(run)
                try:
                    norm = heatmap.axes[7].collections[0].norm
                    self.assertIsInstance(norm, SymLogNorm if energy_mode else LogNorm)
                    if energy_mode:
                        self.assertEqual(profile.axes[4].get_yscale(), "symlog")
                        np.testing.assert_array_equal(
                            profile.axes[5].lines[0].get_ydata(), [2, -1]
                        )
                    profile.savefig(run / "profile.png")
                    heatmap.savefig(run / "heatmap.png")
                finally:
                    plt.close(profile)
                    plt.close(heatmap)


if __name__ == "__main__":
    unittest.main()
