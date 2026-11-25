import tempfile
from pathlib import Path
import unittest

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_output import load_snapshots, plot_output


class HDF5PlotTests(unittest.TestCase):
    def write_snapshots(self, run, populations=()):
        logs = run / "LOGS"
        logs.mkdir()
        for index in range(2):
            with h5py.File(logs / f"snapshot_{index:08d}.h5", "w") as source:
                source["frame/time"] = [index * 10.]
                source["gas_grid/radial_grid"] = [1., 2., 3.]
                source["gas_grid/surface_density"] = np.array([1., 2., 3.]) * 1e-5
                for component in populations:
                    group = "dust_grid" if component == "dust" else "micron_grid"
                    source[group + "/surface_density"] = [1e-7, 0, 3e-7]
                if "dust" in populations:
                    source["particles/index"] = [0, 7, 9]
                    source["particles/position"] = [0, 1.5, 2.5]
                    source["particles/size"] = [0, .1, .2]

    def test_gas_single_and_two_populations(self):
        for populations in ((), ("dust",), ("dust", "micron")):
            with self.subTest(populations=populations), tempfile.TemporaryDirectory() as directory:
                run = Path(directory)
                self.write_snapshots(run, populations)
                snapshots = load_snapshots(run)
                self.assertEqual(set(snapshots[0].density), {"gas", *populations})
                np.testing.assert_array_equal(snapshots[0].radius, [1, 2, 3])
                if populations:
                    np.testing.assert_array_equal(snapshots[0].particles["dust"][:, 0], [7, 9])
                for kind in ("profiles", "heatmap"):
                    if "micron" in populations and kind == "profiles":
                        with self.assertWarnsRegex(UserWarning, "not saved"):
                            figure = plot_output(run, kind=kind)
                    else:
                        figure = plot_output(run, kind=kind)
                    try:
                        figure.savefig(run / (kind + ".png"))
                        self.assertGreater((run / (kind + ".png")).stat().st_size, 0)
                    finally:
                        plt.close(figure)

    def test_shape_and_ambiguous_format(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            path = run / "LOGS" / "snapshot_00000000.h5"
            with h5py.File(path, "a") as source:
                del source["gas_grid/surface_density"]
                source["gas_grid/surface_density"] = [1, 2]
            with self.assertRaisesRegex(ValueError, "shape"):
                load_snapshots(run)
            (run / "LOGS" / "density_profile_00000000.dat").write_text("# Initial gas profile\n")
            with self.assertRaisesRegex(ValueError, "Both ASCII and HDF5"):
                load_snapshots(run)

    def test_single_snapshot_heatmap_error(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            path = run / "LOGS" / "snapshot_00000000.h5"
            with self.assertRaisesRegex(ValueError, "at least two"):
                plot_output(path, kind="heatmap")

    def test_two_population_initial_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run, ("dust", "micron"))
            with h5py.File(run / "LOGS" / "snapshot_00000001.h5", "a") as source:
                source["dust_grid/surface_density"][:] = [2e-7, 1e-7, 6e-7]
                source["micron_grid/surface_density"][:] = [5e-8, 0, 1.5e-7]
            with self.assertWarnsRegex(UserWarning, "undefined ratios"):
                figure = plot_output(run, kind="heatmap", components=["dust", "micron"],
                                     normalize_initial=True)
            try:
                for axis, ratio in zip(figure.axes[:2], (2, .5)):
                    values = axis.collections[0].get_array()
                    np.testing.assert_allclose(values[[0, 2]], [[1, ratio], [1, ratio]])
                    self.assertTrue(np.all(np.ma.getmaskarray(values)[1]))
                figure.canvas.draw()
            finally:
                plt.close(figure)


if __name__ == "__main__":
    unittest.main()
