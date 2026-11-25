import tempfile
from pathlib import Path
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_output import plot_output, load_snapshots


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
                    figure = plot_output(run, linear=linear)
                    try:
                        self.assertEqual(figure.axes[0].get_yscale(), "linear" if linear else "log")
                        points = figure.axes[2].collections
                        np.testing.assert_array_equal(points[0].get_array(), [7, 2])
                        np.testing.assert_array_equal(points[1].get_offsets()[:, 0], [10.25, 10.25])
                        np.testing.assert_array_equal(points[1].get_offsets()[:, 1], [1.5, 1.8])
                        self.assertEqual(figure.axes[2].get_xlabel(), "t [yr]")
                        self.assertEqual(figure.axes[2].get_ylabel(), "r [AU]")
                        self.assertEqual(figure.axes[2].get_xscale(), "linear")
                        self.assertEqual(figure.axes[2].get_yscale(), "linear" if linear else "log")
                        self.assertEqual(figure.axes[0].get_xlabel(), "r [AU]")
                        self.assertFalse(figure.axes[0].get_shared_x_axes().joined(
                            figure.axes[0], figure.axes[2]))
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
                figure = plot_output(run)
            plt.close(figure)
            path = run / "LOGS" / "dust_density_profile_00000010.dat"
            np.savetxt(path, [[10, .9, 1e-7], [10, 1, 2e-7]])
            with self.assertRaisesRegex(ValueError, "grids differ"):
                plot_output(run)

    def test_missing_time_is_error(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            path = run / "LOGS" / "density_profile_00000000.dat"
            np.savetxt(path, [[1, 1, 1e-5, 1e-8, -1e-9],
                              [2, 2, 2e-5, 2e-8, -1e-9]])
            with self.assertRaisesRegex(ValueError, "missing snapshot time"):
                plot_output(run)

    def test_initial_header_and_gas_only(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for pattern in ("dust_density_profile_*.dat", "dust_size_evolution_*.dat"):
                for path in (run / "LOGS").glob(pattern):
                    path.unlink()
            path = run / "LOGS" / "density_profile_00000000.dat"
            path.write_text(path.read_text().replace("Gas at 0 years", "Initial gas profile"))
            figure = plot_output(run)
            try:
                self.assertEqual(len(figure.axes), 2)
                self.assertIn("Initial: 0 yr", figure.axes[0].lines[-2].get_label())
                figure.savefig(run / "gas.png")
            finally:
                plt.close(figure)

    def test_twopop_ascii_heatmap_and_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            logs = run / "LOGS"
            for suffix in ("00000000.dat", "00000010.dat"):
                np.savetxt(logs / ("dust_micron_density_profile_" + suffix),
                           [[0, 1, 1e-8], [0, 2, 2e-8]])
                np.savetxt(logs / ("micron_dust_size_evolution_" + suffix),
                           [[0, 1, 1.5, 1e-4]])
            snapshots = load_snapshots(run)
            self.assertEqual(set(snapshots[0].density), {"gas", "dust", "micron"})
            for kind in ("profiles", "heatmap"):
                figure = plot_output(run, kind=kind, components=["micron"], stride=3)
                try:
                    figure.savefig(run / (kind + ".png"))
                    if kind == "profiles":
                        self.assertEqual(figure.axes[1].collections[0].get_array()[0], 1)
                finally:
                    plt.close(figure)
            with self.assertRaisesRegex(ValueError, "unavailable"):
                plot_output(run, components=["unknown"])
            with self.assertRaisesRegex(ValueError, "stride"):
                plot_output(run, stride=0)

    def test_missing_population_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            (run / "LOGS" / "dust_density_profile_00000010.dat").unlink()
            with self.assertRaises(OSError):
                load_snapshots(run)

    def test_two_population_profiles_use_two_by_two_layout(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for label in (0, 10):
                logs = run / "LOGS"
                np.savetxt(logs / f"dust_micron_density_profile_{label:08d}.dat",
                           [[label, 1, 1e-8], [label, 2, 2e-8]])
                np.savetxt(logs / f"micron_dust_size_evolution_{label:08d}.dat",
                           [[label, 1, 1.5, 1e-4]])
            figure = plot_output(run, components=["micron", "dust"])
            try:
                figure.canvas.draw()
                for index, axis in enumerate(figure.axes[:4]):
                    spec = axis.get_subplotspec()
                    self.assertEqual(spec.get_gridspec().get_geometry(), (2, 2))
                    self.assertEqual(spec.rowspan.start, index // 2)
                    self.assertEqual(spec.colspan.start, index % 2)
                self.assertEqual(figure.axes[0].get_title(), "Dust surface density")
                self.assertIn("Small-grain", figure.axes[1].get_title())
                self.assertEqual(figure.axes[2].get_title(), "Dust positions (colour = ID)")
                self.assertIn("Small-grain", figure.axes[3].get_title())
                self.assertEqual(figure.axes[2].get_xlabel(), "t [yr]")
                self.assertEqual(figure.axes[3].get_ylabel(), "r [AU]")
            finally:
                plt.close(figure)

    def test_static_gas_ascii(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            config = run / "config"
            config.mkdir()
            initial = run / "LOGS" / "density_profile_00000000.dat"
            (config / "initial_gas_profile.dat").write_text(initial.read_text())
            for path in (run / "LOGS").glob("density_profile_*.dat"):
                path.unlink()
            for label, time in ((0, 0), (10, 10.25)):
                path = run / "LOGS" / f"dust_density_profile_{label:08d}.dat"
                np.savetxt(path, [[label, 1, 0], [label, 2, 1e-7]],
                           header="Initial dust profile" if label == 0 else f"Dust at {time} years")
            snapshots = load_snapshots(run)
            self.assertEqual([snapshot.time for snapshot in snapshots], [0, 10.25])
            np.testing.assert_array_equal(snapshots[0].density["gas"],
                                          snapshots[1].density["gas"])

    def test_duplicate_identical_and_conflicting_times(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            logs = run / "LOGS"
            for prefix in ("density_profile", "dust_density_profile", "dust_size_evolution"):
                original = logs / (prefix + "_00000010.dat")
                duplicate = logs / (prefix + "_00000020.dat")
                duplicate.write_text(original.read_text())
            with self.assertWarnsRegex(UserWarning, "Duplicate identical"):
                self.assertEqual(len(load_snapshots(run)), 2)
            path = logs / "density_profile_00000020.dat"
            gas = np.loadtxt(path)
            gas[0, 2] *= 2
            np.savetxt(path, gas, header="Gas at 10.25 years")
            with self.assertRaisesRegex(ValueError, "same time"):
                load_snapshots(run)

    def test_zero_and_constant_density_heatmaps(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for value in (0.0, 1e-7):
                for label in (0, 10):
                    path = run / "LOGS" / f"dust_density_profile_{label:08d}.dat"
                    np.savetxt(path, [[label, 1, value], [label, 2, value]])
                if value == 0:
                    with self.assertWarnsRegex(UserWarning, "plotting it in black"):
                        figure = plot_output(run, kind="heatmap", components=["dust"],
                                             show_particles=False)
                else:
                    figure = plot_output(run, kind="heatmap", components=["dust"],
                                         show_particles=False)
                figure.canvas.draw()
                values = figure.axes[0].collections[0].get_array()
                if value == 0:
                    self.assertTrue(np.all(np.ma.getmaskarray(values)))
                else:
                    self.assertFalse(np.any(np.ma.getmaskarray(values)))
                    np.testing.assert_array_equal(values, value)
                plt.close(figure)

    def test_heatmap_axes_and_black_empty_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            for linear in (False, True):
                figure = plot_output(run, kind="heatmap", components=["dust"], linear=linear)
                try:
                    figure.canvas.draw()
                    axis = figure.axes[0]
                    self.assertEqual(axis.get_xlabel(), "t [yr]")
                    self.assertEqual(axis.get_ylabel(), "r [AU]")
                    self.assertEqual(axis.get_xscale(), "linear")
                    self.assertEqual(axis.get_yscale(), "linear" if linear else "log")
                    np.testing.assert_allclose(axis.get_xlim(), [0, 10.25])
                    np.testing.assert_allclose(axis.get_ylim(), [1, 2])
                    image = axis.collections[0]
                    np.testing.assert_array_equal(np.ma.getmaskarray(image.get_array()),
                                                  [[True, True], [False, False]])
                    np.testing.assert_allclose(image.get_facecolors()[:2],
                                               [[0, 0, 0, 1], [0, 0, 0, 1]])
                    np.testing.assert_allclose(axis.get_facecolor(), [0, 0, 0, 1])
                finally:
                    plt.close(figure)

    def test_initial_normalization_and_absolute_default(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            self.write_snapshots(run)
            path = run / "LOGS" / "density_profile_00000010.dat"
            gas = np.loadtxt(path)
            gas[:, 2] *= [2, .5]
            np.savetxt(path, gas, header="Gas at 10.25 years")
            for kind in ("profiles", "heatmap"):
                for linear in (False, True):
                    figure = plot_output(run, kind=kind, components=["gas"],
                                         linear=linear, normalize_initial=True)
                    try:
                        if kind == "heatmap":
                            np.testing.assert_allclose(figure.axes[0].collections[0].get_array(),
                                                       [[1, 2], [1, .5]])
                            self.assertIn("Sigma(r,0)", figure.axes[1].get_ylabel())
                        else:
                            np.testing.assert_allclose(figure.axes[0].lines[0].get_ydata(), [1, 1])
                            np.testing.assert_allclose(figure.axes[0].lines[1].get_ydata(), [2, .5])
                        figure.canvas.draw()
                    finally:
                        plt.close(figure)
            figure = plot_output(run, kind="heatmap", components=["gas"])
            np.testing.assert_allclose(figure.axes[0].collections[0].get_array(),
                                       [[1e-5, 2e-5], [2e-5, 1e-5]])
            plt.close(figure)
            with self.assertWarnsRegex(UserWarning, "undefined ratios"):
                figure = plot_output(run, kind="heatmap", components=["dust"],
                                     normalize_initial=True)
            np.testing.assert_array_equal(
                np.ma.getmaskarray(figure.axes[0].collections[0].get_array()),
                [[True, True], [False, False]])
            plt.close(figure)
            path = run / "LOGS" / "density_profile_00000000.dat"
            path.write_text(path.read_text().replace("Gas at 0 years", "Gas at 1 years"))
            with self.assertRaisesRegex(ValueError, "t=0"):
                plot_output(run, normalize_initial=True)


if __name__ == "__main__":
    unittest.main()
