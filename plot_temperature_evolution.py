"""Plot radial temperature profiles and their percentage evolution.

Usage: python plot_temperature_evolution.py output_thermal_sg [-o temperature_evolution.png]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np

from plot_gas import load_snapshot


def plot_temperature_evolution(run):
    paths = sorted((Path(run) / "LOGS").glob("temperature_profile_*.dat"))
    if len(paths) < 2:
        raise ValueError(f"{run}: need at least two thermal snapshots.")
    times, profiles = [], []
    radius = None
    for path in paths:
        time, data = load_snapshot(path)
        if (data.ndim != 2 or data.shape[1] != 9 or len(data) < 2 or
                not np.all(np.isfinite(data[:, :2])) or
                np.any(data[:, :2] <= 0) or np.any(np.diff(data[:, 0]) <= 0)):
            raise ValueError(f"{path}: expected nine columns with positive, finite radius and T.")
        if radius is None:
            radius = data[:, 0]
        elif data.shape[0] != len(radius) or not np.allclose(
                data[:, 0], radius, rtol=1e-6, atol=0):
            raise ValueError(f"{path}: radial grid differs from the initial snapshot.")
        times.append(time)
        profiles.append(data[:, 1])
    times = np.asarray(times)
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError(f"{run}: snapshot times must be finite and strictly increasing.")
    profiles = np.asarray(profiles)
    changes = 100.0 * (profiles / profiles[0] - 1.0)
    norm = Normalize(vmin=times[0], vmax=times[-1])
    cmap = plt.get_cmap("plasma")
    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                constrained_layout=True)
    for index, time in enumerate(times):
        color = cmap(norm(time))
        for axis, values in zip(axes, (profiles, changes)):
            axis.plot(radius, values[index], color=color, linewidth=1,
                      alpha=0.8)
    for axis, values in zip(axes, (profiles, changes)):
        axis.plot(radius, values[0], color="black", linestyle="--",
                  linewidth=1.5, label=f"Initial: {times[0]:g} yr")
        axis.plot(radius, values[-1], color=cmap(1.0), linewidth=2,
                  label=f"Final: {times[-1]:g} yr")
        axis.set_xscale("log")
        axis.grid(True, which="both", alpha=0.2)
        axis.legend()
    axes[0].set(title="Midplane temperature", ylabel="T [K]", yscale="log")
    axes[1].set(title="Change from initial temperature profile",
                ylabel=r"$100\,[T(r,t)/T(r,t_0)-1]$ [%]", xlabel="r [AU]")
    figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes, label="t [yr]")
    figure.suptitle(f"{Path(run).name} - temperature evolution")
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="run directory containing thermal snapshots in LOGS")
    parser.add_argument("-o", "--output", help="save instead of opening a window")
    args = parser.parse_args()
    try:
        figure = plot_temperature_evolution(args.run)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    if args.output:
        figure.savefig(args.output, dpi=180)
        plt.close(figure)
        print(f"Saved {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
