"""Quick ASCII overview: gas/dust surface density and ID-coloured particle positions.

Usage: python plot_dust.py RUN_DIR [-o overview.png] [--linear]
Density curves are coloured by time; particle colours identify individual grains.
"""

import argparse
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np

from plot_gas import load_snapshot


def read_table(path, columns, allow_empty=False):
    with path.open() as stream:
        has_data = any(line.strip() and not line.lstrip().startswith("#") for line in stream)
    if not has_data:
        if allow_empty:
            return np.empty((0, columns))
        raise ValueError(f"{path}: empty density snapshot.")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != columns or not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: expected {columns} finite numeric columns.")
    return data


def plot_dust(run, linear=False):
    run = Path(run)
    logs = run / "LOGS"
    paths = sorted(logs.glob("density_profile_*.dat"))
    if not paths:
        raise ValueError(f"{logs}: no ASCII density_profile_*.dat snapshots found.")
    snapshots = []
    for path in paths:
        suffix = path.name.removeprefix("density_profile_")
        time, _ = load_snapshot(path)
        gas = read_table(path, 5)
        dust = read_table(logs / f"dust_density_profile_{suffix}", 3)
        particles = read_table(logs / f"dust_size_evolution_{suffix}", 4, allow_empty=True)
        if (not np.isfinite(time) or time < 0 or np.any(gas[:, 1] <= 0) or
                np.any(np.diff(gas[:, 1]) <= 0) or np.any(gas[:, 2] < 0)):
            raise ValueError(f"{path}: invalid time, gas radius or density.")
        if (len(dust) != len(gas) or not np.allclose(
                dust[:, 1], gas[:, 1], rtol=1e-5, atol=1e-5)):
            raise ValueError(
                f"{suffix}: dust and gas grids differ. Regenerate snapshots with the "
                "corrected dust-density writer; old files may have shifted/truncated rows."
            )
        if np.any(dust[:, 2] < 0):
            raise ValueError(f"{suffix}: negative dust density.")
        ids = particles[:, 1]
        if (np.any(ids < 0) or np.any(ids != np.floor(ids)) or
                len(np.unique(ids)) != len(ids)):
            raise ValueError(f"{suffix}: particle IDs must be unique nonnegative integers.")
        # Particle headers historically contain t/(2*pi); use the matching gas time.
        particles = particles[(particles[:, 2] >= gas[0, 1]) &
                              (particles[:, 2] <= gas[-1, 1])]
        snapshots.append((time, gas, dust, particles))
    snapshots.sort(key=lambda snapshot: snapshot[0])
    times = np.array([snapshot[0] for snapshot in snapshots])
    if len(times) > 1 and np.any(np.diff(times) <= 0):
        raise ValueError(f"{run}: snapshot times must be strictly increasing.")
    norm = Normalize(times[0], times[-1] if times[-1] > times[0] else times[0] + 1)
    cmap = plt.get_cmap("viridis")
    figure, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True,
                                constrained_layout=True)
    max_id = max((np.max(p[:, 1]) for _, _, _, p in snapshots if len(p)), default=1)
    id_norm = Normalize(0, max(1, max_id))
    for time, gas, dust, particles in snapshots:
        for axis, data in zip(axes[:2], (gas, dust)):
            values = np.ma.masked_less_equal(data[:, 2], 0) if not linear else data[:, 2]
            axis.plot(data[:, 1], values, color=cmap(norm(time)), linewidth=1, alpha=.8)
        if len(particles):
            axes[2].scatter(particles[:, 2], np.full(len(particles), time),
                            c=particles[:, 1], cmap="turbo", norm=id_norm,
                            s=5, linewidths=0, rasterized=True)
    for axis, column, title, label in (
        (axes[0], 1, "Gas surface density", r"$\Sigma_g$ [M$_\odot$/AU$^2$]"),
        (axes[1], 2, "Dust surface density", r"$\Sigma_d$ [M$_\odot$/AU$^2$]"),
    ):
        for index, style, name in ((0, "--", "Initial"), (-1, "-", "Final")):
            data = snapshots[index][column]
            values = np.ma.masked_less_equal(data[:, 2], 0) if not linear else data[:, 2]
            axis.plot(data[:, 1], values, style, color="black",
                      linewidth=1.5, label=f"{name}: {times[index]:g} yr")
        axis.set(title=title, ylabel=label, yscale="linear" if linear else "log")
        axis.legend(fontsize=8)
    if not any(len(p) for _, _, _, p in snapshots):
        warnings.warn("No active particles in the saved snapshots.", stacklevel=2)
        axes[2].text(.5, .5, "No active particles", transform=axes[2].transAxes, ha="center")
    axes[2].set(title="Particle positions (colour = ID)", xlabel="r [AU]", ylabel="t [yr]")
    for axis in axes:
        axis.set_xscale("linear" if linear else "log")
        axis.grid(alpha=.2, which="both")
    figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes[:2], label="t [yr]")
    figure.colorbar(ScalarMappable(norm=id_norm, cmap="turbo"), ax=axes[2], label="Particle ID")
    figure.suptitle(run.name)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="run directory containing ASCII snapshots in LOGS")
    parser.add_argument("-o", "--output", help="save instead of opening a window")
    parser.add_argument("--linear", action="store_true", help="use linear radius and density axes")
    args = parser.parse_args()
    try:
        figure = plot_dust(args.run, linear=args.linear)
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
