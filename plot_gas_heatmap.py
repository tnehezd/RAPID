"""Plot time-radius heatmaps of normalized gas density and pressure gradient.

Usage:
    python plot_gas_heatmap.py RUN_DIR [RUN_DIR ...] [--ref REF_DIR] [-o out.png]

Each directory must contain LOGS/density_profile_*.dat in RAPID's five-column
format. Density is divided by each run's earliest saved profile. Columns show
separate runs, with shared color scales. The pressure gradient is read directly
from the snapshots. White solid contours show dP/dr=0; black dots mark radial
positive-to-negative gradient crossings (pressure maxima), at the outer sample
of each crossing pair.
"""
import argparse
import os

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np

from plot_gas import normalized_profile, snapshots


def load_heatmap(run):
    snaps = sorted(snapshots(run), key=lambda snap: snap[0])
    times = np.array([t for t, _ in snaps])
    initial = snaps[0][1]
    if initial.ndim != 2 or initial.shape[1] != 5 or len(initial) < 2:
        raise SystemExit(f"{run}: expected at least two rows and five columns.")
    radii = initial[:, 1]
    if (len(times) < 2 or not np.all(np.isfinite(times)) or
            np.any(np.diff(times) <= 0)):
        raise SystemExit(f"{run}: need at least two distinct, increasing snapshot times.")
    if (not np.all(np.isfinite(radii)) or np.any(radii <= 0) or
            np.any(np.diff(radii) <= 0)):
        raise SystemExit(f"{run}: radii must be finite, positive and increasing.")
    density = np.array([normalized_profile(data, initial, 2, run) for _, data in snaps]).T
    gradient = np.array([data[:, 4] for _, data in snaps]).T
    if (not np.all(np.isfinite(density)) or np.any(density <= 0) or
            not np.all(np.isfinite(gradient))):
        raise SystemExit(f"{run}: need positive density ratios and finite pressure gradients.")
    return times, radii, density, gradient


def plot_heatmaps(runs, linthresh=1e-12):
    if not np.isfinite(linthresh) or linthresh <= 0:
        raise SystemExit("linthresh must be finite and positive.")
    datasets = [load_heatmap(run) for run in runs]
    density_min = min(d.min() for _, _, d, _ in datasets)
    density_max = max(d.max() for _, _, d, _ in datasets)
    if density_min == density_max:
        density_min /= 1.01
        density_max *= 1.01
    density_norm = LogNorm(vmin=density_min, vmax=density_max)
    gradient_limit = max(linthresh, max(np.abs(g).max() for _, _, _, g in datasets))
    gradient_norm = SymLogNorm(linthresh=linthresh, vmin=-gradient_limit,
                               vmax=gradient_limit, base=10)
    fig, axes = plt.subplots(2, len(runs), figsize=(6 * len(runs), 9),
                             squeeze=False, sharex=True, sharey=True,
                             constrained_layout=True)
    for j, (run, (times, radii, density, gradient)) in enumerate(zip(runs, datasets)):
        density_mesh = axes[0, j].pcolormesh(
            times, radii, density, norm=density_norm, cmap="viridis", shading="nearest")
        gradient_mesh = axes[1, j].pcolormesh(
            times, radii, gradient, norm=gradient_norm, cmap="RdBu_r", shading="nearest")
        if gradient.min() < 0 < gradient.max():
            axes[1, j].contour(times, radii, gradient, levels=[0], colors="white",
                               linestyles="-", linewidths=0.8)
            radial_indices, time_indices = np.nonzero(
                (gradient[:-1, :] > 0) & (gradient[1:, :] < 0))
            axes[1, j].plot(times[time_indices], radii[radial_indices + 1],
                            ".", color="black", markersize=4, zorder=4)
        else:
            axes[1, j].text(0.03, 0.97, "No pressure-gradient sign change",
                            transform=axes[1, j].transAxes, va="top", fontsize=8,
                            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
        name = os.path.basename(os.path.normpath(run))
        axes[0, j].set_title(name + "\nNormalized gas surface density")
        axes[1, j].set_title("Pressure gradient\nWhite: $dP/dr=0$; black dots: maxima")
        axes[1, j].set_xlabel("t [yr]")
        for ax in axes[:, j]:
            ax.set_yscale("log")
            ax.set_ylim(min(r.min() for _, r, _, _ in datasets),
                        max(r.max() for _, r, _, _ in datasets))
    axes[0, 0].set_ylabel("$r$ [AU]")
    axes[1, 0].set_ylabel("$r$ [AU]")
    fig.colorbar(density_mesh, ax=axes[0].tolist(), label=r"$\Sigma_g(r,t)/\Sigma_g(r,t_0)$")
    fig.colorbar(gradient_mesh, ax=axes[1].tolist(), label="$dP/dr$ (code units; symlog)")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+", help="run directories containing LOGS")
    parser.add_argument("--ref", help="reference run, displayed as a separate column")
    parser.add_argument("--linthresh", type=float, default=1e-12,
                        help="linear threshold for the pressure-gradient symlog scale (default: 1e-12)")
    parser.add_argument("-o", "--output", help="save figure instead of opening a window")
    args = parser.parse_args()
    runs = args.runs + ([args.ref] if args.ref else [])
    fig = plot_heatmaps(runs, args.linthresh)
    if args.output:
        fig.savefig(args.output, dpi=300)
        print(f"Saved {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
