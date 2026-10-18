"""Plot a nine-column equilibrium or time-dependent thermal profile."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_temperature(path):
    with open(path) as profile:
        energy_mode = "Radiative energy evolution" in profile.readline()
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] != 9 or len(data) < 2:
        raise ValueError(f"{path}: expected at least two rows and nine columns.")
    if (not np.all(np.isfinite(data)) or
            np.any(data[:, [0, 1, 2, 3, 5, 6]] <= 0) or
            np.any(data[:, [4, 7]] < 0)):
        raise ValueError(f"{path}: expected positive radius, T, H/r, density, opacity and tau; "
                         "viscosity and viscous heating must be nonnegative and all data finite.")
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f"{path}: radii must be increasing.")

    radius, temperature, aspect, _, _, opacity, tau, heating, cooling = data.T
    figure, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True,
                                constrained_layout=True)
    profiles = (
        (axes[0, 0], temperature, "Midplane temperature", "T [K]"),
        (axes[0, 1], aspect, "Aspect ratio", "H/r"),
        (axes[0, 2], opacity, "Opacity", r"$\kappa$ [cm$^2$ g$^{-1}$]"),
        (axes[1, 0], tau, "Optical depth", r"$\tau$"),
    )
    for axis, values, title, label in profiles:
        axis.plot(radius, values, color="tab:purple")
        axis.set(title=title, ylabel=label, yscale="log")
    axes[1, 0].axhline(1.0, color="gray", linestyle="--", label=r"$\tau=1$")
    axes[1, 0].legend()

    axes[1, 1].plot(radius, heating, label="Viscous heating")
    axes[1, 1].plot(radius, cooling, "--", label="Net radiative cooling")
    axes[1, 1].set(title="Energy balance", ylabel=r"Flux [erg cm$^{-2}$ s$^{-1}$]",
                   yscale="log" if np.any(heating > 0) or np.any(cooling > 0) else "linear")
    if np.any(cooling < 0):
        limit = float(np.max(np.abs(np.stack((heating, cooling)))))
        axes[1, 1].set_yscale("symlog", linthresh=limit * 1e-3)
    axes[1, 1].legend()

    absolute_source = energy_mode or np.any((heating == 0) & (cooling != 0))
    if absolute_source:
        values = heating - cooling
        title = "Net thermal source (excluding transport/work)"
        label = r"$Q^+ - Q^-$ [erg cm$^{-2}$ s$^{-1}$]"
    else:
        values = np.divide(cooling - heating, heating, out=np.zeros_like(heating),
                           where=heating > 0)
        title = f"Max relative imbalance: {np.max(np.abs(values)):.2e}"
        label = r"$(Q^- - Q^+)/Q^+$"
    axes[1, 2].plot(radius, values, color="tab:blue")
    axes[1, 2].axhline(0.0, color="gray", linewidth=0.8)
    axes[1, 2].set(title=title, ylabel=label)
    axes[1, 2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for axis in axes.flat:
        axis.set_xscale("log")
        axis.grid(True, which="both", alpha=0.2)
    for axis in axes[1]:
        axis.set_xlabel("Radius [AU]")
    figure.suptitle(f"Radiative {'energy evolution' if energy_mode else 'equilibrium'} - {Path(path).name}")
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", nargs="?", default="temperature_profile.dat")
    parser.add_argument("-o", "--output", help="Save figure instead of opening a window.")
    args = parser.parse_args()
    try:
        figure = plot_temperature(args.profile)
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
