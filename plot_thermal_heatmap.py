"""Plot nine gas/thermal time-radius heatmaps from a radiative ASCII run.

Usage: python plot_thermal_heatmap.py output_thermal [-o heatmaps.png]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import numpy as np

from plot_gas import load_snapshot, normalized_profile


def load_profiles(run):
    logs = Path(run) / "LOGS"
    gas_files = sorted(logs.glob("density_profile_*.dat"))
    if len(gas_files) < 2:
        raise ValueError(f"{logs}: need at least two gas snapshots.")
    times, gas_profiles, thermal_profiles = [], [], []
    time_unit_s = 365.242199 * 86400.0 / (2.0 * np.pi)
    pressure_unit_cgs = 1.989e33 / (1.496e13 * time_unit_s**2)
    initial = None
    for gas_file in gas_files:
        thermal_file = logs / gas_file.name.replace("density_profile_", "temperature_profile_", 1)
        time, gas = load_snapshot(gas_file)
        thermal_time, thermal = load_snapshot(thermal_file)
        if gas.ndim != 2 or gas.shape[1] != 5 or len(gas) < 2:
            raise ValueError(f"{gas_file}: expected at least two rows and five columns.")
        if thermal.shape != (len(gas), 9):
            raise ValueError(f"{thermal_file}: expected matching rows and nine columns.")
        if not np.isclose(time, thermal_time, rtol=1e-5, atol=1e-8):
            raise ValueError(f"{gas_file}: gas and thermal snapshot times disagree.")
        if (not np.all(np.isfinite(gas)) or not np.all(np.isfinite(thermal)) or
                np.any(gas[:, 1:4] <= 0) or np.any(np.diff(gas[:, 1]) <= 0) or
                np.any(thermal[:, [1, 2, 3, 5, 6]] <= 0) or
                np.any(thermal[:, [4, 7]] < 0)):
            raise ValueError(f"{gas_file}: invalid radii, density or thermal quantities.")
        if not np.allclose(gas[:, 1], thermal[:, 0], rtol=1e-6, atol=0):
            raise ValueError(f"{thermal_file}: radial grid differs from gas snapshot.")
        if initial is None:
            initial = gas
        gas_profiles.append((
            normalized_profile(gas, initial, 2, str(run)), gas[:, 4],
            np.sqrt(gas[:, 3] * pressure_unit_cgs / thermal[:, 3]) / 1e5,
        ))
        thermal_profiles.append(thermal)
        times.append(thermal_time)
    times = np.asarray(times)
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError(f"{logs}: snapshot times must be finite and strictly increasing.")
    gas_data = np.asarray(gas_profiles)
    thermal_data = np.asarray(thermal_profiles)
    fields = [gas_data[:, 0].T, gas_data[:, 1].T]
    fields.extend(thermal_data[:, :, column].T for column in (1, 2, 5, 6, 7, 8))
    fields[2] = 100.0 * (fields[2] / fields[2][:, :1] - 1.0)
    fields.append(gas_data[:, 2].T)
    fields[8] = 100.0 * (fields[8] / fields[8][:, :1] - 1.0)
    return times, initial[:, 1], fields


def positive_norm(values):
    minimum, maximum = float(values.min()), float(values.max())
    if minimum == maximum:
        if minimum == 0:
            return Normalize(vmin=0, vmax=1)
        minimum /= 1.01
        maximum *= 1.01
    if minimum == 0:
        return Normalize(vmin=minimum, vmax=maximum)
    return LogNorm(vmin=minimum, vmax=maximum)


def plot_heatmaps(run, linthresh=1e-12):
    if not np.isfinite(linthresh) or linthresh <= 0:
        raise ValueError("linthresh must be finite and positive.")
    times, radii, fields = load_profiles(run)
    titles = (
        "Normalized gas surface density", "Pressure gradient",
        "Midplane temperature change", "Aspect ratio", "Opacity", "Optical depth",
        "Viscous heating", "Net radiative cooling", "Sound speed change",
    )
    labels = (
        r"$\Sigma_g(r,t)/\Sigma_g(r,t_0)$", r"$dP/dr$ (code units; symlog)",
        r"$100\,[T(r,t)/T(r,t_0)-1]$ [%]", "H/r", r"$\kappa$ [cm$^2$ g$^{-1}$]", r"$\tau$",
        r"$Q^+$ [erg cm$^{-2}$ s$^{-1}$]",
        r"$Q^-$ [erg cm$^{-2}$ s$^{-1}$]",
        r"$100\,[c_s(r,t)/c_s(r,t_0)-1]$ [%]",
    )
    flux_values = np.stack(fields[6:8])
    signed_flux = np.any(flux_values < 0)
    if signed_flux:
        limit = float(np.abs(flux_values).max())
        flux_norm = SymLogNorm(linthresh=limit * 1e-3, vmin=-limit, vmax=limit, base=10)
    else:
        flux_norm = positive_norm(flux_values)
    fig, axes = plt.subplots(3, 3, figsize=(20, 11), sharex=True, sharey=True,
                             constrained_layout=True)
    for index, (axis, values, title, label) in enumerate(
            zip(axes.flat, fields, titles, labels)):
        norm = positive_norm(values) if index not in (1, 2, 6, 7, 8) else None
        cmap = "viridis"
        if index == 1:
            limit = max(linthresh, float(np.abs(values).max()))
            norm = SymLogNorm(linthresh=linthresh, vmin=-limit, vmax=limit, base=10)
            cmap = "RdBu_r"
        elif index in (2, 8):
            norm = Normalize(vmin=-0.1, vmax=0.1, clip=True)
            cmap = "RdBu_r"
        elif index in (6, 7):
            norm = flux_norm
            cmap = "RdBu_r" if signed_flux else "magma"
        elif index >= 2:
            cmap = "plasma"
        mesh = axis.pcolormesh(times, radii, values, norm=norm, cmap=cmap,
                               shading="nearest")
        if index == 1 and values.min() < 0 < values.max():
            axis.contour(times, radii, values, levels=[0], colors="white",
                         linewidths=0.8)
            radial_indices, time_indices = np.nonzero(
                (values[:-1] > 0) & (values[1:] < 0))
            axis.plot(times[time_indices], radii[radial_indices + 1], ".",
                      color="black", markersize=4, zorder=4)
        axis.set(title=title, yscale="log", ylim=(radii.min(), radii.max()))
        fig.colorbar(mesh, ax=axis, label=label,
                     extend="both" if index in (2, 8) else "neither")
    for axis in axes[-1]:
        axis.set_xlabel("t [yr]")
    for axis in axes[:, 0]:
        axis.set_ylabel("r [AU]")
    fig.suptitle(f"{Path(run).name} - gas and radiative thermodynamics\n"
                 "Pressure gradient: white zero contours; black dots: pressure maxima")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="radiative run directory containing ASCII snapshots in LOGS")
    parser.add_argument("--linthresh", type=float, default=1e-12)
    parser.add_argument("-o", "--output", help="save instead of opening a window")
    args = parser.parse_args()
    try:
        fig = plot_heatmaps(args.run, args.linthresh)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    if args.output:
        fig.savefig(args.output, dpi=180)
        plt.close(fig)
        print(f"Saved {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
