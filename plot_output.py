"""Standalone RAPID output plotter (ASCII or HDF5, gas and dust populations).

python3 plot_output.py RUN_DIR [-o overview.png]
python3 plot_output.py RUN_DIR --kind heatmap --components gas dust micron
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import warnings

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
import numpy as np


COMPONENTS = ("gas", "dust", "micron")
LABELS = {"gas": "Gas", "dust": "Dust", "micron": "Small-grain dust"}
GROUPS = {"gas": "gas_grid", "dust": "dust_grid", "micron": "micron_grid"}
ASCII_DENSITY = {"dust": "dust_density_profile", "micron": "dust_micron_density_profile"}
ASCII_PARTICLES = {"dust": "dust_size_evolution", "micron": "micron_dust_size_evolution"}


@dataclass
class Snapshot:
    time: float
    radius: np.ndarray
    density: dict
    particles: dict


def read_snapshot_time(path):
    with path.open() as stream:
        for line in stream:
            if not line.lstrip().startswith("#"):
                break
            match = re.search(r"at\s+([0-9.eE+-]+)\s+years", line)
            if match:
                return float(match.group(1))
            if line.strip() in ("# Initial gas profile", "# Initial dust profile",
                                "# Initial micron-sized dust profile"):
                return 0.0
    raise ValueError(f"{path}: missing snapshot time in years.")


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


def load_ascii(path, populations, static_gas=None):
    gas = read_table(static_gas or path, 5)
    prefix = "dust_density_profile" if static_gas else "density_profile"
    suffix = path.name[len(prefix):]
    density = {"gas": gas[:, 2]}
    particles = {}
    for component in populations:
        dust_path = path.parent / (ASCII_DENSITY[component] + suffix)
        dust = read_table(dust_path, 3)
        if len(dust) != len(gas) or not np.allclose(
                dust[:, 1], gas[:, 1], rtol=1e-5, atol=1e-5):
            raise ValueError(
                f"{dust_path}: dust and gas grids differ. Regenerate old "
                "shifted/truncated output with the corrected writer."
            )
        density[component] = dust[:, 2]
        # Use gas time: older particle headers contain t/(2*pi).
        particles[component] = read_table(
            path.parent / (ASCII_PARTICLES[component] + suffix), 4, allow_empty=True
        )[:, 1:]
    return Snapshot(read_snapshot_time(path), gas[:, 1], density, particles)


def load_hdf5(path):
    try:
        import h5py
    except ImportError as exc:
        raise ValueError("HDF5 input requires h5py: python3 -m pip install h5py") from exc
    with h5py.File(path, "r") as source:
        required = ("gas_grid/radial_grid", "gas_grid/surface_density", "frame/time")
        for key in required:
            if key not in source:
                raise ValueError(f"{path}: missing HDF5 dataset {key}.")
        time = np.asarray(source["frame/time"]).reshape(-1)
        if time.size != 1:
            raise ValueError(f"{path}: frame/time must contain one value.")
        radius = np.asarray(source["gas_grid/radial_grid"])
        density = {}
        for component, group in GROUPS.items():
            key = group + "/surface_density"
            if group in source:
                if key not in source:
                    raise ValueError(f"{path}: missing HDF5 dataset {key}.")
                density[component] = np.asarray(source[key])
        particles = {}
        for component, group in (("dust", "particles"), ("micron", "micron_particles")):
            if group not in source:
                continue
            keys = [group + "/" + name for name in ("index", "position", "size")]
            if any(key not in source for key in keys):
                raise ValueError(f"{path}: incomplete HDF5 particle group {group}.")
            arrays = [np.asarray(source[key]) for key in keys]
            if any(array.ndim != 1 or array.shape != arrays[0].shape for array in arrays):
                raise ValueError(f"{path}: inconsistent particle arrays in {group}.")
            particles[component] = np.column_stack(arrays)
    return Snapshot(float(time[0]), radius, density, particles)


def validate_snapshot(snapshot, path):
    radius = snapshot.radius
    if (not np.isfinite(snapshot.time) or snapshot.time < 0 or radius.ndim != 1
            or len(radius) < 2 or not np.all(np.isfinite(radius))
            or np.any(radius <= 0) or np.any(np.diff(radius) <= 0)):
        raise ValueError(f"{path}: invalid time or radial grid.")
    for component, density in snapshot.density.items():
        if (density.shape != radius.shape or not np.all(np.isfinite(density))
                or np.any(density < 0)):
            raise ValueError(f"{path}: invalid {component} density or grid shape.")
    for component, particles in snapshot.particles.items():
        if not np.all(np.isfinite(particles)):
            raise ValueError(f"{path}: nonfinite {component} particle data.")
        ids = particles[:, 0]
        if (np.any(ids < 0) or np.any(ids != np.floor(ids))
                or len(np.unique(ids)) != len(ids) or np.any(particles[:, 2] < 0)):
            raise ValueError(f"{path}: invalid {component} particle IDs or sizes.")
        active = ((particles[:, 1] >= radius[0]) & (particles[:, 1] <= radius[-1])
                  & (particles[:, 2] > 0))
        snapshot.particles[component] = particles[active]


def load_snapshots(source, input_format="auto"):
    source = Path(source)
    static_gas = None
    if input_format not in ("auto", "ascii", "hdf5"):
        raise ValueError("format must be auto, ascii or hdf5.")
    if source.is_file():
        candidates = [source]
        detected = "hdf5" if source.suffix in (".h5", ".hdf5") else "ascii"
        if input_format != "auto" and input_format != detected:
            raise ValueError("File extension disagrees with --format.")
        input_format = detected
    else:
        logs = source / "LOGS" if (source / "LOGS").is_dir() else source
        ascii_paths = sorted(logs.glob("density_profile_*.dat"))
        if not ascii_paths and input_format != "hdf5":
            dust_paths = sorted(logs.glob("dust_density_profile_*.dat"))
            if dust_paths:
                static_gas = logs.parent / "config" / "initial_gas_profile.dat"
                if not static_gas.is_file():
                    raise ValueError("No gas snapshots or saved static initial gas profile found.")
                ascii_paths = dust_paths
        hdf5_paths = sorted(logs.glob("snapshot_*.h5")) + sorted(logs.glob("snapshot_*.hdf5"))
        if input_format == "auto":
            if ascii_paths and hdf5_paths:
                raise ValueError("Both ASCII and HDF5 snapshots found; select --format.")
            input_format = "hdf5" if hdf5_paths else "ascii"
        candidates = hdf5_paths if input_format == "hdf5" else ascii_paths
    if not candidates:
        raise ValueError(f"{source}: no {input_format} snapshots found.")
    populations = []
    if input_format == "ascii":
        for component in ("dust", "micron"):
            directory = candidates[0].parent
            has_density = any(directory.glob(ASCII_DENSITY[component] + "_*.dat"))
            has_particles = any(directory.glob(ASCII_PARTICLES[component] + "_*.dat"))
            if has_density != has_particles:
                raise ValueError(f"{directory}: incomplete {component} output.")
            if has_density:
                populations.append(component)
    snapshots = []
    for path in candidates:
        snapshot = (load_hdf5(path) if input_format == "hdf5"
                    else load_ascii(path, populations, static_gas))
        validate_snapshot(snapshot, path)
        snapshots.append(snapshot)
    snapshots.sort(key=lambda snapshot: snapshot.time)
    unique = []
    for snapshot in snapshots:
        if unique and snapshot.time == unique[-1].time:
            previous = unique[-1]
            if (snapshot.density.keys() != previous.density.keys()
                    or snapshot.particles.keys() != previous.particles.keys()
                    or not np.array_equal(snapshot.radius, previous.radius)
                    or any(not np.array_equal(snapshot.density[key], previous.density[key])
                           for key in snapshot.density)
                    or any(not np.array_equal(snapshot.particles[key], previous.particles[key])
                           for key in snapshot.particles)):
                raise ValueError(f"{source}: different snapshots have the same time.")
            warnings.warn(f"Duplicate identical snapshot at t={snapshot.time:g} yr; "
                          "plotting it once.", stacklevel=2)
        else:
            unique.append(snapshot)
    snapshots = unique
    first = snapshots[0]
    for previous, current in zip(snapshots, snapshots[1:]):
        if current.time <= previous.time:
            raise ValueError(f"{source}: snapshot times must be strictly increasing.")
        if current.density.keys() != first.density.keys():
            raise ValueError(f"{source}: populations change between snapshots.")
        if current.particles.keys() != first.particles.keys():
            raise ValueError(f"{source}: particle datasets change between snapshots.")
        if (current.radius.shape != first.radius.shape
                or not np.allclose(current.radius, first.radius, rtol=1e-5, atol=1e-5)):
            raise ValueError(f"{source}: radial grids change between snapshots.")
    return snapshots


def plot_output(source, kind="profiles", components=None, input_format="auto",
                linear=False, cmap="viridis", stride=1, show_particles=True,
                normalize_initial=False):
    if kind not in ("profiles", "heatmap") or stride < 1:
        raise ValueError("kind must be profiles or heatmap; stride must be positive.")
    colourmap = plt.get_cmap(cmap).copy()
    if kind == "heatmap":
        colourmap.set_bad("black")
    snapshots = load_snapshots(source, input_format)
    if normalize_initial and snapshots[0].time != 0:
        raise ValueError("Initial normalization requires the t=0 snapshot.")
    available = snapshots[0].density
    components = list(available) if components is None else list(dict.fromkeys(components))
    if not components or any(component not in available for component in components):
        raise ValueError(f"Requested components unavailable; choose from {', '.join(available)}.")
    if kind == "heatmap" and len(snapshots) < 2:
        raise ValueError("A heatmap requires at least two snapshots.")
    selected = snapshots[::stride]
    if selected[-1] is not snapshots[-1]:
        selected.append(snapshots[-1])
    particle_components = []
    if show_particles and kind == "profiles":
        for component in components:
            if component == "gas":
                continue
            if component in snapshots[0].particles:
                particle_components.append(component)
            else:
                warnings.warn(f"{LABELS[component]} particle positions are not saved in this "
                              "output; only its density is plotted.", stacklevel=2)
    rows = len(components) + len(particle_components)
    dust_grid = (kind == "profiles" and set(components) == {"dust", "micron"}
                 and set(particle_components) == {"dust", "micron"})
    if dust_grid:
        components = ["dust", "micron"]
        particle_components = ["dust", "micron"]
    figsize = (10, 7) if dust_grid else (
        (7, 4.5 * rows) if kind == "heatmap" else (9, 3 * rows))
    figure, axes = plt.subplots(2 if dust_grid else rows, 2 if dust_grid else 1,
                                figsize=figsize, squeeze=False,
                                sharex=kind == "heatmap", constrained_layout=True)
    axes = axes.ravel()
    times = np.array([snapshot.time for snapshot in selected])
    time_norm = Normalize(times[0], times[-1] if times[-1] > times[0] else times[0] + 1)
    for axis, component in zip(axes, components):
        densities = np.array([snapshot.density[component] for snapshot in selected])
        density_label = r"$\Sigma$ [M$_\odot$/AU$^2$]"
        if normalize_initial:
            initial = snapshots[0].density[component]
            if np.any(initial == 0):
                warnings.warn(f"{LABELS[component]} initial density is zero at some radii; "
                              "undefined ratios are masked.", stacklevel=2)
            with np.errstate(over="raise", invalid="raise"):
                try:
                    ratios = np.divide(densities, initial, out=np.zeros_like(densities),
                                       where=initial > 0)
                except FloatingPointError as exc:
                    raise ValueError(f"{LABELS[component]} initial density ratios overflow.") from exc
            densities = np.ma.array(ratios, mask=np.broadcast_to(initial == 0, ratios.shape))
            density_label = r"$\Sigma(r,t)/\Sigma(r,0)$"
        values = densities if linear else np.ma.masked_less_equal(densities, 0)
        if kind == "heatmap":
            values = np.ma.masked_less_equal(densities, 0)
            positive = values.compressed()
            if not positive.size:
                warnings.warn(f"{LABELS[component]} has no positive density values; "
                              "plotting it in black.", stacklevel=2)
            norm = LogNorm(positive.min(), positive.max()) if not linear and positive.size else None
            image = axis.pcolormesh(times, snapshots[0].radius, values.T, shading="nearest",
                                    cmap=colourmap, norm=norm, rasterized=True)
            figure.colorbar(image, ax=axis, label=density_label)
            axis.set_ylabel("r [AU]")
            axis.set_facecolor("black")
            axis.set_xlim(times[0], times[-1])
            axis.set_ylim(snapshots[0].radius[0], snapshots[0].radius[-1])
        else:
            for snapshot, density in zip(selected, values):
                axis.plot(snapshot.radius, density, color=colourmap(time_norm(snapshot.time)),
                          linewidth=1, alpha=.8)
            for index, style, label in ((0, "--", "Initial"), (-1, "-", "Final")):
                axis.plot(selected[index].radius, values[index], style, color="black",
                          label=f"{label}: {times[index]:g} yr", linewidth=1.5)
            axis.set_ylabel(density_label)
            axis.set_xlabel("r [AU]")
            axis.set_yscale("linear" if linear else "log")
            axis.legend(fontsize=9)
        axis.set_title(LABELS[component] + " surface density"
                       + (" / initial" if normalize_initial else ""))
    if kind == "profiles":
        figure.colorbar(ScalarMappable(norm=time_norm, cmap=colourmap),
                       ax=list(axes[:len(components)]), label="t [yr]")
    for axis, component in zip(axes[len(components):], particle_components):
        particles = [snapshot.particles[component] for snapshot in selected]
        max_id = max((np.max(points[:, 0]) for points in particles if len(points)), default=1)
        id_norm = Normalize(0, max(1, max_id))
        for snapshot, points in zip(selected, particles):
            if len(points):
                axis.scatter(np.full(len(points), snapshot.time), points[:, 1],
                             c=points[:, 0], cmap="turbo", norm=id_norm,
                             s=5, linewidths=0, rasterized=True)
        if not any(len(points) for points in particles):
            warnings.warn(f"No active particles for {component}.", stacklevel=2)
            axis.text(.5, .5, "No active particles", transform=axis.transAxes, ha="center")
        axis.set(title=LABELS[component] + " positions (colour = ID)",
                 xlabel="t [yr]", ylabel="r [AU]")
        axis.set_yscale("linear" if linear else "log")
        axis.set_ylim(snapshots[0].radius[0], snapshots[0].radius[-1])
        if times[-1] > times[0]:
            axis.set_xlim(times[0], times[-1])
        figure.colorbar(ScalarMappable(norm=id_norm, cmap="turbo"), ax=axis, label="Particle ID")
    for index, axis in enumerate(axes):
        if kind == "heatmap":
            axis.set_yscale("linear" if linear else "log")
            axis.grid(False)
        else:
            axis.set_xscale("log" if not linear and index < len(components) else "linear")
            axis.grid(alpha=.2, which="both")
    if kind == "heatmap":
        axes[-1].set_xlabel("t [yr]")
    figure.suptitle(Path(source).name)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="run directory, LOGS directory or one snapshot")
    parser.add_argument("-o", "--output", help="save instead of opening a window")
    parser.add_argument("--kind", choices=("profiles", "heatmap"), default="profiles")
    parser.add_argument("--format", choices=("auto", "ascii", "hdf5"), default="auto")
    parser.add_argument("--components", nargs="+", choices=COMPONENTS)
    parser.add_argument("--linear", action="store_true", help="linear radius and density scales")
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--stride", type=int, default=1, help="plot every Nth snapshot, plus final")
    parser.add_argument("--no-particles", action="store_true", help="only surface densities")
    parser.add_argument("--normalize-initial", action="store_true",
                        help="divide each density by its local t=0 value")
    args = parser.parse_args()
    try:
        figure = plot_output(args.run, kind=args.kind, components=args.components,
                             input_format=args.format, linear=args.linear, cmap=args.cmap,
                             stride=args.stride, show_particles=not args.no_particles,
                             normalize_initial=args.normalize_initial)
        if args.output:
            figure.savefig(args.output, dpi=180)
            plt.close(figure)
            print(f"Saved {args.output}")
        else:
            plt.show()
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
