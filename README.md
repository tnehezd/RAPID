# RAPID: Representative Approach for Particle-Integrated Disks

RAPID is a one-dimensional gas and dust evolution code for protoplanetary
disks. It combines an Eulerian gas grid with Lagrangian representative dust
particles to model viscous gas evolution, dust drift, grain growth and
fragmentation, and dust trapping near dead-zone pressure maxima.

## Documentation

The user guide, configuration parameters, numerical methods and API reference
are available on **Read the Docs**:

**https://rapiddocs.readthedocs.io/en/latest/**

Check that the documentation version matches the code used for your run.

## Development status

The current development code identifies itself as 2.3.0. Two-population dust,
self-gravity, dust
backreaction and radiative/energy thermal modes are disabled pending physical
validation; locally isothermal calculations remain available.
Two-population dust is paused pending coagulation mass-transfer validation;
single-population drift and growth remain available. Old two-population
output can still be plotted.

Recent corrections affect gas evolution, dust transport, growth and output
coordinates. Earlier results are not corrected retroactively: affected runs
must be regenerated, not merely replotted. Record the code commit and
configuration with each calculation.

## Build

Build requirements include a C compiler, GNU Make, HDF5 and an OpenMP runtime.
The YAML wrapper requires Python 3.9 or newer and PyYAML.

From the repository root:

```bash
make all
python3 -m pip install PyYAML
```

The executable is created at `bin/simulation`.

## Run

Run the supplied gas-only configuration from the repository root:

```bash
python3 run_simulation.py --config configs/config_pure_gas.yaml
```

More examples are in [configs/](configs/). For your own setup, replace the
configuration path with your YAML file. The wrapper does not create
`config.yaml` automatically.

For direct execution, `./bin/simulation --help` lists the available options.
Configuration details and parameter descriptions are in the
[online documentation](https://rapiddocs.readthedocs.io/en/latest/).

Alternatively, install the pre-built Python package:

```bash
python3 -m pip install rapidsim
rapidsim --config /path/to/config.yaml
```

The published package may contain an older RAPID binary than this development
repository.

## Plot results

For an overview of ASCII or HDF5 gas/dust snapshots:

```bash
python3 plot_output.py RUN_DIR -o overview.png
python3 plot_output.py RUN_DIR --kind heatmap -o heatmap.png
```

The standalone plotter requires NumPy and Matplotlib, plus h5py for HDF5.
It detects available gas and dust populations. Use `--components gas dust micron`
to select fields, `--no-particles` for densities only, and `--help` for other
options. Gas-only runs produce a single gas panel. Current HDF5 files store
small-grain density but not small-grain particle trajectories.

## References

If you use RAPID in your research, cite:

- Regály, Juhász & Nehéz (2017), *Interpreting Brightness Asymmetries in
  Transition Disks: Vortex at Dead Zone or Planet-carved Gap Edges?*
  [doi:10.3847/1538-4357/aa9a3f](https://doi.org/10.3847/1538-4357/aa9a3f).
- Tarczay-Nehéz (2026), *Trajectory-based dust evolution in disks: first
  results from the RAPID simulation code*.
  [doi:10.1007/s10569-026-10278-2](https://doi.org/10.1007/s10569-026-10278-2).
