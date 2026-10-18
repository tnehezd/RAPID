# RAPID: Representative Approach for Particle-Integrated Disks

This repository contains the source code for `RAPID` (Representative Approach for Particle-Integrated Disks), a 1D numerical model developed as part of my PhD research. The model is designed to simulate the dynamics and evolution of dust particles in protoplanetary disks, including dust traps at the edges of an embedded deadzone.

The core of the code was originally written in C during my PhD between 2012 and 2015. It has since been restructured into a more modular form to improve its functionality and maintainability.

#### ⚠️ 2.3.0 Version: Stable Structure – Ongoing Feature Development

The core structure of the codebase has now been finalized, and no further refactoring is planned. The architecture, naming conventions, and overall organization are stable and consistent with the [Coding Standard](CodingStandard.md).

However, the project remains under active development in terms of features and capabilities. Please always make sure to `git pull` and stay up to date with the latest release before running simulations.

Upcoming and ongoing extensions include:

- Additional numerical integrators
- Dust settling implementation  
  (transition from 1D to 1+1D structure)
- General validation and benchmark test problems  
  (e.g. shock tube and other standard test cases)
  
Stay tuned for regular feature updates.

## Model Description

RAPID is a trajectory-based Lagrangian model that tracks the motion of representative dust particles in the radial direction within an evolving 1D Eulerian gaseous disk. A key feature of this model is that the pressure maximum, which acts as a dust trap, is not a static element but a dynamically evolving phenomenon arising naturally in the vicinity of the embedded dead zone's inner and outer edges.

This approach offers computational efficiency due to its 1D nature while providing an effective approximation of particle trajectories. Using RAPID, we can conduct parameter studies to investigate dust mass growth, particle trajectories, and the evolution of feeding zones around pressure bumps, aiming to understand the conditions that facilitate planetesimal formation in protoplanetary disks.

This work is based on the methods described in Regály et al. (2017) and Tarczay-Nehéz (2026).

## Building and Running the Simulation

### Prerequisites
To build the executable, you need a C compiler (e.g., GCC or Clang) and GNU Make.

### Build Instructions

To compile the source code, open your terminal and run the following commands:

```Bash
make clean
make all
```

`make all` defaults to the optimized **release** build (`-O2 -g`), without
`-ffast-math`. For a non-optimized debug build, use `make debug-build`
(`-O0 -g`); `make release` switches back. Equivalently, use
`make BUILD=debug all` or `make BUILD=release all`.

Objects and executables are kept separately under `obj/debug`, `obj/release`,
`bin/debug` and `bin/release`, so switching profiles does not require
`make clean`. Each build selects its executable at `bin/simulation`, preserving
the Python wrapper's default path. Use the profile-specific executables to run
both versions without switching. Do not select different profiles concurrently
in the same worktree. `make debug` builds the debug profile and opens LLDB.
`make test-thermal` uses the selected `BUILD` profile as well.

### Running a Simulation

There are two primary ways to run a simulation:

Using a configuration file:
The `run_simulation.py` Python3 wrapper script can be used to run a series of simulations based on parameters defined in a `config.yaml` file. Simply execute the script in your terminal:

```Bash
python run_simulation.py
```

Using command-line arguments:
Alternatively, you can run the executable directly from the command line by specifying the parameters. A detailed list of all available command-line flags and their corresponding parameters will be provided here shortly.

Gas profile snapshots (`LOGS/density_profile_*.dat`) write radius, surface density,
pressure, and pressure gradient with 10 significant digits. Rebuild and rerun
simulations to obtain this precision; existing snapshots are unchanged.

### Runtime comparison matrix

`test_examples/generate_configs.py` generates 80 independent YAMLs: SG on/off,
photoevaporation on/off (Picogna), gas-only/one growing dust population,
linear/logarithmic grids, and five thermal choices (lociso, radiative with
prescribed/geometric irradiation, radiative_energy with prescribed/geometric
irradiation). Each uses 500 gas cells, 1000 particles when dust is active,
4 OpenMP threads, 500,000 years and 10,000-year output cadence.
Energy runs use outflow boundaries; other modes use zero-gradient boundaries,
so these comparisons do not isolate numerical cost at identical physics.

```bash
python test_examples/generate_configs.py --generate-only
# Generate and run all cases, at most 10 concurrent runs / 40 OpenMP threads:
python test_examples/generate_configs.py
# Reduce concurrency:
python test_examples/generate_configs.py --jobs 2
```

Build the release binary first (`make all`). Each case has its own output
directory; per-case console logs and `timings.tsv` are written to
`benchmark_logs/` (overwritten on a new batch). Failed cases are reported and
the runner exits nonzero after the batch. Concurrent wall times include resource
contention; use `--jobs 1` for isolated timing comparisons.

### Initial dust drift diagnostic

Dust runs write `drift_timecale.dat` (historical filename retained) in the
output root, for both ASCII and HDF5 output. It contains the initial population
and particle IDs, distance in AU, grain radius in cm, Stokes number, signed
radial velocity in AU/year, local `r/abs(v_r)` in years, and direction
(-1 inward, 0 stationary, +1 outward). A stationary particle has infinite
timescale. Velocity includes gas advection: this is not a time-integrated
arrival time at the star. Both populations are listed when enabled.

Growth/fragmentation updates the primary grain size once per step, after dust
density deposition and before RK4 position transport. Gas and deposited dust
fields are interpolated at the particle radius on their gas-grid coordinates;
the dust density's zero-based physical-cell array excludes gas ghost cells.
RK4 holds that size fixed. The secondary micron population is not grown by
this update. This split scheme is not a fourth-order coupled growth solver.
`make test-dust` tests linear/log grids with particle counts below, equal to
and above the gas-cell count.

Gas viscous evolution advances `nu * Sigma` and divides by the same local
viscosity to recover `Sigma` before applying photoevaporation and the density
floor. Results generated without this conversion must be rerun.

----
### References

If you use this code in your research, please cite the following papers, which provides the theoretical background for the core numerical methods:

```bibtex
@ARTICLE{2017ApJ...851...89R,
      author = {{Reg{\'a}ly}, Zs. and {Juh{\'a}sz}, A. and {Neh{\'e}z}, D.},
      title = "{Interpreting Brightness Asymmetries in Transition Disks: Vortex at Dead Zone or Planet-carved Gap Edges?}",
    journal = {\apj},
   keywords = {accretion, accretion disks, hydrodynamics, methods: numerical, protoplanetary disks, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
      year = 2017,
      month = dec,
     volume = {851},
     number = {2},
       eid = {89},
      pages = {89},
        doi = {10.3847/1538-4357/aa9a3f},
archivePrefix = {arXiv},
     eprint = {1711.03548},
primaryClass = {astro-ph.EP},
     adsurl = {[https://ui.adsabs.harvard.edu/abs/2017ApJ...851...89R](https://ui.adsabs.harvard.edu/abs/2017ApJ...851...89R)},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@article{TarczayNehez2026,
    author = {Tarczay-Neh{\'e}z, D{\'o}ra},
    date = {2026/02/03},
    date-added = {2026-02-03 12:16:24 +0100},
    date-modified = {2026-02-03 12:16:24 +0100},
    doi = {10.1007/s10569-026-10278-2},
    id = {Tarczay-Neh{\'e}z2026},
    isbn = {1572-9478},
    journal = {Celestial Mechanics and Dynamical Astronomy},
    number = {1},
    pages = {6},
    title = {Trajectory-based dust evolution in disks: first results from the RAPID simulation code},
    url = {https://doi.org/10.1007/s10569-026-10278-2},
    volume = {138},
    year = {2026},
    bdsk-url-1 = {https://doi.org/10.1007/s10569-026-10278-2}}
}
