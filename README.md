# RAPID: Representative Approach for Particle-Integrated Disks

This repository contains the source code for `RAPID` (Representative Approach for Particle-Integrated Disks), a 1D numerical model developed as part of my PhD research. The model is designed to simulate the dynamics and evolution of dust particles in protoplanetary disks, including dust traps at the edges of an embedded deadzone.

The core of the code was originally written in C during my PhD between 2012 and 2015. It has since been restructured into a more modular form to improve its functionality and maintainability.

#### ⚠️ 2.2.3 Version: Stable Structure – Ongoing Feature Development

The core structure of the codebase has now been finalized, and no further refactoring is planned. The architecture, naming conventions, and overall organization are stable and consistent with the [Coding Standard](docs/CodingStandard.md).

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

### Running a Simulation

There are two primary ways to run a simulation:

Using a configuration file:
The `run_simulation.py` Python3 wrapper script can be used to run a series of simulations based on parameters defined in a `config.yaml` file. Simply execute the script in your terminal:

```Bash
python run_simulation.py
```

Using command-line arguments:
Alternatively, you can run the executable directly from the command line by specifying the parameters. A detailed list of all available command-line flags and their corresponding parameters will be provided here shortly.

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

