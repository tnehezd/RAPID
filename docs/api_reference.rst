.. _api-reference:

API Reference
=============

This section provides a comprehensive, automatically generated API reference for the core C code components of the simulation. It details the functions, variables, and data structures used throughout the project, facilitating understanding and development.

----

Configuration Module
--------------------

This section details the API for the `config.h`, `parser.h` and the `simulation_types.h` header files. Together, these files include declarations for global simulation parameters, file pointers, and constant definitions for the setup and management of the simulation environment.

config.h
^^^^^^^^

This section presents the detailed API for the `config.h` header file, including global variables and file pointers related to simulation setup and configuration.

.. doxygenfile:: config.h
   :project: rapid


parser.h
^^^^^^^^

This section presents the detailed API for the `parser.h` header file, including global variables and file pointers related to simulation setup and configuration.

.. doxygenfile:: parser.h
   :project: rapid



simulation_types.h
^^^^^^^^^^^^^^^^^^^

This section presents the detailed API for the `simulation_types.h` header file, including global variables and file pointers related to simulation setup and configuration.

.. doxygenfile:: simulation_types.h
   :project: rapid


----


Initialize Disk
---------------

This section describes the initialization of disk and the pyhsical paramters of the simumation in the `init_tool_module.h` header file.

init_tool_module.h
^^^^^^^^^^^^^^^^^^

.. doxygenfile:: init_tool_module.h
   :project: rapid

----

Disk Physics Module
-------------------

This section details the API for the ``disk_model.h``, ``dust_physics.h``, ``gas_physics.h`` and the ``particle_data`` header files. These files provide functions for constructing the initial physical state of the protoplanetary gas disk and the containing dusty material.


disk_model.h
^^^^^^^^^^^^
This section details the API for the ``disk_model.h`` header file. It provide functions for constructing the initial physical state of the protoplanetary gas disk.
These routines generate the radial grid, surface density profile, pressure, pressure gradient, and gas velocity fields used throughout the simulation.


.. doxygenfile:: disk_model.h
   :project: rapid


----


dust_physics.h
^^^^^^^^^^^^^^

This section details the API for the ``dust_physics.h`` header file. It contains routines that compute the physical evolution of dust particles in the disk, including radial drift, aerodynamic coupling to the gas, particle growth barriers.


.. doxygenfile:: dust_physics.h
   :project: rapid

----


gas_physics.h
^^^^^^^^^^^^^

This section details the API for the ``gas_physics.h`` header file. It provides functions that compute the physical evolution of the gas disk, including viscosity, pressure scale height, gas velocity, and the time‑dependent update of the gas surface density.


.. doxygenfile:: gas_physics.h
   :project: rapid


particle_data.h
^^^^^^^^^^^^^^^

This section details the API for the ``particle_data.h`` header file. This module initializes and loads arrays connected to the dust particles.

.. doxygenfile:: particle_data.h
   :project: rapid



----

Boundary Conditions Module
--------------------------


boundary_conditions.h
^^^^^^^^^^^^^^^^^^^^^

This section details the API for the ``boundary_conditions.h`` header file. It provides functions for managing ghost-cell extrapolation and enforcing physical boundary conditions at the inner and outer edges of the computational domain to ensure numerical stability.


.. doxygenfile:: boundary_conditions.h
   :project: rapid

----


Simulation Core Module
----------------------

integrator.h
^^^^^^^^^^^^

.. doxygenfile:: integrator.h
   :project: rapid


simulation_core.h
^^^^^^^^^^^^^^^^^

.. doxygenfile:: simulation_core.h
   :project: rapid

----



Utils Module
------------


I/O Utils
^^^^^^^^^

io_utils.h
^^^^^^^^^^

.. doxygenfile:: io_utils.h
   :project: rapid



utils.h
^^^^^^^

.. doxygenfile:: utils.h
   :project: rapid
