.. _api-reference:

API Reference
=============

This section provides a comprehensive, automatically generated API reference for the core C code components of the simulation. It details the functions, variables, and data structures used throughout the project, facilitating understanding and development.

----

Configuration Module
--------------------

This section details the API for the `config.h` header file. Together, these files include declarations for global simulation parameters, file pointers, and constant definitions for the setup and management of the simulation environment.

config.h
^^^^^^^^

This section presents the detailed API for the `config.h` header file, including global variables and file pointers related to simulation setup and configuration.

.. doxygenfile:: config.h
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

This section details the API for the ``disk_model.h`` header file. It provides functions for constructing the initial physical state of the protoplanetary gas disk.
These routines generate the radial grid, surface density profile, pressure, pressure gradient, and gas velocity fields used throughout the simulation.

disk_model.h
^^^^^^^^^^^^

.. doxygenfile:: disk_model.h
   :project: rapid


----

This section details the API for the ``dust_physics.h`` header file. It contains routines that compute the physical evolution of dust particles in the disk, including radial drift, aerodynamic coupling to the gas, particle growth barriers.


dust_physics.h
^^^^^^^^^^^^^^

.. doxygenfile:: dust_physics.h
   :project: rapid

----

This section details the API for the ``gas_physics.h`` header file. It provides functions that compute the physical evolution of the gas disk, including viscosity, pressure scale height, gas velocity, and the time‑dependent update of the gas surface density.



gas_physics.h
^^^^^^^^^^^^^

.. doxygenfile:: gas_physics.h
   :project: rapid



----

Utils Module
------------

I/O Utils
^^^^^^^^^