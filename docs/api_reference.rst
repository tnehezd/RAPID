.. _api-reference:

API Reference
=============

This section provides a comprehensive, automatically generated API reference for the core C code components of the simulation. It details the functions, variables, and data structures used throughout the project, facilitating understanding and development.

----

Configuration Module
--------------------

This section details the API for the `config.h` and `globals.h` header files. Together, these files include declarations for global simulation parameters, file pointers, and constant definitions for the setup and management of the simulation environment.

config.h
^^^^^^^^

This section presents the detailed API for the `config.h` header file, including global variables and file pointers related to simulation setup and configuration.

.. doxygenfile:: config.h
   :project: rapid

----

globals.h
^^^^^^^^^

This section presents the detailed API for the `globals.h` header file that contains global variables for conversions related to simulation setup and configuration.

.. doxygenfile:: globals.h
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

This section details the API for the `disk_model.h` header file. It includes declarations for functions that initialize, and evolve the physical properties of the gas disk within the simulation, such as surface density, pressure, velocity, and pressure gradient.

disk_model.h
^^^^^^^^^^^^

.. doxygenfile:: disk_model.h
   :project: rapid


This section details the API for the `dust_physics.h` header file. It includes declarations for functions that initialize, and evolve the physical properties of the dust particles within the simulation, such as radial distance, friction due to the headwind, size and the size constraining barriers.

dust_physics.h
^^^^^^^^^^^^^^

.. doxygenfile:: dust_physics.h
   :project: rapid



----

Utils Module
------------

I/O Utils
^^^^^^^^^