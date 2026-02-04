.. Your C Project Name documentation master file, created by
   sphinx-quickstart.

Welcome to the Documentation of the RAPID Simulation Code
=========================================================

**RAPID** (Representative Approach for Particle-Integrated Disks) is a high-performance 1D numerical framework 
designed to simulate the co-evolution of gas and dust in protoplanetary disks using a hybrid Euler-Lagrange approach. 

The code specifically targets the formation of pressure traps at dead zone edges, enabling the investigation of 
dust accumulation and planetesimal formation nurseries with high computational efficiency.

----


.. toctree::
   :maxdepth: 2
   :caption: Contents

   self          
   intro.rst
   getting_started.rst
   physical_background.rst

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   api_reference


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


----


How to Cite
-----------

If you use the **RAPID** code for your research, publications, or presentations, we kindly ask you to cite the following papers to support the development of this project.

Primary Reference for the Code
------------------------------
The main reference describing the **RAPID** model, the trajectory-based Lagrangian approach, and its implementation is:

* **Tarczay-Nehéz, D. (2026)**: *Trajectory-based dust evolution in disks: first results from the RAPID simulation code*, Celestial Mechanics and Dynamical Astronomy, 138(1), 6. `https://doi.org/10.1007/s10569-026-10278-2`

Physical Foundation and Dead Zone Model
----------------------------------------
For the physical background regarding the dead zone edges, viscosity transitions, and the vortex-vortex interaction models implemented in RAPID, please cite:

* **Regály, Zs., Juhász, A., & Nehéz, D. (2017)**: *Interpreting Brightness Asymmetries in Transition Disks: Vortex at Dead Zone or Planet-carved Gap Edges?*, The Astrophysical Journal, 851(2), 89. `https://doi.org/10.3847/1538-4357/aa9a3f`

BibTeX
------
For your convenience, you can use the following BibTeX entries:

.. code-block:: bibtex

    @article{TarczayNehez2026,
        author = {Tarczay-Neh{\'e}z, D{\'o}ra},
        title = {Trajectory-based dust evolution in disks: first results from the RAPID simulation code},
        journal = {Celestial Mechanics and Dynamical Astronomy},
        year = {2026},
        volume = {138},
        doi = {10.1007/s10569-026-10278-2}
    }

    @article{Regaly2017,
        author = {{Reg{\'a}ly}, Zs. and {Juh{\'a}sz}, A. and {Neh{\'e}z}, D.},
        title = "{Interpreting Brightness Asymmetries in Transition Disks: Vortex at Dead Zone or Planet-carved Gap Edges?}",
        journal = {The Astrophysical Journal},
        year = 2017,
        volume = {851},
        doi = {10.3847/1538-4357/aa9a3f}
    }