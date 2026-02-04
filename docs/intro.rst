Introduction
============

RAPID is a 1D numerical code designed to simulate the co-evolution of gas and dust in protoplanetary disks. 
The code employs a **hybrid Euler-Lagrange approach** to capture the different physical natures of the two components:

* **Gas Dynamics:** The gas phase is treated as a viscous fluid on a fixed **Eulerian grid**. It evolves according to the viscous evolution equation, accounting for gas pressure, Keplerian rotation, and optionally, disk winds or photoevaporation.
* **Dust Dynamics:** The solid phase is modeled using a **Lagrangian approach** with **representative particles**. 

Representative Particle Model
-----------------------------
Unlike global fluid models for dust, RAPID tracks discrete particles (or "super-particles"). 
Each Lagrangian particle represents a specific amount of dust mass based on its initial location. 
Throughout the simulation, the particle carries this mass and its physical properties (like grain size) 
as it moves radially due to gas drag and pressure gradients. This allows for a more precise tracking 
of dust-to-gas ratios and local accumulations.

Governing Equations
-------------------
To understand the underlying physics, the following equations are implemented:

**1. Gas Surface Density Evolution:**
The viscous evolution of the gas surface density :math:`\Sigma_g` is governed by:

.. math::

   \frac{\partial \Sigma_g}{\partial t} = \frac{3}{r} \frac{\partial}{\partial r} \left[ \sqrt{r} \frac{\partial}{\partial r} (\nu \Sigma_g \sqrt{r}) \right]

where :math:`\nu` is the kinematic viscosity.

**2. Dust Radial Velocity:**
Dust particles drift radially due to their coupling with the gas, described by the stopping time :math:`t_s` (or Stokes number :math:`St`):

.. math::

   v_{r,dust} = \frac{1}{1 + St^2} v_{r,gas} + \frac{St}{1 + St^2} \frac{\eta v_K}{St}

where :math:`\eta` is the pressure support parameter.