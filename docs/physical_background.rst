Physical Background
===================

The RAPID code (Representative Approach for Particle-Integrated Disks) is designed to simulate the coupled evolution of gas and dust in protoplanetary disks. It uses a 1D (radial) framework assuming an axisymmetric, vertically thin disk.

Numerical Approach: Hybrid Euler-Lagrange
-----------------------------------------
RAPID combines two distinct numerical descriptions:
1. **Eulerian Grid (Gas):** The gas surface density is evolved on a fixed radial grid using a finite-difference scheme.
2. **Lagrangian Particles (Dust):** The solid phase is modeled using an ensemble of $N_p$ representative particles. Each particle represents the dust mass of a specific annulus and is tracked individually through the evolving gas disk.

Gas Dynamics
------------
The evolution of the gas surface density :math:`\Sigma_g` is governed by the viscous evolution equation:

.. math::
   \frac{\partial \Sigma_g}{\partial t} = \frac{3}{r} \frac{\partial}{\partial r} \left[ r^{1/2} \frac{\partial}{\partial r} \left(\nu \Sigma_g r^{1/2}\right) \right]

The kinematic viscosity :math:`\nu` follows the Shakura-Sunyaev $\alpha$-prescription:
:math:`\nu = \alpha H c_s`.



Dead Zone Implementation
^^^^^^^^^^^^^^^^^^^^^^^^
A key feature of RAPID is the explicit implementation of an accretionally inactive **Dead Zone**. The viscosity is reduced by a factor :math:`\delta_\alpha` using a combined hyperbolic tangent profile:

.. math::
   \delta_{\alpha} = 1 - \frac{1}{2}\left ( 1 - \alpha_{\mathrm{mod}}  \right ) \left [  \tanh{\left ( \frac{r - r_{\mathrm{dze,i}}}{\Delta r_{\mathrm{dze,i}}} \right ) }+  \tanh{\left ( \frac{r_{\mathrm{dze,o}} - r}{\Delta r_{\mathrm{dze,o}}} \right )} \right]

This creates sharp viscosity gradients at the inner (:math:`r_{dze,i}`) and outer (:math:`r_{dze,o}`) edges, leading to the formation of **pressure traps**.

Dust Dynamics and Growth
------------------------
Dust particles experience a drag force due to the velocity difference between the sub-Keplerian gas and the Keplerian dust.

Radial Drift
^^^^^^^^^^^^
The radial velocity of a dust particle :math:`u_{dust,r}` is given by:

.. math::
   u_{dust,r} = \frac{u_{gas,r}}{1+St^2} + \frac{2}{St + St^{-1}}u_{drift}

where :math:`St` is the Stokes number, representing the coupling strength between dust and gas.



Particle Size Evolution (Growth Barriers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RAPID follows the model of Birnstiel et al. (2012), where particle growth is limited by three physical barriers:

1. **Fragmentation Barrier:** Caused by turbulent relative velocities.
2. **Drift Barrier:** Occurs when the drift timescale is shorter than the growth timescale.
3. **Drift-induced Fragmentation:** Relevant in low-turbulence regions like the Dead Zone.

The actual particle size :math:`a` is the minimum of the exponential growth value and these limiting barriers:

.. math::
   a = \min(a_{frag}, a_{drift}, a_{df}, a_0 e^{t/\tau_{grow}})