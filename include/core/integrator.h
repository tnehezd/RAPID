/**
 * @file integrator.h
 * @brief Numerical integration routines for particle trajectories and growth.
 *
 * This header defines the integrators used to evolve the radial position and 
 * physical size of Lagrangian dust particles over time.
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "disk_model.h"      
#include "simulation_types.h"

/**
 * @brief Performs a fourth-order Runge-Kutta (RK4) integration step for a single dust particle.
 *
 * Evolves radial position at fixed grain size using interpolated gas properties.
 * Growth/fragmentation is applied once by updateParticleSizes before transport.
 *
 * @param[in]  particle_radius        Current physical radius (grain size) of the dust particle.
 * @param[in]  actual_timestep        The time step size (dt) for the integration.
 * @param[in]  particle_distance      Current radial distance of the particle from the central star.
 * @param[out] particle_distance_new  Pointer to store the updated radial distance after the RK4 step.
 * @param[out] particle_radius_new    Unchanged grain size.
 * @param[in]  disk_params            Pointer to the structure containing physical disk parameters and gas arrays.
 *
 * @note The function uses four intermediate slope calculations (k1, k2, k3, k4) to ensure 
 * fourth-order accuracy in the radial trajectory.
 */
void integrateParticleRungeKutta4(double particle_radius,
								  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, 
								  const DiskParameters *disk_params);

#endif // INTEGRATOR_H
