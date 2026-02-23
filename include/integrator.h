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
 * This function evolves both the radial position of the particle and its physical radius (size)
 * over a given timestep. It accounts for the aerodynamic drag (radial drift) and, 
 * if enabled, the dust growth/fragmentation processes by interpolating local 
 * gas properties from the Eulerian grid.
 *
 * @param[in]  actual_time            Current simulation time.
 * @param[in]  particle_radius        Current physical radius (grain size) of the dust particle.
 * @param[in]  dust_surfacedensity    Pointer to the array containing the dust surface density distribution.
 * @param[in]  particle_distance_grid Pointer to the grid used for particle-related spatial calculations.
 * @param[in]  actual_timestep        The time step size (dt) for the integration.
 * @param[in]  particle_distance      Current radial distance of the particle from the central star.
 * @param[out] particle_distance_new  Pointer to store the updated radial distance after the RK4 step.
 * @param[out] particle_radius_new    Pointer to store the updated grain size after the growth step.
 * @param[in]  disk_params            Pointer to the structure containing physical disk parameters and gas arrays.
 * @param[in]  simulation_options     Pointer to the simulation configuration (e.g., enabling growth models).
 *
 * @note The function uses four intermediate slope calculations (k1, k2, k3, k4) to ensure 
 * fourth-order accuracy in the radial trajectory.
 */
void integrateParticleRungeKutta4(double actual_time, double particle_radius, const double *dust_surfacedensity, const double *particle_distance_grid, 
								  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, 
								  const DiskParameters *disk_params, const SimulationOptions *simulation_options);

void integrateParticleRungeKutta4Structured(double actual_time, double particle_radius, const double *dust_surfacedensity, const double *particle_distance_grid, 
                                  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, 
                                  const DiskParameters *disk_params, const SimulationOptions *simulation_options);

#endif // INTEGRATOR_H

