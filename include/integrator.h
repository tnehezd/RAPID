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

void integrateParticleRungeKutta4Structured(double actual_time, double particle_radius,
                                  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, 
                                  const DiskParameters *disk_params, const SimulationOptions *simulation_options);

#endif // INTEGRATOR_H

