/**
 * @file simulation_core.h
 * @brief Core routines for advancing the disk–dust system in time.
 *
 * Contains the main physical update functions used during the simulation:
 * dust drift calculation, timestep determination, and the full time
 * integration loop that evolves gas and dust according to the selected
 * snapshot mode.
 */

#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

#include "disk_model.h"        
#include "simulation_types.h"  

/**
 * @brief Computes the radial drift velocity of dust particles in 1D.
 *
 * This function calculates the drift velocity of a dust particle based on
 * its size, the local pressure gradient, gas surface density, gas velocity,
 * and radial position. The result is written into @p drift_velocity.
 *
 * @param particle_radius Radius of the dust particle.
 * @param pressure_gradient Local radial pressure gradient.
 * @param gas_surface_density Gas surface density at the particle location.
 * @param gas_velocity Radial gas velocity.
 * @param radial_distance Radial distance of the particle from the star.
 * @param drift_velocity Output pointer for the computed drift velocity.
 * @param disk_params Pointer to the disk parameter structure.
 */
void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params);

/**
 * @brief Computes the simulation time step based on disk parameters.
 *
 * The time step is typically determined by stability constraints such as
 * Courant conditions or drift/growth timescales.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 * @return The computed time step.
 */
double calculateTimeStep(const DiskParameters *disk_params);

/**
 * @brief Performs the full time integration loop for the simulation.
 *
 * This function advances the system in time according to the selected
 * snapshot mode, updating disk quantities, dust evolution, and writing
 * output files when necessary.
 *
 * @param mode Snapshot writing mode (e.g., none, periodic, final only).
 * @param disk_params Pointer to the DiskParameters structure to update.
 * @param sim_opts Pointer to the SimulationOptions structure.
 * @param output_files Pointer to the OutputFiles structure for writing results.
 */
void timeIntegrationForTheSystem(SnapshotMode mode, DiskParameters *disk_params, SimulationOptions *sim_opts, OutputFiles *output_files);

#endif // SIMULATION_CORE_H