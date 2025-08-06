/**
 * @file dust_physics.h
 * @brief This header file contains function prototypes and structures
 * related to the physical processes of dust particles in the RAPID model,
 * including dynamics, density calculations, and growth.
 */

#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h>
#include "simulation_types.h"
#include "particle_data.h"


/**
 * @brief Calculates the dust surface density profile on an Eulerian grid.
 *
 * This function performs a Cloud-in-Cell (CIC) mapping from Lagrangian particles
 * onto an Eulerian grid to compute the dust surface density.
 *
 * @param output_sigma_d_grid Output: The calculated dust surface density (per grid cell).
 * @param output_r_grid_centers Output: Radial positions of the grid cell centers.
 * @param particle_radius_au_array Input: A 1D array of particle radial positions (const).
 * @param particle_mass_array Input: A 1D array of particle masses (const).
 * @param n_particles Input: The total number of particles.
 * @param n_grid_cells Input: The number of grid cells for density calculation.
 * @param disk_params Input: A pointer to the disk parameters structure (const).
 *  
 * THERE'S A BUG!!!!
 */
void calculate_dust_surface_density_profile(double *output_sigma_d_grid, double *output_r_grid_centers, const double *particle_radius_au_array, const double *particle_mass_array, int n_particles, int n_grid_cells, const disk_t *disk_params);


/**
 * @brief Calculates the Stokes number for a given particle.
 * @param pradius The radius of the particle.
 * @param sigma The gas surface density.
 * @param r The radial distance from the star.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The calculated Stokes number.
 */
double calculate_stokes_number(double pradius, double sigma, double r, const disk_t *disk_params);

/**
 * @brief Calculates the maximum particle size limited by radial drift.
 *
 * This function determines the maximum particle size based on the drift barrier
 * described in Birnstiel et al. (2012).
 *
 * @param sigmad The dust surface density.
 * @param r The radial distance.
 * @param p The gas surface density power-law exponent.
 * @param dp The pressure gradient parameter.
 * @param rho_p The intrinsic density of the dust particle.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The maximum particle size from the drift barrier.
 */
double calculate_max_size_from_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/**
 * @brief Calculates the maximum particle size limited by turbulence.
 *
 * This function determines the maximum particle size based on the turbulence barrier
 * described in Birnstiel et al. (2012).
 *
 * @param sigma The gas surface density.
 * @param r The radial distance.
 * @param rho_p The intrinsic density of the dust particle.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The maximum particle size from the turbulence barrier.
 */
double calculate_max_size_from_turbulence(double sigma, double r, double rho_p, const disk_t *disk_params);

/**
 * @brief Calculates the maximum particle size limited by combined drift and fragmentation.
 *
 * This function determines the maximum particle size based on the fragmentation barrier
 * and radial drift, as described in Birnstiel et al. (2012).
 *
 * @param sigma The gas surface density.
 * @param r The radial distance.
 * @param p The gas surface density power-law exponent.
 * @param dp The pressure gradient parameter.
 * @param rho_p The intrinsic density of the dust particle.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The maximum particle size from the drift-fragmentation barrier.
 */
double calculate_max_size_from_drift_fragmentation(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/**
 * @brief Calculates the timescale for particle growth.
 * @param r The radial distance.
 * @param eps The dust-to-gas ratio.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The growth timescale.
 */
double calculate_growth_timescale(double r, double eps, const disk_t *disk_params);

/**
 * @brief Updates the particle size based on the Birnstiel et al. (2012) paper.
 * @param prad The current particle radius.
 * @param pdens The particle density.
 * @param sigma The gas surface density.
 * @param sigmad The dust surface density.
 * @param y The gas-to-dust ratio.
 * @param p The gas surface density power-law exponent.
 * @param dpress_val The pressure gradient value.
 * @param dt The time step.
 * @param disk_params A pointer to the disk parameters structure.
 * @return The new particle size.
 */
double update_particle_size(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params);

/**
 * @brief Calculates the dust disk density on a grid.
 *
 * This function retrieves the ParticleData_t structure from tIntegrate and
 * stores the calculated `sigmad` and `sigmadm` values in the disk_params structure.
 *
 * @param p_data A pointer to the ParticleData_t structure (const).
 * @param disk_params A pointer to the disk parameters structure.
 * @param sim_opts A pointer to the simulation options structure (const).
 */
void calculate_dust_density_grid(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts);

/**
 * @brief Updates the radial positions of the dust particles.
 *
 * This function calculates and stores the new distances of the dust particles.
 *
 * @param particles_array A pointer to the array of dust particles (either Pop1 or Pop2).
 * @param num_particles The number of particles in the array.
 * @param deltat The time step.
 * @param t The current simulation time.
 * @param sim_opts A pointer to the simulation options structure (const).
 * @param disk_params A pointer to the disk parameters structure (const).
 */
void update_particle_positions(dust_particle_t *particles_array, int num_particles, double deltat, double t, const simulation_options_t *sim_opts, const disk_t *disk_params);


#endif // DUST_PHYSICS_H