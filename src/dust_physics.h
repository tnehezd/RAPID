// dust_physics.h

#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h>
#include "simulation_types.h" // Required for the simulation_options_t and disk_t structures
#include "particle_data.h"    // This include must also be added for dust_particle_t!


// Function prototype for calculating dust surface density profile
// This function performs a Cloud-in-Cell (CIC) mapping from Lagrangian particles
// onto an Eulerian grid.
void calculate_dust_surface_density_profile(
    double *output_sigma_d_grid,      // Output: The calculated dust surface density (per grid cell)
    double *output_r_grid_centers,    // Output: Radial positions of the grid cell centers
    const double radin[][2],          // Input: Particle radial positions (const)
    const double *massin,             // Input: Particle masses (const)
    int n_particles,                  // Input: Number of particles
    int n_grid_cells,                 // Input: Number of grid cells for density calculation
    const disk_t *disk_params         // Input: The disk parameters (const)
);

/* Calculates the turbulent alpha parameter */
double calculate_turbulent_alpha(double r, const disk_t *disk_params);

/* Calculates the Stokes number for a given particle */
double calculate_stokes_number(double pradius, double sigma, double r, const disk_t *disk_params);

// Functions to determine the maximum particle size based on various physical constraints (Birnstiel et al. 2012)
double calculate_max_size_from_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params);
double calculate_max_size_from_turbulence(double sigma, double r, double rho_p, const disk_t *disk_params);
double calculate_max_size_from_drift_fragmentation(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/* The timescale for particle growth */
double calculate_growth_timescale(double r, double eps, const disk_t *disk_params);

/* Calculates the new particle size at a given location, based on the Birnstiel paper */
double update_particle_size(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params);

// Function to calculate the dust disk density
// This function gets the ParticleData_t structure from tIntegrate.
// The sigmad and sigmadm values will be stored in disk_params.
void calculate_dust_density_grid(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts);

/* Function to store the new distances of the dust particles */
// NEW PROTOTYPE:
// - dust_particle_t *particles_array: The actual particle array (either Pop1 OR Pop2)
// - int num_particles: The number of particles in the array
// - double deltat, double t: timestep and current time
// - const simulation_options_t *sim_opts, const disk_t *disk_params: options and disk parameters
void update_particle_positions(dust_particle_t *particles_array, int num_particles, double deltat, double t,
                 const simulation_options_t *sim_opts, const disk_t *disk_params);


#endif // DUST_PHYSICS_H