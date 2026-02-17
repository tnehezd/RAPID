#ifndef VERTICAL_SETTLING_H
#define VERTICAL_SETTLING_H

#include "simulation_types.h"
#include "particle_data.h"


/**
 * @brief Computes the vertical dust distribution for a given radius and particle size.
 *
 * This function calculates the dust scale height and generates the vertical 
 * dust density profile assuming a Gaussian distribution.
 *
 * @param radial_distance The radial distance from the star [same units as in DiskParameters].
 * @param particle_radius The radius of the dust particle [cm].
 * @param dust_surface_density The dust surface density Σ_d at the given radius [g/cm^2].
 * @param disk_params Pointer to the DiskParameters structure containing disk properties.
 * @param Hd Pointer to a double where the calculated dust scale height will be stored.
 * @param rho_z_array Preallocated array of size `z_points` where the vertical dust density profile will be stored.
 * @param z_points Number of points in the vertical (z) direction.
 * @param N_Hd The number of scale heights above/below the midplane to compute (e.g., 4 → ±4 Hd).
 *
 * @note The function internally computes the dust scale height based on turbulence 
 *       and Stokes number, and then fills rho_z_array with a normalized Gaussian.
 *
 * @warning The rho_z_array must be preallocated with at least `z_points` elements.
 */
void calculateVerticalDistribution(double radial_distance,
                                   double particle_radius,
                                   double dust_surface_density,
                                   const DiskParameters *disk_params,
                                   double *Hd,
                                   double *rho_z_array,
                                   int z_points,
                                   double N_Hd);


/*
 * Apply vertical settling + diffusion for all radial and vertical points.
 * Updates the existing StructuredParticleData in place.
 *
 * Parameters:
 *  sdata           : pointer to structured particle data
 *  disk_params     : disk structure (contains gas density, scale height, etc.)
 *  particle_radius : radius of the dust particles [AU or cm consistent with Stokes calc]
 *  dt              : timestep [internal units]
 *  alpha_turb      : turbulence parameter (dimensionless)
 */
void applyVerticalSettling(StructuredParticleData *sdata, const DiskParameters *disk_params, double dt);


#endif
