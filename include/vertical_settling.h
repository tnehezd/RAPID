#ifndef VERTICAL_SETTLING_H
#define VERTICAL_SETTLING_H

#include "simulation_types.h"
#include "config.h"
#include "gas_physics.h"

/**
 * @brief Calculate the vertical scale height (Hd) of a dust particle.
 * 
 * This function computes the dust scale height based on the local gas scale height,
 * the turbulence parameter alpha, and the particle's Stokes number.
 * 
 * @param radial_distance Radial distance from the star (AU or code units)
 * @param particle_radius Radius of the dust particle (cm or code units)
 * @param disk_params Pointer to disk parameters structure
 * @return Vertical dust scale height Hd
 */
double calculateDustScaleHeight(double radial_distance, double particle_radius, const DiskParameters *disk_params);

/**
 * @brief Compute the vertical density profile of dust particles.
 * 
 * Assumes a Gaussian vertical distribution based on the scale height Hd.
 * 
 * @param Hd Dust scale height
 * @param dust_surface_density Dust surface density at this radial location
 * @param rho_z Output array containing the vertical density profile (size z_points)
 * @param z_points Number of vertical grid points
 * @param z_max Maximum vertical distance from the midplane (positive value)
 */
void calculateVerticalDustProfile(double Hd, double dust_surface_density, double *rho_z, int z_points, double z_max);

/**
 * @brief Convenience function to calculate both Hd and vertical density profile.
 * 
 * Computes the dust scale height and fills the rho_z array with the vertical
 * Gaussian distribution of dust density.
 * 
 * @param radial_distance Radial distance from the star
 * @param particle_radius Radius of the dust particle
 * @param dust_surface_density Dust surface density
 * @param disk_params Pointer to disk parameters
 * @param Hd Output: computed dust scale height
 * @param rho_z_array Output array for vertical density profile
 * @param z_points Number of vertical grid points
 * @param z_max Maximum vertical distance from the midplane
 */
void calculateVerticalDistribution(double radial_distance, double particle_radius, double dust_surface_density, const DiskParameters *disk_params, double *Hd, double *rho_z_array, int z_points, double z_max);

#endif // VERTICAL_SETTLING_H
