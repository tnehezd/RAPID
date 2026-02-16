#ifndef VERTICAL_SETTLING_H
#define VERTICAL_SETTLING_H

#include "simulation_types.h"
#include "config.h"
#include "gas_physics.h"

/**
 * @brief Calculate the vertical scale height (Hd) of a dust particle.
 * 
 * Computes the dust scale height based on the local gas scale height,
 * turbulence parameter alpha, and the particle's Stokes number.
 * 
 * @param radial_distance Radial distance from the star (AU or code units)
 * @param particle_radius Radius of the dust particle (cm or code units)
 * @param disk_params Pointer to the disk parameters structure
 * @return Vertical dust scale height Hd
 */
double calculateDustScaleHeight(double radial_distance, double particle_radius, const DiskParameters *disk_params);

/**
 * @brief Compute the vertical density profile of dust particles.
 * 
 * Assumes a Gaussian vertical distribution based on the dust scale height Hd.
 * 
 * @param Hd Dust scale height
 * @param dust_surface_density Dust surface density at this radial location
 * @param rho_z Output array containing the vertical density profile (size z_points)
 * @param z_points Number of vertical grid points
 * @param z_max Maximum vertical distance from the midplane (positive value)
 */
void calculateVerticalDustProfile(double Hd, double dust_surface_density, double *rho_z, int z_points, double z_max);

/**
 * @brief Compute both dust scale height and vertical density profile.
 * 
 * This convenience function calculates Hd and fills the rho_z array with a
 * Gaussian vertical distribution of dust density. The vertical range is
 * automatically scaled by N_Hd * Hd.
 * 
 * @param radial_distance Radial distance from the star (AU or code units)
 * @param particle_radius Radius of the dust particle (cm or code units)
 * @param dust_surface_density Dust surface density at this radial location
 * @param disk_params Pointer to the disk parameters structure
 * @param Hd Output: computed dust scale height
 * @param rho_z_array Output array for the vertical density profile (size z_points)
 * @param z_points Number of vertical grid points
 * @param N_Hd Maximum vertical distance in units of Hd (e.g., 5 → ±5*Hd)
 */
void calculateVerticalDistribution(double radial_distance, double particle_radius, double dust_surface_density, const DiskParameters *disk_params, double *Hd, double *rho_z_array, int z_points, double N_Hd);

#endif // VERTICAL_SETTLING_H
