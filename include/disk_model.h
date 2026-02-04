/**
 * @file disk_model.h
 * @brief Initialization routines for the 1D protoplanetary gas disk model.
 *
 * This module provides the functions required to construct the initial physical
 * state of the gas disk used throughout the simulation. It generates the radial
 * grid, computes the initial gas surface density, pressure, pressure gradient,
 * and gas velocity fields. These quantities define the baseline disk structure
 * from which both gas evolution and dust dynamics proceed.
 *
 * The module also includes a helper routine for mapping representative dust
 * particle masses onto the gas grid, enabling consistent coupling between the
 * dust and gas components.
 *
 * All functions operate on the DiskParameters structure, which stores the
 * grid‑based physical quantities and global disk configuration parameters.
 */

#ifndef DISK_MODEL_H
#define DISK_MODEL_H

#include "config.h"
#include "simulation_types.h" 

/**
 * @brief Reads disk parameters from configuration and initializes the DiskParameters structure.
 *
 * This function loads the global disk configuration (radial bounds, grid size,
 * physical constants, and initial parameter values) into the DiskParameters
 * structure. It must be called before any grid‑based initialization routines.
 *
 * @param disk_params Pointer to the DiskParameters structure to be filled.
 */
void readDiskParameters(DiskParameters *disk_params);

/**
 * @brief Constructs the radial grid of the gas disk.
 *
 * The radial grid defines the spatial discretization used for all gas and dust
 * calculations. This function fills the radial_grid array and computes the
 * uniform grid spacing \f$\Delta r\f$.
 *
 * @param disk_params Pointer to the DiskParameters structure containing grid arrays.
 */
void createRadialGrid(DiskParameters *disk_params);

/**
 * @brief Computes the initial gas surface density profile.
 *
 * This routine initializes the gas surface density \f$\Sigma(r)\f$ according to the chosen
 * disk model (e.g., power‑law or tapered profile). The resulting array is used
 * by all subsequent gas and dust physics modules.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 */
void createInitialGasSurfaceDensity(DiskParameters *disk_params);

/**
 * @brief Computes the initial gas pressure profile.
 *
 * Using the initialized surface density and temperature (or scale height)
 * structure, this function evaluates the gas pressure \f$P(r)\f$ across the grid.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 */
void createInitialGasPressure(DiskParameters *disk_params);

/**
 * @brief Computes the initial radial pressure gradient \f$dP/dr\f$\f$.
 *
 * The pressure gradient is a key quantity for dust drift, gas velocity, and
 * stability analysis. This function computes the discrete derivative of the
 * pressure profile across the radial grid.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 */
void createInitialGasPressureGradient(DiskParameters *disk_params);

/**
 * @brief Computes the initial gas radial velocity profile.
 *
 * The gas velocity is typically derived from viscous accretion theory and
 * depends on the pressure gradient and viscosity. This function initializes
 * the gas_velocity array in the DiskParameters structure.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 */
void createInitialGasVelocity(DiskParameters *disk_params);

/**
 * @brief Maps representative dust particle masses onto the gas grid.
 *
 * This function reconstructs the dust surface density profile from a set of
 * representative particle radii and masses. It distributes particle mass onto
 * the gas grid to produce a grid‑based dust surface density array.
 *
 * @param input_dust_radii_array     		2D array of particle radii and positions.
 * @param input_mass_array            		Array of particle masses.
 * @param output_dust_surfacedensity_array 	2D array to store the resulting dust Σ(r).
 * @param particle_number             		Number of representative dust particles.
 * @param disk_params                 		ßPointer to the DiskParameters structure.
 */
void calculateDustSurfaceDensityFromRepresentativeMass(double input_dust_radii_array[][2], double *input_mass_array, double output_dust_surfacedensity_array[][3], int particle_number, const DiskParameters *disk_params);

#endif // DISK_MODEL_H
