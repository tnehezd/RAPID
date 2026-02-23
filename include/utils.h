/**
 * @file utils.h
 * @brief Helper functions for interpolation, sorting, particle handling, and trap analysis.
 *
 * Contains numerical utility routines used across the simulation, including
 * interpolation, derivative coefficients, pressure‑trap detection, particle
 * merging and rounding, histogram updates, and radius‑range calculations.
 */


#ifndef UTILS_H
#define UTILS_H

#include "particle_data.h"
#include "simulation_types.h" 

/**
 * @brief Performs linear interpolation on a radial grid.
 *
 * Computes the interpolated value of a quantity defined on @p radial_grid
 * at an arbitrary position @p actual_position. The result is written to @p interpolated_value.
 *
 * @param input_value Input vector containing values defined on the radial grid.
 * @param radial_grid Radial grid array.
 * @param actual_position Position at which interpolation is required.
 * @param interpolated_value Output pointer for the interpolated value.
 * @param grid_spacing Grid spacing or auxiliary radial parameter.
 * @param disk_params Pointer to disk parameters.
 */
void linearInterpolation(double *input_value, double *radial_grid, double actual_position, double *interpolated_value, double grid_spacing, const DiskParameters *disk_params);

/**
 * @brief Returns the minimum of three values.
 *
 * @param value_1 First value.
 * @param value_2 Second value.
 * @param value_3 Third value.
 * @return The smallest of the three inputs.
 */
double findMinimumOfAnArray(double value_1, double value_2, double value_3);

/**
 * @brief Computes the FTCS second‑derivative coefficient at a given radius.
 *
 * Used in finite‑difference schemes for diffusion or pressure calculations.
 *
 * @param radial_distance Radial position.
 * @param disk_params Pointer to disk parameters.
 * @return The second‑derivative coefficient.
 */
double ftcsSecondDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the FTCS first‑derivative coefficient at a given radius.
 *
 * @param radial_distance Radial position.
 * @param disk_params Pointer to disk parameters.
 * @return The first‑derivative coefficient.
 */
double ftcsFirstDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Identifies pressure traps in the disk.
 *
 * A pressure trap is a location where the radial pressure gradient changes sign.
 *
 * @param disk_params Pointer to disk parameters.
 * @param traps Output array of detected traps.
 * @param max_traps Maximum number of traps that can be stored.
 * @return Number of traps found.
 */
int identifyPressureTraps(const DiskParameters *disk_params, PressureTrap *traps, int max_traps);


/**
 * @brief Sorts a 2D array with 3 columns by the first column.
 *
 * @param array 2D array of size n×3.
 * @param number_of_rows Number of rows.
 */
void sortAnArray(double array[][3],int number_of_rows);

/**
 * @brief Rounds particle radii to the nearest grid cell.
 *
 * @param particle_data 2D array containing particle radii and metadata.
 * @param particle_number Number of particles.
 * @param disk_params Pointer to disk parameters.
 */
void roundParticleRadii(double particle_data[][3], int particle_number, const DiskParameters *disk_params);

/**
 * @brief Merges particles that fall within the same radial bin.
 *
 * Used to reduce particle count by combining nearby particles.
 *
 * @param particle_data 2D array of particle data.
 * @param grid_step Radial bin width.
 * @param particle_number Number of particles.
 * @param disk_params Pointer to disk parameters.
 */
void mergeParticlesByRadius(double particle_data[][3], double grid_step, int particle_number, const DiskParameters *disk_params); 


void calculateMassInSpecificTrapStructured(PressureTrap *trap, const StructuredParticleData *data, const SimulationOptions *sim_opts);


void updateStructuredParticleGridIndices(StructuredParticleData *sd,  const DiskParameters *disk_params);


void updateDustSurfaceDensityStructured(StructuredParticleData *data, const DiskParameters *disk_params);

void updateDustSurfaceDensityEulerianCIC(StructuredParticleData *data, const DiskParameters *disk_params);

void updateDustSurfaceDensityEulerian(StructuredParticleData *data, const DiskParameters *disk_params);

void updateDustSurfaceDensityEulerianTSC(StructuredParticleData *data, const DiskParameters *disk_params);

void updateDustSurfaceDensitySmart(StructuredParticleData *data, const DiskParameters *disk_params);


#endif // UTILS_H