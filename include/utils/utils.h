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
 * @brief Computes the dust mass contained within a specific pressure trap.
 *
 * Integrates the dust mass inside the trap boundaries using particle data.
 *
 * @param trap Pointer to the trap structure to update.
 * @param particle_data Pointer to particle data.
 * @param particle_number Number of particles.
 * @param sim_opts Pointer to simulation options.
 */
void calculateMassInSpecificTrap(PressureTrap *trap, const ParticleData *particle_data, int particle_number, const SimulationOptions *sim_opts);

/**
 * @brief Updates the grid index of each particle based on its radial position.
 *
 * @param particle_data Pointer to particle data.
 * @param actual_time Current simulation time.
 * @param particle_number Number of particles.
 * @param disk_params Pointer to disk parameters.
 * @param is_twopop Flag indicating if the two-population model is enabled.
 */
void updateParticleGridIndices(ParticleData *particle_data, double actual_time, int particle_number, const DiskParameters *disk_params, int is_twopop);

/**
 * @brief Computes the minimum and maximum particle radii in the system.
 *
 * Optionally includes micron‑sized particles if enabled.
 *
 * @param particle_data Pointer to particle data.
 * @param particle_number Number of particles.
 * @param has_secondary_population Non‑zero if micron dust is included.
 * @param min_radius Output: minimum radius.
 * @param max_radius Output: maximum radius.
 */
void computeParticleRadiusRange(
    const ParticleData *particle_data,
    int particle_number,
    int has_secondary_population,
    double *min_radius,
    double *max_radius
);

/**
 * @brief Checks if the secondary dust population is enabled.
 *
 * @param mode Snapshot mode.
 * @return Non-zero if the secondary dust population is enabled.
 */
int isSecondaryPopulationEnabled(SnapshotMode mode);

/**
 * @brief Checks if gas evolution is enabled.
 *
 * @param mode Snapshot mode.
 * @return Non-zero if gas evolution is enabled.
 */
int isGasEvolutionEnabled(SnapshotMode mode);

/**
 * @brief Checks if dust is enabled.
 *
 * @param mode Snapshot mode.
 * @return Non-zero if dust is enabled.
 */
int isDustEnabled(SnapshotMode mode);


/**
 * @brief Checks if dust growth is enabled.
 *
 * @param mode Snapshot mode.
 * @return Non-zero if dust growth is enabled.
 */
int isDustGrowthEnabled(SnapshotMode mode);

/**
 * @brief Checks if dust drift is enabled.
 *
 * @param mode Snapshot mode.
 * @return Non-zero if dust drift is enabled.
 */
int isDustDriftEnabled(SnapshotMode mode);

#endif // UTILS_H