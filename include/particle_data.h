/**
 * @file particle_data.h
 * @brief Data structures and utilities for managing dust particle arrays.
 *
 * Defines the ParticleData container used to store particle positions,
 * masses, and surface densities, along with allocation and cleanup
 * routines used by the simulation.
 */

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <stddef.h>


/**
 * @brief Represents a single dust particle in the protoplanetary disk.
 *
 * Each particle stores its unique index, radial and vertical positions,
 * and its mass. Future extensions may include velocity, charge, or
 * other physical properties.
 */
typedef struct {
    size_t index;      /**< Unique identifier for the particle. */
    double r_au;       /**< Radial distance from the star in AU. */
    double z_au;       /**< Vertical height above the disk midplane in AU. */
    double mass_g;     /**< Particle mass in grams. */
    double radius;     /**< Particle radius in cm. */
    double grid_index; /**< Initial particle index. */
} DustParticle;

/**
 * @brief Container for dust particles on a 2D grid (radial x vertical).
 *
 * Stores particles as a dynamically allocated 2D array [n_r][n_z].
 */
typedef struct {
    DustParticle **particles; /**< 2D array of DustParticle structs. */
    size_t n_r;               /**< Number of radial grid points. */
    size_t n_z;               /**< Number of vertical cells per radial point. */
} StructuredParticleData;

void freeStructuredParticleData(StructuredParticleData *structured_data);

#endif // PARTICLE_DATA_H