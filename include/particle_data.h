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
} DustParticle;

/**
 * @brief Container for all particle‑related arrays used in the simulation.
 *
 * This structure stores dynamically allocated arrays that describe the
 * spatial distribution, mass distribution, and surface densities of both
 * regular dust particles and micron‑sized dust particles.  
 *
 * The arrays are allocated based on the number of particles and whether
 * the two‑population dust model is enabled.
 */
typedef struct {
    double (*particle_distance_array)[3];          /**< Radial positions of dust particles (2 columns: r and dr). */
    double (*micron_particle_distance_array)[3];   /**< Radial positions of micron‑sized dust particles. */
    double *dust_particle_mass_grid;               /**< Mass grid for dust particles. */
    double *massmicradial_grid;                    /**< Mass grid for micron‑sized dust particles. */
    double (*dust_particle_mass_array)[5];         /**< Dust particle mass evolution array (5 columns per particle). */
    double (*micron_dust_particle_mass_array)[5];  /**< Micron dust mass evolution array. */
    double *dust_surfacedensity;                   /**< Dust surface density profile. */
    double *micron_dust_surfacedensity;            /**< Micron dust surface density profile. */
    double *particle_distance_grid;                /**< Radial grid for dust particles. */
    double *micron_particle_distance_grid;         /**< Radial grid for micron dust particles. */
    double *particle_z_array;                      /**< NEW: Vertical z positions of each particle */
    size_t allocated_particle_number;              /**< Number of particles allocated in memory. */
} ParticleData;

/**
 * @brief Allocates memory for particle‑related arrays.
 *
 * This function initializes all arrays inside a ParticleData structure
 * based on the number of particles and whether the two‑population dust
 * model is enabled.
 *
 * @param particle_data Pointer to the ParticleData structure to initialize.
 * @param particle_count Number of dust particles to allocate.
 * @param is_twopop_enabled Non‑zero if the micron‑dust population should also be allocated.
 * @return 0 on success, non‑zero on allocation failure.
 */
int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled);

/**
 * @brief Frees all memory allocated inside a ParticleData structure.
 *
 * This function releases all dynamically allocated arrays and resets
 * the structure fields to NULL or zero.
 *
 * @param particle_data Pointer to the ParticleData structure to free.
 */
void freeParticleData(ParticleData *particle_data);

#endif // PARTICLE_DATA_H