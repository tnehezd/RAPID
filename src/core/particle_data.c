/**
 * @file particle_data.c
 * @brief Memory management for dust particle data structures.
 */

#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Required for memset
#include "logger.h"

/**
 * @brief Allocates memory for dust particle arrays based on population settings.
 * 
 * Uses memset to ensure all pointers are initialized to NULL before allocation.
 * If two-population mode is disabled, secondary arrays remain NULL.
 */
int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled) {

    if (particle_data == NULL) {
        LOG_ERROR("ParticleData pointer is NULL.\n");
        return 1; 
    }

    // Initialize all pointers to NULL and numeric fields to 0
    memset(particle_data, 0, sizeof(ParticleData));

    if (particle_count == 0) {
        LOG_DEBUG("Particle count is 0. No particle arrays allocated.\n");
        return 0; 
    }

    // 1. Primary population arrays
    particle_data->particle_distance_array     = calloc(particle_count, sizeof(*particle_data->particle_distance_array));
    particle_data->dust_particle_mass_grid     = calloc(particle_count, sizeof(double));
    particle_data->dust_particle_mass_array    = calloc(particle_count, sizeof(*particle_data->dust_particle_mass_array));
    particle_data->dust_surfacedensity         = calloc(particle_count, sizeof(double));
    particle_data->particle_distance_grid      = calloc(particle_count, sizeof(double));

    // 2. Secondary population arrays (only if enabled)
    if (is_twopop_enabled) {
        particle_data->micron_particle_distance_array = calloc(particle_count, sizeof(*particle_data->micron_particle_distance_array));
        particle_data->massmicradial_grid            = calloc(particle_count, sizeof(double));
        particle_data->micron_dust_particle_mass_array = calloc(particle_count, sizeof(*particle_data->micron_dust_particle_mass_array));
        particle_data->micron_dust_surfacedensity      = calloc(particle_count, sizeof(double));
        particle_data->micron_particle_distance_grid   = calloc(particle_count, sizeof(double));
    }

    // 3. Error checking for primary arrays
    if (!particle_data->particle_distance_array || !particle_data->dust_particle_mass_grid || 
        !particle_data->dust_particle_mass_array || !particle_data->dust_surfacedensity || 
        !particle_data->particle_distance_grid) {
        LOG_ERROR("Primary particle array allocation failed!\n");
        freeParticleData(particle_data); 
        return 1; 
    }

    // 4. Error checking for secondary arrays
    if (is_twopop_enabled) {
        if (!particle_data->micron_particle_distance_array || !particle_data->massmicradial_grid || 
            !particle_data->micron_dust_particle_mass_array || !particle_data->micron_dust_surfacedensity || 
            !particle_data->micron_particle_distance_grid) {
            LOG_ERROR("Secondary particle array allocation failed!\n");
            freeParticleData(particle_data);
            return 1;
        }
    }

    particle_data->allocated_particle_number = particle_count;
    LOG_DEBUG("Particle arrays allocated (%zu particles, Two-Pop: %s).\n", 
              particle_count, is_twopop_enabled ? "ON" : "OFF");
    
    return 0;
}

/**
 * @brief Frees all allocated particle arrays.
 * 
 * Safely calls free() on all members. free(NULL) is a safe no-op in C,
 * so this handles both single and two-population modes seamlessly.
 */
void freeParticleData(ParticleData *particle_data) {

    if (particle_data == NULL) {
        return; 
    }

    // Free primary population
    free(particle_data->particle_distance_array);
    free(particle_data->dust_particle_mass_grid);
    free(particle_data->dust_particle_mass_array);
    free(particle_data->dust_surfacedensity);
    free(particle_data->particle_distance_grid);

    // Free secondary population
    free(particle_data->micron_particle_distance_array);
    free(particle_data->massmicradial_grid);
    free(particle_data->micron_dust_particle_mass_array);
    free(particle_data->micron_dust_surfacedensity);
    free(particle_data->micron_particle_distance_grid);

    // Reset structure to a clean state
    memset(particle_data, 0, sizeof(ParticleData));

    LOG_DEBUG("Particle arrays freed.\n");
}