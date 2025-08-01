// src/particle_data.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memset
#include "particle_data.h" // Includes dust_particle.h

/**
 * @brief Allocates memory for the ParticleData_t structure based on particle counts.
 *
 * This function creates contiguous memory blocks for Pop1 and (if enabled) Pop2 particles.
 * It also initializes the particle counts.
 *
 * @param p_data Pointer to the ParticleData_t structure to be allocated.
 * @param num_particles_pop1 The total number of particles for Population 1.
 * @param num_particles_pop2 The total number of particles for Population 2.
 * @param twopop_enabled Flag indicating if a second particle population is active.
 */
void allocate_particle_data(ParticleData_t *p_data, int num_particles_pop1, int num_particles_pop2, int twopop_enabled) {
    if (p_data == NULL) {
        fprintf(stderr, "ERROR [allocate_particle_data]: p_data pointer is NULL.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize all pointers to NULL for safety before allocation
    p_data->particles_pop1 = NULL;
    p_data->particles_pop2 = NULL;

    p_data->num_particles_pop1 = num_particles_pop1;
    p_data->num_particles_pop2 = num_particles_pop2;

    fprintf(stderr, "DEBUG [allocate_particle_data]: Allocating for Lagrangian particles. Pop1: %d, Pop2: %d, twopop_enabled: %d\n",
            num_particles_pop1, num_particles_pop2, twopop_enabled);

    // Allocate for Population 1 particles
    if (num_particles_pop1 > 0) {
        p_data->particles_pop1 = (dust_particle_t *)malloc(num_particles_pop1 * sizeof(dust_particle_t));
        if (p_data->particles_pop1 == NULL) {
            fprintf(stderr, "ERROR [allocate_particle_data]: Memory allocation failed for particles_pop1 (size: %zu bytes).\n",
                    num_particles_pop1 * sizeof(dust_particle_t));
            exit(EXIT_FAILURE);
        }
        // Initialize allocated memory to zeros (good practice)
        memset(p_data->particles_pop1, 0, num_particles_pop1 * sizeof(dust_particle_t));
        fprintf(stderr, "DEBUG [allocate_particle_data]: Allocated %d particles for Pop1.\n", num_particles_pop1);
    }

    // Allocate for Population 2 particles if twopop is enabled
    if (twopop_enabled && num_particles_pop2 > 0) {
        p_data->particles_pop2 = (dust_particle_t *)malloc(num_particles_pop2 * sizeof(dust_particle_t));
        if (p_data->particles_pop2 == NULL) {
            fprintf(stderr, "ERROR [allocate_particle_data]: Memory allocation failed for particles_pop2 (size: %zu bytes).\n",
                    num_particles_pop2 * sizeof(dust_particle_t));
            // Clean up Pop1 if Pop2 allocation fails
            free(p_data->particles_pop1);
            p_data->particles_pop1 = NULL;
            exit(EXIT_FAILURE);
        }
        memset(p_data->particles_pop2, 0, num_particles_pop2 * sizeof(dust_particle_t));
        fprintf(stderr, "DEBUG [allocate_particle_data]: Allocated %d particles for Pop2.\n", num_particles_pop2);
    }

    fprintf(stderr, "DEBUG [allocate_particle_data]: ParticleData_t structure allocated.\n");
}

/**
 * @brief Frees memory allocated for the ParticleData_t structure.
 *
 * Sets pointers to NULL after freeing to prevent dangling pointers.
 *
 * @param p_data Pointer to the ParticleData_t structure to be freed.
 */
void free_particle_data(ParticleData_t *p_data) {
    if (p_data == NULL) {
        return; // Nothing to free
    }

    fprintf(stderr, "DEBUG [free_particle_data]: Freeing Lagrangian particle data.\n");

    if (p_data->particles_pop1 != NULL) {
        free(p_data->particles_pop1);
        p_data->particles_pop1 = NULL;
    }
    if (p_data->particles_pop2 != NULL) {
        free(p_data->particles_pop2);
        p_data->particles_pop2 = NULL;
    }
    
    // Reset particle counts for clarity
    p_data->num_particles_pop1 = 0;
    p_data->num_particles_pop2 = 0;

    fprintf(stderr, "DEBUG [free_particle_data]: ParticleData_t structure freed.\n");
}