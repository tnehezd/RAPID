// src/particle_data.c
#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h> 


/**
 * @brief Allocates memory for a 2D radial x vertical dust particle grid.
 *
 * @param structured_data Pointer to the StructuredParticleData to initialize.
 * @param n_r Number of radial grid points.
 * @param n_z Number of vertical cells per radial point.
 * @return 0 on success, non-zero on failure.
 */
int allocateStructuredParticleData(StructuredParticleData *structured_data, size_t n_r, size_t n_z) {
    if (!structured_data || n_r == 0 || n_z == 0) {
        fprintf(stderr, "ERROR [allocateStructuredParticleData]: Invalid arguments.\n");
        return 1;
    }

    structured_data->n_r = n_r;
    structured_data->n_z = n_z;

    // allocate rows
    structured_data->particles = malloc(n_r * sizeof(DustParticle*));
    if (!structured_data->particles) {
        fprintf(stderr, "ERROR [allocateStructuredParticleData]: Failed to allocate particle rows.\n");
        return 1;
    }

    // allocate columns for each row
    for (size_t i = 0; i < n_r; i++) {
        structured_data->particles[i] = malloc(n_z * sizeof(DustParticle));
        if (!structured_data->particles[i]) {
            fprintf(stderr, "ERROR [allocateStructuredParticleData]: Failed to allocate particle row %zu.\n", i);
            // free previous rows
            for (size_t j = 0; j < i; j++) free(structured_data->particles[j]);
            free(structured_data->particles);
            structured_data->particles = NULL;
            return 1;
        }
    }

    // optional: initialize particles with zeros
    for (size_t i = 0; i < n_r; i++)
        for (size_t j = 0; j < n_z; j++) {
            structured_data->particles[i][j].index = 0;
            structured_data->particles[i][j].r_au = 0.0;
            structured_data->particles[i][j].z_au = 0.0;
            structured_data->particles[i][j].mass_g = 0.0;
        }

    fprintf(stderr, "DEBUG [allocateStructuredParticleData]: Allocated %zu x %zu dust particles.\n", n_r, n_z);
    return 0;
}

/**
 * @brief Frees memory of a 2D radial x vertical dust particle grid.
 *
 * @param structured_data Pointer to the StructuredParticleData to free.
 */
void freeStructuredParticleData(StructuredParticleData *structured_data) {
    if (!structured_data || !structured_data->particles) return;

    for (size_t i = 0; i < structured_data->n_r; i++)
        free(structured_data->particles[i]);
    free(structured_data->particles);

    structured_data->particles = NULL;
    structured_data->n_r = 0;
    structured_data->n_z = 0;

    fprintf(stderr, "DEBUG [freeStructuredParticleData]: Freed structured dust particle grid.\n");
}