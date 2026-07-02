// src/particle_data.c
#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h> 
#include "logger.h"

int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled) {

    if (particle_data == NULL) {
        LOG_ERROR("ParticleData pointer is NULL.\n");
        return 1; 
    }

    particle_data->particle_distance_array = NULL;
    particle_data->micron_particle_distance_array = NULL;
    particle_data->dust_particle_mass_grid = NULL;
    particle_data->massmicradial_grid = NULL;
    particle_data->dust_particle_mass_array = NULL;
    particle_data->micron_dust_particle_mass_array = NULL;
    particle_data->dust_surfacedensity = NULL;
    particle_data->micron_dust_surfacedensity = NULL;
    particle_data->particle_distance_grid = NULL;
    particle_data->micron_particle_distance_grid = NULL;
    particle_data->allocated_particle_number = 0;

    if (particle_count == 0) {
        LOG_DEBUG("Particle count is 0. No particle arrays allocated.\n");
        return 0; 
    }

    particle_data->particle_distance_array = malloc(particle_count * sizeof(*particle_data->particle_distance_array));
    particle_data->micron_particle_distance_array = malloc(particle_count * sizeof(*particle_data->micron_particle_distance_array));
    particle_data->dust_particle_mass_grid = malloc(particle_count * sizeof(double));
    particle_data->massmicradial_grid = malloc(particle_count * sizeof(double));
    particle_data->dust_particle_mass_array = malloc(particle_count * sizeof(*particle_data->dust_particle_mass_array));
    particle_data->micron_dust_particle_mass_array = malloc(particle_count * sizeof(*particle_data->micron_dust_particle_mass_array));
    particle_data->dust_surfacedensity = malloc(particle_count * sizeof(double));
    particle_data->micron_dust_surfacedensity = malloc(particle_count * sizeof(double));
    particle_data->particle_distance_grid = malloc(particle_count * sizeof(double));
    particle_data->micron_particle_distance_grid = malloc(particle_count * sizeof(double));

    if (!particle_data->particle_distance_array || !particle_data->micron_particle_distance_array  || !particle_data->dust_particle_mass_grid || !particle_data->massmicradial_grid ||
        !particle_data->dust_particle_mass_array || !particle_data->micron_dust_particle_mass_array || !particle_data->dust_surfacedensity || !particle_data->micron_dust_surfacedensity ||
        !particle_data->particle_distance_grid || !particle_data->micron_particle_distance_grid) {
        LOG_ERROR("Primary particle array allocation failed!\n");
        freeParticleData(particle_data); 
        return 1; 
    }

    if (is_twopop_enabled) {
        LOG_INFO("Two-population model is enabled, but something is missing. Need to implement...\n");
    } else {
        LOG_DEBUG("Two-population model is OFF. Secondary particle arrays not allocated.\n");
    }

    particle_data->allocated_particle_number = particle_count;

    LOG_DEBUG("Particle arrays allocated for %zu particles.\n", particle_count);
    return 0;
}

void freeParticleData(ParticleData *particle_data) {

    if (particle_data == NULL) {
        return; 
    }

    free(particle_data->particle_distance_array);
    free(particle_data->micron_particle_distance_array);
    free(particle_data->dust_particle_mass_grid);
    free(particle_data->massmicradial_grid);
    free(particle_data->dust_particle_mass_array);
    free(particle_data->micron_dust_particle_mass_array);
    free(particle_data->dust_surfacedensity);
    free(particle_data->micron_dust_surfacedensity);
    free(particle_data->particle_distance_grid);
    free(particle_data->micron_particle_distance_grid);

    particle_data->particle_distance_array = NULL;
    particle_data->micron_particle_distance_array = NULL;
    particle_data->dust_particle_mass_grid = NULL;
    particle_data->massmicradial_grid = NULL;
    particle_data->dust_particle_mass_array = NULL;
    particle_data->micron_dust_particle_mass_array = NULL;
    particle_data->dust_surfacedensity = NULL;
    particle_data->micron_dust_surfacedensity = NULL;
    particle_data->particle_distance_grid = NULL;
    particle_data->micron_particle_distance_grid = NULL;
    particle_data->allocated_particle_number = 0;

    LOG_DEBUG("Particle arrays freed.\n");
}