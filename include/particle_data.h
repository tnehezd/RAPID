// src/particle_data.h

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <stddef.h>

typedef struct {
    double (*particle_distance_array)[2];
    double (*micron_particle_distance_array)[2];
    double *dust_particle_mass_grid;
    double *massmicradial_grid;
    double (*dust_particle_mass_array)[5];
    double (*micron_dust_particle_mass_array)[5];
    double *dust_surfacedensity;
    double *micron_dust_surfacedensity;
    double *particle_distance_grid;
    double *micron_particle_distance_grid;
    size_t allocated_particle_number;
} ParticleData;

int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled);
void freeParticleData(ParticleData *particle_data);

#endif // PARTICLE_DATA_H