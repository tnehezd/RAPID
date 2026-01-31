// src/particle_data.h

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <stddef.h> // size_t-hoz

// Struktúra a dinamikusan allokált részecskeadatok tárolására
typedef struct {
    double (*particle_distance_array)[2];
    double (*micron_particle_distance_array)[2];
    double *particle_mass_array;
    double *massmicradial_grid;
    double (*partmassind)[5];
    double (*partmassmicrind)[5];

    double *sigmad;
    double *sigmadm;
    double *rdvec;
    double *rmicvec;
    size_t allocated_particle_number; // Az aktuálisan allokált részecskék száma
} ParticleData;

// Függvény a memóriafoglalásra
int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled);

// Függvény a memória felszabadítására
void freeParticleData(ParticleData *particle_data);

#endif // PARTICLE_DATA_H