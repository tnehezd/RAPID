// src/particle_data.h

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <stddef.h> // size_t-hoz

// Struktúra a dinamikusan allokált részecskeadatok tárolására
typedef struct {
    double (*radius)[2];
    double (*radiusmicr)[2];
    double (*radius_rec)[2]; // Temp array for inverse radii calculations
    double *massvec;
    double *massmicrvec;
    double (*partmassind)[5];
    double (*partmassmicrind)[5];
    double *sigmad;
    double *sigmadm;
    double *rdvec;
    double *rmicvec;
    size_t allocated_particle_number; // Az aktuálisan allokált részecskék száma
} ParticleData_t;

// Függvény a memóriafoglalásra
int allocate_particle_data(ParticleData_t *p_data, size_t particle_count, int is_twopop_enabled);

// Függvény a memória felszabadítására
void free_particle_data(ParticleData_t *p_data);

#endif // PARTICLE_DATA_H