#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <stddef.h> // size_t-hoz

// Struktúra a dinamikusan allokált részecskeadatok tárolására
typedef struct {
    long double (*radius)[2]; // Changed to long double
    long double (*radiusmicr)[2]; // Changed to long double
    long double (*radius_rec)[2]; // Temp array for inverse radii calculations // Changed to long double
    long double *massvec; // Changed to long double
    long double *massmicrvec; // Changed to long double
    long double (*partmassind)[5]; // Changed to long double
    long double (*partmassmicrind)[5]; // Changed to long double
    long double *sigmad; // Changed to long double
    long double *sigmadm; // Changed to long double
    long double *rdvec; // Changed to long double
    long double *rmicvec; // Changed to long double
    size_t allocated_particle_number; // Az aktuálisan allokált részecskék száma
} ParticleData_t;

// Függvény a memóriafoglalásra
int allocate_particle_data(ParticleData_t *p_data, size_t particle_count, int is_twopop_enabled);

// Függvény a memória felszabadítására
void free_particle_data(ParticleData_t *p_data);

#endif // PARTICLE_DATA_H