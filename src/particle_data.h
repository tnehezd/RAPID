// PARTICLE_DATA_H
#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include "dust_particle.h" // Include the specific dust particle definition

/**
 * @brief This structure stores the dynamically allocated arrays for particle populations.
 *
 * It acts as a container for the main and micron-sized dust particle arrays,
 * which hold the evolving properties of each particle.
 * This is the ONLY place where particle data arrays are managed globally.
 */
typedef struct {
    // ÚJ és kizárólagos deklarációk az egységesített részecskestruktúra használatával
    // Ezek a tömbök fogják tartalmazni az összes részecske aktuális állapotát.
    dust_particle_t *particles_pop1;      ///< Pointer to an array of dust_particle_t for population 1 (cm-sized)
    dust_particle_t *particles_pop2;      ///< Pointer to an array of dust_particle_t for population 2 (micron-sized, if twopop)

    int num_particles_pop1;             ///< Number of particles in population 1
    int num_particles_pop2;             ///< Number of particles in population 2 (0 if twopop is disabled)


} ParticleData_t;


// Function prototypes
void allocate_particle_data(ParticleData_t *p_data, int num_particles_pop1, int num_particles_pop2, int twopop_enabled);
void free_particle_data(ParticleData_t *p_data);

#endif // PARTICLE_DATA_H