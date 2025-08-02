// dust_particle.h
#ifndef DUST_PARTICLE_H
#define DUST_PARTICLE_H

/**
 * @brief Represents a single Lagrangian dust particle with its properties.
 *
 * Each instance stores the particle's ID, current distance, current size,
 * its initial (fixed) mass, and the reciprocal of its current size.
 */
typedef struct {
    int id;                   ///< Unique identifier for this particle
    double distance_au;       ///< Radial distance from the star in AU
    double distance_au_reciprocal;
    double current_size_cm;   ///< Particle size in centimeters
    double initial_mass_msun; ///< Initial (fixed) representative mass in Solar Masses
    double size_reciprocal;   ///< 1 / current_size_cm (for efficiency in some calculations)
    double drdt;

} dust_particle_t;

#endif // DUST_PARTICLE_H