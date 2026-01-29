// src/particle_data.c

#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h> // malloc, free, exit

int allocateParticleData(ParticleData_t *p_data, size_t particle_count, int is_twopop_enabled) {
    if (p_data == NULL) {
        fprintf(stderr, "ERROR [allocateParticleData]: ParticleData_t pointer is NULL.\n");
        return 1; // Hiba
    }

    // Inicializálás NULL-ra
    p_data->radius = NULL;
    p_data->radiusmicr = NULL;
    p_data->radius_rec = NULL;
    p_data->massvec = NULL;
    p_data->massmicradial_grid = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->allocated_particle_number = 0;

    if (particle_count == 0) {
        fprintf(stderr, "DEBUG [allocateParticleData]: Particle count is 0. No particle arrays allocated.\n");
        return 0; // Sikeres, de nincs allokáció
    }

    // Fő részecske tömbök
    p_data->radius = malloc(particle_count * sizeof(*p_data->radius));
    p_data->radiusmicr = malloc(particle_count * sizeof(*p_data->radiusmicr));
    p_data->radius_rec = malloc(particle_count * sizeof(*p_data->radius_rec));
    p_data->massvec = malloc(particle_count * sizeof(double));
    p_data->massmicradial_grid = malloc(particle_count * sizeof(double));
    p_data->partmassind = malloc(particle_count * sizeof(*p_data->partmassind));
    p_data->partmassmicrind = malloc(particle_count * sizeof(*p_data->partmassmicrind));
    p_data->sigmad = malloc(particle_count * sizeof(double));
    p_data->sigmadm = malloc(particle_count * sizeof(double));
    p_data->rdvec = malloc(particle_count * sizeof(double));
    p_data->rmicvec = malloc(particle_count * sizeof(double));

    // Ellenőrzés
    if (!p_data->radius || !p_data->radiusmicr || !p_data->radius_rec || !p_data->massvec || !p_data->massmicradial_grid ||
        !p_data->partmassind || !p_data->partmassmicrind || !p_data->sigmad || !p_data->sigmadm ||
        !p_data->rdvec || !p_data->rmicvec) {
        fprintf(stderr, "ERROR [allocateParticleData]: Primary particle array allocation failed!\n");
        freeParticleData(p_data); // Felszabadítás, ha valami elszállt
        return 1; // Hiba
    }

    // Secondary particles (csak ha twopop engedélyezve van, feltételezve, hogy a growth ehhez kapcsolódik)
    // A 4-szeres méretet a calculateDustDistance függvényben látottak alapján vettem.
    if (is_twopop_enabled) {

    } else {
        fprintf(stderr, "DEBUG [allocateParticleData]: Two-population model is OFF. Secondary particle arrays not allocated.\n");
    }

    p_data->allocated_particle_number = particle_count;
    fprintf(stderr, "DEBUG [allocateParticleData]: Particle arrays allocated for %zu particles.\n", particle_count);
    return 0; // Sikeres allokáció
}

void freeParticleData(ParticleData_t *p_data) {
    if (p_data == NULL) {
        return; // Nincs mit felszabadítani
    }

    free(p_data->radius);
    free(p_data->radiusmicr);
    free(p_data->radius_rec);
    free(p_data->massvec);
    free(p_data->massmicradial_grid);
    free(p_data->partmassind);
    free(p_data->partmassmicrind);
    free(p_data->sigmad);
    free(p_data->sigmadm);
    free(p_data->rdvec);
    free(p_data->rmicvec);

    // Fontos: a pointereket NULL-ra állítjuk felszabadítás után, hogy elkerüljük a dangling pointereket
    p_data->radius = NULL;
    p_data->radiusmicr = NULL;
    p_data->radius_rec = NULL;
    p_data->massvec = NULL;
    p_data->massmicradial_grid = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->allocated_particle_number = 0;

    fprintf(stderr, "DEBUG [freeParticleData]: Particle arrays freed.\n");
}