// src/particle_data.c

#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h> // malloc, free, exit

int allocateParticleData(ParticleData *particle_data, size_t particle_count, int is_twopop_enabled) {
    if (particle_data == NULL) {
        fprintf(stderr, "ERROR [allocateParticleData]: ParticleData pointer is NULL.\n");
        return 1; // Hiba
    }

    // Inicializálás NULL-ra
    particle_data->radius = NULL;
    particle_data->radiusmicr = NULL;
    particle_data->massvec = NULL;
    particle_data->massmicradial_grid = NULL;
    particle_data->partmassind = NULL;
    particle_data->partmassmicrind = NULL;
    particle_data->sigmad = NULL;
    particle_data->sigmadm = NULL;
    particle_data->rdvec = NULL;
    particle_data->rmicvec = NULL;
    particle_data->allocated_particle_number = 0;

    if (particle_count == 0) {
        fprintf(stderr, "DEBUG [allocateParticleData]: Particle count is 0. No particle arrays allocated.\n");
        return 0; // Sikeres, de nincs allokáció
    }

    // Fő részecske tömbök
    particle_data->radius = malloc(particle_count * sizeof(*particle_data->radius));
    particle_data->radiusmicr = malloc(particle_count * sizeof(*particle_data->radiusmicr));
    particle_data->massvec = malloc(particle_count * sizeof(double));
    particle_data->massmicradial_grid = malloc(particle_count * sizeof(double));
    particle_data->partmassind = malloc(particle_count * sizeof(*particle_data->partmassind));
    particle_data->partmassmicrind = malloc(particle_count * sizeof(*particle_data->partmassmicrind));
    particle_data->sigmad = malloc(particle_count * sizeof(double));
    particle_data->sigmadm = malloc(particle_count * sizeof(double));
    particle_data->rdvec = malloc(particle_count * sizeof(double));
    particle_data->rmicvec = malloc(particle_count * sizeof(double));

    // Ellenőrzés
    if (!particle_data->radius || !particle_data->radiusmicr  || !particle_data->massvec || !particle_data->massmicradial_grid ||
        !particle_data->partmassind || !particle_data->partmassmicrind || !particle_data->sigmad || !particle_data->sigmadm ||
        !particle_data->rdvec || !particle_data->rmicvec) {
        fprintf(stderr, "ERROR [allocateParticleData]: Primary particle array allocation failed!\n");
        freeParticleData(particle_data); // Felszabadítás, ha valami elszállt
        return 1; // Hiba
    }

    // Secondary particles (csak ha twopop engedélyezve van, feltételezve, hogy a growth ehhez kapcsolódik)
    // A 4-szeres méretet a calculateDustDistance függvényben látottak alapján vettem.
    if (is_twopop_enabled) {

    } else {
        fprintf(stderr, "DEBUG [allocateParticleData]: Two-population model is OFF. Secondary particle arrays not allocated.\n");
    }

    particle_data->allocated_particle_number = particle_count;
    fprintf(stderr, "DEBUG [allocateParticleData]: Particle arrays allocated for %zu particles.\n", particle_count);
    return 0; // Sikeres allokáció
}

void freeParticleData(ParticleData *particle_data) {
    if (particle_data == NULL) {
        return; // Nincs mit felszabadítani
    }

    free(particle_data->radius);
    free(particle_data->radiusmicr);
    free(particle_data->massvec);
    free(particle_data->massmicradial_grid);
    free(particle_data->partmassind);
    free(particle_data->partmassmicrind);
    free(particle_data->sigmad);
    free(particle_data->sigmadm);
    free(particle_data->rdvec);
    free(particle_data->rmicvec);

    // Fontos: a pointereket NULL-ra állítjuk felszabadítás után, hogy elkerüljük a dangling pointereket
    particle_data->radius = NULL;
    particle_data->radiusmicr = NULL;
    particle_data->massvec = NULL;
    particle_data->massmicradial_grid = NULL;
    particle_data->partmassind = NULL;
    particle_data->partmassmicrind = NULL;
    particle_data->sigmad = NULL;
    particle_data->sigmadm = NULL;
    particle_data->rdvec = NULL;
    particle_data->rmicvec = NULL;
    particle_data->allocated_particle_number = 0;

    fprintf(stderr, "DEBUG [freeParticleData]: Particle arrays freed.\n");
}