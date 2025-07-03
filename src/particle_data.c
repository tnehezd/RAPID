// src/particle_data.c

#include "particle_data.h"
#include <stdio.h>
#include <stdlib.h> // malloc, free, exit

int allocate_particle_data(ParticleData_t *p_data, size_t particle_count, int is_twopop_enabled) {
    if (p_data == NULL) {
        fprintf(stderr, "ERROR [allocate_particle_data]: ParticleData_t pointer is NULL.\n");
        return 1; // Hiba
    }

    // Inicializálás NULL-ra
    p_data->radius = NULL;
    p_data->radiusmicr = NULL;
    p_data->radiussec = NULL;
    p_data->radius_rec = NULL;
    p_data->massvec = NULL;
    p_data->massmicrvec = NULL;
    p_data->masssecvec = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->partmasssecind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->sigmads = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->rsvec = NULL;
    p_data->allocated_particle_number = 0;

    if (particle_count == 0) {
        fprintf(stderr, "DEBUG [allocate_particle_data]: Particle count is 0. No particle arrays allocated.\n");
        return 0; // Sikeres, de nincs allokáció
    }

    // Fő részecske tömbök
    p_data->radius = malloc(particle_count * sizeof(*p_data->radius));
    p_data->radiusmicr = malloc(particle_count * sizeof(*p_data->radiusmicr));
    p_data->radius_rec = malloc(particle_count * sizeof(*p_data->radius_rec));
    p_data->massvec = malloc(particle_count * sizeof(double));
    p_data->massmicrvec = malloc(particle_count * sizeof(double));
    p_data->partmassind = malloc(particle_count * sizeof(*p_data->partmassind));
    p_data->partmassmicrind = malloc(particle_count * sizeof(*p_data->partmassmicrind));
    p_data->sigmad = malloc(particle_count * sizeof(double));
    p_data->sigmadm = malloc(particle_count * sizeof(double));
    p_data->rdvec = malloc(particle_count * sizeof(double));
    p_data->rmicvec = malloc(particle_count * sizeof(double));

    // Ellenőrzés
    if (!p_data->radius || !p_data->radiusmicr || !p_data->radius_rec || !p_data->massvec || !p_data->massmicrvec ||
        !p_data->partmassind || !p_data->partmassmicrind || !p_data->sigmad || !p_data->sigmadm ||
        !p_data->rdvec || !p_data->rmicvec) {
        fprintf(stderr, "ERROR [allocate_particle_data]: Primary particle array allocation failed!\n");
        free_particle_data(p_data); // Felszabadítás, ha valami elszállt
        return 1; // Hiba
    }

    // Secondary particles (csak ha twopop engedélyezve van, feltételezve, hogy a growth ehhez kapcsolódik)
    // A 4-szeres méretet a Get_Radius függvényben látottak alapján vettem.
    if (is_twopop_enabled) {
        p_data->radiussec = malloc(4 * particle_count * sizeof(*p_data->radiussec));
        p_data->masssecvec = malloc(4 * particle_count * sizeof(double));
        p_data->partmasssecind = malloc(4 * particle_count * sizeof(*p_data->partmasssecind));
        p_data->sigmads = malloc(4 * particle_count * sizeof(double));
        p_data->rsvec = malloc(4 * particle_count * sizeof(double));

        if (!p_data->radiussec || !p_data->masssecvec || !p_data->partmasssecind || !p_data->sigmads || !p_data->rsvec) {
            fprintf(stderr, "ERROR [allocate_particle_data]: Secondary particle array allocation failed!\n");
            free_particle_data(p_data); // Felszabadítás, ha valami elszállt
            return 1; // Hiba
        }
    } else {
        fprintf(stderr, "DEBUG [allocate_particle_data]: Two-population model is OFF. Secondary particle arrays not allocated.\n");
    }

    p_data->allocated_particle_number = particle_count;
    fprintf(stderr, "DEBUG [allocate_particle_data]: Particle arrays allocated for %zu particles.\n", particle_count);
    return 0; // Sikeres allokáció
}

void free_particle_data(ParticleData_t *p_data) {
    if (p_data == NULL) {
        return; // Nincs mit felszabadítani
    }

    free(p_data->radius);
    free(p_data->radiusmicr);
    free(p_data->radiussec);
    free(p_data->radius_rec);
    free(p_data->massvec);
    free(p_data->massmicrvec);
    free(p_data->masssecvec);
    free(p_data->partmassind);
    free(p_data->partmassmicrind);
    free(p_data->partmasssecind);
    free(p_data->sigmad);
    free(p_data->sigmadm);
    free(p_data->sigmads);
    free(p_data->rdvec);
    free(p_data->rmicvec);
    free(p_data->rsvec);

    // Fontos: a pointereket NULL-ra állítjuk felszabadítás után, hogy elkerüljük a dangling pointereket
    p_data->radius = NULL;
    p_data->radiusmicr = NULL;
    p_data->radiussec = NULL;
    p_data->radius_rec = NULL;
    p_data->massvec = NULL;
    p_data->massmicrvec = NULL;
    p_data->masssecvec = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->partmasssecind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->sigmads = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->rsvec = NULL;
    p_data->allocated_particle_number = 0;

    fprintf(stderr, "DEBUG [free_particle_data]: Particle arrays freed.\n");
}