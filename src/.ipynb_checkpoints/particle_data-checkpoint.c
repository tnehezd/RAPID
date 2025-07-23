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
    p_data->radius_rec = NULL;
    p_data->massvec = NULL;
    p_data->massmicrvec = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->allocated_particle_number = 0;

    if (particle_count == 0) {
        fprintf(stderr, "DEBUG [allocate_particle_data]: Particle count is 0. No particle arrays allocated.\n");
        return 0; // Sikeres, de nincs allokáció
    }

    // Fő részecske tömbök
    // Changed to long double allocation
    p_data->radius = malloc(particle_count * sizeof(long double));
    p_data->radiusmicr = malloc(particle_count * sizeof(long double));
    p_data->radius_rec = malloc(particle_count * sizeof(long double));

    p_data->massvec = malloc(particle_count * sizeof(long double)); 
    p_data->massmicrvec = malloc(particle_count * sizeof(long double)); 
    p_data->partmassind = malloc(particle_count * sizeof(long double)); // Assuming partmassind holds long double
    p_data->partmassmicrind = malloc(particle_count * sizeof(long double)); // Assuming partmassmicrind holds long double

    // Changed to long double allocation
    p_data->sigmad = malloc(particle_count * sizeof(long double));
    p_data->sigmadm = malloc(particle_count * sizeof(long double));
    p_data->rdvec = malloc(particle_count * sizeof(long double));
    p_data->rmicvec = malloc(particle_count * sizeof(long double));

    // Ellenőrzés
    if (!p_data->radius || !p_data->radiusmicr || !p_data->radius_rec || !p_data->massvec || !p_data->massmicrvec ||
        !p_data->partmassind || !p_data->partmassmicrind || !p_data->sigmad || !p_data->sigmadm ||
        !p_data->rdvec || !p_data->rmicvec) {
        fprintf(stderr, "ERROR [allocate_particle_data]: Primary particle array allocation failed!\n");
        free_particle_data(p_data); // Felszabadítás, ha valami elszállt
        return 1; // Hiba
    }

    // Secondary particles (csak ha twopop engedélyezve van, feltételezve, hogy a growth ehhez kapcsolódik)
    // A kódban lévő "if (is_twopop_enabled) {} else {}" blokk üres, így a debug üzenet is csak az else ágon jelenik meg.
    // Ez valószínűleg egy helyőrző, és nincs szükség allokációra itt.
    if (is_twopop_enabled) {
        // Ha lenne külön allokáció a twopop-hoz, az itt történne.
    } else {
        fprintf(stderr, "DEBUG [allocate_particle_data]: Two-population model is OFF. Secondary particle arrays (radiusmicr, massmicrvec, partmassmicrind, sigmadm, rmicvec) are allocated but might not be fully utilized.\n");
        // Megjegyzés: A radiusmicr, massmicrvec, partmassmicrind, sigmadm, rmicvec már allokálva van a fő blokkban.
        // Ez a debug üzenet kicsit félrevezető lehet, de a funkcionalitást nem befolyásolja.
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
    free(p_data->radius_rec);
    free(p_data->massvec);
    free(p_data->massmicrvec);
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
    p_data->massmicrvec = NULL;
    p_data->partmassind = NULL;
    p_data->partmassmicrind = NULL;
    p_data->sigmad = NULL;
    p_data->sigmadm = NULL;
    p_data->rdvec = NULL;
    p_data->rmicvec = NULL;
    p_data->allocated_particle_number = 0;

    fprintf(stderr, "DEBUG [free_particle_data]: Particle arrays freed.\n");
}