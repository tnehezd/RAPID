// dust_physics.h

#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h>
#include "simulation_types.h" // Szükséges a simulation_options_t és disk_t struktúrákhoz
#include "particle_data.h"    // Ezt az include-ot is hozzá kell adni a dust_particle_t miatt!


// A calculate_dust_surface_density_profile prototípusa, ha az radin[][2]-t vár
void calculate_dust_surface_density_profile(
    double *output_sigma_d_grid,       // Kimenet: A kiszámított por felületi sűrűség (gridenként)
    double *output_r_grid_centers,     // Kimenet: A grid cellák középpontjainak radiális pozíciói
    const double radin[][2],           // Bemenet: Részecskék radiális pozíciója (const)
    const double *massin,              // Bemenet: Részecskék tömege (const)
    int n_particles,                   // Bemenet: A részecskék száma
    int n_grid_cells,                  // Bemenet: A grid cellák száma, amire a sűrűséget számoljuk
    const disk_t *disk_params          // Bemenet: A korong paraméterei (const)
);

/* alpha turbulens paraméter kiszámolása */
double calculate_turbulent_alpha(double r, const disk_t *disk_params);

/* kiszamolja az adott reszecskehez tartozo Stokes szamot  */
double Stokes_Number(double pradius, double sigma, disk_t *disk_params);


double a_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params);
double a_turb(double sigma, double r, double rho_p, const disk_t *disk_params);
double a_df(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/* a reszecskek novekedesenek idoskalaja   */
double tauGr(double r, double eps, const disk_t *disk_params);

/* kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK  */
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params);

// Porkorong sűrűségének számítása
// Ezt a függvényt úgy fogjuk hívni a tIntegrate-ből, hogy a ParticleData_t struktúrát kapja meg.
// Az *rd és *rmic paraméterek valószínűleg nem kellenek kimenetként, ha a disk_params->rvec-et használjuk.
// A sigmad és sigmadm is a disk_params-ba kerül majd.
void Get_Sigmad(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts);

/* Fuggveny a porszemcsek uj tavolsaganak elraktarozasara      */
// ÚJ PROTOTÍPUS:
// - dust_particle_t *particles_array: A tényleges részecske tömb (Pop1 VAGY Pop2)
// - int num_particles: A tömbben lévő részecskék száma
// - double deltat, double t: időlépés és aktuális idő
// - const simulation_options_t *sim_opts, const disk_t *disk_params: opciók és korong paraméterek
void Get_Radius(dust_particle_t *particles_array, int num_particles, double deltat, double t,
                const simulation_options_t *sim_opts, const disk_t *disk_params);


#endif // DUST_PHYSICS_H