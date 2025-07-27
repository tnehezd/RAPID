#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h> // Szükséges lehet FILE, fprintf stb.
#include "simulation_types.h" // Szükséges a simulation_options_t és disk_t struktúrákhoz

// Extern globális változók deklarációi, ha ezeket a dust_physics.c fájl használja
// Ezeknek a definíciói valószínűleg a config.h-ban vagy a main.c-ben vannak.
// Csak akkor hagyd meg őket itt, ha a dust_physics.c fájl közvetlenül használja őket,
// és nincs beinclude-olva egy másik fájl, ami már deklarálja (pl. config.h)
// Ha a config.h includolva van a dust_physics.c-ben, és ott vannak definiálva,
// akkor ezeket az extern deklarációkat elhagyhatod innen.

/**
 * @brief Reads and calculates derived disk parameters.
 *
 * This function initializes and sets up the fundamental parameters of the disk
 * within the `disk_t` structure. It may involve converting input parameters
 * to the appropriate internal units or calculating derived dimensionless quantities.
 *
 * @param disk_params Pointer to the `disk_t` structure where parameters are stored.
 */
void read_disk_parameters(disk_t *disk_params);


void initial_dust_surface_density_profile(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params);

/*	alpha turbulens paraméter kiszámolása	*/
double calculate_turbulent_alpha(double r, const disk_t *disk_params);

/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
double Stokes_Number(double pradius, double sigma, disk_t *disk_params);

// GetMass függvény prototípusa
void GetMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const simulation_options_t *sim_opts);

double a_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params);
double a_turb(double sigma, double r, double rho_p, const disk_t *disk_params);
double a_df(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/*	a reszecskek novekedesenek idoskalaja	*/
double tauGr(double r, double eps, const disk_t *disk_params);

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params);

// Porkorong sűrűségének számítása
void Get_Sigmad(double max_param, double min_param, double rad[][2], double radmicr[][2],
                double *sigma_d, double *sigma_dm, double *massvec,
                double *massmicrvec,double *rd, double *rmic,
                const simulation_options_t *sim_opts, const disk_t *disk_params);

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
                double deltat, double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params);



#endif // DUST_PHYSICS_H