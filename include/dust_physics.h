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



/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
double calculateStokesNumber(double pradius, double sigma, disk_t *disk_params);

// calculateParticleMass függvény prototípusa
void calculateParticleMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const simulation_options_t *sim_opts);

double calculateRadialDriftBarrier(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params);
double calculateTurbulentFragmentationBarrier(double sigma, double r, double rho_p, const disk_t *disk_params);
double calculateDriftInducedFragmentationBarrier(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params);

/*	a reszecskek novekedesenek idoskalaja	*/
double calculateGrowthTimescale(double r, double eps, const disk_t *disk_params);

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double calculateDustParticleSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params);

// Porkorong sűrűségének számítása
void calculateDustSurfaceDensity(double max_param, double min_param, double rad[][2], double radmicr[][2],
                double *sigma_d, double *sigma_dm, double *massvec,
                double *massmicrvec,double *rd, double *rmic,
                const simulation_options_t *sim_opts, const disk_t *disk_params);

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void calculateDustDistance(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
                double deltat, double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params);



#endif // DUST_PHYSICS_H