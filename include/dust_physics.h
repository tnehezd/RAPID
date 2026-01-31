#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h> // Szükséges lehet FILE, fprintf stb.
#include "simulation_types.h" // Szükséges a SimulationOptions és DiskParameters struktúrákhoz
#include "particle_data.h"

// Extern globális változók deklarációi, ha ezeket a dust_physics.c fájl használja
// Ezeknek a definíciói valószínűleg a config.h-ban vagy a main.c-ben vannak.
// Csak akkor hagyd meg őket itt, ha a dust_physics.c fájl közvetlenül használja őket,
// és nincs beinclude-olva egy másik fájl, ami már deklarálja (pl. config.h)
// Ha a config.h includolva van a dust_physics.c-ben, és ott vannak definiálva,
// akkor ezeket az extern deklarációkat elhagyhatod innen.



/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
double calculateStokesNumber(double pradius, double sigma, DiskParameters *disk_params);

// calculateParticleMass függvény prototípusa
void calculateParticleMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const SimulationOptions *sim_opts);

double calculateRadialDriftBarrier(double sigmad, double r, double p, double dp, double rho_p, const DiskParameters *disk_params);
double calculateTurbulentFragmentationBarrier(double sigma, double r, double rho_p, const DiskParameters *disk_params);
double calculateDriftInducedFragmentationBarrier(double sigma, double r, double p, double dp, double rho_p, const DiskParameters *disk_params);

/*	a reszecskek novekedesenek idoskalaja	*/
double calculateGrowthTimescale(double r, double eps, const DiskParameters *disk_params);

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double calculateDustParticleSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const DiskParameters *disk_params);

// Porkorong sűrűségének számítása
void calculateDustSurfaceDensity(double max_param, double min_param, const ParticleData *particle_data, const SimulationOptions *sim_opts, const DiskParameters *disk_params);

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void calculateDustDistance(const char *nev, int opt, ParticleData *particle_data, double deltat, double t, int n, const SimulationOptions *sim_opts, const DiskParameters *disk_params);



#endif // DUST_PHYSICS_H