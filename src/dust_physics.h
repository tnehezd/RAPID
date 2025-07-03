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


/*	alpha turbulens paraméter kiszámolása	*/
double calculate_turbulent_alpha(double r, const disk_t *disk_params);

/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
double Stokes_Number(double pradius, double sigma, disk_t *disk_params);

/*	Lokalis viszkozitas erteke	*/
double visc(double r, const disk_t *disk_params);

/*	local scale height	*/
double scale_height(double r, const disk_t *disk_params);

/*	lokális kepleri sebesség	*/
double v_kep(double r, const disk_t *disk_params);

/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r, const disk_t *disk_params);

/*	local sound speed		*/
double c_sound(double r, const disk_t *disk_params);

/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r, const disk_t *disk_params);

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r, const disk_t *disk_params);

/*	a nyomas derivaltja	*/
void dpress(disk_t *disk_params);

/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
double Coeff_3(double sigma, double r); // Feltételezve, hogy Coeff_3 a dust_physics.c-ben van

/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(disk_t *disk_params); // disk_params nem const, mert módosítva van az ugvec tagja

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
void Get_Sigmad(double max_param, double min_param, double rad[][2], double radmicr[][2], double radsec[][2],
                double *sigma_d, double *sigma_dm, double *sigma_ds, double *massvec,
                double *massmicrvec, double *masssecvec, double *rd, double *rmic, double *rs,
                const simulation_options_t *sim_opts, const disk_t *disk_params);

/*	Fuggveny a sigma, p, dp kiszamolasara	*/
// disk_params nem const, mert módosítva van a sigmavec, pressvec, dpressvec tagjai
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params);

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
                double deltat, double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params);



#endif // DUST_PHYSICS_H