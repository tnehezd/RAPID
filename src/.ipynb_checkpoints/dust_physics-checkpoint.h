#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h> // Szükséges lehet FILE, fprintf stb.
#include "simulation_types.h" // Szükséges a simulation_options_t és disk_t struktúrákhoz

// Extern globális változók deklarációi, ha ezeket a dust_physics.c fájl használja
// Ha a config.h include-olva van a dust_physics.c-ben, és ott vannak definiálva,
// akkor ezeket az extern deklarációkat elhagyhatod innen.
// Ha itt kellenek, akkor legyenek long double-ök, ha a config.h-ban is azok.
// extern FILE *fout2; // Ha a fout2 itt van extern-ként deklarálva és nem a main.c-ben/config.h-ban

/*	alpha turbulens paraméter kiszámolása	*/
long double calculate_turbulent_alpha(long double r, const disk_t *disk_params);

/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
long double Stokes_Number(long double pradius, long double sigma, disk_t *disk_params);

/*	Lokalis viszkozitas erteke	*/
long double visc(long double r, const disk_t *disk_params);

/*	local scale height	*/
long double scale_height(long double r, const disk_t *disk_params);

/*	lokális kepleri sebesség	*/
long double v_kep(long double r, const disk_t *disk_params);

/*	lokalis kepleri korfrekvencia	*/
long double kep_freq(long double r, const disk_t *disk_params);

/*	local sound speed		*/
long double c_sound(long double r, const disk_t *disk_params);

/*	Suruseg a midplane-ben	*/
long double rho_mp(long double sigma, long double r, const disk_t *disk_params);

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
long double press(long double sigma, long double r, const disk_t *disk_params);

/*	a nyomas derivaltja	*/
void dpress(disk_t *disk_params); // disk_params nem const, mert módosítva van a dpressvec tagja

/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
long double Coeff_3(long double sigma, long double r); // Feltételezve, hogy Coeff_3 a dust_physics.c-ben van

/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(disk_t *disk_params); // disk_params nem const, mert módosítva van az ugvec tagja

// GetMass függvény prototípusa - MINDEN DOUBLE LONG DOUBLE-RE MÓDOSÍTVA
void GetMass(int n, long double (*partmassind)[5], int indii, int indio, int indoi, int indoo, int indoi_buff, int indoo_buff, long double *massiout, long double *massoout, long double *massoout_buff, const simulation_options_t *sim_opts);

long double a_drift(long double sigmad, long double r, long double p, long double dp, long double rho_p, const disk_t *disk_params);
long double a_turb(long double sigma, long double r, long double rho_p, const disk_t *disk_params);
long double a_df(long double sigma, long double r, long double p, long double dp, long double rho_p, const disk_t *disk_params);

/*	a reszecskek novekedesenek idoskalaja	*/
long double tauGr(long double r, long double eps, const disk_t *disk_params);

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
long double getSize(long double prad, long double pdens, long double sigma, long double sigmad, long double y, long double p, long double dpress_val, long double dt, const disk_t *disk_params);

// Porkorong sűrűségének számítása - MINDEN DOUBLE LONG DOUBLE-RE MÓDOSÍTVA
void Get_Sigmad(long double max_param, long double min_param,
                long double rad[][2], long double radmicr[][2],
                long double *massvec,
                long double *massmicrvec,
                long double *sigma_d_out,
                long double *sigma_dm_out,
                long double *rd_vec, long double *rmic_vec,
                const simulation_options_t *sim_opts, const disk_t *disk_params);

/*	Fuggveny a sigma, p, dp kiszamolasara	*/
// disk_params nem const, mert módosítva van a sigmavec, pressvec, dpressvec tagjai
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params);

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
// MINDEN DOUBLE LONG DOUBLE-RE MÓDOSÍTVA
void Get_Radius(const char *nev, int opt, long double radius[][2], const long double *sigmad, const long double *rdvec,
                long double deltat, long double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params);


#endif // DUST_PHYSICS_H