#ifndef GAS_PHYSICS_H
#define GAS_PHYSICS_H

#include "simulation_types.h" 

/*	Calculates turbulent alpha parameter	*/
double calculateTurbulentAlpha(double r, const disk_t *disk_params);

/*	Lokalis viszkozitas erteke	*/
double calculateKinematicViscosity(double r, const disk_t *disk_params);

/*	local scale height	*/
double calculatePressureScaleHeight(double r, const disk_t *disk_params);

/*	lokális kepleri sebesség	*/
double calculateKeplerianVelocity(double r, const disk_t *disk_params);

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

/*	Fuggveny a sigma, p, dp kiszamolasara	*/
// disk_params nem const, mert módosítva van a sigmavec, pressvec, dpressvec tagjai
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params);

#endif // GAS_PHYSICS_H
