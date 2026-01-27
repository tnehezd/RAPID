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
double calculateKeplerianFrequency(double r, const disk_t *disk_params);

/*	local sound speed		*/
double c_sound(double r, const disk_t *disk_params);

/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r, const disk_t *disk_params);

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double calculateGasPressure(double sigma, double r, const disk_t *disk_params);

/*	a nyomas derivaltja	*/
void calculateGasPressureGradient(disk_t *disk_params);

/*	calculateGasRadialVelocity kiszamolasahoz eltarolt koefficiens	*/
double coefficientForGasRadialVelocity(double sigma, double r); // Feltételezve, hogy coefficientForGasRadialVelocity a dust_physics.c-ben van

/*	calculateGasRadialVelocity = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void calculateGasRadialVelocity(disk_t *disk_params); // disk_params nem const, mert módosítva van az ugvec tagja

/*	Fuggveny a sigma, p, dp kiszamolasara	*/
// disk_params nem const, mert módosítva van a sigmavec, pressvec, dpressvec tagjai
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params);

#endif // GAS_PHYSICS_H
